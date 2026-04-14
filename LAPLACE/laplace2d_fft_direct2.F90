program laplace2d_fft_direct_batched
  use iso_c_binding
  use cudafor
  use cufft
  implicit none

  ! --------------------------------------------------------------------------
  ! PARAMETERS
  ! --------------------------------------------------------------------------
  integer, parameter :: Nx  = 4096 + 2
  integer, parameter :: Ny  = 4096 + 2
  integer, parameter :: nxi = Nx - 2
  integer, parameter :: nyi = Ny - 2
  integer, parameter :: nxe = 2 * (nxi + 1)
  integer, parameter :: nye = 2 * (nyi + 1)

  real(8), parameter :: pi = acos(-1.0d0)

  ! --------------------------------------------------------------------------
  ! VARIABLES
  ! --------------------------------------------------------------------------
  real(8) :: dx, dy, x, y, start, finish
  real(8) :: maxErr, exactVal, scale2d
  real(8), allocatable :: u(:,:), spec(:,:)
  complex(8), allocatable :: workx(:,:), worky(:,:)

  integer(c_int) :: planx, plany, ierr
  integer(c_int), dimension(1) :: n1, inembed1, onembed1
  integer :: i, j

  print *, "--------------------------------------------"
  print *, " Direct Laplace (batched 1D DST-I via cuFFT)"
  print *, "--------------------------------------------"

  dx = 1.d0 / dble(Nx - 1)
  dy = 1.d0 / dble(Ny - 1)
  scale2d = 4.d0 / dble((nxi + 1) * (nyi + 1))

  allocate(u(Nx,Ny), spec(nxi,nyi))
  u    = 0.d0
  spec = 0.d0

  call cpu_time(start)

  ! --------------------------------------------------------------------------
  ! STEP 1: BUILD INTERIOR RHS FROM BC
  ! --------------------------------------------------------------------------
  print *, "[1] Build RHS from BC"

  do i = 1, nxi
    x = dble(i) * dx
    spec(i,1)   = spec(i,1)   - sin(pi*x) / dy**2
    spec(i,nyi) = spec(i,nyi) - sin(pi*x) * exp(-pi) / dy**2
  end do

  !$acc enter data copyin(spec)

  ! --------------------------------------------------------------------------
  ! STEP 2: DST IN X (batched over rows)
  ! --------------------------------------------------------------------------
  print *, "[2] DST in x"

  allocate(workx(nxe, nyi))
  !$acc enter data create(workx)

  !$acc parallel loop collapse(2) present(workx)
  do j = 1, nyi
    do i = 1, nxe
      workx(i,j) = (0.d0, 0.d0)
    end do
  end do

  ! Build odd extension along x on device:
  ! workx(2:nxi+1,j) = spec(1:nxi,j)
  ! workx(nxe-i+1,j) = -spec(i,j)
  !$acc parallel loop collapse(2) present(workx,spec)
  do j = 1, nyi
    do i = 1, nxi
      workx(i+1, j)         = cmplx(spec(i,j), 0.d0, kind=8)
      workx(nxe-i+1, j)     = cmplx(-spec(i,j), 0.d0, kind=8)
    end do
  end do

  n1(1)       = nxe
  inembed1(1) = nxe
  onembed1(1) = nxe

  ierr = cufftPlanMany(planx, 1, n1, inembed1, 1, nxe, onembed1, 1, nxe, CUFFT_Z2Z, nyi)
  if (ierr /= 0) stop "cufftPlanMany(planx) failed"

  !$acc host_data use_device(workx)
    ierr = cufftExecZ2Z(planx, workx, workx, CUFFT_FORWARD)
  !$acc end host_data
  if (ierr /= 0) stop "FFT x forward failed"

  ! Extract DST-I coefficients back into spec:
  ! DST = -0.5 * imag(FFT coefficient)
  !$acc parallel loop collapse(2) present(workx,spec)
  do j = 1, nyi
    do i = 1, nxi
      spec(i,j) = -0.5d0 * aimag(workx(i+1, j))
    end do
  end do

  ierr = cufftDestroy(planx)
  if (ierr /= 0) stop "cufftDestroy(planx) failed"

  !$acc exit data delete(workx)
  deallocate(workx)

  ! --------------------------------------------------------------------------
  ! STEP 3: DST IN Y (batched over columns)
  ! --------------------------------------------------------------------------
  print *, "[3] DST in y"

  allocate(worky(nye, nxi))
  !$acc enter data create(worky)

  !$acc parallel loop collapse(2) present(worky)
  do i = 1, nxi
    do j = 1, nye
      worky(j,i) = (0.d0, 0.d0)
    end do
  end do

  ! Build odd extension along y on device:
  !$acc parallel loop collapse(2) present(worky,spec)
  do i = 1, nxi
    do j = 1, nyi
      worky(j+1, i)         = cmplx(spec(i,j), 0.d0, kind=8)
      worky(nye-j+1, i)     = cmplx(-spec(i,j), 0.d0, kind=8)
    end do
  end do

  n1(1)       = nye
  inembed1(1) = nye
  onembed1(1) = nye

  ierr = cufftPlanMany(plany, 1, n1, inembed1, 1, nye, onembed1, 1, nye, CUFFT_Z2Z, nxi)
  if (ierr /= 0) stop "cufftPlanMany(plany) failed"

  !$acc host_data use_device(worky)
    ierr = cufftExecZ2Z(plany, worky, worky, CUFFT_FORWARD)
  !$acc end host_data
  if (ierr /= 0) stop "FFT y forward failed"

  !$acc parallel loop collapse(2) present(worky,spec)
  do i = 1, nxi
    do j = 1, nyi
      spec(i,j) = -0.5d0 * aimag(worky(j+1, i))
    end do
  end do

  ! --------------------------------------------------------------------------
  ! STEP 4: SPECTRAL DIVISION ON GPU
  ! --------------------------------------------------------------------------
  print *, "[4] Spectral solve on GPU"

  !$acc parallel loop collapse(2) present(spec)
  do j = 1, nyi
    do i = 1, nxi
      spec(i,j) = spec(i,j) / &
        ( 2.d0 * (cos(pi * dble(i) / dble(nxi + 1)) - 1.d0) / dx**2 + &
          2.d0 * (cos(pi * dble(j) / dble(nyi + 1)) - 1.d0) / dy**2 )
    end do
  end do

  ! --------------------------------------------------------------------------
  ! STEP 5: INVERSE DST IN Y
  ! --------------------------------------------------------------------------
  print *, "[5] Inverse DST in y"

  !$acc parallel loop collapse(2) present(worky)
  do i = 1, nxi
    do j = 1, nye
      worky(j,i) = (0.d0, 0.d0)
    end do
  end do

  !$acc parallel loop collapse(2) present(worky,spec)
  do i = 1, nxi
    do j = 1, nyi
      worky(j+1, i)         = cmplx(spec(i,j), 0.d0, kind=8)
      worky(nye-j+1, i)     = cmplx(-spec(i,j), 0.d0, kind=8)
    end do
  end do

  !$acc host_data use_device(worky)
    ierr = cufftExecZ2Z(plany, worky, worky, CUFFT_FORWARD)
  !$acc end host_data
  if (ierr /= 0) stop "FFT y inverse/DST failed"

  !$acc parallel loop collapse(2) present(worky,spec)
  do i = 1, nxi
    do j = 1, nyi
      spec(i,j) = -0.5d0 * aimag(worky(j+1, i))
    end do
  end do

  ierr = cufftDestroy(plany)
  if (ierr /= 0) stop "cufftDestroy(plany) failed"

  !$acc exit data delete(worky)
  deallocate(worky)

  ! --------------------------------------------------------------------------
  ! STEP 6: INVERSE DST IN X
  ! --------------------------------------------------------------------------
  print *, "[6] Inverse DST in x"

  allocate(workx(nxe, nyi))
  !$acc enter data create(workx)

  !$acc parallel loop collapse(2) present(workx)
  do j = 1, nyi
    do i = 1, nxe
      workx(i,j) = (0.d0, 0.d0)
    end do
  end do

  !$acc parallel loop collapse(2) present(workx,spec)
  do j = 1, nyi
    do i = 1, nxi
      workx(i+1, j)         = cmplx(spec(i,j), 0.d0, kind=8)
      workx(nxe-i+1, j)     = cmplx(-spec(i,j), 0.d0, kind=8)
    end do
  end do

  n1(1)       = nxe
  inembed1(1) = nxe
  onembed1(1) = nxe

  ierr = cufftPlanMany(planx, 1, n1, inembed1, 1, nxe, onembed1, 1, nxe, CUFFT_Z2Z, nyi)
  if (ierr /= 0) stop "cufftPlanMany(planx 2) failed"

  !$acc host_data use_device(workx)
    ierr = cufftExecZ2Z(planx, workx, workx, CUFFT_FORWARD)
  !$acc end host_data
  if (ierr /= 0) stop "FFT x inverse/DST failed"

  !$acc parallel loop collapse(2) present(workx,spec)
  do j = 1, nyi
    do i = 1, nxi
      spec(i,j) = -0.5d0 * aimag(workx(i+1, j))
    end do
  end do

  ierr = cufftDestroy(planx)
  if (ierr /= 0) stop "cufftDestroy(planx 2) failed"

  !$acc exit data delete(workx)
  deallocate(workx)

  ! --------------------------------------------------------------------------
  ! STEP 7: COPY BACK INTERIOR SOLUTION
  ! --------------------------------------------------------------------------
  print *, "[7] Copy back interior"

  !$acc update self(spec)
  !$acc exit data delete(spec)

  do j = 1, nyi
    do i = 1, nxi
      u(i+1,j+1) = spec(i,j) * scale2d
    end do
  end do

  ! --------------------------------------------------------------------------
  ! STEP 8: APPLY BOUNDARY CONDITIONS
  ! --------------------------------------------------------------------------
  do i = 1, Nx
    x = dble(i-1) * dx
    u(i,1)  = sin(pi*x)
    u(i,Ny) = sin(pi*x) * exp(-pi)
  end do
  u(1,:)  = 0.d0
  u(Nx,:) = 0.d0

  call cpu_time(finish)

  ! --------------------------------------------------------------------------
  ! STEP 9: VALIDATION
  ! --------------------------------------------------------------------------
  print *, "[8] Validation"

  maxErr = 0.d0
  do j = 1, Ny
    y = dble(j-1) * dy
    do i = 1, Nx
      x = dble(i-1) * dx
      exactVal = sin(pi*x) * exp(-pi*y)
      maxErr = max(maxErr, abs(u(i,j) - exactVal))
    end do
  end do

  print *, "--------------------------------------------"
  print *, "Max difference from exact solution = ", maxErr
  print *, "Potential at (x=0.5,y=0.5) ~ ", u((Nx-1)/2, (Ny-1)/2)
  print *, "--------------------------------------------"
  print *, "Total time = ", finish - start

  ! Optional output
   open(unit=10, file='laplace_solution.dat', status='replace')
   do j = 1, Ny
      do i = 1, Nx
         y = (j-1)*dy
         x = (i-1)*dx
         write(10,'(3E16.8)') x, y, u(i,j)
      end do
      write(10,*)
   end do
   close(10)

  deallocate(u, spec)

end program laplace2d_fft_direct_batched
