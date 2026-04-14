program laplace2d_fft_direct
  use iso_c_binding
  use cudafor
  use cufft
  implicit none

  ! --------------------------------------------------------------------------
  ! PARAMETERS
  ! --------------------------------------------------------------------------
  integer, parameter :: Nx  = 2048 + 2
  integer, parameter :: Ny  = 2048 + 2
  integer, parameter :: nxi = Nx - 2
  integer, parameter :: nyi = Ny - 2
  integer, parameter :: nxe = 2 * (nxi + 1)
  integer, parameter :: nye = 2 * (nyi + 1)

  real(8), parameter :: pi = acos(-1.0d0)

  ! --------------------------------------------------------------------------
  ! VARIABLES
  ! --------------------------------------------------------------------------
  real(8) :: dx, dy, x, y, start, finish
  real(8) :: maxErr, exactVal, scale
  real(8) :: lamx, lamy, denom

  real(8),    allocatable :: u(:,:)
  complex(8), allocatable :: ext(:,:)

  integer(c_int) :: plan, ierr
  integer :: i, j

  print *, "--------------------------------------------"
  print *, " Direct Laplace (DST-I via cuFFT, GPU solve)"
  print *, "--------------------------------------------"

  dx    = 1.d0 / dble(Nx - 1)
  dy    = 1.d0 / dble(Ny - 1)
  scale = 1.d0 / dble(nxe * nye)

  allocate(u(Nx,Ny), ext(nxe,nye))

  u   = 0.d0
  ext = (0.d0, 0.d0)

  call cpu_time(start)

  ! --------------------------------------------------------------------------
  ! STEP 1: BUILD FULL ODD-ODD EXTENSION DIRECTLY
  ! --------------------------------------------------------------------------
  print *, "[1] Build full antisymmetric extension"

  ! Interior RHS comes from the known Dirichlet BCs
  do i = 1, nxi
    x = dble(i) * dx

    ! j = 1 interior row gets contribution from bottom boundary
    ext(i+1, 2) = cmplx(-sin(pi*x) / dy**2, 0.d0, kind=8)

    ! j = nyi interior row gets contribution from top boundary
    ext(i+1, nyi+1) = cmplx(-sin(pi*x) * exp(-pi) / dy**2, 0.d0, kind=8)
  end do

  ! Odd reflection in x
  do j = 1, nye
    do i = 1, nxi
      ext(nxe - i + 1, j) = -ext(i + 1, j)
    end do
  end do

  ! Odd reflection in y
  do j = 1, nyi
    do i = 1, nxe
      ext(i, nye - j + 1) = -ext(i, j + 1)
    end do
  end do

  !$acc enter data copyin(ext) create(u)

  ! --------------------------------------------------------------------------
  ! STEP 2: FFT FORWARD
  ! --------------------------------------------------------------------------
  print *, "[2] FFT forward to wavenumber space"

  ierr = cufftPlan2d(plan, nxe, nye, CUFFT_Z2Z)
  if (ierr /= 0) stop "cufftPlan2d failed"

  !$acc host_data use_device(ext)
  ierr = cufftExecZ2Z(plan, ext, ext, CUFFT_FORWARD)
  !$acc end host_data
  if (ierr /= 0) stop "FFT forward failed"

  ! --------------------------------------------------------------------------
  ! STEP 3: SPECTRAL SOLVE ON GPU
  ! --------------------------------------------------------------------------
  print *, "[3] Spectral solve on GPU"

  ! Zero only the DST null lines (where solution is given by BCs)
  !$acc parallel loop present(ext)
  do i = 1, nxe
    ext(i,1)     = (0.d0, 0.d0)
    ext(i,nyi+2) = (0.d0, 0.d0)
  end do

  !$acc parallel loop present(ext)
  do j = 1, nye
    ext(1,j)     = (0.d0, 0.d0)
    ext(nxi+2,j) = (0.d0, 0.d0)
  end do

  ! Divide the independent block by the discrete Dirichlet eigenvalues (wavenumbers)
  !$acc parallel loop collapse(2) private(lamx,lamy,denom) present(ext)
  do j = 1, nyi
    do i = 1, nxi
      lamx  = 2.d0 * (cos(pi * dble(i) / dble(nxi + 1)) - 1.d0) / dx**2
      lamy  = 2.d0 * (cos(pi * dble(j) / dble(nyi + 1)) - 1.d0) / dy**2
      denom = lamx + lamy
      ext(i+1,j+1) = ext(i+1,j+1) / denom
    end do
  end do

  ! Rebuild full odd-odd symmetry in spectral space
  !$acc parallel loop collapse(2) present(ext)
  do j = 1, nyi
    do i = 1, nxi
      ext(nxe - i + 1, j + 1)     = -ext(i + 1, j + 1)
      ext(i + 1, nye - j + 1)     = -ext(i + 1, j + 1)
      ext(nxe - i + 1, nye - j + 1) = ext(i + 1, j + 1)
    end do
  end do

  ! --------------------------------------------------------------------------
  ! STEP 4: FFT INVERSE
  ! --------------------------------------------------------------------------
  print *, "[4] FFT inverse to physical space"

  !$acc host_data use_device(ext)
  ierr = cufftExecZ2Z(plan, ext, ext, CUFFT_INVERSE)
  !$acc end host_data
  if (ierr /= 0) stop "FFT inverse failed"

  ! --------------------------------------------------------------------------
  ! STEP 5: EXTRACT SOLUTION ON GPU
  ! --------------------------------------------------------------------------
  print *, "[5] Extract solution"

  !$acc parallel loop collapse(2) present(ext,u)
  do j = 1, nyi
    do i = 1, nxi
      u(i+1,j+1) = real(ext(i+1,j+1), kind=8) * scale
    end do
  end do

  !$acc update self(u)

  !$acc exit data delete(ext,u)

  ! --------------------------------------------------------------------------
  ! STEP 6: APPLY BOUNDARY CONDITIONS
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
  ! STEP 7: VALIDATION
  ! --------------------------------------------------------------------------
  print *, "[6] Validation"

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

  print *, " Total time = ", finish - start

  ierr = cufftDestroy(plan)
  if (ierr /= 0) stop "cufftDestroy failed"

  ! Write solution to file
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

  deallocate(u, ext)

end program laplace2d_fft_direct
