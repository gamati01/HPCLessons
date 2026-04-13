program laplace_jacobi_gpu_pointer
  implicit none

  !-----------------------------------------------------------------------
  ! Parameters
  !-----------------------------------------------------------------------
  integer, parameter :: Nx = 512+2
  integer, parameter :: Ny = 512+2
  integer, parameter :: maxIter = 500000
  double precision, parameter :: tol = 1.0d-8

  !-----------------------------------------------------------------------
  ! Variables
  !-----------------------------------------------------------------------
  double precision, allocatable, target :: u_buf(:,:), u_new_buf(:,:)
  double precision, pointer, contiguous :: u(:,:), u_new(:,:), tmp(:,:)

  double precision :: dx, dy, x, y, pi
  double precision :: err, diff, exactVal, maxErr
  double precision :: start, finish
  integer :: i, j, iter

  !-----------------------------------------------------------------------
  ! Allocate contiguous storage
  !-----------------------------------------------------------------------
  allocate(u_buf(Nx,Ny), u_new_buf(Nx,Ny))

  ! Associate pointers
  u     => u_buf
  u_new => u_new_buf

  call cpu_time(start)

  !-----------------------------------------------------------------------
  ! Grid
  !-----------------------------------------------------------------------
  pi = acos(-1.0d0)
  dx = 1.0d0 / dble(Nx - 1)
  dy = 1.0d0 / dble(Ny - 1)

  !-----------------------------------------------------------------------
  ! Move data to GPU
  !-----------------------------------------------------------------------
  !$acc enter data create(u_buf, u_new_buf)

  !$acc parallel loop collapse(2) present(u_buf, u_new_buf)
  do j = 1, Ny
    do i = 1, Nx
      u(i,j)     = 0.0d0
      u_new(i,j) = 0.0d0
    end do
  end do

  !-----------------------------------------------------------------------
  ! Boundary conditions
  !-----------------------------------------------------------------------
  !$acc parallel loop present(u_buf, u_new_buf)
  do i = 1, Nx
    x = (i - 1)*dx
    u(i,1)     = sin(pi*x)
    u_new(i,1) = u(i,1)

    u(i,Ny)     = sin(pi*x)*exp(-pi)
    u_new(i,Ny) = u(i,Ny)
  end do

  !-----------------------------------------------------------------------
  ! Jacobi iteration (pointer swap)
  !-----------------------------------------------------------------------
  do iter = 1, maxIter
    err = 0.0d0

    !$acc parallel loop gang vector collapse(2) present(u_buf, u_new_buf) reduction(max:err)
    do j = 2, Ny-1
      do i = 2, Nx-1
        u_new(i,j) = 0.25d0 * ( u(i+1,j) + u(i-1,j) + &
                                u(i,j+1) + u(i,j-1) )
        diff = abs(u_new(i,j) - u(i,j))
        err = max(err, diff)
      end do
    end do

    ! Pointer swap (host-side, zero cost)
    tmp => u
    u   => u_new
    u_new => tmp

    if (err < tol) exit
  end do

  !-----------------------------------------------------------------------
  ! Copy result back
  !-----------------------------------------------------------------------
  !$acc update self(u_buf)

  print *, "Converged in", iter, "iterations with error", err

  !-----------------------------------------------------------------------
  ! Verification
  !-----------------------------------------------------------------------
  maxErr = 0.0d0
  !$acc parallel loop gang vector collapse(2) present(u_buf) reduction(max:maxErr)
  do j = 1, Ny
    do i = 1, Nx
      y = (j - 1)*dy
      x = (i - 1)*dx
      exactVal = sin(pi*x)*exp(-pi*y)
      diff = abs(u(i,j) - exactVal)
      maxErr = max(maxErr, diff)
    end do
  end do

  print *, "Max difference from exact solution = ", maxErr
  print *, "Potential at (x=0.5,y=0.5) ~ ", u((Nx-1)/2,(Ny-1)/2)

  call cpu_time(finish)
  print *, "Elapsed time:", finish - start, "seconds"

  !-----------------------------------------------------------------------
  ! Output
  !-----------------------------------------------------------------------
  open(10,file='laplace_solution.dat',status='replace')
  do j = 1, Ny
    do i = 1, Nx
      y = (j - 1)*dy
      x = (i - 1)*dx
      write(10,'(3E16.8)') x, y, u(i,j)
    end do
  end do
  close(10)

  !-----------------------------------------------------------------------
  ! Cleanup
  !-----------------------------------------------------------------------
  !$acc exit data delete(u_buf, u_new_buf)
  deallocate(u_buf, u_new_buf)

end program laplace_jacobi_gpu_pointer
