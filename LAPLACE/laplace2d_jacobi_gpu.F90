program laplace_jacobi_gpu
  implicit none

  !-----------------------------------------------------------------------
  ! Parameters
  !-----------------------------------------------------------------------
  integer, parameter :: Nx = 512+2       ! # of grid points in x
  integer, parameter :: Ny = 512+2       ! # of grid points in y
  integer, parameter :: maxIter = 500000 ! Max iterations (Jacobi is slow)
  double precision, parameter :: tol = 1.0d-8  ! Convergence tolerance

  !-----------------------------------------------------------------------
  ! Variables
  !-----------------------------------------------------------------------
  double precision, dimension(:,:), allocatable :: u, u_new
  double precision :: dx, dy, x, y, pi
  double precision :: err, diff, exactVal, maxErr, start, finish
  integer :: i, j, iter
  logical :: swap_flag  ! To alternate between u and u_new

  ! Allocate memory
  allocate(u(Nx,Ny), u_new(Nx,Ny))

  ! Start timer
  call cpu_time(start)

  ! Compute pi and uniform grid spacing
  pi = acos(-1.0d0)
  dx = 1.0d0 / dble(Nx - 1)
  dy = 1.0d0 / dble(Ny - 1)

  ! Initialize u and u_new to 0
  !$acc enter data create(u, u_new)
  !$acc kernels
  u = 0.0d0
  u_new = 0.0d0
  !$acc end kernels

  !-----------------------------------------------------------------------
  ! Apply boundary conditions (always applied to u)
  !-----------------------------------------------------------------------
  !$acc parallel loop
  do i = 1, Nx
    x = (i - 1)*dx
    u(i,1)     = sin(pi*x)
    u_new(i,1) = u(i,1)
  end do

  !$acc parallel loop
  do i = 1, Nx
    x = (i - 1)*dx
    u(i,Ny)     = sin(pi*x)*exp(-pi)
    u_new(i,Ny) = u(i,Ny)
  end do

  ! Jacobi iteration
  swap_flag = .false.  ! Start with u as the "current" solution
  do iter = 1, maxIter
    err = 0.0d0

    if (swap_flag) then
      !$acc parallel loop gang vector collapse(2) reduction(max:err)
      do j = 2, Ny-1
        do i = 2, Nx-1
          u(i,j) = 0.25d0 * (u_new(i+1,j) + u_new(i-1,j) + u_new(i,j+1) + u_new(i,j-1))
          diff = abs(u(i,j) - u_new(i,j))
          err = max(err, diff)
        end do
      end do
    else
      !$acc parallel loop gang vector collapse(2) reduction(max:err)
      do j = 2, Ny-1
        do i = 2, Nx-1
          u_new(i,j) = 0.25d0 * (u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1))
          diff = abs(u_new(i,j) - u(i,j))
          err = max(err, diff)
        end do
      end do
    end if

    swap_flag = .not. swap_flag  ! Toggle the flag for next iteration

    ! Convergence check
    if (err < tol) exit
  end do

  ! Ensure final solution is in u (not u_new)
  if (swap_flag) then
    !$acc parallel loop gang vector collapse(2)
    do j = 2, Ny-1
      do i = 2, Nx-1
        u(i,j) = u_new(i,j)
      end do
    end do
  end if

  !$acc update self(u)
  print *, "Converged in", iter, "iterations with error", err

  !-----------------------------------------------------------------------
  ! Compare numerical solution to exact solution
  !-----------------------------------------------------------------------
  maxErr = 0.0d0
  !$acc parallel loop gang vector collapse(2) reduction(max:maxErr)
  do j = 1, Ny
    do i = 1, Nx
      y = (j - 1)*dy
      x = (i - 1)*dx
      exactVal = sin(pi*x)*exp(-pi*y)
      diff = dabs(u(i,j) - exactVal)
      if (diff > maxErr) maxErr = diff
    end do
  end do

  print *, "Max difference from exact solution = ", maxErr
  print *, "Potential at (x=0.5,y=0.5) ~ ", u((Nx-1)/2, (Ny-1)/2)

  ! Stop timer
  call cpu_time(finish)
  print *, "Elapsed time:", finish - start, "seconds"

  !-----------------------------------------------------------------------
  ! Optionally write full solution to file for visualization
  !-----------------------------------------------------------------------
  open(unit=10, file='laplace_solution.dat', status='replace')
  do j = 1, Ny
    y = (j - 1)*dy
    do i = 1, Nx
      x = (i - 1)*dx
      write(10,'(3E16.8)') x, y, u(i,j)
    end do
    write(10,*)
  end do
  close(10)

  ! Cleanup
  !$acc exit data delete(u, u_new)
  deallocate(u, u_new)

end program laplace_jacobi_gpu
