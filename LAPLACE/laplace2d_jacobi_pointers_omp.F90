program laplace_jacobi_pointers
  use omp_lib
  implicit none

  !-----------------------------------------------------------------------
  ! Parameters
  !-----------------------------------------------------------------------
  integer, parameter :: Nx = 512+1       ! # of grid points in x
  integer, parameter :: Ny = 512+1       ! # of grid points in y
  integer, parameter :: maxIter = 500000 ! Max iterations (Jacobi is slow)
  double precision, parameter :: tol = 1.0d-9  ! Convergence tolerance

  !-----------------------------------------------------------------------
  ! Variables
  !-----------------------------------------------------------------------
  double precision, dimension(:,:), pointer :: u, u_new, temp
  double precision :: dx, dy, x, y, pi
  double precision :: err, diff, exactVal, maxErr, start, finish
  double precision :: u_new_ij, u_ij
  integer :: i, j, iter

  ! Allocate memory
  allocate(u(Nx,Ny), u_new(Nx,Ny))
  
  ! Start timer
  start = omp_get_wtime()  ! <-- Replace cpu_time(start)

  ! Compute pi and uniform grid spacing
  pi = acos(-1.0d0)
  dx = 1.0d0 / dble(Nx - 1)
  dy = 1.0d0 / dble(Ny - 1)

  ! Initialize u and u_new to 0 (or with boundary conditions as needed)
  u = 0.0d0
  u_new = 0.0d0

  !-----------------------------------------------------------------------
  ! Apply boundary conditions (parallelized)
  !-----------------------------------------------------------------------
  ! y=0 => u(x,0) = sin(pi*x)
  !$OMP PARALLEL DO PRIVATE(i,x) SHARED(u, u_new) SCHEDULE(static)
  do i = 1, Nx
    x = (i - 1)*dx
    u(i,1)     = sin(pi*x)
    u_new(i,1) = u(i,1)
  end do
  !$OMP END PARALLEL DO

  ! y=1 => u(x,1) = sin(pi*x)*exp(-pi)
  !$OMP PARALLEL DO PRIVATE(i,x) SHARED(u, u_new) SCHEDULE(static)
  do i = 1, Nx
    x = (i - 1)*dx
    u(i,Ny)     = sin(pi*x)*exp(-pi)
    u_new(i,Ny) = u(i,Ny)
  end do
  !$OMP END PARALLEL DO

  ! Jacobi iteration (parallelized)
  do iter = 1, maxIter
    err = 0.0d0
    !$OMP PARALLEL DO PRIVATE(i,j,diff) SHARED(u, u_new) REDUCTION(max:err) SCHEDULE(static)
    do j = 2, Ny-1
      do i = 2, Nx-1
        u_new(i,j) = 0.25d0 * (u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1))
        diff = abs(u_new(i,j) - u(i,j))
        err = max(err, diff)
      end do
    end do
    !$OMP END PARALLEL DO

    ! Swap pointers instead of copying arrays
    temp => u
    u => u_new
    u_new => temp

    ! Convergence check
    if (err < tol) exit
  end do

  print *, "Converged in", iter, "iterations with error", err

  !-----------------------------------------------------------------------
  ! Compare numerical solution to exact solution (parallelized)
  !-----------------------------------------------------------------------
  maxErr = 0.0d0
  !$OMP PARALLEL DO PRIVATE(i, j, x, y, exactVal, diff) SHARED(u) REDUCTION(max:maxErr)
  do j = 1, Ny
    do i = 1, Nx
      y = (j - 1)*dy
      x = (i - 1)*dx
      exactVal = sin(pi*x)*exp(-pi*y)
      diff = dabs(u(i,j) - exactVal)
      maxErr = max(maxErr, diff)
    end do
  end do
  !$OMP END PARALLEL DO

  print *, "Max difference from exact solution = ", maxErr
  print *, "Potential at (x=0.5,y=0.5) ~ ", u((Nx-1)/2, (Ny-1)/2)

  ! Stop timer
  finish = omp_get_wtime()  ! <-- Replace cpu_time(finish)
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

  ! Deallocate memory
  deallocate(u, u_new)

end program laplace_jacobi_pointers
