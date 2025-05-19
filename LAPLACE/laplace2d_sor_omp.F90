program laplace2d_sor_redblack
  use omp_lib
  implicit none

  ! Parameters
  integer, parameter :: Nx = 1024+2
  integer, parameter :: Ny = 1024+2
  double precision, parameter :: tol = 1.0d-9, w = 1.99d0
  integer, parameter :: maxIter = 100000

  ! Variables
  double precision :: u(Nx, Ny), dx, dy, x, y, err, diff, start, finish, pi, maxErr, exactVal
  integer :: i, j, iter, rb

  ! Start timer
  start = omp_get_wtime()

  ! Initialize grid and boundary conditions
  pi = acos(-1.0d0)
  dx = 1.0d0 / dble(Nx - 1)
  dy = 1.0d0 / dble(Ny - 1)
  u = 0.0d0

  ! Boundary conditions
  do i = 1, Nx
    x = (i-1)*dx
    u(i,1)  = sin(pi*x)          ! y=0
    u(i,Ny) = sin(pi*x)*exp(-pi) ! y=1
  end do

  ! Red-Black SOR
  do iter = 1, maxIter
    err = 0.0d0

    ! Phase 1: Update RED cells (i+j even)
    !$OMP PARALLEL DO PRIVATE(i,j,diff) REDUCTION(max:err)
    do j = 2, Ny-1
      do i = 2 + mod(j,2), Nx-1, 2  ! Stride-2 for red cells
        diff = 0.25d0*(u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1)) - u(i,j)
        u(i,j) = u(i,j) + w * diff
        err = max(err, abs(diff))
      end do
    end do
    !$OMP END PARALLEL DO

    ! Phase 2: Update BLACK cells (i+j odd)
    !$OMP PARALLEL DO PRIVATE(i,j,diff) REDUCTION(max:err)
    do j = 2, Ny-1
      do i = 2 + mod(j+1,2), Nx-1, 2  ! Stride-2 for black cells
        diff = 0.25d0*(u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1)) - u(i,j)
        u(i,j) = u(i,j) + w * diff
        err = max(err, abs(diff))
      end do
    end do
    !$OMP END PARALLEL DO

    ! Check convergence
    if (err < tol) then
      print *, "Converged in", iter, "iterations. Error =", err
      exit
    end if
  end do

  ! Verify solution
  maxErr = 0.0d0
  !$OMP PARALLEL DO PRIVATE(i,j,x,y,exactVal,diff) REDUCTION(max:maxErr)
  do j = 1, Ny
    do i = 1, Nx
      x = (i-1)*dx
      y = (j-1)*dy
      exactVal = sin(pi*x)*exp(-pi*y)
      maxErr = max(maxErr, abs(u(i,j) - exactVal))
    end do
  end do
  !$OMP END PARALLEL DO

  print *, "Max difference from exact solution = ", maxErr
  print *, "Potential at (x=0.5,y=0.5) ~ ", u((Nx-1)/2, (Ny-1)/2)

  ! Stop timer
  finish = omp_get_wtime()  ! <-- Replace cpu_time(finish)
  print '("Time = ",f10.4," seconds.")',finish-start

  ! Optional: write the full solution to a file for visualization
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

end program laplace2d_sor_redblack
