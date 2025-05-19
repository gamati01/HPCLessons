program laplace2d_sor
  use omp_lib
  implicit none

  ! --------------------------------------------------------------------------
  ! PARAMETERS
  ! --------------------------------------------------------------------------
  integer, parameter :: Nx = 512+2         ! Number of grid points in x
  integer, parameter :: Ny = 512+2         ! Number of grid points in y
  double precision, parameter :: tol = 1.0d-9
  double precision, parameter :: w   = 1.99d0  ! Over-relaxation factor (1 < w < 2) WARNIG: it is case dependent!
  integer, parameter :: maxIter = 100000

  ! --------------------------------------------------------------------------
  ! VARIABLES
  ! --------------------------------------------------------------------------
  double precision, dimension(Nx, Ny) :: u 
  double precision :: dx, dy, x, y
  double precision :: err, diff, uNew, uOld, start, finish, exactval, maxerr, pi
  integer :: i, j, iter

  ! Start timer
  start = omp_get_wtime()  ! <-- Replace cpu_time(start)

  ! --------------------------------------------------------------------------
  ! COMPUTE GRID SPACING
  ! --------------------------------------------------------------------------
  pi = acos(-1.0d0)
  dx = 1.0d0 / dble(Nx - 1)
  dy = 1.0d0 / dble(Ny - 1)

  ! --------------------------------------------------------------------------
  ! INITIALIZE SOLUTION AND APPLY BOUNDARY CONDITIONS
  ! --------------------------------------------------------------------------
  u = 0.0d0

  ! Along x=0 and x=1 => u=0 (already zero-initialized)
  ! y=0 => u(x,0) = sin(pi*x)
  do i = 0, Nx-1
    x = dble(i) * dx
    u(i+1,1) = sin(pi*x)
  end do

  ! y=1 => u(x,1) = sin(pi*x)*exp(-pi)
  do i = 0, Nx-1
    x = dble(i) * dx
    u(i+1,Ny) = sin(pi*x)*exp(-pi)
  end do

  ! --------------------------------------------------------------------------
  ! SOR ITERATION
  ! --------------------------------------------------------------------------
  do iter = 1, maxIter
    err = 0.0d0

    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,uOld,uNew,diff) SHARED(u) REDUCTION(max:err)
    ! Update interior points (i=2..Nx-1, j=2..Ny-1)
    do j = 2, Ny-1
       do i = 2, Nx-1
          uOld = u(i,j)

          ! Laplace update (5-point stencil, uniform dx=dy)
          uNew = 0.25d0*( u(i+1,j) + u(i-1,j) + &
                          u(i,j+1) + u(i,j-1) )

          ! SOR update
          diff   = uNew - uOld
          u(i,j) = uOld + w*diff

          ! Track maximum update
          err = max(err, abs(diff))
       end do
    end do
    !$OMP END PARALLEL DO

    ! Check for convergence
    if (err < tol) then
       print *, 'Converged after ', iter, ' iterations, err=', err
       exit
    end if

  end do

  if (iter > maxIter) then
     print *, 'Warning: reached maxIter without full convergence.'
  end if

  !-----------------------------------------------------------------------
  ! Compare numerical solution to exact solution: sin(pi*x)*exp(-pi*y)
  !-----------------------------------------------------------------------
  maxErr = 0.0d0
  !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, x, y, exactVal, diff) SHARED(dy, dx, pi, u) REDUCTION(max:maxErr)
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
  print '("Time = ",f10.4," seconds.")',finish-start

  ! Optional: write the full solution to a file for visualization
  open(unit=10, file='laplace_solution.dat', status='replace')
  do j = 1, Ny
     y = (j-1)*dy
     do i = 1, Nx
        x = (i-1)*dx
        write(10,'(3E16.8)') x, y, u(i,j)
     end do
     write(10,*)
  end do
  close(10)

end program laplace2d_sor

