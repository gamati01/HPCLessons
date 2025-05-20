program laplace2d_sor_redblack
  implicit none

  ! --------------------------------------------------------------------------
  ! PARAMETERS
  ! --------------------------------------------------------------------------
  integer, parameter :: Nx = 1024+2         ! Number of grid points in x
  integer, parameter :: Ny = 1024+2         ! Number of grid points in y
  double precision, parameter :: tol = 1.0d-9
  double precision, parameter :: w   = 1.99d0  ! Over-relaxation factor (1 < w < 2)
  integer, parameter :: maxIter = 100000

  ! --------------------------------------------------------------------------
  ! VARIABLES
  ! --------------------------------------------------------------------------
  double precision, dimension(Nx, Ny) :: u
  double precision :: dx, dy, x, y
  double precision :: err, diff, uNew, uOld, start, finish, exactval, maxerr, pi
  integer :: i, j, iter, is_red

  call cpu_time(start)
  
  ! Initialize OpenACC
  !$acc init device_num(0)

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

  ! Copy data to device
  !$acc enter data copyin(u)

  ! --------------------------------------------------------------------------
  ! RED-BLACK SOR ITERATION - GPU ACCELERATED
  ! --------------------------------------------------------------------------
  do iter = 1, maxIter
    err = 0.0d0

    ! Phase 1: Update RED points (where i+j is even)
    !$acc parallel loop gang vector collapse(2) present(u) private(uOld, uNew, diff, is_red) reduction(max:err)
    do j = 2, Ny-1
      do i = 2, Nx-1
        is_red = mod(i+j,2)
        if (is_red == 0) then  ! Only process red points in this phase
          uOld = u(i,j)
          uNew = 0.25d0*( u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1) )
          diff = uNew - uOld
          u(i,j) = uOld + w*diff
          if (abs(diff) > err) err = abs(diff)
        end if
      end do
    end do
    !$acc end parallel loop

    ! Phase 2: Update BLACK points (where i+j is odd)
    !$acc parallel loop gang vector collapse(2) present(u) private(uOld, uNew, diff, is_red) reduction(max:err)
    do j = 2, Ny-1
      do i = 2, Nx-1
        is_red = mod(i+j,2)
        if (is_red == 1) then  ! Only process black points in this phase
          uOld = u(i,j)
          uNew = 0.25d0*( u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1) )
          diff = uNew - uOld
          u(i,j) = uOld + w*diff
          if (abs(diff) > err) err = abs(diff)
        end if
      end do
    end do
    !$acc end parallel loop

    ! Check for convergence
    if (err < tol) then
      print *, 'Converged after ', iter, ' iterations, err=', err
      exit
    end if
  end do

  ! Bring results back to host for analysis
  !$acc update host(u)

  !-----------------------------------------------------------------------
  ! Compare numerical solution to exact solution: sin(pi*x)*exp(-pi*y)
  !-----------------------------------------------------------------------
  maxErr = 0.0d0
  !$acc parallel loop collapse(2) present(u) reduction(max:maxErr) private(x, y, exactVal, diff)
  do j = 1, Ny
    do i = 1, Nx
      y = (j - 1)*dy
      x = (i - 1)*dx
      exactVal = sin(pi*x)*exp(-pi*y)
      diff = dabs(u(i,j) - exactVal)
      if (diff > maxErr) maxErr = diff
    end do
  end do
  !$acc end parallel loop

  print *, "Max difference from exact solution = ", maxErr
  print *, "Potential at (x=0.5,y=0.5) ~ ", u((Nx-1)/2, (Ny-1)/2)

  call cpu_time(finish)
  print '("Time = ",f10.4," seconds.")',finish-start

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

  !$acc exit data delete(u)
  !$acc shutdown

end program laplace2d_sor_redblack
