program laplace_jacobi_correct
  implicit none

  !-----------------------------------------------------------------------
  ! Parameters
  !-----------------------------------------------------------------------
  integer, parameter :: Nx = 512+2       ! # of grid points in x
  integer, parameter :: Ny = 512+2       ! # of grid points in y
  integer, parameter :: maxIter = 500000 ! Max iterations (Jacobi is slow)
  double precision, parameter :: tol = 1.0d-9  ! Convergence tolerance

  !-----------------------------------------------------------------------
  ! Variables
  !-----------------------------------------------------------------------
  double precision, dimension(Nx,Ny) :: u, u_new
  double precision :: dx, dy, x, y, pi
  double precision :: err, diff, exactVal, maxErr, start, finish
  double precision :: u_new_ij, u_ij
  integer :: i, j, iter

  ! Compute pi and uniform grid spacing
  pi = acos(-1.0d0)
  dx = 1.0d0 / dble(Nx - 1)
  dy = 1.0d0 / dble(Ny - 1)

  ! Initialize arrays to zero
  u = 0.0d0
  u_new = 0.0d0

  !-----------------------------------------------------------------------
  ! Apply boundary conditions
  !-----------------------------------------------------------------------
  ! y=0 => u(x,0) = sin(pi*x)
  do i = 1, Nx
    x = (i - 1)*dx
    u(i,1)     = sin(pi*x)
    u_new(i,1) = u(i,1)
  end do

  ! y=1 => u(x,1) = sin(pi*x)*exp(-pi)
  do i = 1, Nx
    x = (i - 1)*dx
    u(i,Ny)     = sin(pi*x)*exp(-pi)
    u_new(i,Ny) = u(i,Ny)
  end do

  ! x=0 => u(0,y)=0, x=1 => u(1,y)=0 (already zero-initialized)

  call cpu_time(start)
  !-----------------------------------------------------------------------
  ! Jacobi iteration
  !-----------------------------------------------------------------------
  do iter = 1, maxIter
    err = 0.0d0

    ! Update only interior points
    do j = 2, Ny - 1
      do i = 2, Nx - 1
        ! Laplace update: average of the four neighbors
        u_ij     = u(i,j)
        u_new_ij = 0.25d0 * (u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1))
        ! Compute difference
        diff = dabs(u_new_ij - u_ij)
        if (diff > err) err = diff
        ! Update
        u_new(i,j) = u_new_ij
      end do
    end do

    ! Check for convergence based on max update
    if (err < tol) then
       print *, "Jacobi converged after ", iter, " iterations; err=", err
       exit
    end if

    ! Copy new interior values back into u
    u(2:Nx-1,2:Ny-1) = u_new(2:Nx-1,2:Ny-1)
  end do
  call cpu_time(finish) 

  print '("Time = ",f10.4," seconds.")',finish-start
                        
  if (iter >= maxIter) then
    print *, "Warning: Jacobi did NOT fully converge after", maxIter, "iterations."
  end if

  !-----------------------------------------------------------------------
  ! Compare numerical solution to exact solution: sin(pi*x)*exp(-pi*y)
  !-----------------------------------------------------------------------
  maxErr = 0.0d0
  do j = 1, Ny
    y = (j - 1)*dy
    do i = 1, Nx
      x = (i - 1)*dx
      exactVal = sin(pi*x)*exp(-pi*y)
      diff = dabs(u(i,j) - exactVal)
      if (diff > maxErr) maxErr = diff
    end do
  end do

  print *, "Max difference from exact solution = ", maxErr
  print *, "Potential at (x=0.5,y=0.5) ~ ", u((Nx-1)/2, (Ny-1)/2)

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

end program laplace_jacobi_correct

