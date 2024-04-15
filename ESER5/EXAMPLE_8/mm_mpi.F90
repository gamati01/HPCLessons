program matrix_matrix_prod
   use mpi
   implicit none
   integer :: n
   integer :: nprocs, myrank, status(MPI_STATUS_SIZE)
   integer :: i, j, k, ierr, jstart, jend
   real*8 :: d
   real*8 :: time1,  time2
   real*8 :: time11, time22
   real*8, dimension(:,:), allocatable   :: a, b, c
   character(len=128) :: command
   character(len=80) :: arg

!mpi stuff (setup)
   call MPI_INIT(ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

!
! reading from standard input the size of the matrix
   if(myrank.eq.0) then
     write(*,*) "MPI version with task = ", nprocs
     call get_command_argument(0,command)
     if (command_argument_count() /= 1) then
       write(0,*) "Usage:", trim(command), "   matrix size"
       call MPI_ABORT(MPI_COMM_WORLD, 911,ierr)
       call MPI_FINALIZE(ierr)
       stop
     else
       call get_command_argument(1,arg)
       read(arg,*) n
     endif
   endif
   call MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!a check
   if (n > 0 ) then
     if (myrank == 0 ) then
       write(*,*) "Matrix size is ", n
     endif
   else
     write(*,*) "Error, matrix size is ", n
     call MPI_FINALIZE(ierr)
     stop
   endif
!
   if(nprocs /= 2) then
     if(myrank == 0) then
       write(*,*) "Error, nprocs=", nprocs, "is not 2!!!!"
       write(*,*) "This (stupid) code works only with 2 task!!!"
     endif
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     call MPI_FINALIZE(ierr)
     stop
   endif
!
!allocation/inizializations
   allocate(a(n,n),stat=ierr)
   if(ierr/=0) STOP "a matrix allocation failed"
   allocate(b(n,n),stat=ierr)
   if(ierr/=0) STOP "b matrix allocation failed"
   allocate(c(n,n),stat=ierr)
   if(ierr/=0) STOP "c matrix allocation failed"
   c = 0.d0

   if(myrank == 0) then
     call random_number(a)
     call random_number(b)
   endif
!
!sending a and b elements
   time1 = MPI_Wtime()
   if(myrank.eq.1) then
     call mpi_recv(a(1,1), n*n, MPI_DOUBLE_PRECISION,0,1,                 &
                             MPI_COMM_WORLD, status,ierr)
     call mpi_recv(b(1,1), n*n, MPI_DOUBLE_PRECISION,0,2,                 &
                             MPI_COMM_WORLD, status,ierr)
   endif
!
   if(myrank.eq.0) then
     call mpi_send(a(1,1), n*n, MPI_DOUBLE_PRECISION,1,1,                   &
                             MPI_COMM_WORLD, ierr)
     call mpi_send(b(1,1), n*n, MPI_DOUBLE_PRECISION,1,2,                   &
                             MPI_COMM_WORLD, ierr)
   endif
   call mpi_barrier(MPI_COMM_WORLD,ierr)
   time11 = MPI_Wtime()
!
   if (myrank == 0) then
     jstart = 1
     jend = n/2
   endif
!
   if (myrank == 1) then
     jstart = n/2+1
     jend = n
   endif
!
   do j=jstart, jend
      do k=1, n
         do i=1, n
            c(i,j) = c(i,j) + a(i,k)*b(k,j)
         end do
      end do
   end do
   time22 = MPI_Wtime()
!
!collecting c elements
   if(myrank == 0) then
     call mpi_recv(c(1,n/2+1), n*n/2, MPI_DOUBLE_PRECISION,1,4,             &
                            MPI_COMM_WORLD, status,ierr)
   endif
!
   if(myrank == 1) then
     call mpi_send(c(1,n/2+1), n*n/2, MPI_DOUBLE_PRECISION,0,4,             &
                            MPI_COMM_WORLD, ierr)
   endif
   call mpi_barrier(MPI_COMM_WORLD,ierr)
   time2 = MPI_Wtime()
! check
   if(myrank == 0) then
     call random_number(d)
     i = int( d*n+1)
     call random_number(d)
     j = int( d*n+1)
     d = 0.d0
     do k=1, n
        d = d + a(i,k)*b(k,j)
     end do
   endif

   if(myrank == 0) then
     write(*,*) "Check on a random element (error):" , abs(d-c(i,j))
     write(*,*) "Elapsed time                      ", time2-time1 ," s"
     write(*,*) "Measured Gflops (no mpi)          ", 2.0*n*n*n/(time22-time11)/1000**3
     write(*,*) "Measured Gflops (with mpi)        ", 2.0*n*n*n/(time2-time1)/1000**3
   endif

   deallocate(a,b,c)

   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   call MPI_FINALIZE(ierr)

   write(6,*) "All done..."

end program matrix_matrix_prod
