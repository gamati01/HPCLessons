! =====================================================================
!     ****** CODICI/mm
!
!     COPYRIGHT
!       (c) 2022 by Giorgio Amati
!     NAME
!       mm
!     DESCRIPTION
!       Matrix Multiplication
!       DO_CONCURRENT parallelism
!     INPUTS
!       none
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!       integer variables used: i,j,k,n
!       real variables used:    randd,a,b,c,mflop
!                               time1, time2
!
!     *****
! =====================================================================

program mm
!
    use precision_module
    use timing_module
!
    implicit none
!
    integer :: i,j,k ! index
    integer :: n                                   ! size of the matrix
!
    real(my_kind), dimension(:,:), allocatable:: a ! matrix (origin)
    real(my_kind), dimension(:,:), allocatable:: b ! matrix (origin)
    real(my_kind), dimension(:,:), allocatable:: c ! matrix (destination)
    real(my_kind), dimension(:,:), allocatable:: d ! matrix (destination)
    real(my_kind):: time1, time2 ! timing
    real(my_kind):: time3, time4 ! timing
    real(my_kind):: mflops ! Mflops
    real(dp_kind):: check
!
    write(6,*) "--------------------------------------"
    write(6,*) " Matrix-Matrix Multiplication         "
    write(6,*) " precision used    ",precision(a(1,1))
    write(6,*) " using do concurrent                        "
    write(6,*) " Which matrix size?                   "
    read(5,*) n
    write(6,*) " Matrix size      =", n
    write(6,*) " Memory size (MB) =", 3*n*n*4/1024/1024 
    write(6,*) "--------------------------------------"
!
    allocate(a(1:n,1:n))
    allocate(b(1:n,1:n))
    allocate(c(1:n,1:n))
    allocate(d(1:n,1:n))
!
    mflops = 2*float(n)*float(n)*float(n)/(1000.0*1000.0)
!
! initialization
    call timing(time1)
    call random_number(a)
    call random_number(b)
    c = 0._my_kind
    d = 0._my_kind
    call timing(time2)
    write(*,*) "initialization", time2-time1
    write(*,*) a(n/2,n/2),b(n/2,n/2),c(n/2,n/2)
!
! main loop
    call system("date       > time.log")
    call timing(time1)
!    
    do concurrent(j=1:n)
       do k = 1, n, 1
          do i = 1, n
             c(i,j) = c(i,j) + a(i,k)*b(k,j)   !&
          enddo
       enddo
    enddo
    call timing(time2)
    call system("date       >> time.log")
!
    write(6,*) "----------------------------"
    write(6,*) "concurrent:time for moltiplication ", time2-time1
    write(6,*) "concurrent:MFLOPS                  ", mflops/(time2-time1)
    write(6,*) "concurrent:check                   ", c(n/2-1,n/2-1)
!
#ifdef VALIDATION
! validation
    call timing(time3)
    do j = 1, n
       do k = 1, n
          do i = 1, n
             d(i,j) = d(i,j) + a(i,k)*b(k,j)
          enddo
       enddo
    enddo
    call timing(time4)
!
    check = 0.0d0
! check
    do j = 1, n
       do i = 1, n
          check = check + (d(i,j)-c(i,j))**2
       enddo
    enddo
!    
    write(6,*) "---------------------------------"
    write(6,*) "CPU:time for moltiplication      ", time4-time3
    write(6,*) "CPU:MFLOPS                       ", mflops/(time4-time3)
    write(6,*) "CPU:check                        ", d(n/2-1,n/2-1)
    write(6,*) "---------------------------------"
    write(6,*) "Difference CPU-GPU               ", check
    write(6,*) "---------------------------------"
#endif
!
end program mm

