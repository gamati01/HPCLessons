! =====================================================================
!     ****** CODICI/mm
!
!     COPYRIGHT
!       (c) 2022 by Giorgio Amati
!     NAME
!       mm
!     DESCRIPTION
!       Matrix Multiplication
!       size of the matrix should be given at standard input
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
    integer :: ii,jj,kk ! index
    integer :: n                                   ! size of the matrix
    integer :: step                                ! blocking size
!
    real(my_kind), dimension(:,:), allocatable:: a ! matrix (origin)
    real(my_kind), dimension(:,:), allocatable:: b ! matrix (origin)
    real(my_kind), dimension(:,:), allocatable:: c ! matrix (destination)
    real(my_kind):: time1, time2 ! timing
    real(my_kind):: time3, time4 ! timing
    real(my_kind):: mflops ! Mflops
    real(my_kind):: temp0, temp1, temp2, temp3     ! exploiting temporal locality
    real(dp_kind):: check
!
    step = 64
    write(6,*) "--------------------------------------"
    write(6,*) " Matrix-Matrix Multiplication         "
    write(6,*) " precision used    ",precision(a(1,1))
    write(6,*) " rel. 4, blocking+unrolling  ", step
    write(6,*) " Which matrix size?                   "
    read(5,*) n
    write(6,*) " Matrix size      =", n
    write(6,*) " Memory size (MB) =", 3*n*n*my_kind/1024/1024 
    write(6,*) "--------------------------------------"
!
    allocate(a(1:n,1:n))
    allocate(b(1:n,1:n))
    allocate(c(1:n,1:n))
!
    mflops = 2*float(n)*float(n)*float(n)/(1000.0*1000.0)
    check = 0.0
!
! initialization
    call timing(time1)
    call random_number(a)
    call random_number(b)
    c = 0._my_kind
    call timing(time2)
    write(*,*) "initialization", time2-time1
    write(*,*) a(n/2,n/2),b(n/2,n/2),c(n/2,n/2)
!
! main loop
    call system("date       > time.log")
    call timing(time1)
!
    do jj = 1, n, step
    do kk = 1, n, step
    do ii = 1, n, step
        do j = jj, jj+step-1, 4
        do k = kk, kk+step-1, 4
        do i = ii, ii+step-1
           temp0 =  a(i,k+0)
           temp1 =  a(i,k+1)
           temp2 =  a(i,k+2)
           temp3 =  a(i,k+3)
           c(i,j+0) = c(i,j+0) + temp0*b(k+0,j+0) + temp1*b(k+1,j+0) + temp2*b(k+2,j+0) + temp3*b(k+3,j+0)
           c(i,j+1) = c(i,j+1) + temp0*b(k+0,j+1) + temp1*b(k+1,j+1) + temp2*b(k+2,j+1) + temp3*b(k+3,j+1)
           c(i,j+2) = c(i,j+2) + temp0*b(k+0,j+2) + temp1*b(k+1,j+2) + temp2*b(k+2,j+2) + temp3*b(k+3,j+2)
           c(i,j+3) = c(i,j+3) + temp0*b(k+0,j+3) + temp1*b(k+1,j+3) + temp2*b(k+2,j+3) + temp3*b(k+3,j+3)
        enddo
        enddo
        enddo
    enddo
    enddo
    enddo
!
    call timing(time2)
    call system("date       >> time.log")
!
    write(6,*) "----------------------------"
    write(6,*) "CPU:time for moltiplication ", time2-time1
    write(6,*) "CPU:MFLOPS                  ", mflops/(time2-time1)
    write(6,*) "CPU:check                   ", c(n/2-1,n/2-1)
!
end program mm

