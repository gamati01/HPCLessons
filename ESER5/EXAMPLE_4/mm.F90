! =====================================================================
!     ****** CODICI/mm
!
!     COPYRIGHT
!       (c) 2022 by Giorgio Amati
!     NAME
!       mm
!     DESCRIPTION
!       Matrix Multiplication
!       CUDAFortran version
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
    use cudafor

!
#ifdef _OPENMP
    use omp_lib
#endif
!
    implicit none
!
    integer :: i,j,k ! index
    integer :: n                                   ! size of the matrix
!
    real(my_kind), dimension(:,:), allocatable:: a ! matrix (origin)
    real(my_kind), dimension(:,:), allocatable:: b ! matrix (origin)
    real(my_kind), dimension(:,:), allocatable:: c ! matrix (destination)
    real(my_kind), dimension(:,:), device, allocatable:: a_gpu ! matrix (origin)
    real(my_kind), dimension(:,:), device, allocatable:: b_gpu ! matrix (origin)
    real(my_kind), dimension(:,:), device, allocatable:: c_gpu ! matrix (destination)

    real(my_kind):: time1, time2 ! timing
    real(my_kind):: time3, time4 ! timing
    real(my_kind):: mflops ! Mflops
    real(dp_kind):: check
!
    write(6,*) "--------------------------------------"
    write(6,*) " Matrix-Matrix Multiplication         "
    write(6,*) " precision used    ",precision(a(1,1))
    write(6,*) " Cuda fortran (CUF)                   "
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
! copy to matrix defined on the device!
    a_gpu = a
    b_gpu = b
    c_gpu = c
!$cuf kernel  do(2) <<<*,*>>>
    do j = 1, n
       do i = 1, n
          do k = 1, n
             c_gpu(i+0,j) = c_gpu(i+0,j) + a_gpu(i+0,k)*b_gpu(k,j)
          enddo
       enddo
    enddo
!
! copy bak the results
    c = c_gpu
!    
    call timing(time2)
    call system("date       >> time.log")
!
    write(6,*) "----------------------------"
    write(6,*) "GPU:time for moltiplication ", time2-time1
    write(6,*) "GPU:MFLOPS                  ", mflops/(time2-time1)
    write(6,*) "GPU:check                   ", c(n/2-1,n/2-1)
!
end program mm

