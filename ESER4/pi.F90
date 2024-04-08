!--------------------------------------------
!  Exercise: Pi                               
!                                             
!  Compute the value of PI using the integral 
!  pi = sqrt(6* sum 1/(x*x))                        
!                                             
!---------------------------------------------

program pigreco
    implicit none

!    integer(selected_int_kind(18)) :: i, intervals 
    integer :: i, elements 
    real :: sum1, sum2
    real :: pi

    real(kind(1.d0)), parameter :: PI25DT = acos(-1.d0)

    real :: time1, time2
    call cpu_time(time1)

    write(6,*) "Elements?"
    read(5,*) elements
    write(6,*) "Number of elements: ", elements
    sum1=0.0
    sum2=0.0

! forward sum...    
    do i=1,elements
        sum1=sum1+(1.0/(float(i)*float(i)))
    end do
!
! backward sum
    do i=elements, 1, -1
        sum2=sum2+(1.0/(float(i)*float(i)))
    end do
!
    call cpu_time(time2)

    write(6,*) ' The True PI =', PI25DT
    write(6,*) ' Computed PI (forward) =', sqrt(6.0*sum1), sqrt(6.0*sum1)/PI25DT
    write(6,*) ' Computed PI (backward)=', sqrt(6.0*sum2), sqrt(6.0*sum2)/PI25DT
!    PRINT '(a13,2x,f30.25)',' The True PI =', PI25DT
    PRINT *, 'Elapsed time ', time2-time1 ,' s' 

end program pigreco

