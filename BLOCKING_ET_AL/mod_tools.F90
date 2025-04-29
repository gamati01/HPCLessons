!
      module precision_module
      integer, parameter :: dp_kind = kind(1.d0)
      integer, parameter :: sp_kind = kind(1.)
#ifdef SINGLEPRECISION
      integer, parameter :: my_kind = sp_kind
#else
      integer, parameter :: my_kind = dp_kind
#endif
      end module precision_module
!
!
      module timing_module
      use precision_module
      contains
      subroutine timing(t)

      real(my_kind), intent(out) :: t
      integer :: time_array(8)

      call date_and_time(values=time_array)
      t = 3600.*time_array(5)+60.*time_array(6)+time_array(7)+time_array(8)/1000.

      end subroutine timing
      end module timing_module

