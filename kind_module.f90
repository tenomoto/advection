module kind_module
  implicit none

  integer, parameter, public:: &
    i8b = selected_int_kind(15), &
    i4b = selected_int_kind(9), & 
    sp = kind(1.0), &
    dp = kind(1.0d0)

end module kind_module

