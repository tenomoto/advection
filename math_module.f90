module math_module
  use kind_module, only: i4b, dp
  implicit none

! Mathematical constants
  real(kind=dp), parameter, public :: &
    math_pi = 3.14159265358979323846264338327950288_dp, &
    math_pir = 1.0_dp/math_pi, &
    math_pi2 = 2.0_dp*math_pi, math_pih = 0.5_dp*math_pi, &
    math_deg2rad = math_pi/180.0_dp

end module math_module
