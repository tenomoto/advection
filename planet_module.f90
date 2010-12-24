module planet_module
  use kind_module, only: i4b, dp
  use math_module, only: pi=>math_pi
  implicit none

! Geophysical constants
  real(kind=dp), public :: &
    planet_radius = 6.371d6, &
    day_in_sec = 86400.0_dp, angular_velocity

  public :: planet_init

contains

  subroutine planet_init()
    implicit none

    namelist /planet/ planet_radius, day_in_sec

    read(unit=5, nml=planet)
    write(unit=6, nml=planet)
    angular_velocity= 2.0_dp*pi/day_in_sec

  end subroutine planet_init

end module planet_module
