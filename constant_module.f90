module constant_module
	implicit none

!	Type
	integer, parameter, public:: &
		i8b = selected_int_kind(15), &
		i4b = selected_int_kind(9), & 
		i2b = selected_int_kind(4), &
		sp = kind(1.0), &
		dp = kind(1.0d0)

! Mathematical constants
	real(kind=dp), parameter :: &
		pi = 3.1415926535897931d0

! Physical constants
	real(kind=dp), parameter, public :: &
		gas_constant = 287.0d0,   &
		specific_heat_at_constant_pressure = 1004.0d0

! Geophysical constants
	real(kind=dp), parameter, public :: &
		planet_radius = 6.371d6,  &
		day_in_sec = 86400.0d0,   &
		hour_in_sec = 3600.0d0,   &
		angular_velocity = 2.d0*pi/day_in_sec
	
end module constant_module

