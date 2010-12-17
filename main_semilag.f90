!! Ritchie (1987) advection test
program advection

  use parameter_module, only: set_parameters, nlon, nlat, ntrunc
  use grid_module, only: grid_init
  use semilag_module, only: semilag_init, semilag_timeint, semilag_clean

  implicit none

! initialization

  call set_parameters()
  call grid_init(nlon,nlat,ntrunc)
  call semilag_init()

! time evolution

  call semilag_timeint()

! clean up

  call semilag_clean()

end program advection
