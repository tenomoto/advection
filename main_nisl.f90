!! Ritchie (1987) advection test
program advection

  use parameter_module, only: set_parameters
  use grid_module, only: grid_init
  use nisl_module, only: nisl_init, nisl_timeint, nisl_clean

  implicit none

! initialization

  call set_parameters()
  call grid_init(nlon,nlat,ntrunc)
  call nisl_init()

! time evolution

  call nisl_timeint()

! clean up

  call nisl_clean()

end program advection
