!! Ritchie (1987) advection test
program advection

  use grid_module, only: grid_init
  use eulerian_module, only: eulerian_init, eulerian_timeint, eulerian_clean

  implicit none

! initialization

  call grid_init()
  call eulerian_init()

! time evolution

  call eulerian_timeint()

! clean up

  call eulerian_clean()

end program advection
