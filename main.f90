!! Ritchie (1987) advection test
program advection

  use planet_module, only: planet_init
  use grid_module, only: grid_init, grid_clean
  use time_module, only: time_init, model
  use eulerian_module, only: eulerian_init, eulerian_timeint, eulerian_clean
  use semilag_module, only: semilag_init, semilag_timeint, semilag_clean
  use nisl_module, only: nisl_init, nisl_timeint, nisl_clean

  implicit none

  call planet_init()
  call grid_init()
  call time_init()
  select case(model)
    case("euler ")
      call eulerian_init()
      call eulerian_timeint()
      call eulerian_clean()
    case("slag  ")
      call semilag_init()
      call semilag_timeint()
      call semilag_clean()
    case("nisl  ")
      call nisl_init()
      call nisl_timeint()
      call nisl_clean()
    case default
      print *, "No matching model for", model
  end select
  call grid_clean()

end program advection
