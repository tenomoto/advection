module time_module
  use kind_module, only: i4b, dp
  implicit none
  private

  integer(kind=i4b), public ::  nstep, hstep
  real(kind=dp), public :: deltat, etf = 0.0d0

  character(len=6), public :: &
    model = "euler ", imethod = "linpol ", imethoduv = "bilin "

  character(len=256), public :: ifile="init.dat", hfile="hist.dat"

  public :: time_init

contains

  subroutine time_init()
    implicit none

    namelist /time/ nstep, hstep, model, deltat, etf, &
      imethod, imethoduv, ifile, hfile

    read(unit=5, nml=time)
    write(unit=6, nml=time)

  end subroutine time_init

end module time_module
