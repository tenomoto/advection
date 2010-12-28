module time_module
  use kind_module, only: i4b, dp
  implicit none
  private

  integer(kind=i4b), public ::  nstep, hstep
  real(kind=dp), public :: deltat, etf = 0.0_dp, kappa = 0.0_dp

  character(len=6), public :: &
    model = "euler ", imethod = "linpol ", imethoduv = "bilin "

  character(len=256), public :: ifile="init.dat", hfile="hist.dat"

  public :: time_init

contains

  subroutine time_init()
    use grid_module, only: ntrunc
    use planet_module, only: d=>day_in_sec
    implicit none

    real(kind=dp) :: tau

    namelist /time/ nstep, hstep, model, deltat, etf, tau, &
      imethod, imethoduv, ifile, hfile

    read(unit=5, nml=time)
    write(unit=6, nml=time)

    if (tau<=0.0_dp) then
      kappa = 0.0_dp
    else
      kappa = 1.0_dp/(tau*d*(ntrunc*(ntrunc+1.0_dp))**2)
    end if
    print *, "tau=", tau, " kappa=", kappa

  end subroutine time_init

end module time_module
