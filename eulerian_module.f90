module eulerian_module

  use constant_module, only: i4b, dp, a=>planet_radius, hour_in_sec
  use parameter_module, only: nlon, nlat, ntrunc, nstep, hstep, deltat
  use grid_module, only: gu, gv, su, sv, sphi, sphi_old, lat
  use legendre_transform_module, only: legendre_analysis, legendre_synthesis, &
                                       legendre_synthesis_dlon, legendre_synthesis_dlat
  use io_module, only: save_data
  private

  real(kind=dp), parameter, public :: time_filter_param = 0.0_dp
  real(kind=dp), dimension(:,:), allocatable, private :: gphi, dgphi

  integer(kind=i4b), private :: nsave = 0
  real(kind=dp), dimension(:), allocatable, private :: cos2
  complex(kind=dp), dimension(:,:), allocatable, private :: sphi1
  character(len=*), parameter, private :: hfile = "history.dat"

  private :: update
  public :: eulerian_init, eulerian_timeint, eulerian_clean

contains

  subroutine eulerian_init()
    implicit none

    integer(kind=i4b) :: j

    allocate(cos2(nlat))
    do j=1, nlat
      cos2(j) = cos(lat(j))*cos(lat(j))
    end do
    allocate(gphi(nlon,nlat), &
      dgphi(nlon,nlat), sphi1(0:ntrunc,0:ntrunc))

    print *, "step=0 hour=0"
    call legendre_synthesis(sphi, gphi)
    call legendre_synthesis(su, gu)
    call legendre_synthesis(sv, gv)
    call save_data(hfile, 1, gphi, "replace")
    call save_data(hfile, 2, gu, "old")
    call save_data(hfile, 3, gv, "old")
    nsave = 1

    print *, "step=1/2", " hour=", real(0.5*deltat/hour_in_sec)
    call update(deltat/2)
    print *, "step=1", " hour=", real(deltat/hour_in_sec)
    call update(deltat)
    if (hstep==1) then
      call save_data(hfile, 3*nsave+1, gphi, "old")
      call save_data(hfile, 3*nsave+2, gu, "old")
      call save_data(hfile, 3*nsave+3, gv, "old")
      nsave = nsave + 1
    end if

  end subroutine eulerian_init

  subroutine eulerian_clean()

    deallocate(gphi, dgphi, sphi1)

  end subroutine eulerian_clean

  subroutine eulerian_timeint()
    implicit none

    integer(kind=i4b) :: i

    do i=2, nstep
      print *, "step=", i, " hour=", real(i*deltat/hour_in_sec)
      call update(2*deltat)
      if (mod(i,hstep)==0) then
        print *, "Saving step=", i
        call save_data(hfile, 3*nsave+1, gphi, "old")
        call save_data(hfile, 3*nsave+2, gu, "old")
        call save_data(hfile, 3*nsave+3, gv, "old")
        nsave = nsave + 1
      end if
    end do

  end subroutine eulerian_timeint

  subroutine update(dt)
    implicit none

    real(kind=dp), intent(in) :: dt

    integer(kind=i4b) :: j, m, n

    call legendre_synthesis(sphi_old, gphi)
! dF/dlon
    call legendre_synthesis_dlon(sphi, dgphi)
    do j=1, nlat
      dgphi(:,j) =  dgphi(:,j)/cos2(j)
    end do
    gphi = gphi - dt*gu*dgphi
! cos(lat)dF/dlat
    call legendre_synthesis_dlat(sphi, dgphi)
    do j=1, nlat
      dgphi(:,j) =  dgphi(:,j)/cos2(j)
    end do
    gphi = gphi - dt*gv*dgphi
    call legendre_analysis(gphi, sphi1)
! update
    do m=0, ntrunc
        sphi_old(m:ntrunc,m) = sphi(m:ntrunc,m) + &
          time_filter_param * (sphi_old(m:ntrunc,m)-2.0_dp*sphi(m:ntrunc,m)+sphi1(m:ntrunc,m))
    end do
    sphi = sphi1
  
  end subroutine update

end module eulerian_module
  
