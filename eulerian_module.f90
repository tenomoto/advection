module eulerian_module

  use constant_module, only: i4b, dp, a=>planet_radius, hour_in_sec
  use parameter_module, only: nlon, nlat, ntrunc, nstep, hstep, deltat
!  use grid_module, only: gu, gv, su, sv, sphi, sphi_old, lon, lat
  use grid_module, only: gu, gv, gphi, sphi, sphi_old, lon, lat
  use legendre_transform_module, only: legendre_analysis, legendre_synthesis, &
                                       legendre_synthesis_dlon, legendre_synthesis_dlat
  use io_module, only: io_save
  private

  real(kind=dp), parameter, public :: time_filter_param = 0.0_dp
  real(kind=dp), dimension(:,:), allocatable, private :: dgphi

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
    allocate(dgphi(nlon,nlat), sphi1(0:ntrunc,0:ntrunc))

    print *, "step=0 hour=0"
    print *, "Saving step=0"
    call io_save(hfile, 1, gphi, "replace")
    call io_save(hfile, 2, gu, "old")
    call io_save(hfile, 3, gv, "old")
    nsave = 1

!    print *, "step=1/2", " hour=", real(0.5*deltat/hour_in_sec)
    print *, "step=1/2", " t=", real(0.5*deltat)
    call update(0.5_dp*deltat,deltat/2)
!    print *, "step=1", " hour=", real(deltat/hour_in_sec)
    print *, "step=1", " t=", real(deltat)
    call update(deltat,deltat)
    if (hstep==1) then
      call io_save(hfile, 3*nsave+1, gphi, "old")
      call io_save(hfile, 3*nsave+2, gu, "old")
      call io_save(hfile, 3*nsave+3, gv, "old")
      nsave = nsave + 1
    end if

  end subroutine eulerian_init

  subroutine eulerian_clean()

    deallocate(dgphi, sphi1)

  end subroutine eulerian_clean

  subroutine eulerian_timeint()
    implicit none

    integer(kind=i4b) :: i, j

    do i=2, nstep
!      print *, "step=", i, " hour=", real(i*deltat/hour_in_sec)
      print *, "step=", i, " t=", real(i*deltat)
      call update(i*deltat,2*deltat)
      if (mod(i,hstep)==0) then
        print *, "Saving step=", i
        do j=1, nlat
          gu(:,j) = gu(:,j) / cos(lat(j)) ! u = U/cos(lat)
          gv(:,j) = gv(:,j) / cos(lat(j)) ! v = V/cos(lat)
        end do
        call io_save(hfile, 3*nsave+1, gphi, "old")
        call io_save(hfile, 3*nsave+2, gu, "old")
        call io_save(hfile, 3*nsave+3, gv, "old")
        nsave = nsave + 1
      end if
    end do

  end subroutine eulerian_timeint

  subroutine update(t,dt)
    use uv_module, only: uv_nodiv
    implicit none

    real(kind=dp), intent(in) :: t, dt

    integer(kind=i4b) :: i,j, m, n

!    call uv_sbody(lon,lat,gu,gv)
    call uv_nodiv(t,lon,lat,gu,gv)
    do j=1, nlat
      do i=1, nlon
        gu(i,j) = gu(i,j) * cos(lat(j)) ! U = u*cos(lat)
        gv(i,j) = gv(i,j) * cos(lat(j)) ! V = v*cos(lat)
      end do
    end do
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
  
