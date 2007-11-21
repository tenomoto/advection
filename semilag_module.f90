module semilag_module

  use constant_module, only: i4b, dp, hour_in_sec, pi, a=>planet_radius
  use parameter_module, only: nlon, nlat, ntrunc, nstep, hstep, deltat
  use glatwgt_module, only: latitudes=>lat
  use init_module, only: su, sv, sphi_old, sphi
  use legendre_transform_module, only: legendre_analysis, legendre_synthesis, &
        legendre_synthesis_dlon, legendre_synthesis_dlat, legendre_synthesis_dlonlat
  use upstream_module, only: find_points
  use interpolate_module, only: interpolate_init, interpolate_clean, &
                          interpolate_set, interpolate_setd, &
                          interpolate_bilinear, interpolate_bicubic, &
                          interpolate_polin2, interpolate_linpol
  use io_module, only: save_data
  private
  
  real(kind=dp), parameter, public :: time_filter_param = 0.000_dp
  integer(kind=i4b), private :: nsave = 0, n = 3
  real(kind=dp), dimension(:,:), allocatable, private :: &
    gu, gv, gphi_old, gphi, gphi1, gphix, gphiy, gphixy, midlon, midlat, deplon, deplat
  complex(kind=dp), dimension(:,:), allocatable, private :: sphi1
  character(len=*), parameter, private :: hfile = "history.dat"

  private :: update
  public :: semilag_init, semilag_timeint, semilag_clean

contains

  subroutine semilag_init()
    implicit none

    integer(kind=i4b) :: i,j

    allocate(sphi1(0:ntrunc,0:ntrunc),gu(nlon,nlat),gv(nlon,nlat), &
             gphi_old(nlon,nlat),gphi(nlon,nlat),gphi1(nlon,nlat), &
             gphix(nlon,nlat),gphiy(nlon,nlat),gphixy(nlon,nlat), &
             midlon(nlon,nlat),midlat(nlon,nlat),deplon(nlon,nlat),deplat(nlon,nlat))
    call interpolate_init(gphi, n)

    print *, "step=0 hour=0"
    call legendre_synthesis(sphi_old,gphi_old)
    gphi = gphi_old
    call legendre_synthesis(su,gu)
    call legendre_synthesis(sv,gv)
    do j=1, nlat
      gu(:,j) = gu(:,j)/cos(latitudes(j))
      gv(:,j) = gv(:,j)/cos(latitudes(j))
    end do
    print *, "umax=", real(maxval(gu)*a), " umin=", real(minval(gu)*a)
    print *, "vmax=", real(maxval(gv)*a), " vmin=", real(minval(gv)*a)
    call save_data(hfile, 1, gphi, "replace")
    call save_data(hfile, 2, gu, "old")
    call save_data(hfile, 3, gv, "old")
    nsave = 1

    do i=1, nlon
      midlon(i,:) = 2.0_dp*pi/nlon*(i-1)
    end do
    do j=1, nlat
      midlat(:,j) = latitudes(j)
    end do

    print *, "step=1", " hour=", real(deltat/hour_in_sec)
    call find_points(gu, gv, 0.5_dp*deltat, midlon, midlat, deplon, deplat)
    call update()
    if (hstep==1) then
      call legendre_synthesis(sphi, gphi)
      call save_data(hfile, 3*nsave+1, gphi, "old")
      call save_data(hfile, 3*nsave+2, gu, "old")
      call save_data(hfile, 3*nsave+3, gv, "old")
      nsave = nsave + 1
    end if
    call find_points(gu, gv, deltat, midlon, midlat, deplon, deplat)

  end subroutine semilag_init

  subroutine semilag_clean()
    implicit none

    deallocate(sphi1,gu,gv,gphi,gphi_old,gphi1,gphix,gphiy,gphixy,midlon,midlat,deplon,deplat)
    call interpolate_clean()

  end subroutine semilag_clean

  subroutine semilag_timeint()
    implicit none

    integer(kind=i4b) :: i

    do i=2, nstep
      print *, "step=", i, " hour=", real(i*deltat/hour_in_sec)
      call update()
      if (mod(i,hstep)==0) then
        print *, "Saving step=", i
! comment below in grid model
!        call legendre_synthesis(sphi, gphi)
        call save_data(hfile, 3*nsave+1, gphi, "old")
        call save_data(hfile, 3*nsave+2, gu, "old")
        call save_data(hfile, 3*nsave+3, gv, "old")
        nsave = nsave + 1
      end if
    end do

  end subroutine semilag_timeint

  subroutine update()
    implicit none

    integer(kind=i4b) :: i, j, m
    real(kind=dp) :: eps, dlonr
    real(kind=dp), dimension(nlon) :: gphitmp

! comment below in grid model
!    call legendre_synthesis(sphi_old,gphi_old)

! calculate spectral derivatives

!    call legendre_synthesis_dlon(sphi_old,gphix)
!    call legendre_synthesis_dlat(sphi_old,gphiy)
!    call legendre_synthesis_dlonlat(sphi_old,gphixy)
!    do j=1, nlat
!      gphiy(:,j) = gphiy(:,j)/cos(latitudes(j))
!      gphixy(:,j) = gphixy(:,j)/cos(latitudes(j))
!    end do

! calculate fd derivatives

!    dlonr = 0.5_dp*nlon/pi
!    gphix(1,:) = dlonr * (gphi_old(2,:) - gphi_old(nlon,:))
!    gphix(nlon,:) = dlonr * (gphi_old(1,:) - gphi_old(nlon-1,:))
!    do i=2, nlon-1
!      gphix(i,:) = dlonr*(gphi_old(i+1,:) - gphi_old(i-1,:))
!    end do
!    eps = 0.5_dp*pi-latitudes(1)
!    gphitmp = cshift(gphi_old(:,1),nlon/2)
!    gphiy(:,1) = (gphitmp-gphi_old(:,2))/(0.5_dp*pi+eps-latitudes(2))
!    gphitmp = cshift(gphix(:,1),nlon/2)
!    gphixy(:,1) = (gphitmp-gphix(:,2))/(0.5_dp*pi+eps-latitudes(2))
!    gphixy(:,1) = (gphix(:,1)-gphix(:,2))/(latitudes(1)-latitudes(2))
!    gphitmp = cshift(gphi_old(:,nlat),nlon/2)
!    gphiy(:,nlat) = (gphitmp-gphi_old(:,nlat-1))/(-0.5_dp*pi-eps-latitudes(nlat-1))
!    gphitmp = cshift(gphix(:,nlat),nlon/2)
!    gphixy(:,nlat) = (gphitmp-gphix(:,nlat-1))/(-0.5_dp*pi-eps-latitudes(nlat-1))
!    do j=2, nlat-1
!      gphiy(:,j) = (gphi_old(:,j+1)-gphi_old(:,j-1))/(latitudes(j+1)-latitudes(j-1))
!      gphixy(:,j) = (gphix(:,j+1)-gphix(:,j-1))/(latitudes(j+1)-latitudes(j-1))
!    end do 
    call interpolate_set(gphi_old)
!    call interpolate_setd(gphix, gphiy, gphixy)
    do j=1, nlat
      do i=1, nlon
!        call interpolate_bilinear(deplon(i,j), deplat(i,j), gphi1(i,j))
!         call interpolate_bicubic(deplon(i,j), deplat(i,j), gphi1(i,j))
!         call interpolate_bicubic(deplon(i,j), deplat(i,j), gphi1(i,j), monotonic=.true.)
!        call interpolate_polin2(deplon(i,j), deplat(i,j), gphi1(i,j))
!        call interpolate_polin2(deplon(i,j), deplat(i,j), gphi1(i,j), monotonic=.true.)
        call interpolate_linpol(deplon(i,j), deplat(i,j), gphi1(i,j))
!        call interpolate_linpol(deplon(i,j), deplat(i,j), gphi1(i,j), monotonic=.true.)
      end do
    end do

! time filter
! spectral
!    call legendre_analysis(gphi1, sphi1)
!    do m=0, ntrunc
!        sphi_old(m:ntrunc,m) = sphi(m:ntrunc,m) + &
!          time_filter_param * (sphi_old(m:ntrunc,m)-2.0_dp*sphi(m:ntrunc,m)+sphi1(m:ntrunc,m))
!    end do
!    sphi = sphi1
! grid
    gphi_old = gphi + time_filter_param * (gphi_old-2.0_dp*gphi+gphi1)
    gphi = gphi1

  end subroutine update

end module semilag_module
