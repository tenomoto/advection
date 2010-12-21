module semilag_module

  use constant_module, only: i4b, dp, hour_in_sec, pi, a=>planet_radius
  use grid_module, only: nlon, nlat, ntrunc, &
    gu, gv, gphi, sphi_old, sphi, latitudes=>lat, lon, coslatr
  use time_module, only: nstep, hstep, deltat, imethod
  use legendre_transform_module, only: legendre_analysis, legendre_synthesis, &
        legendre_synthesis_dlon, legendre_synthesis_dlat, legendre_synthesis_dlonlat
  use upstream_module, only: find_points
  use interpolate_module, only: interpolate_init, interpolate_clean, &
                          interpolate_set, interpolate_setd, &
                          interpolate_bilinear, interpolate_bicubic, &
                          interpolate_polin2, interpolate_linpol, &
                          interpolate_setdx, interpolate_spcher
  use io_module, only: io_save
  private
  
  integer(kind=i4b), private :: nsave = 0, n = 3
  real(kind=dp), dimension(:,:), allocatable, private :: &
    gphi_old, gphi1, gphix, gphiy, gphixy, midlon, midlat, deplon, deplat
  complex(kind=dp), dimension(:,:), allocatable, private :: sphi1
  character(len=*), parameter, private :: &
    ifile = "init.dat", hfile = "history.dat"

  character(len=6), dimension(7), parameter, private :: methods = &
    (/"bilin ", "polin2", "linpol", "fd    ", "sph   ", "fdy   ", "spcher"/)
  logical, private :: spectral = .true., fmono = .false.

  private :: update
  public :: semilag_init, semilag_timeint, semilag_clean

contains

  subroutine semilag_init()
    implicit none

    integer(kind=i4b) :: i,j

    namelist /semilag/ spectral, fmono

    read(unit=5, nml=semilag)
    write(unit=6, nml=semilag)

    allocate(sphi1(0:ntrunc,0:ntrunc), &
             gphi_old(nlon,nlat),gphi1(nlon,nlat), &
             gphix(nlon,nlat),gphiy(nlon,nlat),gphixy(nlon,nlat), &
             midlon(nlon,nlat),midlat(nlon,nlat),deplon(nlon,nlat),deplat(nlon,nlat))
    call interpolate_init(gphi)

!    print *, "step=0 hour=0"
    print *, "step=0 t=0"
    print *, "Saving step=0"
    call legendre_synthesis(sphi_old,gphi_old)
    gphi = gphi_old
!    print *, "umax=", real(maxval(gu)*a), " umin=", real(minval(gu)*a)
!    print *, "vmax=", real(maxval(gv)*a), " vmin=", real(minval(gv)*a)
    print *, "umax=", real(maxval(gu)), " umin=", real(minval(gu))
    print *, "vmax=", real(maxval(gv)), " vmin=", real(minval(gv))
    call io_save(ifile, 1, gphi, "replace")
    call io_save(ifile, 2, gu, "old")
    call io_save(ifile, 3, gv, "old")
    call io_save(hfile, 1, gphi, "replace")
    call io_save(hfile, 2, gu, "old")
    call io_save(hfile, 3, gv, "old")
    nsave = 1

    do i=1, nlon
      midlon(i,:) = lon(i)
    end do
    do j=1, nlat
      midlat(:,j) = latitudes(j)
    end do

    call update(0.5d0*deltat)
!    print *, "step=1", " hour=", real(deltat/hour_in_sec)
    print *, "step=1", " hour=", real(deltat)
    if (hstep==1) then
      call legendre_synthesis(sphi, gphi)
      call io_save(hfile, 3*nsave+1, gphi, "old")
      call io_save(hfile, 3*nsave+2, gu, "old")
      call io_save(hfile, 3*nsave+3, gv, "old")
      nsave = nsave + 1
    end if

  end subroutine semilag_init

  subroutine semilag_clean()
    implicit none

    deallocate(sphi1,gphi_old,gphi1,gphix,gphiy,gphixy,midlon,midlat,deplon,deplat)
    call interpolate_clean()

  end subroutine semilag_clean

  subroutine semilag_timeint()
    implicit none

    integer(kind=i4b) :: i

    do i=2, nstep+1
!      print *, "step=", i, " hour=", real(i*deltat/hour_in_sec)
      print *, "step=", i, " t=", real(i*deltat)
      call update((i-1)*deltat)
      if (mod(i-1,hstep)==0) then
        print *, "Saving step=", i-1
        if (spectral) then
          call legendre_synthesis(sphi, gphi)
        end if
        call io_save(hfile, 3*nsave+1, gphi, "old")
        call io_save(hfile, 3*nsave+2, gu, "old")
        call io_save(hfile, 3*nsave+3, gv, "old")
        nsave = nsave + 1
      end if
    end do

  end subroutine semilag_timeint

  subroutine update(t)
    use time_module, only: etf
    use uv_module, only: uv_nodiv
    implicit none

    real(kind=dp), intent(in) :: t

    integer(kind=i4b) :: i, j, m
    real(kind=dp) :: eps, dlonr
    real(kind=dp), dimension(nlon) :: gphitmp

    call uv_nodiv(t,lon,latitudes,gu,gv)
    call find_points(gu, gv, deltat, midlon, midlat, deplon, deplat)

    if (spectral) then
      call legendre_synthesis(sphi_old,gphi_old)
    end if

! calculate spectral derivatives

    if ((imethod=="sph   ").or.(imethod=="fdy   ").or.(imethod=="spcher")) then
      call legendre_synthesis_dlon(sphi_old,gphix)
    end if
    if (imethod=="sph   ") then
      call legendre_synthesis_dlat(sphi_old,gphiy)
      call legendre_synthesis_dlonlat(sphi_old,gphixy)
      do j=1, nlat
        gphiy(:,j) = gphiy(:,j)*coslatr(j)
        gphixy(:,j) = gphixy(:,j)*coslatr(j)
      end do
    end if

! calculate fd derivatives

    if (imethod=="fd    ") then
! d/dlon
      dlonr = 0.25_dp*nlon/pi
      gphix(1,:) = dlonr * (gphi_old(2,:) - gphi_old(nlon,:))
      gphix(nlon,:) = dlonr * (gphi_old(1,:) - gphi_old(nlon-1,:))
      do i=2, nlon-1
        gphix(i,:) = dlonr*(gphi_old(i+1,:) - gphi_old(i-1,:))
      end do
    end if
    if ((imethod=="fd    ").or.(imethod=="fdy   ")) then
! d/dphi
      eps = 0.5_dp*pi-latitudes(1)
      gphitmp = cshift(gphi_old(:,1),nlon/2)
      gphiy(:,1) = (gphitmp-gphi_old(:,2))/(0.5_dp*pi+eps-latitudes(2))
      gphitmp = cshift(gphix(:,1),nlon/2)
      gphixy(:,1) = (gphitmp-gphix(:,2))/(0.5_dp*pi+eps-latitudes(2))
      gphitmp = cshift(gphi_old(:,nlat),nlon/2)
      gphiy(:,nlat) = (gphitmp-gphi_old(:,nlat-1))/(-0.5_dp*pi-eps-latitudes(nlat-1))
      gphitmp = cshift(gphix(:,nlat),nlon/2)
      gphixy(:,nlat) = (gphitmp-gphix(:,nlat-1))/(-0.5_dp*pi-eps-latitudes(nlat-1))
      do j=2, nlat-1
        gphiy(:,j) = (gphi_old(:,j+1)-gphi_old(:,j-1))/(latitudes(j+1)-latitudes(j-1))
        gphixy(:,j) = (gphix(:,j+1)-gphix(:,j-1))/(latitudes(j+1)-latitudes(j-1))
      end do 
    end if

! set grids
    call interpolate_set(gphi_old)
    if (imethod=="spcher") then
      call interpolate_setdx(gphix)
    end if
    if ((imethod=="fd    ").or.(imethod=="sph   ").or. &
        (imethod=="fdy   ")) then
      call interpolate_setd(gphix, gphiy, gphixy)
    end if
    do j=1, nlat
      do i=1, nlon
        select case (imethod)
          case ("bilin ") ! monotonicity is guranteed
            call interpolate_bilinear(deplon(i,j), deplat(i,j), gphi1(i,j))
          case ("fd    ", "sph   ", "fdy   ")
            call interpolate_bicubic(deplon(i,j), deplat(i,j), gphi1(i,j), monotonic=fmono)
          case ("polin2")
            call interpolate_polin2(deplon(i,j), deplat(i,j), gphi1(i,j), monotonic=fmono)
          case ("linpol")
            call interpolate_linpol(deplon(i,j), deplat(i,j), gphi1(i,j), monotonic=fmono)
          case ("spcher")
            call interpolate_spcher(deplon(i,j), deplat(i,j), gphi1(i,j), monotonic=fmono)
        end select
      end do
    end do

! time filter
! spectral
    if (spectral) then
      call legendre_analysis(gphi1, sphi1)
      do m=0, ntrunc
          sphi_old(m:ntrunc,m) = sphi(m:ntrunc,m) + &
            etf * (sphi_old(m:ntrunc,m)-2.0_dp*sphi(m:ntrunc,m)+sphi1(m:ntrunc,m))
      end do
      sphi = sphi1
    else
! grid
      gphi_old = gphi + etf * (gphi_old-2.0_dp*gphi+gphi1)
      gphi = gphi1
    end if

  end subroutine update

end module semilag_module
