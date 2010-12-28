module semilag_module

  use kind_module, only: i4b, dp
  use grid_module, only: nlon, nlat, ntrunc, &
    gu, gv, gphi, sphi_old, sphi, latitudes=>lat, lon, coslatr, wind, wgt
  private
  
  integer(kind=i4b), private :: nsave = 0, n = 3
  real(kind=dp), dimension(:,:), allocatable, private :: &
    gphi_old, gphix, gphiy, gphixy, midlon, midlat, deplon, deplat

  character(len=6), dimension(7), parameter, private :: methods = &
    (/"bilin ", "polin2", "linpol", "fd    ", "sph   ", "fdy   ", "spcher"/)
  logical, private :: spectral = .true., fmono = .false., conserve = .false.

  private :: update
  public :: semilag_init, semilag_timeint, semilag_clean

contains

  subroutine semilag_init()
    use time_module, only: hstep, deltat, ifile, hfile
    use planet_module, only: a=>planet_radius
    use interpolate_module, only: interpolate_init
    use io_module, only: io_save
    use legendre_transform_module, only: legendre_synthesis
    implicit none

    integer(kind=i4b) :: i,j

    namelist /semilag/ spectral, fmono, conserve

    read(unit=5, nml=semilag)
    write(unit=6, nml=semilag)

    allocate(gphi_old(nlon,nlat), gphix(nlon,nlat),gphiy(nlon,nlat),gphixy(nlon,nlat), &
             midlon(nlon,nlat),midlat(nlon,nlat),deplon(nlon,nlat),deplat(nlon,nlat))
    call interpolate_init(gphi)

    print *, "Saving initial value"
    call io_save(ifile, 1, gphi, "replace")
    call io_save(ifile, 2, gu*a, "old")
    call io_save(ifile, 3, gv*a, "old")
    call legendre_synthesis(sphi_old,gphi_old)
    gphi = gphi_old
    print *, "umax=", real(maxval(gu)*a), " umin=", real(minval(gu)*a)
    print *, "vmax=", real(maxval(gv)*a), " vmin=", real(minval(gv)*a)
    print *, "step=0 t=0"
    print *, "Saving step=0"
    call io_save(hfile, 1, gphi, "replace")
    call io_save(hfile, 2, gu*a, "old")
    call io_save(hfile, 3, gv*a, "old")
    nsave = 1

    do i=1, nlon
      midlon(i,:) = lon(i)
    end do
    do j=1, nlat
      midlat(:,j) = latitudes(j)
    end do

  end subroutine semilag_init

  subroutine semilag_clean()
    use interpolate_module, only: interpolate_clean
    implicit none

    deallocate(gphi_old,gphix,gphiy,gphixy,midlon,midlat,deplon,deplat)
    call interpolate_clean()

  end subroutine semilag_clean

  subroutine semilag_timeint()
    use time_module, only: nstep, hstep, deltat, hfile
    use legendre_transform_module, only: legendre_synthesis
    use io_module, only: io_save
    implicit none

    integer(kind=i4b) :: i

    do i=1, nstep
      print *, "step=", i, " t=", real(i*deltat)
      call update((i-0.5d0)*deltat,0.5_dp*deltat)
      if (mod(i,hstep)==0) then
        print *, "Saving step=", i
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

  subroutine update(t,dt)
    use math_module, only: &
      pi=>math_pi, pir=>math_pir, pih=>math_pih
    use upstream_module, only: find_points
    use time_module, only: etf, imethod, deltat
    use uv_module, only: uv_sbody, uv_nodiv
    use interpolate_module, only: &
      interpolate_set, interpolate_setd, interpolate_setdx, &
      interpolate_bilinear, interpolate_bicubic, interpolate_polin2, &
      interpolate_linpol, interpolate_setdx, interpolate_spcher
    use legendre_transform_module, only: legendre_analysis, legendre_synthesis, &
        legendre_synthesis_dlon, legendre_synthesis_dlat, legendre_synthesis_dlonlat
    implicit none

    real(kind=dp), intent(in) :: t, dt

    integer(kind=i4b) :: i, j, m
    real(kind=dp) :: eps, dlonr
    real(kind=dp), dimension(nlon) :: gphitmp

    select case(wind)
!      case("sbody ")
!        call uv_sbody(lon,latitudes,gu,gv)
      case("nodiv ")
        call uv_nodiv(t,lon,latitudes,gu,gv)
    end select
    call find_points(gu, gv, dt, midlon, midlat, deplon, deplat)

    if (spectral) then
      call legendre_synthesis(sphi_old,gphi_old)
    end if

    if (conserve) then
      do j=1, nlat
        gphi_old(:,j) = gphi_old(:,j)*wgt(j)*coslatr(j)
      end do
      if ((imethod=="sph   ").or.(imethod=="fdy   ").or.(imethod=="spcher")) then
        call legendre_analysis(gphi_old,sphi_old)
      end if
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
      dlonr = 0.25_dp*nlon*pir
      gphix(1,:) = dlonr * (gphi_old(2,:) - gphi_old(nlon,:))
      gphix(nlon,:) = dlonr * (gphi_old(1,:) - gphi_old(nlon-1,:))
      do i=2, nlon-1
        gphix(i,:) = dlonr*(gphi_old(i+1,:) - gphi_old(i-1,:))
      end do
    end if
    if ((imethod=="fd    ").or.(imethod=="fdy   ")) then
! d/dphi
      eps = pih-latitudes(1)
      gphitmp = cshift(gphi_old(:,1),nlon/2)
      gphiy(:,1) = (gphitmp-gphi_old(:,2))/(pih+eps-latitudes(2))
      gphitmp = cshift(gphix(:,1),nlon/2)
      gphixy(:,1) = (gphitmp-gphix(:,2))/(pih+eps-latitudes(2))
      gphitmp = cshift(gphi_old(:,nlat),nlon/2)
      gphiy(:,nlat) = (gphitmp-gphi_old(:,nlat-1))/(-pih-eps-latitudes(nlat-1))
      gphitmp = cshift(gphix(:,nlat),nlon/2)
      gphixy(:,nlat) = (gphitmp-gphix(:,nlat-1))/(-pih-eps-latitudes(nlat-1))
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
            call interpolate_bilinear(deplon(i,j), deplat(i,j), gphi(i,j))
          case ("fd    ", "sph   ", "fdy   ")
            call interpolate_bicubic(deplon(i,j), deplat(i,j), gphi(i,j), monotonic=fmono)
          case ("polin2")
            call interpolate_polin2(deplon(i,j), deplat(i,j), gphi(i,j), monotonic=fmono)
          case ("linpol")
            call interpolate_linpol(deplon(i,j), deplat(i,j), gphi(i,j), monotonic=fmono)
          case ("spcher")
            call interpolate_spcher(deplon(i,j), deplat(i,j), gphi(i,j), monotonic=fmono)
        end select
      end do
    end do

    if (conserve) then
      do j=1, nlat
        gphi(:,j) = gphi(:,j)/(wgt(j)*coslatr(j))
      end do
    end if

! spectral
    if (spectral) then
      call legendre_analysis(gphi,sphi)
      sphi_old(:,:) = sphi(:,:)
    else
! grid
      gphi_old(:,:) = gphi(:,:)
    end if

  end subroutine update

end module semilag_module
