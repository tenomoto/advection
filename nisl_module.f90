module nisl_module

  use kind_module, only: i4b, dp
  use grid_module, only: nlon, nlat, ntrunc, &
    gu, gv, gphi, sphi_old, sphi, longitudes=>lon, latitudes=>lat, coslatr
  private
  
  integer(kind=i4b), private :: nsave = 0
  integer(kind=i4b), allocatable, private :: p(:,:), q(:,:)
  real(kind=dp), dimension(:,:), allocatable, private :: &
    gphi_old, dgphi, dgphim, gphim,&
    midlon, midlat, deplon, deplat, gum, gvm
  complex(kind=dp), dimension(:,:), allocatable, private :: sphi1

  private :: update
  public :: nisl_init, nisl_timeint, nisl_clean

contains

  subroutine nisl_init()
    use time_module, only: hstep, deltat, ifile, hfile
    use planet_module, only: a=>planet_radius
    use interpolate_module, only: interpolate_init
    use legendre_transform_module, only: legendre_synthesis
    use io_module, only: io_save
    implicit none

    integer(kind=i4b) :: i,j

    allocate(sphi1(0:ntrunc,0:ntrunc),gphi_old(nlon,nlat), &
             gphim(nlon,nlat),dgphi(nlon,nlat),dgphim(nlon,nlat), &
             midlon(nlon,nlat),midlat(nlon,nlat), &
             deplon(nlon,nlat),deplat(nlon,nlat), p(nlon,nlat), q(nlon,nlat), &
             gum(nlon,nlat),gvm(nlon,nlat))
    call interpolate_init(gphi)

    print *, "Saving initial value"
    call io_save(ifile, 1, gphi, "replace")
    call io_save(ifile, 2, gu, "old")
    call io_save(ifile, 3, gv, "old")
    call legendre_synthesis(sphi_old,gphi_old)
    gphi = gphi_old
    print *, "step=0 t=0"
    print *, "Saving step=0"
    print *, "umax=", real(maxval(gu)*a), " umin=", real(minval(gu)*a)
    print *, "vmax=", real(maxval(gv)*a), " vmin=", real(minval(gv)*a)
    call io_save(hfile, 1, gphi, "replace")
    call io_save(hfile, 2, gu, "old")
    call io_save(hfile, 3, gv, "old")
    nsave = 1

    do i=1, nlon
      midlon(i,:) = longitudes(i)
    end do
    do j=1, nlat
      midlat(:,j) = latitudes(j)
    end do

    print *, "step=1/2", " t=", real(0.5d0*deltat)
    call update(0.25d0*deltat,0.5*deltat)
    print *, "step=1", " t=", real(deltat)
    call update(0.5d0*deltat,deltat)
    if (hstep==1) then
      call legendre_synthesis(sphi, gphi)
      call io_save(hfile, 3*nsave+1, gphi, "old")
      call io_save(hfile, 3*nsave+2, gu, "old")
      call io_save(hfile, 3*nsave+3, gv, "old")
      nsave = nsave + 1
    end if

  end subroutine nisl_init

  subroutine nisl_clean()
    use interpolate_module, only: interpolate_clean
    implicit none

    deallocate(sphi1,gphi_old,gphim,dgphi,dgphim,gum,gvm, &
      midlon,midlat,deplon,deplat,p,q)
    call interpolate_clean()

  end subroutine nisl_clean

  subroutine nisl_timeint()
    use time_module, only: nstep, hstep, deltat, hfile
    use legendre_transform_module, only: legendre_synthesis
    use io_module, only: io_save
    implicit none

    integer(kind=i4b) :: i

    do i=2, nstep
      print *, "step=", i, " t=", real(i*deltat)
      call update((i-1)*deltat,deltat)
      if (mod(i,hstep)==0) then
        print *, "Saving step=", i
        call legendre_synthesis(sphi, gphi)
        call io_save(hfile, 3*nsave+1, gphi, "old")
        call io_save(hfile, 3*nsave+2, gu, "old")
        call io_save(hfile, 3*nsave+3, gv, "old")
        nsave = nsave + 1
      end if
    end do

  end subroutine nisl_timeint

  subroutine update(t,dt)
    use time_module, only: etf, imethod
    use uv_module, only: uv_nodiv
    use upstream_module, only: find_points
    use legendre_transform_module, only: legendre_analysis, legendre_synthesis, &
        legendre_synthesis_dlon, legendre_synthesis_dlat, legendre_synthesis_dlonlat
    use interpolate_module, only: &
      interpolate_set, interpolate_bilinear, interpolate_polin2
    implicit none

    integer(kind=i4b) :: i, j, m
    real(kind=dp), intent(in) :: t, dt

    call uv_nodiv(t,longitudes,latitudes,gu,gv)
    call find_points(gu, gv, dt, midlon, midlat, deplon, deplat)
    call calc_niuv(dt)

    call legendre_synthesis(sphi_old, gphi_old)
    do j=1, nlat
      do i=1, nlon
        gphi(i,j) = gphi_old(p(i,j),q(i,j))
      end do
    end do

    ! dF/dlon
    call legendre_synthesis_dlon(sphi, dgphi)
    call interpolate_set(dgphi)
    do j=1, nlat
      do i=1, nlon
        if (imethod=="polin2") then
          call interpolate_polin2(midlon(i,j), midlat(i,j), dgphim(i,j))
        else
          call interpolate_bilinear(midlon(i,j), midlat(i,j), dgphim(i,j))
        end if
      end do
      gphim(:,j) = gum(:,j)*coslatr(j)*dgphim(:,j) ! gum: -u'
    end do

    ! cos(lat)dF/dlat
    call legendre_synthesis_dlat(sphi, dgphi) 
    call interpolate_set(dgphi)
    do j=1, nlat
      do i=1, nlon
        if (imethod=="polin2") then
          call interpolate_polin2(midlon(i,j), midlat(i,j), dgphim(i,j))
        else
          call interpolate_bilinear(midlon(i,j), midlat(i,j), dgphim(i,j))
        end if
      end do
      gphim(:,j) = gphim(:,j) + gvm(:,j)*coslatr(j)*dgphim(:,j) ! gvm: -v'
    end do

    gphi = gphi + dt*gphim

! time filter
    call legendre_analysis(gphi, sphi1)
    do m=0, ntrunc
      sphi_old(m:ntrunc,m) = sphi(m:ntrunc,m) + &
        etf * (sphi_old(m:ntrunc,m)-2.0_dp*sphi(m:ntrunc,m)+sphi1(m:ntrunc,m))
      sphi(m:ntrunc,m) = sphi1(m:ntrunc,m)
    end do

  end subroutine update

  subroutine calc_niuv(dt)
    use math_module, only: pir=>math_pir, pi2=>math_pi2
    use time_module, only: imethoduv
    use sphere_module, only: xyz2uv, lonlat2xyz
    use interpolate_module, only: &
       interpolate_setuv, interpolate_bilinearuv, interpolate_polin2uv
    implicit none

    real(kind=dp), intent(in) :: dt

    integer(kind=i4b) :: i,j,ii
    real(kind=dp) :: xg, yg, zg, xr, yr, zr, xm, ym, zm, xd, yd, zd, &
      lon, lat, lonr, latr, lonm, latm, u, v, dlonr, b

    dlonr = 0.5_dp*nlon*pir
    call interpolate_setuv(gu,gv)
    do j=1, nlat
      lat = latitudes(j)
      do i=1, nlon
! find grid points near departure points
        p(i,j) = anint(deplon(i,j)*dlonr+1.0_dp)
        if (p(i,j)>nlon) then
          p(i,j) = p(i,j)-nlon
        end if
! lat = (J+1-2j)pi/(2J+1)
        q(i,j) = anint(0.5_dp*(nlat+1-(2.0_dp*nlat+1.0_dp)*deplat(i,j)*pir))
        lonr = longitudes(p(i,j))
        latr = latitudes(q(i,j))  
        call lonlat2xyz(lonr,latr,xr,yr,zr)
! arrival points
        lon = longitudes(i)
        call lonlat2xyz(lon,lat,xg,yg,zg)
! calculate midpoints between integer departure points and arrival points
        b = 1.0_dp/sqrt(2.0_dp*(1.0_dp+(xg*xr+yg*yr+zg*zr)))
        xm = b*(xg + xr)
        ym = b*(yg + yr)
        zm = b*(zg + zr)
        midlon(i,j) = modulo(atan2(ym,xm)+pi2,pi2)
        midlat(i,j) = asin(zm)
!       print *, real(lon*180/pi), real(midlon(i,j)*180/pi), real(lonr*180/pi), &
!         real(lat*180/pi), real(midlat(i,j)*180/pi), real(latr*180/pi)
! calculate integer velocities at midpoints
        xd = (xg-xr)/dt
        yd = (yg-yr)/dt
        zd = (zg-zr)/dt
        call xyz2uv(xd, yd, zd, midlon(i,j), midlat(i,j), u, v)
!       print *, real(lonr*180/pi), real(latr*180/pi), real(u*a), real(v*a)
!       print *, real(zg), real(zr), real(zd), real(v*a)
        gum(i,j) = u
        gvm(i,j) = v
! calculate velocity at midpoints and -residual velocities
        if (imethoduv=="polin2") then
          call interpolate_polin2uv(midlon(i,j), midlat(i,j), u, v)
        else
          call interpolate_bilinearuv(midlon(i,j), midlat(i,j), u, v)
        end if
        gum(i,j) = gum(i,j) - u
        gvm(i,j) = gvm(i,j) - v
      end do
    end do
        
  end subroutine  calc_niuv

end module nisl_module
