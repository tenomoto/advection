module nisl_module

	use constant_module, only: i4b, dp, hour_in_sec, pi, a=>planet_radius
	use parameter_module, only: nlon, nlat, ntrunc, nstep, hstep, deltat
	use glatwgt_module, only: latitudes=>lat
	use init_module, only: su, sv, sphi_old, sphi
	use legendre_transform_module, only: legendre_analysis, legendre_synthesis, &
	      legendre_synthesis_dlon, legendre_synthesis_dlat, legendre_synthesis_dlonlat
	use upstream_module, only: find_points
  use interpolate_module, only: interpolate_init, interpolate_clean, &
                                interpolate_set, interpolate_bilinear, &
                                interpolate_setuv, interpolate_bilinearuv
	use io_module, only: save_data
	use sphere_module, only: xyz2uv, lonlat2xyz, xy2lon, lat2j
	private
	
	real(kind=dp), parameter, public :: time_filter_param = 0.000_dp
	integer(kind=i4b), private :: nsave = 0
	integer(kind=i4b), allocatable, private :: p(:,:), q(:,:)
	real(kind=dp), dimension(:,:), allocatable, private :: &
		gu, gv, gphi_old, gphi, dgphi, dgphim, gphim,&
		midlon, midlat, deplon, deplat, gum, gvm
	complex(kind=dp), dimension(:,:), allocatable, private :: sphi1
	character(len=*), parameter, private :: hfile = "history.dat"

	private :: update
	public :: nisl_init, nisl_timeint, nisl_clean

contains

	subroutine nisl_init()
		implicit none

		integer(kind=i4b) :: i,j

		allocate(sphi1(0:ntrunc,0:ntrunc),gu(nlon,nlat),gv(nlon,nlat),gphi_old(nlon,nlat), &
		         gphi(nlon,nlat),gphim(nlon,nlat),dgphi(nlon,nlat),dgphim(nlon,nlat), &
		         midlon(nlon,nlat),midlat(nlon,nlat), &
		         deplon(nlon,nlat),deplat(nlon,nlat), p(nlon,nlat), q(nlon,nlat), &
		         gum(nlon,nlat),gvm(nlon,nlat))
    call interpolate_init(gphi)

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
		call calc_niuv(deltat)
		call update(deltat)
		if (hstep==1) then
			call legendre_synthesis(sphi, gphi)
			call save_data(hfile, 3*nsave+1, gphi, "old")
			call save_data(hfile, 3*nsave+2, gu, "old")
			call save_data(hfile, 3*nsave+3, gv, "old")
			nsave = nsave + 1
		end if
		call find_points(gu, gv, deltat, midlon, midlat, deplon, deplat)
		call calc_niuv(2.0_dp*deltat)

	end subroutine nisl_init

	subroutine nisl_clean()
		implicit none

		deallocate(sphi1,gu,gv,gphi,gphi_old,gphim,dgphi,dgphim,gum,gvm, &
			midlon,midlat,deplon,deplat,p,q)
    call interpolate_clean()

	end subroutine nisl_clean

	subroutine nisl_timeint()
		implicit none

		integer(kind=i4b) :: i

		do i=2, nstep
			print *, "step=", i, " hour=", real(i*deltat/hour_in_sec)
			call update(2.0_dp*deltat)
			if (mod(i,hstep)==0) then
				print *, "Saving step=", i
				call legendre_synthesis(sphi, gphi)
				call save_data(hfile, 3*nsave+1, gphi, "old")
				call save_data(hfile, 3*nsave+2, gu, "old")
				call save_data(hfile, 3*nsave+3, gv, "old")
				nsave = nsave + 1
			end if
		end do

	end subroutine nisl_timeint

	subroutine update(dt)
		implicit none

		integer(kind=i4b) :: i, j, m
		real(kind=dp), intent(in) :: dt

		! dF/dlon
		call legendre_synthesis_dlon(sphi, dgphi)
		do j= 1, nlat
			dgphi(:,j) = dgphi(:,j)/cos(latitudes(j))
		end do
		call interpolate_set(dgphi)
		do j=1, nlat
			do i=1, nlon
				call interpolate_bilinear(midlon(i,j), midlat(i,j), dgphim(i,j))
			end do
		end do
		gphim = gum*dgphim

		! cos(lat)dF/dlat
		call legendre_synthesis_dlat(sphi, dgphi) 
		do j= 1, nlat
			dgphi(:,j) = dgphi(:,j)/cos(latitudes(j))
		end do
		call interpolate_set(dgphi)
		do j=1, nlat
			do i=1, nlon
				call interpolate_bilinear(midlon(i,j), midlat(i,j), dgphim(i,j))
			end do
		end do
		gphim = gphim + gvm*dgphim

		call legendre_synthesis(sphi_old, gphi_old)
		do j=1, nlat
			do i=1, nlon
				gphi(i,j) = gphi_old(p(i,j),q(i,j))
			end do
		end do
		gphi = gphi + dt*gphim

! time filter
		call legendre_analysis(gphi, sphi1)
		do m=0, ntrunc
			sphi_old(m:ntrunc,m) = sphi(m:ntrunc,m) + &
				time_filter_param * (sphi_old(m:ntrunc,m)-2.0_dp*sphi(m:ntrunc,m)+sphi1(m:ntrunc,m))
			sphi(m:ntrunc,m) = sphi1(m:ntrunc,m)
		end do

	end subroutine update

	subroutine calc_niuv(dt)
		implicit none

		real(kind=dp), intent(in) :: dt

		integer(kind=i4b) :: i,j,ii
		real(kind=dp) :: xg, yg, zg, xr, yr, zr, xm, ym, zm, xd, yd, zd, &
		                 lon, lat, lonr, latr, lonm, latm, u, v, dlon, b

		dlon = 2.0_dp*pi/nlon
		call interpolate_setuv(gu,gv)
		do j=1, nlat
			lat = latitudes(j)
			do i=1, nlon
! find grid points near departure points
				p(i,j) = anint(deplon(i,j)/dlon+1.0_dp)
				if (p(i,j)>nlon) then
					p(i,j) = p(i,j)-nlon
				end if
! lat = (J+1-2j)pi/(2J+1)
        q(i,j) = anint(0.5_dp*(nlat+1-(2.0_dp*nlat+1.0_dp)*deplat(i,j)/pi))
				lonr = dlon * (p(i,j)-1)
				latr = latitudes(q(i,j))	
				call lonlat2xyz(lonr,latr,xr,yr,zr)
! arrival points
				lon = dlon * (i-1)
				call lonlat2xyz(lon,lat,xg,yg,zg)
! calculate midpoints between new departure points and arrival points
				b = 1.0_dp/sqrt(2.0_dp*(1.0_dp+(xg*xr+yg*yr+zg*zr)))
				xm = b*(xg + xr)
				ym = b*(yg + yr)
				zm = b*(zg + zr)
				midlon(i,j) = xy2lon(xm,ym)
				midlat(i,j) = asin(zm)
!				print *, real(lon*180/pi), real(midlon(i,j)*180/pi), real(lonr*180/pi), &
!					real(lat*180/pi), real(midlat(i,j)*180/pi), real(latr*180/pi)
! calculate velocities at midpoints
				xd = (xg-xr)/dt
				yd = (yg-yr)/dt
				zd = (zg-zr)/dt
				call xyz2uv(xd, yd, zd, midlon(i,j), midlat(i,j), u, v)
!				print *, real(lonr*180/pi), real(latr*180/pi), real(u*a), real(v*a)
!				print *, real(zg), real(zr), real(zd), real(v*a)
				gum(i,j) = u
				gvm(i,j) = v
				call interpolate_bilinearuv(midlon(i,j), midlat(i,j), u, v)
				gum(i,j) = gum(i,j) - u
				gvm(i,j) = gvm(i,j) - v
			end do
		end do
				
	end subroutine	calc_niuv

end module nisl_module
