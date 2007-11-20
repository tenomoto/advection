module upstream_module

! finds departure and mid-points

! Reference: Ritchie (1987)
! Method: use the Cartesian coordinates with the origin at the centre of the sphere

! Author: T. Enomoto
! History: 
! 26 February 2004     First written

	use constant_module, only: i4b, dp, pi, a=>planet_radius
	use glatwgt_module, only: latitudes=>lat
	use interpolate_module, only: interpolate_setuv, interpolate_bilinearuv
	use sphere_module, only: xy2lon, lonlat2xyz, uv2xyz
	private

	integer(kind=i4b), public :: itermax = 10
	real(kind=dp), public :: small = 1.0e-10

	public :: find_points

contains

	subroutine find_points(u, v, dt, midlon, midlat, deplon, deplat)
		implicit none

		real(kind=dp), dimension(:,:), intent(in) :: u, v
		real(kind=dp), intent(in) :: dt
		real(kind=dp), dimension(:,:), intent(inout) :: midlon, midlat, deplon, deplat

		integer(kind=i4b) :: nx, ny, i, j, k
		real(kind=dp) :: un, vn, &     ! normalised velocity
		                 b, &          ! correction factor
		                 xd, yd, zd, & ! Cartesian velocity
		                 xg, yg, zg, & ! arival point in Cartesian coordinates
		                 x0, y0, z0, & ! present point in Cartesian coordinates
		                 x1, y1, z1, & ! updated point in Cartesian coordinates
		                 lon, lat, err

		nx = size(u,1)
		ny = size(u,2)

		call interpolate_setuv(u,v)
		do j=1, ny
			do i=1, nx
! calculate initial values
				un = u(i,j)
				vn = v(i,j)
				lon = 2.0_dp*pi*(i-1)/nx ! calculate (lon,lat) from (i,j)
				lat = latitudes(j)
				call lonlat2xyz(lon, lat, xg, yg, zg) ! transform into Cartesian coordinates
				! r = g as an initial point for the 1st time step
				call lonlat2xyz(midlon(i,j), midlat(i,j), x0, y0, z0) 
				k = 1
				do 
					call uv2xyz(un,vn,lon,lat,xd,yd,zd) ! normalised Cartesian velocity
					! correction factor
					b = 1.0_dp/sqrt(1.0_dp+dt*dt*(xd*xd+yd*yd+zd*zd)-2.0_dp*dt*(xd*xg+yd*yg+zd*zg))
					x1 =  b*(xg - dt*xd) ! calculate new points
					y1 =  b*(yg - dt*yd)
					z1 =  b*(zg - dt*zd)
					! calculate (lon,lat) from (x,y,z)
					lat = asin(z1)
					lon = xy2lon(x1,y1)
          call interpolate_bilinearuv(lon, lat, un, vn) 
					err = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+(z1-z0)*(z1-z0)) ! calculate error
					x0 = x1 ! save as the current point
					y0 = y1
					z0 = z1
					k = k + 1
					if ((err<small).or.(k>itermax)) then
						exit
					end if
				end do
				midlon(i,j) = lon ! store as the mid-point
				midlat(i,j) = lat
				b = 2.0_dp*(x0*xg+y0*yg+z0*zg) ! calculate the departure point
        x1 = b*x0 - xg
				y1 = b*y0 - yg
				z1 = b*z0 - zg
				deplon(i,j) = xy2lon(x1,y1)
				deplat(i,j) = asin(z1)
			end do
		end do
		
	end subroutine find_points

end module upstream_module

