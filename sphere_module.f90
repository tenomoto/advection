module sphere_module
! utility for a spherical topology
	use constant_module, only: i4b, dp, pi
	private

	public :: xy2lon, lonlat2xyz, uv2xyz, xyz2uv, lon2i, lat2j, orthodrome

contains

	subroutine lonlat2xyz(lon, lat, x, y, z)
		implicit none

		real(kind=dp), intent(in) :: lon, lat
		real(kind=dp), intent(out) :: x, y, z

    x = cos(lon)*cos(lat)
		y = sin(lon)*cos(lat)
		z = sin(lat)

	end subroutine lonlat2xyz

	subroutine uv2xyz(u, v, lon, lat, xd, yd, zd)
		implicit none

		real(kind=dp), intent(in) :: u, v, lon, lat
		real(kind=dp), intent(out) :: xd, yd, zd

		xd = -u*sin(lon) - v*cos(lon)*sin(lat)
		yd =  u*cos(lon) - v*sin(lon)*sin(lat)
		zd =  v*cos(lat)

	end subroutine uv2xyz

	subroutine xyz2uv(xd, yd, zd, lon, lat, u, v)
		implicit none

		real(kind=dp), intent(in) :: xd, yd, zd, lon, lat
		real(kind=dp), intent(out) :: u, v

		u = cos(lon)*yd - sin(lon)*xd
		if (abs(lat) == pi/2) then
			v = -cos(lon)*xd-sin(lon)*yd ! omitted division by sin(pi/2) = 1
		else
			v = zd/cos(lat)
		end if

	end subroutine xyz2uv

	function xy2lon(x,y) result(lon)
		implicit none

		real(kind=dp), intent(in) :: x, y
		real(kind=dp) :: lon

		if (x==0.0_dp) then
			if (y>=0.0_dp) then
				lon = pi/2.0_dp
			else
				lon = -pi/2.0_dp
			end if
		else
			lon = atan(y/x)
			if (x<0.0_dp) then
				lon = lon + pi
			else if (y<0.0_dp) then ! x1 > 0.0
				lon = lon + 2.0_dp*pi
			end if
		end if

	end function xy2lon

	function lon2i(lon,nx) result(i)
	! returns the closest longitudinal point i, not exceeding lon
  ! 1 <= return value <= nx
		implicit none

		real(kind=dp), intent(in) :: lon
		integer(kind=i4b), intent(in) :: nx

		integer(kind=i4b) :: i
		real(kind=dp) :: dlon

		dlon = 2.0_dp*pi/nx
		i = floor(lon/dlon+1.0_dp) ! lon = 2pi/nx*(i-1)=dlon*(i-1)

	end function lon2i

	function lat2j(lat,ny) result(j)
	! returns the closest Gaussian point j using approximation
	! lat varies from NP to SP
  ! 1 <= return value <= ny
		implicit none

		real(kind=dp), intent(in) :: lat
		integer(kind=i4b) :: ny

		integer(kind=i4b) :: j

		j = anint(0.5_dp*(ny+1-(2.0_dp*ny+1.0_dp)*lat/pi)) ! lat = (J+1-2j)pi/(2J+1)

	end function lat2j

	function orthodrome(x1, y1, x2, y2) result(l)
	! returns the shortest distance between two points on the unit sphere
		implicit none

		real(kind=dp), intent(in) :: x1, y1, x2, y2
		real(kind=dp) :: l

		l = acos(cos(x1-x2)*cos(y1)*cos(y2)+sin(y1)*sin(y2))

	end function orthodrome

end module sphere_module
