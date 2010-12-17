module init_module

  use constant_module, only: i4b, dp, pi

  integer, parameter, private :: nc = 2
  real(kind=dp), dimension(nc), private :: &
    lonc = (/5.0_dp*pi/6.0_dp, 7.0_dp*pi/6.0_dp/), &
    latc = (/0.0_dp, 0.0_dp/)

  public :: init_ghill, init_ghill2

contains

! Ritche 1987
  subroutine init_ghill(lon,lat,gphi)
    use constant_module, only: a=>planet_radius, deg2rad
    use sphere_module, only: orthodrome
    implicit none

    real(kind=dp), dimension(:), intent(in) :: lon, lat
    real(kind=dp), dimension(:,:), intent(inout) :: gphi

    real(kind=dp), parameter :: &
      phimax = 100.0_dp, xi = 0.0_dp, yi = 0.0_dp, L = 2.5e6_dp ! Gaussian hill

    integer(kind=i4b) :: i, j, nx, ny
    real(kind=dp) :: r, loni, lati

    nx = size(lon)
    ny = size(lat)
    loni = xi*deg2rad
    lati = yi*deg2rad
    do j=1, ny
      do i=1, nx
        r = a*orthodrome(lon(i),lat(j),loni,lati)
        gphi(i,j) = phimax * exp(-(pi*r/L)**2)
      end do
    end do

  end subroutine init_ghill

  subroutine init_ghill2(lon,lat,gphi)
    use sphere_module, only: lonlat2xyz
    implicit none

    real(kind=dp), dimension(:), intent(in) :: lon, lat
    real(kind=dp), dimension(:,:), intent(inout) :: gphi

    real(kind=dp), parameter :: hmax = 0.95_dp, b0 = 5.0_dp

    integer(kind=i4b) :: i, j, k, nx, ny
    real(kind=dp) :: x0, y0, z0, x, y, z

    nx = size(lon)
    ny = size(lat)
    gphi(:,:) = 0.0_dp
    do k=1, nc
      call lonlat2xyz(lonc(k),latc(k),x0,y0,z0)
      do j=1, ny
        do i=1, nx
          call lonlat2xyz(lon(i),lat(j),x,y,z)
          gphi(i,j) = gphi(i,j) + hmax*exp(-b0*((x-x0)**2+(y-y0)**2+(z-z0)**2))
        end do
      end do
    end do

  end subroutine init_ghill2

end module init_module
