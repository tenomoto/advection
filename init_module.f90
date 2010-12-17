module init_module

  use constant_module, only: i4b, dp

  public :: init_ghill

contains

  subroutine init_ghill(lon,lat,gphi)
    use constant_module, only: a=>planet_radius, pi, deg2rad
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

end module init_module
