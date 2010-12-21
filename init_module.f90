module init_module

  use constant_module, only: i4b, dp, pi

  integer, parameter, private :: nc = 2
  real(kind=dp), dimension(nc), private :: &
    lonc = (/5.0_dp*pi/6.0_dp, 7.0_dp*pi/6.0_dp/), &
    latc = (/0.0_dp, 0.0_dp/)

  public :: init_ghill, init_ghill2, init_cbell2, init_scyli2, init_ccbel2

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


  subroutine init_cbell2(lon,lat,gphi)
    use sphere_module, only: orthodrome
    implicit none

    real(kind=dp), dimension(:), intent(in) :: lon, lat
    real(kind=dp), dimension(:,:), intent(inout) :: gphi

    real(kind=dp), parameter :: &
      hmax = 1.0_dp, r0 = 0.5_dp, b = 0.1_dp, c = 0.9_dp

    integer(kind=i4b) :: i, j, k, nx, ny
    real(kind=dp) :: hch, pir0r, r

    nx = size(lon)
    ny = size(lat)
    gphi(:,:) = b
    pir0r = pi/r0
    hch = c*0.5_dp*hmax
    do k=1, nc
      do j=1, ny
        do i=1, nx
          r = orthodrome(lon(i),lat(j),lonc(k),latc(k))
          if (r<r0) then
            gphi(i,j) = gphi(i,j) + hch*(1.0_dp+cos(pir0r*r))
          end if
        end do
      end do
    end do

  end subroutine init_cbell2

  subroutine init_scyli2(lon,lat,gphi)
    use sphere_module, only: orthodrome
    implicit none

    real(kind=dp), dimension(:), intent(in) :: lon, lat
    real(kind=dp), dimension(:,:), intent(inout) :: gphi

    real(kind=dp), parameter :: &
      hmax = 1.0_dp, r0 = 0.5_dp, b = 0.1_dp, c = 1.0_dp

    integer(kind=i4b) :: i, j, nx, ny
    real(kind=dp) :: &
      r1, r2, dlon0, dlon1, dlon2, dlat0, dlat1, dlat2

    nx = size(lon)
    ny = size(lat)
    gphi(:,:) = b
    dlon0 = r0/6.0_dp
    dlat0 = 5.0_dp/12.0_dp*r0
    do j=1, ny
      do i=1, nx
        r1 = orthodrome(lon(i),lat(j),lonc(1),latc(1))
        r2 = orthodrome(lon(i),lat(j),lonc(2),latc(2))
        dlon1 = abs(lon(i)-lonc(1))
        dlon2 = abs(lon(i)-lonc(2))
        dlat1 = lat(j)-latc(1)
        dlat2 = lat(j)-latc(2)
        if (((r1<=r0).and.(dlon1>=dlon0)).or. &
            ((r2<=r0).and.(dlon2>=dlon0)).or. &
            ((r1<=r0).and.(dlon1<dlon0).and.dlat1<-dlat0).or. &
            ((r2<=r0).and.(dlon2<dlon0).and.dlat2>dlat0)) then
          gphi(i,j) = c
        end if
      end do
    end do

  end subroutine init_scyli2

  subroutine init_ccbel2(lon,lat,gphi)
    use sphere_module, only: orthodrome
    implicit none

    real(kind=dp), dimension(:), intent(in) :: lon, lat
    real(kind=dp), dimension(:,:), intent(inout) :: gphi

    real(kind=dp), parameter :: a = -0.8, b = 0.9

    call init_cbell2(lon,lat,gphi)
    gphi(:,:) = a*gphi(:,:)**2 + b

  end subroutine init_ccbel2

end module init_module
