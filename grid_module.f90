module grid_module
  use constant_module, only: i4b, dp, pi
  implicit none
  private

  complex(kind=dp), dimension(:,:), allocatable, public :: sphi, sphi_old, su, sv
  real(kind=dp), dimension(:,:), allocatable, public :: gphi, gu, gv
  real(kind=dp), dimension(:), allocatable, public :: lon, lat

  public :: grid_init

contains

  subroutine grid_init(nlon,nlat,ntrunc)
    implicit none

    integer(kind=i4b), intent(in) :: nlon, nlat, ntrunc

    integer(kind=i4b) :: i
    real(kind=dp) :: dlon

    allocate(lon(nlon), lat(nlat), gphi(nlon,nlat), gu(nlon,nlat), gv(nlon,nlat), &
      su(0:ntrunc,0:ntrunc), sv(0:ntrunc,0:ntrunc), &
      sphi(0:ntrunc,0:ntrunc), sphi_old(0:ntrunc,0:ntrunc))

    dlon = 2.0_dp*pi/nlon
    do i=1, nlon
      lon(i) = dlon * (i-1.0_dp)
    end do

  end subroutine grid_init

  subroutine grid_clean
    implicit none

    deallocate(lon, gphi, gu, gv, su, sv, sphi, sphi_old)

  end subroutine grid_clean

end module grid_module
