module grid_module
  use constant_module, only: i4b, dp, pi
  implicit none
  private

  complex(kind=dp), dimension(:,:), allocatable, public :: sphi, sphi_old!, su, sv
  real(kind=dp), dimension(:,:), allocatable, public :: gphi, gu, gv
  real(kind=dp), dimension(:), allocatable, public :: lon, lat

  public :: grid_init

contains

  subroutine grid_init(nlon,nlat,ntrunc)
    use constant_module, only: a=>planet_radius
    use legendre_transform_module, only: &
      legendre_init, legendre_analysis
    use io_module, only: io_save
    use init_module, only: init_ghill2
    use uv_module, only: uv_sbody, uv_nodiv
    implicit none

    integer(kind=i4b), intent(in) :: nlon, nlat, ntrunc

    integer(kind=i4b) :: i, j, m, n
    real(kind=dp) :: dlon

    allocate(lon(nlon), lat(nlat), gphi(nlon,nlat), gu(nlon,nlat), gv(nlon,nlat), &
      sphi(0:ntrunc,0:ntrunc), sphi_old(0:ntrunc,0:ntrunc))
!      su(0:ntrunc,0:ntrunc), sv(0:ntrunc,0:ntrunc), &
 
    dlon = 2.0_dp*pi/nlon
    do i=1, nlon
      lon(i) = dlon * (i-1.0_dp)
    end do
    call legendre_init(nlon,nlat,ntrunc,lat)

    call init_ghill2(lon,lat,gphi)
    call legendre_analysis(gphi, sphi)
    sphi_old = sphi
    print *, "first few elements of sphi"
    do m=0, 5
      do n=m, 5
        print *, m, n, sphi(n, m)
      end do
    end do

!    call uv_sbody(lon,lat,gu,gv)
    call uv_nodiv(0.0_dp,lon,lat,gu,gv)
!    print *, "Umax=", real(maxval(gu)*a), " Umin=", real(minval(gu)*a)
!    print *, "Vmax=", real(maxval(gv)*a), " Vmin=", real(minval(gv)*a)
    print *, "Umax=", real(maxval(gu)), " Umin=", real(minval(gu))
    print *, "Vmax=", real(maxval(gv)), " Vmin=", real(minval(gv))

    call io_save("init.dat", 1, gphi, "replace")
    call io_save("init.dat", 2, gu, "old")
    call io_save("init.dat", 3, gv, "old")

  end subroutine grid_init

  subroutine grid_clean
    implicit none

!    deallocate(lon, gphi, gu, gv, su, sv, sphi, sphi_old)
    deallocate(lon, gphi, gu, gv, sphi, sphi_old)

  end subroutine grid_clean

end module grid_module
