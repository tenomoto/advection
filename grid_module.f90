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
    use constant_module, only: a=>planet_radius
    use legendre_transform_module, only: &
      legendre_init, legendre_analysis, legendre_synthesis
    use io_module, only: io_save
    use init_module, only: init_ghill
    use uv_module, only: uv_sbody
    implicit none

    integer(kind=i4b), intent(in) :: nlon, nlat, ntrunc

    integer(kind=i4b) :: i, j, m, n
    real(kind=dp) :: dlon

    allocate(lon(nlon), lat(nlat), gphi(nlon,nlat), gu(nlon,nlat), gv(nlon,nlat), &
      su(0:ntrunc,0:ntrunc), sv(0:ntrunc,0:ntrunc), &
      sphi(0:ntrunc,0:ntrunc), sphi_old(0:ntrunc,0:ntrunc))
 
    dlon = 2.0_dp*pi/nlon
    do i=1, nlon
      lon(i) = dlon * (i-1.0_dp)
    end do
    call legendre_init(nlon,nlat,ntrunc,lat)

    call init_ghill(lon,lat,gphi)
    call legendre_analysis(gphi, sphi)
    sphi_old = sphi
    print *, "first few elements of sphi"
    do m=0, 5
      do n=m, 5
        print *, m, n, sphi(n, m)
      end do
    end do

    call uv_sbody(lon,lat,gu,gv)
    do j=1, nlat
      do i=1, nlon
        gu(i,j) = gu(i,j) * cos(lat(j)) ! U = u*cos(lat)
        gv(i,j) = gv(i,j) * cos(lat(j)) ! V = v*cos(lat)
      end do
    end do
    call legendre_analysis(gu, su)
    call legendre_analysis(gv, sv)
    print *, "Umax=", real(maxval(gu)*a), " Umin=", real(minval(gu)*a)
    print *, "Vmax=", real(maxval(gv)*a), " Vmin=", real(minval(gv)*a)

    call io_save("init.dat", 1, gphi, "replace")
    call io_save("init.dat", 2, gu, "old")
    call io_save("init.dat", 3, gv, "old")

  end subroutine grid_init

  subroutine grid_clean
    implicit none

    deallocate(lon, gphi, gu, gv, su, sv, sphi, sphi_old)

  end subroutine grid_clean

end module grid_module
