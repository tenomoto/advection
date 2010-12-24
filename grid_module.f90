module grid_module
  use kind_module, only: i4b, dp
  implicit none
  private

  character(len=6), private :: phi0 = "ghill2"
  character(len=6), public :: wind = "nodiv "
  integer(kind=i4b), public ::  ntrunc, nlon, nlat

  complex(kind=dp), dimension(:,:), allocatable, public :: sphi, sphi_old
  real(kind=dp), dimension(:,:), allocatable, public :: gphi, gu, gv
  real(kind=dp), dimension(:), allocatable, public :: lon, lat, coslat, coslatr, wgt

  public :: grid_init, grid_clean

contains

  subroutine grid_init()
    use math_module, only: pi2=>math_pi2
    use planet_module, only: a=>planet_radius
    use legendre_transform_module, only: &
      legendre_init, legendre_analysis
    use init_module, only: &
      init_ghill, init_ghill2, init_cbell2, init_scyli2, init_ccbel2
    use uv_module, only: uv_sbody, uv_nodiv
    implicit none

    namelist /grid/ ntrunc, nlon, nlat, phi0, wind

    integer(kind=i4b) :: i, j, m, n
    real(kind=dp) :: dlon

    read(unit=5, nml=grid)
    write(unit=6, nml=grid)

    allocate(lon(nlon), lat(nlat), coslat(nlat), coslatr(nlat), wgt(nlat), &
      gphi(nlon,nlat), gu(nlon,nlat), gv(nlon,nlat), &
      sphi(0:ntrunc,0:ntrunc), sphi_old(0:ntrunc,0:ntrunc))
 
    dlon = pi2/nlon
    do i=1, nlon
      lon(i) = dlon * (i-1.0_dp)
    end do
    call legendre_init(nlon,nlat,ntrunc,lat,wgt)
    coslat(:) = cos(lat(:))
    coslatr(:) = 1.0d0/coslat(:)

    select case(phi0)
      case("ghill ")
        call init_ghill(lon,lat,gphi)
      case("ghill2")
        call init_ghill2(lon,lat,gphi)
      case("cbell2")
        call init_cbell2(lon,lat,gphi)
      case("scyli2")
        call init_scyli2(lon,lat,gphi)
      case("ccbel2")
        call init_ccbel2(lon,lat,gphi)
      case default
        print *, "No matching initial phi"
        stop
    end select 
    call legendre_analysis(gphi, sphi)
    sphi_old = sphi
    print *, "first few elements of sphi"
    do m=0, 5
      do n=m, 5
        print *, m, n, sphi(n, m)
      end do
    end do

    select case(wind)
      case("sbody ")
        call uv_sbody(lon,lat,gu,gv)
      case("nodiv ")
        call uv_nodiv(0.0_dp,lon,lat,gu,gv)
      case default
        print *, "No matching initial wind"
        stop
    end select 
      
    print *, "Umax=", real(maxval(gu)*a), " Umin=", real(minval(gu)*a)
    print *, "Vmax=", real(maxval(gv)*a), " Vmin=", real(minval(gv)*a)

  end subroutine grid_init

  subroutine grid_clean
    implicit none

    deallocate(lon, lat, coslat, coslatr, gphi, gu, gv, sphi, sphi_old)

  end subroutine grid_clean

end module grid_module
