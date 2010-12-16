module init_module

  use constant_module, only: i4b, dp, a=>planet_radius, pi, deg2rad
  use parameter_module, only: nlon, nlat, ntrunc
  use grid_module, only: lon, lat, gphi, gu, gv, sphi, sphi_old, su, sv, grid_init
  use legendre_transform_module, only: legendre_init, legendre_analysis, legendre_synthesis
  use sphere_module, only: orthodrome
  use uv_module, only: uv_sbody
  use io_module, only: save_data

  public :: init

contains

  subroutine init()
    implicit none

    real(kind=dp), parameter :: &
      phimax = 100.0_dp, xi = 0.0_dp, yi = 0.0_dp, L = 2.5e6_dp ! Gaussian hill

    integer(kind=i4b) :: i, j, n, m
    real(kind=dp) :: r, loni, lati

    loni = xi*deg2rad
    lati = yi*deg2rad
    do j=1, nlat
      do i=1, nlon
        r = a*orthodrome(lon(i),lat(j),loni,lati)
        gphi(i,j) = phimax * exp(-(pi*r/L)**2)
!       print *, i, j, real(r), real(gphi(i,j))
      end do
    end do
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

    call save_data("init.dat", 1, gphi, "replace")
    call save_data("init.dat", 2, gu, "old")
    call save_data("init.dat", 3, gv, "old")

   print *, "Umax=", real(maxval(gu)*a), " Umin=", real(minval(gu)*a)
   print *, "Vmax=", real(maxval(gv)*a), " Vmin=", real(minval(gv)*a)

  end subroutine init 

end module init_module
