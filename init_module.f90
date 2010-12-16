module init_module

  use constant_module, only: i4b, dp, a=>planet_radius, d=>day_in_sec, pi
  use parameter_module, only: nlon, nlat, ntrunc
  use glatwgt_module, only: lat
  use alf_module, only: pnm
  use legendre_transform_module, only: legendre_init, legendre_analysis, legendre_synthesis, w
  use sphere_module, only: orthodrome
  use io_module, only: save_data

  complex(kind=dp), dimension(:,:), allocatable, public :: sphi, sphi_old, su, sv

  public :: init

contains

  subroutine init()
    implicit none

    real(kind=dp), parameter :: &
      phimax = 100.0_dp, xi = 0.0_dp, yi = 0.0_dp, L = 2.5e6_dp, & ! Gaussian hill
      x0 = 0.0_dp, y0 = 45.0_dp, period = 20.0_dp                ! Rotation

    integer(kind=i4b) :: i, j, n, m
    real(kind=dp) :: r, lon, loni, lati, lon0, lat0, omg
    real(kind=dp), dimension(:,:), allocatable :: gphi, gu, gv

    call legendre_init()
    print *, "first few elements of Pnm"
    do m=0, 0
      do n=m, 5
        print *, m, n, pnm(nlat/4,n,m)
      end do
    end do

    allocate(gphi(nlon,nlat), gu(nlon,nlat), gv(nlon,nlat), &
      su(0:ntrunc,0:ntrunc), sv(0:ntrunc,0:ntrunc), &
      sphi(0:ntrunc,0:ntrunc), sphi_old(0:ntrunc,0:ntrunc))

    lon0 = x0*pi/180.0_dp
    lat0 = y0*pi/180.0_dp
    loni = xi*pi/180.0_dp
    lati = yi*pi/180.0_dp
    omg = (2.0_dp*pi) / (period*d)

    do j=1, nlat
      do i=1, nlon
        lon = 2*pi/nlon * (i-1)
        r = a*orthodrome(lon,lat(j),loni,lati)
        gphi(i,j) = phimax * exp(-(pi*r/L)**2)
!       print *, i, j, real(r), real(gphi(i,j))
      end do
    end do
    call legendre_analysis(gphi, sphi)
    sphi_old = sphi
    print *, "first few elements of sphi"
    do m=0, 5
!     print *, "m=",m, w(m, j0)
      do n=m, 5
        print *, m, n, sphi(n, m)
      end do
    end do
    call save_data("w.dat", 1, real(w,kind=dp), "replace")
    call save_data("w.dat", 2, aimag(w), "old")

    do j=1, nlat
      do i=1, nlon
        lon = 2*pi/nlon * (i-1)
        gu(i,j) = cos(lat(j))*sin(lat0)-cos(lon-lon0)*sin(lat(j))*cos(lat0)
        gv(i,j) = sin(lon-lon0)*cos(lat0)
        gu(i,j) = omg * gu(i,j) * cos(lat(j)) ! U = u*cos(lat)/a
        gv(i,j) = omg * gv(i,j) * cos(lat(j)) ! V = v*cos(lat)/a
      end do
    end do
    call legendre_analysis(gu, su)
    call legendre_analysis(gv, sv)

    call save_data("init.dat", 1, gphi, "replace")
    call save_data("init.dat", 2, gu, "old")
    call save_data("init.dat", 3, gv, "old")

!   print *, "Umax=", real(maxval(gu)*a), " Umin=", real(minval(gu)*a)
!   print *, "Vmax=", real(maxval(gv)*a), " Vmin=", real(minval(gv)*a)

  end subroutine init 

end module init_module
