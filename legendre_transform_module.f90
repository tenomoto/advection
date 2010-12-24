module legendre_transform_module
  use kind_module, only: i4b, dp
  use glatwgt_module, only: glatwgt_calc
  use alf_module, only: pnm, alf_calc, alf_enm
  use fft_module, only: fft_init, fft_analysis, fft_synthesis
  use io_module, only: io_save
  implicit none
  private

  integer(kind=i4b), private :: nlon, nlat, ntrunc
  complex(kind=dp), dimension(:,:), allocatable, private :: w
  real(kind=dp), dimension(:), allocatable, private :: glat, gwgt

  public :: legendre_init, legendre_analysis, legendre_synthesis, &
            legendre_synthesis_dlon, legendre_synthesis_dlat, legendre_synthesis_dlonlat

contains

  subroutine legendre_init(nx,ny,nt,lat,wgt)
    implicit none

    integer(kind=i4b), intent(in) :: nx, ny, nt
    real(kind=dp), dimension(:), intent(inout) :: lat, wgt

    integer(kind=i4b) :: j, n, m
    real(kind=dp), dimension(:,:), allocatable :: g

    nlon = nx
    nlat = ny
    ntrunc = nt
    allocate(glat(nlat),gwgt(nlat))
    call glatwgt_calc(glat,gwgt)
    lat(:) = asin(glat(:))
    wgt(:) = gwgt(:)
    print *, "gaussian latitudes and weights"
    do j=1, nlat/2
      print *, j, lat(j), gwgt(j)
    end do
    call alf_calc(lat, ntrunc+1)
    print *, "first few elements of Pnm"
    do m=0, 5
      do n=m, 5
        print *, m, n, pnm(nlat/4,n,m)
      end do
    end do

    allocate(g(nlon,nlat), w(0:nlon/2,nlat))
    call fft_init(g, w)
    deallocate(g)
    call io_save("w.dat", 1, real(w,kind=dp), "replace")
    call io_save("w.dat", 2, aimag(w), "old")

  end subroutine legendre_init

  subroutine legendre_clean
    implicit none

    deallocate(glat,gwgt,w)

  end subroutine legendre_clean

  subroutine legendre_analysis(g,s)
    implicit none

    real(kind=dp), dimension(:,:), intent(in) :: g
    complex(kind=dp), dimension(0:,0:), intent(inout) :: s
    integer(kind=i4b) :: n, m, nlat2, flg 
    complex(kind=dp) :: tmp

    call fft_analysis(g, w)

    nlat2 = nlat/2

    do m=0, ntrunc
      flg = 1
      do n=m, ntrunc
        tmp = (0.0_dp, 0.0_dp)
        s(n, m) = sum( gwgt(1:nlat2)  * pnm(1:nlat2, n, m) * &
          (w(m, 1:nlat2) + flg * w(m, nlat:nlat2+1:-1)) )
        flg = -flg
      end do
      flg = -flg
    end do

  end subroutine legendre_analysis

  subroutine legendre_synthesis(s,g)
    implicit none

    complex(kind=dp), dimension(0:,0:), intent(in) :: s
    real(kind=dp), dimension(:,:), intent(inout) :: g

    integer(kind=i4b) :: n, m, j, jr, nlat2, flg
    real(kind=dp) :: nh, sh
    complex(kind=dp) :: ntmp, stmp

    nlat2 = nlat/2

    w = (0.0_dp, 0.0_dp)
    do j=1, nlat2
      jr = nlat - j + 1
      do m=0, ntrunc
        ntmp = (0.0_dp, 0.0_dp)
        stmp = ntmp
        flg = 1
        do n=m, ntrunc
          nh = pnm(j, n, m)
          sh = nh * flg
          ntmp = ntmp + nh * s(n, m)
          stmp = stmp + sh * s(n, m)
          flg = -flg
        end do
        w(m, j) = ntmp
        w(m, jr) = stmp
        flg = -flg
      end do
    end do

    call fft_synthesis(w,g)

  end subroutine legendre_synthesis

  subroutine legendre_synthesis_dlon(s,g)
    implicit none

    complex(kind=dp), dimension(0:,0:), intent(in) :: s
    real(kind=dp), dimension(:,:), intent(inout) :: g

    integer(kind=i4b) :: n, m, j, jr, nlat2, flg
    real(kind=dp) :: nh, sh
    complex(kind=dp) :: ntmp, stmp, im

    nlat2 = nlat/2

    w = (0.0_dp, 0.0_dp)
    do j=1, nlat2
      jr = nlat - j + 1
      do m=0, ntrunc
        im = dcmplx(0.0_dp, m)
        ntmp = (0.0_dp, 0.0_dp)
        stmp = ntmp
        flg = 1
        do n=m, ntrunc
          nh = pnm(j, n, m)
          sh = nh * flg
          ntmp = ntmp + im * nh * s(n, m)
          stmp = stmp + im * sh * s(n, m)
          flg = -flg
        end do
        w(m, j) = ntmp
        w(m, jr) = stmp
        flg = -flg
      end do
    end do

    call fft_synthesis(w,g)

  end subroutine legendre_synthesis_dlon

  subroutine legendre_synthesis_dlat(s,g)
    implicit none

    complex(kind=dp), dimension(0:,0:), intent(in) :: s
    real(kind=dp), dimension(:,:), intent(inout) :: g

    integer(kind=i4b) :: n, m, j, jr, nlat2, flg
    real(kind=dp) :: nh, sh
    complex(kind=dp) :: ntmp, stmp

    nlat2 = nlat/2

    w = (0.0_dp, 0.0_dp)
    do j=1, nlat2
      jr = nlat - j + 1
      do m=0, ntrunc
        ntmp = (0.0_dp, 0.0_dp)
        stmp = ntmp
        flg = -1
        do n=m, ntrunc
          nh = (n+1)*alf_enm(n,m)*pnm(j,n-1,m)-n*alf_enm(n+1,m)*pnm(j,n+1,m)
          sh = nh * flg
          ntmp = ntmp + nh * s(n, m)
          stmp = stmp + sh * s(n, m)
          flg = -flg
        end do
        w(m, j) = ntmp
        w(m, jr) = stmp
        flg = -flg
      end do
    end do

    call fft_synthesis(w,g)

  end subroutine legendre_synthesis_dlat

  subroutine legendre_synthesis_dlonlat(s,g)
    implicit none

    complex(kind=dp), dimension(0:,0:), intent(in) :: s
    real(kind=dp), dimension(:,:), intent(inout) :: g

    integer(kind=i4b) :: n, m, j, jr, nlat2, flg
    real(kind=dp) :: nh, sh
    complex(kind=dp) :: ntmp, stmp, im

    nlat2 = nlat/2

    w = (0.0_dp, 0.0_dp)
    do j=1, nlat2
      jr = nlat - j + 1
      do m=0, ntrunc
        im = dcmplx(0.0_dp, m)
        ntmp = (0.0_dp, 0.0_dp)
        stmp = ntmp
        flg = -1
        do n=m, ntrunc
          nh = (n+1)*alf_enm(n,m)*pnm(j,n-1,m)-n*alf_enm(n+1,m)*pnm(j,n+1,m)
          sh = nh * flg
          ntmp = ntmp + im * nh * s(n, m)
          stmp = stmp + im * sh * s(n, m)
          flg = -flg
        end do
        w(m, j) = ntmp
        w(m, jr) = stmp
        flg = -flg
      end do
    end do

    call fft_synthesis(w,g)

  end subroutine legendre_synthesis_dlonlat

end module legendre_transform_module
