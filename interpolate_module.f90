module interpolate_module

! interpolate in a stencil

! References:
!   - Ritchie (1987) describes selection of points across the pole
!   - Interpolation schemes from Numerical Recipes
! Author: T. Enomoto
! History: 
! 2007-11-14 simplified stencil finding
! 2004-09-10 some simplification
! 2004-03

  use constant_module, only: i4b, dp, pi
  use glatwgt_module, only: latitudes=>lat
  use sphere_module, only: lon2i, lat2j
  private

  integer(kind=i4b), private :: nx, ny
  integer(kind=i4b), dimension(4), private :: is, js
  real(kind=dp), private :: u, t, dlon
  real(kind=dp), dimension(:), allocatable, private :: lonf, latf
  real(kind=dp), dimension(:,:), allocatable, private :: ff, ffx, ffy, ffxy, fu, fv

  private :: find_stencil
  public :: interpolate_init, interpolate_clean, &
            interpolate_set, interpolate_setuv, interpolate_setd, &
            interpolate_bilinear, interpolate_bilinearuv, &
            interpolate_bicubic, interpolate_polin2

contains

  subroutine interpolate_init(f)
    implicit none

    real(kind=dp), dimension(:,:) :: f

    integer(kind=i4b) :: i, j

    nx = size(f,1)
    ny = size(f,2)

    allocate(lonf(0:nx+2), latf(-1:ny+2), ff(0:nx+2,-1:ny+2), &
             ffx(0:nx+2,-1:ny+2), ffy(0:nx+2,-1:ny+2), ffxy(0:nx+2,-1:ny+2), &
             fu(0:nx+2,-1:ny+2),fv(0:nx+2,-1:ny+2))

    dlon = 2.0_dp*pi/nx
    do i=0, nx+2
      lonf(i) = dlon*(i-1)
    end do
    latf(-1)   = 0.5_dp*pi + (0.5_dp*pi - latitudes(2))
    latf(0)    = 0.5_dp*pi + (0.5_dp*pi - latitudes(1))
    latf(1:ny) = latitudes
    latf(ny+1) = -0.5_dp*pi + (-0.5_dp*pi - latitudes(ny))
    latf(ny+2) = -0.5_dp*pi + (-0.5_dp*pi - latitudes(ny-1))

  end subroutine interpolate_init

  subroutine interpolate_clean()
    implicit none

    deallocate(lonf, latf, ff, ffx, ffy, ffxy, fu, fv)

  end subroutine  interpolate_clean

  subroutine interpolate_bilinear(lon, lat, fi)
    implicit none

    real(kind=dp), intent(in) :: lon, lat
    real(kind=dp), intent(out) :: fi

    real(kind=dp), dimension(4) :: fs
    integer(kind=i4b) :: k 

    call find_stencil(lon, lat)
    do k=1, 4
      fs(k) = ff(is(k),js(k))
    end do

! Due to geographical location weight t and 1-t should be swapped.
! However, not so sensitive to t or 1-t at points 3 and 4,
! even slightly worse.  Leave it simple.

!    if (abs(lat)<latf(1)) then
      fi = (1.0_dp-u)*((1.0_dp-t)*fs(1)+t*fs(2)) + u*(t*fs(3)+(1.0_dp-t)*fs(4))
!    else
!       fi = (1.0_dp-u)*((1.0_dp-t)*fs(1)+t*fs(2)) + u*(t*fs(4)+(1.0_dp-t)*fs(3))
!    end if

  end subroutine interpolate_bilinear

  subroutine interpolate_bilinearuv(lon, lat, fiu, fiv)
    implicit none

    real(kind=dp), intent(in) :: lon, lat
    real(kind=dp), intent(out) :: fiu, fiv

    real(kind=dp), dimension(4) :: fsu, fsv
    integer(kind=i4b) :: k

    call find_stencil(lon, lat)
    do k=1, 4
      fsu(k) = fu(is(k),js(k))
      fsv(k) = fv(is(k),js(k))
    end do
!    if (abs(lat)<latf(1)) then
      fiu = (1.0_dp-u)*((1.0_dp-t)*fsu(1)+t*fsu(2)) + u*(t*fsu(3)+(1.0_dp-t)*fsu(4))
      fiv = (1.0_dp-u)*((1.0_dp-t)*fsv(1)+t*fsv(2)) + u*(t*fsv(3)+(1.0_dp-t)*fsv(4))
!    else
!      fiu = (1.0_dp-u)*((1.0_dp-t)*fsu(1)+t*fsu(2)) + u*(t*fsu(4)+(1.0_dp-t)*fsu(3))
!      fiv = (1.0_dp-u)*((1.0_dp-t)*fsv(1)+t*fsv(2)) + u*(t*fsv(4)+(1.0_dp-t)*fsv(3))
!    end if

  end subroutine interpolate_bilinearuv

  subroutine interpolate_bicubic(lon, lat, fi, monotonic)
    use bicubic_module, only: bcucof, bcuint, bcuintp
    implicit none

    real(kind=dp), intent(in) :: lon, lat
    real(kind=dp), intent(out) :: fi
    logical, optional, intent(in) :: monotonic

    real(kind=dp), dimension(4) :: z, zx, zy, zxy
    integer(kind=i4b) :: k
    real(kind=dp) :: dlat
    real(kind=dp), dimension(4) :: c

    call find_stencil(lon, lat)
    dlat = latf(js(4)) - latf(js(1))
    do k=1, 4
      z(k) = ff(is(k),js(k))
      zx(k) = ffx(is(k),js(k))
      zy(k) = ffy(is(k),js(k))
      zxy(k) = ffxy(is(k),js(k))
    end do

    call bcucof(z,zx,zy,zxy,dlon,dlat)
!    if (abs(lat)<latf(1)) then
      fi = bcuint(t,u)
!    else
!      fi = bcuintp(t,u)
!    end if

! Bermijo and Staniforth 1992
    if (present(monotonic).and.(monotonic)) then
      fi = min(fi,maxval(z))
      fi = max(fi,minval(z))
    end if

  end subroutine interpolate_bicubic

  subroutine interpolate_polin2(lon, lat, fi, monotonic)
    use polint_module, only : polin2
    implicit none

    real(kind=dp), intent(in) :: lon, lat
    real(kind=dp), intent(out) :: fi
    logical, optional, intent(in) :: monotonic

    real(kind=dp) :: dfi

    integer(kind=i4b) :: i1, i2, j1, j2

    call find_stencil(lon, lat)
    i1 = is(1) - 1
    i2 = is(1) + 2
    j1 = js(1) - 1
    j2 = js(1) + 2
    call polin2(lonf(i1:i2), latf(j1:j2), ff(i1:i2,j1:j2), lon, lat, fi, dfi)

! Bermijo and Staniforth 1992
    if (present(monotonic).and.(monotonic)) then
      fi = min(fi,maxval(ff(i1:i2,j1:j2)))
      fi = max(fi,minval(ff(i1:i2,j1:j2)))
    end if

  end subroutine interpolate_polin2

  subroutine interpolate_set(f)
    implicit none

    real(kind=dp), dimension(:,:), intent(in) :: f

    ff(1:nx,1:ny) = f
    ff(0,1:ny) = f(nx,:)
    ff(nx+1,1:ny) = f(1,:)
    ff(nx+2,1:ny) = f(2,:)
    ff(:,-1) = cshift(ff(:,2),nx/2)
    ff(:,0) = cshift(ff(:,1),nx/2)
    ff(:,ny+1) = cshift(ff(:,ny),nx/2)
    ff(:,ny+2) = cshift(ff(:,ny-1),nx/2)

  end subroutine interpolate_set

  subroutine interpolate_setuv(gu,gv)
    implicit none

    real(kind=dp), dimension(:,:), intent(in) :: gu, gv

    fu(1:nx,1:ny) = gu
    fu(0,1:ny) = gu(nx,:)
    fu(nx+1,1:ny) = gu(1,:)
    fu(nx+2,1:ny) = gu(2,:)
! direction of u is reversed beyond poles
    fu(:,-1) = -cshift(fu(:,2),nx/2)
    fu(:,0) = -cshift(fu(:,1),nx/2)
    fu(:,ny+1) = -cshift(fu(:,ny),nx/2)
    fu(:,ny+2) = -cshift(fu(:,ny-1),nx/2)

    fv(1:nx,1:ny) = gv
    fv(0,1:ny) = gv(nx,:)
    fv(nx+1,1:ny) = gv(1,:)
    fv(nx+2,1:ny) = gv(2,:)
! direction of v is reversed beyond poles
    fv(:,-1) = -cshift(fv(:,2),nx/2)
    fv(:,0) = -cshift(fv(:,1),nx/2)
    fv(:,ny+1) = -cshift(fv(:,ny),nx/2)
    fv(:,ny+2) = -cshift(fv(:,ny-1),nx/2)

  end subroutine interpolate_setuv

  subroutine interpolate_setd(fx,fy,fxy)
    implicit none

    real(kind=dp), dimension(:,:), intent(in) :: fx, fy, fxy

    ffx(1:nx,1:ny) = fx
    ffx(0,1:ny) = fx(nx,:)
    ffx(nx+1,1:ny) = fx(1,:)
    ffx(nx+2,1:ny) = fx(2,:)
! direction of d/dx is reversed beyond poles
    ffx(:,-1) = -cshift(ffx(:,2),nx/2)
    ffx(:,0) = -cshift(ffx(:,1),nx/2)
    ffx(:,ny+1) = -cshift(ffx(:,ny),nx/2)
    ffx(:,ny+2) = -cshift(ffx(:,ny-1),nx/2)

    ffy(1:nx,1:ny) = fy
    ffy(0,1:ny) = fy(nx,:)
    ffy(nx+1,1:ny) = fy(1,:)
    ffy(nx+2,1:ny) = fy(2,:)
! direction of d/dy is reversed beyond poles
    ffy(:,-1) = -cshift(ffy(:,2),nx/2)
    ffy(:,0) = -cshift(ffy(:,1),nx/2)
    ffy(:,ny+1) = -cshift(ffy(:,ny),nx/2)
    ffy(:,ny+2) = -cshift(ffy(:,ny-1),nx/2)

    ffxy(1:nx,1:ny) = fxy
    ffxy(0,1:ny) = fxy(nx,:)
    ffxy(nx+1,1:ny) = fxy(1,:)
    ffxy(nx+2,1:ny) = fxy(2,:)
    ffxy(:,-1) = cshift(ffxy(:,2),nx/2)
    ffxy(:,0) = cshift(ffxy(:,1),nx/2)
    ffxy(:,ny+1) = cshift(ffxy(:,ny),nx/2)
    ffxy(:,ny+2) = cshift(ffxy(:,ny-1),nx/2)

  end subroutine interpolate_setd

  subroutine find_stencil(lon, lat)
    implicit none

    real(kind=dp), intent(in) :: lon, lat
 
    integer(kind=i4b) :: j

    is(1) = lon2i(lon,nx)
    is(2) = is(1) + 1
    t = lon/dlon - is(1) + 1.0_dp ! t = (lon - dlon*(i-1))/dlon
    is(3:4) = is(2:1:-1)

    j = lat2j(lat,ny) 
    if (lat>latf(j)) then 
      j = j - 1
    end if
    js(1:2) = j
    js(3:4) = j + 1
    u = (lat-latf(j))/(latf(j+1)-latf(j))

  end subroutine find_stencil

end module interpolate_module
