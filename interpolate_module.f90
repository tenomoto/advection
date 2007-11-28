module interpolate_module

! interpolate in a stencil

! References:
!   - Ritchie (1987) describes selection of points across the pole
!   - Interpolation schemes from Numerical Recipes
! Author: T. Enomoto
! History: 
! 2007-11-20 higher order Lagrange interpolation
! 2007-11-14 simplified stencil finding
! 2004-09-10 some simplification
! 2004-03

  use constant_module, only: i4b, dp, pi
  use glatwgt_module, only: latitudes=>lat
  use sphere_module, only: lon2i, lat2j
  private

  integer(kind=i4b), private :: nx, ny, n=3, nh, nx1, nx2, ny1, ny2 
  integer(kind=i4b), dimension(4), private :: is, js
  real(kind=dp), private :: u, t, dlon
  real(kind=dp), dimension(:), allocatable, private :: lonf, latf
  real(kind=dp), dimension(:,:), allocatable, private :: ff, ffx, ffy, ffxy, fu, fv

  private :: find_stencil
  public :: interpolate_init, interpolate_clean, &
            interpolate_set, interpolate_setuv, interpolate_setd, &
            interpolate_bilinear, interpolate_bilinearuv, &
            interpolate_bicubic, interpolate_polin2, interpolate_linpol, &
            interpolate_setdx, interpolate_spcher

contains

  subroutine interpolate_init(f,k)
    implicit none

    real(kind=dp), dimension(:,:) :: f
    integer(kind=i4b), optional :: k

    integer(kind=i4b) :: i, j

    nx = size(f,1)
    ny = size(f,2)
    if (present(k)) then
      n = k
    end if
    nh = n/2
    nx1 = 1 - nh
    nx2 = nx + nh + 1
    ny1 = 1 - nh - 1
    ny2 = ny + nh + 1

    allocate(lonf(nx1:nx2), latf(ny1:ny2), ff(nx1:nx2,ny1:ny2), &
             ffx(nx1:nx2,ny1:ny2), ffy(nx1:nx2,ny1:ny2), ffxy(nx1:nx2,ny1:ny2), &
             fu(nx1:nx2,ny1:ny2),fv(nx1:nx2,ny1:ny2))

    dlon = 2.0_dp*pi/nx
    do i=nx1, nx2
      lonf(i) = dlon*(i-1)
    end do
    latf(1:ny) = latitudes
    do j=1, nh+1
      latf(1-j)   = 0.5_dp*pi + (0.5_dp*pi - latitudes(j))
      latf(ny+j)   = -0.5_dp*pi + (-0.5_dp*pi - latitudes(ny-j+1))
    end do

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

! Bermejo and Staniforth 1992
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

    integer(kind=i4b) :: i0, i1, i2, j0, j1, j2
    real(kind=dp) :: dfi

    call find_stencil(lon, lat)
    i0 = is(1)
    i1 = i0 - nh
    i2 = i0 + nh + 1
    j0 = js(1)
    j1 = j0 - nh
    j2 = j0 + nh + 1
    call polin2(lonf(i1:i2), latf(j1:j2), ff(i1:i2,j1:j2), lon, lat, fi, dfi)

! Bermejo and Staniforth 1992
    if (present(monotonic).and.(monotonic)) then
      fi = min(fi,maxval(ff(i0:i0+1,j0:j0+1)))
      fi = max(fi,minval(ff(i0:i0+1,j0:j0+1)))
    end if

  end subroutine interpolate_polin2

  subroutine interpolate_linpol(lon, lat, fi, monotonic)
    use polint_module, only : polint
    implicit none

    real(kind=dp), intent(in) :: lon, lat
    real(kind=dp), intent(out) :: fi
    logical, optional, intent(in) :: monotonic

    integer(kind=i4b) :: i0, j0, j1, j2, j
    real(kind=dp) :: dfi
    real(kind=dp), dimension(n+1) :: ytmp

    call find_stencil(lon, lat)
    i0 = is(1)
    j0 = js(1)
    j1 = j0 - nh
    j2 = j0 + nh + 1
    do j=j1, j2
      ytmp(j-j1+1) =  (1.0_dp-t)*ff(i0,j)+t*ff(i0+1,j)
    end do
    call polint(latf(j1:j2), ytmp, lat, fi, dfi)

! Bermejo and Staniforth 1992
    if (present(monotonic).and.(monotonic)) then
      fi = min(fi,maxval(ff(i0:i0+1,j0:j0+1)))
      fi = max(fi,minval(ff(i0:i0+1,j0:j0+1)))
    end if

  end subroutine interpolate_linpol

  subroutine interpolate_spcher(lon, lat, fi, monotonic)
    use cubicspline_module, only : cubicspline_interpolate
    use polint_module, only : polint
    implicit none

    real(kind=dp), intent(in) :: lon, lat
    real(kind=dp), intent(out) :: fi
    logical, optional, intent(in) :: monotonic

    integer(kind=i4b) :: i0, j0, j1, j2, j
    real(kind=dp) :: dfi
    real(kind=dp), dimension(4) :: fs
    real(kind=dp), dimension(n+1) :: ytmp

    call find_stencil(lon, lat)
    i0 = is(1)
    j0 = js(1)
    j1 = j0 - nh
    j2 = j0 + nh + 1
    do j=j1, j2
      fs(1) = ff(i0,j)
      fs(2) = ff(i0+1,j)
      fs(3) = ffx(i0,j)
      fs(4) = ffx(i0+1,j)
      call cubicspline_interpolate(t, fs, ytmp(j-j1+1))
    end do
    call polint(latf(j1:j2), ytmp, lat, fi, dfi)

! Bermejo and Staniforth 1992
    if (present(monotonic).and.(monotonic)) then
      fi = min(fi,maxval(ff(i0:i0+1,j0:j0+1)))
      fi = max(fi,minval(ff(i0:i0+1,j0:j0+1)))
    end if

  end subroutine interpolate_spcher

  subroutine interpolate_set(f)
    implicit none

    real(kind=dp), dimension(:,:), intent(in) :: f

    integer(kind=i4b) :: i, j

    ff(1:nx,1:ny) = f
    do i=1, nh
      ff(1-i,1:ny) = f(nx-(i-1),:)
    end do
    do i=1, nh+1
      ff(nx+i,1:ny) = f(1+(i-1),:)
    end do
    do j=1, nh+1
      ff(:,1-j) = cshift(ff(:,j),nx/2)
      ff(:,ny+j) = cshift(ff(:,ny-(j-1)),nx/2)
    end do

  end subroutine interpolate_set

  subroutine interpolate_setuv(gu,gv)
    implicit none

    real(kind=dp), dimension(:,:), intent(in) :: gu, gv

    integer(kind=i4b) :: i, j

    fu(1:nx,1:ny) = gu
    fv(1:nx,1:ny) = gv
    do i=1, nh
      fu(1-i,1:ny) = gu(nx-(i-1),:)
      fv(1-i,1:ny) = gv(nx-(i-1),:)
    end do
    do i=1, nh+1
      fu(nx+i,1:ny) = gu(1+(i-1),:)
      fv(nx+i,1:ny) = gv(1+(i-1),:)
    end do
! direction of u, v is reversed beyond poles
    do j=1, nh+1
      fu(:,1-j) = -cshift(fu(:,j),nx/2)
      fu(:,ny+j) = -cshift(fu(:,ny-(j-1)),nx/2)
      fv(:,1-j) = -cshift(fv(:,j),nx/2)
      fv(:,ny+j) = -cshift(fv(:,ny-(j-1)),nx/2)
    end do

  end subroutine interpolate_setuv

  subroutine interpolate_setd(fx,fy,fxy)
    implicit none

    real(kind=dp), dimension(:,:), intent(in) :: fx, fy, fxy

    integer(kind=i4b) :: i, j

    ffx(1:nx,1:ny) = fx
    ffy(1:nx,1:ny) = fy
    ffxy(1:nx,1:ny) = fxy
    do i=1, nh
      ffx(1-i,1:ny) = fx(nx-(i-1),:)
      ffy(1-i,1:ny) = fy(nx-(i-1),:)
      ffxy(1-i,1:ny) = fxy(nx-(i-1),:)
    end do
    do i=1, nh+1
      ffx(nx+i,1:ny) = fx(1+(i-1),:)
      ffy(nx+i,1:ny) = fy(1+(i-1),:)
      ffxy(nx+i,1:ny) = fxy(1+(i-1),:)
    end do
! directions of d/dx and d/dy are reversed beyond poles
    do j=1, nh+1
      ffx(:,1-j) = -cshift(ffx(:,j),nx/2)
      ffx(:,ny+j) = -cshift(ffx(:,ny-(j-1)),nx/2)
      ffy(:,1-j) = -cshift(ffy(:,j),nx/2)
      ffy(:,ny+j) = -cshift(ffy(:,ny-(j-1)),nx/2)
      ffxy(:,1-j) = cshift(ffxy(:,j),nx/2)
      ffxy(:,ny+j) = cshift(ffxy(:,ny-(j-1)),nx/2)
    end do

  end subroutine interpolate_setd

  subroutine interpolate_setdx(fx)
    implicit none

    real(kind=dp), dimension(:,:), intent(in) :: fx

    integer(kind=i4b) :: i, j

    ffx(1:nx,1:ny) = fx
    do i=1, nh
      ffx(1-i,1:ny) = fx(nx-(i-1),:)
    end do
    do i=1, nh+1
      ffx(nx+i,1:ny) = fx(1+(i-1),:)
    end do
! direction of d/dx is reversed beyond poles
    do j=1, nh+1
      ffx(:,1-j) = -cshift(ffx(:,j),nx/2)
      ffx(:,ny+j) = -cshift(ffx(:,ny-(j-1)),nx/2)
    end do
    ffy(:,:) = 0.d0
    ffxy(:,:) = 0.d0

  end subroutine interpolate_setdx

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
