module interpolate_module
! interpolate in a stencil
  use kind_module, only: i4b, dp
  use grid_module, only: latitudes=>lat, coslat, wgt
  use sphere_module, only: lon2i, lat2j
  implicit none
  private

! References:
!   - Ritchie (1987) describes selection of points across the pole
!   - Interpolation schemes from Numerical Recipes
! Author: T. Enomoto
! History: 
! 2007-11-20 higher order Lagrange interpolation
! 2007-11-14 simplified stencil finding
! 2004-09-10 some simplification
! 2004-03

  integer(kind=i4b), private :: nx, ny, n=3, nh, nx1, nx2, ny1, ny2 
  integer(kind=i4b), dimension(4), private :: is, js
  real(kind=dp), private :: u, t, dlon
  real(kind=dp), dimension(:), allocatable, private :: lonf, latf
  real(kind=dp), dimension(:,:), allocatable, private :: &
    ff, ffx, ffy, ffxy, fu, fv, ffxl, ffyl

  private :: find_stencil, quasimonotone_filter, quasimonotone_filter_bicubic
  public :: interpolate_init, interpolate_clean, &
            interpolate_set, interpolate_setuv, &
            interpolate_setd, interpolate_setdx, &
            interpolate_bilinear, interpolate_bilinearuv, &
            interpolate_polin2, interpolate_polin2uv, &
            interpolate_bicubic, interpolate_linpol, &
            interpolate_spcher, interpolate_diff

contains

  subroutine interpolate_init(f)
    use math_module, only: pi2=>math_pi2, pih=>math_pih
    implicit none

    real(kind=dp), dimension(:,:) :: f

    integer(kind=i4b) :: i, j

    namelist /interpolate/ n

    read(unit=5, nml=interpolate)
    write(unit=6, nml=interpolate)

    nx = size(f,1)
    ny = size(f,2)
    nh = n/2
    nx1 = 1 - nh
    nx2 = nx + nh + 1
    ny1 = 1 - nh - 1
    ny2 = ny + nh + 1

    allocate(lonf(nx1:nx2), latf(ny1:ny2), ff(nx1:nx2,ny1:ny2), &
             ffx(nx1:nx2,ny1:ny2), ffy(nx1:nx2,ny1:ny2), ffxy(nx1:nx2,ny1:ny2), &
             ffxl(nx1:nx2,ny1:ny2), ffyl(nx1:nx2,ny1:ny2), &
             fu(nx1:nx2,ny1:ny2),fv(nx1:nx2,ny1:ny2))

    dlon = pi2/nx
    do i=nx1, nx2
      lonf(i) = dlon*(i-1)
    end do
    latf(1:ny) = latitudes
    do j=1, nh+1
      latf(1-j)   = pih + (pih - latitudes(j))
      latf(ny+j)   = -pih + (-pih - latitudes(ny-j+1))
    end do

  end subroutine interpolate_init

  subroutine interpolate_clean()
    implicit none

    deallocate(lonf, latf, ff, ffx, ffy, ffxy, fu, fv, ffxl, ffyl)

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

    real(kind=dp) :: coslatr
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
!      coslatr = 1.0_dp/cos(lat)
!      fiu = fiu*coslatr
!      fiv = fiv*coslatr
!    else
!      fiu = (1.0_dp-u)*((1.0_dp-t)*fsu(1)+t*fsu(2)) + u*(t*fsu(4)+(1.0_dp-t)*fsu(3))
!      fiv = (1.0_dp-u)*((1.0_dp-t)*fsv(1)+t*fsv(2)) + u*(t*fsv(4)+(1.0_dp-t)*fsv(3))
!    end if

  end subroutine interpolate_bilinearuv

  subroutine interpolate_bicubic(lon, lat, fi, monotonic, minmax)
    use bicubic_module, only: bcucof, bcuint, bcuintp
    implicit none

    real(kind=dp), intent(in) :: lon, lat
    real(kind=dp), intent(out) :: fi
    logical, optional, intent(in) :: monotonic
    real(kind=dp), optional, dimension(2), intent(in) :: minmax

    real(kind=dp), dimension(4) :: z, zx, zy, zxy
    integer(kind=i4b) :: k
    real(kind=dp) :: dlat

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

    if (present(monotonic).and.(monotonic)) then
!      call quasimonotone_filter_bicubic(fi,is(1),js(1))
      call quasimonotone_filter(fi,is(1),js(1))
    end if
    if (present(minmax).and.(minmax(1)/=minmax(2))) then
! Ostiguy and Laprise 1990
      fi = max(minmax(1),min(minmax(2),fi))
    end if
   
  end subroutine interpolate_bicubic

  subroutine interpolate_polin2(lon, lat, fi, monotonic, minmax)
    use polint_module, only : polin2
    implicit none

    real(kind=dp), intent(in) :: lon, lat
    real(kind=dp), intent(out) :: fi
    logical, optional, intent(in) :: monotonic
    real(kind=dp), optional, dimension(2), intent(in) :: minmax

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


  end subroutine interpolate_polin2

  subroutine interpolate_polin2uv(lon, lat, fiu, fiv, monotonic)
    use polint_module, only : polin2
    implicit none

    real(kind=dp), intent(in) :: lon, lat
    real(kind=dp), intent(out) :: fiu, fiv
    logical, optional, intent(in) :: monotonic

    integer(kind=i4b) :: i0, i1, i2, j0, j1, j2
    real(kind=dp) :: dfi, coslatr

    call find_stencil(lon, lat)
    i0 = is(1)
    i1 = i0 - nh
    i2 = i0 + nh + 1
    j0 = js(1)
    j1 = j0 - nh
    j2 = j0 + nh + 1
    call polin2(lonf(i1:i2), latf(j1:j2), fu(i1:i2,j1:j2), lon, lat, fiu, dfi)
    call polin2(lonf(i1:i2), latf(j1:j2), fv(i1:i2,j1:j2), lon, lat, fiv, dfi)

! Bermejo and Staniforth 1992
    if (present(monotonic).and.(monotonic)) then
      fiu = min(fiu,maxval(fu(i0:i0+1,j0:j0+1)))
      fiu = max(fiu,minval(fu(i0:i0+1,j0:j0+1)))
      fiv = min(fiv,maxval(fv(i0:i0+1,j0:j0+1)))
      fiv = max(fiv,minval(fv(i0:i0+1,j0:j0+1)))
    end if
!    coslatr = 1.0_dp/cos(lat)
!    fiu = fiu * coslatr
!    fiv = fiv * coslatr

  end subroutine interpolate_polin2uv

  subroutine interpolate_linpol(lon, lat, fi, monotonic, minmax)
    use polint_module, only : polint
    implicit none

    real(kind=dp), intent(in) :: lon, lat
    real(kind=dp), intent(out) :: fi
    logical, optional, intent(in) :: monotonic
    real(kind=dp), optional, dimension(2), intent(in) :: minmax

    integer(kind=i4b) :: i0, i1, i2, j0, j1, j2, j
    real(kind=dp) :: dfi
    real(kind=dp), dimension(n+1) :: ytmp

    call find_stencil(lon, lat)
    i0 = is(1)
    j0 = js(1)
    i1 = i0 - nh
    i2 = i0 + nh + 1
    j1 = j0 - nh
    j2 = j0 + nh + 1
    do j=j1, j0-1
      ytmp(j-j1+1) =  (1.0_dp-t)*ff(i0,j)+t*ff(i0+1,j)
    end do
    do j=j0, j0+1
      call polint(lonf(i1:i2), ff(i1:i2,j), lon, ytmp(j-j1+1), dfi)
    end do
    do j=j0+2, j2
      ytmp(j-j1+1) =  (1.0_dp-t)*ff(i0,j)+t*ff(i0+1,j)
    end do
    call polint(latf(j1:j2), ytmp, lat, fi, dfi)

    if (present(monotonic).and.(monotonic)) then
      call quasimonotone_filter(fi,is(1),js(1))
    end if
    if (present(minmax).and.(minmax(1)/=minmax(2))) then
! Ostiguy and Laprise 1990
      fi = max(minmax(1),min(minmax(2),fi))
    end if

  end subroutine interpolate_linpol

  subroutine interpolate_spcher(lon, lat, fi, monotonic, minmax)
    use cubicspline_module, only : cubicspline_interpolate
    use polint_module, only : polint
    implicit none

    real(kind=dp), intent(in) :: lon, lat
    real(kind=dp), intent(out) :: fi
    logical, optional, intent(in) :: monotonic
    real(kind=dp), optional, dimension(2), intent(in) :: minmax

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
      call cubicspline_interpolate(t, dlon, fs, ytmp(j-j1+1))
    end do
    call polint(latf(j1:j2), ytmp, lat, fi, dfi)

    if (present(monotonic).and.(monotonic)) then
      call quasimonotone_filter(fi,is(1),js(1))
    end if
    if (present(minmax).and.(minmax(1)/=minmax(2))) then
! Ostiguy and Laprise 1990
      fi = max(minmax(1),min(minmax(2),fi))
    end if

  end subroutine interpolate_spcher

  subroutine interpolate_set(f)
    implicit none

    real(kind=dp), dimension(:,:), intent(in) :: f

    integer(kind=i4b) :: i, j

    ff(1:nx,1:ny) = f
    do j=1, nh+1
      ff(1:nx,1-j) = cshift(ff(1:nx,j),nx/2)
      ff(1:nx,ny+j) = cshift(ff(1:nx,ny-(j-1)),nx/2)
    end do
    do i=1, nh
      ff(1-i,:) = ff(nx-(i-1),:)
    end do
    do i=1, nh+1
      ff(nx+i,:) = ff(1+(i-1),:)
    end do

  end subroutine interpolate_set

  subroutine interpolate_setuv(gu,gv)
    implicit none

    real(kind=dp), dimension(:,:), intent(in) :: gu, gv

    integer(kind=i4b) :: i, j

    do j=1, ny
      fu(1:nx,j) = gu(:,j)!*coslat(j)
      fv(1:nx,j) = gv(:,j)!*coslat(j)
    end do
! direction of u, v is reversed beyond poles
    do j=1, nh+1
      fu(1:nx,1-j) = -cshift(fu(1:nx,j),nx/2)
      fu(1:nx,ny+j) = -cshift(fu(1:nx,ny-(j-1)),nx/2)
      fv(1:nx,1-j) = -cshift(fv(1:nx,j),nx/2)
      fv(1:nx,ny+j) = -cshift(fv(1:nx,ny-(j-1)),nx/2)
    end do
    do i=1, nh
      fu(1-i,:) = fu(nx-(i-1),:)
      fv(1-i,:) = fv(nx-(i-1),:)
    end do
    do i=1, nh+1
      fu(nx+i,:) = fu(1+(i-1),:)
      fv(nx+i,:) = fv(1+(i-1),:)
    end do

  end subroutine interpolate_setuv

  subroutine interpolate_setd(fx,fy,fxy)
    implicit none

    real(kind=dp), dimension(:,:), intent(in) :: fx, fy, fxy

    integer(kind=i4b) :: i, j

    ffx(1:nx,1:ny) = fx
    ffy(1:nx,1:ny) = fy
    ffxy(1:nx,1:ny) = fxy
! directions of d/dx and d/dy are reversed beyond poles
    do j=1, nh+1
      ffx(1:nx,1-j) = -cshift(ffx(1:nx,j),nx/2)
      ffx(1:nx,ny+j) = -cshift(ffx(1:nx,ny-(j-1)),nx/2)
      ffy(1:nx,1-j) = -cshift(ffy(1:nx,j),nx/2)
      ffy(1:nx,ny+j) = -cshift(ffy(1:nx,ny-(j-1)),nx/2)
      ffxy(1:nx,1-j) = cshift(ffxy(1:nx,j),nx/2)
      ffxy(1:nx,ny+j) = cshift(ffxy(1:nx,ny-(j-1)),nx/2)
    end do
    do i=1, nh
      ffx(1-i,:) = ffx(nx-(i-1),:)
      ffy(1-i,:) = ffy(nx-(i-1),:)
      ffxy(1-i,:) = ffxy(nx-(i-1),:)
    end do
    do i=1, nh+1
      ffx(nx+i,:) = ffx(1+(i-1),:)
      ffy(nx+i,:) = ffy(1+(i-1),:)
      ffxy(nx+i,:) = ffxy(1+(i-1),:)
    end do

  end subroutine interpolate_setd

  subroutine interpolate_setdx(fx)
    implicit none

    real(kind=dp), dimension(:,:), intent(in) :: fx

    integer(kind=i4b) :: i, j

    ffx(1:nx,1:ny) = fx
! direction of d/dx is reversed beyond poles
    do j=1, nh+1
      ffx(1:nx,1-j) = -cshift(ffx(1:nx,j),nx/2)
      ffx(1:nx,ny+j) = -cshift(ffx(1:nx,ny-(j-1)),nx/2)
    end do
    do i=1, nh
      ffx(1-i,:) = fx(nx-(i-1),:)
    end do
    do i=1, nh+1
      ffx(nx+i,:) = fx(1+(i-1),:)
    end do
    ffy(:,:) = 0.d0
    ffxy(:,:) = 0.d0

  end subroutine interpolate_setdx

  subroutine interpolate_diff()
    implicit none

    integer(kind=i4b) :: i, j

    forall(i=1:nx, j=1:ny)
      ffxl(i,j) = ff(i+1,j) - ff(i,j)
      ffyl(i,j) = ff(i,j+1) - ff(i,j)
    end forall

  end subroutine interpolate_diff

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

  subroutine quasimonotone_filter(fi,i0,j0)
    implicit none

    real(kind=dp), intent(inout) :: fi
    integer(kind=i4b), intent(in) :: i0, j0
    logical :: lmono

! apply filter to a cell that should be monotonic
    lmono = .false.
    if (n/2>=2) then
      lmono = .not. ( &! Nair et al. 1999
        ((ffxl(i0-2,j0)  *ffxl(i0-1,j0)  >0.0_dp).and. &
         (ffxl(i0-1,j0)  *ffxl(i0+1,j0)  <0.0_dp).and. &
         (ffxl(i0+1,j0)  *ffxl(i0+2,j0)  >0.0_dp)).or. &
        ((ffxl(i0-2,j0+1)*ffxl(i0-1,j0+1)>0.0_dp).and. &
         (ffxl(i0-1,j0+1)*ffxl(i0+1,j0+1)<0.0_dp).and. &
         (ffxl(i0+1,j0+1)*ffxl(i0+2,j0+1)>0.0_dp)).or. &
        ((ffyl(i0,  j0-2)*ffyl(i0,  j0-1)>0.0_dp).and. &
         (ffyl(i0,  j0-1)*ffyl(i0,  j0+1)<0.0_dp).and. &
         (ffyl(i0,  j0+1)*ffyl(i0,  j0+2)>0.0_dp)).or. &
        ((ffyl(i0+1,j0-2)*ffyl(i0+1,j0-1)>0.0_dp).and. &
         (ffyl(i0+1,j0-1)*ffyl(i0+1,j0+1)<0.0_dp).and. &
         (ffyl(i0+1,j0+1)*ffyl(i0+1,j0+2)>0.0_dp)) )
    else
      lmono = & ! Sun et al. 1996
        (ffxl(i0-1,j0)*ffxl(i0+1,j0)>=0.0_dp).and. &
        (ffxl(i0-1,j0+1)*ffxl(i0+1,j0+1)>=0.0_dp).and. &
        (ffyl(i0,j0-1)*ffyl(i0,j0+1)>=0.0_dp).and. &
        (ffyl(i0+1,j0-1)*ff(i0+1,j0+1)>=0.0_dp)
    end if
    if (lmono) then
! Bermejo and Staniforth 1992
      fi = max(min(fi,maxval(ff(i0:i0+1,j0:j0+1))),minval(ff(i0:i0+1,j0:j0+1)))
    end if

  end subroutine quasimonotone_filter

  subroutine quasimonotone_filter_bicubic(fi,i0,j0)
    implicit none

    real(kind=dp), intent(inout) :: fi
    integer(kind=i4b), intent(in) :: i0, j0
    logical :: lmono

! apply filter to a cell that should be monotonic
    lmono = .false.
    if (n/2>=2) then
      lmono = .not. ( &! Nair et al. 1999
        ((ffxl(i0-2,j0)  *ffxl(i0-1,j0)  >0.0_dp).and. &
         (ffxl(i0-1,j0)  *ffxl(i0+1,j0)  <0.0_dp).and. &
         (ffxl(i0+1,j0)  *ffxl(i0+2,j0)  >0.0_dp).and. &
         (ffxl(i0-1,j0)  *ffx( i0,  j0)  >0.0_dp).and. &
         (ffxl(i0+1,j0)  *ffx( i0+1,j0)  >0.0_dp)).or. &
        ((ffxl(i0-2,j0+1)*ffxl(i0-1,j0+1)>0.0_dp).and. &
         (ffxl(i0-1,j0+1)*ffxl(i0+1,j0+1)<0.0_dp).and. &
         (ffxl(i0+1,j0+1)*ffxl(i0+2,j0+1)>0.0_dp).and. &
         (ffxl(i0-1,j0+1)*ffx( i0,  j0+1)>0.0_dp).and. &
         (ffxl(i0+1,j0+1)*ffx( i0+1,j0+1)>0.0_dp)).or. &
        ((ffyl(i0,  j0-2)*ffyl(i0,  j0-1)>0.0_dp).and. &
         (ffyl(i0,  j0-1)*ffyl(i0,  j0+1)<0.0_dp).and. &
         (ffyl(i0,  j0+1)*ffyl(i0,  j0+2)>0.0_dp).and. &
         (ffyl(i0,  j0-1)*ffy( i0,  j0)  >0.0_dp).and. &
         (ffyl(i0,  j0+1)*ffy( i0,  j0+1)>0.0_dp)).or. &
        ((ffyl(i0+1,j0-2)*ffyl(i0+1,j0-1)>0.0_dp).and. &
         (ffyl(i0+1,j0-1)*ffyl(i0+1,j0+1)<0.0_dp).and. &
         (ffyl(i0+1,j0+1)*ffyl(i0+1,j0+2)>0.0_dp).and. &
         (ffyl(i0+1,j0-1)*ffy( i0+1,j0)  >0.0_dp).and. &
         (ffyl(i0+1,j0+1)*ffy( i0+1,j0+1)>0.0_dp)) )
    else
      lmono = .not. ( &! Sun et al. 1996
        ((ffxl(i0-1,j0)  *ffxl(i0+1,j0)  <0.0_dp).and. &
         (ffxl(i0-1,j0)  *ffx( i0,  j0)  >0.0_dp).and. &
         (ffxl(i0+1,j0)  *ffx( i0+1,j0)  >0.0_dp)).or. &
        ((ffxl(i0-1,j0+1)*ffxl(i0+1,j0+1)<0.0_dp).and. &
         (ffxl(i0-1,j0+1)*ffx( i0,  j0+1)>0.0_dp).and. &
         (ffxl(i0+1,j0+1)*ffx( i0+1,j0+1)>0.0_dp)).or. &
        ((ffyl(i0,  j0-1)*ffyl(i0,  j0+1)<0.0_dp).and. &
         (ffyl(i0,  j0-1)*ffy( i0,  j0)  >0.0_dp).and. &
         (ffyl(i0,  j0+1)*ffy( i0,  j0+1)>0.0_dp)).or. &
        ((ffyl(i0+1,j0-1)*ffyl(i0+1,j0+1)<0.0_dp).and. &
         (ffyl(i0+1,j0-1)*ffy( i0+1,j0)  >0.0_dp).and. &
         (ffyl(i0+1,j0+1)*ffy( i0+1,j0+1)>0.0_dp)) )
    end if
    if (lmono) then
! Bermejo and Staniforth 1992
      fi = max(min(fi,maxval(ff(i0:i0+1,j0:j0+1))),minval(ff(i0:i0+1,j0:j0+1)))
    end if

  end subroutine quasimonotone_filter_bicubic

end module interpolate_module
