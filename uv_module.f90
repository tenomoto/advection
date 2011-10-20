module uv_module
  use kind_module, only: dp, i4b
  implicit none
  private

  public :: uv_sbody, uv_nodiv, uv_div

contains

  subroutine uv_sbody(lon,lat,gu,gv)
    use planet_module, only: d=>day_in_sec
    use math_module, only: pi2=>math_pi2, deg2rad=>math_deg2rad
    implicit none

    real(kind=dp), dimension(:), intent(in) :: lon, lat
    real(kind=dp), dimension(:,:), intent(inout) :: gu, gv

    real(kind=dp), parameter :: &
      x0 = 0.0_dp, y0 = 45.0_dp, period = 20.0_dp ! Rotation

    real(kind=dp) :: lon0, lat0, omg, l
    integer(kind=i4b) :: i, j, nx, ny

    nx = size(lon)
    ny = size(lat)

    lon0 = x0*deg2rad
    lat0 = y0*deg2rad
    omg = pi2/(period*d)
    do j=1, ny
      do i=1, nx
        gu(i,j) = omg*(cos(lat(j))*sin(lat0)-cos(lon(i)-lon0)*sin(lat(j))*cos(lat0))
        gv(i,j) = omg*(sin(lon(i)-lon0)*cos(lat0))
      end do
    end do

  end subroutine uv_sbody

  subroutine uv_nodiv(t,lon,lat,gu,gv)
    use math_module, only: pi=>math_pi, pi2=>math_pi2
    implicit none

    real(kind=dp), intent(in) :: t
    real(kind=dp), dimension(:), intent(in) :: lon, lat
    real(kind=dp), dimension(:,:), intent(inout) :: gu, gv

    real(kind=dp), parameter :: t1 = 5.0_dp, kappa = 2.0_dp

    real(kind=dp) :: lambda1, pt1, ptt1, kcosptt1
    integer(kind=i4b) :: i, j, nx, ny

    nx = size(lon)
    ny = size(lat)

    pt1 = pi2/t1
    ptt1 = pt1*t
    kcosptt1 = kappa*cos(pi*t/t1)
    do j=1, ny
      do i=1, nx
        lambda1 = lon(i) - ptt1
        gu(i,j) = sin(lambda1)**2*sin(2.0_dp*lat(j))*kcosptt1 + pt1*cos(lat(j))
        gv(i,j) = sin(2.0_dp*lambda1)*cos(lat(j))*kcosptt1
      end do
    end do

  end subroutine uv_nodiv

  subroutine uv_div(t,lon,lat,gu,gv)
    use math_module, only: pi=>math_pi, pi2=>math_pi2
    implicit none

    real(kind=dp), intent(in) :: t
    real(kind=dp), dimension(:), intent(in) :: lon, lat
    real(kind=dp), dimension(:,:), intent(inout) :: gu, gv

    real(kind=dp), parameter :: t1 = 5.0_dp, kappa = 1.0_dp

    real(kind=dp) :: lambda1, pt1, ptt1, kcosptt1
    integer(kind=i4b) :: i, j, nx, ny

    nx = size(lon)
    ny = size(lat)

    pt1 = pi2/t1
    ptt1 = pt1*t
    kcosptt1 = kappa*cos(pi*t/t1)
    do j=1, ny
      do i=1, nx
        lambda1 = lon(i) - ptt1
        gu(i,j) = -sin(0.5_dp*lambda1)**2*sin(2.0_dp*lat(j))*cos(lat(j))**2*kcosptt1 + pt1*cos(lat(j))
        gv(i,j) = 0.5_dp*sin(lambda1)*cos(lat(j))**3*kcosptt1
      end do
    end do

  end subroutine uv_div

end module uv_module
