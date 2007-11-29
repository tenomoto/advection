module cubicspline_module

! Cubic spline interpolation
! Source: Wikipedia, www.cubic.org

! Author: Takeshi Enomoto

! History:
! 2007-11-21: First version

  use constant_module, only : i4b, dp
  implicit none
  private

  integer(kind=i4b), parameter, public :: &
    cubicspline_hermite = 1, cubicspline_bezier = 2
  integer(kind=i4b), private :: m = cubicspline_hermite

  integer(kind=i4b), dimension(4,4,2) :: coeff
  data  coeff / 2, -2,  1,  1,  & ! hermite
               -3,  3, -2, -1,  &
                0,  0,  1,  0,  &
                1,  0,  0,  0,  &
               -1,  3, -3,  1,  & ! bezier
                3, -6,  3,  0,  &
               -3,  3,  0,  0,  &
                1,  0,  0,  0/
  public :: cubicspline_interpolate

contains

  subroutine cubicspline_interpolate(t, dx, f, fi, method)
    implicit none

    integer(kind=i4b), parameter :: n = 4

    real(kind=dp), intent(in) :: t, dx
    real(kind=dp), dimension(n), intent(in) :: f ! f1, f2, df1, df2
    real(kind=dp), intent(out) :: fi
    integer(kind=i4b), optional, intent(in) :: method

    real(kind=dp), dimension(n) :: tt, c, ff
    integer(kind=i4b) :: i

    if (present(method)) then
      m = method
    end if

    tt = (/t*t*t, t*t, t, 1.0_dp/)
    do i=1, n
      c(i) = sum(tt*coeff(i,:,m))
    end do
		ff(1:2) = f(1:2)
    ff(3:4) = f(3:4)*dx
    fi = sum(c*ff)
    
  end subroutine cubicspline_interpolate

end module cubicspline_module
