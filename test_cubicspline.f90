program test_cubicspline
  use cubicspline_module, only : cubicspline_interpolate
  implicit none

  integer, parameter :: n = 2, m = 10
  real*8, dimension(n) :: xa
  real*8, dimension(2*n) :: ya
  real*8 :: dx, ddx, pi, x, y, t
  integer :: i

  pi = acos(-1.0d0)
  dx = 0.5d0*pi/(2*n-1.d0)
  do i=1, n
    xa(i) = dx * i 
    ya(i) = sin(xa(i))
    ya(i+n) = cos(xa(i))
    print *, i, xa(i), ya(i), ya(i+n)
  end do

  ddx = dx/(m-1.d0)
  do i=1, m
    x = xa(1) + ddx*(i-1)
    t = (x-xa(1))/dx
    call cubicspline_interpolate(t, dx, ya, y)
    print *, i, t, x, sin(x), y
  end do
  
end program test_cubicspline
