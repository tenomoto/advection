program test_polint
  use polint_module, only : polint
  implicit none

  integer, parameter :: n = 4, m = 10
  real*8, dimension(n) :: xa, ya
  real*8 :: dx, ddx, pi, x, y, dy
  integer :: i

  pi = acos(-1.0d0)
  dx = 0.5d0*pi/(n-1.d0)
  do i=1, n
    xa(i) = dx * (i-1)
    ya(i) = sin(xa(i))
    print *, i, xa(i), ya(i)
  end do

  ddx = dx/(m-1.d0)
  do i=1, m
    x = xa(2) + ddx*(i-1)
    call polint(xa, ya, x, y, dy)
    print *, i, x, sin(x), y
  end do
  
end program test_polint
