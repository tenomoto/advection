module bicubic_module

! bicubic interpolation
!  call bcucof to calculate coefficents
!  function bcuint returns the interpolated value
! Source: Numerical Recipes

! Author: T. Enomoto
! History: 3 March 2004

  use constant_module, only : i4b, dp
  private

  real(kind=dp), dimension(4,4), private :: c
  public :: bcucof, bcuint, bcuintp

contains

  subroutine bcucof(f, fx, fy, fxy, dx, dy)
    implicit none

    real(kind=dp), dimension(4), intent(in) :: f, fx, fy, fxy
    real(kind=dp), intent(in) :: dx, dy

    integer(kind=i4b) :: i, j, l
    real(kind=dp), dimension(16) :: x, cl
    real(kind=dp), dimension(16,16) :: wt
    real(kind=dp) :: dxdy
    data  wt/&
      1, 0,-3, 2, 0, 0, 0, 0,-3, 0, 9,-6, 2, 0,-6, 4, &
      0, 0, 0, 0, 0, 0, 0, 0, 3, 0,-9, 6,-2, 0, 6,-4, &
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9,-6, 0, 0,-6, 4, &
      0, 0, 3,-2, 0, 0, 0, 0, 0, 0,-9, 6, 0, 0, 6,-4, &
      0, 0, 0, 0, 1, 0,-3, 2,-2, 0, 6,-4, 1, 0,-3, 2, &
      0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 3,-2, 1, 0,-3, 2, &
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 2, 0, 0, 3,-2, &
      0, 0, 0, 0, 0, 0, 3,-2, 0, 0,-6, 4, 0, 0, 3,-2, &
      0, 1,-2, 1, 0, 0, 0, 0, 0,-3, 6,-3, 0, 2,-4, 2, &
      0, 0, 0, 0, 0, 0, 0, 0, 0, 3,-6, 3, 0,-2, 4,-2, &
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 2,-2, &
      0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 3,-3, 0, 0,-2, 2, &
      0, 0, 0, 0, 0, 1,-2, 1, 0,-2, 4,-2, 0, 1,-2, 1, &
      0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 2,-1, 0, 1,-2, 1, &
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0,-1, 1, &
      0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 2,-2, 0, 0,-1, 1/

    dxdy = dx*dy
    do i=1, 4
      x(i) = f(i)
      x(i+4) = fx(i)*dx
      x(i+8) = fy(i)*dy
      x(i+12) = fxy(i)*dxdy
    end do
    do i=1, 16
      cl(i) = sum(wt(i,:)*x)
    end do
    l = 0
    do i=1, 4
      do j=1, 4
        l = l + 1
        c(i,j) = cl(l)
      end do
    end do

  end subroutine bcucof

  function bcuint(t,u) result (fi)
    implicit none

    real(kind=dp), intent(in) :: t, u
    real(kind=dp) :: fi

    integer(kind=i4b) :: i

    fi = 0.0_dp
    do i=4,1,-1
      fi = t*fi+((c(i,4)*u+c(i,3))*u+c(i,2))*u+c(i,1)
    end do

  end function bcuint
      
  function bcuintp(t,u) result (fi)
    implicit none

    real(kind=dp), intent(in) :: t, u
    real(kind=dp) :: fi

    integer(kind=i4b) :: i
    real(kind=dp), dimension(4) :: tt

    tt = (/1.0_dp-t, 1.0_dp-t, t, t/)
    fi = 0.0_dp
    do i=4,1,-1
      fi = tt(i)*fi+((c(i,4)*u+c(i,3))*u+c(i,2))*u+c(i,1)
    end do

  end function bcuintp
      
end module bicubic_module
