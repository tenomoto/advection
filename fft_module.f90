module fft_module
  use kind_module, only: i4b, i8b, dp
  implicit none
  include "fftw3.f"
  private

  integer(kind=i8b), private :: plan_forward, plan_backward

  public :: fft_init, fft_clean, fft_analysis, fft_synthesis

contains

  subroutine fft_init(g, w)
    implicit none

    real(kind=dp), dimension(:,:), intent(in) :: g 
    complex(kind=dp), dimension(0:,:), intent(in) :: w 
    integer(kind=i4b) :: i, j

    i = size(g,1)
    j = size(g,2)

    call dfftw_plan_dft_r2c_1d(plan_forward, i, g(:,1), w(:,1), FFTW_ESTIMATE)
    call dfftw_plan_dft_c2r_1d(plan_backward, i, w(:,1), g(:,1), FFTW_ESTIMATE)

  end subroutine fft_init

  subroutine fft_clean()
    implicit none

    call dfftw_destroy_plan(plan_forward)
    call dfftw_destroy_plan(plan_backward)

  end subroutine fft_clean

  subroutine fft_analysis(g, w)
    implicit none

    real(kind=dp), dimension(:,:), intent(in) :: g 
    complex(kind=dp), dimension(0:,:), intent(inout) :: w 
    integer(kind=i4b) :: j

    do j=1, size(g,2)
      call dfftw_execute_dft_r2c(plan_forward, g(:,j), w(:,j))
    end do
    w = w/size(g,1)

  end subroutine fft_analysis

  subroutine fft_synthesis(w,g)
    implicit none

    complex(kind=dp), dimension(0:,:), intent(in) :: w 
    real(kind=dp), dimension(:,:), intent(inout) :: g 
    integer(kind=i4b) :: j

    do j=1, size(g,2)
      call dfftw_execute_dft_c2r(plan_backward, w(:,j), g(:,j))
    end do

  end subroutine fft_synthesis

end module fft_module
