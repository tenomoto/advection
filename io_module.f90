module io_module

  use constant_module, only: i4b, dp
  private

  public save_data

contains

  subroutine save_data(f, r, g, s)
    implicit none

    character(len=*), intent(in) :: f, s
    integer(kind=i4b), intent(in) :: r
    real(kind=dp), dimension(:,:), intent(in) :: g

    integer(kind=i4b), parameter :: u = 31

    open(unit=u, file=f, access="direct", status=s, &
      action="write", recl=size(g,1)*size(g,2)*8)
    write(unit=u, rec=r) g
    close(unit=u)

  end subroutine save_data

end module io_module
