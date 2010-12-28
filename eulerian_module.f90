module eulerian_module
  use kind_module, only: i4b, dp
  use planet_module, only: a=>planet_radius
  use grid_module, only: &
    nlon, nlat, ntrunc, &
    gu, gv, gphi, sphi, sphi_old, lon, lat, coslatr
  use time_module, only: nstep, hstep, deltat, ifile, hfile
  use legendre_transform_module, only: &
    legendre_analysis, legendre_synthesis, &
    legendre_synthesis_dlon, legendre_synthesis_dlat
  use io_module, only: io_save
  implicit none
  private

  real(kind=dp), dimension(:,:), allocatable, private :: dgphi

  integer(kind=i4b), private :: nsave = 0
  complex(kind=dp), dimension(:,:), allocatable, private :: sphi1

  private :: update
  public :: eulerian_init, eulerian_timeint, eulerian_clean

contains

  subroutine eulerian_init()
    implicit none

    integer(kind=i4b) :: j

    allocate(dgphi(nlon,nlat), sphi1(0:ntrunc,0:ntrunc))

    print *, "Saving initial value"
    call io_save(ifile, 1, gphi, "replace")
    call io_save(ifile, 2, gu, "old")
    call io_save(ifile, 3, gv, "old")
    call legendre_synthesis(sphi,gphi)
    print *, "step=0 t=0"
    print *, "Saving step=0"
    call io_save(hfile, 1, gphi, "replace")
    call io_save(hfile, 2, gu, "old")
    call io_save(hfile, 3, gv, "old")
    nsave = 1
    print *, "step=1/2", " t=", real(0.5*deltat)
    call update(0.0d0*deltat,0.5d0*deltat)
    print *, "step=1", " t=", real(deltat)
    call update(0.5d0*deltat,deltat)
    if (hstep==1) then
      call io_save(hfile, 3*nsave+1, gphi, "old")
      call io_save(hfile, 3*nsave+2, gu, "old")
      call io_save(hfile, 3*nsave+3, gv, "old")
      nsave = nsave + 1
    end if

  end subroutine eulerian_init

  subroutine eulerian_clean()

    deallocate(dgphi, sphi1)

  end subroutine eulerian_clean

  subroutine eulerian_timeint()
    implicit none

    integer(kind=i4b) :: i, j

    do i=2, nstep
      print *, "step=", i, " t=", real(i*deltat)
      call update((i-1)*deltat,2.0d0*deltat)
      if (mod(i,hstep)==0) then
        print *, "Saving step=", i
        call io_save(hfile, 3*nsave+1, gphi, "old")
        call io_save(hfile, 3*nsave+2, gu, "old")
        call io_save(hfile, 3*nsave+3, gv, "old")
        nsave = nsave + 1
      end if
    end do

  end subroutine eulerian_timeint

  subroutine update(t,dt)
    use time_module, only: etf, kappa
    use uv_module, only: uv_nodiv
    implicit none

    real(kind=dp), intent(in) :: t, dt

    real(kind=dp) :: knt
    integer(kind=i4b) :: i, j, m, n

    call uv_nodiv(t,lon,lat,gu,gv)
    call legendre_synthesis(sphi_old, gphi)
! dF/dlon
    call legendre_synthesis_dlon(sphi, dgphi)
    do j=1, nlat
      dgphi(:,j) =  dgphi(:,j)*coslatr(j)
    end do
    gphi = gphi - dt*gu*dgphi
! cos(lat)dF/dlat
    call legendre_synthesis_dlat(sphi, dgphi)
    do j=1, nlat
      dgphi(:,j) =  dgphi(:,j)*coslatr(j)
    end do
    gphi = gphi - dt*gv*dgphi
    call legendre_analysis(gphi, sphi1)
! dissipation
    do n=1, ntrunc
      knt = kappa*dt*(n*(n+1.0_dp))**2
      do m=0, n
        sphi1(n,m) = (1.0_dp-knt)*sphi1(n,m)
      end do
    end do
! update
    do m=0, ntrunc
        sphi_old(m:ntrunc,m) = sphi(m:ntrunc,m) + &
          etf * (sphi_old(m:ntrunc,m)-2.0_dp*sphi(m:ntrunc,m)+sphi1(m:ntrunc,m))
    end do
    do m=0, ntrunc
      do n=m, ntrunc
        sphi(n,m) = sphi1(n,m)
      end do
    end do

  end subroutine update

end module eulerian_module
  
