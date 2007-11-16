module alf_module

! calculates normalised associated Legendre functions

! Source: Based on Swartztrauber (2003)
! Author: T. Enomoto
! Usage:
!   Calculates the values of normalized associated Legendre polynomials
!   at latitudes lat
! NB:
!  alf_calc takes latitudes, alf_calc_old takes glat = sin(lat)
!  normalised to 1. factor (-1)**m is not included.
! History:
!   TE 2004-01-09: Fixed bug in zero clear (does not affect AFES)
!                  Added alf_calc_old for testing purpose
!   TE 2003-12-02: Removed AFES dependent code
!   TE 2003-05-30: Some comments are added.
!   TE 2003-05-23: Created.

	use constant_module, only: i4b, dp
	implicit none
	private

	real(kind=dp), public, dimension(:,:,:), allocatable :: pnm

	public :: alf_calc, alf_calc_old, alf_clean, alf_enm

contains

	function alf_enm(n, m) result(e)
		implicit none

		integer(kind=i4b), intent(in) :: n, m
		real(kind=dp) :: e, nn, mm

		nn = n*n
		mm = m*m
		e =  sqrt((nn - mm)/(4.0_dp*nn - 1.0_dp))

	end function alf_enm

	subroutine alf_calc_old(glat, mmax)
		implicit none

		real(kind=dp), dimension(:), intent(in) :: glat
		integer(kind=i4b), intent(in) :: mmax

		integer(kind=i4b) :: n, m, j, jmax, jmax2, ierr

		jmax = size(glat)
		jmax2 = jmax/2
		allocate(pnm(jmax2, 0:mmax, 0:mmax), stat = ierr)
		if (ierr > 0) then
			print *, "Allocation err in alf_calc"
			stop
		end if
		pnm = 0.0_dp

		pnm(:,0,0) = sqrt(0.5_dp)
		pnm(:,1,0) = sqrt(1.5_dp)*glat(1:jmax2)

		do n=2, mmax
			pnm(:,n,0) = (glat(1:jmax2)*pnm(:,n-1,0)-alf_enm(n-1,0)*pnm(:,n-2,0))/alf_enm(n,0)
		end do
		do m=1, mmax
			pnm(:,m,m) = sqrt(0.5_dp*(2.0_dp*m+1.0_dp)/m*(1.0_dp-glat(1:jmax2)**2))*pnm(:,m-1,m-1)
			pnm(:,m+1,m) = sqrt(2.0_dp*m+3.0_dp)*glat(1:jmax2)*pnm(:,m,m)
			do n=m+2, mmax
				pnm(:,n,m) = (glat(1:jmax2)*pnm(:,n-1,m) - alf_enm(n-1,m)*pnm(:,n-2,m)) / alf_enm(n,m)
			end do
		end do

	end subroutine alf_calc_old
	
	subroutine alf_calc(lat,mmax)

		implicit none

		real(kind=dp), dimension(:), intent(in) :: lat
		integer(kind=i4b), intent(in) :: mmax

		integer(kind=i4b) :: j, m, n, l, k, n2, nmod, jmax, ierr
		real(kind=dp) :: sqrt_nn1_rev, theta, pi2
		real(kind=dp), dimension(:,:), allocatable :: ank
		
		jmax = size(lat)
		allocate(pnm(jmax/2, 0:mmax, 0:mmax), stat = ierr)
		if (ierr > 0) then
			print *, "Allocation err in alf_calc"
			stop
		end if
		pnm = 0.0_dp
		allocate(ank(2:mmax, 0:mmax/2))
		ank = 0.0_dp

! calculate fourier coefficients for Pn
		ank(2,1) = 0.75_dp*sqrt(2.5_dp) 
		do n=3, mmax
			ank(n,n/2) = &
				sqrt(1.0_dp-1.0_dp/(4.0_dp*n*n)) * ank(n-1,(n-1)/2)
		end do
		do n=2, mmax
			n2  = n/2
			do k=1, n2
 		  	l = 2*k
		  	ank(n,n2-k) = (l-1.0_dp)*(2.0_dp*n-l+2.0_dp)/&
					(l*(2.0_dp*n-l+1.0_dp)) * ank(n,n2-k+1)
			end do
			if (n==n2*2) then
				ank(n,0) = 0.5_dp*ank(n,0)
			end if
		end do

! calculate Pnm
		pi2 = 0.5_dp * acos(-1.0_dp)

		do j=1, jmax/2

			theta = pi2 - lat(j)
			pnm(:, 0, 0) = 1.0_dp/sqrt(2.0_dp)

! Pmm and Pm,n=m+1
			do m=1, mmax
				pnm(j, m, m) = sqrt(1.0_dp + 0.5_dp/m) * sin(theta) * pnm(j, m-1, m-1)
			end do
			do m=1, mmax-1
				pnm(j, m+1, m) = sqrt(2.0_dp * m + 3.0_dp) * cos(theta) * pnm(j, m, m)
			end do

! m = 0
			pnm(j, 1, 0) = sqrt(1.5_dp)*cos(theta)
			do n=2, mmax
				nmod = n - n/2*2
				do l=0, n/2
					k = 2*l + nmod ! n even: k=2*l, n odd: k=2*l+1
					pnm(j, n, 0) = pnm(j, n, 0) + ank(n,l)*cos(k*theta)
				end do
			end do

! m = 1
			do n=3, mmax
				nmod = n - n/2*2
				sqrt_nn1_rev = 1.0_dp/sqrt(n*(n+1.0_dp))
				do l=0, n/2
					k = 2*l + nmod ! n even: k=2*l, n odd: k=2*l+1
					pnm(j, n, 1) = pnm(j, n, 1) + ank(n,l)*k*sqrt_nn1_rev*sin(k*theta)
				end do
			end do

! m > 1
			do m=2, mmax-2
				do n=m+2, mmax
					pnm(j, n, m) = ( &
						sqrt((2.0_dp*n+1.0_dp)/(2.0_dp*n-3.0_dp))* &
							( sqrt((n+m-2.0_dp)*(n+m-3.0_dp)) * pnm(j, n-2, m-2) + &
						 	   sqrt((n-m)*(n-m-1.0_dp)) * pnm(j, n-2, m) ) &
						- sqrt((n-m+1.0_dp)*(n-m+2.0_dp)) * pnm(j, n, m-2) &
						) / sqrt((n+m-1.0_dp)*(n+m))
				end do
			end do

		end do! j

		deallocate(ank)

	end subroutine alf_calc

	subroutine alf_clean()
		implicit none

		deallocate(pnm)

	end subroutine alf_clean

end module alf_module
