!! provides experiment parameters
module parameter_module
	use constant_module, only : i4b, dp
	implicit none

	integer(kind=i4b), public ::  ntrunc, nlon, nlat, nstep, hstep
	real(kind=dp), public :: deltat

	public :: set_parameters

contains

	subroutine set_parameters()
		implicit none

		print *, "Enter truncation wavenumber:"
		read *, ntrunc

		nlon = (ntrunc+1)*3
		nlat = nlon/2
		print *, "Horizontal resolution:"
		print *, "ntrunc =", ntrunc, " nlon =", nlon, " nlat =", nlat

		print *, "Enter number of steps:"
		read *, nstep
		print *, "Number of steps = ", nstep

		print *, "Enter delta t:"
		read *, deltat
		print *, "Delta t =", deltat

		print *, "Enter history interval in steps:"
		read *, hstep
		print *, "History interval = ", hstep

	end subroutine set_parameters


end module parameter_module
