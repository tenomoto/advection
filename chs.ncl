load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
n = 10
pi = acos(-1.)
function func(x)
begin
	return sin(x)
end

function dfunc(x)
begin
	return 0.5*pi*cos(x)
end

begin
	p = new(4, "float")
	p(0) = func(0.)
	p(1) = func(0.5*pi)
	p(2) = dfunc(0.)
	p(3) = dfunc(0.5*pi)

	t = fspan(0.,1.,n)
	h = new((/4,n/), "float")
	h(0,:) = 2*t^3-3*t^2+1
	h(1,:) = -2*t^3+3*t^2
  h(2,:) = t^3-2*t^2+t
	h(3,:) = t^3-t^2
	f = new((/2,n/), "float")
	do i=0, n-1
		x = (i*1.)/(n-1.)*0.5*pi
	  f(0,i) = func(x)
		f(1,i) = h(0,i)*p(0)+h(1,i)*p(1)+h(2,i)*p(2)+h(3,i)*p(3)
		print(f(0,i)+" "+f(1,i))
	end do

	wks = gsn_open_wks("x11", "chs")
	res = True
	res@xyDashPattern = 0
	res@xyMarkLineMode = "MarkLines"
	res@xyMarkers      = (/4,5/)
	res@xyMarkerColors = (/"blue","red"/)

	plot = gsn_xy(wks, t, f, res)
end
