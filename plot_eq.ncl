load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
ntrunc = 39
nx = 3*(ntrunc+1)
ny = nx/2
nt =   20
;dir = (/"eulerian",    "eulerian",  "semi-lag",    "semi-lag", \
;        "semi-lag",    "semi-lag",  "semi-lag", "nisl"/)
;run = (/"2500km",      "2500km",    "2500km_bilin","2500km_fd", \
;        "2500km_fdlim","2500km_sph","2500km_fdy"/) 
;dir = (/"eulerian",    "eulerian",  "semi-lag",    "semi-lag", "semi-lag"/)
;run = (/"2500km",      "2500km",    "2500km_bilin","2500km_fd","2500km_fdlimg"/)
dir = (/"eulerian",   "nisl",  "semi-lag",  "semi-lag", "semi-lag"/)
run = (/"2500km",     "2500km","2500km_sph","2500km_fdy", "2500km_polin2limg"/) 
begin
  n = dimsizes(run)
	lon = fspan(-180, 180.-360./nx, nx)
	lon!0 = "lon"
	lon&lon = lon
	lon@units = "degrees_east"
	gau_info = gaus(ny/2)
	lat = gau_info(ny-1:0,0)
	lat!0 = "lat"
	lat&lat = lat
	lat@units = "degrees_north"
	gwgt = gau_info(ny-1:0,1)

	x = new((/n,nx/),"double")
	x!1 = "lon"
	x&lon = lon

	wks = gsn_open_wks("pdf","eq")
	res = True
	res@gsnMaximize = True
	res@trYMinF = -20.
	res@trYMaxF = 125.
	res@trXMinF = -45.
	res@trXMaxF =  45.
	res@xyLineThicknessF = 5
	res@xyMonoDashPattern = True
	res@xyLineColors = (/"black","blue","red","purple","orange"/)

	fname = dir(0)+"/"+"history_L"+run(0)+".dat"
 	buf = fbindirread(fname, 0, (/ny, nx/), "double")
 	x(0,0:nx/2-1) = (/0.5*(buf(ny/2-1,nx/2:nx-1)+buf(ny/2,nx/2:nx-1))/)
 	x(0,nx/2:nx-1) = (/0.5*(buf(ny/2-1,0:nx/2-1)+buf(ny/2,0:nx/2-1))/)
	do i=1, n-1
	 	fname = dir(i)+"/"+"history_L"+run(i)+".dat"
		buf = fbindirread(fname, 3*nt, (/ny, nx/), "double")
 		x(i,0:nx/2-1) = (/0.5*(buf(ny/2-1,nx/2:nx-1)+buf(ny/2,nx/2:nx-1))/)
		x(i,nx/2:nx-1) = (/0.5*(buf(ny/2-1,0:nx/2-1)+buf(ny/2,0:nx/2-1))/)
	end do
 	plot = gsn_csm_xy(wks, lon, x, res)
end
