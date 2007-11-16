ntrunc = 39
nx = 3*(ntrunc+1)
ny = nx/2
nt =   20
dir = (/"eulerian",    "eulerian", "semi-lag",    "semi-lag", \
        "semi-lag",    "nisl",     "semi-lag",    "semi-lag", "semi-lag"/)
run = (/"2500km",      "2500km",   "2500km_bilin","2500km_fd", \
        "2500km_fdlimg","2500km",   "2500km_sph",  "2500km_fdy", "2500km_polin2limg"/) 
begin
	n = dimsizes(run)
	lon = fspan(0, 360.-360./nx, nx)
	lon!0 = "lon"
	lon&lon = lon
	lon@units = "degrees_east"
	gau_info = gaus(ny/2)
	lat = doubletofloat(gau_info(ny-1:0,0))
	lat!0 = "lat"
	lat&lat = lat
	lat@units = "degrees_north"
	gwgt = doubletofloat(gau_info(ny-1:0,1))

	fname = dir(0)+"/"+"history_L"+run(0)+".dat"
	buf0 = fbindirread(fname, 0, (/ny, nx/), "double")
	print(fname+" min="+min(buf0)+" max="+max(buf0))
 	x0 = wgt_areaave(buf0, gwgt, 1., 0)

	do i=1, n-1
	  fname = dir(i)+"/"+"history_L"+run(i)+".dat"
		buf = fbindirread(fname, 3*nt, (/ny, nx/), "double")
 		x = wgt_areaave(buf-buf0, gwgt, 1., 0)
  	print(fname+" min="+min(buf)+" max="+max(buf))
		d = x/x0*100.
 		dp = wgt_areaave((buf-buf0)>0, gwgt, 1., 0)/x0*100.
 		dn = wgt_areaave((buf-buf0)<0, gwgt, 1., 0)/x0*100.
		print("diff total="+d+" diff pos="+dp+" diff neg="+dn)
	end do
end
