load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

	ntrunc = 39
	nx = 3*(ntrunc+1)
	ny = nx/2
	nt =   20
	tstep = 5
	dev = "pdf"
	dir = "nisl"
	run = "2500km"
;	cint =   5
;	cmin = -30
;	cmax =  30
;	cint =  2.
;	cmin = -10
;	cmax =  10
;	cint =  0.5
;	cmin = -3
;	cmax =  3

	fname = dir+"/"+"history_L"+run+".dat"
;	fname = "../history.dat"
	a = 6.371e6

	lon = fspan(0, 360.-360./nx, nx)
	lon!0 = "lon"
	lon&lon = lon
	lon@units = "degrees_east"
	gau_info = gaus(ny/2)
	lat = doubletofloat(gau_info(ny-1:0,0))
	lat!0 = "lat"
	lat&lat = lat
	lat@units = "degrees_north"

	x = new((/ny, nx/), "double")
	x!0 = "lat"
	x&lat = lat
	x!1 = "lon"
	x&lon = lon
	u = x
	v = x
;	fname1="semi-lag/history_L1e4km_bicubic.dat"
;	u = fbindirread(fname1, 1, (/ny, nx/), "double")*a
;	v = fbindirread(fname1, 2, (/ny, nx/), "double")*a

	wks = gsn_open_wks(dev, "history_L"+run)

	res = True
	res@gsnMaximize = True
	res@mpProjection = "satellite"
	res@mpCenterLonF =  0.
	res@mpCenterLatF = 45.
	res@tiMainString = dir+" "+run
	rese = res
	res@cnLevelSelectionMode = "ManualLevels"
	res@cnLevelSpacingF = 10.
	res@cnMinLevelValF = 10.
	res@cnMaxLevelValF = 90.
  res@cnMonoLineThickness = False
	res@cnLineThicknesses = (/1.,2.,1.,2.,1.,2.,1.,2.,1.,2.,1.,2./)
;	res@gsnScalarContour = True
;	res@vcMinDistanceF = 0.017
;	res@vcRefMagnitudeF = 5.0
;	res@vcRefLengthF = 0.01

	do i=0, nt, tstep
		res@gsnCenterString = "day="+i
		x = fbindirread(fname, 3*i, (/ny, nx/), "double")
		print("day="+i+" max phi="+max(x)+" min phi="+min(x))
    res@gsnLeftString = "max="+max(x)
    res@gsnRightString = "min="+min(x)
		plot = gsn_csm_contour_map(wks, x, res)
;		plot = gsn_csm_vector_scalar_map(wks, u, v, x, res)
	end do

; plot error 
	rese@cnMaxLevelCount = 12
  rese@cnMonoLineThickness = False
	rese@cnLineThicknesses = (/1.,2.,1.,2.,1.,2.,1.,2.,1.,2.,1.,2./)
	rese@gsnCenterString = "Simulated-Analytic"
;	rese@cnLevelSpacingF = cint
;	rese@cnMinLevelValF =  cmin
;	rese@cnMaxLevelValF =  cmax
	rese@gsnDraw = False
	rese@gsnFrame = False
	x0 = fbindirread(fname, 0, (/ny, nx/), "double")
	x = x-x0
	print("error max phi="+max(x)+" min phi="+min(x))
  rese@gsnLeftString = "max="+max(x)
  rese@gsnRightString = "min="+min(x)
	plot = gsn_csm_contour_map(wks, x, rese)
	plot = ZeroGoneNegDashLineContour(plot)
	draw(plot)
	frame(wks)

end

