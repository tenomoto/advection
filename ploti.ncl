load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

  ntrunc = 39
  nx = 120
  ny =  60
  nt =  20 
  dt = 5./nt
  tstep =1 
  dev = "x11"

  fname = "init.dat"
;  a = 6.371e6

  lon = fspan(0, 360.-360./nx, nx)
  lon!0 = "lon"
  lon&lon = lon
  lon@units = "degrees_east"
  gau_info = gaus(ny/2)
  lat = doubletofloat(gau_info(ny-1:0,0))
  wgty = doubletofloat(gau_info(:,1))
  wgtx = 1.0
  lat!0 = "lat"
  lat&lat = lat
  lat@units = "degrees_north"

  x = new((/ny, nx/), "double")
  x!0 = "lat"
  x&lat = lat
  x!1 = "lon"
  x&lon = lon

  wks = gsn_open_wks(dev, "init")

  res = True
  res@gsnMaximize = True
  res@mpCenterLonF = 180.
  res@mpFillOn = False
  res@mpOutlineBoundarySets = "NoBoundaries"
  res@tiMainString = fname
  rese = res
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels = ispan(10,90,10)*0.01
  res@cnMonoLineThickness = False
  res@cnLineThicknesses = (/1.,2.,1.,2.,1.,2.,1.,2.,1.,2.,1.,2./)

  res@gsnCenterString = "init"
  x = fbindirread(fname, 0, (/ny, nx/), "double")
  print("max phi="+max(x)+" min phi="+min(x))
  res@gsnLeftString = "max="+max(x)
  res@gsnRightString = "min="+min(x)
  plot = gsn_csm_contour_map(wks, x, res)

; plot error 
  rese@cnMaxLevelCount = 12
  rese@cnMonoLineThickness = False
  rese@cnLineThicknesses = (/1.,2.,1.,2.,1.,2.,1.,2.,1.,2.,1.,2./)
  rese@gsnCenterString = "Simulated-Analytic"
  rese@gsnDraw = False
  rese@gsnFrame = False
  rese@gsnContourZeroLineThicknessF = 0.0
  rese@gsnContourNegLineDashPattern = 2
  x0 = fbindirread("init1.dat", 0, (/ny, nx/), "double")
  x = x-x0
  print("error max phi="+max(x)+" min phi="+min(x))
  print("global mean error="+wgt_areaave(x,wgty,wgtx,0))
  rese@gsnLeftString = "max="+max(x)
  rese@gsnRightString = "min="+min(x)
  plot = gsn_csm_contour_map(wks, x, rese)
  draw(plot)
  frame(wks)

end

