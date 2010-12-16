load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

  nt = 39
  nx = 3*(nt+1)
  ny = nx/2

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

  wks = gsn_open_wks("x11", "init")

  res = True
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnLevelSpacingF = 10.
  res@cnMinLevelValF = 10.
  res@cnMaxLevelValF = 90.

  resv = True
; resv@gsnPolar = "NH"

  x = fbindirread("init.dat", 0, (/ny, nx/), "double")
  print("max phi="+max(x))
  print("min phi="+min(x))
  plot = gsn_csm_contour_map(wks, x, res)
  u = fbindirread("init.dat", 1, (/ny, nx/), "double")
  v = fbindirread("init.dat", 2, (/ny, nx/), "double")
; plot = gsn_csm_contour_map_polar(wks, u, resv)
; plot = gsn_csm_contour_map_polar(wks, v, resv)
  plot = gsn_csm_streamline_map_ce(wks, u, v, resv)

end

