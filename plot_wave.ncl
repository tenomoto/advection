load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

  nt = 19
  nx = 3*(nt+1)
  ny = nx/2
  m = nx/2 + 1 

  wks = gsn_open_wks("x11", "init")
  gsn_define_colormap(wks, "BlueRed")
  res = True
  res@gsnSpreadColors = True
  res@cnFillOn = True
  res@cnLinesOn = False
  res@cnLineLabelsOn = False

  x = fbindirread("w.dat", 0, (/ny, m/), "double")
  y = fbindirread("w.dat", 1, (/ny, m/), "double")
  x = sqrt(x*x + y*y)
  plot = gsn_contour(wks, x, res)

end

