load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

  ntrunc = 39
  nx = 3*(ntrunc+1)
  ny = nx/2
  nt =  20 
  dt = 5./nt
  tstep =1 
  dev = "x11"
  fname = "hist.dat"
;  fname = "init.dat"
  a = 6.371e6
  pi = acos(-1.)
  deg2rad = pi/180.

  lon = fspan(0, 360.-360./nx, nx)
  lon!0 = "lon"
  lon&lon = lon
  lon@units = "degrees_east"
  gau_info = gaus(ny/2)
  lat = doubletofloat(gau_info(ny-1:0,0))
  lat!0 = "lat"
  lat&lat = lat
  lat@units = "degrees_north"

  u = new((/ny, nx/), "double")
  u!0 = "lat"
  u&lat = lat
  u!1 = "lon"
  u&lon = lon
  v = u

  wks = gsn_open_wks(dev, fname)

  res = True
  res@gsnMaximize = True
  res@mpCenterLonF = 180.
  res@mpFillOn = False
  res@mpOutlineBoundarySets = "NoBoundaries"
  res@tiMainString = fname
  res@vcMinDistanceF = 0.017
  res@vcRefMagnitudeF = 2.0
  res@vcRefLengthF = 0.01

  do i=0, nt, tstep
    res@gsnCenterString = "day="+i*dt
    u = fbindirread(fname, 3*i+1, (/ny, nx/), "double")
    v = fbindirread(fname, 3*i+2, (/ny, nx/), "double")
    print("day="+(i*dt)+" umax="+max(u)+" umin="+min(u)+" vmax="+max(v)+" vmin="+min(v))
    res@gsnLeftString = "u "+max(u)+"/"+min(u)
    res@gsnRightString = "v "+max(v)+"/"+min(v)
   plot = gsn_csm_vector_map(wks, u, v, res)
  end do

end

