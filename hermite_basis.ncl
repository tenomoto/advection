load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
n = 100
begin
  wks = gsn_open_wks("x11", "hermite_basis")
  t = fspan(0., 1., n)
  h = new((/4,n/), "float")
  h(0,:) = 2*t^3-3*t^2+1
  h(1,:) = -2*t^3+3*t^2
  h(2,:) = t^3-2*t^2+t
  h(3,:) = t^3-t^2
  plot = gsn_xy(wks, t, h, False)
end
