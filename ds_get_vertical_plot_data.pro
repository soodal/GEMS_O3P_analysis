pro ds_get_vertical_plot_data, data, lon, lat, pres, idx, $
  z_xdir, x_xdir, y_xdir, z_ydir, x_ydir, y_ydir

sz = size(pres)
dim1 = sz[1]
dim2 = sz[2]
dim3 = sz[3]

if sz[0] ne 3 then begin
  print, 'input data must be 3-dimension'
  return
endif

if n_elements(idx) eq 1 then begin
  site_idx = array_indices(lat, idx)
endif else begin
  site_idx = idx
endelse

; for x-dir

x = lon[*, site_idx[1]]
x = rebin(x, [dim1, dim3])
xnanidx = where(x lt -990, /null)
x[xnanidx] = !values.f_nan

y = reform(pres[*, site_idx[1], *])
ynanidx = where(y lt -990, /null)
y[ynanidx] = !values.f_nan

z = reform(data[*, site_idx[1], *])
znanidx = where(z lt -990, /null)
z[znanidx] = !values.f_nan

x_xdir = x
y_xdir = y
z_xdir = z

; for x-dir

x = reform(lat[site_idx[0], *])
x = rebin(x, [dim2, dim3])
xnanidx = where(x lt -990, /null)
x[xnanidx] = !values.f_nan

y = reform(pres[site_idx[0], *, *])
ynanidx = where(y lt -990, /null)
y[ynanidx] = !values.f_nan

z = reform(data[site_idx[0], *, *])
znanidx = where(z lt -990, /null)
z[znanidx] = !values.f_nan

x_ydir = x
y_ydir = y
z_ydir = z

return
end
