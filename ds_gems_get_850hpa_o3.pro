pro ds_gems_get_850hpa_o3, o3p, o3_850hpa_number_density

pres = o3p.pressure
o3 = o3p.o3

temp = o3p.temperature

;;==============================================================================
sz = size(temp)
;temp_center = fltarr(sz[1], sz[2], sz[3]-1)
;temp_center = !values.f_nan
;for i=0, do begin
  ;temp_center = (temp[*, *, i] + temp[*,*,i+1])/2.
;endfor
;;==============================================================================

temp850 = fltarr(sz[1], sz[2])
for ix = 0, sz[1]-2 do begin
  for iy = 0, sz[1]-2 do begin
    temp850[ix, iy] = interpol(temp[ix, iy, *], pres[ix, iy, *], 850)
  endfor
endfor

temp_center = fltarr(sz[1], sz[2], sz[3]-1)
pres_center = fltarr(sz[1], sz[2], sz[3]-1)

for ix = 0, sz[1]-1 do begin
  for iy = 0, sz[2]-1 do begin
    for iz = 0, sz[3]-2 do begin
      pres_center[ix, iy, iz] = $
        ;exp((alog(pres[ix, iy, iz]) + alog(pres[ix, iy, iz+1]))/2.)
        ((pres[ix, iy, iz] + pres[ix, iy, iz+1])/2.)
    endfor
    temp_center[ix, iy, *] = interpol(temp[ix, iy, *], pres[ix, iy, *], $
      pres_center[ix, iy, *])
  endfor
endfor

o3_boundary = fltarr(sz[1], sz[2], sz[3])

o3_850 = fltarr(sz[1], sz[2])

for ix = 0, sz[1]-2 do begin
  for iy = 0, sz[1]-2 do begin
    o3_850[ix, iy] = interpol(o3[ix, iy, *], pres_center[ix, iy, *], 850)
  endfor
endfor

o3_850hpa_number_density = o3_850*0.4462 ; mmol/m^2
return
end
