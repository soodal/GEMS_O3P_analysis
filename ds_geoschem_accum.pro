pro ds_geoschem_accum, data, o3accum, hpa=hpa, height=height


;;----------------------------
;; GEMS vertical column layer for under 300 hPa
;;----------------------------
if keyword_set(hpa) then begin
  o3 = data.o3
  o3size = size(o3, /dim)
  dim1 = o3size[0]
  dim2 = o3size[1]
  gemspres = data.Pressure
  gemslayero3 = fltarr(o3size[0], o3size[1])
  ;gemslayero3[*,*] = !values.f_nan
  for iy=0, dim2-1 do begin
  for ix=0, dim1-1 do begin
    for ilayer=0, 46, 1 do begin
      ;ilevel = ilayer+1
      if gemspres[ilayer, ix, iy] gt hpa $
          and gemspres[ilayer+1, ix, iy] gt hpa then begin
        gemslayero3[ix, iy] = gemslayero3[ix, iy] + o3[ix, iy, ilayer]

      endif else if gemspres[ilayer, ix, iy] gt hpa $
          and gemspres[ilayer+1, ix, iy] le hpa then begin

        gemslayero3[ix, iy] = gemslayero3[ix, iy] + o3[ix, iy, ilayer]*$
          (abs(hpa-gemspres[ilayer, ix, iy]))/$
          (abs(gemspres[ilayer+1, ix, iy]-gemspres[ilayer, ix, iy]))
        ;print, ilayer
        ;stop
        break
      endif
    endfor
  endfor
  endfor
  gemslayero3[where(gemslayero3 lt 0)] = !values.f_nan
  o3accum = gemslayero3
  return
endif

; GEMS vertical column layer for under 10 km
if keyword_set(height) then begin
  o3 = data.o3
  o3size = size(data.o3, /dim)
  altitude = data.Altitude
  nanidx = where(altitude < 0, /null)
  altitude[nanidx] = !values.f_nan
  o3size = size(data.o3, /dim)
  gemslayero3 = fltarr(o3size[0:1])
  dim1 = o3size[0]
  dim2 = o3size[1]
  ;gemslayero3[*,*] = !values.f_nan
  for iy=0, dim2-1 do begin
  for ix=0, dim1-1 do begin
    for ilevel=24, 0, -1 do begin
      ilayer = ilevel-1
      if altitude[ix, iy, ilevel] lt height $
          and altitude[ix, iy, ilevel-1] lt height then begin
        gemslayero3[ix, iy] = gemslayero3[ix, iy] + o3[ix, iy, ilayer]
        ;print, ilayer

      endif else if altitude[ix, iy, ilevel] lt height $
          and altitude[ix, iy, ilevel-1] ge height then begin

        gemslayero3[ix, iy] = gemslayero3[ix, iy] + o3[ix, iy, ilayer]*$
          (height-altitude[ix, iy, ilevel])/$
          (altitude[ix, iy, ilevel-1]-altitude[ix, iy, ilevel])
        ;print, ilayer
        break
      endif
    endfor
  endfor
  endfor
  gemslayero3[where(gemslayero3 le 0)] = !values.f_nan
  o3accum = gemslayero3
  return
endif

end
