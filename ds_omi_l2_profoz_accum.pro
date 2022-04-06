pro ds_omi_l2_profoz_accum, data, o3accum, hpa=pressure_limit, height=height_limit

;if not keyword_set(hpa) then begin
  ;pressure_limit=0
;ENDIF
;if not keyword_set(height_limit) then begin
  ;height_limit=0
;ENDIF


omiecf = data.EffectiveCloudFraction
omialts = data.ProfileLevelAltitude
omipress = data.ProfileLevelPressure

omio3 = data.O3RetrievedProfile
nanidx = where(omio3 lt -1.e30, /null)
omio3size = size(omio3, /dim)
if omio3size[1] eq 30 then begin ; 2021-12-01 OMPROFOZ updated and size changed
  omio3[nanidx] = !values.f_nan
  omio3[*,20:23, *] = !values.f_nan
ENDIF
;omio3 = reform(omio3, [omio3size[0], omio3size[1]*omio3size[2]])
omio3[where(omio3 lt 0, /null)]=!values.f_nan

dim = size(omiecf, /dim)
omio3accum = fltarr(dim[0], dim[1])

;;----------------------------
;; OMI vertical column Layer for under specific pressure
;;----------------------------
if keyword_set(pressure_limit) then begin
  omiecfsize = size(omiecf, /dim)

  for ix=0,dim[0]-1 do begin
  for iy=0,dim[1]-1 do BEGIN

    for ilevel=24, 0, -1 do BEGIN
      ilayer = ilevel - 1
      if omipress[ilevel, ix, iy] gt pressure_limit $
          and omipress[ilevel-1, ix, iy] gt pressure_limit then BEGIN
        omio3accum[ix, iy] = omio3accum[ix, iy] + omio3[ilayer, ix, iy]
        
      ENDIF else if omipress[ilevel, ix, iy] gt pressure_limit $
        and omipress[ilevel-1, ix, iy] le pressure_limit then begin 
        ; upper bound is higher than 300hPa 

        omio3accum[ix, iy] = omio3accum[ix, iy] + omio3[ilayer, ix, iy]*$
          (omipress[ilevel, ix, iy]-pressure_limit)/$
          (omipress[ilevel, ix, iy]-omipress[ilevel-1, ix, iy])
        break
      ENDIF
    endfor
  ENDFOR
  ENDFOR
ENDIF

;;----------------------------
;; OMI vertical column Layer for under specific altitude
;;----------------------------

if keyword_set(height_limit) then begin
  omiecfsize = size(omiecf, /dim)

  altnanidx = where(omialts lt -1.E30, /null)
  omialts[altnanidx] = !values.f_nan

  for ix=0,dim[0]-1 do begin
  for iy=0,dim[1]-1 do BEGIN
    for ilevel=24, 0, -1 do BEGIN
      ilayer = ilevel - 1
      if omialts[ilevel, ix, iy] lt height_limit $
          and omialts[ilevel-1, ix, iy] lt height_limit then BEGIN
        omio3accum[ix, iy] = omio3accum[ix, iy] + omio3[ilayer, ix, iy]
        
      ENDIF else if omialts[ilevel, ix, iy] lt height_limit $
        and omialts[ilevel-1, ix, iy] ge height_limit then begin 

        omio3accum[ix, iy] = omio3accum[ix, iy] + omio3[ilayer, ix, iy]*$
          (height_limit - omialts[ilevel, ix, iy])/$
          (omialts[ilevel-1, ix, iy] - omialts[ilevel, ix, iy])
        break
      ENDIF else BEGIN
        continue
      ENDELSE
    endfor
  ENDFOR
  ENDFOR
ENDIF

o3accum  = omio3accum
return
end
