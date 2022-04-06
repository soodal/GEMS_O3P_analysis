function search_closest_pixel, arrayx, arrayy, pixx, pixy, maxlimit=maxlimit
if not keyword_set(maxlimit) then begin
  maxlimit = 0.2
endif
;sz = size(arrayx, /dimension)

;xnumsmall = round(sz[0]/10.)
;ynumsmall = round(sz[1]/10.)

;apparr = fltarr(round(sz[0]/10.), round(sz[1]/10.))
;apparr[*] = !values.f_nan

;for iy=0, ynumsmall-1 do begin
  ;for ix=0, xnumsmall-1 do begin
    ;apparr[ix, iy] = sqrt((arrayx[ix*10, iy*10]-x)^2 - (arrayy[ix*10, iy*10]-pixy)^2)
  ;endfor
;endfor

;minval = min(apparr, minidx)
;xy = array_indices(apparr, minidx)

;newarr = fltarr(sz[0], sz[1])
;newarr[*] = !values.f_nan
;for iy=round((xy[1]-1)*10.), round((xy[1]+1)*10.) do begin
  ;for ix=round((xy[0]-1)*10.), round((xy[0]+1)*10.) do begin
    ;if ix ge 0 and ix le sz[0]-1 and iy ge 0 and iy le sz[1]-1 then begin
      ;newarr[ix, iy] = sqrt((arrayx[ix, iy]-x)^2 + (arrayy[ix, iy]-pixy)^2)
    ;endif
  ;endfor
;endfor

;minval = min(newarr, newidx, /nan)
;newxyidx = array_indices(newarr, newidx)
;stop
if n_elements(pixx) ne 1 and n_elements(pixy) ne 1 and n_elements(pixx) eq n_elements(pixy) then begin
  results = []
  for i = 0, n_elements(pixx)-1 do begin
    arr = sqrt((arrayx-pixx[i])^2 + (arrayy-pixy[i])^2)
    minval = min(arr, minidx, /nan)
    ;newxyidx = array_indices(arrayx, minidx)
    if minval le maxlimit then begin
      result = minidx
    endif else begin
      result = -999
    endelse
    results = [results, result]
  endfor
  result = results
endif else begin
  arr = sqrt((arrayx-pixx)^2 + (arrayy-pixy)^2)
  minval = min(arr, minidx, /nan)
  ;newxyidx = array_indices(arrayx, minidx)
  if minval le maxlimit then begin
    result = minidx
  endif else begin
    result = -999
  endelse
endelse

return, result

end
