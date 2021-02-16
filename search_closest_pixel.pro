function search_closest_pixel, inputx, inputy, x, y
;sz = size(inputx, /dimension)

;xnumsmall = round(sz[0]/10.)
;ynumsmall = round(sz[1]/10.)

;apparr = fltarr(round(sz[0]/10.), round(sz[1]/10.))
;apparr[*] = !values.f_nan

;for iy=0, ynumsmall-1 do begin
  ;for ix=0, xnumsmall-1 do begin
    ;apparr[ix, iy] = sqrt((inputx[ix*10, iy*10]-x)^2 - (inputy[ix*10, iy*10]-y)^2)
  ;endfor
;endfor

;minval = min(apparr, minidx)
;xy = array_indices(apparr, minidx)

;newarr = fltarr(sz[0], sz[1])
;newarr[*] = !values.f_nan
;for iy=round((xy[1]-1)*10.), round((xy[1]+1)*10.) do begin
  ;for ix=round((xy[0]-1)*10.), round((xy[0]+1)*10.) do begin
    ;if ix ge 0 and ix le sz[0]-1 and iy ge 0 and iy le sz[1]-1 then begin
      ;newarr[ix, iy] = sqrt((inputx[ix, iy]-x)^2 + (inputy[ix, iy]-y)^2)
    ;endif
  ;endfor
;endfor

;minval = min(newarr, newidx, /nan)
;newxyidx = array_indices(newarr, newidx)
;stop

arr = sqrt((inputx-x)^2 + (inputy-y)^2)
minval = min(arr, minidx, /nan)
;newxyidx = array_indices(inputx, minidx)
if minval le 0.1 then begin
  result = minidx
endif else begin
  result = -999
endelse


return, result

end
