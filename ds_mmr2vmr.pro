function ds_mmr2vmr, mmr, ppbv=ppbv, ppmv=ppmv
; ppbv
if keyword_Set(ppbv) then begin
  vmr = 28.9644 / 47.9982 * 1e9 * MMR
endif else if keyword_set(ppmv) then begin
  vmr = 28.9644 / 47.9982 * 1e6 * MMR
endif else begin
  print, 'please choose /ppbv or /ppmv'
  return, !values.f_nan
endelse
return, vmr
end
