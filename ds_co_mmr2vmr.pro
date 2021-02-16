function ds_co_mmr2vmr, mmr, ppbv=ppbv, ppmv=ppmv
; ppbv
if keyword_Set(ppbv) then begin
  vmr = 28.9644 / 47.9982 * 1e9 * MMR
endif else if keyword_set(ppmv) then begin
  vmr = 28.9644 / 47.9982 * 1e6 * MMR
endif else begin
  print, 'Not choosed /ppmv or /ppbv.'
  print, 'It will calculate ppbv'
  vmr = 28.9644 / 47.9982 * 1e9 * MMR
endelse
return, vmr
end
