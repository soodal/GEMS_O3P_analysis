function ds_o3mmr2vmr, mmr, ppbv=ppbv, ppmv=ppmv, kg_kg=kg_kg
; ppbv
if keyword_Set(ppbv) then begin
  vmr = 28.9644 / 47.9982 * 1e9 * MMR
endif else if keyword_set(ppmv) then begin
  vmr = 28.9644 / 47.9982 * 1e6 * MMR
endif else if keyword_Set(kg_kg) then begin
  vmr = 28.9644 / 47.9982 * mmr
endif else begin
  print, 'Not choosed /ppmv or /ppbv.'
  print, 'It will calculate ppbv'
  vmr = 28.9644 / 47.9982 * 1e9 * MMR
endelse
return, vmr
end
