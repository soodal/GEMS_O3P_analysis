function ds_mmr2vmr, mmr, molarmass=molarmass, ppbv=ppbv, ppmv=ppmv

if not keyword_Set(molarmass) then begin
  print, 'keyword molarmass[g/mol] not set. Using O3 molarmass 47.9982g/mol.'
  molarmass = 47.9982
endif

if keyword_Set(ppbv) then begin
  vmr = 28.9644 / molarmass * 1e9 * MMR
endif else if keyword_set(ppmv) then begin
  vmr = 28.9644 / molarmass * 1e6 * MMR
endif else begin
  vmr = 28.9644 / molarmass * MMR
endelse
return, vmr
end
