; ===================================================================
; Ozone (ppbv) = convfac * ozone (DU)  / delta P (in atm)
; Note: alt (km); pres (mb)
; ===================================================================
; conversion factor at 1 atmosphere: 
; 0.00787868 DU*ppmv^(-1)*Pa^(-1) = 0.798307 DU *  ppbv^(-1) * atm^-1
; 1 / conversion factor = 1.25265
; ===================================================================

pro get_du2ppb_factor, convfac, alt=alt, pres=pres

r0 = 6367.45   ; earth radius
c0 = 1.25265

if arg_present(alt) ge 0 and n_elements(alt) gt 0 then begin
    convfac = c0 * (r0/ (r0 + alts))^2
endif else if arg_present(pres) ge 0 and n_elements(pres) gt 0 then begin
    alts = -alog10(pres / 1000.) * 16.  ; convert to pressure altitude
    convfac = c0 * (r0/ (r0 + alts))^2
endif else begin
    convfac = c0   ; work pretty well in the troposphere
endelse

return
end


; Integrate ozone column from surface to burst altitude
;pro get_intoz, nl, p, o3, intoz
pro convert_mpa2du, nl, p, o3, intoz, cols   ; geun added
; o3: mpa
; p : hPa, bottom-up

IF P(1) gt P(0) then begin
  print, 'Its not bottom-up in pressure. It should be bottom-up'
  stop
ENDIF

IF p(1) eq P(2) and p(2) eq p(3) then begin
  intoz =  !values.f_nan
  return
endif

get_du2ppb_factor, convfac ,pres=p

convfac = (convfac(0:nl-2) + convfac(1:nl-1))/2.0
tmp = (p(0:nl-2) - p(1:nl-1)) / 1013.25 / convfac

cols=fltarr(nl) 
for i = 1, nl-1 do begin
    cols(i) = cols(i-1) + (o3(i-1)/p(i-1) + o3(i)/p(i)) * 5.0E3 * tmp(i-1)
endfor
intoz = cols(nl-1)
cols = cols(1:nl-1) - cols(0:nl-2)
return
end
