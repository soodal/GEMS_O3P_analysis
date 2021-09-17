

; Integrate ozone column from surface to burst altitude
pro get_intoz, nl, p, o3, intoz

; o3: mpa
; p : hPa, bottom-up
IF P(1) gt P(0) then begin
  print, 'check beottom-up in pressure'
  stop
ENDIF
IF p(1) eq P(2) and p(2) eq p(3) then begin
  intoz =  !values.f_nan
  return
endif


du2ppb_factor, pres=p, convfac
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
