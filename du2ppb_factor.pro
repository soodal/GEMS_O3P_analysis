; ===================================================================
; Ozone (ppbv) = convfac * ozone (DU)  / delta P (in atm)
; Note: alt (km); pres (mb)
; ===================================================================
; conversion factor at 1 atmosphere: 
; 0.00787868 DU*ppmv^(-1)*Pa^(-1) = 0.798307 DU *  ppbv^(-1) * atm^-1
; 1 / conversion factor = 1.25265
; ===================================================================

pro du2ppb_factor, alt=alt, pres=pres, convfac

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

