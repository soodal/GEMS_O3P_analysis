;===================================================================================================
; PROGRAM NAME:
;  vmr2du_geun
;
; PURPOSE:
;  convert volum mixing ratio(VMR) to dobson unit(DU)
;
; PROCESS:
;
; REFERENCE:
;
; INPUT: (Top2Bot)
;  O3(nlevel, DU), press(nlevel, hPa), Temp(nlevel, K), Alt(nlevel, km) 
;
; OUTPUT:
;  O3_vmr(nlevel-1/nlevel, mol/mol) 
;
; DEVELOPER:
;  Kanghyun Baek (BK)
;  Daegeun Shin (geun)
;  Satellite Remote Sensing Laboratory
;  Division of Earth Environmental System
;  College of Natural Science
;  Pusan National University
;                  Tel: +82-51-510-2172
;                  E-mail: iamdaegeun@gmail.com
;
;  Copyright (C) 2019 Satellite Remote Sensing Laboratory, PNU
;  All right reserved.
;
;===================================================================================================

;+---------------------------------------------------------------------------+
; Function Name : cum_total
;+---------------------------------------------------------------------------+
function cum_total,y

    if (n_elements(y) eq 0) then return,0.
    if (n_elements(y) eq 1) then return,y[0]

    result = dblarr(n_elements(y))
    result[0] = y[0]

    for i=1,n_elements(y)-1 do result[i] = result[i-1] + y[i]
    return,result
end



FUNCTION vmr2du_geun, O3_vmr, Press, Temp, alt=alt, interpol=interp

  nlevel = n_elements(press)
  IF ~keyword_set(alt) THEN BEGIN
    G = 9.80665
    R = 287.053
    GMR = G/R
    T = reverse(temp) & P= reverse(press)
    H = fltarr(nlevel)
    FOR i=1,nlevel-1 DO $
        H[i] = 1./1000.*(1/GMR * (T[i]+T[i-1])/2 * alog(P[i-1]/P[i]) ) + H[i-1]
    alt = reverse(H)
  ENDIF

  nl=N_ELEMENTS(press)
  o3p_du=FLTARR(nl-1)
  zmid=FLTARR(nl-1)
  FOR ilay=0, nl-2 DO BEGIN
    rho1=press[ilay]  /temp[ilay]
    rho2=press[ilay+1]/temp[ilay+1]
    hcon=0.5*1.0D+5*2.68675D+19*273.15D/1013.15D*(alt[ilay]-alt[ilay+1])
                                                                              ;; 0.5*10^5[hPa->mPa]*Na*T0/P0*dz
    zmid[ilay]=(alt[ilay]+alt[ilay+1])/2.
    gascolumns1=hcon*rho1*o3_vmr[ilay]
    gascolumns2=hcon*rho2*o3_vmr[ilay+1]
    xx=(gascolumns1+gascolumns2)/(2.68675D+16)  ; ozon[ppmv] -> ozone[du]
    o3p_du[ilay]=xx
  ENDFOR

  IF keyword_set(interp) THEN BEGIN ; geun modified output option
    cumo3 = cum_total(o3p_du)
    tmpcum =  interpol(cumo3,zmid,alt)
    o3p_du=[tmpcum[0], tmpcum(1:nl-1) - tmpcum(0:nl-2)]
    IF o3p_du[0] LE 0 THEN o3p_du[0]=o3p_du[1]/10.  
  ENDIF ; interp
  RETURN, o3p_du

END
