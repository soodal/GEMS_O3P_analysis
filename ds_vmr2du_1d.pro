;===================================================================================================
; PROGRAM NAME:
;  ds_vmr2du
;
; PURPOSE:
;  convert volum mixing ratio(VMR) to dobson unit(DU)
;
; PROCESS:
;
; REFERENCE:
;
; INPUT: 
;   vmr_o3(nx, ny, nz) [volume mixing ratio]
;   pressure(nx, ny, nz) [hPa]
;   temperature(nx, ny, nz) [K]
;   altitude(nx, ny, nz) [m]
;
; OUTPUT:
;   o3_du_layer(nx, ny, nz-1) [DU]
;
; DEVELOPER:
;   Dae Sung Choi(Choi)
;   Satellite Remote Sensing Laboratory
;   Department of Atmospheric Sciences
;   College of Natural Science
;   Pusan National University
;                  Tel: +82-51-510-2172
;                  E-mail: daesungchoi@pusan.ac.kr
;
;  Copyright (C) 2020 Satellite Remote Sensing Laboratory, PNU
;  All right reserved.
;
;===================================================================================================
pro ds_vmr2du_1d, vmr_o3, pressure, temperature, altitude, o3_du_layer, $
  incresing_height=incheight

if not keyword_set(incheight) then begin
  print, 'increasing_height keyword not set.'
  print, 'the vertical axis of altitude index 0 for the bottom n-1 for top.'
  ;print, 'if this keyword as set to 0 then decreasing order.'
  ;pritn, 'if not set this keyword, default value is 1 to increasing order.'
  incheight = 1
endif


print, 'ds_vmr2du, o3vmr,pressure, temperature, altitude'

vmr_o3_size = n_elements(vmr_o3)

nz = vmr_o3_size


o3_du_layer = fltarr(nz) ;nlayer

;vmr to N_o3 (number density , molec/cm3)
k = 1.3807 * 10.^(-19); J K-1 molc-1 ; 10^4 * kg * cm^2 s^-2

; pV = NkT

;O3 number density
N_o3 = vmr_o3 * pressure/(temperature*k)

;Air number density
N_air = 2.69*10.^(16) ; molec/cm ^3

du_tmp1 = fltarr(nz)
dh      = fltarr(nz)

FOR ih=0,nz-2 DO BEGIN
  dh = altitude[ih+1] - altitude[ih]  ; meter
  dh = dh *  100 ; centimeter
  o3_du_layer[ih] = (N_o3[ih]+N_o3[ih+1])/2/N_air*dh ;cm thickness
ENDFOR

return

;tco = total(o3_du_layer,/nan,3)

;ip = value_locate(p,300)
;sco =  total(o3_du_layer[*,*,0:ip],/nan,3)
;trco =  total(o3_du_layer[*,*,ip+1:nz-2],/nan,3)

;;-----------------------------------------------
;NAN = !values.F_nan
;olon=REPLICATE(NaN,Sz[0],Sz[1])
;FOR iy=0,Sz[1]-1 DO BEGIN
;olon[*,iy]=lon[*]
;ENDFOR

;olat=REPLICATE(NaN,Sz[0],Sz[1])
;FOR ix=0,Sz[0]-1 DO BEGIN
;olat[ix,*]=lat[*]
;ENDFOR

END
