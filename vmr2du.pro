function vmr3du, vmr, 
nl = 42

for ih=0,nl-1 do BEGIN
  Pm[*,*,ih] = P[ih]
endfor

mmr_o3 = reform(o3[*,*,*,it])

;convert mass mxing ratio into vmr
vmr_O3 = 28.9644/47.99 * reform(o3[*,*,*,it])

;vmr to No3 (number density , molec/cm3)
k = 1.3807 * 10.^(-19); J K-1 molc-1 ; 10^4 * kg * cm^2 s^-2

;O3 number density
NO3 = vmr_o3 * Pm/(T*k)

;Air number density
Nair = 2.69*10.^(16) ; molec/cm ^2

du_tmp1 = fltarr(nx,ny,nl-1)
dh      = fltarr(nx,ny,nl-1)

FOR ih=0,nl-2 DO BEGIN
dh = ( h[*,*,ih+1]-h[*,*,ih] )
o3prfs[*,*,ih] =  (no3[*,*,ih]+no3[*,*,ih+1])/2/nair*dh * 100
;m -> cm thickness
ENDFOR
tco = total(o3prfs,/nan,3)


end
