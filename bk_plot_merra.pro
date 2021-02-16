 file = '/data2/MERRA2/MERRA2_400.inst3_3d_asm_Np.20200616.nc4'

 o3 = h5read(file,'O3');*1e9 ; unit mass mixing ratio
 lat = h5read(file,'lat')
 lon = h5read(file,'lon')
 P = h5read(file,'lev')
 alt = h5read(file,'H')
 Temp = h5read(file,'T')

;time (UTC)
 ;03, 06, 09, 12, 15, 18, 21, 24
 IF ~keyword_set(time) then it =0 else it=time

 sz = size(o3,/dim)
 nx = sz[0] & ny = sz[1] & nz = sz[2] & nt = sz[3]
;
 TCO = fltarr(nx,ny)
 o3prfs = fltarR(nx,ny,nz-1) ;nlayer
 Pm = fltarr(nx,ny,nz)

;filter
 o3[where(o3 ge 9.99999e+14)] = !values.F_nan
 Temp[where(Temp ge 9.99999e+14)] = !values.F_nan
 alt[where(alt ge 9.99999e+14)] = !values.F_nan

 it =0
 T = reform(Temp[*,*,*,it])
 H = reform(alt[*,*,*,it])

 for ih=0,41 do Pm[*,*,ih] = P[ih]
 mmr_o3 = reform(o3[*,*,*,it])

 o3molarmass = 47.9982
;convert mass mxing ratio into vmr
 vmr_O3 = 28.9644/47.9982 * reform(o3[*,*,*,it])
 vmr_o3_ds = ds_mmr2vmr(reform(o3[*,*,*,it]), molarmass = o3molarmass)
 print, mean(vmr_o3-vmr_o3_ds, /nan)
 stop


;vmr to No3 (number density , molec/cm3)
 k = 1.3807 * 10.^(-19); J K-1 molc-1 ; 10^4 * kg * cm^2 s^-2

; pV = NkT

;O3 number density
 NO3 = vmr_o3 * Pm/(T*k)

;Air number density
 Nair = 2.69*10.^(16) ; molec/cm ^2

 du_tmp1 = fltarr(nx,ny,nz-1)
 dh      = fltarr(nx,ny,nz-1)

 FOR ih=0,nz-2 DO BEGIN
   dh = h[*,*,ih+1] - h[*,*,ih] 
   o3prfs[*,*,ih] = (no3[*,*,ih]+no3[*,*,ih+1])/2/nair*dh * 100 ;m -> cm thickness
 ENDFOR

 tco = total(o3prfs,/nan,3)

 ip = value_locate(p,300)
 sco =  total(o3prfs[*,*,0:ip],/nan,3)
trco =  total(o3prfs[*,*,ip+1:nz-2],/nan,3)

;-----------------------------------------------
 NAN = !values.F_nan
 olon=REPLICATE(NaN,Sz[0],Sz[1])
 FOR iy=0,Sz[1]-1 DO BEGIN
     olon[*,iy]=lon[*]
 ENDFOR

olat=REPLICATE(NaN,Sz[0],Sz[1])
 FOR ix=0,Sz[0]-1 DO BEGIN
     olat[ix,*]=lat[*]
 ENDFOR

end
