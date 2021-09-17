pro ds_merra2_get_o3underp, data, underp, itime=itime

if keyword_set(itime) then BEGIN
  print, 'keyword itime is not set. itime is set to 0'
  itime = 0
endif

data = {EPV:EPV, $
  H:H, $ 
  lat:lat, $
  lev:lev, $
  lon:lon, $
  O3:O3, $
  OMEGA:OMEGA, $
  PHIS:PHIS, $
  PS:PS, $
  QI:QI, $
  QL:QL, $
  QV:QV, $
  RH:RH, $
  SLP:SLP, $
  T:T, $
  time:time, $
  U:U, $
  V:V}

EPV = data.EPV
H = data.H
lat = data.lat
lev = data.lev
lon = data.lon
O3 = data.o3
OMEGA = data.OMEGA
PHIS = data.PHIS
PS = data.PS
QI = data.QI
QL = data.QL
RH = data.RH
SLP = data.SLP
T = data.T
time = data.time
U = data.U
v = data.V

sz = size(o3)
dim1 = sz[1]
dim2 = sz[2]

itime = 1

pres = 0.0
o3under300 = fltarr(sz[1], sz[2])
o3tmp = reform(o3[*, *, *, 1])
o3total = total(o3tmp, 3)

nl = 42

levm = fltarr(sz[1], sz[2], nl)
for ih=0,nl-1 do BEGIN
  levm[*,*,ih] = lev[ih]
endfor

mmr_o3 = reform(o3[*,*,*,itime])

;convert mass mxing ratio into vmr
vmr_O3 = 28.9644/47.99 * reform(o3[*,*,*,itime])

vmr_o3_ds = ds_mmr2vmr(o3[*, *, *, itime], /ppmv)

;vmr to No3 (number density , molec/cm3)
k = 1.3807 * 10.^(-19); J K-1 molc-1 ; 10^4 * kg * cm^2 s^-2

;O3 number density
NO3 = vmr_o3 * levm/(T*k)

;Air number density
Nair = 2.69*10.^(16) ; molec/cm ^2

du_tmp1 = fltarr(sz[1],sz[2],nl-1)
dh      = fltarr(sz[1],sz[2],nl-1)

o3prfs = fltarr(sz[1], sz[2], nl-1)
FOR ih=0,nl-2 DO BEGIN
  dh = ( h[*,*,ih+1]-h[*,*,ih] )
  o3prfs[*,*,ih] =  (no3[*,*,ih]+no3[*,*,ih+1])/2/nair*dh * 100
  ;m -> cm thickness
ENDFOR

a = o3prfs

o3prfs[where(o3prfs GE 1000 or o3prfs lt 0)] = !values.f_nan
tco = total(o3prfs[*, *, 0:nl-2],/nan,3)

ip = value_locate(lev,300)
trco =  total(o3prfs[*,*,0:ip],/nan,3)

sco =  total(o3prfs[*,*,ip+1:nl-2],/nan,3)

ds_xymake2d, lon, lat, lon2d, lat2d

underp = trco
return
end
