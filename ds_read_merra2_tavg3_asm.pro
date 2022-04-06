function ds_read_merra2_tavg3_asm, fn, itime=itime
if not keyword_set(itime) then BEGIN
  print, 'keyword itime is not set. itime is set to 1'
  itime = 1
endif
;fn = '/data/MODEL/MERRA2/2020/10/MERRA2_400.tavg3_3d_asm_Nv.' + yyyy + mm + dd + '.nc4'
;fn = '/data/MODEL/MERRA2/2020/10/MERRA2_400.tavg3_3d_asm_Nv.' + yyyy + mm + dd + '.nc4'
delp = ds_read_nc(fn, 'DELP')
EPV = ds_read_nc(fn, 'EPV')
H = ds_read_nc(fn, 'H')
lat = ds_read_nc(fn, 'lat')
lev = ds_read_nc(fn, 'lev')
lon = ds_read_nc(fn, 'lon')
O3 = ds_read_nc(fn, 'O3')
OMEGA = ds_read_nc(fn, 'OMEGA')
PHIS = ds_read_nc(fn, 'PHIS')
PS = ds_read_nc(fn, 'PS')
PL = ds_read_nc(fn, 'PL')
QI = ds_read_nc(fn, 'QI')
QL = ds_read_nc(fn, 'QL')
QV = ds_read_nc(fn, 'QV')
RH = ds_read_nc(fn, 'RH')
SLP = ds_read_nc(fn, 'SLP')
T = ds_read_nc(fn, 'T')
time = ds_read_nc(fn, 'time')
U = ds_read_nc(fn, 'U')
V = ds_read_nc(fn, 'V')

sz = size(o3)
dim1 = sz[1]
dim2 = sz[2]
dim3 = sz[3]

pres = 0.0
o3under300 = fltarr(sz[1], sz[2])
o3tmp = reform(o3[*, *, *, 1])
o3total = total(o3tmp, 3)

mmr_o3 = reform(o3[*,*,*,itime])

;convert mass mxing ratio into vmr
vmr_O3 = 28.9644/47.99 * reform(o3[*,*,*,itime])

vmr_o3_ds = ds_mmr2vmr(o3[*, *, *, itime], /ppmv)

;vmr to No3 (number density , molec/cm3)
k = 1.3807 * 10.^(-19); J K-1 molc-1 ; 10^4 * kg * cm^2 s^-2

;O3 number density
NO3 = vmr_o3 * h/(T*k)

;Air number density
Nair = 2.69*10.^(16) ; molec/cm ^2

o3du = fltarr(sz[1], sz[2], sz[3])

dh      = fltarr(sz[1],sz[2],sz[3])
dh[*, *, 0:70] = h[*, *,0:70] - h[*, *, 1:71] 
dh[*, *, 71] = h[*, *, 70] - h[*, *, 71]

o3du = double(no3/nair*dh * 100) ;m -> cm thickness

;o3du[where(o3du GE 1000 or o3du lt 0)] = !values.f_nan

ds_xymake2d, lon, lat, lon2d, lat2d
result = {delp:delp, EPV:EPV, H:H, lon:lon, lat:lat, lev:lev, O3:O3, OMEGA:OMEGA, PHIS:PHIS, $
  PS:PS, PL:PL, QI:QI, QL:QL, QV:QV, RH:RH, SLP:SLP, T:T, time:time, U:U, V:V, $
  lat2d:lat2d, lon2d:lon2d, o3du:o3du}
return, result
end
