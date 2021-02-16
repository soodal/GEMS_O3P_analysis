;merra2o3fn = '/data2/MERRA2/MERRA2_400.inst3_3d_chm_Nv.20200616.nc4'
;AIRDENS = ds_read_nc(merra2o3fn, 'AIRDENS')
;CO = ds_read_nc(merra2o3fn, 'CO')
;delp = ds_read_nc(merra2o3fn, 'DELP') ; [Pa]

;lat = ds_read_nc(merra2o3fn, 'lat')
;lev = ds_read_nc(merra2o3fn, 'lev')
;lon = ds_read_nc(merra2o3fn, 'lon')
;O3 = ds_read_nc(merra2o3fn, 'O3')
;PS = ds_read_nc(merra2o3fn, 'PS')
;time = ds_read_nc(merra2o3fn, 'time')

merra2asmfn = '/data2/MERRA2/MERRA2_400.inst3_3d_asm_Np.20200616.nc4'
EPV = ds_read_nc(merra2asmfn, 'EPV')
H = ds_read_nc(merra2asmfn, 'H')
lat = ds_read_nc(merra2asmfn, 'lat')
lev = ds_read_nc(merra2asmfn, 'lev')
lon = ds_read_nc(merra2asmfn, 'lon')
O3 = ds_read_nc(merra2asmfn, 'O3')
OMEGA = ds_read_nc(merra2asmfn, 'OMEGA')
PHIS = ds_read_nc(merra2asmfn, 'PHIS')
PS = ds_read_nc(merra2asmfn, 'PS')
QI = ds_read_nc(merra2asmfn, 'QI')
QL = ds_read_nc(merra2asmfn, 'QL')
QV = ds_read_nc(merra2asmfn, 'QV')
RH = ds_read_nc(merra2asmfn, 'RH')
SLP = ds_read_nc(merra2asmfn, 'SLP')
T = ds_read_nc(merra2asmfn, 'T')
time = ds_read_nc(merra2asmfn, 'time')
U = ds_read_nc(merra2asmfn, 'U')
V = ds_read_nc(merra2asmfn, 'V')


sz = size(o3)
dim1 = sz[1]
dim2 = sz[2]
fl = 20;n_elements(lev)

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

plot_sat_proj, trco, lon2d, lat2d, $
  title='MERRA2 O3 Column SfcTo300hPa', $
  cb_title='', $
  range=[20, 60], $
  pngfile='./plot/merra_o3column_under300hpa0616.png', /scp_send

plot_sat_proj, tco, lon2d, lat2d, $
  title='MERRA2 O3 Total Column', $
  cb_title='MERRA2 O3 Column SfcTo300hPa', $
  range=[230, 340], $
  pngfile='./plot/merra_o3total_0616.png', /scp_send

plot_sat_proj, sco, lon2d, lat2d, $
  title='MERRA2 O3 Column 300hPaToTOA', $
  range=[190, 320], $
  pngfile='./plot/merra_o3strato_0616.png', /scp_send

end
