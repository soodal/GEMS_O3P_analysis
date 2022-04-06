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

m2mm = 1000.
mm2du = 100.

;merra2asmfn = '/data2/MERRA2/MERRA2_400.inst3_3d_asm_Np.20200616.nc4'
merra2asmfn = '/data/private/soodal/MERRA2_collocated_on_gems/merra2_tavg3_3d_asm_collocated_on_gemsl1c_rad_20210329_0345.nc4'
;EPV = ds_read_nc(merra2asmfn, 'EPV')
H = ds_read_nc(merra2asmfn, 'H')
lat = ds_read_nc(merra2asmfn, 'lat')
lev = ds_read_nc(merra2asmfn, 'lev')
lon = ds_read_nc(merra2asmfn, 'lon')
O3 = ds_read_nc(merra2asmfn, 'O3')
;OMEGA = ds_read_nc(merra2asmfn, 'OMEGA')
;PHIS = ds_read_nc(merra2asmfn, 'PHIS')
PS = ds_read_nc(merra2asmfn, 'PS')
PL = ds_read_nc(merra2asmfn, 'PL')
;QI = ds_read_nc(merra2asmfn, 'QI')
;QL = ds_read_nc(merra2asmfn, 'QL')
;QV = ds_read_nc(merra2asmfn, 'QV')
;RH = ds_read_nc(merra2asmfn, 'RH')
;SLP = ds_read_nc(merra2asmfn, 'SLP')
T = ds_read_nc(merra2asmfn, 'T')
;time = ds_read_nc(merra2asmfn, 'time')
;U = ds_read_nc(merra2asmfn, 'U')
;V = ds_read_nc(merra2asmfn, 'V')


sz = size(o3)
dim1 = sz[1]
dim2 = sz[2]
;fl = 20;n_elements(lev)

itime = 1

pres = 0.0
o3under300 = fltarr(sz[1], sz[2])
o3tmp = reform(o3[*, *, *])
o3total = total(o3tmp, 3)

nl = 72

pressure_mid = fltarr(sz[1], sz[2], nl)
for ih=0,nl-1 do BEGIN
  pressure_mid[*,*,ih] = pl[*, *, ih] ; pl: [Pa]
endfor

mmr_o3 = reform(o3[*,*,*])

;convert mass mxing ratio into vmr
vmr_o3 = ds_o3mmr2vmr(o3,/kg_kg)

;vmr to Number densityNo3 (number density , molec/cm3)

; boltzmann constant k per cm^2
k = 1.3807 * 10.^(-23);[J K^-1] ; ; 10^4 * kg * cm^2 s^-2

;number of O3 molecular
number_o3 = vmr_o3 * pressure_mid/(T*k)

;Air number density
Nair = 2.69*10.^(25) ; molec/m ^3

;du_tmp1 = fltarr(sz[1],sz[2],nl)
dh      = fltarr(sz[1],sz[2],nl)

o3prfs = fltarr(sz[1], sz[2], nl)

h1 = fltarr(sz[1], sz[2], nl+1)

; calculate the grid boundary height
for ih=0, nl-2 do BEGIN
  h1[*, *, ih] = h[*, *, ih] + (h[*, *, ih] - h[*, *, ih+1])/2
endfor
h1[*, *, nl-1] = h[*, *, nl-1] + (h[*, *, nl-2] - h[*, *, nl-1])/2
h1[*, *, nl] = h[*, *, nl-1] - (h[*, *, nl-2] - h[*, *, nl-1])/2

FOR ih=0,nl-1 DO BEGIN
  dh = h1[*, *, ih] - h1[*, *, ih+1]
  o3prfs[*,*,ih] =  number_o3[*,*,ih]/nair*dh * m2mm * mm2du
ENDFOR

o3prfs[where(o3prfs GE 1000 or o3prfs lt 0)] = !values.f_nan
tco = total(o3prfs[*, *, 0:nl-1],/nan,3)

trco = fltarr(sz[1], sz[2])
sco = fltarr(sz[1], sz[2])
for ix=0, sz[1]-1 do begin
  print, ix
  for iy=0, sz[2]-1 do begin
    if total(finite(pl[ix, iy, *])) ne 0 then begin
      ip = value_locate(reform(PL[ix, iy, *]) ,30000)
      trco[ix, iy] = total(o3prfs[ix,iy,ip+1:nl-18],/nan,3)

      sco[ix, iy] = total(o3prfs[ix, iy,0:ip],/nan,3)
      ;print, ix, iy, trco[ix, iy,
    endif
  endfor
endfor

;ds_xymake2d, lon, lat, lon2d, lat2d

plot_sat_proj, trco, lon, lat, $
  title='MERRA2 O3 Sfc-300hPa', $
  cb_title='[DU]', $
  range=[0, 60], $
  pngfile='./plot/merra2_collocated/20210329/merra_o3column_under300hpa20210329.png';, /scp_send

plot_sat_proj, tco, lon, lat, $
  title='MERRA2 O3 Total Column', $
  cb_title='[DU]', $
  range=[220, 430], $
  pngfile='./plot/merra2_collocated/20210329/merra_o3total_20210329.png';, /scp_send
;showme, tco
stop

plot_sat_proj, reform(o3prfs[*, *, nl-1]), lon, lat, $
  title='MERRA2 O3 Layer', $
  cb_title='MERRA2 O3 [DU]', $
  range=[min(o3prfs[*, *, nl-1], /nan), max(o3prfs[*, *, nl-1], /nan)], $
  pngfile='./plot/merra2_collocated/20210621/merra_o3_layer_20210621.png';, /scp_send

;plot_sat_proj, sco, lon2d, lat2d, $
  ;title='MERRA2 O3 Column 300hPaToTOA', $
  ;range=[190, 320], $
  ;pngfile='./plot/merra2_collocated/20210621/merra_o3strato_20210621.png';, /scp_send

gemso3pfn = '/home/soodal/data/FNL/2x_ozmin_ozmax_fnl_daily_310_340/GK2_GEMS_L2_O3P_20210621_0345_ML_prec000000_maxit10_ozminmax2x_BIN4x4.nc'
gemsl2o3p = ds_read_gems_l2_o3p(gemso3pfn)

ds_gems_l2o3p_accum, gemsl2o3p, gemso3accum, hpa=300

plot_sat_proj, gemso3accum - trco, lon, lat, $
  title='GEMS O3P Under 300 hPa - MERRA2 O3 Under 300 hPa', $
  cb_title='O3[DU]', $
  colortable=72, $
  /ctreverse, $
  range=[-20, 20], $
  pngfile='./plot/merra2_collocated/20210621/gemso3p_310340_FNL_ML-merra_o3column_under300hpa20210621.png';, /scp_send
end
