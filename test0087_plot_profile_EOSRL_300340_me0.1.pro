

pos = [0.12, 0.1, 0.9, 0.85]

write_flag = 0
;=============================================================================
; start
;=============================================================================
jd_list = timegen(start=julday(3, 29, 2021, 3, 45), $
  final=julday(3, 29, 2021, 3, 45), units='Hours')

;jd_list = timegen(start=julday(6, 21, 2021, 3, 45), $
  ;final=julday(6, 21, 2021, 3, 45), units='Hours')

itime = 0
caldat, jd_list[itime], month, day, year, hour, minute
yyyy = string(year, format='(i04)')
mm = string(month, format='(i02)')
dd = string(day, format='(i02)')
hh = string(hour, format='(i02)')
mi = string(minute, format='(i02)')

caldat, jd_list[itime]-1, month_y, day_y, year_y, hour_y, minute_y
yyyy_y = string(year_y, format='(i04)')
mm_y = string(month_y, format='(i02)')
dd_y = string(day_y, format='(i02)')
hh_y = string(hour_y, format='(i02)')
mi_y = string(minute_y, format='(i02)')

datetime_str = yyyy + mm + dd + '_' + hh + mi
datetime_y_str = yyyy_y + mm_y + dd_y
;=============================================================================
; GEMS O3P NIER
;=============================================================================
search_str = '/home/soodal/data/ln/GEMS/L2O3P/GK2_GEMS_L2_O3P_' + $
  datetime_str + '_4x4.nc'

fl = file_search(search_str)
nier_data = ds_read_gems_l2_o3p(fl[-1])

nier_lat = nier_data.latitude
nier_lon = nier_data.longitude

nier_o3p = nier_data.o3
nier_p = nier_data.pressure

;=============================================================================
; GEMS O3P EOSRL
;=============================================================================
search_str = '/data/private/soodal/EOSRL/' + $
  'GK2_GEMS_O3P_' + datetime_str + '_*300340_me0.1.nc.debug.nc'


fl = file_search(search_str)
eosrl_data = ds_read_gems_l2_o3p(fl[-1])

;lat = eosrl_data.latitude
;lon = eosrl_data.longitude

eosrl_o3p = eosrl_data.o3
eosrl_p = eosrl_data.pressure

;=============================================================================
; GEMS O3P Synthetic
;=============================================================================
search_str = '/data/private/soodal/softcal_test/corrected/model_300_340/GK2_GEMS_L2_O3P_20210329_0345_winliminit300_prec000000_climML_b4x4_o20211221T225814Z_clima_maxiter10_ecf0.nc'

fl = file_search(search_str)
synthetic_o3p_data = ds_read_gems_l2_o3p(fl[-1])


synthetic_o3p = synthetic_o3p_data.o3
synthetic_o3p_p = synthetic_o3p_data.pressure

;=============================================================================
; SYNTHETIC FITSPEC
;=============================================================================
search_str = '/home/soodal/data/merra2_residual/residuals/' + $
  'GK2_GEMS_L2_O3P_' + datetime_str + '_wli300_prec000000_climML_b4x4_*.nc'

fl = file_search(search_str)
synthetic = ds_read_gems_l2_o3p(fl[-1])

corrected_rad_all = synthetic.corrected_rad_all
corrected_irrad_all = synthetic.corrected_irrad_all
div_sun_all = synthetic.div_sun_all
div_rad_all = synthetic.div_rad_all

; synthetic is the output of O3P.
; It has 300~340 wavelength range
synthetic_simrad = synthetic.simulatedradiances[*, *, 0:203]

syntheticwav = synthetic.wavelengths[*, 0:203]


simrad_sz = size(synthetic_simrad, /dimension)
synthetic_fitspec = dblarr(simrad_sz[0], simrad_sz[1], simrad_sz[2])

for iw = 0, 203 do begin
  synthetic_fitspec[*, *, iw] = corrected_rad_all[*, *, iw] / corrected_irrad_all[*, *, iw] * div_rad_all / div_sun_all
endfor


; solar zenith angle
sza = dblarr(simrad_sz[0], simrad_sz[1])
for iw = 0, 203 do begin
    sza[*, *] = synthetic.solarzenithangle
endfor

sza = synthetic.solarzenithangle
vza = synthetic.viewingzenithangle

lon = synthetic.longitude
lat = synthetic.latitude

nanidx = where(lon lt -1e29, /null)
lon[nanidx] = !values.f_nan

nanidx = where(lat lt -1e29, /null)
lat[nanidx] = !values.f_nan

;=============================================================================
; GEMS Observation
;=============================================================================
;l1cradfn = '/data2/L1C_GEMS/L1C/4x4/GK2_GEMS_L1C_20210621_0345_NOR_694_4x4.nc'
;l1cradfn = '/data/private/soodal/ln/GEMS/L1C/GK2_GEMS_L1C_20210621_0345_4x4.nc'
nier_l1cradfn_search = '/data/private/soodal/ln/GEMS/L1C/GK2_GEMS_L1C_' + datetime_str + '_4x4.nc'
nier_l1cradfn_list = file_search(nier_l1cradfn_search)
nier_l1cradfn = nier_l1cradfn_list[0]

eosrl_l1cradfn_search = '/data2/L1C_GEMS/L1C/eosrl/nier/L1C/202103/29/GK2_GEMS_L1C_20210329_0345_CAPO_NOR_694_BIN4x4.nc'
eosrl_l1cradfn_list = file_search(eosrl_l1cradfn_search)
eosrl_l1cradfn = eosrl_l1cradfn_list[0]

;irrfn = '/data2/L1C_GEMS/L1C/4x4/GK2_GEMS_IRR_20210620_4x4.nc'
irrfn = '/data/private/soodal/ln/GEMS/IRR/GK2_GEMS_IRR_20210620_4x4.nc'
nier_irrfn = '/data/private/soodal/ln/GEMS/IRR/GK2_GEMS_IRR_' + datetime_y_str + '_4x4.nc'
eosrl_irrfn = '/data2/L1C_GEMS/L1C/eosrl/nier/IRR/202103/28/GK2_GEMS_IRR_20210328_BIN4x4.nc'

;cldfn = '/data/nier_ftp/CLOUD/V03/202106/21/GK2_GEMS_L2_20210621_0345_CLOUD_FW_DPRO_BIN4x4.nc'
;cldfn = '/data/private/soodal/ln/GEMS/L2CLD/GK2_GEMS_L2_CLD_20210621_0345_4x4.nc'
cldfn = '/data/private/soodal/ln/GEMS/L2CLD/GK2_GEMS_L2_CLD_' + datetime_str + '_4x4.nc'

nier_gemsl1crad = ds_read_gems_l1c_4x4(nier_l1cradfn)
nier_gemsirr = ds_read_gems_l1c_4x4(nier_irrfn)
nier_gemscld = ds_read_gems_l2_o3p(cldfn)

eosrl_gemsl1crad = ds_read_gems_l1c_4x4(eosrl_l1cradfn)
eosrl_gemsirr = ds_read_gems_l1c_4x4(eosrl_irrfn)


nier_startidxarr = lonarr(simrad_sz[0], simrad_sz[1])

for ix = 0, simrad_sz[0]-1 do begin
    for iy = 0, simrad_sz[1]-1 do begin
        if syntheticwav[iy, 0] lt -1e29 then begin
            nier_startidxarr[ix, iy] = 0
            continue
        endif

        idx = where(nier_gemsl1crad.wavelength[iy, *] eq syntheticwav[iy, 0], /null)
        if n_elements(idx) eq 0 then begin
            stop
        endif
        nier_startidxarr[ix, iy] = idx
    endfor
endfor

nier_l1cradwav = dblarr(simrad_sz[0], simrad_sz[1], simrad_sz[2])
for ix = 0, simrad_sz[0]-1 do begin
    for iy = 0, simrad_sz[1]-1 do begin
        nier_l1cradwav[ix, iy, *] = nier_gemsl1crad.wavelength[iy, nier_startidxarr[ix, iy]:nier_startidxarr[ix, iy]+203]
    endfor
endfor

eosrl_startidxarr = lonarr(simrad_sz[0], simrad_sz[1])

for ix = 0, simrad_sz[0]-1 do begin
    for iy = 0, simrad_sz[1]-1 do begin
        if syntheticwav[iy, 0] lt -1e29 then begin
            eosrl_startidxarr[ix, iy] = 0
            continue
        endif

        idx = where(nier_gemsl1crad.wavelength[iy, *] eq syntheticwav[iy, 0], /null)
        if n_elements(idx) eq 0 then begin
            stop
        endif
        eosrl_startidxarr[ix, iy] = idx
    endfor
endfor
        
eosrl_l1cradwav = dblarr(simrad_sz[0], simrad_sz[1], simrad_sz[2])
for ix = 0, simrad_sz[0]-1 do begin
    for iy = 0, simrad_sz[1]-1 do begin
        eosrl_l1cradwav[ix, iy, *] = eosrl_gemsl1crad.wavelength[iy, eosrl_startidxarr[ix, iy]:eosrl_startidxarr[ix, iy]+203]
    endfor
endfor

;----------------------------------------------------------------------------
nier_rad = dblarr(simrad_sz[0], simrad_sz[1], simrad_sz[2])
for ix = 0, simrad_sz[0]-1 do begin
    for iy = 0, simrad_sz[1]-1 do begin
        nier_rad[ix, iy, 0:simrad_sz[2]-1] = nier_gemsl1crad.image_pixel_values[ix, iy, eosrl_startidxarr[ix, iy]:eosrl_startidxarr[ix, iy]+simrad_sz[2]-1]
    endfor
endfor

syntheticwav[where(syntheticwav le -1e29, /null)] = !values.f_nan

temp = fltarr(simrad_sz[0], simrad_sz[1])
for ix = 0, simrad_sz[0]-1 do begin
  temp[ix, *] = total(finite(syntheticwav), 2)
endfor

;gems_rad_norm = total(nier_rad, 3)/198.
nier_gems_rad_norm = total(nier_rad, 3)/temp

for iw=0, simrad_sz[2]-1 do begin 
  nier_rad[0:simrad_sz[0]-1, 0:simrad_sz[1]-1, iw] = $
      nier_rad[0:simrad_sz[0]-1, 0:simrad_sz[1]-1, iw] / nier_gems_rad_norm
endfor
;----------------------------------------------------------------------------
eosrl_rad = dblarr(simrad_sz[0], simrad_sz[1], simrad_sz[2])
for ix = 0, simrad_sz[0]-1 do begin
    for iy = 0, simrad_sz[1]-1 do begin
        eosrl_rad[ix, iy, 0:simrad_sz[2]-1] = eosrl_gemsl1crad.image_pixel_values[ix, iy, eosrl_startidxarr[ix, iy]:eosrl_startidxarr[ix, iy]+simrad_sz[2]-1]
    endfor
endfor

syntheticwav[where(syntheticwav le -1e29, /null)] = !values.f_nan

temp = fltarr(simrad_sz[0], simrad_sz[1])
for ix = 0, simrad_sz[0]-1 do begin
  temp[ix, *] = total(finite(syntheticwav), 2)
endfor

;gems_rad_norm = total(eosrl_rad, 3)/198.
eosrl_gems_rad_norm = total(eosrl_rad, 3)/temp

for iw=0, simrad_sz[2]-1 do begin 
  eosrl_rad[0:simrad_sz[0]-1, 0:simrad_sz[1]-1, iw] = $
      eosrl_rad[0:simrad_sz[0]-1, 0:simrad_sz[1]-1, iw] / eosrl_gems_rad_norm
endfor

;----------------------------------------------------------------------------
; NIER get irr more wider wavelength range for interpolation
;----------------------------------------------------------------------------

nier_irr0 = dblarr(simrad_sz[1], simrad_sz[2]+2)
nier_irr = dblarr(simrad_sz[1], simrad_sz[2])
nier_l1cirrwav0 = dblarr(simrad_sz[1], simrad_sz[2]+2)
nier_l1cirrwav = dblarr(simrad_sz[1], simrad_sz[2])

ix = floor(simrad_sz[0]/2.)
for iy = 0, simrad_sz[1]-1 do begin
    if nier_startidxarr[ix, iy] ne 0 then begin
        nier_irr0[iy, *] = nier_gemsirr.image_pixel_values[iy, nier_startidxarr[ix, iy]-1:nier_startidxarr[ix, iy]+simrad_sz[2]]
        nier_l1cirrwav0[iy, 0:simrad_sz[2]+1] = $
            nier_gemsirr.wavelength[iy, nier_startidxarr[ix, iy]-1:nier_startidxarr[ix, iy]+simrad_sz[2]]
    endif else begin 
        nier_irr0[iy, *] = nier_gemsirr.image_pixel_values[iy, 0:simrad_sz[2]+1]
        nier_l1cirrwav0[iy, 0:simrad_sz[2]+1] = $
            nier_gemsirr.wavelength[iy, 0:simrad_sz[2]+1]
    endelse
    nier_irr[iy, *] = interpol( $
        nier_irr0[iy, 0:simrad_sz[2]+1], nier_l1cirrwav0[iy, 0:simrad_sz[2]+1], nier_l1cradwav[ix, iy, *])
endfor



;gems_irrad_norm = total(irr, 2)/198.
nier_gems_irrad_norm = total(nier_irr, 2)/total(finite(syntheticwav), 2)

for iw=0, simrad_sz[2]-1 do begin 
  nier_irr[0:simrad_sz[1]-1, iw] = nier_irr[0:simrad_sz[1]-1, iw] / nier_gems_irrad_norm
endfor

irr3d = fltarr(simrad_sz[0], simrad_sz[1], simrad_sz[2])

for i=0, simrad_sz[0]-1 do begin
  irr3d[i, *, 0:simrad_sz[2]-1] = nier_irr
endfor

nier_gems_fitspec = fltarr(simrad_sz[0], simrad_sz[1], simrad_sz[2])

for ix=0,simrad_sz[0]-1 do begin
  for iw=0, simrad_sz[2]-1 do begin
    nier_gems_fitspec[ix, *, iw] = nier_rad[ix, *, iw]/irr3d[ix, *, iw]*nier_gems_rad_norm[ix, *]/nier_gems_irrad_norm[*]
  endfor
endfor
;----------------------------------------------------------------------------
; EOSRL get irr more wider wavelength range for interpolation
;----------------------------------------------------------------------------

eosrl_irr0 = dblarr(simrad_sz[1], simrad_sz[2]+2)
eosrl_irr = dblarr(simrad_sz[1], simrad_sz[2])
eosrl_l1cirrwav0 = dblarr(simrad_sz[1], simrad_sz[2]+2)
eosrl_l1cirrwav = dblarr(simrad_sz[1], simrad_sz[2])

ix = floor(simrad_sz[0]/2.)
for iy = 0, simrad_sz[1]-1 do begin
    if eosrl_startidxarr[ix, iy] ne 0 then begin
        eosrl_irr0[iy, *] = eosrl_gemsirr.image_pixel_values[iy, eosrl_startidxarr[ix, iy]-1:eosrl_startidxarr[ix, iy]+simrad_sz[2]]
        eosrl_l1cirrwav0[iy, 0:simrad_sz[2]+1] = $
            eosrl_gemsirr.wavelength[iy, eosrl_startidxarr[ix, iy]-1:eosrl_startidxarr[ix, iy]+simrad_sz[2]]
    endif else begin 
        eosrl_irr0[iy, *] = eosrl_gemsirr.image_pixel_values[iy, 0:simrad_sz[2]+1]
        eosrl_l1cirrwav0[iy, 0:simrad_sz[2]+1] = $
            eosrl_gemsirr.wavelength[iy, 0:simrad_sz[2]+1]
    endelse
    eosrl_irr[iy, *] = interpol( $
        eosrl_irr0[iy, 0:simrad_sz[2]+1], eosrl_l1cirrwav0[iy, 0:simrad_sz[2]+1], eosrl_l1cradwav[ix, iy, *])
endfor



;gems_irrad_norm = total(irr, 2)/198.
eosrl_gems_irrad_norm = total(eosrl_irr, 2)/total(finite(syntheticwav), 2)

for iw=0, simrad_sz[2]-1 do begin 
  eosrl_irr[0:simrad_sz[1]-1, iw] = eosrl_irr[0:simrad_sz[1]-1, iw] / eosrl_gems_irrad_norm
endfor

irr3d = fltarr(simrad_sz[0], simrad_sz[1], simrad_sz[2])

for i=0, simrad_sz[0]-1 do begin
  irr3d[i, *, 0:simrad_sz[2]-1] = eosrl_irr
endfor

eosrl_gems_fitspec = fltarr(simrad_sz[0], simrad_sz[1], simrad_sz[2])

for ix=0,simrad_sz[0]-1 do begin
  for iw=0, simrad_sz[2]-1 do begin
    eosrl_gems_fitspec[ix, *, iw] = eosrl_rad[ix, *, iw]/irr3d[ix, *, iw]*eosrl_gems_rad_norm[ix, *]/eosrl_gems_irrad_norm[*]
  endfor
endfor

;=============================================================================
; read L2 CLD
;=============================================================================


ecf = nier_gemscld.EffectiveCloudFraction
ecf02nanidx = where(ecf lt -1e29, /null)
_ecf = ecf
_ecf[ecf02nanidx] = !values.f_nan

;for ix=30, 30, 30 do begin
  ;for iy=30, 512, 50 do begin 
    ;outfilename = 'x' + strtrim(string(ix, format='(i03)'), 2) + $
      ;'_y' + strtrim(string(iy, format='(i03)'), 2) + '.png'

    ;plot_softcal_rad, syntheticwav, l1cradwav, $
      ;eosrl_gems_fitspec, synthetic_simrad, $
      ;ix, iy, suffix='adjusted'
      
  ;endfor
;endfor

nanidx = where(synthetic_simrad lt -1e+29, /null)
synthetic_simrad[nanidx] = !values.f_nan

diff = alog(eosrl_gems_fitspec) - alog(synthetic_simrad[*, *, 0:simrad_sz[2]-1])
gems_sz = size(diff, /dimension)

_diff = diff
diff_sz = size(_diff, /dimension)

_diff = reform(_diff, Long(gems_sz[0]) * gems_sz[1], gems_sz[2])
_diff[ecf02nanidx, *] = !values.f_nan
_diff = reform(_diff, gems_sz[0], gems_sz[1], gems_sz[2])

_diff_mean = mean(_diff, dimension=1, /nan)
diff_mean = fltarr(diff_sz[0], diff_sz[1], diff_sz[2])
for ix=0, diff_sz[0]-1 do begin
    diff_mean[ix, *, *] =  _diff_mean
endfor

_diff_median = median(_diff, dimension=1)
diff_median = fltarr(diff_sz[0], diff_sz[1], diff_sz[2])
for ix=0, diff_sz[0]-1 do begin
    diff_median[ix, *, *] =  _diff_median
endfor

diff_raw = diff
diff_raw[nanidx] = -1e30


; save netcdf
outncfn_raw = '/data/private/soodal/softcal_test/residual/gems_merra2_20210621_0345_raw.nc'

if file_test(outncfn_raw) then begin
  file_delete, outncfn_raw
endif

;=============================================================================
; READ MERRA2
;=============================================================================

merra2 = ds_read_merra2_tavg3_3d_asm_nv_collocated_on_gems_l1c($
  '/data/private/soodal/MERRA2_collocated_on_gems/' + $
  'merra2_tavg3_3d_asm_collocated_on_gemsl1c_rad_' + datetime_str + '.nc4')

pres = merra2.pl

itime = 1
o3mmr = reform(merra2.o3)
pres = reform(merra2.PL)
temp = reform(merra2.T)
alt = reform(merra2.H)

lon = merra2.lon
lat = merra2.lat

;convert mass mxing ratio into vmr
o3vmr = ds_mmr2vmr(o3mmr)
merra2_hpa = pres/100. ; [174, 512, 72]
;merra2_hpa_mid = merra2_hpa[*, *, 1:-1] -merra2_hpa[*, *, 0:-2]

ds_vmr2du_3d, o3vmr, merra2_hpa, temp, alt, merra2_o3, increasing_height=0
merra2_sz = size(merra2_o3, /dimension)

merra2_o3_accum = fltarr(merra2_sz)
merra2_o3_accum[*, *, 0] = merra2_o3[*, *, 0]
for i=1, merra2_sz[2]-1 do begin
  merra2_o3_accum[*, *, i] = merra2_o3_accum[*, *, i-1] + merra2_o3[*, *, i]
endfor

merra2_o3_on_gems_pres = fltarr(size(nier_data.pressure, /dimension))
merra2_o3_on_gems_o3 = fltarr(size(nier_data.o3, /dimension))

print, 'before the merra2 interpol loop.'
for iy=0, merra2_sz[1]-1 do begin
  for ix=0, merra2_sz[0]-1 do begin
    merra2_o3_on_gems_pres[ix, iy, *] = $
      interpol(merra2_o3_accum[ix, iy, *], merra2_hpa[ix, iy, 1:-1], nier_data.pressure[ix, iy, *])
  endfor
endfor

merra2_o3_on_gems_o3 = merra2_o3_on_gems_pres[*, *, 1:-1] - merra2_o3_on_gems_pres[*, *, 0:-2]

;===============================================================================
; write start raw
;===============================================================================
if write_flag then begin 
  id = ncdf_create(outncfn_raw)
  _ximage = ncdf_dimdef(id, 'image', diff_sz[0])
  _yspatial = ncdf_dimdef(id, 'spatial', diff_sz[1])
  _wavelength = ncdf_dimdef(id, 'wavelengths', diff_sz[2])

  ; define variables
  _fitspec_difference_from_merra2 = ncdf_vardef(id, 'fitspec_difference_from_merra2', $
      [_ximage, _yspatial, _wavelength], /double)

  ncdf_attput, id, _fitspec_difference_from_merra2, 'units', '1'
  ncdf_attput, id, _fitspec_difference_from_merra2, 'long name', 'fitspec difference from merra2 ozone profile fitspec'

  _wavelength_data = ncdf_vardef(id, 'wavelengths', $
      [_yspatial, _wavelength], /double)

  ncdf_attput, id, _wavelength_data, 'units', '[nm]'
  ncdf_attput, id, _wavelength_data, 'long name', 'wavelengths'

  _sza_data = ncdf_vardef(id, 'solarzenithangle', $
      [_ximage, _yspatial], /double)

  ncdf_attput, id, _sza_data, 'units', 'degrees'
  ncdf_attput, id, _sza_data, 'long name', 'solarzenithangle'

  ; put file in data mode
  ncdf_control, id, /endef

  ; input data
  ncdf_varput, id, _fitspec_difference_from_merra2, diff_raw
  ncdf_varput, id, _wavelength_data, l1cradwav[*, 0:simrad_sz[2]-1] 
  ncdf_varput, id, _sza_data, sza
  ncdf_close, id
endif

ypos = 400
ds_loadct, 33

ctable = colortable([[255, 0, 0], $
  [0, 255, 0], $
  [0, 0, 255]], ncolors=diff_sz[0], /transpose)

nanidx = where(sza lt -1e29, /null)
sza[nanidx] = !values.f_nan
sza_min = min(sza[*, ypos], /nan)
sza_max = max(sza[*, ypos], /nan)
sza_max_cap = 40

nanidx = where(vza lt -1e29, /null)
vza[nanidx] = !values.f_nan
vza_min = min(vza[*, ypos], /nan)
vza_max = max(vza[*, ypos], /nan)
vza_max_cap = 40

relative_difference = (eosrl_gems_fitspec - synthetic_simrad)/synthetic_simrad * 100.

sza_list = fltarr(diff_sz[0])
vza_list = fltarr(diff_sz[0])

sza_list[*] = !values.f_nan
vza_list[*] = !values.f_nan

lon_list = fltarr(diff_sz[0])
lat_list = fltarr(diff_sz[0])

lon_list[*] = !values.f_nan
lat_list[*] = !values.f_nan

ix_list = intarr(diff_sz[0])
ix_list[*] = -999

coloridx_list = intarr(diff_sz[0])
ix_list[*] = 0

alb_list = fltarr(diff_sz[0])
alb_list[*] = !values.f_nan

cld_threshold = 0.2
cld_thres_str = string(cld_threshold, format='(f4.2)')

; multiple plot

cities_name = ['Singapore', 'Kuala Lumpur', 'Hanoi', 'Hong Kong', 'Naha', $
  'Pohang']

cities_lon = [103.9, 101.70, 105.8, 114.10, 127.7, 129.2]
cities_lat = [1.30, 2.7, 21.0, 22.3, 26.2, 36.0]

print, 'start plotting'

project = 'EOSRL'

sub = '310_340'
;if not keyword_set(outputpath) then begin
  outputpath = './plot/'
;endif

;if keyword_set(project) then begin
  outputpath = outputpath + project + '/'
;endif else begin
  ;project = ''
;endelse

;if keyword_set(sub) then begin
  outputpath = outputpath + sub + '/'
;endif else begin
  ;sub = ''
;endelse

outputpath = outputpath + datetime_str + '/'
outputpath = outputpath + 'point_profile/'

plot_cities = 1


pixindices = [65-1, 317-1]
;if plot_cities then begin
  ;for icity=0, n_elements(cities_lat) -1 do begin
    ;print, cities_name[icity]
    ;pixidx = search_closest_pixel(nier_lon, nier_lat, cities_lon[icity], cities_lat[icity], maxlimit=0.2 )

    ;if pixidx ne -999 then begin 
        ;pixindices = array_indices(nier_lon, pixidx)

      ; ===========================================================================
      ; plotting a pixel
      ; ===========================================================================

        nier_ozprof     = nier_o3p[pixindices[0], pixindices[1], *]
        o3apriori   = nier_data.O3Apriori[pixindices[0], pixindices[1], *]

        plot_margin = [0.18, 0.10, 0.10, 0.15]
        plot_xrange = [-10,60]
        plot_yrange = [1000, 1]

        nier_pres = reform(nier_data.pressure[pixindices[0], pixindices[1], *])
        nier_pres = exp((alog(nier_pres[0:-2]) + alog(nier_pres[1:-1]))/2.)

        eosrl_pres = reform(eosrl_data.pressure[pixindices[0], pixindices[1], *])
        eosrl_pres = exp((alog(eosrl_pres[0:-2]) + alog(eosrl_pres[1:-1]))/2.)


        print, 'p1'
        p1 = plot(nier_ozprof, nier_pres, /buffer, /ylog, dim=[500, 600], $
          axis_style=0, $
          margin=plot_margin, $
          xrange=plot_xrange, $
          yrange=plot_yrange, $
          ytitle='Pressure [hPa]', $
          xtitle='O3 [DU]');color=[0, 0, 0], linestyle=0, name='GEMS O3P')
        ;p1.title = cities_name[icity] + ' '+ datetime_str
        ;p1.color = [255, 0, 0]
        p1.color = [228, 160, 0]
        ;p1.yrange= plot_yrange
        p1.linestyle = 0
        ;p1.symbol= '+'
        p1.name = 'NIER O3 profile'

        print, 'p2'
        p2 = plot(o3apriori, nier_pres[0:23], /buffer, $
          ;/ylog, dim=[500, 600], $
          ;axis_style=0, $
          ;margin=plot_margin, $
          ;xrange=plot_xrange, $
          ;yrange=plot_yrange, $
          ;ytitle='Pressure [hPa]', $
          ;xtitle='O3 [DU]')
          /overplot);, $ color='#e4a000', linestyle=1, name='Ozonesonde with GEMS O3P AVGK', /overplot)
        ;p2.title = cities_name[icity] + ' ' + datetime_str
        ;p2.title.font_size=16
        ;p2.yrange = plot_yrange
        ;p2.color = [228, 160, 0]
        p2.color = [0, 0, 0]
        p2.linestyle = 0
        ;p2.symbol= 's'
        p2.name = 'A priori'

        eosrl_ozprof = eosrl_o3p[pixindices[0], pixindices[1], *]
        print, 'p3'
        p3 = plot(eosrl_ozprof, eosrl_pres[0:23], /buffer, /overplot);, $ color='#e4a000', linestyle=1, name='Ozonesonde with GEMS O3P AVGK', /overplot)
        p3.color = [86, 180, 232] ; skyblue
        ;p3.color = [0, 0, 255]
        p3.linestyle = 0
        ;p3.symbol= 's'
        p3.name = 'EOSRL O3 Profile'
        
        ;MERRA2 Profile
        ;p4 = plot(reform(merra2_o3_on_gems_o3[pixindices[0], pixindices[1], *]), $
          ;;merra2_o3_on_gems_pres[pixindices[0], pixindices[1], 0:-2], $
          ;nier_pres, $
          ;/overplot, /buffer)
        ;;p4.color = [0, 159, 115] ; bluish green
        ;p4.color = [255, 0, 0] ; red
        ;;p4.linestyle = 1
        ;p4.name = 'MERRA2 O3 Profile'

        synthetic_pres = reform(synthetic_o3p_p[pixindices[0], pixindices[1], *])
        synthetic_pres = exp((alog(synthetic_pres[0:-2]) + alog(synthetic_pres[1:-1]))/2.)
        ;MERRA2 Synthetic Profile
        ;p5 = plot(reform(synthetic_o3p[pixindices[0], pixindices[1], *]), $
          ;;synthetic_o3p_p[pixindices[0], pixindices[1], 0:23], $
          ;synthetic_pres, $
          ;/overplot, /buffer)
        ;;p5.color = [0, 250, 150]
        ;;p5.color = [240, 228, 66] ; yellow
        ;p5.color = [0, 0, 255] ; blue
        ;;p5.linestyle = 2
        ;p5.name = 'GEMS O3 from Synthetic'


        p_ = plot(fltarr(500), indgen(500)*3-100, 'gray', /buffer, /overplot, $
          axis_style=0, $
          color='gray')

        a1 = axis('x', $
          target = p2, $
          ;textpos=0, $
          major=4, $
          minor=8, $
          ;tickdir = 0, $
          tickvalues=[-10, 0, 10, 20, 30, 40, 50, 60], $
          title='O3 [DU]', $
          location='bottom') ;[min(p2.xrange), 0, 0])

        a2 = axis('y', $
          target = p2, $
          ;textpos=0, $
          major=5, $
          minor=7, $
          tickvalues=[1, 10, 100, 300, 500, 700, 1000], $
          ;tickdir = 0, $
          ;tickunits='Scientific', $
          title='Pressure [hPa]', $
          location='left');[0, min(p2.yrange), 0] );, $

        a3 = axis('x', $
          target = p2, $
          ;textpos=1, $
          major=4, $
          minor=8, $
          tickvalues=[-10, 0, 10, 20, 30, 40, 50, 60], $
          ;tickdir = 1, $
          location='top');[max(p2.xrange), 0, 0])

        a4 = axis('y', $
          target = p2, $
          ;textpos=1, $
          major=5, $
          minor=7, $
          showtext=0, $
          ;tickvalues=[1, 10, 100, 500, 1000], $
          ;tickdir = 1, $
          location='right');[0, max(p2.yrange), 0] );, $
          ;title = 'Temperature')


        ;TWMO, -0.002, temp, 0, 100000., pres*100, pres_tropo, temp_tropo, alt_tropo, 0

        print, 'leg'
        leg = legend(target=[p1, p2, p3], position=[59, 1.5],/data)

        print, 't1'
        t1 = text(0.5, 0.95, 'GEMS O3 Profile '+ datetime_str, $
          font_size=16, $
          /normal, $
          alignment=0.5, $
          vertical_alignment=0.5)

        lat_t = text(0.25, 0.78, 'LAT:' +string(nier_lat[pixindices[0], pixindices[1]], format='(f6.2)'), /normal)
        lon_t = text(0.25, 0.75, 'LON:' +string(nier_lon[pixindices[0], pixindices[1]], format='(f6.2)'), /normal)

        print, 'sav'
        outputpath = './plot/EOSRL/300_340/20210329_0345/point_profile/'
        ;p2.save, outputpath + '/x' + $
          ;string(pixindices[0], format='(i03)') + $
          ;'y' + string(pixindices[1], format='(i03)') + '.png'
        if not file_test(outputpath) then begin
          file_mkdir, outputpath
        endif
        p2.save, outputpath + '300340_me0.1_eosrl_x' + string(pixindices[0], format='(i03)') +  $
            '_y' + string(pixindices[1], format='(i03)') + '.png'
        p2.close


        cldstr = string(ecf[pixindices[0], pixindices[1]], format='(f4.2)')

        ;plot_softcal_rad, syntheticwav, l1cradwav, $
          ;eosrl_gems_fitspec, synthetic_simrad, $
          ;nier_gems_fitspec, nier_l1cradwav, $
          ;pixindices[0], pixindices[1], cldstr, suffix='_310340_adjusted'

        print, 'plot a file' 
    ;endif
  ;endfor
;endif

stop
;=====================================================
; Plot pixels ecf < 0.05
;=====================================================
; plot multiple pixel's Normalized Radiance for each ECF interval

icf_final = 0.05
for icf = 0.05, icf_final, 0.05 do begin
  ecfnanidx = where(ecf lt -1e29, /null)
  ecf[ecfnanidx] = !values.f_nan

  idx = where(ecf lt icf and ecf ge icf - 0.05 and $
    finite(reform(nier_gems_fitspec[*, *, 0])) eq 1 and $
    finite(reform(eosrl_gems_fitspec[*, *, 0])) eq 1, /null)

  if n_elements(idx) ge 1 then begin
    ;for ip = 0, n_elements(idx)-1 do begin
    for ip = 0, n_elements(idx)-1 do begin

      pixindices = array_indices(nier_lon, idx[ip])
      if pixindices[0] eq 91 then begin

    ; ===========================================================================
    ; plotting a pixel
    ; ===========================================================================

      nier_ozprof     = nier_o3p[pixindices[0], pixindices[1], *]
      o3apriori   = nier_data.O3Apriori[pixindices[0], pixindices[1], *]

      plot_margin = [0.18, 0.10, 0.10, 0.15]
      plot_xrange = [-10,60]
      plot_yrange = [1000, 1]

      nier_pres = reform(nier_data.pressure[pixindices[0], pixindices[1], *])
      nier_pres = exp((alog(nier_pres[0:-2]) + alog(nier_pres[1:-1]))/2.)

      eosrl_pres = reform(eosrl_data.pressure[pixindices[0], pixindices[1], *])
      eosrl_pres = exp((alog(eosrl_pres[0:-2]) + alog(eosrl_pres[1:-1]))/2.)


      print, 'p1'
      p1 = plot(nier_ozprof, nier_pres, /buffer, /ylog, dim=[500, 600], $
        axis_style=0, $
        margin=plot_margin, $
        xrange=plot_xrange, $
        ytitle='Pressure [hPa]', $
        xtitle='O3 [DU]');color=[0, 0, 0], linestyle=0, name='GEMS O3P')
      ;p1.title = cities_name[icity] + ' '+ datetime_str
      ;p1.color = [255, 0, 0]
      p1.color = [228, 160, 0]
      p1.yrange= plot_yrange
      p1.linestyle = 0
      ;p1.symbol= '+'
      p1.name = 'NIER O3 profile'

      print, 'p2'
      p2 = plot(o3apriori, nier_pres[0:23], /buffer, $
        ;/ylog, dim=[500, 600], $
        ;axis_style=0, $
        ;margin=plot_margin, $
        ;xrange=plot_xrange, $
        ;ytitle='Pressure [hPa]', $
        ;xtitle='O3 [DU]')
        /overplot);, $ color='#e4a000', linestyle=1, name='Ozonesonde with GEMS O3P AVGK', /overplot)
      t2 = text(0.5, 0.9, 'x: ' + string(pixindices[0], format='(i03)') + $
        ' y:' + string(pixindices[1], format='(i03)') + ' ' + datetime_str, $
        alignment=0.5, font_size=18)
      p2.yrange = plot_yrange
      ;p2.color = [228, 160, 0]
      p2.color = [0, 0, 0]
      p2.linestyle = 0
      ;p2.symbol= 's'
      p2.name = 'A priori'

      eosrl_ozprof = eosrl_o3p[pixindices[0], pixindices[1], *]
      print, 'p3'
      p3 = plot(eosrl_ozprof, eosrl_pres[0:23], /buffer, /overplot);, $ color='#e4a000', linestyle=1, name='Ozonesonde with GEMS O3P AVGK', /overplot)
      ;p3.color = [228, 160, 0]
      p3.color = [0, 0, 255]
      p3.linestyle = 0
      ;p3.symbol= 's'
      p3.name = 'EOSRL O3 Profile'
      
      ;MERRA2 Profile
      p4 = plot(reform(merra2_o3_on_gems_o3[pixindices[0], pixindices[1], *]), $
        ;merra2_o3_on_gems_pres[pixindices[0], pixindices[1], 0:-2], $
        nier_pres, $
        /overplot, /buffer)
      p4.color = [255, 0, 0]
      p4.linestyle = 1
      p4.name = 'MERRA2 O3 Profile'

      synthetic_pres = reform(synthetic_o3p_p[pixindices[0], pixindices[1], *])
      synthetic_pres = exp((alog(synthetic_pres[0:-2]) + alog(synthetic_pres[1:-1]))/2.)
      ;MERRA2 Synthetic Profile
      p5 = plot(reform(synthetic_o3p[pixindices[0], pixindices[1], *]), $
        ;synthetic_pres[pixindices[0], pixindices[1], *], $
        synthetic_pres, $
        /overplot, /buffer)
      p5.color = [0, 250, 150]
      p5.linestyle = 2
      p5.name = 'GEMS O3 from Synthetic'

      ;eosrl_ozprof = eosrl_o3p[pixindices[0], pixindices[1], *]
      ;p3 = plot(eosrl_ozprof, eosrl_pres[0:23], /buffer, /overplot);, $ color='#e4a000', linestyle=1, name='Ozonesonde with GEMS O3P AVGK', /overplot)
      ;;p3.color = [228, 160, 0]
      ;p3.color = [0, 0, 255]
      ;p3.linestyle = 0
      ;;p3.symbol= 's'
      ;p3.name = 'A priori'

      ;p3 = plot(omi_csonprof, nier_pres, /buffer, /overplot);, $ color='#e4a000', linestyle=1, name='Ozonesonde with GEMS O3P AVGK', /overplot)
      ;p3.color = [86, 180, 232]
      ;p3.linestyle = 2
      ;p3.symbol= '*'
      ;p3.name = 'Ozonesonde with GEMS O3P AVGK'

      ;tpres1 = [tpres, tpres]
      ;p4 = plot([0, 100], tpres1, /buffer, /overplot);, $ color='#56b4e8', linestyle=2, name='Ozonesonde', /overplot)
      ;p4.color = [86, 180, 232]
      ;p4.linestyle = 0
      ;;p4.symbol= 'D'
      ;p4.name = 'Tropopause'

      ;temp = reform(nier_data.temperature[xidx[ipix], yidx[ipix], *])
      ;p5 = plot(temp, nier_pres, /buffer, /current, $
        ;axis_style=0, $
        ;margin=plot_margin, $
        ;/ylog);, $ color='#56b4e8', linestyle=2, name='Ozonesonde', /overplot)
      ;p5.color = [0, 159, 115]
      ;p5.yrange= plot_yrange
      ;p5.xrange=[150, 310]
      ;p5.linestyle = 0
      ;p5.symbol= 'D'
      ;p5.name = 'Temperature'
      p_ = plot(fltarr(500), indgen(500)*3-100, 'gray', /buffer, /overplot, $
        axis_style=0, $
        color='gray')

        ;ytitle='Pressure [hPa]', $
        ;xtitle='O3 [DU]');color=[0, 0, 0], linestyle=0, name='GEMS O3P')
      
      a1 = axis('x', $
        target = p2, $
        ;textpos=0, $
        major=4, $
        minor=8, $
        ;tickdir = 0, $
        tickvalues=[-10, 0, 10, 20, 30, 40, 50, 60], $
        title='O3 [DU]', $
        location='bottom') ;[min(p2.xrange), 0, 0])

      a2 = axis('y', $
        target = p2, $
        ;textpos=0, $
        major=5, $
        minor=7, $
        tickvalues=[1, 10, 100, 300, 500, 700, 1000], $
        ;tickdir = 0, $
        ;tickunits='Scientific', $
        title='Pressure [hPa]', $
        location='left');[0, min(p2.yrange), 0] );, $

      a3 = axis('x', $
        target = p2, $
        ;textpos=1, $
        major=4, $
        minor=8, $
        tickvalues=[-10, 0, 10, 20, 30, 40, 50, 60], $
        ;tickdir = 1, $
        location='top');[max(p2.xrange), 0, 0])

      a4 = axis('y', $
        target = p2, $
        ;textpos=1, $
        major=5, $
        minor=7, $
        showtext=0, $
        ;tickvalues=[1, 10, 100, 500, 1000], $
        ;tickdir = 1, $
        location='right');[0, max(p2.yrange), 0] );, $
        ;title = 'Temperature')


      ;TWMO, -0.002, temp, 0, 100000., pres*100, pres_tropo, temp_tropo, alt_tropo, 0

      print, 'leg'
      leg = legend(target=[p1, p2, p3, p4, p5], position=[59, 1.5],/data)

      print, 't1'
      ;t1 = text(0.5, 0.95, 'GEMS O3 Profile '+ datetime_str, $
        ;font_size=16, $
        ;/normal, $
        ;alignment=0.5, $
        ;vertical_alignment=0.5)

      lat_t = text(0.25, 0.78, 'LAT:' +string(nier_lat[pixindices[0], pixindices[1]], format='(f6.2)'), /normal)
      lon_t = text(0.25, 0.75, 'LON:' +string(nier_lon[pixindices[0], pixindices[1]], format='(f6.2)'), /normal)

      print, 'sav'
      ;p2.save, outputpath + '/x' + $
        ;string(pixindices[0], format='(i03)') + $
        ;'y' + string(pixindices[1], format='(i03)') + '.png'

  ;if not keyword_set(outputpath) then begin
    ;outputpath = './plot/'
  ;endif

  ;if keyword_set(project) then begin
    ;outputpath = outputpath + project + '/'
  ;endif else begin
    ;project = ''
  ;endelse

  ;if keyword_set(sub) then begin
    ;outputpath = outputpath + sub + '/'
  ;endif else begin
    ;sub = ''
  ;endelse

  ;outputpath = outputpath + datetime_str + '/'
  ;outputpath = outputpath + 'point_profile/'

      p2.save, outputpath + 'profile_synthetic_310340_x' + string(pixindices[0], format='(i03)') +  $
          '_y' + string(pixindices[1], format='(i03)')  + '.png'
      p2.close


      cldstr = string(ecf[pixindices[0], pixindices[1]], format='(f4.2)')

      plot_softcal_rad, syntheticwav, l1cradwav, $
        eosrl_gems_fitspec, synthetic_simrad, $
        nier_gems_fitspec, nier_l1cradwav, $
        pixindices[0], pixindices[1], cldstr, suffix='_310340_adjusted'

      print, 'plot a file' 
      endif
    endfor
  endif
endfor


; plot multiple pixel's Normalized Radiance for each ECF interval
icf_final = 0.05
for icf = 0.05, icf_final, 0.05 do begin
  ecfnanidx = where(ecf lt -1e29, /null)
  ecf[ecfnanidx] = !values.f_nan

  idx = where(ecf lt icf and ecf ge icf - 0.05 and $
    finite(reform(nier_gems_fitspec[*, *, 0])) eq 1 and $
    finite(reform(eosrl_gems_fitspec[*, *, 0])) eq 1, /null)

  if n_elements(idx) ge 1 then begin
    ;for ip = 0, n_elements(idx)-1 do begin
    for ip = 0, n_elements(idx)-1 do begin
      print, 'ip:', ip
      ;if ip gt 20 then begin
        ;break
      ;endif
      ;if ip eq 0 then begin
        ;op = 0
      ;endif else begin
        ;op = 1
      ;endelse
      op = 0
      indices = array_indices(ecf, idx[ip])
      if indices[0] eq 91 then begin
        yrange = [0, 0.25]
        ;p1 = plot(l1cradwav[ypos, *], eosrl_fitspec[xpos, ypos, *], /buffer, $

        ; Plot 1, EOSRL
        ;command = "p1" + string(ip, format='(i03)') + " = plot(" + $
        command = "p1" + " = plot(" + $
          "syntheticwav[indices[1], *], eosrl_gems_fitspec[indices[0], indices[1], *], /buffer, " + $
          "overplot=op, " + $
          "axis_style=0, " + $
          "xtitle='Wavelength[nm]', " + $
          "ytitle='Radiance[W/cm^2/cm/sr]', " + $
          "title='GEMS L1C Radiances with Synthetic Radiances', " + $
          ;symbol='s',  + $
          "yrange=yrange, " + $
          "xrange=[300, 340], " + $
          "name='GEMS EOSRL L1C Radiances', " + $
          "color='blue', " + $
          "transparency=50, " + $
          "position=pos)"
        print, command
        dummy = execute(command)

        ;command = "p1['axis1'].color='black'"
        ;dummy = execute(command)

        ;;plot simulated radiance
        ;command = "p2" + string(ip, format='(i03)') + " = plot(" + $
        command = "p2" + " = plot(" + $
          ;l1cradwav[ypos, *], synthetic_simrad[xpos, ypos, *], " + $
          "syntheticwav[indices[1], *], synthetic_simrad[indices[0], indices[1], *], " + $
          "/buffer, " + $
          "axis_style=0, " + $
          "name='MERRA2 Synthetic Radiance', " + $
          ;symbol='s', + $
          ;yrange=yrange, + $
          ;"yrange=yrange, " + $
          ;xrange=[300, 340], + $
          "position=pos, " + $
          "color='black', " + $
          "transparency=50, " + $
          "/overplot)"
          ;/current)
        dummy = execute(command)

        ; Plot 3, NIER 
        ;command = "p3" + string(ip, format='(i03)') + " = plot(" + $
        command = "p3" + " = plot(" + $
            "reform(nier_l1cradwav[indices[0], indices[1], *]), reform(nier_gems_fitspec[indices[0], indices[1], *]), " + $
            "/buffer, " + $
            "axis_style=0, " + $
            ;"yrange=yrange, " + $
            "name='NIER Normalized Radiance', " + $
            "color='red', " + $
            "/overplot)"
        print, command
        dummy = execute(command)


        ;yaxis2 = axis('Y', LOCATION='right', yrange=p2.yrange, TARGET=p2, $
          ;color='blue')
        ;xaxis2 = axis('X', LOCATION='TOP', xrange=p1.xrange, TARGET=p2, $
          ;color='black', $
          ;tickname=['', '', '', '', ''])

        ; ratio plot radiance

        eosrl_data = (eosrl_gems_fitspec[indices[0], indices[1], *] - $
          synthetic_simrad[indices[0], indices[1], *])/$
            synthetic_simrad[indices[0], indices[1], *]*100.

        nier_data = (nier_gems_fitspec[indices[0], indices[1], *] - $
          synthetic_simrad[indices[0], indices[1], *])/$
            synthetic_simrad[indices[0], indices[1], *]*100

        p_ = plot(indgen(500), fltarr(500), 'gray', /buffer, /current, $
          axis_style=0, $
          position=pos, $
          yrange=yrange, $
          xrange=[300, 340], $
          color='gray')

        ;idx = where(data lt 0, /null)
        ;data[idx] = 0
        p4 = barplot($
          syntheticwav[indices[1], *], eosrl_data, $
          /buffer, $
          axis_style=0, $
          name='GEMS EOSRL fitspec - Synthetic Ratio[%]', $
          ;symbol='s', $
          yrange=[-100, 100], $
          xrange=[300, 340], $
          position=pos, $
          fill_color='blue', $
          transparency=70, $
          outline=0, $
          /current)

        ;data = (eosrl_fitspec[xpos, ypos, *] - synthetic_simrad[xpos, ypos, *])/$
            ;synthetic_simrad[xpos, ypos, *]*100
        ;idx = where(data gt 0, /null)
        ;data[idx] = 0
        p5 = barplot($
          syntheticwav[indices[1], *], nier_data, $
          /buffer, $
          axis_style=0, $
          name='GEMS NIER fitspec - Synthetic Ratio[%]', $
          ;symbol='s', $
          yrange=[-100, 100], $
          xrange=[300, 340], $
          position=pos, $
          fill_color='red', $
          transparency=70, $
          outline=0, $
          /current)

        yaxis1 = axis('Y', LOCATION='left', yrange=yrange, TARGET=p1, $
          color='black', title='Normalized Radiance')
        xaxis1 = axis('X', LOCATION='BOTTOM', xrange=p1.xrange, TARGET=p1, $
          color='black', title='Wavelength[nm]');, $
          ;tickname=['', '', '', '', ''])

        yaxis2 = axis('Y', LOCATION='right', yrange=[-100, 100], TARGET=p4, $
          color='red', $
          title='Relative Difference [%]')
        xaxis2 = axis('X', LOCATION='TOP', xrange=p1.xrange, TARGET=p4, $
          color='black', $
          tickname=['', '', '', '', ''])

        ;p4 = plot($
          ;o3abs_wav, refabs1, '.', $
          ;/buffer, $
          ;axis_style=0, $
          ;name='O3 Absorption Cross section', $
          ;xrange=[300, 340], $
          ;;yrange=[0, 0.01], $
          ;yrange=[min(refabs1), max(refabs1)/20.], $
          ;position=pos, $
          ;color='black', $
          ;/current)

        leg = legend(target=[p1, p2, p3, p4, p5], position=[0.90, 0.85], $
          /normal, /auto_text_color)

        sIcf0 = String(icf-0.05, format='(f4.2)')
        sIcf1 = String(icf, format='(f4.2)')

        outfilepath = './plot/synthetic_radiance_diff/'
        pngfile = outfilepath + 'gems_merra2_synthetic_simrad_ecf' + sicf0 +'-' + sicf1 + $
          '_x' + string(indices[0], format='(i03)') + $
          '_y' + string(indices[1], format='(i03)') + '.png'
        ;p1.errorbar_color='indian_red'
        ;p1.errorbar_capsize=0.1
        ;p1.title.font_size=16
        p_.save, pngfile
        p_.close
      endif
    endfor
  endif

      ;;plot simulated radiance
  ;command = "p2" + string(ip, format='(i03)') + " = plot(" + $
    ;;l1cradwav[ypos, *], synthetic_simrad[xpos, ypos, *], " + $
    ;"syntheticwav[indices[1], *], synthetic_simrad[indices[0], indices[1], *], " + $
    ;"/buffer, " + $
    ;;";axis_style=0, " + $
    ;"name='MERRA2 Synthetic Radiance', " + $
    ;;symbol='s', + $
    ;;yrange=yrange, + $
    ;"yrange=yrange, " + $
    ;"xrange=[300, 340], " + $
    ;"position=pos, " + $
    ;"color='black', " + $
    ;"/overplot)"
    ;;/current)
  ;dummy = execute(command)

  ;p_ = plot(indgen(500), fltarr(500), 'gray', /buffer, /current, $
    ;axis_style=0, $
    ;position=pos, $
    ;yrange=yrange, $
    ;xrange=[300, 340], $
    ;color='gray')

  ;leg = legend(target=[p1, p2, p3], position=[325, 0.24], $
    ;/data, /auto_text_color)
  ;ecftext = text(330, 0.22, 'ECF = ' + cldstr, /data)

  ;p3 = plot(wav[0:fitresdims[0]-1], stdres, 'r', /buffer, $
    ;/over, name='Stddev')

  ;t1 = text(0.25, 0.75, 'Effective Cloud Fraction < 0.2', /norm)
  ;t2 = text(0.25, 0.70, 'Latitude < 20', /norm)
  ;leg = legend(target=[p0, p1, p2, p3, pap], pos=leg_pos)

  ;sxidx = string(xpos, format='(i03)')
  ;syidx = string(ypos, format='(i03)')
  ;sIcf0 = String(icf-0.05, format='(f4.2)')
  ;sIcf1 = String(icf, format='(f4.2)')

  ;outfilepath = './plot/synthetic_radiance_diff/'
  ;pngfile = outfilepath + 'gems_merra2_synthetic_simrad_ecf' + sicf0 +'-' + sicf1 + '.png'
  ;;p1.errorbar_color='indian_red'
  ;;p1.errorbar_capsize=0.1
  ;;p1.title.font_size=16
  ;p_.save, pngfile
  ;p_.close

endfor
      ; send image to pc
      ;scp_dest = 'soodal@164.125.38.179:/home/soodal/works/plot/'
      ;spawn, 'scp -P18742 -p ' + pngfile +  ' ' + scp_dest


;=========================================================================
; Plot for Clear pixel only MERRA2 with apriori
;=========================================================================

icf_final = 0.0
for icf = 0.05, icf_final, 0.05 do begin
  ecfnanidx = where(ecf lt -1e29, /null)
  ecf[ecfnanidx] = !values.f_nan

  idx = where(ecf lt icf and ecf ge icf - 0.05 and $
    finite(reform(nier_gems_fitspec[*, *, 0])) eq 1 and $
    finite(reform(eosrl_gems_fitspec[*, *, 0])) eq 1, /null)

  if n_elements(idx) ge 1 then begin
    ;for ip = 0, n_elements(idx)-1 do begin
    for ip = 0, n_elements(idx)-1 do begin
      print, 'ip:', ipk
      ;if ip gt 20 then begin
        ;break
      ;endif
      ;if ip eq 0 then begin
        ;op = 0
      ;endif else begin
        ;op = 1
      ;endelse
      op = 0
      indices = array_indices(ecf, idx[ip])
      yrange = [0, 0.25]
      ;p1 = plot(l1cradwav[ypos, *], eosrl_fitspec[xpos, ypos, *], /buffer, $

      ; Plot 1, EOSRL
      ;command = "p1" + string(ip, format='(i03)') + " = plot(" + $
      command = "p1" + " = plot(" + $
        "syntheticwav[indices[1], *], eosrl_gems_fitspec[indices[0], indices[1], *], /buffer, " + $
        "overplot=op, " + $
        "axis_style=0, " + $
        "xtitle='Wavelength[nm]', " + $
        "ytitle='Radiance[W/cm^2/cm/sr]', " + $
        "title='GEMS L1C Radiances with Synthetic Radiances', " + $
        ;symbol='s',  + $
        "yrange=yrange, " + $
        "xrange=[300, 340], " + $
        "name='GEMS EOSRL L1C Radiances', " + $
        "color='blue', " + $
        "transparency=50, " + $
        "position=pos)"
      print, command
      dummy = execute(command)

      ;command = "p1['axis1'].color='black'"
      ;dummy = execute(command)

      ;;plot simulated radiance
      ;command = "p2" + string(ip, format='(i03)') + " = plot(" + $
      command = "p2" + " = plot(" + $
        ;l1cradwav[ypos, *], synthetic_simrad[xpos, ypos, *], " + $
        "syntheticwav[indices[1], *], synthetic_simrad[indices[0], indices[1], *], " + $
        "/buffer, " + $
        "axis_style=0, " + $
        "name='MERRA2 Synthetic Radiance', " + $
        ;symbol='s', + $
        ;yrange=yrange, + $
        ;"yrange=yrange, " + $
        ;xrange=[300, 340], + $
        "position=pos, " + $
        "color='black', " + $
        "transparency=50, " + $
        "/overplot)"
        ;/current)
      dummy = execute(command)

      ; Plot 3, NIER 
      ;command = "p3" + string(ip, format='(i03)') + " = plot(" + $
      command = "p3" + " = plot(" + $
          "reform(nier_l1cradwav[indices[0], indices[1], *]), reform(nier_gems_fitspec[indices[0], indices[1], *]), " + $
          "/buffer, " + $
          "axis_style=0, " + $
          ;"yrange=yrange, " + $
          "name='NIER Normalized Radiance', " + $
          "color='red', " + $
          "/overplot)"
      print, command
      dummy = execute(command)


      ;yaxis2 = axis('Y', LOCATION='right', yrange=p2.yrange, TARGET=p2, $
        ;color='blue')
      ;xaxis2 = axis('X', LOCATION='TOP', xrange=p1.xrange, TARGET=p2, $
        ;color='black', $
        ;tickname=['', '', '', '', ''])

      ; ratio plot radiance

      eosrl_data = (eosrl_gems_fitspec[indices[0], indices[1], *] - $
        synthetic_simrad[indices[0], indices[1], *])/$
          synthetic_simrad[indices[0], indices[1], *]*100.

      nier_data = (nier_gems_fitspec[indices[0], indices[1], *] - $
        synthetic_simrad[indices[0], indices[1], *])/$
          synthetic_simrad[indices[0], indices[1], *]*100

      p_ = plot(indgen(500), fltarr(500), 'gray', /buffer, /current, $
        axis_style=0, $
        position=pos, $
        yrange=yrange, $
        xrange=[300, 340], $
        color='gray')

      ;idx = where(data lt 0, /null)
      ;data[idx] = 0
      p4 = barplot($
        syntheticwav[indices[1], *], eosrl_data, $
        /buffer, $
        axis_style=0, $
        name='GEMS EOSRL fitspec - Synthetic Ratio[%]', $
        ;symbol='s', $
        yrange=[-100, 100], $
        xrange=[300, 340], $
        position=pos, $
        fill_color='blue', $
        transparency=70, $
        outline=0, $
        /current)

      ;data = (eosrl_fitspec[xpos, ypos, *] - synthetic_simrad[xpos, ypos, *])/$
          ;synthetic_simrad[xpos, ypos, *]*100
      ;idx = where(data gt 0, /null)
      ;data[idx] = 0
      p5 = barplot($
        syntheticwav[indices[1], *], nier_data, $
        /buffer, $
        axis_style=0, $
        name='GEMS NIER fitspec - Synthetic Ratio[%]', $
        ;symbol='s', $
        yrange=[-100, 100], $
        xrange=[300, 340], $
        position=pos, $
        fill_color='red', $
        transparency=70, $
        outline=0, $
        /current)

      yaxis1 = axis('Y', LOCATION='left', yrange=yrange, TARGET=p1, $
        color='black', title='Normalized Radiance')
      xaxis1 = axis('X', LOCATION='BOTTOM', xrange=p1.xrange, TARGET=p1, $
        color='black', title='Wavelength[nm]');, $
        ;tickname=['', '', '', '', ''])

      yaxis2 = axis('Y', LOCATION='right', yrange=[-100, 100], TARGET=p4, $
        color='red', $
        title='Relative Difference [%]')
      xaxis2 = axis('X', LOCATION='TOP', xrange=p1.xrange, TARGET=p4, $
        color='black', $
        tickname=['', '', '', '', ''])

      ;p4 = plot($
        ;o3abs_wav, refabs1, '.', $
        ;/buffer, $
        ;axis_style=0, $
        ;name='O3 Absorption Cross section', $
        ;xrange=[300, 340], $
        ;;yrange=[0, 0.01], $
        ;yrange=[min(refabs1), max(refabs1)/20.], $
        ;position=pos, $
        ;color='black', $
        ;/current)

      leg = legend(target=[p1, p2, p3, p4, p5], position=[0.90, 0.85], $
        /normal, /auto_text_color)

      sIcf0 = String(icf-0.05, format='(f4.2)')
      sIcf1 = String(icf, format='(f4.2)')

      outfilepath = './plot/synthetic_radiance_diff/'
      pngfile = outfilepath + 'gems_merra2_synthetic_simrad_ecf' + sicf0 +'-' + sicf1 + $
        '_x' + string(indices[0], format='(i03)') + $
        '_y' + string(indices[1], format='(i03)') + '.png'
      ;p1.errorbar_color='indian_red'
      ;p1.errorbar_capsize=0.1
      ;p1.title.font_size=16
      p_.save, pngfile
      p_.close
    endfor
  endif

      ;;plot simulated radiance
  ;command = "p2" + string(ip, format='(i03)') + " = plot(" + $
    ;;l1cradwav[ypos, *], synthetic_simrad[xpos, ypos, *], " + $
    ;"syntheticwav[indices[1], *], synthetic_simrad[indices[0], indices[1], *], " + $
    ;"/buffer, " + $
    ;;";axis_style=0, " + $
    ;"name='MERRA2 Synthetic Radiance', " + $
    ;;symbol='s', + $
    ;;yrange=yrange, + $
    ;"yrange=yrange, " + $
    ;"xrange=[300, 340], " + $
    ;"position=pos, " + $
    ;"color='black', " + $
    ;"/overplot)"
    ;;/current)
  ;dummy = execute(command)

  ;p_ = plot(indgen(500), fltarr(500), 'gray', /buffer, /current, $
    ;axis_style=0, $
    ;position=pos, $
    ;yrange=yrange, $
    ;xrange=[300, 340], $
    ;color='gray')

  ;leg = legend(target=[p1, p2, p3], position=[325, 0.24], $
    ;/data, /auto_text_color)
  ;ecftext = text(330, 0.22, 'ECF = ' + cldstr, /data)

  ;p3 = plot(wav[0:fitresdims[0]-1], stdres, 'r', /buffer, $
    ;/over, name='Stddev')

  ;t1 = text(0.25, 0.75, 'Effective Cloud Fraction < 0.2', /norm)
  ;t2 = text(0.25, 0.70, 'Latitude < 20', /norm)
  ;leg = legend(target=[p0, p1, p2, p3, pap], pos=leg_pos)

  ;sxidx = string(xpos, format='(i03)')
  ;syidx = string(ypos, format='(i03)')
  ;sIcf0 = String(icf-0.05, format='(f4.2)')
  ;sIcf1 = String(icf, format='(f4.2)')

  ;outfilepath = './plot/synthetic_radiance_diff/'
  ;pngfile = outfilepath + 'gems_merra2_synthetic_simrad_ecf' + sicf0 +'-' + sicf1 + '.png'
  ;;p1.errorbar_color='indian_red'
  ;;p1.errorbar_capsize=0.1
  ;;p1.title.font_size=16
  ;p_.save, pngfile
  ;p_.close

endfor

stop

;==============================================================================
; plot wav with fitspec
;==============================================================================

for icity = 0, n_elements(cities_name)-1 do begin 
  print, icity
  ;data = reform((gems_fitspec[ix, ypos, *] - synthetic_simrad[ix, ypos, *])/$
      ;synthetic_simrad[ix, ypos, *]*100)

  data = reform(relative_difference[ix, ypos, *])

  ;if sza[ix, ypos] lt sza_max_cap and vza[ix, ypos] lt vza_max and _ecf[ix, ypos] lt cld_threshold then begin
  ;if max(abs(data[50:203])) lt 5. then begin 
    ;print, 'plot for ix:', ix
    ;print, data
    ;print, syntheticwav[ypos, *]
    if total(finite(lon[ix-1:ix, ypos])) eq 2 and total(finite(lat[ix-1:ix, ypos])) eq 2 then begin

      ix_list[ix] = 1

      print, sza[ix, ypos], vza[ix, ypos]
      sza_list[ix] = sza[ix, ypos]
      vza_list[ix] = vza[ix, ypos]

      lon_list[ix] = lon[ix, ypos]
      lat_list[ix] = lat[ix, ypos]

      edgelons = fltarr(2)
      edgelats = fltarr(2)

      edgelons[0] = lon[ix, ypos]
      edgelons[1] = lon[ix-1, ypos]

      edgelats[0] = lat[ix, ypos]
      edgelats[1] = lat[ix-1, ypos]

      get_omler_specific_time_loc, omler, month, day, edgelons, edgelats, alb
      alb_list[ix] = alb

      ; 0~174
      ;coloridx = floor((sza[ix, ypos] - sza_min)/(sza_max - sza_min) *(diff_sz[0]-1))
      lat_min = min(lat[*, ypos], /nan)
      lat_max = max(lat[*, ypos], /nan)
      lon_min = min(lon[*, ypos], /nan)
      lon_max = max(lon[*, ypos], /nan)
      coloridx = floor((lon[ix, ypos] - lon_min)/(lon_max - lon_min) *(diff_sz[0]-1))
      coloridx_list[ix] = coloridx
      command = $
      "p3" + string(ix, format='(i03)') + " = plot(" + $
        "reform(syntheticwav[ypos, *]), reform(data), " + $
        "/buffer, " + $
        "axis_style=2, " + $
        ;"title=datetime_str, " + $
        "name='GEMS fitspec - Synthetic Ratio[%]', " + $
        "xtitle='Wavelength[nm]', " + $
        "ytitle='Normalized Radiance Relative Difference[%]', " + $
        ;symbol='s',  + $
        "yrange=[-100, 100], " + $
        "xrange=[300, 340], " + $
        "position=pos, " + $
        ;fill_color='red', + $
        ;"color=ctable[*, ix], " + $
        "color=ctable[*, coloridx], " + $
        "transparency=50, " + $
        "overplot=ix) "
      dummy = execute(command)
        ;outline=0, " + $
      ;"p3" + string(ix, format='(i03)') + " = plot(" + 
        ;reform(syntheticwav[ypos, *]), reform(data), $
        ;;/buffer, $
        ;axis_style=2, $
        ;name='GEMS fitspec - Synthetic Ratio[%]', $
        ;xtitle='Wavelength[nm]', $
        ;ytitle='Normalized Radiance Relative Difference', $
        ;;symbol='s', $
        ;yrange=[-200, 200], $
        ;xrange=[300, 340], $
        ;position=pos, $
        ;;fill_color='red', $
        ;transparency=50, $
        ;overplot=ix);, $
        ;;outline=0, $
      ;endif
    endif
  ;endif

endfor


data = reform((gems_fitspec[173, ypos, *] - synthetic_simrad[173, ypos, *])/$
    synthetic_simrad[173, ypos, *]*100)
;if max(abs(data[50:203])) lt 5. then begin 

if sza[ix, ypos] lt sza_max_cap and vza[ix, ypos] lt vza_max_cap and _ecf[ix, ypos] < cld_threshold then  begin
    ix_list[ix] = 1

    print, sza[ix, ypos], vza[ix, ypos]
    sza_list[ix] = sza[ix, ypos]
    vza_list[ix] = vza[ix, ypos]

  p3_1 = plot($
    reform(syntheticwav[ypos, *]), reform(data), $
    /buffer, $
    ;axis_style=0, $
    name='GEMS fitspec - Synthetic Ratio[%]', $
    xtitle='Wavelength[nm]', $
    ytitle='Normalized Radiance Relative Difference[%]', $
    ;symbol='s', $
    yrange=[-100, 100], $
    xrange=[300, 340], $
    position=pos, $
    ;fill_color='red', $
    transparency=50, $
    /overplot);, $
    ;outline=0, $
endif

p_ = plot(indgen(500), fltarr(500), 'gray', /buffer, /overplot, $
  ;axis_style=0, $
  title='Clear Case', $
  position=pos, $
  yrange=[-100, 100], $
  xrange=[300, 340], $
  color='gray')

t = text(320, 65, /data, 'spatial_index=' + string(ypos))

    ;coloridx = floor((sza[ix, ypos] - sza_min)/(sza_max - sza_min) *diff_sz[0]-1)

;for isza=0, 10 do begin
  ;sza_str = string((sza_max - sza_min)* isza/10. + sza_min, format='(f7.3)')
  ;coloridx = floor(isza/10. *(diff_sz[0]-1))
  ;t = text(0.82, isza*0.05+0.1,color=ctable[*, coloridx], /norm, 'sza=' + sza_str)
;endfor

;==============================================================================
; write text to the plot
;==============================================================================
;ix_idx = where(ix_list eq 1, ix_idx_num, /null)
;for ix_idx_i = 0, ix_idx_num-1 do begin
  ;sza_str = string(sza[ix_idx[ix_idx_i], ypos], format='(f7.3)')
  ;vza_str = string(vza[ix_idx[ix_idx_i], ypos], format='(f7.3)')
  ;alb_str = string(alb_list[ix_idx[ix_idx_i]], format='(f5.2)')
  ;ecf_str = string(_ecf[ix_idx[ix_idx_i], ypos], format='(f4.2)')
  ;t_sza = text(305, ix_idx_i*10-95, color=ctable[*, coloridx_list[ix_idx[ix_idx_i]]], /data, 'sza=' + sza_str)
  ;t_vza = text(314, ix_idx_i*10-95, color=ctable[*, coloridx_list[ix_idx[ix_idx_i]]], /data, 'vza=' + vza_str)
  ;t_alb = text(323, ix_idx_i*10-95, color=ctable[*, coloridx_list[ix_idx[ix_idx_i]]], /data, 'alb=' + alb_str)
  ;t_alb = text(332, ix_idx_i*10-95, color=ctable[*, coloridx_list[ix_idx[ix_idx_i]]], /data, 'ecf=' + ecf_str)
;endfor

;==============================================================================


;cbpos = [0.92, 0.1, 0.95, 0.85]
;cb = colorbar(target=p_, orientation=1, position=cbpos)

t1 = text(320, 75, 'Effective Cloud Fraction < ' + cld_thres_str, /data)
;t2 = text(0.25, 0.70, 'Latitude < 20', /norm)
;leg = legend(target=[p0, p1, p2, p3, pap], pos=leg_pos)
;sxidx = string(ix, format='(i03)')
syidx = string(ypos, format='(i03)')
outfilepath = '/home/soodal/works/GEMS_O3P_analysis/plot/synthetic_radiance_diff/'
pngfile = outfilepath + 'gems_merra2_synthetic_simrad_y' +sYidx + '.png'
print, pngfile
;p1.errorbar_color='indian_red'
;p1.errorbar_capsize=0.1

;p3.title.font_size=16
p_.save, pngfile
p_.close


end
