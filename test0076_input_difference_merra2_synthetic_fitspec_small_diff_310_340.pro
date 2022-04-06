
; MERRA2 FITSPEC : 300.04 ~ 
; L1C, IRR : 297.~
; The Actually fitting range: 300.4~ 

pro plot_softcal_rad, syntheticwav, l1cradwav, $
  obs_fitspec, synthetic_simrad, $
  xpos, ypos, suffix=suffix
  
  if not keyword_set(suffix) then begin
      suffix = ''
  endif

  read_o3_abs, o3abs_wav, refabs1, refabs2, refabs3

  pos = [0.12, 0.1, 0.9, 0.85]
  leg_pos = [0.85, 0.8]

  outfilepath = './plot/synthetic_radiance_diff/'
  if not file_test(outfilepath) then begin
    file_mkdir, outfilepath
  endif

  yrange = [0, 0.25]
  ;p1 = plot(l1cradwav[ypos, *], obs_fitspec[xpos, ypos, *], /buffer, $
  p1 = plot(syntheticwav[ypos, *], obs_fitspec[xpos, ypos, *], /buffer, $
    axis_style=1, $
    xtitle='Wavelength[nm]', $
    ytitle='Radiance[W/cm^2/cm/sr]', $
    title='GEMS L1C Radiances with Synthetic Radiances', $
    ;symbol='s', $
    yrange=yrange, $
    xrange=[300, 340], $
    name='GEMS L1C Radiances', $
    color='red', $
    position=pos)
  p1['axis1'].color='red'

  ;plot simulated radiance
  p2 = plot($
    ;l1cradwav[ypos, *], synthetic_simrad[xpos, ypos, *], $
    syntheticwav[ypos, *], synthetic_simrad[xpos, ypos, *], $
    /buffer, $
    axis_style=0, $
    name='MERRA2 Synthetic Radiance', $
    ;symbol='s', $
    yrange=yrange, $
    xrange=[300, 340], $
    position=pos, $
    color='blue', $
    /current)

  ;yaxis2 = axis('Y', LOCATION='right', yrange=p2.yrange, TARGET=p2, $
    ;color='blue')
  xaxis2 = axis('X', LOCATION='TOP', xrange=p1.xrange, TARGET=p2, $
    color='black', $
    tickname=['', '', '', '', ''])

  ; ratio plot radiance

  data = (obs_fitspec[xpos, ypos, *] - synthetic_simrad[xpos, ypos, *])/$
      synthetic_simrad[xpos, ypos, *]*100

  idx = where(data lt 0, /null)

  data[idx] = 0

  p3 = barplot($
    syntheticwav[ypos, *], data, $
    /buffer, $
    axis_style=0, $
    name='GEMS fitspec - Synthetic Ratio[%]', $
    ;symbol='s', $
    yrange=[-50, 50], $
    xrange=[300, 340], $
    position=pos, $
    fill_color='red', $
    transparency=50, $
    outline=0, $
    /current)

  data = (obs_fitspec[xpos, ypos, *] - synthetic_simrad[xpos, ypos, *])/$
      synthetic_simrad[xpos, ypos, *]*100
  idx = where(data gt 0, /null)
  data[idx] = 0
  p3 = barplot($
    syntheticwav[ypos, *], data, $
    /buffer, $
    axis_style=0, $
    name='GEMS fitspec - Synthetic Ratio[%]', $
    ;symbol='s', $
    yrange=[-50, 50], $
    xrange=[300, 340], $
    position=pos, $
    fill_color='blue', $
    transparency=50, $
    outline=0, $
    /current)

  yaxis2 = axis('Y', LOCATION='right', yrange=[-2, 2], TARGET=p3, $
    color='red')

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

  p_ = plot(indgen(500), fltarr(500), 'gray', /buffer, /current, $
    axis_style=0, $
    position=pos, $
    yrange=p3.yrange, $
    xrange=p1.xrange, $
    color='gray')

  leg = legend(target=[p1, p2, p3], position=[325, 0.24], $
    /data, /auto_text_color)

  ;p3 = plot(wav[0:fitresdims[0]-1], stdres, 'r', /buffer, $
    ;/over, name='Stddev')

  ;t1 = text(0.25, 0.75, 'Effective Cloud Fraction < 0.2', /norm)
  ;t2 = text(0.25, 0.70, 'Latitude < 20', /norm)
  ;leg = legend(target=[p0, p1, p2, p3, pap], pos=leg_pos)
  sxidx = string(xpos, format='(i03)')
  syidx = string(ypos, format='(i03)')
  pngfile = outfilepath + 'gems_merra2_synthetic_simrad_x' +sXidx +'_' + suffix + '.png'
  p1.errorbar_color='indian_red'
  p1.errorbar_capsize=0.1
  p1.title.font_size=16
  p1.save, pngfile
  p1.close

  ; send image to pc
  ;scp_dest = 'soodal@164.125.38.179:/home/soodal/works/plot/'
  ;spawn, 'scp -P18742 -p ' + pngfile +  ' ' + scp_dest
end

ds_read_omler, omler

pos = [0.12, 0.1, 0.9, 0.85]

write_flag = 0
;=============================================================================
; start
;=============================================================================
jd_list = timegen(start=julday(8, 2, 2020, 3, 45), $
  final=julday(8, 2, 2020, 3, 45), units='Hours')

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
; SYNTHETIC
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
synthetic_simrad = synthetic.simulatedradiances[*, *, 0:197]

syntheticwav = synthetic.wavelengths[*, 0:197]


simrad_sz = size(synthetic_simrad, /dimension)
synthetic_fitspec = dblarr(simrad_sz[0], simrad_sz[1], simrad_sz[2])

for iw = 0, 197 do begin
  synthetic_fitspec[*, *, iw] = corrected_rad_all[*, *, iw] / corrected_irrad_all[*, *, iw] * div_rad_all / div_sun_all
endfor


; solar zenith angle
;sza = dblarr(simrad_sz[0], simrad_sz[1])
;for iw = 0, 197 do begin
    ;sza[*, *] = synthetic.solarzenithangle
;endfor

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
l1cradfn = '/data/private/soodal/ln/GEMS/L1C/GK2_GEMS_L1C_20210621_0345_4x4.nc'
l1cradfn = '/data/private/soodal/ln/GEMS/L1C/GK2_GEMS_L1C_' + datetime_str + '_4x4.nc'

;irrfn = '/data2/L1C_GEMS/L1C/4x4/GK2_GEMS_IRR_20210620_4x4.nc'
irrfn = '/data/private/soodal/ln/GEMS/IRR/GK2_GEMS_IRR_20210620_4x4.nc'
irrfn = '/data/private/soodal/ln/GEMS/IRR/GK2_GEMS_IRR_' + datetime_y_str + '_4x4.nc'

;cldfn = '/data/nier_ftp/CLOUD/V03/202106/21/GK2_GEMS_L2_20210621_0345_CLOUD_FW_DPRO_BIN4x4.nc'
cldfn = '/data/private/soodal/ln/GEMS/L2CLD/GK2_GEMS_L2_CLD_20210621_0345_4x4.nc'
cldfn = '/data/private/soodal/ln/GEMS/L2CLD/GK2_GEMS_L2_CLD_' + datetime_str + '_4x4.nc'

gemsl1crad = ds_read_gems_l1c_4x4(l1cradfn)
gemsirr = ds_read_gems_l1c_4x4(irrfn)
gemscld = ds_read_gems_l2_o3p(cldfn)

startidxarr = lonarr(simrad_sz[0], simrad_sz[1])

for ix = 0, simrad_sz[0]-1 do begin
    for iy = 0, simrad_sz[1]-1 do begin
        if syntheticwav[iy, 0] lt -1e29 then begin
            startidxarr[ix, iy] = 0
            continue
        endif

        idx = where(gemsl1crad.wavelength[iy, *] eq syntheticwav[iy, 0], /null)
        if n_elements(idx) eq 0 then begin
            stop
        endif
        startidxarr[ix, iy] = idx
    endfor
endfor

;l1cradwav = double(gemsl1crad.wavelength[*, 0+11+2:197+11+2])
;l1cradwav0 = double(gemsl1crad.wavelength[*, 0+11:197+11])
l1cradwav = dblarr(simrad_sz[0], simrad_sz[1], simrad_sz[2])
for ix = 0, simrad_sz[0]-1 do begin
    for iy = 0, simrad_sz[1]-1 do begin
        l1cradwav[ix, iy, *] = gemsl1crad.wavelength[iy, startidxarr[ix, iy]:startidxarr[ix, iy]+197]
    endfor
endfor
        


;rad = reform(gemsl1crad.image_pixel_values[*, *, 0+11+2:197+11+2])
rad = dblarr(simrad_sz[0], simrad_sz[1], simrad_sz[2])
for ix = 0, simrad_sz[0]-1 do begin
    for iy = 0, simrad_sz[1]-1 do begin
        rad[ix, iy, 0:simrad_sz[2]-1] = gemsl1crad.image_pixel_values[ix, iy, startidxarr[ix, iy]:startidxarr[ix, iy]+simrad_sz[2]-1]
    endfor
endfor

gems_rad_norm = total(rad, 3)/198.

for iw=0, simrad_sz[2]-1 do begin 
  rad[0:simrad_sz[0]-1, 0:simrad_sz[1]-1, iw] = $
      rad[0:simrad_sz[0]-1, 0:simrad_sz[1]-1, iw] / gems_rad_norm
endfor

;----------------------------------------------------------------------------
; get irr more wider wavelength range for interpolation
;----------------------------------------------------------------------------
;irr = reform(gemsirr.image_pixel_values[*, 0+11+2:197+11+2])


irr0 = dblarr(simrad_sz[1], simrad_sz[2]+2)
irr = dblarr(simrad_sz[1], simrad_sz[2])
l1cirrwav0 = dblarr(simrad_sz[1], simrad_sz[2]+2)
l1cirrwav = dblarr(simrad_sz[1], simrad_sz[2])

ix = floor(simrad_sz[0]/2.)
for iy = 0, simrad_sz[1]-1 do begin
    if startidxarr[ix, iy] ne 0 then begin
        irr0[iy, *] = gemsirr.image_pixel_values[iy, startidxarr[ix, iy]-1:startidxarr[ix, iy]+simrad_sz[2]]
        l1cirrwav0[iy, 0:simrad_sz[2]+1] = $
            gemsirr.wavelength[iy, startidxarr[ix, iy]-1:startidxarr[ix, iy]+simrad_sz[2]]
    endif else begin 
        irr0[iy, *] = gemsirr.image_pixel_values[iy, 0:simrad_sz[2]+1]
        l1cirrwav0[iy, 0:simrad_sz[2]+1] = $
            gemsirr.wavelength[iy, 0:simrad_sz[2]+1]
    endelse
    ;l1cirrwav[iy, *] = interpol(irr0[iy, *], l1cirrwav0[iy, *], l1cradwav[ix, iy, *])
    irr[iy, *] = interpol( $
        irr0[iy, 0:simrad_sz[2]+1], l1cirrwav0[iy, 0:simrad_sz[2]+1], l1cradwav[ix, iy, *])
endfor
;l1cirrwav0 = gemsirr.wavelength[*, 0+11:197+11]
;l1cirrwav = gemsirr.wavelength[*, 0+11+2:197+11+2]

gems_irrad_norm = total(irr, 2)/198.

for iw=0, simrad_sz[2]-1 do begin 
  irr[0:simrad_sz[1]-1, iw] = irr[0:simrad_sz[1]-1, iw] / gems_irrad_norm
endfor

irr3d = fltarr(simrad_sz[0], simrad_sz[1], simrad_sz[2])

for i=0, simrad_sz[0]-1 do begin
  irr3d[i, *, 0:simrad_sz[2]-1] = irr
endfor

gems_fitspec = fltarr(simrad_sz[0], simrad_sz[1], simrad_sz[2])

for ix=0,simrad_sz[0]-1 do begin
  for iw=0, simrad_sz[2]-1 do begin
    gems_fitspec[ix, *, iw] = rad[ix, *, iw]/irr3d[ix, *, iw]*gems_rad_norm[ix, *]/gems_irrad_norm[*]
  endfor
endfor


;=============================================================================


ecf = gemscld.EffectiveCloudFraction
;ecf02idx = where(ecf ge 0 and ecf le 0.2, /null)
;ecf02nanidx = where(ecf lt 0 or ecf gt 0.2, /null)
ecf02nanidx = where(ecf lt -1e29, /null)
_ecf = ecf
_ecf[ecf02nanidx] = !values.f_nan

;for ix=30, 174, 30 do begin
for ix=30, 30, 30 do begin
  for iy=30, 512, 50 do begin 
    outfilename = 'x' + strtrim(string(ix, format='(i03)'), 2) + $
      '_y' + strtrim(string(iy, format='(i03)'), 2) + '.png'

    ;plot_softcal_rad, syntheticwav, l1cradwav, $
      ;gems_fitspec, synthetic_simrad, $
      ;ix, iy, suffix='adjusted'
      
  endfor
endfor

nanidx = where(synthetic_simrad lt -1e+29, /null)
synthetic_simrad[nanidx] = !values.f_nan

diff = alog(gems_fitspec) - alog(synthetic_simrad[*, *, 0:simrad_sz[2]-1])
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

;diff_sz = size(diff_mean, /dimension)
; save netcdf
outncfn_mean = '/data/private/soodal/softcal_test/residual/gems_merra2_20210621_0345_mean.nc'
outncfn_median = '/data/private/soodal/softcal_test/residual/gems_merra2_20210621_0345_median.nc'
outncfn_raw = '/data/private/soodal/softcal_test/residual/gems_merra2_20210621_0345_raw.nc'

if file_test(outncfn_mean) then begin
  file_delete, outncfn_mean
endif
if file_test(outncfn_median) then begin
  file_delete, outncfn_median
endif
if file_test(outncfn_raw) then begin
  file_delete, outncfn_raw
endif

;===============================================================================
; write start mean
;===============================================================================
;id = ncdf_create(outncfn_mean)
;_ximage = ncdf_dimdef(id, 'image', diff_sz[0])
;_yspatial = ncdf_dimdef(id, 'spatial', diff_sz[1])
;_wavelength = ncdf_dimdef(id, 'wavelength', diff_sz[2])

;; define variables
;_fitspec_difference_from_merra2 = ncdf_vardef(id, 'fitspec_difference_from_merra2', $
    ;[_ximage, _yspatial, _wavelength], /double)

;ncdf_attput, id, _fitspec_difference_from_merra2, 'units', '1'
;ncdf_attput, id, _fitspec_difference_from_merra2, 'long name', 'fitspec difference from merra2 ozone profile fitspec'

;_wavelength_data = ncdf_vardef(id, 'wavelength', $
    ;[_yspatial, _wavelength], /double)

;ncdf_attput, id, _wavelength_data, 'units', '[nm]'
;ncdf_attput, id, _wavelength_data, 'long name', 'wavelength'

;_sza_data = ncdf_vardef(id, 'solarzenithangle', $
    ;[_ximage, _yspatial], /double)

;ncdf_attput, id, _sza_data, 'units', 'degrees'
;ncdf_attput, id, _sza_data, 'long name', 'solarzenithangle'

;; put file in data mode
;ncdf_control, id, /endef

;; input data
;ncdf_varput, id, _fitspec_difference_from_merra2, diff_mean
;ncdf_varput, id, _wavelength_data, l1cradwav[*, 0:simrad_sz[2]-1]
;ncdf_varput, id, _sza_data, sza
;ncdf_close, id

;===============================================================================
; write start median
;===============================================================================
;id = ncdf_create(outncfn_median)
;_ximage = ncdf_dimdef(id, 'image', diff_sz[0])
;_yspatial = ncdf_dimdef(id, 'spatial', diff_sz[1])
;_wavelength = ncdf_dimdef(id, 'wavelengths', diff_sz[2])

;; define variables
;_fitspec_difference_from_merra2 = ncdf_vardef(id, 'fitspec_difference_from_merra2', $
    ;[_ximage, _yspatial, _wavelength], /double)
;ncdf_attput, id, _fitspec_difference_from_merra2, 'units', '1'
;ncdf_attput, id, _fitspec_difference_from_merra2, 'long name', 'fitspec difference from merra2 ozone profile fitspec'

;_wavelength_data = ncdf_vardef(id, 'wavelengths', $
    ;[_yspatial, _wavelength], /double)

;ncdf_attput, id, _wavelength_data, 'units', '[nm]'
;ncdf_attput, id, _wavelength_data, 'long name', 'wavelengths'

;_sza_data = ncdf_vardef(id, 'solarzenithangle', $
    ;[_ximage, _yspatial], /double)

;ncdf_attput, id, _sza_data, 'units', 'degrees'
;ncdf_attput, id, _sza_data, 'long name', 'solarzenithangle'

;; put file in data mode
;ncdf_control, id, /endef

;; input data
;ncdf_varput, id, _fitspec_difference_from_merra2, diff_median
;ncdf_varput, id, _wavelength_data, l1cradwav[*, 0:simrad_sz[2]-1] 
;ncdf_varput, id, _sza_data, sza
;ncdf_close, id


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

; multiple plot

relative_difference = (gems_fitspec - synthetic_simrad)/synthetic_simrad * 100.

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

for ix = 0, diff_sz[0]-2 do begin 
  ; ratio plot radiance
  print, ix
  ;data = reform((gems_fitspec[ix, ypos, *] - synthetic_simrad[ix, ypos, *])/$
      ;synthetic_simrad[ix, ypos, *]*100)

  data = reform(relative_difference[ix, ypos, *])

  cap_310 = 5.
  ;if sza[ix, ypos] lt sza_max_cap and vza[ix, ypos] lt vza_max and _ecf[ix, ypos] lt 0.05 then begin
  ;if max(abs(data[50:197])) lt cap_310 and _ecf[ix, ypos] lt 0.2 then begin 
  if max(abs(data[50:197])) lt cap_310 then begin 
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
      coloridx = floor((sza[ix, ypos] - sza_min)/(sza_max - sza_min) *(diff_sz[0]-1))
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
        ;symbol='s', " + $
        "yrange=[-100, 100], " + $
        "xrange=[300, 340], " + $
        "position=pos, " + $
        ;fill_color='red', " + $
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
  endif

endfor


data = reform((gems_fitspec[173, ypos, *] - synthetic_simrad[173, ypos, *])/$
    synthetic_simrad[173, ypos, *]*100)
;if max(abs(data[50:197])) lt 5. then begin 
if sza[ix, ypos] and sza_max_cap and vza[ix, ypos] lt vza_max_cap and _ecf[ix, ypos] < 0.05 then  begin
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
  title='Small Relative Difference Case', $
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

;cbpos = [0.92, 0.1, 0.95, 0.85]
;cb = colorbar(target=p_, orientation=1, position=cbpos)

t1 = text(320, 75, 'Effective Cloud Fraction < 0.2', /data)
;t2 = text(0.25, 0.70, 'Latitude < 20', /norm)
;leg = legend(target=[p0, p1, p2, p3, pap], pos=leg_pos)
;sxidx = string(ix, format='(i03)')
syidx = string(ypos, format='(i03)')
outfilepath = '/home/soodal/works/GEMS_O3P_analysis/plot/synthetic_radiance_diff/'
pngfile = outfilepath + 'gems_merra2_synthetic_simrad_y' +sYidx + '_small_diff_notext.png'
print, pngfile
;p1.errorbar_color='indian_red'
;p1.errorbar_capsize=0.1

;p3.title.font_size=16
p_.save, pngfile
p_.close


end
