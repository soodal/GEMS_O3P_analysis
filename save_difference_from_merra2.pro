
;=============================================================================
; start
;=============================================================================
jd_list = timegen(start=julday(3, 29, 2021, 3, 45), $
  final=julday(3, 29, 2021, 3, 45), units='Hours')

for itime = 0, n_elements(jd_list)-1 do begin
  caldat, jd_list[itime], month, day, year, hour, minute
  yyyy = string(year, format='(i04)')
  mm = string(month, format='(i02)')
  dd = string(day, format='(i02)')
  hh = string(hour, format='(i02)')
  mi = string(minute, format='(i02)')

  datetime_str = yyyy + mm + dd + '_' + hh + mi

  caldat, jd_list[itime]-1, month_y, day_y, year_y, hour_y, minute_y
  yyyy_y = string(year_y, format='(i04)')
  mm_y = string(month_y, format='(i02)')
  dd_y = string(day_y, format='(i02)')
  hh_y = string(hour_y, format='(i02)')
  mi_y = string(minute_y, format='(i02)')

  datetime_str = yyyy + mm + dd + '_' + hh + mi
  print, datetime_str

  ;=============================================================================
  ; SYNTHETIC
  ;=============================================================================
  search_str = '/home/soodal/data/merra2_residual/residuals/' + $
    'GK2_GEMS_L2_O3P_' + datetime_str + '_wli300_prec000000_climML_b4x4_*_ecf0.nc'

  fl = file_search(search_str)
  if strlen(file_basename(fl[0])) gt 1 then begin 
    print, fl[0]
    gems_exist = 1
  endif else begin
    gems_exist = 0
  endelse

  ;=============================================================================
  ; GEMS Observation
  ;=============================================================================
  ;l1cradfn_search_str = '/data2/L1C_GEMS/L1C/4x4/GK2_GEMS_L1C_' + datetime_str + '_*_4x4.nc'
  l1cradfn_search_str = '/data/private/soodal/ln/GEMS/L1C/GK2_GEMS_L1C_' + datetime_str + '_4x4.nc'
  l1cradfn_list = file_search(l1cradfn_search_str)
  l1cradfn = l1cradfn_list[0]

  if strlen(file_basename(l1cradfn)) gt 1 then begin 
    print, l1cradfn[0]
    l1crad_exist = 1
  endif else begin
    l1crad_exist = 0
  endelse

  ;irrfn_search_str = '/data2/L1C_GEMS/L1C/4x4/GK2_GEMS_IRR_' + yyyy_y + mm_y + dd_y + '_4x4.nc'
  irrfn_search_str = '/data/private/soodal/ln/GEMS/IRR/GK2_GEMS_IRR_' + yyyy_y + mm_y + dd_y + '_4x4.nc'
  irrfn_list = file_search(irrfn_search_str)
  irrfn = irrfn_list[0]

  if strlen(file_basename(irrfn)) gt 1 then begin 
    print, irrfn
    irrfn_exist = 1
  endif else begin
    irrfn_exist = 0
  endelse

  ;cldfn_search_str = '/data/nier_ftp/CLOUD/V03/' + yyyy + mm + '/' + dd $
  cldfn_search_str = '/data/private/soodal/ln/GEMS/L2CLD/GK2_GEMS_L2_CLD_' + datetime_str + '_4x4.nc'
  cldfn_list = file_search(cldfn_search_str)
  cldfn = cldfn_list[0]

  if strlen(file_basename(cldfn)) gt 1 then begin 
    print, cldfn
    cldfn_exist = 1
  endif else begin
    cldfn_exist = 0
  endelse

  if gems_exist and l1crad_exist and irrfn_exist and cldfn_exist then begin 
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
    ;sza = dblarr(simrad_sz[0], simrad_sz[1])
    ;for iw = 0, 203 do begin
        ;sza[*, *] = synthetic.solarzenithangle
    ;endfor

    sza = synthetic.solarzenithangle
    vza = synthetic.viewingzenithangle


    gemsl1crad = ds_read_gems_l1c_4x4(l1cradfn)
    gemsirr = ds_read_gems_l1c_4x4(irrfn)
    gemscld = ds_read_gems_l2_o3p(cldfn)

    startidxarr = lonarr(simrad_sz[1])

    ;for ix = 0, simrad_sz[0]-1 do begin
        for iy = 0, simrad_sz[1]-1 do begin
            if syntheticwav[iy, 0] lt -1e29 then begin
                startidxarr[iy] = 0
                continue
            endif

            idx = where(gemsl1crad.wavelength[iy, *] eq syntheticwav[iy, 0], /null)
            if n_elements(idx) eq 0 then begin
                stop
            endif
            startidxarr[iy] = idx
        endfor
    ;endfor

    ;l1cradwav = double(gemsl1crad.wavelength[*, 0+11+2:203+11+2])
    ;l1cradwav0 = double(gemsl1crad.wavelength[*, 0+11:203+11])
    l1cradwav = dblarr(simrad_sz[1], simrad_sz[2])
    for iy = 0, simrad_sz[1]-1 do begin
        l1cradwav[iy, *] = gemsl1crad.wavelength[iy, startidxarr[iy]:startidxarr[iy]+203]
    endfor
            


    ;rad = reform(gemsl1crad.image_pixel_values[*, *, 0+11+2:203+11+2])
    rad = dblarr(simrad_sz[0], simrad_sz[1], simrad_sz[2])
    for ix = 0, simrad_sz[0]-1 do begin
        for iy = 0, simrad_sz[1]-1 do begin
            rad[ix, iy, 0:simrad_sz[2]-1] = gemsl1crad.image_pixel_values[ix, iy, startidxarr[iy]:startidxarr[iy]+simrad_sz[2]-1]
        endfor
    endfor

  syntheticwav[where(syntheticwav le -1e29, /null)] = !values.f_nan

  temp = fltarr(simrad_sz[0], simrad_sz[1])
  for ix = 0, simrad_sz[0]-1 do begin
    temp[ix, *] = total(finite(syntheticwav), 2)
  endfor

  ;gems_rad_norm = total(rad, 3)/198.
  gems_rad_norm = total(rad, 3)/temp

    for iw=0, simrad_sz[2]-1 do begin 
      rad[0:simrad_sz[0]-1, 0:simrad_sz[1]-1, iw] = $
          rad[0:simrad_sz[0]-1, 0:simrad_sz[1]-1, iw] / gems_rad_norm
    endfor

    ;----------------------------------------------------------------------------
    ; get irr more wider wavelength range for interpolation
    ;----------------------------------------------------------------------------
    ;irr = reform(gemsirr.image_pixel_values[*, 0+11+2:203+11+2])


    irr0 = dblarr(simrad_sz[1], simrad_sz[2]+2)
    irr = dblarr(simrad_sz[1], simrad_sz[2])
    l1cirrwav0 = dblarr(simrad_sz[1], simrad_sz[2]+2)
    l1cirrwav = dblarr(simrad_sz[1], simrad_sz[2])

    ix = floor(simrad_sz[0]/2.)
    for iy = 0, simrad_sz[1]-1 do begin
        if startidxarr[iy] ne 0 then begin
            irr0[iy, *] = gemsirr.image_pixel_values[iy, startidxarr[iy]-1:startidxarr[iy]+simrad_sz[2]]
            l1cirrwav0[iy, 0:simrad_sz[2]+1] = $
                gemsirr.wavelength[iy, startidxarr[iy]-1:startidxarr[iy]+simrad_sz[2]]
        endif else begin 
            irr0[iy, *] = gemsirr.image_pixel_values[iy, 0:simrad_sz[2]+1]
            l1cirrwav0[iy, 0:simrad_sz[2]+1] = $
                gemsirr.wavelength[iy, 0:simrad_sz[2]+1]
        endelse
        ;l1cirrwav[iy, *] = interpol(irr0[iy, *], l1cirrwav0[iy, *], l1cradwav[iy, *])
        irr[iy, *] = interpol( $
            irr0[iy, 0:simrad_sz[2]+1], l1cirrwav0[iy, 0:simrad_sz[2]+1], l1cradwav[iy, *])
    endfor

    ;gems_irrad_norm = total(irr, 2)/198.
    gems_irrad_norm = total(irr, 2)/total(finite(syntheticwav), 2)

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
    ecf02idx = where(ecf ge 0 and ecf le 0.2, /null)
    ecf02nanidx = where(ecf lt 0 or ecf gt 0.2, /null)
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
    outncfn_mean = '/data/private/soodal/softcal_test/residual/gems_merra2_' + datetime_str + '_mean.nc'
    outncfn_median = '/data/private/soodal/softcal_test/residual/gems_merra2_' + datetime_str + '_median.nc'
    outncfn_raw = '/data/private/soodal/softcal_test/residual/gems_merra2_' + datetime_str + '_raw.nc'

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
    ;ncdf_varput, id, _fitspec_difference_from_merra2, diff_mean
    ;ncdf_varput, id, _wavelength_data, l1cradwav
    ;ncdf_varput, id, _sza_data, sza
    ;ncdf_close, id

    ;===============================================================================
    ; write start median
    ;===============================================================================
    ;id = ncdf_create(outncfn_median)
    ;_ximage = ncdf_dimdef(id, 'image', diff_sz[0])
    ;_yspatial = ncdf_dimdef(id, 'spatial', diff_sz[1])
    ;_wavelength = ncdf_dimdef(id, 'wavelength', diff_sz[2])

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
    ;ncdf_varput, id, _wavelength_data, l1cradwav
    ;ncdf_varput, id, _sza_data, sza
    ;ncdf_close, id


    ;===============================================================================
    ; write start raw
    ;===============================================================================
    id = ncdf_create(outncfn_raw)
    _ximage = ncdf_dimdef(id, 'image', diff_sz[0])
    _yspatial = ncdf_dimdef(id, 'spatial', diff_sz[1])
    _wavelength = ncdf_dimdef(id, 'wavelength', diff_sz[2])

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

    _vza_data = ncdf_vardef(id, 'viewingzenithangle', $
        [_ximage, _yspatial], /double)

    ncdf_attput, id, _vza_data, 'units', 'degrees'
    ncdf_attput, id, _vza_data, 'long name', 'viewingzenithangle'

    _ecf_data = ncdf_vardef(id, 'effectivecloudfraction', $
        [_ximage, _yspatial], /double)

    ncdf_attput, id, _ecf_data, 'units', 'ratio'
    ncdf_attput, id, _ecf_data, 'long name', 'effectivecloudfraction'

    ; put file in data mode
    ncdf_control, id, /endef

    ; input data
    ncdf_varput, id, _fitspec_difference_from_merra2, diff_raw
    ncdf_varput, id, _wavelength_data, l1cradwav
    ncdf_varput, id, _sza_data, sza
    ncdf_varput, id, _vza_data, vza
    ncdf_varput, id, _ecf_data, ecf
    ncdf_close, id
  endif
endfor


end
