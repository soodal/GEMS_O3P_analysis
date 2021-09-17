pro collocate_merra2_chm_inst3_3d_on_gemsl2o3p, year, month, day, hour, minute $
  , method=method


; read GEMS
; read MERRA2
; collocate merra2 on gems
; the time interpolation for the merra2 to the GEMS time
; save the collocated merra2 to the nc file

if not keyword_set(method) then begin
  method = 'nearest'
endif

limit=[-10, 80, 60, 160]

yyyy = string(year, format='(i04)')
mm = string(month, format='(i02)')
dd = string(day, format='(i02)')
hh = string(hour, format='(i02)')
mi = string(minute, format='(i02)')

savpath = './collocate_merra2_on_gemsl2o3p/'
if not file_test(savpath) then begin
  file_mkdir, savpath
endif

;savfile = savpath + 'collocated_merra2_on_gems_'+yyyy+mm+dd+'_'+hh+mi+'.sav'

;savefile =  file_test(savfile) 
savefile = 0

jday = julday(month, day, year, 0, 0)
nday = jday+1
caldat, nday, month_n, day_n, year_n

; set read file name

; GEMS filename
;gemsfn = '../GEMS_O3P_Yonsei/out/310340_EOSRL/GK2B_GEMS_L2_20200616_0345_O3P_ND_DPRO_first.4x4_2008301839.nc4'
gemsfp = '/data2/L2_GEMS/val_1008/'
gemsfn = 'GK2_GEMS_O3P_' $
  + yyyy + mm + dd + '_' + hh + mi + '_v1.ba.nc4'

; MERRA2 filename
;merra2asmfn = '/data/MERRA2/' + yyyy + '/' + mm + '/MERRA2_400.tavg3_3d_asm_Nv.' $
  ;+ yyyy + mm + dd + '.nc4'
merra2chmfn = '/data/MERRA2/' + yyyy + '/' + mm + '/MERRA2_400.inst3_3d_chm_Nv.' $
  + yyyy + mm + dd + '.nc4'

; read file

if file_test(gemsfp+gemsfn) and file_test(merra2chmfn) then begin
  gemsvars = ds_read_gems_l2_o3p(gemsfp+gemsfn)

  ;merra2asm = ds_read_merra2_tavg3_3d_asm_nv(merra2asmfn)
  merra2chm = ds_read_merra2_inst3_3d_chm_nv(merra2chmfn)

  IF hour ge 20 and hour le 27 then BEGIN ; for 2345 UTC
    merra2asmfn_n = '/data/MERRA2/2020/08/MERRA2_400.tavg3_3d_asm_Nv.' $
      + yyyy_n + mm_n + dd_n + '.nc4'
    merra2chmfn_n = '/data/MERRA2/2020/08/MERRA2_400.inst3_3d_chm_Nv.' $
      + yyyy_n + mm_n + dd_n + '.nc4'
    merra2asm_n = ds_read_merra2_tavg3_3d_asm_nv(merra2asmfn_n)
    merra2chm_n = ds_read_merra2_inst3_3d_chm_nv(merra2chmfn_n)
  ENDIF

  ; collocation

  gems_sz = size(gemsvars.longitude, /dim)

  nearest_indices = lonarr(gems_sz[0], gems_sz[1], 2)
  nearest_indices[*] = -999

  ; lon 70-160 => 400-544, lat -15-55 => 150-290
  merra2chm_o3_roi = merra2chm.o3[400:544, 160:290, *, *]
  merra2chm_airdens_roi = merra2chm.airdens[400:544, 160:290, *, *]
  merra2chm_co_roi = merra2chm.co[400:544, 160:290, *, *]
  merra2chm_delp_roi = merra2chm.delp[400:544, 160:290, *, *]
  merra2chm_ps_roi = merra2chm.ps[400:544, 160:290, *, *]

  merra2chm_lon = merra2chm.lon # Replicate(1, n_elements(merra2chm.lat))
  merra2chm_lon_roi = merra2chm_lon[400:544, 160:290]
  merra2chm_lat = replicate(1, n_elements(merra2chm.lon)) # merra2chm.lat
  merra2chm_lat_roi = merra2chm_lat[400:544, 160:290]

  ;gems_roi_idx = where(merra2chm_lon ge 70 and merra2chm_lon lt 160 $
    ;and merra2chm_lat ge -15 and merra2chm_lat lt 55 )

  ;gems_roi_indices = array_indices(merra2chm_lon, gems_roi_idx)

  if method eq 'nearest' then BEGIN
    for ilon = 0, gems_sz[0]-1 do begin
      for ilat = 0, gems_sz[1]-1 do begin
        ; merra2 nearest index from the gems pixel
        merra2_nearest_index = search_closest_pixel( $
          merra2chm_lon_roi, merra2chm_lat_roi, $
          gemsvars.longitude[ilon, ilat], gemsvars.latitude[ilon, ilat], $
          maxlimit=0.5)
        if merra2_nearest_index eq -999 then BEGIN
          _nearest_indices = [-999, -999]
        endif else begin
          _nearest_indices = array_indices(merra2chm_lon_roi, merra2_nearest_index) 
        endelse

        ; merra2 nearest indices from the gems pixel
        nearest_indices[ilon, ilat, *] = _nearest_indices
      ENDFOR
    ENDFOR
  ;else if method eq 'bilinear' then BEGIN
  ENDIF

  ; MERRA2 chm time interpolation to GEMS time
  gemstime_jday = ds_gemsl2o3p_time2julday(gemsvars.time)
  caldat, gemstime_jday, gems_month, gems_year, gems_day, $
    gems_hour, gems_minute

  merra2chm_jday = julday(month, day, year, 0, 0) + merra2chm.time/60./24.
  merra2chm_jday_extended = merra2chm_jday


  ;  GEMS 2020-08-09T23:45 UTC

  time_frac = dblarr(gems_sz[0], 2)
  time_frac[*] = !values.d_nan

  merra2timestep = intarr(gems_sz[0])
  merra2timestep[*] = !values.d_nan

  for ix=0, gems_sz[0] - 1 do begin
    if finite(gemstime_jday[ix]) eq 1 then begin
      if gemstime_jday[ix] gt merra2chm_jday[-1] then BEGIN
        merra2chm_n_jday = julday(month_n, day_n, year_n, 0, 0) + merra2chm_n.time/60./24.
        merra2chm_jday_extended = [merra2chm_jday_extended, merra2chm_n_jday]
      ENDIF

      for imerra2time = 0, n_elements(merra2chm_jday_extended) - 2 do begin
        if gemstime_jday[ix] ge merra2chm_jday_extended[imerra2time] $
            and gemstime_jday[ix] lt merra2chm_jday_extended[imerra2time+1] then BEGIN
          aft_step = gemstime_jday[ix] - merra2chm_jday_extended[imerra2time]
          bef_step = merra2chm_jday_extended[imerra2time+1] - gemstime_jday[ix]
          after_frac = aft_step / (merra2chm_jday_extended[imerra2time+1] - merra2chm_jday_extended[imerra2time])
          before_frac = bef_step / (merra2chm_jday_extended[imerra2time+1] - merra2chm_jday_extended[imerra2time])

          time_frac[ix, 0] = after_frac
          time_frac[ix, 1] = before_frac

          merra2timestep[ix] = imerra2time 
          ;print, 'merra2timestep[ix]', ix, imerra2time
          break
        ENDIF
          
      endfor
    endif
  endfor

  merra2chm_o3_on_gems = dblarr(gems_sz[0], gems_sz[1], 72)
  merra2chm_o3_on_gems[*] = !values.d_nan

  merra2chm_airdens_on_gems = dblarr(gems_sz[0], gems_sz[1], 72)
  merra2chm_airdens_on_gems[*] = !values.d_nan
  
  merra2chm_co_on_gems = dblarr(gems_sz[0], gems_sz[1], 72)
  merra2chm_co_on_gems[*] = !values.d_nan

  merra2chm_delp_on_gems = dblarr(gems_sz[0], gems_sz[1], 72)
  merra2chm_delp_on_gems[*] = !values.d_nan

  merra2chm_ps_on_gems = dblarr(gems_sz[0], gems_sz[1])
  merra2chm_ps_on_gems[*] = !values.d_nan


  FOR ix=0, gems_sz[0] - 1 DO BEGIN
    FOR iy=0, gems_sz[1] -1 DO BEGIN
      if nearest_indices[ix, iy, 0] ne -999 then begin
        merra2chm_o3_on_gems[ix, iy, *] = $
          reform(time_frac[ix, 0] * merra2chm_o3_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]]) $
          + reform(time_frac[ix, 1] * merra2chm_o3_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]+1])

        merra2chm_airdens_on_gems[ix, iy, *] = $
          reform(time_frac[ix, 0] * merra2chm_airdens_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]]) $
          + reform(time_frac[ix, 1] * merra2chm_airdens_roi[$
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]+1])

        merra2chm_co_on_gems[ix, iy, *] = $
          reform(time_frac[ix, 0] * merra2chm_co_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]]) $
          + reform(time_frac[ix, 1] * merra2chm_co_roi[$
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]+1])

        merra2chm_delp_on_gems[ix, iy, *] = $
          reform(time_frac[ix, 0] * merra2chm_delp_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]]) $
          + reform(time_frac[ix, 1] * merra2chm_delp_roi[$
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]+1])

        merra2chm_ps_on_gems[ix, iy] = $
          reform(time_frac[ix, 0] * merra2chm_ps_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], $
            merra2timestep[ix]]) $
          + reform(time_frac[ix, 1] * merra2chm_ps_roi[$
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], $
            merra2timestep[ix]+1])
      endif else begin
        merra2chm_o3_on_gems[ix, iy, *] = !values.d_nan
        merra2chm_airdens_on_gems[ix, iy, *] = !values.d_nan
        merra2chm_co_on_gems[ix, iy, *] = !values.d_nan
        merra2chm_delp_on_gems[ix, iy, *] = !values.d_nan
        merra2chm_ps_on_gems[ix, iy] = !values.d_nan
      endelse

    ENDFOR
  ENDFOR


  ; save netcdf
  outncfp = '/data/MERRA2_collocated_on_gems/' 
  outncfn = gemsfn + '.merra2.nc4'
  if file_test(outncfp+outncfn) then begin
    file_delete, outncfp+outncfn
  endif

  id = ncdf_create(outncfp + outncfn)
  _ximage = ncdf_dimdef(id, 'image', gems_sz[0])
  _yspatial = ncdf_dimdef(id, 'spatial', gems_sz[1])
  _zheight = ncdf_dimdef(id, 'lev', 72)

  ; define variables
  o3mmr = ncdf_vardef(id, 'O3 Mass mixing ratio', [_ximage, _yspatial, _zheight], /double)
  lev = ncdf_vardef(id, 'lev', [_zheight])
  gems_lon = ncdf_vardef(id, 'longitude', [_ximage, _yspatial])
  gems_lat = ncdf_vardef(id, 'latitude', [_ximage, _yspatial])
  _airdens = ncdf_vardef(id, 'airdensity', [_ximage, _yspatial, _zheight])
  _co = ncdf_vardef(id, 'co', [_ximage, _yspatial, _zheight])
  _delp = ncdf_vardef(id, 'delp', [_ximage, _yspatial, _zheight])
  _ps = ncdf_vardef(id, 'ps', [_ximage, _yspatial])


  ncdf_attput, id, o3mmr, 'units', 'kg kg-1'
  ncdf_attput, id, lev, 'units', 'layer'
  ncdf_attput, id, gems_lon, 'units', 'degrees'
  ncdf_attput, id, gems_lat, 'units', 'degrees'
  ncdf_attput, id, _airdens, 'units', 'kg m-3'
  ncdf_attput, id, _co, 'units', 'mol mol-1'
  ncdf_attput, id, _delp, 'units', 'Pa'
  ncdf_attput, id, _ps, 'units', 'Pa'

  ; put file in data mode
  ncdf_control, id, /endef

  ; input data
  ncdf_varput, id, o3mmr, merra2chm_o3_on_gems
  ncdf_varput, id, lev, merra2chm.lev
  ncdf_varput, id, gems_lon, gemsvars.longitude
  ncdf_varput, id, gems_lat, gemsvars.latitude
  ncdf_varput, id, _airdens, merra2chm_airdens_on_gems
  ncdf_varput, id, _co, merra2chm_co_on_gems
  ncdf_varput, id, _delp, merra2chm_delp_on_gems
  ncdf_varput, id, _ps, merra2chm_ps_on_gems
  ncdf_close, id

  return
endif
end
