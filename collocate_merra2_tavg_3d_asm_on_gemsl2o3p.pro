pro collocate_merra2_tavg3_3d_asm_on_gemsl2o3p, year, month, day, hour, minute $
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

jday = julday(month, day, year, 0, 0)
nday = jday+1
caldat, nday, month_n, day_n, year_n

; set read file name

; GEMS filename
gemsfp = '/data2/L2_GEMS/val_1008/'
gemsfn = 'GK2_GEMS_O3P_' $
  + yyyy + mm + dd + '_' + hh + mi + '_v1.ba.nc4'

; MERRA2 filename
merra2_tavg3_3d_asm_fn = '/data/MERRA2/' + yyyy + '/' + mm + '/MERRA2_400.tavg3_3d_asm_Nv.' $
  + yyyy + mm + dd + '.nc4'
merra2fn = '/data/MERRA2/' + yyyy + '/' + mm + '/MERRA2_400.inst3_3d_chm_Nv.' $
  + yyyy + mm + dd + '.nc4'

; read file

if file_test(gemsfp+gemsfn) and file_test(merra2fn) then begin
  gemsvars = ds_read_gems_l2_o3p(gemsfp+gemsfn)

  merra2_tavg3_3d_asm = ds_read_merra2_tavg3_3d_asm_nv(merra2_tavg3_3d_asm_fn)
  ;merra2 = ds_read_merra2_inst3_3d_chm_nv(merra2fn)

  IF hour ge 20 and hour le 27 then BEGIN ; for 2345 UTC
    merra2_tavg3_3d_asm_fn_n = '/data/MERRA2/2020/08/MERRA2_400.tavg3_3d_asm_Nv.' $
      + yyyy_n + mm_n + dd_n + '.nc4'
    ;merra2fn_n = '/data/MERRA2/2020/08/MERRA2_400.inst3_3d_chm_Nv.' $
      ;+ yyyy_n + mm_n + dd_n + '.nc4'
    merra2_tavg3_3d_asm_n = ds_read_merra2_tavg3_3d_asm_nv(merra2_tavg3_3d_asm_fn_n)
    ;merra2_n = ds_read_merra2_inst3_3d_chm_nv(merra2fn_n)
  ENDIF

  ; collocation

  gems_sz = size(gemsvars.longitude, /dim)

  nearest_indices = lonarr(gems_sz[0], gems_sz[1], 2)
  nearest_indices[*] = -999

  ; lon 70-160 => 400-544, lat -15-55 => 150-290
  merra2_cloud_roi = merra2_tavg3_3d_asm.cloud[400:544, 160:290, *, *]

  merra2_delp_roi = merra2_tavg3_3d_asm.delp[400:544, 160:290, *, *]

  merra2_epv_roi = merra2_tavg3_3d_asm.epv[400:544, 160:290, *, *]

  merra2_h_roi = merra2_tavg3_3d_asm.h[400:544, 160:290, *, *]

  merra2_lat = replicate(1, n_elements(merra2_tavg3_3d_asm.lon)) # merra2_tavg3_3d_asm.lat
  merra2_lat_roi = merra2_lat[400:544, 160:290]

  merra2_lon = merra2_tavg3_3d_asm.lon # Replicate(1, n_elements(merra2_tavg3_3d_asm.lat))
  merra2_lon_roi = merra2_lon[400:544, 160:290]

  merra2_o3_roi = merra2_tavg3_3d_asm.o3[400:544, 160:290, *, *]

  merra2_omega_roi = merra2_tavg3_3d_asm.omega[400:544, 160:290, *, *]

  merra2_phis_roi = merra2_tavg3_3d_asm.phis[400:544, 160:290, *] ; surface geopotential height

  merra2_pl_roi = merra2_tavg3_3d_asm.pl[400:544, 160:290, *, *] ; mid level pressure

  merra2_ps_roi = merra2_tavg3_3d_asm.ps[400:544, 160:290, *] ; surface pressure

  merra2_qi_roi = merra2_tavg3_3d_asm.qi[400:544, 160:290, *, *] ; mass fraction of cloud ice water

  merra2_ql_roi = merra2_tavg3_3d_asm.ql[400:544, 160:290, *, *] ; mass fraction of cloud liquid water

  merra2_qv_roi = merra2_tavg3_3d_asm.qv[400:544, 160:290, *, *] ; specific humidity

  merra2_rh_roi = merra2_tavg3_3d_asm.rh[400:544, 160:290, *, *] ; relative humidity after moist

  merra2_slp_roi = merra2_tavg3_3d_asm.slp[400:544, 160:290, *]

  merra2_t_roi = merra2_tavg3_3d_asm.t[400:544, 160:290, *, *]

  merra2_u_roi = merra2_tavg3_3d_asm.u[400:544, 160:290, *, *] ; eastward wind
  merra2_v_roi = merra2_tavg3_3d_asm.v[400:544, 160:290, *, *] ; northward wind

  if method eq 'nearest' then BEGIN
    for ilon = 0, gems_sz[0]-1 do begin
      for ilat = 0, gems_sz[1]-1 do begin
        ; merra2 nearest index from the gems pixel
        merra2_nearest_index = search_closest_pixel( $
          merra2_lon_roi, merra2_lat_roi, $
          gemsvars.longitude[ilon, ilat], gemsvars.latitude[ilon, ilat], $
          maxlimit=0.5)
        if merra2_nearest_index eq -999 then BEGIN
          _nearest_indices = [-999, -999]
        endif else begin
          _nearest_indices = array_indices(merra2_lon_roi, merra2_nearest_index) 
        endelse

        ; merra2 nearest indices from the gems pixel
        nearest_indices[ilon, ilat, *] = _nearest_indices
      ENDFOR
    ENDFOR
  ;else if method eq 'bilinear' then BEGIN
  ENDIF

  ; MERRA2 asm time interpolation to GEMS time
  gemstime_jday = ds_gemsl2o3p_time2julday(gemsvars.time)
  caldat, gemstime_jday, gems_month, gems_year, gems_day, $
    gems_hour, gems_minute

  merra2asm_jday = julday(month, day, year, 0, 0) + merra2_tavg3_3d_asm.time/60./24.
  merra2asm_jday_extended = merra2asm_jday


  time_frac = dblarr(gems_sz[0], 2)
  time_frac[*] = !values.d_nan

  merra2timestep = intarr(gems_sz[0])
  merra2timestep[*] = !values.d_nan

  for ix=0, gems_sz[0] - 1 do begin
    if finite(gemstime_jday[ix]) eq 1 then begin
      if gemstime_jday[ix] gt merra2asm_jday[-1] then BEGIN
        merra2_n_jday = julday(month_n, day_n, year_n, 0, 0) + merra2_tavg3_3d_asm_n.time/60./24.
        merra2asm_jday_extended = [merra2asm_jday_extended, merra2_n_jday]
      ENDIF

      for imerra2time = 0, n_elements(merra2asm_jday_extended) - 2 do begin
        if gemstime_jday[ix] ge merra2asm_jday_extended[imerra2time] $
            and gemstime_jday[ix] lt merra2asm_jday_extended[imerra2time+1] then BEGIN
          aft_step = gemstime_jday[ix] - merra2asm_jday_extended[imerra2time]
          bef_step = merra2asm_jday_extended[imerra2time+1] - gemstime_jday[ix]
          after_frac = aft_step / (merra2asm_jday_extended[imerra2time+1] - merra2asm_jday_extended[imerra2time])
          before_frac = bef_step / (merra2asm_jday_extended[imerra2time+1] - merra2asm_jday_extended[imerra2time])

          time_frac[ix, 0] = after_frac
          time_frac[ix, 1] = before_frac

          merra2timestep[ix] = imerra2time 
          break
        ENDIF
          
      endfor
    endif
  endfor
  
  merra2_cloud_on_gems = dblarr(gems_sz[0], gems_sz[1], 72)
  merra2_cloud_on_gems[*] = !values.d_nan

  merra2_delp_on_gems = dblarr(gems_sz[0], gems_sz[1], 72)
  merra2_delp_on_gems[*] = !values.d_nan

  merra2_epv_on_gems = dblarr(gems_sz[0], gems_sz[1], 72)
  merra2_epv_on_gems[*] = !values.d_nan

  merra2_h_on_gems = dblarr(gems_sz[0], gems_sz[1], 72)
  merra2_h_on_gems[*] = !values.d_nan

  merra2_o3_on_gems = dblarr(gems_sz[0], gems_sz[1], 72)
  merra2_o3_on_gems[*] = !values.d_nan

  merra2_omega_on_gems = dblarr(gems_sz[0], gems_sz[1], 72)
  merra2_omega_on_gems[*] = !values.d_nan

  merra2_phis_on_gems = dblarr(gems_sz[0], gems_sz[1])
  merra2_phis_on_gems[*] = !values.d_nan

  merra2_pl_on_gems = dblarr(gems_sz[0], gems_sz[1], 72)
  merra2_pl_on_gems[*] = !values.d_nan

  merra2_ps_on_gems = dblarr(gems_sz[0], gems_sz[1])
  merra2_ps_on_gems[*] = !values.d_nan

  merra2_qi_on_gems = dblarr(gems_sz[0], gems_sz[1], 72)
  merra2_qi_on_gems[*] = !values.d_nan

  merra2_ql_on_gems = dblarr(gems_sz[0], gems_sz[1], 72)
  merra2_ql_on_gems[*] = !values.d_nan

  merra2_qv_on_gems = dblarr(gems_sz[0], gems_sz[1], 72)
  merra2_qv_on_gems[*] = !values.d_nan

  merra2_rh_on_gems = dblarr(gems_sz[0], gems_sz[1], 72)
  merra2_rh_on_gems[*] = !values.d_nan

  merra2_slp_on_gems = dblarr(gems_sz[0], gems_sz[1])
  merra2_slp_on_gems[*] = !values.d_nan

  merra2_t_on_gems = dblarr(gems_sz[0], gems_sz[1], 72)
  merra2_t_on_gems[*] = !values.d_nan

  merra2_u_on_gems = dblarr(gems_sz[0], gems_sz[1], 72)
  merra2_u_on_gems[*] = !values.d_nan

  merra2_v_on_gems = dblarr(gems_sz[0], gems_sz[1], 72)
  merra2_v_on_gems[*] = !values.d_nan


  FOR ix=0, gems_sz[0] - 1 DO BEGIN
    FOR iy=0, gems_sz[1] -1 DO BEGIN
      if nearest_indices[ix, iy, 0] ne -999 then begin
        merra2_cloud_on_gems[ix, iy, *] = $
          reform(time_frac[ix, 0] * merra2_cloud_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]]) $
          + reform(time_frac[ix, 1] * merra2_cloud_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]+1])
        merra2_delp_on_gems[ix, iy, *] = $
          reform(time_frac[ix, 0] * merra2_delp_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]]) $
          + reform(time_frac[ix, 1] * merra2_delp_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]+1])
        merra2_epv_on_gems[ix, iy, *] = $
          reform(time_frac[ix, 0] * merra2_epv_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]]) $
          + reform(time_frac[ix, 1] * merra2_epv_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]+1])
        merra2_h_on_gems[ix, iy, *] = $
          reform(time_frac[ix, 0] * merra2_h_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]]) $
          + reform(time_frac[ix, 1] * merra2_h_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]+1])
        merra2_o3_on_gems[ix, iy, *] = $
          reform(time_frac[ix, 0] * merra2_o3_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]]) $
          + reform(time_frac[ix, 1] * merra2_o3_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]+1])
        merra2_omega_on_gems[ix, iy, *] = $
          reform(time_frac[ix, 0] * merra2_omega_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]]) $
          + reform(time_frac[ix, 1] * merra2_omega_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]+1])
        merra2_phis_on_gems[ix, iy] = $
          reform(time_frac[ix, 0] * merra2_phis_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], $
            merra2timestep[ix]]) $
          + reform(time_frac[ix, 1] * merra2_phis_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], $
            merra2timestep[ix]+1])
        merra2_pl_on_gems[ix, iy, *] = $
          reform(time_frac[ix, 0] * merra2_pl_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]]) $
          + reform(time_frac[ix, 1] * merra2_pl_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]+1])
        merra2_ps_on_gems[ix, iy] = $
          reform(time_frac[ix, 0] * merra2_ps_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], $
            merra2timestep[ix]]) $
          + reform(time_frac[ix, 1] * merra2_ps_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], $
            merra2timestep[ix]+1])
        merra2_qi_on_gems[ix, iy, *] = $
          reform(time_frac[ix, 0] * merra2_qi_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]]) $
          + reform(time_frac[ix, 1] * merra2_qi_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]+1])
        merra2_ql_on_gems[ix, iy, *] = $
          reform(time_frac[ix, 0] * merra2_ql_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]]) $
          + reform(time_frac[ix, 1] * merra2_ql_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]+1])
        merra2_qv_on_gems[ix, iy, *] = $
          reform(time_frac[ix, 0] * merra2_qv_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]]) $
          + reform(time_frac[ix, 1] * merra2_qv_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]+1])
        merra2_rh_on_gems[ix, iy, *] = $
          reform(time_frac[ix, 0] * merra2_rh_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]]) $
          + reform(time_frac[ix, 1] * merra2_rh_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]+1])
        merra2_slp_on_gems[ix, iy] = $
          reform(time_frac[ix, 0] * merra2_slp_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], $
            merra2timestep[ix]]) $
          + reform(time_frac[ix, 1] * merra2_slp_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], $
            merra2timestep[ix]+1])
        merra2_t_on_gems[ix, iy, *] = $
          reform(time_frac[ix, 0] * merra2_t_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]]) $
          + reform(time_frac[ix, 1] * merra2_t_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]+1])
        merra2_u_on_gems[ix, iy, *] = $
          reform(time_frac[ix, 0] * merra2_u_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]]) $
          + reform(time_frac[ix, 1] * merra2_u_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]+1])
        merra2_v_on_gems[ix, iy, *] = $
          reform(time_frac[ix, 0] * merra2_v_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]]) $
          + reform(time_frac[ix, 1] * merra2_v_roi[ $
            nearest_indices[ix, iy, 0], nearest_indices[ix, iy, 1], *, $
            merra2timestep[ix]+1])
      endif else begin
        merra2_cloud_on_gems[ix, iy, *] = !values.d_nan
        merra2_delp_on_gems[ix, iy, *] = !values.d_nan
        merra2_epv_on_gems[ix, iy, *] = !values.d_nan
        merra2_h_on_gems[ix, iy, *] = !values.d_nan
        merra2_o3_on_gems[ix, iy, *] = !values.d_nan
        merra2_omega_on_gems[ix, iy, *] = !values.d_nan
        merra2_phis_on_gems[ix, iy] = !values.d_nan
        merra2_pl_on_gems[ix, iy, *] = !values.d_nan
        merra2_ps_on_gems[ix, iy] = !values.d_nan
        merra2_qi_on_gems[ix, iy, *] = !values.d_nan
        merra2_ql_on_gems[ix, iy, *] = !values.d_nan
        merra2_qv_on_gems[ix, iy, *] = !values.d_nan
        merra2_rh_on_gems[ix, iy, *] = !values.d_nan
        merra2_slp_on_gems[ix, iy] = !values.d_nan
        merra2_t_on_gems[ix, iy, *] = !values.d_nan
        merra2_u_on_gems[ix, iy, *] = !values.d_nan
        merra2_v_on_gems[ix, iy, *] = !values.d_nan
      endelse

    ENDFOR
  ENDFOR

  ; save netcdf
  outncfp = '/data/MERRA2_collocated_on_gems/' 
  outncfn = gemsfn + '.merra2_tavg3_3d_asm.nc4'
  if file_test(outncfp+outncfn) then begin
    file_delete, outncfp+outncfn
  endif

  id = ncdf_create(outncfp + outncfn)
  _ximage = ncdf_dimdef(id, 'image', gems_sz[0])
  _yspatial = ncdf_dimdef(id, 'spatial', gems_sz[1])
  _zheight = ncdf_dimdef(id, 'lev', 72)

  ; define variables
  _cloud = ncdf_vardef(id, 'CLOUD', [_ximage, _yspatial, _zheight], /double)
  _delp = ncdf_vardef(id, 'DELP', [_ximage, _yspatial, _zheight], /double)
  _epv = ncdf_vardef(id, 'EPV', [_ximage, _yspatial, _zheight], /double)
  _h = ncdf_vardef(id, 'H', [_ximage, _yspatial, _zheight], /double)
  _lev = ncdf_vardef(id, 'LEV', [_zheight], /double)
  _o3 = ncdf_vardef(id, 'O3', [_ximage, _yspatial, _zheight], /double)
  _omega = ncdf_vardef(id, 'OMEGA', [_ximage, _yspatial, _zheight], /double)
  _phis = ncdf_vardef(id, 'PHIS', [_ximage, _yspatial], /double)
  _pl = ncdf_vardef(id, 'PL', [_ximage, _yspatial, _zheight], /double)
  _ps = ncdf_vardef(id, 'PS', [_ximage, _yspatial], /double)
  _qi = ncdf_vardef(id, 'QI', [_ximage, _yspatial, _zheight], /double)
  _ql = ncdf_vardef(id, 'QL', [_ximage, _yspatial, _zheight], /double)
  _qv = ncdf_vardef(id, 'QV', [_ximage, _yspatial, _zheight], /double)
  _rh = ncdf_vardef(id, 'RH', [_ximage, _yspatial, _zheight], /double)
  _slp = ncdf_vardef(id, 'SLP', [_ximage, _yspatial], /double)
  _t = ncdf_vardef(id, 'T', [_ximage, _yspatial, _zheight], /double)
  _u = ncdf_vardef(id, 'U', [_ximage, _yspatial, _zheight], /double)
  _v = ncdf_vardef(id, 'V', [_ximage, _yspatial, _zheight], /double)

  _lon = ncdf_vardef(id, 'lon', [_ximage, _yspatial])
  _lat = ncdf_vardef(id, 'lat', [_ximage, _yspatial])

  ncdf_attput, id, _cloud, 'units', '1'
  ncdf_attput, id, _delp, 'units', 'Pa'
  ncdf_attput, id, _epv, 'units', 'K m^+2 kg^-1 s^-1'
  ncdf_attput, id, _h, 'units', 'm'
  ncdf_attput, id, _lev, 'units', 'layer'
  ncdf_attput, id, _o3, 'units', 'kg kg^-1'
  ncdf_attput, id, _omega, 'units', 'Pa s^-1kg kg-1'
  ncdf_attput, id, _phis, 'units', 'm^2 s^-2'
  ncdf_attput, id, _pl, 'units', 'Pa'
  ncdf_attput, id, _ps, 'units', 'Pa'
  ncdf_attput, id, _qi, 'units', 'kg kg^-1'
  ncdf_attput, id, _ql, 'units', 'kg kg^-1'
  ncdf_attput, id, _qv, 'units', 'kg kg^-1'
  ncdf_attput, id, _rh, 'units', '1'
  ncdf_attput, id, _slp, 'units', 'Pa'
  ncdf_attput, id, _t, 'units', 'K'
  ncdf_attput, id, _u, 'units', 'm s^-1'
  ncdf_attput, id, _v, 'units', 'm s^-1'

  ncdf_attput, id, _lon, 'units', 'degrees'
  ncdf_attput, id, _lat, 'units', 'degrees'

  ncdf_attput, id, _cloud, 'long name', 'cloud fraction for radiation'
  ncdf_attput, id, _delp, 'long name', 'pressure thickness'
  ncdf_attput, id, _epv, 'long name', 'ertels potential vorticity'
  ncdf_attput, id, _h, 'long name', 'mid layer heights'
  ncdf_attput, id, _lev, 'long name', 'vertical level'
  ncdf_attput, id, _o3, 'long name', 'ozone mass mixing ratio'
  ncdf_attput, id, _omega, 'long name', 'vertical pressure velocity'
  ncdf_attput, id, _phis, 'long name', 'surface geopotential height'
  ncdf_attput, id, _pl, 'long name', 'mid level pressure'
  ncdf_attput, id, _ps, 'long name', 'surface pressure'
  ncdf_attput, id, _qi, 'long name', 'mass fraction of cloud ice water'
  ncdf_attput, id, _ql, 'long name', 'mass fraction of cloud liquid water'
  ncdf_attput, id, _qv, 'long name', 'specific humidity'
  ncdf_attput, id, _rh, 'long name', 'relative humidity'
  ncdf_attput, id, _slp, 'long name', 'sea level pressure'
  ncdf_attput, id, _t, 'long name', 'air temperature'
  ncdf_attput, id, _u, 'long name', 'eastward wind'
  ncdf_attput, id, _v, 'long name', 'northward wind'

  ncdf_attput, id, _lon, 'long name', 'longitude'
  ncdf_attput, id, _lat, 'long name', 'latitude'

  ; put file in data mode
  ncdf_control, id, /endef

  ; input data
  ncdf_varput, id, _cloud, merra2_cloud_on_gems
  ncdf_varput, id, _delp, merra2_delp_on_gems
  ncdf_varput, id, _epv, merra2_epv_on_gems
  ncdf_varput, id, _h, merra2_h_on_gems
  ncdf_varput, id, _lev, merra2_tavg3_3d_asm.lev
  ncdf_varput, id, _o3, merra2_o3_on_gems
  ncdf_varput, id, _omega, merra2_omega_on_gems
  ncdf_varput, id, _phis, merra2_phis_on_gems
  ncdf_varput, id, _pl, merra2_pl_on_gems
  ncdf_varput, id, _ps, merra2_ps_on_gems
  ncdf_varput, id, _qi, merra2_qi_on_gems
  ncdf_varput, id, _ql, merra2_ql_on_gems
  ncdf_varput, id, _qv, merra2_qv_on_gems
  ncdf_varput, id, _rh, merra2_rh_on_gems
  ncdf_varput, id, _slp, merra2_slp_on_gems
  ncdf_varput, id, _t, merra2_t_on_gems
  ncdf_varput, id, _u, merra2_u_on_gems
  ncdf_varput, id, _v, merra2_v_on_gems

  ncdf_varput, id, _lon, gemsvars.longitude
  ncdf_varput, id, _lat, gemsvars.latitude
  ncdf_close, id

  return
endif
end