
nx_monthly = 91
ny_monthly = 61

monthly_lon = findgen(nx_monthly) + 70
monthly_lon = rebin(monthly_lon, nx_monthly, ny_monthly)

monthly_lat = findgen(ny_monthly) - 10
monthly_lat = transpose(rebin(monthly_lat, ny_monthly, nx_monthly))
saveflag = 0

for ih = 0, 23 do begin
  print, 'ih', ih
  ; E-W 70~160, 1x1 degree 
  ; S-N -10~50, 1x1 degree 
  monthly_tropo_o3 = fltarr(90, 60)
  monthly_tropo_o3[*] = !values.f_nan

  monthly_tropo_o3_num = fltarr(90, 60)
  monthly_tropo_o3_num[*] = 0

  monthly_tropo_o3_stddev = fltarr(90, 60)
  monthly_tropo_o3_stddev[*] = 0

  monthly_300hpa_o3 = fltarr(90, 60)
  monthly_300hpa_o3[*] = !values.f_nan

  monthly_300hpa_o3_num = fltarr(90, 60)
  monthly_300hpa_o3_num[*] = 0

  monthly_300hpa_o3_stddev = fltarr(90, 60)
  monthly_300hpa_o3_stddev[*] = 0


  for imonth = 6, 6 do begin
    jd_list = timegen(start=julday(imonth, 1, 2021, ih, 45), $
      final=julday(imonth+1, 1, 2021, ih, 45)-1, units='Days')
      ;final=julday(imonth+1, 1, 2020, ih, 45)-1, units='Days')

    ndays = n_elements(jd_list)

    o3_all = fltarr(174, 512, ndays)
    o3_all[*] = !values.f_nan

    o3_300hpa_all = fltarr(174, 512, ndays)
    o3_300hpa_all[*] = !values.f_nan

    lat_all = fltarr(174, 512, ndays)
    lat_all[*] = !values.f_nan
    lon_all = fltarr(174, 512, ndays)
    lon_all[*] = !values.f_nan

    for iday=0, ndays-1 do begin
      print, 'iday', iday
      caldat, jd_list[iday], month, day, year, hour, minute
      yyyy = string(year, format='(i04)')
      mm = string(month, format='(i02)')
      dd = string(day, format='(i02)')
      hh = string(hour, format='(i02)')
      mi = string(minute, format='(i02)')

      datetime_str = yyyy + mm + dd + '_' + hh + mi

      path = '/data/nier_ftp/O3P/V03/' + yyyy + mm + '/' + dd + '/'
      filepattern = 'GK2_GEMS_L2_' + datetime_str + '_O3P_*_*_BIN4x4.nc' 

      ;path = '/data2/L2_GEMS.nier/monthly_val/O3P/'
      ;filepattern = 'GK2_GEMS_L2_' + datetime_str + 'O3P_*_*_BIN4x4.nc' 

      ;path = '/data2/L2_GEMS/val_1008/'
      ;filepattern = 'GK2_GEMS_O3P_' + datetime_str + '.nc4' 

      filelist = file_search(path + filepattern)
      print, path + filepattern
      print, filelist
      if n_elements(filelist) eq 1 then begin
        if strlen(filelist) ne 0 then begin 
          saveflag = 1
        
          gemsl2o3p = ds_read_gems_l2_o3p(filelist[0])

          ds_gems_l2o3p_accum, gemsl2o3p, o3under300hpa, hpa=300.

          casz = size(gemsl2o3p.ColumnAmountO3, /dimension)
          if casz[0] eq 3 then begin
            tropo_o3 = reform(gemsl2o3p.ColumnAmountO3[2, *, *])
          endif else if casz[2] eq 3 then begin
            tropo_o3 = reform(gemsl2o3p.ColumnAmountO3[*, *, 3])
          endif

          sz = size(tropo_o3, /dimension)

          lat = gemsl2o3p.latitude
          lon = gemsl2o3p.longitude

          ecf = gemsl2o3p.EffectiveCloudFractionUV
          nanidx = where(ecf lt 0.2, /null)
          lat[nanidx] = !values.f_nan
          lon[nanidx] = !values.f_nan
          tropo_o3[nanidx] = !values.f_nan

          ;faf = gemsl2o3p.FinalAlgorithmFlags
          ;nanidx = where(faf ne 0, /null)
          ;lat[nanidx] = !values.f_nan
          ;lon[nanidx] = !values.f_nan
          ;tropo_o3[nanidx] = !values.f_nan

          pqf = gemsl2o3p.ProcessingQualityFlags
          nanidx = where(pqf ne 0, /null)
          lat[nanidx] = !values.f_nan
          lon[nanidx] = !values.f_nan
          tropo_o3[nanidx] = !values.f_nan

          o3_all[0:sz[0]-1, 0:sz[1]-1, iday] = tropo_o3
          lon_all[0:sz[0]-1, 0:sz[1]-1, iday] = lon
          lat_all[0:sz[0]-1, 0:sz[1]-1, iday] = lat

          o3_300hpa_all[0:sz[0]-1, 0:sz[1]-1, iday] = o3under300hpa
        endif
      endif
    endfor


    if saveflag then begin
      for ix=0, nx_monthly-2 do begin
        for iy=0, ny_monthly-2 do begin
          idx = where(lon_all GE monthly_lon[ix, iy] and lon_all LT monthly_lon[ix+1, iy] and $
            lat_all GE monthly_lat[ix, iy] and lat_all LT monthly_lat[ix, iy+1], ncount, /null)

          if ncount ne 0 then begin
            monthly_tropo_o3[ix, iy] = mean(o3_all[idx], /nan)
            monthly_tropo_o3_num[ix, iy] = ncount
            monthly_tropo_o3_stddev[ix, iy] = stddev(o3_all[idx], /nan)

            monthly_300hpa_o3[ix, iy] = mean(o3_300hpa_all[idx], /nan)
            monthly_300hpa_o3_num[ix, iy] = ncount
            monthly_300hpa_o3_stddev[ix, iy] = stddev(o3_300hpa_all[idx], /nan)

          endif else begin
            monthly_tropo_o3[ix, iy] = !values.f_nan
            monthly_tropo_o3_num[ix, iy] = 0
            monthly_tropo_o3_stddev[ix, iy] = !values.f_nan

            monthly_300hpa_o3[ix, iy] = !values.f_nan
            monthly_300hpa_o3_num[ix, iy] = 0
            monthly_300hpa_o3_stddev[ix, iy] = !values.f_nan
          endelse

        endfor
      endfor

      
      savepath = '~/data/monthly_tropo_o3/'

      save, filename = savepath + yyyy + mm + '_' + hh + 'utc.sav', $
        monthly_tropo_o3, monthly_tropo_o3_num, monthly_tropo_o3_stddev, monthly_lon, monthly_lat, $
        monthly_300hpa_o3, monthly_300hpa_o3_num, monthly_300hpa_o3_stddev
    endif
    saveflag = 0
  endfor
endfor
end

