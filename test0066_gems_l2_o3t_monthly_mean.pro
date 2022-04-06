
nx_monthly = 361
ny_monthly = 241

monthly_lon = findgen(nx_monthly)/4 + 70
monthly_lon = rebin(monthly_lon, nx_monthly, ny_monthly)

monthly_lat = findgen(ny_monthly)/4 - 10
monthly_lat = transpose(rebin(monthly_lat, ny_monthly, nx_monthly))
saveflag = 0

ny_max = 2048
nx_max = 696
for ih = 0, 23 do begin
  print, 'ih', ih
  ; E-W 70~160, 1x1 degree 
  ; S-N -10~50, 1x1 degree 
  monthly_o3 = fltarr(nx_monthly-1, ny_monthly-1)
  monthly_o3[*] = !values.f_nan

  monthly_o3_num = fltarr(nx_monthly-1, ny_monthly-1)
  monthly_o3_num[*] = 0

  monthly_o3_stddev = fltarr(nx_monthly-1, ny_monthly-1)
  monthly_o3_stddev[*] = 0

  for imonth = 10, 10 do begin
    jd_list = timegen(start=julday(imonth, 1, 2020, ih, 45), $
      final=julday(imonth+1, 1, 2020, ih, 45)-1, units='Days')
      ;final=julday(imonth+1, 1, 2020, ih, 45)-1, units='Days')

    ndays = n_elements(jd_list)

    o3_all = fltarr(nx_max, ny_max, ndays)
    o3_all[*] = !values.f_nan

    lat_all = fltarr(nx_max, ny_max, ndays)
    lat_all[*] = !values.f_nan

    lon_all = fltarr(nx_max, ny_max, ndays)
    lon_all[*] = !values.f_nan

    ecf_all = fltarr(nx_max, ny_max, ndays)
    ecf_all[*] = !values.f_nan

    for iday=0, ndays-1 do begin
      print, 'iday', iday
      caldat, jd_list[iday], month, day, year, hour, minute
      yyyy = string(year, format='(i04)')
      mm = string(month, format='(i02)')
      dd = string(day, format='(i02)')
      hh = string(hour, format='(i02)')
      mi = string(minute, format='(i02)')

      datetime_str = yyyy + mm + dd + '_' + hh + mi

      ;path = '/data/nier_ftp/O3T/V03/' + yyyy + mm + '/' + dd + '/'
      ;filepattern = 'GK2_GEMS_L2_' + datetime_str + '_O3T_*_*_BIN4x4.nc' 

      ;path = '/data2/L2_GEMS.nier/monthly_val/O3P/'
      ;filepattern = 'GK2_GEMS_L2_' + datetime_str + 'O3P_*_*_BIN4x4.nc' 

      path = '/data2/L2_GEMS/val_1008/'
      filepattern = 'GEMS_O3T_' + datetime_str + '.nc' 
      filelist = file_search(path + filepattern)

      cldpath = '/data2/L2_GEMS/val_1008/'
      cldfilepattern = 'GK2_GEMS_CLD_' + datetime_str + '_ver1030.nc' 
      cldfilelist = file_search(cldpath + cldfilepattern)
      ;print, path + filepattern
      print, filelist
      ;print, cldpath + cldfilepattern
      print, cldfilelist
      if n_elements(filelist) eq 1 and n_elements(cldfilelist) eq 1 then begin
        if strlen(filelist) ne 0 and strlen(filelist) ne 0 then begin 
          saveflag = 1
        
          gemsl2o3t = ds_read_gems_l2_o3t(filelist[0])

          cao3 = gemsl2o3t.ColumnAmountO3
          nanidx = where(cao3 lt -990, /null)
          cao3[nanidx] = !values.f_nan
          cao3sz = size(gemsl2o3t.ColumnAmountO3, /dimension)
          ;ecf = gemsl2o3t.EffectiveCloudFraction
          nanidx = where(cao3 lt -990, /null)

          lon = gemsl2o3t.Longitude
          lat = gemsl2o3t.latitude

          ;nanidx = where(ecf lt 0.2, /null)


          ;lat[nanidx] = !values.f_nan
          ;lon[nanidx] = !values.f_nan
          ;cao3[nanidx] = !values.f_nan

          ;nanidx = where(pqf ne 0, /null)
          ;lat[nanidx] = !values.f_nan
          ;lon[nanidx] = !values.f_nan
          ;cao3[nanidx] = !values.f_nan

          o3_all[0:cao3sz[0]-1, 0:cao3sz[1]-1, iday] = cao3
          lon_all[0:cao3sz[0]-1, 0:cao3sz[1]-1, iday] = lon
          lat_all[0:cao3sz[0]-1, 0:cao3sz[1]-1, iday] = lat

          gemsl2cld= ds_read_gems_l2_cld(cldfilelist[0])
          ecf = gemsl2cld.EffectiveCloudFraction


          ecf_all[0:cao3sz[0]-1, 0:cao3sz[1]-1, iday] = ecf

        endif
      endif
    endfor


    if saveflag then begin
      for ix=0, nx_monthly-2 do begin
        for iy=0, ny_monthly-2 do begin
          idx = where(lon_all GE monthly_lon[ix, iy] and lon_all LT monthly_lon[ix+1, iy] and $
            lat_all GE monthly_lat[ix, iy] and lat_all LT monthly_lat[ix, iy+1] and $ 
            ecf_all LT 0.2 , ncount, /null)

          if ncount ne 0 then begin
            monthly_o3[ix, iy] = mean(o3_all[idx], /nan)
            monthly_o3_num[ix, iy] = ncount
            monthly_o3_stddev[ix, iy] = stddev(o3_all[idx], /nan)
          endif else begin
            monthly_o3[ix, iy] = !values.f_nan
            monthly_o3_num[ix, iy] = 0
            monthly_o3_stddev[ix, iy] = !values.f_nan
          endelse

        endfor
      endfor

      savepath = '~/data/monthly_tropo_o3/'

      save, filename = savepath + 'gems_o3t_toz_' + yyyy + mm + '_' + hh + 'utc.sav', $
        monthly_o3, monthly_o3_num, monthly_o3_stddev, monthly_lon, monthly_lat
    endif
    saveflag = 0

    ;save, filename = savepath + 'gems_o3t_cao3_lon_lat_' + $
      ;yyyy + mm + '_' + hh + 'utc.sav', $
      ;o3_all, lon_all, lat_all
  endfor
endfor
end

