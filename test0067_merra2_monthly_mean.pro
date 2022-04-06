

;nx_monthly = 91
;ny_monthly = 61

;monthly_lon = findgen(nx_monthly) + 70
;monthly_lon = rebin(monthly_lon, nx_monthly, ny_monthly)

;monthly_lat = findgen(ny_monthly) - 10
;monthly_lat = transpose(rebin(monthly_lat, ny_monthly, nx_monthly))
saveflag = 0

for ih = 3, 3 do begin
  ;print, 'ih', ih
  ; E-W 70~160, 1x1 degree 
  ; S-N -10~50, 1x1 degree 

  monthly_300hpa_o3 = fltarr(90, 60)
  monthly_300hpa_o3[*] = !values.f_nan

  monthly_300hpa_o3_num = fltarr(90, 60)
  monthly_300hpa_o3_num[*] = 0

  monthly_300hpa_o3_stddev = fltarr(90, 60)
  monthly_300hpa_o3_stddev[*] = 0


  for imonth = 10, 10 do begin
    jd_list = timegen(start=julday(imonth, 1, 2020, ih, 45), $
      ;final=julday(imonth+1, 1, 2021, 23, 45)-1, units='Days')
      final=julday(imonth+1, 1, 2020, ih, 45)-1, units='Days')

    ndays = n_elements(jd_list)

    o3_all = fltarr(576, 361, ndays)
    o3_all[*] = !values.f_nan
    dim1 = 576
    dim2 = 361
    dim3 = 72

    o3_300hpa_all = fltarr(576, 361, ndays)
    o3_300hpa_all[*] = !values.f_nan

    lat_all = fltarr(576, 361, ndays)
    lat_all[*] = !values.f_nan
    lon_all = fltarr(576, 361, ndays)
    lon_all[*] = !values.f_nan

    for iday=0, ndays-1 do begin
      print, 'iday', iday
      caldat, jd_list[iday], month, day, year, hour, minute
      yyyy = string(year, format='(i04)')
      mm = string(month, format='(i02)')
      dd = string(day, format='(i02)')
      hh = string(hour, format='(i02)')
      mi = string(minute, format='(i02)')

      merra2asmfn = '/data/MODEL/MERRA2/2020/10/MERRA2_400.tavg3_3d_asm_Nv.' + yyyy + mm + dd + '.nc4'
      merra2 = ds_read_merra2_tavg3_asm(merra2asmfn)

      datetime_str = yyyy + mm + dd + '_' + hh + mi

      saveflag = 1
        
      ds_merra2_accum, merra2, merra2_o3_300hpa, itime=1, hpa=300

      o3_all[0:dim1-1, 0:dim2-1, iday] = merra2_o3_300hpa
      lon_all[0:dim1-1, 0:dim2-1, iday] = merra2.lon2d
      lat_all[0:dim1-1, 0:dim2-1, iday] = merra2.lat2d

      o3_300hpa_all[0:dim1-1, 0:dim2-1, iday] = merra2_o3_300hpa
    endfor

    if saveflag then begin
      monthly_300hpa_o3 = mean(o3_300hpa_all, dimension=3, /nan)
      ;monthly_300hpa_o3_num = ncount
      monthly_300hpa_o3_stddev = stddev(o3_300hpa_all, dimension=3, /nan)
      
      savepath = '~/data/monthly_tropo_o3/'

      lon2d = merra2.lon2d
      lat2d = merra2.lat2d
      save, filename = savepath + yyyy + mm + '_12utc.sav', $
        lon2d, lat2d, $
        monthly_300hpa_o3, monthly_300hpa_o3_stddev
    endif
    saveflag = 0
  endfor
endfor
end

