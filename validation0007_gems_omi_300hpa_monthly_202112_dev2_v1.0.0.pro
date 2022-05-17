
total_omivals = []
total_gemsvals = []

for imon = 12, 12 do begin
  jd_list = timegen(start=julday(imon, 2, 2021, 00, 45), $
                    final=julday(imon, 14, 2021, 23, 45)-1, units='HOURS')
  monthly_omivals = []
  monthly_gemsvals = []
  ;flag_dailyplot = 0


  for itime = 0, n_elements(jd_list)-1 do begin 
    caldat, jd_list[itime], month, day, year, hour, minute
    yyyy = string(year, format='(i04)')
    mm = string(month, format='(i02)')
    dd = string(day, format='(i02)')
    hh = string(hour, format='(i02)')
    mi = string(minute, format='(i02)')

    datetime_str = yyyy + mm + dd + '_' + hh + mi
    smon = string(imon, format='(I02)')


    gemso3p_outfilelist = file_search($
      '/data/nier_ftp/O3P/V1.0/' + yyyy+mm +'/' + dd + '/' + $
      'GK2_GEMS_L2_' + datetime_str + '*_BIN4x4.nc')

    datetime_str_pos = stregex(gemso3p_outfilelist, '[0-9]{8}_[0-9]{4}')
    gemsyyyymmdd = strmid(gemso3p_outfilelist, datetime_str_pos, 8)
    gemshhmi = strmid(gemso3p_outfilelist, datetime_str_pos + 9, 4)
    gemsyyyymmdd_hhmi = strmid(gemso3p_outfilelist, datetime_str_pos, 13)

    if strlen(gemso3p_outfilelist[0]) gt 3 then begin 
      for igemsfile=0, 0 do begin

        omivals = []
        gemsvals = []
        collocate_gemsl2o3p_omil2profoz, year, month, day, hour, mi, gemsvals, omivals, $
          hpa=300, $ ; omi profile level altitud dimension = [km]
          gemsfile=gemso3p_outfilelist[igemsfile]
          ;height=10. ; omi profile level altitud dimension = [km]

        if n_elements(gemsvals) ge 2 and n_elements(omivals) ge 2 then begin
          ;if flag_dailyplot then begin
            ;plot_gems_validation, gemsvals, omivals, $
              ;filename='./plot/gems_l2_o3p_val_with_omi_' + gemsyyyymmdd[igemsfile] + $
                ;'T' + gemshhmi[igemsfile]+'_under300hpa_ecf0.2.png', $
                ;range=[0, 100], $
                ;delta=2.0
          ;endif
          monthly_omivals = [monthly_omivals, omivals]
          monthly_gemsvals = [monthly_gemsvals, gemsvals]
          total_omivals = [total_omivals, omivals]
          total_gemsvals = [total_gemsvals, gemsvals]
        endif
      endfor
    endif
  endfor

  if n_elements(monthly_gemsvals) gt 1 and $
    n_elements(monthly_omivals) gt 1 then begin
  ; plot tropospheric column ozone for tropics lat < 10
    validx = where(monthly_omivals ge 0.01 and monthly_gemsvals ge 0.01, valnum, /null)
    if valnum ge 1 then begin
      plot_gems_validation, monthly_omivals[validx], monthly_gemsvals[validx], $
        filename='./plot/gems_monthly_validation/gems_l2_o3p_val_with_omi_' + $
          strmid(datetime_str[0], 0, 6) +'_under300hpa_ecf0.2_dev2_20211202_20211214_v1.0.0.png', $
        title='GEMS TropO3(Sfc - 300 hPa)' , $
        ytitle='GEMS Tropospheric O3[DU]', $
        xtitle='OMI Tropospheric O3[DU]', $
        cblim=[0, 2000], $
        range=[0, 60], $
        delta=2.0
    endif else begin
      plot_gems_validation, monthly_omivals, monthly_gemsvals, $
        filename='./plot/gems_monthly_validation/gems_l2_o3p_val_with_omi_' + $
          strmid(datetime_str[0], 0, 6) +'_under300hpa_ecf0.2_dev2_20211202_20211214_v1.0.0.png', $
        title='GEMS TropO3(Sfc - 300 hPa)' , $
        ytitle='GEMS Tropospheric O3[DU]', $
        xtitle='OMI Tropospheric O3[DU]', $
        cblim=[0, 2000], $
        range=[0, 60], $
        delta=2.0
    endelse
  endif
endfor

end
