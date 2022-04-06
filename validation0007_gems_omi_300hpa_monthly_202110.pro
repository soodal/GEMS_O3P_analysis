
total_omivals = []
total_gemsvals = []

for imon = 10, 10 do begin
  jd_list = timegen(start=julday(imon, 1, 2021, 00, 45), $
                    final=julday(imon+1, 1, 2021, 23, 45)-1, units='HOURS')
  monthly_omivals = []
  monthly_gemsvals = []
  flag_dailyplot = 0


  for itime = 0, n_elements(jd_list)-1 do begin 
    caldat, jd_list[itime], month, day, year, hour, minute
    yyyy = string(year, format='(i04)')
    mm = string(month, format='(i02)')
    dd = string(day, format='(i02)')
    hh = string(hour, format='(i02)')
    mi = string(minute, format='(i02)')

    datetime_str = yyyy + mm + dd + '_' + hh + mi
    smon = string(imon, format='(I02)')


    ;if year eq 2021 and imon lt 6 then begin
      ;gemso3p_outfilelist = file_search($
        ;'/data2/L2_GEMS.nier/monthly_val/O3P/GK2_GEMS_L2_' + datetime_str + '_O3P*BIN4x4.nc')
    ;endif else if year eq 2021 and imon ge 6 then begin
      ;gemso3p_outfilelist = file_search($
        ;'/data/nier_ftp/O3P/V03/' + yyyy + mm + '/' + dd + '/' + 'GK2_GEMS_L2_' + datetime_str + '_O3P*BIN4x4.nc')
		;endif else begin
      ;gemso3p_outfilelist = file_search($
        ;'/data2/L2_GEMS/val_1008/' + 'GK2_GEMS_O3P_' + datetime_str + '_*.nc4')
    ;endelse
    gemso3p_outfilelist = file_search($
      '/data/private/soodal/ln/GEMS/L2O3P/GK2_GEMS_L2_O3P_' + datetime_str + '_4x4.nc')


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
          if flag_dailyplot then begin
            plot_gems_validation, gemsvals, omivals, $
              filename='./plot/gems_l2_o3p_val_with_omi_' + gemsyyyymmdd[igemsfile] + $
                'T' + gemshhmi[igemsfile]+'_under300hpa_ecf0.2.png', $
                range=[0, 100], $
                delta=2.0
          endif
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
    plot_gems_validation, monthly_omivals, monthly_gemsvals, $
      filename='./plot/gems_monthly_validation/gems_l2_o3p_val_with_omi_' + $
        strmid(datetime_str[0], 0, 6) +'_under300hpa_ecf0.2.png', $
      title='GEMS O3(Sfc to 300 hPa) Monthly Validation', $
      ytitle='GEMS Tropospheric O3[DU]', $
      xtitle='OMI Tropospheric O3[DU]', $
      cblim=[0, 2000], $
      range=[0, 60], $
      delta=2.0
  endif
endfor
stop

  ; plot tropospheric column ozone for tropics lat < 10
if n_elements(total_gemsvals) gt 1 and $
  n_elements(total_omivals) gt 1 then begin
  plot_gems_validation, total_omivals, total_gemsvals, $
    filename='./plot/gems_l2_o3p_val_with_omi_' + $
      '202104-202108_under300hpa_ecf0.2.png', $
    ytitle='GEMS Tropospheric O3', $
    xtitle='OMI Tropospheric O3', $
    cblim=[0, 2000], $
    range=[0, 60], $
    delta=2.0
endif
end
