
pro comparison_gems_with_omi, jd_total

total_omivals = []
total_gemsvals = []


caldat, jd_total[0], month, day, year, hour, minute
yyyy = string(year, format='(i04)')
mm = string(month, format='(i02)')
dd = string(day, format='(i02)')
hh = string(hour, format='(i02)')
mi = string(minute, format='(i02)')

start_datetime_str = yyyy+mm + dd + '_' + hh+mi
for ijd = 0, n_elements(jd_total)-1 do begin
  caldat, jd_total[ijd], month, day, year, hour, minute
  yyyy = string(year, format='(i04)')
  mm = string(month, format='(i02)')
  dd = string(day, format='(i02)')
  hh = string(hour, format='(i02)')
  mi = string(minute, format='(i02)')

  datetime_str = yyyy+mm + dd + '_' + hh+mi
  path = '/data2/L2_GEMS.nier/monthly_val/O3P/'
  filepattern = 'GK2_GEMS_L2_' + datetime_str + '_O3P_*.nc' 
  filelist = file_search(path + filepattern)

  ;dev2path = '/data/nier_ftp/DEV2/NEW/L2/O3P/' + yyyy + mm + '/' + dd + '/'
  ;dev2filepattern = 'GK2_GEMS_L2_' + datetime_str + '_O3P_*_*_BIN4x4.nc' 
  ;dev2filelist = file_search(dev2path + dev2filepattern)
  ;gemso3p_outfilelist = file_search($
    ;'/data2/L2_GEMS/val_1008/GK2_GEMS_O3P_2020' + smon + '??_??45.nc4')

  ;monthly_omivals = []
  ;monthly_gemsvals = []
  flag_dailyplot = 0
  for igemsfile=0, 0 do begin ;n_elements(gemso3p_outfilelist)-1 do begin
    omivals = []
    gemsvals = []
    if strlen(filelist[0]) GT 1 then begin 
      collocate_gemsl2o3p_omil2profoz, year, month, day, hour, minute, gemsvals, omivals, $
        hpa=300, $ ; omi profile level altitud dimension = [km]
        gemsfile=filelist[0], $
        savepath = './collocate_gemsl2o3p_omiprofoz_old/'
        ;height=10. ; omi profile level altitud dimension = [km]

      if n_elements(gemsvals) ge 2 and n_elements(omivals) ge 2 then begin
        if flag_dailyplot then begin
          plot_gems_validation, gemsvals, omivals, $
            filename='./plot/gems_l2_o3p_val_with_omi_' + gemsyyyymmdd[igemsfile] + $
              'T' + gemshhmi[igemsfile]+'_under300hpa_ecf1.png', $
              range=[0, 100], $
              delta=2.0
        endif
        ;monthly_omivals = [monthly_omivals, omivals]
        ;monthly_gemsvals = [monthly_gemsvals, gemsvals]
        total_omivals = [total_omivals, omivals]
        total_gemsvals = [total_gemsvals, gemsvals]
      endif
    endif
  endfor

  ;if n_elements(monthly_gemsvals) gt 1 and $
    ;n_elements(monthly_omivals) gt 1 then begin
  ;; plot tropospheric column ozone for tropics lat < 10
    ;plot_gems_validation, monthly_gemsvals, monthly_omivals, $
      ;filename='./plot/gems_l2_o3p_val_with_omi_' + $
        ;strmid(gemsyyyymmdd[0], 0, 6) +'_under300hpa_ecf1.png', $
      ;xtitle='GEMS Tropospheric O3', $
      ;ytitle='OMI Tropospheric O3', $
      ;cblim=[0, 60], $
      ;range=[0, 60], $
      ;delta=2.0
  ;endif
endfor

nanidx = where(total_gemsvals lt 0.1, /null)
total_gemsvals[nanidx] = !values.f_nan
nanidx = where(total_omivals lt 0.1, /null)
total_omivals[nanidx] = !values.f_nan
  ; plot tropospheric column ozone for tropics lat < 10
if n_elements(total_gemsvals) gt 1 and $
  n_elements(total_omivals) gt 1 then begin
  plot_gems_validation, total_omivals, total_gemsvals, $
    filename='./plot/202109_update_check/gems_l2_o3p_old_val_with_omi_' + $
      '202109update_under300hpa_ecf02_' + start_datetime_str + '.png', $
      ;'202109update_under300hpa_ecf02_20210415.png', $
    xtitle='OMI Tropospheric O3 under 300 hPa [DU]', $
    ytitle='GEMS Tropospheric O3 under 300 hPa [DU]', $
    cblim=[0, 1000], $
    range=[0, 60], $
    delta=2.0
endif
end


jd_total = []

jd = timegen(start=julday(3, 28, 2021, 0, 45), $
  final=julday(3, 31, 2021, 23, 45), units='Hour')
jd_total=[jd_total, jd]

comparison_gems_with_omi, jd_total

jd_total = []

jd = timegen(start=julday(4, 15, 2021, 0, 45), $
  final=julday(4, 19, 2021, 23, 45), units='Hour')
jd_total=[jd_total, jd]

comparison_gems_with_omi, jd_total

jd_total = []

jd = timegen(start=julday(4, 27, 2021, 0, 45), $
  final=julday(4, 29, 2021, 23, 45), units='Hour')
jd_total=[jd_total, jd]

comparison_gems_with_omi, jd_total

jd_total = []

jd = timegen(start=julday(5, 6, 2021, 0, 45), $
  final=julday(5, 9, 2021, 23, 45), units='Hour')
jd_total=[jd_total, jd]

comparison_gems_with_omi, jd_total

end
