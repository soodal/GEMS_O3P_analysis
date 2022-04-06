
old_accum_total = []
new_accum_total = []

flag_dailyplot = 0

old_ecfs = []
new_ecfs = []


jd = timegen(start=julday(3, 28, 2021, 0, 45), $
  final=julday(3, 31, 2021, 23, 45), units='Hour')

jd = timegen(start=julday(4, 15, 2021, 0, 45), $
  final=julday(4, 19, 2021, 23, 45), units='Hour')

jd = timegen(start=julday(4, 27, 2021, 0, 45), $
  final=julday(4, 29, 2021, 23, 45), units='Hour')

jd = timegen(start=julday(5, 6, 2021, 0, 45), $
  final=julday(5, 9, 2021, 23, 45), units='Hour')

savfile = './update_check_20210427-20210429.sav'
savfile = './update_check_20210506-20210509.sav'

if not file_test(savfile) then begin 
  for iday=0, n_elements(jd)-1 do begin
    caldat, jd[iday], month, day, year, hour, minute
    yyyy = string(year, format='(i04)')
    mm = string(month, format='(i02)')
    dd = string(day, format='(i02)')
    hh = string(hour, format='(i02)')
    mi = string(minute, format='(i02)')

    datetime_str = yyyy + mm + dd + '_' + hh + mi


    path = '/data2/L2_GEMS.nier/monthly_val/O3P/'
    filepattern = 'GK2_GEMS_L2_' + datetime_str + '_O3P_*.nc' 
    filelist = file_search(path + filepattern)

    dev2path = '/data/nier_ftp/DEV2/NEW/L2/O3P/' + yyyy + mm + '/' + dd + '/'
    dev2filepattern = 'GK2_GEMS_L2_' + datetime_str + '_O3P_*_*_BIN4x4.nc' 
    dev2filelist = file_search(dev2path + dev2filepattern)
    if strlen(filelist[0]) ge 10 and strlen(dev2filelist[0]) ge 10 then begin

      old_gemsval = ds_read_gems_l2_o3p(filelist[0])
      old_ecf = old_gemsval.EffectiveCloudFractionUV
      old_ecfs = [old_ecfs, old_ecf]

      new_gemsval = ds_read_gems_l2_o3p(dev2filelist[0])
      new_ecf = new_gemsval.EffectiveCloudFractionUV
      new_ecfs = [new_ecfs, new_ecf]

      ds_gems_l2o3p_accum, old_gemsval, old_accum, hpa=300
      old_accum_total = [old_accum_total, old_accum]
      ds_gems_l2o3p_accum, new_gemsval, new_accum, hpa=300
      new_accum_total = [new_accum_total, new_accum]

      if flag_dailyplot then begin
        plot_gems_validation, new_accum, old_accum_total, $
          filename='./plot/gems_l2_o3p_val_with_omi_' + $
            yyyy + mm + dd + $
            'T' + hh + mi + '_under10km_ecf0.2_30min.png', $
            range=[0, 60], $
            delta=2.0
      endif
    endif

    ;validx = where(new_ecf gt 0.2 and old_ecf gt 0.2 and new_accum gt 5. and old_accum gt 5., /null)
    ;validx = where(old_accum ge 0 and old_accum lt 1000 $
    ;and new_accum ge 0 and new_accum lt 1000 $


    ;plot_gems_validation, new_accum[validx], old_accum[validx], $
      ;filename='./plot/202109_update_comparison/update_check_gems_l2_o3p_val_with_new_' + $
        ;yyyy + mm + dd + 'T' + hh + mi + 'Z_under300hpa_ecf0.2_30min.png', $
      ;xtitle='GEMS new Tropospheric O3', $
      ;ytitle='GEMS old Tropospheric O3', $
      ;cblim=[0, 200], $
      ;range=[0, 60], $
      ;delta=2.0
  endfor
  save, filename=savfile, new_accum_total, old_accum_total, new_ecfs
endif else begin
  restore, savfile
ENDELSE

; plot tropospheric column ozone 


validx = where(new_ecfs lt 0.2 and old_ecfs lt 0.2 and new_accum_total gt 5. and old_accum_total gt 5., /null)
plot_gems_validation, new_accum_total[validx], old_accum_total[validx], $
  ;filename='./plot/202109_update_comparison/update_check_' + $
    ;'20210427T0045-20210429T2345_under300hpa_ecf0.2_30min.png', $
  filename='./plot/202109_update_comparison/update_check_' + $
    '20210506T0045-20210509T2345_under300hpa_ecf0.2_30min.png', $
  xtitle='GEMS NEW 300 hPa Accumulation O3[DU]', $
  ytitle='GEMS OLD 300 hPa Accumulation O3[DU]', $
  cblim=[0, 2000], $
  range=[0, 60], $
  ;cbtitle='Ozone[DU]', $
  delta=2.0
end
