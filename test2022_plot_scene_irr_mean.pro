
jd_list = timegen(start=julday(3, 22, 2021, 0, 45), $
  final=julday(3, 22, 2021, 3, 45), units='Hours')

for i=0, n_elements(jd_list)-1 do begin
  jd = jd_list[i]
  caldat, jd, month, day, year, hour, minute
  yyyy = string(year, format='(i04)')
  mm = string(month, format='(i02)')
  dd = string(day, format='(i02)')
  hh = string(hour, format='(i02)')
  mi = string(minute, format='(i02)')

  datetime_str = yyyy + mm + dd+ '_' +  hh+ mi

  caldat, jd-1, month_yesterday, day_yesterday, year_yesterday, hour_yesterday, minute_yesterday
  yyyy_y = string(year_yesterday, format='(i04)')
  mm_y = string(month_yesterday, format='(i02)')
  dd_y = string(day_yesterday, format='(i02)')
  hh_y = string(hour_yesterday, format='(i02)')
  mi_y = string(minute_yesterday, format='(i02)')

  ; plot for daily irr
  fn = '/home/soodal/data/IRR_mean/o3p/GK2_GEMS_L2_O3P_' + datetime_str + '_IRR_mean_BIN4x4.nc'
  filelist = file_search(fn)
  for i = 0, n_elements(filelist)-1 do begin
    ds_plot_gemso3p_scene, filelist[i], $
        outputpath='./plot/', project='IRR_mean', sub='daily'
  endfor

  ; plot for monthly irr
  fn = '/home/soodal/data/IRR_mean/o3p/IRR_mean/GK2_GEMS_L2_O3P_' + datetime_str + '_IRR_mean_BIN4x4.nc'
  filelist = file_search(fn)
  for i = 0, n_elements(filelist)-1 do begin
    ds_plot_gemso3p_scene, filelist[i], $
        outputpath='./plot/', project='IRR_mean', sub='monthly
  endfor
endfor

end
