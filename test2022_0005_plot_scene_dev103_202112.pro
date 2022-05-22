
jd_list = timegen(start=julday(11, 30, 2021, 19, 45), $
  final=julday(12, 31, 2021, 10, 45), units='Hours')

for i=0, n_elements(jd_list)-1 do begin
  jd = jd_list[i]
  caldat, jd, month, day, year, hour, minute
  yyyy = string(year, format='(i04)')
  mm = string(month, format='(i02)')
  dd = string(day, format='(i02)')
  hh = string(hour, format='(i02)')
  mi = string(minute, format='(i02)')
  ;print, yyyy, mm, dd,'_', hh, mi

  datetime_str = yyyy + mm + dd+ '_' +  hh+ mi


  fn = '/data/nier_ftp/DEV2/O3P/V1.0.3/' + yyyy + mm + '/' + dd + '/GK2_GEMS_L2_' + datetime_str + '_O3P*BIN4x4.nc'
  filelist = file_search(fn)

  for j = 0, n_elements(filelist)-1 do begin
    ds_plot_gemso3p_scene, filelist[j], $
        outputpath='./plot/', project='DEV1.0.3'
  endfor
endfor

end
