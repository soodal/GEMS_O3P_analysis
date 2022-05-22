
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

  datetime_str = yyyy + mm + dd+ '_' +  hh+ mi


  fn100 = '/data/nier_ftp/O3P/V1.0/202103/' + dd + '/GK2_GEMS_L2_' + datetime_str + '_O3P*BIN4x4.nc'
  filelist100 = file_search(fn100)

  fn103 = '/data/nier_ftp/DEV2/O3P/V1.0.3/202103/' + dd + '/GK2_GEMS_L2_' + datetime_str + '_O3P*BIN4x4.nc'
  filelist103 = file_search(fn103)

  if file_test(filelist100[0]) and file_test(filelist103[0]) then begin
    ds_plot_gemso3p_scene_diff, filelist100[0], filelist103[0], $
        outputpath='./plot/', project='DIFF_v1.0.0_v1.0.3'
  endif
endfor

end
