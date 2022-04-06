
jd_list = timegen(start=julday(7, 27, 2021, 3, 45), $
  final=julday(7, 27, 2021, 6, 45), units='Hours')

for i=0, n_elements(jd_list)-1 do begin
  jd = jd_list[i]
  caldat, jd, month, day, year, hour, minute
  yyyy = string(year, format='(i04)')
  mm = string(month, format='(i02)')
  dd = string(day, format='(i02)')
  hh = string(hour, format='(i02)')
  mi = string(minute, format='(i02)')

  datetime_str = yyyy + mm + dd+ '_' +  hh+ mi


  for j = 0, 1 do begin
    folders = ['v1.0.0', 'v1.0.3']
    fn = '/home/soodal/data/update_comparison/' + folders[j] + '/O3P/GK2_GEMS_L2_' + datetime_str + '_O3P_FW_DPRO_BIN4x4.nc'

    filelist = file_search(fn)

    for ifile = 0, n_elements(filelist)-1 do begin
      ds_plot_gemso3p_scene, filelist[ifile], $
          outputpath='./plot/', project='update_comparison', sub=folders[j]
    endfor
  endfor
endfor

end
