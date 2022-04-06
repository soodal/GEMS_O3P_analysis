
jd = timegen(start=julday(3, 28, 2021, 0, 45), $
  final=julday(3, 31, 2021, 23, 45), units='Hour')

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


    ds_plot_gemso3p_scene, filelist[0], $
      outputpath='./plot/', $
      project='202109_update_check', sub='old'

    ds_plot_gemso3p_scene, dev2filelist[0], $
      outputpath='./plot/', $
      project='202109_update_check', sub='new'
  endif
endfor

jd = timegen(start=julday(4, 15, 2021, 0, 45), $
  final=julday(4, 19, 2021, 23, 45), units='Hour')

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


    ds_plot_gemso3p_scene, filelist[0], $
      outputpath='./plot/', $
      project='202109_update_check', sub='old'

    ds_plot_gemso3p_scene, dev2filelist[0], $
      outputpath='./plot/', $
      project='202109_update_check', sub='new'
  endif
endfor

jd = timegen(start=julday(4, 27, 2021, 0, 45), $
  final=julday(4, 29, 2021, 23, 45), units='Hour')


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


    ds_plot_gemso3p_scene, filelist[0], $
      outputpath='./plot/', $
      project='202109_update_check', sub='old'

    ds_plot_gemso3p_scene, dev2filelist[0], $
      outputpath='./plot/', $
      project='202109_update_check', sub='new'
  endif
endfor

jd = timegen(start=julday(5, 6, 2021, 0, 45), $
  final=julday(5, 9, 2021, 23, 45), units='Hour')

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


    ds_plot_gemso3p_scene, filelist[0], $
      outputpath='./plot/', $
      project='202109_update_check', sub='old'

    ds_plot_gemso3p_scene, dev2filelist[0], $
      outputpath='./plot/', $
      project='202109_update_check', sub='new'
  endif
endfor
end
