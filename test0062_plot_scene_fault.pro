
clima = ['ML', 'TB']
atmos = ['FNL', 'UM']

jd_list = timegen(start=julday(4, 6, 2021, 0, 45), $
  final=julday(4, 8, 2021, 23, 45), units='Hours')

;jd_list = timegen(start=julday(5, 5, 2021, 0, 45), $
  ;final=julday(5, 8, 2021, 23, 45), units='Hours')

;jd_list = timegen(start=julday(8, 19, 2021, 0, 45), $
  ;final=julday(8, 27, 2021, 23, 45), units='Hours')

for itime=0, n_elements(jd_list)-1 do begin
  caldat, jd_list[itime], month, day, year, hour, minute
  yyyy = string(year, format='(i04)')
  mm = string(month, format='(i02)')
  dd = string(day, format='(i02)')
  hh = string(hour, format='(i02)')
  mi = string(minute, format='(i02)')

  datetime_str = yyyy + mm + dd + '_' + hh + mi

  ;202106~
  ;search_str = '/data/private/soodal/fault_check/2x_ozmin_ozmax_umdaily/' + $
    ;'GK2_GEMS_L2_O3P_' + datetime_str + '_*BIN4x4.nc'

  search_str = '/data2/L2_GEMS.nier/monthly_val/O3P/' + $
    'GK2_GEMS_L2_' + datetime_str + '_O3P_*BIN4x4.nc'

  fl = file_search(search_str)
  ;print, search_str

  if strlen(fl) gt 0 then begin
    ds_plot_gemso3p_scene, fl[-1], $
      outputpath='./plot/', $
      project='check_fault';, sub='2x_ozmin_ozmax'
  endif
endfor
end
