
jd_list = timegen(start=julday(6, 21, 2021, 3, 45), $
  final=julday(6, 21, 2021, 3, 45), units='Hours')

itime = 0
caldat, jd_list[itime], month, day, year, hour, minute
yyyy = string(year, format='(i04)')
mm = string(month, format='(i02)')
dd = string(day, format='(i02)')
hh = string(hour, format='(i02)')
mi = string(minute, format='(i02)')

datetime_str = yyyy + mm + dd + '_' + hh + mi

search_str = '/home/soodal/data/softcal_test/raw/model_300_340/' + $
  'GK2_GEMS_L2_O3P_' + datetime_str + '_winliminit300_prec000000_climML_b4x4_o20211002T153843Z_clima_maxiter10_ecf0.nc'

fl = file_search(search_str)

if strlen(fl) gt 0 then begin
  ds_plot_gemso3p_scene, fl[-1], $
    outputpath='./plot/', $
    project='softcal_test', sub='merra2_300_340'
endif
end
