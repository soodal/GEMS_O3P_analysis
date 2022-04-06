
jd_list = timegen(start=julday(3, 29, 2021, 3, 45), $
  final=julday(3, 29, 2021, 3, 45), units='Hours')

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


  fn = '/home/soodal/data/softcal_test/corrected/model_310_340/GK2_GEMS_L2_O3P_20210329_0345_winliminit310_prec000000_climML_b4x4_o20211222T021345Z_clima_maxiter10_ecf0.nc'

  ;filelist = file_search('/data/private/soodad/ln/GEMS/L2CLD/GK2_GEMS_L2_CLD_*_4x4.nc')
  filelist = file_search(fn)

  for i = 0, n_elements(filelist)-1 do begin

    ds_plot_gemso3p_scene, filelist[i], $
        outputpath='./plot/', project='merra2_synthetic', sub='310_340'
  endfor
endfor

end
