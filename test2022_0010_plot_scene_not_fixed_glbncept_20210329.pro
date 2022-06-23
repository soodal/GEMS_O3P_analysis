
jd_list = timegen(start=julday(3, 28, 2021, 22, 45), $
  final=julday(3, 29, 2021, 7, 45), units='Hours')

for i=0, n_elements(jd_list)-1 do begin
  jd = jd_list[i]
  caldat, jd, month, day, year, hour, minute
  yyyy = string(year, format='(i04)')
  mm = string(month, format='(i02)')
  dd = string(day, format='(i02)')
  hh = string(hour, format='(i02)')
  mi = string(minute, format='(i02)')

  datetime_str = yyyy + mm + dd+ '_' +  hh+ mi

  cities_name = ['Singapore', 'Kuala Lumpur', 'Hanoi', 'Hong Kong', 'Naha', 'Pohang']

  cities_lon = [103.9, 101.70, 105.8, 114.10, 127.7, 129.2]
  cities_lat = [1.30, 2.7, 21.0, 22.3, 26.2, 36.0]

  clima_list = ['TB', 'ML']
  ;wli_list = ['300', '305', '310']
  wli_list = ['310']
  wle = '340'
  for ic = 0, n_elements(clima_List)-1 do begin
    for iwli = 0, n_elements(wli_list)-1 do begin

      file_str = '/data/private/soodal/test/not_fixed_glbncept/GK2_GEMS_L2_' + $
        yyyy + mm + dd + '_' + hh + mi + '_O3P_winlim' + wli_list[iwli] + wle +'_' + clima_list[ic] + '_prec000000_BIN4x4.nc'

      file_list = file_search(file_str)
      print, file_str

      if file_test(file_list[0]) then begin
        ds_plot_gemso3p_scene, file_list[0], troprange=[20, 60], $
            outputpath='./plot/', project='not_fixed_glbncept', $
            suffix=clima_list[ic]+'_'+wli_list[iwli]+wle
      endif
    endfor
  endfor
endfor

end
