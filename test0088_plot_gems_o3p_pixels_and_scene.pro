
  cities_name = ['Singapore', 'Kuala Lumpur', 'Hanoi', 'Hong Kong', 'Naha', $
    'Pohang']

  cities_lon = [103.9, 101.70, 105.8, 114.10, 127.7, 129.2]
  cities_lat = [1.30, 2.7, 21.0, 22.3, 26.2, 36.0]

;ds_plot_gemso3p_scene, '/data/private/soodal/test/which_slit6_gemsslit/' + $
  ;'GK2_GEMS_L2_O3P_20210329_0345_wli310_ML_prec000000_0.5_ozminmax2x_BIN4x4.nc', $ 
  ;outputpath='./plot/', project='test', sub='which_slit6_gemsslit'
ds_plot_gemso3p_pixels_with_latlon, '/data/private/soodal/test/which_slit6_gemsslit/' + $
  'GK2_GEMS_L2_O3P_20210329_0345_wli310_ML_prec000000_0.5_ozminmax2x_BIN4x4.nc', $
  cities_lon, cities_lat, $
  outputpath='./plot/', project='test', sub='which_slit6_gemsslit', name=cities_name
end
