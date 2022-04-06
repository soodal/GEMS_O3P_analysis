fn = '/home/soodal/data/softcal_test/raw/model_310_340/GK2_GEMS_L2_O3P_20210621_0345_winliminit310_prec000000_climML_b4x4_o20210924T011905Z_clima_maxiter10.nc'

for iy = 10, 510, 50 do begin
for ix = 30, 140, 50 do begin
  ds_plot_gemso3p_pixel, fn, ix, iy, outputpath = './plot/', project='softcal_test', sub='310_340'
endfor
endfor
end
