
clima = ['ML', 'TB']
atmos = ['FNL']

for iwli=310, 310, 5 do begin 
  for iclima=0, 1 do begin
    for iatmos=0, 0 do begin
      wli_str = strtrim(string(iwli), 2)

      ds_plot_gemso3p_scene, '/data/private/soodal/FNL/2x_ozmin_ozmax_fnl_daily_' + wli_str + '_340/' + $
        'GK2_GEMS_L2_O3P_20210621_0345_' + clima[iclima] + $
        '_prec000000_maxit10_ozminmax2x_BIN4x4.nc', $
        outputpath='./plot/', project='clima' + clima[iclima] + '_' + atmos[iatmos] + '_' + wli_str + '-340'
      ds_plot_gemso3p_pixel, '/data/private/soodal/FNL/2x_ozmin_ozmax_fnl_daily_' + wli_str + '_340/' + $
        'GK2_GEMS_L2_O3P_20210621_0345_' + clima[iclima] + $
        '_prec000000_maxit10_ozminmax2x_BIN4x4.nc', [84, 84, 84, 84, 84, 84, 84], [24, 44, 84, 133, 350, 400, 450], $
        outputpath='./plot/', project='clima' + clima[iclima] + '_' + atmos[iatmos] + '_' + wli_str + '-340'
    endfor
  endfor
endfor
end
