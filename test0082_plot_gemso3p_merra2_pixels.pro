
clima = ['ML']
atmos = ['FNL']

for iwli=310, 310, 5 do begin 
  for iclima=0, 0 do begin
    for iatmos=0, 0 do begin
      wli_str = strtrim(string(iwli), 2)

      ;ds_plot_gemso3p_scene, '/data/private/soodal/FNL/2x_ozmin_ozmax_fnl_daily_' + wli_str + '_340/' + $
        ;'GK2_GEMS_L2_O3P_20210621_0345_' + clima[iclima] + $
        ;'_prec000000_maxit10_ozminmax2x_BIN4x4.nc', $
        ;outputpath='./plot/', project='clima' + clima[iclima] + '_' + atmos[iatmos] + '_' + wli_str + '-340'
      ds_plot_gemso3p_pixel, '/data/private/soodal/FNL/2x_ozmin_ozmax_fnl_daily_' + wli_str + '_340/' + $
        'GK2_GEMS_L2_O3P_20210621_0345_' + clima[iclima] + $
        '_prec000000_maxit10_ozminmax2x_BIN4x4.nc', $
        [15, 25, 35, 45, 55, 65], [400, 400, 400, 400, 400, 400], $
        outputpath='./plot/', project='gems_merra2_profile_comp_' + clima[iclima] + '_' + atmos[iatmos] + '_' + wli_str + '-340'

      ds_plot_merra2_collocated_pixel, '/data/private/soodal/MERRA2_collocated_on_gems/' + $
        'merra2_tavg3_3d_asm_collocated_on_gemsl1c_rad_20210621_0345.nc4', $
        [15, 25, 35, 45, 55, 65], [400, 400, 400, 400, 400, 400], $
        outputpath='./plot/', project='gems_merra2_profile_comp_' + clima[iclima] + '_' + atmos[iatmos] + '_' + wli_str + '-340'
    endfor
  endfor
endfor
end
