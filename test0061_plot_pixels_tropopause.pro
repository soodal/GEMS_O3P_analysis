
clima = ['ML', 'TB']
atmos = ['FNL', 'UM']

for iclima=0, 1 do begin
  for iatmos=0, 1 do begin

    ds_plot_gemso3p_scene, '/data/private/soodal/vartest/clima' + clima[iclima] + $
      '_' + atmos[iatmos] + '/GK2_GEMS_L2_O3P_20210621_0345_' + clima[iclima] + $
      '_BIN4x4.nc', outputpath='./plot/', project='clima' + clima[iclima] + '_' + atmos[iatmos] 
    ds_plot_gemso3p_pixel, '/data/private/soodal/vartest/clima' + clima[iclima] + $
      '_' + atmos[iatmos] + '/GK2_GEMS_L2_O3P_20210621_0345_' + clima[iclima] + $
      '_BIN4x4.nc', [84, 84, 84, 84], [24, 44, 84, 133], $
      outputpath='./plot/', project='clima' + clima[iclima] + '_' + atmos[iatmos]
  endfor
endfor
end
