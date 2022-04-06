
fn = '/home/soodal/data/softcal_test/corrected/model_310_340/GK2_GEMS_L2_O3P_20210621_0345_winliminit310_prec000000_climML_b4x4_o20211116T130448Z_clima_maxiter10_ecf0.nc'

filelist = file_search('/data/private/soodal/ln/GEMS/L2CLD/GK2_GEMS_L2_CLD_20210329_*_4x4.nc')

for i = 0, n_elements(filelist)-1 do begin

  ds_plot_gemsl2cld_scene, filelist[i], $
      outputpath='./plot/', project='GEMS_L2_CLD_ecf'
endfor

end
