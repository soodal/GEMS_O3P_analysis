
;fn = '/home/soodal/data/softcal_test/corrected/model_310_340/GK2_GEMS_L2_O3P_20210621_0345_winliminit310_prec000000_climML_b4x4_o20211116T130448Z_clima_maxiter10_ecf0.nc'

;ds_plot_gemso3p_scene, fn, $
    ;outputpath='./plot/', project='corrected_310_340'

;fn = '/home/soodal/data/softcal_test/raw/model_300_340/GK2_GEMS_L2_O3P_20210621_0345_winliminit300_prec000000_climML_b4x4_o20210926T114307Z_clima_maxiter10.nc'
;fn = '/home/soodal/data/softcal_test/corrected/model_300_340/GK2_GEMS_L2_O3P_20210329_0345_winliminit300_prec000000_climML_b4x4_o20211122T154430Z_clima_maxiter10_ecf0.nc'
fn = '/home/soodal/data/softcal_test/corrected/model_300_340/GK2_GEMS_L2_O3P_20210329_0345_winliminit300_prec000000_climML_b4x4_o20211129T141817Z_clima_maxiter10_ecf0.nc'

ds_plot_gemso3p_scene, fn, $
    outputpath='./plot/', project='softcal_test', sub='merra2_300_340'

fn = '/home/soodal/data/softcal_test/corrected/model_310_340/GK2_GEMS_L2_O3P_20210329_0345_winliminit310_prec000000_climML_b4x4_o20211130T000233Z_clima_maxiter10_ecf0.nc'

ds_plot_gemso3p_scene, fn, $
    outputpath='./plot/', project='softcal_test', sub='merra2_310_340'

end
