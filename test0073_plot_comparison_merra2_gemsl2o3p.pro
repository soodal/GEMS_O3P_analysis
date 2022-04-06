
merra2gemso3pfn = "/data/private/soodal/merra2_residual/residuals/GK2_GEMS_L2_O3P_20210621_0345_wli300_prec000000_climML_b4x4_o20211002T123142Z_ecf0.nc"
;merra2gemso3pfn = "/data/private/soodal/softcal_test/raw/model_300_340/GK2_GEMS_L2_O3P_20210621_0345_winliminit300_prec000000_climML_b4x4_o20211002T153843Z_clima_maxiter10_ecf0.nc"

merra2gemsl2o3p = ds_read_gems_l2_o3p(merra2gemso3pfn)


gemso3pfn = '/home/soodal/data/FNL/2x_ozmin_ozmax_fnl_daily_310_340/GK2_GEMS_L2_O3P_20210621_0345_ML_prec000000_maxit10_ozminmax2x_BIN4x4.nc'

gemsl2o3p = ds_read_gems_l2_o3p(gemso3pfn)

merra2_toz = merra2gemsl2o3p.merra2_toz
nanidx = where(merra2_toz lt -1.e29, /null)
merra2_toz[nanidx] = !values.f_nan


cao3 = gemsl2o3p.columnAmountO3
gems_toz = reform(cao3[*, *, 0])
nanidx = where(gems_toz lt -1.e29, /null)
gems_toz[nanidx] = !values.f_nan

showme, merra2_toz - gems_toz


outputpath = './plot/gemso3p-merra2/'

datetime_str = '20210621_0345'
plot_sat_proj, gems_toz - merra2_toz, gemsl2o3p.longitude, gemsl2o3p.latitude, $
    title='GEMS O3P Total ozone - MERRA2 TOZ' + datetime_str, $
    colortable=72, $
    /ctreverse, $
    range=[-30, 30], $
    pngfile=outputpath + '01_totalozone_difference.png';, $

end
