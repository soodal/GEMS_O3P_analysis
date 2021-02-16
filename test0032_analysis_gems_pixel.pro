
;gemsfn = '/data1/gems/o3p/works/GEMS_O3P_Yonsei/out/GK2B_GEMS_L2_O3P_20200616_0345_snr_climML_winliminit310_2020-12-29T223235.nc4.debug.nc4'
;gemsfn = '/app/gemsl2_2020_0714/src/o3p/v1.b/out/GK2B_GEMS_L2_20200616_0345_O3P_winliminit310_prec0.005_climLLM_4x4_2020-12-23T2200KST.nc4.debug.nc4'

; 20200616_0345 which_prec = 0
gemsfn = '/data1/gems/o3p/works/GEMS_O3P_Yonsei/out/GK2B_GEMS_L2_O3P_20200616_0345_fixed0.005_climML_winliminit310_2020-12-30T134832.nc4.debug.nc4'

; 20200616_0345 which_prec = 10
gemsfn = '/data1/gems/o3p/works/GEMS_O3P_Yonsei/out/GK2B_GEMS_L2_O3P_20200616_0345_snr_climML_winliminit310_2020-12-30T134930.nc4.debug.nc4'


; 20200616_0345 which_prec = 10
; winlim 305340
gemsfn = '/data1/gems/o3p/works/GEMS_O3P_Yonsei/out/snr/GK2B_GEMS_L2_O3P_20200616_0345_snr_climML_winliminit305_2020-12-30T140605.nc4.debug.nc4'


xidx = 89 
;yidx = 299;427
yidx = 427
plot_o3p_profile, gemsfn, xidx, yidx
plot_fitweights, gemsfn, xidx, yidx
plot_fitres, gemsfn, xidx, yidx
plot_gems_snr_pix, gemsfn, xidx, yidx

end
