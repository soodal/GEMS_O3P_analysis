fn = '/data/nier_ftp/DEV2/O3P/V1.0.3/202103/29/GK2_GEMS_L2_20210329_0345_O3P_FW_DPRO_BIN4x4.nc'
o3p = ds_read_gems_l2_o3p(fn)

plot_sat_proj, o3p.TropopausePressure, o3p.longitude, o3p.latitude, $
  pngfile='./plot/fnl_tp/um_tp_20210329.png', $
  range=[50, 350], $
  title='UM calculated Tropopause Pressure'
end
