
fn = '/data/nier_ftp/O3T/V03/202106/21/GK2_GEMS_L2_20210621_0345_O3T_FW_DPRO_BIN4x4.nc'

filelist = file_search(fn)

for i = 0, n_elements(filelist)-1 do begin

  ds_plot_gems_l2_o3t_scene, filelist[i], $
      outputpath='./plot/', project='GEMS_L2_O3T'
endfor

end
