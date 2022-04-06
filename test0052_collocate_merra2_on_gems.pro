

;times = timegen(start=julday(12, 21, 2020, 0, 45), final=julday(5, 19, 2020, 5, 45))
times = timegen(start=julday(8, 1, 2020, 0, 45), final=julday(6, 30, 2021, 7, 45), units='Hour')
for idate=0, n_elements(times)-1 do begin
  ;merra2asm_fn = '/data/MERRA2/2020/08/MERRA2_400.tavg3_3d_asm_Nv.20200810.nc4'
  caldat, times[idate], m1, d1, y1, h1, mi1
  m1s = string(m1, format='(i02)')
  d1s = string(d1, format='(i02)')
  y1s = string(y1, format='(i04)')
  h1s = string(h1, format='(i02)')
  mi1s = string(mi1, format='(i02)')

  ;merra2chm_fn = '/data/MERRA2/' + y1s + '/' + m1s + '/MERRA2_400.inst3_3d_chm_Nv.' + y1s + m1s + d1s '.nc4'
  ;gems_path = '/data/private/soodal/merra2_residual/residuals/'
  ;gems_path = '/data2/L1C_GEMS/L1C/4x4/'
  ;gems_path = '/data2/L1C_GEMS/L1C_nier/L1C/202103/29/GK2_GEMS_L1C_20210329_0045_CAPO_NOR_625_BIN4x4.nc/4x4/'
  gems_path = '/data2/L1C_GEMS/L1C_nier/L1C/' + y1s + m1s + '/' + d1s + '/'
  gems_path = '/home/soodal/data/ln/GEMS/L1C/'

  gems_fn = 'GK2_GEMS_L1C_' $
    + y1s + m1s + d1s + '_' + h1s + mi1s + '_4x4.nc'
  collocate_merra2_tavg3_3d_asm_on_gemsl1c_rad, y1, m1, d1, h1, mi1, gems_path =gems_path, gems_fn_pattern=gems_fn
ENDFOR



  ;merra2asm = ds_read_merra2_tavg3_3d_asm_nv(merra2asm_fn)
  ;merra2chm = ds_read_merra2_inst3_3d_chm_nv(merra2chm_fn)

  ;merra2chm_lon = merra2chm.lon # Replicate(1, n_elements(merra2chm.lat))
  ;merra2chm_lat = replicate(1, n_elements(merra2chm.lon)) # merra2chm.lat

  ;gems_roi_idx = where(merra2chm_lon ge 70 and merra2chm_lon lt 160 $
    ;merra2chm_lat ge -15 and merra2chm_lat lt 55 )

  ;gems_roi_indices = array_indices(merra2chm_lon, gems_roi_idx)
end
