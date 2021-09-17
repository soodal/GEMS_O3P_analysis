function ds_read_gems_l1c_4x4, file

;"/data1/L1C_GEMS/float/1x1/nopolc_eosrl/GK2_GEMS_L1C_20200616_0345_NOR_0694.4x4.nc4"
fid=ncdf_open(file)

varids = ncdf_varidsinq(fid)

;struct = {l1c_4x4, $
  ;bad_pixel_mask:fltarr(174, 512, 1033), $
  ;exposure_time:fltarr(174), $
  ;ground_pixel_quality_flag:intarr(174, 512), $
  ;image_acquisition_time:fltarr(174), $
  ;image_pixel_values:fltarr(174, 512, 1033), $
  ;image_pixel_values_bfpc:fltarr(174, 512, 1033), $
  ;nav_average_att_roll:dblarr(174), $
  ;nav_average_att_pitch:dblarr(174), $
  ;nav_average_att_yaw:dblarr(174), $
  ;nav_inr_performance_time:dblarr(174), $
  ;nav_max_att_pitch:dblarr(174), $
  ;nav_max_att_roll:dblarr(174), $
  ;nav_max_att_yaw:dblarr(174), $
  ;nav_min_att_roll:dblarr(174), $
  ;nav_min_att_pitch:dblarr(174), $
  ;nav_min_att_yaw:dblarr(174), $
  ;nav_stddev_att_roll:dblarr(174), $
  ;nav_stddev_att_pitch:dblarr(174), $
  ;nav_stddev_att_yaw:dblarr(174), $
  ;number_of_inr_performance:dblarr(174), $
  ;pixel_longitude:fltarr(174, 512), $
  ;pixel_latitude:fltarr(174, 512), $
  ;sc_azimuth_angle:fltarr(174, 512), $
  ;sc_zenith_angle:fltarr(174, 512), $
  ;snow_index:intarr(174, 512), $
  ;sun_azimuth_angle:fltarr(174, 512), $
  ;sun_zenith_angle:fltarr(174, 512), $
  ;terrain_height:intarr(174, 512), $
  ;wavelength:fltarr(512, 1033), $
  ;wavelength_bfcal:fltarr(512, 1033), $
  ;xtrack_quality_flag:fltarr(174, 512)}

for i=0, n_elements(varids)-1 do begin
  _var = ncdf_varinq(fid, i)
  name = _var.name
  type = _var.datatype
  ncdf_varget, fid, i, var1 
  sz = size(var1, /dim)

  info_tmp=ncdf_varinq(fid,varids[i])
  _varname=info_tmp.name

  if n_elements(struct) eq 0 then begin
    print,'Reading ... ', info_tmp.name
    struct= create_struct(name, var1)
  endif else begin
    print,'Reading ... ', name 
    struct= create_struct(struct, name, var1)
  ENDELSE

  print, name, type, sz

endfor
ncdf_close, fid

return,struct

end
