function ds_read_gems_l1c, file

;/data2/L1C_GEMS/L1C/GK2_GEMS_L1C_20200907_0745_NOR_694.nc'
fid=ncdf_open(file)

varids = ncdf_varidsinq(fid)

varid = ncdf_varid(fid, 'pixel_latitude')
ncdf_varget, fid, varid, lat

varid = ncdf_varid(fid, 'pixel_longitude')
ncdf_varget, fid, varid, lon

nanidx = where(lon le -400, /null)
lon[nanidx] = !values.f_nan
nanidx = where(lat le -400, /null)
lon[nanidx] = !values.f_nan

struct.pixel_longitude = lon
struct.pixel_latitude = lat

return,struct
;struct = {l1c, image_pixel_values: $
  ;pixel_longitude:fltarr(695, 2048), $
  ;pixel_latitude:fltarr(695, 2048)}

;for i=0, n_elements(varids)-1 do begin
  ;struct = ncdf_varinq(fid, i)
  ;name = struct.name
  ;type = struct.datatype
  ;ncdf_varget, fid, i, var1 
  ;sz = size(var1, /dim)
  ;dummy = execute('struct.' + name + '=var1')
  ;print, name, type, sz

;endfor

;return,struct

end
