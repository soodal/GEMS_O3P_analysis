function ds_read_gems_l1c_4x4_lonlat, file

;"/data1/L1C_GEMS/float/1x1/nopolc_eosrl/GK2_GEMS_L1C_20200616_0345_NOR_0694.4x4.nc4"
fid=ncdf_open(file)

varids = ncdf_varidsinq(fid)

struct = {pixel_longitude:fltarr(174, 512), $
  pixel_latitude:fltarr(174, 512) }


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

;for i=0, n_elements(varids)-1 do begin
  ;get_struct = ncdf_varinq(fid, i)
  ;name = get_struct.name
  ;type = get_struct.datatype
  ;ncdf_varget, fid, i, var1 
  ;sz = size(var1, /dim)
  ;if name eq 'pixel_latitude' then begin
    ;dummy = execute('struct.' + name + '=var1')
  ;endif else if name eq 'pixel_longitue' then begin
    ;dummy = execute('struct.' + name + '=var1')
  ;endif
  ;print, name, type, sz

;endfor

;return,struct

end
