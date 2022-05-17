function ds_read_fnl_nc, ncfile, varlist=varlist
;ncfile = '/data1/gems/o3p/ds/GEMS_O3P_analysis/fnl_20200616_00_00.grib2.nc'
if not keyword_set(varlist) then begin
  varlist = ['lv_ISBL0', 'lon_0', 'lat_0', 'TMP_P0_L100_GLL0']
endif

fid=ncdf_open(ncfile)

varid=ncdf_varidsinq(fid)

for vn=0,n_elements(varid)-1 do begin
  info_tmp=ncdf_varinq(fid,varid[vn])
  _varname=info_tmp.name
  if n_elements(_varname) ne 0 eq n_elements(where(_varname eq varlist, /null)) ne 0 then begin
    if n_elements(fnl) eq 0 then begin
      print,'Reading ... ',info_tmp.name
      ncdf_varget,fid,varid[vn],data
      fnl = create_struct(_varname, data)
    endif else begin
      print,'Reading ... ',info_tmp.name
      ncdf_varget,fid,varid[vn],data
      fnl = create_struct(fnl, _varname, data)
    ENDELSE
  endif
ENDFOR  ; variable

ncdf_close, fid
return, fnl

end
