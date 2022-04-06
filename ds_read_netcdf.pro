function ds_read_netcdf, file, varlist=varlist
if not keyword_set(varlist) then begin
  varlist = []
endif

fid=ncdf_open(file)

varid=ncdf_varidsinq(fid)

for vn=0,n_elements(varid)-1 do begin
  info_tmp=ncdf_varinq(fid,varid[vn])
  _varname=info_tmp.name

  if n_elements(_varname) ne 0 or n_elements(where(_varname eq varlist, /null)) ne 0 then begin
    if n_elements(gems) eq 0 then begin
      print,'Reading ... ',info_tmp.name
      ncdf_varget,fid,varid[vn],data
      gems = create_struct(_varname, data)
    endif else begin
      print,'Reading ... ',info_tmp.name
      ncdf_varget,fid,varid[vn],data
      gems = create_struct(gems, _varname, data)
    ENDELSE
  endif
ENDFOR  ; variable

ncdf_close, fid
return, gems

end
