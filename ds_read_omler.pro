pro ds_read_omler, omler, path=path
fn = '/OPER/SYSTEMS/ALGRTH/GEMS/L2_2021/data/o3p/ATMOS/ATMOS/KNMI_OMIALB/OMI-Aura_L3-OMLER_2005m01-2009m12_v003-2010m0503t063707.nc'
if keyword_set(path) then begin
  fn = path
endif

fid=ncdf_open(fn)
varids = ncdf_varidsinq(fid)
for i=0, n_elements(varids)-1 do begin
  _var = ncdf_varinq(fid, i)
  name = _var.name
  type = _var.datatype
  ncdf_varget, fid, i, var1 
  sz = size(var1, /dim)

  info_tmp=ncdf_varinq(fid,varids[i])
  _varname=info_tmp.name

  if n_elements(omler) eq 0 then begin
    print,'Reading ... ', info_tmp.name
    omler= create_struct(name, var1)
  endif else begin
    print,'Reading ... ', name 
    omler= create_struct(omler, name, var1)
  ENDELSE
  print, name, type, sz
endfor
ncdf_close, fid
return
end
