function ds_read_gems_l2_cld, file, varlist=varlist
varlist2get = []
if not keyword_set(varlist) then begin
  varlist = []
endif

fid=ncdf_open(file)

p1id=ncdf_groupsinq(fid)
p1name=ncdf_groupname(p1id[0])

if p1name eq 'NETCDF4' then begin
  p2id=ncdf_groupsinq(p1id[0])
  p3id=ncdf_groupsinq(p2id[0])
  p4id=ncdf_groupsinq(p3id[0])

ENDIF else begin
  p4id=p1id
ENDELSE

;ngroup = n_elements(p4id)
for dn=0,1 do begin
  id=p4id[dn]
  varid=ncdf_varidsinq(id)
  geopath=ncdf_fullgroupname(p4id[0])
  datpath=ncdf_fullgroupname(p4id[1])

  
  
  if dn eq 0 THEN begin ; geopath
    for vn=0,n_elements(varid)-1 do begin
      info_tmp=ncdf_varinq(id,varid[vn])
      _varname=info_tmp.name

      if n_elements(_varname) ne 0 or n_elements(where(_varname eq varlist, /null)) ne 0 then begin
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(_varname, data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, _varname, data)
        ENDELSE
      endif
    ENDFOR  ; variable
  ENDIF ELSE BEGIN ; datpath 
    FOR vn=0,n_elements(varid)-1 DO BEGIN
      info_tmp=ncdf_varinq(id,varid[vn])
      _varname=info_tmp.name

      if n_elements(_varname) ne 0 or n_elements(where(_varname eq varlist, /null)) ne 0 then begin
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(_varname, data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, _varname, data)
        ENDELSE
      endif
    ENDFOR
  ENDELSE
ENDFOR  ; group

ncdf_close, fid
return, gems

end
