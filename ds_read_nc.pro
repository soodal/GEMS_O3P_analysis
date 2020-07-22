function ds_read_nc, fn, subpath

nid=NCDF_OPEN(fn)
ncdf_varget, nid, subpath, var
return, var

;file_quiry=NCDF_INQUIRE(nid)
;nvars=file_quiry.nvars
;FOR i=0,nvars-1 DO BEGIN
  ;_var_struct=NCDF_VARINQ(nid,i)
  ;_var_name=_var_struct.name
  ;_var_ndim=_var_struct.ndims
  ;dum=EXECUTE('NCDF_VARGET,nid,_var_name,_var')

  ;dum=execute('data.'+_var_name + '=_var')
;ENDFOR

;print,data
;stop

;FOR ip=0,np-1 DO BEGIN
  ;var_name=param[ip]
  ;NCDF_VARGET,nid,var_name,var
  ;NCDF_ATTGET,nid,var_name,'_FillValue',fill
  ;NCDF_ATTGET,nid,var_name,'units',dumunit
  ;bad=WHERE(var EQ fill, nbad)
  ;IF(nbad GT 0)THEN begin
  ;var[bad]=!values.f_nan
  ;units[ip]=STRING(dumunit)
  ;dum=EXECUTE(param[ip]+'=FLOAT(var)')
  ;var=!null
;ENDFOR
                                                             

end

