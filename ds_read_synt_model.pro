function ds_read_synt_model, file, varlist=varlist
varlist2get = []
if not keyword_set(varlist) then begin
  varlist = ['Altitude', 'AOD', $
             'CldFrac', 'CldP', 'ColAmt' , $
             'O3', $
             'ColAmtNO2', 'ColAmtOz','ColAmtSO2', $
             'OzRet','Pressure','SurfaceAlbedo','Temperature']
endif

fid=ncdf_open(file)

n_var = n_elements(varlist)

varid=ncdf_varidsinq(fid)

for vn=0,n_elements(varid)-1 do begin
  info_tmp=ncdf_varinq(fid,varid[vn])
  varname=info_tmp.name

  if varname eq 'Altitude' $
    and n_elements(where(varname eq varlist, /null)) ne 0 then begin
    if n_elements(gems) eq 0 then begin
        print,'Reading ... ',info_tmp.name
        ncdf_varget,fid,varid[vn],data
        gems = create_struct('Altitude', data)
    endif else begin
      print,'Reading ... ',info_tmp.name
      ncdf_varget,fid,varid[vn],data
      gems = create_struct(gems, 'Altitude', data)
    ENDELSE
  endif else if varname eq 'AOD' $
    and n_elements(where(varname eq varlist, /null)) ne 0 then begin
    if n_elements(gems) eq 0 then begin
      print,'Reading ... ',info_tmp.name
      ncdf_varget,fid,varid[vn],data
      nanidx = where(data lt -990, /null)
      data[nanidx] = !values.f_nan
      gems = create_struct('AOD', data)
    endif else begin
      print,'Reading ... ',info_tmp.name
      ncdf_varget,fid,varid[vn],data
      nanidx = where(data lt -990, /null)
      data[nanidx] = !values.f_nan
      gems = create_struct(gems, 'AOD', data)
    ENDELSE
  endif else if varname eq 'CldFrac' $
    and n_elements(where(varname eq varlist, /null)) ne 0 then begin
    if n_elements(gems) eq 0 then begin
      print,'Reading ... ',info_tmp.name
      ncdf_varget,fid,varid[vn],data
      gems = create_struct('CldFrac', data)
    endif else begin
      print,'Reading ... ',info_tmp.name
      ncdf_varget,fid,varid[vn],data
      gems = create_struct(gems, 'CldFrac', data)
    ENDELSE
  endif else if varname eq 'CldP' $
    and n_elements(where(varname eq varlist, /null)) ne 0 then begin
    if n_elements(gems) eq 0 then begin
      print,'Reading ... ',info_tmp.name
      ncdf_varget,fid,varid[vn],data
      nanidx = where(data lt -990, /null)
      data[nanidx] = !values.f_nan
      gems = create_struct('CldP', data)
    endif else begin
      print,'Reading ... ',info_tmp.name
      ncdf_varget,fid,varid[vn],data
      nanidx = where(data lt -990, /null)
      data[nanidx] = !values.f_nan
      gems = create_struct(gems, 'CldP', data)
    ENDELSE
  endif else if varname eq 'ColAmt' $
    and n_elements(where(varname eq varlist, /null)) ne 0 then begin
    if n_elements(gems) eq 0 then begin
      print,'Reading ... ',info_tmp.name
      ncdf_varget,fid,varid[vn],data
      gems = create_struct('ColAmt', data)
    endif else begin
      print,'Reading ... ',info_tmp.name
      ncdf_varget,fid,varid[vn],data
      gems = create_struct(gems, 'ColAmt', data)
    ENDELSE
  endif else if varname eq 'ColAmtNO2' $
    and n_elements(where(varname eq varlist, /null)) ne 0 then begin
    if n_elements(gems) eq 0 then begin
      print,'Reading ... ',info_tmp.name
      ncdf_varget,fid,varid[vn],data
      gems = create_struct('ColAmtNO2', data)
    endif else begin
      print,'Reading ... ',info_tmp.name
      ncdf_varget,fid,varid[vn],data
      gems = create_struct(gems, 'ColAmtNO2', data)
    ENDELSE
  endif else if varname eq 'ColAmtOz' $
    and n_elements(where(varname eq varlist, /null)) ne 0 then begin
    if n_elements(gems) eq 0 then begin
      print,'Reading ... ',info_tmp.name
      ncdf_varget,fid,varid[vn],data
      gems = create_struct('ColAmtOz', data)
    endif else begin
      print,'Reading ... ',info_tmp.name
      ncdf_varget,fid,varid[vn],data
      gems = create_struct(gems, 'ColAmtOz', data)
    ENDELSE
  endif else if varname eq 'ColAmtSO2' $
    and n_elements(where(varname eq varlist, /null)) ne 0 then begin 
    if n_elements(gems) eq 0 then begin
      print,'Reading ... ',info_tmp.name
      ncdf_varget,fid,varid[vn],data
      gems = create_struct('ColAmtSO2', data)
    endif else begin
      print,'Reading ... ',info_tmp.name
      ncdf_varget,fid,varid[vn],data
      gems = create_struct(gems, 'ColAmtSO2', data)
    ENDELSE
  endif else if varname eq 'O3' $
    and n_elements(where(varname eq varlist, /null)) ne 0 then begin
    if n_elements(gems) eq 0 then begin
      print,'Reading ... ',info_tmp.name
      ncdf_varget,fid,varid[vn],data
      gems = create_struct('O3', data)
    endif else begin
      print,'Reading ... ',info_tmp.name
      ncdf_varget,fid,varid[vn],data
      gems = create_struct(gems, 'O3', data)
    ENDELSE
  endif else if varname eq 'OzRet' $
    and n_elements(where(varname eq varlist, /null)) ne 0 then begin
    if n_elements(gems) eq 0 then begin
      print,'Reading ... ',info_tmp.name
      ncdf_varget,fid,varid[vn],data
      gems = create_struct('OzRet', data)
    endif else begin
      print,'Reading ... ',info_tmp.name
      ncdf_varget,fid,varid[vn],data
      gems = create_struct(gems, 'OzRet', data)
    ENDELSE
   endif else if varname eq 'Pressure' $
     and n_elements(where(varname eq varlist, /null)) ne 0 then begin
    if n_elements(gems) eq 0 then begin
      print,'Reading ... ',info_tmp.name
      ncdf_varget,fid,varid[vn],data
      gems = create_struct('Pressure', data)
    endif else begin
      print,'Reading ... ',info_tmp.name
      ncdf_varget,fid,varid[vn],data
      gems = create_struct(gems, 'Pressure', data)
    ENDELSE
  endif else if varname eq 'SurfaceAlbedo' $
    and n_elements(where(varname eq varlist, /null)) ne 0 then begin
    if n_elements(gems) eq 0 then begin
      print,'Reading ... ',info_tmp.name
      ncdf_varget,fid,varid[vn],data
      gems = create_struct('SurfaceAlbedo', data)
    endif else begin
      print,'Reading ... ',info_tmp.name
      ncdf_varget,fid,varid[vn],data
      gems = create_struct(gems, 'SurfaceAlbedo', data)
    ENDELSE
  endif else if varname eq 'Temperature' $
    and n_elements(where(varname eq varlist, /null)) ne 0 then begin
    if n_elements(gems) eq 0 then begin
      print,'Reading ... ',info_tmp.name
      ncdf_varget,fid,varid[vn],data
      gems = create_struct('Temperature', data)
    endif else begin
      print,'Reading ... ',info_tmp.name
      ncdf_varget,fid,varid[vn],data
      gems = create_struct(gems, 'Temperature', data)
    ENDELSE
  endif

ENDFOR  ; variable

ncdf_close, fid
return, gems

end
