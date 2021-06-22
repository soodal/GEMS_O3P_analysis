;====================================================================================;
; This function is to read GEOSChem NetCDF file from SNU. 
;
; PROGRAM NAME:
;  read_geoschem_gas_nc.pro
;
; PROCESS:
; 
; REFERENCE:
;
; INPUT:
;  NetCDF file path and Array of parameter name
;
; OUTPUT: 
;
; DEVELOPER:
;  MI JEONG KIM
;  Satellite Remote Sensing Laboratory
;  Division of Earth Environmental System
;  College of Natural Science
;  Pusan National University
;                  Tel: +82-51-510-2172
;                  E-mail: mjkellykim@pusan.ac.kr
;
; REVISION HISTORY:
;  ver1. 2019/12/12
;
;====================================================================================;
FUNCTION read_geoschem_gas_nc, file=file, param=param, onlyp=onlyp

;;=============================================================================
;; 1. NCDF File Open and Get Information
;;=============================================================================
  nid=NCDF_OPEN(file)

  file_struct=NCDF_INQUIRE(nid)
  nvars=file_struct.nvars
  
  FOR i=0,nvars-1 DO BEGIN
    var_struct=NCDF_VARINQ(nid,i)
    var_name=var_struct.name
    var_ndim=var_struct.ndims
    IF var_ndim EQ 1 THEN BEGIN
      dum=EXECUTE('NCDF_VARGET,nid,var_name,'+var_name)
    ENDIF 
  ENDFOR

  ;; Read Info
  nx = N_ELEMENTS(lon)
  ny = N_ELEMENTS(lat)
  nl = N_ELEMENTS(lev)
  nt = N_ELEMENTS(time)
  np = N_ELEMENTS(param)

  ;; Time 
  CASE STRMID(file,6,3,/REVERSE) OF
    'Jan' : mon=1
    'Apr' : mon=4
    'Jul' : mon=7
    'Oct' : mon=10
  ENDCASE

  j0=JULDAY(mon,1,2016,0,0)
  time=time/24D/60D + j0
  CALDAT,time $
        ,month,day,year,hour,minute
  time=STRING(year,F='(i04)')+STRING(month,F='(i02)')$
      +STRING(day, F='(i02)')+STRING(hour, F='(i02)')+STRING(minute, F='(i02)')


  ;; Lon/Lat
  lon = REBIN(          lon ,nx,ny,/SAMPLE)
  lat = REBIN(TRANSPOSE(lat),nx,ny,/SAMPLE)  

  ;; Deallocate!
  j0=!null
  year=!null
  month=!null
  day=!null
  hour=!null

;;=============================================================================
;; 2. Get Parameters and Close
;;=============================================================================
  units=STRARR(np)
  ;; true_data = add_offset + (vlaue)*scale_factor
  FOR ip=0,np-1 DO BEGIN
    var_name=param[ip]
    NCDF_VARGET,nid,var_name,var
      NCDF_ATTGET,nid,var_name,'_FillValue',fill
      NCDF_ATTGET,nid,var_name,'units',dumunit
        bad=WHERE(var EQ fill, nbad)
        IF(nbad GT 0)THEN var[bad]=!values.f_nan
      units[ip]=STRING(dumunit)
      dum=EXECUTE(param[ip]+'=FLOAT(var)')
      var=!null
  ENDFOR
 
  ;; Close
  NCDF_CLOSE,nid

;;=============================================================================
;; 3. Return Output
;;=============================================================================

  outdat=FLTARR(nx,ny,nl,nt,np)
  FOR ip=0,np-1 DO BEGIN
    dum=EXECUTE('outdat[*,*,*,*,ip]='+param[ip])
  ENDFOR
  outdat = REFORM(outdat)
 
  IF ~KEYWORD_SET(onlyp) THEN BEGIN
    output = CREATE_STRUCT('lon', lon, 'lat', lat, 'lev', lev, $
                           'time', time, 'dat', outdat, 'unit', units)
  ENDIF ELSE BEGIN
    output = CREATE_STRUCT('dat', outdat, 'unit', units)
  ENDELSE 

  RETURN, output
END
