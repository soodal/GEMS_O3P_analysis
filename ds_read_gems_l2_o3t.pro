function ds_read_gems_l2_o3t, file, varlist=varlist

if not keyword_set(varlist) then begin
  varlist = ['EffectiveCloudFractionUV', 'ProcessingQualityFlags', $
             'RootMeanSquareErrorOfFit', 'AveragingKernel', 'O3' , $
             'O3Apriori', 'O3AprioriError','ColumnAmountO3', $
             'CloudPressure','ResidualsOfFit','NumberOfIterations','Runtime', $
             'SimulatedRadiances', 'Latitude' ,'Longitude', 'SolarZenithAngle', $
             'ViewingZenithAngle', 'Time','Altitude' ,    $
             'Pressure','Pix','Line','TropopausePressure', 'Temperature', $
             'TerrainPressure', 'RelativeAzimuthAngle', $
             'FitWeights', 'SignalToNoiseRatio', 'Wavelengths', $
             'WavelengthsWholeRange']
endif
fid=ncdf_open(file)

p1id=ncdf_groupsinq(fid)
p1name=ncdf_groupname(p1id[0])

if p1name eq 'NETCDF4' then begin
  p2id=ncdf_groupsinq(p1id[0])
  p3id=ncdf_groupsinq(p2id[0])
  p4id=ncdf_groupsinq(p3id[0])

  lat_id  = ncdf_varid(p4id[0],'Latitude')
  ;res_id  = ncdf_varid(p4id[1],'ResidualsOfFit')
  ;o3solerr_id = ncdf_varid(p4id[1], 'O3SolutionError')

  ;ncdf_varget,p4id[0],lat_id,Latitude
  ;ncdf_varget,p4id[1],res_id,ResidualsOfFit
  ;ncdf_varget,p4id[1],o3solerr_id,O3SolutionError
ENDIF else begin ; Data Fields or Geolocation Fields
  p4id=p1id

  lat_id  = ncdf_varid(p4id[1],'Latitude')
  lon_id = ncdf_varid(p4id[1],'Longitude')
  raa_id = ncdf_varid(p4id[1],'RelativeAzimuthAngle')
  sza_id = ncdf_varid(p4id[1],'SolarZenithAngle')
  vza_id = ncdf_varid(p4id[1],'ViewingZenithAngle')

ENDELSE

;latdims = size(Latitude, /dim)
resdims = size(ResidualsOfFit,/dim)
;o3solerrdims = size(O3SolutionError,/dim)
nlayer2=11
nwl = 7
nsig = 3
spatial=2048
image=701



;nres = resdims[0]
;npix = resdims[1] 
;nline = resdims[2]


;data_list = ['EffectiveCloudFractionUV', 'ProcessingQualityFlags', $
             ;'RootMeanSquareErrorOfFit', 'AveragingKernel', 'O3' , $
             ;'O3Apriori', 'O3AprioriError','ColumnAmountO3', $
             ;'CloudPressure','ResidualsOfFit','NumberOfIterations','Runtime',
             ;'SimulatedRadiances'] ; Added by Dae Sung 2020-08-29

;geo_list  = ['Latitude' ,'Longitude', 'SolarZenithAngle', $
             ;'ViewingZenithAngle', 'Time','Altitude' ,    $
             ;'Pressure','Pix','Line','TropopausePressure']

;struct = { $
    ;AveragingKernel:fltarr(nlayer2, nlayer2), $
    ;CloudPressure:0.0, $
    ;ColumnAmountO3:0.0, $
    ;DegreesOfFreedomForSignal:0.0, $
    ;dNdR:fltarr(2), $
    ;dR_dl:0.0, $
    ;EstimatedError:0.0, $
    ;FinalAlgorithmflags:0.0, $
    ;LayerEfficiency:fltarr(nlayer2), $
    ;Nvalue:fltarr(nwl), $
    ;O3BelowCloud:0.0, $
    ;Reflectivity340:0.0, $
    ;Reflectivity380:0.0, $
    ;Residue:fltarr(nsig), $
    ;StepOneO3:0.0, $
    ;StepTwoO3:0.0, $
    ;TerrainPressure:0.0, $
    ;Latitude:0.0, $
    ;Longitude:0.0, $
    ;RelativeAzimuthAngle:0.0, $
    ;SolarZenithAngle:0.0, $
    ;ViewingZenithAngle:0.0}

;gems=replicate(struct,image,spatial)


for dn=0,1 do begin
  id=p4id[dn]
  varid=ncdf_varidsinq(id)
  ;datpath=ncdf_fullgroupname(p4id[1])
  ;geopath=ncdf_fullgroupname(p4id[0])

  
  if dn eq 1 THEN begin ; Data Fields
    FOR vn=0,n_elements(varid)-1 DO BEGIN
      info_tmp=ncdf_varinq(id,varid[vn])
      print,'Reading ... ',info_tmp.name
      varname=info_tmp.name
      ncdf_varget,id,varid[vn],data

      if varname eq 'AveragingKernel' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then BEGIN
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('AveragingKernel', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'AveragingKernel', data)
        ENDELSE
      ENDIF else if varname eq 'CloudPressure' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then BEGIN
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('CloudPressure', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'CloudPressure', data)
        ENDELSE
      ENDIF else if varname eq 'ColumnAmountO3' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then BEGIN
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('ColumnAmountO3', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'ColumnAmountO3', data)
        ENDELSE
      endif else if varname eq 'DegreesOfFreedomForSignal' $
          and n_elements(where(varname eq varlist, /null)) ne 0 then BEGIN
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('DegreesOfFreedomForSignal', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'DegreesOfFreedomForSignal', data)
        ENDELSE
      ENDIF else if varname eq 'dNdR' $
          and n_elements(where(varname eq varlist, /null)) ne 0 then BEGIN
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('dNdR', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'dNdR', data)
        ENDELSE
      endif else if varname eq 'dR_dl' $
          and n_elements(where(varname eq varlist, /null)) ne 0 then BEGIN
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('dR_dl', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'dR_dl', data)
        ENDELSE
      endif else if varname eq 'EstimatedError' $
          and n_elements(where(varname eq varlist, /null)) ne 0 then BEGIN
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('EstimatedError', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'EstimatedError', data)
        ENDELSE
      endif else if varname eq 'FinalAlgorithmFlags' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then BEGIN
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('FinalAlgorithmFlags', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'FinalAlgorithmFlags', data)
        ENDELSE
      endif else if varname eq 'LayerEfficiency' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then BEGIN
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('LayerEfficiency', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'LayerEfficiency', data)
        ENDELSE
      ENDIF else if varname eq 'Nvalue' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then BEGIN
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('Nvalue', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'Nvalue', data)
        ENDELSE
      ENDif else if varname eq 'O3BelowCloud' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then BEGIN
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('O3BelowCloud', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'O3BelowCloud', data)
        ENDELSE
      ENDif else if varname eq 'Reflectivity340' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then BEGIN
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('Reflectivity340', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'Reflectivity340', data)
        ENDELSE
      ENDif else if varname eq 'Reflectivity380' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then BEGIN
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('Reflectivity380', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'Reflectivity380', data)
        ENDELSE
      ENDif else if varname eq 'Residue' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then BEGIN
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('Residue', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'Residue', data)
        ENDELSE
      ENDif else if varname eq 'StepOneO3' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then BEGIN
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('StepOneO3', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'StepOneO3', data)
        ENDELSE
      endif else if varname eq 'StepTwoO3' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then BEGIN
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('StepTwoO3', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'StepTwoO3', data)
        ENDELSE
      endif else if varname eq 'TerrainPressure' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then BEGIN
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('TerrainPressure', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'TerrainPressure', data)
        ENDELSE
      endif
    ENDFOR
  ENDIF ELSE BEGIN ; GeolocationFields
    for vn=0,n_elements(varid)-1 do begin
      info_tmp=ncdf_varinq(id,varid[vn])
      print,'Reading ... ',info_tmp.name
      varname=info_tmp.name
      ncdf_varget,id,varid[vn],data

      if varname eq 'Longitude' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then BEGIN
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('Longitude', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'Longitude', data)
        ENDELSE
      endif else if varname eq 'Latitude' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then BEGIN
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('Latitude', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'Latitude', data)
        ENDELSE
      endif else if varname eq 'RelativeAzimuthAngle' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then BEGIN
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('RelativeAzimuthAngle', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'RelativeAzimuthAngle', data)
        ENDELSE
      endif else if varname eq 'SolarZenithAngle' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then BEGIN
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('SolarZenithAngle', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'SolarZenithAngle', data)
        ENDELSE
      endif else if varname eq 'ViewingZenithAngle' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then BEGIN
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('ViewingZenithAngle', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'ViewingZenithAngle', data)
        ENDELSE
      ENDIF

    ENDFOR  ; variable
  ENDELSE
ENDFOR  ; group

ncdf_close, fid
return, gems

end
