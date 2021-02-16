function ds_read_gems_l2_o3p, file, varlist=varlist
varlist2get = []
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

  ;lat_id  = ncdf_varid(p4id[0],'Latitude')
  ;res_id  = ncdf_varid(p4id[1],'ResidualsOfFit')
  ;o3solerr_id = ncdf_varid(p4id[1], 'O3SolutionError')

  ;ncdf_varget,p4id[0],lat_id,Latitude
  ;ncdf_varget,p4id[1],res_id,ResidualsOfFit
  ;ncdf_varget,p4id[1],o3solerr_id,O3SolutionError
ENDIF else begin
  p2id=ncdf_groupsinq(p1id[0])
  ;p3id=ncdf_groupsinq(p2id[0])
  ;p4id=ncdf_groupsinq(p3id[0])
  p4id=p1id

  ;lat_id  = ncdf_varid(p4id[0],'Latitude')
  ;res_id  = ncdf_varid(p4id[1],'ResidualsOfFit')
  ;o3solerr_id = ncdf_varid(p4id[1], 'O3SolutionError')

  ;ncdf_varget,p4id[0],lat_id,Latitude
  ;ncdf_varget,p4id[1],res_id,ResidualsOfFit
  ;ncdf_varget,p4id[1],o3solerr_id,O3SolutionError
ENDELSE



;data_list = ['EffectiveCloudFractionUV', 'ProcessingQualityFlags', $
             ;'RootMeanSquareErrorOfFit', 'AveragingKernel', 'O3' , $
             ;'O3Apriori', 'O3AprioriError','ColumnAmountO3', $
             ;'CloudPressure','ResidualsOfFit','NumberOfIterations','Runtime',
             ;'SimulatedRadiances'] ; Added by Dae Sung 2020-08-29

;geo_list  = ['Latitude' ,'Longitude', 'SolarZenithAngle', $
             ;'ViewingZenithAngle', 'Time','Altitude' ,    $
             ;'Pressure','Pix','Line','TropopausePressure']

;struct = {Altitude:fltarr(nl+1), Pressure:fltarr(nl+1), $
  ;AveragingKernel:fltarr(nl, nl), O3:fltarr(nl), O3Apriori:fltarr(nl), $
  ;Runtime:0.0, O3AprioriError:fltarr(nl), ColumnAmountO3:fltarr(3),CloudPressure:0.0, $ 
  ;Time:0.0d, Latitude:0.0, Longitude:0.0, SolarZenithAngle:0.0, $
  ;Temperature:fltarr(nl+1), $
  ;TerrainPressure:0.0, $
  ;RelativeAzimuthAngle: 0.0, $
  ;TerrainReflectivityUV: 0.0, $
  ;Exitval: 0.0, $
  ;ViewingZenithAngle:0.0, EffectiveCloudFractionUV:0.0, ProcessingQualityFlags:0, $
  ;RootMeanSquareErrorOfFit:0.0, Line:0, Pix:0,TropopausePressure:0.0, $
  ;ResidualsOfFit:fltarr(nres), SimulatedRadiances:fltarr(nres), O3SolutionError:fltarr(nl), $
  ;O3RandomNoiseError:fltarr(nl), DegreesOfFreedomForSignal:!values.d_nan, $
  ;NumberOfIterations:0, FinalAlgorithmFlags:0}

;gems=replicate(struct,npix,nline)



for dn=0,1 do begin
  id=p4id[dn]
  varid=ncdf_varidsinq(id)
  geopath=ncdf_fullgroupname(p4id[0])
  datpath=ncdf_fullgroupname(p4id[1])

  
  if dn eq 0 THEN begin ; geopath
    for vn=0,n_elements(varid)-1 do begin
      info_tmp=ncdf_varinq(id,varid[vn])
      varname=info_tmp.name

      if varname eq 'Altitude' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then begin
        if n_elements(gems) eq 0 then begin
            print,'Reading ... ',info_tmp.name
            ncdf_varget,id,varid[vn],data
            gems = create_struct('Altitude', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'Altitude', data)
        ENDELSE
      endif else if varname eq 'Latitude' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then begin
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          nanidx = where(data lt -990, /null)
          data[nanidx] = !values.f_nan
          gems = create_struct('Latitude', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          nanidx = where(data lt -990, /null)
          data[nanidx] = !values.f_nan
          gems = create_struct(gems, 'Latitude', data)
        ENDELSE
      endif else if varname eq 'Line' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then begin
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('Line', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'Line', data)
        ENDELSE
      endif else if varname eq 'Longitude' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then begin
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          nanidx = where(data lt -990, /null)
          data[nanidx] = !values.f_nan
          gems = create_struct('Longitude', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          nanidx = where(data lt -990, /null)
          data[nanidx] = !values.f_nan
          gems = create_struct(gems, 'Longitude', data)
        ENDELSE
      endif else if varname eq 'Pix' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then begin
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('Pix', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'Pix', data)
        ENDELSE
      endif else if varname eq 'Pressure' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then begin
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('Pressure', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'Pressure', data)
        ENDELSE
      endif else if varname eq 'RelativeAzimuthAngle' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then begin
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
        and n_elements(where(varname eq varlist, /null)) ne 0 then begin
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('SolarZenithAngle', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'SolarZenithAngle', data)
        ENDELSE
      endif else if varname eq 'Temperature' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then begin
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('Temperature', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'Temperature', data)
        ENDELSE
      endif else if varname eq 'TerrainPressure' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then begin 
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('TerrainPressure', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'TerrainPressure', data)
        ENDELSE
      endif else if varname eq 'Time' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then begin
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('Time', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'Time', data)
        ENDELSE
       endif else if varname eq 'TropopausePressure' $
         and n_elements(where(varname eq varlist, /null)) ne 0 then begin
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('TropopausePressure', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'TropopausePressure', data)
        ENDELSE
      endif else if varname eq 'ViewingZenithAngle' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then begin
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('ViewingZenithAngle', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'ViewingZenithAngle', data)
        ENDELSE
      endif

    ENDFOR  ; variable
  ENDIF ELSE BEGIN ; datpath 
    FOR vn=0,n_elements(varid)-1 DO BEGIN
      info_tmp=ncdf_varinq(id,varid[vn])
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
      ENDif else if varname eq  'EffectiveCloudFractionUV' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then BEGIN
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('EffectiveCloudFractionUV', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'EffectiveCloudFractionUV', data)
        ENDELSE
      endif else if varname eq 'Exitval' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then BEGIN
        ;gems.Exitval = data
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('Exitval', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'Exitval', data)
        ENDELSE
      ENDif else if varname eq 'FinalAlgorithmFlags' $
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
      endif else if varname eq  'NumberOfIterations' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then begin
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('NumberOfIterations', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'NumberOfIterations', data)
        ENDELSE
      endif else if varname eq 'O3' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then begin
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('O3', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'O3', data)
        ENDELSE
      ENDif else if varname eq 'O3Apriori' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then BEGIN
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('O3Apriori', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'O3Apriori', data)
        ENDELSE
      ENDif else if varname eq 'O3AprioriError' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then BEGIN
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('O3AprioriError', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'O3AprioriError', data)
        ENDELSE
      ENDif else if varname eq 'O3RandomNoiseError' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then BEGIN
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('O3RandomNoiseError', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'O3RandomNoiseError', data)
        ENDELSE
      ENDif else if varname eq 'O3SolutionError' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then BEGIN
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('O3SolutionError', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'O3SolutionError', data)
        ENDELSE
      ENDif else if varname eq  'ProcessingQualityFlags' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then begin
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('ProcessingQualityFlags', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'ProcessingQualityFlags', data)
        ENDELSE
      endif else if varname eq 'ResidualsOfFit' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then begin
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('ResidualsOfFit', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'ResidualsOfFit', data)
        ENDELSE
      endif else if varname eq 'RootMeanSquareErrorOfFit' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then begin
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('RootMeanSquareErrorOfFit', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'RootMeanSquareErrorOfFit', data)
        ENDELSE
      endif else if varname eq 'Runtime' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then begin 
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('Runtime', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'Runtime', data)
        ENDELSE
      endif else if varname eq 'TerrainReflectivityUV' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then begin
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('TerrainReflectivityUV', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'TerrainReflectivityUV', data)
        ENDELSE
      endif else if varname eq 'SimulatedRadiances' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then begin 
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('SimulatedRadiances', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'SimulatedRadiances', data)
        ENDELSE
      endif else if varname eq 'FitWeights' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then begin 
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('Fitweights', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'FitWeights', data)
        ENDELSE
      endif else if varname eq 'SignalToNoiseRatio' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then begin 
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('SignalToNoiseRatio', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'SignalToNoiseRabnewtio', data)
        ENDELSE
      ENDIF else if varname eq 'Wavelengths' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then begin
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('Wavelengths', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'Wavelengths', data)
        ENDELSE
      endif else if varname eq 'WavelengthsWholeRange' $
        and n_elements(where(varname eq varlist, /null)) ne 0 then begin
        if n_elements(gems) eq 0 then begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct('WavelengthsWholeRange', data)
        endif else begin
          print,'Reading ... ',info_tmp.name
          ncdf_varget,id,varid[vn],data
          gems = create_struct(gems, 'WavelengthsWholeRange', data)
        ENDELSE
      endif


    ENDFOR
  ENDELSE
ENDFOR  ; group

ncdf_close, fid
return, gems

end
