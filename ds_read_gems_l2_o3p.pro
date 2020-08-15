function ds_read_gems_l2_o3p, file
nl=24
nres=20
data_list = ['EffectiveCloudFractionUV', 'ProcessingQualityFlags', $
             'RootMeanSquareErrorOfFit', 'AveragingKernel', 'O3' , $
             'O3Apriori', 'O3AprioriError','ColumnAmountO3', $
             'CloudPressure','ResidualsOfFit','NumberOfIterations','Runtime']

geo_list  = ['Latitude' ,'Longitude', 'SolarZenithAngle', $
             'ViewingZenithAngle', 'Time','Altitude' ,    $
             'Pressure','Pix','Line','TropopausePressure']

str = {Altitude:fltarr(nl+1), Pressure:fltarr(nl+1), $
  AveragingKernel:fltarr(nl, nl), O3:fltarr(nl), O3Apriori:fltarr(nl), $
  Runtime:0.0, O3AprioriError:fltarr(nl), ColumnAmountO3:fltarr(3),CloudPressure:0.0, $ 
  Time:0.0d, Latitude:0.0, Longitude:0.0, SolarZenithAngle:0.0, $
  ViewingZenithAngle:0.0, EffectiveCloudFractionUV:0.0, ProcessingQualityFlags:0, $
  RootMeanSquareErrorOfFit:0.0, Line:0, Pix:0,TropopausePressure:0.0, $
  ResidualsOfFit:fltarr(nres), O3SolutionError:fltarr(nl), $
  O3RandomNoiseError:fltarr(nl), DegreesOfFreedomForSignal:fltarr(3), $
  NumberOfIterations:0, FinalAlgorithmFlags:0}


fid=ncdf_open(file)
p1id=ncdf_groupsinq(fid)
p2id=ncdf_groupsinq(p1id[0])
p3id=ncdf_groupsinq(p2id[0])
p4id=ncdf_groupsinq(p3id[0])

tmp_id  = ncdf_varid(p4id[0],'Latitude')
ncdf_varget,p4id[0],tmp_id,Latitude
datadims= size(Latitude,/dim)
npix=datadims[0] & nline=datadims[1]
gems=replicate(str,npix,nline)

tmpdata = fltarr(npix,nline)

for dn=0,1 do begin
  id=p4id[dn]
  varid=ncdf_varidsinq(id)
  datpath=ncdf_fullgroupname(p4id[1])
  geopath=ncdf_fullgroupname(p4id[0])

  
  if dn eq 0 then begin ; geopath
    for vn=0,n_elements(varid)-1 do begin
      info_tmp=ncdf_varinq(id,varid[vn])
      print,'Reading ... ',info_tmp.name
      varname=info_tmp.name
      ncdf_varget,id,varid[vn],data
      IF varname eq 'Latitude'           then gems.Latitude  = data
      IF varname eq 'Longitude'           then gems.Longitude = data
      IF varname eq 'Line'                then gems.Line = data
      IF varname eq 'Pix'                 then gems.Pix  = data
      IF varname eq 'SolarZenithAngle'    then gems.SolarZenithAngle  = data
      IF varname eq 'ViewingZenithAngle'  then gems.ViewingZenithAngle  = data
      IF varname eq 'Pressure'            then BEGIN
        FOR il = 0, nl DO BEGIN
          tmpdata = data[*,*,il]
          gems.Pressure[il] = tmpdata
        ENDFOR
      ENDIF
      IF varname eq 'TropopausePressure'  then gems.TropopausePressure   = data
      ; should be hpa
      IF varname eq 'Altitude'            THEN BEGIN
        for il = 0, nl do begin
          tmpdata = data[*,*,il]
          gems.Altitude[il] = tmpdata
        endfor
      ENDIF
      ; should be km
      IF varname eq 'Time'                then begin
         for j = 0 , npix -1 do begin
           gems[j,*].Time = transpose(data)
         endfor
      endif

    ENDFOR  ; variable
  ENDIF ELSE BEGIN ; datpath 
    FOR vn=0,n_elements(varid)-1 DO BEGIN
      info_tmp=ncdf_varinq(id,varid[vn])
      print,'Reading ... ',info_tmp.name
      varname=info_tmp.name
      ncdf_varget,id,varid[vn],data
      IF varname eq 'O3'                       THEN BEGIN
        FOR il = 0, nl-1 DO BEGIN
          tmpdata = data[*,*,il]
          gems.O3[il] = tmpdata
        ENDFOR
      ENDIF
      IF varname eq 'ColumnAmountO3'           then gems.ColumnAmountO3    = data
      IF varname eq 'O3Apriori'                THEN BEGIN
        FOR il = 0, nl-1 DO BEGIN
          tmpdata = data[*,*,il]
          gems.O3Apriori[il] = tmpdata
        ENDFOR
      ENDIF
      IF varname eq 'O3AprioriError'           then BEGIN
        FOR il = 0, nl-1 DO BEGIN
          tmpdata = data[*,*,il]
          gems.O3AprioriError[il] = tmpdata
        ENDFOR
      ENDIF
      IF varname eq 'O3RandomNoiseError'       then BEGIN
        FOR il = 0, nl-1 DO BEGIN
          tmpdata = data[*,*,il]
          gems.O3RandomNoiseError[il] = tmpdata
        ENDFOR
      ENDIF
      IF varname eq 'O3SolutionError'          then BEGIN
        FOR il = 0, nl-1 DO BEGIN
          tmpdata = data[*,*,il]
          gems.O3SolutionError[il] = tmpdata
        ENDFOR
      ENDIF
      IF varname eq 'RootMeanSquareErrorOfFit' then gems.RootMeanSquareErrorOfFit = data
      IF varname eq 'DegreesOfFreedomForSignal'then gems.DegreesOfFreedomForSignal= data
      IF varname eq 'AveragingKernel'          then BEGIN
        FOR il1 = 0, nl-1 DO BEGIN
        FOR il2 = 0, nl-1 DO BEGIN
          tmpdata = data[*,*,il1,il2]
          gems.AveragingKernel[il1,il2] = tmpdata
        ENDFOR
        ENDFOR
      ENDIF
      IF varname eq 'CloudPressure'            then gems.CloudPressure    = data
      IF varname eq 'ResidualsOfFit'           THEN gems.ResidualsOfFit    = data
      IF varname eq 'EffectiveCloudFractionUV' then gems.EffectiveCloudFractionUV  = data
      IF varname eq 'ProcessingQualityFlags'   then gems.ProcessingQualityFlags  = data
      IF varname eq 'FinalAlgorithmFlags'      then gems.FinalAlgorithmFlags = data
      IF varname eq 'NumberOfIterations'       then gems.NumberOfIterations = data
      IF varname eq 'Runtime'                  then gems.Runtime = data
    ENDFOR
  ENDELSE
ENDFOR  ; group

;datid=p4id[1]
;geoid=p4id[0]

;id=datid
;if field eq 'Geolocation Fields' then id=geoid

;n=0

;varname='O3'

;for n=0,varid[-1] do begin
  ;info_tmp=ncdf_varinq(id,varid[n])
  ;;print,info_tmp.name
  ;if info_tmp.name eq varname then begin
    ;ncdf_varget,id,varid[n],data
    ;info=info_tmp
  ;endif
;endfor

ncdf_close, fid
return, gems

end
