function read_gems_l2_o3p, filename
nl=24
nres=20
data_list = ['EffectiveCloudFractionUV', 'ProcessingQualityFlags', $
             'RootMeanSquareErrorOfFit', 'AveragingKernel', 'O3' , $
             'O3Apriori', 'O3AprioriError','ColumnAmountO3','CloudPressure','ResidualsOfFit','NumberOfIterations','Runtime']

geo_list  = ['Latitude' ,'Longitude', 'SolarZenithAngle', $
             'ViewingZenithAngle', 'Time','Altitude' ,    $
             'Pressure','Pix','Line','TropopausePressure']

str = {alt:fltarr(nl+1), pres:fltarr(nl+1), avgk:fltarr(nl, nl), $
       o3:fltarr(nl), ao3:fltarr(nl),runtime:0.0,$
       ao3e:fltarr(nl), co3:fltarr(3),ctp:0.0,time:0.0d, $ ;wasp
       lat:0.0, lon:0.0, sza:0.0, vza:0.0, cfrac:0.0, exval:0, rms:0.0,$
       line:0, pix:0,tp:0.0,res:fltarr(nres),solerr:fltarr(nl),randerr:fltarr(nl),dfs:fltarr(3),numiter:0,finalqf:0}


fid=ncdf_open(file)
p1id=ncdf_groupsinq(fid)
p2id=ncdf_groupsinq(p1id[0])
p3id=ncdf_groupsinq(p2id[0])
p4id=ncdf_groupsinq(p3id[0])

tmp_id  = ncdf_varid(p4id[0],'Latitude')
ncdf_varget,p4id[0],tmp_id,lat
datadims= size(lat,/dim)
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
      IF varname eq 'Latitude'           then gems.lat  = data
      IF varname eq 'Longitude'           then gems.lon  = data
      IF varname eq 'Line'                then gems.line = data
      IF varname eq 'Pix'                 then gems.pix  = data
      IF varname eq 'SolarZenithAngle'    then gems.sza  = data
      IF varname eq 'ViewingZenithAngle'  then gems.vza  = data
      IF varname eq 'Pressure'            then BEGIN
        FOR il = 0, nl DO BEGIN
          tmpdata = data[*,*,il]
          gems.pres[il] = tmpdata
        ENDFOR
      ENDIF
      IF varname eq 'TropopausePressure'  then gems.tp   = data
      ; should be hpa
      IF varname eq 'Altitude'            THEN BEGIN
        for il = 0, nl do begin
          tmpdata = data[*,*,il]
          gems.alt[il] = tmpdata
        endfor
      ENDIF
      ; should be km
      IF varname eq 'Time'                then begin
         for j = 0 , npix -1 do begin
           gems[j,*].time = transpose(data)
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
          gems.o3[il] = tmpdata
        ENDFOR
      ENDIF
      IF varname eq 'ColumnAmountO3'           then gems.co3    = data
      IF varname eq 'O3Apriori'                THEN BEGIN
        FOR il = 0, nl-1 DO BEGIN
          tmpdata = data[*,*,il]
          gems.ao3[il] = tmpdata
        ENDFOR
      ENDIF
      IF varname eq 'O3AprioriError'           then BEGIN
        FOR il = 0, nl-1 DO BEGIN
          tmpdata = data[*,*,il]
          gems.ao3e[il] = tmpdata
        ENDFOR
      ENDIF
      IF varname eq 'O3RandomNoiseError'       then BEGIN
        FOR il = 0, nl-1 DO BEGIN
          tmpdata = data[*,*,il]
          gems.randerr[il] = tmpdata
        ENDFOR
      ENDIF
      IF varname eq 'O3SolutionError'          then BEGIN
        FOR il = 0, nl-1 DO BEGIN
          tmpdata = data[*,*,il]
          gems.solerr[il] = tmpdata
        ENDFOR
      ENDIF
      IF varname eq 'RootMeanSquareErrorOfFit' then gems.rms    = data
      IF varname eq 'DegreesOfFreedomForSignal'then gems.dfs    = data
      IF varname eq 'AveragingKernel'          then BEGIN
        FOR il1 = 0, nl-1 DO BEGIN
        FOR il2 = 0, nl-1 DO BEGIN
          tmpdata = data[*,*,il1,il2]
          gems.avgk[il1,il2] = tmpdata
        ENDFOR
        ENDFOR
      ENDIF
      IF varname eq 'CloudPressure'            then gems.ctp    = data
      IF varname eq 'ResidualsOfFit'           THEN gems.res    = data
      IF varname eq 'EffectiveCloudFractionUV' then gems.cfrac  = data
      IF varname eq 'ProcessingQualityFlags'   then gems.exval  = data
      IF varname eq 'FianlAlgorithmFlags'      then gems.finalqf = data
      IF varname eq 'NumberOfIterations'       then gems.numiter = data
      IF varname eq 'Runtime'                  then gems.runtime = data
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
;save,file=xdrfn,gems,/xdr
endif else restore,xdrfn
return gems

end
