function ds_read_gems_l2_o3t, file

fid=ncdf_open(file)

p1id=ncdf_groupsinq(fid)
p1name=ncdf_groupname(p1id[0])

if p1name eq 'NETCDF4' then begin
  p2id=ncdf_groupsinq(p1id[0])
  p3id=ncdf_groupsinq(p2id[0])
  p4id=ncdf_groupsinq(p3id[0])

  lat_id  = ncdf_varid(p4id[0],'Latitude')
  res_id  = ncdf_varid(p4id[1],'ResidualsOfFit')
  ;o3solerr_id = ncdf_varid(p4id[1], 'O3SolutionError')

  ;ncdf_varget,p4id[0],lat_id,Latitude
  ncdf_varget,p4id[1],res_id,ResidualsOfFit
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
image=695



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

struct = { $
    AveragingKernel:fltarr(nlayer2, nlayer2), $
    CloudPressure:0.0, $
    ColumnAmountO3:0.0, $
    DegreesOfFreedomForSignal:0.0, $
    dNdR:fltarr(2), $
    dR_dl:0.0, $
    EstimatedError:0.0, $
    FinalAlgorithmflags:0.0, $
    LayerEfficiency:fltarr(nlayer2), $
    Nvalue:fltarr(nwl), $
    O3BelowCloud:0.0, $
    Reflectivity340:0.0, $
    Reflectivity380:0.0, $
    Residue:fltarr(nsig), $
    StepOneO3:0.0, $
    StepTwoO3:0.0, $
    TerrainPressure:0.0, $
    Latitude:0.0, $
    Longitude:0.0, $
    RelativeAzimuthAngle:0.0, $
    SolarZenithAngle:0.0, $
    ViewingZenithAngle:0.0}

gems=replicate(struct,image,spatial)


for dn=0,1 do begin
  id=p4id[dn]
  varid=ncdf_varidsinq(id)
  datpath=ncdf_fullgroupname(p4id[1])
  geopath=ncdf_fullgroupname(p4id[0])

  
  if dn eq 0 THEN begin ; Data Fields
    FOR vn=0,n_elements(varid)-1 DO BEGIN
      info_tmp=ncdf_varinq(id,varid[vn])
      print,'Reading ... ',info_tmp.name
      varname=info_tmp.name
      ncdf_varget,id,varid[vn],data

      CASE varname of
        'AveragingKernel': BEGIN
          FOR il1 = 0, nlayer2-1 DO BEGIN
            FOR il2 = 0, nlayer2-1 DO BEGIN
              tmpdata = reform(data[*,*,il1,il2])
              gems.AveragingKernel[il1,il2] = tmpdata
            ENDFOR
          ENDFOR
          ;gems.AveragingKernel = data
          END
        'CloudPressure': gems.CloudPressure = data
        'ColumnAmountO3': gems.ColumnAmountO3 = data
        'DegreesOfFreedomForSignal': BEGIN
          sz = size(data)
          if sz[0] eq 3 then begin
            data = reform(data[0, *, *])
          ENDIF
          gems.DegreesOfFreedomForSignal = data
        END
        'dNdR': gems.dNdR = data
        'dR_dl': gems.dR_dl = data
        'EstimatedError': gems.EstimatedError = data
        'FinalAlgorithmFlags': gems.FinalAlgorithmFlags = data
        'LayerEfficiency': BEGIN
          for il = 0, nlayer2-1 do begin
            tmpdata = data[*, *, il]
            gems.LayerEfficiency[il] = tmpdata
          endfor
        END
        'Nvalue': BEGIN
          FOR il = 0, nwl-1 DO BEGIN
            tmpdata = reform(data[il, *,*])
            gems.Nvalue[il] = tmpdata
          ENDFOR
          ;gems.Nvalue = data
          END
        'O3BelowCloud': BEGIN
          ;FOR il = 0, nlayer2-1 DO BEGIN
            ;tmpdata = data[*,*,il]
            ;gems.O3BelowCloud = tmpdata
          ;ENDFOR
            gems.O3BelowCloud = data
          END
        'Reflectivity340': BEGIN
          ;FOR il = 0, nlayer2-1 DO BEGIN
            ;tmpdata = reform(data[*,*,il])
            ;gems.Reflectivity340[il] = tmpdata
          ;ENDFOR
          gems.Reflectivity340 = data
          END
        'Reflectivity380': BEGIN
          ;FOR il = 0, nlayer2-1 DO BEGIN
            ;tmpdata = reform(data[*,*,il])
            ;gems.Reflectivity380[il] = tmpdata
          ;ENDFOR
          gems.Reflectivity380 = data
          END
        'Residue': BEGIN
          ;FOR il = 0, nsig-1 DO BEGIN
            ;tmpdata = reform(data[*,*,il])
            ;gems.Residue[il] = tmpdata
          ;ENDFOR
            gems.Residue = data

          END
        'StepOneO3': gems.StepOneO3 = data
        'StepTwoO3': gems.StepTwoO3 = data
        'TerrainPressure': gems.TerrainPressure = data
      ENDCASE

    ENDFOR
  ENDIF ELSE BEGIN ; GeolocationFields
    for vn=0,n_elements(varid)-1 do begin
      info_tmp=ncdf_varinq(id,varid[vn])
      print,'Reading ... ',info_tmp.name
      varname=info_tmp.name
      ncdf_varget,id,varid[vn],data

      CASE varname of
        'Longitude': gems.Longitude = data
        'Latitude': gems.Latitude = data
        'RelativeAzimuthAngle': gems.RelativeAzimuthAngle = data
        'SolarZenithAngle': gems.SolarZenithAngle  = data
        'ViewingZenithAngle': gems.ViewingZenithAngle = data
      ENDCASE


    ENDFOR  ; variable
  ENDELSE
ENDFOR  ; group

ncdf_close, fid
return, gems

end
