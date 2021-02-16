pro ds_read_omi_l2_profoz, file,result
               
 dg = '/HDFEOS/SWATHS/OMI Vertical Ozone Profile/Data Fields'
 gg = '/HDFEOS/SWATHS/OMI Vertical Ozone Profile/Geolocation Fields'

  ; Geolocation Fields
  GroundPixelQualityFlags = h5read(file, gg+'/GroundPixelQualityFlags')
  nanidx = where(GroundPixelQualityFlags lt -1.e30, /null)
  GroundPixelQualityFlags[nanidx] = !values.f_nan

  Latitude = h5read(file,gg+'/Latitude')
  nanidx = where(Latitude lt -1.e30, /null)
  Latitude[nanidx] = !values.f_nan

  LatitudePixelCorner = h5read(file,gg+'/LatitudePixelCorner')
  nanidx = where(LatitudePixelCorner lt -1.e30, /null)
  LatitudePixelCorner[nanidx] = !values.f_nan

  Longitude = h5read(file,gg+'/Longitude')
  nanidx = where(Longitude lt -1.e30, /null)
  Longitude[nanidx] = !values.f_nan

  LongitudePixelCorner = h5read(file,gg+'/LongitudePixelCorner')
  nanidx = where(LongitudePixelCorner lt -1.e30, /null)
  LongitudePixelCorner[nanidx] = !values.f_nan

  RelativeAzimuthAngle = h5read(file, gg+'/RelativeAzimuthAngle')
  nanidx = where(RelativeAzimuthAngle lt -1.e30, /null)
  RelativeAzimuthAngle[nanidx] = !values.f_nan

  SecondsInDay = h5read(file,gg+'/SecondsInDay')
  nanidx = where(SecondsInDay lt -1.e30, /null)
  SecondsInDay[nanidx] = !values.f_nan

  SolarZenithAngle = h5read(file,gg+'/SolarZenithAngle')
  nanidx = where(SolarZenithAngle lt -1.e30, /null)
  SolarZenithAngle[nanidx] = !values.f_nan

  SpacecraftAltitude = h5read(file, gg+'/SpacecraftAltitude')
  nanidx = where(SpacecraftAltitude lt -1.e30, /null)
  SpacecraftAltitude[nanidx] = !values.f_nan

  SpacecraftLatitude = h5read(file, gg+'/SpacecraftLatitude')
  nanidx = where(SpacecraftLatitude lt -1.e30, /null)
  SpacecraftLatitude[nanidx] = !values.f_nan

  SpacecraftLongitude = h5read(file, gg+'/SpacecraftLongitude')
  nanidx = where(SpacecraftLongitude lt -1.e30, /null)
  SpacecraftLongitude[nanidx] = !values.f_nan

  ;terrainHeight = h5read(file, gg+'/TerrainHeight')
  time = h5read(file, gg+'/Time')
  nanidx = where(time lt -1.e30, /null)
  time[nanidx] = !values.f_nan

  ;VAA = h5read(file, gg+'/ViewingAzimuthAngle')
  ViewingZenithAngle = h5read(file,gg+'/ViewingZenithAngle')
  nanidx = where(ViewingZenithAngle lt -1.e30, /null)
  ViewingZenithAngle[nanidx] = !values.f_nan


  ; Data Fields
  AerosolIndex = h5read(file,dg+'/AerosolIndex')
  nanidx = where(AerosolIndex lt -1.e30, /null)
  AerosolIndex[nanidx] = !values.f_nan


  AverageResiduals = h5read(file,dg+'/AverageResiduals')
  nanidx = where(AverageResiduals lt -1.e30, /null)
  AverageResiduals[nanidx] = !values.f_nan

  CloudFlag = h5read(file,dg+'/CloudFlag')
  nanidx = where(CloudFlag lt -1.e30, /null)
  CloudFlag[nanidx] = !values.f_nan

  EffectiveCloudFraction = h5read(file, dg+'/EffectiveCloudFraction')
  nanidx = where(EffectiveCloudFraction lt -1.e30, /null)
  EffectiveCloudFraction[nanidx] = !values.f_nan

  EffectiveCloudPressure = h5read(file, dg+'/EffectiveCloudPressure')
  nanidx = where(EffectiveCloudPressure lt -1.e30, /null)
  EffectiveCloudFraction[nanidx] = !values.f_nan

  ExitStatus = h5read(file, dg+'/ExitStatus')
  nanidx = where(ExitStatus lt -1.e30, /null)
  ExitStatus[nanidx] = !values.f_nan

  GlintProbability = h5read(file, dg+'/GlintProbability')
  nanidx = where(GlintProbability lt -1.e30, /null)
  GlintProbability[nanidx] = !values.f_nan

  MeasurementQualityFlags = h5read(file, dg+'/MeasurementQualityFlags')
  nanidx = where(MeasurementQualityFlags lt -1.e30, /null)
  MeasurementQualityFlags[nanidx] = !values.f_nan

  nChannelWavel = h5read(file, dg+'/nChannelWavel')
  nanidx = where(nChannelWavel lt -1.e30, /null)
  nChannelWavel[nanidx] = !values.f_nan

  nFittingWavel = h5read(file, dg+'/nFittingWavel')
  nanidx = where(nFittingWavel lt -1.e30, /null)
  nFittingWavel[nanidx] = !values.f_nan

  nIteration = h5read(file, dg+'/nIteration')
  nanidx = where(nIteration lt -1.e30, /null)
  nIteration[nanidx] = !values.f_nan

  NonGasParameterAPriori = h5read(file, dg+'/NonGasParameterAPriori')
  nanidx = where(NonGasParameterAPriori lt -1.e30, /null)
  NonGasParameterAPriori[nanidx] = !values.f_nan

  NonGasParameterAPrioriError = h5read(file, dg+'/NonGasParameterAPrioriError')
  nanidx = where(NonGasParameterAPrioriError lt -1.e30, /null)
  NonGasParameterAPrioriError[nanidx] = !values.f_nan

  NonGasParameterRetrieved = h5read(file, dg+'/NonGasParameterRetrieved')
  nanidx = where(NonGasParameterRetrieved lt -1.e30, /null)
  NonGasParameterRetrieved[nanidx] = !values.f_nan

  NonGasParameterRetrievedPrecision = h5read(file, dg+'/NonGasParameterRetrievedPrecision')
  nanidx = where(NonGasParameterRetrievedPrecision lt -1.e30, /null)
  NonGasParameterRetrievedPrecision[nanidx] = !values.f_nan

  NonGasParameterRetrievedSolutionError = h5read(file, dg+'/NonGasParameterRetrievedSolutionError')
  nanidx = where(NonGasParameterRetrievedSolutionError lt -1.e30, /null)
  NonGasParameterRetrievedSolutionError[nanidx] = !values.f_nan

  nSmallPixelColumns = h5read(file, dg+'/nSmallPixelColumns')
  nanidx = where(nSmallPixelColumns lt -1.e30, /null)
  nSmallPixelColumns[nanidx] = !values.f_nan

  O3APrioriProfile = h5read(file, dg+'/O3APrioriProfile')
  nanidx = where(O3APrioriProfile lt -1.e30, /null)
  O3APrioriProfile[nanidx] = !values.f_nan

  O3APrioriProfileError = h5read(file, dg+'/O3APrioriProfileError')
  nanidx = where(O3APrioriProfileError lt -1.e30, /null)
  O3APrioriProfileError[nanidx] = !values.f_nan

  O3AveragingKernel = h5read(file, dg+'/O3AveragingKernel')
  nanidx = where(O3AveragingKernel lt -1.e30, /null)
  O3AveragingKernel[nanidx] = !values.f_nan

  O3InformationContent = h5read(file, dg+'/O3InformationContent')
  nanidx = where(O3InformationContent lt -1.e30, /null)
  O3InformationContent[nanidx] = !values.f_nan

  O3NoiseCorrelationMatrix = h5read(file, dg+'/O3NoiseCorrelationMatrix')
  nanidx = where(O3NoiseCorrelationMatrix lt -1.e30, /null)
  O3NoiseCorrelationMatrix[nanidx] = !values.f_nan

  O3RetrievedProfile = h5read(file, dg+'/O3RetrievedProfile')
  nanidx = where(O3RetrievedProfile lt -1.e30, /null)
  O3RetrievedProfile[nanidx] = !values.f_nan

  O3RetrievedProfilePrecision = h5read(file, dg+'/O3RetrievedProfilePrecision')
  nanidx = where(O3RetrievedProfilePrecision lt -1.e30, /null)
  O3RetrievedProfilePrecision[nanidx] = !values.f_nan

  O3RetrievedProfileSolutionError = h5read(file, dg+'/O3RetrievedProfileSolutionError')
  nanidx = where(O3RetrievedProfileSolutionError lt -1.e30, /null)
  O3RetrievedProfileSolutionError[nanidx] = !values.f_nan

  O3StratosphericColumn = h5read(file, dg+'/O3StratosphericColumn')
  nanidx = where(O3StratosphericColumn lt -1.e30, /null)
  O3StratosphericColumn[nanidx] = !values.f_nan

  O3StratosphericColumnPrecision = h5read(file, dg+'/O3StratosphericColumnPrecision')
  nanidx = where(O3StratosphericColumnPrecision lt -1.e30, /null)
  O3StratosphericColumnPrecision[nanidx] = !values.f_nan

  O3StratosphericColumnSolutionError = h5read(file, dg+'/O3StratosphericColumnSolutionError')
  nanidx = where(O3StratosphericColumnSolutionError lt -1.e30, /null)
  O3StratosphericColumnSolutionError[nanidx] = !values.f_nan

  O3TotalColumn = h5read(file, dg+'/O3TotalColumn')
  nanidx = where(O3TotalColumn lt -1.e30, /null)
  O3TotalColumn[nanidx] = !values.f_nan

  O3TotalColumnPrecision = h5read(file, dg+'/O3TotalColumnPrecision')
  nanidx = where(O3TotalColumnPrecision lt -1.e30, /null)
  O3TotalColumnPrecision[nanidx] = !values.f_nan

  O3TotalColumnSolutionError = h5read(file, dg+'/O3TotalColumnSolutionError')
  nanidx = where(O3TotalColumnSolutionError lt -1.e30, /null)
  O3TotalColumnSolutionError[nanidx] = !values.f_nan

  O3TroposphericColumn = h5read(file, dg+'/O3TroposphericColumn')
  nanidx = where(O3TroposphericColumn lt -1.e30, /null)
  O3TroposphericColumn[nanidx] = !values.f_nan

  O3TroposphericColumnPrecision = h5read(file, dg+'/O3TroposphericColumnPrecision')
  nanidx = where(O3TroposphericColumnPrecision lt -1.e30, /null)
  O3TroposphericColumnPrecision[nanidx] = !values.f_nan

  O3TroposphericColumnSolutionError = h5read(file, dg+'/O3TroposphericColumnSolutionError')
  nanidx = where(O3TroposphericColumn lt -1.e30, /null)
  O3TroposphericColumn[nanidx] = !values.f_nan

  OtherGasAPrioriColumnDensity = h5read(file, dg+'/OtherGasAPrioriColumnDensity')
  nanidx = where(OtherGasAPrioriColumnDensity lt -1.e30, /null)
  OtherGasAPrioriColumnDensity[nanidx] = !values.f_nan

  OtherGasAPrioriColumnDensityError = h5read(file, dg+'/OtherGasAPrioriColumnDensityError')
  nanidx = where(OtherGasAPrioriColumnDensityError lt -1.e30, /null)
  OtherGasAPrioriColumnDensityError[nanidx] = !values.f_nan

  OtherGasRetrievedVerticalColumnDensity = h5read(file, dg+'/OtherGasRetrievedVerticalColumnDensity')
  nanidx = where(OtherGasRetrievedVerticalColumnDensity lt -1.e30, /null)
  OtherGasRetrievedVerticalColumnDensity[nanidx] = !values.f_nan

  OtherGasRetrievedVerticalColumnDensityPrecision = h5read(file, dg+'/OtherGasRetrievedVerticalColumnDensityPrecision')
  nanidx = where(OtherGasRetrievedVerticalColumnDensityPrecision lt -1.e30, /null)
  OtherGasRetrievedVerticalColumnDensityPrecision[nanidx] = !values.f_nan

  OtherGasRetrievedVerticalColumnDensitySolutionError = h5read(file, dg+'/OtherGasRetrievedVerticalColumnDensitySolutionError')
  nanidx = where(OtherGasRetrievedVerticalColumnDensitySolutionError lt -1.e30, /null)
  OtherGasRetrievedVerticalColumnDensitySolutionError[nanidx] = !values.f_nan

  ProfileLevelAltitude = h5read(file, dg+'/ProfileLevelAltitude')
  nanidx = where(ProfileLevelAltitude lt -1.e30, /null)
  ProfileLevelAltitude[nanidx] = !values.f_nan

  ProfileLevelPressure = h5read(file, dg+'/ProfileLevelPressure')
  nanidx = where(ProfileLevelPressure lt -1.e30, /null)
  ProfileLevelPressure[nanidx] = !values.f_nan

  ProfileLevelTemperature = h5read(file, dg+'/ProfileLevelTemperature')
  nanidx = where(ProfileLevelTemperature lt -1.e30, /null)
  ProfileLevelTemperature[nanidx] = !values.f_nan

  RMS = h5read(file, dg+'/RMS')
  nanidx = where(RMS lt -1.e30, /null)
  RMS[nanidx] = !values.f_nan

  SurfaceAlbedo = h5read(file, dg+'/SurfaceAlbedo')
  nanidx = where(SurfaceAlbedo lt -1.e30, /null)
  SurfaceAlbedo[nanidx] = !values.f_nan

  TropopauseIndex = h5read(file, dg+'/TropopauseIndex')
  nanidx = where(TropopauseIndex lt -1.e30, /null)
  TropopauseIndex[nanidx] = !values.f_nan


  sz = SIZE(SolarZenithAngle,/dim)

  tmp = h5read(file,gg+'/Time')
  ;timearr = dblarr(60,n_elements(tmp))
  timearr = dblarr(30,n_elements(tmp))
  ;for i=0,59 do timearr[i,*] = tmp
  for i=0,29 do timearr[i,*] = tmp
  timearr[where(timearr lt 0 , /null)] = !values.f_nan

  caldat,(julday(1,1,1993,0,0,0)*86400.0d0+timearr)/86400.0d0,mon,day,year,hour,min,sec

  date = strarr(sz[0],sz[1])
  FOR iy=0, sz[1]-1 DO BEGIN
  FOR ix=0, sz[0]-1 DO BEGIN
    date[ix,iy] =  string(year[ix,iy],mon[ix,iy],day[ix,iy],format='(I4,I02,I02)')+'T'+$
                   string(hour[ix,iy],min[ix,iy],sec[ix,iy],format='(I02,I02,I02)')+'Z'
  ENDFOR 
  ENDFOR
   
  ;@TODO
  ;jul = julday(mon,day,year,hour,min,sec)
 
  xtracks = fltarr(sz[0],sz[1]) 
  lines = fltarr(sz[0],sz[1])
  FOR I=0, sz[0]-1 DO BEGIN
  FOR J=0, sz[1]-1 DO BEGIN
     Xtracks[i,j] = i
     Lines[i,j]   = j
  ENDFOR
  ENDFOR 

;Quality Flags
;Bits 0 to 3 together contain several output error flags:
;  0 - good sample
;  1 - glint contamination (corrected)
;  2 - sza > 84 (degree)
;  3 - 360 residual > threshold
;  4 - residual at unused ozone wavelength > 4 sigma
;  5 - SOI > 4 sigma (SO2 present)
;  6 - non-convergence
;  7 - abs(residual) > 16.0 (fatal)
;  8 - row anomaly error (same as bit 6 in this field)
 ;gflags = convert_bit2flag(GroundPixelQualityFlags,16)
 ;gqflag= gFLAGS(*,*,0)*1 + gFLAGS(*,*,1)*2 + gFLAGS(*,*,2)*4 + gFLAGS(*,*,3)*8

;row anomaly affect = 1 
 ;rowAnomaly = gFLAGS(*,*,6)

; Use climatological cloud pressure is used
  ;useClimPcld = gFlags(*,*,7)

;Measurment Quality Flags
; The measurement quality flag associated with each "scan" line (Bit value
;is 0 for not set and 1 for set):
 ;mflags = convert_bit2flag(MeasurementQualityFlags,8)
 ;mflag= mFLAGS(*,0)*1 ;+ mFLAGS(*,1)*2  + mFLAGS(*,2)*4

dim = size(EffectiveCloudFraction, /dim)

omixtrackidx = rebin(indgen(dim[0]), dim[0], dim[1])
omialongtrackidx = rebin(transpose(indgen(dim[1])),  dim[0], dim[1])
omilatpixcor = fltarr(dim[0], dim[1], 4)
omilonpixcor = fltarr(dim[0], dim[1], 4)

for ix=0, dim[0]-1 do begin
  for iy=0, dim[1]-1 do begin
    omilatpixcor[ix, iy, *]  = [ $
      LatitudePixelCorner[ix, iy], $ 
      LatitudePixelCorner[ix+1, iy], $ 
      LatitudePixelCorner[ix+1, iy+1], $ 
      LatitudePixelCorner[ix, iy+1]]
    omilonpixcor[ix, iy, *]  = [ $
      LongitudePixelCorner[ix, iy], $ 
      LongitudePixelCorner[ix+1, iy], $ 
      LongitudePixelCorner[ix+1, iy+1], $ 
      LongitudePixelCorner[ix, iy+1]]
  ENDFOR
ENDFOR

 result ={AerosolIndex:AerosolIndex, $
   AverageResiduals:AverageResiduals, $
   CloudFlag:CloudFlag, $
   EffectiveCloudFraction:EffectiveCloudFraction, $
   EffectiveCloudPressure:EffectiveCloudPressure, $
   ExitStatus:ExitStatus, $
   GlintProbability:GlintProbability, $
   MeasurementQualityFlags:MeasurementQualityFlags, $
   nChannelWavel:nChannelWavel, $
   nFittingWavel:nFittingWavel, $
   nIteration:nIteration, $
   NonGasParameterAPriori:NonGasParameterAPriori, $
   NonGasParameterAPrioriError:NonGasParameterAPrioriError, $
   NonGasParameterRetrieved:NonGasParameterRetrieved, $
   NonGasParameterRetrievedPrecision:NonGasParameterRetrievedPrecision, $
   NonGasParameterRetrievedSolutionError:NonGasParameterRetrievedSolutionError,$
   nSmallPixelColumns:nSmallPixelColumns, $
   O3APrioriProfile:O3APrioriProfile, $
   O3APrioriProfileError:O3APrioriProfileError, $
   O3AveragingKernel:O3AveragingKernel, $
   O3InformationContent:O3InformationContent, $
   O3NoiseCorrelationMatrix:O3NoiseCorrelationMatrix, $
   O3RetrievedProfile:O3RetrievedProfile, $
   O3RetrievedProfilePrecision:O3RetrievedProfilePrecision, $
   O3RetrievedProfileSolutionError:O3RetrievedProfileSolutionError, $
   O3StratosphericColumn:O3StratosphericColumn, $
   O3StratosphericColumnPrecision:O3StratosphericColumnPrecision, $
   O3StratosphericColumnSolutionError:O3StratosphericColumnSolutionError, $
   O3TotalColumn:O3TotalColumn, $
   O3TotalColumnPrecision:O3TotalColumnPrecision, $
   O3TotalColumnSolutionError:O3TotalColumnSolutionError, $
   O3TroposphericColumn:O3TroposphericColumn, $
   O3TroposphericColumnPrecision:O3TroposphericColumnPrecision, $
   O3TroposphericColumnSolutionError:O3TroposphericColumnSolutionError, $
   OtherGasAPrioriColumnDensity:OtherGasAPrioriColumnDensity, $
   OtherGasAPrioriColumnDensityError:OtherGasAPrioriColumnDensityError, $
   OtherGasRetrievedVerticalColumnDensity:OtherGasRetrievedVerticalColumnDensity, $
   OtherGasRetrievedVerticalColumnDensityPrecision:OtherGasRetrievedVerticalColumnDensityPrecision, $
   OtherGasRetrievedVerticalColumnDensitySolutionError:OtherGasRetrievedVerticalColumnDensitySolutionError, $
   ProfileLevelAltitude:ProfileLevelAltitude, $
   ProfileLevelPressure:ProfileLevelPressure, $
   ProfileLevelTemperature:ProfileLevelTemperature, $
   RMS:RMS, $
   SurfaceAlbedo:SurfaceAlbedo, $
   TropopauseIndex:TropopauseIndex,$
   GroundPixelQualityFlags:GroundPixelQualityFlags, $
   Latitude:Latitude, $
   LatitudePixelCorner:LatitudePixelCorner, $
   LatitudePixCor:omilatpixcor, $
   Longitude:Longitude, $
   LongitudePixelCorner:LongitudePixelCorner, $
   LongitudePixCor:omilonpixcor, $
   XtrackIndex:omixtrackidx, $
   AlongtrackIndex:omialongtrackidx, $
   RelativeAzimuthAngle:RelativeAzimuthAngle, $
   SecondsInDay:SecondsInDay, $
   SolarZenithAngle:SolarZenithAngle, $
   SpacecraftAltitude:SpacecraftAltitude, $
   SpacecraftLatitude:SpacecraftLatitude, $
   SpacecraftLongitude:SpacecraftLongitude, $
   ;TerrainHeight:TerrainHeight, $
   Time:Time, $
   ViewingAzimuthAngle:ViewingZenithAngle};, $
   ;rowAnomaly:rowAnomaly, $
   ;mflag:mflag,$
   ;useClimpcld:useClimpcld}
 return
end
