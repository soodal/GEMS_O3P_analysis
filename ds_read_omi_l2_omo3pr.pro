pro ds_read_omi_l2_omo3pr, file,result
               
 dg = '/HDFEOS/SWATHS/O3Profile/Data Fields'
 gg = '/HDFEOS/SWATHS/O3Profile/Geolocation Fields'

  alt = h5read(file,gg+'/Altitude')
  Groundpqf = h5read(file, gg+'/GroundPixelQualityFlags')
  lat = h5read(file,gg+'/Latitude')
  lon = h5read(file,gg+'/Longitude')
  pres = h5read(file, gg+'/Pressure')
  SAA = h5read(file,gg+'/SolarAzimuthAngle')
  SZA = h5read(file,gg+'/SolarZenithAngle')
  spacecraftalt = h5read(file, gg+'/SpacecraftAltitude')
  spacecraftLat = h5read(file, gg+'/SpacecraftLatitude')
  spacecraftLon = h5read(file, gg+'/SpacecraftLongitude')
  terrainHeight = h5read(file, gg+'/TerrainHeight')
  time = h5read(file, gg+'/Time')
  VAA = h5read(file, gg+'/ViewingAzimuthAngle')
  VZA = h5read(file,gg+'/ViewingZenithAngle')

  aod = h5read(file, dg+'/AerosolOpticalThickness')
  APrioricov = h5read(file, dg+'/APrioriCovarianceMatrix')
  ak = h5read(file, dg+'/AveragingKernel')
  cloudPressure = h5read(file,dg+'/CloudPressure')
  ColumnAmountO3 = h5read(file, dg+'/ColumnAmountO3')
  toz = h5read(file,dg+'/ColumnAmountO3')  
  covmat = h5read(file, dg+'/CovarianceMatrix')
  dfs = h5read(file, dg+'/DegreesOfFreedomForSignal')
  ecf1 = h5read(file, dg+'/EffectiveCloudFractionUV1')
  ecf2 = h5read(file, dg+'/EffectiveCloudFractionUV2')
  ici = h5read(file, dg+'/InstrumentConfigurationId')
  mqf = h5read(file, dg+'/MeasurementQualityFlags')
  no2 = h5read(file, dg+'/NO2ColumnAmount')
  numiter = h5read(file, dg+'/NumberOfIterations')
  o3 = h5read(file, dg+'/O3')
  o3ap = h5read(file, dg+'/O3Apriori')
  o3aperr = h5read(file, dg+'/O3AprioriError')
  o3prec = h5read(file, dg+'/O3Precision')
  pqf = h5read(file, dg+'/ProcessingQualityFlags')
  refcf = h5read(file, dg+'/ReflectanceCostFunction')
  res = h5read(file, dg+'/ResidualsOfFit')
  rmserroffit = h5read(file, dg+'/RootMeanSquareErrorOfFit')
  sprad = h5read(file, dg+'/SmallPixelRadiance')
  spradpnt = h5read(file, dg+'/SmallPixelRadiancePointer')
  spradvar = h5read(file, dg+'/SmallPixelRadianceVariance')
  so2 = h5read(file, dg+'/SO2ColumnAmount')
  statecost = h5read(file, dg+'/StateCostFunction')
  statevectorspecies = h5read(file, dg+'/StateVectorSpecies')
  temp = h5read(file, dg+'/Temperature')
  terrainref1 = h5read(file, dg+'/TerrainReflectivityUV1')
  terrainref2 = h5read(file, dg+'/TerrainReflectivityUV2')

  sz = SIZE(sza,/dim)

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
 gflags = convert_bit2flag(groundpqf,16)
 gqflag= gFLAGS(*,*,0)*1 + gFLAGS(*,*,1)*2 + gFLAGS(*,*,2)*4 + gFLAGS(*,*,3)*8

;row anomaly affect = 1 
 rowAnomaly = gFLAGS(*,*,6)

; Use climatological cloud pressure is used
  useClimPcld = gFlags(*,*,7)

;Measurment Quality Flags
; The measurement quality flag associated with each "scan" line (Bit value
;is 0 for not set and 1 for set):
 mflags = convert_bit2flag(mqf,8)
 mflag= mFLAGS(*,0)*1 ;+ mFLAGS(*,1)*2  + mFLAGS(*,2)*4

 result ={Altitude:alt, $
   GroundPixelQualityFlags:groundpqf, $
   Latitude:lat, Longitude:lon, $
   Pressure:pres, $
   SolarAzimuthAngle:saa, SolarZenithAngle:sza, $
   SpacecraftAlt:spacecraftAlt, $
   SpacecraftLat:spacecraftLat, $
   SpacecraftLon:spacecraftLon, $
   TerrainHeight:TerrainHeight, $
   time:timearr, $
   ViewingAzimuthAngle:vaa, ViewingZenithAngle:vza, $
   AerosolOpticalThickness:aod, $
   APrioricov:APrioricov, $
   AveragingKernel:ak, $
   CloudPressure:cloudPressure, $
   ColumnAmountO3:ColumnAmountO3, $
   CovarianceMatrix:covmat, $
   DegreesOfFreedomForSignal:dfs, $
   EffectiveCloudFractionUV1:ecf1, $
   EffectiveCloudFractionUV2:ecf2, $
   InstrumentConfigurationId:ici, $
   MeasurementQualityFlags:mqf, $
   NO2ColumnAmount:no2, $
   NumberOfIterations:numiter, $
   O3:o3, $
   o3Apriori:o3ap, o3aprioriError:o3aperr, O3Precision:o3prec, $
   ProcessingQualityFlags:pqf, $
   ReflectanceCostFunction:refcf, $
   ResidualsOfFit:res, $
   RootMeanSquareErrorOfFit:rmserroffit, $
   SmallPixelRadiance:sprad, $
   SmallPixelRadiancePointer:spradpnt, $
   SmallPixelRadianceVariance:spradvar, $
   SO2ColumnAmount:so2, $
   StateCostFunction:statecost, $
   StateVectorSpecies:statevectorspecies, $
   Temperature:temp, $
   TerrainReflectivityUV1:terrainref1, $
   TerrainReflectivityUV2:terrainref2, $
   mflag:mflag,$
   rowAnomaly:rowAnomaly, $
   uClimpcld:useClimpcld,$
   toz:toz,$
   xtracks:xtracks, lines:lines, date:date};, $
end
