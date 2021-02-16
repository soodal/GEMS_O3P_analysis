pro read_l2_v8, file,result
               
 dg = '/HDFEOS/SWATHS/OMI Column Amount O3/Data Fields'
 gg = '/HDFEOS/SWATHS/OMI Column Amount O3/Geolocation Fields'

  lat = h5read(file,gg+'/Latitude')
  lon = h5read(file,gg+'/Longitude')
  SZA = h5read(file,gg+'/SolarZenithAngle')
  VZA = h5read(file,gg+'/ViewingZenithAngle')
  RAA = h5read(file,gg+'/RelativeAzimuthAngle')

  toz = h5read(file,dg+'/ColumnAmountO3')  
  blo3 = h5read(file,dg+'/O3BelowCloud')  
  ozn1 = h5read(file,dg+'/StepOneO3')
  ozn2 = h5read(file,dg+'/StepTwoO3')
  qf =  h5read(file,dg+'/QualityFlags') 
  EL  = h5read(file,dg+'/LayerEfficiency')
  dndt = h5read(file,dg+'/dN_dT')
  dndr = h5read(file,dg+'/dN_dR')
  
  ref331 = h5read(file,dg+'/Reflectivity331')
  ref360 = h5read(file,dg+'/Reflectivity360')
  pter = h5read(file,dg+'/TerrainPressure')
  pc = h5read(file,dg+'/CloudPressure')
  ai = h5read(file,dg+'/UVAerosolIndex' )
  cf = h5read(file,dg+'/fc')

 xqf =  h5read(file,gg+'/XTrackQualityFlags')
 gqf =  h5read(file,gg+'/GroundPixelQualityFlags')
 mqf =  h5read(file,dg+'/MeasurementQualityFlags')

  sz = SIZE(sza,/dim)

  tmp = h5read(file,gg+'/Time')
  time = dblarr(60,n_elements(tmp))
  for i=0,59 do time[i,*] = tmp
  caldat,(julday(1,1,1993,0,0,0)*86400.0d0+time)/86400.0d0,mon,day,year,hour,min,sec

  date = strarr(sz[0],sz[1])
  FOR iy=0, sz[1]-1 DO BEGIN
  FOR ix=0, sz[0]-1 DO BEGIN
    date[ix,iy] =  string(year[ix,iy],mon[ix,iy],day[ix,iy],format='(I4,I02,I02)')+'T'+$
                   string(hour[ix,iy],min[ix,iy],sec[ix,iy],format='(I02,I02,I02)')+'Z'
  ENDFOR & ENDFOR
   
  jul = julday(mon,day,year,hour,min,sec)
 
  xtracks = fltarr(sz[0],sz[1]) & lines = fltarr(sz[0],sz[1])
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
 flags = convert_bit2flag(qf,16)
 qflag= FLAGS(*,*,0)*1 + FLAGS(*,*,1)*2 + FLAGS(*,*,2)*4 + FLAGS(*,*,3)*8

;row anomaly affect = 1 
 rowAnomaly = FLAGS(*,*,6)

; Use climatological cloud pressure is used
  useClimPcld = Flags(*,*,7)


;Xtrack Flags
;The cross track quality flags assigned to each pixel in
;OMI L1B data. Flags indicate detection of the OMI row
;anomaly by analysis of the L1B data.
  ;0 - Not affected
  ;1 - Affected, Not corrected, do not use
  ;2 - Slightly affected, not corrected, use with caution
  ;3 - Affected, corrected, use with caution
  ;4 - Affected, corrected, use pixel
 flags = convert_bit2flag(xqf,8)
 xflag= FLAGS(*,*,0)*1 + FLAGS(*,*,1)*2  + FLAGS(*,*,2)*4


;Measurment Quality Flags
; The measurement quality flag associated with each "scan" line (Bit value
;is 0 for not set and 1 for set):
 flags = convert_bit2flag(mqf,8)
 mflag= FLAGS(*,0)*1 ;+ FLAGS(*,1)*2  + FLAGS(*,2)*4


 result ={lat:lat, lon:lon,sza:sza,vza:vza,$
          qf:qflag,xflag:xflag,mflag:mflag,$
          rowAnomaly:rowAnomaly, uClimpcld:useClimpcld,$
          toz:toz,$
          dndt:dndt, dndr:dndr, EL:EL, raa:raa, pc:pc, pter:pter,$
          ref331:ref331, ref360:ref360,$
          xtracks:xtracks, lines:lines,date:date,jul:jul,ai:ai,blo3:blo3,cf:cf,ozn1:ozn1,ozn2:ozn2}
 
  end
