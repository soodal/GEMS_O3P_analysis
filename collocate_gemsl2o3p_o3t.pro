pro collocate_gemsl2o3p_omo3pr, year, month, day, hour, minute, gemsvals, $
  omivals, on_under_300=on_under_300

if not keyword_Set(on_under_300) then begin
  print, 'Total Ozone'
  on_under_300 = 0
endif

limit=[-10, 80, 60, 160]

yyyy = string(year, format='(i04)')
mm = string(month, format='(i02)')
dd = string(day, format='(i02)')
hh = string(hour, format='(i02)')
mi = string(minute, format='(i02)')

savpath = './collocate_gemsl2o3p_omo3pr/'
savfile = savpath + 'col_gems_omi_'+yyyy+mm+dd+'_'+hh+mi+'.sav'

savefile =  file_test(savfile) 
;savefile = 0
if not savefile then BEGIN

  ; initialize

  omialts = []  
  omigroundpqfs = []
  omilons = []  
  omilats = []  
  omipress = []
  omisaas = []
  omiszas = []
  omitimes = []  
  omivaas = []
  omivzas = []
  omiaods = []
  omiaprioricovs = []
  omiaks = []
  omipclds =  []
  omitozs  =  []
  omidfss = []
  omiecf1s = []
  omiecf2s = []
  omiicis = []
  omimqfs = []
  omino2s = []
  omio3s = fltarr(18, 1)
  ominumiters = []

  omixtracks = []  
  omilines = []
  ;pixlons = []  
  ;pixlats = []  
  ;omir331s = [] 
  omirowanomalys = []
  omirefcfs = []

  ; path for GEMS
  gemsfn = '/data2/L2_GEMS/val_1008/GK2_GEMS_O3P_' $
    + yyyy + mm + dd + '_' + hh + mi + '.nc4'

  ; read GEMS O3P data
  gemsvars = ds_read_gems_l2_o3p(gemsfn)

;;----------------------------
;; GEMS vertical column layer for under 300 hPa
;;----------------------------
  o3 = gemsvars.o3
  o3size = size(gemsvars.o3, /dim)
  ;altitude = gemsvars.Altitude
  gemspres = gemsvars.Pressure
  gemslayero3 = fltarr(o3size[1:2])
  dim1 = o3size[1]
  dim2 = o3size[2]
  ;gemslayero3[*,*] = !values.f_nan
  for iy=0, dim2-1 do begin
  for ix=0, dim1-1 do begin
    for ilevel=24, 0, -1 do begin
      ilayer = ilevel-1
      if gemspres[ilevel, ix, iy] gt 300. $
          and gemspres[ilevel-1, ix, iy] gt 300. then begin
        gemslayero3[ix, iy] = gemslayero3[ix, iy] + o3[ilayer, ix, iy]

      endif else if gemspres[ilevel, ix, iy] gt 300. $
          and gemspres[ilevel-1, ix, iy] le 300. then begin

        gemslayero3[ix, iy] = gemslayero3[ix, iy] + o3[ilayer, ix, iy]*$
          (gemspres[ilevel, ix, iy]-300.)/$
          (gemspres[ilevel, ix, iy]-gemspres[ilevel-1, ix, iy])
        break
      endif
    endfor
  endfor
  endfor
  gemslayero3[where(gemslayero3 lt 0)] = !values.f_nan

  ; GEMS vertical column layer for under 10 km
  ;for iy=0, dim2-1 do begin
  ;for ix=0, dim1-1 do begin
    ;for ilevel=24, 0, -1 do begin
      ;if altitude[ilevel, ix, iy] lt 10. $
          ;and altitude[ilevel-1, ix, iy] lt 10. then begin
        ;gemslayero3[ix, iy] = gemslayero3[ix, iy] + o3[ilevel-1, ix, iy]

      ;endif else if altitude[ilevel, ix, iy] lt 10. $
          ;and altitude[ilevel-1, ix, iy] ge 10. then begin

        ;gemslayero3[ix, iy] = gemslayero3[ix, iy] + o3[ilevel-1, ix, iy]*$
          ;(10.-altitude[ilevel, ix, iy])/$
          ;(altitude[ilevel-1, ix, iy]-altitude[ilevel, ix, iy])
        ;break
      ;endif
    ;endfor
  ;endfor
  ;endfor
  ;gemslayero3[where(gemslayero3 le 0)] = !values.f_nan

  ; path for OMI

  omipath = '/data2/OMI/gdata/' + yyyy + '/L2-OMO3PR/' 
  omipixpath = '/data2/OMI/gdata/' + yyyy + '/L2-OMPIXCOR/' 

  omifiles=FILE_SEARCH(omipath+'OMI-Aura_L2-OMO3PR_'+yyyy+'m'+mm+dd+'*.he5', count=nfiles)

  ; read OMI data
  for ifile=0, nfiles-1 do begin
    fi = file_info(omifiles[ifile])
    if fi.size gt 1500 then begin
      ds_READ_omi_L2_omo3pr, omifiles[ifile], omiresult

      tmp=WHERE(omiresult.longitude GE limit[1] AND omiresult.longitude LE limit[3] AND $
                omiresult.latitude GE limit[0] AND omiresult.latitude LE limit[2], ntmp)

      print, '  READ omifile : ', omifiles, ntmp
      IF ntmp NE 0 THEN BEGIN 
        header='OMI-Aura_L2-OMPIXCOR_'
        time=STRMID(omifiles[ifile],46+4,14,/rev)
        print,time

        ;pixfile=FILE_SEARCH(omipixpath + header+time+'*.he5', count=npixfile)
        ;IF npixfile NE 1 THEN BEGIN
          ;PRINT, '  No Pixfile'
          ;STOP
        ;ENDIF ; npixfile

        ;pixfi = file_info(omifiles[ifile])
        ;if pixfi.size gt 1500 then begin

          ;ds_READ_omi_L2_PIXCOR, pixfile, pixcor
          ;pixdim=SIZE(pixcor.corlon1,/dim)
          ;nx=pixdim[0]
          ;ny=pixdim[1]
          ;pixlon=REFORM(pixcor.corlon1,[nx*ny,4])
          ;pixlat=REFORM(pixcor.corlat1,[nx*ny,4])     
          ;pixlats  =  [pixlats, pixlat]
          ;pixlons  =  [pixlons, pixlon]

          ;pixlon=pixlon[tmp,*]
          ;pixlat=pixlat[tmp,*]
               ;omitime=pixcor.jtime1[tmp]

          omialts = [omialts, omiresult.altitude[tmp]]
          omigroundpqfs   =  [omigroundpqfs,   omiresult.GroundPixelQualityFlags[tmp]     ]
          omilats  =  [omilats, omiresult.Latitude[tmp] ]
          omilons  =  [omilons, omiresult.Longitude[tmp] ]
          omipres = omiresult.Pressure
          omipressize = size(omipres, /dim)
          omipres = reform(omipres, [omipressize[0], omipressize[1]*omipressize[2]])
          omipress = [[omipress], [omipres[*, tmp]]]
          omisaas = [omisaas, omiresult.SolarAzimuthAngle[tmp]]
          omiszas = [omiszas, omiresult.SolarZenithAngle[tmp]]
          omitimes = [omitimes, omiresult.time[tmp]]
          omivaas = [omivaas, omiresult.ViewingAzimuthAngle[tmp]]
          omivzas = [omivzas, omiresult.ViewingZenithAngle[tmp]]
          omiaods = [omiaods, omiresult.AerosolOpticalThickness[tmp]]
          omiaprioricovs = [omiaprioricovs, omiresult.APrioricov[tmp]]
          omiaks = [omiaks, omiresult.AveragingKernel[tmp]]
          omipclds =  [omipclds, omiresult.CloudPressure[tmp]     ]
          omitozs  =  [omitozs, omiresult.ColumnAmountO3[tmp] ]
          omidfss = [omidfss, omiresult.DegreesOfFreedomForSignal[tmp]]
          omiecf1s = [omiecf1s, omiresult.EffectiveCloudFractionUV1[tmp]]
          omiecf2s = [omiecf2s, omiresult.EffectiveCloudFractionUV2[tmp]]
          omiicis = [omiicis, omiresult.InstrumentConfigurationId[tmp]]
          omimqfs = [omimqfs, omiresult.MeasurementQualityFlags[tmp]]
          omino2s = [omino2s, omiresult.NO2ColumnAmount[tmp]]
          ominumiters = [ominumiters, omiresult.NumberOfIterations[tmp]]
          omio3 = omiresult.o3
          omio3size = size(omio3, /dim)
          omio3 = reform(omio3, [omio3size[0], omio3size[1]*omio3size[2]])
          omio3[where(omio3 lt 0, /null)]=!values.f_nan
          omio3s = [[omio3s], [omio3[*, tmp]]]
          
          ;omisize = size(omilon, /dim)
          ;omio3 = reform(o3result.o3, [18, omisize[0]*omisize[1]]) 
          ;omio3s = [omio3s, omiresult.O3[tmp]]

  ;o3ap = h5read(file, dg+'/O3Apriori')
  ;o3aperr = h5read(file, dg+'/O3AprioriError')
  ;o3prec = h5read(file, dg+'/O3Precision')
  ;pqf = h5read(file, dg+'/ProcessingQualityFlags')
  ;refcf = h5read(file, dg+'/ReflectanceCostFunction')
  ;res = h5read(file, dg+'/ResidualsOfFit')
  ;rmserroffit = h5read(file, dg+'/RootMeanSquareErrorOfFit')
  ;sprad = h5read(file, dg+'/SmallPixelRadiance')
  ;spradpnt = h5read(file, dg+'/SmallPixelRadiancePointer')
  ;spradvar = h5read(file, dg+'/SmallPixelRadianceVariance')
  ;so2 = h5read(file, dg+'/SO2ColumnAmount')
  ;statecost = h5read(file, dg+'/StateCostFunction')
  ;statevectorspecies = h5read(file, dg+'/StateVectorSpecies')
  ;temp = h5read(file, dg+'/Temperature')
  ;terrainref1 = h5read(file, dg+'/TerrainReflectivityUV1')
  ;terrainref2 = h5read(file, dg+'/TerrainReflectivityUV2')

          omirefcfs = [omirefcfs, omiresult.ReflectanceCostFunction[tmp]]

          omirowanomalys = [omirowanomalys, omiresult.rowanomaly[tmp]]
          
          omixtracks  =  [omixtracks,omiresult.xtracks[tmp]]
          omilines =  [omilines,omiresult.lines[tmp]]

        ;endif
      ENDIF
    ENDIF  ; ntmp
  ENDFOR  ; ifile


;;----------------------------
;; OMI vertical column Layer for under 300 hPa 
;;----------------------------

  omitozsize = size(omitozs, /dim)
  omilayero3 = fltarr(omitozsize[0])

  for ip=0,omitozsize[0]-1 do BEGIN
    for ilevel=18, 0, -1 do BEGIN
      ilayer = ilevel - 1
      if omipress[ilevel, ip] gt 300. $
          and omipress[ilevel-1, ip] gt 300. then BEGIN
        omilayero3[ip] = omilayero3[ip] + omio3s[ilayer, ip]
        
      ENDIF else if omipress[ilevel, ip] gt 300. $
        and omipress[ilevel-1, ip] le 300. then begin 
        ; upper bound is higher than 300hPa 

        omilayero3[ip] = omilayero3[ip] + omio3s[ilayer, ip]*$
          (omipress[ilevel, ip]-300.)/$
          (omipress[ilevel, ip] - omipress[ilevel-1, ip])
        break
      ENDIF
    endfor
  ENDFOR

  omivars={lon:omilons, lat:omilats, time:omitimes, toz:omitozs,$
    SolarZenithAngle:omiszas, $
    SolarAzimuthAngle:omisaas, $
    ReflectanceCostFunction:omirefcfs, $
    EffectiveCloudFractionUV1:omiecf1s, $
    EffectiveCloudFractionUV2:omiecf2s, $
    omio3Under300hPa:omilayero3, $
    ;pixlat:pixlats, pixlon:pixlons, pixs:omixtracks, lines:omilines,$
    xtracks:omixtracks, $
    lines:omilines, $
    rowanomaly:omirowanomalys, $
    ;ref331:omir331s, $
    ;cf:omicfs, $
    refcf:omirefcfs, $
    pcld:omipclds, qf:omigroundpqfs};,bloz:omiblozs}


  ; collocation

  omilon = omivars.lon
  omilat = omivars.lat
  omitime = omivars.time
  omitoz = omivars.toz
  omi300 = omivars.omio3Under300hPa
  ;omipixlat = omivars.pixlat
  ;omipixlon = omivars.pixlon
  omixtracks = omivars.xtracks
  omilines = omivars.lines
  omipcld = omivars.pcld
  omiqf = omivars.qf
  omirowanomaly = omivars.rowAnomaly


  nanidx = where(omitime le 0, /null)
  omitime[nanidx] = !values.f_nan
  caldat, julday(1,1,1993, 0, 0) + omitime/60./60./24., omimon, omiday, $
    omiyear, omihour, omiminute, omisec

  omimon[nanidx] = 1
  omiday[nanidx] = 1
  omiyear[nanidx] = 1
  omihour[nanidx] = 0
  omiminute[nanidx] = 0
  omisec[nanidx] = 0
  omijulday = julday(omimon, omiday, omiyear, omihour, omiminute, omisec) 
  
  omijulday[nanidx] = !values.f_nan
  
  gemslon = gemsvars.Longitude
  gemslat = gemsvars.Latitude
  ;gemstime = gemsvars.Time
  gemslon[where(gemslon lt -180, /null)] = !values.f_nan
  gemslat[where(gemslat lt -90, /null)] = !values.f_nan
  ;gemstime[where(gemstime lt 0, /null)] = !values.f_nan
  gemstime = julday(month, day, year, hour, minute)
  gemssize = size(gemslon, /dim)
  
  yn_omi_cross_time = intarr(n_elements(omilon))
  yn_omi_cross_time[*] = 0
  yn_omi_cross_time[where(abs(omijulday-gemstime) lt 1./24./2., npix, /null)] = 1

  ;closest_idx = lonarr(gemssize[0], gemssize[1])
  ;closest_idx[*] = -999
  omisize = size(omilon, /dim)
  closest_idx = lonarr(omisize[0])
  closest_idx[*] = -999

  ;for iy = 0, gemssize[1]-1 do BEGIN
    ;for ix = 0, gemssize[0]-1 do BEGIN
      ;x = gemslon[ix, iy]
      ;y = gemslat[ix, iy]


      ;if finite(x) and finite(y) then begin
        ;result = search_closest_pixel(omilon, omilat, x, y)
        ;;result = array_indices(gemslon, result)
        ;closest_idx[ix, iy] = result
      ;ENDIF
    ;ENDFOR
  ;ENDFOR

  for ix = 0, omisize[0]-1 do BEGIN
    x = omilon[ix]
    y = omilat[ix]

    if finite(x) and finite(y) then begin
      result = search_closest_pixel(gemslon, gemslat, x, y)
      ;result = array_indices(gemslon, result)
      closest_idx[ix] = result
    ENDIF
  ENDFOR
  
  save, filename=savfile, gemsvars, gemslayero3, omivars, closest_idx, yn_omi_cross_time

  ;gemstoz = reform(gemsvars.ColumnAmountO3[0, *, *])
  ;gemsvals = []
  ;omivals = []
  ;for iy=0, gemssize[1]-1 do BEGIN
    ;for ix=0, gemssize[0]-1 do BEGIN
      ;x = gemslon[ix, iy]
      ;y = gemslat[ix, iy]

      ;if closest_idx[ix, iy] GE 0 then begin
        ;omival = omitoz[closest_idx[ix, iy]]
        ;gemsval = gemstoz[ix, iy]

        ;if omival ge 0 and omival le 1000 AND $
            ;gemsval ge 0 and gemsval le 1000 then begin
          ;omivals = [omivals, omival]
          ;gemsvals = [gemsvals, gemsval]
        ;endif

      ;endif
    ;ENDFOR
  ;ENDFOR

endif else BEGIN
  RESTORE, savfile
ENDELSE
  gemstoz = reform(gemsvars.ColumnAmountO3[0, *, *])
  gemsvals = []
  omivals = []
  omisize = size(omivars.lon, /dim)
  gemsecf = gemsvars.effectivecloudfractionuv
  for ix=0, omisize[0]-1 do BEGIN

    if closest_idx[ix] GE 0 and yn_omi_cross_time[ix] gt 0 then begin
      if omivars.ReflectanceCostFunction[ix] lt 30 $
          and omivars.SolarZenithAngle[ix] lt 5000 $
          ;and omivars.xtracks[ix] ge 3 $
          ;and omivars.xtracks[ix] le 25 $
          ;and omivars.EffectiveCloudFractionUV1[ix] le 0.2 $
          and gemsecf[closest_idx[ix]] le 0.2 $
        then begin

        if on_under_300 then begin
          omival = omivars.omio3under300hpa[ix]
          gemsval = gemslayero3[closest_idx[ix]]
        endif else begin
          omival = omivars.toz[ix] 
          gemsval = gemstoz[closest_idx[ix]]
        endelse

        ;print,omival, gemsval

        if omival ge 0 and omival le 1000 AND $
            gemsval ge 0 and gemsval le 1000 then begin
          omivals = [omivals, omival]
          gemsvals = [gemsvals, gemsval]
        endif
      endif

    endif
  ENDFOR
  ;plot_omi_satproj, omivars.lon, omivars.lat, omivars.toz,$
    ;filename='plot/omi_toz_20200806T0345.png', $
    ;title='OMI TOZ', range=[220, 350]
  ;plot_omi_satproj, omivars.lon, omivars.lat, omivars.rowanomaly,$
    ;filename='plot/OMI_rowanomaly_20200806T0345.png', $
    ;title='OMI Row Anomaly', range=[0, 1.]
  return
end
