pro collocate_gemsl2o3p_omil2profoz, year, month, day, hour, minute, $
  gemsvals, $
  omivals, $
  gemsidx=gemsidx, $
  hpa=pressure_limit, height=height_limit

if not keyword_Set(pressure_limit) then begin
  pressure_limit = 0
endif

if not keyword_Set(height_limit) then begin
  height_limit = 0
endif

limit=[-10, $ ;minimum latitude
  60, $ ;minimum longitude
  ;60, $ ;maximum latitude
  60, $ ;maximum latitude
  160] ;maximum longitude

yyyy = string(year, format='(i04)')
mm = string(month, format='(i02)')
dd = string(day, format='(i02)')
hh = string(hour, format='(i02)')
mi = string(minute, format='(i02)')

savpath = './collocate_gemsl2o3p_profoz/'
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
  ;omisaas = []
  omiszas = []
  omitimes = []  
  ;omivaas = []
  ;omivzas = []
  omiaods = []
  omiaprioricovs = []
  omiaks = []
  omipclds =  []
  omitozs  =  []
  omidfss = []
  omiecf = []
  omiecp = []
  omiicis = []
  omimqfs = []
  omino2s = []
  omio3s = fltarr(24, 1)
  xtrackindices = []
  alongtrackindices = []
  ominumiters = []

  omixtracks = []  
  omilines = []
  ;pixlons = []  
  ;pixlats = []  
  ;omir331s = [] 
  omirowanomalys = []
  omirefcfs = []
  omilayero3s = []

  ; path for GEMS
  gemsfn = ds_get_gems_o3p_filename(yyyy, mm, dd, hh, mi)
  ;gemsfn = '/data2/L2_GEMS/val_1008/GK2_GEMS_O3P_' $
    ;+ yyyy + mm + dd + '_' + hh + mi + '.nc4'
  if gemsfn ne !null then begin

    ; read GEMS O3P data
    print, gemsfn

    varlist = ['EffectiveCloudFractionUV', 'ProcessingQualityFlags', $
             'O3' , $
             ;'O3Apriori', 'O3AprioriError', $
             ;'CloudPressure', $
             ;'O3Apriori', 'O3AprioriError',$
             'ColumnAmountO3', $
             'SolarZenithAngle', $
             ;'SimulatedRadiances', $
             'Latitude' ,'Longitude', $
             'Time','Altitude' ,    $
             ;'Pressure', 'TropopausePressure', $
             ;'Wavelengths', $
             'WavelengthsWholeRange']
    gemsvars = ds_read_gems_l2_o3p(gemsfn, varlist=varlist)

  ;;----------------------------
  ;; GEMS vertical column layer for under specific height
  ;;----------------------------

    ds_gems_l2o3p_accum, gemsvars, gemslayero3, height=10.


    ; path for OMI

    omipath = '/data2/OMI/gdata/' + yyyy + '/L2-PROFOZ/' 

    omifiles=FILE_SEARCH(omipath+'OMI-Aura_L2-PROFOZ_'+yyyy+'m'+mm+dd+'*.he5')
    ;omifiles = omifiles[$
      ;where($
        ;strmatch(omifiles, '*OMI-Aura_L2-PROFOZ_'+yyyy+'m'+mm+dd+'t0[1-6]*.he5') eq 1)]

    nfiles = n_elements(omifiles)
    ; read OMI data
    for ifile=0, nfiles-1 do begin
    ;for ifile=0, 5 do begin
      fi = file_info(omifiles[ifile])

      omiobstime = strmid(omifiles[ifile], 46, 14, /reverse)
      omiorbitnumb = strmid(omifiles[ifile], 31, 6, /reverse)
      if fi.size gt 1500 then begin
        print, omifiles[ifile]
        ds_READ_omi_L2_profoz, omifiles[ifile], profozresult

        roiidx=WHERE(profozresult.longitude GE limit[1] AND $
                  profozresult.longitude LE limit[3] AND $
                  profozresult.latitude GE limit[0] AND $
                  profozresult.latitude LE limit[2] and $
                  profozresult.EffectiveCloudFraction LE 0.2, ntmp, /null)
        print, ntmp

        if ntmp ge 1 then begin 
          roiindices = array_indices(profozresult.latitude, roiidx)
          ds_omi_l2_profoz_accum, profozresult, omilayero3, height=10.
          dim = size(profozresult.EffectiveCloudFraction, /dim)

          omixtrackidx = rebin(indgen(dim[0]), dim[0], dim[1])
          omialongtrackidx = rebin(transpose(indgen(dim[1])),  dim[0], dim[1])
          omilatpixcor = fltarr(dim[0], dim[1], 4)
          omilonpixcor = fltarr(dim[0], dim[1], 4)

          for ix=0, dim[0]-1 do begin
            for iy=0, dim[1]-1 do begin
              omilatpixcor[ix, iy, *]  = [ $
                profozresult.LatitudePixelCorner[ix, iy], $ 
                profozresult.LatitudePixelCorner[ix+1, iy], $ 
                profozresult.LatitudePixelCorner[ix+1, iy+1], $ 
                profozresult.LatitudePixelCorner[ix, iy+1]]
              omilonpixcor[ix, iy, *]  = [ $
                profozresult.LongitudePixelCorner[ix, iy], $ 
                profozresult.LongitudePixelCorner[ix+1, iy], $ 
                profozresult.LongitudePixelCorner[ix+1, iy+1], $ 
                profozresult.LongitudePixelCorner[ix, iy+1]]
            ENDFOR
          ENDFOR

          xtrackindex = profozresult.XtrackIndex
          alongtrackindex = profozresult.AlongtrackIndex

          xtrackindices = [xtrackindices, xtrackindex[roiidx]]
          alongtrackindices = [alongtrackindices, alongtrackindex[roiidx]]

          time=STRMID(omifiles[ifile],46,14,/rev)

          sz = size(profozresult.O3RetrievedProfile, /dim)
          omiprofilelevelaltitude = reform(profozresult.ProfileLevelAltitude, [sz[0]+1, sz[1] * sz[2]])
          omialts = [[omialts], [omiProfileLevelAltitude[*, roiidx]]]

          omipres = profozresult.ProfileLevelPressure
          omipressize = size(omipres, /dim)
          omipres = reform(omipres, [omipressize[0], omipressize[1]*omipressize[2]])
          omipress = [[omipress], [omipres[*, roiidx]]]

          omigroundpqfs   =  [omigroundpqfs, profozresult.GroundPixelQualityFlags[roiidx]]
          omilats  =  [omilats, profozresult.Latitude[roiidx] ]
          omilons  =  [omilons, profozresult.Longitude[roiidx] ]
          ;omisaas = [omisaas, profozresult.SolarAzimuthAngle[roiidx]]
          omiszas = [omiszas, profozresult.SolarZenithAngle[roiidx]]
          omitimes = [omitimes, profozresult.time[roiidx]]
          ;omivaas = [omivaas, profozresult.ViewingAzimuthAngle[roiidx]]
          ;omivzas = [omivzas, profozresult.ViewingZenithAngle[roiidx]]
          ;omiaods = [omiaods, profozresult.AerosolOpticalThickness[roiidx]]
          ;omiaprioricovs = [omiaprioricovs, profozresult.APrioricov[roiidx]]
          o3ak = reform(profozresult.O3AveragingKernel, [24, 24, long(sz[1]) * sz[2]])
          omiaks = [[[omiaks]], [[o3ak[*, *, roiidx]]]]
          omiecf =  [omiecf, profozresult.EffectiveCloudFraction[roiidx]     ]
          omiecp =  [omiecp, profozresult.EffectiveCloudPressure[roiidx]     ]
          ;omitozs  =  [omitozs, profozresult.ColumnAmountO3[roiidx] ]
          ;omidfss = [omidfss, profozresult.DegreesOfFreedomForSignal[roiidx]]
          ;omiecf1s = [omiecf1s, profozresult.EffectiveCloudFractionUV1[roiidx]]
          ;omiecf2s = [omiecf2s, profozresult.EffectiveCloudFractionUV2[roiidx]]
          ;omiicis = [omiicis, profozresult.InstrumentConfigurationId[roiidx]]
          mqf = profozresult.MeasurementQualityFlags
          mqf = rebin(mqf, sz[2], sz[1])
          mqf = transpose(mqf)
          omimqfs = [omimqfs, mqf[roiidx]]
          ;omino2s = [omino2s, profozresult.NO2ColumnAmount[roiidx]]
          ominumiters = [ominumiters, profozresult.nIteration[roiidx]]
          omio3 = profozresult.O3RetrievedProfile
          nanidx = where(omio3 lt -1.e30, /null)
          omio3[nanidx] = !values.f_nan
          omio3[*,20:23, *] = !values.f_nan
          omio3size = size(omio3, /dim)
          ;omio3 = reform(omio3, [omio3size[0], omio3size[1]*omio3size[2]])
          omio3[where(omio3 lt 0, /null)]=!values.f_nan
          for i = 0, n_elements(roiidx)-1 do begin
            omio3s = [omio3s, omio3[*, roiindices[0, i], roiindices[1, i]]]
          ENDFOR

          omilayero3s = [omilayero3s, omilayero3[roiidx]]
          
          ;omisize = size(omilon, /dim)
          ;omio3 = reform(o3result.o3, [18, omisize[0]*omisize[1]]) 
          ;omio3s = [omio3s, profozresult.O3[roiidx]]

          ;omirefcfs = [omirefcfs, profozresult.ReflectanceCostFunction[roiidx]]

          ;omirowanomalys = [omirowanomalys, omo3prresult.rowanomaly[roiidx]]
          
          ;omixtracks  =  [omixtracks,profozresult.xtracks[roiidx]]
          ;omilines =  [omilines,profozresult.lines[roiidx]]

          ;endif
        ENDIF
      ENDIF  ; ntmp
    ENDFOR  ; ifile


    ;omivars={lon:omilons, lat:omilats, time:omitimes, $
      ;;toz:omitozs,$
      ;SolarZenithAngle:omiszas, $
      ;;SolarAzimuthAngle:omisaas, $
      ;;ReflectanceCostFunction:omirefcfs, $
      ;EffectiveCloudFraction:omiecf, $
      ;EffectiveCloudPressure:omiecp, $
      ;omilayero3:omilayero3, $
      ;;pixlat:pixlats, pixlon:pixlons, pixs:omixtracks, lines:omilines,$
      ;;xtracks:omixtracks, $
      ;;lines:omilines, $
      ;;rowanomaly:omirowanomalys, $
      ;;ref331:omir331s, $
      ;;cf:omicfs, $
      ;;refcf:omirefcfs, $
      ;;pcld:omipclds, $
      ;qf:omigroundpqfs};,bloz:omiblozs}


    ; collocation start

    ;omilon = omivars.Lon
    ;omilat = omivars.Lat
    ;omitime = omivars.Time
    ;;omitoz = omivars.toz
    ;;omilayero3 = omivars.omilayero3
    ;;omipixlat = omivars.pixlat
    ;;omipixlon = omivars.pixlon
    ;;omixtracks = omivars.xtracks
    ;;omilines = omivars.lines
    ;omiecf = omivars.EffectiveCloudFraction
    ;omiecp = omivars.EffectiveCloudPressure
    ;omiqf = omivars.qf
    ;omirowanomaly = omivars.rowAnomaly


    nanidx = where(omitimes lt -5.0E9, /null)
    omitimes[nanidx] = !values.f_nan
    caldat, julday(1,1,1993, 0, 0) + omitimes/60./60./24., omimons, omidays, $
      omiyears, omihours, omiminutes, omiseconds

    ;omimons[nanidx] = 1
    ;omidays[nanidx] = 1
    ;omiyears[nanidx] = 1
    ;omihours[nanidx] = 0
    ;omiminutes[nanidx] = 0
    ;omiseconds[nanidx] = 0
    omijulday = julday(omimons, omidays, omiyears, omihours, omiminutes, omiseconds) 
    
    omijulday[nanidx] = !values.f_nan
    
    gemslon = gemsvars.Longitude
    gemslat = gemsvars.Latitude
    
    ;ds_get_pixcor, gemslon, gemslat, gemslonpixcor, gemslatpixcor
    gems_size = size(gemslon, /dim)

    ; TODO check the input time from gems_l2_o3p output
    ;gemstime = gemsvars.Time


    ;gemstime_sz = size(gemstime)
    ;if gemstime_sz[0] ne 1 or (gemstime_sz[0] eq 1 and gemstime_sz[1] ne 174) then begin 
      ;print, 'size of the variable gemstime is not matched'
      ;gemsvals = []
      ;omivals = []
      ;return
    ;ENDIF
    ;#TODO time variable of gems l2o3p output file is -1.0E30 

    gemstime = fltarr(174)
    gemstime[*] = (julday(month, day, year, hour, minute) - $
      julday(1,1,2000,12,0))*24.*60.*60

    gemstime = rebin(gemstime, [gems_size[0], gems_size[1]])
    ;gemstime = transpose(gemstime)
    ;gemslon[where(gemslon lt -360, /null)] = !values.f_nan
    ;gemslat[where(gemslat lt -90, /null)] = !values.f_nan
    gemslon = gemsvars.longitude
    gemslat = gemsvars.latitude
    ;gemstime[where(gemstime lt 0, /null)] = !values.f_nan
    ;gemstime = julday(month, day, year, hour, minute)
    gemssize = size(gemslon, /dim)
    
    nanidx = where(gemstime lt 0, /null)
    gemstime[nanidx] = !values.f_nan

    caldat, julday(1,1,2000, 12, 0) + gemstime/60./60./24., $
      gemsmons, gemsdays, gemsyears, gemshours, gemsminutes, gemsseconds

    gemsmons[nanidx] = -1
    gemsdays[nanidx] = -1
    gemsyears[nanidx] = -1
    gemshours[nanidx] = -1
    gemsminutes[ nanidx] = -1
    gemsseconds[nanidx] = !values.f_nan

    gemsjulday = julday(gemsmons, gemsdays, gemsyears, gemshours, gemsminutes, gemsseconds)
    gemsjulday[nanidx] = !values.f_nan

    ;closest_idx = lonarr(gemssize[0], gemssize[1])
    ;closest_idx[*] = -999
    ;omisize = size(omilons, /dim)

    ncolloc = n_elements(omitimes)

    closest_idx = lonarr(ncolloc)
    closest_idx[*] = -999

    yn_omi_cross_time = intarr(ncolloc)
    yn_omi_cross_time[*] = 0

    for ip = 0, ncolloc-1 do BEGIN
      x = omilons[ip]
      y = omilats[ip]

      if finite(x) eq 1  and finite(y) eq 1 then begin
        result = search_closest_pixel(gemslon, gemslat, x, y)
        ;result = array_indices(gemslon, result)
        closest_idx[ip] = result
      ENDIF

      ; time difference check
      if closest_idx[ip] ge 0 then begin
        if abs(omijulday[ip] - gemsjulday[closest_idx[ip]]) lt 1./24./2. then BEGIN
          yn_omi_cross_time[ip] = 1
        ENDIF
      endif

    ENDFOR
    
    print, total(yn_omi_cross_time)

    ;save, filename=savfile, gemsvars, gemslayero3, omivars, omilayero3, closest_idx, yn_omi_cross_time
  endif else BEGIN
    RESTORE, savfile
  ENDELSE

  ;gemstoz = reform(gemsvars.ColumnAmountO3[0, *, *])
  gemsvals = []
  omivals = []
  ;omisize = size(omivars.lon, /dim)
  gemsecf = gemsvars.effectivecloudfractionuv
  gemssza = gemsvars.SolarZenithAngle
  gemssza[where(gemssza lt -1E29 , /null)] = !values.f_nan
  nanidx = where(gemsecf lt 0 , /null)
  gemsecf[nanidx] = !values.f_nan
  gemsidx = []

  for ip=0, ncolloc-1 do BEGIN

    if closest_idx[ip] GE 0 and yn_omi_cross_time[ip] eq 1 then begin
      ;if omivars.ReflectanceCostFunction[ix] lt 30 $
          ;and omivars.SolarZenithAngle[ix] lt 5000 $
      if omiszas[ip] lt 90. $
          and xtrackindices[ip] ge 3 $
          and xtrackindices[ip] le 25 $
          ;and omivars.EffectiveCloudFractionUV1[ix] le 0.2 $
            and gemsecf[closest_idx[ip]] le 0.2 $
            and gemsecf[closest_idx[ip]] gt 0 $
            and gemssza[closest_idx[ip]] lt 90. $
            ;and abs(gemssza[closest_idx[ix, iy]] - omivars.SolarZenithAngle[ix, iy]) lt 10 $
            then begin

        if pressure_limit ne 0 or height_limit ne 0 then begin
          omivals = [omivals, omilayero3s[ip]]
          gemsvals = [gemsvals , gemslayero3[closest_idx[ip]]]
          gemsidx = [gemsidx, closest_idx[ip]]
        endif else begin
          omivals = [omivals, omilayero3s[ip]]
          gemsvals = [gemsvals, gemstoz[closest_idx[ip]]]
          gemsidx = [gemsidx, closest_idx[ip]]
        endelse

      endif
    endif

  ENDFOR

  print, n_elements(gemsvals)

  _omival = fltarr(174, 512)
  _omival[*] = !values.f_nan
  _omival[gemsidx] = omivals

  plot_gems_satproj_data, gemslon, gemslat, _omival, $
    filename='./plot/omi_on_gems_grid_'+yyyy+mm+dd+'T'+hh+mi+'.png', $
    title='OMI collocated on the GEMS pixels', range=[20, 60]

  plot_gems_satproj_data, gemslon, gemslat, gemslayero3, $
    filename='./plot/gems_'+yyyy+mm+dd+'T'+hh+mi + '.png', $
    title='GEMS', range=[20, 60]

  ;plot_omi_satproj, omivars.lon, omivars.lat, omivars.rowanomaly,$
    ;filename='plot/OMI_rowanomaly_20200806T0345.png', $
    ;title='OMI Row Anomaly', range=[0, 1.]
  return
endif else begin
ENDELSE

end
