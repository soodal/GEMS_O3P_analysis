;+
; This procedure collocate GEMS L2 O3P with OMI PROFOZ on the GEMS GRID.
;-
;

pro collocate_gemsl2o3p_omil2profoz_on_gems_grid, year, month, day, hour, minute, $
  output, gemsfile=gemsfile

if not keyword_Set(pressure_limit) then begin
  pressure_limit = 0
endif

if not keyword_Set(height_limit) then begin
  height_limit = 0
endif

limit=[-10, $ ;minimum latitude
  60, $ ;minimum longitude
  60, $ ;maximum latitude
  160] ;maximum longitude

yyyy = string(year, format='(i04)')
mm = string(month, format='(i02)')
dd = string(day, format='(i02)')
hh = string(hour, format='(i02)')
mi = string(minute, format='(i02)')

savpath = './collocate_gemsl2o3p_profoz/'
savfile = savpath + 'col_omi_on_gems_grid_'+yyyy+mm+dd+'_'+hh+mi+'.sav'

savefile =  file_test(savfile) 

savefile = 0
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
  omio3s = []
  omio3aps = []
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

  ; GEMS file name set
  if keyword_Set(gemsfile) then begin
    gemsfn = gemsfile
  endif else begin
    if year eq 2020 and month eq 06 and day eq 16 then begin
      gemsfn = '/data2/L2_GEMS/val_1008/GK2_GEMS_O3P_' $
        + yyyy + mm + dd + '_' + hh + mi + '.nc4'
    endif else begin 
      gemsfn = ds_get_gems_o3p_filename(yyyy, mm, dd, hh, mi)
    endelse
  endelse

  if gemsfn ne !null then begin

    ; read GEMS O3P data

    print, gemsfn

    varlist = ['EffectiveCloudFractionUV', 'ProcessingQualityFlags', $
             'O3' , $
             'O3Apriori', $
             ;'O3AprioriError', $
             ;'CloudPressure', $
             ;'O3Apriori', 'O3AprioriError',$
             'ColumnAmountO3', $
             'SolarZenithAngle', $
             ;'SimulatedRadiances', $
             'Latitude' ,'Longitude', $
             'Time','Altitude' ,    $
             'Pressure', $
             ;'TropopausePressure', $
             ;'Wavelengths', $
             'WavelengthsWholeRange']
    gemsvars = ds_read_gems_l2_o3p(gemsfn, varlist=varlist)

    gems_o3 = gemsvars.o3
    gems_o3ap = gemsvars.O3Apriori
    gems_o3_size = size(gems_o3, /dimension)

    ; create output variables

    if gems_o3_size[0] eq 174 then begin
      colloc_omi_o3p_on_gems = fltarr(gems_o3_size[0], gems_o3_size[1], gems_o3_size[2])
      colloc_omi_o3p_on_gems[*] = !values.f_nan
      colloc_omi_xtrack_on_gems = fltarr(gems_o3_size[0], gems_o3_size[1], gems_o3_size[2])
      colloc_omi_xtrack_on_gems[*] = !values.f_nan
      colloc_omi_atrack_on_gems = fltarr(gems_o3_size[0], gems_o3_size[1], gems_o3_size[2])
      colloc_omi_atrack_on_gems[*] = !values.f_nan
      colloc_omi_alt_on_gems = fltarr(gems_o3_size[0], gems_o3_size[1], gems_o3_size[2]+1)
      colloc_omi_alt_on_gems[*] = !values.f_nan
      colloc_omi_pres_on_gems = fltarr(gems_o3_size[0], gems_o3_size[1], gems_o3_size[2]+1)
      colloc_omi_pres_on_gems[*] = !values.f_nan
      colloc_omi_ecf_on_gems = fltarr(gems_o3_size[0], gems_o3_size[1])
      colloc_omi_ecf_on_gems[*] = !values.f_nan
    ENDIF else if gems_o3_size[0] eq 24 then begin
      colloc_omi_o3p_on_gems = fltarr(gems_o3_size[2], gems_o3_size[0], gems_o3_size[1])
      colloc_omi_o3p_on_gems[*] = !values.f_nan
      colloc_omi_xtrack_on_gems = fltarr(gems_o3_size[2], gems_o3_size[0], gems_o3_size[1])
      colloc_omi_xtrack_on_gems[*] = !values.f_nan
      colloc_omi_atrack_on_gems = fltarr(gems_o3_size[2], gems_o3_size[0], gems_o3_size[1])
      colloc_omi_atrack_on_gems[*] = !values.f_nan
      colloc_omi_alt_on_gems = fltarr(gems_o3_size[2]+1, gems_o3_size[0], gems_o3_size[1])
      colloc_omi_alt_on_gems[*] = !values.f_nan
      colloc_omi_pres_on_gems = fltarr(gems_o3_size[2]+1, gems_o3_size[0], gems_o3_size[1])
      colloc_omi_pres_on_gems[*] = !values.f_nan
      colloc_omi_ecf_on_gems = fltarr(gems_o3_size[0], gems_o3_size[1])
      colloc_omi_ecf_on_gems[*] = !values.f_nan
    ENDIF
    ; path for OMI

    omipath = '/data2/OMI/gdata/' + yyyy + '/L2-PROFOZ/' 

    omifiles=FILE_SEARCH(omipath+'OMI-Aura_L2-PROFOZ_'+yyyy+'m'+mm+dd+'*.he5')
    ;omifiles = omifiles[$
      ;where($
        ;strmatch(omifiles, '*OMI-Aura_L2-PROFOZ_'+yyyy+'m'+mm+dd+'t0[1-6]*.he5') eq 1)]

    ; read OMI data

    nfiles = n_elements(omifiles)

    for ifile=0, nfiles-1 do begin
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
          ;ds_omi_l2_profoz_accum, profozresult, omilayero3, height=10.

          omidim = size(profozresult.EffectiveCloudFraction, /dim)

          omixtrackidx = rebin(indgen(omidim[0]), omidim[0], omidim[1])
          omialongtrackidx = rebin($
            transpose(indgen(omidim[1])),  omidim[0], omidim[1])
          omilatpixcor = fltarr(omidim[0], omidim[1], 4)
          omilonpixcor = fltarr(omidim[0], omidim[1], 4)

          for ix=0, omidim[0]-1 do begin
            for iy=0, omidim[1]-1 do begin
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

          ;time=STRMID(omifiles[ifile],46,14,/rev)

          omio3sz = size(profozresult.O3RetrievedProfile, /dim)
          omiprofilelevelaltitude = reform($
            profozresult.ProfileLevelAltitude, $
            [omio3sz[0]+1, omio3sz[1] * omio3sz[2]])

          omialts = [[omialts], [omiProfileLevelAltitude[*, roiidx]]]

          omipres = profozresult.ProfileLevelPressure
          omipressize = size(omipres, /dim)
          omipres = reform(omipres, [omipressize[0], omipressize[1]*omipressize[2]])
          omipress = [[omipress], [omipres[*, roiidx]]]

          omigroundpqfs = [omigroundpqfs, profozresult.GroundPixelQualityFlags[roiidx]]
          omilats = [omilats, profozresult.Latitude[roiidx] ]
          omilons = [omilons, profozresult.Longitude[roiidx] ]
          ;omisaas = [omisaas, profozresult.SolarAzimuthAngle[roiidx]]
          omiszas = [omiszas, profozresult.SolarZenithAngle[roiidx]]

          profoztime = fltarr(omio3sz[1], omio3sz[2])
          profoztime[*] = !values.f_nan
          for ix=0, omio3sz[1]-1 do begin
            profoztime[ix, *] =profozresult.time 
          ENDFOR
          
          omitimes = [omitimes, profoztime[roiidx]]
          ;omivaas = [omivaas, profozresult.ViewingAzimuthAngle[roiidx]]
          ;omivzas = [omivzas, profozresult.ViewingZenithAngle[roiidx]]
          ;omiaods = [omiaods, profozresult.AerosolOpticalThickness[roiidx]]
          ;omiaprioricovs = [omiaprioricovs, profozresult.APrioricov[roiidx]]

          o3ak = reform(profozresult.O3AveragingKernel, $
            [24, 24, long(omio3sz[1]) * omio3sz[2]])
          omiaks = [[[omiaks]], [[o3ak[*, *, roiidx]]]]

          omiecf =  [omiecf, profozresult.EffectiveCloudFraction[roiidx]]
          omiecp =  [omiecp, profozresult.EffectiveCloudPressure[roiidx]]
          ;omitozs  =  [omitozs, profozresult.ColumnAmountO3[roiidx] ]
          ;omidfss = [omidfss, profozresult.DegreesOfFreedomForSignal[roiidx]]
          ;omiecf1s = [omiecf1s, profozresult.EffectiveCloudFractionUV1[roiidx]]
          ;omiecf2s = [omiecf2s, profozresult.EffectiveCloudFractionUV2[roiidx]]
          ;omiicis = [omiicis, profozresult.InstrumentConfigurationId[roiidx]]

          mqf = profozresult.MeasurementQualityFlags
          mqf = rebin(mqf, omio3sz[2], omio3sz[1])
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
            omio3s = [[omio3s], $
              [reform(omio3[*, roiindices[0, i], roiindices[1, i]])]]
          ENDFOR

          omio3ap = profozresult.O3Aprioriprofile
          nanidx = where(omio3ap lt -1.e30, /null)
          omio3ap[nanidx] = !values.f_nan
          omio3ap[*,20:23, *] = !values.f_nan
          omio3apsize = size(omio3ap, /dim)
          omio3ap[where(omio3ap lt 0, /null)]=!values.f_nan

          for i = 0, n_elements(roiidx)-1 do begin
            omio3aps = [[omio3aps], $
              [reform(omio3ap[*, roiindices[0, i], roiindices[1, i]])]]
          ENDFOR
          
        ENDIF
      ENDIF  ; ntmp
    ENDFOR  ; ifile

; collocation start

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
    
    ;ds_get_pixcor, gemsvars.longitude, gemsvars.latitude, gemslonpixcor, gemslatpixcor
    gemssize = size(gemsvars.longitude, /dim)

    ; TODO check the input time from gems_l2_o3p output
    gemstime = gemsvars.Time

    gemstime_sz = size(gemstime)
    if gemstime_sz[0] eq 1 then begin
      if gemstime_sz[1] eq 512 then BEGIN
        print, 'Size of the variable gemstime is 512.'
        print, 'May be this is for the date of 20200616'
        print, 'Actual values are in [0:173] of the varaible'
        gemstime = gemstime[0:173]
      endif else if gemstime_sz[1] ne 174 then begin 
        print, 'size of the variable gemstime is not matched'
        gems_o3s = []
        omi_o3s = []
        stop
        return
      ENDIF
    ENDIF

    ;TODO time variable of gems l2o3p output file is -1.0E30 

    ;gemstime = fltarr(174)
    ;gemstime[*] = (julday(month, day, year, hour, minute) - $
      ;julday(1,1,2000,12,0))*24.*60.*60
    gemstime = rebin(gemstime, [gemssize[0], gemssize[1]])
    
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

    ;close2gemspixidx = lonarr(gemssize[0], gemssize[1])
    ;close2gemspixidx[*] = -999
    ;omisize = size(omilons, /dim)

    ncolloc = n_elements(omitimes)

    close2gemspixidx = lonarr(ncolloc)
    close2gemspixidx[*] = -999

    yn_omi_cross_time = intarr(ncolloc)
    yn_omi_cross_time[*] = 0
    if gems_o3_size[0] eq 174 then BEGIN
      yn_omi_cross_time_on_gems = intarr( gems_o3_size[0], gems_o3_size[1])
      yn_omi_cross_time_on_gems[*] = 0
    ENDIF else if gems_o3_size[0] eq 24 then BEGIN
      yn_omi_cross_time_on_gems = intarr( gems_o3_size[1], gems_o3_size[2])
      yn_omi_cross_time_on_gems[*] = 0
    ENDIF

    for ip = 0, ncolloc-1 do BEGIN
      x = omilons[ip]
      y = omilats[ip]

      if finite(x) eq 1  and finite(y) eq 1 then begin
        result = search_closest_pixel(gemsvars.longitude, gemsvars.latitude, x, y)
        if result ge 0 then begin 
          result2 = array_indices(gemsvars.longitude, result)
          close2gemspixidx[ip] = result
          

          if gems_o3_size[0] eq 174 then begin 
            colloc_omi_o3p_on_gems[result2[0], result2[1], *] =  omio3s[*, ip]
          ENDIF else if gems_o3_size[0] eq 24 then begin
            colloc_omi_o3p_on_gems[*, result2[0], result2[1]] =  omio3s[*, ip]
          ENDIF

          if gems_o3_size[0] eq 174 then begin 
            colloc_omi_alt_on_gems[result2[0], result2[1], *] =  omialts[*, ip]
          ENDIF else if gems_o3_size[0] eq 24 then begin
            colloc_omi_alt_on_gems[*, result2[0], result2[1]] =  omialts[*, ip]
          ENDIF

          if gems_o3_size[0] eq 174 then begin 
            colloc_omi_pres_on_gems[result2[0], result2[1], *] =  omipress[*, ip]
          ENDIF else if gems_o3_size[0] eq 24 then begin
            colloc_omi_pres_on_gems[*, result2[0], result2[1]] =  omipress[*, ip]
          ENDIF

          if gems_o3_size[0] eq 174 then begin 
            colloc_omi_ecf_on_gems[result2[0], result2[1], *] =  omiecf[ip]
          ENDIF else if gems_o3_size[0] eq 24 then begin
            colloc_omi_ecf_on_gems[*, result2[0], result2[1]] =  omiecf[ip]
          ENDIF


          ; time difference check

          if abs(omijulday[ip] - gemsjulday[close2gemspixidx[ip]]) lt 1./24./2. then BEGIN
            yn_omi_cross_time[ip] = 1
            yn_omi_cross_time_on_gems[result2[0], result2[1]] = 1
          ENDIF
        ENDIF
      ENDIF
    ENDFOR
  endif

  ;gemstoz = reform(gemsvars.ColumnAmountO3[0, *, *])
  gems_o3s = []
  gems_o3aps = []

  omi_o3s = []
  omi_o3aps = []
  omi_idx = []
  omi_effectivecloudfraction = []
  omi_solarzenithangle = []
  ;omisize = size(omivars.lon, /dim)

  gems_ecf = gemsvars.effectivecloudfractionuv
  gems_sza = gemsvars.SolarZenithAngle
  gems_alts = []
  gems_press = []
  gems_effectivecloudfraction = []
  gems_solarzenithangle = []
  gems_latitude = []
  gems_longitude = []

  gems_sza[where(gems_sza lt -1E29 , /null)] = !values.f_nan
  nanidx = where(gems_ecf lt 0 , /null)
  gems_ecf[nanidx] = !values.f_nan
  gems_idx = []
  gems_altitude = gemsvars.Altitude
  gems_pressure = gemsvars.Pressure

  valididx = where(close2gemspixidx ge 0, /null)
  ; close2gemspixindices has different origin with close2gemspixidx.
  close2gemspixindices = array_indices(gemsvars.longitude, close2gemspixidx[valididx])
  valididx_sz = size(valididx, /dim)

  for ip=0, valididx_sz[0]-1 do BEGIN
    if close2gemspixidx[valididx[ip]] GE 0 and yn_omi_cross_time[valididx[ip]] eq 1 then begin
      ;if omivars.ReflectanceCostFunction[ix] lt 30 $
          ;and omivars.SolarZenithAngle[ix] lt 5000 $
      if omiszas[valididx[ip]] lt 90. $
          and xtrackindices[valididx[ip]] ge 3 $
          and xtrackindices[valididx[ip]] le 25 $
          ;and omivars.EffectiveCloudFractionUV1[ix] le 0.2 $
            ;and gems_ecf[close2gemspixidx[valididx[ip]]] le 0.2 $
            ;and gems_ecf[close2gemspixidx[valididx[ip]]] gt 0 $
            and gems_sza[close2gemspixidx[valididx[ip]]] lt 90. $
            ;and abs(gems_sza[close2gemspixidx[ix, iy]] - omivars.SolarZenithAngle[ix, iy]) lt 10 $
            then begin

        omi_o3s = [[omi_o3s], [reform(omio3s[*, valididx[ip]])]]
        omi_o3aps = [[omi_o3aps], [reform(omio3aps[*, valididx[ip]])]]
        omi_idx = [omi_idx, valididx[ip]]
        omi_effectivecloudfraction = [omi_effectivecloudfraction, $
          omiecf[valididx[ip]]]
        omi_solarzenithangle = [omi_solarzenithangle, omiszas[valididx[ip]]]


      endif
    endif

  ENDFOR

  if n_elements(where(close2gemspixidx ge 0 and $
      yn_omi_cross_time eq 1, /null)) gt 0 then begin
    output = create_struct('gemsl2o3p', gemsvars $
      , 'colloc_omi_o3p_on_gems', colloc_omi_o3p_on_gems $
      , 'colloc_omi_xtrack_on_gems', colloc_omi_xtrack_on_gems $
      , 'colloc_omi_atrack_on_gems', colloc_omi_atrack_on_gems $
      , 'colloc_omi_alt_on_gems', colloc_omi_alt_on_gems $
      , 'colloc_omi_pres_on_gems', colloc_omi_pres_on_gems $
      , 'colloc_omi_ecf_on_gems', colloc_omi_ecf_on_gems $
      )
  endif else begin
    output = !null
  ENDELSE

  ;_omival = fltarr(174, 512)
  ;_omival[*] = !values.f_nan
  ;_omival[gems_idx] = omio3s

  ;plot_gems_satproj_data, gemsvars.longitude, gemsvars.latitude, _omival, $
    ;filename='./plot/omi_on_gems_grid_'+yyyy+mm+dd+'T'+hh+mi+'.png', $
    ;title='OMI collocated on the GEMS pixels', range=[20, 60]

  ;plot_gems_satproj_data, gemsvars.longitude, gemsvars.latitude, gemslayero3, $
    ;filename='./plot/gems_'+yyyy+mm+dd+'T'+hh+mi + '.png', $
    ;title='GEMS', range=[20, 60]

  ;plot_omi_satproj, omivars.lon, omivars.lat, omivars.rowanomaly,$
    ;filename='plot/OMI_rowanomaly_20200806T0345.png', $
    ;title='OMI Row Anomaly', range=[0, 1.]
  save, filename=savfile, output
endif else begin
  RESTORE, savfile
ENDELSE
return

end
