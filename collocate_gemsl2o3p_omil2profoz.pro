;+
; This procedure collocate GEMS L2 O3P with OMI PROFOZ on the GEMS GRID.
;-
pro collocate_gemsl2o3p_omil2profoz, year, month, day, hour, minute, $
  gemsout, omiout, gemsfile=gemsfile

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

savpath = '/data/private/soodal/collocated_gemsl2o3p_omil2profoz/'
savfile = savpath + 'col_gems_omi_'+yyyy+mm+dd+'_'+hh+mi+'.sav'

isSaved =  file_test(savfile) 

if not isSaved then BEGIN

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

  omio3layers = []
  omio3aps = []
  xtrackindices = []
  alongtrackindices = []
  glintprobabilities= []
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

  ;omipath = '~/data/OMPROFOZ/' + yyyy + '/' + mm + '/' + dd + '/'
  omipath = '/data/SAT/OMI/OMPROFOZ/' +  yyyy + '/' + mm + '/' + dd + '/'
  omifiles=FILE_SEARCH(omipath+'OMI-Aura_L2-PROFOZ_'+yyyy+'m'+mm+dd+'*.he5')
  ;omifiles = omifiles[$
    ;where($
      ;strmatch(omifiles, '*OMI-Aura_L2-PROFOZ_'+yyyy+'m'+mm+dd+'t0[1-6]*.he5') eq 1)]

  ; read OMI data
  nfiles = n_elements(omifiles)
  if nfiles eq 1 and strlen(omifiles[0]) eq 0 then begin
    ;message, 'omi he5 file is not found at ', omipath
    print, 'omi he5 file is not found at ', omipath
    omiout = []
    gemsout = []
    return
  endif



  if gemsfn ne !null and (not (nfiles eq 1 and strlen(omifiles[0]))) then begin

    ; read GEMS O3P data
    ;print, gemsfn

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

    ds_gems_l2o3p_accum, gemsvars, gemso3accum, height=10


    omi_o3p_on_gems_grid = fltarr(gems_o3_size[0], gems_o3_size[1], gems_o3_size[2])
    omi_o3p_on_gems_grid[*] = !values.f_nan
    omi_xtrack_on_gems_grid = fltarr(gems_o3_size[0], gems_o3_size[1], gems_o3_size[2])
    omi_xtrack_on_gems_grid[*] = !values.f_nan
    omi_atrack_on_gems_grid = fltarr(gems_o3_size[0], gems_o3_size[1], gems_o3_size[2])
    omi_atrack_on_gems_grid[*] = !values.f_nan
    omi_alt_on_gems_grid = fltarr(gems_o3_size[0], gems_o3_size[1], gems_o3_size[2])
    omi_alt_on_gems_grid[*] = !values.f_nan
    omi_pres_on_gems_grid = fltarr(gems_o3_size[0], gems_o3_size[1], gems_o3_size[2])
    omi_pres_on_gems_grid[*] = !values.f_nan
    ;omi_on_gems_grid = fltarr(gems_o3_size[0], gems_o3_size[1], gems_o3_size[2])
    ;omi_on_gems_grid[*] = !values.f_nan
    ;omi_on_gems_grid = fltarr(gems_o3_size[0], gems_o3_size[1], gems_o3_size[2])
    ;omi_on_gems_grid[*] = !values.f_nan
    omi_o3p_accum_on_gems_grid = fltarr(gems_o3_size[0], gems_o3_size[1])
    omi_o3p_accum_on_gems_grid[*] = !values.f_nan

    ; path for OMI

    ;omipath = '~/data/OMPROFOZ/' + yyyy + '/' + mm + '/' + dd + '/'
    ;omipath = '/data/SAT/OMI/OMPROFOZ/' +  yyyy + '/' + mm + '/' + dd + '/'
    omifiles=FILE_SEARCH(omipath+'OMI-Aura_L2-PROFOZ_'+yyyy+'m'+mm+dd+'*.he5')
    ;omifiles = omifiles[$
      ;where($
        ;strmatch(omifiles, '*OMI-Aura_L2-PROFOZ_'+yyyy+'m'+mm+dd+'t0[1-6]*.he5') eq 1)]

    ; read OMI data
    nfiles = n_elements(omifiles)
    if nfiles eq 1 and strlen(omifiles[0]) eq 0 then begin
      ;message, 'omi he5 file is not found at ', omipath
      print, 'omi he5 file is not found at ', omipath
      omiout = []
      gemsout = []
      return
    endif

    for ifile=0, nfiles-1 do begin
      fi = file_info(omifiles[ifile])
      omiobstime = strmid(omifiles[ifile], 46, 14, /reverse)
      omiorbitnumb = strmid(omifiles[ifile], 31, 6, /reverse)
      if fi.size gt 1500 then begin
        ;print, omifiles[ifile]
        ds_READ_omi_L2_profoz, omifiles[ifile], profozresult


        roiidx=WHERE(profozresult.longitude GE limit[1] AND $
                  profozresult.longitude LE limit[3] AND $
                  profozresult.latitude GE limit[0] AND $
                  profozresult.latitude LE limit[2] and $
                  profozresult.EffectiveCloudFraction LE 0.2, ntmp, /null)
        print, ntmp

        if ntmp ge 1 then begin 
          roiindices = array_indices(profozresult.latitude, roiidx)
          ds_omi_l2_profoz_accum, profozresult, omio3accum, height=10.

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
          glintprobabilities = [glintprobabilities, profozresult.glintprobability[roiidx]]


          omio3sz = size(profozresult.O3RetrievedProfile, /dim)
          omiprofilelevelaltitude = reform($
            profozresult.ProfileLevelAltitude, $
            [omio3sz[0]+1, omio3sz[1] * omio3sz[2]])

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

          profoztime = fltarr(omio3sz[1], omio3sz[2])
          profoztime[*] = !values.f_nan
          for ix=0, omio3sz[1]-1 do begin
            profoztime[ix, *] =profozresult.time 
          ENDFOR
          
          omitimes = [omitimes, profoztime[roiidx]]
          o3ak = reform(profozresult.O3AveragingKernel, $
            [24, 24, long(omio3sz[1]) * omio3sz[2]])
          omiaks = [[[omiaks]], [[o3ak[*, *, roiidx]]]]
          omiecf =  [omiecf, profozresult.EffectiveCloudFraction[roiidx]]
          omiecp =  [omiecp, profozresult.EffectiveCloudPressure[roiidx]]
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

            omio3layers = [[omio3layers], omio3accum[roiindices[0, i], roiindices[1, i]]]
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

    if n_elements(omitimes) gt 0 then begin
      nanidx = where(omitimes lt -5.0E9, /null)
      omitimes[nanidx] = !values.f_nan
      caldat, julday(1,1,1993, 0, 0) + omitimes/60./60./24., omimons, omidays, $
        omiyears, omihours, omiminutes, omiseconds

      omijulday = julday(omimons, omidays, omiyears, omihours, omiminutes, omiseconds) 
      
      omijulday[nanidx] = !values.f_nan
      
      gemslon = gemsvars.Longitude
      gemslat = gemsvars.Latitude
      
      sz = size(gemsvars.solarzenithangle, /dimension)
      nx = sz[0]

      ;ds_get_pixcor, gemslon, gemslat, gemslonpixcor, gemslatpixcor
      gemssize = size(gemslon, /dim)

      ; TODO check the input time from gems_l2_o3p output
      gemstime = gemsvars.Time

      gemstime_sz = size(gemstime)
      if gemstime_sz[0] eq 1 then begin
        if gemstime_sz[1] eq 512 then BEGIN
          print, 'Size of the variable gemstime is 512.'
          print, 'May be this is for the date of 20200616'
          print, 'Actual values are in [0:173] of the varaible'
          gemstime = gemstime[0:173]
        endif else if gemstime_sz[1] ne nx then begin 
          print, 'size of the variable gemstime is not matched'
          gems_o3s = []
          omi_o3s = []
          stop
          return
        ENDIF
      ENDIF

      ;#TODO time variable of gems l2o3p output file is -1.0E30 

      gemstime = fltarr(nx)
      gemstime[*] = (julday(month, day, year, hour, minute) - $
        julday(1,1,2000,12,0))*24.*60.*60
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

      omiongemsnum = n_elements(omitimes)

      gemsindices = lonarr(omiongemsnum)
      gemsindices[*] = -999

      yn_omi_cross_time = intarr(omiongemsnum)
      yn_omi_cross_time[*] = 0

      ; loop for every collocation point in a scene for the GEMS
      for ipix = 0, omiongemsnum-1 do BEGIN
        x = omilons[ipix]
        y = omilats[ipix]

        if finite(x) eq 1  and finite(y) eq 1 then begin
          gemsidx = search_closest_pixel(gemslon, gemslat, x, y)
          if gemsidx ge 0 then begin
            indices = array_indices(gemslon, gemsidx)
            gemsindices[ipix] = gemsidx
            omi_o3p_on_gems_grid[indices[0], indices[1], *] =  omio3s[*, ipix]

            omi_o3p_accum_on_gems_grid[indices[0], indices[1]] = omio3layers[ipix]
          endif
        ENDIF

        ; time difference check
        if gemsindices[ipix] ge 0 then begin
          if abs(omijulday[ipix] - gemsjulday[gemsindices[ipix]]) lt 1./24./2. then BEGIN
            yn_omi_cross_time[ipix] = 1
          ENDIF
        endif
      ENDFOR ; ipix
    endif
  endif

  gems_o3s = []
  gems_o3aps = []

  omi_o3s = []
  omi_o3layers = []

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

  colgemsidx = where(gemsindices ge 0, /null)
  ; close2gemspixindices has different origin with gemsindices.
  close2gemspixindices = array_indices(gemslon, gemsindices[colgemsidx])
  colgemsidx_sz = size(colgemsidx, /dim)

  for ipix=0, colgemsidx_sz[0]-1 do BEGIN
    if gemsindices[colgemsidx[ipix]] GE 0 and yn_omi_cross_time[colgemsidx[ipix]] eq 1 then begin
      ;if omivars.ReflectanceCostFunction[ix] lt 30 $
          ;and omivars.SolarZenithAngle[ix] lt 5000 $
      if omiszas[colgemsidx[ipix]] lt 70. $
          and xtrackindices[colgemsidx[ipix]] ge 3 $
          and xtrackindices[colgemsidx[ipix]] le 25 $
          and xtrackindices[colgemsidx[ipix]] ne 22 $
          and glintprobabilities[colgemsidx[ipix]] eq 0.0 $
          ;and omivars.EffectiveCloudFractionUV1[ix] le 0.2 $
            ;and gems_ecf[gemsindices[colgemsidx[ipix]]] le 0.2 $
            ;and gems_ecf[gemsindices[colgemsidx[ipix]]] gt 0 $
            and gems_sza[gemsindices[colgemsidx[ipix]]] lt 70. $
            ;and abs(gems_sza[gemsindices[ix, iy]] - omivars.SolarZenithAngle[ix, iy]) lt 10 $
            then begin

        omi_o3s = [[omi_o3s], [reform(omio3s[*, colgemsidx[ipix]])]]

        omi_o3layers = [[omi_o3layers], omio3accum[colgemsidx[ipix]]]

        omi_o3aps = [[omi_o3aps], [reform(omio3aps[*, colgemsidx[ipix]])]]
        omi_idx = [omi_idx, colgemsidx[ipix]]
        omi_effectivecloudfraction = [omi_effectivecloudfraction, $
          omiecf[colgemsidx[ipix]]]
        omi_solarzenithangle = [omi_solarzenithangle, omiszas[colgemsidx[ipix]]]

        gems_o3s = [[gems_o3s] , $
          [reform(gems_o3[close2gemspixindices[0, ipix], $
            close2gemspixindices[1, ipix], *])]]
        gems_o3aps = [[gems_o3aps] , $
          [reform(gems_o3ap[close2gemspixindices[0, ipix], $
            close2gemspixindices[1, ipix], *])]]

        gems_idx = [gems_idx, gemsindices[colgemsidx[ipix]]]
        gems_alts = [[gems_alts], $
          [reform(gems_altitude[close2gemspixindices[0, ipix], $
            close2gemspixindices[1, ipix], *])]]
        gems_press = [[gems_press], $
          [reform(gems_pressure[close2gemspixindices[0, ipix], $
            close2gemspixindices[1, ipix], *])]]

        gems_effectivecloudfraction = [gems_effectivecloudfraction, $
          gems_ecf[gemsindices[colgemsidx[ipix]]]]

        gems_solarzenithangle = [gems_solarzenithangle, $
          gems_sza[gemsindices[colgemsidx[ipix]]]]

        gems_latitude = [gems_latitude, $
          gemslat[gemsindices[colgemsidx[ipix]]]]

        gems_longitude = [gems_longitude, $
          gemslon[gemsindices[colgemsidx[ipix]]]]

      endif
    endif

  ENDFOR ; ipix

  ;if n_elements(where(gemsindices ge 0 and $
      ;yn_omi_cross_time eq 1, /null)) gt 0 then begin

  if n_elements(omi_o3s) gt 0 then begin
    output = create_struct('omi_o3', omi_o3s, $
      'omi_o3_apriori', omi_o3aps, $
      'omi_pressure', omipress[*, omi_idx], $
      'omi_altitue',omialts[*, omi_idx], $
      'omi_longitude', omilons[omi_idx], $
      'omi_latitude', omilats[omi_idx], $
      'omi_effectivecloudfraction', omi_effectivecloudfraction, $
      'omi_solarzenithangle', omi_solarzenithangle, $
      'gems_o3', gems_o3s, $
      'gems_o3_apriori', gems_o3aps, $
      'gems_altitude', gems_alts, $
      'gems_pressure', gems_press, $
      'gems_effectivecloudfraction', gems_effectivecloudfraction, $
      'gems_solarzenithangle', gems_solarzenithangle, $
      'gems_latitude', gems_latitude, $
      'gems_longitude', gems_longitude, $
      'gems_index', gems_idx)
  endif else begin
    output = !null
  ENDELSE



  ;_omival = fltarr(nx, 512)
  ;_omival[*] = !values.f_nan
  ;_omival[gems_idx] = omio3s

  ;plot_gems_satproj_data, gemslon, gemslat, _omival, $
    ;filename='./plot/omi_on_gems_grid_'+yyyy+mm+dd+'T'+hh+mi+'.png', $
    ;title='OMI collocated on the GEMS pixels', range=[20, 60]

  ;plot_gems_satproj_data, gemslon, gemslat, gemslayero3, $
    ;filename='./plot/gems_'+yyyy+mm+dd+'T'+hh+mi + '.png', $
    ;title='GEMS', range=[20, 60]

  ;plot_omi_satproj, omivars.lon, omivars.lat, omivars.rowanomaly,$
    ;filename='plot/OMI_rowanomaly_20200806T0345.png', $
    ;title='OMI Row Anomaly', range=[0, 1.]
  save, filename=savfile, omi_o3p_accum_on_gems_grid, gemso3accum
endif else begin
  RESTORE, savfile
ENDELSE

colidx = where(finite(omi_o3p_accum_on_gems_grid) eq 1 and finite(gemso3accum) eq 1, /null)

omiout = omi_o3p_accum_on_gems_grid[colidx]
gemsout = gemso3accum[colidx]

return

end
