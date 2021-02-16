pro plot_daily_omi_collocated_with_gems_daily, year, month, day, hour, minute, $
  omivals, $
  hpa=pressure_limit, height=height_limit

if not keyword_Set(pressure_limit) then begin
  pressure_limit = 0
endif

if not keyword_Set(height_limit) then begin
  height_limit = 0
endif

limit=[-10, $ ;minimum latitude
  80, $ ;minimum longitude
  ;60, $ ;maximum latitude
  10, $ ;maximum latitude
  160] ;maximum longitude

yyyy = string(year, format='(i04)')
mm = string(month, format='(i02)')
dd = string(day, format='(i02)')
hh = string(hour, format='(i02)')
mi = string(minute, format='(i02)')

savpath = './collocate_gemsl2o3p_omo3pr/'
savfile = savpath + 'col_gems_omi_'+yyyy+mm+dd+'_'+hh+mi+'.sav'

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
  omio3s = fltarr(24, 1)
  ominumiters = []
  omixtrackidxs = []
  omialongtrackidxs = []

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
;; GEMS vertical column layer for under specific height
;;----------------------------

  ds_gems_l2o3p_accum, gemsvars, gemslayero3, height=10.


  ; path for OMI

  omipath = '/data2/OMI/gdata/' + yyyy + '/L2-PROFOZ/' 
  omipixpath = '/data2/OMI/gdata/' + yyyy + '/L2-OMPIXCOR/' 

  omifiles=FILE_SEARCH(omipath+'OMI-Aura_L2-PROFOZ_'+yyyy+'m'+mm+dd+'*.he5', count=nfiles)
  ; read OMI data
  for ifile=0, nfiles-1 do begin
    fi = file_info(omifiles[ifile])

    omiobstime = strmid(omifiles[ifile], 46, 14, /reverse)
    omiorbitnumb = strmid(omifiles[ifile], 31, 6, /reverse)
    if fi.size gt 1500 then begin
      ds_READ_omi_L2_profoz, omifiles[ifile], profozresult

      roiidx=WHERE(profozresult.longitude GE limit[1] AND $
                profozresult.longitude LE limit[3] AND $
                profozresult.latitude GE limit[0] AND $
                profozresult.latitude LE limit[2] $
                profozresult.time
                , ntmp, /null)

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
          omio3size = size(omio3, /dim)
          omio3 = reform(omio3, [omio3size[0], omio3size[1]*omio3size[2]])
          omio3[where(omio3 lt 0, /null)]=!values.f_nan
          omio3s = [[omio3s], [omio3[*, roiidx]]]
          
          ;omirefcfs = [omirefcfs, profozresult.ReflectanceCostFunction[roiidx]]

          ;omirowanomalys = [omirowanomalys, omo3prresult.rowanomaly[roiidx]]
          
          ;omixtracks  =  [omixtracks,profozresult.xtracks[roiidx]]
          ;omilines =  [omilines,profozresult.lines[roiidx]]

        ;endif
      ENDIF
  ENDFOR  ; ifile


;;----------------------------
;; OMI vertical column Layer for under specific pressure
;;----------------------------
if pressure_limit ne 0 then begin
  ;omitozsize = size(omitozs, /dim)
  omiecfsize = size(omiecf, /dim)
  omilayero3 = fltarr(omiecfsize[0])

  altnanidx = where(omialts lt -1.E30, /null)
  omialts[altnanidx] = !values.f_nan

  for ip=0,omiecfsize[0]-1 do BEGIN
    for ilevel=24, 0, -1 do BEGIN
      ilayer = ilevel - 1
      if omipress[ilevel, ip] gt pressure_limit $
          and omipress[ilevel-1, ip] gt pressure_limit then BEGIN
        omilayero3[ip] = omilayero3[ip] + omio3s[ilayer, ip]
        
      ENDIF else if omipress[ilevel, ip] gt pressure_limit $
        and omipress[ilevel-1, ip] le pressure_limit then begin 
        ; upper bound is higher than 300hPa 

        omilayero3[ip] = omilayero3[ip] + omio3s[ilayer, ip]*$
          (omipress[ilevel, ip]-pressure_limit)/$
          (omipress[ilevel, ip]-omipress[ilevel-1, ip])
        break
      ENDIF
    endfor
  ENDFOR
ENDIF


omilayero3[altnanidx] = !values.f_nan

;;----------------------------
;; OMI vertical column Layer for under specific altitude
;;----------------------------

if height_limit ne 0 then begin
  ;omitozsize = size(omitozs, /dim)
  omiecfsize = size(omiecf, /dim)
  omilayero3 = fltarr(omiecfsize[0])

  altnanidx = where(omialts lt -1.E30, /null)
  omialts[altnanidx] = !values.f_nan

  for ip=0,omiecfsize[0]-1 do BEGIN
    for ilevel=24, 0, -1 do BEGIN
      ilayer = ilevel - 1
      if omialts[ilevel, ip] lt height_limit $
          and omialts[ilevel-1, ip] lt height_limit then BEGIN
        omilayero3[ip] = omilayero3[ip] + omio3s[ilayer, ip]
        
      ENDIF else if omialts[ilevel, ip] lt height_limit $
        and omialts[ilevel-1, ip] ge height_limit then begin 

        omilayero3[ip] = omilayero3[ip] + omio3s[ilayer, ip]*$
          (height_limit - omialts[ilevel, ip])/$
          (omialts[ilevel-1, ip] - omialts[ilevel, ip])
        break
      ENDIF else BEGIN
        continue
      ENDELSE
    endfor
  ENDFOR
ENDIF



wrongidx = where(omilayero3 eq 0, /null)
omilayero3[wrongidx] = !values.f_nan
;if height_limit eq 0 and pressure_limit eq 0 then BEGIN
;ENDIF

  omivars={lon:omilons, lat:omilats, time:omitimes, $
    ;toz:omitozs,$
    SolarZenithAngle:omiszas, $
    ;SolarAzimuthAngle:omisaas, $
    ;ReflectanceCostFunction:omirefcfs, $
    EffectiveCloudFraction:omiecf, $
    EffectiveCloudPressure:omiecp, $
    omilayero3:omilayero3, $
    ;pixlat:pixlats, pixlon:pixlons, pixs:omixtracks, lines:omilines,$
    ;xtracks:omixtracks, $
    ;lines:omilines, $
    ;rowanomaly:omirowanomalys, $
    ;ref331:omir331s, $
    ;cf:omicfs, $
    ;refcf:omirefcfs, $
    ;pcld:omipclds, $
    qf:omigroundpqfs};,bloz:omiblozs}


  ; collocation start

  omilon = omivars.Lon
  omilat = omivars.Lat
  omitime = omivars.Time
  ;omitoz = omivars.toz
  omilayero3 = omivars.omilayero3
  ;omipixlat = omivars.pixlat
  ;omipixlon = omivars.pixlon
  ;omixtracks = omivars.xtracks
  ;omilines = omivars.lines
  omiecf = omivars.EffectiveCloudFraction
  omiecp = omivars.EffectiveCloudPressure
  omiqf = omivars.qf
  ;omirowanomaly = omivars.rowAnomaly


  nanidx = where(omitime lt -5.0E9, /null)
  omitime[nanidx] = !values.f_nan
  caldat, julday(1,1,1993, 0, 0) + omitime/60./60./24., omimons, omidays, $
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
  ds_get_pixcor, gemslon, gemslat, gemslonpixcor,gemslatpixcor
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
  gemslon[where(gemslon lt -180, /null)] = !values.f_nan
  gemslat[where(gemslat lt -90, /null)] = !values.f_nan
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

  for ip = 0, omisize[0]-1 do BEGIN
    x = omilon[ip]
    y = omilat[ip]

    if finite(x) and finite(y) then begin
      result = search_closest_pixel(gemslon, gemslat, x, y)
      ;result = array_indices(gemslon, result)
      closest_idx[ip] = result
    ENDIF
  ENDFOR
  nanidx = where(closest_idx lt -990, /null)
  
  yn_omi_cross_time = intarr(n_elements(omilon))
  yn_omi_cross_time[*] = 0
  for ip = 0, n_elements(omilon)-1 do BEGIN
    if abs(omijulday[ip] - gemsjulday[closest_idx[ip]]) lt 1./24./1. then BEGIN
      yn_omi_cross_time[ip] = 1
    ENDIF
  endfor

  save, filename=savfile, gemsvars, gemslayero3, omivars, omilayero3, closest_idx, yn_omi_cross_time

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
gemssza = gemsvars.SolarZenithAngle
gemssza[where(gemssza lt -1E29 , /null)] = !values.f_nan
nanidx = where(gemsecf lt 0 , /null)
gemsecf[nanidx] = !values.f_nan
gemsidx = []
for ix=0, omisize[0]-1 do BEGIN

  if closest_idx[ix] GE 0 and yn_omi_cross_time[ix] eq 1 then begin
    ;if omivars.ReflectanceCostFunction[ix] lt 30 $
        ;and omivars.SolarZenithAngle[ix] lt 5000 $
    ;if omivars.SolarZenithAngle[ix] lt 30. $
        ;;and omivars.xtracks[ix] ge 3 $
        ;;and omivars.xtracks[ix] le 25 $
        ;;and omivars.EffectiveCloudFractionUV1[ix] le 0.2 $
          ;and gemsecf[closest_idx[ix]] le 0.2 $
          ;and gemsecf[closest_idx[ix]] gt 0 $
          ;and gemssza[closest_idx[ix]] lt 30. $
          ;and abs(gemssza[closest_idx[ix]] - omivars.SolarZenithAngle[ix]) lt 10 then begin

      if pressure_limit ne 0 or height_limit ne 0 then begin
        ;omival = omivars.omilayero3[ix]
        omival = omilayero3[ix]
        gemsval = gemslayero3[closest_idx[ix]]
      endif else begin
        omival = omivars.omilayero3[ix] 
        gemsval = gemstoz[closest_idx[ix]]
      endelse

      ;print,omival, gemsval

      if omival ge 0 and omival le 1000 AND $
          gemsval ge 0 and gemsval le 1000 then begin
        omivals = [omivals, omival]
        gemsvals = [gemsvals, gemsval]
        gemsidx = [gemsidx, closest_idx[ix]]
      endif
    ;endif

  endif
ENDFOR

  mondate = yyyy+mm+dd

  outputpath ='./plot/'
  plot_sat_proj, omilayero3, gemsvars.longitude, gemsvars.latitude, $
    title='GEMS L2 O3P for ' + mondate, $
    range=[20, 60], $
    pngfile=outputpath + 'omi_profoz_collocated_gems_' + mondate + '_under300hpa_tropo_wl310340_me0.5.png', $
    /scp_send

stop
return
end
