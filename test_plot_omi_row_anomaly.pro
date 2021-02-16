limit=[-10, 80, 60, 160]

year = 2020
month = 8

day = 11
hour = 3
minute = 45

yyyy = string(year, format='(i04)')
mm = string(month, format='(i02)')
dd = string(day, format='(i02)')
hh = string(hour, format='(i02)')
mi = string(minute, format='(i02)')

savpath = './collocate_gemsl2o3p_omo3pr/'
savfile = savpath + 'col_gems_omi_'+yyyy+mm+dd+hh+mi+'.sav'


;if not file_test(savfile) then BEGIN


  ;read_omo3prof

  ; initialize

  omilons = []  
  omilats = []  
  omitozs = []  
  omitimes = []  
  omipixs = []  
  omilines = []
  pixlons = []  
  pixlats = []  
  omir331s = [] 
  omipclds = [] 
  omicfs = []  
  omiblozs = []
  omigroudpqfs = []
  omirowanomalys = []
  omiszas = []

  ; path for GEMS
  gemsfn = '/data2/L2_GEMS/val_1008/GK2_GEMS_O3P_20200806_0345.nc4'

  ; read GEMS O3P data
  gemsvars = ds_read_gems_l2_o3p(gemsfn)

  ; path for OMI

  omipath = '/data2/OMI/gdata/' + yyyy + '/L2-OMO3PR/' 
  omipixpath = '/data2/OMI/gdata/' + yyyy + '/L2-OMPIXCOR/' 

  omifiles=FILE_SEARCH(omipath+'OMI-Aura*'+yyyy+'m'+mm+dd+'*.he5', count=nfiles)
  omifile = '/data2/OMI/gdata/2020/L2-OMO3PR/OMI-Aura_L2-OMO3PR_2020m0806t0345-o85422_v003-2020m0807t120244.SUB.he5'



  ; read OMI data
  ;for ifile=0, nfiles-1 do begin
    fi = file_info(omifile)
    if fi.size gt 1500 then begin
      ds_READ_omi_L2_omo3pr, omifile, omiresult

      tmp=WHERE(omiresult.longitude GE limit[1] AND omiresult.longitude LE limit[3] AND $
                omiresult.latitude GE limit[0] AND omiresult.latitude LE limit[2], ntmp)

      print, '  READ omifile : ', omifile, ntmp
      IF ntmp NE 0 THEN BEGIN 
        header='OMI-Aura_L2-OMPIXCOR_'
        time=STRMID(omifile,46+4,14,/rev)
        print,time

        pixfile=FILE_SEARCH(omipixpath + header+time+'*.he5', count=npixfile)
        IF npixfile NE 1 THEN BEGIN
          PRINT, '  No Pixfile'
          STOP
        ENDIF ; npixfile

        pixfi = file_info(pixfile)
        if pixfi.size gt 1500 then begin

          ds_READ_omi_L2_PIXCOR, pixfile, pixcor
          pixdim=SIZE(pixcor.corlon1,/dim)
          nx=pixdim[0]
          ny=pixdim[1]
          pixlon=REFORM(pixcor.corlon1,[nx*ny,4])
          pixlat=REFORM(pixcor.corlat1,[nx*ny,4])     
          ;pixlon=pixlon[tmp,*]
          ;pixlat=pixlat[tmp,*]

          pixlats  =  [pixlats, pixlat]
          pixlons  =  [pixlons, pixlon]


        ;        omitime=pixcor.jtime[tmp]
          omilons  =  [omilons, omiresult.longitude[tmp] ]
          omilats  =  [omilats, omiresult.latitude[tmp] ]
          omitozs  =  [omitozs, omiresult.toz[tmp] ]
          omitimes =  [omitimes,pixcor.jtime1[tmp]]
          omiszas = [omiszas, omiresult.SolarZenithAngle[tmp]]

          omipixs  =  [omipixs,omiresult.xtracks[tmp]]
          omilines =  [omilines,omiresult.lines[tmp]]

          ;omicfs =    [omicfs,   omiresult.cf[tmp]   ]
          ;omiblozs =  [omiblozs,   omiresult.blo3[tmp]   ]
          omipclds =  [omipclds, omiresult.pc[tmp]     ]
          ;omir331s =  [omir331s, omiresult.ref331[tmp] ]
          omigroudpqfs   =  [omigroudpqfs,   omiresult.groundpqf[tmp]     ]
          omirowanomalys = [omirowanomalys, omiresult.rowAnomaly[tmp]]
        endif
      ENDIF
    ENDIF  ; ntmp
  ;ENDFOR  ; ifile




  ; collocation

  omilon = omivars.longitude
  omilat = omivars.latitude
  omitime = omivars.time
  omitoz = omivars.toz
  omipixlat = omivars.pixlat
  omipixlon = omivars.pixlon
  omipixs = omivars.pixs
  omilines = omivars.lines
  omipcld = omivars.pcld
  omiqf = omivars.qf
  omirowanomaly = omivars.rowAnomaly


  caldat, julday(1,1,1993, 0, 0) + omitime/60./60./24., omimon, omiday, $
    omiyear, omihour, omiminute, omisec

  omijulday = julday(omimon, omiday, omiyear, omihour, omiminute, omisec) 
  


  gemslon = gemsvars.Longitude
  gemslat = gemsvars.Latitude
  gemstime = gemsvars.Time
  gemslon[where(gemslon lt -180, /null)] = !values.f_nan
  gemslat[where(gemslat lt -90, /null)] = !values.f_nan
  ;gemstime[where(gemstime lt 0, /null)] = !values.f_nan
  gemstime = julday(month, day, year, hour, minute)
  gemssize = size(gemslon, /dim)
  
  yn_omi_cross_time = intarr(n_elements(omilon))
  yn_omi_cross_time[where(abs(omijulday-gemstime) lt 1./24., /null)] = 1

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
      result = search_closest_pixel(omilon, omilat, x, y)
      ;result = array_indices(gemslon, result)
      closest_idx[ix] = result
    ENDIF
  ENDFOR
  save, filename=savfile, gemsvars, omivars, closest_idx, yn_omi_cross_time
  

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

;endif else BEGIN
  ;RESTORE, savfile
;ENDELSE
  gemstoz = reform(gemsvars.ColumnAmountO3[0, *, *])
  gemsvals = []
  omivals = []
  omisize = size(omivars.longitude, /dim)
  for ix=0, omisize[0]-1 do BEGIN

    if closest_idx[ix] GE 0 and yn_omi_cross_time[ix] gt 0 then begin
      omival = omitoz[ix]
      gemsval = gemstoz[closest_idx[ix]]

      if omival ge 0 and omival le 1000 AND $
          gemsval ge 0 and gemsval le 1000 then begin
        omivals = [omivals, omival]
        gemsvals = [gemsvals, gemsval]
      endif

    endif
  ENDFOR
  plot_gems_validation, gemsvals, omivals, filename='./plot/gems_l2_o3p_val_with_omi_20200806T0345.png'

end
