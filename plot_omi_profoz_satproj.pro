pro plot_omi_satproj, corlon, corlat, data, $
  filename=filename, $
  pngfile=pngfile, $
  scp_send=scp, $
  scp_dest=scp_dest1, $
  range=range, $
  title=pic_title, $
  cb_title=cb_title


;;----------------------------
;; call retrieval value
;;----------------------------

outpath = '/data1/gems/o3p/ds/plot/'
;files = file_search(path+'GEMS_ASIAAOP_AODSSA_'+yyyy+'m'+mm+dd+utc+'.sav',count = files_num)
;restore, files(0)  ;; elon, elat, aodres


;;----------------------------
;; set dimension
;;----------------------------

;dim1 = 1199
;dim2 = 899
sz = size(lon, /dim)
dim1 = sz[0]
;dim2 = sz[1]

;lon = replicate(!values.f_nan, dim1, dim2)
;lat = replicate(!values.f_nan, dim1, dim2)
;pic_aod = replicate(!values.f_nan, dim1, dim2)

;lon = reform(elon, dim1, dim2)
;lat = reform(elat, dim1, dim2)
;pic_aod = reform(aodres, dim1, dim2)
;data = 





;;----------------------------
;; Column Amount O3 plot procedure
;;----------------------------
;cao3 = data.ColumnAmountO3

;for icao3 = 1, 2 do begin
  ;layercao3 = reform(cao3[icao3-1, *, *])
  ;layercao3[where(layercao3 lt 0, /null)] = !values.f_nan

  ;if icao3 eq 1 then begin
    ;columnname = 'Total'
    ;o3range = [220, 320]
    
  ;endif else if icao3 eq 2 then begin
    ;columnname = 'Row Anomaly'
    ;o3range = [0, 1]
  ;endif else if icao3 eq 3 then begin
    ;columnname = 'Tropospheric'
    ;o3range = [20, 60]
  ;endif

  ;pic_title =yyyy+'-'+mm+'-'+dd+'_'+utc+'UTC'+ columnname
  pic_title = 'OMI_TOZ_20200806T0345'

  ;x_para = elon           
  ;X_name = 'Longitude'
  ;y_para = elat           
  ;Y_name = 'Latitude'

  Slat = -5.
  Nlat = 60.
  Llon = 75.
  Rlon = 150.

  cgPS_Open, 'gems.ps'
  cgDisplay, aspect=0.7

  pos = [0.1,0.15,0.9,0.90]
  charsize = (!D.Name EQ 'PS') ? cgDefCharsize()*0.7 : cgDefCharsize()*0.75
  thick = (!D.Name EQ 'PS') ? 6 :3
  londel = 10
  latdel = 10

  cgMap_set,26,128, Limit=[-60, 0, 60, 360],$
                pos=pos,charsize=charsize,$
                /satellite,/iso,$
                /horizon, $
                title=pic_title

  ;LoadCT, 22, Ncolors=254,bottom=1
  LoadCT, 33, Ncolors=254,bottom=1

  ;for j = 1L, dim2-2 do begin
    ;for i = 1L, dim1-2 do begin
      ;if (lat[i,j] ge Slat) and (lat[i,j] le Nlat) and $
          ;(lon[i,j] ge Llon) and (lon[i,j] le Rlon) then begin
        ;;if data[i,j] gt 0 and data[i, j] le 500 then begin
          ;if keyword_Set(corlon) and keyword_set(corlat) then begin
            ;xbox = corlon[i, j, *]
            ;ybox = corlat[i, j, *]
          ;endif else begin
            ;xbox = [lon[i,j], lon[i+1,j], lon[i+1,j], lon[i,j]]
            ;ybox = [lat[i,j], lat[i,j], lat[i,j+1], lat[i,j+1]]
          ;ENDELSE
          ;polyfill, xbox, ybox, $
            ;color=bytscl(data[i,j], $
              ;min=min(data), $
              ;max=max(data))
        ;;endif
      ;endif
    ;endfor
  ;endfor

    for i = 1L, dim1-2 do begin
      if (lat[i] ge Slat) and (lat[i] le Nlat) and $
          (lon[i] ge Llon) and (lon[i] le Rlon) then begin
          if keyword_Set(corlon) and keyword_set(corlat) then begin
            xbox = corlon[i, j, *]
            ybox = corlat[i, j, *]
          endif else begin
            xbox = [lon[i]-.1, lon[i]+.1, lon[i]+.1, lon[i]-.1]
            ybox = [lat[i]-.1, lat[i]-.1, lat[i]+.1, lat[i]+.1]
            ;xbox = [lon[i,j], lon[i+1,j], lon[i+1,j], lon[i,j]]
            ;ybox = [lat[i,j], lat[i,j], lat[i,j+1], lat[i,j+1]]
          ENDELSE
          polyfill, xbox, ybox, $
            color=bytscl(data[i], $
              min=range[0], $
              max=range[1])
        ;endif
      endif
    endfor

  cgMap_Grid, Color='black', GLinestyle=1,londel=londel, latdel=latdel, charsize=1 $
                   ,LONLAB=20 ,LATLAB=160, LABEL=1

  cgMap_continents, color='black'

  if not keyword_Set(range) then begin
    range = [min(data), max(data)]
  endif
  cgColorbar, format='(f6.1)', range=range, $
    Position=[pos[0]+0.065, pos[1]-0.04, pos[2]-0.065, pos[1]-0.02], $
    color=0,$
    OOB_High=254, $
    OOB_low=1,$
    divisions=5, $
    bottom=1, $
    ncolor=256, $
    minor=1,  $
    charsize=1.5,$
    XTicklen=1,  $
    XMinor=0, $
    AnnotateColor='black',  $
    Title=title


  ;-----------------------------------------------------------------------
  ; * colorbar title for each product
  ;
  ; * NO2 --> title =  'GEMS NO!d2 !n [x 10!e15!n molecule/cm !e2 !n]',
  ; * TropNO2 --> title =  'GEMS Trop NO!d2 !n [x 10!e15!n molecule/cm !e2 !n]',
  ; * SO2 --> title =  'GEMS SO!d2 !n [DU]',
  ; * HCHO -->  title =  'GEMS HCHO [x 10!e16!n molecule/cm !e2 !n]',
  ; * CHOCHO --> title =  'GEMS CHOCHO [x 10!e16!n molecule/cm !e2 !n]',
  ; * O3T--> title= 'GEMS O3T !n [DU]',
  ; * Trop O3T--> title= 'GEMS TropO3 !n [DU]',
  ; * AOD,SSA, ECF, CRF -->  title= 'GEMS AOD 443nm',title= 'GEMS ECF'
  ; * ALH, AEH -->  title= 'GEMS AEH [km]',
  ; * Surface Reflectance -->  title= 'GEMS SR',
  ; * ECP --> title =  'GEMS ECP [mb]'
  ; * UVI --> title =  'GEMS UVI'
  ;----------------------------------------------------------------

  cgPS_Close


  pngfile = outpath+pic_title + '.png'

  cgPS2Raster,'gems.ps', filename, /png, density = 1000

  scp_dest = 'soodal@164.125.38.179:/home/soodal/WindowsHome/Pictures/scp'

  ; send image to pc
  spawn, 'scp -P18742 -p ' + filename + $
    ' ' + scp_dest

;endfor



end
