pro plot_gems_l2_o3p_satproj, filename, scp_dest=scp_dest1, height=height, $
  presunder=presunder, $
  pngfile=pngfile


if not keyword_Set(height) then begin
  height=0
endif
if not keyword_Set(pressureunder) then begin
  pressureunder=0
endif
varlist = ['EffectiveCloudFractionUV', 'ProcessingQualityFlags', $
           'O3' , $
           ;'O3Apriori', 'O3AprioriError', $
           ;'CloudPressure', $
           ;'O3Apriori', 'O3AprioriError',$
           'ColumnAmountO3', $
           ;'SimulatedRadiances', $
           'Latitude' ,'Longitude', $
           'Time','Altitude' ,    $
           ;'Pressure', 'TropopausePressure', $
           ;'Wavelengths', $
           'WavelengthsWholeRange']

data = ds_read_gems_l2_o3p(filename, varlist=varlist)

basename = file_basename(filename, '.nc4') 
lat = data.latitude
lon = data.longitude



datetime_offset = 33
;yyyy = strmid(filename, 16 + datetime_offset, 4, /reverse)
;mm = strmid(filename, 12 + datetime_offset, 2, /reverse)
;dd = strmid(filename, 10 + datetime_offset, 2, /reverse)
;utc = strmid(filename, 7 + datetime_offset, 4, /reverse)

;wloffset = 55
;precoffset = 50
;initwl = strmid(filename, 16 + wloffset, 3, /reverse)
;prec = '0.00' + strmid(filename, 16 + wloffset-6, 1, /reverse)
;print, yyyy, mm, dd, utc, initwl, prec


do_plot_cao3 = [1, 1, 1]
;do_plot_p300 = 0
;do_plot_10km = 1


;;----------------------------
;; call retrieval value
;;----------------------------

outpath = './plot/'
;files = file_search(path+'GEMS_ASIAAOP_AODSSA_'+yyyy+'m'+mm+dd+utc+'.sav',count = files_num)
;restore, files(0)  ;; elon, elat, aodres


;;----------------------------
;; set dimension
;;----------------------------

;dim1 = 1199
;dim2 = 899
sz = size(lon)
dim1 = sz[1]
dim2 = sz[2]

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
cao3 = data.ColumnAmountO3
on_cao3 = 0
if on_cao3 eq 1 then begin
  for icao3 = 1, -1 do begin
    if do_plot_cao3[icao3-1] then begin
      layercao3 = reform(cao3[icao3-1, *, *])
      layercao3[where(layercao3 lt 0, /null)] = !values.f_nan

      if icao3 eq 1 then begin
        columnname = 'Total'
        o3range = [250, 400]
      endif else if icao3 eq 2 then begin
        columnname = 'Stratospheric'
        o3range = [80, 350]
      endif else if icao3 eq 3 then begin
        columnname = 'Tropospheric'
        o3range = [20, 60]
      endif

      ;pic_title =yyyy+'-'+mm+'-'+dd+'_'+utc+'UTC'+ columnname +initwl + '-340nm_' +prec 
      pic_title =yyyy+'-'+mm+'-'+dd+'_'+utc+'UTC'+ columnname 

      ;x_para = elon           
      ;X_name = 'Longitude'
      ;y_para = elat           
      ;Y_name = 'Latitude'

      Slat = -5.
      Nlat = 60.
      Llon = 60.
      Rlon = 160.

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
      dsloadct, 33

      for j = 1L, dim2-2 do begin
        for i = 1L, dim1-2 do begin
          if (lat[i,j] ge Slat) and (lat[i,j] le Nlat) and $
              (lon[i,j] ge Llon) and (lon[i,j] le Rlon) then begin
            if layercao3[i,j] gt 0 and layercao3[i, j] le 500 then begin
              xbox = [lon[i,j], lon[i+1,j], lon[i+1,j], lon[i,j]]
              ybox = [lat[i,j], lat[i,j], lat[i,j+1], lat[i,j+1]]
              if total(finite(xbox) + finite(ybox)) eq 8 then begin
                polyfill, xbox, ybox, $
                  color=bytscl(layercao3[i,j], $
                  min=o3range[0], $
                  max=o3range[1], $
                  top=253)
              endif
            endif
          endif
        endfor
      endfor

      cgMap_Grid, Color='black', GLinestyle=1,londel=londel, latdel=latdel, charsize=1 $
                       ,LONLAB=20 ,LATLAB=160, LABEL=1

      cgMap_continents, color='black'

      cgColorbar, format='(f6.1)', range=o3range, $
        Position=[pos[0]+0.065, pos[1]-0.04, pos[2]-0.065, pos[1]-0.02], $
        Color='black',$
        OOB_High=253, $
        OOB_low=0,$
        divisions=5, $
        bottom=0, $
        ncolor=254, $
        minor=1,  $
        charsize=1.5,$
        XTicklen=1,  $
        XMinor=0, $
        AnnotateColor='black',  $
        Title='GEMS ' + columnname + ' Column O3'


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


      ;pngfile = outpath+pic_title + '.png'
      pngfile = outpath+filename + '.png'
      cgPS2Raster,'gems.ps', pngfile, /png, density = 1000



      if not keyword_Set(scp_dest1) then begin
        scp_dest = 'soodal@164.125.38.179:/home/soodal/works/plot'
      endif else begin
        scp_dest = scp_dest1
      endelse


      ; send image to pc
      spawn, 'scp -P18742 -p ' + pngfile + $
        ' ' + scp_dest

    endif
  endfor
endif


;;----------------------------
;; vertical 10km column Layer
;;----------------------------

if keyword_Set(height) then begin
  ds_gems_l2o3p_accum, data, layero3, height=height
  layero3[where(layero3 le 0)] = !values.f_nan

  columnname = 'SfcTo10KM'
  o3range = [20, 60]
  ;o3range = [20, 200]

  ;pic_title =yyyy+'-'+mm+'-'+dd+'_'+utc+'UTC'+ columnname +initwl + '-340nm_' +prec 
  ;pic_title =yyyy+'-'+mm+'-'+dd+'_'+utc+'UTC'+ columnname
  pic_title = '2020-06-16 0345UTC SfcTo10KM'

  ;x_para = elon           
  ;X_name = 'Longitude'
  ;y_para = elat           
  ;Y_name = 'Latitude'

  Slat = -5.
  Nlat = 60.
  Llon = 60.
  Rlon = 160.

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
  dsloadct, 33

  for j = 1L, dim2-2 do begin
    for i = 1L, dim1-2 do begin
      if (lat[i,j] ge Slat) and (lat[i,j] le Nlat) and $
          (lon[i,j] ge Llon) and (lon[i,j] le Rlon) then begin
        if layero3[i,j] gt 0 and layero3[i, j] le 500 then begin
          xbox = [lon[i,j], lon[i+1,j], lon[i+1,j], lon[i,j]]
          ybox = [lat[i,j], lat[i,j], lat[i,j+1], lat[i,j+1]]
          if total(finite(xbox) + finite(ybox)) eq 8 then begin
            polyfill, xbox, ybox, $
              color=bytscl(layero3[i,j], $
              min=o3range[0], $
              max=o3range[1], $
              top=253)
          endif
        endif
      endif
    endfor
  endfor

  cgMap_Grid, Color='black', GLinestyle=1,londel=londel, latdel=latdel, charsize=1 $
                   ,LONLAB=20 ,LATLAB=160, LABEL=1

  cgMap_continents, color='black'

  cgColorbar, format='(f6.1)', range=o3range, $
    Position=[pos[0]+0.065, pos[1]-0.04, pos[2]-0.065, pos[1]-0.02], $
    Color='black',$
    OOB_High=253, $
    OOB_low=0,$
    divisions=5, $
    bottom=0, $
    ncolor=254, $
    minor=1,  $
    charsize=1.5,$
    XTicklen=1,  $
    XMinor=0, $
    AnnotateColor='black',  $
    Title='GEMS ' + columnname + ' Column O3'


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


  ;pngfile = outpath+pic_title + '.png'
  ;pngfile = outpath+pngfile

  cgPS2Raster,'gems.ps', outpath+pngfile, /png, density = 1000

    if not keyword_Set(scp_dest1) then begin
      scp_dest = 'soodal@164.125.38.179:/home/soodal/works/plot'
    endif else begin
      scp_dest = scp_dest1
    endelse

  ; send image to pc
  spawn, 'scp -P18742 -p ' + outpath+pngfile + $
    ' ' + scp_dest

endif




;;----------------------------
;; vertical 300 hPa column Layer
;;----------------------------

if keyword_Set(presunder) then begin
  ds_gems_l2o3p_accum, data, layero3, hpa=300
  layero3[where(layero3 le 0)] = !values.f_nan

  columnname = 'SfcTo300hPa'
  o3range = [20, 60]
  ;o3range = [20, 200]

  ;pic_title =yyyy+'-'+mm+'-'+dd+'_'+utc+'UTC'+ columnname +initwl + '-340nm_' +prec 
  ;pic_title =yyyy+'-'+mm+'-'+dd+'_'+utc+'UTC'+ columnname

  ;x_para = elon           
  ;X_name = 'Longitude'
  ;y_para = elat           
  ;Y_name = 'Latitude'

  Slat = -5.
  Nlat = 60.
  Llon = 60.
  Rlon = 160.

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
  dsloadct, 33

  for j = 1L, dim2-2 do begin
    for i = 1L, dim1-2 do begin
      if (lat[i,j] ge Slat) and (lat[i,j] le Nlat) and $
          (lon[i,j] ge Llon) and (lon[i,j] le Rlon) then begin
        if layero3[i,j] gt 0 and layero3[i, j] le 500 then begin
          xbox = [lon[i,j], lon[i+1,j], lon[i+1,j], lon[i,j]]
          ybox = [lat[i,j], lat[i,j], lat[i,j+1], lat[i,j+1]]
          if total(finite(xbox) + finite(ybox)) eq 8 then begin
            polyfill, xbox, ybox, $
              color=bytscl(layero3[i,j], $
              min=o3range[0], $
              max=o3range[1], $
              top=253)
          endif
        endif
      endif
    endfor
  endfor

  cgMap_Grid, Color='black', GLinestyle=1,londel=londel, latdel=latdel, charsize=1 $
                   ,LONLAB=20 ,LATLAB=160, LABEL=1

  cgMap_continents, color='black'

  cgColorbar, format='(f6.1)', range=o3range, $
    Position=[pos[0]+0.065, pos[1]-0.04, pos[2]-0.065, pos[1]-0.02], $
    Color='black',$
    OOB_High=253, $
    OOB_low=0,$
    divisions=5, $
    bottom=0, $
    ncolor=254, $
    minor=1,  $
    charsize=1.5,$
    XTicklen=1,  $
    XMinor=0, $
    AnnotateColor='black',  $
    Title='GEMS ' + columnname + ' Column O3'


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
  pngfile = outpath+pngfile
  cgPS2Raster,'gems.ps', pngfile, /png, density = 1000

    if not keyword_Set(scp_dest1) then begin
      scp_dest = 'soodal@164.125.38.179:/home/soodal/works/plot'
    endif else begin
      scp_dest = scp_dest1
    endelse

  ; send image to pc
  spawn, 'scp -P18742 -p ' + pngfile + $
    ' ' + scp_dest

endif
end
