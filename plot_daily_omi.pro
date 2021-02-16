pro plot_daily_omi, year, month, day, hour, minute, $
  omivals, $
  hpa=pressure_limit, height=height_limit, $
  range=range


if not keyword_Set(pressure_limit) then begin
  pressure_limit = 0
endif

if not keyword_Set(height_limit) then begin
  height_limit = 0
endif


limit=[-10, $ ;minimum latitude
  80, $ ;minimum longitude
  10, $ ;maximum latitude
  160] ;maximum longitude

yyyy = string(year, format='(i04)')
mm = string(month, format='(i02)')
dd = string(day, format='(i02)')
hh = string(hour, format='(i02)')
mi = string(minute, format='(i02)')

pic_title = 'OMI_daily_'+yyyy+mm+dd
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

LoadCT, 33, Ncolors=254,bottom=1

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
    ds_READ_omi_L2_profoz, omifiles[ifile], profoz

    
    sz = size(profoz.latitude)
    dim1 = sz[1]
    dim2 = sz[2]
    time=STRMID(omifiles[ifile],46,14,/rev)

    ;endif
    ;;----------------------------
    ;; OMI vertical column Layer for under specific pressure
    ;;----------------------------
    if pressure_limit ne 0 then begin
      ;omitozsize = size(omitozs, /dim)
      omiecfsize = size(profoz.EffectiveCloudFraction, /dim)
      omilayero3 = fltarr(omiecfsize[0], omiecfsize[1])
      omialts = profoz.ProfileLevelAltitude
      omipress = profoz.ProfileLevelPressure
      

      altnanidx = where(omialts lt -1.E30, /null)
      omialts[altnanidx] = !values.f_nan
      o3rp = profoz.O3RetrievedProfile
      nanidx = where(o3rp lt -1.0e30, /null)
      o3rp[nanidx] = !values.f_nan
      o3rp[*,20:23, *] = !values.f_nan

      for ix=0,omiecfsize[0]-1 do BEGIN
      for iy=0,omiecfsize[1]-1 do BEGIN
        for ilevel=24, 0, -1 do BEGIN
          ilayer = ilevel - 1
          if omipress[ilevel, ix, iy] gt pressure_limit $
              and omipress[ilevel-1, ix, iy] gt pressure_limit then BEGIN
            omilayero3[ix, iy] = omilayero3[ix, iy] + o3rp[ilayer, ix, iy]
            
          ENDIF else if omipress[ilevel, ix, iy] gt pressure_limit $
            and omipress[ilevel-1, ix, iy] le pressure_limit then begin 
            ; upper bound is higher than 300hPa 

            omilayero3[ix, iy] = omilayero3[ix, iy] + o3rp[ilayer, ix, iy]*$
              (omipress[ilevel, ix, iy]-pressure_limit)/$
              (omipress[ilevel, ix, iy]-omipress[ilevel-1, ix, iy])
            break
          ENDIF
        endfor
      ENDFOR
      ENDFOR
    ENDIF
    nanidx = where(omilayero3 le 0 or omilayero3 gt 100, /null)
    omilayero3[nanidx] = !values.f_nan

    ;omilayero3[altnanidx] = !values.f_nan

    ;;----------------------------
    ;; OMI vertical column Layer for under specific altitude
    ;;----------------------------

    if height_limit ne 0 then begin
      ;omitozsize = size(omitozs, /dim)
      omiecfsize = size(profoz.EffectiveCloudFraction, /dim)
      omilayero3 = fltarr(omiecfsize[0], omiecfsize[1])
      omialts = profoz.ProfileLevelAltitude
      omipress = profoz.ProfileLevelPressure

      altnanidx = where(omialts lt -1.E30, /null)
      omialts[altnanidx] = !values.f_nan

      for ip=0,long(omiecfsize[0])*omiecfsize[1]-1 do BEGIN
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

    for j = 0L, dim2-1 do begin
      for i = 0L, dim1-1 do begin
          if omilayero3[i,j] gt 0 and omilayero3[i, j] le 500 then begin
            xbox = [profoz.longitudepixelcorner[i,j], $
              profoz.longitudepixelcorner[i+1,j], $
              profoz.longitudepixelcorner[i+1,j], $
              profoz.longitudepixelcorner[i,j]]
            ybox = [profoz.latitudepixelcorner[i,j], $
              profoz.latitudepixelcorner[i,j], $
              profoz.latitudepixelcorner[i,j+1], $
              profoz.latitudepixelcorner[i,j+1]]
            polyfill, xbox, ybox, $
              color=bytscl(omilayero3[i,j], $
                min=range[0], $
                max=range[1], $
                top=253)
          endif
        ;endif
      endfor
    endfor

  ENDIF
  idx = where(omilayero3 gt 0 and omilayero3 lt 100, ntmp, /null)
  print, ntmp
ENDFOR  ; ifile

cgMap_Grid, Color='black', GLinestyle=1,londel=londel, latdel=latdel, charsize=1 $
                 ,LONLAB=20 ,LATLAB=160, LABEL=1

cgMap_continents, color='black'

if not keyword_Set(range) then begin
  range = [min(omilayero3, /nan), max(omilayero3, /nan)]
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

cgPS_Close


outpath = './plot/'
pngfile = outpath+pic_title + '_10km.png'

cgPS2Raster,'gems.ps', pngfile, /png, density = 1000

scp_dest = 'soodal@164.125.38.179:/home/soodal/works/plot'

; send image to pc
spawn, 'scp -P18742 -p ' + pngfile + $
  ' ' + scp_dest
return
end
