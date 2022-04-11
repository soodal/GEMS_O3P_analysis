pro plot_sat_proj, data, lon, lat, $
  pngfile=pngfile, $
  scp_send=scp_send1, $
  scp_dest=scp_dest1, $
  range=range, $
  title=pic_title, $
  colortable=ct, $
  cornerpixel=cornerpixel, $
  ctreverse=ctreverse, $
  slat=slat, $
  nlat=nlat, $
  llon=llon, $
  rlon=rlon, $
  cb_title=cb_title
if not keyword_Set(pic_title) then begin
  pic_title = 'test'
endif

if not keyword_Set(cb_title) then begin
  cb_title = pic_title
endif

if not keyword_Set(range) then begin
  range = minmax(data)
endif

if not keyword_Set(ct) then begin
  ct = 33
endif

if not keyword_Set(ctreverse) then begin
  ctreverse = 0
endif

if not keyword_Set(slat) then begin
  Slat = -5.
endif
if not keyword_Set(nlat) then begin
  Nlat = 60.
endif
if not keyword_Set(llon) then begin
  Llon = 60.
endif
if not keyword_Set(rlon) then begin
  Rlon = 170.
endif

if not keyword_Set(cornerpixel) then begin
  cornerpixel = 0
endif

outdir = file_dirname(pngfile)
if not file_test(outdir) then begin
  file_mkdir, outdir + '/'
endif

out_basename = file_basename(pngfile, '.png')

cgPS_Open, outdir + '/' + out_basename + '.ps'
cgDisplay, aspect=8./9.

pos = [0.1,0.15,0.9,0.90]
charsize = (!D.Name EQ 'PS') ? cgDefCharsize()*0.7 : cgDefCharsize()*0.75
thick = (!D.Name EQ 'PS') ? 6 :3
londel = 10
latdel = 10

cgMap_set,20,120, Limit=[-60, 0, 60, 360],$
              pos=pos,charsize=charsize,$
              /satellite,/iso,$
              /horizon, $
              title=pic_title

;LoadCT, 22, Ncolors=254,bottom=1
if ctreverse then  begin
  dsloadct, ct, /ctreverse
endif ELSE begin
  dsloadct, ct
ENDELSE

if cornerpixel then begin
  sz = size(data)
  dim1 = sz[1]
  for ip = 0L, dim1-1 do begin
    print, ip
    print, lat[ip, *]
    print, lon[ip, *]
    if (min(lat[ip, *]) ge Slat) and (max(lat[ip, *]) le Nlat) and $
        (min(lon[ip, *]) ge Llon) and (max(lon[ip, *]) le Rlon) then begin
      print, data[ip]
      if finite(data[ip]) eq 1 or data[ip] gt -999 then begin
        xbox = lon[ip, *]
        ybox = lat[ip, *]
        xbox = xbox[where(finite(xbox) eq 1, /null)]
        ybox = ybox[where(finite(ybox) eq 1, /null)]
        print, xbox, ybox
        print, data[ip]
        polyfill, xbox, ybox, $
          color=bytscl(data[ip], $
            min=range[0], $
            max=range[1], $
            top=252)
      endif
    endif
  endfor
endif else begin
  sz = size(data)
  dim1 = sz[1]
  dim2 = sz[2]
  for jy = 0L, dim2-3 do begin
    for ix = 0L, dim1-3 do begin
      if (lat[ix,jy] ge Slat) and (lat[ix,jy] le Nlat) and $
          (lon[ix,jy] ge Llon) and (lon[ix,jy] le Rlon) then begin
        if finite(data[ix,jy]) gt 0 or data[ix, jy] gt -999 then begin
          xbox = [lon[ix,jy], lon[ix+1,jy], lon[ix+1,jy], lon[ix,jy]]
          ybox = [lat[ix,jy], lat[ix,jy], lat[ix,jy+1], lat[ix,jy+1]]
          xbox = xbox[where(finite(xbox) eq 1, /null)]
          ybox = ybox[where(finite(ybox) eq 1, /null)]
          polyfill, xbox, ybox, $
            color=bytscl(data[ix,jy], $
              min=range[0], $
              max=range[1], $
              top=252)
        endif
      endif
    endfor
  endfor
ENDELSE

cgMap_Grid, Color='black', GLinestyle=1,londel=londel, latdel=latdel, charsize=1 $
                 ,LONLAB=20 ,LATLAB=160, LABEL=1

cgMap_continents, color='black'

cgColorbar, format='(f6.1)', range=range, $
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
  Title=cb_title

cgPS_Close

if not keyword_Set(pngfile) then begin
  pngfile = 'test.png'
endif

cgPS2Raster, outdir + '/' + out_basename + '.ps', pngfile, /png, density = 1000

file_delete, outdir + '/' + out_basename + '.ps'

if not keyword_Set(scp_dest1) then begin
  scp_dest1 = 'soodal@164.125.38.179:/home/soodal/works/plot'
endif 

; send image to pc
if keyword_set(scp_send1) then begin
  spawn, 'scp -P18742 -p ' + pngfile + $
  ' ' + scp_dest1
endif
end
