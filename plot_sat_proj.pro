pro plot_sat_proj, data, lon, lat, $
  pngfile=pngfile, $
  scp_send=scp, $
  scp_dest=scp_dest1, $
  range=range, $
  title=pic_title, $
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
Slat = -5.
Nlat = 60.
Llon = 60.
Rlon = 170.

cgPS_Open, 'gems.ps'
cgDisplay, aspect=0.7

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
dsloadct, 33

sz = size(data)
dim1 = sz[1]
dim2 = sz[2]

for j = 0L, dim2-3 do begin
  for i = 0L, dim1-3 do begin
    if (lat[i,j] ge Slat) and (lat[i,j] le Nlat) and $
        (lon[i,j] ge Llon) and (lon[i,j] le Rlon) then begin
      if data[i,j] gt 0 and data[i, j] le 500 then begin
        xbox = [lon[i,j], lon[i+1,j], lon[i+1,j], lon[i,j]]
        ybox = [lat[i,j], lat[i,j], lat[i,j+1], lat[i,j+1]]
        polyfill, xbox, ybox, $
          color=bytscl(data[i,j], $
            min=range[0], $
            max=range[1], $
            top=253)
      endif
    endif
  endfor
endfor

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

cgPS2Raster,'gems.ps', pngfile, /png, density = 1000

if not keyword_Set(scp_dest1) then begin
  scp_dest = 'soodal@164.125.38.179:/home/soodal/works/plot'
endif else begin
  scp_dest = scp_dest1
endelse

; send image to pc
if keyword_set(scp) then begin
  spawn, 'scp -P18742 -p ' + pngfile + $
  ' ' + scp_dest
endif

end
