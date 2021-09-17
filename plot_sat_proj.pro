pro plot_sat_proj, data, lon, lat, $
  pngfile=pngfile, $
  scp_send=scp_send1, $
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

outdir = file_dirname(pngfile)
if not file_test(outdir) then begin
  file_mkdir, outdir + '/'
endif

out_basename = file_basename(pngfile, '.png')

cgPS_Open, outdir + '/' + out_basename + '.ps'
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

for jy = 0L, dim2-3 do begin
  for ix = 0L, dim1-3 do begin
    if (lat[ix,jy] ge Slat) and (lat[ix,jy] le Nlat) and $
        (lon[ix,jy] ge Llon) and (lon[ix,jy] le Rlon) then begin
      if data[ix,jy] gt 0 and data[ix, jy] le 500 then begin
        xbox = [lon[ix,jy], lon[ix+1,jy], lon[ix+1,jy], lon[ix,jy]]
        ybox = [lat[ix,jy], lat[ix,jy], lat[ix,jy+1], lat[ix,jy+1]]
        xbox = xbox[where(finite(xbox) eq 1, /null)]
        ybox = ybox[where(finite(ybox) eq 1, /null)]
        polyfill, xbox, ybox, $
          color=bytscl(data[ix,jy], $
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
