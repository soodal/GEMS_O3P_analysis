pro  plot_gems_field_of_regards, pngfile=pngfile, scp_send=scp, $
  scp_dest=scp_dest1, title=pic_title

nominal = "/data1/L1C_GEMS/float/1x1/nopolc_eosrl/GK2_GEMS_L1C_20200616_0245_NOR_0694.nc"
full_central = "/data1/L1C_GEMS/float/1x1/nopolc_eosrl/GK2_GEMS_L1C_20200616_0545_NOR_0694.nc"
half_korea = "/data2/L1C_GEMS/L1C/GK2_GEMS_L1C_20201002_0045_NOR_693_v2.nc"
full_west = "/data2/L1C_GEMS/L1C/GK2_GEMS_L1C_20201002_0445_NOR_694.nc"
full_west1 = "/data2/L1C_GEMS/L1C/GK2_GEMS_L1C_20201002_0645_NOR_624_v2.nc"
full_west2 = "/data2/L1C_GEMS/L1C/GK2_GEMS_L1C_20201002_0545_NOR_440_v2.nc"
;half_east = "/data2/L1C_GEMS/L1C/GK2_GEMS_L1C_20201002_0"

fns = [nominal, $
  full_central, $
  half_korea, $
  full_west, $
  full_west1, $
  full_west2];, half_east]

colortxt = ['orange', $
  'sky blue', $
  'teal', $
  'gold', $
  'blu6', $
  'org6', $
  'hot pink']

boxcolor = [cgcolor('red'), $
  cgcolor('blue'), $
  cgcolor('magenta'),$
  cgcolor('yellow'), $
  cgcolor('gold'), $
  cgcolor('goldenrod'), $
  cgcolor('black')]

titles = ['Nominal', $
  'Full central', $
  'Half korea', $
  'Full west', $
  'Full west1', $
  'Full west2'];, 'Half east']

if not keyword_Set(pic_title) then begin
  pic_title = 'GEMS field of regards'
endif

device, decomposed=1
cgPS_Open, 'gems.ps';, font=1, TT_Font='Times'
;cgDisplay, aspect=0.7

pos = [0.1,0.10,0.9,0.90]
charsize = (!D.Name EQ 'PS') ? cgDefCharsize()*0.7 : cgDefCharsize()*0.75
thick = (!D.Name EQ 'PS') ? 6 :3
londel = 20
latdel = 20

;dsloadct, 33
latmin = -40
lonmin = 60
latmax = 90
lonmax = 180

bluemarblefn = './data/world.topo.bathy.200407.3x5400x2700.png'

xstart = -180
ystart = -90
xsize = 5400
ysize = 2700

; read image file
bm = image(bluemarblefn, $
  grid_units='degrees', $
  image_dimensions=[5400, 2700], $
  image_location=[-180, -90], $
  xrange=[-180, 180], $
  yrange=[-90, 90], $
  margin=0, $
  /buffer)

; get data from image
bm.getdata, bm_data, x, y

; set the map projection to plot
center_lon = 126.0
center_lat = 0.0
cgMap_set,center_lat,center_lon, $
  ;Limit=[latmin, lonmin, latmax, lonmax],$
  ;Limit=[-90, -180, 90, 180],$
  pos=pos, $
  charsize=charsize,$
  /satellite, $
  sat_p=[5.6, 17.0, 0], $
  /iso,$
  /horizon, $
  /noborder, $
  title=pic_title

warped_image1 = map_image(reform(bm_data[0, *, *]), $
  xstart, ystart, $
  xsize, ysize, $
  latmin=-90, $
  latmax=90, $
  lonmin=-180, $
  lonmax=180, $
  scale=1.0, $ ; this keyword for high resolution image
  mask=mask, $
  compress=1)
warped_image2 = map_image(reform(bm_data[1, *, *]), $
  xstart, ystart, $
  xsize, ysize, $
  latmin=-90, $
  latmax=90, $
  lonmin=-180, $
  lonmax=180, $
  scale=1.0, $ ; this keyword for high resolution image
  mask=mask, $
  compress=1)
warped_image3 = map_image(reform(bm_data[2, *, *]), $
  xstart, ystart, $
  xsize, ysize, $
  latmin=-90, $
  latmax=90, $
  lonmin=-180, $
  lonmax=180, $
  scale=1.0, $ ; this keyword for high resolution image
  mask=mask, $
  compress=1)


nanidx = where(mask eq 0, /null)
warped_image1[nanidx] = 255
warped_image2[nanidx] = 255
warped_image3[nanidx] = 255

sz = size(warped_image1, /dimension)
a = bytarr(3, sz[0], sz[1])
a[0, *, *] = warped_image1
a[1, *, *] = warped_image2
a[2, *, *] = warped_image3

tv, a, xstart, ystart, xsize=xsize, ysize=ysize, true=1

for i =0, n_elements(fns)-1 do begin
  lonlat = ds_read_gems_l1c_lonlat(fns[i])
  lon = lonlat.pixel_longitude
  lat = lonlat.pixel_latitude

  sz = size(lon)
  dim1 = sz[1]
  dim2 = sz[2]

  xbox = [reform(lon[0:dim1-1, 0]), reform(lon[dim1-1, 0:dim2-1]), $
    reverse(reform(lon[0:dim1-1,dim2-1])), reverse(reform(lon[0, 0:dim2-1]))]
  ybox = [reform(lat[0:dim1-1, 0]), reform(lat[dim1-1, 0:dim2-1]), $
    reverse(reform(lat[0:dim1-1,dim2-1])), reverse(reform(lat[0, 0:dim2-1]))]

  idx = where(finite(xbox) eq 1 and finite(ybox) eq 1, /null)

  xbox = xbox[idx]
  ybox = ybox[idx]

  case i of
    0: oplot, xbox, ybox, color=cgcolor(colortxt[i]), thick=5
    1: oplot, xbox, ybox, color=cgcolor(colortxt[i]), thick=5
    2: oplot, xbox, ybox, color=cgcolor(colortxt[i]), thick=5
    3: oplot, xbox, ybox, color=cgcolor(colortxt[i]), thick=5
    4: oplot, xbox, ybox, color=cgcolor(colortxt[i]), thick=5
    5: oplot, xbox, ybox, color=cgcolor(colortxt[i]), thick=5
  endcase
endfor

plots, 128.2, 0, psym=2, symsize=1.2, color=cgcolor(colortxt[6]), /data, thick=6

xyouts, 118.2, 2, 'GEMS', charsize=1.2, color=cgcolor(colortxt[6])

cgMap_Grid, Color='white', Linestyle=1,londel=londel, latdel=latdel, charsize=1 $
   ,LONLAB=-10 ,LATLAB=160 $
   ,LATNAMES=['', '', '45S', '20S', '5S', '0', '20N', '45N', '', ''] $
   ,LATS=[-80, -60, -45, -20, -5, 0, 20, 45, 60, 80] $
   ,LONNAMES=['', '80E', '90E', '100E', '120E', '140E', '160E', '180E', ''] $
   ,LONS=[60, 80, 90, 100, 120, 140, 160, 180, 200] $
   ,LABEL=1

map_continents, color=cgcolor('white'), linestyle=0, /hires, thick=2
;map_continents, color=cgcolor('white'), linestyle=0, /hires, thick=2
;map_continents, color=cgcolor('black'), /horizon, /hires, thick=1, 

; again for horizon
cgMap_Grid, Color='black', Linestyle=6, /horizon $
   ,LONLAB=-90 ,LATLAB=center_lon-90 $
   ,LATNAMES=['', ''] $
   ,LATS=[-90, 90] $
   ,LONNAMES=['', ''] $
   ,LONS=[center_lon-90, center_lon+90] $
   ,LABEL=1

;boxcolor = [cgcolor('red'), cgcolor('blue'), cgcolor('magenta'),$
  ;cgcolor('yellow'), cgcolor('gold'), cgcolor('goldenrod')] ;, cgcolor('black')]
cglegend, /box, title=titles, $
  color=[colortxt[0], $
    colortxt[1], $
    colortxt[2], $
    colortxt[3], $
    colortxt[4], $
    colortxt[5]], $
  /background, $
  bg_color='almond', $
  psym=[-0, $
    -0, $
    -0, $
    -0, $
    -0, $
    -0], $
  linestyles=[0, $
    0, $
    0, $
    0, $
    0, $
    0], $
  thick=10, $
  length=0.03, $
  charsize=1., $
  location=[0.30, 0.40]

cgPS_Close

if not keyword_Set(pngfile) then begin
  pngfile = 'fig1_gems_for.png'
endif

cgPS2Raster,'gems.ps', pngfile, /png, density = 1000

if not keyword_Set(scp_dest1) then begin
  scp_dest = 'soodal@164.125.38.179:/home/soodal/works/plot'
endif else begin
  scp_dest = scp_dest1
endelse

; send image to pc
if keyword_set(scp) then begin
  spawn, 'scp -P18742 ' + pngfile + $
  ' ' + scp_dest
endif

end
