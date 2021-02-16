pro  plot_gems_field_of_regards, lon, lat, pngfile=pngfile, scp_send=scp, $
  scp_dest=scp_dest1, title=pic_title

;fn = "/data1/L1C_GEMS/float/1x1/nopolc_eosrl/GK2_GEMS_L1C_20200616_0345_NOR_0694.4x4.nc4"
nominal = "/data1/L1C_GEMS/float/1x1/nopolc_eosrl/GK2_GEMS_L1C_20200616_0245_NOR_0694.nc"
full_central = "/data1/L1C_GEMS/float/1x1/nopolc_eosrl/GK2_GEMS_L1C_20200616_0545_NOR_0694.nc"
;half_east = "/data1/L1C_GEMS/float/1x1/nopolc/GK2_GEMS_L1C_20200803_2345_NOR_694.nc"

fns = [nominal]


if not keyword_Set(pic_title) then begin
  pic_title = 'Field of Regards'
endif


cgPS_Open, 'gems.ps'
cgDisplay, aspect=0.7

pos = [0.1,0.15,0.9,0.90]
charsize = (!D.Name EQ 'PS') ? cgDefCharsize()*0.7 : cgDefCharsize()*0.75
thick = (!D.Name EQ 'PS') ? 6 :3
londel = 10
latdel = 10

cgMap_set,20,120, Limit=[-60, 0, 60, 360],$
              pos=pos, $
              charsize=charsize,$
              /satellite, $
              sat_p=[2.2, 20, 0], $
              /iso,$
              /horizon, $
              title=pic_title
bluemarblefn = './data/july_blue_marble_new_generation_w_topo_world.topo.200407.3x5400x2700.jpg'

xstart = -180
ystart = -90
xsize = 5400
ysize = 2700
bm = image(bluemarblefn, $
  image_dimensions=[5400, 2700], $
  image_location=[-180, -90], $
  xrange=[-180, 180], $
  yrange=[-90, 90], $
  dimensions=[5400, 2700], $
  margin=0, $
  /overplot, $
  /buffer)
bm.getdata, bm_data

warped_image1 = Map_image(reform(bm_data[0, *, *]), xstart, ystart, $;-180, -90, $
  latmin=-90, $
  latmax=90, $
  lonmin=-180, $
  lonmax=180)
warped_image2 = Map_image(reform(bm_data[1, *, *]), xstart, ystart, $;-180, -90, $
  latmin=-90, $
  latmax=90, $
  lonmin=-180, $
  lonmax=180)
warped_image3 = Map_image(reform(bm_data[2, *, *]), xstart, ystart, $;-180, -90, $
  latmin=-90, $
  latmax=90, $
  lonmin=-180, $
  lonmax=180)

device=1
sz = size(warped_image1, /dimension)
a = fltarr(3, sz[0], sz[1])
a[0, *, *] = warped_image1
a[1, *, *] = warped_image2
a[2, *, *] = warped_image3
tv, a, true=3, xsize=xsize, ysize=ysize
;tv, warped_image2, 1, true=3, xsize=xsize, ysize=ysize
;tv, warped_image3, 2, true=3, xsize=xsize, ysize=ysize
map_continents, /coast, /hires, color=cgcolor('white')

;LoadCT, 22, Ncolors=254,bottom=1
dsloadct, 33



for i =0, n_elements(fns)-1 do begin
  lonlat = ds_read_gems_l1c_lonlat(fns[i])
  lon = lonlat.pixel_longitude
  lat = lonlat.pixel_latitude

  sz = size(lon)
  dim1 = sz[1]
  dim2 = sz[2]

  ;where(finite
  ;if nnan then begin
  xbox = [lon[0:dim1-1, 0], lon[dim1-1, 0:dim2-1], lon[dim1-1:0,dim2-1], lon[0:dim1-1, dim2-1]]
  ybox = [lat[0:dim1-1, 0], lat[dim1-1, 0:dim2-1], lat[dim1-1:0,dim2-1], lat[0:dim1-1, dim2-1]]

  oplot, xbox, ybox, color=cgcolor('r')
endif

cgMap_Grid, Color='white', GLinestyle=1,londel=londel, latdel=latdel, charsize=1 $
                 ,LONLAB=20 ,LATLAB=160, LABEL=1

cgMap_continents, color='white'


;cgColorbar, format='(f6.1)', range=range, $
  ;Position=[pos[0]+0.065, pos[1]-0.04, pos[2]-0.065, pos[1]-0.02], $
  ;Color='black',$
  ;OOB_High=253, $
  ;OOB_low=0,$
  ;divisions=5, $
  ;bottom=0, $
  ;ncolor=254, $
  ;minor=1,  $
  ;charsize=1.5,$
  ;XTicklen=1,  $
  ;XMinor=0, $
  ;AnnotateColor='black',  $
  ;Title=cb_title

cgPS_Close

if not keyword_Set(pngfile) then begin
  pngfile = 'gems_for.png'
endif

cgPS2Raster,'gems.ps', pngfile, /png, density = 1000

if not keyword_Set(scp_dest1) then begin
  scp_dest = 'soodal@164.125.38.179:/home/soodal/WindowsHome/Pictures/scp'
endif else begin
  scp_dest = scp_dest1
endelse

; send image to pc
if keyword_set(scp) then begin
  spawn, 'scp -P18742 -p ' + pngfile + $
  ' ' + scp_dest
endif



end
