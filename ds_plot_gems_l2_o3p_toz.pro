pro ds_plot_gems_l2_o3p_toz, gemsl2fn, outputpath = outputpath, scppath = scppath

; Set keyword
if not keyword_set(outputpath) then begin
  outputpath = './plot/'
endif

if not keyword_set(outputpath) then begin
  scppath = ''
endif
scppath = '/home/soodal/WindowsHome/Pictures/scp/' + scppath


; Set parameters
pos = [0.10, 0.1, 0.8, 0.85]
leg_pos = [0.85, 0.3]
limit=[-5, 85, 55, 155]
nlayer = 24

basename = file_basename(gemsl2fn)

data = ds_read_gems_l2_o3p(gemsl2fn)

cao3 = data.ColumnAmountO3
o3total = reform(cao3[0, *, *], 174, 512)
o3strato = reform(cao3[1, *, *], 174, 512)
o3tropo = reform(cao3[2, *, *], 174, 512)

latitude = data.Latitude
longitude = data.Longitude
dim = n_elements(latitude)

lat = reform(latitude, dim)
lon = reform(longitude, dim)

o3total_geolocidx = where(lon gt 0 and lat gt 0 and o3total gt 0, /null)
o3tropo_geolocidx = where(lon gt 0 and lat gt 0 and o3tropo gt 0, /null)
o3strato_geolocidx = where(lon gt 0 and lat gt 0 and o3strato gt 0, /null)

; for CAO3

map = MAP('Orthographic', /BUFFER, $
  limit=limit, $
  CENTER_LATITUDE=0, $
  CENTER_LONGITUDE=130)

zrange = [250, 380]
levels = findgen(zrange[1] - zrange[0] + 1, start=zrange[0])

ct = colortable(74, /reverse)
mag = bytscl(o3total[o3total_geolocidx], min=zrange[0], max=zrange[1])
p1 = scatterplot(lon[o3total_geolocidx], lat[o3total_geolocidx], $
  magnitude=mag, /buffer, $
  overplot=map, $
  rgb_table=ct, $
  clip=0, $
  position=pos, $
  symbol='dot', $
  sym_size=0.7, $
  title='Total Column Amount O3 KARI L1C 310-340 nm fitting', $
  axis_style=0)

grid = map.MAPGRID
grid.box_axes=1
grid.LINESTYLE = "dashed"
grid.LABEL_POSITION = 0
map['latitude'].LABEL_ANGLE=80
map['longitude'].LABEL_ANGLE=0
grid.longitude_min=limit[1]
grid.latitude_min=limit[0]
grid.grid_Latitude=10
grid.grid_Longitude=10
grid.FONT_SIZE=10
 
m1 = MAPCONTINENTS(/hires)
     
tickvalues = linspace(zrange[0], zrange[1], 10)
cb = colorbar( $
  orientation=1, $
  position=[0.93, 0.1, 0.97, 0.85], $
  tickvalues=tickvalues,$
  range=[zrange[0], zrange[1]], $
  color='black', $
  title='O3[DU]', $
  taper=1)

caldat, systime(/julian), mo, dd, yyyy, hh, mi, ss
isodate = string(yyyy, format='(I04)') + '-' +string(mo, format='(i02)') $
  + '-' + string(dd, format='(i02)') + 'T' + string(hh, format='(i02)') $
  + ':' + string(mi, format='(i02)') + ':' + string(ss, format='(i02)')

t = text(0.4, 0.025, 'Dae Sung Choi, ' + isodate)
pngfile = outputpath + basename + '_toz.png'

p1.title.font_size=16
p1.save, pngfile
p1.close
; send image to pc
spawn, 'scp -P18742 -p ' + pngfile + $
  ' soodal@164.125.38.179:'+scppath

end
