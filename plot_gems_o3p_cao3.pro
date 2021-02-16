
fn = '../GEMS_O3P_Yonsei/out/310340_KARI/GK2B_GEMS_L2_20160115_0300_O3P_ND_DPRO_first.4x4_2008102141.nc4'

data = ds_read_gems_l2_o3p(fn)
;hs_read_gems2003, fn



pos = [0.10, 0.1, 0.8, 0.85]
leg_pos = [0.85, 0.3]
limit=[-5, 85, 55, 155]


cao3 = data.ColumnAmountO3
o3total = reform(cao3[0, *, *])
o3strato = reform(cao3[1, *, *])
o3tropo = reform(cao3[2, *, *])

lat = data.Latitude
lon = data.Longitude


nanidx = where(lon < 0, /null)
lon[nanidx] = !values.f_nan
nanidx = where(lat < 0, /null)
lat[nanidx] = !values.f_nan

nanidx = where(o3total < 0, /null)
o3total[nanidx] = !values.f_nan

nanidx = where(o3tropo < 0, /null)
o3tropo[nanidx] = !values.f_nan

nanidx = where(o3strato < 0, /null)
o3strato[nanidx] = !values.f_nan

; for CAO3

map = MAP('Orthographic', /BUFFER, $
  limit=limit, $
  CENTER_LATITUDE=25, $
  CENTER_LONGITUDE=125)

levels = findgen(380-250 + 1, start=250)

ct = colortable(74)
p1 = contour(o3total, lon, lat, /buffer, $
  /fill, $
  overplot=map, $
  grid_units='degrees', $
  rgb_table=ct, $
  c_value=levels, $
  clip=0, $
  position=pos, $
  title='Total Column Amount O3 KARI L1C 300-340 nm fitting', $
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
     
cb = colorbar(target=[p1], orientation=1, $
  position=[0.93, 0.1, 0.97, 0.85], $
  ;tickdir=1, $
  ;tickinterval=10, $
  ;ticklen=0.5, $
  ;border_on=1, $
  color='black', $
  title='O3[DU]', $
  major=10, $
  ;minor=0, $
  taper=1)

caldat, systime(/julian), mo, dd, yyyy, hh, mi, ss
isodate = string(yyyy, format='(I04)') + '-' +string(mo, format='(i02)') $
  + '-' + string(dd, format='(i02)') + 'T' + string(hh, format='(i02)') $
  + ':' + string(mi, format='(i02)') + ':' + string(ss, format='(i02)')

t = text(0.6, 0.05, 'Dae Sung Choi, ' + isodate)

pngfile = './plot/gems_l2_cao3_300340_KARI.png'

p1.title.font_size=16
p1.save, pngfile
p1.close
; send image to pc
spawn, 'scp -P18742 -p ' + pngfile + $
  ' soodal@164.125.38.179:/home/soodal/WindowsHome/Pictures/scp'

; for Stratospheric O3

map = MAP('Orthographic', /BUFFER, $
  limit=limit, $
  CENTER_LATITUDE=25, $
  CENTER_LONGITUDE=125)

levels = findgen(320-190+1, start=190)

ct = colortable(74)
p1 = contour(o3strato, lon, lat, /buffer, $
  /fill, $
  overplot=map, $
  grid_units='degrees', $
  rgb_table=ct, $
  c_value=levels, $
  clip=0, $
  position=pos, $
  title='Stratospheric Column Amount O3 KARI L1C 300-340 nm fitting', $
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
     
cb = colorbar(target=[p1], orientation=1, $
  position=[0.93, 0.1, 0.97, 0.85], $
  ;tickdir=1, $
  ;tickinterval=10, $
  ;ticklen=0.5, $
  ;border_on=1, $
  color='black', $
  title='O3[DU]', $
  major=10, $
  ;minor=0, $
  taper=1)


caldat, systime(/julian), mo, dd, yyyy, hh, mi, ss
isodate = string(yyyy, format='(I04)') + '-' +string(mo, format='(i02)') $
  + '-' + string(dd, format='(i02)') + 'T' + string(hh, format='(i02)') $
  + ':' + string(mi, format='(i02)') + ':' + string(ss, format='(i02)')

t = text(0.6, 0.05, 'Dae Sung Choi, ' + isodate)
pngfile = './plot/gems_l2_o3strato_300340_KARI.png'

p1.title.font_size=16
p1.save, pngfile
p1.close
; send image to pc
spawn, 'scp -P18742 -p ' + pngfile + $
  ' soodal@164.125.38.179:/home/soodal/WindowsHome/Pictures/scp'

; for Tropospheric O3

map = MAP('Orthographic', /BUFFER, $
  limit=limit, $
  CENTER_LATITUDE=25, $
  CENTER_LONGITUDE=125)

levels = findgen(60-20+1, start=20)

ct = colortable(74)
p1 = contour(o3tropo, lon, lat, /buffer, $
  /fill, $
  overplot=map, $
  grid_units='degrees', $
  rgb_table=ct, $
  c_value=levels, $
  clip=0, $
  position=pos, $
  title='Troposheric Column Amount O3 KARI L1C 300-340 nm fitting', $
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
     
cb = colorbar(target=[p1], orientation=1, $
  position=[0.93, 0.1, 0.97, 0.85], $
  color='black', $
  title='O3[DU]', $
  major=10, $
  ;minor=0, $
  taper=1)


caldat, systime(/julian), mo, dd, yyyy, hh, mi, ss
isodate = string(yyyy, format='(I04)') + '-' +string(mo, format='(i02)') $
  + '-' + string(dd, format='(i02)') + 'T' + string(hh, format='(i02)') $
  + ':' + string(mi, format='(i02)') + ':' + string(ss, format='(i02)')

t = text(0.6, 0.05, 'Dae Sung Choi, ' + isodate)
pngfile = './plot/gems_l2_o3tropo_300340_KARI.png'

p1.title.font_size=16
p1.save, pngfile
p1.close
; send image to pc
spawn, 'scp -P18742 -p ' + pngfile + $
  ' soodal@164.125.38.179:/home/soodal/WindowsHome/Pictures/scp'


end
