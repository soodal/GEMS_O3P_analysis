
fn = '/data1/L2_GEMS/o3p/4x4/nopolc/eosrl/310-340/GEMS_O3P_20200616_0345_nopolc_2008111857.nc4'
basename = file_basename(fn)
data = ds_read_gems_l2_o3p(fn)

; for levels

O3 = data.O3
pres = data.PRESSURE

latitude = data.Latitude
longitude = data.Longitude

dim = n_elements(latitude)

lat = reform(latitude, dim)
lon = reform(longitude, dim)

pos = [0.10, 0.1, 0.8, 0.85]
leg_pos = [0.85, 0.3]
limit=[-5, 85, 55, 155]

nlayer = 24

outputpath = './plot/'
scppath = '/home/soodal/WindowsHome/Pictures/scp/'

for i=0, nlayer-1 do begin
  iPres = reform(pres[i, *, *])

  pres_geolocidx = where(lon gt 0 and lat gt 0 and iPres gt 0, /null)
  nanidx = where(iPres lt 0, /null)

  iPres[nanidx] = !values.f_nan
  meanpres = mean(iPres, /nan)
  strpres = string(meanpres, format='(f07.2)')

  map = MAP('Orthographic', /BUFFER, $
    limit=limit, $
    CENTER_LATITUDE=0, $
    CENTER_LONGITUDE=130)

  zrange = [min(iPres, /nan), max(iPres, /nan)]
  levels = findgen(zrange[1] - zrange[0] + 1, start=zrange[0])

  ct = colortable(74, /reverse)
  mag = bytscl(iPres[pres_geolocidx], min=zrange[0], max=zrange[1])
  p1 = scatterplot(lon[pres_geolocidx], lat[pres_geolocidx], $
    magnitude=mag, /buffer, $
    overplot=map, $
    rgb_table=ct, $
    clip=0, $
    position=pos, $
    symbol='dot', $
    sym_size=0.7, $
    title='Pressure for ' + strpres + 'hPa (mean)', $
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
  pngfile = outputpath + basename + '_1' + string(i, format='(i02)') + '_total_o3.png'

  p1.title.font_size=16
  p1.save, pngfile
  p1.close
  ; send image to pc
  spawn, 'scp -P18742 -p ' + pngfile + $
    ' soodal@164.125.38.179:'+scppath
endfor
end
