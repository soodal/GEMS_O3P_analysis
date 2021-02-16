
o3 = ds_read_nc('../data/era5_global.nc', 'o3')
a = total(o3, 3)
lat = ds_read_nc('../data/era5_global.nc', 'latitude')
lon = ds_read_nc('../data/era5_global.nc', 'longitude')

ds_xymake2d, lon, lat, lon2d, lat2d

map = map('Orthographic', CENTER_LONGITUDE=125, CENTER_LATITUDE=20., /buffer)
grid = map.MAPGRID
grid.LINESTYLE='dotted'
grid.LABEL_POSITION=0
;grid.FONT_SIZE=14

 
ct = colortable(34)
;levels = findgen(51)/50.*(60-20) + 20
c = contour(a, lon2d, lat2d, /fill, /buffer, $
  ;c_value=levels, $
  overplot=map,$
  n_levels=100, $
  grid_units='degrees', rgb_table=ct)

m1 = mapcontinents()

pngfile = './plot/ERA5_o3_column_under300hpa.png'
;cb = colorbar(TITLE='Ozone Mass Mixing Ratio under 300 hPa[kg/kg]', $
  ;position=[0.2, 0.05, 0.8, 0.10], $
  ;;tickinterval=
  ;target=c)
c.save, pngfile

scp_dest = 'soodal@164.125.38.179:/home/soodal/WindowsHome/Pictures/scp'
spawn, 'scp -P18742 -p ' + pngfile + $
  ' ' + scp_dest
end
