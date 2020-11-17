o3 = ds_read_nc('../data/era5.nc', 'o3')
a = total(o3, 3)
lat = ds_read_nc('../data/era5.nc', 'latitude')
lon = ds_read_nc('../data/era5.nc', 'longitude')

ds_xymake2d, lon, lat, lon2d, lat2d

c = contour(lon2d, lat2d, a, /fill, /buffer)
pngfile = './plot/ERA5_o3_column_under300hpa.png'
c.save, pngfile

scp_dest = 'soodal@164.125.38.179:/home/soodal/WindowsHome/Pictures/scp'
spawn, 'scp -P18742 -p ' + pngfile + $
  ' ' + scp_dest
end
