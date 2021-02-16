image = bytscl(sin(dist(400)/10))

TV, image
latmin = -65
latmax = 65
lonmin = 160
lonmax = -70 + 360

map_set, 0, -140, /orthographic, /isotropic, $
  limit=[latmin, lonmin, latmax, lonmax]

result = map_image(image, startx, starty, compress=1, $
  latmin=latmin, lonmin=lonmin, $
  latmax=latmax, lonmax=lonmax)

TV, result, startx, starty
map_grid, latdel=10, londel=10, /label, /horizon
map_continents, /coasts

end
