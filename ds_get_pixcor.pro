pro ds_get_pixcor, lon, lat, lonpixcor, latpixcor, nanval=nanval

dim = size(lon, /dim)
dimlat = size(lat, /dim)
if dim[0] ne dimlat[0] or dim[1] ne dimlat[1] then begin
  print, 'dimension of the longitude is different with the dimension of the latitude'
  stop
endif

lonpixcor = fltarr(dim[0], dim[1], 4)
latpixcor = fltarr(dim[0], dim[1], 4)

for ix=1, dim[0]-2 do begin
for iy=1, dim[1]-2 do begin
  
  lonll = (lon[ix, iy-1] + lon[ix-1, iy-1])/2.
  lonlr = (lon[ix+1, iy-1] + lon[ix, iy-1])/2.
  lonul = (lon[ix, iy+1] + lon[ix-1, iy+1])/2.
  lonur = (lon[ix+1, iy+1] + lon[ix, iy+1])/2.

  latll = (lat[ix, iy-1] + lat[ix-1, iy-1])/2.
  latlr = (lat[ix+1, iy-1] + lat[ix, iy-1])/2.
  latul = (lat[ix, iy+1] + lat[ix-1, iy+1])/2.
  latur = (lat[ix+1, iy+1] + lat[ix, iy+1])/2.
 
  latpixcor[ix, iy, *] = [latll, latlr, latur, latul]
  lonpixcor[ix, iy, *] = [lonll, lonlr, lonur, lonul]

ENDFOR
ENDFOR

; ix = 0
for iy=1, dim[1]-2 do begin
  lonll = lon[0, iy-1] - (lon[0, iy-1] - lon[1, iy-1])/2.
  lonlr = (lon[1, iy-1] + lon[0, iy-1])/2.
  lonul = lon[0, iy+1] - (lon[0, iy+1] - lon[1, iy+1])/2.
  lonur = (lon[1, iy+1] + lon[0, iy+1])/2.

  latll = lat[0, iy-1] - (lat[0, iy-1]-lat[1, iy-1])/2.
  latlr = (lat[0+1, iy-1] + lat[0, iy-1])/2.
  latul = lat[0, iy+1] - (lat[0, iy+1]-lat[1, iy-1])/2.
  latur = (lat[0+1, iy+1] + lat[0, iy+1])/2.

  latpixcor[ix, iy, *] = [latll, latlr, latur, latul]
  lonpixcor[ix, iy, *] = [lonll, lonlr, lonur, lonul]
endfor

; ix = -1
for iy=1, dim[1]-2 do begin
  lonll = (lon[-1, iy-1] + lon[-2, iy-1])/2.
  lonlr = (lon[-1, iy-1] + (lon[-2, iy-1] + lon[-1, iy-1])/2.
  lonul = (lon[-1, iy+1] + lon[-2, iy+1])/2.
  lonur = (lon[-1, iy+1] + (lon[-2, iy+1] + lon[-1, iy+1])/2.

  latll = (lat[-1, iy-1] + lat[-2, iy-1])/2.
  latlr = (lat[-1, iy-1] + (lat[-2, iy-1] + lat[-1, iy-1])/2.
  latul = (lat[-1, iy+1] + lat[-2, iy+1])/2.
  latur = (lat[-1, iy+1] + (lat[-2, iy+1] + lat[-1, iy+1])/2.

  latpixcor[ix, iy, *] = [latll, latlr, latur, latul]
  lonpixcor[ix, iy, *] = [lonll, lonlr, lonur, lonul]
endfor

; iy = 0
for ix = 1 dim[0]-2 do begin
  lonll = (lon[ix, 0] + lon[ix-1, 0])/2.
  lonlr = (lon[ix+1, 0] + lon[ix, 0])/2.
  lonul = (lon[ix, 1] + lon[ix-1, 1])/2.
  lonur = (lon[ix+1, 1] + lon[ix, 1])/2.

  latll = (lat[ix, 0] + lat[ix-1, 0])/2.
  latlr = (lat[ix+1, 0] + lat[ix, 0])/2.
  latul = (lat[ix, 1] + lat[ix-1, 1])/2.
  latur = (lat[ix+1, 1] + lat[ix, 1])/2.

  latpixcor[ix, iy, *] = [latll, latlr, latur, latul]
  lonpixcor[ix, iy, *] = [lonll, lonlr, lonur, lonul]
endfor
 
; iy = -1
for ix = 1 dim[0]-2 do begin
  lonll = (lon[ix, -2] + lon[ix-1, -2])/2.
  lonlr = (lon[ix+1, -2] + lon[ix, -2])/2.
  lonul = (lon[ix, -2] + lon[ix-1, -2])/2.
  lonur = (lon[ix+1, -2] + lon[ix, -2])/2.

  latll = (lat[ix, -2] + lat[ix-1, -2])/2.
  latlr = (lat[ix+1, -2] + lat[ix, -2])/2.
  latul = (lat[ix, -2] + lat[ix-1, -2])/2.
  latur = (lat[ix+1, -2] + lat[ix, -2])/2.

  latpixcor[ix, iy, *] = [latll, latlr, latur, latul]
  lonpixcor[ix, iy, *] = [lonll, lonlr, lonur, lonul]
endfor

; ix = 0, iy = 0
ix = 0
iy = 0
lonll = lon[0, 0]
lonlr = (lon[1, 0] + lon[0, 0])/2.
lonul = lon[0, 1]
lonur = (lon[1, 1] + lon[0, 1])/2.

latll = lat[0, 0]
latlr = (lat[1, 0] + lat[0, 0])/2.
latul = lat[0, 1]
latur = (lat[1, 1] + lat[0, 1])/2.
latpixcor[ix, iy, *] = [latll, latlr, latur, latul]
lonpixcor[ix, iy, *] = [lonll, lonlr, lonur, lonul]


; ix = 0, iy = -1
ix = 0
iy = -1
lonll = lon[0, -2]
lonlr = (lon[1, -2] + lon[0, -2])/2.
lonul = lon[0, -1]
lonur = (lon[1, -1] + lon[0, -1])/2.

latll = lat[0, -2]
latlr = (lat[1, -2] + lat[0, -2])/2.
latul = lat[0, -1]
latur = (lat[1, -1] + lat[0, -1])/2.

latpixcor[ix, iy, *] = [latll, latlr, latur, latul]
lonpixcor[ix, iy, *] = [lonll, lonlr, lonur, lonul]

; ix = -1, iy = 0
ix = -1
iy = 0
lonll = (lon[-1, 0] + lon[-2, 0])/2.
lonlr = lon[-1, 0] 
lonul = (lon[-1, 1] + lon[-2, 1])/2.
lonur = lon[-1, 1]

latll = (lat[-1, 0] + lat[-2, 0])/2.
latlr = lat[-1, 0] 
latul = (lat[-1, 1] + lat[-2, 1])/2.
latur = lat[-1, 1]
 
latpixcor[ix, iy, *] = [latll, latlr, latur, latul]
lonpixcor[ix, iy, *] = [lonll, lonlr, lonur, lonul]

; ix = -1, iy = -1
ix = -1
iy = -1
lonll = (lon[-1, -2] + lon[-2, -2])/2.
lonlr = (lon[-1, -2] + lon[-1, -2])/2.
lonul = (lon[-1, -1] + lon[-2, -1])/2.
lonur = (lon[-1, -1] + lon[-1, -1])/2.

latll = (lat[-1, -2] + lat[-2, -2])/2.
latlr = (lat[-1, -2] + lat[-1, -2])/2.
latul = (lat[-1, -1] + lat[-2, -1])/2.
latur = (lat[-1, -1] + lat[-1, -1])/2.
 
latpixcor[ix, iy, *] = [latll, latlr, latur, latul]
lonpixcor[ix, iy, *] = [lonll, lonlr, lonur, lonul]

return
end
