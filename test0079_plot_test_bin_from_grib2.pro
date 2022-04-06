fn2016 = '/data/private/soodal/dist2016.bin'
openr, lun, fn2016, /get_lun

a2016 = fltarr(1024, 769)
readu, lun, a2016
free_lun, lun

lon = findgen(1024) / 1024 * 360
lat = findgen(769)/ 769. * 180 - 90

ds_xymake2d, lon, lat, lon2d, lat2d

;c = contour(a, lon2d, lat2d, /fill)

fn2017 = '/data/private/soodal/dist2017.bin'
openr, lun, fn2017, /get_lun

a2017 = fltarr(1024, 769)
readu, lun, a2017
free_lun, lun

showme, a2017-a2016
end

