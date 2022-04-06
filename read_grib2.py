import pygrib

fn = '/data/private/soodal/g512_v070_ergl_unis_h024.2016101200.gb2'

grib = pygrib.open(fn)
for g in grib:
    print(g)
