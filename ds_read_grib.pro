pro ds_read_grib, filename
grib_list, filename
nr = grib_count(filename)
print, nr

return
data = grib_getdata(filename, rnum)

end
