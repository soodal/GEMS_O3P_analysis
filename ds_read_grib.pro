pro ds_read_grib, filename, data
grib_list, filename
nr = grib_count(filename)
print, nr


for ir = 1, nr do begin
  data = grib_getdata(filename, ir)
  help, data
endfor

return
end
