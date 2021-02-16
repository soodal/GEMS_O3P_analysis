Function h5read, filename, path, tr=tr
fi = file_info(filename)
if fi.size le 2000 then BEGIN
  print, 'File size is too small. Please check the file.'
  return, 0
endif
fid = h5f_open(filename)
did = h5d_open(fid,path)
data = h5d_read(did)
h5d_close, did
h5f_close, fid
if keyword_set(tr) then BEGIN
  data = transpose(data) 
endif
return, data
end
