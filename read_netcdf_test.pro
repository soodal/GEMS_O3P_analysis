fid=ncdf_open(file)

groupsl1id=ncdf_groupsinq(fid)

ngroup = n_elements(groupsl1id)

for i =0, ngroup-1
  gid=groupsl1id[i]
  varids=ncdf_varidsinq(gid)
  for vn=0, n_elements(varids) - 1 do begin
      ncdf_varget,fid,varids[vn], data
  ENDFOR
ENDFOR

end
