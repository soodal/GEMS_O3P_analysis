function ds_read_merra2_inst3_3d_chm_nv, fn

merra2asmfn = fn

AIRDENS = ds_read_nc(merra2asmfn, 'AIRDENS')
CO = ds_read_nc(merra2asmfn, 'CO')
DELP = ds_read_nc(merra2asmfn, 'DELP')
O3 = ds_read_nc(merra2asmfn, 'O3')
PS = ds_read_nc(merra2asmfn, 'PS')
lat = ds_read_nc(merra2asmfn, 'lat')
lev = ds_read_nc(merra2asmfn, 'lev')
lon = ds_read_nc(merra2asmfn, 'lon')
time = ds_read_nc(merra2asmfn, 'time')


data = {AIRDENS:AIRDENS, $
  CO:CO, $ 
  DELP:DELP, $
  O3:O3, $
  PS:PS, $
  lat:lat, $
  lev:lev, $
  lon:lon, $
  time:time}

return, data

end
