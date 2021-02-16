pro ds_read_merra2_inst3_3d_asm_np, fn, data

merra2asmfn = fn

EPV = ds_read_nc(merra2asmfn, 'EPV')
H = ds_read_nc(merra2asmfn, 'H')
lat = ds_read_nc(merra2asmfn, 'lat')
lev = ds_read_nc(merra2asmfn, 'lev')
lon = ds_read_nc(merra2asmfn, 'lon')
O3 = ds_read_nc(merra2asmfn, 'O3')
OMEGA = ds_read_nc(merra2asmfn, 'OMEGA')
PHIS = ds_read_nc(merra2asmfn, 'PHIS')
PS = ds_read_nc(merra2asmfn, 'PS')
QI = ds_read_nc(merra2asmfn, 'QI')
QL = ds_read_nc(merra2asmfn, 'QL')
QV = ds_read_nc(merra2asmfn, 'QV')
RH = ds_read_nc(merra2asmfn, 'RH')
SLP = ds_read_nc(merra2asmfn, 'SLP')
T = ds_read_nc(merra2asmfn, 'T')
time = ds_read_nc(merra2asmfn, 'time')
U = ds_read_nc(merra2asmfn, 'U')
V = ds_read_nc(merra2asmfn, 'V')


data = {EPV:EPV, $
  H:H, $ 
  lat:lat, $
  lev:lev, $
  lon:lon, $
  O3:O3, $
  OMEGA:OMEGA, $
  PHIS:PHIS, $
  PS:PS, $
  QI:QI, $
  QL:QL, $
  QV:QV, $
  RH:RH, $
  SLP:SLP, $
  T:T, $
  time:time, $
  U:U, $
  V:V}

return

end
