function ds_read_merra2_tavg3_3d_asm_nv_collocated_on_gems_l1c, fn, verbose=verbose
merra2asmfn = fn

if not keyword_set(verbose) then begin
  verbose = 0
endif


if verbose then begin
  print, 'reading CLOUD from the netcdf file.'
endif
CLOUD = ds_read_nc(merra2asmfn, 'CLOUD')
if verbose then begin
  print, 'reading DELP from the netcdf file.'
endif
DELP = ds_read_nc(merra2asmfn, 'DELP')
if verbose then begin
  print, 'reading EPV from the netcdf file.'
endif
EPV = ds_read_nc(merra2asmfn, 'EPV')
if verbose then begin
  print, 'reading H from the netcdf file.'
endif
H = ds_read_nc(merra2asmfn, 'H')
if verbose then begin
  print, 'reading O3 from the netcdf file.'
endif
O3 = ds_read_nc(merra2asmfn, 'O3')
if verbose then begin
  print, 'reading OMEGA from the netcdf file.'
endif
OMEGA = ds_read_nc(merra2asmfn, 'OMEGA')
if verbose then begin
  print, 'reading PHIS from the netcdf file.'
endif
PHIS = ds_read_nc(merra2asmfn, 'PHIS')
if verbose then begin
  print, 'reading PL from the netcdf file.'
endif
PL = ds_read_nc(merra2asmfn, 'PL')
if verbose then begin
  print, 'reading PSD from the netcdf file.'
endif
PS = ds_read_nc(merra2asmfn, 'PS')
if verbose then begin
  print, 'reading QI from the netcdf file.'
endif
QI = ds_read_nc(merra2asmfn, 'QI')
if verbose then begin
  print, 'reading QL from the netcdf file.'
endif
QL = ds_read_nc(merra2asmfn, 'QL')
if verbose then begin
  print, 'reading QV from the netcdf file.'
endif
QV = ds_read_nc(merra2asmfn, 'QV')
if verbose then begin
  print, 'reading RH from the netcdf file.'
endif
RH = ds_read_nc(merra2asmfn, 'RH')
if verbose then begin
  print, 'reading SLP from the netcdf file.'
endif
SLP = ds_read_nc(merra2asmfn, 'SLP')
if verbose then begin
  print, 'reading T from the netcdf file.'
endif
T = ds_read_nc(merra2asmfn, 'T')
if verbose then begin
  print, 'reading U from the netcdf file.'
endif
U = ds_read_nc(merra2asmfn, 'U')
if verbose then begin
  print, 'reading V from the netcdf file.'
endif
V = ds_read_nc(merra2asmfn, 'V')
if verbose then begin
  print, 'reading lat from the netcdf file.'
endif
lat = ds_read_nc(merra2asmfn, 'lat')
if verbose then begin
  print, 'reading lev from the netcdf file.'
endif
lev = ds_read_nc(merra2asmfn, 'lev')
if verbose then begin
  print, 'reading lon from the netcdf file.'
endif
lon = ds_read_nc(merra2asmfn, 'lon')
if verbose then begin
  print, 'reading time from the netcdf file.'
endif
;time = ds_read_nc(merra2asmfn, 'time')


data = {CLOUD:CLOUD, $
  DELP:DELP, $
  EPV:EPV, $
  H:H, $ 
  O3:O3, $
  OMEGA:OMEGA, $
  PHIS:PHIS, $
  PL:PL, $
  PS:PS, $
  QI:QI, $
  QL:QL, $
  QV:QV, $
  RH:RH, $
  SLP:SLP, $
  T:T, $
  U:U, $
  V:V, $
  lat:lat, $
  lev:lev, $
  lon:lon};, $
  ;time:time}

if verbose then begin
  print, 'returning result struct.'
endif
return, data

end
