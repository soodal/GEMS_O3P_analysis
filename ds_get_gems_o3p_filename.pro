function ds_get_gems_o3p_filename, year, month, day, hour, minute
ysz = size(year)
if ysz[1] eq 2 then begin
  yyyy = string(year, format='(i04)')
endif else begin
  yyyy = year
endelse

msz = size(month)
if msz[1] eq 2 then begin
  mm = string(month, format='(i02)')
endif else begin
  mm = month
endelse

dsz = size(day)
if dsz[1] eq 2 then begin
  dd = string(day, format='(i02)')
endif else begin
  dd = day
endelse

hsz = size(hour)
if hsz[1] eq 2 then begin
  hh = string(hour, format='(i02)')
endif else begin
  hh = hour
endelse

misz = size(minute)
if misz[1] eq 2 then begin
  mi = string(minute, format='(i02)')
endif else begin
  mi = minute
endelse

gemsfn = '/data2/L2_GEMS/val_1008/GK2_GEMS_O3P_' $
  + yyyy + mm + dd + '_' + hh + mi + '.nc4'

gemsfn_v1ba = '/data2/L2_GEMS/val_1008/GK2_GEMS_O3P_' $
  + yyyy + mm + dd + '_' + hh + mi + '_v1.ba.nc4'

if year eq 2021 and month ge 4 and month le 6 then begin
  gemsfn = '/data2/L2_GEMS.nier/monthly_val/O3P/GK2_GEMS_L2_' $
    + yyyy + mm + dd + '_' + hh + mi + '_O3P_*.nc'
  fl = file_search(gemsfn)
  gemsfn = fl[0]
  return, gemsfn
endif

if year eq 2021 and month ge 7 then begin
  gemsfn = '/data/nier_ftp/O3P/V03/'+ yyyy + mm + '/' + dd + '/' + $
    'GK2_GEMS_L2_' + yyyy + mm + dd + '_' + hh + mi + '_O3P_*BIN4x4.nc'
  fl = file_search(gemsfn)
  gemsfn = fl[0]
  return, gemsfn
endif


if file_test(gemsfn_v1ba) then begin
  return, gemsfn_v1ba
endif else if file_test(gemsfn) then begin
  return, gemsfn
endif else begin
  return, !null
endelse



end
