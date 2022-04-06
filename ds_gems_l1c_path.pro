pro ds_gems_l1c_path, year, month, day, hour, minute, path, filename

yyyy = string(year, format='(i04)')
mm = string(month, format='(i02)')
dd = string(day, format='(i02)')
hh = string(hour, format='(i02)')
mi = string(minute, format='(i02)')

if year eq 2021 then begin 
  gemsfp = '/data2/L1C_GEMS/L1C_nier/L1C/' + yyyy + mm + '/' + dd + '/'
  gemsfn = 'GK2_GEMS_O3P_' $
    + yyyy + mm + dd + '_' + hh + mi + '_v1.ba.nc4'
  filelist = file_search(gemsfp+gemsfn)

endif


path = gemsfp


end

