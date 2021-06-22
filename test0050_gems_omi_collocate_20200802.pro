
imonth = 8

smonth = string(imonth, format='(i02)')

year = 2020
month = 8
day = 16
hour = 03
minute = 45

;gemspath = '/data1/gems/o3p/works/GEMS_O3P/out/310340_KARI/'

pos = [0.10, 0.1, 0.95, 0.85]
leg_pos = [0.90, 0.80]

wininit = [300, 305, 310]
wininit = [310]

for iwininit=0, n_elements(wininit)-1 do begin
  winlim_str = string(wininit[iwininit], format='(i03)')
  gemspath = '/data2/L2_GEMS/val_1008/'
  gemspath = '/data1/gems/o3p/works/GEMS_O3P/out/paper/'
  for iday = 2, 2 do begin
    yyyy = string(year, format='(i04)')
    mm = string(month, format='(i02)')
    dd = string(iday, format='(i02)')
    hh = string(hour, format='(i02)')
    mi = string(minute, format='(i02)')
    gemsfn = 'GK2B_GEMS_L2_O3P_'+yyyy+mm+dd+'_'+hh+mi+'_winliminit' + $
      winlim_str + '_prec000000_climTB_b4x4_*.nc4'
    gemsfile = file_search(gemspath + gemsfn)
    if n_elements(gemsfile) ne 1 then begin
      stop
    endif
    colloc_omil2profoz_on_gems, year, month, iday, hour, minute, $
      output, gemsfile = gemsfile[0], /writegemsl2o3p
  endfor
endfor

end
