pro read_ozonesonde_pohang, filename, result


openr, lun, filename, /get_lun
ss = ''

readf, lun, ss, format='(A41)'
operator = ss

readf, lun, ss, format='(A48)'
title = ss

readf, lun, ss, format='(A24)'
station_name = strmid(ss, 19, 6)

readf, lun, ss, format='(A26)'
latitude = float(strmid(ss, 18, 8))

readf, lun, ss, format='(A27)'
longitude = float(strmid(ss, 18, 9))

readf, lun, ss, format='(A28)'
launch_year = fix(strmid(ss, 18, 4))
;print, strmid(ss, 23, 2)
launch_month = fix(strmid(ss, 23, 2))
;print, strmid(ss, 26, 2)
launch_day = fix(strmid(ss, 26, 2))

readf, lun, ss, format='(A30)'
launch_hour = fix(strmid(ss, 19, 2))
launch_minute = fix(strmid(ss, 22, 2))
launch_second = fix(strmid(ss, 25, 2))

readf, lun, ss, format='(A77)'
dummyline = ss
readf, lun, ss, format='(A76)'
column_name = strsplit(ss, ', ', /extract)
readf, lun, ss, format='(A77)'
dummyline = ss

free_lun, lun

readcol, filename, $
  time, pres, altitude, geop_h, temp, rh ,ozone, $
  skipline=10, $
  format='A,F,F,F,F,F,F'

;Korea Meteorological Administration (KMA)
;Ozone Sonde Vertical Profile Data - BIRM CY 1803
;Station         : Pohang
;Latitude        : 36.03259
;Longitude       : 129.37963
;Launch Date     : 2020.08.19
;Launch Time     : 05:16:00 UTC
;-----------------------------------------------------------------------------
;TIME[m:s], PRES[hPa], Altitude[m], GeoP_H[gpm], TEMP(degC], RH[%], OZON[mPa]
;-----------------------------------------------------------------------------

nline = n_elements(time)
minarr = []
secarr = []
ms = strsplit(time, ':', /extract)
for i=0, nline-1 do begin
  minute = fix(ms[i, 0])
  second = fix(ms[i, 1])
  minarr = [minarr, minute]
  secarr = [secarr, second]
endfor

jd = julday(launch_month, launch_day, launch_year, $
  launch_minute + minarr, launch_second + secarr)

result = create_struct('time_julday', jd, $
  'PRES_hPa', pres, $
  'Altitude_m', altitude, $
  'GeoP_H_gpm', geop_h, $
  'TEMP_degC', temp, $
  'RH_perc', rh, $
  'OZON_mPa', ozone, $
  'launch_latitude', latitude, $
  'launch_longitude', longitude)

return
end
