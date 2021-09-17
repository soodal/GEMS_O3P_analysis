pro read_ozonesonde_pohang, filename, result


openr, lun, filename, /get_lun
ss = ''

readf, lun, ss, format='(A41)'
operator = ss
if strmid(ss, 0, 41) ne 'Korea Meteorological Administration (KMA)' then begin
  message, 'Error: Header format is different from expected.'
end

readf, lun, ss, format='(A48)'
title = ss
if strmid(ss, 0, 33) ne 'Ozone Sonde Vertical Profile Data' then begin
  message, 'Error: Dataname header format is different from expected.'
end

readf, lun, ss, format='(A24)'
station_name = strmid(ss, 19, 6)

if strmid(ss, 0, 7) ne 'Station' then begin
  message, 'Error: Station header format is different from expected.'
end

readf, lun, ss, format='(A26)'
latitude = float(strmid(ss, 18, 8))
if strmid(ss, 0, 8) ne 'Latitude' then begin
  message, 'Error: Latitude header format is different from expected.'
end

readf, lun, ss, format='(A27)'
longitude = float(strmid(ss, 18, 9))
if strmid(ss, 0, 9) ne 'Longitude' then begin
  message, 'Error: Longitude header format is different from expected.'
end

readf, lun, ss, format='(A28)'
if strmid(ss, 0, 11) ne 'Launch Date' then begin
  message, 'Error: Launch Date header format is different from expected.'
end
launch_year = fix(strmid(ss, 18, 4))

;print, strmid(ss, 23, 2)
launch_month = fix(strmid(ss, 23, 2))
;print, strmid(ss, 26, 2)
launch_day = fix(strmid(ss, 26, 2))

readf, lun, ss, format='(A30)'
if strmid(ss, 0, 11) ne 'Launch Time' then begin
  message, 'Error: Launch Time header format is different from expected.'
end
launch_hour = fix(strmid(ss, 18, 2))
launch_minute = fix(strmid(ss, 21, 2))
launch_second = fix(strmid(ss, 24, 2))

readf, lun, ss, format='(A77)'
if strmid(ss, 0, 11) ne '-----------' then begin
  message, 'Error: ----------- header format is different from expected.'
end
dummyline = ss

readf, lun, ss, format='(A76)'
if strmid(ss, 0, 9) ne 'TIME[m:s]' then begin
  message, 'Error: TIME[m:s] header format is different from expected.'
end

column_name = strsplit(ss, ', ', /extract)
readf, lun, ss, format='(A77)'
if strmid(ss, 0, 11) ne '-----------' then begin
  message, 'Error: ----------- header format is different from expected.'
end
dummyline = ss

free_lun, lun

readcol, filename, $
  time, pres, geop_h, temp, relhum ,ozone, $
  skipline=10, $
  delimiter=' ', $
  format='A,F,F,F,F,F'

;Korea Meteorological Administration (KMA)
;Ozone Sonde Vertical Profile Data - BIRM CY 1803
;Station         : Pohang
;Latitude        : 36.03259
;Longitude       : 129.37963
;Launch Date     : 2020.08.19
;Launch Time     : 05:16:00 UTC
;-----------------------------------------------------------------------------
;TIME[m:s], PRES[hPa], Altitude[m], GeoP_H[gpm], TEMP(degC], relhum[%], OZON[mPa]
;-----------------------------------------------------------------------------

nline = n_elements(time)
minarr = []
secarr = []
for i=0, nline-1 do begin
  ms = strsplit(time[i], ':', /extract)
  minute = fix(ms[0])
  second = fix(ms[1])
  minarr = [minarr, minute]
  secarr = [secarr, second]
endfor

jd = julday(launch_month, launch_day, launch_year, $
  launch_hour, launch_minute + minarr, launch_second + secarr)

result = create_struct('time_julday', jd, $
  'year', launch_year, $
  'month', launch_month, $
  'day', launch_day, $
  'hour', launch_hour, $
  'minute', launch_minute, $

  'PRES_hPa', pres, $
  ;'Altitude_m', altitude, $
  'GeoP_H_gpm', geop_h, $
  'TEMP_degC', temp, $
  'RH_perc', relhum, $
  'OZON_mPa', ozone, $
  'launch_latitude', latitude, $
  'launch_longitude', longitude)

return
end
