pro ds_plot_gems_l2_o3t_scene, gemsfile, outputpath=outputpath, project=project, sub=sub
; plotting every scene aug to oct for gems l2 o3p
;
if not keyword_set(outputpath) then begin
  outputpath = './plot/'
endif

if keyword_set(project) then begin
  outputpath = outputpath + project + '/'
endif

if keyword_set(sub) then begin
  outputpath = outputpath + sub + '/'
endif
  
;scp_dest = 'soodal@164.125.38.179:/home/soodal/works/plot/paper/2021_initial_report'

;initial variable setting

;path = '/data/private/soodal/vartest/climaTB_FNL/'

;varlist = ['EffectiveCloudFraction', $
           ;'O3' , $
           ;;'O3Apriori', 'O3AprioriError', $
           ;;'CloudPressure', $
           ;;'SimulatedRadiances', $
           ;;'O3Apriori', 'O3AprioriError',$
           ;'ColumnAmountO3', $
           ;'Latitude' ,'Longitude', $
           ;'Time','Altitude' ,    $
           ;;'Pressure', 'TropopausePressure', $
           ;;'Wavelengths', $
           ;'WavelengthsWholeRange']

; run

;jd_list = timegen(start=julday(6, 21, 2021, 3, 45), $
  ;final=julday(6, 21, 2021, 3, 45), units='Hours')


datetimepos = stregex(gemsfile, '[[:digit:]]{8}_[[:digit:]]{4}')

print, gemsfile
print, datetimepos
year = fix(strmid(gemsfile, datetimepos, 4))
print, strmid(gemsfile, datetimepos, 4)
month = fix(strmid(gemsfile, datetimepos+4, 2))
print, strmid(gemsfile, datetimepos+4, 2)
day = fix(strmid(gemsfile, datetimepos+6, 2))
print, strmid(gemsfile, datetimepos+6, 2)
hour = fix(strmid(gemsfile, datetimepos+9, 2))
print, strmid(gemsfile, datetimepos+9, 2)
minute = fix(strmid(gemsfile, datetimepos+11, 2))
print, strmid(gemsfile, datetimepos+11, 2)

;caldat, jd_list[idate], month, day, year, hour, minute
yyyy = string(year, format='(i04)')
mm = string(month, format='(i02)')
dd = string(day, format='(i02)')
hh = string(hour, format='(i02)')
mi = string(minute, format='(i02)')

datetime_str = yyyy + mm + dd + '_' + hh + mi
;search_str = 'GK2_GEMS_L2_O3P_' + datetime_str + '_TB_BIN4x4.nc'
;print, path+search_str

outputpath = outputpath + datetime_str + '/'

if not file_test(outputpath) then begin
  file_mkdir, outputpath
endif

fileok = file_test(gemsfile)


if fileok then begin 
  gemsl2o3t = ds_read_gems_l2_o3t(gemsfile)

  ;ecf = gemsl2o3t.EffectiveCloudFraction

  ;ecf02idx = where(ecf gt 0.2, /null)
  ;ecfnanidx = where(ecf lt -1e29, /null)

  longitude = gemsl2o3t.longitude
  latitude = gemsl2o3t.latitude


  nanidx = where(longitude lt -500, /null)
  longitude[nanidx] = !values.f_nan
  latitude[nanidx] = !values.f_nan

  nanidx = where(latitude lt -500, /null)
  longitude[nanidx] = !values.f_nan
  latitude[nanidx] = !values.f_nan

  cao3 = gemsl2o3t.ColumnAmountO3

  rangetotal = [220, 430]

  sz = size(cao3, /dim)

  if not file_test(outputpath) then begin
    file_mkdir, outputpath
  endif

  nanidx = where(cao3 lt -1.0e29, /null)
  cao3[nanidx] = !values.f_nan
  plot_sat_proj, cao3, longitude, latitude, $
    title='GEMS O3T Total ozone ' + datetime_str, $
    range=rangetotal, $
    pngfile=outputpath + '01_totalozone.png';, $
    ;/scp_send, $
    ;scp_dest=scp_dest

endif else begin
  print, 'Can not find any file.'
endelse


end
