pro ds_plot_gemso3p_scene, gemsfile, outputpath=outputpath, project=project, sub=sub
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

varlist = ['EffectiveCloudFractionUV', 'ProcessingQualityFlags', $
           'O3' , $
           ;'O3Apriori', 'O3AprioriError', $
           ;'CloudPressure', $
           ;'SimulatedRadiances', $
           ;'O3Apriori', 'O3AprioriError',$
           'ColumnAmountO3', $
           'Latitude' ,'Longitude', $
           'Time','Altitude' ,    $
           ;'Pressure', 'TropopausePressure', $
           ;'Wavelengths', $
           'WavelengthsWholeRange']

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
  o3p = ds_read_gems_l2_o3p(gemsfile, varlist=varlist)

  ecf = o3p.EffectiveCloudFractionUV

  ecf02idx = where(ecf gt 0.2, /null)
  ecfnanidx = where(ecf lt 0, /null)

  ds_gems_l2o3p_accum, o3p, o3under300hpa, hpa=300.
  ds_gems_l2o3p_accum, o3p, o3under500hpa, hpa=500.

  longitude = o3p.longitude
  latitude = o3p.latitude

  nanidx = where(o3under300hpa lt 0, /null)

  ;o3under300hpa[ecf02idx] = !values.f_nan
  o3under300hpa[ecfnanidx] = !values.f_nan
  o3under300hpa[nanidx] = !values.f_nan

  nanidx = where(longitude lt -500, /null)
  longitude[nanidx] = !values.f_nan
  latitude[nanidx] = !values.f_nan

  nanidx = where(latitude lt -500, /null)
  longitude[nanidx] = !values.f_nan
  latitude[nanidx] = !values.f_nan

  trppsprs = o3p.TropopausePressure



  cao3 = o3p.ColumnAmountO3
  ;if mondate eq '0806' then begin
    ;rangetotal = [250, 310]
  ;endif else begin
    ;rangetotal = [250, 345]
  ;ENDELSE

  rangetotal = [220, 430]

  sz = size(cao3, /dim)
  if sz[0] eq 3 then begin
    cao3 = transpose(cao3, [1, 2, 0])
  endif

  o3 = o3p.o3
  sz = size(o3, /dim)
  if sz[0] eq 24 then begin
    o3 = transpose(o3, [1, 2, 0])
  endif
    

  if not file_test(outputpath) then begin
    file_mkdir, outputpath
  endif

  nanidx = where(cao3 lt -1.0e20, /null)
  cao3[nanidx] = !values.f_nan
  plot_sat_proj, reform(cao3[*, *, 0]), longitude, latitude, $
    title='GEMS O3P Total ozone ' + datetime_str, $
    range=rangetotal, $
    pngfile=outputpath + '01_totalozone.png';, $
    ;/scp_send, $
    ;scp_dest=scp_dest

  plot_sat_proj, reform(cao3[*, *, 2]), longitude, latitude, $
    title='GEMS O3P Tropospheric ozone ' + datetime_str, $
    ;range=[20, 100], $
    range=[0, 60], $
    pngfile=outputpath + '03_tropsphericozone.png';, $

  nanidx = where(o3under300hpa lt -1.0e20, /null)
  o3under300hpa[nanidx] = !values.f_nan
  plot_sat_proj, o3under300hpa, longitude, latitude, $
    title='GEMS O3P Under 300 hPa ' + datetime_str, $
    ;range=[20, 100], $
    range=[0, 60], $
    pngfile=outputpath + '04_under300hPa.png';, $

  nanidx = where(o3under500hpa lt -1.0e20, /null)
  o3under500hpa[nanidx] = !values.f_nan
  plot_sat_proj, o3under500hpa, longitude, latitude, $
    title='GEMS O3P Under 500 hPa ' + datetime_str, $
    range=[20, 100], $
    pngfile=outputpath + '05_under500hPa.png';, $
    ;/scp_send, $
    ;scp_dest=scp_dest

  plot_sat_proj, trppsprs, longitude, latitude, $
    title='GEMS O3P TropopausePressure ' + datetime_str, $
    range=[50, 300], $
    pngfile=outputpath + '10_tropopausepressure.png';, $

  ;if not file_test(outputpath) then begin
    ;file_mkdir, outputpath
  ;endif

  ; plotting layers
  ;for i = 0, 23 do begin
    ;plot_sat_proj, reform(o3[*, *, i]), longitude, latitude, $
      ;title='GEMS O3P Layer' + string(i, format='(i02)') + ' ' + datetime_str, $
      ;range=[0, 40], $
      ;pngfile=outputpath + 'gems_l2_o3p_' + datetime_str + '_layer' + string(i, format='(i02)') + '_wl310340_tb_fnl.png';, $
      ;;/scp_send
  ;endfor

endif else begin
  print, 'Can not find any file.'
endelse


end
