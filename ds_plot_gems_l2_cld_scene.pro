pro ds_plot_gems_l2_cld_scene, gemsfile, outputpath=outputpath, project=project, sub=sub
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
  gemsl2cld = ds_read_gems_l2_cld(gemsfile, varlist=varlist)

  ecf = gemsl2cld.EffectiveCloudFraction

  ;ecf02idx = where(ecf gt 0.2, /null)
  ecfnanidx = where(ecf lt -1e29, /null)
  ecf[ecfnanidx] = !values.f_nan

  longitude = gemsl2cld.longitude
  latitude = gemsl2cld.latitude


  ccp = gemsl2cld.CloudCentroidPressure
  nanidx = where(ccp lt -1e29, /null)
  ccp[nanidx] = !values.f_nan

  crf = gemsl2cld.CloudRadianceFraction
  nanidx = where(crf lt -1e29, /null)
  crf[nanidx] = !values.f_nan

  tp = gemsl2cld.TerrainPressure
  nanidx = where(tp lt -1e29, /null)
  tp[nanidx] = !values.f_nan

  tr = gemsl2cld.TerrainReflectivity
  nanidx = where(tr lt -1e29, /null)
  tr[nanidx] = !values.f_nan

  sza = gemsl2cld.solarzenithangle
  nanidx = where(sza lt -1e29, /null)
  sza[nanidx] = !values.f_nan

  saa = gemsl2cld.solarazimuthangle
  nanidx = where(saa lt -1e29, /null)
  saa[nanidx] = !values.f_nan

  vza = gemsl2cld.viewingzenithangle
  nanidx = where(vza lt -1e29, /null)
  vza[nanidx] = !values.f_nan

  vaa = gemsl2cld.viewingazimuthangle
  nanidx = where(vaa lt -1e29, /null)
  vaa[nanidx] = !values.f_nan

  ;raa = gemsl2cld.RelativeAzimuthAngle
  ;nanidx = where(raa lt -1e29, /null)
  ;raa[nanidx] = !values.f_nan

  range = [0, 1]

  sz = size(ecf, /dim)

  if not file_test(outputpath) then begin
    file_mkdir, outputpath
  endif

  plot_sat_proj, ecf, longitude, latitude, $
    title='GEMS L2 CLD EffectiveCloudFraction ' + datetime_str, $
    range=range, $
    pngfile=outputpath + '01_ecf.png';, $

  plot_sat_proj, sza, longitude, latitude, $
    title='GEMS L2 CLD SolarZenithAngle ' + datetime_str, $
    range=[0, 90], $
    pngfile=outputpath + '02_sza.png';, $

  plot_sat_proj, vza, longitude, latitude, $
    title='GEMS L2 CLD ViewingZenithAngle ' + datetime_str, $
    range=[0, 90], $
    pngfile=outputpath + '03_vza.png';, $

  ;plot_sat_proj, , longitude, latitude, $
    ;title='GEMS L2 CLD ' + datetime_str, $
    ;range=range, $
    ;pngfile=outputpath + '01_ecf.png';, $


endif else begin
  print, 'Can not find any file.'
endelse


end
