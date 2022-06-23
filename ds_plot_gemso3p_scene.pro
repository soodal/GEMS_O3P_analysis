pro ds_plot_gemso3p_scene, gemsfile, outputpath=outputbasepath, project=project, $
  suffix=suffix, troprange=troprange, plotnum=plotnum

if not keyword_set(plotnum) then begin
  plotnum = indgen(20)
endif


if not keyword_set(outputbasepath) then begin
  outputbasepath = './plot/'
endif

if not keyword_set(troprange) then begin
  troprange = [15, 55]
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

if keyword_set(project) then begin
  project_subdir = [project, datetime_str]
  outputpath = filepath('', root_dir=outputbasepath, subdirectory=project_subdir)
endif

if not file_test(outputpath) then begin
  file_mkdir, outputpath
endif

fileok = file_test(gemsfile)

if fileok then begin 
  o3p = ds_read_gems_l2_o3p(gemsfile, varlist=varlist)

  ecf = o3p.EffectiveCloudFractionUV
  ecp = o3p.CloudPressure

  ecf02idx = where(ecf gt 0.2, /null)
  ecfnanidx = where(ecf lt 0, /null)

  ds_gems_l2o3p_accum, o3p, o3under300hpa, hpa=300.
  ds_gems_l2o3p_accum, o3p, o3under500hpa, hpa=500.
  ds_gems_l2o3p_apriori_accum, o3p, o3under300hpa_apriori, hpa=300.

  longitude = o3p.longitude
  latitude = o3p.latitude

  nanidx = where(o3under300hpa lt 0, /null)

  ;o3under300hpa[ecf02idx] = !values.f_nan
  o3under300hpa[ecfnanidx] = !values.f_nan
  o3under300hpa[nanidx] = !values.f_nan
  o3under300hpa_apriori[ecfnanidx] = !values.f_nan
  o3under300hpa_apriori[nanidx] = !values.f_nan

  nanidx = where(longitude lt -500, /null)
  longitude[nanidx] = !values.f_nan
  latitude[nanidx] = !values.f_nan

  nanidx = where(latitude lt -500, /null)
  longitude[nanidx] = !values.f_nan
  latitude[nanidx] = !values.f_nan

  trppsprs = o3p.TropopausePressure
  dfs = o3p.DegreesOfFreedomForSignal
  residualmean = reform(o3p.ResidualsOfFit[*, *, 0])
  numiter = o3p.NumberOfIterations




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

  failedidx = where(cao3 lt -1.0e20, /null)

  cao3[failedidx] = !values.f_nan


  if n_elements(where(1 eq plotnum, /null)) then begin 
    if keyword_set(suffix) then begin
      outfilename = '01_totalozone_' + suffix +'.png'
    endif else begin
      outfilename = '01_totalozone.png'
    endelse
    plot_sat_proj, reform(cao3[*, *, 0]), longitude, latitude, $
      title='GEMS O3P Total ozone ' + datetime_str, $
      colortable=33, $
      range=rangetotal, $
      cb_title='[DU]', $
      pngfile=outputpath + outfilename;, $
      ;/scp_send, $
      ;scp_dest=scp_dest
  endif

  if n_elements(where(3 eq plotnum, /null)) then begin 
    if keyword_set(suffix) then begin
      outfilename = '03_troposphericozone_' + suffix +'.png'
    endif else begin
      outfilename = '03_troposphericozone.png'
    endelse
    plot_sat_proj, reform(cao3[*, *, 2]), longitude, latitude, $
      title='GEMS O3P Tropospheric ozone ' + datetime_str, $
      colortable=33, $
      ;range=[20, 100], $
      range=troprange, $
      cb_title='[DU]', $
      pngfile=outputpath + outfilename;, $
  endif

  if n_elements(where(4 eq plotnum, /null)) then begin 
    o3under300hpa[failedidx] = !values.f_nan
    nanidx = where(o3under300hpa lt -1.0e20, /null)
    o3under300hpa[nanidx] = !values.f_nan
    if keyword_set(suffix) then begin
      outfilename = '04_under300hPa_' + suffix +'.png'
    endif else begin
      outfilename = '04_under300hPa.png'
    endelse
    plot_sat_proj, o3under300hpa, longitude, latitude, $
      title='GEMS O3P Under 300 hPa ' + datetime_str, $
      colortable=33, $
      ;range=[20, 100], $
      range=troprange, $
      cb_title='[DU]', $
      pngfile=outputpath + outfilename;, $
  endif

  if n_elements(where(8 eq plotnum, /null)) then begin 
    o3under300hpa_apriori[failedidx] = !values.f_nan
    nanidx = where(o3under300hpa_apriori lt -1.0e20, /null)
    o3under300hpa_apriori[nanidx] = !values.f_nan
    if keyword_set(suffix) then begin
      outfilename = '08_A_Priori_under300hPa_' + suffix +'.png'
    endif else begin
      outfilename = '08_A_Priori_under300hPa.png'
    endelse
    plot_sat_proj, o3under300hpa_apriori, longitude, latitude, $
      title='GEMS O3P O3 A Priori Under 300 hPa ' + datetime_str, $
      colortable=33, $
      ;range=[20, 100], $
      range=troprange, $
      cb_title='[DU]', $
      pngfile=outputpath + outfilename;, $
  endif

  if n_elements(where(5 eq plotnum, /null)) then begin 
    ecf[failedidx] = !values.f_nan
    nanidx = where(ecf lt -1.0e20, /null)
    ecf[nanidx] = !values.f_nan
    if keyword_set(suffix) then begin
      outfilename = '05_effectivecloudfraction_' + suffix +'.png'
    endif else begin
      outfilename = '05_effectivecloudfraction.png'
    endelse
    plot_sat_proj, ecf, longitude, latitude, $
      title='GEMS_L2_O3P Effective Cloud Fraction', $
      colortable=33, $
      ;range=[20, 100], $
      range=[0.0, 1.0], $
      cb_title='Cloud Fraction', $
      pngfile=outputpath + outfilename;, $
  endif

  ;cp[failedidx] = !values.f_nan
  ;nanidx = where(cp lt -1.0e20, /null)
  ;cp[nanidx] = !values.f_nan
  ;if keyword_set(suffix) then begin
    ;outfilename = '06_cloudpressure_' + suffix +'.png'
  ;endif else begin
    ;outfilename = '06_cloudpressure.png'
  ;endelse
  ;plot_sat_proj, cp, longitude, latitude, $
    ;title='GEMS_L2_O3P Cloud Pressure', $
    ;colortable=22, $
    ;;range=[20, 100], $
    ;range=[1000.0, 0.0], $
    ;pngfile=outputpath + outfilename;, $

  ;cld_ccp[failedidx] = !values.f_nan
  ;nanidx = where(cld_ccp lt -1.0e20, /null)
  ;cld_ccp[nanidx] = !values.f_nan
  ;if keyword_set(suffix) then begin
    ;outfilename = '07_cloudcentroidpressure_' + suffix +'.png'
  ;endif else begin
    ;outfilename = '07_cloudcentroidpressure.png'
  ;endelse
  ;plot_sat_proj, cld_ccp, longitude, latitude, $
    ;title='GEMS_L2_CLD Cloud Centroid Pressure', $
    ;colortable=22, $
    ;;range=[20, 100], $
    ;range=[1000.0, 0.0], $
    ;pngfile=outputpath + outfilename;, $

  ;nanidx = where(o3under500hpa lt -1.0e20, /null)
  ;o3under500hpa[nanidx] = !values.f_nan
  ;plot_sat_proj, o3under500hpa, longitude, latitude, $
    ;title='GEMS O3P Under 500 hPa ' + datetime_str, $
    ;range=[20, 100], $
    ;pngfile=outputpath + '05_under500hPa.png';, $
    ;;/scp_send, $
    ;;scp_dest=scp_dest


  if n_elements(where(10 eq plotnum, /null)) then begin 
    if keyword_set(suffix) then begin
      outfilename = '10_tropopausepressure_' + suffix +'.png'
    endif else begin
      outfilename = '10_tropopausepressure.png'
    endelse
    plot_sat_proj, trppsprs, longitude, latitude, $
      title='GEMS O3P TropopausePressure ' + datetime_str, $
      colortable=33, $
      range=[10, 400], $
      cb_title='[hPa]', $
      pngfile=outputpath + outfilename;, $
  endif

  if n_elements(where(11 eq plotnum, /null)) then begin 
    if keyword_set(suffix) then begin
      outfilename = '11_dfs_' + suffix +'.png'
    endif else begin
      outfilename = '11_dfs.png'
    endelse
    plot_sat_proj, dfs, longitude, latitude, $
      title='GEMS O3P DegreesOfFreedomForSignal ' + datetime_str, $
      colortable=33, $
      range=[1, 3], $
      cb_title='unitless', $
      pngfile=outputpath + outfilename;, $
  endif

  if n_elements(where(12 eq plotnum, /null)) then begin 
    if keyword_set(suffix) then begin
      outfilename = '12_residualsOfFit_mean_' + suffix +'.png'
    endif else begin
      outfilename = '12_ResidualsOfFit_mean.png'
    endelse
    plot_sat_proj, residualmean, longitude, latitude, $
      title='GEMS O3P ResidualsOfFit ' + datetime_str, $
      colortable=33, $
      range=[0.5, 1.0], $
      cb_title='unitless', $
      pngfile=outputpath + outfilename;, $

  endif

  if n_elements(where(17 eq plotnum, /null)) then begin 
    if keyword_set(suffix) then begin
      outfilename = '17_numberofiterations_' + suffix +'.png'
    endif else begin
      outfilename = '17_numberofiterations.png'
    endelse
    plot_sat_proj, numiter, longitude, latitude, $
      title='GEMS O3P NumberOfIterations ' + datetime_str, $
      colortable=33, $
      range=[0, 10], $
      cb_title='unitless', $
      pngfile=outputpath + outfilename;, $
  endif

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
