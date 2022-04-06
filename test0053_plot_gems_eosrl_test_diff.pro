; plotting every scene aug to oct for gems l2 o3p
;

scp_dest = 'soodal@164.125.38.179:/home/soodal/works/plot/paper/2021_initial_report'

;initial variable setting
outputpath = './plot/sol_test/'

nier_path = '/data2/L2_GEMS/L1C_test/o3p/4x4/NIER/'
eosrl_path = '/data2/L2_GEMS/L1C_test/o3p/4x4/EOSRL/'

varlist = ['EffectiveCloudFractionUV', 'ProcessingQualityFlags', $
           'O3' , $
           ;'O3Apriori', 'O3AprioriError', $
           ;'CloudPressure', $
           ;'SimulatedRadiances', $
           ;'O3Apriori', 'O3AprioriError',$
           'ColumnAmountO3', $
           'Latitude' ,'Longitude', $
           'Time','Altitude' ,    $
           'Pressure', $
           ;'TropopausePressure', $
           ;'Wavelengths', $
           'WavelengthsWholeRange']

; run

juldays = timegen(start=julday(3, 29, 2021, 0, 45), final=julday(5, 2, 2021, 23, 45)

sub = ['EOSRL']
for i=0, n_elements(juldays)-1 do begin
  jd = juldays[i]
  caldat, jd, month, day, year, hour, minute
  yyyy = string(year, format='(i04)')
  mm = string(month, format='(i02)')
  dd = string(day, format='(i02)')
  hh = string(hour, format='(i02)')
  mi = string(minute, format='(i02)')
  print, yyyy, mm, dd, hh, mi
  datetime_str = yyyy + mm + dd+ '_' +  hh+ mi

  search_str1 = 'GK2_GEMS_L2_O3P_' + datetime_str + '_' + sub[0] + '_BIN4x4.nc'
  filelist1 = file_search(path + search_str1)
  nfile1 = n_elements(filelist1)

  search_str2 = 'GK2_GEMS_L2_O3P_' + datetime_str + '_' + sub[1] + '_BIN4x4.nc'
  filelist2 = file_search(path + search_str2)
  nfile2 = n_elements(filelist2)

  if total(strlen(filelist1)) ne 0 then begin 
    o3p_avg = ds_read_gems_l2_o3p(filelist1[0], varlist=varlist)
    o3p_odd = ds_read_gems_l2_o3p(filelist2[0], varlist=varlist)
    utchhmm = strmid(filelist1[0], 7, 4, /reverse)

    ecf_avg = o3p_avg.EffectiveCloudFractionUV

    ecf02idx = where(ecf_avg gt 0.2, /null)
    ecfnanidx = where(ecf_avg lt 0, /null)

    ds_gems_l2o3p_accum, o3p_avg, o3under300hpa_avg, hpa=300.
    ds_gems_l2o3p_accum, o3p_odd, o3under300hpa_odd, hpa=300.

    longitude = o3p_avg.longitude
    latitude = o3p_avg.latitude

    nanidx = where(o3under300hpa_avg lt 0 and o3under300hpa_odd lt 0, /null)

    o3under300hpa_avg[ecf02idx] = !values.f_nan
    o3under300hpa_avg[ecfnanidx] = !values.f_nan
    o3under300hpa_avg[nanidx] = !values.f_nan

    nanidx = where(longitude lt -500, /null)
    longitude[nanidx] = !values.f_nan
    latitude[nanidx] = !values.f_nan

    nanidx = where(latitude lt -500, /null)
    longitude[nanidx] = !values.f_nan
    latitude[nanidx] = !values.f_nan

    plot_sat_proj, o3under300hpa_avg - o3under300hpa_odd, longitude, latitude, $
      title='GEMS L2 O3P for ' + datetime_str + 'UTC', $
      range=[-0.5, 0.5], $
      pngfile=outputpath + 'gems_l2_o3p_' + datetime_str + '_under300hpa_sol_avg_odd_diff.png', $
      /scp_send, $
      scp_dest=scp_dest


    cao3_avg = o3p_avg.ColumnAmountO3
    cao3_odd = o3p_odd.ColumnAmountO3
    ;if mondate eq '0806' then begin
      ;rangetotal = [250, 310]
    ;endif else begin
      ;rangetotal = [250, 345]
    ;ENDELSE

    rangetotal = [250, 500]
    sz = size(cao3_avg, /dim)
    if sz[0] ne 3 then begin
      data = cao3_avg[*, *, 0]
    endif else begin
      data = cao3_avg[0, *, *]
    endelse

    ;plot_sat_proj, reform(data), longitude, latitude, $
      ;title='GEMS L2 O3P Total ozone for ' + datetime_str + 'UTC', $
      ;range=rangetotal, $
      ;pngfile=outputpath + 'gems_l2_o3p_date_' + datetime_str +'_toz_' + sub[isub] + '.png', $
      ;/scp_send, $
      ;scp_dest=scp_dest
    ;plot_sat_proj, reform(cao3[2, *, *]), longitude, latitude, $
      ;title='GEMS L2 O3P Tropospheric ozone for ' + mondate + 'T' + utchhmm + 'UTC', $
      ;range=[20, 100], $
      ;pngfile=outputpath + 'gems_l2_o3p_2020' + mondate + utchhmm + '_tropocolumn_wl310340_me0.5.png', $
      ;/scp_send
    ;plot_sat_proj, reform(cao3[1, *, *]), longitude, latitude, $
      ;title='GEMS L2 O3P Stratospheric ozone for ' + mondate + 'T' + utchhmm + 'UTC', $
      ;range=[100, 200], $
      ;pngfile=outputpath + 'gems_l2_o3p_2020' + mondate + utchhmm + '_stratocolumn_wl310340_me0.5.png', $
      ;/scp_send
  ;endfor
  endif else begin
    print, 'Cannot find any file.'
  endelse
endfor
end
