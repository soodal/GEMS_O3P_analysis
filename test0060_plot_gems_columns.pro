; plotting every scene aug to oct for gems l2 o3p
;

scp_dest = 'soodal@164.125.38.179:/home/soodal/works/plot/paper/2021_initial_report'

;initial variable setting


path = '/data/private/soodal/vartest/climaTB_FNL/'

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

jd_list = timegen(start=julday(6, 21, 2021, 3, 45), $
  final=julday(6, 21, 2021, 3, 45), units='Hours')

for idate = 0, n_elements(jd_list)-1 do begin
  caldat, jd_list[idate], mon, day, year, hour, minute
  yyyy = string(year, format='(i04)')
  mm = string(mon, format='(i02)')
  dd = string(day, format='(i02)')
  hh = string(hour, format='(i02)')
  mi = string(minute, format='(i02)')

  datetime_str = yyyy + mm + dd + '_' + hh + mi
  search_str = 'GK2_GEMS_L2_O3P_' + datetime_str + '_TB_BIN4x4.nc'
  print, path+search_str
  filelist = file_search(path + search_str)
  nfile = n_elements(filelist)

  outputpath = './plot/vartest/climaTB_FNL/'
  if total(strlen(filelist)) gt 0 then begin 
    for inum = 0, nfile-1 do begin ; gems l2o3p file
      o3p = ds_read_gems_l2_o3p(filelist[inum], varlist=varlist)
      ;mondate = strmid(filelist[inum], 12, 4, /reverse)
      utchhmm = strmid(filelist[inum], 7, 4, /reverse)

      ecf = o3p.EffectiveCloudFractionUV

      ecf02idx = where(ecf gt 0.2, /null)
      ecfnanidx = where(ecf lt 0, /null)

      ds_gems_l2o3p_accum, o3p, o3under10km, height=10.

      longitude = o3p.longitude
      latitude = o3p.latitude

      nanidx = where(o3under10km lt 0, /null)

      ;o3under10km[ecf02idx] = !values.f_nan
      o3under10km[ecfnanidx] = !values.f_nan
      o3under10km[nanidx] = !values.f_nan

      nanidx = where(longitude lt -500, /null)
      longitude[nanidx] = !values.f_nan
      latitude[nanidx] = !values.f_nan

      nanidx = where(latitude lt -500, /null)
      longitude[nanidx] = !values.f_nan
      latitude[nanidx] = !values.f_nan

      trppsprs = o3p.TropopausePressure

      if not file_test(outputpath) then begin
        file_mkdir, outputpath
      endif

      plot_sat_proj, o3under10km, longitude, latitude, $
        title='GEMS L2 O3P Under 10 KM ' + datetime_str, $
        range=[20, 100], $
        pngfile=outputpath + 'gems_l2_o3p_' + datetime_str + '_under10km_wl310340_tb_fnl.png';, $
        ;/scp_send, $
        ;scp_dest=scp_dest


      cao3 = o3p.ColumnAmountO3
      ;if mondate eq '0806' then begin
        ;rangetotal = [250, 310]
      ;endif else begin
        ;rangetotal = [250, 345]
      ;ENDELSE

      rangetotal = [250, 500]

      sz = size(cao3, /dim)
      if sz[0] eq 3 then begin
        cao3 = transpose(cao3, [1, 2, 0])
      endif

      o3 = o3p.o3
      sz = size(o3, /dim)
      if sz[0] eq 24 then begin
        o3 = transpose(o3, [1, 2, 0])
      endif
        

      plot_sat_proj, reform(cao3[*, *, 0]), longitude, latitude, $
        title='GEMS L2 O3P Total ozone ' + datetime_str, $
        range=rangetotal, $
        pngfile=outputpath + 'gems_l2_o3p_' + datetime_str +'_toz_wl310340_tb_fnl.png';, $
        ;/scp_send, $
        ;scp_dest=scp_dest

      plot_sat_proj, reform(cao3[*, *, 2]), longitude, latitude, $
        title='GEMS L2 O3P Tropospheric ozone ' + datetime_str, $
        range=[20, 100], $
        pngfile=outputpath + 'gems_l2_o3p_' + datetime_str + '_tropocolumn_wl310340_tb_fnl.png';, $

      plot_sat_proj, reform(cao3[*, *, 2]), longitude, latitude, $
        title='GEMS L2 O3P Tropospheric ozone ' + datetime_str, $
        range=[20, 100], $
        pngfile=outputpath + 'gems_l2_o3p_' + datetime_str + '_tropocolumn_wl310340_tb_fnl.png';, $

      plot_sat_proj, trppsprs, longitude, latitude, $
        title='GEMS L2 O3P TropopausePressure ' + datetime_str, $
        range=[20, 300], $
        pngfile=outputpath + 'gems_l2_o3p_' + datetime_str + '_tropopausepressure_fnl.png';, $

      ; plotting layers
      ;for i = 0, 23 do begin
        ;plot_sat_proj, reform(o3[*, *, i]), longitude, latitude, $
          ;title='GEMS L2 O3P Layer' + string(i, format='(i02)') + ' ' + datetime_str, $
          ;range=[0, 40], $
          ;pngfile=outputpath + 'gems_l2_o3p_' + datetime_str + '_layer' + string(i, format='(i02)') + '_wl310340_tb_fnl.png';, $
          ;;/scp_send
      ;endfor
      ;plot_sat_proj, reform(cao3[1, *, *]), longitude, latitude, $
        ;title='GEMS L2 O3P Stratospheric ozone for ' + mondate + 'T' + utchhmm + 'UTC', $
        ;range=[100, 200], $
        ;pngfile=outputpath + 'gems_l2_o3p_2020' + mondate + utchhmm + '_stratocolumn_wl310340_me0.5.png';, $
        ;/scp_send
    endfor
  endif else begin
    print, 'Can not find any file.'
  endelse
endfor


end
