; plotting every scene aug to oct for gems l2 o3p
;

scp_dest = 'soodal@164.125.38.179:/home/soodal/works/plot/paper/2021_initial_report'

;initial variable setting
outputpath = './plot/2021_initial_report/'

path = '/data2/L2_GEMS/nier_L1C/'

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
for imon = 3, 3 do begin
  mon_str = string(imon, format='(i02)')
  for iday = 22, 22 do begin
    day_str = string(iday, format='(i02)')
    for ihour = 00, 23 do begin
      hour_str = string(ihour, format='(i02)')
      datetime_str = '2021' + mon_str + day_str + '_' + hour_str + '45_4x4'
      search_str = 'GK2B_GEMS_L2_O3P_' + datetime_str + '.nc4'
      filelist = file_search(path + search_str)
      nfile = n_elements(filelist)

      if total(strlen(filelist)) ne 0 then begin 
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

          plot_sat_proj, o3under10km, longitude, latitude, $
            title='GEMS L2 O3P for ' + datetime_str + 'UTC', $
            range=[20, 100], $
            pngfile=outputpath + 'gems_l2_o3p_date_' + datetime_str + '_under10km_tropo_wl310340_me0.5.png', $
            /scp_send, $
            scp_dest=scp_dest


          cao3 = o3p.ColumnAmountO3
          ;if mondate eq '0806' then begin
            ;rangetotal = [250, 310]
          ;endif else begin
            ;rangetotal = [250, 345]
          ;ENDELSE

          rangetotal = [250, 500]

          plot_sat_proj, reform(cao3[0, *, *]), longitude, latitude, $
            title='GEMS L2 O3P Total ozone for ' + datetime_str + 'UTC', $
            range=rangetotal, $
            pngfile=outputpath + 'gems_l2_o3p_date_' + datetime_str +'_totalozone_wl310340_me0.5.png', $
            /scp_send, $
            scp_dest=scp_dest
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
        endfor
      endif else begin
        print, 'Can not find any file.'
      endelse
    endfor
  endfor
endfor


end
