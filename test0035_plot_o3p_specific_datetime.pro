outputpath = './plot/'

path = '/data2/L2_GEMS/val_1008/'

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


datetime = '1006'
searchstr = 'GK2_GEMS_O3P_2020' + datetime + '_0345.nc4'
filelist = file_search(path + searchstr)
nfile = n_elements(filelist)

for inum = 0, nfile-1 do begin
  o3p = ds_read_gems_l2_o3p(filelist[inum], varlist=varlist)
  mondate = strmid(filelist[inum], 12, 4, /reverse)
  utchhmm = strmid(filelist[inum], 7, 4, /reverse)
  ecf = o3p.EffectiveCloudFractionUV

  ecf02idx = where(ecf gt 0.2, /null)
  ecfnanidx = where(ecf lt 0, /null)

  ds_gems_l2o3p_accum, o3p, o3under10km, height=10.

  if inum eq 0 then begin
    longitude = o3p.longitude
    latitude = o3p.latitude
  ENDIf


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
    title='GEMS L2 O3P for ' + mondate, $
    range=[20, 60], $
    pngfile=outputpath + 'gems_l2_o3p_date_' + mondate + 'T' + utchhmm + '_under10km_tropo_wl310340_me0.5.png', $
    /scp_send

  cao3 = o3p.ColumnAmountO3
  if mondate eq '0806' then begin
    rangetotal = [250, 310]
  endif else begin
    rangetotal = [250, 345]
  ENDELSE

  plot_sat_proj, reform(cao3[0, *, *]), longitude, latitude, $
    title='GEMS L2 O3P Total ozone for ' + mondate + 'T' + utchhmm + 'UTC', $
    range=rangetotal, $
    pngfile=outputpath + 'gems_l2_o3p_2020' + mondate + utchhmm +'_totalozone_wl310340_me0.5.png', $
    /scp_send

  ;plot_sat_proj, reform(cao3[1, *, *]), longitude, latitude, $
    ;title='GEMS L2 O3P Tropospheric ozone for ' + mondate + 'T' + utchhmm + 'UTC', $
    ;range=[20, 100], $
    ;pngfile=outputpath + 'gems_l2_o3p_2020' + mondate + utchhmm + '_tropocolumn_wl310340_me0.5.png', $
    ;/scp_send

  ;plot_sat_proj, reform(cao3[2, *, *]), longitude, latitude, $
    ;title='GEMS L2 O3P Stratospheric ozone for ' + mondate + 'T' + utchhmm + 'UTC', $
    ;range=[100, 200], $
    ;pngfile=outputpath + 'gems_l2_o3p_2020' + mondate + utchhmm + '_stratocolumn_wl310340_me0.5.png', $
    ;/scp_send
endfor


end
