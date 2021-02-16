outputpath = './plot/'

path = '/data2/L2_GEMS/val_1008/'
path = '/data1/L2_GEMS/group/o3p/4x4/nopolc/v6/'

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

searchstr = 'GEMS_O3P_20200616_*'
filelist = file_search(path + searchstr)
nfile = n_elements(filelist)
fn = outputpath + 'gems_l2_o3p_20200616_hourly_under10km_tropo.gif'

closeflag = 0
for inum = 0, nfile-1 do begin
  o3p = ds_read_gems_l2_o3p(filelist[inum], varlist=varlist)
  mondate = strmid(filelist[inum], 12+18, 4, /reverse)
  utchhmm = strmid(filelist[inum], 7+18, 4, /reverse)
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

  outfile = outputpath + 'gems_l2_o3p_' + mondate + 'T' + utchhmm + $
      '_hourly_under10km_tropo.png'
  plot_gems_satellite, o3under10km, longitude, latitude, $
    ;title='GEMS L2 O3P for ' + mondate+'T'+utchhmm, $
    title=' ', $
    range=[20, 60], $
    outfile=outfile, $
    /scp_send
  
  cao3 = o3p.ColumnAmountO3
  if mondate eq '0806' then begin
    rangetotal = [250, 310]
  endif else begin
    rangetotal = [250, 345]
  ENDELSE

  ;plot_gems_satellite, reform(cao3[0, *, *]), longitude, latitude, $
    ;title='GEMS L2 O3P Total ozone for ' + mondate + 'T' + utchhmm + 'UTC', $
    ;range=rangetotal, $
    ;pngfile=outputpath + 'gems_l2_o3p_2020' + mondate + utchhmm +'_totalozone_wl310340_me0.5.png', $
    ;/scp_send
  ;plot_gems_satellite, reform(cao3[1, *, *]), longitude, latitude, $
    ;title='GEMS L2 O3P Tropospheric ozone for ' + mondate + 'T' + utchhmm + 'UTC', $
    ;range=[20, 100], $
    ;pngfile=outputpath + 'gems_l2_o3p_2020' + mondate + utchhmm + '_tropocolumn_wl310340_me0.5.png', $
    ;/scp_send
  ;plot_gems_satellite, reform(cao3[2, *, *]), longitude, latitude, $
    ;title='GEMS L2 O3P Stratospheric ozone for ' + mondate + 'T' + utchhmm + 'UTC', $
    ;range=[100, 200], $
    ;pngfile=outputpath + 'gems_l2_o3p_2020' + mondate + utchhmm + '_stratocolumn_wl310340_me0.5.png', $
    ;/scp_send
endfor

spawn, 'convert -delay 200 ' + $
  './plot/gems_l2_o3p_0616T*_hourly_under10km_tropo.png ' + $
  '-loop 0 '+ fn

if not keyword_Set(scp_dest1) then begin
  scp_dest = 'soodal@164.125.38.179:/home/soodal/works/plot'
endif else begin
  scp_dest = scp_dest1
endelse

; send image to pc
  spawn, 'scp -P18742 -p ' + fn + $
  ' ' + scp_dest

end
