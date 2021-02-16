outputpath = './plot/'

path = '/data2/L2_GEMS/val_1008/'

varlist = ['EffectiveCloudFractionUV', 'ProcessingQualityFlags', $
           'O3' , $
           'O3Apriori', 'O3AprioriError', $
           'CloudPressure', $
           'SimulatedRadiances', 'Latitude' ,'Longitude', $
           'Time','Altitude' ,    $
           'Pressure', 'TropopausePressure', $
           'Wavelengths', $
           'WavelengthsWholeRange']
for imon = 8, 10 do begin
  smon = string(imon, format='(i02)')
  searchstr = 'GK2_GEMS_O3P_2020' + smon + '??_0345.nc4'
  filelist = file_search(path + searchstr)
  nfile = n_elements(filelist)

  o3p_total = fltarr(174, 512, nfile)
  o3pclear_total = fltarr(174, 512, nfile)

  for inum = 0, nfile-1 do begin
    o3p = ds_read_gems_l2_o3p(filelist[inum], varlist=varlist)
    ecf = o3p.EffectiveCloudFractionUV

    ecf02idx = where(ecf gt 0.2, /null)
    ecfnanidx = where(ecf lt 0, /null)

    ds_gems_l2o3p_accum, o3p, o3under10k, height=10.

    o3under10k_clear = o3under10k

    if inum eq 0 then begin
      longitude = o3p.longitude
      latitude = o3p.latitude
    ENDIf

    nanidx = where(o3under10k lt 0, /null)

    o3under10k_clear[ecf02idx] = !values.f_nan
    o3under10k_clear[ecfnanidx] = !values.f_nan
    o3under10k_clear[nanidx] = !values.f_nan

    ;o3under10k[ecf02idx] = !values.f_nan
    o3under10k[ecfnanidx] = !values.f_nan
    o3under10k[nanidx] = !values.f_nan

    o3p_total[*, *, inum] = o3under10k
    o3pclear_total[*, *, inum] = o3under10k_clear

  endfor

  o3p_mean = mean(o3p_total, dim=3, /nan)
  o3pclear_mean = mean(o3pclear_total, dim=3, /nan)


  plot_sat_proj, o3p_mean, longitude, latitude, $
    title='GEMS L2 O3P Monthly mean for ' + smon, $
    range=[20, 60], $
    pngfile=outputpath + 'gems_l2_o3p_monthly_' + smon + 'tropo_wl310340_me0.5.png', $
    /scp_send

  plot_sat_proj, o3pclear_mean, longitude, latitude, $
    title='GEMS L2 O3P Monthly mean for ' + smon, $
    range=[20, 60], $
    pngfile=outputpath + 'gems_l2_o3p_monthly_' + smon + 'tropo_wl310340_me0.5_ecf0.2.png', $
    /scp_send
endfor


end
