o3p_fn = 'GK2B_GEMS_L2_O3P_20200616_0345_winliminit310_prec000_climML_4x4_2021-04-02T102342.nc4'
;o3p_fn = 'GK2B_GEMS_L2_O3P_20200616_0345_winliminit310_prec000_climTB_4x4_2021-04-02T222159.nc4'

o3p_path = '/data2/L2_GEMS/val_1008/'

o3p_path = '/data1/gems/o3p/works/GEMS_O3P_Yonsei/out/paper/'

o3p_fn = 'GK2B_GEMS_L2_O3P_20200616_0345_winliminit300_prec200_climTB_4x4_2021-04-07T092422.nc4'
fl = file_search(o3p_path, 'GK2B_GEMS_L2_O3P_20200616_0345_winliminit???_prec000500_climTB*.nc4')

outpath = '/data1/gems/o3p/works/GEMS_O3P_Yonsei/out/paper/'

varlist = ['EffectiveCloudFractionUV', 'ProcessingQualityFlags', $
            'AveragingKernel', $
            'O3' , $
            'O3Apriori', 'O3AprioriError', $
            'CloudPressure', $
            'SimulatedRadiances', 'Latitude' ,'Longitude', $
            'O3RandomNoiseError', $
            'O3SolutionError', $
            'Time','Altitude' ,    $
            'Temperature', $
            'Pressure', 'TropopausePressure', $
            'Wavelengths', $
            'WavelengthsWholeRange']
;basename = file_basename(o3p_filepath)
;nc_dirname = file_dirname(o3p_filepath) + '/'
;print, nc_dirname

pohang_lon = 129.37963
pohang_lat = 36.03259

for i=0, n_elements(fl)-1 do begin
  o3p_fn = file_basename(fl[i])
  o3p = ds_read_gems_l2_o3p(o3p_path + o3p_fn, varlist=varlist)


  idx = collocate_gemsl2o3p_ozonesonde(o3p.longitude, o3p.latitude, pohang_lon, pohang_lat)

  pohang_pixlon = o3p.longitude[idx]
  pohang_pixlat = o3p.latitude[idx]

  ; plot vertical profiles
  ds_get_vertical_plot_data, o3p.o3, o3p.longitude, o3p.latitude, o3p.pressure, $
    idx, zx, xx, yx, zy, xy, yy

  plot_gems_l2o3p_vertical_slice_profile, zx, xx, yx, zy, xy, yy, outpath, o3p_fn, $
    pohang_pixlon, pohang_pixlat, $
    /do_contour, title='GEMS L2 Ozone profile', range=[0, 50.0]

  ; plot ozone apriori
  ds_get_vertical_plot_data, o3p.o3apriori, o3p.longitude, o3p.latitude, o3p.pressure, $
    idx, zx, xx, yx, zy, xy, yy

  plot_gems_l2o3p_vertical_slice_profile, zx, xx, yx, zy, xy, yy, $
    outpath, o3p_fn+'.o3apriori', $
    pohang_pixlon, pohang_pixlat, $
    /do_contour, title='GEMS L2 Ozone a priori', range=[0, 50.0]

  ; plot random noise error
  ds_get_vertical_plot_data, o3p.O3RandomNoiseError, o3p.longitude, o3p.latitude, o3p.pressure, $
    idx, zx, xx, yx, zy, xy, yy

  plot_gems_l2o3p_vertical_slice_profile, zx, xx, yx, zy, xy, yy, $
    outpath, o3p_fn+'.RNE', $
    pohang_pixlon, pohang_pixlat, $
    /do_contour, title='GEMS L2 Ozone profile random noise error'

  ; plot solution noise error
  ds_get_vertical_plot_data, o3p.O3SolutionError, o3p.longitude, o3p.latitude, o3p.pressure, $
    idx, zx, xx, yx, zy, xy, yy

  plot_gems_l2o3p_vertical_slice_profile, zx, xx, yx, zy, xy, yy, $
    outpath, o3p_fn+'.SE', $
    pohang_pixlon, pohang_pixlat, $
    /do_contour, title='GEMS L2 Ozone profile solution error'

  ; plot Residuals of fit
  ds_get_vertical_plot_data, o3p.O3SolutionError, o3p.longitude, o3p.latitude, o3p.pressure, $
    idx, zx, xx, yx, zy, xy, yy

  plot_gems_l2o3p_vertical_slice_profile, zx, xx, yx, zy, xy, yy, $
    outpath, o3p_fn+'.residualsoffit', $
    pohang_pixlon, pohang_pixlat, $
    /do_contour, title='GEMS L2 Ozone profile residuals of fit'
  ; plot spatial distribution


  ds_gems_l2o3p_accum, o3p, o3accum, hpa=300
  plot_sat_proj, o3accum, o3p.longitude, o3p.latitude, $
    title='GEMS L2 O3P Tropospheric Ozone for under 300hPa', $
    range=[20, 80], $
    pngfile='/data1/gems/o3p/works/GEMS_O3P_Yonsei/out/paper/'+o3p_fn+'.under300hpa.png', $
    /scp_send, $
    scp_dest = 'soodal@164.125.38.179:/home/soodal/works/plot/paper'

  ds_gems_get_850hpa_o3, o3p, o3_850hpa_number_density
  ds_plot_gems_l2_o3p_850hpa, o3_850hpa_number_density, $
    o3p.longitude, o3p.latitude, $
    title='GEMS L2 O3P Tropospheric Ozone for 850hPa', $
    outfile= '/data1/gems/o3p/works/GEMS_O3P_Yonsei/out/paper/'+o3p_fn+'.850hpa.png', $
    scppath='soodal@164.125.38.179:/home/soodal/works/plot/paper'


endfor

end
