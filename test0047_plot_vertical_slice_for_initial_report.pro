scp_dest = 'soodal@164.125.38.179:/home/soodal/works/plot/paper/2021_initial_report/ds_radiance_yame_cal'

o3p_path = '/data2/L2_GEMS/nier_L1C/'
o3p_path = '/data1/gems/o3p/works/GEMS_O3P/out/paper/'
png_path = './plot/2021_initial_report/ds_radiance_yame_cal/'

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

pohang_lon = 129.37963
pohang_lat = 36.03259

; run
for imon = 3, 3 do begin
  mon_str = string(imon, format='(i02)')
  for iday = 22, 22 do begin
    day_str = string(iday, format='(i02)')
    for ihour = 00, 23 do begin
      hour_str = string(ihour, format='(i02)')
      datetime_str = '2021' + mon_str + day_str + '_' + hour_str + '45'
      search_str = 'GK2B_GEMS_L2_O3P_' + datetime_str + '_winlim*.nc4'
      filelist = file_search(o3p_path + search_str)
      nfile = n_elements(filelist)

      if total(strlen(filelist)) ne 0 then begin 
        for inum = 0, nfile-1 do begin ; gems l2o3p file

          o3p_fn = file_basename(filelist[inum])

          o3p = ds_read_gems_l2_o3p(o3p_path + o3p_fn, varlist=varlist)

          idx = collocate_gemsl2o3p_ozonesonde(o3p.longitude, o3p.latitude, pohang_lon, pohang_lat)

          pohang_pixlon = o3p.longitude[idx]
          pohang_pixlat = o3p.latitude[idx]

          ; plot vertical profiles
          ds_get_vertical_plot_data, o3p.o3, o3p.longitude, o3p.latitude, o3p.pressure, $
            idx, zx, xx, yx, zy, xy, yy

          plot_gems_l2o3p_vertical_slice_profile, zx, xx, yx, zy, xy, yy, png_path, o3p_fn, $
            pohang_pixlon, pohang_pixlat, $
            /do_contour, title='GEMS L2 Ozone profile', range=[0, 50.0], $
            /scp_send, $
            scp_dest=scp_dest

          ; plot ozone apriori
          ds_get_vertical_plot_data, o3p.o3apriori, o3p.longitude, o3p.latitude, o3p.pressure, $
            idx, zx, xx, yx, zy, xy, yy

          plot_gems_l2o3p_vertical_slice_profile, zx, xx, yx, zy, xy, yy, $
            png_path, o3p_fn+'.o3apriori', $
            pohang_pixlon, pohang_pixlat, $
            /do_contour, title='GEMS L2 Ozone a priori', range=[0, 50.0], $
            /scp_send, $
            scp_dest=scp_dest

          ; plot random noise error
          ds_get_vertical_plot_data, o3p.O3RandomNoiseError, o3p.longitude, o3p.latitude, o3p.pressure, $
            idx, zx, xx, yx, zy, xy, yy

          plot_gems_l2o3p_vertical_slice_profile, zx, xx, yx, zy, xy, yy, $
            png_path, o3p_fn+'.RNE', $
            pohang_pixlon, pohang_pixlat, $
            /do_contour, title='GEMS L2 Ozone profile random noise error', $
            /scp_send, $
            scp_dest=scp_dest

          ; plot solution noise error
          ds_get_vertical_plot_data, o3p.O3SolutionError, o3p.longitude, o3p.latitude, o3p.pressure, $
            idx, zx, xx, yx, zy, xy, yy

          plot_gems_l2o3p_vertical_slice_profile, zx, xx, yx, zy, xy, yy, $
            png_path, o3p_fn+'.SE', $
            pohang_pixlon, pohang_pixlat, $
            /do_contour, title='GEMS L2 Ozone profile solution error', $
            /scp_send, $
            scp_dest=scp_dest

          ; plot Residuals of fit
          ds_get_vertical_plot_data, o3p.O3SolutionError, o3p.longitude, o3p.latitude, o3p.pressure, $
            idx, zx, xx, yx, zy, xy, yy

          plot_gems_l2o3p_vertical_slice_profile, zx, xx, yx, zy, xy, yy, $
            png_path, o3p_fn+'.residualsoffit', $
            pohang_pixlon, pohang_pixlat, $
            /do_contour, title='GEMS L2 Ozone profile residuals of fit', $
            /scp_send, $
            scp_dest=scp_dest
          ; plot spatial distribution


          ds_gems_l2o3p_accum, o3p, o3accum, hpa=300
          plot_sat_proj, o3accum, o3p.longitude, o3p.latitude, $
            title='GEMS L2 O3P Tropospheric Ozone for under 300hPa', $
            range=[20, 80], $
            pngfile=png_path+o3p_fn+'.under300hpa.png', $
            /scp_send, $
            scp_dest = scp_dest

          ds_gems_get_850hpa_o3, o3p, o3_850hpa_number_density
          ds_plot_gems_l2_o3p_850hpa, o3_850hpa_number_density, $
            o3p.longitude, o3p.latitude, $
            title='GEMS L2 O3P Tropospheric Ozone for 850hPa', $
            outfile= png_path+o3p_fn+'.850hpa.png', $
            scppath= scp_dest

        endfor
      endif else begin
        print, 'Can not find any file.'
      endelse
    endfor
  endfor
endfor

end
