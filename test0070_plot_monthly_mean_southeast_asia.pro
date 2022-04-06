

outputpath = './plot/monthly_mean/'
for ih = 0, 23 do begin
  yyyymm = '202106'
  utc = string(ih, format='(i02)')

  fn = '~/data/monthly_tropo_o3/gems_o3p_tropo_o3_' + yyyymm + '_' + utc + 'utc.sav'
  if file_test(fn) then begin
    restore, fn
    ;plot_south_east_asia, monthly_tropo_o3[0:89, 0:31], monthly_lon[0:89, 0:31], monthly_lat[0:89, 0:31], $
            ;title='GEMS L2 O3P Tropospheric O3 Monthly Mean', $
            ;range=[20, 45], $
            ;colortable=22, $
            ;pngfile=outputpath + 'gems_l2_o3p_' + yyyymm + '_' + utc + 'utc_tropo.png';, $

    ;plot_south_east_asia, monthly_300hpa_o3[0:89, 0:41], monthly_lon[0:89, 0:41], monthly_lat[0:89, 0:41], $
            ;title='GEMS L2 O3P Under 300 hPa Monthly Mean 2021 06 ' + utc + 'UTC', $
            ;cb_title='O3[DU]', $
            ;range=[0, 60], $
            ;colortable=22, $
            ;pngfile=outputpath + 'gems_l2_o3p_' + yyyymm + '_' + utc + 'utc_300hpa.png';, $

    zonal_anomaly = monthly_300hpa_o3[0:89, 0:41] - $
      transpose(rebin(mean(monthly_300hpa_o3[0:89, 0:41], dimension=1, /nan), 42, 90))
    plot_south_east_asia, zonal_anomaly, monthly_lon[0:89, 0:41], monthly_lat[0:89, 0:41], $
            title='GEMS L2 O3P Under 300 hPa Monthly Mean 2021 06 ' + utc + 'UTC', $
            cb_title='O3[DU]', $
            range=[-10, 10], $
            colortable=70, $
            /ctreverse, $
            pngfile=outputpath + 'gems_l2_o3p_' + yyyymm + '_' + utc + 'utc_300hpa_zonal_anomaly.png';, $
  endif
endfor
end
