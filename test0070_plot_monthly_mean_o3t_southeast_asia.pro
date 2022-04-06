

outputpath = './plot/monthly_mean_ecf02/'


fn_03utc = '~/data/monthly_gems_o3t/gems_o3t_toz_202010_03utc.sav'
restore, fn_03utc
monthly_o3_03utc = monthly_o3
sz = size(monthly_o3_03utc, /dimension)
zonal_anomaly_03utc = monthly_o3_03utc[0:sz[0]-1, 0:sz[1]-1] - $
  transpose(rebin(mean(monthly_o3_03utc[0:sz[0]-1, 0:sz[1]-1], dimension=1, /nan), sz[1], sz[0]))

for ih = 0, 23 do begin
  yyyymm = '202010'
  utc = string(ih, format='(i02)')

  fn = '~/data/monthly_gems_o3t_ecf02/gems_o3t_toz_' + yyyymm + '_' + utc + 'utc.sav'
  if file_test(fn) then begin
    restore, fn
    ;==========================================================================
    ; O3P
    ;==========================================================================
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

    ;==========================================================================
    ; O3T
    ;==========================================================================
    ; Satellite projection
    sz = size(monthly_o3, /dimension)
    zonal_anomaly = monthly_o3[0:sz[0]-1, 0:sz[1]-1] - $
      transpose(rebin(mean(monthly_o3[0:sz[0]-1, 0:sz[1]-1], dimension=1, /nan), sz[1], sz[0]))
    ;plot_south_east_asia, monthly_o3[0:sz[0]-1, 0:sz[1]-1], $
      ;monthly_lon[0:sz[0]-1, 0:sz[1]-1], $
      ;monthly_lat[0:sz[0]-1, 0:sz[1]-1], $
      ;title='GEMS L2 O3T Monthly Mean 2020 10 ' + utc + 'UTC', $
      ;cb_title='O3[DU]', $
      ;range=[25, 50], $
      ;colortable=22, $
      ;;/ctreverse, $
      ;pngfile=outputpath + 'gems_l2_o3t_' + yyyymm + '_' + utc + 'utc_southeast_asia.png';, $

    range = [5, 10, 15, 20]
    for irange=0, n_elements(range)-1 do begin 
      range_str = strtrim(string(range[irange]), 2)
      plot_south_east_asia, zonal_anomaly, $
        monthly_lon[0:sz[0]-1, 0:sz[1]-1], $
        monthly_lat[0:sz[0]-1, 0:sz[1]-1], $
        title='GEMS L2 O3T Monthly Mean Zonal Anomaly 2020 10 ' + utc + 'UTC', $
        cb_title='O3[DU]', $
        range=[-1*range[irange], range[irange]], $
        colortable=70, $
        /ctreverse, $
        pngfile=outputpath + range_str + '/' + 'gems_l2_o3t_' + yyyymm + '_' + utc + $
        'utc_zonal_anomaly_southeast_asia_-' + range_str + '_' + range_str + '.png';, $
    endfor

    plot_south_east_asia, zonal_anomaly-zonal_anomaly_03utc, $
      monthly_lon[0:sz[0]-1, 0:sz[1]-1], $
      monthly_lat[0:sz[0]-1, 0:sz[1]-1], $
      title='GEMS L2 O3T Monthly Mean Zonal Anomaly from 03UTC 2020 10 ' + utc + 'UTC', $
      cb_title='O3[DU]', $
      range=[-3, 3], $
      colortable=70, $
      /ctreverse, $
      pngfile=outputpath + 'gems_l2_o3t_' + yyyymm + '_' + utc + $
      'utc_zonal_anomaly_from_03utc_southeast_asia_-3_3.png';, $
    ; South easth asia region
    ;plot_sat_proj, monthly_o3, monthly_lon[0:sz[0]-1, 0:sz[1]-1], monthly_lat[0:sz[0]-1, 0:sz[1]-1], $
            ;title='GEMS L2 O3T Monthly Mean 2021 06 ' + utc + 'UTC', $
            ;cb_title='O3[DU]', $
            ;range=[260, 280], $
            ;colortable=22, $
            ;;/ctreverse, $
            ;pngfile=outputpath + 'gems_l2_o3t_' + yyyymm + '_' + utc + 'utc.png';, $

    ;zonal_anomaly = monthly_o3[0:sz[0]-1, 0:sz[1]-1] - $
      ;transpose(rebin(mean(monthly_o3[0:sz[0]-1, 0:sz[1]-1], dimension=1, /nan), sz[1], sz[0]))
    ;plot_sat_proj, zonal_anomaly, monthly_lon[0:sz[0]-1, 0:sz[1]-1], monthly_lat[0:sz[0]-1, 0:sz[1]-1], $
            ;title='GEMS L2 O3T Monthly Mean 2021 06 ' + utc + 'UTC', $
            ;cb_title='O3[DU]', $
            ;range=[-5, 5], $
            ;colortable=70, $
            ;/ctreverse, $
            ;pngfile=outputpath + 'gems_l2_o3t_' + yyyymm + '_' + utc + 'utc_zonal_anomaly.png';, $

  endif
endfor
end
