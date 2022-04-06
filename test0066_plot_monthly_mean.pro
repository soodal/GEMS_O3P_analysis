

outputpath = './plot/monthly_mean/'
for ih = 3, 4 do begin
  utc = string(ih, format='(i02)')

  restore, '~/data/monthly_tropo_o3/202010_' + utc + 'utc.sav'
  plot_sat_proj, monthly_tropo_o3[0:89, 0:31], monthly_lon[0:89, 0:31], monthly_lat[0:89, 0:31], $
          title='GEMS L2 O3P Tropospheric O3 Monthly Mean', $
          range=[20, 45], $
          colortable=22, $
          pngfile=outputpath + 'gems_l2_o3p_202010_' + utc + 'utc_tropo.png';, $

  plot_sat_proj, monthly_300hpa_o3[0:89, 0:31], monthly_lon[0:89, 0:31], monthly_lat[0:89, 0:31], $
          title='GEMS L2 O3P Under 300 hPa Monthly Mean', $
          range=[20, 45], $
          colortable=22, $
          pngfile=outputpath + 'gems_l2_o3p_202010_' + utc + 'utc_300hpa.png';, $
endfor
end
