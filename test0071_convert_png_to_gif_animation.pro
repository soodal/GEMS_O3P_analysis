;fn = './plot/monthly_mean/gems_l2_o3p_202106_300hpa.gif'
;spawn, 'convert -delay 100 ' + $
  ;'./plot/monthly_mean/gems_l2_o3p_202106_*utc_300hpa.png ' + $
  ;'-loop 0 '+ fn

;fn = './plot/monthly_mean/gems_l2_o3p_202106_300hpa_zonal_anomaly.gif'
;spawn, 'convert -delay 100 ' + $
  ;'./plot/monthly_mean/gems_l2_o3p_202106_*utc_300hpa_zonal_anomaly.png ' + $
  ;'-loop 0 '+ fn

fn = './plot/monthly_mean/gems_l2_o3t_202010.gif'
spawn, 'convert -delay 100 ' + $
  './plot/monthly_mean/gems_l2_o3p_202010_*utc.png ' + $
  '-loop 0 '+ fn

fn = './plot/monthly_mean/gems_l2_o3t_202010_zonal_anomaly.gif'
spawn, 'convert -delay 100 ' + $
  './plot/monthly_mean/gems_l2_o3t_202010_*utc_zonal_anomaly.png ' + $
  '-loop 0 '+ fn
end
