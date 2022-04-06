
restore, '~/data/monthly_tropo_o3/202010_03utc.sav'
outputpath = './plot/monthly_mean/'
idx = where(abs(lat2d) gt 20, /null)
monthly_300hpa_o3[idx] = !values.f_nan
plot_sat_proj, monthly_300hpa_o3, lon2d, lat2d, $
        title='MERRA2 O3 under 300 hPa Monthly Mean', $
        range=[15, 25], $
        colortable=22, $
        pngfile=outputpath + 'merra2_o3_under300hpa_202010_03utc_tropical_0_35.png';, $
END

