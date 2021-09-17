;afile = '/data2/L2_GEMS/L1C_test/o3p/4x4'

total_sol_gemsvals = []
total_oper_gemsvals = []

juldays = [julday(3, 29, 2021, 0, 45), julday(5, 2, 2021, 23, 45)]

sol_accum_total = []
oper_accum_total = []
flag_dailyplot = 0

sol_ecfs = []
oper_ecfs = []

for ifile=0, n_elements(juldays)-1 do begin

  jd = juldays[ifile]

  caldat, jd, mon, day, year, hour, minute, second
  yyyy = string(year, format='(i04)')
  mm = string(mon, format='(i02)')
  dd = string(day, format='(i02)')
  hh = string(hour, format='(i02)')
  mi = string(minute, format='(i02)')
  ss = string(second, format='(i02)')

  datetime = yyyy + mm + dd + '_' + hh + mi
  print, datetime
  sol_avg_path = '/data2/L2_GEMS/L1C_test/o3p/4x4/SOL/GK2_GEMS_L2_O3P_' + datetime + '_SOL_avg_BIN4x4.nc'
  oper_gems_path = '/data2/L2_GEMS/L1C_test/o3p/4x4/NIER/GK2_GEMS_L2_O3P_' + datetime + '_NIER_BIN4x4.nc'


  sol_gemsval = ds_read_gems_l2_o3p(sol_avg_path)

  ds_gems_l2o3p_accum, sol_gemsval, sol_accum, hpa=500
  sol_ecf = sol_gemsval.EffectiveCloudFractionUV

  oper_gemsval = ds_read_gems_l2_o3p(oper_gems_path)

  ds_gems_l2o3p_accum, oper_gemsval, oper_accum, hpa=500
  oper_ecf = oper_gemsval.EffectiveCloudFractionUV

  sol_ecfs = [sol_ecfs, sol_ecf]
  sol_accum_total = [sol_accum_total, sol_accum]

  oper_ecfs = [oper_ecfs, oper_ecf]
  oper_accum_total = [oper_accum_total, oper_accum]


  if flag_dailyplot then begin
    plot_gems_validation, oper_accum, sol_accum_total, $
      filename='./plot/gems_l2_o3p_val_with_omi_' + $
        yyyy + mm + dd + $
        'T' + hh + mi + '_under10km_ecf0.2_30min.png', $
        range=[0, 60], $
        delta=2.0
  endif

  ;if n_elements(oper_accum_total) gt 1 and $
    ;n_elements(sol_accum_total) gt 1 then begin
  ; plot tropospheric column ozone for tropics lat < 10
    validx = where(oper_ecf gt 0.2 and sol_ecf gt 0.2 and oper_accum gt 5. and sol_accum gt 5., /null)
    plot_gems_validation, oper_accum[validx], sol_accum[validx], $
      filename='./plot/sol_test/sol_avg_gems_l2_o3p_val_with_oper_' + $
        yyyy + mm + dd + 'T' + hh + mi + 'Z_under500hpa_ecf0.2_30min.png', $
      xtitle='GEMS OPER Tropospheric O3', $
      ytitle='GEMS SOL Tropospheric O3', $
      cblim=[0, 200], $
      range=[0, 60], $
      delta=2.0
  ;endif
endfor

  ; plot tropospheric column ozone 
validx = where(oper_ecfs gt 0.2 and sol_ecfs gt 0.2 and oper_accum_total gt 5. and sol_accum_total gt 5., /null)
plot_gems_validation, oper_accum_total[validx], sol_accum_total[validx], $
  filename='./plot/sol_test/sol_avg_' + $
    '20210329T0045-20210502T2345_under500hpa_ecf0.2_30min.png', $
  xtitle='GEMS OPER Tropospheric O3', $
  ytitle='GEMS SOL Tropospheric O3', $
  cblim=[0, 200], $
  range=[0, 60], $
  delta=2.0
end
