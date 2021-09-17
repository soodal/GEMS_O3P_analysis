;afile = '/data2/L2_GEMS/L1C_test/o3p/4x4'

total_test_gemsvals = []
total_oper_gemsvals = []

juldays = [julday(3, 29, 2021, 0, 45), julday(5, 2, 2021, 23, 45)]

project = 'BTDF'
sub = 'POLY'
test_gemsvals_total = []
oper_gemsvals_total = []
flag_dailyplot = 0

test_ecfs = []
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
  test_avg_path = '/data2/L2_GEMS/L1C_test/o3p/4x4/' + project + '/GK2_GEMS_L2_O3P_' + datetime + '_' + sub + '_BIN4x4.nc'
  oper_gems_path = '/data2/L2_GEMS/L1C_test/o3p/4x4/NIER/GK2_GEMS_L2_O3P_' + datetime + '_NIER_BIN4x4.nc'


  test_gemsval = ds_read_gems_l2_o3p(test_avg_path)

  ds_gems_l2o3p_accum, test_gemsval, test_accum, hpa=300
  test_ecf = test_gemsval.EffectiveCloudFractionUV

  oper_gemsval = ds_read_gems_l2_o3p(oper_gems_path)

  ds_gems_l2o3p_accum, oper_gemsval, oper_accum, hpa=300
  oper_ecf = oper_gemsval.EffectiveCloudFractionUV

  test_ecfs = [test_ecfs, test_ecf]
  test_gemsvals_total = [test_gemsvals_total, test_accum]

  oper_ecfs = [oper_ecfs, oper_ecf]
  oper_gemsvals_total = [oper_gemsvals_total, oper_accum]


  if flag_dailyplot then begin
    plot_gems_validation, oper_gemsvals_total, test_gemsvals_total, $
      filename='./plot/gems_l2_o3p_val_with_omi_' + $
        yyyy + mm + dd + $
        'T' + hh + mi + '_under10km_ecf0.2_30min.png', $
        range=[0, 60], $
        delta=2.0
  endif

  ;if n_elements(oper_gemsvals_total) gt 1 and $
    ;n_elements(test_gemsvals_total) gt 1 then begin
  ; plot tropospheric column ozone for tropics lat < 10
    ecfidx = where(oper_ecf > 0.2 and test_ecf > 0.2 and oper_accum > 0.1 and test_accum > 0.1, /null)
    plot_gems_validation, oper_accum, test_accum, $
      filename='./plot/' + project + '/' + project + '_' + sub + '_gems_l2_o3p_val_with_oper_' + $
        yyyy + mm + dd + 'T' + hh + mi + 'Z_under500hpa_ecf0.2_30min.png', $
      xtitle='GEMS OPER Tropospheric O3', $
      ytitle='GEMS ' + sub + ' Tropospheric O3', $
      cblim=[0, 60], $
      range=[0, 60], $
      delta=2.0
  ;endif
endfor

  ; plot tropospheric column ozone 
plot_gems_validation, oper_gemsvals_total, test_gemsvals_total, $
  filename='./plot/' + project + '/' + project + '_' + sub + '_test_' + $
    '20210329T0045-20210502T2345_under500hpa_ecf0.2_30min.png', $
  xtitle='GEMS OPER Tropospheric O3', $
  ytitle='GEMS ' + sub + ' Tropospheric O3', $
  cblim=[0, 60], $
  range=[0, 60], $
  delta=2.0
end
