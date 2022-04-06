; plotting every scene aug to oct for gems l2 o3p
;

scp_dest = 'soodal@164.125.38.179:/home/soodal/works/plot/paper/2021_initial_report'

;initial variable setting
outputpath = './plot/eosrl_test/'

path = '/data2/L2_GEMS/L1C_test/o3p/4x4/EOSRL/'
path_nier = '/data2/L2_GEMS/L1C_test/o3p/4x4/NIER/'

varlist = ['EffectiveCloudFractionUV', 'ProcessingQualityFlags', $
           'O3' , $
           ;'O3Apriori', 'O3AprioriError', $
           ;'CloudPressure', $
           ;'SimulatedRadiances', $
           ;'O3Apriori', 'O3AprioriError',$
           'ColumnAmountO3', $
           'Latitude' ,'Longitude', $
           'Time','Altitude' ,    $
           'Pressure', $
           ;'TropopausePressure', $
           ;'Wavelengths', $
           'WavelengthsWholeRange']

; run

;juldays = [julday(3, 29, 2021, 0, 45), julday(5, 2, 2021, 23, 45)]
juldays = timegen(start=julday(3, 29, 2021, 0, 45), $
  final=julday(5, 2, 2021, 23, 45), units='Hours')

juldays1 = timegen(start=julday(3, 29, 2021, 0, 45), $
  final=julday(3, 29, 2021, 23, 45), units='Hours')
juldays2 = timegen(start=julday(5, 2, 2021, 0, 45), $
  final=julday(5, 2, 2021, 23, 45), units='Hours')

juldays = [juldays1, juldays2]


total_omivals = []
total_gemsvals = []
;for imon = 8, 10 do begin
monthly_omivals = []
monthly_gemsvals = []
for itime = 0, n_elements(juldays)-1 do begin
  jd = juldays[itime]
  caldat, jd, month, day, year, hour, minute
  yyyy = string(year, format='(i04)')
  mm = string(month, format='(i02)')
  dd = string(day, format='(i02)')
  hh = string(hour, format='(i02)')
  mi = string(minute, format='(i02)')
  ;print, yyyy, mm, dd, hh, mi

  gemso3p_outfilelist = file_search($
    '/data2/L2_GEMS/L1C_test/o3p/4x4/EOSRL/' $
    + 'GK2_GEMS_L2_O3P_' + yyyy + mm + dd + '_' + hh + mi + '_EOSRL_BIN4x4.nc')

  gemsyyyymmdd = strmid(gemso3p_outfilelist, 16, 8, /reverse)
  gemshhmi = strmid(gemso3p_outfilelist, 7, 4, /reverse)
  gemsyyyymmdd_hhmi = strmid(gemso3p_outfilelist, 16, 13, /reverse)

  flag_dailyplot = 1
  print, gemso3p_outfilelist[0]
  if strlen(gemso3p_outfilelist[0]) ne 0 then begin
    ;for igemsfile=0, 0 do begin ;n_elements(gemso3p_outfilelist)-1 do begin
    ;for igemsfile=0, n_elements(gemso3p_outfilelist)-1 do begin

      omivals = []
      gemsvals = []
      collocate_gemsl2o3p_omil2profoz, year, month, day, hour, minute, gemsvals, omivals, $
        hpa=300, $ 
        gemsfile=gemso3p_outfilelist[0]
        ;height=10. ; omi profile level altitud dimension = [km]

      if n_elements(gemsvals) ge 2 and n_elements(omivals) ge 2 then begin
        if flag_dailyplot then begin
          plot_gems_validation, gemsvals, omivals, $
            filename='./plot/gems_l2_o3p_val_with_omi_' + yyyy + mm + dd + $
              'T' + hh + mi + '_under300hpa_ecf1.png', $
              range=[0, 100], $
              delta=2.0
        endif
        monthly_omivals = [monthly_omivals, omivals]
        monthly_gemsvals = [monthly_gemsvals, gemsvals]
        total_omivals = [total_omivals, omivals]
        total_gemsvals = [total_gemsvals, gemsvals]
      endif
    ;endfor
  endif

  ;if n_elements(monthly_gemsvals) gt 1 and $
    ;n_elements(monthly_omivals) gt 1 then begin
  ;; plot tropospheric column ozone for tropics lat < 10
    ;plot_gems_validation, monthly_gemsvals, monthly_omivals, $
      ;filename='./plot/gems_l2_o3p_val_with_omi_' + $
        ;yyyy + mm + dd + '_' + hh + mi + '_under300hpa_ecf1.png', $
      ;xtitle='GEMS Tropospheric O3', $
      ;ytitle='OMI Tropospheric O3', $
      ;cblim=[0, 60], $
      ;range=[0, 60], $
      ;delta=2.0
  ;endif
endfor

  ; plot tropospheric column ozone for tropics lat < 10
if n_elements(total_gemsvals) gt 1 and $
  n_elements(total_omivals) gt 1 then begin
  validx = where(total_gemsvals ge 1 and total_omivals ge 1, /null)

  plot_gems_validation, total_gemsvals[validx], total_omivals[validx], $
    filename='./plot/eosrl_test/gems_l2_o3p_eosrl_test_with_omi_' + $
      '20210329_20210502_under300hpa_ecf1.png', $
    xtitle='GEMS Tropospheric O3 under 300hPa', $
    ytitle='OMI Tropospheric O3 under 300hPa', $
    cblim=[0, 200], $
    range=[0, 60], $
    delta=2.0
endif
end
