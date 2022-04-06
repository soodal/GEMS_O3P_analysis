
total_omivals = []
total_gemsvals = []
for imon = 8, 10 do begin
  smon = string(imon, format='(I02)')
  gemso3p_outfilelist = file_search($
    '/data2/L2_GEMS/val_1008/GK2_GEMS_O3P_2020' + smon + '??_??45.nc4')

  gemsyyyymmdd = strmid(gemso3p_outfilelist, 16, 8, /reverse)
  gemshhmi = strmid(gemso3p_outfilelist, 7, 4, /reverse)
  gemsyyyymmdd_hhmi = strmid(gemso3p_outfilelist, 16, 13, /reverse)

  monthly_omivals = []
  monthly_gemsvals = []
  flag_dailyplot = 0
  for igemsfile=0, 0 do begin ;n_elements(gemso3p_outfilelist)-1 do begin
  ;for igemsfile=0, n_elements(gemso3p_outfilelist)-1 do begin
    year = fix(strmid(gemsyyyymmdd[igemsfile], 0, 4))
    mon = fix(strmid(gemsyyyymmdd[igemsfile], 4, 2))
    day = fix(strmid(gemsyyyymmdd[igemsfile], 6, 2)) ; utc
    hour = fix(strmid(gemshhmi[igemsfile], 0, 2)) ; utc
    mi = fix(strmid(gemshhmi[igemsfile], 2, 2)) ; utc

    omivals = []
    gemsvals = []
    collocate_gemsl2o3p_omil2profoz, year, mon, day, hour, mi, gemsvals, omivals, $
      hpa=300 ; omi profile level altitud dimension = [km]
      ;height=10. ; omi profile level altitud dimension = [km]

    if n_elements(gemsvals) ge 2 and n_elements(omivals) ge 2 then begin
      if flag_dailyplot then begin
        plot_gems_validation, gemsvals, omivals, $
          filename='./plot/gems_l2_o3p_val_with_omi_' + gemsyyyymmdd[igemsfile] + $
            'T' + gemshhmi[igemsfile]+'_under300hpa_ecf1.png', $
            range=[0, 100], $
            delta=2.0
      endif
      monthly_omivals = [monthly_omivals, omivals]
      monthly_gemsvals = [monthly_gemsvals, gemsvals]
      total_omivals = [total_omivals, omivals]
      total_gemsvals = [total_gemsvals, gemsvals]
    endif
  endfor

  if n_elements(monthly_gemsvals) gt 1 and $
    n_elements(monthly_omivals) gt 1 then begin
  ; plot tropospheric column ozone for tropics lat < 10
    plot_gems_validation, monthly_omivals, monthly_gemsvals, $
      filename='./plot/gems_l2_o3p_val_with_omi_' + $
        strmid(gemsyyyymmdd[0], 0, 6) +'_under300hpa_ecf1.png', $
      ytitle='GEMS Tropospheric O3', $
      xtitle='OMI Tropospheric O3', $
      cblim=[0, 60], $
      range=[0, 60], $
      delta=2.0
  endif
endfor

  ; plot tropospheric column ozone for tropics lat < 10
if n_elements(total_gemsvals) gt 1 and $
  n_elements(total_omivals) gt 1 then begin
  plot_gems_validation, total_omivals, total_gemsvals, $
    filename='./plot/gems_l2_o3p_val_with_omi_' + $
      '202008-202010_under300hpa_ecf1.png', $
    ytitle='GEMS Tropospheric O3', $
    xtitle='OMI Tropospheric O3', $
    cblim=[0, 60], $
    range=[0, 60], $
    delta=2.0
endif
end
