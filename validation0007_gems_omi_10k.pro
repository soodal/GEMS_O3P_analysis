
total_omivals = []
total_gemsvals = []
for imon = 10, 10 do begin
  smon = string(imon, format='(I02)')
  gemso3p_outfilelist = file_search($
    '/data2/L2_GEMS/val_1008/GK2_GEMS_O3P_2020' + smon + '??_????_v1.ba.nc4')

  gemsyyyymmdd = strmid(gemso3p_outfilelist, 16+6, 8, /reverse)
  gemshhmi = strmid(gemso3p_outfilelist, 7+6, 4, /reverse)
  gemsyyyymmdd_hhmi = strmid(gemso3p_outfilelist, 16+6, 13, /reverse)

  monthly_omivals = []
  monthly_gemsvals = []
  flag_dailyplot = 0
  ;for ifile=0, 3 do begin
  for ifile=0, n_elements(gemso3p_outfilelist)-1 do begin
    year = fix(strmid(gemsyyyymmdd[ifile], 0, 4))
    mon = fix(strmid(gemsyyyymmdd[ifile], 4, 2))
    day = fix(strmid(gemsyyyymmdd[ifile], 6, 2)) ; utc
    hour = fix(strmid(gemshhmi[ifile], 0, 2)) ; utc
    mi = fix(strmid(gemshhmi[ifile], 2, 2)) ; utc

    omivals = []
    gemsvals = []
    collocate_gemsl2o3p_omil2profoz, year, mon, day, hour, mi, gemsvals, omivals, $
      height=10. ; omi profile level altitud dimension = [km]

    if n_elements(gemsvals) ge 2 and n_elements(omivals) ge 2 then begin
      if flag_dailyplot then begin
        plot_gems_validation, gemsvals, omivals, $
          filename='./plot/gems_l2_o3p_val_with_omi_' + gemsyyyymmdd[ifile] + $
            'T' + gemshhmi[ifile]+'_under10km_ecf1_30min.png', $
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
    plot_gems_validation, monthly_gemsvals, monthly_omivals, $
      filename='./plot/gems_l2_o3p_val_with_omi_' + $
        strmid(gemsyyyymmdd[0], 0, 6) +'_under10km_ecf1_30min.png', $
      xtitle='GEMS Tropospheric O3', $
      ytitle='OMI Tropospheric O3', $
      cblim=[0, 60], $
      range=[0, 60], $
      delta=2.0
  endif
endfor

  ; plot tropospheric column ozone for tropics lat < 10
if n_elements(total_gemsvals) gt 1 and $
  n_elements(total_omivals) gt 1 then begin
  plot_gems_validation, total_gemsvals, total_omivals, $
    filename='./plot/gems_l2_o3p_val_with_omi_' + $
      '202008-202010_under10k_ecf02.png', $
    xtitle='GEMS Tropospheric O3', $
    ytitle='OMI Tropospheric O3', $
    cblim=[0, 60], $
    range=[0, 60], $
    delta=2.0
endif
end
