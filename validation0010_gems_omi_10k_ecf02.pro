
total_omivals = []
total_gemsvals = []
for imon = 4, 6 do begin
  smon = string(imon, format='(I02)')
  ;gemso3p_outfilelist = file_search($
    ;'/data2/L2_GEMS/val_1008/GK2_GEMS_O3P_2020' + smon + '??_????_v1.ba.nc4')

  hours = timegen(start=julday(imon, 1, 2021, 0, 45), $
    final=julday(imon+1, 0, 2021, 0, 45), $
    units='H')

  ;gemsyyyymmdd = strmid(gemso3p_outfilelist, 16+6, 8, /reverse)
  ;gemshhmi = strmid(gemso3p_outfilelist, 7+6, 4, /reverse)
  ;gemsyyyymmdd_hhmi = strmid(gemso3p_outfilelist, 16+6, 13, /reverse)
  caldat, hours, mon, day, year, hour, minute, second
  yyyy = string(year, format='(i04)')
  mm = string(mon, format='(i02)')
  dd = string(day, format='(i02)')
  hh = string(hour, format='(i02)')
  mi = string(minute, format='(i02)')
  ss = string(second, format='(i02)')

  monthly_omivals = []
  monthly_gemsvals = []
  flag_dailyplot = 0
  ;for ifile=0, 3 do begin

  for i=0, n_elements(hours)-1 do begin
  ;for ifile=0, n_elements(gemso3p_outfilelist)-1 do begin
    ;year = fix(strmid(gemsyyyymmdd[ifile], 0, 4))
    ;mon = fix(strmid(gemsyyyymmdd[ifile], 4, 2))
    ;day = fix(strmid(gemsyyyymmdd[ifile], 6, 2)) ; utc
    ;hour = fix(strmid(gemshhmi[ifile], 0, 2)) ; utc
    ;mi = fix(strmid(gemshhmi[ifile], 2, 2)) ; utc

    gemso3p_hourfilelist = file_search($
      '/data2/L2_GEMS.nier/monthly_val/O3P/GK2_GEMS_L2_' + yyyy $
      + mm[i] + dd[i]+'_' +hh[i]+ mi[i]+'_O3P_*.nc', $
      count=nfile)
    print, gemso3p_hourfilelist
    
    if nfile ne 0 then begin
      print, gemso3p_hourfilelist
      omivals = []
      gemsvals = []
      collocate_gemsl2o3p_omil2profoz_on_gems_grid, $
        year[i], mon[i], day[i],$ hour[i], minute[i], $
        colloc, $
        height=10., $; omi profile level altitud dimension = [km]
        gemsfile=gemsfile


      if n_elements(gemsvals) ge 2 and n_elements(omivals) ge 2 then begin
        if flag_dailyplot then begin
          plot_gems_validation, gemsvals, omivals, $
            filename='./plot/gems_l2_o3p_val_with_omi_' + $
              yyyy[i] + mm[i] + dd[i] + $
              'T' + hh[i] + mi[i] + '_under10km_ecf0.2_30min.png', $
              range=[0, 60], $
              delta=2.0
        endif
        monthly_omivals = [monthly_omivals, omivals]
        monthly_gemsvals = [monthly_gemsvals, gemsvals]
        total_omivals = [total_omivals, omivals]
        total_gemsvals = [total_gemsvals, gemsvals]
      endif
    endif
  endfor

  if n_elements(monthly_gemsvals) gt 1 and $
    n_elements(monthly_omivals) gt 1 then begin
  ; plot tropospheric column ozone for tropics lat < 10
    plot_gems_validation, monthly_gemsvals, monthly_omivals, $
      filename='./plot/gems_l2_o3p_val_with_omi_' + $
        yyyy[0] + mm[0] + '_under10km_ecf0.2_30min.png', $
      xtitle='GEMS Tropospheric O3', $
      ytitle='OMI Tropospheric O3', $
      cblim=[0, 60], $
      range=[0, 60], $
      delta=2.0
  endif
endfor

  ; plot tropospheric column ozone 
if n_elements(total_gemsvals) gt 1 and $
  n_elements(total_omivals) gt 1 then begin
  plot_gems_validation, total_gemsvals, total_omivals, $
    filename='./plot/gems_l2_o3p_val_with_omi_' + $
      '202008-202010_under10k_ecf0.2_30min.png', $
    xtitle='GEMS Tropospheric O3', $
    ytitle='OMI Tropospheric O3', $
    cblim=[0, 60], $
    range=[0, 60], $
    delta=2.0
endif
end
