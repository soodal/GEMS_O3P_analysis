
o3p_outfilelist = file_search('/data2/L2_GEMS/val_1008/GK2_GEMS_O3P_202008??_0345.nc4')

gemsyyyymmdd = strmid(o3p_outfilelist, 16, 8, /reverse)
gemshhmi = strmid(o3p_outfilelist, 7, 4, /reverse)
gemsyyyymmdd_hhmi = strmid(o3p_outfilelist, 16, 13, /reverse)

total_omivals = []
total_gemsvals = []
flag_dailyplot = 0
for ifile=0, n_elements(o3p_outfilelist)-1 do begin
  year = fix(strmid(gemsyyyymmdd[ifile], 0, 4))
  mon = fix(strmid(gemsyyyymmdd[ifile], 4, 2))
  day = fix(strmid(gemsyyyymmdd[ifile], 6, 2))
  hour = fix(strmid(gemshhmi[ifile], 0, 2))
  mi = fix(strmid(gemshhmi[ifile], 2, 2))

  omivals = []
  gemsvals = []
  collocate_gemsl2o3p_omil2profoz, year, mon, day, hour, mi, gemsvals, omivals, $
    hpa=300

  if n_elements(gemsvals) ge 2 and n_elements(omivals) ge 2 then begin
    if flag_dailyplot then begin
      plot_gems_validation, gemsvals, omivals, $
        filename='./plot/gems_l2_o3p_val_with_omi_' + gemsyyyymmdd[ifile] + $
          'T' + gemshhmi[ifile]+'_under300hpa.png_ecf02', $
          range=[0, 100], $
          delta=2.0
    endif
    total_omivals = [total_omivals, omivals]
    total_gemsvals = [total_gemsvals, gemsvals]
  endif
endfor

; plot tropospheric column ozone
plot_gems_validation, total_gemsvals, total_omivals, $
  filename='./plot/gems_l2_o3p_val_with_omi_' + $
    strmid(gemsyyyymmdd[0], 0, 6) +'_under300hpa_ecf02.png', $
  cblim=[0, 100], $
  range=[0, 100], $
  delta=2.0

end
