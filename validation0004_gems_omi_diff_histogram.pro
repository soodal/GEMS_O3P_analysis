
o3p_outfilelist = file_search('/data2/L2_GEMS/val_1008/GK2_GEMS_O3P_202008??_0345.nc4')

gemsyyyymmdd = strmid(o3p_outfilelist, 16, 8, /reverse)
gemshhmi = strmid(o3p_outfilelist, 7, 4, /reverse)
gemsyyyymmdd_hhmi = strmid(o3p_outfilelist, 16, 13, /reverse)

total_omivals = []
total_gemsvals = []

yn_daily = 0

for ifile=0, n_elements(o3p_outfilelist)-1 do begin
  year = fix(strmid(gemsyyyymmdd[ifile], 0, 4))
  mon = fix(strmid(gemsyyyymmdd[ifile], 4, 2))
  day = fix(strmid(gemsyyyymmdd[ifile], 6, 2))
  hour = fix(strmid(gemshhmi[ifile], 0, 2))
  mi = fix(strmid(gemshhmi[ifile], 2, 2))

  omivals = []
  gemsvals = []
  collocate_gemsl2o3p_omo3pr, year, mon, day, hour, mi, gemsvals, omivals, $
    hpa=300

  if n_elements(gemsvals) ge 2 and n_elements(omivals) ge 2 then begin
    if yn_daily then begin
    plot_gems_validation, gemsvals, omivals, $
      filename='./plot/gems_l2_o3p_val_with_omi_' + gemsyyyymmdd[ifile] + $
        'T' + gemshhmi[ifile]+'.png'
    endif
    total_omivals = [total_omivals, omivals]
    total_gemsvals = [total_gemsvals, gemsvals]
  endif
endfor

pdf = histogram(total_gemsvals-total_omivals, locations=xbin)
phist = barplot(xbin, pdf, xtitle='GEMS-OMI', ytitle='Frequency',$
  title='Histogram', $
  xrange=[-40, 40], $
  width=1.0, $
  /buffer)

pngfile = './histogram_gems_omi_diff.png'
phist.save, pngfile
scp_dest = 'soodal@164.125.38.179:/home/soodal/WindowsHome/Pictures/scp'
spawn, 'scp -P18742 -p ' + pngfile + $
  ' ' + scp_dest


end
