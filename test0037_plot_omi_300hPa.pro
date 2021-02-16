total_omivals = []
total_gemsvals = []
for imon = 8, 10 do begin
  for iday = 1, 31 do begin
    plot_daily_omi, 2020, imon, iday, 4, 45, omivals, height=10., range=[0, 60]
  endfor
endfor

end
