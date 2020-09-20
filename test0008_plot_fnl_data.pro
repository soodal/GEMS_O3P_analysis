
; for 2020
filename = '/data1/app/gemsl2_2020.v0.2.i/data/o3p/ATMOS/fnl13.75LST/fnltemp/fnltemp_20200616.dat'
fnltemp = ds_read_fnl_dat(filename, 360, 180, 26)

x=findgen(360, start=-180)
y=findgen(180)-90

xx = rebin(x, [360, 180])
yy = rebin(transpose(y), [360, 180])

limit=[-90, -180, 90, 180]
;limit=[-90, 0, 90, 360] ; my fnl
ds_plot_geographic, fnltemp[*, *, 25], xx, yy, limit=limit, $
  title='FNL Daily', $
  cbtitle='hPa', $
  range=[250, 350], $
  pngfile='./plot/fnltemp_20200616.png'


; for monthly
filename = '/data1/app/gemsl2_2020.v0.2.i/data/o3p/ATMOS/fnl13.75LST/fnltemp/fnltempavg06.dat'
fnlavgtemp = ds_read_fnl_dat(filename, 360, 180, 26)

x=findgen(360, start=-180)
y=findgen(180)-90

xx = rebin(x, [360, 180])
yy = rebin(transpose(y), [360, 180])

limit=[-90, -180, 90, 180]
;limit=[-90, 0, 90, 360] ; my fnl
ds_plot_geographic, fnlavgtemp[*, *, 25], xx, yy, limit=limit, $
  title='FNL Monthly', $
  cbtitle='hPa', $
  range=[250, 350], $
  pngfile='./plot/fnltempavg06.png'

ds_plot_geographic, fnltemp[*, *, 25] - fnlavgtemp[*, *, 25], xx, yy, limit=limit, $
  title='FNL daily - Monthly', $
  cbtitle='hPa', $
  range=[-20, 20], $
  pngfile='./plot/tropopause_diff_monthly_daily.png'
end
