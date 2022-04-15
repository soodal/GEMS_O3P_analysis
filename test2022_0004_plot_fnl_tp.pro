;fnlfile = '/OPER/SYSTEMS/ALGRTH/GEMS/L2/data/o3p/ATMOS/fnl13.75LST/fnltp/fnltp/fnltp_20050703.dat'

;openr, 10, fnlfile

;x = fltarr(360)
;tp = fltarr(360, 180)

;i = 0
;while not eof(10) do begin
  ;readf, 10, x, format='(360i3)'
  ;tp[*, i] = x
  ;i = i+1
;endwhile

lon = findgen(360)-180
lon = rebin(lon, [360, 180])

lat = findgen(180)-90
lat = rebin(lat, [180, 360])
lat = transpose(lat)


ds_convert_fnl_to_dat_for_gemso3p, '20210329', '/data/MODEL/FNL/2021/',fnlout, outdat='/OPER/SYSTEMS/ALGRTH/GEMS/L2/data/o3p/ATMOS/fnl13.75LST/fnltp/fnltp/fnltp_20210329.dat'

plot_sat_proj, fnlout.tp, lon, lat, $
  pngfile='./plot/fnl_tp/fnl_20210329.png', $
  range=[50, 350], $
  title='FNL Tropopause Pressure'
end
