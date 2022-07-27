
timelist = timegen(start=julday(3, 1, 2021, 0), final=julday(4, 1, 2021, 0), $
  step_size=6, units='Hours')

for it = 0, n_elements(timelist)-1 do begin
  caldat, timelist[it], month, day, year, hour
  ds_julday2yyyymmddhhmi, timelist[it], yyyy, mm, dd, hh, mi

  fn = '/data/MODEL/FNL/2021/fnl_' + yyyy + mm + dd + '_' + hh + '_00.grib2.nc'
  fnl = ds_read_fnl_nc(fn)

  ds_xymake2d, fnl.lon_0, fnl.lat_0, lon2d, lat2d

  plot_sat_proj, fnl.PRES_P0_L7_GLL0/100., lon2d, lat2d, $
    title='FNL Tropopause Pressure ' + yyyy + mm + dd + ' ' + hh + 'UTC', $
    range=[0, 400], $
    pngfile='./plot/FNL/fnl_tropopause_' + yyyy + mm + dd + '_' + hh + 'UTC.png';, $
  endfor

END
