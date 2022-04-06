

jd_list = timegen(start=julday(6, 21, 2021, 3, 45), $
  final=julday(6, 21, 2021, 3, 45), units='Hour')

jd = jd_list[0]

caldat, jd, month, day, year, hour, minute
yyyy = string(year, format='(i04)')
mm = string(month, format='(i02)')
dd = string(day, format='(i02)')
hh = string(hour, format='(i02)')
mi = string(minute, format='(i02)')
;print, yyyy, mm, dd,'_', hh, mi

datetime_str = yyyy + mm + dd+ '_' +  hh+ mi


yyyymmdd = yyyy + mm + dd

fn = '/home/soodal/data/softcal_test/residual/gems_merra2_' + datetime_str + '_raw.nc'

o3p = ds_read_gems_l2_o3p('/home/soodal/data/ln/GEMS/L2O3P/GK2_GEMS_L2_O3P_' + datetime_str + '_4x4.nc')
longitude = o3p.longitude
latitude = o3p.latitude

fitspec_diff = ds_read_nc(fn, 'fitspec_difference_from_merra2') 

diff = reform(fitspec_diff[*, *, 100])
nandix = where(diff lt -1e29, /null)
diff[nandix] = !values.f_nan

nandix = where(diff lt -1, /null)
diff[nandix] = !values.f_nan

rangetotal = minmax(diff, /nan)

plot_sat_proj, diff, longitude, latitude, $
  title='Synthetic Normalized Radiance Difference' + datetime_str, $
  range=[-1, 1], $
  colortable=70, $
  pngfile='./plot/synthetic_radiance_diff/' + '01_fitspec_diff_' + datetime_str + '.png';, $

end
