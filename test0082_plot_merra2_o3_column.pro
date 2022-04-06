
jd_list = timegen(start=julday(6, 21, 2021, 3, 0), $
  final=julday(6, 21, 2021, 3, 0), units='Hour')

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

merra2fn = '/data/MODEL/MERRA2/tavg3_3d_asm_Nv/' + yyyy + '/' + mm + '/' + $
    'MERRA2_400.tavg3_3d_asm_Nv.' + yyyy + mm + dd + '.nc4'
merra2 = ds_read_merra2_tavg3_3d_asm_nv(merra2fn)

pres = merra2.pl

itime = 1
o3mmr = reform(merra2.o3[*, *, *, itime])
pres = reform(merra2.PL[*, *, *, itime]) ; mid pressure
temp = reform(merra2.T[*, *, *, itime])
alt = reform(merra2.H[*, *, *, itime])

lon = merra2.lon
lat = merra2.lat

sz = size(o3mmr, /dimension)

plevel = fltarr(sz[0], sz[1], sz[2]+1)

plevel[*, *, 1:71] = exp((alog(pres[*, *, 0:-2]) + alog(pres[*, *, 1:-1]))/2.)
plevel[*, *, 72] = exp(alog(pres[*, *, -1]) + (alog(pres[*,*, -1] - alog(pres[*, *, -2])))/2.)
plevel[*, *, 0] = exp(alog(pres[*, *, 0]) - (alog(pres[*,*, 1] - alog(pres[*, *, 0])))/2.)



ds_xymake2d, lon, lat, lon2d, lat2d

;convert mass mxing ratio into vmr
o3vmr = ds_mmr2vmr(o3mmr)

ds_vmr2du_3d, o3vmr, pres/100., temp, alt, o3, increasing_height=0

toz = total(o3, 3)


plot_sat_proj, toz, lon2d, lat2d, $
  title='MERRA2 O3 Total Column', $
  cb_title='[DU]', $
  range=[220, 430], $
  slat=-40, $
  nlat=60, $
  llon=60, $
  rlon=180, $
  pngfile='./plot/merra2/merra2_o3total_' + yyyymmdd + '.png';, /scp_send

plot_sat_proj, toz, lon2d, lat2d, $
  title='MERRA2 O3 Total Column', $
  cb_title='[DU]', $
  range=[220, 430], $
  ;slat=-40, $
  slat=-5, $
  ;nlat=60, $
  nlat=45, $
  ;llon=60, $
  llon=70, $
  ;rlon=180, $
  rlon=132, $
  pngfile='./plot/merra2/merra2_o3total_gemsroi_' + yyyymmdd + '.png';, /scp_send


;accum = total(o3[*, *, 30:-1], 3)

;plot_sat_proj, accum, lon2d, lat2d, $
  ;title='MERRA2 O3 SFC-300hPa ' + datetime_str, $
  ;cb_title='[DU]', $
  ;range=[0, 60], $
  ;slat=-40, $
  ;nlat=60, $
  ;llon=60, $
  ;rlon=180, $
  ;pngfile='./plot/merra2/merra2_o3_sfc-300hPa_20210329.png'

hpa = 300

merra2o3accum = fltarr(sz[0], sz[1])
  for iy=0, sz[1]-1 do begin
  for ix=0, sz[0]-1 do begin

    for ilevel=sz[2]-1, 0, -1 do begin
      ilayer = ilevel - 1

      if plevel[ix, iy, ilevel] gt hpa*100. $
          and plevel[ix, iy,ilevel-1] gt hpa*100. then begin
        merra2o3accum[ix, iy] = merra2o3accum[ix, iy] + o3[ix, iy, ilayer]

      endif else if plevel[ix, iy, ilevel] gt hpa*100. $
          and plevel[ix, iy, ilevel-1] le hpa*100. then begin

        merra2o3accum[ix, iy] = merra2o3accum[ix, iy] + o3[ix, iy, ilayer]*$
          (plevel[ix, iy, ilevel]-hpa*100.)/$
          (plevel[ix, iy, ilevel]-plevel[ix, iy, ilevel-1])
        break
      endif
    endfor
  endfor
  endfor
  merra2o3accum[where(merra2o3accum lt 0)] = !values.f_nan

plot_sat_proj, merra2o3accum, lon2d, lat2d, $
  title='MERRA2 O3 Sfc-300hPa', $
  cb_title='[DU]', $
  range=[0, 60], $
  ;slat=-40, $
  slat=-5, $
  ;nlat=60, $
  nlat=45, $
  ;llon=60, $
  llon=70, $
  ;rlon=180, $
  rlon=132, $
  pngfile='./plot/merra2/merra2_o3_gemsroi_sfc-300hPa_' + yyyymmdd + '.png'


plot_sat_proj, toz, lon2d, lat2d, $
  title='MERRA2 O3 Total Column', $
  cb_title='[DU]', $
  range=[220, 430], $
  slat=-40, $
  nlat=60, $
  llon=60, $
  rlon=180, $
  pngfile='./plot/merra2/merra2_o3total_' + yyyymmdd + '.png';, /scp_send


end
