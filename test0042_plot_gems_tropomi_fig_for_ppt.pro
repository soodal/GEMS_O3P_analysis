
outputpath = './plot/'

path = '/home/o3p_geun/GEMS/o3p/dat/out/gems_oldout/tropomi_test_20180419t034816/'

searchstr = 'TROPOMI_TEST_L2_O3P_*o02664_mpi8.h5'
filelist = file_search(path + searchstr)
nfile = n_elements(filelist)
fn = outputpath + 'GEMS_O3P_from_TROPOMI_DATA_20180419t034816.png'

closeflag = 0
for inum = 0, nfile-1 do begin
  tropomi = ds_read_tropomi_test_for_gems_o3p(filelist[inum])

  ecf = tropomi.EffectiveCloudFractionUV

  ecf02idx = where(ecf gt 0.2, /null)
  ecfnanidx = where(ecf lt 0, /null)

  ;ds_gems_l2o3p_accum, o3p, o3under10km, height=10.

  o3 = tropomi.o3
  lon = tropomi.longitude
  lat = tropomi.latitude
  alt = tropomi.altitude ; KM
  pres = tropomi.pressure ; hPa

  o3 = transpose(o3, [1, 2, 0])
  alt = transpose(alt, [1, 2, 0])
  pres = transpose(pres, [1, 2, 0])

  if n_elements(o3total) eq 0 then begin
    o3total = o3
    lontotal = lon
    lattotal = lat
    alttotal = alt
    prestotal = pres
  endif else begin
    o3total = [o3total, o3]
    lontotal = [lontotal, lon]
    lattotal = [lattotal, lat]
    alttotal = [alttotal, alt]
    prestotal = [prestotal, pres]
  ENDELSE
  
endfor
;o3total = transpose(o3total, [2, 0, 1])
;alttotal = transpose(alttotal, [2, 0, 1])
;prestotal = transpose(prestotal, [2, 0, 1])

tropomi_allx = create_struct('O3', o3total)
tropomi_allx = create_struct(tropomi_allx, 'Altitude', alttotal)
tropomi_allx = create_struct(tropomi_allx, 'Pressure', prestotal)
nanidx = where(lattotal lt -900 or lon lt -900 or o3total lt -900, /null)
o3total[nanidx] = !values.f_nan

ds_gems_from_tropomi_l2o3p_accum, tropomi_allx, o3accum, height=10

outfile = './plot/gems_from_tropomi_20180419t034816.png'
plot_gems_satellite, o3accum, lontotal, lattotal, $
  ;title='GEMS L2 O3P for ' + mondate+'T'+utchhmm, $
  title=' ', $
  range=[20, 60], $
  outfile=outfile, $
  /scp_send
  


end
