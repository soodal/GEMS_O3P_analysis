
;cities_name = ['Beijing', 'Chungching', 'Shanghai', 'Gwangjou', 'Seoul', 'Busan']
;xidxs = [105, 144, 93, 123, 74, 68]
;yidxs = [38, 129, 111, 189, 57, 76]
cities_name = '112x_447y'
xidxs = [112-1]
yidxs = [447-1]

pos = [0.15, 0.1, 0.8, 0.85]
leg_pos = [0.85, 0.8]




path ='/data1/gems/o3p/works/GEMS_O3P_Yonsei/' 
outputpath = path + 'out/'
projectpath = outputpath + 'snr/' 

filelist = file_search(projectpath + $
  'GK2B_GEMS_L2_O3P_20200616_0345_winliminit3??_prec0200_climML*.nc4')
print, filelist


L1Cmaker = 'EOSRL'
runtime = strmid(filelist, 13+9+8, 10+7, /reverse)
print, runtime

;runtimeidx = sort(runtime)
;recentrunfile = filelist[runtimeidx[-1]]

varlist = ['EffectiveCloudFractionUV', $
  'ProcessingQualityFlags', $
  'AveragingKernel', $
  'O3' , $
  'O3Apriori', 'O3AprioriError','ColumnAmountO3', $
  'CloudPressure', $
  'Latitude' ,'Longitude', $
  'Time','Altitude' ,    $
  'Pressure', 'TropopausePressure', $
  'Wavelengths', $
  'FinalAlgorithmFlags', $
  'WavelengthsWholeRange']

for i=0, n_elements(filelist) - 1 do begin
  fn = file_basename(filelist[i])
  print, fn
  winliminitpos = strpos(fn, 'winliminit')
  fitrange = strmid(fn, winliminitpos+10, 3)+'340'
  fit_range = strmid(fitrange, 0, 3) + '-' + strmid(fitrange, 3, 3)

  precpos = strpos(fn, 'prec')
  prec = strmid(fn, precpos+4, 4)
  ;prec = string(fix(prec), format='(i04)')

  date = strmid(fn, 13, 13)

  ;scp_dest = 'soodal@164.125.38.179:/home/soodal/WindowsHome/Pictures/scp/' $
    ;+ 'ml_clima/' + fitrange + '_' + L1Cmaker 

  limit=[-5, 85, 55, 155]
  nlayer = 24

  zrange1 = [0.15, 0.15, 0.25, 0.5, 0.9, 1.6, 2.9, 5, 8, 12, $
    15, 14, 8, 6, 14, 12, 2, 0, 0, 0, 0, 0, 4, 5]

  zrange2 = [0.3, 0.2, 0.35, 0.7, 1.2, 2.2, 3.8, 7, 11, 20, $
    25, 40, 45, 45, 35, 35, 30, 40, 40, 20, 10, 15, 20, 25]

  data = ds_read_gems_l2_o3p(projectpath + fn, varlist=varlist)


  o3 = data.O3
  o3ap = data.O3Apriori
  pressure = data.Pressure


  cao3 = data.ColumnAmountO3
  o3total = reform(cao3[0, *, *], 174, 512)
  o3strato = reform(cao3[1, *, *], 174, 512)
  o3tropo = reform(cao3[2, *, *], 174, 512)


  _latitude = data.Latitude
  _longitude = data.Longitude
  dim = n_elements(_latitude)

  lat = reform(_latitude, dim)
  lon = reform(_longitude, dim)

  o3total_geolocidx = where(lon gt 0 and lat gt 0 and o3total gt 0, /null)
  o3tropo_geolocidx = where(lon gt 0 and lat gt 0 and o3tropo gt 0, /null)
  o3strato_geolocidx = where(lon gt 0 and lat gt 0 and o3strato gt 0, /null)

  faf = data.finalalgorithmflags
  cp = data.cloudpressure
  ecf = data.effectivecloudfractionuv
  altitude = data.Altitude

  for icity = 0, n_elements(cities_name)-1 do begin

    ;color = ['red', 'blue', 'green', 'cyan'] 
    name = cities_name[icity]

    o3size = size(o3, /dimension)
    if o3size[0] eq 24 then begin
      o3_site = o3[*, xidxs[icity], yidxs[icity]] 
      o3ap_site = o3ap[*, xidxs[icity], yidxs[icity]]
      p = pressure[*, xidxs[icity], yidxs[icity]]
      alt = altitude[*, xidx[icity], yidxs[icity]]

    endif else if o3size[0] eq 174 then begin
      o3_site = o3[xidxs[icity], yidxs[icity], *] 
      o3ap_site = o3ap[xidxs[icity], yidxs[icity], *]
      p = pressure[xidxs[icity], yidxs[icity], *]
      alt = altitude[xidxs[icity], yidxs[icity], *]

    endif


    ;command = "p" + string(ime, format='(i1)') + $
      ;" = plot(o3, p, color=color[ime], overplot=ime, /buffer, " + $
      ;" title='GEMS O3PROFILE 2020-06-16T0345', " + $
      ;" xtitle='O3[DU]', " + $
      ;" ytitle='Pressure[hPa]', " + $
      ;" yrange=[1000, 0], /ylog, pos=pos, name=name[ime])"
    ;print, command
    ;dumm = execute(command)

    p1 = plot(o3_site, p[0:23], color='black', /overplot, /buffer, $
      title='GEMS O3PROFILE 2020-06-16T0345 ' + fitrange + ' ' + L1Cmaker + ' ' + name , $
      xtitle='O3[DU]', $
      ytitle='Pressure[hPa]', $
      xrange=[0, 50], $
      yrange=[1000, 0], /ylog, pos=pos, name=name)

    pap = plot(o3ap_site, p[0:23], color='red', /overplot, /buffer, $
      yrange=[1000, 0], /ylog, pos=pos, name='A priori')
    print, o3_site
    print, o3ap_site
      

    
    t0 = text(0.45, 0.25, 'FinalAlgorithmFlags:' + string(faf[xidxs[icity], yidxs[icity]], format='(I1)'))
    ;t1 = text(0.45, 0.25, 'Pressure Quality Flag:' + string(info[1], format='(I1)'))
    t2 = text(0.45, 0.20, 'Cloud Pressure:' + string(cp[xidxs[icity], yidxs[icity]]))
    t3 = text(0.45, 0.15, 'Effective Cloud Pressure:' + string(ecf[xidxs[icity], yidxs[icity]], format='(F7.2)'))

    leg = legend(target=[p1, pap], pos=leg_pos)
    pngfile = './plot/gems_l2_o3_profile_ml_clima_' +fitrange+'_prec'+ prec + $
      '_' + name + '.png'
    p1.save, pngfile
    p1.close
    ; send image to pc
    spawn, 'scp -P18742 -p ' + pngfile + $
      ' soodal@164.125.38.179:/home/soodal/works/plot/112x_447y'
    stop
  endfor
endfor

end

