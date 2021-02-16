
cities_name = ['Beijing', 'Chungching', 'Shanghai', 'Gwangjou', 'Seoul', 'Busan']
xidxs = [105, 144, 93, 123, 74, 68]
yidxs = [38, 129, 111, 189, 57, 76]

pos = [0.10, 0.1, 0.8, 0.85]
leg_pos = [0.85, 0.8]


fitrange = '300340'
fit_range = strmid(fitrange, 0, 3) + '-' + strmid(fitrange, 3, 3)
L1Cmaker = 'NIER'

path ='/data1/gems/o3p/ds/GEMS_O3P_Yonsei/' 
outputpath = path + 'out/'
projectpath = outputpath + 'ml_clima/' + fitrange + '_' + L1Cmaker + '/'

filelist = file_search(projectpath + '*.nc4')
print, filelist

runtime = strmid(filelist, 13, 10, /reverse)
runtimeidx = sort(runtime)
recentrunfile = filelist[runtimeidx[-1]]

fn = file_basename(recentrunfile)

date = strmid(fn, 13, 13)

scp_dest = 'soodal@164.125.38.179:/home/soodal/WindowsHome/Pictures/scp/' $
  + 'ml_clima/' + fitrange + '_' + L1Cmaker 

limit=[-5, 85, 55, 155]
nlayer = 24

zrange1 = [0.15, 0.15, 0.25, 0.5, 0.9, 1.6, 2.9, 5, 8, 12, $
  15, 14, 8, 6, 14, 12, 2, 0, 0, 0, 0, 0, 4, 5]

zrange2 = [0.3, 0.2, 0.35, 0.7, 1.2, 2.2, 3.8, 7, 11, 20, $
  25, 40, 45, 45, 35, 35, 30, 40, 40, 20, 10, 15, 20, 25]

data = ds_read_gems_l2_o3p(projectpath + fn)


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

for icity = 0, n_elements(cities_name)-1 do begin

  ;color = ['red', 'blue', 'green', 'cyan'] 
  name = cities_name[icity]

  o3_site = o3[*, xidxs[icity], yidxs[icity]] 
  o3ap_site = o3ap[*, xidxs[icity], yidxs[icity]]
  p = pressure[*, xidxs[icity], yidxs[icity]]


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
    yrange=[1000, 0], /ylog, pos=pos, name=name)

  pap = plot(o3ap_site, p[0:23], color='red', /overplot, /buffer, $
    yrange=[1000, 0], /ylog, pos=pos, name='A priori')
    

  
  t0 = text(0.45, 0.25, 'Final Quality Flag:' + string(faf[xidxs[icity], yidxs[icity]], format='(I1)'))
  ;t1 = text(0.45, 0.25, 'Pressure Quality Flag:' + string(info[1], format='(I1)'))
  t2 = text(0.45, 0.20, 'Cloud Pressure:' + string(cp[xidxs[icity], yidxs[icity]]))
  t3 = text(0.45, 0.15, 'Effective Cloud Pressure:' + string(ecf[xidxs[icity], yidxs[icity]], format='(F7.2)'))

  leg = legend(target=[p1, pap], pos=leg_pos)
  pngfile = './plot/gems_l2_o3_profile_ml_clima_' +fitrange+'_'+L1Cmaker +'_'+ name + '.png'
  p1.save, pngfile
  p1.close
  ; send image to pc
  spawn, 'scp -P18742 -p ' + pngfile + $
    ' soodal@164.125.38.179:/home/soodal/WindowsHome/Pictures/scp/profiles'
endfor
end




