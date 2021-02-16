
pos = [0.10, 0.1, 0.8, 0.85]
leg_pos = [0.85, 0.8]

fitrange = '300340'
fit_range = strmid(fitrange, 0, 3) + '-' + strmid(fitrange, 3, 3)
L1Cmaker = 'NIER'

path ='/data1/gems/o3p/ds/GEMS_O3P_Yonsei/' 
outputpath = path + 'out/'
projectpath = outputpath + 'softcal/' + fitrange + '_' + L1Cmaker $
  + '_image_mean'

filelist = file_search(projectpath + '*.nc4')
print, filelist

runtime = strmid(filelist, 13, 10, /reverse)
runtimeidx = sort(runtime)
recentrunfile = filelist[runtimeidx[-1]]

fn = file_basename(recentrunfile)

date = strmid(fn, 13, 13)

scp_dest = 'soodal@164.125.38.179:/home/soodal/WindowsHome/Pictures/scp/' $
  + 'softcal/' + fitrange + '_' + L1Cmaker + '_image_mean'

limit=[-5, 85, 55, 155]
nlayer = 24

zrange1 = [0.15, 0.15, 0.25, 0.5, 0.9, 1.6, 2.9, 5, 8, 12, $
  15, 14, 8, 6, 14, 12, 2, 0, 0, 0, 0, 0, 4, 5]

zrange2 = [0.3, 0.2, 0.35, 0.7, 1.2, 2.2, 3.8, 7, 11, 20, $
  25, 40, 45, 45, 35, 35, 30, 40, 40, 20, 10, 15, 20, 25]


;fn = '/data1/L2_GEMS/o3p/4x4/nopolc/eosrl/310-340/GEMS_O3P_20200616_0345_nopolc_2008111857.nc4'

data = ds_read_gems_l2_o3p(projectpath + fn)

cao3 = data.ColumnAmountO3
o3total = reform(cao3[0, *, *], 174, 512)
o3strato = reform(cao3[1, *, *], 174, 512)
o3tropo = reform(cao3[2, *, *], 174, 512)

latitude = data.Latitude
longitude = data.Longitude
dim = n_elements(latitude)

lat = reform(latitude, dim)
lon = reform(longitude, dim)

o3total_geolocidx = where(lon gt 0 and lat gt 0 and o3total gt 0, /null)
o3tropo_geolocidx = where(lon gt 0 and lat gt 0 and o3tropo gt 0, /null)
o3strato_geolocidx = where(lon gt 0 and lat gt 0 and o3strato gt 0, /null)


for ici = 0, n_elements(ci)-1 do begin
  for icj = 0, n_elements(cs)-1 do begin
    cics = 'ci' + strtrim(string(ci[ici], format='(i3)'), 2) + 'cs' + string(cs[icj], format='(i3)')

    color = ['red', 'blue', 'green', 'cyan'] 
    name = ['ME 0.01', '0.001', '0.005', '0.04']
    o3p_total = fltarr(4, 24)

    for ime = 0, 3 do begin
      fp = path + cics + '/' + string(ime, format='(i03)') + '/'  
      o3fn = 'output_l2o3p_ozret.txt'
      presfn = 'output_l2o3p_p.txt'
      extrafn = 'output_l2o3p_extra_info.txt'
      ozapfn = 'output_l2o3p_ozap.txt'

      readcol, fp + o3fn, o3, format='F', numline=24
      readcol, fp + presfn, p, format='F', numline=24
      readcol, fp + extrafn, info, numline=4
      readcol, fp + ozapfn, ozap, format='F', numline=24

      o3p_total[ime,*] = o3
      if ime eq 0 then begin
        command = "p" + string(ime, format='(i1)') + $
          " = plot(o3, p, color=color[ime], overplot=ime, /buffer, " + $
          " title='GEMS O3PROFILE 2020-06-16T0345', " + $
          " xtitle='O3[DU]', " + $
          " ytitle='Pressure[hPa]', " + $
          " yrange=[1000, 0], /ylog, pos=pos, name=name[ime])"
        print, command
      endif else begin
        command = "p" + string(ime, format='(i1)') + $
          " = plot(o3, p, color=color[ime], overplot=ime, /buffer, " + $
          " yrange=[1000, 0], /ylog, pos=pos, name=name[ime])"
      endelse
      dumm = execute(command)

    endfor

    pap = plot(o3, p, color='black', /overplot, /buffer, $
      yrange=[1000, 0], /ylog, pos=pos, name='A priori')
      
    t0 = text(0.45, 0.3, 'Final Quality Flag:' + string(info[0], format='(I1)'))
    t1 = text(0.45, 0.25, 'Pressure Quality Flag:' + string(info[1], format='(I1)'))
    t2 = text(0.45, 0.20, 'Cloud Pressure:' + string(info[2]))
    t3 = text(0.45, 0.15, 'Effective Cloud Pressure:' + string(info[3], format='(F7.2)'))

    leg = legend(target=[p0, p1, p2, p3, pap], pos=leg_pos)
    pngfile = './plot/gems_l2_o3_profile_' + cics + '.png'
    p1.save, pngfile
    p1.close
    ; send image to pc
    spawn, 'scp -P18742 -p ' + pngfile + $
      ' soodal@164.125.38.179:/home/soodal/WindowsHome/Pictures/scp/profiles'
    print, o3p_total
  endfor
endfor
end




