; INDEX
; 310340_KARI
;   for cao3
;   for strato
;   for tropo
; 310340_EOSRL
;   for cao3
;   for strato
;   for tropo
; 305340_KARI
;   for cao3
;   for strato
;   for tropo
; 305340_EOSRL
;   for cao3
;   for strato
;   for tropo

; output file calculated by me
; 
;fn = '../GEMS_O3P_Yonsei/out/310340_KARI/GK2B_GEMS_L2_20160115_0300_O3P_ND_DPRO_first.4x4_2008102141.nc4'
;fn = '../GEMS_O3P_Yonsei/out/310340_EOSRL/GK2B_GEMS_L2_20160115_0300_O3P_ND_DPRO_first.4x4_2008101648.nc4'
;fn = '../GEMS_O3P_Yonsei/out/305340_KARI/GK2B_GEMS_L2_20160115_0300_O3P_ND_DPRO_first.4x4_2008111012.nc4'
;fn = '../GEMS_O3P_Yonsei/out/305340_EOSRL/GK2B_GEMS_L2_20160115_0300_O3P_ND_DPRO_first.4x4_2008112047.nc4'

; output file calculated by hs
;fn = '/data1/L2_GEMS/o3p/4x4/nopolc/eosrl/310-340/GEMS_O3P_20200616_0345_nopolc_2008111857.nc4'
;fn = '/data1/L2_GEMS/o3p/4x4/nopolc/eosrl/310-340/GEMS_O3P_20200616_0345_nopolc_2008111857.nc4'
;fn = '../GEMS_O3P_Yonsei/out/310340_EOSRL/GK2B_GEMS_L2_20160115_0300_O3P_ND_DPRO_first.4x4_2008101648.nc4'
;fn = '../GEMS_O3P_Yonsei/out/305340_KARI/GK2B_GEMS_L2_20160115_0300_O3P_ND_DPRO_first.4x4_2008111012.nc4'
;fn = '../GEMS_O3P_Yonsei/out/305340_EOSRL/GK2B_GEMS_L2_20160115_0300_O3P_ND_DPRO_first.4x4_2008112047.nc4'


; Configuration

pos = [0.10, 0.1, 0.8, 0.85]
leg_pos = [0.85, 0.8]

fitrange = '305340'
fit_range = strmid(fitrange, 0, 3) + '-' + strmid(fitrange, 3, 3)
L1Cmaker = 'EOSRL'

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

; for CAO3

map = MAP('Orthographic', /BUFFER, $
  limit=limit, $
  CENTER_LATITUDE=25, $
  CENTER_LONGITUDE=125)

zrange = [250, 380]
levels = findgen(zrange[1] - zrange[0] + 1, start=zrange[0])

ct = colortable(74, /reverse)
mag = bytscl(o3total[o3total_geolocidx], min=zrange[0], max=zrange[1])
p1 = scatterplot(lon[o3total_geolocidx], lat[o3total_geolocidx], $
  magnitude=mag, /buffer, $
  overplot=map, $
  rgb_table=ct, $
  clip=0, $
  position=pos, $
  symbol='dot', $
  sym_size=0.7, $
  title='Total Column Amount O3 EOSRL L1C ' + fit_range + ' nm fitting', $
  axis_style=0)

grid = map.MAPGRID
grid.box_axes=1
grid.LINESTYLE = "dashed"
grid.LABEL_POSITION = 0
map['latitude'].LABEL_ANGLE=80
map['longitude'].LABEL_ANGLE=0
grid.longitude_min=limit[1]
grid.latitude_min=limit[0]
grid.grid_Latitude=10
grid.grid_Longitude=10
grid.FONT_SIZE=10
 
m1 = MAPCONTINENTS(/hires)
     
tickvalues = linspace(zrange[0], zrange[1], 10)
cb = colorbar( $
  orientation=1, $
  position=[0.93, 0.1, 0.97, 0.85], $
  tickvalues=tickvalues,$
  range=[zrange[0], zrange[1]], $
  color='black', $
  title='O3[DU]', $
  taper=1)

caldat, systime(/julian), mo, dd, yyyy, hh, mi, ss
isodate = string(yyyy, format='(I04)') + '-' +string(mo, format='(i02)') $
  + '-' + string(dd, format='(i02)') + 'T' + string(hh, format='(i02)') $
  + ':' + string(mi, format='(i02)') + ':' + string(ss, format='(i02)')

t = text(0.4, 0.025, 'Dae Sung Choi, ' + isodate)
pngfile = './plot/gems_l2_cao3_' + fitrange + '_' + L1Cmaker +'.png'

p1.title.font_size=16
p1.save, pngfile
p1.close
; send image to pc
spawn, 'scp -P18742 -p ' + pngfile + $
  ' ' + scp_dest 

; for Stratospheric O3

map = MAP('Orthographic', /BUFFER, $
  limit=limit, $
  CENTER_LATITUDE=25, $
  CENTER_LONGITUDE=125)

zrange = [190, 320]
levels = findgen(zrange[1] - zrange[0] + 1, start=zrange[0])

ct = colortable(74, /reverse)
mag = bytscl(o3strato[o3strato_geolocidx], min=zrange[0], max=zrange[1])
p1 = scatterplot(lon[o3strato_geolocidx], lat[o3strato_geolocidx], $
  magnitude=mag, /buffer, $
  overplot=map, $
  rgb_table=ct, $
  clip=0, $
  position=pos, $
  symbol='dot', $
  sym_size=0.7, $
  title='Stratospheric Column Amount O3 EOSRL L1C ' + fit_range + ' nm fitting', $
  axis_style=0)

grid = map.MAPGRID
grid.box_axes=1
grid.LINESTYLE = "dashed"
grid.LABEL_POSITION = 0
map['latitude'].LABEL_ANGLE=80
map['longitude'].LABEL_ANGLE=0
grid.longitude_min=limit[1]
grid.latitude_min=limit[0]
grid.grid_Latitude=10
grid.grid_Longitude=10
grid.FONT_SIZE=10
 
m1 = MAPCONTINENTS(/hires)
     
tickvalues = linspace(zrange[0], zrange[1], 10)
cb = colorbar( $
  orientation=1, $
  position=[0.93, 0.1, 0.97, 0.85], $
  tickvalues=tickvalues,$
  range=[zrange[0], zrange[1]], $
  color='black', $
  title='O3[DU]', $
  taper=1)

caldat, systime(/julian), mo, dd, yyyy, hh, mi, ss
isodate = string(yyyy, format='(I04)') + '-' +string(mo, format='(i02)') $
  + '-' + string(dd, format='(i02)') + 'T' + string(hh, format='(i02)') $
  + ':' + string(mi, format='(i02)') + ':' + string(ss, format='(i02)')

t = text(0.4, 0.025, 'Dae Sung Choi, ' + isodate)
pngfile = './plot/gems_l2_o3strato_' + fitrange + '_' + L1Cmaker + '.png'

p1.title.font_size=16
p1.save, pngfile
p1.close
; send image to pc
spawn, 'scp -P18742 -p ' + pngfile + $
  ' ' + scp_dest

; for Tropospheric O3

map = MAP('Orthographic', /BUFFER, $
  limit=limit, $
  CENTER_LATITUDE=25, $
  CENTER_LONGITUDE=125)

zrange = [20, 60]
levels = findgen(zrange[1] - zrange[0] + 1, start=zrange[0])

ct = colortable(74, /reverse)
mag = bytscl(o3tropo[o3tropo_geolocidx], min=zrange[0], max=zrange[1])
p1 = scatterplot(lon[o3tropo_geolocidx], lat[o3tropo_geolocidx], $
  magnitude=mag, /buffer, $
  overplot=map, $
  rgb_table=ct, $
  clip=0, $
  position=pos, $
  symbol='dot', $
  sym_size=0.7, $
  title='Troposheric Column Amount O3 EOSRL L1C ' + fit_range + ' nm fitting', $
  axis_style=0)

grid = map.MAPGRID
grid.box_axes=1
grid.LINESTYLE = "dashed"
grid.LABEL_POSITION = 0
map['latitude'].LABEL_ANGLE=80
map['longitude'].LABEL_ANGLE=0
grid.longitude_min=limit[1]
grid.latitude_min=limit[0]
grid.grid_Latitude=10
grid.grid_Longitude=10
grid.FONT_SIZE=10
 
m1 = MAPCONTINENTS(/hires)
     
tickvalues = linspace(zrange[0], zrange[1], 10)
cb = colorbar( $
  orientation=1, $
  position=[0.93, 0.1, 0.97, 0.85], $
  tickvalues=tickvalues,$
  range=[zrange[0], zrange[1]], $
  color='black', $
  title='O3[DU]', $
  taper=1)

caldat, systime(/julian), mo, dd, yyyy, hh, mi, ss
isodate = string(yyyy, format='(I04)') + '-' +string(mo, format='(i02)') $
  + '-' + string(dd, format='(i02)') + 'T' + string(hh, format='(i02)') $
  + ':' + string(mi, format='(i02)') + ':' + string(ss, format='(i02)')

t = text(0.4, 0.025, 'Dae Sung Choi, ' + isodate)
pngfile = './plot/gems_l2_o3tropo_' + fitrange + '_' + L1Cmaker + '.png'

p1.title.font_size=16
p1.save, pngfile
p1.close
; send image to pc
spawn, 'scp -P18742 -p ' + pngfile + $
  ' ' + scp_dest

; for levels

O3 = data.O3
pres = data.PRESSURE

;for i=0, nlayer-1 do begin
  ;iO3 = reform(O3[i, *, *], dim)
  ;o3layer_geolocidx = where(lon gt 0 and lat gt 0 and iO3 gt 0, /null)
  
  ;iPres = reform(pres[i, *, *], dim)
  ;nanidx = where(iPres lt 0, /null)
  ;iPres[nanidx] = !values.f_nan
  ;meanpres = mean(iPres, /nan)
  ;strpres = string(meanpres, format='(f07.2)')

  ;map = MAP('Orthographic', /BUFFER, $
    ;limit=limit, $
    ;CENTER_LATITUDE=25, $
    ;CENTER_LONGITUDE=125)

  ;zrange = [zrange1[i], zrange2[i]]
  ;levels = findgen(zrange[1] - zrange[0] + 1, start=zrange[0])

  ;ct = colortable(74, /reverse)
  ;mag = bytscl(iO3[o3layer_geolocidx], min=zrange[0], max=zrange[1])
  ;p1 = scatterplot(lon[o3layer_geolocidx], lat[o3layer_geolocidx], $
    ;magnitude=mag, /buffer, $
    ;overplot=map, $
    ;rgb_table=ct, $
    ;clip=0, $
    ;position=pos, $
    ;symbol='dot', $
    ;sym_size=0.7, $
    ;title='O3 ' + strpres + 'hPa (mean) EOSRL L1C ' + fit_range + ' nm fitting', $
    ;axis_style=0)

  ;grid = map.MAPGRID
  ;grid.box_axes=1
  ;grid.LINESTYLE = "dashed"
  ;grid.LABEL_POSITION = 0
  ;map['latitude'].LABEL_ANGLE=80
  ;map['longitude'].LABEL_ANGLE=0
  ;grid.longitude_min=limit[1]
  ;grid.latitude_min=limit[0]
  ;grid.grid_Latitude=10
  ;grid.grid_Longitude=10
  ;grid.FONT_SIZE=10
   
  ;m1 = MAPCONTINENTS(/hires)
       
  ;tickvalues = linspace(zrange[0], zrange[1], 10)
  ;cb = colorbar( $
    ;orientation=1, $
    ;position=[0.93, 0.1, 0.97, 0.85], $
    ;tickvalues=tickvalues,$
    ;range=[zrange[0], zrange[1]], $
    ;color='black', $
    ;title='O3[DU]', $
    ;taper=1)

  ;caldat, systime(/julian), mo, dd, yyyy, hh, mi, ss
  ;isodate = string(yyyy, format='(I04)') + '-' +string(mo, format='(i02)') $
    ;+ '-' + string(dd, format='(i02)') + 'T' + string(hh, format='(i02)') $
    ;+ ':' + string(mi, format='(i02)') + ':' + string(ss, format='(i02)')

  ;t = text(0.4, 0.025, 'Dae Sung Choi, ' + isodate)
  ;pngfile = './plot/gems_l2_cao3_310340_EOSRL_'+strpres+'.png'

  ;p1.title.font_size=16
  ;p1.save, pngfile
  ;p1.close
  ;; send image to pc
  ;spawn, 'scp -P18742 -p ' + pngfile + $
    ;' ' + scp_dest
;endfor
end
