; Plot Ozone Spatial distribution from GEMS L2 O3P output file
pro ds_plot_gems_l2_o3p_all, input, title=title, $
  outputpath = outputpath, scppath = scppath, use_filename=use_filename, $
  use_structure=use_structure, basename=basename

; Set keyword
if not keyword_set(use_filename) then begin
  use_filename=0
endif
if not keyword_set(use_structure) then begin
  use_structure=0
endif
if not keyword_set(outputpath) then begin
  outputpath = './plot/'
endif

if not keyword_set(scppath) then begin
  scppath = 'soodal@164.125.38.179:/home/soodal/works/plot'
endif

; Set parameters
pos = [0.10, 0.1, 0.8, 0.85]
leg_pos = [0.85, 0.3]
limit=[-5, 80, 55, 160]
nlayer = 24

zrange1 = [0.15, 0.15, 0.25, 0.5, 0.9, 1.6, 2.9, 5, 8, 12, $
  15, 14, 8, 6, 14, 12, 2, 0, 0, 0, 0, 0, 4, 5]

zrange2 = [0.3, 0.2, 0.35, 0.7, 1.2, 2.2, 3.8, 7, 11, 20, $
  25, 40, 45, 45, 35, 35, 30, 40, 40, 20, 10, 15, 20, 25]


; read structure from nc4 file
if use_filename or (not use_structure) then begin
  basename = file_basename(input)
  data = ds_read_gems_l2_o3p(input)
endif

if use_structure then begin
  if not keyword_set(basename) then begin
    message, 'You need to set basename keyword.'
  endif
  data = input
endif


cao3 = data.ColumnAmountO3
o3total = reform(cao3[0, *, *], 174, 512)
o3strato = reform(cao3[1, *, *], 174, 512)
o3tropo = reform(cao3[2, *, *], 174, 512)

ds_gems_l2o3p_accum, data, o3under300, hpa=300.

latitude = data.Latitude
longitude = data.Longitude
dim = n_elements(latitude)

lat = reform(latitude, dim)
lon = reform(longitude, dim)

o3total_geolocidx = where(lon gt 0 and lat gt 0 and o3total gt 0, /null)
o3tropo_geolocidx = where(lon gt 0 and lat gt 0 and o3tropo gt 0, /null)
o3strato_geolocidx = where(lon gt 0 and lat gt 0 and o3strato gt 0, /null)

plot_sat_proj, o3total, longitude, latitude, $
  title=strmid(basename, 0, 24), $
  range=[230, 340], $
  pngfile=outputpath + basename+'_01total.png', $
  /scp_send, scp_dest=scppath

plot_sat_proj, o3strato, longitude, latitude, $
  title=strmid(basename, 0, 24), $
  range=[190, 320], $
  pngfile=outputpath + basename+'_02strato.png', $
  /scp_send, scp_dest=scppath

plot_sat_proj, o3tropo, longitude, latitude, $
  title=strmid(basename, 0, 24), $
  range=[20, 80], $
  pngfile=outputpath + basename+'_03tropo.png', $
  /scp_send, scp_dest=scppath

plot_sat_proj, o3under300, longitude, latitude, $
  title=strmid(basename, 0, 24), $
  range=[20, 60], $
  pngfile=outputpath + basename+'_04under300.png', $
  /scp_send, scp_dest=scppath

; for CAO3

orthographic = 0
if orthographic then begin
  map = MAP('Orthographic', /BUFFER, $
    limit=limit, $
    CENTER_LATITUDE=0, $
    CENTER_LONGITUDE=130)

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
    title='Total Column Amount O3', $
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
  pngfile = outputpath + basename + '_01_total_o3.png'

  p1.title.font_size=16
  p1.save, pngfile
  p1.close
  ; send image to pc
  spawn, 'scp -P18742 -p ' + pngfile + $
    ' '+scppath

  ; for Stratospheric Column O3

    map = MAP('Orthographic', /BUFFER, $
      limit=limit, $
      CENTER_LATITUDE=0, $
      CENTER_LONGITUDE=130)

    zrange = [190, 320]
    levels = findgen(zrange[1] - zrange[0] + 1, start=zrange[0])

    ct = colortable(74, /reverse)
    mag = bytscl(o3strato[o3strato_geolocidx], min=zrange[0], max=zrange[1])
    p1 = scatterplot(lon[o3strato_geolocidx], lat[o3strato_geolocidx],$
      magnitude=mag, /buffer, $
      overplot=map, $
      rgb_table=ct, $
      clip=0, $
      position=pos, $
      symbol='dot', $
      sym_size=0.7, $
      title='Stratospheric Column Amount O3', $
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
    pngfile = outputpath + basename + '_02_strato_o3.png'

    p1.title.font_size=16
    p1.save, pngfile
    p1.close
    ; send image to pc
    spawn, 'scp -P18742 -p ' + pngfile + $
      ' '+scppath

    ; for Tropospheric Column O3

    map = MAP('Orthographic', /BUFFER, $
      limit=limit, $
      CENTER_LATITUDE=0, $
      CENTER_LONGITUDE=130)

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
      title='Troposheric Column Amount O3', $
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
    pngfile = outputpath + basename + '_03_tropo_o3.png'

    p1.title.font_size=16
    p1.save, pngfile
    p1.close
    ; send image to pc
    spawn, 'scp -P18742 -p ' + pngfile + $
      ' '+scppath

    ; for levels

    ;O3 = data.O3
    ;pres = data.PRESSURE

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
        ;CENTER_LATITUDE=0, $
        ;CENTER_LONGITUDE=130)

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
        ;title='O3 ' + strpres + 'hPa (mean)', $
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
      ;pngfile = outputpath + basename + '_1' + string(i, format='(i02)') + '_total_o3.png'

      ;p1.title.font_size=16
      ;p1.save, pngfile
      ;p1.close
      ;; send image to pc
      ;spawn, 'scp -P18742 -p ' + pngfile + $
        ;' '+scppath
    ;endfor
endif
end
