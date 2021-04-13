pro plot_gems_l2o3p_vertical_slice_profile, $
  z_xdir, x_xdir, y_xdir, z_ydir, x_ydir, y_ydir, $
  outpath, basename, $
  pixlon, pixlat, $
  do_polyfill=do_polyfill, do_contour=do_contour, title=title, range=range

if (not keyword_Set(do_polyfill)) and (not keyword_Set(do_contour)) then begin
  print, 'do_polyfill and do_contour keyword are not set'
  print, 'plotting with polyfill'
  do_polyfill = 1
endif

if do_contour then begin
  do_polyfill = 0
endif

if not keyword_Set(title) then begin
  title = 'GEMS L2 Ozone profile'
endif

; collocate
;cities_name = ['pohang']
;pohang_lon = 129.37963
;pohang_lat = 36.03259

;o3 = o3p.o3
;pres = o3p.pressure

if not keyword_Set(range) then begin
  du_min = min([z_xdir, z_ydir])
  du_max = max([z_xdir, z_ydir])
  range = [du_min, du_max]
endif
du_min = min(range)
du_max = max(range)

szx = size(z_xdir, /dimension)
szy = size(z_ydir, /dimension)

dim1 = szx[0]
dim2 = szy[0]
dim3 = szx[1]

; ==============================================================================
; plotting for x direction
; ==============================================================================

xrange = [min(x_xdir), max(x_xdir)]
yrange = [1000, 0.1]
n_levels = 100
levels = findgen(n_levels)*max(range)/n_levels

pngfile = outpath + basename + '.ph_lat' + $
  string(pixlat, format='(F6.3)') + '.png'

pos = [0.13, 0.10, 0.85, 0.90]
cbpos = [0.87, 0.10, 0.90, 0.90]

if do_polyfill then begin
  cgPS_Open, 'gems.ps'
  cgDisplay, aspect=0.7


  charsize = (!D.Name EQ 'PS') ? cgDefCharsize()*0.7 : cgDefCharsize()*0.75
  thick = (!D.Name EQ 'PS') ? 6 :3
  londel = 10
  latdel = 10

  dsloadct, 33

  device, decomposed=0
  !x.range=xrange
  !y.range=yrange
  !x.style=1
  !y.style=1
  contour, z_xdir, x_xdir[*, 0:-2], y_xdir[*, 0:-2],/fill, /ylog, $
    position=pos, $
    xtitle='Longitude', $
    xtickv=findgen(6)*10+100, $
    xtickname=['100E', '110E', '120E', '130E', '140E', '150E'], $
    ytitle='Pressure[hPa]', $
    title=title, $
    xrange=xrange, yrange=yrange, zrange=[0, 50], levels =levels, /nodata

  for iy=0, dim3-1-1 do begin
    for ix=0, dim1-1-1 do begin
      xbox = [x_xdir[ix,iy], x_xdir[ix+1,iy], x_xdir[ix+1,iy], x_xdir[ix,iy]]
      ybox = [y_xdir[ix,iy], y_xdir[ix,iy], y_xdir[ix,iy+1], y_xdir[ix,iy+1]]
      oob_idx = where(ybox lt 0.1, /null)
      ybox[oob_idx] = 0.1

      polyfill, xbox, ybox, $
        color=bytscl(z[ix, iy], min=du_min, max=du_max, top=253), $
        /data
    endfor
  endfor

  cgColorbar, format='(f6.1)', range=range, $
    Position=[pos[2]+0.1, pos[1], pos[2]+0.13, pos[3]], $
    Color='black',$
    OOB_High=253, $
    OOB_low=0,$
    divisions=5, $
    bottom=0, $
    tlocation='LEFT', $
    Ticknames=['0', '10', '20', '30', '40', '50'], $
    ncolor=254, $
    minor=1,  $
    charsize=1.5,$
    XTicklen=1,  $
    XMinor=0, $
    /vertical, $
    AnnotateColor='black',  $
    Title='Ozone[DU]'

  cgPS_Close

  cgPS2Raster, 'gems.ps', pngfile , /png, density = 1000
endif

if do_contour then begin
  ct = colortable(72, /reverse)
  c = contour(z_xdir, x_xdir[*, 0:-2], y_xdir[*, 0:-2], /fill, /ylog, $
    position=pos, $
    xtitle='Longitude', $
    xtickv=findgen(6)*10+100, $
    xtickname=['100E', '110E', '120E', '130E', '140E', '150E'], $
    ytitle='Pressure[hPa]', $
    title=title, $
    rgb_table=ct, $
    xrange=xrange, yrange=yrange, zrange=range, c_value=levels, $
    /buffer)
  cb = colorbar(target=c, /taper, /orientation, position=cbpos, /textpos)
    ;/border)
  p = plot([pixlon], [900],/overplot, SYMBOL='Star', /sym_filled, $
    sym_fill_color='Purple')
  c.save, pngfile
  c.close

endif

if not keyword_Set(scp_dest1) then begin
  scp_dest = 'soodal@164.125.38.179:/home/soodal/works/plot/paper'
endif else begin
  scp_dest = scp_dest1
endelse

; send image to pc
;if keyword_set(scp) then begin
  spawn, 'scp -P18742 -p ' + pngfile + $
  ' ' + scp_dest
;endif


; ==============================================================================
; plotting for y direction
; ==============================================================================

xrange = [min(x_ydir), max(x_ydir)]
yrange = [1000, 0.1]
n_levels = 100
levels = findgen(n_levels)*max(range)/n_levels

pngfile = outpath + basename + '.ph_lon' + $
  string(pixlon, format='(F7.3)') + '.png'

pos = [0.13, 0.10, 0.85, 0.90]
cbpos = [0.87, 0.10, 0.90, 0.90]

if do_polyfill then begin
  cgPS_Open, 'gems.ps'
  cgDisplay, aspect=0.7


  charsize = (!D.Name EQ 'PS') ? cgDefCharsize()*0.7 : cgDefCharsize()*0.75
  thick = (!D.Name EQ 'PS') ? 6 :3
  londel = 10
  latdel = 10

  dsloadct, 33

  device, decomposed=0
  !x.range=xrange
  !y.range=yrange
  !x.style=1
  !y.style=1
  contour, z_ydir, x_ydir[*, 0:-2], y_ydir[*, 0:-2],/fill, /ylog, $
    position=pos, $
    xtitle='Latitude', $
    xtickv=findgen(5)*10, $
    xtickname=['0N', '10N', '20N', '30N', '40N'], $
    ytitle='Pressure[hPa]', $
    title=title, $
    xrange=xrange, yrange=yrange, zrange=[0, 50], levels =levels, /nodata

  for iy=0, dim3-1-1 do begin
    for ix=0, dim2-1-1 do begin
      xbox = [x_ydir[ix,iy], x_ydir[ix+1,iy], x_ydir[ix+1,iy], x_ydir[ix,iy]]
      ybox = [y_ydir[ix,iy], y_ydir[ix,iy], y_ydir[ix,iy+1], y_ydir[ix,iy+1]]
      oob_idx = where(ybox lt 0.1, /null)
      ybox[oob_idx] = 0.1

      polyfill, xbox, ybox, $
        color=bytscl(z[ix, iy], min=du_min, max=du_max, top=253), $
        /data
    endfor
  endfor

  cgColorbar, format='(f6.1)', range=range, $
    Position=[pos[2]+0.1, pos[1], pos[2]+0.13, pos[3]], $
    Color='black',$
    OOB_High=253, $
    OOB_low=0,$
    divisions=5, $
    bottom=0, $
    tlocation='LEFT', $
    Ticknames=['0', '10', '20', '30', '40', '50'], $
    ncolor=254, $
    minor=1,  $
    charsize=1.5,$
    XTicklen=1,  $
    XMinor=0, $
    /vertical, $
    AnnotateColor='black',  $
    Title='Ozone[DU]'

  cgPS_Close
  cgPS2Raster, 'gems.ps', pngfile , /png, density = 1000

endif

if do_contour then begin
  ct = colortable(72, /reverse)
  c = contour(z_ydir, x_ydir[*, 0:-2], y_ydir[*, 0:-2], /fill, /ylog, $
    position=pos, $
    xtitle='Latitude', $
    xtickv=findgen(5)*10, $
    xtickname=['0N', '10N', '20N', '30N', '40N'], $
    ytitle='Pressure[hPa]', $
    title=title, $
    rgb_table=ct, $
    xrange=xrange, yrange=yrange, zrange=range, c_value=levels, $
    /buffer)
  cb = colorbar(target=c, /taper, /orientation, position=cbpos, /textpos)
    ;/border)
  p = plot([pixlat], [900],/overplot, SYMBOL='Star', /sym_filled, $
    sym_fill_color='Purple')
  c.save, pngfile
  c.close

endif


if not keyword_Set(scp_dest1) then begin
  scp_dest = 'soodal@164.125.38.179:/home/soodal/works/plot/paper'
endif else begin
  scp_dest = scp_dest1
endelse

; send image to pc
;if keyword_set(scp) then begin
  spawn, 'scp -P18742 -p ' + pngfile + $
  ' ' + scp_dest
;endif
end
