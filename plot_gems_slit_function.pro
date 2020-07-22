read_gems_srf_interpolated_on = 1
plot_gems_srf_interploated_on = 1

; read EWHA interpolated SRF
if read_gems_srf_interpolated_on then begin
  restore, 'data/SRF/Slit function/GEMS_region_bpdata_all_polyfit.sav'
endif
yidx = 1023
centeridx = 30

nwl = 1033

ncolors = 256
ct = colortable([[100, 0, 100], [0, 0, 100], [0, 0, 255]], $
  ncolors=levels,$
  /transpose)

c1 = contour(reform(gems_bp_peaknorm[1023, *, *]), $
  title='GEMS Slit Function', $
  pos=[0.15, 0.2, 0.9, 0.7], $
  /fill, $
  xtitle='Slit', $
  ytitle='Wavelength idx', $
  n_levels=100, $
  rgb_table=10, $
  /buffer)

pngfile = 'fig/gems_slit_function.png'
c1.save, pngfile
c1.close

; send image to pc
spawn, 'scp -P18742 ' + pngfile + $
  ' soodal@164.125.38.179:/home/soodal/WindowsHome/Pictures/scp'

if plot_gems_srf_interploated_on then begin
  pos=[0.1, 0.1, 0.9, 0.9]

  coloridx = fix(jwl/1033.*256.)
  for jwl=0, nwl-1 do begin
    x = gems_laserwv[*, jwl]
    y = fltarr(61)
    y[*] = gems_laserwv[30, jwl]
    z = reform(gems_bp_peaknorm[yidx, *, jwl])
    p=plot3d(x,y,z, color=coloridx, xrange=[300, 500], yrange=[300, 500], $
      title='GEMS Slit Function', $
      axis_style=2, margin=[0.2, 0.2, 0.2, 0], $
      xminor=0, yminor=0, zminor=0, $
      rgb_table=ct, $
      ;vert_colors=bytscl(t)
      ;shadow_color='white', $
      xy_shadow=1, yz_shadow=1, xz_shadow=1, $
      xtitle='Wavelength', $
      ytitle='Center Wavelength', $
      overplot=jwl, $
      /buffer)


    ;p2=barplot(gems_wl_mean, (gems_irr_mean-ref_GEMSs)/ref_GEMSs, /buffer, $
      ;pos=pos, $
      ;color='red', $
      ;name='Difference', $
      ;axis_style=0, $
      ;/current)
  endfor
  
  ;yaxis2=axis('Y', LOCATION='right', yrange=p2.yrange, TARGET=p2, $
      ;color='red')
  ;xaxis2=axis('X', LOCATION='top', xrange=p2.xrange, TARGET=p2)

  pngfile = 'fig/gems_slit_function_3d.png'
  p.save, pngfile
  p.close

  ; send image to pc
  spawn, 'scp -P18742 ' + pngfile + $
    ' soodal@164.125.38.179:/home/soodal/WindowsHome/Pictures/scp'
end




end
