
pro plot_softcal_rad, wav, $
  obs_rad, synthetic_rad, $
  obs_irrad, synthetic_irrad, $
  xpos, ypos

  pos = [0.12, 0.1, 0.9, 0.85]
  leg_pos = [0.85, 0.8]

  outfilepath = './plot/synthetic_radiance_diff/'
  if not file_test(outfilepath) then begin
    file_mkdir, outfilepath
  endif

  yrange = [0, 3.]
  p1 = plot(wav[xpos, ypos, *], obs_rad[xpos, ypos, *], /buffer, $
    axis_style=1, $
    xtitle='Wavelength[nm]', $
    ytitle='Radiance[W/cm^2/cm/sr]', $
    title='GEMS L1C Radiances with Synthetic Radiances', $
    ;symbol='s', $
    yrange=yrange, $
    xrange=[300, 340], $
    name='GEMS L1C Radiances', $
    color='red', $
    position=pos)
  p1['axis1'].color='red'

  ;plot simulated radiance
  p2 = plot($
    wav[xpos, ypos, *], synthetic_rad[xpos, ypos, *], $
    /buffer, $
    axis_style=0, $
    name='MERRA2 Synthetic Radiance', $
    ;symbol='s', $
    yrange=yrange, $
    xrange=[300, 340], $
    position=pos, $
    color='blue', $
    /current)
  ;yaxis2 = axis('Y', LOCATION='right', yrange=p2.yrange, TARGET=p2, $
    ;color='blue')
  xaxis2 = axis('X', LOCATION='TOP', xrange=p1.xrange, TARGET=p2, $
    color='black', $
    tickname=['', '', '', '', ''])

  ;plot observed radiance
  p3 = plot($
    wav[xpos, ypos, *], $
    (obs_rad[xpos, ypos, *] - synthetic_rad[xpos, ypos, *])/synthetic_rad[xpos, ypos, *], $
    /buffer, $
    axis_style=0, $
    name='Ratio', $
    ;symbol='s', $
    yrange=[-2, 2], $
    xrange=[300, 340], $
    position=pos, $
    color='black', $
    /current)
  yaxis2 = axis('Y', LOCATION='right', yrange=[-2, 2], TARGET=p3, $
    color='red')

  p_ = plot(indgen(500), fltarr(500), 'gray', /buffer, /current, $
    axis_style=0, $
    position=pos, $
    yrange=p3.yrange, $
    xrange=p1.xrange, $
    color='gray')
  leg = legend(target=[p1, p2, p3], position=[340, 0.8], $
    /data, /auto_text_color)

  ;p3 = plot(wav[0:fitresdims[0]-1], stdres, 'r', /buffer, $
    ;/over, name='Stddev')

  ;t1 = text(0.25, 0.75, 'Effective Cloud Fraction < 0.2', /norm)
  ;t2 = text(0.25, 0.70, 'Latitude < 20', /norm)
  ;leg = legend(target=[p0, p1, p2, p3, pap], pos=leg_pos)
  sxidx = string(xpos, format='(i03)')
  syidx = string(ypos, format='(i03)')
  pngfile = outfilepath + 'gems_merra2_synthetic_radiance_x' +sXidx +'_y' +sYidx +'.png'
  p1.errorbar_color='indian_red'
  p1.errorbar_capsize=0.1
  p1.title.font_size=16
  p1.save, pngfile
  p1.close

  ; PLOT for IRRAD

  p1 = plot(wav[xpos, ypos, *], obs_irrad[xpos, ypos, *], /buffer, $
    axis_style=1, $
    xtitle='Wavelength[nm]', $
    ytitle='Irradiance[W/cm^2/cm/sr]', $
    title='GEMS L1C Irradiances with Synthetic Irradiances', $
    ;symbol='s', $
    yrange=yrange, $
    xrange=[300, 340], $
    name='GEMS L1C Radiances', $
    color='red', $
    position=pos)
  p1['axis1'].color='red'

  ;plot simulated irradiance
  p2 = plot($
    wav[xpos, ypos, *], synthetic_irrad[xpos, ypos, *], $
    /buffer, $
    axis_style=0, $
    name='MERRA2 Synthetic Irradiance', $
    ;symbol='s', $
    yrange=yrange, $
    xrange=[300, 340], $
    position=pos, $
    color='blue', $
    /current)
  ;yaxis2 = axis('Y', LOCATION='right', yrange=p2.yrange, TARGET=p2, $
    ;color='blue')
  xaxis2 = axis('X', LOCATION='TOP', xrange=p1.xrange, TARGET=p2, $
    color='black', $
    tickname=['', '', '', '', ''])

  ;plot observed radiance
  p3 = plot($
    wav[xpos, ypos, *], $
    (obs_irrad[xpos, ypos, *] - synthetic_irrad[xpos, ypos, *])/synthetic_irrad[xpos, ypos, *], $
    /buffer, $
    axis_style=0, $
    name='Ratio', $
    ;symbol='s', $
    yrange=[-2, 2], $
    xrange=[300, 340], $
    position=pos, $
    color='black', $
    /current)
  yaxis2 = axis('Y', LOCATION='right', yrange=[-2, 2], TARGET=p3, $
    color='red')

  p_ = plot(indgen(500), fltarr(500), 'gray', /buffer, /current, $
    axis_style=0, $
    position=pos, $
    yrange=p3.yrange, $
    xrange=p1.xrange, $
    color='gray')
  leg = legend(target=[p1, p2, p3], position=[340, 0.8], $
    /data, /auto_text_color)

  ;p3 = plot(wav[0:fitresdims[0]-1], stdres, 'r', /buffer, $
    ;/over, name='Stddev')

  ;t1 = text(0.25, 0.75, 'Effective Cloud Fraction < 0.2', /norm)
  ;t2 = text(0.25, 0.70, 'Latitude < 20', /norm)
  ;leg = legend(target=[p0, p1, p2, p3, pap], pos=leg_pos)
  sxidx = string(xpos, format='(i03)')
  syidx = string(ypos, format='(i03)')
  pngfile = outfilepath + 'gems_merra2_synthetic_irradiance_x' +sXidx +'_y' +sYidx +'.png'
  p1.errorbar_color='indian_red'
  p1.errorbar_capsize=0.1
  p1.title.font_size=16
  p1.save, pngfile
  p1.close



  ; send image to pc
  ;scp_dest = 'soodal@164.125.38.179:/home/soodal/works/plot/'
  ;spawn, 'scp -P18742 -p ' + pngfile +  ' ' + scp_dest
end




jd_list = timegen(start=julday(6, 21, 2021, 3, 45), $
  final=julday(6, 21, 2021, 3, 45), units='Hours')

itime = 0
caldat, jd_list[itime], month, day, year, hour, minute
yyyy = string(year, format='(i04)')
mm = string(month, format='(i02)')
dd = string(day, format='(i02)')
hh = string(hour, format='(i02)')
mi = string(minute, format='(i02)')

datetime_str = yyyy + mm + dd + '_' + hh + mi

; SYNTHETIC
search_str = '/home/soodal/data/merra2_residual/residuals/' + $
  'GK2_GEMS_L2_O3P_' + datetime_str + '_wli300_prec000000_climML_b4x4_o20211002T123142Z_ecf0.nc'

fl = file_search(search_str)
gemso3p = ds_read_gems_l2_o3p(fl[0])
corrected_rad_all = gemso3p.corrected_rad_all
corrected_irrad_all = gemso3p.corrected_irrad_all
div_sun_all = gemso3p.div_sun_all
div_rad_all = gemso3p.div_rad_all

wav = gemso3p.WavelengthsWholeRange[*, *, 0:203]

; GEMS Observation
l1cradfn = '/data2/L1C_GEMS/L1C/4x4/GK2_GEMS_L1C_20210621_0345_NOR_694_4x4.nc'
irrfn = '/data2/L1C_GEMS/L1C/4x4/GK2_GEMS_IRR_20210620_4x4.nc'

gemsrad = ds_read_gems_l1c_4x4(l1cradfn)
gemsirr = ds_read_gems_l1c_4x4(irrfn)

rad = reform(gemsrad.image_pixel_values[*, *, 0:203])
sz = size(rad, /dimension)

gems_rad_norm = total(rad, 3)/204.
for iw=0, 203 do begin 
  rad[0:sz[0]-1, 0:sz[1]-1, iw] = rad[0:sz[0]-1, 0:sz[1]-1, iw] / gems_rad_norm
endfor

irr = reform(gemsirr.image_pixel_values[*, 0:203])
gems_irrad_norm = total(irr, 2)/204.
for iw=0, 203 do begin 
  irr[0:sz[1]-1, iw] = irr[0:sz[1]-1, iw] / gems_irrad_norm
endfor

irr3d = fltarr(sz[0], sz[1], sz[2])

for i=0, sz[0]-1 do begin
  irr3d[i, *, 0:203] = irr
endfor


for ix=30, 174, 30 do begin
  for iy=30, 512, 50 do begin 
    outfilename = 'x' + strtrim(string(ix, format='(i03)'), 2) + $
      '_y' + strtrim(string(iy, format='(i03)'), 2) + '.png'

    plot_softcal_rad, wav, $
      rad, corrected_rad_all, $
      irr3d, corrected_irrad_all, $
      ix, iy
      
  endfor
endfor

end
