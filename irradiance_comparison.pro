;=============================================================================;
compile_opt IDL2


read_gems_irr_on = 1
plot_gems_irr_on = 0

read_ref_solar_irr_on = 1
plot_ref_solar_irr_on = 1

read_gems_srf_on = 1
plot_gems_srf_on = 0

plot_bp_srf_on = 0

read_gems_srf_interpolated_on = 1

cal_gems_ref_irr_on = 1
plot_gems_ref_irr_on = 1
plot_gems_ref_irr_300_340_on = 1

; read GEMS irradiance data
if read_gems_irr_on then begin
  gems_irr = ds_read_nc('./GK2_GEMS_IRR_20200616.nc', $
    'image_pixel_values')
  gems_wl = ds_read_nc('./GK2_GEMS_IRR_20200616.nc', $
    'wavelength')
endif

; plot GEMS irradiance data
gems_wl_pixel = gems_wl[1023, *]
gems_irr_pixel = gems_wl[1023, *]
if plot_gems_irr_on then begin

  p = plot(gems_wl_pixel, gems_irr_pixel, /buffer, $
    title='GEMS irradiance 20200616 [W/cm2/cmgt/sr]')

  pngfile = 'gems_irr.png'
  p.save, pngfile
  p.close

  ; send image to pc
  spawn, 'scp -P18742 ' + pngfile + $
    ' soodal@164.125.38.179:/home/soodal/WindowsHome/Pictures/scp'
endif

; read Reference Solar irradiance data
c = 299792458. ; [m/s] speed of light
h = 6.62607015E-34 ; [j*s] planck constant
if read_ref_solar_irr_on then begin
  refsolspecfn='./OMSAO_NKPNO_SolarSpec_ReCal.dat'
  restore, './stemplate.sav'
  data = READ_ASCII(refsolspecfn, $
    TEMPLATE = sTemplate)
  flux = data.field2 * (h * c / data.field1) * $ ; data.field1 : [nm]
    1.0E14 * 1.E9 * 1E6 ; [E9j/s] 
  ; @TODO 1E6??
  refsolar = {wavelength:data.field1, photons:data.field2, irradiance:flux}
endif

; plot Reference Solar irradiance data
if plot_ref_solar_irr_on then begin
  p = plot(refsolar.wavelength, refsolar.photons, $
    title='OMSAO Reference Solar Spectrum, unit=[photons/s/cm2/nm]', $
    /buffer)

  pngfile = 'reference_solar_irradiance.png'
  p.save, pngfile
  p.close

  ; send image to pc
  spawn, 'scp -P18742 ' + pngfile + $
    ' soodal@164.125.38.179:/home/soodal/WindowsHome/Pictures/scp'
endif

; read SRF data
if read_gems_srf_on then begin
  restore, './SRF/Slit function/GEMS_bpdata.sav'
endif

; plot Slit Response Function for GEMS 7 channel
if plot_gems_srf_on then begin
  clevels = 7
  ctable = COLORTABLE([[0,0,100],$
                      [0,0,255],$
                      [0,255,255],$
                      [0,255,100]], $
                      NCOLORS = clevels, /TRANSPOSE)

  for ichan=0, 6 do begin
    if ichan ne 0 then overplot=1
    p = plot(laserwv[*, ichan], bp_peaknorm[0, *, 0], $
      title='Normalized Slit Response Function', $
      color=ctable[*, ichan], /buffer, overplot=overplot)
  endfor

  pngfile = 'srf.png'
  p.save, pngfile
  p.close

  ; send image to pc
  spawn, 'scp -P18742 ' + pngfile + $
    ' soodal@164.125.38.179:/home/soodal/WindowsHome/Pictures/scp'
endif


if plot_bp_srf_on then begin
  peaknorm = mean(bp_peaknorm, dim=1)
  ; yield solar irradiance using slit function with reference irradiance
  wlcenter = [301.8, 330.0, 365.0, 390.0, 435.0, 470.0, 498.2]

  satirr = fltarr(61, 7)
  for ichan=0, 6 do begin
    wls = []
    for i=0, 60 do begin
      wl = wlcenter[ichan] - 1.8 + 0.06*i
      wls = [wls, wl]
      refidx = long(round(wlcenter[ichan]*100))-24000L-181 + i*6 
      ;print, wl, refsolar.wavelength[refidx]
      satirr[i, ichan] = refsolar.irradiance[refidx] * peaknorm[i, ichan]
    endfor
    if ichan ne 0 then overplot=1
    p = plot(wls, satirr[*, ichan], $
      title='Reference Solar irradiance for 7 channel', $
      color=ctable[*, ichan], /buffer, overplot=overplot)
  endfor

  pngfile = 'satellite_irradiance.png'
  p.save, pngfile
  p.close

  ; send image to pc
  spawn, 'scp -P18742 ' + pngfile + $
    ' soodal@164.125.38.179:/home/soodal/WindowsHome/Pictures/scp'
endif

; read EWHA interpolated SRF
if read_gems_srf_interpolated_on then begin
  restore, './SRF/Slit function/GEMS_region_bpdata_all_polyfit.sav'
endif

; calculate gems reference irradiance
if cal_gems_ref_irr_on then begin
  GEMS_bp_peaknorm_mean = mean(GEMS_bp_peaknorm, dim=1)
  ref_GEMSs = []

  for iwlcenter = 0, n_elements(GEMS_laserwv[0, *])-1 do begin
    ; interpolate reference irradiance wavelength to gems wavelength
    refsolar_irr_gems = interpol(refsolar.irradiance, $
      refsolar.wavelength, GEMS_laserwv[*, iwlcenter])

    convolved = refsolar_irr_gems * GEMS_bp_peaknorm_mean[*, iwlcenter]
    ;print, total(convolved, /nan)
    ref_GEMSs = [ref_GEMSs, total(convolved, /nan)]
  endfor ; iwlcenter
endif

; plot interpolated GEMS SRF applied to reference irradiance spectrum
if plot_gems_ref_irr_on then begin
  pos=[0.1, 0.1, 0.9, 0.9]
  p=plot(gems_wl_pixel, ref_GEMSs, $
    pos=pos, $
    title='Reference GEMS Spectrum', $
    name='GEMS reference irradiance', $
    /buffer)

  p1=plot(gems_wl_pixel, gems_irr_pixel, /buffer, $
    pos=pos, $
    axis_style=1, $
    color='blue', $
    name='GEMS Irradiance', $
    /overplot)

  p2=barplot(gems_wl_pixel, (gems_irr_pixel-ref_GEMSs)/ref_GEMSs, /buffer, $
    pos=pos, $
    color='red', $
    name='Difference', $
    axis_style=0, $
    /current)
  
  yaxis2=axis('Y', LOCATION='right', yrange=p2.yrange, TARGET=p2, $
      color='red')
  xaxis2=axis('X', LOCATION='top', xrange=p2.xrange, TARGET=p2)

  pngfile = 'ref_GEMS_irradiance.png'
  p.save, pngfile
  p.close

  ; send image to pc
  spawn, 'scp -P18742 ' + pngfile + $
    ' soodal@164.125.38.179:/home/soodal/WindowsHome/Pictures/scp'
end

; plot interpolated GEMS SRF applied to reference irradiance spectrum
if plot_gems_ref_irr_300_340_on then begin
  pos=[0.1, 0.1, 0.9, 0.85]
  p=plot(gems_wl_pixel, ref_GEMSs, $
    pos=pos, $
    title='Reference GEMS Spectrum', $
    name='GEMS reference irradiance', $
    /buffer, $
    yrange=[-100, 1300], $
    xrange=[300, 340])

  p1=plot(gems_wl_pixel, gems_irr_pixel, /buffer, $
    pos=pos, $
    axis_style=1, $
    color='blue', $
    name='GEMS Irradiance', $
    xrange=[300, 340], $
    /overplot)

  p2=barplot(gems_wl_pixel, (gems_irr_pixel-ref_GEMSs)/ref_GEMSs, /buffer, $
    pos=pos, $
    color='red', $
    name='Difference', $
    axis_style=0, $
    yrange=[-0.5, 2.], $
    xrange=[300, 340], $
    /current)
  leg = legend(target=[p, p1, p2], position=[320, 1200], /data)
  
  yaxis2=axis('Y', LOCATION='right', yrange=p2.yrange, TARGET=p2, $
      color='red')
  xaxis2=axis('X', LOCATION='top', xrange=p2.xrange, TARGET=p2)
  pngfile = 'ref_GEMS_irradiance_300-340.png'
  p.save, pngfile
  p.close

  ; send image to pc
  spawn, 'scp -P18742 ' + pngfile + $
    ' soodal@164.125.38.179:/home/soodal/WindowsHome/Pictures/scp'
endif

end
