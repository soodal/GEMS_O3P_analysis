pro plot_gems_ak, o3p_file, xidx, yidx


;------------------------------------------------------------------------------
; Configuration
;------------------------------------------------------------------------------
pos = [0.12, 0.1, 0.9, 0.85]
leg_pos = [0.85, 0.8]

fitrange = '310340'
fit_range = strmid(fitrange, 0, 3) + '-' + strmid(fitrange, 3, 3)
L1Cmaker = 'EOSRL'

path ='/data1/gems/o3p/ds/GEMS_O3P_Yonsei/' 
outputpath = path + 'out/'
projectpath = outputpath ;+ 'softcal/' + fitrange + '/'

filelist = file_search(projectpath + '*.nc4')


runtime = strmid(filelist, 13, 10, /reverse)
runtimeidx = sort(runtime)
recentrunfile = filelist[runtimeidx[-1]]

filefullpath = '/data1/gems/o3p/ds/GEMS_O3P_Yonsei/out/GK2B_GEMS_L2_O3P_20200616_0345_climML_winliminit310_2020-12-05T0626.nc4'


;fn = file_basename(recentrunfile)
fn = file_basename(filefullpath)

date = strmid(fn, 17, 13)

scp_dest = 'soodal@164.125.38.179:/home/soodal/windowshome/works/plot/'

; read GEMS_L2_O3P
data = ds_read_gems_l2_o3p(o3p_file)

; get variable from GEMS_L2_O3P

ak = data.AveragingKernel
akdim = size(ak, /dimension)

pressure = data.Pressure
pbdim = size(pressure, /dimension) ; number of boundary of pressure
; xidx=89, yidx=43

ct = colortable([[255, 0, 0], $
  [127, 127, 0], $
  [0, 255, 0], $
  [0, 127, 127],$
  [0, 0, 255], $
  [0, 127, 127], $
  [255, 0, 255]], $
  /transpose, $
  ncolor=24) 

w = window(/buffer, $
  dimension=[500, 700])
for ilayer=0, akdim[0]-1 do begin
  ak_ = ak[ilayer, *, xidx, yidx]
  pressure_boundary = pressure[*, xidx, yidx] ; boundary
  pressure_center = fltarr(pbdim[0]-1)
  for jp = 0, pbdim[0]-2 do begin
    pressure_center[jp] = (pressure_boundary[jp] + pressure_boundary[jp+1])/2.
  endfor
  dummy = execute($
    " p" + string(ilayer, format='(i02)')+ " = " + $
    "plot(ak_, reverse(pressure_center), " + $
    "/buffer, " + $
    "color=ct[*, ilayer], " + $
    "/ylog, " + $
    "yrange=[1000, 0.1], " + $
    "name='layer'+string(ilayer, format='(i02)'), " + $
    "title='Averaging Kernel', " + $
    "ytitle='Pressure[hPa]', " + $
    "overplot=ilayer) ")

endfor
leg = legend(target=[p00, p01, p02, p03,p04, p05, p06, p07, p08, $
  p09, p10, p11, p12, p13, p14, p15, p16, p17, p18, $
  p19, p20, p21, p22, p23], $
  font_size=10)
p23.save, 'fig3_AveragingKernelfor_x'+$
  strmid(string(xidx), 2) + 'y'+strmid(string(yidx), 2)+'.png'
p23.close
;for ilayer=0, akdim[0]-1 do begin
  ;ak_ = ak[*, ilayer, xidx, yidx]
  ;pressure_boundary = pressure[*, xidx, yidx] ; boundary
  ;pressure_center = fltarr(pbdim[0]-1)
  ;for jp = 0, pbdim[0]-2 do begin
    ;pressure_center[jp] = (pressure_boundary[jp] + pressure_boundary[jp+1])/2.
  ;endfor
  ;p = plot(ak_, reverse(pressure_center), $
    ;/buffer, $
    ;color=ct[*, ilayer], $
    ;/ylog, $
    ;yrange=[1000, 0.1], $
    ;title='Averaging Kernel', $
    ;xtitle='Pressure[hPa]', $
    ;overplot=ilayer)
;endfor
;p.save, 'test2.png'
;p.close
stop



ecf = data.EffectiveCloudFractionUV
ecfnanidx = where(ecf lt -990, /null)
ecf[ecfnanidx] = !values.f_nan

lat = data.latitude
latnanidx = where(lat lt -990, /null)
lat[latnanidx] = !values.f_nan

lon = data.longitude
lonnanidx = where(lon lt -990, /null)
lon[lonnanidx] = !values.f_nan

fitres = data.ResidualsOfFit
nanidx = where(fitres lt -990)
fitres[nanidx] = !values.f_nan
zeroidx = where(fitres eq 0)
fitres[zeroidx] = !values.f_nan
fitresdims = size(fitres, /dim)

simrad = data.SimulatedRadiances
nanidx = where(simrad lt -990)
simrad[nanidx] = !values.f_nan
zeroidx = where(simrad eq 0)
simrad[zeroidx] = !values.f_nan
simraddims = size(simrad, /dim)

wavpath = strmid(fn, 13, 10, /reverse) + '/'
wavpath = ''
wavfn = projectpath + wavpath + 'waves_rad_310_340.txt'
readcol, wavfn, wav

;==============================================================================
; Plot 
;------------------------------------------------------------------------------
; for latitude < 20
;------------------------------------------------------------------------------
latclr_idx = where(lat lt 20 and ecf lt 0.2)
fitres = reform(fitres, fitresdims[0], 174L*512)
simrad = reform(simrad, simraddims[0], 174L*512)
residuals = fltarr(fitresdims[0])
residuals[*] = !values.f_nan
stdres = fltarr(fitresdims[0])
stdres[*] = !values.f_nan

npix = n_elements(latclr_idx)
; pick 10 pixel from 1/20, 3/20, 5/20, ..., 19/20  from the whole idx
selected_idx = (indgen(10) * 2 + 1) * npix/20

for in = 0, n_elements(selected_idx)-1 do begin 
  ;for iwav = 0, fitresdims[0]-1 do begin
    ;_fitres = reform(fitres[iwav, latclr_idx[in]])
    ;_simrad = reform(simrad[iwav, latclr_idx[in]])

    ;residuals[iwav] = mean(_fitres, /nan)/mean(_simrad, /nan)*100
    ;stdres[iwav] = stddev(_fitres, /nan)/mean(_simrad, /nan)*100
  ;endfor
  residuals = fitres[*, latclr_idx[selected_idx[in]]]$
    /simrad[*, latclr_idx[selected_idx[in]]] * 100

  wav0idx = where(wav eq 0, /null)
  wav[wav0idx] = !values.f_nan
  residuals[wav0idx] = !values.f_nan

  mm = minmax(residuals, /nan)
  yrange = [-0.3, 0.3]
  ;IF mm[1] GE 0.3 THEN BEGIN
    ;yrange[1] = mm[1]
  ;ENDIf
  ;IF mm[0] LE -0.3 THEN BEGIN
    ;yrange[0] = mm[0]
  ;ENDIF
  p1 = plot(wav[0:fitresdims[0]-1],residuals, /buffer, $
    xtitle='Wavelength[nm]', $
    ytitle='Residuals[%]', $
    title='GEMS L2 O3P ' + fit_range + 'nm Fitting Residual', $
    symbol='s', $
    yrange=yrange, $
    xrange=[300, 340], $
    name='Residuals', $
    position=pos)
  p_ = plot(indgen(500), fltarr(500), 'gray', /buffer, /over)
  ;p2 = plot(wav[0:fitresdims[0]-1], stdres, 'r', /buffer, $
    ;/over, name='Stddev')

  t1 = text(0.25, 0.75, 'Effective Cloud Fraction < 0.2', /norm)
  t2 = text(0.25, 0.70, 'Latitude < 20', /norm)
  ;leg = legend(target=[p0, p1, p2, p3, pap], pos=leg_pos)
  pngfile = './plot/fitting_residual_' + date + '_' + fitrange $
    + '_latlt20_' + strtrim(string(latclr_idx[selected_idx[in]]),2) + '.png'
  p1.errorbar_color='indian_red'
  p1.errorbar_capsize=0.1
  p1.title.font_size=16
  p1.save, pngfile
  p1.close
  ; send image to pc
  ;spawn, 'scp -P18742 -p ' + pngfile +  ' ' + scp_dest
endfor


;-----------------------------------------------------------------------------
; for 20 <= latitude < 40
;-----------------------------------------------------------------------------
latclr_idx = where(lat ge 20  and lat lt 40 and ecf lt 0.2)
fitres = reform(fitres, fitresdims[0], 174L*512)
simrad = reform(simrad, simraddims[0], 174L*512)
residuals = fltarr(fitresdims[0])
residuals[*] = !values.f_nan
stdres = fltarr(fitresdims[0])
stdres[*] = !values.f_nan

npix = n_elements(latclr_idx)
; pick 10 pixel from 1/20, 3/20, 5/20, ..., 19/20  from the whole idx
selected_idx = (indgen(10) * 2 + 1) * npix/20

for in = 0, n_elements(selected_idx)-1 do begin 
  ;for i = 0, fitresdims[0]-1 do begin
    ;_fitres = reform(fitres[i, latclr_idx])
    ;_simrad = reform(simrad[i, latclr_idx])
    ;residuals[i] = mean(_fitres, /nan)/mean(_simrad, /nan)*100
    ;stdres[i] = stddev(_fitres, /nan)/mean(_simrad, /nan)*100
  ;endfor
  residuals = fitres[*, latclr_idx[selected_idx[in]]]$
    /simrad[*, latclr_idx[selected_idx[in]]] * 100

  wav0idx = where(wav eq 0, /null)
  wav[wav0idx] = !values.f_nan
  residuals[wav0idx] = !values.f_nan

  mm = minmax(residuals, /nan)
  yrange = [-0.3, 0.3]
  ;IF mm[1] GE 0.3 THEN BEGIN
    ;yrange[1] = mm[1]
  ;ENDIf
  ;IF mm[0] LE -0.3 THEN BEGIN
    ;yrange[0] = mm[0]
  ;ENDIF
  p1 = plot(wav[0:fitresdims[0]-1],residuals, /buffer, $
    xtitle='Wavelength[nm]', $
    ytitle='Residuals[%]', $
    title='GEMS L2 O3P ' + fit_range + 'nm Fitting Residual', $
    symbol='s', $
    yrange=yrange, $
    xrange=[300, 340], $
    name='Residuals', $
    position=pos)
  p_ = plot(indgen(500), fltarr(500), 'gray', /buffer, /over)
  ;p2 = plot(wav[0:fitresdims[0]-1], stdres, 'r', /buffer, $
    ;/over, name='Stddev')

  t1 = text(0.25, 0.75, 'Effective Cloud Fraction < 0.2', /norm)
  t2 = text(0.25, 0.70, '20 < Latitude < 40', /norm)
  ;leg = legend(target=[p0, p1, p2, p3, pap], pos=leg_pos)
  pngfile = './plot/fitting_residual_' + date + '_' + fitrange + $
    '_20gelatlt40_' + strtrim(string(latclr_idx[selected_idx[in]]),2) + '.png'
  p1.errorbar_color='indian_red'
  p1.errorbar_capsize=0.1
  p1.title.font_size=16
  p1.save, pngfile
  p1.close
  ; send image to pc
  ;spawn, 'scp -P18742 -p ' + pngfile +  ' ' + scp_dest
endfor

;------------------------------------------------------------------------------
; for latitude < 20 mae
;------------------------------------------------------------------------------
latclr_idx = where(lat lt 20 and ecf lt 0.2)
fitres = reform(fitres, fitresdims[0], 174L*512)
simrad = reform(simrad, simraddims[0], 174L*512)
residuals = fltarr(fitresdims[0])
residuals[*] = !values.f_nan
stdres = fltarr(fitresdims[0])
stdres[*] = !values.f_nan

npix = n_elements(latclr_idx)
; pick 10 pixel from 1/20, 3/20, 5/20, ..., 19/20  from the whole idx
selected_idx = (indgen(10) * 2 + 1) * npix/20

for in = 0, n_elements(selected_idx)-1 do begin 
  ;for i = 0, fitresdims[0]-1 do begin
    ;_fitres = reform(fitres[i, latclr_idx])
    ;_simrad = reform(simrad[i, latclr_idx])
    ;residuals[i] = mean(abs(_fitres), /nan)/mean(abs(_simrad), /nan)*100
    ;stdres[i] = stddev(abs(_fitres), /nan)/mean(abs(_simrad), /nan)*100
  ;endfor
  residuals = fitres[*, latclr_idx[selected_idx[in]]]$
    /simrad[*, latclr_idx[selected_idx[in]]] * 100

  wav0idx = where(wav eq 0, /null)
  wav[wav0idx] = !values.f_nan
  residuals[wav0idx] = !values.f_nan

  mm = minmax(residuals, /nan)
  yrange = [0, 0.3]
  ;IF mm[1] GE 0.3 THEN BEGIN
    ;yrange[1] = mm[1]
  ;ENDIf
  ;IF mm[0] LE -0.3 THEN BEGIN
    ;yrange[0] = mm[0]
  ;ENDIF

  p1 = plot(wav[0:fitresdims[0]-1],residuals, /buffer, $
    xtitle='Wavelength[nm]', $
    ytitle='Residuals[%]', $
    title='GEMS L2 O3P ' + fit_range + 'nm Fitting Residual(MAE)', $
    symbol='s', $
    yrange=yrange, $
    xrange=[300, 340], $
    name='Residuals', $
    position=pos)
  p_ = plot(indgen(500), fltarr(500), 'gray', /buffer, /over)
  ;p2 = plot(wav[0:fitresdims[0]-1], stdres, 'r', /buffer, $
    ;/over, name='Stddev')

  t1 = text(0.25, 0.75, 'Effective Cloud Fraction < 0.2', /norm)
  t2 = text(0.25, 0.70, 'Latitude < 20', /norm)
  ;leg = legend(target=[p0, p1, p2, p3, pap], pos=leg_pos)
  pngfile = './plot/fitting_residual_' + date + '_' + fitrange $
    + '_latlt20_mae_' + strtrim(string(latclr_idx[selected_idx[in]]),2) + '.png'
  p1.errorbar_color='indian_red'
  p1.errorbar_capsize=0.1
  p1.title.font_size=16
  p1.save, pngfile
  p1.close
  ; send image to pc
  ;spawn, 'scp -P18742 -p ' + pngfile +  ' ' + scp_dest
endfor

;-----------------------------------------------------------------------------
; for 20 <= latitude < 40 mae
;-----------------------------------------------------------------------------
latclr_idx = where(lat ge 20  and lat lt 40 and ecf lt 0.2)
fitres = reform(fitres, fitresdims[0], 174L*512)
simrad = reform(simrad, simraddims[0], 174L*512)
residuals = fltarr(fitresdims[0])
residuals[*] = !values.f_nan
stdres = fltarr(fitresdims[0])
stdres[*] = !values.f_nan

npix = n_elements(latclr_idx)
; pick 10 pixel from 1/20, 3/20, 5/20, ..., 19/20  from the whole idx
selected_idx = (indgen(10) * 2 + 1) * npix/20

for in = 0, n_elements(selected_idx)-1 do begin 
  ;for i = 0, fitresdims[0]-1 do begin
    ;_fitres = reform(fitres[i, latclr_idx])
    ;_simrad = reform(simrad[i, latclr_idx])
    ;residuals[i] = mean(abs(_fitres), /nan)/mean(abs(_simrad), /nan)*100
    ;stdres[i] = stddev(_fitres, /nan)/mean(_simrad, /nan)*100
  ;endfor
  residuals = fitres[*, latclr_idx[selected_idx[in]]]$
    /simrad[*, latclr_idx[selected_idx[in]]] * 100

  wav0idx = where(wav eq 0, /null)
  wav[wav0idx] = !values.f_nan
  residuals[wav0idx] = !values.f_nan

  mm = minmax(residuals, /nan)
  yrange = [0, 0.3]
  ;IF mm[1] GE 0.3 THEN BEGIN
    ;yrange[1] = mm[1]
  ;ENDIf
  ;IF mm[0] LE -0.3 THEN BEGIN
    ;yrange[0] = mm[0]
  ;ENDIF
  p1 = plot(wav[0:fitresdims[0]-1],residuals, /buffer, $
    xtitle='Wavelength[nm]', $
    ytitle='Residuals[%]', $
    title='GEMS L2 O3P ' + fit_range + 'nm Fitting Residual(MAE)', $
    symbol='s', $
    yrange=yrange, $
    xrange=[300, 340], $
    name='Residuals', $
    position=pos)
  p_ = plot(indgen(500), fltarr(500), 'gray', /buffer, /over)
  ;p2 = plot(wav[0:fitresdims[0]-1], 'r', /buffer, $
    ;/over, name='Stddev')

  t1 = text(0.25, 0.75, 'Effective Cloud Fraction < 0.2', /norm)
  t2 = text(0.25, 0.70, '20 < Latitude < 40', /norm)
  ;leg = legend(target=[p0, p1, p2, p3, pap], pos=leg_pos)
  pngfile = './plot/fitting_residual_' + date + '_' + fitrange $
    + '_20gelatlt40_mae_' + strtrim(string(latclr_idx[selected_idx[in]]),2) + '.png'
  p1.errorbar_color='indian_red'
  p1.errorbar_capsize=0.1
  p1.title.font_size=16
  p1.save, pngfile
  p1.close
  ; send image to pc
  ;spawn, 'scp -P18742 -p ' + pngfile +  ' ' + scp_dest
endfor
  
;------------------------------------------------------------------------------
; for latitude < 20 rmse
;------------------------------------------------------------------------------
latclr_idx = where(lat lt 20 and ecf lt 0.2)
fitres = reform(fitres, fitresdims[0], 174L*512)
simrad = reform(simrad, simraddims[0], 174L*512)
residuals = fltarr(fitresdims[0])
residuals[*] = !values.f_nan
stdres = fltarr(fitresdims[0])
stdres[*] = !values.f_nan

npix = n_elements(latclr_idx)
; pick 10 pixel from 1/20, 3/20, 5/20, ..., 19/20  from the whole idx
selected_idx = (indgen(10) * 2 + 1) * npix/20

for in = 0, n_elements(selected_idx)-1 do begin 
  ;for i = 0, fitresdims[0]-1 do begin
    ;_fitres = reform(fitres[i, latclr_idx])
    ;_simrad = reform(simrad[i, latclr_idx])
    ;residuals[i] = sqrt(total(_fitres^2, 0, /nan)/total(finite(_fitres)))$
      ;/mean(_simrad)*100
    ;stdres[i] = stddev(_fitres, /nan)/mean(_simrad, /nan)*100
  ;endfor
  residuals = fitres[*, latclr_idx[selected_idx[in]]]$
    /simrad[*, latclr_idx[selected_idx[in]]] * 100

  wav0idx = where(wav eq 0, /null)
  wav[wav0idx] = !values.f_nan
  residuals[wav0idx] = !values.f_nan

  mm = minmax(residuals, /nan)
  yrange = [0, 0.3]
  ;IF mm[1] GE 0.3 THEN BEGIN
    ;yrange[1] = mm[1]
  ;ENDIf
  ;IF mm[0] LE -0.3 THEN BEGIN
    ;yrange[0] = mm[0]
  ;ENDIF
  p1 = plot(wav[0:fitresdims[0]-1],residuals, /buffer, $
    xtitle='Wavelength[nm]', $
    ytitle='Residuals[%]', $
    title='GEMS L2 O3P ' + fit_range + 'nm Fitting Residual(RMSE)', $
    symbol='s', $
    yrange=yrange, $
    xrange=[300, 340], $
    name='Residuals', $
    position=pos)
  p_ = plot(indgen(500), fltarr(500), 'gray', /buffer, /over)
  ;p2 = plot(wav[0:fitresdims[0]-1], stdres, 'r', /buffer, $
    ;/over, name='Stddev')

  t1 = text(0.25, 0.75, 'Effective Cloud Fraction < 0.2', /norm)
  t2 = text(0.25, 0.70, 'Latitude < 20', /norm)
  ;leg = legend(target=[p0, p1, p2, p3, pap], pos=leg_pos)
  pngfile = './plot/fitting_residual_' + date + '_' + fitrange $
    + '_latlt20_rmse_' + strtrim(string(latclr_idx[selected_idx[in]]),2) + '.png'
  p1.errorbar_color='indian_red'
  p1.errorbar_capsize=0.1
  p1.title.font_size=16
  p1.save, pngfile
  p1.close
  ; send image to pc
  ;spawn, 'scp -P18742 -p ' + pngfile +  ' ' + scp_dest
endfor

;-----------------------------------------------------------------------------
; for 20 <= latitude < 40 rmse
;-----------------------------------------------------------------------------
latclr_idx = where(lat ge 20  and lat lt 40 and ecf lt 0.2)
fitres = reform(fitres, fitresdims[0], 174L*512)
simrad = reform(simrad, simraddims[0], 174L*512)
residuals = fltarr(fitresdims[0])
residuals[*] = !values.f_nan
stdres = fltarr(fitresdims[0])
stdres[*] = !values.f_nan

npix = n_elements(latclr_idx)
; pick 10 pixel from 1/20, 3/20, 5/20, ..., 19/20  from the whole idx
selected_idx = (indgen(10) * 2 + 1) * npix/20

for in = 0, n_elements(selected_idx)-1 do begin 
  ;for i = 0, fitresdims[0]-1 do begin
    ;_fitres = reform(fitres[i, latclr_idx])
    ;_simrad = reform(simrad[i, latclr_idx])
    ;residuals[i] = sqrt(total(_fitres^2, 0, /nan)/total(finite(_fitres)))$
      ;/mean(_simrad)*100
    ;stdres[i] = stddev(_fitres, /nan)/mean(_simrad, /nan)*100
  ;endfor
  residuals = fitres[*, latclr_idx[selected_idx[in]]]$
    /simrad[*, latclr_idx[selected_idx[in]]] * 100

  wav0idx = where(wav eq 0, /null)
  wav[wav0idx] = !values.f_nan
  residuals[wav0idx] = !values.f_nan

  yrange = [0, 0.3]
  ;IF mm[1] GE 0.3 THEN BEGIN
    ;yrange[1] = mm[1]
  ;ENDIf
  ;IF mm[0] LE -0.3 THEN BEGIN
    ;yrange[0] = mm[0]
  ;ENDIF
  p1 = plot(wav[0:fitresdims[0]-1],residuals, /buffer, $
    xtitle='Wavelength[nm]', $
    ytitle='Residuals[%]', $
    title='GEMS L2 O3P ' + fit_range + 'nm Fitting Residual(RMSE)', $
    symbol='s', $
    yrange=yrange, $
    xrange=[300, 340], $
    name='Residuals', $
    position=pos)
  p_ = plot(indgen(500), fltarr(500), 'gray', /buffer, /over)
  ;p2 = plot(wav[0:fitresdims[0]-1], stdres, 'r', /buffer, $
    ;/over, name='Stddev')

  t1 = text(0.25, 0.75, 'Effective Cloud Fraction < 0.2', /norm)
  t2 = text(0.25, 0.70, '20 < Latitude < 40', /norm)
  ;leg = legend(target=[p0, p1, p2, p3, pap], pos=leg_pos)
  pngfile = './plot/fitting_residual_' + date + '_' + fitrange $
    + '_20gelatlt40_rmse_' + strtrim(string(latclr_idx[selected_idx[in]]),2) + '.png'
  p1.errorbar_color='indian_red'
  p1.errorbar_capsize=0.1
  p1.title.font_size=16
  p1.save, pngfile
  p1.close
  ; send image to pc
  ;spawn, 'scp -P18742 -p ' + pngfile +  ' ' + scp_dest
endfor

print, 'Done.'
end
