;------------------------------------------------------------------------------
; Configuration
;------------------------------------------------------------------------------
pos = [0.12, 0.1, 0.9, 0.85]
leg_pos = [0.85, 0.8]

fitrange = '305340'
fit_range = strmid(fitrange, 0, 3) + '-' + strmid(fitrange, 3, 3)
L1Cmaker = 'EOSRL'

path ='/data1/gems/o3p/ds/GEMS_O3P_Yonsei/' 
outputpath = path + 'out/'
projectpath = outputpath + 'softcal/' + fitrange + '/'

filelist = file_search(projectpath + '*.nc4')

runtime = strmid(filelist, 13, 10, /reverse)
runtimeidx = sort(runtime)
recentrunfile = filelist[runtimeidx[-1]]

fn = file_basename(recentrunfile)

date = strmid(fn, 13, 13)

scp_dest = 'soodal@164.125.38.179:/home/soodal/WindowsHome/Pictures/scp/softcal'

; read GEMS_L2_O3P
data = ds_read_gems_l2_o3p(projectpath + fn)

; get variable from GEMS_L2_O3P
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
wavfn = projectpath + wavpath + 'waves_rad_.txt'
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
    + '_latlt20_' + strtrim(string(latclr_idx[selected_idx[in]]),2) $
    + '_softcal.png'
  p1.errorbar_color='indian_red'
  p1.errorbar_capsize=0.1
  p1.title.font_size=16
  p1.save, pngfile
  p1.close
  ; send image to pc
  spawn, 'scp -P18742 -p ' + pngfile +  ' ' + scp_dest
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
  pngfile = './plot/fitting_residual_' + date + '_' + fitrange $
    + '_20gelatlt40_' + strtrim(string(latclr_idx[selected_idx[in]]),2) $
    + '_softcal.png'
  p1.errorbar_color='indian_red'
  p1.errorbar_capsize=0.1
  p1.title.font_size=16
  p1.save, pngfile
  p1.close
  ; send image to pc
  spawn, 'scp -P18742 -p ' + pngfile +  ' ' + scp_dest
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
    + '_latlt20_mae_' + strtrim(string(latclr_idx[selected_idx[in]]),2) $
    + '_softcal.png'
  p1.errorbar_color='indian_red'
  p1.errorbar_capsize=0.1
  p1.title.font_size=16
  p1.save, pngfile
  p1.close
  ; send image to pc
  spawn, 'scp -P18742 -p ' + pngfile +  ' ' + scp_dest
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
    + '_20gelatlt40_mae_' + strtrim(string(latclr_idx[selected_idx[in]]),2) $
    + '_softcal.png'
  p1.errorbar_color='indian_red'
  p1.errorbar_capsize=0.1
  p1.title.font_size=16
  p1.save, pngfile
  p1.close
  ; send image to pc
  spawn, 'scp -P18742 -p ' + pngfile +  ' ' + scp_dest
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
    + '_latlt20_rmse_' + strtrim(string(latclr_idx[selected_idx[in]]),2) $
    + '_softcal.png'
  p1.errorbar_color='indian_red'
  p1.errorbar_capsize=0.1
  p1.title.font_size=16
  p1.save, pngfile
  p1.close
  ; send image to pc
  spawn, 'scp -P18742 -p ' + pngfile +  ' ' + scp_dest
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
    + '_20gelatlt40_rmse_' + strtrim(string(latclr_idx[selected_idx[in]]),2) $
    + '_softcal.png'
  p1.errorbar_color='indian_red'
  p1.errorbar_capsize=0.1
  p1.title.font_size=16
  p1.save, pngfile
  p1.close
  ; send image to pc
  spawn, 'scp -P18742 -p ' + pngfile +  ' ' + scp_dest
endfor

print, 'Done.'
end
