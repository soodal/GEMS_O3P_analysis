pro plot_fitres, filename, xidx, yidx

sXidx = string(xidx, format='(i03)')
sYidx = string(yidx, format='(i03)')

;------------------------------------------------------------------------------
; Configuration
;------------------------------------------------------------------------------
pos = [0.12, 0.1, 0.9, 0.85]
leg_pos = [0.85, 0.8]

winliminitpos = strpos(filename, 'winliminit')
fitrange = strmid(filename, winliminitpos+10, 3) + '340'
fit_range = strmid(fitrange, 0, 3) + '-' + strmid(fitrange, 3, 3)
L1Cmaker = 'EOSRL'
stop

path ='/data1/gems/o3p/ds/GEMS_O3P_Yonsei/' 
outputpath = path + 'out/'
projectpath = outputpath ;+ 'softcal/' + fitrange + '/'

filelist = file_search(projectpath + '*.nc4')


runtime = strmid(filelist, 13, 10, /reverse)
runtimeidx = sort(runtime)
recentrunfile = filelist[runtimeidx[-1]]

;filefullpath = '/data1/gems/o3p/ds/GEMS_O3P_Yonsei/out/GK2B_GEMS_L2_O3P_20200616_0345_climML_winliminit310_2020-12-05T0626.nc4'

;snrnc = '/data1/gems/o3p/ds/GEMS_O3P_Yonsei/out/GK2B_GEMS_L2_O3P_20200616_0345_snr__climML_winliminit310_2020-12-16T1724.nc4


;fn = file_basename(recentrunfile)
;fn = file_basename(filefullpath)

datepos = stregex(filename,'[0-9]{8}_[0-9]{4}', length=len)
date = strmid(filename, datepos, len)

; xidx=89, yidx=43

;scp_dest = 'soodal@164.125.38.179:/home/soodal/windowshome/works/plot/'

; read GEMS_L2_O3P
;data = ds_read_gems_l2_o3p(filename)
data = ds_read_gems_l2_o3p(filename, $
  varlist=['EffectiveCloudFractionUV', $
  'ResidualsOfFit', 'Latitude', 'Longitude', $
  'SimulatedRadiances', 'FitWeights', 'SignalToNoiseRatio', 'WavelengthsWholeRange'])
basename = file_basename(filename)

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

;wavpath = strmid(fn, 13, 10, /reverse) + '/'
;wavpath = ''
;wavfn = projectpath + wavpath + 'waves_rad_310_340.txt'
;readcol, wavfn, wav
wav = data.WavelengthsWholeRange

;==============================================================================
; Plot 
;------------------------------------------------------------------------------
; for latitude < 20
;------------------------------------------------------------------------------
latclr_idx = where(lat lt 20 and ecf lt 0.2)
;fitres = reform(fitres, fitresdims[0], 174L*512)
;simrad = reform(simrad, simraddims[0], 174L*512)
residuals = fltarr(fitresdims[0])
residuals[*] = !values.f_nan
stdres = fltarr(fitresdims[0])
stdres[*] = !values.f_nan

npix = n_elements(latclr_idx)
; pick 10 pixel from 1/20, 3/20, 5/20, ..., 19/20  from the whole idx
;selected_idx = (indgen(10) * 2 + 1) * npix/20

;for in = 0, n_elements(selected_idx)-1 do begin 
  ;for iwav = 0, fitresdims[0]-1 do begin
    ;_fitres = reform(fitres[iwav, latclr_idx[in]])
    ;_simrad = reform(simrad[iwav, latclr_idx[in]])

    ;residuals[iwav] = mean(_fitres, /nan)/mean(_simrad, /nan)*100
    ;stdres[iwav] = stddev(_fitres, /nan)/mean(_simrad, /nan)*100
  ;endfor
  residuals = fitres[*, xidx, yidx]$
    /simrad[*, xidx, yidx] * 100

  wav0idx = where(wav eq 0, /null)
  wav[wav0idx] = !values.f_nan
  residuals[wav0idx] = !values.f_nan

  mm = minmax(residuals, /nan)
  yrange = [-0.5, 0.5]
  ;IF mm[1] GE 0.3 THEN BEGIN
    ;yrange[1] = mm[1]
  ;ENDIf
  ;IF mm[0] LE -0.3 THEN BEGIN
    ;yrange[0] = mm[0]
  ;ENDIF
  p1 = plot(wav[0:fitresdims[0]-1, xidx, yidx],residuals, /buffer, $
    axis_style=1, $
    xtitle='Wavelength[nm]', $
    ytitle='Residuals[%]', $
    title='GEMS L2 O3P ' + fit_range + 'nm Fitting Residual', $
    ;symbol='s', $
    yrange=yrange, $
    xrange=[300, 340], $
    name='Residuals', $
    color='red', $
    position=pos)
  p1['axis1'].color='red'

  ;plot simulated radiance
  p2 = plot($
    wav[0:fitresdims[0]-1, xidx, yidx], simrad[0:fitresdims[0]-1, xidx, yidx], $
    /buffer, $
    axis_style=0, $
    name='Simulated', $
    ;symbol='s', $
    ;yrange=yrange, $
    xrange=[300, 340], $
    position=pos, $
    color='blue', $
    /current)
  yaxis2 = axis('Y', LOCATION='right', yrange=p2.yrange, TARGET=p2, $
    color='blue')
  xaxis2 = axis('X', LOCATION='TOP', xrange=p1.xrange, TARGET=p2, $
    color='black', $
    tickname=['', '', '', '', ''])

  ;plot observed radiance
  p3 = plot($
    wav[0:fitresdims[0]-1, xidx, yidx], $
    simrad[0:fitresdims[0]-1, xidx, yidx] + fitres[0:fitresdims[0]-1, xidx, yidx], $
    /buffer, $
    axis_style=0, $
    name='Observed', $
    ;symbol='s', $
    yrange=p2.yrange, $
    xrange=[300, 340], $
    position=pos, $
    color='black', $
    /current)
  ;yaxis2 = axis('Y', LOCATION='right', yrange=p3.yrange, TARGET=p3, $
    ;color='red')

  p_ = plot(indgen(500), fltarr(500), 'gray', /buffer, /current, $
    axis_style=0, $
    position=pos, $
    yrange=p1.yrange, $
    xrange=p1.xrange, $
    color='gray')
  leg = legend(target=[p1, p2, p3], position=[318, 0.4], $
    /data, /auto_text_color)
  ;p3 = plot(wav[0:fitresdims[0]-1], stdres, 'r', /buffer, $
    ;/over, name='Stddev')

  ;t1 = text(0.25, 0.75, 'Effective Cloud Fraction < 0.2', /norm)
  ;t2 = text(0.25, 0.70, 'Latitude < 20', /norm)
  ;leg = legend(target=[p0, p1, p2, p3, pap], pos=leg_pos)
  pngfile = './plot/' +basename + '_fitres_x' +sXidx +'_y' +sYidx +'.png'
  p1.errorbar_color='indian_red'
  p1.errorbar_capsize=0.1
  p1.title.font_size=16
  p1.save, pngfile
  p1.close
  ; send image to pc
  scp_dest = 'soodal@164.125.38.179:/home/soodal/works/plot/'
  spawn, 'scp -P18742 -p ' + pngfile +  ' ' + scp_dest
;endfor

print, 'Done.'
end
