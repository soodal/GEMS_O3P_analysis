savfile = './test0043_variables.sav'
if not file_test(savfile) then begin

  m_param = ['Met_PMID', 'Met_U', 'Met_V', 'Met_T']
  g_param = ['SpeciesConc_O3']

  homepath = '/data1/syndata_20190326/'
  synt_fn = 'GEMS_2016011503_Trueinput.nc'

  synt = ds_read_synt_model(homepath+synt_fn)

  synt_gems_o3t_path = '/data1/syndata_20190326/L2/'
  synt_gems_o3t_fn = 'GEMS_SYN_L2_O3T_2016011503_Cloud_gaussLUT.nc'

  ;synt_alt_true_4x4 = fltarr(176, 512, 56);congrid(synt.altitude, 176, 512, 55)
  ;synt_pres_true_4x4 = fltarr(176, 512, 56);congrid(synt.pressure, 176, 512, 55)
  ;synt_temp_4x4 = fltarr(176,512,55);congrid(synt.temperature, 176, 512, 55)

  ;synt_o3_true_4x4 = fltarr(176, 512, 55)
  ;for ix=0, 175 do begin
    ;for iy = 0, 511 do begin
      ;for iz = 0, 54 do begin 
        ;if ix ne 175 then begin
          ;synt_o3_true_4x4[ix, iy, iz] = mean(synt.o3[$
            ;ix*4:ix*4+3, iy*4:iy*4+3, iz])
        ;endif else begin
          ;synt_o3_true_4x4[ix, iy, iz] = mean(synt.o3[$
            ;ix*4, iy*4:iy*4+3, iz])
        ;endelse
      ;endfor
    ;endfor
  ;endfor
  synt_o3_true_4x4 = synt.o3

  ;synt_alt_true_4x4 = fltarr(176, 512, 56)
  ;for ix=0, 175 do begin
    ;for iy = 0, 511 do begin
      ;for iz = 0, 55 do begin 
        ;if ix ne 175 then begin
          ;synt_alt_true_4x4[ix, iy, iz] = mean(synt.altitude[$
            ;ix*4:ix*4+3, iy*4:iy*4+3, iz])
        ;endif else begin
          ;synt_alt_true_4x4[ix, iy, iz] = mean(synt.altitude[$
            ;ix*4, iy*4:iy*4+3, iz])
        ;endelse
      ;endfor
    ;endfor
  ;endfor
  synt_alt_true_4x4 = synt.altitude

  ;synt_pres_true_4x4 = fltarr(176, 512, 56)
  ;for ix=0, 175 do begin
    ;for iy = 0, 511 do begin
      ;for iz = 0, 55 do begin 
        ;if ix ne 175 then begin
          ;synt_pres_true_4x4[ix, iy, iz] = mean(synt.pressure[$
            ;ix*4, iy*4:iy*4+3, iz])
        ;endif else begin 
          ;synt_pres_true_4x4[ix, iy, iz] = mean(synt.pressure[$
            ;ix*4, iy*4:iy*4+3, iz])
        ;endelse
      ;endfor
    ;endfor
  ;endfor
  synt_pres_true_4x4 = synt.pressure
  save, filename=savfile
endif else begin
  restore, filename=savfile
endelse

varlist = [ $;'EffectiveCloudFractionUV', 'ProcessingQualityFlags', $
           ;'RootMeanSquareErrorOfFit', 'AveragingKernel', $
           'ColumnAmountO3' , $
           ;'O3Apriori', 'O3AprioriError','ColumnAmountO3', $
           ;'CloudPressure','ResidualsOfFit','NumberOfIterations','Runtime', $
           ;'SimulatedRadiances', $
           'Latitude' ,'Longitude'];, $
           ;'SolarZenithAngle', $
           ;'ViewingZenithAngle', 'Time','Altitude' ,    $
           ;'Pressure','Pix','Line','TropopausePressure', 'Temperature', $
           ;'TerrainPressure', 'RelativeAzimuthAngle', $
           ;'FitWeights', 'SignalToNoiseRatio', 'Wavelengths', $
           ;'WavelengthsWholeRange']
synt_gems_o3t = ds_read_gems_l2_o3t(synt_gems_o3t_path + synt_gems_o3t_fn, varlist=varlist)
;synt_gems_o3t_lon_4x4 = fltarr(176, 512)
;for ix=0, 175 do begin
  ;for iy = 0, 511 do begin
    ;if ix ne 175 then begin
      ;synt_gems_o3t_lon_4x4[ix, iy] = mean(synt_gems_o3t.longitude[$
        ;ix*4:ix*4+3, iy*4:iy*4+3])
    ;endif else begin
      ;synt_gems_o3t_lon_4x4[ix, iy] = mean(synt_gems_o3t.longitude[$
        ;ix*4, iy*4:iy*4+3])
    ;endelse
  ;endfor
;endfor
synt_gems_o3t_lon_4x4 = synt_gems_o3t.longitude

;synt_gems_o3t_lat_4x4 = fltarr(176, 512)
;for ix=0, 175 do begin
  ;for iy = 0, 511 do begin
    ;if ix ne 175 then begin
      ;synt_gems_o3t_lat_4x4[ix, iy] = mean(synt_gems_o3t.latitude[$
        ;ix*4:ix*4+3, iy*4:iy*4+3])
    ;endif else begin
      ;synt_gems_o3t_lat_4x4[ix, iy] = mean(synt_gems_o3t.latitude[$
        ;ix*4, iy*4:iy*4+3])
    ;endelse
  ;endfor
;endfor
synt_gems_o3t_lat_4x4 = synt_gems_o3t.latitude

;synt_gems_o3t_cao3_4x4 = fltarr(176, 512)
;for ix=0, 175 do begin
  ;for iy = 0, 511 do begin
    ;if ix ne 175 then begin
      ;synt_gems_o3t_cao3_4x4[ix, iy] = mean(synt_gems_o3t.ColumnAmountO3[$
        ;ix*4:ix*4+3, iy*4:iy*4+3])
    ;endif else begin
      ;synt_gems_o3t_cao3_4x4[ix, iy] = mean(synt_gems_o3t.ColumnAmountO3[$
        ;ix*4, iy*4:iy*4+3])
    ;endelse
  ;endfor
;endfor
synt_gems_o3t_cao3_4x4 = synt_gems_o3t.ColumnAmountO3


n_air = 2.69E19
synt_o3_true_4x4_du = fltarr(176, 512, 55)

;for ix = 0, 175 do begin
  ;for iy = 0, 511 do begin
    for iz = 0, 54 do begin
      synt_o3_true_4x4_du[*, *, iz] = synt_o3_true_4x4[*, *, iz]*$
        (synt_alt_true_4x4[*, *, iz] - $
        synt_alt_true_4x4[*, *, iz+1])*1000.*1000.*100./n_air/200000.
    endfor
  ;endfor
;endfor


o3_true_du = synt_o3_true_4x4/n_air*1000

synt_o3_true_4x4_du = o3_true_du


nanidx = where(synt_o3_true_4x4 lt -990, /null)
synt_o3_true_4x4_du[nanidx] = !values.f_nan

synt_o3_true_4x4_new = create_struct('O3', synt_o3_true_4x4_du)
synt_o3_true_4x4_new = create_struct($
  synt_o3_true_4x4_new, 'Altitude', synt_alt_true_4x4)
synt_o3_true_4x4_new = create_struct($
  synt_o3_true_4x4_new, 'Pressure', synt_pres_true_4x4)

ds_model_syntest_accum, synt_o3_true_4x4_new, $
  synt_o3_true_4x4_accum, height=10

;m = map('Robinson', /buffer)
c = contour($
  synt_o3_true_4x4_accum, synt_gems_o3t_lon_4x4, synt_gems_o3t_lat_4x4, $
  /fill, $
  ;overplot=m, $
  n_levels=50, /buffer)
cb = colorbar()
pngfile = 'synt_o3_tropo.png'
c.save, pngfile
c.close

scp_dest = 'soodal@164.125.38.179:/home/soodal/works/plot'
spawn, 'scp -P18742 -p ' + pngfile + $
  ' ' + scp_dest

;m = map('Robinson', /buffer)
c = contour(total(synt_o3_true_4x4_du, 3), synt_gems_o3t_lon_4x4, synt_gems_o3t_lat_4x4, /fill, $
  ;overplot=m, $
  n_levels=50, /buffer)
cb = colorbar()
pngfile = 'synt_o3_true_toz.png'
c.save, pngfile
c.close

scp_dest = 'soodal@164.125.38.179:/home/soodal/works/plot'
spawn, 'scp -P18742 -p ' + pngfile + $
  ' ' + scp_dest

;m = map('Robinson', /buffer)
c = contour(synt_gems_o3t_cao3_4x4, synt_gems_o3t_lon_4x4, synt_gems_o3t_lat_4x4, /fill, $
  ;overplot=m, $
  n_levels=50, /buffer)
cb = colorbar()
pngfile = 'synt_o3t_cao3.png'
c.save, pngfile
c.close

scp_dest = 'soodal@164.125.38.179:/home/soodal/works/plot'
spawn, 'scp -P18742 -p ' + pngfile + $
  ' ' + scp_dest

; scatterplot x: SYNT_GEMS_O3T_CAO3 y: SYNT_GEMS_O3P_CAO3
plot_gems_validation, synt_gems_o3t_cao3_4x4, total(synt_o3_true_4x4_du,3), $
    filename='./plot/gems_synt_o3p_true_toz_synt_o3t_cao3_validation_' + $
    ;strmid(gemsyyyymmdd[0], 0, 6) +'_TOZ_ecf02.png', $
    '20160115t0300utc.png', $
    cblim=[0, 100], $
    range=[250, 450], $
    scp_send=1, $
    delta=2.0

o3size = size(synt_o3_true_4x4, /dim)

;;----------------------------------------------------------------
;; GEMS_O3P_from SYNTHETIC
;;----------------------------------------------------------------

varlist = ['EffectiveCloudFractionUV', 'ProcessingQualityFlags', $
       ;'RootMeanSquareErrorOfFit', 'AveragingKernel', $
         'O3' , $
           ;'O3Apriori', 'O3AprioriError','ColumnAmountO3', $
           ;'CloudPressure','ResidualsOfFit','NumberOfIterations','Runtime', $
           ;'SimulatedRadiances', $
           'Latitude' ,'Longitude', $
           ;'SolarZenithAngle', $
           'ViewingZenithAngle', 'Time','Altitude' ,    $
           ;'Pressure','Pix','Line','TropopausePressure', 'Temperature', $
           ;'TerrainPressure', 'RelativeAzimuthAngle', $
           ;'FitWeights', 'SignalToNoiseRatio', 'Wavelengths', $
           'WavelengthsWholeRange']

;synt_gems_o3p_path = '/data1/gems/o3p/GEMS/o3p/dat/out/'
;synt_gems_o3p_fn = 'GEM_TEST_L2_O3P_20160115_0300_v201907_radiance_cloud.nc'
;synt_gems_o3p = ds_read_gems_l2_o3p(synt_gems_o3p_path+synt_gems_o3p_fn)

outfile = './plot/synt_model_20160616t0345utc.png'
plot_gems_satellite, synt_o3_true_4x4_accum, $
  synt_gems_o3t.longitude, synt_gems_o3t.latitude, $
  ;title='GEMS L2 O3P for ' + mondate+'T'+utchhmm, $
  title=' ', $
  ;range=[20, 60], $
  outfile=outfile, $
  /scp_send
  ;scp_send=0
  
outfile = './plot/synt_model_20160616t0345utc_latitude.png'
plot_gems_satellite, synt_gems_o3t.latitude, $
  synt_gems_o3t.longitude, synt_gems_o3t.latitude, $
  ;title='GEMS L2 O3P for ' + mondate+'T'+utchhmm, $
  title=' ', $
  ;range=[20, 60], $
  outfile=outfile, $
  /scp_send
  ;scp_send=0
  
outfile = './plot/synt_model_20160616t0345utc_longitude.png'
plot_gems_satellite, synt_gems_o3t.longitude, $
  synt_gems_o3t.longitude, synt_gems_o3t.latitude, $
  ;title='GEMS L2 O3P for ' + mondate+'T'+utchhmm, $
  title=' ', $
  ;range=[20, 60], $
  outfile=outfile, $
  /scp_send
  ;scp_send=0

;;--------------------------------------------------
;; Make grid for OMI 0.5 grid
;;--------------------------------------------------

Slat = -5.
Nlat = 45.
Llon = 75.
Rlon = 145.

;; * MAKE grid

res = 0.5
nx = (Rlon-Llon)/res+1.
ny = (Nlat-Slat)/res+1.
LAT = fltarr(nx,ny)
LON = fltarr(nx,ny)
for j = 0, ny-1 do lat(*,j) = Nlat-res*j
for i = 0, nx-1 do lon(i,*) = Llon+res*i

;;--------------------------------------------------

;;----------------------------------------
;; Synthetic O3P True input for 20160616
;;----------------------------------------

outputpath = './plot/'


synt_o3p_stacking = replicate(0., nx, ny)
synt_o3p_stack_num = replicate(0, nx, ny)
synt_o3p_accum_mean = replicate(!values.f_nan, nx, ny)


; synt not have lon lat.
synt_roiidx = where( (synt_gems_o3t_lat_4x4 ge Slat) and $
  (synt_gems_o3t_lat_4x4 le Nlat) and $
  (synt_gems_o3t_lon_4x4 ge Llon) and $
  (synt_gems_o3t_lon_4x4 le Rlon), roipixnum, /null)


;if roipixnum gt 0 then begin
  ;for jj = 0, ny-1 do begin
  ;yidx = where(abs(synt_gems_o3t_lat_4x4-(Nlat-jj*res)) lt res*0.5,ynum)
    ;if ynum gt 0 then begin
    ;for ii = 0, nx-1 do begin
      ;xidx = where(abs(synt_gems_o3t_lon_4x4[yidx]-(Llon+ii*res)) lt res*0.5,xnum)
      ;if xnum gt 0 then begin

        ;data_o3p = total(synt_o3_true_4x4_accum[yidx[xidx]], /nan)
        ;pixnum = total(finite(synt_o3_true_4x4_accum[yidx[xidx]]))
        ;synt_o3p_stacking[ii, jj] = synt_o3p_stacking[ii, jj] + data_o3p
        ;synt_o3p_stack_num[ii, jj] = synt_o3p_stack_num[ii, jj] + pixnum
      ;endif
    ;endfor
    ;endif
  ;endfor
;endif

;synt_o3p_accum_mean = synt_o3p_stacking/synt_o3p_stack_num
;numzeroidx = where(synt_o3p_stack_num eq 0, /null)
;;notzeroidx = where(synt_o3p_stack_num gt 0, /null)
;synt_o3p_accum_mean[numzeroidx] = !values.f_nan


;;----------------------------------------------------------------
;; GEMS_O3P_from SYNTHETIC
;;----------------------------------------------------------------

;synt_gems_o3p_path = '/data1/gems/o3p/GEMS/o3p/dat/out/'
;synt_gems_o3p_fn = 'GEM_TEST_L2_O3P_20160115_0300_v201903_radiance_clear_test_update.nc'

varlist = ['EffectiveCloudFractionUV', 'ProcessingQualityFlags', $
           'O3' , $
           ;'O3Apriori', 'O3AprioriError', $
           ;'CloudPressure', $
           ;'SimulatedRadiances', $
           ;'O3Apriori', 'O3AprioriError',$
           'ColumnAmountO3', $
           'Latitude' ,'Longitude', $
           'Time','Altitude' ,    $
           'Pressure', $
           ;'TropopausePressure', $
           ;'Wavelengths', $
           'WavelengthsWholeRange']


synt_gems_o3p_path = '/data1/gems/o3p/works/'
synt_gems_o3p_fn = 'GK2B_GEMS_L2_20160115_0300_O3P_flipped_bugfix2.4x4.nc4'
synt_gems_o3p = ds_read_gems_l2_o3p(synt_gems_o3p_path + synt_gems_o3p_fn, varlist=varlist)

ds_gems_l2o3p_accum, synt_gems_o3p, synt_gems_o3p_under10km, height=10.

;synt_gems_o3p_accum_stacking = replicate(0., nx, ny)
;synt_gems_o3p_accum_stack_num = replicate(0, nx, ny)
;synt_gems_o3p_accum_mean = replicate(!values.f_nan, nx, ny)


;synt_gems_o3p_roiidx = where( (synt_gems_o3t_lat_4x4 ge Slat) and $
  ;(synt_gems_o3t_lat_4x4 le Nlat) and $
  ;(synt_gems_o3t_lon_4x4 ge Llon) and $
  ;(synt_gems_o3t_lon_4x4 le Rlon), $
  ;synt_gems_o3t_roiidx_num, /null)

;if synt_gems_o3t_roiidx_num gt 0 then begin
  ;for jj = 0, ny-1 do begin
  ;yidx = where(abs(synt_gems_o3t_lat_4x4-(Nlat-jj*res)) lt res*0.5,ynum)
    ;if ynum gt 0 then begin
    ;for ii = 0, nx-1 do begin
      ;xidx = where(abs(synt_gems_o3t_lon_4x4[yidx]-(Llon+ii*res)) lt $
        ;res*0.5,xnum)
      ;if xnum gt 0 then begin
        ;data_o3p = total(synt_gems_o3p_under10km[yidx[xidx]], /nan)
        ;pixnum = total(finite(synt_gems_o3p_under10km[yidx[xidx]]))
        ;synt_gems_o3p_accum_stacking[ii, jj] = $
          ;synt_gems_o3p_accum_stacking[ii, jj] + data_o3p
        ;synt_gems_o3p_accum_stack_num[ii, jj] = $
          ;synt_gems_o3p_accum_stack_num[ii, jj] + pixnum
      ;endif
    ;endfor
    ;endif
  ;endfor
;endif

;synt_gems_o3p_accum_mean = synt_gems_o3p_accum_stacking/synt_gems_o3p_accum_stack_num
;numzeroidx = where(synt_gems_o3p_accum_stack_num eq 0, /null)
;notzeroidx = where(gemsomi_stack_num gt 0, /null)
;synt_gems_o3p_accum_mean[numzeroidx] = !values.f_nan

;valididx = where(finite(synt_gems_o3p_accum_mean) eq 1 and finite(synt_o3p_accum_mean) eq 1, /null)

flipped_synt_o3_true_4x4_accum = reverse(reverse(synt_o3_true_4x4_accum), 2)

valididx = where(finite(synt_gems_o3p_under10km) eq 1 and finite(flipped_synt_o3_true_4x4_accum) eq 1, /null)

plot_gems_validation, synt_gems_o3p_under10km[valididx], flipped_synt_o3_true_4x4_accum[valididx], $
  filename='./plot/gems_synt_o3p_validation_' + $
    ;strmid(gemsyyyymmdd[0], 0, 6) +'_TOZ_ecf02.png', $
    '20160115t0300utc.png', $
  xtitle='SYNTHETIC GEMS O3P Tropospheric Column', $
  ytitle='SYNTHETIC O3P Tropospheric Column', $
  cblim=[0, 2000], $
  range=[0, 60], $
  scp_send=1, $
  delta=2.0

outfile = './plot/synt_gems_o3p_accum_mean_20160616t0345utc.png'
plot_gems_satellite, synt_gems_o3p_under10km, $
  synt_gems_o3p.longitude, synt_gems_o3p.latitude, $
  ;title='GEMS L2 O3P for ' + mondate+'T'+utchhmm, $
  title=' ', $
  range=[20, 60], $
  outfile=outfile, $
  /scp_send
  ;scp_send=0

end
