; not complete

gemsfile1 = '/data2/L2_GEMS/L1C_test/o3p/EOSRL/GK2_GEMS_L2_O3P_20210329_0345_EOSRL_BIN4x4.nc'
gemsfile2 = '/data2/L2_GEMS/L1C_test/o3p/EOSRL/GK2_GEMS_L2_O3P_20210329_0345_EOSRL_BIN4x4.nc'

gems_basename1 = file_basename(gemsfile1)
gems_basename2 = file_basename(gemsfile2)


project = 'EOSRL'

sub = '310_340'
;if not keyword_set(outputpath) then begin
  outputpath = './plot/'
;endif

;if keyword_set(project) then begin
  outputpath = outputpath + project + '/'
;endif else begin
  ;project = ''
;endelse

;if keyword_set(sub) then begin
  outputpath = outputpath + sub + '/'
;endif else begin
  ;sub = ''
;endelse


;==================================
;(1)  read gems l2 o3p file
;==================================  
datetimepos1 = stregex(gemsfile1, '[[:digit:]]{8}_[[:digit:]]{4}')
datetimepos2 = stregex(gemsfile2, '[[:digit:]]{8}_[[:digit:]]{4}')

;print, gemsfile
;print, datetimepos
year1 = fix(strmid(gemsfile1, datetimepos1, 4))
print, strmid(gemsfile1, datetimepos1, 4)
month1 = fix(strmid(gemsfile1, datetimepos1+4, 2))
print, strmid(gemsfile1, datetimepos1+4, 2)
day1 = fix(strmid(gemsfile1, datetimepos1+6, 2))
print, strmid(gemsfile1, datetimepos1+6, 2)
hour1 = fix(strmid(gemsfile1, datetimepos1+9, 2))
print, strmid(gemsfile1, datetimepos1+9, 2)
minute1 = fix(strmid(gemsfile1, datetimepos1+11, 2))
print, strmid(gemsfile1, datetimepos1+11, 2)

year2 = fix(strmid(gemsfile2, datetimepos2, 4))
print, strmid(gemsfile2, datetimepos2, 4)
month2 = fix(strmid(gemsfile2, datetimepos2+4, 2))
print, strmid(gemsfile2, datetimepos2+4, 2)
day2 = fix(strmid(gemsfile2, datetimepos2+6, 2))
print, strmid(gemsfile2, datetimepos2+6, 2)
hour2 = fix(strmid(gemsfile2, datetimepos2+9, 2))
print, strmid(gemsfile2, datetimepos2+9, 2)
minute2 = fix(strmid(gemsfile2, datetimepos2+11, 2))
print, strmid(gemsfile2, datetimepos2+11, 2)

;caldat, jd_list[idate], month, day, year, hour, minute
yyyy1 = string(year1, format='(i04)')
mm1 = string(month1, format='(i02)')
dd1 = string(day1, format='(i02)')
hh1 = string(hour1, format='(i02)')
mi1 = string(minute1, format='(i02)')

datetime_str1= yyyy1 + mm1 + dd1 + '_' + hh1 + mi1

outputpath1 = outputpath1 + '/'  + datetime_str1 + '/point_profile/'

yyyy2 = string(year2, format='(i04)')
mm2 = string(month2, format='(i02)')
dd2 = string(day2, format='(i02)')
hh2 = string(hour2, format='(i02)')
mi2 = string(minute2, format='(i02)')

datetime_str2= yyyy2 + mm2 + dd2 + '_' + hh2 + mi2

outputpath2 = outputpath2 + '/'  + datetime_str2 + '/point_profile/'


if not file_test(outputpath) then begin
  file_mkdir, outputpath 
endif

varlist = ['EffectiveCloudFractionUV', 'ProcessingQualityFlags', $
           'O3' , $
           ;'O3Apriori', 'O3AprioriError', $
           ;'CloudPressure', $
           ;'SimulatedRadiances', $
           ;'O3Apriori', 'O3AprioriError',$
           'ColumnAmountO3', $
           'Latitude' ,'Longitude', $
           'Time','Altitude' ,    $
           ;'Pressure', 'TropopausePressure', $
           ;'Wavelengths', $
           'WavelengthsWholeRange']

gemsvar1 = ds_read_gems_l2_o3p(gemsfile1, varlist=varlist)
gemsvar2 = ds_read_gems_l2_o3p(gemsfile2, varlist=varlist)

cldvar1 = ds_read_gems_l2_cld('/data/private/soodal/ln/GEMS/L2CLD/GK2_GEMS_L2_CLD_' + datetime_str1 + '_4x4.nc')
cldvar2 = ds_read_gems_l2_cld('/data/private/soodal/ln/GEMS/L2CLD/GK2_GEMS_L2_CLD_' + datetime_str2 + '_4x4.nc')

;=============================================================================
; READ MERRA2
;=============================================================================

merra2_1 = ds_read_merra2_tavg3_3d_asm_nv_collocated_on_gems_l1c('/data/private/soodal/MERRA2_collocated_on_gems/merra2_tavg3_3d_asm_collocated_on_gemsl1c_rad_' + datetime_str1 + '.nc4')
merra2_2 = ds_read_merra2_tavg3_3d_asm_nv_collocated_on_gems_l1c('/data/private/soodal/MERRA2_collocated_on_gems/merra2_tavg3_3d_asm_collocated_on_gemsl1c_rad_' + datetime_str2 + '.nc4')


itime = 1
o3mmr1 = reform(merra2_1.o3)
pres1 = reform(merra2_1.PL)
temp1 = reform(merra2_1.T)
alt1 = reform(merra2_1.H)

lon1 = merra2_1.lon
lat1 = merra2_1.lat

itime = 1
o3mmr2 = reform(merra2_2.o3)
pres2 = reform(merra2_2.PL)
temp2 = reform(merra2_2.T)
alt2 = reform(merra2_2.H)

lon2 = merra2_2.lon
lat2 = merra2_2.lat

ds_xymake2d, lon1, lat1, lon2d1, lat2d1
ds_xymake2d, lon2, lat2, lon2d2, lat2d2

;convert mass mxing ratio into vmr
o3vmr1 = ds_mmr2vmr(o3mmr1)
merra2_hpa1 = pres1/100.
merra2_hpa_mid1 = merra2_hpa1[*, *, 1:-1] -merra2_hpa1[*, *, 0:-2]

ds_vmr2du_3d, o3vmr1, merra2_hpa1, temp1, alt1, merra2_o3_1, increasing_height=0
merra2_sz1 = size(merra2_o3_1, /dimension)

merra2_o3_accum1 = fltarr(merra2_sz1)
merra2_o3_accum1[*, *, 0] = merra2_o3_1[*, *, 0]
for i=1, merra2_sz1[2]-1 do begin
  merra2_o3_accum1[*, *, i] = merra2_o3_accum1[*, *, i-1] + merra2_o3_1[*, *, i]
endfor

merra2_o3_on_gems_pres1 = fltarr(size(gemsvar1.pressure, /dimension))
merra2_o3_on_gems_o3_1 = fltarr(size(gemsvar1.o3, /dimension))

for iy=0, merra2_sz1[1]-1 do begin
  for ix=0, merra2_sz1[0]-1 do begin
    merra2_o3_on_gems_pres1[ix, iy, *] = $
      interpol(merra2_o3_accum1[ix, iy, *], merra2_hpa1[ix, iy, 1:-1], gemsvar1.pressure[ix, iy, *])
  endfor
endfor


merra2_o3_on_gems_o3_1 = merra2_o3_on_gems_pres1[*, *, 1:-1] - merra2_o3_on_gems_pres1[*, *, 0:-2]

;convert mass mxing ratio into vmr
o3vmr2 = ds_mmr2vmr(o3mmr2)
merra2_hpa2 = pres2/100.
merra2_hpa_mid2 = merra2_hpa2[*, *, 1:-1] -merra2_hpa2[*, *, 0:-2]

ds_vmr2du_3d, o3vmr2, merra2_hpa2, temp2, alt2, merra2_o3_2, increasing_height=0
merra2_sz2 = size(merra2_o3_2, /dimension)

merra2_o3_accum2 = fltarr(merra2_sz2)
merra2_o3_accum2[*, *, 0] = merra2_o3_2[*, *, 0]
for i=1, merra2_sz2[2]-1 do begin
  merra2_o3_accum2[*, *, i] = merra2_o3_accum2[*, *, i-1] + merra2_o3_2[*, *, i]
endfor

merra2_o3_on_gems_pres2 = fltarr(size(gemsvar2.pressure, /dimension))
merra2_o3_on_gems_o3_2 = fltarr(size(gemsvar2.o3, /dimension))

for iy=0, merra2_sz2[1]-1 do begin
  for ix=0, merra2_sz2[0]-1 do begin
    merra2_o3_on_gems_pres2[ix, iy, *] = $
      interpol(merra2_o3_accum2[ix, iy, *], merra2_hpa2[ix, iy, 1:-1], gemsvar2.pressure[ix, iy, *])
  endfor
endfor


merra2_o3_on_gems_o3_2 = merra2_o3_on_gems_pres2[*, *, 1:-1] - merra2_o3_on_gems_pres2[*, *, 0:-2]


;=============================================================================
ecf1 = cldvar1.EffectiveCloudFraction
ecf2 = cldvar2.EffectiveCloudFraction

xidx = [50, 18, 22, 105]
yidx = [30, 63, 106, 262]

xidx = [50, 18, 22, 22, 22, 22, 22, 22, 22, 22, $
  22, 22, 22, 22, 22, 22, 22, 22, 22, 22, $
  105]
yidx = [30, 63, 106, 130, 145, 160, 180, 200, 220, 240, $
  260, 280, 300, 320, 340, 360, 380, 400, 420, 440, $
  262]

for i = 0, 3 do begin
  print, ecf1[xidx[i], yidx[i]]
endfor

; apply basic data filtering
nl = 24
;sel  = where(gemsvar.SolarZenithAngle le 88 and gemsvar.ProcessingQualityFlags eq 0 and gemsvar.FinalAlgorithmFlags eq 0, nprof, /null) 

if n_elements(xidx) eq n_elements(yidx) then begin
  okay = 1
endif else begin
  okay = 0
  print, 'Number of xidx and yidx is not same.'
  return
endelse

if okay then begin
  for ipix = 0, n_elements(xidx)-1 do begin

    arrsize = 10
    minidx = where(ecf[xidx[ipix]-arrsize:xidx[ipix]+arrsize, yidx[ipix]-arrsize:yidx[ipix]+arrsize] eq $
      min(ecf[xidx[ipix]-arrsize:xidx[ipix]+arrsize, yidx[ipix]-arrsize:yidx[ipix]+arrsize]))


    minindices = array_indices(fltarr(arrsize*2+1, arrsize*2+1), minidx)

    xidx_fixed = xidx[ipix] - arrsize + minindices[0]
    yidx_fixed = yidx[ipix] - arrsize + minindices[1]

    print, xidx_fixed, yidx_fixed
    print, ecf[xidx_fixed, yidx_fixed]
    ozprof     = gemsvar.O3[xidx_fixed, yidx_fixed, *]
    o3apriori   = gemsvar.O3Apriori[xidx_fixed, yidx_fixed, *]

    tpres = gemsvar.TropopausePressure[xidx_fixed, yidx_fixed]
    lat = gemsvar.latitude[xidx_fixed, yidx_fixed]
    lon = gemsvar.longitude[xidx_fixed, yidx_fixed]

    ;gems_avgk = fltarr(nl, nl)
    ;gems_avgk[*] = !values.f_nan

    ;gemsavgk = reform(gemsvar.AveragingKernel[xidx_fixed, yidx_fixed, *, *])
    ;;if finite(sonde_o3_total_du) eq 0 then continue
    ;IF n_elements(gemsavgk) ge 24 then begin
      ;;print, 'ifasdfafdasd'
      ;avgk       = reform(gemsavgk)

      ;for ip = 0 , nl-1  do begin
        ;;omi(sidx).avgk(*,ip) = avgk(*,ip)
        ;gems_avgk[*,ip] = avgk[*,ip]
      ;endfor
      ;gems_avgk = avgk
    ;endif else begin
          ;csontco = 0 & csontco200 = 0 & csonsco = 0 & csonsco200 = 0 
          ;csono3du0 = fltarr(nl)
    ;endelse

    merra2_ozprof = merra2_o3_on_gems_o3[xidx_fixed, yidx_fixed, *]

    ; plotting

    plot_margin = [0.18, 0.15, 0.10, 0.15]
    plot_xrange = [0,60]
    plot_yrange = [1000, 1]

    pres = reform(gemsvar.pressure[xidx_fixed, yidx_fixed, *])

    p1 = plot(ozprof, pres, /buffer, /overplot, /ylog, dim=[500, 600], $
      axis_style=1, $
      margin=plot_margin, $
      xrange=plot_xrange, $
      ytitle='Pressure [hPa]', $
      xtitle='O3 [DU]');color=[0, 0, 0], linestyle=0, name='GEMS O3P')
    ;p1.title = 'GEMS O3 Profile '+ datetime_str
    p1.color = [0, 0, 0]
    p1.yrange= plot_yrange
    p1.linestyle = 0
    p1.symbol= '+'
    p1.name = 'O3 profile'

    ;p1_1 = plot(merra2_ozprof, pres, /buffer, /overplot)
    ;p1_1.color = [0, 114, 177]
    ;p1_1.color = [255, 0, 0] ; for professor
    ;p1_1.linestyle = 2
    ;p1_1.symbol= '*'
    ;p1_1.name = 'MERRA2 O3 profile'

    p2 = plot(o3apriori, pres[0:23], /buffer, /overplot);, $ color='#e4a000', linestyle=1, name='Ozonesonde with GEMS O3P AVGK', /overplot)
    p2.color = [228, 160, 0]
    p2.color = [0, 0, 255] ; for professor
    p2.linestyle = 1
    p2.symbol= 's'
    p2.name = 'A priori'

    ;p3 = plot(omi_csonprof, pres, /buffer, /overplot);, $ color='#e4a000', linestyle=1, name='Ozonesonde with GEMS O3P AVGK', /overplot)
    ;p3.color = [86, 180, 232]
    ;p3.linestyle = 2
    ;p3.symbol= '*'
    ;p3.name = 'Ozonesonde with GEMS O3P AVGK'

    ; ========================================================================
    ; TROPOPAUSE PRESSURE PLOT
    ; ========================================================================
    ;tpres1 = [tpres, tpres]
    ;p4 = plot([0, 100], tpres1, /buffer, /overplot);, $ color='#56b4e8', linestyle=2, name='Ozonesonde', /overplot)
    ;p4.color = [86, 180, 232]
    ;p4.linestyle = 0
    ;;p4.symbol= 'D'
    ;p4.name = 'Tropopause'

    ; ========================================================================
    ; TEMPERATURE PLOT
    ; ========================================================================
    ;temp = reform(gemsvar.temperature[xidx_fixed, yidx_fixed, *])
    ;p5 = plot(temp, pres, /buffer, /current, $
      ;axis_style=0, $
      ;margin=plot_margin, $
      ;/ylog);, $ color='#56b4e8', linestyle=2, name='Ozonesonde', /overplot)
    ;p5.color = [0, 159, 115]
    ;p5.yrange= plot_yrange
    ;p5.xrange=[150, 310]
    ;p5.linestyle = 0
    ;p5.symbol= 'D'
    ;p5.name = 'Temperature'
    
    a5 = axis('y', $
      target = p1, $
      textpos=1, $
      major=4, $
      minor=8, $
      ;tickdir = 1, $
      location=[max(p1.xrange), 0, 0])

    a5 = axis('x', $
      target = p1, $
      textpos=1, $
      major=5, $
      minor=7, $
      tickdir = 1, $
      location=[0, max(p1.yrange), 0], $
      title = 'Temperature')

    ;TWMO, -0.002, temp, 0, 100000., pres*100, pres_tropo, temp_tropo, alt_tropo, 0

    leg = legend($
      target=[p1, $
        ;p1_1, $
        p2], $
        ;p4, $
        ;p5], $
      position=[59, 1.5],/data)

    t1 = text(0.5, 0.95, 'GEMS O3 Profile '+ datetime_str, $
      font_size=16, $
      /normal, $
      alignment=0.5, $
      vertical_alignment=0.5)

    lat_t = text(0.3, 0.78, 'LAT:' +string(lat, format='(f6.2)'), /normal)
    lat_t = text(0.3, 0.75, 'LON:' +string(lon, format='(f6.2)'), /normal)

    p1.save, outputpath + '/x' + $
      string(xidx_fixed, format='(i03)') + $
      'y' + string(yidx_fixed, format='(i03)') + '.png'
    p1.close
  endfor
ENDif

end


