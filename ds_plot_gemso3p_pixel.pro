PRO calc_cum ,use_spline, pres0, cum0, Ptop, Pbot, out
   pres = pres0
   cum  = [0, cum_total(cum0)]
   out = !values.f_nan

   temp = interpol(cum, (pres), ([ptop, pbot]))
   out = temp(1) - temp(0)
END


;============================
; MAIN PROGRAM 
; (1) read omil2file
; (2) screen sonpxl with no collocation with OMI  if colist are gvien
; (3) find collocation between OMI and SONDE within 0.5 for lon, 1.5 for lat, for 12 hours
;pro load_gems_sonde, sondes, gemsfile,  gems,  nprof, show=show, colist=colist,gemssta=gemssta

pro ds_plot_gemso3p_pixel, gemsfile, xidx, yidx, outputpath=outputpath, project=project, sub=sub, suffix=suffix


if not keyword_set(suffix) then begin
  put_suffix = 0
  suffix = ''
endif else begin
  put_suffix = 1
endelse

gems_basename = file_basename(gemsfile)

print, 'This procedures use the Index start with 1 not 0.'

if not keyword_set(outputpath) then begin
  outputpath = './plot/'
endif

if keyword_set(project) then begin
  outputpath = outputpath + project + '/'
endif else begin
  project = ''
endelse

if keyword_set(sub) then begin
  outputpath = outputpath + sub + '/'
endif else begin
  sub = ''
endelse


;==================================
;(1)  read gems l2 o3p file
;==================================  
datetimepos = stregex(gemsfile, '[[:digit:]]{8}_[[:digit:]]{4}')

print, gemsfile
print, datetimepos
year = fix(strmid(gemsfile, datetimepos, 4))
print, strmid(gemsfile, datetimepos, 4)
month = fix(strmid(gemsfile, datetimepos+4, 2))
print, strmid(gemsfile, datetimepos+4, 2)
day = fix(strmid(gemsfile, datetimepos+6, 2))
print, strmid(gemsfile, datetimepos+6, 2)
hour = fix(strmid(gemsfile, datetimepos+9, 2))
print, strmid(gemsfile, datetimepos+9, 2)
minute = fix(strmid(gemsfile, datetimepos+11, 2))
print, strmid(gemsfile, datetimepos+11, 2)

;caldat, jd_list[idate], month, day, year, hour, minute
yyyy = string(year, format='(i04)')
mm = string(month, format='(i02)')
dd = string(day, format='(i02)')
hh = string(hour, format='(i02)')
mi = string(minute, format='(i02)')

datetime_str = yyyy + mm + dd + '_' + hh + mi

outputpath = outputpath + '/'  + datetime_str + '/point_profile/'

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

gemsvar = ds_read_gems_l2_o3p(gemsfile, varlist=varlist)

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
  fnlncfile = '/data/MODEL/FNL/' + yyyy + '/fnl_' + yyyy + mm + dd + '_' + '00_00.grib2.nc'

  print, fnlncfile
  if file_test(fnlncfile) then begin
    print, fnlncfile, ' is exist.'
    fnl = ds_read_fnl_nc(fnlncfile) 
  endif
  for ipix = 0, n_elements(xidx)-1 do begin

    ozprof     = gemsvar.O3[xidx[ipix]-1, yidx[ipix]-1, *]
    o3apriori   = gemsvar.O3Apriori[xidx[ipix]-1, yidx[ipix]-1, *]

    tpres = gemsvar.TropopausePressure[xidx[ipix]-1, yidx[ipix]-1]
    lat = gemsvar.latitude[xidx[ipix]-1, yidx[ipix]-1]
    lon = gemsvar.longitude[xidx[ipix]-1, yidx[ipix]-1]
    ecf = gemsvar.EffectiveCloudFractionUV[xidx[ipix]-1, yidx[ipix]-1]
    cp = gemsvar.CloudPressure[xidx[ipix]-1, yidx[ipix]-1]
    sza = gemsvar.SolarZenithAngle[xidx[ipix]-1, yidx[ipix]-1]

    ;gems_avgk = fltarr(nl, nl)
    ;gems_avgk[*] = !values.f_nan

    ;gemsavgk = reform(gemsvar.AveragingKernel[xidx[ipix], yidx[ipix], *, *])
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


    ; plotting

    plot_margin = [0.18, 0.15, 0.10, 0.15]
    plot_xrange = [0,60]
    plot_yrange = [1000, 0.08746]

    pres = reform(gemsvar.pressure[xidx[ipix]-1, yidx[ipix]-1, *])
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

    p2 = plot(o3apriori, pres[0:23], /buffer, /overplot);, $ color='#e4a000', linestyle=1, name='Ozonesonde with GEMS O3P AVGK', /overplot)
    p2.color = [228, 160, 0]
    p2.linestyle = 1
    p2.symbol= 's'
    p2.name = 'A priori'

    ;p3 = plot(omi_csonprof, pres, /buffer, /overplot);, $ color='#e4a000', linestyle=1, name='Ozonesonde with GEMS O3P AVGK', /overplot)
    ;p3.color = [86, 180, 232]
    ;p3.linestyle = 2
    ;p3.symbol= '*'
    ;p3.name = 'Ozonesonde with GEMS O3P AVGK'

    tpres1 = [tpres, tpres]
    p4 = plot([0, 100], tpres1, /buffer, /overplot);, $ color='#56b4e8', linestyle=2, name='Ozonesonde', /overplot)
    p4.color = [86, 180, 232]
    p4.linestyle = 0
    ;p4.symbol= 'D'
    p4.name = 'Tropopause'

    temp = reform(gemsvar.temperature[xidx[ipix]-1, yidx[ipix]-1, *])
    p5 = plot(temp, pres, /buffer, /current, $
      axis_style=0, $
      margin=plot_margin, $
      /ylog);, $ color='#56b4e8', linestyle=2, name='Ozonesonde', /overplot)
    p5.color = [0, 159, 115]
    p5.yrange= plot_yrange
    p5.xrange=[150, 310]
    p5.linestyle = 0
    p5.symbol= 'D'
    p5.name = 'Temperature'
    
    a5 = axis('y', $
      target = p5, $
      textpos=1, $
      major=4, $
      minor=8, $
      ;tickdir = 1, $
      location=[max(p5.xrange), 0, 0])

    a5 = axis('x', $
      target = p5, $
      textpos=1, $
      major=5, $
      minor=7, $
      tickdir = 1, $
      location=[0, max(p5.yrange), 0], $
      title = 'Temperature')

    ;TWMO, -0.002, temp, 0, 100000., pres*100, pres_tropo, temp_tropo, alt_tropo, 0

    leg_list = [p1, p2, p4, p5]
    if file_test(fnlncfile) then begin
      TMP_P0_L100_GLL0 = fnl.TMP_P0_L100_GLL0
      fnl_lon_0 = fnl.lon_0
      fnl_lat_0 = fnl.lat_0
      ds_xymake2d, fnl_lon_0, fnl_lat_0, fnl_lon, fnl_lat
      idx = search_closest_pixel(fnl_lon, fnl_lat, lon, lat, maxlimit=1.0)
      if idx ne -999 then begin
        indices = array_indices(fnl_lon, idx)
        fnl_pres = fnl.lv_ISBL0/100.
        p6 = plot(TMP_P0_L100_GLL0[indices[0], indices[1], *], fnl_pres, /buffer, /current, $
          axis_style=0, $
          margin=plot_margin, $
          /ylog);, $ color='#56b4e8', linestyle=2, name='Ozonesonde', /overplot)
        p6.color = [240, 228, 66]
        p6.yrange= plot_yrange
        p6.xrange=[150, 310]
        p6.linestyle = 0
        p6.symbol= 'tu'
        p6.name = 'FNL Temperature'
        leg_list = [p1, p2, p4, p5, p6]
      ENDif
    ENDif

    leg = legend(target=leg_list, position=[59, 30],/data)

    t1 = text(0.5, 0.95, 'GEMS O3 Profile '+ datetime_str, $
      font_size=16, $
      /normal, $
      alignment=0.5, $
      vertical_alignment=0.5)

    lat_t = text(0.4, 0.78, 'LAT:' +string(lat, format='(f6.2)'), /normal)
    lon_t = text(0.4, 0.75, 'LON:' +string(lon, format='(f6.2)'), /normal)
    ecf_t = text(0.4, 0.72, 'ECF:' +string(ecf, format='(f6.2)'), /normal)
    cp_t = text(0.4, 0.69, 'CP:' +string(cp, format='(f6.2)'), /normal)
    sza_t = text(0.4, 0.66, 'SZA:' +string(sza, format='(f6.2)'), /normal)

    if put_suffix then begin
      outfn = outputpath + '/x' + $
        string(xidx[ipix], format='(i03)') + $
        'y' + string(yidx[ipix], format='(i03)') + '_' + suffix + '.png'
    ENDif else begin
      outfn = outputpath + '/x' + $
        string(xidx[ipix], format='(i03)') + $
        'y' + string(yidx[ipix], format='(i03)') + '.png'
    endelse

    p1.save, outfn
    p1.close
  endfor
ENDif
END
