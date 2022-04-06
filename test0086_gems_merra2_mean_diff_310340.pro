
;=============================================================================
; start
;=============================================================================
jd_list = timegen(start=julday(3, 29, 2021, 3, 45), $
  final=julday(3, 29, 2021, 3, 45), units='Hours')

;jd_list = timegen(start=julday(6, 21, 2021, 3, 45), $
  ;final=julday(6, 21, 2021, 3, 45), units='Hours')

itime = 0
caldat, jd_list[itime], month, day, year, hour, minute
yyyy = string(year, format='(i04)')
mm = string(month, format='(i02)')
dd = string(day, format='(i02)')
hh = string(hour, format='(i02)')
mi = string(minute, format='(i02)')

caldat, jd_list[itime]-1, month_y, day_y, year_y, hour_y, minute_y
yyyy_y = string(year_y, format='(i04)')
mm_y = string(month_y, format='(i02)')
dd_y = string(day_y, format='(i02)')
hh_y = string(hour_y, format='(i02)')
mi_y = string(minute_y, format='(i02)')

datetime_str = yyyy + mm + dd + '_' + hh + mi
datetime_y_str = yyyy_y + mm_y + dd_y

;=============================================================================
; READ GEMS O3P(MERRA2)
;=============================================================================
search_str = '/data/private/soodal/softcal_test/corrected/model_310_340/GK2_GEMS_L2_O3P_20210621_0345_winliminit310_prec000000_climML_b4x4_o20211116T130448Z_clima_maxiter10_ecf0.nc'

fl = file_search(search_str)
synthetic_o3p_data = ds_read_gems_l2_o3p(fl[-1])


synthetic_o3p = synthetic_o3p_data.o3
synthetic_o3p_apriori = synthetic_o3p_data.o3apriori
synthetic_o3p_p = synthetic_o3p_data.pressure

;=============================================================================
; read MERRA2
;=============================================================================
merra2 = ds_read_merra2_tavg3_3d_asm_nv_collocated_on_gems_l1c($
'/data/private/soodal/MERRA2_collocated_on_gems/' + $
'merra2_tavg3_3d_asm_collocated_on_gemsl1c_rad_' + datetime_str + '.nc4')

pres = merra2.pl

itime = 1
o3mmr = reform(merra2.o3)
pres = reform(merra2.PL)
temp = reform(merra2.T)
alt = reform(merra2.H)

lon = merra2.lon
lat = merra2.lat

;convert mass mxing ratio into vmr
o3vmr = ds_mmr2vmr(o3mmr)
merra2_hpa = pres/100. ; [174, 512, 72]
;merra2_hpa_mid = merra2_hpa[*, *, 1:-1] -merra2_hpa[*, *, 0:-2]

ds_vmr2du_3d, o3vmr, merra2_hpa, temp, alt, merra2_o3, increasing_height=0
merra2_sz = size(merra2_o3, /dimension)

merra2_o3_accum = fltarr(merra2_sz)
merra2_o3_accum[*, *, 0] = merra2_o3[*, *, 0]
for i=1, merra2_sz[2]-1 do begin
merra2_o3_accum[*, *, i] = merra2_o3_accum[*, *, i-1] + merra2_o3[*, *, i]
endfor

merra2_o3_on_gems_pres = fltarr(size(synthetic_o3p_data.pressure, /dimension))
merra2_o3_on_gems_o3 = fltarr(size(synthetic_o3p_data.o3, /dimension))

print, 'before the merra2 interpol loop.'
for iy=0, merra2_sz[1]-1 do begin
for ix=0, merra2_sz[0]-1 do begin
  merra2_o3_on_gems_pres[ix, iy, *] = $
    interpol(merra2_o3_accum[ix, iy, *], merra2_hpa[ix, iy, 1:-1], synthetic_o3p_data.pressure[ix, iy, *])
endfor
endfor

merra2_o3_on_gems_o3 = merra2_o3_on_gems_pres[*, *, 1:-1] - merra2_o3_on_gems_pres[*, *, 0:-2]

;=============================================================================
; read L2 CLD
;=============================================================================


cldfn = '/data/private/soodal/ln/GEMS/L2CLD/GK2_GEMS_L2_CLD_' + datetime_str + '_4x4.nc'
nier_gemscld = ds_read_gems_l2_o3p(cldfn)
ecf = nier_gemscld.EffectiveCloudFraction
ecf02nanidx = where(ecf lt -1e29, /null)
_ecf = ecf
_ecf[ecf02nanidx] = !values.f_nan
ecf = _ecf

;=============================================================================
; PLOTTING
;=============================================================================
pos = [0.12, 0.1, 0.9, 0.85]
plot_margin = [0.18, 0.10, 0.10, 0.15]
plot_xrange = [-10.,60.]
plot_yrange = [1000., 0.1]

outputpath = '/home/soodal/works/GEMS_O3P_analysis/plot/merra2_synthetic/'

for ilat=0., 40., 10. do begin
  ilon = 120
  valididx = where(lon ge ilon and lon lt ilon+10. and $
    lat ge ilat and lat lt ilat + 10. and $
    ecf lt 0.05 and ecf ge 0.0, /null)
  
  lat_str = string(ilat, format='(f05.1)')
  lon_str = string(ilon, format='(f06.1)')
  for ipix = 0 , n_elements(valididx)-1 do begin
    print, 'ipix:', ipix
    pixindices = array_indices(lon, valididx[ipix])
    ;if pixindices[0] eq 25 then begin

      merra2_o3_profile = reform(merra2_o3_on_gems_o3[pixindices[0], pixindices[1], *])
      gems_o3_profile = reform(synthetic_o3p[pixindices[0], pixindices[1], *])
      gems_o3_ap = reform(synthetic_o3p_apriori[pixindices[0], pixindices[1], *])
      
      if ipix ge 10 then begin
        break
      endif

      ;if n_elements(merra2_o3_profiles) eq 0 then begin
      if ipix eq 0 then begin
        merra2_o3_profiles = merra2_o3_profile
      endif else begin
        merra2_o3_profiles = [[merra2_o3_profiles], $
          [merra2_o3_profile]]
      endelse

      ;if n_elements(gems_o3_profiles) eq 0 then begin
      if ipix eq 0 then begin
        gems_o3_profiles = gems_o3_profile
      endif else begin
        gems_o3_profiles = [[gems_o3_profiles], $
          [gems_o3_profile]]
      endelse

      if ipix eq 0 then begin
        gems_o3_aps = gems_o3_ap
      endif else begin
        gems_o3_aps = [[gems_o3_aps], $
          [gems_o3_ap]]
      endelse

      ;if n_elements(gems_pressure_grids) eq 0 then begin 
      if ipix eq 0 then begin 
        gems_pressure_grids = reform(synthetic_o3p_data.pressure[pixindices[0], pixindices[1], *])
      endif else begin
        gems_pressure_grids = [[gems_pressure_grids], $
          [reform(synthetic_o3p_data.pressure[pixindices[0], pixindices[1], *])]]
      endelse

      gems_pressure_mid = exp((alog(synthetic_o3p_data.pressure[pixindices[0], pixindices[1], 0:-2]) + $
        alog(synthetic_o3p_data.pressure[pixindices[0], pixindices[1], 1:-1]))/2.)
      if n_elements(gems_pressure_mids) eq 0 then begin
        gems_pressure_mids = reform(gems_pressure_mid)
      endif else begin
        gems_pressure_mids = [[gems_pressure_mids], [reform(gems_pressure_mid)]]
      endelse

      ; PLOT GEMS O3 Profile
      if ipix eq 0 then begin
        command = "p1_" + string(ipix, format='(i03)') + " = plot(" + $
          "gems_o3_profiles[*, ipix], gems_pressure_mid, " + $
          "/buffer, " + $
          ;"overplot=ipix, " + $
          "/ylog, " + $
          "axis_style=0, " + $
          "xtitle='O3 [DU]', " + $
          "ytitle='Pressure [hPa]', " + $
          ;"title='GEMS O3P from Synthetic with MERRA2', " + $
          ;symbol='s',  + $
          "yrange=plot_yrange, " + $
          "xrange=plot_xrange, " + $
          "name='x' + string(pixindices[0], format='(i03)') + '_' + " + $
          "'y' + string(pixindices[1], format='(i03)'), " + $
          "color=[0, 0, 255], " + $
          ;"color='blue', " + $
          "transparency=50, " + $
          "position=pos)"
        dummy = execute(command)
      endif else begin
        command = "p1_" + string(ipix, format='(i03)') + " = plot(" + $
          "gems_o3_profiles[*, ipix], gems_pressure_mid, " + $
          "/buffer, " + $
          "/overplot, " + $
          "name='x' + string(pixindices[0], format='(i03)') + '_' + " + $
          "'y' + string(pixindices[1], format='(i03)'), " + $
          ;"'y' + string(pixindices[1], format='(i03)'))"
          ;"color=[0, 0, 255])";, " + $
          "color=[0, 0, 255], " + $
          ;"color='blue')";, " + $
          "transparency=50)"; + $
          ;"position=pos)"
        dummy = execute(command)
      endelse

      ; PLOT MERRA2 O3 Profile
      if ipix eq 0 then begin
        command = "p2_" + string(ipix, format='(i03)') + " = plot(" + $
          "merra2_o3_profiles[*, ipix], gems_pressure_mid, " + $
          "/buffer, " + $
          "/overplot, " + $
          "name='x' + string(pixindices[0], format='(i03)') + '_' + " + $
          "'y' + string(pixindices[1], format='(i03)'), " + $
          ;"'y' + string(pixindices[1], format='(i03)'))"
          ;"color=[255, 0, 0])";, " + $
          "color=[255, 0, 0], " + $
          ;"color='red')";, " + $
          "transparency=50)"; + $
          ;"position=pos)"
        dummy = execute(command)
      endif else begin
        command = "p2_" + string(ipix, format='(i03)') + " = plot(" + $
          "merra2_o3_profiles[*, ipix], gems_pressure_mid, " + $
          "/buffer, " + $
          "/overplot, " + $
          "name='x' + string(pixindices[0], format='(i03)') + '_' + " + $
          "'y' + string(pixindices[1], format='(i03)'), " + $
          ;"'y' + string(pixindices[1], format='(i03)'))"
          ;"color=[255, 0, 0])";, " + $
          "color=[255, 0, 0], " + $
          ;"color='red')";, " + $
          "transparency=50)"; + $
          ;"position=pos)"
        dummy = execute(command)
      endelse

      ; PLOT GEMS A-priori O3 Profile
      if ipix eq 0 then begin
        command = "pap_" + string(ipix, format='(i03)') + " = plot(" + $
          "gems_o3_aps[*, ipix], gems_pressure_mid, " + $
          "/buffer, " + $
          "/overplot, " + $
          "name='x' + string(pixindices[0], format='(i03)') + '_' + " + $
          "'y' + string(pixindices[1], format='(i03)'), " + $
          ;"'y' + string(pixindices[1], format='(i03)'))"
          ;"color=[100, 100, 100])";, " + $
          "color=[100, 100, 100], " + $
          ;"color='red')";, " + $
          "transparency=50)"; + $
          ;"position=pos)"
        dummy = execute(command)
      endif else begin
        command = "pap_" + string(ipix, format='(i03)') + " = plot(" + $
          "gems_o3_aps[*, ipix], gems_pressure_mid, " + $
          "/buffer, " + $
          "/overplot, " + $
          "name='x' + string(pixindices[0], format='(i03)') + '_' + " + $
          "'y' + string(pixindices[1], format='(i03)'), " + $
          ;"'y' + string(pixindices[1], format='(i03)'))"
          ;"color=[100, 100, 100])";, " + $
          "color=[100, 100, 100], " + $
          ;"color='red')";, " + $
          "transparency=50)"; + $
          ;"position=pos)"
        dummy = execute(command)
      endelse
    ;endif

  endfor

  gems_merra2_diff = gems_o3_profiles - merra2_o3_profiles
  gems_merra2_diff_mean = reform(mean(gems_merra2_diff, dim=2))
  gems_merra2_diff_stdev = reform(stddev(gems_merra2_diff, dim=2))
  gems_pressure_mid_mean = reform(exp(mean(alog(gems_pressure_mids), dim=2)))

  merra2_mean = reform(mean(merra2_o3_profiles, dim=2))

  p_ = plot([0, 0], [0.1, 1000], $
    /buffer, $
    color=[150, 150, 150], $
    /overplot)

  ;p3 = barplot(reform(exp(mean(alog(gems_pressure_mids), dim=2))), $
    ;reform(mean(gems_merra2_diff, dim=2)), $
    ;/buffer, $
    ;/ylog, $
    ;;index=0, nbar=1, $
    ;;fill_color='red', $
    ;/horizontal, $
    ;name='GEMS O3P - MERRA2 O3P', $
    ;/overplot);, $
    ;;yrange=plot_yrange, $
    ;;xrange=[-10, 60], $
    ;;position=pos)

  p3 = plot($
    gems_merra2_diff_mean, $
    gems_Pressure_mid_mean, $
    /buffer, $
    ;/ylog, $
    ;index=0, nbar=1, $
    ;fill_color='red', $
    ;/horizontal, $
    name='GEMS O3P - MERRA2 O3P', $
    /overplot);, $
    ;yrange=plot_yrange, $
    ;xrange=[-10, 60], $
    ;position=pos)

  a1 = axis('x', $
    target = p3, $
    ;textpos=0, $
    major=5, $
    minor=7, $
    tickvalues=[-10, 0, 10, 20, 30, 40, 50, 60], $
    tickname=[-10, 0, 10, 20, 30, 40, 50, 60], $
    title='O3 [DU]', $
    ;tickdir = 0, $
    ;tickunits='Scientific', $
    location='bottom');[0, min(p2.yrange), 0] );, $
  a2 = axis('y', $
    target = p3, $
    ;textpos=0, $
    major=5, $
    minor=7, $
    tickvalues=[0.1, 1.0, 10.0, 100.0, 300.0, 500.0, 700.0, 1000.0], $
    tickname=[0.1, 1, 10, 100, 300, 500, 700, 1000], $
    ;tickdir = 0, $
    ;tickunits='Scientific', $
    title='Pressure [hPa]', $
    location='left');[0, min(p2.yrange), 0] );, $
  a3 = axis('x', $
    target = p3, $
    ;textpos=0, $
    major=5, $
    minor=7, $
    tickvalues=[-10, 0, 10, 20, 30, 40, 50, 60], $
    tickname=['', '', '', '', '', '', '', ''], $
    ;tickdir = 0, $
    ;tickunits='Scientific', $
    ;title='Pressure [hPa]', $
    location='top');[0, min(p2.yrange), 0] );, $
  a4 = axis('y', $
    target = p3, $
    ;textpos=0, $
    major=5, $
    minor=7, $
    tickvalues=[0.1, 1.0, 10.0, 100.0, 300.0, 500.0, 700.0, 1000.0], $
    tickname=['', '', '', '', '', '', '', ''], $
    ;tickdir = 0, $
    ;tickunits='Scientific', $
    ;title='Pressure [hPa]', $
    location='right');[0, min(p2.yrange), 0] );, $

  leg = legend(target=[p1_000, p1_001, p1_002, p1_004, p1_005, p1_005, $
          p1_007, p1_008, p1_009], $
        position=[55, 1], /data)

  plot_xrange2=[-250, 250]
  p_center = plot([0, 0], plot_yrange, $
    axis_style=0, $
    /current, $
    xrange=plot_xrange2, $
    yrange=plot_yrange, $
    color=[150, 150, 150], $
    pos=pos)

  ;a3 = axis('x', $
    ;target = p_center, $
    ;;textpos=0, $
    ;major=5, $
    ;minor=7, $
    ;tickvalues=[-250, -200, -150, -100, -50, 0, 50, 100, 150, 200, 250], $
    ;;tickdir = 0, $
    ;;tickunits='Scientific', $
    ;title='DIFF RATIO [%]', $
    ;location='top');[0, min(p2.yrange), 0] );, $

  title_text = text(0.5, 0.95, 'GEMS O3P from Synthetic with MERRA2', $
    alignment=0.5, font_size=18)

  ; MEAN DIFF BARPLOT for log axis
  for ilevel = 0, 23 do begin
    x = gems_merra2_diff_mean[ilevel]
    x0 = pos[0] + (-plot_xrange[0]) / (plot_xrange[1]-(plot_xrange[0])) * (pos[2] - pos[0]) 
    x1 = pos[0] + (x-plot_xrange[0]) / (plot_xrange[1]-(plot_xrange[0])) * (pos[2] - pos[0])

    y = gems_Pressure_mid_mean[ilevel]
    y0 = pos[3] - alog(y) / (alog(plot_yrange[0]) - alog(plot_yrange[1])) * (pos[3] - pos[1])
    y0 = pos[3] - (alog(y) - alog(plot_yrange[1])) / (alog(plot_yrange[0]) - alog(plot_yrange[1])) * (pos[3] - pos[1])

    x_list = [x0, x1, x1, x0]
    y_list = [y0-0.005, y0-0.005, y0+0.005, y0+0.005]

    if x ge 0 then begin
      dummy = polygon(x_list, y_list, /current, $
        /norm, $
        transparency=50, $
        fill_color=[255, 100, 100])
    endif else begin
      dummy = polygon(x_list, y_list, /current, $
        /norm, $
        transparency=50, $
        fill_color=[100, 100, 255])
    endelse
  endfor

  ; MEAN DIFF RATIO BARPLOT for log axis
  ;for ilevel = 0, 23 do begin
    ;xstart = 0.
    ;x = gems_merra2_diff_mean[ilevel] / merra2_mean[ilevel] * 100.
    ;x0 = pos[0] + (xstart-plot_xrange2[0]) / (plot_xrange2[1]-(plot_xrange2[0])) * (pos[2] - pos[0]) 
    ;x1 = pos[0] + (x-plot_xrange2[0]) / (plot_xrange2[1]-(plot_xrange2[0])) * (pos[2] - pos[0])

    ;y = gems_Pressure_mid_mean[ilevel]
    ;y0 = pos[3] - alog(y) / (alog(plot_yrange[0]) - alog(plot_yrange[1])) * (pos[3] - pos[1])
    ;y0 = pos[3] - (alog(y) - alog(plot_yrange[1])) / (alog(plot_yrange[0]) - alog(plot_yrange[1])) * (pos[3] - pos[1])

    ;x_list = [x0, x1, x1, x0]
    ;y_list = [y0-0.007, y0-0.007, y0+0.007, y0+0.007]

    ;if x ge 0 then begin
      ;dummy = polygon(x_list, y_list, /current, $
        ;/norm, $
        ;transparency=50, $
        ;fill_color=[255, 150, 150])
    ;endif else begin
      ;dummy = polygon(x_list, y_list, /current, $
        ;/norm, $
        ;transparency=50, $
        ;fill_color=[150, 150, 255])
    ;endelse
  ;endfor

  ;p3.title = cities_name[icity] + ' ' + datetime_str
  ;p3.title.font_size=16
  ;p3.yrange = plot_yrange
  ;;p2.color = [228, 160, 0]
  ;p2.color = [0, 0, 0]
  ;p2.linestyle = 0
  ;;p2.symbol= 's'
  ;p2.name = 'A priori'
  p_.save, outputpath + 'gems_merra2_diff_profile_310340_lat' + lat_str + '_lon' + lon_str + '.png'
  p_.close
endfor



end

