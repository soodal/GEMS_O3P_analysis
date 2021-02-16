
o3p_outfilelist = file_search('/data2/L2_GEMS/val_1008/GK2_GEMS_O3P_202008??_0345.nc4')
;o3p_outfilelist = file_search('/data2/L2_GEMS/val_1008/GEMS_O3T')


gemsyyyymmdd = strmid(o3p_outfilelist, 16, 8, /reverse)
gemshhmi = strmid(o3p_outfilelist, 7, 4, /reverse)
gemsyyyymmdd_hhmi = strmid(o3p_outfilelist, 16, 13, /reverse)

total_o3pvals = []
total_o3tvals = []
total_lats = []
total_lons = []
total_ecfs = []

yn_daily = 0


savfile = './data/gems_m08_o3p_o3t_comparison.sav'

nlocs = []
flags = []
if not file_test(savfile) then begin
;if 1 then begin
  for ifile=0, n_elements(o3p_outfilelist)-1 do begin
    year = fix(strmid(gemsyyyymmdd[ifile], 0, 4))
    mon = fix(strmid(gemsyyyymmdd[ifile], 4, 2))
    day = fix(strmid(gemsyyyymmdd[ifile], 6, 2))
    hour = fix(strmid(gemshhmi[ifile], 0, 2))
    mi = fix(strmid(gemshhmi[ifile], 2, 2))

    omivals = []
    gemsvals = []
    o3p = ds_read_gems_l2_o3p(o3p_outfilelist[ifile])
    o3p_cao3 = reform(o3p.ColumnAmountO3[0, *, *])
    nanidx = where(o3p_cao3 lt 0, /null)
    o3p_cao3[nanidx] = !values.f_nan
    
    o3p_ecf = o3p.EffectiveCloudFractionUV
    o3p_ecf[where(o3p_ecf lt 0)] =!values.f_nan

    o3p_lon = o3p.Longitude
    o3p_lat = o3p.Latitude


    o3tfn = strmid(o3p_outfilelist[ifile], 0, 24) + 'GEMS_O3T' + $
      strmid(o3p_outfilelist[ifile], 36, 17)
    o3t = ds_read_gems_l2_o3t(o3tfn)
    o3tcao3 = o3t.ColumnAmountO3 
    o3tlat = o3t.latitude
    o3tlon = o3t.longitude
    ;nanidx = where(o3tlon lt -180 or o3tlon gt 360, /null)
    o3tlon[nanidx]=!values.f_nan
    o3tlat[nanidx]=!values.f_nan
    nanidx = where(o3tlat lt -90 or o3tlon gt 90, /null)
    o3tlon[nanidx]=!values.f_nan
    o3tlat[nanidx]=!values.f_nan
    

    ;pngfile = 'o3t_daily.png'
    ;c = contour(o3tcao3, o3tlon, o3tlat, /fill, /buffer, $
      ;zrange=[200, 360])

    ;c.save, pngfile
    ;c.close

    ;scp_dest = 'soodal@164.125.38.179:/home/soodal/WindowsHome/Pictures/scp'
    ;spawn, 'scp -P18742 -p ' + pngfile + $
      ;' ' + scp_dest

    notroiidx = where(o3tlat lt -5. or o3tlat gt 45, /null)
    o3tcao3[notroiidx] = !values.f_nan
    nanidx = where(o3tcao3 lt 0 or finite(o3tcao3) ne 1, /null)
    o3tcao3[nanidx] = !values.f_nan
    

    flag = o3t.finalalgorithmflags
    flags = [flags, flag]
    nanidx = where(flag ne 0)
    ;o3tcao3[nanidx] = !values.f_nan
    

    collocate_gems_l2_o3t2o3p, o3tcao3, o3t_cao3_4x4

    ;o3t_cao3_4x4 = congrid(o3tcao3, 512, 174)
    ;o3t_cao3_4x4 = transpose(o3t_cao3_4x4)

    ;pngfile = 'o3t_4x4_daily.png'
    ;c = contour(o3tcao3, o3tlon, o3tlat, /fill, /buffer, $
      ;zrange=[200, 360])

    ;c.save,pngfile 
    ;c.close
    ;scp_dest = 'soodal@164.125.38.179:/home/soodal/WindowsHome/Pictures/scp'
    ;spawn, 'scp -P18742 -p ' + pngfile + $
      ;' ' + scp_dest

    o3tval = o3t_cao3_4x4
    o3pval = o3p_cao3

    idx = where(o3p_lon ge 80 and o3p_lon le 150 and o3p_lat ge -5 and $
      o3p_lat le 45 and finite(o3t_cao3_4x4) eq 1 and $
      o3t_cao3_4x4 ge 50 and o3p_ecf le 0.1, nloc)
    nlocs = [nlocs, nloc]
    print, nloc

    ;plot_sat_proj, o3t_Cao3_4x4, o3p.longitude, o3p.latitude, range=[200, 400], /scp_send


    if n_elements(o3tval) ge 2 and n_elements(o3pval) ge 2 then begin
      if yn_daily then begin
        plot_gems_validation, o3pval, o3tval, $
          filename='./plot/gems_l2_o3p_val_o3t_' + gemsyyyymmdd[ifile] + $
            'T' + gemshhmi[ifile]+'.png', $
          xtitle='GEMS_O3P', ytitle='GEMS_O3T', $
          range=[200, 350], delta=2, $
          cblim=[0, 400]
      endif
      total_o3pvals = [total_o3pvals, o3pval]
      total_o3tvals = [total_o3tvals, o3tval]
      total_lats = [total_lats, o3p_lat]
      total_lons = [total_lons, o3p_lon]
      total_ecfs = [total_ecfs, o3p_ecf]

    endif
  endfor
  save, filename=savfile
endif else begin
  restore, savfile
endelse

roiidx = where(finite(total_o3pvals) eq 1 and finite(total_o3tvals) eq 1 and $
  total_o3tvals ge 0 and $
  total_o3pvals ge 0 and $
  total_ecfs lt 0.1 and $
  total_lons gt 80 and $
  total_lons lt 150 and $
  total_lats gt -5 and $
  total_lats lt 45 $
  , /null)

plot_gems_validation, total_o3pvals[roiidx], total_o3tvals[roiidx], $
  filename='./plot/gems_l2_o3p_val_o3t.png', $
  ct=22, $
  xtitle='GEMS_O3P', ytitle='GEMS_O3T', $
  range=[200, 350], delta=2, $
  cblim=[0, 4000]

diffratio = (total_o3pvals[roiidx]-total_o3tvals[roiidx])/total_o3tvals[roiidx]*100
pdf = histogram(diffratio, locations=xbin, binsize=0.1)

phist = barplot(xbin, pdf, $
  xtitle='GEMS_O3P - GEMS_O3T[%]', $
  ytitle='Frequency',$
  fill_color='pink', $
  linestyle=6, $
  title='Histogram', $
  xrange=[-5, 5], $
  width=1.0, $
  /buffer)

yfit = gaussfit(xbin, pdf, coeff)

p1 = plot(xbin, yfit, 'r', /overplot, /buffer)
t1 = text(0.25, 0.83, $
  'Center: ' + string(coeff[1], format='(F7.3)'), 'r', $
  /normal)
t2 = text(0.25, 0.80, $
  'Width: ' + string(coeff[2], format='(F7.3)'), 'r', $
  /normal)
t3 = text(0.25, 0.77, $
  '#: ' + string(n_elements(total_o3tvals[roiidx]), format='(i7)'), $
  'r', /normal)

p2 = plot(fltarr(10), findgen(10)*10000, '--g', /overplot, /buffer)
p1.yrange=[0, max(pdf)*1.1]


pngfile = './histogram_gems_o3p_o3t_diff_ratio.png'
phist.save, pngfile
scp_dest = 'soodal@164.125.38.179:/home/soodal/WindowsHome/Pictures/scp'
spawn, 'scp -P18742 -p ' + pngfile + $
  ' ' + scp_dest

toc

end
