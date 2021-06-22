
imonth = 8

smonth = string(imonth, format='(i02)')

gemso3p_outfilelist = file_search('/data2/L2_GEMS/val_1008/GK2_GEMS_O3P_2020' + smonth + '01_0345.nc4')
gemso3p_outfilelist = file_search('/data2/L2_GEMS/nier_L1C/GK2B_GEMS_L2_O3P_20210322_*0045_4x4.nc4')

year = 2020
month = 8
day = 1
hour = 03
minute = 45

for iday = 4, 31 do begin
collocate_gemsl2o3p_omil2profoz_on_gems_grid, year, month, iday, hour, minute, $
  output

if n_elements(output) gt 0 then begin
  latlt10idx = where(output.omi_latitude lt 10 and $
    output.omi_effectivecloudfraction lt 0.2 and $
    output.omi_effectivecloudfraction ge 0 and $
    output.gems_latitude lt 10 and $
    output.gems_effectivecloudfraction lt 0.2 and $
    output.gems_effectivecloudfraction ge 0, numidx, /null)

  if numidx ge 1 then begin 
    omi_o3 = output.omi_o3[*, latlt10idx[0]]
    omi_o3ap = output.omi_o3_apriori[*, latlt10idx[0]]
    omipres = output.omi_pressure[*, latlt10idx[0]]
    gems_o3 = output.gems_o3[*, latlt10idx[0]]
    gems_o3ap = output.gems_o3_apriori[*, latlt10idx[0]]
    gemspres = output.gems_pressure[*, latlt10idx[0]]

    pos = [0.10, 0.1, 0.8, 0.85]
    leg_pos = [0.85, 0.8]

    p1 = plot(omi_o3, omipres[0:23], color='black', /buffer, $
      title='GEMS O3PROFILE 2020-08-01T0345', $
      xtitle='O3[DU]', $
      ytitle='Pressure[hPa]', $
      symbol='*', $
      font_size=10, $
      yrange=[1000, 0], /ylog, pos=pos, name='OMI O3 profile')

    p2 = plot(omi_o3ap, omipres[0:23], color='grey', /overplot, /buffer, $
      symbol='diamond', $
      yrange=[1000, 0], /ylog, pos=pos, $
      linestyle=1, $
      name='OMI O3 a priori profile')

    p3 = plot(gems_o3, gemspres[0:23], color='blue', /overplot, /buffer, $
      symbol='+', $
      yrange=[1000, 0], /ylog, pos=pos, name='GEMS O3 profile')

    p4 = plot(gems_o3ap, gemspres[0:23], color='sky blue', /overplot, /buffer, $
      symbol='triangle', $
      linestyle=1, $
      yrange=[1000, 0], /ylog, pos=pos, name='GEMS O3 a priori profile')

    leg = legend(target=[p1, p2, p3, p4], pos=leg_pos)

    yyyy = string(year, format='(i04)')
    mm = string(month, format='(i02)')
    dd = string(iday, format='(i02)')
    hh = string(hour, format='(i04)')
    mi = string(minute, format='(i04)')

    pngfile = './plot/gems_omi_profile_comparison_' + $
      yyyy + mm + dd + hh + mi + 'z.png'

    gemslat1 = output.gems_latitude[latlt10idx[0]]
    gemslon1 = output.gems_longitude[latlt10idx[0]]
    gemsecf1 = output.gems_effectivecloudfraction[latlt10idx[0]]

    ;t1 = text(0.45, 0.25, 'Pressure Quality Flag:' + string(info[1], format='(I1)'))
    ;t2 = text(0.45, 0.20, 'Cloud Pressure:' + string(cp[xidx, yidx]))
    t1 = text(0.45, 0.30, 'Longitude:' + string(gemslon1, format='(f7.2)'), font_size=10)
    t1 = text(0.45, 0.25, 'Latitude:' + string(gemslat1, format='(f7.2)'), font_size=10)
    t2 = text(0.45, 0.20, 'Effective Cloud Fraction:' + string(gemsecf1, format='(F7.2)'))

    p1.save, pngfile
    p1.close
    ; send image to pc
    scp_dest = 'soodal@164.125.38.179:/home/soodal/works/plot/gems_omi_comp'
    spawn, 'scp -P18742 -p ' + pngfile +  ' ' + scp_dest
  endif

  endif
endfor
stop


gemso3p = ds_read_gems_l2_o3p(gemso3p_outfilelist[0], varlist=varlist)

gems_ak = gemso3p.AveragingKernel

ds_gems_l2o3p_accum, gemso3p, o3under10k_gems, height=10.

gems_o3 = [gems_o3, o3under10k_gems[gems_pohang_idx[0], gems_pohang_idx[1]]]
;=============================================================================
; Ozonesonde
;=============================================================================

o3sonde_path = '/data1/gems/o3p/works/GEMS_O3P_analysis/data/ozonesonde/pohang/'
o3sonde_filelist = file_search(o3sonde_path + '2020' + smonth  + '*_L0.txt')

read_ozonesonde_pohang, o3sonde_filelist[0], sonde

vmr = sonde.ozon_mpa/1000/100/sonde.pres_hpa

ds_vmr2du_1d, vmr, sonde.pres_hpa, sonde.temp_degc + 273.15, sonde.geop_h_gpm, o3_du_layer

nsonde = n_elements(vmr)
o3under10k_sonde = 0.
for i=0, nsonde-1 do begin
  if sonde.geop_h_gpm[i] lt 10000 then begin
    o3under10k_sonde += o3_du_layer[i]
  endif
endfor
 
sonde_o3 = [sonde_o3, o3under10k_sonde]

stop



gemsyyyymmdd = strmid(gemso3p_outfilelist, 16, 8, /reverse)
gemshhmi = strmid(gemso3p_outfilelist, 7, 4, /reverse)
gemsyyyymmdd_hhmi = strmid(gemso3p_outfilelist, 16, 13, /reverse)

total_omivals = []
total_gemsvals = []
on_dailyplot = 0
for ifile=0, n_elements(gemso3p_outfilelist)-1 do begin
  year = fix(strmid(gemsyyyymmdd[ifile], 0, 4))
  mon = fix(strmid(gemsyyyymmdd[ifile], 4, 2))
  day = fix(strmid(gemsyyyymmdd[ifile], 6, 2))
  hour = fix(strmid(gemshhmi[ifile], 0, 2))
  mi = fix(strmid(gemshhmi[ifile], 2, 2))

  omivals = []
  gemsvals = []
  collocate_gemsl2o3p_omo3pr, year, mon, day, hour, mi, gemsvals, omivals, $
    on_under_300=1

  if n_elements(gemsvals) ge 2 and n_elements(omivals) ge 2 then begin
    if on_dailyplot then begin
      plot_gems_validation, gemsvals, omivals, $
        filename='./plot/gems_l2_o3p_val_with_omi_' + gemsyyyymmdd[ifile] + $
          'T' + gemshhmi[ifile]+'_under300hpa.png_ecf02', $
          range=[0, 100], $
          delta=2.0
    endif
    total_omivals = [total_omivals, omivals]
    total_gemsvals = [total_gemsvals, gemsvals]
  endif
endfor

; plot tropospheric column ozone
plot_gems_validation, total_gemsvals, total_omivals, $
  filename='./plot/gems_l2_o3p_val_with_omi_' + $
    strmid(gemsyyyymmdd[0], 0, 6) +'_under300hpa_ecf02.png', $
  cblim=[0, 100], $
  range=[0, 100], $
  delta=2.0

end
