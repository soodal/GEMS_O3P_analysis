
imonth = 8

smonth = string(imonth, format='(i02)')

gemso3p_outfilelist = file_search('/data2/L2_GEMS/val_1008/GK2_GEMS_O3P_2020' + smonth + '01_0345.nc4')
gemso3p_outfilelist = file_search('/data2/L2_GEMS/nier_L1C/GK2B_GEMS_L2_O3P_20210322_*0045_4x4.nc4')

year = 2020
month = 6
day = 16
hour = 03
minute = 45

;gemspath = '/data1/gems/o3p/works/GEMS_O3P/out/310340_KARI/'

pos = [0.10, 0.1, 0.95, 0.85]
leg_pos = [0.90, 0.80]


for wininit=300, 310, 5 do begin
  winlim_str = string(wininit, format='(i03)')
  gemspath = '/data1/gems/o3p/works/GEMS_O3P/out/'+winlim_str+'340_KARI/'
  gemsfile = file_search(gemspath + 'GK2B_GEMS_L2_*.nc4')
  for iday = 16, 16 do begin
    if n_elements(gemsfile) ne 1 then begin
      stop
    endif
    colloc_omil2profoz_on_gems, year, month, iday, hour, minute, $
      output, gemsfile = gemsfile[0]

    datetime_str = string(year, format='(i04)') + '-' $
      + string(month, format='(i02)') + '-' $
      + string(iday, format='(i02)') + 'T' $
      + string(hour, format='(i02)') + ':' $
      + string(minute, format='(i02)')


    if n_elements(output) gt 0 then begin
      idx = where($
        output.gemsl2o3p.longitude ge 127 and $
        output.gemsl2o3p.longitude lt 129 and $
        output.colloc_omi_ecf_on_gems lt 0.2 and $
        output.colloc_omi_ecf_on_gems ge 0.0 $
        , numidx, /null)

      indices = array_indices(output.gemsl2o3p.longitude, idx)

      gems_o3_sz = size(output.gemsl2o3p.o3)
      if numidx ge 1 then begin 

        ; each pixel

        for ipix=0, numidx-1 do begin
          gems_pres_pix1 = fltarr(24)
          gems_pres_pix1[*] = !values.f_nan

          omi_pres_pix1 = fltarr(24)
          omi_pres_pix1[*] = !values.f_nan

          if gems_o3_sz[1] eq 174 then begin

            ; get gems variables
            gems_o3 = output.gemsl2o3p.o3[indices[0, ipix], indices[1, ipix], *]
            gems_o3ap = output.gemsl2o3p.O3Apriori[indices[0, ipix], indices[1, ipix], *]
            gems_pres = output.gemsl2o3p.pressure[indices[0, ipix], indices[1, ipix], *]
            gems_pres_pix = output.gemsl2o3p.pressure[$
              indices[0, ipix], indices[1, ipix], *]
            for il = 0, 23 do begin
              gems_pres1 = gems_pres[*, *, il] + gems_pres[*, *, il+1]
              gems_pres_pix1[il] = gems_pres_pix[il] + gems_pres_pix[il+1]
            endfor

            ; get omi variables
            omi_o3 = output.colloc_omi_o3p_on_gems[indices[0, ipix], indices[1, ipix], *]
            omi_o3ap = output.colloc_omi_o3ap_on_gems[indices[0, ipix], indices[1, ipix], *]
            omi_pres = output.colloc_omi_pres_on_gems[indices[0, ipix], indices[1, ipix], *]
            omi_pres_pix = output.colloc_omi_pres_on_gems[$
              indices[0, ipix], indices[1, ipix], *]
            for il = 0, 23 do begin
              omi_pres_pix1[il] = omi_pres_pix[il] + omi_pres_pix[il+1]
            endfor
          endif else if gems_o3_sz[1] eq 24 then begin
            ; get gems variables
            gems_o3 = output.gemsl2o3p.o3[*, indices[0, ipix], indices[1, ipix]]
            gems_o3ap = output.gemsl2o3p.O3Apriori[*, indices[0, ipix], indices[1, ipix]]
            gems_pres = output.gemsl2o3p.pressure[*, indices[0, ipix], indices[1, ipix]]
            gems_pres_pix = output.gemsl2o3p.pressure[$
              *, indices[0, ipix], indices[1, ipix]]
            for il = 0, 23 do begin
              gems_pres1 = gems_pres[il, *, *] + gems_pres[il, *, *]
              gems_pres_pix1[il] = gems_pres_pix[il] + gems_pres_pix[il+1]
            endfor

            ; get omi variables
            omi_o3 = output.colloc_omi_pres_on_gems[*, indices[0, ipix], indices[1, ipix]]
            omi_o3ap = output.colloc_omi_o3ap_on_gems[*, indices[0, ipix], indices[1, ipix]]
            omi_pres = output.colloc_omi_pres_on_gems[*, indices[0, ipix], indices[1, ipix]]
            omi_pres_pix = output.colloc_omi_pres_on_gems[$
              *, indices[0, ipix], indices[1, ipix]]
            for il = 0, 23 do begin
              omi_pres_pix1[il] = omi_pres_pix[il] + omi_pres_pix[il+1]
            endfor
          endif

          p1 = plot(omi_o3, omi_pres_pix1, color='black', /buffer, $
            title='GEMS O3PROFILE ' + datetime_str, $
            xtitle='O3[DU]', $
            ytitle='Pressure[hPa]', $
            symbol='*', $
            font_size=14, $
            dimension=[500, 700], $
            yrange=[1000, 0], /ylog, pos=pos, name='OMI O3 profile')

          p2 = plot(omi_o3ap, omi_pres_pix1, color='grey', /overplot, /buffer, $
            symbol='diamond', $
            yrange=[1000, 0], /ylog, pos=pos, $
            linestyle=1, $
            name='OMI O3 a priori profile')

          p3 = plot(gems_o3, gems_pres_pix1, color='blue', /overplot, /buffer, $
            symbol='+', $
            yrange=[1000, 0], /ylog, pos=pos, name='GEMS O3 profile')

          p4 = plot(gems_o3ap, gems_pres_pix1, color='sky blue', /overplot, /buffer, $
            symbol='triangle', $
            linestyle=1, $
            yrange=[1000, 0], /ylog, pos=pos, name='GEMS O3 a priori profile')

          leg = legend(target=[p1, p2, p3, p4], pos=leg_pos)

          yyyy = string(year, format='(i04)')
          mm = string(month, format='(i02)')
          dd = string(iday, format='(i02)')
          hh = string(hour, format='(i02)')
          mi = string(minute, format='(i02)')

          pngfile = './plot/gems_omi_profile_comparison_' $
            + yyyy + mm + dd + hh + mi + 'z_winlim' +  winlim_str + '340' $
            + '_ci' + string(indices[0, ipix], format='(i03)') $
            + '_cs' + string(indices[1, ipix], format='(i03)') + '.png'

          gemslat1 = output.gemsl2o3p.latitude[indices[0, ipix], indices[1, ipix]]
          gemslon1 = output.gemsl2o3p.longitude[indices[0, ipix], indices[1, ipix]]
          omiecf1 = output.colloc_omi_ecf_on_gems[indices[0, ipix], indices[1, ipix]]

          ;t1 = text(0.45, 0.25, 'Pressure Quality Flag:' + string(info[1], format='(I1)'))
          ;t2 = text(0.45, 0.20, 'Cloud Pressure:' + string(cp[xidx, yidx]))
          t1 = text(0.45, 0.62, $
            'Longitude:' + string(gemslon1, format='(f7.2)'), font_size=12)
          t1 = text(0.45, 0.59, $
            'Latitude:' + string(gemslat1, format='(f7.2)'), font_size=12)
          t2 = text(0.45, 0.56, $
            'Effective Cloud Fraction:' + string(omiecf1, format='(F7.2)'), font_size=12)

          p1.save, pngfile
          p1.close
          ; send image to pc
          scp_dest = 'soodal@164.125.38.179:/home/soodal/works/plot/gems_omi_comp'
          spawn, 'scp -P18742 -p ' + pngfile +  ' ' + scp_dest
        endfor
      endif

    endif
  endfor
endfor

end
