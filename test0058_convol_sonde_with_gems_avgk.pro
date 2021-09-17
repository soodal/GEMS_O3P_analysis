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

sonde_fl = []

timeseries = timegen(start=julday(7, 1, 2021), final=julday(7, 31, 2021), step_size = 1, units='Day')
for itime = 0, n_elements(timeseries)-1 do begin
  caldat, timeseries[itime], month, day, year, hour, minute
  yyyy = string(year, format='(i04)')
  mm = string(month, format='(i02)')
  dd = string(day, format='(i02)')
  ;hh = string(hour, format='(i02)')
  ;mi = string(minute, format='(i02)')

  ;ozonesonde_fn = '/data/ozonesonde_pohang/2021/20210729D_L1.txt'
  ozonesonde_fp = '/data/ozonesonde_pohang/' + yyyy + '/' 
  ozonesonde_fn = + yyyy + mm + dd +'D_L1.txt'

  fl = file_search(ozonesonde_fp + ozonesonde_fn)
  if n_elements(fl) ge 1 and strlen(fl[0]) gt 0 then begin
    print, fl[0]
    sonde_fl = [sonde_fl, fl[0]]
  endif

ENDFOR



sonde_o3prof_avgk_total = fltarr(n_elements(sonde_fl), 24)
sonde_o3prof_avgk_total[*] = !values.f_nan
sonde_o3prof_total = fltarr(n_elements(sonde_fl), 24)
sonde_o3prof_total[*] = !values.f_nan
gems_o3prof_total = fltarr(n_elements(sonde_fl), 24)
gems_o3prof_total[*] = !values.f_nan

sonde_o3prof_avgk_tco_total = fltarr(n_elements(sonde_fl))
sonde_o3prof_avgk_tco_total[*] = !values.f_nan
sonde_o3prof_tco_total = fltarr(n_elements(sonde_fl))
sonde_o3prof_tco_total[*] = !values.f_nan
gems_o3prof_tco_total = fltarr(n_elements(sonde_fl))
gems_o3prof_tco_total[*] = !values.f_nan

for ifile = 0, n_elements(sonde_fl) -1 do begin
  print, sonde_fl[ifile]

  read_ozonesonde_pohang, sonde_fl[ifile], ozonesonde_data

  launch_hour = ozonesonde_data.hour
  launch_minute = ozonesonde_data.minute

  sonde_basename = file_basename(sonde_fl[ifile], '.txt')

  yyyy = strmid(sonde_basename, 0, 4)
  mm = strmid(sonde_basename, 4, 2)
  dd = strmid(sonde_basename, 6, 2)

  gemsfp = '/data/nier_ftp/O3P/V03/' + yyyy + mm + '/' + dd + '/'
  
  ; TODO find the nearest time 
  gemsfn = 'GK2_GEMS_L2_' + yyyy + mm + dd + '_' + string(launch_hour, format='(i02)') + '*.nc'

  gemsfl = file_search(gemsfp + gemsfn)

  ;==================================
  ;(1)  read gems l2 o3p file
  ;==================================  
  if n_elements(gemsfl) ne 0 and strlen(gemsfl[0]) ne 0 then begin
    gemsvar = ds_read_gems_l2_o3p(gemsfl[0])

    ;read_omil2_he5,omifile, nl, nf, nalb, ngas, naer, ominprof, $ ---------------------------------------
      ;omilon, omilat, omisza, omivza, omiaza, omirms, omiavgres, $
       ;omicfrac, omictp, omicldflg, omiai, omiutc, omintp, $
       ;ominw, omisaa, omiexval, ominiter, ominspike, omiglint, omildwt, omisnow, $
       ;omimon, omiday, omiyear, omipix, omiline, omitaod, omitsca, omisaod, omialb, atmos, $
       ;ozprofs, omicol, omitrace,  omifitvar, omiavgk, omicorrel, omicovar, omicontri, omifitspec, $
       ;omisimrad, omiwaves, omiclmrad, omiactrad, omiwf, omisnr, omiring, $ ---------------------------------------
       ; omiptrp = omiptrp, omiztrp = omiztrp, omiorb=omiorb, omijul=omijul, /quiet
       ;nprof = ominprof, /get_avgk, ;

       
    ;omiptrp = fltarr(ominprof) 
    ;omiztrp = fltarr(ominprof)
    ;FOR i = 0 , ominprof -1 do begin
         ;omiptrp(i) = atmos(i, 0, omintp(i))
         ;omiztrp(i) = atmos(i, 1, omintp(i))
    ;ENDFOR

    ;numwin = n_elements( omirms(0,*)) -1

    ; apply basic data filtering
    nl = 24

    sel  = where(gemsvar.SolarZenithAngle le 88 and gemsvar.ProcessingQualityFlags eq 0 and gemsvar.FinalAlgorithmFlags eq 0, nprof, /null) 
    ;if keyword_set(omista) then begin
      ;the_sta = strmid(strcompress(sondes(0).station, /remo), 0, 4)
      ;sel  = where( omiday ne 0 and omisza le 88 and omiexval gt 0  and omiexval lt 100 and omista eq the_sta and omirms le 5, nprof) 
    ;endif

    ;if nprof eq 0 then return 
    ;omijul=omiutc 
    ;omijul(sel)= julday(omimon(sel),omiday(sel),omiyear(sel),omiutc(sel))  
    ;IF size(omicfrac, /n_dim) eq 1 then begin
       ;ncfrac = 1
    ;ENDIF ELSE BEGIN
       ;ncfrac = n_elements(omicfrac(0,*))
    ;ENDELSE

    ;==================================
    ;(2) collocation
    ; find one omipix per sonpix
    ;==================================
    ;ss = where( sondes.year ge min(omiyear(sel)) and sondes.year le max(omiyear(sel)), nsonde)
    ;if nsonde eq 0 then return
    ;sondes = sondes(ss)
    closeidx = search_closest_pixel(gemsvar.longitude, gemsvar.latitude, $
      ozonesonde_data.launch_longitude, ozonesonde_data.launch_latitude, maxlimit=1)

    if closeidx ge 0 then begin
      closeindices = array_indices(gemsvar.longitude, closeidx)

      ;mind     = fltarr(nsonde) & mind(*) = -1
      ;omi      = replicate(define_omison_str(nl, nalb, ncfrac), nsonde)
      ;sonjul   = julday( sondes.mon, sondes.day, sondes.year, sondes.utc)
      ;onejul   = julday(1, 1, 2000,2) - julday(1, 1, 2000,1)
      
      ; find collocation among all omi pixels 

      ;djuls    = abs(sonjul(sidx)     - omijul(sel) )
      ;dutcs    = djuls/onejul
      ;dlats    = abs(sondes(sidx).lat - omilat(sel,4))
      ;dlons    = abs(sondes(sidx).lon - omilon(sel,4))
      ;da = where( djuls le 0.6 and dlats le 1. and dlons le 1.5, nda)
      ;da = where( djuls le 0.6 and dlats le 1.5 and dlons le 3, nda)
      
      ;if nda eq 0 then continue
      
      ;find_overpass, nda, omilon(sel(da)), omilat(sel(da)), sondes(sidx(da)).lon, sondes(sidx).lat, midx, /closest, dis=dis
      ;da = da(midx(0)) 
      ;nda = 1 
      ;oidx       =   sel(da) 
      ;;;;;;;;;;;;;;;;;;;;;;;;;  

      ;nsl        = sondes(sidx).nlvl-1
      nsl        = n_elements(ozonesonde_data.ozon_mpa) - 1
      sonpres    = ozonesonde_data.pres_hpa
      if max(sonpres) ge 2000 then continue
      sono3      = ozonesonde_data.ozon_mpa
      sonT       = ozonesonde_data.temp_degc
      ord = reverse(sort(sonpres))

      sonpres    = sonpres(ord) 
      sont = sont(ord) 
      sono3 = sono3(ord)

      diff = [1, sonpres(0:nsl-2) - sonpres(1:nsl-1) ]
      ss   = where( diff gt 0 and sono3 gt  0, nsl)
      if nsl le 10 then continue

      sonpres = sonpres(ss) 
      sont = sont(ss) 
      sono3 = sono3(ss)

      TWMO, -0.002, sonT, 0, 100000., sonpres*100, sonptrp, sonttrp, sonztrp, 0

      sonztrp    = -16. * alog10(sonptrp / 101325.)
      sonptrp    = sonptrp*0.01

      atmos = fltarr(3, 25)
      atmos(0, *) = gemsvar.pressure[closeindices[0], closeindices[1], *]
      atmos(1, *) = gemsvar.Altitude[closeindices[0], closeindices[1], *]
      atmos(2, *) = gemsvar.Temperature[closeindices[0], closeindices[1], *]

      ;omipres    = reform(atmos(oidx, 0, *))
      ;omiz       = reform(atmos(oidx, 1, *))
      ;ts       = reform(atmos(oidx, 2, *))
      
      omipres    = reform(atmos(0, *))
      omiz       = reform(atmos(1, *))
      ts       = reform(atmos(2, *))

      tmid      = zmid(ts)
      ;ozprof     = reform(ozprofs(oidx, 2, *))
      ;apozprof   = reform(ozprofs(oidx, 0, *))

      ozprof     = gemsvar.O3[closeindices[0], closeindices[1], *]
      apozprof   = gemsvar.O3Apriori[closeindices[0], closeindices[1], *]

      ;ntp        = omintp(oidx)
      tpres = gemsvar.TropopausePressure[closeindices[0], closeindices[1]]

      ;nwave      = ominw(oidx)
      convert_mpa2du, nsl, sonpres, sono3, sonde_o3_total_du, sono3du ;, cols=sono3du ,toc200=t200, toc500=t500, toc750=t750

      omi_avgk = fltarr(nl, nl)
      omi_avgk[*] = !values.f_nan

      omiavgk = gemsvar.AveragingKernel[closeindices[0], closeindices[1], *, *]
      if finite(sonde_o3_total_du) eq 0 then continue
      IF n_elements(omiavgk) ge 24 then begin
        print, 'ifasdfafdasd'
        avgk       = reform(omiavgk)

        ;convol_sonde_with_gems_ak, sonpres, sono3du, omipres, apozprof, ozprof, avgk, ntp, $ ; inputs
        convol_sonde_with_gems_ak, sonpres, sono3du, omipres, apozprof, ozprof, avgk, -1, $ ; inputs
               sonstlvl, csono3du0, csonsco, csontco, csonsco200, csontco200, $
               omisco, omitco, omisco200, omitco200, $
               omiasco, omiatco, omiasco200, omiatco200, psontop, $
               convol=1, cconvol = 1, use_spline=1, append_ret = 1, print_result = 1, midpres=200, insert_ptrp=tpres
        for ip = 0 , nl-1  do begin
          ;omi(sidx).avgk(*,ip) = avgk(*,ip)
          omi_avgk[*,ip] = avgk[*,ip]
        endfor
        omi_avgk = avgk
      endif else begin
            csontco = 0 & csontco200 = 0 & csonsco = 0 & csonsco200 = 0 
            csono3du0 = fltarr(nl)
      endelse
        ;convol_sonde_with_gems_ak, sonpres, sono3du, omipres, apozprof, ozprof, avgk, ntp, $ ; inputs
        convol_sonde_with_gems_ak, sonpres, sono3du, omipres, apozprof, ozprof, avgk, -1, $ ; inputs
               sonstlvl, sono3du0, sonsco, sontco, sonsco200, sontco200, $
               omisco,  omitco,  omisco200, omitco200, $
               omiasco, omiatco, omiasco200, omiatco200, psontop, $
               convol=0, cconvol = 0, use_spline=1, append_ret = 1, print_result = 0, midpres=200,  $
              insert_ptrp = sonptrp

      ;omi(sidx).alt       = omiz  
      omi_alt       = omiz
      ;omi(sidx).tmid       = tmid
      omi_tmid       = tmid
      ;omi(sidx).pres      = omipres  
      omi_pres      = omipres  
      ;omi(sidx).aozprofer  = transpose(reform(ozprofs(oidx,1,*)))
      omi_aozprofer  = transpose(reform(gemsvar.O3AprioriError[closeindices[0], closeindices[1], *]))
      ;omi(sidx).solutioner = transpose(reform(ozprofs(oidx,3,*)))
      omi_solutioner = transpose(reform(gemsvar.O3SolutionError[closeindices[0], closeindices[1], *]))
      ;omi(sidx).randomer   = transpose(reform(ozprofs(oidx,4,*)))
      omi_randomer   = transpose(reform(gemsvar.O3RandomNoiseError[closeindices[0], closeindices[1], *]))

      ; ozone profile
      ;omi(sidx).ozprof     = ozprof   
      omi_ozprof     = ozprof   

      ; ozone a-priori profile
      ;omi(sidx).aozprof    = apozprof        
      omi_aozprof    = apozprof        

      ; ozonesonde profile
      ;omi(sidx).sonprof    = sono3du0  
      omi_sonprof    = sono3du

      ; ozone_sonde profile 
      ;omi(sidx).csonprof   = csono3du0  
      omi_csonprof   = csono3du0  

      ; tropopause height
      ;ztrp = omiztrp(oidx)
      ztrp = tpres

      calc_cum ,0, omiz, apozprof, ztrp, ztrp-3, out
      omiatco = out
      bad = where( finite(sono3du0) eq 0)
      sono3du0(bad) = ozprof(bad)
      calc_cum ,0, omiz, sono3du0, ztrp, ztrp-5, out
      sontco = out
      ;omi(sidx).tco    = [omitco, omiatco, sontco, csontco]
      ;omi(sidx).sco    = [omisco, omiasco, sonsco, csonsco]
      ;omi(sidx).tco200 = [omitco200, omiatco200, sontco200,csontco200]
      ;omi(Sidx).sco200 = [omisco200, omiasco200, sonsco200,csonsco200]
      ;omi(sidx).sonztrp = sonztrp
      ;omi(Sidx).sonptrp = sonptrp
      ;omi(Sidx).psontop = psontop

      omi_tco    = [omitco, omiatco, sontco, csontco]
      omi_sco    = [omisco, omiasco, sonsco, csonsco]
      omi_tco200 = [omitco200, omiatco200, sontco200, csontco200]
      omi_sco200 = [omisco200, omiasco200, sonsco200, csonsco200]
      omi_sonztrp = sonztrp
      omi_sonptrp = sonptrp
      omi_psontop = psontop
     
      ;omi(sidx).dlat = dlats(da)
      ;omi(sidx).dlon = dlons(da)
      ;omi(Sidx).dutc = dutcs(da)
      ;omi(Sidx).djul = djuls(da)
      ;omi(sidx).nwave = nwave
      ;omi(sidx).dis   = dis
      ;mind (sidx) = oidx

      ;omi_dlat = dlats(da)
      ;omi_dlon = dlons(da)
      ;omi_dutc = dutcs(da)
      ;omi_djul = djuls(da)
      ;omi_nwave = nwave
      ;omi_dis   = dis
      ;mind = oidx

      ;if abs(omi(sidx).tco(2) - omi(sidx).tco(3)) ge 9 then stop
      ;print , omi(sidx).tco, omi(sidx).psontop, omiptrp(da), omicfrac(da), format='(100f7.2)'
      ;if omi(sidx).tco(2) le 8 then stop

      ; plotting
      ;for ifile = 0, n_elements(sonde_fl) do begin
      pres = reform(gemsvar.pressure[closeindices[0], closeindices[1], *])
      p1 = plot(ozprof, pres, /buffer, /overplot, /ylog, dim=[500, 700]);color=[0, 0, 0], linestyle=0, name='GEMS O3P')
      p1.title = sonde_basename
      p1.color = [0, 0, 0]
      p1.yrange= [1000, 1]
      p1.linestyle = 0
      p1.symbol= '+'
      p1.name = 'GEMS O3P'

      p2 = plot(omi_aozprof, pres, /buffer, /overplot);, $ color='#e4a000', linestyle=1, name='Ozonesonde with GEMS O3P AVGK', /overplot)
      p2.color = [228, 160, 0]
      p2.linestyle = 1
      p2.symbol= 's'
      p2.name = 'GEMS O3P A-priori'

      p3 = plot(omi_csonprof, pres, /buffer, /overplot);, $ color='#e4a000', linestyle=1, name='Ozonesonde with GEMS O3P AVGK', /overplot)
      p3.color = [86, 180, 232]
      p3.linestyle = 2
      p3.symbol= '*'
      p3.name = 'Ozonesonde with GEMS O3P AVGK'

      p4 = plot(sono3du0, pres, /buffer, /overplot);, $ color='#56b4e8', linestyle=2, name='Ozonesonde', /overplot)
      p4.color = [0, 159, 115]
      p4.linestyle = 3
      p4.symbol= 'D'
      p4.name = 'Ozonesonde'

      leg = legend(target=[p1, p2, p3, p4])

      p1.save, '/data/soodal_data/collocated_gemsl2o3p_ozonesonde/' + $
        sonde_basename + '.png'
      p1.close

      ;endfor


      ; ------------
      sonde_o3prof_avgk_total[ifile, *] = omi_csonprof
      sonde_o3prof_total[ifile, *] = sono3du0
      gems_o3prof_total[ifile, *] = ozprof

      sonde_o3prof_avgk_tco_total[ifile] = csontco
      sonde_o3prof_tco_total[ifile] = sontco
      gems_o3prof_tco_total[ifile] = omitco

    ENDif
  ENDif
endfor

  sonde_finiteidx = finite(sonde_o3prof_avgk_total)
  gems_finiteidx = finite(gems_o3prof_total)

  both_finiteidx = where(sonde_finiteidx and gems_finiteidx, /null)



  print, correlate(sonde_o3prof_avgk_total[both_finiteidx], gems_o3prof_total[both_finiteidx])
  
  result = linfit(sonde_o3prof_avgk_total[both_finiteidx], gems_o3prof_total[both_finiteidx], $
    yfit=yfit)

  plot_gems_validation, gems_o3prof_total[both_finiteidx], sonde_o3prof_avgk_total[both_finiteidx], $
    filename='/data/soodal_data/collocated_gemsl2o3p_ozonesonde/gems_l2_o3p_val_with_ozonesonde_' + $
      '202107_pohang.png', $
    xtitle='GEMS O3', $
    ytitle='Ozonesonde O3', $
    cblim=[0, 5], $
    range=[0, 60], $
    delta=2.0
  stop
  da = where( mind ne -1, nprof)  &  if nprof eq 0 then return
  omi(da).sta  = sondes(da).station
  omi(da).sonlat = sondes(da).lat
  omi(da).sonz0  = sondes(da).z0
  omi(da).sontype   = sondes(da).type
  omi(da).sonsource = sondes(da).source
  omi(da).nson     = nsonde
  omi(da).nprof    = nprof
  omi(da).sonjul   = sonjul(da)
  omi(da).sonutc   = sondes(da).utc
  omi(da).corr     = sondes(da).corr
  omi(da).ok       = sondes(da).ok
  omi(da).day  = omiday(mind(da))   
  omi(da).mon  = omimon(mind(da))
  omi(da).year = omiyear(mind(da))
  omi(da).utc  = omiutc(mind(da))
  omi(da).lat  = omilat(mind(da),4)
  omi(da).lon  = omilon(mind(da),4)
  omi(da).sza  = omisza(mind(da))
  omi(da).vza  = omivza(mind(da))
  omi(da).pix  = omipix(mind(da))
  omi(da).orb  = omiorb(mind(da))
  omi(da).line = omiline(mind(da))
  omi(da).ptrp = omiptrp(mind(da))
  omi(da).ztrp = omiztrp(mind(da))

  omi(da).ctp  = omictp(mind(da))
  omi(da).ai   = omiai(mind(da))
  omi(da).rms    = omirms(mind(da),0)
  omi(da).rmsUv1    = omirms(mind(da),1)
  omi(da).spike    = ominspike(mind(da))
  omi(da).saa      = omisaa(mind(da))
  ;omi(da).xflag      = omixflag(mind(da))
  IF numwin ge 2 then omi(da).rmsUV2    = omirms(mind(da),2)
  omi(da).etco      = omicol(mind(da), 2,1)
  omi(da).esco      = omicol(mind(da), 1,1)
  omi(da).avgres = omiavgres(mind(da),0)
  omi(da).avgresuv1 = omiavgres(mind(da),1)

  IF numwin ge 2 then omi(da).avgresuv2 = omiavgres(mind(da),2)
  omi(da).alb(0:nalb-1)    = transpose(omialb(mind(da),0:nalb-1))
  omi(da).cfrac(0:ncfrac-1)  = transpose(omicfrac(mind(da),0:ncfrac-1)   )
  files = omifile(mind(da))
  omi(da).jul    = omijul(mind(da))
  omi = omi(da)
  omi = omi(sort(omi.sonjul))
END
