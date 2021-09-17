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
pro load_gems_sonde, sondes, gemsfile,  gems,  nprof, show=show, colist=colist,gemssta=gemssta

  IF keyword_set(colist) THEN begin 
     ncolist = n_elements(colist.year)
     mind = fltarr(ncolist) & mind(*) = -999
     sonjul = julday( sondes.mon, sondes.day, sondes.year, sondes.utc)
     FOR i = 0 , ncolist -1 do begin
         da = where( strpos( strmid(sondes.station, 0, 5), colist(i).sta) ne -1 and  $
                    abs(sonjul - colist(i).sonjul) le 0.4, nda)
         if nda ge 0 then begin
            da = da(0) & nda = 1
         endif
         if nda ne 1 then begin 
              print, i, nda
               stop
         endif
         mind(i) = da
     ENDFOR
     mind = mind(where(mind ge 0, nsonde))
     sondes = sondes(mind)
    ; print,'nsonde:', nsonde
  ENDIF
  ;==================================
  ;(1)  read omil2file
  ;==================================  
  IF strpos(omifile(0), '.out') ne -1 then begin 
  read_omil2_file,omifile, nl, nf, nalb, ngas, naer, $
     ominprof, omilon, omilat, omisza, omivza, omiaza, omirms, omiavgres, $
     omicfrac, omictp, omicldflg, omiai, omiutc,  $
     omintp, ominw, omisaa, omiexval, ominiter, ominspike, omiglint, omildwt, omisnow, $
     omimon, omiday, omiyear, omipix, omiline, omitaod, omitsca, omisaod, omialb, atmos, $
     ozprofs, omicol, omitrace,  omifitvar, omiavgk, omicorrel, omicovar, omicontri, omifitspec, $
     omisimrad, omiwaves, omiclmrad, omiactrad, omiwf, omisnr, omiring, /get_avgk, $
     OMIptrp = omiptrp, omiztrp = omiztrp, omiorb=omiorb,omijul=omijul,  /quiet
     nprof = ominprof
  ENDIF ELSE IF strpos(omifile(0), '.he5') ne -1 then begin
  read_omil2_he5,omifile, nl, nf, nalb, ngas, naer, $
     ominprof, omilon, omilat, omisza, omivza, omiaza, omirms, omiavgres, $
     omicfrac, omictp, omicldflg, omiai, omiutc,  $
     omintp, ominw, omisaa, omiexval, ominiter, ominspike, omiglint, omildwt, omisnow, $
     omimon, omiday, omiyear, omipix, omiline, omitaod, omitsca, omisaod, omialb, atmos, $
     ozprofs, omicol, omitrace,  omifitvar, omiavgk, omicorrel, omicovar, omicontri, omifitspec, $
     omisimrad, omiwaves, omiclmrad, omiactrad, omiwf, omisnr, omiring, /get_avgk, $
     omiptrp = omiptrp, omiztrp = omiztrp, omiorb=omiorb, omijul=omijul, /quiet
     nprof = ominprof
  ENDIF ELSE IF strpos(omifile(0), '.xdr') ne -1 THEN BEGIN
    restore, omifile 
    xprofoz = 0
    IF xprofoz eq 1 then begin
     nalb = 3
     omialb=fltarr(ominprof, 3)
     omijul= julday(omimon,omiday,omiyear,omiutc)  
     omivza = omiday
     omiptrp = omiday & omiztrp = omiday 
     FOR i = 0 , ominprof -1 do begin
       omiptrp = atmos(i, 0, omintp(i))
       omiztrp = atmos(i, 1, omintp(i))
     ENDFOR
   ENDIF
  ENDIF
  IF ominprof eq 0 then return
     
  omiptrp = fltarr(ominprof) & omiztrp = fltarr(ominprof)
  FOR i = 0 , ominprof -1 do begin
       omiptrp(i) = atmos(i, 0, omintp(i))
       omiztrp(i) = atmos(i, 1, omintp(i))
  ENDFOR

  numwin = n_elements( omirms(0,*)) -1
  ; apply basic data filtering
  sel  = where( omiday ne 0 and omisza le 88 and omiexval gt 0 and omiexval lt 100  , nprof) 
  if keyword_set(omista) then begin
    the_sta = strmid(strcompress(sondes(0).station, /remo), 0, 4)
    sel  = where( omiday ne 0 and omisza le 88 and omiexval gt 0  and omiexval lt 100 and   omista eq the_sta and omirms le 5, nprof) 
  endif
  if nprof eq 0 then return 
  omijul=omiutc 
  omijul(sel)= julday(omimon(sel),omiday(sel),omiyear(sel),omiutc(sel))  
  IF size(omicfrac, /n_dim) eq 1 then begin
     ncfrac = 1
  ENDIF ELSE BEGIN
     ncfrac = n_elements(omicfrac(0,*))
  ENDELSE
  ;==================================
  ;(2) collocation
  ; find one omipix per sonpix
  ;==================================
  ss = where( sondes.year ge min(omiyear(sel)) and sondes.year le max(omiyear(sel)), nsonde)
  if nsonde eq 0 then return
  sondes = sondes(ss)
  mind     = fltarr(nsonde) & mind(*) = -1
  omi      = replicate(define_omison_str(nl, nalb, ncfrac), nsonde)
  sonjul   = julday( sondes.mon, sondes.day, sondes.year, sondes.utc)
  onejul   = julday(1, 1, 2000,2) - julday(1, 1, 2000,1)
  
  FOR sidx = 0, nsonde -1 Do BeGIN
      ; find collocation among all omi pixels 

      djuls    = abs(sonjul(sidx)     - omijul(sel) )
      dutcs    = djuls/onejul
      dlats    = abs(sondes(sidx).lat - omilat(sel,4))
      dlons    = abs(sondes(sidx).lon - omilon(sel,4))
      da = where( djuls le 0.6 and dlats le 1. and dlons le 1.5, nda)
      da = where( djuls le 0.6 and dlats le 1.5 and dlons le 3, nda)
      
      if nda eq 0 then continue
      
      find_overpass, nda, omilon(sel(da)), omilat(sel(da)), sondes(sidx(da)).lon, sondes(sidx).lat, midx, /closest, dis=dis
      da = da(midx(0)) & nda = 1 
      oidx       =   sel(da) 
      ;;;;;;;;;;;;;;;;;;;;;;;;;  

      nsl        = sondes(sidx).nlvl-1
      sonpres    = sondes(sidx).p(0:nsl-1) & if max(sonpres) ge 2000 then continue
      sono3      = sondes(sidx).o3(0:nsl-1)
      sonT       = sondes(sidx).t(0:nsl-1)
      ord = reverse(sort(sonpres))
      sonpres    = sonpres(ord) & sont = sont(ord) & sono3=sono3(ord)
      diff = [1, sonpres(0:nsl-2) - sonpres(1:nsl-1) ]
      ss   = where( diff gt 0 and sono3 gt  0, nsl)
      if nsl le 10 then continue
      sonpres = sonpres(ss) &  sont = sont(ss) & sono3 = sono3(ss)
      TWMO, -0.002, sonT,0, 100000., sonpres*100, sonptrp, sonttrp, sonztrp,0
      sonztrp    = -16. * alog10(sonptrp / 101325.)
      sonptrp    = sonptrp*0.01

      omipres    = reform(atmos(oidx, 0, *))
      omiz       = reform(atmos(oidx, 1, *))
      ts       = reform(atmos(oidx, 2, *))
      tmid      = zmid(ts)
      ozprof     = reform(ozprofs(oidx, 2, *))
      apozprof   = reform(ozprofs(oidx, 0, *))
      ntp        = omintp(oidx)
      nwave      = ominw(oidx)
      get_intoz, nsl, sonpres, sono3, aintoz , cols=sono3du ,toc200=t200, toc500=t500, toc750=t750
      if finite(aintoz) eq 0 then continue
      IF n_elements(omiavgk) ge nl then begin

      avgk       = reform(omiavgk(oidx, *,*))
      convol_sonde_with_omiak,sonpres,  sono3du, omipres, apozprof, ozprof, avgk, ntp,  $
             sonstlvl,csono3du0, csonsco, csontco, csonsco200, csontco200, $
             omisco, omitco, omisco200, omitco200, $
             omiasco, omiatco, omiasco200, omiatco200, psontop, $
             convol=1, cconvol = 1, use_spline=1, append_ret = 1, print_result = 0, midpres=200
      for ip = 0 , nl-1  do omi(sidx).avgk(*,ip) = avgk(*,ip)

      endif else begin
            csontco = 0 & csontco200 = 0 & csonsco = 0 & csonsco200 = 0 
            csono3du0 = fltarr(nl)
      endelse
      convol_sonde_with_omiak,sonpres,  sono3du, omipres, apozprof, ozprof, avgk, ntp,  $
             sonstlvl, sono3du0, sonsco, sontco, sonsco200, sontco200, $
             omisco,  omitco,  omisco200, omitco200, $
             omiasco, omiatco, omiasco200, omiatco200, psontop, $
             convol=0, cconvol = 0, use_spline=1, append_ret = 1, print_result = 0, midpres=200 ;, $
            ;insert_ptrp = sonptrp

      omi(sidx).alt       = omiz  
      omi(sidx).tmid       = tmid
      omi(sidx).pres      = omipres  
      omi(sidx).aozprofer  = transpose(reform(ozprofs(oidx,1,*)))
      omi(sidx).solutioner = transpose(reform(ozprofs(oidx,3,*)))
      omi(sidx).randomer   = transpose(reform(ozprofs(oidx,4,*)))

      omi(sidx).ozprof     = ozprof   
      omi(sidx).aozprof    = apozprof        
      omi(sidx).sonprof    = sono3du0  
      omi(sidx).csonprof   = csono3du0  
;      ztrp = omiztrp(oidx)
;      calc_cum ,0, omiz, apozprof, ztrp, ztrp-3, out
;      omiatco = out
;      bad = where( finite(sono3du0) eq 0)
;      sono3du0(bad) = ozprof(bad)
;      calc_cum ,0, omiz, sono3du0, ztrp, ztrp-5, out
;      sontco = out
      omi(sidx).tco    = [omitco, omiatco, sontco, csontco]
      omi(sidx).sco    = [omisco, omiasco, sonsco, csonsco]
      omi(sidx).tco200 = [omitco200, omiatco200, sontco200,csontco200]
      omi(Sidx).sco200 = [omisco200, omiasco200, sonsco200,csonsco200]
      omi(sidx).sonztrp = sonztrp
      omi(Sidx).sonptrp = sonptrp
      omi(Sidx).psontop = psontop
     
      omi(sidx).dlat = dlats(da)
      omi(sidx).dlon = dlons(da)
      omi(Sidx).dutc = dutcs(da)
      omi(Sidx).djul = djuls(da)
      omi(sidx).nwave = nwave
      omi(sidx).dis   = dis
      mind (sidx) = oidx
    ;  if abs(omi(sidx).tco(2) - omi(sidx).tco(3)) ge 9 then stop
     ;print , omi(sidx).tco, omi(sidx).psontop, omiptrp(da), omicfrac(da), format='(100f7.2)'
;      if omi(sidx).tco(2) le 8 then stop
     ENDFOR
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
