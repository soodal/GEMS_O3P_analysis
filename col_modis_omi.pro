;===================================================================================================
; PROGRAM NAME:
;  col_modis_omi
;
; PURPOSE:
;  mapping collocated OMI TO3 and MODIS CTP
;
; PROCESS:
;
; REFERENCE:
;
; REQUIRED:
;  read_modis_myd06.pro 
;  read_l2_v8.pro
;  convert_bit2flag.pro
;  h5read.pro
;
; REQUIRED 
;
;
; DEVELOPER:
;  Daegeun Shin (geun)
;  Kanghyun Back
;
;  Satellite Remote Sensing Laboratory
;  Division of Earth Environmental System
;  College of Natural Science
;  Pusan National University
;                  Tel: +82-51-510-2172
;                  E-mail: daegeun@pusan.ac.kr
;
; REVISION HISTORY:
;
;  Copyright (C) 2017 Satellite Remote Sensing Laboratory, PNU
;  All right reserved.
;
;===================================================================================================

;+---------------------------------------------------------------------------+
; Main Procedure 
;+---------------------------------------------------------------------------+
;PRO col_modis_omi

  ; user control

start = systime(2)

dates = '2007'+string( findgen(30) + 401 ,f='(I4.4)')
FOR id=1,1 DO BEGIN 
  date=dates[id]
  call_xdr=1  ; 0: run without xdr,  1: xdr generation, 2: use xdr

  year=STRMID(date,0,4) 
  mon=STRMID(date,4,2)
  day=STRMID(date,6,2)

  limit=[-20,-20,5,10]

  xdrfile='col_modis_omi_'+date+'.xdr'
  IF call_xdr LE 1 THEN BEGIN
    omipath='/home/Data/OMI/2_OML2TO/'
    omifiles=FILE_SEARCH(omipath+year+'/'+'OMI-Aura*'+year+'m'+mon+day+'*.he5', count=nfiles)
    omilons=[]  &  omilats=[]  & omitozs=[]  & omitimes=[]  & omipixs=[]  & omilines=[]
    pixlons=[]  &  pixlats=[]  
    omir331s=[ ] & omipclds = [ ] & omicfs = [ ]  & omiblozs = [ ]
    omiqfs = [ ]
    FOR ifile=0,nfiles-1 DO BEGIN
      omifile=omifiles[ifile]
      READ_L2_V8, omifile, omires
      tmp=WHERE(omires.lon GE limit[1] AND omires.lon LE limit[3] AND $
                omires.lat GE limit[0] AND omires.lat LE limit[2], ntmp)
      IF ntmp NE 0 THEN BEGIN 
        print, '  READ omifile : ', omifile, ntmp
        header='OMI-Aura_L2-OMPIXCOR_'
        time=STRMID(omifile,46,14,/rev)
        pixfile=FILE_SEARCH('/home/Data/OMI/2_OML2PIXCOR/'+year+'/'+header+time+'*.he5', count=npixfile)
        IF npixfile NE 1 THEN BEGIN
          PRINT, '  No Pixfile'
          STOP
        ENDIF ; npixfile
        dim=SIZE(omires.lon,/dim)
        nx=dim[0]
        ny=dim[1]

        READ_L2_PIXCOR, pixfile, pixres
;        omilon=omires.lon[tmp]
;        omilat=omires.lat[tmp]
;        omitoz=omires.toz[tmp]
;        omitime=omires.date[tmp]
;        omipix=omires.xtracks[tmp]
;        omiline=omires.lines[tmp]
        pixlon=REFORM(pixres.corlon,[nx*ny,4])
        pixlat=REFORM(pixres.corlat,[nx*ny,4])     
        pixlon=pixlon[tmp,*]
        pixlat=pixlat[tmp,*]
;        omitime=pixres.jtime[tmp]
        omilons  =  [omilons, omires.lon[tmp] ]
        omilats  =  [omilats, omires.lat[tmp] ]
        omitozs  =  [omitozs, omires.toz[tmp] ]
        omitimes =  [omitimes,pixres.jtime[tmp]]
        pixlats  =  [pixlats, pixlat]
        pixlons  =  [pixlons, pixlon]

        omipixs  =  [omipixs,omires.xtracks[tmp]]
        omilines =  [omilines,omires.lines[tmp]]

        omicfs =    [omicfs,   omires.cf[tmp]   ]
        omiblozs =  [omiblozs,   omires.blo3[tmp]   ]
        omipclds =  [omipclds, omires.pc[tmp]     ]
        omir331s =  [omir331s, omires.ref331[tmp] ]
        omiqfs   =  [omiqfs,   omires.qf[tmp]     ]
      ENDIF  ; ntmp
    ENDFOR  ; ifile

    ; modis filename
    julmod=STRING(JULDAY(mon,day,year)-JULDAY('01','01',year)+1,F='(I03)')
    modpath='/home/Data2/5_MODIS/MYD06/'
    modfiles=FILE_SEARCH(modpath+'MYD06_L2.A'+year+julmod+'*.hdf', count=nmodfile)

    mod03path='/home/Data2/5_MODIS/MYD03/'
    mod03files=FILE_SEARCH(mod03path+'MYD03.A'+year+julmod+'*.hdf', count=nmod03file)

    IF nmodfile EQ 0 or nmod03file eq 0 THEN BEGIN
      PRINT, '  No Modfile'
      STOP
    ENDIF  ; nmodfile

    IF nmodfile ne nmod03file then begin
      times06 =  strmid(modfiles,strlen(modpath)+18,4)
      times03 =  strmid(mod03files,strlen(mod03path)+15,4)

      loc06 = cmset_op(times06,'AND',times03,/index) 
      loc03 = cmset_op(times03,'AND',times06,/index) 

      modfiles = modfiles[loc06]
      mod03files = mod03files[loc03]

      nmodfile = n_elements(modfiles)
    ENDIF

    ; initialize modis variables

    modlons=[]  
    modlats=[]  
    modtimes=[] 
    mod03lons=[]  
    mod03lats=[]  
    mod03times=[] 
    modcfs=[]
    modctps=[]
    modcphases=[]
    mod03cots = [ ]

    ; read n modis files
    FOR ifile=0,nmodfile-1 DO BEGIN
      modfile=modfiles[ifile]
      mod03file=mod03files[ifile]
      modres=READ_MODIS_MYD06(modfile)
      mod03res=READ_MODIS_MOD03(mod03file)

      mlon  = modres.lon
      mlat = modres.lat
      hlon  = mod03res.lon 
      hlat = mod03res.lat

      tmp=WHERE(mlon GE limit[1] AND mlon LE limit[3] AND $
                mlat GE limit[0] AND mlat LE limit[2], ntmp)

      tmp1=WHERE(hlon GE limit[1] AND hlon LE limit[3] AND $
                 hlat GE limit[0] AND hlat LE limit[2], ntmp1)

      IF ntmp NE 0 THEN BEGIN
        print, '  READ modfile : ', modfile, ntmp
        modlon=modres.lon[tmp]
        modlat=modres.lat[tmp]
        modtime=modres.jtime[tmp]
        modctp=modres.ctp[tmp]
        modcf=modres.cf[tmp]
        modcphase=modres.cphase[tmp]     
 
        modlons       = [modlons, modlon]
        modlats       = [modlats, modlat]
        modtimes      = [modtimes, modtime]
        modctps       = [modctps, modctp]
        modcfs        = [modcfs, modcf]
        modcphases    = [modcphases, modcphase]
      ENDIF  ; ntmp

      IF ntmp1 NE 0 THEN BEGIN
        mod03lons       = [mod03lons, mod03res.lon[tmp1]]
        mod03lats       = [mod03lats, mod03res.lat[tmp1]]
        mod03times       = [mod03times, mod03res.jtime[tmp1]]
        mod03cots       = [mod03cots, modres.cot[tmp1]]
      ENDIF
    ENDFOR  ; ifile

    omivars={lon:omilons, lat:omilats, time:omitimes, toz:omitozs,$
             pixlat:pixlats, pixlon:pixlons, pixs:omipixs, lines:omilines,$
             ref331:omir331s, cf:omicfs, pcld:omipclds, qf:omiqfs, bloz:omiblozs}
    modvars={lon:modlons, lat:modlats, time:modtimes, ctp:modctps, cf:modcfs, cphase:modcphases,$
             hlon:mod03lons,hlat:mod03lats,htime:mod03times,cot:mod03cots}
    IF call_xdr EQ 1 THEN BEGIN
      SAVE, file=xdrfile, omivars, modvars, /xdr
    ENDIF ; call_xdr

  ENDIF ELSE BEGIN  ; call_xdr

    RESTORE, xdrfile

  ENDELSE  ; call_xdr

;  ; collocation
  npo=N_ELEMENTS(omivars.lon)
  time_bo=120*60  ; (unit: second)
  dist_bo=20   ; distance limit for nearest modis center from omi center (unit : km)
  omictp=FLTARR(npo)
  omicf=FLTARR(npo)
  omicphase=FLTARR(npo)
  omicot = fltarr(npo)

  omitime = omivars.time

;5km modis value
  mlon = modvars.lon &  mlat = modvars.lat
  mtime = modvars.time
;1km 
  hlon=modvars.hlon &  hlat=modvars.hlat
  htime=modvars.htime &  hcot=modvars.cot
 
  FOR ipo=0,npo-1 DO BEGIN
    print, 'ipo: ', ipo
    print,  '  '

    lon_min=MIN(omivars.pixlon[ipo,*]) & lon_max=MAX(omivars.pixlon[ipo,*])
    lat_min=MIN(omivars.pixlat[ipo,*]) & lat_max=MAX(omivars.pixlat[ipo,*])

  ;5km spatial resolution
    tmp=WHERE(mlon GE lon_min AND mlon LE lon_max AND $
              mlat GE lat_min AND mlat LE lat_max AND $
              ABS(mtime-omitime[ipo]) LE time_bo, ntmp) 
    IF ntmp NE 0 THEN BEGIN
      ; ctp collocation
      ctps=modvars.ctp[tmp]
      val=WHERE(ctps NE -999, nval)
      avgctp=(nval NE 0) ? MEAN(ctps[val]) : -999
      omictp[ipo]=avgctp     

      ; cf collocation
      cfs=modvars.cf[tmp]
      val=WHERE(cfs NE -999, nval)
      avgcf=(nval NE 0) ? MEAN(cfs[val]) : -999
      omicf[ipo]=avgcf     

      ; cf collocation
      cphases=modvars.cphase[tmp]
      val=WHERE(cphases NE -999, nval)
      IF nval NE 0 THEN BEGIN
        dists=ABS(omivars.lon[ipo]-modvars.lon[tmp[val]])+ABS(omivars.lat[ipo]-modvars.lat[tmp[val]])
        mindist=WHERE(dists EQ MIN(dists), nmindist)
        mindist=mindist[0]
        mdist=MAP_2POINTS(omivars.lon[ipo], omivars.lat[ipo], $
                          modvars.lon[tmp[val[mindist]]], modvars.lat[tmp[val[mindist]]], /meters)/1000.
        sel_cphase=(mdist LE dist_bo) ? cphases[val[mindist]] : -999
        
      ENDIF ELSE BEGIN  ; nval
        sel_cphase=-999
      ENDELSE
      omicphase[ipo]=sel_cphase

    ENDIF ELSE BEGIN ; ntmp
      omictp[ipo] =-999
      omicf[ipo]  =-999
      omicphase[ipo] =-999
    ENDELSE ; ntmp

; IF keyword_Set(include_COT) THEN BEGIN 
  ;1km spatial resolution
;    tmp1=WHERE(modvars.hlon GE lon_min AND modvars.hlon LE lon_max AND $
;               modvars.hlat GE lat_min AND modvars.hlat LE lat_max AND $
;              ABS(modvars.htime-omivars.time[ipo]) LE time_bo, ntmp1) 
    tmp1=WHERE(hlon GE lon_min AND hlon LE lon_max AND $
               hlat GE lat_min AND hlat LE lat_max AND $
              ABS(htime-omivars.time[ipo]) LE time_bo, ntmp1) 
    IF ntmp1 ne 0 then begin 
      ; ctp collocation
      cots=modvars.cot[tmp1]
      val=WHERE(cots NE -999, nval)
      avgcot=(nval NE 0) ? MEAN(cots[val]) : -999
      omicot[ipo]=avgcot    
    ENDIF ELSE BEGIN
      omicot[ipo]= -999
    ENDELSE
;  ENDIF
    ; check point
    ;CALDAT,(JULDAY(1,1,1993,0,0,0)*86400.0d0+ABS(omivars.time[ipo]-modvars.time[tmp]))/86400.0d0,mon,day,year,hour,min,sec
    ;time_diff='h'+STRING(hour,F='(I02)')+'-m'+ STRING(min,F='(I02)')+'-s'+STRING(sec,F='(I02)')
    ;print, ipo, omicphase[ipo], mdist
    ;print, ipo, time_diff[0]

  ENDFOR  ; ipix
  res={omi:omivars,$
       omictp:omictp, omicf:omicf, omicphase:omicphase,omicot:omicot}

  file = 'final_col_modis_omi_'+date+'.xdr'
  save,file=file,/xdr,res
  print, "saving : ",file
ENDFOR

final = systime(2)

time = (final - start)/60.
print, 'Time: ' , time

END
