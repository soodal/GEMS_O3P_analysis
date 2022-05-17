;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 2013 4 29  by Jbak
; (1) Main read_fnl_nc, date, utc, fnl, fnl_new
; (2) sub procedure
;     rearray_fnl, fnl, fnl_new
;     function : endday(year, mon)
;     read_fnl_nc, date, utc, fnl, fnl_new, get_rearray=get_rearray
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;+---------------------------------------------------------------------------+
; Sub Function 1 : endday(year, mon)
;+---------------------------------------------------------------------------+
function enday, year, mon
  yr=FIX(year)
  mm=FIX(mon)-1
  IF (((yr mod 4) eq 0) and ((yr mod 100) ne 0)) or ((yr mod 400) eq 0) THEN BEGIN
    out_date=[31,29,31,30,31,30,31,31,30,31,30,31]
  ENDIF ELSE BEGIN
    out_date=[31,28,31,30,31,30,31,31,30,31,30,31]
  ENDELSE
  return, out_date[mm]
end


;+---------------------------------------------------------------------------+
; Sub Procedure 1 : rearray_fnl 
;+---------------------------------------------------------------------------+
pro sub1_rearray_fnl , fnl, fnl_new

  nutc  = size(fnl, /di)
  nlon = 360
  nlat = 180
  nl   = 26
  var1 = fltarr(nlon, nlat)
  var2 = fltarr(nlon, nlat, nl)
  str  = {sp:var1, tp:var1,st:var1, temp:var2, lon:fltarr(nlon), lat:fltarr(nlat), pres:fltarr(nl)}
  fnl_new   = replicate(str,nutc(0))

  lat0  = fnl(0).lat ;[90, -90]
  lon0  = fnl(0).lon ;[0, 359 ]
  pres0 = fnl(0).pres ;[10, 1000]


  FOR k = 0 , nutc(0)-1 do begin
   sp0=fnl(k).sp
   tp0=fnl(k).tp
   st0=fnl(k).st
   temp0=fnl(k).temp

;(1) rearray latitude from [90, -90] to [-90, 90]
   ord   = sort(lat0)
   lat   = lat0(ord)
   sp    = sp0(*,ord)
   tp    = tp0(*,ord)
   st    = st0(*,ord)
   temp  = temp0(*,ord,*)

  ;(2) rearray latitude from [0, 359] to [-179, 180]
   lon   = lon0
   thewest = where( lon0 gt 180 and lon0 lt 360, nwest)
   lon(thewest) = lon(thewest)-360
   ord = sort(lon)
   lon = lon(ord)
   sp    = sp(ord,*)
   tp    = tp(ord,*)
   st    = st(ord,*)
   temp  = temp(ord,*,*)

  ;(3) grid to center [0.5]
   lat   = (lat (0:179)+ lat(1:180))*0.5
   sp    = exp((alog(sp (*,0:179))+ alog(sp(*,1:180)))*0.5)
   tp    = exp((alog(tp (*,0:179))+ alog(tp(*,1:180)))*0.5)
   st    = (st (*,0:179)+ st(*,1:180))*0.5
   temp  = (temp (*,0:179,*)+ temp(*,1:180,*))*0.5

   nlon = n_elements(lon)
   lon1  = [nlon-1, indgen(nlon-1)]
   lon2  = [indgen(nlon)]
   lon   =  (lon(lon1)+lon(lon2))/2.0
   lon(0)= -179.5
   sp    =exp((alog(sp(lon1, *))+ alog(sp(lon2,*)))*0.5)
   tp    =exp((alog(tp (lon1,*))+ alog(tp(lon2,*)))*0.5)
   st    =(st (lon1, *)+ st(lon2,*))*0.5
   temp  =(temp (lon1, *,*)+ temp(lon2,*,*))*0.5
   fnl_new(k).sp = sp
   fnl_new(k).tp = tp
   fnl_new(k).st = st
   fnl_new(k).temp = temp
   fnl_new(k).lon  = lon
   fnl_new(k).lat  = lat
   fnl_new(k).pres = pres0
  ENDFOR
  print , 'rearray FNL'

END


;===================================================================================================
; PROGRAM NAME:
;  read_fnl_nc
;
; PURPOSE:
;  Read NCEP FNL Operational Model Global Tropospheric Analyses,
;  continuing from July 1999mapping column ozone with wind vector data.
;
; PROCESS:
;
; REFERENCE:
;
; INPUT:
;
; OUTPUT:
;
; DEVELOPER:
;  Daegeun Shin (geun)
;  Satellite Remote Sensing Laboratory
;  Division of Earth Environmental System
;  College of Natural Science
;  Pusan National University
;                  Tel: +82-51-510-2172
;                  E-mail: daegeun@pusan.ac.kr
;
; REVISION HISTORY:
;   original sorce from Daegeun Shin
;   Updated by Dae Sung Choi 2020-09-08
;   Updated by Dae Sung Choi 2022-04-12 stop -> message
;     Satellite Remote Sensing Laboratory
;     Division of Earth Environmental System
;     College of Natural Science
;     Pusan National University
;                  Tel: +82-51-510-2172
;                  E-mail: daesungchoi@pusan.ac.kr
;
; Copyright (C) 2020 Satellite Remote Sensing Laboratory, PNU
; All right reserved.
;
;===================================================================================================
PRO read_fnl_nc, date, utc, fnl, fnl_new, dir, get_rearray=get_rearray

  year = fix(strmid(date, 0, 4)) 
  month = fix(strmid(date, 4, 2))
  day = fix(strmid(date, 6, 2))

  ;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; obtain files for oneday
  ;;;;;;;;;;;;;;;;;;;;;;;;;;
  sel   = where ( fix(utc) ge 0 and fix(utc) le 18)
  da    = fix(utc(sel))/6
  file  = findfile(dir + 'fnl_'+date + '_'+'*'+'_*.nc', count=nfile) 
  if nfile ne 4 then begin
    message, 'FNL files for each 6 hours are not exists.'
  endif

  file  = file(da)
  sel1  = where( fix(utc) eq -6, nsel1)
  sel2  = where( fix(utc) eq 24, nsel2)

  if nsel1 eq 1 then begin
        day = day-1
     if day eq 0 then begin
        day = enday(year,month-1) & month = month - 1
     endif
     if month eq 0 then begin
        month = 12 & year = year - 1
     endif
     adate = string(year, format='(i4.4)') + string(month, format='(i2.2)')+ string(day, format='(i2.2)')
     afile = findfile(dir +'fnl_'+adate + '_18_*.nc', count=nfile) & if nfile ne 4 then stop
     file  = [afile(0), file]
  endif

  year = fix(strmid(date, 0, 4)) &  month = fix(strmid(date, 4, 2)) & day = fix(strmid(date, 6, 2))

  if nsel2 eq 1 then begin
        day = day+1
     if day gt enday(year,month) then begin
        day = 1 & month = month + 1
     endif
     if month gt 12 then begin
        month = 1 & year = year + 1
        dir   = '/data/MODEL/FNL/'+string(year,format='(i4.4)') + '/'
     endif
     adate = string(year, format='(i4.4)') + string(month, format='(i2.2)')+ string(day, format='(i2.2)')

     afile = findfile(dir +'fnl_'+adate + '_00_*.nc', count=nfile)
     if nfile ne 1 then begin
        print , 'No file for next'
        stop
     endif
     file = [file, afile(0)]
  endif
  nfile  = n_elements(file)

  varname =  ['PRES_3_SFC',  'PRES_3_TRO','TMP_3_SFC', 'TMP_3_ISBL']
  nlon = 360
  nlat = 181
  nl   = n_elements( ncread(file(0), 'lv_ISBL0'  ))
  var1 = fltarr(nlon, nlat)
  var2 = fltarr(nlon, nlat, nl)

  str   = {sp:var1, tp:var1,st:var1, temp:var2, lon:fltarr(nlon), lat:fltarr(nlat), pres:fltarr(nl), $
           uwind:var2, vwind:var2, omega:var2, pv:var2}
  fnl   = replicate(str,nfile)
  ;ncdf_cat, file(0)

  FOR k = 0 , nfile-1 do begin
      fname = file(k)
      surfp = ncread(fname, 'PRES_P0_L1_GLL0')
      tropp = ncread(fname, 'PRES_P0_L7_GLL0')
      surft = ncread(fname, 'TMP_P0_L1_GLL0')
      temp  = ncread(fname, 'TMP_P0_L100_GLL0')

      uwind  = ncread(fname, 'UGRD_P0_L100_GLL0')
      vwind  = ncread(fname, 'VGRD_P0_L100_GLL0')
      omega  = ncread(fname, 'VVEL_P0_L100_GLL0') 
      pv  = ncread(fname, 'lv_PVL3')
      lat = ncread(fname, 'lat_0')
      lon = ncread(fname, 'lon_0')
      pres  = ncread(fname,  'lv_ISBL0')

      fnl(k).sp   =  surfp  ;pa
      fnl(k).tp   =  tropp  ;pa
      fnl(k).st   =  surft ; (K)
      fnl(k).temp =  temp ; (K)
      fnl(k).lat = lat
      fnl(k).lon = lon
      fnl(k).pres= pres
      fnl(k).uwind= uwind
      fnl(k).vwind= vwind
      fnl(k).omega= omega
  ENDFOR

  IF keyword_set (get_rearray) then begin
     sub1_rearray_fnl , fnl, fnl_new
  ENDIF

END



;+---------------------------------------------------------------------------+
; Main Procedure : prepare_fnl_geun
;   Input  : date, lsttime
;   Output : fnlout 
; 
;+---------------------------------------------------------------------------+

PRO ds_convert_fnl_to_dat_for_gemso3p, date, dir, fnlout, outdat=outdat, lsttime=lsttime

  ;date    = '20120612'
  IF ~KEYWORD_SET(lsttime) THEN BEGIN
    lsttime = 13.75
  ENDIF

  fnltime = [0, 6, 12, 18]

  ;1) loading fnl data per day
  print, dir
  read_fnl_nc, date, fnltime, fnl1, fnlnew, dir

  fnl     = fnl1
  surfp0  = fnl.sp ; surface pressure
  tropp0  = fnl.tp ; tropopause pressure
  surft0  = fnl.st ; surface temperature
  temp0   = fnl.temp ; temperature profile
  uwnd0   = fnl.uwind ; uwind component
  vwnd0   = fnl.vwind ; vwind component
  omega0  = fnl.omega ; vertical veolocity
  lat0    = fnl(0).lat
  lon0    = fnl(0).lon
  pres0   = fnl(0).pres
  nlon0   = n_elements(lon0)
  nlat0   = n_elements(lat0)
  nl0     = n_elements(pres0)

  ;2) convert lsttime to utctime using longitude

  da   = where(lon0 gt 180.0, nda)

  if nda ne 0 then lon0(da) = lon0(da) - 360.

  ; OMI가 극궤도 주기로 지구를 도는 경우 적도를 지나는 시각 1345
  ; 이 시각에 맞추어 경도에 따라 다른 비율로 interpolation
  ;utctime = -lon0/15. + lsttime
  
  ; gems 는 한씬에 전 영역을 모두 관측하므로
  ; -lon0 대신 -125를 넣어서 일괄 계산
  tmp = fltarr(360)
  tmp[*] = 125
  utctime = - tmp /15. + lsttime

  da = where(utctime ge 24.0, nda)
  if nda gt 0 then utctime(da) = utctime(da) -24.
   da = where(utctime lt 0, nda)
  if nda gt 0 then begin
   utctime(da) = utctime(da) + 24.
   print, ' negative time'
  endif


  ;3) interpolate to lsttime

  surfp = fltarr(nlon0, nlat0)
  tropp = fltarr(nlon0, nlat0)
  surft = fltarr(nlon0, nlat0)
  tempp  = fltarr(nlon0, nlat0, nl0)
  uwnd  = fltarr(nlon0, nlat0, nl0)
  vwnd  = fltarr(nlon0, nlat0, nl0)
  omegap  = fltarr(nlon0, nlat0, nl0)

  utctime = utctime / 6.
  fidx    = floor(utctime)
  lidx    = fidx + 1
  frac2   = utctime-fidx
  frac1   = 1.0 - frac2

  FOR ilon = 0 , nlon0 - 1 do begin
    surfp(ilon,*)= surfp0[ilon, *,fidx[ilon]] *frac1[ilon] + surfp0[ilon, *,lidx[ilon]] *frac2[ilon]
    tropp(ilon,*)= tropp0[ilon, *,fidx[ilon]] *frac1[ilon] + Tropp0[ilon, *,lidx[ilon]] *frac2[ilon]
    surft(ilon,*)= surft0[ilon, *,fidx[ilon]] *frac1[ilon] + surft0[ilon, *,lidx[ilon]] *frac2[ilon]
    tempp(ilon,*,*)= temp0[ilon, *,*,fidx[ilon]] *frac1[ilon] + temp0[ilon, *,*,lidx[ilon]] *frac2[ilon]
    uwnd(ilon,*,*)= uwnd0[ilon, *,*,fidx[ilon]] *frac1[ilon] + uwnd0[ilon, *,*,lidx[ilon]] *frac2[ilon]
    vwnd(ilon,*,*)= vwnd0[ilon, *,*,fidx[ilon]] *frac1[ilon] + vwnd0[ilon, *,*,lidx[ilon]] *frac2[ilon]
    omegap(ilon,*,*)= omega0[ilon, *,*,fidx[ilon]] *frac1[ilon] + omega0[ilon, *,*,lidx[ilon]] *frac2[ilon]
  ENDFOR


  ;4) rearray

   pres  = pres0
  ;(1) rearray latitude from [90, -90] to [-90, 90]
   ord   = sort(lat0)
   lat   = lat0(ord)
   sp    = surfp(*,ord)
   tp    = tropp(*,ord)
   st    = surft(*,ord)
   temp  = tempp(*,ord,*)
   uwind  = uwnd(*,ord,*)
   vwind  = vwnd(*,ord,*)
   omega  = omegap(*,ord,*)


  ;(2) rearray latitude from [0, 359] to [-179, 180]
   lon   = lon0
   thewest = where( lon0 gt 180 and lon0 lt 360, nwest)
   if nwest ne 0 then lon(thewest) = lon(thewest)-360
   ord   = sort(lon)
   lon   = lon(ord)
   sp    = sp(ord,*)
   tp    = tp(ord,*)
   st    = st(ord,*)
   temp  = temp(ord,*,*)
   uwind  = uwind(ord,*,*)
   vwind  = vwind(ord,*,*)
   omega  = omega(ord,*,*)

   if nlat0 eq 181 then begin
  ;(3) grid to center [0.5]
     lat   = (lat (0:179)+ lat(1:180))*0.5
     sp    = exp((alog(sp (*,0:179))+ alog(sp(*,1:180)))*0.5)
     tp    = exp((alog(tp (*,0:179))+ alog(tp(*,1:180)))*0.5)
     st    = (st (*,0:179)+ st(*,1:180))*0.5
     temp  = (temp (*,0:179,*)+ temp(*,1:180,*))*0.5
     uwind  = (uwind (*,0:179,*)+ uwind(*,1:180,*))*0.5
     vwind  = (vwind (*,0:179,*)+ vwind(*,1:180,*))*0.5
     omega  = (omega (*,0:179,*)+ omega(*,1:180,*))*0.5

     nlon = n_elements(lon)
     lon1  = [nlon-1, indgen(nlon-1)]
     lon2  = [indgen(nlon)]
     lon   = (lon(lon1)+lon(lon2))/2.0
     lon(0)= -179.5
     sp    =exp((alog(sp(lon1, *))+ alog(sp(lon2,*)))*0.5)
     tp    =exp((alog(tp (lon1,*))+ alog(tp(lon2,*)))*0.5)
     st    =(st (lon1, *)+ st(lon2,*))*0.5
     temp  =(temp (lon1, *,*)+ temp(lon2,*,*))*0.5
     uwind  =(uwind (lon1, *,*)+ uwind(lon2,*,*))*0.5
     vwind  =(vwind (lon1, *,*)+ vwind(lon2,*,*))*0.5
     omega  =(omega (lon1, *,*)+ omega(lon2,*,*))*0.5
   endif

   fnlout = {sp:sp/100. , tp:tp/100., st:st, temp:temp, lon:lon, lat:lat, pres:pres/100., uwind:uwind, vwind:vwind, omega:omega}

   tp_size = size(tp, /dimension)

   openw, lun, outdat, /get_lun
   
   for iy=0, tp_size[1] - 1 do BEGIN
     printf, lun, tp[*, iy], format='(360i3)'
   ENDFOR
   free_lun, lun
END



