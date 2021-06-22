PRO load_SONDE_v10 , sondes, list, uncorr=uncorr, sonfile=sonfile, sta=sta, quiet=quiet

IF keyword_set ( sonfile ) THEN BEGIN
    restore, sonfile
ENDIF ELSE BEGIN
    restore,  homedir()+'Data/SONDE/'+'all_sondes_2004-2014_v10.xdr'
ENDELSE

sonde_toznorm = 1
IF keyword_set(uncorr) then sonde_toznorm = 0
mncorr  = 0.85   &  mxcorr = 1.15

IF keyword_set (sta) THEN BEGIN
   da = where( strpos(sondes.station, sta) ne -1, nda)
   sondes = sondes(da)
ENDIF

if sonde_toznorm eq 1 then begin
    da = where(sondes.corr eq 0, ct)
    if ct gt 0 then sondes(da).corr = 1.0 ; no correction factor is provided
    nsonde = n_elements(sondes)
    i = 0l
    while i lt nsonde do begin
        nlvl = sondes(i).nlvl
        corr = abs(sondes(i).corr)
        if corr ge mncorr and corr le mxcorr then begin
           IF sondes(i).corr lt 0 then sondes(i).o3(0:nlvl-1) = sondes(i).o3(0:nlvl-1) * corr
        endif else begin
        ;   IF sondes(i).corr gt 0 then sondes(i).o3(0:nlvl-1) = sondes(i).o3(0:nlvl -1) /corr
        ;   sondes(i).ok = 2
        endelse
        i = i + 1
    endwhile
endif else if sonde_toznorm eq 0 then begin
    da = where( sondes.ok eq 1 and sondes.corr gt 0, nda)
    i = 0L
    while i lt nda do begin
       corr = sondes(da(i)).corr
       nlvl = sondes(da(i)).nlvl
      ; sondes(da(i)).o3(0:nlvl-1) = sondes(da(i)).o3(0:nlvl-1)/corr
      ; sondes(da(i)).ok = 0
       i = i + 1
    endwhile
endif

;---------------------------------------------------------------------------------
IF not keyword_set (quiet) then begin
list = get_sonlist( sondes)
ENDIF ELSE BEGIN
list = get_sonlist( sondes, /quiet)
ENDELSE
return
END
