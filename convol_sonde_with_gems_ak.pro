PRO calc_o3cum ,use_spline, pres0, cum0, Ptop, Pbot, out

  pres = pres0
  cum  = cum0

  ; bottom-up to top-down
  if pres(0) gt pres(1) then begin
    pres = reverse(pres) 
    cum = reverse(cum) 
  endif

  out = !values.f_nan

  ; input error
  if ptop ge pbot then begin
    print, 'ptop ge pbot'
    return
  endif

  if ptop lt min(pres) or pbot gt max(pres) then begin
    print, 'ptop lt min(pres) or pbot gt max(pres)'

    return
  endif

  ; spline 
  if use_spline eq 1 then begin
    temp = spline(alog(pres), cum, alog([ptop, pbot]))
  endif else begin
    temp = interpol(cum, alog(pres), alog([ptop, pbot]))
  endelse
  out = temp(1) - temp(0)
END

;==============================================================================
; Note ozonesonde input profiles: bottom up
;                 output profiles:  top down 
; input
;   sonpres     : sonde pressure profile (hPa)
;   sono3       : sonde partial ozone column (DU)
;   omipres     : OMI pressure profile
;   apozprof    : OMI a priori profile
;   ozprof      : OMI retrieved profile
;   avgk        : OMI avgk
;   ntp         : OMI tropopause pressure level if ntp eq -1 then define 
;                 tropopause with insert_ptrp
;
; output 
;   sonstlvl    : the level in omipres corresponding to the top level of sonpres
;   psontop     : psontop = omipres(sonstlvl) 

; option
;   use_spline  : (present) apply spline interpolation (Not. present) apply liear interpolatino
;   convol      : convolved sonprof with omiavgk is calculated
;   cconvol     : convolved column amount is calculated
;   append_ret  : partial column above sonde top is filled with retrieved ozone 
;                 with a priori ozone if not present. 
;   insert_ptrp : if present as 500 hPa, omisco is for partial column between 
;                 500 hPa to sontop omitco is for partial column between 
;                 bottom to 500 hPa
;   midpres     : if present as 100 hPa, omisco200 is for partial column 
;                 between 100hpa to sontop omitco200 is for partial column 
;                 between bottom to 100hpa

pro convol_sonde_with_gems_ak, $
  sonpres0, sono30, omipres, apozprof, ozprof, avgk, ntp, $ ; inputs
  sonstlvl, sondu, sonsco, sontco, sonsco200, sontco200, omisco, $
  omitco, omisco200, omitco200, omiasco, omiatco, omiasco200, omiatco200, $
  psontop, convol=convol, cconvol = cconvol, use_spline=use_spline, $
  append_ret = append_ret, print_result = print_result, midpres=midpres, $
  insert_ptrp=insert_ptrp

; (1) define variable
missing = !values.f_nan
if not(keyword_set(convol))       then convol   = 0
if not(keyword_set(cconvol))      then cconvol   = 0
if not(keyword_set(use_spline))   then use_spline   = 0
if not(keyword_set(print_result)) then print_result = 0
if not(keyword_set(append_ret))   then append_ret = 0
if not(keyword_set(midpres))      then midpres = 200.

nl     = n_elements(ozprof)
sonsco    = missing ; from omipres(sonstlvl) to omipres(ntp)
sontco    = missing ; from son(nsl) to omipres(ntp)
sonsco200 = missing ; from omipres(sonstlvl) to 200
sontco200 = missing ; from son(nsl) to 200
omitco    = missing
omitco200    = missing
omiatco200 = missing
omiatco   = missing
omisco    = missing ; from omipres(nl) to 200
omisco200 = missing ; from omipres(sonstlvl) to omipres(ntp) 
omiasco    = missing
omiasco200 = missing
sondu0 = fltarr(nl)  
sondu0(*) = missing  ; original sonde on OMI grid
sondu  = sondu0  ; initialized with nan, this variable will filled with convolved sonde data

; (2) cumulated column ozone : cumson, cuma, cum
sonpres = sonpres0  
sono3   = sono30
nsl   = n_elements(sono3)
ptop  = min(sonpres)

; Change it to top-down if it is bottom-up 
if sonpres(0) gt sonpres(1) then begin
 sonpres = reverse(sonpres)  
 sono3   = reverse(sono3)
endif
cumson = [0, cum_total(sono3)]


cuma = [0, cum_total(apozprof)]
cum  = [0, cum_total(ozprof)]

;-----------------------------------------------------------------------------
; Adjust ozonesonde data if OMI surface pressure is larger
;-----------------------------------------------------------------------------
if omipres(nl) gt sonpres(nsl)   then begin
    tmp = (cumson(nsl)-cumson(nsl-1)) * (omipres(nl) - sonpres(nsl-1)) $
                                      / (sonpres(nsl) - sonpres(nsl-1))
    cumson(nsl)  = cumson(nsl-1) + tmp
    sonpres(nsl) = omipres(nl)
endif

da = where(omipres ge sonpres(0))
sonstlvl = da(0)
if sonstlvl eq -1 then begin
  print, 'sonstlvl eq -1'
  stop
  return
endif
psontop = omipres(sonstlvl)   ; top pressure used in comparison
if ntp ge 0 then begin
  ptrop   = omipres(ntp)
endif
if keyword_set (insert_ptrp) then  ptrop   = insert_ptrp

ntemp = nl - sonstlvl 
if ntemp le 3 then begin
  print, 'ntemp le 3'
  return
endif


; OMI retrieval
    calc_o3cum ,use_spline, omipres, cum, psontop, midpres,     omisco200
    calc_o3cum ,use_spline, omipres, cum, midpres, omipres(nl), omitco200
    calc_o3cum ,use_spline, omipres, cum, psontop, ptrop,       omisco
    calc_o3cum ,use_spline, omipres, cum, ptrop,   omipres(nl), omitco 
; OMI a priori
    calc_o3cum ,use_spline, omipres, cuma, psontop, midpres,     omiasco200
    calc_o3cum ,use_spline, omipres, cuma, midpres, omipres(nl), omiatco200
    calc_o3cum ,use_spline, omipres, cuma, psontop, ptrop, omiasco
    calc_o3cum ,use_spline, omipres, cuma, ptrop, omipres(nl), omiatco 

;-----------------------------------------------------------------------------
; Interpolate ozonesonde profile to OMI altitude grid
;-----------------------------------------------------------------------------
sondu0 = fltarr(nl)  
sondu0(*) = missing ; original sonde on OMI grid
sondu  = sondu0  ; initialized with nan, this variable will filled with convolved sonde data
if use_spline eq 1 then begin
    tmpcum = spline(alog(sonpres), cumson, alog(omipres(sonstlvl:nl)))
    da = where ( finite(tmpcum) eq 0 , nda)
    if nda ne 0 then begin
      print , 'spline error' 
       stop
      tmpcum = interpol(cumson, alog(sonpres), alog(omipres(sonstlvl:nl)))
    endif
endif else begin
    tmpcum = interpol(cumson, alog(sonpres), alog(omipres(sonstlvl:nl)))
endelse
sondu0(sonstlvl:nl-1)= tmpcum(1:ntemp) - tmpcum(0:ntemp-1)

;-----------------------------------------------------------------------------
; Apply OMI averaging kernels
;-----------------------------------------------------------------------------
if convol eq 1 then begin
    if append_ret eq 1 then begin
        sondu0(0:sonstlvl-1) = ozprof(0:sonstlvl-1)
    endif else begin
        sondu0(0:sonstlvl-1) = apozprof(0:sonstlvl-1)  ; essentially do nothing
    endelse

    diff = sondu0 - apozprof

    for i = sonstlvl, nl - 1 do begin
        sondu(i) = apozprof(i) + total(avgk(*, i) * diff)
    endfor
endif else begin
    sondu = sondu0
endelse

;-----------------------------------------------------------------------------
; Get convolved ozone columns if necessary
;-----------------------------------------------------------------------------
    tmpres = omipres(sonstlvl:nl)
    ntemp  = nl - sonstlvl 
    tmpcum = [0, cum_total(sondu(sonstlvl:nl-1))]   
    calc_o3cum ,use_spline, tmpres,  tmpcum, psontop, midpres, sonsco200
    calc_o3cum ,use_spline, tmpres,  tmpcum, midpres, omipres(nl), sontco200
    calc_o3cum ,use_spline, tmpres,  tmpcum, psontop, ptrop, sonsco
    calc_o3cum ,use_spline, tmpres,  tmpcum, ptrop,   omipres(nl), sontco    
; OMI retrieval
    calc_o3cum ,use_spline, omipres, cum, psontop, midpres,     omisco200
if print_result eq 1 then begin
    print, ''
    print, 'Top pressure level below burst & tropopause pressure: ', psontop, ptrop
    print, 'Surface Pressure: ', sonpres0(0), sonpres(nsl), omipres(nl)
    print, 'Ozonesonde profile before convolution: '
    print, sondu0, format='(12F8.2)'
    print, 'Ozonesonde profile after convolution:'
    print, sondu, format='(12F8.2)'

    print, 'Sonde SCO, TCO, SCO200, TCO200: ', sonsco, sontco, sonsco200, sontco200, format='(A, 4F8.2)'
    print, 'OMI   SCO, TCO, SCO200, TCO200: ', omisco, omitco, omisco200, omitco200, format='(A, 4F8.2)'
    print, 'OMIap SCO, TCO, SCO200, TCO200: ', omiasco, omiatco, omiasco200, omiatco200, format='(A, 4F8.2)'

endif


return
end
