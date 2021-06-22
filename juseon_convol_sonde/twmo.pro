;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; WMO definition in calculating thermal tropopause
; calculate pressure, temperature, and height of tropopause using thermal criterium
; by Thomas Reichler
; usage:
; twmo, -0.002, input temperature, surfacealt, surfacepres, inputpress, outputptrp, outputtrp, outputztrop
; TWMO, -0.002, temp(from top to bottom),surface altitue (in meter), surface pressure (in Pa), vertical pressure (in Pa), $
;      ptrp (output press in Pa), ttrp (output temperaure), ztrp (output altitude in meter),0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro  TWMO, ZGWMO, T0_INI, ZS, PS, PFULL0, $   ; input variables
           PTRP, TTRP, ZTRP, PL         ; output variables

; convert the unit of temperature to kelvin
T0 = T0_INI
IF max(T0_INI) < 100 THEN  BEGIN
    T0 = T0_INI +273.15
ENDIF

T = T0
PFULL = PFULL0
IF PFULL(0) gt PFULL(1) THEN BEGIN ; sort by top-down
   PFULL = reverse(PFULL)
   T     = reverse(T)  
ENDIF

level 	= n_elements(PFULL); input pressures


plimu 	= 45000.    
plimo 	= 5000.    

kap    	= 0.286
ka1    	= kap - 1.0
faktor 	= -9.81/287.0
deltaz 	= 2000.0
dtdz0  	= 0.					; data from level below
pmk0   	= 0.					
pm0   	= 0.
tm0   	= 0.
a0	= 0.
b0	= 0.

PTRP  	= .0
TTRP  	= .0
ZTRP  	= .0

;;;;;;;;;;;;;;;;;;;;;;;;;
J=LEVEL-1
REPEAT BEGIN
;;;;;;;;;;;;;;;;;;;;;;;;;

pmk = .5 * (PFULL(j-1)^kap + PFULL(j)^kap)      ; p mitte ^kappa
pm = pmk^(1./kap)                               ; p mitte



a = (t(j-1)-t(j))/(PFULL(j-1)^kap-PFULL(j)^kap)
b = t(j)-(a*PFULL(j)^kap)
tm = a * pmk + b                                ; T mitte
dtdp = a * kap * (pm^ka1)
dtdz = faktor*dtdp*pm/tm



if j ne level-1 and pm le plimu and dtdz gt zgwmo and dtdz0 le zgwmo then begin

  ag = (dtdz-dtdz0) / (pmk-pmk0)      
  bg = dtdz0 - (ag * pmk0)           
  ptph = exp(alog((zgwmo-bg)/ag)/kap)

  ; calculate temperature at this ptph assuming linear gamma interpolation
  prs=pm0
  inttmp0=(bg*alog(prs)+ag*(prs^kap)/kap)/faktor*t(j)
  prs=ptph
  inttmp=(bg*alog(prs)+ag*(prs^kap)/kap)/faktor*t(j)
  ttph=tm0+inttmp-inttmp0
  ;if ptph lt PFULL(j) then begin
  ;   ttph = a * ptph^kap + b                 
  ;endif else begin
  ;   ttph = a0 * ptph^kap + b0                  
  ;endelse

  if PL then begin

    ; plot temperature profile assuming linear gamma interpolation
    dpn=10
    dp=(pm0-pm)/dpn
    for prs=pm0,pm,-dp do begin
      inttmp=(bg*alog(prs)+ag*(prs^kap)/kap)/faktor*t(j)
      tmp=tm0+inttmp-inttmp0
;      print,tmp,prs/100.
      oplot, [tmp], [prs]/100., psym=1
    endfor
   
    oplot, [ttph], [ptph]/100., psym=2, color=240 ; tropopause
  endif

  if ptph ge plimo and ptph le plimu then begin

    ; dt/dz above 2km must not be lower then -2K/km
    p2km = ptph + deltaz*(pm/tm)*faktor     	; p at ptph + 2km
    asum = .0                               	; dtdz above
    icount = 0                              	; number of levels above
    PTRP=0.
    
    JJ=J
    repeat begin
      pmk2 = .5 * (PFULL(jj-1)^kap+PFULL(jj)^kap)   	; p mitte ^kappa
      pm2 = pmk2^(1./kap)                   	; p mitte
      if pm2 lt p2km then begin
        PTRP=ptph                           	; ptph valid
      endif else begin
          a2 = (t(jj-1)-t(jj))              	; a
          a2 = a2/(PFULL(jj-1)^kap-PFULL(jj)^kap)
          b2 = t(jj)-(a2*PFULL(jj)^kap)         	; b
          tm2 = a2 * pmk2 + b2              	; T mitte
          dtdp2 = a2 * kap * (pm2^(kap-1))  	; dt/dp
          dtdz2 = faktor*dtdp2*pm2/tm2      	; gamma
          asum = asum+dtdz2
          icount = icount+1
          aquer = asum/icount               	; dt/dz mean
          if aquer le zgwmo  then begin     	; dt/dz above < 2K/1000
             PTRP=-1.                       	; to skip loop
          endif
      endelse
    JJ=JJ-1
    endrep until PTRP ne 0. or JJ eq 1

    if PTRP ne 0. and PTRP ne -1. then begin
      PTRP=ptph 
      TTRP=ttph 
    endif else begin
      PTRP=0. 
      TTRP=0. 
    endelse
  endif

endif
 
pm0 = pm
tm0 = tm
pmk0  = pmk
dtdz0 = dtdz
a0 = a
b0 = b

;;;;;;;;;;;;;;;;;;;;;;;;;
J=J-1
ENDREP UNTIL J eq 1 or PTRP ne 0.
;;;;;;;;;;;;;;;;;;;;;;;;;

; calculate geometric height of tropopause by integrating hypsometric equation
; from tropopause pressure down to surface pressure
; Ztp - Zs = R/g Int(T)dlnp]

;;;;;;;;;;;;;;;;;;;;;;;;;
if PTRP ne 0. and PTRP lt PS then begin		; only if valid tropopause found
;;;;;;;;;;;;;;;;;;;;;;;;;

R=287.
g=9.81

if 1 then begin	; Note: both methods lead to identical results

;;;;;;;;;;;;;;
; method I
;;;;;;;;;;;;;;

J = LEVEL-1					; bottom
DLNP = alog(PS/PFULL(j))			; from surface pressure
TM = T(j)					; take TM of lowest level (better: extrapolate)
TDLNP = TM*DLNP

while J ge 1 and PTRP lt PFULL(j-1) do begin 
  DLNP = alog(PFULL(j)/PFULL(j-1))
  TM = 0.5 * (T(j) + T(j-1))
  TDLNP = TDLNP + TM*DLNP
  J=J-1 
endwhile

DLNP = alog(PFULL(j)/PTRP)			; up to tropopause pressure
TM = 0.5 * (T(j) + TTRP)			; use TTRP to get TM of this level 
TDLNP = TDLNP + TM*DLNP

ZTRP = ZS + TDLNP*R/g
;print, ZTRP

endif else begin

;;;;;;;;;;;;;;
; method II, interpolate T mean in p^kappa
;;;;;;;;;;;;;;

J = LEVEL-1                                     ; bottom
DLNP = alog(PS/PFULL(j))                        ; from surface pressure
TM = T(j)                                       ; take TM of lowest level (better: extrapolate)
TDLNP = TM*DLNP

while J ge 1 and PTRP lt PFULL(j-1) do begin 
  DLNP = alog(PFULL(j)/PFULL(j-1))
  PMK = .5 * (PFULL(j-1)^kap + PFULL(j)^kap) 	; p mean ^ kappa
  a = (t(j-1)-t(j))/(PFULL(j-1)^kap-PFULL(j)^kap)
  b = t(j)-(a*PFULL(j)^kap)
  TM = a * PMK + b                      
  TDLNP = TDLNP + TM*DLNP
  J=J-1 
endwhile

DLNP = alog(PFULL(j)/PTRP)			; up to tropopause pressure
PMK = .5 * (PTRP^kap + PFULL(j)^kap) 		; p mean ^ kappa
a = (TTRP-t(j))/(PTRP^kap-PFULL(j)^kap)		; now use PTRP and TTRP
b = t(j)-(a*PFULL(j)^kap)
TM = a * PMK + b                      
TDLNP = TDLNP + TM*DLNP

ZTRP = ZS + TDLNP*R/g
print, ZTRP

endelse

;;;;;;;;;;;;;;;;;;;;;;;;;
endif
;;;;;;;;;;;;;;;;;;;;;;;;;

if PL eq 1 then begin
  XM=[3,3]
  YM=[0,0]

  PRANGE=[130,70]	; tropics
  TRANGE=[180,220]

  PRANGE=[280,200]	; extratropics
  TRANGE=[200,260]

  PRANGE=[500,100]	; antarctica
  TRANGE=[180,230]

  PRANGE=[500,100]	; antarctica
  TRANGE=[200,260]
  yrange = [1000,10]
  xrange =[min(t0), max(t0)]
  set_plot, 'ps'
  mywindow
  plot,findgen(10), findgen(10), $
        yrange=yRANGE, ystyle=1, $
        xrange=xRANGE, xstyle=1
  oplot, T0, PFULL0/100., linestyle=0
  plots, xrange,[ptrp, ptrp]/100, color=200

endif
return
end

