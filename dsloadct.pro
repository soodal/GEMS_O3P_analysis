PRO dsloadct, index, NCOLORS=ncolors, XYTHICK=xyThick, ctreverse=ctreverse, wb=wb

;!X.MINOR=1 & !Y.MINOR=1
;!X.STYLE=1 & !Y.STYLE=1

;IF NOT Keyword_Set(index) THEN BEGIN
  ;index=33
;ENDIF

IF NOT Keyword_Set(ctreverse) THEN BEGIN
  ctreverse=0
ENDIF

IF NOT Keyword_Set(wb) THEN BEGIN
  wb=0
ENDIF

IF Keyword_Set(xyThick) THEN BEGIN
  IF KeyWord_Set(xyThick) EQ 1 THEN BEGIN
    SET_PLOT, 'x'
    DEVICE, DECOMPOSED = 0
    !X.THICK=xyThick & !Y.THICK=xyThick
  ENDIF ELSE IF N_Elements(xyThick) EQ 2 THEN BEGIN
    SET_PLOT, 'x'
    DEVICE, DECOMPOSED = 0
    !X.THICK=xyThick[0] & !Y.THICK=xyThick[1]
  ENDIF
ENDIF else begin
  !X.THICK=1 & !Y.THICK=1
ENDELSE

;IF index EQ 20 THEN BEGIN
  ;TVLCT,r,g,b,/GET
  ;r_=CONGRID(r[100:255],255)
  ;g_=CONGRID(g[100:255],255)
  ;b_=CONGRID(b[100:255],255)
  ;TVLCT,r_,g_,b_
;ENDIF

;@@ define black and white color
If Keyword_Set(ncolors) Then BEGIN
  nc=ncolors
ENDIF Else BEGIN
  nc=256
ENDelse
LOADCT, index, NCOLORS=nc, /SI

IF ctreverse THEN BEGIN
  tvlct, r, g, b, /GET
  tvlct, reverse(r), reverse(g), reverse(b)
ENDIF

IF wb THEN BEGIN
  TVLCT, 0,0,0, 255 ; Drawing color.
  TVLCT, 255,255,255, 254 ; Background color.
  !p.color = 255
  !p.background = 254
ENDIF

RETURN
END

