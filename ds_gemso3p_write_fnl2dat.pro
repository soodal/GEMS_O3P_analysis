; 2013.4.30 JBAK
; subroutines required to run write_fnl_for_sao.pro
; 1) read_fnl_nc     ; load original data from FNLNCEP files
; 2) prepare_fnl_sao ; interpolation to local time 
;                      rearray lon/lat
;                      grid to center
; ~ 2011.09 by Xiong / after then by Juseon
; 2021.09 by Dae Sung Choi

pro ds_gemso3p_write_fnl2dat, date, lsttime, fnl, outpath=dir

nl  = n_elements(fnl.pres)
nlat= n_elements(fnl.lat)
nlon= n_elements(fnl.lon)
surfp=fnl.sp
tropp=fnl.tp
surft=fnl.st
temp =fnl.temp

;dir = '/home/geun/2_O3PR/ATMOS/fnl'+string(lsttime, format='(f5.2)')+'LST/'
if not keyword_set(dir) then begin
  dir = '/home/geun/2_O3PR/tool_for_geun/fnl/fnlout/'
  dir = '/OPER/SYSTEMS/ALGRTH/GEMS/L2/data/o3p/ATMOS/fnl13.75LST/'
endif

; write surface pressure
close,1 &  close,2 & close,3 & close,4 
openw, 1, dir + 'fnlsp/fnlsp_'+date+'.dat'
openw, 2, dir + 'fnltp/fnltp_'+date+'.dat'
openw, 3, dir + 'fnlst/fnlst_'+date+'.dat'
openw, 4, dir + 'fnltemp/fnltemp_'+date+'.dat'

FOR i = 0 , nlat-1 do begin
 printf,1, round(surfp(*,i)), format='(360i4)'
 printf,2, round(tropp(*,i)), format='(360i3)'
 printf,3, round(surft(*,i)), format='(360i3)'
ENDFOR

;IF nl NE 26 THEN BEGIN
  ;temp=REFORM(temp[*,*,5:nl-1])
  ;IF nl ne 31 THEN STOP
  ;nl=26
;ENDIF
fnl_layers = [10., 20., 30., 50., 70., 100., $
  150., 200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 700., $
  750., 800., 850., 900., 925., 950., 975., 1000.]
sz = size(temp, /dimension)

temp_new = fltarr(sz[0], sz[1], 26)
case nl of 
  31: begin
    temp_new=reform(temp[*, *, 5:nl-1])
    nl = 26
  end
  41: begin
    for il=0, n_elements(fnl_layers)-1 do begin
      layer_idx = where(fnl.pres eq fnl_layers[il], /null)
      if n_elements(layer_idx) eq 0 then begin
        STOP
      endif
      temp_layer = temp[*, *, layer_idx]
      temp_new[*, *, il] = temp_layer
    ENDFOR
    nl = 26
  end
else: stop
endcase


FOR j = 0 , nl -1 do begin
    FOR i = 0 , nlat-1 do begin
        printf,4, round(temp_new(*,i,j)), format='(360i3)'
    ENDFOR
ENDFOR

close,1 &  close,2 & close,3 & close,4 
return

END 
