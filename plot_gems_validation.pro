
;============================
;
; 2018/05/23
; GEMS validation plot code
;
;============================

pro plot_gems_validation, xval, yval, filename=filename, cblim=cblim,$
  title=title, $
  xtitle = xtitle, $
  ytitle = ytitle, $
  ct = ctnum, $
  background = bgcolor, $
  range=range, delta=delta


if not keyword_Set(cblim) then begin
  cblim = [0., 100.]
endif
; set plotting parameters

result_path = './plot/'

if not keyword_Set(range) then begin
  minv = min([xval, yval]) 							;; set minimum value of x, y axis
  maxv = max([xval, yval]) 							;; set maximum value of x, y axis
  delta = 2.							;; set histogram interval of x, y axis
  range=[(floor(minv/delta)-1)*delta, (ceil(maxv/delta)+1)*delta]
endif

  
if not keyword_set(ctnum) then begin
  ctnum = 33
endif
cgPS_Open, 'test.ps', /landscape

cgDisplay, aspect=1
dsloadct, ctnum

position = [0.15, 0.15, 0.85, 0.9]
thick = (!D.Name EQ 'PS') ? 6 :3

x_para = xval 
y_para = yval 
if not keyword_Set(xtitle) then begin
  xtitle = 'GEMS Tropospheric O3'
endif
if not keyword_Set(ytitle) then begin
  ytitle = 'OMI Tropospheric O3'
endif
if not keyword_Set(bgcolor) then begin
  bgcolor = 'white'
endif

cgPlot, x_para, y_para, /NoData, $
  ;background=bgcolor, $
  XTitle =  xtitle, YTitle = ytitle, $
  XRange = range, YRange = range, $
  Position=position, Charsize=2.0, $
  XThick=3.0, YThick=3.0, Thick=6.0, /NoErase

xbox = [range[0], range[1], range[1], range[0]]
ybox = [range[0], range[0], range[1], range[1]]
polyfill, xbox, ybox, color=cgcolor(bgcolor)

;linfit_xdata=findgen(101)/100.*(range[1]-range[0])-minv
;back_eq1 = linfit_xdata*1.
;oplot, linfit_xdata, back_eq1, color=0, linestyle=0

if not keyword_set(title) then begin
  title = title
endif

cgPlot, x_para, y_para, /NoData, $
  ;background=bgcolor, $
  title = title, $
  XTitle =  xtitle, YTitle = ytitle, $
  XRange = range, YRange = range, $
  Position=position, Charsize=2.0, $
  XThick=3.0, YThick=3.0, Thick=6.0, /NoErase

;============================
; plot histogram
;============================
for j = range[0], range[1]-delta, delta  do begin
  for i = range[0], range[1]-delta, delta  do begin
  ;; delta 구간내에 있는 data 갯수 찾기
    idx_histo= where(y_para ge j $
      and y_para lt j+delta $
      and x_para ge i $
      and x_para lt i+delta, idx_histo_num)

    xbox = [i, i+delta, i+delta, i]
    ybox = [j, j, j+delta, j+delta]

    if idx_histo_num gt 0.0 then begin
      ;print, i, j, idx_histo_num
      polyfill, xbox, ybox, $
        color = bytscl(idx_histo_num, min=cblim[0] ,max=cblim[1], top=253)
    endif
  endfor
endfor

x_para = xval 
y_para = yval 
idx = where(finite(xval) eq 1 and finite(yval) eq 1, /null)
corr = correlate(x_para[idx], y_para[idx])

x_para = xval
y_para = yval 
coeff = linfit(x_para[idx], y_para[idx])

centerlinex = findgen(101)*(range[1]-range[0])/100.+range[0]
centerliney = findgen(101)*(range[1]-range[0])/100.+range[0]

x_para = xval 
y_para = yval 
linfitx = findgen(101)*(max([x_para, y_para])-min([x_para, y_para]))/100 + min([x_para, y_para])

x_para = xval 
y_para = yval 
linfity = coeff[0] +coeff[1]*linfitx

x_para = xval 
y_para = yval 
rmse = sqrt(mean((y_para-x_para)^2, /nan))

x_para = xval 
y_para = yval 
mae = mean(abs(x_para-y_para), /nan)

x_para = xval 
y_para = yval 
mbe = mean(x_para-y_para, /nan)

dsloadct, 0
oplot, centerlinex, centerliney, color=200, linestyle=0

dsloadct, 34
oplot, linfitx, linfity, color=253, linestyle=0
dsloadct, 0

; upper left box
;polyfill, [(range[1]-range[0])*0.15+range[0], $
            ;(range[1]-range[0])*0.45+range[0], $
            ;(range[1]-range[0])*0.45+range[0], $
            ;(range[1]-range[0])*0.15+range[0]], $
            ;[(range[1]-range[0])*0.45+range[0], $
            ;(range[1]-range[0])*0.45+range[0], $
            ;(range[1]-range[0])*0.95+range[0], $
            ;(range[1]-range[0])*0.95+range[0]], color=200

;polyfill, [(range[1]-range[0])*0.45+range[0], $
            ;(range[1]-range[0])*0.95+range[0], $
            ;(range[1]-range[0])*0.95+range[0], $
            ;(range[1]-range[0])*0.45+range[0]], $
            ;[(range[1]-range[0])*0.05+range[0], $
            ;(range[1]-range[0])*0.05+range[0], $
            ;(range[1]-range[0])*0.45+range[0], $
            ;(range[1]-range[0])*0.45+range[0]], color=200

xyouts, (range[1]-range[0])*0.05+range[0], (range[1]-range[0])*0.95+range[0], $
  'N(#) ='+string(n_elements(xval)), /data

xyouts, (range[1]-range[0])*0.05+range[0], (range[1]-range[0])*0.90+range[0], $
  'CORR ='+string(corr), /data

If coeff[0] lt 0 then begin
  interceptstr = ' - ' + strtrim(string(abs(coeff[0]), format='(f9.3)'), 2)
endif else begin
  interceptstr = ' + ' + strtrim(string(abs(coeff[0]), format='(f9.3)'), 2)
ENDELSE

xyouts, (range[1]-range[0])*0.05+range[0], (range[1]-range[0])*0.85+range[0], $
  'Y =   '+strtrim(string(coeff[1], format='(f7.3)'), 2) + $
  'X ' + interceptstr


;lower right box
xyouts, (range[1]-range[0])*0.55+range[0], (range[1]-range[0])*0.12+range[0], $
  'RMSE = '+string(rmse), /data

xyouts, (range[1]-range[0])*0.55+range[0], (range[1]-range[0])*0.07+range[0], $
  'MAE =  '+string(mae), /data

xyouts, (range[1]-range[0])*0.55+range[0], (range[1]-range[0])*0.02+range[0], $
  'MBE =  '+string(mbe), /data

;xyouts, (range[1]-range[0])*0.45+range[0], (range[1]-range[0])*0.50+range[0], $
  ;'RMSE='+string(corr), /data

;xyouts, (range[1]-range[0])*0.45+range[0], (range[1]-range[0])*0.45+range[0], $
  ;'Y='+strtrim(string(coeff[1]),2) + 'x+'+strtrim(string(coeff[0]), 2), /data

dsloadct, ctnum

cgColorbar,format='(i4)', $
  Color='black',$
  range=cblim, $
  position=[0.87, 0.15, 0.9, 0.9],$
  /vertical,$
  /right, $ 
  divisions=5, $
  bottom=0, $
  ncolor=253, $
  minor=1, $
  Charsize=2
cgPS_Close

pngfile = filename
cgPS2Raster, 'test.ps', pngfile, density =1000, /png
spawn, 'convert ' + pngfile + $
  ' -background "rgba(0,0,0,0.5)" -rotate -90 ' + pngfile


scp_dest = 'soodal@164.125.38.179:/home/soodal/works/plot'

; send image to pc
spawn, 'scp -P18742 -p ' + pngfile + $
  ' ' + scp_dest


end
