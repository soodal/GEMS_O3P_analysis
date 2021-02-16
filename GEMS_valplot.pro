;pro GEMS_valplot

;============================
;
; 2018/05/23
; GEMS validation plot code
;
;============================

;============================
; read file
;============================


	result_path = '/home/myid0121/2016GEMS/'
	result_files = file_search(result_path+'outfile.txt', count=file_num)

	for f = 0L, file_num-1  do begin

		lines = file_lines(result_files(f))
		if lines lt 5 then goto, file_out

		data = fltarr(3,lines)
		openr,1,result_files(f)
		readf,1,data
		close,1

		xval = fltarr(lines) & xval(*) = reform(data[1,*])
		yval = fltarr(lines) & yval(*) = reform(data[2,*])

	file_out :
	endfor


;============================
; plot val data
;============================


minv = 0.0 							;; set minimum value of x, y axis
maxv = 3.0 							;; set maximum value of x, y axis
delta = 0.03 							;; set histogram interval of x, y axis



cgPS_Open, 'test.ps'

cgDisplay, aspect=1
LoadCT, 33, Ncolors=254, Bottom=1

position = [0.15, 0.15, 0.85, 0.9]
thick = (!D.Name EQ 'PS') ? 6 :3

x_para = xval           & X_name = 'AERONET_AOD [443nm]'
y_para = yval           & Y_name = 'GEMS_AOD [443nm]'


cgPlot, x_para, y_para, /NoData, $
			XTitle =  X_name, YTitle = Y_name, $
			XRange = [minv, maxv], YRange = [minv, maxv], $
			Position=position, Charsize=2.0,XThick=3.0, YThick=3.0, Thick=6.0, /NoErase

linfit_xdata=findgen(200)-100
back_eq1 = linfit_xdata*1.
oplot, linfit_xdata, back_eq1, color=0, linestyle=1


;============================
; plot histogram
;============================

			for j = minv, maxv, delta  do begin
				for i = minv, maxv, delta  do begin

				;; delta 구간내에 있는 data 갯수 찾기
				idx_histo= where( yval ge j    and    yval lt j+delta    and    xval ge i    and    xval lt i+delta,    idx_histo_num)

				xbox = [i, i+delta, i+delta, i]
				ybox = [j, j, j+delta, j+delta]

				if idx_histo_num gt 0.0 then begin

				print, i, j, idx_histo_num
				polyfill, xbox, ybox, color = bytscl(idx_histo_num, min = 0.0 , max = 5.0,top = 254)    ;; set color bar min / max

				endif

				endfor
			endfor

		cgColorbar,format='(i4)',range=[0.0,5.0], position=[0.87, 0.15, 0.9, 0.9],color=0, /vertical,/right, $      ;; set color bar min / max
			                divisions=5,bottom=1,ncolor=254,minor=1, Charsize=2

cgPS_Close



cgPS2Raster, 'test.ps', result_path+ 'GEMSAODval_20050101-20051231_histo.png',density =1000,/png




stop

end





















