pro plot_gems_satellite, data, lon, lat, $
  outfile=outfile, $
  scp_send=scp, $
  scp_dest=scp_dest1, $
  range=range, $
  title=pic_title, $
  cb_title=cb_title

if not keyword_Set(pic_title) then begin
  pic_title = 'test'
endif

if not keyword_Set(cb_title) then begin
  cb_title = pic_title
endif

if not keyword_Set(range) then begin
  range = minmax(data)
endif

;;----------------------------
;; call retrieval value
;;----------------------------


;;----------------------------
;; set dimension
;;----------------------------

sz = size(data)
dim1 = sz[1]
dim2 = sz[2]

;;----------------------------
;; plot procedure
;;----------------------------

;x_para = elon           & X_name = 'Longitude'
;y_para = elat           & Y_name = 'Latitude'

Slat = -3.
Nlat = 50.
Llon = 75.
Rlon = 160.

cgPS_Open, 'gems.ps'
cgDisplay, aspect=0.7
dsLoadCT, 22;bottom=1

pos = [0.1,0.15,0.9,0.97]
charsize = (!D.Name EQ 'PS') ? cgDefCharsize()*0.7 : cgDefCharsize()*0.75
thick = (!D.Name EQ 'PS') ? 6 :3
londel = 10
latdel = 10


cgMap_set,36,128, Limit=[-60,70 , 60, 160],$
		           /noerase, /advance,pos=pos,charsize=charsize,$
		          /satellite,/noborder,/iso,$
		          title=pic_title

for j = 1L, dim2-3 do begin
  for i = 1L, dim1-3 do begin

    if (lat(i,j) ge Slat) and $
        (lat(i,j) le Nlat) and $
        (lon(i,j) ge Llon) and $
        (lon(i,j) le Rlon) then begin

      if data(i,j) gt 0.01 then begin
              xbox = [lon(i,j), lon(i+1,j), lon(i+1,j), lon(i,j)]
              ybox = [lat(i,j), lat(i,j), lat(i,j+1), lat(i,j+1)]
      polyfill, xbox, ybox, $
        ;color = bytscl(data(i,j), $
        color=bytscl(data[i,j], $
          min=range[0], $
          max=range[1], $
          top=253)
        ;min = 0.000001, $
        ;max = 1.5, $
        ;top = 255)
      endif
    endif


  endfor
endfor

cgMap_Grid, Color='black', GLinestyle=1,londel=londel,latdel=latdel,charsize=1 $
			           , LONLAB=20, LATLAB=160, LABEL=1

cgMap_continents,color='black'

cgColorbar, format='(f6.1)', range=range, $
  Position=[pos[0]+0.065, pos[1]-0.04, pos[2]-0.065, pos[1]-0.02], $
  color='black', $
	OOB_High=253, $
  OOB_low=0, $
	divisions=5, $
  bottom=0, $
  ncolor=254, $
  minor=1, $
  charsize=1.5, $
	XTicklen=1, $
  XMinor=0, $
 	AnnotateColor='black', $
  Title=cb_title


;-----------------------------------------------------------------------
; * colorbar title for each product
;
; * NO2 --> title =  'GEMS NO!d2 !n [x 10!e15!n molecule/cm !e2 !n]',
; * TropNO2 --> title =  'GEMS Trop NO!d2 !n [x 10!e15!n molecule/cm !e2 !n]',
; * SO2 --> title =  'GEMS SO!d2 !n [DU]',
; * HCHO -->  title =  'GEMS HCHO [x 10!e16!n molecule/cm !e2 !n]',
; * CHOCHO --> title =  'GEMS CHOCHO [x 10!e16!n molecule/cm !e2 !n]',
; * O3T--> title= 'GEMS O3T !n [DU]',
; * Trop O3T--> title= 'GEMS TropO3 !n [DU]',
; * AOD,SSA, ECF, CRF -->  title= 'GEMS AOD 443nm',title= 'GEMS ECF'
; * ALH, AEH -->  title= 'GEMS AEH [km]',
; * Surface Reflectance -->  title= 'GEMS SR',
; * ECP --> title =  'GEMS ECP [mb]'
; * UVI --> title =  'GEMS UVI'
;----------------------------------------------------------------

cgPS_Close

if not keyword_Set(outfile) then begin
  outfile = 'test.png'
endif

cgPS2Raster,'gems.ps', outfile, /png, density = 1000

if not keyword_Set(scp_dest1) then begin
  scp_dest = 'soodal@164.125.38.179:/home/soodal/works/plot'
endif else begin
  scp_dest = scp_dest1
endelse

; send image to pc
if keyword_set(scp) then begin
  spawn, 'scp -P18742 -p ' + outfile + $
  ' ' + scp_dest
endif


end
