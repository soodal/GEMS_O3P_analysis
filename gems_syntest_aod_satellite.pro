;pro gems_syntest_aod

sjul = julday(01, 15, 2016)
ejul = julday(01, 15, 2016)

utc = '04'

for jul = sjul, ejul do begin

caldat, jul, mm, dd, yy
yyst = string(yy,'(I4.4)')
mmst = string(mm,'(I2.2)')
ddst = string(dd,'(I2.2)')

;;----------------------------
;; call retrieval value
;;----------------------------

path = '/home/myid0121/test_syn/vector/'
files = file_search(path+'GEMS_ASIAAOP_AODSSA_'+yyst+'m'+mmst+ddst+utc+'.sav',count = files_num)
restore, files(0)  ;; elon, elat, aodres


;;----------------------------
;; set dimension
;;----------------------------

	dim1 = 1199
	dim2 = 899

lon = replicate(!values.f_nan, dim1, dim2)
lat = replicate(!values.f_nan, dim1, dim2)
pic_aod = replicate(!values.f_nan, dim1, dim2)

lon = reform(elon, dim1, dim2)
lat = reform(elat, dim1, dim2)
pic_aod = reform(aodres, dim1, dim2)

;;----------------------------
;; plot procedure
;;----------------------------

pic_title =yyst+mmst+ddst+utc +  ' AOD'

x_para = elon           & X_name = 'Longitude'
y_para = elat           & Y_name = 'Latitude'

Slat = -5.
Nlat = 50.
Llon = 75.
Rlon = 150.

cgPS_Open, 'gems.ps'
cgDisplay, aspect=0.7
		LoadCT, 22, Ncolors=254,bottom=1

pos = [0.1,0.15,0.9,0.97]
charsize = (!D.Name EQ 'PS') ? cgDefCharsize()*0.7 : cgDefCharsize()*0.75
thick = (!D.Name EQ 'PS') ? 6 :3
londel = 10
latdel = 10


cgMap_set,36,128, Limit=[-60,70 , 60, 160],$
		           /noerase, /advance,pos=pos,charsize=charsize,$
		          /satellite,/noborder,/iso,$
		          title=title



		for j = 1L, dim2-2 do begin
			for i = 1L, dim1-2 do begin

				if (lat(i,j) ge Slat) and (lat(i,j) le Nlat) and (lon(i,j) ge Llon) and (lon(i,j) le Rlon)  then begin

				if pic_aod(i,j) gt 0.01 then begin
              	xbox = [lon(i,j), lon(i+1,j), lon(i+1,j), lon(i,j)]
                ybox = [lat(i,j), lat(i,j), lat(i,j+1), lat(i,j+1)]
				polyfill, xbox, ybox, color = bytscl(pic_aod(i,j), min = 0.000001 , max = 1.5,top = 255)
				endif
				endif


			endfor
		endfor



cgMap_Grid, Color='black', GLinestyle=1,londel=londel,latdel=latdel,charsize=1 $
			           , LONLAB=20, LATLAB=160, LABEL=1


cgMap_continents,color='black'


cgColorbar,format='(f6.1)',range=[0.0,1.5], Position=[pos[0]+0.065, pos[1]-0.04, pos[2]-0.065, pos[1]-0.02],color=0,$
							OOB_High=254,OOB_low=1,$
	                        divisions=5,bottom=1,ncolor=254,minor=1, charsize=1.5,$
	                        XTicklen=1, XMinor=0, $
 	     				    AnnotateColor='black', Title='GEMS AOD'


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
cgPS2Raster,'gems.ps', path+pic_title+'_global.png', /png, density = 1000


endfor ;; julday



stop

end
