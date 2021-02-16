pro plot_gems_satproj_data, lon, lat, data, $
  range=range, filename=filename, $
  scp_dest=scp_dest1, $
  title=title


datetime_offset = 33
yyyy = strmid(filename, 16 + datetime_offset, 4, /reverse)
mm = strmid(filename, 12 + datetime_offset, 2, /reverse)
dd = strmid(filename, 10 + datetime_offset, 2, /reverse)
utc = strmid(filename, 7 + datetime_offset, 4, /reverse)

;wloffset = 55
;precoffset = 50
;initwl = strmid(filename, 16 + wloffset, 3, /reverse)
;prec = '0.00' + strmid(filename, 16 + wloffset-6, 1, /reverse)
;print, yyyy, mm, dd, utc, initwl, prec


do_plot_cao3 = [1, 1, 1]
do_plot_p300 = 0
do_plot_10km = 1


;;----------------------------
;; call retrieval value
;;----------------------------
outpath = './plot/'
;files = file_search(path+'GEMS_ASIAAOP_AODSSA_'+yyyy+'m'+mm+dd+utc+'.sav',count = files_num)
;restore, files(0)  ;; elon, elat, aodres


;;----------------------------
;; set dimension
;;----------------------------

;dim1 = 1199
;dim2 = 899
sz = size(lon)
dim1 = sz[1]
dim2 = sz[2]

;lon = replicate(!values.f_nan, dim1, dim2)
;lat = replicate(!values.f_nan, dim1, dim2)
;pic_aod = replicate(!values.f_nan, dim1, dim2)

;lon = reform(elon, dim1, dim2)
;lat = reform(elat, dim1, dim2)
;pic_aod = reform(aodres, dim1, dim2)





;;----------------------------
;; 
;;----------------------------
columnname = 'Tropospheric'
if not keyword_set(range) then begin
  range = [20, 60]
endif

;pic_title =yyyy+'-'+mm+'-'+dd+'_'+utc+'UTC'+ columnname +initwl + '-340nm_' +prec 
if not keyword_Set(title) then begin
  title =yyyy+'-'+mm+'-'+dd+'_'+utc+'UTC'+ columnname 
endif
  

Slat = -5.
Nlat = 60.
Llon = 80.
Rlon = 150.

cgPS_Open, 'gems.ps'
cgDisplay, aspect=0.7

pos = [0.1,0.15,0.9,0.90]
charsize = (!D.Name EQ 'PS') ? cgDefCharsize()*0.7 : cgDefCharsize()*0.75
thick = (!D.Name EQ 'PS') ? 6 :3
londel = 10
latdel = 10

cgMap_set,26,128, Limit=[-60, 0, 60, 360],$
              pos=pos,charsize=charsize,$
              /satellite,/iso,$
              /horizon, $
              title=title

;LoadCT, 22, Ncolors=254,bottom=1
dsloadct, 33

;for j = 1L, dim2-2 do begin
  ;for i = 1L, dim1-2 do begin
    ;if (lat[i,j] ge Slat) and (lat[i,j] le Nlat) and $
        ;(lon[i,j] ge Llon) and (lon[i,j] le Rlon) then begin
      ;if data[i,j] gt 0 and data[i, j] le 500 then begin
        ;xbox = [lon[i,j], lon[i+1,j], lon[i+1,j], lon[i,j]]
        ;ybox = [lat[i,j], lat[i,j], lat[i,j+1], lat[i,j+1]]
        ;polyfill, xbox, ybox, $
          ;color=bytscl(data[i,j], $
            ;min=range[0], $
            ;max=range[1], $
            ;top=253)
      ;endif
    ;endif
  ;endfor
;endfor
nanidx = where(finite(data, /nan) eq 1, /null)
lon[nanidx] = !values.f_nan
lat[nanidx] = !values.f_nan
plots, lon, lat, color=bytscl(data, min=range[0], max=range[1], top=253), psym=6, symsize=0.5

cgMap_Grid, Color='black', GLinestyle=1,londel=londel, latdel=latdel, charsize=1 $
                 ,LONLAB=20 ,LATLAB=160, LABEL=1

cgMap_continents, color='black'

cgColorbar, format='(f6.1)', range=range, $
  Position=[pos[0]+0.065, pos[1]-0.04, pos[2]-0.065, pos[1]-0.02], $
  Color='black',$
  OOB_High=253, $
  OOB_low=0,$
  divisions=5, $
  bottom=0, $
  ncolor=254, $
  minor=1,  $
  charsize=1.5,$
  XTicklen=1,  $
  XMinor=0, $
  AnnotateColor='black',  $
  Title='GEMS ' + columnname + ' Column O3'


cgPS_Close


pngfile = filename
cgPS2Raster,'gems.ps', pngfile, /png, density = 1000



if not keyword_Set(scp_dest1) then begin
  scp_dest = 'soodal@164.125.38.179:/home/soodal/works/plot'
endif else begin
  scp_dest = scp_dest1
endelse

; send image to pc
spawn, 'scp -P18742 -p ' + pngfile + $
  ' ' + scp_dest

end
