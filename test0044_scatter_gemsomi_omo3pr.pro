
outputpath = './plot/'

gemsomipath = '/home/geun/2_O3PR/tool_for_geun/gemso3p_analysis/monthly_yonform/total_out/'
omo3prpath = '/home/Data/OMI/6_OML2O3P/2005/'
profozpath = '/home/Data/OMI/7_OML2PROFOZ/2005/'


searchstr = 'GEM_TEST_L2_O3P_2005m????_*.h5'
filelist = file_search(gemsomipath + searchstr)
nfile = n_elements(filelist)
fn = outputpath + 'GEMS_O3P_from_OMI_DATA_2005_yearly_scatterplot.png'

;;--------------------------------------------------
;; Make grid for OMI 0.5 grid
;;--------------------------------------------------

	Slat = -5.
	Nlat = 45.
	Llon = 75.
	Rlon = 145.

;; * MAKE grid

	res = 0.5
	nx = (Rlon-Llon)/res+1.
	ny = (Nlat-Slat)/res+1.
	LAT = fltarr(nx,ny)
	LON = fltarr(nx,ny)
	for j = 0, ny-1 do lat(*,j) = Nlat-res*j
	for i = 0, nx-1 do lon(i,*) = Llon+res*i

;;--------------------------------------------------

;;----------------------------------------
;; GEMS O3P for 2005
;;----------------------------------------

  gemsomi_o3p_col = []
  omo3pr_o3p_col = []



  savefile = './gems_omi_yearly_scatter.sav'
  if not file_test(savefile) then begin
    for ifile = 0L, nfile-1. do begin
      print, ifile
      sdate = strmid(filelist[ifile], 91,9)
      stime = strmid(filelist[ifile], 101, 6)


      ;profozfilelist = file_search(profozpath + '*'+sdate+ '*' + stime + '*')
      ;help, profozfilelist
      ;if n_elements(profozfilelist) eq 1 then begin
        ;ds_read_omi_l2_profoz, profozfilelist[0], profoz
      ;endif
      gemsomi_o3p = ds_read_gems_from_omi_data(filelist[ifile])

      gemsomi_stacking = replicate(0., nx, ny)
      gemsomi_stack_num = replicate(0, nx, ny)
      gemsomi_orbit = replicate(!values.f_nan, nx, ny)

      gemsomi_o3 = gemsomi_o3p.o3
      gemsomi_lon = gemsomi_o3p.longitude
      gemsomi_lat = gemsomi_o3p.latitude
      gemsomi_alt = gemsomi_o3p.altitude ; KM
      gemsomi_pres = gemsomi_o3p.pressure ; hPa
      gemsomi_o3 = transpose(gemsomi_o3, [1, 2, 0])
      gemsomi_alt = transpose(gemsomi_alt, [1, 2, 0])
      gemsomi_pres = transpose(gemsomi_pres, [1, 2, 0])
      gemsomi_o3p_1 = create_struct('O3', gemsomi_o3)
      gemsomi_o3p_1 = create_struct(gemsomi_o3p_1, 'Altitude', gemsomi_alt)
      gemsomi_o3p_1 = create_struct(gemsomi_o3p_1, 'Pressure', gemsomi_pres)
      nanidx = where(gemsomi_lat lt -900 or gemsomi_lon lt -900 or gemsomi_o3 lt -900, /null)
      gemsomi_o3[nanidx] = !values.f_nan

      roiidx = where( (gemsomi_lat ge Slat) and (gemsomi_lat le Nlat) and $
        (gemsomi_lon ge Llon) and (gemsomi_lon le Rlon), roipixnum, /null)
      ds_gems_from_omi_l2o3p_accum, gemsomi_o3p_1, o3p_accum, height=10
      if roipixnum gt 0 then begin
          for jj = 0, ny-1 do begin
          yidx = where(abs(gemsomi_lat-(Nlat-jj*res)) lt res*0.5,ynum)
            if ynum gt 0 then begin
            for ii = 0, nx-1 do begin
              xidx = where(abs(gemsomi_lon[yidx]-(Llon+ii*res)) lt res*0.5,xnum)
              if xnum gt 0 then begin

                data_o3p = total(o3p_accum[yidx[xidx]], /nan)
                pixnum = total(finite(o3p_accum[yidx[xidx]]))
                gemsomi_stacking[ii, jj] = gemsomi_stacking[ii, jj] + data_o3p
                gemsomi_stack_num[ii, jj] = gemsomi_stack_num[ii, jj] + pixnum
              endif
            endfor
            endif
          endfor
      endif
      
      omo3prfilelist = file_search(omo3prpath + '*'+sdate+ '*' + stime + '*')
      help, omo3prfilelist
      if n_elements(omo3prfilelist) eq 1 then begin
        ds_read_omi_l2_omo3pr, omo3prfilelist[0], omo3pr
      endif

      omo3pr_stacking = replicate(0., nx, ny)
      omo3pr_stack_num = replicate(0, nx, ny)
      omo3pr_orbit = replicate(!values.f_nan, nx, ny)

      omo3pr_o3 = omo3pr.o3
      omo3pr_lon = omo3pr.longitude
      omo3pr_lat = omo3pr.latitude
      omo3pr_alt = omo3pr.altitude ; KM
      omo3pr_pres = omo3pr.pressure ; hPa
      omo3pr_o3 = transpose(omo3pr_o3, [1, 2, 0])
      omo3pr_alt = transpose(omo3pr_alt, [1, 2, 0])
      omo3pr_pres = transpose(omo3pr_pres, [1, 2, 0])
      omo3pr_o3p_1 = create_struct('O3', omo3pr_o3)
      omo3pr_o3p_1 = create_struct(omo3pr_o3p_1, 'Altitude', omo3pr_alt)
      omo3pr_o3p_1 = create_struct(omo3pr_o3p_1, 'Pressure', omo3pr_pres)
      nanidx = where(omo3pr_lat lt -900 or omo3pr_lon lt -900 or omo3pr_o3 lt -900, /null)
      omo3pr_o3[nanidx] = !values.f_nan
      roiidx = where( (omo3pr_lat ge Slat) and (omo3pr_lat le Nlat) and $
        (omo3pr_lon ge Llon) and (omo3pr_lon le Rlon), roipixnum, /null)
      ds_omo3pr_accum, omo3pr_o3p_1, o3p_accum, height=10
      if roipixnum gt 0 then begin
        for jj = 0, ny-1 do begin
          yidx = where(abs(omo3pr_lat-(Nlat-jj*res)) lt res*0.5,ynum)
          if ynum gt 0 then begin
            for ii = 0, nx-1 do begin
              xidx = where(abs(omo3pr_lon[yidx]-(Llon+ii*res)) lt res*0.5,xnum)
              if xnum gt 0 then begin

                data_o3p = total(o3p_accum[yidx[xidx]], /nan)
                pixnum = total(finite(o3p_accum[yidx[xidx]]))
                omo3pr_stacking[ii, jj] = omo3pr_stacking[ii, jj] + data_o3p
                omo3pr_stack_num[ii, jj] = omo3pr_stack_num[ii, jj] + pixnum
              endif
            endfor
          endif
        endfor
      endif

      gemsomi_orbit = gemsomi_stacking/gemsomi_stack_num
      numzeroidx = where(gemsomi_stack_num eq 0, /null)
      ;notzeroidx = where(gemsomi_stack_num gt 0, /null)
      gemsomi_orbit[numzeroidx] = !values.f_nan

      omo3pr_orbit = omo3pr_stacking/omo3pr_stack_num
      numzeroidx = where(omo3pr_stack_num eq 0, /null)
      ;notzeroidx = where(omo3pr_stack_num gt 0, /null)
      omo3pr_orbit[numzeroidx] = !values.f_nan

      notzeroidx = where(omo3pr_stack_num gt 0 and gemsomi_stack_num gt 0, /null)

      gemsomi_o3p_col = [gemsomi_o3p_col, gemsomi_orbit[notzeroidx]]
      omo3pr_o3p_col = [omo3pr_o3p_col, omo3pr_orbit[notzeroidx]]

    endfor  ;for nfile-1


    
    save, filename=savefile, $
      omo3pr_o3p_col, gemsomi_o3p_col
  endif else begin
    restore,savefile
  endelse

  stop
;;----------------------------------------
;; * plot procedure

pic_title = 'GEMS TropO3 [DU]'

x_para = LON           & X_name = 'Longitude'
y_para = LAT           & Y_name = 'Latitude'

cgPS_Open, 'gems1.ps'
cgDisplay, aspect=0.7
LoadCT, 33, Ncolors=254, Bottom=1

position = [0.15, 0.10, 0.85, 0.80]
thick = (!D.Name EQ 'PS') ? 6 :3

cgPlot, LON, LAT, /NoData, $
  XRange = [Llon, Rlon], YRange = [slat, nlat], $
  xstyle=1, ystyle=1, xminor=1, yminor=1, $
  Position=position, Charsize=1.5,XThick=3.0, YThick=3.0, Thick=6.0, /NoErase, $
  ytickname = ['0'+cgsymbol('deg')+'N','10'+cgsymbol('deg')+'N',$
  '20'+cgsymbol('deg')+'N','30'+cgsymbol('deg')+'N','40'+cgsymbol('deg')+'N'], $
  xtickname =  ['80'+cgsymbol('deg')+'E','90'+cgsymbol('deg')+'E',$
  '100'+cgsymbol('deg')+'E','110'+cgsymbol('deg')+'E',$
  '120'+cgsymbol('deg')+'E', '130'+cgsymbol('deg')+'E','140'+cgsymbol('deg')+'E']


for j = 1L, ny-2 do begin
for i = 1L, nx-2 do begin

  xbox = [lon(i,j), lon(i+1,j), lon(i+1,j), lon(i,j)]
  ybox = [lat(i,j), lat(i,j), lat(i,j+1), lat(i,j+1)]
  polyfill, xbox, ybox, color = bytscl(gemsomi_orbit[i,j], min = 0.0 , max = 60.0,top = 255)   ;;set color bar

endfor
endfor


cgMap_Set, /NoErase, /Cylindrical, Limit=[Slat, Llon, Nlat, Rlon], /NoBorder, Color=0, $
  POSITION = [!X.WINDOW[0], !Y.WINDOW[0], !X.WINDOW[1], !Y.WINDOW[1]]

cgMap_Continents, /Coasts, Color=0
cgMap_Grid, charsize=1.5

cgColorbar,format='(f6.2)',range=[0.0,60.0], position=[0.18, 0.87, 0.82, 0.89],color=0,$    ;;set color bar
          divisions=5,bottom=1,ncolor=254,minor=1, charsize=1.6, title = 'GEMS TropO3', Tcharsize=1.6, Tlocation = 'top'
                     ;divisions=5,bottom=1,ncolor=254,minor=1, charsize=1.6, title = 'GEMS SO!d2 !n [x 10!e16!n molecule/cm !e2 !n]', Tcharsize=1.6, Tlocation = 'top'

;-----------------------------------------------------------------------
; * colorbar title
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
; * ECP --> title =  'GEMS ECP [mb]'
; * UVI --> title =  'GEMS UVI'

;----------------------------------------------------------------



cgPS_Close

outfile = './plot/gems_from_omi_2005_yearly.png'

cgPS2Raster,'gems1.ps', outfile, /png, density = 1000


if not keyword_Set(scp_dest1) then begin
scp_dest = 'soodal@164.125.38.179:/home/soodal/works/plot'
endif else begin
scp_dest = scp_dest1
endelse

; send image to pc
spawn, 'scp -P18742 -p ' + outfile + $
' ' + scp_dest



end
