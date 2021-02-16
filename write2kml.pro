;file ='/data1/gems/o3t/mywork/val_omi_cross-section_2020-0616_bands.xdr'
file = '/data2/OMI/val_omi_cross-section_2020-0616_bands.xdr'
restore, file

nm_bands = ['310-340','300-340','305-340']
iband=2
case iband of
 0:begin
   o3prfs = data.bd 
   end

 1:begin
   o3prfs = data.bd1 
   end
 2:begin
   o3prfs = data.bd2
   end
endcase


 apri = data.apri
 gemsPlev = data.pres
 omiPlev = data.omiP

 sz = size(data.omip,/dim)  & xsz=sz[0] & psz=sz[1]
 pmid_omi = fltarr(xsz,psz-1)
 pmid_gems = fltarr(xsz,psz-1)

 FOR ic=0, xsz-1 DO BEGIN
   omip = reform(omiplev[ic,*]) 
   gemsp = reform(gemsPLev[ic,*]) 
   for iz=0,psz-2 do pmid_omi[ic,iz] =  exp( (alog(omip[iz])+alog(omip[iz+1]))/2 )
   for iz=0,psz-2 do pmid_gems[ic,iz] =  exp( (alog(gemsp[iz])+alog(gemsp[iz+1]))/2 )
 ENDFOR

 lon1 = 121 & lon2 = 155
 pos =[0.15,0.15,0.8,0.95]
 
 xtickv = [string(lon1,f='(I3)'), data.omilon[0:*:5], string(lon2,f='(I3)')]
 nxtic = n_elements(xtickv)-1
 nm_lonlat = string(data.omilon[0:*:5],f='(F5.1)') +'!C' + string(data.omilat[0:*:5],f='(F5.1)')
 xtickn = [' Lon: !C'+'Lat:',nm_lonlat, ' ']


 ylon = data.omilon[0]+(lon2-data.omilon[0])/7.*findgen(7) 
 ylat = interpol(data.omilat,data.omilon, ylon )
 xtickv = [string(lon1,f='(I3)'), ylon, string(lon2,f='(I3)')]
 nxtic = n_elements(xtickv)-1
 nm_lonlat = string(ylon,f='(F5.1)') +'!C' + string(ylat,f='(F5.1)')
 xtickn = [' lon: !C'+'lat:',nm_lonlat, ' ']

 cgname ='OMI_cross_section_2020_0616'

 ;cgps_open,'tmp.ps'
 ;cgdisplay,800,600
    ;title = 'OMI SAO'
    ;levels = findgen(21)*2
    ;;pomi = plot(fltarr(10,10, 10),/nodata,yr=[1000,0.1],/ylog,$
      ;/buffer, $
             ;;xr=[lon1,lon2],title=title,$
             ;;ytitle='Pressure [hPa]',xtitle='!CLongitude',$
             ;;ytickv = [1000,700,500,200,100,50,10,1,0.1],yticks=8,pos=pos,$
             ;;xtickv = xtickv,$
             ;;xtickname=xtickn, $
             ;;xticks=nxtic
           ;;)
    ;cgloadct,33,ncolors=22,bottom=0
    ;nlevels =21
    ;title = 'OMI SAO'
    ;levels = findgen(21)*2

    ;cgloadct,33,ncolors=22,bottom=0
    ;nlevels =21
    ;c = contour(data.omi,data.omilon,pmid_omi,/overplot,c_colors=indgen(21),$;bytscl(levels,min=levels[0],m
      ;/buffer, $
               ;levels=levels,$
              ;/cell_fill  )
 ;;colorbar

  ;cbpos = [0.82, 0.15, 0.85,0.95]
  ;cbxname = string( levels[indgen(10)*2],f='(I3)' )
  ;cfmt = '(I3)'
  ;charsize=1.3
  ;cb = Colorbar(Range=[min(levels),Max(levels)], $
  ;Divisions=10, XTicklen=1, XMinor=0, $
  ;AnnotateColor='black', NColors=nlevels, $
  ;Position=cbpos, $
  ;/buffer, $
  ;Title='O3 [DU]',Charsize=charsize,/right,/vertical,$
  ;format=cfmt,tickname=cbxname,TCHARSIZE=1,tlocation='top')
  
 ;cgps_close,/png
  ;rm_psfile,cgname  
 
 cgname ='GEMS_cross_section_2020_0616_'+nm_bands[iband]
 ;cgps_open,'tmp.ps'
 ;cgdisplay,800,600
    title = 'GEMS [' +nm_bands[iband]+ 'nm] '
    levels = findgen(21)*2
    pgems = plot(fltarr(10,10), /buffer);/nodata,yr=[1000,0.1],/ylog,$
             ;xr=[lon1,lon2],title=title,$
             ;ytitle='Pressure [hPa]',xtitle='!CLongitude',$
             ;ytickv = [1000,700,500,200,100,50,10,1,0.1],yticks=8,pos=pos,$
             ;xtickv = xtickv,xtickname=xtickn, xticks=nxtic
    cgloadct,33,ncolors=22,bottom=0
    m = map('Mollweide', limit=[80, -5, 150, 45])
    nlevels =21
    cgems = contour(o3prfs, data.omilon,pmid_gems, /buffer,$
      c_colors=indgen(21),$;bytscl(levels,min=levels[0],m
               levels=levels,$
              /cell_fill ) 
 ;colorbar
  cbpos = [0.82, 0.15, 0.85,0.95]
  cbxname = string( levels[indgen(10)*2],f='(I3)' )
  cfmt = '(I3)'
  charsize=1.3
  cbgems = Colorbar(Range=[min(levels),Max(levels)], $
  Divisions=10, XTicklen=1, XMinor=0, $
  AnnotateColor='black', $
  NColors=nlevels, $
  Position=cbpos, $
  /buffer, $
  Title='O3 [DU]',Charsize=charsize,/right,/vertical,$
  format=cfmt,tickname=cbxname,TCHARSIZE=1,tlocation='top')
  
cgems.save, 'gems.kml'
cgems.close
stop
 
 cgname ='Apriori_cross_section_2020_0616'
 ;cgps_open,'tmp.ps'
 ;cgdisplay,800,600
    title = 'ML Climatology'
    levels = findgen(21)*2
    pap = plot(fltarr(10,10),/nodata,yr=[1000,0.1],/ylog,$
      /buffer, $
             xr=[lon1,lon2],title=title);,$
             ;ytitle='Pressure [hPa]',xtitle='!CLongitude',$
             ;ytickv = [1000,700,500,200,100,50,10,1,0.1],yticks=8,pos=pos,$
             ;xtickv = xtickv,xtickname=xtickn, xticks=nxtic)
    cgloadct,33,ncolors=22,bottom=0
    nlevels =21
    cap = contour(data.apri,data.omilon,pmid_gems,/overplot,c_colors=indgen(21),$;bytscl(levels,min=levels[0],m
      /buffer, $
               levels=levels,$
              /cell_fill  )
 ;colorbar
  cbpos = [0.82, 0.15, 0.85,0.95]
  cbxname = string( levels[indgen(10)*2],f='(I3)' )
  cfmt = '(I3)'
  charsize=1.3
  cbap = Colorbar(Range=[min(levels),Max(levels)], $
  Divisions=10, XTicklen=1, XMinor=0, $
  AnnotateColor='black', NColors=nlevels, $
  Position=cbpos, $
  /buffer, $
  Title='O3 [DU]',Charsize=charsize,/right,/vertical,$
  format=cfmt,tickname=cbxname,TCHARSIZE=1,tlocation='top')
  
 ;cgps_close,/png
  ;rm_psfile,cgname  
  pap.save, 'apriori.kml'
  pap.close

end
