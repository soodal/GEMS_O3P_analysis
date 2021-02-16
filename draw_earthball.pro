

;x_para = elon           
;X_name = 'Longitude'
;y_para = elat           
;Y_name = 'Latitude'

Slat = -5.
Nlat = 60.
Llon = 75.
Rlon = 150.

cgPS_Open, 'gems.ps'
cgDisplay, aspect=0.7
;LoadCT, 22, Ncolors=254,bottom=1
LoadCT, 33, Ncolors=254,bottom=1

pos = [0.1,0.15,0.9,0.97]
charsize = (!D.Name EQ 'PS') ? cgDefCharsize()*0.7 : cgDefCharsize()*0.75
thick = (!D.Name EQ 'PS') ? 6 :3
londel = 10
latdel = 10

cgMap_set, /satellite, sat_p=[2., 0., 0.], 26, 128, $
              Limit=[-90, 0, 90, 360],$
              pos=pos, $
              ;/noborder, $
              charsize=charsize,$
              /iso,$
              /horizon, $
              e_horizon={fill:0, thick:10}
              ;title=title


cgMap_Grid, Color='black', GLinestyle=1,londel=londel, latdel=latdel, charsize=1 $
                 ,LONLAB=20 ,LATLAB=160, LABEL=0

cgMap_continents, color='black', thick=10
cgPS_Close


pngfile = 'earthball.png'
cgPS2Raster,'gems.ps', pngfile, /png, density = 1000

scp_dest = 'soodal@164.125.38.179:/home/soodal/WindowsHome/Pictures/scp'

; send image to pc
spawn, 'scp -P18742 -p ' + pngfile + $
  ' ' + scp_dest

end
