pro read_l2_pixcor, file,$
               result
               
 dg = '/HDFEOS/SWATHS/OMI Ground Pixel Corners UV-2/Data Fields'
 gg = '/HDFEOS/SWATHS/OMI Ground Pixel Corners UV-2/Geolocation Fields'

  lat = h5read(file,gg+'/Latitude')
  lon = h5read(file,gg+'/Longitude')

  corlat = h5read(file,dg+'/FoV75CornerLatitude')  
  corlon = h5read(file,dg+'/FoV75CornerLongitude')  
  
  sz = SIZE(lon,/dim)

  tmp = h5read(file,gg+'/Time')
  time = dblarr(60,n_elements(tmp))
  for i=0,59 do time[i,*] = tmp
  caldat,(julday(1,1,1993,0,0,0)*86400.0d0+time)/86400.0d0,mon,day,year,hour,min,sec

  date = strarr(sz[0],sz[1])
  FOR iy=0, sz[1]-1 DO BEGIN
  FOR ix=0, sz[0]-1 DO BEGIN
    date[ix,iy] =  string(year[ix,iy],mon[ix,iy],day[ix,iy],format='(I4,I02,I02)')+'T'+$
                   string(hour[ix,iy],min[ix,iy],sec[ix,iy],format='(I02,I02,I02)')+'Z'
  ENDFOR & ENDFOR
   
 result ={date:date, lat:lat, lon:lon, corlat:corlat, corlon:corlon, jtime:time}
 
  end
