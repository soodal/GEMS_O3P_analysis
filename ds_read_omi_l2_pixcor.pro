pro ds_read_omi_l2_pixcor, file, result
  dg1 = '/HDFEOS/SWATHS/OMI Ground Pixel Corners UV-1/Data Fields'
  gg1 = '/HDFEOS/SWATHS/OMI Ground Pixel Corners UV-1/Geolocation Fields'
             
  dg2 = '/HDFEOS/SWATHS/OMI Ground Pixel Corners UV-2/Data Fields'
  gg2 = '/HDFEOS/SWATHS/OMI Ground Pixel Corners UV-2/Geolocation Fields'

  dgvis = '/HDFEOS/SWATHS/OMI Ground Pixel Corners VIS/Data Fields'
  ggvis = '/HDFEOS/SWATHS/OMI Ground Pixel Corners VIS/Geolocation Fields'

  lat1 = h5read(file,gg1+'/Latitude')
  lon1 = h5read(file,gg1+'/Longitude')

  lat2 = h5read(file,gg2+'/Latitude')
  lon2 = h5read(file,gg2+'/Longitude')

  latvis = h5read(file,ggvis+'/Latitude')
  lonvis = h5read(file,ggvis+'/Longitude')

  corlat1 = h5read(file,dg1+'/FoV75CornerLatitude')  
  corlon1 = h5read(file,dg1+'/FoV75CornerLongitude')  

  corlat2 = h5read(file,dg2+'/FoV75CornerLatitude')  
  corlon2 = h5read(file,dg2+'/FoV75CornerLongitude')  
  
  corlatvis = h5read(file,dgvis+'/FoV75CornerLatitude')  
  corlonvis = h5read(file,dgvis+'/FoV75CornerLongitude')  

  sz1 = SIZE(lon1,/dim)
  sz2 = SIZE(lon2,/dim)
  szvis = SIZE(lonvis,/dim)

  tmp1 = h5read(file,gg1+'/Time')

  tmp2 = h5read(file,gg2+'/Time')

  tmpvis = h5read(file,ggvis+'/Time')

  time1 = dblarr(60,n_elements(tmp1))

  time2 = dblarr(60,n_elements(tmp2))

  timevis = dblarr(60,n_elements(tmpvis))

  for i=0,59 do time1[i,*] = tmp1

  for i=0,59 do time2[i,*] = tmp2

  for i=0,59 do timevis[i,*] = tmpvis

  caldat,(julday(1,1,1993,0,0,0)*86400.0d0+time1)/86400.0d0,mon1,day1,year1,hour1,min1,sec1

  caldat,(julday(1,1,1993,0,0,0)*86400.0d0+time2)/86400.0d0,mon2,day2,year2,hour2,min2,sec2

  caldat,(julday(1,1,1993,0,0,0)*86400.0d0+timevis)/86400.0d0,monvis,dayvis,yearvis,hourvis,minvis,secvis

  date1 = strarr(sz1[0],sz1[1])
  date2 = strarr(sz2[0],sz2[1])
  datevis = strarr(szvis[0],szvis[1])

  FOR iy=0, sz1[1]-1 DO BEGIN
    FOR ix=0, sz1[0]-1 DO BEGIN
      date1[ix,iy] =  string(year1[ix,iy],mon1[ix,iy],day1[ix,iy],format='(I4,I02,I02)')+'T'+$
                     string(hour1[ix,iy],min1[ix,iy],sec1[ix,iy],format='(I02,I02,I02)')+'Z'
    ENDFOR
  ENDFOR

  FOR iy=0, sz2[1]-1 DO BEGIN
    FOR ix=0, sz2[0]-1 DO BEGIN
      date2[ix,iy] =  string(year2[ix,iy],mon2[ix,iy],day2[ix,iy],format='(I4,I02,I02)')+'T'+$
                     string(hour2[ix,iy],min2[ix,iy],sec2[ix,iy],format='(I02,I02,I02)')+'Z'
    ENDFOR
  ENDFOR

  FOR iy=0, szvis[1]-1 DO BEGIN
    FOR ix=0, szvis[0]-1 DO BEGIN
      datevis[ix,iy] =  string(yearvis[ix,iy],monvis[ix,iy],dayvis[ix,iy],format='(I4,I02,I02)')+'T'+$
                     string(hourvis[ix,iy],minvis[ix,iy],secvis[ix,iy],format='(I02,I02,I02)')+'Z'
    ENDFOR
  ENDFOR
   
  result ={date1:date1, lat1:lat1, lon1:lon1, corlat1:corlat1, corlon1:corlon1, jtime1:time1, $
    date2:date2, lat2:lat2, lon2:lon2, corlat2:corlat2, corlon2:corlon2, jtime2:time2, $
    datevis:datevis, latvis:latvis, lonvis:lonvis, corlatvis:corlatvis, corlonvis:corlonvis, jtimevis:timevis}
 
end
