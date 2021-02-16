
FUNCTION read_modis_mod03, file

  ;file = '/home/Data2/5_MODIS/MYD06_L2.A2017213.0325.006.2017213164720.hdf'
;  file ='/home/Data2/5_MODIS/MYD03/MYD03.A2007109.1330.006.2012073205151.hdf'
  fid = eos_sw_open(file,/read)
  sid = eos_sw_attach(fid, 'MODIS_Swath_Type_GEO') 
  status = eos_sw_readfield(sid,'Latitude', lat)
  status = eos_sw_readfield(sid,'Longitude', lon)
  status = eos_sw_readfield(fid,'SD Start Time', time)
  status = eos_sw_detach(sid)
  status = eos_sw_close(fid)

  fid = HDF_SD_START(file,/READ)
  icf =  HDF_SD_NAMETOINDEX(fid, 'SD start time')
  sid=HDF_SD_SELECT(fid, icf)
  HDF_SD_GETDATA, sid,time
  HDF_SD_ENDACCESS, sid

  sz = SIZE(lon,/dim)
  nx = sz[0]
  ny = sz[1]
  date = strarr(sz[0],sz[1])
  times = dblarr(sz[0],sz[1])
  caldat,(julday(1,1,1993,0,0,0)*86400.0d0+time)/86400.0d0,mon,day,year,hour,min,sec

  FOR iy=0, n_elements(time)-1 DO BEGIN
    dat  =  string(year[iy],mon[iy],day[iy],format='(I4,I02,I02)')+'T'+$
            string(hour[iy],min[iy],sec[iy],format='(I02,I02,I02)')+'Z'
    date[*,10*iy:iy*10+10 -1] = dat
    times[*,10*iy:iy*10+10-1] = time[iy]
  ENDFOR 

   result = {lon:lon,lat:lat,jtime:times}  
   return,result
END
