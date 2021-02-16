
FUNCTION read_modis_myd06, file

  ;file = '/home/Data2/5_MODIS/MYD06_L2.A2017213.0325.006.2017213164720.hdf'

  varnames=['Cloud_Fraction', 'Cloud_Optical_Thickness', 'Cloud_Phase_Infrared', 'Cloud_Top_Pressure']
  vars=['cf', 'cot', 'cphase', 'ctp']
  nvar=N_ELEMENTS(varnames)

  fid=HDF_SD_START(file, /READ)
  FOR ivar=0,nvar-1 DO BEGIN 
    icf = HDF_SD_NAMETOINDEX(fid, varnames[ivar])
    sid=HDF_SD_SELECT(fid, icf)
    fillvalue_index=HDF_SD_ATTRFIND(sid, '_FillValue')
    valid_range_index=HDF_SD_ATTRFIND(sid, 'valid_range')
    scalefactor_index=HDF_SD_ATTRFIND(sid, 'scale_factor')
    offset_index=HDF_SD_ATTRFIND(sid, 'add_offset')
    HDF_SD_ATTRINFO, sid, fillvalue_index, data=fillvalue
    HDF_SD_ATTRINFO, sid, valid_range_index, data=valid_range
    HDF_SD_ATTRINFO, sid, scalefactor_index, DATA=scale_factor
    HDF_SD_ATTRINFO, sid, offset_index, data=offset
    HDF_SD_GETDATA, sid,data 
    HDF_SD_ENDACCESS, sid
    nan=WHERE(data EQ fillvalue[0], nnan)
    data=data*scale_factor[0]+offset[0]
    IF nnan NE 0 THEN data[nan]=-999 
    res=execute(vars[ivar]+'=data')
  ENDFOR  ; ivar

  fid = eos_sw_open(file,/read)
  sid = eos_sw_attach(fid, 'mod06') 
  status = eos_sw_readfield(sid,'Latitude', lat)
  status = eos_sw_readfield(sid,'Longitude', lon)
  status = eos_sw_readfield(sid,'Scan_Start_Time', time)
  status = eos_sw_detach(sid)
  status = eos_sw_close(fid)
  sz = SIZE(lon,/dim)
  nx = sz[0]
  ny = sz[1]
  date = strarr(sz[0],sz[1])
  caldat,(julday(1,1,1993,0,0,0)*86400.0d0+time)/86400.0d0,mon,day,year,hour,min,sec
  FOR iy=0, ny-1 DO BEGIN
  FOR ix=0, nx-1 DO BEGIN
    date[ix,iy] =  string(year[ix,iy],mon[ix,iy],day[ix,iy],format='(I4,I02,I02)')+'T'+$
                   string(hour[ix,iy],min[ix,iy],sec[ix,iy],format='(I02,I02,I02)')+'Z'
  ENDFOR & ENDFOR

  result = {nxtrack:nx, nline:ny, lon:lon, lat:lat, cot:cot, cf:cf, cphase:cphase, ctp:ctp, jtime:time}
  
  RETURN, result

END
