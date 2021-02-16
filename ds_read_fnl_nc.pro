
file = '/data1/gems/o3p/ds/GEMS_O3P_analysis/fnl_20200616_00_00.grib2.nc'
fid=ncdf_open(file)
ncdf_list, file, /variables, /dimensions, /gatt, /vatt

tid = ncdf_varid(fid, 'TMP_P0_L100_GLL0')
ncdf_varget, fid, tid, fnltemp


lvid = ncdf_varid(fid, 'lv_ISBL0')
ncdf_varget, fid, lvid, fnllevel


;;< this case for which_tprof = 2 
;pold0 = [1013.25, 1000., 950., 925., 850., $
  ;700., 600., 500., 400., 300., $
  ;250., 200., 150., 100., 70., $
  ;50., 30., 20., 10., 7., $
  ;5., 3., 2., 1., 0.50, $
  ;0.40, 0.25, 0.175, 0.125, 0.0874604]

pold0 = [1013.25, 1000., 975., 950., 925., 900., 850., $
        800., 750., 700., 650., 600., 550., 500., 450., $
        400., 350., 300., 250., 200., 150., 100., 70., 50., $
        30., 20., 10., 7., 5., 3., 2., 1., 0.70, 0.35, 0.25, $
        0.175, 0.125, 0.0874604]
              

fnltemp_interpol = fltarr(360, 180, 26)
FOR ix=0, 359 DO BEGIN
  FOR iy=0, 179 DO BEGIN
    fnltemp_interpol[ix, iy, *] = interpol(reform(fnltemp[ix, iy, *]), $
      fnllevel[*]/100., pold0[25:0:-1]);pold0[-5:-30:-1])
  ENDFOR
ENDFOR



filename = '/data1/app/gemsl2_2020.v0.2.i/data/o3p/ATMOS/fnl13.75LST/fnltemp/fnltemp_20200616.dat'
IF file_test(filename) THEN BEGIN
  file_delete, filename
ENDIF
openw, 1, filename

FOR ilayer=0, 25 DO BEGIN
  FOR iy=179, 0, -1 DO BEGIN
    printf, 1, fnltemp_interpol[*, iy, ilayer], format='(360i3)'
    ;print,mean(fnltemp_interpol[*, iy, ilayer])
  ENDFOR
  print, mean(fnltemp_interpol[*, *, ilayer])
ENDFOR
close, 1
                                                                                                                                     
stop

p1id=ncdf_groupsinq(fid)
p2id=ncdf_groupsinq(p1id[0])
p3id=ncdf_groupsinq(p2id[0])
p4id=ncdf_groupsinq(p3id[0])

lat_id  = ncdf_varid(p4id[0],'Latitude')
res_id  = ncdf_varid(p4id[1],'ResidualsOfFit')
o3solerr_id = ncdf_varid(p4id[1], 'O3SolutionError')

;ncdf_varget,p4id[0],lat_id,Latitude
ncdf_varget,p4id[1],res_id,ResidualsOfFit
ncdf_varget,p4id[1],o3solerr_id,O3SolutionError

;latdims = size(Latitude, /dim)
resdims = size(ResidualsOfFit,/dim)
o3solerrdims = size(O3SolutionError,/dim)

nl = o3solerrdims[2]

nres = resdims[0]
npix = resdims[1] 
nline = resdims[2]


;data_list = ['EffectiveCloudFractionUV', 'ProcessingQualityFlags', $
             ;'RootMeanSquareErrorOfFit', 'AveragingKernel', 'O3' , $
             ;'O3Apriori', 'O3AprioriError','ColumnAmountO3', $
             ;'CloudPressure','ResidualsOfFit','NumberOfIterations','Runtime',
             ;'SimulatedRadiances'] ; Added by Dae Sung 2020-08-29

;geo_list  = ['Latitude' ,'Longitude', 'SolarZenithAngle', $
             ;'ViewingZenithAngle', 'Time','Altitude' ,    $
             ;'Pressure','Pix','Line','TropopausePressure']

str = {Altitude:fltarr(nl+1), Pressure:fltarr(nl+1), $
  AveragingKernel:fltarr(nl, nl), O3:fltarr(nl), O3Apriori:fltarr(nl), $
  Runtime:0.0, O3AprioriError:fltarr(nl), ColumnAmountO3:fltarr(3),CloudPressure:0.0, $ 
  Time:0.0d, Latitude:0.0, Longitude:0.0, SolarZenithAngle:0.0, $
  ViewingZenithAngle:0.0, EffectiveCloudFractionUV:0.0, ProcessingQualityFlags:0, $
  RootMeanSquareErrorOfFit:0.0, Line:0, Pix:0,TropopausePressure:0.0, $
  ResidualsOfFit:fltarr(nres), SimulatedRadiances:fltarr(nres), O3SolutionError:fltarr(nl), $
  O3RandomNoiseError:fltarr(nl), DegreesOfFreedomForSignal:fltarr(3), $
  NumberOfIterations:0, FinalAlgorithmFlags:0}

gems=replicate(str,npix,nline)

tmpdata = fltarr(npix,nline)



end
