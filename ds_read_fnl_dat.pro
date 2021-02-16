function ds_read_fnl_dat, filename, nx, ny, nl

;filename = '/data1/app/gemsl2_2020.v0.2.i/data/o3p/ATMOS/fnl13.75LST/fnltemp/fnltempavg06.dat'
openr, lun, filename, /get_lun

;x=findgen(360, start=-180)
;y=findgen(180)-90

;xx = rebin(x, [360, 180])
;yy = rebin(transpose(y), [360, 180])

fnltemp = fltarr(nx, ny, nl)
line = fltarr(360)
FOR ilayer=0, nl-1 DO BEGIN
  FOR iy=0, ny-1 DO BEGIN
    readf, lun, line, format='('+strtrim(string(nx),2)+'i3)'
    fnltemp[*, iy, ilayer] = line
  ENDFOR
ENDFOR
free_lun, lun
return, fnltemp

end
