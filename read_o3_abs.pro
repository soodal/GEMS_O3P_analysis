pro read_o3_abs, wav, refabs1, refabs2, refabs3, filename = filename



if not keyword_set(filename) then begin
  path = '/OPER/SYSTEMS/ALGRTH/GEMS/L2/data/o3p/tbl/'
  fn = "o3abs_brion_uv_vacfinal.dat"
  filename = path + fn
endif 

if not file_test(filename) then begin
  print, 'o3 abs file is not exist. check the file path and set the keyword filename.'
endif

;openr, lun, /get_lun
;readf, lun
;outdata = 


readcol, filename, wav, refabs1, refabs2, refabs3, skipline=10

end
