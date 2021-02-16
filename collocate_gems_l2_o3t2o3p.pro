pro collocate_gems_l2_o3t2o3p, o3t, result

;sz = size(o3p)
dim1 = 174
dim2 = 512

result = fltarr(dim1, dim2)


for j = 0, dim2-1 do begin ; spatial 512
for i = 0, dim1-1 do begin ; image 174
  if i eq dim1-1 then begin
    ;if total(finite(o3t[i*4:(i+1)*4-2, j*4:(j+1)*4-1], /nan)) ge 10 then begin
      ;result[i, j] = !values.f_nan
    ;endif else begin
      result[i, j] = mean(o3t[i*4:(i+1)*4-2, j*4:(j+1)*4-1], /nan)
    ;endelse
  endif else begin
    ;if total(finite(o3t[i*4:(i+1)*4-1, j*4:(j+1)*4-1], /nan)) ge 10 then begin
      ;result[i, j] = !values.f_nan
    ;endif else begin
      result[i, j] = mean(o3t[i*4:(i+1)*4-1, j*4:(j+1)*4-1], /nan)
    ;endelse
  endelse

endfor
endfor

return

end
