function ds_bytscl, input, min=minvalue, max=maxvalue, top=topvalue
x = input

if not keyword_set(minvalue) then begin
  minvalue = min(x, /nan)
endif

if not keyword_set(maxvalue) then begin
  maxvalue = max(x, /nan)
endif

if not keyword_set(topvalue) then begin
  topvalue = 255
endif

maxidx = where(x gt maxvalue, nmaxidx, /null)
if nmaxidx gt 0 then begin
  x[maxidx] = maxvalue
endif

minidx = where(x lt minvalue, nminidx, /null)
if nminidx gt 0 then begin
  x[minidx] = minvalue
endif

result = bytscl(x, min=minvalue, max=maxvalue, top=topvalue) 

return, result

end

