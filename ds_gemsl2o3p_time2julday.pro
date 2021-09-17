function ds_gemsl2o3p_time2julday, gemsl2o3p_time

compile_opt idl2

jday = double(gemsl2o3p_time)/3600./24. + julday(1,1,2000, 12, 0)

nanidx = where(jday LT 0)
jday[nanidx] = !values.d_nan

return, jday
end
