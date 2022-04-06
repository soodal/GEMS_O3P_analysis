function angle_minus_twopi, gamma0, pival, result

if ( gamma0 > pival ) then begin
  result = gamma0 - 2.0 * pival ;sign(2.0_r8*pival - gamma0, gamma0)
endif else if ( gamma0 lt -pival ) then begin
  result = gamma0 + 2.0 * pival
endif else begin
  result = gamma0
end if

return result
end 
