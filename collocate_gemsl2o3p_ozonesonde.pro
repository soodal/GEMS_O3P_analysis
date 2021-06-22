function collocate_gemsl2o3p_ozonesonde, lon, lat, sondelon, sondelat

distance = sqrt((lon-sondelon)^2 + (lat-sondelat)^2)

closest_idx = where(distance eq min(distance), /null)

return, closest_idx
end

