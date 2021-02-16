
for i =0, 18 do begin
  o3pval = total_o3pvals[174*i:174*(i+1)-1, *]
  o3tval = total_o3tvals[174*i:174*(i+1)-1, *]

  o3p_lon = total_lons[174*i:174*(i+1)-1, *]
  o3p_lat = total_lats[174*i:174*(i+1)-1, *]

  o3p_ecf = total_ecfs[174*i:174*(i+1)-1, *]

  oc = where( o3p_lon ge 80 and $
              o3p_lon le 150 and $
             o3p_lat ge -5 and $
             o3p_lat le 45 and $
            finite(o3pval) eq 1 and  $
            o3pval ge 50 and $
            o3p_ecf le 0.1, nloc )

  oc = where( o3p_lon ge 80 and $
              o3p_lon le 150 and $
             o3p_lat ge -5 and $
             o3p_lat le 45 and $
            finite(o3tval) eq 1 and  $
            o3tval ge 50 and $
            o3p_ecf le 0.1, nloct)

  oc1 = where(o3p_lon ge 80 and $
              o3p_lon le 150 and $
             o3p_lat ge -5 and $
             o3p_lat le 45 and $
            finite(o3pval) eq 1 and  $
            o3pval ge 50 and $
            o3p_ecf le 0.1, nloc1 )

  oc2 = where( o3p_lon ge 80 and $
              o3p_lon le 150 and $
             o3p_lat ge -5 and $
             o3p_lat le 45 and $
            finite(o3tval) eq 1 and  $
            o3tval ge 50 and $
            o3p_ecf le 0.1, nloct1)

            print, nloc, nloct, size(o3pval)
if i eq 1 then stop
endfor
end



