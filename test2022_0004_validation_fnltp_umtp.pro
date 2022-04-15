
fn = '/data/nier_ftp/DEV2/O3P/V1.0.3/202103/29/GK2_GEMS_L2_20210329_0345_O3P_FW_DPRO_BIN4x4.nc'
o3p = ds_read_gems_l2_o3p(fn)

o3p_lat = o3p.latitude
o3p_lon = o3p.longitude

nanidx = where(o3p_lat lt -1e29, /null)
o3p_lat[nanidx] = !values.f_nan

nanidx = where(o3p_lon lt -1e29, /null)
o3p_lon[nanidx] = !values.f_nan

ds_convert_fnl_to_dat_for_gemso3p, '20210329', '/data/MODEL/FNL/2021/',fnlout, outdat='/OPER/SYSTEMS/ALGRTH/GEMS/L2/data/o3p/ATMOS/fnl13.75LST/fnltp/fnltp/fnltp_20210329.dat'

fnl_lon = findgen(360) - 180
fnl_lon = rebin(fnl_lon, [360, 180])

fnl_lat = findgen(180)-90
fnl_lat = rebin(fnl_lat, [180, 360])
fnl_lat = transpose(fnl_lat)


gemsnum = n_elements(o3p_lat)

gemsRoiIdx = lonarr(gemsnum)
gemsRoiIdx[*] = -999

;omi_o3p_on_gems_grid = fltarr(gems_o3_size[0], gems_o3_size[1], gems_o3_size[2])
;omi_o3p_on_gems_grid[*] = !values.f_nan
;omi_o3p_accum_on_gems_grid = fltarr(gems_o3_size[0], gems_o3_size[1])
;omi_o3p_accum_on_gems_grid[*] = !values.f_nan


fnl_size = size(fnl_lon, /dimension)

gems_tp_on_fnl_grid = fltarr(fnl_size[0], fnl_size[1])
gems_tp_on_fnl_grid[*] = !values.f_nan

for ipix = 0, gemsnum-1 do BEGIN
  x = o3p_lon[ipix]
  y = o3p_lat[ipix]

  if finite(x) eq 1 and finite(y) eq 1 then begin
    GemsPixonFNLidx = search_closest_pixel(fnl_lon, fnl_lat, x, y, maxlimit=0.5)
    if n_elements(GemsPixonFNLidx) gt 0 and GemsPixonFNLidx ne -999 then begin
      ;GemsPixOnFNLIndices = array_indices(fnl_lon, GemsPixonFNLidx)

      gemsRoiIdx[ipix] = GemsPixonFNLidx

      ;omi_o3p_on_gems_grid[OmiPixOnGemsIndices[0], OmiPixOnGemsIndices[1], *] =  omio3s[*, ipix]
      gems_tp_on_fnl_grid[gemspixonFNLidx] =  o3p.TropopausePressure[ipix]

    endif
  ENDIF

  ; time difference check
  ;if gemsRoiIdx[ipix] ge 0 then begin
    ;if abs(omijulday[ipix] - gemsjulday[gemsRoiIdx[ipix]]) lt 1./24./2. then BEGIN
      ;yn_omi_cross_time[ipix] = 1
    ;ENDIF
  ;endif
ENDFOR

plot_gems_validation, gems_tp_on_fnl_grid, fnlout.tp, $
  filename='./plot/fnl_tp/validation_tp_fnl_um_20210329.png', $
  xtitle='UM calculated Tropopause Pressure', $
  ytitle='FNL Tropopause Pressure', $
  range=[50, 400], $
  delta=2.0


end

