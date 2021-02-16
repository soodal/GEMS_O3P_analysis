pro collocate_gemsl2o3p_omo3pr, gems, merra2, $
  gemsvals, merra2vals, on_under_300=on_under_300

if not keyword_Set(on_under_300) then begin
  print, 'Total Ozone'
  on_under_300 = 0
endif

limit=[-10, 80, 60, 160]

yyyy = string(year, format='(i04)')
mm = string(month, format='(i02)')
dd = string(day, format='(i02)')
hh = string(hour, format='(i02)')
mi = string(minute, format='(i02)')

savpath = './collocate_gemsl2o3p_merra2/'
if not file_test(savpath) then begin
  file_mkdir, savpath
endif

savfile = savpath + 'col_gems_merra2_'+yyyy+mm+dd+'_'+hh+mi+'.sav'

savefile =  file_test(savfile) 
;savefile = 0
if not savefile then BEGIN

  ; path for GEMS
  gemsfn = '../GEMS_O3P_Yonsei/out/310340_EOSRL/GK2B_GEMS_L2_20200616_0345_O3P_ND_DPRO_first.4x4_2008301839.nc4'

  ; read GEMS O3P data
  gemsvars = ds_read_gems_l2_o3p(gemsfn)

;;----------------------------
;; GEMS vertical column layer for under 300 hPa
;;----------------------------
  o3 = gemsvars.o3
  o3size = size(gemsvars.o3, /dim)
  ;altitude = gemsvars.Altitude
  gemspres = gemsvars.Pressure
  gemslayero3 = fltarr(o3size[1:2])
  dim1 = o3size[1]
  dim2 = o3size[2]
  ;gemslayero3[*,*] = !values.f_nan
  for iy=0, dim2-1 do begin
  for ix=0, dim1-1 do begin
    for ilevel=24, 0, -1 do begin
      ilayer = ilevel-1
      if gemspres[ilevel, ix, iy] gt 300. $
          and gemspres[ilevel-1, ix, iy] gt 300. then begin
        gemslayero3[ix, iy] = gemslayero3[ix, iy] + o3[ilayer, ix, iy]

      endif else if gemspres[ilevel, ix, iy] gt 300. $
          and gemspres[ilevel-1, ix, iy] le 300. then begin

        gemslayero3[ix, iy] = gemslayero3[ix, iy] + o3[ilayer, ix, iy]*$
          (gemspres[ilevel, ix, iy]-300.)/$
          (gemspres[ilevel, ix, iy]-gemspres[ilevel-1, ix, iy])
        break
      endif
    endfor
  endfor
  endfor
  gemslayero3[where(gemslayero3 lt 0)] = !values.f_nan

  ; GEMS vertical column layer for under 10 km
  ;for iy=0, dim2-1 do begin
  ;for ix=0, dim1-1 do begin
    ;for ilevel=24, 0, -1 do begin
      ;if altitude[ilevel, ix, iy] lt 10. $
          ;and altitude[ilevel-1, ix, iy] lt 10. then begin
        ;gemslayero3[ix, iy] = gemslayero3[ix, iy] + o3[ilevel-1, ix, iy]

      ;endif else if altitude[ilevel, ix, iy] lt 10. $
          ;and altitude[ilevel-1, ix, iy] ge 10. then begin

        ;gemslayero3[ix, iy] = gemslayero3[ix, iy] + o3[ilevel-1, ix, iy]*$
          ;(10.-altitude[ilevel, ix, iy])/$
          ;(altitude[ilevel-1, ix, iy]-altitude[ilevel, ix, iy])
        ;break
      ;endif
    ;endfor
  ;endfor
  ;endfor
  ;gemslayero3[where(gemslayero3 le 0)] = !values.f_nan

  ; collocation

  merra2lon = merra2.lon
  merra2lat = merra2.lat
  merra2time = merra2.time

  ds_merra2_get_o3underp, data, merra2under300, itime=1

  
  gemslon = gemsvars.Longitude
  gemslat = gemsvars.Latitude
  ;gemstime = gemsvars.Time
  gemslon[where(gemslon lt -180, /null)] = !values.f_nan
  gemslat[where(gemslat lt -90, /null)] = !values.f_nan
  ;gemstime[where(gemstime lt 0, /null)] = !values.f_nan
  ;gemstime = julday(month, day, year, hour, minute)
  gemssize = size(gemslon, /dim)
  
  ;yn_omi_cross_time = intarr(n_elements(omilon))
  ;yn_omi_cross_time[*] = 0
  ;yn_omi_cross_time[where(abs(omijulday-gemstime) lt 1./24./2., npix, /null)] = 1

  ;closest_idx = lonarr(gemssize[0], gemssize[1])
  ;closest_idx[*] = -999

  nmerrapix = n_elements(merra2lon)

  omisize = size(omilon, /dim)
  closest_idx = lonarr(omisize[0])
  closest_idx[*] = -999

  ;for iy = 0, gemssize[1]-1 do BEGIN
    ;for ix = 0, gemssize[0]-1 do BEGIN
      ;x = gemslon[ix, iy]
      ;y = gemslat[ix, iy]

      ;if finite(x) and finite(y) then begin
        ;result = search_closest_pixel(omilon, omilat, x, y)
        ;;result = array_indices(gemslon, result)
        ;closest_idx[ix, iy] = result
      ;ENDIF
    ;ENDFOR
  ;ENDFOR

  for ix = 0, nmerrapix-1 do BEGIN
    x = merra2lon[ix]
    y = merra2lat[ix]

    if finite(x) and finite(y) then begin
      result = search_closest_pixel(gemslon, gemslat, x, y)
      closest_idx[ix] = result
    ENDIF
  ENDFOR
  
  save, filename=savfile, gemsvars, gemslayero3, merra2, merra2under300, $
    closest_idx

  ;gemstoz = reform(gemsvars.ColumnAmountO3[0, *, *])
  ;gemsvals = []
  ;omivals = []
  ;for iy=0, gemssize[1]-1 do BEGIN
    ;for ix=0, gemssize[0]-1 do BEGIN
      ;x = gemslon[ix, iy]
      ;y = gemslat[ix, iy]

      ;if closest_idx[ix, iy] GE 0 then begin
        ;omival = omitoz[closest_idx[ix, iy]]
        ;gemsval = gemstoz[ix, iy]

        ;if omival ge 0 and omival le 1000 AND $
            ;gemsval ge 0 and gemsval le 1000 then begin
          ;omivals = [omivals, omival]
          ;gemsvals = [gemsvals, gemsval]
        ;endif

      ;endif
    ;ENDFOR
  ;ENDFOR

endif else BEGIN
  RESTORE, savfile
ENDELSE

  gemstoz = reform(gemsvars.ColumnAmountO3[0, *, *])
  gemsvals = []
  gemsecf = gemsvars.effectivecloudfractionuv
  for ix=0, omisize[0]-1 do BEGIN

    if closest_idx[ix] GE 0 then begin
          and gemsecf[closest_idx[ix]] le 0.2 $
        then begin

        merra2val = merra2under300[ix] 
        gemsval = gemstoz[closest_idx[ix]]

        ;print,merra2val, gemsval

        if merra2val ge 0 and merra2val le 1000 AND $
            gemsval ge 0 and gemsval le 1000 then begin
          merra2vals = [merra2vals, merra2val]
          gemsvals = [gemsvals, gemsval]
        endif
      endif

    endif
  ENDFOR
  ;plot_omi_satproj, omivars.lon, omivars.lat, omivars.toz,$
    ;filename='plot/omi_toz_20200806T0345.png', $
    ;title='OMI TOZ', range=[220, 350]
  ;plot_omi_satproj, omivars.lon, omivars.lat, omivars.rowanomaly,$
    ;filename='plot/OMI_rowanomaly_20200806T0345.png', $
    ;title='OMI Row Anomaly', range=[0, 1.]
  return
end
