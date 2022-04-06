; plotting every scene aug to oct for gems l2 o3p
;

scp_dest = 'soodal@164.125.38.179:/home/soodal/works/plot/paper/2021_initial_report'

;initial variable setting
outputpath = './plot/eosrl_test/'

path = '/data2/L2_GEMS/L1C_test/o3p/'
;path_nier = '/data2/L2_GEMS/L1C_test/o3p/4x4/NIER/'
path_nier = '/data/private/soodal/ln/GEMS/L2O3P/'

;varlist = ['EffectiveCloudFractionUV', 'ProcessingQualityFlags', $
           ;'O3' , $
           ;;'O3Apriori', 'O3AprioriError', $
           ;;'CloudPressure', $
           ;;'SimulatedRadiances', $
           ;;'O3Apriori', 'O3AprioriError',$
           ;'ColumnAmountO3', $
           ;'Latitude' ,'Longitude', $
           ;'Time','Altitude' ,    $
           ;'Pressure', $
           ;;'TropopausePressure', $
           ;;'Wavelengths', $
           ;'WavelengthsWholeRange']

; run

juldays = timegen(start=julday(3, 29, 2021, 1, 45), $
  final=julday(5, 2, 2021, 23, 45), units='Hours')

juldays = [julday(3, 29, 2021, 0, 45), julday(3, 29, 2021, 3, 45)]

sub = ['EOSRL']

total_o3accum_eosrl_500 = []
total_o3accum_nier_500 = []
total_o3accum_eosrl_300 = []
total_o3accum_nier_300 = []

for i=0, n_elements(juldays)-1 do begin

  jd = juldays[i]
  caldat, jd, month, day, year, hour, minute
  yyyy = string(year, format='(i04)')
  mm = string(month, format='(i02)')
  dd = string(day, format='(i02)')
  hh = string(hour, format='(i02)')
  mi = string(minute, format='(i02)')
  ;print, yyyy, mm, dd,'_', hh, mi

  datetime_str = yyyy + mm + dd+ '_' +  hh+ mi


  search_str1 = 'GK2_GEMS_L2_O3P_' + datetime_str + '_' + sub[0] + '_BIN4x4.nc'
  filelist_eosrl = file_search(path + sub[0] + '/' + search_str1)
  print, filelist_eosrl
  nfile1 = n_elements(filelist_eosrl)

  ;search_str2 = 'GK2_GEMS_L2_O3P_' + datetime_str + '_' + sub[1] + '_BIN4x4.nc'
  ;filelist_poly = file_search(path + sub[1] + '/' + search_str2)
  ;nfile2 = n_elements(filelist_poly)

  search_str3 = 'GK2_GEMS_L2_O3P_' + datetime_str + '_4x4.nc'
  filelist_nier = file_search(path_nier + '/' + search_str3)
  print, filelist_nier
  nfile3 = n_elements(filelist_nier)

  if total(strlen(filelist_eosrl)) ne 0 $
    and total(strlen(filelist_nier)) ne 0 then begin 
  ;if total(strlen(filelist_eosrl)) ne 0 then begin
    ;for inum = 0, nfile-1 do begin ; gems l2o3p file
      o3p_eosrl = ds_read_gems_l2_o3p(filelist_eosrl[0])
      ;o3p_poly = ds_read_gems_l2_o3p(filelist_poly[0])
      o3p_nier = ds_read_gems_l2_o3p(filelist_nier[0])
      ;mondate = strmid(filelist[inum], 12, 4, /reverse)
      utchhmm = strmid(filelist_eosrl[0], 7, 4, /reverse)

      ecf_eosrl = o3p_eosrl.EffectiveCloudFractionUV

      ecf02idx = where(ecf_eosrl gt 0.2, /null)
      ecfnanidx = where(ecf_eosrl lt 0, /null)

      ds_gems_l2o3p_accum, o3p_eosrl, o3under500hpa_eosrl, hpa=500.
      ;ds_gems_l2o3p_accum, o3p_poly, o3under500hpa_poly, hpa=500.
      ;ds_gems_l2o3p_accum, o3p_nier, o3under500hpa_nier, hpa=500.

      ds_gems_l2o3p_accum, o3p_eosrl, o3under300hpa_eosrl, hpa=300.
      ;ds_gems_l2o3p_accum, o3p_poly, o3under300hpa_poly, hpa=300.
      ds_gems_l2o3p_accum, o3p_nier, o3under300hpa_nier, hpa=300.

      longitude_eosrl = o3p_eosrl.longitude
      latitude_eosrl = o3p_eosrl.latitude

      ;longitude_poly = o3p_poly.longitude
      ;latitude_poly = o3p_poly.latitude
      
      longitude_nier = o3p_nier.longitude
      latitude_nier = o3p_nier.latitude

      ;nanidx = where(o3under500hpa_eosrl lt 0 and o3under500hpa_poly lt 0, /null)

      ;o3under500hpa_eosrl[ecf02idx] = !values.f_nan
      ;o3under500hpa_eosrl[ecfnanidx] = !values.f_nan
      ;o3under500hpa_eosrl[nanidx] = !values.f_nan

      ;nanidx = where(longitude_eosrl lt -500, /null)
      ;longitude_eosrl[nanidx] = !values.f_nan
      ;latitude_eosrl[nanidx] = !values.f_nan

      ;nanidx = where(latitude_poly lt -500, /null)
      ;longitude_poly[nanidx] = !values.f_nan
      ;latitude_poly[nanidx] = !values.f_nan

      ;nanidx = where(latitude_nier lt -500, /null)
      ;longitude_nier[nanidx] = !values.f_nan
      ;latitude_nier[nanidx] = !values.f_nan

      ;plot_sat_proj, o3under500hpa_eosrl - o3under500hpa_poly,
        ;longitude_eosrl, latitude_eosrl, $
        ;title='GEMS L2 O3P BTDF eosrl POLY DIFF SFC-500hPa ' + datetime_str + 'UTC', $
        ;range=[-10., 10.], $
        ;pngfile=outputpath + 'gems_l2_o3p_' + datetime_str +
        ;'_under500hpa_btdf_eosrl_poly_diff.png';, $
        ;/scp_send, $
        ;scp_dest=scp_dest

      ;valididx = where(o3under500hpa_nier eq 0 or o3under500hpa_eosrl eq 0, /null)
      valididx = where(o3under500hpa_eosrl eq 0, /null)
      o3under500hpa_eosrl[valididx] = !values.f_nan
      ;o3under500hpa_nier[valididx] = !values.f_nan

      valididx = where(o3under300hpa_nier eq 0 or o3under300hpa_eosrl eq 0, /null)
      o3under300hpa_eosrl[valididx] = !values.f_nan
      o3under300hpa_nier[valididx] = !values.f_nan

     
      ;plot_sat_proj, o3under500hpa_eosrl - o3under500hpa_nier, longitude_eosrl, latitude_eosrl, $
        ;title='GEMS L2 O3P EOSRL NIER DIFF SFC-500hPa ' + datetime_str + 'UTC', $
        ;range=[-10., 10.], $
        ;pngfile=outputpath + 'gems_l2_o3p_' + datetime_str + '_under500hpa_eosrl_nier_diff.png', $
        ;cb_title='[DU]', $
        ;/diff

      ;plot_sat_proj, o3under300hpa_eosrl - o3under300hpa_nier, longitude_eosrl, latitude_eosrl, $
        ;title='GEMS L2 O3P EOSRL NIER DIFF SFC-300hPa ' + datetime_str + 'UTC', $
        ;range=[-10., 10.], $
        ;pngfile=outputpath + 'gems_l2_o3p_' + datetime_str + '_under300hpa_eosrl_nier_diff.png', $
        ;cb_title='[DU]', $
        ;/diff
        ;/scp_send, $
        ;scp_dest=scp_dest

      ;plot_sat_proj, o3under500hpa_poly - o3under500hpa_nier, longitude_poly,
        ;latitude_eosrl, $
        ;title='GEMS L2 O3P POLY NIER DIFF SFC-500hPa ' + datetime_str + 'UTC', $
        ;range=[-10., 10.], $
        ;pngfile=outputpath + 'gems_l2_o3p_' + datetime_str + '_under500hpa_poly_nier_diff.png';, $
        ;/scp_send, $
        ;scp_dest=scp_dest

      ;plot_sat_proj, o3under500hpa_eosrl, longitude_eosrl, latitude_eosrl, $
        ;title='GEMS L2 O3P EOSRL SFC-500hPa ' + sub[0] + ' ' + datetime_str + 'UTC', $
        ;range=[20, 60], $
        ;pngfile=outputpath + 'gems_l2_o3p_' + datetime_str + '_under500hpa_eosrl.png', $
        ;cb_title='[DU]'
        ;/scp_send, $
        ;scp_dest=scp_dest

      plot_sat_proj, o3under300hpa_eosrl, longitude_eosrl, latitude_eosrl, $
        title='GEMS L2 O3P EOSRL SFC-300hPa ' + sub[0] + ' ' + datetime_str + 'UTC', $
        range=[0, 60], $
        pngfile=outputpath + 'gems_l2_o3p_' + datetime_str + '_under300hpa_eosrl.png', $
        cb_title='[DU]'
        ;/scp_send, $
        ;scp_dest=scp_dest
      if i eq 0 then begin
        eosrl_o3_300hPa_00utc = o3under300hpa_eosrl
        lon_00utc = longitude_eosrl
        lat_00utc = latitude_eosrl
      endif

      if i eq 1 then begin

        sz03 = size(o3under300hpa_eosrl, /dimension)
        collocated_00_on_03 = fltarr(sz03[0], sz03[1])
        collocated_00_on_03 = !values.f_nan

        sz00 = size(eosrl_o3_300hpa_00utc, /dimension)


        for iy=0, sz00[1]-1 do begin
          for ix=0, sz00[0]-1 do begin
            x = lon_00utc[ix, iy]
            y = lat_00utc[ix, iy]
            result = search_closest_pixel(longitude_eosrl, latitude_eosrl, x, y)
            collocated_00_on_03[result] = eosrl_o3_300hPa_00utc[ix, iy]
          endfor
        endfor

        stop
      endif




      ;plot_sat_proj, o3under500hpa_poly, longitude_poly, latitude_poly, $
        ;title='GEMS L2 O3P BTDF SFC-500hPa ' + sub[1] + ' ' + datetime_str + 'UTC', $
        ;range=[20, 60], $
        ;pngfile=outputpath + 'gems_l2_o3p_' + datetime_str + '_under500hpa_btdf_poly.png';, $
        ;/scp_send, $
        ;scp_dest=scp_dest

      ;plot_sat_proj, o3under500hpa_nier, longitude_nier, latitude_nier, $
        ;title='GEMS L2 O3P NIER SFC-500hPa ' + ' ' + datetime_str + 'UTC', $
        ;range=[20, 60], $
        ;pngfile=outputpath + 'gems_l2_o3p_' + datetime_str + '_under500hpa_nier.png', $
        ;cb_title='[DU]'
        ;;/scp_send, $
        ;;scp_dest=scp_dest

      plot_sat_proj, o3under300hpa_nier, longitude_nier, latitude_nier, $
        title='GEMS L2 O3P NIER SFC-300hPa ' + ' ' + datetime_str + 'UTC', $
        range=[0, 60], $
        pngfile=outputpath + 'gems_l2_o3p_' + datetime_str + '_under300hpa_nier.png', $
        cb_title='[DU]'
        ;/scp_send, $
        ;scp_dest=scp_dest

      total_o3accum_eosrl_500 = [total_o3accum_eosrl_500, o3under500hpa_eosrl]
      ;total_o3accum_nier_500 = [total_o3accum_nier_500, o3under500hpa_nier]

      total_o3accum_eosrl_300 = [total_o3accum_eosrl_300, o3under300hpa_eosrl]
      ;total_o3accum_nier_300 = [total_o3accum_nier_300, o3under300hpa_nier]
      ;cao3_eosrl = o3p_eosrl.ColumnAmountO3
      ;cao3_poly = o3p_poly.ColumnAmountO3
      ;if mondate eq '0806' then begin
        ;rangetotal = [250, 310]
      ;endif else begin
        ;rangetotal = [250, 345]
      ;ENDELSE

      ;rangetotal = [250, 500]
      ;sz = size(cao3_eosrl, /dim)
      ;if sz[0] ne 3 then begin
        ;data = cao3_eosrl[*, *, 0]
      ;endif else begin
        ;data = cao3_eosrl[0, *, *]
      ;endelse

      ;plot_sat_proj, reform(data), longitude, latitude, $
        ;title='GEMS L2 O3P Total ozone for ' + datetime_str + 'UTC', $
        ;range=rangetotal, $
        ;pngfile=outputpath + 'gems_l2_o3p_date_' + datetime_str +'_toz_' + sub[isub] + '.png', $
        ;/scp_send, $
        ;scp_dest=scp_dest
      ;plot_sat_proj, reform(cao3[2, *, *]), longitude, latitude, $
        ;title='GEMS L2 O3P Tropospheric ozone for ' + mondate + 'T' + utchhmm + 'UTC', $
        ;range=[20, 100], $
        ;pngfile=outputpath + 'gems_l2_o3p_2020' + mondate + utchhmm + '_tropocolumn_wl310340_me0.5.png', $
        ;/scp_send
      ;plot_sat_proj, reform(cao3[1, *, *]), longitude, latitude, $
        ;title='GEMS L2 O3P Stratospheric ozone for ' + mondate + 'T' + utchhmm + 'UTC', $
        ;range=[100, 200], $
        ;pngfile=outputpath + 'gems_l2_o3p_2020' + mondate + utchhmm + '_stratocolumn_wl310340_me0.5.png', $
        ;/scp_send
    ;endfor
  endif else begin
    print, 'Cannot find any file.'
  endelse
endfor

;plot_gems_validation, total_o3accum_eosrl_500, total_o3accum_nier_500, $
  ;filename='./plot/eosrl_test/gems_l2_o3p_eosrl_test_with_omi_' + $
    ;'20210329_20210502.png', $
  ;xtitle='GEMS EOSRL Tropospheric O3 under 500hPa', $
  ;ytitle='GEMS NIER Tropospheric O3 under 500hPa', $
  ;cblim=[0, 10000], $
  ;range=[0, 60], $
  ;delta=2.0

;plot_gems_validation, total_o3accum_eosrl_300, total_o3accum_nier_300, $
  ;filename='./plot/eosrl_test/gems_l2_o3p_eosrl_test_with_omi_' + $
    ;'20210329_20210502.png', $
  ;xtitle='GEMS EOSRL Tropospheric O3 under 300hPa', $
  ;ytitle='GEMS NIER Tropospheric O3 under 300hPa', $
  ;cblim=[0, 10000], $
  ;range=[0, 60], $
  ;delta=2.0

end
