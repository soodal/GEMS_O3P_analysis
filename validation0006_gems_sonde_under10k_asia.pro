
;gems_pohang_idx = [66,71]

stn = ['pohang', 'stn014', "stn344"]
stnName = ['Pohang', 'Tateno(Tsukuba)', "King's Park(Hong Kong)"]

stns = [] 
months = []

;=============================================================================
; Ozonesonde
;=============================================================================

ecfs = []

for istn = 0, 2 do begin
  o3sonde_path = '/data1/gems/o3p/works/GEMS_O3P_analysis/data/ozonesonde/' $
    + stn[istn] + '/'
  sonde_o3 = []
  gems_o3 = []
  
  for imonth = 8, 10 do begin
    sonde_o3_monthly = []
    gems_o3_monthly = []
    smonth = string(imonth, format='(i02)')

    if stn[istn] eq 'pohang' then begin
      o3sonde_filelist = file_search(o3sonde_path + '2020' + smonth  + '*.txt')
    endif else if stn[istn] eq 'stn014' then begin
      o3sonde_filelist = file_search(o3sonde_path + '2020' + smonth  + '*.csv')
      print, o3sonde_filelist
    endif else if stn[istn] eq 'stn344' then begin
      o3sonde_filelist = file_search(o3sonde_path + '20' + smonth  + '*.csv')
      print, o3sonde_filelist
    endif
    nobs = n_elements(o3sonde_filelist)

    

    if nobs eq 1 then begin
      if strlen(o3sonde_filelist) eq 0 then begin
        nofile = 1
      endif else begin
        nofile = 0
      endelse
    endif else begin
      nofile = 0
    endelse


    if not nofile then begin
      for isonde=0, nobs-1 do begin
        if stn[istn] eq 'pohang' then begin 
          sdate = strmid(o3sonde_filelist[isonde], /reverse, 8, 2)
        endif else if stn[istn] eq 'stn014' then begin 
          sdate = strmid(o3sonde_filelist[isonde], /reverse, 24, 2)
        endif else if stn[istn] eq 'stn344' then begin 
          sdate = strmid(o3sonde_filelist[isonde], /reverse, 7, 2)
        endif

        ;print, o3sonde_filelist[isonde], ' date: ', sdate
        if stn[istn] eq 'pohang' then begin
          sonde = !null
          read_ozonesonde_pohang, o3sonde_filelist[isonde], sonde
          vmr = sonde.ozon_mpa/1000/100/sonde.pres_hpa
          ds_vmr2du_1d, vmr, sonde.pres_hpa, sonde.temp_degc + 273.15, sonde.geop_h_gpm, o3_du_layer
          height = sonde.geop_h_gpm
          oslat = sonde.launch_latitude
          oslon = sonde.launch_longitude
        endif else begin
          sonde = !null
          read_ozonesonde_woudc, o3sonde_filelist[isonde], sonde
          vmr = sonde.O3PartialPressure/1000/100/sonde.pressure
          ds_vmr2du_1d, vmr, sonde.pressure, sonde.temperature + 273.15, sonde.GPHeight, o3_du_layer
          height = sonde.GPHeight
          oslat = sonde.launch_latitude
          oslon = sonde.launch_longitude
        endelse

        nvertical = n_elements(vmr)
        o3under10k_sonde = 0.
        for ilayer=0, nvertical-1 do begin
          if height[ilayer] lt 10000 then begin
            o3under10k_sonde += o3_du_layer[ilayer]
            ;print, 'ilayer:', ilayer, ' height:', height[ilayer]
          endif
        endfor

;=============================================================================
; gems
;=============================================================================
        ;o3p_outfilelist = file_search('/data2/L2_GEMS/val_1008/GK2_GEMS_O3P_2020' + smonth + sdate + '_0345.nc4')
        ;o3p_outfile = '/data2/L2_GEMS/val_1008/GK2_GEMS_O3P_2020' + smonth + sdate + '_0345.nc4'
        ;print, 'date: ', smonth, sdate
        ;print, o3p_outfile
        year = 2020
        month = imonth

        day = fix(sdate)
        hour = 3
        minute = 45

        print, year, month, day, hour, minute
        o3p_outfile = ds_get_o3p_filename(year, month, day, hour, minute)
        if o3p_outfile eq !NULL then begin
          o3p_outfile = ''
          print, 'there is no gems file'
        endif

        if file_test(o3p_outfile) then begin
          varlist = ['EffectiveCloudFractionUV', 'ProcessingQualityFlags', $
                      'AveragingKernel', $
                      'O3' , $
                       'O3Apriori', 'O3AprioriError','ColumnAmountO3', $
                      'CloudPressure', $
                      'Latitude' ,'Longitude', $
                      'Time','Altitude' ,    $
                      'Pressure', 'TropopausePressure', $
                      'Wavelengths', $
                      'WavelengthsWholeRange']
          o3p = ds_read_gems_l2_o3p(o3p_outfile, varlist=varlist)
          gems_ak = o3p.AveragingKernel
          ds_gems_l2o3p_accum, o3p, o3under10k_gems, height=10.

          gemslat = o3p.latitude
          gemslon = o3p.longitude
          nanidx = where(gemslat lt 500, /null)
          gemslat[nanidx] = !values.f_nan
          nanidx = where(gemslon lt 500, /null)
          gemslon[nanidx] = !values.f_nan

          distance = (o3p.latitude - oslat)^2 + (o3p.longitude - oslon)^2
          nearest_idx = where(distance le 2 and distance eq min(distance, $
            /nan), /null)
          if n_elements(nearest_idx) ne 0 then begin
            nearest_indices = array_indices(o3p.latitude, nearest_idx)
            print, nearest_indices
            print, oslon, oslat

            ;print, 'nearest_idx', nearest_idx
            ;print, 'oslon, oslat', oslon, oslat

            sonde_o3 = [sonde_o3, o3under10k_sonde]
            sonde_o3_monthly = [sonde_o3_monthly, o3under10k_sonde]
            if nearest_idx eq -1 then begin
              gems_o3 = [gems_o3, !values.f_nan]
              gems_o3_monthly = [gems_o3_monthly, !values.f_nan]
              ecfs = [ecfs, !values.f_nan]
            endif else begin
              gems_o3 = [gems_o3, o3under10k_gems[nearest_idx]]
              gems_o3_monthly = [gems_o3_monthly, o3under10k_gems[nearest_idx]]
              ecfs = [ecfs, o3p.EffectiveCloudFractionUV[nearest_idx]]
            endelse
              ;mean(o3under10k_gems[gems_pohang_idx[0]-2:gems_pohang_idx[0]+2, $
                ;gems_pohang_idx[1]-2:gems_pohang_idx[1]+2], /nan)]
            stns = [stns, stn[istn]]
            months = [months, imonth]
          endif


          ;if n_elements(o3p.o3[nearest_indices[0], nearest_indices[1], *]) ne 1 then begin
            ;print, n_elements(o3p.o3[nearest_indices[0],nearest_indices[1], *])
            ;p1 = plot(vmr, height, color='k', $
              ;name='Ozonesonde', /buffer, $
              ;xrange=[0, 50], $
              ;yrange=[0, 10000])
            ;p2 = plot(o3p.o3[nearest_indices[0], nearest_indices[1], *], $
                ;o3p.altitude[nearest_indices[0], nearest_indices[1], *]*1000, $
              ;color='b', $
              ;name='Satellite profile', $
              ;/Buffer, /overplot)
            ;p3 = plot(o3p.o3apriori[nearest_indices[0], nearest_indices[1], *], $
                ;o3p.altitude[nearest_indices[0], nearest_indices[1], *]*1000, $
              ;color='b', $
              ;name='Satellite profile', $
              ;/Buffer, /overplot)
            ;pngfile = 'profile_' + stn[istn] + '_' +smonth+sdate+'0345.png' 

            ;p1.save, pngfile
            ;p1.close

            ;scp_dest = 'soodal@164.125.38.179:/home/soodal/works/plot/profile'

            ;; send image to pc
            ;spawn, 'scp -P18742 -p ' + pngfile + $
              ;' ' + scp_dest
          ;endif


        endif
      endfor
    endif
    if n_elements(gems_o3_monthly) ge 1 and n_elements(gems_o3_monthly) ge 1 then begin
      idx = where(finite(gems_o3_monthly) eq 1 and finite(sonde_o3_monthly) eq 1, /null)
      plot_gems_validation, gems_o3_monthly[idx], sonde_o3_monthly[idx], $
        filename='./plot/validation0006_gems_sonde_' + stn[istn] + $
          '_' + '2020' + smonth+ '.png', $
        title='GEMS O3P under 10 KM', $
        xtitle = 'GEMS O3P Tropospheric column ozone', $
        ytitle = stnName[istn] + ' Tropospheric Column Ozone', $
        ct = ctnum, $
        background = bgcolor, $
        range=range, delta=delta
    endif
  endfor
  idx = where(finite(gems_o3) eq 1 and finite(sonde_o3) eq 1, /null)
  plot_gems_validation, gems_o3[idx], sonde_o3[idx], $
    filename='./plot/validation0006_gems_sonde_' + stn[istn] + '.png', $
    title='GEMS O3P under 10 KM', $
    xtitle = 'GEMS O3P Tropospheric column ozone', $
    ytitle = stnName[istn] + ' Tropospheric Column Ozone', $
    ct = ctnum, $
    background = bgcolor, $
    range=range, delta=delta

endfor

idx = where(finite(gems_o3) eq 1 and finite(sonde_o3) eq 1 and $
  stns eq 'pohang' or  stns eq 'sn344')

print, correlate(gems_o3[idx], sonde_o3[idx])


stop

;gemsyyyymmdd = strmid(o3p_outfilelist, 16, 8, /reverse)
;gemshhmi = strmid(o3p_outfilelist, 7, 4, /reverse)
;gemsyyyymmdd_hhmi = strmid(o3p_outfilelist, 16, 13, /reverse)

for ifile=0, n_elements(o3p_outfilelist)-1 do begin

  omivals = []
  gemsvals = []
  collocate_gemsl2o3p_omo3pr, year, mon, day, hour, mi, gemsvals, omivals, $
    on_under_300=1

  if n_elements(gemsvals) ge 2 and n_elements(omivals) ge 2 then begin
    if on_dailyplot then begin
      plot_gems_validation, gemsvals, omivals, $
        filename='./plot/gems_l2_o3p_val_with_sonde_' + gemsyyyymmdd[ifile] + $
          'T' + gemshhmi[ifile]+'_under10km.png_ecf02', $
          range=[0, 100], $
          delta=2.0
    endif
    total_omivals = [total_omivals, omivals]
    total_gemsvals = [total_gemsvals, gemsvals]
  endif
endfor

; plot tropospheric column ozone
plot_gems_validation, total_gemsvals, total_omivals, $
  filename='./plot/gems_l2_o3p_val_with_omi_' + $
    strmid(gemsyyyymmdd[0], 0, 6) +'_under300hpa_ecf02.png', $
  title='GEMS O3P under 10 KM', $
  xtitle = 'GEMS O3P Tropospheric column ozone', $
  ytitle = stn[istn] + ' Tropospheric Column Ozone', $
  cblim=[0, 100], $
  range=[0, 100], $
  delta=2.0

end
