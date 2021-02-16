
imonth = 8

smonth = string(imonth, format='(i02)')

sonde_o3 = []
gems_o3 = []
;=============================================================================
; gems
;=============================================================================
o3p_outfilelist = file_search('/data2/L2_GEMS/val_1008/GK2_GEMS_O3P_2020' + smonth + '??_0345.nc4')

gems_pohang_idx = [66,71]
pohang_latlon = [36.03259, 129.37963]
nearest_idx = [67, 70]


varlist = ['EffectiveCloudFractionUV', 'ProcessingQualityFlags', $
            'AveragingKernel', $
            'O3' , $
            'O3Apriori', 'O3AprioriError', $
            'CloudPressure', $
            'SimulatedRadiances', 'Latitude' ,'Longitude', $
            'Time','Altitude' ,    $
            'Pressure', 'TropopausePressure', $
            'Wavelengths', $
            'WavelengthsWholeRange']

o3p = ds_read_gems_l2_o3p(o3p_outfilelist[0], varlist=varlist)

gems_ak = o3p.AveragingKernel

ds_gems_l2o3p_accum, o3p, o3under10k_gems, height=10.

gems_o3 = [gems_o3, o3under10k_gems[gems_pohang_idx[0], gems_pohang_idx[1]]]
;=============================================================================
; Ozonesonde
;=============================================================================

o3sonde_path = '/data1/gems/o3p/works/GEMS_O3P_analysis/data/ozonesonde/pohang/'
o3sonde_filelist = file_search(o3sonde_path + '2020' + smonth  + '*_L0.txt')

read_ozonesonde_pohang, o3sonde_filelist[0], sonde

vmr = sonde.ozon_mpa/1000/100/sonde.pres_hpa

ds_vmr2du_1d, vmr, sonde.pres_hpa, sonde.temp_degc + 273.15, sonde.geop_h_gpm, o3_du_layer

nsonde = n_elements(vmr)
o3under10k_sonde = 0.
for i=0, nsonde-1 do begin
  if sonde.geop_h_gpm[i] lt 10000 then begin
    o3under10k_sonde += o3_du_layer[i]
  endif
endfor
 
sonde_o3 = [sonde_o3, o3under10k_sonde]

stop



gemsyyyymmdd = strmid(o3p_outfilelist, 16, 8, /reverse)
gemshhmi = strmid(o3p_outfilelist, 7, 4, /reverse)
gemsyyyymmdd_hhmi = strmid(o3p_outfilelist, 16, 13, /reverse)

total_omivals = []
total_gemsvals = []
on_dailyplot = 0
for ifile=0, n_elements(o3p_outfilelist)-1 do begin
  year = fix(strmid(gemsyyyymmdd[ifile], 0, 4))
  mon = fix(strmid(gemsyyyymmdd[ifile], 4, 2))
  day = fix(strmid(gemsyyyymmdd[ifile], 6, 2))
  hour = fix(strmid(gemshhmi[ifile], 0, 2))
  mi = fix(strmid(gemshhmi[ifile], 2, 2))

  omivals = []
  gemsvals = []
  collocate_gemsl2o3p_omo3pr, year, mon, day, hour, mi, gemsvals, omivals, $
    on_under_300=1

  if n_elements(gemsvals) ge 2 and n_elements(omivals) ge 2 then begin
    if on_dailyplot then begin
      plot_gems_validation, gemsvals, omivals, $
        filename='./plot/gems_l2_o3p_val_with_omi_' + gemsyyyymmdd[ifile] + $
          'T' + gemshhmi[ifile]+'_under300hpa.png_ecf02', $
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
  cblim=[0, 100], $
  range=[0, 100], $
  delta=2.0

end
