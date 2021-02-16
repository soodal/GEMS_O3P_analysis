
omifn = '/data2/OMI/gdata/2020/L2OMO3PR/OMI-Aura_L2-OMO3PR_2020m0801t0009-o85347_v003-2020m0802t082803.SUB.he5'

filelist = file_search('/data2/L2_GEMS/val_1008/GK2_GEMS_O3P_*.nc4')

for i = 0, n_elements(filelist)-1 do begin
  plot_gems_satproj,filelist[i]
endfor
end
