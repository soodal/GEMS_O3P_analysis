
pngpath ='/data1/gems/o3p/ds/GEMS_O3P_Yonsei/out/winlim_test/'
for i = 310, 307, -1 do begin
  testpath =  string(i, format='(i03)') + '340_EOSRL/'
  filelist = file_search(pngpath + testpath + '*.nc4')
  for j=0, n_elements(filelist) - 1 do begin
    ds_plot_gems_l2_o3p_all, filelist[j], scppath ='winlim_test/' $
      + strmid(testpath, 0, strlen(testpath)-1)
  endfor
endfor
end
