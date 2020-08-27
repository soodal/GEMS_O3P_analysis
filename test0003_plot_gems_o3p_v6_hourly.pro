filelist = file_search('/data1/L2_GEMS/group/o3p/4x4/nopolc/v6/*')
for i =0, n_elements(filelist) - 1 do begin
  ds_plot_gems_l2_o3p_all, filelist[i], scppath ='hourly/'
endfor
end
