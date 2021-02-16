

filelist = file_search('/data2/app/gemsl2_2020_0714/src/o3p/v1.b/out/winlim/*/*0345*.nc4')
filelist = filelist[sort(filelist)]

for i=0, n_elements(filelist)-1 do BEGIN
  plot_gems_satproj,filelist[i]
endfor

end
