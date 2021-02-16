; plot for 305340, 310340, kari and eosrl
path = '/data1/gems/o3p/works/GEMS_O3P_Yonsei/out/'
maker = ['KARI', 'EOSRL']


for iws = 305, 310, 5 do begin
  for imaker=0, 1 do begin

    subpath = string(iws, format='(i03)') + '340_'+maker[imaker] + '/'
    fl = file_search(path+'/'+subpath+'/'+'*.nc4')

    pngfile = string(iws, format='(i03)') + '340_' +maker[imaker] + $
      '_2020-06-16_0345UTC.png'
    plot_gems_l2_o3p_satproj, fl[0], height=10, pngfile=pngfile
  endfor
endfor


end
