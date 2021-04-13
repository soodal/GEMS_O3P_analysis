
pohang_lon = 129.37963
pohang_lat = 36.03259

;o3p_fn = file_basename(fl[i])
;o3p = ds_read_gems_l2_o3p(o3p_path + o3p_fn, varlist=varlist)


ncpath ='/data1/gems/o3p/works/GEMS_O3P_Yonsei/out/winlim_test/'

scp_dest = 'soodal@164.125.38.179:/home/soodal/works/plot/winlim_test'
for i = 310, 305, -1 do begin
  subpath =  string(i, format='(i03)') + '340_EOSRL/'
  filelist = file_search(ncpath + subpath + '*.nc4')
  ;stop
  for j=0, n_elements(filelist) - 1 do begin
    basename = file_basename(filelist[j])
    o3p = ds_read_gems_l2_o3p(ncpath + subpath + basename, $
      varlist=varlist)

    ;ds_plot_gems_l2_o3p_all, o3p, $
      ;/use_structure, $
      ;basename=basename, $
      ;scppath ='soodal@164.125.38.179:/home/soodal/works/plot/winlim_test/' $
        ;+ strmid(subpath, 0, strlen(subpath)-1)

    idx = collocate_gemsl2o3p_ozonesonde(o3p.longitude, o3p.latitude, $
      pohang_lon, pohang_lat)
    indices = array_indices(o3p.longitude, idx)

    plot_gems_ak, o3p, indices[0], indices[1], $
      /use_structure, $
      scp_dest='soodal@164.125.38.179:/home/soodal/works/plot/winlim_test/' $
        + strmid(subpath, 0, strlen(subpath)-1), $
      basename=basename

  endfor
endfor
end
