filename = '/data1/gems/o3p/works/GEMS_O3P_analysis/data/ozonesonde/stn014/20200916.ECC.1Z.1Z35891.JMA.csv'
;read_ozonesonde_woudc, filename, result
filename = '/data1/gems/o3p/works/GEMS_O3P_analysis/data/ozonesonde/stn344/20080705.csv'
read_ozonesonde_woudc, filename, result
print, n_elements(result.temperature), n_elements(result.gpheight)
end
