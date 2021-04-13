; Plot Ozone Spatial distribution from GEMS L2 O3P output file
pro ds_plot_gems_l2_o3p_850hpa, o3_850hpa_number_density, o3p_lon, o3p_lat, title=title, $
  outfile = outfile, scppath=scppath

; Set keyword
if not keyword_set(outfile) then begin
  outfile = './plot/gems.png'
endif

if not keyword_set(scppath) then begin
  scppath = '/home/soodal/works/plot/'
endif
;scppath = '/home/soodal/works/plot/' + scppath


plot_sat_proj, o3_850hpa_number_density, o3p_lon, o3p_lat, $
  title=title, $
  range=[0, 10], $
  pngfile=outfile, $
  /scp_send, $
  scp_dest = scppath
;print, outfile
;print, scppath

end
