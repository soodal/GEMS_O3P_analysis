
fn = 'OMIO3PROF_onlyuv2_2040_120125-o05139_L0754-0959_X02-33.out'
path = '/data1/gems/o3p/works/2_O3PR/OZBOREAS-OMI/src/'

fnames = path + fn

omi_result = ds_read_omi_profoz_ascii_out(fnames)

outputpath = '/data1/gems/o3p/works/GEMS_O3P_analysis/plot/' + 'omi_uv2/'
plot_sat_proj, reform(omi_result.omicol[*, 2, 4]), $
  reform(omi_result.omilon[*, 0:3]), $
  reform(omi_result.omilat[*, 0:3]), $
  /cornerpixel, $
  title='OMI L2 O3P SFC-TROPOPAUSE', $
  range=[20, 60], $
  pngfile=outputpath + 'omi_uv2_tropospheric_o3.png', $
  cb_title='[DU]'

end
