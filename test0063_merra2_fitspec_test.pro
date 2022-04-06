varlist = ['EffectiveCloudFractionUV', 'ProcessingQualityFlags', $
           'O3' , $
           ;'O3Apriori', 'O3AprioriError', $
           ;'CloudPressure', $
           ;'SimulatedRadiances', $
           ;'O3Apriori', 'O3AprioriError',$
           'ColumnAmountO3', $
           'Latitude' ,'Longitude', $
           'Time','Altitude' ,    $
           ;'Pressure', 'TropopausePressure', $
           ;'Wavelengths', $
           'WavelengthsWholeRange', $
           'Fitspec', $
           'SimulatedRadiances', $
           'div_rad_all', $
           'div_rad_sun', $
           'corrected_rad', $
           'corrected_irrad']

residualfn = '/data/private/soodal/merra2_residual/residuals/GK2_GEMS_L2_O3P_20210621_0345_winliminit300_prec000000_climML_b4x4_o20210917T205155Z.nc'

nsub = 198


merra2_gems = ds_read_gems_l2_o3p(residualfn, varlist=varlist)

merra2_corrected_rad_all = merra2_gems.corrected_rad_all[*, *, 0:nsub-1]
nanidx = where(merra2_corrected_rad_all lt -1.e29, /null)
merra2_corrected_rad_all[nanidx] = !values.f_nan

merra2_corrected_irrad_all = merra2_gems.corrected_irrad_all[*, *, 0:nsub-1]
nanidx = where(merra2_corrected_irrad_all lt -1.e29, /null)
merra2_corrected_irrad_all[nanidx] = !values.f_nan

merra2_div_sun_all = merra2_gems.div_sun_all
nanidx = where(merra2_div_sun_all lt -1.e29, /null)
merra2_div_sun_all[nanidx] = !values.f_nan

merra2_div_rad_all = merra2_gems.div_rad_all 
nanidx = where(merra2_div_rad_all lt -1.e29, /null)
merra2_div_rad_all[nanidx] = !values.f_nan

merra2_fitspec = merra2_gems.fitspec 
nanidx = where(merra2_fitspec lt -1.e29, /null)
merra2_fitspec[nanidx] = !values.f_nan


;-----------------------------------------------------------------------------

l1cfn = '/data2/L1C_GEMS/L1C/4x4/GK2_GEMS_L1C_20210621_0345_NOR_694_4x4.nc'
l1crad = ds_read_gems_l1c_4x4(l1cfn) 

gems_rad_spec = l1crad.image_pixel_values[*, *, 0:nsub-1]
gems_rad_norm = total(gems_rad_spec, 3)/nsub
div_rad = gems_rad_norm
sz = size(div_rad, /dimension)

subspec = gems_rad_spec
gems_rad_spec = subspec / gems_rad_norm
curr_rad_spec = gems_rad_spec
corrected_rad = curr_rad_spec
;---------

irrfn = '/data2/L1C_GEMS/L1C/4x4/GK2_GEMS_IRR_20210620_4x4.nc'
l1cirr = ds_read_gems_l1c_4x4(irrfn)

l1birr_irr_all = l1cirr.image_pixel_values[*, 0:nsub-1]
irr_norm = total(l1birr_irr_all, 2)/nsub

gems_irrad_spec = l1birr_irr_all
subspec = gems_irrad_spec
gems_irrad_norm = total(subspec, 2)  / nsub

div_sun = gems_irrad_norm
div_sun_2d = transpose(rebin(div_sun, sz[1], sz[0]))


gems_irrad_spec = subspec / gems_irrad_norm
curr_irrad_spec = gems_irrad_spec

curr_sol_spec = gems_irrad_spec
spline_sun = curr_sol_spec
database = spline_sun
sunspec_ss = database
corrected_irrad = sunspec_ss

;-----------------------------------------------------------------------------
diff_corrected_rad = (merra2_corrected_rad_all - corrected_rad)/corrected_rad
diff_corrected_irrad = (merra2_corrected_irrad_all - corrected_irrad)/corrected_irrad
diff_div_sun = (merra2_div_sun_all - div_sun_2d)/div_sun_2d
diff_div_rad = (merra2_div_rad_all - div_rad)/div_rad


gems_fitspec = corrected_rad / corrected_irrad * gems_rad_norm / irr_norm

diff_fitspec = (merra2_fitspec - gems_fitspec)/gems_fitspec
end

