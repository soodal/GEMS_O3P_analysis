import os
import glob

fl = os.listdir('/data2/L2_GEMS/val_1008/softcal/')
fl_nc4 = [file for file in fl if file.endswith('.nc4')]

fp = fl_nc4[0]
print(fp)
quit()

o3p = ds_read_gems_l2_o3p(fp)


# cldnanidx = where(o3p.EffectiveCloudFractionUV lt -1e29, /null)
# cldnanindices = array_indices(o3p.EffectiveCloudFractionUV, cldnanidx)

# cldnanidx = where(o3p.EffectiveCloudFractionUV lt 0 or o3p.EffectiveCloudFractionUV gt 0.15, /null)
# cldnanindices = array_indices(o3p.EffectiveCloudFractionUV, cldnanidx)
# cnisz = size(cldnanindices, /dimension)

# cloudyidx = where(o3p.EffectiveCloudFractionUV gt 0.15, /null)
# cloudyindices = array_indices(o3p.EffectiveCloudFractionUV, cloudyidx)


# clearidx = where(o3p.EffectiveCloudFractionUV le 0.15 and $
  # o3p.EffectiveCloudFractionUV ge 0, /null)
# clearindices = array_indices(o3p.EffectiveCloudFractionUV, clearidx)
# ;cisz = size(clearindices, /dimension)

# allresfit = o3p.all_wavelength_residualsoffit

# for i = 0, cnisz[1]-1 do begin
  # allresfit[*, cldnanindices[0, i], cldnanindices[1, i]] = !values.f_nan
# ENDFOR

# nanidx = where(allresfit lt -1e29, /null)

# allresfit[nanidx] = !values.f_nan


# residuals_mean = fltarr(204, 512)
# residuals_stddev = fltarr(204, 512)
# for iw=0, 203 do begin
# FOR iy=0, 511 DO BEGIN
  # if total(finite(allresfit[iw, *, iy])) gt 1 then begin
    # residuals_mean[iw, iy] = mean(allresfit[iw, *, iy], /nan)
    # residuals_stddev[iw, iy] = stddev(allresfit[iw, *, iy], /nan)
  # endif
# ENDFOR
# endfor

# c = contour(residuals_mean


# END
