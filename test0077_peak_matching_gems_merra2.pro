
compile_opt idl2

plotpath = './plot/peak_fitting/'
savefn ='./peak_fitting/gems_merra2_20210621_0345.sav' 

drawplot = 1

if not file_test(savefn) then begin
	;=============================================================================
	; read MERRA2 residual
	;=============================================================================
	merra2gemso3pfn = "/data/private/soodal/merra2_residual/residuals/GK2_GEMS_L2_O3P_20210621_0345_wli300_prec000000_climML_b4x4_o20211002T123142Z_ecf0.nc"
	merra2gemsl2o3p = ds_read_gems_l2_o3p(merra2gemso3pfn)

	merra2_fitspec = merra2gemsl2o3p.fitspec

	;=============================================================================
	; read L2 CLD
	;=============================================================================
	gemscldfn = '/data2/L2_GEMS.nier/monthly_val/CLOUD/GK2_GEMS_L2_20210621_0345_CLOUD_FW_DPRO_BIN4x4.nc'
	gemsl2clddata = ds_read_gems_l2_o3p(gemscldfn)

	gemsl2cld_ecf = gemsl2clddata.EffectiveCloudFraction
	nanidx = where(gemsl2cld_ecf lt -1e29, /null)
	gemsl2cld_ecf[nanidx] = !values.f_nan

	;=============================================================================
	; read L2 O3P with MERRA radiance 
	;=============================================================================
	gemso3pfn = '/home/soodal/data/softcal_test/raw/model_300_340/GK2_GEMS_L2_O3P_20210621_0345_winliminit300_prec000000_climML_b4x4_o20211002T153843Z_clima_maxiter10_ecf0.nc'
	gemsl2o3p = ds_read_gems_l2_o3p(gemso3pfn)

	;=============================================================================
	; read L1C RAD, IRR
	;=============================================================================
	l1cradfn = '/data2/L1C_GEMS/L1C/4x4/GK2_GEMS_L1C_20210621_0345_NOR_694_4x4.nc'
	irrfn = '/data2/L1C_GEMS/L1C/4x4/GK2_GEMS_IRR_20210620_4x4.nc'

	gemsrad = ds_read_gems_l1c_4x4(l1cradfn)
	gemsirr = ds_read_gems_l1c_4x4(irrfn)

	rad = reform(gemsrad.image_pixel_values[*, *, 0:203])
	sz = size(rad, /dimension)

	gems_rad_norm = total(rad, 3)/204.
	for iw=0, 203 do begin 
		rad[0:sz[0]-1, 0:sz[1]-1, iw] = rad[0:sz[0]-1, 0:sz[1]-1, iw] / gems_rad_norm
	endfor

	irr = reform(gemsirr.image_pixel_values[*, 0:203])
	gems_irrad_norm = total(irr, 2)/204.
	for iw=0, 203 do begin 
		irr[0:sz[1]-1, iw] = irr[0:sz[1]-1, iw] / gems_irrad_norm
	endfor

	irr3d = fltarr(sz[0], sz[1], sz[2])

	for i=0, sz[0]-1 do begin
		irr3d[i, *, 0:203] = irr
	endfor

	gems_fitspec = fltarr(sz[0], sz[1], sz[2])
	for ix=0,sz[0]-1 do begin
		for iw=0, 203 do begin
			gems_fitspec[ix, *, iw] = rad[ix, *, iw]/irr3d[ix, *, iw]*gems_rad_norm[ix, *]/gems_irrad_norm[*]
		endfor
	endfor


	save, filename=savefn
endif else begin
	restore, savefn
endelse

kernel = [0.1, 0.2, 0.4, 0.2, 0.1]
kernel = [0.1, 0.8, 0.1]

initialize = 0

;=============================================================================
; pick a target pixel
;=============================================================================

;for iy = 450, 512, 100 do begin
for iy = 0, sz[1]-1 do begin
  ;for ix = 0, sz[0]-1 do begin 
  for ix = 130, 130 do begin 
    ;xidx = 130
    ;yidx = 418

    xidx = ix
    yidx = iy

    xidx_str = string(xidx, format='(i03)')
    yidx_str = string(yidx, format='(i03)')
    ;=============================================================================
    ; plot for a merra2
    ;=============================================================================
    merra2_wwr = merra2gemsl2o3p.WavelengthsWholeRange
    merra2_xvals = reform(merra2_wwr[xidx, yidx, *])
    merra2_yvals = reform(merra2_fitspec[xidx, yidx, *])

    ; nan value
    nanidx = where(merra2_yvals lt -1e29, /null)
    merra2_yvals[nanidx] = !values.f_nan


    ;merra2_yvals_5 = convol(merra2_yvals, kernel, $
    merra2_yvals_3 = convol(merra2_yvals, kernel, $
      /nan, $
      /normalize)

    ;if total(finite(merra2_xvals[0:196], /nan)) + total(finite(merra2_yvals[0:196], /nan)) ne 0 then begin
    ;if total(finite(merra2_xvals[0:196], /nan)) + total(finite(merra2_yvals_5[0:196], /nan)) ne 0 then begin
    if total(finite(merra2_xvals[0:196], /nan)) + total(finite(merra2_yvals_3[0:196], /nan)) ne 0 then begin
      continue
    endif

    ;p = plot(merra2_xvals, merra2_yvals,$
    ;p = plot(merra2_xvals, merra2_yvals_5,$
    if drawplot then begin
      p = plot(merra2_xvals, merra2_yvals_3,$
          'k', $
          title = 'Normalized Radiance with Local Max/Min', $
          xtitle = 'Wavelength[nm]', ytitle = 'Normalized Radiance',$
          ;yrange = [0, 0.2], $
          name = 'MERRA2 Radiance', $
          /buffer, $
          yminor = 3, xminor = 3)
    endif
        

    merra2_local_extrema_index = []
    ;find the local max index
    ;merra2_local_max_index = local_max_finder(merra2_xvals, merra2_yvals)
    ;merra2_local_max_index = local_max_finder(merra2_xvals, merra2_yvals_5)
    merra2_local_max_index = local_max_finder(merra2_xvals, merra2_yvals_3)

    merra2_local_extrema_index = [merra2_local_extrema_index, merra2_local_max_index]

    ;extract the x/y pairings of the local max/min
    merra2_x_maxima = merra2_xvals[merra2_local_max_index]
    ;merra2_y_maxima = merra2_yvals[merra2_local_max_index]
    ;merra2_y_maxima = merra2_yvals_5[merra2_local_max_index]
    merra2_y_maxima = merra2_yvals_3[merra2_local_max_index]

    ;overplot the max/min on the existing plot
    if drawplot then begin
      p2 = scatterplot(merra2_x_maxima,  merra2_y_maxima, /current, /overplot, $
          /buffer, $
          symbol = 'o', sym_color = 'r', sym_thick = 2)
    endif


    ;find the local min index with /minima keyword
    ;local_min_index = local_max_finder(merra2_xvals, merra2_yvals, /minima)
    ;local_min_index = local_max_finder(merra2_xvals, merra2_yvals_5, /minima)
    local_min_index = local_max_finder(merra2_xvals, merra2_yvals_3, /minima)

    merra2_local_extrema_index = [merra2_local_extrema_index, local_min_index]

    order = sort(merra2_local_extrema_index)
    merra2_local_extrema_index = merra2_local_extrema_index[order]

    merra2_huggins_extrema_idx = where(merra2_xvals[merra2_local_extrema_index] gt 313.5 and $
              merra2_xvals[merra2_local_extrema_index] lt 330, /null)

    ;extract the x/y pairings of the local max/min
    merra2_x_minima = merra2_xvals[local_min_index]
    ;merra2_y_minima = merra2_yvals[local_min_index]
    ;merra2_y_minima = merra2_yvals_5[local_min_index]
    merra2_y_minima = merra2_yvals_3[local_min_index]

    merra2_x_extrema_vals = merra2_xvals[merra2_local_extrema_index[merra2_huggins_extrema_idx]]
    ;merra2_y_extrema_vals = merra2_yvals[merra2_local_extrema_index[merra2_huggins_extrema_idx]]
    ;merra2_y_extrema_vals = merra2_yvals_5[merra2_local_extrema_index[merra2_huggins_extrema_idx]]
    merra2_y_extrema_vals = merra2_yvals_3[merra2_local_extrema_index[merra2_huggins_extrema_idx]]

    ;overplot the min on the existing plot
    if drawplot then begin
      p3 = scatterplot(merra2_x_minima, merra2_y_minima, /current, /overplot, $
          /buffer, $
          symbol = 'o', sym_color = 'b', sym_thick = 2)
    endif

    ;=============================================================================
    ; plot for a gems
    ;=============================================================================
    ;gems_wwr = merra2gemsl2o3p.WavelengthsWholeRange
    gems_xvals = reform(merra2_wwr[xidx, yidx, *])
    gems_yvals = reform(gems_fitspec[xidx, yidx, *])

    ; nan value
    nanidx = where(gems_yvals lt -1e29, /null)
    gems_yvals[nanidx] = !values.f_nan

    ;gems_yvals_5 = convol(gems_yvals, kernel, $
    gems_yvals_3 = convol(gems_yvals, kernel, $
      /nan, $
      /normalize)

    ;p_gems = plot(gems_xvals, gems_yvals, 'b-.', $
    ;p_gems = plot(gems_xvals, gems_yvals_5, 'b-.', $
    if drawplot then begin
      p_gems = plot(gems_xvals, gems_yvals_3, 'b-.', $
          title = 'Normalized Radiance with Local Max/Min', $
          xtitle = 'Wavelength[nm]', ytitle = 'Normalized Radiance',$
          ;yrange = [0, 0.2], $
          name='GEMS Radiance', $
          /buffer, $
          yminor = 3, xminor = 3, /overplot)
    endif
        

    gems_local_extrema_index = []
    ;find the local max index
    ;gems_local_max_index = local_max_finder(gems_xvals, gems_yvals)
    ;gems_local_max_index = local_max_finder(gems_xvals, gems_yvals_5)
    gems_local_max_index = local_max_finder(gems_xvals, gems_yvals_3)
    gems_local_extrema_index = [gems_local_extrema_index, gems_local_max_index]

    ;extract the x/y pairings of the local max/min
    gems_x_maxima = gems_xvals[gems_local_max_index]
    ;gems_y_maxima = gems_yvals[gems_local_max_index]
    ;gems_y_maxima = gems_yvals_5[gems_local_max_index]
    gems_y_maxima = gems_yvals_3[gems_local_max_index]

    ;overplot the max/min on the existing plot
    if drawplot then begin
      p_gems2 = scatterplot(gems_x_maxima, gems_y_maxima, /current, /overplot, $
          /buffer, $
          symbol = 'o', sym_color = 'm', sym_thick = 2)
    endif


    ;find the local min index with /minima keyword
    ;gems_local_min_index = local_max_finder(gems_xvals, gems_yvals, /minima)
    ;gems_local_min_index = local_max_finder(gems_xvals, gems_yvals_5, /minima)
    gems_local_min_index = local_max_finder(gems_xvals, gems_yvals_3, /minima)
    gems_local_extrema_index = [gems_local_extrema_index, gems_local_min_index]

    order = sort(gems_local_extrema_index)
    gems_local_extrema_index = gems_local_extrema_index[order]

    gems_huggins_extrema_idx = where(gems_xvals[gems_local_extrema_index] gt 313.5 and $
              gems_xvals[gems_local_extrema_index] lt 330, /null)

    ;extract the x/y pairings of the local max/min
    gems_x_minima = gems_xvals[gems_local_min_index]
    ;gems_y_minima = gems_yvals[gems_local_min_index]
    ;gems_y_minima = gems_yvals_5[gems_local_min_index]
    gems_y_minima = gems_yvals_3[gems_local_min_index]

    gems_x_extrema_vals = gems_xvals[gems_local_extrema_index[gems_huggins_extrema_idx]]
    ;gems_y_extrema_vals = gems_yvals[gems_local_extrema_index[gems_huggins_extrema_idx]]
    ;gems_y_extrema_vals = gems_yvals_5[gems_local_extrema_index[gems_huggins_extrema_idx]]
    gems_y_extrema_vals = gems_yvals_3[gems_local_extrema_index[gems_huggins_extrema_idx]]



    print, gems_local_extrema_index[gems_huggins_extrema_idx]

    print, gems_x_extrema_vals
    print, merra2_x_extrema_vals
    print, gems_x_extrema_vals - merra2_x_extrema_vals
    ;overplot the min on the existing plot
    if drawplot then begin
      p_gems3 = scatterplot(gems_x_minima, gems_y_minima, /current, /overplot, $
          /buffer, $
          symbol = 'o', sym_color = 'c', sym_thick = 2)
    endif



    ;idx = where(finite(gems_yvals) eq 1 and finite(merra2_yvals) eq 1, /null)
    ;idx = where(finite(gems_yvals) eq 1 and finite(merra2_yvals_5) eq 1, /null)
    ;idx = where(finite(gems_yvals_5) eq 1 and finite(merra2_yvals_5) eq 1, /null)
    idx = where(finite(gems_yvals_3) eq 1 and finite(merra2_yvals_3) eq 1, /null)
    ;basic_err = total(abs(gems_yvals[idx] - merra2_yvals[idx]))
    ;basic_err = total(abs(gems_yvals[idx] - merra2_yvals_5[idx]))
    ;basic_err = total(abs(gems_yvals_5[idx] - merra2_yvals_5[idx]))
    basic_err = total(abs(gems_yvals_3[idx] - merra2_yvals_3[idx]))
    ; basic_err = 1.5766745303389320


    ;------------------------------------------------------------------------------
    ; GEMS peak fitting to RTM
    ;------------------------------------------------------------------------------
    new_gems_xvals = dblarr(1033)
    new_gems_yvals = dblarr(1033)


    merra2_extrema_num = n_elements(merra2_x_extrema_vals)
    gems_extrema_num = n_elements(gems_x_extrema_vals)
 
    if initialize eq 0 then begin
      merra2_start_index_array = fltarr(sz[0], sz[1], merra2_extrema_num)
      initialize = 1
    endif

    if drawplot then begin 
      if merra2_extrema_num ne gems_extrema_num then begin
        p.save, plotpath+'gems_merra2_radiance_x' + xidx_str + '_y' + yidx_str + '.png'
        p.close
        continue
      endif
    endif

    for iex=0, merra2_extrema_num-1 do begin 
      if iex eq 0 then begin
        merra2_interval = merra2_x_extrema_vals[0] - merra2_xvals[0]
        ;print, 'merra2_interval:', merra2_interval
        merra2_start_index = 0
      endif else begin 
        merra2_interval = merra2_x_extrema_vals[iex] - merra2_x_extrema_vals[iex-1]
        ;print, 'merra2_interval:', merra2_interval
        merra2_start_index = merra2_local_extrema_index[merra2_huggins_extrema_idx[iex-1]]
      endelse
      merra2_final_index = merra2_local_extrema_index[merra2_huggins_extrema_idx[iex]]
      ;print, 'merra2_xvals[iex]:', merra2_x_extrema_vals[iex]

      merra2_start_index_array[xidx, yidx, iex] = merra2_start_index


      if iex eq 0 then begin
        gems_interval = gems_x_extrema_vals[0] - gems_xvals[0]
        ;print, 'gems_interval:', gems_interval
        gems_start_index = 0
      endif else begin
        gems_interval = gems_x_extrema_vals[iex] - gems_x_extrema_vals[iex-1]
        ;print, 'gems_interval:', gems_interval
        gems_start_index =	gems_local_extrema_index[gems_huggins_extrema_idx[iex-1]]
      endelse
      gems_final_index =	gems_local_extrema_index[gems_huggins_extrema_idx[iex]]
      ;print, 'gems_xvals[iex]:', gems_x_extrema_vals[iex]

      ;gems_num = gems_local_extrema_index[gems_huggins_extrema_idx[iex]] 

      ;extrema_diff = gems_x_extrema_vals[iex] - merra2_x_extrema_vals[iex]


      new_gems_xvals[gems_start_index:gems_final_index] = $
        (gems_xvals[gems_start_index:gems_final_index] - gems_xvals[gems_start_index]) * $
                    merra2_interval / gems_interval + $
                    merra2_xvals[merra2_start_index]

        
      ;print, gems_xvals[0:gems_local_extrema_index[gems_huggins_extrema_idx[0]]]
      ;print, new_gems_xvals[0:gems_final_index]

      new_gems_yvals[merra2_start_index:merra2_final_index] = $
        ;interpol(gems_yvals[gems_start_index:gems_final_index], $
        ;interpol(gems_yvals_5[gems_start_index:gems_final_index], $
        interpol(gems_yvals_3[gems_start_index:gems_final_index], $
          new_gems_xvals[gems_start_index:gems_final_index], $
          merra2_xvals[merra2_start_index:merra2_final_index])
      ;if iex eq 0 then begin
        ;p_gems_m = plot(merra2_xvals[0:merra2_final_index], new_gems_yvals[0:merra2_final_index], $
            ;/overplot)
      ;endif else begin
        ;p_gems_m = plot(merra2_xvals[merra2_start_index:merra2_final_index], $
          ;new_gems_yvals[merra2_start_index:merra2_final_index], $
            ;/overplot)
      ;endelse

    endfor

      merra2_interval = merra2_xvals[196] - merra2_x_extrema_vals[merra2_extrema_num-1]
      merra2_start_index = merra2_local_extrema_index[merra2_huggins_extrema_idx[merra2_extrema_num-1]]
      merra2_final_index = 196

      print, 'merra2_interval:', merra2_interval
      print, 'merra2_xvals[iex]:', merra2_x_extrema_vals[merra2_extrema_num-1]

      gems_interval = gems_xvals[196] - gems_x_extrema_vals[gems_extrema_num-1]
      gems_start_index =	gems_local_extrema_index[gems_huggins_extrema_idx[gems_extrema_num-1]]
      gems_final_index =	196
      print, 'gems_interval:', gems_interval
      print, 'gems_xvals[iex]:', gems_x_extrema_vals[gems_extrema_num-1]

      ;gems_num = gems_local_extrema_index[gems_huggins_extrema_idx[iex]] 

      ;extrema_diff = gems_x_extrema_vals[iex] - merra2_x_extrema_vals[iex]


      new_gems_xvals[gems_start_index:gems_final_index] = $
        (gems_xvals[gems_start_index:gems_final_index] - gems_xvals[gems_start_index]) * $
                    merra2_interval / gems_interval + $
                    merra2_xvals[merra2_start_index]

        
      ;print, gems_xvals[0:gems_local_extrema_index[gems_huggins_extrema_idx[0]]]
      ;print, new_gems_xvals[0:gems_final_index]

      new_gems_yvals[merra2_start_index:merra2_final_index] = $
        ;interpol(gems_yvals[gems_start_index:gems_final_index], $
        ;interpol(gems_yvals_5[gems_start_index:gems_final_index], $
        interpol(gems_yvals_3[gems_start_index:gems_final_index], $
          new_gems_xvals[gems_start_index:gems_final_index], $
          merra2_xvals[merra2_start_index:merra2_final_index])

    if drawplot then begin
      p_gems_m = plot(merra2_xvals[0:196], $
        new_gems_yvals[0:196], $
        'r--', $
          /overplot, $
          name='Adjusted GEMS Radiance')

      leg =legend(target=[p, p_gems, p_gems_m], /normal, position=[0.94, 0.3])

      p.save, plotpath+'gems_merra2_radiance_x' + xidx_str + '_y' + yidx_str + '_gems_merra2smoothed3.png'
      p.close

      ;p_diff = plot(merra2_xvals[0:196], gems_yvals[0:196] - merra2_yvals[0:196], $
      ;p_diff = plot(merra2_xvals[0:196], gems_yvals[0:196] - merra2_yvals_5[0:196], $
      ;p_diff = plot(merra2_xvals[0:196], gems_yvals_5[0:196] - merra2_yvals_5[0:196], $
      p_diff = plot(merra2_xvals[0:196], gems_yvals_3[0:196] - merra2_yvals_3[0:196], $
        'k', $
        title='GEMS - MERRA2 Normalized Radiance', $
        name='GEMS MERRA2 Difference', $
        /buffer)
      ;p_diff_modified = plot(merra2_xvals[0:196], new_gems_yvals[0:196] - merra2_yvals[0:196], $
      p_diff_modified = plot(merra2_xvals[0:196], new_gems_yvals[0:196] - merra2_yvals_3[0:196], $
        'b', $
        name='GEMS(Modified) MERRA2 Difference', $
        /overplot, $
        /buffer)
      p_diff.yrange=[-0.05, 0.05]

      leg = legend(target=[p_diff, p_diff_modified], /normal, position=[0.94, 0.9])
      p_diff.save, plotpath+'gems_merra2_radiance_difference_modified_x' + xidx_str + '_y' + yidx_str + '_gems_merra2smoothed3.png'
      p_diff.close
    endif

  endfor
endfor
stop


; filtering the pixels

ecf = gemsl2cld.EffectiveCloudFraction
nanidx = where(ecf lt -1e29, /null)
ecf[nanidx] = !values.f_nan

nanidx = where(ecf lt 0.2, /null)


gems_to3[nanidx] = !values.f_nan
nanidx = where(ecf gt 0.2, /null)
gems_to3[nanidx] = !values.f_nan



end

