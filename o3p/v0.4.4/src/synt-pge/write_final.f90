  ! ***************** Modification History ******************
  ! xiong liu, July 2003
  ! 1. Add an option for writing ozone retrieval results in 
  !    SUBROUTINE write_intermed
  ! *********************************************************

 
SUBROUTINE write_final ( fitted_col, rmsavg, davg, drelavg, n_spec )

  ! **********************************
  !
  !   Final WRITE of fitting results
  !
  ! **********************************

  USE OMSAO_precision_module
  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  INTEGER,        INTENT (IN) :: n_spec
  REAL (KIND=dp), INTENT (IN) :: fitted_col, rmsavg, davg, drelavg

  ! Write out the average fitting statistics
  WRITE (*,*)
  WRITE (*,'(A, 1PE13.5)') '                      Avg Col = ', fitted_col/n_spec
  WRITE (*,'(A, 1PE13.5)') '                      Avg RMS = ', rmsavg/n_spec
  WRITE (*,'(A, 1PE13.5)') '                     Avg dCol = ', davg/n_spec
  WRITE (*,'(A, 1PE13.5)') ' Avg relative Col uncertainty = ', drelavg/n_spec

  RETURN
END SUBROUTINE write_final

SUBROUTINE omi_write_intermed (founit, fitcol, dfitcol, rms, exval )
  
  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: maxloc, maxwin,max_fit_pts
  USE OMSAO_indices_module,    ONLY: ring_idx, shi_idx, squ_idx,&
       wvl_idx, spc_idx, sig_idx, hwe_idx, hwr_idx, vgr_idx, vgl_idx,  &
       asy_idx, hwl_idx
  USE ozprof_data_module,      ONLY: ozprof, ozprof_std, ozprof_init, &
       ozprof_ap, ozprof_apstd, eff_alb, eff_alb_init, nlay, ozdfs, ozinfo, num_iter, &
       covar, contri, avg_kernel, use_lograd, use_oe, nalb, atmosprof, ntp, nlay_fit, &
       ozprof_start_index, ozprof_end_index, ozfit_start_index, ozfit_end_index, &
       start_layer, end_layer, the_ctp, the_cfrac, the_cod, lambcld_refl, do_lambcld, &
       which_cld, which_alb, ozprof_nstd, strataod, stratsca, tropaod, tropsca, &
       aerwavs, maxawin, actawin, ozwrtcorr, ozwrtcovar, ozwrtcontri, the_ai, the_cld_flg,  &
       ozwrtres, ozwrtavgk, ozwrtvar, fgasidxs, ngas, nfgas, gaswrt, tracegas, saa_flag, &
       nsaa_spike, ozwrtfavgk, favg_kernel, radcalwrt, nsfc, ozwrtwf, ozwrtsnr, &
       weight_function, do_simu, glintprob, trace_prof, trace_avgk, trace_contri, which_clima,use_flns
  USE OMSAO_variables_module, ONLY : the_sza_atm, the_vza_atm, the_aza_atm, the_lons, the_lats, nloc, &
       the_utc, fitvar_rad, mask_fitvar_rad, n_fitvar_rad, fitvar_rad_std, n_rad_wvl, fitvar_rad_apriori, &
       fitspec_rad, fitres_rad, fitwavs, fitvar_rad_str, fitvar_rad_nstd, currline, currpix, currloop, &
       simspec_rad, clmspec_rad, actspec_rad, fitweights, database, refidx, numwin, nradpix, curr_x, curr_y,&
       wavcal, yn_varyslit
  USE SYNT_data_module, ONLY : land_water_flg, glint_flg, snow_ice_flg 


  IMPLICIT NONE

  ! ===============
  ! Input Variables
  ! ===============
  INTEGER,        INTENT (IN) :: founit, exval
  REAL (KIND=dp), INTENT (IN) :: rms
  REAL (KIND=dp), DIMENSION(3), INTENT (IN)    :: fitcol
  REAL (KIND=dp), DIMENSION(3, 2), INTENT (IN) :: dfitcol


  ! ===============
  ! Local variables
  ! ===============
  INTEGER                                  :: i, j, id, fidx, lidx
  REAL (KIND=dp), DIMENSION (2*maxloc)     :: latlon
  REAL (KIND=dp)                           :: avgres
  REAL (KIND=dp), DIMENSION (n_fitvar_rad) :: correl
  REAL (KIND=dp), DIMENSION (maxwin)       :: allrms, allavgres
  CHARACTER(1) :: ff



  !REAL (KIND=dp), DIMENSION (n_rad_wvl)    :: simrad
  
  WRITE(founit, '(A11,2I5,1X,A28)') 'Line/XPix: ',curr_y, curr_x, the_utc

  DO i = 1, nloc
     latlon(2*i-1) = the_lats(i); latlon(2*i) = the_lons(i)
  END DO
  WRITE(founit,'(13F8.2)') latlon(1:2*nloc), the_sza_atm,the_vza_atm,the_aza_atm
  WRITE(founit, '(A)') &
       'Exit Status, # Iterations, SAA, # SAA Spike, Land/Water, Glint, Snow/Ice Glint Probability'
  WRITE(founit, '(2I4, L4, 4I4, F5.2)') exval, num_iter, saa_flag, nsaa_spike,   &
       land_water_flg(currpix, currloop), glint_flg(currpix, currloop), &
       snow_ice_flg(currpix, currloop), glintprob

  IF (exval >= 0) THEN

     fidx = 1
     DO i = 1, numwin
        lidx = fidx + nradpix(i) - 1
        allrms(i) = SQRT(SUM((fitres_rad(fidx:lidx)/fitweights(fidx:lidx))**2) / nradpix(i))
        fidx = lidx + 1
     ENDDO
     simspec_rad(1:n_rad_wvl) = fitspec_rad(1:n_rad_wvl) - fitres_rad(1:n_rad_wvl)
     IF (use_lograd) THEN
        fitspec_rad(1:n_rad_wvl) = EXP(fitspec_rad(1:n_rad_wvl))
        simspec_rad(1:n_rad_wvl) = EXP(simspec_rad(1:n_rad_wvl))
        fitres_rad(1:n_rad_wvl) = fitspec_rad(1:n_rad_wvl) -  simspec_rad(1:n_rad_wvl)
     END IF
     avgres = SQRT(SUM((ABS(fitres_rad(1:n_rad_wvl)) / &
          fitspec_rad(1:n_rad_wvl))**2.0)/n_rad_wvl)*100.0   

     fidx = 1
     DO i = 1, numwin
        lidx = fidx + nradpix(i) - 1
        allavgres(i) = SQRT(SUM((ABS(fitres_rad(fidx:lidx)) / &
             fitspec_rad(fidx:lidx))**2.0)/nradpix(i))*100.0
        fidx = lidx + 1
     ENDDO
 
     WRITE(founit, '(A)')    'rms, avgres, dfs, info'
     WRITE(founit, '(14f8.3)') rms, avgres, ozdfs, ozinfo, allrms(1:numwin), allavgres(1:numwin)

     WRITE(founit, '(A8, I3, 1P28E10.2)') 'Albedo: ', nalb, eff_alb_init(1:nalb), eff_alb(1:nalb)
     !print * , 'alb0',eff_alb_init(1:2)
     !print * , 'alb ', eff_alb(1:2)
     WRITE(founit, '(A8,3F10.3, I3, L3, F10.3)') 'Cloud :  ', the_cfrac, the_cod, the_ctp, &
          the_cld_flg, do_lambcld, lambcld_refl
     WRITE(founit, '(A8, I3, 1P36E10.2)') 'Aerosol:  ', actawin, aerwavs(1:actawin),       &
          tropaod(1:actawin), tropsca(1:actawin), strataod(1:actawin), stratsca(1:actawin), the_ai
     
     WRITE(founit, '(A31, 3I5)') 'Atmosphere and ozone profiles: ', nlay, ntp, nsfc
     WRITE(founit, '(A)') '  #    P(mb)   Z(km)  T(K)     O3-ap O3-apstd   O3     STD     NSTD (DU)'
     WRITE(founit, '(A2, 1x, F9.3, 2F8.3)') ' 0', atmosprof(1:3, 0)
     DO i = 1, nlay 
        WRITE(founit, '(I2, 1X, F9.3, 7F8.3)') i, atmosprof(1:3, i), ozprof_ap(i), &
             ozprof_apstd(i), ozprof(i), ozprof_std(i), ozprof_nstd(i)
     END DO
     WRITE(founit, '(A24,4x,F8.3,8X,3F8.3)') ' Total Ozone: ', &
          SUM(ozprof_ap(1:nsfc)), fitcol(1), dfitcol(1, 1), dfitcol(1, 2)
     IF (ntp > 0) THEN
        WRITE(founit, '(A24,4x,F8.3,8X,3F8.3)') ' Stratospheric Ozone: ',&
             SUM(ozprof_ap(1:ntp)), fitcol(2), dfitcol(2, 1), dfitcol(2, 2)
        WRITE(founit, '(A24,4x,F8.3,8X, 3F8.3)') ' Tropospheric Ozone: ',&
             SUM(ozprof_ap(ntp+1:nsfc)), fitcol(3), dfitcol(3, 1), dfitcol(3, 2)          
     ENDIF

     IF (gaswrt) THEN
        WRITE(founit, '(A, I5)') 'Fitted trace gases and uncertainty: ', nfgas
        WRITE(founit, '(A)')  &
             '    Var    Initial   A Priori A Priori Std  VCD      STD       NSTD      AMF      ACFRAC    AVGK(1)   AVGK(2)'
        DO i = 1, ngas
           IF (fgasidxs(i) > 0) THEN
              j = mask_fitvar_rad(fgasidxs(i))
              WRITE(founit, '(I2,1x,A6,1P10d10.2)') fgasidxs(i), fitvar_rad_str(j), tracegas(i, 1:10)
           ENDIF
        ENDDO
        
        ! xliu, 08/11/2010, Add averaging kernels/ a priori profile (shape)
        IF ( ozwrtavgk ) THEN
           WRITE(founit, '(A, 2I5)') 'Trace gas averaging kernels and a priori: ', nfgas, nlay
           DO i = 1, ngas
              IF (fgasidxs(i) > 0) THEN
                 WRITE(founit, '(1P100d12.4)') trace_prof(i, 1:nlay)
                 WRITE(founit, '(100F9.4)') trace_avgk(i, 1:nlay)
              ENDIF
           ENDDO
        ENDIF

        IF ( ozwrtcontri) THEN
           WRITE(founit, '(A, 2I5)') 'Trace gas contribution function: ', nfgas, n_rad_wvl
           DO i = 1, ngas
              IF (fgasidxs(i) > 0) THEN
                 WRITE(founit, '(1P1000d12.4)') trace_contri(i, 1:n_rad_wvl)
              ENDIF
           ENDDO
        ENDIF                  
     ENDIF
     
     IF (ozwrtvar) THEN                      
        WRITE(founit, '(A, I5)') 'Fitted variables and uncertainty: ', n_fitvar_rad
        DO i = 1, n_fitvar_rad
           WRITE(founit, '(I2, 1x,A6,1P3d10.2)') i, fitvar_rad_str(mask_fitvar_rad(i)), &
                fitvar_rad(mask_fitvar_rad(i)), fitvar_rad_std(mask_fitvar_rad(i)), &
                fitvar_rad_nstd(mask_fitvar_rad(i))
        END DO
     ENDIF
     
     IF (ozwrtcorr) THEN
        WRITE(founit, '(A, I5)') 'Correlation matrix: ', n_fitvar_rad
        DO i = 1, n_fitvar_rad
           correl(i) = 1.0
           DO j = 1, i - 1
              correl(j) = covar(i, j) / SQRT(covar(i, i) * covar(j, j))
           END DO
           WRITE(founit, '(I2,1X,100f6.2)') i, correl(1:i)
        END DO
     ENDIF
        
     IF (ozwrtcovar) THEN
        WRITE(founit, '(A, 3I5)') 'Covariance matrix: ', nlay_fit, start_layer, end_layer
        DO i = ozfit_start_index, ozfit_end_index
           WRITE(founit, '(1p60d10.2)') covar(i, ozfit_start_index:ozfit_end_index)
        END DO
     ENDIF
     
     IF (ozwrtavgk) THEN
        WRITE(founit, '(A, 3I5)') 'Average kernel: ', nlay_fit, start_layer, end_layer
        DO i = ozfit_start_index, ozfit_end_index
           WRITE(founit, '(1p60d10.2)') avg_kernel(i, ozfit_start_index:ozfit_end_index)
        END DO
     ENDIF

     IF (ozwrtfavgk) THEN
        WRITE(founit, '(A, 3I5)') 'Average kernel1: ', nlay_fit, start_layer, end_layer
        DO i = ozfit_start_index, ozfit_end_index
           WRITE(founit, '(1p60d10.2)') favg_kernel(i, ozfit_start_index:ozfit_end_index)
        END DO
     ENDIF
     
     IF (ozwrtcontri) THEN
        WRITE(founit, '(A, I5)') 'Contribution Function: ', n_rad_wvl
        DO i = 1, n_rad_wvl
           WRITE(founit, '(1p60d10.2)') contri(ozfit_start_index:ozfit_end_index, i)
        ENDDO
     ENDIF

     IF (ozwrtwf) THEN
        WRITE(founit, '(A, 2I5)') 'Weighting Function: ', n_rad_wvl, n_fitvar_rad
        DO i = 1, n_rad_wvl
           WRITE(founit, '(1p100d11.3)') weight_function(i, 1:n_fitvar_Rad)
        ENDDO
     ENDIF

     IF (ozwrtsnr) THEN
        WRITE(founit, '(A, I5)') 'Measurement Error: ', n_rad_wvl
        WRITE(founit, '(1p100d11.3)') fitweights(1:n_rad_wvl)
        !WRITE(*, '(1p100d11.3)') fitweights(1:n_rad_wvl)
     ENDIF
     
     IF (ozwrtres) THEN
        WRITE(founit, '(A, I5)') 'Fit residual: ', n_rad_wvl
        DO i = 1, n_rad_wvl
           WRITE(founit, '(f9.4, 1p3d12.4)') & 
                    fitwavs(i), fitspec_rad(i), simspec_rad(i),actspec_rad(i)
        END DO

        !   write(*,'(10f8.2)') fitwavs(1:3), fitwavs(n_rad_wvl-3:n_rad_wvl)
         !    print * , '-------------'
         !  write(*,'(10f8.3)')  (fitspec_rad(1:3) -simspec_rad(1:3))*100/simspec_rad(1:3) , &
         !   (fitspec_rad(n_rad_wvl-3:n_rad_wvl) -simspec_rad(n_rad_wvl-3:n_rad_wvl))*100/simspec_rad(n_rad_wvl-3:n_rad_wvl)
     ENDIF

     !IF (ozwrtres) THEN
     !   WRITE(founit, '(A, I5)') 'Ring+Fit residual: ', n_rad_wvl
     !   DO i = 1, n_rad_wvl
     !      WRITE(founit, '(f9.4, 1p3d12.4)') fitwavs(i), fitspec_rad(i), &
     !      simspec_rad(i), database(ring_idx, refidx(i))
     !   END DO
     !ENDIF

     IF (radcalwrt .AND. .NOT. do_simu) THEN
        WRITE(founit, '(A, I5)') 'Radiance Calibration: ', n_rad_wvl
        DO i = 1, n_rad_wvl
           WRITE(founit, '(f9.4, 1p4d12.4)') fitwavs(i), fitspec_rad(i), simspec_rad(i), &
                clmspec_rad(i), actspec_rad(i)          
        END DO
     ENDIF
  END IF
 

 RETURN
END SUBROUTINE omi_write_intermed





!SUBROUTINE gome_write_intermed (founit, npix, nsub, fitcol, dfitcol, rms, amf, amfgeo, &
     !sol_zen_eff, nang, sza_atm, vza_atm, ngeo, lat, lon, exval )
  
  !USE OMSAO_precision_module
  !USE OMSAO_indices_module,    ONLY: no2_t1_idx, so2_idx, bro_idx, hcho_idx
  !USE ozprof_data_module,      ONLY: ozprof_flag, ozprof, ozprof_std, ozprof_init, &
       !ozprof_ap, ozprof_apstd, eff_alb, eff_alb_init, nlay, ozdfs, ozinfo, num_iter, &
       !covar, contri, avg_kernel, use_lograd, use_oe, nalb, atmosprof, ntp, nlay_fit, &
       !ozprof_start_index, ozprof_end_index, ozfit_start_index, ozfit_end_index, &
       !start_layer, end_layer, the_ctp, the_cfrac, the_cod, the_orig_cfr, the_orig_ctp, &
       !the_cld_flg, the_orig_cod, the_ai, lambcld_refl, do_lambcld, &
       !which_cld, which_alb, ozprof_nstd, strataod, stratsca, tropaod, tropsca, &
       !aerwavs, maxawin, actawin, ozwrtcorr, ozwrtcovar, ozwrtcontri,  &
       !ozwrtres, ozwrtavgk, ozwrtvar, fgasidxs, ngas, nfgas, gaswrt, tracegas, &
       !saa_flag, nsaa_spike, ozwrtfavgk, favg_kernel, nsfc
  !USE OMSAO_variables_module,  ONLY : the_sza_atm, the_vza_atm, the_aza_atm, chisq, &
       !fitvar_rad, mask_fitvar_rad, n_fitvar_rad, fitvar_rad_std, n_rad_wvl, fitvar_rad_apriori, &
       !fitspec_rad, fitres_rad, fitwavs, fitvar_rad_str, fitvar_rad_nstd
  !USE OMSAO_gome_data_module,  ONLY : gome_orbc, gome_pixdate, gome_angles_wrtn, &
       !orbnum, gome_stpix, gome_endpix, gome_npix
  
  !IMPLICIT NONE

  !! ===============
  !! Input Variables
  !! ===============
  !INTEGER,        INTENT (IN) :: founit, npix, nsub, nang, ngeo, exval
  !REAL (KIND=dp), INTENT (IN) :: rms, amf, amfgeo, sol_zen_eff
  !REAL (KIND=dp), DIMENSION(3), INTENT (IN)    :: fitcol
  !REAL (KIND=dp), DIMENSION(3, 2), INTENT (IN) :: dfitcol
  !REAL (KIND=dp), DIMENSION (nang),INTENT (IN) :: sza_atm, vza_atm
  !REAL (KIND=dp), DIMENSION (ngeo),INTENT (IN) :: lat, lon

  !! ===============
  !! Local variables
  !! ===============
  !INTEGER                                  :: i, j, id
  !REAL (KIND=dp), DIMENSION (2*ngeo)       :: latlon
  !REAL (KIND=dp)                           :: avgres
  !REAL (KIND=dp), DIMENSION (n_fitvar_rad) :: correl
  !REAL (KIND=dp), DIMENSION (n_rad_wvl)    :: simrad

  !DO i = 1, ngeo
     !latlon(2*i-1) = lat(i); latlon(2*i) = lon(i)
  !END DO

  !IF (.NOT. ozprof_flag) THEN
     !WRITE (founit, '(I5,I3,1P8E12.4,I8, 2X, 0P16F7.2)') &
          !npix, nsub, fitcol, dfitcol, fitcol/amf, dfitcol/amf, rms, amf, amfgeo, &
          !sol_zen_eff, exval, sza_atm(1:nang), vza_atm(1:nang), latlon(1:2*ngeo)
  !ELSE     
     !WRITE(founit, '(A12,I5,A5,2I5,1x,A24, 2I5)') 'GOME Pixel: ', orbnum, &
          !gome_orbc, npix, nsub, gome_pixdate, gome_stpix, gome_npix
     !WRITE(founit,'(3F7.2)') the_sza_atm,the_vza_atm,the_aza_atm
     !WRITE(founit,'(10F7.2)') latlon(1:2*ngeo)
     !WRITE(founit, '(A13,I5, A15, I5, A6, L5, I5)') 'Exit Status: ', exval, '# Iterations: ', &
          !num_iter, ' SAA: ', saa_flag, nsaa_spike
     
    !IF (exval > 0) THEN
       !simrad = fitspec_rad(1:n_rad_wvl) - fitres_rad(1:n_rad_wvl)
       !IF (use_lograd) THEN
          !fitspec_rad(1:n_rad_wvl) = EXP(fitspec_rad(1:n_rad_wvl))
          !simrad(1:n_rad_wvl) = EXP(simrad(1:n_rad_wvl))
          !fitres_rad(1:n_rad_wvl) = fitspec_rad(1:n_rad_wvl) -  simrad
       !END IF
       !avgres = SQRT(SUM((ABS(fitres_rad(1:n_rad_wvl)) / &
            !fitspec_rad(1:n_rad_wvl))**2.0)/n_rad_wvl)*100.0
       
       !WRITE(founit, '(4(A9, f8.2))') 'rms = ', rms, ' avgres = ', &
            !avgres, ' dfs = ', ozdfs, ' info. = ', ozinfo
       !WRITE(founit, '(A8,2I3,1P28E12.4)') 'Albedo: ', which_alb, nalb, eff_alb_init(1:nalb), &
            !eff_alb(1:nalb)
       !WRITE(founit, '(A8,2I3, 1P5E12.4)') 'Cloud:  ', which_cld, the_cld_flg, the_cfrac,     &
            !the_ctp, the_orig_cfr, the_orig_ctp, lambcld_refl
       !WRITE(founit, '(A8, I3, 1P36E10.3)') 'Aerosol:  ', actawin, the_ai, aerwavs(1:actawin), &
            !tropaod(1:actawin), tropsca(1:actawin), strataod(1:actawin), stratsca(1:actawin)
     
       !WRITE(founit, '(A31, 2I5)') 'Atmosphere and ozone profiles: ', nlay, ntp, nsfc
       !WRITE(founit, '(A)') '  #    P(mb)   Z(km)  T(K)     O3-ap O3-apstd   O3     STD     NSTD (DU)'
       !WRITE(founit, '(A2, 1x, F9.3, 2F8.3)') ' 0', atmosprof(1:3, 0)
       !DO i = 1, nlay 
          !WRITE(founit, '(I2, 1X, F9.3, 7F8.3)') i, atmosprof(1:3, i), ozprof_ap(i), &
               !ozprof_apstd(i), ozprof(i), ozprof_std(i), ozprof_nstd(i)
       !END DO
       !WRITE(founit, '(A24,4x,F8.3,8X,3F8.3)') ' Total Ozone: ', &
            !SUM(ozprof_ap(1:nsfc)), fitcol(1), dfitcol(1, 1), dfitcol(1, 2)
       !IF (ntp > 0) THEN
          !WRITE(founit, '(A24,4x,F8.3,8X,3F8.3)') ' Stratospheric Ozone: ',&
               !SUM(ozprof_ap(1:ntp)), fitcol(2), dfitcol(2, 1), dfitcol(2, 2)
          !WRITE(founit, '(A24,4x,F8.3,8X, 3F8.3)') ' Tropospheric Ozone: ',&
               !SUM(ozprof_ap(ntp+1:nsfc)), fitcol(3), dfitcol(3, 1), dfitcol(3, 2)          
       !ENDIF

       !IF (gaswrt) THEN
          !WRITE(founit, '(A, I5)') 'Fitted trace gases and uncertainty: ', nfgas
          !WRITE(founit, '(A)')  &
               !'    Var    Initial   A Priori A Priori Std  VCD      STD       NSTD      AMF      ACFRAC    AVGK(1)   AVGK(2)'
          !DO i = 1, ngas
              !IF (fgasidxs(i) > 0) THEN
                 !j = mask_fitvar_rad(fgasidxs(i))
                 !WRITE(founit, '(I2,1x,A6,1P10d10.2)') fgasidxs(i), fitvar_rad_str(j), tracegas(i, 1:10)
              !ENDIF
          !ENDDO
       !ENDIF

       !IF (ozwrtvar) THEN                      
           !WRITE(founit, '(A, I5)') 'Fitted variables and uncertainty: ', n_fitvar_rad
           !DO i = 1, n_fitvar_rad
              !WRITE(founit, '(I2, 1x,A6,1P3d10.2)') i, fitvar_rad_str(mask_fitvar_rad(i)), &
                   !fitvar_rad(mask_fitvar_rad(i)), fitvar_rad_std(mask_fitvar_rad(i)), &
                   !fitvar_rad_nstd(mask_fitvar_rad(i))
           !END DO
       !ENDIF
       
       !IF (ozwrtcorr) THEN
           !WRITE(founit, '(A, I5)') 'Correlation matrix: ', n_fitvar_rad
           !DO i = 1, n_fitvar_rad
              !correl(i) = 1.0
              !DO j = 1, i - 1
                 !correl(j) = covar(i, j) / SQRT(covar(i, i) * covar(j, j))
              !END DO
              !WRITE(founit, '(I2,1X,100f6.2)') i, correl(1:i)
           !END DO
       !ENDIF

       !IF (ozwrtcovar) THEN
          !WRITE(founit, '(A, 3I5)') 'Covariance matrix: ', nlay_fit, start_layer, end_layer
          !DO i = ozfit_start_index, ozfit_end_index
             !WRITE(founit, '(1p60d10.2)') covar(i, ozfit_start_index:ozfit_end_index)
          !END DO
       !ENDIF
       
       !IF (ozwrtavgk) THEN
           !WRITE(founit, '(A, 3I5)') 'Average kernel: ', nlay_fit, start_layer, end_layer
           !DO i = ozfit_start_index, ozfit_end_index
              !WRITE(founit, '(1p60d10.2)') avg_kernel(i, ozfit_start_index:ozfit_end_index)
           !END DO
       !ENDIF

       !IF (ozwrtfavgk) THEN
           !WRITE(founit, '(A, 3I5)') 'Average kernel1: ', nlay_fit, start_layer, end_layer
           !DO i = ozfit_start_index, ozfit_end_index
              !WRITE(founit, '(1p60d10.2)') favg_kernel(i, ozfit_start_index:ozfit_end_index)
           !END DO
       !ENDIF
       
       !IF (ozwrtcontri) THEN
          !WRITE(founit, '(A, I5)') 'Contribution Function: ', n_rad_wvl
          !DO i = 1, n_rad_wvl
             !WRITE(founit, '(1p60d10.2)') contri(ozfit_start_index:ozfit_end_index, i)
          !ENDDO
       !ENDIF

       !IF (ozwrtres) THEN
          !WRITE(founit, '(A, I5)') 'Fit residual: ', n_rad_wvl
          !DO i = 1, n_rad_wvl
             !WRITE(founit, '(f9.4, 1p2d12.4)') fitwavs(i), fitspec_rad(i), simrad(i)
          !END DO
       !ENDIF
             
    !END IF
    !WRITE(founit, *)
 !END IF
 
 !RETURN
!END SUBROUTINE gome_write_intermed


