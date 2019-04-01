  ! ***************************************************************************
  ! Author:  xiong liu
  ! Date  :  July 23, 2003
  ! Purpose: read input ozone profile variables and their bounds and add
  !          those variable to fitting control variables for radiance fit
  ! ***************************************************************************

SUBROUTINE read_ozprof_input (fit_ctrl_unit, fit_ctrl_file, pge_error_status )     

  USE OMSAO_precision_module
  USE ozprof_data_module,        ONLY: nlay, nlay_fit, ozstr, albstr, othstr,  &
       ozprof_start_index, ozprof_end_index, ozfit_start_index,                &
       ozfit_end_index, start_layer, end_layer, atmos_prof_fname,              &
       lcurve_write, lcurve_fname, ozwrtint, ozwrtint_fname,                   &
       lcurve_gcv, ptr_order, ptr_w0, ptr_w1, ptr_w2, ozabs_fname,             &
       aerosol, strat_aerosol, use_lograd, polcorr, t_fidx, t_lidx, tf_fidx,   &
       tf_lidx, nt_fit, nalb, nfalb, albfidx, albidx, use_oe, cloud, use_flns, &
       useasy, ndiv, albmax, albmin, do_multi_vza, do_lambcld,                 &
       ring_on_line, ring_convol, fit_atanring, degcorr, do_subfit, fgassidxs, &
       fgaspos, which_clima, which_alb, which_cld, use_logstate, radcalwrt,    &
       smooth_ozbc, do_ch2reso, coadd_after_b1ab, atmos_prof_fname, ozwrtcorr, &
       ozwrtcovar, ozwrtcontri, ozwrtres, ozwrtvar, ozwrtavgk,  ozwrtfavgk,    &
       degfname, biascorr, biasfname, do_tracewf, fgasidxs, ngas, gasidxs,     &
       algorithm_name, algorithm_version, atmwrt, gaswrt, nfgas, do_radinter,  &
       osind, osfind, slind, slfind, shind, shfind, rnind, rnfind, dcind,      &
       dcfind, isind, isfind, irind, irfind, oswins, slwins, shwins, rnwins,   &
       dcwins, iswins, irwins, nos, nsl, nsh, nrn, ndc, nis, nir, nothgrp,     &
       use_reg_presgrid, presgrid_fname, use_tropopause, adjust_trop_layer,    &
       pst0, fixed_ptrop, ntp0, which_biascorr, which_caloz, caloz_fname,      &
       ozwrtwf, ozwrtsnr, which_aerosol, scale_aod, scaled_aod, do_simu,       &
       pos_alb, toms_fwhm, ozcrs_alb_fname, nmom, maxmom, ngksec, maxgksec,    &
       maxgkmatc, ngkmatc, wrtring, wrtozcrs, ncldaer, cldaerstr, ecfrind,     &
       ecodind, ectpind, taodind, twaeind, saodind,  ecfrfind, ecodfind,      &
       ectpfind, taodfind, twaefind, saodfind, lambcld_initalb, scacld_initcod,&
       which_aperr, loose_aperr, min_serr, min_terr, which_toz, norm_tropo3,  &
       sprsind, sprsfind, nwfc, wfcstr, nfwfc, wfcmax, wfcmin, wfcfpix, wfclpix,&
       wfcidx, wfcfidx, so2zind, so2zfind, do_bothstep, do_twostep, &
       use_large_so2_aperr, use_effcrs, radc_msegsr, radc_nsegsr, radc_samprate,&
       radc_lambnd, hres_samprate, thealbidx, thewfcidx, do_simu_rmring,        &
       update_o3, update_sao3
       
  USE OMSAO_parameters_module,   ONLY: maxlay,  maxchlen, maxwin         
  USE OMSAO_indices_module,      ONLY: max_rs_idx,  max_calfit_idx, mxs_idx,   &
       n_max_fitpars, shi_idx, maxalb, solar_idx, so2_idx, bro_idx, hcho_idx,  &
       no2_t1_idx, maxoth, maxgrp, shift_offset, maxcldaer, maxwfc, so2v_idx
  USE OMSAO_variables_module,    ONLY: fitvar_rad_init,  fitvar_rad_saved,     &
       mask_fitvar_rad, n_fitvar_rad, lo_radbnd, up_radbnd, fitvar_rad_str,    &
       rad_identifier, refspec_fname, outdir, numwin, refdbdir, scnwrt,        &
       rmask_fitvar_rad, fitvar_rad_init_saved, winlim, fothvarpos, fitvar_rad_unit
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER,           INTENT (IN) :: fit_ctrl_unit
  CHARACTER (LEN=*), INTENT (IN) :: fit_ctrl_file

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER,           INTENT (OUT):: pge_error_status

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER                  :: i, j, k, iw, fidx, lidx, file_read_stat, idx, &
       swin, ewin, ntemp, nord, ntotp, thewin, theord, np, fstlay, lstlay, &
       fstlayT, lstlayT
  INTEGER, DIMENSION (maxoth, 2)      :: tmpwins
  INTEGER, DIMENSION (maxwin, maxoth) :: tmpind, tmpfind
  CHARACTER (LEN=maxchlen)            :: tmpchar
  CHARACTER (LEN=3)                   :: idxchar
  CHARACTER (LEN=6)                   :: idxchar1
  REAL      (KIND=dp)                 :: vartmp, lotmp, uptmp, mnoz, minoz, maxoz, mnT, minT, maxT
  LOGICAL                             :: vart

  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=17), PARAMETER :: modulename = 'read_ozprof_input'

  ! ========================
  ! Error handling variables
  ! ========================
  INTEGER :: errstat

  ! =================================
  ! External OMI and Toolkit routines
  ! =================================
  INTEGER :: OMI_SMF_setmsg

  use_oe = .true.
  ! -------------------------
  ! Open fitting control file
  ! -------------------------
  OPEN ( UNIT=fit_ctrl_unit, FILE=TRIM(ADJUSTL(fit_ctrl_file)), &
       STATUS='OLD', IOSTAT=errstat)

  IF ( errstat /= pge_errstat_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_open_fitctrl_file, &
          TRIM(ADJUSTL(fit_ctrl_file)), modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  END IF

  READ (fit_ctrl_unit, *) 
  READ (fit_ctrl_unit, '(A)') algorithm_name
  READ (fit_ctrl_unit, '(A)') algorithm_version
  algorithm_name = TRIM(ADJUSTL(algorithm_name))
  algorithm_version = TRIM(ADJUSTL(algorithm_version))
  
  ! -------------------------------
  ! Read general control variables
  ! -------------------------------
  READ (fit_ctrl_unit, *) do_multi_vza
  IF (do_multi_vza) THEN
     WRITE(*, *) 'Only effective viewing geometry is computed for a ground pixel!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF
  IF (do_ch2reso) do_multi_vza = .FALSE.

  READ (fit_ctrl_unit, *) do_radinter
  READ (fit_ctrl_unit, *) do_simu, do_simu_rmring
  IF ( .NOT. do_simu ) do_simu_rmring = .FALSE.
  READ (fit_ctrl_unit, *) do_twostep, do_bothstep, use_large_so2_aperr
  READ (fit_ctrl_unit, *) do_tracewf 
  READ (fit_ctrl_unit, *) 
  READ (fit_ctrl_unit, *) atmwrt, scnwrt, ozwrtint, gaswrt, ozwrtvar, ozwrtcorr, &
       ozwrtcovar, ozwrtavgk, ozwrtfavgk, ozwrtcontri, ozwrtres, ozwrtwf, ozwrtsnr, &
       wrtring, wrtozcrs
  ozwrtint_fname = TRIM(ADJUSTL(outdir)) // 'inter_' // rad_identifier // '.dat'
  READ (fit_ctrl_unit, *) 
  READ (fit_ctrl_unit, '(A)') atmos_prof_fname

  !  Calibration options
  READ (fit_ctrl_unit, *) 
  READ (fit_ctrl_unit, *) radcalwrt, which_caloz
  READ (fit_ctrl_unit, '(A)') caloz_fname
  IF (radcalwrt .AND. .NOT. do_simu) ozwrtres = .FALSE.
  IF (radcalwrt .AND. do_simu) THEN
     ozwrtavgk = .FALSE.; ozwrtfavgk = .FALSE. ; ozwrtint = .FALSE.
     ozwrtcorr = .FALSE.; ozwrtcovar = .FALSE.;  ozwrtcontri = .FALSE.
     ozwrtwf   = .FALSE.; ozwrtsnr   = .FALSE.;  ozwrtres = .TRUE.; ozwrtvar = .FALSE.
  ENDIF

  READ (fit_ctrl_unit, *) use_lograd
  READ (fit_ctrl_unit, *) use_logstate
  READ (fit_ctrl_unit, *) use_flns
  READ (fit_ctrl_unit, *) ring_on_line, ring_convol, fit_atanring
  READ (fit_ctrl_unit, *) smooth_ozbc
  READ (fit_ctrl_unit, *) biascorr, degcorr  
  READ (fit_ctrl_unit, *) which_biascorr
  READ (fit_ctrl_unit, '(A)') biasfname
  READ (fit_ctrl_unit, '(A)') degfname
  READ (fit_ctrl_unit, *) do_subfit
  READ (fit_ctrl_unit, *) lcurve_gcv
  READ (fit_ctrl_unit, *) lcurve_write
  lcurve_fname = TRIM(ADJUSTL(outdir)) // 'lcurve_' // rad_identifier // '.dat'
  READ (fit_ctrl_unit, *) ptr_order, ptr_w0, ptr_w1, ptr_w2
  READ (fit_ctrl_unit, *) which_clima
  READ (fit_ctrl_unit, *) which_aperr
  READ (fit_ctrl_unit, *) which_toz
  READ (fit_ctrl_unit, *) update_o3
  READ (fit_ctrl_unit, *) update_sao3
 
  IF (which_clima > 14 .OR. which_clima <= 0) THEN
     WRITE(*, *) modulename, ' No such ozone profile climatology!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF

  IF (which_aperr > 14 .OR. which_aperr <= 0 .or. which_aperr  ==7 .or. which_aperr ==8) THEN
     WRITE(*, *) modulename, ' No such ozone profile a priori error!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF

  IF (which_toz > 2 .OR. which_toz < 0) THEN
     WRITE(*, *) modulename, ' No such total ozone field!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF

  IF ( which_toz == 0) THEN
     IF ( which_clima == 5 .or. which_clima == 6 .or. which_aperr == 5 ) THEN 
     WRITE(*, *) modulename, ' Total ozone is needed to use this climatology!!!'
     pge_error_status = pge_errstat_error; RETURN
     ENDIF
     
     IF ( update_o3 == .TRUE. .or. update_sao3 == .TRUE. ) THEN
     WRITE(*, *) modulename, ' Total ozone is needed to update o3 or sao3 !!!'
     pge_error_status = pge_errstat_error; RETURN 
     ENDIF
  ELSE IF ( which_toz /= 0 ) THEN 
  ! toz might be need when doing radiometric calibration and we need total ozone to normalize the profile
  !   IF ( (which_clima /= 5 .and. which_clima /= 6) .and. (which_aperr /= 5) ) THEN 
  !         which_toz = 0
  !        WRITE(*, *) modulename, ' which_toz is set to be 0 because clima is not IUP or V* !!!'      
  !        IF ( update_o3 == .TRUE. .or. update_sao3 == .TRUE. ) THEN
  !             update_o3 = .FALSE. ; update_sao3 = .FALSE.
  !        WRITE(*, *) modulename, ' updateo3 is set to .FASE. because clima is not IUP or V* !!!'      
  !        ENDIF
  !   ENDIF     
  ENDIF 

      
  READ (fit_ctrl_unit, *) loose_aperr, min_serr, min_terr
  READ (fit_ctrl_unit, *) norm_tropo3
  READ (fit_ctrl_unit, *) which_alb

  IF (which_alb > 5) THEN
     WRITE(*, *) modulename, ' No such albedo database!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF

  READ (fit_ctrl_unit, *) which_cld
  IF (which_cld > 4) THEN
     WRITE(*, *) modulename, ' No such cloud option!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF

 
  READ (fit_ctrl_unit, *) aerosol, strat_aerosol
  IF( aerosol ) print *,'aerosol is considered !!!' 
  IF (.NOT. aerosol) strat_aerosol = .FALSE.
  READ (fit_ctrl_unit, *)
  READ (fit_ctrl_unit, *) which_aerosol, scale_aod, scaled_aod
  READ (fit_ctrl_unit, *) cloud, do_lambcld, lambcld_initalb, scacld_initcod
  IF (.NOT. cloud) do_lambcld = .FALSE.
  !IF (.NOT. do_lambcld) THEN
  !   WRITE(*, *) modulename, ': Clouds must be assummed Lambertian!!!'
  !   pge_error_status = pge_errstat_error; RETURN
  !ENDIF
  READ (fit_ctrl_unit, *) useasy
  READ (fit_ctrl_unit, *) nmom
  IF (nmom > maxmom) THEN
     WRITE(*, *) modulename, ': Input # of phase moments exceeds maxmom: ', maxmom
     WRITE(*, *) modulename, ': Reduce it or increase maxmom!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF
  IF (.NOT. aerosol .AND. (.NOT. cloud .OR. do_lambcld) ) THEN
     nmom = 2
  ENDIF
  READ (fit_ctrl_unit, *) use_effcrs
  IF (.NOT. use_effcrs) do_radinter = .FALSE.
  READ (fit_ctrl_unit, *) hres_samprate
  IF (hres_samprate < 0.01) THEN
     WRITE(*, *) modulename, ':hres_samprate must be > 0.01!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF

  IF ( MOD(hres_samprate * 100, 1.0) /= 0) THEN
     WRITE(*, *) modulename, ': hres_samprate must be multiples of 0.01!!!'
     hres_samprate = NINT(hres_samprate * 100) / 100.
     WRITE(*, *) modulename, ': hres_samprate is reset to: ', hres_samprate
  ENDIF
       
  ! Read parameters for specifying sampling rate for specified # of spectral region
  ! Only valid if use_effcrs is set to false
  READ (fit_ctrl_unit, *) radc_nsegsr
  IF (radc_nsegsr > radc_msegsr) THEN
     WRITE(*, *) modulename, ': Need to increase radc_msegsr!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF
  IF (radc_nsegsr <= 0) THEN
     WRITE(*, *) modulename, ': radc_nsegsr must be > 0!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF
  READ (fit_ctrl_unit, *) (radc_lambnd(i), i = 1, radc_nsegsr)
  DO i = 2, radc_nsegsr
     IF (radc_lambnd(i) <= radc_lambnd(i-1)) THEN
        WRITE(*, *) modulename, ': Need to be increasing order!!!'
        pge_error_status = pge_errstat_error; RETURN
     ENDIF
  ENDDO
  IF (radc_lambnd(1) > winlim(1, 1) - 5) THEN
     radc_lambnd(1) = winlim(1, 1) - 5
  ENDIF
  READ (fit_ctrl_unit, *) (radc_samprate(i), i = 1, radc_nsegsr)
  IF ( MINVAL(radc_samprate(1:radc_nsegsr)) < hres_samprate ) THEN
     WRITE(*, *) modulename, ': radc_samprate must be > hres_samprate!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF
  
  READ (fit_ctrl_unit, *) polcorr 
  IF (polcorr > 5 .OR. polcorr == 2) THEN
     WRITE(*, *) modulename, ': No such polarization correction option!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF
  IF (polcorr == 0 .OR. polcorr >= 3 ) THEN
     useasy = .FALSE.
     ngksec = maxgksec; ngkmatc = maxgkmatc
  ELSE
     ngksec = 1; ngkmatc = 1
  ENDIF

  ! -------------------------------------
  ! Read ozone profile control variables
  ! -------------------------------------
  REWIND ( fit_ctrl_unit )
  CALL skip_to_filemark ( fit_ctrl_unit, ozstr, tmpchar, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, &
          TRIM(ADJUSTL(fit_ctrl_file)), modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  END IF

  !xliu (02/08/2007): modify the way of reading input for atmospheric profiles
  !                   add more options for atmospheric layering scheme
  READ (fit_ctrl_unit, '(A)') ozabs_fname
  ozabs_fname = TRIM(ADJUSTL(refdbdir)) // ozabs_fname
  READ (fit_ctrl_unit, *) nlay
  IF (nlay > maxlay) THEN
     WRITE(*, *) modulename, ' : # of layers exceed maximum # of layers allowed!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF 
  READ (fit_ctrl_unit, *) ndiv
  READ (fit_ctrl_unit, *) fstlay, lstlay
  READ (fit_ctrl_unit, *) mnoz, minoz, maxoz
  READ (fit_ctrl_unit, *) fstlayT, lstlayT
  READ (fit_ctrl_unit, *) mnT, minT, maxT
  READ (fit_ctrl_unit, *) use_reg_presgrid
  READ (fit_ctrl_unit, '(A)') presgrid_fname
  READ (fit_ctrl_unit, *) use_tropopause
  READ (fit_ctrl_unit, *) fixed_ptrop, pst0, ntp0
  READ (fit_ctrl_unit, *) adjust_trop_layer  
  
  ! -----------------------------------------------
  ! Read initial ozone profile variables
  ! -----------------------------------------------
  idx  = shift_offset + max_rs_idx
  ozprof_start_index = idx + 1
  ozprof_end_index = idx + nlay
  t_fidx = ozprof_start_index + maxlay
  t_lidx = ozprof_end_index + maxlay
  fitvar_rad_unit(ozprof_start_index:ozprof_end_index) = 'DU'
  fitvar_rad_unit(t_fidx:t_lidx) = 'K'
  
  fitvar_rad_init(ozprof_start_index:ozprof_end_index) = mnoz
  lo_radbnd(ozprof_start_index:ozprof_end_index)       = minoz
  up_radbnd(ozprof_start_index:ozprof_end_index)       = maxoz
  DO i = 1, fstlay-1
     lo_radbnd(idx+i) = mnoz
     up_radbnd(idx+i) = mnoz
  ENDDO
  DO i = lstlay+1, nlay
     lo_radbnd(idx+i) = mnoz
     up_radbnd(idx+i) = mnoz
  ENDDO

  fitvar_rad_init(t_fidx:t_lidx) = mnT
  lo_radbnd(t_fidx:t_lidx)       = minT
  up_radbnd(t_fidx:t_lidx)       = maxT
  DO i = 1, fstlayT-1
     lo_radbnd(t_fidx+i-1) = mnT
     up_radbnd(t_fidx+i-1) = mnT
  ENDDO
  DO i = lstlayT+1, nlay
     lo_radbnd(t_fidx+i-1) = mnT
     up_radbnd(t_fidx+i-1) = mnT
  ENDDO
  DO i = 1, nlay
     WRITE(fitvar_rad_str(idx+i), '(A2, I2.2)') 'oz', i
  ENDDO
  DO i = fstlayT, lstlayT
     WRITE(fitvar_rad_str(t_fidx+i-1), '(A2, I2.2)') 'tt', i
  ENDDO
    
!  ozpfpars: DO i = 1, nlay
!     READ (fit_ctrl_unit, *, IOSTAT=errstat) vartmp, lotmp, uptmp, vart
!     IF ( errstat /= pge_errstat_ok ) THEN
!        errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, &
!             TRIM(ADJUSTL(fit_ctrl_file)), modulename, 0)
!        WRITE(*, *) modulename, ' : Error in reading initial ozone variables!!!'
!        pge_error_status = pge_errstat_error; RETURN
!     END IF
!     ! ---------------------------------------------------------
!     ! Check for consitency of bounds and adjust where necessary
!     ! ---------------------------------------------------------
!     IF ( lotmp > vartmp .OR. uptmp < vartmp ) THEN
!        lotmp = vartmp ; uptmp = vartmp
!     END IF
!     IF ( lotmp == uptmp .AND. lotmp /= vartmp ) THEN
!        uptmp = vartmp ; lotmp = vartmp
!     END IF
!     
!     fitvar_rad_init(idx + i) = vartmp
!     WRITE(fitvar_rad_str(idx+i), '(A2, I2.2)') 'oz', i
!     
!     lo_radbnd(idx + i) = lotmp
!     up_radbnd(idx + i) = uptmp  
!     
!     ! determine whether to vary temperature
!     IF (vart .AND. lotmp < uptmp) THEN
!        fitvar_rad_init(idx + i + maxlay) = 200.0  ! rechange after preparing atmos.
!        lo_radbnd(idx + i + maxlay) = 0.
!        up_radbnd(idx + i + maxlay) = 400.
!        WRITE(fitvar_rad_str(idx+i+maxlay), '(A2, I2)') 'tmp', i
!
!     ELSE
!        fitvar_rad_init(idx + i+ maxlay) = 0. 
!        lo_radbnd(idx + i + maxlay) = 0.
!        up_radbnd(idx + i + maxlay) = 0.
!     END IF
!
!  END DO ozpfpars

  ! -------------------------------------------------------------
  ! Add unfixed ozone profile variables to the variable list 
  ! -------------------------------------------------------------
  ntemp = n_fitvar_rad
  ozfit_start_index = n_fitvar_rad + 1
  DO i = idx + 1, idx + nlay
     IF ( lo_radbnd(i) < up_radbnd(i) ) THEN
        n_fitvar_rad = n_fitvar_rad + 1
        mask_fitvar_rad(n_fitvar_rad) = i
        rmask_fitvar_rad(i) = n_fitvar_rad        
     END IF
     IF (n_fitvar_rad == ntemp + 1) THEN
        start_layer = i - idx
     END IF
  END DO
  ozfit_end_index = n_fitvar_rad
  nlay_fit = ozfit_end_index - ozfit_start_index + 1
  end_layer = start_layer + nlay_fit - 1
  
  ntemp = n_fitvar_rad 
  tf_fidx = 0; tf_lidx = 0; nt_fit = 0
  
  ! Add unfixed temperature variables to the variable list
  DO i = idx + 1 + maxlay, idx + nlay + maxlay
     IF ( lo_radbnd(i) < up_radbnd(i) ) THEN
        n_fitvar_rad = n_fitvar_rad + 1
        mask_fitvar_rad(n_fitvar_rad) = i
        rmask_fitvar_rad(i) = n_fitvar_rad  
        IF (tf_fidx == 0) tf_fidx = n_fitvar_rad        
     END IF

  END DO
  IF (tf_fidx > 0) THEN
     tf_lidx = n_fitvar_rad
     nt_fit = tf_lidx - tf_fidx + 1
  ENDIF
  
  ! -----------------------------------------------
  ! Read wavelength dependent surface albedo terms
  ! -----------------------------------------------
  REWIND ( fit_ctrl_unit )
  CALL skip_to_filemark ( fit_ctrl_unit, albstr, tmpchar, file_read_stat)
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, &
          TRIM(ADJUSTL(fit_ctrl_file)), modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  END IF

  READ (fit_ctrl_unit, *) pos_alb, toms_fwhm
  READ (fit_ctrl_unit, '(A)') ozcrs_alb_fname
  READ (fit_ctrl_unit, *) nalb
  IF (nalb > maxalb) THEN
     WRITE(*, *) modulename, ' : # of albedo terms exceeds allowed  ', maxalb
     pge_error_status = pge_errstat_error; RETURN
  ELSE IF (nalb < 1) THEN
     WRITE(*, *) modulename, ' : Need to specify at least 1 albedo term!!!'
     pge_error_status = pge_errstat_error; RETURN     
  ENDIF
  
  idx = shift_offset + max_rs_idx + maxlay * 2 
  albidx = idx + 1
  albmin= 0.0; albmax = 0.0
  albpars: DO i = 1, nalb
     READ (fit_ctrl_unit, *, IOSTAT=errstat) idxchar1, vartmp, &
          lotmp, uptmp, albmin(i), albmax(i)
     IF ( errstat /= pge_errstat_ok ) THEN
		errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, &
             TRIM(ADJUSTL(fit_ctrl_file)), modulename, 0)
		WRITE(*, *) modulename, ' : Error in reading initial albedo variables!!!'
		pge_error_status = pge_errstat_error; RETURN
     END IF

     ! ---------------------------------------------------------
     ! Check for consitency of bounds and adjust where necessary
     ! ---------------------------------------------------------
     IF ( lotmp > vartmp .OR. uptmp < vartmp ) THEN
		lotmp = vartmp ; uptmp = vartmp
     END IF
     IF ( lotmp == uptmp .AND. lotmp /= vartmp ) THEN
		uptmp = vartmp ; lotmp = vartmp
     END IF
    
     fitvar_rad_init(idx + i) = vartmp
     fitvar_rad_str(idx+i) = TRIM(ADJUSTL(idxchar1))
     lo_radbnd(idx + i) = lotmp
     up_radbnd(idx + i) = uptmp
  END DO albpars

  ! -------------------------------------------------------------
  ! Add albedo variable to the variable list
  ! -------------------------------------------------------------
  albfidx = 0; thealbidx = 0
  DO i = idx + 1, idx + nalb
     IF ( lo_radbnd(i) < up_radbnd(i) ) THEN
        n_fitvar_rad = n_fitvar_rad + 1
        mask_fitvar_rad(n_fitvar_rad) = i
        rmask_fitvar_rad(i) = n_fitvar_rad  
        IF (albfidx == 0) albfidx = n_fitvar_rad
        IF ( fitvar_rad_str(i)(4:4)== '0') thealbidx = n_fitvar_rad
     ENDIF
  ENDDO
  IF (albfidx > 0) THEN
     nfalb = n_fitvar_rad - albfidx + 1
     thealbidx = thealbidx - albfidx + 1
  ENDIF

  ! ------------------------------------------------
  ! Read wavelength dependent cloud fraction terms
  ! ------------------------------------------------
  REWIND ( fit_ctrl_unit )
  CALL skip_to_filemark ( fit_ctrl_unit, wfcstr, tmpchar, file_read_stat)
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, &
          TRIM(ADJUSTL(fit_ctrl_file)), modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  END IF

  READ (fit_ctrl_unit, *) nwfc
  IF (nwfc > maxwfc) THEN
     WRITE(*, *) modulename, ' : # of wavelength dependent terms exceed allowed  ', maxwfc
     pge_error_status = pge_errstat_error; RETURN
  !ELSE IF (nwfc < 1) THEN
  !   WRITE(*, *) modulename, ' : Need to specify at least 1 cloud fraction term!!!'
  !   pge_error_status = pge_errstat_error; RETURN     
  ELSE IF (nwfc < 1) THEN
     nwfc = 0
  ENDIF

  idx = idx + maxalb 
  wfcidx = idx + 1
  wfcmin= 0.0; wfcmax = 0.0
  wfcpars: DO i = 1, nwfc
     READ (fit_ctrl_unit, *, IOSTAT=errstat) idxchar1, vartmp, &
          lotmp, uptmp, wfcmin(i), wfcmax(i)
     IF ( errstat /= pge_errstat_ok ) THEN
		errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, &
             TRIM(ADJUSTL(fit_ctrl_file)), modulename, 0)
		WRITE(*, *) modulename, ' : Error in reading initial cloud fraction variables!!!'
		pge_error_status = pge_errstat_error; RETURN
     END IF

     ! ---------------------------------------------------------
     ! Check for consitency of bounds and adjust where necessary
     ! ---------------------------------------------------------
     IF ( lotmp > vartmp .OR. uptmp < vartmp ) THEN
		lotmp = vartmp ; uptmp = vartmp
     END IF
     IF ( lotmp == uptmp .AND. lotmp /= vartmp ) THEN
		uptmp = vartmp ; lotmp = vartmp
     END IF
    
     fitvar_rad_init(idx + i) = vartmp
     fitvar_rad_str(idx+i) = TRIM(ADJUSTL(idxchar1))
     lo_radbnd(idx + i) = lotmp
     up_radbnd(idx + i) = uptmp
  END DO wfcpars

  ! ---------------------------------------------------------------------
  ! Add wavelength-dependent cloud fraction variables to the variable list
  ! ----------------------------------------------------------------------
  wfcfidx = 0; thewfcidx=0
  DO i = idx + 1, idx + nwfc
     IF ( lo_radbnd(i) < up_radbnd(i) ) THEN
        n_fitvar_rad = n_fitvar_rad + 1
        mask_fitvar_rad(n_fitvar_rad) = i
        rmask_fitvar_rad(i) = n_fitvar_rad  
        IF (wfcfidx == 0) wfcfidx = n_fitvar_rad
        IF ( fitvar_rad_str(i)(4:4)== '0') thewfcidx = n_fitvar_rad
     END IF
  END DO
  IF (wfcfidx > 0) THEN
     nfwfc = n_fitvar_rad - wfcfidx + 1
     thewfcidx = thewfcidx - albfidx + 1
  ENDIF

  ! -------------------------------
  ! Read cloud/aerosol variables
  ! -------------------------------
  REWIND ( fit_ctrl_unit )
  CALL skip_to_filemark ( fit_ctrl_unit, cldaerstr, tmpchar, file_read_stat)
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, &
          TRIM(ADJUSTL(fit_ctrl_file)), modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  END IF
  READ (fit_ctrl_unit, *) ncldaer
  IF (ncldaer > maxcldaer) THEN
     WRITE(*, *) modulename, ' : Increase maxcldaer in OMSAO_indices...!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF

  idx = wfcidx + maxwfc - 1
   
  DO i = 1, ncldaer
     READ (fit_ctrl_unit, *, IOSTAT=errstat) idxchar1, vartmp, lotmp, uptmp
     IF ( errstat /= pge_errstat_ok ) THEN
		errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, &
             TRIM(ADJUSTL(fit_ctrl_file)), modulename, 0)
		WRITE(*, *) modulename, ' : Error in reading initial albedo variables!!!'
		pge_error_status = pge_errstat_error; RETURN
     END IF

     ! ---------------------------------------------------------
     ! Check for consitency of bounds and adjust where necessary
     ! ---------------------------------------------------------
     IF ( lotmp > vartmp .OR. uptmp < vartmp ) THEN
		lotmp = vartmp ; uptmp = vartmp
     END IF
     IF ( lotmp == uptmp .AND. lotmp /= vartmp ) THEN
		uptmp = vartmp ; lotmp = vartmp
     END IF

     IF ( .NOT. cloud .AND. i > 0 .AND. i < 4) THEN   ! Cannot fit any cloud variables
        vartmp = 0.0; lotmp = 0.0; uptmp = 0.0
     ENDIF
     ! Disalble this cloud fraction if wavelength dependent cloud fraction is selected
     IF (nwfc > 0 .AND. i == 1) THEN  
        vartmp = 0.0; lotmp = 0.0; uptmp = 0.0
     ENDIF

     IF ( i == 8 ) THEN ! For SO2Z, SO2V needs to be selected
        j = max_calfit_idx + (so2v_idx - 1) * mxs_idx + 2  ! Index for SO2V
        IF (rmask_fitvar_rad(j) == 0) THEN
           vartmp = 0.0; lotmp = 0.0; uptmp = 0.0          ! Disable SO2Z if SO2V is not selected
        ELSE                                               
           IF (vartmp == 0.0 .AND. lotmp >= uptmp) THEN 
              vartmp = 5.0; lotmp = 5.0; uptmp = 5.0       ! Initialize to 5 km if not initialized or set to zero
           ENDIF
        ENDIF
     ENDIF

     IF ( do_lambcld .AND. i == 2 ) THEN              ! Cannot fit cloud optical thickness
        vartmp = 0.0; lotmp = 0.0; uptmp = 0.0
     ENDIF
     IF ( .NOT. aerosol .AND. i > 3 .AND. i < 7) THEN ! Cannot fit aerosol variables
        vartmp = 0.0; lotmp = 0.0; uptmp = 0.0
     ENDIF
     IF ( .NOT. strat_aerosol .AND. i == 6) THEN      ! Cannot fit stratospheric aerosols
        vartmp = 0.0; lotmp = 0.0; uptmp = 0.0
     ENDIF
     ! Disable fitting cloud top pressure (not implemented)
     ! Cloud optical thickness is directly fitted either from longer wavelengths, 
     ! which is enabled when using scattering clouds, or a cloud fraction is directly fitted 
     ! and when cloud fraction is greater than 1 (need to exchange fitting variable between cloud 
     ! fraction and cloud optical thickness)
     IF ( i == 2 .OR. i == 3) THEN                   
        vartmp = 0.0; lotmp = 0.0; uptmp = 0.0
     ENDIF
    
     fitvar_rad_init(idx + i) = vartmp
     fitvar_rad_str(idx+i) = TRIM(ADJUSTL(idxchar1))
     lo_radbnd(idx + i) = lotmp
     up_radbnd(idx + i) = uptmp          
  ENDDO
  ecfrind = idx + 1; ecodind = idx + 2; ectpind = idx + 3
  taodind = idx + 4; twaeind = idx + 5; saodind = idx + 6
  sprsind = idx + 7; so2zind = idx + 8
  
  ecfrfind = 0; ecodfind = 0; ectpfind = 0
  taodfind = 0; twaefind = 0; saodfind = 0
  sprsfind = 0; so2zfind = 0
  ! -------------------------------------------------------------
  ! Add Cloud/Aerosol variables to the variable list
  ! -------------------------------------------------------------
  DO i = idx + 1, idx + ncldaer
     j = i -idx
     IF ( lo_radbnd(i) < up_radbnd(i) ) THEN
        n_fitvar_rad = n_fitvar_rad + 1
        mask_fitvar_rad(n_fitvar_rad) = i
        rmask_fitvar_rad(i) = n_fitvar_rad 
        IF      (j == 1 ) THEN
           ecfrfind = n_fitvar_rad
        ELSE IF (j == 2) THEN
           ecodfind = n_fitvar_rad
        ELSE IF (j == 3) THEN
           ectpfind = n_fitvar_rad
           fitvar_rad_unit(i) = 'mb'
        ELSE IF (j == 4) THEN
           taodfind = n_fitvar_rad
        ELSE IF (j == 5) THEN
           twaefind = n_fitvar_rad
        ELSE IF (j == 6) THEN
           saodfind = n_fitvar_rad
        ELSE IF (j == 7) THEN
           sprsfind = n_fitvar_rad
           fitvar_rad_unit(i) = 'mb'
        ELSE IF (j == 8) THEN
           so2zfind = n_fitvar_rad
           fitvar_rad_unit(i) = 'km'
        ENDIF
     END IF
  END DO
  
  ! Read other parameters (for multiple window)
  REWIND ( fit_ctrl_unit )
  CALL skip_to_filemark ( fit_ctrl_unit, othstr, tmpchar, file_read_stat)
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, &
          TRIM(ADJUSTL(fit_ctrl_file)), modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  END IF    
  READ (fit_ctrl_unit, *) nothgrp 
  IF (nothgrp > maxgrp) THEN
     WRITE(*, *) modulename, ' : Increase maxgrp in OMSAO_indices...!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF
   
  k = wfcidx + maxwfc + maxcldaer - 1 
  IF (do_subfit) THEN
     ntotp = numwin * maxoth; np = numwin
  ELSE
     ntotp = maxoth; np = 1
  ENDIF

  DO i = 1, nothgrp        ! for each group of parameters

     nord = 0; tmpind = 0; tmpfind = 0; tmpwins = 0

     DO j = 1, maxoth   ! for each order of parameters

        READ (fit_ctrl_unit, *, IOSTAT=errstat) idxchar1, vartmp, lotmp, uptmp, swin, ewin
        IF ( errstat /= pge_errstat_ok ) THEN
           errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, &
                TRIM(ADJUSTL(fit_ctrl_file)), modulename, 0)
           WRITE(*, *) modulename, ' : Error in reading other variables!!!'
           pge_error_status = pge_errstat_error; RETURN
        END IF

        ! special conditions for fit_atanring
        IF (i == 4 .AND. fit_atanring) THEN
           swin = 1;   ewin = 1
           IF (j > 3) THEN
              lotmp = 0.0; uptmp = 0.0
           ELSE
              lotmp = -1.0D99; uptmp = 1.0D99
           ENDIF
           IF (j == 1) THEN
              vartmp = 0.4
           ELSE IF (j == 2) THEN
              vartmp = 303.
           ELSE IF (j == 3) THEN
              vartmp = 3.0
           ELSE
              vartmp = 0.0
           ENDIF
        ENDIF

        IF (swin > ewin) THEN
           ntemp = swin; swin = ewin; ewin = ntemp
        ENDIF
        IF (swin < 1 .AND. ewin < 1) THEN
           swin = 1; ewin = numwin
        ENDIF
        IF (swin < 1) swin = 1
        IF (ewin > numwin) ewin = numwin       
        tmpwins(j, 1) = swin; tmpwins(j, 2) = ewin

        ! ---------------------------------------------------------
        ! Check for consitency of bounds and adjust where necessary
        ! ---------------------------------------------------------
        IF (lotmp > vartmp .OR. uptmp < vartmp ) THEN
           lotmp = vartmp ; uptmp = vartmp
        END IF
        IF ( lotmp == uptmp .AND. lotmp /= vartmp ) THEN
           uptmp = vartmp ; lotmp = vartmp
        END IF
        
        IF (do_subfit) THEN
           lo_radbnd (k + 1 : k + numwin) = 0.0
           up_radbnd (k + 1 : k + numwin) = 0.0
           fitvar_rad_init (k + 1 : k + numwin) = 0.0
           IF (i == 1 .OR. i == 3) fitvar_rad_unit(k+1:k+numwin) = 'nm'
                      
           fitvar_rad_init(k + swin : k + ewin) = vartmp
           fitvar_rad_str (k + swin : k + ewin) = TRIM(ADJUSTL(idxchar1))
           lo_radbnd      (k + swin : k + ewin) = lotmp
           up_radbnd      (k + swin : k + ewin) = uptmp
           IF (i == 1 .OR. i == 3) fitvar_rad_unit(k+swin:k+ewin) = 'nm'

           tmpind(1:numwin, j) = k + (/(idx, idx = 1, numwin)/)
           k = k + numwin

        ELSE
           fitvar_rad_init(k + 1) = vartmp
           fitvar_rad_str (k + 1) = TRIM(ADJUSTL(idxchar1))
           IF (i == 1 .OR. i == 3) fitvar_rad_unit(k + 1) = 'nm'
           lo_radbnd(k + 1) = lotmp
           up_radbnd(k + 1) = uptmp
           k = k + 1; tmpind(1, j) = k

           ! The windows must be the same for all the orders
           IF (j > 1) THEN
              tmpwins(j, :) = tmpwins(1, :)
           ENDIF
        ENDIF

     END DO

     ! If lower order, then no higher orders
     DO j = 1, maxoth-1
        DO iw = 1, np
           fidx = tmpind(iw, j)
           IF (lo_radbnd(fidx) == up_radbnd(fidx)) THEN
              lo_radbnd(tmpind(iw, j+1:maxoth)) = 0.0
              up_radbnd(tmpind(iw, j+1:maxoth)) = 0.0
              fitvar_rad_init(tmpind(iw, j+1:maxoth)) = 0.0
           ENDIF
        ENDDO
     ENDDO

     ntemp = n_fitvar_rad; idx = tmpind(1, 1) - 1
     DO j = 1, ntotp
        IF (lo_radbnd(idx + j) < up_radbnd(idx + j)) THEN
           n_fitvar_rad = n_fitvar_rad + 1
           mask_fitvar_rad(n_fitvar_rad) = idx + j
           rmask_fitvar_rad(idx + j) = n_fitvar_rad  
           theord = CEILING(1.0 * j / np); thewin = j - (theord - 1) * np; 
           tmpfind(thewin, theord) = n_fitvar_rad 
        ENDIF
     ENDDO
     IF (n_fitvar_rad > ntemp) nord =  theord
     
     IF (i == 1) THEN
        nos = nord; oswins = tmpwins; osind = tmpind; osfind = tmpfind
     ELSE IF ( i == 2) THEN
        nsl = nord; slwins = tmpwins; slind = tmpind; slfind = tmpfind
     ELSE IF (i == 3)  THEN
        nsh = nord; shwins = tmpwins; shind = tmpind; shfind = tmpfind
     ELSE IF (i == 4)  THEN
        nrn = nord; rnwins = tmpwins; rnind = tmpind; rnfind = tmpfind
     ELSE IF (i == 5 ) THEN 
        ndc = nord; dcwins = tmpwins; dcind = tmpind; dcfind = tmpfind
     ELSE IF (i == 6)  THEN
        nis = nord; iswins = tmpwins; isind = tmpind; isfind = tmpfind
     ELSE IF (i == 7) THEN
        nir = nord; irwins = tmpwins; irind = tmpind; irfind = tmpfind
     ENDIF

     !WRITE(*, *) i, nord
     !WRITE(*, '(8I5)') tmpind(1:np,  1:maxoth)
     !WRITE(*, '(8I5)') tmpfind(1:np, 1:maxoth)
  ENDDO
  !DO i = 1, n_fitvar_rad
  !   print *, i, fitvar_rad_str(mask_fitvar_rad(i))
  !ENDDO
  !WRITE(*, '(7I5)') nos, nsl, nsh, nrn, ndc, nis, nir

  ! get indices for auxiliary variables in the final fitted array
  fgasidxs = 0; nfgas = 0
  DO i = 1, ngas
     fidx =  max_calfit_idx + (gasidxs(i) - 1) * mxs_idx + 1; lidx = fidx + 2
     fitvar_rad_unit(fidx:lidx) = 'molecumes cm^-2'
     fgasidxs(i) = MAXVAL(rmask_fitvar_rad(fidx:lidx))
     fgassidxs(i) = rmask_fitvar_rad(shift_offset + gasidxs(i))
     IF (fgasidxs(i) > 0) THEN
        nfgas = nfgas + 1; fgaspos(nfgas) = i
     ENDIF
  ENDDO

  ! Find indices of variables (other than trace gases and ozone)
  j = 1
  DO i = 1, n_fitvar_rad
     !WRITE(*, '(I5, A10, A20)') i, fitvar_rad_str(mask_fitvar_rad(i)), &
     !     fitvar_rad_unit(mask_fitvar_rad(i))
     IF (nfgas > 0) THEN
        IF (i >= fgasidxs(fgaspos(1)) .AND. i <= fgasidxs(fgaspos(nfgas))) CYCLE
     ENDIF
     IF (i >= ozfit_start_index .AND. i <= ozfit_end_index) CYCLE
     fothvarpos (j) = i; j = j + 1
  ENDDO

  ! -----------------------------------------------
  ! Close fitting control file, report SUCCESS read
  ! -----------------------------------------------
  CLOSE ( UNIT=fit_ctrl_unit )

  fitvar_rad_saved = fitvar_rad_init
  fitvar_rad_init_saved = fitvar_rad_init

  IF (scnwrt) WRITE(*, '(A, I8, A, I8)') ' n_fitvar_rad = ', n_fitvar_rad, '      nlayer = ', nlay

  RETURN
END SUBROUTINE read_ozprof_input
