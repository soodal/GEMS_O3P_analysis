  ! *********************** Modification History **********************
  ! Xiong Liu; July, 2003  (!xliu)
  ! 1. Read whether to do ozone profile retrieval and set a flag 
  ! 2. Read additional fitting control variables if ozprof_flag is set
  ! 3. Read an option about variable slit width, add yn_varyslit variable 
  !    n_slit_itnerval, slit_fname, slit_redo, wavcal_redo, wavcal_fname, 
  !    in USE OMSAO_variables_module
  ! 4. Add 1 nm more for winwav_min, winwav_max to avoid interpolation
  !    out of bounds for solar spectrum calibration
  ! 5. Read option use_meas_sig
  ! 6. Read option use_pixel_bin
  ! *******************************************************************

SUBROUTINE gems_o3p_read_ctr_file (fit_ctrl_unit, fit_ctrl_file,  pge_error_status)     

  ! ***********************************************************
  !
  !   Read fitting control parameters from input control file
  !
  ! ***********************************************************

  USE OMSAO_precision_module
  USE OMSAO_indices_module, ONLY: &
       max_rs_idx, calfit_strings, max_calfit_idx, radfit_strings, &
       mns_idx, mxs_idx, ad1_idx, bro_idx, lbe_idx, n_max_fitpars, refspec_strings, &
       genline_str, socline_str, racline_str, rafline_str, rspline_str,  &
       molline_str, eoi3str, solar_idx, amf_idx, us1_idx, us2_idx, &
       shift_offset, comm_idx, com1_idx, comfidx, cm1fidx, comvidx, cm1vidx
  USE OMSAO_parameters_module,   ONLY: maxchlen, maxwin, max_fit_pts,max_mol_fit 
  USE OMSAO_variables_module,    ONLY: use_backup, use_solcomp, avg_solcomp, avgsol_allorb,use_meas_sig, &
       fitcol_idx, fincol_idx, n_fincol_idx, n_mol_fit, &
       weight_sun, max_itnum_sol,weight_rad,max_itnum_rad,renorm, szamax, zatmos,  &
       which_slit, slit_trunc_limit, yn_varyslit, wavcal, wavcal_sol,smooth_slit,slit_rad, &
       slit_fit_pts,n_slit_step,slit_redo,wavcal_fit_pts, n_wavcal_step, wavcal_redo,   &
       yn_smooth, yn_doas,pm_one,  tol,  epsrel,  epsabs,  epsx,radwavcal_freq, phase, &
       n_fitvar_sol, fitvar_sol_init, fitvar_sol_saved, mask_fitvar_rad,fitvar_rad_str,      &
       n_fitvar_rad, fitvar_rad_init, fitvar_rad_saved, mask_fitvar_sol,fitvar_rad_unit,rmask_fitvar_rad,&  
       lo_sunbnd, up_sunbnd, lo_sunbnd_init,up_sunbnd_init, lo_radbnd, up_radbnd, &          
       static_input_fnames,refspec_fname, &     
       outdir, atmdbdir, refdbdir, &
       reduce_resolution,reduce_slit, redsampr, redlam, redfixwav,use_redfixwav,redfixwav_fname,nredfixwav, &
       rm_mgline,numwin, do_bandavg,wcal_bef_coadd,band_selectors,  n_band_avg, n_band_samp,& 
       winlim,winwav_min, winwav_max, &   
       have_amftable,have_undersampling,   &
       scnwrt, database_indices, radnhtrunc, refnhextra, & 
       pixnum_lim, linenum_lim,  l1b_rad_filename,l1b_irrad_filename,l2_cld_filename, l2_filename, geo_nc_filename
       
  USE OMSAO_errstat_module
  USE GEMS_O3P_gemsdata_module, ONLY: mchannel, nchannel,  ncoadd, gems_redslw, &
                                      upper_wvls, lower_wvls, retlbnd, retubnd ,&
				       do_xbin, do_ybin, nxbin, nybin

  USE ozprof_data_module,   ONLY: ozprof_str, ozprof_flag, ozprof_input_fname,         &
                   fullorb, do_ch2reso, l1l2inp_unit, nos, nsh, nsl, do_simu, radcalwrt

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER,           INTENT (IN)     :: fit_ctrl_unit
  CHARACTER (LEN=*), INTENT (IN)     :: fit_ctrl_file
  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER,           INTENT (OUT)   :: pge_error_status
  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER                  :: i, j, k, file_read_stat, sidx, ridx, idx, nidx, ntsh
  CHARACTER (LEN=maxchlen) :: tmpchar, l1l2_file
  CHARACTER (LEN=3)        :: idxchar, xbinchar, ybinchar
  CHARACTER (LEN=5)        :: idxchar1
  LOGICAL                  :: yn_eoi, l1l2_here
  REAL      (KIND=dp)      :: vartmp, lotmp, uptmp

  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=30), PARAMETER :: modulename = 'read_fitting_control_file'

  ! ========================
  ! Error handling variables
  ! ========================
  INTEGER :: errstat, version, ios

  CHARACTER (LEN=25), PARAMETER :: lm_instrument = 'Satellite instrument name'
  CHARACTER (LEN=20), PARAMETER :: lm_l1l2inputs   = 'Level 1/2 input files'
  CHARACTER (LEN=30), PARAMETER :: lm_atmdb      = 'Atmospheric database directory'
  CHARACTER (LEN=27), PARAMETER :: lm_refdb      = 'Reference spectra directory'
  CHARACTER (LEN=16), PARAMETER :: lm_outdb      = 'Output directory'
  CHARACTER (LEN=30), PARAMETER :: lm_bandselect = 'GEMS radiance bands to be used'
  CHARACTER (LEN=30), PARAMETER :: lm_coadding   = 'GEMS spatial coadding'
  CHARACTER (LEN=27), PARAMETER :: lm_reduceres  = 'Reduce spectral resolution'

  ! =================================
  ! External OMI and Toolkit routines
  ! =================================
  INTEGER :: OMI_SMF_setmsg

  pge_error_status = pge_errstat_ok
  fullorb=.TRUE.
  do_ch2reso=.FALSE.
  have_amftable = .FALSE.
  ! -----------------------------------------------------------
  ! Initialize array with reference spectrum names to Zero_Spec
  ! -----------------------------------------------------------
  DO j = solar_idx, max_rs_idx
     refspec_fname(j) = 'OMSAO_Zero_Spec.dat'
  END DO

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


  !xliu, 09/23/05 Add direcotry, remove hard code directory
  ! ----------------------------------------------------------
  ! Position cursor to read database directory
  ! ----------------------------------------------------------  
  REWIND ( fit_ctrl_unit )
  CALL skip_to_filemark ( fit_ctrl_unit, lm_atmdb, tmpchar, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, lm_atmdb, modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  ELSE
     READ (fit_ctrl_unit, '(A)') atmdbdir
  ENDIF

  REWIND ( fit_ctrl_unit)
  CALL skip_to_filemark ( fit_ctrl_unit, lm_refdb, tmpchar, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, lm_refdb, modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  ELSE
     READ (fit_ctrl_unit, '(A)') refdbdir
  ENDIF
 
  REWIND ( fit_ctrl_unit)
  CALL skip_to_filemark ( fit_ctrl_unit, lm_outdb, tmpchar, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file,lm_outdb, modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  ELSE
     READ (fit_ctrl_unit, '(A)') outdir
  ENDIF

  !xliu: add the following block 
  ! ----------------------------------------------------------
  ! Position cursor to read whether to retrieve ozone profile
  ! ----------------------------------------------------------
  REWIND ( fit_ctrl_unit )
  CALL skip_to_filemark ( fit_ctrl_unit, ozprof_str, tmpchar, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) THEN
     ozprof_flag = .FALSE.
     WRITE(*, *) 'This algorithm is only for ozone profile retrieval!!!'
     pge_error_status = pge_errstat_error; RETURN
  ELSE 
     READ (fit_ctrl_unit, *) ozprof_flag
     IF (ozprof_flag)  THEN
        READ (fit_ctrl_unit, '(A)') ozprof_input_fname
     END IF
  END IF
  
  ! ----------------------------------------------------------      
  !xliu, 01/03/2007, read options to degrade spectral resolution
  ! ----------------------------------------------------------

  REWIND ( fit_ctrl_unit )
  CALL skip_to_filemark ( fit_ctrl_unit, lm_reduceres, tmpchar, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, lm_reduceres, modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  END IF
  READ (fit_ctrl_unit, *) reduce_resolution
  READ (fit_ctrl_unit, *) reduce_slit
  READ (fit_ctrl_unit, *) gems_redslw(1:mchannel)
  IF ( reduce_slit == 1 ) gems_redslw(1:mchannel) = gems_redslw(1:mchannel) / 1.66511  ! convert from FWHM to hw1e 
  READ (fit_ctrl_unit, *) use_redfixwav
  IF (.NOT. reduce_resolution) use_redfixwav = .FALSE.
  READ (fit_ctrl_unit, '(A)') redfixwav_fname
  IF (use_redfixwav) THEN
     OPEN(UNIT=l1l2inp_unit, FILE=TRIM(ADJUSTL(redfixwav_fname)), STATUS='OLD', IOSTAT=errstat)
     IF ( errstat /= pge_errstat_ok ) THEN
        errstat = OMI_SMF_setmsg (omsao_e_open_fitctrl_file, &
             TRIM(ADJUSTL(redfixwav_fname)), modulename, 0)
        pge_error_status = pge_errstat_error; RETURN
     ELSE
        READ(l1l2inp_unit, *) nredfixwav
        IF (nredfixwav > max_fit_pts) THEN
           pge_error_status = pge_errstat_error; RETURN
        ENDIF
        READ(l1l2inp_unit, *) (redfixwav(i), i = 1, nredfixwav)
        CLOSE(UNIT=l1l2inp_unit) 
     END IF
  ENDIF
  READ (fit_ctrl_unit, *) redsampr 
  READ (fit_ctrl_unit, *) redlam

  ! -----------------------------------------------------
  ! Position cursor to read OMI channels used for fitting
  ! -----------------------------------------------------

  REWIND ( fit_ctrl_unit )
  CALL skip_to_filemark ( fit_ctrl_unit, lm_bandselect, tmpchar, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, lm_bandselect, modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  END IF

  READ (fit_ctrl_unit, *) rm_mgline
  READ (fit_ctrl_unit, *) numwin, do_bandavg, wcal_bef_coadd !, winwav_min, winwav_max
  IF (reduce_resolution) THEN
     do_bandavg = .FALSE.; rm_mgline = .FALSE.
  ENDIF
  IF (numwin > maxwin .OR. numwin < 1) THEN
     WRITE(*, *) 'Number of windows exceeds maxwin or less than 1!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF

  retlbnd = 1000.0; retubnd = 0.0
  DO i = 1, numwin
     READ(fit_ctrl_unit, *) band_selectors(i), winlim(i, 1), winlim(i, 2), &
          n_band_avg(i), n_band_samp(i)
     
     IF ((band_selectors(i) < 0) .OR. (band_selectors(i) > mchannel)) THEN
        WRITE(*, *) 'No such bands exist !!!'
        pge_error_status = pge_errstat_error; RETURN
     ENDIF
    
     IF (numwin > 1) THEN
        IF (winlim(i, 1) < lower_wvls(band_selectors(i)) .OR. winlim(i, 2) > upper_wvls(band_selectors(i))) THEN
           WRITE(*, *) 'Specified fitting windows does not make sense!!!'
           pge_error_status = pge_errstat_error; RETURN
        ENDIF
     ELSE
        ! Allow 2 extra nm for radiance calibration
        IF (winlim(i, 1) < lower_wvls(band_selectors(i)) - 1.0 .OR. winlim(i, 2) > upper_wvls(band_selectors(i)) + 1.0) THEN
            print*,lower_wvls(band_selectors(i))-1.0,upper_wvls(band_selectors(i))+1.0
           WRITE(*, *) 'Specified fitting windows does not make sense!!!'
           pge_error_status = pge_errstat_error; RETURN
        ENDIF
     ENDIF

     IF (i > 1) THEN
        IF  (band_selectors(i) < band_selectors(i-1) .OR. winlim(i, 1) < winlim(i-1, 2))  THEN
           WRITE(*, *) 'Incorrect band selection (must be in increasing wavelength) !!!'
           pge_error_status = pge_errstat_error; RETURN
        ENDIF
     ENDIF

     IF (winlim(i, 1) < retlbnd(band_selectors(i))) retlbnd(band_selectors(i)) = winlim(i, 1)
     IF (winlim(i, 2) > retubnd(band_selectors(i))) retubnd(band_selectors(i)) = winlim(i, 2)
  END DO 
  IF (MAXVAL(n_band_avg(1:numwin)) == 1) do_bandavg = .FALSE.

  DO i = 1, mchannel
     IF (retlbnd(i) == 1000.0) retlbnd(i) = lower_wvls(i)
     IF (retubnd(i) == 0.0)    retubnd(i) = upper_wvls(i)
  ENDDO
  
  IF (do_bandavg) THEN
     IF (ANY(n_band_avg(1:numwin) < 1)) THEN
        WRITE(*, *) 'Number of points for averaging must >= 1!!!'
        pge_error_status = pge_errstat_error; RETURN
     ENDIF
     
     IF (ANY(n_band_samp(1:numwin) < 1)) THEN
        WRITE(*, *) 'Number of points for sampling must >= 1!!!'
        pge_error_status = pge_errstat_error; RETURN
     ENDIF
  ENDIF
  winwav_min = winlim(1, 1) - 10.0
  winwav_max = winlim(numwin, 2) + 10.0

  ! ------------------------------------------------
  ! Position cursor to read spatial coadding option parameter
  ! ------------------------------------------------  
  REWIND ( fit_ctrl_unit )
  CALL skip_to_filemark ( fit_ctrl_unit, lm_coadding, tmpchar, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, lm_coadding, modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  END IF
  READ (fit_ctrl_unit, *) do_xbin, nxbin
  READ (fit_ctrl_unit, *) do_ybin, nybin
  IF (nxbin == 1) do_xbin = .false.
  IF (nybin == 1) do_ybin = .false.
  READ(fit_ctrl_unit, *) pixnum_lim
  READ(fit_ctrl_unit, *) linenum_lim
  READ(fit_ctrl_unit, *) l1l2_here
  READ(fit_ctrl_unit, '(A)') l1l2_file
  IF (l1l2_here) THEN
      OPEN( UNIT=l1l2inp_unit, FILE=TRIM(ADJUSTL(l1l2_file)),STATUS='OLD',IOSTAT=errstat)
      IF (errstat /= pge_errstat_ok) THEN
        errstat = OMI_SMF_setmsg (omsao_e_open_fitctrl_file,TRIM(ADJUSTL(l1l2_file)) ,modulename, 0)
        pge_error_status = pge_errstat_error; RETURN
      ELSE
        READ(l1l2inp_unit, '(A)') l1b_irrad_filename
        READ(l1l2inp_unit, '(A)') l1b_rad_filename
        READ(l1l2inp_unit, '(A)') l2_cld_filename
        READ(l1l2inp_unit, '(A)') l2_filename
        READ(l1l2inp_unit, *)  linenum_lim, pixnum_lim
      ENDIF
  ENDIF
  ! ------------------------------------------------
  ! Position cursor to read general input parameters
  ! ------------------------------------------------
  REWIND (fit_ctrl_unit)
  CALL skip_to_filemark ( fit_ctrl_unit, genline_str, tmpchar, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, genline_str, modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  END IF

  READ (fit_ctrl_unit, *) use_backup, use_solcomp, avg_solcomp, avgsol_allorb
  READ (fit_ctrl_unit, *) which_slit
  IF (reduce_resolution .AND. which_slit /= 0 .AND. which_slit /= 3) THEN
     WRITE(*, *) 'Have to use consistent slit function!!!'
     IF (reduce_slit == 1) which_slit = 0
     IF (reduce_slit == 2) which_slit = 3
  ENDIF

  READ (fit_ctrl_unit, *) slit_trunc_limit
  READ (fit_ctrl_unit, *) yn_varyslit, wavcal, wavcal_sol, smooth_slit, slit_rad
  IF (reduce_resolution .AND. yn_varyslit) THEN
     yn_varyslit = .FALSE.
     WRITE(*, *) 'Could not use variable slit function when to reduce resolution!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF
  IF (.NOT. wavcal) wavcal_sol = .FALSE.
  READ (fit_ctrl_unit, *) slit_fit_pts, n_slit_step, slit_redo
  READ (fit_ctrl_unit, *) wavcal_fit_pts, n_wavcal_step, wavcal_redo
  READ (fit_ctrl_unit, *) yn_smooth
  READ (fit_ctrl_unit, *) yn_doas
  READ (fit_ctrl_unit, *) use_meas_sig 
  READ (fit_ctrl_unit, *) tol
  READ (fit_ctrl_unit, *) epsrel
  READ (fit_ctrl_unit, *) epsabs
  READ (fit_ctrl_unit, *) epsx

  ! ----------------------------------------------------------
  ! Position cursor to read solar calibration input parameters
  ! ----------------------------------------------------------
  REWIND (fit_ctrl_unit)
  CALL skip_to_filemark ( fit_ctrl_unit, socline_str, tmpchar, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, socline_str, modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  END IF

  ! ---------------------------------------------------
  ! First thing to read is the Solar Reference Spectrum
  ! ---------------------------------------------------
  READ (fit_ctrl_unit, '(A)') refspec_fname(solar_idx)  
  refspec_fname(solar_idx) = TRIM(ADJUSTL(refdbdir)) // refspec_fname(solar_idx)  
  READ (fit_ctrl_unit, *) weight_sun
  READ (fit_ctrl_unit, *) max_itnum_sol

  n_fitvar_sol = 0;  fitvar_sol_init = 0.0
  solpars: DO i = 1, max_calfit_idx

     READ (fit_ctrl_unit, *) idxchar, vartmp, lotmp, uptmp
     ! ---------------------------------------------------------
     ! Check for consitency of bounds and adjust where necessary
     ! ---------------------------------------------------------
     IF ( lotmp > vartmp .OR. uptmp < vartmp ) THEN
        lotmp = vartmp ; uptmp = vartmp
     END IF
     IF ( lotmp == uptmp .AND. lotmp /= vartmp ) THEN
        uptmp = vartmp ; lotmp = vartmp
     END IF

     IF ( idxchar == eoi3str ) EXIT solpars

     CALL string2index ( calfit_strings, max_calfit_idx, idxchar, sidx )
     IF ( sidx > 0 ) THEN
        fitvar_sol_init(sidx) = vartmp
        lo_sunbnd(sidx) = lotmp ; up_sunbnd(sidx) = uptmp
        IF ( lotmp < uptmp ) THEN
           n_fitvar_sol = n_fitvar_sol + 1
           mask_fitvar_sol(n_fitvar_sol) = i
          ! print *, 'fitsol:', n_fitvar_sol,i, idxchar, vartmp
        ENDIF
     END IF
  END DO solpars
  fitvar_sol_saved = fitvar_sol_init
  lo_sunbnd_init = lo_sunbnd; up_sunbnd_init = up_sunbnd

  ! -------------------------------------------------------------
  ! Position cursor to read radiance calibration input parameters
  ! -------------------------------------------------------------
  REWIND (fit_ctrl_unit)
  CALL skip_to_filemark ( fit_ctrl_unit, racline_str, tmpchar, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, racline_str, modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  END IF

  READ (fit_ctrl_unit, *) renorm
  READ (fit_ctrl_unit, *) weight_rad
  READ (fit_ctrl_unit, *) max_itnum_rad
  READ (fit_ctrl_unit, *) radwavcal_freq
  READ (fit_ctrl_unit, *) szamax
  READ (fit_ctrl_unit, *) zatmos
  READ (fit_ctrl_unit, *) phase

  fitvar_rad_init = 0.0
  fitvar_rad_str = '      '
  radpars: DO i = 1, max_calfit_idx

     READ (fit_ctrl_unit, *) idxchar, vartmp, lotmp, uptmp
     ! ---------------------------------------------------------
     ! Check for consitency of bounds and adjust where necessary
     ! ---------------------------------------------------------
     IF ( lotmp > vartmp .OR. uptmp < vartmp ) THEN
        lotmp = vartmp ; uptmp = vartmp
     END IF
     IF ( lotmp == uptmp .AND. lotmp /= vartmp ) THEN
        uptmp = vartmp ; lotmp = vartmp
     END IF

     IF ( idxchar == eoi3str ) EXIT radpars
     CALL string2index ( calfit_strings, max_calfit_idx, idxchar, sidx )
     IF ( sidx > 0 ) THEN
        fitvar_rad_init(sidx) = vartmp
        fitvar_rad_str (sidx) = TRIM(ADJUSTL(idxchar))
        lo_radbnd(sidx) = lotmp ; up_radbnd(sidx) = uptmp
     END IF
  END DO radpars
   
  ! ---------------------------------------------------------
  ! Position cursor to read radiance fitting input parameters
  ! ---------------------------------------------------------
  REWIND (fit_ctrl_unit)
  CALL skip_to_filemark ( fit_ctrl_unit, rafline_str, tmpchar, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, rafline_str, modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  END IF

  ! By default we set the undersampling spectrum to FALSE. Only if we select
  ! it to be included in the fitting does it become trues. This way we save
  ! computation time for cases where we don't include the undersampling.
  have_undersampling = .FALSE.

  ! -------------------------------------------------------------
  ! Now keep reading spectrum blocks until EOF. This, obviously, 
  ! has to be the last READ action performed from the input file.
  ! -------------------------------------------------------------
  comvidx = 0; cm1vidx = 0; comfidx = 0; cm1fidx = 0
  fitvar_rad_unit = 'NoUnits'
  getpars: DO j = 1, max_rs_idx

     ! Read the spectrum identification string (SIS)
     READ (UNIT=fit_ctrl_unit, FMT='(A)', IOSTAT=errstat) tmpchar
     IF ( errstat /= file_read_ok ) THEN
        errstat = OMI_SMF_setmsg ( &
             omsao_e_read_fitctrl_file, 'radiance fitting parameters', modulename, 0)
        pge_error_status = pge_errstat_error; RETURN
     END IF
     CALL check_for_endofinput ( TRIM(ADJUSTL(tmpchar)), yn_eoi )
     IF ( yn_eoi ) EXIT getpars

     ! Convert SIS to index
     CALL string2index ( refspec_strings, max_rs_idx, tmpchar, ridx )

     ! the name of the corresponding reference spectrum
     READ (UNIT=fit_ctrl_unit, FMT='(A)', IOSTAT=errstat) refspec_fname(ridx)
     refspec_fname(ridx) = TRIM(ADJUSTL(refdbdir)) // refspec_fname(ridx)  

     ! Read the block of fitting parameters for current reference spectrum
     DO k = 1, mxs_idx

        READ (fit_ctrl_unit, *) idxchar, vartmp, lotmp, uptmp
        ! ---------------------------------------------------------
        ! Check for consitency of bounds and adjust where necessary
        ! ---------------------------------------------------------
        IF ( lotmp > vartmp .OR. uptmp < vartmp ) THEN
           lotmp = vartmp ; uptmp = vartmp
        END IF
        IF ( lotmp == uptmp .AND. lotmp /= vartmp ) THEN
           uptmp = vartmp ; lotmp = vartmp
        END IF

        CALL string2index ( radfit_strings, mxs_idx, idxchar, sidx )
        IF ( sidx > 0 ) THEN
           i = max_calfit_idx + (ridx-1)*mxs_idx + sidx
           fitvar_rad_init(i) = vartmp
           fitvar_rad_str (i) = TRIM(ADJUSTL(tmpchar))
           lo_radbnd (i) = lotmp ; up_radbnd (i) = uptmp
              
           IF ( (ridx == us1_idx .OR. ridx == us2_idx) .AND. &
                ANY ( (/ vartmp,lotmp,uptmp /) /= 0.0 ) ) have_undersampling = .TRUE.

           IF ( ridx == comm_idx .AND. lotmp < uptmp ) THEN
              comvidx = i
           ENDIF

           IF ( ridx == com1_idx .AND. lotmp < uptmp ) THEN
              cm1vidx = i
           ENDIF
        END IF
     END DO

     ! read shift parameter
     READ (fit_ctrl_unit, *) idxchar1, vartmp, lotmp, uptmp
     IF ( lotmp > vartmp .OR. uptmp < vartmp ) THEN
        lotmp = vartmp ; uptmp = vartmp
     END IF
     IF ( lotmp == uptmp .AND. lotmp /= vartmp ) THEN
        uptmp = vartmp ; lotmp = vartmp
     END IF

     i =  max_calfit_idx + (ridx-1)*mxs_idx + 1
     IF  (ALL(lo_radbnd(i:i+2) - up_radbnd(i:i+2) >= 0.0)) THEN
        vartmp = 0.0; uptmp = 0.0; lotmp = 0.0
     END IF

     i =  shift_offset + ridx
     fitvar_rad_init(i) = vartmp
     fitvar_rad_str(i) = TRIM(ADJUSTL(idxchar1))
     fitvar_rad_unit(i) = 'nm'
     
     lo_radbnd (i) = lotmp ; up_radbnd (i) = uptmp        
  END DO getpars

  ! -----------------------------------------------------
  ! For safety, copy REFSPEC_FNAME to STATIC_INPUT_FNAMES
  ! (in the OMI branch this is the other way around since
  !  there we get file names from the PCF)
  ! -----------------------------------------------------
  DO j = solar_idx, max_rs_idx
     static_input_fnames(j) = refspec_fname(j)
  END DO

  ! -------------------------------------------------------------
  ! Find the indices of those variables that are actually varied
  ! during the fitting, and save those in MASK_FITVAR_RAD. Save
  ! the number of varied parameters in N_FITVAR_RAD.
  !
  ! In addition, we need to determine the fitting indices that
  ! will make up the final fitted column of the molecule(s) in
  ! question. For this we have to jump through a double loop:
  ! Since we are compressing the fitting parameter array to 
  ! include only the varied parameters, the final covariance 
  ! matrix, which is crucial for determining the uncertainties,
  ! only knows the compressed indices. Therefore we have to have
  ! an "index of an index" type array, that remembers the index
  ! position of indices of the final molecule(s), AS THEY APPEAR
  ! IN THE COMPRESSED FITTING PARAMETER LIST.
  !
  ! For the latter task it is easier to split the loops into
  ! calibration parameters (unrelated to the final fitted column)
  ! and reference spectra parameters.
  ! -------------------------------------------------------------
  n_fitvar_rad = 0 ; mask_fitvar_rad = 0
  rmask_fitvar_rad = 0; database_indices = 0
  ! --------------------------------
  ! First the calibration parameters
  ! --------------------------------
  DO i = 1, max_calfit_idx
     IF ( lo_radbnd(i) < up_radbnd(i) ) THEN
        n_fitvar_rad = n_fitvar_rad + 1
        mask_fitvar_rad(n_fitvar_rad) = i
        rmask_fitvar_rad(i) = n_fitvar_rad      
     END IF
  END DO

  ! ------------------------------------
  ! Now the reference spectra parameters
  ! ------------------------------------
  n_fincol_idx = 0
  DO i = 1, max_rs_idx
     idx = max_calfit_idx + (i-1) * mxs_idx
     DO j = mns_idx, mxs_idx        
        ! ----------------------------------------------
        ! Assign only entries that are varied in the fit
        ! ----------------------------------------------
        IF ( lo_radbnd(idx+j) < up_radbnd(idx+j) ) THEN
           n_fitvar_rad = n_fitvar_rad + 1
           mask_fitvar_rad(n_fitvar_rad) = idx+j
           rmask_fitvar_rad(idx+j) = n_fitvar_rad
           database_indices(n_fitvar_rad) = i

           ! -------------------------------------------------
           ! And here the loop over the final column molecules.
           ! We have to match the FITCOL_IDX with the current
           ! molecule index, and then remember the position of
           ! the fitting index in the MASK_FITVAR_RAD array.
           ! The second index remembers the reference spectrum
           ! that is associated with this molecule, so that we
           ! can easily access its normalization factor.
           ! -------------------------------------------------
           IF (.NOT. ozprof_flag) THEN  !xliu
              getfincol: DO k = 1, n_mol_fit
                 IF ( fitcol_idx(k) == i ) THEN
                    n_fincol_idx = n_fincol_idx + 1
                    fincol_idx (1,n_fincol_idx) = n_fitvar_rad
                    fincol_idx (2,n_fincol_idx) = i
                    EXIT getfincol
                 END IF
              END DO getfincol
           ENDIF   !xliu

           IF (idx + j == comvidx) comfidx = n_fitvar_rad
           IF (idx + j == cm1vidx) cm1fidx = n_fitvar_rad  


        END IF
     END DO
  END DO

  ! For shift
  ntsh = 0
  DO i = 1, max_rs_idx
     j = shift_offset + i
     IF ( lo_radbnd(j) < up_radbnd(j) ) THEN
        n_fitvar_rad = n_fitvar_rad + 1
        mask_fitvar_rad(n_fitvar_rad) = j
        rmask_fitvar_rad(j) = n_fitvar_rad
        ntsh = ntsh + 1
        print *,'shift' , n_fitvar_rad,idx+j,fitvar_rad_str(idx+j),fitvar_rad_init(idx+j), fitvar_rad_unit(idx+j)
     END IF
  END DO

  ! -----------------------------------------------
  ! Close fitting control file, report SUCCESS read
  ! -----------------------------------------------
  CLOSE ( UNIT=fit_ctrl_unit )

  !xliu: add the following block
  ! ------------------------------------------------------------------------------
  ! Read fitting conrol parameters from input file for ozone profile variables
  IF (ozprof_flag) THEN 
     CALL read_ozprof_input ( &
          fit_ctrl_unit, ozprof_input_fname, pge_error_status )
     IF ( pge_error_status >= pge_errstat_error ) RETURN 
  END IF

  ! refnhextra must >= 1 and radnhtrunc > refnhextra, if interpolation is performed
  ! radnhtrunc should be 
  IF ( (ntsh == 0 .AND. nsh == 0 .AND. nos == 0 .AND. nsl == 0) &
       .OR. (do_simu .AND. .NOT. radcalwrt)) THEN
     radnhtrunc = 3; refnhextra = 2
  !ELSE IF (reduce_resolution .AND. use_redfixwav) THEN
  !   radnhtrunc = 2; refnhextra = 1
  ELSE
     radnhtrunc = 2; refnhextra = 1
  ENDIF
     

  ! ------------------------------------------------------------------------------
  fitvar_rad_saved = fitvar_rad_init

  IF ( yn_doas ) THEN
     pm_one     = -1.D0
  ELSE
     pm_one     = 1.D0
  END IF

  RETURN
END SUBROUTINE gems_o3p_read_ctr_file


