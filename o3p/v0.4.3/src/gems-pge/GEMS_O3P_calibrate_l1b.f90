SUBROUTINE gems_o3p_rad_cross_calibrate (first_pix, last_pix,ny,offline, pge_error_status)
  USE OMSAO_precision_module
  USE OMSAO_indices_module,         ONLY: wvl_idx, spc_idx, sig_idx, max_calfit_idx, hwe_idx
  USE OMSAO_parameters_module,      ONLY: downweight, normweight, max_fit_pts
  USE OMSAO_variables_module,       ONLY: currpix, currpixchar,currline, currtrack, n_rad_wvl, curr_rad_spec,  &
                                          numwin, nradpix ,scnwrt,use_meas_sig, yn_varyslit, which_slit,&
                                          wincal_wav,wavcal,    nslit_rad,  radslitfit,slitwav_rad,radwinfit,     radwavfit,slit_rad,     &
                                          wavcal_sol, nslit,      slitfit,slitwav ,  solwinfit,     nwavcal_rad, sswav_rad, slitdis
          
  USE ozprof_data_module,          ONLY: calunit
  USE OMSAO_errstat_module
  USE GEMS_O3P_gemsdata_module,    only:gems_rad, & !INPUT
                                        gems_solcal, gems_radcal, gems_wincal_wav !OUTPUT
  IMPLICIT NONE

  ! ================
  ! Output variables
  ! ================
  INTEGER, INTENT (IN)        :: first_pix, last_pix, ny, offline
  INTEGER, INTENT (OUT)       :: pge_error_status

  ! ===============
  ! Local variables
  ! ===============
  REAL (KIND=dp), DIMENSION(max_fit_pts) :: diff, corr, wav
  REAL (KIND=dp)                         :: dwvl
  INTEGER  :: i, j, ix, iline, mline, inter_errstat, np, fidx, lidx
  LOGICAL  :: error
  CHARACTER (LEN=23), PARAMETER :: modulename = 'gems_o3p_rad_cross_calibrate'

  pge_error_status = pge_errstat_ok; error   = .FALSE.
  mline = NINT(ny / 2.0) - 1

  DO ix = first_pix, last_pix
     currpix = ix;    WRITE(currpixchar, '(A1, I2.2)') 'x', ix
     np = MAXVAL(gems_rad%nwav(ix, 1:ny))
     IF ( ALL( gems_rad%errstat(ix, 1:ny) == pge_errstat_error )) CYCLE
     ! Only process a scan line in the middle 
     DO i = mline, 1, -1  ! first half      
        IF (gems_rad%errstat(ix, i) == pge_errstat_ok .AND. gems_rad%nwav(ix, i) == np) EXIT
     ENDDO

     IF ( i < 0 ) THEN    ! second half
        DO i = mline + 1, ny ! need more consideration
           IF (gems_rad%errstat(ix, i) == pge_errstat_ok .AND. gems_rad%nwav(ix, i) == np) EXIT
        ENDDO 
     ENDIF
  
     IF (i > ny ) CYCLE  ! no valid pixels for this xtrack position
     iline = i
 
    ! Spectrum
     n_rad_wvl = gems_rad%nwav(ix, iline)  
     curr_rad_spec(wvl_idx, 1:n_rad_wvl) = gems_rad%wavl(1:n_rad_wvl, ix, iline)
     curr_rad_spec(spc_idx, 1:n_rad_wvl) = gems_rad%spec(1:n_rad_wvl, ix, iline)
   
     IF (use_meas_sig) THEN
        curr_rad_spec(sig_idx, 1:n_rad_wvl) = gems_rad%prec(1:n_rad_wvl, ix, iline)
     ELSE
        curr_rad_spec(sig_idx, 1:n_rad_wvl) = normweight
     ENDIF

     nradpix(1:numwin) = gems_rad%npix(1:numwin, ix, iline)   
     ! Perform calibration
     IF (scnwrt) WRITE(*, '(A, 3I5)') 'Performing radiance wavelength calibration: ', currline, currtrack
     IF (which_slit < 4) THEN
        IF (yn_varyslit .AND. slit_rad ) THEN
           CALL rad_fit_vary(calunit, n_rad_wvl, curr_rad_spec(wvl_idx:sig_idx, 1:n_rad_wvl), error)
           gems_radcal%nslit (ix) = nslit_rad
           gems_radcal%slitfit(1:nslit_rad, 1:max_calfit_idx, :, ix) = radslitfit(1:nslit_rad, 1:max_calfit_idx, :)
           gems_radcal%slitwav(1:nslit_rad, ix) = slitwav_rad(1:nslit_rad)         
        ELSE IF (yn_varyslit ) THEN
           nslit            = gems_solcal%nslit(ix)
           slitwav(1:nslit) = gems_solcal%slitwav(1:nslit, ix)
           slitfit          = gems_solcal%slitfit(:, :, :, ix)  !max_fit_pts, max_calfit_idx, 2, nxtrack_max
        ELSE 
           solwinfit = gems_solcal%winfit(:, :, :, ix)          !maxwin,       max_calfit_idx,2, nxtrack_max
        ENDIF
     ENDIF 

     IF (wavcal) THEN
        IF (yn_varyslit .AND. (wavcal_sol .OR. which_slit == 4) ) THEN
           CALL rad_wavcal_vary (calunit, n_rad_wvl, curr_rad_spec(wvl_idx:sig_idx, 1:n_rad_wvl), error)
           gems_radcal%nwavcal(ix) = nwavcal_rad
           gems_radcal%wavfit(1:nwavcal_rad, 1:max_calfit_idx, :, ix) = radwavfit(1:nwavcal_rad, 1:max_calfit_idx, :)
           gems_radcal%sswav(1:nwavcal_rad, ix) = sswav_rad(1:nwavcal_rad)
        ELSE       
           CALL radiance_wavcal (n_rad_wvl, curr_rad_spec(wvl_idx:sig_idx, 1:n_rad_wvl), error)
           gems_wincal_wav(1:numwin, ix) = wincal_wav(1:numwin)
           gems_radcal%winfit(1:numwin, 1:max_calfit_idx, :, ix) = radwinfit(1:numwin, 1:max_calfit_idx, :)
        END IF       
     ENDIF

     IF (error) THEN
        pge_error_status = pge_errstat_warning
     ENDIF
     
     !----------------------------
     ! Save calibrated wavelengths
     !---------------------------- 
     diff(1:n_rad_wvl)  = curr_rad_spec(wvl_idx, 1:n_rad_wvl) - gems_rad%wavl(1:n_rad_wvl, ix, iline)

     ! Check wavelength shifts, if it is too large, do not apply
     fidx = 1
     DO i = 1, numwin
        lidx = fidx + nradpix(i) - 1
        dwvl =  curr_rad_spec(wvl_idx, fidx+1) - curr_rad_spec(wvl_idx, fidx)
        IF (ANY(ABS(diff(fidx:lidx)) > dwvl)) diff(fidx:lidx) = 0.0
        fidx = lidx + 1
     ENDDO
     IF (MAXVAL(ABS(diff(1:n_rad_wvl))) == 0.0) RETURN

     wav (1:n_rad_wvl)  = gems_rad%wavl(1:n_rad_wvl, ix, iline)
     DO i = 1, ny
        IF (gems_rad%errstat(ix, i) == pge_errstat_ok) THEN
           np = gems_rad%nwav(ix, i)

           IF (np == n_rad_wvl) THEN
             gems_rad%wavl(1:np, ix, i) = gems_rad%wavl(1:np, ix, i) + diff(1:np)
              IF ( ANY(gems_rad%wavl(2:np, ix, i) - gems_rad%wavl(1:np-1, ix, i) < 0.0) ) THEN
                gems_rad%wavl(1:np, ix, i) = gems_rad%wavl(1:np, ix, i) - diff(1:np)
                gems_rad%errstat(ix, i) = pge_errstat_warning
              ENDIF
           ELSE  ! Perform interpolation
              CALL INTERPOL(wav(1:n_rad_wvl), diff(1:n_rad_wvl), n_rad_wvl,  &
                   gems_rad%wavl(1:np, ix, i), corr(1:np), np, inter_errstat)
              IF (inter_errstat < 0)  THEN
                 gems_rad%errstat(ix, i) = pge_errstat_warning ! wavelength not calibrated
              ELSE
                 gems_rad%wavl(1:np, ix, i) = gems_rad%wavl(1:np, ix, i) + corr(1:np)
                 IF ( ANY(gems_rad%wavl(2:np, ix, i) -gems_rad%wavl(1:np-1, ix, i) < 0.0) ) THEN
                    gems_rad%wavl(1:np, ix, i) = gems_rad%wavl(1:np, ix, i) - corr(1:np)
                    gems_rad%errstat(ix, i) = pge_errstat_warning
                 ENDIF
              ENDIF
           ENDIF
        ENDIF

     ENDDO
  ENDDO

  RETURN
END SUBROUTINE gems_o3p_rad_cross_calibrate

SUBROUTINE GEMS_O3P_SUB3_Irrad_cross_calibrate (first_pix, last_pix,  pge_error_status)

  USE OMSAO_precision_module
  USE OMSAO_indices_module,         ONLY: wvl_idx, spc_idx, sig_idx, max_calfit_idx, hwe_idx, vgr_idx, shi_idx
  USE OMSAO_parameters_module,      ONLY: downweight, normweight
  USE OMSAO_variables_module,       ONLY: currpix, currpixchar, band_selectors, numwin,  &
                                          wavcal, wavcal_sol, reduce_resolution, yn_varyslit,scnwrt, &
                                          which_slit, use_meas_sig, nsolpix,n_irrad_wvl, curr_sol_spec, &
                                          nsol_ring, sol_spec_ring, &
                                          wincal_wav, solwinfit, &
                                          nslit_sol, solslitfit, slitwav_sol, &
                                          nwavcal_sol, solwavfit, sswav_sol
  USE ozprof_data_module,          ONLY: calunit
  USE OMSAO_errstat_module
  USE GEMS_O3P_gemsdata_module,    ONLY: nxbin, gems_irrad, gems_ring, &  !INPUT
                                         gems_wincal_wav,gems_solcal !output  

  IMPLICIT NONE

  ! ================
  ! Output variables
  ! ================
  INTEGER, INTENT (IN)        :: first_pix, last_pix
  INTEGER, INTENT (OUT)       :: pge_error_status

  ! ===============
  ! Local variables
  ! ===============
  INTEGER                     :: i, ix, fidx, lidx
  REAL (KIND=dp)              :: dwvl
  LOGICAL                     :: error
  CHARACTER (LEN=23), PARAMETER :: modulename = 'gems_o3p_irrad_cross_calibrate'

  pge_error_status = pge_errstat_ok

  DO ix = first_pix, last_pix

     ! Bad solar, no calibration
     IF (gems_irrad%errstat(ix) == pge_errstat_error) CYCLE
     error   = .FALSE.

     currpix = ix
     WRITE(currpixchar, '(A1, I2.2)') 'x', ix

     ! Spectrum
     n_irrad_wvl = gems_irrad%nwav(ix)
     curr_sol_spec(wvl_idx,1:n_irrad_wvl) = gems_irrad%wavl(1:n_irrad_wvl, ix)
     curr_sol_spec(spc_idx,1:n_irrad_wvl) = gems_irrad%spec(1:n_irrad_wvl, ix)

     IF (use_meas_sig) THEN
        curr_sol_spec(sig_idx, 1:n_irrad_wvl) = gems_irrad%prec(1:n_irrad_wvl, ix)
     ELSE
        curr_sol_spec(sig_idx, 1:n_irrad_wvl) = normweight
     ENDIF
     nsolpix(1:numwin) = gems_irrad%npix(1:numwin, ix)


     ! Ring Spectrum
     nsol_ring = gems_ring%sol_uin(ix) - gems_ring%sol_lin(ix) + 1
     sol_spec_ring(1, 1:nsol_ring) = gems_ring%solwavl(gems_ring%sol_lin(ix):gems_ring%sol_uin(ix), ix)

     ! Perform calibration
     IF (scnwrt) WRITE(*, '(A,I5)') 'Performing solar wavelength calibration: ', (ix - 1) * nxbin + 1

     IF (yn_varyslit) THEN !F
        IF (which_slit <= 3 ) THEN
           CALL solar_fit_vary (calunit, error )
           gems_solcal%nslit(ix) = nslit_sol
           gems_solcal%slitfit(1:nslit_sol, 1:max_calfit_idx, :, ix) = solslitfit(1:nslit_sol, 1:max_calfit_idx, :)
           gems_solcal%slitwav(1:nslit_sol, ix) = slitwav_sol(1:nslit_sol)
        ENDIF

        IF (wavcal .AND. (wavcal_sol .OR. which_slit == 4)) THEN
           CALL solar_wavcal_vary(calunit, error)
           gems_solcal%nwavcal(ix) = nwavcal_sol
           gems_solcal%wavfit(1:nwavcal_sol, 1:max_calfit_idx, :, ix) = solwavfit(1:nwavcal_sol, 1:max_calfit_idx, :)
           gems_solcal%sswav(1:nwavcal_sol, ix) = sswav_sol(1:nwavcal_sol)
        ENDIF
     ELSE ! yn_varyslit is False
        IF (wavcal .OR. which_slit /= 4) THEN ! F, 0
           CALL solar_fit ( error )
           gems_wincal_wav(1:numwin, ix) = wincal_wav(1:numwin) ! averaged lamda over each fitting window

           IF (reduce_resolution) THEN ! F
              WRITE(*,*) 'No process for reduce_resolution in ', modulename
           ENDIF
           gems_solcal%winfit(1:numwin, 1:max_calfit_idx, :, ix) = solwinfit(1:numwin, 1:max_calfit_idx, :)

        ENDIF
     END IF

     ! Check wavelength shifts, if it is too large, do not apply
     fidx = 1
     DO i = 1, numwin
        lidx = fidx + nsolpix(i) - 1
        dwvl =  curr_sol_spec(wvl_idx, fidx+1) - curr_sol_spec(wvl_idx, fidx)

        IF ( ANY(ABS(curr_sol_spec(wvl_idx, fidx:lidx) - gems_irrad%wavl(fidx:lidx, ix)) > dwvl) .OR. &
             ANY(curr_sol_spec(wvl_idx, fidx+1:lidx) - curr_sol_spec(wvl_idx, fidx:lidx-1) < 0.0) ) THEN
           error = .TRUE.

        ELSE
           gems_irrad%wavl(fidx:lidx, ix) = curr_sol_spec(wvl_idx, fidx:lidx)
        ENDIF

        fidx = lidx + 1
     ENDDO

     IF (error) THEN
        gems_irrad%errstat(ix) = pge_errstat_warning; CYCLE
     ENDIF


     gems_ring%solwavl(gems_ring%sol_lin(ix):gems_ring%sol_uin(ix), ix)  = sol_spec_ring(1, 1:nsol_ring)
  ENDDO

  RETURN
END SUBROUTINE GEMS_O3P_SUB3_Irrad_cross_calibrate



SUBROUTINE solwavcal_coadd(wcal_bef_coadd, nspec, ncoadd, allspec, wshis, wsqus, error)
  
  USE OMSAO_precision_module
  USE OMSAO_indices_module,     ONLY: max_calfit_idx, shi_idx, squ_idx, &
       wvl_idx, spc_idx, sig_idx, hwe_idx, asy_idx, hwr_idx, vgl_idx
  USE OMSAO_variables_module,   ONLY: n_fitvar_sol, fitvar_sol,   &
       mask_fitvar_sol, lo_sunbnd, up_sunbnd, fixslitcal, fitwavs, &
       fitweights, currspec, fitvar_sol_saved, which_slit
  USE OMSAO_errstat_module
       
  IMPLICIT NONE

  ! ================================
  ! Input and Output variables
  ! =================================
  INTEGER, INTENT(IN)  :: nspec, ncoadd
  LOGICAL, INTENT(IN)  :: wcal_bef_coadd
  REAL (KIND=dp), DIMENSION(ncoadd, sig_idx, nspec), INTENT(INOUT) :: allspec
  REAL (KIND=dp), DIMENSION(ncoadd),                 INTENT(OUT) :: wshis, wsqus
  LOGICAL,                                           INTENT(OUT) :: error

  ! ===============
  ! Local variables
  ! ===============
  REAL (KIND=dp), DIMENSION(max_calfit_idx, 2) :: tmp_fitvar
  REAL (KIND=dp)                               :: tmpwave, solar_norm, dwvl         
  INTEGER       :: i, j, solfit_exval, errstat
  INTEGER, SAVE :: slit_unit
  LOGICAL, SAVE :: wrt_to_screen, wrt_to_file, slitcal, first = .TRUE.

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=15), PARAMETER :: modulename = 'solwavcal_coadd'

  ! ------------------
  ! External functions
  ! ------------------
  INTEGER :: OMI_SMF_setmsg

  error = .FALSE.
  wshis  = 0.0; wsqus = 0.0

  IF (wcal_bef_coadd) THEN

     ! find the locations of actually used fitting variables
     IF (first) THEN
        fixslitcal = .TRUE.    ; slitcal = .TRUE.
        wrt_to_screen = .FALSE.; wrt_to_file = .FALSE. ; slit_unit = 1000

        IF (which_slit == 4) THEN
           fitvar_sol(hwe_idx:asy_idx) = 0_dp
           lo_sunbnd(hwe_idx:asy_idx)  = 0_dp
           up_sunbnd(hwe_idx:asy_idx)  = 0_dp
           
           fitvar_sol(vgl_idx:hwr_idx) = 0_dp
           lo_sunbnd(vgl_idx:hwr_idx)  = 0_dp
           up_sunbnd(vgl_idx:hwr_idx)  = 0_dp
           slitcal = .FALSE.; fixslitcal = .FALSE.
        ENDIF
        
        n_fitvar_sol = 0
        DO i = 1, max_calfit_idx
           IF (lo_sunbnd(i) < up_sunbnd(i) ) THEN
              n_fitvar_sol =  n_fitvar_sol + 1
              mask_fitvar_sol(n_fitvar_sol) = i
           END IF
        ENDDO
        first = .FALSE.
     ENDIF
     
     solar_norm = SUM(allspec(1, spc_idx, 1:nspec)) / nspec
     DO i = 1, ncoadd
        fitwavs   (1:nspec)  = allspec(i, wvl_idx, 1:nspec)
        currspec  (1:nspec)  = allspec(i, spc_idx, 1:nspec) / solar_norm
        fitweights(1:nspec)  = allspec(i, sig_idx, 1:nspec) / solar_norm
        
        fitvar_sol = fitvar_sol_saved
        CALL cal_fit_one (nspec, n_fitvar_sol, wrt_to_screen, wrt_to_file,&
             slitcal, slit_unit, tmpwave, tmp_fitvar, solfit_exval)
        
        IF (solfit_exval < 0) THEN
           WRITE(*, *) 'solwavcal_coadd: calibration not converge for pixel: ', i
           error = .TRUE.; wshis(i) = 0.; wsqus(i) = 0.
        END IF

        ! Shift and squeeze earthshine spectrum
        allspec(i, wvl_idx, 1:nspec) = (fitwavs(1:nspec) - fitvar_sol(shi_idx)) / (1.0 + fitvar_sol(squ_idx)) 

        ! Make sure that the correction is less than a pixel
        dwvl = fitwavs(2) - fitwavs(1)
        IF ( ANY( ABS(allspec(i, wvl_idx, 1:nspec) - fitwavs(1:nspec)) > dwvl ) ) THEN
           allspec(i, wvl_idx, 1:nspec) = fitwavs(1:nspec) ! Roll back
           wshis(i) = 0.0; wsqus(i) = 0.0
        ELSE             
           wshis(i) = fitvar_sol(shi_idx); wsqus(i) = fitvar_sol(squ_idx)
        ENDIF
     ENDDO
  ENDIF

  IF (.NOT. wcal_bef_coadd) THEN  ! simple coadding due to unsuccessful calibration
     DO i = 2, ncoadd
        allspec(1, :, :)    = allspec(1, :, :) + allspec(i, :, :)
     ENDDO
     allspec(1, wvl_idx, :) = allspec(1, wvl_idx, :) / ncoadd  
     allspec(1, spc_idx, :) = allspec(1, spc_idx, :) / ncoadd  ! reduce S/N
     allspec(1, sig_idx, :) = allspec(1, sig_idx, :) / ncoadd  / SQRT(1.0 * ncoadd)  ! reduce S/N
  ELSE     
     ! Interpolate every other spectrum to the wavelength grids of first spectrum
     DO i = 2, ncoadd
        CALL interpolation (nspec, allspec(i, wvl_idx, 1:nspec), allspec(i, spc_idx, 1:nspec), &
             nspec - 2,  allspec(1, wvl_idx, 2:nspec-1), allspec(i, spc_idx, 2:nspec-1), errstat )
        IF ( errstat > pge_errstat_warning ) THEN
           errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; STOP 1
        END IF
        
        allspec(1, wvl_idx, 1) = allspec(1, wvl_idx, 1) + allspec(i, wvl_idx, 1)
        allspec(1, wvl_idx, nspec) = allspec(1, wvl_idx, nspec) + allspec(i, wvl_idx, nspec)
        allspec(1, spc_idx, :) = allspec(1, spc_idx, :) + allspec(i, spc_idx, :)
        allspec(1, sig_idx, :) = allspec(1, sig_idx, :) + allspec(i, sig_idx, :)
     ENDDO
     allspec(1, wvl_idx, 1) = allspec(1, wvl_idx, 1) / ncoadd  
     allspec(1, wvl_idx, nspec) = allspec(1, wvl_idx, nspec) / ncoadd  
     allspec(1, spc_idx, :) = allspec(1, spc_idx, :) / ncoadd  ! reduce S/N
     allspec(1, sig_idx, :) = allspec(1, sig_idx, :) / ncoadd  / SQRT(1.0 * ncoadd)  ! reduce S/N    
  ENDIF
  
  RETURN
  
END SUBROUTINE solwavcal_coadd


SUBROUTINE radwavcal_coadd(wcal_bef_coadd, wavcal, iwin, ix, nspec, ncoadd, allspec, wshis, wsqus, error)
  
  USE OMSAO_precision_module
  USE OMSAO_indices_module,     ONLY: max_calfit_idx, shi_idx, squ_idx, &
       wvl_idx, spc_idx, sig_idx, hwe_idx, asy_idx, vgl_idx, hwr_idx 
  USE OMSAO_variables_module,   ONLY: n_fitvar_sol, fitvar_sol,   &
       mask_fitvar_sol, lo_sunbnd, up_sunbnd, fixslitcal, fitwavs, &
       fitweights, currspec, fitvar_sol_saved, which_slit, fitvar_sol_init, &
       yn_varyslit, nslit, slitwav, slitfit, lo_sunbnd_init, up_sunbnd_init
  USE OMSAO_errstat_module
  
  IMPLICIT NONE
  
  ! ================================
  ! Input and Output variables
  ! =================================
  INTEGER, INTENT(IN) :: nspec, ncoadd, iwin, ix
  REAL (KIND=dp), DIMENSION(ncoadd, sig_idx, nspec), INTENT(INOUT)   :: allspec
  REAL (KIND=dp), DIMENSION(ncoadd),                   INTENT(INOUT) :: wshis, wsqus
  LOGICAL,                                             INTENT( IN)   :: wcal_bef_coadd, wavcal
  LOGICAL,                                             INTENT(OUT)   :: error

  ! ===============
  ! Local variables
  ! ===============
  REAL (KIND=dp), DIMENSION(max_calfit_idx, 2) :: tmp_fitvar
  REAL (KIND=dp)                               :: tmpwave, solar_norm, dwvl         
  INTEGER         :: i, n_fit_pts, solfit_exval, errstat
  LOGICAL, SAVE   :: wrt_to_screen, wrt_to_file, slitcal, first = .TRUE.
  INTEGER, SAVE   :: slit_unit

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=14), PARAMETER :: modulename = 'radwavcal_coadd'

  ! ------------------
  ! External functions
  ! ------------------
  INTEGER :: OMI_SMF_setmsg

  error = .FALSE.
  IF (wcal_bef_coadd .AND. wavcal) THEN

     ! find the locations of actually used fitting variables
     !IF (first) THEN
        fixslitcal  = .TRUE. ; slitcal = .TRUE.;  wrt_to_screen = .FALSE.
        wrt_to_file = .FALSE.; slit_unit = 1000

        fitvar_sol = fitvar_sol_init
        lo_sunbnd  = lo_sunbnd_init
        up_sunbnd  = up_sunbnd_init

        IF (which_slit == 4) THEN
           ! fix variables for hwe, asy, vgr, vgl, hwr, hwl
           fitvar_sol(hwe_idx:asy_idx) = 0.0
           lo_sunbnd(hwe_idx:asy_idx) = 0.0
           up_sunbnd(hwe_idx:asy_idx) = 0.0
           fitvar_sol(vgl_idx:hwr_idx) = 0.0
           lo_sunbnd(vgl_idx:hwr_idx) = 0.0
           up_sunbnd(vgl_idx:hwr_idx) = 0.0
        ENDIF
        
        n_fitvar_sol = 0
        DO i = 1, max_calfit_idx
           IF (lo_sunbnd(i) < up_sunbnd(i) ) THEN
              n_fitvar_sol  =  n_fitvar_sol + 1
              mask_fitvar_sol(n_fitvar_sol) = i
           END IF
        ENDDO
     !   first = .FALSE.
     !ENDIF

     !IF (which_slit /= 4) THEN
     !   IF ( .NOT. yn_varyslit)  THEN
     !      fitvar_sol(hwe_idx:hwr_idx) = omi_solwinfit(iwin, hwe_idx:hwr_idx, 1, ix)
     !      print *, fitvar_sol(hwe_idx : hwr_idx)
     !   ELSE
     !      nslit = omi_nslit_sol(ix)
     !      slitwav(1:nslit) = omi_slitwav_sol(1:nslit, ix)
     !      slitfit(1:nslit, :, :) = omi_solslitfit(1:nslit, :, :, ix)
     !   ENDIF
     !ENDIF
         
     solar_norm = SUM(allspec(1, spc_idx, 1:nspec)) / nspec
     n_fit_pts = nspec
     
     DO i = 1, ncoadd
        fitwavs   (1:nspec)  = allspec(i, wvl_idx, 1:nspec)
        currspec  (1:nspec)  = allspec(i, spc_idx, 1:nspec) / solar_norm
        fitweights(1:nspec)  = allspec(i, sig_idx, 1:nspec) / solar_norm
        
        CALL cal_fit_one (n_fit_pts, n_fitvar_sol, wrt_to_screen, wrt_to_file,&
             slitcal, slit_unit, tmpwave, tmp_fitvar, solfit_exval)
        
        IF (solfit_exval < 0) THEN
           WRITE(*, *) 'radwavcal_coadd: calibration not converge for pixel: ', i
           error = .TRUE.;           wshis(i) = 0.0; wsqus(i) = 0.0
        END IF

        ! Shift and squeeze earthshine spectrum
        allspec(i, wvl_idx, 1:nspec) = (fitwavs(1:nspec) - fitvar_sol(shi_idx)) / (1.0 + fitvar_sol(squ_idx)) 

        ! Make sure that the correction is less than a pixel
        dwvl = fitwavs(2) - fitwavs(1)
        IF ( ANY( ABS(allspec(i, wvl_idx, 1:nspec) - fitwavs(1:nspec)) > dwvl ) ) THEN
           allspec(i, wvl_idx, 1:nspec) = fitwavs(1:nspec) ! Roll back
           wshis(i) = 0.0; wsqus(i) = 0.0
        ELSE             
           wshis(i) = fitvar_sol(shi_idx); wsqus(i) = fitvar_sol(squ_idx)
        ENDIF
     ENDDO
  ELSE IF (wcal_bef_coadd .AND. .NOT. wavcal) THEN
     DO i = 1, ncoadd
        allspec(i, wvl_idx, 1:nspec) = (allspec(i, wvl_idx, 1:nspec) - wshis(i)) / (1.0 + wsqus(i)) 
     ENDDO        
  ENDIF

  IF (.NOT. wcal_bef_coadd) THEN  ! simple coadding due to unsuccessful calibration
     wshis = 0.0; wsqus = 0.0
     DO i = 2, ncoadd
        allspec(1, :, :)    = allspec(1, :, :) + allspec(i, :, :)
     ENDDO
     allspec(1, wvl_idx, :) = allspec(1, wvl_idx, :) / ncoadd  
     allspec(1, spc_idx, :) = allspec(1, spc_idx, :) / ncoadd  ! reduce S/N
     allspec(1, sig_idx, :) = allspec(1, sig_idx, :) / ncoadd  / SQRT(1.0 * ncoadd)  ! reduce S/N
  ELSE     
     ! Interpolate every other spectrum to the wavelength grids of first spectrum
     DO i = 2, ncoadd
        CALL interpolation (nspec, allspec(i, wvl_idx, 1:nspec), allspec(i, spc_idx, 1:nspec), &
             nspec - 2,  allspec(1, wvl_idx, 2:nspec-1), allspec(i, spc_idx, 2:nspec-1), errstat )
        IF ( errstat > pge_errstat_warning ) THEN
           errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; STOP 1
        END IF
        
        allspec(1, wvl_idx, 1) = allspec(1, wvl_idx, 1) + allspec(i, wvl_idx, 1)
        allspec(1, wvl_idx, nspec) = allspec(1, wvl_idx, nspec) + allspec(i, wvl_idx, nspec)
        allspec(1, spc_idx, :) = allspec(1, spc_idx, :) + allspec(i, spc_idx, :)
        allspec(1, sig_idx, :) = allspec(1, sig_idx, :) + allspec(i, sig_idx, :)
     ENDDO
     allspec(1, wvl_idx, 1) = allspec(1, wvl_idx, 1) / ncoadd  
     allspec(1, wvl_idx, nspec) = allspec(1, wvl_idx, nspec) / ncoadd  
     allspec(1, spc_idx, :) = allspec(1, spc_idx, :) / ncoadd  ! reduce S/N
     allspec(1, sig_idx, :) = allspec(1, sig_idx, :) / ncoadd / SQRT(1.0 * ncoadd)  ! reduce S/N    
  ENDIF
  
  RETURN
  
END SUBROUTINE radwavcal_coadd


! Properly align several cross-track spectra (to be coadded) to within one pixel 
SUBROUTINE prespec_align(nw, nspec, wavl, spec, prec, qflg)

  USE OMSAO_precision_module
  USE OMSAO_errstat_module
  
  IMPLICIT NONE
  
  ! ================================
  ! Input and Output variables
  ! =================================
  INTEGER, INTENT(IN)                                    :: nspec, nw
  REAL (KIND=r4),    DIMENSION(nw, nspec), INTENT(INOUT) :: wavl, spec, prec
  INTEGER (KIND=i2), DIMENSION(nw, nspec), INTENT(INOUT) :: qflg

  ! Local variables
  REAL (KIND=dp)                                         :: dwvl, mx_fwav 
  REAL (KIND=dp), DIMENSION(nw)                          :: diff
  INTEGER                                                :: fidx, i, j, nsel 

  mx_fwav = MAXVAL(wavl(1, 1:nspec))
 
  DO i = 1, nspec
     diff(1:nw) = ABS(wavl(1:nw, i) - mx_fwav)
     fidx       = MINVAL(MINLOC(diff(1:nw)))
     nsel = nw - fidx + 1
         
     ! Shift the spectrum
     wavl(1:nsel, i) = wavl(fidx:nw, i)
     spec(1:nsel, i) = spec(fidx:nw, i); spec(nsel+1:nw, i) = 0.0
     prec(1:nsel, i) = prec(fidx:nw, i) 
     qflg(1:nsel, i) = qflg(fidx:nw, i)  

     ! Make sure that wavelengths are increasing
     dwvl = wavl(2, i) - wavl(1, i)
     DO j = nsel + 1, nw
        wavl(j, i) = wavl(j-1, i) + dwvl
     ENDDO   
  ENDDO

  RETURN
END SUBROUTINE prespec_align
