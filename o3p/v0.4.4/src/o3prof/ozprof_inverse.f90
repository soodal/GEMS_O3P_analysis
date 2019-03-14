! ***********************************************************************************
! Author: Xiong Liu
! Date:   July 24, 2003
! Purpose: Ozone profile retrieval using PTR with GSVD-Lcurve/GCV
!		     or Optimal Estimation
!
! Modification History
! 1. xliu, jan 6, 2004, interface with optimal estimation
!    a. Remove convergence determination here (done in gsvd or oem)
!       The convergence criterion is determined just using relative
!       chisq change.
!    b. Remove writing fitting results for each iteration
!    c. Modify setting up return value
!    d. Deal with negative values in retrieved ozone (initialization
!       is too difficult from the actual)
! Exit status
! exval = 0, not converge
!       = 1, using absolute radiance converge in OE/PTR (delchi < epsrel)
!       = 2, converge in LIDORT (delchi < epsrel): 
!       = 4, ozone parameters change less than epsrel
!       =+10, used modified a  priori (e.g., ozone hole condition)
!       =+20, used UV2 retrievals (maybe modified a priori)
!       = +100, negative values occur for the last iteration 
!       = 3, 1 + 2      
!       = 5, 1 + 4     
!       = 6, 2 + 4
!       = 7, 4 + 2 + 1 
!       = -1, failure due to set_cldalb (specfit_ozprof.f90)
!       = -2, failure due to regular matrix with Tikhonov (specfit_ozprof.f90) or making atmosphere
!       = -3, failure due to out of bounds
!       = -4, failure in get_raman
!       = -5, failure in pseudo_model.f90
!       = -6, failure due to NaN values
! *********************************************************************************

SUBROUTINE ozprof_inverse (nf, varname, fitvar, fitvarap, lowbnd, upbnd,  &
     ns, np, sa, bb, nchisq, fitspec, fitres, exval)

  USE OMSAO_precision_module
  USE OMSAO_parameters_module,  ONLY: maxlay 
  USE ozprof_data_module,       ONLY: ffidx=>ozfit_start_index, flidx=>ozfit_end_index,     &
       ozwrtint_unit, ozwrtint, num_iter, avg_kernel, contri, covar, ncovar, ozdfs, ozinfo, &
       nlay, pfidx=>ozprof_start_index, plidx=>ozprof_end_index, ring_on_line,      &
       albfidx, nfalb, use_logstate, tf_fidx, tf_lidx, nt_fit, lcurve_unit, fgasidxs,       &
       gasidxs, ngas, saa_flag, tracegas, ozwrtfavgk, favg_kernel, radcalwrt, use_lograd,   &
       which_caloz, start_layer, end_layer, atmosprof, ozwrtwf, weight_function, do_simu, &
       wrtring, wfcfidx, nfwfc, ecfrind, ecfrfind, so2zfind, so2valts, do_twostep, do_bothstep, &
       use_large_so2_aperr, use_effcrs, do_simu_rmring, trace_avgk, trace_contri, trace_profwf, &
       ozwrtcontri, update_o3, update_sao3, which_toz
  USE OMSAO_variables_module,   ONLY: fitvar_rad, mask_fitvar_rad, epsrel, the_nspc,        &
       fitwavs, fitweights, fitvar_rad_apriori, maxit=>max_itnum_rad,  clmspec_rad,         &
       fitvar_rad_init, nradpix, numwin, npix_fitted, the_sza_atm, the_vza_atm, scnwrt,     &
       currpix, currline, currloop, the_surfalt, band_selectors
  USE OMSAO_indices_module,     ONLY: no2_t1_idx, so2_idx, bro_idx, hcho_idx, so2v_idx, &
       bro2_idx, o2o2_idx
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  INTEGER, INTENT (IN)                              :: ns, nf, np
  REAL (KIND=dp), DIMENSION(nf), INTENT (IN)        :: lowbnd, upbnd
  CHARACTER (LEN=6), DIMENSION(nf), INTENT (IN)     :: varname
  REAL (KIND=dp), DIMENSION(nf), INTENT (INOUT)     :: fitvarap
  REAL (KIND=dp), DIMENSION(nf, nf), INTENT(INOUT)  :: sa
  REAL (KIND=dp), DIMENSION(nf, nf), INTENT(IN)     :: bb

  ! ==================
  ! Modified variables  
  ! ==================
  REAL (KIND=dp), DIMENSION (nf), INTENT (INOUT) :: fitvar

  ! ================
  ! Output variables
  ! ================
  INTEGER,        INTENT (OUT)                   :: exval
  REAL (KIND=dp), INTENT (OUT)                   :: nchisq
  REAL (KIND=dp), DIMENSION (ns), INTENT (OUT)   :: fitspec, fitres

  ! ================
  ! Local Variables
  ! ================
  LOGICAL        :: refl_only, proceed, do_sa_diagonal, conv, varconv, negval, &
       last_iter, so2aperr_update, decorrelate, use_uv2init, correct_merr
  INTEGER        :: i, j, k, errstat, so2zfidx, uv2fy, uv2ly, nuv2, fidx, lidx, &
       cmerr_niter, uv12_retflg
  REAL (KIND=dp) :: ochisq, lchisq, sumn, sumo, avgres, delchi, sumdfs, &
       nradrms, oradrms, readout_noise
  REAL (KIND=dp), DIMENSION(nf)       :: delta_x, xold, xap, std, nstd, aperr, fitvarap0
  REAL (KIND=dp), DIMENSION(ns)       :: gspec, sig, fitspec1, fitres1, selidx, gspec_new, simrad
  REAL (KIND=dp), DIMENSION(ns, nf)   :: dyda
  REAL (KIND=dp), DIMENSION(ns+1)     :: sig1
  REAL (KIND=dp), DIMENSION(ns+1, nf) :: dyda1
  REAL (KIND=dp), DIMENSION(nf, ns+1) :: contri1
  REAL (KIND=dp), DIMENSION(nf, nf)   :: correl
  REAL (KIND=dp), DIMENSION(maxlay)   :: tozprof

  LOGICAL, SAVE                       :: first = .TRUE.
  INTEGER, SAVE                       :: no2fidx, so2vfidx, so2fidx, brofidx, hchofidx, o4fidx
  
  ! Variables for adjusting measurement errors based on fitting residuals
  LOGICAL                           :: adjust_merr
  INTEGER, PARAMETER                :: nreg = 5
  INTEGER, DIMENSION(nreg)          :: regfidxs, reglidxs, regnpts
  REAL (KIND=dp), DIMENSION(nreg)   :: reg_rms, reg_res
  REAL (KIND=dp), DIMENSION(0:nreg) :: reg_waves = &
       (/260.0, 290.0, 300.0, 310.0, 325.0, 350./)

  ! For detect a NaN
  INTEGER, PARAMETER :: DBPRECISION = SELECTED_INT_KIND(PRECISION(1.0d0))
  INTEGER (DBPRECISION), PARAMETER :: NAN = Z"7FF8000000000000"
  INTEGER, PARAMETER :: DPSB = BIT_SIZE(NAN) - 1

  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=14), PARAMETER :: modulename = 'ozprof_inverse'

  IF (first) THEN  
     DO i = 1, ngas
        IF (gasidxs(i) == no2_t1_idx) no2fidx  = fgasidxs(i)
        IF (gasidxs(i) == so2v_idx)   so2vfidx = fgasidxs(i)  
        IF (gasidxs(i) == so2_idx)    so2fidx  = fgasidxs(i)  
        IF (gasidxs(i) == bro_idx)    brofidx  = fgasidxs(i) 
        IF (gasidxs(i) == hcho_idx)   hchofidx = fgasidxs(i)    
        IF (gasidxs(i) == o2o2_idx)   o4fidx   = fgasidxs(i)  
     ENDDO
     
     first = .FALSE.
  ENDIF

  use_uv2init = .FALSE.
  IF (numwin >= 2 ) THEN
     fidx = 1
     DO i = 1, numwin
        lidx = fidx + nradpix(i) - 1 
        IF (band_selectors(i) == 2) THEN
           uv2fy = fidx; uv2ly = lidx; nuv2 = nradpix(i)
           use_uv2init = .TRUE.; EXIT
        ENDIF        
        fidx = lidx + 1
     ENDDO
  ENDIF
  uv12_retflg = 0 ! 0: uv1+uv2 1: uv1+uv2+modified a priori 2: uv2 retrieval only

  num_iter  = 0  
  refl_only = .FALSE.       ! need both radiances and weighting function
  proceed   = .TRUE.;    conv      = .FALSE.; varconv = .FALSE.
  do_sa_diagonal = .FALSE.; negval = .FALSE.
  sig = 1.0                 ! measurement error included in dyda and gspec (i.e., normalized to 1)

  exval = 0
  covar = 0.0; ncovar = 0.0 ! covariance matrix 
  fitvarap0 = fitvarap
  decorrelate = .FALSE.
  correct_merr = .FALSE. 
  cmerr_niter = 0 !maxit + 1
  readout_noise=1.0

  ! After retrievals are done with currently assumed measurement errors, preform retrievals again
  ! with adjusted measurement errors based on fitting residuals in several different spectral regions
  adjust_merr = .FALSE.
  IF (adjust_merr) THEN
     correct_merr = .FALSE. 
     fidx = 1;     i = fidx
     DO j = 1, nreg
        DO WHILE (fitwavs(i) < reg_waves(j) .AND. i <= ns)
           i = i + 1
        ENDDO
        lidx = i - 1
        regfidxs(j) = fidx; reglidxs(j) = lidx; regnpts(j) = lidx - fidx + 1
        fidx = lidx + 1
     ENDDO
  ENDIF

  ! Calculate ring spectrum
  IF (ring_on_line .AND. .NOT. (do_simu .AND. radcalwrt .AND. .NOT. wrtring .AND. .NOT. do_simu_rmring) ) THEN
     CALL GET_RAMAN(nlay, fitvar_rad(pfidx:plidx), errstat)
     IF (errstat == pge_errstat_error) THEN
        exval = -4; RETURN
     ENDIF
  ENDIF
    
  IF (ozwrtint) WRITE(ozwrtint_unit, '(A,I5,A10,I5, A10, I5)')  'Line = ', &
       currline, ' XPix = ', currpix, ' Loop = ', currloop

  ochisq = 10.0D20;  oradrms = 100.0
  inverse: DO WHILE (proceed)

     ! compute radiance, spectrum, and weighting function
     ! xliu: 0/1/28/1010: add use_effcrs here because use_hres always 
     ! needs weighting function for correction
     IF (radcalwrt .AND. do_simu .AND. .NOT. do_simu_rmring .AND. use_effcrs) THEN
        refl_only = .TRUE.
     ENDIF
     IF (num_iter == cmerr_niter .AND. correct_merr) fitweights(1:ns) = &
          fitweights(1:ns) / SQRT(1.0d0 * the_nspc) * readout_noise
!write(*,*) 'hello wasp! is it zero?' ,fitspec ! wasp : yes
!write(*,*) 'hello wasp! ozprof_inverse gspec here!' ,gspec
    CALL pseudo_model(num_iter, refl_only, ns, nf, fitvar, fitvarap, dyda, gspec, &
          fitres, fitspec, nchisq, nradrms, errstat)
!stop
!write(*,*) 'hello wasp! is it zero?' ,fitspec ! wasp : no
!write(*,*) 'hello wasp! pseudo model done!'
!stop
!write(*,*) 'hello wasp! fitvar here ozprof_inverse! ',fitvar

     IF (errstat == pge_errstat_error) THEN
        proceed = .FALSE.; exval = -5; CYCLE
     ENDIF
     IF ( radcalwrt .AND. do_simu .AND. .NOT. do_simu_rmring) THEN
        exval = 0; RETURN
     ENDIF

     xold = fitvar 
    
     IF (use_logstate) THEN
        xold(ffidx:flidx) = LOG(xold(ffidx:flidx)) 
        DO i = ffidx, flidx
           dyda(:, i) = dyda(:, i) * fitvar(i)
        ENDDO
     ENDIF
     xap = fitvarap
     IF (use_logstate) xap(ffidx:flidx) = LOG(xap(ffidx:flidx))

     delchi = ABS(nradrms - oradrms) / SQRT(oradrms)  
     IF (IEOR(IBCLR(TRANSFER(delchi, NAN), DPSB), NAN) == 0) THEN  ! check for NAN
        proceed = .FALSE.; exval = -6; CYCLE
     ENDIF

     ochisq = nchisq; oradrms = nradrms
     IF (delchi < epsrel) THEN     ! converge, exit
        exval = 1
        !proceed = .FALSE.; CYCLE
     ENDIF

     num_iter = num_iter + 1
     so2aperr_update = .FALSE.

     DO 
        last_iter = .TRUE.
        IF (.NOT. do_twostep) THEN  ! wasp : if(T)
           CALL oe_inversion (do_sa_diagonal, ozwrtint, ozwrtint_unit, epsrel,        &
                last_iter, num_iter, ns, nf, gspec, sig, dyda, xap, xold, sa, &
                varname, ffidx, flidx, delta_x, covar(1:nf, 1:nf), ncovar(1:nf, 1:nf),&
                conv, avg_kernel(1:nf, 1:nf), contri(1:nf, 1:ns), ozdfs, ozinfo, lchisq, gspec_new)
        !write(*,*) 'hello wasp! oe_inversion done!'
        !stop
        ELSE
           CALL twostep_inversion (do_sa_diagonal, ozwrtint, ozwrtint_unit, epsrel,        &
                last_iter, num_iter, ns, nf, gspec, sig, dyda, xap, xold, lowbnd, upbnd, sa, &
                varname, ffidx, flidx, delta_x, covar(1:nf, 1:nf), ncovar(1:nf, 1:nf),&
                conv, avg_kernel(1:nf, 1:nf), contri(1:nf, 1:ns), ozdfs, ozinfo, lchisq)  
        ENDIF

        IF ( radcalwrt .AND. do_simu .AND. do_simu_rmring) THEN
           ! Remove Ring effect and I/F wavelength shift
           simrad(1:ns) = fitspec(1:ns) - fitres(1:ns)
           !DO i = 1, ns
           !   WRITE(90, '(4D14.6)') fitwavs(i), fitspec(i), simrad(i), fitres(i)
           !ENDDO
           
           ! Ring effect, assume ring effect last variable, one single band
           gspec(1:ns) = gspec(1:ns) - dyda(1:ns, nf) * delta_x(nf)   

           ! I/F wavelength shift (should not since wavelength scale is not changed)
           !gspec(1:ns) = gspec(1:ns) - dyda(1:ns, nf-1) * delta_x(nf-1)
           fitres(1:ns) = gspec(1:ns) * fitweights(1:ns)
           fitspec(1:ns) = fitres(1:ns) + simrad(1:ns)
     
           exval = 0; RETURN
        ENDIF
        
        !write(*,*) 'hello wasp! ozprof_inverse fitvar here! ', fitvar(ffidx:flidx)
        !write(*,*) 'hello wasp! ozprof_inverse xold here! ', xold
        !write(*,*) 'hello wasp! ozprof_inverse delta_x here! ', delta_x
        fitvar = delta_x + xold
        !write(*,*) 'hello wasp! ozprof_inverse fitvar here! ', fitvar(ffidx:flidx)
!stop
        IF (IEOR(IBCLR(TRANSFER(lchisq, NAN), DPSB), NAN) == 0) THEN  ! check for NAN
           proceed = .FALSE.; exval = -6; EXIT
        ENDIF
        
        ! Special treatment for SO2
        IF (do_twostep .OR. (so2vfidx <= 0 .AND. so2fidx <= 0) .OR. use_large_so2_aperr) EXIT
        so2aperr_update = .FALSE.
        IF ( so2vfidx > 0) THEN
           IF ( ABS(fitvar(so2vfidx)) > SQRT(sa(so2vfidx, so2vfidx)) * 0.5) THEN
              sa(so2vfidx, so2vfidx) = 4.0 * fitvar(so2vfidx) ** 2.0 
              so2aperr_update  = .TRUE.
           ENDIF
        ENDIF
        
        IF ( so2fidx > 0 ) THEN
           IF ( ABS(fitvar(so2fidx)) > SQRT(sa(so2fidx, so2fidx)) * 0.5) THEN
              sa(so2fidx, so2fidx) = 4.0 * fitvar(so2fidx) ** 2.0 
              so2aperr_update = .TRUE.
           ENDIF
        ENDIF
        IF (.NOT. so2aperr_update ) EXIT
     ENDDO

     IF (use_logstate) THEN
        fitvar(ffidx:flidx) = EXP(fitvar(ffidx:flidx)) 
        xold(ffidx:flidx) = EXP(xold(ffidx:flidx)) 
        delta_x(ffidx:flidx) = fitvar(ffidx:flidx) - xold(ffidx:flidx)
     ENDIF
     !WRITE(*, '(I3,1X,A6,3d14.6)') ((i, varname(i), xold(i), &
     !     delta_x(i), fitvar(i)), i=1, nf) 
     !WRITE(*, *) SUM(xold(ffidx:flidx)), SUM(fitvar(ffidx:flidx)), SUM(xap(ffidx:flidx))

     uv12_retflg = 0
     IF (ANY (fitvar(ffidx:flidx) <= lowbnd(ffidx:flidx)) .OR. &
          ANY (fitvar(ffidx:flidx) >= upbnd(ffidx:flidx))) THEN                             
        IF (use_uv2init) THEN
           CALL negativeo3_inversion (uv2fy, uv2ly, nuv2, do_sa_diagonal,  &
                ozwrtint, ozwrtint_unit, epsrel, last_iter, num_iter, ns, nf, gspec, &
                sig, dyda, xap, xold, lowbnd, upbnd, sa, varname, ffidx, flidx,      &
                delta_x, covar(1:nf, 1:nf), ncovar(1:nf, 1:nf), conv, avg_kernel(1:nf, 1:nf), &
                contri(1:nf, 1:ns), ozdfs, ozinfo, lchisq, gspec_new, uv12_retflg)

           fitvar = delta_x + xold
           fitvarap(ffidx:flidx) = xap(ffidx:flidx)
        ENDIF
     END IF
     !WRITE(*, '(I3,1X,A6,3d14.6)') ((i, varname(i), xold(i), &
     !     delta_x(i), fitvar(i)), i=1, nf) 
     !WRITE(*, *) SUM(xold(ffidx:flidx)), SUM(fitvar(ffidx:flidx)), SUM(xap(ffidx:flidx))
     !write(*,*) 'hello wasp! ozprof_inverse fitvar here! ', fitvar(ffidx:flidx)
     !write(*,*) 'hello wasp! ozprof_inverse fitvar here! ', lowbnd(ffidx:flidx)
     !write(*,*) 'hello wasp! ozprof_inverse fitvar here! ', upbnd(ffidx:flidx) 

     IF (ANY (fitvar(ffidx:flidx) <= lowbnd(ffidx:flidx)) .OR. &
          ANY (fitvar(ffidx:flidx) >= upbnd(ffidx:flidx))) THEN                             
        IF (so2aperr_update) THEN
           ! Redo the retrieval the using original state except with updated SO2 fields
           IF (so2vfidx > 0) THEN
              xold(so2vfidx) = fitvar(so2vfidx); fitvar = xold
              fitvar_rad(mask_fitvar_rad(1:nf)) = fitvar
              fitvarap(so2vfidx) = fitvar(so2vfidx)
           ENDIF
        ELSE
           WRITE(*, *) modulename, ': Retrieved ozone values out of bounds!!!'
           proceed = .FALSE.;     exval = -3; CYCLE
        ENDIF        
     ELSE    
        ! Unnecessary here, it is the a priori error that matters
        !IF (so2aperr_update) THEN
        !   ! Update the a priori for the next iteration
        !   fitvarap(so2vfidx) = fitvar(so2vfidx)
        !ENDIF

        ! Update a priori value for ozone at every iterations
        !IF (num_iter >= 2) xap(ffidx:flidx) = fitvar(ffidx:flidx)

        ! update the uncondensed fitting variables
        fitvar_rad(mask_fitvar_rad(1:nf)) = fitvar
     END IF

     IF ( ALL(ABS(delta_x(ffidx:flidx) / fitvar(ffidx:flidx)) <= epsrel) ) THEN
        varconv = .TRUE.
     ENDIF

     IF (num_iter >= maxit .AND. .NOT. conv)  THEN
        proceed = .FALSE.;      exval = 0; CYCLE
     END IF

!     IF ((conv .OR. varconv)) THEN
     IF ( conv ) THEN
        IF (conv) exval = exval + 2
        IF (varconv) exval = exval + 4
        proceed = .FALSE.; CYCLE
     ENDIF

     IF (ring_on_line .AND. ABS(SUM(delta_x(ffidx:flidx))) > 20.0D0) THEN
        CALL GET_RAMAN(nlay, fitvar_rad(pfidx:plidx), errstat) 
        IF (errstat == pge_errstat_error) THEN 
           exval = -4; proceed = .FALSE.
        ENDIF
     ENDIF

     ! Check if need to update the SO2V profile shape
     IF (so2zfind > 0 ) THEN
        IF (fitvar(so2zfind) /= xold(so2zfind)) THEN
           IF (fitvar(so2zfind) > the_surfalt) THEN
              so2valts(0) = fitvar(so2zfind)
              so2valts(-1) = so2valts(0) - 1.0
              IF (so2valts(-1) < the_surfalt) so2valts(-1) = (the_surfalt + so2valts(0)) / 2.0
              so2valts(1)  = so2valts(0) + 1.0  
              ! Update the a priori value for SO2 plume height for the next iteration 
              ! (since the apriori value is an arbitrary guess)
              fitvarap(so2zfind) = fitvar(so2zfind)
           ELSE
              fitvar(so2zfind) = xold(so2zfind)
              fitvar_rad(mask_fitvar_rad(so2zfind)) = fitvar(so2zfind)
           ENDIF
           CALL ADJUST_SO2VPLUMEZ(errstat)
        ENDIF        
     ENDIF
 
     ! Updating O3 and Sao3 at first iteration
     ! When total ozone dependent climatology is selected.

       
     IF ( do_twostep == .False. .and. num_iter == 1) THEN 
       IF ( update_o3 == .TRUE. .or. update_sao3 ==.TRUE. .and. which_toz /= 0) THEN 
         print * , 'Update apriori information based on retrieved total o3'
         CALL update_o3_sao3 ( fitvar(ffidx:flidx), fitvarap(ffidx:flidx), sa(ffidx:flidx,ffidx:flidx ))
         
       ENDIF
     ENDIF

  END DO inverse

  DO i = 1, nf
     aperr(i) = SQRT(sa(i, i))
  ENDDO

  DO k = 1, ngas
     i = fgasidxs(k)
     IF (i > 0) THEN
        tracegas(k, 9) = avg_kernel(i, i)  
        ! consider interference from others (???)
        tracegas(k, 10) = SUM(avg_kernel(i, 1:nf) * aperr(1:nf) / aperr(i) )  
     ENDIF
  ENDDO
  !IF (ozwrtfavgk) favg_kernel(1:nf, 1:nf) = avg_kernel(1:nf, 1:nf)

  ! xliu: 03/19/2010
  ! Averaging kernels have already been calculated for each iteration
  ! No need to call oe_inversion unless:
  ! a. Decorrelate ozone with other varaibles
  ! b. Recalculate averaging kernels with corrected measurement errors
  ! 05/26/2010
  ! c. Adjust measurement error to reflect actual fitting residuals
  
  IF ( exval >= 0 .AND. (decorrelate .OR. &
       (correct_merr .AND. cmerr_niter == maxit + 1) .OR. adjust_merr) ) THEN   

     IF (correct_merr .AND. cmerr_niter == maxit + 1) THEN
        gspec(1:ns) = gspec(1:ns) * SQRT(1.0d0 * the_nspc) / readout_noise
        DO i = 1, nf
           dyda(1:ns, i) = dyda(1:ns, i) * SQRT(1.0d0 * the_nspc) / readout_noise
        ENDDO
        fitweights(1:ns) = fitweights(1:ns) / SQRT(1.0d0 * the_nspc) * readout_noise
        nchisq = nchisq  * the_nspc / readout_noise / readout_noise
     ENDIF

     IF (adjust_merr) THEN
        DO i = 1, nreg
           fidx = regfidxs(i); lidx = reglidxs(i)
           reg_rms(i) = SQRT(SUM((fitres(fidx:lidx)/fitweights(fidx:lidx))**2) / regnpts(i))
           reg_res(i) = SQRT(SUM(fitres(fidx:lidx)**2) / regnpts(i))

           fitweights(fidx:lidx) = fitweights(fidx:lidx) * reg_rms(i)
           gspec(fidx:lidx) = gspec(fidx:lidx) / reg_rms(i)
           dyda(fidx:lidx, 1:nf) = dyda(fidx:lidx, 1:nf) / reg_rms(i)
           nchisq  = SUM((fitres(1:ns) / fitweights(1:ns))**2.0)
           !WRITE(*, '(2F10.4,3I5,2D14.5)') reg_waves(i-1), reg_waves(i), &
           !     fidx, lidx, regnpts(i), reg_rms(i), reg_res(i)
        ENDDO
     ENDIF

     ! save standard deviations for those variables
     DO i = 1, nf
        std(i) = covar(i, i); nstd(i) = ncovar(i, i)
     ENDDO

     IF ( decorrelate )  THEN
        DO i = 1, nf
           IF (i < ffidx .OR. i > flidx ) THEN
              dyda(:, i)  = 0.0D0
              fitvarap(i) = fitvar(i)
           ENDIF
        ENDDO
     ENDIF
     !xold = fitvar  
     xap = fitvarap 

     IF (use_logstate) THEN
        xold(ffidx:flidx) = LOG(xold(ffidx:flidx)) 
        xap(ffidx:flidx)  = LOG(xap(ffidx:flidx)) 
     ENDIF

     last_iter = .TRUE.
     IF (.NOT. do_twostep) THEN
           CALL oe_inversion (do_sa_diagonal, ozwrtint, ozwrtint_unit, epsrel,        &
                last_iter, num_iter, ns, nf, gspec, sig, dyda, xap, xold, sa, &
                varname, ffidx, flidx, delta_x, covar(1:nf, 1:nf), ncovar(1:nf, 1:nf),&
               conv, avg_kernel(1:nf, 1:nf), contri(1:nf, 1:ns), ozdfs, ozinfo, lchisq, gspec_new)            
     ELSE
           CALL twostep_inversion (do_sa_diagonal, ozwrtint, ozwrtint_unit, epsrel,        &
                last_iter, num_iter, ns, nf, gspec, sig, dyda, xap, xold, lowbnd, upbnd, sa, &
                varname, ffidx, flidx, delta_x, covar(1:nf, 1:nf), ncovar(1:nf, 1:nf),&
                conv, avg_kernel(1:nf, 1:nf), contri(1:nf, 1:ns), ozdfs, ozinfo, lchisq) 
     ENDIF
        
     DO i = 1, nf
           contri(i, 1:ns) = contri(i, 1:ns) / fitweights(1:ns)
     ENDDO
        
     fitvar = xold  + delta_x

     !WRITE(*, '(I3,1X,A6,3d14.6)') ((i, varname(i), xold(i), delta_x(i), fitvar(i)), i=1, nf)  

     IF (use_logstate) THEN
        fitvar(ffidx:flidx) = EXP(fitvar(ffidx:flidx)) 
        xold(ffidx:flidx) = EXP(xold(ffidx:flidx)) 
        delta_x(ffidx:flidx) = fitvar(ffidx:flidx) - xold(ffidx:flidx)
     ENDIF

     IF (ANY (fitvar(ffidx:flidx) <= lowbnd(ffidx:flidx)) .OR. &
          ANY (fitvar(ffidx:flidx) >= upbnd(ffidx:flidx))) THEN

        WRITE(*, *) modulename, ': Retrieved ozone values out of bounds!!!'
        exval = -3
     ELSE       
        ! update the uncondensed fitting variables
        fitvar_rad(mask_fitvar_rad(1:nf)) = fitvar
     END IF

     DO i = 1, nf
        IF (i < ffidx .OR. i > flidx) THEN
           covar(i, i) = std(i); ncovar(i, i) = nstd(i) 
        ENDIF
     ENDDO
  ENDIF

  IF (ozwrtwf .AND. exval >= 0) THEN
     DO i = 1, nf  
        weight_function(1:ns, i) = dyda(1:ns, i) * fitweights(1:ns)
     ENDDO
  ENDIF

  DO i = 1, nf  
     contri(i, 1:ns) = contri(i, 1:ns) / fitweights(1:ns)
  ENDDO

  ! xliu: 08/06/2010, Add trace gas averaging kernel (dc / dx)
  ! dc/dx = dc/dY * dY/dx
  ! dc/dY: contribution function for VCD of a trace gas
  ! dY/dx: trace gas profile weighting function
  DO k = 1, ngas
     i = fgasidxs(k)
     IF (i > 0) THEN
        DO j = 1, nlay
           trace_avgk(k, j) = SUM(contri(i, 1:ns) * trace_profwf(k, 1:ns, j))
        ENDDO
      
        !WRITE(*, *) k, i, avg_kernel(i, i), SUM(trace_avgk(k, 1:nlay))
        !WRITE(*, '(12F8.5)') trace_avgk(k, 1:nlay)
     ENDIF
  ENDDO

  ! This is incorrect, and need to be changed
  IF  (use_logstate) THEN
     DO i = ffidx, flidx
        covar(i, i) = covar(i, i) * fitvar(i) ** 2.0
        ncovar(i, i) = ncovar(i, i) * fitvar(i) ** 2.0
        DO j = i+1, flidx
           covar(i, j) =  covar(i, j) * fitvar(i) * fitvar(j); covar(j, i) = covar(i,j) 
           ncovar(i, j) = ncovar(i, j) * fitvar(i) * fitvar(j); ncovar(j, i) = ncovar(i,j)
        ENDDO
     ENDDO
  ENDIF

  ! For exval  > 1, ideally still need to calculate final spectrum
  ! but this step is unnecssary, since retrieval is already done.
  ! On the other hand, it can save computation by one iteration (significant)
  ! for exval == 1, final spectrum already calculated
  ! IF (exval > 1) THEN  
  !   refl_only = .TRUE.
  !   CALL pseudo_model(num_iter, refl_only, ns, nf, fitvar, fitvarap, dyda, gspec, &
  !        fitres, fitspec, nchisq, nradrms, errstat)     
  !   delchi = ABS(nradrms - oradrms) / oradrms   ! converge, exit
  !   IF (delchi < epsrel)  exval = exval + 1
  !END IF

  IF (exval >= 0 .AND. radcalwrt) THEN
     refl_only = .TRUE.
     xold = fitvar     

     !Save retrievals  
     !Use a priori albedo and ozone, the other are the same
     IF (which_caloz == 1) THEN
        fitvar(ffidx:flidx) = fitvarap(ffidx:flidx) 
     ELSE
        tozprof(1:nlay)                = fitvar_rad(pfidx:plidx)
        tozprof(start_layer:end_layer) = fitvarap(ffidx:flidx)
        CALL get_caloz(nlay, atmosprof(1, 0:nlay), tozprof(1:nlay))
        fitvar(ffidx:flidx) = tozprof(start_layer:end_layer)         
     ENDIF

     IF (nfalb > 0) fitvar(albfidx:albfidx+nfalb-1) = fitvarap0(albfidx:albfidx+nfalb-1)
     IF (nfwfc > 0) fitvar(wfcfidx:wfcfidx+nfwfc-1) = fitvarap0(wfcfidx:wfcfidx+nfwfc-1)
     IF (ecfrfind > 0) fitvar(ecfrind) = fitvarap0(ecfrind)
     fitvar_rad(mask_fitvar_rad(1:nf)) = fitvar    

     CALL pseudo_model(num_iter, refl_only, ns, nf, fitvar, fitvarap, dyda, gspec, &
          fitres1, fitspec1, nchisq, nradrms, errstat)

     ! Restore the retrieved variables
     fitvar = xold
     fitvar_rad(mask_fitvar_rad(1:nf)) = fitvar

     IF (.NOT. use_lograd) THEN
        clmspec_rad(1:ns)  = fitspec1(1:ns) - fitres1(1:ns)
     ELSE
        clmspec_rad(1:ns)  = EXP(fitspec1(1:ns) - fitres1(1:ns))
     ENDIF
  ENDIF
  IF (exval > 0) THEN
     IF (uv12_retflg == 1) exval = exval + 10
     IF (uv12_retflg == 2) exval = exval + 20
     IF ( ANY(fitvar(ffidx:flidx) <= 0.0)) negval = .TRUE.
     IF (negval) exval = exval + 100 
  ENDIF


  RETURN

END SUBROUTINE ozprof_inverse

SUBROUTINE update_o3_sao3 (first_oz, ozprof, sao3)

  USE OMSAO_precision_module 
  USE OMSAO_parameters_module, ONLY: p0
  USE OMSAO_variables_module,  ONLY: the_month, the_day, the_lat
  USE ozprof_data_module,      ONLY: update_o3, update_sao3, which_clima, which_aperr, &
                                    atmosprof, ps0,nlay, nsfc , the_cfrac, &
                                    ozone_above60km, &
                                    smooth_ozbc, sacldscl0, pv811

  IMPLICIT NONE
  ! ======================
  ! Input/Output variables
  ! ======================

  REAL (KIND=dp), DIMENSION(nlay), INTENT(IN)    :: first_oz
  REAL (KIND=dp), DIMENSION(nlay), INTENT(INOUT) :: ozprof
  REAL (KIND=dp), DIMENSION(nlay, nlay), INTENT(INOUT) :: sao3

  ! ======================
  ! Local variables
  ! ======================
  REAL (KIND=dp) ::  toz, ozone_belowsfc

  INTEGER, PARAMETER ::  nref=71, nmpref=61, nmipas=121, nv8=11
  REAL (KIND=dp), DIMENSION(1:nmipas)       :: mipasp, mipast, mipaso3
  REAL (KIND=dp), DIMENSION(0:nref)         :: ozref, refp, refpg
  REAL (KIND=dp), DIMENSION(0:nlay) :: cum, cum1, ps, psg
  REAL (KIND=dp), DIMENSION(0:nv8) ::  pv8, pv8g, v8oz
  REAL (KIND=dp) ::  tmp, tmpsa
  INTEGER :: i , errstat

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=8), PARAMETER    :: modulename = 'update_o3 sao3'


  ps(0:nlay) = atmosprof(1, 0:nlay)
! calc toz
  
  IF (which_clima == 5 ) then
      toz = sum(first_oz) -ozone_above60km 
      
  ELSE IF ( which_clima == 6) then 
      toz = toz
  ENDIF
      
  IF (p0 > ps0 ) THEN ! ex) surface  = 5 km 
      ozone_belowsfc = first_oz(nsfc)*(p0- ps0)/(ps(nsfc) - ps(nsfc-1))
      toz = toz + ozone_belowsfc
  ENDIF

  !----------------------------
  ! updating ozone profile
  !----------------------------
  IF (update_o3 == .TRUE.) THEN 
      CALL REVERSE(ps(0:nlay), nlay +1 )       
      psg   = LOG(ps)
   
      IF (which_clima == 5 ) THEN 

         refp(0:nref) = (/(p0 * 10.0D0 ** (-REAL(i, KIND=dp)/16.D0), i=0, nref )/)

        ! Get a priori climatology for 60-70 km from MIPAS climatology
         CALL GET_MIPASIG2O3(the_month, the_day, the_lat, mipasp, mipaso3)
         ozref(nmpref:nref-1) = mipaso3(nmpref:nref-1)
        ! Ozone above 70 km and one more layer
         ozref(nref) = SUM(mipaso3(nref:nmipas))  
         CALL GET_IUPPROF (toz, ozref(1:nmpref-1))

         IF (ps0 > p0 ) THEN 
            tmp = ( ps0 - p0)/(refp(0)-refp(1))
            ozref(1) = ozref(1)*(1+tmp)
            refp(0) = ps0
         ENDIF

         IF (refp(nref) > ps(nlay) ) refp(nref) = ps(nlay)

         ozref(0) = 0.0
         DO i = 1, nref         
            ozref(i) = ozref(i-1) + ozref(i)
         ENDDO

         refpg = LOG(refp)  
         
    
         CALL BSPLINE(refpg(0:nref), ozref(0:nref), nref+1, psg(0:nlay), cum (0:nlay), nlay+1, errstat)
 
         IF (errstat < 0) THEN
            WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat; stop
         ENDIF
         cum = cum - cum(0) !; toz = umkoz(numk)
         ozprof(1:nlay) = cum(1:nlay) - cum(0:nlay-1)  

      ELSE IF (which_clima == 6 ) THEN 
         pv8 = pv811
         CALL GET_V8PROF(toz,  v8oz(1:nv8))
            v8oz(0) = 0.0   ! get cumulative ozone amount
         DO i = 1, nv8
            v8oz(i) = v8oz(i) + v8oz(i-1)
         ENDDO

      ! May have negative values for upper layers or lead to inproper interpolation 
      ! when interpolated from the 11 layer climatology and the partial column ozone 
      ! for these layers are very small, need readjustment using the McPeters Clima
         cum1 = cum

         IF (ps0 > p0 ) THEN 
            tmp = ( ps0 - p0)/(pv8(0)-pv8(1))
            v8oz(1) = v8oz(1)*(1+tmp)
            pv8(0) = ps0
         ENDIF

         IF (pv8(nv8) > ps(nlay) ) pv8(nv8) = ps(nlay)
 
         pv8g = LOG(pv8)
         CALL BSPLINE(pv8g(0:nv8), v8oz(0:nv8), nv8+1, psg(0:nlay), cum (0:nlay), nlay+1, errstat)
         IF (errstat < 0) THEN
         WRITE(*, *) modulename, ': BSPLINE error 1, errstat = ', errstat ; stop
         ENDIF
         
        DO i = nlay, 0, -1
           IF (psg(i) >= pv8g(nv8-1)) EXIT
        ENDDO
        tmp = (cum(nlay)-cum(i)) / (cum1(nlay)-cum1(i))

        cum = cum - cum(0); cum(1:nlay) = cum(1:nlay) - cum(0:nlay-1)  
        cum1 = cum1 - cum1(0); cum1(1:nlay) = cum1(1:nlay) - cum1(0:nlay-1)  
        cum(i+1:nlay) =cum1(i+1:nlay) * tmp


      ENDIF

      ! Reverse data on the retrieval grid
       CALL REVERSE(ozprof(1:nlay), nlay)  
      


  ENDIF

 
  !-----------------------------
  ! updating std profile 
  !-----------------------------
  IF (update_sao3 == .TRUE. .and. which_aperr == 5) THEN

     call get_apriori_covar(toz, ozprof, sao3)
     IF (nsfc < nlay) THEN
        sao3(nsfc+1:nlay, :) = 0.0; sao3(:, nsfc+1:nlay) = 0.0
     ENDIF  
  !xliu, 08/29/05 scaling a priori for layers below clouds to avoid smoothing even
  !for full cloudy conditions
     IF (.NOT. smooth_ozbc) THEN
     sacldscl0 = sacldscl0 * the_cfrac + (1.0 - the_cfrac )
     DO i = 1, nsfc
        IF (sacldscl0(i) < 1.0) THEN 
           tmpsa = sao3(i, i)
           sao3(i, 1:nlay) = sao3(i, 1:nlay) * sacldscl0(i)
           sao3(1:nlay, i) = sao3(1:nlay, i) * sacldscl0(i)
           sao3(i, i) = tmpsa
        ENDIF
     ENDDO
    ENDIF

  ENDIF

 RETURN

END SUBROUTINE update_o3_sao3

SUBROUTINE get_caloz (nl, pres, ozprof)

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY : mflay, rearth
  USE ozprof_data_module,      ONLY : caloz_fname, profunit
  IMPLICIT NONE
  
  ! ===============
  ! Input variables
  ! ===============
  INTEGER, INTENT (IN)                            :: nl
  REAL (KIND=dp), DIMENSION(0:nl), INTENT (INOUT) :: pres
  REAL (KIND=dp), DIMENSION(nl),   INTENT (INOUT) :: ozprof

  ! ==================
  ! Logical variables
  ! ==================
  INTEGER :: i, nz, oztyp, sidx, eidx, nlay, errstat
  REAL (KIND=dp)                                  :: alt, gcorr
  REAL (KIND=dp), DIMENSION(0:mflay)              :: zs
  REAL (KIND=dp), DIMENSION(0:nl)                 :: cozprof
  
  LOGICAL                                         :: first = .true.
  REAL (KIND=dp), DIMENSION(0:mflay), SAVE        :: ps, ozs, cozs


  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=9), PARAMETER :: modulename = 'get_caloz'

  IF ( first ) THEN 
     OPEN(profunit, FILE=TRIM(ADJUSTL(caloz_fname)), STATUS='old')
     READ(profunit, *) nz, oztyp
     IF (nz > mflay) THEN
        WRITE(*, *) modulename, ': Need to increase mflay!!!'; STOP
     ENDIF
     READ (profunit, *) 
     IF (oztyp <= 2) THEN
        READ (profunit, *) ozs(1:nz)
     ELSE
        READ (profunit, *) ozs(0:nz)
     ENDIF
     
     IF ( oztyp == 1 .AND. nz /= nl ) THEN
        WRITE(*, *) modulename, ': Number of layers are inconsistent!!!'; STOP 
     ELSE
        READ (profunit, *) 
        READ (profunit, *) ps(0:nz)  ! mb
     ENDIF
     CLOSE (profunit)
     
     IF (oztyp >= 2) THEN
        ! make profiles top-down
        IF ( ps(0) > ps(1) ) THEN
           CALl REVERSE(ozs(1:nz), nz)
           CALL REVERSE(ps(0:nz), nz+1)
        ENDIF
        
        ! Convert ozone from mixing ratio to partial column ozone
        IF (oztyp == 2) THEN        ! Get cumulative ozone profile
           cozs(0) = 0.0
           DO i = 1, nz
              cozs(i) = cozs(i-1) + ozs(i)
           ENDDO
        ELSE IF (oztyp == 3) THEN   ! Integrate from ppbv to DU (with gravity correction)
           zs = - 16.0 * LOG10( ps / 1013.25)
           cozs(0) = 0.0
           DO i = 1, nz             ! 2533.12 = 1.25 / 0.5 * 1013.25 
              alt     = ( zs(i-1) + zs(i) ) / 2.0
              gcorr   = ( rearth / (rearth + alt) ) ** 2.0 * 2533.125 !* 1.25 / 0.5 * 1013.25
              cozs(i) = cozs(i-1) + ( ozs(i-1) + ozs(i) ) * (ps(i) - ps(i-1)) / gcorr
           ENDDO
        ENDIF
     ENDIF
     
     first = .FALSE.
  ENDIF
     
  IF (oztyp == 1 ) THEN
     ozprof(1:nl) = ozs(1:nz) 
  ELSE
     ! Interpolate to the retrieval grid
     sidx = MINVAL(MINLOC(pres, MASK=(pres >= ps(0 )))) - 1
     eidx = MINVAL(MAXLOC(pres, MASK=(pres <= ps(nz)))) - 1
     nlay = eidx - sidx

     ps = LOG(ps); pres = LOG(pres)     
     CALL BSPLINE(ps, cozs, nz+1, pres(sidx:eidx), cozprof(sidx:eidx), nlay+1, errstat)
     ozprof(sidx+1:eidx) = cozprof(sidx+1:eidx) - cozprof(sidx:eidx-1)  
     ps = EXP(ps); pres = EXP(pres) 
     
  ENDIF

  RETURN

END SUBROUTINE get_caloz


SUBROUTINE twostep_inversion (do_sa_diagonal, ozwrtint, ozwrtint_unit, epsrel,  &
     last_iter, num_iter, ns, nf, gspec, sig, dyda, xap, xold, lowbnd, upbnd, sa, &
     varname, ffidx, flidx, delta_x, covar, ncovar, conv, avg_kernel, &
     contri, ozdfs, ozinfo, lchisq) 
  
  USE OMSAO_precision_module
  USE OMSAO_parameters_module,  ONLY: maxlay, max_spec_pts, elsunc_np, elsunc_nw 
  USE ozprof_data_module,       ONLY: fgasidxs, gasidxs, ngas, do_bothstep
  USE OMSAO_variables_module,   ONLY: fitvar_rad, mask_fitvar_rad, tol, epsabs,  epsx, &
       max_itnum_rad, step2_y, step2_dyda, fitwavs
  USE OMSAO_indices_module,     ONLY: elsunc_userdef
  USE bounded_nonlin_LS,        ONLY: elsunc
  USE OMSAO_errstat_module


  IMPLICIT NONE

  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER, INTENT (IN)                             :: ns, nf, ozwrtint_unit, num_iter, ffidx, flidx
  REAL (KIND=dp), INTENT (IN)                      :: epsrel
  REAL (KIND=dp), INTENT (OUT)                     :: ozdfs, ozinfo, lchisq
  REAL (KIND=dp), DIMENSION(nf), INTENT (IN)       :: lowbnd, upbnd, xap, xold
  REAL (KIND=dp), DIMENSION(nf), INTENT (OUT)      :: delta_x
  REAL (KIND=dp), DIMENSION(ns), INTENT (IN)       :: gspec, sig
  CHARACTER (LEN=6), DIMENSION(nf), INTENT (IN)    :: varname
  REAL (KIND=dp), DIMENSION(ns, nf), INTENT (INOUT):: dyda
  REAL (KIND=dp), DIMENSION(nf, ns), INTENT (OUT)  :: contri
  REAL (KIND=dp), DIMENSION(nf, nf), INTENT(INOUT) :: sa
  REAL (KIND=dp), DIMENSION(nf, nf), INTENT(OUT)   :: covar, ncovar, avg_kernel
  LOGICAL, INTENT(IN)                              :: do_sa_diagonal, ozwrtint, last_iter
  LOGICAL, INTENT(OUT)                             :: conv
  
  ! =======================
  ! Local variables
  ! =======================
  INTEGER                                          :: i, j, idx, n2f, elbnd, exval
  INTEGER, DIMENSION(nf)                           :: step2idxs
  REAL (KIND=dp), DIMENSION(nf, nf)                :: sa_sav
  REAL (KIND=dp), DIMENSION(nf)                    :: step2_fitvar
  CHARACTER (LEN=6), DIMENSION(nf)                 :: step2_varname
  REAL (KIND=dp), DIMENSION(ns)                    :: gspec1

  INTEGER,        DIMENSION (elsunc_np)            :: p
  REAL (KIND=dp), DIMENSION (elsunc_nw)            :: w
  REAL (KIND=dp), DIMENSION (nf)                   :: blow, bupp !, stderr
  REAL (KIND=dp), DIMENSION (ns)                   :: step2_fitres
  REAL (KIND=dp), DIMENSION (ns, nf)               :: dfda
  REAL (KIND=dp), DIMENSION (nf, nf)               :: correl

  EXTERNAL step2_specfit
  
  ! Freeze trace gas variables other than ozone by either
  ! 1. Setting weighting function to zero
  ! 2. Setting the a priori covariance matrix to zero
  sa_sav = sa

  n2f = 0
  DO i = 1, ngas
     IF (fgasidxs(i) > 0) THEN
        idx = fgasidxs(i)
        IF (.NOT. do_bothstep) sa(idx, idx) = 0.0
        n2f = n2f + 1
        step2idxs(n2f) = idx
     ENDIF
  ENDDO
  
  !DO i = flidx + 1, nf
  !   IF (varname(i)(3:4) == 'a1' .OR. varname(i)(3:4) == 'a2' .OR. varname(i)(3:4) == 'a3') THEN
  !      IF (.NOT. do_bothstep) sa(i, i) = 0.0
  !      n2f = n2f + 1
  !      step2idxs(n2f) = i
  !   ENDIF
  !ENDDO

  IF (n2f > 0) THEN
     step2_fitvar(1:n2f)  = xold(step2idxs(1:n2f))
     step2_dyda(1:ns, 1:n2f) = dyda(:, step2idxs(1:n2f))
     step2_varname(1:n2f) = varname(step2idxs(1:n2f))
  ENDIF

  delta_x = 0.0

  ! Call optimal estimation routine for those un-frozen parameters
  CALL oe_inversion (do_sa_diagonal, ozwrtint, ozwrtint_unit, epsrel,         &
       last_iter, num_iter, ns, nf, gspec, sig, dyda, xap, xold, sa,          &
       varname, ffidx, flidx, delta_x, covar(1:nf, 1:nf), ncovar(1:nf, 1:nf), &
       conv, avg_kernel(1:nf, 1:nf), contri(1:nf, 1:ns), ozdfs, ozinfo, lchisq, gspec1) 
  !WRITE(*, '(I3,1X,A6,2d14.6)') ((i, varname(i), xold(i), &
  !     delta_x(i)), i=1, nf)

  IF (n2f <= 0 ) RETURN

  ! Second step: use NLLS (ELSUNC, but linear here, i.e., 1 iterations) to
  ! fit those freezed trace gases 
  elbnd = 0   ! Unconstrained
  exval = 0         

  p  = -1   ;          p(1)   = 0;     p(3) = 1   !max_itnum_rad
  w  = -1.0 ;          w(1:4) = (/ tol,  epsrel,  epsabs,  epsx /)
  blow(1:n2f) = -1.0D99;   bupp(1:n2f) = 1.0D99
  blow(1:n2f) = -1.0D99;   bupp(1:n2f) = 1.0D99
  !blow(2) = 0.0; bupp(2) = 1.0
  !blow(5) = 0.0; bupp(5) = 1.0

  step2_fitres = 0.0;  dfda(1:ns, 1:n2f) = 0.0
  step2_y(1:ns) = gspec1(1:ns)

  CALL elsunc ( step2_fitvar(1:n2f), n2f, ns, ns, step2_specfit,          &
       elbnd, blow(1:n2f),  bupp(1:n2f), p, w, exval, step2_fitres(1:ns), &
       dfda(1:ns, 1:n2f) )

  ! ---------------------------------------------------------------
  ! Compute fitting residual
  ! FITRES is the negative of the returned function F = Model-Data.
  ! -------------------------------------opl--------------------------
  step2_fitres(1:ns) = -step2_fitres(1:ns)

  !WRITE(90, '(3D14.6)') ((gspec(i), gspec1(i), step2_fitres(i)), i=1, ns)

  ! Fitting RMS and CHI**2
  ! ----------------------
  lchisq = SQRT (SUM (step2_fitres(1:ns)**2 ) / ns )

  !! compute standard deviation for each variable
  !DO i = 1, n2f
  !   stderr(i) =  SQRT(dfda(i, i) * ns / (ns - n2f)) !* lchisq
  !END DO

  delta_x(step2idxs(1:n2f)) = delta_x(step2idxs(1:n2f)) + step2_fitvar(1:n2f)
  
  DO i = 1, n2f
     DO j = 1, n2f
        covar(step2idxs(i), step2idxs(j)) = covar(step2idxs(i), step2idxs(j)) + dfda(i, j)
        ncovar(step2idxs(i), step2idxs(j)) = ncovar(step2idxs(i), step2idxs(j)) + dfda(i, j)
     ENDDO
  ENDDO

  DO i = 1, ngas
     IF (fgasidxs(i) > 0) THEN
        idx = fgasidxs(i)
        sa(idx, idx) = sa_sav(idx, idx)
        !IF (ABS(xold(idx) + delta_x(idx)) > SQRT(sa(idx, idx))) &
        !     sa(idx, idx)= (xold(idx) + delta_x(idx))**2 
     ENDIF
  ENDDO

  !WRITE(*, '(I3,1X,A6,2d14.6)') ((i, varname(i), xold(i), delta_x(i)), i=1, nf)
  !print *, lchisq

  RETURN
END SUBROUTINE twostep_inversion

SUBROUTINE step2_specfit ( a, na, y, m, ctrl, dyda, mdy )

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY : max_spec_pts
  USE OMSAO_variables_module,  ONLY : step2_y, step2_dyda

  IMPLICIT NONE

  ! Input parameters
  ! ================
  INTEGER,                         INTENT (IN)  :: na, m, mdy
  REAL (KIND=dp), DIMENSION (na),  INTENT (IN)  :: a

  ! Modified parameters
  ! ===================
  INTEGER, INTENT (INOUT) :: ctrl

  ! Output parameters
  ! =================
  REAL (KIND=dp), DIMENSION (m),     INTENT (OUT) :: y
  REAL (KIND=dp), DIMENSION (m, na), INTENT (OUT) :: dyda

  ! Local variables
  ! ===============
  REAL (KIND=dp), DIMENSION (m) :: y0
  INTEGER                       :: i

  y0 = 0.0_dp
  DO i = 1, na
     y0 =  y0 + a(i) * (step2_dyda(1:m, i))
  END DO
 
  SELECT CASE ( ABS(ctrl) )
  CASE ( 1 )
     y  =  y0 - step2_y(1:m)
  CASE ( 2 )
     dyda(1:m, 1:na) = step2_dyda(1:m, 1:na)
  CASE ( 3 )
     ! This CASE is included to get the complete fitted spectrum
     y  = y0
  CASE DEFAULT
     WRITE (*, '(A,I3)') "Don't know how to handle CTRL = ", ctrl
  END SELECT

  RETURN
END SUBROUTINE step2_specfit


SUBROUTINE negativeo3_inversion (uv2fy, uv2ly, nuv2, do_sa_diagonal,  &
     ozwrtint, ozwrtint_unit, epsrel, last_iter, num_iter, ns, nf, gspec, &
     sig, dyda, xap, xold, lowbnd, upbnd, sa, varname, ffidx, flidx,      &
     delta_x, covar, ncovar, conv, avg_kernel, contri, ozdfs, ozinfo,     &
     lchisq, gspec_new, uv12_retflg) 
  
  USE OMSAO_precision_module
  USE OMSAO_parameters_module,  ONLY: maxlay, max_spec_pts
  USE ozprof_data_module,       ONLY: atmosprof, nlay
  USE OMSAO_variables_module,   ONLY: fitvar_rad, mask_fitvar_rad, fitwavs, &
       the_lat, the_month, the_day
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER, INTENT (IN)                             :: uv2fy, uv2ly, nuv2, &
       ns, nf, ozwrtint_unit, num_iter, ffidx, flidx
  INTEGER, INTENT (OUT)                            :: uv12_retflg
  REAL (KIND=dp), INTENT (IN)                      :: epsrel
  REAL (KIND=dp), INTENT (OUT)                     :: ozdfs, ozinfo, lchisq
  REAL (KIND=dp), DIMENSION(nf), INTENT (IN)       :: lowbnd, upbnd, xold
  REAL (KIND=dp), DIMENSION(nf), INTENT (OUT)      :: delta_x
  REAL (KIND=dp), DIMENSION(ns), INTENT (IN)       :: gspec, sig
  CHARACTER (LEN=6), DIMENSION(nf), INTENT (IN)    :: varname
  REAL (KIND=dp), DIMENSION(ns, nf), INTENT (INOUT):: dyda
  REAL (KIND=dp), DIMENSION(nf, ns), INTENT (OUT)  :: contri
  REAL (KIND=dp), DIMENSION(nf, nf), INTENT(INOUT) :: sa
  REAL (KIND=dp), DIMENSION(nf, nf), INTENT(OUT)   :: covar, ncovar, avg_kernel
  REAL (KIND=dp), DIMENSION(ns), INTENT (OUT)      :: gspec_new
  REAL (KIND=dp), DIMENSION(nf), INTENT (OUT)      :: xap
  LOGICAL, INTENT(IN)                              :: do_sa_diagonal, ozwrtint, last_iter
  LOGICAL, INTENT(OUT)                             :: conv
  
  ! =======================
  ! Local variables
  ! =======================
  INTEGER                         :: errstat
  REAL (KIND=dp)                  :: newtoz, oldtoz, aptoz
  REAL (KIND=dp), DIMENSION(nlay) :: newapoz, apoz
  REAL (KIND=dp), DIMENSION(nf) :: delta_x1, tmp_fitvar

  ! First try inversion with UV2 only
  uv12_retflg = 2
  CALL oe_inversion (do_sa_diagonal, ozwrtint, ozwrtint_unit, epsrel,         &
       last_iter, num_iter, nuv2, nf, gspec(uv2fy:uv2ly), sig(uv2fy:uv2ly),   &
       dyda(uv2fy:uv2ly, :), xap, xold, sa, varname, ffidx, flidx, delta_x,   &
       covar(1:nf, 1:nf), ncovar(1:nf, 1:nf), conv, avg_kernel(1:nf, 1:nf),   &
       contri(1:nf, uv2fy:uv2ly), ozdfs, ozinfo, lchisq, gspec_new(uv2fy:uv2ly)) 

  aptoz  = SUM(xap(ffidx:flidx))
  apoz   = xap(ffidx:flidx)
  oldtoz = SUM(xold(ffidx:flidx))
  newtoz = SUM(delta_x(ffidx:flidx)) + oldtoz

  !print *, 'Inside'
  !WRITE(*, '(4F8.3)') aptoz, oldtoz, newtoz, SUM(delta_x(ffidx:flidx))

  ! If UV2 total ozone is larger than a priori by 50 DU, then 
  ! use total ozone dependent TOMS V8 climatology to replace xap
  IF (ABS(newtoz - aptoz) >= 50.0 .AND. newtoz >= 100. .AND. newtoz <= 600.) THEN
     CALL get_tomsv8_clima(the_month, the_day, the_lat, newtoz, &
          nlay, atmosprof(1, 0:nlay), apoz(1:nlay), newapoz(1:nlay), errstat)

     IF (errstat /= pge_errstat_error) THEN
        xap(ffidx:flidx) = newapoz(1:nlay)

        ! Try rertievals with updated a priori and with both retrievals, accpet results
        ! if it is succesful, otherwise, still use UV2 retrievals
        CALL oe_inversion (do_sa_diagonal, ozwrtint, ozwrtint_unit, epsrel, &
             last_iter, num_iter, ns, nf, gspec, sig, dyda, xap, xold, sa,  &
             varname, ffidx, flidx, delta_x1, covar(1:nf, 1:nf), &
             ncovar(1:nf, 1:nf), conv, avg_kernel(1:nf, 1:nf),  &
             contri(1:nf, 1:ns), ozdfs, ozinfo, lchisq, gspec_new)
        
        tmp_fitvar = xold + delta_x1

        IF (ALL (tmp_fitvar(ffidx:flidx) > lowbnd(ffidx:flidx)) .AND. &
             ALL (tmp_fitvar(ffidx:flidx) < upbnd(ffidx:flidx))) THEN
           delta_x = delta_x1
           uv12_retflg = 1
           !print *, 'Use modified UV1 + UV2'
           !WRITE(*, '(4F8.3)') aptoz, oldtoz, oldtoz+SUM(delta_x(ffidx:flidx)), SUM(delta_x(ffidx:flidx)
        ELSE
           CALL oe_inversion (do_sa_diagonal, ozwrtint, ozwrtint_unit, epsrel,         &
                last_iter, num_iter, nuv2, nf, gspec(uv2fy:uv2ly), sig(uv2fy:uv2ly),   &
                dyda(uv2fy:uv2ly, :), xap, xold, sa, varname, ffidx, flidx, delta_x,   &
                covar(1:nf, 1:nf), ncovar(1:nf, 1:nf), conv, avg_kernel(1:nf, 1:nf),   &
                contri(1:nf, uv2fy:uv2ly), ozdfs, ozinfo, lchisq, gspec_new(uv2fy:uv2ly)) 
           !print *, 'Still use UV2'
        ENDIF
     ENDIF
  ENDIF

  RETURN
END SUBROUTINE negativeo3_inversion



