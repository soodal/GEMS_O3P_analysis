  ! *********************************************************************************
  ! Author: Xingu Liu
  ! Date:   July 24, 2003
  ! Purpose: calculated simulated reflectance and weighting function/first derivative
  ! for all the fitting variables. Basically, call the lidort to calculate radiance
  ! and ozone, albedo weighting functions, where species other than ozone are not
  ! taken into account. Then calibrate the measured radiance/solar spectra to obtain
  ! the measured reflectance and use the finite difference to obtain the first 
  ! derivative for all the other variables.
  ! *********************************************************************************
  ! Need to do more about the albedo 

SUBROUTINE pseudo_model (num_iter, refl_only, ns, nf, fitvar, fitvarap, dyda, gspec,    &
     fitres, fitspec, chisq, relrms, errstat)
  
  USE OMSAO_precision_module
  USE OMSAO_parameters_module,ONLY : maxlay
  USE OMSAO_indices_module,   ONLY : shi_idx, squ_idx, max_calfit_idx, ring_idx,  &
       max_rs_idx, mxs_idx, maxalb, maxoth, solar_idx, bro_idx, bro2_idx, so2_idx,   &
       no2_t1_idx, hcho_idx, shift_offset, us1_idx, us2_idx, maxwfc, so2v_idx, o2o2_idx
!      comm_idx, com1_idx, fsl_idx, rsl_idx, com2_idx, com3_idx
  USE OMSAO_variables_module, ONLY : fitwavs, fitweights, sza => the_sza_atm, &
       vza => the_vza_atm, aza => the_aza_atm, sca=>the_sca_atm, fitvar_rad,  &
       mask_fitvar_rad, yn_varyslit, the_lat, the_lon, rad_wav_avg,           &
       lo_radbnd, up_radbnd, database_shiwf, slwf, npix_fitted, database,     &
       fitvar_rad_str, numwin, nradpix, band_selectors, refidx, npix_fitting, &
       scnwrt, refspec_norm, rmask_fitvar_rad, database_indices,              &
       fitvar_rad_apriori, actspec_rad
  USE ozprof_data_module,     ONLY : nlay, use_lograd,                       &
       ozf_fidx => ozfit_start_index, ozf_lidx => ozfit_end_index,           &
       ozp_fidx=>ozprof_start_index, ozp_lidx => ozprof_end_index,           &
       stlay => start_layer, endlay => end_layer, albfidx, nalb, nfalb,      &
       albidx, albmin, albmax, albfpix, alblpix, t_fidx, t_lidx, tf_fidx,    &
       tf_lidx, nt_fit, do_subfit, saa_flag, fgasidxs, tracegas, ngas,       &
       gasidxs, osind, osfind, slind, slfind, shind, shfind, rnind, rnfind,  &
       dcind, dcfind, isind, isfind, irind, irfind, oswins, slwins, shwins,  &
       rnwins, dcwins, iswins, irwins, nos, nsl, nsh, nrn, ndc, nis, nir,    &
       nothgrp, the_cfrac, polcorr, radcalwrt, do_simu, ecfrfind, ecodfind,  &
       ectpfind, taodfind, twaefind, saodfind, ecfrind, ecodind, ectpind,    &
       taodind, twaeind, saodind, sprsind, sprsfind, wfcfidx, nwfc, nfwfc,   &
       wfcidx, wfcmin, wfcmax, wfcfpix, wfclpix, so2zind, so2zfind, fit_atanring, &
       use_effcrs, ncalcp, do_simu_rmring, nsaa_spike, the_ai
  USE OMSAO_errstat_module

  IMPLICIT NONE
  
  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER, INTENT(IN)                             :: ns, nf, num_iter
  INTEGER, INTENT(OUT)                            :: errstat
  LOGICAL, INTENT(INOUT)                          :: refl_only
  REAL (KIND=dp), INTENT(INOUT), DIMENSION (nf)   :: fitvar, fitvarap
  REAL (KIND=dp), INTENT(OUT), DIMENSION (ns)     :: gspec, fitres, fitspec
  REAL (KIND=dp), INTENT(OUT), DIMENSION (ns, nf) :: dyda
  REAL (KIND=dp), INTENT(OUT)                     :: chisq, relrms

  ! ===============
  ! Local variables
  ! ===============
  INTEGER, PARAMETER :: MSTKS = 4
  INTEGER :: n0alb, n0wfc, i, j, k, iw, ReturnStatus, ridx, sidx, fidx, lidx, &
       idx, albord, min_ssa_iter, swin, ewin, ig, nord, nostk, wfcord, ntmp
  INTEGER, DIMENSION(maxalb)               :: albpmax, albpmin
  INTEGER, DIMENSION(maxwfc)               :: wfcpmax, wfcpmin
  INTEGER, DIMENSION(numwin, maxoth)       :: tmpind, tmpfind
  INTEGER, DIMENSION(maxoth, 2)            :: tmpwins
  INTEGER, DIMENSION(nf)                   :: newind
  REAL (KIND=dp), DIMENSION(maxalb)        :: albarr 
  REAL (KIND=dp), DIMENSION(maxwfc)        :: wfcarr 
  REAL (KIND=dp), DIMENSION (ns)           :: delpos, waves, meas1, meas2, &
       simrad, simrad1, fitspec1, diff, temporwf
  REAL (KIND=dp), DIMENSION(ns,nlay,MSTKS) :: ozwf, tmpwf
  REAL (KIND=dp), DIMENSION(ns, 4)         :: albothwf, wfcothwf
  REAL (KIND=dp), DIMENSION(ns, MSTKS)     :: o3shiwf, cfracwf, albwf, fsimrad, &
       ctpwf, codwf, saodwf, taodwf, twaewf, sprswf, so2zwf
  REAL (KIND=dp), DIMENSION(numwin, maxoth):: o3shi
  REAL (KIND=dp), DIMENSION(nlay)          :: tprof, ozprof, ozadj, ozaprof
  !REAL (KIND=dp), DIMENSION (n_radwvl_sav) :: waves_sav, delpos_sav
  REAL (KIND=dp), DIMENSION (nf)           :: fitvar_saved
  REAL (KIND=dp)                           :: rms, radrms, wavavg, cfrac, newoz
  REAL (KIND=dp), DIMENSION (numwin)       :: allrms, allchisq, allradrms

  LOGICAL :: do_ozwf, do_albwf, do_o3shi, do_tmpwf, do_shiwf, do_taodwf, do_twaewf, &
       do_saodwf, do_cfracwf, do_ctpwf, do_codwf, negval, do_sprswf, do_so2zwf
  LOGICAL, DIMENSION (nlay)     :: ozvary
  REAL (KIND=dp), DIMENSION(ns) :: walb0s 

  ! xliu, 08/10/2010
  ! Current VLIDORT calculation is based on single surface albedo (per channel). 
  ! The wavelength dependence of surface albedo on radiance is corrected through
  ! weighting function. Howver, there are not accounted for in the calculation of
  ! weighting functions. 
  LOGICAL :: vary_sfcalb                 

  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=12), PARAMETER :: modulename = 'pseudo_model'

  ! To compute fine and radiance calculaiton wavelngth grids when using high resolution wavelength grids
  LOGICAL, SAVE             :: first = .TRUE.

  errstat = pge_errstat_ok

  IF (first .AND. .NOT. use_effcrs) THEN
     CALL get_hres_radcal_waves(errstat)
     IF (errstat == pge_errstat_error) THEN
        WRITE(*, *) modulename, ': Errors in getting fine & radiance calculation wavelength grids!!!'
        RETURN
     ENDIF
     first = .FALSE.
  ENDIF
 
  ! ================ Determine flags for linearization ======================
  
  do_albwf  = .TRUE.; do_ozwf = .TRUE.; do_o3shi = .TRUE.; do_tmpwf = .TRUE.; do_cfracwf = .TRUE.
  IF (nfalb <= 0) do_albwf = .FALSE.
  IF (nfwfc <= 0) do_cfracwf = .FALSE.
  IF (refl_only .AND. .NOT. use_effcrs) THEN
     do_ozwf = .FALSE.; do_o3shi = .FALSE.; do_tmpwf = .FALSE.
  END IF
  IF (nos <= 0)    do_o3shi = .FALSE.
  IF (nt_fit <= 0) do_tmpwf = .FALSE. 
  ozvary = .FALSE.; ozvary(stlay:endlay) = .TRUE.
  
  do_ctpwf  = .FALSE.;  do_codwf  = .FALSE.
  do_taodwf  = .FALSE.; do_twaewf = .FALSE.;  do_saodwf = .FALSE.
  do_sprswf = .FALSE. ; do_so2zwf = .FALSE. 
  IF (.NOT. refl_only) THEN
     IF (ecfrfind > 0) do_cfracwf = .TRUE.
     IF (ecodfind > 0) do_codwf   = .TRUE.
     IF (ectpfind > 0) do_codwf   = .TRUE.
     IF (taodfind > 0) do_taodwf  = .TRUE.
     IF (saodfind > 0) do_saodwf  = .TRUE.
     IF (twaefind > 0) do_twaewf  = .TRUE.
     IF (sprsfind > 0) do_sprswf  = .TRUE.
     IF (so2zfind > 0) do_so2zwf  = .TRUE.
  ENDIF

  ! Update cloud fraction or disable fitting cloud fraction when fc<=1.0E-3 or fc>=0.999
  IF ( ecfrfind > 0 ) THEN
     the_cfrac = fitvar_rad(ecfrind)
     IF (the_cfrac <= 1.0E-3 ) THEN
        the_cfrac = 0.0D0;  fitvar(ecfrfind) = 0.0D0; fitvar_rad(ecfrind) = 0.0D0
        fitvarap(ecfrfind) = 0.0D0;    fitvar_rad_apriori(ecfrind) = 0.0D0
        do_cfracwf = .FALSE.      
     ELSE IF (the_cfrac >= 0.999 ) THEN
        the_cfrac = 1.0D0;  fitvar(ecfrfind) = 1.0D0; fitvar_rad(ecfrind) = 1.0D0
        fitvarap(ecfrfind) = 1.0D0;    fitvar_rad_apriori(ecfrind) = 1.0D0
        do_cfracwf = .FALSE. 
     
        ! When using Lambertian cloud model, surface albedo is already fitted, so nothing needs to be done
        ! But when using scattering cloud model, need to fit an optical thickness and a surface albedo
        ! To be implemented !!!
     ENDIF
  ENDIF
  ! ======= Set up ozone, temperature, trace gases, albedo, lamda for LIDORT ============
  ozprof(1:nlay)  = fitvar_rad (ozp_fidx:ozp_lidx)
  ozaprof(1:nlay) = fitvar_rad_apriori(ozp_fidx:ozp_lidx)
  tprof(1:nlay)   = fitvar_rad(t_fidx:t_lidx)
   
  !WRITE(*, '(A)') 'Initial Guess (pseudo_model-1): '
  !WRITE(*, '(12F8.3)') fitvar_rad(ozp_fidx:ozp_lidx), SUM(fitvar_rad(ozp_fidx:ozp_lidx))
  !WRITE(*, '(A)') 'Initial Guess (pseudo_model-2): '
  !WRITE(*, '(12F8.3)') fitvar(ozf_fidx:ozf_lidx), SUM(fitvar(ozf_fidx:ozf_lidx))

  !xliu (02/01/2007): adjust ozone profile for negative ozone values
  !                   radiances will be corrected using ozone weighting function
  ozadj(1:nlay) = 0.0; negval = .FALSE.
  DO i = 1, nlay
     IF (ozprof(i) <= 0.0) THEN
        newoz  = MIN(0.5, ozaprof(i))
        negval = .TRUE. ; ozadj(i)  = newoz - ozprof(i); ozprof(i) = newoz
        do_ozwf = .TRUE.; ozvary(i) = .TRUE.
     ENDIF
  ENDDO

 DO k = 1, ngas
    i = fgasidxs(k)
    IF (i > 0) THEN
       j = mask_fitvar_rad(i)
       tracegas(k, 4) = fitvar_rad(j) !/ refspec_norm(gasidxs(k)) ! trace gas column in molecules cm-2
    ENDIF
   
 ENDDO
 waves = fitwavs(1:ns)
 
 vary_sfcalb = .FALSE. !.TRUE.
 n0alb = 0
 
 DO i = 1, nalb
     j = albidx - 1 + i
!xliu, 02/08/2012, add albord and **albord
     READ(fitvar_rad_str(j)(4:4), '(I1)') albord

     IF (albord == 0) THEN
!     IF (fitvar_rad_str(j)(4:4) == '0') THEN
        n0alb = n0alb + 1
        
        !IF (fitvar_rad(j) > 1.0) THEN
        !   fitvar_rad(j) = 1.0
        !   IF (rmask_fitvar_rad(j) > 0) fitvar(rmask_fitvar_rad(j)) = 1.0
        !ENDIF
        albarr(n0alb)  = fitvar_rad(j)
        albpmax(n0alb) = alblpix(i); albpmin(n0alb) = albfpix(i)
        IF (vary_sfcalb) walb0s(albpmin(n0alb):albpmax(n0alb)) = albarr(n0alb)
     ELSE
        IF (vary_sfcalb) THEN
           fidx=albfpix(i); lidx=alblpix(i)
           wavavg = SUM(waves(fidx:lidx)/(1.0+lidx-fidx))
           walb0s(fidx:lidx) = walb0s(fidx:lidx) + fitvar_rad(j) * (waves(fidx:lidx) - wavavg)**albord
        ENDIF
     ENDIF        
  ENDDO
   
  n0wfc = 0
  DO i = 1, nwfc
     j = wfcidx - 1 + i
     IF (fitvar_rad_str(j)(4:4) == '0') THEN
        n0wfc = n0wfc + 1; wfcarr(n0wfc) = fitvar_rad(j)
        wfcpmax(n0wfc) = wfclpix(i); wfcpmin(n0wfc) = wfcfpix(i)
     ENDIF      
  ENDDO
  
  IF (do_subfit) THEN
     DO i = 1, maxoth
        o3shi(1:numwin, i) = fitvar_rad(osind(1:numwin, i))
     ENDDO
  ELSE
     o3shi(1, 1:maxoth)    = fitvar_rad(osind(1, 1:maxoth))
  ENDIF
     
 !Iraddaince/radiance shift is done when interpolating solar reference to wavelength grid, not here
 !IF (nsh > 0) THEN
 !   IF (do_subfit) THEN
 !      fidx = 1
 !      DO j = 1, numwin
 !         lidx = fidx + nradpix(j) - 1
 !         delpos(fidx:lidx) =  waves(fidx:lidx) - (waves(fidx) + waves(lidx)) / 2.0
 !         IF (shfind(j, 1) > 0)  waves(fidx:lidx) = waves(fidx:lidx) + fitvar_rad(shind(j, 1)) 
 !
 !         DO i = 2, nsh            
 !            IF (shfind(j, i) > 0)  waves(fidx:lidx) = waves(fidx:lidx) + &
 !                 fitvar_rad(shind(j, i)) * (delpos(fidx:lidx) ** (i-1))
 !         ENDDO
 !         fidx = lidx + 1
 !      ENDDO
 !   ELSE
 !      IF (shwins(1, 1) == 1) THEN
 !         fidx = 1
 !      ELSE
 !         fidx = SUM(nradpix(1: shwins(1, 1)-1)) + 1
 !      ENDIF
 !      lidx = SUM(nradpix(1: shwins(1, 2)))
 !      delpos(fidx:lidx) =  waves(fidx:lidx) - (waves(fidx) + waves(lidx)) / 2.0
 !      IF (shfind(1, 1) > 0) waves(fidx:lidx)  = waves(fidx:lidx) + fitvar_rad(shind(1, 1))
 !
 !      DO i = 2, nsh  
 !         IF (shfind(1, i) > 0) THEN
 !            waves(fidx:lidx)  = waves(fidx:lidx) + fitvar_rad(shind(1, i)) * (delpos(fidx:lidx) ** (i-1))
 !         ENDIF
 !      ENDDO
 !   ENDIF
 !ENDIF

  ! === Call LIDORT, polarization correction, and additional wf calc =====
  !the_cfrac = 0.
  nostk = 1

  IF (use_effcrs) THEN
     CALL LIDORT_PROF_ENV(do_ozwf, do_albwf, do_tmpwf, do_o3shi, ozvary,       &
          do_taodwf, do_twaewf, do_saodwf, do_cfracwf, do_ctpwf, do_codwf, do_sprswf, do_so2zwf, &
          ns, waves, maxoth, o3shi, sza, vza, aza, nlay, ozprof, tprof, n0alb, &
          albarr, albpmin, albpmax, vary_sfcalb, walb0s, n0wfc, wfcarr, wfcpmin, wfcpmax, nostk, &
          albwf(1:ns, 1:nostk), ozwf(1:ns, 1:nlay, 1:nostk), tmpwf(1:ns, 1:nlay, 1:nostk),           &
          o3shiwf(1:ns, 1:nostk), cfracwf(1:ns, 1:nostk), codwf(1:ns, 1:nostk), &
          ctpwf(1:ns, 1:nostk), taodwf(1:ns, 1:nostk), twaewf(1:ns, 1:nostk), &
          saodwf(1:ns, 1:nostk), sprswf(1:ns, 1:nostk), so2zwf(1:ns, 1:nostk), fsimrad(1:ns, 1:nostk), errstat)
     IF (errstat == pge_errstat_error) &
          WRITE(*, *) modulename, ': Errors in calling LIDORT_PROF_ENV!!!'
  ELSE    
       !      print * , 'o3shi', o3shi(:,1)
       !      print * , 'albarr',albarr(1:2)
       !      write(*,'(a10, 10e15.7)')  'walb0s',walb0s(201:210)
       !      write(*,'(a10, 10e15.7)')  'albwf', albwf(201:210,1)

     ntmp = MAX(ncalcp, ns) ! 94, 218
     CALL HRES_RADCALC_ENV(ntmp, do_ozwf, do_albwf, do_tmpwf, do_o3shi, ozvary,            &
          do_taodwf, do_twaewf, do_saodwf, do_cfracwf, do_ctpwf, do_codwf, do_sprswf,do_so2zwf, &
          ns, waves, maxoth, o3shi, sza, vza, aza, nlay, ozprof, tprof, n0alb,  &
          albarr, albpmin, albpmax, vary_sfcalb, walb0s, n0wfc, wfcarr, wfcpmin, wfcpmax, nostk, &
          albwf(1:ns, 1:nostk), ozwf(1:ns, 1:nlay, 1:nostk), tmpwf(1:ns, 1:nlay, 1:nostk), &
          o3shiwf(1:ns, 1:nostk), cfracwf(1:ns, 1:nostk), codwf(1:ns, 1:nostk),                      &
          ctpwf(1:ns, 1:nostk), taodwf(1:ns, 1:nostk), twaewf(1:ns, 1:nostk), saodwf(1:ns, 1:nostk), &
          sprswf(1:ns, 1:nostk), so2zwf(1:ns, 1:nostk), fsimrad(1:ns, 1:nostk), errstat)
          
        !     write(*,'(a10, 10e15.7)')  'albwf', albwf(201:210,1)
     IF (errstat == pge_errstat_error) &
          WRITE(*, *) modulename, ': Errors in calling HRES_RADCALC_ENV!!!'
  ENDIF
  IF (errstat == pge_errstat_error) RETURN
 

  !xliu (02/01/2007): correct radiances based on ozone weighting function to deal with negative ozone values
  IF (negval) THEN
     DO i = 1, nlay 
        IF (ozadj(i) > 0) THEN
           fsimrad(1:ns, 1:nostk) = fsimrad(1:ns, 1:nostk) - ozadj(i) * ozwf(1:ns, i, 1:nostk) 
        ENDIF
     ENDDO
  ENDIF
  simrad = fsimrad(1:ns, 1)

  !IF (radcalwrt .AND. do_simu) THEN
  !   albwf(1:ns, 1) = 0.0
  !   o3shiwf(1:ns, 1) = 0.0
  !   ozwf(1:ns, :, 1) = 0.0
  !   database(so2_idx, refidx(1:ns)) = 0.0
  !   database(bro_idx, refidx(1:ns)) = 0.0
  !   database(bro2_idx, refidx(1:ns)) = 0.0
  !   database(hcho_idx, refidx(1:ns)) = 0.0
  !   database(no2_t1_idx, refidx(1:ns)) = 0.0    
  !ENDIF
    
  !PRINT *, polcorr, ns, nlay
  !print *, sza, aza, vza
  !DO i = 1, ns
  ! WRITE(77, '(f12.5, 43d14.6)') fitwavs(i), simrad(i), ozwf(i, 1:nlay, 1), albwf(i, 1)!, taodwf(i, 1) !&
  !        !, o3shiwf(i, 1) !* ozprof(1:nlay)
  !ENDDO
  !CLOSE(77)
  !STOP
!!$
!!$  do i = 1,ns
!!$
!!$  write(*,*)  i, so2zwf(i, 1)
!!$  enddo

  ! get dlnI/dx
  IF (use_lograd) THEN
     IF (do_albwf)   albwf(1:ns, 1)   = albwf(1:ns, 1)   / simrad     
     IF (do_o3shi)   o3shiwf(1:ns, 1) = o3shiwf(1:ns, 1) / simrad
     IF (do_codwf)   codwf(1:ns, 1)   = codwf(1:ns, 1)   / simrad
     IF (do_sprswf)  sprswf(1:ns, 1)  = sprswf(1:ns, 1)  / simrad
     IF (do_so2zwf)  so2zwf(1:ns, 1)  = so2zwf(1:ns, 1)  / simrad
     IF (do_ctpwf)   ctpwf(1:ns, 1)   = ctpwf(1:ns, 1)   / simrad
     IF (do_cfracwf) cfracwf(1:ns, 1) = cfracwf(1:ns, 1) / simrad
     IF (do_taodwf)  taodwf(1:ns, 1)  = taodwf(1:ns, 1)  / simrad
     IF (do_saodwf)  saodwf(1:ns, 1)  = saodwf(1:ns, 1)  / simrad
     IF (do_twaewf)  twaewf(1:ns, 1)  = twaewf(1:ns, 1)  / simrad
      
     IF (do_ozwf) THEN
        DO i = stlay, endlay
           ozwf (:, i, 1) = ozwf(:, i, 1) / simrad
           IF (nt_fit > 0) tmpwf(:, i, 1) = tmpwf(:, i, 1) / simrad
        END DO
     ENDIF

     simrad = LOG(simrad)           ! get dlnI     
  END IF
 
  !WRITE(77, *) 'Ozone weighting function D(lnI)/D(lnx)'
  !DO i = 1, ns
  !  WRITE(77, '(f8.4, 43d14.6)') fitwavs(i), ozwf(i, 1:nlay) !* ozprof(1:nlay)
  !ENDDO

  ! correct for linear/quardratic wavelength dependent in albedo
  albothwf = 0.0
  DO i = 1, nalb
     j = albidx + i - 1

     READ(fitvar_rad_str(j)(4:4), '(I1)') albord
     IF (albord == 0) CYCLE

     fidx=albfpix(i); lidx=alblpix(i)
     wavavg = SUM(waves(fidx:lidx)/(1.0+lidx-fidx))
     !albothwf(fidx:lidx,albord) = albwf(fidx:lidx, 1)*(waves(fidx:lidx) / wavavg)**albord
     albothwf(fidx:lidx,albord) = albwf(fidx:lidx, 1)*(waves(fidx:lidx) - wavavg)**albord ! much better than above
     IF (.NOT. vary_sfcalb) simrad(fidx:lidx) = simrad(fidx:lidx) +  albothwf(fidx:lidx, albord) * fitvar_rad(j)         
     
  ENDDO
         
  IF (nwfc > 0) THEN
     wfcothwf = 0.0
     DO i = 1, nwfc
        j = wfcidx + i - 1
        
        READ(fitvar_rad_str(j)(4:4), '(I1)') wfcord
        IF (wfcord == 0) CYCLE
        
        fidx=wfcfpix(i); lidx=wfclpix(i)
        wavavg = SUM(waves(fidx:lidx)/(1.0+lidx-fidx))
        wfcothwf(fidx:lidx,wfcord) = cfracwf(fidx:lidx, 1)*(waves(fidx:lidx) - wavavg)**wfcord 
        simrad(fidx:lidx) = simrad(fidx:lidx) +  wfcothwf(fidx:lidx, wfcord) * fitvar_rad(j)        
     ENDDO
  ENDIF

  IF (radcalwrt .AND. do_simu .AND. .NOT. do_simu_rmring) THEN
     fitspec = actspec_rad(1:ns)
     IF ( use_lograd ) fitspec = LOG( fitspec )
     fitres = fitspec - simrad
     chisq  = SUM((fitres / fitweights(1:ns))**2.0)
     RETURN
  ENDIF
    
  ! get calibrated reflectance, correcting for ring, undersampling, and trace gases
  ! calculate weighting function for shift parameter if required, just need to do it once
  IF (num_iter == 0) THEN 
     do_shiwf = .TRUE.
  ELSE
     do_shiwf = .FALSE.
  ENDIF
  CALL spectra_reflectance (ns, nf, fitvar, do_shiwf, fitspec, errstat)
  IF (errstat == pge_errstat_error) THEN
     WRITE(*, *) modulename, ': Errors in spectra_reflectance!!!'; RETURN
  ENDIF
  IF (num_iter == 0) THEN
     CALL UV1_SPIKE_DETECT(ns, fitspec, simrad, nsaa_spike)
  ENDIF
           
  ! get residual between measured and simulated reflectance
  fitres = fitspec - simrad

  ! compute chi-square difference
  chisq  = SUM((fitres / fitweights(1:ns))**2.0)
  rms =    SQRT(chisq  / REAL(ns, KIND=dp))
  relrms = 100.D0 * SQRT(SUM(ABS((simrad-fitspec) / fitspec)**2.0) &
          / REAL(ns, KIND=dp))

  IF (scnwrt) THEN
     fidx = 1
     DO i = 1, numwin
        lidx = fidx + nradpix(i) - 1
        allchisq(i) = SUM((fitres(fidx:lidx) / fitweights(fidx:lidx))**2.0)
        allrms(i)   = SQRT(allchisq(i) / nradpix(i))
        allradrms(i) = 100.D0 * SQRT(SUM(ABS((simrad(fidx:lidx)-fitspec(fidx:lidx)) &
             / fitspec(fidx:lidx))**2) / nradpix(i)) 
        fidx = lidx + 1
     ENDDO

     ! Relative rms difference between calculated and simulated log-radiances
     ! /radiances depending on the flag use_lograd
     IF (use_lograd) THEN     
        simrad1 = EXP(simrad); fitspec1 = EXP(fitspec)
        radrms = 100.D0 * SQRT(SUM(ABS((simrad1-fitspec1) / fitspec1)**2) &
             / REAL(ns, KIND=dp))   
        
        fidx = 1
        DO i = 1, numwin
           lidx = fidx + nradpix(i) - 1
           allradrms(i) = 100.D0 * SQRT(SUM(ABS((simrad1(fidx:lidx)-fitspec1(fidx:lidx)) &
                / fitspec1(fidx:lidx))**2) / nradpix(i))
           fidx = lidx + 1
        ENDDO
     ELSE 
        radrms = relrms
     ENDIF
     
     WRITE (*, '(I5, 4(A10,1pd11.3))') num_iter, ' Chi = ', chisq, ' rms = ', rms, &
          ' rms(%) = ', relrms, ' Irms(%) = ', radrms
     DO i = 1, numwin
        WRITE (*, '(A13, I2, 2(A14, 1pd11.3))') 'Win ', i, ': allrms = ', &
             allrms(i), 'allIrms(%) = ', allradrms(i)
     ENDDO

  ENDIF

  !DO i = 1, ns
  !   WRITE(90, '(F8.3, 2D14.6)') fitwavs(i), EXP(fitspec(i)), EXP(simrad(i))
  !ENDDO
  !CLOSE (90)
  !print *, sza, vza, aza
  !STOP

  IF (npix_fitted == 0) THEN
     min_ssa_iter = 2
  ELSE
     min_ssa_iter = 1
  ENDIF
   
  ! South Atlantic Anomaly Correction, use 10 channel 2 pixels as reference because they are
  ! less subject to South Atlantic Anomaly
  IF (num_iter == min_ssa_iter .AND. num_iter <=min_ssa_iter+2 .AND. .NOT. refl_only .AND. &
       band_selectors(1) == 1 .AND. saa_flag  ) THEN 
     CALL spike_detect_correct(nradpix(1)+10, fitspec(1:nradpix(1)+10), simrad(1:nradpix(1)+10))  
  ENDIF
  fitres = fitspec - simrad
    
  IF (.NOT. refl_only) THEN 
     dyda = 0.0
    
     ! albedo weighting function
     DO i = 1, nfalb
        j = albfidx + i - 1
        k = mask_fitvar_rad(j) 
        fidx = albfpix(k -albidx + 1); lidx=alblpix(k - albidx + 1)

        READ(fitvar_rad_str(k)(4:4), '(I1)') albord

        IF (albord == 0) THEN
           dyda(fidx:lidx, j)=albwf(fidx:lidx, 1)
        ELSE 
           dyda(fidx:lidx, j)=albothwf(fidx:lidx, albord)
        ENDIF
              
        ! xliu: 07/01/2010, compute aerosol index (defined with relative to 20 nm distance)
        ! defined similar to TOMS aerosol index
        IF (albord == 1) THEN
           the_ai = (dyda(lidx, j) - dyda(fidx, j)) * fitvar(j) * 100. / (waves(lidx)-waves(fidx)) * 20.
        
        ENDIF
     ENDDO  
     ! wavelength-dependent cloud fraction weighting function
     DO i = 1, nfwfc
        j = wfcfidx + i - 1
        k = mask_fitvar_rad(j) 
        fidx = wfcfpix(k -wfcidx + 1); lidx=wfclpix(k - wfcidx + 1)

        READ(fitvar_rad_str(k)(4:4), '(I1)') wfcord

        IF (wfcord == 0) THEN
           dyda(fidx:lidx, j) = cfracwf(fidx:lidx, 1)
        ELSE 
           dyda(fidx:lidx, j) = wfcothwf(fidx:lidx, wfcord)
        ENDIF
     ENDDO  

     ! ozone profile and temperature weighting function
     dyda(:, ozf_fidx:ozf_lidx) = ozwf(:, stlay:endlay, 1) 
     IF (nt_fit > 0) THEN
        fidx = stlay + mask_fitvar_rad(tf_fidx)-maxlay-ozp_fidx
        lidx = endlay + mask_fitvar_rad(tf_lidx)-maxlay-ozp_lidx
        dyda(:, tf_fidx:tf_lidx) =tmpwf(1:ns, fidx:lidx, 1) 
     ENDIF

     ! Aerosol/cloud parameters 
     IF (taodfind > 0) dyda(:, taodfind) = taodwf(:, 1)
     IF (saodfind > 0) dyda(:, saodfind) = saodwf(:, 1)
     IF (twaefind > 0) dyda(:, twaefind) = twaewf(:, 1)
     IF (ecfrfind > 0) dyda(:, ecfrfind) = cfracwf(:, 1)
     IF (ecodfind > 0) dyda(:, ecodfind) = codwf(:, 1)
     IF (sprsfind > 0) dyda(:, sprsfind) = sprswf(:, 1)
     IF (so2zfind > 0) dyda(:, so2zfind) = so2zwf(:, 1)
     IF (ectpfind > 0) dyda(:, ectpfind) = ctpwf(:, 1)     

     ! get 1st derivative for calibration and reference parameters
     ! wfs for undersampling, common mode, other gases and shift terms
   
     fitvar_saved(1:nf) = fitvar(1:nf); do_shiwf = .FALSE.
     DO i = 1, ozf_fidx - 1     
        ! shift parameter, already calculate dR/dS
      
        IF (mask_fitvar_rad(i) > shift_offset) CYCLE 

        ! check for shift indices
     
        ridx = database_indices(i) ! indices in the database
        sidx = shift_offset + ridx ! shift indice in fitvar_rad
          
        IF (ridx == so2_idx .OR. ridx == hcho_idx .OR. ridx == no2_t1_idx .OR. &
           ridx == bro_idx .OR. ridx == bro2_idx .OR. ridx == so2v_idx .OR. ridx == o2o2_idx) THEN
        !IF (ridx /= us1_idx .AND. ridx /= us2_idx .AND. ridx /= comm_idx .AND. ridx /= com1_idx &
        !     .AND. ridx /= com2_idx .AND. ridx /= com3_idx .AND. ridx /= fsl_idx .AND. ridx /= rsl_idx) THEN
           ! for trace gases, wfs are just the cross sections (or amf * cross sections)
           dyda(:, i) = -database(ridx, refidx(1:ns))
           IF (.NOT. use_lograd)  dyda(:, i) =  dyda(:, i) * simrad
        ! xliu, 11/01/2011, the following for undersampling is incorrect
        !ELSE  IF (ridx == us1_idx .OR. ridx == us2_idx) THEN
        !   dyda(:, i) = - EXP(fitspec(1:ns)) * database(ridx, refidx(1:ns))   
        ELSE
           fitvar(i) = fitvar_saved(i) * 1.001           
           CALL spectra_reflectance (ns, nf, fitvar, do_shiwf, meas1, errstat)
           IF (errstat == pge_errstat_error) THEN
              WRITE(*, *) modulename, ': Errors in spectra_reflectance!!!'
              fitvar(i) = fitvar_saved(i); RETURN
           ENDIF
           fitvar(i) = fitvar_saved(i) * 0.999
           CALL spectra_reflectance (ns, nf, fitvar, do_shiwf, meas2, errstat)  
           IF (errstat == pge_errstat_error) THEN
              WRITE(*, *) modulename, ': Errors in spectra_reflectance!!!'
              fitvar(i) = fitvar_saved(i); RETURN
           ENDIF      
           dyda(:, i) = -(meas1 - meas2) / (0.002 * fitvar_saved(i))
           fitvar(i) = fitvar_saved(i)
        ENDIF 

        ! check for shifting
        IF (rmask_fitvar_rad(sidx) > 0) THEN
           dyda(:, rmask_fitvar_rad(sidx)) = dyda(:, i) * database_shiwf(ridx, refidx(1:ns))  
        ENDIF
     ENDDO
     
     ! Weighting function when use fit_atanring
     ! 1. analytical  2. finite difference (use 2 to validate 1)
       
     IF (fit_atanring) THEN
        do_shiwf = .FALSE.
        DO i = rnfind(1, 1), rnfind(1, 3)
           fitvar(i) = fitvar_saved(i) * 1.001           
           CALL spectra_reflectance (ns, nf, fitvar, do_shiwf, meas1, errstat)
           IF (errstat == pge_errstat_error) THEN
              WRITE(*, *) modulename, ': Errors in spectra_reflectance!!!'
              fitvar(i) = fitvar_saved(i); RETURN
           ENDIF
           fitvar(i) = fitvar_saved(i) * 0.999
           CALL spectra_reflectance (ns, nf, fitvar, do_shiwf, meas2, errstat)  
           IF (errstat == pge_errstat_error) THEN
              WRITE(*, *) modulename, ': Errors in spectra_reflectance!!!'
              fitvar(i) = fitvar_saved(i); RETURN
           ENDIF      
           dyda(:, i) = -(meas1 - meas2) / (0.002 * fitvar_saved(i))
           fitvar(i) = fitvar_saved(i)
        ENDDO

        !! Analytical weighting function
        !delpos = (fitwavs-fitvar(rnfind(1, 2))) / fitvar(rnfind(1, 3))
        !dyda(:, rnfind(1, 1)) = database(ring_idx, refidx(1:ns)) * (atan(delpos)+1.54223)
        !dyda(:, rnfind(1, 2)) = -database(ring_idx, refidx(1:ns)) / (1.0 + delpos**2) * &
        !     fitvar(rnfind(1, 1)) / fitvar(rnfind(1, 3))
        !dyda(:, rnfind(1, 3)) = dyda(:, rnfind(1, 2)) * delpos             
     ENDIF

     ! wfs for ozcrs shift, slit_shift, wavelength shift, Ring effect, degradation correction,
     ! internal scattering in irradiance, internal scattering in radiance
     DO ig = 1, nothgrp
        IF (ig == 1) THEN            ! ozone cross section  
           nord = nos; tmpind = osind; tmpfind = osfind; tmpwins = oswins
           temporwf = -o3shiwf(1:ns, 1)
        ELSE IF ( ig == 2) THEN      ! Radiance/Irradiance Slit Difference
           nord = nsl; tmpind = slind; tmpfind = slfind; tmpwins = slwins
           temporwf = -slwf(1:ns)
        ELSE IF (ig == 3)  THEN      ! Radiance/irradince wavelength shift
           nord = nsh; tmpind = shind; tmpfind = shfind; tmpwins = shwins
        ELSE IF (ig == 4)  THEN      ! Ring effect
           IF (fit_atanring) CYCLE
           nord = nrn; tmpind = rnind; tmpfind = rnfind; tmpwins = rnwins
           temporwf = -database(ring_idx, refidx(1:ns))
           IF (.NOT. use_lograd) temporwf(1:ns) = temporwf * simrad
        ELSE IF (ig == 5 ) THEN      ! Degradation correction
           nord = ndc; tmpind = dcind; tmpfind = dcfind; tmpwins = dcwins
        ELSE IF (ig == 6)  THEN      ! Internal scattering in irradinace
           nord = nis; tmpind = isind; tmpfind = isfind; tmpwins = iswins
        ELSE IF (ig == 7) THEN       ! Internal scattering in radiance
           nord = nir; tmpind = irind; tmpfind = irfind; tmpwins = irwins
        ENDIF

        ! Use finite differences to get the zero-order weighting functions
        IF (nord > 0 .AND. ig >= 3 .AND. ig <= nothgrp .AND. ig /= 4) THEN
           temporwf = 0.0

           IF (do_subfit) THEN
              swin = tmpwins(1, 1); ewin = tmpwins(1, 2)
           ELSE
              swin = 1; ewin = 1
           ENDIF

           IF (ig == 3) THEN  ! irradiance/radiance wavelength shift
              fitvar(tmpfind(swin:ewin, 1)) = 0.001 
              CALL spectra_reflectance (ns, nf, fitvar, do_shiwf, meas1, errstat) 
              IF (errstat == pge_errstat_error) THEN
                 WRITE(*, *) modulename, ': Errors in spectra_reflectance!!!'
                 fitvar(tmpfind(swin:ewin, 1)) = fitvar_saved(tmpfind(swin:ewin, 1)); RETURN
              ENDIF
              fitvar(tmpfind(swin:ewin, 1)) = -0.001 
              CALL spectra_reflectance (ns, nf, fitvar, do_shiwf, meas2, errstat)
              IF (errstat == pge_errstat_error) THEN
                 WRITE(*, *) modulename, ': Errors in spectra_reflectance!!!'
                 fitvar(tmpfind(swin:ewin, 1)) = fitvar_saved(tmpfind(swin:ewin, 1)); RETURN
              ENDIF
              temporwf = -(meas1 - meas2) / 0.002
              fitvar(tmpfind(swin:ewin, 1)) = fitvar_saved(tmpfind(swin:ewin, 1))
           ELSE
              fitvar(tmpfind(swin:ewin, 1)) = fitvar_saved(tmpfind(swin:ewin, 1)) * 1.001 
              CALL spectra_reflectance (ns, nf, fitvar, do_shiwf, meas1, errstat) 
              IF (errstat == pge_errstat_error) THEN
                 WRITE(*, *) modulename, ': Errors in spectra_reflectance!!!'
                 fitvar(tmpfind(swin:ewin, 1)) = fitvar_saved(tmpfind(swin:ewin, 1)); RETURN
              ENDIF
              fitvar(tmpfind(swin:ewin, 1))  = fitvar_saved(tmpfind(swin:ewin, 1)) * 0.999
              CALL spectra_reflectance (ns, nf, fitvar, do_shiwf, meas2, errstat)
              IF (errstat == pge_errstat_error) THEN
                 WRITE(*, *) modulename, ': Errors in spectra_reflectance!!!'
                 fitvar(tmpfind(swin:ewin, 1)) = fitvar_saved(tmpfind(swin:ewin, 1)); RETURN
              ENDIF
              fitvar(tmpfind(swin:ewin, 1))  = fitvar_saved(tmpfind(swin:ewin, 1))

              IF (swin == ewin) THEN
                 IF (fitvar_saved(tmpfind(swin, 1)) /= 0.0) THEN
                    temporwf = -(meas1 - meas2) / (0.002 * fitvar_saved(tmpfind(swin, 1)))
                 ELSE
                    temporwf = 0.0d0
                 ENDIF
              ELSE
                 IF (swin == 1) THEN
                    fidx = 1
                 ELSE
                    fidx = SUM(nradpix(1: swin-1)) + 1
                 ENDIF

                 DO iw = swin, ewin
                    lidx = fidx + nradpix(iw) - 1
                    IF (fitvar_saved(tmpfind(iw, 1)) /= 0.0) THEN
                       temporwf(fidx:lidx) = -(meas1(fidx:lidx) - meas2(fidx:lidx)) / &
                            (0.002 * fitvar_saved(tmpfind(iw, 1)))
                    ELSE
                       temporwf(fidx:lidx) = 0.d0
                    ENDIF
                    fidx = lidx + 1
                 ENDDO
              ENDIF
           ENDIF
        ENDIF

        IF (nord > 0) THEN
           IF (do_subfit) THEN
              fidx = 1
              DO j = 1, numwin
                 lidx = fidx + nradpix(j) - 1
                 IF ( tmpfind(j, 1) > 0 ) dyda(fidx:lidx, tmpfind(j, 1)) = temporwf(fidx:lidx)

                 delpos(fidx:lidx) = fitwavs(fidx:lidx) - (fitwavs(fidx) + fitwavs(lidx)) * 0.5
                 DO i = 2, nord
                    IF ( tmpfind(j, i) > 0 ) THEN
                       dyda(fidx:lidx, tmpfind(j, i)) = dyda(fidx:lidx, tmpfind(j, i-1)) * delpos(fidx:lidx)
                    ENDIF
                 ENDDO
                 fidx = lidx + 1
              ENDDO
           ELSE 
              IF (tmpwins(1, 1) == 1) THEN
                 fidx = 1
              ELSE
                 fidx = SUM(nradpix(1: tmpwins(1, 1) - 1)) + 1
              ENDIF
              lidx = SUM(nradpix(1: tmpwins(1, 2)))

              dyda(fidx:lidx, tmpfind(1, 1)) = temporwf(fidx:lidx)
              delpos(fidx:lidx) = fitwavs(fidx:lidx) - (fitwavs(fidx) + fitwavs(lidx)) * 0.5
              
              DO i = 2, nord
                 IF ( tmpfind(1, i) > 0 ) THEN
                    dyda(fidx:lidx, tmpfind(1, i)) = dyda(fidx:lidx, tmpfind(1, i-1)) * delpos(fidx:lidx)
                 ENDIF
              ENDDO
           ENDIF
        ENDIF
     ENDDO

    !WRITE(92, *) ns, nf
    !DO i = 1, ns
    !   WRITE(92, '(f10.4, 80d14.6)') fitwavs(i), dyda(i, 1:nf)
    !ENDDO
    !errstat = pge_errstat_error; RETURN
    !
    !CLOSE(92)
    !STOP
     
     DO i = 1, nf  
        dyda(:, i) = dyda(:, i) / fitweights(1:ns)
     END DO
     
     ! finnally obtain the new spectrum to be fitted in the GSVD
     gspec(1:ns) = fitres(1:ns) / fitweights(1:ns)
     
     ! Restore the unperturbated fitting variables
     fitvar(1:nf) = fitvar_saved(1:nf)
     fitvar_rad(mask_fitvar_rad(1:nf)) = fitvar(1:nf)  
  END IF
  

  RETURN
  
END SUBROUTINE pseudo_model
