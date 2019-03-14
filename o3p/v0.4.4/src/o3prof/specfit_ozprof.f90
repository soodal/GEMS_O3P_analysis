  ! ***************************************************************
  ! Author:  Xiong Liu
  ! Date:    July 24, 2003
  ! Purpose: Ozone profile retrieval using PTR or OE
  ! 
  ! Modification History:
  ! 1.  xliu, Jan 19, 2004
  !     a. Remove using ELSUNC 
  !     b. Add a interface for preparing atmosphere using tomsv8 
  !        profile, EP total ozone, temperature profiles, surface 
  !        pressure
  !     c. Add a interface to prepare a priori covariance  
  !     d. Calculate error in total ozone here
  !     e. Move get_reg_matrix here
  !     f. Prepare measurment error and measurement vector
  ! ***************************************************************

SUBROUTINE specfit_ozprof (initval, fitcol, dfitcol, rms, exval)
  
  USE OMSAO_precision_module
  USE OMSAO_parameters_module,   ONLY: missing_value_dp, maxlay, du2mol
  USE OMSAO_indices_module,      ONLY: n_max_fitpars, wvl_idx, spc_idx,         &
       sig_idx, maxalb, no2_t1_idx, so2_idx, so2v_idx, bro_idx, hcho_idx, us1_idx,        &
       us2_idx, max_calfit_idx, max_rs_idx, mxs_idx, maxoth, shift_offset, maxwfc, &
       comvidx, cm1vidx, comfidx, cm1fidx, bro2_idx, o2o2_idx, solar_idx ! wasp
  USE OMSAO_variables_module,    ONLY: curr_rad_spec, rad_wav_avg, fitwavs,     &
       currspec, fitweights, fitvar_rad, fitvar_rad_apriori, fitvar_rad_saved,  &
       fitvar_rad_init, fitvar_rad_str, lo_radbnd, up_radbnd, n_fitvar_rad, &
       n_rad_wvl,mask_fitvar_rad,fitspec_rad, fitres_rad, use_meas_sig,        &
       weight_rad, yn_doas, fitvar_rad_std, the_lon, the_lat, the_month,        &
       the_year, the_day, npix_fitted, edgelons, edgelats,  currloop,           &
       fitvar_rad_nstd, numwin, nradpix, refspec_norm, scnwrt, ozabs_convl, the_surfalt, &
       fitvar_rad_init_saved, the_lons, the_lats, the_surfalt, nloc, fitvar_rad_aperror, &
       the_sza_atm, currline, refwvl, n_refwvl, reduce_resolution,curr_sol_spec,n_irrad_wvl,&
       database, database_shiwf, the_sza_atm, the_vza_atm, the_aza_atm, atmdbdir ! atmdbdir added : geun  
  USE ozprof_data_module,        ONLY: ozprof_start_index, ozprof_end_index,    &
       ozfit_start_index, ozfit_end_index, covar, ozprof_std, ozprof_ap,& 
       which_clima, ncovar, ozprof_apstd, ozprof_init, ozprof, start_layer,     &
       end_layer, eff_alb, eff_alb_init, nlay, nlay_fit, nflay, ptr_order,      &
       aerosol, cloud, ntp, nsfc, atmosprof, ndiv, ps0, pst, useasy, nup2p, nfalb, &
       nalb, tf_fidx, tf_lidx,t_fidx,t_lidx, albidx, albfidx, nt_fit, pos_alb,  &
       the_cfrac, the_cod, the_ctp, the_cld_flg, do_lambcld, which_cld,         &
       ozprof_nstd, maxawin, use_logstate, ozinfo, avg_kernel, contri,  &
       num_iter, smooth_ozbc, colprof, actawin, aerwavs, do_tracewf, mgasprof,  &
       fgasidxs, gasidxs, tracegas, ngas, nflay, do_subfit, osind,     &
       rnind, dcind, isind, irind, slind, shind, nos, nsh, nsl, nrn, ndc, nir,  &
       nis, shfind, osfind, use_tropopause, pst0, nsfc, nfsfc, radcalwrt, do_simu,   &
       tropaod, tropsca, tropwaer, strataod, stratsca, taodind, taodfind, twaeind,   &
       saodfind, ecfrind, ecfrfind, ecodind, ecodfind, ectpind, ectpfind, has_glint, &
       twaefind, saodind, glintprob, which_toz, sprsind, sprsfind, wfcidx, wfcfidx,  &
       nwfc, nfwfc, eff_wfc, eff_wfc_init, so2zind, so2zfind, fit_atanring, &
       use_large_so2_aperr, ozwrtcontri, ozwrtwf, weight_function, contri, &
       trace_profwf, trace_contri, trace_prof, trace_avgk , sacldscl0, &
       which_spres, which_sfct, which_tprof, ncep_fname,div_sun,div_rad ! geun wasp
  
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! =============================
  !  Input / Output Variables
  ! =============================
  INTEGER,            INTENT(IN)  :: initval
  INTEGER,            INTENT(OUT) :: exval
  REAL     (KIND=dp), INTENT(OUT) :: rms
  REAL     (KIND=dp), DIMENSION(3), INTENT(OUT)    :: fitcol
  REAL     (KIND=dp), DIMENSION(3, 2), INTENT(OUT) :: dfitcol  ! smooth+noise, noise

  ! ===============
  ! Local variables
  ! ===============
  INTEGER  :: i, j, nump, errstat, k1, npoints, k, u1idx, u2idx, nsub, nord, is, fidx, lidx
  REAL (KIND=dp) :: asum, ssum, chisq, tmpsa, aodscl, waerscl, salbedo, ncepreso_z0, omi_z0, fdum ! wasp
  REAL (KIND=dp), DIMENSION(n_max_fitpars, n_max_fitpars)    :: bb, sa
  REAL (KIND=dp), DIMENSION(nlay, nlay)                   :: sao3
  REAL (KIND=dp), DIMENSION (n_max_fitpars)               :: lowbond, upbond, fitvar, &
       stderr, fitvarap, stderr1

  REAL (KIND=dp) :: toz, albfc_aperr, albfc_aperr1, albfc_aperr2
  CHARACTER (LEN=6), DIMENSION(n_max_fitpars)             :: varname

  REAL (KIND=dp), DIMENSION(maxlay)                       :: ozprof_std_sav
  REAL (KIND=dp), DIMENSION(n_max_fitpars, n_max_fitpars) :: covar_sav

  REAL (KIND=dp), DIMENSION(nloc)                         :: fine_z0
  REAL (kind=dp), DIMENSION(3)                            :: dfitcol_xa

  LOGICAL, SAVE :: first = .TRUE.
  INTEGER, SAVE :: ozp_fidx,  ozp_lidx, ozf_fidx, ozf_lidx, nf
  
  INTEGER :: nold, nbatm, dum,e  ! geun,wasp
  CHARACTER (LEN=2)               :: monc, dayc   ! geun
  CHARACTER (LEN=4)               :: yrc          ! geun
  LOGICAL                         :: file_exist   ! geun

  !xliu: 09/03/05, add sacldscal, scaling factor for scaling a priori covariance below clouds
  REAL (KIND=dp), DIMENSION(nlay)           :: sacldscl
  INTEGER :: idx_jbak1, idx_jbak2
  INTEGER :: nlat_atm, nlon_atm  !geun 
  REAL (KIND=dp) :: longrid_atm, latgrid_atm !geun

  ! ==============================
  ! Name of this module/subroutine
  ! ==============================

  CHARACTER (LEN=14), PARAMETER :: modulename = 'specfit_ozprof'

  ! Initialize variables for convenience
  IF (first) THEN  ! only need to be initialized once
     nf = n_fitvar_rad;   fitvar_rad = 0.0
     ozp_fidx = ozprof_start_index; ozp_lidx = ozprof_end_index
     ozf_fidx = ozfit_start_index; ozf_lidx = ozfit_end_index

     sa = 0.0; fitvar_rad_apriori = 0.0;   fitvar_rad_std = 0.0;  fitvar_rad_nstd = 0.0
     ozprof_std = 0.0;  ozprof_nstd = 0.0  

     fitvar = 0.0 ; lowbond = 0.0 ; upbond = 0.0 
     varname(1:nf) = fitvar_rad_str(mask_fitvar_rad(1:nf))

  ! tprof files check -----------------------------------------------------------------------------------------
     WRITE(monc, '(I2.2)') the_month          ! from 9 to '09'
     WRITE(dayc, '(I2.2)') the_day            ! from 9 to '09'
     WRITE(yrc,  '(I4.4)') the_year

     IF ( which_tprof == 0 ) THEN  ! geun
        ncep_fname = TRIM(ADJUSTL(atmdbdir)) // 'fnl13.75LST/fnltemp/fnltemp_' // yrc // monc // dayc // '.dat'
        print *, 'TPROF is taken from FNL daily'
     ELSE IF (which_tprof == 1 ) THEN   ! added by geun
        ncep_fname = TRIM(ADJUSTL(atmdbdir)) // 'fnl13.75LST/fnltemp/fnltempavg' // monc // '.dat'
        print *, 'TPROF is taken from FNL monthly'
     ELSE IF (which_tprof == 2 ) THEN   ! added by geun
        ncep_fname = TRIM(ADJUSTL(atmdbdir)) // 'umatmos/umtemp/umtemp_' // yrc // monc // dayc // '.dat'
        print *, 'TPROF is taken from UM daily'
     ENDIF  ! geun

     ! Determine if file exists or not
     INQUIRE (FILE= ncep_fname, EXIST= file_exist)
     IF (.NOT. file_exist) THEN
        WRITE(*, *) 'Warning: no T profile file found, use monthly mean!!!'
        ncep_fname = TRIM(ADJUSTL(atmdbdir)) // 'fnl13.75LST/fnltemp/fnltempavg' // monc // '.dat'
        which_tprof = 1
     ENDIF
  ! -----------------------------------------------------------------------------------------------------------

     first = .FALSE.
  ENDIF

!+---------------------------------------------------------------------------+
! Check radiation input 
!+---------------------------------------------------------------------------+
  !print*, n_irrad_wvl, n_rad_wvl
  !do i=1, n_rad_wvl
    !write(*,'(1x,f13.8,1x,f10.7,1x,f13.8,1x,f13.8)'),curr_rad_spec(wvl_idx,i),checkspec(i),curr_rad_spec(spc_idx,i),&
                                                    !(checkspec(i) - curr_rad_spec(spc_idx,i))/curr_rad_spec(spc_idx,i)
    !print*,curr_rad_spec(wvl_idx,i),curr_rad_spec(spc_idx,i)
  !enddo

  !do i=1,n_irrad_wvl
    !print*,curr_sol_spec(wvl_idx,i),curr_sol_spec(spc_idx,i)
    !print*, curr_sol_spec(wvl_idx,i)/curr_rad_spec(wvl_idx,i),curr_sol_spec(spc_idx,i)/curr_rad_spec(spc_idx,i)
  !enddo
  !stop
!+---------------------------------------------------------------------------+

  ! Initialize variables
  errstat = pge_errstat_ok
  npoints = n_rad_wvl 

  ! use previous fitting results except T, albedo, cloud will be updated
  ! use previous ozone will speed the convergence (could even double)
  fitvar_rad_init = fitvar_rad_saved
  
  ! ===================================================================
  !	         Set up measurement vector and measurement error
  ! ===================================================================


  fitwavs   (1:npoints) = curr_rad_spec(wvl_idx,1:npoints)
  currspec  (1:npoints) = curr_rad_spec(spc_idx,1:npoints)
  fitweights(1:npoints) = curr_rad_spec(sig_idx,1:npoints)
  !print * , npoints !wasp : 204
!do i=1,npoints
  !write(*,*) fitwavs(i),currspec(i) ! wasp: 270~330, normalized, 1e-2~1e0
!enddo
!stop

  IF (ozabs_convl) THEN
     ! For aerosol properties
     actawin = numwin + 2
     IF (actawin > maxawin) STOP 'Need to increase maxawin!!!'
     fidx = 1
     DO i = 1, numwin
        lidx = fidx + nradpix(i) - 1
        aerwavs(i) = fitwavs(fidx)
        fidx = lidx + 1
     ENDDO
     aerwavs(numwin + 1) = fitwavs(lidx)
     i = COUNT(mask = (aerwavs(1:numwin+1) < pos_alb))
     aerwavs(i+2:numwin+2) = aerwavs(i+1:numwin+1)
     aerwavs(i+1) = pos_alb 
     
     ! calculate approximate average wavelength for the window
    
     IF ( weight_rad ) THEN !wasp : F
        asum = SUM ( fitwavs(1:n_rad_wvl) / fitweights(1:n_rad_wvl)**2 )
        ssum = SUM ( 1.D0 / fitweights(1:n_rad_wvl)**2 )
        rad_wav_avg = asum / ssum
     ELSE
        rad_wav_avg = (fitwavs(n_rad_wvl) + fitwavs(1)) / 2.0
     END IF
  ENDIF
  ! ======================================================================
  !       Set up atmospheric cloud properties, albedo and atmosphere
  ! ======================================================================
  ! Spres is provided at NCEP resolution: 2.5 x 2.5
  ! To get the spres for OMI pixel:
  ! 1. get spres at ncep/ncar reso
  ! 2. get z0 at ncep/ncar reso
  ! 3. get z0 at omi spatial resolution
  ! 4. get spres at omi spatial resolution
  IF (which_spres == 0 .OR. which_spres == 1) THEN  ! geun
    nlat_atm = 180  ;  nlon_atm = 360 
    longrid_atm = 1.0  ;  latgrid_atm = 1.0 
  ELSE IF (which_spres == 2) THEN
    nlat_atm = 769  ;  nlon_atm = 1024 
    longrid_atm = 0.351562  ;  latgrid_atm = 0.234375 
  ENDIF  ! geun

  CALL GET_SPRES(the_year, the_month, the_day, the_lon, the_lat, ps0, nlon_atm, nlat_atm, longrid_atm, latgrid_atm)
  CALL get_ncepreso_surfalt(the_lon, the_lat, ncepreso_z0)
  DO i = 1, nloc
     CALL get_finereso_surfalt(the_lons(i), the_lats(i), fine_z0(i))
  ENDDO

  omi_z0 = (SUM(fine_z0(1:4)) + fine_z0(5) * 4.) / 8.

  ! Adjust surface pressure ! dpres (p_omi-p_ncep) = dz(z_omi - z_ncep)
       
  ps0 = ps0 + 1013.25 * (10.**(-omi_z0/16.) - 10.**(-ncepreso_z0/16.))
  the_surfalt = omi_z0
  IF (use_tropopause) THEN 
     IF (which_spres == 0 .OR. which_spres == 1) THEN  ! geun
       nlat_atm = 180  ;  nlon_atm = 360 
       longrid_atm = 1.0  ;  latgrid_atm = 1.0 
     ELSE IF (which_spres == 2) THEN
       nlat_atm = 769  ;  nlon_atm = 1024 
       longrid_atm = 0.351562  ;  latgrid_atm = 0.234375 
     ENDIF  ! geun
     CALL GET_TPRES(the_year, the_month, the_day, the_lon, the_lat, pst, nlon_atm, nlat_atm, longrid_atm, latgrid_atm)
  ELSE
     pst = pst0
  ENDIF
  
  IF (which_toz == 1 ) THEN
     CALL GET_TOZ(the_year, the_month, the_day, the_lon, the_lat, toz)
  !ELSE IF (which_toz == 2) THEN
  !   CALL GET_ZMTOZ(the_year, the_month, the_day, the_lat, toz)
  ELSE
     toz = 0.0
  ENDIF
      
  IF (scnwrt) THEN 
     WRITE(*, '(6(A,F8.2))') '   spres =', ps0, ' tpres =', pst, ' toz = ', toz, ' ctp =', the_ctp
     WRITE(*, '(6(A,F8.3))') '   Lat =', the_lat, 'Lon =', the_lon,' SZA =',the_sza_atm, ' VZA =', the_vza_atm, ' aza = ', the_aza_atm
  ENDIF
  ! ====================================================================
  !	                 Set up atmospheric profiles
  ! ====================================================================
  ! 1. Monthly mean total ozone (EP)
  ! 2. Surface and Tropopause Pressure (also used for getting albedo)
  ! 3. Aerosols (SAGE + GOCART)
  ! 4. Clouds   (GOMECAT)  
  ! 5. Albedo Database (If use clouds)
  ! 6. Temperature (ECMWF)
  ! 7. Ozone, BrO, SO2, NO2, HCHO
  ! 8. Atmos. Profiles
  !    a. TOMS V7 with TOMSV8 temperature profiles (T, P, h, O3) (deleted)
  !    b. GOME workgroup a priori (O3, P, h, T.)  (deleted)
  !    c. ECMWF Temperature and TOMS V8/McPeters Profiles (best conditions)
  ! 9. Get apriori covariance (ozone and non-ozone paramters)
  ! 10.Get initial albedo
  ! ======================================================================
  ! note: the returned ozprof,atmosprof and nup2p are counted bottom up
  !the_sza_atm = 45.
  !the_ctp = 205.819
  !the_cfrac = 1.0
  !the_ctp  = the_ctp + 50 
 
   IF (which_tprof == 0 .OR. which_tprof == 1) THEN ! geun
     nbatm=26 ! nlfnl
     nold=37  ! nlecm(31)+6
   ELSE IF (which_tprof == 2) THEN 
     nbatm=25 ! nlum
     nold=29 ! nlum+4
   ENDIF  !geun
   
   CALL make_atm(the_year, the_month, the_day, ndiv, &
       the_cod, the_cfrac, the_ctp, nlay, toz, ps0, pst, atmosprof(:,0:nlay),    &
       ozprof(1:nlay), nup2p(0:nlay), sacldscl, nbatm, nold, errstat) ! nbatm, nold added : geun 
  IF (errstat == pge_errstat_error)  THEN
     exval = -2; RETURN
  ENDIF

 

  ! second step, use channel 1 retrieval results
  !IF (curr_rad_spec(1, 1) > 307.2) THEN
  !   CALL GET_FIRST_RETRIEVAL(the_month, the_lat, nlay, atmosprof(1, 0:nlay), &
  !        atmosprof(2, 0:nlay), atmosprof(2, 0:nlay), ozprof(1:nlay), sao3)
  !   STOP 'NO Two-Step Retrieval!!!'
  !ENDIF

  ! Add ozone and trace gases into initialized array for the first retrieval
  ! Then use the previous retrieval as the initial
  ! NO2 and HCHO: always using GEOS-CHEM fields with 100% error
  ! For BrO:      a prioir from model fields but with 1.0E-14 error globally (enough information)
  ! For SO2:      a priori from model fields but with dynamic a priori error 
  !               to deal with volcanic eruption (implemented in ozone_reverse.f90)
  IF (initval == 0 .OR. ANY(fitvar_rad_init(ozp_fidx:ozp_lidx) <= 0.0))  THEN  !geun  off -> fix
     fitvar_rad_init(ozp_fidx:ozp_lidx) = ozprof(1:nlay)  
  ENDIF !geun
  IF (nsfc < nlay) fitvar_rad_init(ozp_lidx + 1 - nlay + nsfc : ozp_lidx) = ozprof(nsfc+1:nlay)
  fitvar_rad_apriori(ozp_fidx:ozp_lidx) = ozprof(1:nlay)
  
  ! Always use climatological temperature profiles
  fitvar_rad_init(t_fidx:t_lidx) = &
       (atmosprof(3, 1:nlay) + atmosprof(3, 0:nlay-1)) / 2.0
   
  ! For aerosols, a priori is based on model, a priori error is assumed as 100%, 
  ! and initial value is from previous retrievals.
  fitvar_rad_apriori(taodind) = tropaod(actawin)
  fitvar_rad_apriori(twaeind) = tropwaer(actawin)
  fitvar_rad_apriori(saodind) = strataod(actawin)  
  IF (initval == 0 .OR. taodfind == 0) fitvar_rad_init(taodind) = tropaod(actawin)
  IF (initval == 0 .OR. twaefind == 0) fitvar_rad_init(twaeind) = tropwaer(actawin)
  IF (initval == 0 .OR. saodfind == 0) fitvar_rad_init(saodind) = strataod(actawin)

  ! Set up albedo and cloud fraction in the retrieval
  ! albedo and cloud fraction can be adjusted based on 370.2 nm reflectance
   
  CALL SET_CLDALB(npoints, fitwavs, the_cod, the_ctp, the_cfrac, salbedo, errstat)
  IF (errstat == pge_errstat_error)  THEN
     exval = -1; RETURN
  ENDIF

  ! For clouds, initial ctp, cod is based on assumed input (e.g., 20/10) or from other products, 
  ! which maybe re-adjusted using longer wavelengths
  fitvar_rad_init(ecfrind) = the_cfrac; fitvar_rad_apriori(ecfrind) = the_cfrac
  fitvar_rad_init(ecodind) = the_cod;   fitvar_rad_apriori(ecodind) = the_cod
  fitvar_rad_init(ectpind) = the_ctp;   fitvar_rad_apriori(ectpind) = the_ctp
  IF (so2zfind > 0) THEN
     fitvar_rad_init(so2zind) = fitvar_rad_init_saved(so2zind)
     lo_radbnd(so2zind) = the_surfalt; up_radbnd(so2zind) = 30.
     fitvar_rad_apriori(so2zind)  = fitvar_rad_init(so2zind)
  ENDIF
  
  ! For surface pressure
  fitvar_rad_init(sprsind) = ps0;       fitvar_rad_apriori(sprsind) = ps0 
  IF (sprsfind > 0) THEN
     lo_radbnd(sprsind) = atmosprof(1, nlay) - (atmosprof(1, nlay)-atmosprof(1, nlay-1)) * 0.5
     up_radbnd(sprsind) = atmosprof(1, nlay) + (atmosprof(1, nlay)-atmosprof(1, nlay-1)) * 0.5
  ENDIF

  DO k = 1, ngas
     i = fgasidxs(k)
     IF (i > 0) THEN
        j = mask_fitvar_rad(i)          
        
        ! tracegas(k, 8) = 1.0 - the_cfrac * tracegas(k, 8)  
        ! fitvar_rad_apriori(j) = mgasprof(k, nflay + 1) * refspec_norm(gasidxs(k)) !* tracegas(k, 8)
        ! This is incorrect, since everything is taken into account implicitly through weighting function       
        fitvar_rad_apriori(j) = mgasprof(k, nflay + 1) * refspec_norm(gasidxs(k)) 
        !IF ( npix_fitted == 0 .OR. gasidxs(k) == so2_idx .OR. gasidxs(k) == so2v_idx &
        !     .OR. initval == 0 .OR. fitvar_rad_init(j) < 0.0) &

        fitvar_rad_init(j) =  fitvar_rad_apriori(j) 
        IF (gasidxs(k) == so2_idx .OR. gasidxs(k) == so2v_idx) fitvar_rad_init(j) = 0.0 
        IF (gasidxs(k) == bro2_idx) fitvar_rad_apriori(j) = fitvar_rad_apriori(j) * 5. / 2.
        
        ! initial values for trace gas shift parameters needs to be fixed for every pixel
        ! otherwise increasing from North to South to unreasonably large
        j = shift_offset + gasidxs(k)
        fitvar_rad_init(j) = 0.00
          
     ENDIF
  
  ENDDO

  IF ( initval == 0 ) THEN
     DO i = 1, numwin
        fitvar_rad_init (osind(i, 1:maxoth)) = 0.0D0
        fitvar_rad_init (shind(i, 1:maxoth)) = 0.0D0
     ENDDO
  ENDIF

  DO i = 1, numwin
     fitvar_rad_init (osind(i, 1:maxoth)) = 0.0D0
     fitvar_rad_init (shind(i, 1:maxoth)) = 0.0D0
  ENDDO

  DO i = albidx, albidx + nalb - 1
     IF (fitvar_rad_str(i)(4:4) /= '0') fitvar_rad_init(i) = 0.D0
  ENDDO
  fitvar_rad_init(irind(1, 1))    = -1.0E-5    ! non zero
  !fitvar_rad_init(irind(1, 1)+1) = 1.0e-7    !geun  on : in0 uv2 fix
  !fitvar_rad_init(rnind(1, 1))   = -1.87    !geun  on ring fix
  !fitvar_rad_init(rnind(1, 1)+1) = -1.87    !geun  on ring fix
 
  IF (do_subfit) THEN
     nsub = numwin
  ELSE
     nsub = 1
  ENDIF
             
  ! ======================================================================
  !	 Set up state vector, a priori state vector and covariance matrix
  ! ======================================================================
  ! set up fitting variables
  fitvar_rad = fitvar_rad_init    
  IF (start_layer /= 1) THEN
     lo_radbnd(ozp_fidx:ozp_fidx+start_layer-2) = &
          fitvar_rad(ozp_fidx:ozp_fidx+start_layer-2)
     up_radbnd(ozp_fidx:ozp_fidx+start_layer-2) = &
          fitvar_rad(ozp_fidx:ozp_fidx+start_layer-2)     
  ENDIF
    
  IF (end_layer /= nlay) THEN
     i = nlay - end_layer - 1
     lo_radbnd(ozp_lidx-i:ozp_lidx) = fitvar_rad(ozp_lidx-i:ozp_lidx)
     up_radbnd(ozp_lidx-i:ozp_lidx) = fitvar_rad(ozp_lidx-i:ozp_lidx)    
  ENDIF
   
  ! Get a priori ozone covariance matrix
  CALL GET_APRIORI_COVAR(toz, ozprof(1:nlay), sao3)
  IF (nsfc < nlay) THEN
        sao3(nsfc+1:nlay, :) = 0.0; sao3(:, nsfc+1:nlay) = 0.0
  ENDIF     

  !xliu, 08/29/05 scaling a priori for layers below clouds to avoid smoothing even
  !for full cloudy conditions
     sacldscl0 = sacldscl ! used in update_o3_sao3
  IF (.NOT. smooth_ozbc) THEN
     sacldscl = sacldscl * the_cfrac + (1.0 - the_cfrac )
     DO i = 1, nsfc
        IF (sacldscl(i) < 1.0) THEN 
           tmpsa = sao3(i, i)
           sao3(i, 1:nlay) = sao3(i, 1:nlay) * sacldscl(i)
           sao3(1:nlay, i) = sao3(1:nlay, i) * sacldscl(i)
           sao3(i, i) = tmpsa
        ENDIF
     ENDDO
  ENDIF
  
  ! Set up a priori state vector and covariance matrix
  ! use a priori for O3, T, 0th albedo, trace gas, Ring effect
  ! Zero for others

     fitvar_rad_apriori(ozp_fidx:ozp_lidx) = ozprof(1:nlay)
     fitvar_rad_apriori(t_fidx:t_lidx)     = fitvar_rad(t_fidx:t_lidx)
     DO i = albidx, albidx + nalb - 1
        IF (fitvar_rad_str(i)(4:4) == '0') fitvar_rad_apriori(i) = fitvar_rad(i)
     ENDDO 
     DO i = wfcidx, wfcidx + nwfc - 1
        IF (fitvar_rad_str(i)(4:4) == '0') fitvar_rad_apriori(i) = fitvar_rad(i)
     ENDDO 
     IF (fit_atanring) THEN        
        fitvar_rad_apriori(rnind(1, 1) : rnind(1, 3) + nsub - 1) = &
             fitvar_rad(rnind(1, 1) : rnind(1, 3) + nsub - 1) 
     ELSE
        fitvar_rad_apriori(rnind(1, 1) : rnind(1, 1) + nsub - 1) = &
             fitvar_rad(rnind(1, 1) : rnind(1, 1) + nsub - 1) 
     ENDIF
     IF (comfidx > 0) fitvar_rad_apriori(comvidx) = 1.0
     IF (cm1fidx > 0) fitvar_rad_apriori(cm1vidx) = 1.0

     u1idx = max_calfit_idx + (us1_idx - 1) * mxs_idx 
     u2idx = max_calfit_idx + (us2_idx - 1) * mxs_idx 
     ! The following covariance matrix are fixed
     IF ( npix_fitted == 0) THEN       

        DO i = 1, nf 
           j = mask_fitvar_rad(i)
           k = j - shift_offset
           
           ! 2500% error for parameters unless specified  (xliu, 03/21/2006)
           IF (k > 0 .AND. k < max_rs_idx) THEN  ! shift parameters specified in BOREAS.inp
              sa(i, i) = 4.0E-1
           ELSE IF (j > u1idx .AND. j <= u1idx + mxs_idx) THEN
              sa(i, i) = 0.5
           ELSE IF (j > u2idx .AND. j <= u2idx + mxs_idx) THEN
              sa(i, i) = 0.5
               
           ELSE IF (j >= rnind(1, 1) .AND. j <= rnind(nsub, 1)) THEN
              IF (fit_atanring) THEN 
                 sa(i, i) = 0.5
              ELSE
                 sa(i, i) = 1.0
              ENDIF
           ELSE IF (j >= rnind(1, 2) .AND. j <= rnind(nsub, 2)) THEN
              IF (fit_atanring) THEN 
                 sa(i, i) = 9.
              ELSE
                 sa(i, i) = 2.0E-3
              ENDIF
         
           ELSE IF (j >= rnind(1, 3) .AND. j <= rnind(nsub, 3) ) THEN
              IF (fit_atanring) THEN 
                 sa(i, i) = 4.
              ELSE
                 sa(i, i) = 1.0E-4
              ENDIF
           
           ELSE IF (j >= rnind(1, 4) .AND. j <= rnind(nsub, 4) ) THEN
              sa(i, i) = 5.0E-5   
           ELSE IF (j >= isind(1, 1) .AND. j <= isind(nsub, 1)) THEN
              sa(i, i) = 1.0E-1 
           ELSE IF (j >= isind(2, 1) .AND. j <= isind(numwin, maxoth)) THEN
              sa(i, i) = 1.0E-2 
           ELSE IF (j >= irind(1, 1) .AND. j <= irind(nsub, 1)) THEN
              sa(i, i) = 1.0E-8 
           ELSE IF (j >= irind(2, 1) .AND. j <= irind(numwin, maxoth)) THEN
              sa(i, i) = 1.0E-10
           ELSE IF (j >= dcind(1, 1) .AND. j <= dcind(nsub, 1) ) THEN
              sa(i, i) = 0.1  
           ELSE IF (j >= dcind(1, 2) .AND. j <= dcind(nsub, maxoth)) THEN
              sa(i, i) = 0.02 
           ELSE IF (i == comfidx .OR. i == cm1fidx) THEN
              sa(i, i) = 1.0 
           ELSE IF (i < ozf_fidx .OR. i > ozf_lidx) THEN
              sa(i, i) = (fitvar_rad(j))**2.0 * 25.0  
           ENDIF
        ENDDO
        
        IF (nt_fit > 0) THEN   ! 5 K std. deviation
           DO i = tf_fidx, tf_lidx
              sa(i, i) =  25.0
           ENDDO
        ENDIF

        ! zero-order error +/-0.05 nm
        ! error decreases by a factor of 10 when the order increases by 1
        IF (nos > 0) THEN
           DO i = 1, nos
              DO j = 1, nsub
                 k = osfind(j, i)
                 IF (k > 0) sa(k, k) = 5.0E-4 * (10.0 ** (-(i - 1) * 2.0)) 
              ENDDO
           ENDDO
        ENDIF
        
        ! zero-order error +/-0.01 nm
        ! error decreases by a factor of 10 when the order increases by 1
        IF (nsh > 0) THEN
           DO i = 1, nsh
              DO j = 1, nsub
                 k = shfind(j, i)
                 IF (k > 0) sa(k, k) = 5.0E-4  * (10.0 ** (-(i - 1) * 2.0))
              ENDDO
           ENDDO
        ENDIF        
     ENDIF

     ! need to update o3 a priori covariance matrix for each retrieval
     sa(ozf_fidx:ozf_lidx, ozf_fidx:ozf_lidx) = &
          sao3(start_layer:end_layer,start_layer:end_layer) 
!print*, nfalb >0 .OR. nfwfc >0 .OR. ecfrfind > 0
     IF (nfalb > 0 .OR. nfwfc > 0 .OR. ecfrfind > 0) THEN
!print*, albfc_aperr
        albfc_aperr = 0.05**2.0
!print*, albfc_aperr
!print*, tropaod(1) >= 0.25 .AND. taodfind == 0 .AND. twaefind == 0
!print*, has_glint
       IF (tropaod(1) >= 0.25 .AND. taodfind == 0 .AND. twaefind == 0) THEN
           albfc_aperr1 = albfc_aperr * (1.0 + 8.0 * tropaod(1) - 2.0) 
        ELSE IF (has_glint) THEN  ! Assume a priori error of 0.2 instead of 0.05 for 100% sun glint
           albfc_aperr1 = albfc_aperr * ( 1.0 + 15.0 * glintprob) 
        ELSE
           albfc_aperr1 = albfc_aperr 
        ENDIF
        albfc_aperr2    = albfc_aperr * (1.0 + 15.0 * salbedo)
!print*, albfc_aperr1,albfc_aperr2,salbedo
!print*,'end check'
        albfc_aperr     = MAX(albfc_aperr1, albfc_aperr2)
     ENDIF

     IF (nfalb > 0) THEN
        DO i = albfidx, albfidx + nfalb - 1
           ! The a priori std. for non-zero albedo terms are based on retrievals 
           ! (1.6E-3, 1.0E-5, 1.0E-7, 1.0E-9)
           ! (3.0E-2, 4.0E-4, 1.6E-5, 6.4E-7)
           READ (fitvar_rad_str(mask_fitvar_rad(i))(4:4), '(I1.1)') nord
           sa(i, i) = albfc_aperr * (5.0 ** ( - nord * 2.0))
           IF (salbedo > 0.6 .AND. nord >= 1) sa(i, i) = 0.0
        ENDDO
     ENDIF

     IF (nfwfc > 0) THEN
        DO i = wfcfidx, wfcfidx + nfwfc - 1
           READ (fitvar_rad_str(mask_fitvar_rad(i))(4:4), '(I1.1)') nord
           sa(i, i) = albfc_aperr * (5.0 ** ( - nord * 2.0))
        ENDDO
     ENDIF
     
     ! Use 50% (NO2 and HCHO) or 100% (SO2 and BrO) error for other minor trace gases
     ! The apiori of these trace gases are determined by climatology (NO2, HCHO)
   
     DO i = 1, ngas
        j = fgasidxs(i)
        IF (j > 0) THEN
           sa(j, j) = (fitvar_rad_apriori(mask_fitvar_rad(j)))**2.0
           IF (gasidxs(i) == hcho_idx .OR. gasidxs(i) == no2_t1_idx ) sa(j, j) = sa(j, j) * 0.25  !50% error
           IF (gasidxs(i) == o2o2_idx) sa(j, j) = sa(j, j) * 0.25       !50% error
           IF ((gasidxs(i) == so2_idx .OR. gasidxs(i) == so2v_idx) .AND. use_large_so2_aperr)  &
                sa(j, j) = (1.0E17 * refspec_norm(gasidxs(i))) ** 2.0
        ENDIF
     ENDDO

     ! A priori covariance matrix for aerosols and clouds
     IF (taodfind > 0) THEN  ! 100%
        sa(taodfind, taodfind) = fitvar_rad_apriori(taodind)**2   !* 0.25
     ENDIF
     IF (saodfind > 0) THEN  ! 50%
        sa(saodfind, saodfind) = fitvar_rad_apriori(saodind)**2 * 0.25
     ENDIF
     IF (twaefind > 0) THEN  ! aerosol single scattering albedo change by 0.05
        sa(twaefind, twaefind) = 2.5E-3
     ENDIF
     IF (ecfrfind > 0) THEN  
        sa(ecfrfind, ecfrfind) = albfc_aperr
     ENDIF
     IF (ecodfind > 0) THEN  ! 10%
        sa(ecodfind, ecodfind) = (fitvar_rad_apriori(ecodind) * 0.1)**2
     ENDIF
     IF (ectpfind > 0) THEN  ! 100 mb
        sa(ectpfind, ectpfind) = 100.**2
     ENDIF    

     IF (sprsfind > 0) THEN  ! 100 mb
        sa(sprsfind, sprsfind) = 30.**2
     ENDIF

     IF (so2zfind > 0) THEN  ! 1.0 km
        sa(so2zfind, so2zfind) = 1.0 ** 2.
     ENDIF
        
    
  ! Initialize exit value to be zero (missing)
  exval = 0  
    
  ! Create a condensed array of fitting variables that are varied. 
  ! This considerably reduces the execution time of the fitting routine.
  IF (radcalwrt .AND. do_simu) fitvar_rad = fitvar_rad_apriori ! F - wasp!
  fitvar(1:nf)    = fitvar_rad(mask_fitvar_rad(1:nf))
  fitvarap(1:nf)  = fitvar_rad_apriori(mask_fitvar_rad(1:nf))
  lowbond(1:nf)   = lo_radbnd(mask_fitvar_rad(1:nf))
  upbond(1:nf)    = up_radbnd(mask_fitvar_rad(1:nf))
 
  !WRITE(*, '(A)') 'Initial Guess: '
  !WRITE(*, '(12F8.3)') fitvar_rad(ozp_fidx:ozp_lidx), SUM(fitvar_rad(ozp_fidx:ozp_lidx))
  IF  ( reduce_resolution) then ! Jbak's modification to avoide errors in interpolating solar spectrum with the fitwavs
      do i = 1, n_refwvl-1  
    if ( refwvl(i) > refwvl(i+1)) then
    
     print *, refwvl(i), 'decrease', refwvl(i+1)
 refwvl(i) = refwvl(i+1) - (refwvl(i+2)-refwvl(i+1)) 
     endif
   enddo 
  ENDIF 

  CALL ozprof_inverse (nf, varname(1:nf), fitvar(1:nf), fitvarap(1:nf), &
       lowbond(1:nf), upbond(1:nf), npoints, nump, sa(1:nf,1:nf), bb(1:nf,1:nf), &
       chisq, fitspec_rad(1:npoints), fitres_rad(1:npoints), exval)
  fitvar_rad(mask_fitvar_rad(1:nf)) = fitvar(1:nf)  ! for safe 
  fitvar_rad_apriori(mask_fitvar_rad(1:nf)) = fitvarap(1:nf)  ! Some a priori values can be changed
  DO i = 1, nf
     fitvar_rad_aperror(i) = SQRT(sa(i, i))
  ENDDO

  !WRITE(*, '(A)') 'Retrieval: '
  !WRITE(*, '(12F8.3)') fitvar_rad(ozp_fidx:ozp_lidx), SUM(fitvar_rad(ozp_fidx:ozp_lidx))
  
  ! reset the wavelength shifts to zero if retrievals are not successful
  ! because the failure of retrievals are due to too large wavelength shifts most of the time
  IF (exval <= 0) THEN
     DO i = 1, numwin
        fitvar_rad_init(osind(i, 1:maxoth)) = 0.0;   fitvar_rad_init(shind(i, 1:maxoth)) = 0.0
     ENDDO
  ENDIF
    
  IF (exval < 0) THEN     ! Terminate the whole retrieval
     fitcol  = missing_value_dp
     dfitcol = missing_value_dp    
     fitvar_rad_saved = fitvar_rad_init
     RETURN
  ELSE
     IF (exval == 0 ) THEN
        fitvar_rad_saved = fitvar_rad_init 
     ELSE
        fitvar_rad_saved = fitvar_rad 
     ENDIF
  ENDIF
  
  ! calculate rms difference between measurements and calculations
  ! If using measurement error, rms is ideally to be 1, if it > 1, suggesting
  ! the estiamted measurement error is too small, and vice versa too large
  rms = SQRT(chisq / REAL(npoints, KIND=dp))
  
  ! need to multiply rms to get the actual retrieval random noise error
  ! no matter whether use measurement error or not (sig=1.0)
  ! because measurement error is not reliable, a rms of > 1 indicates
  ! underestimated measurement error and vice versa
  ! but for smoothing error, do not do it
  ! covar(1:nf,1:nf)=rms**2.0*covar(1:nf, 1:nf) !*npoints/(npoints-nf)
  
  ozprof(1:nlay) = fitvar_rad(ozp_fidx:ozp_lidx)
  DO i = 1, nf
     stderr(i) = SQRT(covar(i, i)); stderr1(i) = SQRT(ncovar(i, i))
  END DO


!print*,ecfrfind,ecfrind,n_max_fitpars,nf !,sa(ecfrfind,ecfrfind),albfc_aperr
!print*,'fitvar_rad(ecfrind),fitvar_rad_init(ecfrind),fitvar_rad_apriori(ecfrind)'
!print*,fitvar_rad(ecfrind),fitvar_rad_init(ecfrind),fitvar_rad_apriori(ecfrind)
!print*,'sa(ecfrfind,ecfrfind),albfc_aperr'
!print*,sa(ecfrfind,ecfrfind),albfc_aperr
!print*,'stderr(1:nf)'
!print*,stderr(1:nf)
!print*,'stderr1(1:nf)'
!print*,stderr1(1:nf)
 !stop 

  fitvar_rad_std(mask_fitvar_rad(1:nf))  = stderr(1:nf)
  fitvar_rad_nstd(mask_fitvar_rad(1:nf)) = stderr1(1:nf)
  
  ! compute error in total ozone, stratospheric ozone, tropospheric ozone
  ! -------------------------------------------------------------
  ! FITCOL:       Used for total o3, strat o3, and trop o3
  ! DFITCOL:      Uncertainty for corresponding columns
  ! NLAY_FIT:     Number of ozone layers that are varied
  ! STATRT_LAYER: first layer that are varied among all layers
  ! END_LAYER:    last layer that are varied among all layers
  ! -------------------------------------------------------------
  fitcol = 0.0  ;  dfitcol = 0.0 ; dfitcol_xa = 0.0
  fitcol(1) = SUM (ozprof(1:nsfc)) 
  j = start_layer 
  DO i = ozf_fidx, ozf_lidx
     ozprof_std(j) = stderr(i); ozprof_nstd(j) = stderr1(i); j = j + 1
  END DO

  k1 = ozf_fidx - start_layer
  DO is = 1, 2
     IF (is == 1 ) THEN
        ozprof_std_sav(ozf_fidx:ozf_lidx) = ozprof_std(ozf_fidx:ozf_lidx)
        covar_sav(ozf_fidx:ozf_lidx, ozf_fidx:ozf_lidx) = covar(ozf_fidx:ozf_lidx, ozf_fidx:ozf_lidx)
     ELSE
        ozprof_std_sav(ozf_fidx:ozf_lidx) = ozprof_nstd(ozf_fidx:ozf_lidx)
        covar_sav(ozf_fidx:ozf_lidx, ozf_fidx:ozf_lidx) = ncovar(ozf_fidx:ozf_lidx, ozf_fidx:ozf_lidx)
     ENDIF
            
     ! remove correlated error to get error in total ozone
     DO i = start_layer, nsfc
        IF (ozprof_std_sav(i) > 0.0) THEN
          
           dfitcol(1, is) = dfitcol(1, is) + ozprof_std_sav(i)  ** 2.0
!           IF (is ==1) dfitcol_xa(1) =dfitcol_xa(1) + sao3(i,i)
           DO j = start_layer, i - 1
                 
              dfitcol(1, is) = dfitcol(1, is) + 2.0 * ozprof_std_sav(i) * ozprof_std_sav(j) * &
                    covar_sav(k1+i, k1+j) / SQRT(covar_sav(k1+i, k1+i) * covar_sav(k1+j, k1+j))          
!              if (is ==1) dfitcol_xa(1) = dfitcol_xa(1) + 2.0 * sqrt(sao3(i,i)) * sqrt(sao3(j,j)) * &
!                    sao3(k1+i, k1+j) / SQRT(sao3(k1+i, k1+i) * sao3(k1+j, k1+j))          
           END DO
            
        ENDIF
     ENDDO
     IF (ntp > 0) THEN     
        ! stratospheric ozone
        DO i = start_layer, ntp
           IF (ozprof_std_sav(i) > 0.0) THEN
              dfitcol(2, is) = dfitcol(2, is) + ozprof_std_sav(i)**2.0
              
!           IF (is ==1) dfitcol_xa(2) =dfitcol_xa(2) + sao3(i,i)
              DO j = start_layer, i - 1
                 dfitcol(2, is) = dfitcol(2, is) + 2.0 * ozprof_std_sav(i) * ozprof_std_sav(j) * &
                      covar_sav(k1+i, k1+j) / SQRT(covar_sav(k1+i, k1+i) * covar_sav(k1+j, k1+j))
              
!              if (is ==1) dfitcol_xa(2) = dfitcol_xa(2) + 2.0 * sqrt(sao3(i,i)) * sqrt(sao3(j,j)) * &
!                    sao3(k1+i, k1+j) / SQRT(sao3(k1+i, k1+i) * sao3(k1+j, k1+j))          
              ENDDO
           ENDIF
        ENDDO

        ! tropospheric ozone
        DO i = ntp + 1, nsfc
           IF (ozprof_std_sav(i) > 0.0) THEN
              dfitcol(3, is) = dfitcol(3, is) + ozprof_std_sav(i)**2.0

  !         IF (is ==1) dfitcol_xa(3) =dfitcol_xa(3) + sao3(i,i)
              DO j = ntp + 1, i - 1
                 dfitcol(3, is) = dfitcol(3, is) + 2.0 * ozprof_std_sav(i) * ozprof_std_sav(j) *  &
                      covar_sav(k1+i, k1+j) / SQRT(covar_sav(k1+i, k1+i) * covar_sav(k1+j, k1+j))

 !             if (is ==1) dfitcol_xa(3) = dfitcol_xa(3) + 2.0 * sqrt(sao3(i,i)) * sqrt(sao3(j,j)) * &
 !                   sao3(k1+i, k1+j) / SQRT(sao3(k1+i, k1+i) * sao3(k1+j, k1+j))          
                 !print *, ozprof_std_sav(i), ozprof_std_sav(j), sqrt(covar_sav(k1+i, k1+i)), sqrt(covar_sav(k1+j, k1+j))
              ENDDO
           ENDIF
        ENDDO
     ENDIF
  ENDDO

!            print * , sqrt(dfitcol_xa(:))
  dfitcol(1, :)    = SQRT (dfitcol(1, :))
  IF (ntp > 0) THEN
     fitcol(2)     = SUM (ozprof(1:ntp)) 
     fitcol(3)     = SUM (ozprof(ntp+1:nsfc))
     dfitcol(2, :) = SQRT (dfitcol(2, :))
     dfitcol(3, :) = SQRT (dfitcol(3, :))
  ENDIF

  IF (nwfc > 0) THEN
     DO i = wfcidx + nwfc - 1, wfcidx, -1
        IF (fitvar_rad_str(i)(4:4) == '0') EXIT
     ENDDO
     the_cfrac = fitvar_rad(i)   ! Use UV2/last channel cloud fraction
  ELSE IF ( ecfrfind > 0) THEN
     the_cfrac = fitvar_rad(ecfrind )
  ENDIF
  
  !IF (the_cfrac <= 1.0D-3 .AND. ecfrind > 0) THEN
  !   the_cod = 0.d0; the_ctp = 0.d0
  !ENDIF
  
  ozprof(1:nlay) = fitvar_rad(ozp_fidx:ozp_lidx)
  ozprof_init(1:nlay) = fitvar_rad_init(ozp_fidx:ozp_lidx)
  ozprof_ap(1:nlay) = fitvar_rad_apriori(ozp_fidx:ozp_lidx) 
  ozprof_apstd = 0.0
  j = start_layer
  IF (use_logstate) THEN
     DO i = ozf_fidx, ozf_lidx
        ozprof_apstd(j) = SQRT(sa(i, i))*ozprof_ap(j); j = j + 1
     ENDDO
  ELSE
     DO i = ozf_fidx, ozf_lidx
        ozprof_apstd(j) = SQRT(sa(i, i)); j = j + 1
     ENDDO
  ENDIF
     

  IF (nsfc < nlay) THEN
     ozprof(nsfc+1:nlay) = -99.99;     ozprof_ap(nsfc+1:nlay) = -99.99
  ENDIF

  eff_alb_init = fitvar_rad_init(albidx:albidx+maxalb-1)
  eff_alb      = fitvar_rad(albidx:albidx+maxalb-1) 
  eff_wfc_init = fitvar_rad_init(wfcidx:wfcidx+maxwfc-1)
  eff_wfc      = fitvar_rad(wfcidx:wfcidx+maxwfc-1) 

 ! trace gases
 DO k = 1, ngas
    i = fgasidxs(k)
    IF (i > 0) THEN
       j = mask_fitvar_rad(i)
       tracegas(k, 1) = fitvar_rad_init(j)
       tracegas(k, 2) = fitvar_rad_apriori(j)
       tracegas(k, 3) = SQRT(sa(i, i))
       tracegas(k, 4) = fitvar_rad(j)
       tracegas(k, 5) = fitvar_rad_std(j) 
       tracegas(k, 6) = fitvar_rad_nstd(j) 
       tracegas(k, 1:6) = tracegas(k, 1:6) / refspec_norm(gasidxs(k))  
       fitvar_rad(j) = tracegas(k, 4)
       fitvar_rad_std(j) = tracegas(k, 5); fitvar_rad_nstd(j) = tracegas(k, 6)

       ! xliu: 07/01/2010, 08/09/2010
       ! Change weighting function and contribution function for trace gas variables 
       ! wrt to the reported unit instead of normalized quantities
       IF (ozwrtwf) THEN
          weight_function(1:npoints, i) = weight_function(1:npoints, i) * refspec_norm(gasidxs(k)) 
          trace_profwf(k, 1:npoints, 1:nlay) = trace_profwf(k, 1:npoints, 1:nlay) * refspec_norm(gasidxs(k))
       ENDIF
       IF (ozwrtcontri) THEN
          contri(i, 1:npoints) = contri(i, 1:npoints) / refspec_norm(gasidxs(k))  
          trace_contri(k, 1:npoints) = contri(i, 1:npoints)
       ENDIF

       ! Add trace gas profile (a priori profile/shape)
       DO j = 1, nlay
          fidx = nup2p(j - 1) + 1; lidx = nup2p(j)
          trace_prof(k, j) = SUM(mgasprof(k, fidx:lidx))
       ENDDO    
    ENDIF
 ENDDO

 IF ( taodfind > 0) THEN
    aodscl = fitvar_rad(taodind) / tropaod(actawin)
    tropaod(1:actawin) = tropaod(1:actawin) * aodscl
    
    IF ( twaefind == 0 ) THEN  ! Single scattering albedo does not change
       tropsca(1:actawin) = tropsca(1:actawin) * aodscl
    ENDIF
 ELSE
    aodscl = 1.0
 ENDIF

IF ( twaefind > 0 ) THEN
    waerscl = fitvar_rad(twaeind) / tropwaer(actawin)   ! Scale single scattering albedo
    tropwaer(1:actawin) = tropwaer(1:actawin) * waerscl
    tropsca(1:actawin)  = tropsca(1:actawin) * waerscl * aodscl
 ENDIF

 IF ( saodfind > 0 ) THEN
    aodscl = fitvar_rad(saodind) / strataod(actawin)
    strataod(1:actawin) = strataod(1:actawin) * aodscl    
    stratsca(1:actawin) = stratsca(1:actawin) * aodscl
 ENDIF

 IF ( ecfrfind > 0) THEN
    the_cfrac = fitvar_rad(ecfrind)
 ENDIF

 IF (sprsfind > 0) THEN
    atmosprof(1, nsfc) = fitvar_rad(sprsind)
 ENDIF


 RETURN
END SUBROUTINE specfit_ozprof
