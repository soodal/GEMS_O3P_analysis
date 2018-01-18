SUBROUTINE gems_o3p_adj_solar_data (pge_error_status)

  USE OMSAO_precision_module
  USE OMSAO_indices_module,         ONLY: wvl_idx, spc_idx, sig_idx, hwe_idx
  USE OMSAO_parameters_module,      ONLY: normweight, mrefl
  USE OMSAO_variables_module,       ONLY: curr_sol_spec, n_irrad_wvl, use_meas_sig, &
       numwin, nsol_ring, sol_spec_ring, nsolpix, yn_varyslit, slit_rad, solwinfit, &
       nslit, slitwav, slitfit, sring_fidx, sring_lidx,  currpix, slitdis, which_slit
  USE GEMS_O3P_gemsdata_module,         ONLY: gems_rad, gems_irrad, &
                                              gems_refl,gems_ring, gems_solcal
                                            ! omi_nslit_rad, omi_slitwav_rad, omi_radslitfit, omi_nslit_sol,  &
                                            ! omi_slitwav_sol, omi_solslitfit
  USE ozprof_data_module,           ONLY: div_sun, sun_posr, sun_specr, nrefl
  USE OMSAO_errstat_module 

  IMPLICIT NONE
  
  ! =================
  ! Output variables
  ! =================
  INTEGER, INTENT (OUT) :: pge_error_status

  INTEGER :: i

  pge_error_status = pge_errstat_ok

  ! Solar Spectrum
  n_irrad_wvl = gems_irrad%nwav(currpix) 
  div_sun     = gems_irrad%norm(currpix)

  curr_sol_spec(wvl_idx, 1:n_irrad_wvl) = gems_irrad%wavl(1:n_irrad_wvl, currpix) 
  curr_sol_spec(spc_idx, 1:n_irrad_wvl) = gems_irrad%spec(1:n_irrad_wvl, currpix) 
  IF (use_meas_sig ) THEN
     curr_sol_spec(sig_idx, 1:n_irrad_wvl) = gems_irrad%prec(1:n_irrad_wvl, currpix)
  ELSE
     curr_sol_spec(sig_idx, 1:n_irrad_wvl) = normweight
  ENDIF
  nsolpix(1:numwin) = gems_irrad%npix(1:numwin, currpix) 

  ! Solar Spectrum for Ring Calculation
  nsol_ring = gems_ring%nsol(currpix)
  sol_spec_ring(1, 1:nsol_ring) = gems_ring%solwavl(1:nsol_ring, currpix)
  sol_spec_ring(2, 1:nsol_ring) = gems_ring%solspec(1:nsol_ring, currpix)
  sring_fidx = gems_ring%sol_lin(currpix)
  sring_lidx = gems_ring%sol_uin(currpix)
  
  !WRITE(90, '(F12.6, D14.6)') (sol_spec_ring(1:2, i), i = 1, nsol_ring)

  ! Reflectance spectrum at ~370 nm
  sun_posr (1:nrefl) = gems_refl%solwavl(1:nrefl, currpix)
  sun_specr(1:nrefl) = gems_refl%solspec(1:nrefl, currpix)

  ! Load slit calibration parameters
  IF (which_slit < 4) THEN
     IF (yn_varyslit) THEN
        IF (slit_rad) THEN
           !nslit = omi_nslit_rad(currpix)
           !slitwav(1:nslit) = omi_slitwav_rad(1:nslit, currpix)
           !slitfit = omi_radslitfit(:, :, :, currpix)
        ELSE
           !nslit = gems_solcal%nslit(currpix)
           !slitwav(1:nslit) = gems_solcal%slitwav(1:nslit, currpix)
           !slitfit = omi_solslitfit(:, :, :, currpix) 
        ENDIF
     ELSE 
        solwinfit = gems_solcal%winfit(:, :, :, currpix)
     END IF
  ENDIF
  
  RETURN

END SUBROUTINE gems_o3p_adj_solar_data 

SUBROUTINE gems_o3p_adj_earthshine_data (theline, nlines_max, pge_error_status)

  USE OMSAO_precision_module
  USE OMSAO_indices_module,         ONLY: wvl_idx, spc_idx, sig_idx, n_max_fitpars, &
       solar_idx, fsl_idx, rsl_idx, comm_idx, com1_idx
  USE OMSAO_parameters_module,      ONLY: normweight, mrefl, max_fit_pts, deg2rad, maxchlen
  USE OMSAO_variables_module,       ONLY: curr_rad_spec, curr_sol_spec, n_rad_wvl, &
       use_meas_sig, numwin, nradpix, the_sza_atm, the_vza_atm, the_aza_atm, the_sca_atm,  &
       the_month, the_year, the_day, the_lon, the_lat, the_lats, the_lons, edgelons, edgelats, &
       the_surfalt, nview, nloc, the_utc, n_radwvl_sav, radwvl_sav,  nradpix_sav, saa_minlat, &
       saa_maxlat, saa_minlon, saa_maxlon, saa_minlat1, saa_maxlat1, saa_minlon1, saa_maxlon1, &
       do_bandavg, refidx, fitvar_rad_saved, n_fitvar_rad, radwavcal_freq, currpix, currloop, &
        n_irrad_wvl, nsolpix, actspec_rad, database, band_selectors, &
       refidx_sav, mask_fitvar_rad, lo_radbnd, up_radbnd, radnhtrunc, refnhextra, refdbdir, the_time

  USE ozprof_data_module,           ONLY: div_rad, div_sun, rad_posr, rad_specr, nsaa_spike, saa_flag, &
       the_cfrac, the_ctp, the_cld_flg, which_cld, the_cod, the_orig_cfr, the_orig_ctp, scacld_initcod,&
       the_orig_cod, the_ai, radcalwrt, biasfname, biascorr, l1l2inp_unit, which_biascorr, nrefl,      &
       aerosol, which_aerosol, scale_aod, scaled_aod, do_simu, the_fixalb, do_lambcld, lambcld_refl, &
       has_glint, glintprob, ozprof_start_index, ozprof_end_index, sun_posr, sun_specr, pos_alb, nir, &
       the_snowice, the_landwater_flag, the_glint_flag
  USE GEMS_O3P_gemsdata_module, ONLY:mchannel, nchannel,nxtrack, nfxtrack,nxtrack_max, &
                                     gems_clouds,gems_rad, gems_irrad, gems_refl, gems_gflag,&
				     gems_time, gems_sza, gems_vza,gems_aza, gems_sca,  &  ! omi_esca, omi_eaza
				     gems_lon, gems_lat,gems_clon, gems_clat, gems_elon, gems_elat, gems_height, orbnumsol, & 
                                     gems_initval, gems_exitval, gems_fitvar
  USE OMSAO_errstat_module 

  IMPLICIT NONE
  
  ! =================
  ! In/Out variables
  ! =================
  INTEGER, INTENT (IN)        :: theline, nlines_max
  INTEGER, INTENT (OUT)       :: pge_error_status

  ! =================
  ! Local variables
  ! ================= 
  REAL (KIND=r8), PARAMETER:: pi         = 3.14159265358979_r8  ! 2*ASIN(1.0_r8)
  REAL (KIND=r8), PARAMETER:: rad2deg    = 180.0_r8 / pi

  INTEGER                     :: hour, minute, fidx, lidx, i, j, west_idx, south_idx, idxoff, &
       nhtrunc, ntrunc, ntrunc1, errstat, ntempx, nch, ix, nord, ch, nw, is, nsub, idum, iw
  INTEGER (KIND=i4)           :: estat
  REAL (KIND=r8)              :: second, finit
  REAL (KIND=dp), DIMENSION (n_max_fitpars) :: fitvar
  LOGICAL                     :: redo_database
  
  ! xliu (02/03/2007): variables for correcting across-track dependent biases
  INTEGER, PARAMETER                                             :: maxord = 12
  INTEGER, DIMENSION (mchannel), SAVE                              :: corr_npars, nxcorr, nxwav
  LOGICAL, SAVE                                                  :: first = .TRUE.
  REAL (KIND=dp), DIMENSION(mchannel, nxtrack_max, 0:maxord), SAVE :: corrpars, offset_pars, slope_pars
  REAL (KIND=dp), DIMENSION(mchannel), SAVE                        :: corr_woffset
  REAL (KIND=dp), DIMENSION(max_fit_pts)                         :: corr, offset, slope, slope1, xoffset
  REAL (KIND=dp), DIMENSION(nxtrack_max, max_fit_pts), SAVE      :: allcorr, alloffset
  REAL (KIND=dp), DIMENSION(mchannel, nxtrack_max), SAVE           :: xcorr
  REAL (KIND=dp), DIMENSION(mchannel,nxtrack_max,max_fit_pts),SAVE :: xwcorr, xwslp, xwoff, &
       cldclrdf_xwcorr, xwslp1, xwxoff
  REAL (KIND=dp), DIMENSION(nxtrack_max), SAVE                   :: avg347clr
  REAL (KIND=dp), DIMENSION(mchannel, max_fit_pts), SAVE           :: xwavs
  INTEGER, SAVE                                                  :: nxgascorr, nxw2corr
  REAL (KIND=dp), DIMENSION(nxtrack_max, max_fit_pts, 2), SAVE   :: gascorr, xw2corr
  INTEGER, DIMENSION(nxtrack_max, 3), SAVE                       :: gascorr_npts, xw2corr_npts
  REAL (KIND=dp), DIMENSION(1:maxord, max_fit_pts)               :: del
  REAL (KIND=dp), DIMENSION (2, max_fit_pts)                     :: strayspec
  REAL (KIND=dp)                                                 :: woffset, rad347, irad347
  CHARACTER (LEN=maxchlen)                                       :: cldclrdf_biasfname, gascorr_fname, xw2corr_fname

    
  ! ================================
  !   External functions
  ! ================================
  INTEGER (KIND=i4), EXTERNAL :: PGS_TD_TAItoUTC  
  INTEGER :: OMI_SMF_setmsg

  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=26), PARAMETER :: modulename = 'gems_adjust_earthshine_data'

  pge_error_status = pge_errstat_ok

  IF (first .AND. biascorr) THEN
    !T  F   bias correction, time-dependent degradation correction
    !7      which biascorr (1: direct corr. 2: para. corr 3: dirrect corr (offset & slope) 4: para. (offset & slope) 
    !       5. same as 3 but LOG(I) 6. xtrack correction (theta) 7 xtrack and wave (don't need to match positions like 1) 
    !       8. xw, relative difference I/F diff. vs. rad347  9.  xw, relative I/F difference vs. rad347 nm, but with option 7 correction as offset
    !       It is much better to use direct correction instead of parameterized correction
   IF ( which_biascorr == 7) THEN ! need xw_bias*.dat
        OPEN(UNIT=l1l2inp_unit, FILE=TRIM(ADJUSTL(biasfname)), STATUS='OLD', IOSTAT=errstat)
        IF ( errstat /= pge_errstat_ok ) THEN
           errstat = OMI_SMF_setmsg (omsao_e_open_fitctrl_file, &
           TRIM(ADJUSTL(biasfname)), modulename, 0)
           pge_error_status = pge_errstat_error; RETURN
        ENDIF
        xwcorr = 1.0 ! Initialize to one
        DO is = 1, mchannel
           READ (l1l2inp_unit, *) nxcorr(is), nxwav(is)
           DO iw = 1, nxwav(is)
              READ (l1l2inp_unit, *) xwavs(is, iw), xwcorr(is, 1:nxcorr(is), iw)
           ENDDO
        
           IF (nxcorr(is) > nfxtrack) THEN
              nsub = nxcorr(is) / nfxtrack
              DO ix = 1, nfxtrack
                 fidx = (ix - 1) * nsub + 1
                 lidx = fidx + nsub - 1
                 DO iw = 1, nxwav(is)
                    xwcorr(is, ix, iw) = SUM(xwcorr(is, fidx:lidx, iw)) / nsub
                 ENDDO
              ENDDO
           ENDIF
        ENDDO
        CLOSE(UNIT=l1l2inp_unit) 

        ! Correction for trace gases
        gascorr_fname = ADJUSTL(TRIM(refdbdir)) // 'OMI/OMIO3PROF_corrgas_hres_o10582.dat'
     
        OPEN (UNIT=l1l2inp_unit, FILE=TRIM(ADJUSTL(gascorr_fname)), STATUS='UNKNOWN', IOSTAT=errstat)
        IF ( errstat /= pge_errstat_ok ) THEN
           WRITE(*, '(2A)') modulename, ': Cannot open trace gas correction file!!!'
           errstat = pge_errstat_error; RETURN
        END IF
        
        READ (l1l2inp_unit, *); READ (l1l2inp_unit, *); READ(l1l2inp_unit, *) 
        READ (l1l2inp_unit, *) nxgascorr
        gascorr(:, :, 2) = 1.d0
        DO ix = 1, nxgascorr
           READ (l1l2inp_unit, *) idum, gascorr_npts(ix, 1:3)
           DO i = 1, gascorr_npts(ix, 1)
              READ (l1l2inp_unit, *) gascorr(ix, i, 1:2)
              gascorr(ix, i, 2) = EXP(gascorr(ix, i, 2))
           ENDDO
        ENDDO
        CLOSE(UNIT=l1l2inp_unit) 

        ! Additional correction for x-track dependent biases
        xw2corr_fname = ADJUSTL(TRIM(refdbdir)) // 'OMI/OMIO3PROF_hres_xwcorr-2006m071116.dat'
     
        OPEN (UNIT=l1l2inp_unit, FILE=TRIM(ADJUSTL(xw2corr_fname)), STATUS='UNKNOWN', IOSTAT=errstat)
        IF ( errstat /= pge_errstat_ok ) THEN
           WRITE(*, '(2A)') modulename, ': Cannot open additional x-track dependent correction file!!!'
           errstat = pge_errstat_error; RETURN
        END IF
        
        READ (l1l2inp_unit, *)
        READ (l1l2inp_unit, *) nxw2corr
        xw2corr(:, :, 2) = 0.0d0
        DO ix = 1, nxw2corr
           READ (l1l2inp_unit, *) idum, xw2corr_npts(ix, 1:3)
           DO i = 1, xw2corr_npts(ix, 1)
              READ (l1l2inp_unit, *) xw2corr(ix, i, 1:2)
           ENDDO
        ENDDO
        CLOSE(UNIT=l1l2inp_unit) 
   ELSE
     WRITE(*,'(A)') 'Which_bias should be 7 in'//modulename
   ENDIF
     first = .FALSE.
  ENDIF

  ! Radiance Spectrum
  n_rad_wvl = gems_rad%nwav(currpix, currloop) 
  div_rad   = gems_rad%norm(currpix, currloop)

 
  curr_rad_spec(wvl_idx, 1:n_rad_wvl) = gems_rad%wavl(1:n_rad_wvl, currpix, currloop) 
  curr_rad_spec(spc_idx, 1:n_rad_wvl) = gems_rad%spec(1:n_rad_wvl, currpix, currloop)   
  IF (use_meas_sig) THEN
     curr_rad_spec(sig_idx, 1:n_rad_wvl) = gems_rad%prec(1:n_rad_wvl, currpix, currloop)
  ELSE
     curr_rad_spec(sig_idx, 1:n_rad_wvl) = normweight
  ENDIF
  nradpix(1:numwin) = gems_rad%npix(1:numwin, currpix, currloop)     ! Solar Spectrum


  IF ( biascorr .AND. which_biascorr == 7 ) THEN
     fidx = 1
     DO i = 1, numwin
        lidx = fidx + nradpix(i) - 1; ch = band_selectors(i) 
        !print *, ch, currpix, fidx, lidx, nradpix(i), nxwav(ch)
        CALL INTERPOL(xwavs(ch, 1:nxwav(ch)), xwcorr(ch, currpix, 1:nxwav(ch)), nxwav(ch), &
             curr_rad_spec(wvl_idx, fidx:lidx),  corr(1:nradpix(i)), nradpix(i), errstat)
        IF (errstat < 0) THEN
           WRITE(*, *) modulename, ': INTERPOL error, errstat = ', errstat
           errstat = pge_errstat_error; RETURN
        ENDIF     
        curr_rad_spec(spc_idx, fidx:lidx) = curr_rad_spec(spc_idx, fidx:lidx) / corr(1:nradpix(i))
        fidx = lidx + 1
     ENDDO
  ENDIF

  ! IRRadiance Spectrum
  n_irrad_wvl = n_rad_wvl
  div_sun     = gems_irrad%norm(currpix)
  curr_sol_spec(wvl_idx, 1:n_irrad_wvl) = gems_irrad%wavl(gems_rad%wind(1:n_rad_wvl, currpix, currloop), currpix) 
  curr_sol_spec(spc_idx, 1:n_irrad_wvl) = gems_irrad%spec(gems_rad%wind(1:n_rad_wvl, currpix, currloop), currpix) 

  IF (use_meas_sig .AND. orbnumsol /= 99999) THEN
     curr_sol_spec(sig_idx, 1:n_irrad_wvl) = gems_irrad%prec(gems_rad%wind(1:n_rad_wvl, currpix, currloop), currpix)
  ELSE IF (orbnumsol == 99999) THEN
     curr_sol_spec(sig_idx, 1:n_irrad_wvl) = 0.0  ! Ignore error in solar irradiance
  ELSE
     curr_sol_spec(sig_idx, 1:n_irrad_wvl) = normweight
  ENDIF
  nsolpix(1:numwin) = nradpix(1:numwin)  
  

  ! Reflectance spectrum at ~370 nm ==> 345-348
  rad_posr (1:nrefl) = gems_refl%radwavl(1:nrefl, currpix, currloop)
  rad_specr(1:nrefl) = gems_refl%radspec(1:nrefl, currpix, currloop)

  IF ( biascorr ) THEN
      IF ( which_biascorr == 7 ) THEN
        IF ( xwavs(mchannel, nxwav(mchannel)) < rad_posr(1) ) THEN
           rad_specr(1:nrefl) = rad_specr(1:nrefl) / xwcorr(mchannel, currpix, nxwav(mchannel))
        ELSE
           fidx = MINVAL( MINLOC( xwavs(mchannel, 1:nxwav(mchannel)), MASK = &
                (xwavs(mchannel, 1:nxwav(mchannel)) > rad_posr(1) )))
           lidx = MINVAL(MAXLOC( xwavs(mchannel, 1:nxwav(mchannel)), MASK = &
                (xwavs(mchannel, 1:nxwav(mchannel)) < rad_posr(nrefl) )))
           IF (fidx > lidx) THEN
              idum = fidx; fidx = lidx; lidx = idum
           ENDIF
           rad_specr(1:nrefl) = rad_specr(1:nrefl) / &
                (SUM(xwcorr(mchannel, currpix, fidx:lidx)) / (lidx - fidx + 1))
         ENDIF
      ENDIF
  ENDIF

  nview       = 1
  the_sza_atm = gems_sza  (currpix, currloop)
  the_vza_atm = gems_vza  (currpix, currloop)
  the_aza_atm = gems_aza  (currpix, currloop)

  ! compute scattering angle for polarization correction
  the_sca_atm = gems_sca (currpix, currloop)
  the_lon     = gems_lon (currpix, currloop)
  the_lat     = gems_lat (currpix, currloop)
  the_surfalt = gems_height   (currpix, currloop) / 1000.
 
  nloc          = 5
  the_lons(1) = gems_clon  (currpix-1, currloop)
  the_lons(2) = gems_clon  (currpix-1, currloop + 1)
  the_lons(3) = gems_clon  (currpix,   currloop + 1)
  the_lons(4) = gems_clon  (currpix,   currloop)
  the_lons(5) = the_lon
  the_lats(1) = gems_clat  (currpix-1, currloop)
  the_lats(2) = gems_clat  (currpix-1, currloop + 1)
  the_lats(3) = gems_clat  (currpix,   currloop + 1)
  the_lats(4) = gems_clat  (currpix,   currloop)
  the_lats(5) = the_lat
  edgelons(1) = gems_elon  (currpix-1, currloop)
  edgelons(2) = gems_elon  (currpix,   currloop)
  edgelats(1) = gems_elat  (currpix-1, currloop)
  edgelats(2) = gems_elat  (currpix,   currloop)

  the_time = gems_time(currloop)
  estat = PGS_TD_TAItoUTC(gems_time(currloop), the_utc)
  READ (the_utc, '(I4, 1x, I2, 1x, I2, 1x, I2, 1x, I2, 1x, F9.6)') the_year, the_month, &
       the_day, hour, minute, second

  ! Chechk cloud fractions !now which_cld = 3 (OMCLDO2 + TOMS) ELSE 2: only TOMS
  IF (which_cld /= 2) THEN
     the_cld_flg = gems_clouds%qflags(currpix, currloop)
     IF (the_cld_flg /= 10) THEN  ! Bad clouds for 10
        the_cfrac = gems_clouds%cfr(currpix, currloop)
        the_ctp   = gems_clouds%ctp(currpix, currloop)
     ELSE
        the_cfrac = 0.0;      the_ctp = 0.0
     ENDIF
  ENDIF

  IF (which_cld == 2 .OR. (the_cld_flg == 10 .AND. which_cld >= 3 ))  THEN
     the_cld_flg = 2  ! Derived based on OMCLDRR
     CALL GET_TOMSV8_CTP(the_month, the_day, the_lon, the_lat, the_ctp, pge_error_status)
     the_cfrac = 0.5  ! will be updated anyway at longer wavelength
     the_ai = -999.0
  ENDIF

  ! Special treatments for sea glint
  has_glint = .FALSE.; glintprob = 0.0
  ! Land-water flag=1: >=8 not used, else contain water
  IF  (gems_gflag%land_water(currpix, currloop) /= 1 .AND. gems_gflag%land_water(currpix, currloop) < 8) THEN 
     IF (gems_gflag%glint(currpix, currloop) == 1) THEN
        has_glint = .TRUE.
        CALL SUNGLINT_PROBABILITY (the_sza_atm, the_vza_atm, the_aza_atm, glintprob)
        !PRINT *, 'Glint Probability: ', glintprob, the_cfrac
        !IF (the_cfrac < 0.30 * glintprob) the_cfrac = 0.0
     ENDIF
  ENDIF

  ! Snow/ice flag
  ! 00: Snow-free land
  ! 1-100: sea ice concentration
  ! 101: permanent ice
  ! 102: not used
  ! 103: snow
  ! 104: ocean
  !105-123: Reserved
  !124: mixed pixels
  !125: suspect ice value
  !126: corners
  !17: Error

  the_snowice   = gems_gflag%snow_ice(currpix, currloop)
  the_landwater_flag = gems_gflag%land_water(currpix, currloop)
  the_glint_flag     = gems_gflag%glint(currpix, currloop)

  IF ( do_lambcld ) THEN
     the_cod = 0.0
  ELSE
     ! Pixel-independent approximation: cloudy scence with an effective COD 20.0 
     ! (cloud thickness 100 mb) and clear-sky scene. If cloud fraction is 20, 
     ! then rederive the effective COD.  Since CTP from OMI products are based on
     ! Lambertian cloud model, it is better to assume thin cloud layer (e.g., 100 mb)
     ! even for thick clouds
     the_cod = scacld_initcod
  ENDIF

  IF (do_simu .AND. .NOT. radcalwrt) THEN
     OPEN(UNIT=l1l2inp_unit, FILE='INP/sim.inp', STATUS='unknown')
     READ(l1l2inp_unit, *) the_sza_atm, the_vza_atm, the_aza_atm, the_fixalb, the_surfalt, &
          the_cfrac, the_ctp, the_cod, the_lon, the_lat, the_month, the_day, which_aerosol, &
          scaled_aod, do_lambcld, lambcld_refl
     IF (which_aerosol < 0 ) THEN
        aerosol = .FALSE.
        scale_aod = .FALSE.
        scaled_aod = 0.0
     ELSE
        aerosol = .TRUE.
        scale_aod = .TRUE.
     ENDIF

     IF (the_cfrac == 0.0 .OR. the_cod == 0) THEN
        the_ctp = 0.0; the_cod = 0.0; the_cfrac = 0
     ENDIF
     
     IF (do_lambcld ) THEN
        the_cod = 0.0
     ENDIF

     CLOSE (l1l2inp_unit)
  ENDIF

  ! These properties may be slightly modified later (save them)
  the_orig_cfr = the_cfrac
  the_orig_ctp = the_ctp
  the_orig_cod = the_cod

  ! Detect Spikes over the South Atlantic Anomaly region 
  IF (.NOT. do_simu .AND. .NOT. radcalwrt) THEN
     IF ( ((the_lat > saa_minlat  .AND. the_lat < saa_maxlat    .AND. &
          the_lon  > saa_minlon  .AND. the_lon < saa_maxlon)   .OR.  & 
          (the_lat > saa_minlat1 .AND. the_lat < saa_maxlat1   .AND. &
          the_lon  > saa_minlon1 .AND. the_lon < saa_maxlon1)) .AND. .NOT. do_simu )  THEN   
        CALL ROUGH_SPIKE_DETECT(n_rad_wvl, curr_rad_spec(wvl_idx, 1:n_rad_wvl), &
             curr_rad_spec(spc_idx, 1:n_rad_wvl), curr_sol_spec(spc_idx, 1:n_rad_wvl), nsaa_spike)
        saa_flag = .FALSE. !; nsaa_spike = 0
     ELSE
        saa_flag = .FALSE.;     nsaa_spike = 0
     ENDIF
  ELSE
     saa_flag = .FALSE.;     nsaa_spike = 0
  ENDIF
 
  ! Obtain measurement error in term sun-normalized radiance
  CALL gems_o3p_adj_sig (curr_rad_spec(wvl_idx:sig_idx, 1:n_rad_wvl),  &
       curr_sol_spec(wvl_idx:sig_idx, 1:n_rad_wvl), n_rad_wvl)
    
  ! Make sure that reference spectra has  more wavelengths than
  ! irradiance interpolation and shifting      
  nhtrunc = radnhtrunc
  ntrunc = nhtrunc * 2; ntrunc1 = ntrunc + 1

  fidx = 1
  DO i = 1, numwin
     lidx = fidx + nradpix(i) - ntrunc1
     curr_rad_spec(1:sig_idx, fidx:lidx) = curr_rad_spec(1:sig_idx, fidx + nhtrunc : lidx + nhtrunc)
     IF (lidx  < n_rad_wvl - ntrunc1 ) THEN
        curr_rad_spec(1:sig_idx, lidx+1:n_rad_wvl - ntrunc) =  &
             curr_rad_spec(1:sig_idx, lidx+ntrunc1:n_rad_wvl)
     ENDIF    
     nradpix(i) = nradpix(i) - ntrunc; fidx = lidx + 1; n_rad_wvl = n_rad_wvl - ntrunc
  ENDDO
 
  ! save the original grids for later obtain ozone cross section
  n_radwvl_sav = n_rad_wvl; nradpix_sav = nradpix
  radwvl_sav(1:n_rad_wvl) = curr_rad_spec(wvl_idx, 1:n_rad_wvl)  

  ! redo_database = .FALSE.
  redo_database = .TRUE.  !wasp

  IF (theline >= 1) THEN

     IF (gems_rad%nwav(currpix, currloop) /= gems_rad%nwav(currpix, currloop-1)) THEN
        redo_database = .TRUE.
     ELSE 
        IF (ANY(gems_rad%wind(1:n_rad_wvl, currpix, currloop) - &
             gems_rad%wind(1:n_rad_wvl, currpix, currloop-1) /= 0)) redo_database = .TRUE.
     ENDIF
  ENDIF

  IF ( MOD (theline, radwavcal_freq) == 0 .OR. redo_database) THEN 
     ! --------------------------------------------------------------
     ! Spline data bases, compute undersampling spectrum, and prepare
     ! reference spectra for fitting.
     ! ---- ----------------------------------------------------------
    
     CALL prepare_databases ( n_rad_wvl, curr_rad_spec(wvl_idx,1:n_rad_wvl), pge_error_status )
     IF ( pge_error_status >= pge_errstat_error ) RETURN

     !CALL avg_band_spec(curr_rad_spec(wvl_idx, 1:n_rad_wvl), strayspec(1, 1:n_rad_wvl), &
     !     n_rad_wvl, idum, errstat)
     !CALL avg_band_spec(curr_rad_spec(wvl_idx, 1:n_rad_wvl), strayspec(2, 1:n_rad_wvl), &
     !     n_rad_wvl, idum, errstat)
  ENDIF 

  !WRITE original wavelengths and solar irradiance spectra
  !DO i = 1, n_radwvl_sav
  !   WRITE(96, *) radwvl_sav(i), i0sav(refidx_sav(i)) * div_sun
  !ENDDO

  ! average and subsampling for selected bands and update nradpix
  IF (do_bandavg) THEN 
     CALL avg_band_radspec (curr_rad_spec(wvl_idx:sig_idx, 1:n_rad_wvl), &
          n_rad_wvl, pge_error_status)
     IF ( pge_error_status >= pge_errstat_error ) RETURN
  ENDIF

  fidx = 1
  DO i = 1, numwin
     lidx = fidx + nradpix(i) - 1
     idxoff = refnhextra + (i - 1) * 2 * refnhextra
     refidx(fidx:lidx) = (/(j, j = fidx + idxoff, lidx + idxoff)/)
     fidx = lidx  + 1
  ENDDO

  IF (radcalwrt) THEN
     actspec_rad(1:n_rad_wvl) = curr_rad_spec(spc_idx, 1:n_rad_wvl) / &
          database(solar_idx, refidx(1:n_rad_wvl)) * div_rad / div_sun
  ENDIF

  ! load databases for common modes
        
  IF ( MOD (theline, radwavcal_freq) == 0 .OR. redo_database) THEN
     CALL load_omi_comres(pge_error_status)
      
     IF ( pge_error_status >= pge_errstat_error ) RETURN
  ENDIF
   

  
  IF (biascorr) THEN
     IF ( which_biascorr == 7) THEN
	    IF (curr_rad_spec(1,1) > 271) then 
        curr_rad_spec(spc_idx, 1:n_rad_wvl)  = curr_rad_spec(spc_idx, 1:n_rad_wvl) * gascorr(currpix, 1:n_rad_wvl, 2)
        curr_rad_spec(spc_idx, 1:n_rad_wvl)  = curr_rad_spec(spc_idx, 1:n_rad_wvl) * &
                                               (1.0d0 + xw2corr(currpix, 1:n_rad_wvl, 2) / 100.)
         !print * , 'which_biascorr=7 is blocked for different fitting window' ! jbak 
        ENDIF
     ENDIF
  ENDIF
  
  ! Initialized fitted variables from valid western and southern neighbors 
  ! since the retrievals are performed from west to east and from south to north
  west_idx  = currpix  - 1; south_idx = currloop - 1
  IF (south_idx < 0) south_idx = nlines_max - 1

  gems_initval(currpix, currloop) = 0.0; fitvar = 0.0; finit = 0.0
  IF (west_idx > 0) THEN
     IF (gems_exitval(west_idx, currloop) > 0) THEN  ! Western pixel (success retrieval)
        fitvar(1:n_fitvar_rad) = fitvar(1:n_fitvar_rad) + &
             gems_fitvar(west_idx, currloop, 1:n_fitvar_rad)
        finit = finit + 1.0
     ENDIF
     
     IF (south_idx >= 0 .AND. south_idx /= nlines_max - 1) THEN
        IF (gems_exitval(west_idx, south_idx) > 0) THEN ! Southwestern pixel (success retrieval)
           fitvar(1:n_fitvar_rad) = fitvar(1:n_fitvar_rad) &
                + gems_fitvar(west_idx, south_idx, 1:n_fitvar_rad) * 0.5
           finit = finit + 0.5
        ENDIF
     ENDIF
  ENDIF
  
  IF ( south_idx >= 0 ) THEN
     IF (gems_exitval(currpix, south_idx) > 0) THEN     ! Southern pixel (success retrieval)
        fitvar(1:n_fitvar_rad) = fitvar(1:n_fitvar_rad) + &
             gems_fitvar(currpix, south_idx, 1:n_fitvar_rad) 
        finit = finit + 1.0
     ENDIF
  ENDIF

  IF (finit > 0) THEN
     fitvar_rad_saved(mask_fitvar_rad(1:n_fitvar_rad)) = fitvar(1:n_fitvar_rad) / finit
     gems_initval(currpix, currloop) = 1
  ENDIF

  RETURN
END SUBROUTINE gems_o3p_adj_earthshine_data


SUBROUTINE gems_o3p_adj_sig (radspec, solspec, np)

  USE OMSAO_precision_module
  USE OMSAO_indices_module,    ONLY: wvl_idx, spc_idx, sig_idx
  USE OMSAO_parameters_module, ONLY: maxwin
  USE OMSAO_variables_module,  ONLY: numwin, nradpix, band_selectors
  USE ozprof_data_module,      ONLY: div_rad, div_sun, use_lograd, use_flns
  
  IMPLICIT NONE
  
  ! ====================
  ! In/Output variables
  ! ====================
  INTEGER, INTENT(IN)                             :: np
  REAL (KIND=dp), DIMENSION(3, np), INTENT(IN)    :: solspec
  REAL (KIND=dp), DIMENSION(3, np), INTENT(INOUT) :: radspec
  
  ! ====================
  ! Local variables
  ! ====================
  INTEGER, PARAMETER                :: nreg = 3
  INTEGER                           :: i, j, iwin, fidx, lidx
  REAL(KIND=dp)                     :: dw1, dw2
  REAL (KIND=dp), DIMENSION(np)     :: relsig, normrad, sig
  REAL (KIND=dp), DIMENSION(maxwin) :: floor_noise =  &
       (/0.004, 0.002, 0.001, 0.001, 0.001/)
  REAL (KIND=dp), DIMENSION(nreg) :: reg_noise =  &
       (/0.004, 0.004, 0.002/)
  REAL (KIND=dp), DIMENSION(0:nreg) :: reg_waves = &
       (/260.0, 300.0, 310.0, 350./)
  
  ! I/F
  normrad = radspec(spc_idx,:) / solspec(spc_idx,:) * (div_rad / div_sun) 

  ! Relative photon noise (I)
  relsig = radspec(sig_idx,:) / radspec(spc_idx,:)
  relsig = SQRT( relsig ** 2.0 + (solspec(sig_idx, :) / solspec(spc_idx,:)) ** 2.0)

  IF (use_flns) THEN
     fidx = 1;     i = fidx
     DO j = 1, nreg
        DO WHILE (i <= np )
           IF (radspec(wvl_idx, i) < reg_waves(j)) THEN
              i = i + 1
           ELSE
              exit
           ENDIF
        ENDDO
        lidx = i - 1
        WHERE (relsig(fidx:lidx) < reg_noise(j))              
           relsig(fidx:lidx) = reg_noise(j)
        ENDWHERE
!          print * , relsig(fidx:lidx)
!          print *,'=='
        fidx = lidx + 1
     ENDDO
  ENDIF
  IF (use_lograd) THEN        ! error in logarithmic radiance
     sig = LOG(1.0 + relsig)  ! ~relative error
  ELSE
     sig = relsig * normrad   ! absolute measurement error in I/F
  ENDIF
  
  radspec(sig_idx, :) = sig

  !DO i = 1, np - 1
  !   WRITE(92, *) radspec(1, i), sig(i)
  !ENDDO
  !STOP
  
  RETURN
END SUBROUTINE gems_o3p_adj_sig



! Correction for wavelength registration at 1:67 and 498:557
SUBROUTINE corruv2wav(nw, nx, waves)
  USE OMSAO_precision_module
  
  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                  INTENT (IN)    :: nw, nx
  REAL (KIND=r4), DIMENSION (nw, nx), INTENT (INOUT) :: waves

  ! ----------------
  ! Local variables
  ! ---------------
  INTEGER        :: ix
  REAL (KIND=r4) :: del, ndel, delp1, delp2, delm1, delm2, ndelp1, ndelm1, sh1, sh2

  !RETURN

  !print *, nw, nx
  DO ix = 1, nx   
     ! At position 67
     del    = waves(68, ix) - waves(67, ix)
     delm1  = waves(67, ix) - waves(66, ix)
     delp1  = waves(69, ix) - waves(68, ix)
     delp2  = waves(70, ix) - waves(69, ix)
     ndel   = delp1 * 2 - delp2
     ndelm1 = ndel * 2  - delp1
     
     ! Shifts
     sh1 = (ndel - del)
     sh2 = sh1 + (ndelm1 - delm1)
     
     !print *, 67, sh1, sh2
     waves(67, ix)   = waves(67, ix) - sh1
     waves(1:66, ix) = waves(1:66, ix) - sh2
    
     ! At position 498
     delm2 = waves(496, ix) - waves(495, ix)
     delm1 = waves(497, ix) - waves(496, ix)
     del   = waves(498, ix) - waves(497, ix)
     delp1 = waves(499, ix) - waves(498, ix)
     ndel = delm1 * 2 - delm2
     ndelp1 = ndel * 2 - delm1

     sh1 = (ndel - del)
     sh2 = sh1 + (ndelp1 - delp1)
     !print *, 498, sh1, sh2
     
     waves(498, ix) = waves(498, ix) + sh1
     waves(499:nw, ix) = waves(499:nw, ix) + sh2
 ENDDO

RETURN
END SUBROUTINE corruv2wav



SUBROUTINE convert_gpqualflag_info ( &
     nxtrack, omi_geoflg, land_water_flg, glint_flg, snow_ice_flg )

  USE OMSAO_precision_module
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                      INTENT (IN) :: nxtrack
  INTEGER (KIND=i2), DIMENSION (nxtrack), INTENT (IN) :: omi_geoflg

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER (KIND=i2), DIMENSION (nxtrack), INTENT (OUT) :: land_water_flg, glint_flg, snow_ice_flg

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i4),                PARAMETER      :: nbyte = 16
  INTEGER (KIND=i2), DIMENSION (7), PARAMETER      :: seven_byte = (/ 1, 2, 4, 8, 16, 32, 64 /)
  INTEGER (KIND=i2)                                :: i
  INTEGER (KIND=i2), DIMENSION (nxtrack)           :: tmp_flg
  INTEGER (KIND=i2), DIMENSION (nxtrack,0:nbyte-1) :: tmp_bytes

  ! ----------------------------
  ! Initialize output quantities
  ! ----------------------------
  land_water_flg = 0 ; glint_flg = 0 ; snow_ice_flg = 0

  ! -----------------------------------------------
  ! Save input variable in TMP_FLG for modification
  ! -----------------------------------------------
  tmp_flg(1:nxtrack) = omi_geoflg(1:nxtrack)  ;  tmp_bytes = 0

  CALL convert_2bytes_to_16bits ( &
       nbyte, nxtrack, tmp_flg(1:nxtrack), tmp_bytes(1:nxtrack,0:nbyte-1) )
  
  ! ------------------------------
  ! The Glint flag is easy: Byte 4
  ! ------------------------------
  glint_flg(1:nxtrack) = tmp_bytes(1:nxtrack,4)

  ! ------------------------------------------------------------------
  ! Land/Water and Ice require a bit more work. The BIT slices must be
  ! multiplied with the corresponding powers of 2. The sum over this
  ! product is the information we seek.
  ! ------------------------------------------------------------------
  DO i = 1, nxtrack
     land_water_flg(i) = SUM(tmp_bytes(i,0:3 )*seven_byte(1:4))
     snow_ice_flg  (i) = SUM(tmp_bytes(i,8:14)*seven_byte(1:7))
  END DO

  RETURN
END SUBROUTINE convert_gpqualflag_info

SUBROUTINE convert_2bytes_to_16bits ( nbits, ndim, byte_num, bit_num )

  ! ==========================================================
  ! Takes an NDIM dimensional 2Byte integer BYTE_NUM and
  ! converts it into an NDIM x 16 dimensional interger BIT_NUM
  ! ==========================================================

  USE OMSAO_precision_module
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                   INTENT (IN) :: ndim, nbits
  INTEGER (KIND=i2), DIMENSION (ndim), INTENT (IN) :: byte_num

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER (KIND=i2), DIMENSION (ndim,0:nbits-1), INTENT (OUT) :: bit_num

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i2)                   :: i, j
  REAL (KIND=dp)                      :: powval
  INTEGER (KIND=i2), DIMENSION (ndim) :: tmp_byte

  ! ----------------------------
  ! Initialize output quantities
  ! ----------------------------
  bit_num = 0_i2

  ! ------------------------------------------------
  ! Save input variable in TMP_BYTE for modification
  ! ------------------------------------------------
  tmp_byte(1:ndim) = byte_num(1:ndim)

  WHERE ( tmp_byte(1:ndim) < 0)
     tmp_byte(1:ndim) = tmp_byte(1:ndim) + 65536
  ENDWHERE

  ! -------------------------------------------------------------------
  ! Starting with the highest power NBYTES-1, subtract powers of 2 and
  ! assign 1 whereever the power fits in the flag number. At the end we
  ! arrive at a 16 BIT binary representation from which we can extract
  ! the surface information.
  ! -------------------------------------------------------------------
  DO i = nbits-1, 0, -1
     powval = 2.0 ** i
     IF ( powval > 0 ) THEN
        WHERE ( tmp_byte(1:ndim) >= powval )
           bit_num(1:ndim,i) = 1_i2
           tmp_byte(1:ndim) = tmp_byte(1:ndim) - powval
        ENDWHERE
     END IF
  END DO

  RETURN
END SUBROUTINE convert_2bytes_to_16bits

SUBROUTINE convert_16bits_to_2bytes ( nbits, ndim, bit_num, byte_num )

  ! ==========================================================
  ! Takes an NDIM dimensional 2Byte integer BYTE_NUM and
  ! converts it into an NDIM x 16 dimensional interger BIT_NUM
  ! ==========================================================

  USE OMSAO_precision_module
  IMPLICIT NONE

  ! ----------------------
  ! Input/Out variables
  ! ----------------------
  INTEGER (KIND=i4),                   INTENT (IN)           :: ndim, nbits
  INTEGER (KIND=i2), DIMENSION (ndim,0:nbits-1), INTENT (IN) :: bit_num
  INTEGER (KIND=i2), DIMENSION (ndim), INTENT (OUT)          :: byte_num

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i2) :: i

  ! ----------------------------
  ! Initialize output quantities
  ! ----------------------------
  byte_num = 0

  DO i = 0, nbits - 1 
     byte_num = byte_num + bit_num(:, i) * 2 ** i
  ENDDO

  WHERE (byte_num > 32767)
     byte_num = byte_num - 65536
  ENDWHERE 

  RETURN
END SUBROUTINE convert_16bits_to_2bytes


SUBROUTINE coadd_2bytes_qflgs(nbits, ndim, qflg1, qflg2) 
  
  USE OMSAO_precision_module
  IMPLICIT NONE

  ! --------------------------
  ! Input/Output variables
  ! --------------------------
  INTEGER (KIND=i4),                    INTENT (IN)  :: nbits, ndim
  INTEGER (KIND=i2), DIMENSION(ndim), INTENT (INOUT) :: qflg1
  INTEGER (KIND=i2), DIMENSION(ndim), INTENT (IN)    :: qflg2

  INTEGER (KIND=i2), DIMENSION (ndim, 0:nbits-1) :: bit_num1, bit_num2

  CALL convert_2bytes_to_16bits ( nbits, ndim, qflg1, bit_num1 )
  CALL convert_2bytes_to_16bits ( nbits, ndim, qflg2, bit_num2 )
  bit_num1 = bit_num1 + bit_num2

  WHERE(bit_num1 > 1) 
     bit_num1 = 1
  ENDWHERE
  
  ! Convert 16 bits to 2bytes
  CALL convert_16bits_to_2bytes (nbits, ndim, bit_num1, qflg1)

  RETURN
  
END SUBROUTINE coadd_2bytes_qflgs



subroutine timestamp (curtime )
!
!*******************************************************************************
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none
!
  character (len=24), INTENT(OUT) :: curtime
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 3 ), parameter, dimension(12) :: month = (/ &
    'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC' /)
  integer n
  integer s
  character ( len = 10 )  time
  integer values(8)
  integer y
  character ( len = 5 ) zone
!
  call date_and_time ( date, time, zone, values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  !if ( h < 12 ) then
  !  ampm = 'AM'
  !else if ( h == 12 ) then
  !  if ( n == 0 .and. s == 0 ) then
  !    ampm = 'Noon'
  !  else
  !    ampm = 'PM'
  !  end if
  !else
  !  h = h - 12
  !  if ( h < 12 ) then
  !    ampm = 'PM'
  !  else if ( h == 12 ) then
  !    if ( n == 0 .and. s == 0 ) then
  !      ampm = 'Midnight'
  !    else
  !      ampm = 'AM'
  !    end if
  !  end if
  !end if

  write ( curtime, '(a3,1x,i2.2,1x,i4,1x,i2.2,a1,i2.2,a1,i2.2,a1,i3.3)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm

  return
end subroutine timestamp


SUBROUTINE GET_DOY(theyear, themon,  theday, thedoy)
  USE OMSAO_precision_module
  IMPLICIT NONE

  ! --------------------------
  ! Input/Output variables
  ! --------------------------
  INTEGER, INTENT (IN)   :: theyear, themon, theday
  INTEGER, INTENT (OUT)  :: thedoy

  ! Local variables
  INTEGER, DIMENSION(12) :: ndays = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

  IF (MOD(theyear, 4) == 0) ndays(2) = 29
  IF (themon == 1) THEN
     thedoy = theday
  ELSE 
     thedoy = SUM(ndays(1:themon-1)) + theday
  ENDIF
  
  RETURN
END SUBROUTINE GET_DOY

! xliu, 03/25/2011, add several subroutines for bit-based operation of 8-bit unsigned integer
SUBROUTINE convert_byte_to_8bits ( nbits, ndim, byte_num, bit_num )

  ! ==========================================================
  ! Takes an NDIM dimensional unsigned Byte integer BYTE_NUM and
  ! converts it into an NDIM x 8 dimensional interger BIT_NUM
  ! ==========================================================

  USE OMSAO_precision_module
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                   INTENT (IN) :: ndim, nbits
  INTEGER (KIND=i1), DIMENSION (ndim), INTENT (IN) :: byte_num

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER (KIND=i1), DIMENSION (ndim,0:nbits-1), INTENT (OUT) :: bit_num

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i2)                   :: i, j
  REAL (KIND=dp)                      :: powval
  INTEGER (KIND=i2), DIMENSION (ndim) :: tmp_byte

  ! ----------------------------
  ! Initialize output quantities
  ! ----------------------------
  bit_num = 0_i1

  ! ------------------------------------------------
  ! Save input variable in TMP_BYTE for modification
  ! ------------------------------------------------
  ! Note: fortran does not have unsigned integer, so copy it to 2 bytes
  ! and 256 for negative nunmbers
  tmp_byte(1:ndim) = byte_num(1:ndim)
  WHERE ( tmp_byte(1:ndim) < 0 )
     tmp_byte(1:ndim) = tmp_byte(1:ndim) + 256
  ENDWHERE

  ! -------------------------------------------------------------------
  ! Starting with the highest power NBYTES-1, subtract powers of 2 and
  ! assign 1 whereever the power fits in the flag number. At the end we
  ! arrive at a 8 BIT binary representation from which we can extract
  ! the surface information.
  ! -------------------------------------------------------------------
  DO i = nbits-1, 0, -1
     powval = 2.0 ** i
     IF ( powval > 0 ) THEN
        WHERE ( tmp_byte(1:ndim) >= powval )
           bit_num(1:ndim,i) = 1_i1
           tmp_byte(1:ndim) = tmp_byte(1:ndim) - powval
        ENDWHERE
     END IF
  END DO

  RETURN
END SUBROUTINE convert_byte_to_8bits

SUBROUTINE convert_8bits_to_byte ( nbits, ndim, bit_num, byte_num )

  ! ==========================================================
  ! Takes an NDIM dimensional 2Byte integer BYTE_NUM and
  ! converts it into an NDIM x 16 dimensional interger BIT_NUM
  ! ==========================================================

  USE OMSAO_precision_module
  IMPLICIT NONE

  ! ----------------------
  ! Input/Out variables
  ! ----------------------
  INTEGER (KIND=i4),                   INTENT (IN)           :: ndim, nbits
  INTEGER (KIND=i1), DIMENSION (ndim,0:nbits-1), INTENT (IN) :: bit_num
  INTEGER (KIND=i1), DIMENSION (ndim), INTENT (OUT)          :: byte_num

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i2) :: i

  ! ----------------------------
  ! Initialize output quantities
  ! ----------------------------
  byte_num = 0

  DO i = 0, nbits - 1 
     byte_num = byte_num + bit_num(:, i) * 2 ** i
  ENDDO

  WHERE (byte_num > 127)
     byte_num = byte_num - 256
  ENDWHERE

  RETURN
END SUBROUTINE convert_8bits_to_byte



SUBROUTINE coadd_byte_qflgs(nbits, ndim, qflg1, qflg2) 
  
  USE OMSAO_precision_module
  IMPLICIT NONE

  ! --------------------------
  ! Input/Output variables
  ! --------------------------
  INTEGER (KIND=i4),                    INTENT (IN)  :: nbits, ndim
  INTEGER (KIND=i1), DIMENSION(ndim), INTENT (INOUT) :: qflg1
  INTEGER (KIND=i1), DIMENSION(ndim), INTENT (IN)    :: qflg2

  INTEGER (KIND=i1), DIMENSION (ndim, 0:nbits-1) :: bit_num1, bit_num2

  CALL convert_byte_to_8bits ( nbits, ndim, qflg1, bit_num1 )
  CALL convert_byte_to_8bits ( nbits, ndim, qflg1, bit_num2 )
  bit_num1 = bit_num1 + bit_num2

  WHERE(bit_num1 > 1) 
     bit_num1 = 1
  ENDWHERE
  
  ! Convert 8 bits to 1 byte
  CALL convert_8bits_to_byte (nbits, ndim, bit_num1, qflg1)

  RETURN
  
END SUBROUTINE coadd_byte_qflgs


SUBROUTINE convert_xtrackqfag_info ( nxtrack, omi_xtrackqflg, &
     rowanomaly_flg, waveshift_flg, blockage_flg, straysun_flg, strayearth_flg )

  USE OMSAO_precision_module
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                      INTENT (IN) :: nxtrack
  INTEGER (KIND=i1), DIMENSION (nxtrack), INTENT (IN) :: omi_xtrackqflg

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER (KIND=i1), DIMENSION (nxtrack), INTENT (OUT) :: rowanomaly_flg, waveshift_flg, &
       blockage_flg, straysun_flg, strayearth_flg

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i2),                PARAMETER      :: nbyte = 8
  INTEGER (KIND=i2), DIMENSION (7), PARAMETER      :: seven_byte = (/ 1, 2, 4, 8, 16, 32, 64 /)
  INTEGER (KIND=i2)                                :: i
  INTEGER (KIND=i1), DIMENSION (nxtrack)           :: tmp_flg
  INTEGER (KIND=i1), DIMENSION (nxtrack,0:nbyte-1) :: tmp_bytes

  ! ----------------------------
  ! Initialize output quantities
  ! ----------------------------
  rowanomaly_flg = 0; waveshift_flg = 0; blockage_flg = 0
  straysun_flg = 0; strayearth_flg = 0

  ! -----------------------------------------------
  ! Save input variable in TMP_FLG for modification
  ! -----------------------------------------------
  tmp_flg(1:nxtrack) = omi_xtrackqflg(1:nxtrack)  ;  tmp_bytes = 0

  CALL convert_byte_to_8bits (nbyte, nxtrack, tmp_flg(1:nxtrack), tmp_bytes(1:nxtrack,0:nbyte-1))

  waveshift_flg(1:nxtrack)  = tmp_bytes(1:nxtrack,4)
  blockage_flg(1:nxtrack)   = tmp_bytes(1:nxtrack,5)
  straysun_flg(1:nxtrack)   = tmp_bytes(1:nxtrack,6)
  strayearth_flg(1:nxtrack) = tmp_bytes(1:nxtrack,7)
  
  ! ------------------------------------------------------------------
  ! Row anomaly require a bit more work. The BIT slices must be
  ! multiplied with the corresponding powers of 2. The sum over this
  ! product is the information we seek.
  ! ------------------------------------------------------------------
  DO i = 1, nxtrack
     rowanomaly_flg(i) = SUM(tmp_bytes(i,0:2)*seven_byte(1:3))
  ENDDO

  RETURN
END SUBROUTINE convert_xtrackqfag_info



SUBROUTINE ROUGH_SPIKE_DETECT(ns, waves, rad, sol, nspike)

  USE OMSAO_precision_module
  IMPLICIT NONE

  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER, INTENT(IN)                             :: ns
  INTEGER, INTENT(OUT)                            :: nspike
  REAL (KIND=dp), INTENT(INOUT), DIMENSION (ns)   :: rad
  REAL (KIND=dp), INTENT(IN), DIMENSION (ns)      :: sol, waves

  ! =======================
  ! local variables
  ! =======================
  INTEGER                       :: i, j
  REAL (KIND=dp), PARAMETER     :: thresh = 0.20
  REAL (KIND=dp)                :: diff, oldratio, newratio
  REAL (KIND=dp), DIMENSION(ns) :: normrad

  !WRITE(90, '(f8.3, 2d14.6)') ((waves(i), rad(i), sol(i)), i=1, ns)

  normrad = rad / sol
  nspike  = 0           
  DO i = ns - 3, 1, - 1      ! only detect spike below 305 nm

     IF (waves(i) > 305.0D0) CYCLE

     oldratio = normrad(i) / normrad(i+1)
     diff = oldratio-1.0
     IF (diff > thresh ) THEN   !  This pixel got large error
        nspike = nspike + 1
        newratio = normrad(i+1) / normrad(i+2)
        
        IF (ABS(newratio - 1.0) < 0.5 * thresh) THEN
           normrad(i) = normrad(i+1) * newratio
           rad(i) = rad(i) * newratio / oldratio
        ELSE
           rad(i+1) =  rad(i+1) * normrad(i+2) / normrad(i+3) / newratio
           newratio =  normrad(i+2) / normrad(i+3)
           rad(i) = rad(i) * newratio / oldratio
        ENDIF
     ENDIF
  ENDDO
  !WRITE(*, *) 'Number of spikes = ', nspike
  
  RETURN
END SUBROUTINE ROUGH_SPIKE_DETECT


SUBROUTINE SPIKE_DETECT_CORRECT(ns, fitspec, simrad)
  
  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY : fitwavs, fitweights, currspec
  USE ozprof_data_module,     ONLY : use_lograd
  IMPLICIT NONE
  
  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER, INTENT(IN)                             :: ns
  REAL (KIND=dp), INTENT(INOUT), DIMENSION (ns)   :: simrad
  REAL (KIND=dp), INTENT(INOUT), DIMENSION (ns)   :: fitspec

  ! =======================
  ! local variables
  ! =======================
  INTEGER                        :: i, approach
  REAL (KIND=dp), DIMENSION(ns)  :: mratio, sratio
  REAL (KIND=dp)                 :: fitavg, simavg, fitavg3, &
       simavg3, fitavg2, simavg2

  IF (use_lograd) THEN
     simrad = EXP(simrad); fitspec = EXP(fitspec)
  ENDIF  
  !WRITE(90, '(f8.3, 2d14.6)') ((fitwavs(i), fitspec(i), simrad(i)), i=1, ns)

  ! use approach 1 is better
  approach = 1

  IF (approach == 1) THEN

     ! Assume the average of last 10 pixels are without errors
     ! Also assume the ratio of adjacent pixels is less dependent on ozone profile
     sratio(1:ns-1) = simrad(1:ns-1)  /  simrad(2:ns) ; sratio(ns) = sratio(ns-1)
     mratio(1:ns-1) = fitspec(1:ns-1) / fitspec(2:ns);  mratio(ns) = mratio(ns-1) 

     fitavg = SUM(fitspec(ns-9:ns)) / 10.0D0
     simavg = SUM(simrad(ns-9:ns))  / 10.0D0
     
     sratio(ns-10)   = simrad(ns-10)  / simavg
     mratio(ns-10)   = fitspec(ns-10) / fitavg
     fitspec(ns-10)  = fitspec(ns-10)  * sratio(ns-10) / mratio(ns-10)
     currspec(ns-10) = currspec(ns-10) * sratio(ns-10) / mratio(ns-10)  
     
     DO i = ns - 11, 1, - 1 
        mratio(i)  = fitspec(i)  / fitspec(i + 1)
        fitspec(i) = fitspec(i)  * sratio(i) / mratio(i)
        currspec(i)= currspec(i) * sratio(i) / mratio(i)  
     ENDDO
  ELSE

     ! Assume the average of last 20 pixels are without errors
     ! Also assume the ratio of adjacent pixels (first * third / second^2) &
     ! is less dependent on ozone profile

     sratio(1:ns-2) = simrad(1:ns-2)  * simrad(3:ns)  / (simrad(2:ns-1)**2.0)
     mratio(1:ns-2) = fitspec(1:ns-2) * fitspec(3:ns) / (fitspec(2:ns-1)**2.0)

     fitavg3= SUM(fitspec(ns-9:ns)) / 10.0D0
     simavg3= SUM(simrad(ns-9:ns))  / 10.0D0
     fitavg2= SUM(fitspec(ns-19:ns-10)) / 10.0D0
     simavg2= SUM(simrad(ns-19:ns-10))  / 10.0D0

     sratio(ns-20)  = simrad(ns-20)  * simavg3 / (simavg2**2.0)
     mratio(ns-20)  = fitspec(ns-20) * fitavg3 / (fitavg2**2.0)  
     fitspec (ns-20)= fitspec(ns-20) * sratio(ns-20) / mratio(ns-20)
     currspec(ns-20)= currspec(ns-20)* sratio(ns-20) / mratio(ns-20)  

     sratio(ns-21)  = simrad(ns-21)  * simavg2 / (simrad(ns-20)**2.0)
     mratio(ns-21)  = fitspec(ns-21) * fitavg2 / (fitspec(ns-20)**2.0) 
     fitspec (ns-21)= fitspec(ns-21) * sratio(ns-21) / mratio(ns-21)
     currspec(ns-21)= currspec(ns-21)* sratio(ns-21) / mratio(ns-21) 

     DO i = ns - 22, 1, - 1 
        mratio(i)  = fitspec(i)  * fitspec(i+2) / (fitspec(i + 1)**2.0)
        fitspec(i) = fitspec(i)  * sratio(i)    /  mratio(i)
        currspec(i)= currspec(i) * sratio(i)    /  mratio(i)  
     ENDDO
  ENDIF

  !WRITE(91, '(f8.3, 2d14.6)') ((fitwavs(i), fitspec(i), simrad(i)), i=1, ns) 
  IF (use_lograd) THEN
     simrad = LOG(simrad); fitspec = LOG(fitspec)
  ENDIF

  RETURN

END SUBROUTINE SPIKE_DETECT_CORRECT

! xliu, 05/23/2010
! Identify spikes (e.g. due to emission lines, protons) in UV-1
! And modify measurement error
SUBROUTINE UV1_SPIKE_DETECT(ns, fitspec, simrad, nspike)

  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY : fitwavs, fitweights, currspec, nradpix
  USE ozprof_data_module,     ONLY : use_lograd
  
  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER, INTENT(IN)                             :: ns
  REAL (KIND=dp), INTENT(INOUT), DIMENSION (ns)   :: simrad
  REAL (KIND=dp), INTENT(INOUT), DIMENSION (ns)   :: fitspec
  INTEGER, INTENT(OUT)                            :: nspike
  
  ! =======================
  ! local variables
  ! =======================
  INTEGER                        :: iw, j
  REAL (KIND=dp)                 :: rms, thresh, diff
  
  IF (use_lograd) THEN
     simrad = EXP(simrad); fitspec = EXP(fitspec)
  ENDIF
  
  nspike = 0
  j = nradpix(1)
  DO iw =  nradpix(1) - 1, 1, -1
     fratio = fitspec(iw) / fitspec(j)
     sratio = simrad(iw) / simrad(j)
     
     rms = SQRT(fitweights(iw)**2 + fitweights(j)**2) * 3.0 * SQRT(( 1.0 * j - iw))
     thresh = MAX(rms, 0.06)
     diff   = ABS(fratio - sratio)
     IF ( diff > thresh) THEN
        nspike = nspike + 1
        fitweights(iw) = diff / 3.0
        !WRITE(*, '(I5, 4F10.4)') nspike, fitwavs(iw), rms, thresh, diff
     ELSE
        j = iw
     ENDIF

  ENDDO       
  !WRITE(*, *) 'Number of spikes = ', nspike
  
  IF (use_lograd) THEN
     simrad = LOG(simrad); fitspec = LOG(fitspec)
  ENDIF

  !DO iw = 1, ns
  !   WRITE(90, *) fitwavs(iw), simrad(iw), fitspec(iw)
  !ENDDO
  !STOP
  
  RETURN
  
END SUBROUTINE UV1_SPIKE_DETECT
