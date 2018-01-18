! ******************************************************************************************
! When use_effcrs is set to true, outputs (radiance and weighting function) correspond to 
! input wavelengths. 
! When use_effcrs is set to false, i.e., use high resolution cross 
! section in the radiative transfer calculation and convolve with slot functions after
! radiative transfer calculation, waves do not correspond to the output rad and weighting 
! function, waves(1:nw) store the wavelengths (waves(1:ncalcp) where radiance calculation 
! are done (waves(ncalp+1:nw) = 0.0. And outputs corresponds to measurents (1:ns).
! Here nw is max(ncalcp, ns).
! ******************************************************************************************

SUBROUTINE LIDORT_PROF_ENV (do_ozwf, do_albwf, do_tmpwf, do_o3shi, ozvary, &
     do_taodwf, do_twaewf, do_saodwf, do_cfracwf, do_ctpwf, do_codwf, do_sprswf,  do_so2zwf, &
     nw, waves, nos, o3shi, sza, vza, aza, nl, ozprof, tprof, n0alb, albarr, &
     albpmin, albpmax, vary_sfcalb, walb0s, n0wfc, wfcarr, wfcpmin, wfcpmax, nostk, albwf, ozwf,  &
     tmpwf, o3shiwf, cfracwf, codwf, ctpwf, taodwf, twaewf, saodwf, sprswf, so2zwf, rad, errstat)

  USE OMSAO_precision_module
  USE OMSAO_parameters_module,ONLY  : maxlay, max_fit_pts, du2mol, rearth, max_spec_pts
  USE OMSAO_variables_module, ONLY  : rad_wav_avg, vza_atm, aza_atm,           &
       radwvl_sav, n_radwvl_sav, fitvar_rad, rmask_fitvar_rad, fitvar_rad_init, &
       numwin, nradpix, nradpix_sav, avgsza, avgvza, avgaza, b1ab_div_wav,     &
       database, database_shiwf, database_save, refidx, band_selectors,        &
       refwvl, n_refwvl, refspec_norm, the_surfalt, currloop, refspec_orig_data, &
       n_refspec_pts, winlim, n_rad_wvl
  USE OMSAO_indices_module,   ONLY : so2_idx, so2v_idx, o2o2_idx
  USE ozprof_data_module,     ONLY  : nt_fit, atmosprof, nup2p, nflay,         &
       start_layer, end_layer, do_multi_vza, do_lambcld, lambcld_refl, ntp,    &
       do_subfit, polcorr, the_cod, the_cfrac, maxawin, do_ch2reso, maxmom,    &
       nmom, fts, fps, fzs, fozs, frhos, gaext, gasca, gaasy, gamoms, gcq,     &
       gcw, gcasy, gcmoms, the_cbeta, aerwavs, has_clouds, aerosol, useasy,    &
       aerwavs, actawin, ncbp, nctp, num_iter, do_tracewf, mgasprof, ngas,     &
       tracegas, fgasidxs, gasidxs, nallgas, raycof, depol, do_radinter,       &
       gassidxs, fgassidxs, strat_aerosol, osfind, oswins, mflay, nfsfc, nsfc, &
       do_simu, maxgksec, ngksec, radcalwrt, wrtozcrs, &
       tropsca, tropaod, tropwaer, strataod, stratsca, taodind, taodfind,  &
       saodind, saodfind,  twaeind, twaefind, sprsind, nwfc, use_effcrs, &
       radcidxs, ncalcp, VlidortNstream
  USE OMSAO_errstat_module
  
  IMPLICIT NONE
  
  !===============================  Define Variables ===========================
  ! Include files of dimensions and numbers
  INCLUDE 'VLIDORT.PARS'
  
  ! Include files of input variables
  INCLUDE 'VLIDORT_INPUTS.VARS'
  INCLUDE 'VLIDORT_SETUPS.VARS'
  INCLUDE 'VLIDORT_L_INPUTS.VARS'
  INCLUDE 'VLIDORT_BOOKKEEP.VARS'

  ! Include files of result variables
  INCLUDE 'VLIDORT_RESULTS.VARS'
  INCLUDE 'VLIDORT_L_RESULTS.VARS'
  
  !INTEGER, PARAMETER :: ngas=1 
  
  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER, INTENT(IN)  :: nw, nl, nos, n0alb, nostk, n0wfc
  LOGICAL, INTENT(IN)  :: do_ozwf, do_albwf, do_tmpwf, do_o3shi, do_taodwf, vary_sfcalb, &
       do_twaewf, do_saodwf, do_cfracwf, do_codwf, do_ctpwf, do_sprswf, do_so2zwf
  INTEGER, INTENT(OUT) :: errstat
  INTEGER, DIMENSION(n0alb), INTENT(IN)        :: albpmax, albpmin
  INTEGER, DIMENSION(n0wfc), INTENT(IN)        :: wfcpmax, wfcpmin
  LOGICAL, DIMENSION(nl), INTENT(IN)           :: ozvary
  REAL (KIND=dp), DIMENSION(nw),  INTENT(IN)   :: waves, walb0s
  REAL (KIND=dp), DIMENSION(nw, nostk),  INTENT(OUT)   :: rad, albwf, cfracwf, o3shiwf, &
       codwf, ctpwf, taodwf, twaewf, saodwf, sprswf, so2zwf
  REAL (KIND=dp), DIMENSION(numwin, nos), INTENT(IN)   :: o3shi
  REAL (KIND=dp), DIMENSION(nl),  INTENT(IN)           :: ozprof, tprof
  REAL (KIND=dp), DIMENSION(nw, nl, nostk),INTENT(OUT) :: ozwf, tmpwf
  REAL (KIND=dp), DIMENSION(n0alb), INTENT(IN)         :: albarr 
  REAL (KIND=dp), DIMENSION(n0wfc), INTENT(IN)         :: wfcarr 
  REAL (KIND=dp), INTENT(IN)                           :: sza, vza, aza
  
  ! =======================
  ! Local variables
  ! =======================
  LOGICAL :: problems, ex, do_clouds, do_fozwf, do_faerwf, do_fraywf, openfileflag
  INTEGER :: status_inputread, status_inputcheck, status_calculation, ncheckmessages, nreadmessages, &
       nf, na, ic, iw, i, j, k, ii, kk, jj, jk, idum, low, hgh, fidx, lidx, nvza, naza, nz1, tempi, &
       nfgas, nstep, istk, npolmod, ipol, nsprs, npts, nw0
  INTEGER :: ozwfidx, aodwfidx, twaewfidx, codwfidx, sprswfidx, raywfidx
  LOGICAL, DIMENSION(nflay)                    :: cldmsk, varyprof!, aermsk
  REAL (KIND=dp)                               :: lamda, xg, frac, toz, temp, aodscl, waerscl
  REAL (KIND=dp), DIMENSION(0:nflay)           :: ozs, delps
  REAL (KIND=dp), DIMENSION(nflay)             :: cldsca, cldext, cldasy, eta, cldext0, aersca, aerext, aerasy
  REAL (KIND=dp), DIMENSION(0:nmom,maxgksec,nflay) :: cldmoms, aermoms
  REAL (KIND=dp), DIMENSION(0:nl)              :: cumoz
  REAL (KIND=dp), DIMENSION(nl)                :: oztmpwf, tmptwf
  REAL (KIND=dp), DIMENSION(nw)                :: albs, tmpalbs, gshiwf, wfcs
  REAL (KIND=dp), DIMENSION(nw, 2, nostk)      :: polerr
  REAL (KIND=dp), DIMENSION(nw, nostk)         :: radclr, radcld
  REAL (KIND=dp), DIMENSION(n_radwvl_sav)      :: delpos, swaves, delshi
  INTEGER,        DIMENSION(3)                 :: vind
  REAL (KIND=dp), DIMENSION(2, nostk)          :: radclrcld
  REAL (KIND=dp), DIMENSION(nw, nflay, nostk)  :: fozwf, faerwf, faerswf, fcodwf, fsprswf, fraywf
  REAL (KIND=dp), DIMENSION(nw, nflay)         :: pfozwf, pfaerwf, pfaerswf, pfcodwf, pfsprswf, pfraywf
  REAL (KIND=dp), DIMENSION(nw, nl)            :: pozwf, ptmpwf
  REAL (KIND=dp), DIMENSION(nw)                :: prad, palbwf, pctpwf, pcfracwf
  REAL (KIND=dp), DIMENSION (max_spec_pts)     :: tmprefwav, tmprefspec

  INTEGER, DIMENSION (nallgas)                  :: gasin   
  INTEGER, DIMENSION (5)                        :: tmp_gasidxs = (/1, 4, 5, 7, 3/)
  REAL (KIND=dp), DIMENSION(nw, nallgas, nflay) :: allcrs
  REAL (KIND=dp), DIMENSION(nallgas, nflay)     :: alleta, allcol  
  REAL (KIND=dp), DIMENSION(nw, nflay) :: deltau, delsca, delray, delo3abs, delabs
  
!  Exception handling for VLIDORT Model Calculation. New code, 13 October 2010
  CHARACTER(LEN=maxchlen), DIMENSION(0:MAX_MESSAGES) :: checkmessages, checkactions, &
       readmessages, readactions
  CHARACTER(LEN=maxchlen) :: message, trace_1, trace_2, trace_3

  LOGICAL, SAVE                                 :: first = .TRUE.
  INTEGER, SAVE                                 :: nz, faer_lvl, npolcorr, nradcal
  REAL (KIND=dp), DIMENSION(0:mflay),                       SAVE :: ts
  !REAL (KIND=dp), DIMENSION(max_fit_pts, mflay),            SAVE :: aersca, aerext, aerasy 
  !REAL (KIND=dp), DIMENSION(max_fit_pts,0:maxmom,maxgksec,mflay), SAVE :: aermoms
  REAL (KIND=dp), DIMENSION(3, max_fit_pts),                SAVE :: abscrs_qtdepen
  REAL (KIND=dp), DIMENSION(max_fit_pts, mflay),            SAVE :: abscrs, dads, dadt
  LOGICAL, DIMENSION(mflay),                                SAVE :: aermsk
  LOGICAL, DIMENSION(max_fit_pts),                          SAVE :: do_radcals, do_polcorrs
  INTEGER, DIMENSION(max_fit_pts),                          SAVE :: polcorr_idxs, radcal_idxs
  REAL (KIND=dp), DIMENSION(max_fit_pts, mflay),            SAVE :: so2crs
  INTEGER,                                                  SAVE :: so2idx, so2vidx, o4idx, so2crsidx

  ! Will become an input parameter later in ozprof.inp
  LOGICAL :: do_ssfullb295 = .FALSE.
  LOGICAL :: use_so2dtcrs  = .true.  ! Use SO2 cross section at different temperatures

  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=15), PARAMETER :: modulename = 'LIDORT_PROF_ENV'

  errstat = pge_errstat_ok

  IF (first) THEN
     status_inputcheck  = vlidort_success
     status_calculation = vlidort_success
     status_inputread   = vlidort_success

     ! ======================= Read LIDORT Control Input ==========================
     !xliu: 03/07/2011, switch VLIDORT from vv2p4 to vv2p4RTC
     !CALL VLIDORT_L_INPUT_MASTER ('INP/vlidort_control.inp', &
     !     'o3prof_lidort_error', status_inputread)

     openfileflag   = .FALSE.
     CALL VLIDORT_L_INPUT_MASTER ( '../run/conf/INP/vlidort_control_vv2p4RTC.inp', &
           status_inputread, nreadmessages, readmessages, readactions )

     IF (status_inputread .NE. vlidort_success) THEN
        OPEN(vlidort_errunit, file='o3prof_lidort_input_error.log', status = 'unknown')
        WRITE(vlidort_errunit, *)' FATAL:   Wrong input from VLIDORT input file-read'
        WRITE(vlidort_errunit, *)'  ------ Here are the messages and actions '
        WRITE(vlidort_errunit,'(A,I3)')'    ** Number of messages = ', nreadmessages
        DO i = 1, nreadmessages
           nf = LEN(readmessages(i))
           na = LEN(readactions(i))
           WRITE(vlidort_errunit,'(A,I3,A,A)')'Message # ', i,' : ',readmessages(i)(1:nf)
           WRITE(vlidort_errunit,'(A,I3,A,A)')'Action  # ', i,' : ',readactions(i)(1:na)
        ENDDO
        CLOSE(vlidort_errunit)
               
        WRITE(*, *) modulename, ': Problems encountered with input read!!!'
        errstat = pge_errstat_error; RETURN
     ENDIF

     ngreek_moments_input = nmom
     earth_radius = rearth
     
     IF (.NOT. aerosol .AND. (.NOT. has_clouds .OR. do_lambcld)) THEN
        DO_SSCORR_TRUNCATION = .FALSE.; DO_DELTAM_SCALING = .FALSE.
        DO_RAYLEIGH_ONLY = .TRUE.; DO_SOLUTION_SAVING = .FALSE.; DO_BVP_TELESCOPING = .FALSE.
     ENDIF

     DO i = 1, ngas
        IF (gasidxs(i) == so2v_idx)   so2vidx = i  
        IF (gasidxs(i) == so2_idx)    so2idx  = i    
        IF (gasidxs(i) == o2o2_idx)    o4idx  = i
     ENDDO

     IF (fgasidxs(so2vidx) <= 0 .AND. fgasidxs(so2idx) <= 0) use_so2dtcrs = .FALSE.

     ! Could not handle two different wavelength shifts for the same cross sections
     IF (fgassidxs(so2vidx) > 0 .AND. fgassidxs(so2idx) > 0 .AND. use_so2dtcrs) use_so2dtcrs = .FALSE.

     first = .FALSE.
  ENDIF
 
  ! ============= Overridden some control and atmospheric variables ============== 
   IF (num_iter == 0 ) THEN
     n_szangles = 1; szangles(1) = sza 
     IF (sza >= 90.0 .OR. sza < 0) THEN
        WRITE(*, *) modulename, ' : SZA is >= 90 or < 0 !!!'
        errstat = pge_errstat_error; RETURN
     ENDIF
     
     n_user_vzangles = 1; n_user_relazms = 1
     user_vzangles(1) = vza; user_relazms(1) = aza
     
     !IF (.NOT. do_multi_vza) THEN
     !ELSE
     !   n_user_streams = 3; n_user_relazms = 3
     !   user_angles_input(1:3) = ABS(vza_atm); user_relazms(1:3) = aza_atm
     !   ! vza will be sorted automatically in LIDORT, need to get order index
     !   vind(1) = MINVAL(MINLOC(ABS(vza_atm)))
     !   vind(3) = MINVAL(MAXLOC(ABS(vza_atm)))
     !   vind(2) = 6 - vind(1) - vind(3)
     !   IF (vind(2) < 1 .OR. vind(2) > 3) THEN
     !      WRITE(*, *) modulename, ' : Sth. wrong in view geometry!!!'
     !      errstat = pge_errstat_error; RETURN
     !   ENDIF
     !
     !   vind((/vind(1), vind(2), vind(3)/)) = (/1, 2, 3/)
     !   naza =3; nvza = 3
     !ENDIF
     nz = nflay;   nlayers = nz

     IF (nz > maxlayers) THEN
        WRITE(*, *) modulename, ' : # of layers exceeded allowed !!!',nz, maxlayers

        errstat = pge_errstat_error; RETURN
     ENDIF
     IF (nl > nlayers) THEN
        WRITE(*, *) modulename, ' : Coarse grids cannot be finer than fine grids!!!'
        errstat = pge_errstat_error; RETURN
     ENDIF
       
     ! Set height grid for doing Chapman Function Calculation 
     height_grid(0:nlayers) = fzs(0:nz)
     geometry_specheight  = the_surfalt
     !print *, the_surfalt, height_grid(nlayers), fzs(nz)

     ! set the first aerosol layer (from TOA down) 1:faer_lvl-1 (without aerosols)
     faer_lvl = 1;   IF (.NOT. strat_aerosol) faer_lvl = nup2p(ntp) + 1
    
     IF (aerosol) THEN 
        aermsk = .TRUE.
        IF (.NOT. strat_aerosol) aermsk(1:faer_lvl-1) = .FALSE.
     ELSE
        aermsk = .FALSE.
     ENDIF
         
     ts(1:nz) = (fts(1:nz) + fts(0:nz-1)) / 2.0
     !WRITE(*, '(12F8.2)') ts(1:nz)
     !CALL get_efft(nz, fzs(0:nz), fozs(1:nz), fts(0:nz), ts(1:nz), errstat)
     !WRITE(*, *)
     !WRITE(*, '(12F8.2)') ts(1:nz)
  ENDIF

  ! =================== Interpolate Ozone, T to fine grids =======================
  IF (nz /= nl) THEN
     DO i = 1, nl
        varyprof(nup2p(i-1)+1:nup2p(i)) = ozvary(i)
        ozs(nup2p(i-1)+1:nup2p(i)) = fozs(nup2p(i-1)+1:nup2p(i)) * ozprof(i) / &
             SUM(fozs(nup2p(i-1)+1:nup2p(i)))
     ENDDO

     IF (nt_fit > 0) THEN
        
        DO i = 1, nl
           atmosprof(3, i) = tprof(i) * 2.0 - atmosprof(3,  i-1)
        ENDDO
        
        CALL BSPLINE(atmosprof(2, 0:nl), atmosprof(3, 0:nl), &
             nl+1, fzs(0:nz), ts(0:nz), nz+1, errstat)
        IF (errstat < 0) THEN
           WRITE(*, *) modulename, ' : BSPLINE error, errstat = ', errstat; RETURN
        ENDIF
        ts(1:nz) = (ts(1:nz) + ts(0:nz-1)) / 2.0
     ENDIF
 ELSE
    ozs(1:nz) = ozprof(1:nl);    varyprof(1:nz) = ozvary(1:nl)
 ENDIF
 
 nz1 = nfsfc - 1 
 ! Update aerosol fields: first AOD
 IF (do_taodwf .AND. num_iter > 0) THEN
    aodscl = fitvar_rad(taodind) / tropaod(actawin)
    tropaod(1:actawin) = tropaod(1:actawin) * aodscl
    gaext(1:actawin, nup2p(ntp)+1:nz1) = gaext(1:actawin, nup2p(ntp)+1:nz1) * aodscl
    
    IF (.NOT. do_twaewf) THEN ! Single scattering albedo does not change
       tropsca(1:actawin) = tropsca(1:actawin) * aodscl
       gasca(1:actawin, nup2p(ntp)+1:nz1) = gasca(1:actawin, nup2p(ntp)+1:nz1) * aodscl
    ENDIF
 ELSE
    aodscl = 1.0
 ENDIF

 IF (do_twaewf .AND. num_iter > 0) THEN
    waerscl = fitvar_rad(twaeind) / tropwaer(actawin)   ! Scale single scattering albedo
    tropwaer(1:actawin) = tropwaer(1:actawin) * waerscl
    tropsca(1:actawin)  = tropsca(1:actawin) * waerscl * aodscl
    gasca(1:actawin, nup2p(ntp)+1:nz1) = gasca(1:actawin, nup2p(ntp)+1:nz1) * waerscl * aodscl
 ENDIF

 IF (do_saodwf .AND. num_iter > 0) THEN
    aodscl = fitvar_rad(saodind) / strataod(actawin)
    strataod(1:actawin) = strataod(1:actawin) * aodscl
    stratsca(1:actawin) = stratsca(1:actawin) * aodscl
    gaext(1:actawin, 1:nup2p(ntp)) = gaext(1:actawin, 1:nup2p(ntp)) * aodscl
    gasca(1:actawin, 1:nup2p(ntp)) = gasca(1:actawin, 1:nup2p(ntp)) * aodscl
 ENDIF

 IF (do_sprswf .AND. num_iter > 0) THEN
    temp = (fitvar_rad(sprsind) - fps(nup2p(nsfc-1)))/ (fps(nz1) - fps(nup2p(nsfc-1)))
    frhos(nup2p(nsfc-1)+1:nz1) = frhos(nup2p(nsfc-1)+1:nz1) * temp
    delps(nup2p(nsfc-1)+1:nz1) = (fps(nup2p(nsfc-1)+1:nz1)-fps(nup2p(nsfc-1):nz1-1))  * temp
    DO i = nup2p(nsfc-1)+1, nz1
       fps(i) = fps(i-1) + delps(i)
    ENDDO
 ENDIF
  
 ! albedo array
 IF (.NOT. vary_sfcalb) THEN
    DO i = 1, n0alb
       albs(albpmin(i):albpmax(i)) = albarr(i)
    ENDDO
 ELSE
    albs(1:nw) = walb0s(1:nw)
 ENDIF
            
 DO i = 1, n0wfc
    wfcs(wfcpmin(i):wfcpmax(i)) = wfcarr(i)
 ENDDO

  ! Determinine atmospheric weighting functions to be calculated
  IF (do_ozwf .OR.  do_tmpwf .OR. do_o3shi .OR. do_taodwf .OR. &
       do_twaewf .OR. do_saodwf .OR. do_codwf .OR. do_sprswf )  do_atmos_linearization = .TRUE.
  do_surface_linearization = do_albwf    
  IF (do_atmos_linearization .OR. do_albwf) THEN
     do_simulation_only =    .FALSE.;    do_linearization = .TRUE.
  ELSE 
     do_simulation_only =    .TRUE.;     do_linearization = .FALSE.
  ENDIF
  
  i = 0;  ozwfidx = 0; aodwfidx = 0; twaewfidx = 0; codwfidx = 0; sprswfidx = 0
  layer_vary_flag(1:nz) = .FALSE.
  layer_vary_number(1:nz) = 0
  IF (do_ozwf .OR. do_tmpwf .OR. do_o3shi) THEN
     i = i + 1;  ozwfidx = i
     profilewf_names(i) = 'ozone volume mixing ratio------'

     layer_vary_flag(1:nz) = varyprof(1:nz)
     WHERE (varyprof(1:nz) )
        layer_vary_number(1:nz) = layer_vary_number(1:nz) + 1
     ENDWHERE
  ENDIF
  IF ( do_taodwf .OR. do_saodwf) THEN
     i = i + 1; aodwfidx = i
     profilewf_names(i) = 'aerosol extinction coefficient-'

     IF (do_taodwf) THEN
        layer_vary_flag(nup2p(ntp)+1:nz1) = .TRUE.
        layer_vary_number(nup2p(ntp)+1:nz1) = layer_vary_number(nup2p(ntp)+1:nz1) + 1
     ENDIF

     IF (do_saodwf) THEN
        layer_vary_flag(1:nup2p(ntp)) = .TRUE.
        layer_vary_number(1:nup2p(ntp)) = layer_vary_number(1:nup2p(ntp)) + 1
     ENDIF
  ENDIF
  IF ( do_twaewf) THEN
     i = i + 1; twaewfidx = i
     profilewf_names(i) = 'aerosol scattering coefficient-'
     layer_vary_flag(nup2p(ntp)+1:nz1) = .TRUE.
     layer_vary_number(nup2p(ntp)+1:nz1) = layer_vary_number(nup2p(ntp)+1:nz1) + 1
  ENDIF
  IF ( do_codwf) THEN
     i = i + 1; codwfidx = i
     profilewf_names(i) = 'cloud extinction coefficient---'
     layer_vary_flag(nctp:ncbp) = .TRUE.
     layer_vary_number(nctp:ncbp) = layer_vary_number(nctp:ncbp) + 1
  ENDIF
  do_fraywf = .FALSE.
  IF ( do_sprswf  .OR. (.NOT. use_effcrs .AND. nw > 1) ) THEN
     ! Need to use jacobians wrt rayleigh OD to perform interpolation

     i = i + 1; sprswfidx = i; raywfidx = i
     profilewf_names(i) = 'rayleigh optical thickness-----'
     IF (.NOT. use_effcrs) THEN
        do_fraywf = .TRUE.
        layer_vary_flag(1:nz1) = .TRUE.
        layer_vary_number(1:nz1) = layer_vary_number(1:nz1) + 1
     ELSE
        layer_vary_flag(nup2p(nsfc-1)+1:nz1) = .TRUE.
        layer_vary_number(nup2p(nsfc-1)+1:nz1) = layer_vary_number(nup2p(nsfc-1)+1:nz1) + 1
     ENDIF
  ENDIF

  n_totalatmos_wfs = i
  ! xliu, 03/8/11, has to be initialized in v2p4RTC
  n_surface_wfs = 1   

  IF (n_totalatmos_wfs > 1) THEN
     WHERE(layer_vary_number(1:nz) > 0) 
        layer_vary_number(1:nz) = n_totalatmos_wfs
     ENDWHERE
  ENDIF

  do_fozwf = .FALSE.
  IF (do_ozwf .OR. do_tmpwf .OR. do_o3shi .OR. (.NOT. use_effcrs .AND. nw > 1)) do_fozwf = .TRUE.
  do_faerwf = .FALSE.
  IF (do_taodwf .OR. do_saodwf) do_faerwf = .TRUE.
      
!  print *, n_totalatmos_wfs, nz
!  print *, layer_vary_flag(1:nz)
!  print *, layer_vary_number(1:nz)
  
  ! ==================== Get Ozone Absorption Cross Section ====================  
  
  IF (nw > 1) THEN
     allcol(1, 1:nz1) = ozs(1:nz1) * du2mol
     nfgas = 1; gasin(1) = 1
     
     DO k = 1, ngas
        IF (fgasidxs(k) > 0) THEN
           nfgas = nfgas + 1; gasin(nfgas) = nfgas
           IF (use_so2dtcrs) so2crsidx = nfgas
           
           ! molecules cm^-2, but normalized by refspec_norm(gasidxs(k))
           ! allcol / refspec_norm will be molecules cm^-2
           allcol(nfgas, 1:nz1) = mgasprof(k, 1:nz1) * tracegas(k, 4) / mgasprof(k, nz+1)
        ENDIF
     ENDDO
              
     IF (use_effcrs) THEN
        IF (nos > 0) THEN
           delshi = 0.0

           IF (do_subfit) THEN
              fidx = 1
              DO j = 1, numwin
                 lidx = fidx + nradpix_sav(j) - 1
                 delpos(fidx:lidx) =  radwvl_sav(fidx:lidx) - (radwvl_sav(fidx) + radwvl_sav(lidx)) / 2.0
                 IF (osfind(j, 1) > 0) delshi(fidx:lidx) =  o3shi(j, 1) 
                 
                 DO i = 2, nos
                    IF (osfind(j, i) > 0) delshi(fidx:lidx) = delshi(fidx:lidx)  + &
                         o3shi(j, i) * delpos(fidx:lidx) ** (i-1)
                 ENDDO
                 fidx = lidx + 1
              ENDDO
           ELSE
              IF (oswins(1, 1) == 1) THEN
                 fidx = 1
              ELSE
                 fidx = SUM(nradpix(1: oswins(1, 1)-1)) + 1
              ENDIF
              lidx = SUM(nradpix(1: oswins(1, 2)))
              
              delpos(fidx:lidx) =  radwvl_sav(fidx:lidx) - (radwvl_sav(fidx) + radwvl_sav(lidx)) / 2.0
              IF (osfind(1, 1) > 0) delshi(fidx:lidx) =  + o3shi(1, 1) 
              
              DO i = 2, nos  
                 IF (osfind(1, i) > 0) delshi(fidx:lidx) = delshi(fidx:lidx)  + &
                      o3shi(1, i) * delpos(fidx:lidx) ** (i-1)
              ENDDO
           ENDIF
        
           swaves = radwvl_sav - delshi        
        ELSE
           swaves = radwvl_sav
        ENDIF

        ! Get ozone absorption coefficients
        IF (ANY(ABS(delshi) >= 2.0)) THEN
           WRITE(*, *) modulename, ' : Ozone wavelength shifts are too large!!!'
           errstat = pge_errstat_error; RETURN
        ENDIF

        ! Get temperature-dependent ozone cross section at instrumental spectral resolution
        IF (num_iter == 0 .OR. (do_o3shi .AND. nw > 1) ) THEN
           CALL GETABS_CRS(swaves, n_radwvl_sav, nw, 1, nz1, ts(1:nz1), abscrs(1:nw, 1:nz1), &
                do_o3shi, do_tmpwf, dads(1:nw, 1:nz1), dadt(1:nw, 1:nz1), problems, &
                abscrs_qtdepen(1:3, 1:nw))
           
           IF (problems) THEN
              WRITE(*, *) modulename, ' : Problems in reading ozone absorption !!!'
              errstat = pge_errstat_error; RETURN
           ENDIF
        ENDIF      
        allcrs(1:nw, 1, 1:nz1) = abscrs(1:nw, 1:nz1)
        
        ! Get Rayleigh scattering coefficients and deplolarizaiton factor
        !IF (num_iter == 0) THEN 
        !   CALL GET_ALL_RAYCOF_DEPOL(nw, waves, raycof(1:nw), depol(1:nw))
        IF (num_iter == 0) THEN 
           CALL GET_ALL_RAYCOF_DEPOL1(n_radwvl_sav, radwvl_sav, nw, raycof(1:nw), depol(1:nw), problems)
        ENDIF
        
        IF ( num_iter == 0 .AND. (fgasidxs(so2idx) > 0 .OR. fgasidxs(so2vidx) > 0) .AND. use_so2dtcrs) THEN
           CALL GETSO2_CRS(radwvl_sav, n_radwvl_sav, nw, nz1, ts(1:nz1), so2crs(1:nw, 1:nz1), problems)
           IF (problems) THEN
              WRITE(*, *) modulename, ' : Problems in reading SO2 absorption !!!'
              errstat = pge_errstat_error; RETURN
           ENDIF
        ENDIF
        
        nfgas = 1
        DO k = 1, ngas
           IF (fgasidxs(k) > 0) THEN
              nfgas = nfgas + 1; gasin(nfgas) = nfgas
                            
              IF ((k /= so2idx .AND. k /= so2vidx) .OR. .NOT. use_so2dtcrs) THEN
                 IF (fgassidxs(k) > 0 ) THEN
                    npts = n_refspec_pts(gasidxs(k))
                    tmprefwav(1:npts) = refspec_orig_data(gasidxs(k), 1:npts, 1)
                    tmprefspec(1:npts) = refspec_orig_data(gasidxs(k), 1:npts, 3)
                    fidx = MINVAL(MINLOC(waves(1:nw), MASK=(waves(1:nw) >= &
                         tmprefwav(1) + 0.1 .AND. waves(1:nw) <= tmprefwav(npts) - 0.1)))
                    lidx = MINVAL(MAXLOC(waves(1:nw), MASK=(waves(1:nw) >= &
                         tmprefwav(1) + 0.1 .AND. waves(1:nw) <= tmprefwav(npts) - 0.1)))
                    
                    IF (lidx > fidx .AND. lidx > 0 .AND. fidx > 0) THEN 
                       temp = fitvar_rad(rmask_fitvar_rad(fgassidxs(k)))
                       CALL BSPLINE1(  tmprefwav(1:npts) - temp, tmprefspec(1:npts), npts, &
                            waves(fidx:lidx), allcrs(fidx:lidx, nfgas, 1),  gshiwf(fidx:lidx), lidx-fidx+1, errstat) 
                       
                       !print *, k, gasidxs(k), fgasidxs(k), gassidxs(k), fgassidxs(k), fitvar_rad(gassidxs(k)), npts
                       !WRITE(90, '(F10.4, 2D14.6)') ((waves(i), database_shiwf(gasidxs(k), refidx(i)), gshiwf(i)), i=fidx, lidx)
                       
                       database_shiwf(gasidxs(k), refidx(fidx:lidx)) = gshiwf(fidx:lidx)
                       database(gasidxs(k), refidx(fidx:lidx)) = allcrs(fidx:lidx, nfgas, 1)               
                       IF (errstat < 0) THEN
                          WRITE(*, *) modulename, ' : BSPLINE error, errstat = ', errstat; RETURN
                       ENDIF
                    ENDIF
                 ELSE
                    allcrs(1:nw, nfgas, 1) = database_save(gasidxs(k), refidx(1:nw))
                 ENDIF
                 
                 DO i = 2, nz1
                    allcrs(1:nw, nfgas, i) = allcrs(1:nw, nfgas, 1)
                 ENDDO
              ELSE
                 allcrs(1:nw, nfgas, 1:nz1) = so2crs(1:nw, 1:nz1) / refspec_norm(gasidxs(k))
              ENDIF
           ENDIF
        ENDDO
        
        IF (do_radinter) THEN
           IF (num_iter == 0) THEN
              do_radcals(1:nw) = .FALSE.
              fidx = 1; k = nz1 / 2
              
              DO iw = 1, numwin
                 lidx = fidx + nradpix(iw) - 1
                 
                 ! Do radiative transfer calculations at local maxima/minima
                 ! always do radiative transfer calculation at end points
                 IF (band_selectors(iw) == 2) THEN
                    DO i = fidx + 1, lidx - 1
                       IF (abscrs(i, k) > abscrs(i-1, k) .AND. abscrs(i, k) > abscrs(i+1, k)) &
                            do_radcals(i) = .TRUE.
                       IF (abscrs(i, k) < abscrs(i-1, k) .AND. abscrs(i, k) < abscrs(i+1, k)) &
                            do_radcals(i) = .TRUE.
                       !IF (abscrs(i, nz1) > abscrs(i-1, nz1) .AND. abscrs(i, nz1) > abscrs(i+1, nz1)) &
                       !    do_radcals(i) = .TRUE.
                       !IF (abscrs(i, nz1) < abscrs(i-1, nz1) .AND. abscrs(i, nz1) < abscrs(i+1, nz1)) &
                       !    do_radcals(i) = .TRUE.
                    ENDDO
                    
                    nstep = 5; do_radcals(fidx) = .TRUE.; do_radcals(lidx-1:lidx) = .TRUE.                
                 ELSE IF (band_selectors(iw) == 1) THEN
                    nstep = 3; do_radcals(fidx) = .TRUE.; do_radcals(lidx) = .TRUE.
                 ELSE
                    nstep = 1
                 ENDIF
                 
                 DO i = fidx + 1, lidx - 1, nstep
                    do_radcals(i) = .TRUE.
                 ENDDO
                 
                 fidx = lidx + 1        
              ENDDO
              
              nradcal = 0
              DO i = 1, nw
                 IF (do_radcals(i)) THEN
                    nradcal = nradcal + 1; radcal_idxs(nradcal) = i
                 ENDIF
              ENDDO
           ENDIF
        ELSE
           IF (num_iter == 0) THEN
              do_radcals(1:nw) = .TRUE.; nradcal = nw
              DO i = 1, nw
                 radcal_idxs(i) = i
              ENDDO
           ENDIF
        ENDIF
        
        !WRITE(77, *) nw, nz1
        !WRITE(77, '(30F8.3)') ts(1:nz1)
        !WRITE(77, '(30F8.3)') ozs(1:nz1)
        !DO i = 1, nw
        !   WRITE(77, '(60D14.6)') waves(i), abscrs(i, 12) !, dads(i, 1:nz1)
        !ENDDO  
        !CLOSE (77)
        
        !WRITE(91, *) nz1
        !WRITE(91, '(30D14.6)')  frhos(1:nz1)
        !WRITE(91, '(30D14.6)')  fzs(0:nz1)
        !WRITE(91, *) nw
        !DO i = 1, nw
        !   WRITE(91, '(F10.4, D14.6)') waves(i), raycof(i)
        !ENDDO
     ELSE   ! .NOT. use_effcrs
        ! O3/SO2 (use_so2dtcrs=.TRUE.) cross section: if do_tmpwf = .FALSE. and do_o3shi is false, 
        ! just need to get once for each retrieval
        ! Other trace gas cross section: just need to get it once for all the retrievals if no shifts 
        CALL GET_HRES_GASCRS_RAY(nw, waves, nz1, ts(1:nz1), do_o3shi, o3shi, do_tmpwf, nfgas, &
             use_so2dtcrs, num_iter, allcrs(1:nw, 1:nfgas, 1:nz1), allcol(1:nfgas, 1:nz1), frhos(1:nz1), &
             abscrs_qtdepen(1:3, 1:nw), raycof(1:nw), depol(1:nw), errstat)
          
        abscrs(1:ncalcp, 1:nz1) = allcrs(1:ncalcp, 1, 1:nz1)
        IF (use_so2dtcrs) so2crs(1:ncalcp, 1:nz1) = allcrs(1:ncalcp, so2crsidx, 1:nz1)
     
        IF (num_iter == 0) THEN
           do_radcals(1:nw) = .FALSE.
           nradcal = ncalcp
           do_radcals(1:ncalcp) = .TRUE.
           radcal_idxs(1:ncalcp) = (/(i, i=1, ncalcp)/)
        ENDIF
                 
     ENDIF

     !IF (wrtozcrs) THEN
     !   WRITE(91, *) nw
     !   DO iw = 1, nw 
     !      WRITE(91, '(F10.4, 6D14.6)') waves(iw), abscrs_qtdepen(1:3, iw), raycof(iw), depol(iw) !,  &
     !           !database(6, refidx(iw))*refspec_norm(6), database(9, refidx(iw))*refspec_norm(9), &
     !           !database(10, refidx(iw))*refspec_norm(10), database(12, refidx(iw))*refspec_norm(12)
     !   ENDDO
     !   
     !   WRITE(91, *) nz1
     !   WRITE(91, '(2D14.6)') fps(0), fzs(0)
     !   DO i = 1, nz1 
     !      WRITE(91, '(10D14.6)') fps(i), fzs(i), ts(i), frhos(i), allcol(1, i), allcol(2, i)/refspec_norm(6), &
     !           allcol(3, i)/refspec_norm(9), allcol(4, i)/refspec_norm(10), allcol(5, i)/refspec_norm(12) 
     !   ENDDO
     !ENDIF
     
 ELSE
     do_radcals(1) = .TRUE.; nradcal = 1; radcal_idxs(1) = 1

     ! o3 absorption coefficient at 370.2 nm with TOMS FWHM
     dads(1, 1:nz1) = 0.0; dadt(1, 1:nz1) = 0.0

     ! Weighted by solar flux
     !abscrs(1, 1:nz) = 9.1231787D-24 + (ts(1:nz) - 273.15) * 1.9005502D-25 + &
     !     (ts(1:nz) - 273.15)**2.0 * 1.2275286D-27     
     !raycof(1) = 2.3184501D-26; depol(1) = 0.030247913D0

     nfgas = 7
     CALL GET_ALB_OZCRS_RAY(nz1, ts(1:nz1), nfgas, allcrs(1, 1:nfgas, 1:nz1), raycof(1), depol(1), problems)
     IF (problems) THEN
        WRITE(*, *) modulename, ' : Problems in reading O3 XSec for determining Fc!!!'
        errstat = pge_errstat_error; RETURN
     ENDIF

     !CALL GET_ALL_RAYCOF_DEPOL(1, 360.D0, raycof(1), depol(1))
     !print *, raycof(1), depol(1)
 
     ! Ozone
     nfgas = 6;  gasin(1) = 1
     allcol(1, 1:nz1) = ozs(1:nz1) * du2mol

     DO i = 1, 4
        allcol(i+1, 1:nz1) = mgasprof(tmp_gasidxs(i), 1:nz1)
        gasin(i+1) = i + 1
     ENDDO
     
     ! O2-O2 Concentration (Oxygen: 20.95%)
     allcol(6, 1:nz1) = ( frhos(1:nz1) * 0.2095 ) ** 2.0 / (fzs(0:nz1-1) - fzs(1:nz1)) / 1.0D5
     gasin(6) = 6
  ENDIF
    
  ! Initialize aerosol/cloud property profiles
  IF (num_iter == 0) THEN
     aersca(1:nz1) = 0.0; aerext(1:nz1) = 0.0
     aerasy(1:nz1) = 0.0;  aermoms(0:nmom, 1:ngksec, 1:nz1) = 0.0
  ENDIF

  cldmsk = .FALSE.; IF (nctp /= 0) cldmsk(nctp:ncbp)=.TRUE.
  cldsca=0.0; cldext=0.0; cldasy=0.0;  cldmoms = 0.0
  
  ! Initialize output variables  
  rad = 0.0; radclr = 0.0; radcld = 0.0
  cfracwf = 0.0
  IF ( do_fozwf )  fozwf   = 0.0
  IF ( do_albwf)   albwf   = 0.0
  IF ( do_codwf )  fcodwf  = 0.0
  IF ( do_sprswf ) fsprswf  = 0.0
  IF ( do_fraywf ) fraywf  = 0.0
  IF ( do_ctpwf )  ctpwf   = 0.0
  IF ( do_faerwf)  faerwf  = 0.0
  IF ( do_twaewf ) faerswf = 0.0
  
  ! senstitivity of absorption cross section to temperature, used 
  ! for calculating temperature wf directly with LIDORT
  !eta = 0.0             ! dummy variable here
  alleta = 0.0

  polerr = 1.0
  IF ( polcorr == 2) THEN 
     IF (nwfc > 0) THEN
        STOP 'Polarization correction + variable fc: not implemented!!!'
     ENDIF
     tmpalbs = albs 
     ! clear-sky
     IF (the_cfrac /= 1.0) THEN
        toz = SUM(ozs(1:nz1))
        
        WHERE(tmpalbs >= 1.0)
           tmpalbs = 0.999
        ENDWHERE
        
        WHERE(tmpalbs < 0.0)
           tmpalbs = 0.001
        ENDWHERE
        
!        CALL lup_polerror(toz, fps(nz1), nw, waves, tmpalbs, polerr(1:nw, 1, 1), errstat)
        IF (errstat == pge_errstat_error) RETURN
        polerr(1:nw, 1, 1) = 1.0 + polerr(1:nw, 1, 1) / 100.0
     ENDIF
     
     ! cloudy part
     IF (the_cfrac /= 0.0) THEN
        toz = SUM(ozs(1:nctp-1))
        IF (do_lambcld) THEN
           tmpalbs = lambcld_refl
        ELSE
           ! approxmation of cod as cloud albedo for polarization correction
           ! Follows Kokhanovsky et al., 2003
           tmpalbs = 1.0 - 1.0 / (1.072 + 0.1125 * the_cod)
        ENDIF
        
        WHERE(tmpalbs >= 1.0)
           tmpalbs = 0.999
        ENDWHERE
        
        WHERE(tmpalbs < 0.0)
           tmpalbs = 0.001
        ENDWHERE
        
!        CALL lup_polerror(toz, fps(nctp-1), nw, waves, tmpalbs, polerr(1:nw, 2, 1), errstat)
        IF (errstat == pge_errstat_error) RETURN
        polerr(1:nw, 2, 1) = 1.0 + polerr(1:nw, 2, 1) / 100.0
     ENDIF
  ENDIF

  ! Determine wavelengths where exact polarization correction (NSTOKES: 4 vs 1) is calculated
  ! In UV1 (or between 270 and 310 nm): ~292 nm, ~298 nm, ~300 nm, ~302 nm, ~304 nm, ~306 nm, last wavelength 
  ! In UV2 (or between 310 and 340 nm): first, 1/4, middle and last wavelength
  ! So exact vector LIDORT calculation is done at 11 wavelengths.
  ! This option works when radiance interpolation option is turned on
  IF (num_iter == 0) THEN
     do_polcorrs(1:nw) = .FALSE.; npolcorr = 0
     IF ( (polcorr >= 3 .AND. polcorr <= 5) .AND. nw > 1 ) THEN
        fidx = 1; idum = 0
        DO iw = 1, numwin
           IF (iw == numwin) THEN
              temp = winlim(iw, 2)
           ELSE
              temp = (winlim(iw, 2) + winlim(iw + 1, 1)) / 2.
           ENDIF
           lidx = MINVAL(MAXLOC(waves(1:nw), MASK=(waves(1:nw) < temp .AND. waves(1:nw) > 0))) 
           !lidx = fidx + nradpix(iw) - 1
           IF ( waves(lidx) <= 312.0 ) THEN
              ! Error from using single scattering is about 0.2% at 270 nm, it needs o be corrected
              !IF (do_ssfullb295) do_polcorrs(fidx) = .TRUE. 
              IF (idum == 0) idum = 1
              DO i = fidx + 1, lidx - 1
                 IF (do_radcals(i) ) THEN
                    IF ( (waves(idum) < 290. .AND. waves(i) >= 290.0) .OR.  &
                         (waves(idum) < 295. .AND. waves(i) >= 295.0) .OR.  &
                         (waves(idum) < 299. .AND. waves(i) >= 299.0) .OR.  &
                         (waves(idum) < 301. .AND. waves(i) >= 301.0) .OR.  &
                         (waves(idum) < 303. .AND. waves(i) >= 303.) .OR.  &
                         !(waves(idum) < 304. .AND. waves(i) >= 304.0) .OR.  &
                         (waves(idum) < 305. .AND. waves(i) >= 305.0) )  do_polcorrs(i) = .TRUE.
                    idum = i
                 ENDIF
              ENDDO
              do_polcorrs(lidx) = .TRUE.
           ELSE
              IF (idum == 0) idum = 1
              do_polcorrs(fidx) = .TRUE.; do_polcorrs(lidx) = .TRUE.     
              j = (fidx + lidx ) / 2; k = fidx + (j - fidx) / 3  

              jj = (j + lidx )  / 2; kk = fidx + (k - fidx) / 3
              jk = (j + k ) / 2.
                          
              DO i = fidx + 1, lidx - 1
                 IF (do_radcals(i) ) THEN
                    !IF ( (waves(idum) < waves(j) .AND. waves(i) >= waves(j)) .OR. &
                    !     (waves(idum) < waves(k) .AND. waves(i) >= waves(k)) .OR. &
                    !     (waves(idum) < waves(jj) .AND. waves(i) >= waves(jj)) .OR. &
                    !     (waves(idum) < waves(kk) .AND. waves(i) >= waves(kk)) .OR. &
                    !     (waves(idum) < waves(jk) .AND. waves(i) >= waves(jk)))  do_polcorrs(i) = .TRUE.
                    IF ( (waves(idum) < waves(j) .AND. waves(i) >= waves(j)) .OR. &
                         (waves(idum) < waves(k) .AND. waves(i) >= waves(k)))  do_polcorrs(i) = .TRUE.
                    idum = i
                 ENDIF
              ENDDO
           ENDIF
           
           fidx = lidx + 1
        ENDDO
        DO i = 1, nw
           IF ( do_polcorrs(i) ) THEN
              !IF (do_ssfullb295 .AND. npolcorr == 2 .AND. waves(1) < 290.0) THEN
              !   ! Use single scattering (ss), scalar below 295 nm
              !   ! Use ss + multiple scattering, scalar above 295 nm
              !   ! Correction: 270 (ss + scalr), 290  (ss + scalar),
              !   !             first lamda before 295 (ss + scalar)
              !   !             first lamda after  295 (scalar)
              !   do_polcorrs(i-1) = .TRUE.; npolcorr = npolcorr + 2
              !   polcorr_idxs(3) = i-1; polcorr_idxs(4) = i
              !ELSE
                 npolcorr = npolcorr + 1; polcorr_idxs(npolcorr) = i
              !ENDIF
           ENDIF
        ENDDO
        
        !print *, nradcal, npolcorr
        !print *, polcorr_idxs(1:npolcorr)
        !print *, do_radcals(polcorr_idxs(1:npolcorr))
        !print *, waves(polcorr_idxs(1:npolcorr))
        !STOP
     ENDIF
  ENDIF
  
  ! ====================== Call LIDORT and Do Post Processing =====================
  DO iw = 1, nw 
     
     IF (.NOT. do_radcals(iw) ) CYCLE   ! Radiances and weighting functions will be interpolated 
     lamda = waves(iw) 
     
     IF (nwfc > 0) the_cfrac = wfcs(iw)
      
     ! NSTOKES = 4 when
     ! (1) Vector LIDORT or 
     ! (2) on-line polarization correction (polcorr=3) at cloud wavelength or wavelengths 
     !     where exact radiances are calcualted or
     ! (3) on-line polarization correction (polcorr=4) at cloud wavelength or wavelengths 
     !     where exact radiances are calcualted in the first iteration
     IF (polcorr == 0 .OR. ((polcorr == 3 .OR. polcorr == 5) .AND. (nw == 1 .OR. do_polcorrs(iw)))  &
          .OR. (polcorr == 4 .AND. (nw == 1 .OR. (do_polcorrs(iw) .AND. num_iter == 0))) ) THEN
        NSTOKES = 3 ; NSTREAMS = 4
     ELSE
        NSTOKES = 1 ; NSTREAMS = 4
     ENDIF
     VlidortNstream = NSTREAMS
  
     IF ( lamda < 295.0 .AND. NSTOKES == 1 .AND. nw > 1 .AND. do_ssfullb295) THEN
        DO_SSFULL = .TRUE.;  DO_SSCORR_TRUNCATION = .FALSE.; DO_DELTAM_SCALING = .FALSE.
     ELSE
        DO_SSFULL = .FALSE.
        IF (aerosol .OR. (has_clouds .AND. .NOT. do_lambcld)) THEN
           DO_SSCORR_TRUNCATION = .TRUE.; DO_DELTAM_SCALING = .TRUE.
        ENDIF
     ENDIF
          
     IF ( aerosol ) THEN  !.AND. num_iter == 0 ) THEN
        ! Interpolate/extrapoalte for aerosol properties
        hgh = actawin
        DO i = 1, actawin
           IF (lamda < aerwavs(i)) THEN
              hgh = i; EXIT
           ENDIF
        ENDDO

        IF (hgh==1) hgh = 2    ! extrapolation
        low = hgh - 1
        xg = (lamda - aerwavs(low)) / (aerwavs(hgh) - aerwavs(low))
        
        aersca(faer_lvl:nfsfc-1) = gasca(low, faer_lvl:nfsfc-1) * (1.0 - xg) +  gasca(hgh, faer_lvl:nfsfc-1) * xg
        aerext(faer_lvl:nfsfc-1) = gaext(low, faer_lvl:nfsfc-1) * (1.0 - xg) +  gaext(hgh, faer_lvl:nfsfc-1) * xg
        
        !IF (num_iter == 0) THEN  ! Don't need to be updated
        IF (useasy) THEN
           aerasy(faer_lvl:nfsfc-1) = gaasy(low, faer_lvl:nfsfc-1) * (1.0 - xg) + &
                gaasy(hgh, faer_lvl:nfsfc-1) * xg
        ELSE        
           DO i = 0, nmom
              DO j = 1, ngksec
                 aermoms(i, j, faer_lvl:nfsfc-1) = &
                      gamoms(low, faer_lvl:nfsfc-1, i, j) * (1.0 - xg) +  &
                      gamoms(hgh, faer_lvl:nfsfc-1, i, j) * xg
              ENDDO
           ENDDO
        ENDIF
        !ENDIF
     ENDIF
     
     IF (has_clouds .AND. .NOT. do_lambcld) THEN
        cldext(nctp:ncbp) = (gcq(low) * (1.0 - xg) + gcq(hgh) * xg) * &
             (fzs(nctp-1 : ncbp-1) - fzs(nctp:ncbp)) * the_cbeta
        cldsca(nctp:ncbp) = cldext(nctp:ncbp)  ! w0 = 1.0
        IF (iw == 1) cldext0(nctp:ncbp) = cldext(nctp:ncbp) * gcq(actawin) / (gcq(low) * (1.0 - xg) + gcq(hgh) * xg) 
        
        IF (useasy) THEN
           cldasy(nctp:ncbp) = gcasy(low) * (1.0 - xg) + gcasy(hgh) * xg
        ELSE
           DO i = 0, nmom
              DO j = 1, ngksec
                 cldmoms(i, j, nctp:ncbp) = gcmoms(low, i, j) &
                      * (1.0 - xg) + gcmoms(hgh, i, j) * xg 
              ENDDO
           ENDDO
        ENDIF
     ENDIF
     
     IF ( do_polcorrs(iw) .AND. ((polcorr == 3 .OR. polcorr == 5) &
          .OR. (polcorr == 4 .AND. num_iter == 0) ) ) THEN
        npolmod = 2   ! Twice, one vector and one scalar          
     ELSE 
        npolmod = 1   ! Only once either scalar or vector
     ENDIF
     
     ! When polcorr = 5 is selected, only calculate weighting function for iteration 1
     ! iteration 0 if 1st pixel of a x-track position is being retrieved
     IF (polcorr == 5 .AND. nw > 1 .AND. (num_iter > 1 .OR. (num_iter == 0 .AND. currloop /= 0))) THEN
        do_simulation_only = .TRUE.
        do_linearization   = .FALSE.
     ENDIF
    
     DO ipol = 1, npolmod
        IF (ipol == 2) THEN
           NSTOKES = 1 ; NSTREAMS = 4 ! Always scalar for second mode
           
           IF ( lamda < 295.0 .AND. nw > 1 .AND. do_ssfullb295) THEN
              DO_SSFULL = .TRUE.;  DO_SSCORR_TRUNCATION = .FALSE.; DO_DELTAM_SCALING = .FALSE.
           ELSE
              DO_SSFULL = .FALSE.; DO_SSCORR_TRUNCATION = .TRUE.; DO_DELTAM_SCALING = .TRUE.
              !DO_SSCORR_TRUNCATION = .FALSE.; DO_DELTAM_SCALING = .FALSE.
           ENDIF
           IF (.NOT. aerosol) THEN
              DO_SSCORR_TRUNCATION = .FALSE.; DO_DELTAM_SCALING = .FALSE.
           ENDIF
           
           ! Save VECTOR LIDORT results
           prad(iw) = rad(iw, 1)
           rad(iw, 1:nostk) = 0.d0
           IF (do_cfracwf) pcfracwf(iw) = cfracwf(iw, 1)
           IF (do_cfracwf) cfracwf(iw, 1:nostk) =  0.d0
           
           IF ( do_linearization ) THEN
              IF (do_albwf)  palbwf(iw)             = albwf(iw, 1)
              IF (do_fozwf)  pfozwf(iw, :)          = fozwf(iw, :, 1)
              IF (do_codwf)  pfcodwf(iw, nctp:ncbp) = fcodwf(iw, nctp:ncbp, 1)
              IF (do_sprswf) pfsprswf(iw, nup2p(nsfc-1)+1:nfsfc-1) = fsprswf(iw, nup2p(nsfc-1)+1:nfsfc-1, 1)
              IF (do_fraywf) pfraywf(iw, :)         = fraywf(iw, :, 1)
              IF (do_ctpwf)  pctpwf(iw)             = ctpwf(iw, 1)
              IF (do_faerwf) pfaerwf(iw, faer_lvl:nfsfc-1)  = faerwf(iw, faer_lvl:nfsfc-1, 1)           
              IF (do_twaewf) pfaerswf(iw, faer_lvl:nfsfc-1) = faerwf(iw, faer_lvl:nfsfc-1, 1)
              
              ! Initalize those variables to zero again          
              IF (do_albwf) albwf(iw, 1:nostk) = 0.d0
              IF (do_fozwf) fozwf(iw, :, 1:nostk) = 0.d0 
              IF (do_codwf) fcodwf(iw, nctp:ncbp, 1:nostk) = 0.d0
              IF (do_sprswf) fsprswf(iw, nup2p(nsfc-1)+1:nfsfc-1, 1:nostk) = 0.d0
              IF (do_fraywf) fraywf(iw, :, 1:nostk) = 0.d0
              IF (do_ctpwf)  ctpwf(iw, 1:nostk)       = 0.d0
              IF (do_faerwf) faerwf(iw, faer_lvl:nfsfc-1, 1:nostk)  = 0.d0
              IF (do_twaewf) faerswf(iw, faer_lvl:nfsfc-1, 1:nostk) = 0.d0

           ENDIF
        ENDIF
        
        !print *, iw, ipol, do_simulation_only, do_linearization, do_atmos_linearization, do_surface_linearization
        
        radclrcld = 0.0
        DO ic = 1, 2               ! for clear and cloud

           IF (ic == 1) THEN
              do_clouds = .FALSE.; frac = 1.0 - the_cfrac
           ELSE
              do_clouds = .TRUE. ; frac = the_cfrac
           ENDIF
           !WRITE(*, *) ic, frac, do_clouds, do_lambcld, lambcld_refl, do_multi_vza
           IF (frac == 0.0) CYCLE  ! No clear/cloudy part
           
           ! Note convert ozone from DU to molecule/cm^2  here
           IF ((ic == 1) .OR. (.NOT. do_lambcld))  THEN
              lambertian_albedo = albs(iw)
              
              ! Reset up number of layers since Lambertian 
              ! cloudy scene got different layers
              nlayers = nfsfc - 1; nz1 = nlayers   
              
              ! zero O3 weighting function below surface (This is necessary)           
              profilewf(1:n_totalatmos_wfs, nz1+1:nz, 1, 1, 1:NSTOKES, 1) = 0.0
              
              IF (ipol == 1) THEN
                 !WRITE(91, '(4F10.4,3I5)') sza, vza, aza, lambertian_albedo, 1, nlayers, ngreek_moments_input
                 !WRITE(91, '(40D14.6)') height_grid(0:nlayers)
                 !WRITE(91, '(2D14.6)') lamda, depol(iw)
                 
                 CALL LIDORT_PROF_PREP(lamda, raycof(iw), depol(iw), fzs(0:nz1), frhos(1:nz1), &
                      varyprof(1:nz1), nfgas, gasin(1:nfgas), allcrs(iw, 1:nfgas, 1:nz1), allcol(1:nfgas, 1:nz1), &
                      alleta(1:nfgas, 1:nz1), useasy, nmom, aerosol, aersca(1:nz1),      &
                      aerext(1:nz1), aerasy(1:nz1), aermoms(0:nmom, 1:maxgksec, 1:nz1), aermsk(1:nz1), &
                      do_clouds, cldsca(1:nz1), cldext(1:nz1), cldasy(1:nz1), &
                      cldmoms(0:nmom, 1:maxgksec, 1:nz1), cldmsk(1:nz1), problems, &
                      deltau(iw, 1:nz1), delsca(iw, 1:nz1), delo3abs(iw, 1:nz1), delray(iw, 1:nz1))
                 IF (problems) then
                    WRITE(*, *) modulename, ' : Problems encountered in lidort preparation!!!'
                    errstat = pge_errstat_error; RETURN
                 END IF
              ENDIF
           ELSE                 ! lambertian clouds
              nlayers = nctp-1  ! from cloud top to TOA
              nz1 = nlayers;    do_clouds = .FALSE.
              
              IF (the_cfrac == 1.0 .AND. nw /= 1) THEN
                 lambertian_albedo = albs(iw); lambcld_refl = lambertian_albedo
              ELSE
                 lambertian_albedo = lambcld_refl ! use 80% (could be adjusted when the_cfrac gt 0.90)
              ENDIF
              
              ! zero O3 weighting function (This is necessary)           
              profilewf(1:n_totalatmos_wfs, nz1+1:nz, 1, 1, 1:NSTOKES, 1) = 0.0
              
              IF (frac == 1.0 .AND. ipol == 1) THEN
                 CALL LIDORT_PROF_PREP(lamda, raycof(iw), depol(iw), fzs(0:nz1), frhos(1:nz1), &
                      varyprof(1:nz1), nfgas, gasin(1:nfgas), allcrs(iw, 1:nfgas, 1:nz1), allcol(1:nfgas, 1:nz1), &
                      alleta(1:nfgas, 1:nz1), useasy, nmom, aerosol, aersca(1:nz1),      &
                      aerext(1:nz1), aerasy(1:nz1), aermoms(0:nmom, 1:maxgksec, 1:nz1), aermsk(1:nz1), &
                      do_clouds, cldsca(1:nz1), cldext(1:nz1), cldasy(1:nz1), &
                      cldmoms(0:nmom, 1:maxgksec, 1:nz1), cldmsk(1:nz1), problems , &
                      deltau(iw, 1:nz1), delsca(iw, 1:nz1), delo3abs(iw, 1:nz1), delray(iw, 1:nz1))
                 
                 IF (problems) THEN
                    WRITE(*, *) modulename, ' : Problems encountered in lidort preparation!!!'
                    errstat = pge_errstat_error; RETURN
                 END IF
              ENDIF
           ENDIF
           
           ! No Need to reset viewing geometry for each LIDORT call
           !geometry_specheight = height_grid(nlayers)
           !szangles(1) = sza; user_vzangles(1) = vza; user_relazms(1) = aza  

           !xliu: 03/07/2011, switch VLIDORT version (from vv2p4 to vv2p4RTC)
           ! Used in v2p4
           !CALL VLIDORT_L_MASTER_LAMBERTIAN (status_inputcheck, status_calculation) 
           !CALL VLIDORT_STATUS ( status_inputcheck, status_calculation )

           ! used in v2p4RTC
           CALL VLIDORT_L_MASTER ( status_inputcheck, ncheckmessages, checkmessages, &
                checkactions, status_calculation, message, trace_1, trace_2, trace_3)
           CALL VLIDORT_STATUS ( 'o3prof_lidort_error', vlidort_errunit, openfileflag, &
                status_inputcheck,  ncheckmessages, checkmessages, checkactions, &
                status_calculation, message, trace_1, trace_2, trace_3 )

           IF (status_inputcheck   /= vlidort_success) THEN
              WRITE(*, *) modulename, ' : Problems encountered in lidort input check!!!'
              errstat = pge_errstat_error; RETURN
           ENDIF
           IF (status_calculation  /= vlidort_success) THEN
              WRITE(*, *) modulename, ' : Problems encountered in lidort calculation!!!'
              errstat = pge_errstat_error; RETURN
           ENDIF
           
           ! Pixel-independent approximation
           radclrcld(ic, 1:nostk) = stokes(1, 1, 1:nostk, 1) * polerr(iw, ic, 1:nostk)
           rad(iw, 1:nostk)       = rad(iw, 1:nostk) + radclrcld(ic, 1:nostk) * frac
           
           IF (do_linearization ) THEN
              
              ! weighting function per Dobson Unit
              IF ( do_ozwf ) THEN
                 DO istk = 1, nostk
                    fozwf(iw, 1:nfsfc-1, istk) = fozwf(iw, 1:nfsfc-1, istk) + profilewf(ozwfidx, 1:nfsfc-1, 1, 1, istk, 1) &
                         / ozs(1:nfsfc-1)  * polerr(iw, ic, istk) * frac 
                 ENDDO
              ENDIF
              
              ! Weighting function with respect to aerosol/cloud optical depth at the last wavelength
              ! so as to keep the scaling since we are fitting the aod at that wavelength
              IF ( do_taodwf .OR. do_saodwf ) THEN
                 !print *, ipol, faer_lvl, nfsfc-1, frac
                 DO istk = 1, nostk
                    faerwf(iw, faer_lvl:nfsfc-1, istk)  = faerwf(iw, faer_lvl:nfsfc-1, istk) + &
                         profilewf(aodwfidx, faer_lvl:nfsfc-1, 1, 1, istk, 1) / &
                         gaext(actawin, faer_lvl:nfsfc-1) * polerr(iw, ic, istk) *  frac
                    !print *, faerwf(iw, faer_lvl:nfsfc-1, 1)
                 ENDDO
              ENDIF
              IF ( do_twaewf ) THEN
                 DO istk = 1, nostk
                    faerswf(iw, faer_lvl:nfsfc-1, istk)  = faerswf(iw, faer_lvl:nfsfc-1, istk) + &
                         profilewf(twaewfidx, faer_lvl:nfsfc-1, 1, 1, istk, 1) / &
                         gasca(actawin, faer_lvl:nfsfc-1) * gaext(actawin, faer_lvl:nfsfc-1) * &
                         polerr(iw, ic, istk) *  frac
                 ENDDO
              ENDIF
              IF ( do_codwf ) THEN
                 DO istk = 1, nostk
                    fcodwf(iw, nctp:ncbp, istk)  = fcodwf(iw, nctp:ncbp, istk) + &
                         profilewf(codwfidx, nctp:ncbp, 1, 1, istk, 1) / &
                         cldext0(nctp:ncbp) * polerr(iw, ic, istk) *  frac
                 ENDDO
              ENDIF
              IF ( do_sprswf ) THEN
                 DO istk = 1, nostk 
                    fsprswf(iw, nup2p(nsfc-1)+1:nfsfc-1, istk)  = fsprswf(iw, nup2p(nsfc-1)+1:nfsfc-1, istk) + &
                         profilewf(sprswfidx, nup2p(nsfc-1)+1:nfsfc-1, 1, 1, istk, 1) / &
                         (fps(nup2p(nsfc-1)+1:nfsfc-1) - fps(nup2p(nsfc-1):nfsfc-2)) * polerr(iw, ic, istk) *  frac
                 ENDDO
              ENDIF
              
              IF ( do_fraywf ) THEN
                 DO istk = 1, nostk
                    fraywf(iw, 1:nfsfc-1, istk) = fraywf(iw, 1:nfsfc-1, istk) + profilewf(raywfidx, 1:nfsfc-1, 1, 1, istk, 1) &
                         / delray(iw, 1:nfsfc-1)  * polerr(iw, ic, istk) * frac 
                 ENDDO
              ENDIF
              
              ! Non-lambertian clouds: albedo wf is from both clear and cloud
              ! Lamberitan clouds:     if the_cfrac  < 1, albwf from clear only
              !                        if the_cfrac == 1, albwf from cloud only

              
              ! xliu, 03/08/11: the surface albedo weighting functions from v2p4RTC are un-normalized.
              ! Also, we need to initialize n_surface_wfs = 1 for calculating lambertian surface weighting function
              IF (do_albwf) THEN 
                 IF (.NOT. do_lambcld) THEN
                    albwf(iw, 1:nostk) = albwf(iw, 1:nostk) + surfacewf(1, 1, 1, 1:nostk, 1) &
                         * polerr(iw, ic, 1:nostk) * frac !/ lambertian_albedo
                 ELSE IF (the_cfrac < 1.0) THEN
                    IF (ic == 1) albwf(iw, 1:nostk) = albwf(iw, 1:nostk) + &
                            surfacewf(1, 1, 1, 1:nostk, 1) * polerr(iw, ic, 1:nostk) * frac ! / lambertian_albedo 
                 ELSE IF (the_cfrac == 1.0) THEN
                    IF (ic == 2) albwf(iw, 1:nostk) = albwf(iw, 1:nostk) + &
                         surfacewf(1, 1, 1, 1:nostk, 1) * polerr(iw, ic, 1:nostk) * frac !/ lambertian_albedo
                 ENDIF
              ENDIF
              
           ENDIF
           
        ENDDO ! end clear/cloudy scene loop
        radclr(iw, 1:nostk) = radclrcld(1, 1:nostk); radcld(iw, 1:nostk) = radclrcld(2, 1:nostk)
        IF (do_cfracwf) cfracwf(iw, 1:nostk) = radcld(iw, 1:nostk) - radclr(iw, 1:nostk)
     ENDDO    ! end scalar/vector modes
  ENDDO       ! end wavelength loop 
  
  ! On-line polarization correction
  IF ( the_cfrac == 1.0 .AND. do_lambcld) THEN
     nz1 = nctp - 1
  ELSE
     nz1 = nfsfc - 1
  ENDIF

  IF ( (polcorr >= 3 .AND. polcorr <= 5 ) .AND. nw > 1) THEN
     IF (use_effcrs) THEN
        nw0 = nw
     ELSE
        nw0 = ncalcp
     ENDIF
     
     nsprs = nup2p(nsfc-1)+1
     delabs(1:nw0, 1:nz1) = deltau(1:nw0, 1:nz1) - delsca(1:nw0, 1:nz1)
     !xliu, 11/02/2011, add abscrs in the variables
     CALL polcorr_online(num_iter, polcorr, nw0, nz1, nctp, ncbp, nsprs, faer_lvl,        &
          npolcorr, polcorr_idxs(1:npolcorr), do_fozwf, do_albwf, do_faerwf,            &
          do_twaewf, do_codwf, do_sprswf, do_fraywf, do_cfracwf, waves(1:nw0), rad(1:nw0, 1), &
          prad(1:nw0), delabs(1:nw0, 1:nz1), abscrs(1:nw0, 1:nz1), ozs(1:nz1), albwf(1:nw0, 1), palbwf(1:nw0), &
          fozwf(1:nw0, 1:nz1, 1), pfozwf(1:nw0, 1:nz1), faerwf(1:nw0, 1:nz1, 1),         &
          pfaerwf(1:nw0, 1:nz1), faerswf(1:nw0, 1:nz1, 1), pfaerswf(1:nw0, 1:nz1),       &
          fcodwf(1:nw0, 1:nz1, 1), pfcodwf(1:nw0, 1:nz1), fsprswf(1:nw0, 1:nz1, 1),      &
          pfsprswf(1:nw0, 1:nz1), fraywf(1:nw0, 1:nz1, 1), pfraywf(1:nw0, 1:nz1),        &
          cfracwf(1:nw0, 1), pcfracwf(1:nw0))
  ENDIF
  
  ! Radiance Interpolation
  
  IF (nw > 1 .AND. do_radinter ) THEN
     nsprs = nup2p(nsfc-1)+1
    
     CALL radwf_interpol(nw, nz1, nctp, ncbp, nsprs, faer_lvl, do_radcals(1:nw),             &
          do_fozwf, do_albwf, do_faerwf, do_twaewf, do_codwf, do_sprswf, do_cfracwf,         &
          waves, abscrs(1:nw, 1:nz1), ozs(1:nz1), rad(1:nw, 1), fozwf(1:nw, 1:nz1, 1),       &
          albwf(1:nw, 1), cfracwf(1:nw, 1), faerwf(1:nw, 1:nz1, 1), faerswf(1:nw, 1:nz1, 1), &
          fcodwf(1:nw, 1:nz1, 1), fsprswf(1:nw, 1:nz1, 1), errstat)
    
     IF (errstat == pge_errstat_error) RETURN
  ENDIF
  
  IF (nw > 1 .AND. .NOT. use_effcrs) THEN
     IF (do_fraywf .AND. do_fozwf) THEN 
        nsprs = nup2p(nsfc-1)+1
          
        CALL hres_radwf_inter_convol(nw, nz1, nctp, ncbp, nsprs, faer_lvl, do_albwf, &
             do_faerwf, do_twaewf, do_codwf, do_sprswf, do_cfracwf, do_tracewf, use_so2dtcrs, &
             do_o3shi, do_tmpwf, waves, ozs(1:nz1), rad(1:nw, 1), fozwf(1:nw, 1:nz1, 1), &
             albwf(1:nw, 1), cfracwf(1:nw, 1), faerwf(1:nw, 1:nz1, 1), faerswf(1:nw, 1:nz1, 1), &
             fcodwf(1:nw, 1:nz1, 1), fsprswf(1:nw, 1:nz1, 1), fraywf(1:nw, 1:nz1, 1), &
             dads(1:nw, 1:nz1), dadt(1:nw, 1:nz1), abscrs(1:nw, 1:nz), so2crs(1:nw, 1:nz), errstat)

        IF (errstat == pge_errstat_error) RETURN
     ELSE
        WRITE(*, *) 'Must have O3/rayleigh weighting function to do radiance interpolation!!!'
        errstat = pge_errstat_error; RETURN
     ENDIF
  ENDIF
  
  
  !WRITE(71, *) nw, nflay
  !DO iw = 1, nw
  !   IF (do_radcals(iw) ) THEN
  !      write(71, *) '1 '
  !   ELSE
  !      write(71, *) '0 '
  !   ENDIF
  !ENDDO
  !    
  !DO iw = 1, nw 
  !   WRITE(71, '(1000D16.7)') waves(iw), rad(iw, 1), fozwf(iw, 1:nflay, 1), albwf(iw, 1), &
  !        deltau(iw, 1:nflay), delsca(iw, 1:nflay), delo3abs(iw, 1:nflay)
  !ENDDO
  !WRITE(71, '(10D16.7)') ozs(1:nflay)
  !print *, fozwf(100, 1:nflay, 1)
  !
  !STOP

  IF (nw > 1 .AND. .NOT. use_effcrs) THEN
     nw0 = n_rad_wvl
  ELSE
     nw0 = nw
  ENDIF

  ! Calculate desired weighting functions at the end after applying all the correction
  IF ( do_ozwf ) THEN  
     DO i = 1, nl
        fidx = nup2p(i - 1) + 1; lidx = nup2p(i)         
        DO iw = 1, nw0
           DO istk = 1, nostk
              ozwf(iw, i, istk) = SUM(fozwf(iw, fidx:lidx, istk) * ozs(fidx:lidx)) / SUM(ozs(fidx:lidx))
           ENDDO
        ENDDO
     ENDDO
  ENDIF
  
  IF ( do_tmpwf ) THEN  
     DO i = 1, nl
        fidx = nup2p(i - 1) + 1; lidx = nup2p(i)         
        DO iw = 1, nw0
           DO istk = 1, nostk
              tmpwf(iw, i, istk) = SUM(fozwf(iw, fidx:lidx, istk) * ozs(fidx:lidx) &
                   * dadt(iw, fidx:lidx) * (fzs(fidx-1:lidx-1)-fzs(fidx:lidx))) / (fzs(fidx-1) - fzs(lidx))
           ENDDO
        ENDDO
     ENDDO
  ENDIF
  
  IF (do_o3shi) THEN
     DO iw = 1, nw0
        DO istk = 1, nostk
           o3shiwf(iw, istk) = SUM(fozwf(iw, 1:nz1, istk) * ozs(1:nz1) * dads(iw, 1:nz1))
        ENDDO
     ENDDO
  ENDIF

  IF (do_taodwf) THEN
     DO iw = 1, nw0
        DO istk = 1, nostk
           taodwf(iw, istk) = SUM( faerwf(iw, nup2p(ntp)+1:nz1, istk) * &
                gaext(actawin, nup2p(ntp)+1:nz1) ) / SUM (gaext(actawin, nup2p(ntp)+1:nz1))
        ENDDO
     ENDDO
  ENDIF

  IF (do_twaewf) THEN
     DO iw = 1, nw0
        DO istk = 1, nostk
           twaewf(iw, istk) = SUM( faerswf(iw, nup2p(ntp)+1:nz1, istk) * &
                gasca(actawin, nup2p(ntp)+1:nz1) / gaext(actawin, nup2p(ntp)+1:nz1) ) &
                / SUM (gasca(actawin, nup2p(ntp)+1:nz1) / gaext(actawin, nup2p(ntp)+1:nz1))
        ENDDO
     ENDDO
  ENDIF

  IF (do_saodwf) THEN
     DO iw = 1, nw0
        DO istk = 1, nostk
           saodwf(iw, istk) = SUM( faerwf(iw, 1:nup2p(ntp), istk) * &
                gaext(actawin, 1:nup2p(ntp)) ) / SUM (gaext(actawin, 1:nup2p(ntp)))
        ENDDO
     ENDDO
  ENDIF

  IF (do_codwf) THEN
     DO iw = 1, nw0
        DO istk = 1, nostk
           codwf(iw, istk) = SUM( fcodwf(iw, nctp:ncbp, istk) * &
                cldext0(nctp:ncbp) ) / SUM (cldext0(nctp:ncbp))
        ENDDO
     ENDDO
  ENDIF

  IF (do_sprswf) THEN
     DO iw = 1, nw0
        DO istk = 1, nostk
           sprswf(iw, istk) = SUM( fsprswf(iw, nup2p(nsfc-1)+1:nz1, istk) * &
                (fps(nup2p(nsfc-1)+1:nz1) - fps(nup2p(nsfc-1):nz1-1))) / (fps(nz1)-fps(nup2p(nsfc-1)))
        ENDDO
     ENDDO
  ENDIF
   
  !IF (nw > 1) THEN
  !   DO iw = 1, nw
  !      WRITE(78, '(4D16.7)') waves(iw), albwf(iw, 1), sprswf(iw, 1), cfracwf(iw, 1)
  !   ENDDO
  !ENDIF
  
  IF (nw > 1 .AND. do_ozwf .AND. do_tracewf ) &
       CALL GET_TRACEGAS_WF (fozwf(1:nw0, 1:nz, 1), abscrs(1:nw0, 1:nz), so2crs(1:nw0, 1:nz), use_so2dtcrs, &
       rad(1:nw0, 1), nw0, nz, nz1, ozs(1:nz), waves(1:nw0), do_so2zwf, so2zwf(1:nw0, 1))
  
  IF (nw > 1 .AND. do_simu .AND. .NOT. radcalwrt) THEN
     WRITE(78, *) 2, nz1
     WRITE(78, *) 'Profile: Z, P, T, O3, SO2'
     DO i = 1, nz1
        WRITE(78, '(F10.4, 4D16.7)') fzs(i), fps(i), fts(i), fozs(i), allcol(3, i)/refspec_norm(9)
     ENDDO
     WRITE(78, '(A)') 'Wavelength, radiance, albedo wf, ozone wf, aerosol wf, ozcrs, so2crs'
     DO i = 1, nw !160, 369, 209
        WRITE(78, '(F10.4, 500D16.7)') waves(i), rad(i, 1), albwf(i, 1), fozwf(i, 1:nz1, 1), & !, faerwf(i, 1:nz1, 1), &
             allcrs(i, 1, 1:nz1), allcrs(i, 3, nz1)*refspec_norm(9)
     ENDDO
     
     STOP
  ENDIF
 
  RETURN
END SUBROUTINE LIDORT_PROF_ENV

