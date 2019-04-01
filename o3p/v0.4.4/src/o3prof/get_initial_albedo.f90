  ! ***********************************************************
  ! Author:  xiong liu
  ! Date  :  July 24, 2003
  ! Purpose: compute the albedo at 370 nm using a TOMS look-up
  !          table from the convolved reflectance at 370 nm 
  !          triangular slit function with fwhm = 1.132 nm
  !          Assume terrain pressure is fixed at 1.0 atm
  ! ***********************************************************

  ! The albedo from this routine need not be very accurate and 
  ! the solar and radiance spectrum here are not well calibrated, 
  ! think about this more

SUBROUTINE get_initial_albedo (noalb, albedo, pge_error_status)
  
  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY : the_sza_atm, the_vza_atm, &
                                     the_aza_atm, reduce_resolution, use_redfixwav
  USE ozprof_data_module,     ONLY : pos_alb, toms_fwhm, rad_posr, rad_specr, &
                                     sun_posr, sun_specr, ps0, measref, nrefl
  USE OMSAO_errstat_module

  IMPLICIT NONE
 
  ! =================
  ! Output variables
  ! =================
  LOGICAL, INTENT (IN)        :: noalb
  REAL (KIND=dp), INTENT(OUT) :: albedo
  INTEGER, INTENT (OUT)       :: pge_error_status
 
  ! ===============
  ! Local variables
  ! ===============
  INTEGER, PARAMETER            :: nw = 45, midin = 23    ! 45 * 0.05 nm ~ 2.25 nm
  REAL (KIND=dp), DIMENSION(nw) :: wave_arr, rad_arr, irrad_arr, weight
  INTEGER                       :: i, errstat, naw
  REAL (KIND=dp)                :: calc_albedo, wav_interval
 
  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=18), PARAMETER    :: modulename = 'get_initial_albedo'
 
  pge_error_status = pge_errstat_ok
 
  wav_interval = toms_fwhm * 2 / (nw - 1)

  IF (.NOT. reduce_resolution) THEN 
    ! check if the provided wavelength range is enough

    !    print * ,nrefl, rad_posr(1:nrefl)

    IF ( (rad_posr(nrefl) <= pos_alb + toms_fwhm) .OR. &
         (rad_posr(1) >= pos_alb - toms_fwhm)) THEN
      WRITE(*, *) modulename, ' : Raidance does not cover: ', pos_alb, ' nm'
      pge_error_status = pge_errstat_error; RETURN
    END IF
    
    IF ((sun_posr(nrefl) <= pos_alb + toms_fwhm) .OR. &
         (sun_posr(1) >= pos_alb - toms_fwhm)) THEN
      WRITE(*, *) modulename, ': Solar does not cover: ', pos_alb, ' nm'
      pge_error_status = pge_errstat_error; RETURN
    END IF

    ! get wavelength positions
    DO i = 1, nw
      wave_arr(i) = pos_alb + REAL(i - midin, KIND=dp) * wav_interval
    END DO
    weight(1:nw) = 1.0 - ABS(wave_arr(1:nw) - pos_alb) / toms_fwhm
    naw = nw
  ELSE
    ! check if the provided wavelength range is enough
    IF ( rad_posr(nrefl) < pos_alb .OR. rad_posr(1) > pos_alb) THEN
      WRITE(*, *) modulename, ' : Raidance does not cover: ', pos_alb, ' nm'
      pge_error_status = pge_errstat_error; RETURN
    END IF
    
    IF (sun_posr(nrefl) < pos_alb .OR. sun_posr(1) > pos_alb ) THEN
      WRITE(*, *) modulename, ': Solar does not cover: ', pos_alb, ' nm'
      pge_error_status = pge_errstat_error; RETURN
    END IF
    
    naw = 1; wave_arr(1) = pos_alb; weight(1) = 1.0
  
  ENDIF
  ! wave_r ~ 346nm
   
  IF ( use_redfixwav .AND. nrefl == 1) THEN
    measref = rad_specr(1) / sun_specr(1)
  ELSE
    ! interpolate solar and radiance spectra to get spectra at the above positions
    CALL BSPLINE(rad_posr(1:nrefl), rad_specr(1:nrefl), nrefl, &
         wave_arr(1:naw), rad_arr(1:naw), naw, errstat)
    IF (errstat < 0) THEN
      WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat; STOP
    ENDIF
    
    CALL BSPLINE(sun_posr(1:nrefl), sun_specr(1:nrefl), nrefl,  wave_arr(1:naw), &
         irrad_arr(1:naw), naw, errstat)
    IF (errstat < 0) THEN
      WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat; STOP
    ENDIF
    
    ! compute the reflectance and calculate albedo
    measref = SUM(rad_arr(1:naw) / irrad_arr(1:naw) * weight(1:naw)) / SUM (weight(1:naw)) 
  ENDIF
  
  ! noalb is TRUE' 
  
  IF (noalb) RETURN
  
  albedo = calc_albedo(measref, ps0, the_sza_atm, the_vza_atm, the_aza_atm)
  ! WRITE(*, '(A, d12.4)') ' The initial albedo is: ', albedo

  IF (albedo <= -0.1 .OR. albedo >= 1.2) THEN
    WRITE(*, *) modulename, ' : Surface albedo out of bounds!!!'
    pge_error_status = pge_errstat_error; RETURN
  END IF
  
  RETURN
  
END SUBROUTINE get_initial_albedo


! ==============================================================
!  Calculate surface albedo for a certain spres, sza, vza, aza 
!  using a reflectance look-up table calculated using TOMRAD
! ==============================================================

FUNCTION calc_albedo(refl, spres, sza, vza, aza) RESULT (albedo)
  
  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY : deg2rad
  USE OMSAO_variables_module, ONLY  : refdbdir
  USE ozprof_data_module,      ONLY : alb_tbl_fname, atmos_unit
  IMPLICIT NONE
  
  ! =================
  ! Input variables
  ! =================
  REAL (KIND=dp), INTENT(IN)       :: refl, sza, vza, aza, spres

  ! =================
  ! Result variables
  ! =================
  REAL (KIND=dp)                   :: albedo
  
  ! =================
  ! local variables
  ! =================
  INTEGER, PARAMETER              :: mp = 12, msza = 13, mvza = 7   ! look-up table dimension
  INTEGER                         :: i, j, k, pin, szain, vzain
  INTEGER, SAVE                   ::  nsza, nvza, np
  REAL (KIND=dp), SAVE, DIMENSION(msza) :: sza_arr
  REAL (KIND=dp), SAVE, DIMENSION(mvza) :: vza_arr
  REAL (KIND=dp), SAVE, DIMENSION(mp)   :: parr, sb
  REAL (KIND=dp), SAVE, DIMENSION(mp, msza, mvza) :: rad0, rad1, rad2, t  
  REAL (KIND=dp) :: q1, q2, the_i0, the_i1, the_i2, the_t, isurf, sza1, vza1,&
       aza1, spres1, pfrac, sfrac, vfrac, the_sb
  LOGICAL, SAVE  :: first = .TRUE.
  
  IF (first) THEN
     alb_tbl_fname = TRIM(ADJUSTL(refdbdir)) // 'toms_370nm_refl_new.dat'
     OPEN (UNIT=atmos_unit, file = alb_tbl_fname, status = 'old')
     READ (atmos_unit, *) np, nsza, nvza
     READ (atmos_unit, *) parr(1:np)
     READ (atmos_unit, *) sza_arr(1:nsza)
     READ (atmos_unit, *) vza_arr(1:nvza)
     READ (atmos_unit, *) (((rad0(k, i, j), j=1, nvza), i=1, nsza), k=1, np)
     READ (atmos_unit, *) (((rad1(k, i, j), j=1, nvza), i=1, nsza), k=1, np)
     READ (atmos_unit, *) (((rad2(k, i, j), j=1, nvza), i=1, nsza), k=1, np)
     READ (atmos_unit, *) (((t(k, i, j), j=1, nvza), i=1, nsza), k=1, np)
     READ (atmos_unit, *) sb
     CLOSE (UNIT=atmos_unit)
     first = .FALSE.
  ENDIF

  spres1 = spres / 1013.25;

  ! find index for spres, sza, vza
  pin = MINVAL(MINLOC(parr, MASK =(parr >= spres1)))
  szain = MINVAL(MAXLOC(sza_arr, MASK = (sza_arr <= sza)))
  vzain = MINVAL(MAXLOC(vza_arr, MASK = (vza_arr <= vza)))

  IF (pin > np) pin = np - 1
  IF (szain > nsza) szain = szain - 1
  IF (vzain > nvza) vzain = vzain - 1

  pfrac = (spres1 - parr(pin)) / (parr(pin+1)-parr(pin))
  sfrac = (COS(sza*deg2rad) - COS(sza_arr(szain)*deg2rad)) / &
       (COS(sza_arr(szain+1)*deg2rad) - COS(sza_arr(szain)*deg2rad))
  vfrac = (COS(vza*deg2rad) - COS(vza_arr(vzain)*deg2rad)) / &
       (COS(vza_arr(vzain+1)*deg2rad) - COS(vza_arr(vzain)*deg2rad))

  ! use 3-D inteprolation (pressure, cos(sza), vza(sza))
  CALL BLEND_103(pfrac, sfrac, vfrac, rad0(pin, szain, vzain), &
       rad0(pin, szain, vzain+1), rad0(pin, szain+1, vzain),   &
       rad0(pin, szain+1, vzain+1), rad0(pin+1, szain, vzain), &
       rad0(pin+1, szain, vzain+1), rad0(pin+1, szain+1, vzain),&
       rad0(pin+1,szain+1, vzain+1), the_i0)
  
  CALL BLEND_103(pfrac, sfrac, vfrac, rad1(pin, szain, vzain), &
       rad1(pin, szain, vzain+1), rad1(pin, szain+1, vzain),   &
       rad1(pin, szain+1, vzain+1), rad1(pin+1, szain, vzain), &
       rad1(pin+1, szain, vzain+1), rad1(pin+1, szain+1, vzain),&
       rad1(pin+1,szain+1, vzain+1), the_i1)  

  CALL BLEND_103(pfrac, sfrac, vfrac, rad2(pin, szain, vzain), &
       rad2(pin, szain, vzain+1), rad2(pin, szain+1, vzain),   &
       rad2(pin, szain+1, vzain+1), rad2(pin+1, szain, vzain), &
       rad2(pin+1, szain, vzain+1), rad2(pin+1, szain+1, vzain), &
       rad2(pin+1,szain+1, vzain+1), the_i2)  
  
  CALL BLEND_103(pfrac, sfrac, vfrac, t(pin, szain, vzain), &
       t(pin, szain, vzain+1), t(pin, szain+1, vzain),      &
       t(pin, szain+1, vzain+1), t(pin+1, szain, vzain),    &
       t(pin+1, szain, vzain+1), t(pin+1, szain+1, vzain),  &
       t(pin+1,szain+1, vzain+1), the_t) 

  ! use 1-D interpolation for sb
  CALL BLEND_101(pfrac, sb(pin), sb(pin+1), the_sb)

  ! Calculate surface albedo
  sza1 = sza * deg2rad; vza1 = vza * deg2rad; aza1 = aza * deg2rad
  
  q1 = -0.375D0 * COS(sza1) *  SIN(sza1) * SIN(vza1)
  q2 = 0.09375D0 * (SIN(sza1)**2) * (SIN(vza1)**2) / COS(vza1)
  the_i1 = the_i1 * q1 * COS(aza1)
  the_i2 = the_i2 * q2 * COS(2.0 * aza1)
  
  isurf = refl - the_i0 - the_i1 - the_i2
  albedo = isurf / (the_t + isurf * the_sb)
  
  RETURN
END FUNCTION calc_albedo

! =============================================================
! Obtain GOME MLER surface albedo developed by Kolemeijer/TOMS
! =============================================================
SUBROUTINE get_gome_alb(month, elons, elats, region, albarr, nalb)
  
  USE OMSAO_precision_module 
  USE OMSAO_variables_module, ONLY: atmdbdir
  USE ozprof_data_module,     ONLY: atmos_unit
  IMPLICIT NONE
  
  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)         :: month, region, nalb
  REAL (KIND=dp), DIMENSION(2), INTENT(IN)  :: elons, elats  
  REAL (KIND=dp), INTENT(OUT), DIMENSION(nalb) :: albarr

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER :: numy = 8, nlat=180, nlon=360
  REAL (KIND=dp), PARAMETER :: longrid = 1.0, latgrid = 1.0
  CHARACTER (LEN=2), DIMENSION(12),   PARAMETER :: monstr = (/'01', '02', &
       '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'/)
  CHARACTER (LEN=3), DIMENSION(4), PARAMETER    :: vistr = (/'495', '555',&
       '610', '670'/)
  CHARACTER (LEN=130)            :: alb_fname
  CHARACTER (LEN=14)             :: tmpstr
  INTEGER, SAVE, DIMENSION(5, nlon, nlat) :: glbalb
  INTEGER                        :: i, j, k, ialb, latin, lonin,  npix, nact
  LOGICAL                        :: file_exist
  LOGICAL, SAVE                  :: first = .TRUE.
  REAL (KIND=dp)                 :: sumalb, lon, lat

  ! read albedo at 380 nm  (closest to 370 nm)
  IF (first) THEN
     IF (region /= 2) THEN
        alb_fname = TRIM(ADJUSTL(atmdbdir)) // 'albdb/alb' // monstr(month) // '380.dat'
        
        ! Determine if file exists or not
        INQUIRE (FILE= alb_fname, EXIST= file_exist)
        IF (.NOT. file_exist) alb_fname = TRIM(ADJUSTL(atmdbdir)) // 'albdb/alb' // '380.dat'
        
        OPEN (UNIT = atmos_unit, file=alb_fname , status = 'unknown')     
        DO i = 1, 3
           READ (atmos_unit, '(A)')
        END DO
        DO i = 1, nlat 
           READ (atmos_unit,'(14(25i3/),10i3,a14)') (glbalb(1, j, i), j=1, nlon), tmpstr
        ENDDO
        CLOSE (atmos_unit)
     ENDIF

     IF (region /= 1) THEN
        DO k = 1, 4
           alb_fname = TRIM(ADJUSTL(atmdbdir)) // 'albdb/alb' // monstr(month) // vistr(k) // '.dat'
           
           ! Determine if file exists or not
           INQUIRE (FILE= alb_fname, EXIST= file_exist)
           IF (.NOT. file_exist) alb_fname = TRIM(ADJUSTL(atmdbdir)) // 'albdb/alb' // vistr(k) // '.dat'
           OPEN (UNIT = atmos_unit, file=alb_fname , status = 'unknown')     
           DO i = 1, 3
              READ (atmos_unit, '(A)')
           END DO
           DO i = 1, nlat 
              READ (atmos_unit,'(14(25i3/),10i3,a14)') (glbalb(k+1, j, i), j=1, nlon), tmpstr
           ENDDO
           CLOSE (atmos_unit)
        ENDDO
     ENDIF
           
     first = .FALSE.
  ENDIF

  ialb = 0
  npix = NINT((elons(2) - elons(1)) / longrid)
  IF (npix == 0) npix = 1  

  IF (region /= 2) THEN
     sumalb = 0.0D0; nact=0
     ialb = ialb + 1
     DO i = 1, npix
        IF (npix > 1) THEN
           lon = elons(1)  + (i - 1 + 0.5) * longrid
           lat = elats(1)  + (lon - elons(1)) / (elons(2)-elons(1)) * (elats(2) - elats(1))
        ELSE
           lon = (elons(1) + elons(2) ) / 2.0
           lat = (elats(1) + elats(2) ) / 2.0
        ENDIF
        
        lonin = INT((lon + 180.0) / longrid) + 1
        latin = INT((lat + 90.0)  / latgrid) + 1
        lonin = MOD(lonin, nlon)
        IF (latin > nlat) latin = nlat
        
        IF (glbalb(1, lonin, latin) >= 0.0) THEN
           sumalb = sumalb + glbalb(1, lonin, latin)
           nact = nact + 1
        ENDIF
     ENDDO

     IF (nact > 0) THEN
        albarr(ialb) = sumalb / nact / 1000.0
     ELSE
        albarr(ialb) = 0.10
     ENDIF
  ENDIF

  ! read albedo in visible: 495, 555, 610, 670 nm 
  IF (region /= 1) THEN
     DO k = 1, 4
        ialb = ialb + 1
        sumalb = 0.0D0; nact=0
        
        DO i = 1, npix   
           IF (npix > 1) THEN
              lon = elons(1)  + (i - 1 + 0.5) * longrid
              lat = elats(1)  + (lon - elons(1)) / (elons(2)-elons(1)) * (elats(2) - elats(1))
           ELSE
              lon = (elons(1) + elons(2) ) / 2.0
              lat = (elats(1) + elats(2) ) / 2.0
           ENDIF
           
           lonin = INT((lon + 180.0) / longrid) + 1
           latin = INT((lat + 90.0)  / latgrid) + 1
           lonin = MOD(lonin, nlon)
           IF (latin > nlat) latin = nlat
           
           IF (glbalb(k+1, lonin, latin) >= 0.0) THEN
              sumalb = sumalb + glbalb(k+1, lonin, latin)
              nact = nact + 1
           ENDIF
        ENDDO
        IF (nact > 0) THEN
           albarr(ialb) = sumalb / nact / 1000.0
        ELSE
           albarr(ialb) = 0.10
        ENDIF
     ENDDO
  ENDIF

  RETURN
END SUBROUTINE get_gome_alb

  
SUBROUTINE adj_albcfrac(albedo, cfrac, ctau, errstat)
  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY : sza => the_sza_atm, fitvar_rad_apriori, &
       vza => the_vza_atm, aza => the_aza_atm, fitvar_rad_init, scnwrt
  USE ozprof_data_module,     ONLY : nlay, ozp_fidx => ozprof_start_index, &
       ozp_lidx => ozprof_end_index, t_fidx, t_lidx, pos_alb, measref, &
       do_lambcld, lambcld_refl, the_cfrac, num_iter, the_cbeta, ncbp, &
       nctp, fzs, lambcld_initalb, taodfind, twaefind, scacld_initcod
 USE OMSAO_errstat_module

  IMPLICIT NONE

  ! Modified variables
  INTEGER, INTENT(OUT)            :: errstat
  REAL (KIND=dp), INTENT(INOUT)   :: albedo, cfrac, ctau

  ! Local variables
  INTEGER, PARAMETER              :: ns = 1, nos = 1, nalb=1, nostk=1, nwfc=1
  INTEGER, DIMENSION(nalb)        :: albpmax, albpmin
  INTEGER, DIMENSION(nalb)        :: wfcpmax, wfcpmin
  INTEGER                         :: i
  REAL (KIND=dp), DIMENSION(nlay) :: tprof, ozprof, ozadj, ozaprof
  REAL (KIND=dp), DIMENSION(ns)   :: waves, albwf, cfracwf, simrad, codwf, &
       ctpwf, taodwf, twaewf, saodwf, sprswf, so2zwf, walb0s
  REAL (KIND=dp), DIMENSION(ns)   :: o3shiwf
  REAL (KIND=dp), DIMENSION(ns, nlay) :: ozwf, tmpwf
  REAL (KIND=dp), DIMENSION(nalb) :: albarr
  REAL (KIND=dp), DIMENSION(nalb) :: wfcarr
  REAL (KIND=dp), DIMENSION(nos)  :: o3shi
  REAL (KIND=dp)                  :: delta_cfrac, delta_alb, delta_cod, &
       simrad1, simrad2, spres, newoz, initalb
  LOGICAL                  :: do_albwf, do_tmpwf, do_ozwf, do_o3shi, negval, vary_sfcalb, &
       do_taodwf, do_twaewf, do_saodwf, do_cfracwf, do_ctpwf, do_codwf, do_sprswf, do_so2zwf
  LOGICAL, DIMENSION(nlay) :: ozvary

  do_tmpwf  = .FALSE.; do_ozwf   = .FALSE.; do_o3shi = .FALSE.; ozvary   = .FALSE.
  do_taodwf = .FALSE.; do_saodwf = .FALSE.; do_twaewf = .FALSE.; do_ctpwf = .FALSE.
  do_sprswf = .FALSE.; do_so2zwf = .FALSE.
  IF (cfrac == 1.0) cfrac = 0.95 ! Calculate cloud fraction weighting function

  ! ======= Set up ozone, temperature, albedo, lamda for LIDORT ============
  ozprof(1:nlay) =  fitvar_rad_init(ozp_fidx:ozp_lidx)
  ozaprof(1:nlay) = fitvar_rad_apriori(ozp_fidx:ozp_lidx)
  tprof(1:nlay)  =  fitvar_rad_init(t_fidx:t_lidx)
  !xliu (02/01/2007): adjust ozone profile for negative ozone values
  !                   radiances will be corrected using ozone weighting function
  ozadj(1:nlay) = 0.0; negval = .FALSE.
  DO i = 1, nlay
     IF (ozprof(i) <= 0.0) THEN
        newoz  = MIN(0.5d0, ozaprof(i))
        negval = .TRUE. ; ozadj(i)  = newoz - ozprof(i); ozprof(i) = newoz
        do_ozwf = .TRUE.; ozvary(i) = .TRUE.
     ENDIF
  ENDDO

  o3shi(1) = 0.0  ; waves(1) = pos_alb
  albpmax(1) = 1  ; albpmin(1)=1 
  wfcpmax(1) = 1  ; wfcpmin(1)=1 

  num_iter = 0
  IF (scnwrt) WRITE(*, '(A20,4d14.5)') '  Rs, Rc, Fc, Tc: ', albedo, lambcld_refl, cfrac, ctau

  initalb = albedo
  DO 
     ! Weighting function needs to be calculated
     do_albwf = .FALSE.; do_cfracwf = .FALSE.; do_codwf = .FALSE.
     IF (cfrac <= 0.0) THEN
        do_albwf  = .TRUE. 
     ELSE IF (cfrac < 1.0) THEN
        do_cfracwf = .TRUE.
     ELSE IF (do_lambcld) THEN
        do_albwf = .TRUE.
     ELSE
        do_codwf = .TRUE.
     ENDIF

     albarr(1) = albedo; vary_sfcalb = .FALSE.; walb0s(1) = albedo
     the_cfrac = cfrac;  wfcarr(1) = cfrac
     
     CALL LIDORT_PROF_ENV(do_ozwf, do_albwf, do_tmpwf, do_o3shi, ozvary,    &
          do_taodwf, do_twaewf, do_saodwf, do_cfracwf, do_ctpwf, do_codwf,  do_sprswf, do_so2zwf, &
          ns, waves, nos, o3shi, sza, vza, aza, nlay, ozprof, tprof, nalb, albarr, albpmin, &
          albpmax, vary_sfcalb, walb0s, nwfc, wfcarr, wfcpmin, wfcpmax, nostk, albwf, ozwf, tmpwf, o3shiwf, cfracwf, &
          codwf, ctpwf, taodwf, twaewf, saodwf, sprswf, so2zwf, simrad, errstat)


     !xliu (02/01/2007): correct radiances based on ozone weighting function to deal with negative ozone values
     IF (negval) THEN
        DO i = 1, nlay 
           IF (ozadj(i) > 0) THEN
              simrad(1:ns) = simrad(1:ns) - ozadj(i) * ozwf(1:ns, i) 
           ENDIF
        ENDDO
     ENDIF
     
     IF (errstat == pge_errstat_error) RETURN 

     ! Initial delta values to zero
     delta_cfrac = 0.0; delta_alb = 0.0; delta_cod = 0.0
     IF (cfrac > 0.0 .AND. cfrac < 1.0D0) THEN  
        delta_cfrac = (measref - simrad(1)) / cfracwf(1)
        cfrac = cfrac + delta_cfrac
        IF (cfrac > 0.0 .AND. cfrac < 1.0D0) THEN
           IF (scnwrt) WRITE(*, '(A20,4d14.5)') '  Rs, Rc, Fc, Tc: ', albedo, lambcld_refl, cfrac, ctau
           EXIT
        ENDIF
     ELSE IF (cfrac >= 1.0D0 .AND. do_lambcld) THEN   ! Adjust the lambertian cloud albedo
        ! Here albwf is the albedo wf for a fully cloudy conditions
        delta_alb = (measref - simrad(1)) / albwf(1)
        lambcld_refl = lambcld_refl + delta_alb
        
        IF (lambcld_refl < lambcld_initalb) THEN  ! cfracwf < 0, surface contribution comparable to cloud
           lambcld_refl = lambcld_initalb; cfrac = 0.0
           delta_cfrac = -1.0; delta_alb = 0.0
        ENDIF
     ELSE IF (cfrac >= 1.0D0 .AND. .NOT. do_lambcld) THEN  ! Cloud optical thickness
        delta_cod = (measref - simrad(1)) / codwf(1)
        ctau = ctau + delta_cod

        IF (ctau < scacld_initcod) THEN  ! cfracwf < 0, surface contribution comparable to cloud
           ctau = scacld_initcod; cfrac = 0.0
           delta_cfrac = -1.0; delta_cod = 0.0
        ENDIF
        the_cbeta  = ctau /(fzs(nctp-1) - fzs(ncbp))
     ELSE
        delta_alb = (measref - simrad(1)) / albwf(1)
        albedo = albedo + delta_alb
     ENDIF

     IF (scnwrt) WRITE(*, '(A20,4d14.5)') '  Rs, Rc, Fc, Tc: ', albedo, lambcld_refl, cfrac, ctau
     IF (cfrac < 0.00 .OR. albedo < 0.00 .OR. cfrac > 1.0D0 .OR. albedo > 1.0D0) THEN
        IF (cfrac < 0.0)        THEN
           cfrac = 0.0
        ELSE IF (cfrac > 1.0D0) THEN
           cfrac = 1.0D0

        ! Quality Flags can be added for the following two cases
        ELSE IF (albedo < 0.0)  THEN ! Suggest aerosols are under/over estimated (absorbing/nonabsorbing)
           albedo = 0.001
           EXIT  
        ELSE                         ! Suggest cloud exists (polar regions), but need cloud-top pressure
           albedo = 1.0;     EXIT    
        ENDIF
     ENDIF

     ! Exit if the change in albedo or cloud fraction is smaller than 0.001
     IF ( ABS(delta_cfrac) <= 0.001 .AND. ABS(delta_alb) <= 0.001 .AND. ABS(delta_cod) <= 0.001 ) EXIT
     
     num_iter = num_iter + 1
     IF (num_iter >= 10) EXIT
  ENDDO

  IF (taodfind > 0 .OR. twaefind> 0) THEN
     albedo = initalb
  ENDIF

  RETURN
END SUBROUTINE adj_albcfrac

! ==============================================================
!  Calculate backscattered radiance for a certain spres, sza, 
!  vza, aza, surface albedo using a reflectance look-up table 
!  calculated using TOMRAD
! ==============================================================

SUBROUTINE get_refrad (nm, rs, spres, sza, vza, aza, refl)
  
  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY : deg2rad
  USE OMSAO_variables_module,  ONLY : refdbdir
  USE ozprof_data_module,      ONLY : alb_tbl_fname, atmos_unit
  IMPLICIT NONE
  
  ! =================
  ! Input variables
  ! =================
  INTEGER, INTENT (IN)                      :: nm
  REAL (KIND=dp), INTENT(IN)                :: sza, vza, aza, spres
  REAL (KIND=dp), DIMENSION(nm), INTENT(IN) :: rs

  ! =================
  ! Output variables
  ! =================
  REAL (KIND=dp), DIMENSION(nm), INTENT(OUT) :: refl
  
  ! =================
  ! local variables
  ! =================
  INTEGER, PARAMETER              :: mp = 12, msza = 13, mvza = 7   ! look-up table dimension
  INTEGER                         :: i, j, k, pin, szain, vzain
  INTEGER, SAVE                   ::  nsza, nvza, np
  REAL (KIND=dp), SAVE, DIMENSION(msza) :: sza_arr
  REAL (KIND=dp), SAVE, DIMENSION(mvza) :: vza_arr
  REAL (KIND=dp), SAVE, DIMENSION(mp)   :: parr, sb
  REAL (KIND=dp), SAVE, DIMENSION(mp, msza, mvza) :: rad0, rad1, rad2, t  
  REAL (KIND=dp) :: q1, q2, the_i0, the_i1, the_i2, the_t, sza1, vza1,&
       aza1, pfrac, sfrac, vfrac, the_sb
  LOGICAL, SAVE  :: first = .TRUE.
  
  IF (first) THEN
     alb_tbl_fname = TRIM(ADJUSTL(refdbdir)) // 'toms_370nm_refl_new.dat'
     OPEN (UNIT=atmos_unit, file = alb_tbl_fname, status = 'old')
     READ (atmos_unit, *) np, nsza, nvza
     READ (atmos_unit, *) parr(1:np)
     READ (atmos_unit, *) sza_arr(1:nsza)
     READ (atmos_unit, *) vza_arr(1:nvza)
     READ (atmos_unit, *) (((rad0(k, i, j), j=1, nvza), i=1, nsza), k=1, np)
     READ (atmos_unit, *) (((rad1(k, i, j), j=1, nvza), i=1, nsza), k=1, np)
     READ (atmos_unit, *) (((rad2(k, i, j), j=1, nvza), i=1, nsza), k=1, np)
     READ (atmos_unit, *) (((t(k, i, j), j=1, nvza), i=1, nsza), k=1, np)
     READ (atmos_unit, *) sb
     CLOSE (UNIT=atmos_unit)
     first = .FALSE.
  ENDIF

  ! find index for spres, sza, vza
  pin = MINVAL(MINLOC(parr, MASK =(parr >= spres)))
  szain = MINVAL(MAXLOC(sza_arr, MASK = (sza_arr <= sza)))
  vzain = MINVAL(MAXLOC(vza_arr, MASK = (vza_arr <= vza)))

  IF (pin > np) pin = np - 1
  IF (szain > nsza) szain = szain - 1
  IF (vzain > nvza) vzain = vzain - 1

  pfrac = (spres - parr(pin)) / (parr(pin+1)-parr(pin))
  sfrac = (COS(sza*deg2rad) - COS(sza_arr(szain)*deg2rad)) / &
       (COS(sza_arr(szain+1)*deg2rad) - COS(sza_arr(szain)*deg2rad))
  vfrac = (COS(vza*deg2rad) - COS(vza_arr(vzain)*deg2rad)) / &
       (COS(vza_arr(vzain+1)*deg2rad) - COS(vza_arr(vzain)*deg2rad))

  ! use 3-D inteprolation (pressure, cos(sza), vza(sza))
  CALL BLEND_103(pfrac, sfrac, vfrac, rad0(pin, szain, vzain), &
       rad0(pin, szain, vzain+1), rad0(pin, szain+1, vzain),   &
       rad0(pin, szain+1, vzain+1), rad0(pin+1, szain, vzain), &
       rad0(pin+1, szain, vzain+1), rad0(pin+1, szain+1, vzain),&
       rad0(pin+1,szain+1, vzain+1), the_i0)
  
  CALL BLEND_103(pfrac, sfrac, vfrac, rad1(pin, szain, vzain), &
       rad1(pin, szain, vzain+1), rad1(pin, szain+1, vzain),   &
       rad1(pin, szain+1, vzain+1), rad1(pin+1, szain, vzain), &
       rad1(pin+1, szain, vzain+1), rad1(pin+1, szain+1, vzain),&
       rad1(pin+1,szain+1, vzain+1), the_i1)  

  CALL BLEND_103(pfrac, sfrac, vfrac, rad2(pin, szain, vzain), &
       rad2(pin, szain, vzain+1), rad2(pin, szain+1, vzain),   &
       rad2(pin, szain+1, vzain+1), rad2(pin+1, szain, vzain), &
       rad2(pin+1, szain, vzain+1), rad2(pin+1, szain+1, vzain), &
       rad2(pin+1,szain+1, vzain+1), the_i2)  
  
  CALL BLEND_103(pfrac, sfrac, vfrac, t(pin, szain, vzain), &
       t(pin, szain, vzain+1), t(pin, szain+1, vzain),      &
       t(pin, szain+1, vzain+1), t(pin+1, szain, vzain),    &
       t(pin+1, szain, vzain+1), t(pin+1, szain+1, vzain),  &
       t(pin+1,szain+1, vzain+1), the_t) 

  ! use 1-D interpolation for sb
  CALL BLEND_101(pfrac, sb(pin), sb(pin+1), the_sb)

  ! Calculate surface albedo
  sza1 = sza * deg2rad; vza1 = vza * deg2rad; aza1 = aza * deg2rad
  
  q1 = -0.375D0 * COS(sza1) *  SIN(sza1) * SIN(vza1)
  q2 = 0.09375D0 * (SIN(sza1)**2) * (SIN(vza1)**2) / COS(vza1)
  the_i1 = the_i1 * q1 * COS(aza1)
  the_i2 = the_i2 * q2 * COS(2.0 * aza1)

  refl = the_i0 + the_i1 + the_i2 + rs * the_t / (1.0 - rs * the_sb)

  !isurf = refl - the_i0 - the_i1 - the_i2
  !albedo = isurf / (the_t + isurf * the_sb)
  
  RETURN
END SUBROUTINE get_refrad

! This version calculate radiance from a look-up table, which have treated polarization 
! and performed spherical correction for the outgoing beam
SUBROUTINE adj_albcfrac1(spres, ctp, albedo, cfrac)
  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY : sza => the_sza_atm, &
       vza => the_vza_atm, aza => the_aza_atm, scnwrt
  USE ozprof_data_module,     ONLY : measref, do_lambcld, lambcld_refl, the_cfrac, num_iter

  IMPLICIT NONE

  ! Modified variables
  REAL (KIND=dp), INTENT(INOUT)   :: albedo, cfrac, spres, ctp

  ! Local variables
  INTEGER, PARAMETER  :: nalb = 2
  INTEGER             :: ic
  REAL (KIND=dp)      :: wave, albwf, cfracwf, simrad, delta_cfrac, &
       delta_alb, delta_cod, simrad1, simrad2
  REAL (KIND=dp), DIMENSION(nalb)    :: refl, albs
  REAL (KIND=dp), DIMENSION(2)       :: radcldclr, albwfs

  spres = spres / 1013.25; ctp = ctp / 1013.25
 
  ! Initialize Lambertian cloud albedo to be 80%
  IF (do_lambcld .AND. cfrac > 0) lambcld_refl = 0.80D0

  num_iter = 0
  IF (scnwrt) WRITE(*, '(A20,3d14.5)') '  Rs, Rc, Fc: ', albedo, lambcld_refl, cfrac

  DO 
     IF (cfrac /= 1.0) THEN
        albs(1) = albedo; albs(2) = albedo + 0.01
        CALL get_refrad (nalb, albs, spres, sza, vza, aza, refl)
        albwfs(1) = (refl(2) - refl(1)) / 0.01
        radcldclr(1) = refl(1)
     ELSE
        albwfs(1) = 0.0; radcldclr(1) = 0.0
     ENDIF

     IF (cfrac /= 0.0) THEN
        albs(1) = lambcld_refl; albs(2) = lambcld_refl + 0.01
        CALL get_refrad (nalb, albs, ctp, sza, vza, aza, refl)
        albwfs(2) = (refl(2) - refl(1)) / 0.01
        radcldclr(2) = refl(1)
     ELSE
        albwfs(2) = 0.0; radcldclr(2) = 0.0
     ENDIF

     simrad = radcldclr(1) * (1.0 - cfrac) + radcldclr(2) * cfrac
     albwf = albwfs(1) * (1.0 - cfrac) + albwfs(2) * cfrac
     cfracwf = radcldclr(2) - radcldclr(1)
     
     !print *, measref, simrad
     !print *, albwf, cfracwf

     IF (cfrac > 0.0 .AND. cfrac < 1.0D0) THEN  
        delta_cfrac = (measref - simrad) / cfracwf
        cfrac = cfrac + delta_cfrac
        delta_alb = 1.0
     ELSE IF (cfrac >= 1.0D0) THEN  ! adjust the lambertian cloud albedo
        ! here the albwf is the albedo wf for a complete cloud
        delta_alb = (measref - simrad) / albwf
        lambcld_refl = lambcld_refl + delta_alb
        delta_cfrac = 1.0
     ELSE
        delta_alb = (measref - simrad) / albwf
        albedo = albedo + delta_alb
        delta_cfrac = 1.0
     ENDIF

     IF (scnwrt) WRITE(*, '(A20,3d14.5)') ' Rs, Rc, Fc: ', albedo,  lambcld_refl, cfrac

     IF (cfrac < 0.00 .OR. albedo < 0.00 .OR. cfrac > 1.0D0 .OR. albedo > 1.0D0) THEN
        IF (cfrac < 0.0) THEN
           cfrac = 0.0
        ELSE IF (cfrac > 1.0D0) THEN
           cfrac = 1.0D0
        ELSE IF (albedo < 0.0) THEN
           albedo = 0.00001;  EXIT
        ELSE
           albedo = 1.0; EXIT
        ENDIF
     ENDIF   

     ! Exit if the change in albedo or cloud fraction is smaller than 0.001
     IF ( ABS(delta_cfrac) <= 0.001 .OR. ABS(delta_alb) <= 0.001 ) EXIT
     
     num_iter = num_iter + 1

  ENDDO

  spres = spres * 1013.25; ctp = ctp * 1013.25

  RETURN
END SUBROUTINE adj_albcfrac1
  
! =============================================================
! Obtain TOMS MLER surface albedo
! =============================================================
SUBROUTINE get_toms_alb(month, elons, elats, albedo)
  
  USE OMSAO_precision_module 
  USE OMSAO_variables_module, ONLY: atmdbdir
  USE ozprof_data_module,     ONLY: atmos_unit
  IMPLICIT NONE
  
  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)         :: month
  REAL (KIND=dp), DIMENSION(2), INTENT(IN)  :: elons, elats  
  REAL (KIND=dp), INTENT(OUT) :: albedo

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER        :: nlat=180, nlon=288
  REAL (KIND=dp), PARAMETER :: longrid = 1.25, latgrid = 1.0
  CHARACTER (LEN=2), DIMENSION(12),   PARAMETER :: monstr = (/'01', '02', &
       '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'/)
  CHARACTER (LEN=130)            :: alb_fname
  CHARACTER (LEN=14)             :: tmpstr
  INTEGER, SAVE, DIMENSION(nlon, nlat) :: glbalb
  INTEGER                        :: i, j, k, latin, lonin, nact, npix
  LOGICAL                        :: file_exist
  LOGICAL, SAVE                  :: first = .TRUE.
  REAL (KIND=dp)                 :: sumalb, lat, lon
  
  IF (first) THEN
     ! determine lon and lat index
     alb_fname = TRIM(ADJUSTL(atmdbdir)) // 'albdb/nrsclim' // monstr(month) // '.dat'
     
     ! Determine if file exists or not
     INQUIRE (FILE= alb_fname, EXIST= file_exist)
     IF (.NOT. file_exist) THEN
        STOP 'Albedo file not exist: alb_fname'
     ENDIF
     
     OPEN (UNIT = atmos_unit, file=alb_fname , status = 'unknown')     
     DO i = 1, 5
        READ (atmos_unit, '(A)')
     END DO
     DO i = 1, nlat 
        READ (atmos_unit,'(11(25i3/),13i3)') (glbalb(j, i), j=1, nlon)
     ENDDO
     CLOSE (atmos_unit)
     first = .FALSE.
  ENDIF
     
  npix = NINT((elons(2) - elons(1)) / longrid)
  IF (npix == 0) npix = 1
  sumalb = 0.0D0; nact=0
  DO i = 1, npix
     
     IF (npix > 1) THEN
        lon = elons(1)  + (i - 1 + 0.5) * longrid
        lat = elats(1)  + (lon - elons(1)) / (elons(2)-elons(1)) * (elats(2) - elats(1))
     ELSE
        lon = (elons(1) + elons(2) ) / 2.0
        lat = (elats(1) + elats(2) ) / 2.0
     ENDIF

     lonin = INT((lon + 180.0) / longrid) + 1
     latin = INT((lat + 90.0)  / latgrid) + 1
     lonin = MOD(lonin, nlon)
     IF (latin > nlat) latin = nlat

     IF (glbalb(lonin, latin) >= 0.0) THEN
        sumalb = sumalb + glbalb(lonin, latin)
        nact = nact + 1
     ENDIF
  ENDDO

  IF (nact > 0) THEN
     albedo = sumalb / nact / 1000.0
  ELSE
     albedo = 0.10
  ENDIF
   
  RETURN
END SUBROUTINE get_toms_alb

! ========================================================================
! Obtain OMI MLER surface albedo developed by Omar Torres and Changwoo Ahn
! ========================================================================
SUBROUTINE get_omi_alb(month, day, elons, elats, albedo)
  
  USE OMSAO_precision_module 
  USE OMSAO_variables_module, ONLY: atmdbdir
  USE ozprof_data_module,     ONLY: atmos_unit, pos_alb
  IMPLICIT NONE
  
  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)         :: month, day
  REAL (KIND=dp), DIMENSION(2), INTENT(IN)  :: elons, elats  
  REAL (KIND=dp), INTENT(OUT) :: albedo

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER        :: nlat=180, nlon=360
  REAL (KIND=dp), PARAMETER :: longrid = 1.0, latgrid = 1.0
  CHARACTER (LEN=2), DIMENSION(12),   PARAMETER :: monstr = (/'01', '02', &
       '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'/)
  CHARACTER (LEN=3), DIMENSION(4), PARAMETER :: cwavs = (/'331', '340', '360', '380'/)
  REAL (KIND=dp),    DIMENSION(4), PARAMETER :: wavs  = (/331.0, 340.0, 360.0, 380.0/)            
  CHARACTER (LEN=130)             :: alb_fname
  CHARACTER (LEN=14)              :: tmpstr
  INTEGER, DIMENSION(2)           :: wavin, monin
  INTEGER                         :: i, j, k, ialb, latin, lonin,  npix, nact, nm, nw
  REAL (KIND=dp), DIMENSION(2)    :: wavfrac, monfrac
  LOGICAL                         :: file_exist
  REAL (KIND=dp)                  :: sumalb, lon, lat, frac
  INTEGER, DIMENSION(nlon)        :: tmpalb

  LOGICAL, SAVE                  :: first = .TRUE.
  REAL (KIND=dp), SAVE, DIMENSION(nlon, nlat) :: glbalb

  ! Read albedo database
  IF (first) THEN
     
     IF (day <= 15) THEN
        monin(1) = month - 1
        IF (monin(1) == 0) monin(1) = 12
        monin(2) = month
        monfrac(1) = (15.0 - day) / 30.0
        monfrac(2) = 1.0 - monfrac(1)
     ELSE 
        monin(2) = month + 1
        IF (monin(2) == 13) monin(2) = 1
        monin(1) = month
        monfrac(2) = (day - 15) / 30.0
        monfrac(1) = 1.0 - monfrac(2)
     ENDIF
     nm = 2

     DO i = 1, 4
        IF (wavs(i) >= pos_alb) EXIT
     ENDDO
     IF (i == 1) THEN
        nw = 1; wavin(1) = 1; wavfrac(1) = 1.0
     ELSE IF ( i == 4 .AND. pos_alb >= wavs(4) ) THEN
        nw = 1; wavin(1) = 4; wavfrac(1) = 1.0
     ELSE 
        nw = 2; wavin(1) = i - 1; wavin(2) = i 
        wavfrac(2) = (pos_alb - wavs(i-1)) / (wavs(i)-wavs(i-1))
        wavfrac(1) = 1.0 - wavfrac(2)
     END IF

     glbalb = 0.0

     DO i = 1, nw
        DO j = 1, nm
           alb_fname = TRIM(ADJUSTL(atmdbdir)) // 'omialb/omialb' // monstr(monin(j)) // cwavs(wavin(i)) // '.dat'
           
           ! Determine if file exists or not
           INQUIRE (FILE= alb_fname, EXIST= file_exist)
           IF (.NOT. file_exist) THEN
              WRITE(*, *) 'GET_OMI_ALB: albedo file does not exist!!!'; STOP
           ENDIF
           
           OPEN (UNIT = atmos_unit, file=alb_fname , status = 'unknown')     
           DO k = 1, 3
              READ (atmos_unit, '(A)')
           ENDDO
           frac = monfrac(j) * wavfrac(i)
           DO k = 1, nlat 
              READ (atmos_unit,'(14(25i4/),10i4,a14)') tmpalb, tmpstr
              glbalb(:, k) = glbalb(:, k) + tmpalb * frac
           ENDDO
           CLOSE (atmos_unit)
        ENDDO
     ENDDO

     first = .FALSE.
  ENDIF

  npix = NINT((elons(2) - elons(1)) / longrid)
  IF (npix == 0) npix = 1
  sumalb = 0.0D0; nact=0
  DO i = 1, npix
     IF (npix > 1) THEN
        lon = elons(1)  + (i - 1 + 0.5) * longrid
        lat = elats(1)  + (lon - elons(1)) / (elons(2)-elons(1)) * (elats(2) - elats(1))
     ELSE
        lon = (elons(1) + elons(2) ) / 2.0
        lat = (elats(1) + elats(2) ) / 2.0
     ENDIF
     
     lonin = INT((lon + 180.0) / longrid) + 1
     latin = INT((lat + 90.0)  / latgrid) + 1
     lonin = MOD(lonin, nlon)
     IF (latin > nlat) latin = nlat
     
     IF (glbalb(lonin, latin) >= 0.0) THEN
        sumalb = sumalb + glbalb(lonin, latin)
        nact = nact + 1
     ENDIF
  ENDDO

  IF (nact > 0) THEN
     albedo = sumalb / nact / 10000.0
  ELSE
     albedo = 0.10
  ENDIF

  RETURN
END SUBROUTINE get_omi_alb


! ========================================================================
! Obtain OMLER (OMI product from Kleipool et al 2008)
! ========================================================================
SUBROUTINE get_omler_alb(month, day, elons, elats, albedo)
  
  USE OMSAO_precision_module 
  USE OMSAO_variables_module, ONLY: atmdbdir
  USE ozprof_data_module,     ONLY: atmos_unit, pos_alb
  USE OMSAO_he5_module
  IMPLICIT NONE
  
  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)         :: month, day
  REAL (KIND=dp), DIMENSION(2), INTENT(IN)  :: elons, elats  
  REAL (KIND=dp), INTENT(OUT) :: albedo

  ! ======================
  ! Local variables
  ! ======================
  INTEGER (KIND=8), PARAMETER     :: nlat = 360, nlon = 720
  REAL (KIND=dp),   PARAMETER     :: longrid = 0.5, latgrid = 0.5
  INTEGER (KIND=4), PARAMETER     :: nwvl = 23, nmon = 12
  REAL (KIND=4), DIMENSION(nlat)  :: lats
  REAL (KIND=4), DIMENSION(nlon)  :: lons
  REAL (KIND=4), DIMENSION(nwvl)  :: wvls           
  CHARACTER (LEN=130)             :: alb_fname
  CHARACTER (LEN=14)              :: tmpstr
  INTEGER, DIMENSION(2)           :: wavin, monin
  INTEGER                         :: i, j, k, ialb, latin, lonin,  npix, nact, nm, nw
  REAL (KIND=dp), DIMENSION(2)    :: wavfrac, monfrac
  LOGICAL                         :: file_exist
  REAL (KIND=dp)                  :: sumalb, lon, lat, frac
  INTEGER (KIND=2), DIMENSION(nlon, nlat, nwvl, nmon)  :: tmpalb
  LOGICAL, SAVE                   :: first = .TRUE.
  REAL (KIND=dp), SAVE, DIMENSION(nlon, nlat) :: glbalb

  INTEGER :: alb_fid, grid_id, status
  INTEGER (KIND=8)               :: start = 0, stride = 1
  INTEGER (KIND=8), DIMENSION(4) :: start_alb = 0, stride_alb = 1
  INTEGER (KIND=8), DIMENSION(4) :: edge_alb = (/nlon, nlat, nwvl, nmon/)

  REAL (KIND=dp), PARAMETER      :: scale_factor = 0.001
 
  alb_fname = TRIM(ADJUSTL(atmdbdir)) // 'KNMI_OMIALB/OMI-Aura_L3-OMLER_2004m10-2007m10_v003-2008m0910t100324.he5'

  nw = 2
  nm = 2

  ! Read albedo database
  IF (first) THEN
     
     IF (day <= 15) THEN
        monin(1) = month - 1
        IF (monin(1) == 0) monin(1) = 12
        monin(2) = month
        monfrac(1) = (15.0 - day) / 30.0
        monfrac(2) = 1.0 - monfrac(1)
     ELSE 
        monin(2) = month + 1
        IF (monin(2) == 13) monin(2) = 1
        monin(1) = month
        monfrac(2) = (day - 15) / 30.0
        monfrac(1) = 1.0 - monfrac(2)
     ENDIF
     nm = 2

     ! Determine if file exists or not
     INQUIRE (FILE = alb_fname, EXIST = file_exist)
     IF (.NOT. file_exist) THEN
        WRITE(*, *) 'GET_OMLER_ALB: albedo file does not exist!!!'; STOP
     ENDIF

     alb_fid = HE5_GDopen(TRIM(ADJUSTL(alb_fname)), he5f_acc_rdonly)

     grid_id = HE5_GDattach(alb_fid, TRIM(ADJUSTL('EarthSurfaceReflectanceClimatology')))

     status = HE5_GDrdfld(grid_id, 'Latitude',  start, stride, nlat, lats)
     status = HE5_GDrdfld(grid_id, 'Longitude', start, stride, nlon, lons)
     status = HE5_GDrdfld(grid_id, 'Wavelength', start, stride, nwvl, wvls)

     DO i = 1, nwvl
        IF (wvls(i) >= pos_alb) EXIT
     ENDDO
     IF (i == 1) THEN
        nw = 1; wavin(1) = 1; wavfrac(1) = 1.0
     ELSE IF ( i == nwvl .AND. pos_alb >= wvls(nwvl) ) THEN
        nw = 1; wavin(1) = nwvl; wavfrac(1) = 1.0
     ELSE 
        nw = 2; wavin(1) = i - 1; wavin(2) = i 
        wavfrac(2) = (pos_alb - wvls(i-1)) / (wvls(i)-wvls(i-1))
        wavfrac(1) = 1.0 - wavfrac(2)
     ENDIF
          
     glbalb = 0.0
     status = HE5_GDrdfld(grid_id, 'MonthlySurfaceReflectance', &
          start_alb, stride_alb, edge_alb, tmpalb)

     DO i = 1, nw
        DO j = 1, nm
           frac = monfrac(j) * wavfrac(i)
           glbalb = glbalb + tmpalb(:,:,wavin(i),monin(j)) * scale_factor * frac
        ENDDO
     ENDDO

     status = HE5_GDdetach(grid_id)
     status = HE5_GDclose(alb_fid)

     first = .FALSE.
  ENDIF

  npix = NINT((elons(2) - elons(1)) / longrid)
  IF (npix == 0) npix = 1
  sumalb = 0.0D0; nact=0
  DO i = 1, npix
     IF (npix > 1) THEN
        lon = elons(1)  + (i - 1 + 0.5) * longrid
        lat = elats(1)  + (lon - elons(1)) / (elons(2)-elons(1)) * (elats(2) - elats(1))
     ELSE
        lon = (elons(1) + elons(2) ) / 2.0
        lat = (elats(1) + elats(2) ) / 2.0
     ENDIF
     
     lonin = INT((lon + 180.0) / longrid) + 1
     latin = INT((lat + 90.0)  / latgrid) + 1
     lonin = MOD(lonin, nlon)
     IF (latin > nlat) latin = nlat
     
     IF (glbalb(lonin, latin) >= 0.0) THEN
        sumalb = sumalb + glbalb(lonin, latin)
        nact = nact + 1
     ENDIF
  ENDDO

  IF (nact > 0) THEN
     albedo = sumalb / nact
  ELSE
     albedo = 0.10
  ENDIF

  RETURN
END SUBROUTINE get_omler_alb



! ========================================================================
! Read synthetic surface albedo
! ========================================================================
SUBROUTINE get_synt_alb(albedo)

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: maxchlen
  USE OMSAO_variables_module, ONLY: atmdbdir, currpix, currline
  USE ozprof_data_module,     ONLY: pos_alb
  USE SYNT_read_l1b,          ONLY: read_synt_alb, synt_db, synt_date, synt_rad_wavl, nwavel, &
                                    nxtrack, nytrack, synt_num
  IMPLICIT NONE

  !! ======================
  !! Input/Output variables
  !! ======================
  REAL (KIND=dp), INTENT(OUT) :: albedo

  !! ======================
  !! Local variables
  !! ======================
  REAL (KIND=4), DIMENSION(nwavel)  :: wvls
  CHARACTER (LEN=maxchlen)             :: alb_fname
  INTEGER, DIMENSION(2)           :: wavin
  INTEGER                         :: i, nw  !, j, k, ialb, latin, lonin,  npix, nact, nm, nw
  REAL (KIND=dp), DIMENSION(2)    :: wavfrac
  REAL                            :: alb
  LOGICAL                         :: file_exist
  REAL (KIND=dp)                  :: frac
  REAL , DIMENSION(nxtrack, nytrack, nwavel)  :: synt_alb
  LOGICAL, SAVE                   :: first = .TRUE.

  INTEGER :: status !alb_fid, grid_id

  !alb_fname = TRIM(ADJUSTL(synt_db)) // TRIM(synt_date) // '_00UTC/GEMS_20130715_surface_albedo.nc3'
  !alb_fname = TRIM(ADJUSTL(synt_db))//'aux/'//'GEMS_'//synt_date //'_albedoPresInfo.nc'  ! total data
  alb_fname = TRIM(ADJUSTL(synt_db))//'aux/'//'GEMS_'//synt_date //'_albedoPresInfo_'//synt_num//'.nc'  ! sliced data

! Read albedo database
  IF (first) THEN

! Determine if file exists or not
    INQUIRE (FILE = TRIM(alb_fname), EXIST = file_exist)
    IF (.NOT. file_exist) THEN
      WRITE(*, *) 'GET_OMLER_ALB: albedo file does not exist!!!'; STOP
    ENDIF
    CALL read_synt_alb (alb_fname, synt_alb, status)
    wvls = synt_rad_wavl (:, currpix, currline)
    DO i = 1, nwavel
      IF (wvls(i) >= pos_alb) EXIT
    ENDDO

    IF (i == 1) THEN
      nw = 1; wavin(1) = 1; wavfrac(1) = 1.0
    ELSE IF ( i == nwavel .AND. pos_alb >= wvls(nwavel) ) THEN
      nw = 1; wavin(1) = nwavel; wavfrac(1) = 1.0
    ELSE
      nw = 2; wavin(1) = i - 1; wavin(2) = i
      wavfrac(2) = (pos_alb - wvls(i-1)) / (wvls(i)-wvls(i-1))
      wavfrac(1) = 1.0 - wavfrac(2)
    ENDIF

    alb = 0.0
    DO i = 1, nw
      frac = wavfrac(i)
      alb = alb + synt_alb(currpix, currline, wavin(i)) * frac
    ENDDO

    first = .FALSE.
  ENDIF

  IF (alb >= 0.0) THEN
    albedo = alb
  ELSE
    albedo = 0.10
  ENDIF

  RETURN
END SUBROUTINE get_synt_alb
