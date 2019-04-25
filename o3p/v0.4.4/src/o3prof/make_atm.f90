! =======================================================================
! Author: Xiong Liu                                                      !
! Date: Feb 9, 2004                                                      !
! Purpose : Make atmosphere profiles used in LIDORT calculation for      ! 
!           GOME ozone profile retrieval                                 !
! 1. ECMWF temperature profiles                                          !
! 2. TOMS EP monthly mean total ozone                                    !
! 3. TOMS V8 ozone profiles                                              !
! 4. NCAR-NCEP surface pressure and terrain pressure                     !
! 5. SAGE-II stratospheric aerosols and GOCART tropospheric aerosols     !
! 6. GOMECAT cloud information                                           !
! Note:                                                                  !
! Ozone retrieval is done at different grid from LIDORT calculation      !
! Retrieved o3 needs to be interpolated into fine grids at each iteration!
! Other information other than ozone are fixed                           !    
! =======================================================================

! xliu (Dec 06 2006): 
!     1. Added missing ozone above ~60 km based on AFGL profiles
!     2. Added one more layer above ~60 km
!     3. Organize the codes and made some minor changes 
! xliu (Dec 26 2006):
!     1. Get surface pressure from the temperature profile and surface altitude (from L1B)
!        Previously, surface pressure was from NCEP-NCAR data (2.5x2.5 spatial resolution,
!        which is too coarse for OMI data).

!xliu, 03/09/11
! Notes about the number of layers and top of the atmosphere
!  Need to have 50 layers to make the calculation accurate (at least for single scattering?)
!  Need to divide the first (from top) layer to more layers (maybe only for single scattering?)
!  May need to increase the top layer (currently 0.087). To change it, just replace 13.5 with
!  a larger number, and may also need to increase nref to a larger number.

SUBROUTINE MAKE_ATM ( year, month, day, ndiv, tauc, cfrac,                       &
     ctp, numk, toz, spres, tpres, atmosprof, ozprof, nup2p, sacldscl, errstat)  ! nbatm, nold added : geun
  
  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: p0, boltz, xmair, accgrav, ugc, avo, rearth, du2mol
  USE OMSAO_variables_module, ONLY : sza => the_sza_atm, database, n_refspec_pts, &
       the_surfalt, the_lon, the_lat, fitvar_rad_init, currline, atmdbdir
  USE OMSAO_indices_module,   ONLY : so2_idx, hcho_idx, bro_idx, bro2_idx, no2_t1_idx, so2v_idx, o2o2_idx
  USE ozprof_data_module,     ONLY : atmos_prof_fname, profunit, do_lambcld, &
       which_clima, strataod, stratsca, tropaod, tropsca, tropwaer, stratwaer, useasy, mflay, &
       aerosol, cloud, has_clouds, ncbp, nctp, ntp, aerwavs, nw=>actawin, fts, nmom,          &
       fps, fzs, fozs, frhos, gaext, gasca, gaasy, gamoms, gcq, gcw, gcasy,                   &
       gcmoms, the_cbeta, nflay, mgasprof, do_tracewf, gasidxs, fgasidxs, atmwrt,             &
       ngas, the_cld_flg, tracegas, maxlay, use_reg_presgrid, presgrid_fname,                 &
       adjust_trop_layer, fixed_ptrop, ntp0, pst0, nsfc, nfsfc, use_tropopause,               &
       scaled_aod, scale_aod, ngksec,  maxgksec, do_simu, radcalwrt, which_toz, norm_tropo3,  &
       so2zind, so2vprofn1p1, so2valts, atmos_unit, ozone_above60km, pv811, trpz

  USE OMSAO_errstat_module
  IMPLICIT NONE
  
  ! xliu:, 03/08/2011
  ! Switch from using 17 2.5x2 NCEP Reanalysis temperture to 26 level 1x1 NCEP FNL data
  !INTEGER, PARAMETER :: nv8=11, nold=28, nlecm=22, nref=71, nmpref=61, nmipas=121
  INTEGER, PARAMETER :: nv8=11, nlfnl=26, nlecm=31, nold=nlecm+6, nref=71, nmpref=61, nmipas=121
  
  ! Input variables
  INTEGER, INTENT(IN)          :: year, month, day, ndiv, numk
  INTEGER, INTENT(OUT)         :: errstat
  REAL (KIND=dp), INTENT(INOUT):: tauc, cfrac
  REAL (KIND=dp), INTENT(INOUT):: ctp, toz, spres, tpres
  
  ! Output variables
  INTEGER, DIMENSION(0:numk), INTENT(OUT)           :: nup2p
  REAL (KIND=dp), DIMENSION(numk), INTENT(OUT)      :: ozprof
  REAL (KIND=dp), DIMENSION(3, 0:numk), INTENT(OUT) :: atmosprof

  !xliu, 09/03/05 scaling a priori for layers below clouds to avoud smoothing even
  !for full cloudy conditions
  REAL (KIND=dp), DIMENSION(numk), INTENT(OUT)      :: sacldscl
  
  ! Local variables
  INTEGER                    :: i, j, k, n, ncld, ncldinc, np, nftp, ntemp, mnorstd, tmpntp
  INTEGER, DIMENSION(0:numk) :: indarr
  REAL (KIND=dp)             :: dlgp, cbp, mt, accr, mindiff, fndiv, tmp, lon, lat, tmpscl, &
       so2v_fwhm, so2v_z0, so2v_vcd, sfct, halfdz

  ! xliu:, 03/08/2011
  ! Switch from using 17 2.5x2 NCEP Reanalysis temperture to 26 level 1x1 NCEP FNL data
!  REAL (KIND=dp), DIMENSION (0:nold), PARAMETER :: pold0 = (/1013.25, 1000., 925., 850., &
!       700., 600., 500., 400., 300., 250., 200., 150., 100., 70., 50., 30., 20., 10., 7., &
!       5., 3., 2., 1., 0.70, 0.35, 0.25, 0.175, 0.125, 0.0874604/)
!  REAL (KIND=dp), DIMENSION (0:nold), PARAMETER :: told0 = (/288.2, 287.5, 283.3, 278.7, &
!       268.6, 261., 252., 241.4, 228.6, 223.1, 220.5, 216.5, 216.6, 216.6, 216.6, 220.5,  &
!       223.1, 227.7, 232.4, 239.3, 249.5, 257.9, 271.3, 270.7, 260.8, 250., 247.0, 242., 236.0/)
  REAL (KIND=dp), DIMENSION (0:nold), PARAMETER :: pold0 = (/1013.25, 1000., 975., 950., 925., 900., 850., &
       800., 750., 700., 650., 600., 550., 500., 450., 400., 350., 300., 250., 200., 150., 100., 70., 50., &
       30., 20., 10., 7., 5., 3., 2., 1., 0.70, 0.35, 0.25, 0.175, 0.125, 0.0874604/)
  REAL (KIND=dp), DIMENSION (0:nold), PARAMETER :: told0 = (/288.2, 287.5, 286.1, 284.7, 283.2, 281.8, &
       278.7, 275.5, 272.2, 268.6, 264.8, 260.8, 256.6, 251.9, 246.9, 241.4, 235.4, 228.6, 220.5, &
       216.5, 216.7, 216.7, 216.7, 217.2, 220.5, 223.1, 227.7, 232.4, 239.3, 249.5, 257.9, 271.3, &
       269.3, 256.8, 249.7, 242.5, 235.9, 229.1/)

  REAL (KIND=dp), DIMENSION (0:nold)         :: pold, told, zold
  REAL (KIND=dp), DIMENSION (1:nmipas)       :: mipasp, mipast, mipaso3
  REAL (KIND=dp), DIMENSION (4)              :: ptemp    
  REAL (KIND=dp), DIMENSION (0:nv8)          :: pv8, v8oz
  REAL (KIND=dp), DIMENSION (nv8)            :: v8temp
  REAL (KIND=dp), DIMENSION (0:numk)         :: umkp_save, umkp, umkt, umkz, umkoz, umkoz1
  REAL (KIND=dp), DIMENSION (0:nref)         :: ozref, refp
  REAL (KIND=dp), DIMENSION (0:mflay)        :: oznref
  REAL (KIND=dp), DIMENSION (2)              :: afgltmp
  REAL (KIND=dp), DIMENSION (0:24)           :: fixed_p, fixed_cumoz
  REAL (KIND=dp), DIMENSION (1:24)           :: fixed_oz
  REAL (KIND=dp), DIMENSION (1:60)           :: a1, a2, a3
  LOGICAL, SAVE                              :: first = .TRUE., first1 = .TRUE., is_pgrid = .true., isinc=.TRUE.
  REAL (KIND=dp), DIMENSION (0:nv8), SAVE    :: pv80
  REAL (KIND=dp), DIMENSION (0:maxlay), SAVE :: umkp0, umkz0

  LOGICAL :: use_input_spres = .TRUE., fixed_temp = .FALSE.

  !xliu, 09/23/05, indicator for a layer if it is above a cloud or not
  !If above a cloud then 1 else 0
  REAL (KIND=dp), DIMENSION (mflay)            :: acld  
  
  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=8), PARAMETER    :: modulename = 'make_atm'

  errstat = pge_errstat_ok
  
  ! pressure grid for temperature profile
  pold = pold0; pold(nold) = p0 * 2.0D0 ** (-13.5D0) 
 
    
  ! xliu, 03/08/11, switch from NCEP Reanalysis temperature fields to NCEP FNL temperature fields
  ! Use NCEP FNL for up to 10 mb and ECMWFT average between 7 and 1 mb
  ! NCEP FNL: 26 layers, 1x1 degree, ECMWFT: 2.5 x 2.5 degree
  ! So need to merge NCEP FNL and ECMWF data from different spatial grid
  ! Get temperature profiles
! CALL GET_NCEPT(year, month, day, the_lon, the_lat, told(1:nlecm))
  CALL GET_NCEPFNLT(year, month, day, the_lon, the_lat, told(1:nlfnl))
  CALL GET_ECMWFAVGT(month, day, the_lon, the_lat, told(nlfnl+1:nlecm))

  IF (ALL(told(1:nlecm) == 0.0)) told = told0
 
  ! Use TOMS V8 temperature climatology for 0.70, 0.35 mb 
  CALL GET_V8TEMP(month, day, the_lat, v8temp)
  told(nlecm+1) = v8temp(10); told(nlecm+2) = v8temp(11)

  ! Use MIPAS Temperature cimatology for the upper 4 layers)
  CALL GET_MIPASIG2T(month, day, the_lat, mipasp, mipast)
  mipasp = LOG(mipasp); ptemp = LOG(pold(nlecm+3:nlecm+6))
  CALL BSPLINE(mipasp, mipast, nmipas, ptemp, told(nlecm+3:nlecm+6), 4, errstat)

  IF (errstat < 0) THEN
     WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
     errstat = pge_errstat_error; RETURN
  ENDIF
  if (fixed_temp) told = told0

  ! xliu, 03/09/11, get surface temperature
  CALL GET_SFCT(year, month, day, the_lon, the_lat, sfct)
  ! xliu (12/05/2006): Change from 7 K / 90 mb to 5.65 K / 100 mb (~6.5 K / km)
  !told(0) = told(1) + (pold(0) - pold(1)) / 100.0 * 5.65
  ! xliu (08/17/2008): Use extrapolation 
  IF (use_input_spres) THEN
    IF (spres > pold(0) ) THEN
       pold(0) = spres; told(0) = sfct
    ELSE
       told(0) = told(1) + (told(1)-told(2))/(LOG(pold(1))-LOG(pold(2))) * (LOG(pold(0))-LOG(pold(1)))
    ENDIF
  ELSE
     told(0) = told(1) + (told(1)-told(2))/(LOG(pold(1))-LOG(pold(2))) * (LOG(pold(0))-LOG(pold(1)))
  ENDIF

  ! Compute altitude grids for this temperature profile
  ! xliu, 03/09/11: Slightly modify the way of computing altitude grids (based on 
  !       surface pressure and surface altitude). Also use surface temperature 
  !       from met. data.
  pold = LOG(pold); zold(0) = 0; ntemp = 0
  IF (use_input_spres) THEN
     IF (spres == pold(0)) THEN
        zold(0) = the_surfalt
        ntemp = 0
     ELSE
        !CALL BSPLINE(pold, told, nold+1, LOG(spres), tmp, 1, errstat)
        !IF (errstat < 0) THEN
        !   WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
        !   errstat = pge_errstat_error; RETURN
        !ENDIF

        ! Use surface T and T at P(0) to determine the altitude difference
        ! Maynot be very accurate for high moutainous regions

        !mt = (told(0) + tmp)/ 2.0
        !dlgp = pold(0) - LOG(spres)
        !accr = (rearth / (rearth + the_surfalt/2.))**2
        !zold(0) = the_surfalt - boltz/(xmair/avo)/(accgrav*accr)*mt*dlgp

        DO i = 1, nold
           IF ( pold(i) <= LOG(spres) ) EXIT
        ENDDO
        ntemp = i - 1
        
        mt = (told(ntemp) + sfct)/ 2.0
        dlgp = pold(ntemp) - LOG(spres)
        halfdz = -(dlgp/2.30259) * 8.0 ! Approximate delta z using zs=-16. * alog10(Ps/P0)
        accr = (rearth / (rearth + the_surfalt + halfdz))**2
        zold(ntemp) = the_surfalt - boltz/(xmair/avo)/(accgrav*accr)*mt*dlgp
        
        DO i = ntemp-1, 0, -1
           mt = (told(i) + told(i+1))/ 2.0
           dlgp = pold(i) - pold(i+1)
           halfdz = -(dlgp/2.30259) * 8.0 ! Approximate delta z using zs=-16. * alog10(Ps/P0)
           accr = (rearth / (rearth + zold(i+1) + halfdz))**2
           zold(i) = the_surfalt - boltz/(xmair/avo)/(accgrav*accr)*mt*dlgp
        ENDDO
     ENDIF
  ENDIF
        
  DO i = ntemp+1, nold
     mt=(told(i)+told(i-1))/ 2.0
     dlgp= pold(i) - pold(i-1)
     halfdz = -(dlgp/2.30259) * 8.0 ! Approximate delta z using zs=-16. * alog10(Ps/P0)
     accr=(rearth / (rearth + zold(i-1) + halfdz))**2 	  
     zold(i)=zold(i-1)-boltz/(xmair/avo)/(accgrav*accr)*mt*dlgp
  ENDDO

  IF ( first1 ) THEN
     ! ================================ Generate Ozone Retrieval Layers ===========================
     IF (use_reg_presgrid) THEN
        ! 12 layers from 1.0 atm to 0.24 atm (equal log-pressure) as a basis except for the last layer
        fndiv = 12.0 / numk
        DO i = 0, numk
           umkp0(i) = 2.D0 ** (-fndiv * i)
        ENDDO
        umkp0(numk) = 2.0D0 ** (-13.5D0)  ! 65 km xliu (Dec 06 2006)
        umkp0(0:numk) = LOG(umkp0(0:numk) * p0)
     ELSE
        OPEN (UNIT=profunit, FILE=TRIM(ADJUSTL(presgrid_fname)), status='old', IOSTAT=errstat)
        IF ( errstat /= pge_errstat_ok ) THEN
           WRITE(*, *) modulename, ': Error in opening pressure grid file!!!'; STOP
        ENDIF
        READ (profunit, *) ntemp, is_pgrid
        IF (ntemp /= numk + 1) THEN
           WRITE(*, *) modulename, ': The number of layers is inconsistent!!!'; STOP
        ENDIF
        IF (is_pgrid) THEN 
           READ (profunit, *) umkp0(0:numk) 
        ELSE 
           READ (profunit, *) umkz0(0:numk)
        ENDIF
        CLOSE (profunit)
        
        IF (is_pgrid) THEN
           IF (umkp0(numk) > p0 * 2.0D0 ** (-13.5D0)) umkp0(numk) = p0 * 2.0D0 ** (-13.5D0) 
           umkp0(0:numk) = LOG(umkp0(0:numk))
        ENDIF
     ENDIF
     
     ! V8 TOMS Climatology
     IF (which_clima == 6) THEN
        DO i = 0, nv8 - 1
           pv80(i) = 2.0D0 ** (-i)
        ENDDO
        pv80(nv8) = umkp0(numk)
        pv80(0:nv8) = LOG(pv80(0:nv8) * p0)
     ENDIF
     
     IF (.NOT. fixed_ptrop .AND. .NOT. use_tropopause .AND. is_pgrid) THEN
        pst0 = EXP(umkp0(ntp0)); tpres = pst0
     ENDIF
        
     first1 = .FALSE.
  ENDIF
  pv8 = pv80
  pv811 = EXP(pv8) ! used in update_o3_sao3
  IF ( .NOT. is_pgrid) THEN
     CALL INTERPOL(zold, pold, nold+1, umkz0(0:numk), umkp0(0:numk), numk+1, errstat)
     IF (errstat < 0) THEN
        WRITE(*, *) modulename, ': INTERPOL error, errstat = ', errstat
        errstat = pge_errstat_error; RETURN
     ENDIF
     IF ( EXP(umkp0(numk)) > p0 * 2.0D0 ** (-13.5D0) ) umkp0(numk) = LOG(p0 * 2.0D0 ** (-13.5D0))

     IF (.NOT. fixed_ptrop .AND. .NOT. use_tropopause) THEN
        pst0 = EXP(umkp0(ntp0)); tpres = pst0
     ENDIF
  ENDIF
  umkp(0:numk) = umkp0(0:numk)

  ! Find surface pressure, adjust P, Z ,T fields if terrain height is below zero (e.g. dead sea)
  IF (.NOT. use_input_spres) THEN
     IF (the_surfalt < 0.0) THEN
        zold(0) = the_surfalt
        pold(0) = pold(1) + (pold(0) - pold(1)) * (zold(1) - zold(0)) / zold(1)
        !told(0) = told(1) + (EXP(pold(0)) - EXP(pold(1))) / 100.0 * 5.65
        told(0) = told(1) + (told(1)-told(2))/(LOG(pold(1))-LOG(pold(2))) * (LOG(pold(0))-LOG(pold(1)))
        spres = EXP(pold(0))
     ELSE IF (the_surfalt == 0.0) THEN
        spres = EXP(pold(0))
     ELSE
        DO i = 1,  nold 
           IF (the_surfalt <= zold(i)) THEN
              spres = pold(i) + (pold(i-1) - pold(i)) *  (zold(i) - the_surfalt) / (zold(i) - zold(i-1))
              spres = EXP(spres); EXIT
           ENDIF
        ENDDO
     ENDIF
  ENDIF

  ! Find the level of tropopause
  ! ntp = 3 means there are three layers below tropopause
  IF (.NOT. use_tropopause .AND. .NOT. fixed_ptrop) THEN
     ntp = ntp0
  ELSE
     ! Replace the closest level (ranges from 85 mb to 450 mb) with tropopause pressure 
     ntp = MAXVAL(MINLOC(ABS(LOG(tpres)-umkp))) - 1;  umkp(ntp) = LOG(tpres)
     IF ( ntp == 0 ) ntp = 1 ! should not replace surface (xliu: 08/09/12) 
  ENDIF
  
  IF (adjust_trop_layer) THEN
     ! Replace first layer with surface pressure 
     ! Redivide below-tropopause layers to equal-altitude layers
     umkp(0) = LOG(spres)
     dlgp = (umkp(0) - umkp(ntp)) / REAL(ntp, KIND=dp)
     DO i = 1, ntp-1
        umkp(i) = umkp(i-1) - dlgp
     ENDDO
     nsfc = 0
  ELSE
     ! Find the layer that are closest to surface pressure
     !nsfc = MAXVAL(MINLOC(ABS(LOG(spres)-umkp(0:ntp)))) - 1;  umkp(nsfc) = LOG(spres)

     ! Find the first layer that are below surface (used in TES)
     IF (umkp(0) < LOG(spres)) THEN
        nsfc = 0
     ELSE
        DO i = 1, ntp
           IF (umkp(i) < LOG(spres)) EXIT
        ENDDO
        IF (1.05 * spres / EXP(umkp(i-1)) > EXP(umkp(i)) / spres) THEN
           nsfc = i - 1
        ELSE
           nsfc = i
        ENDIF
     ENDIF
     
     ! Replace that layer with surface pressure
     ! mark number of layers below surface pressure
     umkp(nsfc) = LOG(spres)
  ENDIF


 
  !=================================================================================!
  !  GET ozone profile in retrieval grid
  !=================================================================================!

  refp(0:nref) = (/(p0 * 10.0D0 ** (-REAL(i, KIND=dp)/16.D0), i=0, nref )/)
  ! Get a priori climatology for 60-70 km from MIPAS climatology
  CALL GET_MIPASIG2O3(month, day, the_lat, mipasp, mipaso3)
  ozref(nmpref:nref-1) = mipaso3(nmpref:nref-1)
  ! Ozone above 70 km and one more layer
  ozref(nref) =  SUM(mipaso3(nref:nmipas))
  
  ! Get a priori climatology for 0-60 km (pressure altitude)
  ozone_above60km = SUM( ozref(nmpref:nref))   
  toz = toz - ozone_above60km
  ! loading troppause height in km

  CALL INTERPOL(pold, zold, nold+1,umkp(ntp), trpz, 1, errstat)

  IF (errstat < 0) THEN
     WRITE(*, *) modulename, ': INTERPOL error, errstat = ', errstat
     errstat = pge_errstat_error; RETURN
  ENDIF 
! 1 = surface   
  IF (which_clima == 1 .or. which_clima == 6) THEN 
         CALL GET_MCPROF (ozref(1:nmpref-1)) 
  ELSE IF  (which_clima >= 2 .and. which_clima <=4) THEN
         CALL GET_TBPROF (ozref(1:nmpref-1))
  ELSE  IF (which_clima == 5) THEN
         CALL GET_IUPPROF(toz, ozref(1:nmpref-1))   
  ELSE  IF (which_clima == 12) THEN 
         CALL GET_MLprof(ozref(1:nmpref-1), 1)
  ELSE  IF (which_clima == 13) THEn 
         CALL GET_TJprof(ozref(1:nmpref-1), 1)
  ELSE 
         CALL GET_MCPROF (ozref(1:nmpref-1)) 
  ENDIF
!    CALL GET_MCPROF (a1(1:60)) ;CALL GET_TBPROF (a2(1:60))  ;CALL GET_IUPPROF(toz, a3(1:60))
!           DO i = 1, 60 
!            print * , a1(i), a2(i), a3(i)
!            enddo
!           print * , sum(a1), sum(a2), sum(a3)

  ! Bondary layer correction
  IF (spres > p0 ) THEN !sfc
      tmp = ( spres - p0)/(refp(0)-refp(1))
      ozref(1) = ozref(1)*(1+tmp)
      refp(0) = spres
  ENDIF
  IF (refp(nref) > EXP(umkp(numk)) ) refp(nref) = EXP(umkp(numk)) !top

  ozref(0) = 0.0
  DO i = 1, nref
     ozref(i) = ozref(i-1) + ozref(i)
  ENDDO

  refp = LOG(refp)

  CALL BSPLINE(refp, ozref, nref+1, umkp(0:numk), umkoz(0:numk), numk+1, errstat)
  IF (errstat < 0) THEN
     WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat ;stop
  ENDIF

  IF (which_clima == 6) THEN
   
  CALL get_v8prof(toz, v8oz(1:nv8)) 
     v8oz(0) = 0.0   ! get cumulative ozone amount
     DO i = 1, nv8
        v8oz(i) = v8oz(i) + v8oz(i-1)
     ENDDO

     ! get cumulative ozone amount at modified grids
     IF (spres > p0 ) THEN 
        tmp = ( spres - p0)/(exp(pv8(0))-exp(pv8(1)))
        v8oz(1) = v8oz(1)*(1+tmp)
        pv8(0) = LOG(spres)
     ENDIF
        pv8(nv8) = umkp(numk)

     ! May have negative values for upper layers or lead to inproper interpolation 
     ! when interpolated from the 11 layer climatology and the partial column ozone 
     ! for these layers are very small, need readjustment using the McPeters Clima
     umkoz1 = umkoz

     CALL BSPLINE(pv8, v8oz, nv8+1, umkp, umkoz, numk+1, errstat)
     IF (errstat < 0) THEN
        WRITE(*, *) modulename, ': BSPLINE error 1, errstat = ', errstat
        errstat = pge_errstat_error; RETURN
     ENDIF
         
     DO i = numk, 0, -1
        IF (umkp(i) >= pv8(nv8-1)) EXIT
     ENDDO
     tmp = (umkoz(numk)-umkoz(i)) / (umkoz1(numk)-umkoz1(i))

     umkoz = umkoz - umkoz(0); umkoz(1:numk) = umkoz(1:numk) - umkoz(0:numk-1)  
     umkoz1 = umkoz1 - umkoz1(0); umkoz1(1:numk) = umkoz1(1:numk) - umkoz1(0:numk-1)  
     umkoz(i+1:numk) = umkoz1(i+1:numk) * tmp
  ELSE 
     umkoz = umkoz - umkoz(0) !; toz = umkoz(numk)
     umkoz(1:numk) = umkoz(1:numk) - umkoz(0:numk-1)  
  ENDIF
    
  ! ====================== Generate Fine Layers for Radiance Calculation =================
  ! Generate fine grids for ozone retrieval grids according to following rules
  ! Clear-sky: each ozone retrieval layer is broken to ndiv layers
  ! In presence of clouds and ctp/cbp are not at bounaries, ctp/cbp are added as levels
  ! [ncbp, nctp] clouds is located between ncbp and nctp levels
 
  ! xliu (Dec 06 2006): Add one more layer for the top layer  
  IF (numk * ndiv + 3 > mflay) THEN
     STOP 'MAKE_ATM: need to increas mflay!!!'
  ENDIF
  
  fndiv = REAL(ndiv, KIND=dp)
  fps(0) = umkp(0); np = 0
  DO i = 1, numk
     ! xliu, 03/08/11
     ! Need to add at least 3 more layers to reduce the error to below 0.1% at 270 nm
     ! If add one more layer, the error will be ~0.4% at 270 nm and decreases to 0.08% at 300 nm
     IF ( i == numk) fndiv = fndiv + 1  ! Add one more layers for the top layer
     dlgp = (umkp(i-1) - umkp(i)) / fndiv
     
     DO j = 1, fndiv
        np = np + 1
        fps(np) = umkp(i-1) - dlgp * j
    ENDDO
  ENDDO


  ! Determine CBP by assumming COD of 20 per 100 hPa
  ! several heuristic assumptions about clouds are implemented
  ! due to inperfect or inadequate cloud information
  ! xliu: 10/09/2005:
  !       1. change 200 m (20 hpa) to 30 m (3 hpa) if the given CTP is below surface
  !       2. change from 15 per 100 hPa to 20 Per 100 hPa (i.e., 6.6 to 5)
  has_clouds = .FALSE.; nctp = 0; ncbp = 0; cbp =0.0
  !ctp = 900 !wasp
  IF (cloud .AND. ctp > 0.0) THEN 
     fps = EXP(fps)
     has_clouds = .TRUE.

     IF ( do_lambcld ) THEN            ! Lambertian clouds (not need for cbp)
        IF (ctp >= spres  ) THEN
           CALL GET_TOMSV8_CTP(month, day, the_lon, the_lat, ctp, errstat) !wasp
           cfrac = 0.5  ! will be updated anyway at longer wavelength
           IF (errstat == pge_errstat_error) THEN
              WRITE(*, *) modulename, ': Error in getting OMI Cloud Climatology ', errstat
              RETURN
           ENDIF
           
           IF (ctp >= spres) ctp = spres - 20.0
           the_cld_flg = 2            ! Set up flag for adjusting cloud-top pressure
        ENDIF
     ELSE                              ! For scattering clouds
        cbp = ctp + tauc * 5.0       
        IF (cbp >= spres - 20.0 ) THEN ! extreme case when ctp > Fps
           cbp = spres - 20.0          ! Cloud base height: ~ 30 m
           ctp = cbp - tauc * 5.0      ! COD 20 per 100 mb
           the_cld_flg = 3             ! Set up flag for adjusting cloud-top pressure
        ENDIF
     ENDIF

     ! Insert cloud-bottom pressure only necessary when scattering clouds   
     IF (.NOT. do_lambcld) THEN     
        DO i = 0, np-1
           IF (ABS(fps(i)-cbp) <= 1.0D-3*cbp) THEN
              ncbp = i + 1; cbp = fps(i); EXIT
           ELSE IF (ABS(fps(i+1)-cbp) <= 1.0D-3*cbp) THEN
              ncbp = i + 2; cbp = fps(i+1); EXIT
           ELSE IF (fps(i) > cbp .AND. fps(i+1) < cbp) THEN  ! CBP here, add 1 layer
              fps(i+2:np+1) = fps(i+1:np)
              np = np + 1;  ncbp = i + 2; fps(i+1) = cbp; EXIT
           ENDIF
        ENDDO
     ELSE
        ncbp = 1
     ENDIF
     
     ! Insert cloud top pressure
     k = 0
     IF (.NOT. do_lambcld) k = ncbp - 1
     DO i = k, np-1
        IF (ABS(fps(i)-ctp) <= 1.0D-3 * ctp ) THEN
           nctp = i ;  ctp = fps(i); EXIT
        ELSE IF (ABS(fps(i+1)-ctp) <= 1.0D-3 * ctp ) THEN
           nctp = i+1 ; ctp = fps(i+1); the_cld_flg = 3; EXIT
        ELSE IF (fps(i) > ctp .AND. fps(i+1) < ctp) THEN     ! CTP here, add 1 layer
           fps(i+2:np+1) = fps(i+1:np)
           np = np + 1;  nctp = i + 1; fps(i+1) = ctp; EXIT
        ENDIF
     ENDDO
     
     fps = LOG(fps)  ! convert fps back to logarithmic for convenience

     ! need to break cloud to more layers ( tau = 2 / per layer) for scattering clouds
     IF (.NOT. do_lambcld) THEN
        !IF (tauc >= 50.0) THEN
        !   ncld = INT(tauc / 2.5)
        !ELSE IF (tauc >= 40.0) THEN
        !   ncld = INT(tauc / 2.0)
        !ELSE IF (tauc >= 30.0) THEN
        !   ncld = INT(tauc / 1.5)
        !ELSE IF (tauc >= 20.0) THEN
        !   ncld = INT(tauc / 1.0)
        !ELSE IF (tauc >= 10.0) THEN
        !   ncld = INT(tauc / 0.5)
        !ELSE
        !   ncld = INT(tauc / 0.25)
        !ENDIF
        ncld = 10
        IF (do_simu .AND. .NOT. radcalwrt) ncld = 20

        ncldinc = ncld - (nctp-ncbp+1)
        fps(nctp+ncldinc:np+ncldinc) = fps(nctp:np)
        nctp = nctp + ncldinc
        np = np + ncldinc

        ! divide those cloud layers to equal log-pressure layers
        dlgp = (fps(ncbp-1) - fps(nctp)) / REAL(nctp-ncbp+1, KIND=dp)
        DO i = ncbp, nctp-1
           fps(i) = fps (i-1) - dlgp
        ENDDO

        ! read just layers in umkp that is between clouds
        DO i = 0, numk
           IF (umkp(i) > LOG(cbp)) THEN
              CYCLE
           ELSE IF (umkp(i) < LOG(ctp)) THEN
              EXIT
           ELSE IF (umkp(i) >= LOG(ctp) .AND. umkp(i) <= LOG(cbp)) THEN
              mindiff = 1000.0
              DO j = ncbp-1, nctp
                 IF (ABS(umkp(i) - fps(j)) < mindiff) THEN
                    mindiff = ABS(umkp(i) - fps(j))
                    k = j
                 ENDIF
              ENDDO
              umkp(i) = fps(k)
           ENDIF
        ENDDO
     ENDIF
  ENDIF ! cloud 
 
  ! Get indices for ozone retrieval levels in fine grids
  j = 1; nup2p = 0
  DO i = 1, np
     IF (ABS(fps(i) - umkp(j)) <= ABS(umkp(j)*1.0D-7)) THEN
		nup2p(j) = i; j = j + 1
     ENDIF
  ENDDO
  nfsfc = nup2p(nsfc)
  
  ! ================================ Get T, P, Z, aerosols ============================== 
  ! Interpolate temperature and altitude to FINE grids
  CALL BSPLINE(pold, told, nold+1, fps(0:np), fts(0:np), np+1, errstat)
  IF (errstat < 0) THEN
     WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
     errstat = pge_errstat_error; RETURN
  ENDIF

  CALL INTERPOL(pold, zold, nold+1, fps(0:np), fzs(0:np), np+1, errstat)
  IF (errstat < 0) THEN
     WRITE(*, *) modulename, ': INTERPOL error, errstat = ', errstat
     errstat = pge_errstat_error; RETURN
  ENDIF
  fzs(nsfc) = the_surfalt
 
  ! Calculate air column density and altitude grid 
  ! Assuming hydrostatic and Z (MAX(1atm,spres)) = 0.0 km
  !IF (spres > p0) THEN
  !   fzs(0) = 0.0
  !ELSE
  !  mt = (told(0) + fts(0)) / 2.0
  !   dlgp= pold(0) - fps(0)    !fps is already in log()
  !   fzs(0) = boltz/(xmair/avo)/(accgrav)*mt*dlgp
  !END IF
  ! use terrain height database
  !CALL GET_SURFALT(the_lon, the_lat, fzs(0))
  !
  !fzs(0) = the_surfalt
  !DO i=1, np
  !   mt=(fts(i)+fts(i-1))/ 2.0
  !   dlgp= fps(i) - fps(i-1)
  !   accr=(rearth / (rearth + fzs(i-1)))**2 	  
  !   fzs(i)=fzs(i-1)-boltz/(xmair/avo)/(accgrav*accr)*mt*dlgp
  !ENDDO
   
  frhos(0:np) = EXP(fps(0:np)) * 100.0 / ugc / fts(0:np) * avo * 1.0D-6  ! number density
   
  frhos(1:np) = (frhos(1:np) + frhos(0:np-1)) / 2.0 * (fzs(1:np)-fzs(0:np-1)) * 1.0D5  
  ! Add molecules about fps(0): molecule ~ pressure
  frhos(np) = frhos(np) * EXP(fps(np-1)) / (EXP(fps(np-1))-EXP(fps(np)))
  
  umkt = fts(nup2p); umkz = fzs(nup2p)     !  T and Z grid for ozone retrieval layers

 
  ! ======================== Interpolate Ozone to LIDORT Grid ======================
  CALL BSPLINE(refp, ozref, nref+1, fps(0:np), oznref(0:np), np+1, errstat)
  IF (errstat < 0) THEN
     WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
     errstat = pge_errstat_error; RETURN
  ENDIF
  oznref(1:np) = oznref(1:np) - oznref(0:np-1); oznref(0) = 0.0  

  ! convert back to mbar for pressure
  umkp = EXP(umkp); fps(0:np) = EXP(fps(0:np))
  atmosprof(1, 0:numk) = umkp; atmosprof(2, 0:numk) = umkz; atmosprof(3, 0:numk)=umkt
  ozprof = umkoz(1:numk)

  ! Loading minor trace gas profile
  IF (do_tracewf .OR. first) THEN
     DO i = 1, ngas
       IF (fgasidxs(i) > 0) THEN
           IF (gasidxs(i) == no2_t1_idx) THEN
              CALL GET_NO2(year, month, the_lon, the_lat, fzs(nfsfc:np), fps(nfsfc:np), sza, mgasprof(i, nfsfc+1:np), np-nfsfc)
           ELSE IF (gasidxs(i) == bro_idx) THEN     
              CALL GET_BRO(month, the_lat, fzs(nfsfc:np), fps(nfsfc:np), sza, mgasprof(i, nfsfc+1:np), np-nfsfc)
              !print *, 'pressure: '
              !WRITE(*, '(10F10.4)') fps(nfsfc:np)
              !print *, 'BrO: '
              !WRITE(*, '(10D12.4)') mgasprof(i, nfsfc+1:np)
              !WRITE(*, '(D12.4)')  SUM(mgasprof(i, nfsfc+1:np))
              
           ELSE IF (gasidxs(i) == bro2_idx) THEN     
              !Put BrO in the first layer
              mgasprof(i, nfsfc+1:np) = 0.0
              mgasprof(i, nfsfc+1:nfsfc+nup2p(1)) = (fps(nfsfc:nfsfc+nup2p(1)-1) - fps(nfsfc+1:nfsfc+nup2p(1))) / &
                   (fps(nfsfc) - fps(nfsfc+nup2p(1)))  * 2.0E13
              !print *, 'BrO2: '
              !WRITE(*, '(10D12.4)') mgasprof(i, nfsfc+1:np)
              !WRITE(*, '(D12.4)')  SUM(mgasprof(i, nfsfc+1:np))
              !print *, nup2p(1)
              !STOP
              
           ELSE IF (gasidxs(i) == hcho_idx) THEN 
              CALL GET_GEOSCHEM_HCHO(month, the_lon, the_lat, fps(nfsfc:np), mgasprof(i, nfsfc+1:np), np-nfsfc) 
           ELSE IF ( gasidxs(i) == so2_idx) THEN 
              nftp = nup2p(ntp) - nfsfc
              !write(90, *) currline 
              CALL GET_GEOSCHEM_SO2(year, month, the_lon, the_lat, fps(nfsfc:np), mgasprof(i, nfsfc+1:np), nftp, np-nsfc) 
              !errstat = pge_errstat_error; RETURN
           ELSE IF (gasidxs(i) == so2v_idx) THEN
              so2valts(0)  = fitvar_rad_init(so2zind)
              IF (so2valts(0) < the_surfalt) THEN
                 fitvar_rad_init(so2zind) = the_surfalt + 2.5
                 so2valts(0) = fitvar_rad_init(so2zind)
              ELSE IF (so2valts(0) < the_surfalt + 1.0) THEN
                 so2valts(0) = so2valts(0) + 1.0
              ENDIF         
              so2valts(-1) = so2valts(0) - 1.0
              so2valts(1)  = so2valts(0) + 1.0  

              ! Initializing: a Gaussian plume with 0.5 km FWHM and 0.1 DU SO2
              so2v_fwhm = 1.0; so2v_vcd = 0.1  
              CALL INSERT_GAUSSIAN_PLUME(so2v_fwhm, so2v_vcd, so2valts(0), np-nfsfc, &
                   fzs(nfsfc:np), isinc, mgasprof(i, nfsfc+1:np), errstat)
              mgasprof(i, nfsfc+1:np) = mgasprof(i, nfsfc+1:np) * du2mol
              
              CALL INSERT_GAUSSIAN_PLUME(so2v_fwhm, so2v_vcd, so2valts(-1), np-nfsfc, &
                   fzs(nfsfc:np), isinc, so2vprofn1p1(nfsfc+1:np, 1), errstat)
              so2vprofn1p1(nfsfc+1:np, 1) = so2vprofn1p1(nfsfc+1:np, 1) * du2mol
              so2vprofn1p1(np+1, 1) = SUM(so2vprofn1p1(nfsfc+1:np, 1))
              CALL REVERSE(so2vprofn1p1(1:np, 1), np)

              CALL INSERT_GAUSSIAN_PLUME(so2v_fwhm, so2v_vcd, so2valts(1), np-nfsfc, &
                   fzs(nfsfc:np), isinc, so2vprofn1p1(nfsfc+1:np, 2), errstat)
              so2vprofn1p1(nfsfc+1:np, 2) = so2vprofn1p1(nfsfc+1:np, 2) * du2mol
              so2vprofn1p1(np+1, 2) = SUM(so2vprofn1p1(nfsfc+1:np, 2))  
              CALL REVERSE(so2vprofn1p1(1:np, 2), np)     
           ELSE IF (gasidxs(i) == o2o2_idx) THEN
              ! Assume 20.95% for O2
              ! Here Divide 1.0E10, and absorption coefficients multiply by 1.D10
              mgasprof(i, 1:np) = (frhos(1:np)*0.2095)**2 / (fzs(1:np)-fzs(0:np-1)) / 1.0D5/1.D10
           ELSE
              print *, 'Not Implemented!!!'
              errstat = pge_errstat_error; RETURN
           ENDIF

           mgasprof(i, 1:nfsfc) = 0.0
           mgasprof(i, np+1)    = SUM(mgasprof(i, nfsfc+1:np))
           IF (do_lambcld) THEN
              tracegas(i, 8)    = SUM(mgasprof(i, nfsfc+1:nctp)) / mgasprof(i, np+1)
           ELSE
              tracegas(i, 8)    = 0.0
           ENDIF
           CALL REVERSE(mgasprof(i, 1:np), np)


        ENDIF
     ENDDO
     IF (first) first = .FALSE.
  ENDIF
  
  IF (which_clima == 11 ) THEN
     mnorstd = 1
     CALL get_mlso3prof_single(year, month, day, the_lat, numk, mnorstd, umkp, umkz, ozprof, tmpntp, errstat)
     IF (errstat < 0) THEN
        WRITE(*, *) modulename, ': Error in getting MLS ozone profiles!!!'
        errstat = pge_errstat_error; RETURN
     ENDIF
  ELSE IF (which_clima == 10 ) THEN
     mnorstd = 1
     CALL get_mlso3prof(year, month, day, the_lat, numk, mnorstd, umkp, umkz, ozprof, tmpntp, errstat)
     IF (errstat < 0) THEN
        WRITE(*, *) modulename, ': Error in getting MLS ozone profiles!!!'
        errstat = pge_errstat_error; RETURN
     ENDIF
  ELSE IF (which_clima == 9) THEN      ! 72 x 46, Logan clima only
     CALL GET_LOGAN_CLIMA(month, the_lon, the_lat, umkp, ozprof, numk, ntp)  
  ELSE IF (which_clima == 8) THEN      ! 144 x 91, profile only
     CALL GET_GEOSCHEM_O31(month, the_lon, the_lat, umkp, ozprof, numk, ntp)
  ELSE IF (which_clima == 7) THEN      ! 18 x 12, profile
     CALL GET_GEOSCHEM_O3MEAN(month, the_lon, the_lat, umkp, ozprof, numk, ntp)
  ENDIF

!  ! Use fixed a priori 
!  OPEN(UNIT = atmos_unit, FILE = TRIM(ADJUSTL(atmdbdir)) // 'mpclima/fixed_apriori.dat', status='old')
!  READ(atmos_unit, *) 
!  READ(atmos_unit, *) fixed_p
!  READ(atmos_unit, *) fixed_oz
!  CLOSE (atmos_unit)
!  
!  IF (fixed_p(0) < umkp(0)) THEN
!     fixed_oz(1) = fixed_oz(1) * (umkp(0) - fixed_p(1)) / (fixed_p(0) - fixed_p(1))
!     fixed_p(0) = umkp(0)
!  ENDIF
!  
!  fixed_cumoz = 0.0
!  DO i = 1, 24
!     fixed_cumoz(i) = fixed_cumoz(i-1) + fixed_oz(i)
!  ENDDO
!
!  CALL BSPLINE(fixed_p, fixed_cumoz, 25, umkp(0:numk), umkoz(0:numk), numk+1, errstat)
!  IF (errstat < 0) THEN
!     WRITE(*, *) modulename, ': Error in interpolating fixed a priori!!!'
!     errstat = pge_errstat_error; RETURN
!  ENDIF
!  ozprof(1:numk) = umkoz(1:numk) - umkoz(0:numk-1) 
!  ozprof(1:numk)=(/0.1920, 0.1560, 0.2780, 0.5060, 0.9120, 1.5850, 2.5755,3.8795, &
!       5.4425,7.4515,10.1115,12.5915,14.0555,12.7130,6.1330,-0.2220,-0.9420,4.2495,&
!       11.6435, 4.4670, 2.8970, 2.4810,2.6700,3.0605/)
!  CALL REVERSE(ozprof(1:numk), numk)
  
  ! Need to normalize ozone profile (Strat. and trop.) or tropospheric ozone profile with input total ozone
  IF (which_toz /= 0 .AND. which_clima /= 6 .and. which_clima /= 5) THEN
     IF (which_clima /= 10 .AND. which_clima /= 11) tmpntp = ntp
     CALL GET_NORMTOZ(year, month, day, the_lat, toz, numk, tmpntp, umkp, umkz, ozprof, errstat)  
     IF (errstat < 0) THEN
        WRITE(*, *) modulename, ': Error in getting zonal mean total ozone!!!'
        errstat = pge_errstat_error; RETURN
     ENDIF
  ENDIF
  
  DO i = 1, np
     IF (fps(i) >= ctp-1.0E-7) THEN
        acld(i) = 0.0 
     ELSE 
        acld(i) = 1.0
     ENDIF
  ENDDO
  
  !xliu, 08/29/05 scaling a priori for layers below clouds to avoid smoothing even
  !for full cloudy conditions
  umkoz(1:numk) = ozprof(1:numk)
  DO i = 1, numk
     fozs(nup2p(i-1) + 1:nup2p(i)) = oznref(nup2p(i-1) + 1 : nup2p(i)) * umkoz(i)     &
          / SUM(oznref(nup2p(i-1) + 1 : nup2p(i)))
     
     sacldscl(i) = SUM(fozs(nup2p(i-1) + 1:nup2p(i)) * acld(nup2p(i-1) + 1:nup2p(i)))  &
          /  SUM(fozs(nup2p(i-1) + 1:nup2p(i)))
  ENDDO
    
  ! ============================ Get Aerosols ======================================
  !gaext=0.0; gasca=0.0; gaasy=0.0; gamoms=0.0 
  IF (aerosol) THEN
       CALL READ_AEROSOL_PROF(year, month, the_lon, the_lat, fps(nfsfc:np), fzs(nfsfc:np), &
            nup2p(ntp)-nfsfc, np-nfsfc, aerwavs(1:nw), nw, useasy, gaext(1:nw, nfsfc+1:np),       &
            gasca(1:nw, nfsfc+1:np), gaasy(1:nw, nfsfc+1:np), &
            ngksec, gamoms(1:nw, nfsfc+1:np, 0:nmom, 1:ngksec), errstat)
       IF (errstat == pge_errstat_error) RETURN
       tropwaer(1:nw) = 0.0; stratwaer(1:nw) = 0.0
       DO i = 1, nw
          !gaext(i, nfsfc+1:5) = 0.0
          !gaext(i, 7:nup2p(ntp)) = 0.0
          !gasca(i, nfsfc+1:5) = 0.0
          !gasca(i, 7:nup2p(ntp)) = 0.0

          tropaod(i)  = SUM(gaext(i, nfsfc+1:nup2p(ntp)))
          tropsca(i)  = SUM(gasca(i, nfsfc+1:nup2p(ntp)))
          strataod(i) = SUM(gaext(i, nup2p(ntp)+1:np))
          stratsca(i) = SUM(gasca(i, nup2p(ntp)+1:np))
          IF (tropaod(i)  /= 0.0) tropwaer(i)  = tropsca(i)  / tropaod(i)
          IF (strataod(i) /= 0.0) stratwaer(i) = stratsca(i) / strataod(i)
       ENDDO

       IF (scale_aod)  THEN
          tmpscl  = scaled_aod / tropaod(nw)
          tropaod = tropaod * tmpscl
          tropsca = tropsca * tmpscl
          gaext(1:nw, nfsfc+1:nup2p(ntp)) = gaext(1:nw, nfsfc+1:nup2p(ntp)) * tmpscl
          gasca(1:nw, nfsfc+1:nup2p(ntp)) = gasca(1:nw, nfsfc+1:nup2p(ntp)) * tmpscl
       ENDIF
    ELSE 
       gaext=0.0; gasca=0.0; gaasy=0.0; gamoms=0.0 
       strataod = 0.0; stratsca = 0.0; tropaod = 0.0; tropsca = 0.0
    ENDIF
  
  ! ========================= Get Clouds(Homo) Properties ===========================
  IF (has_clouds .AND. .NOT. do_lambcld) THEN
     !WRITE(*, *) nctp, ncbp, fzs(nctp), fzs(ncbp-1), fps(nctp), fps(ncbp-1), tauc
     the_cbeta = tauc /(fzs(nctp) - fzs(ncbp-1))
     CALL GET_CLOUD_MIPROP(aerwavs(1:nw), nw, useasy, gcq(1:nw), &
          gcw(1:nw), gcasy(1:nw), ngksec, gcmoms(1:nw, 0:nmom, 1:ngksec), errstat) 
     IF (errstat == pge_errstat_error) RETURN
  ELSE
     gcq =0.0; gcw =0.0; gcasy=0.0; gcmoms=0.0
  ENDIF

  !WRITE(90, '(3D14.6)') fts(0), fps(0), fzs(0)
  !DO i = 1, np
  !   WRITE(90, '(9D14.6)') fts(i), fps(i), fzs(i), &
  !        fozs(i), frhos(i), mgasprof(1, i), mgasprof(4, i), mgasprof(5, i), mgasprof(7, i)
  !ENDDO
  !STOP
  
  ! =============== Write to atmospheric profiles as LIDORT Input====================
  IF (atmwrt) THEN
     OPEN (UNIT=profunit, FILE=TRIM(ADJUSTL(atmos_prof_fname)), status='unknown')
     WRITE(profunit, '(3L3, 2I3)') aerosol, has_clouds, useasy, nmom, nw
     WRITE(profunit, '(10F10.2)')  aerwavs(1:nw)
     
     ! write atmosphere
     WRITE(profunit, '(3X, F7.2, 21x, D14.6, F9.4)') fts(0), fps(0), fzs(0)
     DO i = 1, np
        IF (useasy) THEN
           WRITE(profunit, '(I3, 2F7.2,2D14.6,F9.4,100D11.3)') &
                i, fts(i), fozs(i), frhos(i), fps(i), fzs(i), &
                (gaext(j,i),j=1,nw), (gasca(j,i),j=1,nw), (gaasy(j,i),j=1,nw)
        ELSE
           WRITE(profunit, '(I3, 2F7.2,2D14.6,F9.4,150D11.3)') i, fts(i), fozs(i), &
                frhos(i), fps(i), fzs(i), (gaext(j,i),j=1,nw), (gasca(j,i),j=1,nw), &
                (((gamoms(j, i, k, n), n=1, ngksec), k=0, nmom), j = 1, nw)
        ENDIF
     ENDDO
     
     ! Write cloud information
     IF (has_clouds) THEN 
        WRITE(profunit, '(3I3, 2F10.4)') nfsfc, ncbp, nctp, the_cbeta, cfrac
        WRITE(profunit, '(20F11.4)') (gcq(i),i=1, nw)
        WRITE(profunit, '(20F11.7)') (gcw(i),i=1, nw)
        WRITE(profunit, '(20F11.4)') (gcasy(i),i=1,nw)
        WRITE(profunit, '(150F11.4)')  (((gcmoms(i, j, k), k = 1, ngksec), j = 0, nmom), i=1, nw)
     ENDIF
     CLOSE(profunit)
  ENDIF

  nflay = np 
  ncbp = np - ncbp + 1; nctp = np - nctp + 1; nfsfc = np - nfsfc + 1


  ! Reverse (from bottom-up to top-down)
  CALL REVERSE(fts(0:np), np+1)
  CALL REVERSE(fps(0:np), np+1)
  CALL REVERSE(fzs(0:np), np+1)
  CALL REVERSE(fozs(1:np), np)
  CALL REVERSE(frhos(1:np), np)

  DO i = 1, nw
     CALL REVERSE(gaext(i, 1:np), np)
     CALL REVERSE(gasca(i, 1:np), np)

     IF (useasy) THEN
        CALL REVERSE(gaasy(i, 1:np), np)
     ELSE
        DO j = 0, nmom
           DO k = 1, ngksec 
              CALL REVERSE(gamoms(i, 1:np, j, k), np)
           ENDDO
        ENDDO
     ENDIF

     ! No need to revert cloud properties because they are homogenuous (i.e., same at each layer)
  ENDDO

  ! Reverse data on the retrieval grid
  CALL REVERSE(ozprof(1:numk), numk)
  DO i = 1, 3 
     CALL REVERSE(atmosprof(i, 0:numk), numk+1)
  ENDDO
  CALL REVERSE(sacldscl(1:numk), numk)
  indarr = (/((i), i=0, numk)/)
  nup2p(0:numk) = MAXVAL(nup2p(0:numk)) - nup2p(numk - indarr) 
  ntp = numk - ntp; nsfc = numk - nsfc 
  
  RETURN

END SUBROUTINE MAKE_ATM


SUBROUTINE INSERT_GAUSSIAN_PLUME(fwhm, vcd, z0, nz, zs, isinc, gprof, errstat)

  !fwhm:  Full with at half maximum (km)
  !vcd:   Vertical column density (DU)
  !z0:    Center altitude
  !nz:    Number of layers
  !zs:    Altitudes at each level from surface to top (km) (either increasing or decreasing)
  !gprof: Partial column at each provided layer
  
  USE OMSAO_precision_module
  USE OMSAO_errstat_module
  IMPLICIT NONE
  
  ! Input/output variables
  INTEGER, INTENT(IN)                          :: nz
  INTEGER, INTENT(OUT)                         :: errstat
  LOGICAL, INTENT(IN)                          :: isinc
  REAL (KIND=dp), INTENT (IN)                  :: fwhm, vcd, z0
  REAL (KIND=dp), DIMENSION (0:nz), INTENT(IN) :: zs
  REAL (KIND=dp), DIMENSION (nz), INTENT(OUT)  :: gprof

  ! Local variables
  INTEGER, PARAMETER                 :: nlvl = 101
  INTEGER                            :: i
  REAL (KIND=dp), DIMENSION (nlvl)   :: x, y
  REAL (KIND=dp), DIMENSION (nlvl-1) :: ycol
  REAL (KIND=dp), DIMENSION (nlvl+2) :: x2, cum
  REAL (KIND=dp), DIMENSION (nlvl+1) :: ycol2
  REAL (KIND=dp), DIMENSION (0:nz)   :: ocum
  REAL (KIND=dp)                     :: dx, hw1e, xsum, temp

  ! -------------------------------
  ! Name of this subroutine/module
  ! --------------------------------
  CHARACTER (LEN=21), PARAMETER      :: modulename = 'INSERT_GAUSSIAN_PLUME'

  errstat = 0
  hw1e = fwhm / 1.66551
  xsum = fwhm * 3.0
  dx   =   xsum / (nlvl - 1.0)
  DO i = 1, nlvl 
     x(i)  = -xsum * 0.5 + ( i - 1 ) * dx   ! DU / KM
     y(i) = EXP( -(x(i) / hw1e)**2 )
  ENDDO
  
  DO i = 1, nlvl - 1
     ycol(i) = (y(i) + y(i+1)) / 2.0 * dx
  ENDDO
  ycol = ycol * vcd / SUM(ycol)    ! partial column 
  x    = z0 + x 

  IF (.NOT. isinc) CALL REVERSE(x(1:nlvl), nlvl)

  ! Interpolate (linear) the Gaussian plume onto input grid
  x2(2:nlvl+1) = x(1:nlvl);  ycol2(2: nlvl) = ycol; ycol2(1) = 0.0; ycol2(nlvl+1) = 0.0
  IF (isinc) THEN
     ycol2(2: nlvl) = ycol; ycol2(1) = 0.0; ycol2(nlvl+1) = 0.0
     x2(1)        = x2(2) - 1.0
     x2(nlvl + 2) = x2(nlvl+1) + 1.0
     IF (zs(0) < x(1) ) x2(1) = zs(0)
     IF (zs(nz) > x(nlvl)) x2(nlvl+2) = zs(nz)
  ELSE
     x2(1)        = x2(2) + 1.0
     x2(nlvl + 2) = x2(nlvl+1) - 1.0
     IF (zs(0)  > x(1) ) x2(1) = zs(0)
     IF (zs(nz) < x(nlvl)) x2(nlvl+2) = zs(nz)
  ENDIF
  
  cum(1) = 0.0
  DO i = 2, nlvl + 2
     cum(i) = cum(i-1) + ycol2(i-1)
  ENDDO

  CALL INTERPOL(x2(1:nlvl+2), cum(1:nlvl+2), nlvl+2, zs(0:nz), ocum(0:nz), nz+1, errstat)
  IF (errstat < 0) THEN
     WRITE(*, *) modulename, ': INTERPOL error, errstat = ', errstat
     errstat = pge_errstat_error; RETURN
  ENDIF
  gprof = ocum(1:nz) - ocum(0:nz-1)

  temp = MAXVAL(gprof)
  WHERE (ABS(gprof) < 1.0E-7*temp)
     gprof = 0.0
  ENDWHERE
  temp = SUM(gprof)

  IF (temp /= 0.0 .AND. temp /= vcd ) gprof = gprof * vcd / temp

  !WRITE(*, *) errstat
  !WRITE(*, '(10F12.4)')  zs(0:nz)
  !WRITE(*, '(10D12.4)')  gprof 
  !WRITE(*, *) SUM(gprof(1:nz))

  RETURN 

END SUBROUTINE INSERT_GAUSSIAN_PLUME

SUBROUTINE ADJUST_SO2VPLUMEZ(errstat)
  
  USE OMSAO_parameters_module, ONLY : mflay, du2mol
  USE OMSAO_indices_module,    ONLY : so2v_idx
  USE ozprof_data_module,      ONLY : mgasprof, so2vprofn1p1, fzs, nfsfc, so2valts, nflay, ngas, gasidxs
  USE OMSAO_precision_module
  USE OMSAO_errstat_module
  IMPLICIT NONE
  
  ! Input/output variables
  INTEGER, INTENT(OUT)                   :: errstat

  ! Local variables
  INTEGER                                :: iz, i, so2idx
  LOGICAL                                :: isinc = .FALSE.
  REAL (KIND=dp)                         :: so2v_fwhm, so2v_vcd
  REAL (KIND=dp), DIMENSION(mflay, -1:1) :: gprofs

  ! -------------------------------
  ! Name of this subroutine/module
  ! --------------------------------
  CHARACTER (LEN=17), PARAMETER      :: modulename = 'ADJUST_SO2VPLUMEZ'

  errstat = 0
  so2v_fwhm = 1.0;  so2v_vcd  = 0.1

  DO iz = -1, 1
     CALL INSERT_GAUSSIAN_PLUME(so2v_fwhm, so2v_vcd, so2valts(iz), nfsfc-1, &
          fzs(0:nfsfc-1), isinc, gprofs(1:nfsfc-1, iz), errstat)
  ENDDO

  DO i = 1, ngas
     IF (gasidxs(i) == so2v_idx) so2idx = i
  ENDDO

  mgasprof(so2idx, 1:nfsfc-1) = gprofs(1:nfsfc-1, 0)  * du2mol
  so2vprofn1p1(1:nfsfc, 1)    = gprofs(1:nfsfc-1, -1) * du2mol
  so2vprofn1p1(1:nfsfc, 2)    = gprofs(1:nfsfc-1, 1)  * du2mol
     
  RETURN 

END SUBROUTINE ADJUST_SO2VPLUMEZ

