!===================================================================================
! Author: Xiong Liu
!
! Date  : Created on Feb 10, 2004
!         Modified from Randall's routine in amfv4.sd/main.F for reading
!         GEOS-CHEM aerosols
!
! Modificcation History:
! Nov 27, 2004: (before, use 9612-9711 for other years; now: 9609-9708)
! Jan 26, 2005: (use GOES-4(96-97,99-00) fields and GEOS-3 (98( fields for individual years
!               and GOES-4 average fields for other years
!
! Purpose:
! This program reads aerosol extinction profiles:
! (a): Stratospheric aerosols from SAGE (Jan 1995- Aug 1999)
!      Reference: Bauman et al., 2003a, b, c
!       4 wavelengths, 385, 453, 525, 1020 nm (extrapolated to 0.26 and 0.35um based on Q)
!      28 latitude 5 degree bands from -70s to 70n
!      Altitude from 10 km to 40 km at every 1 km (valid values about tropopause)
! (b): GEOS-CHEM monthly mean aerosol extinction profiles (i.e., GOCART aerosols) for
!      Mineral dust, tropospheric sulfate, black carbon,  organic carbon, sea salt
!      (accumulation mode and coarse mode)
!      Reference:
!      Chin, M., P. Ginoux, S. Kinne, O. Torres, B. Holben, B.N. Duncan, 
!      R.V. Martin, J.A. Logan, A. Higurashi, and T. Nakajima, Tropospheric 
!      aerosol optical thickness from the GOCART model and comparisons with 
!      satellite and sunphotometer measurements, J. Atmos. Sci., 59, 461-483, 2002.
!
!      Dimension: GEOS-STRAT, GRID2x2.5x26 (2.5lon x 2lat x 26 layers)
!      INTEGER, PARAMETER :: IGLOB  = 144 (-181.25 -- 178.75) pixel center
!      INTEGER, PARAMETER :: JGLOB  = 91  (-89.0 --- 89.0)
!      INTEGER, PARAMETER :: LGLOB  = 26 
!
!      Sigma Pressure Coordinate
!      Levels
!      DATA MSIGMAE
!     &   / 1.0d0, .987871d0, .954730d0, .905120d0, .845000d0, .78d0,
!     &         .710000d0, .639000d0, .570000d0, .503000d0, .440000d0,
!     &         .380000d0, .325000d0, .278000d0, .237954d0, .202593d0,
!     &         .171495d0, .144267d0, .121347d0, .102098d0, .085972d0,
!     &         .072493d0, .061252d0, .051896d0, .037692d0, .019958d0,
!     &         .000000d0 / 

! Optical properties:
! (a) Tropospheric aerosols: 
! The optical properties (Q, w, phase moments) for each of these aerosols 
! are computed by R.V. Martin at four wave lengths (300, 400, 600, 999)
! Note: for wavelengths less than 300., extrapolated
! (b) Stratospheric aerosols
! Since it consists of sulfate (70.85%H2SO4), single scattering albedo=1
! Asymmetric factors and phase moments are computed for a range of reff and sigma
! at four different wavelengths, 0.260, 0.337, 0.525, 1.0230
!===================================================================================
     
! This subroutine is modified from Randall's routine in amfv4.sd/main.F for reading
! GEOS-CHEM aerosols
! Note in gasmoms, first moments are not stored, since they are always 1

SUBROUTINE READ_AEROSOL_PROF(year, month, lon, lat, ps, zs, tp, nz, &
     waves, nwave, useasy, gaext, gasca, gaasy, ngksec, gamoms, errstat)
  
  USE OMSAO_precision_module 
  USE OMSAO_variables_module, ONLY: atmdbdir, refdbdir
  USE ozprof_data_module,     ONLY: atmos_unit, nmom, maxgksec, &
       strat_aerosol, which_aerosol, maxmom
  USE OMSAO_errstat_module

  IMPLICIT NONE
  
  ! ======================
  ! Parameter Variables
  ! ======================
  INTEGER, PARAMETER :: NAER    = 6, MWL=10, NWL=4, NSLAYER=31, NLAT=91, NSLAT=28, NLON=144
  !INTEGER, PARAMETER :: NTLAYER = 26 (for using v4.26)
  INTEGER, PARAMETER :: NTLAYER = 21  !(for using GEOS-3 and 4 fields)
  REAL (KIND=dp), PARAMETER :: latgrid=2.0, latsgrid=5.0, longrid=2.5, lon0=-181.25, lat0=-91.0
  ! Tracer numbers for each aerosol type
  INTEGER, DIMENSION(NAER)              :: IND = (/4, 6, 9, 12, 15, 18/)
  CHARACTER (LEN=8)                     :: CATEGORY_IN = 'OD-MAP-$'
  CHARACTER (LEN=3), DIMENSION(12)      :: monstr = (/'jan', 'feb', &
       'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/) 
  
  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)  :: month, year, nz, tp, nwave, ngksec
  INTEGER, INTENT(OUT) :: errstat
  LOGICAL, INTENT(IN)  :: useasy
  REAL (KIND=dp), INTENT (IN)                             :: lon, lat
  REAL (KIND=dp), DIMENSION(nwave), INTENT(IN)            :: waves
  REAL (KIND=dp), DIMENSION(0:nz),  INTENT(IN)            :: zs, ps  ! index 0 for TOA
  REAL (KIND=dp), DIMENSION(nwave, nz), INTENT(OUT)       :: gaext, gasca, gaasy
  REAL (KIND=dp), DIMENSION(nwave, nz, 0:nmom, ngksec), INTENT(OUT) :: gamoms
  
  ! ======================
  ! LOCAL variables
  ! ======================
  CHARACTER (LEN=130)                                       :: aername
  REAL (KIND=dp), SAVE, DIMENSION(NLON, NLAT, NTLAYER, NAER):: tarsl
  REAL (KIND=sp), DIMENSION(NLON, NLAT, NTLAYER)            :: array
  REAL (KIND=dp), SAVE, DIMENSION(NWL, NLAT, NSLAYER)       :: sext, sasy
  REAL (KIND=dp), SAVE, DIMENSION(NWL, NLAT, NSLAYER, 0:maxmom, maxgksec) :: smoms
  REAL (KIND=dp), DIMENSION(NWL, NSLAYER)                   :: nsext, nsasy
  REAL (KIND=dp), DIMENSION(NWL, NSLAYER, 0:nmom, ngksec)   :: nsmoms
  REAL (KIND=dp), DIMENSION(0:NSLAYER)                      :: zstrat
  REAL (KIND=dp), SAVE, DIMENSION(MWL, NAER)                :: wl
  REAL (KIND=dp), SAVE, DIMENSION(NWL)                      :: wls
  REAL (KIND=dp), SAVE, DIMENSION(MWL, NAER)                :: raa, qext, assa, qasy
  REAL (KIND=dp), SAVE, DIMENSION(maxgksec, 0:maxmom, MWL, NAER)  :: phfcn
  INTEGER,        SAVE, DIMENSION(MWL)                      :: nwls
  REAL (KIND=dp), DIMENSION(ngksec, 0:nmom, NAER)           :: phsmoms
  REAL (KIND=dp), DIMENSION(NTLAYER, NAER)                  :: tprof
  REAL (KIND=dp), DIMENSION(0:NTLAYER, NAER)                :: ctprof
  REAL (KIND=dp), DIMENSION(0:nz, NAER)                     :: ntprof
  REAL (KIND=dp), DIMENSION(NSLAYER)                        :: sprof, sg
  REAL (KIND=dp), DIMENSION(NSLAYER, 0:nmom, ngksec)        :: sph
  REAL (KIND=dp), DIMENSION(0:NSLAYER)                      :: csprof
  REAL (KIND=dp), DIMENSION(0:nz)                           :: nsprof
  REAL (KIND=dp), DIMENSION(NAER)                           :: ext, waer, ext400, asy
  REAL (KIND=dp)                                            :: tau_in, xg, xg1, divalt, slope, frac
  REAL (KIND=dp), DIMENSION(2) :: latfrac, lonfrac
  INTEGER       , DIMENSION(2) :: lonin, latin
  INTEGER  :: i, j, k, n, m, tracer_in, low, high, ntlvl, nslow, nshigh, ib, nblat, nblon, faer, laer
  LOGICAL  :: file_exist

!     ! sigma coordinate for GEOS-STRAT
!     REAL (KIND=dp), DIMENSION (0: NTLAYER), PARAMETER :: sigma0 =      &
!          (/1.0d0, .987871d0, .954730d0, .905120d0, .845000d0, .78d0,& 
!          .710000d0, .639000d0, .570000d0, .503000d0, .440000d0,     &
!          .380000d0, .325000d0, .278000d0, .237954d0, .202593d0,     &
!          .171495d0, .144267d0, .121347d0, .102098d0, .085972d0,     &
!          .072493d0, .061252d0, .051896d0, .037692d0, .019958d0,.0d0 /)

  ! Correct coordinates (for geos3 fields)
  REAL (KIND=DP), DIMENSION(0:NTLAYER), PARAMETER:: pres3 = (/ &
       1.00000D+00, 9.97095D-01, 9.91200D-01, 9.81500D-01, 9.67100D-01, 9.46800D-01, &
       9.19500D-01, 8.84000D-01, 8.39000D-01, 7.83000D-01, 7.18200D-01, 6.47600D-01, &
       5.74100D-01, 5.00000D-01, 4.27800D-01, 3.59500D-01, 2.97050D-01, 2.41950D-01, &
       1.94640D-01, 1.55000D-01, 1.22680D-01, 9.69D-02/)

  ! for geos4 fields
  REAL (KIND=DP), DIMENSION(0:NTLAYER), PARAMETER:: ap4    = (/  &
       0.000000D0,   0.000000D0,   12.704939D0,  35.465965D0,  66.098427D0,  101.671654D0,  &
       138.744400D0, 173.403183D0, 198.737839D0, 215.417526D0, 223.884689D0, 224.362869D0, &
       216.864929D0, 201.192093D0, 176.929993D0, 150.393005D0, 127.837006D0, 108.663429D0, &
       92.365662D0,  78.512299D0,  56.387939D0, 40.175419D0/)
  REAL (KIND=DP), DIMENSION(0:NTLAYER), PARAMETER:: bp4    = (/  &
       1.000000D0,   0.985110D0,   0.943290D0,   0.867830D0, 0.764920D0,  0.642710D0,  &
       0.510460D0,   0.378440D0,   0.270330D0,   0.183300D0, 0.115030D0,  0.063720D0,  &
       0.028010D0,   0.006960D0,   0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0, &
       0.000000D0,   0.000000D0,   0.000000D0, 0.000000D0/)
  !p(L) = AP4(L) + BP4(L) * ps

  REAL (KIND=dp), DIMENSION (0: NTLAYER)    :: zsig, sigma 
  CHARACTER (LEN=2)  :: monc, yrc 
  LOGICAL, SAVE      :: first = .TRUE.
  INTEGER, SAVE      :: ntlayer1

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=17), PARAMETER    :: modulename = 'read_aerosol_prof'
  
  ! --------------------------------------------------------------
  ! Initialize output variables
  ! --------------------------------------------------------------
  errstat = pge_errstat_ok
  gaext = 0.0; gasca = 0.0; gaasy=0.0; gamoms = 0.0

  IF (which_aerosol == 0) THEN
     faer = 1; laer = NAER
  ELSE
     faer = which_aerosol; laer = which_aerosol
  ENDIF
  
  IF (first) THEN
     !---------------------------------------------------------------
     ! Read tropospheric aerosol profiles from the binary punch file
     !---------------------------------------------------------------
     WRITE(monc, '(I2.2)') month          ! from 9 to '09'   
     WRITE(yrc, '(I2.2)') MOD(year, 100)  ! from 1997 to '97'

     !   aername = TRIM(ADJUSTL(atmdbdir)) // 'tropaer/gocart_' // monstr(month) // yrc //'.dat'
     !   
     !   ! Determine if file exists or not
     !   INQUIRE (FILE= aername, EXIST= file_exist)
     !   IF (.NOT. file_exist .AND. month >= 9) THEN 
     !      WRITE(*,*) 'Warning: no trop. aerosol file found, use 96-97 aerosols!!!'
     !      aername = TRIM(ADJUSTL(atmdbdir)) // 'tropaer/gocart_' // monstr(month) // '96.dat' ! 9-12 in 96
     !   ELSE IF (.NOT. file_exist) THEN
     !      WRITE(*,*) 'Warning: no trop. aerosol file found, use 96-97 aerosols!!!'
     !      aername = TRIM(ADJUSTL(atmdbdir)) // 'tropaer/gocart_' // monstr(month) // '97.dat' ! 1-8 in 97
     !   ENDIF
     !   ntlayer1 = NTLAYER

     ! CRN, 19-Jan-2011: Have got new aerosols from M. Fu (via Thomas Kurosu) for 2005. Use them 
     ! for year 2005.
     ! CRN, 29-Feb-2011: Use aerosols from Chulkyu Lee for 2006 and later.
     IF (year == 2004) yrc = '05'
     IF (year >= 2006) yrc = '06'
      
     aername = TRIM(ADJUSTL(atmdbdir)) // 'gcaer/aer_' // yrc // monc // '.dat'
     ! Determine if file exists or not
     INQUIRE (FILE= aername, EXIST= file_exist)
     IF (.NOT. file_exist) THEN
        WRITE(*, *) 'Warning: no aerosol profile file found, use 96-97, 99-00 average!!!'
        aername = TRIM(ADJUSTL(atmdbdir)) // 'gcaer/aer_avg' // monc // '.dat'
     ENDIF

     OPEN (UNIT=atmos_unit, FILE=aername, STATUS='old')
     DO i = 1, 5
        READ(atmos_unit, *) 
     ENDDO
     IF (year == 1998) THEN
        ntlayer1 = NTLAYER
     ELSE
        ntlayer1  = 17
     ENDIF
        
     READ(atmos_unit, '(144E8.2)') ((((tarsl(i,j,k,n),i=1,NLON),j=1,NLAT),k=1,NTLAYER),n=1,NAER)
     CLOSE(atmos_unit)
      
     !---------------------------------------------------------------
     ! Read stratospheric aerosol profiles from files
     !---------------------------------------------------------------
     IF (strat_aerosol) THEN
        aername = TRIM(ADJUSTL(atmdbdir)) // 'strataer/strataerv_' // yrc // '_' // monc // '.dat'
        
        ! Determine if file exists or not
        INQUIRE (FILE= aername, EXIST= file_exist)
        IF (.NOT. file_exist) THEN
           WRITE(*,*) 'Warning: no strat. aerosol file found, use 98-99 aerosols!!!'
           aername = TRIM(ADJUSTL(atmdbdir)) // 'strataer/strataerv_bk_' // monc // '.dat'
        ENDIF
        
        wls = (/260.0, 337.0, 525.0, 1020./)
        OPEN (UNIT=atmos_unit, FILE=aername, STATUS='old')
        DO i = 1, NSLAT 
           DO j = 1, NSLAYER
              READ(atmos_unit, *) (sext(k, i, j), k=1, NWL), &
                   (((smoms(k, i, j, n, m), k=1, NWL),  m = 1, maxgksec), n=0, 64)
           ENDDO
        ENDDO
        CLOSE (atmos_unit)
        smoms(:, :, :, 65:maxmom, :) = 0.0
        sext = 10.D0**(sext)
        IF (useasy) sasy = smoms(:, :, :, 1, 1) / 3.0
     ENDIF
     
     !------------------------------------------------------------------------
     ! Read tropospheric aerosol parameters
     !------------------------------------------------------------------------
!     OPEN (UNIT=atmos_unit, FILE = TRIM(ADJUSTL(refdbdir)) // 'jv_spec.dat', &
!          FORM='formatted', STATUS='old')
!     OPEN (UNIT=atmos_unit, FILE = TRIM(ADJUSTL(refdbdir)) // 'aer_vproperty_kirchstetter1000.dat', &
!          FORM='formatted',STATUS='old')
     OPEN (UNIT=atmos_unit, FILE = TRIM(ADJUSTL(refdbdir)) // 'aer_vpropertynew1000.dat', &
          FORM='formatted',STATUS='old')
     READ (atmos_unit, *)      ! Header
     DO n = 1, NAER 
        READ(atmos_unit, *)    ! Aerosol label
        READ(atmos_unit, *) nwls(n)
        IF (nwls(n) > MWL) THEN
           WRITE(*, *) modulename, ' # Aerosol input wavelengths exceeds MWL. Increase MWL!!!'
           errstat = pge_errstat_error; RETURN           
        ENDIF
        DO i = 1, nwls(n)
           READ(atmos_unit, *) wl(i,n), qext(i,n), raa(i,n), assa(i,n), &
                (phfcn(1, j, i, n), j = 0, nmom)
           DO k = 2, maxgksec
              READ(atmos_unit, *) (phfcn(k, j, i, n), j = 0, nmom)
           ENDDO
        ENDDO
     ENDDO
     IF (useasy) qasy = phfcn(1, 1, :, :) / 3.0
     CLOSE ( atmos_unit ) 
     
     first = .FALSE.
  ENDIF
  
  
  !-----------------------------------------------------------------------
  ! Convert to weighted aerosol properties: w, taer, and phase function
  !-----------------------------------------------------------------------
  CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
       lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
  tprof(:, faer:laer) = 0.0
  DO i = 1, nblon
     DO j = 1, nblat 
        tprof(:, faer:laer) = tprof(:, faer:laer) + &
             tarsl(lonin(i), latin(j), :, faer:laer) * lonfrac(i) * latfrac(j)
     ENDDO
  ENDDO 
  
  ! The optical thickness for tropospheric aerosols is scaled 
  ! using extinction coefficient for tropospheric aerosol to provided
  ! optical thickness at 400 nm
  DO j = faer, laer
     DO k = 1, nwls(j)
        IF (wl(k, j) >= 400.0) EXIT
     ENDDO
     IF (k == 1) k = 2
     low = k-1; high = k
     
     xg     = (400.0-wl(low, j))/(wl(high, j)-wl(low, j))
     ext400(j) = qext(low, j) * (1.0 - xg) + qext(high, j) * xg
  ENDDO
  
  ! Get height coordinate for GEOS-CHEM
  ! sigma = sigma0 * ps(0); 
  IF (year == 1998) THEN
     sigma = pres3 * ps(0)
  ELSE
     sigma = ap4 + bp4 * ps(0)
  ENDIF
  
  !zsig(NTLAYER) = zs(nz)
  CALL BSPLINE(LOG(ps), zs, nz+1, LOG(sigma(0:NTLAYER)), zsig(0:NTLAYER),&
       NTLAYER+1, errstat)
  IF (errstat < 0) THEN
     WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat 
     errstat = pge_errstat_error; RETURN
  ENDIF
  
  ! Find the level just above tropopause from GEOS CHEM grid
  ntlvl = 1; divalt = MAX(zs(tp), 10.0)
  DO WHILE (zsig(ntlvl) < divalt .AND. ntlvl <= NTLAYER-1)
     ntlvl = ntlvl + 1
  ENDDO

  IF (strat_aerosol) THEN
     ! Index for stratospheric aerosols
     IF (lat <= -67.5) THEN
        nblat = 1; latin(1) = 1; latfrac(1) = 1.0
     ELSE IF (lat >= 67.5) THEN
        nblat = 1; latin(1) = NSLAT; latfrac(1) = 1.0
     ELSE
        nblat = 2; frac = (lat + 67.5) / latsgrid + 1
        latin(1) = INT(frac); latin(2) = latin(1) + 1
        latfrac(1) = latin(2) - frac; latfrac(2) = 1.0 - latfrac(1)
     ENDIF
     
     nsext = 0.0; nsasy = 0.0; nsmoms = 0.0
     DO ib = 1, nblat
        nsext  =  nsext + sext(:, latin(ib), :) * latfrac(ib)
        nsasy  =  nsasy + sasy(:, latin(ib), :) * latfrac(ib)
        nsmoms(1:NWL, 1:NSLAYER, 0:nmom, 1:ngksec) =  nsmoms(1:NWL, 1:NSLAYER, 0:nmom, 1:ngksec) &
             + smoms(1:NWL, latin(ib), 1:NSLAYER, 0:nmom, 1:ngksec) * latfrac(ib)
     ENDDO

     ! Find the level that is just below tropopause
     DO i = 0, NSLAYER
        zstrat(i) = 10.0 + i
     ENDDO
  
     nslow = 0
     DO WHILE (zstrat(nslow) <= divalt)
        nslow = nslow + 1
     ENDDO
     nslow = nslow - 1     ! if nslow < 0, then has gap in it (extrapolation)
     
     nshigh = tp
     DO WHILE (zs(nshigh) < zstrat(NSLAYER-1))
        nshigh = nshigh + 1
     ENDDO
     nshigh = nshigh - 1
  ENDIF
  
  DO i = 1, nwave
     DO j = faer, laer
        DO k = 1, nwls(j)
           IF (waves(i) <= wl(k, j)) EXIT
        ENDDO
        IF (k == 1) k = 2
        high = k; low = k - 1
     
        ! handling tropospheric part    
        xg      = (waves(i) - wl(low, j)) / (wl(high, j)-wl(low, j))
        ext(j)  = qext(low, j) * (1.0 - xg) + qext(high, j) * xg
        waer(j) = assa(low, j) * (1.0 - xg) + assa(high, j) * xg
     
        IF (useasy) THEN
           asy(j) = qasy(low, j) * (1.0 - xg) + qasy(high, j) * xg
        ELSE
           phsmoms(1:ngksec, 0:nmom, j) = phfcn(1:ngksec, 0:nmom, low, j)*(1.0 - xg) + &
                phfcn(1:ngksec, 0:nmom, high, j) * xg
        ENDIF
     ENDDO
     
     ! scaling aerosol optical thickness to 400 nm values
     ctprof(0, faer:laer) = 0.0    
     DO j = 1, ntlvl
        ctprof(j, faer:laer) = ctprof(j-1, faer:laer) + ext(faer:laer) &
             / ext400(faer:laer) * tprof(j, faer:laer)
     ENDDO
     
     DO j = faer, laer
        CALL INTERPOL(zsig(0:ntlvl), ctprof(0:ntlvl, j), ntlvl+1, zs(0:tp), &
             ntprof(0:tp, j), tp+1, errstat)
        IF (errstat < 0) THEN
           WRITE(*, *) modulename, ': INTERPOL error, errstat = ', errstat
           errstat = pge_errstat_error; RETURN
        ENDIF
     ENDDO
     ntprof(1:tp, faer:laer) = ntprof(1:tp, faer:laer) - ntprof(0:tp-1, faer:laer)
     
     DO j = 1, tp
        gaext(i, j) = SUM(ntprof(j, faer:laer)) 
        gasca(i, j) = SUM(ntprof(j, faer:laer) * waer(faer:laer)) 
        
        IF (useasy) THEN
           gaasy(i, j) = SUM(ntprof(j, faer:laer) * waer(faer:laer) * asy(faer:laer)) / gasca(i, j)
        ELSE        
           DO k = 0, nmom
              DO n = 1, ngksec
                 gamoms(i, j, k, n) = SUM(phsmoms(n, k, faer:laer) * ntprof(j, faer:laer) &
                      * waer(faer:laer)) / gasca(i, j)
              ENDDO
           ENDDO
        ENDIF   
     ENDDO
    
     IF (strat_aerosol) THEN
        ! handling stratospheric part
        IF (waves(i) > 525.0) THEN    
           low = 3; high = 4
        ELSE IF(waves(i) > 337.0) THEN
           low=2; high=3
        ELSE
           low=1; high=2
        ENDIF
       
        xg = (waves(i)-wls(low))/(wls(high)-wls(low))
        sprof  = nsext(low, :) * (1.0 - xg) + nsext(high, :) * xg
        
        IF (useasy) THEN
           sg  = nsasy(low, :) * (1.0 - xg) + nsasy(high, :) * xg
        ELSE
           sph = nsmoms(low, :, :, :) * (1.0 - xg) + nsmoms(high, :, :, :) * xg
        ENDIF
     
        csprof(0) = 0.0
        DO j = 1, NSLAYER
           csprof(j) = csprof(j-1) + sprof(j)
        ENDDO
     
        CALL INTERPOL(zstrat, csprof, NSLAYER+1, zs(tp:nshigh), &
             nsprof(tp:nshigh), nshigh-tp + 1, errstat) 
        IF (errstat < 0) THEN
           WRITE(*, *) modulename, ': INTERPOL error, errstat = ', errstat 
           errstat = pge_errstat_error; RETURN
        ENDIF
        gaext(i, tp+1:nshigh) = nsprof(tp+1:nshigh) - nsprof(tp:nshigh-1)
        gasca(i, tp+1:nshigh) = gaext(i, tp+1:nshigh)  ! single albedo 1 for sulfate
        
        IF (useasy) THEN	
           CALL INTERPOL(zstrat(0:NSLAYER-1), sg, NSLAYER, &
                zs(tp+1:nshigh), gaasy(i, tp+1:nshigh), nshigh-tp, errstat)
           IF (errstat < 0) THEN
              WRITE(*, *) modulename, ': INTERPOL error, errstat = ', errstat
              errstat = pge_errstat_error; RETURN
           ENDIF
        ELSE
           DO k = 0, nmom
              DO n = 1, ngksec 
                 CALL INTERPOL(zstrat(0:NSLAYER-1), sph(:, k, n), NSLAYER, &
                      zs(tp+1:nshigh), gamoms(i, tp+1:nshigh, k, n), nshigh-tp, errstat)
                 IF (errstat < 0) THEN
                    WRITE(*, *) modulename, ': INTERPOL error, errstat = ', errstat
                    errstat = pge_errstat_error; RETURN
                 ENDIF
              ENDDO
           ENDDO
        ENDIF
        
        ! Handling mesophere (i.e., above 41 km)
        slope = (LOG10(gaext(i, nshigh))-LOG10(gaext(i, nshigh-1))) / &
             (zs(nshigh)-zs(nshigh-1))
        DO j = nshigh+1, nz
           gaext(i, j) = 10.0D0**(slope * (zs(j) - zs(nshigh-1)) + LOG10(gaext(i, nshigh-1))) 
        ENDDO
        
        gasca(i, nshigh+1:nz) = gaext(i, nshigh+1:nz)
        IF (useasy) THEN
           gaasy(i, nshigh+1:nz) = gaasy(i, nshigh)
        ELSE
           DO j = 0, nmom
              DO k = 1, ngksec
                 gamoms(i, nshigh+1:nz, j, k) = gamoms(i, nshigh, j, k) 
              ENDDO
           ENDDO
        ENDIF
     ENDIF
  ENDDO   ! end wavelength loop
     
  RETURN
  
END SUBROUTINE READ_AEROSOL_PROF

