! Cloud optical thickness is derived at approximately 760 nm
! gcq: ratio of extinction efficiency at wavelengths to that at 760 nm
! gcw: single scattering albedo = 1.0
! gcasy: asymmetric factor
! gcmoms: phase moments

SUBROUTINE GET_CLOUD_MIPROP(waves, nwave, useasy, gcq, gcw, gcasy, ngksec, gcmoms, errstat) 

  USE OMSAO_precision_module 
  USE OMSAO_variables_module, ONLY: refdbdir
  USE ozprof_data_module,     ONLY: atmos_unit, nmom, maxgksec, pos_alb, maxmom
  USE OMSAO_errstat_module
  
  IMPLICIT NONE
  
  ! ======================
  ! Parameter Variables
  ! ======================
  INTEGER, PARAMETER :: MWL = 10
  
  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)  :: nwave, ngksec
  INTEGER, INTENT(OUT) :: errstat
  LOGICAL, INTENT(IN)  :: useasy
  REAL (KIND=dp), DIMENSION(nwave), INTENT(IN)  :: waves
  REAL (KIND=dp), DIMENSION(nwave), INTENT(OUT) :: gcq, gcw, gcasy
  REAL (KIND=dp), DIMENSION(nwave, 0:nmom, ngksec), INTENT(OUT) :: gcmoms
  
  ! ======================
  ! LOCAL variables
  ! ======================
  REAL (KIND=dp), SAVE, DIMENSION(MWL)                     :: wl
  REAL (KIND=dp), SAVE, DIMENSION(MWL)                     :: raa, qext, assa, qasy
  REAL (KIND=dp), SAVE, DIMENSION(MWL, 0:maxmom, 1:maxgksec) :: phfcn
  LOGICAL,        SAVE                                     :: first = .TRUE.          
  INTEGER                                                  :: i, j, k, low, high, nwl
  REAL (KIND=dp)                                           :: extcld, xg
  
  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=17), PARAMETER    :: modulename = 'get_cloud_miprop'

  ! --------------------------------------------------------------
  ! Initialize to out variables
  ! --------------------------------------------------------------
  gcq = 0.0; gcw = 1.0; gcasy=0.0; gcmoms = 0.0
  errstat = pge_errstat_ok
  
  IF (first) THEN
     ! --------------------------------------------------------------
     ! Read cloud optical properties
     ! --------------------------------------------------------------
     OPEN (UNIT=atmos_unit, FILE=TRIM(ADJUSTL(refdbdir)) // 'cld_vprop1000.dat',&
          FORM='formatted',STATUS='old')
     READ(atmos_unit, *); READ(atmos_unit, *)      ! Header and cloud label
     READ(atmos_unit, *) nwl
        IF (nwl > MWL) THEN
           WRITE(*, *) modulename, ' # cloud input wavelengths exceeds MWL. Increase MWL!!!'
           errstat = pge_errstat_error; RETURN           
        ENDIF
     DO i = 1, NWL
        READ(atmos_unit, *) wl(i), qext(i), raa(i), assa(i), &
             (phfcn(i, j, 1), j = 0, nmom)
        assa(i) = 1.0  ! Assume a single scattering albedo of 1.0
        DO k = 2, maxgksec
           READ(atmos_unit, *) (phfcn(i, j, k), j = 0, nmom)
        ENDDO
     ENDDO

     CLOSE ( atmos_unit ) 
     IF (useasy) qasy = phfcn(:, 1, 1) / 3.0
     first = .FALSE.
  ENDIF
 
  ! --------------------------------------------------------------
  ! Interpolated to desired wavelengths
  ! --------------------------------------------------------------
  DO k = 1, nwl
     IF (pos_alb <= wl(k)) EXIT
  ENDDO
  IF (k == 1) k = 2
  high   = k; low = k - 1
  xg     = (pos_alb-wl(low))/(wl(high)-wl(low))
  extcld = qext(low) * (1.0 - xg) + qext(high) * xg
  
  DO i =1, nwave 
     DO k = 1, nwl
        IF (waves(i) <= wl(k)) EXIT
     ENDDO
     IF (k == 1) k = 2
     high = k; low = k - 1
        
     xg   = (waves(i)-wl(low))/(wl(high)-wl(low))
     gcq(i)  = qext(low) * (1_dp - xg) + qext(high) * xg
     gcw(i)  = assa(low) * (1_dp - xg) + assa(high)  * xg
     
     IF (useasy) THEN
        gcasy(i) = qasy(low) * (1.0 - xg) + qasy(high) * xg
     ELSE
        DO j = 0, nmom 
           DO k = 1, ngksec
              gcmoms(i, j, k) = phfcn(low, j, k) * (1.0 - xg) + phfcn(high, j, k) * xg
           ENDDO
        END DO
     ENDIF
  ENDDO
  
  gcq = gcq / extcld
  !gcw = gcw / extcld
  !WRITE(*, *) extcld, gcq
  !WRITE(*, *) gcasy
  !WRITE(*, *) gcw
  !STOP
  
  RETURN
END SUBROUTINE GET_CLOUD_MIPROP


! Get cloud macrophysical properties (CFRAC, COD, CTP) for Ch1(1.2s) and Ch2b(0.375s)
! For properties in Ch1, need to average over 6 pixels
! For properties in Ch2, use individual cld info for forward scan but average
! over three pixels for backward scan
! Variables
! year, mon, day, orb: year, month, day, which orbit in that day (e.g. 103)
! pix:  start pixels to be read
! cldflg: determines what cloud info is needed
!   = 0, channel 1, average over 6 pixels
!   = 1, individual pixel
!   = 2, backscan, average over three pixels
!   = 3 (i.e., 0 + 1)   (not implemented)
!   = 4 (i.e., 0 + 2)   (not implemented)
!   = 5 same as 0 except there are gome_npix = 40 
! sndpix: used for cldflg = 0, 3, 4
! sndpix = 1, 2, 3, 5, 6, 7, 8, 16 (8 for sum of 1, 2, 3;    16 for sum of 4, 5, 6)
! taucs, cfracs, ctps: returned 1/2 cloud optical thickness, fraction and 
! cloud top pressure (cloud optical thickness is scaled down to 1/4 for 
! consistency with ISCCP)


SUBROUTINE GET_CLOUD_MAPROP (year, month, day, orb, stpix, sndpix, cldflg, &
     ctaus, cfracs, ctps, errstat)
  
  USE OMSAO_precision_module 
  USE OMSAO_variables_module, ONLY: sol_identifier,atmdbdir
  USE ozprof_data_module,     ONLY: cldunit
  USE OMSAO_errstat_module

  IMPLICIT NONE
  !INTEGER, PARAMETER :: dp = KIND(1.0D0)
  
  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)           :: year, month, day, cldflg, stpix
  INTEGER, INTENT(OUT)          :: errstat
  INTEGER, INTENT(INOUT)        :: sndpix
  CHARACTER (LEN=3), INTENT(IN) :: orb
  REAL (KIND=dp), DIMENSION(2), INTENT(OUT) :: ctaus, cfracs, ctps
  
  
  ! ======================
  ! LOCAL variables
  ! ======================
  CHARACTER (LEN=3), DIMENSION(12),   PARAMETER :: monstr = (/'jan', 'feb', &
       'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/) 
  CHARACTER (LEN=4) :: yrc
  CHARACTER (LEN=2) :: dayc
  CHARACTER (LEN=3) :: actorb
  CHARACTER(LEN=130):: cldname
  CHARACTER(LEN=14) :: idstr
  LOGICAL           :: found, iflag, sflag, nothere, file_exist
  LOGICAL, SAVE     :: first = .TRUE.
  INTEGER           :: i, j, nump, stat, apix, pixtyp, nactp, nactp1, actyr, actmon
  REAL(KIND=dp)     :: pmdcfrac, oxcfrac, atau, actp, ifrac, sfrac, alb, &
       sumtau, sumctp, sumfrac, actz
  
  ! Initialize variables
  errstat = pge_errstat_ok
  ctaus = 0.0; ctps = 0.0; cfracs = 0.0

  IF (first) THEN
     ! get month, day, year from level 1 filename. Even though the actual date 
     ! could be different from the date on level 1 data (e.g. last orbit)
     
     ! Construct cloud file name to be searched
     READ(sol_identifier, '(I1,I2,A2,A3)') actyr, actmon, dayc, actorb
     
     IF (actyr >= 5) THEN
        actyr = actyr + 1990
     ELSE
        actyr = actyr + 2000
     ENDIF
     WRITE(yrc, '(I4)') actyr
     
     cldname = TRIM(ADJUSTL(atmdbdir)) // 'gomecat/' // yrc // '/' // monstr(actmon) &
          // '/GOME-' // yrc // '-' // monstr(month) // '-' // dayc // '-' // actorb &
          // '_R2_clouds.dat_short'  
     
     ! search file 
     INQUIRE (FILE= cldname, EXIST= file_exist)
     IF (.NOT. file_exist) THEN
        !don't use clouds, effective surface albedo derived
        WRITE(*, *) 'Warning: Cldfile not found!!!'; 
        errstat = pge_errstat_error; RETURN
     ENDIF
     OPEN (UNIT=cldunit, FILE=cldname, STATUS='old')
     
     first = .FALSE.
  ENDIF

  IF (cldflg == 1 .OR. cldflg == 2 .OR. cldflg == 5) sndpix = 0

  ! Determine number of pixels to be read
  IF (cldflg == 1) THEN
     nump = 1
  ELSE IF (cldflg == 2) THEN
     nump = 3
  ELSE IF (cldflg == 0) THEN
     nump = 6
  ELSE IF (cldflg == 5) THEN
     nump = 30
  ELSE
     WRITE(*, *) 'GET_CLOUD_MAPROP: Should not happen!!!'
     errstat = pge_errstat_error; RETURN
  ENDIF
    
  found = .FALSE.; stat=0; nothere=.FALSE.
  
  ! search until find the start pixel or EOF or error
  DO WHILE (.NOT. found .OR. stat /= 0 .OR. nothere)  
     
     READ (cldunit, FMT='(A14,2I6,8F10.4,1x,2L3)', IOSTAT=stat) idstr, apix,   &
          pixtyp, pmdcfrac, oxcfrac, actz, actp, atau, alb, ifrac, sfrac, &
          iflag, sflag
     
     IF (apix > stpix) THEN        ! no such pixel, since pixels are arranged in order
		nothere = .TRUE.
        BACKSPACE (cldunit)
        EXIT	
     ELSE IF (apix == stpix) THEN  ! pixel found
		found = .TRUE.
        
        IF (nump == 1) THEN
           IF (atau > 0) THEN
              ctaus(1) = atau
              cfracs(1) = oxcfrac
              ctps(1) = actp	
           ENDIF
		ELSE  ! multiple pixels
		   IF (atau > 0) THEN
              sumtau = atau
              sumfrac = oxcfrac
              sumctp = actp
              nactp = 1
           ELSE 
              sumtau =   0.0
              sumfrac =  0.0
              sumctp =   0.0
              nactp =    0
           ENDIF
           
           IF (sndpix == 1) THEN
              ctaus(2) = sumtau
              cfracs(2) = sumfrac
              ctps(2) = sumctp	
           ENDIF
           
           DO i = 2, nump
              READ(cldunit, FMT='(A14,2I6,8F10.4,1x,2L3)', IOSTAT=stat) idstr, &
                   apix, pixtyp, pmdcfrac, oxcfrac, actz, actp, atau, alb, &
                   ifrac, sfrac, iflag, sflag
              IF (stat /= 0) THEN
                 errstat = pge_errstat_error
                 CLOSE (cldunit)
                 RETURN   ! not eough cloud info, just exist assuming no clouds
              ENDIF

              IF (atau > 0) THEN
                 sumtau = sumtau + atau
                 sumfrac = sumfrac + oxcfrac
                 sumctp = sumctp + actp
                 nactp = nactp + 1
              ENDIF
              
              IF (sndpix == i .AND. atau > 0) THEN
                 ctaus(2) =  atau
                 cfracs(2) = oxcfrac
                 ctps(2) =   actp
              ENDIF
              
              IF (i == 3 .AND. sndpix > 6 .AND. nactp >= 1) THEN
                 ctaus(2) =  sumtau  / REAL(nactp, KIND=dp)
                 cfracs(2) = sumfrac / REAL(nactp, KIND=dp)
                 ctps(2) =   sumctp  / REAL(nactp, KIND=dp)
                 nactp1 = nactp
              ENDIF
           ENDDO
           
           IF (nactp > 0) THEN
              ctaus(1) =  sumtau  / REAL(nactp, KIND=dp)
              cfracs(1) = sumfrac / REAL(nactp, KIND=dp)
              ctps(1) =   sumctp  / REAL(nactp, KIND=dp)
           ENDIF
           
           IF (sndpix == 16) THEN
              IF (nactp > nactp1)  THEN
                 ctaus(2) =  (sumtau   - ctaus(2)  * nactp1) / &
                      REAL(nactp - nactp1, KIND=dp)
                 cfracs(2) = (sumfrac - cfracs(2) * nactp1) / &
                      REAL(nactp - nactp1, KIND=dp)
                 ctps(2) =   (sumctp   - ctps(2)   * nactp1) / &
                      REAL(nactp -nactp1, KIND=dp)
              ELSE
                 ctaus(2)  = ctaus(1)
                 cfracs(2) = cfracs(1)
                 ctps(2)   = ctps(1)
              ENDIF
           ENDIF
		ENDIF ! end second if
        
		EXIT ! for the case pixel is found
     ENDIF   ! end first if
  ENDDO      ! end while loop
  
  IF (ALL(ctaus <= 0.0D0)) THEN
     ctps = 0.0D0; cfracs = 0.0D0
  ELSE
     ! optical thickness is only used during the change of cloud-top height
     ! when the cloud-top height/cloud optical thickness does not look reasonable
     cfracs = cfracs / 100.0      ! convert from % 
  ENDIF

  ! special handling
  IF (ctaus(1) > 0.0 .AND. cfracs(1) < 0.0) cfracs(1) = 0.5
      
  RETURN
  
END SUBROUTINE GET_CLOUD_MAPROP

! =============================================================
! Obtain TOMSV7/ISCCP monthly mean cloud-top pressure
! =============================================================
SUBROUTINE GET_ISCCP_CTP(month, elons, elats, ctp, errstat)
  
  USE OMSAO_precision_module 
  USE OMSAO_variables_module, ONLY: atmdbdir
  USE ozprof_data_module,     ONLY: atmos_unit
  USE OMSAO_errstat_module

  IMPLICIT NONE
  
  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)         :: month
  INTEGER, INTENT(OUT)        :: errstat
  REAL (KIND=dp), DIMENSION(2), INTENT(IN)  :: elons, elats  
  REAL (KIND=dp), INTENT(OUT) :: ctp

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER :: numy = 8, nlat=72, nlon=144
  REAL (KIND=dp), PARAMETER :: longrid = 2.5, latgrid = 2.5
  CHARACTER (LEN=3), DIMENSION(12),   PARAMETER :: monstr = (/'JAN', 'FEB', &
       'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'/)
  CHARACTER (LEN=130)            :: isccp_fname
  CHARACTER (LEN=14)             :: tmpstr
  INTEGER, SAVE, DIMENSION(nlon, nlat) :: glbctp
  INTEGER                        :: i, j, k, latin, lonin, nact, npix
  LOGICAL                        :: file_exist
  LOGICAL, SAVE                  :: first = .TRUE.
  REAL (KIND=dp)                 :: sumctp, lat, lon

  errstat = pge_errstat_ok

  IF (first) THEN
     ! determine lon and lat index
     isccp_fname = TRIM(ADJUSTL(atmdbdir)) // 'isccp/CLDPRES' // '.' // monstr(month)
     
     ! Determine if file exists or not
     INQUIRE (FILE= isccp_fname, EXIST= file_exist)
     IF (.NOT. file_exist) THEN
        WRITE(*, *) 'ISCCP file not exist: ', isccp_fname
        errstat = pge_errstat_error; RETURN
     ENDIF
     
     OPEN (UNIT = atmos_unit, file=isccp_fname , status = 'unknown')     
     DO i = 1, 3
        READ (atmos_unit, '(A)')
     END DO
     DO i = 1, nlat 
        READ (atmos_unit,'(144I3)') glbctp(:, i)
     ENDDO
     CLOSE (atmos_unit)
     first = .FALSE.
  ENDIF

  npix = NINT((elons(2) - elons(1)) / longrid)
  IF (npix == 0) npix = npix + 1
  sumctp = 0.0D0; nact=0
  DO i = 1, npix

     IF (npix > 1) THEN
        lon = elons(1)  + (i - 1 + 0.5) * longrid
        lat = elats(1)  + (lon - elons(1)) / (elons(2)-elons(1)) * (elats(2) - elats(1))
     ELSE
        lon = (elons(1) + elons(2))/2.0
        lat = (elats(1) + elats(2))/2.0
     ENDIF

     lonin = INT((lon + 180.0) / longrid) + 1
     latin = INT((lat + 90.0)  / latgrid) + 1
     lonin = MOD(lonin, nlon)
     IF (latin > nlat) latin = nlat

     IF (glbctp(lonin, latin) > 0.0) THEN
        sumctp = sumctp + glbctp(lonin, latin)
        nact = nact + 1
     ENDIF
  ENDDO

  ctp = sumctp / nact
  
  RETURN
END SUBROUTINE GET_ISCCP_CTP


! =============================================================
! Obtain TOMSV7/ISCCP monthly mean cloud-top pressure
! =============================================================
SUBROUTINE GET_TOMSV8_CTP(month, day, lon, lat, ctp, errstat)
  
  USE OMSAO_precision_module 
  USE OMSAO_variables_module, ONLY: atmdbdir
  USE ozprof_data_module,     ONLY: atmos_unit
  USE OMSAO_errstat_module

  IMPLICIT NONE
  
  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)         :: month, day
  INTEGER, INTENT(OUT)        :: errstat
  REAL (KIND=dp), INTENT(IN)  :: lon, lat 
  REAL (KIND=dp), INTENT(OUT) :: ctp

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER             :: nlat=180, nlon=360, nmon=12
  REAL (KIND=dp), PARAMETER      :: longrid = 1.0, latgrid = 1.0, mongrid=1.0, &
       lon0=-180.0, lat0=-90.0, mon0=0.0
  CHARACTER (LEN=130)            :: isccp_fname
  CHARACTER (LEN=14)             :: tmpstr
  INTEGER, DIMENSION(2)          :: latin, lonin, monin 
  REAL (KIND=dp), DIMENSION(2)   :: latfrac, lonfrac, monfrac

  REAL (KIND=dp), SAVE, DIMENSION(nmon, nlon, nlat) :: glbctp
  INTEGER                        :: i, j, k, nblat, nblon, nbmon
  REAL (KIND=dp)                 :: mon
  LOGICAL                        :: file_exist
  LOGICAL, SAVE                  :: first = .TRUE.

  errstat = pge_errstat_ok
  IF (first) THEN

     ! determine lon and lat index
     !isccp_fname = TRIM(ADJUSTL(atmdbdir)) // 'isccp/CLDPRES_V8.txt'
     isccp_fname = TRIM(ADJUSTL(atmdbdir)) // 'isccp/omcldrr_clim_rev.dat'

     ! Determine if file exists or not
     INQUIRE (FILE= isccp_fname, EXIST= file_exist)
     IF (.NOT. file_exist) THEN
        WRITE(*, *) 'OMCLDRR file not exist: ', isccp_fname
        errstat = pge_errstat_error; RETURN
     ENDIF
     
     OPEN (UNIT = atmos_unit, file=isccp_fname , status = 'unknown')    
     DO i = 1, nmon
        DO j = nlat, 1, -1
           READ (atmos_unit, *) glbctp(i, :, j)
        ENDDO
     ENDDO
     CLOSE (atmos_unit)
     first = .FALSE.
  ENDIF

  mon = (month - 1.0) + day / 31.0

  CALL get_gridfrac1(nlon, nlat, nmon, longrid, latgrid, mongrid, lon0, lat0, mon0, &
     lon, lat, mon, nblon, nblat, nbmon, lonfrac, latfrac, monfrac, lonin, latin, monin)
  ctp = 0.0
  DO i = 1, nblon
     DO j = 1, nblat 
        DO k = 1, nbmon
           ctp = ctp + glbctp(monin(k), lonin(i), latin(j)) * lonfrac(i) * latfrac(j) * monfrac(k)
        ENDDO
     ENDDO
  ENDDO
  
  RETURN
END SUBROUTINE GET_TOMSV8_CTP
