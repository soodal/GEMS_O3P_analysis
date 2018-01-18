SUBROUTINE gems_o3p_prep_cld ( ny, offline, pge_error_status)

  USE OMSAO_precision_module
  USE GEMS_O3P_gemsdata_module, ONLY: nxtrack, ntimes, gems_clouds, ncoadd, nfxtrack, nxbin, nybin
  USE OMSAO_errstat_module
  IMPLICIT NONE

  ! ================
  ! Output variables
  ! ================
    INTEGER, INTENT(IN) :: ny, offline
  ! ================
  ! Output variables
  ! ================
  INTEGER, INTENT (OUT)       :: pge_error_status

  ! ================
  ! Local variables
  ! ================
  INTEGER ::j,k, ix,iy, iix,iiy, nbin, sline,eline, nbx,nx,nt, errstat
  INTEGER (KIND=i4), PARAMETER  :: nbit = 16
  REAL :: tmpsum, count1,count0, scfr
  REAL (KIND=r4),    DIMENSION (nxtrack, ntimes) :: cfr, ctp
  INTEGER (KIND=i2), DIMENSION (nxtrack, ntimes) :: ctp1
  INTEGER (KIND=i2), DIMENSION (nxtrack, ntimes) :: qflag
  INTEGER (KIND=i2), DIMENSION (nxtrack, ntimes,0:nbit-1) :: flgbits
  CHARACTER (LEN=23), PARAMETER :: modulename = 'gems_o3p_read_l2_cloud'

  pge_error_status = pge_errstat_ok 
  nx  = nxtrack
  nt  = ntimes
  sline = offline + 1
  eline = offline + ny*nybin
  !--------------------------------
  ! cloud data before coading, ctp, cfr
  !--------------------------------
  CALL Gems_o3p_share_l2_cld(nx, nt, cfr, ctp1, qflag)

  ctp = ctp1*1.0
  
  
  ! convert 2bytes to 16 its
  Do ix = 1, nx
    CALL convert_2bytes_to_16bits (nbit, nt, qflag(ix, 1:nt), flgbits(ix, 1:nt, 0:nbit-1))
  ENDDO
  qflag(1:nx, :)  = flgbits(1:nx, :, 1) + flgbits(1:nx, :, 2) &
         + flgbits(1:nx, :, 5) + flgbits(1:nx, :, 7)  &
         + flgbits(1:nx, :, 9) + flgbits(1:nx, :, 11) 
  ! Fill in cloud top pressure values for bad pixels
  CALL fill_in_ctp (nx,nt, ctp, qflag)   


  !--------------------------------
  ! cloud data after coading, gems_cloud
  !--------------------------------
  
  nbin = nxbin*ncoadd
  nbx = nxtrack/nbin

  gems_clouds%cfr    (1:nbx, 1:ny) = 0.0 
  gems_clouds%ctp    (1:nbx, 1:ny) = 0.0 
  gems_clouds%qflags (1:nbx, 1:ny) = 1 
  
  DO ix = 1, nbx 
    DO iy = 1, ny

      iix = (ix - 1) * nbin + 1 
      iiy = offline + 1+ (iy - 1) * nybin 
      scfr = 0.0; count1 = 0.0; count0 = 0.0; tmpsum = 0.0

      DO j = iix, iix + nbin - 1      
          DO k = iiy, iiy + nybin - 1
             IF (cfr(j, k) >= 0.0) THEN
                 gems_clouds%cfr(ix, iy) = gems_clouds%cfr(ix, iy) + cfr(j, k)
                 count1 = count1 + 1.0
             ENDIF
             IF ( ctp(j, k) > 0.0 .AND. cfr(j, k) >= 0.0  ) THEN                 
                 gems_clouds%ctp(ix, iy) = gems_clouds%ctp(ix, iy) + LOG(ctp(j, k)) * cfr(j, k)
                 tmpsum = tmpsum + LOG(ctp(j, k)) 
                 scfr = scfr + cfr(j, k)
                 count0 = count0 + 1.0
                 ! print * , j,k, ctp(j, k), cfr(j, k)
             ENDIF
           ENDDO
      ENDDO
    
      IF (scfr /= 0.0) THEN       ! Weighted by Cloud Fraction
        gems_clouds%ctp(ix, iy) = EXP(gems_clouds%ctp(ix, iy) / scfr)
      ELSE IF (count0 > 0.0 ) THEN ! Simple average if cloud fraction is all zero
        gems_clouds%ctp(ix, iy) = EXP(tmpsum / count0)      
      ELSE
        gems_clouds%ctp(ix, iy) = 0.0
        gems_clouds%qflags(ix,iy) = 10 ! bad results (should not be used)
      ENDIF
          
      IF (count1 /= 0.0) THEN
        gems_clouds%cfr(ix, iy) = gems_clouds%cfr(ix, iy) / count1
      ELSE
        gems_clouds%cfr(ix, iy) = 0.0
      ENDIF

       print * , scfr, count0, gems_clouds%ctp(ix, iy), ctp1(ix, iy)
     ENDDO ! Loop of Y-track
  ENDDO ! Loop of X-track

 
  RETURN
  END SUBROUTINE gems_o3p_prep_cld

  
  

 SUBROUTINE fill_in_ctp(nx, nt, ctp, qflag)
   USE OMSAO_precision_module
   ! ----------------------
    ! Input/Output variables
    ! ----------------------  
    INTEGER (KIND=i4), INTENT (IN)                      :: nx, nt 
    INTEGER (KIND=i2), DIMENSION (nx, nt), INTENT(INOUT) :: qflag
    REAL    (KIND=r4), DIMENSION (nx, nt), INTENT(INOUT) :: ctp
    
    ! -----------------------
    ! Local variables
    ! -----------------------
    INTEGER   (KIND=i4)      :: ix, i, j, k, nline, fidx, lidx, sidx, eidx
    REAL      (KIND=r4)      :: frac, dis




    ! If cloud top pressure is too small, then flag that retrieval
    WHERE (ctp < 90.0 ) 
       qflag = qflag + 1
    ENDWHERE
    
    ! If cloud top pressure is too large, then reset it to 1013 mb (i.e., surface pressure)
    WHERE (ctp > 1013.0 ) 
       ctp = 1013.0
    ENDWHERE
    
    ! Reset all flagged pixels to zero 
    WHERE (qflag >= 1) 
       ctp = 0.0
    ENDWHERE
    !PRINT *, 1.D0 * COUNT(ctp == 0.0) / (1.D0 * nx * nt * 1.0)

    ! Fill in flagged pixels
    DO i = 2, nt - 1   ! Fill in along the track
       DO ix = 1, nx 
          IF ( ctp(ix, i) == 0.0 .AND. qflag(ix, i-1) == 0 .AND. qflag(ix, i+1) == 0 ) &
               ctp(ix, i) = ( ctp(ix, i-1) + ctp(ix, i+1) ) / 2.0
       ENDDO
    ENDDO
    
    DO i = 1, nt    ! Fill in across the track 
       DO ix = 2, nx - 1 
          IF ( ctp(ix, i) == 0.0 .AND. qflag(ix - 1, i) == 0 .AND. qflag(ix + 1, i) == 0 ) &
               ctp(ix, i) = ( ctp(ix - 1, i) + ctp(ix + 1, i) ) / 2.0
       ENDDO
    ENDDO
    !PRINT *, 1.D0 * COUNT(ctp == 0.0) / (1.D0 * nx * nt * 1.0)
    
    ! Linear interpolation along the track
    DO ix = 1, nx      
       DO i = 1, nt
          IF (ctp(ix, i) > 0.0) EXIT
       ENDDO
       fidx = i
       
       DO i = nt, 1, -1
          IF (ctp(ix, i) > 0.0) EXIT
       ENDDO
       lidx = i
       
       IF (fidx >= lidx ) CYCLE
       
       i = fidx + 1
       DO WHILE ( i <= lidx )
          IF (ctp(ix, i) == 0.0 ) THEN
             sidx = i - 1; i = i + 1
             
             eidx = sidx - 1
             DO WHILE (i <= lidx) 
                IF ( ctp(ix, i) > 0.0 ) THEN
                   eidx = i; i = i + 1; EXIT
                ELSE
                   i = i + 1
                ENDIF
             ENDDO
             
             dis = REAL((eidx - sidx), KIND=dp)
             IF (dis <= 12) THEN
                DO j = sidx + 1, eidx - 1
                   frac = 1.0 - REAL((j - sidx), KIND=dp) / dis
                   ctp(ix, j) = frac * ctp(ix, sidx) + (1.0 - frac) * ctp(ix, eidx)
                ENDDO
             ENDIF
          ELSE
             i = i + 1
          ENDIF
       ENDDO
    ENDDO
    !print *, 1.D0 * COUNT(ctp == 0.0) / (1.D0 * nx * nt * 1.0)

    ! Linear interpolation across the track
    DO i = 1, nt      
       DO ix = 1, nx
          IF (ctp(ix, i) > 0.0) EXIT
       ENDDO
       fidx = ix

       DO ix = nx, 1, -1
          IF (ctp(ix, i) > 0.0) EXIT
       ENDDO
       lidx = ix

       IF (fidx >= lidx ) CYCLE
       
       ix = fidx + 1
       DO WHILE ( ix <= lidx )
          IF (ctp(ix, i) == 0.0 ) THEN
             sidx = ix - 1; ix = ix + 1
             
             eidx = sidx  - 1
             DO WHILE (ix <= lidx) 
                IF ( ctp(ix, i) > 0.0 ) THEN
                   eidx = ix; ix = ix + 1; EXIT
                ELSE
                   ix = ix + 1
                ENDIF
             ENDDO
             
             dis = REAL((eidx - sidx), KIND=dp)
             IF (dis <= 6) THEN
                DO j = sidx + 1, eidx - 1
                   frac = 1.0 - REAL((j - sidx), KIND=dp) / dis
                   ctp(j, i) = frac * ctp(sidx, i) + (1.0 - frac) * ctp(eidx, i)
                ENDDO
             ENDIF
          ELSE
             ix = ix + 1
          ENDIF
       ENDDO
    ENDDO
    !print *, 1.D0 * COUNT(ctp == 0.0) / (1.D0 * nx * nt)  

    RETURN
END SUBROUTINE fill_in_ctp
