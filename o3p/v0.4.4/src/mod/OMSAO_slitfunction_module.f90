MODULE OMSAO_slitfunction_module

  USE OMSAO_precision_module
  USE OMSAO_variables_module,     ONLY: refdbdir, coadd_uv2, currpix, band_selectors, numwin, winlim
  USE GEMS_O3P_gemsdata_module,   ONLY: nxtrack_max, nfxtrack, nxbin, ncoadd
  USE ozprof_data_module,         ONLY: calunit
  USE OMSAO_errstat_module
  IMPLICIT NONE

  INTEGER, PARAMETER :: NP  = 5, & ! # of parameterized parameters 
                        NP1 = 6, & ! # of actual slit parameters 
                        NO  = 7, & ! order of polynomials for parameterizing each parameter
                        NCH = 3    ! UV1, UV2, VIS
  REAL (KIND=r8), DIMENSION(NCH), PARAMETER:: omiscale = (/0.01, 0.001, 0.001/)
  REAL (KIND=r8), DIMENSION(NCH), PARAMETER:: omishi   = (/0.0,  350.0, 450.0/) 
  REAL (KIND=r8), DIMENSION(NCH, nxtrack_max, NP, NO)  :: omi_slitpars

CONTAINS

  SUBROUTINE load_slitpars (pge_error_status)

    ! Input/Output variables
    INTEGER, INTENT(OUT)          :: pge_error_status

    CHARACTER(LEN=3), DIMENSION(NCH), PARAMETER :: sltables = (/'uv1', 'uv2', 'vis'/)
    CHARACTER (LEN=maxchlen)                    :: slitpar_fname
    INTEGER                                     :: i, j, ic, ix, nx, errstat
    CHARACTER (LEN=12), PARAMETER               :: modulename = 'load_slitpars'

!    ! For testing the slit parameters
!    INTEGER :: ch, xpix, nw
!    REAL (KIND=r8) :: lamda
!    REAL (KIND=r8), DIMENSION(401) :: wave, slit
    
    DO ic = 1, NCH
       slitpar_fname = TRIM(ADJUSTL(refdbdir)) // 'OMI/omiopf_' // sltables(ic) // 'slit_xtrack.dat' 
       OPEN (UNIT=calunit, FILE=TRIM(ADJUSTL(slitpar_fname)), STATUS='UNKNOWN', IOSTAT=errstat)
       IF ( errstat /= pge_errstat_ok ) THEN
          WRITE(*, '(2A)') modulename, ': Cannot open slit parameterization file!!!'
          pge_error_status = pge_errstat_error; RETURN
       END IF

       DO i = 1, 10
          READ(calunit, *)
       ENDDO
       READ (calunit, *) nx
       IF (nx /= nxtrack_max .AND. nx /= nxtrack_max / 2) THEN
          WRITE(*, '(2A)') modulename, ': # of across track positions does not match!!!'
          pge_error_status = pge_errstat_error; RETURN
       ENDIF
       
       DO ix = 1, nx
          READ(calunit, *)
          READ(calunit, *) ((omi_slitpars(ic, ix, i, j), j=1, NO), i=1, NP)
       ENDDO
       CLOSE(calunit)    
    ENDDO

!    ! Test the slit parameters
!    ch = 3;   nw = 401; xpix = 60
!    lamda = 450.0
!    wave(:) = (/((j * 0.01 - 2. + lamda), j = 0, 400)/)
!    CALL compute_slitprofile (ch, xpix, lamda, wave, nw, slit)

    ! Binning UV-2/VIS slit parameters
    IF (coadd_uv2) THEN
       DO ix = 1, nfxtrack
          i = (ix - 1) * ncoadd + 1
          omi_slitpars(2:3, ix, :, :) = omi_slitpars(2:3, i, :, :) 
          
          DO j = 1, ncoadd-1      
             omi_slitpars(2:3, ix, :, :) = omi_slitpars(2:3, ix, :, :) + omi_slitpars(2:3, i + j, :, :)
          ENDDO
          omi_slitpars(2:3, ix, :, :) = omi_slitpars(2:3, ix, :, :) / ncoadd
       ENDDO
    ENDIF

    ! furthur binning across the track position
    IF (nxbin > 1) THEN
       DO ix = 1, nfxtrack / nxbin
          i = (ix - 1) * nxbin + 1
          omi_slitpars(:, ix, :, :) = omi_slitpars(:, i, :, :) 
          
          DO j = 1, nxbin - 1      
             omi_slitpars(:, ix, :, :) = omi_slitpars(:, ix, :, :) + omi_slitpars(:, i + j, :, :)
          ENDDO
          omi_slitpars(:, ix, :, :) = omi_slitpars(:, ix, :, :) / nxbin
       ENDDO
    ENDIF

    RETURN

  END SUBROUTINE load_slitpars

  SUBROUTINE compute_slitprofile (ch, lamda, wave, nw, slit)

    ! Input/Output variables   
    INTEGER, INTENT(IN)                        :: ch, nw
    REAL (KIND=r8), INTENT(IN)                 :: lamda
    REAL (KIND=r8), DIMENSION(nw), INTENT(IN)  :: wave
    REAL (KIND=r8), DIMENSION(nw), INTENT(OUT) :: slit
  
    ! Local variables
    INTEGER                         :: i                 
    REAL (KIND=r8)                  :: lam, tmp
    REAL (KIND=r8), DIMENSION (NO)  :: lamarr
    REAL (KIND=r8), DIMENSION (NP)  :: pars1
    REAL (KIND=r8), DIMENSION (NP1) :: pars

    lam = omiscale(ch) * (lamda - omishi(ch))
    lamarr(1) = 1.0
    DO i = 2, NO
       lamarr(i) = lamarr(i-1) * lam
    ENDDO

    IF (ch == 1) THEN
      pars(3)  = SUM(lamarr(1:6)  * omi_slitpars(ch, currpix, 3, 1:6) )
      slit = EXP(- ((wave-lamda)/pars(3))**2.0 )
   ELSE
      pars1(1) = SUM(lamarr(1:6)  * omi_slitpars(ch, currpix, 1, 1:6) )     ! 5th
      pars1(2) = SUM(lamarr(1:NO) * omi_slitpars(ch, currpix, 2, 1:NO))     ! 6th
      DO i = 3, 5
         pars1(i) = SUM(lamarr(1:6) * omi_slitpars(ch, currpix, i, 1:6))    ! 5th
      ENDDO
      tmp = pars1(2) / (pars1(1) + pars1(4))
      pars(1) = pars1(1)
      pars(2) = lamda !+ pars1(4) * tmp
      pars(3) = pars1(3)
      pars(4) = pars1(4)
      pars(5) = lamda !- pars1(1) * tmp
      pars(6) = pars1(5) 
      slit = pars(1) * EXP(- ((wave-pars(2))/pars(3))**2.0 ) + &
           pars(4) * EXP(- ((wave-pars(5))/pars(6))**4.0 )          
   ENDIF
   !WRITE(*, '(6D16.8)') pars
   !WRITE(*, '(15D12.4)') slit(1:nw) / MAXVAL(slit)
  
   RETURN

  END SUBROUTINE compute_slitprofile

! Use an average slit function for each fitting window instead of using variable slit widths
SUBROUTINE omislit_multi (wvlarr, specarr, specmod, npoints)

  ! ===============
  ! Input variables
  ! ===============
  INTEGER,                             INTENT (IN) :: npoints
  REAL (KIND=dp), DIMENSION (npoints), INTENT (IN) :: wvlarr, specarr
  
  ! ================
  ! Output variables
  ! ================
  REAL (KIND=dp), DIMENSION (npoints), INTENT (OUT) :: specmod
  
  ! ===============
  ! Local variables
  ! ===============
  INTEGER :: i, j, iw, ch, fidx, lidx, sidx, eidx, nslit, lpnt, fpnt, nhalf
  REAL (KIND=dp)                      :: delwvl, bndwav, slitsum
  REAL (KIND=dp), DIMENSION (npoints) :: locsli, locsli1, locsli2, locwvl
  INTEGER, DIMENSION (npoints)        :: idxs
    
  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=12), PARAMETER :: modulename = 'omislit_vary'
  
  fidx = 1
  DO iw = 1, numwin
     IF (iw < numwin) THEN
        bndwav = (winlim(iw, 2) + winlim(iw+1, 1) ) / 2.0
        lidx = MINVAL(MAXLOC(wvlarr(fidx:npoints), MASK=(wvlarr(fidx:npoints) < bndwav))) + fidx - 1
     ELSE
        lidx = npoints
     ENDIF
     IF (wvlarr(lidx) >= wvlarr(npoints) - 2.0) lidx = npoints
     IF (wvlarr(lidx)-wvlarr(fidx) <= 2.0) CYCLE

     ch = band_selectors(iw)   
     delwvl = wvlarr(fidx + 1) - wvlarr(fidx)
     IF (ch == 1) THEN
        nhalf = 1.0 / delwvl  ! 2 nm for UV-1 (truncate these < ~0.001)
     ELSE 
        nhalf = 0.5 / delwvl  ! 1 nm for UV-2 and VIS (truncate those < ~0.001)
     ENDIF
     nslit = nhalf * 2 + 1

     sidx = MAX(fidx, nhalf+1)
     eidx = MIN(lidx, npoints - nhalf)

     ! compute the average slit functions
     i = INT((fidx + lidx) / 2)

     !locwvl(1:nslit) = wvlarr(i-nhalf:i+nhalf) - wvlarr(i)     
     !CALL compute_slitprofile (ch, wvlarr(fidx), wvlarr(fidx)+locwvl(1:nslit), nslit, locsli1(1:nslit))
     !CALL compute_slitprofile (ch, wvlarr(lidx), wvlarr(lidx)+locwvl(1:nslit), nslit, locsli2(1:nslit))
     !locsli(1:nslit) = locsli1(1:nslit) + locsli2(1:nslit)

     CALL compute_slitprofile (ch, wvlarr(i), wvlarr(i-nhalf:i+nhalf), nslit, locsli(1:nslit))
     slitsum = SUM(locsli(1:nslit))

     ! safe convolution mode 
     DO i = sidx, eidx
        specmod(i) = DOT_PRODUCT(locsli(1:nslit), specarr(i-nhalf:i+nhalf)) / slitsum
     ENDDO

     ! left wrap mode
     DO i = fidx, sidx - 1
        idxs(1:nslit) = (/(j, j = i - nhalf, i+nhalf)/)
        lpnt = nhalf - i + 1
        idxs(1 : lpnt) = ABS(idxs(1 : lpnt)) + 2
        specmod(i) = DOT_PRODUCT(locsli(1:nslit), specarr(idxs(1:nslit))) / slitsum
     ENDDO

     ! right wrap mode
     DO i = eidx + 1, lidx
        idxs(1:nslit) = (/(j, j = i - nhalf, i+nhalf)/)
        fpnt = nhalf + (npoints - i) + 1
        idxs(fpnt:nslit) =  npoints * 2 - idxs(fpnt:nslit)
        specmod(i) = DOT_PRODUCT(locsli(1:nslit), specarr(idxs(1:nslit))) / slitsum
     ENDDO

     IF (lidx == npoints)   EXIT
     fidx = lidx + 1
  END DO

  RETURN
END SUBROUTINE omislit_multi

! Use variable slit function for during the convolution
SUBROUTINE omislit_vary (wvlarr, specarr, specmod, npoints)

  ! ===============
  ! Input variables
  ! ===============
  INTEGER,                             INTENT (IN) :: npoints
  REAL (KIND=dp), DIMENSION (npoints), INTENT (IN) :: wvlarr, specarr
  
  ! ================
  ! Output variables
  ! ================
  REAL (KIND=dp), DIMENSION (npoints), INTENT (OUT) :: specmod
  
  ! ===============
  ! Local variables
  ! ===============
  INTEGER  :: i, j, iw, ch, fidx, lidx, sidx, eidx, nslit, lpnt, fpnt, nhalf
  REAL (KIND=dp)                      :: delwvl, bndwav
  REAL (KIND=dp), DIMENSION (npoints) :: locsli
  INTEGER, DIMENSION (npoints)        :: idxs
    
  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=12), PARAMETER :: modulename = 'omislit_vary'
  
  fidx = 1
  DO iw = 1, numwin
     IF (iw < numwin) THEN
        bndwav = (winlim(iw, 2) + winlim(iw+1, 1) ) / 2.0
        lidx = MINVAL(MAXLOC(wvlarr(fidx:npoints), MASK=(wvlarr(fidx:npoints) <= bndwav))) + fidx - 1
     ELSE
        lidx = npoints
     ENDIF
     IF (wvlarr(lidx) >= wvlarr(npoints) - 2.0) lidx = npoints
     IF (wvlarr(lidx)-wvlarr(fidx) <= 2.0) CYCLE

     ch = band_selectors(iw)
     
     delwvl = wvlarr(fidx + 1) - wvlarr(fidx)
     IF (ch == 1) THEN
        nhalf = 1.0 / delwvl  ! 2 nm for UV-1 (truncate these < ~0.001)
     ELSE 
        nhalf = 0.5 / delwvl  ! 1 nm for UV-2 and VIS (truncate those < ~0.001)
     ENDIF
     nslit = nhalf * 2 + 1

     sidx = MAX(fidx, nhalf + 1)
     eidx = MIN(lidx, npoints - nhalf)

     ! safe convolution mode 
     DO i = sidx, eidx
        CALL compute_slitprofile (ch, wvlarr(i), wvlarr(i-nhalf:i+nhalf), nslit, locsli(1:nslit))
        specmod(i) = DOT_PRODUCT(locsli(1:nslit), specarr(i-nhalf:i+nhalf)) / SUM(locsli(1:nslit))
     ENDDO

     ! left wrap mode
     DO i = fidx, sidx - 1
        idxs(1:nslit) = (/(j, j = i - nhalf, i+nhalf)/)
        lpnt = nhalf - i + 1
        idxs(1 : lpnt) = ABS(idxs(1 : lpnt)) + 2
        CALL compute_slitprofile (ch, wvlarr(i), wvlarr(idxs(1:nslit)), nslit, locsli(1:nslit))
        specmod(i) = DOT_PRODUCT(locsli(1:nslit), specarr(idxs(1:nslit))) / SUM(locsli(1:nslit))
     ENDDO
     
     ! right wrap mode
     DO i = eidx + 1, lidx
        idxs(1:nslit) = (/(j, j = i - nhalf, i+nhalf)/)
        fpnt = nhalf + (npoints - i) + 1
        idxs(fpnt:nslit) =  npoints * 2 - idxs(fpnt:nslit)
        CALL compute_slitprofile (ch, wvlarr(i), wvlarr(idxs(1:nslit)), nslit, locsli(1:nslit))
        specmod(i) = DOT_PRODUCT(locsli(1:nslit), specarr(idxs(1:nslit))) / SUM(locsli(1:nslit))
     ENDDO
     
     IF (lidx == npoints)       EXIT
     fidx = lidx + 1
  END DO

  RETURN
END SUBROUTINE omislit_vary

!xliu, 10/29/2009, add the capability to convole multiple spectrum simultaneously
SUBROUTINE omislit_f2c (fwave, fspec, nf, nspec, cwave, cspec, nc)

  ! ===============
  ! Input variables
  ! ===============
  INTEGER,                        INTENT (IN)         :: nc, nf, nspec
  REAL (KIND=dp), DIMENSION (nf), INTENT (IN)         :: fwave
  REAL (KIND=dp), DIMENSION (nf, nspec), INTENT (IN)  :: fspec
  REAL (KIND=dp), DIMENSION (nc), INTENT (IN)         :: cwave
  REAL (KIND=dp), DIMENSION (nc, nspec), INTENT (OUT) :: cspec

  ! ===============
  ! Local variables
  ! ===============
  INTEGER        :: i, j, iwin, iw, fidx, fidxc, lidx, lidxc, midx, sidx, eidx, ch, imid
  REAL (KIND=dp) :: temp, dhalf, ssum
  REAL (KIND=dp), DIMENSION (nf) :: slit, tmpwav
  
  !slit(:) = 0.0; cspec(:, :)= 0.0
  fidx = 1; fidxc = 1
  DO iwin = 1, numwin

     IF (iwin == numwin) THEN
        lidx = nf; lidxc = nc
     ELSE
        temp = (winlim(iwin, 2) + winlim(iwin + 1, 1)) / 2.0
        lidx =  MINVAL(MAXLOC(fwave, MASK=(fwave <= temp)))
        lidxc = MINVAL(MAXLOC(cwave, MASK=(cwave <= temp)))
     ENDIF    

     ch = band_selectors(iwin)
     IF (ch == 1) THEN
        dhalf = 2.0
     ELSE
        dhalf = 1.0
     ENDIF

     ! compute average slit functions
     imid = INT((fidxc + lidxc) / 2)

     DO i = fidxc, lidxc
        ! Find the closest pixel
        midx = MINVAL(MAXLOC(fwave(fidx:lidx), MASK=(fwave(fidx:lidx) <= cwave(i)))) + fidx

        DO j = midx, fidx, -1
           IF ( ABS(fwave(j)-cwave(i)) >= dhalf .OR. j == 1) EXIT
        ENDDO
        sidx = j

        DO j = midx, lidx
           IF ( ABS(fwave(j)-cwave(i)) >= dhalf .OR. j == 1) EXIT
        ENDDO
        eidx = j
        tmpwav(sidx:eidx) = fwave(sidx:eidx) + (cwave(imid) - cwave(i))
        CALL compute_slitprofile (ch, cwave(imid), tmpwav(sidx:eidx), eidx-sidx+1, slit(sidx:eidx))
        
        ssum = SUM(slit(sidx:eidx))
        DO j = 1, nspec
           cspec(i, j) = SUM(fspec(sidx:eidx, j) * slit(sidx:eidx)) / ssum
        ENDDO
     ENDDO
     
     fidx = lidx + 1; fidxc = lidxc + 1
  ENDDO

  RETURN

END SUBROUTINE omislit_f2c

!xliu, 10/29/2009, add the capability to convole multiple spectrum simultaneously
SUBROUTINE omislit_vary_f2c (fwave, fspec, nf, nspec, cwave, cspec, nc)

  ! ===============
  ! Input variables
  ! ===============
  INTEGER,                        INTENT (IN)         :: nc, nf, nspec
  REAL (KIND=dp), DIMENSION (nf), INTENT (IN)         :: fwave
  REAL (KIND=dp), DIMENSION (nf, nspec), INTENT (IN)  :: fspec
  REAL (KIND=dp), DIMENSION (nc), INTENT (IN)         :: cwave
  REAL (KIND=dp), DIMENSION (nc, nspec), INTENT (OUT) :: cspec

  ! ===============
  ! Local variables
  ! ===============
  INTEGER        :: i, j, iwin, iw, fidx, fidxc, lidx, lidxc, midx, sidx, eidx, ch
  REAL (KIND=dp) :: temp, dhalf, ssum
  REAL (KIND=dp), DIMENSION (nf) :: slit
  
  !slit(:) = 0.0; cspec(:, :)= 0.0
  fidx = 1; fidxc = 1
  DO iwin = 1, numwin

     IF (iwin == numwin) THEN
        lidx = nf; lidxc = nc
     ELSE
        temp = (winlim(iwin, 2) + winlim(iwin + 1, 1)) / 2.0
        lidx =  MINVAL(MAXLOC(fwave, MASK=(fwave <= temp)))
        lidxc = MINVAL(MAXLOC(cwave, MASK=(cwave <= temp)))
     ENDIF    

     ch = band_selectors(iwin)
     IF (ch == 1) THEN
        dhalf = 2.0
     ELSE
        dhalf = 1.0
     ENDIF

     DO i = fidxc, lidxc
        ! Find the closest pixel
        midx = MINVAL(MAXLOC(fwave(fidx:lidx), MASK=(fwave(fidx:lidx) <= cwave(i)))) + fidx

        DO j = midx, fidx, -1
           IF ( ABS(fwave(j)-cwave(i)) >= dhalf .OR. j == 1) EXIT
        ENDDO
        sidx = j

        DO j = midx, lidx
           IF ( ABS(fwave(j)-cwave(i)) >= dhalf .OR. j == 1) EXIT
        ENDDO
        eidx = j
        CALL compute_slitprofile (ch, cwave(i), fwave(sidx:eidx), eidx-sidx+1, slit(sidx:eidx))
        
        ssum = SUM(slit(sidx:eidx))
        DO j = 1, nspec
           cspec(i, j) = SUM(fspec(sidx:eidx, j) * slit(sidx:eidx)) / ssum
        ENDDO
     ENDDO
     
     fidx = lidx + 1; fidxc = lidxc + 1
  ENDDO

  RETURN

END SUBROUTINE omislit_vary_f2c
  
END MODULE OMSAO_slitfunction_module
