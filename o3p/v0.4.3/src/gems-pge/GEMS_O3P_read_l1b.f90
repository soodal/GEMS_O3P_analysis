SUBROUTINE GEMS_O3P_read_l1b_irrad (nxcoadd, first_pix, last_pix, pge_error_status)
! OMI irradiance is uploaded from not daily measurements, but averaged irradiance data set

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: maxchlen, mrefl, max_ring_pts
  USE OMSAO_indices_module,    ONLY: wvl_idx, spc_idx, sig_idx, maxwin
  USE OMSAO_variables_module,  ONLY: refdbdir,band_selectors, numwin, winpix, winlim, scnwrt, &
                                     currpix, rm_mgline,wcal_bef_coadd
  USE ozprof_data_module,     ONLY: lun=>calunit, pos_alb, toms_fwhm, nrefl
  USE OMSAO_errstat_module
  USE GEMS_O3P_gemsdata_module, ONLY: mchannel, nwavel_max, nxtrack, nchannel,ncoadd, nxbin, gemsraddate, chs, &
                                      gems_irrad,gems_refl, gems_ring
  IMPLICIT NONE
  !-----------------------
  ! Input/Output varialbes
  !-----------------------
  INTEGER, INTENT (IN)  :: nxcoadd, first_pix, last_pix
  INTEGER, INTENT (OUT) :: pge_error_status
 
  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i4), PARAMETER            :: nbits = 16
  INTEGER (KIND=i2), DIMENSION(0:nbits-1) :: mflgbits
  INTEGER (KIND=i2), DIMENSION(nxcoadd, nwavel_max, 0:nbits-1) :: flgbits
  INTEGER (KIND=i2), DIMENSION(nwavel_max):: flgmsks
  INTEGER :: theyear, themon, theday, thedoy
  INTEGER :: i, j, is, ix, idx, iix, iw, irefl,ic
  INTEGER :: fidx, lidx, ch
  INTEGER :: nsolbin,nwavel,nbin, nx,  nsub, nring, noff1, noff2
 
  INTEGER (KIND=i4)             :: errstat
  INTEGER, DIMENSION (mchannel) :: spos, epos
  INTEGER, DIMENSION (mchannel) :: nwls
  INTEGER, DIMENSION (maxwin)   :: nwbin 
  REAL (KIND = dp)              :: normsc, wcenter
  REAL (KIND = dp), DIMENSION (maxwin, nxcoadd) :: wshis, wsqus  
  REAL (KIND = dp), DIMENSION (nxcoadd*2, sig_idx, nwavel_max) :: subspec
  REAL (KIND = dp), DIMENSION (spc_idx, max_ring_pts)          :: subrsol 
   
  LOGICAL :: error
  CHARACTER (LEN=maxchlen)                :: bkfname
  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=24), PARAMETER :: modulename = 'GEMS_O3P_read_l1b_irrad' 


  ! ----------------------------
  ! Initialize variables
  ! ---------------------------- 
  pge_error_status = pge_errstat_ok
  errstat  = pge_errstat_ok 
  gems_irrad%spec (1:nwavel_max,1:nxtrack) = 0.0
  gems_irrad%prec (1:nwavel_max,1:nxtrack) = 0.0
  gems_irrad%wavl (1:nwavel_max,1:nxtrack) = 0.0
  gems_irrad%qflg (1:nwavel_max,1:nxtrack) = 0
  gems_irrad%errstat(1:nxtrack) = pge_errstat_ok

  nsolbin = 1
  !IF (.NOT. use_backup) THEN
  ! DO is = 1, nchannel
  !     ch = chs(is)
  !     IF ch == 1. AND. nx == nfxtrack*2 ) nsolbin = 2 ???
  ! ENDDO
  !ENDIF

  IF (nsolbin == 2) THEN
    pge_error_status  = pge_errstat_error
    WRITE(*,*) ' Not process for zoom-in mode'
    return
  ENDIF

  !-----------------------------------------
  ! Determine sun-earth distance correction
  !-----------------------------------------
  READ(gemsraddate, '(I4,1X,2I2)') theyear, themon, theday
  CALL GET_DOY(theyear, themon,  theday, thedoy)
  IF (thedoy == 366) thedoy = 365
 
  OPEN (UNIT=lun, FILE= ADJUSTL(TRIM(refdbdir)) // 'solar-distance.dat', &
        STATUS='UNKNOWN', IOSTAT=errstat)
  IF ( errstat /= pge_errstat_ok ) THEN
      WRITE(*, '(2A)') modulename, ': Cannot open Sun-Earth Distance datafile!!!'
      pge_error_status = pge_errstat_error; RETURN
  END IF
  DO i = 1, 12
    READ(LUN, *)
  ENDDO
  DO i = 1, thedoy
    READ(LUN, *) normsc, normsc
  ENDDO
  CLOSE(LUN)
     
  normsc = 1.0 / normsc ** 2  ! solar energy is inversely proportional to square distance
    
  !-----------------------------------------
  ! load irradiance dataset to gems_irrad array
  !-----------------------------------------
  bkfname = ADJUSTL(TRIM(refdbdir)) // 'OMI/omisol_v003_avg_nshi_backup.dat'
  OPEN (UNIT=lun, FILE=TRIM(ADJUSTL(bkfname)), STATUS='UNKNOWN', IOSTAT=errstat)
  IF ( errstat /= pge_errstat_ok ) THEN
    WRITE(*, '(2A)') modulename, ': Cannot open solar backup file!!!'
    pge_error_status = pge_errstat_error; RETURN
  END IF

  nwavel = 0
  DO is = 1, 2
    READ(lun, *) nx, nwls(is)
    spos(is) = nwavel + 1; epos(is) = nwavel + nwls(is)
    DO i = 1, nx
      READ(lun, *) 
      DO j = 1, nwls(is)
        READ(lun, *) gems_irrad%wavl(nwavel + j, i), gems_irrad%spec(nwavel + j, i), &
                     gems_irrad%prec(nwavel + j, i), nsub, nsub
        IF (nsub > 0) gems_irrad%prec(nwavel + j, i) = gems_irrad%prec(nwavel + j, i) &
                     / SQRT( REAL(nsub, KIND=dp) )
      ENDDO
    ENDDO
    gems_irrad%spec(spos(is):epos(is), 1:nx) = gems_irrad%spec(spos(is):epos(is), 1:nx) * normsc
    gems_irrad%prec(spos(is):epos(is), 1:nx) = gems_irrad%prec(spos(is):epos(is), 1:nx) * normsc
    nwavel = epos(is)
     
    !print * , spos(is), epos(is), nwls(is), nwavel ! 1 159 159  159 ! 160 716 557 716
  ENDDO

  IF (nchannel == 1) THEN
    IF (band_selectors(1) == 1) THEN
       nwavel = nwls(1)         
    ELSE
      nwavel = nwls(2)
      gems_irrad%wavl(1:nwavel, :) = gems_irrad%wavl(spos(2):epos(2), :)
      gems_irrad%spec(1:nwavel, :) = gems_irrad%spec(spos(2):epos(2), :) 
      gems_irrad%prec(1:nwavel, :) = gems_irrad%prec(spos(2):epos(2), :)
      spos(2) = spos(2) - nwls(1); epos(2) = epos(2) - nwls(1) 
    ENDIF
  ENDIF
  CLOSE(LUN)


  ! dwavmax = ( omi_irradiance_wavl(2, 1) - omi_irradiance_wavl(1, 1) ) * 1.1
  ! IF (reduce_resolution) THEN
  ! ENDIF

  IF (nwavel > nwavel_max) THEN
    WRITE(*, '(A, 2i5)') "Need to increase nwavel_max!!!", nwavel, nwavel_max
    pge_error_status = pge_errstat_error; RETURN
  ENDIF
  ! Determine number of wavelengths to be read for deteriming cloud fraction
  fidx = MAXVAL ( MINLOC ( gems_irrad%wavl(1:nwavel, 1), MASK = &
       (gems_irrad%wavl(1:nwavel, 1) > pos_alb - toms_fwhm * 1.4) ))
  lidx = MAXVAL ( MAXLOC ( gems_irrad%wavl(1:nwavel, 1), MASK = &
       (gems_irrad%wavl(1:nwavel, 1) < pos_alb + toms_fwhm * 1.4) ))
  IF (fidx <1 .OR. lidx > nwavel) THEN
     WRITE(*, '(2A)') modulename, ': Need to change pos_alb/toms_fwhm!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF
  nrefl = lidx - fidx + 1 

  IF (nrefl > mrefl .or. nrefl == 1) THEN
    WRITE(*, '(2A)') modulename, ': check of nref, mrefl!!!'
    pge_error_status = pge_errstat_error; RETURN
  ENDIF

  ! Determine number of binning for different fitting windows
  DO iw = 1, numwin     
    ch = band_selectors(iw)
    nwbin(iw) = nxbin 
    IF (ch == 2 .and. ncoadd .ne. 1 ) THEN 
      nwbin(iw) = nxbin * ncoadd 
    ENDIF
    nwbin(iw) = nwbin(iw) * nsolbin
  ENDDO
  !-----------------------------------------
  ! subset and coadding
  !-----------------------------------------

  DO ix = first_pix, last_pix
    ! print *  , 'read irrad ix', ix
    currpix = ix
    ! Indices for UV-1 are from 1 to nfxtrack
    ! Indices for UV-2 are from 1 to nxtrack (nfxtrack * 2)
    ! UV2 indices corresponding to UV-1 pixel ix are ix * 2 -1 & ix * 2 (or iix+1, iix+2), respectively
    ! If additional across track coadding (e.g., nxbin) is performed, then for a particular ix
    ! UV1: nbin = nxbin; UV-2: nbin = nxbin * 2
    ! Coadded original across track pixels are: iix + 1 : iix + nbin (iix = (ix-1) * nbin)
    ! convert flags

    ! Get quality flag bits, coadd flags 
    flgmsks = 0 
    DO is = 1, nchannel 
      ch = chs(is)
      ! not use uwbin here because nwbin for different fitting window, the flags just for different channel
      IF (is == 1 ) THEN 
        nbin = nxbin
      ELSE 
        nbin = nxbin*ncoadd
      ENDIF
        nbin = nbin*nsolbin ! zoom mode UV1 solar irradiance got 60 positions
        iix = (ix-1)*nbin
        ! shift the position by 15
      IF (ch == 2 .AND. nsolbin == 2 ) THEN 
        ! iix = iix - (zoom_p1 - 1) * nsolbin
        pge_error_status  = pge_errstat_error
        WRITE (*,*) 'No process for zoom_p1!!! in read irradiance'
        return
      ENDIF

      IF (nbin  > 2) THEN
	    ! properly align cross track positions to be coadded (should be within one pixel)
        CALL prespec_align(nwls(ch), nbin, &
              gems_irrad%wavl(spos(ch):epos(ch), iix+1:iix+nbin), &
              gems_irrad%spec(spos(ch):epos(ch), iix+1:iix+nbin), &
              gems_irrad%prec(spos(ch):epos(ch), iix+1:iix+nbin), &       
              gems_irrad%qflg(spos(ch):epos(ch), iix+1:iix+nbin))	
	    !if (is == 2) print * , nbin, iix, gems_irrad%wavl(spos(ch):epos(ch), iix+nbin:iix+nbin)
      ENDIF
      
      DO ic = 1, nbin
        call convert_2bytes_to_16bits (nbits, nwls(ch), gems_irrad%qflg(spos(ch):epos(ch), iix+ic) , &
                                      flgbits(ic, spos(ch):epos(ch), 0:nbits-1))
             flgmsks(spos(ch):epos(ch)) = + flgbits(ic,spos(ch):epos(ch), 0)                &   ! Missing
                                          + flgbits(ic,spos(ch):epos(ch), 1)                &   ! Bad 
                                          + flgbits(ic,spos(ch):epos(ch), 2)                    ! Processing error
	
      ENDDO	
      ! More process for reduce_resolution
  ENDDO ! loop of nchannel
   
  ! Subset valid spectrum
  nsub =0; subspec = 0.0
  Do iw = 1, numwin
    ch   = band_selectors (iw)
    nbin = nwbin(iw);  iix = (ix - 1) * nbin
      
    winpix(iw, 1) = MINVAL ( MINLOC ( gems_irrad%wavl(spos(ch):epos(ch),iix + 1), &
              MASK = gems_irrad%wavl(spos(ch):epos(ch),iix + 1) >= winlim(iw, 1)) )
    winpix(iw, 2) = MAXVAL ( MAXLOC ( gems_irrad%wavl(spos(ch):epos(ch),iix + 1), &
             MASK = gems_irrad%wavl(spos(ch):epos(ch),iix + 1) <= winlim(iw, 2)) )


    gems_irrad%npix  (iw, ix) = nsub
    fidx = winpix(iw, 1) + spos(ch) - 1
    lidx = winpix(iw, 2) + spos(ch) - 1
    gems_irrad%winpix(iw,ix,1:2) = 0

    IF (rm_mgline .AND. winlim(iw, 1) < 286.0 .AND. winlim(iw, 2) > 286.0) THEN
      DO i = fidx, lidx
        IF (ALL(gems_irrad%spec(i, iix+1:iix+nbin)  > 0.0)   .AND. &
            ALL(gems_irrad%spec(i, iix+1:iix+nbin)  < 4.0E14) .AND. flgmsks(i) == 0 .AND. &
            (ALL(gems_irrad%wavl(i, iix+1:iix+nbin) < 278.8)   .OR. &
            ALL(gems_irrad%wavl(i, iix+1:iix+nbin)  > 281.0))  .AND. &
            (ALL(gems_irrad%wavl(i, iix+1:iix+nbin) < 284.7)   .OR. &
            ALL(gems_irrad%wavl(i, iix+1:iix+nbin)  > 285.7))) THEN
          nsub = nsub + 1
          subspec(1:nbin, wvl_idx, nsub) = gems_irrad%wavl(i, iix+1:iix+nbin)
          subspec(1:nbin, spc_idx, nsub) = gems_irrad%spec(i, iix+1:iix+nbin)
          subspec(1:nbin, sig_idx, nsub) = gems_irrad%prec(i, iix+1:iix+nbin)
        IF ( gems_irrad%winpix(iw, ix, 1) == 0)   gems_irrad%winpix(iw, ix, 1) = i
             gems_irrad%winpix(iw, ix, 2) = i
             gems_irrad%wind(nsub, ix) = i
        ENDIF
	    ENDDO
    ELSE
      DO i = fidx, lidx
        IF (ALL(gems_irrad%spec(i, iix+1:iix+nbin)  > 0.0)   .AND. &
            ALL(gems_irrad%spec(i, iix+1:iix+nbin)  < 4.0E14) .AND. flgmsks(i) == 0 ) THEN
          nsub = nsub + 1
          subspec(1:nbin, wvl_idx, nsub) = gems_irrad%wavl(i, iix+1:iix+nbin)
          subspec(1:nbin, spc_idx, nsub) = gems_irrad%spec(i, iix+1:iix+nbin)
          subspec(1:nbin, sig_idx, nsub) = gems_irrad%prec(i, iix+1:iix+nbin)
          IF ( gems_irrad%winpix(iw, ix, 1) == 0)   gems_irrad%winpix(iw, ix, 1) = i
                 gems_irrad%winpix(iw, ix, 2) = i
                 gems_irrad%wind(nsub, ix) = i
        ENDIF          
      ENDDO    
    ENDIF
    gems_irrad%npix(iw, ix) = nsub - gems_irrad%npix(iw, ix)
  ENDDO ! loop of numwin
  gems_irrad%nwav(ix) = nsub

  ! Coadding spectrum when necessary
  fidx = 1
  DO iw = 1, numwin       
    ch = band_selectors(iw)
    nbin = nwbin(iw)
    lidx = fidx + gems_irrad%npix(iw, ix) - 1 
    IF (nbin > 1) THEN ! wcal_bef_coadd F
      !WRITE(*,'(A)') ' Perform solwavcal_coadd'
      !print * , subspec(1:nbin,1, fidx:fidx), sum(subspec(1:nbin,1, fidx:fidx))/nbin
      CALL solwavcal_coadd(wcal_bef_coadd, gems_irrad%npix(iw, ix), nbin, &
           subspec(1:nbin, :, fidx:lidx), wshis(iw, 1:nbin), wsqus(iw, 1:nbin), error)
      ! print * , subspec(1:nbin,1, fidx:fidx)
	      
      IF (error) THEN
        WRITE(*, '(A)') 'No solar wavelength calibration before coadding!!!'
        gems_irrad%errstat(ix) = pge_errstat_warning
      ENDIF
    ENDIF
    fidx = lidx + 1    
  ENDDO
 
  ! GET solar spectrum for ring effect
  subrsol = 0.0

  !ring (a) iw=1 ring(13:126) = 269.13-308.817 Same As SubSPectra
  ch = band_selectors(1)
  IF (ch == 1) THEN 
    noff1 = 12
  ELSE
    noff1 = 25
  ENDIF

  nring = gems_irrad%npix(1,ix) + noff1
  subrsol(1:spc_idx, noff1+1: nring) = subspec(1, 1:spc_idx, 1:gems_irrad%npix(1,ix))
  gems_ring%sol_ndiv(ix) = 0

  ! ring (b) iw=1, ring(1:12) = 264.825-268.78
  ! add extra spectra before first window (uncoadded)
  ! if unavailable, needed to ammened with solar reference spectrum
  noff1 = noff1 + 1 ;     nbin = nwbin(1); iix = (ix - 1) * nbin

  DO i = winpix(1, 1) - 1, 1, -1
    IF (ALL(gems_irrad%spec(i, iix+1:iix+nbin) > 0.0) .AND. &
        ALL(gems_irrad%spec(i, iix+1:iix+nbin) < 4.0E14) .AND. flgmsks(i) == 0) THEN
        noff1 = noff1 - 1
        subrsol(wvl_idx, noff1) = SUM(gems_irrad%wavl(i, iix+1:iix+nbin)) / nbin
        subrsol(spc_idx, noff1) = SUM(gems_irrad%spec(i, iix+1:iix+nbin)) / nbin
        	
        IF (noff1 == 1) EXIT
    ENDIF
  ENDDO

  DO iw = 2, numwin
    ch = band_selectors(iw)
    nbin = nwbin(iw - 1) ; iix = (ix - 1) * nbin           

     
    IF (ch == band_selectors(iw - 1) ) THEN   
      DO i = winpix(iw-1, 2)+1, winpix(iw, 1)-1 
        IF (ALL(gems_irrad%spec(i, iix+1:iix+nbin) > 0.0) .AND. &
            ALL(gems_irrad%spec(i, iix+1:iix+nbin) < 4.0E14) .AND. flgmsks(i) == 0) THEN
            nring = nring + 1

            subrsol(wvl_idx, nring) = SUM(gems_irrad%wavl(i, iix+1:iix+nbin)) / nbin
            subrsol(spc_idx, nring) = SUM(gems_irrad%spec(i, iix+1:iix+nbin)) / nbin
                 
        ENDIF
      ENDDO
    ELSE  ! first channel 1 and second channel 2
      wcenter = (winlim(iw-1, 2) + winlim(iw, 1)) / 2.0 
      idx = MAXVAL ( MAXLOC ( gems_irrad%wavl(spos(ch-1):epos(ch-1), iix+1), &
      MASK =gems_irrad%wavl(spos(ch-1):epos(ch-1), iix+1) < wcenter ) ) + spos(ch-1) - 1

!(c) ring(127:131) : 309.114~310.29              
      DO i = winpix(iw-1, 2)+1, idx 
        IF (ALL(gems_irrad%spec(i, iix+1:iix+nbin) > 0.0) .AND. &
            ALL(gems_irrad%spec(i, iix+1:iix+nbin) < 4.0E14) .AND. flgmsks(i) == 0 .AND. &
            ALL(gems_irrad%wavl(i, iix+1:iix+nbin) > subrsol(wvl_idx, nring)) ) THEN
            nring = nring + 1

            subrsol(wvl_idx, nring) = SUM(gems_irrad%wavl(i, iix+1:iix+nbin)) / nbin
            subrsol(spc_idx, nring) = SUM(gems_irrad%spec(i, iix+1:iix+nbin)) / nbin                   
!print * , i, nring, subrsol(1, nring)
        ENDIF
      ENDDO
      !(c) ring(132:141) : 310.58-312.01
      gems_ring%sol_ndiv(ix) = nring               ! 1:nring is from the same channel
      nbin = nwbin(iw) ; iix = (ix - 1) * nbin
              
      idx = MAXVAL ( MINLOC ( gems_irrad%wavl(spos(ch):epos(ch), iix+1),   &
      MASK = gems_irrad%wavl(spos(ch):epos(ch), iix+1) > wcenter ) ) + spos(ch) - 1

      DO i = idx, winpix(iw, 1) - 2 + spos(ch)                 
        IF (ALL(gems_irrad%spec(i, iix+1:iix+nbin) > 0.0) .AND. &
            ALL(gems_irrad%spec(i, iix+1:iix+nbin) < 4.0E14) .AND. flgmsks(i) == 0 .AND. &
            ALL(gems_irrad%wavl(i, iix+1:iix+nbin) > subrsol(wvl_idx, nring)) ) THEN
          nring = nring + 1

          subrsol(wvl_idx, nring) = SUM(gems_irrad%wavl(i, iix+1:iix+nbin)) / nbin
          subrsol(spc_idx, nring) = SUM(gems_irrad%spec(i, iix+1:iix+nbin)) / nbin                   
		   ! print * , i, nring, subrsol(1, nring), nbin
        ENDIF
      ENDDO
    ENDIF
           
    idx = SUM(gems_irrad%npix(1:iw-1, ix))
    subrsol(wvl_idx:spc_idx, nring+1:nring+gems_irrad%npix(iw, ix)) = subspec(1, wvl_idx:spc_idx, &
            idx+1:idx+gems_irrad%npix(iw, ix))
    nring = nring + gems_irrad%npix(iw, ix)
  ENDDO

  ! Add extra spectra after fitting window
  noff2 = nring
  IF (ch == 2) THEN
    nring = nring + 25
  ELSE
    nring = nring + 12
  ENDIF
        
  DO i = spos(ch) + winpix(numwin, 2), epos(ch)
    IF (ch == 1 ) THEN
      nbin = nxbin
    ELSE
      nbin = nxbin * ncoadd
    ENDIF
      iix = (ix - 1) * nbin
           
    IF (ALL(gems_irrad%spec(i, iix+1:iix+nbin) > 0.0) .AND. &
        ALL(gems_irrad%spec(i, iix+1:iix+nbin) < 4.0E14) .AND. flgmsks(i) == 0) THEN
      noff2 = noff2 + 1
    
      subrsol(wvl_idx, noff2) = SUM(gems_irrad%wavl(i, iix+1:iix+nbin)) / nbin
      subrsol(spc_idx, noff2) = SUM(gems_irrad%spec(i, iix+1:iix+nbin)) / nbin
      ! write(*,'(3i5,3f10.4)') i, noff2,iix, subrsol(wvl_idx, noff2),gems_irrad%wavl(i, iix+1:iix+nbin)
    ENDIF
  
    IF (noff2 == nring) EXIT
  ENDDO     
  gems_ring%nsol(ix) = nring; gems_ring%sol_lin(ix) = noff1; gems_ring%sol_uin(ix) = noff2
   
  ! Get data for surface albedo & cloud fraction at 370.2 nm +/- 15 pixels
  irefl = 0; gems_refl%solwinpix(ix, 1:2) = 0
  nbin = nxbin * ncoadd
  iix = (ix - 1) * nbin
     
  idx = MAXVAL ( MINLOC ( gems_irrad%wavl(1:nwavel, iix+1), &
                          MASK = (gems_irrad%wavl(1:nwavel, iix+1) > pos_alb - toms_fwhm * 1.4) ))
  
  DO i  = idx, nwavel
    IF (ALL(gems_irrad%spec(i, iix+1:iix+nbin) > 0.0) .AND. &
        ALL(gems_irrad%spec(i, iix+1:iix+nbin) < 4.0E14) .AND. flgmsks(i) == 0) THEN
      irefl = irefl + 1
      gems_refl%solwavl(irefl, ix) = SUM(gems_irrad%wavl(i, iix+1:iix+nbin)) / nbin
      gems_refl%solspec(irefl, ix) = SUM(gems_irrad%spec(i, iix+1:iix+nbin)) / nbin
      !write(*,'(3i5,3f10.3)') irefl,nbin, iix+1,gems_refl%solwavl(irefl, ix), gems_irrad%wavl(i, iix+1:iix+nbin)
      IF ( gems_refl%solwinpix(ix, 1) == 0)  gems_refl%solwinpix(ix, 1) = i
        gems_refl%solwinpix(ix, 2) = i
      ENDIF
      IF (irefl == nrefl) EXIT
  ENDDO
  IF (irefl /= nrefl) THEN
    WRITE(*, *) 'Could not get enough irradiance points for cloud fraction!!!'
    gems_irrad%errstat(ix) = pge_errstat_error; CYCLE
  ENDIF
      
  IF (scnwrt) THEN
    WRITE(*, *) 'End Of Reading Irradiance Spectrum: ', ix
    DO i = 1, numwin
      WRITE(*,'(A10,I4,2f8.3,I4)') 'win = ', i, winlim(i,1), winlim(i,2),gems_irrad%npix(i, ix)
      IF (gems_irrad%npix(i, ix) < 4) THEN
        WRITE(*, '(A,f8.3,A3,f8.3)') ' Not enough points (>=4)  in window: ', winlim(i,1), ' - ', winlim(i,2)
        pge_error_status = pge_errstat_error
      ENDIF
    ENDDO
  ENDIF
  
  gems_irrad%norm(ix) = SUM ( subspec(1, spc_idx, 1:nsub) ) / nsub
     
  IF (gems_irrad%norm(ix) <= 0.0 ) THEN 
    gems_irrad%errstat(ix) = pge_errstat_error; CYCLE
  ENDIF
  
  gems_irrad%wavl(1:nsub, ix)  = subspec(1, wvl_idx, 1:nsub)
  gems_irrad%spec(1:nsub, ix)  = subspec(1, spc_idx, 1:nsub) / gems_irrad%norm(ix)
  gems_irrad%prec(1:nsub, ix)  = subspec(1, sig_idx, 1:nsub) / gems_irrad%norm(ix)    
  
  gems_ring%solwavl(1:nring, ix) = subrsol(1, 1:nring)
  gems_ring%solspec(1:nring, ix) = subrsol(2, 1:nring) / gems_irrad%norm(ix) 
  !print * , gems_irrad%norm(ix)
  !print * , gems_irrad%wavl(1:5, ix)
  !print * , gems_irrad%spec(1:5, ix)
  !print * , gems_irrad%prec(1:5, ix)
 
  ENDDO ! xtrack
RETURN
END SUBROUTINE GEMS_O3P_read_l1b_irrad


SUBROUTINE GEMS_O3P_read_l1b_rad (nxcoadd, first_pix, last_pix,ny, offline, pge_error_status)
! OMI irradiance is uploaded from not daily measurements, but averaged irradiance data set

  USE OMSAO_precision_module
  USE OMSAO_indices_module,    ONLY: wvl_idx, spc_idx, sig_idx, maxwin
  USE OMSAO_variables_module,  ONLY: refdbdir,band_selectors, numwin, winpix, winlim, scnwrt, &
                                     currpix, rm_mgline,wcal_bef_coadd, szamax
  USE OMSAO_parameters_module, ONLY: maxchlen, mrefl
  USE OMSAO_errstat_module
  USE GEMS_O3P_gemsdata_module, ONLY: nchannel, mchannel, nwavel_max, nxtrack_max, &
                                      nxtrack, nfxtrack, ntimes, chs, ncoadd, nxbin, nybin, &
                                     gemsl1btype, gems_uv1, gems_uv2,  &
                                     gems_rad, gems_irrad,  &
                                     gems_mqflg, gems_saa,  gems_xflag, &
                                     gems_sza, gems_refl, gems_aza, gems_vza
  USE ozprof_data_module,       ONLY:  nrefl
  IMPLICIT NONE	
  !-----------------------
  ! INput/Output varialbes
  !-----------------------
  INTEGER, INTENT (IN)  :: nxcoadd, first_pix, last_pix, ny, offline
  INTEGER, INTENT (OUT) :: pge_error_status
  !-----------------------
  ! Local varialbes
  !-----------------------
  INTEGER, PARAMETER    :: nbits = 16
  INTEGER (KIND=i2), DIMENSION(0:nbits-1)   :: mflgbits , tmp_mflgbits 
  INTEGER (KIND=i2), DIMENSION(nxcoadd, nwavel_max, 0:nbits-1) :: flgbits
  INTEGER (KIND=i2), DIMENSION(nwavel_max)                     :: flgmsks
  REAL   (KIND= dp)                         :: tmpNinteg
  INTEGER, DIMENSION (nwavel_max)           :: idxs
  INTEGER   (KIND=4), DIMENSION (nchannel)  :: nwls
  INTEGER   (KIND=i2), DIMENSION(nchannel)  :: spos, epos
  INTEGER, DIMENSION (maxwin)   :: nwbin 
  REAL (KIND = dp), DIMENSION (nxcoadd, sig_idx, nwavel_max) :: subspec
  LOGICAL, DIMENSION (maxwin, nxtrack_max)  :: wavcals 
  REAL (KIND = dp), DIMENSION (maxwin, nxtrack_max, nxcoadd) :: wshis, wsqus 
  INTEGER :: i,j,k, is, ix,iix, iy,iiy,iw, ch, nwavel,nwl, nx, nbin, ic, nsub, fidx, lidx, ii, error, irefl
  LOGICAL :: correct_merr = .TRUE.
  TYPE (gemsl1btype) :: gems_uv
  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=23), PARAMETER :: modulename = 'omi_read_radiance_lines'

  
  pge_error_status = pge_errstat_ok 


  !--------------------------------
  ! Coadding Y-direction
  !--------------------------------
  gems_rad%wavl(:,:,:) = 0.0
  gems_rad%spec(:,:,:) = 0.0
  gems_rad%prec(:,:,:) = 0.0
  gems_rad%qflg(:,:,:) = pge_errstat_ok
  gems_rad%errstat(1:nxtrack, 1:ny) = pge_errstat_ok
  gems_rad%line_errstat( 1:ny) = pge_errstat_ok
  gems_rad%nwav(:,:) = 0
 
  j = 1  ; nwavel = 0 ;  wavcals = .TRUE.

  DO is = 1, nchannel

    ch = chs(is)
    IF (ch == 1) THEN
	    gems_uv    = gems_uv1

    ELSE IF (ch == 2) THEN
      gems_uv    = gems_uv2

    ENDIF
    nwls(ch) = gems_uv%nwavel
    nx       = gems_uv%nxtrack
    DO iy = 1, ny

      mflgbits = 0
      iiy = (iy-1)*nybin + offline 
      Do i = 1,  nybin   
        iiy =  iiy +1

        CALL convert_2bytes_to_16bits ( nbits, 1, gems_uv%mqflg(iiy), tmp_mflgbits(0:nbits-1))
        mflgbits = mflgbits + tmp_mflgbits(0:nbits-1)

        IF (correct_merr) THEN
          tmpNinteg = 2.d0 / gems_uv%ExposureTime(iiy)
        ENDIF      

        nwl = j + nwls(ch) - 1
            
        gems_rad%spec(j:nwl, 1:nx, iy) = gems_rad%spec(j:nwl, 1:nx, iy)+ gems_uv%spec(1:nwls(ch),1:nx, iiy)
        gems_rad%wavl(j:nwl, 1:nx, iy) = gems_rad%wavl(j:nwl, 1:nx, iy)+ gems_uv%wavl(1:nwls(ch),1:nx, iiy)
        ! if (is == 2) write(*,'(1i5, 11f10.5)') iiy, gems_uv%wavl(1:11,14, iiy)

        IF (correct_merr) THEN
          gems_rad%prec(j:nwl, 1:nx, iy) = gems_rad%prec(j:nwl, 1:nx, iy) + &
                                           gems_uv%prec(1:nwls(ch),1:nx, iiy) / SQRT( tmpNinteg ) 
        ELSE
          gems_rad%prec(j:nwl, 1:nx, iy) = gems_rad%prec(j:nwl, 1:nx, iy) + &
                                           gems_uv%prec(1:nwls(ch),1:nx, iiy)
        ENDIF

        DO ix = 1, nxtrack
          CALL coadd_2bytes_qflgs(nbits, nwls(ch), gems_rad%qflg(j:nwl, ix,iy), gems_uv%qflg(1:nwls(ch), ix,iiy)) 
        ENDDO
	   
      ENDDO ! End of Y-Coadding
            
      IF (mflgbits(0) >= 1 .OR. mflgbits(1) >= 1 .OR. mflgbits(3) >= 1 .OR. mflgbits(12) >= 1) THEN
        WRITE(*, *) 'All radiances could not be used: line ', iiy, ' Swath ', is
        gems_rad%line_errstat(iy) = pge_errstat_error      
      ELSE IF (ANY(mflgbits >= 1)) THEN
        !WRITE(*, *) 'Warning set on all radiances: line',  blockline, ' Swath ', is
        IF ( gems_rad%line_errstat(iy) /= pge_errstat_error) &
             gems_rad%line_errstat(iy) = pge_errstat_warning      
        ! Over SAA region
        IF (mflgbits(10) >= 1) gems_saa(iy) = 1
      ENDIF

	
      gems_rad%spec(j:nwl, 1:nx, iy) = gems_rad%spec(j:nwl, 1:nx,iy) / nybin
      gems_rad%wavl(j:nwl, 1:nx, iy) = gems_rad%wavl(j:nwl, 1:nx,iy) / nybin     
     IF ( correct_merr ) then !jbak
        gems_rad%prec(j:nwl, 1:nx, iy) = gems_rad%prec(j:nwl, 1:nx,iy) / nybin / SQRT(1.0D0 * nybin)
      ELSE
        gems_rad%prec(j:nwl, 1:nx, iy) = gems_rad%prec(j:nwl, 1:nx,iy) / nybin
      ENDIF
    ENDDO ! end iloop

    ! Sort data in wavelength increasing order   
    nwavel = nwavel + nwls(ch); spos(ch) = j; j = nwavel + 1; epos(ch) = nwavel

    IF (gems_rad%wavl(spos(ch), 1, 1) > gems_rad%wavl(epos(ch), 1, 1)) THEN
      idxs(spos(ch):epos(ch)) = (/ (i, i = epos(ch), spos(ch), -1) /)
      gems_rad%wavl(spos(ch):epos(ch), 1:nx, 1:ny) = gems_rad%wavl(idxs(spos(ch):epos(ch)), 1:nx, 1:ny)
      gems_rad%spec(spos(ch):epos(ch), 1:nx, 1:ny) = gems_rad%spec(idxs(spos(ch):epos(ch)), 1:nx, 1:ny)
      gems_rad%prec(spos(ch):epos(ch), 1:nx, 1:ny) = gems_rad%prec(idxs(spos(ch):epos(ch)), 1:nx, 1:ny)
      gems_rad%qflg(spos(ch):epos(ch), 1:nx, 1:ny) = gems_rad%qflg(idxs(spos(ch):epos(ch)), 1:nx, 1:ny)     
    ENDIF 

  ENDDO ! end  channel loop


! Process for reduced resoution (not included for gems)

  IF (nwavel > nwavel_max) THEN
    WRITE(*, *) "Need to increase nwavel_max!!!"
    pge_error_status = pge_errstat_error; RETURN
  ENDIF

  !--------------------------------
  ! SubSet & Coadding X-direction
  !--------------------------------
  ! Determine number of binning for different fitting windows
  DO iw = 1, numwin
    ch = band_selectors(iw)
    IF (ch == 1 ) THEN
      nwbin(iw) = nxbin
    ELSE IF (ch == 2) THEN 
      nwbin(iw) = nxbin * ncoadd
    ENDIF
  ENDDO

  DO iy = 1, ny 
    IF (gems_rad%line_errstat(iy) == pge_errstat_error) THEN
      gems_rad%errstat(first_pix:last_pix, iy) = pge_errstat_error; CYCLE
    ENDIF
    DO ix = first_pix, last_pix
      currpix = ix
      IF (gems_sza(ix, iy) > szamax .OR. gems_sza(ix, iy) < 0 ) THEN
        gems_rad%errstat(ix, iy) = pge_errstat_error; CYCLE
      ENDIF

      IF (gems_xflag%rowanomaly(ix, iy) == 1) THEN
        print * , 'rowanomaly_flg'
        gems_rad%errstat(ix, iy) = pge_errstat_error; CYCLE
      ENDIF

      ! Get quality flag bits
      ! Coadd uv-2 flags if necessary to avoid coadding inconsistent # of pixels 
      flgmsks = 0 
      DO is = 1, nchannel
        ch = chs(is)
        IF (ch == 1) THEN 
          nbin = nxbin
        ELSE IF (ch == 2) THEN
          nbin = nxbin*ncoadd           
        ENDIF
        iix = (ix - 1) * nbin 
      
        ! properly align cross track positions to be coadded (should be within one pixel)

        IF (nbin > 2) CALL prespec_align(nwls(ch), nbin,                                       &
                                         gems_rad%wavl(spos(ch):epos(ch), iix+1:iix+nbin, iy), &
                                         gems_rad%spec(spos(ch):epos(ch), iix+1:iix+nbin, iy), &
                                         gems_rad%prec(spos(ch):epos(ch), iix+1:iix+nbin, iy), &       
                                         gems_rad%qflg(spos(ch):epos(ch), iix+1:iix+nbin, iy))
              
        DO ic = 1, nbin
          CALL convert_2bytes_to_16bits ( nbits, nwls(ch),gems_rad%qflg(spos(ch):epos(ch), &
                                          iix + ic, iy), flgbits(ic, spos(ch):epos(ch), 0:nbits-1))
                   flgmsks(spos(ch):epos(ch))                       &
                 = flgmsks(spos(ch):epos(ch))                       &
                 + flgbits(ic, spos(ch):epos(ch), 0)                &   ! Missing
                 + flgbits(ic, spos(ch):epos(ch), 1)                &   ! Bad 
                 + flgbits(ic, spos(ch):epos(ch), 2)                &   ! Processing error
                 + flgbits(ic, spos(ch):epos(ch), 4)                &   ! RTS_Pixel_Warning Flag
                 + flgbits(ic, spos(ch):epos(ch), 5)                &   ! Saturation Possibility Flag
                 + flgbits(ic, spos(ch):epos(ch), 7)                    ! Dark Current Warning Flag
        ENDDO
      ENDDO ! End of nchannel
           
      ! Subset valid data
      nsub = 0 ; subspec=0.0  ; fidx=1
         
      DO iw = 1, numwin
        ch = band_selectors(iw)
        gems_rad%npix(iw, ix, iy) = nsub
        nbin = nwbin(iw)
        iix = (ix - 1) * nbin
        !fidx = gems_irrad%winpix(iw, ix, 1) 
        !lidx = gems_irrad%winpix(iw, ix, 2) 
        lidx = fidx + gems_irrad%npix(iw, ix) - 1
        
        DO ii = fidx, lidx
          i = gems_irrad%wind(ii, ix)
          IF (ALL( gems_rad%spec(i, iix+1:iix+nbin, iy) > 0.0) .AND. &
              ALL( gems_rad%spec(i, iix+1:iix+nbin, iy) < 4.0E14) .AND. flgmsks(i) == 0 ) THEN
            nsub = nsub + 1
            subspec(1:nbin, wvl_idx, nsub) =  gems_rad%wavl(i, iix+1:iix+nbin, iy)
            subspec(1:nbin, spc_idx, nsub) =  gems_rad%spec(i, iix+1:iix+nbin, iy)
            subspec(1:nbin, sig_idx, nsub) =  gems_rad%prec(i, iix+1:iix+nbin, iy)
            gems_rad%wind(nsub, ix, iy) = ii
            !print * , i,ii,subspec(1, 1, nsub),subspec(1, sig_idx, nsub)
          ENDIF
        ENDDO
         
        fidx = lidx + 1
        gems_rad%npix(iw, ix, iy) = nsub - gems_rad%npix(iw, ix, iy)

! If the # of wavelengths is <= 75% of the # of irradiances, stop processing this pixel
        IF (gems_rad%npix(iw, ix, iy) <= gems_irrad%npix(iw, ix) * 0.8 ) THEN   !geun  0.9 -> 0.8
          WRITE(*, '(A,5I5,F9.2)') 'Too fewer radiance points: ', ix, iy, iw, &
          gems_rad%npix(iw, ix, iy), gems_irrad%npix(iw, ix), gems_sza(ix, iy)
          gems_rad%errstat(ix, iy) = pge_errstat_error; CYCLE  ! geun modified from EXIT to CYCLE
        ENDIF
      ENDDO ! End of numwin
     
      IF ( gems_rad%errstat(ix, iy) == pge_errstat_error) CYCLE   ! This pixel will not be processed.     
      gems_rad%nwav(ix, iy) = nsub

! Perform coadding if UV-2 is selected with UV-1
      fidx = 1     
      DO iw = 1, numwin
        ch = band_selectors(iw);  nbin = nwbin(iw)
        lidx = fidx + gems_rad%npix(iw, ix, iy) - 1 
        IF (nbin > 1) THEN
          CALL radwavcal_coadd(wcal_bef_coadd, wavcals(iw, ix), iw, ix, gems_rad%npix(iw, ix, iy), nbin, &
               subspec(1:nbin, :, fidx:lidx), wshis(iw, ix, 1:nbin), wsqus(iw, ix, 1:nbin), error)
               wavcals(iw, ix) = .FALSE.
          IF (error) THEN
            WRITE(*, '(A)') 'No radiance wavelength calibration before coadding!!!'
            pge_error_status = pge_errstat_warning
          ENDIF
        ENDIF
        fidx = lidx + 1
      ENDDO

! Get data for surface albedo & cloud fraction at 370.2 nm +/- 20 pixels
      irefl = 0; fidx = gems_refl%solwinpix(ix, 1)
      nbin = nxbin * ncoadd  
      iix = (ix - 1) * nbin
      DO i  = fidx, nwavel
        IF ( ALL(gems_rad%spec(i, iix+1:iix+nbin, iy) > 0.0) .AND. &
             ALL(gems_rad%spec(i, iix+1:iix+nbin, iy) < 4.0E14) .AND. flgmsks(i) == 0) THEN
        irefl = irefl + 1
        gems_refl%radwavl(irefl, ix, iy) = SUM(gems_rad%wavl(i, iix+1:iix+nbin, iy)) / nbin
        gems_refl%radspec(irefl, ix, iy) = SUM(gems_rad%spec(i, iix+1:iix+nbin, iy)) / nbin
        !write(*,'(3i5,3f10.3)')  irefl, ix, iix+1, gems_refl%radwavl(irefl, ix, iy) ,gems_rad%wavl(i, iix+1:iix+nbin, iy)
        ENDIF	
        IF (irefl == nrefl) EXIT
      ENDDO
      IF (irefl /= nrefl) THEN
        WRITE(*, '(A, 2I5, F9.2)') 'Number of rad/sol points (cloud fraction) do not match: ', &
        ix, iy, gems_sza(ix, iy)
        gems_rad%errstat(ix, iy) = pge_errstat_error  
      ENDIF

  
      gems_rad%norm(ix, iy) = SUM (subspec(1, spc_idx, 1:nsub) ) / nsub
      !gems_rad%norm(ix,iy) = 1.0E11 
      IF ( gems_rad%norm(ix, iy) <= 0.0 ) THEN 
        write(*,'(a)') 'gems_rad%normal(ix,iy) <0'
        print * , gems_rad%norm(ix, iy) 
        pge_error_status = pge_errstat_error; RETURN
      ENDIF
      gems_rad%wavl(1:nsub, ix, iy) = subspec(1, wvl_idx, 1:nsub)
      gems_rad%spec(1:nsub, ix, iy) = subspec(1, spc_idx, 1:nsub) / gems_rad%norm(ix, iy)
      gems_rad%prec(1:nsub, ix, iy) = subspec(1, sig_idx, 1:nsub) / gems_rad%norm(ix, iy)       
      !write(*,*) 'hello wasp !  ',nsub !gems_rad%spec(1:nsub,ix,iy)
    ENDDO ! END of x-direction
  ENDDO ! End of y-direction
!DO i=1,nsub
  !write(*,*) 'hello wasp !  ',gems_rad%wavl(i,15,245),gems_rad%spec(i,15,245)
!ENDDO
  !write(*,*) 'hello wasp !  ',gems_rad%spec(1:nsub,ix,iy)
!stop
END SUBROUTINE GEMS_O3P_read_l1b_rad
