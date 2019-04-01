SUBROUTINE SYNT_O3P_read_l1b_irrad (nxcoadd, first_pix, last_pix, pge_error_status)
! OMI irradiance is uploaded from not daily measurements, but averaged irradiance data set

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: maxchlen, mrefl, max_ring_pts
  USE OMSAO_indices_module,    ONLY: wvl_idx, spc_idx, sig_idx, maxwin
  USE OMSAO_variables_module,  ONLY: refdbdir,band_selectors, numwin, winpix, winlim, scnwrt, &
                                     currpix, rm_mgline,wcal_bef_coadd
  USE ozprof_data_module,     ONLY: lun=>calunit, pos_alb, toms_fwhm, nrefl
  USE OMSAO_errstat_module
  USE GEMS_O3P_gemsdata_module, ONLY: mchannel, nwavel_max, nchannel, ncoadd, nxbin, gemsraddate, chs, &
                                      gems_irrad,gems_refl, gems_ring
  USE SYNT_read_l1b,            ONLY: synt_irrad_spec, synt_irrad_wavl
  USE SYNT_data_module,         ONLY: nwavel, lower_spec, upper_spec, nxtrack

  IMPLICIT NONE
 !-----------------------
 ! INput/Output varialbes
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
  INTEGER :: nsolbin,nbin, nx,  nsub, nring, noff1, noff2

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
  CHARACTER (LEN=24), PARAMETER :: modulename = 'SYNT_O3P_read_l1b_irrad' 


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

  IF (nsolbin == 2) THEN
      pge_error_status  = pge_errstat_error
      WRITE(*,*) ' NOt process for zoom-in mode'
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
    
  ! one channel in SYNT
  spos(1) = 1; epos(1) = nwavel 


  ! load irradiance data 
  gems_irrad%wavl=synt_irrad_wavl
  gems_irrad%spec=synt_irrad_spec
  gems_irrad%prec=1.0   ! precision for test

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
  nwbin(1) = nxbin 
  nwbin(1) = nwbin(1) * nsolbin
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
  ! not use uwbin here because nwbin for different fitting window, the flags just for different channel
  nbin = nxbin
  nbin = nbin*nsolbin ! zoom mode UV1 solar irradiance got 60 positions
  iix = (ix-1)*nbin

  IF (nbin  > 2) THEN
  ! properly align cross track positions to be coadded (should be within one pixel)
    CALL prespec_align(nwavel, nbin, &
         gems_irrad%wavl(1:nwavel, iix+1:iix+nbin), &
         gems_irrad%spec(1:nwavel, iix+1:iix+nbin), &
         gems_irrad%prec(1:nwavel, iix+1:iix+nbin), &       
         gems_irrad%qflg(1:nwavel, iix+1:iix+nbin))	
   !if (is == 2) print * , nbin, iix, gems_irrad%wavl(1:nwavel, iix+nbin:iix+nbin)
  ENDIF

  DO ic = 1, nbin
    call convert_2bytes_to_16bits (nbits, nwavel, gems_irrad%qflg(1:nwavel, iix+ic) , &
      flgbits(ic, 1:nwavel, 0:nbits-1))
      flgmsks(1:nwavel) = + flgbits(ic,1:nwavel, 0)                &   ! Missing
                                   + flgbits(ic,1:nwavel, 1)                &   ! Bad 
                                   + flgbits(ic,1:nwavel, 2)                    ! Processing error
  ENDDO

  ! Subset valid spectrum
  nsub =0; subspec = 0.0
  nbin = nwbin(1);  iix = (ix - 1) * nbin
  
  winpix(1, 1) = MINVAL ( MINLOC ( gems_irrad%wavl(1:nwavel,iix + 1), &
         MASK = gems_irrad%wavl(1:nwavel,iix + 1) >= winlim(1, 1)) )
  winpix(1, 2) = MAXVAL ( MAXLOC ( gems_irrad%wavl(1:nwavel,iix + 1), &
         MASK = gems_irrad%wavl(1:nwavel,iix + 1) <= winlim(1, 2)) )

  gems_irrad%npix(1, ix) = nsub
  fidx = winpix(1, 1) 
  lidx = winpix(1, 2) 
  gems_irrad%winpix(1,ix,1:2) = 0

  IF (rm_mgline .AND. winlim(1, 1) < 286.0 .AND. winlim(1, 2) > 286.0) THEN
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
          IF ( gems_irrad%winpix(1, ix, 1) == 0)   gems_irrad%winpix(1, ix, 1) = i
          gems_irrad%winpix(1, ix, 2) = i
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
          IF ( gems_irrad%winpix(1, ix, 1) == 0)   gems_irrad%winpix(1, ix, 1) = i
          gems_irrad%winpix(1, ix, 2) = i
          gems_irrad%wind(nsub, ix) = i
      ENDIF          
    ENDDO    
  ENDIF
  gems_irrad%npix(1, ix) = nsub - gems_irrad%npix(1, ix)
  gems_irrad%nwav(ix) = nsub

  ! Coadding spectrum when necessary
  fidx = 1
  nbin = nwbin(1)
  lidx = fidx + gems_irrad%npix(1, ix) - 1 
  IF (nbin > 1) THEN ! wcal_bef_coadd F
    !WRITE(*,'(A)') ' Perform solwavcal_coadd'
     !print * , subspec(1:nbin,1, fidx:fidx), sum(subspec(1:nbin,1, fidx:fidx))/nbin
    CALL solwavcal_coadd(wcal_bef_coadd, gems_irrad%npix(1, ix), nbin, &
         subspec(1:nbin, :, fidx:lidx), wshis(1, 1:nbin), wsqus(1, 1:nbin), error)
    ! print * , subspec(1:nbin,1, fidx:fidx)
    
    IF (error) THEN
       WRITE(*, '(A)') 'No solar wavelength calibration before coadding!!!'
       gems_irrad%errstat(ix) = pge_errstat_warning
    ENDIF
  ENDIF
  fidx = lidx + 1    
 
 ! GET solar spectrum for ring effect
  subrsol = 0.0

  !ring (a) iw=1 ring(13:126) = 269.13-308.817 Same As SubSPectra
  noff1 = 12
  nring = gems_irrad%npix(1,ix) + noff1
  subrsol(1:spc_idx, noff1+1: nring) = subspec(1, 1:spc_idx, 1:gems_irrad%npix(1,ix))
  gems_ring%sol_ndiv(ix) = 0

  ! ring (b) iw=1, ring(1:12) = 264.825-268.78
  ! add extra spectra before first window (uncoadded)
  ! if unavailable, needed to ammened with solar reference spectrum
   noff1 = noff1 + 1 ;     nbin = nwbin(1); iix = (ix - 1) * nbin


   DO i = winpix(1, 1) - 1, 1, -1
      IF (ALL(gems_irrad%spec(i, iix+1:iix+nbin) > lower_spec) .AND. &
          ALL(gems_irrad%spec(i, iix+1:iix+nbin) < upper_spec) .AND. flgmsks(i) == 0) THEN
          noff1 = noff1 - 1
          subrsol(wvl_idx, noff1) = SUM(gems_irrad%wavl(i, iix+1:iix+nbin)) / nbin
          subrsol(spc_idx, noff1) = SUM(gems_irrad%spec(i, iix+1:iix+nbin)) / nbin
	
           IF (noff1 == 1) EXIT
       ENDIF
   ENDDO

   ! Add extra spectra after fitting window
   noff2 = nring
   nring = nring + 12

   DO i = 1 + winpix(1, 2), nwavel    ! wavelength range outside fitting window

       nbin = nxbin
       iix = (ix - 1) * nbin
           
       IF (ALL(gems_irrad%spec(i, iix+1:iix+nbin) > lower_spec) .AND. &
           ALL(gems_irrad%spec(i, iix+1:iix+nbin) < upper_spec) .AND. flgmsks(i) == 0) THEN
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
   
   idx = MAXVAL ( MINLOC ( gems_irrad%wavl(1:nwavel, iix+1), MASK = &
        (gems_irrad%wavl(1:nwavel, iix+1) > pos_alb - toms_fwhm * 1.4) ))

   DO i  = idx, nwavel
      IF (ALL(gems_irrad%spec(i, iix+1:iix+nbin) > lower_spec) .AND. &
           ALL(gems_irrad%spec(i, iix+1:iix+nbin) < upper_spec) .AND. flgmsks(i) == 0) THEN
         irefl = irefl + 1
         gems_refl%solwavl(irefl, ix) = SUM(gems_irrad%wavl(i, iix+1:iix+nbin)) / nbin
         gems_refl%solspec(irefl, ix) = SUM(gems_irrad%spec(i, iix+1:iix+nbin)) / nbin
       !  write(*,'(3i5,3f10.3)') irefl,nbin, iix+1,gems_refl%solwavl(irefl, ix), gems_irrad%wavl(i, iix+1:iix+nbin)
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
END SUBROUTINE SYNT_O3P_read_l1b_irrad


SUBROUTINE SYNT_O3P_read_l1b_rad (nxcoadd, first_pix, last_pix,ny, offline, pge_error_status)
! OMI irradiance is uploaded from not daily measurements, but averaged irradiance data set

  USE OMSAO_precision_module
  USE OMSAO_indices_module,    ONLY: wvl_idx, spc_idx, sig_idx, maxwin
  USE OMSAO_variables_module,  ONLY: refdbdir,band_selectors, numwin, winpix, winlim, scnwrt, &
                                     currpix, rm_mgline,wcal_bef_coadd, szamax
  USE OMSAO_parameters_module, ONLY: maxchlen, mrefl
  USE OMSAO_errstat_module
  USE GEMS_O3P_gemsdata_module, ONLY: mchannel, nwavel_max, nchannel, nxtrack_max,nfxtrack, ntimes, ncoadd, nxbin, nybin, &
                                     gemsl1btype, gems_uv1, gems_uv2,  &
                                     gems_rad,gems_irrad, &
                                     gems_mqflg,gems_saa,  gems_xflag, &
                                     gems_sza, gems_refl
  USE ozprof_data_module,       ONLY:  nrefl
  !USE SYNT_read_l1b,            ONLY: synt_rad_spec, synt_rad_wavl
  USE SYNT_data_module,         ONLY: nwavel, synt_rad, rad_group, lower_spec, upper_spec,synt_rad_spec, synt_rad_wavl, &
                                      synt_rad_snr, nxtrack


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
  INTEGER   (KIND=i2), DIMENSION(nchannel)  :: spos, epos
  INTEGER, DIMENSION (maxwin)   :: nwbin 
  REAL (KIND = dp), DIMENSION (nxcoadd, sig_idx, nwavel_max) :: subspec
  LOGICAL, DIMENSION (maxwin, nxtrack_max)  :: wavcals 
  REAL (KIND = dp), DIMENSION (maxwin, nxtrack_max, nxcoadd) :: wshis, wsqus 
  INTEGER :: i,j,k, is, ix,iix, iy,iiy,iw, nx, nbin, ic, nsub, fidx, lidx, ii, error, irefl
  LOGICAL :: correct_merr = .TRUE.
  TYPE (rad_group) :: gems_uv
  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=23), PARAMETER :: modulename = 'omi_read_radiance_lines'

  pge_error_status = pge_errstat_ok 

  synt_rad%spec = synt_rad_spec
  synt_rad%wavl = synt_rad_wavl
  !synt_rad%prec = 1.0
  !synt_rad%prec = synt_rad_spec*0.01
  synt_rad%prec = synt_rad_spec*(1./synt_rad_snr)   ! read_snr data
  !synt_rad%prec = synt_rad_spec*(1./720)   ! snr720
  !synt_rad%nwavel = nwavel
  !synt_rad%nxtrack = nxtrack

  !--------------------------------
  ! Coadding Y-direction
  !--------------------------------
  gems_rad%wavl(:,:,:) = 0.0
  gems_rad%spec(:,:,:) = 0.0
  gems_rad%prec(:,:,:) = 0.0
  gems_rad%qflg(:,:,:) = pge_errstat_ok
  gems_rad%errstat(1:nxtrack, 1:ny) = pge_errstat_ok
  gems_rad%line_errstat( 1:ny) = pge_errstat_ok
  gems_rad%nwav(:,:) = nwavel 

  wavcals = .TRUE.

  gems_uv = synt_rad
 
  nx       = nxtrack
  DO iy = 1, ny

  mflgbits = 0
  iiy = (iy-1)*nybin + offline 
  Do i = 1,  nybin   
    iiy =  iiy +1

    !CALL convert_2bytes_to_16bits ( nbits, 1, gems_uv%mqflg(iiy), tmp_mflgbits(0:nbits-1))
    !mflgbits = mflgbits + tmp_mflgbits(0:nbits-1)

    !IF (correct_merr) THEN
       !tmpNinteg = 2.d0 / gems_uv%ExposureTime(iiy)
    !ENDIF      

    gems_rad%spec(1:nwavel, 1:nx, iy) = gems_rad%spec(1:nwavel, 1:nx, iy)+ gems_uv%spec(1:nwavel,1:nx, iiy)
    gems_rad%wavl(1:nwavel, 1:nx, iy) = gems_rad%wavl(1:nwavel, 1:nx, iy)+ gems_uv%wavl(1:nwavel,1:nx, iiy)

    !IF (correct_merr) THEN
       !gems_rad%prec(1:nwavel, 1:nx, iy) = gems_rad%prec(1:nwavel, 1:nx, iy)+gems_uv%prec(1:nwavel,1:nx, iiy) / SQRT( tmpNinteg ) 
    !ELSE
       !gems_rad%prec(1:nwavel, 1:nx, iy) = gems_rad%prec(1:nwavel, 1:nx, iy)+gems_uv%prec(1:nwavel,1:nx, iiy)
    !ENDIF
    gems_rad%prec(1:nwavel, 1:nx, iy) = gems_rad%prec(1:nwavel, 1:nx, iy)+gems_uv%prec(1:nwavel,1:nx, iiy)

    !DO ix = 1, nxtrack
      !CALL coadd_2bytes_qflgs(nbits, nwavel, gems_rad%qflg(1:nwavel, ix,iy), gems_uv%qflg(1:nwavel, ix,iiy)) 
    !ENDDO
  ENDDO ! End of Y-Coadding
          
  !IF (mflgbits(0) >= 1 .OR. mflgbits(1) >= 1 .OR. mflgbits(3) >= 1 .OR. mflgbits(12) >= 1) THEN
      !WRITE(*, *) 'All radiances could not be used: line ', iiy, ' Swath ', is
      !gems_rad%line_errstat(iy) = pge_errstat_error      
  !ELSE IF (ANY(mflgbits >= 1)) THEN
     !!WRITE(*, *) 'Warning set on all radiances: line',  blockline, ' Swath ', is
  !IF ( gems_rad%line_errstat(iy) /= pge_errstat_error) &
       !gems_rad%line_errstat(iy) = pge_errstat_warning      
     !! Over SAA region
      !IF (mflgbits(10) >= 1) gems_saa(iy) = 1
  !ENDIF


  ! load radiance data
  gems_rad%spec(1:nwavel, 1:nx, iy) = gems_rad%spec(1:nwavel, 1:nx,iy) / nybin
  gems_rad%wavl(1:nwavel, 1:nx, iy) = gems_rad%wavl(1:nwavel, 1:nx,iy) / nybin     
  gems_rad%prec(1:nwavel, 1:nx, iy) = gems_rad%prec(1:nwavel, 1:nx,iy) / nybin

  !IF ( correct_merr ) then !jbak
  !gems_rad%prec(1:nwavel, 1:nx, iy) = gems_rad%prec(1:nwavel, 1:nx,iy) / nybin / SQRT(1.0D0 * nybin)
  !ELSE
  !gems_rad%prec(1:nwavel, 1:nx, iy) = gems_rad%prec(1:nwavel, 1:nx,iy) / nybin
  !ENDIF

  ! Sort data in wavelength increasing order   
  IF (gems_rad%wavl(1, 1, iy) > gems_rad%wavl(nwavel, 1, iy)) THEN
    idxs(1:nwavel) = (/ (i, i = nwavel, 1, -1) /)
    gems_rad%wavl(1:nwavel, 1:nx, iy) = gems_rad%wavl(idxs(1:nwavel), 1:nx, iy)
    gems_rad%spec(1:nwavel, 1:nx, iy) = gems_rad%spec(idxs(1:nwavel),1:nx, iy)
    gems_rad%prec(1:nwavel, 1:nx, iy) = gems_rad%prec(idxs(1:nwavel),1:nx, iy)
    gems_rad%qflg(1:nwavel, 1:nx, iy) = gems_rad%qflg(idxs(1:nwavel), 1:nx, iy)     
  ENDIF 

  ENDDO ! end iloop


 ! Process for reduced resoution (not included for gems)
  IF (nwavel > nwavel_max) THEN
     WRITE(*, *) "Need to increase nwavel_max!!!"
     pge_error_status = pge_errstat_error; RETURN
  ENDIF

 !--------------------------------
 ! SubSet & Coadding X-direction
 !--------------------------------

  ! Determine number of binning for different fitting windows

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
      nbin = nxbin
      iix = (ix - 1) * nbin 
     ! properly align cross track positions to be coadded (should be within one pixel)

      IF (nbin > 2) CALL prespec_align(nwavel, nbin, gems_rad%wavl(1:nwavel,&
           iix+1:iix+nbin, iy), gems_rad%spec(1:nwavel, iix+1:iix+nbin, iy), &
           gems_rad%prec(1:nwavel, iix+1:iix+nbin, iy), &       
           gems_rad%qflg(1:nwavel, iix+1:iix+nbin, iy))
      
      !DO ic = 1, nbin
         !CALL convert_2bytes_to_16bits ( nbits, nwavel,gems_rad%qflg(1:nwavel, &
              !iix + ic, iy), flgbits(ic, 1:nwavel, 0:nbits-1))
         !flgmsks(1:nwavel) = flgmsks(1:nwavel) &
              !+ flgbits(ic, 1:nwavel, 0)                &   ! Missing
              !+ flgbits(ic, 1:nwavel, 1)                &   ! Bad 
              !+ flgbits(ic, 1:nwavel, 2)                &   ! Processing error
              !+ flgbits(ic, 1:nwavel, 4)                &   ! RTS_Pixel_Warning Flag
              !+ flgbits(ic, 1:nwavel, 5)                &   ! Saturation Possibility Flag
              !+ flgbits(ic, 1:nwavel, 7)                     ! Dark Current Warning Flag
      !ENDDO
           
      ! Subset valid data
      nsub = 0 ; subspec=0.0;fidx=1
       
      gems_rad%npix(1, ix, iy) = nsub
      nbin = nxbin 
      iix = (ix - 1) * nbin

      !fidx = gems_irrad%winpix(1, ix, 1) 
      !lidx = gems_irrad%winpix(1, ix, 2) 
      lidx = fidx + gems_irrad%npix(1, ix) - 1

      DO ii = fidx, lidx
        i = gems_irrad%wind(ii, ix)
        IF (ALL( gems_rad%spec(i, iix+1:iix+nbin, iy) > lower_spec) .AND. &
             ALL( gems_rad%spec(i, iix+1:iix+nbin, iy) < upper_spec) .AND. flgmsks(i) == 0 ) THEN
           nsub = nsub + 1
           subspec(1:nbin, wvl_idx, nsub) =  gems_rad%wavl(i, iix+1:iix+nbin, iy)
           subspec(1:nbin, spc_idx, nsub) =  gems_rad%spec(i, iix+1:iix+nbin, iy)
           subspec(1:nbin, sig_idx, nsub) =  gems_rad%prec(i, iix+1:iix+nbin, iy)
           gems_rad%wind(nsub, ix, iy) = ii
           !print * , i,ii,subspec(1, 1, nsub),subspec(1, sig_idx, nsub)
        ENDIF
      ENDDO
       
      fidx = lidx + 1
      gems_rad%npix(1, ix, iy) = nsub - gems_rad%npix(1, ix, iy)

      ! If the # of wavelengths is <= 75% of the # of irradiances, stop processing this pixel
      IF (gems_rad%npix(1, ix, iy) <= gems_irrad%npix(1, ix) * 0.8 ) THEN   !geun  0.9 -> 0.8
        WRITE(*, '(A,5I5,F9.2)') 'Too fewer radiance points: ', ix, iy, 1, &
              gems_rad%npix(1, ix, iy), gems_irrad%npix(1, ix), gems_sza(ix, iy)
        !gems_rad%errstat(ix, iy) = pge_errstat_error; EXIT
        gems_rad%errstat(ix, iy) = pge_errstat_error
stop; CYCLE  ! geun modified from EXIT to CYCLE
      ENDIF
     
      IF ( gems_rad%errstat(ix, iy) == pge_errstat_error) CYCLE   ! This pixel will not be processed.     
      gems_rad%nwav(ix, iy) = nsub

      ! Perform coadding if UV-2 is selected with UV-1
      fidx = 1     
      nbin = nxbin 
      lidx = fidx + gems_rad%npix(1, ix, iy) - 1 
      IF (nbin > 1) THEN
         CALL radwavcal_coadd(wcal_bef_coadd, wavcals(1, ix), 1, ix, gems_rad%npix(1, ix, iy), nbin, &
              subspec(1:nbin, :, fidx:lidx), wshis(1, ix, 1:nbin), wsqus(1, ix, 1:nbin), error)
         wavcals(1, ix) = .FALSE.
         IF (error) THEN
            WRITE(*, '(A)') 'No radiance wavelength calibration before coadding!!!'
            pge_error_status = pge_errstat_warning
         ENDIF
      ENDIF
      fidx = lidx + 1

      ! Get data for surface albedo & cloud fraction at 370.2 nm +/- 20 pixels
      irefl = 0; fidx = gems_refl%solwinpix(ix, 1)
      nbin = nxbin * ncoadd  
      iix = (ix - 1) * nbin

      DO i  = fidx, nwavel
        IF ( ALL(gems_rad%spec(i, iix+1:iix+nbin, iy) > lower_spec) .AND. &
           ALL(gems_rad%spec(i, iix+1:iix+nbin, iy) < upper_spec) .AND. flgmsks(i) == 0) THEN
           irefl = irefl + 1
           gems_refl%radwavl(irefl, ix, iy) = SUM(gems_rad%wavl(i, iix+1:iix+nbin, iy)) / nbin
           gems_refl%radspec(irefl, ix, iy) = SUM(gems_rad%spec(i, iix+1:iix+nbin, iy)) / nbin
         ! write(*,'(3i5,3f10.3)')  irefl, ix, iix+1, gems_refl%radwavl(irefl, ix, iy) ,gems_rad%wavl(i, iix+1:iix+nbin, iy)
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
      gems_rad%prec(1:nsub, ix, iy) = subspec(1, sig_idx, 1:nsub) /gems_rad%norm(ix, iy)       

    ENDDO ! END of x-direction
  ENDDO ! End of y-direction

 END SUBROUTINE SYNT_O3P_read_l1b_rad
