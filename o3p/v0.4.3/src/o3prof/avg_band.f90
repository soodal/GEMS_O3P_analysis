! Average radiance spectra and also update nradpix

SUBROUTINE avg_band_radspec(currspec, npts, errstat)

  USE OMSAO_precision_module
  USE OMSAO_indices_module,   ONLY: wvl_idx, spc_idx, sig_idx 
  USE OMSAO_variables_module, ONLY: numwin, nradpix, n_band_avg, n_band_samp, dwavmax
  USE OMSAO_errstat_module
  
  IMPLICIT NONE
  
  ! ========================
  !Input/Output variables
  ! ========================
  INTEGER, INTENT(INOUT)                                   :: npts
  INTEGER, INTENT(OUT)                                     :: errstat
  REAL (KIND=dp), DIMENSION(sig_idx, npts), INTENT (INOUT) :: currspec
  
  ! ===============
  ! Local variables
  ! ===============
  REAL (KIND=dp), DIMENSION(sig_idx, npts)  :: currspec_temp
  REAL (KIND=dp)                            :: sumsol, sumsig, sumwav
  INTEGER  :: i, iwin, j, k, npix, npts_out, fidx, lidx, navg, nsamp, ntemp
 
  errstat = pge_errstat_ok
  IF (npts /= SUM(nradpix(1:numwin))) THEN
     WRITE(*, *) 'avg_band_radspec: Number of points inconsistent!!!'
     errstat = pge_errstat_ok; RETURN
  ENDIF

  currspec_temp = 0.0    !initalize spectra to be returned

  fidx = 1; npts_out = 0
  DO iwin = 1, numwin
     lidx = fidx + nradpix(iwin) - 1
     navg = n_band_avg(iwin); nsamp = n_band_samp(iwin)
     
     IF ( navg == 1 .AND. nsamp == 1) THEN   ! copy all data
        currspec_temp(:, npts_out+1:npts_out+nradpix(iwin)) = & 
             currspec(:, fidx:lidx)
        npts_out = npts_out + nradpix(iwin) 
     ELSE IF (navg == 1 .AND. nsamp > 1) THEN ! copy data
        ntemp = npts_out
        DO  i = fidx, lidx, nsamp
           npts_out = npts_out + 1
           currspec_temp(:, npts_out) = currspec(:, i)
        ENDDO
        nradpix(iwin) = npts_out - ntemp
     ELSE                                     ! Need more averaging
        ntemp = npts_out
        DO  i = fidx, lidx, nsamp 
           npix = 0; sumsol=0.0; sumsig=0.0; sumwav=0.0
    
           DO k = i - navg/2, i + navg/2        
              ! IF navg points are not continuous, then do not average
              IF (k >= fidx .AND. k <= lidx) THEN 
                 IF (ABS(currspec(wvl_idx, k) - currspec(wvl_idx,i)) <= ABS(k-i) * dwavmax ) THEN
                    sumsol = sumsol + currspec(spc_idx, k) 
                    sumsig = sumsig + currspec(sig_idx, k)
                    sumwav = sumwav + currspec(wvl_idx, k) 
                    npix = npix + 1
                 END IF
              ENDIF
           ENDDO
           
           IF (npix >=  navg / 2 ) THEN
              npts_out = npts_out + 1
              currspec_temp(wvl_idx, npts_out) = sumwav / npix 
              currspec_temp(spc_idx, npts_out) = sumsol / npix
              currspec_temp(sig_idx, npts_out) = sumsig / npix     &
                  / SQRT(1.0D0 * navg )  !coadding reduces error
           END IF
        ENDDO
        nradpix(iwin) = npts_out - ntemp
     ENDIF
     fidx = lidx + 1

  ENDDO

  npts = npts_out
  currspec(:, 1:npts) = currspec_temp(:, 1:npts)

  RETURN

END SUBROUTINE avg_band_radspec


! NOTE: Used to average reference spectrum, which is already 
! convolved and interpolated to the same radiance grid
SUBROUTINE avg_band_refspec(pos, spec, npts_in, npts_out, errstat)

  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY: numwin, nradpix, n_band_avg, &
       n_band_samp, n_radwvl_sav, nradpix_sav, dwavmax, refnhextra
  USE OMSAO_errstat_module

  IMPLICIT NONE
 
  ! ========================
  !Input/Output variables
  ! ========================
  INTEGER, INTENT(IN)                                :: npts_in
  INTEGER, INTENT(OUT)                               :: errstat
  REAL (KIND=dp), DIMENSION(npts_in), INTENT (IN)    :: pos
  INTEGER, INTENT(OUT)                               :: npts_out
  REAL (KIND=dp), DIMENSION(npts_in), INTENT (INOUT) :: spec

  ! ===============
  ! Local variables
  ! ===============
  REAL (KIND=dp), DIMENSION(npts_in) :: spec_temp
  REAL (KIND=dp)                     :: sumsol
  INTEGER, DIMENSION(numwin)         :: winpts
  INTEGER  :: i, iwin, j, k, npix, fidx, lidx, navg, nsamp, nextra

  errstat = pge_errstat_ok

  nextra = refnhextra * 2
  IF (npts_in == n_radwvl_sav + numwin * nextra) THEN
     winpts(1:numwin) = nradpix_sav(1:numwin) + nextra
  ELSE IF (npts_in == SUM(nradpix(1:numwin)) + numwin * nextra) THEN
     npts_out = npts_in;      RETURN
  ELSE
     WRITE(*, *) 'avg_band_refspec: Number of points inconsistent!!!'
     WRITE(*, *) npts_in, n_radwvl_sav, SUM(nradpix(1:numwin))
     errstat = pge_errstat_ok; RETURN
  ENDIF
 
  spec_temp = 0.0 !initalize spectra to be returned

  fidx = 1; npts_out = 0
  DO iwin = 1, numwin
     lidx = fidx + winpts(iwin) - 1
     navg = n_band_avg(iwin); nsamp = n_band_samp(iwin)
     
     IF ( navg == 1 .AND. nsamp == 1) THEN    ! copy all data
        spec_temp(npts_out+1:npts_out+winpts(iwin)) = spec(fidx:lidx)
        npts_out = npts_out + winpts(iwin)       
     ELSE IF (navg == 1 .AND. nsamp > 1) THEN ! copy data
        DO i = fidx, fidx + refnhextra - 1
           npts_out = npts_out + 1; spec_temp(npts_out) = spec(i)
        ENDDO

        DO  i = fidx+refnhextra, lidx-refnhextra, nsamp
           npts_out = npts_out + 1; spec_temp(npts_out) = spec(i)
        ENDDO

        DO i = lidx - refnhextra + 1, lidx
           npts_out = npts_out + 1; spec_temp(npts_out) = spec(i)
        ENDDO
     ELSE                                     ! Need more averaging
        
        DO i = fidx, fidx + refnhextra - 1
           npts_out = npts_out + 1; spec_temp(npts_out) = spec(i)
        ENDDO

        DO  i = fidx+refnhextra, lidx-refnhextra, nsamp 
           npix = 0; sumsol=0.0
    
           DO k = i - navg/2, i + navg/2        
              ! IF navg points are not continuous, then do not average
              IF (k >= fidx + refnhextra .AND. k <= lidx-refnhextra) THEN
                 IF ( ABS(pos(k)-pos(i)) <= ABS(k-i) * dwavmax ) THEN
                    sumsol = sumsol + spec(k);  npix = npix + 1
                 ENDIF
              END IF
           END DO
          
           IF (npix >= navg / 2 ) THEN
              npts_out = npts_out + 1; spec_temp(npts_out) = sumsol / npix
           END IF
        ENDDO

        DO i = lidx - refnhextra + 1, lidx
           npts_out = npts_out + 1; spec_temp(npts_out) = spec(i)
        ENDDO

     ENDIF
     fidx = lidx + 1
  ENDDO

  spec(1:npts_out) = spec_temp(1:npts_out)

  RETURN

END SUBROUTINE avg_band_refspec

! NOTE: Used to average reference spectrum, which is already 
! convolved and interpolated to the same radiance grid
SUBROUTINE avg_band_effrefspec(pos, spec, npts_in, npts_out, errstat)

  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY: numwin, nradpix, n_band_avg, &
       n_band_samp, n_radwvl_sav, nradpix_sav, i0sav, dwavmax, refnhextra
  USE OMSAO_errstat_module

  IMPLICIT NONE
 
  ! ========================
  !Input/Output variables
  ! ========================
  INTEGER, INTENT(IN)                                :: npts_in
  INTEGER, INTENT(OUT)                               :: errstat
  REAL (KIND=dp), DIMENSION(npts_in), INTENT (IN)    :: pos
  INTEGER, INTENT(OUT)                               :: npts_out
  REAL (KIND=dp), DIMENSION(npts_in), INTENT (INOUT) :: spec

  ! ===============
  ! Local variables
  ! ===============
  REAL (KIND=dp), DIMENSION(npts_in) :: spec_temp
  REAL (KIND=dp)                     :: sumsol, sumf
  INTEGER, DIMENSION(numwin)         :: winpts
  INTEGER  :: i, iwin, j, k, npix, fidx, lidx, navg, nsamp, nextra

  errstat = pge_errstat_ok

  nextra = refnhextra * 2
  IF (npts_in == n_radwvl_sav + numwin * nextra) THEN
     winpts(1:numwin) = nradpix_sav(1:numwin) + nextra
  ELSE IF (npts_in == SUM(nradpix(1:numwin)) + numwin * nextra) THEN
     npts_out = npts_in;      RETURN
  ELSE
     WRITE(*, *) 'avg_band_effrefspec: Number of points inconsistent!!!'
     WRITE(*, *) npts_in, n_radwvl_sav, SUM(nradpix(1:numwin))
     errstat = pge_errstat_ok; RETURN
  ENDIF
 
  spec_temp = 0.0 !initalize spectra to be returned

  fidx = 1; npts_out = 0
  DO iwin = 1, numwin
     lidx = fidx + winpts(iwin) - 1
     navg = n_band_avg(iwin); nsamp = n_band_samp(iwin)
     
     IF ( navg == 1 .AND. nsamp == 1) THEN    ! copy all data
        spec_temp(npts_out+1:npts_out+winpts(iwin)) = spec(fidx:lidx)
        npts_out = npts_out + winpts(iwin)       
     ELSE IF (navg == 1 .AND. nsamp > 1) THEN ! copy data
        DO i = fidx, fidx + refnhextra - 1
           npts_out = npts_out + 1; spec_temp(npts_out) = spec(i)
        ENDDO

        DO  i = fidx+refnhextra, lidx-refnhextra, nsamp
           npts_out = npts_out + 1; spec_temp(npts_out) = spec(i)
        ENDDO

        DO i = lidx - refnhextra + 1, lidx
           npts_out = npts_out + 1; spec_temp(npts_out) = spec(i)
        ENDDO
     ELSE                                     ! Need more averaging
        
        DO i = fidx, fidx + refnhextra - 1
           npts_out = npts_out + 1; spec_temp(npts_out) = spec(i)
        ENDDO

        DO  i = fidx+refnhextra, lidx-refnhextra, nsamp 
           npix = 0; sumsol=0.0; sumf = 0.0
    
           DO k = i - navg/2, i + navg/2        
              ! IF navg points are not continuous, then do not average
              IF (k >= fidx + refnhextra .AND. k <= lidx-refnhextra) THEN
                 IF ( ABS(pos(k)-pos(i)) <= ABS(k-i) * dwavmax ) THEN
                    sumsol = sumsol + spec(k) * i0sav(k);  npix = npix + 1
                    sumf = sumf + i0sav(k)
                 ENDIF
              END IF
           END DO
          
           IF (npix >= navg / 2 ) THEN
              npts_out = npts_out + 1; spec_temp(npts_out) = sumsol / sumf
           END IF
        ENDDO

        DO i = lidx - refnhextra + 1, lidx
           npts_out = npts_out + 1; spec_temp(npts_out) = spec(i)
        ENDDO

     ENDIF
     fidx = lidx + 1
  ENDDO

  spec(1:npts_out) = spec_temp(1:npts_out)

  RETURN

END SUBROUTINE avg_band_effrefspec



! NOTE: Used to average reference spectrum, which is already 
! convolved and interpolated to the same radiance grid
SUBROUTINE avg_band_effozcrs(pos, spec, npts_in, npts_out, errstat)

  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY: numwin, nradpix, n_band_avg, &
       n_band_samp, n_radwvl_sav, nradpix_sav, i0sav, refidx_sav, dwavmax
  USE OMSAO_errstat_module

  IMPLICIT NONE
 
  ! ========================
  !Input/Output variables
  ! ========================
  INTEGER, INTENT(IN)                                :: npts_in
  INTEGER, INTENT(OUT)                               :: errstat
  REAL (KIND=dp), DIMENSION(npts_in), INTENT (IN)    :: pos
  INTEGER, INTENT(OUT)                               :: npts_out
  REAL (KIND=dp), DIMENSION(npts_in), INTENT (INOUT) :: spec

  ! ===============
  ! Local variables
  ! ===============
  REAL (KIND=dp), DIMENSION(npts_in) :: spec_temp
  REAL (KIND=dp)                     :: sumsol, sumf
  INTEGER, DIMENSION(numwin)         :: winpts
  INTEGER  :: i, iwin, j, k, npix, fidx, lidx, navg, nsamp

  errstat = pge_errstat_ok

  IF (npts_in == n_radwvl_sav ) THEN
     winpts = nradpix_sav(1:numwin)
  ELSE IF (npts_in == SUM(nradpix(1:numwin))) THEN
     npts_out = npts_in;      RETURN
  ELSE
     WRITE(*, *) 'avg_band_effozcrs: Number of points inconsistent!!!'
     WRITE(*, *) npts_in, n_radwvl_sav, SUM(nradpix(1:numwin))
     errstat = pge_errstat_ok; RETURN
  ENDIF
 
  spec_temp = 0.0 !initalize spectra to be returned

  fidx = 1; npts_out = 0
  DO iwin = 1, numwin
     lidx = fidx + winpts(iwin) - 1
     navg = n_band_avg(iwin); nsamp = n_band_samp(iwin)
     
     IF ( navg == 1 .AND. nsamp == 1) THEN    ! copy all data
        spec_temp(npts_out+1:npts_out+winpts(iwin)) = spec(fidx:lidx)
        npts_out = npts_out + winpts(iwin)       
     ELSE IF (navg == 1 .AND. nsamp > 1) THEN ! copy data
        DO  i = fidx, lidx, nsamp
           npts_out = npts_out + 1; spec_temp(npts_out) = spec(i)
        ENDDO
     ELSE                                     ! Need more averaging
        DO  i = fidx, lidx, nsamp 
           npix = 0; sumsol=0.0; sumf = 0.0
    
           DO k = i - navg/2, i + navg/2        
              ! IF navg points are not continuous, then do not average
              IF (k >= fidx .AND. k <= lidx) THEN
                 IF ( ABS(pos(k)-pos(i)) <= ABS(k-i) * dwavmax ) THEN
                    sumsol = sumsol + spec(k) * i0sav(refidx_sav(k))
                    sumf = sumf + i0sav(refidx_sav(k));  npix = npix + 1
                 ENDIF 
              END IF
           END DO
           
           IF (npix >= navg / 2 ) THEN
              npts_out = npts_out + 1; spec_temp(npts_out) = sumsol / sumf
           END IF
        ENDDO
     ENDIF
     fidx = lidx + 1
  ENDDO

  spec(1:npts_out) = spec_temp(1:npts_out)

  RETURN

END SUBROUTINE avg_band_effozcrs


! NOTE: Used to average reference spectrum, which is already 
! convolved and interpolated to the same radiance grid
SUBROUTINE avg_band_ozcrs(pos, spec, npts_in, npts_out, errstat)

  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY: numwin, nradpix, n_band_avg, &
       n_band_samp, n_radwvl_sav, nradpix_sav, i0sav, refidx_sav, dwavmax
  USE OMSAO_errstat_module

  IMPLICIT NONE
 
  ! ========================
  !Input/Output variables
  ! ========================
  INTEGER, INTENT(IN)                                :: npts_in
  INTEGER, INTENT(OUT)                               :: errstat
  REAL (KIND=dp), DIMENSION(npts_in), INTENT (IN)    :: pos
  INTEGER, INTENT(OUT)                               :: npts_out
  REAL (KIND=dp), DIMENSION(npts_in), INTENT (INOUT) :: spec

  ! ===============
  ! Local variables
  ! ===============
  REAL (KIND=dp), DIMENSION(npts_in) :: spec_temp
  REAL (KIND=dp)                     :: sumsol, sumf
  INTEGER, DIMENSION(numwin)         :: winpts
  INTEGER  :: i, iwin, j, k, npix, fidx, lidx, navg, nsamp

  errstat = pge_errstat_ok

  IF (npts_in == n_radwvl_sav ) THEN
     winpts = nradpix_sav(1:numwin)
  ELSE IF (npts_in == SUM(nradpix(1:numwin))) THEN
     npts_out = npts_in;      RETURN
  ELSE
     WRITE(*, *) 'avg_band_ozcrs: Number of points inconsistent!!!'
     WRITE(*, *) npts_in, n_radwvl_sav, SUM(nradpix(1:numwin))
     errstat = pge_errstat_ok; RETURN
  ENDIF
 
  spec_temp = 0.0 !initalize spectra to be returned

  fidx = 1; npts_out = 0
  DO iwin = 1, numwin
     lidx = fidx + winpts(iwin) - 1
     navg = n_band_avg(iwin); nsamp = n_band_samp(iwin)
     
     IF ( navg == 1 .AND. nsamp == 1) THEN    ! copy all data
        spec_temp(npts_out+1:npts_out+winpts(iwin)) = spec(fidx:lidx)
        npts_out = npts_out + winpts(iwin)       
     ELSE IF (navg == 1 .AND. nsamp > 1) THEN ! copy data
        DO  i = fidx, lidx, nsamp
           npts_out = npts_out + 1; spec_temp(npts_out) = spec(i)
        ENDDO
     ELSE                                     ! Need more averaging
        DO  i = fidx, lidx, nsamp 
           npix = 0; sumsol=0.0; sumf = 0.0
    
           DO k = i - navg/2, i + navg/2        
              ! IF navg points are not continuous, then do not average
              IF (k >= fidx .AND. k <= lidx) THEN
                 IF ( ABS(pos(k)-pos(i)) <= ABS(k-i) * dwavmax ) THEN
                    sumsol = sumsol + spec(k) !* i0sav(refidx_sav(k))
                    sumf = sumf + i0sav(refidx_sav(k));  npix = npix + 1
                 ENDIF 
              END IF
           END DO
           
           IF (npix >= navg / 2 ) THEN
              npts_out = npts_out + 1; spec_temp(npts_out) = sumsol / npix !sumf
           END IF
        ENDDO
     ENDIF
     fidx = lidx + 1
  ENDDO

  spec(1:npts_out) = spec_temp(1:npts_out)

  RETURN

END SUBROUTINE avg_band_ozcrs



SUBROUTINE avg_band_subrefspec(pos, spec, npts_in, iwin, nedge, npts_out)

  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY: numwin, nradpix, n_band_avg, &
       n_band_samp, n_radwvl_sav, nradpix_sav, dwavmax
  USE OMSAO_errstat_module

  IMPLICIT NONE
 
  ! ========================
  !Input/Output variables
  ! ========================
  INTEGER, INTENT(IN)                                :: npts_in, iwin, nedge
  REAL (KIND=dp), DIMENSION(npts_in), INTENT (IN)    :: pos
  INTEGER, INTENT(OUT)                               :: npts_out
  REAL (KIND=dp), DIMENSION(npts_in), INTENT (INOUT) :: spec

  ! ===============
  ! Local variables
  ! ===============
  REAL (KIND=dp), DIMENSION(npts_in) :: spec_temp
  REAL (KIND=dp)                     :: sumsol
  INTEGER, DIMENSION(numwin)         :: winpts
  INTEGER  :: i, j, k, npix, fidx, lidx, navg, nsamp

  spec_temp = 0.0 !initalize spectra to be returned

  fidx = 1; lidx = npts_in; npts_out = 0
  navg = n_band_avg(iwin); nsamp = n_band_samp(iwin)
     
  IF ( navg == 1 .AND. nsamp == 1) THEN    ! copy all data
     spec_temp(npts_out+1:npts_out+winpts(iwin)) = spec(fidx:lidx)
     npts_out = npts_out + winpts(iwin)       
  ELSE IF (navg == 1 .AND. nsamp > 1) THEN ! copy data
     DO i = fidx, fidx + nedge-1
        npts_out = npts_out + 1; spec_temp(npts_out) = spec(i)
     ENDDO
     
     DO  i = fidx + nedge, lidx - nedge, nsamp
        npts_out = npts_out + 1; spec_temp(npts_out) = spec(i)
     ENDDO
     
     DO i = lidx - nedge + 1, lidx
        npts_out = npts_out + 1; spec_temp(npts_out) = spec(i)
     ENDDO
  ELSE                                     ! Need more averaging
     
     DO i = fidx, fidx + nedge-1
        npts_out = npts_out + 1; spec_temp(npts_out) = spec(i)
     ENDDO
     
     DO  i = fidx + nedge, lidx - nedge, nsamp 
        npix = 0; sumsol=0.0
        
        DO k = i - navg/2, i + navg/2        
           ! IF navg points are not continuous, then do not average
           IF (k >= fidx + nedge .AND. k <= lidx - nedge) THEN
              IF ( ABS(pos(k)-pos(i)) <= ABS(k-i) * dwavmax ) THEN
                 sumsol = sumsol + spec(k);  npix = npix + 1
              ENDIF
           END IF
        END DO
        
        IF (npix >= navg / 2 ) THEN
           npts_out = npts_out + 1; spec_temp(npts_out) = sumsol / npix
        END IF
     ENDDO
     
     DO i = lidx - nedge + 1, lidx
        npts_out = npts_out + 1; spec_temp(npts_out) = spec(i)
     ENDDO     
  ENDIF

  spec(1:npts_out) = spec_temp(1:npts_out)

  RETURN

END SUBROUTINE avg_band_subrefspec

! NOTE: Used to average reference spectrum, which is already 
! convolved and interpolated to the same radiance grid
SUBROUTINE avg_band_spec(pos, spec, npts_in, npts_out, errstat)

  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY: numwin, nradpix_sav, n_band_avg, &
       n_band_samp, dwavmax
  USE OMSAO_errstat_module

  IMPLICIT NONE
 
  ! ========================
  !Input/Output variables
  ! ========================
  INTEGER, INTENT(IN)                                :: npts_in
  INTEGER, INTENT(OUT)                               :: errstat
  REAL (KIND=dp), DIMENSION(npts_in), INTENT (IN)    :: pos
  INTEGER, INTENT(OUT)                               :: npts_out
  REAL (KIND=dp), DIMENSION(npts_in), INTENT (INOUT) :: spec

  ! ===============
  ! Local variables
  ! ===============
  REAL (KIND=dp), DIMENSION(npts_in) :: spec_temp
  REAL (KIND=dp)                     :: sumsol, sumf
  INTEGER, DIMENSION(numwin)         :: winpts
  INTEGER  :: i, iwin, j, k, npix, fidx, lidx, navg, nsamp

  errstat = pge_errstat_ok

  IF (npts_in /= SUM(nradpix_sav(1:numwin))) THEN
     WRITE(*, *) 'avg_band_spec: Number of points inconsistent!!!'
     WRITE(*, *) npts_in, SUM(nradpix_sav(1:numwin))
     errstat = pge_errstat_ok; RETURN
  ENDIF
 
  spec_temp = 0.0 !initalize spectra to be returned

  fidx = 1; npts_out = 0
  DO iwin = 1, numwin
     lidx = fidx + nradpix_sav(iwin) - 1
     navg = n_band_avg(iwin); nsamp = n_band_samp(iwin)
    
     IF ( navg == 1 .AND. nsamp == 1) THEN    ! copy all data
        spec_temp(npts_out+1:npts_out+nradpix_sav(iwin)) = spec(fidx:lidx)
        npts_out = npts_out + nradpix_sav(iwin)       
     ELSE IF (navg == 1 .AND. nsamp > 1) THEN ! copy data
        DO  i = fidx, lidx, nsamp
           npts_out = npts_out + 1; spec_temp(npts_out) = spec(i)
        ENDDO
     ELSE                                     ! Need more averaging
        DO  i = fidx, lidx, nsamp 
           npix = 0; sumsol=0.0
    
           DO k = i - navg/2, i + navg/2        
              ! IF navg points are not continuous, then do not average
              IF (k >= fidx .AND. k <= lidx) THEN
                 IF ( ABS(pos(k)-pos(i)) <= ABS(k-i) * dwavmax ) THEN
                    sumsol = sumsol + spec(k); npix = npix + 1
                 ENDIF 
              END IF
           END DO
           
           IF (npix >= navg / 2 ) THEN
              npts_out = npts_out + 1; spec_temp(npts_out) = sumsol / npix
           END IF
        ENDDO
     ENDIF
     fidx = lidx + 1
  ENDDO

  spec(1:npts_out) = spec_temp(1:npts_out)

  RETURN

END SUBROUTINE avg_band_spec
