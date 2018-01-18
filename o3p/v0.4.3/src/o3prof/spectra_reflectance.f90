  ! ******************************************************************************
  ! Author:  xiong liu
  ! Date  :  July 17, 2003
  ! Purpose: 1. Calibrate solar and radiance spectrum to compute reflectance
  !          for comparing with the LIDORT-calculated radiance
  !          2. Also correct absorption by other species other than ozone because
  !          these species are not modelled in the LIDORT  
  ! ******************************************************************************

SUBROUTINE spectra_reflectance (ns, nf, fitvar, do_shiwf, fitspec, errstat)

  USE OMSAO_precision_module
  USE OMSAO_indices_module, ONLY    : &
       max_rs_idx, max_calfit_idx, solar_idx, ring_idx, ad1_idx, &
       lbe_idx, ad2_idx, mxs_idx, wvl_idx, spc_idx, comm_idx, &
       sin_idx, hwe_idx, asy_idx, shi_idx, squ_idx, us1_idx, us2_idx, &
       com1_idx, no2_t1_idx, no2_t2_idx, &
       so2_idx, so2v_idx, bro_idx, bro2_idx, shift_offset, fsl_idx, rsl_idx
  USE OMSAO_variables_module,  ONLY : fitvar_rad, mask_fitvar_rad, rad_wav_avg, &
       database, currspec, fitwavs, n_refwvl, refwvl, lo_radbnd, refnhextra, &
       up_radbnd, database_shiwf, slwf, numwin, nradpix, refspec_norm, refidx
  USE ozprof_data_module,      ONLY : div_rad, div_sun, use_lograd, do_subfit, &
       osind, osfind, slind, slfind, shind, shfind, rnind, rnfind, dcind,      &
       dcfind, isind, isfind, irind, irfind, oswins, slwins, shwins, &
       rnwins, dcwins, iswins, irwins, nos, nsl, nsh, nrn, ndc, nis, nir, &
       fit_atanring
  USE OMSAO_errstat_module
  
  IMPLICIT NONE


  ! ========================
  ! Input/Output variables
  ! ========================
  INTEGER,  INTENT (IN)  :: ns, nf
  INTEGER,  INTENT (OUT) :: errstat
  LOGICAL,  INTENT (IN)  :: do_shiwf
  REAL (KIND=dp), DIMENSION (nf), INTENT (IN)  :: fitvar
  REAL (KIND=dp), DIMENSION (ns), INTENT (OUT) :: fitspec
  
  ! =================================
  ! External OMI and Toolkit routines
  ! =================================
  INTEGER :: OMI_SMF_setmsg

  ! ===============
  ! Local variables
  ! ===============
  INTEGER  :: i, j, k, idx, fidx, lidx, nextra
  REAL (KIND=dp), DIMENSION (ns)      :: del, sunspec_ss, locdbs, dfdsl, gshiwf, tempsum
  REAL (KIND=dp), DIMENSION (n_refwvl):: sunpos_ss, refpos_ss, delref
  REAL (KIND=dp), DIMENSION (ns)      :: corrected_irrad, corr, corrected_rad
  REAL (KIND=dp)                      :: wavg
  LOGICAL                             :: cal_shiwf

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=19), PARAMETER      :: modulename = 'spectra_reflectance'

  errstat = pge_errstat_ok

  ! Obtain the uncondensed array of variables
  fitvar_rad(mask_fitvar_rad(1:nf)) = fitvar(1:nf)

  sunpos_ss(1:n_refwvl) =  refwvl(1:n_refwvl)
  IF (nsh > 0) THEN
     nextra = 2 * refnhextra
     IF (do_subfit) THEN
        fidx = 1
        DO i = 1, numwin
           lidx =  fidx + nradpix(i) + nextra - 1
           delref(fidx:lidx) = refwvl(fidx:lidx) - (refwvl(fidx)+refwvl(lidx))*0.5
           IF (shfind(i, 1) > 0) sunpos_ss(fidx:lidx) = sunpos_ss(fidx:lidx) + fitvar_rad(shind(i, 1)) 
  
           DO j = 2, nsh
              IF (shfind(i, j) > 0) sunpos_ss(fidx:lidx) = sunpos_ss(fidx:lidx) +  &
                   fitvar_rad(shind(i, j)) * delref(fidx:lidx)**(j-1)
           ENDDO
           fidx = lidx + 1
        ENDDO
     ELSE
        IF (shwins(1, 1) == 1) THEN
           fidx = 1
        ELSE
           fidx = SUM(nradpix(1: shwins(1, 1)-1)) + 1 + (shwins(1, 1)-1) * nextra
        ENDIF       
        lidx = SUM(nradpix(1: shwins(1, 2))) +  shwins(1, 2) * nextra
        delref(fidx:lidx) = refwvl(fidx:lidx) - (refwvl(fidx)+refwvl(lidx))*0.5
        IF (shfind(1, 1) > 0) sunpos_ss(fidx:lidx) = sunpos_ss(fidx:lidx) + fitvar_rad(shind(1, 1)) 
        
        DO j= 2, nsh
           IF (shfind(1, j) > 0) sunpos_ss(fidx:lidx) = sunpos_ss(fidx:lidx) +  &
                fitvar_rad(shind(1, j)) * delref(fidx:lidx)**(j-1)
        ENDDO
     ENDIF
  ENDIF

  ! Re-sample the solar reference spectrum to the radiance grid
  CALL interpolation (n_refwvl, sunpos_ss(1:n_refwvl), database(solar_idx, 1:n_refwvl), &
       ns, fitwavs(1:ns), sunspec_ss(1:ns), errstat)  
  IF ( errstat >= pge_errstat_warning ) THEN
     !errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0)
     WRITE(*, '(2A,I3)') modulename, ': interpolation out of bounds or not ascending: ', errstat
     errstat = pge_errstat_error; RETURN
  END IF

  IF (do_shiwf .AND. nsl > 0) THEN
     CALL interpolation(n_refwvl, sunpos_ss(1:n_refwvl), &
          database_shiwf(solar_idx,1:n_refwvl), ns, fitwavs(1:ns), dfdsl(1:ns), errstat)

     IF ( errstat >= pge_errstat_warning ) THEN
        !errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0)
        WRITE(*, '(2A,I3)') modulename, ': interpolation out of bounds or not ascending: ', errstat
        errstat = pge_errstat_error; RETURN
     END IF
     
     IF (use_lograd) THEN
        slwf(1:ns) = -dfdsl(1:ns) / sunspec_ss(1:ns)
     ELSE
        slwf(1:ns) = -dfdsl(1:ns)*currspec(1:ns)/(sunspec_ss(1:ns)**2)
     END IF
  END IF

  ! correct for irradiance due to using slit width different from radiance spectrum
  IF ( nsl > 0 ) THEN
     IF ( do_subfit ) THEN
        fidx = 1
        DO i = 1, numwin
           lidx =  fidx + nradpix(i) - 1
           del(fidx:lidx) = fitwavs(fidx:lidx) - (fitwavs(fidx)+fitwavs(lidx)) * 0.5
           IF (slfind(i, 1) > 0) sunspec_ss(fidx:lidx) = sunspec_ss(fidx:lidx) + &
                   fitvar_rad(slind(i, 1)) 
           DO j = 2, nsl
              IF (slfind(i, j) > 0) sunspec_ss(fidx:lidx) = sunspec_ss(fidx:lidx) + &
                   fitvar_rad(slind(i, j)) * del(fidx:lidx)**(j-1)
           ENDDO
           fidx = lidx + 1
        ENDDO
     ELSE
        IF (slwins(1, 1) == 1) THEN
           fidx = 1
        ELSE
           fidx = SUM(nradpix(1: slwins(1, 1)-1)) + 1 
        ENDIF
        lidx = SUM(nradpix(1: slwins(1, 2))) 
        del(fidx:lidx) = fitwavs(fidx:lidx) - (fitwavs(fidx)+fitwavs(lidx)) * 0.5 
        IF (slfind(1, 1) > 0) sunspec_ss(fidx:lidx) = sunspec_ss(fidx:lidx) + &
             fitvar_rad(slind(1, 1)) 
        DO j = 2, nsl
           IF (slfind(1, j) > 0) sunspec_ss(fidx:lidx) = sunspec_ss(fidx:lidx) + &
                fitvar_rad(slind(1, j)) * del(fidx:lidx)**(j-1)
        ENDDO
     ENDIF
  ENDIF

  ! Degradation Correction
  IF ( ndc > 0 ) THEN
     corr = 0.0

     IF (do_subfit) THEN
        fidx = 1
        DO i = 1, numwin
           lidx =  fidx + nradpix(i) - 1
           del(fidx:lidx) = fitwavs(fidx:lidx) - (fitwavs(fidx) + fitwavs(lidx)) * 0.5

           IF (dcfind(i, 1) > 0) corr(fidx:lidx) = corr(fidx:lidx) + fitvar_rad(dcind(i, 1)) 

           DO j = 2, ndc
              IF (dcfind(i, j) > 0) corr(fidx:lidx) = corr(fidx:lidx) + &
                   fitvar_rad(dcind(i, j)) * del(fidx:lidx)**(j-1)
           ENDDO
           fidx = lidx + 1
        ENDDO
     ELSE
        IF (dcwins(1, 1) == 1) THEN
           fidx = 1
        ELSE
           fidx = SUM(nradpix(1: dcwins(1, 1)-1)) + 1 
        ENDIF
        lidx = SUM(nradpix(1: dcwins(1, 2))) 
        del(fidx:lidx) = fitwavs(fidx:lidx) - (fitwavs(fidx) + fitwavs(lidx)) * 2.0
        IF (dcfind(1, 1) > 0) corr(fidx:lidx) = corr(fidx:lidx) + fitvar_rad(dcind(1, 1)) 
        
        DO j = 2, ndc
           IF (dcfind(1, j) > 0) corr(fidx:lidx) = corr(fidx:lidx) + &
                fitvar_rad(dcind(1, j)) * del(fidx:lidx)**(j-1)
        ENDDO
     ENDIF
     sunspec_ss = sunspec_ss * EXP (corr)
  ENDIF
  !WRITE(91, '(F10.4, D14.6)') ((fitwavs(i), sunspec_ss(i)), i = 1, ns)

  ! Internal Scattering in Irradiance
  IF ( nis > 0 ) THEN
     IF (do_subfit) THEN
        fidx = 1
        DO i = 1, numwin
           lidx =  fidx + nradpix(i) - 1
           wavg = (fitwavs(fidx) + fitwavs(lidx)) * 0.5
           del(fidx:lidx) = fitwavs(fidx:lidx) - wavg
           tempsum(fidx:lidx) = SUM(sunspec_ss(fidx:lidx)) / wavg * fitwavs(fidx:lidx)

           IF (isfind(i, 1) > 0) sunspec_ss(fidx:lidx) = sunspec_ss(fidx:lidx) &
                + fitvar_rad(isind(i, 1)) * tempsum(fidx:lidx)
           DO j = 2, nis
              IF (isfind(i, j) > 0) sunspec_ss(fidx:lidx) = sunspec_ss(fidx:lidx) &
                   + fitvar_rad(isind(i, j)) * tempsum(fidx:lidx) * del(fidx:lidx) ** (j-1) 
           ENDDO
           fidx = lidx + 1
        ENDDO
     ELSE
        IF (iswins(1, 1) == 1) THEN
           fidx = 1
        ELSE
           fidx = SUM(nradpix(1: iswins(1, 1)-1)) + 1 
        ENDIF
        lidx = SUM(nradpix(1: iswins(1, 2)))
 
        wavg = (fitwavs(fidx) + fitwavs(lidx)) * 0.5
        del(fidx:lidx) = fitwavs(fidx:lidx) - wavg
        tempsum(fidx:lidx) = SUM(sunspec_ss(fidx:lidx)) / wavg * fitwavs(fidx:lidx)

        IF (isfind(1, 1) > 0) sunspec_ss(fidx:lidx) = sunspec_ss(fidx:lidx) &
                + fitvar_rad(isind(1, 1)) * tempsum(fidx:lidx)
        DO j = 2, nis
           IF (isfind(1, j) > 0) sunspec_ss(fidx:lidx) = sunspec_ss(fidx:lidx) &
                + fitvar_rad(isind(1, j)) * tempsum(fidx:lidx) * del(fidx:lidx) ** (j-1)
        ENDDO
     ENDIF
  ENDIF
  corrected_irrad = sunspec_ss 

!  ! Internal Scattering in Irradiance (use stray-light spectra)
!  IF ( nis > 0 ) THEN
!     IF (do_subfit) THEN
!        fidx = 1
!        DO i = 1, numwin
!           lidx =  fidx + nradpix(i) - 1
!           wavg = (fitwavs(fidx) + fitwavs(lidx)) * 0.5
!           del(fidx:lidx) = fitwavs(fidx:lidx) - wavg
!
!           IF (isfind(i, 1) > 0) sunspec_ss(fidx:lidx) = sunspec_ss(fidx:lidx) &
!                + fitvar_rad(isind(i, 1)) * database(fsl_idx, refidx(fidx:lidx))
!           DO j = 2, nis
!              IF (isfind(i, j) > 0) sunspec_ss(fidx:lidx) = sunspec_ss(fidx:lidx) + &
!                   fitvar_rad(isind(i, j)) * database(fsl_idx, refidx(fidx:lidx)) * del(fidx:lidx) ** (j-1) 
!           ENDDO
!           fidx = lidx + 1
!        ENDDO
!     ELSE
!        IF (iswins(1, 1) == 1) THEN
!           fidx = 1
!        ELSE
!           fidx = SUM(nradpix(1: iswins(1, 1)-1)) + 1 
!        ENDIF
!        lidx = SUM(nradpix(1: iswins(1, 2)))
! 
!        wavg = (fitwavs(fidx) + fitwavs(lidx)) * 0.5
!        del(fidx:lidx) = fitwavs(fidx:lidx) - wavg
!
!        IF (isfind(1, 1) > 0) sunspec_ss(fidx:lidx) = sunspec_ss(fidx:lidx) &
!                + fitvar_rad(isind(1, 1)) * database(fsl_idx, refidx(fidx:lidx))
!        DO j = 2, nis
!           IF (isfind(1, j) > 0) sunspec_ss(fidx:lidx) = sunspec_ss(fidx:lidx) + &
!                fitvar_rad(isind(1, j)) * database(fsl_idx, refidx(fidx:lidx)) * del(fidx:lidx) ** (j-1)
!        ENDDO
!     ENDIF
!  ENDIF
!  corrected_irrad = sunspec_ss 

  ! Initial add-on contributions to irradiance: undersampling
  DO j = no2_t1_idx, max_rs_idx
     i = max_calfit_idx + (j-1) * mxs_idx + ad1_idx
     IF (fitvar_rad(i) == 0.0) CYCLE                !  no such parameter
        
     k = shift_offset + j  ! corresponding shift parameter
     IF (fitvar_rad(k) == 0.0) THEN
        locdbs(1:ns) = database(j, refidx(1:ns))
     ELSE 
        refpos_ss(1:n_refwvl) = refwvl(1:n_refwvl) + fitvar_rad(k)          
        cal_shiwf = .FALSE.   
        IF (lo_radbnd(k) < up_radbnd(k) .AND. do_shiwf) cal_shiwf = .TRUE.
        CALL bspline2(refpos_ss(1:n_refwvl)-fitvar_rad(k), database(j,1:n_refwvl), n_refwvl, &
             cal_shiwf,fitwavs(1:ns), locdbs(1:ns),  gshiwf, ns, errstat)
        database_shiwf(j, refidx(1:ns)) = gshiwf
        IF (errstat < 0) THEN
           WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
           errstat = pge_errstat_error; RETURN
        ENDIF
     END IF
        
     corrected_irrad(1:ns) = corrected_irrad(1:ns) + fitvar_rad(i) * locdbs(1:ns)
  END DO
  !DO i = 1, ns
  !   WRITE(*, '(F10.4,D14.7)') fitwavs(i), corrected_irrad(i) * div_sun
  !ENDDO
  !STOP

  ! Internal Scattering in Radiance
  corrected_rad = currspec(1:ns)
  IF ( nir > 0 ) THEN
     IF (do_subfit) THEN
        fidx = 1
        DO i = 1, numwin
           lidx =  fidx + nradpix(i) - 1

           wavg = (fitwavs(fidx) + fitwavs(lidx)) * 0.5
           del(fidx:lidx) = fitwavs(fidx:lidx) - wavg
           tempsum(fidx:lidx) = SUM(corrected_rad(fidx:lidx)) / wavg * fitwavs(fidx:lidx)
           !tempsum(fidx:lidx) = 1.0

           IF (irfind(i, 1) > 0) corrected_rad(fidx:lidx) = corrected_rad(fidx:lidx) &
                   + fitvar_rad(irind(i, 1)) * tempsum(fidx:lidx)
           DO j = 2, nir
              IF (irfind(i, j) > 0) corrected_rad(fidx:lidx) = corrected_rad(fidx:lidx) &
                   + fitvar_rad(irind(i, j)) * tempsum(fidx:lidx) * del(fidx:lidx) ** (j-1)  
           ENDDO
           fidx = lidx + 1
        ENDDO
     ELSE
        IF (irwins(1, 1) == 1) THEN
           fidx = 1
        ELSE
           fidx = SUM(nradpix(1: irwins(1, 1)-1)) + 1 
        ENDIF
        lidx = SUM(nradpix(1: irwins(1, 2))) 
        wavg = (fitwavs(fidx) + fitwavs(lidx)) * 0.5
        del(fidx:lidx) = fitwavs(fidx:lidx) - wavg
        tempsum(fidx:lidx) = SUM(corrected_rad(fidx:lidx)) / wavg * fitwavs(fidx:lidx)
        
        IF (irfind(1, 1) > 0) corrected_rad(fidx:lidx) = corrected_rad(fidx:lidx) &
                   + fitvar_rad(irind(1, 1)) * tempsum(fidx:lidx) 
        
        DO j = 2, nir
           IF (irfind(1, j) > 0) corrected_rad(fidx:lidx) = corrected_rad(fidx:lidx) &
                + fitvar_rad(irind(1, j)) * tempsum(fidx:lidx) * del(fidx:lidx) ** (j-1)
        ENDDO
     ENDIF
  ENDIF
  fitspec = corrected_rad / corrected_irrad * div_rad / div_sun 

! ! Internal Scattering in Radiance  (fit straylight spectrum)
!  corrected_rad = currspec(1:ns)
!  IF ( nir > 0 ) THEN
!     IF (do_subfit) THEN
!        fidx = 1
!        DO i = 1, numwin
!           lidx =  fidx + nradpix(i) - 1
!     
!           wavg = (fitwavs(fidx) + fitwavs(lidx)) * 0.5
!           del(fidx:lidx) = fitwavs(fidx:lidx) - wavg
!     
!           IF (irfind(i, 1) > 0) corrected_rad(fidx:lidx) = corrected_rad(fidx:lidx) * &
!                   (1.0 - fitvar_rad(irind(i, 1)) * database(rsl_idx, refidx(fidx:lidx)))
!           DO j = 2, nir
!              IF (irfind(i, j) > 0) corrected_rad(fidx:lidx) = corrected_rad(fidx:lidx) *  &
!                (1.0 - fitvar_rad(irind(i, j)) * database(rsl_idx, refidx(fidx:lidx)) * del(fidx:lidx) ** (j-1))  
!           ENDDO
!           fidx = lidx + 1
!        ENDDO
!     ELSE
!        IF (irwins(1, 1) == 1) THEN
!           fidx = 1
!        ELSE
!           fidx = SUM(nradpix(1: irwins(1, 1)-1)) + 1 
!        ENDIF
!        lidx = SUM(nradpix(1: irwins(1, 2))) 
!        wavg = (fitwavs(fidx) + fitwavs(lidx)) * 0.5
!        del(fidx:lidx) = fitwavs(fidx:lidx) - wavg
!        
!        IF (irfind(1, 1) > 0) corrected_rad(fidx:lidx) = corrected_rad(fidx:lidx) * &
!                  (1.0 - fitvar_rad(irind(1, 1)) * database(rsl_idx, refidx(fidx:lidx)))
!        
!        DO j = 2, nir
!           IF (irfind(1, j) > 0) corrected_rad(fidx:lidx) = corrected_rad(fidx:lidx) * &
!                (1.0 - fitvar_rad(irind(1, j)) * database(rsl_idx, refidx(fidx:lidx)) * del(fidx:lidx) ** (j-1))
!        ENDDO
!     ENDIF
!  ENDIF
!  fitspec = corrected_rad / corrected_irrad * div_rad / div_sun  
  
  ! Beer's law contributions for species other than ozone profile, because 
  ! measured radiance was already contaminated by absorption other than ozone, 
  ! we need to multiply exp(+...) in order to match the modeled reflectance.
  ! Now it is being treated in LIDORT directly
  !DO j = no2_t1_idx, max_rs_idx
  !   IF (j == comm_idx .OR. j == com1_idx) THEN
  !      i = max_calfit_idx + (j-1) * mxs_idx + lbe_idx
  !      
  !      IF (fitvar_rad(i) == 0.0) CYCLE !no such parameter
  !      
  !      k = shift_offset + j
  !      IF (fitvar_rad(k) == 0.0) THEN
  !         locdbs(1:ns) = database(j, refidx(1:ns))
  !      ELSE 
  !         refpos_ss(1:n_refwvl) = refwvl(1:n_refwvl) + fitvar_rad(k)       
  !         cal_shiwf = .FALSE.           
  !         IF (lo_radbnd(k) < up_radbnd(k) .AND. do_shiwf) cal_shiwf = .TRUE.
  !         CALL bspline2(refpos_ss(1:n_refwvl)-fitvar_rad(k), database(j,1:n_refwvl), n_refwvl, &
  !              cal_shiwf, fitwavs(1:ns), locdbs(1:ns), gshiwf, ns, errstat)
  !         database_shiwf(j, refidx(1:ns)) = gshiwf
  !         IF (errstat < 0) THEN
  !            WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
  !            errstat = pge_errstat_error; RETURN
  !         ENDIF
  !      END IF
  !      
  !      fitspec(1:ns) = fitspec(1:ns) * EXP( fitvar_rad(i) * locdbs(1:ns))
  !   ENDIF
  !END DO

  ! Ring effect (using Beer-law contribution)
  IF (.NOT. fit_atanring) THEN
     IF ( nrn > 0 ) THEN
        corr = 0.0
        
        IF (do_subfit) THEN
           fidx = 1
           DO i = 1, numwin
              lidx =  fidx + nradpix(i) - 1
              
              del(fidx:lidx) = fitwavs(fidx:lidx) - ( (fitwavs(fidx) + fitwavs(lidx)) * 0.5)
              
              IF (rnfind(i, 1) > 0) corr(fidx:lidx) = corr(fidx:lidx) + fitvar_rad(rnind(i, 1))
              
              DO j = 2, nrn
                 IF (rnfind(i, j) > 0) corr(fidx:lidx) = corr(fidx:lidx) + &
                      fitvar_rad(rnind(i, j)) * del(fidx:lidx) ** (j-1)
              ENDDO
              
              fidx = lidx + 1
           ENDDO
        ELSE
           IF (rnwins(1, 1) == 1) THEN
              fidx = 1
           ELSE
              fidx = SUM(nradpix(1: rnwins(1, 1)-1)) + 1 
           ENDIF
           lidx = SUM(nradpix(1: rnwins(1, 2))) 
           
           del(fidx:lidx)   = fitwavs(fidx:lidx) - (fitwavs(fidx) + fitwavs(lidx)) * 0.5
           
           IF (rnfind(1, 1) > 0) corr(fidx:lidx) = corr(fidx:lidx) + fitvar_rad(rnind(1, 1))
           
           DO j = 2, nrn
              IF (rnfind(1, j) > 0) corr(fidx:lidx) = corr(fidx:lidx) + &
                   fitvar_rad(rnind(1, j)) * del(fidx:lidx) ** (j-1)
           ENDDO
        ENDIF
        
        fitspec = fitspec * EXP(database(ring_idx, refidx(1:ns)) * corr)
     ENDIF
  ELSE
     ! Fitting Ring effect using y = -[a0 * [atan((x-a1)/a2) + 1.54223] + 1]
     corr = atan( (fitwavs - fitvar_rad(rnind(1, 2))) / fitvar_rad(rnind(1, 3)) )
     corr = -(fitvar_rad(rnind(1, 1)) * (corr - corr(1)) + 1.0)
     fitspec = fitspec * EXP(database(ring_idx, refidx(1:ns)) * corr)
  ENDIF

  ! Second add-on contributions: add to the Sun-normalized radiance 
  ! (e.g. commomd mode, ring filling in)
  DO j = no2_t1_idx, max_rs_idx
     i = max_calfit_idx + (j-1) * mxs_idx + ad2_idx
     IF (fitvar_rad(i) == 0.0) CYCLE ! no such parameter
     
     k = shift_offset + j
     IF (fitvar_rad(k) == 0.0) THEN
        locdbs(1:ns) = database(j, refidx(1:ns))
     ELSE 
        refpos_ss(1:n_refwvl) = refwvl(1:n_refwvl) + fitvar_rad(k)
       
        cal_shiwf = .FALSE.
        IF (lo_radbnd(k) < up_radbnd(k) .AND. do_shiwf) cal_shiwf = .TRUE.
        CALL bspline2(refpos_ss(1:n_refwvl)-fitvar_rad(k), database(j,1:n_refwvl), n_refwvl, &
             cal_shiwf, fitwavs(1:ns), locdbs(1:ns), gshiwf, ns, errstat)
        database_shiwf(j, refidx(1:ns)) = gshiwf
        IF (errstat < 0) THEN
           WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
           errstat = pge_errstat_error; RETURN
        ENDIF
     END IF
     
     IF (j == comm_idx .OR. j == com1_idx) THEN ! locdbs is rel % diff. btw. measured and simulated
        fitspec(1:ns) = fitspec(1:ns)  * (1.0 - fitvar_rad(i) * locdbs(1:ns))
     ELSE
        fitspec(1:ns) = fitspec(1:ns)  + fitvar_rad(i) * locdbs(1:ns)
     END IF
  ENDDO

  ! Use logarithmic of the reflectance 
  IF (use_lograd)  fitspec(1:ns) = LOG(fitspec(1:ns))

  RETURN
END SUBROUTINE spectra_reflectance
