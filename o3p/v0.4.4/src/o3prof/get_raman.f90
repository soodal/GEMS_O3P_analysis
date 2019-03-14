SUBROUTINE GET_RAMAN(nl, ozprof, errstat)
  
  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY : du2mol, deg2rad, max_ring_pts
  USE OMSAO_variables_module,  ONLY : database, database_shiwf, n_refwvl,   &
       refwvl, the_sca_atm, the_sza_atm,the_vza_atm, refwvl_sav, lo_radbnd, &
       up_radbnd, nsol_ring, n_refwvl_sav, do_bandavg, fitvar_rad,curr_rad_spec,n_rad_wvl,  &
       sol_spec_ring, ozabs_convl, refdbdir, the_aza_atm, n_rad_wvl, refidx
  USE ozprof_data_module,      ONLY : nflay, atmos_prof_fname, nos, mflay,  &
       atmosprof, actawin, aerwavs, gaext, fts, fzs, fozs, frhos, num_iter, wrtring
  USE OMSAO_indices_module,    ONLY : solar_idx, ring_idx, ring1_idx
  USE OMSAO_errstat_module      
  IMPLICIT NONE 

  ! Input/output variables
  INTEGER, INTENT(IN)                       :: nl
  INTEGER, INTENT(OUT)                      :: errstat
  REAL (KIND=dp), DIMENSION(nl), INTENT(IN) :: ozprof
  
  ! Local Variables
  INTEGER, PARAMETER :: maxnu = 15000
  LOGICAL            :: do_o3shi, do_tmpwf, problems, do_bandavg_sav
  INTEGER            :: i, j,ntemp, low, hgh, nnref, fidx, lidx
  REAL (KIND=dp)     :: scl, xg, deltlam
  REAL (KIND=dp), DIMENSION(0:nflay) :: tauin, ts, ozs
  REAL (KIND=dp), DIMENSION(0:nl)    :: cumoz
  REAL (KIND=dp), DIMENSION(nflay)   :: ext
  REAL (KIND=dp), DIMENSION(nsol_ring, nflay) :: strans, vtrans, dads, dadt
  REAL (KIND=dp), DIMENSION(3, nsol_ring)     :: abscrs_qtdepen
  REAL (KIND=dp), DIMENSION(maxnu, nflay)     :: st, vt
  REAL (KIND=dp), DIMENSION(maxnu)            :: ring
  REAL (KIND=dp), DIMENSION(n_refwvl_sav)    :: newring

  ! Saved variables
  INTEGER,                                        SAVE :: nuhi, nulo, nu
  REAL (KIND=dp),                                 SAVE :: avgt, cosvza, cossza
  REAL (KIND=dp), DIMENSION(max_ring_pts),        SAVE :: swavs, raycof
  REAL (KIND=dp), DIMENSION(maxnu),               SAVE :: ramanwav
  REAL (KIND=dp), DIMENSION(max_ring_pts, mflay), SAVE :: abscrs, aerext
  REAL (KIND=dp), DIMENSION(maxnu)                     :: tmpring

  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=9), PARAMETER :: modulename = 'GET_RAMAN'
  
  errstat = pge_errstat_ok 
  
  IF (ozabs_convl) THEN   ! Check this for each cross track position

     ! Get position for raman calculation
     swavs(1:nsol_ring) = sol_spec_ring(1, 1:nsol_ring) ! 264-334
     IF (swavs(1) > refwvl(1)-2.0 .OR. swavs(nsol_ring) < refwvl(n_refwvl) + 2.0) THEN
        WRITE(*, *) modulename, ': Increase wavelength range for solar spectra in ring calculation!!!'
        print * , swavs(1), refwvl(1)-2.0, swavs(nsol_ring), refwvl(n_refwvl)+2.0
        errstat = pge_errstat_error; RETURN
     ENDIF
   
     nuhi = INT(1.0D7 / swavs(1)) ! maximum of frequency, 37745
     nulo = INT(1.0D7 / swavs(nsol_ring)) + 1 ! minimum of frequency, 29977
     nu = 0
     DO i = nulo, nuhi
        nu = nu + 1
        ramanwav(nu) = 1.0D7 / REAL(i, KIND=dp) ! very high resolution
     ENDDO
     ramanwav(nu) = swavs(1); ramanwav(1) = swavs(nsol_ring)
     IF (nuhi - nulo + 1 > maxnu) THEN
        WRITE(*, *) nuhi, nulo, refwvl(1), refwvl(n_refwvl)
        WRITE(*, *) modulename, ': Need to increase maxnu!!!'
        errstat = pge_errstat_error; RETURN
     ELSE IF (nuhi <= nulo) THEN
        WRITE(*, *) modulename, ': nulo>=nuhi, should never happen!!!'
        errstat = pge_errstat_error; RETURN
     ENDIF

     ! Get Rayleigh scattering coefficients once
     CALL GET_ALL_RAYCOF(nsol_ring, swavs(1:nsol_ring), raycof(1:nsol_ring))
     
     ! ozabs_convl will be changed after calling get_abscrs routine
  ENDIF

  IF (num_iter == 0) THEN
     ts(1:nflay) = (fts(1:nflay) + fts(0:nflay-1)) / 2.d0

     ! Use effective temperature to speed up the calculation
     avgt = SUM(ts(1:nflay) * fozs(1:nflay)) / SUM (fozs(1:nflay))

     cossza = COS(the_sza_atm * deg2rad); cosvza = COS(the_vza_atm * deg2rad)

     do_o3shi = .FALSE.; do_tmpwf = .FALSE.
     do_bandavg_sav = do_bandavg; do_bandavg = .FALSE.
     CALL GETABS_CRS(swavs(1:nsol_ring), nsol_ring, nsol_ring, 1, nflay, ts(1:nflay), &
          abscrs(1:nsol_ring, 1:nflay), do_o3shi,  do_tmpwf, dads, dadt, problems, &
          abscrs_qtdepen(1:3, 1:nsol_ring) )

     IF (problems) THEN
        print *, swavs(1), swavs(nsol_ring)
        WRITE(*, *) modulename, ': Problems in reading trace gas absorption!!!'
        do_bandavg = do_bandavg_sav; errstat = pge_errstat_error; RETURN
     ENDIF
     do_bandavg = do_bandavg_sav

     DO i = 1, nsol_ring
        low = COUNT(MASK=(aerwavs(1:actawin) < swavs(i)))
        IF (low == 0) THEN
           low = low + 1
        ELSE IF (low == actawin) THEN 
           low = low - 1
        ENDIF
        hgh = low + 1
        
        xg = (swavs(i) - aerwavs(low)) / (aerwavs(hgh) - aerwavs(low))
        aerext(i, 1:nflay) = gaext(low, 1:nflay) * (1.0 - xg) +  gaext(hgh, 1:nflay) * xg
     ENDDO
  ENDIF
 
  ! Get ozone profile
  IF (num_iter == 0) THEN
     ozs(1:nflay) = fozs(1:nflay)
  ELSE
     IF (nflay /= nl) THEN
        cumoz(0) = 0.0
        DO i = 1, nl
           cumoz(i) = ozprof(i) + cumoz(i-1)
        ENDDO
        
        CALL INTERPOL(atmosprof(2, 0:nl), cumoz, nl+1, fzs(0:nflay), ozs(0:nflay), nflay+1, errstat)
        IF (errstat < 0) THEN
           WRITE(*, *) modulename, ': INTERPOL error, errstat = ', errstat ;  RETURN
        ENDIF
        ozs(1:nflay) = (ozs(1:nflay) - ozs(0:nflay-1))
     ELSE
        ozs(1:nflay) = ozprof(1:nl)
     ENDIF
  ENDIF
  ozs(1:nflay) = ozs(1:nflay) * du2mol  ! convert to molecules

  !WRITE(77, '(2I5, 4D16.7)') nsol_ring, nflay, the_sza_atm, the_vza_atm, the_aza_atm, the_sca_atm
  !WRITE(77, '(100F14.7)') ts(1:nflay)
  !WRITE(77, '(100D16.7)') frhos(1:nflay)
    
     
  ! Compute optical depth
  DO i = 1, nsol_ring
     ext = raycof(i) * frhos(1:nflay) + abscrs(i, 1:nflay) * ozs(1:nflay) + aerext(i, 1:nflay) 
     
     tauin = 0.0
     DO j = 1, nflay
        tauin(j) = tauin(j-1) + ext(j)
     ENDDO
     
     ! Plane parallel
     scl = sol_spec_ring(2, i) * ((swavs(i) / swavs(1)) ** (-4.0))
     strans(i, 1:nflay) = scl * EXP(-tauin(1:nflay) / cossza)
     vtrans(i, 1:nflay) = EXP(-tauin(1:nflay) / cosvza)
     
     !! Spherical geometry (unnecessary)
     !CALL GET_SLANT_TAU (nflay, zs, tauin, sza, tauout)
     !strans(i, 1:nflay) = scl * EXP(-tauout(1:nflay))
     !CALL GET_SLANT_TAU (nflay, zs, tauin, vza, tauout)
     !vtrans(i, 1:nflay) = EXP(-tauout(1:nflay))
     !WRITE(77, '(F12.5, 100D16.7)') swavs(i), sol_spec_ring(2, i), ext(1:nflay)
  ENDDO
  
  ! Interpolate to raman grid in wavenumber
  DO i = 1, nflay
     CALL BSPLINE(swavs(1:nsol_ring), strans(1:nsol_ring, i), nsol_ring, &
          ramanwav(1:nu), st(1:nu, i), nu, errstat)
     IF (errstat < 0) THEN
        WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat; RETURN
     ENDIF
     
     CALL BSPLINE(swavs(1:nsol_ring), vtrans(1:nsol_ring, i), nsol_ring, &
          ramanwav(1:nu), vt(1:nu, i), nu, errstat)
     IF (errstat < 0) THEN
        WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat; RETURN
     ENDIF
  ENDDO
  
  ! Call raman program
  CALL RAMAN(refdbdir, nulo, nuhi, nu, nflay, avgt, the_sca_atm, cossza, &
       st(1:nu,:), vt(1:nu,:), frhos(1:nflay), ring(1:nu))

  !CALL BSPLINE(ramanwav(1:nu), ring(1:nu), nu, swavs(1:nsol_ring), &
  !     tmpring(1:nsol_ring), nsol_ring, errstat) 
  !WRITE(78, '(F12.5, D16.7)') ((swavs(i), tmpring(i)), i = 1, nsol_ring)
  
  ! Interpolate calculated ring back to gome radiance grids
  CALL BSPLINE(ramanwav(1:nu), ring(1:nu), nu, refwvl_sav(1:n_refwvl_sav), &
       newring(1:n_refwvl_sav), n_refwvl_sav, errstat) 

  
  IF (errstat < 0) THEN
     WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat; RETURN
  ENDIF

  IF (do_bandavg) THEN
     CALL avg_band_refspec(refwvl_sav(1:n_refwvl_sav), newring(1:n_refwvl_sav), &
          n_refwvl_sav, ntemp, errstat)
     IF ( errstat /= 0 .OR. ntemp /= n_refwvl) THEN
        WRITE(*, *) modulename, ': Ring Spectra Averaging Error ', n_refwvl_sav, ntemp
        errstat = pge_errstat_error; RETURN
     ENDIF
  ENDIF
  
  ! put the Ring spectrum into database
  database(ring_idx,  1:n_refwvl)  = newring(1:n_refwvl)
  !database(ring1_idx,  1:n_refwvl) = newring(1:n_refwvl)

  IF (wrtring) THEN
     WRITE(92, *) n_rad_wvl
     DO i = 1, n_rad_wvl
        WRITE(92, '(f8.4, D14.5)') refwvl(refidx(i)), newring(refidx(i))
     ENDDO
  ENDIF
  RETURN
  
END SUBROUTINE GET_RAMAN



