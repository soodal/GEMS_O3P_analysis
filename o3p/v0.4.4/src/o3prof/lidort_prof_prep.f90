! ====================================================================
!  Author: Xiong Liu
!  Date: Jan. 30, 2004
!  Purpose: Prepare input (tau, single scattering albedo, phase
!		moments, their variation, surface BDRF, spherical attenutaion
!		for LIDORT calculation
! ====================================================================
! Reminder:
! Lambertian Surface model is implemented
! Refractive atmosphere is not implemented
! Modication History:
! ====================================================================

! Description of Auguments
! lamda:   wavelength 
! zsgrid:  altitude at each level from TOA to BOS in km, nlayers+1 levels
! airgrid: air column density for each layer, molecules / cm^2 
! varyprof:arrays of linearization flags for each layer, 0:no 1:yes

! Currently, number of gases to be allowed is one, i.e., O3
! Number of gases possibly to be considered later includes
! 1: O3 2: NO2 3: O4 4: BrO 5: SO2 6: HCHO  7: OCLO  8: O2  9:H2O
! These species will be added when needed, just need to add corresponding
! cross section database and read them
 
! ngas  :  number of gases
! gasin :  pointer to gases that are used
! abscrs:  Input/output
!		   If get_crs is set, then it refers to absorption cross section 
!		   at each layer and for each species for each molecule 
!          On return, it gives the absorption od for each species at each layer
! gascol : column density for each species at each layer, molecules/cm^2

! useasy    :  use asymmetric factor/phase moments for clouds/aerosols
! nmoms     :  number of phase moments
! do_aerosols: include aerosols
! aersca    :  aerosol scattering coefficients at each layer
! aerext    :  aerosol extinction coefficients at each layer
! aerasy    :  aerosol asymmetric factor  at each layer
! aermsk    :  Indicator of aerosols for each layer, 1: with aerosol, 0: no aerosols
! aermoms   :  Aerosol moments at layers with aerosols

! do_clouds : include clouds 
! cldsca    :  cloud scattering coefficients at each layer
! cldext    :  cloud extinction coefficients at each layer
! cldasy    :  cloud asymmetric factor  at each layer
! cldmsk    :  Indicator of clouds for each layer, 1: with cloud, 0: no clouds
! cldmoms   :  cloud moments at layers with clouds

SUBROUTINE LIDORT_PROF_PREP (lamda, raycof, depol, zsgrid, airgrid,  varyprof, &
     ngas, gasin, abscrs, gascol, eta, useasy, nmoms, &          
     do_aerosols, aersca, aerext, aerasy, aermoms, aermsk, &
     do_clouds, cldsca, cldext, cldasy, cldmoms, cldmsk, problems, &
     deltau, delsca, delo3abs, delray)

  USE OMSAO_precision_module
  USE ozprof_data_module, ONLY : maxgksec, maxgkmatc, ngksec, ngkmatc
  IMPLICIT NONE  
   
  !===============================  Define Variables ===========================
  ! Include files of dimensions and numbers
  INCLUDE 'VLIDORT.PARS'  
  ! Include files of input variables
  INCLUDE 'VLIDORT_INPUTS.VARS'
  INCLUDE 'VLIDORT_SETUPS.VARS'
!  INCLUDE 'VLIDORT_REFLECTANCE.VARS'
  INCLUDE 'VLIDORT_L_INPUTS.VARS'
  INCLUDE 'VLIDORT_BOOKKEEP.VARS'

  ! Input variables
  INTEGER, INTENT(IN)                     :: ngas, nmoms
  INTEGER, INTENT(IN), DIMENSION(ngas)    :: gasin
  LOGICAL, INTENT(IN), DIMENSION(nlayers) :: cldmsk, aermsk, varyprof
  LOGICAL, INTENT(IN)                     :: useasy

  REAL (KIND=dp), INTENT(IN)              :: raycof, depol, lamda
  REAL (KIND=dp), DIMENSION(0:nlayers), INTENT(IN) :: zsgrid
  REAL (KIND=dp), DIMENSION(nlayers), INTENT(IN)   :: airgrid,  &
       aersca, aerext, aerasy, cldsca, cldext, cldasy
  REAL (KIND=dp), DIMENSION(0:nmoms, maxgksec, nlayers), INTENT(IN) :: aermoms
  REAL (KIND=dp), DIMENSION(0:nmoms, maxgksec, nlayers), INTENT(IN) :: cldmoms
  REAL (KIND=dp), DIMENSION(ngas, nlayers), INTENT(IN)  :: abscrs, gascol, eta

  ! Optional output
  REAL (KIND=dp), DIMENSION(nlayers), INTENT(OUT) :: deltau, delsca, delo3abs, delray
  
  
  ! Output variables
  LOGICAL, INTENT(OUT)   :: problems

  ! Modified variables
  LOGICAL, INTENT(INOUT) :: do_aerosols, do_clouds

  ! Local variables
  INTEGER, PARAMETER     :: maxngas = 7, maxscatter=3, allngas = 9
  INTEGER, DIMENSION(maxgkmatc), PARAMETER :: &
       greekmat_idxs = (/1, 2, 5, 6, 11, 12, 15, 16/), phasmoms_idxs = (/1, 5, 5, 2, 3, 6, 6, 4/)

  INTEGER					  :: ui, i, j, k, q, nscatter, idx, cldidx, aeridx, nactgksec, nactgkmatc
  INTEGER, DIMENSION(allngas) :: absin
  REAL (KIND=dp)		      :: scaco_r, absco_r, &  ! raycof, depol, 
		extco_r, extco, scaco, pvar, extco_a, scaco_a, extco_c, scaco_c, j0, j1
  REAL (KIND=dp), DIMENSION(maxscatter)              :: scaco_input
  REAL (KIND=dp), DIMENSION(nlayers)                 :: extconf
  REAL (KIND=dp), DIMENSION(ngas, nlayers)           :: absod
  REAL (KIND=dp), DIMENSION(0:maxmoments_input, 1:maxgksec, maxscatter),     SAVE :: phasmoms_input
  REAL (KIND=dp), DIMENSION(0:maxmoments_input, 1:maxgksec),                 SAVE :: phasmoms_total_input
  REAL (KIND=dp), DIMENSION(0:max_atmoswfs, 0:maxmoments_input, 1:maxgksec), SAVE :: l_phasmoms_total_input
  LOGICAL, SAVE                                                             :: first = .TRUE.
  
  ! ========================== Check for Input ==================================
  problems = .FALSE.
  !IF (MAXVAL(gasin) > allngas) THEN
  !   WRITE(*, *) 'Not weighting functions are implemented for all gases!!!'
  !   problems = .TRUE.; RETURN
  !ENDIF

  IF (first) THEN
     IF (.NOT. useasy .AND. nmoms > maxmoments_input) THEN
        WRITE(*, *) 'Need to increase maxmoments_input for aerosols/clouds!!!'
        problems = .TRUE.; RETURN
     ENDIF
     
     IF (.NOT. useasy .AND. nmoms < ngreek_moments_input) THEN
        WRITE(*, *) 'Need to increase input moments for aerosols/clouds!!!'
        problems = .TRUE.; RETURN
     ENDIF

     ! This only needs to be initialized once
     phasmoms_input        = ZERO
     phasmoms_total_input  = ZERO
     greekmat_total_input  = ZERO
				   
     l_deltau_vert_input   = ZERO
     l_omega_total_input   = ZERO
     l_greekmat_total_input= ZERO 
     l_phasmoms_total_input= ZERO
     
     first =.FALSE.
  ENDIF
  
  IF (NSTOKES == 1) THEN
     nactgksec = 1;  nactgkmatc = 1
  ELSE
     nactgksec = ngksec; nactgkmatc = ngkmatc
  ENDIF

  !WRITE(*, *) nmoms, nmoments, ngreek_moments_input, maxmoments 
  !WRITE(*, *) SUM(gascol), SUM(airgrid)
  !WRITE(*, *) varyprof(1), varyprof(nlayers)
  !WRITE(*, *) abscrs(1, 1), abscrs(1, 30)

  ! 1: O3 2: NO2  3:O2  4: O4 5: BrO 6: H2O 7 SO2 8: HCHO  9: OCLO
  absin(:) = 0
  DO i =1, ngas
		absin(gasin(i)) = i
  ENDDO

  ! Disable clouds and aerosols if for Rayleigh scattering atmosph:q!ere
  !IF (do_rayleigh_only) THEN
  !   do_clouds = .FALSE.; do_aerosols = .FALSE.
  !ENDIF

  ! Enable delta-M scaling for clouds or aerosols
  !IF (do_clouds .OR. do_aerosols) do_deltam_scaling = .TRUE.

  ! Start layer loop
  taugrid_input(0) = ZERO 

  ! Get rayleigh scattering phase function moments (Same for each layer)
  ! unassigned elements have already initialized to zero
  phasmoms_input(0, 1, 1) = ONE
  phasmoms_input(2, 1, 1) = (ONE - depol) / (TWO + depol)  
  IF (nactgksec == 6) THEN
     phasmoms_input(2, 2, 1) = 6.0D0 * phasmoms_input(2, 1, 1)
     phasmoms_input(2, 5, 1) = -SQRT(6.0D0) * phasmoms_input(2, 1, 1)
     phasmoms_input(1, 4, 1) = 3.0D0 * (ONE - 2.0D0 * depol) / (TWO + depol)
  ENDIF
 
  DO i = 1, nlayers   
     ! Rayleigh scattering
     scaco_r = raycof * airgrid(i)
     delray(i) = scaco_r

     ! Gas absorption
     absod(1:ngas, i) = abscrs(1:ngas, i) * gascol(1:ngas, i)
     absco_r = SUM(absod(1:ngas, i))   

     IF (absco_r <= 0.0) THEN
        problems  = .TRUE.
        print *, 'Negative total absorption: ', lamda, i, ngas
        print *, abscrs(1:ngas, i)
        print *, gascol(1:ngas, i)       
        RETURN
     ENDIF
     extco_r = absco_r + scaco_r
     scaco_input(1) = scaco_r
     
     ! Aerosols and clouds
     extco = extco_r
     nscatter = 1; cldidx = 0; aeridx = 0
     extco_a = ZERO; scaco_a = ZERO; extco_c = ZERO; scaco_c = ZERO; 

     !IF (.NOT. do_rayleigh_only) THEN        
        IF (do_clouds .AND. cldmsk(i)) THEN
           nscatter = nscatter + 1
           extco = extco + cldext(i)
           extco_c = cldext(i);   scaco_c = cldsca(i)
           scaco_input(nscatter) = cldsca(i)
           cldidx = nscatter
           
           ! get phase moments for clouds
           IF (.NOT. useasy) THEN
              phasmoms_input(0:nmoms, 1:nactgksec, nscatter) = cldmoms(0:nmoms, 1:nactgksec, i)                         
           ELSE ! use H_G function
              phasmoms_input(0, 2, nscatter) = ONE
              j0 = ONE             
              DO j = 1, ngreek_moments_input
                 j1 = REAL(2*j+1, KIND=dp)
                 phasmoms_input(j, 2, nscatter) = (j1/j0) * cldasy(i) * phasmoms_input(j-1, 2, nscatter)
                 j0 = j1
              ENDDO
           ENDIF
           
        ENDIF  ! end clouds
        
        IF (do_aerosols .AND. aermsk(i)) THEN
           nscatter = nscatter + 1
           extco = extco + aerext(i)
           scaco_input(nscatter) = aersca(i)
           extco_a = aerext(i);  scaco_a = aersca(i)
           aeridx = nscatter

           ! get phase moments for aerosols
           IF (.NOT. useasy) THEN  
              phasmoms_input(0:nmoms, 1:nactgksec, nscatter) = aermoms(0:nmoms, 1:nactgksec, i)  
           ELSE ! use H_G function
              phasmoms_input(0, 2, nscatter) = ONE
              j0 = ONE             
              DO j = 1, ngreek_moments_input
                 j1 = REAL(2*j+1, KIND=dp)
                 phasmoms_input(j, 2, nscatter) = (j1/j0) * aerasy(i) * phasmoms_input(j-1, 2, nscatter)
                 j0 = j1
              ENDDO
           ENDIF
        ENDIF  ! end aerosols
     !ENDIF     ! end non-rayleigh
     
     ! setup LIDORT input for tau and omega
     scaco = SUM(scaco_input(1:nscatter))
     omega_total_input(i) = scaco / extco
     
     IF (omega_total_input(i) < OMEGA_SMALLNUM) omega_total_input(i) = OMEGA_SMALLNUM 
     IF (omega_total_input(i) > 1.0 - OMEGA_SMALLNUM)  &
          omega_total_input(i) = 1.0 - OMEGA_SMALLNUM
     taugrid_input(i) = taugrid_input(i-1) + extco
     deltau_vert_input(i) = extco
     !extconf(i) = extco / (zsgrid(i-1) - zsgrid(i))   ! extinction coefficients
     !IF (i > 25) THEN
     !   print *, i, aermsk(i), taugrid_input(i), omega_total_input(i)
     !   print *, extco, absco_r, scaco_r, extco_a, scaco_a
     !ENDIF
  
     ! sum up phase moments as required in LIDORT
     DO j = 0, ngreek_moments_input
        DO k = 1, nactgksec
           phasmoms_total_input(j, k) = SUM(phasmoms_input(j, k, 1:nscatter) &
                * scaco_input(1:nscatter)) / scaco
        ENDDO
     ENDDO
     !phasmoms_total_input(ngreek_moments_input+1:maxmoments, 1:maxgksec) = 0.0  
     
     ! Set up greek scattering matrix for each moment 
     !greekmat_total_input(0:ngreek_moments_input, i, 1)  = phasmoms_total_input(0:ngreek_moments_input, 1)
     !greekmat_total_input(0:ngreek_moments_input, i, 2)  = phasmoms_total_input(0:ngreek_moments_input, 5)
     !greekmat_total_input(0:ngreek_moments_input, i, 5)  = phasmoms_total_input(0:ngreek_moments_input, 5)
     !greekmat_total_input(0:ngreek_moments_input, i, 6)  = phasmoms_total_input(0:ngreek_moments_input, 2)
     !greekmat_total_input(0:ngreek_moments_input, i, 11) = phasmoms_total_input(0:ngreek_moments_input, 3)
     !greekmat_total_input(0:ngreek_moments_input, i, 12) = phasmoms_total_input(0:ngreek_moments_input, 6)
     !greekmat_total_input(0:ngreek_moments_input, i, 15) = -phasmoms_total_input(0:ngreek_moments_input, 6)
     !greekmat_total_input(0:ngreek_moments_input, i, 16) = phasmoms_total_input(0:ngreek_moments_input, 4)
     !greekmat_total_input(ngreek_moments_input+1:maxmoments, i, 1:MAXSTOKES_SQ) = 0.0
     greekmat_total_input(0:ngreek_moments_input, i, greekmat_idxs(1:nactgkmatc)) = &
          phasmoms_total_input(0:ngreek_moments_input, phasmoms_idxs(1:nactgkmatc))
     IF ( nactgkmatc > 1 ) greekmat_total_input(0:ngreek_moments_input, i, 15) &
          = -greekmat_total_input(0:ngreek_moments_input, i, 15)

     ! This should always be 1, but may be slightly different due to numerical truncation
     greekmat_total_input(0, i, 1) = 1.0  

     !IF  (i == 26) THEN
     !WRITE (91, '(I5, 7D24.12, I5)') i, extco, scaco, scaco_a, absco_r, scaco_r, &
     !     omega_total_input(i), greekmat_total_input(0, i, 1), nscatter
     !DO k = 1, 16
     !   WRITE (91, '(1000D24.12)') (greekmat_total_input(j, i, k), j=0, ngreek_moments_input)
     !ENDDO
     !ENDIF

     !IF (do_simulation_only .OR. .NOT. varyprof(i)) THEN   ! no linearition
     !   ! zero out quantity for safety
     !   layer_vary_flag(i) = .FALSE.
     !   layer_vary_number(i) = 0	
	 !   	   
     !   !l_deltau_vert_input(:, i) =   ZERO
     !   !l_omega_total_input(:, i) = ZERO
     !   !l_greekmat_total_input(:, : , :, i) = ZERO        
     !
     !ELSE
     !   layer_vary_flag(i) = .TRUE.
     !   layer_vary_number(i) = n_totalatmos_wfs
     !The above part has been taken care of in routine: lidort_prof_env.f90
        
     DO q = 1, n_totalatmos_wfs        
        !  w.r.t ozone volume mixing ratio: 1
        !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IF ( profilewf_names(q) == 'ozone volume mixing ratio------' ) THEN
           idx = absin(1)
           IF (idx < 1) THEN
              WRITE(*, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.; RETURN
           ENDIF
           l_omega_total_input(q, i) = - absod(idx, i) / extco
           l_deltau_vert_input(q, i) = + absod(idx, i) / extco
           !l_greekmat_total_input(q, 0:maxmoments , i, 1:16) = ZERO
           
           !  w.r.t NO2 volume mixing ratio: 2
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ELSE IF ( profilewf_names(q) == 'NO2 volume mixing ratio------' ) THEN
           idx = absin(2)
           IF (idx < 1) THEN
              WRITE(*, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.; RETURN
           ENDIF
           l_omega_total_input(q, i) = - absod(idx, i) / extco
           l_deltau_vert_input(q, i) = + absod(idx, i) / extco
           !l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t O2 volume mixing ratio: 8
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ELSE IF ( profilewf_names(q) == 'O2 volume mixing ratio------' ) THEN
           idx = absin(8)
           IF (idx < 1) THEN
              WRITE(*, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.; RETURN
           ENDIF
           l_omega_total_input(q, i) = - absod(idx, i) / extco
           l_deltau_vert_input(q, i) = + absod(idx, i) / extco
           !l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t O4 volume mixing ratio: 3
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ELSE IF ( profilewf_names(q) == 'O4 volume mixing ratio------' ) THEN
           idx = absin(3)
           IF (idx < 1) THEN
              WRITE(*, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.;  RETURN
           ENDIF
           l_omega_total_input(q, i) = - absod(idx, i) / extco
           l_deltau_vert_input(q, i) = + absod(idx, i) / extco
           !l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t BrO volume mixing ratio: 4
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ELSE IF ( profilewf_names(q) == 'BrO volume mixing ratio------' ) THEN
           idx = absin(4)
           IF (idx < 1) THEN
              WRITE(*, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.; RETURN
           ENDIF
           l_omega_total_input(q, i) = - absod(idx, i) / extco
           l_deltau_vert_input(q, i) = + absod(idx, i) / extco
           !l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t H2O volume mixing ratio: 9
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ELSE IF ( profilewf_names(q) == 'H2O volume mixing ratio------' ) THEN
           idx = absin(9)
           IF (idx < 1) THEN
              WRITE(*, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.; RETURN
           ENDIF
           l_omega_total_input(q, i) = - absod(idx, i) / extco
           l_deltau_vert_input(q, i)         = + absod(idx, i) / extco
           !l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t SO2 volume mixing ratio: 5
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ELSE IF ( profilewf_names(q) == 'SO2 volume mixing ratio------' ) THEN
           idx = absin(5)
           IF (idx < 1) THEN
              WRITE(*, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.; RETURN
           ENDIF
           l_omega_total_input(q, i) = - absod(idx, i) / extco
           l_deltau_vert_input(q, i) = + absod(idx, i) / extco
           !l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t HCHO volume mixing ratio: 6
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ELSE IF ( profilewf_names(q) == 'HCHO volume mixing ratio------' ) THEN
           idx = absin(6)
           IF (idx < 1) THEN
              WRITE(*, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.; RETURN
           ENDIF
           l_omega_total_input(q, i) = - absod(idx, i) / extco
           l_deltau_vert_input(q, i)         = + absod(idx, i) / extco
           !l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t OCLO volume mixing ratio: 7
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ELSE IF ( profilewf_names(q) == 'OCLO volume mixing ratio------' ) THEN
           idx = absin(7)
           IF (idx < 1) THEN
              WRITE(*, *) idx, 'This gas is not modeled. No WF can be done!!!'
              problems = .TRUE.; RETURN
           ENDIF
           l_omega_total_input(q, i) = - absod(idx, i) / extco
           l_deltau_vert_input(q, i)         = + absod(idx, i) / extco
           !l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t average temperature of layer
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           !  no variation of phase functions
           !  Assume no effects on air density              
        ELSE IF ( profilewf_names(q) == 'average temperature of layer---' ) THEN
           l_omega_total_input(q, i) = - SUM(absod(:, i) * eta(:, i)) / extco
           l_deltau_vert_input(q, i) = + SUM(absod(:, i) * eta(:, i)) / extco
           !l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t average pressure of layer
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           !  no variation of phase functions
        ELSE IF ( profilewf_names(q) == 'average pressure of layer------' ) THEN
           
           pvar = extco_r/ extco
           l_omega_total_input(q, i) = ((ONE - pvar) * scaco_input(1) - &
                pvar * (scaco - scaco_input(1))) / scaco
           l_deltau_vert_input(q,i) = extco_r / extco
           !l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t rayleigh soptical thickness
           ! xliu: August 12, 2008
           !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
        ELSE IF ( profilewf_names(q) == 'rayleigh optical thickness-----' ) THEN
           pvar = scaco_r / extco	       
           l_omega_total_input(q,i) = (1.0 - omega_total_input(i)) * scaco_r / extco / omega_total_input(i)
           l_deltau_vert_input(q,i) = pvar
           l_greekmat_total_input(q, 0:ngreek_moments_input, i, :) = ZERO
           DO j = 0, ngreek_moments_input
             DO k = 1, nactgksec
                IF (phasmoms_total_input(j, k) /= 0.0) THEN
                   l_phasmoms_total_input(q, j, k) = ( phasmoms_input(j, k, 1) - phasmoms_total_input(j, k) ) &
                        / phasmoms_total_input(j, k) * scaco_a / scaco
                ELSE
                   l_phasmoms_total_input(q, j, k) = 0.0
                ENDIF
             ENDDO
          ENDDO
          l_greekmat_total_input(q, 0:ngreek_moments_input, i, greekmat_idxs(1:nactgkmatc)) = &
               l_phasmoms_total_input(q, 0:ngreek_moments_input, phasmoms_idxs(1:nactgkmatc))
          IF ( nactgkmatc > 1 )  l_greekmat_total_input(q, 0:ngreek_moments_input, i, 15) &
               = - l_greekmat_total_input(q, 0:ngreek_moments_input, i, 15)          
        ELSE IF ( profilewf_names(q) == 'rayleigh scattering coefficient' ) THEN
           !  xliu: April 13, 2007 
           !  Still need to consider the variation in phase function  
           pvar = scaco_r / extco	       
           l_omega_total_input(q, i) = ((ONE - pvar) * scaco_input(1) - &
                pvar * (scaco - scaco_input(1)) ) / scaco
           l_omega_total_input(q,i) = pvar
           l_deltau_vert_input(q,i)     = scaco_r / extco
           !l_greekmat_total_input(q, : , i, :) = ZERO
           
           !  w.r.t aerosol extinction coefficient / aerosol optical thickness
           !  aerosol scattering albedo does not change   
           !  xliu: April 13, 2007 (consider the variation in phase function)
        ELSE IF ( profilewf_names(q) == 'aerosol extinction coefficient-' ) THEN
           IF (aeridx > 0) THEN
              l_deltau_vert_input(q,i) = + extco_a / extco
              l_omega_total_input(q,i) = (scaco_a  / extco_a - omega_total_input(i) ) &
                   / omega_total_input(i)  * extco_a / extco
              l_greekmat_total_input(q, 0:ngreek_moments_input, i, :) = ZERO
              DO j = 0, ngreek_moments_input
                 DO k = 1, nactgksec
                    IF (phasmoms_total_input(j, k) /= 0.0) THEN
                       l_phasmoms_total_input(q, j, k) = ( phasmoms_input(j, k, aeridx) &
                            - phasmoms_total_input(j, k) ) / phasmoms_total_input(j, k) &
                            * scaco_a / scaco
                    ELSE
                       l_phasmoms_total_input(q, j, k) = 0.0
                    ENDIF
                 ENDDO
              ENDDO
              l_greekmat_total_input(q, 0:ngreek_moments_input, i, greekmat_idxs(1:nactgkmatc)) = &
                   l_phasmoms_total_input(q, 0:ngreek_moments_input, phasmoms_idxs(1:nactgkmatc))
              IF ( nactgkmatc > 1 )  l_greekmat_total_input(q, 0:ngreek_moments_input, i, 15) &
                   = - l_greekmat_total_input(q, 0:ngreek_moments_input, i, 15)
              !print *, maxval(l_greekmat_total_input), minval(l_greekmat_total_input)
              !print *, i, l_deltau_vert_input(q,i), l_omega_total_input(q,i)
              !WRITE(*, '(6D14.6)') (l_phasmoms_total_input(1, j, 1:nactgksec), j = 0, ngreek_moments_input)
              !STOP
           ENDIF
           !  w.r.t  aerosol scattering coefficient / single scattering albedo
           !  aerosol optical thickness will not change
        ELSE IF ( profilewf_names(q) == 'aerosol scattering coefficient-' ) THEN
           IF (aeridx > 0) THEN
              l_deltau_vert_input(q,i) = ZERO
              l_omega_total_input(q,i) = scaco_a  / omega_total_input(i) / extco
              l_greekmat_total_input(q, 0:ngreek_moments_input, i, :) = ZERO
              DO j = 0, ngreek_moments_input
                 DO k = 1, nactgksec
                    IF (phasmoms_total_input(j, k) /= 0.0) THEN
                       l_phasmoms_total_input(q, j, k) = ( phasmoms_input(j, k, aeridx) &
                            - phasmoms_total_input(j, k) ) / phasmoms_total_input(j, k) &
                            * scaco_a / scaco
                    ELSE
                       l_phasmoms_total_input(q, j, k) = 0.0
                    ENDIF
                 ENDDO
              ENDDO
              l_greekmat_total_input(q, 0:ngreek_moments_input, i, greekmat_idxs(1:nactgkmatc)) = &
                   l_phasmoms_total_input(q, 0:ngreek_moments_input, phasmoms_idxs(1:nactgkmatc))
              IF ( nactgkmatc > 1 )  l_greekmat_total_input(q, 0:ngreek_moments_input, i, 15) &
                   = - l_greekmat_total_input(q, 0:ngreek_moments_input, i, 15)
           ENDIF
           !   w.r.t cloud extinction coefficient / optical thickness
        ELSE IF ( profilewf_names(q) == 'cloud extinction coefficient---' ) THEN
           IF (cldidx > 0) THEN
              l_deltau_vert_input(q,i) = + extco_c / extco
              l_omega_total_input(q,i) = (scaco_c  / extco_c - omega_total_input(i) ) &
                   / omega_total_input(i)  * extco_c / extco
              l_greekmat_total_input(q, 0:ngreek_moments_input, i, :) = ZERO
              DO j = 0, ngreek_moments_input
                 DO k = 1, nactgksec
                    IF (phasmoms_total_input(j, k) /= 0.0) THEN
                       l_phasmoms_total_input(q, j, k) = ( phasmoms_input(j, k, cldidx) &
                            - phasmoms_total_input(j, k) ) / phasmoms_total_input(j, k) &
                            * scaco_a / scaco
                    ELSE
                       l_phasmoms_total_input(q, j, k) = 0.0
                    ENDIF
                 ENDDO
              ENDDO
              l_greekmat_total_input(q, 0:ngreek_moments_input, i, greekmat_idxs(1:nactgkmatc)) = &
                   l_phasmoms_total_input(0, 1:ngreek_moments_input, phasmoms_idxs(1:nactgkmatc))
              IF ( nactgkmatc > 1 )  l_greekmat_total_input(q, 0:ngreek_moments_input, i, 15) &
                   = - l_greekmat_total_input(q, 0:ngreek_moments_input, i, 15)
           ENDIF
           
           !  w.r.t clouds scattering coefficient
        ELSE IF ( profilewf_names(q) == 'cloud scattering coefficient---' ) THEN
           IF (cldidx > 0) THEN
              l_deltau_vert_input(q,i) = ZERO
              l_omega_total_input(q,i) = scaco_c  / omega_total_input(i) / extco
              l_greekmat_total_input(q, 0:ngreek_moments_input, i, :) = ZERO
              DO j = 0, ngreek_moments_input
                 DO k = 1, nactgksec
                    IF (phasmoms_total_input(j, k) /= 0.0) THEN
                       l_phasmoms_total_input(q, j, k) = ( phasmoms_input(j, k, cldidx) &
                            - phasmoms_total_input(j, k) ) / phasmoms_total_input(j, k) &
                            * scaco_a / scaco
                    ELSE
                       l_phasmoms_total_input(q, j, k) = 0.0
                    ENDIF
                 ENDDO
              ENDDO
              l_greekmat_total_input(q, 0:ngreek_moments_input, i, greekmat_idxs(1:nactgkmatc)) = &
                   l_phasmoms_total_input(q, 0:ngreek_moments_input, phasmoms_idxs(1:nactgkmatc))
              IF ( nactgkmatc > 1 )  l_greekmat_total_input(q, 0:ngreek_moments_input, i, 15) &
                   = - l_greekmat_total_input(q, 0:ngreek_moments_input, i, 15)
           ENDIF
        ENDIF	     ! end selection of weighting function 
     ENDDO           ! n_totalatmos_wfs loop
  !ENDIF             ! end of do_linearization    
  ENDDO              ! layer loop
  
  deltau(1:nlayers) = deltau_vert_input(1:nlayers)
  delsca(1:nlayers) = deltau_vert_input(1:nlayers) * omega_total_input(1:nlayers)
  delo3abs(1:nlayers) = absod(1, 1:nlayers)
    
    
!  ! Prepare for surface albedo  (unnecessary)
!  IF (do_lambertian_surface) THEN
!     DO i = 1, nstreams
!        bireflec_0 (0, i, 1) = ONE
!        emissivity (i) = ONE - lambertian_albedo
!        DO j = 1, nstreams
!           bireflec (0, i, j) = ONE
!        ENDDO
!     ENDDO
!     DO ui = 1, n_user_streams
!        user_bireflec_0 (0, ui, 1) = ONE
!        user_emissivity (ui) = ONE - lambertian_albedo
!        DO j = 1, nstreams
!           user_bireflec (0, ui, j) = ONE
!        ENDDO
!     ENDDO
!  ENDIF
!   
!  ! Prepare for spherical attenuation
!  IF ( .NOT. do_chapman_function ) THEN
!     CALL PREPARE_SPHERICAL (nlayers, do_plane_parallel, beam_szas(1), earth_radius, &
!          extconf, zsgrid, deltau_slant_input(1:nlayers, 1:nlayers, 1), &
!          sza_local_input(1:nlayers, 1))
!  ENDIF
  
END SUBROUTINE LIDORT_PROF_PREP

! ===============================================================
! Modified from provided routine in LODORT V23 by Rob (F77->F90)
! This routine is not consistent with CHAPMAN_FUNCTION, maybe sth.
! is wrong
! ===============================================================

SUBROUTINE PREPARE_SPHERICAL (nlayers, do_plane_parallel,  input_sunzen, &
     re, ext, z_grid, tauthick_input, sunlocal_input )
  
  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY : deg2rad
  IMPLICIT NONE
  
  !  Input arguments
  INTEGER, INTENT(IN)		   :: nlayers
  LOGICAL, INTENT(IN)		   :: do_plane_parallel
  REAL (KIND=dp), INTENT(IN)   :: input_sunzen, re
  REAL (KIND=dp), DIMENSION (0:nlayers), INTENT(IN) :: z_grid
  REAL (KIND=dp), DIMENSION (nlayers), INTENT(IN)   :: ext
  
  !  Output
  REAL (KIND=dp), DIMENSION (nlayers, nlayers), INTENT(OUT) :: tauthick_input
  REAL (KIND=dp), DIMENSION (nlayers), INTENT(OUT)          :: sunlocal_input
  
  !  local variables
  INTEGER        :: n, j, m
  REAL (KIND=dp) :: gm_toa, th_toa, th0, th1, gm0, gm1
  REAL (KIND=dp) :: h(0:nlayers), delz, taup, mu_toa
  REAL (KIND=dp) :: z_dIff, z_0, z, const0
  REAL (KIND=dp) :: x, xd, cumdep, s2, hf, deltm, delt
  

  !  get spherical optical depths
  !  ----------------------------
  !  prepare spherical attenuation (shell geometry)
  
  IF ( .NOT. do_plane_parallel ) THEN
     
     mu_toa = COS ( input_sunzen * deg2rad )
     gm_toa = SQRT ( 1.0d0 - mu_toa * mu_toa )
     th_toa = ASIN (gm_toa)
     h(0:nlayers) = z_grid(0:nlayers) + re
     const0 = gm_toa / h(0)
     cumdep = 0.0d0
     
     DO n = 1, nlayers
	    delz = z_grid(n-1)-z_grid(n)
	    z_diff   = delz
	    z_0 = z_grid(n-1)
	    z  = z_0
	    DO j = 1, 1
           z = z - z_dIFf
           x = z_0 - z
           xd = x + cumdep
           hf = h(0) - xd
           gm0 = const0 * hf	  
           th0 = ASIN ( gm0 )
           taup = 0.0d0
           DO m = 1, n-1
              gm1 = h(m-1) * gm0 / h(m)
              th1 = ASIN(gm1)
              s2 = h(m-1) * SIN(th1-th0) / gm1
              IF ( j == 1 ) tauthick_input(n,m) = ext(m) * s2
              taup = taup + ext(m) * s2
              th0 = th1
              gm0 = gm1        
           ENDDO
           s2 = h(n-1)* SIN(th_toa-th0) / gm_toa
           taup = taup + ext(n) * s2
           IF ( j== 1 ) tauthick_input(n,n) = ext(n) * s2
	    ENDDO
	    cumdep = xd
	    sunlocal_input(n) = mu_toa
     ENDDO
  ELSE
     
     mu_toa = COS( input_sunzen * deg2rad )
     DO n = 1, nlayers
	    sunlocal_input(n) = mu_toa
	    delz = z_grid(n-1)-z_grid(n)
	    delt = ext(n) * delz
	    tauthick_input(n,n) = delt/mu_toa
	    DO m = 1, n-1
           deltm = ext(m) * delz
           tauthick_input(n,m) = deltm / mu_toa
	    ENDDO
     ENDDO
     
  ENDIF

END SUBROUTINE PREPARE_SPHERICAL


! ==============================================================
!	Routine to calculate Rayleight scattering coefficients
!	  and molecular depolarization factor
!   LAMDA: Wavelength in nm
!   RAYCOF:Rayleight scattering coefficients cm-2/molecule
!   DEPOL: Molecular depolarization ratio
! ==============================================================

!SUBROUTINE GET_RAYCOF_DEPOL(lamda, raycof, depol)
!  
!  USE OMSAO_precision_module
!  IMPLICIT NONE
! 
!  !	Input/Output
!  REAL (KIND=dp), INTENT(IN)  :: lamda
!  REAL (KIND=dp), INTENT(OUT) :: raycof, depol
!  
!  !	Local variables
!  REAL (KIND=dp) :: sig, sig2, sig2p, sig4, fk_n2, fk_o2, fking
!  REAL (KIND=dp), PARAMETER :: abod = 1.0455996d0, bbod = -341.29061d0, &
!       cbod = -0.90230850d0, dbod = 0.0027059889d0, ebod = -85.968563d0
! 
!  !	Rayleigh coefficient
!  ! Using bodhaine et al, j. atm. oceanic tech. 16, 1854-1861, 1999.
!  sig =    1.0d3 / lamda
!  sig2 =   sig * sig
!  sig2p =  1.d0 / sig2
!  sig4 =   sig2 * sig2
!  raycof = (abod + bbod * sig2 + cbod * sig2p) &
!       / (1.d0 + dbod * sig2 + ebod * sig2p) * 1.d-28
!  
!  !	Derivation of depolarization factor d from king factors for air,
!  !	fking = (6 + 3.depol) / (6 - 7.depol) 
!  !	bodhaine et al., 370 ppmv co2
!  fk_n2 = 1.034d0 + 3.17d-4 * sig2
!  fk_o2 = 1.096d0 + 1.385d-3 * sig2 + 1.448d-4 * sig4
!  !fk_ar = 1.00
!  !fk_co2 = 1.15d0
!  !fking = (78.084d0 * fk_n2 + 20.946d0 * fk_o2 + 0.934d0 * fk_ar + &
!  !     0.037d0 * fk_co2) / (78.084d0 + 20.946d0 + 0.934d0 + 0.037d0)
!  fking = (78.084d0 * fk_n2 + 20.946d0 * fk_o2 + 0.97655d0) / 100.001d0
!  depol = 6.d0 * (fking - 1.d0) / (3.d0 + 7.d0 * fking)
!  
!  RETURN
!  
!END SUBROUTINE GET_RAYCOF_DEPOL


SUBROUTINE GET_ALL_RAYCOF_DEPOL(nw, waves, raycof, depol)
  
  USE OMSAO_precision_module
  IMPLICIT NONE

  !	Input/Output
  INTEGER, INTENT(IN)                        :: nw
  REAL (KIND=dp), DIMENSION(nw), INTENT(IN)  :: waves
  REAL (KIND=dp), DIMENSION(nw), INTENT(OUT) :: raycof, depol
  
  !	Local variables
  REAL (KIND=dp), DIMENSION(nw) :: sig, sig2, sig2p, sig4, fk_n2, &
       fk_o2, fking
  REAL (KIND=dp), PARAMETER     :: abod = 1.0455996d0, bbod = -341.29061d0, &
       cbod = -0.90230850d0, dbod = 0.0027059889d0, ebod = -85.968563d0
 
  !	Rayleigh coefficient
  ! Using bodhaine et al, j. atm. oceanic tech. 16, 1854-1861, 1999.
  sig =    1.0d3 / waves
  sig2 =   sig * sig
  sig2p =  1.d0 / sig2
  sig4 =   sig2 * sig2
  raycof = (abod + bbod * sig2 + cbod * sig2p) &
       / (1.d0 + dbod * sig2 + ebod * sig2p) * 1.d-28
  
  !	Derivation of depolarization factor d from king factors for air,
  !	fking = (6 + 3.depol) / (6 - 7.depol) 
  !	bodhaine et al., 370 ppmv co2
  fk_n2 = 1.034d0 + 3.17d-4 * sig2
  fk_o2 = 1.096d0 + 1.385d-3 * sig2 + 1.448d-4 * sig4
  fking = (78.084d0 * fk_n2 + 20.946d0 * fk_o2 + 0.97655d0) / 100.001d0
  depol = 6.d0 * (fking - 1.d0) / (3.d0 + 7.d0 * fking)
  
  RETURN
  
END SUBROUTINE GET_ALL_RAYCOF_DEPOL

SUBROUTINE GET_ALL_RAYCOF(nw, waves, raycof)
  
  USE OMSAO_precision_module
  IMPLICIT NONE
  
  !	Input/Output
  INTEGER, INTENT(IN)                        :: nw
  REAL (KIND=dp), DIMENSION(nw), INTENT(IN)  :: waves
  REAL (KIND=dp), DIMENSION(nw), INTENT(OUT) :: raycof
  
  !	Local variables
  REAL (KIND=dp), DIMENSION(nw) :: sig, sig2, sig2p, sig4
  REAL (KIND=dp), PARAMETER     :: abod = 1.0455996d0, bbod = -341.29061d0, &
       cbod = -0.90230850d0, dbod = 0.0027059889d0, ebod = -85.968563d0
  
  !	Rayleigh coefficient
  ! Using bodhaine et al, j. atm. oceanic tech. 16, 1854-1861, 1999.
  sig =    1.0d3 / waves
  sig2 =   sig * sig
  sig2p =  1.d0 / sig2
  sig4 =   sig2 * sig2
  raycof = (abod + bbod * sig2 + cbod * sig2p) &
       / (1.d0 + dbod * sig2 + ebod * sig2p) * 1.d-28
  
  RETURN
  
END SUBROUTINE GET_ALL_RAYCOF

SUBROUTINE GET_ALL_RAYCOF_DEPOL1(nw, waves, nw1, raycof, depol, problems)

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY : max_spec_pts
  USE OMSAO_variables_module,  ONLY : refdbdir, n_refspec_pts, refspec_orig_data, do_bandavg

  IMPLICIT NONE

  !     Input/Output
  INTEGER, INTENT(IN)                        :: nw, nw1
  REAL (KIND=dp), DIMENSION(nw), INTENT(IN)  :: waves
  REAL (KIND=dp), DIMENSION(nw1), INTENT(OUT):: raycof, depol
  LOGICAL, INTENT(OUT)                       :: problems

  ! Local variables
  !REAL (KIND=dp), DIMENSION(nw1)             ::raycof1, depol1
  REAL (KIND=dp), DIMENSION(max_spec_pts)     :: lowresi0
  REAL (KIND=dp), DIMENSION(nw)               :: raycof1, depol1
  LOGICAL                                     :: get_lowresi0

  INTEGER, SAVE                                 :: nref
  REAL (KIND=dp), DIMENSION(max_spec_pts), SAVE :: ray, dep, refwavs
  REAL (KIND=dp),                          SAVE :: rnorm, dnorm
  LOGICAL, SAVE                                 :: first = .TRUE.

  INTEGER                               :: i, errstat, ni0, ntemp
  REAL (KIND=dp)                        :: scalex
  
  problems = .FALSE.

  IF (first) THEN

     ni0 = n_refspec_pts(1); nref = ni0
     refwavs(1:nref) = refspec_orig_data(1, 1:ni0, 1)

     CALL GET_ALL_RAYCOF_DEPOL(nref, refwavs, ray(1:nref), dep(1:nref))
     rnorm = 1.0E-25; dnorm = 1.0E-2
     ray(1:nref) = ray(1:nref) / rnorm; dep(1:nref) = dep(1:nref) / dnorm
     
     scalex = 1.0  ! dummy variable here
     get_lowresi0 = .FALSE.
     CALL CORRECT_I0EFFECT(refwavs(1:nref), ray(1:nref), nref, &
          refspec_orig_data(1, 1:ni0, 1), refspec_orig_data(1, 1:ni0, 2), ni0, &
          scalex, get_lowresi0, errstat, lowresi0(1:nref))
     IF ( errstat /= 0 ) THEN
        WRITE(*, *) 'Error in correcting I0 effect for raycof!!!'
        problems = .TRUE.; RETURN
     ENDIF
     CALL CORRECT_I0EFFECT(refwavs(1:nref), dep(1:nref), nref, &
          refspec_orig_data(1, 1:ni0, 1), refspec_orig_data(1, 1:ni0, 2), ni0, &
          scalex, get_lowresi0, errstat, lowresi0(1:nref))
     IF ( errstat /= 0 ) THEN
        WRITE(*, *) 'Error in correcting I0 effect for depol!!!'
        problems = .TRUE.; RETURN
     ENDIF

     first = .FALSE.
  ENDIF

  ! interpolation
  CALL BSPLINE(refwavs(1:nref), ray(1:nref), nref, waves, raycof1, nw, errstat)
  IF (errstat < 0) THEN
     WRITE(*, *) 'Error in interpolating raycof!!!'
     problems = .TRUE.; RETURN
  ENDIF
  
  CALL BSPLINE(refwavs(1:nref), dep(1:nref), nref, waves, depol1,  nw, errstat)
  IF (errstat < 0) THEN
     WRITE(*, *) 'Error in interpolating depol!!!'
     problems = .TRUE.; RETURN
  ENDIF

  IF (do_bandavg) THEN
     CALL avg_band_effozcrs(waves, raycof1, nw, ntemp, errstat)
     IF ( errstat /= 0 .OR. ntemp /= nw1) THEN
        WRITE(*, *) 'Raycof Spectra Averaging Error: ', nw, nw1, ntemp
        problems = .TRUE.; RETURN
     ENDIF

     CALL avg_band_effozcrs(waves, depol1, nw, ntemp, errstat)
     IF ( errstat /= 0 .OR. ntemp /= nw1) THEN
        WRITE(*, *) 'Depol Spectra Averaging Error: ', nw, nw1, ntemp
        problems = .TRUE.; RETURN
     ENDIF   
  ENDIF
  
  raycof = raycof1(1:nw1) * rnorm
  depol  = depol1(1:nw1)  * dnorm

  RETURN

END SUBROUTINE GET_ALL_RAYCOF_DEPOL1

