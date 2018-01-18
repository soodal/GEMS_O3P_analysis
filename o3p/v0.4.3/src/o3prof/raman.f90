!! The original model is developed by C.E. Sioris (Sioris and Evans, 2001)

! It is converted to FORTRAN 90 and slightly modified by x.liu
! 1. Converted to FORTRAN 90
! 2. Read input parameters only once
! 3. Optimize some computation
!
! This has been swtiched to a forward model, but R is for elastic and I is for inelastic
! Christopher E. Sioris, Centre for Research in Earth and Space Science
! to appear in Can.J. Phys. 2001.  
! 
! CAVEATS: 1) apply only to obs at moderate spectral res (>=1cm^-1)
! 2) valid only for 1000nm>lambda>200nm
! 3) valid only in weakly absorbing regions or where absorption
! does not exhibit fine structure (i.e. Chappuis). Thus avoid 
! Hartley/Huggins, O2 and H2O bands        (not true anymore)
! 4) assumes no aerosols, no water vapour  (not true anymore)
! 5) assumes all scattering is in line-of-sight (i.e. ignores
! phase, polarization of multiply scattered radiation)
! 6) assumes an effective temperature for the entire atmosphere
! 
! newer versions of Raman scattering model address 1-6
! 
! Using this model
! 
! 1) Change nuhi and nulo to the appropriate upper and lower 
! wavenumber boundaries.
! 2) Change TH to the number of measured spectra
! 3) Save the temperatures(T) and single scattering angles(SSA) 
! as the file tempSSA.txt in the following format:
! T1 (in Kelvin)    SSA1 (in degrees)
! T2                SSA2
! 99.9999 (or any real number to prevent EOF error)
! SEE SAMPLE tempSSA.txt
! 4) Save the measured radiances as follows:   
! spectrum1 (1 data point/line, in ascending order wrt lambda) 
! spectrum2(don't skip any lines between spectra, see sample Int.txt) 

      SUBROUTINE raman(refdir, nulo, nuhi, nline, nz, T, sca, cossza,&
            R, tran,rhos, ring)
         
      USE OMSAO_precision_module
      USE OMSAO_parameters_module, ONLY : maxchlen, pi, o2mix, N2mix, CO2mix
      IMPLICIT NONE

! ========================
! Input/output variables
! ========================
      CHARACTER (LEN=maxchlen), INTENT(IN)                 :: refdir
      INTEGER, INTENT (IN)                                 :: nulo, nuhi, nline, nz
      REAL (KIND=dp), DIMENSION(nulo:nuhi, nz), INTENT(IN) :: R, tran
      REAL (KIND=dp), DIMENSION(nline),        INTENT(OUT) :: ring
      REAL (KIND=dp), DIMENSION(nz),         INTENT(INOUT) :: rhos
      REAL (KIND=dp), INTENT(IN)                           :: sca, cossza, T

! ========================
! Local variables
! ========================
     REAL (KIND=dp), PARAMETER :: c1=1.438769, NL=2.686763D19
     INTEGER,        PARAMETER :: N2Jmax=28, O2maxJ=53, O2max=94, maxpos=218, pixelno=6521 
     
     INTEGER        :: nu, iz, j, k, fidx, lidx
     REAL (KIND=dp) :: ZN2, ZO2, temp, temp1, temp2, phasefnc, pi3, avgt
     REAL (KIND=dp), DIMENSION (nz)                         :: tempz
     REAL (KIND=dp), DIMENSION (nulo:nuhi)                  :: gammaN2, gammaO2
     REAL (KIND=dp), DIMENSION(nulo+maxpos:nuhi-maxpos)     :: RaylP, nr, Raylro, N2so, O2so, I_tot, &
       R_tot, e, Raylcsec
     REAL (KIND=dp), DIMENSION(nulo+maxpos:nuhi-maxpos, nz) :: I, N2sumin, O2sumin, diff
     REAL (KIND=dp), DIMENSION(0:N2Jmax)                    :: N2pop
     REAL (KIND=dp), DIMENSION(0:2*N2Jmax-3)                :: N2csec
     REAL (KIND=dp), DIMENSION(0:O2maxJ)                    :: O2popz
     REAL (KIND=dp), DIMENSION(0:2*O2max-7)                 :: O2csec, O2pop
     CHARACTER (len=maxchlen) :: N2En, N2pos, O2En, O2pos, O2EnfZ, O2J, O2JfZ, O2PT, N2PT 

    LOGICAL,                                 SAVE :: first = .TRUE.
    REAL (KIND=dp), DIMENSION(0:N2Jmax),     SAVE :: N2E
    REAL (KIND=dp), DIMENSION(0:2*N2Jmax-3), SAVE :: N2b
    REAL (KIND=dp), DIMENSION(0:O2maxJ),     SAVE :: O2EnZ
    REAL (KIND=dp), DIMENSION(0:2*O2max-7),  SAVE :: O2E,  O2b
    INTEGER, DIMENSION (0:O2maxJ),           SAVE :: O2JZ
    INTEGER, DIMENSION (0:2*O2max-7),        SAVE :: O2J2, O2shift
    INTEGER, DIMENSION (0:2*N2Jmax-3),       SAVE :: N2shift  
    REAL (kind=dp), Dimension(nz) :: tmpn2, tmpo2, tmpdiff,dum1, dum2 , dum   ! jbak   
    REAL (kind=dp) :: dum_I_toc ,dum_R_toc, dum_ring, dum_one
! Real input parameters needed for each calculation
  IF (first) THEN
     N2En   = TRIM(ADJUSTL(refdir)) // 'raman/N2En.txt'
     N2pos  = TRIM(ADJUSTL(refdir)) // 'raman/N2pos.txt'     
     O2En   = TRIM(ADJUSTL(refdir)) // 'raman/O2En.txt'
     O2pos  = TRIM(ADJUSTL(refdir)) // 'raman/O2pos.txt'
     O2EnfZ = TRIM(ADJUSTL(refdir)) // 'raman/O2EnfZ.txt'
     O2J    = TRIM(ADJUSTL(refdir)) // 'raman/O2J.txt'
     O2JfZ  = TRIM(ADJUSTL(refdir)) // 'raman/O2JfZ.txt'
     O2PT   = TRIM(ADJUSTL(refdir)) // 'raman/O2PT.txt' 
     N2PT   = TRIM(ADJUSTL(refdir)) // 'raman/N2PT.txt'

     ! Read in the line frequencies and strengths
     OPEN (48, file = N2PT,  status='old');  OPEN (49, file = N2pos, status='old')
     DO j = 0, 2 * N2Jmax - 3
        READ (48, *) N2b(j);  READ (49, *) N2shift(j)
     ENDDO
     CLOSE(48); CLOSE(49) 
     
     ! calculate partitioning of N2 
     OPEN(48, file = N2En, status = 'old')     	
     DO j = 0, N2Jmax
        READ (48, *) N2E(j)    
     ENDDO
     CLOSE(48)
     
     ! calculate state sum for O2 
     OPEN(48, file = O2JfZ,  status='old');     OPEN(49, file = O2EnfZ, status='old') 
     do k = 0, O2maxJ
        READ (48,*) O2JZ(k);        READ (49,*) O2EnZ(k)
     ENDDO
     CLOSE(48);  CLOSE(49)
     
     ! O2 cross sections
     OPEN(48, file = O2J,   status = 'old');   OPEN(49, file = O2En,  status = 'old')  
     DO k = 0, 2 * O2max-7
        READ(48, *) O2J2(k); READ(49, *) O2E(k)
     ENDDO
     CLOSE(48);  CLOSE(49)
     
     OPEN(48, file = O2PT,  status = 'old');  OPEN(49, file = O2pos, status = 'old')  
     DO k = 0, 2 * O2max-7
        READ(48, *) O2b(k);  READ(49, *) O2shift(k)
     ENDDO
     CLOSE(48);  CLOSE(49)

     first = .FALSE.
     ENDIF
   
    fidx = nulo + maxpos; lidx = nuhi - maxpos
! calculate dynamic optical parameters
     DO nu = nulo, nuhi
     temp = (1d0 * nu) ** 2 
     gammaN2(nu) = -6.01466E-25 + 2.38557E-14 / (1.86099E10 -temp)
     gammaO2(nu) = 7.149E-26    + 4.59364E-15 / (4.82716E9 - temp)  
     ENDDO
	
     pi3 = pi ** 3
     DO nu = fidx, lidx 
        temp   = (nu / 1d4) ** 2
        nr(nu) = 1d-4 * (0.7041 + 315.9 / (157.39 - temp ) )
        nr(nu) = nr(nu) + 8.4127d-4 / (50.429 - temp )
     
        e(nu) = O2mix * (0.096 + 0.001385 * temp + 1.448d-4 * temp ** 2)
        e(nu) = e(nu) + (N2mix * (0.034 + 0.000317 * temp) + CO2mix * 0.15 )
        e(nu) = e(nu) * 4.5d0 !9.0d0 /2.0d0

! code doesn't work below ~316 nm because nu**2 > maxint, use real         
       Raylcsec(nu) = 32.0d0 * (REAL(nu))**4 * pi3 * (nr(nu))**2 * (1.0d0 + e(nu) / 4.5d0) / 3.d0 
       Raylcsec(nu) = Raylcsec(nu) / NL / NL     
       Raylro(nu)   = 6.0d0 * e(nu) / (45.d0 + 7.d0 * e(nu))         
     ENDDO

! calculate Rayleigh and Raman scattering phase functions
      temp = COS(sca * pi / 180.d0) ** 2
      phasefnc = (13.d0 + temp) * 3.d0 / 40.d0  
      RaylP(fidx:lidx) = 1.d0 + Raylro(fidx:lidx) + (1.d0 - Raylro(fidx:lidx) ) * temp
      RaylP(fidx:lidx) = 3.d0 / (4.d0 + 2.d0 * Raylro(fidx:lidx) ) * RaylP(fidx:lidx)        
  
      temp = 256 * pi ** 5 / 27

! The following loop could be avoided by just using an effective temperature
! Will have small effect on the computed Ring effect spectrum (the effect on
! troospheric column ozone for one orbit is (0.014+/-0.06 DU)
! Use effective temperature, replace T(iz) with T
      
!DO iz = 1,nz	     
!  calculate partitioning of N2  	
!iz = 12
      ZN2=0
      DO j = 0, N2Jmax
         IF ( ifix(10* (j / 2.0 - ifix( j / 2.0) ) ) == 5) THEN
            N2pop(j) = 3 * ( 2 * j + 1.0d0) * EXP(-c1 * N2E(j) / T) 
         ELSE
            N2pop(j) = 6 * ( 2 * j + 1.0d0) * EXP(-c1 * N2E(j) / T)
         ENDIF
         ZN2= ZN2 + N2pop(j)  
      ENDDO
! calculate static part of Anti-Stokes cross sections for N2
       DO j = 0, N2Jmax - 2 
             N2csec(j) = temp * N2pop(j + 2) * N2b(j) /  ZN2
      ENDDO

! calculate static part of Stokes cross sections for N2	
      DO j = N2Jmax - 1, 2 * N2Jmax - 3 
         N2csec(j) = temp * N2pop(j - N2Jmax + 1) * N2b(j) / ZN2
      ENDDO
  
! calculate state sum for O2 
      ZO2 = 0
      DO k = 0, O2maxJ
        O2popz(k) = ( 2 * O2JZ(k) + 1 ) * EXP(- c1 * O2EnZ(k) / T )
        ZO2 = ZO2 + O2popz(k)
      ENDDO
       
! O2 cross sections
       DO k = 0, 2 * O2max-7
        O2pop(k) = (2 * O2J2(k) + 1 ) * EXP(-c1 * O2E(k) / T )
        O2csec(k)= temp * O2pop(k) * O2b(k) / ZO2
      ENDDO


! set arrays to zero initially
       N2so(fidx:lidx) = 0.0; N2sumin(fidx:lidx, 1:nz) = 0.0
       O2so(fidx:lidx) = 0.0; O2sumin(fidx:lidx, 1:nz) = 0.0

! calculate relative amounts of light shifted in/out of a given nu     
     DO nu = fidx, lidx
        temp1 = gammaN2(nu) * gammaN2(nu)
        temp2 = (REAL(nu))**4
        tmpn2(1:nz) = 0.d0  
          DO j = 0, 2 * N2Jmax-3
            N2so(nu) = N2so(nu) + N2csec(j) * (REAL(nu - N2shift(j)))**4 * temp1
            tmpn2(1:nz) =  tmpn2(1:nz) + R(nu + N2shift(j), 1:nz) / R(nu, 1:nz)*N2csec(j) * gammaN2( nu + N2shift(j)) ** 2 * temp2                
          ENDDO
          N2sumin(nu, 1:nz) = tmpn2(1:nz)    
       
       
         temp1 = gammaO2(nu) * gammaO2(nu)
         tmpo2(:) =  0.d0      
         DO k = 0, 2 * O2max-7           
          O2so(nu) = O2so(nu) + O2csec(k) * (REAL(nu - O2shift(k)))**4 * temp1
          tmpo2(1:nz) =  tmpo2(1:nz) + R(nu + O2shift(k), 1:nz) / R(nu, 1:nz) * O2csec(k) &
             * gammaO2(nu + O2shift(k)) ** 2 * temp2
        
        ENDDO
        O2sumin(nu, 1:nz) = tmpo2(1:nz)

        tmpdiff(1:nz) = N2mix * (tmpN2(1:nz) - N2so(nu)) + O2mix * (tmpO2(1:nz) - O2so(nu))
        diff(nu, 1:nz) = tmpdiff(1:nz) 
        I(nu, 1:nz) = R(nu, 1:nz)  * ( 1.0 + phasefnc * tmpdiff(1:nz) / (Raylcsec(nu) * RaylP(nu)))
     
      
         tempz(1:nz) =  rhos(1:nz) *  tran(nu, 1:nz)   
         dum1(1:nz) = R(nu, 1:nz)  * ( 1.0 + phasefnc * tmpdiff(1:nz) / (Raylcsec(nu) * RaylP(nu))) ! I(nu,1:nz)
    !     dum1(1:nz)  =  I(nu, 1:nz)
         dum2(1:nz) = R(nu, 1:nz)
         I_tot(nu)  = sum(dum1(1:nz)*tempz(1:nz)) !  I_tot(nu) = sum(I*tempz)
         R_tot(nu)  = sum(dum2(1:nz)*tempz(1:nz)) ! I_tot(nu) = sum(r*tempz)
         ring(nu - nulo + 1) = I_tot(nu)/R_tot(nu)  - 1.0d0
    ENDDO
 
  
       ring(1:maxpos) = 0.d0
       ring(nline - maxpos + 1 : nline) = 0.d0

     

      RETURN
  
      END SUBROUTINE raman
