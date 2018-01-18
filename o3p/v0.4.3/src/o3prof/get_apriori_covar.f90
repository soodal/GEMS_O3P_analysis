! ==============================================================
! Construct a priori covariance for ozone based on
! ozone standard deviation of Fortuin Climatology 
! Diagonal elements are directly from this climatology
! Non-diagonal elements are calculated by assuming a 
! correlation length (5 km for now)
! ==============================================================
SUBROUTINE get_apriori_covar(toz, ozprof, sao3)
  
  USE OMSAO_precision_module 
  USE OMSAO_parameters_module, ONLY: p0
  USE ozprof_data_module,     ONLY: atmos_unit, which_aperr,  &
       min_serr, min_terr, loose_aperr, use_logstate, atmosprof, nlay, ntp
  USE OMSAO_variables_module, ONLY: atmdbdir, the_year, the_month, the_day, the_lon, the_lat
  IMPLICIT NONE
  
  ! ======================
  ! Input/Output variables
  ! ======================
  REAL (KIND=dp),INTENT(IN) :: toz
  REAL (KIND=dp), DIMENSION(nlay, nlay), INTENT(OUT) :: sao3
  REAL (KIND=dp), DIMENSION(nlay),      INTENT(IN) :: ozprof
  
  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER :: mref=60
  INTEGER            :: errstat
  REAL (KIND=dp), PARAMETER       :: corrlen=6.0     ! changed from 6 km to 8 km (more uniform in a priori influence)
  
  INTEGER :: nz ! nz = nlay
  REAL (KIND=dp), DIMENSION(0:nlay)       :: ps, zs, ts
  REAL (KIND=dp), DIMENSION(nlay)       :: zmid
  REAL (KIND=dp), DIMENSION(0:nlay)     :: pslg, nstd, nstd1, ps1, zs1
  REAL (KIND=dp) :: tmp
  REAL (KIND=dp), DIMENSION(mref)       :: astd, a1, a2, a3 
  REAL (KIND=dp), DIMENSION(0: mref)    :: cumastd, preslg, pres

  INTEGER                               :: i, j, k,mnorstd, tmpntp, nref


  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=17), PARAMETER :: modulename = 'get_apriori_covar'  

  nz = nlay
  ps(0:nz) = atmosprof(1, 0:nz)
  zs(0:nz) = atmosprof(2, 0:nz)
  ts(0:nz) = atmosprof(3, 0:nz)

  ! ==============================
  ! get astd
  ! ==============================   
  sao3 = 0.0; astd = 0.0

  nref = 60
  pres(1:60) = (/(1013.25*10.0**(-1.0*i/16.0), i = 59, 0, -1)/)
  pres(0) = 0.05  ! about 70 km

 IF ( which_aperr == 1 )  call get_mpstd(astd(1:nref))
 IF ( which_aperr == 12 ) call get_mlprof(astd(1:nref), 2)
 IF ( which_aperr == 13 ) call get_tjprof(astd(1:nref), 2)  

 IF ( which_aperr >= 2 .and. which_aperr <= 4 )  call get_tbstd (astd(1:nref))

 IF ( which_aperr == 6) then  ! Fourtuine 
     nref = 19 ! need more '
     pres(0:nref) = (/0.05, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 20., &
             30., 50., 70., 100., 150., 200., 300., 500., 700., 1000.0/)  
     call get_fortstd (astd(1:nref ))
 ENDIF 

  ! ==============================
  ! get nstd
  ! ==============================   
   
  ! Bondary layer correction
    IF (ps(nz) > p0) then !sfc
        tmp = ( ps(nz) - p0)/(pres(nref)-pres(nref-1))
        astd(nref) = astd(nref)*(1+tmp)
        pres(nref) = ps(nz)
    ENDIF
    IF ( ps(0) < pres(0) ) pres(0) = ps(0) !top

  ! convert  partial column to accumulate
    cumastd(0) = 0.0
    DO i = 1, nref
       cumastd(i) = cumastd(i-1) + astd(i) 
    ENDDO   

    preslg = LOG(pres); pslg = LOG(ps)
    CALL BSPLINE(preslg(0:nref), cumastd(0:nref),nref+1, pslg(0:nz),&
         nstd(0:nz), nz+1, errstat)
    IF (errstat < 0) THEN
       WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat; STOP
    ENDIF

  ! Contruct the full covariance matrix for ozone (in Dobson units)
  nstd(1:nz) = nstd(1:nz) - nstd(0:nz-1)
  !print *, SUM(nstd(ntp+1:nz)) / SUM(ozprof(ntp+1:nz))
  !nstd(1:nz) =  ozprof(1:nz) * 0.5
 
  IF (which_aperr == 9) THEN
     ps1(0) = ps(nz)
     DO i = 1, nz
        ps1(i) = ps(nz-i); nstd1(i) = nstd(nz-i+1)
     ENDDO
     CALL GET_GEOSCHEM_O3STD(the_month, the_lon, the_lat, ps1, nstd1(1:nz), nz, nz-ntp)  
     DO i = 1, nz
        nstd(i) = nstd1(nz-i+1)
     ENDDO
  ELSE IF (which_aperr == 10) THEN
     ps1(0) = ps(nz)
     DO i = 1, nz
        ps1(i) = ps(nz-i); nstd1(i) = nstd(nz-i+1)
     ENDDO
  
     mnorstd = 2
     CALL get_mlso3prof(the_year, the_month, the_day, the_lat, nz, mnorstd, ps1(0:nz), zs1(0:nz), nstd1(1:nz), tmpntp, errstat)
     IF (errstat < 0) THEN
        WRITE(*, *) modulename, ': Error in getting MLS ozone variabilities!!!'; STOP
     ENDIF
     DO i = 1, nz
        nstd(i) = nstd1(nz-i+1)
     ENDDO
  ELSE IF (which_aperr == 11) THEN
     ps1(0) = ps(nz)
     DO i = 1, nz
        ps1(i) = ps(nz-i); nstd1(i) = nstd(nz-i+1)
     ENDDO
  
     mnorstd = 2
     CALL get_mlso3prof_single(the_year, the_month, the_day, the_lat, nz, mnorstd, ps1(0:nz), zs1(0:nz), nstd1(1:nz), tmpntp, errstat)
     IF (errstat < 0) THEN
        WRITE(*, *) modulename, ': Error in getting MLS ozone variabilities!!!'; STOP
     ENDIF
     DO i = 1, nz
        nstd(i) = nstd1(nz-i+1)
     ENDDO
  ENDIF

  ! Loose a priori constraint (because those from climatology are sometimes too small)
  IF (loose_aperr) THEN !jbak
     DO i = 1, ntp-1 
!        IF (nstd(i) / ozprof(i) < min_serr) THEN
           nstd(i) = ozprof(i) * min_serr
!        ENDIF
     ENDDO
     
     DO i = ntp, nz 
!        IF (nstd(i) / ozprof(i) < min_terr) THEN
           nstd(i) = ozprof(i) * min_terr
!        ENDIF
     ENDDO     
  ENDIF
  
  IF (use_logstate) nstd(1:nz) = nstd(1:nz)/ozprof(1:nz)

  DO i = 1, nz
     sao3(i, i)= nstd(i) ** 2.0 
  ENDDO
  
  ! This is based on retrieval stastistics 
  zmid = (zs(0:nz-1) + zs(1:nz)) / 2.0  
  DO i = 1, nz
     DO j = 1, i - 1
      
        sao3(i, j) = SQRT(sao3(i,i) * sao3(j, j)) * &
             EXP(- ABS((zmid(i)-zmid(j)) / corrlen)**2 )
        sao3(j, i) = sao3(i, j) 
     ENDDO
  ENDDO
 
    
  RETURN
END SUBROUTINE get_apriori_covar





! ===============================================================
! Obtain McPeter std profiles (12 month, 18 latitude bands, 61 levels)
! ===============================================================
SUBROUTINE get_mpstd(std)
  
  USE OMSAO_precision_module 
  USE ozprof_data_module,     ONLY: atmos_unit, which_aperr
  USE OMSAO_variables_module, ONLY: atmdbdir, the_month, the_day, the_lat
  IMPLICIT NONE

  INTEGER, PARAMETER                             ::  nref = 60
  ! ======================
  ! Input/Output variables
  ! ======================

  REAL (KIND=dp), DIMENSION(nref), INTENT(OUT)   :: std
  
  ! ======================
  ! Local variables
  ! ======================
  CHARACTER (LEN=130)           :: apfname

  INTEGER, PARAMETER :: nlat=18, nmon=12, nlay=61
  REAL (KIND=dp), DIMENSION(nlay) :: std0, pres
  INTEGER, DIMENSION(2)           :: latin, monin
  REAL (KIND=dp)                  :: frac
  REAL (KIND=dp), DIMENSION(2)    :: latfrac, monfrac
  INTEGER                         :: i, j, k, nband, nm, idum

  ! ======================
  ! SAVE variables
  ! ======================
  LOGICAL, SAVE                 :: first = .TRUE.
  REAL (KIND=dp), SAVE, DIMENSION(nmon, nlat, nlay) :: stds 

  
! ** load std profiles ** !
!  IF (which_aperr == 1) print *,'Sao3 set to be LLM'
  IF (first) THEN
        apfname = TRIM(ADJUSTL(atmdbdir)) // 'mpclima/llmclima_std.dat'
        OPEN (UNIT = atmos_unit, file=apfname, status = 'unknown')
        READ (atmos_unit, '(A)') ;  READ(atmos_unit, '(A)') 
        DO i = 1, nmon   
           READ(atmos_unit, '(A)') ;  READ(atmos_unit, '(A)')  ! read month label
           DO k = nlay, 1, -1
              READ(atmos_unit, *) idum, (stds(i, j, k), j=1, nlat) ! ppmv
           ENDDO
        ENDDO
        CLOSE(atmos_unit)
        first = .FALSE.
 ENDIF
         
! ** interpolation for lat, mon** ! 
  IF (the_day <= 15) THEN
     monin(1) = the_month - 1
     IF (monin(1) == 0) monin(1) = 12
     monin(2) = the_month
     monfrac(1) = (15.0 - the_day) / 30.0
     monfrac(2) = 1.0 - monfrac(1)
  ELSE 
     monin(2) = the_month + 1
     IF (monin(2) == 13) monin(2) = 1
     monin(1) = the_month
     monfrac(2) = (the_day - 15) / 30.0
     monfrac(1) = 1.0 - monfrac(2)
  ENDIF
  nm = 2


  IF (the_lat <= -85.0) THEN
        nband = 1; latin(1) = 1; latfrac(1) = 1.0
  ELSE IF (the_lat >= 85.0) THEN
        nband = 1; latin(1) = nlat; latfrac(1) = 1.0
  ELSE
        nband = 2     ; frac = (the_lat + 85.0) / 10.0 + 1
        latin(1) = INT(frac); latin(2) = latin(1) + 1
        latfrac(1) = latin(2) - frac; latfrac(2) = 1.0 - latfrac(1)
  ENDIF

  std0 =0.0
  DO i = 1, nband
        DO j = 1, nm
           std0 =  std0 + stds(monin(j), latin(i), :) * monfrac(j) * latfrac(i)
        ENDDO
  ENDDO
     
! ** convert ppm into DU ** ! 
  pres(1) = 0.05  ! about 70 km
  pres(2:61) = (/(1013.25*10.0**(-1.0*i/16.0), i = 59, 0, -1)/)
  DO i = 1, nref
     std(i) = (std0(i+1) + std0(i))*0.5 *(pres(i+1) - pres(i))/ 1.267
  ENDDO

  RETURN

END SUBROUTINE get_mpstd

! ===============================================================
! Obtain TB hybrid variance
! drive variance by varialbe shifht or non shifht depending on
!      the number of sources with trpz below 14 km and latitude
! And then merging with LLM variance
! ** variable shift region    between tropz - 5 km and tropz + 5 km
! ** non shift region         below tropz tropz-5 km
! ** merging with LLM region  between tropz+5km and tropz+10 km
! ** LLM                      above tropz + 10 km 
! 2011.6.15 Jbak
! ===============================================================
SUBROUTINE get_tbstd(std)

  USE OMSAO_precision_module 
  USE ozprof_data_module,     ONLY: atmos_unit,trpz, which_aperr
  USE OMSAO_variables_module, ONLY: atmdbdir, the_month, the_day, the_lat
  IMPLICIT NONE


  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, PARAMETER                             :: nref = 60
  REAL (KIND=dp), DIMENSION(nref), INTENT(OUT)   :: std
  
  ! ======================
  ! Local variables
  ! ======================
  INTEGER            :: i , which_offset
  REAL (KIND=dp), DIMENSION(nref):: std1, std2,llm,ml,ab,refz,tmp
  REAL (KIND=dp) :: meg,meg1, meg2, weight,weight1, weight2, trpz1, trpz2, del, del1, del2
 trpz1 = 13
 trpz2 = 15
 meg1  = 5
 meg2  = 5

 del1 = meg1 + 5
 del2 = meg2 + 5

 if (trpz < 10 ) then 
   meg2 = 3; meg2=3
 endif

 refz(1:nref) = (/(i*1.0, i= 59, 0,-1 )/)

! ** basic clima

 CALL get_mpstd(llm)
!CALL get_mlprof(ml,2) ;  llm(:) =ml(:)
 
IF (which_aperr == 4) then
 call get_ab(tmp,std, 1)
ELSE IF (which_aperr == 3) then 
 call get_tb (tmp, std1,1)  
 call get_tb (tmp, std2,2) 
 call get_ab (tmp, ab,1)   
!      std2(:) = AB(:)
ENDIF


! Vertical mixing.
DO i = 1, nref


   weight1 = 1-(abs( refz(i) - trpz )-meg1)/(del1-meg1)
   weight2 = 1-(abs( refz(i) - trpz )-meg2)/(del2-meg2)
  if ( weight1 < 0 )  weight1 = 0
  if ( weight1 > 1 )  weight1 = 1
  if ( weight2 < 0 )  weight2 = 0
  if ( weight2 > 1 )  weight2 = 1
  
  IF ( refz(i) >= trpz ) then ! above the tropopause
      IF ( which_aperr == 4 ) then 
       std(i) = std(i)*weight1 +LLM(i)*(1-weight1)
      elseif ( which_aperr == 3 ) then 
       std1(i) = std1(i)*weight1 +LLM(i)*(1-weight1)
       std2(i) = std2(i)*weight1 +LLM(i)*(1-weight1)
      ENDIF
  endif
  IF ( which_aperr == 3 .and. refz(i) < trpz) then ! below the tropopause 
       std1(i) = std1(i)*weight2 +ab(i)*(1-weight2)
       std2(i) = std2(i)*weight2 +ab(i)*(1-weight2)
  endif           
ENDDO 
! spaceing mixing 
IF ( which_aperr == 3 ) then  
 weight = (trpz2-trpz)/(trpz2-trpz1)
 if (trpz <= trpz1) then  !100 % TB
          weight = 1
 else if (trpz > trpz2) then !100 % AB  	 
          weight = 0
 endif 
 STD(:) = std1(:)*weight + std2(:)*(1-weight)      
ENDIF


IF (any(std(:) < 0)) then ; print * , 'error at get_tbprof' ; stop ; ENDIF

RETURN

END SUBROUTINE get_tbstd



! ===============================================================
! Obtain IUP std profiles (2 season, 6 latitude bands,
!   14 profiles with total ozone at a step of 50 DU
! ===============================================================
SUBROUTINE get_iupstd(toz, std)

  USE OMSAO_precision_module 
  USE OMSAO_variables_module, ONLY: atmdbdir ,the_month, the_lat
  USE ozprof_data_module,     ONLY: atmos_unit ,pst,which_aperr
  IMPLICIT NONE


  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, PARAMETER                           ::  nref = 60
  REAL (KIND=dp), INTENT(IN)                   :: toz
  REAL (KIND=dp), DIMENSION(nref), INTENT(OUT) :: std

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER :: nlat=6, maxprof=14, nmon=2
  REAL (KIND=dp), DIMENSION(nref)   :: iupstd, llmstd, weight,refz
  REAL (KIND=dp), DIMENSION(14)                    :: temp1, temp2
  REAL (KIND=dp)                                   :: frac, fdum, maxoz, minoz,meg
  REAL (KIND=dp), DIMENSION(2)                     :: latfrac
  INTEGER,        DIMENSION(2)                     :: latin
  INTEGER                                          :: season

  CHARACTER (LEN=130)                              :: apfname
  INTEGER :: i, im, ib, il, k,iprof, profin, nprof, nband
  character(10), dimension(6) :: bandname =['90S-60S','60S-30S','30S-0S', '0N-30N', '30N-60N', '60N-90N'] 
  ! ======================
  ! saved variables
  INTEGER,        SAVE, DIMENSION(nmon, nlat)              :: nprofs
  INTEGER,        SAVE, DIMENSION(nmon, nlat, maxprof)     :: tozindex
  REAL (KIND=dp), SAVE, DIMENSION(nmon, nlat, maxprof,nref):: stds 
  LOGICAL,        SAVE                                     :: first = .TRUE.
  ! ======================

! ** load std profiles ** !
IF (which_aperr == 5) print *,'Sao3 set to be IUP with toz',toz  

IF (first) THEN
    apfname = TRIM(ADJUSTL(atmdbdir)) // 'iupclima/iupclima_o3du_sd.dat'
    OPEN (UNIT = atmos_unit, file=  apfname, status = 'unknown')
        
    ! read loop        
    DO i = 1, 6 ; READ (atmos_unit, '(A)') ;ENDDO

    DO im = 1, nmon
        Do ib = 1,  nlat
           READ(atmos_unit, '(A)')   ! read month label         
           READ(atmos_unit,*) (temp1(iprof), iprof=1, maxprof)! read total ozone label    

           DO il = 1, nref ! read bottom -> top
           READ(atmos_unit, *) fdum, (temp2(iprof),iprof = 1, maxprof ) ! du  
   
            nprof = 1         
            DO k = 1, maxprof             
                   IF (temp1(k) /= 0.0 ) then                   
                   stds(im,ib,nprof, il) = temp2(k) 
                   tozindex(im, ib, nprof) = temp1(k)                       
                   nprof = nprof + 1
                   ENDIF
            ENDDO
                   nprofs (im,ib) = nprof-1      
           ENDDO
        ENDDO
     ENDDO       

    CLOSE (atmos_unit)  
   first = .FALSE.
ENDIF 

! ** interpolation for lat, mon** ! 

  !find latitude fraction  
     IF (the_lat <= -75.0) THEN
        nband = 1; latin(1) = 1; latfrac(1) = 1.0
     ELSE IF (the_lat >= 75.0) THEN
        nband = 1; latin(1) = nlat; latfrac(1) = 1.0
     ELSE
        nband = 2     ; frac = (the_lat + 75.0) / 30.0 + 1
        latin(1) = INT(frac); latin(2) = latin(1) + 1
        latfrac(1) = latin(2) - frac; latfrac(2) = 1.0 - latfrac(1)
     ENDIF
  !find total ozone fraction
     
   season = 1  
   IF (the_month > 5 .and. the_month < 12 ) season = 2 
      

  !find total ozone fraction
      iupstd =0.

     DO ib = 1, nband

           nprof = nprofs(season, latin(ib))
           minoz = tozindex(season, latin(ib), 1)
           maxoz = tozindex(season, latin(ib), nprof) 

                      
           IF (toz < minoz) THEN
              WRITE(*,*), 'Warning: no a priori profile available!!!'
               iupstd  = iupstd + stds(season, latin(ib), 1, :)  * latfrac(ib)
           ELSE IF (toz > maxoz) THEN
              WRITE(*,*), 'Warning: no a priori profile available!!!'
               iupstd = iupstd + stds(season, latin(ib), nprof, :)  * latfrac(ib)
           ELSE

               profin = INT ((toz - minoz ) / 30.0)+1

              IF (profin == 0) THEN 
                 profin = 1
              ELSE IF (profin == nprof) THEN
                 profin = profin -1
              ENDIF
              
              frac = 1.0 - (toz - (minoz + (profin-1) * 30.0)) / 30.0

              iupstd = iupstd + latfrac(ib) * (frac * stds(season, latin(ib), profin, :) &
                   + (1.0 - frac) * stds(season, latin(ib), profin+1, :))
           ENDIF
     
     ENDDO
               
  ! merging with LLM Clim
  refz(1:nref) = (/(i*1.0+0.5, i= 59, 0,-1 )/)

  meg =25
  CALL get_mpstd(llmstd)
  DO i = 1, nref 
     IF ( refz(i) <= meg ) THEN 
      weight(i) = 1.0
     ELSE IF ( refz(i) > meg+5 ) THEN 
      weight(i) = 0.0
     ELSE 
      weight(i) = (meg+5-refz(i))/5.
     ENDIF 
      std(i) = iupstd(i)*weight(i) + llmstd(i)*(1-weight(i)) 
      ! print * , refz(i), weight(i), std(i)
  ENDDO

  RETURN
END SUBROUTINE get_iupstd




! ===============================================================
! Obtain fortstd profiles (12 month, 17 latitude bands, 19 levels)
! ===============================================================
SUBROUTINE get_fortstd( std)
  
  USE OMSAO_precision_module 
  USE ozprof_data_module,     ONLY: atmos_unit
  USE OMSAO_variables_module, ONLY: atmdbdir, the_month, the_day, the_lat
  IMPLICIT NONE

  INTEGER, PARAMETER                            ::  nref = 19
  ! ======================
  ! Input/Output variables
  ! ======================
  REAL (KIND=dp), DIMENSION(nref), INTENT(OUT)  :: std
  
  ! ======================
  ! Local variables
  ! ======================
  CHARACTER (LEN=130)           :: apfname
  INTEGER, PARAMETER :: nlat=17, nm=12, nlay=20
  
  REAL (KIND=dp), DIMENSION(nlay) :: std0, pres
  INTEGER, DIMENSION(2)           :: latin, monin
  REAL (KIND=dp)                  :: frac
  REAL (KIND=dp), DIMENSION(2)    :: latfrac, monfrac
  INTEGER                         :: i, j, k, nband, nmon

  LOGICAL, SAVE                   :: first = .TRUE.
  REAL (KIND=dp), SAVE, DIMENSION(nm, nlat, nlay) :: stds 

  
! ** load std profiles ** !
  IF (first) THEN
  print * , 'fk error'
        apfname = TRIM(ADJUSTL(atmdbdir)) // 'fkclima/fortuin_o3_sdev.dat'
        
        OPEN (UNIT = atmos_unit, file=apfname, status = 'unknown')
        DO i = 1, nm 
           READ(atmos_unit, '(A)')  ! read month label
           READ(atmos_unit, *) ((stds(i, j, k), j=1, nlat), k=nlay, 1, -1) ! ppmv
        ENDDO
        CLOSE(atmos_unit)
        first = .FALSE.
 ENDIF

! ** interpolation for lat, mon** ! 
  IF (the_day <= 15) THEN
     monin(1) = the_month - 1
     IF (monin(1) == 0) monin(1) = 12
     monin(2) = the_month
     monfrac(1) = (15.0 - the_day) / 30.0
     monfrac(2) = 1.0 - monfrac(1)
  ELSE 
     monin(2) = the_month + 1
     IF (monin(2) == 13) monin(2) = 1
     monin(1) = the_month
     monfrac(2) = (the_day - 15) / 30.0
     monfrac(1) = 1.0 - monfrac(2)
  ENDIF
  nmon = 2


  IF (the_lat <= -80.0) THEN
     nband = 1; latin(1) = 1; latfrac(1) = 1.0
  ELSE IF (the_lat >= 80.0) THEN
     nband = 1; latin(1) = nlat; latfrac(1) = 1.0
  ELSE
     nband = 2     ; frac = (the_lat + 80.0) / 10.0 + 1
     latin(1) = INT(frac); latin(2) = latin(1) + 1
     latfrac(1) = latin(2) - frac; latfrac(2) = 1.0 - latfrac(1)
  ENDIF

  std0 =0.0
  DO i = 1, nband
     DO j = 1, nmon
           std0 =  std0 + stds(monin(j), latin(i), :) * monfrac(j) * latfrac(i)
     ENDDO
  ENDDO


  pres(1:nlay) = (/0.05, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 20., &
             30., 50., 70., 100., 150., 200., 300., 500., 700., 1000.0/) 
  DO i = 1, nref
     std(i) = (std0(i+1) + std0(i))*0.5 *(pres(i+1) - pres(i))/ 1.267
  ENDDO

  RETURN

END SUBROUTINE get_fortstd 

! The following routine is outdated
! Two step approach does not work

!SUBROUTINE get_first_retrieval(month, lat, nz, ps, zs, ts, ozprof, sao3)
!
!  USE OMSAO_precision_module 
!  USE OMSAO_variables_module, ONLY : ch1_out_file 
!  IMPLICIT NONE
!  

!  ! ======================
!  ! Input/Output variables
!  ! ======================
!  INTEGER, INTENT(IN)                            :: nz, month
!  REAL (KIND=dp), INTENT(IN)                     :: lat
!  REAL(KIND=dp),  DIMENSION(nz), INTENT(INOUT)   :: ozprof
!  REAL (KIND=dp), DIMENSION(nz, nz), INTENT(OUT) :: sao3
!  REAL (KIND=dp), DIMENSION(0:nz), INTENT(INOUT) :: zs, ps, ts
!
!  INTEGER                           :: i, j, idum
!  REAL (KIND=dp)                    :: fdum
!  REAL (KIND=dp), DIMENSION(nz)     :: zmid, ozs
!  REAL (KIND=dp), PARAMETER         :: corrlen=6.0
!  CHARACTER (LEN=18), PARAMETER     :: modulename = 'get_first_retrieval' 
!  CHARACTER (LEN=10)                :: atmos_str
!  CHARACTER (LEN=10)                :: covar_str
!  CHARACTER (LEN=130)               :: fname
!
!  fname = ch1_out_file
!  fname = '/home/xliu/OzoneFit/OZBOREAS-IMPROVE/testout/lv2_70409103_0avg_b1a.out'
!  OPEN(UNIT=55, FILE=fname, STATUS = 'OLD')
!  DO 
!     READ(55, '(A10)') atmos_str
!     IF (atmos_str == 'Atmosphere') EXIT
!  ENDDO
!  READ(55, '(A)')
!  READ(55, '(A)')
!
!  DO i = 1, nz
!     READ(55, *) idum, fdum, fdum, fdum, fdum, fdum, ozs(i)
!  ENDDO
!  ozprof(1:nz) = ozs(1:nz)
!
!  DO 
!     READ(55, '(A10)') covar_str
!     IF (covar_str == 'Covariance') EXIT
!  ENDDO
!  READ(55, *) sao3
!  CLOSE (UNIT=55)
!
!  zmid = (zs(0:nz-1) + zs(1:nz)) / 2.0  
!  DO i = 1, nz
!     DO j = 1, i - 1
!        sao3(i, j) = SQRT(sao3(i,i) * sao3(j, j)) * &
!             EXP(- ABS((zmid(i)-zmid(j)) / corrlen)**2.0 )
!        sao3(j, i) = sao3(i, j) 
!     ENDDO
!  ENDDO
!  
!  RETURN
!
!END SUBROUTINE GET_FIRST_RETRIEVAL

