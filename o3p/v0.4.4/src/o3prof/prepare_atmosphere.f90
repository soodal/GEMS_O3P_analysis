  ! ************************************************************************
  ! Author:  xiong liu
  ! Date  :  July 24, 2003
  ! Purpose: Routine to read atmospheric surface and tropopause pressure, 
  !          temperature, ozone, trace gases, surface altitdue and so on.
  ! adding : ML (2011), TB (2012), TJ(2013)  climaotlogy by jbak
  ! ************************************************************************


! =====================================================================
! Obtain TB hybrid oz profiles
! drive ozone profile by variable shift or non shift depending on
!      the number of profile with trpz below 14 km and latitude
! And then merging with LLM profile 
! ** variable shift region    between tropz - 5 km and tropz + 5 km
! ** non shift region         below tropz tropz-5 km
! ** merging with LLM region  between tropz+5km and tropz+10 km
! ** LLM                      above tropz + 10 km 
! 2011.6.15 Jbak
! ======================================================================
SUBROUTINE get_tbprof (ozref)

  USE OMSAO_precision_module 
  USE ozprof_data_module,     ONLY: atmos_unit, trpz,which_clima
  USE OMSAO_variables_module, ONLY: atmdbdir, the_month, the_day, the_lat
  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, PARAMETER                             ::  nref = 60
  REAL (KIND=dp), DIMENSION(nref), INTENT(OUT)   ::  ozref

  ! ======================
  ! Local variables
  ! ======================
  REAL (KIND=dp), DIMENSION(nref):: ozref1, ozref2, llm,ml,ab, tmp, refz
  REAL (KIND=dp)                 :: weight, weight1, weight2,meg, trpz1, trpz2, del, del1, del2, meg1, meg2
  INTEGER :: i,which_offset
  trpz1 = 13 ; trpz2 = 15
  meg1   = 5 ; meg2   = 5
 
  del1 = 5 + meg1
  del2 = 5 + meg2
  if ( trpz < 10 ) then 
    meg2 = 3 ; del2 = 3 + meg2
  endif
  refz(1:nref) = (/(i*1.0, i= 0, 59 )/)

! Load basic clima 

  CALL get_mcprof(llm)
 
  IF (which_clima == 4) then 
    call get_AB (ozref,tmp, 1) 
  ELSE IF (which_clima == 3 ) THEN   
    call get_tb (ozref1, tmp,1)    
    call get_tb (ozref2, tmp,2)
    call get_ab (AB, tmp,1) 
  ENDIF

! Vertical mixing. 
  DO i = 1, nref
    weight1 = 1-(abs( refz(i) - trpz )-meg1)/(del1-meg1)
    weight2 = 1-(abs( refz(i) - trpz )-meg2)/(del2-meg2)
    if ( weight1 < 0 )  weight1 = 0
    if ( weight1 > 1 )  weight1 = 1
    if ( weight2 < 0 )  weight2 = 0
    if ( weight2 > 1 )  weight2 = 1
 
  
    IF ( refz(i) >= trpz ) then 
      IF ( which_clima == 4 ) then 
        ozref(i) = ozref(i)*weight1 +LLM(i)*(1-weight1)
      elseif ( which_clima == 3 ) then 
        ozref1(i) = ozref1(i)*weight1 +LLM(i)*(1-weight1)
        ozref2(i) = ozref2(i)*weight1 +LLM(i)*(1-weight1)
      ENDIF
    endif
    if ( which_clima == 3 .and. refz(i) < trpz ) then 
      ozref1(i) = ozref1(i)*weight2 +AB(i)*(1-weight2)
      ozref2(i) = ozref2(i)*weight2 +AB(i)*(1-weight2)   
    ENDIF
  ENDDO 

  IF ( which_clima == 3 ) THEN 
    weight = (trpz2-trpz)/(trpz2-trpz1)
    IF (trpz <= trpz1) THEN  !100 % TB
      weight =1
    ELSE IF (trpz > trpz2) THEN !100 % AB
      weight =0
    ENDIF   
    ozref(:) = ozref1(:)*weight + ozref2(:)*(1-weight)
    write(*,'(2(a7, i6),6(a7,f6.2), f6.2)') 'MON:', the_month,' DAY:',the_day,' LAT:', the_lat, ' TPH', trpz,'=TBx',weight,' del: ',del1,' tbllm:',meg,' smooth:',trpz1,trpz2
  ENDIF

  IF (any(ozref(:) < 0)) then ; print * , 'error at get_tbprof' ; stop ; ENDIF

  RETURN

END SUBROUTINE get_tbprof


! ===============================================================
! Obtain TB-based oz profiles (12 month, 18 latitude bands, 80 layers : ppb)
! Variable shifht
! 2011.6.2 Jbak
! ===============================================================

SUBROUTINE get_tb(ozref,std,which_tb)

  USE OMSAO_precision_module 
  USE ozprof_data_module,     ONLY: atmos_unit, trpz, mzt
  USE OMSAO_variables_module, ONLY: atmdbdir, the_month, the_day, the_lat
  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, PARAMETER                             ::  nref = 60
  INTEGER, INTENT(IN)                            ::  which_tb
  REAL (KIND=dp), DIMENSION(nref), INTENT(OUT)   ::  ozref,std
  
  ! ======================
  ! Local variables
  ! ======================
  CHARACTER (LEN=130)           :: apfname

  INTEGER, PARAMETER :: nlat=18, nmon=12, nlay=80
  REAL (KIND=dp), PARAMETER         :: lat0=-90., latgrid=10.
  REAL (KIND=dp), DIMENSION(nlay)   :: ozref0,std0 ! orignal profile
  REAL (KIND=dp), DIMENSION(0:nlay) :: cum0,cums0, refz0, zstar, tb0
  REAL (KIND=dp), DIMENSION(0:nref) :: cum,cums,refz, offset, tb

  INTEGER                           :: i, j, k,fidx, lidx, nband, nm, errstat
  INTEGER, DIMENSION(2)             :: latin, monin
  REAL (KIND=dp)                    :: frac,fdum
  REAL (KIND=dp), DIMENSION(2)      :: latfrac, monfrac
  REAL (KIND=dp)                    :: meg
  REAL (KIND=dp)                    :: gravity_correct ! used for converting unit


  LOGICAL, SAVE                     :: first = .TRUE.
  REAL (KIND=dp), SAVE, DIMENSION(nmon, nlat, nlay) ::ozrefs,ozrefs1, ozrefs2
  REAL (KIND=dp), SAVE, DIMENSION(nmon, nlat, nlay) ::stds,stds1, stds2
  REAL (KIND=dp), SAVE, DIMENSION(nmon, nlat)       ::mtropz, mtropz1, mtropz2
  REAL (KIND=dp), SAVE, DIMENSION(nlat)             ::lats
  REAL (KIND=dp), SAVE, DIMENSION(nlay) :: z0

  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=17), PARAMETER :: modulename = 'get_TB_VS_OZ' 

! ** load oz profiles ** !
  IF (first) THEN
          
        !apfname = '../../ATMOS/tbclima/TB14L-5.vs' ! jbak
        apfname = TRIM(ADJUSTL(atmdbdir)) //'tbclima/TB14L-5.vs' ! geun
        OPEN (UNIT = atmos_unit, file=apfname, status = 'unknown')
        READ (atmos_unit, '(A)') ;  READ(atmos_unit, '(A)')             
        DO i = 1, nmon
          DO j = nlat, 1, -1
           READ(atmos_unit, *) fdum, lats(j), mtropz1(i,j) ! nsample, lat, mean ztrop
           READ(atmos_unit, *) (z0(k), ozrefs1(i, j, k),stds1(i, j, k), k = nlay, 1, -1) ! ppb                      
          ENDDO
        ENDDO
        CLOSE(atmos_unit)


        apfname = TRIM(ADJUSTL(atmdbdir)) //'tbclima/TB14H-5.vs' ! geun
        !apfname = '../../ATMOS/tbclima/TB14H-5.vs' ! jbak

        OPEN (UNIT = atmos_unit, file=apfname, status = 'unknown')
        READ (atmos_unit, '(A)') ;  READ(atmos_unit, '(A)')             
        DO i = 1, nmon
          DO j = nlat, 1, -1
           READ(atmos_unit, *) fdum, lats(j),mtropz2(i,j) 
           READ(atmos_unit, *) (z0(k), ozrefs2(i, j, k),stds2(i, j, k), k = nlay, 1, -1) ! ppb                      
         ENDDO
        ENDDO
        CLOSE(atmos_unit)
       
        ! extratropical TB: fill 5, 15, 25 
        fidx=minval(minloc(lats,mask = (lats(1:nlat) < 35 .and. lats(1:nlat) > 0)))
        lidx=minval(maxloc(lats,mask = (lats(1:nlat) < 35 .and. lats(1:nlat) > 0)))
        
        DO i = fidx, lidx
          ozrefs1(:,i, :) =ozrefs1(:, lidx+1, :)
          stds1(:, i, :)  =stds1(:, lidx+1, :)
          mtropz1(:,i)    =mtropz1(:,lidx+1)
        ENDDO   

        fidx=minval(minloc(lats,mask = (lats(1:nlat) < 0 .and. lats(1:nlat) > -35 )))
        lidx=minval(maxloc(lats,mask = (lats(1:nlat) < 0 .and. lats(1:nlat) > -35)))
        
        DO i = fidx, lidx
          ozrefs1(:,i, :) =ozrefs1(:, fidx-1, :)
          stds1(:, i, :)  =stds1(:, fidx-1, :)
          mtropz1(:,i)    =mtropz1(:,fidx-1)
        ENDDO   
      
  
        ! tropical TB: fill -85~-35 with -25, fill 35~85 with 35
        DO j = 1, 6
              ozrefs2(:, j, :) =ozrefs2(:, 7, :) !at -25
              stds2(:, j, :)   =stds2(:, 7, :)
              mtropz2(:,j)     =mtropz2(:,7)               

              ozrefs2(:, j+12, :)   =ozrefs1(:, 13, :) ! at 35
              stds2(:, j+12, :)     =stds2(:, 13, :)
              mtropz2(:,j+12)       =mtropz2(:,13)
        ENDDO
        IF (any(ozrefs1 < 0) .or. any(ozrefs <0) ) then
           print *, 'TB clima contain -999'  ; stop          
        ENDIF
        first = .FALSE.
ENDIF
        
if (which_tb == 1 ) then 
   ozrefs (:,:,:) = ozrefs1(:,:,:)
   stds (:,:,:) = stds1(:,:,:)
ELSE
   ozrefs (:,:,:) = ozrefs2(:,:,:)
   stds (:,:,:) = stds2(:,:,:)
ENDIF

! ** interpolation for lat, mon** ! 
  CALL get_monfrac(nmon, the_month, the_day, nm, monfrac, monin)
  CALL get_latfrac(nlat,latgrid, lat0,the_lat, nband, latfrac, latin)

 
  ozref0 =0.0 ; std0 = 0.0 ; mzt = 0.0
  DO i = 1, nband
        DO j = 1, nm
            ozref0 =  ozref0+ ozrefs(monin(j), latin(i), :) * monfrac(j) * latfrac(i)
            std0   =  std0+ stds(monin(j), latin(i), :) * monfrac(j) * latfrac(i)
            mzt    =  mzt+mtropz(monin(j), latin(i)) * monfrac(j) * latfrac(i)
        ENDDO
  ENDDO

!** convert tb reference into regular reference
!   zs0 = [ -20, 60] is data grid
!   tb  = reg - offset 
!   reg = tb + offset 
!   offset = (tropz - mzt)*(1-|reg(i)-tropz|/5)+mzt
!   reg(i) = tb(i) +tropz+(mtz-tropz)|reg(i)-tropz|/5   
!   convert PPB to DU [ here, just use constant shift ]
  tb0(0:nlay-1)   = z0-0.5 ; tb0(nlay)= 60 
  refz0(0:nlay)   = tb0(0:nlay)+trpz
  zstar(0:nlay)   = 1.0/(10**((refz0)/16.0))
  cum0(:) = 0.0 ; cums0(:)=0.0

  DO i = 1, nlay
   gravity_correct = (6367. / (6367. + refz0(i)+0.5 ))**2.
   ozref0(i) = ozref0(i)*( zstar(i-1)-zstar(i)) / ( 1.25 * gravity_correct)
   cum0(i)   = cum0(i-1) + ozref0(i)
   std0(i)   = std0(i)*( zstar(i-1)-zstar(i)) / ( 1.25 * gravity_correct)
   cums0(i)   = cums0(i-1) + std0(i)
  ENDDO

! derive variable shifht from LLM grid algitude covering 0 to 60 km
  refz(0:nref)   = (/(i*1.0, i= 0, 60 )/)   
  tb(0:nref)     = refz(0:nref)- trpz

!  DO i = 0, nref
!     IF ( abs(refz(i)-trpz ) <= 6. ) then 
!	 
!        offset(i) = (trpz-mzt)*(1-abs(refz(i)-trpz)/6.)+mzt
!		write(*,'(a5,10f8.4)')'v:',refz(i), offset(i), refz(i)-offset(i)
!     ELSE
!        offset(i) = mzt
!     ENDIF
!	 
!  ENDDO
!  tb(0:nref) = refz(0:nref)-offset(0:nref) 
!ENDIF

  IF (tb(0) < tb0(0) .or. tb(nref) > tb0(nlay) ) then 
      print * , 'check boundary condition in TB clim' 
      print * , TB(0), tb0(0), tb(nref), tb0(nlay) ; stop
  ENDIF

  CALL BSPLINE(tb0, cum0, nlay+1, tb, cum, nref+1, errstat)
  CALL BSPLINE(tb0, cums0, nlay+1, tb, cums, nref+1, errstat)
  IF (errstat < 0) THEN
    WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat ; stop
  ENDIF

  ozref(1:nref) = cum(1:nref)-cum(0:nref-1)
  std(1:nref)   = cums(1:nref)-cums(0:nref-1)
  IF (any(ozref(:) < 0)) then
     ozref(:) = -999 ; std(:) = -999 ; print *, 'TB <0' ; return
  endif
  CALL REVERSE(STD(1:nref), nref)

RETURN
END SUBROUTINE get_tb


SUBROUTINE get_ab (ozref,std, which_ab)

  USE OMSAO_precision_module 
  USE ozprof_data_module,     ONLY: atmos_unit
  USE OMSAO_variables_module, ONLY: atmdbdir, the_month, the_day, the_lat
  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, PARAMETER                             ::  nref = 60
  REAL (KIND=dp), DIMENSION(nref), INTENT(OUT)   ::  ozref,std
  INTEGER, INTENT(IN)                            ::  which_ab 
  ! ======================
  ! Local variables
  ! ======================
  CHARACTER (LEN=130)           :: apfname

  INTEGER, PARAMETER :: nlat=18, nmon=12, nlay=60
  REAL (KIND=dp), parameter::  latgrid=10., lat0=-90
  REAL (KIND=dp), DIMENSION(nlay)   :: ozref0,std0
  REAL (KIND=dp), DIMENSION(0:nlay) :: refz, zstar
  REAL (KIND=dp)                    :: gravity_correct

  REAL (KIND=dp)                :: frac,fdum
  REAL (KIND=dp), DIMENSION(2)  :: latfrac, monfrac
  INTEGER                       :: i, j, k, nband, nm, errstat
  INTEGER, DIMENSION(2)         :: latin, monin

  LOGICAL, SAVE                 :: first = .TRUE.
  REAL (KIND=dp), SAVE, DIMENSION(nmon, nlat, nlay) ::ozrefs1, ozrefs2, ozrefs
  REAL (KIND=dp), SAVE, DIMENSION(nmon, nlat, nlay) ::stds1, stds2, stds
  REAL (KIND=dp), SAVE, DIMENSION(nlay) :: zs0

  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=17), PARAMETER :: modulename = 'get_TB_NS_OZ'  
  

! ** load std profiles ** !
  IF (first) THEN
! AB clima. with all sonde profiles
    apfname = TRIM(ADJUSTL(atmdbdir)) // 'tbclima/ABall.ns' ! jbak
!    apfname = '/home/jbak/OzoneFit/OZBOREAS-OMI/src/tbclim/ABall.ns' ! jbak
    OPEN (UNIT = atmos_unit, file=apfname, status = 'unknown')
    READ (atmos_unit, '(A)') ;  READ(atmos_unit, '(A)') 
     
    DO i = 1, nmon
      DO j = nlat, 1, -1
        READ(atmos_unit, *) fdum, fdum
        READ(atmos_unit, *) (zs0(k), ozrefs1(i, j, k),stds1(i, j, k), k=nlay, 1, -1)! ppb
      ENDDO
    ENDDO
    CLOSE(atmos_unit)

! AB clima. with tropopause > 14.5   
    apfname = TRIM(ADJUSTL(atmdbdir)) // 'tbclima/AB14H.ns' ! jbak
!    apfname = '/home/jbak/OzoneFit/OZBOREAS-OMI/src/tbclim/AB14.5H.ns' ! jbak

    OPEN (UNIT = atmos_unit, file=apfname, status = 'unknown')
    READ (atmos_unit, '(A)') ;  READ(atmos_unit, '(A)') 
     
    DO i = 1, nmon
      DO j = nlat, 1, -1
        READ(atmos_unit, *) fdum, fdum
        READ(atmos_unit, *) (zs0(k), ozrefs2(i, j, k),stds2(i, j, k), k= nlay,1,-1)! ppb
      ENDDO
    ENDDO
    CLOSE(atmos_unit)

    DO j = 1, 6 ! filling mid/high with -/+25
      ozrefs2(:, j, :)   =ozrefs2(:, 7, :) ! -25
      stds2(:, j, :)     =stds2(:, 7, :)
    ENDDO
         
    DO j = 13, 18
      ozrefs2(1:5, j, :) =ozrefs2(1:5, 12, :) !25
      stds2(1:5, j, :)   =stds2(1:5, 12, :)
      ozrefs2(:, j, :) =ozrefs2(:, 13, :) !35
      stds2(:, j, :)   =stds2(:, 13, :)
    ENDDO
    IF (any(ozrefs2 < 0) .or. any(stds2 < 0)) then
      print *, 'AB clima contain -999'
    ENDIF

    first = .FALSE.
  ENDIF

  if (which_ab == 1 ) then 
    ozrefs (:,:,:) = ozrefs1(:,:,:)
    stds (:,:,:) = stds1(:,:,:)
  ELSE
    ozrefs (:,:,:) = ozrefs2(:,:,:)
    stds (:,:,:) = stds2(:,:,:)
  ENDIF

! ** interpolation for lat, mon** ! 

  CALL get_monfrac(nmon, the_month, the_day, nm, monfrac, monin)
  CALL get_latfrac(nlat,latgrid, lat0,the_lat, nband, latfrac, latin)

  ozref0(:) =0.0 ; std0(:)=0.0
  DO i = 1, nband
    DO j = 1, nm
      ozref0=  ozref0+ ozrefs(monin(j), latin(i), :) * monfrac(j) * latfrac(i)
      std0=  std0+ stds(monin(j), latin(i), :) * monfrac(j) * latfrac(i)
    ENDDO
  ENDDO

! ** convert ppb into DU ** ! 
  refz(0:nref) = (/(i*1.0, i= 0, 60 )/)
  zstar(0:nref) = 1.0/(10**((refz)/16.0))
  DO i = 1, nlay
    gravity_correct = (6367. / (6367. + refz(i)+0.5 ))**2.
    ozref0(i) = ozref0(i)*( zstar(i-1)-zstar(i)) / ( 1.25 * gravity_correct)
    std0(i) = std0(i)*( zstar(i-1)-zstar(i)) / ( 1.25 * gravity_correct)
  ENDDO
  
  ozref(:) = ozref0(:) 
  std(:)   = std0(:) 
  CALL REVERSE(STD(1:nref), nref)
 
RETURN
END SUBROUTINE get_ab

SUBROUTINE get_tjprof(out,index_out)

  USE OMSAO_precision_module 
  USE ozprof_data_module,     ONLY: atmos_unit,trpz
  USE OMSAO_variables_module, ONLY: atmdbdir, the_month, the_day, the_lat, the_lon
  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, PARAMETER                             ::  nref = 60
  INTEGER :: index_out
  REAL (KIND=dp), DIMENSION(nref), INTENT(OUT)   ::  out

  ! ======================
  ! Local variables
  ! ======================
  CHARACTER (LEN=130)           :: apfname

  INTEGER, PARAMETER :: nlat=60, nlon=72, nmon=12, nlay=20
  REAL (KIND=dp), PARAMETER         :: longrid = 5, latgrid = 3,lon0=-180.0, lat0=-90.0
  REAL (KIND=dp), DIMENSION(nref)   :: ozref0,std0, ozref, std
  REAL (KIND=dp), DIMENSION(0:nref) :: refz, zstar
  REAL (KIND=dp)                    :: gravity_correct, weight

  INTEGER                       :: i, j, k,im, nblat, nblon,nbmon, errstat
  REAL (KIND=dp), DIMENSION(2)  :: latfrac,lonfrac, monfrac
  INTEGER, DIMENSION(2)         :: latin,lonin, monin

  REAL (kind=dp), DIMENSION(nref) :: llm, llmstd
  LOGICAL, SAVE                 :: first = .TRUE.
  REAL (KIND=dp), SAVE, DIMENSION(nmon, nlon, nlat, nlay) ::ozrefs
  REAL (KIND=dp), SAVE, DIMENSION(nmon, nlon, nlat, nlay) ::stds

  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=17), PARAMETER :: modulename = 'get_tjprof'  

  IF (first) THEN
     apfname = TRIM(ADJUSTL(atmdbdir)) // 'tjclima/' //'tjclima_sealevel_o3ppb_mn.dat'
     OPEN (UNIT = atmos_unit, file = apfname, status = 'unknown')
     DO im = 1, nmon 
     READ (atmos_unit, '(72f8.3)') (((ozrefs(im, i, j,k), i=1, nlon), j=1, nlat), k=1,nlay)
     ENDDO
     CLOSE (atmos_unit)
     apfname = TRIM(ADJUSTL(atmdbdir)) // 'tjclima/' //'tjclima_sealevel_o3ppb_sd.dat'
     OPEN (UNIT = atmos_unit, file = apfname, status = 'unknown')
     DO im = 1, nmon 
      READ (atmos_unit, '(72f8.3)') (((stds(im, i, j,k), i=1, nlon), j=1, nlat), k=1,nlay)
     ENDDO
     CLOSE (atmos_unit)
     first = .FALSE.
  ENDIf

  CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
       the_lon, the_lat, nblon, nblat, lonfrac, latfrac, lonin, latin)

  CALL get_monfrac(nmon, the_month, the_day, nbmon, monfrac, monin)

  ozref0=0.0 ; std0 = 0.0
  DO i = 1, nblon
    DO j = 1, nblat
      DO im = 1, nbmon
          ozref0(1:nlay) =  ozref0(1:nlay)+ ozrefs(monin(im),lonin(i), latin(j), :) * monfrac(im) * latfrac(j)*lonfrac(i)
          std0(1:nlay)   =  std0(1:nlay)+ stds(monin(im),lonin(i), latin(j), :) * monfrac(im) * latfrac(j)*lonfrac(i)
      ENDDO
    ENDDO
  ENDDO

! ** convert ppm into DU ** ! 
  refz(0:nref) = (/(i*1.0, i= 0, nref )/)
  zstar(0:nlay) = 1.0/(10**((refz(0:nlay))/16.0))
  DO i = 1, nlay
   gravity_correct = (6367. / (6367. + refz(i)+0.5 ))**2.
   ozref0(i) = ozref0(i)*( zstar(i-1)-zstar(i)) / ( 1.25 * gravity_correct)
  ENDDO

  std0(1:nlay) = ozref0(1:nlay)*std0(1:nlay)/100.
  call get_mcprof(llm) 
  call get_mpstd(llmstd) 
  CALL REVERSE(llmstd(1:nref), nref)

! merging above trpz
  DO i = 1, nref 
     IF (refz(i) <= trpz-2.5) then 
         weight = 1.
     else if (refz(i) > trpz+2.5) then 
         weight = 0.
     else 
         weight = 1 - abs(refz(i) - trpz)/5.
     endif         
     if (ozref0(i) == 0.0) weight = 0.0
     ozref(i) = ozref0(i)*weight + LLM(i)*(1-weight)
     std0(i) = ozref(i)*0.1
     std(i) = std0(i)*weight + LLMstd(i)*(1-weight)
     !write(*, '(i2,5f8.2)') i,weight, ozref(i), ozref0(i), llm(i)
  ENDDO
  CALL REVERSE(std(1:nref), nref)

  IF ( index_out == 1 ) out(:) = ozref(:)
  IF ( index_out == 2 ) out(:) = std(:)
  RETURN
END SUBROUTINE get_tjprof

SUBROUTINE get_monfrac(nmon, mon, day, nbmon, monfrac, monin)

  USE OMSAO_precision_module
  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)                       :: nmon, mon, day
  INTEGER, INTENT(OUT)                      :: nbmon
  INTEGER, DIMENSION(2), INTENT(OUT)        :: monin
  REAL (KIND=dp), DIMENSION(2), INTENT(OUT) :: monfrac

  IF (day <= 15) THEN
     monin(1) = mon - 1
     IF (monin(1) == 0) monin(1) = 12
     monin(2) = mon
     monfrac(1) = (15.0 - day) / 30.0
     monfrac(2) = 1.0 - monfrac(1)
  ELSE 
     monin(2) = mon + 1
     IF (monin(2) == 13) monin(2) = 1
     monin(1) = mon
     monfrac(2) = (day - 15) / 30.0
     monfrac(1) = 1.0 - monfrac(2)
  ENDIF
     nbmon=2
END SUBROUTINE get_monfrac

SUBROUTINE get_latfrac( nlat, latgrid, lat0, lat,  nblat, latfrac, latin)

  USE OMSAO_precision_module
  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)                       :: nlat
  REAL (KIND=dp), INTENT(IN)                :: lat0, lat, latgrid
  INTEGER, INTENT(OUT)                      :: nblat
  INTEGER, DIMENSION(2), INTENT(OUT)        :: latin
  REAL (KIND=dp), DIMENSION(2), INTENT(OUT) :: latfrac
  
  ! ======================
  ! Local variables
  ! ======================  
  REAL (KIND=dp) :: frac, lat_offset
  
  lat_offset   = lat0 + latgrid / 2.0
  nblat = 2; frac = (lat - lat_offset) / latgrid + 1
  latin(1) = INT(frac); latin(2) = latin(1) + 1
  latfrac(1) = latin(2) - frac; latfrac(2) = 1.0 - latfrac(1)
  IF (latin(1) == 0)   THEN 
     latin(1) = 1;    latfrac(1) = 1.0; nblat = 1
  ENDIF
  
  IF (latin(2) > nlat) THEN
     latin(1) = nlat; latfrac(1) = 1.0; nblat = 1
  ENDIF


END SUBROUTINE get_latfrac

SUBROUTINE get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
     lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)

  USE OMSAO_precision_module
  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)                       :: nlon, nlat
  REAL (KIND=dp), INTENT(IN)                :: lon0, lat0, lat, lon, longrid, latgrid
  INTEGER, INTENT(OUT)                      :: nblon, nblat
  INTEGER, DIMENSION(2), INTENT(OUT)        :: latin, lonin
  REAL (KIND=dp), DIMENSION(2), INTENT(OUT) :: latfrac, lonfrac
  
  ! ======================
  ! Local variables
  ! ======================  
  REAL (KIND=dp) :: frac, lat_offset, lon_offset
  
  lat_offset   = lat0 + latgrid / 2.0
  lon_offset   = lon0 + longrid  / 2.0
  
  nblat = 2; frac = (lat - lat_offset) / latgrid + 1
  latin(1) = INT(frac); latin(2) = latin(1) + 1
  latfrac(1) = latin(2) - frac; latfrac(2) = 1.0 - latfrac(1)
  IF (latin(1) == 0)   THEN 
     latin(1) = 1;    latfrac(1) = 1.0; nblat = 1
  ENDIF
  
  IF (latin(2) > nlat) THEN
     latin(1) = nlat; latfrac(1) = 1.0; nblat = 1
  ENDIF
  
  ! Circular in longitude direction
  nblon = 2; frac = (lon - lon_offset) / longrid + 1
  lonin(1) = INT(frac); lonin(2) = lonin(1) + 1
  lonfrac(1) = lonin(2) - frac; lonfrac(2) = 1.0 - lonfrac(1)
  IF (lonin(1) == 0)   lonin(1) = nlon
  IF (lonin(2) > nlon) lonin(2) = 1
  
  RETURN
  
END SUBROUTINE get_gridfrac


SUBROUTINE get_gridfrac1(nlon, nlat, nmon, longrid, latgrid, mongrid, lon0, lat0, mon0, &
     lon, lat, mon, nblon, nblat, nbmon, lonfrac, latfrac, monfrac, lonin, latin, monin)

  USE OMSAO_precision_module
  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)                       :: nlon, nlat, nmon
  REAL (KIND=dp), INTENT(IN)                :: lon0, lat0, mon0, lat, lon, mon, longrid, latgrid, mongrid
  INTEGER, INTENT(OUT)                      :: nblon, nblat, nbmon
  INTEGER, DIMENSION(2), INTENT(OUT)        :: latin, lonin, monin
  REAL (KIND=dp), DIMENSION(2), INTENT(OUT) :: latfrac, lonfrac, monfrac
  
  ! ======================
  ! Local variables
  ! ======================  
  REAL (KIND=dp) :: frac, lat_offset, lon_offset, mon_offset
  
  lat_offset   = lat0   + latgrid / 2.0
  lon_offset   = lon0   + longrid / 2.0
  mon_offset   = mon0   + mongrid / 2.0
  
  nblat = 2; frac = (lat - lat_offset) / latgrid + 1
  latin(1) = INT(frac); latin(2) = latin(1) + 1
  latfrac(1) = latin(2) - frac; latfrac(2) = 1.0 - latfrac(1)
  IF (latin(1) == 0)   THEN 
     latin(1) = 1;    latfrac(1) = 1.0; nblat = 1
  ENDIF
  
  IF (latin(2) > nlat) THEN
     latin(1) = nlat; latfrac(1) = 1.0; nblat = 1
  ENDIF
  
  ! Circular in longitude direction
  nblon = 2; frac = (lon - lon_offset) / longrid + 1
  lonin(1) = INT(frac); lonin(2) = lonin(1) + 1
  lonfrac(1) = lonin(2) - frac; lonfrac(2) = 1.0 - lonfrac(1)
  IF (lonin(1) == 0)   lonin(1) = nlon
  IF (lonin(2) > nlon) lonin(2) = 1

  ! Circular in year
  nbmon = 2; frac = (mon - mon_offset) / mongrid + 1
  monin(1) = INT(frac); monin(2) = monin(1) + 1
  monfrac(1) = monin(2) - frac; monfrac(2) = 1.0 - monfrac(1)
  IF (monin(1) == 0)   monin(1) = nmon
  IF (monin(2) > nmon) monin(2) = 1
  
  RETURN
  
END SUBROUTINE get_gridfrac1

SUBROUTINE get_geoschem_o3mean(month, lon, lat, ps, ozprof, nz, ntp)
  
  USE OMSAO_precision_module 
  USE OMSAO_variables_module, ONLY: atmdbdir
  USE ozprof_data_module,     ONLY: atmos_unit
  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)         :: month, nz, ntp
  REAL (KIND=dp), INTENT(IN)  :: lat, lon
  REAL (KIND=dp), DIMENSION(0:nz), INTENT(IN)     :: ps
  REAL (KIND=dp), DIMENSION(nz),   INTENT(INOUT)  :: ozprof

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER               :: nlat=18, nlon=12, nalt=19
  INTEGER                          :: errstat, i, j, k, nblat, nblon, nalt0, ntp0
  REAL (KIND=dp), PARAMETER        :: latgrid=10.0, longrid=30.0, lon0=-180.0, lat0=-90.0

  INTEGER, DIMENSION(2)            :: latin, lonin
  REAL (KIND=dp), DIMENSION(2)     :: latfrac, lonfrac
  REAL (KIND=dp), DIMENSION(nalt)  :: gprof
  REAL (KIND=dp), DIMENSION(0:nz)  :: tempoz
  
  ! Saved variables
  REAL (KIND=dp), DIMENSION(nlon, nlat, nalt), SAVE :: geosoz
  LOGICAL                                    , SAVE :: first = .TRUE.

  REAL (KIND=dp), DIMENSION(0:nalt)           :: geospres, cumoz
  CHARACTER (LEN=3), DIMENSION(12)            :: months = (/'jan', 'feb',&
       'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/)
  CHARACTER (LEN=130)                         :: geosfile
  
  ! Correct coordinates
  REAL (KIND=DP), DIMENSION(0:nalt), PARAMETER:: pres = (/1.0d0,          &
       .987871d0, .954730d0, .905120d0, .845000d0, .78d0,    .710000d0,   &
       .639000d0, .570000d0, .503000d0, .440000d0,.380000d0, .325000d0,   &
       .278000d0, .237954d0, .202593d0, .171495d0, .144267d0, .121347d0,  &
       .102098d0/)
  
  IF (first) THEN
     geosfile = TRIM(ADJUSTL(atmdbdir)) // 'geoschem_tropclima/' // months(month) // '_o3_mean.dat'
     
     OPEN(UNIT = atmos_unit, FILE = geosfile, status='old')
     READ(atmos_unit, *) (((geosoz(i, j, k), k = 1, nalt), j = 1, nlat), i = 1, nlon)
     CLOSE (atmos_unit)
     first = .FALSE.
  ENDIF 

  CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
       lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)

  gprof = 0.0
  DO i = 1, nblon
     DO j = 1, nblat 
        gprof = gprof + geosoz(lonin(i), latin(j), :) * lonfrac(i) * latfrac(j)
     ENDDO
  ENDDO
     
  geospres = pres * ps(0)
  
  cumoz = 0.0
  DO i = 1, nalt
     cumoz(i) = cumoz(i-1) + gprof(i) * 1000.0 / 1.25 * &
          (geospres(i-1) - geospres(i)) / 1013.25 
     IF (ANY(geosoz(lonin(1:nblon), latin(1:nblat), i) <= 0.0)) THEN
        j = i - 1; EXIT
     ELSE
        j = i
     ENDIF
  ENDDO
  nalt0 = j

  DO i = 1, ntp
     IF (ps(i) < geospres(nalt0)) THEN
        ntp0 = i - 1; EXIT
     ENDIF
  ENDDO
     
  CALL BSPLINE(geospres, cumoz, nalt0+1, ps(0:ntp0), tempoz(0:ntp0), ntp0+1, errstat)
  tempoz(1:ntp0) = tempoz(1:ntp0) - tempoz(0:ntp0-1)     
  ozprof(1:ntp0) =  tempoz(1:ntp0) !* SUM(ozprof(1:ntp)) / SUM(tempoz(1:ntp)) *
  
  RETURN  
END SUBROUTINE get_geoschem_o3mean

SUBROUTINE get_geoschem_o3std(month, lon, lat, ps, ozprof, nz, ntp)
  
  USE OMSAO_precision_module 
  USE OMSAO_variables_module, ONLY: atmdbdir
  USE ozprof_data_module,     ONLY: atmos_unit
  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)         :: month, nz, ntp
  REAL (KIND=dp), INTENT(IN)  :: lat, lon
  REAL (KIND=dp), DIMENSION(0:nz), INTENT(IN)     :: ps
  REAL (KIND=dp), DIMENSION(nz),   INTENT(INOUT)  :: ozprof

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER               :: nlat=18, nlon=12, nalt=19
  INTEGER                          :: errstat, i, j, k, nblat, nblon, nalt0, ntp0
  REAL (KIND=dp), PARAMETER        :: latgrid=10.0, longrid=30.0, lon0=-180.0, lat0=-90.0

  INTEGER, DIMENSION(2)            :: latin, lonin
  REAL (KIND=dp), DIMENSION(2)     :: latfrac, lonfrac
  REAL (KIND=dp), DIMENSION(nalt)  :: gprof
  REAL (KIND=dp), DIMENSION(0:nz)  :: tempoz
  
  ! Saved variables
  REAL (KIND=dp), DIMENSION(nlon, nlat, nalt), SAVE :: geosoz
  LOGICAL                                    , SAVE :: first = .TRUE.

  REAL (KIND=dp), DIMENSION(0:nalt)           :: geospres, cumoz
  CHARACTER (LEN=3), DIMENSION(12)            :: months = (/'jan', 'feb',&
       'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/)
  CHARACTER (LEN=130)                         :: geosfile
  
  ! Correct coordinates
  REAL (KIND=DP), DIMENSION(0:nalt), PARAMETER:: pres = (/1.0d0,          &
       .987871d0, .954730d0, .905120d0, .845000d0, .78d0,    .710000d0,   &
       .639000d0, .570000d0, .503000d0, .440000d0,.380000d0, .325000d0,   &
       .278000d0, .237954d0, .202593d0, .171495d0, .144267d0, .121347d0,  &
       .102098d0/)
  
  IF (first) THEN
     geosfile = TRIM(ADJUSTL(atmdbdir)) // 'geoschem_tropclima/' // months(month) // '_o3_std.dat'
     
     OPEN(UNIT = atmos_unit, FILE = geosfile, status='old')
     READ(atmos_unit, *) (((geosoz(i, j, k), k = 1, nalt), j = 1, nlat), i = 1, nlon)
     CLOSE (atmos_unit)
     first = .FALSE.
  ENDIF

  CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
       lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)

  gprof = 0.0
  DO i = 1, nblon
     DO j = 1, nblat 
        gprof = gprof + geosoz(lonin(i), latin(j), :) * lonfrac(i) * latfrac(j)
     ENDDO
  ENDDO     
  geospres = pres * ps(0)
  
  cumoz = 0.0
  DO i = 1, nalt
     cumoz(i) = cumoz(i-1) + gprof(i) * 1000.0 / 1.25 * &
          (geospres(i-1) - geospres(i)) / 1013.25 
     IF (ANY(geosoz(lonin(1:nblon), latin(1:nblat), i) <= 0.0)) THEN
        j = i - 1; EXIT
     ELSE
        j = i
     ENDIF
  ENDDO
  nalt0 = j

  DO i = 1, ntp
     IF (ps(i) < geospres(nalt0)) THEN
        ntp0 = i - 1; EXIT
     ENDIF
  ENDDO
  
  CALL BSPLINE(geospres, cumoz, nalt+1, ps(0:ntp0), tempoz(0:ntp0), ntp0+1, errstat)
  tempoz(1:ntp0) = tempoz(1:ntp0) - tempoz(0:ntp0-1)    
  ozprof(1:ntp0) =  tempoz(1:ntp0) 
  
  RETURN  
END SUBROUTINE get_geoschem_o3std

! ==================================================================
! Obtain NCAR/NCEP 12pm surface pressure (mb) 
!   for  each 2.5 by 2.5 region
! If no data is available, then use mean surface
!    pressure from all years
! xliu: 03/08/11, switch from NCEP to 1x1 NCEP FNL (archived at NCAR)
! ==================================================================
SUBROUTINE get_spres(year, month, day, lon, lat, spres)

  USE OMSAO_precision_module 
  USE OMSAO_variables_module, ONLY: atmdbdir
  USE ozprof_data_module,     ONLY: atmos_unit
  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)           :: month, year, day
  REAL (KIND=dp), INTENT(IN)    :: lon, lat
  REAL (KIND=dp), INTENT(OUT)   :: spres

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER             :: nlat=180, nlon=360
  REAL (KIND=dp), PARAMETER      :: longrid = 1.0, latgrid = 1.0, lon0=-180.0, lat0=-90.0
  INTEGER                        :: i, j,  nblat, nblon
  LOGICAL                        :: file_exist
  CHARACTER (LEN=2)              :: monc, dayc
  CHARACTER (LEN=4)              :: yrc
  CHARACTER (LEN=130)            :: spres_fname
  INTEGER, DIMENSION(2)          :: latin, lonin
  REAL (KIND=dp), DIMENSION(2)   :: latfrac, lonfrac

  INTEGER, SAVE, DIMENSION(nlon, nlat) :: glbspres
  LOGICAL, SAVE                        :: first = .TRUE.

  IF (first) THEN
  
    WRITE(monc, '(I2.2)') month          ! from 9 to '09'
    WRITE(dayc, '(I2.2)') day ; WRITE(yrc,  '(I4.4)') year

    !spres_fname =TRIM(ADJUSTL(atmdbdir)) // 'nspres/spres' // yrc // monc // dayc // '.dat'
    spres_fname =TRIM(ADJUSTL(atmdbdir)) // 'fnl13.75LST/fnlsp/fnlsp_' // yrc // monc // dayc // '.dat'

! Determine if file exists or not
    INQUIRE (FILE= spres_fname, EXIST= file_exist)
    IF (.NOT. file_exist) THEN
      WRITE(*, *) 'Warning: no surface pressure file found, use monthly mean!!!'
      !spres_fname = TRIM(ADJUSTL(atmdbdir)) // 'nspres/spresavg' // monc // '.dat'
      spres_fname = TRIM(ADJUSTL(atmdbdir)) // 'fnl13.75LST/fnlsp/fnlspavg' // monc // '.dat'
    ENDIF
    OPEN (UNIT = atmos_unit, file = spres_fname, status = 'unknown')
    !READ (atmos_unit, '(144I4)') ((glbspres(i, j), i=1, nlon), j=1, nlat)
    READ (atmos_unit, '(360I4)') ((glbspres(i, j), i=1, nlon), j=1, nlat)
    CLOSE (atmos_unit)
    first = .FALSE.
  ENDIF

  CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
       lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
  spres = 0.0
  DO i = 1, nblon
     DO j = 1, nblat 
        spres = spres + glbspres(lonin(i), latin(j)) * lonfrac(i) * latfrac(j)
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE get_spres

SUBROUTINE get_sfct(year, month, day, lon, lat, sfct)

  USE OMSAO_precision_module 
  USE OMSAO_variables_module, ONLY: atmdbdir
  USE ozprof_data_module,     ONLY: atmos_unit
  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)           :: month, year, day
  REAL (KIND=dp), INTENT(IN)    :: lon, lat
  REAL (KIND=dp), INTENT(OUT)   :: sfct

  ! ======================
  ! Local variables
  ! ======================
  REAL (KIND=dp), PARAMETER      :: longrid = 1.0, latgrid = 1.0, lon0=-180.0, lat0=-90.0
  INTEGER, PARAMETER             :: nlon=180, nlat=360
  INTEGER                        :: i, j,  nblat, nblon
  LOGICAL                        :: file_exist
  CHARACTER (LEN=2)              :: monc, dayc
  CHARACTER (LEN=4)              :: yrc
  CHARACTER (LEN=130)            :: sfct_fname
  INTEGER, DIMENSION(2)          :: latin, lonin
  REAL (KIND=dp), DIMENSION(2)   :: latfrac, lonfrac

  INTEGER, SAVE, DIMENSION(nlon,nlat) :: glbsfct 
  LOGICAL, SAVE                        :: first = .TRUE.

  IF (first) THEN
    WRITE(monc, '(I2.2)') month          ! from 9 to '09' 
    WRITE(dayc, '(I2.2)') day            ! from 9 to '09'     
    WRITE(yrc,  '(I4.4)') year
     
    sfct_fname =TRIM(ADJUSTL(atmdbdir)) // 'fnl13.75LST/fnlst/fnlst_' // yrc // monc // dayc // '.dat' 

! Determine if file exists or not
    INQUIRE (FILE= sfct_fname, EXIST= file_exist)
    IF (.NOT. file_exist) THEN
      WRITE(*, *) 'Warning: no surface temperature file found, use monthly mean!!!'
      sfct_fname = TRIM(ADJUSTL(atmdbdir)) // 'fnl13.75LST/fnlst/fnlstavg' // monc // '.dat'
    ENDIF
    
    OPEN (UNIT = atmos_unit, file = sfct_fname, status = 'unknown')
    READ (atmos_unit, '(360I3)') ((glbsfct(i, j), i=1, nlon), j=1, nlat)
    CLOSE (atmos_unit)
    first = .FALSE.
  ENDIF
  CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
       lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
  sfct = 0.0
  DO i = 1, nblon
     DO j = 1, nblat 
        sfct = sfct + glbsfct(lonin(i), latin(j)) * lonfrac(i) * latfrac(j)
     ENDDO
  ENDDO
      
  RETURN
END SUBROUTINE get_sfct
 
! ====================================================================
! Obtain NCAR/NCEP 12pm tropopause pressure (mb) 
!   for  each 2.5 by 2.5 region
! If no data is available, then use mean surface
!    pressure from all years
! xliu: 03/08/11, switch from NCEP to 1x1 NCEP FNL (archived at NCAR)
 ! =================================================================
SUBROUTINE get_tpres(year, month, day, lon, lat, tpres)

  USE OMSAO_precision_module 
  USE OMSAO_variables_module, ONLY: atmdbdir
  USE ozprof_data_module,     ONLY: atmos_unit
  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)           :: month, year, day
  REAL (KIND=dp), INTENT(IN)    :: lon, lat
  REAL (KIND=dp), INTENT(OUT)   :: tpres

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER             :: nlat=180, nlon=360
  REAL (KIND=dp), PARAMETER      :: longrid = 1.0, latgrid = 1.0, lon0=-180.0, lat0=-90.0
  INTEGER                        :: i, j,  nblat, nblon
  LOGICAL                        :: file_exist
  CHARACTER (LEN=2)              :: monc, dayc
  CHARACTER (LEN=4)              :: yrc
  CHARACTER (LEN=130)            :: tpres_fname
  INTEGER, DIMENSION(2)          :: latin, lonin
  REAL (KIND=dp), DIMENSION(2)   :: latfrac, lonfrac


  INTEGER, SAVE, DIMENSION(nlon,nlat)  :: glbtpres
  LOGICAL, SAVE                        :: first = .TRUE.

  IF (first) THEN
    WRITE(monc, '(I2.2)') month          ! from 9 to '09' 
    WRITE(dayc, '(I2.2)') day            ! from 9 to '09'     
    WRITE(yrc,  '(I4.4)') year 

    tpres_fname =TRIM(ADJUSTL(atmdbdir)) // 'fnl13.75LST/fnltp/fnltp_' // yrc // monc // dayc // '.dat'
    

    ! Determine if file exists or not
    INQUIRE (FILE= tpres_fname, EXIST= file_exist)
    IF (.NOT. file_exist) THEN
      WRITE(*, *) 'Warning: no tropopause pressure file found, use monthly mean!!!'
      !tpres_fname = TRIM(ADJUSTL(atmdbdir)) // 'ntpres/tpresavg' // monc // '.dat'
      tpres_fname = TRIM(ADJUSTL(atmdbdir)) // 'fnl13.75LST/fnltp/fnltpavg' // monc // '.dat'
    ENDIF

    OPEN (UNIT = atmos_unit, file = tpres_fname, status = 'unknown')
    READ (atmos_unit, '(360I3)') ((glbtpres(i, j), i=1, nlon), j=1, nlat)
    CLOSE (atmos_unit)
    first = .FALSE.
  ENDIF
    
  CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
       lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
  tpres = 0.0
  DO i = 1, nblon
     DO j = 1, nblat 
        tpres = tpres + glbtpres(lonin(i), latin(j)) * lonfrac(i) * latfrac(j)
     ENDDO
  ENDDO
   
  RETURN
END SUBROUTINE get_tpres


! ===================================================
! Obtain EP TOMS monthly mean total ozone (DU) for
!    each 1.25 by 1 region
! If no data is available, then use mean total ozone 
!    from all years
! ====================================================
SUBROUTINE get_toz(year, month, day, lon, lat, toz)

  USE OMSAO_precision_module 
  USE OMSAO_variables_module, ONLY: atmdbdir
  USE ozprof_data_module,     ONLY: atmos_unit
  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)         :: month, year, day
  REAL (KIND=dp),INTENT(IN)   :: lon, lat
  REAL (KIND=dp), INTENT(OUT) :: toz

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER           :: nlat=180, nlon=288
  REAL (KIND=dp), PARAMETER    :: longrid = 1.25, latgrid = 1.0, lon0=-180.0, lat0=-90.0
  CHARACTER (LEN=2)            :: monc, yrc, dayc
  CHARACTER (LEN=130)          :: toz_fname

  INTEGER                      :: i, j, nblat, nblon
  LOGICAL                      :: file_exist
  INTEGER, DIMENSION(2)        :: latin, lonin
  REAL (KIND=dp), DIMENSION(2) :: latfrac, lonfrac

  INTEGER, SAVE, DIMENSION(nlon, nlat) :: glbtoz
  LOGICAL, SAVE                        :: first = .TRUE.

  IF (first) THEN
     WRITE(dayc, '(I2.2)') day             ! from 9 to '09'  
     WRITE(monc, '(I2.2)') month           ! from 9 to '09'  
     WRITE(yrc,  '(I2.2)') MOD(year, 100)  ! from 1997 to '97'

     toz_fname = TRIM(ADJUSTL(atmdbdir)) // 'eptoz/ep' // yrc // monc // '.dat'

     ! Determine if file exists or not
     INQUIRE (FILE= toz_fname, EXIST= file_exist)
     IF (.NOT. file_exist) THEN
        WRITE(*, *) 'Warning: no EP O3 file found, use monthly mean!!!'
        toz_fname = TRIM(ADJUSTL(atmdbdir)) // 'eptoz/avgep' // monc // '.dat'
     ENDIF

     OPEN (UNIT = atmos_unit, file=toz_fname, status = 'unknown')
     DO i = 1, 3
        READ (atmos_unit, '(A)')
     END DO
     READ (atmos_unit, *) ((glbtoz(i, j), i=1, nlon), j=1, nlat)
     CLOSE (atmos_unit)
     first = .FALSE.
  ENDIF

  CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
       lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
  toz = 0.0
  DO i = 1, nblon
     DO j = 1, nblat
        IF (glbtoz(lonin(i), latin(j)) > 0) &
             toz = toz + glbtoz(lonin(i), latin(j)) * lonfrac(i) * latfrac(j)
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE get_toz


! Obtain ECMWF temperature profile
SUBROUTINE get_ecmwft(year, month, day, lon, lat, ecmwft)

  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY: atmdbdir
  USE ozprof_data_module,     ONLY: atmos_unit
  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, PARAMETER          :: nlecm = 23
  INTEGER, INTENT(IN)         :: month, year, day
  REAL (KIND=dp), INTENT(IN)  :: lon, lat
  REAL (KIND=dp), DIMENSION(nlecm), INTENT(OUT) :: ecmwft

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER           :: nlat=72, nlon=144, nalt=23
  REAL (KIND=dp), PARAMETER    :: longrid = 2.5, latgrid = 2.5, lon0=-180.0, lat0=-90.0
  CHARACTER (LEN=2)            :: yrc, monc, dayc
  CHARACTER (LEN=130)          :: ecmwft_fname, ncep_fname
  INTEGER                      :: i, j, k, nblat, nblon
  INTEGER, DIMENSION(2)        :: latin, lonin
  REAL (KIND=dp), DIMENSION(2) :: latfrac, lonfrac
  LOGICAL                      :: file_exist

  INTEGER, SAVE, DIMENSION(nlon, nlat, nalt) :: glbecmwft
  LOGICAL, SAVE                              :: first = .TRUE.

  IF (first) THEN
     WRITE(monc, '(I2.2)') month          ! from 9 to '09' 
     WRITE(dayc, '(I2.2)') day            ! from 9 to '09'     
     WRITE(yrc,  '(I2.2)') MOD(year, 100) ! from 1997 to '97'

     ! use ECMWF
     IF (year <= 2001) THEN      
        ecmwft_fname = TRIM(ADJUSTL(atmdbdir)) // 'ecmwft/ecmwft' // yrc // monc // dayc // '.dat'      
        ! Determine if file exists or not
        INQUIRE (FILE= ecmwft_fname, EXIST= file_exist)
        IF (.NOT. file_exist) THEN
           WRITE(*, *) 'Warning: no T profile file found, use monthly mean!!!'
           ecmwft_fname = TRIM(ADJUSTL(atmdbdir)) // 'ecmwft/ecmwftavg' // monc // '.dat'
        ENDIF
        OPEN (UNIT = atmos_unit, file = ecmwft_fname, status = 'unknown')
        READ (atmos_unit, '(144i3)') (((glbecmwft(i, j, k), i=1, nlon), j=1, nlat), k=1, nalt)
        CLOSE(atmos_unit)
     ELSE  ! Use NCEP for up to 10 mb and ECMWFT average for up between 10 and 1 mb
        ! ECMWFT average between 10mb and 1mb (7, 5, 3, 2, 1)', other layers will be overlapped if no more data
        ecmwft_fname = TRIM(ADJUSTL(atmdbdir)) // 'ecmwft/ecmwftavg' // monc // '.dat'
        OPEN (UNIT = atmos_unit, file = ecmwft_fname, status = 'unknown')
        READ (atmos_unit, '(144i3)') (((glbecmwft(i, j, k), i=1, nlon), j=1, nlat), k=1, nalt)
        CLOSE(atmos_unit) 
        
        ncep_fname = TRIM(ADJUSTL(atmdbdir)) // 'ecmwft/ncep' // yrc // monc // dayc // '.dat'      
        ! Determine if file exists or not
        INQUIRE (FILE= ncep_fname, EXIST= file_exist)
        IF (.NOT. file_exist) THEN
           WRITE(*, *) 'Warning: no T profile file found, use monthly mean!!!'
           ! already read the monthly mean above
        ELSE
           OPEN (UNIT = atmos_unit, file = ncep_fname, status = 'unknown')
           ! NCEP misses the 775 level, which is shown in ECMWFT, other levels are the same
           READ(atmos_unit, '(144I3)') (((glbecmwft(i, j, k), i = 1, nlon), j = 1, nlat), k = 1, 3)
           READ(atmos_unit, '(144I3)') (((glbecmwft(i, j, k), i = 1, nlon), j = 1, nlat), k = 5, 18)
           glbecmwft(:, :, 4) = (glbecmwft(:, :, 3) + glbecmwft(:, :, 5)) / 2.0
        ENDIF
     ENDIF

     first = .FALSE.
  ENDIF

  CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
       lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
  ecmwft = 0.0
  DO i = 1, nblon
     DO j = 1, nblat
        ecmwft = ecmwft + glbecmwft(lonin(i), latin(j), :) * lonfrac(i) * latfrac(j)
     ENDDO
  ENDDO
  
  RETURN
END SUBROUTINE get_ecmwft

! Obtain NCEP temperature profile appended with ECMWF temperature profiles above 10 mb
SUBROUTINE get_ncept(year, month, day, lon, lat, ncept)

  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY: atmdbdir
  USE ozprof_data_module,     ONLY: atmos_unit
  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, PARAMETER                            :: nlecm=22
  INTEGER, INTENT(IN)                           :: month, year, day
  REAL (KIND=dp), INTENT(IN)                    :: lon, lat
  REAL (KIND=dp), DIMENSION(nlecm), INTENT(OUT) :: ncept

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER              :: nlat=72, nlon=144, nalt=22
  REAL (KIND=dp), PARAMETER       :: longrid = 2.5, latgrid = 2.5, lon0=-180.0, lat0=-90.0
  CHARACTER (LEN=2)               :: yrc, monc, dayc
  CHARACTER (LEN=130)             :: ecmwft_fname, ncep_fname
  INTEGER                         :: i, j, k, nblat, nblon
  INTEGER, DIMENSION(2)           :: latin, lonin
  REAL (KIND=dp), DIMENSION(2)    :: latfrac, lonfrac
  LOGICAL                         :: file_exist


  INTEGER, SAVE, DIMENSION(nlon, nlat, nalt) :: glbncept
  LOGICAL, SAVE                              :: first = .TRUE.

  IF (first) THEN
     WRITE(monc, '(I2.2)') month          ! from 9 to '09' 
     WRITE(dayc, '(I2.2)') day            ! from 9 to '09'     
     WRITE(yrc,  '(I2.2)') MOD(year, 100) ! from 1997 to '97'

     ! Use NCEP for up to 10 mb and ECMWFT average for up between 10 and 1 mb
     ! ECMWFT average between 10mb and 1mb (7, 5, 3, 2, 1)', 
     ! NCEP: 17 layers ECMWFT: 23 layers (including 7, 5, 3, 2, 1, 775 mb)
     ecmwft_fname = TRIM(ADJUSTL(atmdbdir)) // 'ecmwft/ecmwftavg' // monc // '.dat'
     OPEN (UNIT = atmos_unit, file = ecmwft_fname, status = 'unknown')
     ! nalt + 1 = 23
     READ (atmos_unit, '(144I3)') (((glbncept(i, j, k), i=1, nlon), j=1, nlat), k=1, 1)
     READ (atmos_unit, '(144I3)') (((glbncept(i, j, k), i=1, nlon), j=1, nlat), k=1, nalt)
     CLOSE(atmos_unit) 
        
     ncep_fname = TRIM(ADJUSTL(atmdbdir)) // 'nncept/ncep' // yrc // monc // dayc // '.dat'      
     ! Determine if file exists or not
     INQUIRE (FILE= ncep_fname, EXIST= file_exist)
     IF (.NOT. file_exist) THEN
        WRITE(*, *) 'Warning: no T profile file found, use monthly mean!!!'
        ncep_fname = TRIM(ADJUSTL(atmdbdir)) // 'nncept/ncepavg' // monc // '.dat'
     ENDIF
     OPEN (UNIT = atmos_unit, file = ncep_fname, status = 'unknown')
     READ(atmos_unit, '(144I3)') (((glbncept(i, j, k), i = 1, nlon), j = 1, nlat), k = 1, 17)

     first = .FALSE.
  ENDIF

  CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
       lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
  ncept = 0.0
  DO i = 1, nblon
     DO j = 1, nblat
        ncept = ncept + glbncept(lonin(i), latin(j), :) * lonfrac(i) * latfrac(j)
     ENDDO
  ENDDO
  
  RETURN
END SUBROUTINE get_ncept

! Obtain monthly mean ECMWF temperature for 7, 5, 3, 2, 1 mb 
SUBROUTINE get_ecmwfavgt(month, day, lon, lat, ecmwft)

  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY: atmdbdir
  USE ozprof_data_module,     ONLY: atmos_unit
  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, PARAMETER                            :: nl=5
  INTEGER, INTENT(IN)                           :: month, day
  REAL (KIND=dp), INTENT(IN)                    :: lon, lat
  REAL (KIND=dp), DIMENSION(nl), INTENT(OUT)    :: ecmwft

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER              :: nlat=72, nlon=144
  REAL (KIND=dp), PARAMETER       :: longrid = 2.5, latgrid = 2.5, lon0=-180.0, lat0=-90.0
  CHARACTER (LEN=2)               :: yrc, monc, dayc
  CHARACTER (LEN=130)             :: ecmwft_fname
  INTEGER                         :: il, i, j, k, nblat, nblon
  INTEGER, DIMENSION(2)           :: latin, lonin
  REAL (KIND=dp), DIMENSION(2)    :: latfrac, lonfrac
  LOGICAL                         :: file_exist


  INTEGER, SAVE, DIMENSION(nlon, nlat, nl) :: glbecmwft
  LOGICAL, SAVE                            :: first = .TRUE.

  IF (first) THEN
     WRITE(monc, '(I2.2)') month          ! from 9 to '09' 
     WRITE(dayc, '(I2.2)') day            ! from 9 to '09'     

     ! There are 23 layers in the data, but only read the last 
     ! five layers from ECMWF, i.e., at 7, 5, 3, 2, 1 mb, respectively
     ecmwft_fname = TRIM(ADJUSTL(atmdbdir)) // 'ecmwft/ecmwftavg' // monc // '.dat'
     OPEN (UNIT = atmos_unit, file = ecmwft_fname, status = 'unknown')
     DO il = 1, 18
        READ (atmos_unit, '(144I3)') (((glbecmwft(i, j, k), i=1, nlon), j=1, nlat), k=1, 1)
     ENDDO
     READ (atmos_unit, '(144I3)') (((glbecmwft(i, j, k), i=1, nlon), j=1, nlat), k=1, nl)
     CLOSE(atmos_unit) 
        
     first = .FALSE.
  ENDIF

  CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
       lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
  ecmwft = 0.0
  DO i = 1, nblon
     DO j = 1, nblat
        ecmwft = ecmwft + glbecmwft(lonin(i), latin(j), :) * lonfrac(i) * latfrac(j)
     ENDDO
  ENDDO
  
  RETURN
END SUBROUTINE get_ecmwfavgt

SUBROUTINE get_ncepfnlt(year, month, day, lon, lat, ncept)

  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY: atmdbdir
  USE ozprof_data_module,     ONLY: atmos_unit
  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, PARAMETER                            :: nl=26 
  INTEGER, INTENT(IN)                           :: month, year, day
  REAL (KIND=dp), INTENT(IN)                    :: lon, lat
  REAL (KIND=dp), DIMENSION(nl), INTENT(OUT)    :: ncept

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER              :: nlat=180, nlon=360
  REAL (KIND=dp), PARAMETER       :: longrid = 1.0, latgrid = 1.0, lon0=-180.0, lat0=-90.0
  CHARACTER (LEN=2)               :: monc, dayc
  CHARACTER (LEN=4)               :: yrc
  CHARACTER (LEN=130)             :: ncep_fname
  INTEGER                         :: i, j, k, nblat, nblon
  INTEGER, DIMENSION(2)           :: latin, lonin
  REAL (KIND=dp), DIMENSION(2)    :: latfrac, lonfrac
  LOGICAL                         :: file_exist

  INTEGER, SAVE, DIMENSION(nlon, nlat, nl)  :: glbncept  ! geun
  LOGICAL, SAVE                            :: first = .TRUE.

  IF (first) THEN
    WRITE(monc, '(I2.2)') month          ! from 9 to '09' 
    WRITE(dayc, '(I2.2)') day            ! from 9 to '09'     
    WRITE(yrc,  '(I4.4)') year           
 
    ncep_fname = TRIM(ADJUSTL(atmdbdir)) // 'fnl13.75LST/fnltemp/fnltemp_' // yrc // monc // dayc // '.dat'
      
    ! Determine if file exists or not
    INQUIRE (FILE= ncep_fname, EXIST= file_exist)
    IF (.NOT. file_exist) THEN
      WRITE(*, *) 'Warning: no T profile file found, use monthly mean!!!'
      ncep_fname = TRIM(ADJUSTL(atmdbdir)) // 'fnl13.75LST/fnltemp//fnltempavg' // monc // '.dat'
    ENDIF


     ! NCEP FNL: 26 layers (top down from 10 to 1000 mb), but data will be bottom up after being read     
     OPEN (UNIT = atmos_unit, file = ncep_fname, status = 'unknown')
     READ(atmos_unit, '(360I3)') (((glbncept(i, j, k), i = 1, nlon), j = 1, nlat), k = nl, 1, -1)

     first = .FALSE.
  ENDIF

  CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
       lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
  ncept = 0.0
  DO i = 1, nblon
     DO j = 1, nblat
        ncept = ncept + glbncept(lonin(i), latin(j), :) * lonfrac(i) * latfrac(j)
     ENDDO
  ENDDO
  
  RETURN
END SUBROUTINE get_ncepfnlt


SUBROUTINE get_mcprof(ozref)

  USE OMSAO_precision_module 
  USE OMSAO_variables_module, ONLY: atmdbdir, the_month, the_day, the_lat
  USE ozprof_data_module,     ONLY: atmos_unit,which_clima
  IMPLICIT NONE

  INTEGER, PARAMETER                           :: nref = 60
  ! ======================
  ! Input/Output variables
  ! ======================
  REAL (KIND=dp), DIMENSION(nref), INTENT(OUT) :: ozref
  
  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER :: nlat=18, nmon=12

  CHARACTER (LEN=130)                                :: ozprof_fname
  CHARACTER (LEN=200)                                :: line

  REAL (KIND=dp)                                           :: frac
  REAL (KIND=dp), DIMENSION(2)                             :: latfrac, monfrac
  INTEGER,        DIMENSION(2)                             :: latin, monin
  INTEGER :: i, j, im,ib,nband, nm 

  ! saved variables
  REAL (KIND=dp), SAVE, DIMENSION(nmon, nlat, nref)        :: ozrefs
  LOGICAL,        SAVE                                     :: first = .TRUE.

!  IF (which_clima == 1) print *,'A-priori set to be LLM'

! ** load oz profiles ** !
  IF (first) THEN
     ozprof_fname = TRIM(ADJUSTL(atmdbdir)) // 'mpclima/llmclima_prof.dat'
     OPEN (UNIT = atmos_unit, file= ozprof_fname, status = 'unknown')
     DO im = 1, nmon
        READ (atmos_unit, *)
        DO i = nref, 1, -1
           READ (atmos_unit, *) ozrefs(im, :, i)
        ENDDO
     ENDDO
     CLOSE (atmos_unit)     
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
  
  ozref = 0.0
  DO im = 1, nm
     DO ib = 1, nband
        ozref = ozref + latfrac(ib) * monfrac(im) * ozrefs(monin(im), latin(ib), :) 
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE get_mcprof

!  DU table 
!  lat [-85, 85]
!  mon [1, 12]
!  lat [0-1, 64-65, 66-90]
SUBROUTINE get_mlprof(out, index_out)

  USE OMSAO_precision_module 
  USE OMSAO_variables_module, ONLY: atmdbdir, the_month, the_day, the_lat
  USE ozprof_data_module,     ONLY: atmos_unit,which_clima
  IMPLICIT NONE

  INTEGER, PARAMETER                           :: nref = 60
  ! ======================
  ! Input/Output variables
  ! ======================

  INTEGER, INTENT(IN) :: index_out
  REAL (KIND=dp), DIMENSION(nref), INTENT(OUT) :: out
  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER :: nlat=18, nmon=12, nlay=66
  REAL (KIND=dp), parameter :: latgrid=10, lat0=-90

  CHARACTER (LEN=130)                                      :: ozprof_fname
  CHARACTER (LEN=200)                                      :: line
  CHARACTER (LEN=10)                                       :: cdum
  REAL (KIND=dp)                                           :: frac
  REAL (KIND=dp), DIMENSION(2)                             :: latfrac, monfrac
  INTEGER,        DIMENSION(2)                             :: latin, monin
  INTEGER :: i, j, im,ib,nband, nm  
  REAL (KIND=dp),DIMENSION(nlay)                            :: ozref0,std0, pres
  REAL (KIND=dp),DIMENSION(nref)                            :: std, ozref
  ! saved variables
  REAL (KIND=dp), SAVE, DIMENSION(nmon, nlat, nlay)        :: ozrefs, stds
  LOGICAL,        SAVE                                     :: first = .TRUE.
! ** load oz profiles ** !
  IF (first) THEN
      IF (which_clima == 12) print *,'A-priori set to be ML'
     ! LOAD ozone DU table
     ozprof_fname = TRIM(ADJUSTL(atmdbdir)) // 'MLclima/ML_du_table.dat'
     OPEN (UNIT = atmos_unit, file= ozprof_fname, status = 'unknown')

     DO ib = 1, nlat
        READ (atmos_unit, *) ;READ (atmos_unit, *)
        READ (atmos_unit, *) ;READ (atmos_unit, *)
        DO i =  1,nlay
           READ (atmos_unit, '(a10, 12f7.3)') cdum, ozrefs(:, ib, i) ! bottom-top
        ENDDO
     ENDDO
     CLOSE (atmos_unit)   

     ! LOAD STD ppmv table
     ozprof_fname = TRIM(ADJUSTL(atmdbdir)) //'MLclima/ML_ppmv_stats.dat'
     OPEN (UNIT = atmos_unit, file= ozprof_fname, status = 'unknown')
      READ (atmos_unit, *) ;READ (atmos_unit, *)
     DO ib = 1, nlat
        READ (atmos_unit, *) ;READ (atmos_unit, *)
        READ (atmos_unit, *)
        DO i =  nlay,1, -1
           READ (atmos_unit, '(a4, 12f7.3)') cdum, stds(:, ib, i) ! top-bottom
        ENDDO
     ENDDO
     CLOSE (atmos_unit)   
     first = .FALSE.    
  ENDIF

! ** interpolation for lat, mon** ! 

  CALL get_monfrac(nmon, the_month, the_day, nm, monfrac, monin)
  CALL get_latfrac(nlat,latgrid, lat0,the_lat, nband, latfrac, latin)

  ozref0 = 0.0; std0=0.0 
  DO im = 1, nm
     DO ib = 1, nband
        ozref0 = ozref0 + latfrac(ib) * monfrac(im) * ozrefs(monin(im), latin(ib), :) 
        std0 = std0 + latfrac(ib) * monfrac(im) * stds(monin(im), latin(ib), :) 
     ENDDO
  ENDDO
  ozref(:) = ozref0(1:nref)

! convert ppmb to DU for std profile
  pres(1) = 0.05  ! about 70 km
  pres(2:nlay) = (/(1013.25*10.0**(-1.0*i/16.0), i = nlay-2, 0, -1)/)
  
  DO i = nlay-nref, nlay-1
    std(i-5)= (std0(i+1) + std0(i))*0.5 *(pres(i+1) - pres(i))/ 1.267
  ENDDO

  IF ( index_out == 1 ) out(:) = ozref(:)
  IF ( index_out == 2 ) out(:) = std(:)
  

  RETURN
END SUBROUTINE get_mlprof



SUBROUTINE get_iupprof(toz, ozref)

  USE OMSAO_precision_module 
  USE OMSAO_variables_module, ONLY: atmdbdir, the_month, the_lat
  USE ozprof_data_module,     ONLY: atmos_unit,which_clima
  IMPLICIT NONE


  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, PARAMETER                           :: nref = 60
  REAL (KIND=dp), INTENT(IN)                   :: toz
  REAL (KIND=dp), DIMENSION(nref), INTENT(OUT) :: ozref

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER :: nlat=6, maxprof=14, nmon=2

  REAL (KIND=dp), DIMENSION(nref)   :: iupozref, llmozref, weight,refz
  REAL (KIND=dp), DIMENSION(14)                    :: tozindex, temp
  REAL (KIND=dp)                                   :: frac, fdum, maxoz, minoz, meg
  REAL (KIND=dp), DIMENSION(2)                     :: latfrac
  INTEGER,        DIMENSION(2)                     :: latin
  INTEGER                                          :: season

  CHARACTER (LEN=130)                              :: ozprof_fname
  INTEGER :: i, im, ib, il, k, iprof, profin, nprof, nband
  character(10), dimension(6) :: bandname =['90S-60S','60S-30S','30S-0S', '0N-30N', '30N-60N', '60N-90N'] 
    
  ! ======================
  ! saved variables
  INTEGER,        SAVE, DIMENSION(nmon, nlat)              :: nprofs
  REAL (KIND=dp), SAVE, DIMENSION(nmon, nlat,maxprof, nref):: ozrefs
  LOGICAL,        SAVE                                     :: first = .TRUE.
  ! ======================


! ** load oz profiles ** !
IF (which_clima == 5) print *,'A-priori set to be IUP with toz',toz     
IF (first) THEN

    ozprof_fname = TRIM(ADJUSTL(atmdbdir)) // 'iupclima/iupclima_o3du_mn.dat'
    OPEN (UNIT = atmos_unit, file= ozprof_fname, status = 'unknown')
        
    ! read loop        
    DO i = 1, 6 
         READ (atmos_unit, '(A)') ! read header
    ENDDO
        
    DO im = 1, nmon
        Do ib = 1,  nlat
           READ(atmos_unit, '(A)')   ! read month label         
           READ(atmos_unit,*) (tozindex(iprof), iprof=1, maxprof)! read total ozone label    

           DO il = nref, 1, -1 ! read bottom -> top
           READ(atmos_unit, *) fdum, (temp(iprof),iprof = 1, maxprof ) ! du     
           nprof = 1         
           DO k = 1, maxprof             
                  IF (tozindex(k) /= 0.0 ) then                   
                  ozrefs(im,ib,nprof, il) = temp(k)                        
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

  !!!!!! find latitude fraction  
     IF (the_lat <= -75.0) THEN
        nband = 1; latin(1) = 1; latfrac(1) = 1.0
     ELSE IF (the_lat >= 75.0) THEN
        nband = 1; latin(1) = nlat; latfrac(1) = 1.0
     ELSE
        nband = 2     ; frac = (the_lat + 75.0) / 30.0 + 1
        latin(1) = INT(frac); latin(2) = latin(1) + 1
        latfrac(1) = latin(2) - frac; latfrac(2) = 1.0 - latfrac(1)
     ENDIF
  !!!!!! find season
   season = 1                                    ! Dec-May
   IF (the_month > 5 .and. the_month < 12 ) season =2   ! July-Nov   

  !!!!!! find total ozone fraction
   iupozref(:) =0.   

   DO ib = 1, nband

           nprof = nprofs(season, latin(ib))
           minoz = sum( ozrefs(season, latin(ib), 1, :))
           maxoz = sum( ozrefs(season, latin(ib), nprof, :)) 
         
           IF (toz < minoz) THEN
              WRITE(*,*), 'Warning: no a priori profile available!!!'
               iupozref  = iupozref + ozrefs(season, latin(ib),1, :) * toz / minoz * latfrac(ib)
           ELSE IF (toz > maxoz) THEN

              WRITE(*,*), 'Warning: no a priori profile available!!!'
               iupozref =  iupozref + ozrefs(season, latin(ib), nprof, :) * toz / maxoz * latfrac(ib)
           ELSE

             profin = INT ((toz - minoz ) / 30.0) +1
                
             IF (profin == 0) THEN 
                 profin = 1
             ELSE IF (profin == nprof) THEN
                 profin = profin -1
             ENDIF
             frac = 1.0 - (toz - (minoz + (profin-1) * 30.0)) / 30.0

             iupozref = iupozref + latfrac(ib) * (frac * ozrefs(season, latin(ib),profin, :) &
                   + (1.0 - frac) * ozrefs(season, latin(ib),profin+1, : ))
           ENDIF
   ENDDO

  
  ! merging with LLM Clim
  refz(1:nref) = (/(i*1.0+0.5, i= 0, 59 )/)
  CALL get_mcprof(llmozref)
  meg = 25
  DO i = 1, nref 
        weight(i) = (meg+5-refz(i))/5.
     IF ( refz(i) <= meg )  weight(i) = 1.0
     IF ( refz(i) > meg+5 ) weight(i) = 0.0  
      ozref(i) = iupozref(i)*weight(i) + llmozref(i)*(1-weight(i)) 
  ENDDO

  RETURN      
END SUBROUTINE get_iupprof



! ===============================================================
! Obtain TOMS V8 ozone profiles (12 month, 18 latitude bands,
!   3-10 profiles with total ozone at a step of 50 DU
! ===============================================================
SUBROUTINE get_v8prof(toz, oz)

  USE OMSAO_precision_module 
  USE OMSAO_variables_module, ONLY: atmdbdir, the_month, the_day, the_lat
  USE ozprof_data_module,     ONLY: atmos_unit
  IMPLICIT NONE

  INTEGER, PARAMETER                           :: nl = 11
  ! ======================
  ! Input/Output variables
  ! ======================

  REAL (KIND=dp), INTENT(INOUT)                :: toz
  REAL (KIND=dp), DIMENSION(nl), INTENT(OUT)   :: oz

  
  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER :: nlat=18, maxprof=10, nmon=12

  CHARACTER (LEN=130)                                :: ozprof_fname
  CHARACTER (LEN=200)                                :: line

  REAL (KIND=dp)                                           :: frac, fdum, maxoz, minoz
  REAL (KIND=dp), DIMENSION(2)                             :: latfrac, monfrac
  INTEGER,        DIMENSION(2)                             :: latin, monin
  INTEGER :: i, j, ib, profin, nprof, nband, nm, im

  ! saved variables
  REAL (KIND=dp), SAVE, DIMENSION(nmon, nlat, maxprof, nl) :: ozprofs
  INTEGER,        SAVE, DIMENSION(nmon, nlat)              :: nprofs
  LOGICAL,        SAVE                                     :: first = .TRUE.




! ** load oz profiles ** !

  IF (first) THEN
  print * , 'EP+V8 profile '
  ozprof_fname = TRIM(ADJUSTL(atmdbdir)) // 'v8clima/tomsv8_ozone_clima.dat'
  OPEN (UNIT = atmos_unit, file= ozprof_fname, status = 'unknown')
        
        ! Read until the target month        
         DO im = 1, nmon
           DO i = 1, nlat 
              READ(atmos_unit, *) 
              nprof = 1
              DO j = 1, maxprof
                 READ (atmos_unit, '(A)') line;  READ (line, *) fdum

                 IF (fdum < 999.0) THEN
                    READ (line, *) fdum, ozprofs(im, i, nprof, :)
                    nprof = nprof + 1
                 ENDIF
              ENDDO
              nprofs(im, i) = nprof - 1              
           ENDDO
        ENDDO

  CLOSE (atmos_unit)
     
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
  

  oz = 0.0
  DO im = 1, nm
        DO ib = 1, nband   
           nprof = nprofs(monin(im), latin(ib))
           minoz = SUM(ozprofs(monin(im), latin(ib), 1, :))
           maxoz = SUM(ozprofs(monin(im), latin(ib), nprof, :))
                      
           IF (toz < minoz) THEN
              WRITE(*,*), 'Warning: no a priori profile available!!!'
              oz  = oz + ozprofs(monin(im), latin(ib), 1, :) * toz / minoz * latfrac(ib)
           ELSE IF (toz > maxoz) THEN
              WRITE(*,*), 'Warning: no a priori profile available!!!'
              oz = oz + ozprofs(monin(im), latin(ib), nprof, :) * toz / maxoz * latfrac(ib)
           ELSE
              profin = INT ((toz - minoz ) / 50.0)+1
              IF (profin == 0) THEN 
                 profin = 1
              ELSE IF (profin == nprof) THEN
                 profin = profin - 1
              ENDIF
              
              frac = 1.0 - (toz - (minoz + (profin-1) * 50.0)) / 50.0
              oz = oz + latfrac(ib) * monfrac(im) * (frac * ozprofs(monin(im), latin(ib), profin, :) &
                   + (1.0 - frac) * ozprofs(monin(im), latin(ib), profin+1, :))           
           ENDIF
        ENDDO
  ENDDO

  RETURN
END SUBROUTINE get_v8prof


! ===============================================================================
! Obtain TOMS V8 temperatire profiles (11, levels, 12 months, 18 latitude bands)
! ===============================================================================
SUBROUTINE get_v8temp(month, day, lat, v8temp)

  USE OMSAO_precision_module 
  USE OMSAO_variables_module, ONLY: atmdbdir
  USE ozprof_data_module,     ONLY: atmos_unit
  IMPLICIT NONE
  
  INTEGER, PARAMETER                              :: nl = 11
  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)                             :: month, day
  REAL (KIND=dp), INTENT(IN)                      :: lat
  REAL (KIND=dp), DIMENSION(nl), INTENT(OUT)      :: v8temp
  
  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER                              :: nlat=18, nmon=12
  CHARACTER (LEN=130)                             :: tfname

  ! saved variables
  REAL (KIND=dp), SAVE, DIMENSION(nl, nlat, nmon) :: tprofs
  LOGICAL,        SAVE                            :: first = .TRUE.
  REAL (KIND=dp)                                  :: frac
  REAL (KIND=dp), DIMENSION(2)                    :: latfrac, monfrac
  INTEGER,        DIMENSION(2)                    :: latin, monin
  INTEGER                                         :: ib, nb, nm, im

  IF (first) THEN
     ! read the reference profile for climatology
     tfname = TRIM(ADJUSTL(atmdbdir)) // 'v8clima/tv8_temp.dat'
     OPEN (UNIT = atmos_unit, file= tfname, status = 'unknown')
     READ  (atmos_unit, *) tprofs
     CLOSE (atmos_unit)     
     first = .FALSE.
  ENDIF
     
  IF (day <= 15) THEN
     monin(1) = month - 1
     IF (monin(1) == 0) monin(1) = 12
     monin(2) = month
     monfrac(1) = (15.0 - day) / 30.0
     monfrac(2) = 1.0 - monfrac(1)
  ELSE 
     monin(2) = month + 1
     IF (monin(2) == 13) monin(2) = 1
     monin(1) = month
     monfrac(2) = (day - 15) / 30.0
     monfrac(1) = 1.0 - monfrac(2)
  ENDIF
  nm = 2

  IF (lat <= -85.0) THEN
     nb = 1; latin(1) = 1; latfrac(1) = 1.0
  ELSE IF (lat >= 85.0) THEN
     nb = 1; latin(1) = nlat; latfrac(1) = 1.0
  ELSE
     nb = 2     ; frac = (lat + 85.0) / 10.0 + 1
     latin(1) = INT(frac); latin(2) = latin(1) + 1
     latfrac(1) = latin(2) - frac; latfrac(2) = 1.0 - latfrac(1)
  ENDIF
  
  v8temp = 0.0
  DO im = 1, nm
     DO ib = 1, nb
        v8temp = v8temp + latfrac(ib) * monfrac(im) * tprofs(:, latin(ib), monin(im)) 
     ENDDO
  ENDDO
  
  RETURN
END SUBROUTINE get_v8temp


! Use MIPAS IG2 Temperature Profile cimatology 
! 121 levels (pressre altitude from 120 km to 0 km), 4 months (1,4,7,10)
! and 6 latitude bands (-75, -45, -10, 10, 45, 75)
SUBROUTINE GET_MIPASIG2T(month, day, lat, xx, yy)
  
  USE OMSAO_precision_module 
  USE OMSAO_variables_module, ONLY: atmdbdir
  USE ozprof_data_module,     ONLY: atmos_unit
  IMPLICIT NONE
  
  INTEGER, PARAMETER                              :: nl = 121
  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)                             :: month, day
  REAL (KIND=dp), INTENT(IN)                      :: lat
  REAL (KIND=dp), DIMENSION(nl), INTENT(OUT)      :: xx, yy
  
  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER                              :: nlat=6, nmon=4
  REAL (KIND=dp), DIMENSION(1:nmon), PARAMETER    :: mons = (/0.5, 3.5, 6.5, 9.5/)
  REAL (KIND=dp), DIMENSION(1:nlat), PARAMETER    :: lats = (/-75.0, -45.0, -10.0, 10.0, 45.0, 75.0/)

  CHARACTER (LEN=130)                             :: fname

  ! saved variables
  REAL (KIND=dp), SAVE, DIMENSION(nl, nlat, nmon) :: profs
  REAL (KIND=dp), SAVE, DIMENSION(nl)             :: pres0
  LOGICAL,        SAVE                            :: first = .TRUE.

  REAL (KIND=dp), DIMENSION(0:nlat)               :: temp
  REAL (KIND=dp)                                  :: frac, fmon
  REAL (KIND=dp), DIMENSION(2)                    :: latfrac, monfrac
  INTEGER,        DIMENSION(2)                    :: latin, monin
  INTEGER                                         :: ib, nb, nm, im, i, nheader

  IF (first) THEN
     fname = TRIM(ADJUSTL(atmdbdir)) // 'mipasprof/MIPAS_IG2_Tclima.dat'
     nheader = 8
     
     OPEN (UNIT = atmos_unit, file= fname, status = 'unknown')
     DO i = 1, nheader
        READ (atmos_unit, *) 
     ENDDO
     
     DO im = 1, nmon
        READ (atmos_unit, *)
        DO i = nl, 1, -1
           READ (atmos_unit, *) temp(0:nlat)
           profs(i, 1:nlat, im) = temp(1:nlat)
        ENDDO
        READ (atmos_unit, *)        
     ENDDO

     DO i = 1, nl
        pres0(i) = 1013.25 * 10. ** (- (i - 1.0) / 16. )
     ENDDO

     CLOSE (atmos_unit)     
     first = .FALSE.
  ENDIF

  fmon = month - 1.0 + 1.0 * day / 31.
  IF (fmon <= mons(1)) THEN
     monin(1) = nmon; monin(2) = 1
     frac = 1.0 - (fmon + 2.5) / 3.0
     monfrac(1) = frac; monfrac(2) = 1.0 - frac
  ELSE IF (fmon >= mons(nmon) ) THEN
     monin(1) = nmon; monin(2) = 1
     frac = 1.0 - (fmon - mons(nmon)) / 3.0
     monfrac(1) = frac; monfrac(2) = 1.0 - frac
  ELSE
     DO i = 2, nmon
        IF (fmon < mons(i)) EXIT
     ENDDO

     monin(1) = i - 1; monin(2) = i
     frac     = 1.0 - (fmon - mons(i-1)) / (mons(i) - mons(i-1))
     monfrac(1) = frac; monfrac(2) = 1.0 - frac
  ENDIF
  nm = 2

  IF (lat <= lats(1)) THEN
     nb = 1; latin(1) = 1; latfrac(1) = 1.0
  ELSE IF (lat >= lats(nlat)) THEN
     nb = 1; latin(1) = nlat; latfrac(1) = 1.0
  ELSE
     DO i = 2, nlat
        IF (lat < lats(i)) EXIT
     ENDDO

     nb = 2; latin(1) = i - 1; latin(2) = i
     frac     = 1.0 - (lat - lats(i-1)) / (lats(i) - lats(i-1))
     latfrac(1) = frac; latfrac(2) = 1.0 - frac
  ENDIF
  
  xx = pres0; yy = 0.0
  DO im = 1, nm
     DO ib = 1, nb
        yy = yy + latfrac(ib) * monfrac(im) * profs(:, latin(ib), monin(im)) 
     ENDDO
  ENDDO

  RETURN

END SUBROUTINE GET_MIPASIG2T

! Use MIPAS IG2 Temperature Profile cimatology 
! 121 levels (pressre altitude from 120 km to 0 km), 4 months (1,4,7,10)
! and 6 latitude bands (-75, -45, -10, 10, 45, 75)

SUBROUTINE GET_MIPASIG2O3(month, day, lat, xx, yy)
  
  USE OMSAO_precision_module 
  USE OMSAO_variables_module, ONLY: atmdbdir
  USE ozprof_data_module,     ONLY: atmos_unit
  IMPLICIT NONE
  
  INTEGER, PARAMETER                              :: nl = 121
  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)                             :: month, day
  REAL (KIND=dp), INTENT(IN)                      :: lat
  REAL (KIND=dp), DIMENSION(nl), INTENT(OUT)      :: xx, yy
  
  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER                              :: nlat=6, nmon=4
  REAL (KIND=dp), DIMENSION(1:nmon), PARAMETER    :: mons = (/0.5, 3.5, 6.5, 9.5/)
  REAL (KIND=dp), DIMENSION(1:nlat), PARAMETER    :: lats = (/-75.0, -45.0, -10.0, 10.0, 45.0, 75.0/)

  CHARACTER (LEN=130)                             :: fname

  ! saved variables
  REAL (KIND=dp), SAVE, DIMENSION(nl, nlat, nmon) :: profs
  REAL (KIND=dp), SAVE, DIMENSION(nl)             :: pres0
  LOGICAL,        SAVE                            :: first = .TRUE.

  REAL (KIND=dp), DIMENSION(0:nlat)               :: temp
  REAL (KIND=dp)                                  :: frac, fmon
  REAL (KIND=dp), DIMENSION(2)                    :: latfrac, monfrac
  INTEGER,        DIMENSION(2)                    :: latin, monin
  INTEGER                                         :: ib, nb, nm, im, i, nheader

  IF (first) THEN
     fname = TRIM(ADJUSTL(atmdbdir)) // 'mipasprof/MIPAS_IG2_O3clima.dat'
     nheader = 9
     
     OPEN (UNIT = atmos_unit, file= fname, status = 'unknown')
     DO i = 1, nheader
        READ (atmos_unit, *) 
     ENDDO
     
     DO im = 1, nmon
        READ (atmos_unit, *)
        DO i = nl, 1, -1
           READ (atmos_unit, *) temp(0:nlat)
           profs(i, 1:nlat, im) = temp(1:nlat)
        ENDDO
        READ (atmos_unit, *)        
     ENDDO

     DO i = 1, nl
        pres0(i) = 1013.25 * 10. ** (- (i - 1.0) / 16. )
     ENDDO

     CLOSE (atmos_unit)     
     first = .FALSE.
  ENDIF

  fmon = month - 1.0 + 1.0 * day / 31.
  IF (fmon <= mons(1)) THEN
     monin(1) = nmon; monin(2) = 1
     frac = 1.0 - (fmon + 2.5) / 3.0
     monfrac(1) = frac; monfrac(2) = 1.0 - frac
  ELSE IF (fmon >= mons(nmon) ) THEN
     monin(1) = nmon; monin(2) = 1
     frac = 1.0 - (fmon - mons(nmon)) / 3.0
     monfrac(1) = frac; monfrac(2) = 1.0 - frac
  ELSE
     DO i = 2, nmon
        IF (fmon < mons(i)) EXIT
     ENDDO

     monin(1) = i - 1; monin(2) = i
     frac     = 1.0 - (fmon - mons(i-1)) / (mons(i) - mons(i-1))
     monfrac(1) = frac; monfrac(2) = 1.0 - frac
  ENDIF
  nm = 2

  IF (lat <= lats(1)) THEN
     nb = 1; latin(1) = 1; latfrac(1) = 1.0
  ELSE IF (lat >= lats(nlat)) THEN
     nb = 1; latin(1) = nlat; latfrac(1) = 1.0
  ELSE
     DO i = 2, nlat
        IF (lat < lats(i)) EXIT
     ENDDO

     nb = 2; latin(1) = i - 1; latin(2) = i
     frac     = 1.0 - (lat - lats(i-1)) / (lats(i) - lats(i-1))
     latfrac(1) = frac; latfrac(2) = 1.0 - frac
  ENDIF
  
  xx = pres0; yy = 0.0
  DO im = 1, nm
     DO ib = 1, nb
        yy = yy + latfrac(ib) * monfrac(im) * profs(:, latin(ib), monin(im)) 
     ENDDO
  ENDDO

  RETURN

END SUBROUTINE GET_MIPASIG2O3

SUBROUTINE get_surfalt(lon, lat, z0)

  USE OMSAO_precision_module 
  USE OMSAO_variables_module, ONLY: atmdbdir
  USE ozprof_data_module,     ONLY: atmos_unit
  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  REAL (KIND=dp), INTENT(IN)     :: lon, lat
  REAL (KIND=dp), INTENT(OUT)    :: z0

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER             :: nlat=360, nlon=720
  REAL (KIND=dp), PARAMETER      :: longrid = 0.5, latgrid = 0.5, lon0=-180.0, lat0=-90.0

  INTEGER                        :: i, j, nblat, nblon
  LOGICAL                        :: file_exist
  CHARACTER (LEN=130)            :: surfalt_fname
  INTEGER, DIMENSION(2)          :: latin, lonin
  REAL (KIND=dp), DIMENSION(2)   :: latfrac, lonfrac

  INTEGER, SAVE, DIMENSION(nlon, nlat) :: glbz
  LOGICAL, SAVE                        :: first = .TRUE.

  IF (first) THEN
      surfalt_fname = TRIM(ADJUSTL(atmdbdir)) // 'terrain_height/tomsv7_terrain.dat'
      
      ! Determine if file exists or not
      INQUIRE (FILE= surfalt_fname, EXIST= file_exist)
      IF (.NOT. file_exist) THEN
         STOP 'No Terrain Elevation datafile found!!!'
      ENDIF
      
      OPEN (UNIT = atmos_unit, file = surfalt_fname, status = 'unknown')
      DO i = 1, 4
         READ(atmos_unit, *)
      ENDDO
      
      READ (atmos_unit, '(720I4)') ((glbz(i, j), i=1, nlon), j=1, nlat)
      CLOSE (atmos_unit)
      first = .FALSE.
  ENDIF

  CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
       lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
  z0 = 0.0
  DO i = 1, nblon
     DO j = 1, nblat 
        z0 = z0 + glbz(lonin(i), latin(j)) * lonfrac(i) * latfrac(j)
     ENDDO
  ENDDO
  z0 = z0 / 1000.0  ! convert tp km
  
  RETURN
END SUBROUTINE get_surfalt

! xliu: 03/08/11, switch from NCEP to 1x1 NCEP FNL
SUBROUTINE get_ncepreso_surfalt(lon, lat, z0)

  USE OMSAO_precision_module 
  USE OMSAO_variables_module, ONLY: atmdbdir
  USE ozprof_data_module,     ONLY: atmos_unit
  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  REAL (KIND=dp), INTENT(IN)     :: lon, lat
  REAL (KIND=dp), INTENT(OUT)    :: z0

  ! ======================
  ! Local variables
  ! ======================
  !INTEGER, PARAMETER             :: nlat=72, nlon=144
  !REAL (KIND=dp), PARAMETER      :: longrid = 2.5, latgrid = 2.5, lon0=-180.0, lat0=-90.0
  INTEGER, PARAMETER             :: nlat=180, nlon=360
  REAL (KIND=dp), PARAMETER      :: longrid = 1.0, latgrid = 1.0, lon0=-180.0, lat0=-90.0

  INTEGER                        :: i, j, nblat, nblon
  LOGICAL                        :: file_exist
  CHARACTER (LEN=130)            :: surfalt_fname
  INTEGER, DIMENSION(2)          :: latin, lonin
  REAL (KIND=dp), DIMENSION(2)   :: latfrac, lonfrac

  INTEGER, SAVE, DIMENSION(nlon, nlat) :: glbz
  LOGICAL, SAVE                        :: first = .TRUE.

  IF (first) THEN
      !surfalt_fname = TRIM(ADJUSTL(atmdbdir)) // 'terrain_height/dem2.5x2.5.dat'
      surfalt_fname = TRIM(ADJUSTL(atmdbdir)) // 'terrain_height/fnlsh1x1.dat'
      
      ! Determine if file exists or not
      INQUIRE (FILE= surfalt_fname, EXIST= file_exist)
      IF (.NOT. file_exist) THEN
         STOP 'No Terrain Elevation datafile found!!!'
      ENDIF
      
      OPEN (UNIT = atmos_unit, file = surfalt_fname, status = 'unknown')
      DO i = 1, 4
         READ(atmos_unit, *)
      ENDDO
      
      !READ (atmos_unit, '(144I4)') ((glbz(i, j), i=1, nlon), j=1, nlat)
      READ (atmos_unit, '(360I4)') ((glbz(i, j), i=1, nlon), j=1, nlat)
      CLOSE (atmos_unit)
      first = .FALSE.
  ENDIF

  CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
       lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
  z0 = 0.0
  DO i = 1, nblon
     DO j = 1, nblat 
        z0 = z0 + glbz(lonin(i), latin(j)) * lonfrac(i) * latfrac(j)
     ENDDO
  ENDDO
  z0 = z0 / 1000.0  ! convert tp km
  
  RETURN
END SUBROUTINE get_ncepreso_surfalt

SUBROUTINE get_finereso_surfalt(lon, lat, z0)

  USE OMSAO_precision_module 
  USE OMSAO_variables_module, ONLY: atmdbdir
  USE ozprof_data_module,     ONLY: atmos_unit
  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  REAL (KIND=dp), INTENT(IN)     :: lon, lat
  REAL (KIND=dp), INTENT(OUT)    :: z0

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER             :: nlat=1800, nlon=3600
  REAL (KIND=dp), PARAMETER      :: longrid = 0.1, latgrid = 0.1, lon0=-180.0, lat0=-90.0

  INTEGER                        :: i, j, nblat, nblon
  LOGICAL                        :: file_exist
  CHARACTER (LEN=130)            :: surfalt_fname
  INTEGER, DIMENSION(2)          :: latin, lonin
  REAL (KIND=dp), DIMENSION(2)   :: latfrac, lonfrac

  INTEGER, SAVE, DIMENSION(nlon, nlat) :: glbz
  LOGICAL, SAVE                        :: first = .TRUE.

  IF (first) THEN
      surfalt_fname = TRIM(ADJUSTL(atmdbdir)) // 'terrain_height/dem0.1x0.1.dat'
      
      ! Determine if file exists or not
      INQUIRE (FILE= surfalt_fname, EXIST= file_exist)
      IF (.NOT. file_exist) THEN
         STOP 'No Terrain Elevation datafile found!!!'
      ENDIF
      
      OPEN (UNIT = atmos_unit, file = surfalt_fname, status = 'unknown')
      DO i = 1, 4
         READ(atmos_unit, *)
      ENDDO
      
      READ (atmos_unit, '(3600I4)') ((glbz(i, j), i=1, nlon), j=1, nlat)
      CLOSE (atmos_unit)
      first = .FALSE.
  ENDIF

  CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
       lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
  z0 = 0.0
  DO i = 1, nblon
     DO j = 1, nblat 
        z0 = z0 + glbz(lonin(i), latin(j)) * lonfrac(i) * latfrac(j)
     ENDDO
  ENDDO
  z0 = z0 / 1000.0  ! convert tp km
  
  RETURN
END SUBROUTINE get_finereso_surfalt


SUBROUTINE get_geoschem_o31(month, lon, lat, ps, ozprof, nz, ntp)  
  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY: atmdbdir 
  USE ozprof_data_module,     ONLY: atmos_unit

  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)         :: month, nz, ntp
  REAL (KIND=dp), INTENT(IN)  :: lon, lat
  REAL (KIND=dp), DIMENSION(0:nz), INTENT(IN)     :: ps
  REAL (KIND=dp), DIMENSION(nz),   INTENT(INOUT)  :: ozprof

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER               :: nlat=91, nlon=144, nalt=19
  REAL (KIND=dp), PARAMETER        :: longrid = 2.5, latgrid = 2.0, lon0=-181.25, lat0=-91.0
  INTEGER                          :: errstat, i, j, k, nblat, nblon, ntp0, nalt0

  REAL (KIND=dp), DIMENSION(nalt)  :: gprof
  REAL (KIND=dp), DIMENSION(0:nz)  :: tempoz
  INTEGER, DIMENSION(2)            :: latin, lonin
  REAL (KIND=dp), DIMENSION(2)     :: latfrac, lonfrac

  REAL (KIND=dp), SAVE, DIMENSION(nlon, nlat, nalt) :: geosoz
  LOGICAL, SAVE                                     :: first = .TRUE.

  ! Correct coordinates
  REAL (KIND=DP), DIMENSION(0:nalt), PARAMETER:: pres = (/1.0d0,          &
       .987871d0, .954730d0, .905120d0, .845000d0, .78d0, .710000d0,      &
       .639000d0, .570000d0, .503000d0, .440000d0,.380000d0, .325000d0,   &
       .278000d0, .237954d0, .202593d0, .171495d0, .144267d0, .121347d0,  &
       .102098d0/)

  REAL (KIND=dp), DIMENSION(0:nalt)           :: geospres, cumoz
  CHARACTER (LEN=3), DIMENSION(12) :: months = (/'jan', 'feb','mar', 'apr', &
       'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/)

  CHARACTER (LEN=130)      :: geosfile
  
  IF (first) THEN
     geosfile = TRIM(ADJUSTL(atmdbdir)) // 'geoschem_tropclima/' // months(month) // '_o3_avg.dat'
     OPEN(UNIT = atmos_unit, FILE = geosfile, status='old')
     READ(atmos_unit, *) (((geosoz(i, j, k), k = 1, nalt), j = 1, nlat), i = 1, nlon)
     CLOSE (atmos_unit)
     first = .FALSE.
  ENDIF

  CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
       lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
  gprof = 0.0
  DO i = 1, nblon
     DO j = 1, nblat 
        gprof = gprof + geosoz(lonin(i), latin(j), :) * lonfrac(i) * latfrac(j)
     ENDDO
  ENDDO 
  geospres = pres * ps(0)

  ! Integrate from ppm to DU  
  cumoz = 0.0
  DO i = 1, nalt  ! 1000.0 / 1.25 / 1013.25 = 1.0 / 1.2665625
     cumoz(i) = cumoz(i-1) + gprof(i) * (geospres(i-1) - geospres(i)) / 1.266525
     IF (ANY(geosoz(lonin(1:nblon), latin(1:nblat), i) <= 0.0)) THEN
        j = i - 1; EXIT
     ELSE
        j = i
     ENDIF
  ENDDO
  nalt0 = j

  DO i = 1, ntp
     IF (ps(i) < geospres(nalt0)) THEN
        ntp0 = i - 1; EXIT
     ENDIF
  ENDDO
  
  CALL BSPLINE(geospres, cumoz, nalt0+1, ps(0:ntp0), tempoz(0:ntp0), ntp0+1, errstat)
  tempoz(1:ntp0) = tempoz(1:ntp0) - tempoz(0:ntp0-1)     
  ozprof(1:ntp0) =  tempoz(1:ntp0) !* SUM(ozprof(1:ntp)) / SUM(tempoz(1:ntp)) *
  ! use profile shape only
  !ozprof(1:ntp) =  tempoz(1:ntp) * SUM(ozprof(1:ntp)) / SUM(tempoz(1:ntp)) 
  
  RETURN  
END SUBROUTINE GET_GEOSCHEM_O31

SUBROUTINE get_logan_clima(month, lon, lat, ps, ozprof, nz, ntp)  
  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY: atmdbdir 
  USE ozprof_data_module,     ONLY: atmos_unit
  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)                             :: month, nz, ntp
  REAL (KIND=dp), INTENT(IN)                      :: lon, lat
  REAL (KIND=dp), DIMENSION(0:nz), INTENT(IN)     :: ps
  REAL (KIND=dp), DIMENSION(nz),   INTENT(INOUT)  :: ozprof

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER        :: nlat=46, nlon=72, nalt=13
  REAL (KIND=dp), PARAMETER :: longrid = 5.0, latgrid = 4.0, lon0=-180.0, lat0=-92.0
  INTEGER                   :: errstat, i, j, k, nblat, nblon, ntp0

  REAL (KIND=dp), DIMENSION(nalt)             :: gprof
  REAL (KIND=dp), DIMENSION(0:nz)             :: tempoz
  INTEGER, DIMENSION(2)                       :: latin, lonin
  REAL (KIND=dp), DIMENSION(2)                :: latfrac, lonfrac

  REAL (KIND=dp), SAVE, DIMENSION(nlon, nlat, nalt) :: geosoz
  LOGICAL, SAVE                                     :: first = .TRUE.

  ! Correct coordinates
  REAL (KIND=DP), DIMENSION(1:nalt), PARAMETER:: pres = (/1000., 900., &
       800., 700., 600., 500., 400., 300., 250., 200., 150., 125., 100./)

  REAL (KIND=dp), DIMENSION(1:nalt)           :: cumoz, presmod
  CHARACTER (LEN=3), DIMENSION(12)            :: months = (/'jan', 'feb', &
       'mar', 'apr',  'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/)
  CHARACTER (LEN=2)                           :: monc
  CHARACTER (LEN=130)                         :: geosfile
  
  IF (first) THEN
     WRITE(monc, '(I2.2)') month
     geosfile = TRIM(ADJUSTL(atmdbdir)) // 'logan_clima/ozone.13.4x5.' // monc

     OPEN(UNIT = atmos_unit, FILE = geosfile, status='old')
     READ(atmos_unit, '(9E10.3)') geosoz
     CLOSE (atmos_unit)
     first = .FALSE.
  ENDIF

  CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
       lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
  gprof = 0.0
  DO i = 1, nblon
     DO j = 1, nblat 
        gprof = gprof + geosoz(lonin(i), latin(j), :) * lonfrac(i) * latfrac(j)
     ENDDO
  ENDDO 
  
  ! Integrate from ppb to DU
  cumoz = 0.0
  DO i = 2, nalt  ! 2533.125 = 2 * 1.25 * 1013.25
     cumoz(i) = cumoz(i-1) + (gprof(i-1) + gprof(i)) * (pres(i-1) - pres(i)) / 2533.125 
  ENDDO

  presmod = pres
  DO i = 1, ntp
     IF (ps(i) < presmod(nalt)) THEN
        ntp0 = i - 1; EXIT
     ENDIF
  ENDDO
  IF (presmod(1) < ps(0))  presmod(1) = ps(0)
  
  CALL BSPLINE(presmod, cumoz, nalt, ps(0:ntp0), tempoz(0:ntp0), ntp0+1, errstat)
  tempoz(1:ntp0) = tempoz(1:ntp0) - tempoz(0:ntp0-1)    
  ozprof(1:ntp0) =  tempoz(1:ntp0)  ! use actual profile shape

  ! use profile shape only
  ! ozprof(1:ntp) =  tempoz(1:ntp) * SUM(ozprof(1:ntp)) / SUM(tempoz(1:ntp)) 
 
  RETURN  
END SUBROUTINE GET_LOGAN_CLIMA

! ps: bottom up
! read GEOS-STRAT (V6.13) from May Fu
SUBROUTINE GET_GEOSCHEM_HCHO(month, lon, lat, ps, hcho, nz)  
  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: du2mol
  USE OMSAO_variables_module, ONLY: atmdbdir 
  USE ozprof_data_module,     ONLY: atmos_unit

  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)                          :: month, nz
  REAL (KIND=dp), INTENT(IN)                   :: lon, lat
  REAL (KIND=dp), DIMENSION(0:nz), INTENT(IN)  :: ps
  REAL (KIND=dp), DIMENSION(nz),   INTENT(OUT) :: hcho  ! in Dobson Units

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER        :: nlat=91, nlon=144, nalt=19
  REAL (KIND=dp), PARAMETER :: longrid = 2.5, latgrid = 2.0, lon0=-181.25, lat0=-91.0
  INTEGER                   :: errstat, i, j, k, nblat, nblon, ntp

  REAL (KIND=dp), DIMENSION(nalt)             :: gprof
  REAL (KIND=dp), DIMENSION(0:nalt)           :: geospres, cumhcho
  REAL (KIND=dp), DIMENSION(0:nz)             :: temphcho
  INTEGER, DIMENSION(2)                       :: latin, lonin
  REAL (KIND=dp), DIMENSION(2)                :: latfrac, lonfrac

  REAL (KIND=dp), SAVE, DIMENSION(nlon, nlat, nalt) :: geoshcho
  LOGICAL, SAVE                                     :: first = .TRUE.

  ! Correct coordinates
  REAL (KIND=DP), DIMENSION(0:nalt), PARAMETER:: pres = (/1.0d0,          &
       .987871d0, .954730d0, .905120d0, .845000d0, .78d0, .710000d0,      &
       .639000d0, .570000d0, .503000d0, .440000d0,.380000d0, .325000d0,   &
       .278000d0, .237954d0, .202593d0, .171495d0, .144267d0, .121347d0,  &
       .102098d0/)
  CHARACTER (LEN=3), DIMENSION(12)  :: months = (/'jan', 'feb','mar', 'apr', &
       'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/)
  CHARACTER (LEN=130)               :: geosfile
  
  IF (first) THEN
     geosfile = TRIM(ADJUSTL(atmdbdir)) // 'geoschem_hcho/' // months(month) // '_hcho_avg.dat'
     OPEN(UNIT = atmos_unit, FILE = geosfile, status='old')
     DO i = 1, 5
        READ(atmos_unit, *) 
     ENDDO
     READ(atmos_unit, *) (((geoshcho(i, j, k), k = 1, nalt), j = 1, nlat), i = 1, nlon)
     CLOSE (atmos_unit)
     first = .FALSE.
  ENDIF

  CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
       lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
  gprof = 0.0
  DO i = 1, nblon
     DO j = 1, nblat 
        gprof = gprof + geoshcho(lonin(i), latin(j), :) * lonfrac(i) * latfrac(j)
     ENDDO
  ENDDO 
  geospres = pres * ps(0)

  ! Integrate from ppb to DU  
  cumhcho = 0.0
  DO i = 1, nalt  ! 1266.5625 = 1.25 * 1013.25
     cumhcho(i) = cumhcho(i-1) + gprof(i) * (geospres(i-1) - geospres(i)) / 1266.5625
  ENDDO

  ! MAXLOC (maximum value), the index of maxloc starts from 1
  ntp = MINVAL(MINLOC(ps(0:nz), MASK = (ps(0:nz) >= geospres(nalt)))) - 1
  
  hcho = 0.0
  CALL INTERPOL(geospres, cumhcho, nalt+1, ps(0:ntp), temphcho(0:ntp), ntp+1, errstat)
  hcho(1:ntp) = temphcho(1:ntp) - temphcho(0:ntp-1)   ! DU at each layer 
  
  hcho(1:ntp) = hcho(1:ntp) * du2mol

  RETURN  
END SUBROUTINE GET_GEOSCHEM_HCHO


! zs, ps: bottom up (BOS -> TOA)
! read stratospheric BrO from Chris + 0.2 ppbv in the troposphere (too small)
SUBROUTINE GET_BRO(month, lat, zs, ps, sza, bro, nz)  

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: deg2rad, du2mol
  USE OMSAO_variables_module,  ONLY: atmdbdir 
  USE ozprof_data_module,      ONLY: atmos_unit
  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)                          :: month, nz
  REAL (KIND=dp),                  INTENT(IN)  :: lat
  REAL (KIND=dp), INTENT(IN)                   :: sza
  REAL (KIND=dp), DIMENSION(0:nz), INTENT(IN)  :: zs, ps
  REAL (KIND=dp), DIMENSION(nz),   INTENT(OUT) :: bro     ! In Dobson Units

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER          :: nlat = 18, nalt = 25, maxsza = 17, nmax=601
  REAL (KIND=dp), PARAMETER   :: latgrid = 10.0
  INTEGER                     :: errstat, i, j, nsza, fidx, lidx, nband
  INTEGER, DIMENSION(2)       :: latin
  REAL(KIND=dp), DIMENSION(2) :: latfrac

  INTEGER,        SAVE, DIMENSION(nlat)               :: nszas
  REAL (KIND=dp), SAVE, DIMENSION(nlat, nalt, maxsza) :: allbro
  REAL (KIND=dp), SAVE, DIMENSION(nlat, nalt)         :: alts
  REAL (KIND=dp), SAVE, DIMENSION(nlat, maxsza)       :: szas
  REAL (KIND=dp), SAVE, DIMENSION(nmax)               :: ptmp0, brotmp, ptmp
  INTEGER, SAVE                                       :: ntmp
  LOGICAL, SAVE                                       :: first = .TRUE.

  REAL (KIND=dp), DIMENSION(maxsza)           :: tempszas
  REAL (KIND=dp), DIMENSION(nalt)             :: bprof
  REAL (KIND=dp), DIMENSION(nalt+5)           :: temprof, tempalt      ! Stuff to the surface
  REAL (KIND=dp), DIMENSION(0:nz)             :: tempbro, temp
  REAL (KIND=dp)                              :: frac, csza, airdens

  INTEGER :: which_bro = 0  ! 0: from PRATMO, 1: from SAO  2: from GEOS-5


  CHARACTER (LEN=3), DIMENSION(12)  :: months = (/'jan', 'feb','mar', 'apr', &
       'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/)
  CHARACTER (LEN=130)               :: fname
  
  IF (first) THEN
     fname = TRIM(ADJUSTL(atmdbdir)) // 'BRO/' // months(month) // '_atm_BrO.dat' 
     OPEN(UNIT = atmos_unit, FILE = fname, status='old')
     
     DO i = 1, nlat 
        READ(atmos_unit, '(5X,2I5)') nszas(i), nszas(i)
        ! Only use AM data (from mid-night to noon)
        nszas(i) = nszas(i) / 2 
        READ(atmos_unit, '(5X, 34E13.4)') szas(i, 1:nszas(i)), szas(i, 1:nszas(i))
        DO j = nalt, 1, -1  ! Read top-down (j = 1, lowest level) 
           READ(atmos_unit, *) alts(i, j), allbro(i, j, 1:nszas(i)), &
                allbro(i, j, 1:nszas(i)), airdens
           allbro(i, j, 1:nszas(i)) = allbro(i, j, 1:nszas(i)) / airdens * 1.0E9 ! to ppbv           
        ENDDO
        READ(atmos_unit, *)
     ENDDO     
     CLOSE (atmos_unit)
     
     IF (which_bro == 1) THEN
        fname = TRIM(ADJUSTL(atmdbdir)) // 'BRO/' // 'omsao_bro_oper.dat' 
        OPEN(UNIT = atmos_unit, FILE = fname, status='old')     
        READ(atmos_unit, *) ntmp
        DO i = 1, ntmp
           READ(atmos_unit, *) ptmp0(i), brotmp(i)
        ENDDO
        CLOSE(atmos_unit)
     ELSE IF (which_bro == 2) THEN
        fname = TRIM(ADJUSTL(atmdbdir)) // 'BRO/' // 'geos5_bro_barr_2008m0406.dat' 
        OPEN(UNIT = atmos_unit, FILE = fname, status='old')     
        READ(atmos_unit, *) ntmp
        DO i = 1, ntmp
           READ(atmos_unit, *) ptmp0(i), brotmp(i)
        ENDDO
        CLOSE(atmos_unit)
     ENDIF

     first = .FALSE.
  ENDIF
  
  ! Interpolate over latitude (avoid discontinuity) 
  IF (lat <= -85.0) THEN
     nband = 1; latin(1) = 1; latfrac(1) = 1.0
  ELSE IF (lat >= 85.0) THEN
     nband = 1; latin(1) = nlat; latfrac(1) = 1.0
  ELSE
     nband = 2; frac = (lat + 85.0) / latgrid + 1
     latin(1) = INT(frac); latin(2) = latin(1) + 1
     latfrac(1) = latin(2) - frac; latfrac(2) = 1.0 - latfrac(1)
  ENDIF  
  
  tempalt(1:5) = (/0.0, 2.0, 4.0, 6.0, 8.0/)
  ! temprof(1:5) = 2.0E-4   ! Assume 0.2 pptv in the troposphere (from S.E. Chris)
  ! Now assume constant mixing ratio with lowest level mixing ratio
  ! This will better capture seasonal and latitudinal variation of the tropospheric BrO

  tempbro = 0.0
  csza = COS(sza * deg2rad)
  DO i = 1, nband
     nsza = nszas(latin(i))
     tempszas(1:nsza) = COS(szas(latin(i), 1:nsza) * deg2rad)

     fidx = MINVAL(MAXLOC( tempszas(1:nsza), MASK=(tempszas(1:nsza) <= csza)))
     IF (fidx == nsza) THEN
        bprof = allbro(latin(i), :, nsza)
     ELSE IF (fidx == 0) THEN
        bprof = allbro(latin(i), :, 1)
     ELSE
        lidx = fidx + 1
        frac = 1.0 - (csza - tempszas(fidx)) / (tempszas(lidx) - tempszas(fidx))
        bprof = frac * allbro(latin(i), :, fidx) + (1 - frac) * allbro(latin(i), :, lidx)
     ENDIF
     
     tempalt(6:nalt+5) = alts(latin(i), :);     temprof(6:nalt+5) = bprof
     temprof(1:5) = 0. !temprof(6)  ! Extratoplate to the whole troposphere
     IF (tempalt(nalt+5) < zs(nz)) tempalt(nalt+5) = zs(nz)
     IF (tempalt(1) > zs(0)) tempalt(1) = zs(0)  ! Avoid extrapolation
     
     ! Interpolate to GOME retrieval altitudes
     CALL BSPLINE(tempalt, temprof, nalt+5, zs(0:nz), temp(0:nz), nz+1, errstat)
     IF (temprof(5) == 0.0) THEN
        WHERE (zs(0:nz) < tempalt(6))
           temp(0:nz) = 0.0
        ENDWHERE
     ENDIF   
     tempbro = tempbro + temp * latfrac(i)     
  ENDDO

  IF (which_bro /= 0) THEN
     ptmp(1:ntmp) = ptmp0(1:ntmp)
     IF (ptmp0(1) < ps(0)) ptmp(1) = ps(0)
     IF (ptmp0(ntmp) > ps(nz)) ptmp(ntmp) = ps(nz)  
     CALL BSPLINE(ptmp(1:ntmp), brotmp(1:ntmp), ntmp, ps(0:nz), tempbro(0:nz), nz+1, errstat)
  ENDIF

  !WRITE(90, *) nz + 1
  !WRITE(90, '(2D14.5)') ((ps(i), tempbro(i)), i=0, nz)
  !STOP
 
  ! Integrate from ppbv to DU  
  DO i = 1, nz
     ! accurate to within 1% (2533.125 = 2 * 1.25 * 1013.25
     bro(i) = (tempbro(i) + tempbro(i-1)) * (ps(i-1) - ps(i)) / 2533.125
  ENDDO
  bro(1:nz) = bro(1:nz) * du2mol
     
  RETURN  
END SUBROUTINE GET_BRO


! zs, ps: bottom up (BOS -> TOA)
! read GEOS-STRAT (V6.13) from May Fu + stratospheric NO2 from Chris
! CRN: After 2006, use NO2 from Lok Lamsal produced for SCIA overpass times in 2006 (Originally provided
! in HDF on variable pressure grid, changed to ASCII format and interpolated in log(p) to constant grid
! on GEOS-4 reduced levels)
SUBROUTINE GET_NO2(year, month, lon, lat, zs, ps, sza, no2, nz)  

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: deg2rad, du2mol
  USE OMSAO_variables_module,  ONLY: atmdbdir 
  USE ozprof_data_module,      ONLY: atmos_unit
  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)                          :: year, month, nz 
  REAL (KIND=dp), INTENT(IN)                   :: lat, lon
  REAL (KIND=dp), INTENT(IN)                   :: sza
  REAL (KIND=dp), DIMENSION(0:nz), INTENT(IN)  :: zs, ps
  REAL (KIND=dp), DIMENSION(nz),   INTENT(OUT) :: no2     ! In Dobson Units

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER          :: nlat = 18, nalt = 25, maxsza = 17
  INTEGER, PARAMETER          :: nglat=91, nglon=144, ngalt=19
  REAL (KIND=dp), PARAMETER   :: longrid = 2.5, latgrid = 2.0, lon0=-181.25, lat0=-91.0
  INTEGER                     :: errstat, i, j, k, nsza, fidx, lidx, nband, nblat, nblon, ntp
  INTEGER, DIMENSION(2)       :: latin,   lonin
  REAL(KIND=dp), DIMENSION(2) :: latfrac, lonfrac

  INTEGER,        SAVE, DIMENSION(nlat)                :: nszas
  REAL (KIND=dp), SAVE, DIMENSION(nlat, nalt, maxsza)  :: allno2
  REAL (KIND=dp), SAVE, DIMENSION(nlat, nalt)          :: alts
  REAL (KIND=dp), SAVE, DIMENSION(nlat, maxsza)        :: szas
  REAL (KIND=dp), SAVE, DIMENSION(nglon, nglat, ngalt) :: geosno2
  LOGICAL, SAVE                                        :: first = .TRUE.

  REAL (KIND=dp), DIMENSION(ngalt)            :: gprof
  REAL (KIND=dp), DIMENSION(0:ngalt)          :: geospres, geosalt
  REAL (KIND=dp), DIMENSION(maxsza)           :: tempszas
  REAL (KIND=dp), DIMENSION(nalt)             :: bprof, bpres
  REAL (KIND=dp), DIMENSION(0:nalt+ngalt)     :: cumno2, tempalt      
  REAL (KIND=dp), DIMENSION(0:nz)             :: temp
  REAL (KIND=dp)                              :: frac, csza, airdens

  ! Correct coordinates
  REAL (KIND=DP), DIMENSION(0:ngalt), PARAMETER:: pres = (/1.0d0,              &
       0.987871d0, 0.954730d0, 0.905120d0, 0.845000d0, 0.7800d0,   0.710000d0, &
       0.639000d0, 0.570000d0, 0.503000d0, 0.440000d0, 0.3800d0,   0.325000d0, &
       0.278000d0, 0.237954d0, 0.202593d0, 0.171495d0, 0.144267d0, 0.121347d0, &
       0.102098d0/)
  CHARACTER (LEN=3), DIMENSION(12)  :: months = (/'jan', 'feb','mar', 'apr',    &
       'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/)
  CHARACTER (LEN=130)               :: fname  
  CHARACTER (LEN=2)                 :: monc

  IF (first) THEN
     fname = TRIM(ADJUSTL(atmdbdir)) // 'NO2/' // months(month) // '_strat_NO2.dat' 
     OPEN(UNIT = atmos_unit, FILE = fname, status='old')
     
     DO i = 1, nlat 
        READ(atmos_unit, '(5X,2I5)') nszas(i), nszas(i)
        ! Only use AM data (from mid-night to noon)
        nszas(i) = nszas(i) / 2 
        READ(atmos_unit, '(5X, 34E13.4)') szas(i, 1:nszas(i)), szas(i, 1:nszas(i))
        DO j = nalt, 1, -1  ! Read top-down (j = 1, lowest level) 
           READ(atmos_unit, *) alts(i, j), allno2(i, j, 1:nszas(i)), &
                allno2(i, j, 1:nszas(i)), airdens
           allno2(i, j, 1:nszas(i)) = allno2(i, j, 1:nszas(i)) / airdens * 1.0E9 ! to ppbv           
        ENDDO
        READ(atmos_unit, *)
     ENDDO     
     CLOSE (atmos_unit)

     IF (year >= 2004) THEN
        WRITE(monc,  '(I2.2)') month  
        fname = TRIM(ADJUSTL(atmdbdir)) // 'NO2/gcno2_06' // monc // '.dat'
     ELSE
        fname = TRIM(ADJUSTL(atmdbdir)) // 'NO2/' // months(month) // '_no2_avg.dat'
     ENDIF
     OPEN(UNIT = atmos_unit, FILE = fname, status='old')
     DO i = 1, 5
        READ(atmos_unit, *) 
     ENDDO
     READ(atmos_unit, *) (((geosno2(i, j, k), k = 1, ngalt), j = 1, nglat), i = 1, nglon)
     CLOSE (atmos_unit)

     first = .FALSE.
  ENDIF 

  !  Get GEOS-CHEM profiles (ppbv)  
  CALL get_gridfrac(nglon, nglat, longrid, latgrid, lon0, lat0, &
       lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
  gprof = 0.0
  DO i = 1, nblon
     DO j = 1, nblat 
        gprof = gprof + geosno2(lonin(i), latin(j), :) * lonfrac(i) * latfrac(j)
     ENDDO
  ENDDO 

  geospres = pres * ps(0)  
  CALL BSPLINE(ps, zs, nz+1, geospres,  geosalt, ngalt+1, errstat)

  ! Interpolate over latitude (avoid discontinuity) 
  IF (lat <= -85.0) THEN
     nband = 1; latin(1) = 1; latfrac(1) = 1.0
  ELSE IF (lat >= 85.0) THEN
     nband = 1; latin(1) = nlat; latfrac(1) = 1.0
  ELSE
     nband = 2; frac = (lat + 85.0) / 10.0 + 1
     latin(1) = INT(frac); latin(2) = latin(1) + 1
     latfrac(1) = latin(2) - frac; latfrac(2) = 1.0 - latfrac(1)
  ENDIF  

  ! Get Profile in the stratosphere
  no2 = 0.0
  csza = COS(sza * deg2rad)
  DO i = 1, nband
     nsza = nszas(latin(i))
     tempszas(1:nsza) = COS(szas(latin(i), 1:nsza) * deg2rad)

     fidx = MINVAL(MAXLOC( tempszas(1:nsza), MASK=(tempszas(1:nsza) <= csza)))
     IF (fidx == nsza) THEN
        bprof = allno2(latin(i), :, nsza)
     ELSE IF (fidx == 0) THEN
        bprof = allno2(latin(i), :, 1)
     ELSE
        lidx = fidx + 1
        frac = 1.0 - (csza - tempszas(fidx)) / (tempszas(lidx) - tempszas(fidx))
        bprof = frac * allno2(latin(i), :, fidx) + (1 - frac) * allno2(latin(i), :, lidx)
     ENDIF

     ! Approximation: the layer-averaged profile is used as level (bottom) profile
     ntp = MINVAL(MAXLOC(geosalt(1:ngalt), MASK=(gprof > 0 .AND. geosalt(1:ngalt) < alts(latin(i), 1))))

     tempalt(0:ntp)          = geosalt(0:ntp)
     tempalt(ntp+1:ntp+nalt) = alts(latin(i), :)
     IF (tempalt(0) > zs(0)) tempalt(0) = zs(0)
     IF (tempalt(ntp+nalt) < zs(nz)) tempalt(ntp+nalt) = zs(nz)

     CALL BSPLINE(zs, ps, nz+1, alts(latin(i), :), bpres, nalt, errstat)
     cumno2 = 0.0
     DO j = 1, ntp
        cumno2(j) = cumno2(j-1) + gprof(j) * (geospres(j-1) - geospres(j)) / 1266.5625
     ENDDO
     cumno2(ntp + 1) = cumno2(ntp) + (gprof(ntp) + bprof(1)) * (geospres(ntp) - bpres(1)) / 2533.125
     DO j = 2, nalt
        k = ntp + j
        cumno2(k) = cumno2(k-1) + (bprof(j-1) + bprof(j)) * (bpres(j-1)-bpres(j)) / 2533.125
     ENDDO

     ! Interpolate to GOME retrieval altitudes
     CALL INTERPOL(tempalt(0:nalt+ntp), cumno2(0:nalt+ntp), nalt + ntp + 1, zs(0:nz), temp(0:nz), nz+1, errstat)
     temp(1:nz) = temp(1:nz) - temp(0:nz-1)
     no2 = no2 + temp(1:nz) * latfrac(i)     
  ENDDO

  no2(1:nz) = no2(1:nz) * du2mol
 
  RETURN  
END SUBROUTINE GET_NO2

! ps: bottom up
! read GEOS-4 fields (96-97, 99-00) and GEOS-3 fields (98) of SO2 from Randall Martin and Neil Moore
! use average for other years
! 2006 is provide by Chulkyu Lee on 30 altitude grids (reduced GEOS-4) in HDF5 format
! will be used for after 2004
SUBROUTINE GET_GEOSCHEM_SO2(year, month, lon, lat, ps, so2, ntp, nz)  
  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: du2mol
  USE OMSAO_variables_module, ONLY: atmdbdir 
  USE ozprof_data_module,     ONLY: atmos_unit

  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)                          :: month, nz, year, ntp
  REAL (KIND=dp), INTENT(IN)                   :: lon, lat
  REAL (KIND=dp), DIMENSION(0:nz), INTENT(IN)  :: ps
  REAL (KIND=dp), DIMENSION(nz),   INTENT(OUT) :: so2  ! in Dobson Units

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER        :: nlat = 91, nlon = 144, nalt_pre2006 = 21, nalt_post2006 = 30
  INTEGER, SAVE             :: nalt
  REAL (KIND=dp), PARAMETER :: longrid = 2.5, latgrid = 2.0, lon0=-181.25, lat0=-91.0
  INTEGER                   :: errstat, i, j, k, nblat, nblon, error
  REAL (KIND=dp), DIMENSION(0:nalt_post2006)          :: geospres, cumso2
  REAL (KIND=dp), DIMENSION(nalt_post2006)            :: gprof
  REAL (KIND=dp), DIMENSION(0:nz)                     :: tempso2
  LOGICAL, SAVE                                       :: file_exist
  INTEGER, DIMENSION(2)                               :: latin,   lonin
  REAL(KIND=dp), DIMENSION(2)                         :: latfrac, lonfrac

  REAL (KIND=dp), SAVE, DIMENSION(nlon,nlat,nalt_post2006) :: geosso2
  LOGICAL, SAVE                                       :: first = .TRUE.

  ! Correct coordinates (for geos3 fields)
  REAL (KIND=DP), DIMENSION(0:nalt_pre2006), PARAMETER:: pres3 = (/ &
       1.00000D+00, 9.97095D-01, 9.91200D-01, 9.81500D-01, 9.67100D-01, 9.46800D-01, &
       9.19500D-01, 8.84000D-01, 8.39000D-01, 7.83000D-01, 7.18200D-01, 6.47600D-01, &
       5.74100D-01, 5.00000D-01, 4.27800D-01, 3.59500D-01, 2.97050D-01, 2.41950D-01, &
       1.94640D-01, 1.55000D-01, 1.22680D-01, 9.69000D-02/)

  ! for geos4 fields
  REAL (KIND=DP), DIMENSION(0:nalt_pre2006), PARAMETER:: ap4    = (/  &
       0.000000D0,   0.000000D0,   12.704939D0,  35.465965D0,  66.098427D0,  101.671654D0,  &
       138.744400D0, 173.403183D0, 198.737839D0, 215.417526D0, 223.884689D0, 224.362869D0, &
       216.864929D0, 201.192093D0, 176.929993D0, 150.393005D0, 127.837006D0, 108.663429D0, &
       92.365662D0,  78.512299D0,  56.387939D0,  40.175419D0/)
  REAL (KIND=DP), DIMENSION(0:nalt_pre2006), PARAMETER:: bp4    = (/  &
       1.000000D0,   0.985110D0,   0.943290D0,   0.867830D0, 0.764920D0,  0.642710D0,  &
       0.510460D0,   0.378440D0,   0.270330D0,   0.183300D0, 0.115030D0,  0.063720D0,  &
       0.028010D0,   0.006960D0,   0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0, &
       0.000000D0,   0.000000D0,   0.000000D0,   0.000000D0/)
  !p(L) = AP4(L) + BP4(L) * ps

  ! for GEOS-4, 30 levels
!!$  REAL (KIND=DP), DIMENSION(nalt_post2006), PARAMETER:: etac    = (/  &
!!$       0.9926, 0.9707, 0.9300, 0.8680, 0.7891, 0.6987, 0.6031, 0.5135, 0.4373, 0.3724, &
!!$       0.3171, 0.2701, 0.2299, 0.1956, 0.1663, 0.1414, 0.1202, 0.1021, 0.0868, 0.0685, &
!!$       0.0491, 0.0348, 0.0245, 0.0148, 0.0068, 0.0029, 0.0011, 0.0004, 0.0001, 0.0000/)
  REAL (KIND=DP), DIMENSION(0:nalt_post2006), PARAMETER:: etae    = (/  &
    1.0000,  0.9851, 0.9562, 0.9039, 0.8321, 0.7460, 0.6515, 0.5547, 0.4723, 0.4022, &
    0.3425,  0.2917, 0.2484, 0.2114, 0.1798, 0.1528, 0.1299, 0.1104, 0.0939, 0.0798, &
    0.0573,  0.0408, 0.0288, 0.0201, 0.0094, 0.0041, 0.0017, 0.0006, 0.0002, 0.0001, 0/)



  CHARACTER (LEN=130)               :: geosfile
  CHARACTER (LEN=2)                 :: yearc, monc


  IF (first) THEN
     WRITE(yearc, '(I2.2)') MOD(year, 100)
     WRITE(monc,  '(I2.2)') month  

     IF (year >= 2004) THEN
        ! Currently only have 2006 from Chulkyu, we will use for all years since
        geosfile = TRIM(ADJUSTL(atmdbdir)) // 'gcso2/so2_06' // monc // '.dat'
     ELSE
        geosfile = TRIM(ADJUSTL(atmdbdir)) // 'gcso2/so2_' // yearc // monc // '.dat'
     ENDIF

     ! Determine if file exists or not
     INQUIRE (FILE= geosfile, EXIST= file_exist)
     IF (.NOT. file_exist) THEN
        WRITE(*, *) 'Warning: no SO2 profile file found, use 96-97, 99-00 average!!!'
        geosfile = TRIM(ADJUSTL(atmdbdir)) // 'gcso2/so2_avg' // monc // '.dat'
        nalt = nalt_pre2006
     ELSE 
        IF (year >= 2004) THEN 
           nalt = nalt_post2006
        ELSE
           nalt = nalt_pre2006
        ENDIF
     ENDIF

     OPEN(UNIT = atmos_unit, FILE = geosfile, status='old')
     DO i = 1, 5
        READ(atmos_unit, *) 
     ENDDO
     IF (year >= 2004 .AND. file_exist) THEN
        READ(atmos_unit, '(30E8.2)') (((geosso2(i, j, k), k = 1, nalt), j = 1, nlat), i = 1, nlon)
     ELSE
        READ(atmos_unit, '(21E8.2)') (((geosso2(i, j, k), k = 1, nalt), j = 1, nlat), i = 1, nlon)
     ENDIF
     CLOSE (atmos_unit)
     first = .FALSE.

  ENDIF

  CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
       lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)


  gprof = 0.0
  DO i = 1, nblon
     DO j = 1, nblat 
        gprof = gprof + geosso2(lonin(i), latin(j), :) * lonfrac(i) * latfrac(j)
     ENDDO
  ENDDO


  ! Need special processing for getting pressure profile for GEOS-4
  IF (year == 1998) THEN
     DO i = 0, nalt
        geospres(i) = pres3(i+1) * ps(0)
     ENDDO
  ELSEIF (year >= 2004 .AND. file_exist) THEN
     geospres = etae * ps(0)
  ELSE
     DO i = 0, nalt
        !geospres(i) = ap4(i+1) + bp4(i+1) * ps(0)
        geospres(i) = ap4(i) + bp4(i) * ps(0)
     ENDDO
  ENDIF


  ! Integrate from ppb to DU  
  cumso2 = 0.0
  DO i = 1, nalt  ! 1266.5625 = 1.25 * 1013.25
     cumso2(i) = cumso2(i-1) + gprof(i) * (geospres(i-1) - geospres(i)) / 1266.5625
  ENDDO
  !WRITE(90, *) nalt
  !WRITE(90, '(30D14.6)') (geospres(i), i=0, nalt)
  !WRITE(90, '(30D14.6)') (gprof(i), i=1, nalt)

  ! MAXLOC (maximum value), the index of maxloc starts from 1
  ! Note geos4 fields include some stratospheric part, which should not be used here
  !ntp = MINVAL(MINLOC(ps(0:nz), MASK = (ps(0:nz) >= geospres(nalt)))) - 1

  so2 = 0.0
  CALL INTERPOL(geospres, cumso2, nalt+1, ps(0:ntp), tempso2(0:ntp), ntp+1, errstat)
  so2(1:ntp) = tempso2(1:ntp) - tempso2(0:ntp-1)   ! DU at each layer 

  ! Assume 0.015 ppbv for stratospheric SO2
  ! 0.015 / 1.25 / 1013.25 = 1.1843E-5
  !DO i = ntp+1, nz
  !   so2(i) = (ps(i-1) - ps(i)) * 1.1843E-5
  !ENDDO


  so2(1:nz) = so2(1:nz) * du2mol


  RETURN  
END SUBROUTINE GET_GEOSCHEM_SO2

! =====================================================================
! Obtain AURA MLS zonal mean ozone profiles and its standard deviations
! (quality flags applied) 0.1-215 mb (i.e., 10-64 km), 36 latitude bins
! =====================================================================
SUBROUTINE get_mlso3prof(year, month, day, lat, nz, mnorstd, ps, zs, oz, ntp, errstat)
  USE OMSAO_precision_module 
  USE OMSAO_variables_module, ONLY: atmdbdir
  USE ozprof_data_module,     ONLY: atmos_unit, which_clima
  IMPLICIT NONE

  INTEGER, PARAMETER                           :: ml = 37, mlat=36
  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)                          :: year, month, day, nz, mnorstd
  INTEGER, INTENT(OUT)                         :: errstat, ntp
  REAL (KIND=dp),INTENT(IN)                    :: lat
  REAL (KIND=dp), DIMENSION(0:nz), INTENT(IN)  :: ps, zs
  REAL (KIND=dp), DIMENSION(nz), INTENT(INOUT) :: oz
   
  ! ======================
  ! Local variables
  ! ======================
  CHARACTER (LEN=130)              :: mlsfname
  CHARACTER (LEN=2)                :: monc, dayc
  CHARACTER (LEN=4)                :: yrc
  LOGICAL                          :: file_exist
  INTEGER                          :: i, j, ib, nband, fidx, lidx, ntmpl, sl, el
  REAL (KIND=dp), DIMENSION (ml)   :: tmpoz, tmpozstd, ratio
  REAL (KIND=dp), DIMENSION (0:ml) :: cumoz
  REAL (KIND=dp), DIMENSION (0:nz) :: tmpcumoz, tmps
  REAL (KIND=dp), DIMENSION (2)    :: latfrac
  REAL (KIND=dp)                   :: sumfrac
  INTEGER,        DIMENSION (2)    :: latin

  ! Saved variables
  INTEGER, SAVE                             :: nlat, nl
  REAL (KIND=dp), SAVE, DIMENSION(mlat, ml) :: mlsprofs, mlstds
  REAL (KIND=dp), SAVE, DIMENSION(0:ml)     :: mlsps
  REAL (KIND=dp), SAVE, DIMENSION(mlat)     :: mlslats
  LOGICAL,        SAVE                      :: first = .TRUE.

  errstat = 0
  IF (first) THEN
     WRITE(monc, '(I2.2)') month          ! from 9 to '09' 
     WRITE(dayc, '(I2.2)') day            ! from 9 to '09'     
     WRITE(yrc,  '(I4.4)') year           ! from 1997 to '1997'
     
     ! Check the availablity of MLS ozone profiles
     mlsfname =TRIM(ADJUSTL(atmdbdir)) // 'MLSO3/zm_v02_' // yrc // monc // dayc // '.dat'
    
     ! Determine if file exists or not
     INQUIRE (FILE= mlsfname, EXIST= file_exist)
     IF (.NOT. file_exist) THEN
        WRITE(*, *) 'No MLS ozone profile found!!!'; errstat = -1; RETURN
     ENDIF
     
     OPEN (UNIT = atmos_unit, file = mlsfname, status = 'unknown')
     READ (atmos_unit, *) nl
     !nl = nl-1      ! Only use above 147 mb
     IF (nl > ml) THEN
        WRITE(*, *) 'Need to increase ml from ', ml, ' to ', nl
        errstat = -1; RETURN
     ENDIF
     ! Reading pressure bottom up
     READ (atmos_unit, *); READ (atmos_unit, *) (mlsps(i), i = nl, 0, -1)
     mlsps(0:nl) = LOG(mlsps(0:nl)) 
     
     READ (atmos_unit, *) nlat
     IF (nlat > mlat) THEN
        WRITE(*, *) 'Need to increase mlat from ', mlat, ' to ', nlat
        errstat = -1; RETURN
     ENDIF
     READ (atmos_unit, *); READ (atmos_unit, *) mlslats(1:nlat)

     ! Reading bottom up
     READ (atmos_unit, *)
     DO i = 1, nlat
        READ (atmos_unit, *) (mlsprofs(i, j), j = nl, 1, -1)
     ENDDO
     READ (atmos_unit, *)
     DO i = 1, nlat
        READ (atmos_unit, *) (mlstds(i, j), j = nl, 1, -1)
     ENDDO
     CLOSE (atmos_unit)
    
     first = .FALSE.
  ENDIF

  IF (lat <= mlslats(1)) THEN
     nband = 1; latin(1) = 1; latfrac(1) = 1.0
  ELSE IF (lat >= mlslats(nlat)) THEN
     nband = 1; latin(1) = nlat; latfrac(1) = 1.0
  ELSE
     nband = 2 
     DO i = 2, nlat
        IF ( lat <= mlslats(i) ) THEN
           latin(1) = i - 1; latin(2) = i
           latfrac(2) = (lat - mlslats(i - 1)) / (mlslats(i) - mlslats(i-1))
           latfrac(1) = 1.0d0 - latfrac(2)
           EXIT
        ENDIF
     ENDDO
  ENDIF

  tmpoz(1:nl) = 0.0; tmpozstd(1:nl) = 0.0
  DO ib = 1, nband
     tmpoz(1:nl) = tmpoz(1:nl) + mlsprofs(latin(ib), 1:nl) * latfrac(ib)
     tmpozstd(1:nl) = tmpozstd(1:nl) + mlstds(latin(ib), 1:nl) * latfrac(ib)
  ENDDO
  ratio(1:nl) = tmpozstd(1:nl) / tmpoz(1:nl) * 100.0
 
  ! Only use MLS altitude range where reltative variability is < 50%
  ! Find first MLS layer to be used
  DO i = 1, nl
     IF (ratio(i) <= 50.0) THEN
        sl = i; EXIT
     ENDIF
  ENDDO

  ! Find last MLS layer to be used
  DO i = nl, 1, -1
     IF (ratio(i) <= 50.0) THEN
        el = i; EXIT
     ENDIF
  ENDDO
  
  IF (mnorstd == 2) tmpoz(1:nl) = tmpozstd(1:nz)
   
  ! Get cumulative ozone profile from (215 mb to 0.1 mb)
  cumoz(0) = 0.0
  DO i = 1, nl
     cumoz(i) = cumoz(i-1) + tmpoz(i)
  ENDDO
  tmps = LOG(ps(0:nz))

  fidx = MINVAL(MAXLOC(tmps(0:nz), MASK = (tmps(0:nz) <= mlsps(sl-1))) - 1)
  lidx = MINVAL(MINLOC(tmps(0:nz), MASK = (tmps(0:nz) >= mlsps(el))) - 1)
  ntmpl = lidx - fidx + 1
  !print *, sl, el, fidx, lidx, ntmpl
  !print *, EXP(mlsps(sl-1)), EXP(mlsps(el))
  !print *, EXP(tmps(fidx)), EXP(tmps(lidx))

  !print *, ' ozone before: ', SUM(oz)
  !print *, oz
  CALL BSPLINE(mlsps(0:nl), cumoz(0:nl), nl+1, tmps(fidx:lidx), &
       tmpcumoz(0:ntmpl-1), ntmpl, errstat)
  IF (errstat < 0) THEN
     WRITE(*, *) 'GET_MLSO3PROF: BSPLINE error, errstat = ', errstat; RETURN
  ENDIF
  oz(fidx+1:lidx) = tmpcumoz(1:ntmpl-1) - tmpcumoz(0:ntmpl-2)
  !print *, ' ozone after: ', SUM(oz)
  !print *, fidx+1, lidx, SUM(oz(fidx+1:lidx))
  !print *, oz

  ntp = fidx
  
  RETURN
END SUBROUTINE get_mlso3prof

! =====================================================================
! Obtain AURA MLS zonal mean ozone profiles and its standard deviations
! (quality flags applied) 0.1-215 mb (i.e., 10-64 km), 36 latitude bins
! =====================================================================
SUBROUTINE get_mlso3prof_single(year, month, day, lat, nz, mnorstd, ps, zs, oz, ntp, errstat)
  USE OMSAO_precision_module 
  USE OMSAO_variables_module, ONLY: atmdbdir
  USE ozprof_data_module,     ONLY: atmos_unit, which_clima
  IMPLICIT NONE

  INTEGER, PARAMETER                           :: ml = 37
  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)                          :: year, month, day, nz, mnorstd
  INTEGER, INTENT(OUT)                         :: errstat, ntp
  REAL (KIND=dp),INTENT(IN)                    :: lat
  REAL (KIND=dp), DIMENSION(0:nz), INTENT(IN)  :: ps, zs
  REAL (KIND=dp), DIMENSION(nz), INTENT(INOUT) :: oz
   
  ! ======================
  ! Local variables
  ! ======================
  CHARACTER (LEN=130)              :: mlsfname
  CHARACTER (LEN=2)                :: monc, dayc
  CHARACTER (LEN=4)                :: yrc
  LOGICAL                          :: file_exist
  INTEGER                          :: i, j, fidx, lidx, ntmpl, sl, el, nl, theprof, ios, nm
  REAL (KIND=dp), DIMENSION (ml)   :: tmpoz, ratio, mlsprof, mlstd
  REAL (KIND=dp), DIMENSION (0:ml) :: cumoz, mlsps
  REAL (KIND=dp), DIMENSION (0:nz) :: tmpcumoz, tmps
  REAL (KIND=dp)                   :: tmplon, tmplat, tmpsza, tmptime

  errstat = 0

  WRITE(monc, '(I2.2)') month          ! from 9 to '09' 
  WRITE(dayc, '(I2.2)') day            ! from 9 to '09'     
  WRITE(yrc,  '(I4.4)') year           ! from 1997 to '1997'
  
  ! Check the availablity of MLS ozone profiles
  mlsfname =TRIM(ADJUSTL(atmdbdir)) // 'MLSO3/mlso3_v02_' // yrc // monc // dayc // '.dat'
    
  ! Determine if file exists or not
  INQUIRE (FILE= mlsfname, EXIST= file_exist)
  IF (.NOT. file_exist) THEN
     WRITE(*, *) 'No MLS ozone profile found!!!'; errstat = -1; RETURN
  ENDIF

  OPEN (UNIT = atmos_unit, file = '../run/conf/INP/mlsprof_index.inp', status = 'unknown', IOSTAT=ios)
  IF (ios /= 0) THEN
     WRITE(*, *) 'Do not know which profile to choose!!!'; errstat = -1; RETURN
  ELSE
     READ (atmos_unit, *) theprof; CLOSE (atmos_unit)
  ENDIF
   
  OPEN (UNIT = atmos_unit, file = mlsfname, status = 'unknown')
  READ (atmos_unit, *) nm, nl
  IF (nl > ml) THEN
     WRITE(*, *) 'Need to increase ml from ', ml, ' to ', nl
     errstat = -1; CLOSE(atmos_unit); RETURN
  ENDIF
  IF (theprof > nm - 1) THEN
     WRITE(*, *) 'Do not have this profile!!!'
     errstat = -1; CLOSE(atmos_unit); RETURN
  ENDIF

  ! Reading pressure bottom up
  READ (atmos_unit, *) (mlsps(i), i = nl, 0, -1)
  mlsps(0:nl) = LOG(mlsps(0:nl)) 

  ! Skip profiles until the one we want
  DO i = 1, theprof
     READ (atmos_unit, *); READ (atmos_unit, *); READ (atmos_unit, *)
  ENDDO

  READ (atmos_unit, *) i, tmplon, tmplat, tmpsza, tmptime
  !WRITE(*, '(4F10.4)') tmplon, tmplat, tmpsza, tmptime

  ! Reading bottom up
  READ (atmos_unit, *) (mlsprof(j), j = nl, 1, -1)
  READ (atmos_unit, *) (mlstd(j),   j = nl, 1, -1)
  CLOSE (atmos_unit)
  ratio(1:nl) = mlstd(1:nl) / mlsprof(1:nl) * 100.0
 
  ! Only use MLS altitude range where reltative variability is < 50%
  ! Find first MLS layer to be used
  DO i = 1, nl
     IF (ratio(i) <= 50.0) THEN
        sl = i; EXIT
     ENDIF
  ENDDO

  ! Find last MLS layer to be used
  DO i = nl, 1, -1
     IF (ratio(i) <= 50.0) THEN
        el = i; EXIT
     ENDIF
  ENDDO
  
  IF (mnorstd == 1) THEN
     tmpoz(1:nl) = mlsprof(1:nz)
  ELSE IF (mnorstd == 2) THEN
     tmpoz(1:nl) = mlstd(1:nz)
  ENDIF
       
  ! Get cumulative ozone profile from 
  cumoz(0) = 0.0
  DO i = 1, nl
     cumoz(i) = cumoz(i-1) + tmpoz(i)
  ENDDO
  tmps = LOG(ps(0:nz))

  fidx = MINVAL(MAXLOC(tmps(0:nz), MASK = (tmps(0:nz) <= mlsps(sl-1))) - 1)
  lidx = MINVAL(MINLOC(tmps(0:nz), MASK = (tmps(0:nz) >= mlsps(el))) - 1)
  ntmpl = lidx - fidx + 1
  !print *, sl, el, fidx, lidx, ntmpl
  !print *, EXP(mlsps(sl-1)), EXP(mlsps(el))
  !print *, EXP(tmps(fidx)), EXP(tmps(lidx))

  !print *, ' ozone before: ', SUM(oz)
  !print *, oz
  CALL BSPLINE(mlsps(0:nl), cumoz(0:nl), nl+1, tmps(fidx:lidx), &
       tmpcumoz(0:ntmpl-1), ntmpl, errstat)
  IF (errstat < 0) THEN
     WRITE(*, *) 'GET_MLSO3PROF: BSPLINE error, errstat = ', errstat; RETURN
  ENDIF
  oz(fidx+1:lidx) = tmpcumoz(1:ntmpl-1) - tmpcumoz(0:ntmpl-2)
  !print *, ' ozone after: ', SUM(oz)
  !print *, fidx+1, lidx, SUM(oz(fidx+1:lidx))
  !print *, oz

  ntp = fidx
  
  RETURN
END SUBROUTINE get_mlso3prof_single

SUBROUTINE GET_NORMTOZ(year, month, day, lat, toz, nz, ntp, ps, zs, oz, errstat)

 USE OMSAO_precision_module 
 USE OMSAO_variables_module, ONLY: atmdbdir
 USE ozprof_data_module,     ONLY: atmos_unit, norm_tropo3, which_toz
 IMPLICIT NONE

  INTEGER, PARAMETER                           :: ntlat = 180

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)                          :: year, month, day, nz, ntp
  INTEGER, INTENT(OUT)                         :: errstat
  REAL (KIND=dp),INTENT(IN)                    :: lat
  REAL (KIND=dp),INTENT(INOUT)                 :: toz
  REAL (KIND=dp), DIMENSION(0:nz), INTENT(IN)  :: ps, zs
  REAL (KIND=dp), DIMENSION(nz), INTENT(INOUT) :: oz
   
  ! ======================
  ! Local variables
  ! ======================
  CHARACTER (LEN=130)              :: omto3fname
  CHARACTER (LEN=2)                :: monc, dayc
  CHARACTER (LEN=4)                :: yrc
  LOGICAL                          :: file_exist
  INTEGER                          :: i, j, ib, nband
  REAL (KIND=dp), DIMENSION (2)    :: latfrac
  REAL (KIND=dp)                   :: sumfrac, mnalt, do3
  INTEGER,        DIMENSION (2)    :: latin

  ! Saved variables
  REAL (KIND=dp), SAVE, DIMENSION(ntlat) :: zmto3, tlats, zmalt
  LOGICAL,        SAVE                   :: first = .TRUE.

  errstat = 0
  IF (which_toz == 2) THEN 
     IF (first) THEN
        WRITE(monc, '(I2.2)') month          ! from 9 to '09' 
        WRITE(dayc, '(I2.2)') day            ! from 9 to '09'     
        WRITE(yrc,  '(I4.4)') year           ! from 1997 to '1997'
        
        ! Check the availablity of MLS ozone profiles
        omto3fname =TRIM(ADJUSTL(atmdbdir)) // 'OMTO3/zm_v003_' // yrc // 'm' // monc // dayc // '.dat'
        
        ! Determine if file exists or not
        INQUIRE (FILE= omto3fname, EXIST= file_exist)
        IF (.NOT. file_exist) THEN
           WRITE(*, *) 'No Zonal Mean OMTO3 found!!!'; errstat = -1; RETURN
        ENDIF
        OPEN (UNIT = atmos_unit, file = omto3fname, status = 'unknown')
        DO i = 1, ntlat
           tlats(i) = REAL(i, KIND=dp) - 89.5
        ENDDO
        READ (atmos_unit, *)
        READ (atmos_unit, *) ((zmto3(i), zmalt(i)), i = 1, ntlat)
        CLOSE(atmos_unit)
        
        first = .FALSE.
     ENDIF
     
     IF (lat <= -89.5) THEN
        nband = 1; latin(1) = 1; latfrac(1) = 1.0
     ELSE IF (lat >= 89.5) THEN
        nband = 1; latin(1) = ntlat; latfrac(1) = 1.0
     ELSE
        nband = 2 
        DO i = 2, ntlat
           IF ( lat <= tlats(i)) THEN
              latin(1) = i - 1; latin(2) = i
              latfrac(2) = (lat - tlats(i - 1)) / (tlats(i) - tlats(i-1))
              latfrac(1) = 1.0d0 - latfrac(2)
              EXIT
           ENDIF
        ENDDO
     ENDIF
     
     toz = 0.0; sumfrac = 0.0; mnalt = 0.0
     DO ib = 1, nband
        IF ( zmto3(latin(ib)) > 0.0 ) THEN
           toz = toz + zmto3(latin(ib)) * latfrac(ib)
           mnalt = mnalt + zmalt(latin(ib)) * latfrac(ib)
           sumfrac = sumfrac + latfrac(ib)
        ENDIF
     ENDDO
     toz  = toz / sumfrac
     mnalt = mnalt / sumfrac / 1000.0 
  
     ! Accounting for different terrain height using approximate pressure conversion
     do3 = ( 1013.25 * (10.0**(-mnalt / 16.0)) - ps(0) ) / (ps(0) - ps(1)) * oz(1)
     toz = toz - do3
  ENDIF

  IF (norm_tropo3) THEN
     oz(1:ntp) = oz(1:ntp) * (toz - SUM(oz(ntp+1:nz))) / SUM(oz(1:ntp))
  ELSE
     oz(1:nz) = oz(1:nz) * toz / SUM(oz(1:nz))
  ENDIF
     
  RETURN
END SUBROUTINE GET_NORMTOZ
  

! ===============================================================
! Obtain TOMS V8 ozone profiles (12 month, 18 latitude bands,
!   3-10 profiles with total ozone at a step of 50 DU
! ===============================================================
SUBROUTINE get_tomsv8_clima(month, day, lat, toz, nl, ps, apoz, oz, errstat)

  USE OMSAO_parameters_module, ONLY: p0
  USE OMSAO_precision_module 
  USE OMSAO_variables_module,  ONLY: atmdbdir
  USE ozprof_data_module,      ONLY: atmos_unit
  USE OMSAO_errstat_module
  IMPLICIT NONE

  INTEGER, PARAMETER                          :: nl0 = 11
  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)                          :: month, day, nl
  INTEGER, INTENT(OUT)                         :: errstat
  REAL (KIND=dp),INTENT(IN)                    :: lat
  REAL (KIND=dp), INTENT(IN)                   :: toz
  REAL (KIND=dp), DIMENSION(0:nl), INTENT(IN)  :: ps
  REAL (KIND=dp), DIMENSION(nl), INTENT(IN)    :: apoz
  REAL (KIND=dp), DIMENSION(nl), INTENT(OUT)   :: oz
  
  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER :: nlat=18, maxprof=10, nmon=12
  CHARACTER (LEN=3), DIMENSION(12)  :: months = (/'jan', 'feb','mar', 'apr', &
       'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/)
  CHARACTER (LEN=130)                                :: ozprof_fname
  CHARACTER (LEN=200)                                :: line

  ! saved variables
  REAL (KIND=dp), SAVE, DIMENSION(nmon, nlat, maxprof, nl0) :: ozprofs
  INTEGER,        SAVE, DIMENSION(nmon, nlat)               :: nprofs
  REAL (KIND=dp), SAVE, DIMENSION(0:nl0)                    :: pv80
  LOGICAL,        SAVE                                      :: first = .TRUE.

  REAL (KIND=dp)                                            :: frac, tmp, fdum, maxoz, minoz
  REAL (KIND=dp), DIMENSION(nl0)                            :: oz0
  REAL (KIND=dp), DIMENSION(0:nl0)                          :: cum0,pv80g
  REAL (KIND=dp), DIMENSION(0:nl)                           :: logps, cum
  REAL (KIND=dp), DIMENSION(2)                              :: latfrac, monfrac
  INTEGER,        DIMENSION(2)                              :: latin, monin
  INTEGER :: i, j, ib, profin, nprof, nband, nm, im

  CHARACTER (LEN=16), PARAMETER :: modulename = 'get_tomsv8_clima'

  IF (first) THEN
     ! read the TOMS V8 profiles
     ozprof_fname = TRIM(ADJUSTL(atmdbdir)) // 'v8clima/tomsv8_ozone_clima.dat'
     OPEN (UNIT = atmos_unit, file= ozprof_fname, status = 'unknown')
     
     ! Read until the target month        
     DO im = 1, nmon
        DO i = 1, nlat 
           READ(atmos_unit, *) 
           nprof = 1
           DO j = 1, maxprof
              READ (atmos_unit, '(A)') line;  READ (line, *) fdum
              
              IF (fdum < 999.0) THEN
                 READ (line, *) fdum, ozprofs(im, i, nprof, :)
                 nprof = nprof + 1
              ENDIF
           ENDDO
           nprofs(im, i) = nprof - 1              
        ENDDO
     ENDDO
     CLOSE (atmos_unit)

     pv80(0) = ps(0)
     DO i = 1, nl0
        pv80(i) = p0 * 2.0D0 ** (+ i - nl0)
     ENDDO
     
     
     first = .FALSE.
  ENDIF

  IF (day <= 15) THEN
     monin(1) = month - 1
     IF (monin(1) == 0) monin(1) = 12
     monin(2) = month
     monfrac(1) = (15.0 - day) / 30.0
     monfrac(2) = 1.0 - monfrac(1)
  ELSE 
     monin(2) = month + 1
     IF (monin(2) == 13) monin(2) = 1
     monin(1) = month
     monfrac(2) = (day - 15) / 30.0
     monfrac(1) = 1.0 - monfrac(2)
  ENDIF
  nm = 2

  IF (lat <= -85.0) THEN
     nband = 1; latin(1) = 1; latfrac(1) = 1.0
  ELSE IF (lat >= 85.0) THEN
     nband = 1; latin(1) = nlat; latfrac(1) = 1.0
  ELSE
     nband = 2     ; frac = (lat + 85.0) / 10.0 + 1
     latin(1) = INT(frac); latin(2) = latin(1) + 1
     latfrac(1) = latin(2) - frac; latfrac(2) = 1.0 - latfrac(1)
  ENDIF

  oz0 = 0.0
  DO im = 1, nm
     DO ib = 1, nband   
        nprof = nprofs(monin(im), latin(ib))
        minoz = SUM(ozprofs(monin(im), latin(ib), 1, :))
        maxoz = SUM(ozprofs(monin(im), latin(ib), nprof, :))
        
        IF (toz < minoz) THEN
           !WRITE(*,*), 'Warning: no a priori profile available!!!'
           oz0  = oz0 + ozprofs(monin(im), latin(ib), 1, :) * toz / minoz * latfrac(ib)
        ELSE IF (toz > maxoz) THEN
           !WRITE(*,*), 'Warning: no a priori profile available!!!'
           oz0 = oz0 + ozprofs(monin(im), latin(ib), nprof, :) * toz / maxoz * latfrac(ib)
        ELSE
           profin = INT ((toz - minoz ) / 50.0) + 1
           IF (profin == 0) THEN 
              profin = 1
           ELSE IF (profin == nprof) THEN
              profin = profin - 1
           ENDIF
           
           frac = 1.0 - (toz - (minoz + (profin-1) * 50.0)) / 50.0
           oz0 = oz0 + latfrac(ib) * monfrac(im) * (frac * ozprofs(monin(im), latin(ib), profin, :) &
                + (1.0 - frac) * ozprofs(monin(im), latin(ib), profin+1, :))
        ENDIF
     ENDDO
  ENDDO
  CALL REVERSE(oz0(1:nl0), nl0)

    

  ! Bondary layer correction  
  IF (ps(nl) > p0 ) then !sfc
      tmp = ( ps(nl) - p0)/(pv80(nl0)-pv80(nl0-1))
      pv80(nl0) = ps(nl)
  ENDIF
  IF ( ps(0) < pv80(0) ) pv80(0) = ps(0)

  ! Interpolate ozone profile to the input pressure grid
  cum0(0) = 0.0
  DO i = 1, nl0
     cum0(i) = cum0(i-1) + oz0(i)
  ENDDO
  
  pv80g(0:nl0) = LOG(pv80(0:nl0)) 
  logps        = LOG(ps)
     
  errstat = pge_errstat_ok
  
  CALL BSPLINE(pv80g, cum0, nl0+1, logps, cum, nl+1, errstat)
  IF (errstat < 0) THEN
     WRITE(*, *) modulename, ': INTERPOL error, errstat = ', errstat
     errstat = pge_errstat_error; RETURN
  ENDIF
  oz = cum(1:nl) - cum(0:nl-1)

  ! Correct for top few layers using original xap (based on McPeters Clima)
  DO i = 0, nl
     IF (logps(i) >= pv80g(1)) EXIT
  ENDDO
  oz(1:i) = apoz(1:i) * SUM(oz(1:i)) / SUM(apoz(1:i))
  
  RETURN
END SUBROUTINE get_tomsv8_clima

