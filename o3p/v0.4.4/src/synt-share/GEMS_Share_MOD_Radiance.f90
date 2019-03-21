!-------------------------------------------------------------------------------
!+Module to calculate Radinace

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE Share_MOD_Radiance

!-------------------------------------------------------------------------------
!+Description: 
!     Randiance Calculation Module
!
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2014.12.01 Fisrt Code      (BR Kim    , EWHA WOMAN UNIV.) 
! 0.2     2014.12.30 Standardization (YuGeun Ki , SAEASoft) 
!-------------------------------------------------------------------------------

! Modules used:


!**********************************************************!
IMPLICIT NONE

INTEGER(KIND=4)                     :: nwavel_rad    
INTEGER(KIND=4)                     :: nxtrack_rad   
INTEGER(KIND=4)                     :: ncoef_rad 
INTEGER(KIND=4)                     :: ntimes_rad      
INTEGER, POINTER                    :: l_missing_rad(:, :)

INTERFACE GEMS_Share_Calc_ObservationWavelen
    MODULE PROCEDURE GEMS_Share_Calc_ObservationWavelen_2d
    MODULE PROCEDURE GEMS_Share_Calc_ObservationWavelen_3d
END INTERFACE

INTERFACE GEMS_Share_Calc_RadianceValue
    MODULE PROCEDURE GEMS_Share_Calc_RadianceValue_0d       ! Output Data type REAL*8
    MODULE PROCEDURE GEMS_Share_Calc_RadianceValue_1d       ! Output Data type REAL*8
    MODULE PROCEDURE GEMS_Share_Calc_RadianceValue_2d       ! Output Data type REAL*8       
    MODULE PROCEDURE GEMS_Share_Calc_RadianceValue_3d       ! Output Data type REAL*8
    MODULE PROCEDURE GEMS_Share_Calc_RadianceValue_0d_2     ! Output Data type REAL*4
    MODULE PROCEDURE GEMS_Share_Calc_RadianceValue_1d_2     ! Output Data type REAL*4
    MODULE PROCEDURE GEMS_Share_Calc_RadianceValue_2d_2     ! Output Data type REAL*4
    MODULE PROCEDURE GEMS_Share_Calc_RadianceValue_3d_2     ! Output Data type REAL*4
END INTERFACE

INTERFACE GEMS_Share_init4RadianceCalc
    MODULE PROCEDURE GEMS_Share_init4RadianceCalci16
    MODULE PROCEDURE GEMS_Share_init4RadianceCalci32
END INTERFACE

CONTAINS

!-------------------------------------------------------------------------------
!+Description: 
!        2차원 입력값으로부터 2차원 관측파장을 계산
!
!
! Method: 
!        
!
! Input files:
!       refcol : 
!       coef   : 
!
! Output files:
!       wave   : 
!
!-------------------------------------------------------------------------------
SUBROUTINE GEMS_Share_Calc_ObservationWavelen_2d(refcol, coef, wave)
    IMPLICIT NONE

    INTEGER*2,            INTENT(IN)    :: refcol
    REAL, DIMENSION(:,:), INTENT(IN)    :: coef
    REAL, DIMENSION(:,:), INTENT(OUT)   :: wave

    INTEGER     :: i, j     !, nxtrack_rad, nwavel_rad
    REAL        :: x, w
    REAL, ALLOCATABLE, DIMENSION(:)     :: c

    ALLOCATE (c(ncoef_rad))

    DO i = 1, nwavel_rad
        DO j = 1, nxtrack_rad
            c = coef(:,j)
            x = (i-1)-refcol
            w = c(1)+x*(c(2)+x*(c(3)+x*(c(4)+x*c(5))))

            wave(i,j) = w
        ENDDO
    ENDDO

    DEALLOCATE (c)

END SUBROUTINE GEMS_Share_Calc_ObservationWavelen_2d

!-------------------------------------------------------------------------------
!+Description: 
!        3차원 입력값으로부터 3차원 관측파장을 계산
!
!
! Method:
!
!
! Input files:
!       refcol : 
!       coef   : 
!
! Output files:
!       wave   : 
!
!-------------------------------------------------------------------------------
SUBROUTINE GEMS_Share_Calc_ObservationWavelen_3d(refcol, coef, wave)
    IMPLICIT NONE

    INTEGER*2, DIMENSION(:),     INTENT(IN)     :: refcol
    REAL,      DIMENSION(:,:,:), INTENT(IN)     :: coef
    REAL,      DIMENSION(:,:,:), INTENT(OUT)    :: wave

    INTEGER     :: i, j, k    !, ntimes_rad, nxtrack_rad, nwavel_rad !loop parameter
    REAL        :: x, w
    REAL, ALLOCATABLE, DIMENSION(:) :: c

    ALLOCATE (c(ncoef_rad))

! the slow way to calculate wavelength -- brute force
    DO i = 1, nwavel_rad
        DO j = 1, nxtrack_rad
            DO k = 1, ntimes_rad
                c = coef(:,j,k)
                x = (i-1)-refcol(k)
                w = c(1)+x*(c(2)+x*(c(3)+x*(c(4)+x*c(5))))

                wave(i,j,k) = w
            ENDDO
        ENDDO
    ENDDO

    DEALLOCATE (c)

END SUBROUTINE GEMS_Share_Calc_ObservationWavelen_3d

!-------------------------------------------------------------------------------
!+Description: 
!        scalar 입력값에서 Radiance계산
!
!
! Method:
!
!
! Input files:
!       radman : 
!       radexp : 
!
! Output files:
!       rad    : 
!
!-------------------------------------------------------------------------------
SUBROUTINE GEMS_Share_Calc_RadianceValue_0d(radman, radexp, rad)
! This subroutine calculate radiance using Radiance Mantissa and Exponent value.
    IMPLICIT NONE

    INTEGER*1,  INTENT(IN)  :: radexp
    INTEGER*2,  INTENT(IN)  :: radman
    REAL*8,     INTENT(OUT) :: rad

    rad = float(radman)*10.**float(radexp)
    IF (radexp .lt. 0   .or.  &
        radexp .gt. 128 .or.  &
        radman .lt. 0.          ) THEN

        rad = -999
        !l_missing_rad(:,:)=1

        !print*,'warning'
    ENDIF

END SUBROUTINE GEMS_Share_Calc_RadianceValue_0d

!-------------------------------------------------------------------------------
!+Description: 
!        1차원 배열의 입력값에서 Radiance계산
!
!
! Method:
!
!
! Input files:
!       radman : 
!       radexp : 
!
! Output files:
!       rad    : 
!
!-------------------------------------------------------------------------------
SUBROUTINE GEMS_Share_Calc_RadianceValue_1d(radman, radexp, rad)
    IMPLICIT NONE

    INTEGER*1, DIMENSION(:), INTENT(IN)     :: radexp
    INTEGER*2, DIMENSION(:), INTENT(IN)     :: radman
    REAL*8,    DIMENSION(:), INTENT(OUT)    :: rad

    INTEGER :: nn, mm, ll
    INTEGER ::  i,  j,  k

    nn = size(radexp)

    DO i = 1, nn
        rad(i) = float(radman(i))*10.**float(radexp(i))
        IF (radexp(i) .lt. 0    .or.  &
            radexp(i) .gt. 128  .or.  &
            radman(i) .lt. 0.           ) THEN

            rad(i) = -999
            !l_missing_rad(:,:)=1

            !print*,'warning'
        ENDIF
    END DO

END SUBROUTINE GEMS_Share_Calc_RadianceValue_1d

!-------------------------------------------------------------------------------
!+Description: 
!        2차원 배열의 입력값에서 Radiance계산
!
!
! Method:
!
!
! Input files:
!       radman : 
!       radexp : 
!
! Output files:
!       rad    : 
!
!-------------------------------------------------------------------------------
SUBROUTINE GEMS_Share_Calc_RadianceValue_2d(radman, radexp, rad)
    IMPLICIT NONE

    INTEGER*1,  DIMENSION(:,:), INTENT(IN)  :: radexp
    INTEGER*2,  DIMENSION(:,:), INTENT(IN)  :: radman
    REAL*8,     DIMENSION(:,:), INTENT(OUT) :: rad

    INTEGER :: nn, mm, ll
    INTEGER ::  i,  j,  k

    nn = size(radexp,1)
    mm = size(radexp,2)

    DO j = 1, mm
        DO i = 1, nn
            rad(i,j) = float(radman(i,j))*10.**float(radexp(i,j))
            IF (radexp(i,j) .lt. 0   .or.  &
                radexp(i,j) .gt. 128 .or.  &
                radman(i,j) .lt. 0.           ) THEN

                rad(i,j) = 0.       !1.e+34
                l_missing_rad(j,:)=1
            ENDIF
        END DO
    END DO

END SUBROUTINE GEMS_Share_Calc_RadianceValue_2d

!-------------------------------------------------------------------------------
!+Description:
!        3차원 배열의 입력값에서 Radiance계산
!
!
! Method:
!
!
! Input files:
!       radman : 
!       radexp : 
!
! Output files:
!       rad    : 
!
!-------------------------------------------------------------------------------
SUBROUTINE GEMS_Share_Calc_RadianceValue_3d(radman, radexp, rad)
    IMPLICIT NONE

    INTEGER*1,  DIMENSION(:,:,:), INTENT(IN)    :: radexp
    INTEGER*2,  DIMENSION(:,:,:), INTENT(IN)    :: radman
    REAL*8,     DIMENSION(:,:,:), INTENT(OUT)   :: rad

    INTEGER :: nn, mm, ll
    INTEGER ::  i,  j,  k

    nn = size(radexp,1)
    mm = size(radexp,2)
    ll = size(radexp,3)

    DO k = 1, ll
        DO j = 1, mm
            DO i = 1, nn
                rad(i,j,k) = float(radman(i,j,k))*10.**float(radexp(i,j,k))
                IF (radexp(i,j,k) .lt. 0   .or.  &
                    radexp(i,j,k) .gt. 128 .or.  &
                    radman(i,j,k) .lt. 0.         ) THEN

                    rad(i,j,k) = 0.             !1.e+34
                    l_missing_rad(j,k)=1

                    !print*,'exp missing',radexp(i,j,k)
                ENDIF
            END DO
        END DO
    END DO

END SUBROUTINE GEMS_Share_Calc_RadianceValue_3d

!-------------------------------------------------------------------------------
!+Description: 
!        scalar 입력값에서 Radiance계산
!
!
! Method:
!
!
! Input files:
!       radman : 
!       radexp : 
!
! Output files:
!       rad    : 
!
!-------------------------------------------------------------------------------
SUBROUTINE GEMS_Share_Calc_RadianceValue_0d_2(radman, radexp, rad)
! This subroutine calculate radiance using Radiance Mantissa and Exponent value.
    IMPLICIT NONE

    INTEGER*1,  INTENT(IN)  :: radexp
    INTEGER*2,  INTENT(IN)  :: radman
    REAL*4,     INTENT(OUT) :: rad

    rad = float(radman)*10.**float(radexp)
    IF (radexp .lt. 0   .or.  &
        radexp .gt. 128 .or.  &
        radman .lt. 0.          ) THEN

        rad = -999
        !l_missing_rad(:,:)=1

        !print*,'warning'
    ENDIF

END SUBROUTINE GEMS_Share_Calc_RadianceValue_0d_2

!-------------------------------------------------------------------------------
!+Description: 
!        1차원 배열의 입력값에서 Radiance계산
!
!
! Method:
!
!
! Input files:
!       radman : 
!       radexp : 
!
! Output files:
!       rad    : 
!
!-------------------------------------------------------------------------------
SUBROUTINE GEMS_Share_Calc_RadianceValue_1d_2(radman, radexp, rad)
    IMPLICIT NONE

    INTEGER*1, DIMENSION(:), INTENT(IN)     :: radexp
    INTEGER*2, DIMENSION(:), INTENT(IN)     :: radman
    REAL*4,    DIMENSION(:), INTENT(OUT)    :: rad

    INTEGER :: nn, mm, ll
    INTEGER ::  i,  j,  k

    nn = size(radexp)

    DO i = 1, nn
        rad(i) = float(radman(i))*10.**float(radexp(i))
        IF (radexp(i) .lt. 0    .or.  &
            radexp(i) .gt. 128  .or.  &
            radman(i) .lt. 0.           ) THEN

            rad(i) = -999
            !l_missing_rad(:,:)=1

            !print*,'warning'
        ENDIF
    END DO

END SUBROUTINE GEMS_Share_Calc_RadianceValue_1d_2

!-------------------------------------------------------------------------------
!+Description: 
!        2차원 배열의 입력값에서 Radiance계산
!
!
! Method:
!
!
! Input files:
!       radman : 
!       radexp : 
!
! Output files:
!       rad    : 
!
!-------------------------------------------------------------------------------
SUBROUTINE GEMS_Share_Calc_RadianceValue_2d_2(radman, radexp, rad)
    IMPLICIT NONE

    INTEGER*1,  DIMENSION(:,:), INTENT(IN)  :: radexp
    INTEGER*2,  DIMENSION(:,:), INTENT(IN)  :: radman
    REAL*4,     DIMENSION(:,:), INTENT(OUT) :: rad

    INTEGER :: nn, mm, ll
    INTEGER ::  i,  j,  k

    nn = size(radexp,1)
    mm = size(radexp,2)

    DO j = 1, mm
        DO i = 1, nn
            rad(i,j) = float(radman(i,j))*10.**float(radexp(i,j))
            IF (radexp(i,j) .lt. 0   .or.  &
                radexp(i,j) .gt. 128 .or.  &
                radman(i,j) .lt. 0.           ) THEN

                rad(i,j) = 0.       !1.e+34
                l_missing_rad(j,:)=1
            ENDIF
        END DO
    END DO

END SUBROUTINE GEMS_Share_Calc_RadianceValue_2d_2

!-------------------------------------------------------------------------------
!+Description:
!        3차원 배열의 입력값에서 Radiance계산
!
!
! Method:
!
!
! Input files:
!       radman : 
!       radexp : 
!
! Output files:
!       rad    : 
!
!-------------------------------------------------------------------------------
SUBROUTINE GEMS_Share_Calc_RadianceValue_3d_2(radman, radexp, rad)
    IMPLICIT NONE

    INTEGER*1,  DIMENSION(:,:,:), INTENT(IN)    :: radexp
    INTEGER*2,  DIMENSION(:,:,:), INTENT(IN)    :: radman
    REAL*4,     DIMENSION(:,:,:), INTENT(OUT)   :: rad

    INTEGER :: nn, mm, ll
    INTEGER ::  i,  j,  k

    nn = size(radexp,1)
    mm = size(radexp,2)
    ll = size(radexp,3)

    DO k = 1, ll
        DO j = 1, mm
            DO i = 1, nn
                rad(i,j,k) = float(radman(i,j,k))*10.**float(radexp(i,j,k))
                IF (radexp(i,j,k) .lt. 0   .or.  &
                    radexp(i,j,k) .gt. 128 .or.  &
                    radman(i,j,k) .lt. 0.         ) THEN

                    rad(i,j,k) = 0.             !1.e+34
                    l_missing_rad(j,k)=1

                    !print*,'exp missing',radexp(i,j,k)
                ENDIF
            END DO
        END DO
    END DO

END SUBROUTINE GEMS_Share_Calc_RadianceValue_3d_2


!-------------------------------------------------------------------------------
!+Description: 
!      Radinace계산 서브루틴에서 사용된 전역상수값을 초기화
!
!
! Method:
!
!
! Input files:
!
! Output files:
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2015. .   Fisrt Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------
SUBROUTINE GEMS_Share_init4RadianceCalci16(nwavel, nxtrack, ntimes, ncoef)
    INTEGER(KIND=2),       INTENT(IN) :: nwavel
    INTEGER(KIND=2),       INTENT(IN) :: nxtrack
    INTEGER(KIND=2),       INTENT(IN) :: ntimes
    INTEGER(KIND=2),       INTENT(IN) :: ncoef
   
    nwavel_rad      =   nwavel
    nxtrack_rad     =   nxtrack
    ntimes_rad      =   ntimes
    ncoef_rad       =   ncoef

END SUBROUTINE GEMS_Share_init4RadianceCalci16

SUBROUTINE GEMS_Share_init4RadianceCalci32(nwavel, nxtrack, ntimes, ncoef)
    INTEGER(KIND=4),       INTENT(IN) :: nwavel
    INTEGER(KIND=4),       INTENT(IN) :: nxtrack
    INTEGER(KIND=4),       INTENT(IN) :: ntimes
    INTEGER(KIND=4),       INTENT(IN) :: ncoef
   
    nwavel_rad      =   nwavel
    nxtrack_rad     =   nxtrack
    ntimes_rad      =   ntimes
    ncoef_rad       =   ncoef

END SUBROUTINE GEMS_Share_init4RadianceCalci32

!-------------------------------------------------------------------------------
!+Description: 
!      Radinace계산 서브루틴에서 사용된 전역변수값을 MAIN에서 사용된 변수로 저장
!
!
! Method:
!
!
! Input files:
!
! Output files:
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2015. .   Fisrt Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------
SUBROUTINE saveGvtoLv4RadianceCalc(l_missing)
    INTEGER, INTENT(OUT)    :: l_missing(:,:)

    l_missing       =   l_missing_rad

END SUBROUTINE saveGvtoLv4RadianceCalc
  

!-------------------------------------------------------------------------------
!+Description: 
!      Radiance계산을 위한 전역변수 메모리를 할당해줌
!
!
! Method:
!
!
! Input files:
!
! Output files:
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2015.01.27 Fisrt Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------
SUBROUTINE GEMS_Share_allocMem4Rad

    ALLOCATE( l_missing_rad(nxtrack_rad, ntimes_rad) )

    l_missing_rad   = 0

END SUBROUTINE GEMS_Share_allocMem4Rad

!-------------------------------------------------------------------------------
!+Description: 
!      Radiance계산에 참여한 전역변수 메모리를 해제함
!
!
! Method:
!
!
! Input files:
!
! Output files:
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2015.01.27 Fisrt Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------
SUBROUTINE GEMS_Share_deAllocMem4Rad

    IF (ASSOCIATED(l_missing_rad))   DEALLOCATE(l_missing_rad )

END SUBROUTINE GEMS_Share_deAllocMem4Rad

END MODULE Share_MOD_Radiance

!-------------------------------------------------------------------------------
