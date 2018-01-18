!-------------------------------------------------------------------------------
!+

MODULE O3P_MOD_Declaration
!-------------------------------------------------------------------------------
!+Description: 
!     .
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.0     2016.05.11 Fisrt Code (R&D, SaeaSoft Co., Ltd.) 
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Modules:
! USE O3P_MOD_dummy
 USE O3P_MOD_share_init
 USE OMSAO_precision_module
!**********************************************************!
 IMPLICIT NONE

! ----------
! Parameters
! ----------

!**********************************************************!
CONTAINS

SUBROUTINE GEMS_O3P_SUB1_Declaration(retCode)
    IMPLICIT NONE

    INTEGER(KIND=4), INTENT(OUT)    :: retCode         ! return error code

    !-------------------------------------------
    INTEGER(KIND=4)                 :: errCode         ! error check code

    REAL (KIND=dp), DIMENSION(3)    :: fitcol
    REAL (KIND=dp), DIMENSION(3, 2) :: dfitcol

    errCode = 0

    print*, "+++ 1. Declaration"
    print*, "Here is GEMS_O3P_SUB1_Declaration() subroutine. !!!"
    print*, " "
    print*, " "

    !--- 공통 모듈 사용을 위한 초기화 작업
    CALL GEMS_O3P_MOD_share_init(errCode)
    IF ( errCode /= 0 ) THEN
        print*, "ERROR ! [공통 모듈 사용을 위한 초기화 작업]"
        retCode = errCode
        RETURN
    END IF
    retCode = errCode
    RETURN
END SUBROUTINE GEMS_O3P_SUB1_Declaration
END MODULE O3P_MOD_Declaration
