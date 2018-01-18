!-------------------------------------------------------------------------------
!+

MODULE O3P_MOD_share_init
!-------------------------------------------------------------------------------
!+Description: 
!     .
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.0     2016.05.12 Fisrt Code (Last First, ooo UNIV.) 
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Modules:
                                                ! Added by Saeasoft
 USE Share_MOD_Constants_Variables              ! GEMS Share Module for common variables              
 USE Share_MOD_Env
 USE Share_MOD_Radiance
 USE Share_l2_o3p_mod_write_read

!**********************************************************!
 IMPLICIT NONE

!**********************************************************!
 CONTAINS

!
!--  SUBROUTINE GEMS_O3P_MOD_share_init --
!
SUBROUTINE GEMS_O3P_MOD_share_init(retCode)
    IMPLICIT NONE

    INTEGER(KIND=4), INTENT(OUT)    :: retCode         ! return error code

    !-------------------------------------------
    CHARACTER(LEN=128)              :: cfg_fpath
    INTEGER(KIND=4)                 :: errCode         ! error check code

    errCode = 0

    print*, "+++ Init. Environment & Memory Initialize "
    print*, " "
    print*, " "

    !
    !--------  GEMS Share: Gether Configuration  ---- 
    !
    cfg_fpath = '../../../share/conf/gems.conf'
    CALL GEMS_Share_GetL2CldEnv(cfg_fpath, errCode)
    IF ( errCode /= 0 ) THEN
        WRITE(*, '( A36 )') "[Env O3P] Configuration Read Failed !"
        retCode = -9001
        RETURN
    ENDIF

    !
    !--------  GEMS Share: Global Variables  --------- 
    !
    CALL GEMS_Share_Init_GlobalConstants
!    CALL GEMS_Share_Init_L2O3P_GlobalConstants

    CALL GEMS_Share_init4RadianceCalc(gci_nwavel3, gci_nxtrack2,    &
                                      gci_ntimes, gci_nwavelcoef)     ! for using the module to calculate radiance
    CALL GEMS_Share_allocMem4Rad                                      ! for memory allocation of the module to calculate radiance
 
    !
    !-------- GEMS Share: 동적 메모리 할당  --------- 
    !
    errCode = GEMS_Share_MOD_L2O3P_MemAlloc()
    IF ( errCode /= 0 ) THEN
        WRITE(*, '( A36 )') "[L2 O3P] Memory Allocation Failed !"
        retCode = -9002
        RETURN
    ENDIF
    !

    retCode = errCode

    RETURN
END SUBROUTINE GEMS_O3P_MOD_share_init


END MODULE O3P_MOD_share_init
