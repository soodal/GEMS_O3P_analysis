!-------------------------------------------------------------------------------
!+

MODULE O3P_MOD_share_close
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
 USE Share_MOD_Radiance
 USE Share_l2_O3P_mod_write_read

!**********************************************************!
 IMPLICIT NONE

!**********************************************************!
 CONTAINS

!
!--  SUBROUTINE GEMS_O3P_MOD_share_close --
!
SUBROUTINE GEMS_O3P_MOD_share_close()
    IMPLICIT NONE

    INTEGER(KIND=4)                 :: errCode         ! error check code

    !--------------------
    !---- 메모리 해제 ---
    !--------------------
    print*, "+++ Fin. Share Memory Close "
    print*, " "
    print*, " "

    CALL GEMS_Share_deAllocMem4Rad                  ! memory deallocation of the module to calculate radiance

    errCode = GEMS_Share_MOD_L2O3P_MemDeAlloc()
    IF ( errCode /= 0 ) THEN
        WRITE(*, '( A36 )') "[L2 O3P] Memory DeAllocation Failed !"
        RETURN
    ENDIF

    RETURN
END SUBROUTINE GEMS_O3P_MOD_share_close

END MODULE O3P_MOD_share_close
