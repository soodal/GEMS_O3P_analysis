!-------------------------------------------------------------------------------
!+

MODULE O3P_MOD_close
!-------------------------------------------------------------------------------
!+Description: 
!     .
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.0     2016.05.12 Fisrt Code (R&D, SaeaSoft Co., Ltd.) 
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Modules:

 USE GEMS_O3P_gemsdata_module,    ONLY: deallocate_o3p_l1b, deallocate_o3p_var, &
                                        gems_uv1, gems_uv2

!**********************************************************!
 IMPLICIT NONE

!**********************************************************!
 CONTAINS

!
!--  SUBROUTINE GEMS_O3P_MOD_close --
!
SUBROUTINE GEMS_O3P_SUB6_Proc_close()
    IMPLICIT NONE

    INTEGER(KIND=4)                 :: errCode         ! error check code
    INTEGER :: pge_error_status
    !--------------------
    !---- 메모리 해제 ---
    !--------------------

    print*, "+++ 6. Local Memory Close "
    print*, "Here is GEMS_O3P_SUB6_Proc_close() subroutine. !!!"
    print*, " "
    print*, " "

    CALL deallocate_o3p_l1b (gems_uv1, pge_error_status)
    CALL deallocate_o3p_l1b (gems_uv2, pge_error_status)
    CALL deallocate_o3p_var (pge_error_status)

    RETURN
END SUBROUTINE GEMS_O3P_SUB6_Proc_close

END MODULE O3P_MOD_close
