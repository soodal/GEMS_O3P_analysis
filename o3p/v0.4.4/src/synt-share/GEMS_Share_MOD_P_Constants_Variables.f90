!-------------------------------------------------------------------------------
!+Module to declare global variables and parameters for log processing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE Share_MOD_P_Constants_Variables 

!-------------------------------------------------------------------------------
!+Description: 
!     Module to declare GEMS global constants and variables for log processing.
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2014.12.23 Fisrt Code (YuGeun Ki, Seasoft)  
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Modules:
      ! USE Share_MOD_Types

      IMPLICIT NONE
      SAVE
 
!-------------------------------------------------------------------------------
! Global Constants:

!-------------------------------------------------------------------------------
! Global Variables:

      CHARACTER(LEN=160)    :: ctl_msg
      CHARACTER(LEN=160)    :: disp_msg
      CHARACTER(LEN=160)    :: disp_line

      CHARACTER(LEN=9)      :: LOGMSG
      INTEGER               :: llvl         ! Log Level(-2,-1,0,1,2,3), 
                                            ! -2=ERROR, -1=WARNING, 0=General, 2=INFO, 3=DETAIL INFO 
END MODULE Share_MOD_P_Constants_Variables

!-------------------------------------------------------------------------------
