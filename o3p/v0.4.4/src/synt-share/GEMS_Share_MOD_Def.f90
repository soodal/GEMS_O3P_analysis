!-------------------------------------------------------------------------------
!+Module to define the Default Value of GEMS Algorithm 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE Share_MOD_Def

!-------------------------------------------------------------------------------
!+Description: 
!      Module to define the Default Value of GEMS Algorithm 
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2015.04.03 Fisrt Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Modules:

!**********************************************************!
IMPLICIT NONE

!-------------------------------------------------------------------------------
!+Global Default Value Definition:
    REAL(KIND=8),    PARAMETER      :: DEF_dVAL = -999
    REAL(KIND=4),    PARAMETER      :: DEF_fVAL = -999
    INTEGER(KIND=8), PARAMETER      :: DEF_lVAL = -999
    INTEGER(KIND=4), PARAMETER      :: DEF_iVAL = -999
    INTEGER(KIND=2), PARAMETER      :: DEF_sVAL = -32767
    INTEGER(KIND=1), PARAMETER      :: DEF_cVAL = -127

!--------------------------------
!+External Function Declaration For Using GEMS HDF5 API
!    INTEGER  GEMS_Share_Hdf5ReadData
!    INTEGER  GEMS_Share_Hdf5WriteData
!    INTEGER  GEMS_Share_Hdf5ReadAttr
!    INTEGER  GEMS_Share_Hdf5WriteAttr
!    INTEGER  GEMS_Share_Hdf5CreateL2File
!    INTEGER  GEMS_Share_Hdf5CreateL1BFile
!    external GEMS_Share_Hdf5ReadData
!    external GEMS_Share_Hdf5WriteData
!    external GEMS_Share_Hdf5ReadAttr
!    external GEMS_Share_Hdf5WriteAttr
!    external GEMS_Share_Hdf5CreateL2File
!    external GEMS_Share_Hdf5CreateL1BFile
!--------------------------------

!    INTEGER  GEMS_Share_xxx
!    external GEMS_Share_xxx
!    INTEGER  GEMS_Share_xxx2
!    external GEMS_Share_xxx2

  CONTAINS 

!-------------------------------------------------------------------------------
!+SUBROUTINE:


END MODULE Share_MOD_Def

!-------------------------------------------------------------------------------
