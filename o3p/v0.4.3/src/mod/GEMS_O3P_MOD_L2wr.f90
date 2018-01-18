!-------------------------------------------------------------------------------
!+

MODULE O3P_MOD_L2wr
!-------------------------------------------------------------------------------
!+Description: 
!     . write L2 Data of Cloud Algorithm
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2015.02.24 Fisrt Code (YuGeun Ki, SAEASoft) 
!-------------------------------------------------------------------------------


      
! Declarations:

! Modules used:

USE Share_MOD_Env
USE Share_MOD_Constants_Variables
USE Share_MOD_L2FileName
USE Share_l2_o3p_mod_write_read
USE GEMS_O3P_gemsdata_module, ONLY: gems_nx, gems_ny
USE ozprof_data_module,       ONLY: nlay
USE OMSAO_variables_module, ONLY:L2_filename
!**********************************************************!
IMPLICIT NONE

!**********************************************************!
CONTAINS

SUBROUTINE GEMS_O3P_L2wr(L2O3P_WR, retCode)
    TYPE(L2_o3p),    INTENT(INOUT)  :: L2O3P_WR
    INTEGER(KIND=4), INTENT(INOUT)  :: retCode    

    !-------------------------------------------------------!
    TYPE(O3P_ds)                    :: O3Pds
    CHARACTER(LEN=256)              :: wr_file_path         ! L2 출력파일이 생성될 경로 및 파일명
    INTEGER(KIND=4)                 :: file_path_sz = 256    

    INTEGER(KIND=4)                 :: hdferr       ! HDF5 Function return code 
    INTEGER(KIND=4)                 :: errCode      !    
    INTEGER(KIND=4)                 :: nlayp
    CHARACTER(LEN=128)              :: Lv2FileName  ! L2파일명을 저장

    CHARACTER(LEN=128)              :: nml_fpath    ! namelist파일이 존재하는 경로 및 파일명

    !
    ! External Function Declaration
    !
    INTEGER  GEMS_Share_Hdf5CreateL2File4O3P
    external GEMS_Share_Hdf5CreateL2File4O3P
    INTEGER GEMS_Share_Hdf5InitL2StorageSize
    external GEMS_Share_Hdf5InitL2StorageSize

    LOGMSG  =" PROC MSG"
    !------------------------------------------
    !---- 1. Gether L2 File Name         ---
    !------------------------------------------
    IF ( gds_O3pL2Env%setUp == 0 ) THEN
        llvl    =0
        disp_msg="did not setup L2 O3P Env"
        CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)
        RETURN
    END IF

    nlayp = nlay+1
    
    hdferr = GEMS_Share_Hdf5InitL2StorageSize(                  & ! need to consider GEMS_Share_MOD_Constants_Variables.f90
                gems_ny            ,                            & ! Size of cross track 
                gci_nline2_l2      ,                            & ! Not used for O3P
                gci_nwavel_l2      ,                            & ! Not used for O3P, used for AOD
                gci_nwavel2_l2     ,                            & ! Not used for O3P, used for O3T, SO2
                gci_nwavel3_l2     ,                            & ! O3P
                gci_nwavel4_l2     ,                            & ! Not used for O3P, use for ALBD
                gems_nx            ,                            & ! Size of cross track 
                gci_nxtrack2_l2    ,                            & ! Not Used for O3P
                nlay         ,                                  & ! Ozone Column layer
                nlay          ,                                 & ! Not Used for O3P
                nlayp        ,                                  & ! Atmospheric Profiles
                gci_nmon           ,                            & ! Not Used for O3P
                gci_nlat           ,                            & ! Not Used for O3P
                gci_nlon           ,                            & ! Not Used for O3P
                gci_ncomp          ,                            & ! Not Used for O3P
                gci_nResiduals     ,                            & ! Size of fitting residuals for O3P
                gci_nColumns       ,                            & ! Size of column o3 for O3P
                6                                               & 
      )
  
    CALL GEMS_Share_GetLv2FileName("TEST", "O3P", Lv2FileName)
    Lv2FileName = trim(gds_O3pL2Env%out_lv2_fname)    ! wasp
    wr_file_path  = trim(gds_O3pL2Env%out_lv2_fpath) // trim(Lv2FileName)
    !wr_file_path  = trim(gds_O3pL2Env%out_lv2_fpath) //'MPINO_YBIN01'//trim(Lv2FileName)
    !wr_file_path=trim(adjustl(l2_filename))
    WRITE(*,'(A)') ' => SaveFile:', adjustl(trim(wr_file_path))
    !------------------------------------------
    !---- 2. Create L2 File               ---
    !------------------------------------------
    hdferr = GEMS_Share_Hdf5CreateL2File4O3P(wr_file_path, file_path_sz)
    IF ( hdferr /= 0 ) THEN
        retCode = hdferr
        RETURN
    ENDIF

    nml_fpath  = trim(gds_O3pL2Env%nml_lv2_fpath) //trim(gds_O3pL2Env%nml_lv2_fname)

    !------------------------------------------
    !---- 3. Add values to L2 File          ---
    !------------------------------------------
    IF ( L2O3P_WR%status == 0 ) THEN
        llvl    =0
        disp_msg="did not allocate L2 O3P Memory"
        CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)
        RETURN
    END IF

    errCode = GEMS_Share_MOD_L2O3P_Write2(L2O3P_WR, nml_fpath, wr_file_path)
    IF ( errCode /= 0 ) THEN
        retCode = errCode
        RETURN
    END IF

    retCode = 0

    RETURN

END SUBROUTINE

END MODULE
