!-------------------------------------------------------------------------------
!+Module to get number times

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   MODULE NumTimes

!-------------------------------------------------------------------------------
!+Description: 
!     +Module to get number times in L1B file.
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2015. 3. 2 Fisrt  Code (Chiok  An, Seasoft) 
!-------------------------------------------------------------------------------

      TYPE lineds
        INTEGER(KIND=4) :: n_line_rug
        INTEGER(KIND=4) :: n_line_rvg
        INTEGER(KIND=4) :: n_line_irr
      END TYPE lineds
  
      TYPE(lineds)      :: nlines

      CONTAINS

      SUBROUTINE GEMS_NumTimes_Read( n_line )

!-------------------------------------------------------------------------------
!+Description: 
!      GEMS Share Subroutine to get number times in L1B file.
!
!
! Method:
!
!
! Input files:
!
!
! Output files:
!
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2015. 3. 2 Fisrt  Code (Chiok  An, Seasoft) 
!-------------------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=4) :: n_line

      CHARACTER(LEN=128)            :: ctl_fpath

      CHARACTER*200                 :: l1b_rug_file_path
      CHARACTER*200                 :: ert_uv1_dgrp_path
      CHARACTER*200                 :: ert_uv1_ggrp_path
      CHARACTER*200                 :: ert_uv1_sgrp_path
      CHARACTER*200                 :: ert_uv2_dgrp_path
      CHARACTER*200                 :: ert_uv2_ggrp_path
      CHARACTER*200                 :: ert_uv2_sgrp_path

      CHARACTER*200                 :: NumTimes_DataName

      INTEGER                       :: file_path_sz     ! File Path Size
      INTEGER                       :: grp_path_sz      ! Group Path Size
      INTEGER                       :: data_name_sz     ! Data Name Size
      !---------- HDF Read Various -----------

      CHARACTER*200                 :: file_path        ! File CPath
      CHARACTER*200                 :: grp_path         ! First CGroup Path
      CHARACTER*200                 :: data_name        ! Data Cname
      INTEGER                       :: pstart(3)        ! Start coordinate
      INTEGER                       :: pedge(3)         ! End coorinate
      INTEGER                       :: hdferr           ! status retrun value

      !
      ! External Function Declaration
      !
      INTEGER  GEMS_Share_Hdf5ReadData


      Namelist /L1B_RUG_File_List/l1b_rug_file_path, &
                                  ert_uv1_dgrp_path, &
                                  ert_uv1_ggrp_path, &
                                  ert_uv1_sgrp_path, &
                                  ert_uv2_dgrp_path, &
                                  ert_uv2_ggrp_path, &
                                  ert_uv2_sgrp_path

      Namelist /L1B_File_List_Value_Size/file_path_sz, &
                                         grp_path_sz,  &
                                         data_name_sz      

      Namelist /L1B_RUG_DATA_List_nt/NumTimes_DataName

      !$GEMS/[account]/v0.0/bin
      !$GEMS/share/conf
      ctl_fpath = '../../../share/conf/l1breadmdl.nml' 
      OPEN(10, file=trim(ctl_fpath), status='old')
      READ(10, L1B_RUG_File_List)
      READ(10, L1B_RUG_DATA_List_nt)
      READ(10, L1B_File_List_Value_Size)
      CLOSE(10)
!
!
      pstart(1) = -1
      file_path = l1b_rug_file_path
      grp_path = ert_uv1_sgrp_path
      data_name = NumTimes_DataName
!
!!--------- READ HDF --------------------------
!!
!!---------------------------------------------
!!--------- HDF Swath Attributes --------------
!!---------------------------------------------
!!          EXTERNAL FUNCTION USE
!!          - GEMS_Share_Hdf5ReadData
!!            1. Numtimes Read
!!---------
!
           hdferr = GEMS_Share_Hdf5ReadData(                     &
                                              file_path,         &
                                              file_path_sz,      &
                                              grp_path,          &
                                              grp_path_sz,       &
                                              data_name,         &
                                              data_name_sz,      &
                                              n_line,            &
                                              pstart,            &
                                              pedge               )


      END SUBROUTINE GEMS_NumTimes_Read

      SUBROUTINE GEMS_NumTimes_Read2( n_line , ctl_fpath)

!-------------------------------------------------------------------------------
!+Description: 
!      GEMS Share Subroutine to get number times in L1B file.
!
!
! Method:
!
!
! Input files:
!
!
! Output files:
!
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2016. 9. 27 Fisrt  Code (Chiok  An, Seasoft) 
!-------------------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=4) :: n_line

      CHARACTER(LEN=128),INTENT(IN) :: ctl_fpath

      CHARACTER*200                 :: l1b_rug_file_path
      CHARACTER*200                 :: ert_uv1_dgrp_path
      CHARACTER*200                 :: ert_uv1_ggrp_path
      CHARACTER*200                 :: ert_uv1_sgrp_path
      CHARACTER*200                 :: ert_uv2_dgrp_path
      CHARACTER*200                 :: ert_uv2_ggrp_path
      CHARACTER*200                 :: ert_uv2_sgrp_path

      CHARACTER*200                 :: NumTimes_DataName

      INTEGER                       :: file_path_sz     ! File Path Size
      INTEGER                       :: grp_path_sz      ! Group Path Size
      INTEGER                       :: data_name_sz     ! Data Name Size
      !---------- HDF Read Various -----------

      CHARACTER*200                 :: file_path        ! File CPath
      CHARACTER*200                 :: grp_path         ! First CGroup Path
      CHARACTER*200                 :: data_name        ! Data Cname
      INTEGER                       :: pstart(3)        ! Start coordinate
      INTEGER                       :: pedge(3)         ! End coorinate
      INTEGER                       :: hdferr           ! status retrun value

      !
      ! External Function Declaration
      !
      INTEGER  GEMS_Share_Hdf5ReadData


      Namelist /L1B_RUG_File_List/l1b_rug_file_path, &
                                  ert_uv1_dgrp_path, &
                                  ert_uv1_ggrp_path, &
                                  ert_uv1_sgrp_path, &
                                  ert_uv2_dgrp_path, &
                                  ert_uv2_ggrp_path, &
                                  ert_uv2_sgrp_path

      Namelist /L1B_File_List_Value_Size/file_path_sz, &
                                         grp_path_sz,  &
                                         data_name_sz      

      Namelist /L1B_RUG_DATA_List_nt/NumTimes_DataName

      !$GEMS/[account]/v0.0/bin
      !$GEMS/share/conf
      OPEN(10, file=trim(ctl_fpath), status='old')
      READ(10, L1B_RUG_File_List)
      READ(10, L1B_RUG_DATA_List_nt)
      READ(10, L1B_File_List_Value_Size)
      CLOSE(10)
!
!
      pstart(1) = -1
      file_path = l1b_rug_file_path
      grp_path = ert_uv1_sgrp_path
      data_name = NumTimes_DataName
!
!!--------- READ HDF --------------------------
!!
!!---------------------------------------------
!!--------- HDF Swath Attributes --------------
!!---------------------------------------------
!!          EXTERNAL FUNCTION USE
!!          - GEMS_Share_Hdf5ReadData
!!            1. Numtimes Read
!!---------
!
           hdferr = GEMS_Share_Hdf5ReadData(                     &
                                              file_path,         &
                                              file_path_sz,      &
                                              grp_path,          &
                                              grp_path_sz,       &
                                              data_name,         &
                                              data_name_sz,      &
                                              n_line,            &
                                              pstart,            &
                                              pedge               )


      END SUBROUTINE GEMS_NumTimes_Read2

      SUBROUTINE GEMS_NumTimes_Read3(ctl_fpath)

!-------------------------------------------------------------------------------
!+Description: 
!      GEMS Share Subroutine to get number times in L1B file.
!
!
! Method:
!
!
! Input files:
!
!
! Output files:
!
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2016. 9. 27 Fisrt  Code (Chiok  An, Seasoft) 
!-------------------------------------------------------------------------------

      IMPLICIT NONE


!      INTEGER(KIND=4) :: nline

      CHARACTER(LEN=128),INTENT(IN) :: ctl_fpath

      CHARACTER*200                 :: l1b_rug_file_path
      CHARACTER*200                 :: ert_uv1_dgrp_path
      CHARACTER*200                 :: ert_uv1_ggrp_path
      CHARACTER*200                 :: ert_uv1_sgrp_path
      CHARACTER*200                 :: ert_uv2_dgrp_path
      CHARACTER*200                 :: ert_uv2_ggrp_path
      CHARACTER*200                 :: ert_uv2_sgrp_path
      CHARACTER*200                 :: l1b_rvg_file_path
      CHARACTER*200                 :: ert_vis_dgrp_path
      CHARACTER*200                 :: ert_vis_ggrp_path
      CHARACTER*200                 :: ert_vis_sgrp_path
      CHARACTER*200                 :: l1b_irr_file_path
      CHARACTER*200                 :: sun_uv1_dgrp_path
      CHARACTER*200                 :: sun_uv1_ggrp_path
      CHARACTER*200                 :: sun_uv1_sgrp_path
      CHARACTER*200                 :: sun_uv2_dgrp_path
      CHARACTER*200                 :: sun_uv2_ggrp_path
      CHARACTER*200                 :: sun_uv2_sgrp_path
      CHARACTER*200                 :: sun_vis_dgrp_path
      CHARACTER*200                 :: sun_vis_ggrp_path
      CHARACTER*200                 :: sun_vis_sgrp_path


      CHARACTER*200                 :: NumTimes_DataName

      INTEGER                       :: file_path_sz     ! File Path Size
      INTEGER                       :: grp_path_sz      ! Group Path Size
      INTEGER                       :: data_name_sz     ! Data Name Size
      !---------- HDF Read Various -----------

      CHARACTER*200                 :: file_path        ! File CPath
      CHARACTER*200                 :: grp_path         ! First CGroup Path
      CHARACTER*200                 :: data_name        ! Data Cname
      INTEGER                       :: pstart(3)        ! Start coordinate
      INTEGER                       :: pedge(3)         ! End coorinate
      INTEGER                       :: hdferr           ! status retrun value

      !
      ! External Function Declaration
      !
      INTEGER  GEMS_Share_Hdf5ReadData


      Namelist /L1B_RUG_File_List/l1b_rug_file_path, &
                                  ert_uv1_dgrp_path, &
                                  ert_uv1_ggrp_path, &
                                  ert_uv1_sgrp_path, &
                                  ert_uv2_dgrp_path, &
                                  ert_uv2_ggrp_path, &
                                  ert_uv2_sgrp_path

      Namelist /L1B_RVG_File_List/l1b_rvg_file_path, &
                                  ert_vis_dgrp_path, &
                                  ert_vis_ggrp_path, &
                                  ert_vis_sgrp_path

      Namelist /L1B_IRR_File_List/l1b_irr_file_path, &
                                  sun_uv1_dgrp_path, &
                                  sun_uv1_ggrp_path, &
                                  sun_uv1_sgrp_path, &
                                  sun_uv2_dgrp_path, &
                                  sun_uv2_ggrp_path, &
                                  sun_uv2_sgrp_path, &
                                  sun_vis_dgrp_path, &
                                  sun_vis_ggrp_path, &
                                  sun_vis_sgrp_path

      Namelist /L1B_File_List_Value_Size/file_path_sz, &
                                         grp_path_sz,  &
                                         data_name_sz      

      Namelist /L1B_RUG_DATA_List_nt/NumTimes_DataName

!      Namelist /L1B_RVG_DATA_List_nt/NumTimes_DataName

!      Namelist /L1B_IRR_DATA_List_nt/NumTimes_DataName

      !$GEMS/[account]/v0.0/bin
      !$GEMS/share/conf
      OPEN(10, file=trim(ctl_fpath), status='old')
      READ(10, L1B_RUG_File_List)
      READ(10, L1B_RVG_File_List)
      READ(10, L1B_IRR_File_List)
      READ(10, L1B_RUG_DATA_List_nt)
!      READ(10, L1B_RVG_DATA_List_nt)
!      READ(10, L1B_IRR_DATA_List_nt)
      READ(10, L1B_File_List_Value_Size)
      CLOSE(10)

!
!
      pstart(1) = -1
      file_path = l1b_rug_file_path
      grp_path = ert_uv1_sgrp_path
      data_name = NumTimes_DataName
!
!!--------- READ HDF --------------------------
!!
!!---------------------------------------------
!!--------- HDF Swath Attributes --------------
!!---------------------------------------------
!!          EXTERNAL FUNCTION USE
!!          - GEMS_Share_Hdf5ReadData
!!            1. Numtimes Read
!!---------
!
           hdferr = GEMS_Share_Hdf5ReadData(                     &
                                              file_path,         &
                                              file_path_sz,      &
                                              grp_path,          &
                                              grp_path_sz,       &
                                              data_name,         &
                                              data_name_sz,      &
                                              nlines%n_line_rug,  &
                                              pstart,            &
                                              pedge               )

      file_path = l1b_rvg_file_path
      grp_path  = ert_vis_sgrp_path
!      data_name = NumTimes_DataName
!
!!--------- READ HDF --------------------------
!!
!!---------------------------------------------
!!--------- HDF Swath Attributes --------------
!!---------------------------------------------
!!          EXTERNAL FUNCTION USE
!!          - GEMS_Share_Hdf5ReadData
!!            1. Numtimes Read
!!---------
!
           hdferr = GEMS_Share_Hdf5ReadData(                     &
                                              file_path,         &
                                              file_path_sz,      &
                                              grp_path,          &
                                              grp_path_sz,       &
                                              data_name,         &
                                              data_name_sz,      &
                                              nlines%n_line_rvg,  &
                                              pstart,            &
                                              pedge               )

      file_path = l1b_irr_file_path
      grp_path  = sun_vis_sgrp_path
!      data_name = NumTimes_DataName
!
!!--------- READ HDF --------------------------
!!
!!---------------------------------------------
!!--------- HDF Swath Attributes --------------
!!---------------------------------------------
!!          EXTERNAL FUNCTION USE
!!          - GEMS_Share_Hdf5ReadData
!!            1. Numtimes Read
!!---------
!
           hdferr = GEMS_Share_Hdf5ReadData(                     &
                                              file_path,         &
                                              file_path_sz,      &
                                              grp_path,          &
                                              grp_path_sz,       &
                                              data_name,         &
                                              data_name_sz,      &
                                              nlines%n_line_irr,  &
                                              pstart,            &
                                              pedge               )


      END SUBROUTINE GEMS_NumTimes_Read3

   END MODULE
