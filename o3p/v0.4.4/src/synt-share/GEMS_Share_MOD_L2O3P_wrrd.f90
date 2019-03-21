
!-------------------------------------------------------------------------------
!+Module to write and read GEMS L2File

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE Share_l2_o3p_mod_write_read

!-------------------------------------------------------------------------------
!+Description: 
!     Module to wirte and read GEMS L2O3P File
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2014.12.11 Fisrt Code (HJ Lee, Seasoft) 
! 0.2     2017.01.31 L2데이터정의서 v0.5버전 적용(YuGeun) 
!-------------------------------------------------------------------------------

      USE Share_MOD_Constants_Variables, ONLY:  gci_nline_l2,       &
                                                gci_nxtrack2_l2,    &
                                                gci_nlayer,         &
                                                gci_nlayerp,        &
                                                gci_nwavel3_l2,     &
                                                gci_nResiduals,     &
                                                gci_nColumns
      USE Share_MOD_P_Constants_Variables, ONLY: LOGMSG, llvl, ctl_msg, &
                                                 disp_msg, disp_line
      USE Share_MOD_Log
      USE Share_MOD_Def
!**********************************************************!
      IMPLICIT NONE

      PUBLIC  :: GEMS_Share_MOD_L2O3P_Write
      PUBLIC  :: GEMS_Share_MOD_L2O3P_Read
      PUBLIC  :: GEMS_Share_MOD_L2O3P_Write2
      PUBLIC  :: read_namelist_o3p

      INTEGER(KIND=4), PRIVATE      :: nlevel_o3p
      INTEGER(KIND=4), PRIVATE      :: nxtrack2_o3p
      INTEGER(KIND=4), PRIVATE      :: nline
      INTEGER(KIND=4), PRIVATE      :: nlayer_o3p
      INTEGER(KIND=4), PRIVATE      :: nlayerp_o3p
      INTEGER(KIND=4), PRIVATE      :: nwavel3_o3p
      INTEGER(KIND=4), PRIVATE      :: nResiduals_o3p
      INTEGER(KIND=4), PRIVATE      :: nColumns_o3p

    !----------O3P----------!
     ! <O3Profile>

      !--- Data Fields --- 
      TYPE :: o3p_d_fld
       REAL   (KIND=8) , POINTER    :: AvgK           (:,:,:,:) !AveragingKernel            
       REAL   (KIND=8) , POINTER    :: OzRet          (:,:,:)   !O3                         
       REAL   (KIND=8) , POINTER    :: OzAp           (:,:,:)   !O3Apriori                  
       REAL   (KIND=8) , POINTER    :: OzApErr        (:,:,:)   !O3AprioriError             
       REAL   (KIND=8) , POINTER    :: OzNErr         (:,:,:)   !O3RandomNoiseError              
       REAL   (KIND=8) , POINTER    :: OzSMErr        (:,:,:)   !O3SolutionError              
       REAL   (KIND=8) , POINTER    :: ColAmtOz       (:,:,:)   !ColumnAmountO3             
       REAL   (KIND=8) , POINTER    :: DFS            (:,:,:)   !DegreesOfFreedomForSignal  
       REAL   (KIND=4) , POINTER    :: CldP           (:,:)     !CloudPressure              
       REAL   (KIND=8) , POINTER    :: EffCldFrac     (:,:)     !EffectiveCloudFractionUV  
       REAL   (KIND=8) , POINTER    :: SfcAlb         (:,:)     !TerrainReflectivityUV     
       INTEGER(KIND=4) , POINTER    :: PrsQFlag       (:,:)     !ProcessingQualityFlags     
       INTEGER(KIND=2) , POINTER    :: NumIter        (:,:)     !NumberOfIterations         
       REAL   (KIND=8) , POINTER    :: AvgRes         (:,:,:)   !ResidualsOfFit             
       REAL   (KIND=8) , POINTER    :: RMS            (:,:)     !RootMeanSquareErrorOfFit   
      END TYPE o3p_d_fld

      !--- Geolocation Fields --- 
      TYPE :: o3p_g_fld
       REAL   (KIND=4) , POINTER    :: TrppsP         (:,:)     !TropopausePressure
       REAL   (KIND=4) , POINTER    :: P              (:,:,:)   !Pressure                   
       REAL   (KIND=4) , POINTER    :: Alt            (:,:,:)   !Altitude                   
       REAL   (KIND=4) , POINTER    :: TEMP           (:,:,:)   !Temperature                
       REAL   (KIND=4) , POINTER    :: TerrP          (:,:)     !TerrainPressure
       REAL   (KIND=4) , POINTER    :: LAT            (:,:)     !Latitude                   
       REAL   (KIND=4) , POINTER    :: LON            (:,:)     !Longitude                  
       INTEGER(KIND=4) , POINTER    :: Line_tmp       (:,:)     !Line                ! add 2016.09.23     line -> line_tmp   : geun
       INTEGER(KIND=4) , POINTER    :: Pix_tmp        (:,:)     !Pix                 ! add 2016.09.23     for unknown pointer error  
       REAL   (KIND=4) , POINTER    :: SolZenAng      (:,:)     !SolarZenithAngle           
       REAL   (KIND=4) , POINTER    :: ViewZenAng     (:,:)     !ViewingZenithAngle         
       REAL   (KIND=4) , POINTER    :: RelAziAng      (:,:)     !RelativeAzimuthAngle  
       REAL   (KIND=8) , POINTER    :: Time           (:)       !Time                       
      END TYPE o3p_g_fld

      TYPE :: L2_o3p
       TYPE(o3p_d_fld)              :: o3p_dfld
       TYPE(o3p_g_fld)              :: o3p_gfld
       INTEGER(KIND=1)              :: status   = 0             ! 메모리를 할당 받았는지 여부를 확인
      END TYPE L2_o3p

      TYPE(L2_o3p)                  :: gds_L2O3P

     ! O3P datasets in namelist
      TYPE :: O3P_ds
       CHARACTER*200                 :: l2_o3p_file_path
       CHARACTER*200                 :: o3p_dgrp_path
       CHARACTER*200                 :: o3p_ggrp_path

!      -- NEW DATA --
       CHARACTER*200                 :: AvgKernel_DataName
       CHARACTER*200                 :: O3_DataName
       CHARACTER*200                 :: O3Apriori_DataName
       CHARACTER*200                 :: O3AprioriErr_DataName
       CHARACTER*200                 :: O3RndmNsErr_DataName
       CHARACTER*200                 :: O3SolutionErr_DataName
       CHARACTER*200                 :: ColAmountO3_DataName
       CHARACTER*200                 :: DegOfFreedom_DataName
       CHARACTER*200                 :: CldPressure_DataName
       CHARACTER*200                 :: EffCldFracUV_DataName
       CHARACTER*200                 :: TerrRefltUV_DataName
       CHARACTER*200                 :: ProcessQFlag_DataName
       CHARACTER*200                 :: NumOfIteration_DataName
       CHARACTER*200                 :: ResidualsOfFit_DataName
       CHARACTER*200                 :: RMSErrOfFit_DataName
!      -- NEW GEO --
       CHARACTER*200                 :: TrppsPressure_DataName
       CHARACTER*200                 :: Pressure_DataName
       CHARACTER*200                 :: Altitude_DataName
       CHARACTER*200                 :: Temperature_DataName
       CHARACTER*200                 :: TerrPressure_DataName
       CHARACTER*200                 :: Latitude_DataName
       CHARACTER*200                 :: Longitude_DataName
       CHARACTER*200                 :: Line_DataName         ! add 2016.09.23
       CHARACTER*200                 :: Pix_DataName          ! add 2016.09.23
       CHARACTER*200                 :: SolarZenAng_DataName
       CHARACTER*200                 :: ViewZenAng_DataName
       CHARACTER*200                 :: RelAziAng_DataName
       CHARACTER*200                 :: Time_DataName

       INTEGER                       :: file_path_sz
       INTEGER                       :: grp_path_sz
       INTEGER                       :: data_name_sz 
      END TYPE O3P_ds

      TYPE(O3P_ds)                         :: O3Pdss

CONTAINS

!--------------------------------

!-------------------------------------------------------------------------------
!+Description: 
!      FUNCTION GEMS_Share_MOD_L2_O3P_Write_Read .
!
!
! Method:

! write and read  L2O3P files:
!       L2_o3p    : the Level 2 O3P data structure
!       ctl_fpath : the control file path
!
! output files:
!       status    : the result code of L2O3P write and read
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2014.12.23  Fisrt Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!L2 O3P write function   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Declarations:
      FUNCTION GEMS_Share_MOD_L2O3P_Write(L2O3P, O3Pds, ctl_fpath) RESULT (status)

! input/output parameters:
      CHARACTER(LEN=128),   INTENT(INOUT) :: ctl_fpath
      TYPE(L2_o3p),         INTENT(INOUT) :: L2O3P
      TYPE(O3P_ds),         INTENT(INOUT) :: O3Pds

! Local parameters For L2 Reading:
      CHARACTER*200                 :: data_name
      CHARACTER*200                 :: grp_path

      INTEGER                       :: pstart(3)
      INTEGER                       :: pedge(3)

      INTEGER                       :: ierr        ! Error code 
      INTEGER                       :: hdferr      ! HDF5 Read Error code 
      INTEGER*8                     :: iAddr       ! memory address#

      INTEGER (KIND = 4)            :: status

! Local parameters For Logging:
      INTEGER       :: i, j, k

      INTEGER(1), PARAMETER :: X_FF = INT(Z'FF',1)

!--------------------------------
      !
      ! External Function Declaration
      !
      INTEGER  GEMS_Share_Hdf5ReadData
      INTEGER  GEMS_Share_Hdf5WriteData
      INTEGER  GEMS_Share_Hdf5ReadAttr
      INTEGER  GEMS_Share_Hdf5WriteAttr
      INTEGER  GEMS_Share_Hdf5CreateL2File
      INTEGER  GEMS_Share_Hdf5CreateL1BFile
      external GEMS_Share_Hdf5ReadData
      external GEMS_Share_Hdf5WriteData
      external GEMS_Share_Hdf5ReadAttr
      external GEMS_Share_Hdf5WriteAttr
      external GEMS_Share_Hdf5CreateL2File
      external GEMS_Share_Hdf5CreateL1BFile
!--------------------------------


!--------------------------------
!---- FUNCTION Start -------------
!--------------------------------
      LOGMSG=" PROC MSG"
      !CALL GEMS_Share_MOD_log(llvl, '--- Start L2O3P WRITE FUNCTION ---', 'START MSG')

!
!--------------------------------
!
      pstart(1) = -1                                            ! HDF5 데이터셑의 모든 데이터값을 읽는다는 의미

      CALL read_namelist_o3p(ctl_fpath, O3Pds)

      !----------      L2O3P     ----------!
      !--------------------------------------
      ! <O3 Profile: Data Fields>
      !--------------------------------------

      disp_msg = "L2O3P WRITE FILE PATH="//O3Pds%l2_o3p_file_path
      !CALL GEMS_Share_MOD_log(llvl, disp_msg,  LOGMSG)

      !CALL GEMS_Share_MOD_log(llvl, "O3 Profile :Data Fields",  LOGMSG)

      !-----------------  AveragingKernel  ---------------
      data_name = O3Pds%AvgKernel_DataName
      grp_path  = O3Pds%o3p_dgrp_path
      hdferr = GEMS_Share_Hdf5WriteData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%AvgK, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%AvgK, nline,nxtrack2_o3p,nlayer_o3p,nlayer_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  O3  ---------------
      data_name = O3Pds%O3_DataName
      hdferr = GEMS_Share_Hdf5WriteData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%OzRet, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%OzRet, nline,nxtrack2_o3p,nlayer_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  O3Apriori  ---------------
      data_name = O3Pds%O3Apriori_DataName
      hdferr = GEMS_Share_Hdf5WriteData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%OzAp, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%OzAp, nline,nxtrack2_o3p,nlayer_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  O3AprioriError  ---------------
      data_name = O3Pds%O3AprioriErr_DataName
      hdferr = GEMS_Share_Hdf5WriteData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%OzApErr, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%OzApErr, nline,nxtrack2_o3p,nlayer_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  O3 Random Noise Error  ---------------
      data_name = O3Pds%O3RndmNsErr_DataName
      hdferr = GEMS_Share_Hdf5WriteData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%OzNErr, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%OzNErr, nline,nxtrack2_o3p,nlayer_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  O3 Solution Error  ---------------
      data_name = O3Pds%O3SolutionErr_DataName
      hdferr = GEMS_Share_Hdf5WriteData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%OzSMErr, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%OzSMErr,nline,nxtrack2_o3p,nlayer_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  ColumnAmountO3  ---------------
      data_name = O3Pds%ColAmountO3_DataName
      hdferr = GEMS_Share_Hdf5WriteData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%ColAmtOz, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%ColAmtOz, nline,nxtrack2_o3p,nColumns_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  DegreesOfFreedomForSignal  ---------------
      data_name = O3Pds%DegOfFreedom_DataName
      hdferr = GEMS_Share_Hdf5WriteData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%DFS, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%DFS, nline,nxtrack2_o3p,nColumns_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  CloudPressure  ---------------
      data_name = O3Pds%CldPressure_DataName
      hdferr = GEMS_Share_Hdf5WriteData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%CldP, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%CldP, nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  EffectiveCloudFractionUV  ---------------
      data_name = O3Pds%EffCldFracUV_DataName
      hdferr = GEMS_Share_Hdf5WriteData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%EffCldFrac, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%EffCldFrac, nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  TerrainReflectivityUV  ---------------
      data_name = O3Pds%TerrRefltUV_DataName
      hdferr = GEMS_Share_Hdf5WriteData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%SfcAlb, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%SfcAlb,nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  ProcessingQualityFlags  ---------------
      data_name = O3Pds%ProcessQFlag_DataName
      hdferr = GEMS_Share_Hdf5WriteData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%PrsQFlag, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%PrsQFlag,nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  NumberOfIterations  ---------------
      data_name = O3Pds%NumOfIteration_DataName      
      hdferr = GEMS_Share_Hdf5WriteData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%NumIter, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%NumIter, nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  ResidualsOfFit  ---------------
      data_name = O3Pds%ResidualsOfFit_DataName
      hdferr = GEMS_Share_Hdf5WriteData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%AvgRes, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%AvgRes,nline,nxtrack2_o3p,nResiduals_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  RootMeanSquareErrorOfFit  ---------------
      data_name = O3Pds%RMSErrOfFit_DataName
      hdferr = GEMS_Share_Hdf5WriteData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%RMS, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%RMS,nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !--------------------------------------
      ! <O3Profile: Geolocation Fields>
      !--------------------------------------
      !CALL GEMS_Share_MOD_log(llvl, "O3Profile: Geolocation Fields", LOGMSG)

      !-----------------  TropopausePressure  ---------------
      data_name = O3Pds%TrppsPressure_DataName
      grp_path  = O3Pds%o3p_ggrp_path
      hdferr = GEMS_Share_Hdf5WriteData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%TrppsP, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%TrppsP, nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  Pressure  ---------------
      data_name = O3Pds%Pressure_DataName
      hdferr = GEMS_Share_Hdf5WriteData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%P, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%P, nline,nxtrack2_o3p,nlevel_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  Altitude  ---------------
      data_name = O3Pds%Altitude_DataName
      hdferr = GEMS_Share_Hdf5WriteData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%Alt, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%Alt, nline,nxtrack2_o3p,nlevel_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  Temperature  ---------------
      data_name = O3Pds%Temperature_DataName
      hdferr = GEMS_Share_Hdf5WriteData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%TEMP, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%TEMP,nline,nxtrack2_o3p,nlevel_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  TerrainPressure  ---------------
      data_name = O3Pds%TerrPressure_DataName
      hdferr = GEMS_Share_Hdf5WriteData(O3Pds%l2_o3p_file_path,O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz,L2O3P%o3p_gfld%TerrP, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl,L2O3P%o3p_gfld%TerrP,nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  Latitude  ---------------
      data_name = O3Pds%Latitude_DataName
      hdferr = GEMS_Share_Hdf5WriteData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%LAT, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%LAT, nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  Longitude  ---------------
      data_name = O3Pds%Longitude_DataName
      hdferr = GEMS_Share_Hdf5WriteData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%LON, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%LON, nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  Line       ---------------
      data_name = O3Pds%Line_DataName
      hdferr = GEMS_Share_Hdf5WriteData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%Line_tmp, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%Line_tmp, nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  Pix        ---------------
      data_name = O3Pds%Pix_DataName
      hdferr = GEMS_Share_Hdf5WriteData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%Pix_tmp, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%Pix_tmp, nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SolarZenithAngle  ---------------
      data_name = O3Pds%SolarZenAng_DataName
      hdferr = GEMS_Share_Hdf5WriteData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%SolZenAng, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%SolZenAng, nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  ViewingZenithAngle  ---------------
      data_name = O3Pds%ViewZenAng_DataName
      hdferr = GEMS_Share_Hdf5WriteData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%ViewZenAng, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%ViewZenAng, nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  RelativeAzimuthAngle  ---------------
      data_name = O3Pds%RelAziAng_DataName
      hdferr = GEMS_Share_Hdf5WriteData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%RelAziAng, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%RelAziAng, nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  Time  ---------------
      data_name = O3Pds%Time_DataName
      hdferr = GEMS_Share_Hdf5WriteData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%Time, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%Time, nline, data_name)
      ELSE
            GOTO 2999
      ENDIF

      write(disp_msg,'(" GEMS HDF5 File READ OK, hdferr=",i6)') hdferr
      !CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)

!--------------------------------
!---- END FUNCTION -------------
!--------------------------------

      !CALL GEMS_Share_MOD_log(llvl, '---   END L2 WRITE FUNCTION ---', '  END MSG')
      !CALL GEMS_Share_MOD_log(llvl, '.', LOGMSG)
      status = 0

      RETURN

999   CONTINUE
      write(disp_msg,'( A20, "Memory Allocation Failed, ierr=",i6)') ctl_msg,ierr
      CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)
      status = ierr

      RETURN

2999  CONTINUE
      write(disp_msg,'(" !!! GEMS HDF5 File WRITE ERROR !!! ", "[", A30, "],hdferr=",i6)') data_name, hdferr
      CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)
      status = hdferr

      RETURN

      END FUNCTION GEMS_Share_MOD_L2O3P_Write

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!L2 O3P read function   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Declarations:
      FUNCTION GEMS_Share_MOD_L2O3P_Read(L2O3P, O3Pds, ctl_fpath) RESULT (status)

! input/output parameters:
      CHARACTER(LEN=128),   INTENT(INOUT) :: ctl_fpath
      TYPE(L2_o3p),         INTENT(INOUT) :: L2O3P
      TYPE(O3P_ds),         INTENT(INOUT) :: O3Pds

! Local parameters For L2 Reading:
      CHARACTER*200                 :: data_name
      CHARACTER*200                 :: grp_path

      INTEGER                       :: pstart(3)
      INTEGER                       :: pedge(3)

      INTEGER                       :: ierr        ! Error code 
      INTEGER                       :: hdferr      ! HDF5 Read Error code 
      INTEGER*8                     :: iAddr       ! memory address#

      INTEGER (KIND = 4)            :: status

! Local parameters For Logging:
      INTEGER       :: i, j, k

      INTEGER(1), PARAMETER :: X_FF = INT(Z'FF',1)

!--------------------------------
      !
      ! External Function Declaration
      !
      INTEGER  GEMS_Share_Hdf5ReadData
      INTEGER  GEMS_Share_Hdf5WriteData
      INTEGER  GEMS_Share_Hdf5ReadAttr
      INTEGER  GEMS_Share_Hdf5WriteAttr
      INTEGER  GEMS_Share_Hdf5CreateL2File
      INTEGER  GEMS_Share_Hdf5CreateL1BFile
      external GEMS_Share_Hdf5ReadData
      external GEMS_Share_Hdf5WriteData
      external GEMS_Share_Hdf5ReadAttr
      external GEMS_Share_Hdf5WriteAttr
      external GEMS_Share_Hdf5CreateL2File
      external GEMS_Share_Hdf5CreateL1BFile
!--------------------------------


!--------------------------------

!--------------------------------
!---- FUNCTION Start -------------
!--------------------------------
      LOGMSG=" PROC MSG"
      !CALL GEMS_Share_MOD_log(llvl, '--- Start L2O3P READ FUNCTION  ---', 'STARTMSG')

!--------------------------------
!
      pstart(1) = -1                                            ! HDF5 데이터셑의 모든 데이터값을 읽는다는 의미

      CALL read_namelist_o3p(ctl_fpath, O3Pds)

      !----------      L2O3P     ----------!
      !--------------------------------------
      ! <O3 Profile: Data Fields>
      !--------------------------------------

      !disp_msg = "L2O3P READ FILE PATH="//O3Pds%l2_o3p_file_path
      !CALL GEMS_Share_MOD_log(llvl, disp_msg,  LOGMSG)

      !CALL GEMS_Share_MOD_log(llvl, "O3 Profile :Data Fields",  LOGMSG)

      !-----------------  AveragingKernel  ---------------
      IF (ASSOCIATED(L2O3P%o3p_dfld%AvgK)) DEALLOCATE(L2O3P%o3p_dfld%AvgK)
      data_name = O3Pds%AvgKernel_DataName
      grp_path  = O3Pds%o3p_dgrp_path
      ALLOCATE (L2O3P%o3p_dfld%AvgK(nlayer_o3p,nlayer_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                       data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%AvgK, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%AvgK, nline,nxtrack2_o3p,nlayer_o3p,nlayer_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  O3  ---------------
      IF (ASSOCIATED(L2O3P%o3p_dfld%OzRet)) DEALLOCATE(L2O3P%o3p_dfld%OzRet)
      data_name = O3Pds%O3_DataName
      ALLOCATE (L2O3P%o3p_dfld%OzRet(nlayer_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                       data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%OzRet, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%OzRet, nline,nxtrack2_o3p,nlayer_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  O3Apriori  ---------------
      IF (ASSOCIATED(L2O3P%o3p_dfld%OzAp)) DEALLOCATE(L2O3P%o3p_dfld%OzAp)
      data_name = O3Pds%O3Apriori_DataName
      ALLOCATE (L2O3P%o3p_dfld%OzAp(nlayer_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                       data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%OzAp, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%OzAp,nline,nxtrack2_o3p,nlayer_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  O3AprioriError  ---------------
      IF (ASSOCIATED(L2O3P%o3p_dfld%OzApErr)) DEALLOCATE(L2O3P%o3p_dfld%OzApErr)
      data_name = O3Pds%O3AprioriErr_DataName
      ALLOCATE (L2O3P%o3p_dfld%OzApErr(nlayer_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                       data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%OzApErr, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%OzApErr, nline,nxtrack2_o3p,nlayer_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  O3 Random Noise Error  ---------------
      IF (ASSOCIATED(L2O3P%o3p_dfld%OzNErr)) DEALLOCATE(L2O3P%o3p_dfld%OzNErr)
      data_name = O3Pds%O3RndmNsErr_DataName
      ALLOCATE (L2O3P%o3p_dfld%OzNErr(nlayer_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                       data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%OzNErr, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%OzNErr, nline,nxtrack2_o3p,nlayer_o3p,data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  O3 Solution Error  ---------------
      IF (ASSOCIATED(L2O3P%o3p_dfld%OzSMErr)) DEALLOCATE(L2O3P%o3p_dfld%OzSMErr)
      data_name = O3Pds%O3SolutionErr_DataName
      ALLOCATE (L2O3P%o3p_dfld%OzSMErr(nlayer_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path,O3Pds%grp_path_sz, &
                                       data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%OzSMErr, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%OzSMErr,nline,nxtrack2_o3p,nlayer_o3p,data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  ColumnAmountO3  ---------------
      IF (ASSOCIATED(L2O3P%o3p_dfld%ColAmtOz)) DEALLOCATE(L2O3P%o3p_dfld%ColAmtOz)
      data_name = O3Pds%ColAmountO3_DataName
      ALLOCATE(L2O3P%o3p_dfld%ColAmtOz(nColumns_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                       data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%ColAmtOz, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%ColAmtOz, nline,nxtrack2_o3p,nColumns_o3p,data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  DegreesOfFreedomForSignal  ---------------
      IF (ASSOCIATED(L2O3P%o3p_dfld%DFS)) DEALLOCATE(L2O3P%o3p_dfld%DFS)
      data_name = O3Pds%DegOfFreedom_DataName
      ALLOCATE (L2O3P%o3p_dfld%DFS(nColumns_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                       data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%DFS, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%DFS, nline,nxtrack2_o3p,nColumns_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  CloudPressure  ---------------
      IF (ASSOCIATED(L2O3P%o3p_dfld%CldP)) DEALLOCATE(L2O3P%o3p_dfld%CldP)
      data_name = O3Pds%CldPressure_DataName
      ALLOCATE(L2O3P%o3p_dfld%CldP(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                       data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%CldP, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%CldP, nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  EffectiveCloudFractionUV  ---------------
      IF (ASSOCIATED(L2O3P%o3p_dfld%EffCldFrac)) DEALLOCATE(L2O3P%o3p_dfld%EffCldFrac)
      data_name = O3Pds%EffCldFracUV_DataName
      ALLOCATE (L2O3P%o3p_dfld%EffCldFrac(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                       data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%EffCldFrac, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%EffCldFrac, nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  TerrainReflectivityUV  ---------------
      IF (ASSOCIATED(L2O3P%o3p_dfld%SfcAlb)) DEALLOCATE(L2O3P%o3p_dfld%SfcAlb)
      data_name = O3Pds%TerrRefltUV_DataName
      ALLOCATE (L2O3P%o3p_dfld%SfcAlb(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                       data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%SfcAlb, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%SfcAlb,nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  ProcessingQualityFlags  ---------------
      IF (ASSOCIATED(L2O3P%o3p_dfld%PrsQFlag)) DEALLOCATE(L2O3P%o3p_dfld%PrsQFlag)
      data_name = O3Pds%ProcessQFlag_DataName
      ALLOCATE (L2O3P%o3p_dfld%PrsQFlag(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                       data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%PrsQFlag, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%PrsQFlag,nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  NumberOfIterations  ---------------
      IF (ASSOCIATED(L2O3P%o3p_dfld%NumIter)) DEALLOCATE(L2O3P%o3p_dfld%NumIter)
      data_name = O3Pds%NumOfIteration_DataName
      ALLOCATE (L2O3P%o3p_dfld%NumIter(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                       data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%NumIter, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%NumIter, nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  ResidualsOfFit  ---------------
      IF (ASSOCIATED(L2O3P%o3p_dfld%AvgRes)) DEALLOCATE(L2O3P%o3p_dfld%AvgRes)
      data_name = O3Pds%ResidualsOfFit_DataName
      ALLOCATE (L2O3P%o3p_dfld%AvgRes(nResiduals_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                       data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%AvgRes, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%AvgRes,nline,nxtrack2_o3p,nResiduals_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  RootMeanSquareErrorOfFit  ---------------
      IF (ASSOCIATED(L2O3P%o3p_dfld%RMS)) DEALLOCATE(L2O3P%o3p_dfld%RMS)
      data_name = O3Pds%RMSErrOfFit_DataName
      ALLOCATE (L2O3P%o3p_dfld%RMS(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                       data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%RMS, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%RMS,nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !--------------------------------------
      ! <O3Profile: Geolocation Fields>
      !--------------------------------------
      !CALL GEMS_Share_MOD_log(llvl, "O3Profile: Geolocation Fields", LOGMSG)

      !-----------------  TropopausePressure  ---------------
      IF (ASSOCIATED(L2O3P%o3p_gfld%TrppsP)) DEALLOCATE(L2O3P%o3p_gfld%TrppsP)
      data_name = O3Pds%TrppsPressure_DataName
      grp_path  = O3Pds%o3p_ggrp_path
      ALLOCATE (L2O3P%o3p_gfld%TrppsP(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                       data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%TrppsP, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%TrppsP, nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  Pressure  ---------------
      IF (ASSOCIATED(L2O3P%o3p_gfld%P)) DEALLOCATE(L2O3P%o3p_gfld%P)
      data_name = O3Pds%Pressure_DataName
      ALLOCATE (L2O3P%o3p_gfld%P(nlevel_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                       data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%P, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%P, nline,nxtrack2_o3p,nlevel_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  Altitude  ---------------
      IF (ASSOCIATED(L2O3P%o3p_gfld%Alt)) DEALLOCATE(L2O3P%o3p_gfld%Alt)
      data_name = O3Pds%Altitude_DataName
!      grp_path  = O3Pds%o3p_ggrp_path
      ALLOCATE (L2O3P%o3p_gfld%Alt(nlevel_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                       data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%Alt, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%Alt, nline,nxtrack2_o3p,nlevel_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  Temperature  ---------------
      IF (ASSOCIATED(L2O3P%o3p_gfld%TEMP)) DEALLOCATE(L2O3P%o3p_gfld%TEMP)
      data_name = O3Pds%Temperature_DataName
      ALLOCATE (L2O3P%o3p_gfld%TEMP(nlevel_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                       data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%TEMP, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%TEMP,nline,nxtrack2_o3p,nlevel_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  TerrainPressure  ---------------
      IF (ASSOCIATED(L2O3P%o3p_gfld%TerrP)) DEALLOCATE(L2O3P%o3p_gfld%TerrP)
      data_name = O3Pds%TerrPressure_DataName
      ALLOCATE (L2O3P%o3p_gfld%TerrP(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(O3Pds%l2_o3p_file_path,O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                       data_name, O3Pds%data_name_sz,L2O3P%o3p_gfld%TerrP, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl,L2O3P%o3p_gfld%TerrP, nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  Latitude  ---------------
      IF (ASSOCIATED(L2O3P%o3p_gfld%LAT)) DEALLOCATE(L2O3P%o3p_gfld%LAT)
      data_name = O3Pds%Latitude_DataName
      ALLOCATE (L2O3P%o3p_gfld%LAT(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                       data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%LAT, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%LAT, nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  Longitude  ---------------
      IF (ASSOCIATED(L2O3P%o3p_gfld%LON)) DEALLOCATE(L2O3P%o3p_gfld%LON)
      data_name = O3Pds%Longitude_DataName
      ALLOCATE (L2O3P%o3p_gfld%LON(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                       data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%LON, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%LON, nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  Line       ---------------
      IF (ASSOCIATED(L2O3P%o3p_gfld%Line_tmp)) DEALLOCATE(L2O3P%o3p_gfld%Line_tmp)
      data_name = O3Pds%Line_DataName
      ALLOCATE (L2O3P%o3p_gfld%Line_tmp(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                       data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%Line_tmp, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%Line_tmp, nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  Pix        ---------------
      IF (ASSOCIATED(L2O3P%o3p_gfld%Pix_tmp)) DEALLOCATE(L2O3P%o3p_gfld%Pix_tmp)
      data_name = O3Pds%Pix_DataName
      ALLOCATE (L2O3P%o3p_gfld%Pix_tmp(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                       data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%Pix_tmp, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%Pix_tmp, nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SolarZenithAngle  ---------------
      IF (ASSOCIATED(L2O3P%o3p_gfld%SolZenAng)) DEALLOCATE(L2O3P%o3p_gfld%SolZenAng)
      data_name = O3Pds%SolarZenAng_DataName
      ALLOCATE (L2O3P%o3p_gfld%SolZenAng(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                       data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%SolZenAng, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%SolZenAng, nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  ViewingZenithAngle  ---------------
      IF (ASSOCIATED(L2O3P%o3p_gfld%ViewZenAng)) DEALLOCATE(L2O3P%o3p_gfld%ViewZenAng)
      data_name = O3Pds%ViewZenAng_DataName
      ALLOCATE (L2O3P%o3p_gfld%ViewZenAng(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                       data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%ViewZenAng, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%ViewZenAng, nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  RelativeAzimuthAngle  ---------------
      IF (ASSOCIATED(L2O3P%o3p_gfld%RelAziAng)) DEALLOCATE(L2O3P%o3p_gfld%RelAziAng)
      data_name = O3Pds%RelAziAng_DataName
      ALLOCATE (L2O3P%o3p_gfld%RelAziAng(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                       data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%RelAziAng, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%RelAziAng, nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  Time  ---------------
      IF (ASSOCIATED(L2O3P%o3p_gfld%Time)) DEALLOCATE(L2O3P%o3p_gfld%Time)
      data_name = O3Pds%Time_DataName
      ALLOCATE (L2O3P%o3p_gfld%Time(nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(O3Pds%l2_o3p_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                       data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%Time, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%Time, nline, data_name)
      ELSE
            GOTO 2999
      ENDIF

      write(disp_msg,'(" GEMS HDF5 File READ OK, hdferr=",i6)') hdferr
      !CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)


!--------------------------------
!---- END FUNCTION -------------
!--------------------------------
      !CALL GEMS_Share_MOD_log(llvl, '---   END L2 READ FUNCTION ---', ' END MSG')
      !CALL GEMS_Share_MOD_log(llvl, '.', LOGMSG)
      status = 0

      RETURN

999   CONTINUE
      write(disp_msg,'( A20, "Memory Allocation Failed, ierr=",i6)') ctl_msg, ierr
      CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)
      status = ierr

      RETURN

2999  CONTINUE
      write(disp_msg,'(" !!! GEMS HDF5 File READ ERROR !!! ", "[", A30, "], hdferr=",i6)') data_name, hdferr
      CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)


      RETURN

      END FUNCTION GEMS_Share_MOD_L2O3P_Read



!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Read NameList function   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE read_namelist_o3p(ctl_fpath, O3Pds)

! input/output parameters:
      CHARACTER(LEN=128),   INTENT(INOUT) :: ctl_fpath
      TYPE(O3P_ds),         INTENT(INOUT) :: O3Pds

! Local parameters For L2 Reading:
      CHARACTER*200                 :: l2_o3p_file_path
      CHARACTER*200                 :: o3p_dgrp_path
      CHARACTER*200                 :: o3p_ggrp_path

      CHARACTER*200                 :: AvgKernel_DataName
      CHARACTER*200                 :: O3_DataName
      CHARACTER*200                 :: O3Apriori_DataName
      CHARACTER*200                 :: O3AprioriErr_DataName
      CHARACTER*200                 :: O3RndmNsErr_DataName
      CHARACTER*200                 :: O3SolutionErr_DataName
      CHARACTER*200                 :: ColAmountO3_DataName
      CHARACTER*200                 :: DegOfFreedom_DataName
      CHARACTER*200                 :: CldPressure_DataName
      CHARACTER*200                 :: EffCldFracUV_DataName
      CHARACTER*200                 :: TerrRefltUV_DataName
      CHARACTER*200                 :: ProcessQFlag_DataName
      CHARACTER*200                 :: NumOfIteration_DataName
      CHARACTER*200                 :: ResidualsOfFit_DataName
      CHARACTER*200                 :: RMSErrOfFit_DataName
!      -- NEW GEO --
      CHARACTER*200                 :: TrppsPressure_DataName
      CHARACTER*200                 :: Pressure_DataName
      CHARACTER*200                 :: Altitude_DataName
      CHARACTER*200                 :: Temperature_DataName
      CHARACTER*200                 :: TerrPressure_DataName
      CHARACTER*200                 :: Latitude_DataName
      CHARACTER*200                 :: Longitude_DataName
      CHARACTER*200                 :: Line_DataName         ! add 2016.09.23
      CHARACTER*200                 :: Pix_DataName          ! add 2016.09.23
      CHARACTER*200                 :: SolarZenAng_DataName
      CHARACTER*200                 :: ViewZenAng_DataName
      CHARACTER*200                 :: RelAziAng_DataName
      CHARACTER*200                 :: Time_DataName

      INTEGER                       :: file_path_sz
      INTEGER                       :: grp_path_sz
      INTEGER                       :: data_name_sz

      Namelist /L2_O3P_File_List/l2_o3p_file_path, &
                                 o3p_dgrp_path,    &
                                 o3p_ggrp_path
      Namelist /L2_File_List_Value_Size/file_path_sz, &
                                        grp_path_sz,  &
                                        data_name_sz
      Namelist /L2_O3P_DATA_List/AvgKernel_DataName,     &
                                 O3_DataName,            &
                                 O3Apriori_DataName,     &
                                 O3AprioriErr_DataName,  &
                                 O3RndmNsErr_DataName,   &
                                 O3SolutionErr_DataName, &
                                 ColAmountO3_DataName,   &
                                 DegOfFreedom_DataName,  &
                                 CldPressure_DataName,   &
                                 EffCldFracUV_DataName,  &
                                 TerrRefltUV_DataName,   &
                                 ProcessQFlag_DataName,  &
                                 NumOfIteration_DataName,&
                                 ResidualsOfFit_DataName,&
                                 RMSErrOfFit_DataName,   &
                                 TrppsPressure_DataName, &
                                 Pressure_DataName,      &
                                 Altitude_DataName,      &
                                 Temperature_DataName,   &
                                 TerrPressure_DataName,  &
                                 Latitude_DataName,      &
                                 Longitude_DataName,     &
                                 Line_DataName,          &
                                 Pix_DataName,           &
                                 SolarZenAng_DataName,   &
                                 ViewZenAng_DataName,    &
                                 RelAziAng_DataName,     &
                                 Time_DataName



!--------------------------------
!---     Read Control file    ---
!--------------------------------
      OPEN(10, file=trim(ctl_fpath), status='old')
      READ(10, L2_O3P_File_List)
      READ(10, L2_O3P_DATA_List)
      READ(10, L2_File_List_Value_Size)
      CLOSE(10)

      O3Pds%AvgKernel_DataName     = AvgKernel_DataName
      O3Pds%O3_DataName            = O3_DataName
      O3Pds%O3Apriori_DataName     = O3Apriori_DataName
      O3Pds%O3AprioriErr_DataName  = O3AprioriErr_DataName
      O3Pds%O3RndmNsErr_DataName   = O3RndmNsErr_DataName
      O3Pds%O3SolutionErr_DataName = O3SolutionErr_DataName
      O3Pds%ColAmountO3_DataName   = ColAmountO3_DataName
      O3Pds%DegOfFreedom_DataName  = DegOfFreedom_DataName
      O3Pds%CldPressure_DataName   = CldPressure_DataName
      O3Pds%EffCldFracUV_DataName  = EffCldFracUV_DataName
      O3Pds%TerrRefltUV_DataName   = TerrRefltUV_DataName
      O3Pds%ProcessQFlag_DataName  = ProcessQFlag_DataName
      O3Pds%NumOfIteration_DataName= NumOfIteration_DataName
      O3Pds%ResidualsOfFit_DataName= ResidualsOfFit_DataName
      O3Pds%RMSErrOfFit_DataName   = RMSErrOfFit_DataName
      O3Pds%TrppsPressure_DataName = TrppsPressure_DataName
      O3Pds%Pressure_DataName      = Pressure_DataName
      O3Pds%Altitude_DataName      = Altitude_DataName
      O3Pds%Temperature_DataName   = Temperature_DataName
      O3Pds%TerrPressure_DataName  = TerrPressure_DataName
      O3Pds%Latitude_DataName      = Latitude_DataName
      O3Pds%Longitude_DataName     = Longitude_DataName
      O3Pds%Line_DataName          = Line_DataName
      O3Pds%Pix_DataName           = Pix_DataName
      O3Pds%SolarZenAng_DataName   = SolarZenAng_DataName
      O3Pds%ViewZenAng_DataName    = ViewZenAng_DataName
      O3Pds%RelAziAng_DataName     = RelAziAng_DataName
      O3Pds%Time_DataName          = Time_DataName
      
      O3Pds%l2_o3p_file_path          = l2_o3p_file_path 
      O3Pds%o3p_dgrp_path             = o3p_dgrp_path
      O3Pds%o3p_ggrp_path             = o3p_ggrp_path
      O3Pds%file_path_sz              = file_path_sz
      O3Pds%grp_path_sz               = grp_path_sz
      O3Pds%data_name_sz              = data_name_sz
      END SUBROUTINE read_namelist_o3p

SUBROUTINE GEMS_Share_Init_L2O3P_GlobalConstants
      ! Size for L2O3P dataset
      nline            =    gci_nline_l2
      nxtrack2_o3p     =    gci_nxtrack2_l2
      nlayer_o3p       =    gci_nlayer
      nlevel_o3p       =    gci_nlayer + 1
      nwavel3_o3p      =    gci_nwavel3_l2
      nResiduals_o3p   =    gci_nResiduals
      nColumns_o3p     =    gci_nColumns
END SUBROUTINE GEMS_Share_Init_L2O3P_GlobalConstants

SUBROUTINE GEMS_Share_Print_L2O3P_GlobalConstants
      ! Size for L2O3P dataset
    WRITE(*, '( A15 I5)') "nline         =",nline   
    WRITE(*, '( A15 I5)') "nxtrack2_o3p  =",nxtrack2_o3p
    WRITE(*, '( A15 I5)') "nlayer_o3p    =",nlayer_o3p  
    WRITE(*, '( A15 I5)') "nResiduals_o3p=",nResiduals_o3p
    WRITE(*, '( A15 I5)') "nColumns_o3p  =",nColumns_o3p
    WRITE(*, '( A15 I5)') "nlayer_o3p    =",nlayer_o3p 
    WRITE(*, '( A15 I5)') "nwavel3_o3p   =",nwavel3_o3p 
     
END SUBROUTINE GEMS_Share_Print_L2O3P_GlobalConstants

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!L2 O3P Memory Allocation function   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Declarations:
      FUNCTION GEMS_Share_MOD_L2O3P_MemAlloc() RESULT (status)

! input/output parameters:

! Local parameters For L2 Reading:
      INTEGER                       :: ierr        ! Error code 
      INTEGER (KIND = 4)            :: status

! Local parameters For Logging:


!--------------------------------
!---- FUNCTION Start -------------
!--------------------------------
      llvl  =9
      LOGMSG=" PROC MSG"
      !CALL GEMS_Share_MOD_log(llvl, '--- Start L2 O3P Memory Allocation FUNCTION ---', 'STARTMSG')

      !--------------------------------------
      ! <O3Profile: Data Fields>
      !--------------------------------------

      !-----------------  AveragingKernel  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_dfld%AvgK)) DEALLOCATE(gds_L2O3P%o3p_dfld%AvgK)
      ALLOCATE (gds_L2O3P%o3p_dfld%AvgK(nlayer_o3p,nlayer_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="AveragingKernel"
              GOTO 999
           ENDIF

      !-----------------  O3  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_dfld%OzRet)) DEALLOCATE(gds_L2O3P%o3p_dfld%OzRet)
      ALLOCATE (gds_L2O3P%o3p_dfld%OzRet(nlayer_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="O3"
              GOTO 999
           ENDIF

      !-----------------  O3Apriori  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_dfld%OzAp)) DEALLOCATE(gds_L2O3P%o3p_dfld%OzAp)
      ALLOCATE (gds_L2O3P%o3p_dfld%OzAp(nlayer_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="O3Apriori"
              GOTO 999
           ENDIF

      !-----------------  O3AprioriError  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_dfld%OzApErr)) DEALLOCATE(gds_L2O3P%o3p_dfld%OzApErr)
      ALLOCATE (gds_L2O3P%o3p_dfld%OzApErr(nlayer_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="O3AprioriError"
              GOTO 999
           ENDIF

      !-----------------  O3 Random Noise Error  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_dfld%OzNErr)) DEALLOCATE(gds_L2O3P%o3p_dfld%OzNErr)
      ALLOCATE (gds_L2O3P%o3p_dfld%OzNErr(nlayer_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="O3RandomNoiseError"
              GOTO 999
           ENDIF

      !-----------------  O3 Solution Error  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_dfld%OzSMErr)) DEALLOCATE(gds_L2O3P%o3p_dfld%OzSMErr)
      ALLOCATE (gds_L2O3P%o3p_dfld%OzSMErr(nlayer_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="O3SolutionError"
              GOTO 999
           ENDIF

      !-----------------  ColumnAmountO3  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_dfld%ColAmtOz)) DEALLOCATE(gds_L2O3P%o3p_dfld%ColAmtOz)
      ALLOCATE (gds_L2O3P%o3p_dfld%ColAmtOz(nColumns_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="ColumnAmountO3"
              GOTO 999
           ENDIF

      !-----------------  DegreesOfFreedomForSignal  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_dfld%DFS)) DEALLOCATE(gds_L2O3P%o3p_dfld%DFS)
      ALLOCATE (gds_L2O3P%o3p_dfld%DFS(nColumns_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="DegreesOfFreedomForSignal"
              GOTO 999
           ENDIF

      !-----------------  CloudPressure  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_dfld%CldP)) DEALLOCATE(gds_L2O3P%o3p_dfld%CldP)
      ALLOCATE (gds_L2O3P%o3p_dfld%CldP(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="CloudPressure"
              GOTO 999
           ENDIF

      !-----------------  EffectiveCloudFractionUV  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_dfld%EffCldFrac)) DEALLOCATE(gds_L2O3P%o3p_dfld%EffCldFrac)
      ALLOCATE (gds_L2O3P%o3p_dfld%EffCldFrac(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="EffectiveCloudFractionUV"
              GOTO 999
           ENDIF

      !-----------------  TerrainReflectivityUV  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_dfld%SfcAlb)) DEALLOCATE(gds_L2O3P%o3p_dfld%SfcAlb)
      ALLOCATE (gds_L2O3P%o3p_dfld%SfcAlb(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="TerrainReflectivityUV"
              GOTO 999
           ENDIF

      !-----------------  ProcessingQualityFlags  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_dfld%PrsQFlag)) DEALLOCATE(gds_L2O3P%o3p_dfld%PrsQFlag)
      ALLOCATE (gds_L2O3P%o3p_dfld%PrsQFlag(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="ProcessingQualityFlags"
              GOTO 999
           ENDIF

      !-----------------  NumberOfIterations  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_dfld%NumIter)) DEALLOCATE(gds_L2O3P%o3p_dfld%NumIter)
      ALLOCATE (gds_L2O3P%o3p_dfld%NumIter(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="NumberOfIterations"
              GOTO 999
           ENDIF

      !-----------------  ResidualsOfFit  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_dfld%AvgRes)) DEALLOCATE(gds_L2O3P%o3p_dfld%AvgRes)
      ALLOCATE (gds_L2O3P%o3p_dfld%AvgRes(nResiduals_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="ResidualsOfFit"
              GOTO 999
           ENDIF

      !-----------------  RootMeanSquareErrorOfFit  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_dfld%RMS)) DEALLOCATE(gds_L2O3P%o3p_dfld%RMS)
      ALLOCATE (gds_L2O3P%o3p_dfld%RMS(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="RootMeanSquareErrorOfFit"
              GOTO 999
           ENDIF

      !--------------------------------------
      ! <O3Profile: Geolocation Fields>
      !--------------------------------------
      !CALL GEMS_Share_MOD_log(llvl, "O3Profile: Geolocation Fields", LOGMSG)

      !-----------------  TropopausePressure  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_gfld%TrppsP)) DEALLOCATE(gds_L2O3P%o3p_gfld%TrppsP)
      ALLOCATE (gds_L2O3P%o3p_gfld%TrppsP(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="TropopausePressure"
              GOTO 999
           ENDIF

      !-----------------  Pressure  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_gfld%P)) DEALLOCATE(gds_L2O3P%o3p_gfld%P)
      ALLOCATE (gds_L2O3P%o3p_gfld%P(nlevel_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="Pressure"
              GOTO 999
           ENDIF

      !-----------------  Altitude  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_gfld%Alt)) DEALLOCATE(gds_L2O3P%o3p_gfld%Alt)
      ALLOCATE (gds_L2O3P%o3p_gfld%Alt(nlevel_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="Altitude"
              GOTO 999
           ENDIF

      !-----------------  Temperature  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_gfld%TEMP)) DEALLOCATE(gds_L2O3P%o3p_gfld%TEMP)
      ALLOCATE (gds_L2O3P%o3p_gfld%TEMP(nlevel_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="Temperature"
              GOTO 999
           ENDIF

      !-----------------  TerrainPressure  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_gfld%TerrP)) DEALLOCATE(gds_L2O3P%o3p_gfld%TerrP)
      ALLOCATE (gds_L2O3P%o3p_gfld%TerrP(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="TerrainPressure"
              GOTO 999
           ENDIF

      !-----------------  Latitude  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_gfld%LAT)) DEALLOCATE(gds_L2O3P%o3p_gfld%LAT)
      ALLOCATE (gds_L2O3P%o3p_gfld%LAT(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="Latitude"
              GOTO 999
           ENDIF

      !-----------------  Longitude  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_gfld%LON)) DEALLOCATE(gds_L2O3P%o3p_gfld%LON)
      ALLOCATE (gds_L2O3P%o3p_gfld%LON(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="Longitude"
              GOTO 999
           ENDIF

      !-----------------  Line       ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_gfld%Line_tmp)) DEALLOCATE(gds_L2O3P%o3p_gfld%Line_tmp)
      ALLOCATE (gds_L2O3P%o3p_gfld%Line_tmp(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="Line"
              GOTO 999
           ENDIF

      !-----------------  Pix        ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_gfld%Pix_tmp)) DEALLOCATE(gds_L2O3P%o3p_gfld%Pix_tmp)
      ALLOCATE (gds_L2O3P%o3p_gfld%Pix_tmp(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="Pix"
              GOTO 999
           ENDIF

      !-----------------  SolarZenithAngle  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_gfld%SolZenAng)) DEALLOCATE(gds_L2O3P%o3p_gfld%SolZenAng)
      ALLOCATE (gds_L2O3P%o3p_gfld%SolZenAng(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="SolarZenithAngle"
              GOTO 999
           ENDIF

      !-----------------  ViewingZenithAngle  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_gfld%ViewZenAng)) DEALLOCATE(gds_L2O3P%o3p_gfld%ViewZenAng)
      ALLOCATE (gds_L2O3P%o3p_gfld%ViewZenAng(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="ViewingZenithAngle"
              GOTO 999
           ENDIF

      !-----------------  RelativeAzimuthAngle  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_gfld%RelAziAng)) DEALLOCATE(gds_L2O3P%o3p_gfld%RelAziAng)
      ALLOCATE (gds_L2O3P%o3p_gfld%RelAziAng(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="RelativeAzimuthAngle"
              GOTO 999
           ENDIF

      !-----------------  Time  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_gfld%Time)) DEALLOCATE(gds_L2O3P%o3p_gfld%Time)
      ALLOCATE (gds_L2O3P%o3p_gfld%Time(nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="Time"
              GOTO 999
           ENDIF

      gds_L2O3P%status = 1
      write(disp_msg,'("L2 O3P Memory Allocation OK !")')
      !CALL GEMS_Share_MOD_log(0, disp_msg, LOGMSG)

!--------------------------------
!---- END FUNCTION -------------
!--------------------------------
      !CALL GEMS_Share_MOD_log(llvl, '---   END L2 O3P Memory Allocation FUNCTION ---', LOGMSG)
      !CALL GEMS_Share_MOD_log(llvl, '.', LOGMSG)
      status = 0

      RETURN

999   CONTINUE
      llvl  =9
      write(disp_msg,'( A30, "L2 O3P Memory Allocation Failed, ierr=",i6)') ctl_msg,ierr
      CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)
      status = ierr

      RETURN

      END FUNCTION GEMS_Share_MOD_L2O3P_MemAlloc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!L2 O3P Memory DeAllocation function   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Declarations:
      FUNCTION GEMS_Share_MOD_L2O3P_MemDeAlloc() RESULT (status)

! input/output parameters:

! Local parameters For L2 Reading:
      INTEGER (KIND = 4)            :: status

! Local parameters For Logging:


!--------------------------------
!---- FUNCTION Start -------------
!--------------------------------
      llvl=9
      LOGMSG=" PROC MSG"
      !CALL GEMS_Share_MOD_log(llvl, '--- Start L2 O3P Memory DeAllocation FUNCTION  ---', LOGMSG)

      !--------------------------------------
      ! <O3Profile: Data Fields>
      !--------------------------------------
      !-----------------  AveragingKernel  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_dfld%AvgK)) DEALLOCATE(gds_L2O3P%o3p_dfld%AvgK)

      !-----------------  O3  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_dfld%OzRet)) DEALLOCATE(gds_L2O3P%o3p_dfld%OzRet)

      !-----------------  O3Apriori  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_dfld%OzAp)) DEALLOCATE(gds_L2O3P%o3p_dfld%OzAp)

      !-----------------  O3AprioriError  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_dfld%OzApErr)) DEALLOCATE(gds_L2O3P%o3p_dfld%OzApErr)

      !-----------------  O3 Random Noise Error  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_dfld%OzNErr)) DEALLOCATE(gds_L2O3P%o3p_dfld%OzNErr)

      !-----------------  O3 Solution Error  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_dfld%OzSMErr)) DEALLOCATE(gds_L2O3P%o3p_dfld%OzSMErr)

      !-----------------  ColumnAmountO3  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_dfld%ColAmtOz)) DEALLOCATE(gds_L2O3P%o3p_dfld%ColAmtOz)

      !-----------------  DegreesOfFreedomForSignal  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_dfld%DFS)) DEALLOCATE(gds_L2O3P%o3p_dfld%DFS)

      !-----------------  CloudPressure  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_dfld%CldP)) DEALLOCATE(gds_L2O3P%o3p_dfld%CldP)

      !-----------------  EffectiveCloudFractionUV  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_dfld%EffCldFrac)) DEALLOCATE(gds_L2O3P%o3p_dfld%EffCldFrac)

      !-----------------  TerrainReflectivityUV  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_dfld%SfcAlb)) DEALLOCATE(gds_L2O3P%o3p_dfld%SfcAlb)

      !-----------------  ProcessingQualityFlags  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_dfld%PrsQFlag)) DEALLOCATE(gds_L2O3P%o3p_dfld%PrsQFlag)

      !-----------------  NumberOfIterations  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_dfld%NumIter)) DEALLOCATE(gds_L2O3P%o3p_dfld%NumIter)

      !-----------------  ResidualsOfFit  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_dfld%AvgRes)) DEALLOCATE(gds_L2O3P%o3p_dfld%AvgRes)

      !-----------------  RootMeanSquareErrorOfFit  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_dfld%RMS)) DEALLOCATE(gds_L2O3P%o3p_dfld%RMS)

      !--------------------------------------
      ! <O3Profile: Geolocation Fields>
      !--------------------------------------

      !-----------------  TropopausePressure  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_gfld%TrppsP)) DEALLOCATE(gds_L2O3P%o3p_gfld%TrppsP)

      !-----------------  Pressure  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_gfld%P)) DEALLOCATE(gds_L2O3P%o3p_gfld%P)

      !-----------------  Altitude  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_gfld%Alt)) DEALLOCATE(gds_L2O3P%o3p_gfld%Alt)

      !-----------------  Temperature  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_gfld%TEMP)) DEALLOCATE(gds_L2O3P%o3p_gfld%TEMP)

      !-----------------  TerrainPressure  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_gfld%TerrP)) DEALLOCATE(gds_L2O3P%o3p_gfld%TerrP)

      !-----------------  Latitude  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_gfld%LAT)) DEALLOCATE(gds_L2O3P%o3p_gfld%LAT)

      !-----------------  Longitude  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_gfld%LON)) DEALLOCATE(gds_L2O3P%o3p_gfld%LON)

      !-----------------  Line       ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_gfld%Line_tmp)) DEALLOCATE(gds_L2O3P%o3p_gfld%Line_tmp)

      !-----------------  Pix        ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_gfld%Pix_tmp)) DEALLOCATE(gds_L2O3P%o3p_gfld%Pix_tmp)

      !-----------------  SolarZenithAngle  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_gfld%SolZenAng)) DEALLOCATE(gds_L2O3P%o3p_gfld%SolZenAng)

      !-----------------  ViewingZenithAngle  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_gfld%ViewZenAng)) DEALLOCATE(gds_L2O3P%o3p_gfld%ViewZenAng)

      !-----------------  RelativeAzimuthAngle  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_gfld%RelAziAng)) DEALLOCATE(gds_L2O3P%o3p_gfld%RelAziAng)

      !-----------------  Time  ---------------
      IF (ASSOCIATED(gds_L2O3P%o3p_gfld%Time)) DEALLOCATE(gds_L2O3P%o3p_gfld%Time)

      gds_L2O3P%status = 0
      write(disp_msg,'(" L2 O3P Memory DeAllocation OK !")')
      !CALL GEMS_Share_MOD_log(0, disp_msg, LOGMSG)

!--------------------------------
!---- END FUNCTION -------------
!--------------------------------
      !CALL GEMS_Share_MOD_log(llvl, '---   END L2 O3P Memory DeAllocation FUNCTION ---', LOGMSG)
      !CALL GEMS_Share_MOD_log(llvl, '.', LOGMSG)
      status = 0

      RETURN

      END FUNCTION GEMS_Share_MOD_L2O3P_MemDeAlloc


!!!!!!!!!!!!!!!!!!!!!!!!!!!
!L2 O3P write function2   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION GEMS_Share_MOD_L2O3P_Write2(L2O3P, ctl_fpath, wr_file_path) RESULT (status)

! input/output parameters:
      CHARACTER(LEN=128),   INTENT(INOUT) :: ctl_fpath
      CHARACTER(LEN=128),   INTENT(IN )   :: wr_file_path
      TYPE(L2_o3p),         INTENT(INOUT) :: L2O3P

! Local parameters For L2 Reading:
      TYPE(O3P_ds)                  :: O3Pds

      CHARACTER*200                 :: data_name
      CHARACTER*200                 :: grp_path

      INTEGER                       :: pstart(3)
      INTEGER                       :: pedge(3)

      INTEGER                       :: ierr        ! Error code 
      INTEGER                       :: hdferr      ! HDF5 Read Error code 
      INTEGER*8                     :: iAddr       ! memory address#

      INTEGER (KIND = 4)            :: status


! Local parameters For Logging:

      INTEGER       :: i, j, k

      INTEGER(1), PARAMETER :: X_FF = INT(Z'FF',1)


!--------------------------------
      !
      ! External Function Declaration
      !
      INTEGER  GEMS_Share_Hdf5ReadData
      INTEGER  GEMS_Share_Hdf5WriteData
      INTEGER  GEMS_Share_Hdf5ReadAttr
      INTEGER  GEMS_Share_Hdf5WriteAttr
      INTEGER  GEMS_Share_Hdf5CreateL2File
      INTEGER  GEMS_Share_Hdf5CreateL1BFile
      external GEMS_Share_Hdf5ReadData
      external GEMS_Share_Hdf5WriteData
      external GEMS_Share_Hdf5ReadAttr
      external GEMS_Share_Hdf5WriteAttr
      external GEMS_Share_Hdf5CreateL2File
      external GEMS_Share_Hdf5CreateL1BFile
!--------------------------------


!--------------------------------
!---- FUNCTION Start -------------
!--------------------------------
      LOGMSG=" PROC MSG"
      !CALL GEMS_Share_MOD_log(llvl, '--- Start L2O3P WRITE FUNCTION ver2  ---', 'START MSG')

!
!--------------------------------

     pstart(1) = -1                                            ! HDF5 데이터셑의 모든 데이터값을 읽는다는 의미

#ifdef DEBUG
      !WRITE(disp_msg, '( A25, A128, A1 )') "NameList file path     =[",trim(ctl_fpath),  "]"
      !CALL GEMS_Share_MOD_log(llvl, disp_msg,  LOGMSG)
#endif
      CALL read_namelist_o3p(ctl_fpath, O3Pds)
!      CALL GEMS_Share_Init_L2O3P_GlobalConstants


      !----------      L2O3P     ----------!
      !--------------------------------------
      ! <O3 Profile: Data Fields>
      !--------------------------------------

      llvl  =0
      !disp_msg = "L2O3P WRITE FILE PATH="//wr_file_path
      !CALL GEMS_Share_MOD_log(llvl, disp_msg,  LOGMSG)

      !CALL GEMS_Share_MOD_log(llvl, "O3 Profile :Data Fields",  LOGMSG)

      !-----------------  AveragingKernel  ---------------
      data_name = O3Pds%AvgKernel_DataName
      grp_path  = O3Pds%o3p_dgrp_path
#ifdef DEBUG
      !WRITE(disp_msg, '( A25, A128, A1 )') "L2 O3P write file path =[",trim(wr_file_path)           , "]"
      !CALL GEMS_Share_MOD_log(llvl, disp_msg,  LOGMSG)
      !WRITE(disp_msg, '( A25, I5,   A1 )') "file path size         =[",O3Pds%file_path_sz          , "]"
      !CALL GEMS_Share_MOD_log(llvl, disp_msg,  LOGMSG)
      !WRITE(disp_msg, '( A25, A128, A1 )') "data group path        =[",trim(grp_path)  , "]"
      !CALL GEMS_Share_MOD_log(llvl, disp_msg,  LOGMSG)
      !WRITE(disp_msg, '( A25, I5,   A1 )') "data group path size   =[",O3Pds%grp_path_sz           , "]"
      !CALL GEMS_Share_MOD_log(llvl, disp_msg,  LOGMSG)
      !WRITE(disp_msg, '( A25, A128, A1 )') "dataset name           =[",trim(data_name)              , "]"
      !CALL GEMS_Share_MOD_log(llvl, disp_msg,  LOGMSG)
      !WRITE(disp_msg, '( A25, I5,   A1 )') "data name size         =[",O3Pds%data_name_sz          , "]"
      !CALL GEMS_Share_MOD_log(llvl, disp_msg,  LOGMSG)
      !WRITE(disp_msg, '( A31, I12,   A1 )') "AveragingKernel size =[",SIZE(L2O3P%o3p_dfld%AvgK)  , "]"
      !CALL GEMS_Share_MOD_log(llvl, disp_msg,  LOGMSG)
#endif

      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%AvgK, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%AvgK, nline,nxtrack2_o3p,nlayer_o3p,nlayer_o3p, data_name)
      ELSE
            GOTO 999
      ENDIF

      !-----------------  O3  ---------------
      data_name = O3Pds%O3_DataName
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%OzRet, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%OzRet, nline,nxtrack2_o3p,nlayer_o3p, data_name)
      ELSE
            GOTO 999
      ENDIF

      !-----------------  O3Apriori  ---------------
      data_name = O3Pds%O3Apriori_DataName
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%OzAp, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%OzAp, nline,nxtrack2_o3p,nlayer_o3p, data_name)
      ELSE
            GOTO 999
      ENDIF

      !-----------------  O3AprioriError  ---------------
      data_name = O3Pds%O3AprioriErr_DataName
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%OzApErr, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%OzApErr, nline,nxtrack2_o3p,nlayer_o3p, data_name)
      ELSE
            GOTO 999
      ENDIF

      !-----------------  O3 Random Noise Error  ---------------
      data_name = O3Pds%O3RndmNsErr_DataName
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%OzNErr, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%OzNErr, nline,nxtrack2_o3p,nlayer_o3p, data_name)
      ELSE
            GOTO 999
      ENDIF

      !-----------------  O3 Solution Error  ---------------
      data_name = O3Pds%O3SolutionErr_DataName
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%OzSMErr, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%OzSMErr,nline,nxtrack2_o3p,nlayer_o3p, data_name)
      ELSE
            GOTO 999
      ENDIF

      !-----------------  ColumnAmountO3  ---------------
      data_name = O3Pds%ColAmountO3_DataName
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%ColAmtOz, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%ColAmtOz, nline,nxtrack2_o3p,nColumns_o3p, data_name)
      ELSE
            GOTO 999
      ENDIF

      !-----------------  DegreesOfFreedomForSignal  ---------------
      data_name = O3Pds%DegOfFreedom_DataName
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%DFS, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%DFS, nline,nxtrack2_o3p,nColumns_o3p, data_name)
      ELSE
            GOTO 999
      ENDIF

      !-----------------  CloudPressure  ---------------
      data_name = O3Pds%CldPressure_DataName
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%CldP, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%CldP, nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 999
      ENDIF

      !-----------------  EffectiveCloudFractionUV  ---------------
      data_name = O3Pds%EffCldFracUV_DataName
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%EffCldFrac, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%EffCldFrac, nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 999
      ENDIF

      !-----------------  TerrainReflectivityUV  ---------------
      data_name = O3Pds%TerrRefltUV_DataName
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%SfcAlb, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%SfcAlb,nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 999
      ENDIF

      !-----------------  ProcessingQualityFlags  ---------------
      data_name = O3Pds%ProcessQFlag_DataName
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%PrsQFlag, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%PrsQFlag,nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 999
      ENDIF

      !-----------------  NumberOfIterations  ---------------
      data_name = O3Pds%NumOfIteration_DataName      
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%NumIter, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%NumIter, nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 999
      ENDIF

      !-----------------  ResidualsOfFit  ---------------
      data_name = O3Pds%ResidualsOfFit_DataName
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%AvgRes, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%AvgRes,nline,nxtrack2_o3p,nResiduals_o3p, data_name)
      ELSE
            GOTO 999
      ENDIF

      !-----------------  RootMeanSquareErrorOfFit  ---------------
      data_name = O3Pds%RMSErrOfFit_DataName
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_dfld%RMS, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_dfld%RMS,nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 999
      ENDIF

      !--------------------------------------
      ! <O3Profile: Geolocation Fields>
      !--------------------------------------
      !CALL GEMS_Share_MOD_log(llvl, "O3Profile: Geolocation Fields", LOGMSG)

      !-----------------  TropopausePressure  ---------------
      data_name = O3Pds%TrppsPressure_DataName
      grp_path  = O3Pds%o3p_ggrp_path
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%TrppsP, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%TrppsP, nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 999
      ENDIF

      !-----------------  Pressure  ---------------
      data_name = O3Pds%Pressure_DataName
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%P, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%P, nline,nxtrack2_o3p,nlevel_o3p, data_name)
      ELSE
            GOTO 999
      ENDIF

      !-----------------  Altitude  ---------------
      data_name = O3Pds%Altitude_DataName
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%Alt, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%Alt, nline,nxtrack2_o3p,nlevel_o3p, data_name)
      ELSE
            GOTO 999
      ENDIF

      !-----------------  Temperature  ---------------
      data_name = O3Pds%Temperature_DataName
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%TEMP, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%TEMP,nline,nxtrack2_o3p,nlevel_o3p, data_name)
      ELSE
            GOTO 999
      ENDIF

      !-----------------  TerrainPressure  ---------------
      data_name = O3Pds%TerrPressure_DataName
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%TerrP, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%TerrP,nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 999
      ENDIF

      !-----------------  Latitude  ---------------
      data_name = O3Pds%Latitude_DataName
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%LAT, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%LAT, nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 999
      ENDIF

      !-----------------  Longitude  ---------------
      data_name = O3Pds%Longitude_DataName
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%LON, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%LON, nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 999
      ENDIF

      !-----------------  Line       ---------------
      data_name = O3Pds%Line_DataName
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%Line_tmp, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%Line_tmp, nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 999
      ENDIF

      !-----------------  Pix        ---------------
      data_name = O3Pds%Pix_DataName
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%Pix_tmp, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%Pix_tmp, nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 999
      ENDIF

      !-----------------  SolarZenithAngle  ---------------
      data_name = O3Pds%SolarZenAng_DataName
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%SolZenAng, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%SolZenAng, nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 999
      ENDIF

      !-----------------  ViewingZenithAngle  ---------------
      data_name = O3Pds%ViewZenAng_DataName
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%ViewZenAng, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%ViewZenAng, nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 999
      ENDIF

      !-----------------  RelativeAzimuthAngle  ---------------
      data_name = O3Pds%RelAziAng_DataName
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%RelAziAng, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%RelAziAng, nline,nxtrack2_o3p, data_name)
      ELSE
            GOTO 999
      ENDIF

      !-----------------  Time  ---------------
      data_name = O3Pds%Time_DataName
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, O3Pds%file_path_sz, grp_path, O3Pds%grp_path_sz, &
                                        data_name, O3Pds%data_name_sz, L2O3P%o3p_gfld%Time, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2O3P%o3p_gfld%Time, nline, data_name)
      ELSE
            GOTO 999
      ENDIF

      write(disp_msg,'(" GEMS HDF5 File WRITE OK, hdferr=",i6)') hdferr
      !CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)

!--------------------------------
!---- END FUNCTION -------------
!--------------------------------
      llvl  =9
      !CALL GEMS_Share_MOD_log(llvl, '---   END L2 WRITE FUNCTION ver2 ---', '  END MSG')
      !CALL GEMS_Share_MOD_log(llvl, '.', LOGMSG)
      status = 0

      RETURN

999   CONTINUE
      write(disp_msg,'(" !!! GEMS HDF5 File WRITE ver2 ERROR !!! ", "[", A30, "],hdferr=",i6)') data_name, hdferr
      CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)
      status = hdferr

      RETURN

END FUNCTION GEMS_Share_MOD_L2O3P_Write2

!-------------------------------------------------------------------------------
!+Description: 
!      FUNCTION to initialize L2 O3P Memory.
!
!
! Method:
!
! input  :
!
! input/output :
!
! output :
!       status    : the result code of L2O3P Memory Initialization
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2015.04.03  Fisrt Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!L2 O3P Memory Allocation function  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION GEMS_Share_MOD_L2O3P_MemInit() RESULT (status)

! input/output parameters:

! Local parameters For L2 Reading:
      CHARACTER*200                 :: data_name
      INTEGER                       :: ierr        ! Error code 
      INTEGER (KIND = 4)            :: status

! Local parameters For Logging:
!--------------------------------
!---- FUNCTION Start -------------
!--------------------------------
      llvl  =9
      LOGMSG=" PROC MSG"
      !CALL GEMS_Share_MOD_log(llvl, '--- Start L2 O3P Memory Initialization FUNCTION  ---', LOGMSG)

      IF ( gds_L2O3P%status == 0 ) THEN
          llvl    =0
          disp_msg="did not allocate L2 O3P Memory"
          CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)
          status = -7001  ! MEMOEY ALLOCATION ERROR CODE = -7001
          RETURN
      END IF


      !--------------------------------------
      ! <O3Profile: Data Fields>
      !--------------------------------------
      !-----------------  AveragingKernel  ---------------
      gds_L2O3P%o3p_dfld%AvgK = DEF_fVAL

      !-----------------  O3  ---------------
      gds_L2O3P%o3p_dfld%OzRet = DEF_fVAL

      !-----------------  O3Apriori  ---------------
      gds_L2O3P%o3p_dfld%OzAp = DEF_fVAL

      !-----------------  O3AprioriError  ---------------
      gds_L2O3P%o3p_dfld%OzApErr = DEF_fVAL

      !-----------------  O3 Random Noise Error  ---------------
      gds_L2O3P%o3p_dfld%OzNErr = DEF_fVAL

      !-----------------  O3 Solution Error  ---------------
      gds_L2O3P%o3p_dfld%OzSMErr = DEF_fVAL

      !-----------------  ColumnAmountO3  ---------------
      gds_L2O3P%o3p_dfld%ColAmtOz = DEF_fVAL

      !-----------------  DegreesOfFreedomForSignal  ---------------
      gds_L2O3P%o3p_dfld%DFS = DEF_fVAL

      !-----------------  CloudPressure  ---------------
      gds_L2O3P%o3p_dfld%CldP = DEF_fVAL

      !-----------------  EffectiveCloudFractionUV  ---------------
      gds_L2O3P%o3p_dfld%EffCldFrac = DEF_fVAL

      !-----------------  TerrainReflectivityUV  ---------------
      gds_L2O3P%o3p_dfld%SfcAlb = DEF_sVAL

      !-----------------  ProcessingQualityFlags  ---------------
      gds_L2O3P%o3p_dfld%PrsQFlag = DEF_iVAL

      !-----------------  NumberOfIterations  ---------------
      gds_L2O3P%o3p_dfld%NumIter = DEF_sVAL

      !-----------------  ResidualsOfFit  ---------------
      gds_L2O3P%o3p_dfld%AvgRes = DEF_fVAL

      !-----------------  RootMeanSquareErrorOfFit  ---------------
      gds_L2O3P%o3p_dfld%RMS = DEF_fVAL


      !--------------------------------------
      ! <O3Profile: Geolocation Fields>
      !--------------------------------------

      !-----------------  TropopausePressure  ---------------
      gds_L2O3P%o3p_gfld%TrppsP = DEF_fVAL

      !-----------------  Pressure  ---------------
      gds_L2O3P%o3p_gfld%P = DEF_fVAL

      !-----------------  Altitude  ---------------
      gds_L2O3P%o3p_gfld%Alt = DEF_fVAL

      !-----------------  Temperature  ---------------
      gds_L2O3P%o3p_gfld%TEMP = DEF_fVAL

      !-----------------  TerrainPressure -------------
      gds_L2O3P%o3p_gfld%TerrP = DEF_fVAL

      !-----------------  Latitude  ---------------
      gds_L2O3P%o3p_gfld%LAT = DEF_fVAL

      !-----------------  Longitude  ---------------
      gds_L2O3P%o3p_gfld%LON = DEF_fVAL

      !-----------------  Line  ---------------
      gds_L2O3P%o3p_gfld%Line_tmp = DEF_iVAL

      !-----------------  Pix  ---------------
      gds_L2O3P%o3p_gfld%Pix_tmp = DEF_iVAL

      !-----------------  SolarZenithAngle  ---------------
      gds_L2O3P%o3p_gfld%SolZenAng = DEF_fVAL

      !-----------------  ViewingZenithAngle  ---------------
      gds_L2O3P%o3p_gfld%ViewZenAng = DEF_fVAL

      !-----------------  RelativeAzimuthAngle  ---------------
      gds_L2O3P%o3p_gfld%RelAziAng = DEF_fVAL

      !-----------------  Time  ---------------
      gds_L2O3P%o3p_gfld%Time = DEF_dVAL


      write(disp_msg,'("L2 O3P Memory Initialization OK !")')
      !CALL GEMS_Share_MOD_log(0, disp_msg, LOGMSG)

      status = 0

      RETURN

END FUNCTION GEMS_Share_MOD_L2O3P_MemInit

!-------------------------------------------------------------------------------
!+Description: 
!      FUNCTION ver2 to allocate L2 O3P Memory.
!
!
! Method:
!
! input  :
!
! input/output :
!     local_L2O3P : L2 O3P Data
!
! output :
!       status    : the result code of L2O3P Memory allocation
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2016.06.13  전역변수를 매개변수로 사용하기 위해 생성(YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------
! Declarations:
FUNCTION GEMS_Share_MOD_L2O3P_MemAlloc2(local_L2O3P) RESULT (status)

! input/output parameters:
      TYPE(L2_O3P)                  :: local_L2O3P

! Local parameters For L2 Reading:
      INTEGER                       :: ierr        ! Error code 
      INTEGER (KIND = 4)            :: status

! Local parameters For Logging:


!--------------------------------
!---- FUNCTION Start -------------
!--------------------------------
      llvl  =9
      LOGMSG=" PROC MSG"
      !CALL GEMS_Share_MOD_log(llvl, '--- Start L2 O3P Memory Allocation FUNCTION ver2 ---', 'STARTMSG')


      !--------------------------------------
      ! <O3Profile: Data Fields>
      !--------------------------------------

      !-----------------  AveragingKernel  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_dfld%AvgK)) DEALLOCATE(local_L2O3P%o3p_dfld%AvgK)
      ALLOCATE (local_L2O3P%o3p_dfld%AvgK(nlayer_o3p,nlayer_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="AveragingKernel"
              GOTO 999
           ENDIF

      !-----------------  O3  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_dfld%OzRet)) DEALLOCATE(local_L2O3P%o3p_dfld%OzRet)
      ALLOCATE (local_L2O3P%o3p_dfld%OzRet(nlayer_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="O3"
              GOTO 999
           ENDIF

      !-----------------  O3Apriori  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_dfld%OzAp)) DEALLOCATE(local_L2O3P%o3p_dfld%OzAp)
      ALLOCATE (local_L2O3P%o3p_dfld%OzAp(nlayer_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="O3Apriori"
              GOTO 999
           ENDIF

      !-----------------  O3AprioriError  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_dfld%OzApErr)) DEALLOCATE(local_L2O3P%o3p_dfld%OzApErr)
      ALLOCATE (local_L2O3P%o3p_dfld%OzApErr(nlayer_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="O3AprioriError"
              GOTO 999
           ENDIF

      !-----------------  O3 Random Noise Error  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_dfld%OzNErr)) DEALLOCATE(local_L2O3P%o3p_dfld%OzNErr)
      ALLOCATE (local_L2O3P%o3p_dfld%OzNErr(nlayer_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="O3RandomNoiseError"
              GOTO 999
           ENDIF

      !-----------------  O3 Solution Error  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_dfld%OzSMErr)) DEALLOCATE(local_L2O3P%o3p_dfld%OzSMErr)
      ALLOCATE (local_L2O3P%o3p_dfld%OzSMErr(nlayer_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="O3SolutionError"
              GOTO 999
           ENDIF

      !-----------------  ColumnAmountO3  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_dfld%ColAmtOz)) DEALLOCATE(local_L2O3P%o3p_dfld%ColAmtOz)
      ALLOCATE (local_L2O3P%o3p_dfld%ColAmtOz(nColumns_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="ColumnAmountO3"
              GOTO 999
           ENDIF

      !-----------------  DegreesOfFreedomForSignal  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_dfld%DFS)) DEALLOCATE(local_L2O3P%o3p_dfld%DFS)
      ALLOCATE (local_L2O3P%o3p_dfld%DFS(nColumns_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="DegreesOfFreedomForSignal"
              GOTO 999
           ENDIF

      !-----------------  CloudPressure  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_dfld%CldP)) DEALLOCATE(local_L2O3P%o3p_dfld%CldP)
      ALLOCATE (local_L2O3P%o3p_dfld%CldP(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="CloudPressure"
              GOTO 999
           ENDIF

      !-----------------  EffectiveCloudFractionUV  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_dfld%EffCldFrac)) DEALLOCATE(local_L2O3P%o3p_dfld%EffCldFrac)
      ALLOCATE (local_L2O3P%o3p_dfld%EffCldFrac(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="EffectiveCloudFractionUV"
              GOTO 999
           ENDIF

      !-----------------  TerrainReflectivityUV  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_dfld%SfcAlb)) DEALLOCATE(local_L2O3P%o3p_dfld%SfcAlb)
      ALLOCATE (local_L2O3P%o3p_dfld%SfcAlb(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="TerrainReflectivityUV"
              GOTO 999
           ENDIF

      !-----------------  ProcessingQualityFlags  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_dfld%PrsQFlag)) DEALLOCATE(local_L2O3P%o3p_dfld%PrsQFlag)
      ALLOCATE (local_L2O3P%o3p_dfld%PrsQFlag(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="ProcessingQualityFlags"
              GOTO 999
           ENDIF

      !-----------------  NumberOfIterations  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_dfld%NumIter)) DEALLOCATE(local_L2O3P%o3p_dfld%NumIter)
      ALLOCATE (local_L2O3P%o3p_dfld%NumIter(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="NumberOfIterations"
              GOTO 999
           ENDIF

      !-----------------  ResidualsOfFit  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_dfld%AvgRes)) DEALLOCATE(local_L2O3P%o3p_dfld%AvgRes)
      ALLOCATE (local_L2O3P%o3p_dfld%AvgRes(nResiduals_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="ResidualsOfFit"
              GOTO 999
           ENDIF

      !-----------------  RootMeanSquareErrorOfFit  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_dfld%RMS)) DEALLOCATE(local_L2O3P%o3p_dfld%RMS)
      ALLOCATE (local_L2O3P%o3p_dfld%RMS(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="RootMeanSquareErrorOfFit"
              GOTO 999
           ENDIF

      !--------------------------------------
      ! <O3Profile: Geolocation Fields>
      !--------------------------------------
      !CALL GEMS_Share_MOD_log(llvl, "O3Profile: Geolocation Fields", LOGMSG)


      !-----------------  TropopausePressure  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_gfld%TrppsP)) DEALLOCATE(local_L2O3P%o3p_gfld%TrppsP)
      ALLOCATE (local_L2O3P%o3p_gfld%TrppsP(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="Altitude"
              GOTO 999
           ENDIF

      !-----------------  Pressure  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_gfld%P)) DEALLOCATE(local_L2O3P%o3p_gfld%P)
      ALLOCATE (local_L2O3P%o3p_gfld%P(nlevel_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="Pressure"
              GOTO 999
           ENDIF

      !-----------------  Altitude  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_gfld%Alt)) DEALLOCATE(local_L2O3P%o3p_gfld%Alt)
      ALLOCATE (local_L2O3P%o3p_gfld%Alt(nlevel_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="Altitude"
              GOTO 999
           ENDIF

      !-----------------  Temperature  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_gfld%TEMP)) DEALLOCATE(local_L2O3P%o3p_gfld%TEMP)
      ALLOCATE (local_L2O3P%o3p_gfld%TEMP(nlevel_o3p,nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="Temperature"
              GOTO 999
           ENDIF

      !-----------------  TerrainPressure  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_gfld%TerrP)) DEALLOCATE(local_L2O3P%o3p_gfld%TerrP)
      ALLOCATE (local_L2O3P%o3p_gfld%TerrP(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="TerrainPressure"
              GOTO 999
           ENDIF

      !-----------------  Latitude  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_gfld%LAT)) DEALLOCATE(local_L2O3P%o3p_gfld%LAT)
      ALLOCATE (local_L2O3P%o3p_gfld%LAT(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="Latitude"
              GOTO 999
           ENDIF

      !-----------------  Longitude  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_gfld%LON)) DEALLOCATE(local_L2O3P%o3p_gfld%LON)
      ALLOCATE (local_L2O3P%o3p_gfld%LON(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="Longitude"
              GOTO 999
           ENDIF

      !-----------------  Line       ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_gfld%Line_tmp)) DEALLOCATE(local_L2O3P%o3p_gfld%Line_tmp)
      ALLOCATE (local_L2O3P%o3p_gfld%Line_tmp(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="Line"
              GOTO 999
           ENDIF

      !-----------------  Pix        ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_gfld%Pix_tmp)) DEALLOCATE(local_L2O3P%o3p_gfld%Pix_tmp)
      ALLOCATE (local_L2O3P%o3p_gfld%Pix_tmp(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="Pix"
              GOTO 999
           ENDIF

      !-----------------  SolarZenithAngle  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_gfld%SolZenAng)) DEALLOCATE(local_L2O3P%o3p_gfld%SolZenAng)
      ALLOCATE (local_L2O3P%o3p_gfld%SolZenAng(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="SolarZenithAngle"
              GOTO 999
           ENDIF

      !-----------------  ViewingZenithAngle  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_gfld%ViewZenAng)) DEALLOCATE(local_L2O3P%o3p_gfld%ViewZenAng)
      ALLOCATE (local_L2O3P%o3p_gfld%ViewZenAng(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="ViewingZenithAngle"
              GOTO 999
           ENDIF

      !-----------------  RelativeAzimuthAngle  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_gfld%RelAziAng)) DEALLOCATE(local_L2O3P%o3p_gfld%RelAziAng)
      ALLOCATE (local_L2O3P%o3p_gfld%RelAziAng(nxtrack2_o3p,nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="RelativeAzimuthAngle"
              GOTO 999
           ENDIF

      !-----------------  Time  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_gfld%Time)) DEALLOCATE(local_L2O3P%o3p_gfld%Time)
      ALLOCATE (local_L2O3P%o3p_gfld%Time(nline), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="Time"
              GOTO 999
           ENDIF

      local_L2O3P%status = 1
      !write(disp_msg,'("L2 O3P Memory Allocation ver2 OK !")')
      !CALL GEMS_Share_MOD_log(0, disp_msg, LOGMSG)

!--------------------------------
!---- END FUNCTION -------------
!--------------------------------
      !CALL GEMS_Share_MOD_log(llvl, '---   END L2 O3P Memory Allocation FUNCTION ver2 ---', '  END MSG')
      !CALL GEMS_Share_MOD_log(llvl, '.', LOGMSG)
      status = 0

      RETURN

999   CONTINUE
      llvl  =9
      write(disp_msg,'( A30, "L2 O3P Memory Allocation ver2 Failed, ierr=",i6)')ctl_msg,ierr
      !CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)
      status = ierr

      RETURN

END FUNCTION GEMS_Share_MOD_L2O3P_MemAlloc2


!-------------------------------------------------------------------------------
!+Description: 
!      FUNCTION ver2 to deallocate L2 O3P Memory.
!
!
! Method:
!
! input  :
!
! input/output :
!     local_L2O3P : L2 O3P Data
!
! output :
!       status    : the result code of L2O3P Memory deallocation
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2016.06.13  전역변수를 매개변수로 사용하기 위해 생성(YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------
! Declarations:
FUNCTION GEMS_Share_MOD_L2O3P_MemDeAlloc2(local_L2O3P) RESULT (status)

! input/output parameters:
      TYPE(L2_O3P)                  :: local_L2O3P

! Local parameters For L2 Reading:
      INTEGER (KIND = 4)            :: status

! Local parameters For Logging:


!--------------------------------
!---- FUNCTION Start -------------
!--------------------------------
      llvl=9
      LOGMSG=" PROC MSG"
      !CALL GEMS_Share_MOD_log(llvl, '--- Start L2 O3P Memory DeAllocation FUNCTION ver2 ---', LOGMSG )

      !--------------------------------------
      ! <O3Profile: Data Fields>
      !--------------------------------------

      !-----------------  AveragingKernel  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_dfld%AvgK)) DEALLOCATE(local_L2O3P%o3p_dfld%AvgK)

      !-----------------  O3  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_dfld%OzRet)) DEALLOCATE(local_L2O3P%o3p_dfld%OzRet)

      !-----------------  O3Apriori  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_dfld%OzAp)) DEALLOCATE(local_L2O3P%o3p_dfld%OzAp)

      !-----------------  O3AprioriError  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_dfld%OzApErr)) DEALLOCATE(local_L2O3P%o3p_dfld%OzApErr)

      !-----------------  O3 Random Noise Error  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_dfld%OzNErr)) DEALLOCATE(local_L2O3P%o3p_dfld%OzNErr)

      !-----------------  O3 Solution Error  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_dfld%OzSMErr)) DEALLOCATE(local_L2O3P%o3p_dfld%OzSMErr)

      !-----------------  ColumnAmountO3  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_dfld%ColAmtOz)) DEALLOCATE(local_L2O3P%o3p_dfld%ColAmtOz)

      !-----------------  DegreesOfFreedomForSignal  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_dfld%DFS)) DEALLOCATE(local_L2O3P%o3p_dfld%DFS)

      !-----------------  CloudPressure  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_dfld%CldP)) DEALLOCATE(local_L2O3P%o3p_dfld%CldP)

      !-----------------  EffectiveCloudFractionUV  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_dfld%EffCldFrac)) DEALLOCATE(local_L2O3P%o3p_dfld%EffCldFrac)

      !-----------------  TerrainReflectivityUV  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_dfld%SfcAlb)) DEALLOCATE(local_L2O3P%o3p_dfld%SfcAlb)

      !-----------------  ProcessingQualityFlags  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_dfld%PrsQFlag)) DEALLOCATE(local_L2O3P%o3p_dfld%PrsQFlag)

      !-----------------  NumberOfIterations  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_dfld%NumIter)) DEALLOCATE(local_L2O3P%o3p_dfld%NumIter)

      !-----------------  ResidualsOfFit  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_dfld%AvgRes)) DEALLOCATE(local_L2O3P%o3p_dfld%AvgRes)

      !-----------------  RootMeanSquareErrorOfFit  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_dfld%RMS)) DEALLOCATE(local_L2O3P%o3p_dfld%RMS)

      !--------------------------------------
      ! <O3Profile: Geolocation Fields>
      !--------------------------------------

      !-----------------  TropopausePressure  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_gfld%TrppsP)) DEALLOCATE(local_L2O3P%o3p_gfld%TrppsP)

      !-----------------  Pressure  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_gfld%P)) DEALLOCATE(local_L2O3P%o3p_gfld%P)

      !-----------------  Altitude  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_gfld%Alt)) DEALLOCATE(local_L2O3P%o3p_gfld%Alt)

      !-----------------  Temperature  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_gfld%TEMP)) DEALLOCATE(local_L2O3P%o3p_gfld%TEMP)

      !-----------------  TerrainPressure  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_gfld%TerrP)) DEALLOCATE(local_L2O3P%o3p_gfld%TerrP)

      !-----------------  Latitude  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_gfld%LAT)) DEALLOCATE(local_L2O3P%o3p_gfld%LAT)

      !-----------------  Longitude  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_gfld%LON)) DEALLOCATE(local_L2O3P%o3p_gfld%LON)

      !-----------------  Line       ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_gfld%Line_tmp)) DEALLOCATE(local_L2O3P%o3p_gfld%Line_tmp)

      !-----------------  Pix        ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_gfld%Pix_tmp)) DEALLOCATE(local_L2O3P%o3p_gfld%Pix_tmp)

      !-----------------  SolarZenithAngle  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_gfld%SolZenAng)) DEALLOCATE(local_L2O3P%o3p_gfld%SolZenAng)

      !-----------------  ViewingZenithAngle  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_gfld%ViewZenAng)) DEALLOCATE(local_L2O3P%o3p_gfld%ViewZenAng)

      !-----------------  RelativeAzimuthAngle  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_gfld%RelAziAng)) DEALLOCATE(local_L2O3P%o3p_gfld%RelAziAng)

      !-----------------  Time  ---------------
      IF (ASSOCIATED(local_L2O3P%o3p_gfld%Time)) DEALLOCATE(local_L2O3P%o3p_gfld%Time)

      local_L2O3P%status = 0
      write(disp_msg,'(" L2 O3P Memory DeAllocation ver2 OK !")')
      !CALL GEMS_Share_MOD_log(0, disp_msg, LOGMSG)

!--------------------------------
!---- END FUNCTION -------------
!--------------------------------
      !CALL GEMS_Share_MOD_log(llvl, '---   END L2 O3P Memory DeAllocation FUNCTION ver2 ---', LOGMSG )
      !CALL GEMS_Share_MOD_log(llvl, '.', LOGMSG)
      status = 0

      RETURN

END FUNCTION GEMS_Share_MOD_L2O3P_MemDeAlloc2


!-------------------------------------------------------------------------------
!+Description: 
!      FUNCTION ver2 to initialize L2 O3P Memory.
!
!
! Method:
!
! input  :
!
! input/output :
!     local_L2O3P : L2 O3P Data
!
! output :
!       status    : the result code of L2O3P Memory Initialization
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2016.06.13  전역변수를 매개변수로 사용하기 위해 생성(YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------
FUNCTION GEMS_Share_MOD_L2O3P_MemInit2(local_L2O3P) RESULT (status)

! input/output parameters:
      TYPE(L2_O3P)                  :: local_L2O3P

! Local parameters For L2 Reading:
      INTEGER (KIND = 4)            :: status

! Local parameters For Logging:

!--------------------------------
!---- FUNCTION Start -------------
!--------------------------------
      llvl  =9
      LOGMSG=" PROC MSG"
      !CALL GEMS_Share_MOD_log(llvl, '--- Start L2 O3P Memory Initialization FUNCTION ver2 ---', LOGMSG)
      
      IF ( local_L2O3P%status == 0 ) THEN
          llvl    =0
          disp_msg="did not allocate L2 O3P Memory"
          CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)
          status = -7001  ! MEMOEY ALLOCATION ERROR CODE = -7001
          RETURN
      END IF

      !--------------------------------------
      ! <O3Profile: Data Fields>
      !--------------------------------------
      !-----------------  AveragingKernel  ---------------
      local_L2O3P%o3p_dfld%AvgK = DEF_fVAL

      !-----------------  O3  ---------------
      local_L2O3P%o3p_dfld%OzRet = DEF_fVAL

      !-----------------  O3Apriori  ---------------
      local_L2O3P%o3p_dfld%OzAp = DEF_fVAL

      !-----------------  O3AprioriError  ---------------
      local_L2O3P%o3p_dfld%OzApErr = DEF_fVAL

      !-----------------  O3 Random Noise Error  ---------------
      local_L2O3P%o3p_dfld%OzNErr = DEF_fVAL

      !-----------------  O3 Solution Error  ---------------
      local_L2O3P%o3p_dfld%OzSMErr = DEF_fVAL

      !-----------------  ColumnAmountO3  ---------------
      local_L2O3P%o3p_dfld%ColAmtOz = DEF_fVAL

      !-----------------  DegreesOfFreedomForSignal  ---------------
      local_L2O3P%o3p_dfld%DFS = DEF_fVAL

      !-----------------  CloudPressure  ---------------
      local_L2O3P%o3p_dfld%CldP = DEF_fVAL

      !-----------------  EffectiveCloudFractionUV  ---------------
      local_L2O3P%o3p_dfld%EffCldFrac = DEF_fVAL

      !-----------------  TerrainReflectivityUV  ---------------
      local_L2O3P%o3p_dfld%SfcAlb = DEF_sVAL

      !-----------------  ProcessingQualityFlags  ---------------
      local_L2O3P%o3p_dfld%PrsQFlag = DEF_iVAL

      !-----------------  NumberOfIterations  ---------------
      local_L2O3P%o3p_dfld%NumIter = DEF_sVAL

      !-----------------  ResidualsOfFit  ---------------
      local_L2O3P%o3p_dfld%AvgRes = DEF_fVAL

      !-----------------  RootMeanSquareErrorOfFit  ---------------
      local_L2O3P%o3p_dfld%RMS = DEF_fVAL


      !--------------------------------------
      ! <O3Profile: Geolocation Fields>
      !--------------------------------------

      !-----------------  TropopausePressure  ---------------
      local_L2O3P%o3p_gfld%TrppsP = DEF_fVAL

      !-----------------  Pressure  ---------------
      local_L2O3P%o3p_gfld%P = DEF_fVAL

      !-----------------  Altitude  ---------------
      local_L2O3P%o3p_gfld%Alt = DEF_fVAL

      !-----------------  Temperature  ---------------
      local_L2O3P%o3p_gfld%TEMP = DEF_fVAL

      !-----------------  TerrainPressure -------------
      local_L2O3P%o3p_gfld%TerrP = DEF_fVAL

      !-----------------  Latitude  ---------------
      local_L2O3P%o3p_gfld%LAT = DEF_fVAL

      !-----------------  Longitude  ---------------
      local_L2O3P%o3p_gfld%LON = DEF_fVAL

      !-----------------  Line  ---------------
      local_L2O3P%o3p_gfld%Line_tmp = DEF_iVAL

      !-----------------  Pix  ---------------
      local_L2O3P%o3p_gfld%Pix_tmp = DEF_iVAL

      !-----------------  SolarZenithAngle  ---------------
      local_L2O3P%o3p_gfld%SolZenAng = DEF_fVAL

      !-----------------  ViewingZenithAngle  ---------------
      local_L2O3P%o3p_gfld%ViewZenAng = DEF_fVAL

      !-----------------  RelativeAzimuthAngle  ---------------
      local_L2O3P%o3p_gfld%RelAziAng = DEF_fVAL

      !-----------------  Time  ---------------
      local_L2O3P%o3p_gfld%Time = DEF_dVAL

      write(disp_msg,'("L2 O3P Memory Initialization ver2 OK !")')
      !CALL GEMS_Share_MOD_log(0, disp_msg, LOGMSG)

      status = 0

      RETURN

END FUNCTION GEMS_Share_MOD_L2O3P_MemInit2


END MODULE Share_l2_o3p_mod_write_read

!-------------------------------------------------------------------------------
