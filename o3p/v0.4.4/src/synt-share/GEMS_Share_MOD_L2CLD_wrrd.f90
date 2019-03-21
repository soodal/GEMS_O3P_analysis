!-------------------------------------------------------------------------------
!+Module to write and read GEMS L2File

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE Share_l2_cld_mod_write_read

!-------------------------------------------------------------------------------
!+Description: 
!     Module to wirte and read GEMS L2CLD File
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2014.12.11 Fisrt Code (HJ Lee, Seasoft) 
! 0.1.1   2016.12.13 Add memory management function ver2(YG Ki) 
!-------------------------------------------------------------------------------
      
      USE Share_MOD_Constants_Variables, ONLY:  gci_nline_l2   ,   &
                                                gci_nxtrack_l2 ,   &
                                                gci_ncomp 
      USE Share_MOD_P_Constants_Variables, ONLY: LOGMSG, llvl, ctl_msg, &
                                                 disp_msg, disp_line
      USE Share_MOD_Log
      USE Share_MOD_Def
!**********************************************************!

      IMPLICIT NONE

      PUBLIC  :: GEMS_Share_MOD_L2CLD_Write
      PUBLIC  :: GEMS_Share_MOD_L2CLD_Read
      PUBLIC  :: GEMS_Share_MOD_L2CLD_Write2
      PUBLIC  :: read_namelist_cld

      INTEGER(KIND=4), PRIVATE      :: nline_cld
      INTEGER(KIND=4), PRIVATE      :: nxtrack_cld
!      INTEGER(KIND=4), PRIVATE      :: ncomp_cld
      INTEGER(KIND=4), PARAMETER    :: DS_NAME_SZ = 200

    !----------CLD----------!
     ! <Cloud Fraction and Pressure>

      !--- Data Fields ---
      TYPE :: cld_d_fld
       REAL   (KIND=4) , POINTER    :: CldFrac         (:,:)     ! CloudFraction                  
       REAL   (KIND=4) , POINTER    :: CldFracPrec     (:,:)     ! CloudFractionPrecision         
       REAL   (KIND=4) , POINTER    :: CldP            (:,:)     ! CloudPressure                  
       REAL   (KIND=4) , POINTER    :: CldPPrec        (:,:)     ! CloudPressurePrecision         
       REAL   (KIND=4) , POINTER    :: ContmAtRefWlen  (:,:)     ! ContinuumAtReferenceWavelength 
       INTEGER(KIND=2) , POINTER    :: PrsQFlag        (:,:)     ! ProcessingQualityFlags
       REAL   (KIND=4) , POINTER    :: SltColAmtNO2    (:,:)     ! SlantColumnAmountNO2 
       REAL   (KIND=4) , POINTER    :: SltColAmtO2O2   (:,:)     ! SlantColumnAmountO2O2
       REAL   (KIND=4) , POINTER    :: SltColAmtO3     (:,:)     ! SlantColumnAmountO3          
       INTEGER(KIND=1) , POINTER    :: SnwIceExt       (:,:)     ! SnowIceExtent
       REAL   (KIND=4) , POINTER    :: TerrP           (:,:)     ! TerrainPressure                
       REAL   (KIND=4) , POINTER    :: TerrReflt       (:,:)     ! TerrainReflectivity            
      END TYPE cld_d_fld

      TYPE :: omicld_d_fld
       REAL   (KIND=4) , POINTER    :: CldFrac         (:,:)     ! CloudFraction                  
       REAL   (KIND=4) , POINTER    :: CldFracPrec     (:,:)     ! CloudFractionPrecision         
       INTEGER(KIND=2) , POINTER    :: CldP            (:,:)     ! CloudPressure                  
       INTEGER(KIND=2) , POINTER    :: CldPPrec        (:,:)     ! CloudPressurePrecision         
       REAL   (KIND=4) , POINTER    :: ContmAtRefWlen  (:,:)     ! ContinuumAtReferenceWavelength 
       INTEGER(KIND=2) , POINTER    :: PrsQFlag        (:,:)     ! ProcessingQualityFlags
       REAL   (KIND=4) , POINTER    :: SltColAmtNO2    (:,:)     ! SlantColumnAmountNO2 
       REAL   (KIND=4) , POINTER    :: SltColAmtO2O2   (:,:)     ! SlantColumnAmountO2O2
       REAL   (KIND=4) , POINTER    :: SltColAmtO3     (:,:)     ! SlantColumnAmountO3          
       INTEGER(KIND=1) , POINTER    :: SnwIceExt       (:,:)     ! SnowIceExtent
       INTEGER(KIND=2) , POINTER    :: TerrP           (:,:)     ! TerrainPressure                
       INTEGER(KIND=1) , POINTER    :: TerrReflt       (:,:)     ! TerrainReflectivity            
      END TYPE omicld_d_fld

      TYPE :: cld_g_fld
       REAL   (KIND=4) , POINTER    :: LAT             (:,:)     ! Latitude                       
       REAL   (KIND=4) , POINTER    :: LON             (:,:)     ! Longitude                      
       REAL   (KIND=4) , POINTER    :: SolAziAng       (:,:)     ! SolarAzimuthAngle              
       REAL   (KIND=4) , POINTER    :: SolZenAng       (:,:)     ! SolarZenithAngle               
       REAL   (KIND=4) , POINTER    :: ScAlt           (:)       ! SpacecraftAltitude             
       REAL   (KIND=4) , POINTER    :: ScLat           (:)       ! SpacecraftLatitude             
       REAL   (KIND=4) , POINTER    :: ScLon           (:)       ! SpacecraftLongitude            
       REAL   (KIND=4) , POINTER    :: TerrHgt         (:,:)     ! TerrainHeight                  
       REAL   (KIND=4) , POINTER    :: ViewAziAng      (:,:)     ! ViewingAzimuthAngle            
       REAL   (KIND=4) , POINTER    :: ViewZenAng      (:,:)     ! ViewingZenithAngle             
      END TYPE cld_g_fld

      TYPE :: L2_cld
       TYPE(cld_d_fld)              :: cld_dfld
       TYPE(cld_g_fld)              :: cld_gfld
       INTEGER(KIND=1)              :: status   = 0             ! 메모리를 할당 받았는지 여부를 확인
      END TYPE L2_cld

      TYPE :: OMI_L2CLD
       TYPE(omicld_d_fld)           :: omicld_dfld
       TYPE(cld_g_fld)              :: cld_gfld
       INTEGER(KIND=1)              :: status   = 0
      END TYPE OMI_L2CLD

      TYPE(L2_cld)                  :: gds_L2CLD
      TYPE(OMI_L2CLD)               :: gds_OMI_L2CLD

     ! CLD datasets in namelist
      TYPE :: CLD_ds
       CHARACTER(LEN=DS_NAME_SZ)            :: l2_cld_file_path
       CHARACTER(LEN=DS_NAME_SZ)            :: cld_dgrp_path
       CHARACTER(LEN=DS_NAME_SZ)            :: cld_ggrp_path

       INTEGER                              :: file_path_sz
       INTEGER                              :: grp_path_sz
       INTEGER                              :: data_name_sz      

       CHARACTER(LEN=DS_NAME_SZ)            :: CldFraction_DataName
       CHARACTER(LEN=DS_NAME_SZ)            :: CldFracPrec_DataName
       CHARACTER(LEN=DS_NAME_SZ)            :: CldPressure_DataName
       CHARACTER(LEN=DS_NAME_SZ)            :: CldPPrec_DataName
       CHARACTER(LEN=DS_NAME_SZ)            :: ContmAtRefWlen_DataName
       CHARACTER(LEN=DS_NAME_SZ)            :: PrsQFlag_DataName
       CHARACTER(LEN=DS_NAME_SZ)            :: SltColAmtNO2_DataName
       CHARACTER(LEN=DS_NAME_SZ)            :: SltColAmtO2O2_DataName
       CHARACTER(LEN=DS_NAME_SZ)            :: SltColAmtO3_DataName
       CHARACTER(LEN=DS_NAME_SZ)            :: SnowIceExtent_DataName
       CHARACTER(LEN=DS_NAME_SZ)            :: TerrPressure_DataName
       CHARACTER(LEN=DS_NAME_SZ)            :: TerrReflt_DataName
       CHARACTER(LEN=DS_NAME_SZ)            :: Latitude_DataName
       CHARACTER(LEN=DS_NAME_SZ)            :: Longitude_DataName
       CHARACTER(LEN=DS_NAME_SZ)            :: SolarAziAng_DataName
       CHARACTER(LEN=DS_NAME_SZ)            :: SolarZenAng_DataName
       CHARACTER(LEN=DS_NAME_SZ)            :: ScraftAlt_DataName
       CHARACTER(LEN=DS_NAME_SZ)            :: ScraftLat_DataName
       CHARACTER(LEN=DS_NAME_SZ)            :: ScraftLon_DataName
       CHARACTER(LEN=DS_NAME_SZ)            :: TerrainHgt_DataName
       CHARACTER(LEN=DS_NAME_SZ)            :: ViewAziAng_DataName
       CHARACTER(LEN=DS_NAME_SZ)            :: ViewZenAng_DataName
      END TYPE CLD_ds

      TYPE(CLD_ds)                          :: gds_CLDds


CONTAINS

!--------------------------------

!-------------------------------------------------------------------------------
!+Description: 
!      FUNCTION GEMS_Share_MOD_L2_CLD_Write_Read .
!
!
! Method:

! write and read  L2CLD files:
!       L2_cld    : the Level 2 CLD data structure
!       ctl_fpath : the control file path
!
! output files:
!       status    : the result code of L2CLD write and read
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2014.12.23  Fisrt Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!L2 CLD write function   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Declarations:
      FUNCTION GEMS_Share_MOD_L2CLD_Write(L2CLD, CLDds, ctl_fpath) RESULT (status)

! input/output parameters:
      CHARACTER(LEN=128),   INTENT(INOUT) :: ctl_fpath
      TYPE(L2_cld),         INTENT(INOUT) :: L2CLD
      TYPE(CLD_ds),         INTENT(INOUT) :: CLDds

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
      CALL GEMS_Share_MOD_log(llvl, '--- Start L2CLD WRITE FUNCTION  ---', 'START MSG')

!
!--------------------------------
!
      pstart(1) = -1                                            ! HDF5 데이터셑의 모든 데이터값을 읽는다는 의미

      CALL read_namelist_cld(ctl_fpath, CLDds)

      !----------      L2CLD     ----------!
      !--------------------------------------
      ! <Cloud Fraction And Pressure :Data Fields>
      !--------------------------------------

      disp_msg = "L2CLD WRITE FILE PATH="//CLDds%l2_cld_file_path
      CALL GEMS_Share_MOD_log(llvl, disp_msg,  LOGMSG)

      CALL GEMS_Share_MOD_log(llvl, "Cloud Fraction And Pressure :Data Fields",LOGMSG)

      !-----------------  CloudFraction --------------- 
      data_name = CLDds%CldFraction_DataName
      grp_path  = CLDds%cld_dgrp_path

      hdferr = GEMS_Share_Hdf5WriteData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_dfld%CldFrac, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_dfld%CldFrac, nxtrack_cld,nline_cld,data_name)
      ELSE
            GOTO 2999
      ENDIF
 
    !-----------------  CloudFractionPrecision  ---------------
      data_name = CLDds%CldFracPrec_DataName
      hdferr = GEMS_Share_Hdf5WriteData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_dfld%CldFracPrec, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_dfld%CldFracPrec, nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  CloudPressure  ---------------
      data_name = CLDds%CldPressure_DataName
      hdferr = GEMS_Share_Hdf5WriteData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_dfld%CldP, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_dfld%CldP, nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  CloudPressurePrecision  ---------------
      data_name = CLDds%CldPPrec_DataName           
      hdferr = GEMS_Share_Hdf5WriteData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_dfld%CldPPrec, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_dfld%CldPPrec, nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  ContinuumAtReferenceWavelength  ---------------
      data_name = CLDds%ContmAtRefWlen_DataName           
      hdferr = GEMS_Share_Hdf5WriteData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_dfld%ContmAtRefWlen, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_dfld%ContmAtRefWlen, nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SlantColumnAmountNO2  ---------------
      data_name = CLDds%SltColAmtNO2_DataName
      hdferr = GEMS_Share_Hdf5WriteData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_dfld%SltColAmtNO2, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl,L2CLD%cld_dfld%SltColAmtNO2, nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  ProcessingQualityFlags  ---------------
      data_name = CLDds%PrsQFlag_DataName
      hdferr = GEMS_Share_Hdf5WriteData(CLDds%l2_cld_file_path,CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz,L2CLD%cld_dfld%PrsQFlag, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_dfld%PrsQFlag,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SlantColumnAmountO2O2  ---------------
      data_name = CLDds%SltColAmtO2O2_DataName           
      hdferr = GEMS_Share_Hdf5WriteData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_dfld%SltColAmtO2O2, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_dfld%SltColAmtO2O2, nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SlantColumnAmountO3  ---------------
      data_name = CLDds%SltColAmtO3_DataName
      hdferr = GEMS_Share_Hdf5WriteData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                             data_name, CLDds%data_name_sz, L2CLD%cld_dfld%SltColAmtO3, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl,L2CLD%cld_dfld%SltColAmtO3,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SnowIceExtent  ---------------
      data_name = CLDds%SnowIceExtent_DataName
      hdferr = GEMS_Share_Hdf5WriteData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_dfld%SnwIceExt, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl,L2CLD%cld_dfld%SnwIceExt,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  TerrainPressure  ---------------
      data_name = CLDds%TerrPressure_DataName           
      hdferr = GEMS_Share_Hdf5WriteData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_dfld%TerrP, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_dfld%TerrP, nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  TerrainReflectivity  ---------------
      data_name = CLDds%TerrReflt_DataName           
      hdferr = GEMS_Share_Hdf5WriteData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_dfld%TerrReflt, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_dfld%TerrReflt, nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !--------------------------------------
      ! <Cloud Fraction And Pressure: Geolocation Fields>
      !--------------------------------------
      CALL GEMS_Share_MOD_log(llvl, "Cloud Fraction And Pressure: Geolocation Fields", LOGMSG)


      !-----------------  Latitude  ---------------
      data_name = CLDds%Latitude_DataName
      grp_path  = CLDds%cld_ggrp_path
           
      hdferr = GEMS_Share_Hdf5WriteData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_gfld%LAT, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%LAT, nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  Longitude  ---------------
      data_name = CLDds%Longitude_DataName           
      hdferr = GEMS_Share_Hdf5WriteData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_gfld%LON, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%LON, nxtrack_cld,nline_cld,data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SolarAzimuthAangle  ---------------
      data_name = CLDds%SolarAziAng_DataName           
      hdferr = GEMS_Share_Hdf5WriteData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_gfld%SolAziAng, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%SolAziAng, nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SolarZenithAangle  ---------------
      data_name = CLDds%SolarZenAng_DataName           
      hdferr = GEMS_Share_Hdf5WriteData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_gfld%SolZenAng, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%SolZenAng,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SpacecraftAltitude  ---------------
      data_name = CLDds%ScraftAlt_DataName           
      hdferr = GEMS_Share_Hdf5WriteData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_gfld%ScAlt, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%ScAlt, nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SpacecraftLatitude  ---------------
      data_name = CLDds%ScraftLat_DataName           
      hdferr = GEMS_Share_Hdf5WriteData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_gfld%ScLat, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%ScLat, nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF


      !-----------------  SpacecraftLongitude  ---------------
      data_name = CLDds%ScraftLon_DataName           
      hdferr = GEMS_Share_Hdf5WriteData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_gfld%ScLon, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%ScLon, nline_cld,data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  TerrainHeight  ---------------
      data_name = CLDds%TerrainHgt_DataName           
      hdferr = GEMS_Share_Hdf5WriteData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_gfld%TerrHgt, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%TerrHgt, nxtrack_cld,nline_cld,data_name)
      ELSE
            GOTO 2999
      ENDIF


      !-----------------  ViewingAzimuthAngle  ---------------
      data_name = CLDds%ViewAziAng_DataName           
      hdferr = GEMS_Share_Hdf5WriteData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_gfld%ViewAziAng, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%ViewAziAng,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  ViewingZenithAngle  ---------------
      data_name = CLDds%ViewZenAng_DataName
      hdferr = GEMS_Share_Hdf5WriteData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_gfld%ViewZenAng, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%ViewZenAng,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      write(disp_msg,'(" GEMS HDF5 File WRITE OK, hdferr=",i6)') hdferr
      CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)

!--------------------------------
!---- END FUNCTION -------------
!--------------------------------

      CALL GEMS_Share_MOD_log(llvl, '---   END L2 WRITE FUNCTION ---', '  END MSG')
      CALL GEMS_Share_MOD_log(llvl, '.', LOGMSG)
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

      END FUNCTION GEMS_Share_MOD_L2CLD_Write



!!!!!!!!!!!!!!!!!!!!!!!!!!!
!L2 CLD read function   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Declarations:
      FUNCTION GEMS_Share_MOD_L2CLD_Read(L2CLD, CLDds, ctl_fpath) RESULT (status)

! input/output parameters:
      CHARACTER(LEN=128),   INTENT(INOUT) :: ctl_fpath
      TYPE(L2_cld),         INTENT(INOUT) :: L2CLD
      TYPE(CLD_ds),         INTENT(INOUT) :: CLDds

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
      CALL GEMS_Share_MOD_log(llvl, '--- Start L2CLD READ FUNCTION  ---', 'STARTMSG')


!--------------------------------
!
      pstart(1) = -1                                            ! HDF5 데이터셑의 모든 데이터값을 읽는다는 의미

      CALL read_namelist_cld(ctl_fpath, CLDds)

!      !--- Setting to Global Variable of CLD namelist
!      gds_CLDds = CLDds

      !----------      L2CLD     ----------!
      !--------------------------------------
      ! <Cloud Fraction And Pressure :Data Fields>
      !--------------------------------------

      disp_msg = "L2CLD READ FILE PATH="//CLDds%l2_cld_file_path
      CALL GEMS_Share_MOD_log(llvl, disp_msg,  LOGMSG)

      CALL GEMS_Share_MOD_log(llvl, "Cloud Fraction And Pressure :Data Fields",LOGMSG)


      !-----------------  CloudFraction ---------------
      IF (ASSOCIATED(L2CLD%cld_dfld%CldFrac)) DEALLOCATE(L2CLD%cld_dfld%CldFrac)
      data_name = CLDds%CldFraction_DataName
      grp_path  = CLDds%cld_dgrp_path
      ALLOCATE (L2CLD%cld_dfld%CldFrac(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%cld_dfld%CldFrac, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_dfld%CldFrac, nxtrack_cld,nline_cld,data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  CloudFractionPrecision  ---------------
      IF (ASSOCIATED(L2CLD%cld_dfld%CldFracPrec)) DEALLOCATE(L2CLD%cld_dfld%CldFracPrec)
      data_name = CLDds%CldFracPrec_DataName
      ALLOCATE (L2CLD%cld_dfld%CldFracPrec(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%cld_dfld%CldFracPrec, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_dfld%CldFracPrec, nline_cld, nxtrack_cld,data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  CloudPressure  ---------------
      IF (ASSOCIATED(L2CLD%cld_dfld%CldP)) DEALLOCATE(L2CLD%cld_dfld%CldP)
      data_name = CLDds%CldPressure_DataName
      ALLOCATE (L2CLD%cld_dfld%CldP(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%cld_dfld%CldP, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_dfld%CldP, nxtrack_cld,nline_cld,data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  CloudPressurePrecision  ---------------
      IF (ASSOCIATED(L2CLD%cld_dfld%CldPPrec)) DEALLOCATE(L2CLD%cld_dfld%CldPPrec)
      data_name = CLDds%CldPPrec_DataName
      ALLOCATE (L2CLD%cld_dfld%CldPPrec(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%cld_dfld%CldPPrec, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_dfld%CldPPrec,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  ContinuumAtReferenceWavelength  ---------------
      IF (ASSOCIATED(L2CLD%cld_dfld%ContmAtRefWlen)) DEALLOCATE(L2CLD%cld_dfld%ContmAtRefWlen)
      data_name = CLDds%ContmAtRefWlen_DataName
      ALLOCATE (L2CLD%cld_dfld%ContmAtRefWlen(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%cld_dfld%ContmAtRefWlen, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_dfld%ContmAtRefWlen,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  ProcessingQualityFlags  ---------------
      IF (ASSOCIATED(L2CLD%cld_dfld%PrsQFlag)) DEALLOCATE(L2CLD%cld_dfld%PrsQFlag)
      data_name = CLDds%PrsQFlag_DataName
      ALLOCATE (L2CLD%cld_dfld%PrsQFlag(nxtrack_cld,nline_cld), STAT =ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path,CLDds%file_path_sz, grp_path,CLDds%grp_path_sz, &
                                        data_name,CLDds%data_name_sz,L2CLD%cld_dfld%PrsQFlag, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl,L2CLD%cld_dfld%PrsQFlag,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SlantColumnAmountNO2  ---------------
      IF (ASSOCIATED(L2CLD%cld_dfld%SltColAmtNO2)) DEALLOCATE(L2CLD%cld_dfld%SltColAmtNO2)
      data_name = CLDds%SltColAmtNO2_DataName
      ALLOCATE (L2CLD%cld_dfld%SltColAmtNO2(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%cld_dfld%SltColAmtNO2, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl,L2CLD%cld_dfld%SltColAmtNO2,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SlantColumnAmountO2O2  ---------------
      IF (ASSOCIATED(L2CLD%cld_dfld%SltColAmtO2O2)) DEALLOCATE(L2CLD%cld_dfld%SltColAmtO2O2)
      data_name = CLDds%SltColAmtO2O2_DataName
      ALLOCATE (L2CLD%cld_dfld%SltColAmtO2O2( nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%cld_dfld%SltColAmtO2O2, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_dfld%SltColAmtO2O2, nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SlantColumnAmountO3  ---------------
      IF (ASSOCIATED(L2CLD%cld_dfld%SltColAmtO3)) DEALLOCATE(L2CLD%cld_dfld%SltColAmtO3)
      data_name = CLDds%SltColAmtO3_DataName
      ALLOCATE (L2CLD%cld_dfld%SltColAmtO3(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%cld_dfld%SltColAmtO3, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl,L2CLD%cld_dfld%SltColAmtO3,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SnowIceExtent  ---------------
      IF (ASSOCIATED(L2CLD%cld_dfld%SnwIceExt)) DEALLOCATE(L2CLD%cld_dfld%SnwIceExt)
      data_name = CLDds%SnowIceExtent_DataName
      ALLOCATE (L2CLD%cld_dfld%SnwIceExt(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%cld_dfld%SnwIceExt, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl,L2CLD%cld_dfld%SnwIceExt,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  TerrainPressure  ---------------
      IF (ASSOCIATED(L2CLD%cld_dfld%TerrP)) DEALLOCATE(L2CLD%cld_dfld%TerrP)
      data_name = CLDds%TerrPressure_DataName
      ALLOCATE (L2CLD%cld_dfld%TerrP(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%cld_dfld%TerrP, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_dfld%TerrP, nxtrack_cld,nline_cld,data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  TerrainReflectivity  ---------------
      IF (ASSOCIATED(L2CLD%cld_dfld%TerrReflt)) DEALLOCATE(L2CLD%cld_dfld%TerrReflt)
      data_name = CLDds%TerrReflt_DataName
      ALLOCATE (L2CLD%cld_dfld%TerrReflt(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%cld_dfld%TerrReflt, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_dfld%TerrReflt,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !--------------------------------------
      ! <Cloud Fraction And Pressure: Geolocation Fields>
      !--------------------------------------
      CALL GEMS_Share_MOD_log(llvl, "Cloud Fraction And Pressure: Geolocation Fields", LOGMSG)

      !-----------------  Latitude  ---------------
      IF (ASSOCIATED(L2CLD%cld_gfld%LAT)) DEALLOCATE(L2CLD%cld_gfld%LAT)
      data_name = CLDds%Latitude_DataName
      grp_path  = CLDds%cld_ggrp_path
      ALLOCATE (L2CLD%cld_gfld%LAT(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%cld_gfld%LAT, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%LAT, nxtrack_cld,nline_cld,data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  Longitude  ---------------
      IF (ASSOCIATED(L2CLD%cld_gfld%LON)) DEALLOCATE(L2CLD%cld_gfld%LON)
      data_name = CLDds%Longitude_DataName
      ALLOCATE (L2CLD%cld_gfld%LON(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%cld_gfld%LON, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%LON, nxtrack_cld,nline_cld,data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SolarAzimuthAangle  ---------------
      IF (ASSOCIATED(L2CLD%cld_gfld%SolAziAng)) DEALLOCATE(L2CLD%cld_gfld%SolAziAng)
      data_name = CLDds%SolarAziAng_DataName
      ALLOCATE (L2CLD%cld_gfld%SolAziAng(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%cld_gfld%SolAziAng, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%SolAziAng,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SolarZenithAangle  ---------------
      IF (ASSOCIATED(L2CLD%cld_gfld%SolZenAng)) DEALLOCATE(L2CLD%cld_gfld%SolZenAng)
      data_name = CLDds%SolarZenAng_DataName
      ALLOCATE (L2CLD%cld_gfld%SolZenAng(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%cld_gfld%SolZenAng, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%SolZenAng,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SpacecraftAltitude  ---------------
      IF (ASSOCIATED(L2CLD%cld_gfld%ScAlt)) DEALLOCATE(L2CLD%cld_gfld%ScAlt)
      data_name = CLDds%ScraftAlt_DataName
      ALLOCATE (L2CLD%cld_gfld%ScAlt(nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%cld_gfld%ScAlt, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%ScAlt, nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SpacecraftLatitude  ---------------
      IF (ASSOCIATED(L2CLD%cld_gfld%ScLat)) DEALLOCATE(L2CLD%cld_gfld%ScLat)
      data_name = CLDds%ScraftLat_DataName
      ALLOCATE (L2CLD%cld_gfld%ScLat(nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%cld_gfld%ScLat, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%ScLat, nline_cld,data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SpacecraftLongitude  ---------------
      IF (ASSOCIATED(L2CLD%cld_gfld%ScLon)) DEALLOCATE(L2CLD%cld_gfld%ScLon)
      data_name = CLDds%ScraftLon_DataName
      ALLOCATE (L2CLD%cld_gfld%ScLon(nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%cld_gfld%ScLon, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%ScLon, nline_cld,data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  TerrainHeight  ---------------
      IF (ASSOCIATED(L2CLD%cld_gfld%TerrHgt)) DEALLOCATE(L2CLD%cld_gfld%TerrHgt)
      data_name = CLDds%TerrainHgt_DataName
      ALLOCATE (L2CLD%cld_gfld%TerrHgt(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%cld_gfld%TerrHgt, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%TerrHgt, nxtrack_cld,nline_cld,data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  ViewingAzimuthAngle  ---------------
      IF (ASSOCIATED(L2CLD%cld_gfld%ViewAziAng)) DEALLOCATE(L2CLD%cld_gfld%ViewAziAng)
      data_name = CLDds%ViewAziAng_DataName
      ALLOCATE (L2CLD%cld_gfld%ViewAziAng(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%cld_gfld%ViewAziAng, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%ViewAziAng,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  ViewingZenithAngle  ---------------
      IF (ASSOCIATED(L2CLD%cld_gfld%ViewZenAng)) DEALLOCATE(L2CLD%cld_gfld%ViewZenAng)
      data_name = CLDds%ViewZenAng_DataName
      ALLOCATE (L2CLD%cld_gfld%ViewZenAng(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%cld_gfld%ViewZenAng, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%ViewZenAng,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      write(disp_msg,'(" GEMS HDF5 File READ OK, hdferr=",i6)') hdferr
      CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)

!--------------------------------
!---- END FUNCTION -------------
!--------------------------------
      CALL GEMS_Share_MOD_log(llvl, '---   END L2 READ FUNCTION ---', '  END MSG')
      CALL GEMS_Share_MOD_log(llvl, '.', LOGMSG)
      status = 0

      RETURN

999   CONTINUE
      write(disp_msg,'( A20, "Memory Allocation Failed, ierr=",i6)') ctl_msg,ierr
      CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)
      status = ierr

      RETURN

2999  CONTINUE
      write(disp_msg,'(" !!! GEMS HDF5 File READ ERROR !!! ", "[", A30, "],hdferr=",i6)') data_name, hdferr
      CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)
      status = hdferr

      RETURN

      END FUNCTION GEMS_Share_MOD_L2CLD_Read


!!!!!!!!!!!!!!!!!!!!!!!!!!!
!OMI L2 CLD read function   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Declarations:
      FUNCTION GEMS_Share_MOD_OMI_L2CLD_Read(L2CLD, CLDds, ctl_fpath) RESULT (status)

! input/output parameters:
      CHARACTER(LEN=128),   INTENT(INOUT) :: ctl_fpath
      TYPE(OMI_L2CLD),      INTENT(INOUT) :: L2CLD
      TYPE(CLD_ds),         INTENT(INOUT) :: CLDds

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
      CALL GEMS_Share_MOD_log(llvl, '--- Start L2CLD READ FUNCTION  ---', 'STARTMSG')


!--------------------------------
!
      pstart(1) = -1                                            ! HDF5 데이터셑의 모든 데이터값을 읽는다는 의미

      CALL read_namelist_cld(ctl_fpath, CLDds)

!      !--- Setting to Global Variable of CLD namelist
!      gds_CLDds = CLDds

      !----------      L2CLD     ----------!
      !--------------------------------------
      ! <Cloud Fraction And Pressure :Data Fields>
      !--------------------------------------

      disp_msg = "L2CLD READ FILE PATH="//CLDds%l2_cld_file_path
      CALL GEMS_Share_MOD_log(llvl, disp_msg,  LOGMSG)

      CALL GEMS_Share_MOD_log(llvl, "Cloud Fraction And Pressure :Data Fields",LOGMSG)


      !-----------------  CloudFraction ---------------
      IF (ASSOCIATED(L2CLD%omicld_dfld%CldFrac)) DEALLOCATE(L2CLD%omicld_dfld%CldFrac)
      data_name = CLDds%CldFraction_DataName
      grp_path  = CLDds%cld_dgrp_path
      ALLOCATE (L2CLD%omicld_dfld%CldFrac(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%omicld_dfld%CldFrac, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%omicld_dfld%CldFrac, nxtrack_cld,nline_cld,data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  CloudFractionPrecision  ---------------
      IF (ASSOCIATED(L2CLD%omicld_dfld%CldFracPrec)) DEALLOCATE(L2CLD%omicld_dfld%CldFracPrec)
      data_name = CLDds%CldFracPrec_DataName
      ALLOCATE (L2CLD%omicld_dfld%CldFracPrec(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%omicld_dfld%CldFracPrec, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%omicld_dfld%CldFracPrec, nline_cld, nxtrack_cld,data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  CloudPressure  ---------------
      IF (ASSOCIATED(L2CLD%omicld_dfld%CldP)) DEALLOCATE(L2CLD%omicld_dfld%CldP)
      data_name = CLDds%CldPressure_DataName
      ALLOCATE (L2CLD%omicld_dfld%CldP(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%omicld_dfld%CldP, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%omicld_dfld%CldP, nxtrack_cld,nline_cld,data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  CloudPressurePrecision  ---------------
      IF (ASSOCIATED(L2CLD%omicld_dfld%CldPPrec)) DEALLOCATE(L2CLD%omicld_dfld%CldPPrec)
      data_name = CLDds%CldPPrec_DataName
      ALLOCATE (L2CLD%omicld_dfld%CldPPrec(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%omicld_dfld%CldPPrec, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%omicld_dfld%CldPPrec,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  ContinuumAtReferenceWavelength  ---------------
      IF (ASSOCIATED(L2CLD%omicld_dfld%ContmAtRefWlen)) DEALLOCATE(L2CLD%omicld_dfld%ContmAtRefWlen)
      data_name = CLDds%ContmAtRefWlen_DataName
      ALLOCATE (L2CLD%omicld_dfld%ContmAtRefWlen(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%omicld_dfld%ContmAtRefWlen, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%omicld_dfld%ContmAtRefWlen,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  ProcessingQualityFlags  ---------------
      IF (ASSOCIATED(L2CLD%omicld_dfld%PrsQFlag)) DEALLOCATE(L2CLD%omicld_dfld%PrsQFlag)
      data_name = CLDds%PrsQFlag_DataName
      ALLOCATE (L2CLD%omicld_dfld%PrsQFlag(nxtrack_cld,nline_cld), STAT =ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path,CLDds%file_path_sz, grp_path,CLDds%grp_path_sz, &
                                        data_name,CLDds%data_name_sz,L2CLD%omicld_dfld%PrsQFlag, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl,L2CLD%omicld_dfld%PrsQFlag,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SlantColumnAmountNO2  ---------------
      IF (ASSOCIATED(L2CLD%omicld_dfld%SltColAmtNO2)) DEALLOCATE(L2CLD%omicld_dfld%SltColAmtNO2)
      data_name = CLDds%SltColAmtNO2_DataName
      ALLOCATE (L2CLD%omicld_dfld%SltColAmtNO2(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%omicld_dfld%SltColAmtNO2, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl,L2CLD%omicld_dfld%SltColAmtNO2,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SlantColumnAmountO2O2  ---------------
      IF (ASSOCIATED(L2CLD%omicld_dfld%SltColAmtO2O2)) DEALLOCATE(L2CLD%omicld_dfld%SltColAmtO2O2)
      data_name = CLDds%SltColAmtO2O2_DataName
      ALLOCATE (L2CLD%omicld_dfld%SltColAmtO2O2(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%omicld_dfld%SltColAmtO2O2, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%omicld_dfld%SltColAmtO2O2, nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SlantColumnAmountO3  ---------------
      IF (ASSOCIATED(L2CLD%omicld_dfld%SltColAmtO3)) DEALLOCATE(L2CLD%omicld_dfld%SltColAmtO3)
      data_name = CLDds%SltColAmtO3_DataName
      ALLOCATE (L2CLD%omicld_dfld%SltColAmtO3(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%omicld_dfld%SltColAmtO3, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl,L2CLD%omicld_dfld%SltColAmtO3,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SnowIceExtent  ---------------
      IF (ASSOCIATED(L2CLD%omicld_dfld%SnwIceExt)) DEALLOCATE(L2CLD%omicld_dfld%SnwIceExt)
      data_name = CLDds%SnowIceExtent_DataName
      ALLOCATE (L2CLD%omicld_dfld%SnwIceExt(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%omicld_dfld%SnwIceExt, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl,L2CLD%omicld_dfld%SnwIceExt,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  TerrainPressure  ---------------
      IF (ASSOCIATED(L2CLD%omicld_dfld%TerrP)) DEALLOCATE(L2CLD%omicld_dfld%TerrP)
      data_name = CLDds%TerrPressure_DataName
      ALLOCATE (L2CLD%omicld_dfld%TerrP(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%omicld_dfld%TerrP, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%omicld_dfld%TerrP, nxtrack_cld,nline_cld,data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  TerrainReflectivity  ---------------
      IF (ASSOCIATED(L2CLD%omicld_dfld%TerrReflt)) DEALLOCATE(L2CLD%omicld_dfld%TerrReflt)
      data_name = CLDds%TerrReflt_DataName
      ALLOCATE (L2CLD%omicld_dfld%TerrReflt(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%omicld_dfld%TerrReflt, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%omicld_dfld%TerrReflt,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !--------------------------------------
      ! <Cloud Fraction And Pressure: Geolocation Fields>
      !--------------------------------------
      CALL GEMS_Share_MOD_log(llvl, "Cloud Fraction And Pressure: Geolocation Fields", LOGMSG)

      !-----------------  Latitude  ---------------
      IF (ASSOCIATED(L2CLD%cld_gfld%LAT)) DEALLOCATE(L2CLD%cld_gfld%LAT)
      data_name = CLDds%Latitude_DataName
      grp_path  = CLDds%cld_ggrp_path
      ALLOCATE (L2CLD%cld_gfld%LAT(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%cld_gfld%LAT, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%LAT, nxtrack_cld,nline_cld,data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  Longitude  ---------------
      IF (ASSOCIATED(L2CLD%cld_gfld%LON)) DEALLOCATE(L2CLD%cld_gfld%LON)
      data_name = CLDds%Longitude_DataName
      ALLOCATE (L2CLD%cld_gfld%LON(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%cld_gfld%LON, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%LON, nxtrack_cld,nline_cld,data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SolarAzimuthAangle  ---------------
      IF (ASSOCIATED(L2CLD%cld_gfld%SolAziAng)) DEALLOCATE(L2CLD%cld_gfld%SolAziAng)
      data_name = CLDds%SolarAziAng_DataName
      ALLOCATE (L2CLD%cld_gfld%SolAziAng(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%cld_gfld%SolAziAng, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%SolAziAng,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SolarZenithAangle  ---------------
      IF (ASSOCIATED(L2CLD%cld_gfld%SolZenAng)) DEALLOCATE(L2CLD%cld_gfld%SolZenAng)
      data_name = CLDds%SolarZenAng_DataName
      ALLOCATE (L2CLD%cld_gfld%SolZenAng(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%cld_gfld%SolZenAng, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%SolZenAng,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SpacecraftAltitude  ---------------
      IF (ASSOCIATED(L2CLD%cld_gfld%ScAlt)) DEALLOCATE(L2CLD%cld_gfld%ScAlt)
      data_name = CLDds%ScraftAlt_DataName
      ALLOCATE (L2CLD%cld_gfld%ScAlt(nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%cld_gfld%ScAlt, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%ScAlt, nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SpacecraftLatitude  ---------------
      IF (ASSOCIATED(L2CLD%cld_gfld%ScLat)) DEALLOCATE(L2CLD%cld_gfld%ScLat)
      data_name = CLDds%ScraftLat_DataName
      ALLOCATE (L2CLD%cld_gfld%ScLat(nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%cld_gfld%ScLat, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%ScLat, nline_cld,data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SpacecraftLongitude  ---------------
      IF (ASSOCIATED(L2CLD%cld_gfld%ScLon)) DEALLOCATE(L2CLD%cld_gfld%ScLon)
      data_name = CLDds%ScraftLon_DataName
      ALLOCATE (L2CLD%cld_gfld%ScLon(nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%cld_gfld%ScLon, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%ScLon, nline_cld,data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  TerrainHeight  ---------------
      IF (ASSOCIATED(L2CLD%cld_gfld%TerrHgt)) DEALLOCATE(L2CLD%cld_gfld%TerrHgt)
      data_name = CLDds%TerrainHgt_DataName
      ALLOCATE (L2CLD%cld_gfld%TerrHgt(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%cld_gfld%TerrHgt, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%TerrHgt, nxtrack_cld,nline_cld,data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  ViewingAzimuthAngle  ---------------
      IF (ASSOCIATED(L2CLD%cld_gfld%ViewAziAng)) DEALLOCATE(L2CLD%cld_gfld%ViewAziAng)
      data_name = CLDds%ViewAziAng_DataName
      ALLOCATE (L2CLD%cld_gfld%ViewAziAng(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%cld_gfld%ViewAziAng, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%ViewAziAng,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  ViewingZenithAngle  ---------------
      IF (ASSOCIATED(L2CLD%cld_gfld%ViewZenAng)) DEALLOCATE(L2CLD%cld_gfld%ViewZenAng)
      data_name = CLDds%ViewZenAng_DataName
      ALLOCATE (L2CLD%cld_gfld%ViewZenAng(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(CLDds%l2_cld_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                       data_name, CLDds%data_name_sz, L2CLD%cld_gfld%ViewZenAng, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%ViewZenAng,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      write(disp_msg,'(" GEMS HDF5 File READ OK, hdferr=",i6)') hdferr
      CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)

!--------------------------------
!---- END FUNCTION -------------
!--------------------------------
      CALL GEMS_Share_MOD_log(llvl, '---   END L2 READ FUNCTION ---', '  END MSG')
      CALL GEMS_Share_MOD_log(llvl, '.', LOGMSG)
      status = 0

      RETURN

999   CONTINUE
      write(disp_msg,'( A20, "Memory Allocation Failed, ierr=",i6)') ctl_msg,ierr
      CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)
      status = ierr

      RETURN

2999  CONTINUE
      write(disp_msg,'(" !!! GEMS HDF5 File READ ERROR !!! ", "[", A30, "],hdferr=",i6)') data_name, hdferr
      CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)
      status = hdferr

      RETURN

      END FUNCTION GEMS_Share_MOD_OMI_L2CLD_Read

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Read NameList function   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE read_namelist_cld(ctl_fpath, CLDds)

! input/output parameters:
      CHARACTER(LEN=128),   INTENT(INOUT) :: ctl_fpath
      TYPE(CLD_ds),         INTENT(INOUT) :: CLDds

! Local parameters For L2 Reading:
      CHARACTER*200                       :: l2_cld_file_path
      CHARACTER*200                       :: cld_dgrp_path
      CHARACTER*200                       :: cld_ggrp_path

      CHARACTER*200                       :: CldFraction_DataName
      CHARACTER*200                       :: CldFracPrec_DataName
      CHARACTER*200                       :: CldPressure_DataName
      CHARACTER*200                       :: CldPPrec_DataName
      CHARACTER*200                       :: ContmAtRefWlen_DataName
      CHARACTER*200                       :: PrsQFlag_DataName
      CHARACTER*200                       :: SltColAmtNO2_DataName
      CHARACTER*200                       :: SltColAmtO2O2_DataName
      CHARACTER*200                       :: SltColAmtO3_DataName
      CHARACTER*200                       :: SnowIceExtent_DataName
      CHARACTER*200                       :: TerrPressure_DataName
      CHARACTER*200                       :: TerrReflt_DataName
      CHARACTER*200                       :: Latitude_DataName
      CHARACTER*200                       :: Longitude_DataName
      CHARACTER*200                       :: SolarAziAng_DataName
      CHARACTER*200                       :: SolarZenAng_DataName
      CHARACTER*200                       :: ScraftAlt_DataName
      CHARACTER*200                       :: ScraftLat_DataName
      CHARACTER*200                       :: ScraftLon_DataName
      CHARACTER*200                       :: TerrainHgt_DataName
      CHARACTER*200                       :: ViewAziAng_DataName
      CHARACTER*200                       :: ViewZenAng_DataName

      INTEGER                             :: file_path_sz
      INTEGER                             :: grp_path_sz
      INTEGER                             :: data_name_sz


      Namelist /L2_CLD_File_List/l2_cld_file_path,  &
                                 cld_dgrp_path,     &
                                 cld_ggrp_path
      Namelist /L2_File_List_Value_Size/file_path_sz, &
                                        grp_path_sz,  &
                                        data_name_sz

      Namelist /L2_CLD_DATA_List/CldFraction_DataName,   &
                                 CldFracPrec_DataName,   &
                                 CldPressure_DataName,   &
                                 CldPPrec_DataName,      &
                                 ContmAtRefWlen_DataName,&
                                 PrsQFlag_DataName,      &
                                 SltColAmtNO2_DataName,  &
                                 SltColAmtO2O2_DataName, &
                                 SltColAmtO3_DataName,   &
                                 SnowIceExtent_DataName, &
                                 TerrPressure_DataName,  &
                                 TerrReflt_DataName,     &
                                 Latitude_DataName,      &
                                 Longitude_DataName,     &
                                 SolarAziAng_DataName,   &
                                 SolarZenAng_DataName,   &
                                 ScraftAlt_DataName,     &
                                 ScraftLat_DataName,     &
                                 ScraftLon_DataName,     &
                                 TerrainHgt_DataName,    &
                                 ViewAziAng_DataName,    &
                                 ViewZenAng_DataName 


!--------------------------------
!---     Read Control file    ---
!--------------------------------
      OPEN(10, file=trim(ctl_fpath), status='old')
      READ(10, L2_CLD_File_List)
      READ(10, L2_CLD_DATA_List)
      READ(10, L2_File_List_Value_Size)
      CLOSE(10)

      CLDds%CldFraction_DataName    = CldFraction_DataName
      CLDds%CldFracPrec_DataName    = CldFracPrec_DataName
      CLDds%CldPressure_DataName    = CldPressure_DataName
      CLDds%CldPPrec_DataName       = CldPPrec_DataName
      CLDds%ContmAtRefWlen_DataName = ContmAtRefWlen_DataName
      CLDds%PrsQFlag_DataName       = PrsQFlag_DataName
      CLDds%SltColAmtNO2_DataName   = SltColAmtNO2_DataName
      CLDds%SltColAmtO2O2_DataName  = SltColAmtO2O2_DataName
      CLDds%SltColAmtO3_DataName    = SltColAmtO3_DataName
      CLDds%SnowIceExtent_DataName  = SnowIceExtent_DataName
      CLDds%TerrPressure_DataName   = TerrPressure_DataName
      CLDds%TerrReflt_DataName      = TerrReflt_DataName
      CLDds%Latitude_DataName       = Latitude_DataName
      CLDds%Longitude_DataName      = Longitude_DataName
      CLDds%SolarAziAng_DataName    = SolarAziAng_DataName
      CLDds%SolarZenAng_DataName    = SolarZenAng_DataName
      CLDds%ScraftAlt_DataName      = ScraftAlt_DataName
      CLDds%ScraftLat_DataName      = ScraftLat_DataName
      CLDds%ScraftLon_DataName      = ScraftLon_DataName
      CLDds%TerrainHgt_DataName     = TerrainHgt_DataName
      CLDds%ViewAziAng_DataName     = ViewAziAng_DataName
      CLDds%ViewZenAng_DataName     = ViewZenAng_DataName

      CLDds%l2_cld_file_path      = l2_cld_file_path
      CLDds%cld_dgrp_path         = cld_dgrp_path
      CLDds%cld_ggrp_path         = cld_ggrp_path
      CLDds%file_path_sz          = file_path_sz
      CLDds%grp_path_sz           = grp_path_sz
      CLDds%data_name_sz          = data_name_sz

      END SUBROUTINE read_namelist_cld


      SUBROUTINE GEMS_Share_Init_L2CLD_GlobalConstants
      ! Size for L2CLD dataset
      nline_cld       =    gci_nline_l2
      nxtrack_cld     =    gci_nxtrack_l2
!      ncomp_cld       =    gci_ncomp
      END SUBROUTINE GEMS_Share_Init_L2CLD_GlobalConstants

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!L2 CLD Memory Allocation function   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Declarations:
      FUNCTION GEMS_Share_MOD_L2CLD_MemAlloc() RESULT (status)

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
      CALL GEMS_Share_MOD_log(llvl, '--- Start L2 CLD Memory Allocation FUNCTION  ---', LOGMSG )

      !-----------------  CloudFraction ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_dfld%CldFrac)) DEALLOCATE(gds_L2CLD%cld_dfld%CldFrac)
      ALLOCATE (gds_L2CLD%cld_dfld%CldFrac(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="CloudFraction"
              GOTO 2999
           ENDIF

      !-----------------  CloudFractionPrecision  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_dfld%CldFracPrec)) DEALLOCATE(gds_L2CLD%cld_dfld%CldFracPrec)
      ALLOCATE (gds_L2CLD%cld_dfld%CldFracPrec(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="CloudFractionPrecision"
              GOTO 2999
           ENDIF

      !-----------------  CloudPressure  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_dfld%CldP)) DEALLOCATE(gds_L2CLD%cld_dfld%CldP)
      ALLOCATE (gds_L2CLD%cld_dfld%CldP(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="CloudPressure"
              GOTO 2999
           ENDIF

      !-----------------  CloudPressurePrecision  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_dfld%CldPPrec)) DEALLOCATE(gds_L2CLD%cld_dfld%CldPPrec)
      ALLOCATE (gds_L2CLD%cld_dfld%CldPPrec(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="CloudPressurePrecision"
              GOTO 2999
           ENDIF

      !-----------------  ContinuumAtReferenceWavelength  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_dfld%ContmAtRefWlen)) DEALLOCATE(gds_L2CLD%cld_dfld%ContmAtRefWlen)
      ALLOCATE (gds_L2CLD%cld_dfld%ContmAtRefWlen(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="ContinuumAtReferenceWavelength"
              GOTO 2999
           ENDIF

      !-----------------  ProcessingQualityFlags  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_dfld%PrsQFlag)) DEALLOCATE(gds_L2CLD%cld_dfld%PrsQFlag)
      ALLOCATE (gds_L2CLD%cld_dfld%PrsQFlag(nxtrack_cld,nline_cld), STAT =ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="ProcessingQualityFlags"
              GOTO 2999
           ENDIF

      !-----------------  SlantColumnAmountNO2  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_dfld%SltColAmtNO2)) DEALLOCATE(gds_L2CLD%cld_dfld%SltColAmtNO2)
      ALLOCATE (gds_L2CLD%cld_dfld%SltColAmtNO2(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="SlantColumnAmountNO2"
              GOTO 2999
           ENDIF

      !-----------------  SlantColumnAmountO2O2  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_dfld%SltColAmtO2O2)) DEALLOCATE(gds_L2CLD%cld_dfld%SltColAmtO2O2)
      ALLOCATE (gds_L2CLD%cld_dfld%SltColAmtO2O2(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="SlantColumnAmountO2O2"
              GOTO 2999
           ENDIF

      !-----------------  SlantColumnAmountO3  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_dfld%SltColAmtO3)) DEALLOCATE(gds_L2CLD%cld_dfld%SltColAmtO3)
      ALLOCATE (gds_L2CLD%cld_dfld%SltColAmtO3(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="SlantColumnAmountO3"
              GOTO 2999
           ENDIF

      !-----------------  SnowIceExtent  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_dfld%SnwIceExt)) DEALLOCATE(gds_L2CLD%cld_dfld%SnwIceExt)
      ALLOCATE (gds_L2CLD%cld_dfld%SnwIceExt(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="SnowIceExtent"
              GOTO 2999
           ENDIF

      !-----------------  TerrainPressure  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_dfld%TerrP)) DEALLOCATE(gds_L2CLD%cld_dfld%TerrP)
      ALLOCATE (gds_L2CLD%cld_dfld%TerrP(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="TerrainPressure"
              GOTO 2999
           ENDIF

      !-----------------  TerrainReflectivity  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_dfld%TerrReflt)) DEALLOCATE(gds_L2CLD%cld_dfld%TerrReflt)
      ALLOCATE (gds_L2CLD%cld_dfld%TerrReflt(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="TerrainReflectivity"
              GOTO 2999
           ENDIF

      !--------------------------------------
      ! <Cloud Fraction And Pressure: Geolocation Fields>
      !--------------------------------------

      !-----------------  Latitude  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_gfld%LAT)) DEALLOCATE(gds_L2CLD%cld_gfld%LAT)
      ALLOCATE (gds_L2CLD%cld_gfld%LAT(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="Latitude"
              GOTO 2999
           ENDIF

      !-----------------  Longitude  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_gfld%LON)) DEALLOCATE(gds_L2CLD%cld_gfld%LON)
      ALLOCATE (gds_L2CLD%cld_gfld%LON(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="Longitude"
              GOTO 2999
           ENDIF

      !-----------------  SolarAzimuthAangle  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_gfld%SolAziAng)) DEALLOCATE(gds_L2CLD%cld_gfld%SolAziAng)
      ALLOCATE (gds_L2CLD%cld_gfld%SolAziAng(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="SolarAzimuthAangle"
              GOTO 2999
           ENDIF

      !-----------------  SolarZenithAangle  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_gfld%SolZenAng)) DEALLOCATE(gds_L2CLD%cld_gfld%SolZenAng)
      ALLOCATE (gds_L2CLD%cld_gfld%SolZenAng(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="SolarZenithAangle"
              GOTO 2999
           ENDIF

      !-----------------  SpacecraftAltitude  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_gfld%ScAlt)) DEALLOCATE(gds_L2CLD%cld_gfld%ScAlt)
      ALLOCATE (gds_L2CLD%cld_gfld%ScAlt(nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="SpacecraftAltitude"
              GOTO 2999
           ENDIF

      !-----------------  SpacecraftLatitude  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_gfld%ScLat)) DEALLOCATE(gds_L2CLD%cld_gfld%ScLat)
      ALLOCATE (gds_L2CLD%cld_gfld%ScLat(nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="SpacecraftLatitude"
              GOTO 2999
           ENDIF

      !-----------------  SpacecraftLongitude  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_gfld%ScLon)) DEALLOCATE(gds_L2CLD%cld_gfld%ScLon)
      ALLOCATE (gds_L2CLD%cld_gfld%ScLon(nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="SpacecraftLongitude"
              GOTO 2999
           ENDIF

      !-----------------  TerrainHeight  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_gfld%TerrHgt)) DEALLOCATE(gds_L2CLD%cld_gfld%TerrHgt)
      ALLOCATE (gds_L2CLD%cld_gfld%TerrHgt(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="TerrainHeight"
              GOTO 2999
           ENDIF

      !-----------------  ViewingAzimuthAngle  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_gfld%ViewAziAng)) DEALLOCATE(gds_L2CLD%cld_gfld%ViewAziAng)
      ALLOCATE (gds_L2CLD%cld_gfld%ViewAziAng(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="ViewingAzimuthAngle"
              GOTO 2999
           ENDIF

      !-----------------  ViewingZenithAngle  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_gfld%ViewZenAng)) DEALLOCATE(gds_L2CLD%cld_gfld%ViewZenAng)
      ALLOCATE (gds_L2CLD%cld_gfld%ViewZenAng(nxtrack_cld,nline_cld), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg="ViewingZenithAngle"
              GOTO 2999
           ENDIF

      gds_L2CLD%status = 1
      write(disp_msg,'("L2 CLD Memory Allocation OK !")')
      CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)

!--------------------------------
!---- END FUNCTION -------------
!--------------------------------
      CALL GEMS_Share_MOD_log(llvl, '---   END L2 CLD Memory Allocation FUNCTION ---', '  END MSG')
      CALL GEMS_Share_MOD_log(llvl, '.', LOGMSG)
      status = 0

      RETURN

2999   CONTINUE
      llvl  =9
      write(disp_msg,'( A30, "L2 CLD Memory Allocation Failed, ierr=",i6)') ctl_msg,ierr
      CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)
      status = ierr

      RETURN

      END FUNCTION GEMS_Share_MOD_L2CLD_MemAlloc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!L2 CLD Memory DeAllocation function   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Declarations:
      FUNCTION GEMS_Share_MOD_L2CLD_MemDeAlloc() RESULT (status)

! input/output parameters:

! Local parameters For L2 Reading:
      INTEGER (KIND = 4)            :: status

! Local parameters For Logging:

!--------------------------------
!---- FUNCTION Start -------------
!--------------------------------
      llvl=9
      LOGMSG=" PROC MSG"
      CALL GEMS_Share_MOD_log(llvl, '--- Start L2 CLD Memory DeAllocation FUNCTION  ---', LOGMSG )

      !-----------------  CloudFraction ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_dfld%CldFrac)) DEALLOCATE(gds_L2CLD%cld_dfld%CldFrac)

      !-----------------  CloudFractionPrecision  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_dfld%CldFracPrec)) DEALLOCATE(gds_L2CLD%cld_dfld%CldFracPrec)

      !-----------------  CloudPressure  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_dfld%CldP)) DEALLOCATE(gds_L2CLD%cld_dfld%CldP)

      !-----------------  CloudPressurePrecision  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_dfld%CldPPrec)) DEALLOCATE(gds_L2CLD%cld_dfld%CldPPrec)

      !-----------------  ContinuumAtReferenceWavelength  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_dfld%ContmAtRefWlen)) DEALLOCATE(gds_L2CLD%cld_dfld%ContmAtRefWlen)

      !----------------- ProcessingQualityFlags  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_dfld%PrsQFlag)) DEALLOCATE(gds_L2CLD%cld_dfld%PrsQFlag) 

      !-----------------  SlantColumnAmountNO2  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_dfld%SltColAmtNO2)) DEALLOCATE(gds_L2CLD%cld_dfld%SltColAmtNO2)

      !-----------------  SlantColumnAmountO2O2  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_dfld%SltColAmtO2O2)) DEALLOCATE(gds_L2CLD%cld_dfld%SltColAmtO2O2)

      !-----------------  SlantColumnAmountO3  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_dfld%SltColAmtO3)) DEALLOCATE(gds_L2CLD%cld_dfld%SltColAmtO3)

      !-----------------  SnowIceExtent  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_dfld%SnwIceExt)) DEALLOCATE(gds_L2CLD%cld_dfld%SnwIceExt)

      !-----------------  TerrainPressure  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_dfld%TerrP)) DEALLOCATE(gds_L2CLD%cld_dfld%TerrP)

      !-----------------  TerrainReflectivity  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_dfld%TerrReflt)) DEALLOCATE(gds_L2CLD%cld_dfld%TerrReflt)

      !--------------------------------------
      ! <Cloud Fraction And Pressure: Geolocation Fields>
      !--------------------------------------

      !-----------------  Latitude  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_gfld%LAT)) DEALLOCATE(gds_L2CLD%cld_gfld%LAT)

      !-----------------  Longitude  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_gfld%LON)) DEALLOCATE(gds_L2CLD%cld_gfld%LON)

      !-----------------  SolarAzimuthAangle  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_gfld%SolAziAng)) DEALLOCATE(gds_L2CLD%cld_gfld%SolAziAng)

      !-----------------  SolarZenithAangle  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_gfld%SolZenAng)) DEALLOCATE(gds_L2CLD%cld_gfld%SolZenAng)

      !-----------------  SpacecraftAltitude  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_gfld%ScAlt)) DEALLOCATE(gds_L2CLD%cld_gfld%ScAlt)

      !-----------------  SpacecraftLatitude  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_gfld%ScLat)) DEALLOCATE(gds_L2CLD%cld_gfld%ScLat)

      !-----------------  SpacecraftLongitude  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_gfld%ScLon)) DEALLOCATE(gds_L2CLD%cld_gfld%ScLon)

      !-----------------  TerrainHeight  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_gfld%TerrHgt)) DEALLOCATE(gds_L2CLD%cld_gfld%TerrHgt)

      !-----------------  ViewingAzimuthAngle  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_gfld%ViewAziAng)) DEALLOCATE(gds_L2CLD%cld_gfld%ViewAziAng)

      !-----------------  ViewingZenithAngle  ---------------
      IF (ASSOCIATED(gds_L2CLD%cld_gfld%ViewZenAng)) DEALLOCATE(gds_L2CLD%cld_gfld%ViewZenAng)

      gds_L2CLD%status = 0
      write(disp_msg,'(" L2 CLD Memory DeAllocation OK !")')
      CALL GEMS_Share_MOD_log(0, disp_msg, LOGMSG)

!--------------------------------
!---- END FUNCTION -------------
!--------------------------------
      CALL GEMS_Share_MOD_log(llvl, '---   END L2 CLD Memory DeAllocation FUNCTION ---', LOGMSG)
      CALL GEMS_Share_MOD_log(llvl, '.', LOGMSG)
      status = 0

      RETURN

      END FUNCTION GEMS_Share_MOD_L2CLD_MemDeAlloc


!!!!!!!!!!!!!!!!!!!!!!!!!!!
!L2 CLD write function2   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!

      FUNCTION GEMS_Share_MOD_L2CLD_Write2(L2CLD, ctl_fpath, wr_file_path) RESULT (status)

! input/output parameters:
      CHARACTER(LEN=128),   INTENT(INOUT) :: ctl_fpath
      CHARACTER(LEN=128),   INTENT(IN )   :: wr_file_path
      TYPE(L2_cld),         INTENT(INOUT) :: L2CLD

! Local parameters For L2 Reading:
      TYPE(CLD_ds)                  :: CLDds

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
      CALL GEMS_Share_MOD_log(llvl, '--- Start L2CLD WRITE FUNCTION2  ---', 'START MSG')

!
!--------------------------------
!
      pstart(1) = -1                                            ! HDF5 데이터셑의 모든 데이터값을 읽는다는 의미


#ifdef DEBUG
      WRITE(disp_msg, '( A25, A128, A1 )') "NameList file path     =[",trim(ctl_fpath),  "]"
      CALL GEMS_Share_MOD_log(llvl, disp_msg,  LOGMSG)
#endif
      CALL read_namelist_cld(ctl_fpath, CLDds)

      !----------      L2CLD     ----------!
      !--------------------------------------
      ! <Cloud Fraction And Pressure :Data Fields>
      !--------------------------------------

      llvl  =0
      disp_msg = "L2CLD WRITE FILE PATH="//wr_file_path
      CALL GEMS_Share_MOD_log(llvl, disp_msg,  LOGMSG)

      CALL GEMS_Share_MOD_log(llvl, "Cloud Fraction And Pressure :Data Fields",LOGMSG)

      !-----------------  CloudFraction --------------- 
      data_name = CLDds%CldFraction_DataName
      grp_path  = CLDds%cld_dgrp_path
#ifdef DEBUG
      WRITE(disp_msg, '( A25, A128, A1 )') "L2 CLD write file path =[",trim(wr_file_path)           , "]"
      CALL GEMS_Share_MOD_log(llvl, disp_msg,  LOGMSG)
      WRITE(disp_msg, '( A25, I5,   A1 )') "file path size         =[",CLDds%file_path_sz          , "]"
      CALL GEMS_Share_MOD_log(llvl, disp_msg,  LOGMSG)
      WRITE(disp_msg, '( A25, A128, A1 )') "data group path        =[",trim(grp_path)  , "]"
      CALL GEMS_Share_MOD_log(llvl, disp_msg,  LOGMSG)
      WRITE(disp_msg, '( A25, I5,   A1 )') "data group path size   =[",CLDds%grp_path_sz           , "]"
      CALL GEMS_Share_MOD_log(llvl, disp_msg,  LOGMSG)
      WRITE(disp_msg, '( A25, A128, A1 )') "dataset name           =[",trim(data_name)              , "]"
      CALL GEMS_Share_MOD_log(llvl, disp_msg,  LOGMSG)
      WRITE(disp_msg, '( A25, I5,   A1 )') "data name size         =[",CLDds%data_name_sz          , "]"
      CALL GEMS_Share_MOD_log(llvl, disp_msg,  LOGMSG)
      WRITE(disp_msg, '( A31, f10.3,   A1 )') "CloudFraction =[",SIZE(L2CLD%cld_dfld%CldFrac)  , "]"
      CALL GEMS_Share_MOD_log(llvl, disp_msg,  LOGMSG)
#endif
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_dfld%CldFrac, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_dfld%CldFrac, nxtrack_cld,nline_cld,data_name)
      ELSE
            GOTO 2999
      ENDIF
 
    !-----------------  CloudFractionPrecision  ---------------
      data_name = CLDds%CldFracPrec_DataName
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_dfld%CldFracPrec, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_dfld%CldFracPrec, nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF


      !-----------------  CloudPressure  ---------------
      data_name = CLDds%CldPressure_DataName
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_dfld%CldP, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_dfld%CldP, nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  CloudPressurePrecision  ---------------
      data_name = CLDds%CldPPrec_DataName           
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_dfld%CldPPrec, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_dfld%CldPPrec, nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  ContinuumAtReferenceWavelength  ---------------
      data_name = CLDds%ContmAtRefWlen_DataName           
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, & 
                                        data_name, CLDds%data_name_sz, L2CLD%cld_dfld%ContmAtRefWlen, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_dfld%ContmAtRefWlen,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  ProcessingQualityFlags  ---------------
      data_name = CLDds%PrsQFlag_DataName
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, CLDds%file_path_sz,grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz,L2CLD%cld_dfld%PrsQFlag, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_dfld%PrsQFlag ,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SlantColumnAmountNO2  ---------------
      data_name = CLDds%SltColAmtNO2_DataName            
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_dfld%SltColAmtNO2, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_dfld%SltColAmtNO2,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SlantColumnAmountO2O2  ---------------
      data_name = CLDds%SltColAmtO2O2_DataName
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_dfld%SltColAmtO2O2, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_dfld%SltColAmtO2O2, nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SlantColumnAmountO3  ---------------
      data_name = CLDds%SltColAmtO3_DataName
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_dfld%SltColAmtO3, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_dfld%SltColAmtO3,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SnowIceExtent  ---------------
      data_name = CLDds%SnowIceExtent_DataName           
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_dfld%SnwIceExt, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_dfld%SnwIceExt,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  TerrainPressure  ---------------
      data_name = CLDds%TerrPressure_DataName           
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_dfld%TerrP, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_dfld%TerrP, nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  TerrainReflectivity  ---------------
      data_name = CLDds%TerrReflt_DataName           
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_dfld%TerrReflt, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_dfld%TerrReflt, nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !--------------------------------------
      ! <Cloud Fraction And Pressure: Geolocation Fields>
      !--------------------------------------
      CALL GEMS_Share_MOD_log(llvl, "Cloud Fraction And Pressure: Geolocation Fields", LOGMSG)


      !-----------------  Latitude  ---------------
      data_name = CLDds%Latitude_DataName
      grp_path  = CLDds%cld_ggrp_path
           
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_gfld%LAT, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%LAT, nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  Longitude  ---------------
      data_name = CLDds%Longitude_DataName           
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_gfld%LON, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%LON, nxtrack_cld,nline_cld,data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SolarAzimuthAangle  ---------------
      data_name = CLDds%SolarAziAng_DataName          
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_gfld%SolAziAng, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%SolAziAng, nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SolarZenithAangle  ---------------
      data_name = CLDds%SolarZenAng_DataName
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_gfld%SolZenAng, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%SolZenAng,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SpacecraftAltitude  ---------------
      data_name = CLDds%ScraftAlt_DataName           
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_gfld%ScAlt, pstart, pedge)
      IF (hdferr == 0) THEN
            !CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%ScAlt, nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SpacecraftLatitude  ---------------
      data_name = CLDds%ScraftLat_DataName           
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_gfld%ScLat, pstart, pedge)
      IF (hdferr == 0) THEN
            !CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%ScLat, nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF


      !-----------------  SpacecraftLongitude  ---------------
      data_name = CLDds%ScraftLon_DataName           
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_gfld%ScLon, pstart, pedge)
      IF (hdferr == 0) THEN
            !CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%ScLon, nline_cld,data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  TerrainHeight  ---------------
      data_name = CLDds%TerrainHgt_DataName           
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_gfld%TerrHgt, pstart, pedge)
      IF (hdferr == 0) THEN
            !CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%TerrHgt, nxtrack_cld,nline_cld,data_name)
      ELSE
            GOTO 2999
      ENDIF


      !-----------------  ViewingAzimuthAngle  ---------------
      data_name = CLDds%ViewAziAng_DataName           
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_gfld%ViewAziAng, pstart, pedge)
      IF (hdferr == 0) THEN
            !CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%ViewAziAng,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  ViewingZenithAngle  ---------------
      data_name = CLDds%ViewZenAng_DataName           
      hdferr = GEMS_Share_Hdf5WriteData(wr_file_path, CLDds%file_path_sz, grp_path, CLDds%grp_path_sz, &
                                        data_name, CLDds%data_name_sz, L2CLD%cld_gfld%ViewZenAng, pstart, pedge)
      IF (hdferr == 0) THEN
            !CALL GEMS_Share_MOD_log(llvl, L2CLD%cld_gfld%ViewZenAng,nxtrack_cld,nline_cld, data_name)
      ELSE
            GOTO 2999
      ENDIF

      write(disp_msg,'(" GEMS HDF5 File WRITE OK, hdferr=",i6)') hdferr
      CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)

!--------------------------------
!---- END FUNCTION -------------
!--------------------------------
      llvl  =9
      CALL GEMS_Share_MOD_log(llvl, '---   END L2 WRITE FUNCTION ---', '  END MSG')
      CALL GEMS_Share_MOD_log(llvl, '.', LOGMSG)
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

      END FUNCTION GEMS_Share_MOD_L2CLD_Write2

!-------------------------------------------------------------------------------
!+Description: 
!      FUNCTION to initialize L2 CLD Memory.
!
!
! Method:
!
! input  :
!
! input/output :
!
! output :
!       status    : the result code of L2CLD Memory Initialization
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2015.04.03  Fisrt Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!L2 CLD Memory Allocation function  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GEMS_Share_MOD_L2CLD_MemInit() RESULT (status)

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
      CALL GEMS_Share_MOD_log(llvl, '--- Start L2 CLD Memory Initialization FUNCTION  ---', LOGMSG)

      IF ( gds_L2CLD%status == 0 ) THEN
          llvl    =0
          disp_msg="did not allocate L2 CLD Memory"
          CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)
          status = -7001  ! MEMOEY ALLOCATION ERROR CODE = -7001
          RETURN
      END IF

      !--------------------------------------
      ! <Cloud Fraction And Pressure :Data Fields>
      !--------------------------------------
      !-----------------  CloudFraction ---------------
      gds_L2CLD%cld_dfld%CldFrac = DEF_fVAL

      !-----------------  CloudFractionPrecision  ---------------
      gds_L2CLD%cld_dfld%CldFracPrec = DEF_fVAL

      !-----------------  CloudPressure  ---------------
      gds_L2CLD%cld_dfld%CldP = DEF_fVAL

      !-----------------  CloudPressurePrecision  ---------------
      gds_L2CLD%cld_dfld%CldPPrec = DEF_fVAL

      !-----------------  ContinuumAtReferenceWavelength  ---------------
      gds_L2CLD%cld_dfld%ContmAtRefWlen = DEF_dVAL

      !-----------------  ProcessingQualityFlags  ---------------
      gds_L2CLD%cld_dfld%PrsQFlag = DEF_sVAL

      !-----------------  SlantColumnAmountNO2  ---------------
      gds_L2CLD%cld_dfld%SltColAmtNO2 = DEF_fVAL

      !-----------------  SlantColumnAmountO2O2  ---------------
      gds_L2CLD%cld_dfld%SltColAmtO2O2 = DEF_fVAL

      !-----------------  SlantColumnAmountO3  ---------------
      gds_L2CLD%cld_dfld%SltColAmtO3 = DEF_fVAL

      !-----------------  SnowIceExtent  ---------------
      gds_L2CLD%cld_dfld%SnwIceExt = DEF_cVAL

      !-----------------  TerrainPressure  ---------------
      gds_L2CLD%cld_dfld%TerrP = DEF_fVAL

      !-----------------  TerrainReflectivity  ---------------
      gds_L2CLD%cld_dfld%TerrReflt = DEF_fVAL

      !--------------------------------------
      ! <Cloud Fraction And Pressure: Geolocation Fields>
      !--------------------------------------

      !-----------------  Latitude  ---------------
      gds_L2CLD%cld_gfld%LAT = DEF_fVAL

      !-----------------  Longitude  ---------------
      gds_L2CLD%cld_gfld%LON = DEF_fVAL 

      !-----------------  SolarAzimuthAangle  ---------------
      gds_L2CLD%cld_gfld%SolAziAng = DEF_fVAL

      !-----------------  SolarZenithAangle  ---------------
      gds_L2CLD%cld_gfld%SolZenAng = DEF_fVAL

      !-----------------  SpacecraftAltitude  ---------------
      gds_L2CLD%cld_gfld%ScAlt = DEF_fVAL

      !-----------------  SpacecraftLatitude  ---------------
      gds_L2CLD%cld_gfld%ScLat = DEF_fVAL

      !-----------------  SpacecraftLongitude  ---------------
      gds_L2CLD%cld_gfld%ScLon = DEF_fVAL

      !-----------------  TerrainHeight  ---------------
      gds_L2CLD%cld_gfld%TerrHgt = DEF_fVAL

      !-----------------  ViewingAzimuthAngle  ---------------
      gds_L2CLD%cld_gfld%ViewAziAng = DEF_fVAL

      !-----------------  ViewingZenithAngle  ---------------
      gds_L2CLD%cld_gfld%ViewZenAng = DEF_fVAL


      write(disp_msg,'("L2 CLD Memory Initialization OK !")')
      CALL GEMS_Share_MOD_log(0, disp_msg, LOGMSG)

      status = 0

      RETURN

      END FUNCTION GEMS_Share_MOD_L2CLD_MemInit

!-------------------------------------------------------------------------------
!+Description: 
!      FUNCTION ver2 to allocate L2 CLD Memory.
!
!
! Method:
!
! input  :
!
! input/output :
!     local_L2CLD : L2 CLD Data
!
! output :
!       status    : the result code of L2CLD Memory Initialization
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2016.12.13  전역변수를 매개변수로 사용하기 위해 생성(YuGeun Ki) 
!-------------------------------------------------------------------------------
! Declarations:
FUNCTION GEMS_Share_MOD_L2CLD_MemAlloc2(local_L2CLD) RESULT (status)

! input/output parameters:
    TYPE(L2_cld), INTENT(INOUT)     :: local_L2CLD

! Local parameters For L2 Reading:
    INTEGER                         :: ierr        ! Error code 
    INTEGER (KIND = 4)              :: status

! Local parameters For Logging:


!--------------------------------
!---- FUNCTION Start -------------
!--------------------------------
    llvl  =9
    LOGMSG=" PROC MSG"
    CALL GEMS_Share_MOD_log(llvl, '--- Start L2 CLD Memory Allocation FUNCTION ver2 ---', LOGMSG )

    !-----------------  CloudFraction ---------------
    IF (ASSOCIATED(local_L2CLD%cld_dfld%CldFrac)) DEALLOCATE(local_L2CLD%cld_dfld%CldFrac)
    ALLOCATE (local_L2CLD%cld_dfld%CldFrac(nxtrack_cld,nline_cld), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg="CloudFraction"
            GOTO 2999
         ENDIF

    !-----------------  CloudFractionPrecision  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_dfld%CldFracPrec)) DEALLOCATE(local_L2CLD%cld_dfld%CldFracPrec)
    ALLOCATE (local_L2CLD%cld_dfld%CldFracPrec(nxtrack_cld,nline_cld), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg="CloudFractionPrecision"
            GOTO 2999
         ENDIF

    !-----------------  CloudPressure  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_dfld%CldP)) DEALLOCATE(local_L2CLD%cld_dfld%CldP)
    ALLOCATE (local_L2CLD%cld_dfld%CldP(nxtrack_cld,nline_cld), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg="CloudPressure"
            GOTO 2999
         ENDIF

    !-----------------  CloudPressurePrecision  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_dfld%CldPPrec)) DEALLOCATE(local_L2CLD%cld_dfld%CldPPrec)
    ALLOCATE (local_L2CLD%cld_dfld%CldPPrec(nxtrack_cld,nline_cld), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg="CloudPressurePrecision"
            GOTO 2999
         ENDIF

    !-----------------  ContinuumAtReferenceWavelength  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_dfld%ContmAtRefWlen)) DEALLOCATE(local_L2CLD%cld_dfld%ContmAtRefWlen)
    ALLOCATE (local_L2CLD%cld_dfld%ContmAtRefWlen(nxtrack_cld,nline_cld), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg="ContinuumAtReferenceWavelength"
            GOTO 2999
         ENDIF

    !-----------------  ProcessingQualityFlags  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_dfld%PrsQFlag)) DEALLOCATE(local_L2CLD%cld_dfld%PrsQFlag)
    ALLOCATE (local_L2CLD%cld_dfld%PrsQFlag(nxtrack_cld,nline_cld), STAT =ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg="ProcessingQualityFlags"
            GOTO 2999
         ENDIF

    !-----------------  SlantColumnAmountNO2  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_dfld%SltColAmtNO2)) DEALLOCATE(local_L2CLD%cld_dfld%SltColAmtNO2)
    ALLOCATE (local_L2CLD%cld_dfld%SltColAmtNO2(nxtrack_cld,nline_cld), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg="SlantColumnAmountNO2"
            GOTO 2999
         ENDIF

    !-----------------  SlantColumnAmountO2O2  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_dfld%SltColAmtO2O2)) DEALLOCATE(local_L2CLD%cld_dfld%SltColAmtO2O2)
    ALLOCATE (local_L2CLD%cld_dfld%SltColAmtO2O2(nxtrack_cld,nline_cld), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg="SlantColumnAmountO2O2"
            GOTO 2999
         ENDIF

    !-----------------  SlantColumnAmountO3  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_dfld%SltColAmtO3)) DEALLOCATE(local_L2CLD%cld_dfld%SltColAmtO3)
    ALLOCATE (local_L2CLD%cld_dfld%SltColAmtO3(nxtrack_cld,nline_cld), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg="SlantColumnAmountO3"
            GOTO 2999
         ENDIF

    !-----------------  SnowIceExtent  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_dfld%SnwIceExt)) DEALLOCATE(local_L2CLD%cld_dfld%SnwIceExt)
    ALLOCATE (local_L2CLD%cld_dfld%SnwIceExt(nxtrack_cld,nline_cld), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg="SnowIceExtent"
            GOTO 2999
         ENDIF

    !-----------------  TerrainPressure  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_dfld%TerrP)) DEALLOCATE(local_L2CLD%cld_dfld%TerrP)
    ALLOCATE (local_L2CLD%cld_dfld%TerrP(nxtrack_cld,nline_cld), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg="TerrainPressure"
            GOTO 2999
         ENDIF

    !-----------------  TerrainReflectivity  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_dfld%TerrReflt)) DEALLOCATE(local_L2CLD%cld_dfld%TerrReflt)
    ALLOCATE (local_L2CLD%cld_dfld%TerrReflt(nxtrack_cld,nline_cld), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg="TerrainReflectivity"
            GOTO 2999
         ENDIF

    !--------------------------------------
    ! <Cloud Fraction And Pressure: Geolocation Fields>
    !--------------------------------------

    !-----------------  Latitude  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_gfld%LAT)) DEALLOCATE(local_L2CLD%cld_gfld%LAT)
    ALLOCATE (local_L2CLD%cld_gfld%LAT(nxtrack_cld,nline_cld), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg="Latitude"
            GOTO 2999
         ENDIF

    !-----------------  Longitude  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_gfld%LON)) DEALLOCATE(local_L2CLD%cld_gfld%LON)
    ALLOCATE (local_L2CLD%cld_gfld%LON(nxtrack_cld,nline_cld), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg="Longitude"
            GOTO 2999
         ENDIF

    !-----------------  SolarAzimuthAangle  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_gfld%SolAziAng)) DEALLOCATE(local_L2CLD%cld_gfld%SolAziAng)
    ALLOCATE (local_L2CLD%cld_gfld%SolAziAng(nxtrack_cld,nline_cld), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg="SolarAzimuthAangle"
            GOTO 2999
         ENDIF

    !-----------------  SolarZenithAangle  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_gfld%SolZenAng)) DEALLOCATE(local_L2CLD%cld_gfld%SolZenAng)
    ALLOCATE (local_L2CLD%cld_gfld%SolZenAng(nxtrack_cld,nline_cld), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg="SolarZenithAangle"
            GOTO 2999
         ENDIF

    !-----------------  SpacecraftAltitude  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_gfld%ScAlt)) DEALLOCATE(local_L2CLD%cld_gfld%ScAlt)
    ALLOCATE (local_L2CLD%cld_gfld%ScAlt(nline_cld), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg="SpacecraftAltitude"
            GOTO 2999
         ENDIF

    !-----------------  SpacecraftLatitude  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_gfld%ScLat)) DEALLOCATE(local_L2CLD%cld_gfld%ScLat)
    ALLOCATE (local_L2CLD%cld_gfld%ScLat(nline_cld), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg="SpacecraftLatitude"
            GOTO 2999
         ENDIF

    !-----------------  SpacecraftLongitude  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_gfld%ScLon)) DEALLOCATE(local_L2CLD%cld_gfld%ScLon)
    ALLOCATE (local_L2CLD%cld_gfld%ScLon(nline_cld), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg="SpacecraftLongitude"
            GOTO 2999
         ENDIF

    !-----------------  TerrainHeight  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_gfld%TerrHgt)) DEALLOCATE(local_L2CLD%cld_gfld%TerrHgt)
    ALLOCATE (local_L2CLD%cld_gfld%TerrHgt(nxtrack_cld,nline_cld), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg="TerrainHeight"
            GOTO 2999
         ENDIF

    !-----------------  ViewingAzimuthAngle  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_gfld%ViewAziAng)) DEALLOCATE(local_L2CLD%cld_gfld%ViewAziAng)
    ALLOCATE (local_L2CLD%cld_gfld%ViewAziAng(nxtrack_cld,nline_cld), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg="ViewingAzimuthAngle"
            GOTO 2999
         ENDIF

    !-----------------  ViewingZenithAngle  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_gfld%ViewZenAng)) DEALLOCATE(local_L2CLD%cld_gfld%ViewZenAng)
    ALLOCATE (local_L2CLD%cld_gfld%ViewZenAng(nxtrack_cld,nline_cld), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg="ViewingZenithAngle"
            GOTO 2999
         ENDIF

    local_L2CLD%status = 1
    write(disp_msg,'("L2 CLD Memory Allocation ver2 OK !")')
    CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)

!--------------------------------
!---- END FUNCTION -------------
!--------------------------------
    CALL GEMS_Share_MOD_log(llvl, '---   END L2 CLD Memory Allocation FUNCTION ver2 ---', '  END MSG')
    CALL GEMS_Share_MOD_log(llvl, '.', LOGMSG)
    status = 0

    RETURN

2999   CONTINUE
    llvl  =9
    write(disp_msg,'( A30, "L2 CLD Memory Allocation ver2 Failed, ierr=",i6)') ctl_msg,ierr
    CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)
    status = ierr

    RETURN

END FUNCTION GEMS_Share_MOD_L2CLD_MemAlloc2


!-------------------------------------------------------------------------------
!+Description: 
!      FUNCTION ver2 to deallocate L2 CLD Memory.
!
!
! Method:
!
! input  :
!
! input/output :
!     local_L2CLD : L2 CLD Data
!
! output :
!       status    : the result code of L2CLD Memory DeAllocation
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2016.12.13  전역변수를 매개변수로 사용하기 위해 생성(YG Ki) 
!-------------------------------------------------------------------------------
! Declarations:
FUNCTION GEMS_Share_MOD_L2CLD_MemDeAlloc2(local_L2CLD) RESULT (status)

! input/output parameters:
    TYPE(L2_cld), INTENT(INOUT)     :: local_L2CLD

! Local parameters For L2 Reading:
    INTEGER (KIND = 4)              :: status

! Local parameters For Logging:

!--------------------------------
!---- FUNCTION Start -------------
!--------------------------------
    llvl=9
    LOGMSG=" PROC MSG"
    CALL GEMS_Share_MOD_log(llvl, '--- Start L2 CLD Memory DeAllocation FUNCTION ver2 ---', LOGMSG )

    !-----------------  CloudFraction ---------------
    IF (ASSOCIATED(local_L2CLD%cld_dfld%CldFrac)) DEALLOCATE(local_L2CLD%cld_dfld%CldFrac)

    !-----------------  CloudFractionPrecision  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_dfld%CldFracPrec)) DEALLOCATE(local_L2CLD%cld_dfld%CldFracPrec)

    !-----------------  CloudPressure  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_dfld%CldP)) DEALLOCATE(local_L2CLD%cld_dfld%CldP)

    !-----------------  CloudPressurePrecision  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_dfld%CldPPrec)) DEALLOCATE(local_L2CLD%cld_dfld%CldPPrec)

    !-----------------  ContinuumAtReferenceWavelength  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_dfld%ContmAtRefWlen)) DEALLOCATE(local_L2CLD%cld_dfld%ContmAtRefWlen)

    !----------------- ProcessingQualityFlags  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_dfld%PrsQFlag)) DEALLOCATE(local_L2CLD%cld_dfld%PrsQFlag) 

    !-----------------  SlantColumnAmountNO2  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_dfld%SltColAmtNO2)) DEALLOCATE(local_L2CLD%cld_dfld%SltColAmtNO2)

    !-----------------  SlantColumnAmountO2O2  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_dfld%SltColAmtO2O2)) DEALLOCATE(local_L2CLD%cld_dfld%SltColAmtO2O2)

    !-----------------  SlantColumnAmountO3  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_dfld%SltColAmtO3)) DEALLOCATE(local_L2CLD%cld_dfld%SltColAmtO3)

    !-----------------  SnowIceExtent  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_dfld%SnwIceExt)) DEALLOCATE(local_L2CLD%cld_dfld%SnwIceExt)

    !-----------------  TerrainPressure  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_dfld%TerrP)) DEALLOCATE(local_L2CLD%cld_dfld%TerrP)

    !-----------------  TerrainReflectivity  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_dfld%TerrReflt)) DEALLOCATE(local_L2CLD%cld_dfld%TerrReflt)

    !--------------------------------------
    ! <Cloud Fraction And Pressure: Geolocation Fields>
    !--------------------------------------

    !-----------------  Latitude  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_gfld%LAT)) DEALLOCATE(local_L2CLD%cld_gfld%LAT)

    !-----------------  Longitude  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_gfld%LON)) DEALLOCATE(local_L2CLD%cld_gfld%LON)

    !-----------------  SolarAzimuthAangle  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_gfld%SolAziAng)) DEALLOCATE(local_L2CLD%cld_gfld%SolAziAng)

    !-----------------  SolarZenithAangle  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_gfld%SolZenAng)) DEALLOCATE(local_L2CLD%cld_gfld%SolZenAng)

    !-----------------  SpacecraftAltitude  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_gfld%ScAlt)) DEALLOCATE(local_L2CLD%cld_gfld%ScAlt)

    !-----------------  SpacecraftLatitude  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_gfld%ScLat)) DEALLOCATE(local_L2CLD%cld_gfld%ScLat)

    !-----------------  SpacecraftLongitude  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_gfld%ScLon)) DEALLOCATE(local_L2CLD%cld_gfld%ScLon)

    !-----------------  TerrainHeight  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_gfld%TerrHgt)) DEALLOCATE(local_L2CLD%cld_gfld%TerrHgt)

    !-----------------  ViewingAzimuthAngle  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_gfld%ViewAziAng)) DEALLOCATE(local_L2CLD%cld_gfld%ViewAziAng)

    !-----------------  ViewingZenithAngle  ---------------
    IF (ASSOCIATED(local_L2CLD%cld_gfld%ViewZenAng)) DEALLOCATE(local_L2CLD%cld_gfld%ViewZenAng)

    local_L2CLD%status = 0
    write(disp_msg,'(" L2 CLD Memory DeAllocation ver2 OK !")')
    CALL GEMS_Share_MOD_log(0, disp_msg, LOGMSG)

!--------------------------------
!---- END FUNCTION -------------
!--------------------------------
    CALL GEMS_Share_MOD_log(llvl, '---   END L2 CLD Memory DeAllocation FUNCTION ver2 ---', LOGMSG)
    CALL GEMS_Share_MOD_log(llvl, '.', LOGMSG)
    status = 0

    RETURN

END FUNCTION GEMS_Share_MOD_L2CLD_MemDeAlloc2

!-------------------------------------------------------------------------------
!+Description: 
!      FUNCTION ver2 to initialize L2 CLD Memory.
!
!
! Method:
!
! input  :
!
! input/output :
!     local_L2CLD : L2 CLD Data
!
! output :
!       status    : the result code of L2CLD Memory Initialization
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2016.12.13  전역변수를 매개변수로 사용하기 위해 생성(YG Ki) 
!-------------------------------------------------------------------------------
FUNCTION GEMS_Share_MOD_L2CLD_MemInit2(local_L2CLD) RESULT (status)

! input/output parameters:
    TYPE(L2_cld), INTENT(INOUT) :: local_L2CLD

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
    CALL GEMS_Share_MOD_log(llvl, '--- Start L2 CLD Memory Initialization FUNCTION ver2 ---', LOGMSG)

    IF ( local_L2CLD%status == 0 ) THEN
        llvl    =0
        disp_msg="did not allocate L2 CLD Memory"
        CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)
        status = -7001  ! MEMOEY ALLOCATION ERROR CODE = -7001
        RETURN
    END IF

    !--------------------------------------
    ! <Cloud Fraction And Pressure :Data Fields>
    !--------------------------------------
    !-----------------  CloudFraction ---------------
    local_L2CLD%cld_dfld%CldFrac = DEF_fVAL

    !-----------------  CloudFractionPrecision  ---------------
    local_L2CLD%cld_dfld%CldFracPrec = DEF_fVAL

    !-----------------  CloudPressure  ---------------
    local_L2CLD%cld_dfld%CldP = DEF_fVAL

    !-----------------  CloudPressurePrecision  ---------------
    local_L2CLD%cld_dfld%CldPPrec = DEF_fVAL

    !-----------------  ContinuumAtReferenceWavelength  ---------------
    local_L2CLD%cld_dfld%ContmAtRefWlen = DEF_dVAL

    !-----------------  ProcessingQualityFlags  ---------------
    local_L2CLD%cld_dfld%PrsQFlag = DEF_sVAL

    !-----------------  SlantColumnAmountNO2  ---------------
    local_L2CLD%cld_dfld%SltColAmtNO2 = DEF_fVAL

    !-----------------  SlantColumnAmountO2O2  ---------------
    local_L2CLD%cld_dfld%SltColAmtO2O2 = DEF_fVAL

    !-----------------  SlantColumnAmountO3  ---------------
    local_L2CLD%cld_dfld%SltColAmtO3 = DEF_fVAL

    !-----------------  SnowIceExtent  ---------------
    local_L2CLD%cld_dfld%SnwIceExt = DEF_cVAL

    !-----------------  TerrainPressure  ---------------
    local_L2CLD%cld_dfld%TerrP = DEF_fVAL

    !-----------------  TerrainReflectivity  ---------------
    local_L2CLD%cld_dfld%TerrReflt = DEF_fVAL

    !--------------------------------------
    ! <Cloud Fraction And Pressure: Geolocation Fields>
    !--------------------------------------

    !-----------------  Latitude  ---------------
    local_L2CLD%cld_gfld%LAT = DEF_fVAL

    !-----------------  Longitude  ---------------
    local_L2CLD%cld_gfld%LON = DEF_fVAL 

    !-----------------  SolarAzimuthAangle  ---------------
    local_L2CLD%cld_gfld%SolAziAng = DEF_fVAL

    !-----------------  SolarZenithAangle  ---------------
    local_L2CLD%cld_gfld%SolZenAng = DEF_fVAL

    !-----------------  SpacecraftAltitude  ---------------
    local_L2CLD%cld_gfld%ScAlt = DEF_fVAL

    !-----------------  SpacecraftLatitude  ---------------
    local_L2CLD%cld_gfld%ScLat = DEF_fVAL

    !-----------------  SpacecraftLongitude  ---------------
    local_L2CLD%cld_gfld%ScLon = DEF_fVAL

    !-----------------  TerrainHeight  ---------------
    local_L2CLD%cld_gfld%TerrHgt = DEF_fVAL

    !-----------------  ViewingAzimuthAngle  ---------------
    local_L2CLD%cld_gfld%ViewAziAng = DEF_fVAL

    !-----------------  ViewingZenithAngle  ---------------
    local_L2CLD%cld_gfld%ViewZenAng = DEF_fVAL


    write(disp_msg,'("L2 CLD Memory Initialization ver2 OK !")')
    CALL GEMS_Share_MOD_log(0, disp_msg, LOGMSG)

    status = 0

    RETURN

END FUNCTION GEMS_Share_MOD_L2CLD_MemInit2

END MODULE Share_l2_cld_mod_write_read

!-------------------------------------------------------------------------------
