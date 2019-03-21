
!+Module to read GEMS L1BFile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE Share_MOD_L1B_Read

!-------------------------------------------------------------------------------
!+Description: 
!     Module to read GEMS L1BFile
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1.0   2014.12.11 Fisrt Code (HJ Lee, Seasoft) 
! 0.1.1   2016.09.20 L1BIRR/Sun Volume VIS Swath/Data Fields/
!                    PixelQualityFlags HDF5 Read 함수 추가     (YG Ki, Seasoft) 
!-------------------------------------------------------------------------------

! Modules used:
      USE Share_MOD_Constants_Variables, ONLY:  gci_nline      ,   & 
                                                gci_nline_rug  ,   & 
                                                gci_nline_rvg  ,   & 
                                                gci_nline_irr  ,   & 
                                                gci_nwavel     ,   & 
                                                gci_nwavel2    ,   & 
                                                gci_nwavel3    ,   & 
                                                gci_nxtrack    ,   & 
                                                gci_nxtrack2   ,   & 
                                                gci_nwavelcoef
      USE Share_MOD_P_Constants_Variables, ONLY: LOGMSG, llvl, ctl_msg, &
                                                 disp_msg, disp_line
      USE Share_MOD_Log

!**********************************************************!
      IMPLICIT NONE

      PUBLIC  :: GEMS_Share_MOD_L1B_Read
      PUBLIC  :: GEMS_Share_MOD_L1Brug_Read
      PUBLIC  :: GEMS_Share_MOD_L1Brvg_Read
      PUBLIC  :: GEMS_Share_MOD_L1Birr_Read
      
      INTEGER(KIND=4), PRIVATE      :: nLine_L1B          
      INTEGER(KIND=4), PRIVATE      :: nLine_L1B_rug      
      INTEGER(KIND=4), PRIVATE      :: nLine_L1B_rvg      
      INTEGER(KIND=4), PRIVATE      :: nLine_L1B_irr      
      INTEGER(KIND=4), PRIVATE      :: nWavel_L1B         
      INTEGER(KIND=4), PRIVATE      :: nWavel2_L1B        
      INTEGER(KIND=4), PRIVATE      :: nWavel3_L1B        
      INTEGER(KIND=4), PRIVATE      :: nxtrack_L1B        
      INTEGER(KIND=4), PRIVATE      :: nxtrack2_L1B       
      INTEGER(KIND=4), PRIVATE      :: nWavelCoef_L1B     

      INTEGER(KIND=4), PARAMETER    :: PATH_SZ = 128     
      !
      ! --- filepath and namelist needed to read L1B File 
      !
      CHARACTER(LEN=PATH_SZ)        :: gcc_l1b_ctl_file = "../../../share/conf/l1breadmdl.nml"    

      TYPE :: l1b_file_path
        CHARACTER(LEN=PATH_SZ)      :: rug
        CHARACTER(LEN=PATH_SZ)      :: rvg
        CHARACTER(LEN=PATH_SZ)      :: irr
      END TYPE

      TYPE(l1b_file_path)           :: gds_l1b_fpath

     !----------L1BRUG----------!
     ! <Earth UV-1 Swath and Earth UV-2 Swath>

      TYPE :: brug_d_fld
       REAL   (KIND=4) , POINTER    :: ExpTime      (:)         ! ExposureTime
       INTEGER(KIND=1) , POINTER    :: ExpType      (:)         ! ExposureType
       INTEGER(KIND=2) , POINTER    :: MQFlag       (:)         ! MeasurementQualityFlags
       INTEGER(KIND=1) , POINTER    :: NumSmPixCol  (:)         ! NumberSmallPixelColumns
       INTEGER(KIND=2) , POINTER    :: PQFlag       (:,:,:)     ! PixelQualityFlags
       INTEGER(KIND=1) , POINTER    :: RadExponent  (:,:,:)     ! RadianceExponent
       INTEGER(KIND=2) , POINTER    :: RadMantissa  (:,:,:)     ! RadianceMantissa
       INTEGER(KIND=2) , POINTER    :: RadPrecision (:,:,:)     ! RadiancePrecisionMantissa
       REAL   (KIND=4) , POINTER    :: ReadTime     (:)         ! ReadoutTime
       REAL   (KIND=4) , POINTER    :: WlenCoef     (:,:,:)     ! WavelengthCoefficient
       REAL   (KIND=4) , POINTER    :: WlenPrec     (:,:,:)     ! WavelengthCoefficientPrecision
       INTEGER(KIND=2) , POINTER    :: WlenRefCol   (:)         ! WavelengthReferenceColumn
      END TYPE brug_d_fld

      TYPE :: brug_g_fld
       INTEGER(KIND=2) , POINTER    :: GPQFlag      (:,:)       ! GroundPixelQualityFlags
       REAL   (KIND=4) , POINTER    :: LAT          (:,:)       ! Latitude
       REAL   (KIND=4) , POINTER    :: LON          (:,:)       ! Longitude
       REAL   (KIND=4) , POINTER    :: SecInDay     (:)         ! SecondsInDay
       REAL   (KIND=4) , POINTER    :: SolAziAng    (:,:)       ! SolarAzimuthAngle
       REAL   (KIND=4) , POINTER    :: SolZenAng    (:,:)       ! SolarZenithAngle
       REAL   (KIND=4) , POINTER    :: ScAlt        (:)         ! SpacecraftAltitude
       REAL   (KIND=4) , POINTER    :: ScLat        (:)         ! SpacecraftLatitude
       REAL   (KIND=4) , POINTER    :: ScLon        (:)         ! SpacecraftLongitude
       INTEGER(KIND=2) , POINTER    :: TerrHgt      (:,:)       ! TerrainHeight      
       REAL   (KIND=8) , POINTER    :: Time         (:)         ! Time
       REAL   (KIND=4) , POINTER    :: ViewAziAng   (:,:)       ! ViewingAzimuthAngle 
       REAL   (KIND=4) , POINTER    :: ViewZenAng   (:,:)       ! ViewingZenithAngle
       INTEGER(KIND=1) , POINTER    :: XTQFlag      (:,:)       ! XTrackQualityFlags 
      END TYPE brug_g_fld

      TYPE :: brug_sw_att
       REAL   (KIND=4)              :: ES_Dist                  ! EarthSunDistance
      END TYPE brug_sw_att

      TYPE :: e_uv
       TYPE(brug_d_fld)             :: brug_dfld
       TYPE(brug_g_fld)             :: brug_gfld
       TYPE(brug_sw_att)            :: brug_swatt
      END TYPE e_uv

      TYPE :: L1B_rug
       TYPE(e_uv)                   :: euv1
       TYPE(e_uv)                   :: euv2
      END TYPE L1B_rug

      TYPE(L1B_rug)                 :: gds_brug


     !----------L1BRVG----------!
     ! <Earth VIS Swath>
      
      !--- Data Fields --- 
      TYPE :: brvg_d_fld
       INTEGER(KIND=2) , POINTER    :: PQFlag        (:,:,:)    ! PixelQualityFlags
       INTEGER(KIND=1) , POINTER    :: RadExponent   (:,:,:)    ! RadianceExponent
       INTEGER(KIND=2) , POINTER    :: RadMantissa   (:,:,:)    ! RadianceMantissa
       REAL   (KIND=4) , POINTER    :: WlenCoef      (:,:,:)    ! WavelengthCoefficient
       INTEGER(KIND=2) , POINTER    :: WlenRefCol    (:)        ! WavelengthReferenceColumn
      END TYPE brvg_d_fld

      !--- Geolocation Fields --- 
      TYPE :: brvg_g_fld
       REAL   (KIND=4) , POINTER    :: LAT           (:,:)      ! Longitude           
       REAL   (KIND=4) , POINTER    :: LON           (:,:)      ! Latitude            
       REAL   (KIND=4) , POINTER    :: SolAziAng     (:,:)      ! SolarAzimuthAngle 
       REAL   (KIND=4) , POINTER    :: SolZenAng     (:,:)      ! SolarZenithAngle   
       INTEGER(KIND=2) , POINTER    :: TerrHgt       (:,:)      ! TerrainHeight
       REAL   (KIND=4) , POINTER    :: ViewZenAng    (:,:)      ! ViewingZenithAngle   
       REAL   (KIND=4) , POINTER    :: ViewAziAng    (:,:)      ! ViewingAzimuthAngle 
      END TYPE brvg_g_fld
     
      TYPE :: brvg_sw_att
       REAL   (KIND=4)              :: ES_Dist                  ! EarthSunDistance
      END TYPE brvg_sw_att
 
      TYPE :: L1B_rvg
       TYPE(brvg_d_fld)             :: brvg_dfld
       TYPE(brvg_g_fld)             :: brvg_gfld
       TYPE(brvg_sw_att)            :: brvg_swatt
      END TYPE L1B_rvg

      TYPE(L1B_rvg)                 ::  gds_brvg

     !----------L1BIRR----------!
     ! <Sun Volume UV-2 Swath and Sun Volume VIS Swath >

      !--- Data Fields --- 
      TYPE :: birr_d_fld
       INTEGER(KIND=1) , POINTER    :: IrrRadExponent (:,:)     ! IrradianceExponent 
       INTEGER(KIND=2) , POINTER    :: IrrRadMantissa (:,:)     ! IrradianceMantissa 
       INTEGER(KIND=2) , POINTER    :: IrrRadPrecision(:,:)     ! IrradiancePrecisionMantissa
       INTEGER(KIND=2) , POINTER    :: PQFlag         (:,:)     ! PixelQualityFlags 
       REAL   (KIND=4) , POINTER    :: WlenCoef       (:,:)     ! WavelengthCoefficient 
       INTEGER(KIND=2)              :: WlenRefCol               ! WavelengthReferenceColumn 
      END TYPE birr_d_fld

      !--- Swath Attributes --- 
      TYPE :: birr_sw_att
       REAL   (KIND=4)              :: ES_Dist                   ! EarthSunDistance 
      END TYPE birr_sw_att

      TYPE :: s_uv
       TYPE(birr_d_fld)             :: birr_dfld
       TYPE(birr_sw_att)            :: birr_swatt   
      END TYPE s_uv

      TYPE :: L1B_irr
       TYPE(s_uv)                   :: suv1
       TYPE(s_uv)                   :: suv2
       TYPE(s_uv)                   :: svis
      END TYPE L1B_irr

      TYPE(L1B_irr)                 :: gds_birr

      CONTAINS


!-------------------------------------------------------------------------------
!+Description: 
!      FUNCTION GEMS_Share_MOD_L1B_Read .
!
!
! Method:
!
!
! Input files:
!       L1Brug    : the Level 1B RUG data structure
!       L1Brvg    : the Level 1B RVG data structure
!       L1Birr    : the Level 1B IRR data structure
!       ctl_fpath : the control file path
!
! Output files:
!       status    : the result code of L1B read
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2014.12.23  Fisrt Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------

      
! Declarations:

      FUNCTION GEMS_Share_MOD_L1B_Read(L1Brug, L1Brvg, L1Birr, ctl_fpath) RESULT (status)

        TYPE(L1B_rug),      INTENT(INOUT)     :: L1Brug
        TYPE(L1B_rvg),      INTENT(INOUT)     :: L1Brvg
        TYPE(L1B_irr),      INTENT(INOUT)     :: L1Birr
        CHARACTER(LEN=128), INTENT(INOUT)     :: ctl_fpath

        INTEGER (KIND = 4)                    :: status

        status = GEMS_Share_MOD_L1Brug_Read(L1Brug, ctl_fpath)
        if ( status .ne. 0 ) GOTO 1999

        status = GEMS_Share_MOD_L1Brvg_Read(L1Brvg, ctl_fpath)
        if ( status .ne. 0 ) GOTO 1999

        status = GEMS_Share_MOD_L1Birr_Read(L1Birr, ctl_fpath)
        if ( status .ne. 0 ) GOTO 1999

1999    CONTINUE
        RETURN 

      END FUNCTION GEMS_Share_MOD_L1B_Read

!-------------------------------------------------------------------------------
!+Description: 
!      FUNCTION GEMS_Share_MOD_L1Brug_Read .
!
!
! Method:
!
!
! Input files:
!       L1Brug    : the Level 1B RUG data structure
!       ctl_fpath : the control file path
!
! Output files:
!       status    : the result code of L1B RUG read
!
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2014.12.23  Fisrt Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------

      
! Declarations:

      FUNCTION GEMS_Share_MOD_L1Brug_Read(L1Brug, ctl_fpath) RESULT (status)

! Modules used:

!      USE Share_MOD_Log
!      USE Share_MOD_P_Constants_Variables, ONLY: LOGMSG, llvl, ctl_msg, &
!                                                 disp_msg, disp_line

! input/output parameters:
      CHARACTER(LEN=128),   INTENT(INOUT) :: ctl_fpath
      TYPE(L1B_rug),        INTENT(INOUT) :: L1Brug


! Local parameters For L1B Reading:
      CHARACTER*200                 :: l1b_rug_file_path
      CHARACTER*200                 :: ert_uv1_dgrp_path
      CHARACTER*200                 :: ert_uv1_ggrp_path
      CHARACTER*200                 :: ert_uv1_sgrp_path
      CHARACTER*200                 :: ert_uv2_dgrp_path
      CHARACTER*200                 :: ert_uv2_ggrp_path
      CHARACTER*200                 :: ert_uv2_sgrp_path

      CHARACTER*200                 :: grp_path
      CHARACTER*200                 :: data_name

      CHARACTER*200                 :: ExposureTime_DataName
      CHARACTER*200                 :: ExposureType_DataName
      CHARACTER*200                 :: PixelQFlags_DataName
      CHARACTER*200                 :: MstQFlags_DataName 
      CHARACTER*200                 :: NumsPixelCol_DataName
      CHARACTER*200                 :: RadExponent_DataName
      CHARACTER*200                 :: RadMantissa_DataName
      CHARACTER*200                 :: RadPcMantissa_DataName
      CHARACTER*200                 :: ReadoutTime_DataName
      CHARACTER*200                 :: WlenCoeff_DataName
      CHARACTER*200                 :: WlenCoeffPc_DataName
      CHARACTER*200                 :: WlenRefCol_DataName
      CHARACTER*200                 :: GdPixelQFlags_DataName
      CHARACTER*200                 :: Latitude_DataName
      CHARACTER*200                 :: Longitude_DataName
      CHARACTER*200                 :: SecondsInDay_DataName
      CHARACTER*200                 :: SolAziAng_DataName
      CHARACTER*200                 :: SolZenAng_DataName
      CHARACTER*200                 :: ScraftAlt_DataName
      CHARACTER*200                 :: ScraftLat_DataName
      CHARACTER*200                 :: ScraftLon_DataName
      CHARACTER*200                 :: TerrainHgt_DataName
      CHARACTER*200                 :: Time_DataName
      CHARACTER*200                 :: ViewAziAng_DataName
      CHARACTER*200                 :: ViewZenAng_DataName
      CHARACTER*200                 :: XTrackQFlags_DataName
      CHARACTER*200                 :: EarthSunDist_DataName

      INTEGER                       :: file_path_sz    
      INTEGER                       :: grp_path_sz 
      INTEGER                       :: data_name_sz

      INTEGER                       :: pstart(3)
      INTEGER                       :: pedge(3)

      INTEGER                       :: ierr        ! Error code 
      INTEGER                       :: hdferr      ! HDF5 Read Error code 
      INTEGER*8                     :: iAddr       ! memory address#

      INTEGER (KIND = 4)            :: status

! Local parameters For Logging:

      INTEGER       :: i, j, k

#ifdef INTEL
      INTEGER(1), PARAMETER :: X_FF = INT(Z'FF',1)
#endif

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
      Namelist /L1B_RUG_DATA_List/ExposureTime_DataName, &
                                  ExposureType_DataName, &
                                  PixelQFlags_DataName,  &
                                  MstQFlags_DataName,    &
                                  NumsPixelCol_DataName, &
                                  RadExponent_DataName,  & 
                                  RadMantissa_DataName,  &
                                  RadPcMantissa_DataName,&
                                  ReadoutTime_DataName,  &
                                  WlenCoeff_DataName,    &    
                                  WlenCoeffPc_DataName,  &  
                                  WlenRefCol_DataName,   &   
                                  GdPixelQFlags_DataName,&
                                  Latitude_DataName,     &     
                                  Longitude_DataName,    &    
                                  SecondsInDay_DataName, & 
                                  SolAziAng_DataName,    &    
                                  SolZenAng_DataName,    &    
                                  ScraftAlt_DataName,    &
                                  ScraftLat_DataName,    &    
                                  ScraftLon_DataName,    &    
                                  TerrainHgt_DataName,   &
                                  Time_DataName,         &
                                  ViewAziAng_DataName,   &   
                                  ViewZenAng_DataName,   &   
                                  XTrackQFlags_DataName, & 
                                  EarthSunDist_DataName
  
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
      CALL GEMS_Share_MOD_log(llvl, '--- Start L1BRUG READ FUNCTION  ---', 'START MSG')


!--------------------------------
!---     Read Control file    ---
!--------------------------------
      OPEN(10, file=trim(ctl_fpath), status='old')
      READ(10, L1B_RUG_File_List)
      READ(10, L1B_RUG_DATA_List)
      READ(10, L1B_File_List_Value_Size)
      CLOSE(10)

!--------------------------------
      !write(disp_line,'("",i6)') __LINE__
      !disp_msg = "program name=["//__FILE__//"], current line#=["//trim(disp_line)//"]"

!
!--------------------------------
!
      pstart(1) = -1                                            ! HDF5 데이터셑의 모든 데이터값을 읽는다는 의미

      !--- Setting to Global Variable of L1B File Path
      gds_l1b_fpath%rug = l1b_rug_file_path

      !-------------------------------------------------------------------------------------------------------------------------------
      !----------      L1BRUG     ----------!
      !--------------------------------------
      ! <Earth UV-1 Swath:Data Fields>
      !--------------------------------------
      disp_msg = "L1BRUG READ FILE PATH="//l1b_rug_file_path
      CALL GEMS_Share_MOD_log(llvl, disp_msg,  LOGMSG)
      status = lprint(disp_msg)

      CALL GEMS_Share_MOD_log(llvl, "Earth UV-1 Swath:Data Fields",  LOGMSG)

      !-----------------  ExposureTime  ---------------
      IF (ASSOCIATED(L1Brug%euv1%brug_dfld%ExpTime)) DEALLOCATE(L1Brug%euv1%brug_dfld%ExpTime)
      data_name = ExposureTime_DataName
      ALLOCATE (L1Brug%euv1%brug_dfld%ExpTime(nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv1_dgrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv1%brug_dfld%ExpTime, pstart, pedge)
      llvl = 9   
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv1%brug_dfld%ExpTime, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  ExposureType  ---------------
      IF (ASSOCIATED(L1Brug%euv1%brug_dfld%ExpType)) DEALLOCATE(L1Brug%euv1%brug_dfld%ExpType)
      data_name = ExposureType_DataName 
      ALLOCATE (L1Brug%euv1%brug_dfld%ExpType(nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv1_dgrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv1%brug_dfld%ExpType, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv1%brug_dfld%ExpType, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  MeasurementQualityFlags  ---------------
      IF (ASSOCIATED(L1Brug%euv1%brug_dfld%MQFlag)) DEALLOCATE(L1Brug%euv1%brug_dfld%MQFlag)
      data_name = MstQFlags_DataName
      grp_path  = ert_uv1_dgrp_path
      ALLOCATE (L1Brug%euv1%brug_dfld%MQFlag(nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv1_dgrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv1%brug_dfld%MQFlag, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv1%brug_dfld%MQFlag, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  NumberSmallPixelColumns  ---------------
      IF (ASSOCIATED(L1Brug%euv1%brug_dfld%NumSmPixCol)) DEALLOCATE(L1Brug%euv1%brug_dfld%NumSmPixCol)
      data_name = NumsPixelCol_DataName
      ALLOCATE (L1Brug%euv1%brug_dfld%NumSmPixCol(nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv1_dgrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv1%brug_dfld%NumSmPixCol, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv1%brug_dfld%NumSmPixCol,nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  PixelQualityFlags  ---------------
      IF (ASSOCIATED(L1Brug%euv1%brug_dfld%PQFlag)) DEALLOCATE(L1Brug%euv1%brug_dfld%PQFlag)
      data_name = PixelQFlags_DataName
      ALLOCATE (L1Brug%euv1%brug_dfld%PQFlag(nWavel_L1B,nxtrack_L1B,nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv1_dgrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv1%brug_dfld%PQFlag, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv1%brug_dfld%PQFlag, nWavel_L1B,nxtrack_L1B,nLine_L1B, 1, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  RadianceExponent ---------------
      IF (ASSOCIATED(L1Brug%euv1%brug_dfld%RadExponent)) DEALLOCATE(L1Brug%euv1%brug_dfld%RadExponent)
      data_name = RadExponent_DataName
      ALLOCATE (L1Brug%euv1%brug_dfld%RadExponent(nWavel_L1B,nxtrack_L1B,nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv1_dgrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv1%brug_dfld%RadExponent, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv1%brug_dfld%RadExponent, nWavel_L1B,nxtrack_L1B,nLine_L1B, 1, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  RadianceMantissa  ---------------
      IF (ASSOCIATED(L1Brug%euv1%brug_dfld%RadMantissa)) DEALLOCATE(L1Brug%euv1%brug_dfld%RadMantissa)
      data_name = RadMantissa_DataName
      ALLOCATE (L1Brug%euv1%brug_dfld%RadMantissa(nWavel_L1B,nxtrack_L1B,nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv1_dgrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv1%brug_dfld%RadMantissa, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv1%brug_dfld%RadMantissa, nWavel_L1B,nxtrack_L1B,nLine_L1B, 1, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  RadiancePrecisionMantissa  ---------------
      IF (ASSOCIATED(L1Brug%euv1%brug_dfld%RadPrecision)) DEALLOCATE(L1Brug%euv1%brug_dfld%RadPrecision)
      data_name = RadPcMantissa_DataName
      ALLOCATE (L1Brug%euv1%brug_dfld%RadPrecision(nWavel_L1B,nxtrack_L1B,nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv1_dgrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv1%brug_dfld%RadPrecision, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv1%brug_dfld%RadPrecision, nWavel_L1B,nxtrack_L1B,nLine_L1B, 1, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  ReadoutTime  ---------------
      IF (ASSOCIATED(L1Brug%euv1%brug_dfld%ReadTime)) DEALLOCATE(L1Brug%euv1%brug_dfld%ReadTime)
      data_name = ReadoutTime_DataName
      ALLOCATE (L1Brug%euv1%brug_dfld%ReadTime(nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv1_dgrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv1%brug_dfld%ReadTime, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv1%brug_dfld%ReadTime, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  WavelengthCoefficient  ---------------
      IF (ASSOCIATED(L1Brug%euv1%brug_dfld%WlenCoef)) DEALLOCATE(L1Brug%euv1%brug_dfld%WlenCoef)
      data_name = WlenCoeff_DataName
      ALLOCATE (L1Brug%euv1%brug_dfld%WlenCoef(nWavelCoef_L1B,nxtrack_L1B,nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv1_dgrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv1%brug_dfld%WlenCoef, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv1%brug_dfld%WlenCoef, nWavelCoef_L1B,nxtrack_L1B,nLine_L1B, 1, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !---------------  WavelengthCoefficientPrecision  -----------------
      IF (ASSOCIATED(L1Brug%euv1%brug_dfld%WlenPrec)) DEALLOCATE(L1Brug%euv1%brug_dfld%WlenPrec)
      data_name = WlenCoeffPc_DataName
      ALLOCATE (L1Brug%euv1%brug_dfld%WlenPrec(nWavelCoef_L1B,nxtrack_L1B,nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv1_dgrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv1%brug_dfld%WlenPrec, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv1%brug_dfld%WlenPrec, nWavelCoef_L1B,nxtrack_L1B,nLine_L1B, 1, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  WavelengthReferenceColumn  ---------------
      IF (ASSOCIATED(L1Brug%euv1%brug_dfld%WlenRefCol)) DEALLOCATE(L1Brug%euv1%brug_dfld%WlenRefCol)
      data_name = WlenRefCol_DataName
      ALLOCATE (L1Brug%euv1%brug_dfld%WlenRefCol(nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv1_dgrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv1%brug_dfld%WlenRefCol, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv1%brug_dfld%WlenRefCol, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !--------------------------------------
      ! <Earth UV-1 Swath : Geolocation Fields>
      !--------------------------------------
      CALL GEMS_Share_MOD_log(llvl, "Earth UV-1 Swath:Geolocation Fields",  LOGMSG)

      !-----------------  GroundPixelQualityFlags  ---------------
      IF (ASSOCIATED(L1Brug%euv1%brug_gfld%GPQFlag)) DEALLOCATE(L1Brug%euv1%brug_gfld%GPQFlag)
      data_name = GdPixelQFlags_DataName
      ALLOCATE (L1Brug%euv1%brug_gfld%GPQFlag(nxtrack_L1B, nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv1_ggrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv1%brug_gfld%GPQFlag, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv1%brug_gfld%GPQFlag, nxtrack_L1B, nLine_L1B, 0, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  Latitude  ---------------
      IF (ASSOCIATED(L1Brug%euv1%brug_gfld%LAT)) DEALLOCATE(L1Brug%euv1%brug_gfld%LAT)
      data_name = Latitude_DataName
      ALLOCATE (L1Brug%euv1%brug_gfld%LAT(nxtrack_L1B, nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv1_ggrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv1%brug_gfld%LAT, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv1%brug_gfld%LAT, nxtrack_L1B, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  Longitude  ---------------
      IF (ASSOCIATED(L1Brug%euv1%brug_gfld%LON)) DEALLOCATE(L1Brug%euv1%brug_gfld%LON)
      data_name = Longitude_DataName
      ALLOCATE (L1Brug%euv1%brug_gfld%LON(nxtrack_L1B, nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv1_ggrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv1%brug_gfld%LON, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv1%brug_gfld%LON, nxtrack_L1B, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF


      !-----------------  SecondsInDay  ---------------
      IF (ASSOCIATED(L1Brug%euv1%brug_gfld%SecInDay)) DEALLOCATE(L1Brug%euv1%brug_gfld%SecInDay)
      data_name = SecondsInDay_DataName
      ALLOCATE (L1Brug%euv1%brug_gfld%SecInDay(nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv1_ggrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv1%brug_gfld%SecInDay, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv1%brug_gfld%SecInDay, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SolarAzimuthAngle  ---------------
      IF (ASSOCIATED(L1Brug%euv1%brug_gfld%SolAziAng)) DEALLOCATE(L1Brug%euv1%brug_gfld%SolAziAng)
      data_name = SolAziAng_DataName
      ALLOCATE (L1Brug%euv1%brug_gfld%SolAziAng(nxtrack_L1B, nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv1_ggrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv1%brug_gfld%SolAziAng, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv1%brug_gfld%SolAziAng, nxtrack_L1B, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF
     
     
      !-----------------  SolarZenithAngle  ---------------
      IF (ASSOCIATED(L1Brug%euv1%brug_gfld%SolZenAng)) DEALLOCATE(L1Brug%euv1%brug_gfld%SolZenAng)
      data_name = SolZenAng_DataName
      ALLOCATE (L1Brug%euv1%brug_gfld%SolZenAng(nxtrack_L1B, nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv1_ggrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv1%brug_gfld%SolZenAng, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv1%brug_gfld%SolZenAng, nxtrack_L1B, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF
    
      !-----------------  SpacecraftAltitude  ---------------
      IF (ASSOCIATED(L1Brug%euv1%brug_gfld%ScAlt)) DEALLOCATE(L1Brug%euv1%brug_gfld%ScAlt)
      data_name = ScraftAlt_DataName
      ALLOCATE (L1Brug%euv1%brug_gfld%ScAlt(nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv1_ggrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv1%brug_gfld%ScAlt, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv1%brug_gfld%ScAlt, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SpacecraftLatitude  ---------------
      IF (ASSOCIATED(L1Brug%euv1%brug_gfld%ScLat)) DEALLOCATE(L1Brug%euv1%brug_gfld%ScLat)
      data_name = ScraftLat_DataName
      ALLOCATE (L1Brug%euv1%brug_gfld%ScLat(nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv1_ggrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv1%brug_gfld%ScLat, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv1%brug_gfld%ScLat, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  SpacecraftLongitude  ---------------
      IF (ASSOCIATED(L1Brug%euv1%brug_gfld%ScLon)) DEALLOCATE(L1Brug%euv1%brug_gfld%ScLon)
      data_name = ScraftLon_DataName
      ALLOCATE (L1Brug%euv1%brug_gfld%ScLon(nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv1_ggrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv1%brug_gfld%ScLon, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv1%brug_gfld%ScLon, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  TerrainHeight  ---------------
      IF (ASSOCIATED(L1Brug%euv1%brug_gfld%TerrHgt)) DEALLOCATE(L1Brug%euv1%brug_gfld%TerrHgt)
      data_name = TerrainHgt_DataName
      ALLOCATE (L1Brug%euv1%brug_gfld%TerrHgt(nxtrack_L1B, nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv1_ggrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv1%brug_gfld%TerrHgt, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv1%brug_gfld%TerrHgt, nxtrack_L1B, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  Time  ---------------
      IF (ASSOCIATED(L1Brug%euv1%brug_gfld%Time)) DEALLOCATE(L1Brug%euv1%brug_gfld%Time)
      data_name = Time_DataName
      ALLOCATE (L1Brug%euv1%brug_gfld%Time(nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv1_ggrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv1%brug_gfld%Time, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv1%brug_gfld%Time, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  ViewingAzimuthAngle  ---------------
      IF (ASSOCIATED(L1Brug%euv1%brug_gfld%ViewAziAng)) DEALLOCATE(L1Brug%euv1%brug_gfld%ViewAziAng)
      data_name = ViewAziAng_DataName
      ALLOCATE (L1Brug%euv1%brug_gfld%ViewAziAng(nxtrack_L1B, nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv1_ggrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv1%brug_gfld%ViewAziAng, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv1%brug_gfld%ViewAziAng, nxtrack_L1B, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  ViewingZenithAngle  ---------------
      IF (ASSOCIATED(L1Brug%euv1%brug_gfld%ViewZenAng)) DEALLOCATE(L1Brug%euv1%brug_gfld%ViewZenAng)
      data_name = ViewZenAng_DataName
      ALLOCATE (L1Brug%euv1%brug_gfld%ViewZenAng(nxtrack_L1B, nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv1_ggrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv1%brug_gfld%ViewZenAng, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv1%brug_gfld%ViewZenAng, nxtrack_L1B, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

!      !-----------------  XTrackQualityFlags   ---------------
      IF (ASSOCIATED(L1Brug%euv1%brug_gfld%XTQFlag)) DEALLOCATE(L1Brug%euv1%brug_gfld%XTQFlag)
      data_name = XTrackQFlags_DataName
      ALLOCATE (L1Brug%euv1%brug_gfld%XTQFlag(nxtrack_L1B, nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv1_ggrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv1%brug_gfld%XTQFlag, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv1%brug_gfld%XTQFlag, nxtrack_L1B, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF


      !--------------------------------------
      ! <Earth UV-1 Swath: Swath Attributes>
      !--------------------------------------
      CALL GEMS_Share_MOD_log(llvl, "Earth UV-1 Swath:Swath Attributes",  LOGMSG)

      !-----------------  EarthSunDistance  ---------------
      data_name = EarthSunDist_DataName
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv1_sgrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv1%brug_swatt%ES_Dist, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv1%brug_swatt%ES_Dist, data_name)
      ELSE
            GOTO 2999
      ENDIF


     !-------------------------------------------------------------------------------------------------------------------------------
     !----------L1BRUG----------!
     ! <Earth UV-2 Swath:Data Fileds>
      CALL GEMS_Share_MOD_log(llvl, "Earth UV-2 Swath:Data Fileds",  LOGMSG)

      !-----------------  ExposureTime  ---------------
      IF (ASSOCIATED(L1Brug%euv2%brug_dfld%ExpTime)) DEALLOCATE(L1Brug%euv2%brug_dfld%ExpTime)
      data_name = ExposureTime_DataName
      ALLOCATE (L1Brug%euv2%brug_dfld%ExpTime(nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv2_dgrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv2%brug_dfld%ExpTime, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv2%brug_dfld%ExpTime, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  ExposureType  ---------------
      IF (ASSOCIATED(L1Brug%euv2%brug_dfld%ExpType)) DEALLOCATE(L1Brug%euv2%brug_dfld%ExpType)
      data_name = ExposureType_DataName
      ALLOCATE (L1Brug%euv2%brug_dfld%ExpType(nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv2_dgrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv2%brug_dfld%ExpType, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv2%brug_dfld%ExpType, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  MeasurementQualityFlags  ---------------
      IF (ASSOCIATED(L1Brug%euv2%brug_dfld%MQFlag)) DEALLOCATE(L1Brug%euv2%brug_dfld%MQFlag)
      data_name = MstQFlags_DataName
      grp_path  = ert_uv2_dgrp_path
      ALLOCATE (L1Brug%euv2%brug_dfld%MQFlag(nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv2_dgrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv2%brug_dfld%MQFlag, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv2%brug_dfld%MQFlag, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  NumberSmallPixelColumns  ---------------
      IF (ASSOCIATED(L1Brug%euv2%brug_dfld%NumSmPixCol)) DEALLOCATE(L1Brug%euv2%brug_dfld%NumSmPixCol)
      data_name = NumsPixelCol_DataName
      ALLOCATE (L1Brug%euv2%brug_dfld%NumSmPixCol(nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv2_dgrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv2%brug_dfld%NumSmPixCol, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv2%brug_dfld%NumSmPixCol, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  PixelQualityFlags  ---------------
      IF (ASSOCIATED(L1Brug%euv2%brug_dfld%PQFlag)) DEALLOCATE(L1Brug%euv2%brug_dfld%PQFlag)
      data_name = PixelQFlags_DataName
      ALLOCATE (L1Brug%euv2%brug_dfld%PQFlag(nWavel2_L1B,nxtrack2_L1B,nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv2_dgrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv2%brug_dfld%PQFlag, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv2%brug_dfld%PQFlag, nWavel2_L1B, nxtrack2_L1B, nLine_L1B, 1, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  RadianceExponent ---------------
      IF (ASSOCIATED(L1Brug%euv2%brug_dfld%RadExponent)) DEALLOCATE(L1Brug%euv2%brug_dfld%RadExponent)
      data_name = RadExponent_DataName
      ALLOCATE (L1Brug%euv2%brug_dfld%RadExponent(nWavel2_L1B,nxtrack2_L1B,nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv2_dgrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv2%brug_dfld%RadExponent, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv2%brug_dfld%RadExponent, nWavel2_L1B, nxtrack2_L1B, nLine_L1B, 1, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  RadianceMantissa  ---------------
      IF (ASSOCIATED(L1Brug%euv2%brug_dfld%RadMantissa)) DEALLOCATE(L1Brug%euv2%brug_dfld%RadMantissa)
      data_name = RadMantissa_DataName
      ALLOCATE (L1Brug%euv2%brug_dfld%RadMantissa(nWavel2_L1B,nxtrack2_L1B,nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv2_dgrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv2%brug_dfld%RadMantissa, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv2%brug_dfld%RadMantissa, nWavel2_L1B, nxtrack2_L1B, nLine_L1B, 1, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  RadiancePrecisionMantissa  ---------------
      IF (ASSOCIATED(L1Brug%euv2%brug_dfld%RadPrecision)) DEALLOCATE(L1Brug%euv2%brug_dfld%RadPrecision)
      data_name = RadPcMantissa_DataName
      ALLOCATE (L1Brug%euv2%brug_dfld%RadPrecision(nWavel2_L1B,nxtrack2_L1B,nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv2_dgrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv2%brug_dfld%RadPrecision, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv2%brug_dfld%RadPrecision, nWavel2_L1B, nxtrack2_L1B, nLine_L1B, 1, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  ReadoutTime  ---------------
      IF (ASSOCIATED(L1Brug%euv2%brug_dfld%ReadTime)) DEALLOCATE(L1Brug%euv2%brug_dfld%ReadTime)
      data_name = ReadoutTime_DataName
      ALLOCATE (L1Brug%euv2%brug_dfld%ReadTime(nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv2_dgrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv2%brug_dfld%ReadTime, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv2%brug_dfld%ReadTime, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  WavelengthCoefficient  ---------------
      IF (ASSOCIATED(L1Brug%euv2%brug_dfld%WlenCoef)) DEALLOCATE(L1Brug%euv2%brug_dfld%WlenCoef)
      data_name = WlenCoeff_DataName
      ALLOCATE (L1Brug%euv2%brug_dfld%WlenCoef(nWavelCoef_L1B,nxtrack2_L1B,nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv2_dgrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv2%brug_dfld%WlenCoef, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv2%brug_dfld%WlenCoef, nWavelCoef_L1B,nxtrack2_L1B,nLine_L1B, 1, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  WavelengthReferenceColumn  ---------------
      IF (ASSOCIATED(L1Brug%euv2%brug_dfld%WlenRefCol)) DEALLOCATE(L1Brug%euv2%brug_dfld%WlenRefCol)
      data_name = WlenRefCol_DataName
      ALLOCATE (L1Brug%euv2%brug_dfld%WlenRefCol(nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv2_dgrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv2%brug_dfld%WlenRefCol, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv2%brug_dfld%WlenRefCol, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !--------------------------------------
      ! <Earth UV-2 Swath: Geolocation Fields>
      !--------------------------------------
      CALL GEMS_Share_MOD_log(llvl, "Earth UV-2 Swath:Geolocation Fileds",  LOGMSG)

      !-----------------  GroundPixelQualityFlags  ---------------
      IF (ASSOCIATED(L1Brug%euv2%brug_gfld%GPQFlag)) DEALLOCATE(L1Brug%euv2%brug_gfld%GPQFlag)
      data_name = GdPixelQFlags_DataName
      ALLOCATE (L1Brug%euv2%brug_gfld%GPQFlag(nxtrack2_L1B, nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv2_ggrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv2%brug_gfld%GPQFlag, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv2%brug_gfld%GPQFlag, nxtrack2_L1B, nLine_L1B, 0, data_name)
      ELSE
           GOTO 2999
      ENDIF


      !-----------------  Latitude  ---------------
      IF (ASSOCIATED(L1Brug%euv2%brug_gfld%LAT)) DEALLOCATE(L1Brug%euv2%brug_gfld%LAT)
      data_name = Latitude_DataName
      ALLOCATE (L1Brug%euv2%brug_gfld%LAT(nxtrack2_L1B, nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv2_ggrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv2%brug_gfld%LAT, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv2%brug_gfld%LAT, nxtrack2_L1B, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  Longitude  ---------------
      IF (ASSOCIATED(L1Brug%euv2%brug_gfld%LON)) DEALLOCATE(L1Brug%euv2%brug_gfld%LON)
      data_name = Longitude_DataName
      ALLOCATE (L1Brug%euv2%brug_gfld%LON(nxtrack2_L1B, nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv2_ggrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv2%brug_gfld%LON, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv2%brug_gfld%LON, nxtrack2_L1B, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF


      !-----------------  SolarAzimuthAngle  ---------------
      IF (ASSOCIATED(L1Brug%euv2%brug_gfld%SolAziAng)) DEALLOCATE(L1Brug%euv2%brug_gfld%SolAziAng)
      data_name = SolAziAng_DataName
      ALLOCATE (L1Brug%euv2%brug_gfld%SolAziAng(nxtrack2_L1B, nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv2_ggrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv2%brug_gfld%SolAziAng, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv2%brug_gfld%SolAziAng, nxtrack2_L1B, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF
     
     
      !-----------------  SolarZenithAngle  ---------------
      IF (ASSOCIATED(L1Brug%euv2%brug_gfld%SolZenAng)) DEALLOCATE(L1Brug%euv2%brug_gfld%SolZenAng)
      data_name = SolZenAng_DataName
      ALLOCATE (L1Brug%euv2%brug_gfld%SolZenAng(nxtrack2_L1B, nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv2_ggrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv2%brug_gfld%SolZenAng, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv2%brug_gfld%SolZenAng, nxtrack2_L1B, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF
  
      !-----------------  SpacecraftAltitude  ---------------
      IF (ASSOCIATED(L1Brug%euv2%brug_gfld%ScAlt)) DEALLOCATE(L1Brug%euv2%brug_gfld%ScAlt)
      data_name = ScraftAlt_DataName
      ALLOCATE (L1Brug%euv2%brug_gfld%ScAlt(nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv2_ggrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv2%brug_gfld%ScAlt, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv2%brug_gfld%ScAlt, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF
 
      !-----------------  TerrainHeight  ---------------
      IF (ASSOCIATED(L1Brug%euv2%brug_gfld%TerrHgt)) DEALLOCATE(L1Brug%euv2%brug_gfld%TerrHgt)
      data_name = TerrainHgt_DataName
      ALLOCATE (L1Brug%euv2%brug_gfld%TerrHgt(nxtrack2_L1B, nLine_L1B), STAT =ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz,ert_uv2_ggrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv2%brug_gfld%TerrHgt, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv2%brug_gfld%TerrHgt,nxtrack2_L1B, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  Time  ---------------
      IF (ASSOCIATED(L1Brug%euv2%brug_gfld%Time)) DEALLOCATE(L1Brug%euv2%brug_gfld%Time)
      data_name = Time_DataName
      ALLOCATE (L1Brug%euv2%brug_gfld%Time(nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv2_ggrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv2%brug_gfld%Time, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv2%brug_gfld%Time, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF
 
      !-----------------  ViewingAzimuthAngle  ---------------
      IF (ASSOCIATED(L1Brug%euv2%brug_gfld%ViewAziAng)) DEALLOCATE(L1Brug%euv2%brug_gfld%ViewAziAng)
      data_name = ViewAziAng_DataName
      ALLOCATE (L1Brug%euv2%brug_gfld%ViewAziAng(nxtrack2_L1B, nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv2_ggrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv2%brug_gfld%ViewAziAng, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv2%brug_gfld%ViewAziAng, nxtrack2_L1B, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  ViewingZenithAngle  ---------------
      IF (ASSOCIATED(L1Brug%euv2%brug_gfld%ViewZenAng)) DEALLOCATE(L1Brug%euv2%brug_gfld%ViewZenAng)
      data_name = ViewZenAng_DataName
      ALLOCATE (L1Brug%euv2%brug_gfld%ViewZenAng(nxtrack2_L1B, nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv2_ggrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv2%brug_gfld%ViewZenAng, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv2%brug_gfld%ViewZenAng, nxtrack2_L1B, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  XTrackQualityFlags   ---------------
      IF (ASSOCIATED(L1Brug%euv2%brug_gfld%XTQFlag)) DEALLOCATE(L1Brug%euv2%brug_gfld%XTQFlag)
      data_name = XTrackQFlags_DataName
      ALLOCATE (L1Brug%euv2%brug_gfld%XTQFlag(nxtrack2_L1B, nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv2_ggrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv2%brug_gfld%XTQFlag, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv2%brug_gfld%XTQFlag, nxtrack2_L1B, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !--------------------------------------
      ! <Earth UV-2 Swath: Swath Attributes>
      !--------------------------------------
      CALL GEMS_Share_MOD_log(llvl, "Earth UV-2 Swath:Swath Attributes",LOGMSG)

      !-----------------  EarthSunDistance  ---------------
      data_name = EarthSunDist_DataName
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rug_file_path, file_path_sz, ert_uv2_sgrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brug%euv2%brug_swatt%ES_Dist, pstart, pedge)
      IF (hdferr == 0) THEN 
            CALL GEMS_Share_MOD_log(llvl, L1Brug%euv2%brug_swatt%ES_Dist, data_name)
      ELSE
            GOTO 2999
      ENDIF
      llvl = 0   

      write(disp_msg,'(" OMI HDF5 File READ OK, hdferr=",i6)') hdferr
      CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)

!--------------------------------
!---- END FUNCTION -------------
!--------------------------------
      CALL GEMS_Share_MOD_log(llvl, '---   END L1BRUG READ FUNCTION ---', '  END MSG')
      CALL GEMS_Share_MOD_log(llvl, '.', LOGMSG)
      status = 0

      RETURN 

999   CONTINUE
      llvl = 0   
      write(disp_msg,'( A20, "Memory Allocation Failed, ierr=",i6)') ctl_msg, ierr
      CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)
      status = ierr

      RETURN 

2999  CONTINUE
      llvl = 0   
      write(disp_msg,'(" !!! OMI HDF5 File READ ERROR !!! ", "[", A30, "], hdferr=",i6)') data_name, hdferr
      CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)
      status = hdferr

      RETURN 

      END FUNCTION GEMS_Share_MOD_L1Brug_Read


!-------------------------------------------------------------------------------
!+Description: 
!      FUNCTION GEMS_Share_MOD_L1Brvg_Read .
!
!
! Method:
!
!
! Input files:
!       L1Brvg    : the Level 1B RVG data structure
!       ctl_fpath : the control file path
!
! Output files:
!       status    : the result code of L1B RVG read
!
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2014.12.23  Fisrt Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------

      
! Declarations:

      FUNCTION GEMS_Share_MOD_L1Brvg_Read(L1Brvg, ctl_fpath) RESULT (status)


! Modules used:

!      USE Share_MOD_Log
!      USE Share_MOD_P_Constants_Variables, ONLY: LOGMSG, llvl, ctl_msg, &
!                                                 disp_msg, disp_line

! input/output parameters:
      CHARACTER(LEN=128),   INTENT(INOUT) :: ctl_fpath
      TYPE(L1B_rvg),        INTENT(INOUT) :: L1Brvg


! Local parameters For L1B Reading:
      CHARACTER*200                 :: l1b_rvg_file_path
      CHARACTER*200                 :: ert_vis_dgrp_path
      CHARACTER*200                 :: ert_vis_ggrp_path
      CHARACTER*200                 :: ert_vis_sgrp_path

      CHARACTER*200                 :: grp_path
      CHARACTER*200                 :: data_name

      CHARACTER*200                 :: PixelQFlags_DataName
      CHARACTER*200                 :: RadExponent_DataName
      CHARACTER*200                 :: RadMantissa_DataName
      CHARACTER*200                 :: WlenCoeff_DataName
      CHARACTER*200                 :: WlenRefCol_DataName
      CHARACTER*200                 :: Longitude_DataName
      CHARACTER*200                 :: Latitude_DataName
      CHARACTER*200                 :: SolAziAng_DataName
      CHARACTER*200                 :: SolZenAng_DataName
      CHARACTER*200                 :: TerrainHgt_DataName
      CHARACTER*200                 :: ViewZenAng_DataName
      CHARACTER*200                 :: ViewAziAng_DataName
      CHARACTER*200                 :: EarthSunDist_DataName

      INTEGER                       :: file_path_sz    
      INTEGER                       :: grp_path_sz 
      INTEGER                       :: data_name_sz

      INTEGER                       :: pstart(3)
      INTEGER                       :: pedge(3)

      INTEGER                       :: ierr        ! Error code 
      INTEGER                       :: hdferr      ! HDF5 Read Error code 
      INTEGER*8                     :: iAddr       ! memory address#

      INTEGER (KIND = 4)            :: status

! Local parameters For Logging:

      INTEGER       :: i, j, k

#ifdef INTEL
      INTEGER(1), PARAMETER :: X_FF = INT(Z'FF',1)
#endif

      Namelist /L1B_RVG_File_List/l1b_rvg_file_path, &
                                  ert_vis_dgrp_path, &
                                  ert_vis_ggrp_path, &
                                  ert_vis_sgrp_path
      Namelist /L1B_File_List_Value_Size/file_path_sz, &
                                         grp_path_sz,  &
                                         data_name_sz
      Namelist /L1B_RVG_DATA_List/PixelQFlags_DataName, &
                                  RadExponent_DataName, &
                                  RadMantissa_DataName, &
                                  WlenCoeff_DataName,   &
                                  WlenRefCol_DataName,  &
                                  Latitude_DataName,    &
                                  Longitude_DataName,   &
                                  SolAziAng_DataName,   &
                                  SolZenAng_DataName,   &
                                  TerrainHgt_DataName,  &
                                  ViewAziAng_DataName,  &
                                  ViewZenAng_DataName,  &
                                  EarthSunDist_DataName  

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
      CALL GEMS_Share_MOD_log(llvl, '--- Start L1BRVG READ FUNCTION ---', 'START MSG')


!--------------------------------
!---     Read Control file    ---
!--------------------------------
      OPEN(10, file=trim(ctl_fpath), status='old')
      READ(10, L1B_RVG_File_List)
      READ(10, L1B_RVG_DATA_List)
      READ(10, L1B_File_List_Value_Size)
      CLOSE(10)

!--------------------------------
      !write(disp_line,'("",i6)') __LINE__
      !disp_msg = "program name=["//__FILE__//"], current line#=["//trim(disp_line)//"]"

!
!--------------------------------
!
      pstart(1) = -1                                            ! HDF5 데이터셑의 모든 데이터값을 읽는다는 의미

      !--- Setting to Global Variable of L1B File Path
      gds_l1b_fpath%rvg = l1b_rvg_file_path

      !-------------------------------------------------------------------------------------------------------------------------------
      !----------      L1BRVG     ----------!
      !--------------------------------------
      ! <Earth VIS Swath:Data Fields>
      !--------------------------------------
      disp_msg = "L1BRVG READ FILE PATH="//l1b_rvg_file_path
      CALL GEMS_Share_MOD_log(llvl, disp_msg,  LOGMSG)
      status = lprint(disp_msg)

      CALL GEMS_Share_MOD_log(llvl, "Earth VIS Swath:Data Fields",  LOGMSG)

      !-----------------  PixelQualityFlags ---------------
      IF (ASSOCIATED(L1Brvg%brvg_dfld%PQFlag)) DEALLOCATE(L1Brvg%brvg_dfld%PQFlag)
      data_name = PixelQFlags_DataName
      ALLOCATE(L1Brvg%brvg_dfld%PQFlag(nWavel3_L1B,nxtrack2_L1B,nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rvg_file_path, file_path_sz,ert_vis_dgrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brvg%brvg_dfld%PQFlag, pstart, pedge)
      llvl = 9   
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brvg%brvg_dfld%PQFlag, nWavel3_L1B,nxtrack2_L1B,nLine_L1B, 1, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  RadianceExponent ---------------
      IF (ASSOCIATED(L1Brvg%brvg_dfld%RadExponent)) DEALLOCATE(L1Brvg%brvg_dfld%RadExponent)
      data_name = RadExponent_DataName
      ALLOCATE (L1Brvg%brvg_dfld%RadExponent(nWavel3_L1B,nxtrack2_L1B,nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rvg_file_path, file_path_sz, ert_vis_dgrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brvg%brvg_dfld%RadExponent, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brvg%brvg_dfld%RadExponent, nWavel3_L1B,nxtrack2_L1B,nLine_L1B, 1, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  RadianceMantissa  ---------------
      IF (ASSOCIATED(L1Brvg%brvg_dfld%RadMantissa)) DEALLOCATE(L1Brvg%brvg_dfld%RadMantissa)
      data_name = RadMantissa_DataName
      ALLOCATE (L1Brvg%brvg_dfld%RadMantissa(nWavel3_L1B,nxtrack2_L1B,nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rvg_file_path, file_path_sz, ert_vis_dgrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brvg%brvg_dfld%RadMantissa, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brvg%brvg_dfld%RadMantissa, nWavel3_L1B,nxtrack2_L1B,nLine_L1B, 1, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  WavelengthCoefficient  ---------------
      IF (ASSOCIATED(L1Brvg%brvg_dfld%WlenCoef)) DEALLOCATE(L1Brvg%brvg_dfld%WlenCoef)
      data_name = WlenCoeff_DataName
      ALLOCATE (L1Brvg%brvg_dfld%WlenCoef(nWavelCoef_L1B,nxtrack2_L1B,nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rvg_file_path, file_path_sz, ert_vis_dgrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brvg%brvg_dfld%WlenCoef, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brvg%brvg_dfld%WlenCoef, nWavelCoef_L1B,nxtrack2_L1B,nLine_L1B, 1, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  WavelengthReferenceColumn  ---------------
      IF (ASSOCIATED(L1Brvg%brvg_dfld%WlenRefCol)) DEALLOCATE(L1Brvg%brvg_dfld%WlenRefCol)
      data_name = WlenRefCol_DataName
      ALLOCATE (L1Brvg%brvg_dfld%WlenRefCol(nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rvg_file_path, file_path_sz, ert_vis_dgrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brvg%brvg_dfld%WlenRefCol, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brvg%brvg_dfld%WlenRefCol, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !--------------------------------------
      ! <Earth VIS Swath/Geolocation Fields>
      !--------------------------------------
      CALL GEMS_Share_MOD_log(llvl, "Earth VIS Swath:Geolocation Fields",  LOGMSG)

      print*,"TEST2"

      !-----------------  Latitude  ---------------
      IF (ASSOCIATED(L1Brvg%brvg_gfld%LAT)) DEALLOCATE(L1Brvg%brvg_gfld%LAT)
      data_name = Latitude_DataName
      ALLOCATE (L1Brvg%brvg_gfld%LAT(nxtrack2_L1B, nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rvg_file_path, file_path_sz, ert_vis_ggrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brvg%brvg_gfld%LAT, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brvg%brvg_gfld%LAT, nxtrack2_L1B, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  Longitude  ---------------
      IF (ASSOCIATED(L1Brvg%brvg_gfld%LON)) DEALLOCATE(L1Brvg%brvg_gfld%LON)
      data_name = Longitude_DataName
      ALLOCATE (L1Brvg%brvg_gfld%LON(nxtrack2_L1B, nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rvg_file_path, file_path_sz, ert_vis_ggrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brvg%brvg_gfld%LON, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brvg%brvg_gfld%LON, nxtrack2_L1B, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF


      !-----------------  SolarAzimuthAngle  ---------------
      IF (ASSOCIATED(L1Brvg%brvg_gfld%SolAziAng)) DEALLOCATE(L1Brvg%brvg_gfld%SolAziAng)
      data_name = SolAziAng_DataName
      ALLOCATE (L1Brvg%brvg_gfld%SolAziAng(nxtrack2_L1B, nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rvg_file_path, file_path_sz, ert_vis_ggrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brvg%brvg_gfld%SolAziAng, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brvg%brvg_gfld%SolAziAng, nxtrack2_L1B, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF
     
     
      !-----------------  SolarZenithAngle  ---------------
      IF (ASSOCIATED(L1Brvg%brvg_gfld%SolZenAng)) DEALLOCATE(L1Brvg%brvg_gfld%SolZenAng)
      data_name = SolZenAng_DataName
      ALLOCATE (L1Brvg%brvg_gfld%SolZenAng(nxtrack2_L1B, nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rvg_file_path, file_path_sz, ert_vis_ggrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brvg%brvg_gfld%SolZenAng, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brvg%brvg_gfld%SolZenAng, nxtrack2_L1B, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  TerrainHeight  ---------------
      IF (ASSOCIATED(L1Brvg%brvg_gfld%TerrHgt)) DEALLOCATE(L1Brvg%brvg_gfld%TerrHgt)
      data_name = TerrainHgt_DataName
      ALLOCATE (L1Brvg%brvg_gfld%TerrHgt(nxtrack2_L1B,nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rvg_file_path, file_path_sz, ert_vis_ggrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brvg%brvg_gfld%TerrHgt, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brvg%brvg_gfld%TerrHgt, nxtrack2_L1B,nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF
    
      !-----------------  ViewingAzimuthAngle  ---------------
      IF (ASSOCIATED(L1Brvg%brvg_gfld%ViewAziAng)) DEALLOCATE(L1Brvg%brvg_gfld%ViewAziAng)
      data_name = ViewAziAng_DataName
      ALLOCATE (L1Brvg%brvg_gfld%ViewAziAng(nxtrack2_L1B, nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rvg_file_path, file_path_sz, ert_vis_ggrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brvg%brvg_gfld%ViewAziAng, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brvg%brvg_gfld%ViewAziAng, nxtrack2_L1B, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  ViewingZenithAngle  ---------------
      IF (ASSOCIATED(L1Brvg%brvg_gfld%ViewZenAng)) DEALLOCATE(L1Brvg%brvg_gfld%ViewZenAng)
      data_name = ViewZenAng_DataName
      ALLOCATE (L1Brvg%brvg_gfld%ViewZenAng(nxtrack2_L1B, nLine_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rvg_file_path, file_path_sz, ert_vis_ggrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brvg%brvg_gfld%ViewZenAng, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brvg%brvg_gfld%ViewZenAng, nxtrack2_L1B, nLine_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      print*,"TEST3"
      !--------------------------------------
      ! <Earth VIS Swath: Swath Attributes>
      !--------------------------------------
      CALL GEMS_Share_MOD_log(llvl, "Earth VIS Swath:Swath Attributes",LOGMSG)

      !-----------------  EarthSunDistance  ---------------
      data_name = EarthSunDist_DataName
      hdferr = GEMS_Share_Hdf5ReadData(l1b_rvg_file_path, file_path_sz,ert_vis_sgrp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Brvg%brvg_swatt%ES_Dist, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Brvg%brvg_swatt%ES_Dist,data_name)
      ELSE
            GOTO 2999
      ENDIF
      llvl = 0   

      write(disp_msg,'(" OMI HDF5 File READ OK, hdferr=",i6)') hdferr
      CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)



!--------------------------------
!---- END FUNCTION -------------
!--------------------------------
      CALL GEMS_Share_MOD_log(llvl, '---   END L1BRVG READ FUNCTION ---', '  END MSG')
      CALL GEMS_Share_MOD_log(llvl, '.', LOGMSG)
      status = 0

      RETURN 

999   CONTINUE
      llvl = 0   
      write(disp_msg,'( A20, "Memory Allocation Failed, ierr=",i6)') ctl_msg, ierr
      CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)
      status = ierr

      RETURN 

2999  CONTINUE
      llvl = 0   
      write(disp_msg,'(" !!! OMI HDF5 File READ ERROR !!! ", "[", A30, "], hdferr=",i6)') data_name, hdferr
      CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)
      status = hdferr

      RETURN 

      END FUNCTION GEMS_Share_MOD_L1Brvg_Read

!-------------------------------------------------------------------------------
!+Description: 
!      FUNCTION GEMS_Share_MOD_L1Birr_Read .
!
!
! Method:
!
!
! Input files:
!       L1Birr    : the Level 1B IRR data structure
!       ctl_fpath : the control file path
!
! Output files:
!       status    : the result code of L1B IRR read
!
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2014.12.23  Fisrt Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------

      
! Declarations:

      FUNCTION GEMS_Share_MOD_L1Birr_Read(L1Birr, ctl_fpath) RESULT (status)


! Modules used:

!      USE Share_MOD_Log
!      USE Share_MOD_P_Constants_Variables, ONLY: LOGMSG, llvl, ctl_msg, &
!                                                 disp_msg, disp_line

! input/output parameters:
      CHARACTER(LEN=128),   INTENT(INOUT) :: ctl_fpath
      TYPE(L1B_irr),        INTENT(INOUT) :: L1Birr

! Local parameters For L1B Reading:
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

      CHARACTER*200                 :: grp_path
      CHARACTER*200                 :: data_name

      CHARACTER*200                 :: IrrExponent_DataName
      CHARACTER*200                 :: IrrMantissa_DataName
      CHARACTER*200                 :: IrrPcMantissa_DataName
      CHARACTER*200                 :: PixelQFlags_DataName
      CHARACTER*200                 :: WlenCoeff_DataName
      CHARACTER*200                 :: WlenRefCol_DataName
      CHARACTER*200                 :: EarthSunDist_DataName

      INTEGER                       :: file_path_sz    
      INTEGER                       :: grp_path_sz 
      INTEGER                       :: data_name_sz

      INTEGER                       :: pstart(3)
      INTEGER                       :: pedge(3)

      INTEGER                       :: ierr        ! Error code 
      INTEGER                       :: hdferr      ! HDF5 Read Error code 
      INTEGER*8                     :: iAddr       ! memory address#

      INTEGER (KIND = 4)            :: status

! Local parameters For Logging:

      INTEGER       :: i, j, k

#ifdef INTEL
      INTEGER(1), PARAMETER :: X_FF = INT(Z'FF',1)
#endif

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
      Namelist /L1B_IRR_DATA_List/IrrExponent_DataName,  &
                                  IrrMantissa_DataName,  &
                                  IrrPcMantissa_DataName,& 
                                  PixelQFlags_DataName,  &
                                  WlenCoeff_DataName,    &
                                  WlenRefCol_DataName,   &
                                  EarthSunDist_DataName

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
      CALL GEMS_Share_MOD_log(llvl, '--- Start L1BIRR READ FUNCTION  ---', 'START MSG')


!--------------------------------
!---     Read Control file    ---
!--------------------------------
!      namelist_file = 'l1breadmdl.ctl'
      OPEN(10, file=trim(ctl_fpath), status='old')
      READ(10, L1B_IRR_File_List)
      READ(10, L1B_IRR_DATA_List)
      READ(10, L1B_File_List_Value_Size)
      CLOSE(10)

!
!--------------------------------
!
      pstart(1) = -1                                            ! HDF5 데이터셑의 모든 데이터값을 읽는다는 의미

      !--- Setting to Global Variable of L1B File Path
      gds_l1b_fpath%irr = l1b_irr_file_path

      !-------------------------------------------------------------------------------------------------------------------------------
      !----------      L1BIRR     ----------!
      !--------------------------------------
      ! <Sun Volume UV-1 Swath:Data Fields>
      !--------------------------------------

      disp_msg = "L1BIRR READ FILE PATH="//l1b_irr_file_path
      CALL GEMS_Share_MOD_log(llvl, disp_msg,  LOGMSG)
      status = lprint(disp_msg)

      CALL GEMS_Share_MOD_log(llvl, "Sun Volume UV-1 Swath:Data Fields",  LOGMSG)

      !-----------------  IrradianceExponent ---------------
!      IF (ASSOCIATED(L1Birr%suv1%birr_dfld%IrrRadExponent)) DEALLOCATE(L1Birr%suv1%birr_dfld%IrrRadExponent)
      IF (ASSOCIATED(L1Birr%suv1%birr_dfld%IrrRadExponent)) DEALLOCATE(L1Birr%suv1%birr_dfld%IrrRadExponent)
      data_name = IrrExponent_DataName
      grp_path  = sun_uv1_dgrp_path
!      ALLOCATE (L1Birr%suv1%birr_dfld%IrrRadExponent(nWavel_L1B,nxtrack_L1B), STAT = ierr)
      ALLOCATE (L1Birr%suv1%birr_dfld%IrrRadExponent(nWavel_L1B,nxtrack_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_irr_file_path, file_path_sz, grp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Birr%suv1%birr_dfld%IrrRadExponent, pstart, pedge)
      llvl = 9   
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Birr%suv1%birr_dfld%IrrRadExponent, nWavel_L1B, nxtrack_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  IrradianceMantissa  ---------------
      IF (ASSOCIATED(L1Birr%suv1%birr_dfld%IrrRadMantissa)) DEALLOCATE(L1Birr%suv1%birr_dfld%IrrRadMantissa)
      data_name = IrrMantissa_DataName
      ALLOCATE (L1Birr%suv1%birr_dfld%IrrRadMantissa(nWavel_L1B,nxtrack_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_irr_file_path, file_path_sz, grp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Birr%suv1%birr_dfld%IrrRadMantissa, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Birr%suv1%birr_dfld%IrrRadMantissa, nWavel_L1B, nxtrack_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  WavelengthCoefficient  ---------------
      IF (ASSOCIATED(L1Birr%suv1%birr_dfld%WlenCoef)) DEALLOCATE(L1Birr%suv1%birr_dfld%WlenCoef)
      data_name = WlenCoeff_DataName
      ALLOCATE (L1Birr%suv1%birr_dfld%WlenCoef(nWavelCoef_L1B,nxtrack_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_irr_file_path, file_path_sz, grp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Birr%suv1%birr_dfld%WlenCoef, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Birr%suv1%birr_dfld%WlenCoef, nWavelCoef_L1B,nxtrack_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  WavelengthReferenceColumn  ---------------
      data_name = WlenRefCol_DataName
      hdferr = GEMS_Share_Hdf5ReadData(l1b_irr_file_path, file_path_sz, grp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Birr%suv1%birr_dfld%WlenRefCol, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Birr%suv1%birr_dfld%WlenRefCol, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !--------------------------------------
      ! <Sun Volume UV-1 Swath:Swath Attributes>
      !--------------------------------------
      CALL GEMS_Share_MOD_log(llvl, "Sun Volume UV-1 Swath:Swath Attributes",  LOGMSG)

      !-----------------  EarthSunDistance  ---------------
      data_name = EarthSunDist_DataName
      grp_path  = sun_uv1_sgrp_path
      hdferr = GEMS_Share_Hdf5ReadData(l1b_irr_file_path, file_path_sz, grp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Birr%suv1%birr_swatt%ES_Dist, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Birr%suv1%birr_swatt%ES_Dist, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-------------------------------------------------------------------------------------------------------------------------------
      !----------      L1BIRR     ----------!
      !--------------------------------------
      ! <Sun Volume UV-2 Swath:Data Fields>
      !--------------------------------------

      disp_msg = "L1BIRR READ FILE PATH="//l1b_irr_file_path
      CALL GEMS_Share_MOD_log(llvl, disp_msg,  LOGMSG)
      status = lprint(disp_msg)

      CALL GEMS_Share_MOD_log(llvl, "Sun Volume UV-2 Swath:Data Fields",  LOGMSG)

      !-----------------  IrradianceExponent ---------------
!      IF (ASSOCIATED(L1Birr%suv2%birr_dfld%IrrRadExponent)) DEALLOCATE(L1Birr%suv2%birr_dfld%IrrRadExponent)
      IF (ASSOCIATED(L1Birr%suv2%birr_dfld%IrrRadExponent)) DEALLOCATE(L1Birr%suv2%birr_dfld%IrrRadExponent)
      data_name = IrrExponent_DataName
      grp_path  = sun_uv2_dgrp_path
!      ALLOCATE (L1Birr%suv2%birr_dfld%IrrRadExponent(nWavel2_L1B,nxtrack2_L1B), STAT = ierr)
      ALLOCATE (L1Birr%suv2%birr_dfld%IrrRadExponent(nWavel2_L1B,nxtrack2_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_irr_file_path, file_path_sz, grp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Birr%suv2%birr_dfld%IrrRadExponent, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Birr%suv2%birr_dfld%IrrRadExponent, nWavel2_L1B, nxtrack2_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  IrradianceMantissa  ---------------
      IF (ASSOCIATED(L1Birr%suv2%birr_dfld%IrrRadMantissa)) DEALLOCATE(L1Birr%suv2%birr_dfld%IrrRadMantissa)
      data_name = IrrMantissa_DataName
      ALLOCATE (L1Birr%suv2%birr_dfld%IrrRadMantissa(nWavel2_L1B,nxtrack2_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_irr_file_path, file_path_sz, grp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Birr%suv2%birr_dfld%IrrRadMantissa, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Birr%suv2%birr_dfld%IrrRadMantissa, nWavel2_L1B, nxtrack2_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  IrradiancePrecisionMantissa  ---------------
      IF (ASSOCIATED(L1Birr%suv2%birr_dfld%IrrRadPrecision)) DEALLOCATE(L1Birr%suv2%birr_dfld%IrrRadPrecision)
      data_name = IrrPcMantissa_DataName
      ALLOCATE (L1Birr%suv2%birr_dfld%IrrRadPrecision(nWavel2_L1B,nxtrack2_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_irr_file_path, file_path_sz, grp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Birr%suv2%birr_dfld%IrrRadPrecision, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Birr%suv2%birr_dfld%IrrRadPrecision, nWavel2_L1B, nxtrack2_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  PixelQualityFlags  ---------------
      IF (ASSOCIATED(L1Birr%suv2%birr_dfld%PQFlag)) DEALLOCATE(L1Birr%suv2%birr_dfld%PQFlag)
      data_name = PixelQFlags_DataName
      ALLOCATE (L1Birr%suv2%birr_dfld%PQFlag(nWavel2_L1B,nxtrack2_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_irr_file_path, file_path_sz, grp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Birr%suv2%birr_dfld%PQFlag, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Birr%suv2%birr_dfld%PQFlag, nWavel2_L1B, nxtrack2_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF


      !-----------------  WavelengthCoefficient  ---------------
      IF (ASSOCIATED(L1Birr%suv2%birr_dfld%WlenCoef)) DEALLOCATE(L1Birr%suv2%birr_dfld%WlenCoef)
      data_name = WlenCoeff_DataName
      ALLOCATE (L1Birr%suv2%birr_dfld%WlenCoef(nWavelCoef_L1B,nxtrack2_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_irr_file_path, file_path_sz, grp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Birr%suv2%birr_dfld%WlenCoef, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Birr%suv2%birr_dfld%WlenCoef, nWavelCoef_L1B,nxtrack2_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  WavelengthReferenceColumn  ---------------
      data_name = WlenRefCol_DataName
      hdferr = GEMS_Share_Hdf5ReadData(l1b_irr_file_path, file_path_sz, grp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Birr%suv2%birr_dfld%WlenRefCol, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Birr%suv2%birr_dfld%WlenRefCol, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !--------------------------------------
      ! <Sun Volume UV-2 Swath:Swath Attributes>
      !--------------------------------------
      CALL GEMS_Share_MOD_log(llvl, "Sun Volume UV-2 Swath:Swath Attributes",  LOGMSG)

      !-----------------  EarthSunDistance  ---------------
      data_name = EarthSunDist_DataName
      grp_path  = sun_uv2_sgrp_path
      hdferr = GEMS_Share_Hdf5ReadData(l1b_irr_file_path, file_path_sz, grp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Birr%suv2%birr_swatt%ES_Dist, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Birr%suv2%birr_swatt%ES_Dist, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-------------------------------------------------------------------------------------------------------------------------------
      !----------      L1BIRR     ----------!
      !--------------------------------------
      ! <Sun Volume VIS Swath:Data Fields
      !--------------------------------------

      disp_msg = "L1BIRR READ FILE PATH="//l1b_irr_file_path
      CALL GEMS_Share_MOD_log(llvl, disp_msg,  LOGMSG)

      CALL GEMS_Share_MOD_log(llvl, "Sun Volume VIS Swath:Data Fields",  LOGMSG)

      !-----------------  IrradianceExponent ---------------
      IF (ASSOCIATED(L1Birr%svis%birr_dfld%IrrRadExponent)) DEALLOCATE(L1Birr%svis%birr_dfld%IrrRadExponent)
      data_name = IrrExponent_DataName
      grp_path  = sun_vis_dgrp_path
      ALLOCATE (L1Birr%svis%birr_dfld%IrrRadExponent(nWavel3_L1B,nxtrack2_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_irr_file_path, file_path_sz, grp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Birr%svis%birr_dfld%IrrRadExponent, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Birr%svis%birr_dfld%IrrRadExponent, nWavel3_L1B, nxtrack2_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  IrradianceMantissa  ---------------
      IF (ASSOCIATED(L1Birr%svis%birr_dfld%IrrRadMantissa)) DEALLOCATE(L1Birr%svis%birr_dfld%IrrRadMantissa)
      data_name = IrrMantissa_DataName
      ALLOCATE (L1Birr%svis%birr_dfld%IrrRadMantissa(nWavel3_L1B,nxtrack2_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_irr_file_path, file_path_sz, grp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Birr%svis%birr_dfld%IrrRadMantissa, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Birr%svis%birr_dfld%IrrRadMantissa, nWavel3_L1B, nxtrack2_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  PixelQualityFlags  ---------------
      IF (ASSOCIATED(L1Birr%svis%birr_dfld%PQFlag)) DEALLOCATE(L1Birr%svis%birr_dfld%PQFlag)
      data_name = PixelQFlags_DataName
      ALLOCATE (L1Birr%svis%birr_dfld%PQFlag(nWavel3_L1B,nxtrack2_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_irr_file_path, file_path_sz, grp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Birr%svis%birr_dfld%PQFlag, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Birr%svis%birr_dfld%PQFlag, nWavel3_L1B, nxtrack2_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  WavelengthCoefficient  ---------------
      IF (ASSOCIATED(L1Birr%svis%birr_dfld%WlenCoef)) DEALLOCATE(L1Birr%svis%birr_dfld%WlenCoef)
      data_name = WlenCoeff_DataName
      ALLOCATE (L1Birr%svis%birr_dfld%WlenCoef(nWavelCoef_L1B,nxtrack2_L1B), STAT = ierr)
           IF(ierr .NE. 0) THEN
              ctl_msg=data_name
              GOTO 999
           ENDIF
      hdferr = GEMS_Share_Hdf5ReadData(l1b_irr_file_path, file_path_sz, grp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Birr%svis%birr_dfld%WlenCoef, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Birr%svis%birr_dfld%WlenCoef, nWavelCoef_L1B,nxtrack2_L1B, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !-----------------  WavelengthReferenceColumn  ---------------
      data_name = WlenRefCol_DataName
      hdferr = GEMS_Share_Hdf5ReadData(l1b_irr_file_path, file_path_sz, grp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Birr%svis%birr_dfld%WlenRefCol, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Birr%svis%birr_dfld%WlenRefCol, data_name)
      ELSE
            GOTO 2999
      ENDIF

      !--------------------------------------
      ! <Sun Volume VIS Swath:Swath Attributes>
      !--------------------------------------
      CALL GEMS_Share_MOD_log(llvl, "Sun Volume VIS Swath:Swath Attributes",LOGMSG)

      !-----------------  EarthSunDistance  ---------------
      data_name = EarthSunDist_DataName
      grp_path  = sun_vis_sgrp_path
      hdferr = GEMS_Share_Hdf5ReadData(l1b_irr_file_path, file_path_sz,grp_path, grp_path_sz, data_name, data_name_sz,  &
                         L1Birr%svis%birr_swatt%ES_Dist, pstart, pedge)
      IF (hdferr == 0) THEN
            CALL GEMS_Share_MOD_log(llvl, L1Birr%svis%birr_swatt%Es_Dist,data_name)
      ELSE
            GOTO 2999
      ENDIF
      llvl = 0   

      write(disp_msg,'(" OMI HDF5 File READ OK, hdferr=",i6)') hdferr
      CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)

!--------------------------------
!---- END FUNCTION -------------
!--------------------------------
      CALL GEMS_Share_MOD_log(llvl, '---   END L1B IRR READ FUNCTION ---', '  END MSG')
      CALL GEMS_Share_MOD_log(llvl, '.', LOGMSG)
      status = 0

      RETURN 

999   CONTINUE
      llvl = 0   
      write(disp_msg,'( A20, "Memory Allocation Failed, ierr=",i6)') ctl_msg, ierr
      CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)
      status = ierr

      RETURN 

2999  CONTINUE
      llvl = 0   
      write(disp_msg,'(" !!! OMI HDF5 File READ ERROR !!! ", "[", A30, "], hdferr=",i6)') data_name, hdferr
      CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)
      status = hdferr

      RETURN 

      END FUNCTION GEMS_Share_MOD_L1Birr_Read


!-------------------------------------------------------------------------------
!+Description: 
!      L1B파일을 읽기 위해 사용될 전역상수값을 초기화
!
!
! Method:
!
!
! Input files:
!
! Output files:
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2015.01.28   Fisrt Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------
SUBROUTINE GEMS_Share_Init_L1BGlobalConstants

  ! Size for L1B dataset
    nLine_L1B       =    gci_nline      
    nLine_L1B_rug   =    gci_nline_rug  ! -- add chiok 2016.09.29
    nLine_L1B_rvg   =    gci_nline_rvg  ! -- add chiok 2016.09.29     
    nLine_L1B_irr   =    gci_nline_irr  ! -- add chiok 2016.09.29     
    nWavel_L1B      =    gci_nwavel     
    nWavel2_L1B     =    gci_nwavel2    
    nWavel3_L1B     =    gci_nwavel3    
    nxtrack_L1B     =    gci_nxtrack    
    nxtrack2_L1B    =    gci_nxtrack2   
    nWavelCoef_L1B  =    gci_nwavelcoef 

END SUBROUTINE GEMS_Share_Init_L1BGlobalConstants

!-------------------------------------------------------------------------------
!+Description: 
!      L1B파일을 읽기 위해 사용된 메모리를 해제
!
!
! Method:
!
!
! Input files:
!
! Output files:
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2015.01.28   Fisrt Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------
SUBROUTINE GEMS_Share_DeallocL1BVars(L1Brug, L1Brvg, L1Birr)
    TYPE(L1B_rug),  INTENT(IN)      :: L1Brug           ! structure to save the data of L1BRUG File
    TYPE(L1B_rvg),  INTENT(IN)      :: L1Brvg           ! structure to save the data of L1BRVG File
    TYPE(L1B_irr),  INTENT(IN)      :: L1Birr           ! structure to save the data of L1BIRR File

    CALL GEMS_Share_DeallocL1BrugVars(L1Brug)
    CALL GEMS_Share_DeallocL1BrvgVars(L1Brvg)
    CALL GEMS_Share_DeallocL1BirrVars(L1Birr)

    RETURN
END SUBROUTINE GEMS_Share_DeallocL1BVars

SUBROUTINE GEMS_Share_DeallocL1BrugVars(L1Brug)
    TYPE(L1B_rug)       :: L1Brug

    IF (ASSOCIATED(L1Brug%euv1%brug_dfld%ExpTime     )) DEALLOCATE (L1Brug%euv1%brug_dfld%ExpTime      )
    IF (ASSOCIATED(L1Brug%euv1%brug_dfld%ExpType     )) DEALLOCATE (L1Brug%euv1%brug_dfld%ExpType      )
    IF (ASSOCIATED(L1Brug%euv1%brug_dfld%MQFlag      )) DEALLOCATE (L1Brug%euv1%brug_dfld%MQFlag       )
    IF (ASSOCIATED(L1Brug%euv1%brug_dfld%NumSmPixCol )) DEALLOCATE (L1Brug%euv1%brug_dfld%NumSmPixCol  )
    IF (ASSOCIATED(L1Brug%euv1%brug_dfld%PQFlag      )) DEALLOCATE (L1Brug%euv1%brug_dfld%PQFlag       )
    IF (ASSOCIATED(L1Brug%euv1%brug_dfld%RadExponent )) DEALLOCATE (L1Brug%euv1%brug_dfld%RadExponent  )
    IF (ASSOCIATED(L1Brug%euv1%brug_dfld%RadMantissa )) DEALLOCATE (L1Brug%euv1%brug_dfld%RadMantissa  )
    IF (ASSOCIATED(L1Brug%euv1%brug_dfld%RadPrecision)) DEALLOCATE (L1Brug%euv1%brug_dfld%RadPrecision )
    IF (ASSOCIATED(L1Brug%euv1%brug_dfld%ReadTime    )) DEALLOCATE (L1Brug%euv1%brug_dfld%ReadTime     )
    IF (ASSOCIATED(L1Brug%euv1%brug_dfld%WlenCoef    )) DEALLOCATE (L1Brug%euv1%brug_dfld%WlenCoef     )
    IF (ASSOCIATED(L1Brug%euv1%brug_dfld%WlenPrec    )) DEALLOCATE (L1Brug%euv1%brug_dfld%WlenPrec     )
    IF (ASSOCIATED(L1Brug%euv1%brug_dfld%WlenRefCol  )) DEALLOCATE (L1Brug%euv1%brug_dfld%WlenRefCol   )
    IF (ASSOCIATED(L1Brug%euv1%brug_gfld%GPQFlag     )) DEALLOCATE (L1Brug%euv1%brug_gfld%GPQFlag      )
    IF (ASSOCIATED(L1Brug%euv1%brug_gfld%LAT         )) DEALLOCATE (L1Brug%euv1%brug_gfld%LAT          )
    IF (ASSOCIATED(L1Brug%euv1%brug_gfld%LON         )) DEALLOCATE (L1Brug%euv1%brug_gfld%LON          )
    IF (ASSOCIATED(L1Brug%euv1%brug_gfld%SecInDay    )) DEALLOCATE (L1Brug%euv1%brug_gfld%SecInDay     )
    IF (ASSOCIATED(L1Brug%euv1%brug_gfld%SolAziAng   )) DEALLOCATE (L1Brug%euv1%brug_gfld%SolAziAng    )
    IF (ASSOCIATED(L1Brug%euv1%brug_gfld%SolZenAng   )) DEALLOCATE (L1Brug%euv1%brug_gfld%SolZenAng    )
    IF (ASSOCIATED(L1Brug%euv1%brug_gfld%ScAlt       )) DEALLOCATE (L1Brug%euv1%brug_gfld%ScAlt        )
    IF (ASSOCIATED(L1Brug%euv1%brug_gfld%ScLat       )) DEALLOCATE (L1Brug%euv1%brug_gfld%ScLat        )
    IF (ASSOCIATED(L1Brug%euv1%brug_gfld%ScLon       )) DEALLOCATE (L1Brug%euv1%brug_gfld%ScLon        )
    IF (ASSOCIATED(L1Brug%euv1%brug_gfld%TerrHgt     )) DEALLOCATE (L1Brug%euv1%brug_gfld%TerrHgt      )
    IF (ASSOCIATED(L1Brug%euv1%brug_gfld%Time        )) DEALLOCATE (L1Brug%euv1%brug_gfld%Time         )
    IF (ASSOCIATED(L1Brug%euv1%brug_gfld%ViewAziAng  )) DEALLOCATE (L1Brug%euv1%brug_gfld%ViewAziAng   )
    IF (ASSOCIATED(L1Brug%euv1%brug_gfld%ViewZenAng  )) DEALLOCATE (L1Brug%euv1%brug_gfld%ViewZenAng   )
    IF (ASSOCIATED(L1Brug%euv1%brug_gfld%XTQFlag     )) DEALLOCATE (L1Brug%euv1%brug_gfld%XTQFlag      )

    IF (ASSOCIATED(L1Brug%euv2%brug_dfld%ExpTime     )) DEALLOCATE (L1Brug%euv2%brug_dfld%ExpTime      )
    IF (ASSOCIATED(L1Brug%euv2%brug_dfld%ExpType     )) DEALLOCATE (L1Brug%euv2%brug_dfld%ExpType      )
    IF (ASSOCIATED(L1Brug%euv2%brug_dfld%MQFlag      )) DEALLOCATE (L1Brug%euv2%brug_dfld%MQFlag       )
    IF (ASSOCIATED(L1Brug%euv2%brug_dfld%NumSmPixCol )) DEALLOCATE (L1Brug%euv2%brug_dfld%NumSmPixCol  )
    IF (ASSOCIATED(L1Brug%euv2%brug_dfld%PQFlag      )) DEALLOCATE (L1Brug%euv2%brug_dfld%PQFlag       )
    IF (ASSOCIATED(L1Brug%euv2%brug_dfld%RadExponent )) DEALLOCATE (L1Brug%euv2%brug_dfld%RadExponent  )
    IF (ASSOCIATED(L1Brug%euv2%brug_dfld%RadMantissa )) DEALLOCATE (L1Brug%euv2%brug_dfld%RadMantissa  )
    IF (ASSOCIATED(L1Brug%euv2%brug_dfld%RadPrecision)) DEALLOCATE (L1Brug%euv2%brug_dfld%RadPrecision )
    IF (ASSOCIATED(L1Brug%euv2%brug_dfld%ReadTime    )) DEALLOCATE (L1Brug%euv2%brug_dfld%ReadTime     )
    IF (ASSOCIATED(L1Brug%euv2%brug_dfld%WlenCoef    )) DEALLOCATE (L1Brug%euv2%brug_dfld%WlenCoef     )
    IF (ASSOCIATED(L1Brug%euv2%brug_dfld%WlenRefCol  )) DEALLOCATE (L1Brug%euv2%brug_dfld%WlenRefCol   )
    IF (ASSOCIATED(L1Brug%euv2%brug_gfld%GPQFlag     )) DEALLOCATE (L1Brug%euv2%brug_gfld%GPQFlag      )
    IF (ASSOCIATED(L1Brug%euv2%brug_gfld%LAT         )) DEALLOCATE (L1Brug%euv2%brug_gfld%LAT          )
    IF (ASSOCIATED(L1Brug%euv2%brug_gfld%LON         )) DEALLOCATE (L1Brug%euv2%brug_gfld%LON          )
    IF (ASSOCIATED(L1Brug%euv2%brug_gfld%SolAziAng   )) DEALLOCATE (L1Brug%euv2%brug_gfld%SolAziAng    )
    IF (ASSOCIATED(L1Brug%euv2%brug_gfld%SolZenAng   )) DEALLOCATE (L1Brug%euv2%brug_gfld%SolZenAng    )
    IF (ASSOCIATED(L1Brug%euv2%brug_gfld%ScAlt       )) DEALLOCATE (L1Brug%euv2%brug_gfld%ScAlt        )
    IF (ASSOCIATED(L1Brug%euv2%brug_gfld%TerrHgt     )) DEALLOCATE (L1Brug%euv2%brug_gfld%TerrHgt      )
    IF (ASSOCIATED(L1Brug%euv2%brug_gfld%Time        )) DEALLOCATE (L1Brug%euv2%brug_gfld%Time         )
    IF (ASSOCIATED(L1Brug%euv2%brug_gfld%ViewAziAng  )) DEALLOCATE (L1Brug%euv2%brug_gfld%ViewAziAng   )
    IF (ASSOCIATED(L1Brug%euv2%brug_gfld%ViewZenAng  )) DEALLOCATE (L1Brug%euv2%brug_gfld%ViewZenAng   )
    IF (ASSOCIATED(L1Brug%euv2%brug_gfld%XTQFlag     )) DEALLOCATE (L1Brug%euv2%brug_gfld%XTQFlag      )

    NULLIFY (L1Brug%euv1%brug_dfld%ExpTime      )
    NULLIFY (L1Brug%euv1%brug_dfld%ExpType      )
    NULLIFY (L1Brug%euv1%brug_dfld%MQFlag       )
    NULLIFY (L1Brug%euv1%brug_dfld%NumSmPixCol  )
    NULLIFY (L1Brug%euv1%brug_dfld%PQFlag       )
    NULLIFY (L1Brug%euv1%brug_dfld%RadExponent  )
    NULLIFY (L1Brug%euv1%brug_dfld%RadMantissa  )
    NULLIFY (L1Brug%euv1%brug_dfld%RadPrecision )
    NULLIFY (L1Brug%euv1%brug_dfld%ReadTime     )
    NULLIFY (L1Brug%euv1%brug_dfld%WlenCoef     )
    NULLIFY (L1Brug%euv1%brug_dfld%WlenPrec     )
    NULLIFY (L1Brug%euv1%brug_dfld%WlenRefCol   )
    NULLIFY (L1Brug%euv1%brug_gfld%GPQFlag      )
    NULLIFY (L1Brug%euv1%brug_gfld%LAT          )
    NULLIFY (L1Brug%euv1%brug_gfld%LON          )
    NULLIFY (L1Brug%euv1%brug_gfld%SecInDay     )
    NULLIFY (L1Brug%euv1%brug_gfld%SolAziAng    )
    NULLIFY (L1Brug%euv1%brug_gfld%SolZenAng    )
    NULLIFY (L1Brug%euv1%brug_gfld%ScAlt        )
    NULLIFY (L1Brug%euv1%brug_gfld%ScLat        )
    NULLIFY (L1Brug%euv1%brug_gfld%ScLon        )
    NULLIFY (L1Brug%euv1%brug_gfld%TerrHgt      )
    NULLIFY (L1Brug%euv1%brug_gfld%Time         )
    NULLIFY (L1Brug%euv1%brug_gfld%ViewAziAng   )
    NULLIFY (L1Brug%euv1%brug_gfld%ViewZenAng   )
    NULLIFY (L1Brug%euv1%brug_gfld%XTQFlag      )
                                                    
    NULLIFY (L1Brug%euv2%brug_dfld%ExpTime      )
    NULLIFY (L1Brug%euv2%brug_dfld%ExpType      )
    NULLIFY (L1Brug%euv2%brug_dfld%MQFlag       )
    NULLIFY (L1Brug%euv2%brug_dfld%NumSmPixCol  )
    NULLIFY (L1Brug%euv2%brug_dfld%PQFlag       )
    NULLIFY (L1Brug%euv2%brug_dfld%RadExponent  )
    NULLIFY (L1Brug%euv2%brug_dfld%RadMantissa  )
    NULLIFY (L1Brug%euv2%brug_dfld%RadPrecision )
    NULLIFY (L1Brug%euv2%brug_dfld%ReadTime     )
    NULLIFY (L1Brug%euv2%brug_dfld%WlenCoef     )
    NULLIFY (L1Brug%euv2%brug_dfld%WlenRefCol   )
    NULLIFY (L1Brug%euv2%brug_gfld%GPQFlag      )
    NULLIFY (L1Brug%euv2%brug_gfld%LAT          )
    NULLIFY (L1Brug%euv2%brug_gfld%LON          )
    NULLIFY (L1Brug%euv2%brug_gfld%SolAziAng    )
    NULLIFY (L1Brug%euv2%brug_gfld%SolZenAng    )
    NULLIFY (L1Brug%euv2%brug_gfld%ScAlt        )
    NULLIFY (L1Brug%euv2%brug_gfld%Time         )
    NULLIFY (L1Brug%euv2%brug_gfld%TerrHgt      )
    NULLIFY (L1Brug%euv2%brug_gfld%ViewAziAng   )
    NULLIFY (L1Brug%euv2%brug_gfld%ViewZenAng   )
    NULLIFY (L1Brug%euv2%brug_gfld%XTQFlag      )

    write(disp_msg,'("L1Brug Memory deAllocation OK")')
    CALL GEMS_Share_MOD_log(llvl, TRIM(ADJUSTL(disp_msg)), LOGMSG)

    RETURN
END SUBROUTINE GEMS_Share_DeallocL1BrugVars

SUBROUTINE GEMS_Share_DeallocL1BrvgVars(L1Brvg)
    TYPE(L1B_rvg)       :: L1Brvg

    IF (ASSOCIATED(L1Brvg%brvg_dfld%PQFlag           )) DEALLOCATE (L1Brvg%brvg_dfld%PQFlag            )
    IF (ASSOCIATED(L1Brvg%brvg_dfld%RadExponent      )) DEALLOCATE (L1Brvg%brvg_dfld%RadExponent       )
    IF (ASSOCIATED(L1Brvg%brvg_dfld%RadMantissa      )) DEALLOCATE (L1Brvg%brvg_dfld%RadMantissa       )
    IF (ASSOCIATED(L1Brvg%brvg_dfld%WlenCoef         )) DEALLOCATE (L1Brvg%brvg_dfld%WlenCoef          )
    IF (ASSOCIATED(L1Brvg%brvg_dfld%WlenRefCol       )) DEALLOCATE (L1Brvg%brvg_dfld%WlenRefCol        )
    IF (ASSOCIATED(L1Brvg%brvg_gfld%LAT              )) DEALLOCATE (L1Brvg%brvg_gfld%LAT               )
    IF (ASSOCIATED(L1Brvg%brvg_gfld%LON              )) DEALLOCATE (L1Brvg%brvg_gfld%LON               )
    IF (ASSOCIATED(L1Brvg%brvg_gfld%SolAziAng        )) DEALLOCATE (L1Brvg%brvg_gfld%SolAziAng         )
    IF (ASSOCIATED(L1Brvg%brvg_gfld%SolZenAng        )) DEALLOCATE (L1Brvg%brvg_gfld%SolZenAng         )
    IF (ASSOCIATED(L1Brvg%brvg_gfld%TerrHgt          )) DEALLOCATE (L1Brvg%brvg_gfld%TerrHgt           )
    IF (ASSOCIATED(L1Brvg%brvg_gfld%ViewAziAng       )) DEALLOCATE (L1Brvg%brvg_gfld%ViewAziAng        )
    IF (ASSOCIATED(L1Brvg%brvg_gfld%ViewZenAng       )) DEALLOCATE (L1Brvg%brvg_gfld%ViewZenAng        )

    NULLIFY (L1Brvg%brvg_dfld%PQFlag            )
    NULLIFY (L1Brvg%brvg_dfld%RadExponent       )
    NULLIFY (L1Brvg%brvg_dfld%RadMantissa       )
    NULLIFY (L1Brvg%brvg_dfld%WlenCoef          )
    NULLIFY (L1Brvg%brvg_dfld%WlenRefCol        )
    NULLIFY (L1Brvg%brvg_gfld%LAT               )
    NULLIFY (L1Brvg%brvg_gfld%LON               )
    NULLIFY (L1Brvg%brvg_gfld%SolAziAng         )
    NULLIFY (L1Brvg%brvg_gfld%SolZenAng         )
    NULLIFY (L1Brvg%brvg_gfld%TerrHgt           )
    NULLIFY (L1Brvg%brvg_gfld%ViewAziAng        )
    NULLIFY (L1Brvg%brvg_gfld%ViewZenAng        )

    write(disp_msg,'("L1Brvg Memory deAllocation OK")')
    CALL GEMS_Share_MOD_log(llvl, TRIM(ADJUSTL(disp_msg)), LOGMSG)

    RETURN
END SUBROUTINE GEMS_Share_DeallocL1BrvgVars

SUBROUTINE GEMS_Share_DeallocL1BirrVars(L1Birr)
    TYPE(L1B_irr)       :: L1Birr

    IF (ASSOCIATED(L1Birr%suv1%birr_dfld%IrrRadExponent   )) DEALLOCATE (L1Birr%suv1%birr_dfld%IrrRadExponent    )
    IF (ASSOCIATED(L1Birr%suv1%birr_dfld%IrrRadMantissa   )) DEALLOCATE (L1Birr%suv1%birr_dfld%IrrRadMantissa    )
    IF (ASSOCIATED(L1Birr%suv1%birr_dfld%WlenCoef         )) DEALLOCATE (L1Birr%suv1%birr_dfld%WlenCoef          )

    IF (ASSOCIATED(L1Birr%suv2%birr_dfld%IrrRadExponent   )) DEALLOCATE (L1Birr%suv2%birr_dfld%IrrRadExponent    )
    IF (ASSOCIATED(L1Birr%suv2%birr_dfld%IrrRadMantissa   )) DEALLOCATE (L1Birr%suv2%birr_dfld%IrrRadMantissa    )
    IF (ASSOCIATED(L1Birr%suv2%birr_dfld%IrrRadPrecision  )) DEALLOCATE (L1Birr%suv2%birr_dfld%IrrRadPrecision   )
    IF (ASSOCIATED(L1Birr%suv2%birr_dfld%PQFlag           )) DEALLOCATE (L1Birr%suv2%birr_dfld%PQFlag            )
    IF (ASSOCIATED(L1Birr%suv2%birr_dfld%WlenCoef         )) DEALLOCATE (L1Birr%suv2%birr_dfld%WlenCoef          )

    IF (ASSOCIATED(L1Birr%svis%birr_dfld%IrrRadExponent   )) DEALLOCATE (L1Birr%svis%birr_dfld%IrrRadExponent    )
    IF (ASSOCIATED(L1Birr%svis%birr_dfld%IrrRadMantissa   )) DEALLOCATE (L1Birr%svis%birr_dfld%IrrRadMantissa    )
    IF (ASSOCIATED(L1Birr%svis%birr_dfld%IrrRadPrecision  )) DEALLOCATE (L1Birr%svis%birr_dfld%IrrRadPrecision)
    IF (ASSOCIATED(L1Birr%svis%birr_dfld%PQFlag           )) DEALLOCATE (L1Birr%svis%birr_dfld%PQFlag            )
    IF (ASSOCIATED(L1Birr%svis%birr_dfld%WlenCoef         )) DEALLOCATE (L1Birr%svis%birr_dfld%WlenCoef          )

    NULLIFY (L1Birr%suv1%birr_dfld%IrrRadExponent    ) 
    NULLIFY (L1Birr%suv1%birr_dfld%IrrRadMantissa    ) 
    NULLIFY (L1Birr%suv1%birr_dfld%WlenCoef          ) 
                                                          
    NULLIFY (L1Birr%suv2%birr_dfld%IrrRadExponent    ) 
    NULLIFY (L1Birr%suv2%birr_dfld%IrrRadMantissa    ) 
    NULLIFY (L1Birr%suv2%birr_dfld%IrrRadPrecision   ) 
    NULLIFY (L1Birr%suv2%birr_dfld%PQFlag            ) 
    NULLIFY (L1Birr%suv2%birr_dfld%WlenCoef          ) 
                                                          
    NULLIFY (L1Birr%svis%birr_dfld%IrrRadExponent    ) 
    NULLIFY (L1Birr%svis%birr_dfld%IrrRadMantissa    ) 
    NULLIFY (L1Birr%svis%birr_dfld%IrrRadPrecision   )    
    NULLIFY (L1Birr%svis%birr_dfld%PQFlag            ) 
    NULLIFY (L1Birr%svis%birr_dfld%WlenCoef          ) 

    write(disp_msg,'("L1Birr Memory deAllocation OK")')
    CALL GEMS_Share_MOD_log(llvl, TRIM(ADJUSTL(disp_msg)), LOGMSG)

    RETURN
END SUBROUTINE GEMS_Share_DeallocL1BirrVars

!-------------------------------------------------------------------------------
!+Description: 
!      L1B파일을 읽기 위해 사용된 메모리를 할당
!
!
! Method:
!
!
! Input files:
!
! Output files:
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2015.04.01   Fisrt Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------
SUBROUTINE GEMS_Share_AllocL1BVars(L1Brug, L1Brvg, L1Birr)
    TYPE(L1B_rug),  INTENT(INOUT)      :: L1Brug           ! structure to save the data of L1BRUG File
    TYPE(L1B_rvg),  INTENT(INOUT)      :: L1Brvg           ! structure to save the data of L1BRVG File
    TYPE(L1B_irr),  INTENT(INOUT)      :: L1Birr           ! structure to save the data of L1BIRR File

    CALL GEMS_Share_AllocL1BrugVars
    CALL GEMS_Share_AllocL1BrvgVars
    CALL GEMS_Share_AllocL1BirrVars

    RETURN
END SUBROUTINE GEMS_Share_AllocL1BVars

SUBROUTINE GEMS_Share_AllocL1BrugVars

    INTEGER                         :: ierr        ! Error code 
    CHARACTER*200                   :: data_name = "Earth UV1 Data Field"

    !--------------------------------------
    ! <Earth UV-1 Swath:Data Fields>
    !--------------------------------------
    !-----------------  ExposureTime  ---------------
    IF (ASSOCIATED(gds_brug%euv1%brug_dfld%ExpTime)) DEALLOCATE(gds_brug%euv1%brug_dfld%ExpTime)
    ALLOCATE (gds_brug%euv1%brug_dfld%ExpTime(nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  ExposureType  ---------------
    IF (ASSOCIATED(gds_brug%euv1%brug_dfld%ExpType)) DEALLOCATE(gds_brug%euv1%brug_dfld%ExpType)
    ALLOCATE (gds_brug%euv1%brug_dfld%ExpType(nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  MeasurementQualityFlags  ---------------
    IF (ASSOCIATED(gds_brug%euv1%brug_dfld%MQFlag)) DEALLOCATE(gds_brug%euv1%brug_dfld%MQFlag)
    ALLOCATE (gds_brug%euv1%brug_dfld%MQFlag(nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  NumberSmallPixelColumns  ---------------
    IF (ASSOCIATED(gds_brug%euv1%brug_dfld%NumSmPixCol)) DEALLOCATE(gds_brug%euv1%brug_dfld%NumSmPixCol)
    ALLOCATE (gds_brug%euv1%brug_dfld%NumSmPixCol(nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  PixelQualityFlags  ---------------
    IF (ASSOCIATED(gds_brug%euv1%brug_dfld%PQFlag)) DEALLOCATE(gds_brug%euv1%brug_dfld%PQFlag)
    ALLOCATE (gds_brug%euv1%brug_dfld%PQFlag(nWavel_L1B,nxtrack_L1B,nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  RadianceExponent ---------------
    IF (ASSOCIATED(gds_brug%euv1%brug_dfld%RadExponent)) DEALLOCATE(gds_brug%euv1%brug_dfld%RadExponent)
    ALLOCATE (gds_brug%euv1%brug_dfld%RadExponent(nWavel_L1B,nxtrack_L1B,nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  RadianceMantissa  ---------------
    IF (ASSOCIATED(gds_brug%euv1%brug_dfld%RadMantissa)) DEALLOCATE(gds_brug%euv1%brug_dfld%RadMantissa)
    ALLOCATE (gds_brug%euv1%brug_dfld%RadMantissa(nWavel_L1B,nxtrack_L1B,nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  RadiancePrecisionMantissa  ---------------
    IF (ASSOCIATED(gds_brug%euv1%brug_dfld%RadPrecision)) DEALLOCATE(gds_brug%euv1%brug_dfld%RadPrecision)
    ALLOCATE (gds_brug%euv1%brug_dfld%RadPrecision(nWavel_L1B,nxtrack_L1B,nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  ReadoutTime  ---------------
    IF (ASSOCIATED(gds_brug%euv1%brug_dfld%ReadTime)) DEALLOCATE(gds_brug%euv1%brug_dfld%ReadTime)
    ALLOCATE (gds_brug%euv1%brug_dfld%ReadTime(nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  WavelengthCoefficient  ---------------
    IF (ASSOCIATED(gds_brug%euv1%brug_dfld%WlenCoef)) DEALLOCATE(gds_brug%euv1%brug_dfld%WlenCoef)
    ALLOCATE (gds_brug%euv1%brug_dfld%WlenCoef(nWavelCoef_L1B,nxtrack_L1B,nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !---------------  WavelengthCoefficientPrecision  -----------------
    IF (ASSOCIATED(gds_brug%euv1%brug_dfld%WlenPrec)) DEALLOCATE(gds_brug%euv1%brug_dfld%WlenPrec)
    ALLOCATE (gds_brug%euv1%brug_dfld%WlenPrec(nWavelCoef_L1B,nxtrack_L1B,nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  WavelengthReferenceColumn  ---------------
    IF (ASSOCIATED(gds_brug%euv1%brug_dfld%WlenRefCol)) DEALLOCATE(gds_brug%euv1%brug_dfld%WlenRefCol)
    ALLOCATE (gds_brug%euv1%brug_dfld%WlenRefCol(nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !--------------------------------------
    ! <Earth UV-1 Swath/Geolocation Fields>
    !--------------------------------------
    data_name = "Earth UV1 Geolocation Field"
    !-----------------  GroundPixelQualityFlags  ---------------
    IF (ASSOCIATED(gds_brug%euv1%brug_gfld%GPQFlag)) DEALLOCATE(gds_brug%euv1%brug_gfld%GPQFlag)
    ALLOCATE (gds_brug%euv1%brug_gfld%GPQFlag(nxtrack_L1B, nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  Latitude  ---------------
    IF (ASSOCIATED(gds_brug%euv1%brug_gfld%LAT)) DEALLOCATE(gds_brug%euv1%brug_gfld%LAT)
    ALLOCATE (gds_brug%euv1%brug_gfld%LAT(nxtrack_L1B, nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  Longitude  ---------------
    IF (ASSOCIATED(gds_brug%euv1%brug_gfld%LON)) DEALLOCATE(gds_brug%euv1%brug_gfld%LON)
    ALLOCATE (gds_brug%euv1%brug_gfld%LON(nxtrack_L1B, nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  SecondsInDay  ---------------
    IF (ASSOCIATED(gds_brug%euv1%brug_gfld%SecInDay)) DEALLOCATE(gds_brug%euv1%brug_gfld%SecInDay)
    ALLOCATE (gds_brug%euv1%brug_gfld%SecInDay(nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  SolarAzimuthAngle  ---------------
    IF (ASSOCIATED(gds_brug%euv1%brug_gfld%SolAziAng)) DEALLOCATE(gds_brug%euv1%brug_gfld%SolAziAng)
    ALLOCATE (gds_brug%euv1%brug_gfld%SolAziAng(nxtrack_L1B, nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF
    
    !-----------------  SolarZenithAngle  ---------------
    IF (ASSOCIATED(gds_brug%euv1%brug_gfld%SolZenAng)) DEALLOCATE(gds_brug%euv1%brug_gfld%SolZenAng)
    ALLOCATE (gds_brug%euv1%brug_gfld%SolZenAng(nxtrack_L1B, nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF
    
    !-----------------  SpacecraftAltitude  ---------------
    IF (ASSOCIATED(gds_brug%euv1%brug_gfld%ScAlt)) DEALLOCATE(gds_brug%euv1%brug_gfld%ScAlt)
    ALLOCATE (gds_brug%euv1%brug_gfld%ScAlt(nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  SpacecraftLatitude  ---------------
    IF (ASSOCIATED(gds_brug%euv1%brug_gfld%ScLat)) DEALLOCATE(gds_brug%euv1%brug_gfld%ScLat)
    ALLOCATE (gds_brug%euv1%brug_gfld%ScLat(nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  SpacecraftLongitude  ---------------
    IF (ASSOCIATED(gds_brug%euv1%brug_gfld%ScLon)) DEALLOCATE(gds_brug%euv1%brug_gfld%ScLon)
    ALLOCATE (gds_brug%euv1%brug_gfld%ScLon(nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  TerrainHeight  ---------------
    IF (ASSOCIATED(gds_brug%euv1%brug_gfld%TerrHgt)) DEALLOCATE(gds_brug%euv1%brug_gfld%TerrHgt)
    ALLOCATE (gds_brug%euv1%brug_gfld%TerrHgt(nxtrack_L1B, nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  Time  ---------------
    IF (ASSOCIATED(gds_brug%euv1%brug_gfld%Time)) DEALLOCATE(gds_brug%euv1%brug_gfld%Time)
    ALLOCATE (gds_brug%euv1%brug_gfld%Time(nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  ViewingAzimuthAngle  ---------------
    IF (ASSOCIATED(gds_brug%euv1%brug_gfld%ViewAziAng)) DEALLOCATE(gds_brug%euv1%brug_gfld%ViewAziAng)
    ALLOCATE (gds_brug%euv1%brug_gfld%ViewAziAng(nxtrack_L1B, nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  ViewingZenithAngle  ---------------
    IF (ASSOCIATED(gds_brug%euv1%brug_gfld%ViewZenAng)) DEALLOCATE(gds_brug%euv1%brug_gfld%ViewZenAng)
    ALLOCATE (gds_brug%euv1%brug_gfld%ViewZenAng(nxtrack_L1B, nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  XTrackQualityFlags   ---------------
    IF (ASSOCIATED(gds_brug%euv1%brug_gfld%XTQFlag)) DEALLOCATE(gds_brug%euv1%brug_gfld%XTQFlag)
    ALLOCATE (gds_brug%euv1%brug_gfld%XTQFlag(nxtrack_L1B, nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF


    !----------L1BRUG----------!
    ! <Earth UV-2 Swath:Data Fileds>
    data_name = "Earth UV2 Data Field"

    !-----------------  ExposureTime  ---------------
    IF (ASSOCIATED(gds_brug%euv2%brug_dfld%ExpTime)) DEALLOCATE(gds_brug%euv2%brug_dfld%ExpTime)
    ALLOCATE (gds_brug%euv2%brug_dfld%ExpTime(nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  ExposureType  ---------------
    IF (ASSOCIATED(gds_brug%euv2%brug_dfld%ExpType)) DEALLOCATE(gds_brug%euv2%brug_dfld%ExpType)
    ALLOCATE (gds_brug%euv2%brug_dfld%ExpType(nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  MeasurementQualityFlags  ---------------
    IF (ASSOCIATED(gds_brug%euv2%brug_dfld%MQFlag)) DEALLOCATE(gds_brug%euv2%brug_dfld%MQFlag)
    ALLOCATE (gds_brug%euv2%brug_dfld%MQFlag(nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  NumberSmallPixelColumns  ---------------
    IF (ASSOCIATED(gds_brug%euv2%brug_dfld%NumSmPixCol)) DEALLOCATE(gds_brug%euv2%brug_dfld%NumSmPixCol)
    ALLOCATE (gds_brug%euv2%brug_dfld%NumSmPixCol(nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  PixelQualityFlags  ---------------
    IF (ASSOCIATED(gds_brug%euv2%brug_dfld%PQFlag))      DEALLOCATE(gds_brug%euv2%brug_dfld%PQFlag)
    ALLOCATE (gds_brug%euv2%brug_dfld%PQFlag(nWavel2_L1B,nxtrack2_L1B,nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  RadianceExponent ---------------
    IF (ASSOCIATED(gds_brug%euv2%brug_dfld%RadExponent)) DEALLOCATE(gds_brug%euv2%brug_dfld%RadExponent)
    ALLOCATE (gds_brug%euv2%brug_dfld%RadExponent(nWavel2_L1B,nxtrack2_L1B,nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  RadianceMantissa  ---------------
    IF (ASSOCIATED(gds_brug%euv2%brug_dfld%RadMantissa)) DEALLOCATE(gds_brug%euv2%brug_dfld%RadMantissa)
    ALLOCATE (gds_brug%euv2%brug_dfld%RadMantissa(nWavel2_L1B,nxtrack2_L1B,nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

   !-----------------  RadiancePrecisionMantissa  ---------------
    IF (ASSOCIATED(gds_brug%euv2%brug_dfld%RadPrecision)) DEALLOCATE(gds_brug%euv2%brug_dfld%RadPrecision)
    ALLOCATE (gds_brug%euv2%brug_dfld%RadPrecision(nWavel2_L1B,nxtrack2_L1B,nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  ReadoutTime  ---------------
    IF (ASSOCIATED(gds_brug%euv2%brug_dfld%ReadTime)) DEALLOCATE(gds_brug%euv2%brug_dfld%ReadTime)
    ALLOCATE (gds_brug%euv2%brug_dfld%ReadTime(nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  WavelengthCoefficient  ---------------
    IF (ASSOCIATED(gds_brug%euv2%brug_dfld%WlenCoef)) DEALLOCATE(gds_brug%euv2%brug_dfld%WlenCoef)
    ALLOCATE (gds_brug%euv2%brug_dfld%WlenCoef(nWavelCoef_L1B,nxtrack2_L1B,nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  WavelengthReferenceColumn  ---------------
    IF (ASSOCIATED(gds_brug%euv2%brug_dfld%WlenRefCol)) DEALLOCATE(gds_brug%euv2%brug_dfld%WlenRefCol)
    ALLOCATE (gds_brug%euv2%brug_dfld%WlenRefCol(nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !--------------------------------------
    ! <Earth UV-2 Swath/Geolocation Fields>
    !--------------------------------------
    data_name = "Earth UV2 Geolocation Field"

    !-----------------  GroundPixelQualityFlags  ---------------
    IF (ASSOCIATED(gds_brug%euv2%brug_gfld%GPQFlag)) DEALLOCATE(gds_brug%euv2%brug_gfld%GPQFlag)
    ALLOCATE (gds_brug%euv2%brug_gfld%GPQFlag(nxtrack2_L1B, nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  Latitude  ---------------
    IF (ASSOCIATED(gds_brug%euv2%brug_gfld%LAT)) DEALLOCATE(gds_brug%euv2%brug_gfld%LAT)
    ALLOCATE (gds_brug%euv2%brug_gfld%LAT(nxtrack2_L1B, nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  Longitude  ---------------
    IF (ASSOCIATED(gds_brug%euv2%brug_gfld%LON)) DEALLOCATE(gds_brug%euv2%brug_gfld%LON)
    ALLOCATE (gds_brug%euv2%brug_gfld%LON(nxtrack2_L1B, nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  SolarAzimuthAngle  ---------------
    IF (ASSOCIATED(gds_brug%euv2%brug_gfld%SolAziAng)) DEALLOCATE(gds_brug%euv2%brug_gfld%SolAziAng)
    ALLOCATE (gds_brug%euv2%brug_gfld%SolAziAng(nxtrack2_L1B, nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF
    
    !-----------------  SolarZenithAngle  ---------------
    IF (ASSOCIATED(gds_brug%euv2%brug_gfld%SolZenAng)) DEALLOCATE(gds_brug%euv2%brug_gfld%SolZenAng)
    ALLOCATE (gds_brug%euv2%brug_gfld%SolZenAng(nxtrack2_L1B, nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  SpacecraftAltitude  ---------------
    IF (ASSOCIATED(gds_brug%euv2%brug_gfld%ScAlt)) DEALLOCATE(gds_brug%euv2%brug_gfld%ScAlt)
    ALLOCATE (gds_brug%euv2%brug_gfld%ScAlt(nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF
 
    !-----------------  TerrainHeight  ---------------
    IF (ASSOCIATED(gds_brug%euv2%brug_gfld%TerrHgt)) DEALLOCATE(gds_brug%euv2%brug_gfld%TerrHgt)
    ALLOCATE (gds_brug%euv2%brug_gfld%TerrHgt(nxtrack2_L1B, nLine_L1B), STAT=ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  Time  ---------------
    IF (ASSOCIATED(gds_brug%euv2%brug_gfld%Time)) DEALLOCATE(gds_brug%euv2%brug_gfld%Time)
    ALLOCATE (gds_brug%euv2%brug_gfld%Time(nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  ViewingAzimuthAngle  ---------------
    IF (ASSOCIATED(gds_brug%euv2%brug_gfld%ViewAziAng)) DEALLOCATE(gds_brug%euv2%brug_gfld%ViewAziAng)
    ALLOCATE (gds_brug%euv2%brug_gfld%ViewAziAng(nxtrack2_L1B, nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  ViewingZenithAngle  ---------------
    IF (ASSOCIATED(gds_brug%euv2%brug_gfld%ViewZenAng)) DEALLOCATE(gds_brug%euv2%brug_gfld%ViewZenAng)
    ALLOCATE (gds_brug%euv2%brug_gfld%ViewZenAng(nxtrack2_L1B, nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  XTrackQualityFlags   ---------------
    IF (ASSOCIATED(gds_brug%euv2%brug_gfld%XTQFlag)) DEALLOCATE(gds_brug%euv2%brug_gfld%XTQFlag)
    ALLOCATE (gds_brug%euv2%brug_gfld%XTQFlag(nxtrack2_L1B, nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !----------L1BRUG----------!
    ! <Earth UV-2 Swath: Swath Attributes>

    write(disp_msg,'(" Memory Allocation OK, errCode=",i6)') ierr
    CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)

    RETURN

999 CONTINUE
    write(disp_msg,'( A20, "L1BRUG Memory Allocation Failed, ierr=",i6)') ctl_msg, ierr
    CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)

    RETURN 

END SUBROUTINE GEMS_Share_AllocL1BrugVars

SUBROUTINE GEMS_Share_AllocL1BrvgVars

    INTEGER                         :: ierr        ! Error code 
    CHARACTER*200                   :: data_name = "Earth VIS Data Field"

    !--------------------------------------
    ! <Earth VIS Swath:Data Fields>
    !--------------------------------------
    !-----------------  PixelQualityFlags ---------------
    IF (ASSOCIATED(gds_brvg%brvg_dfld%PQFlag)) DEALLOCATE(gds_brvg%brvg_dfld%PQFlag)
    ALLOCATE(gds_brvg%brvg_dfld%PQFlag(nWavel3_L1B,nxtrack2_L1B,nLine_L1B), STAT= ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF


    !-----------------  RadianceExponent ---------------
    IF (ASSOCIATED(gds_brvg%brvg_dfld%RadExponent)) DEALLOCATE(gds_brvg%brvg_dfld%RadExponent)
    ALLOCATE (gds_brvg%brvg_dfld%RadExponent(nWavel3_L1B,nxtrack2_L1B,nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  RadianceMantissa  ---------------
    IF (ASSOCIATED(gds_brvg%brvg_dfld%RadMantissa)) DEALLOCATE(gds_brvg%brvg_dfld%RadMantissa)
    ALLOCATE (gds_brvg%brvg_dfld%RadMantissa(nWavel3_L1B,nxtrack2_L1B,nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  WavelengthCoefficient  ---------------
    IF (ASSOCIATED(gds_brvg%brvg_dfld%WlenCoef)) DEALLOCATE(gds_brvg%brvg_dfld%WlenCoef)
    ALLOCATE (gds_brvg%brvg_dfld%WlenCoef(nWavelCoef_L1B,nxtrack2_L1B,nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  WavelengthReferenceColumn  ---------------
    IF (ASSOCIATED(gds_brvg%brvg_dfld%WlenRefCol)) DEALLOCATE(gds_brvg%brvg_dfld%WlenRefCol)
    ALLOCATE (gds_brvg%brvg_dfld%WlenRefCol(nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !--------------------------------------
    ! <Earth VIS Swath/Geolocation Fields>
    !--------------------------------------
    data_name = "Earth VIS Geolocation Field"
    !-----------------  Latitude  ---------------
    IF (ASSOCIATED(gds_brvg%brvg_gfld%LAT)) DEALLOCATE(gds_brvg%brvg_gfld%LAT)
    ALLOCATE (gds_brvg%brvg_gfld%LAT(nxtrack2_L1B, nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  Longitude  ---------------
    IF (ASSOCIATED(gds_brvg%brvg_gfld%LON)) DEALLOCATE(gds_brvg%brvg_gfld%LON)
    ALLOCATE (gds_brvg%brvg_gfld%LON(nxtrack2_L1B, nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  SolarAzimuthAngle  ---------------
    IF (ASSOCIATED(gds_brvg%brvg_gfld%SolAziAng)) DEALLOCATE(gds_brvg%brvg_gfld%SolAziAng)
    ALLOCATE (gds_brvg%brvg_gfld%SolAziAng(nxtrack2_L1B, nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF
    
    !-----------------  SolarZenithAngle  ---------------
    IF (ASSOCIATED(gds_brvg%brvg_gfld%SolZenAng)) DEALLOCATE(gds_brvg%brvg_gfld%SolZenAng)
    ALLOCATE (gds_brvg%brvg_gfld%SolZenAng(nxtrack2_L1B, nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  TerrainHeight  ---------------
    IF (ASSOCIATED(gds_brvg%brvg_gfld%TerrHgt)) DEALLOCATE(gds_brvg%brvg_gfld%TerrHgt)
    ALLOCATE (gds_brvg%brvg_gfld%TerrHgt(nxtrack2_L1B,nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF
    
    !-----------------  ViewingAzimuthAngle  ---------------
    IF (ASSOCIATED(gds_brvg%brvg_gfld%ViewAziAng)) DEALLOCATE(gds_brvg%brvg_gfld%ViewAziAng)
    ALLOCATE (gds_brvg%brvg_gfld%ViewAziAng(nxtrack2_L1B, nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  ViewingZenithAngle  ---------------
    IF (ASSOCIATED(gds_brvg%brvg_gfld%ViewZenAng)) DEALLOCATE(gds_brvg%brvg_gfld%ViewZenAng)
    ALLOCATE (gds_brvg%brvg_gfld%ViewZenAng(nxtrack2_L1B, nLine_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    write(disp_msg,'(" Memory Allocation OK, errCode=",i6)') ierr
    CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)

    RETURN 

999 CONTINUE
    write(disp_msg,'( A20, "Memory Allocation Failed, ierr=",i6)') ctl_msg, ierr
    CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)

    RETURN
END SUBROUTINE GEMS_Share_AllocL1BrvgVars

SUBROUTINE GEMS_Share_AllocL1BirrVars

    INTEGER                         :: ierr        ! Error code 
    CHARACTER*200                   :: data_name = "SUN UV1 Data Field"

    !--------------------------------------
    ! <Sun Volume UV-1 Swath:Data Fields>
    !--------------------------------------
    !-----------------  IrradianceExponent ---------------
    IF (ASSOCIATED(gds_birr%suv1%birr_dfld%IrrRadExponent)) DEALLOCATE(gds_birr%suv1%birr_dfld%IrrRadExponent)
    ALLOCATE (gds_birr%suv1%birr_dfld%IrrRadExponent(nWavel_L1B,nxtrack_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  IrradianceMantissa  ---------------
    IF (ASSOCIATED(gds_birr%suv1%birr_dfld%IrrRadMantissa)) DEALLOCATE(gds_birr%suv1%birr_dfld%IrrRadMantissa)
    ALLOCATE (gds_birr%suv1%birr_dfld%IrrRadMantissa(nWavel_L1B,nxtrack_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  WavelengthCoefficient  ---------------
    IF (ASSOCIATED(gds_birr%suv1%birr_dfld%WlenCoef)) DEALLOCATE(gds_birr%suv1%birr_dfld%WlenCoef)
    ALLOCATE (gds_birr%suv1%birr_dfld%WlenCoef(nWavelCoef_L1B,nxtrack_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  WavelengthReferenceColumn  ---------------

    !--------------------------------------
    ! <Sun Volume UV-2 Swath: Swath Attributes>
    !--------------------------------------

    !--------------------------------------
    !--------------------------------------
    ! <Sun Volume UV-2 Swath:Data Fields>
    !--------------------------------------
    data_name = "SUN UV2 Data Field"

    !-----------------  IrradianceExponent ---------------
    IF (ASSOCIATED(gds_birr%suv2%birr_dfld%IrrRadExponent)) DEALLOCATE(gds_birr%suv2%birr_dfld%IrrRadExponent)
    ALLOCATE (gds_birr%suv2%birr_dfld%IrrRadExponent(nWavel2_L1B,nxtrack2_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  IrradianceMantissa  ---------------
    IF (ASSOCIATED(gds_birr%suv2%birr_dfld%IrrRadMantissa)) DEALLOCATE(gds_birr%suv2%birr_dfld%IrrRadMantissa)
    ALLOCATE (gds_birr%suv2%birr_dfld%IrrRadMantissa(nWavel2_L1B,nxtrack2_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  IrradiancePrecisionMantissa  ---------------
    IF (ASSOCIATED(gds_birr%suv2%birr_dfld%IrrRadPrecision)) DEALLOCATE(gds_birr%suv2%birr_dfld%IrrRadPrecision)
    ALLOCATE (gds_birr%suv2%birr_dfld%IrrRadPrecision(nWavel2_L1B,nxtrack2_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  PixelQualityFlags  ---------------
    IF (ASSOCIATED(gds_birr%suv2%birr_dfld%PQFlag)) DEALLOCATE(gds_birr%suv2%birr_dfld%PQFlag)
    ALLOCATE (gds_birr%suv2%birr_dfld%PQFlag(nWavel2_L1B,nxtrack2_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  WavelengthCoefficient  ---------------
    IF (ASSOCIATED(gds_birr%suv2%birr_dfld%WlenCoef)) DEALLOCATE(gds_birr%suv2%birr_dfld%WlenCoef)
    ALLOCATE (gds_birr%suv2%birr_dfld%WlenCoef(nWavelCoef_L1B,nxtrack2_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  WavelengthReferenceColumn  ---------------

    !--------------------------------------
    ! <Sun Volume UV-2 Swath: Swath Attributes>
    !--------------------------------------

    !--------------------------------------
    ! <Sun Volume VIS Swath:Data Fields
    !--------------------------------------
    data_name = "SUN VIS Data Field"
    !-----------------  IrradianceExponent ---------------
    IF (ASSOCIATED(gds_birr%svis%birr_dfld%IrrRadExponent)) DEALLOCATE(gds_birr%svis%birr_dfld%IrrRadExponent)
    ALLOCATE (gds_birr%svis%birr_dfld%IrrRadExponent(nWavel3_L1B,nxtrack2_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  IrradianceMantissa  ---------------
    IF (ASSOCIATED(gds_birr%svis%birr_dfld%IrrRadMantissa)) DEALLOCATE(gds_birr%svis%birr_dfld%IrrRadMantissa)
    ALLOCATE (gds_birr%svis%birr_dfld%IrrRadMantissa(nWavel3_L1B,nxtrack2_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  PixelQualityFlags  ---------------
    IF (ASSOCIATED(gds_birr%svis%birr_dfld%PQFlag)) DEALLOCATE(gds_birr%svis%birr_dfld%PQFlag)
    ALLOCATE (gds_birr%svis%birr_dfld%PQFlag(nWavel3_L1B,nxtrack2_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  WavelengthCoefficient  ---------------
    IF (ASSOCIATED(gds_birr%svis%birr_dfld%WlenCoef)) DEALLOCATE(gds_birr%svis%birr_dfld%WlenCoef)
    ALLOCATE (gds_birr%svis%birr_dfld%WlenCoef(nWavelCoef_L1B,nxtrack2_L1B), STAT = ierr)
         IF(ierr .NE. 0) THEN
            ctl_msg=data_name
            GOTO 999
         ENDIF

    !-----------------  WavelengthReferenceColumn  ---------------

    write(disp_msg,'(" Memory Allocation OK, errCode=",i6)') ierr
    CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)

    RETURN 

999 CONTINUE
    write(disp_msg,'( A20, "Memory Allocation Failed, ierr=",i6)') ctl_msg, ierr
    CALL GEMS_Share_MOD_log(llvl, disp_msg, LOGMSG)

    RETURN
END SUBROUTINE GEMS_Share_AllocL1BirrVars

!-------------------------------------------------------------------------------
!+Description: 
!      SUBROUTINE GEMS_Share_GetL1BInfo_OMI.
!
!
! Method:
!
!
! Input files:
!       l1b_fpath   : File Class
!
!
! Output files:
!       obit        : orbit number of OMI
!       cly         : L1B data classificaion 
!       date        : granule start time and production date
!
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2016.09.09  Fisrt Code (ChiOk An, Seasoft) 
!-------------------------------------------------------------------------------
SUBROUTINE GEMS_Share_GetL1BInfo_OMI (                 &
                                         l1b_fpath,    &
                                         obit,         &
                                         cly,          &
                                         date          &
                                     )

    CHARACTER(LEN=*), INTENT(IN)      :: l1b_fpath
    CHARACTER(LEN=*), INTENT(OUT)     :: obit
    CHARACTER(LEN=*), INTENT(OUT)     :: cly
    CHARACTER(LEN=*), INTENT(OUT)     :: date(2)
!--- local variable
    CHARACTER(LEN=4)        :: year(2)
    CHARACTER(LEN=2)        :: mon(2)
    CHARACTER(LEN=2)        :: day(2)
    CHARACTER(LEN=4)        :: time(2)
    CHARACTER(LEN=1)        :: tmp2
    CHARACTER(LEN=100)      :: tmp3
    INTEGER                 :: i, fstr, estr,cont

    fstr = 1
    estr = 0
    cont = 0
    do i = 1, 100
       tmp2 = l1b_fpath(i:i) 
       if(tmp2 .eq. '/' .or. tmp2 .eq. '-' .or. tmp2 .eq. '_'       &
          .or. tmp2 .eq. '.'                                    ) then      ! the point cond needs

          estr = i-1
          tmp3 = l1b_fpath(fstr:estr)

          if( tmp3(1:1) .eq. 'o') obit = tmp3(2:6) ! obit number 
          if( tmp3(1:1) .eq. 'O') cly  = tmp3(5:8) ! date classify
          if( tmp3(1:1) .eq. '2') then
              cont = cont + 1
              year(cont) = tmp3(1:4) ! year
               mon(cont) = tmp3(6:7) ! month
               day(cont) = tmp3(8:9) ! day
              time(cont) = tmp3(11:14) ! time
              if ( cont .EQ. 1 ) then
                 date(cont) = tmp3(1:14)//'00'
              else
                 date(cont) = tmp3(1:16)
              end if
!              date(cont) = year(cont)//mon(cont)//day(cont)//time(cont)
          end if
        end if
        fstr = estr+2
    end do

    RETURN

END SUBROUTINE GEMS_Share_GetL1BInfo_OMI

END MODULE Share_MOD_L1B_Read

!-------------------------------------------------------------------------------
