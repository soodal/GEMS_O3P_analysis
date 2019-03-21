!-------------------------------------------------------------------------------
!+Module to process Level2 File Name

MODULE Share_MOD_L2FileName
!-------------------------------------------------------------------------------
!+Description: 
!     Module to process Level2 File Name
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2015.02.24 Fisrt Code (YuGeun Ki, SAEAsoft) 
! 0.2     2015.08.10 Modify Product Semantic Descriptor  (HJ Lee, SAEAsoft)
!-------------------------------------------------------------------------------
      
! Declarations:

! Modules used:

!**********************************************************!
IMPLICIT NONE


!**********************************************************!
CONTAINS

!-------------------------------------------------------------------------------
!+Description: 
!      SUBROUTINE GEMS_Share_GetLv2FileName.
!
!
! Method:
!
!
! Input files:
!       mode            : File Class
!                           TEST : internal testing
!                           OGCA : on-ground calibration
!                           GSOV : ground segment overall validation, system level testing
!                           OPER : operational processing
!                           NRTI : near-real time processing
!                           OFFL : offline processing
!                           RPRO : reprocessing
!       prodName        : Level2 Final Production Name
!
! Output files:
!       Lv2FileName    : Level2 File Name
!
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2014.01.01  Fisrt Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------

SUBROUTINE GEMS_Share_GetLv2FileName(mode, prodName, Lv2FileName)

    CHARACTER(LEN=*), INTENT(IN)     :: mode
    CHARACTER(LEN=*), INTENT(IN)     :: prodName
    CHARACTER(LEN=*), INTENT(OUT)    :: Lv2FileName

    CHARACTER(LEN=3)      :: MissionID
    CHARACTER(LEN=4)      :: FileClass
    CHARACTER(LEN=10)     :: FileType
    CHARACTER(LEN=4)      :: FileCategory
    CHARACTER(LEN=6)      :: ProductSemanticDescriptor
    CHARACTER(LEN=63)     :: FileInstanceID

    CHARACTER(LEN=8)      :: yyyy_mm_dd
    CHARACTER(LEN=1)      :: T_CHAR
    CHARACTER(LEN=6)      :: HHMMSS

    CHARACTER(LEN=15)     :: SensingStartTime
    CHARACTER(LEN=15)     :: SensingEndTime
    CHARACTER(LEN=15)     :: WorkingStartTime

    CHARACTER(LEN=5)      :: absoluteOrbitNumber
    CHARACTER(LEN=2)      :: collectionNumber
    CHARACTER(LEN=6)      :: processorVersionNumber
    CHARACTER(LEN=2)      :: fileExt

    INTEGER(KIND=4)       :: idateTime(8)
    CHARACTER(LEN=8)      :: cDate
    CHARACTER(LEN=10)     :: cTime
    CHARACTER(LEN=5)      :: cZone

    MissionID    = "GEM"
    FileClass    = mode
                          ! TEST : internal testing
                          ! OGCA : on-ground calibration
                          ! GSOV : ground segment overall validation, system level testing
                          ! OPER : operational processing
                          ! NRTI : near-real time processing
                          ! OFFL : offline processing
                          ! RPRO : reprocessing

    FileCategory = "L2__"

    SELECT CASE ( prodName )
     CASE ( 'CLD' )
        ProductSemanticDescriptor = "CLOUD_"
     CASE ( 'O3P' )
        ProductSemanticDescriptor = "O3__PR"
     CASE ( 'O3T' )
        ProductSemanticDescriptor = "O3____"
     CASE ( 'NO2' )
        ProductSemanticDescriptor = "NO2___"
     CASE ( 'HCHO' )
        ProductSemanticDescriptor = "HCHO__"
     CASE ( 'SO2' )
        ProductSemanticDescriptor = "SO2___"
     CASE ( 'AOD' )
        ProductSemanticDescriptor = "AER_AI"
     CASE ( 'ALH' )
        ProductSemanticDescriptor = "AER_LH"
     CASE ( 'ALBD' )
        ProductSemanticDescriptor = "ALBD__"
    END SELECT

    FileType  = FileCategory//ProductSemanticDescriptor

    yyyy_mm_dd = "20150113"
    T_CHAR     = "T"
    HHMMSS     = "140000"
    SensingStartTime = yyyy_mm_dd//T_CHAR//HHMMSS

    yyyy_mm_dd = "20150113"
    T_CHAR     = "T"
    HHMMSS     = "140012"
    SensingEndTime = yyyy_mm_dd//T_CHAR//HHMMSS

    absoluteOrbitNumber    = "00001"
    collectionNumber       = "00"
    processorVersionNumber = "000001"     ! 010001 의 의미는 version 1.0.1

    call date_and_time(cDate, cTime, cZone, idateTime)
    yyyy_mm_dd = cDate
    T_CHAR     = "T"
    HHMMSS     = cTime

    WorkingStartTime = yyyy_mm_dd//T_CHAR//HHMMSS   

    FileInstanceID  = SensingStartTime//"_"//SensingEndTime//"_"//      &
                      absoluteOrbitNumber//"_"//collectionNumber//"_"// &
                      processorVersionNumber//"_"//WorkingStartTime 

    fileExt = "h5"

    Lv2FileName = MissionID//"_"//FileClass//"_"//FileType//"_"// &
                  FileInstanceID//"."//fileExt

    RETURN
END SUBROUTINE GEMS_Share_GetLv2FileName

!-------------------------------------------------------------------------------
!+Description: 
!      SUBROUTINE GEMS_Share_GetLv2FileName_OMI.
!
!
! Method:
!
!
! Input files:
!       mode            : File Class
!                           TEST : internal testing
!                           OGCA : on-ground calibration
!                           GSOV : ground segment overall validation, system level testing
!                           OPER : operational processing
!                           NRTI : near-real time processing
!                           OFFL : offline processing
!                           RPRO : reprocessing
!       prodName        : Level2 Final Production Name
!       granuleStTime   : L1B granule start time
!       orbitNo         : orbit number
!       prodDate        : production date
!
! Output files:
!       Lv2FileName    : Level2 File Name
!
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2016.09.09  Fisrt Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------

SUBROUTINE GEMS_Share_GetLv2FileName_OMI(mode, prodName, granuleStTime, orbitNo, prodDate, Lv2FileName)

    CHARACTER(LEN=*), INTENT(IN)     :: mode
    CHARACTER(LEN=*), INTENT(IN)     :: prodName
    CHARACTER(LEN=*), INTENT(IN)     :: granuleStTime
    CHARACTER(LEN=*), INTENT(IN)     :: orbitNo
    CHARACTER(LEN=*), INTENT(IN)     :: prodDate
    CHARACTER(LEN=*), INTENT(OUT)    :: Lv2FileName

    CHARACTER(LEN=3)      :: MissionID
    CHARACTER(LEN=4)      :: FileClass
    CHARACTER(LEN=10)     :: FileType
    CHARACTER(LEN=4)      :: FileCategory
    CHARACTER(LEN=6)      :: ProductSemanticDescriptor
    CHARACTER(LEN=63)     :: FileInstanceID

    CHARACTER(LEN=8)      :: yyyy_mm_dd
    CHARACTER(LEN=1)      :: T_CHAR
    CHARACTER(LEN=6)      :: HHMMSS

    CHARACTER(LEN=15)     :: SensingStartTime
    CHARACTER(LEN=15)     :: SensingEndTime
    CHARACTER(LEN=15)     :: WorkingStartTime

    CHARACTER(LEN=5)      :: absoluteOrbitNumber
    CHARACTER(LEN=2)      :: collectionNumber
    CHARACTER(LEN=6)      :: processorVersionNumber
    CHARACTER(LEN=2)      :: fileExt

    INTEGER(KIND=4)       :: idateTime(8)
    CHARACTER(LEN=8)      :: cDate
    CHARACTER(LEN=10)     :: cTime
    CHARACTER(LEN=5)      :: cZone

    MissionID    = "GEM"
    FileClass    = mode
                          ! TEST : internal testing
                          ! OGCA : on-ground calibration
                          ! GSOV : ground segment overall validation, system level testing
                          ! OPER : operational processing
                          ! NRTI : near-real time processing
                          ! OFFL : offline processing
                          ! RPRO : reprocessing

    FileCategory = "L2__"

    SELECT CASE ( prodName )
     CASE ( 'CLD' )
        ProductSemanticDescriptor = "CLOUD_"
     CASE ( 'O3P' )
        ProductSemanticDescriptor = "O3__PR"
     CASE ( 'O3T' )
        ProductSemanticDescriptor = "O3____"
     CASE ( 'NO2' )
        ProductSemanticDescriptor = "NO2___"
     CASE ( 'HCHO' )
        ProductSemanticDescriptor = "HCHO__"
     CASE ( 'SO2' )
        ProductSemanticDescriptor = "SO2___"
     CASE ( 'AOD' )
        ProductSemanticDescriptor = "AER_AI"
     CASE ( 'ALH' )
        ProductSemanticDescriptor = "AER_LH"
     CASE ( 'ALBD' )
        ProductSemanticDescriptor = "ALBD__"
     CASE ( 'UVidx' )
        ProductSemanticDescriptor = "UVidx_"
    END SELECT

    FileType  = FileCategory//ProductSemanticDescriptor

    yyyy_mm_dd = "20150113"
    T_CHAR     = "T"
    HHMMSS     = "140000"
    SensingStartTime = yyyy_mm_dd//T_CHAR//HHMMSS

    yyyy_mm_dd = "20150113"
    T_CHAR     = "T"
    HHMMSS     = "140012"
    SensingEndTime = yyyy_mm_dd//T_CHAR//HHMMSS

    absoluteOrbitNumber    = "00001"
    collectionNumber       = "00"
    processorVersionNumber = "000001"     ! 010001 의 의미는 version 1.0.1

    call date_and_time(cDate, cTime, cZone, idateTime)
    yyyy_mm_dd = cDate
    T_CHAR     = "T"
    HHMMSS     = cTime

    WorkingStartTime = yyyy_mm_dd//T_CHAR//HHMMSS   

    FileInstanceID  = granuleStTime//"_"//"o"//orbitNo//"_"//       &
                      prodDate//"_"//collectionNumber//"_"//        &
                      processorVersionNumber//"_"//WorkingStartTime 

    fileExt = "h5"

    Lv2FileName = MissionID//"_"//FileClass//"_"//FileType//"_"// &
                  FileInstanceID//"."//fileExt

    RETURN
END SUBROUTINE GEMS_Share_GetLv2FileName_OMI

END MODULE Share_MOD_L2FileName

!-------------------------------------------------------------------------------
