!-------------------------------------------------------------------------------
!+Module to declare and initialize global variables and parameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE Share_MOD_Constants_Variables_Synt 

!-------------------------------------------------------------------------------
!+Description: 
!     Defines GEMS global constants and variables.
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2013.10.31 Fisrt Code (Sung-Hyun Park, Seasoft) 
! 0.2     2014.08.27 remove Cloud Mask Flag  
!                    remove Pixel Availity Flag 
!                    add new variables(Y.K. Ki, Seasoft) 
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Modules:
      ! USE Share_MOD_Types
!      USE Share_MOD_P_Constants_Variables

      IMPLICIT NONE
      SAVE
 
!-------------------------------------------------------------------------------
! Global Constants:

  ! Calculation Constants
      REAL              :: gcf_pi_c               ! pi
      REAL              :: gcf_pi2_c              ! pi*2
      REAL              :: gcf_deg2rad_c          ! Degree to Radian
      REAL              :: gcf_rad2deg_c          ! Radian to Degree

      REAL              :: gcf_planck_c           ! Planck's constant  [J*s]  
      REAL              :: gcf_light_c            ! Speed of light [m/s]
      REAL              :: gcf_boltz_c            ! Boltzmann constant [J/K]

      REAL              :: gcf_solar_c            ! Solar constant     [W/m^2]
      REAL              :: gcf_press_c            ! standard pressure  [hPa]  
      REAL              :: gcf_kelvin_c           ! Absolute Temperature[K]
      
  ! Land/Sea Flag        
      INTEGER           :: gci_flg_ls_sea         ! Land/Sea Flag; Sea
      INTEGER           :: gci_flg_ls_land        ! Land/Sea Flag; Land
      INTEGER           :: gci_flg_ls_coast       ! Land/Sea Flag; Coast
      INTEGER           :: gci_flg_ls_unavail     ! Land/Sea Flag; Unavailable
      INTEGER           :: gci_flg_ls_exosphere   ! Land/Sea Flag; Exospher
      
  ! Cloud Mask Flag

      
  ! Pixel Availity Flag  

!-------------------------------------------------------------------------------
! Global Variables:

  ! Type of Angles
      TYPE angles_type          
        REAL, ALLOCATABLE   :: vza(:,:)           ! Viewing zenith angle
        REAL, ALLOCATABLE   :: sza(:,:)           ! Solar zenith angle
        REAL, ALLOCATABLE   :: raa(:,:)           ! Relative azimuth angle (sun/satellite)
      END TYPE angles_type                             

      TYPE(angles_type)     :: gvf_angles         ! Angles of vza, sza, raa
          
  ! Type of Latitude/Longitude    
      TYPE latlon_type
         REAL, ALLOCATABLE  :: lat(:,:)           ! Latitude
         REAL, ALLOCATABLE  :: lon(:,:)           ! Longitude
      END TYPE latlon_type
      
      TYPE(latlon_type)     :: gvf_latlon         ! Latitude/Longitude 
      
  ! Variables of Spectrum    
      REAL(KIND=8),     POINTER :: gvf_euv1_rad   (:,:,:) ! Earth UV1 Radiance
      REAL(KIND=8),     POINTER :: gvf_euv2_rad   (:,:,:) ! Earth UV2 Radiance
      REAL(KIND=8),     POINTER :: gvf_evis_rad   (:,:,:) ! Earth VIS Radiance
      REAL(KIND=8),     POINTER :: gvf_suv2_irrad (:,:  ) ! Sun UV2 Irradiance
      REAL(KIND=8),     POINTER :: gvf_suv2_prec  (:,:  ) ! Sun UV2 IrradiancePrecision
      REAL(KIND=8),     POINTER :: gvf_euv1_prec  (:,:,:) ! Earth UV1 Radiance Precision
      REAL(KIND=8),     POINTER :: gvf_euv2_prec  (:,:,:) ! Earth UV2 Radiance Precision
      REAL(KIND=8),     POINTER :: gvf_svis_irrad (:,:  ) ! Sun Visible Irradiance
      REAL(KIND=8),     ALLOCATABLE :: gvf_albedo     (:,:,:) ! Surface albedo
      
  ! Variables of Others     
      REAL(KIND=8),     ALLOCATABLE :: gvf_elevation  (:,:)   ! Surface elevation
      REAL(KIND=8),     ALLOCATABLE :: gvf_amfgeo     (:,:)   ! geometric air mass factor
  
  ! Variables of Flag    
      INTEGER,    ALLOCATABLE :: gvi_flg_ls     (:,:)   ! Land/Sea Flag
      INTEGER,    ALLOCATABLE :: gvi_flg_cm     (:,:)   ! Cloud Mask Flag
      INTEGER,    ALLOCATABLE :: gvi_flg_px     (:,:)   ! Pixel Availity Flag

  ! Reflectance
      REAL(KIND=8),     ALLOCATABLE :: gvf_rflctnc            ! 관측된 대기상한 반사도

  ! shift
      REAL(KIND=4),     ALLOCATABLE :: gvf_shft               ! 복사조도와 복사휘도 calibration에 사용되는  shift factor

  ! squaze
      REAL(KIND=4),     ALLOCATABLE :: gvf_sqz                ! 복사조도와 복사휘도 calibration에 사용되는  squaze factor

  ! Measurement time
      CHARACTER(LEN=15),     ALLOCATABLE :: gvc_msr_time      ! 관측 년/월/일/시간

  ! Normalized  radiance
      REAL(KIND=4),     ALLOCATABLE :: gvf_nrad_obs           ! 위성에서 관측된  복사휘도/복사조도 일별 관측값
      REAL(KIND=4),     ALLOCATABLE :: gvf_nrad_avg           ! 위성에서 관측된  복사휘도/복사조도 1년 평균값

  ! Wavelength
      REAL(KIND=4),     ALLOCATABLE :: gvf_wav                ! 복사조도, 복사휘도에 해당되는  파장
      REAL(KIND=4),     POINTER :: gvf_euv1_radwave(:,:,:)    ! Earth UV1 Wavelength          
      REAL(KIND=4),     POINTER :: gvf_euv2_radwave(:,:,:)    ! Earth UV2 Wavelength          
      REAL(KIND=4),     POINTER :: gvf_evis_radwave(:,:,:)    ! Earth VIS Wavelength          
      REAL(KIND=4),     POINTER :: gvf_suv2_irrwave(:,:  )    ! Sun UV2 Wavelength          
      REAL(KIND=4),     POINTER :: gvf_svis_irrwave(:,:  )    ! Sun Visible Irradiance      

  ! terrain pressure
      REAL(KIND=4),     ALLOCATABLE :: gvf_pt                 ! 지형 압력

  ! cloud pressure
      REAL(KIND=8),     ALLOCATABLE :: gvf_pc                 ! 구름 압력

  ! array size of global allocatable variables
      TYPE sz_kind          
        INTEGER(KIND=4)     :: nx
        INTEGER(KIND=4)     :: ny
        INTEGER(KIND=4)     :: nz
      END TYPE sz_kind                             

      TYPE glb_alloc_sz_type          
        TYPE(sz_kind)       :: angles       ! memory size of Angles of vza, sza, raa
        TYPE(sz_kind)       :: latlon       ! memory size of Latitude/Longitude 
        TYPE(sz_kind)       :: euv1_rad     ! memory size of Earth UV1 Radiance       
        TYPE(sz_kind)       :: euv1_radwav  ! memory size of Earth UV1 Wavelength       
        TYPE(sz_kind)       :: euv1_prec    ! memory size of Earth UV1 RadiancePrecision
        TYPE(sz_kind)       :: euv2_rad     ! memory size of Earth UV2 Radiance       
        TYPE(sz_kind)       :: euv2_radwav  ! memory size of Earth UV2 Wavelength    
        TYPE(sz_kind)       :: euv2_prec    ! memory size of Earth UV2 RadiancePrecision   
        TYPE(sz_kind)       :: evis_rad     ! memory size of Earth VIS Radiance       
        TYPE(sz_kind)       :: evis_radwav  ! memory size of Earth VIS Wavelength       
        TYPE(sz_kind)       :: suv1_irrad   ! memory size of SUN UV1 Irradiance     
        TYPE(sz_kind)       :: suv1_irrwav  ! memory size of SUN UV1 Wavelength     
        TYPE(sz_kind)       :: suv2_irrad   ! memory size of SUN UV2 Irradiance     
        TYPE(sz_kind)       :: suv2_prec    ! memory size of SUN UV2 IrradiancePrecision     
        TYPE(sz_kind)       :: suv2_irrwav  ! memory size of SUN UV2 Wavelength     
        TYPE(sz_kind)       :: svis_irrad   ! memory size of SUN VIS Irradiance     
        TYPE(sz_kind)       :: svis_irrwav  ! memory size of SUN VIS Wavelength     
        TYPE(sz_kind)       :: albedo       ! memory size of Surface albedo 
        TYPE(sz_kind)       :: elevation    ! memory size of Surface elevation        
        TYPE(sz_kind)       :: amfgeo       ! memory size of geometric air mass factor
        TYPE(sz_kind)       :: flg_ls       ! memory size of Land/Sea Flag      
        TYPE(sz_kind)       :: flg_cm       ! memory size of Cloud Mask Flag    
        TYPE(sz_kind)       :: flg_px       ! memory size of Pixel Availity Flag

	TYPE(sz_kind)       :: euv1         ! memory size of Earth UV1       
	TYPE(sz_kind)       :: euv2         ! memory size of Earth UV2       
	TYPE(sz_kind)       :: evis         ! memory size of Earth VIS       
	TYPE(sz_kind)       :: suv1         ! memory size of SUN   UV1       
	TYPE(sz_kind)       :: suv2         ! memory size of SUN   UV2       
	TYPE(sz_kind)       :: svis         ! memory size of SUN   VIS       

      END TYPE glb_alloc_sz_type                             

      TYPE(glb_alloc_sz_type)   :: gvs_size

  ! array size being used on L1B Dataset 
      INTEGER(KIND=4)       :: gci_nline      
      INTEGER(KIND=4)       :: gci_nline_rug      
      INTEGER(KIND=4)       :: gci_nline_rvg      
      INTEGER(KIND=4)       :: gci_nline_irr    
      INTEGER(KIND=4)       :: gci_nwavel     
      INTEGER(KIND=4)       :: gci_nwavel2    
      INTEGER(KIND=4)       :: gci_nwavel3    
      INTEGER(KIND=4)       :: gci_nxtrack    
      INTEGER(KIND=4)       :: gci_nxtrack2   
      INTEGER(KIND=4)       :: gci_nwavelcoef 
      INTEGER(KIND=4)       :: gci_ntimes      

  ! array size being used on L2 Dataset
      INTEGER(KIND=4)       :: gci_nline_l2               
      INTEGER(KIND=4)       :: gci_nline_rug_l2               
      INTEGER(KIND=4)       :: gci_nline_rvg_l2               
      INTEGER(KIND=4)       :: gci_nline_irr_l2               
      INTEGER(KIND=4)       :: gci_nline2_l2              
      INTEGER(KIND=4)       :: gci_nwavel_l2              
      INTEGER(KIND=4)       :: gci_nwavel2_l2             
      INTEGER(KIND=4)       :: gci_nwavel3_l2             
      INTEGER(KIND=4)       :: gci_nwavel4_l2             
      INTEGER(KIND=4)       :: gci_nxtrack_l2             
      INTEGER(KIND=4)       :: gci_nxtrack2_l2            
      INTEGER(KIND=4)       :: gci_nlayer                 
      INTEGER(KIND=4)       :: gci_nlayer2               
      INTEGER(KIND=4)       :: gci_nlayerp
      INTEGER(KIND=4)       :: gci_nmon                  
      INTEGER(KIND=4)       :: gci_nlat                  
      INTEGER(KIND=4)       :: gci_nlon                  
      INTEGER(KIND=4)       :: gci_ncomp
      INTEGER(KIND=4)       :: gci_nResiduals
      INTEGER(KIND=4)       :: gci_nColumns 

CONTAINS

!-------------------------------------------------------------------------------
!+Description: 
!      전역상수값을 초기화
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
! 0.1     2015. .   Fisrt Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------
SUBROUTINE GEMS_Share_Init_GlobalConstants_Synt
    USE NumTimes
    IMPLICIT NONE
    INTEGER :: nline
    INTEGER :: rc

    INTEGER  GEMS_Share_Hdf5InitL1BStorageSize
    INTEGER  GEMS_Share_Hdf5InitL2StorageSize
    external GEMS_Share_Hdf5InitL1BStorageSize
    external GEMS_Share_Hdf5InitL2StorageSize

  ! Calculation Constants
    gcf_pi_c                = 3.141592        ! pi
    gcf_pi2_c               = gcf_pi_c*2.     ! pi*2
    gcf_deg2rad_c           = gcf_pi_c/180.   ! Degree to Radian
    gcf_rad2deg_c           = 180./gcf_pi_c   ! Radian to Degree

    gcf_planck_c            = 6.62606876E-34  ! Planck's constant  [J*s]  
    gcf_light_c             = 2.99793E8       ! Speed of light [m/s]
    gcf_boltz_c             = 1.380622E-23    ! Boltzmann constant [J/K]

    gcf_solar_c             = 1367.0          ! Solar constant     [W/m^2]
    gcf_press_c             = 1013.25         ! standard pressure  [hPa]  
    gcf_kelvin_c            = 273.16          ! Absolute Temperature[K]
      
  ! Land/Sea Flag     
    gci_flg_ls_sea          = 0               ! Land/Sea Flag; Sea
    gci_flg_ls_land         = 1               ! Land/Sea Flag; Land
    gci_flg_ls_coast        = 2               ! Land/Sea Flag; Coast
    gci_flg_ls_unavail      = -999            ! Land/Sea Flag; Unavailable
    gci_flg_ls_exosphere    = -999            ! Land/Sea Flag; Exospher


  ! Size for L1B dataset
    !CALL GEMS_NumTimes_Read(nline)    ! geun closed
    gci_nline               = nline
    gci_nwavel              = 159
    gci_nwavel2             = 557
    gci_nwavel3             = 751
    gci_nxtrack             = 30
    gci_nxtrack2            = 60
    gci_nwavelcoef          = 5
    gci_ntimes              = gci_nline

 ! Size for L2 dataset
    !gci_nline_l2               = nline     !변경  1644에서 nline으로 
    gci_nline_l2               = 2048   !geun
    gci_nline2_l2              = 330
    gci_nwavel_l2              = 3
    gci_nwavel2_l2             = 12
    gci_nwavel3_l2             = 20
    gci_nwavel4_l2             = 23
    gci_nxtrack_l2             = 60
    gci_nxtrack2_l2            = 10   ! geun
    gci_nlayer                 = 24
    gci_nlayer2                = 11
!    gci_nlayerp = 777
    gci_nlayerp                = gci_nlayer+1
    gci_nmon                   = 12
    gci_nlat                   = 360
    gci_nlon                   = 720
    gci_ncomp                  = 2
    gci_nResiduals             = 20
    gci_nColumns               = 3

    !
    !--- Size Declaration for Global Allocatable Variables
    !
    gvs_size%angles%nx          =   gci_nxtrack2     ! x coordi memory size of Angles of vza, sza, raa
    gvs_size%angles%ny          =   gci_ntimes       ! y coordi memory size of Angles of vza, sza, raa
    gvs_size%latlon%nx          =   gci_nxtrack2     ! x coordi memory size of Latitude/Longitude 
    gvs_size%latlon%ny          =   gci_ntimes       ! y coordi memory size of Latitude/Longitude 

    !
    !--- Dimension Size for L1B datasets
    !
    gvs_size%euv1%nx            =   gci_nwavel       ! x coordi memory size of EUV1       
    gvs_size%euv1%ny            =   gci_nxtrack      ! y coordi memory size of EUV1       
    gvs_size%euv1%nz            =   gci_ntimes       ! z coordi memory size of EUV1       
			       
    gvs_size%euv2%nx            =   gci_nwavel2      ! x coordi memory size of EUV2       
    gvs_size%euv2%ny            =   gci_nxtrack2     ! y coordi memory size of EUV2       
    gvs_size%euv2%nz            =   gci_ntimes       ! z coordi memory size of EUV2       
			       
    gvs_size%evis%nx            =   gci_nwavel3      ! x coordi memory size of EVIS     
    gvs_size%evis%ny            =   gci_nxtrack2     ! y coordi memory size of EVIS     
    gvs_size%evis%nz            =   gci_ntimes       ! z coordi memory size of EVIS     

    gvs_size%suv1%nx            =   gci_nwavel       ! x coordi memory size of SUV1     
    gvs_size%suv1%ny            =   gci_nxtrack      ! y coordi memory size of SUV1     
    gvs_size%suv1%nz            =   gci_ntimes       ! z coordi memory size of SUV1     
			       
    gvs_size%suv2%nx            =   gci_nwavel2      ! x coordi memory size of SUV2     
    gvs_size%suv2%ny            =   gci_nxtrack2     ! y coordi memory size of SUV2     
    gvs_size%suv2%nz            =   gci_ntimes       ! z coordi memory size of SUV2     
			       
    gvs_size%svis%nx            =   gci_nwavel3      ! x coordi memory size of SVIS     
    gvs_size%svis%ny            =   gci_nxtrack2     ! y coordi memory size of SVIS     
    gvs_size%svis%nz            =   gci_ntimes       ! z coordi memory size of SVIS     
    !---

    !
    !--- Dimension size  of radiance, wavelength, irradiance based on OMI
    !                   nx          ny          nz
    !                  ---         ---        ----
    ! UV1_rad:         159          30        1644    gci_nwavel    gci_nxtrack     gci_ntimes
    ! wave:            159          30        1644
    !--
    ! UV2_rad:         557          60        1644    gci_nwavel2   gci_nxtrack2    gci_ntimes
    ! wave:            557          60        1644
    !--
    ! VIS_rad:         751          60        1644    gci_nwavel3   gci_nxtrack2    gci_ntimes
    ! wave:            751          60        1644
    !--
    ! UV1_sun_rad:     159          30                gci_nwavel    gci_nxtrack
    ! wave             159          30
    !--
    ! UV2_sun_rad:     557          60                gci_nwavel2   gci_nxtrack2   
    ! wave             557          60
    ! prec             557          60
    !--
    ! VIS_sun_rad      751          60                gci_nwavel3   gci_nxtrack2  
    ! wave             751          60
    !
    !------------------------------------------------------------------------------------------

    !
    !--- Dimension Size for Radiance, Wavelength, Irradiance datasets
    !
    gvs_size%euv1_rad%nx        =   gvs_size%euv1%nx ! x coordi memory size of EUV1 Radiance       
    gvs_size%euv1_rad%ny        =   gvs_size%euv1%ny ! y coordi memory size of EUV1 Radiance       
    gvs_size%euv1_rad%nz        =   gvs_size%euv1%nz ! z coordi memory size of EUV1 Radiance       
    gvs_size%euv1_radwav%nx     =   gvs_size%euv1%nx ! x coordi memory size of EUV1 Wavelength       
    gvs_size%euv1_radwav%ny     =   gvs_size%euv1%ny ! y coordi memory size of EUV1 Wavelength       
    gvs_size%euv1_radwav%nz     =   gvs_size%euv1%nz ! z coordi memory size of EUV1 Wavelength       
    gvs_size%euv1_prec%nx       =   gvs_size%euv1%nx ! x coordi memory size of EUV1 RadiancePrecision     
    gvs_size%euv1_prec%ny       =   gvs_size%euv1%ny ! y coordi memory size of EUV1 RadiancePrecision     
    gvs_size%euv1_prec%nz       =   gvs_size%euv1%nz ! z coordi memory size of EUV1 RadiancePrecision

    gvs_size%euv2_rad%nx        =   gvs_size%euv2%nx ! x coordi memory size of EUV2 Radiance       
    gvs_size%euv2_rad%ny        =   gvs_size%euv2%ny ! y coordi memory size of EUV2 Radiance       
    gvs_size%euv2_rad%nz        =   gvs_size%euv2%nz ! z coordi memory size of EUV2 Radiance       
    gvs_size%euv2_radwav%nx     =   gvs_size%euv2%nx ! x coordi memory size of EUV2 Wavelength       
    gvs_size%euv2_radwav%ny     =   gvs_size%euv2%ny ! y coordi memory size of EUV2 Wavelength       
    gvs_size%euv2_radwav%nz     =   gvs_size%euv2%nz ! z coordi memory size of EUV2 Wavelength       
    gvs_size%euv2_prec%nx       =   gvs_size%euv2%nx ! x coordi memory size of EUV2 RadiancePrecision     
    gvs_size%euv2_prec%ny       =   gvs_size%euv2%ny ! y coordi memory size of EUV2 RadiancePrecision     
    gvs_size%euv2_prec%nz       =   gvs_size%euv2%nz ! z coordi memory size of EUV2 RadiancePrecision


    gvs_size%evis_rad%nx        =   gvs_size%evis%nx ! x coordi memory size of EVIS Radiance       
    gvs_size%evis_rad%ny        =   gvs_size%evis%ny ! y coordi memory size of EVIS Radiance       
    gvs_size%evis_rad%nz        =   gvs_size%evis%nz ! z coordi memory size of EVIS Radiance       
    gvs_size%evis_radwav%nx     =   gvs_size%evis%nx ! x coordi memory size of EVIS Wavelength       
    gvs_size%evis_radwav%ny     =   gvs_size%evis%ny ! y coordi memory size of EVIS Wavelength       
    gvs_size%evis_radwav%nz     =   gvs_size%evis%nz ! z coordi memory size of EVIS Wavelength       

    gvs_size%suv1_irrad%nx      =   gvs_size%suv1%nx ! x coordi memory size of SUV1 Irradiance     
    gvs_size%suv1_irrad%ny      =   gvs_size%suv1%ny ! y coordi memory size of SUV1 Irradiance     
    gvs_size%suv1_irrad%nz      =   gvs_size%suv1%nz ! z coordi memory size of SUV1 Irradiance     
    gvs_size%suv1_irrwav%nx     =   gvs_size%suv1%nx ! x coordi memory size of SUV1 Wavelength     
    gvs_size%suv1_irrwav%ny     =   gvs_size%suv1%ny ! y coordi memory size of SUV1 Wavelength     
    gvs_size%suv1_irrwav%nz     =   gvs_size%suv1%nz ! z coordi memory size of SUV1 Wavelength     

    gvs_size%suv2_irrad%nx      =   gvs_size%suv2%nx ! x coordi memory size of SUV2 Irradiance     
    gvs_size%suv2_irrad%ny      =   gvs_size%suv2%ny ! y coordi memory size of SUV2 Irradiance     
    gvs_size%suv2_irrad%nz      =   gvs_size%suv2%nz ! z coordi memory size of SUV2 Irradiance     
    gvs_size%suv2_prec%nx       =   gvs_size%suv2%nx ! x coordi memory size of SUV2 IrradiancePrecision     
    gvs_size%suv2_prec%ny       =   gvs_size%suv2%ny ! y coordi memory size of SUV2 IrradiancePrecision     
    gvs_size%suv2_prec%nz       =   gvs_size%suv2%nz ! z coordi memory size of SUV2 IrradiancePrecision     
    gvs_size%suv2_irrwav%nx     =   gvs_size%suv2%nx ! x coordi memory size of SUV2 Wavelength     
    gvs_size%suv2_irrwav%ny     =   gvs_size%suv2%ny ! y coordi memory size of SUV2 Wavelength     
    gvs_size%suv2_irrwav%nz     =   gvs_size%suv2%nz ! z coordi memory size of SUV2 Wavelength     

    gvs_size%svis_irrad%nx      =   gvs_size%svis%nx ! x coordi memory size of SVIS Irradiance     
    gvs_size%svis_irrad%ny      =   gvs_size%svis%ny ! y coordi memory size of SVIS Irradiance     
    gvs_size%svis_irrad%nz      =   gvs_size%svis%nz ! z coordi memory size of SVIS Irradiance     
    gvs_size%svis_irrwav%nx     =   gvs_size%svis%nx ! x coordi memory size of SVIS Wavelength     
    gvs_size%svis_irrwav%ny     =   gvs_size%svis%ny ! y coordi memory size of SVIS Wavelength     
    gvs_size%svis_irrwav%nz     =   gvs_size%svis%nz ! z coordi memory size of SVIS Wavelength     
    !---

    gvs_size%albedo%nx          =   gci_nxtrack2     ! x coordi memory size of Surface albedo 
    gvs_size%albedo%ny          =   gci_ntimes       ! y coordi memory size of Surface albedo 
    gvs_size%albedo%nz          =   gci_ntimes       ! z coordi memory size of Surface albedo 
    gvs_size%elevation%nx       =   gci_nxtrack2     ! x coordi memory size of Surface elevation        
    gvs_size%elevation%ny       =   gci_ntimes       ! y coordi memory size of Surface elevation        
    gvs_size%amfgeo%nx          =   gci_nxtrack2     ! x coordi memory size of geometric air mass factor
    gvs_size%amfgeo%ny          =   gci_ntimes       ! y coordi memory size of geometric air mass factor
    gvs_size%flg_ls%nx          =   gci_nxtrack2     ! x coordi memory size of Land/Sea Flag      
    gvs_size%flg_ls%ny          =   gci_ntimes       ! y coordi memory size of Land/Sea Flag      
    gvs_size%flg_cm%nx          =   gci_nxtrack2     ! x coordi memory size of Cloud Mask Flag    
    gvs_size%flg_cm%ny          =   gci_ntimes       ! y coordi memory size of Cloud Mask Flag    
    gvs_size%flg_px%nx          =   gci_nxtrack2     ! x coordi memory size of Pixel Availity Flag
    gvs_size%flg_px%ny          =   gci_ntimes       ! y coordi memory size of Pixel Availity Flag

    !
    ! initialize storage size of L1B File for C-Interface
    !
    rc = GEMS_Share_Hdf5InitL1BStorageSize(            &
      gci_nline            ,                           &
      gci_nwavel           ,                           &
      gci_nwavel2          ,                           &
      gci_nwavel3          ,                           &
      gci_nxtrack          ,                           &
      gci_nxtrack2         ,                           &
      gci_nwavelcoef       ,                           &
      gci_ntimes                                       &
    )

    !
    ! initialize storage size of L2 File for C-Interface
    !
    rc = GEMS_Share_Hdf5InitL2StorageSize(             &
              gci_nline_l2       ,                     &
              gci_nline2_l2      ,                     &
              gci_nwavel_l2      ,                     &
              gci_nwavel2_l2     ,                     &
              gci_nwavel3_l2     ,                     &
              gci_nwavel4_l2     ,                     &
              gci_nxtrack_l2     ,                     &
              gci_nxtrack2_l2    ,                     &
              gci_nlayer         ,                     &
              gci_nlayer2        ,                     &
              gci_nlayerp        ,                     &
              gci_nmon           ,                     &
              gci_nlat           ,                     &
              gci_nlon           ,                     &
              gci_ncomp          ,                     &
              gci_nResiduals     ,                     &
              gci_nColumns       ,                     &
              6                                        &   ! Time Column
    )

END SUBROUTINE GEMS_Share_Init_GlobalConstants_Synt

SUBROUTINE GEMS_Share_Init_GlobalConstants2_Synt( ctl_fpath )
    USE NumTimes
    IMPLICIT NONE
    INTEGER :: nline
    INTEGER :: rc
    CHARACTER(LEN=128),INTENT(IN) :: ctl_fpath

    INTEGER  GEMS_Share_Hdf5InitL1BStorageSize
    INTEGER  GEMS_Share_Hdf5InitL2StorageSize
    external GEMS_Share_Hdf5InitL1BStorageSize
    external GEMS_Share_Hdf5InitL2StorageSize

  ! Calculation Constants
    gcf_pi_c                = 3.141592        ! pi
    gcf_pi2_c               = gcf_pi_c*2.     ! pi*2
    gcf_deg2rad_c           = gcf_pi_c/180.   ! Degree to Radian
    gcf_rad2deg_c           = 180./gcf_pi_c   ! Radian to Degree

    gcf_planck_c            = 6.62606876E-34  ! Planck's constant  [J*s]  
    gcf_light_c             = 2.99793E8       ! Speed of light [m/s]
    gcf_boltz_c             = 1.380622E-23    ! Boltzmann constant [J/K]

    gcf_solar_c             = 1367.0          ! Solar constant     [W/m^2]
    gcf_press_c             = 1013.25         ! standard pressure  [hPa]  
    gcf_kelvin_c            = 273.16          ! Absolute Temperature[K]
      
  ! Land/Sea Flag     
    gci_flg_ls_sea          = 0               ! Land/Sea Flag; Sea
    gci_flg_ls_land         = 1               ! Land/Sea Flag; Land
    gci_flg_ls_coast        = 2               ! Land/Sea Flag; Coast
    gci_flg_ls_unavail      = -999            ! Land/Sea Flag; Unavailable
    gci_flg_ls_exosphere    = -999            ! Land/Sea Flag; Exospher


  ! Size for L1B dataset
    CALL GEMS_NumTimes_Read2(nline, ctl_fpath)
    gci_nline               = nline
    gci_nwavel              = 159
    gci_nwavel2             = 557
    gci_nwavel3             = 751
    gci_nxtrack             = 30
    gci_nxtrack2            = 60
    gci_nwavelcoef          = 5
    gci_ntimes              = gci_nline

 ! Size for L2 dataset
    gci_nline_l2               = nline     !변경  1644에서 nline으로 
    gci_nline2_l2              = 330
    gci_nwavel_l2              = 3
    gci_nwavel2_l2             = 12
    gci_nwavel3_l2             = 20
    gci_nwavel4_l2             = 23
    gci_nxtrack_l2             = 60
    gci_nxtrack2_l2            = 30
    gci_nlayer                 = 24
    gci_nlayer2                = 11
!    gci_nlayerp = 777
    gci_nlayerp                = gci_nlayer+1
    gci_nmon                   = 12
    gci_nlat                   = 360
    gci_nlon                   = 720
    gci_ncomp                  = 2
    !
    !--- Size Declaration for Global Allocatable Variables
    !
    gvs_size%angles%nx          =   gci_nxtrack2     ! x coordi memory size of Angles of vza, sza, raa
    gvs_size%angles%ny          =   gci_ntimes       ! y coordi memory size of Angles of vza, sza, raa
    gvs_size%latlon%nx          =   gci_nxtrack2     ! x coordi memory size of Latitude/Longitude 
    gvs_size%latlon%ny          =   gci_ntimes       ! y coordi memory size of Latitude/Longitude 

    !
    !--- Dimension Size for L1B datasets
    !
    gvs_size%euv1%nx            =   gci_nwavel       ! x coordi memory size of EUV1       
    gvs_size%euv1%ny            =   gci_nxtrack      ! y coordi memory size of EUV1       
    gvs_size%euv1%nz            =   gci_ntimes       ! z coordi memory size of EUV1       
			       
    gvs_size%euv2%nx            =   gci_nwavel2      ! x coordi memory size of EUV2       
    gvs_size%euv2%ny            =   gci_nxtrack2     ! y coordi memory size of EUV2       
    gvs_size%euv2%nz            =   gci_ntimes       ! z coordi memory size of EUV2       
			       
    gvs_size%evis%nx            =   gci_nwavel3      ! x coordi memory size of EVIS     
    gvs_size%evis%ny            =   gci_nxtrack2     ! y coordi memory size of EVIS     
    gvs_size%evis%nz            =   gci_ntimes       ! z coordi memory size of EVIS     

    gvs_size%suv1%nx            =   gci_nwavel       ! x coordi memory size of SUV1     
    gvs_size%suv1%ny            =   gci_nxtrack      ! y coordi memory size of SUV1     
    gvs_size%suv1%nz            =   gci_ntimes       ! z coordi memory size of SUV1     
			       
    gvs_size%suv2%nx            =   gci_nwavel2      ! x coordi memory size of SUV2     
    gvs_size%suv2%ny            =   gci_nxtrack2     ! y coordi memory size of SUV2     
    gvs_size%suv2%nz            =   gci_ntimes       ! z coordi memory size of SUV2     
			       
    gvs_size%svis%nx            =   gci_nwavel3      ! x coordi memory size of SVIS     
    gvs_size%svis%ny            =   gci_nxtrack2     ! y coordi memory size of SVIS     
    gvs_size%svis%nz            =   gci_ntimes       ! z coordi memory size of SVIS     
    !---

    !
    !--- Dimension size  of radiance, wavelength, irradiance based on OMI
    !                   nx          ny          nz
    !                  ---         ---        ----
    ! UV1_rad:         159          30        1644    gci_nwavel    gci_nxtrack     gci_ntimes
    ! wave:            159          30        1644
    !--
    ! UV2_rad:         557          60        1644    gci_nwavel2   gci_nxtrack2    gci_ntimes
    ! wave:            557          60        1644
    !--
    ! VIS_rad:         751          60        1644    gci_nwavel3   gci_nxtrack2    gci_ntimes
    ! wave:            751          60        1644
    !--
    ! UV1_sun_rad:     159          30                gci_nwavel    gci_nxtrack
    ! wave             159          30
    !--
    ! UV2_sun_rad:     557          60                gci_nwavel2   gci_nxtrack2   
    ! wave             557          60
    ! prec             557          60
    !--
    ! VIS_sun_rad      751          60                gci_nwavel3   gci_nxtrack2  
    ! wave             751          60
    !
    !------------------------------------------------------------------------------------------

    !
    !--- Dimension Size for Radiance, Wavelength, Irradiance datasets
    !
    gvs_size%euv1_rad%nx        =   gvs_size%euv1%nx ! x coordi memory size of EUV1 Radiance       
    gvs_size%euv1_rad%ny        =   gvs_size%euv1%ny ! y coordi memory size of EUV1 Radiance       
    gvs_size%euv1_rad%nz        =   gvs_size%euv1%nz ! z coordi memory size of EUV1 Radiance       
    gvs_size%euv1_radwav%nx     =   gvs_size%euv1%nx ! x coordi memory size of EUV1 Wavelength       
    gvs_size%euv1_radwav%ny     =   gvs_size%euv1%ny ! y coordi memory size of EUV1 Wavelength       
    gvs_size%euv1_radwav%nz     =   gvs_size%euv1%nz ! z coordi memory size of EUV1 Wavelength       
    gvs_size%euv1_prec%nx       =   gvs_size%euv1%nx ! x coordi memory size of EUV1 RadiancePrecision     
    gvs_size%euv1_prec%ny       =   gvs_size%euv1%ny ! y coordi memory size of EUV1 RadiancePrecision     
    gvs_size%euv1_prec%nz       =   gvs_size%euv1%nz ! z coordi memory size of EUV1 RadiancePrecision

    gvs_size%euv2_rad%nx        =   gvs_size%euv2%nx ! x coordi memory size of EUV2 Radiance       
    gvs_size%euv2_rad%ny        =   gvs_size%euv2%ny ! y coordi memory size of EUV2 Radiance       
    gvs_size%euv2_rad%nz        =   gvs_size%euv2%nz ! z coordi memory size of EUV2 Radiance       
    gvs_size%euv2_radwav%nx     =   gvs_size%euv2%nx ! x coordi memory size of EUV2 Wavelength       
    gvs_size%euv2_radwav%ny     =   gvs_size%euv2%ny ! y coordi memory size of EUV2 Wavelength       
    gvs_size%euv2_radwav%nz     =   gvs_size%euv2%nz ! z coordi memory size of EUV2 Wavelength       
    gvs_size%euv2_prec%nx       =   gvs_size%euv2%nx ! x coordi memory size of EUV2 RadiancePrecision     
    gvs_size%euv2_prec%ny       =   gvs_size%euv2%ny ! y coordi memory size of EUV2 RadiancePrecision     
    gvs_size%euv2_prec%nz       =   gvs_size%euv2%nz ! z coordi memory size of EUV2 RadiancePrecision


    gvs_size%evis_rad%nx        =   gvs_size%evis%nx ! x coordi memory size of EVIS Radiance       
    gvs_size%evis_rad%ny        =   gvs_size%evis%ny ! y coordi memory size of EVIS Radiance       
    gvs_size%evis_rad%nz        =   gvs_size%evis%nz ! z coordi memory size of EVIS Radiance       
    gvs_size%evis_radwav%nx     =   gvs_size%evis%nx ! x coordi memory size of EVIS Wavelength       
    gvs_size%evis_radwav%ny     =   gvs_size%evis%ny ! y coordi memory size of EVIS Wavelength       
    gvs_size%evis_radwav%nz     =   gvs_size%evis%nz ! z coordi memory size of EVIS Wavelength       

    gvs_size%suv1_irrad%nx      =   gvs_size%suv1%nx ! x coordi memory size of SUV1 Irradiance     
    gvs_size%suv1_irrad%ny      =   gvs_size%suv1%ny ! y coordi memory size of SUV1 Irradiance     
    gvs_size%suv1_irrad%nz      =   gvs_size%suv1%nz ! z coordi memory size of SUV1 Irradiance     
    gvs_size%suv1_irrwav%nx     =   gvs_size%suv1%nx ! x coordi memory size of SUV1 Wavelength     
    gvs_size%suv1_irrwav%ny     =   gvs_size%suv1%ny ! y coordi memory size of SUV1 Wavelength     
    gvs_size%suv1_irrwav%nz     =   gvs_size%suv1%nz ! z coordi memory size of SUV1 Wavelength     

    gvs_size%suv2_irrad%nx      =   gvs_size%suv2%nx ! x coordi memory size of SUV2 Irradiance     
    gvs_size%suv2_irrad%ny      =   gvs_size%suv2%ny ! y coordi memory size of SUV2 Irradiance     
    gvs_size%suv2_irrad%nz      =   gvs_size%suv2%nz ! z coordi memory size of SUV2 Irradiance     
    gvs_size%suv2_prec%nx       =   gvs_size%suv2%nx ! x coordi memory size of SUV2 IrradiancePrecision     
    gvs_size%suv2_prec%ny       =   gvs_size%suv2%ny ! y coordi memory size of SUV2 IrradiancePrecision     
    gvs_size%suv2_prec%nz       =   gvs_size%suv2%nz ! z coordi memory size of SUV2 IrradiancePrecision     
    gvs_size%suv2_irrwav%nx     =   gvs_size%suv2%nx ! x coordi memory size of SUV2 Wavelength     
    gvs_size%suv2_irrwav%ny     =   gvs_size%suv2%ny ! y coordi memory size of SUV2 Wavelength     
    gvs_size%suv2_irrwav%nz     =   gvs_size%suv2%nz ! z coordi memory size of SUV2 Wavelength     

    gvs_size%svis_irrad%nx      =   gvs_size%svis%nx ! x coordi memory size of SVIS Irradiance     
    gvs_size%svis_irrad%ny      =   gvs_size%svis%ny ! y coordi memory size of SVIS Irradiance     
    gvs_size%svis_irrad%nz      =   gvs_size%svis%nz ! z coordi memory size of SVIS Irradiance     
    gvs_size%svis_irrwav%nx     =   gvs_size%svis%nx ! x coordi memory size of SVIS Wavelength     
    gvs_size%svis_irrwav%ny     =   gvs_size%svis%ny ! y coordi memory size of SVIS Wavelength     
    gvs_size%svis_irrwav%nz     =   gvs_size%svis%nz ! z coordi memory size of SVIS Wavelength     
    !---

    gvs_size%albedo%nx          =   gci_nxtrack2     ! x coordi memory size of Surface albedo 
    gvs_size%albedo%ny          =   gci_ntimes       ! y coordi memory size of Surface albedo 
    gvs_size%albedo%nz          =   gci_ntimes       ! z coordi memory size of Surface albedo 
    gvs_size%elevation%nx       =   gci_nxtrack2     ! x coordi memory size of Surface elevation        
    gvs_size%elevation%ny       =   gci_ntimes       ! y coordi memory size of Surface elevation        
    gvs_size%amfgeo%nx          =   gci_nxtrack2     ! x coordi memory size of geometric air mass factor
    gvs_size%amfgeo%ny          =   gci_ntimes       ! y coordi memory size of geometric air mass factor
    gvs_size%flg_ls%nx          =   gci_nxtrack2     ! x coordi memory size of Land/Sea Flag      
    gvs_size%flg_ls%ny          =   gci_ntimes       ! y coordi memory size of Land/Sea Flag      
    gvs_size%flg_cm%nx          =   gci_nxtrack2     ! x coordi memory size of Cloud Mask Flag    
    gvs_size%flg_cm%ny          =   gci_ntimes       ! y coordi memory size of Cloud Mask Flag    
    gvs_size%flg_px%nx          =   gci_nxtrack2     ! x coordi memory size of Pixel Availity Flag
    gvs_size%flg_px%ny          =   gci_ntimes       ! y coordi memory size of Pixel Availity Flag

    !
    ! initialize storage size of L1B File for C-Interface
    !
    rc = GEMS_Share_Hdf5InitL1BStorageSize(            &
      gci_nline            ,                           &
      gci_nwavel           ,                           &
      gci_nwavel2          ,                           &
      gci_nwavel3          ,                           &
      gci_nxtrack          ,                           &
      gci_nxtrack2         ,                           &
      gci_nwavelcoef       ,                           &
      gci_ntimes                                       &
    )

    !
    ! initialize storage size of L2 File for C-Interface
    !
    rc = GEMS_Share_Hdf5InitL2StorageSize(             &
              gci_nline_l2       ,                     &
              gci_nline2_l2      ,                     &
              gci_nwavel_l2      ,                     &
              gci_nwavel2_l2     ,                     &
              gci_nwavel3_l2     ,                     &
              gci_nwavel4_l2     ,                     &
              gci_nxtrack_l2     ,                     &
              gci_nxtrack2_l2    ,                     &
              gci_nlayer         ,                     &
              gci_nlayer2        ,                     &
              gci_nlayerp        ,                     &
              gci_nmon           ,                     &
              gci_nlat           ,                     &
              gci_nlon           ,                     &
              gci_ncomp          ,                     &
              gci_nResiduals     ,                     &
              gci_nColumns       ,                     &
              6                                        &   ! Time Column
    )

END SUBROUTINE GEMS_Share_Init_GlobalConstants2_Synt

SUBROUTINE GEMS_Share_Init_GlobalConstants3_Synt( ctl_fpath )
    USE NumTimes
    IMPLICIT NONE

    INTEGER :: rc
    CHARACTER(LEN=128),INTENT(IN) :: ctl_fpath

    INTEGER  GEMS_Share_Hdf5InitL1BStorageSize
    INTEGER  GEMS_Share_Hdf5InitL2StorageSize
    external GEMS_Share_Hdf5InitL1BStorageSize
    external GEMS_Share_Hdf5InitL2StorageSize

  ! Calculation Constants
    gcf_pi_c                = 3.141592        ! pi
    gcf_pi2_c               = gcf_pi_c*2.     ! pi*2
    gcf_deg2rad_c           = gcf_pi_c/180.   ! Degree to Radian
    gcf_rad2deg_c           = 180./gcf_pi_c   ! Radian to Degree

    gcf_planck_c            = 6.62606876E-34  ! Planck's constant  [J*s]  
    gcf_light_c             = 2.99793E8       ! Speed of light [m/s]
    gcf_boltz_c             = 1.380622E-23    ! Boltzmann constant [J/K]

    gcf_solar_c             = 1367.0          ! Solar constant     [W/m^2]
    gcf_press_c             = 1013.25         ! standard pressure  [hPa]  
    gcf_kelvin_c            = 273.16          ! Absolute Temperature[K]
      
  ! Land/Sea Flag     
    gci_flg_ls_sea          = 0               ! Land/Sea Flag; Sea
    gci_flg_ls_land         = 1               ! Land/Sea Flag; Land
    gci_flg_ls_coast        = 2               ! Land/Sea Flag; Coast
    gci_flg_ls_unavail      = -999            ! Land/Sea Flag; Unavailable
    gci_flg_ls_exosphere    = -999            ! Land/Sea Flag; Exospher


  ! Size for L1B dataset
    CALL GEMS_NumTimes_Read3(ctl_fpath)
    gci_nline_rug           = nlines%n_line_rug
    gci_nline_rvg           = nlines%n_line_rvg
    gci_nline_irr           = nlines%n_line_irr
    gci_nwavel              = 159
    gci_nwavel2             = 557
    gci_nwavel3             = 751
    gci_nxtrack             = 30
    gci_nxtrack2            = 60
    gci_nwavelcoef          = 5
    gci_ntimes              = gci_nline

 ! Size for L2 dataset
    gci_nline_rug_l2           = nlines%n_line_rug     !변경  1644에서 nline으로 
    gci_nline_rvg_l2           = nlines%n_line_rvg     !변경  1644에서 nline으로 
    gci_nline_irr_l2           = nlines%n_line_irr     !변경  1644에서 nline으로 
    gci_nline2_l2              = 330
    gci_nwavel_l2              = 3
    gci_nwavel2_l2             = 12
    gci_nwavel3_l2             = 20
    gci_nwavel4_l2             = 23
    gci_nxtrack_l2             = 60
    gci_nxtrack2_l2            = 30
    gci_nlayer                 = 24
    gci_nlayer2                = 11
!    gci_nlayerp = 777
    gci_nlayerp                = gci_nlayer+1
    gci_nmon                   = 12
    gci_nlat                   = 360
    gci_nlon                   = 720
    gci_ncomp                  = 2
    !
    !--- Size Declaration for Global Allocatable Variables
    !
    gvs_size%angles%nx          =   gci_nxtrack2     ! x coordi memory size of Angles of vza, sza, raa
    gvs_size%angles%ny          =   gci_ntimes       ! y coordi memory size of Angles of vza, sza, raa
    gvs_size%latlon%nx          =   gci_nxtrack2     ! x coordi memory size of Latitude/Longitude 
    gvs_size%latlon%ny          =   gci_ntimes       ! y coordi memory size of Latitude/Longitude 

    !
    !--- Dimension Size for L1B datasets
    !
    gvs_size%euv1%nx            =   gci_nwavel       ! x coordi memory size of EUV1       
    gvs_size%euv1%ny            =   gci_nxtrack      ! y coordi memory size of EUV1       
    gvs_size%euv1%nz            =   gci_ntimes       ! z coordi memory size of EUV1       
			       
    gvs_size%euv2%nx            =   gci_nwavel2      ! x coordi memory size of EUV2       
    gvs_size%euv2%ny            =   gci_nxtrack2     ! y coordi memory size of EUV2       
    gvs_size%euv2%nz            =   gci_ntimes       ! z coordi memory size of EUV2       
			       
    gvs_size%evis%nx            =   gci_nwavel3      ! x coordi memory size of EVIS     
    gvs_size%evis%ny            =   gci_nxtrack2     ! y coordi memory size of EVIS     
    gvs_size%evis%nz            =   gci_ntimes       ! z coordi memory size of EVIS     

    gvs_size%suv1%nx            =   gci_nwavel       ! x coordi memory size of SUV1     
    gvs_size%suv1%ny            =   gci_nxtrack      ! y coordi memory size of SUV1     
    gvs_size%suv1%nz            =   gci_ntimes       ! z coordi memory size of SUV1     
			       
    gvs_size%suv2%nx            =   gci_nwavel2      ! x coordi memory size of SUV2     
    gvs_size%suv2%ny            =   gci_nxtrack2     ! y coordi memory size of SUV2     
    gvs_size%suv2%nz            =   gci_ntimes       ! z coordi memory size of SUV2     
			       
    gvs_size%svis%nx            =   gci_nwavel3      ! x coordi memory size of SVIS     
    gvs_size%svis%ny            =   gci_nxtrack2     ! y coordi memory size of SVIS     
    gvs_size%svis%nz            =   gci_ntimes       ! z coordi memory size of SVIS     
    !---

    !
    !--- Dimension size  of radiance, wavelength, irradiance based on OMI
    !                   nx          ny          nz
    !                  ---         ---        ----
    ! UV1_rad:         159          30        1644    gci_nwavel    gci_nxtrack     gci_ntimes
    ! wave:            159          30        1644
    !--
    ! UV2_rad:         557          60        1644    gci_nwavel2   gci_nxtrack2    gci_ntimes
    ! wave:            557          60        1644
    !--
    ! VIS_rad:         751          60        1644    gci_nwavel3   gci_nxtrack2    gci_ntimes
    ! wave:            751          60        1644
    !--
    ! UV1_sun_rad:     159          30                gci_nwavel    gci_nxtrack
    ! wave             159          30
    !--
    ! UV2_sun_rad:     557          60                gci_nwavel2   gci_nxtrack2   
    ! wave             557          60
    ! prec             557          60
    !--
    ! VIS_sun_rad      751          60                gci_nwavel3   gci_nxtrack2  
    ! wave             751          60
    !
    !------------------------------------------------------------------------------------------

    !
    !--- Dimension Size for Radiance, Wavelength, Irradiance datasets
    !
    gvs_size%euv1_rad%nx        =   gvs_size%euv1%nx ! x coordi memory size of EUV1 Radiance       
    gvs_size%euv1_rad%ny        =   gvs_size%euv1%ny ! y coordi memory size of EUV1 Radiance       
    gvs_size%euv1_rad%nz        =   gvs_size%euv1%nz ! z coordi memory size of EUV1 Radiance       
    gvs_size%euv1_radwav%nx     =   gvs_size%euv1%nx ! x coordi memory size of EUV1 Wavelength       
    gvs_size%euv1_radwav%ny     =   gvs_size%euv1%ny ! y coordi memory size of EUV1 Wavelength       
    gvs_size%euv1_radwav%nz     =   gvs_size%euv1%nz ! z coordi memory size of EUV1 Wavelength       
    gvs_size%euv1_prec%nx       =   gvs_size%euv1%nx ! x coordi memory size of EUV1 RadiancePrecision     
    gvs_size%euv1_prec%ny       =   gvs_size%euv1%ny ! y coordi memory size of EUV1 RadiancePrecision     
    gvs_size%euv1_prec%nz       =   gvs_size%euv1%nz ! z coordi memory size of EUV1 RadiancePrecision

    gvs_size%euv2_rad%nx        =   gvs_size%euv2%nx ! x coordi memory size of EUV2 Radiance       
    gvs_size%euv2_rad%ny        =   gvs_size%euv2%ny ! y coordi memory size of EUV2 Radiance       
    gvs_size%euv2_rad%nz        =   gvs_size%euv2%nz ! z coordi memory size of EUV2 Radiance       
    gvs_size%euv2_radwav%nx     =   gvs_size%euv2%nx ! x coordi memory size of EUV2 Wavelength       
    gvs_size%euv2_radwav%ny     =   gvs_size%euv2%ny ! y coordi memory size of EUV2 Wavelength       
    gvs_size%euv2_radwav%nz     =   gvs_size%euv2%nz ! z coordi memory size of EUV2 Wavelength       
    gvs_size%euv2_prec%nx       =   gvs_size%euv2%nx ! x coordi memory size of EUV2 RadiancePrecision     
    gvs_size%euv2_prec%ny       =   gvs_size%euv2%ny ! y coordi memory size of EUV2 RadiancePrecision     
    gvs_size%euv2_prec%nz       =   gvs_size%euv2%nz ! z coordi memory size of EUV2 RadiancePrecision


    gvs_size%evis_rad%nx        =   gvs_size%evis%nx ! x coordi memory size of EVIS Radiance       
    gvs_size%evis_rad%ny        =   gvs_size%evis%ny ! y coordi memory size of EVIS Radiance       
    gvs_size%evis_rad%nz        =   gvs_size%evis%nz ! z coordi memory size of EVIS Radiance       
    gvs_size%evis_radwav%nx     =   gvs_size%evis%nx ! x coordi memory size of EVIS Wavelength       
    gvs_size%evis_radwav%ny     =   gvs_size%evis%ny ! y coordi memory size of EVIS Wavelength       
    gvs_size%evis_radwav%nz     =   gvs_size%evis%nz ! z coordi memory size of EVIS Wavelength       

    gvs_size%suv1_irrad%nx      =   gvs_size%suv1%nx ! x coordi memory size of SUV1 Irradiance     
    gvs_size%suv1_irrad%ny      =   gvs_size%suv1%ny ! y coordi memory size of SUV1 Irradiance     
    gvs_size%suv1_irrad%nz      =   gvs_size%suv1%nz ! z coordi memory size of SUV1 Irradiance     
    gvs_size%suv1_irrwav%nx     =   gvs_size%suv1%nx ! x coordi memory size of SUV1 Wavelength     
    gvs_size%suv1_irrwav%ny     =   gvs_size%suv1%ny ! y coordi memory size of SUV1 Wavelength     
    gvs_size%suv1_irrwav%nz     =   gvs_size%suv1%nz ! z coordi memory size of SUV1 Wavelength     

    gvs_size%suv2_irrad%nx      =   gvs_size%suv2%nx ! x coordi memory size of SUV2 Irradiance     
    gvs_size%suv2_irrad%ny      =   gvs_size%suv2%ny ! y coordi memory size of SUV2 Irradiance     
    gvs_size%suv2_irrad%nz      =   gvs_size%suv2%nz ! z coordi memory size of SUV2 Irradiance     
    gvs_size%suv2_prec%nx       =   gvs_size%suv2%nx ! x coordi memory size of SUV2 IrradiancePrecision     
    gvs_size%suv2_prec%ny       =   gvs_size%suv2%ny ! y coordi memory size of SUV2 IrradiancePrecision     
    gvs_size%suv2_prec%nz       =   gvs_size%suv2%nz ! z coordi memory size of SUV2 IrradiancePrecision     
    gvs_size%suv2_irrwav%nx     =   gvs_size%suv2%nx ! x coordi memory size of SUV2 Wavelength     
    gvs_size%suv2_irrwav%ny     =   gvs_size%suv2%ny ! y coordi memory size of SUV2 Wavelength     
    gvs_size%suv2_irrwav%nz     =   gvs_size%suv2%nz ! z coordi memory size of SUV2 Wavelength     

    gvs_size%svis_irrad%nx      =   gvs_size%svis%nx ! x coordi memory size of SVIS Irradiance     
    gvs_size%svis_irrad%ny      =   gvs_size%svis%ny ! y coordi memory size of SVIS Irradiance     
    gvs_size%svis_irrad%nz      =   gvs_size%svis%nz ! z coordi memory size of SVIS Irradiance     
    gvs_size%svis_irrwav%nx     =   gvs_size%svis%nx ! x coordi memory size of SVIS Wavelength     
    gvs_size%svis_irrwav%ny     =   gvs_size%svis%ny ! y coordi memory size of SVIS Wavelength     
    gvs_size%svis_irrwav%nz     =   gvs_size%svis%nz ! z coordi memory size of SVIS Wavelength     
    !---

    gvs_size%albedo%nx          =   gci_nxtrack2     ! x coordi memory size of Surface albedo 
    gvs_size%albedo%ny          =   gci_ntimes       ! y coordi memory size of Surface albedo 
    gvs_size%albedo%nz          =   gci_ntimes       ! z coordi memory size of Surface albedo 
    gvs_size%elevation%nx       =   gci_nxtrack2     ! x coordi memory size of Surface elevation        
    gvs_size%elevation%ny       =   gci_ntimes       ! y coordi memory size of Surface elevation        
    gvs_size%amfgeo%nx          =   gci_nxtrack2     ! x coordi memory size of geometric air mass factor
    gvs_size%amfgeo%ny          =   gci_ntimes       ! y coordi memory size of geometric air mass factor
    gvs_size%flg_ls%nx          =   gci_nxtrack2     ! x coordi memory size of Land/Sea Flag      
    gvs_size%flg_ls%ny          =   gci_ntimes       ! y coordi memory size of Land/Sea Flag      
    gvs_size%flg_cm%nx          =   gci_nxtrack2     ! x coordi memory size of Cloud Mask Flag    
    gvs_size%flg_cm%ny          =   gci_ntimes       ! y coordi memory size of Cloud Mask Flag    
    gvs_size%flg_px%nx          =   gci_nxtrack2     ! x coordi memory size of Pixel Availity Flag
    gvs_size%flg_px%ny          =   gci_ntimes       ! y coordi memory size of Pixel Availity Flag

    !
    ! initialize storage size of L1B File for C-Interface
    !
    rc = GEMS_Share_Hdf5InitL1BStorageSize(            &
      gci_nline            ,                           &
      gci_nwavel           ,                           &
      gci_nwavel2          ,                           &
      gci_nwavel3          ,                           &
      gci_nxtrack          ,                           &
      gci_nxtrack2         ,                           &
      gci_nwavelcoef       ,                           &
      gci_ntimes                                       &
    )

    !
    ! initialize storage size of L2 File for C-Interface
    !
    rc = GEMS_Share_Hdf5InitL2StorageSize(             &
              gci_nline_l2       ,                     &
              gci_nline2_l2      ,                     &
              gci_nwavel_l2      ,                     &
              gci_nwavel2_l2     ,                     &
              gci_nwavel3_l2     ,                     &
              gci_nwavel4_l2     ,                     &
              gci_nxtrack_l2     ,                     &
              gci_nxtrack2_l2    ,                     &
              gci_nlayer         ,                     &
              gci_nlayer2        ,                     &
              gci_nlayerp        ,                     &
              gci_nmon           ,                     &
              gci_nlat           ,                     &
              gci_nlon           ,                     &
              gci_ncomp          ,                     &
              gci_nResiduals     ,                     &
              gci_nColumns       ,                     &
              6                                        &   ! Time Column
    )

END SUBROUTINE GEMS_Share_Init_GlobalConstants3_Synt

END MODULE Share_MOD_Constants_Variables_Synt

!-------------------------------------------------------------------------------
