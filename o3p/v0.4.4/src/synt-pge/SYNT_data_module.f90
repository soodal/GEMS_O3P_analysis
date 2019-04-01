MODULE SYNT_data_module

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: maxchlen, n_rad_winwav, maxwin, max_fit_pts, max_ring_pts, mrefl
  USE GEMS_O3P_gemsdata_module, ONLY: nxtrack_max
  !USE OMSAO_indices_module,    ONLY: n_max_fitpars, max_rs_idx, max_calfit_idx, o3_t1_idx, o3_t3_idx, spc_idx

  IMPLICIT NONE

  ! Define Dimensions 
  INTEGER, PARAMETER :: mchannel=1, mswath = 1, nytrack_max =899
  !INTEGER, PARAMETER :: nswath=1
  INTEGER :: nxtrack, nytrack, nwavel  
  !CHARACTER (LEN=5)        :: orbc, orbcsol
  !INTEGER (KIND=4)         :: orbnum, orbnumsol
  !INTEGER :: which_omps_cloud, which_omps

  !! ------------------------------------------------------------
  !! Boundary wavelengths (approximate) 
  !! ------------------------------------------------------------
  !REAL (KIND=4), DIMENSION(mswath),  PARAMETER :: lower_wvls = (/290/), upper_wvls = (/510/)
  !REAL (KIND=8), DIMENSION(mswath),  PARAMETER :: lower_spec = (/1E8/), upper_spec = (/4.0E15/) 
  REAL (KIND=8), DIMENSION(mswath),  PARAMETER :: lower_spec = (/0.0/), upper_spec = (/4.0E15/)                               
  !REAL (KIND=8), DIMENSION(mswath)             :: reduce_ubnd, reduce_lbnd, retubnd, retlbnd


 
  ! ------------------------------------------------------------
  ! OMPS origianl Data
  ! ------------------------------------------------------------
  !INTEGER, POINTER, DIMENSION (:)        :: omps_IQFlag
  !REAL, POINTER, DIMENSION (:)           :: omps_exposuretime
  !REAL (KIND=DP), POINTER, DIMENSION (:) :: omps_time
  !INTEGER, POINTER, DIMENSION (:,:)      :: omps_GPQFlag

  !INTEGER, POINTER, DIMENSION (:,:)   :: omps_irrad_qflg
  REAL, POINTER, DIMENSION (:,:)      :: synt_irrad_wavl,synt_irrad_spec !, omps_irrad_prec, omps_rad_coeff

  INTEGER, POINTER, DIMENSION (:,:,:) :: synt_rad_qflg
  REAL, POINTER, DIMENSION (:,:,:) :: synt_rad_spec, synt_rad_wavl, synt_rad_snr, synt_surfalb !,omps_rad_prec, omps_stray_corr, omps_smear_corr


  REAL, POINTER, DIMENSION (:,:)      :: synt_latitude, synt_Longitude,& 
                                         synt_vzenith,  synt_saa, synt_vaa, &
                                         synt_szenith,  synt_hter


  CHARACTER (LEN=maxchlen)            :: synt_db
  CHARACTER (LEN=3)                   :: synt_num
  CHARACTER (LEN=10)                  :: synt_date

  !REAL :: omps_sunearth

!!-----------------------
!! Array for processed geolocation dataset
!!-----------------------
  TYPE geo_group
  REAL (KIND=r8), DIMENSION (:)  , POINTER    :: time
  REAL (KIND=r8), DIMENSION (:,:), POINTER    :: lon, lat, ctp, height
  REAL (KIND=r8), DIMENSION (:,:), POINTER    :: sza, vza, aza,sca
  REAL (KIND=r8), DIMENSION (:,:), POINTER    :: saz, vaz
  REAL (KIND=r8), DIMENSION (:,:), POINTER    :: clat,clon
  REAL (KIND=r8), DIMENSION (:,:), POINTER    :: elat,elon
  REAL (KIND=r8), DIMENSION (:,:), POINTER    :: groundqflg
  END TYPE geo_group
  TYPE (geo_group) :: synt_geo


  !! ------------------------------------------------------------
  !! OMPS processed Data
  !! ------------------------------------------------------------
  !TYPE IRRAD_GROUP
     !INTEGER, DIMENSION(nxtrack_max) :: nwav,errstat
     !INTEGER, DIMENSION(max_fit_pts, nxtrack_max) :: qflg
     !INTEGER (KINd=2), DIMENSION(max_fit_pts,nxtrack_max) :: wind
     !INTEGER, DIMENSION (maxwin, nxtrack_max) :: npix
     !REAL (KIND=8), DIMENSION(nxtrack_max) :: norm
     !REAL (KIND=8), DIMENSION(maxwin, nxtrack_max,2) :: winpix
     !REAL,    DIMENSION(max_fit_pts, nxtrack_max):: spec, prec, wavl     
  !END TYPE IRRAD_GROUP
  !TYPE (IRRAD_GROUP) :: omps_irrad

  TYPE RAD_GROUP
    !INTEGER :: nxtrack, nwavel
    !INTEGER (KIND=2), POINTER, DIMENSION (:)        :: mqflg
    !INTEGER, DIMENSION(nxtrack_max, nytrack_max) :: nwav,errstat
    !INTEGER, DIMENSION(max_fit_pts, 10, nytrack_max) :: qflg
    !INTEGER (KINd=2), DIMENSION(max_fit_pts,10,nytrack_max) :: wind
    !INTEGER, DIMENSION (maxwin, nxtrack_max,nytrack_max) :: npix
    !REAL (KIND=8), DIMENSION(nxtrack_max, nytrack_max) :: norm
    !REAL (KIND=8), DIMENSION(maxwin, nxtrack_max,nytrack_max,2) :: winpix
    !REAL,    DIMENSION(max_fit_pts, nxtrack_max, nytrack_max):: spec, prec, wavl,stray,smear   

    INTEGER (KIND=2), POINTER, DIMENSION (:)       :: mqflg
    INTEGER,          POINTER, DIMENSION (:, :)    :: nwav,errstat
    INTEGER,          POINTER, DIMENSION (:, :, :) :: qflg
    INTEGER (KINd=2), POINTER, DIMENSION (:,:,:)   :: wind
    INTEGER,          POINTER, DIMENSION (:, :,:)  :: npix
    REAL (KIND=8),    POINTER, DIMENSION (:, :)    :: norm
    REAL (KIND=8),    POINTER, DIMENSION (:,:,:,:) :: winpix
    REAL,             POINTER, DIMENSION (:, :, :) :: spec, prec, wavl,stray,smear   
    !INTEGER (KIND=i8)                               :: nxtrack, nwavel   : geun closed
  END TYPE RAD_GROUP
  TYPE (RAD_GROUP) :: synt_rad


  !TYPE ring_group
  !REAL    (KIND=4), DIMENSION ( max_ring_pts, nxtrack_max)   :: solspec, solwavl
  !INTEGER, DIMENSION (nxtrack_max)                           :: sol_lin, sol_uin, nsol, ndiv
  !END TYPE ring_group
  !TYPE (ring_group) :: omps_ring


!!-----------------------
!! Array for reflectance
!!-----------------------
  !TYPE reflectance
     !REAL (KIND=r4), DIMENSION (mrefl, nxtrack_max) :: solspec
     !REAL (KIND=r4), DIMENSION (mrefl, nxtrack_max) :: solwavl
     !REAL (KIND=r4), DIMENSION (mrefl, nxtrack_max,nytrack_max) :: radspec
     !REAL (KINd=r4), DIMENSION (mrefl, nxtrack_max,nytrack_max) :: radwavl
     !REAL (KIND=r4), DIMENSION (nxtrack_max, 2) :: solwinpix
  !END TYPE reflectance
  !TYPE (reflectance) ::  omps_refl
  !REAL (KIND=r4), DIMENSION (spc_idx, mrefl, nxtrack_max) ::omps_solspecr
  !REAL (KIND=r4), DIMENSION (nxtrack_max,2) :: omps_solrwinpix



!!-----------------------
!! Array for ground flags
!!-----------------------
!TYPE gflag_group
  !INTEGER (KIND=I2), DIMENSION (nytrack_max)         :: saa
  !INTEGER (KIND=i2), DIMENSION (nxtrack_max,nytrack_max) :: land_water
  !INTEGER (KIND=i2), DIMENSION (nxtrack_max,nytrack_max) :: glint
  !INTEGER (KIND=i2), DIMENSION (nxtrack_max,nytrack_max) :: snow_ice
!END TYPE gflag_group
!TYPE (gflag_group) :: omps_gflag
  INTEGER (KIND=i2), DIMENSION (nxtrack_max, nytrack_max) :: land_water_flg, snow_ice_flg, glint_flg

!!-----------------------
!! calibration varialbes
!!-----------------------
  !REAL (KIND=dp), DIMENSION(maxwin)                                   :: omps_redslw
  !REAL (KIND=dp), DIMENSION(maxwin, nxtrack_max)                      :: omps_wincal_wav
  !REAL (KIND=dp), DIMENSION(maxwin,max_calfit_idx,2, nxtrack_max)     :: omps_solwinfit

  !REAL (KIND=dp), DIMENSION (max_fit_pts, nxtrack_max) :: omps_slitwav_sol=0.0, omps_sswav_sol=0.0, &
                                                          !omps_sswav_rad=0.0 
  !REAL (KIND=dp), DIMENSION (max_fit_pts, max_calfit_idx, 2, nxtrack_max)  ::omps_solslitfit=0.0,  &
       !omps_solwavfit=0.0, omps_radwavfit=0.0
  !INTEGER, DIMENSION (nxtrack_max) ::omps_nslit_sol, omps_nwavcal_sol, omps_nwavcal_rad

  !REAL (KIND=dp), DIMENSION (:,:,:,:), POINTER  :: omps_cali_sol

  !! radiance calibration 
  !INTEGER, DIMENSION(nxtrack_max, nytrack_max) :: omps_nslit_rad
  !REAL (KIND=DP), DIMENSION(max_fit_pts, nxtrack_max, nytrack_max) :: omps_slitwav_rad
  !REAL (KIND=dp), DIMENSION(maxwin, max_calfit_idx, 2, nxtrack_max, nytrack_max) :: omps_radwinfit
  !REAL (KIND=dp), DIMENSION(:,:,:,:,:), POINTER :: omps_radslitfit
  !REAL (KIND=dp), DIMENSION (:,:,:,:,:), POINTER  :: omps_cali_rad


!-----------------------
! OUTPUT variables
!-----------------------
  !INTEGER, DIMENSION (nxtrack_max, nytrack_max)                          :: omps_exitval, omps_initval
!REAL (KIND=dp),   DIMENSION (nxtrack_max, nytrack_max,  n_max_fitpars)   :: omps_fitvar 

  CONTAINS
  SUBROUTINE allocate_synt_raddata (nx, ny, nw, max_fit_pts, maxwin, status)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nx, ny, nw, max_fit_pts, maxwin
    INTEGER, INTENT(OUT) :: status

!!-----------------------
!! Array for radiance dataset
!!-----------------------
     !INTEGER :: nxtrack, nwavel
     !INTEGER (KIND=2), POINTER, DIMENSION (:)        :: mqflg
     !INTEGER, DIMENSION(nxtrack_max, nytrack_max) :: nwav,errstat
     !INTEGER, DIMENSION(max_fit_pts, 10, nytrack_max) :: qflg
     !INTEGER (KINd=2), DIMENSION(max_fit_pts,10,nytrack_max) :: wind
     !INTEGER, DIMENSION (maxwin, nxtrack_max,nytrack_max) :: npix
     !REAL (KIND=8), DIMENSION(nxtrack_max, nytrack_max) :: norm
     !REAL (KIND=8), DIMENSION(maxwin, nxtrack_max,nytrack_max,2) :: winpix
     !REAL,    DIMENSION(max_fit_pts, nxtrack_max, nytrack_max):: spec, prec, wavl,stray,smear   

   IF (ASSOCIATED (synt_rad%nwav)) DEALLOCATE (synt_rad%nwav)
   ALLOCATE (synt_rad%nwav(nx, ny), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_rad%errstat)) DEALLOCATE (synt_rad%errstat)
   ALLOCATE (synt_rad%errstat(nx, ny), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_rad%qflg)) DEALLOCATE (synt_rad%qflg)
   ALLOCATE (synt_rad%qflg(max_fit_pts, nx, ny), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_rad%wind)) DEALLOCATE (synt_rad%wind)
   ALLOCATE (synt_rad%wind(max_fit_pts, nx, ny), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_rad%npix)) DEALLOCATE (synt_rad%npix)
   ALLOCATE (synt_rad%npix(maxwin, nx, ny), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_rad%norm)) DEALLOCATE (synt_rad%norm)
   ALLOCATE (synt_rad%norm(nx, ny), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_rad%winpix)) DEALLOCATE (synt_rad%winpix)
   ALLOCATE (synt_rad%winpix(maxwin, nx, ny, 2), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_rad%spec)) DEALLOCATE (synt_rad%spec)
   ALLOCATE (synt_rad%spec(max_fit_pts, nx, ny), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_rad%prec)) DEALLOCATE (synt_rad%prec)
   ALLOCATE (synt_rad%prec(max_fit_pts, nx, ny), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_rad%wavl)) DEALLOCATE (synt_rad%wavl)
   ALLOCATE (synt_rad%wavl(max_fit_pts, nx, ny), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_rad%stray)) DEALLOCATE (synt_rad%stray)
   ALLOCATE (synt_rad%stray(max_fit_pts, nx, ny), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_rad%smear)) DEALLOCATE (synt_rad%smear)
   ALLOCATE (synt_rad%smear(max_fit_pts, nx, ny), stat=status); IF (status .ne. 0 ) GOTO 111

   RETURN

   111 continue
   status  = -1
   return  
  END SUBROUTINE allocate_synt_raddata


  SUBROUTINE allocate_synt_data (nx, ny, nw, status)
  IMPLICIT NONE
   INTEGER, INTENT(IN) :: nx, ny, nw
   INTEGER, INTENT(OUT) :: status
  ! ------------------------------------------------------------
  ! OMPS origianl Data
  ! ------------------------------------------------------------


   !IF (ASSOCIATED (omps_IQFlag)) DEALLOCATE (omps_IQFlag)
   !ALLOCATE (omps_IQFlag(ny), stat=status); IF (status .ne. 0 ) GOTO 111

   !IF (ASSOCIATED (omps_exposuretime)) DEALLOCATE (omps_exposuretime)
   !ALLOCATE (omps_exposuretime(ny), stat=status); IF (status .ne. 0 ) GOTO 111

   !IF (ASSOCIATED (omps_time)) DEALLOCATE (omps_time)
   !ALLOCATE (omps_time(ny), stat=status); IF (status .ne. 0 ) GOTO 111

   !IF (ASSOCIATED (omps_GPQFlag)) DEALLOCATE (omps_GPQFlag)
   !ALLOCATE (omps_GPQFlag(nx, ny), stat=status); IF (status .ne. 0 ) GOTO 111
   
   !IF (ASSOCIATED (omps_irrad_qflg)) DEALLOCATE (omps_irrad_qflg)
   !ALLOCATE (omps_irrad_qflg(nw, nx), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_irrad_wavl)) DEALLOCATE (synt_irrad_wavl)
   ALLOCATE (synt_irrad_wavl(nw, nx), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_irrad_spec)) DEALLOCATE (synt_irrad_spec)
   ALLOCATE (synt_irrad_spec(nw, nx), stat=status); IF (status .ne. 0 ) GOTO 111

   !IF (ASSOCIATED (omps_irrad_prec)) DEALLOCATE (omps_irrad_prec)
   !ALLOCATE (omps_irrad_prec(nw, nx), stat=status); IF (status .ne. 0 ) GOTO 111

   !IF (ASSOCIATED (omps_rad_coeff)) DEALLOCATE (omps_rad_coeff)
   !ALLOCATE (omps_rad_coeff(nw, nx), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_rad_qflg)) DEALLOCATE (synt_rad_qflg)
   ALLOCATE (synt_rad_qflg(nw, nx, ny), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_rad_spec)) DEALLOCATE (synt_rad_spec)
   ALLOCATE (synt_rad_spec(nw, nx, ny), stat=status); IF (status .ne. 0 ) GOTO 111

   !IF (ASSOCIATED (synt_rad_prec)) DEALLOCATE (synt_rad_prec)
   !ALLOCATE (synt_rad_prec(nw, nx, ny), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_rad_wavl)) DEALLOCATE (synt_rad_wavl)
   ALLOCATE (synt_rad_wavl(nw, nx, ny), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_rad_snr)) DEALLOCATE (synt_rad_snr)
   ALLOCATE (synt_rad_snr(nw, nx, ny), stat=status); IF (status .ne. 0 ) GOTO 111

   !IF (ASSOCIATED (omps_stray_corr)) DEALLOCATE (omps_stray_corr)
   !ALLOCATE (omps_stray_corr(nw, nx, ny), stat=status); IF (status .ne. 0 ) GOTO 111

   !IF (ASSOCIATED (omps_smear_corr)) DEALLOCATE (omps_smear_corr)
   !ALLOCATE (omps_smear_corr(nw, nx, ny), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_latitude)) DEALLOCATE (synt_latitude)
   ALLOCATE (synt_latitude(nx, ny), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_Longitude)) DEALLOCATE (synt_Longitude)
   ALLOCATE (synt_Longitude(nx, ny), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_hter)) DEALLOCATE (synt_hter)
   ALLOCATE (synt_hter(nx, ny), stat=status); IF (status .ne. 0 ) GOTO 111

   !IF (ASSOCIATED (omps_vazimuth)) DEALLOCATE (omps_vazimuth)
   !ALLOCATE (omps_vazimuth(nx, ny), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_vzenith)) DEALLOCATE (synt_vzenith)
   ALLOCATE (synt_vzenith(nx, ny), stat=status); IF (status .ne. 0 ) GOTO 111

   !IF (ASSOCIATED (omps_sazimuth)) DEALLOCATE (omps_sazimuth)
   !ALLOCATE (omps_sazimuth(nx, ny), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_szenith)) DEALLOCATE (synt_szenith)
   ALLOCATE (synt_szenith(nx, ny), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_saa)) DEALLOCATE (synt_saa)
   ALLOCATE (synt_saa(nx, ny), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_vaa)) DEALLOCATE (synt_vaa)
   ALLOCATE (synt_vaa(nx, ny), stat=status); IF (status .ne. 0 ) GOTO 111

   !IF (ASSOCIATED (omps_ctp)) DEALLOCATE (omps_ctp)
   !ALLOCATE (omps_ctp(nx, ny), stat=status); IF (status .ne. 0 ) GOTO 111

   !IF (ASSOCIATED (omps_snowicefrc)) DEALLOCATE (omps_snowicefrc)
   !ALLOCATE (omps_snowicefrc(nx, ny), stat=status); IF (status .ne. 0 ) GOTO 111

   !IF (ASSOCIATED (omps_albedo)) DEALLOCATE (omps_albedo)
   !ALLOCATE (omps_albedo(nx, ny), stat=status); IF (status .ne. 0 ) GOTO 111

   !IF (ASSOCIATED (omps_spres)) DEALLOCATE (omps_spres)
   !ALLOCATE (omps_spres(nx, ny), stat=status); IF (status .ne. 0 ) GOTO 111

!!-----------------------
!! Array for alb dataset
!!-----------------------

   IF (ASSOCIATED (synt_surfalb)) DEALLOCATE (synt_surfalb)
   ALLOCATE (synt_surfalb(nx, ny, nw), stat=status); IF (status .ne. 0 ) GOTO 111

!!-----------------------
!! Array for processed geolocation dataset
!!-----------------------

   IF (ASSOCIATED (synt_geo%time)) DEALLOCATE (synt_geo%time)
   ALLOCATE (synt_geo%time(ny), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_geo%lon)) DEALLOCATE (synt_geo%lon)
   ALLOCATE (synt_geo%lon(nx,ny), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_geo%lat)) DEALLOCATE (synt_geo%lat)
   ALLOCATE (synt_geo%lat(nx,ny), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_geo%ctp)) DEALLOCATE (synt_geo%ctp)
   ALLOCATE (synt_geo%ctp(nx,ny), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_geo%height)) DEALLOCATE (synt_geo%height)
   ALLOCATE (synt_geo%height(nx,ny), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_geo%sza)) DEALLOCATE (synt_geo%sza)
   ALLOCATE (synt_geo%sza(nx,ny), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_geo%vza)) DEALLOCATE (synt_geo%vza)
   ALLOCATE (synt_geo%vza(nx,ny), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_geo%saz)) DEALLOCATE (synt_geo%saz)
   ALLOCATE (synt_geo%saz(nx,ny), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_geo%vaz)) DEALLOCATE (synt_geo%vaz)
   ALLOCATE (synt_geo%vaz(nx,ny), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_geo%aza)) DEALLOCATE (synt_geo%aza)
   ALLOCATE (synt_geo%aza(nx,ny), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_geo%sca)) DEALLOCATE (synt_geo%sca)
   ALLOCATE (synt_geo%sca(nx,ny), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_geo%clat)) DEALLOCATE (synt_geo%clat)
   ALLOCATE (synt_geo%clat(0:nx,1:ny+1), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_geo%clon)) DEALLOCATE (synt_geo%clon)
   ALLOCATE (synt_geo%clon(0:nx,1:ny+1), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_geo%elon)) DEALLOCATE (synt_geo%elon)
   ALLOCATE (synt_geo%elon(0:nx,ny), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_geo%elat)) DEALLOCATE (synt_geo%elat)
   ALLOCATE (synt_geo%elat(0:nx,ny), stat=status); IF (status .ne. 0 ) GOTO 111

   IF (ASSOCIATED (synt_geo%groundqflg)) DEALLOCATE (synt_geo%groundqflg)
   ALLOCATE (synt_geo%groundqflg(nx,ny), stat=status); IF (status .ne. 0 ) GOTO 111

!!-------------------
!!  calibration data
!!-------------------

   !IF (ASSOCIATED (omps_radslitfit)) DEALLOCATE (omps_radslitfit)
   !ALLOCATE (omps_radslitfit(max_fit_pts, max_calfit_idx,2,nx,ny), stat=status); IF (status .ne. 0 ) GOTO 111

   !IF (ASSOCIATED (omps_cali_sol)) DEALLOCATE (omps_cali_sol)
   !ALLOCATE (omps_cali_sol(max_fit_pts,3, max_fit_pts,nx), stat=status); IF (status .ne. 0 ) GOTO 111
   !IF (ASSOCIATED (omps_cali_rad)) DEALLOCATE (omps_cali_rad)
   !ALLOCATE (omps_cali_rad(30,3, max_fit_pts, nx, ny), stat=status); IF (status .ne. 0 ) GOTO 111

   RETURN

   111 continue
   status  = -1
   return  
  END SUBROUTINE allocate_synt_data
  

END MODULE SYNT_data_module


