MODULE GEMS_O3P_gemsdata_module

USE OMSAO_precision_module
USE OMSAO_parameters_module, ONLY: maxwin, mrefl, max_ring_pts,max_fit_pts
USE OMSAO_indices_module,    ONLY: max_calfit_idx, n_max_fitpars


CHARACTER (LEN=5)         :: orbc, orbcsol
CHARACTER (LEN=9)         :: gemsraddate
INTEGER (KIND=i4)         :: orbnum, orbnumsol !, gemssol_version


! Define Data Maximum Dimension
INTEGER (KIND=i4), PARAMETER :: ntimes_max= 2000, nxtrack_max  = 60, nwavel_max  = 720 !wasp
INTEGER (KIND=i4), PARAMETER :: mchannel=2 ! # of maximum channel , mswath
INTEGER :: nchannel, nxtrack, nfxtrack, ntimes  !nfxtrack = nxtrakc for UV1 if selected
INTEGER, DIMENSION (mchannel):: chs
INTEGER :: nxbin, nybin, ncoadd,first_pix, last_pix, first_line, last_line, offline, gems_nx, gems_ny
INTEGER, DIMENSION (ntimes_max) ::  lineloc
INTEGER, DIMENSION (nxtrack_max) ::  pixloc !wasp
LOGICAL :: do_xbin, do_ybin


! ----------------------------------------------------------

! Boundary wavelengths (approximate) for UV-2 and VIS channels
! ------------------------------------------------------------
REAL (KIND=r4), DIMENSION(mchannel),  PARAMETER :: upper_wvls = (/310.0, 387.0/), &
                                                   lower_wvls = (/260.0, 310.0/)
REAL (KIND=dp), DIMENSION(mchannel)             :: reduce_ubnd, reduce_lbnd, retubnd, retlbnd


! Time varialbees

!-------------------------------
! Gems L1b data
!------------------------------

TYPE gemsl1btype

  INTEGER :: nxtrack, nwavel
  INTEGER (KIND=i2), POINTER, DIMENSION (:)        :: mqflg
  REAL    (KIND=r4), POINTER, DIMENSION (:)        :: ExposureTime
  REAL    (KIND=r8), POINTER, DIMENSION (:)        :: time 
  INTEGER (KIND=i1), POINTER, DIMENSION (:,:)      :: xtrackqflg  !  8bit
  INTEGER (KIND=i2), POINTER, DIMENSION (:,:)      :: height 
  INTEGER (KIND=i2), POINTER, DIMENSION (:,:)      :: groundqflg  ! 16bit
  REAL    (KIND=r4), POINTER, DIMENSION (:,:)      :: lon, lat
  REAL    (KIND=r4), POINTER, DIMENSION (:,:)      :: sza, vza, saz, vaz, aza

  INTEGER (KIND=i2), POINTER, DIMENSION (:,:,:)    :: qflg 
  REAL    (KIND=r4), POINTER, DIMENSION (:,:,:)    :: spec , prec, wavl

END TYPE gemsl1btype
TYPE (gemsl1btype) :: gems_uv1, gems_uv2


!-------------------------------
! Gems processed radiance data
!-------------------------------
 TYPE radtype
    INTEGER (KIND=i2), POINTER, DIMENSION(:,:,:)  :: qflg              ! w, x, y
    REAL    (KIND=r4), POINTER, DIMENSION(:,:,:)  :: spec , prec, wavl ! w, x, y
    INTEGER (KIND=i2), POINTER, DIMENSION(:,:)    :: errstat   ! x, y
    REAL    (KIND=r8), POINTER, DIMENSION(:,:)    :: norm  ! x, y
    INTEGER,           POINTER, DIMENSION(:,:)    :: nwav ! x, y
    INTEGER,           POINTER, DIMENSION(:,:,:)  :: npix ! maxwin, nx, ny
    INTEGER,           POINTER, DIMENSION(:,:,:)  :: wind ! w, nx, ny , the indiced of selected radiance
    INTEGER (KIND=i2), POINTER, DIMENSION(:)      :: line_errstat !y
 END TYPE radtype
 TYPE (radtype) :: gems_rad

!-------------------------------
! Gems processed irradiance data
!-------------------------------
 TYPE irradtype
    INTEGER (KIND=i2), POINTER, DIMENSION (:,:)   :: qflg ! w, x
    REAL    (KIND=r4), POINTER, DIMENSION (:,:)   :: spec, prec, wavl ! w,x
    INTEGER (KIND=i2), POINTER, DIMENSION (:)     :: errstat ! x
    REAL    (KIND=r8), POINTER, DIMENSION (:)     :: norm !x
    INTEGER,           POINTER, DIMENSION (:)     :: nwav !x
    INTEGER,           POINTER, DIMENSION (:,:)   :: npix !maxin, nx
    INTEGER,           POINTER, DIMENSION (:,:)   :: wind !w, nx  ! the indices ofselected irradiance
    REAL    (KIND=r8), POINTER, DIMENSION (:,:,:) :: winpix ! maxwin, x, 2
 END TYPE irradtype
 TYPE (irradtype) :: gems_irrad


!-------------------------------
! Gems processed Geocation data, nx X ny 
!-------------------------------
INTEGER (KIND=i2),     POINTER, DIMENSION (:)     :: gems_mqflg ! ny
INTEGER (KIND=i2),     POINTER, DIMENSION (:)     :: gems_saa !ny
REAL    (KIND=r8),     POINTER, DIMENSION (:)     :: gems_time  ! ny
INTEGER (KIND=i1),     POINTER, DIMENSION (:,:)   :: gems_xtrackQflg  !  8bit
INTEGER (KIND=i2),     POINTER, DIMENSION (:,:)   :: gems_height ! 16bit integer
INTEGER (KIND=i2),     POINTER, DIMENSION (:,:)   :: gems_groundQflg  ! 16bit
REAL    (KIND=r8),     POINTER, DIMENSION (:,:)   :: gems_lon, gems_lat
REAL    (KIND=r8),     POINTER, DIMENSION (:,:)   :: gems_sza, gems_vza, gems_aza, gems_sca
REAL    (KIND=r8),     POINTER, DIMENSION (:,:)   :: gems_clat,gems_clon
REAL    (KIND=r8),     POINTER, DIMENSION (:,:)   :: gems_elat, gems_elon

! + Converted Flags From XtrackQflag
TYPE xflagtype
 INTEGER (KIND=i1),     POINTER, DIMENSION (:,:)  :: rowanomaly, waveshift
 INTEGER (KIND=i1),     POINTER, DIMENSION (:,:)  :: blockage, straysun, strayearth
END TYPE xflagtype
TYPE (xflagtype) :: gems_xflag

! + Converted Flags From Ground Pixel Quality Flags
TYPE groundflagtype
  INTEGER (KIND=i2),   POINTER, DIMENSION (:,:)   ::  geo
  INTEGER (KIND=i2),   POINTER, DIMENSION (:,:)   ::  land_water
  INTEGER (KIND=i2),   POINTER, DIMENSION (:,:)   ::  glint
  INTEGER (KIND=i2),   POINTER, DIMENSION (:,:)   ::  snow_ice
END TYPE groundflagtype
TYPE (groundflagtype) :: gems_gflag

!------------------------
! cloud variables
!-------------------------
TYPE CloudBlock
  REAL (KIND=r4),   POINTER, DIMENSION (:,:)      ::  cfr, ctp
  INTEGER(KIND=i1), POINTER, DIMENSION (:,:)      ::  qflags 
END TYPE CloudBLock
TYPE (cloudBlock):: gems_clouds

!------------------------
! Array for ring spectra
!-------------------------
TYPE ringtype 
    INTEGER,        POINTER, DIMENSION (:)       :: nsol              !nxtrack_max, N of ring spectrum
    INTEGER,        POINTER, DIMENSION (:)       :: sol_lin, sol_uin  !nxtrack_max, Lower/upper index of ring spectrum
    INTEGER,        POINTER, DIMENSION (:)       :: sol_ndiv          !nxtrack_max
    REAL (KIND=r4), POINTER, DIMENSION (:,: )    :: solspec, solwavl  !max_ring_pts,nxtrack_max 
END TYPE ringtype
TYPE (ringtype) gems_ring

!-----------------------
! Array for reflectance
!-----------------------
 TYPE reflectance
     REAL (KIND=r4), POINTER, DIMENSION (:,: )      ::  solspec    !mrefl, nxtrack_max
     REAL (KIND=r4), POINTER, DIMENSION (:,: )      ::  solwavl    !mrefl, nxtrack_max
     REAL (KIND=r4), POINTER, DIMENSION (:,:,:)     ::  radspec    !mrefl, nxtrack_max,ntimes_max
     REAL (KINd=r4), POINTER, DIMENSION (:,:,:)     ::  radwavl    !mrefl, nxtrack_max,ntimes_max
     REAL (KIND=r4), POINTER, DIMENSION (:,:)       ::  solwinpix  !nxtrack_max, 2
 END TYPE reflectance
 TYPE (reflectance) ::  gems_refl


!------------------------
! calibration variables
!-------------------------
 REAL (KIND=dp), DIMENSION (maxwin, nxtrack_max)                         :: gems_wincal_wav

 REAL (KIND=dp), DIMENSION (maxwin)                 :: gems_redslw
 TYPE solcaltype
      INTEGER,       POINTER, DIMENSION (:)         :: nslit           !(nxtrack_max)
      INTEGER,       POINTER, DIMENSION (:)         :: nwavcal         !(nxtrack_max)
      REAL (KIND=dp),POINTER, DIMENSION (:,:,:,:)   :: winfit          !(maxwin,      max_calfit_idx,2, nxtrack_max)
      REAL (KIND=dp),POINTER, DIMENSION (:,:)       :: slitwav, sswav  !(max_fit_pts, nxtrack_max)  
      REAL (KIND=dp),POINTER, DIMENSION (:,:,:,:)   :: slitfit, wavfit !(max_fit_pts, max_calfit_idx, 2, nxtrack_max)
 END TYPE solcaltype
 TYPE (solcaltype) :: gems_solcal

 TYPE radcaltype
      INTEGER,       POINTER, DIMENSION (:)         :: nslit           !(nxtrack_max)
      INTEGER,       POINTER, DIMENSION (:)         :: nwavcal         !(nxtrack_max)
      REAL (KIND=dp),POINTER, DIMENSION (:,:,:,:)   :: winfit          !(maxwin,      max_calfit_idx,2, nxtrack_max)
      REAL (KIND=dp),POINTER, DIMENSION (:,:)       :: slitwav, sswav  !(max_fit_pts, nxtrack_max)  
      REAL (KIND=dp),POINTER, DIMENSION (:,:,:,:)   :: slitfit, wavfit !(max_fit_pts, max_calfit_idx, 2, nxtrack_max)
 END TYPE radcaltype
 TYPE (radcaltype) :: gems_radcal                                                                          

!------------------------
! gems out
!-------------------------
INTEGER,             POINTER, DIMENSION (:,:)       :: gems_exitval, gems_initval !(nxtrack_max, ntimes_max)
REAL (KIND=dp),      POINTER, DIMENSION (:,:,:)     :: gems_fitvar ! (nxtrack_max, ntimes_max,  n_max_fitpars)


CONTAINS 

SUBROUTINE  allocate_o3p_var (ny,  pge_error_status)
 USE OMSAO_errstat_module  

 IMPLICIT NONE
! ----------------
! Output variables
! ----------------
INTEGER, INTENT(IN)   :: ny
INTEGER, INTENT (OUT) :: pge_error_status
! ------------------------------
!  Local variables
! ------------------------------
  INTEGER(KIND=4)         :: status
  INTEGER(KIND=4)         :: nx, nw

 pge_error_status = pge_errstat_ok
 nx  = nxtrack_max
 nw  = nwavel_max


 !-------------------------- GEMS_RAD structure ------------------------------------
 IF (ASSOCIATED(GEMS_Rad%qflg)) DEALLOCATE(GEMS_Rad%qflg)
 ALLOCATE (GEMS_Rad%qflg(nw, nx, ny), STAT = status) ; IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED(GEMS_Rad%wavl)) DEALLOCATE(GEMS_Rad%wavl)
 ALLOCATE (GEMS_Rad%wavl(nw, nx, ny), STAT = status) ; IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED(GEMS_Rad%spec)) DEALLOCATE(GEMS_Rad%spec)
 ALLOCATE (GEMS_Rad%spec(nw, nx, ny), STAT = status) ; IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED(GEMS_Rad%prec)) DEALLOCATE(GEMS_Rad%prec)
 ALLOCATE (GEMS_Rad%prec(nw, nx, ny), STAT = status) ; IF (status.ne. 0)  GOTO 111

 IF (ASSOCIATED(GEMS_Rad%errstat)) DEALLOCATE(GEMS_Rad%errstat)
 ALLOCATE (GEMS_Rad%errstat (nx, ny), STAT = status) ; IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED(GEMS_Rad%norm)) DEALLOCATE(GEMS_Rad%norm)
 ALLOCATE (GEMS_Rad%norm (nx, ny), STAT = status) ; IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED(GEMS_Rad%line_errstat)) DEALLOCATE(GEMS_Rad%line_errstat)
 ALLOCATE (GEMS_Rad%line_errstat (ny), STAT = status) ; IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED(GEMS_Rad%nwav)) DEALLOCATE(GEMS_Rad%nwav)
 ALLOCATE (GEMS_Rad%nwav(nx, ny), STAT = status) ; IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED(GEMS_Rad%npix)) DEALLOCATE(GEMS_Rad%npix)
 ALLOCATE (GEMS_Rad%npix(maxwin, nx, ny), STAT = status) ; IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED(GEMS_Rad%wind)) DEALLOCATE(GEMS_Rad%wind)
 ALLOCATE (GEMS_Rad%wind (nw, nx, ny), STAT = status) ; IF (status .ne. 0)  GOTO 111

 !-------------------------- GEMS_IRRAD structure ------------------------------------
 IF (ASSOCIATED(GEMS_irrad%qflg)) DEALLOCATE(GEMS_irrad%qflg)
 ALLOCATE (GEMS_irrad%qflg(nw, nx), STAT = status) ; IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED(GEMS_irrad%wavl)) DEALLOCATE(GEMS_irrad%wavl)
 ALLOCATE (GEMS_irrad%wavl(nw, nx), STAT = status) ; IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED(GEMS_irrad%spec)) DEALLOCATE(GEMS_irrad%spec)
 ALLOCATE (GEMS_irrad%spec(nw, nx), STAT = status) ; IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED(GEMS_irrad%prec)) DEALLOCATE(GEMS_irrad%prec)
 ALLOCATE (GEMS_irrad%prec(nw, nx), STAT = status) ; IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED(GEMS_irrad%errstat)) DEALLOCATE(GEMS_irrad%errstat)
 ALLOCATE (GEMS_irrad%errstat(nx), STAT = status) ; IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED(GEMS_irrad%norm)) DEALLOCATE(GEMS_irrad%norm)
 ALLOCATE (GEMS_irrad%norm(nx), STAT = status) ; IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED(GEMS_irrad%nwav)) DEALLOCATE(GEMS_irrad%nwav)
 ALLOCATE (GEMS_irrad%nwav(nx), STAT = status) ; IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED(GEMS_irrad%npix)) DEALLOCATE(GEMS_irrad%npix)
 ALLOCATE (GEMS_irrad%npix(maxwin,nx), STAT = status) ; IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED(GEMS_irrad%wind)) DEALLOCATE(GEMS_irrad%wind)
 ALLOCATE (GEMS_irrad%wind(nw, nx), STAT = status) ; IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED(GEMS_irrad%winpix)) DEALLOCATE(GEMS_irrad%winpix)
 ALLOCATE (GEMS_irrad%winpix(maxwin,nx,2), STAT = status) ; IF (status .ne. 0)  GOTO 111


 !-------------------------- GEMS processed GEOLocation data ------------------------------------

 IF (ASSOCIATED(gems_mqflg))       DEALLOCATE(gems_mqflg)
 ALLOCATE (gems_mqflg(ny)         ,STAT = status) ; IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED(gems_saa  ))       DEALLOCATE(gems_saa)
 ALLOCATE (gems_saa(ny)           ,STAT = status) ; IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED(gems_time ))       DEALLOCATE(gems_time)
 ALLOCATE (gems_time(ny)          ,STAT = status) ; IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED(gems_xtrackQflg  ))DEALLOCATE(gems_xtrackQflg)
 ALLOCATE (gems_xtrackQflg(nx, ny),STAT = status) ; IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED(gems_height  ))    DEALLOCATE(gems_height)
 ALLOCATE (gems_height(nx, ny)    ,STAT = status) ; IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED(gems_groundQflg  ))DEALLOCATE(gems_groundQflg)
 ALLOCATE (gems_groundQflg(nx, ny),STAT = status) ; IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED(gems_lon  ))       DEALLOCATE(gems_lon)
 ALLOCATE (gems_lon(nx, ny),STAT = status) ; IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED(gems_lat  ))       DEALLOCATE(gems_lat)
 ALLOCATE (gems_lat(nx, ny),STAT = status) ; IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED( gems_sza  ))       DEALLOCATE( gems_sza)
 ALLOCATE (gems_sza(nx, ny),STAT = status) ; IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED( gems_vza  ))       DEALLOCATE( gems_vza)
 ALLOCATE (gems_vza(nx, ny),STAT = status) ; IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED( gems_aza  ))       DEALLOCATE( gems_aza)
 ALLOCATE (gems_aza(nx, ny),STAT = status) ; IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED( gems_sca  ))       DEALLOCATE( gems_sca)
 ALLOCATE (gems_sca(nx, ny),STAT = status) ; IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED( gems_clon  ))      DEALLOCATE( gems_clon)
 ALLOCATE (gems_clon(0:nx, ny+1),STAT = status); IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED( gems_clat  ))      DEALLOCATE( gems_clat)
 ALLOCATE (gems_clat(0:nx, ny+1),STAT = status); IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED( gems_elon  ))      DEALLOCATE( gems_elon)
 ALLOCATE (gems_elon(0:nx, ny),STAT = status); IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED( gems_elat  ))      DEALLOCATE( gems_elat)
 ALLOCATE (gems_elat(0:nx, ny),STAT = status); IF (status .ne. 0)  GOTO 111

 ! + Converted Flags From XtrackQflag
 IF (ASSOCIATED( gems_xflag%rowanomaly  ))      DEALLOCATE( gems_xflag%rowanomaly)
 ALLOCATE (gems_xflag%rowanomaly(nx, ny),STAT = status); IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED( gems_xflag%waveshift  ))      DEALLOCATE( gems_xflag%waveshift)
 ALLOCATE (gems_xflag%waveshift(nx, ny),STAT = status); IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED( gems_xflag%blockage  ))      DEALLOCATE( gems_xflag%blockage)
 ALLOCATE (gems_xflag%blockage(nx, ny),STAT = status); IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED( gems_xflag%straysun  ))      DEALLOCATE( gems_xflag% straysun)
 ALLOCATE (gems_xflag% straysun(nx, ny),STAT = status); IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED( gems_xflag%strayearth  ))      DEALLOCATE( gems_xflag%strayearth)
 ALLOCATE (gems_xflag%strayearth(nx, ny),STAT = status); IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED( gems_gflag%geo  ))      DEALLOCATE( gems_gflag%geo)
 ALLOCATE (gems_gflag%geo(nx, ny),STAT = status); IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED( gems_gflag%land_water  ))      DEALLOCATE( gems_gflag%land_water)
 ALLOCATE (gems_gflag%land_water(nx, ny),STAT = status); IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED( gems_gflag%glint  ))      DEALLOCATE( gems_gflag%glint)
 ALLOCATE (gems_gflag%glint(nx, ny),STAT = status); IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED( gems_gflag%snow_ice  ))      DEALLOCATE( gems_gflag%snow_ice)
 ALLOCATE (gems_gflag%snow_ice(nx, ny),STAT = status); IF (status .ne. 0)  GOTO 111


 !-------------------------- GEMS cloud data ------------------------------------
 IF (ASSOCIATED( gems_clouds%cfr  ))      DEALLOCATE( gems_clouds%cfr)
 ALLOCATE (gems_clouds%cfr(nx, ny),STAT = status); IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED( gems_clouds%ctp  ))      DEALLOCATE( gems_clouds%ctp)
 ALLOCATE (gems_clouds%ctp(nx, ny),STAT = status); IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED( gems_clouds%qflags   ))      DEALLOCATE( gems_clouds%qflags)
 ALLOCATE (gems_clouds%qflags(nx, ny),STAT = status); IF (status .ne. 0)  GOTO 111

 !-------------------------- ring spectrum  ------------------------------------
 IF (ASSOCIATED( gems_ring%nsol   ))      DEALLOCATE( gems_ring%nsol )
 ALLOCATE (gems_ring%nsol (nx), STAT = status); IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED( gems_ring%sol_lin   ))      DEALLOCATE( gems_ring%sol_lin )
 ALLOCATE (gems_ring%sol_lin (nx), STAT = status); IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED( gems_ring%sol_uin   ))      DEALLOCATE( gems_ring%sol_uin )
 ALLOCATE (gems_ring%sol_uin (nx), STAT = status); IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED( gems_ring%sol_ndiv   ))      DEALLOCATE( gems_ring%sol_ndiv )
 ALLOCATE (gems_ring%sol_ndiv (nx), STAT = status); IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED( gems_ring%solspec   ))      DEALLOCATE( gems_ring%solspec )
 ALLOCATE (gems_ring%solspec (max_ring_pts,nx), STAT = status); IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED( gems_ring%solwavl  ))      DEALLOCATE( gems_ring%solwavl )
 ALLOCATE (gems_ring%solwavl  (max_ring_pts,nx), STAT = status); IF (status .ne. 0)  GOTO 111

 !-------------------------- reflectance spectrum  ------------------------------------
 IF (ASSOCIATED( gems_refl%solwavl  ))      DEALLOCATE( gems_refl%solwavl )
 ALLOCATE (gems_refl%solwavl  (mrefl,nx),   STAT = status); IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED( gems_refl%solspec  ))      DEALLOCATE( gems_refl%solspec )
 ALLOCATE (gems_refl%solspec  (mrefl,nx),   STAT = status); IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED( gems_refl%radwavl  ))      DEALLOCATE( gems_refl%radwavl )
 ALLOCATE (gems_refl%radwavl  (mrefl,nx,ny),STAT = status); IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED( gems_refl%radspec  ))      DEALLOCATE( gems_refl%radspec )
 ALLOCATE (gems_refl%radspec  (mrefl,nx,ny),STAT = status); IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED( gems_refl%solwinpix  ))      DEALLOCATE( gems_refl%solwinpix )
 ALLOCATE (gems_refl%solwinpix (nx, 2),   STAT = status); IF (status .ne. 0)  GOTO 111

 !-------------------------- irradiance calibration data  ------------------------------------
 
 IF (ASSOCIATED(gems_solcal%nslit  ))      DEALLOCATE( gems_solcal%nslit )
 ALLOCATE (gems_solcal%nslit (nx),   STAT = status); IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED(gems_solcal%nwavcal  ))      DEALLOCATE( gems_solcal%nwavcal )
 ALLOCATE (gems_solcal%nwavcal (nx),   STAT = status); IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED(gems_solcal%winfit  ))      DEALLOCATE( gems_solcal%winfit )
 ALLOCATE (gems_solcal%winfit (maxwin, max_calfit_idx, 2, nx),   STAT = status); IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED(gems_solcal%slitwav  ))      DEALLOCATE( gems_solcal%slitwav )
 ALLOCATE (gems_solcal%slitwav (max_fit_pts, nx),   STAT = status); IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED(gems_solcal%sswav  ))      DEALLOCATE( gems_solcal%sswav )
 ALLOCATE (gems_solcal%sswav (max_fit_pts, nx),   STAT = status); IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED(gems_solcal%slitfit  ))      DEALLOCATE( gems_solcal%slitfit )
 ALLOCATE (gems_solcal%slitfit (max_fit_pts, max_calfit_idx, 2, nx),   STAT = status); IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED(gems_solcal%wavfit  ))      DEALLOCATE( gems_solcal%wavfit )
 ALLOCATE (gems_solcal%wavfit (max_fit_pts, max_calfit_idx, 2, nx),   STAT = status); IF (status .ne. 0)  GOTO 111


 !-------------------------- gems fitting variables  ------------------------------------
 IF (ASSOCIATED( gems_exitval  ))     DEALLOCATE( gems_exitval )
 ALLOCATE (gems_exitval(nx, ny),      STAT = status); IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED( gems_initval  ))     DEALLOCATE( gems_initval )
 ALLOCATE (gems_initval(nx, ny),      STAT = status); IF (status .ne. 0)  GOTO 111

 IF (ASSOCIATED( gems_fitvar  ))     DEALLOCATE( gems_fitvar )
 ALLOCATE (gems_fitvar(nx, ny, n_max_fitpars),      STAT = status); IF (status .ne. 0)  GOTO 111


RETURN
111 CONTINUE
    WRITE(*,'(A)') 'Error in '//'Error in gems_o3P_allocation'
    pge_error_status = pge_errstat_error
    RETURN

END SUBROUTINE allocate_o3p_var


SUBROUTINE allocate_o3p_l1b (gems, pge_error_status)
 USE OMSAO_errstat_module  
 IMPLICIT NONE
! ----------------
! Output variables
! ----------------
TYPE (gemsl1btype), INTENT(INOUT)    :: gems
INTEGER, INTENT (OUT) :: pge_error_status

 pge_error_status = pge_errstat_ok
 !-----------------  measurement quality flags  ---------------
 IF (ASSOCIATED(gems%mqflg)) DEALLOCATE(gems%mqflg)
 ALLOCATE (gems%mqflg(ntimes_max), STAT = pge_error_status)
 IF(pge_error_status .NE. 0) GOTO 111

 !-----------------  Exposuretime  ---------------
 IF (ASSOCIATED(gems%ExposureTime)) DEALLOCATE(gems%ExposureTime)
 ALLOCATE (gems%ExposureTime(ntimes_max), STAT = pge_error_status)
 IF(pge_error_status .NE. 0) GOTO 111

 !-----------------  Time  ---------------
 IF (ASSOCIATED(gems%Time)) DEALLOCATE(gems%Time)
 ALLOCATE (gems%Time(ntimes_max), STAT = pge_error_status)
 IF(pge_error_status .NE. 0) GOTO 111

 !-----------------  Lat  ---------------
 IF (ASSOCIATED(gems%Lat)) DEALLOCATE(gems%Lat)
 ALLOCATE (gems%Lat(nxtrack_max, ntimes_max), STAT = pge_error_status)
 IF(pge_error_status .NE. 0) GOTO 111

 !-----------------  Lon  ---------------
 IF (ASSOCIATED(gems%Lon)) DEALLOCATE(gems%Lon)
 ALLOCATE (gems%Lon(nxtrack_max, ntimes_max), STAT = pge_error_status)
 IF(pge_error_status .NE. 0) GOTO 111

 !-----------------  height  ---------------
 IF (ASSOCIATED(gems%height)) DEALLOCATE(gems%height)
 ALLOCATE (gems%height(nxtrack_max, ntimes_max), STAT = pge_error_status)
 IF(pge_error_status .NE. 0) GOTO 111

  !-----------------  groundqflg  ---------------
 IF (ASSOCIATED(gems%groundqflg)) DEALLOCATE(gems%groundqflg)
 ALLOCATE (gems%groundqflg(nxtrack_max, ntimes_max), STAT = pge_error_status)
 IF(pge_error_status .NE. 0) GOTO 111

  !----------------- xtrackqflg  ---------------
 IF (ASSOCIATED(gems%xtrackqflg)) DEALLOCATE(gems%xtrackqflg)
 ALLOCATE (gems%xtrackqflg(nxtrack_max, ntimes_max), STAT = pge_error_status)
 IF(pge_error_status .NE. 0) GOTO 111

  !-----------------  sza  ---------------
 IF (ASSOCIATED(gems%sza)) DEALLOCATE(gems%sza)
 ALLOCATE (gems%sza(nxtrack_max, ntimes_max), STAT = pge_error_status)
 IF(pge_error_status .NE. 0) GOTO 111

  !-----------------  vza  ---------------
 IF (ASSOCIATED(gems%vza)) DEALLOCATE(gems%vza)
 ALLOCATE (gems%vza(nxtrack_max, ntimes_max), STAT = pge_error_status)
 IF(pge_error_status .NE. 0) GOTO 111

  !-----------------  saz  ---------------
 IF (ASSOCIATED(gems%saz)) DEALLOCATE(gems%saz)
 ALLOCATE (gems%saz(nxtrack_max, ntimes_max), STAT = pge_error_status)
 IF(pge_error_status .NE. 0) GOTO 111

  !-----------------  vaz ---------------
 IF (ASSOCIATED(gems%vaz)) DEALLOCATE(gems%vaz)
 ALLOCATE (gems%vaz(nxtrack_max, ntimes_max), STAT = pge_error_status)
 IF(pge_error_status .NE. 0) GOTO 111

 !-----------------  aza ---------------
 IF (ASSOCIATED(gems%aza)) DEALLOCATE(gems%aza)
 ALLOCATE (gems%aza(nxtrack_max, ntimes_max), STAT = pge_error_status)
 IF(pge_error_status .NE. 0) GOTO 111


  !-----------------  spec  ---------------
 IF (ASSOCIATED(gems%spec)) DEALLOCATE(gems%spec)
 ALLOCATE (gems%spec(nwavel_max,nxtrack_max, ntimes_max), STAT = pge_error_status)
 IF(pge_error_status .NE. 0) GOTO 111

  !-----------------  prec  ---------------
 IF (ASSOCIATED(gems%prec)) DEALLOCATE(gems%prec)
 ALLOCATE (gems%prec(nwavel_max,nxtrack_max, ntimes_max), STAT = pge_error_status)
 IF(pge_error_status .NE. 0) GOTO 111

  !-----------------  wavl ---------------
 IF (ASSOCIATED(gems%wavl)) DEALLOCATE(gems%wavl)
 ALLOCATE (gems%wavl(nwavel_max,nxtrack_max, ntimes_max), STAT = pge_error_status)
 IF(pge_error_status .NE. 0) GOTO 111

 !-----------------  qflg ---------------
 IF (ASSOCIATED(gems%qflg)) DEALLOCATE(gems%qflg)
 ALLOCATE (gems%qflg(nwavel_max,nxtrack_max, ntimes_max), STAT = pge_error_status)
 IF(pge_error_status .NE. 0) GOTO 111

RETURN

111 CONTINUE
    WRITE(*,'(A)') 'Error in '//'Error in allocate_o3p_l1b'
    pge_error_status = pge_errstat_error
    RETURN
RETURN
END SUBROUTINE allocate_o3p_l1b


SUBROUTINE deallocate_o3p_l1b (gems, pge_error_status)
 USE OMSAO_errstat_module  
 IMPLICIT NONE
! ----------------
! Output variables
! ----------------
TYPE (gemsl1btype), INTENT(INOUT)    :: gems
INTEGER, INTENT (OUT) :: pge_error_status

 pge_error_status = pge_errstat_ok

 !-----------------  measurement quality flags  ---------------
 IF (ASSOCIATED(gems%mqflg)) DEALLOCATE(gems%mqflg)
 !-----------------  Exposuretime  ---------------
 IF (ASSOCIATED(gems%ExposureTime)) DEALLOCATE(gems%ExposureTime)
 !-----------------  Time  ---------------
 IF (ASSOCIATED(gems%Time)) DEALLOCATE(gems%Time)
 !-----------------  Lat  ---------------
 IF (ASSOCIATED(gems%Lat)) DEALLOCATE(gems%Lat)
 !-----------------  Lon  ---------------
 IF (ASSOCIATED(gems%Lon)) DEALLOCATE(gems%Lon)
 !-----------------  height  ---------------
 IF (ASSOCIATED(gems%height)) DEALLOCATE(gems%height)
  !-----------------  groundqflg  ---------------
 IF (ASSOCIATED(gems%groundqflg)) DEALLOCATE(gems%groundqflg)
  !----------------- xtrackqflg  ---------------
 IF (ASSOCIATED(gems%xtrackqflg)) DEALLOCATE(gems%xtrackqflg)
  !-----------------  sza  ---------------
 IF (ASSOCIATED(gems%sza)) DEALLOCATE(gems%sza)
  !-----------------  vza  ---------------
 IF (ASSOCIATED(gems%vza)) DEALLOCATE(gems%vza)
  !-----------------  saz  ---------------
 IF (ASSOCIATED(gems%saz)) DEALLOCATE(gems%saz)
  !-----------------  vaz ---------------
 IF (ASSOCIATED(gems%vaz)) DEALLOCATE(gems%vaz)
 !-----------------  aza ---------------
 IF (ASSOCIATED(gems%aza)) DEALLOCATE(gems%aza)
  !-----------------  spec  ---------------
 IF (ASSOCIATED(gems%spec)) DEALLOCATE(gems%spec)
  !-----------------  prec  ---------------
 IF (ASSOCIATED(gems%prec)) DEALLOCATE(gems%prec)
  !-----------------  wavl ---------------
 IF (ASSOCIATED(gems%wavl)) DEALLOCATE(gems%wavl)
 !-----------------  qflg ---------------
 IF (ASSOCIATED(gems%qflg)) DEALLOCATE(gems%qflg)

RETURN

111 CONTINUE
    WRITE(*,'(A)') 'Error in '//'Error in allocate_o3p_l1b'
    pge_error_status = pge_errstat_error
    RETURN
RETURN
END SUBROUTINE deallocate_o3p_l1b


SUBROUTINE deallocate_o3p_var ( pge_error_status)
 USE OMSAO_errstat_module  
 IMPLICIT NONE
! ----------------
! Output variables
! ----------------
INTEGER, INTENT (OUT) :: pge_error_status

pge_error_status = pge_errstat_ok

 !-------------------------- GEMS_RAD structure ------------------------------------
 IF (ASSOCIATED(GEMS_Rad%qflg))     DEALLOCATE(GEMS_Rad%qflg)
 IF (ASSOCIATED(GEMS_Rad%wavl))     DEALLOCATE(GEMS_Rad%wavl)
 IF (ASSOCIATED(GEMS_Rad%spec))     DEALLOCATE(GEMS_Rad%spec)
 IF (ASSOCIATED(GEMS_Rad%prec))     DEALLOCATE(GEMS_Rad%prec)
 IF (ASSOCIATED(GEMS_Rad%errstat))  DEALLOCATE(GEMS_Rad%errstat)
 IF (ASSOCIATED(GEMS_Rad%norm))     DEALLOCATE(GEMS_Rad%norm)
 IF (ASSOCIATED(GEMS_Rad%line_errstat)) DEALLOCATE(GEMS_Rad%line_errstat)
 IF (ASSOCIATED(GEMS_Rad%nwav))     DEALLOCATE(GEMS_Rad%nwav)
 IF (ASSOCIATED(GEMS_Rad%npix))     DEALLOCATE(GEMS_Rad%npix)
 IF (ASSOCIATED(GEMS_Rad%wind))     DEALLOCATE(GEMS_Rad%wind)

 !-------------------------- GEMS_IRRAD structure ------------------------------------
 IF (ASSOCIATED(GEMS_irrad%qflg))   DEALLOCATE(GEMS_irrad%qflg)
 IF (ASSOCIATED(GEMS_irrad%wavl))   DEALLOCATE(GEMS_irrad%wavl)
 IF (ASSOCIATED(GEMS_irrad%spec))   DEALLOCATE(GEMS_irrad%spec)
 IF (ASSOCIATED(GEMS_irrad%prec))   DEALLOCATE(GEMS_irrad%prec)
 IF (ASSOCIATED(GEMS_irrad%errstat)) DEALLOCATE(GEMS_irrad%errstat)
 IF (ASSOCIATED(GEMS_irrad%norm))   DEALLOCATE(GEMS_irrad%norm)
 IF (ASSOCIATED(GEMS_irrad%nwav))   DEALLOCATE(GEMS_irrad%nwav)
 IF (ASSOCIATED(GEMS_irrad%npix))   DEALLOCATE(GEMS_irrad%npix)
 IF (ASSOCIATED(GEMS_irrad%wind))   DEALLOCATE(GEMS_irrad%wind)
 IF (ASSOCIATED(GEMS_irrad%winpix)) DEALLOCATE(GEMS_irrad%winpix)
 !-------------------------- GEMS processed GEOLocation data ------------------------------------
 IF (ASSOCIATED(gems_mqflg))        DEALLOCATE(gems_mqflg)
 IF (ASSOCIATED(gems_saa  ))        DEALLOCATE(gems_saa)
 IF (ASSOCIATED(gems_time ))        DEALLOCATE(gems_time)
 IF (ASSOCIATED(gems_xtrackQflg  )) DEALLOCATE(gems_xtrackQflg)
 IF (ASSOCIATED(gems_height  ))     DEALLOCATE(gems_height)
 IF (ASSOCIATED(gems_groundQflg  )) DEALLOCATE(gems_groundQflg)
 IF (ASSOCIATED(gems_lon  ))        DEALLOCATE(gems_lon)
 IF (ASSOCIATED(gems_lat  ))        DEALLOCATE(gems_lat)
 IF (ASSOCIATED( gems_sza  ))       DEALLOCATE( gems_sza)
 IF (ASSOCIATED( gems_vza  ))       DEALLOCATE( gems_vza)
 IF (ASSOCIATED( gems_aza  ))       DEALLOCATE( gems_aza)
 IF (ASSOCIATED( gems_sca  ))       DEALLOCATE( gems_sca)
 IF (ASSOCIATED( gems_clon  ))      DEALLOCATE( gems_clon)
 IF (ASSOCIATED( gems_clat  ))      DEALLOCATE( gems_clat)
 IF (ASSOCIATED( gems_elon  ))      DEALLOCATE( gems_elon)
 IF (ASSOCIATED( gems_elat  ))      DEALLOCATE( gems_elat)
 ! + Converted Flags From XtrackQflag
 IF (ASSOCIATED( gems_xflag%rowanomaly  ))  DEALLOCATE( gems_xflag%rowanomaly)
 IF (ASSOCIATED( gems_xflag%waveshift  ))   DEALLOCATE( gems_xflag%waveshift)
 IF (ASSOCIATED( gems_xflag%blockage  ))    DEALLOCATE( gems_xflag%blockage)
 IF (ASSOCIATED( gems_xflag%straysun  ))    DEALLOCATE( gems_xflag% straysun)
 IF (ASSOCIATED( gems_xflag%strayearth  ))  DEALLOCATE( gems_xflag%strayearth)
 IF (ASSOCIATED( gems_gflag%geo  ))         DEALLOCATE( gems_gflag%geo)
 IF (ASSOCIATED( gems_gflag%land_water  ))  DEALLOCATE( gems_gflag%land_water)
 IF (ASSOCIATED( gems_gflag%glint  ))       DEALLOCATE( gems_gflag%glint)
 IF (ASSOCIATED( gems_gflag%snow_ice  ))    DEALLOCATE( gems_gflag%snow_ice)
 !-------------------------- GEMS cloud data ------------------------------------
 IF (ASSOCIATED( gems_clouds%cfr  ))        DEALLOCATE( gems_clouds%cfr)
 IF (ASSOCIATED( gems_clouds%ctp  ))        DEALLOCATE( gems_clouds%ctp)
 IF (ASSOCIATED( gems_clouds%qflags   ))    DEALLOCATE( gems_clouds%qflags)
 !-------------------------- ring spectrum  ------------------------------------
 IF (ASSOCIATED( gems_ring%nsol   ))        DEALLOCATE( gems_ring%nsol )
 IF (ASSOCIATED( gems_ring%sol_lin   ))     DEALLOCATE( gems_ring%sol_lin )
 IF (ASSOCIATED( gems_ring%sol_uin   ))     DEALLOCATE( gems_ring%sol_uin )
 IF (ASSOCIATED( gems_ring%sol_ndiv   ))    DEALLOCATE( gems_ring%sol_ndiv )
 IF (ASSOCIATED( gems_ring%solspec   ))     DEALLOCATE( gems_ring%solspec )
 IF (ASSOCIATED( gems_ring%solwavl  ))      DEALLOCATE( gems_ring%solwavl )
 !-------------------------- reflectance spectrum  ------------------------------------
 IF (ASSOCIATED( gems_refl%solwavl  ))      DEALLOCATE( gems_refl%solwavl )
 IF (ASSOCIATED( gems_refl%solspec  ))      DEALLOCATE( gems_refl%solspec )
 IF (ASSOCIATED( gems_refl%radwavl  ))      DEALLOCATE( gems_refl%radwavl )
 IF (ASSOCIATED( gems_refl%radspec  ))      DEALLOCATE( gems_refl%radspec )
 IF (ASSOCIATED( gems_refl%solwinpix  ))    DEALLOCATE( gems_refl%solwinpix )
 !-------------------------- irradiance calibration data  ------------------------------------
 
 IF (ASSOCIATED(gems_solcal%nslit  ))      DEALLOCATE( gems_solcal%nslit )
 IF (ASSOCIATED(gems_solcal%nwavcal  ))    DEALLOCATE( gems_solcal%nwavcal )
 IF (ASSOCIATED(gems_solcal%winfit  ))     DEALLOCATE( gems_solcal%winfit )
 IF (ASSOCIATED(gems_solcal%slitwav  ))    DEALLOCATE( gems_solcal%slitwav )
 IF (ASSOCIATED(gems_solcal%sswav  ))      DEALLOCATE( gems_solcal%sswav )
 IF (ASSOCIATED(gems_solcal%slitfit  ))    DEALLOCATE( gems_solcal%slitfit )
 IF (ASSOCIATED(gems_solcal%wavfit  ))     DEALLOCATE( gems_solcal%wavfit )

 !-------------------------- gems fitting variables  ------------------------------------
 IF (ASSOCIATED( gems_exitval  ))    DEALLOCATE( gems_exitval )
 IF (ASSOCIATED( gems_initval  ))    DEALLOCATE( gems_initval )
 IF (ASSOCIATED( gems_fitvar  ))     DEALLOCATE( gems_fitvar )

RETURN

111 CONTINUE
    WRITE(*,'(A)') 'Error in '//'Error in allocate_o3p_l1b'
    pge_error_status = pge_errstat_error
    RETURN
RETURN
END SUBROUTINE deallocate_o3p_var

END MODULE GEMS_O3P_gemsdata_module


