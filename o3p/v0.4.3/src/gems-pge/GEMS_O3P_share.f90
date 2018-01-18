SUBROUTINE gems_o3p_share_l1b (pge_error_status)

! Share Module
  USE Share_MOD_Constants_Variables
  USE Share_MOD_Radiance2
  USE Share_MOD_L1B_Read
  USE Share_MOD_Log
! my Module
  USE OMSAO_variables_module,   ONLY: band_selectors, sol_identifier, rad_identifier, & 
                                      outdir,l1b_rad_filename,l1b_irrad_filename, &
                                      wavcal_fname, swavcal_fname,slit_fname,rslit_fname

  USE OMSAO_errstat_module  
  USE GEMS_O3P_gemsdata_module, ONLY: ntimes_max, nxtrack_max, nwavel_max, nxtrack , nfxtrack, ntimes,nchannel, chs, &
                                      orbnum, orbnumsol, orbc, orbcsol, ncoadd, &
                                      gems_uv1, gems_uv2, gemsraddate, &
                                      allocate_o3p_l1b


  IMPLICIT NONE
! ----------------
! Output variables
! ----------------
  INTEGER, INTENT (OUT) :: pge_error_status
! ------------------------------
!  Share Local variables
! ------------------------------
  INTEGER(KIND=4)         :: status
  INTEGER(KIND=4)         :: nArg1, nArg2, nArg3, nArg4
! ------------------------------
!  Local variables
! ------------------------------
  INTEGER :: i, ch, nx_uv1, nx_uv2, nx, ny, nw_uv1, nw_uv2, nw

  CHARACTER(LEN=128)    :: ctl_fpath
! ------------------------------
! Name of this module/subroutine
! ------------------------------
  CHARACTER (LEN=26), PARAMETER :: modulename = 'gems_o3p_allocate_gemsl1b'


!-------------------------
!-- GEMS L1B read module
!-------------------------
 call GEMS_Share_GetRadiance2
!status  = GEMS_Share_MOD_L1Brug_Read(gds_brug, gcc_l1b_ctl_file)
!status = GEMS_Share_MOD_L1Birr_Read(gds_birr, gcc_l1b_ctl_file)


gds_l1b_fpath%irr=gds_l1b_fpath%rug ! L1Birr 를 읽지 않았으므로, irradiance 이름은 빈공간
l1b_rad_filename    = adjustl(trim(gds_l1b_fpath%rug))
l1b_irrad_filename  = adjustl(trim(gds_l1b_fpath%irr))

! obtain orbit number from irradiance file
i = INDEX(l1b_irrad_filename , '-o') + 2
orbcsol = l1b_irrad_filename (i : i + 5)
READ (orbcsol, *) orbnumsol

! obtain orbit number from radiance file
i = INDEX(l1b_rad_filename, '-o') + 2
orbc = l1b_rad_filename(i : i + 5)
READ (orbc, *) orbnum

! Read date from radiance file (used for correcting sun-earth distance when using backupirradiance)
i = INDEX(l1b_rad_filename, '-o') -14
gemsraddate = l1b_rad_filename(i : i + 8)

! generate identifer for irradiance and radiance spectrum
rad_identifier = 'o' // orbc
sol_identifier = 'o' // orbcsol

slit_fname    = TRIM(ADJUSTL(outdir)) // 'slit_o'    !// orbcsol
swavcal_fname = TRIM(ADJUSTL(outdir)) // 'swavcal_o' !// orbcsol
rslit_fname   = TRIM(ADJUSTL(outdir)) // 'rslit_o' // orbc
wavcal_fname  = TRIM(ADJUSTL(outdir)) // 'wavcal_o'// orbc

 pge_error_status = pge_errstat_ok
 ntimes = 0 ; nxtrack = 0


! set fitting parameter
! nfxtrack = the number of pixels for UV-1 if it is selected
ntimes = gci_ntimes
IF (nchannel == 2) THEN 
   chs(1) = 1 ; chs(2) = 2
   nfxtrack = 30  ;   nxtrack = 60
ELSE  IF (nchannel == 1) THEN
   IF (band_selectors (1) == 1 ) THEN
     chs(1) = 1
     nfxtrack = 30 ; nxtrack = 30 ! this case is impossible due to albedo uv2 channel 
   ELSE
     chs(1) = 2
     nfxtrack = 60 ; nxtrack = 60
   ENDIF
ELSE
   WRITE(*, '(A)') 'Check # of channel selected '
   pge_error_status = pge_errstat_error    
ENDIF


IF (ntimes > ntimes_max) THEN
   pge_error_status = pge_errstat_error
   WRITE(*, '(A,I5)') 'Need to increase ntimes_max >= ', ntimes
ENDIF

IF (nxtrack > nxtrack_max) THEN
   pge_error_status = pge_errstat_error
   WRITE(*, '(A)') 'Need to increase nxtrack_max!!!'
ENDIF

!-----------------------------------------------------------
! Allocate gems l1b rad variables from GEMS_SHARE_MOD
!----------------------------------------------------------
nx_uv1 = gci_nxtrack
nx_uv2 = gci_nxtrack2
nw_uv1 = gci_nwavel
nw_uv2 = gci_nwavel2
ny     = ntimes

! OPEN data block for geolocation fields (UV-1, if both are selected) 
   call allocate_o3p_l1b (gems_uv1, pge_error_status)
   call allocate_o3p_l1b (gems_uv2, pge_error_status)



    nw = nw_uv1 ; nx = nx_uv1

    gems_uv1%nxtrack = nx
    gems_uv1%nwavel  = nw
    gems_uv1%time(1:ny)             = gds_brug%euv1%brug_gfld%Time(1:ny)
    gems_uv1%height(1:nx, 1:ny)     = gds_brug%euv1%brug_gfld%TerrHgt(1:nx, 1:ny)
    gems_uv1%lon(1:nx, 1:ny)        = gds_brug%euv1%brug_gfld%LON(1:nx, 1:ny) 
    gems_uv1%lat(1:nx, 1:ny)        = gds_brug%euv1%brug_gfld%LAT(1:nx, 1:ny)
    gems_uv1%sza(1:nx, 1:ny)        = gds_brug%euv1%brug_gfld%SolZenAng(1:nx, 1:ny)
    gems_uv1%vza(1:nx, 1:ny)        = gds_brug%euv1%brug_gfld%ViewZenAng(1:nx, 1:ny)
    gems_uv1%saz(1:nx, 1:ny)        = gds_brug%euv1%brug_gfld%SolAziAng(1:nx, 1:ny)
    gems_uv1%vaz(1:nx, 1:ny)        = gds_brug%euv1%brug_gfld%ViewAziAng(1:nx, 1:ny)
    gems_uv1%groundqflg(1:nx, 1:ny) = gds_brug%euv1%brug_gfld%GPQFlag(1:nx, 1:ny) 
    gems_uv1%xtrackqflg(1:nx,1:ny)  = gds_brug%euv1%brug_gfld%XTQFlag(1:nx,1:ny)
    gems_uv1%mqflg(1:ny)            = gds_brug%euv1%brug_dfld%MQFlag(1:ny)   
    gems_uv1%exposureTime(1:ny)     = gds_brug%euv1%brug_dfld%ExpTime(1:ny)
    gems_uv1%wavl(1:nw, 1:nx, 1:ny) = gvf_euv1_radwave(1:nw, 1:nx, 1:ny)
    gems_uv1%spec(1:nw, 1:nx, 1:ny) = gvf_euv1_rad(1:nw, 1:nx, 1:ny)
    gems_uv1%prec(1:nw, 1:nx, 1:ny)=  gvf_euv1_prec(1:nw, 1:nx, 1:ny)
    gems_uv1%qflg(1:nw, 1:nx, 1:ny) = gds_brug%euv1%brug_dfld%PQFlag(1:nw, 1:nx, 1:ny)

 nw = nw_uv2 ; nx = nx_uv2
    gems_uv2%nxtrack = nx 
    gems_uv2%nwavel  = nw    
    gems_uv2%height(1:nx, 1:ny)     = gds_brug%euv2%brug_gfld%TerrHgt(1:nx, 1:ny)
    gems_uv2%lon(1:nx, 1:ny)        = gds_brug%euv2%brug_gfld%LON(1:nx, 1:ny)
    gems_uv2%lat(1:nx, 1:ny)        = gds_brug%euv2%brug_gfld%LAT(1:nx, 1:ny)  
    gems_uv2%sza(1:nx, 1:ny)        = gds_brug%euv2%brug_gfld%SolZenAng(1:nx, 1:ny)
    gems_uv2%vza(1:nx, 1:ny)        = gds_brug%euv2%brug_gfld%ViewZenAng(1:nx, 1:ny)
    gems_uv2%saz(1:nx, 1:ny)        = gds_brug%euv2%brug_gfld%SolAziAng(1:nx, 1:ny)
    gems_uv2%vaz(1:nx, 1:ny)        = gds_brug%euv2%brug_gfld%ViewAziAng(1:nx, 1:ny)
    gems_uv2%groundqflg(1:nx, 1:ny) = gds_brug%euv2%brug_gfld%GPQFlag(1:nx, 1:ny)
    gems_uv2%xtrackqflg(1:nx,1:ny)  = gds_brug%euv2%brug_gfld%XTQFlag(1:nx,1:ny)
    gems_uv2%mqflg(1:ny)            = gds_brug%euv2%brug_dfld%MQFlag(1:ny)
    gems_uv2%exposureTime(1:ny)     = gds_brug%euv2%brug_dfld%ExpTime(1:ny)
    gems_uv2%wavl(1:nw, 1:nx, 1:ny) = gvf_euv2_radwave(1:nw, 1:nx, 1:ny)
    gems_uv2%spec(1:nw, 1:nx, 1:ny) = gvf_euv2_rad(1:nw, 1:nx, 1:ny)
    gems_uv2%prec(1:nw, 1:nx, 1:ny) = gvf_euv2_prec(1:nw, 1:nx, 1:ny)
    gems_uv2%qflg(1:nw, 1:nx, 1:ny) = gds_brug%euv2%brug_dfld%PQFlag(1:nw, 1:nx, 1:ny)
  
 CALL GEMS_Share_DeallocMem2

RETURN

END SUBROUTINE gems_o3p_share_l1b


SUBROUTINE Gems_o3p_share_l2_cld (nx, ny, gems_cfr, gems_ctp, gems_qflag)
  USE OMSAO_precision_module
  USE Share_l2_cld_mod_write_read
  USE Share_MOD_Log
  USE Share_MOD_Constants_Variables
  IMPLICIT NONE


  ! ================
  ! Input variables
  ! ================
  INTEGER, INTENT (IN) :: nx, ny
  ! ================
  ! Output variables
  ! ================
  REAL (KIND=r4), DIMENSION (nx, ny), INTENT(OUT) ::  gems_cfr
  INTEGER (KIND=i2), DIMENSION (nx, ny), INTENT(OUT) :: gems_ctp
  INTEGER (KIND=i2), DIMENSION (nx, ny), INTENT(OUT) :: gems_qflag
  ! ================
  ! Local variables
  ! ================
   CHARACTER(LEN=128)    :: ctl_fpath
   TYPE(CLD_ds)          :: CLDds
!  TYPE(L2_cld)          :: L2CLD      ! For GEMS L2 CLD
   TYPE(OMI_L2CLD)       :: L2CLD       ! For OMI L2 CLD
   INTEGER :: i 
   INTEGER (KIND = 4)    :: status


   ctl_fpath = '../../../share/conf/lv2readmdl.nml'
   !-------------------------------
   !-- Initialize Global Constants
   !-------------------------------
   !CALL GEMS_Share_Init_GlobalConstants
   CALL GEMS_Share_Init_L2CLD_GlobalConstants

   !-------------------------------
   !-- L2 CLD read module 
   !-------------------------------
!   status = GEMS_Share_MOD_L2CLD_Read(L2CLD, CLDds, ctl_fpath)     ! For GEMS L2 CLD
   status = GEMS_Share_MOD_OMI_L2CLD_Read(L2CLD, CLDds, ctl_fpath)  ! For OMI L2 CLD
      
!   gems_cfr(1:nx,1:ny) = L2CLD%cld_dfld%CldFrac(1:nx,1:ny)         ! For GEMS L2 CLD
!   gems_ctp(1:nx,1:ny) = L2CLD%cld_dfld%CldP(1:nx,1:ny)
!   gems_qflag(1:nx,1:ny)  = L2CLD%cld_dfld%PrsQFlag(1:nx,1:ny)

   gems_cfr(1:nx,1:ny)    = L2CLD%omicld_dfld%CldFrac(1:nx,1:ny)       ! For OMI L2 CLD
   gems_ctp(1:nx,1:ny)    = L2CLD%omicld_dfld%CldP(1:nx,1:ny)
   gems_qflag(1:nx,1:ny)  = L2CLD%omicld_dfld%PrsQFlag(1:nx,1:ny)


  ! Need Delocation 
   RETURN
  END SUBROUTINE Gems_o3p_share_l2_cld

