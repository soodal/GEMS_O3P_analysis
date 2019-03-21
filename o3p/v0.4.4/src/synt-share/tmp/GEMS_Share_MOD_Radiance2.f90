!--------------------------------------------------------------
!+Module to use radiance.
!
MODULE Share_MOD_Radiance2

!-------------------------------------------------------------------------------
!+Description: 
!     Module to use GEMS Radiance Module
!
!+Version Date       Comment
! ------- ---------- -------
! 0.1.0   2015.03.09 Fisrt  Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------

! Declarations:

! Modules used:
   USE Share_MOD_Log
   USE Share_MOD_P_Constants_Variables, ONLY: LOGMSG, llvl
   USE Share_MOD_Constants_Variables
   USE Share_MOD_DynamicMem
   USE Share_MOD_Radiance
   USE Share_MOD_L1B_Read


!**********************************************************!
   IMPLICIT NONE

CONTAINS

!-------------------------------------------------------------------------------
!+Description: 
!      GEMS Share Module Radiance Getting Program.
!
!
! Method:
!
!
! Input :
!       L1BRUG, L1BRVG, L1BIRR file
!
! Output:
!       Radiance, Irradiance, Wavelength value
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2015.04.02 Fisrt  Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------
SUBROUTINE GEMS_Share_GetRadiance22
    print*, 'xxx'
END SUBROUTINE GEMS_Share_GetRadiance22


SUBROUTINE GEMS_Share_GetRadiance

! Local parameters For L1B Reading:
    CHARACTER(LEN=128)      :: ctl_fpath

    !TYPE(L1B_rug)           :: L1Brug       !-- L1BRUG data structure
    !TYPE(L1B_rvg)           :: L1Brvg       !-- L1BRVG data structure
    !TYPE(L1B_irr)           :: L1Birr       !-- L1BIRR data structure

    INTEGER(KIND=4)         :: status

    INTEGER(KIND=4)         :: nArg1, nArg2, nArg3, nArg4

    !--------------------------------
    !---- GEMS_Share_GetRadiance PROGRAM Start ---
    !--------------------------------
    LOGMSG=" PROC MSG"
    llvl = 9                     ! Log Level
    CALL GEMS_Share_MOD_log(llvl, '--- Start GEMS_Share_GetRadiance PROGRAM ---', 'START MSG')
    
    !-------------------------------
    !-- Initialize Global Constants
    !-------------------------------
    CALL GEMS_Share_Init_GlobalConstants
    CALL GEMS_Share_Init_L1BGlobalConstants

    !-------------------------------
    !-- Memory Allocation
    !-------------------------------
    CALL GEMS_Share_AllocMem2       ! Radiance 전역변수의 메모리를 할당

    !/-------------------------------
    !-- Namelist setting
    !-------------------------------
    ctl_fpath = gcc_l1b_ctl_file

    !-------------------------
    !-- GEMS L1B read module
    !-------------------------
    !print*, 'L1BRUG File READ'
    status = GEMS_Share_MOD_L1Brug_Read(gds_brug, ctl_fpath)
    IF ( status .ne. 0 ) GOTO 1999

    !print*, 'L1BRVG File READ'
    status = GEMS_Share_MOD_L1Brvg_Read(gds_brvg, ctl_fpath)
    IF ( status .ne. 0 ) GOTO 1999

    !print*, 'L1BIRR File READ'
    status = GEMS_Share_MOD_L1Birr_Read(gds_birr, ctl_fpath)
    IF ( status .ne. 0 ) GOTO 1999

    !-------------------------
    !-- Initialize For Earth UV1 Wavelength Calculation-1st
    !-------------------------
    nArg1 = gci_nwavel
    nArg2 = gci_nxtrack
    nArg3 = gci_ntimes
    nArg4 = gci_nwavelcoef
    CALL GEMS_Share_init4RadianceCalc(nArg1, nArg2, nArg3, nArg4) ! for using the share module to calculate radiance

    !-------------------------
    !-- Wavelength Calculation
    !-------------------------
    CALL GEMS_Share_Calc_ObservationWavelen(gds_brug%euv1%brug_dfld%WlenRefCol, gds_brug%euv1%brug_dfld%WlenCoef,               &
                                            gvf_euv1_radwave)
    !print*,'gvf_euv1_radwave(1,1,1)=', gvf_euv1_radwave(1,1,1)

    !-------------------------
    !-- Initialize For Earth UV2 and Sun UV2 Wavelength  Calculation-2nd
    !-------------------------
    nArg1 = gci_nwavel2
    nArg2 = gci_nxtrack2
    CALL GEMS_Share_init4RadianceCalc(nArg1, nArg2, nArg3, nArg4) ! for using the share module to calculate radiance

    CALL GEMS_Share_Calc_ObservationWavelen(gds_brug%euv2%brug_dfld%WlenRefCol, gds_brug%euv2%brug_dfld%WlenCoef,               &
                                            gvf_euv2_radwave)
    !print*,'gvf_euv2_radwave(1,1,1)=', gvf_euv2_radwave(1,1,1)

    CALL GEMS_Share_Calc_ObservationWavelen(gds_birr%suv2%birr_dfld%WlenRefCol, gds_birr%suv2%birr_dfld%WlenCoef,               &
                                            gvf_suv2_irrwave)
    !print*, 'gvf_suv2_irrwave  (1,1)=', gvf_suv2_irrwave(1,1)

    !-------------------------
    !-- Initialize For Earth VIS and Sun VIS Wavelength  Calculation-2nd-1
    !-------------------------
    nArg1 = gci_nwavel3
    CALL GEMS_Share_init4RadianceCalc(nArg1, nArg2, nArg3, nArg4) ! for using the share module to calculate radiance

    CALL GEMS_Share_Calc_ObservationWavelen(gds_brvg%brvg_dfld%WlenRefCol, gds_brvg%brvg_dfld%WlenCoef,                         &
                                            gvf_evis_radwave)
    !print*, 'gvf_evis_radwave  (1,1,1)=', gvf_evis_radwave(1,1,1)

    CALL GEMS_Share_Calc_ObservationWavelen(gds_birr%svis%birr_dfld%WlenRefCol, gds_birr%svis%birr_dfld%WlenCoef,               &
                                            gvf_svis_irrwave)
    !print*, 'gvf_svis_irrwave  (1,1)=', gvf_svis_irrwave(1,1)

    !-------------------------
    !-- Initialize For Earth UV1 Radiance Calculation-3rd
    !-------------------------
    nArg1 = gci_nwavel
    nArg2 = gci_nxtrack
    CALL GEMS_Share_init4RadianceCalc(nArg1, nArg2, nArg3, nArg4) ! for using the share module to calculate radiance
    CALL GEMS_Share_allocMem4Rad                                  ! for memory allocation of the share module to calculate radiance

    !-------------------------
    !-- UV1 Radiance Radiance Calculation
    !-------------------------
    CALL GEMS_Share_Calc_RadianceValue(gds_brug%euv1%brug_dfld%RadMantissa,    gds_brug%euv1%brug_dfld%RadExponent,             &
                                       gvf_euv1_rad)
    CALL GEMS_Share_Calc_RadianceValue(gds_brug%euv1%brug_dfld%RadPrecision, gds_brug%euv1%brug_dfld%RadExponent, &
                                       gvf_euv1_prec)
    !print*, 'gvf_euv1_rad  (1,1,1)=', gvf_euv1_rad(1,1,1)
!stop
    !-------------------------
    !-- Initialize For Earth UV2 Radiance Calculation-4th
    !-------------------------
    nArg1 = gci_nwavel2
    nArg2 = gci_nxtrack2
    CALL GEMS_Share_deAllocMem4Rad       ! for memory deallocation in the share module to calculate radiance
    CALL GEMS_Share_init4RadianceCalc(nArg1, nArg2, nArg3, nArg4) ! for using the share module to calculate radiance
    CALL GEMS_Share_allocMem4Rad                                  ! for memory allocation of the share module to calculate radiance

    !-------------------------
    !-- UV2 Radiance Calculation
    !-------------------------
    CALL GEMS_Share_Calc_RadianceValue(gds_brug%euv2%brug_dfld%RadMantissa, gds_brug%euv2%brug_dfld%RadExponent, &
                                       gvf_euv2_rad)
    CALL GEMS_Share_Calc_RadianceValue(gds_brug%euv2%brug_dfld%RadPrecision, gds_brug%euv2%brug_dfld%RadExponent, &
                                       gvf_euv2_prec)
    !print*, 'gvf_euv2_rad  (1,1,1)=', gvf_euv2_rad(1,1,1)

    !-------------------------
    !-- Initialize For Earth VIS Radiance Calculation-5th
    !-------------------------
    nArg1 = gci_nwavel3
    CALL GEMS_Share_deAllocMem4Rad       ! for memory deallocation in the share module to calculate radiance
    CALL GEMS_Share_init4RadianceCalc(nArg1, nArg2, nArg3, nArg4) ! for using the share module to calculate radiance
    CALL GEMS_Share_allocMem4Rad                                  ! for memory allocation of the share module to calculate radiance

    !-------------------------
    !-- Earth VIS Radiance Calculation
    !-------------------------
    CALL GEMS_Share_Calc_RadianceValue(gds_brvg%brvg_dfld%RadMantissa,    gds_brvg%brvg_dfld%RadExponent,                       &
                                       gvf_evis_rad)
    !print*, 'gvf_evis_rad  (1,1,1)=', gvf_evis_rad(1,1,1)

    !-------------------------
    !-- Initialize For Sun UV2 Radiance Calculation-6th
    !-------------------------
    nArg1 = gci_nwavel2
    CALL GEMS_Share_deAllocMem4Rad       ! for memory deallocation in the share module to calculate radiance
    CALL GEMS_Share_init4RadianceCalc(nArg1, nArg2, nArg3, nArg4) ! for using the share module to calculate radiance
    CALL GEMS_Share_allocMem4Rad                                  ! for memory allocation of the share module to calculate radiance

    !-------------------------
    !-- Sun UV2 Radiance Calculation
    !-------------------------
    CALL GEMS_Share_Calc_RadianceValue(gds_birr%suv2%birr_dfld%IrrRadMantissa,  gds_birr%suv2%birr_dfld%IrrRadExponent,          &
                                       gvf_suv2_irrad)
    CALL GEMS_Share_Calc_RadianceValue(gds_birr%suv2%birr_dfld%IrrRadPrecision, gds_birr%suv2%birr_dfld%IrrRadExponent,          &
                                       gvf_suv2_prec)
    !print*, 'gvf_suv2_irrad(1,1  )=', gvf_suv2_irrad(1,1)

    !-------------------------
    !-- Initialize For Sun VIS Radiance Calculation-7th
    !-------------------------
    nArg1 = gci_nwavel3
    CALL GEMS_Share_deAllocMem4Rad       ! for memory deallocation in the share module to calculate radiance
    CALL GEMS_Share_init4RadianceCalc(nArg1, nArg2, nArg3, nArg4) ! for using the share module to calculate radiance
    CALL GEMS_Share_allocMem4Rad                                  ! for memory allocation of the share module to calculate radiance

    !-------------------------
    !-- Sun VIS Radiance Calculation
    !-------------------------
    CALL GEMS_Share_Calc_RadianceValue(gds_birr%svis%birr_dfld%IrrRadMantissa, gds_birr%svis%birr_dfld%IrrRadExponent,          &
                                       gvf_svis_irrad)
    !print*, 'gvf_svis_irrad(1,1  )=', gvf_svis_irrad(1,1)

    !-------------------------------
    !-- Memory DeAllocation For Normal
    !-------------------------------

    CALL GEMS_Share_deAllocMem4Rad       ! for memory deallocation in the share module to calculate radiance

    !CALL GEMS_Share_DeallocL1BVars(gds_brug, gds_brvg, gds_birr)
    !CALL GEMS_Share_DeallocMem2

    print*, "-- End GEMS_Share_GetRadiance Subroutine --"
    print*, " "

    RETURN

1999  CONTINUE

    !-------------------------------
    !-- Memory DeAllocation For Error
    !-------------------------------

    CALL GEMS_Share_deAllocMem4Rad       ! for memory deallocation in the share module to calculate radiance

    !CALL GEMS_Share_DeallocL1BVars(gds_brug, gds_brvg, gds_birr)
    !CALL GEMS_Share_DeallocMem2

    print*, "-- ERROR GEMS_Share_GetRadiance Subroutine --"
    print*, " "

    RETURN

END SUBROUTINE GEMS_Share_GetRadiance

!-------------------------------------------------------------------------------
!+Description: 
!      GEMS Share Module Radiance Getting Program .
!
!
! Method:
!
!
! Input :
!       L1BRUG file
!
! Output:
!       Radiance, Wavelength value
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2015.04.02 Fisrt  Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------

SUBROUTINE GEMS_Share_GetRadiance2      ! Getting data from L1BRUG

! Local parameters For L1B Reading:
    CHARACTER(LEN=128)      :: ctl_fpath

    !TYPE(L1B_rug)           :: L1Brug       !-- L1BRUG data structure

    INTEGER(KIND=4)         :: status

    INTEGER(KIND=4)         :: nArg1, nArg2, nArg3, nArg4

    !--------------------------------
    !---- GEMS_Share_GetRadiance PROGRAM Start ---
    !--------------------------------
    LOGMSG=" PROC MSG"
    llvl = 9                     ! Log Level
    CALL GEMS_Share_MOD_log(llvl, '--- Start GEMS_Share_GetRadiance2 PROGRAM ---', 'START MSG')
    
    !-------------------------------
    !-- Initialize Global Constants
    !-------------------------------
    CALL GEMS_Share_Init_GlobalConstants
    CALL GEMS_Share_Init_L1BGlobalConstants

    !-------------------------------
    !-- Memory Allocation
    !-------------------------------
    CALL GEMS_Share_AllocMem2       ! Radiance 전역변수의 메모리를 할당

    !-------------------------------
    !-- Namelist setting
    !-------------------------------
    ctl_fpath = gcc_l1b_ctl_file

    !-------------------------
    !-- GEMS L1B read module
    !-------------------------
    !print*, 'L1BRUG File READ'
    status = GEMS_Share_MOD_L1Brug_Read(gds_brug, ctl_fpath)
    IF ( status .ne. 0 ) GOTO 1999

    !-------------------------
    !-- Initialize For Earth UV1 Wavelength Calculation-1st
    !-------------------------
    nArg1 = gci_nwavel
    nArg2 = gci_nxtrack
    nArg3 = gci_ntimes
    nArg4 = gci_nwavelcoef
    CALL GEMS_Share_init4RadianceCalc(nArg1, nArg2, nArg3, nArg4) ! for using the share module to calculate radiance

    !-------------------------
    !-- Earth UV1 Wavelength Calculation
    !-------------------------
    CALL GEMS_Share_Calc_ObservationWavelen(gds_brug%euv1%brug_dfld%WlenRefCol, gds_brug%euv1%brug_dfld%WlenCoef,               &
                                            gvf_euv1_radwave)

    !-------------------------
    !-- Initialize For Earth UV2 Wavelength  Calculation-2nd
    !-------------------------
    nArg1 = gci_nwavel2
    nArg2 = gci_nxtrack2
    CALL GEMS_Share_init4RadianceCalc(nArg1, nArg2, nArg3, nArg4) ! for using the share module to calculate radiance

    !-------------------------
    !-- Earth UV2 Wavelength Calculation
    !-------------------------
    CALL GEMS_Share_Calc_ObservationWavelen(gds_brug%euv2%brug_dfld%WlenRefCol, gds_brug%euv2%brug_dfld%WlenCoef,               &
                                            gvf_euv2_radwave)

    !-------------------------
    !-- Initialize For Earth UV1 Radiance Calculation-3rd
    !-------------------------
    nArg1 = gci_nwavel
    nArg2 = gci_nxtrack
    CALL GEMS_Share_init4RadianceCalc(nArg1, nArg2, nArg3, nArg4) ! for using the share module to calculate radiance
    CALL GEMS_Share_allocMem4Rad                                  ! for memory allocation of the share module to calculate radiance

    !-------------------------
    !-- Earth UV1 Radiance Calculation
    !-------------------------
    CALL GEMS_Share_Calc_RadianceValue(gds_brug%euv1%brug_dfld%RadMantissa, gds_brug%euv1%brug_dfld%RadExponent,  &
                                       gvf_euv1_rad)
    CALL GEMS_Share_Calc_RadianceValue(gds_brug%euv1%brug_dfld%RadPrecision, gds_brug%euv1%brug_dfld%RadExponent, &
                                       gvf_euv1_prec)

    !-------------------------
    !-- Initialize For Earth UV2 Radiance Calculation-4th
    !-------------------------
    nArg1 = gci_nwavel2
    nArg2 = gci_nxtrack2
    CALL GEMS_Share_deAllocMem4Rad       ! for memory deallocation in the share module to calculate radiance
    CALL GEMS_Share_init4RadianceCalc(nArg1, nArg2, nArg3, nArg4) ! for using the share module to calculate radiance
    CALL GEMS_Share_allocMem4Rad                                  ! for memory allocation of the share module to calculate radiance

    !-------------------------
    !-- Earth UV2 Radiance Calculation
    !-------------------------
    CALL GEMS_Share_Calc_RadianceValue(gds_brug%euv2%brug_dfld%RadMantissa, gds_brug%euv2%brug_dfld%RadExponent,  &
                                       gvf_euv2_rad)
    CALL GEMS_Share_Calc_RadianceValue(gds_brug%euv2%brug_dfld%RadPrecision, gds_brug%euv2%brug_dfld%RadExponent, &
                                       gvf_euv2_prec)

    !-------------------------------
    !-- Memory DeAllocation For Normal
    !-------------------------------

    CALL GEMS_Share_deAllocMem4Rad       ! for memory deallocation in the share module to calculate radiance

    !CALL GEMS_Share_DeallocMem2

    print*, "-- End GEMS_Share_GetRadiance2 Subroutine --"
    print*, " "

    RETURN

1999  CONTINUE

    !-------------------------------
    !-- Memory DeAllocation For Error
    !-------------------------------

    CALL GEMS_Share_deAllocMem4Rad       ! for memory deallocation in the share module to calculate radiance

    !CALL GEMS_Share_DeallocMem2

    print*, "-- ERROR GEMS_Share_GetRadiance2 Subroutine --"
    print*, " "

    RETURN

END SUBROUTINE GEMS_Share_GetRadiance2

!-------------------------------------------------------------------------------
!+Description: 
!      GEMS Share Module Radiance Getting Program 3.
!
!
! Method:
!
!
! Input :
!       L1BRUG, L1BIRR file
!
! Output:
!       Radiance, Irradiance, Wavelength value
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2015.04.02 Fisrt  Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------

SUBROUTINE GEMS_Share_GetRadiance3      ! Getting data from L1BRUG and L1BIRR 

! Local parameters For L1B Reading:
    CHARACTER(LEN=128)      :: ctl_fpath

    !TYPE(L1B_rug)           :: L1Brug       !-- L1BRUG data structure
    !TYPE(L1B_irr)           :: L1Birr       !-- L1BIRR data structure

    INTEGER(KIND=4)         :: status

    INTEGER(KIND=4)         :: nArg1, nArg2, nArg3, nArg4

    !--------------------------------
    !---- GEMS_Share_GetRadiance PROGRAM Start ---
    !--------------------------------
    LOGMSG=" PROC MSG"
    llvl = 9                     ! Log Level
    CALL GEMS_Share_MOD_log(llvl, '--- Start GEMS_Share_GetRadiance3 PROGRAM ---', 'START MSG')
    
    !-------------------------------
    !-- Initialize Global Constants
    !-------------------------------
    CALL GEMS_Share_Init_GlobalConstants
    CALL GEMS_Share_Init_L1BGlobalConstants

    !-------------------------------
    !-- Memory Allocation
    !-------------------------------
    CALL GEMS_Share_AllocMem2       ! Radiance 전역변수의 메모리를 할당

    !-------------------------------
    !-- Namelist setting
    !-------------------------------
    ctl_fpath = gcc_l1b_ctl_file

    !-------------------------
    !-- GEMS L1B read module
    !-------------------------
    !print*, 'L1BRUG File READ'
    status = GEMS_Share_MOD_L1Brug_Read(gds_brug, ctl_fpath)
    IF ( status .ne. 0 ) GOTO 1999

    !print*, 'L1BIRR File READ'
    status = GEMS_Share_MOD_L1Birr_Read(gds_birr, ctl_fpath)
    IF ( status .ne. 0 ) GOTO 1999

    !-------------------------
    !-- Initialize For Earth UV1 Wavelength Calculation-1st
    !-------------------------
    nArg1 = gci_nwavel
    nArg2 = gci_nxtrack
    nArg3 = gci_ntimes
    nArg4 = gci_nwavelcoef
    CALL GEMS_Share_init4RadianceCalc(nArg1, nArg2, nArg3, nArg4) ! for using the share module to calculate radiance

    !-------------------------
    !-- Earth UV1 Wavelength Calculation
    !-------------------------
    CALL GEMS_Share_Calc_ObservationWavelen(gds_brug%euv1%brug_dfld%WlenRefCol, gds_brug%euv1%brug_dfld%WlenCoef,               &
                                            gvf_euv1_radwave)

    !-------------------------
    !-- Initialize For Earth UV2 and Sun UV2 Wavelength  Calculation-2nd
    !-------------------------
    nArg1 = gci_nwavel2
    nArg2 = gci_nxtrack2
    CALL GEMS_Share_init4RadianceCalc(nArg1, nArg2, nArg3, nArg4) ! for using the share module to calculate radiance

    !-------------------------
    !-- Earth UV2 Wavelength Calculation
    !-------------------------
    CALL GEMS_Share_Calc_ObservationWavelen(gds_brug%euv2%brug_dfld%WlenRefCol, gds_brug%euv2%brug_dfld%WlenCoef,               &
                                            gvf_euv2_radwave)

    !-------------------------
    !-- Sun UV2 Wavelength Calculation
    !-------------------------
    CALL GEMS_Share_Calc_ObservationWavelen(gds_birr%suv2%birr_dfld%WlenRefCol, gds_birr%suv2%birr_dfld%WlenCoef,               &
                                            gvf_suv2_irrwave)

    !-------------------------
    !-- Initialize For Sun VIS Wavelength  Calculation-3rd
    !-------------------------
    nArg1 = gci_nwavel3
    CALL GEMS_Share_init4RadianceCalc(nArg1, nArg2, nArg3, nArg4) ! for using the share module to calculate radiance

    !-------------------------
    !-- Sun VIS Wavelength Calculation
    !-------------------------
    CALL GEMS_Share_Calc_ObservationWavelen(gds_birr%svis%birr_dfld%WlenRefCol, gds_birr%svis%birr_dfld%WlenCoef,               &
                                            gvf_svis_irrwave)

    !-------------------------
    !-- Initialize For Earth UV1 Radiance Calculation-4th
    !-------------------------
    nArg1 = gci_nwavel
    nArg2 = gci_nxtrack
    CALL GEMS_Share_init4RadianceCalc(nArg1, nArg2, nArg3, nArg4) ! for using the share module to calculate radiance
    CALL GEMS_Share_allocMem4Rad                                  ! for memory allocation of the share module to calculate radiance

    !-------------------------
    !-- Earth UV1 Radiance Calculation
    !-------------------------
    CALL GEMS_Share_Calc_RadianceValue(gds_brug%euv1%brug_dfld%RadMantissa, gds_brug%euv1%brug_dfld%RadExponent,  &
                                       gvf_euv1_rad)
    CALL GEMS_Share_Calc_RadianceValue(gds_brug%euv1%brug_dfld%RadPrecision, gds_brug%euv1%brug_dfld%RadExponent, &
                                       gvf_euv1_prec)

    !-------------------------
    !-- Initialize For Earth UV2 Radiance Calculation-5th
    !-------------------------
    nArg1 = gci_nwavel2
    nArg2 = gci_nxtrack2
    CALL GEMS_Share_deAllocMem4Rad       ! for memory deallocation in the share module to calculate radiance
    CALL GEMS_Share_init4RadianceCalc(nArg1, nArg2, nArg3, nArg4) ! for using the share module to calculate radiance
    CALL GEMS_Share_allocMem4Rad                                  ! for memory allocation of the share module to calculate radiance

    !-------------------------
    !-- Earth UV2 Radiance Calculation
    !-------------------------
    CALL GEMS_Share_Calc_RadianceValue(gds_brug%euv2%brug_dfld%RadMantissa, gds_brug%euv2%brug_dfld%RadExponent,  &
                                       gvf_euv2_rad)
    CALL GEMS_Share_Calc_RadianceValue(gds_brug%euv2%brug_dfld%RadPrecision, gds_brug%euv2%brug_dfld%RadExponent, &
                                       gvf_euv2_prec)

    !-------------------------
    !-- Initialize For Sun UV2 Radiance Calculation-6th
    !-------------------------
    CALL GEMS_Share_deAllocMem4Rad       ! for memory deallocation in the share module to calculate radiance
    CALL GEMS_Share_init4RadianceCalc(nArg1, nArg2, nArg3, nArg4) ! for using the share module to calculate radiance
    CALL GEMS_Share_allocMem4Rad                                  ! for memory allocation of the share module to calculate radiance

    !-------------------------
    !-- Sun UV2 Radiance Calculation
    !-------------------------
    CALL GEMS_Share_Calc_RadianceValue(gds_birr%suv2%birr_dfld%IrrRadMantissa,  gds_birr%suv2%birr_dfld%IrrRadExponent,          &
                                       gvf_suv2_irrad)
    CALL GEMS_Share_Calc_RadianceValue(gds_birr%suv2%birr_dfld%IrrRadPrecision, gds_birr%suv2%birr_dfld%IrrRadExponent,          &
                                       gvf_suv2_prec)

    !-------------------------
    !-- Initialize For Sun VIS Radiance Calculation-7th
    !-------------------------
    nArg1 = gci_nwavel3
    CALL GEMS_Share_deAllocMem4Rad       ! for memory deallocation in the share module to calculate radiance
    CALL GEMS_Share_init4RadianceCalc(nArg1, nArg2, nArg3, nArg4) ! for using the share module to calculate radiance
    CALL GEMS_Share_allocMem4Rad                                  ! for memory allocation of the share module to calculate radiance

    !-------------------------
    !-- Sun VIS Radiance Calculation
    !-------------------------
    CALL GEMS_Share_Calc_RadianceValue(gds_birr%svis%birr_dfld%IrrRadMantissa, gds_birr%svis%birr_dfld%IrrRadExponent,          &
                                       gvf_svis_irrad)

    !-------------------------------
    !-- Memory DeAllocation For Normal
    !-------------------------------

    CALL GEMS_Share_deAllocMem4Rad       ! for memory deallocation in the share module to calculate radiance

    !CALL GEMS_Share_DeallocMem2

    print*, "-- End GEMS_Share_GetRadiance3 Subroutine --"
    print*, " "

    RETURN

1999  CONTINUE

    !-------------------------------
    !-- Memory DeAllocation For Error
    !-------------------------------

    CALL GEMS_Share_deAllocMem4Rad       ! for memory deallocation in the share module to calculate radiance

    !CALL GEMS_Share_DeallocMem2

    print*, "-- ERROR GEMS_Share_GetRadiance3 Subroutine --"
    print*, " "

    RETURN

END SUBROUTINE GEMS_Share_GetRadiance3

!-------------------------------------------------------------------------------
!+Description: 
!      GEMS Share Module Radiance Getting Program 4.
!
!
! Method:
!
!
! Input :
!       L1BRVG, L1BIRR file
!
! Output:
!       Radiance, Irradiance, Wavelength value
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2015.04.02 Fisrt  Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------

SUBROUTINE GEMS_Share_GetRadiance4      ! Getting data from L1BRVG and L1BIRR 

! Local parameters For L1B Reading:
    CHARACTER(LEN=128)      :: ctl_fpath

    !TYPE(L1B_rvg)           :: L1Brvg       !-- L1BRVG data structure
    !TYPE(L1B_irr)           :: L1Birr       !-- L1BIRR data structure

    INTEGER(KIND=4)         :: status

    INTEGER(KIND=4)         :: nArg1, nArg2, nArg3, nArg4

    !--------------------------------
    !---- GEMS_Share_GetRadiance PROGRAM Start ---
    !--------------------------------
    LOGMSG=" PROC MSG"
    llvl = 9                     ! Log Level
    CALL GEMS_Share_MOD_log(llvl, '--- Start GEMS_Share_GetRadiance4 PROGRAM ---', 'START MSG')
    
    !-------------------------------
    !-- Initialize Global Constants
    !-------------------------------
    CALL GEMS_Share_Init_GlobalConstants
    CALL GEMS_Share_Init_L1BGlobalConstants

    !-------------------------------
    !-- Memory Allocation
    !-------------------------------
    CALL GEMS_Share_AllocMem2       ! Radiance 전역변수의 메모리를 할당

    !-------------------------------
    !-- Namelist setting
    !-------------------------------
    ctl_fpath = gcc_l1b_ctl_file

    !-------------------------
    !-- GEMS L1B read module
    !-------------------------
    !print*, 'L1BRVG File READ'
    print*,"TEST--------1"
    print*,ctl_fpath
    status = GEMS_Share_MOD_L1Brvg_Read(gds_brvg, ctl_fpath)
    IF ( status .ne. 0 ) GOTO 1999

    !print*, 'L1BIRR File READ'
    print*,"TEST--------2"
    status = GEMS_Share_MOD_L1Birr_Read(gds_birr, ctl_fpath)
    IF ( status .ne. 0 ) GOTO 1999

    !-------------------------
    !-- Initialize For Sun UV2 Wavelength  Calculation-1st
    !-------------------------
    nArg1 = gci_nwavel2
    nArg2 = gci_nxtrack2
    nArg3 = gci_ntimes
    nArg4 = gci_nwavelcoef
    CALL GEMS_Share_init4RadianceCalc(nArg1, nArg2, nArg3, nArg4) ! for using the share module to calculate radiance

    !-------------------------
    !-- Sun UV2 Wavelength Calculation
    !-------------------------
    CALL GEMS_Share_Calc_ObservationWavelen(gds_birr%suv2%birr_dfld%WlenRefCol, gds_birr%suv2%birr_dfld%WlenCoef,               &
                                            gvf_suv2_irrwave)

    !-------------------------
    !-- Initialize For Earth VIS and Sun VIS Wavelength  Calculation-2nd
    !-------------------------
    nArg1 = gci_nwavel3
    CALL GEMS_Share_init4RadianceCalc(nArg1, nArg2, nArg3, nArg4) ! for using the share module to calculate radiance

    !-------------------------
    !-- Earth VIS Wavelength Calculation
    !-------------------------
    CALL GEMS_Share_Calc_ObservationWavelen(gds_brvg%brvg_dfld%WlenRefCol, gds_brvg%brvg_dfld%WlenCoef, gvf_evis_radwave)

    !-------------------------
    !-- Sun VIS Wavelength Calculation
    !-------------------------
    CALL GEMS_Share_Calc_ObservationWavelen(gds_birr%svis%birr_dfld%WlenRefCol, gds_birr%svis%birr_dfld%WlenCoef,               &
                                            gvf_svis_irrwave)

    !-------------------------
    !-- Initialize For Sun UV2 Radiance Calculation-3rd
    !-------------------------
    nArg1 = gci_nwavel2
    CALL GEMS_Share_deAllocMem4Rad       ! for memory deallocation in the share module to calculate radiance
    CALL GEMS_Share_init4RadianceCalc(nArg1, nArg2, nArg3, nArg4) ! for using the share module to calculate radiance
    CALL GEMS_Share_allocMem4Rad                                  ! for memory allocation of the share module to calculate radiance

    !-------------------------
    !-- Sun UV2 Radiance Calculation
    !-------------------------
    CALL GEMS_Share_Calc_RadianceValue(gds_birr%suv2%birr_dfld%IrrRadMantissa,  gds_birr%suv2%birr_dfld%IrrRadExponent,          &
                                       gvf_suv2_irrad)
    CALL GEMS_Share_Calc_RadianceValue(gds_birr%suv2%birr_dfld%IrrRadPrecision, gds_birr%suv2%birr_dfld%IrrRadExponent,          &
                                       gvf_suv2_prec)

    !-------------------------
    !-- Initialize For Earth VIS Radiance Calculation-4th
    !-------------------------
    nArg1 = gci_nwavel3
    CALL GEMS_Share_init4RadianceCalc(nArg1, nArg2, nArg3, nArg4) ! for using the share module to calculate radiance
    CALL GEMS_Share_allocMem4Rad                                  ! for memory allocation of the share module to calculate radiance

    !-------------------------
    !-- Earth VIS Radiance Calculation
    !-------------------------
    CALL GEMS_Share_Calc_RadianceValue(gds_brvg%brvg_dfld%RadMantissa,    gds_brvg%brvg_dfld%RadExponent,                       &
                                       gvf_evis_rad)

    !-------------------------
    !-- Initialize For Sun VIS Radiance Calculation-5th
    !-------------------------
    CALL GEMS_Share_deAllocMem4Rad       ! for memory deallocation in the share module to calculate radiance
    CALL GEMS_Share_init4RadianceCalc(nArg1, nArg2, nArg3, nArg4) ! for using the share module to calculate radiance
    CALL GEMS_Share_allocMem4Rad                                  ! for memory allocation of the share module to calculate radiance

    !-------------------------
    !-- Sun VIS Radiance Calculation
    !-------------------------
    CALL GEMS_Share_Calc_RadianceValue(gds_birr%svis%birr_dfld%IrrRadMantissa, gds_birr%svis%birr_dfld%IrrRadExponent,          &
                                       gvf_svis_irrad)

    !-------------------------------
    !-- Memory DeAllocation For Normal
    !-------------------------------

    CALL GEMS_Share_deAllocMem4Rad       ! for memory deallocation in the share module to calculate radiance

    !CALL GEMS_Share_DeallocMem2

    print*, "-- End GEMS_Share_GetRadiance4 Subroutine --"
    print*, " "

    RETURN

1999  CONTINUE

    !-------------------------------
    !-- Memory DeAllocation For Error
    !-------------------------------

    CALL GEMS_Share_deAllocMem4Rad       ! for memory deallocation in the share module to calculate radiance

    !CALL GEMS_Share_DeallocMem2

    print*, "-- ERROR GEMS_Share_GetRadiance4 Subroutine --"
    print*, " "

    RETURN

END SUBROUTINE GEMS_Share_GetRadiance4
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------

SUBROUTINE GEMS_Share_GetRadiance5(L1Brug, L1Brvg, L1Birr)

! Local parameters For L1B Reading:
    CHARACTER(LEN=128)      :: ctl_fpath

    TYPE(L1B_rug)           :: L1Brug       !-- L1BRUG data structure
    TYPE(L1B_rvg)           :: L1Brvg       !-- L1BRVG data structure
    TYPE(L1B_irr)           :: L1Birr       !-- L1BIRR data structure

    INTEGER(KIND=4)         :: status

    INTEGER(KIND=4)         :: nArg1, nArg2, nArg3, nArg4

    !--------------------------------
    !---- GEMS_Share_GetRadiance PROGRAM Start ---
    !--------------------------------
    LOGMSG=" PROC MSG"
    llvl = 9                     ! Log Level
    CALL GEMS_Share_MOD_log(llvl, '--- Start GEMS_Share_GetRadiance PROGRAM ---', 'START MSG')
    
    !-------------------------------
    !-- Initialize Global Constants
    !-------------------------------
    CALL GEMS_Share_Init_GlobalConstants
    CALL GEMS_Share_Init_L1BGlobalConstants

    !-------------------------------
    !-- Memory Allocation
    !-------------------------------
    CALL GEMS_Share_AllocMem2       ! Radiance 전역변수의 메모리를 할당

    !/-------------------------------
    !-- Namelist setting
    !-------------------------------
!    ctl_fpath = gcc_l1b_ctl_file

    !-------------------------
    !-- GEMS L1B read module
    !-------------------------
    !print*, 'L1BRUG File READ'
!    status = GEMS_Share_MOD_L1Brug_Read(gds_brug, ctl_fpath)
!    IF ( status .ne. 0 ) GOTO 1999
!
!    !print*, 'L1BRVG File READ'
!    status = GEMS_Share_MOD_L1Brvg_Read(gds_brvg, ctl_fpath)
!    IF ( status .ne. 0 ) GOTO 1999

    !print*, 'L1BIRR File READ'
!    status = GEMS_Share_MOD_L1Birr_Read(gds_birr, ctl_fpath)
!    IF ( status .ne. 0 ) GOTO 1999

    !-------------------------
    !-- Initialize For Earth UV1 Wavelength Calculation-1st
    !-------------------------
    nArg1 = gci_nwavel
    nArg2 = gci_nxtrack
    nArg3 = gci_ntimes
    nArg4 = gci_nwavelcoef
    CALL GEMS_Share_init4RadianceCalc(nArg1, nArg2, nArg3, nArg4) ! for using the share module to calculate radiance

    !-------------------------
    !-- Wavelength Calculation
    !-------------------------
    CALL GEMS_Share_Calc_ObservationWavelen(L1Brug%euv1%brug_dfld%WlenRefCol, L1Brug%euv1%brug_dfld%WlenCoef,               &
                                            gvf_euv1_radwave)
    !print*,'gvf_euv1_radwave(1,1,1)=', gvf_euv1_radwave(1,1,1)

    !-------------------------
    !-- Initialize For Earth UV2 and Sun UV2 Wavelength  Calculation-2nd
    !-------------------------
    nArg1 = gci_nwavel2
    nArg2 = gci_nxtrack2
    CALL GEMS_Share_init4RadianceCalc(nArg1, nArg2, nArg3, nArg4) ! for using the share module to calculate radiance

    CALL GEMS_Share_Calc_ObservationWavelen(L1Brug%euv2%brug_dfld%WlenRefCol, L1Brug%euv2%brug_dfld%WlenCoef,               &
                                            gvf_euv2_radwave)
    !print*,'gvf_euv2_radwave(1,1,1)=', gvf_euv2_radwave(1,1,1)

    CALL GEMS_Share_Calc_ObservationWavelen(L1Birr%suv2%birr_dfld%WlenRefCol, L1Birr%suv2%birr_dfld%WlenCoef,               &
                                            gvf_suv2_irrwave)
    !print*, 'gvf_suv2_irrwave  (1,1)=', gvf_suv2_irrwave(1,1)

    !-------------------------
    !-- Initialize For Earth VIS and Sun VIS Wavelength  Calculation-2nd-1
    !-------------------------
    nArg1 = gci_nwavel3
    CALL GEMS_Share_init4RadianceCalc(nArg1, nArg2, nArg3, nArg4) ! for using the share module to calculate radiance

    CALL GEMS_Share_Calc_ObservationWavelen(L1Brvg%brvg_dfld%WlenRefCol, L1Brvg%brvg_dfld%WlenCoef,                         &
                                            gvf_evis_radwave)
    !print*, 'gvf_evis_radwave  (1,1,1)=', gvf_evis_radwave(1,1,1)

    CALL GEMS_Share_Calc_ObservationWavelen(L1Birr%svis%birr_dfld%WlenRefCol, L1Birr%svis%birr_dfld%WlenCoef,               &
                                            gvf_svis_irrwave)
    !print*, 'gvf_svis_irrwave  (1,1)=', gvf_svis_irrwave(1,1)

    !-------------------------
    !-- Initialize For Earth UV1 Radiance Calculation-3rd
    !-------------------------
    nArg1 = gci_nwavel
    nArg2 = gci_nxtrack
    CALL GEMS_Share_init4RadianceCalc(nArg1, nArg2, nArg3, nArg4) ! for using the share module to calculate radiance
    CALL GEMS_Share_allocMem4Rad                                  ! for memory allocation of the share module to calculate radiance

    !-------------------------
    !-- UV1 Radiance Radiance Calculation
    !-------------------------
    CALL GEMS_Share_Calc_RadianceValue(L1Brug%euv1%brug_dfld%RadMantissa,    L1Brug%euv1%brug_dfld%RadExponent,             &
                                       gvf_euv1_rad)
    CALL GEMS_Share_Calc_RadianceValue(L1Brug%euv1%brug_dfld%RadPrecision, L1Brug%euv1%brug_dfld%RadExponent, &
                                       gvf_euv1_prec)
    !print*, 'gvf_euv1_rad  (1,1,1)=', gvf_euv1_rad(1,1,1)
!stop
    !-------------------------
    !-- Initialize For Earth UV2 Radiance Calculation-4th
    !-------------------------
    nArg1 = gci_nwavel2
    nArg2 = gci_nxtrack2
    CALL GEMS_Share_deAllocMem4Rad       ! for memory deallocation in the share module to calculate radiance
    CALL GEMS_Share_init4RadianceCalc(nArg1, nArg2, nArg3, nArg4) ! for using the share module to calculate radiance
    CALL GEMS_Share_allocMem4Rad                                  ! for memory allocation of the share module to calculate radiance

    !-------------------------
    !-- UV2 Radiance Calculation
    !-------------------------
    CALL GEMS_Share_Calc_RadianceValue(L1Brug%euv2%brug_dfld%RadMantissa, L1Brug%euv2%brug_dfld%RadExponent, &
                                       gvf_euv2_rad)
    CALL GEMS_Share_Calc_RadianceValue(L1Brug%euv2%brug_dfld%RadPrecision, L1Brug%euv2%brug_dfld%RadExponent, &
                                       gvf_euv2_prec)
    !print*, 'gvf_euv2_rad  (1,1,1)=', gvf_euv2_rad(1,1,1)

    !-------------------------
    !-- Initialize For Earth VIS Radiance Calculation-5th
    !-------------------------
    nArg1 = gci_nwavel3
    CALL GEMS_Share_deAllocMem4Rad       ! for memory deallocation in the share module to calculate radiance
    CALL GEMS_Share_init4RadianceCalc(nArg1, nArg2, nArg3, nArg4) ! for using the share module to calculate radiance
    CALL GEMS_Share_allocMem4Rad                                  ! for memory allocation of the share module to calculate radiance

    !-------------------------
    !-- Earth VIS Radiance Calculation
    !-------------------------
    CALL GEMS_Share_Calc_RadianceValue(L1Brvg%brvg_dfld%RadMantissa,    L1Brvg%brvg_dfld%RadExponent,                       &
                                       gvf_evis_rad)
    !print*, 'gvf_evis_rad  (1,1,1)=', gvf_evis_rad(1,1,1)

    !-------------------------
    !-- Initialize For Sun UV2 Radiance Calculation-6th
    !-------------------------
    nArg1 = gci_nwavel2
    CALL GEMS_Share_deAllocMem4Rad       ! for memory deallocation in the share module to calculate radiance
    CALL GEMS_Share_init4RadianceCalc(nArg1, nArg2, nArg3, nArg4) ! for using the share module to calculate radiance
    CALL GEMS_Share_allocMem4Rad                                  ! for memory allocation of the share module to calculate radiance

    !-------------------------
    !-- Sun UV2 Radiance Calculation
    !-------------------------
    CALL GEMS_Share_Calc_RadianceValue(L1Birr%suv2%birr_dfld%IrrRadMantissa,  L1Birr%suv2%birr_dfld%IrrRadExponent,          &
                                       gvf_suv2_irrad)
    CALL GEMS_Share_Calc_RadianceValue(L1Birr%suv2%birr_dfld%IrrRadPrecision, L1Birr%suv2%birr_dfld%IrrRadExponent,          &
                                       gvf_suv2_prec)
    !print*, 'gvf_suv2_irrad(1,1  )=', gvf_suv2_irrad(1,1)

    !-------------------------
    !-- Initialize For Sun VIS Radiance Calculation-7th
    !-------------------------
    nArg1 = gci_nwavel3
    CALL GEMS_Share_deAllocMem4Rad       ! for memory deallocation in the share module to calculate radiance
    CALL GEMS_Share_init4RadianceCalc(nArg1, nArg2, nArg3, nArg4) ! for using the share module to calculate radiance
    CALL GEMS_Share_allocMem4Rad                                  ! for memory allocation of the share module to calculate radiance

    !-------------------------
    !-- Sun VIS Radiance Calculation
    !-------------------------
    CALL GEMS_Share_Calc_RadianceValue(L1Birr%svis%birr_dfld%IrrRadMantissa, L1Birr%svis%birr_dfld%IrrRadExponent,          &
                                       gvf_svis_irrad)
    !print*, 'gvf_svis_irrad(1,1  )=', gvf_svis_irrad(1,1)

    !-------------------------------
    !-- Memory DeAllocation For Normal
    !-------------------------------

    CALL GEMS_Share_deAllocMem4Rad       ! for memory deallocation in the share module to calculate radiance

    !CALL GEMS_Share_DeallocL1BVars(L1Brug, L1Brvg, L1Birr)
    !CALL GEMS_Share_DeallocMem2

    print*, "-- End GEMS_Share_GetRadiance Subroutine --"
    print*, " "

    RETURN

1999  CONTINUE

    !-------------------------------
    !-- Memory DeAllocation For Error
    !-------------------------------

    CALL GEMS_Share_deAllocMem4Rad       ! for memory deallocation in the share module to calculate radiance

    !CALL GEMS_Share_DeallocL1BVars(L1Brug, L1Brvg, L1Birr)
    !CALL GEMS_Share_DeallocMem2

    print*, "-- ERROR GEMS_Share_GetRadiance Subroutine --"
    print*, " "

    RETURN

END SUBROUTINE GEMS_Share_GetRadiance5

END MODULE Share_MOD_Radiance2
