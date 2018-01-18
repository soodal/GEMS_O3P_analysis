!-------------------------------------------------------------------------------
!+

MODULE O3P_MOD_Input
!-------------------------------------------------------------------------------
!+Description: 
!     .
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.0     2016.05.11 Fisrt Code (R&D, SaeaSoft Co., Ltd.) 
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Modules:

 USE O3P_MOD_Declaration
 USE OMSAO_errstat_module
 USE OMSAO_parameters_module,     ONLY:maxchlen
 USE GEMS_O3P_gemsdata_module,    ONLY: do_xbin, do_ybin, nxbin, nybin,  ncoadd, nxtrack, &
                                        gems_ny,gems_nx, ntimes , allocate_o3p_var
 USE OMSAO_variables_module,      ONLY:pixnum_lim, linenum_lim, coadd_uv2
!**********************************************************!
 IMPLICIT NONE

! ---------
! Variables
! ---------

!**********************************************************!
CONTAINS

SUBROUTINE GEMS_O3P_SUB2_Proc_Input(fit_ctrl_file, pge_error_status)
    IMPLICIT NONE
   
    CHARACTER (Len=maxchlen), INTENT(IN)  :: fit_ctrl_file
    INTEGER(KIND=4), INTENT(INOUT)          :: pge_error_status  ! return error code

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER                         :: nxcoadd,i
    INTEGER,             PARAMETER  :: fcunit = 11,  specunit = 12


    !------------------------------------
    !---- Standard MAIN PROGRAM Start ---
    !------------------------------------
    pge_error_status = pge_errstat_ok
 

    print*, "+++ 2. Input "
    print*, "Here is GEMS_O3P_SUB2_Proc_Input() subroutine. !!!"
    print*, " "

    WRITE(*,'(A)') ' => Reading control option'
    CALL gems_o3p_read_ctr_file (fcunit,fit_ctrl_file, pge_error_status )
    IF ( pge_error_status >= pge_errstat_error ) THEN
         WRITE(*,'(A)') ' => Error in Reading fitting control file'
         RETURN
    END IF
    !
    !--
    !
    WRITE(*,'(A)') ' => Reading Reference Spectra'
    CALL read_reference_spectra ( specunit, pge_error_status )
    IF ( pge_error_status >= pge_errstat_error ) THEN
         WRITE(*,'(A)') ' => Error in read_reference_spectra'
         RETURN
    END IF 
    !
    !--
    !
    WRITE(*,'(A)') ' => Reading GEMS l1b Data from GEMS SHARE MODULE'
   ! Allocate & initialize gemsdata variables from gems share module
    CALL gems_o3p_share_l1b (pge_error_status)
    If (pge_error_status >= pge_errstat_error)  THEN
         WRITE(*,'(A)') ' => Error in gems_o3p_share_l1b'
         RETURN
    END IF
    
    !------------------------------------------------------------------------
    ! check pixnum_lim & linenum_lim
    !--------------------------------------------------------------------------
    ! boundary check
    IF ( ALL ( linenum_lim < 0 ) )         linenum_lim(1:2)    = (/ 1, ntimes/)
    IF ( linenum_lim(1) > linenum_lim(2) ) linenum_lim([1, 2]) = linenum_lim([2,1])
    IF ( linenum_lim(1) < 1 )              linenum_lim(1)      = 1
    IF ( linenum_lim(2) > ntimes)          linenum_lim(2)      = ntimes

    IF ( ALL ( pixnum_lim < 0 ) )          pixnum_lim(1:2)    = (/ 1, nxtrack /)
    IF ( pixnum_lim(1) > pixnum_lim(2) )   pixnum_lim([1, 2]) = pixnum_lim([2, 1])
    IF ( pixnum_lim(1) < 1 )               pixnum_lim(1)      = 1
    IF ( pixnum_lim(2) > nxtrack )        pixnum_lim(2)      = nxtrack
    
    ! linenum_lim check
    If (do_ybin .and. nybin > 1 ) THEN 
      i = linenum_lim(2) - linenum_lim(1) +1
      IF (MOD (i, nybin) /= 0 ) THEN 
          linenum_lim(2) = NINT(1.0*i/nybin)*nybin + linenum_lim(1) - 1
          IF (linenum_lim(2) > ntimes) linenum_lim(2) = linenum_lim(2) - nybin
          IF (linenum_lim(1) > linenum_lim(2)) THEN 
              print * , 'check linenum_lim in GEMS_O3P_sub2input' ; stop
          ENDIF
      ENDIF
    ELSE
      nybin = 1
    ENDIF

    ! check pixnum_lim
    ! check for selected across track position (must start from odd positions)
    IF (coadd_uv2)  THEN
     i = pixnum_lim(2)-pixnum_lim(1) + 1
     IF ( MOD(pixnum_lim(1), ncoadd) /= 1 .OR. MOD(i, ncoadd) /= 0 ) THEN
        WRITE(*, '(A,2I4)') 'Incorrect across track positions to be coadded: ',pixnum_lim(1:2)
        pge_error_status = pge_errstat_error; RETURN
     ENDIF
     ! pixnum_lim for UV2
     pixnum_lim = CEILING(1.0 * pixnum_lim / ncoadd)  !nint
     ! pixnum_lim for Uv1 if ncoadd = 2
    ENDIF

   ! must divide and must start from odd coadded positions
   IF (do_xbin .AND. nxbin > 1) THEN
     i = pixnum_lim(2)-pixnum_lim(1) + 1
     IF ( MOD (i, nxbin) /= 0 .OR. MOD(pixnum_lim(1), nxbin) /= 1 ) THEN
        print * , mod(i,nxbin), mod(pixnum_lim(1), nxbin)
        WRITE(*, '(A,2I4)') 'Incorrect across track binning option: ',pixnum_lim(1:2)
        pge_error_status = pge_errstat_error; RETURN
     ENDIF
   ELSE
     nxbin = 1
   ENDIF
    
    !---------------------------------------------------------------------

    nxcoadd    = ncoadd*nxbin            ! Need just for OMI
    gems_nx    = nxtrack/nxcoadd         ! N of total x-pixels after coadding
    gems_ny    = INT ( ntimes*1.0/nybin) ! N of total y-pixels after coadding
    !gems_ny = 10 !wasp
    
    call allocate_o3p_var(ntimes,  pge_error_status)
 
    !------------------------------
    ! Read Irradiance
    !------------------------------
    WRITE(*,'(A,i2,"-",i2)') ' => Reading GEMS irradiance for ', 1, gems_nx
    CALL gems_o3p_read_l1b_irrad ( nxcoadd, 1, gems_nx, pge_error_status)
    If (pge_error_status >= pge_errstat_error) THEN
         RETURN
    END IF

   
    RETURN

END SUBROUTINE GEMS_O3P_SUB2_Proc_Input



END MODULE O3P_MOD_Input
