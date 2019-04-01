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
 USE GEMS_O3P_gemsdata_module,    ONLY: do_xbin, do_ybin, nxbin, nybin,  ncoadd, nchannel, &
                                        gems_ny,gems_nx, ntimes , allocate_o3p_var, nfxtrack

 USE OMSAO_variables_module,      ONLY:pixnum_lim, linenum_lim,  &
                                       l1b_rad_filename
 USE SYNT_data_module             
 USE SYNT_read_l1b,               ONLY: read_synt_l1b,  get_synt_dims


!**********************************************************!
 IMPLICIT NONE

! ---------
! Variables
! ---------

!**********************************************************!
CONTAINS

SUBROUTINE SYNT_O3P_SUB2_Proc_Input(fit_ctrl_file, pge_error_status)
    IMPLICIT NONE
   
    CHARACTER (Len=maxchlen), INTENT(IN)  :: fit_ctrl_file
    INTEGER(KIND=4), INTENT(INOUT)        :: pge_error_status  ! return error code
    INTEGER (KIND=4)                      :: estat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER                         :: nxcoadd, npix
    INTEGER                         :: nx, ny, nw, didx
    INTEGER,             PARAMETER  :: fcunit = 11,  specunit = 12
   
    ! For synt
    INTEGER :: first_line, last_line, first_pix, last_pix
    REAL (KIND=r8) :: time_tai
    CHARACTER (LEN = 28)              :: tmp_utc !geun
    CHARACTER (LEN = 8)               :: deli !geun
    CHARACTER (LEN = maxchlen)        :: snr_filename !geun

    ! Externeal function 
    INTEGER (KIND=i4), EXTERNAL :: PGS_TD_UTCtoTAI, PGS_TC_TAItoUTC !geun

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

    WRITE(*,'(A)') ' => Reading SYNT l1b Data from SYNT MODULE'

    !Allocate & initialize gemsdata variables from gems share module
    CALL get_synt_dims (l1b_rad_filename, nxtrack, nytrack, nwavel)

    ! spatial binning is not implemented for omps
    if (nxtrack .ne. nxtrack_max .and. nytrack .ne. nytrack_max) then
        print *, 'nxtrack=',nxtrack,nxtrack_max,nytrack, nytrack_max, 'nwavel=',nwavel, l1b_rad_filename
        write(*,'(a)') 'this case is not performed '
        pge_error_status = pge_errstat_error ; return
    endif
    nxbin = 1 ;  nybin = 1

    if (linenum_lim(2) >= nytrack) linenum_lim(2) = nytrack

    first_pix  = CEILING(1.0 * pixnum_lim(1) / nxbin)
    last_pix   = NINT(1.0 * pixnum_lim(2) / nxbin )

    ! Line number, starting from zero and keep track of offset
    first_line  = CEILING(1.0 * linenum_lim(1) / nybin)
    last_line   = NINT(1.0 * linenum_lim(2) / nybin )

    WRITE(*,'(A)') '=> Reading radiance / irradiance / Geolocation files'
    call allocate_synt_raddata( nxtrack, nytrack, nwavel, max_fit_pts, maxwin, pge_error_status) 
                                                                                 ! add for using pointer variables : geun
    call allocate_synt_data( nxtrack, nytrack, nwavel, pge_error_status)

    !didx=INDEX(l1b_rad_filename, '.nc')  ;  didx=didx-19  ! total data
    didx=INDEX(l1b_rad_filename, '.nc')  ;  didx=didx-23   ! sliced data
    synt_db=l1b_rad_filename(:didx-10)    
    synt_num=l1b_rad_filename(didx+20:didx+22)
    synt_date=l1b_rad_filename(didx:didx+9)
    !snr_filename=TRIM(synt_db)//'Interpolated_SNR.nc3'  !geun
        
    !CALL read_synt_l1b (l1b_rad_filename, snr_filename, pge_error_status)
    CALL read_synt_l1b (l1b_rad_filename, pge_error_status)

    If (pge_error_status >= pge_errstat_error)  THEN
         !WRITE(*,'(A)') ' => Error in gems_o3p_share_l1b'
         WRITE(*,'(A)') ' => Error in reading l1b files'
         RETURN
    END IF

   
    !------------------------------------------------------------------------
    ! read geolocation data 
    !--------------------------------------------------------------------------
    ntimes = nytrack   ! geun
    nfxtrack = nxtrack  ! geun
    ny = ntimes
    nw = nwavel ; nx = nxtrack 

    tmp_utc=l1b_rad_filename(didx:didx+3)//'-'//l1b_rad_filename(didx+4:didx+5)//'-'//l1b_rad_filename(didx+6:didx+7)
    !tmp_utc=TRIM(tmp_utc)//'T00:00:00.000001Z'
    tmp_utc=TRIM(tmp_utc)//'T'//l1b_rad_filename(didx+8:didx+9)//':00:00.000001Z'
    estat             = PGS_TD_UTCtoTAI(tmp_utc,time_tai)
    synt_geo%time(1:ny)             = time_tai 
    !synt_geo%time(1:ny)             = 6.48006826138914E8    ! 20130715 for atm data 
    synt_geo%height(1:nx, 1:ny)     = synt_hter(1:nx,1:ny)
    synt_geo%lon(1:nx, 1:ny)        = synt_longitude(1:nx,1:ny) 
    synt_geo%lat(1:nx, 1:ny)        = synt_latitude(1:nx,1:ny)
    synt_geo%sza(1:nx, 1:ny)        = synt_szenith(1:nx, 1:ny)
    synt_geo%vza(1:nx, 1:ny)        = synt_vzenith(1:nx, 1:ny)
    synt_geo%saz(1:nx, 1:ny)        = synt_saa(1:nx, 1:ny) 
    synt_geo%vaz(1:nx, 1:ny)        = synt_vaa(1:nx, 1:ny)
    synt_geo%groundqflg(1:nx, 1:ny) = 0.0

    !------------------------------------------------------------------------
    ! check pixnum_lim & linenum_lim
    !--------------------------------------------------------------------------
    ! boundary check
    IF ( ALL ( linenum_lim < 0 ) )         linenum_lim(1:2)    = (/ 1, ntimes/)
    IF ( linenum_lim(1) > linenum_lim(2) ) linenum_lim([1, 2]) = linenum_lim([2,1])
    IF ( linenum_lim(1) < 1 )              linenum_lim(1)      = 1
    IF ( linenum_lim(2) > ntimes)          linenum_lim(2)      = ntimes

    IF ( ALL ( pixnum_lim < 0 ) )          pixnum_lim(1:2)    = (/ 1, nxtrack/)
    IF ( pixnum_lim(1) > pixnum_lim(2) )   pixnum_lim([1, 2]) = pixnum_lim([2,1])
    IF ( pixnum_lim(1) < 1 )               pixnum_lim(1)      = 1
    IF ( pixnum_lim(2) > nxtrack )         pixnum_lim(2)      = nxtrack
    
    ! linenum_lim check
    If (do_ybin .and. nybin > 1 ) THEN 
      npix = linenum_lim(2) - linenum_lim(1) +1
      IF (MOD (npix, nybin) /= 0 ) THEN 
          linenum_lim(2) = NINT(1.0*npix/nybin)*nybin + linenum_lim(1) - 1
          IF (linenum_lim(2) > ntimes) linenum_lim(2) = linenum_lim(2) - nybin
          IF (linenum_lim(1) > linenum_lim(2)) THEN 
              print * , 'check linenum_lim in GEMS_O3P_sub2input' ; stop
          ENDIF
      ENDIF
    ELSE
      nybin = 1
    ENDIF

    !---------------------------------------------------------------------

    gems_nx     = nxtrack                 ! N of total x-pixels after coadding
    gems_ny     = INT ( ntimes*1.0/nybin) ! N of total y-pixels after coadding
    !gems_ny    = 10    !geun
    nxcoadd     = 1
    nchannel    = 1
    ncoadd      = 1
    print*,'nxbin,nybin=',nxbin,nybin
    print*,'nchannel=',nchannel
    print*,'ncoadd=',ncoadd
    print*,'gems_nx=',gems_nx
    print*,'gems_ny=',gems_ny

    call allocate_o3p_var(ntimes,  pge_error_status)
 
    !!------------------------------
    !! Irradiance  -> gems irradiance variable
    !!------------------------------
    WRITE(*,'(A,i2,"-",i2)') ' => Reading SYNT irradiance for ', 1, gems_nx
    !CALL gems_o3p_read_l1b_irrad ( nxcoadd, 1, gems_nx, pge_error_status)  ! read climatological irradiance
    CALL synt_o3p_read_l1b_irrad ( nxcoadd, 1, gems_nx, pge_error_status)  ! read synthtic irradiance

    If (pge_error_status >= pge_errstat_error) THEN
         RETURN
    END IF
    RETURN

END SUBROUTINE SYNT_O3P_SUB2_Proc_Input

END MODULE O3P_MOD_Input
