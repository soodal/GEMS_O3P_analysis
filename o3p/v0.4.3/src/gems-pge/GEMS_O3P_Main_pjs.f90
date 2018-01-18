PROGRAM GEMS_O3P_Main

!--------------------------------------------------------------
!+Description :
!     GEMS Standard Main Program
!     
! Method:
!     Algorithm description file - 5 page
!
! Input files:
!     fit_ctrl_file : GEMS spectrum 자료를 가공하는데 필요한 옵션 및 변수가 저장되어 있는 파일
!
! Output files:
!    
!
!
!+Version  Date        Comment
! -------- ----------- -------
!  0.1.0   2011.       First Code for OMI     (Liu Xiong,  Harvard-SAO)
!  0.2.0   2015.07.01  Modified Code for GEMS (Bak Juseon, PNU UNIV.)
!  0.3.0               Improved Code for GEMS (Bak Juseon, PNU UNIV.)
!--------------------------------------------------------------

  USE OMSAO_errstat_module
  USE GEMS_O3P_gemsdata_module,    ONLY: nxbin, nybin, do_xbin, do_ybin, ncoadd, nfxtrack, &
                                         gems_gemsdata_deallocate, gems_ny, gems_nx, ntimes
 
  IMPLICIT NONE

  ! -------------------------
  ! Name of module/subroutine
  ! -------------------------
  CHARACTER (LEN=13), PARAMETER :: modulename = 'GEMS_O3P_Main'

  ! -------------------------
  ! LOCAL variables
  ! -------------------------
  INTEGER :: errstat, OMI_SMF_setmsg, nxcoadd
  INTEGER, PARAMETER               :: fcunit = 11,  specunit = 12
  CHARACTER (Len=100), PARAMETER   :: fit_ctrl_file='../run/conf/GEMS_O3P.inp'
  INTEGER :: pge_error_status, exit
  INTEGER, DIMENSION(2)           :: pixlim, linelim
  REAL (KIND=dp), DIMENSION(3)    :: fitcol
  REAL (KIND=dp), DIMENSION(3, 2) :: dfitcol
  INTEGER, PARAMETER              :: n_mpi = 3
  INTEGER, DIMENSION (n_mpi, 2)   :: mpi_line


  !------------------------------------
  !---- Standard MAIN PROGRAM Start ---
  !------------------------------------

  do_xbin = .FALSE. ; do_ybin =.FALSE.
  nxbin = 1 ; if (nxbin /= 1 ) do_xbin = .TRUE.
  nybin = 4 ; if (nybin /= 1)  do_ybin = .TRUE.

  pge_error_status = pge_errstat_ok

  CALL unbufferSTDout() ! Make PGE write STD/IO unbuffered

  WRITE(*,'(A)') '_________________________________'
  WRITE(*,'(A)') ' STEP1 : Setting algorithm'
  WRITE(*,'(A)') '_________________________________'

  WRITE(*,'(A)') ' => Reading fitting control file'
  CALL gems_o3p_read_fitting_control_file (fcunit,fit_ctrl_file, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) GOTO 666

  WRITE(*,'(A)') ' => Reading Reference Spectra'
  CALL read_reference_spectra ( specunit, pge_error_status )   
  IF ( pge_error_status >= pge_errstat_error ) GOTO 666
  
  WRITE(*,'(A)') ' => Loading GEMS l1b Data from GEMS SHARE MODULE'
 ! Allocate & initialize gemsdata variables from gems share module
  CALL gems_o3p_share_l1b (pge_error_status)
  If (pge_error_status >= pge_errstat_error)  GOTO 666

  nxcoadd    = ncoadd*nxbin
  gems_nx    = nfxtrack/nxbin          ! N of x-pixels after coadding
  gems_ny    = INT ( ntimes*1.0/nybin) ! N of y-pixels after coadding
  !------------------------------
  ! Read Irradiance
  !------------------------------
  WRITE(*,'(A,i2,"-",i2)') ' => Reading GEMS irradiance for ', 1, gems_nx
  CALL gems_o3p_read_l1b_irrad ( nxcoadd, 1, gems_nx, pge_error_status)
  If (pge_error_status >= pge_errstat_error) GOTO 666

  !------------------------------
  ! Calibrate Irradiance
  !------------------------------
  WRITE(*, '(A)') ' => Calibrating GEMS irradiance'
  CALL gems_o3p_irrad_cross_calibrate (1, gems_nx, pge_error_status)
  IF ( pge_error_status >= pge_errstat_error ) GOTO 666

  !------------------------------
  ! Setting for writing output
  !------------------------------
  CALL gems_o3p_write_h5 (1,fitcol, dfitcol,0)

   print * , ntimes
  Call  Assign_mpi_line   (n_mpi, mpi_line, .true.)

  WRITE(*,'(A)') '_________________________________'
  WRITE(*,'(A)') ' STEP2 : O3 fitting'
  WRITE(*,'(A)') '_________________________________'
  
  
  pixlim  = (/1,2/) ! [1,60]첫번째 픽셀은 반드시 홀수, 총 갯수는 반드시 짝수 


  linelim = (/mpi_line(1,1),mpi_line(1,2)/)! [1,1644]
  CALL Assign_pixline ( linelim, pixlim, pge_error_status)
  IF ( pge_error_status >= pge_errstat_error ) GOTO 666
  CALL gems_o3p_fitting_process (  pge_error_status)
  IF ( pge_error_status >= pge_errstat_error ) GOTO 666
 
  IF (n_mpi > 1 ) THEN 
   linelim = (/mpi_line(2,1),mpi_line(2,2)/)! [1,1644]
   print * , linelim
   CALL Assign_pixline ( linelim, pixlim, pge_error_status)
   IF ( pge_error_status >= pge_errstat_error ) GOTO 666
   CALL gems_o3p_fitting_process (  pge_error_status)
   IF ( pge_error_status >= pge_errstat_error ) GOTO 666
  ENDIF

  call gems_gemsdata_deallocate 
  WRITE(*,'(A)') '_________________________________'
  WRITE(*,'(A)') ' STEP3 : Saving Results to *.h5'
  WRITE(*,'(A)') '_________________________________'


  CALL gems_o3p_write_h5 (3,fitcol, dfitcol,1)
 

  ! ------------------------------------
  ! Write END_OF_RUN message to log file
  ! ------------------------------------
666 SELECT CASE ( pge_error_status )
  CASE ( pge_errstat_ok )
     ! ----------------------------------------------------------------
     ! PGE execution completed successfully. All is well.
     ! ----------------------------------------------------------------
     errstat = OMI_SMF_setmsg ( OMSAO_S_ENDOFRUN, '', modulename, 0 )
     print * , '!!!! algorithm  ended without any error'    
  CASE ( pge_errstat_warning )
     ! ----------------------------------------------------------------
     ! PGE execution raised non-terminal warnings. Nothing serious, we
     ! hope, so execution completed but with a non-zero exit status.
     ! ----------------------------------------------------------------
     errstat = OMI_SMF_setmsg ( OMSAO_W_ENDOFRUN, '', modulename, 0 )
     STOP 1
  CASE ( pge_errstat_error )
     ! ----------------------------------------------------------------
     ! PGE execution encountered an error that lead to termination.
     ! ----------------------------------------------------------------
     errstat = OMI_SMF_setmsg ( OMSAO_E_ENDOFRUN, '', modulename, 0 )
     STOP 1
  CASE DEFAULT
     ! ----------------------------------------------------------------
     ! If we ever reach here, then PGE_ERRSTAT has been set to a funny
     ! value. This should never happen, but we buffer this case anyway.
     ! ----------------------------------------------------------------
     errstat = OMI_SMF_setmsg ( OMSAO_U_ENDOFRUN, '', modulename, 0 )
     STOP
  END SELECT

END PROGRAM GEMS_O3P_Main

SUBROUTINE  Assign_mpi_line  (n_mpi, mpi_line, is_print)
 USE GEMS_O3P_gemsdata_module, ONLY: gems_ny,  nybin, lineloc
 IMPLICIT NONE

 ! IN/OUT variables
 INTEGER, INTENT (IN)            :: n_mpi 
 LOGICAL, INTENT(IN)             :: is_print
 INTEGER , DIMENSION (n_mpi, 2 ), INTENT(OUT) :: mpi_line 
 ! local variables
 INTEGER :: i, j, line1, line2, currline

 line1  = 1
 lineloc (:)   = 0
 mpi_line(1,1) = 1
 DO i = 1, n_mpi
    line2 = line1 + INT (gems_ny*1.0/n_mpi) -1
    IF ( i == n_mpi)  line2 = gems_ny  
    DO  j = line1, line2
        currline = (j-1)*nybin+1
        lineloc(currline) = j	
	print * , i, j, currline
    ENDDO
    IF (is_print) write(*,'("MPI of ",i4,":", i6,"-",i6)') i,line1, line2
    mpi_line(i,1) = (line1-1)*nybin +1 ; mpi_line(i,2) = line2*nybin
    line1 = line2 + 1

 ENDDO

 RETURN
END SUBROUTINE Assign_mpi_line

SUBROUTINE Assign_pixline ( linelim, pixlim , pge_error_status )
 
 USE OMSAO_errstat_module
 USE OMSAO_variables_module,  ONLY: pixnum_lim, linenum_lim,coadd_uv2
 USE GEMS_O3P_gemsdata_module, ONLY: nxtrack_max, ntimes_max,ntimes, nxtrack, &
                                     do_xbin, do_ybin, nxbin, nybin, ncoadd, &
				     offline, first_pix, last_pix, first_line,last_line
 
 IMPLICIT NONE

 !INPUT/OUTPUT variables
 INTEGER, INTENT(IN), DIMENSION(2)    :: pixlim, linelim
 INTEGER, INTENT(OUT) :: pge_error_status
 
 !local variables
 INTEGER :: i

  ! ------------------------------------------------
  ! Check for consistency of pixel limits to process
  ! -----------------------------------------------
  pge_error_status = pge_errstat_ok
 

  linenum_lim =  linelim
  pixnum_lim  =  pixlim

  IF ( ALL ( linenum_lim < 0 ) )         linenum_lim(1:2) = (/ 1, ntimes_max /)
  IF ( linenum_lim(1) > linenum_lim(2) ) linenum_lim([1, 2]) = linenum_lim([2, 1])   
  IF ( linenum_lim(1) < 1 )              linenum_lim(1) = 1
  IF ( linenum_lim(2) > ntimes)          linenum_lim(2) = ntimes

  IF ( ALL ( pixnum_lim < 0 ) )          pixnum_lim(1:2) = (/ 1, nxtrack_max /)
  IF ( pixnum_lim(1) > pixnum_lim(2) )   pixnum_lim([1, 2]) = pixnum_lim([2, 1])   
  IF ( pixnum_lim(1) < 1 )               pixnum_lim(1) = 1
  IF ( pixnum_lim(2) > nxtrack_max )     pixnum_lim(2) = nxtrack_max
  

  ! check for selected across track position (must start from odd positions)
  IF (coadd_uv2)  THEN
     i = pixnum_lim(2)-pixnum_lim(1) + 1
     IF ( MOD(pixnum_lim(1), ncoadd) /= 1 .OR. MOD(i, ncoadd) /= 0 ) THEN
        WRITE(*, '(A,2I4)') 'Incorrect across track positions to be coadded: ', pixnum_lim(1:2)
        pge_error_status = pge_errstat_error; RETURN
     ENDIF
     
     pixnum_lim = CEILING(1.0 * pixnum_lim / ncoadd)  !nint 
  ENDIF

  ! must divide and must start from odd coadded positions
  IF (do_xbin .AND. nxbin > 1) THEN
     i = pixnum_lim(2)-pixnum_lim(1) + 1
     IF ( MOD (i, nxbin) /= 0 .OR. MOD(pixnum_lim(1), nxbin) /= 1 ) THEN
        WRITE(*, '(A,2I4)') 'Incorrect across track binning option: ', pixnum_lim(1:2)
        pge_error_status = pge_errstat_error; RETURN
     ENDIF
  ELSE
     nxbin = 1
  ENDIF

  ! Could start from any positions, adjust line positions if necessary

  IF (do_ybin .AND. nybin > 1)  THEN
     i = linenum_lim(2)-linenum_lim(1) + 1
     IF (MOD(i, nybin) /= 0) THEN
        linenum_lim(2) = NINT(1.0 * i / nybin) * nybin + linenum_lim(1) - 1
	! print * , linenum_lim(2), ntimes
        IF (linenum_lim(2) > ntimes) linenum_lim(2) = linenum_lim(2) - nybin
     ENDIF
  ELSE
     nybin = 1
  ENDIF
  
  IF (nxbin == 1) do_xbin = .FALSE.
  IF (nybin == 1) do_ybin = .FALSE.

  
  first_pix  = CEILING(1.0 * pixnum_lim(1) / nxbin)
  last_pix   = NINT(1.0 * pixnum_lim(2) / nxbin )

  
  offline    = linenum_lim(1) -1
  first_line = 1 
  last_line    =  NINT (1.0 *(linenum_lim(2) - linenum_lim(1) ) / nybin)
  IF (nybin == 1) last_line  =  linenum_lim(2) - linenum_lim(1) + 1


  do i = 1, last_line
  !   print * , i, (i-1)*nybin + offline+1, (i-1)*nybin + offline+1+nybin-1
  enddo
  do i = first_pix, last_pix
   !  print * , i, (i-1)*nxbin+1, (i-1)*nxbin+1+nxbin -1
  enddo

 RETURN
END SUBROUTINE
