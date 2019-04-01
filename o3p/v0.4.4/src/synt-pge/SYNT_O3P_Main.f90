
PROGRAM SYNT_O3P_Main
!--------------------------------------------------------------
!+Description :
!     GEMS Standard Main Program  -> synthetic data
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
!  0.4.0   2016.06.15  MPI applied Code for GEMS (Bak Juseon, PNU UNIV.)
!
!--------------------------------------------------------------

!-------------------------------------------------------------------------------
! Modules 
!-------------------------------------------------------------------------------
! GEMS Shared modules
 USE Share_MOD_Log
 USE Share_l2_o3p_mod_write_read
 USE Share_MOD_MPI
 USE O3P_MOD_Declaration
 USE O3P_MOD_Input
 USE O3P_MOD_Output
 USE O3P_MOD_close
 USE O3P_MOD_share_close
! o3P only modules
 USE OMSAO_errstat_module
 USE OMSAO_variables_module,      ONLY:pixnum_lim, linenum_lim
 USE OMSAO_parameters_module,     ONLY:maxchlen
 USE GEMS_O3P_gemsdata_module,    ONLY:first_line, last_line, first_pix, last_pix, lineloc, pixloc
!**********************************************************!
 IMPLICIT NONE
!-------------------------- 
! GEMS MPI Variable Declaration
!----------------------------
 INTEGER(KIND=4)             :: nprocs   ! N of MPI process
 INTEGER(KIND=4)             :: myrank   ! current ID of MPI process [0, nprocs-1]
 INTEGER(KIND=4)             :: mpi_idx  ! myrank + 1
 CHARACTER(LEN=120)          :: str_attr ! error message
 CHARACTER(LEN=120)          :: log_file_name ! MPI process log file
 INTERFACE
   SUBROUTINE GEMS_Share_MOD_MPI_DataSumProcess(myrank, mpi_pix, mpi_line, L2O3P_wr, retCode)
     USE Share_l2_o3p_mod_write_read
     INTEGER(KIND=4), INTENT(IN)     :: myrank
     INTEGER(KIND=4), INTENT(IN)     :: mpi_pix(:)
     INTEGER(KIND=4), INTENT(IN)     :: mpi_line(:,:)
     TYPE(L2_o3p),    INTENT(OUT)    :: L2O3P_wr
     INTEGER(KIND=4), INTENT(OUT)    :: retCode
   END SUBROUTINE GEMS_Share_MOD_MPI_DataSumProcess
 END INTERFACE

 ! MPI는 y축방향으로 나누어짐 
 ! mpi_line 은 나누어서 계산할 line의 범위를 할당함, mpi_line_out은 계산후 최종
 ! 산출물이 저장될 위치를 할당함
 ! nybin = 1 면 mpi_line = mpi_line_out
 ! case: linenum_lim =[1, 14] , nprocs = 3 
 !                               (a) nybin = 1             (b) nybin = 2
 !       mpi_line(1,*) =[1,4]    o3(1:4) => o3_out(1:4)     o3(1:2)=> o3_out(1:2)
 !       mpi_line(2,*) =[5,8]    o3(1:4) => O3_out(5:8)     o3(1:2)=> o3_out(3:4)
 !       mpi_line(3,*) =[9,14]   o3(1:6) => o3_out(9:14)    o3(1:3)=> o3_out(5:7)
 ! case: linenum_lim=[11,24], nprocs = 3
 !                               (a) nybin = 1             (b) nybin = 2
 !       mpi_line(1,*) =[11,14]    o3(1:4) => o3_out(11:14)     o3(1:2)=> o3_out(6:7)
 !       mpi_line(2,*) =[15,18]    o3(1:4) => O3_out(15:18)     o3(1:2)=> o3_out(8:9)
 !       mpi_line(3,*) =[19,24]    o3(1:6) => o3_out(19:24)     o3(1:3)=> o3_out(10:12)
 INTEGER(KIND=4),ALLOCATABLE :: mpi_line(:,:)    
 INTEGER(KIND=4),ALLOCATABLE :: mpi_line_out(:,:)
 INTEGER(KIND=4)             :: mpi_pix_out(2)
 TYPE(L2_o3p)                :: L2O3P_wr ! variables for output of MPI
!----------------------------------------------
! other local variables
!----------------------------------------------
 CHARACTER (Len=maxchlen), PARAMETER  :: fit_ctrl_file='../run/conf/GEMS_O3P.inp'
 INTEGER(KIND=4)             :: pixlim(2), linelim(2) ! pixline limit for current mpi process
 INTEGER(KIND=4)             :: errstat, OMI_SMF_setmsg
 INTEGER(KIND=4)             :: errCode = 9       ! processing error flag
! -------------------------
! Name of module/subroutine
! -------------------------
 CHARACTER (LEN=13), PARAMETER :: modulename = 'SYNT_O3P_Main'


 !************************************************
 !******** O3P MAIN PROGRAM Start   ************ 
 !************************************************
 !-------------------------------------------------------------------------
 !-------- . MPI Start
 !-------------------------------------------------------------------------
  CALL GEMS_Share_MOD_MPI_Start(myrank, nprocs, errCode)
  IF( errCode /= 0 ) THEN
      str_attr = 'MPI Start Error!!!' ; GOTO 2999
  ENDIF
print*,'##GEMS_Share_MOD_MPI_Start End##'
 !-------------------------------------------------------------------------
 !-------- O3P MAIN PROGRAM Start                            ------------ 
 !-------------------------------------------------------------------------
    llvl   = 9
    LOGMSG ="START "
    CALL GEMS_Share_MOD_log(llvl,"---- START GEMS O3P Main ----" ,LOGMSG) 

    !
    !-- Open MPI LOG File 
    !
    log_file_name = "o3p_mpi_"
    CALL GEMS_Share_MOD_MPI_LogOpen(myrank, log_file_name, errCode)
    IF ( errCode /= 0 ) THEN
        PRINT*, "MPI Log File Open Error !"
        STOP
    ENDIF

    CALL unbufferSTDout() ! Make PGE write STD/IO unbuffered

    WRITE(*,'(A)') '_________________________________'
    WRITE(*,'(A)') ' STEP1 : Setting algorithm'
    WRITE(*,'(A)') '_________________________________'

    !-------------------------------------------------------------------------
    !-------- 1. Global Variable Declaration                      ------------ 
    !-------------------------------------------------------------------------
    ! 공통 모듈 사용을 위한 초기화 작업
    !CALL GEMS_O3P_SUB1_Declaration(errCode)
    !IF( errCode /= 0 ) THEN
        !str_attr = "Program Declaration Error!!!"
        !GOTO 2999
    !ENDIF
    !-------------------------------------------------------------------------
    !-------- 2. Input Processing                                 ------------ 
    !-------------------------------------------------------------------------
    !1. GEMS_O3P.inp 읽음
    !2. tbl 읽음
    !3. 공통모듈로부터 위성자료 읽음
    !4. o3p 알고리즘 변수 할당
    !5. irradiance 자료처리 (valid subset for o3 fitting window, coadding process)
    CALL SYNT_O3P_SUB2_Proc_Input(fit_ctrl_file, errCode)
    IF ( errCode /= 0 ) GOTO 1999
    !-------------------------------------------------------------------------
    !-------- 3. Irrad_Cross_Calibrate
    !-------------------------------------------------------------------------
    ! calibrate slit/wavelength 
    WRITE(*,'(A)') '=> GEMS O3P Irrad Cross Calibration'
    CALL GEMS_O3P_SUB3_Irrad_cross_calibrate(1, gems_nx, errCode)
    IF ( errCode /= 0 ) GOTO 1999
       
    !------------------------------
    ! MPI 설정 및 영역 설정
    !-----------------------------
    ALLOCATE (mpi_line(nprocs, 2), STAT = errCode)
    IF (errCode /= 0) THEN
        str_attr = "mpi_line Memory Allocation Error!!!"
        GOTO 2999
    ENDIF
    ALLOCATE (mpi_line_out(nprocs, 2), STAT = errCode)
    IF (errCode /= 0) THEN
        str_attr = "mpi_line Memory Allocation Error!!!"
        GOTO 2999
    ENDIF

    Call  Assign_mpi_line   (myrank, nprocs, linenum_lim, mpi_line, mpi_line_out)
    
    mpi_idx = myrank + 1
    WRITE(*,'("MPI_IDX:",i3)') mpi_idx
    pixlim  = pixnum_lim
    linelim = (/mpi_line(mpi_idx,1), mpi_line(mpi_idx,2)/) ! [1,1644]
    CALL Assign_pixline (linelim, pixlim, errCode)
    IF ( errCode /= 0 ) GOTO 2999

    mpi_pix_out (1)= pixloc(first_pix)
    mpi_pix_out (2)= pixloc(last_pix)
    IF ( mpi_line_out(mpi_idx,1) /= lineloc(first_line) ) STOP
    IF ( mpi_line_out(mpi_idx,2) /= lineloc(last_line) ) STOP
    
    !
    !---  write the size of line & pix to MPI log file
    !        
    CALL GEMS_Share_MOD_MPI_preLog(myrank, mpi_pix_out, mpi_line_out, errCode)
    IF ( errCode /= 0 ) GOTO 2999

    WRITE(*,'(A)') '_________________________________'
    WRITE(*,'(A)') ' STEP2 : O3 fitting'
    WRITE(*,'(A)') '_________________________________'
  
    !-------------------------------------------------------------------------
    !-------- 4. Fitting_Process
    !-------------------------------------------------------------------------
    CALL GEMS_O3P_SUB4_fitting_process(errCode)
    IF ( errCode /= 0 ) THEN
        print*, "O3 Fitting Error !!! "
        GOTO 1999  
    ENDIF
print*,'GEMS_O3P_SUB4_fitting_process DONE STOPPED.'
stop
    !-------------------------------------------------------------------------
    !-------- . MPI Processing
    !-------------------------------------------------------------------------
    
    CALL GEMS_Share_MOD_MPI_DataSumProcess(myrank, mpi_pix_out, mpi_line_out, L2O3P_wr, errCode)
    IF ( errCode /= 0 ) THEN
        print*, "MPI Running Error !!! "
        GOTO 1999  
    ENDIF

    !
    !-- Close MPI LOG File
    !
    CALL GEMS_Share_MOD_MPI_LogClose

    !-------------------------------------------------------------------------
    !-------- . MPI Finalize
    !-------------------------------------------------------------------------
    CALL GEMS_Share_MOD_MPI_Fin(errCode)
    IF ( errCode /= 0 ) THEN
        PRINT*, "MPI Finalize Error !!! "
        GOTO 1999  
    ENDIF

    !-------------------------------------------------------------------------
    !-------- 5. OutPut Processing
    !-------------------------------------------------------------------------
    IF ( myrank == 0 ) THEN

    WRITE(*,'(A)') '_________________________________'
    WRITE(*,'(A)') ' STEP3 : Saving Results to *.h5'
    WRITE(*,'(A)') '_________________________________'
    
     CALL GEMS_O3P_SUB5_Output(L2O3P_wr, errCode)
     IF ( errCode /= 0 ) GOTO 1999
    ENDIF

    !-------------------------------------------------------------------------
    !-------- 6. Close
    !-------------------------------------------------------------------------
    CALL GEMS_O3P_SUB6_Proc_close()
    errCode = GEMS_Share_MOD_L2O3P_MemDeAlloc2(L2O3P_wr)
    IF ( errCode /= 0 ) THEN
        print*, "Error to memory deallocation of L2 O3P Output Variable "
    ENDIF

    !************************************************
    !******** O3P MAIN PROGRAM END as Normal ******
    !************************************************
    !--- 공통 모듈 사용 후 작업
    CALL GEMS_O3P_MOD_share_close()

    CALL GEMS_Share_MOD_log(llvl,"---- End GEMS O3P Main ----" ,LOGMSG) 

    print*, "End O3P Program OK!"
    STOP

1999 CONTINUE
    !************************************************
    !******** O3P MAIN PROGRAM END as Error  ******
    !************************************************
    !--- 공통 모듈 사용 후 작업
    IF( errCode /= -9001 ) THEN
        IF( errCode == -9002 ) THEN
            CALL GEMS_Share_deAllocMem4Rad
        ELSE
            CALL GEMS_O3P_SUB6_Proc_close()
            errCode = GEMS_Share_MOD_L2O3P_MemDeAlloc2(L2O3P_wr)
            IF ( errCode /= 0 ) THEN
                print*, "Error to memory deallocation of L2 O3P Output Variable "
            ENDIF
            CALL GEMS_O3P_MOD_share_close()
        ENDIF
    ENDIF

    ! ------------------------------------
    ! Write END_OF_RUN message to log file
    ! ------------------------------------
    SELECT CASE ( errCode )
    CASE ( pge_errstat_ok )
       ! ----------------------------------------------------------------
       ! PGE execution completed successfully. All is well.
       ! ----------------------------------------------------------------
       errstat = OMI_SMF_setmsg ( OMSAO_S_ENDOFRUN, '', modulename, 0 )
       print * , '!!!! algorithm  ended without any error'
       STOP
    CASE ( pge_errstat_warning )
       ! ----------------------------------------------------------------
       ! PGE execution raised non-terminal warnings. Nothing serious, we
       ! hope, so execution completed but with a non-zero exit status.
       ! ----------------------------------------------------------------
       errstat = OMI_SMF_setmsg ( OMSAO_W_ENDOFRUN, '', modulename, 0 )
    CASE ( pge_errstat_error )
       ! ----------------------------------------------------------------
       ! PGE execution encountered an error that lead to termination.
       ! ----------------------------------------------------------------
       errstat = OMI_SMF_setmsg ( OMSAO_E_ENDOFRUN, '', modulename, 0 )
    CASE DEFAULT
       ! ----------------------------------------------------------------
       ! If we ever reach here, then PGE_ERRSTAT has been set to a funny
       ! value. This should never happen, but we buffer this case anyway.
       ! ----------------------------------------------------------------
       errstat = OMI_SMF_setmsg ( OMSAO_U_ENDOFRUN, '', modulename, 0 )
    END SELECT

2999 CONTINUE
     print*, "Error O3P Program! => ", "[", str_attr, "]"
     STOP


END PROGRAM SYNT_O3P_Main


!
!--- SUBROUTINE
!
SUBROUTINE  Assign_mpi_line  (myrank,nprocs,linelim, mpi_line, mpi_line_out)
 USE GEMS_O3P_gemsdata_module, ONLY: gems_ny,  nybin
 IMPLICIT NONE

 ! IN/OUT variables
 INTEGER, INTENT (IN)            :: myrank,nprocs
 INTEGER, INTENT (IN)            :: linelim(2)
 INTEGER, INTENT(OUT)            :: mpi_line(nprocs,2), mpi_line_out(nprocs, 2)
 ! local variables
 INTEGER :: i, j, line1, line2
 INTEGER :: offline, nline

 offline = linelim(1) -1
 nline   = (linelim(2) - linelim(1) +1)/nybin
 line1  = 1
 mpi_line(1,1) = 1
 DO i = 1, nprocs
    line2 = line1 + INT (nline*1.0/nprocs) -1
    IF ( i == nprocs)  line2 = nline
    mpi_line(i,1) = (line1-1)*nybin +1 +offline
    mpi_line(i,2)  = line2*nybin+offline
    mpi_line_out(i,1) = line1 + INT(offline/nybin) 
    mpi_line_out(i,2) = line2 + INT(offline/nybin)
    IF (myrank == 0 ) THEN 
        WRITE(*,'(A10,"[",i5,i5,"]",A10,"[",i5,i5,"]")') 'ORI_LINE', mpi_line(i,:), ' SAVE_LINE:[',mpi_line_out(i,:)
    ENDIF
    line1 = line2 + 1
 ENDDO
   
 RETURN
END SUBROUTINE Assign_mpi_line


!
!--- SUBROUTINE
!
SUBROUTINE Assign_pixline ( linenum_lim, pixnum_lim , pge_error_status )
 
 USE OMSAO_errstat_module
 USE GEMS_O3P_gemsdata_module, ONLY: do_xbin, do_ybin, nxbin, nybin, &
                                     offline, first_pix, last_pix, first_line,last_line,&
                                     lineloc, pixloc ! index to be saved in GEMS_final_output 
 
 IMPLICIT NONE

 !INPUT/OUTPUT variables
 INTEGER, INTENT(IN), DIMENSION(2)    :: pixnum_lim, linenum_lim
 INTEGER, INTENT(OUT) :: pge_error_status
 
 !local variables
 INTEGER :: i
 REAL :: currtrack , currline
  ! ------------------------------------------------
  ! Check for consistency of pixel limits to process
  ! -----------------------------------------------
  pge_error_status = pge_errstat_ok

  IF (do_ybin .AND. nybin > 1)  THEN
     i = linenum_lim(2)-linenum_lim(1) + 1
     IF (MOD(i, nybin) /= 0) THEN
            stop
     ENDIF
  ENDIF
  
  first_pix  = CEILING(1.0 * pixnum_lim(1) / nxbin)
  last_pix   = NINT(1.0 * pixnum_lim(2) / nxbin )
  offline    = linenum_lim(1) -1
  first_line = 1 
  last_line  = NINT (1.0 *(linenum_lim(2) - linenum_lim(1) +1) / nybin)
  IF (nybin == 1) last_line  =  linenum_lim(2) - linenum_lim(1) + 1

  do i = first_line, last_line
    currline =   (i-1)*nybin + offline+1
    lineloc(i) = i +  INT(offline/nybin)
    print * , '[line]',i, currline, lineloc(i)
  enddo
  do i = first_pix, last_pix
    currtrack = (i-1)*nxbin + 1
    pixloc(i) = i 
   print * , '[pix ]',i, currtrack, pixloc(i)
  enddo
 RETURN
END SUBROUTINE Assign_pixline


!
!--- SUBROUTINE MPI Data Sum Processing
!
SUBROUTINE GEMS_Share_MOD_MPI_DataSumProcess(myrank, mpi_pix, mpi_line, L2O3P_wr, retCode)
 USE Share_l2_o3p_mod_write_read
 USE Share_MOD_MPI

 IMPLICIT NONE
 !--------------------------------------------
 !INPUT/OUTPUT variables
 !-------------------------------------------
    INTEGER(KIND=4), INTENT(IN)     :: myrank
    INTEGER(KIND=4), INTENT(IN)     :: mpi_pix(:)
    INTEGER(KIND=4), INTENT(IN)     :: mpi_line(:,:)
    TYPE(L2_o3p),    INTENT(OUT)    :: L2O3P_wr
    INTEGER(KIND=4), INTENT(OUT)    :: retCode
 !------------------------------------------------    
 !local variables
 !------------------------------------------------
    INTEGER(KIND=4)                 :: errCode = 0
    INTEGER(KIND=8)             :: iAddr
    !
    !
    !-- Get the memory of L2 Output Variable for write for MPI
    !
    errCode = GEMS_Share_MOD_L2O3P_MemAlloc2(L2O3P_wr)
    IF ( errCode /= 0 ) THEN
        print*, "Error to memory allocation of L2 O3P Output Variable "
        retCode = errCode
        RETURN 
    ENDIF
    
    iAddr = loc(L2O3P_wr%O3P_dfld%ColAmtOz)
    write(unit_num, '(A2, A, I12, A, I12, A, /)'),                              &
            ' ',                                                                &
            ' iAddr after allocating L2O3P_wr%O3P_dfld%ColAmtOz memory   = [',  &
            iAddr, ']', SIZE(L2O3P_wr%O3P_dfld%ColAmtOz), '+'

    errCode = GEMS_Share_MOD_L2O3P_MemInit2(L2O3P_wr)
    IF ( errCode /= 0 ) THEN
        print*, "Error to memory initialization of L2 O3P Output Variable "
        retCode = errCode
        RETURN  
    ENDIF
    
    !
    !--- Sum the value of Data Fields in L2 Ozone Ouput
    !
    CALL GEMS_Share_MOD_MPI_Process(myrank,                        &
                                    gds_L2O3P%o3p_dfld%AvgK,       &
                                    L2O3P_wr%o3p_dfld%AvgK,        &
                                    mpi_pix, mpi_line, errCode)
    IF ( errCode /= 0 ) THEN
        PRINT*, "[AvgK] MPI Processing Error !!!"
        retCode = errCode
        RETURN  
    ENDIF

    CALL GEMS_Share_MOD_MPI_Process(myrank,                        &
                                    gds_L2O3P%o3p_gfld%TEMP,       &
                                    L2O3P_wr%o3p_gfld%TEMP,        &
                                    mpi_pix, mpi_line, errCode)
    IF ( errCode /= 0 ) THEN
        PRINT*, "[TEMP] MPI Processing Error !!!"
        retCode = errCode
        RETURN  
    ENDIF

    CALL GEMS_Share_MOD_MPI_Process(myrank,                        &
                                    gds_L2O3P%o3p_dfld%OzRet,      &
                                    L2O3P_wr%o3p_dfld%OzRet,       &
                                    mpi_pix, mpi_line, errCode)
    IF ( errCode /= 0 ) THEN
        PRINT*, "[OzRet] MPI Processing Error !!!"
        retCode = errCode
        RETURN  
    ENDIF

    CALL GEMS_Share_MOD_MPI_Process(myrank,                        &
                                    gds_L2O3P%o3p_dfld%OzAP,       &
                                    L2O3P_wr%o3p_dfld%OzAP,        &
                                    mpi_pix, mpi_line, errCode)
    IF ( errCode /= 0 ) THEN
        PRINT*, "[OzAP] MPI Processing Error !!!"
        retCode = errCode
        RETURN  
    ENDIF

    CALL GEMS_Share_MOD_MPI_Process(myrank,                        &
                                    gds_L2O3P%o3p_dfld%OzApErr,    &
                                    L2O3P_wr%o3p_dfld%OzApErr,     &
                                    mpi_pix, mpi_line, errCode)
    IF ( errCode /= 0 ) THEN
        PRINT*, "[OzApErr] MPI Processing Error !!!"
        retCode = errCode
        RETURN  
    ENDIF

    CALL GEMS_Share_MOD_MPI_Process(myrank,                        &
                                    gds_L2O3P%o3p_dfld%OzNErr,     &
                                    L2O3P_wr%o3p_dfld%OzNErr,      &
                                    mpi_pix, mpi_line, errCode)
    IF ( errCode /= 0 ) THEN
        PRINT*, "[OzNErr] MPI Processing Error !!!"
        retCode = errCode
        RETURN  
    ENDIF

    CALL GEMS_Share_MOD_MPI_Process(myrank,                        &
                                    gds_L2O3P%o3p_dfld%OzSMErr,    &
                                    L2O3P_wr%o3p_dfld%OzSMErr,     &
                                    mpi_pix, mpi_line, errCode)
    IF ( errCode /= 0 ) THEN
        PRINT*, "[OzSMErr] MPI Processing Error !!!"
        retCode = errCode
        RETURN  
    ENDIF

    CALL GEMS_Share_MOD_MPI_Process(myrank,                        &
                                    gds_L2O3P%o3p_dfld%ColAmtOz,   &
                                    L2O3P_wr%o3p_dfld%ColAmtOz,    &
                                    mpi_pix, mpi_line, errCode)
    IF ( errCode /= 0 ) THEN
        PRINT*, "[ColAmtOz] MPI Processing Error !!!"
        retCode = errCode
        RETURN  
    ENDIF

    CALL GEMS_Share_MOD_MPI_Process(myrank,                        &
                                    gds_L2O3P%o3p_dfld%DFS,        &
                                    L2O3P_wr%o3p_dfld%DFS,         &
                                    mpi_pix, mpi_line, errCode)
    IF ( errCode /= 0 ) THEN
        PRINT*, "[DFS] MPI Processing Error !!!"
        retCode = errCode
        RETURN  
    ENDIF

    CALL GEMS_Share_MOD_MPI_Process(myrank,                        &
                                    gds_L2O3P%o3p_dfld%CldP,       &
                                    L2O3P_wr%o3p_dfld%CldP,        &
                                    mpi_pix, mpi_line, errCode)
    IF ( errCode /= 0 ) THEN
        PRINT*, "[CldP] MPI Processing Error !!!"
        retCode = errCode
        RETURN  
    ENDIF

    CALL GEMS_Share_MOD_MPI_Process(myrank,                        &
                                    gds_L2O3P%o3p_dfld%EffCldFrac, &
                                    L2O3P_wr%o3p_dfld%EffCldFrac,  &
                                    mpi_pix, mpi_line, errCode)
    IF ( errCode /= 0 ) THEN
        PRINT*, "[EffCldFrac] MPI Processing Error !!!"
        retCode = errCode
        RETURN  
    ENDIF

    CALL GEMS_Share_MOD_MPI_Process(myrank,                        &
                                    gds_L2O3P%o3p_dfld%SfcAlb,     &
                                    L2O3P_wr%o3p_dfld%SfcAlb,      &
                                    mpi_pix, mpi_line, errCode)
    IF ( errCode /= 0 ) THEN
        PRINT*, "[SfcAlb] MPI Processing Error !!!"
        retCode = errCode
        RETURN  
    ENDIF

!------ Not Calculated Product
!    CALL GEMS_Share_MOD_MPI_Process(myrank,                        &
!                                    gds_L2O3P%o3p_dfld%TerrP,      &
!                                    L2O3P_wr%o3p_dfld%TerrP,       &
!                                    mpi_pix, mpi_line, errCode)
!    IF ( errCode /= 0 ) THEN
!        PRINT*, "[TerrP] MPI Processing Error !!!"
!        retCode = errCode
!        RETURN  
!    ENDIF

    CALL GEMS_Share_MOD_MPI_Process(myrank,                        &
                                    gds_L2O3P%o3p_dfld%PrsQFlag,   &
                                    L2O3P_wr%o3p_dfld%PrsQFlag,    &
                                    mpi_pix, mpi_line, errCode)
    IF ( errCode /= 0 ) THEN
        PRINT*, "[PrsQFlag] MPI Processing Error !!!"
        retCode = errCode
        RETURN  
    ENDIF

    CALL GEMS_Share_MOD_MPI_Process(myrank,                        &
                                    gds_L2O3P%o3p_dfld%NumIter,    &
                                    L2O3P_wr%o3p_dfld%NumIter,     &
                                    mpi_pix, mpi_line, errCode)
    IF ( errCode /= 0 ) THEN
        PRINT*, "[NumIter MPI Processing Error !!!"
        retCode = errCode
        RETURN  
    ENDIF

    CALL GEMS_Share_MOD_MPI_Process(myrank,                        &
                                    gds_L2O3P%o3p_dfld%AvgRes,     &
                                    L2O3P_wr%o3p_dfld%AvgRes,      &
                                    mpi_pix, mpi_line, errCode)
    IF ( errCode /= 0 ) THEN
        PRINT*, "[AvgRes] MPI Processing Error !!!"
        retCode = errCode
        RETURN  
    ENDIF

    CALL GEMS_Share_MOD_MPI_Process(myrank,                        &
                                    gds_L2O3P%o3p_dfld%RMS,        &
                                    L2O3P_wr%o3p_dfld%RMS,         &
                                    mpi_pix, mpi_line, errCode)
    IF ( errCode /= 0 ) THEN
        PRINT*, "[RMS] MPI Processing Error !!!"
        retCode = errCode
        RETURN  
    ENDIF

    !
    !--- Sum the value of Geolocation Fields in L2 Ozone Ouput
    !
    CALL GEMS_Share_MOD_MPI_Process(myrank,                        &
                                    gds_L2O3P%o3p_gfld%Alt,        &
                                    L2O3P_wr%o3p_gfld%Alt,         &
                                    mpi_pix, mpi_line, errCode)
    IF ( errCode /= 0 ) THEN
        PRINT*, "MPI Processing Error !!!"
        retCode = errCode
        RETURN  
    ENDIF

    CALL GEMS_Share_MOD_MPI_Process(myrank,                        &
                                    gds_L2O3P%o3p_gfld%P,          &
                                    L2O3P_wr%o3p_gfld%P,           &
                                    mpi_pix, mpi_line, errCode)
    IF ( errCode /= 0 ) THEN
        PRINT*, "[P] MPI Processing Error !!!"
        retCode = errCode
        RETURN  
    ENDIF

    CALL GEMS_Share_MOD_MPI_Process(myrank,                        &
                                    gds_L2O3P%o3p_gfld%LAT,        &
                                    L2O3P_wr%o3p_gfld%LAT,         &
                                    mpi_pix, mpi_line, errCode)
    
  
    IF ( errCode /= 0 ) THEN
        PRINT*, "[LAT] MPI Processing Error !!!"
        retCode = errCode
        RETURN  
    ENDIF

    CALL GEMS_Share_MOD_MPI_Process(myrank,                        &
                                    gds_L2O3P%o3p_gfld%LON,        &
                                    L2O3P_wr%o3p_gfld%LON,         &
                                    mpi_pix, mpi_line, errCode)
    IF ( errCode /= 0 ) THEN
        PRINT*, "[LON] MPI Processing Error !!!"
        retCode = errCode
        RETURN  
    ENDIF

    CALL GEMS_Share_MOD_MPI_Process(myrank,                        &
                                    gds_L2O3P%o3p_gfld%SolZenAng,  &
                                    L2O3P_wr%o3p_gfld%SolZenAng,   &
                                    mpi_pix, mpi_line, errCode)
    IF ( errCode /= 0 ) THEN
        PRINT*, "[SolZenAng] MPI Processing Error !!!"
        retCode = errCode
        RETURN  
    ENDIF

    CALL GEMS_Share_MOD_MPI_Process(myrank,                        &
                                    gds_L2O3P%o3p_gfld%ViewZenAng, &
                                    L2O3P_wr%o3p_gfld%ViewZenAng,  &
                                    mpi_pix, mpi_line, errCode)
    IF ( errCode /= 0 ) THEN
        PRINT*, "[ViewZenAng MPI Processing Error !!!"
        retCode = errCode
        RETURN  
    ENDIF

    CALL GEMS_Share_MOD_MPI_Process(myrank,                        &
                                    gds_L2O3P%o3p_gfld%RelAziAng,  &
                                    L2O3P_wr%o3p_gfld%RelAziAng,   &
                                    mpi_pix, mpi_line, errCode)
    IF ( errCode /= 0 ) THEN
        PRINT*, "[RelAziAng MPI Processing Error !!!"
        retCode = errCode
        RETURN  
    ENDIF
    
    CALL GEMS_Share_MOD_MPI_Process(myrank,                        &
                                    gds_L2O3P%o3p_gfld%Time,       &
                                    L2O3P_wr%o3p_gfld%Time,        &
                                    mpi_pix, mpi_line, errCode)
    IF ( errCode /= 0 ) THEN
        PRINT*, "[Time] MPI Processing Error !!!"
        retCode = errCode
        RETURN  
    ENDIF
   
    ! geun added PIX output
    CALL GEMS_Share_MOD_MPI_Process(myrank,                        &
                                    gds_L2O3P%o3p_gfld%Pix_tmp,       &
                                    L2O3P_wr%o3p_gfld%Pix_tmp,        &
                                    mpi_pix, mpi_line, errCode)
    IF ( errCode /= 0 ) THEN
        PRINT*, "[Pix] MPI Processing Error !!!"
        retCode = errCode
        RETURN  
    ENDIF


    ! geun added Line output
    CALL GEMS_Share_MOD_MPI_Process(myrank,                        &
                                    gds_L2O3P%o3p_gfld%Line_tmp,       &
                                    L2O3P_wr%o3p_gfld%Line_tmp,        &
                                    mpi_pix, mpi_line, errCode)
    IF ( errCode /= 0 ) THEN
        PRINT*, "[Line] MPI Processing Error !!!"
        retCode = errCode
        RETURN  
    ENDIF

    retCode = errCode
    RETURN    

END SUBROUTINE GEMS_Share_MOD_MPI_DataSumProcess 
