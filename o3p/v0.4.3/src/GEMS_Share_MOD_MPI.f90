!-------------------------------------------------------------------------------
!+

MODULE Share_MOD_MPI
!-------------------------------------------------------------------------------
!+Description: 
!     . Module to handle Log & MPI function
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.0     2016.06.16 Fisrt Code (R&D, SaeaSoft Co., Ltd.) 
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Modules:
! USE O3P_MOD_dummy

 USE Share_MOD_Log


!**********************************************************!
 IMPLICIT NONE

!
! GEMS MPI Include
!
include 'mpif.h'

!
!--- Interface definition
!
INTERFACE GEMS_Share_MOD_MPI_Process
    MODULE PROCEDURE   GEMS_Share_MOD_MPI_Process_3D
END INTERFACE

!--
INTEGER(KIND=4)             :: istatus(MPI_STATUS_SIZE)
INTEGER(KIND=4)             :: gci_nprocs

CHARACTER(len=MPI_MAX_LIBRARY_VERSION_STRING) version

INTEGER(KIND=4)             :: unit_num


CONTAINS

!
!---
!
SUBROUTINE GEMS_Share_MOD_MPI_Start(myrank, nprocs, retCode)

! argument parameter
    INTEGER(KIND=4), INTENT(OUT)    :: myrank
    INTEGER(KIND=4), INTENT(OUT)    :: nprocs
    INTEGER(KIND=4), INTENT(OUT)    :: retCode

! local variables
    INTEGER(KIND=4)             :: errCode = 0
    INTEGER(KIND=4)             :: i, j

    !-------------------------------------------------------------------------
    !-------- . MPI Start
    !-------------------------------------------------------------------------
    call MPI_INIT(errCode)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, gci_nprocs, errCode)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, errCode)
    IF( errCode /= 0 ) THEN
        PRINT*, 'Error !!! [', errCode, ']'
        retCode = errCode
        STOP
    ENDIF

    nprocs = gci_nprocs
    retCode = errCode
    RETURN
END SUBROUTINE GEMS_Share_MOD_MPI_Start


!
!---
!
SUBROUTINE GEMS_Share_MOD_MPI_LogOpen(myrank, log_file_name, retCode)

! argument parameter
    INTEGER(KIND=4),    INTENT(IN)      :: myrank
    CHARACTER(LEN=120), INTENT(IN)      :: log_file_name
    INTEGER(KIND=4),    INTENT(OUT)     :: retCode

! local variables
    CHARACTER(LEN=120)                  :: mpi_log_file
    INTEGER(KIND=4)                     :: errCode = 0

    !
    !-- Open MPI LOG File 
    !
    unit_num = 8
    WRITE(mpi_log_file, '(A, A, I0.3, A)'), '../../log/', log_file_name, myrank, '.log'
    OPEN(UNIT=unit_num, FILE=mpi_log_file, ACCESS='append', IOSTAT=errCode)
    IF ( errCode /= 0 ) THEN
        retCode = errCode
        PRINT*, "Error!, did not open [ ", TRIM(mpi_log_file), "] => ", errCode
        STOP
    ENDIF

    RETURN
END SUBROUTINE GEMS_Share_MOD_MPI_LogOpen

!
!---
!
SUBROUTINE GEMS_Share_MOD_MPI_LogClose()

    !
    !-- Close MPI LOG File
    !
    CLOSE(UNIT=unit_num)

    RETURN
END SUBROUTINE GEMS_Share_MOD_MPI_LogClose


!
!---
!
SUBROUTINE GEMS_Share_MOD_MPI_preLog(myrank, ix, iyy, retCode)

! argument parameter
    INTEGER(KIND=4),    INTENT(IN)  :: myrank
    INTEGER(KIND=4),    INTENT(IN)  :: ix(:)
    INTEGER(KIND=4),    INTENT(IN)  :: iyy(:, :)
    INTEGER(KIND=4),    INTENT(OUT) :: retCode

! local variables
    INTEGER(KIND=4)                 :: errCode = 0
    INTEGER                         :: n_mpi_tmp

    INTEGER(KIND=4)                 :: i, j

    retCode = errCode
    !
    !-- write the position of line & pix to MPI Log
    !
    IF ( myrank == 0 ) THEN
        DO i=1, gci_nprocs
            
            n_mpi_tmp = i
            WRITE(unit_num, '(A15, I7, 2(A20,I7,A1))'),      &
                    '###  myrank = [',  n_mpi_tmp-1 ,        &
                    '], [iy1]=[',  iyy(n_mpi_tmp,1), ']',    &
                    ', [iy2]=[',    iyy(n_mpi_tmp,2), ']'    
            WRITE(unit_num, '(A15, I7, 2(A20,I7,A1))'),      &
                    '###  myrank = [',  n_mpi_tmp-1,         &
                    '],  [ix1]=[',   ix(1), ']',             &
                    ' , [ix2]=[',    ix(2), ']'

        ENDDO
    ENDIF

    RETURN
END SUBROUTINE GEMS_Share_MOD_MPI_preLog


!
!---
!
SUBROUTINE GEMS_Share_MOD_MPI_Process_3D(myrank, pSend, pRecv, ix, iyy, retCode)

! argument parameter
    INTEGER(KIND=4),       INTENT(IN)   :: myrank
    REAL(KIND=8), POINTER, INTENT(IN)   :: pSend(:, :, :)
    REAL(KIND=8), POINTER, INTENT(IN)   :: pRecv(:, :, :)
    INTEGER(KIND=4),       INTENT(IN)   :: ix(:)
    INTEGER(KIND=4),       INTENT(IN)   :: iyy(:, :)
    INTEGER(KIND=4),       INTENT(OUT)  :: retCode

!    REAL(KIND=4), POINTER       :: pfSend(:, :), pfRecv(:, :)
!    REAL(KIND=4), ALLOCATABLE   :: fSend(:, :), fRecv(:, :)

! local variables
    INTEGER(KIND=4)             :: errCode = 0
    CHARACTER(LEN=120)          :: w_fmt

    INTEGER(KIND=4)             :: iscnt
    INTEGER(KIND=4)             :: iscnt2
    INTEGER, POINTER            :: ircnt(:), idisp(:)

    INTEGER(KIND=4)             :: n_mpi
    INTEGER(KIND=4)             :: iy(2)
    INTEGER(KIND=4)             :: ix_first  
    INTEGER(KIND=4)             :: ix_last  
    INTEGER(KIND=4)             :: iy_size  

    INTEGER(KIND=4)             :: ixTotSize
    INTEGER(KIND=4)             :: iyTotSize 

    INTEGER(KIND=4)             :: i,j

    n_mpi = myrank+1
    iy    = iyy(n_mpi, 1:2)

    retCode = errCode
    IF ( gci_nprocs <= 0 ) THEN
        PRINT*, 'MPI Process Total Count Input Error!!!'
        STOP
    ENDIF

    !
    !--- point the 3D pointer variable to 3D target variable
    !
!    pSend=>gds_L2O3P%o3p_dfld%ColAmtOz(:, :, :)
!    pRecv=>L2O3P_wr%o3p_dfld%ColAmtOz(:, :, :)
    WRITE(unit_num, '(A15, I7, 2(A20,3I7,A1))'),            &
            '###  myrank = [',  myrank,                     &
            '], shape(pSend)=[',   shape(pSend), ']',     &
            '   shape(pRecv)=[',   shape(pRecv), ']'

    ixTotSize = SIZE(pSend(1,:,1))
    iyTotSize = SIZE(pSend(1,1,:))

    WRITE(unit_num, '(A15, I7, 2(A32,I7,A1))'),                &
            '###  myrank = [',  myrank,                        &
            '   SIZE(pSend(1,:,1))=ixTotSize=[',  ixTotSize, ']',  &  
            '], SIZE(pSend(1,1,:))=iyTotSize=[',  iyTotSize, ']'
    
    !
    !-- MPI variable 
    !
    IF (ASSOCIATED(ircnt)) DEALLOCATE(ircnt)
    ALLOCATE (ircnt(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Receive Count Allocation Error!!!'
            RETURN
         ENDIF
    IF (ASSOCIATED(idisp)) DEALLOCATE(idisp)
    ALLOCATE (idisp(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Display Count Allocation Error!!!'
            RETURN
         ENDIF

    !
    !-- Set send count, receive count, display position for MPI
    !
    iscnt   = SIZE(pSend(:, :, iy(1):iy(2)))
  
    IF ( myrank == 0 ) THEN
        iscnt2   = 0
        DO i=1, gci_nprocs
     
            !
            !-- set display position for MPI
            !
            IF ( i == 1 ) THEN
                idisp(i) = 0
            ELSE
                idisp(i) = iscnt2
            ENDIF

!            n_mpi_tmp = i
!            iyy1  = mpi_line(n_mpi_tmp,1)
!            iyy2  = mpi_line(n_mpi_tmp,2)

            iscnt2   = iscnt2 + SIZE(pSend(:, :, iyy(i,1):iyy(i,2)))
            !
            !-- set receive count for MPI
            !
            ircnt(i) = SIZE(pSend(:, :, iyy(i,1):iyy(i,2)))
        ENDDO
    ENDIF

    !
    !-- Write MPI Log
    !
    WRITE(unit_num, '(A15, I7, A12, I9, A1)'),                                      &
            '###  myrank = [', myrank, '], iscnt = [', iscnt, ']'
    
    IF ( myrank == 0 ) THEN
        WRITE(unit_num, '(A15, I7, A1)'),                                           &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                                           &
                '      ircnt  = ['
        DO i=1, SIZE(ircnt)
            WRITE(unit_num, '(A14, I9)'), ' ', ircnt(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
        WRITE(unit_num, '(A15, I7, A1)'),                                           &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                                           &
                '      idisp  = ['
        DO i=1, SIZE(idisp)
            WRITE(unit_num, '(A14, I9)'), ' ', idisp(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
    ENDIF

    WRITE(unit_num, '(A3)') "---"

    !
    !-- Gather the output From other processes
    !
    CALL MPI_GATHERV(pSend(:, :, iy(1)), iscnt, MPI_REAL,            &
                     pRecv,            ircnt, idisp, MPI_REAL,     &
                     0, MPI_COMM_WORLD, errCode)

    retCode = errCode

    !
    !--- set pixel size
    !
    ix_first = 1            
    ix_last  = ixTotSize    

    !
    !
    !-- Write MPI Log
    !
    iy_size = iy(2) - iy(1)
    WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', SIZE(pSend(:,1,1))*(iy_size+1),'(F10.3))'

    IF (myrank==0) THEN
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(:, i, iy(1):iy(2)))
        ENDDO

        DO j=1, gci_nprocs

!            n_mpi_tmp = j
!            iyy1 = mpi_line(n_mpi_tmp,1)
!            iyy2 = mpi_line(n_mpi_tmp,2)

            iy_size = iyy(i,2) - iyy(i,1)
            WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(F10.3))'
            DO i=ix_first, ix_last
                WRITE(unit_num, w_fmt),                                                 &
                        'myrank = [', j-1, '], irecv =',                                &
                        pRecv(:, i, iyy(i,1):iyy(i,2))
            ENDDO
        ENDDO
    ELSE
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(:, i, iy(1):iy(2)))
        ENDDO
    ENDIF

    RETURN
END SUBROUTINE GEMS_Share_MOD_MPI_Process_3D


!
!--- 
!
SUBROUTINE GEMS_Share_MOD_MPI_Fin(retCode)

! argument parameter
    INTEGER(KIND=4),    INTENT(OUT)     :: retCode

    CALL MPI_FINALIZE(retCode)

    RETURN
END SUBROUTINE GEMS_Share_MOD_MPI_Fin

END MODULE Share_MOD_MPI