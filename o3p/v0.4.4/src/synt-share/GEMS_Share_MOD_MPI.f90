!-------------------------------------------------------------------------------
!+

MODULE Share_MOD_MPI
!-------------------------------------------------------------------------------
!+Description: 
!     . Module to handle MPI function and write to Log
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
    MODULE PROCEDURE   GEMS_Share_MOD_MPI_Process_1D_I1
    MODULE PROCEDURE   GEMS_Share_MOD_MPI_Process_1D_I2
    MODULE PROCEDURE   GEMS_Share_MOD_MPI_Process_1D_I4
    MODULE PROCEDURE   GEMS_Share_MOD_MPI_Process_1D_R4
    MODULE PROCEDURE   GEMS_Share_MOD_MPI_Process_1D_R8
    MODULE PROCEDURE   GEMS_Share_MOD_MPI_Process_2D_I1
    MODULE PROCEDURE   GEMS_Share_MOD_MPI_Process_2D_I2
    MODULE PROCEDURE   GEMS_Share_MOD_MPI_Process_2D_I4
    MODULE PROCEDURE   GEMS_Share_MOD_MPI_Process_2D_R4
    MODULE PROCEDURE   GEMS_Share_MOD_MPI_Process_2D_R8
    MODULE PROCEDURE   GEMS_Share_MOD_MPI_Process_3D_I1
    MODULE PROCEDURE   GEMS_Share_MOD_MPI_Process_3D_I2
    MODULE PROCEDURE   GEMS_Share_MOD_MPI_Process_3D_I4
    MODULE PROCEDURE   GEMS_Share_MOD_MPI_Process_3D_R4
    MODULE PROCEDURE   GEMS_Share_MOD_MPI_Process_3D_R8
    MODULE PROCEDURE   GEMS_Share_MOD_MPI_Process_4D_I1
    MODULE PROCEDURE   GEMS_Share_MOD_MPI_Process_4D_I2
    MODULE PROCEDURE   GEMS_Share_MOD_MPI_Process_4D_I4
    MODULE PROCEDURE   GEMS_Share_MOD_MPI_Process_4D_R4
    MODULE PROCEDURE   GEMS_Share_MOD_MPI_Process_4D_R8
END INTERFACE

!--
INTEGER(KIND=4)             :: gci_nprocs               ! 병렬처리 총 개수

INTEGER(KIND=4)             :: istatus(MPI_STATUS_SIZE)

INTEGER(KIND=4)             :: unit_num                 ! 로그파일의 출력디바이스

CHARACTER(len=MPI_MAX_LIBRARY_VERSION_STRING) version

CHARACTER(LEN=120)          :: mpi_log_dir              ! 로그파일의 디렉토리 위치


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

print*,'MPI START Done'
    !-------------------------------------------------------------------------
    !-------- . MPI Start
    !-------------------------------------------------------------------------
    call MPI_INIT(errCode)
print*,'MPI INIT Done'
    call MPI_COMM_SIZE(MPI_COMM_WORLD, gci_nprocs, errCode)
print*,'MPI COMM SIZE Done'
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, errCode)
print*,'MPI COMM Rank Done'
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
    CHARACTER(LEN=256)                  :: mpi_log_file
    INTEGER(KIND=4)                     :: errCode = 0

    mpi_log_dir = '../../log/'

    !
    !-- Open MPI LOG File 
    !
    unit_num = 8
    WRITE(mpi_log_file, '(A, A, I0.3, A)'),                                 &
            TRIM(mpi_log_dir), TRIM(log_file_name), myrank, '.log'
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
            WRITE(unit_num, '(A15, I7, 2(A13,I7,A1))'),         &
                    '###  myrank = [',  n_mpi_tmp-1 ,           &
                    '], [iyy(1)]=[',  iyy(n_mpi_tmp,1), ']',    &
                    ',  [iyy(2)]=[',  iyy(n_mpi_tmp,2), ']'    
            WRITE(unit_num, '(A15, I7, 2(A13,I7,A1))'),         &
                    '###  myrank = [',  n_mpi_tmp-1,            &
                    '],  [ix(1)]=[',   ix(1), ']',              &
                    ' ,  [ix(2)]=[',   ix(2), ']'

        ENDDO
    ENDIF

    RETURN
END SUBROUTINE GEMS_Share_MOD_MPI_preLog

!
!--- 1Dim 8bit 정수형 변수를 병렬처리
!
SUBROUTINE GEMS_Share_MOD_MPI_Process_1D_I1(myrank, pSend, pRecv, ix, iyy, retCode)

! argument parameter
    INTEGER(KIND=4),          INTENT(IN)    :: myrank
    INTEGER(KIND=1), POINTER, INTENT(IN)    :: pSend(:)
    INTEGER(KIND=1), POINTER, INTENT(IN)    :: pRecv(:)
    INTEGER(KIND=4),          INTENT(IN)    :: ix(:)
    INTEGER(KIND=4),          INTENT(IN)    :: iyy(:, :)
    INTEGER(KIND=4),          INTENT(OUT)   :: retCode

! local variables
    INTEGER(KIND=4)             :: errCode = 0
    CHARACTER(LEN=120)          :: w_fmt

    INTEGER(KIND=4)             :: iscnt
    INTEGER(KIND=4)             :: iscnt2
    INTEGER(KIND=4),ALLOCATABLE :: ircnt(:), idisp(:)

    INTEGER(KIND=4)             :: n_mpi
    INTEGER(KIND=4)             :: n_mpi_tmp
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

    WRITE(unit_num, '(/, A15, I7, A, /)'),                 &
            '###  myrank = [',  myrank,                    &
            '] => 1Dim 8bit 정수형 변수를 병렬처리'

    !
    !--- point the 1D pointer variable to 1D target variable
    !
    WRITE(unit_num, '(A15, I7, 2(A17,I7,A1))'),            &
            '###  myrank = [',  myrank,                     &
            '], shape(pSend)=[',   shape(pSend), ']',       &
            '   shape(pRecv)=[',   shape(pRecv), ']'

    ixTotSize = SIZE(pSend(:))
    iyTotSize = SIZE(pSend(:))

    WRITE(unit_num, '(A15, I7, 2(A32,I7,A1))'),                 &
            '###  myrank = [',  myrank,                         &
            '], SIZE(pSend(:))=ixTotSize=[',  ixTotSize, ']',   &  
            '], SIZE(pSend(:))=iyTotSize=[',  iyTotSize, ']'
    
    !
    !-- MPI variable 
    !
    ALLOCATE (ircnt(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Receive Count Allocation Error!!!'
            RETURN
         ENDIF
    ALLOCATE (idisp(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Display Count Allocation Error!!!'
            RETURN
         ENDIF

    !
    !-- Set send count, receive count, display position for MPI
    !
    iscnt   = SIZE(pSend(iy(1):iy(2)))        ! 마스터 프로세스(rank=0)로 전송할 데이터 개수
  
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

            n_mpi_tmp = i
            iscnt2   = iscnt2 + SIZE(pSend(iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))

            !
            !-- set receive count for MPI
            !
            ircnt(i) = SIZE(pSend(iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))
        ENDDO
    ENDIF

    !
    !-- Write MPI Log
    !
    WRITE(unit_num, '(A15, I7, A12, I7, A1, A15, I7, A1, I7, A2)'),     &
            '###  myrank = [', myrank, '], iscnt = [', iscnt, ']',      &
            ' => SIZE(pSend(', iy(1), ':', iy(2), '))'
    
    IF ( myrank == 0 ) THEN
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      ircnt  = ['
        DO i=1, SIZE(ircnt)
            WRITE(unit_num, '(A14, I9)'), ' ', ircnt(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      idisp  = ['
        DO i=1, SIZE(idisp)
            WRITE(unit_num, '(A14, I9)'), ' ', idisp(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
    ENDIF

    WRITE(unit_num, '(A3, I5)') "---iy(1)", iy(1)

    !
    !-- Gather the output From other processes
    !
    CALL MPI_GATHERV(pSend(iy(1))   ,    iscnt, MPI_CHARACTER,            &
                     pRecv,              ircnt, idisp, MPI_CHARACTER,     &
                     0, MPI_COMM_WORLD, errCode)

    retCode = errCode

    !
    !--- set pixel size calculated in defined total pixel size
    !
    ix_first = 1            
    ix_last  = 1   ! ix(2)    

    !
    !
    !-- Write MPI Log
    !
    iy_size = iy(2) - iy(1)
    WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(I10))'

    IF (myrank==0) THEN
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(iy(1):iy(2)))
        ENDDO

        DO j=1, gci_nprocs

            n_mpi_tmp = j
            iy_size = iyy(n_mpi_tmp,2) - iyy(n_mpi_tmp,1)
            WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(I10))'
            DO i=ix_first, ix_last
                WRITE(unit_num, w_fmt),                                                 &
                        'myrank = [', n_mpi_tmp-1, '], irecv =',                        &
                        pRecv(iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2))
            ENDDO
        ENDDO
    ELSE
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(iy(1):iy(2)))
        ENDDO
    ENDIF

    DEALLOCATE(ircnt)
    DEALLOCATE(idisp)

    RETURN
END SUBROUTINE GEMS_Share_MOD_MPI_Process_1D_I1


!
!--- 1Dim 16bit 정수형 변수를 병렬처리
!
SUBROUTINE GEMS_Share_MOD_MPI_Process_1D_I2(myrank, pSend, pRecv, ix, iyy, retCode)

! argument parameter
    INTEGER(KIND=4),          INTENT(IN)    :: myrank
    INTEGER(KIND=2), POINTER, INTENT(IN)    :: pSend(:)
    INTEGER(KIND=2), POINTER, INTENT(IN)    :: pRecv(:)
    INTEGER(KIND=4),          INTENT(IN)    :: ix(:)
    INTEGER(KIND=4),          INTENT(IN)    :: iyy(:, :)
    INTEGER(KIND=4),          INTENT(OUT)   :: retCode

! local variables
    INTEGER(KIND=4)             :: errCode = 0
    CHARACTER(LEN=120)          :: w_fmt

    INTEGER(KIND=4)             :: iscnt
    INTEGER(KIND=4)             :: iscnt2
    INTEGER(KIND=4),ALLOCATABLE :: ircnt(:), idisp(:)

    INTEGER(KIND=4)             :: n_mpi
    INTEGER(KIND=4)             :: n_mpi_tmp
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

    WRITE(unit_num, '(/, A15, I7, A, /)'),                 &
            '###  myrank = [',  myrank,                    &
            '] => 1Dim 16bit 정수형 변수를 병렬처리'

    !
    !--- point the 1D pointer variable to 1D target variable
    !
    WRITE(unit_num, '(A15, I7, 2(A17,I7,A1))'),             &
            '###  myrank = [',  myrank,                     &
            '], shape(pSend)=[',   shape(pSend), ']',       &
            '   shape(pRecv)=[',   shape(pRecv), ']'

    ixTotSize = SIZE(pSend(:))
    iyTotSize = SIZE(pSend(:))

    WRITE(unit_num, '(A15, I7, 2(A32,I7,A1))'),                   &
            '###  myrank = [',  myrank,                           &
            '], SIZE(pSend(:))=ixTotSize=[',  ixTotSize, ']',     &  
            '], SIZE(pSend(:))=iyTotSize=[',  iyTotSize, ']'
    
    !
    !-- MPI variable 
    !
    ALLOCATE (ircnt(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Receive Count Allocation Error!!!'
            RETURN
         ENDIF
    ALLOCATE (idisp(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Display Count Allocation Error!!!'
            RETURN
         ENDIF

    !
    !-- Set send count, receive count, display position for MPI
    !
    iscnt   = SIZE(pSend(iy(1):iy(2)))        ! 마스터 프로세스(rank=0)로 전송할 데이터 개수
  
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

            n_mpi_tmp = i
            iscnt2   = iscnt2 + SIZE(pSend(iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))

            !
            !-- set receive count for MPI
            !
            ircnt(i) = SIZE(pSend(iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))
        ENDDO
    ENDIF

    !
    !-- Write MPI Log
    !
    WRITE(unit_num, '(A15, I7, A12, I7, A1, A15, I7, A1, I7, A2)'),     &
            '###  myrank = [', myrank, '], iscnt = [', iscnt, ']',      &
            ' => SIZE(pSend(', iy(1), ':', iy(2), '))'
    
    IF ( myrank == 0 ) THEN
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      ircnt  = ['
        DO i=1, SIZE(ircnt)
            WRITE(unit_num, '(A14, I9)'), ' ', ircnt(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      idisp  = ['
        DO i=1, SIZE(idisp)
            WRITE(unit_num, '(A14, I9)'), ' ', idisp(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
    ENDIF

    WRITE(unit_num, '(A3, I5)') "---iy(1)", iy(1)

    !
    !-- Gather the output From other processes
    !
    CALL MPI_GATHERV(pSend(iy(1))   ,    iscnt, MPI_SHORT,            &
                     pRecv,              ircnt, idisp, MPI_SHORT,     &
                     0, MPI_COMM_WORLD, errCode)

    retCode = errCode

    !
    !--- set pixel size calculated in defined total pixel size
    !
    ix_first = 1            
    ix_last  = 1   ! ix(2)    

    !
    !
    !-- Write MPI Log
    !
    iy_size = iy(2) - iy(1)
    WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(I10))'

    IF (myrank==0) THEN
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(iy(1):iy(2)))
        ENDDO

        DO j=1, gci_nprocs

            n_mpi_tmp = j
            iy_size = iyy(n_mpi_tmp,2) - iyy(n_mpi_tmp,1)
            WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(I10))'
            DO i=ix_first, ix_last
                WRITE(unit_num, w_fmt),                                                 &
                        'myrank = [', n_mpi_tmp-1, '], irecv =',                        &
                        pRecv(iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2))
            ENDDO
        ENDDO
    ELSE
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(iy(1):iy(2)))
        ENDDO
    ENDIF

    DEALLOCATE(ircnt)
    DEALLOCATE(idisp)

    RETURN
END SUBROUTINE GEMS_Share_MOD_MPI_Process_1D_I2


!
!--- 1Dim 32bit 정수형 변수를 병렬처리
!
SUBROUTINE GEMS_Share_MOD_MPI_Process_1D_I4(myrank, pSend, pRecv, ix, iyy, retCode)

! argument parameter
    INTEGER(KIND=4),          INTENT(IN)    :: myrank
    INTEGER(KIND=4), POINTER, INTENT(IN)    :: pSend(:)
    INTEGER(KIND=4), POINTER, INTENT(IN)    :: pRecv(:)
    INTEGER(KIND=4),          INTENT(IN)    :: ix(:)
    INTEGER(KIND=4),          INTENT(IN)    :: iyy(:, :)
    INTEGER(KIND=4),          INTENT(OUT)   :: retCode

! local variables
    INTEGER(KIND=4)             :: errCode = 0
    CHARACTER(LEN=120)          :: w_fmt

    INTEGER(KIND=4)             :: iscnt
    INTEGER(KIND=4)             :: iscnt2
    INTEGER(KIND=4),ALLOCATABLE :: ircnt(:), idisp(:)

    INTEGER(KIND=4)             :: n_mpi
    INTEGER(KIND=4)             :: n_mpi_tmp
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

    WRITE(unit_num, '(/, A15, I7, A, /)'),                 &
            '###  myrank = [',  myrank,                    &
            '] => 1Dim 32bit 정수형 변수를 병렬처리'

    !
    !--- point the 1D pointer variable to 1D target variable
    !
    WRITE(unit_num, '(A15, I7, 2(A17,I7,A1))'),             &
            '###  myrank = [',  myrank,                     &
            '], shape(pSend)=[',   shape(pSend), ']',       &
            '   shape(pRecv)=[',   shape(pRecv), ']'

    ixTotSize = SIZE(pSend(:))
    iyTotSize = SIZE(pSend(:))

    WRITE(unit_num, '(A15, I7, 2(A32,I7,A1))'),                   &
            '###  myrank = [',  myrank,                           &
            '], SIZE(pSend(:))=ixTotSize=[',  ixTotSize, ']',     &  
            '], SIZE(pSend(:))=iyTotSize=[',  iyTotSize, ']'
    
    !
    !-- MPI variable 
    !
    ALLOCATE (ircnt(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Receive Count Allocation Error!!!'
            RETURN
         ENDIF
    ALLOCATE (idisp(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Display Count Allocation Error!!!'
            RETURN
         ENDIF

    !
    !-- Set send count, receive count, display position for MPI
    !
    iscnt   = SIZE(pSend(iy(1):iy(2)))        ! 마스터 프로세스(rank=0)로 전송할 데이터 개수
  
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

            n_mpi_tmp = i
            iscnt2   = iscnt2 + SIZE(pSend(iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))

            !
            !-- set receive count for MPI
            !
            ircnt(i) = SIZE(pSend(iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))
        ENDDO
    ENDIF

    !
    !-- Write MPI Log
    !
    WRITE(unit_num, '(A15, I7, A12, I7, A1, A15, I7, A1, I7, A2)'),     &
            '###  myrank = [', myrank, '], iscnt = [', iscnt, ']',      &
            ' => SIZE(pSend(', iy(1), ':', iy(2), '))'
    
    IF ( myrank == 0 ) THEN
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      ircnt  = ['
        DO i=1, SIZE(ircnt)
            WRITE(unit_num, '(A14, I9)'), ' ', ircnt(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      idisp  = ['
        DO i=1, SIZE(idisp)
            WRITE(unit_num, '(A14, I9)'), ' ', idisp(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
    ENDIF

    WRITE(unit_num, '(A3, I5)') "---iy(1)", iy(1)

    !
    !-- Gather the output From other processes
    !
    CALL MPI_GATHERV(pSend(iy(1))   ,    iscnt, MPI_INTEGER,            &
                     pRecv,              ircnt, idisp, MPI_INTEGER,     &
                     0, MPI_COMM_WORLD, errCode)

    retCode = errCode

    !
    !--- set pixel size calculated in defined total pixel size
    !
    ix_first = 1            
    ix_last  = 1   ! ix(2)    

    !
    !
    !-- Write MPI Log
    !
    iy_size = iy(2) - iy(1)
    WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(I10))'

    IF (myrank==0) THEN
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(iy(1):iy(2)))
        ENDDO

        DO j=1, gci_nprocs

            n_mpi_tmp = j
            iy_size = iyy(n_mpi_tmp,2) - iyy(n_mpi_tmp,1)
            WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(I10))'
            DO i=ix_first, ix_last
                WRITE(unit_num, w_fmt),                                                 &
                        'myrank = [', n_mpi_tmp-1, '], irecv =',                        &
                        pRecv(iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2))
            ENDDO
        ENDDO
    ELSE
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(iy(1):iy(2)))
        ENDDO
    ENDIF

    DEALLOCATE(ircnt)
    DEALLOCATE(idisp)

    RETURN
END SUBROUTINE GEMS_Share_MOD_MPI_Process_1D_I4

!
!--- 1Dim 32bit 실수형 변수를 병렬처리
!
SUBROUTINE GEMS_Share_MOD_MPI_Process_1D_R4(myrank, pSend, pRecv, ix, iyy, retCode)

! argument parameter
    INTEGER(KIND=4),          INTENT(IN)    :: myrank
    REAL(KIND=4),    POINTER, INTENT(IN)    :: pSend(:)
    REAL(KIND=4),    POINTER, INTENT(IN)    :: pRecv(:)
    INTEGER(KIND=4),          INTENT(IN)    :: ix(:)
    INTEGER(KIND=4),          INTENT(IN)    :: iyy(:, :)
    INTEGER(KIND=4),          INTENT(OUT)   :: retCode

! local variables
    INTEGER(KIND=4)             :: errCode = 0
    CHARACTER(LEN=120)          :: w_fmt

    INTEGER(KIND=4)             :: iscnt
    INTEGER(KIND=4)             :: iscnt2
    INTEGER(KIND=4),ALLOCATABLE :: ircnt(:), idisp(:)

    INTEGER(KIND=4)             :: n_mpi
    INTEGER(KIND=4)             :: n_mpi_tmp
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

    WRITE(unit_num, '(/, A15, I7, A, /)'),                 &
            '###  myrank = [',  myrank,                    &
            '] => 1Dim 32bit 실수형 변수를 병렬처리'

    !
    !--- point the 1D pointer variable to 1D target variable
    !
    WRITE(unit_num, '(A15, I7, 2(A17,I7,A1))'),             &
            '###  myrank = [',  myrank,                     &
            '], shape(pSend)=[',   shape(pSend), ']',       &
            '   shape(pRecv)=[',   shape(pRecv), ']'

    ixTotSize = SIZE(pSend(:))
    iyTotSize = SIZE(pSend(:))

    WRITE(unit_num, '(A15, I7, 2(A32,I7,A1))'),                   &
            '###  myrank = [',  myrank,                           &
            '], SIZE(pSend(:))=ixTotSize=[',  ixTotSize, ']',     &  
            '], SIZE(pSend(:))=iyTotSize=[',  iyTotSize, ']'
    
    !
    !-- MPI variable 
    !
    ALLOCATE (ircnt(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Receive Count Allocation Error!!!'
            RETURN
         ENDIF
    ALLOCATE (idisp(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Display Count Allocation Error!!!'
            RETURN
         ENDIF

    !
    !-- Set send count, receive count, display position for MPI
    !
    iscnt   = SIZE(pSend(iy(1):iy(2)))        ! 마스터 프로세스(rank=0)로 전송할 데이터 개수
  
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

            n_mpi_tmp = i
            iscnt2   = iscnt2 + SIZE(pSend(iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))

            !
            !-- set receive count for MPI
            !
            ircnt(i) = SIZE(pSend(iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))
        ENDDO
    ENDIF

    !
    !-- Write MPI Log
    !
    WRITE(unit_num, '(A15, I7, A12, I7, A1, A15, I7, A1, I7, A2)'),     &
            '###  myrank = [', myrank, '], iscnt = [', iscnt, ']',      &
            ' => SIZE(pSend(', iy(1), ':', iy(2), '))'
    
    IF ( myrank == 0 ) THEN
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      ircnt  = ['
        DO i=1, SIZE(ircnt)
            WRITE(unit_num, '(A14, I9)'), ' ', ircnt(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      idisp  = ['
        DO i=1, SIZE(idisp)
            WRITE(unit_num, '(A14, I9)'), ' ', idisp(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
    ENDIF

    WRITE(unit_num, '(A3, I5)') "---iy(1)", iy(1)

    !
    !-- Gather the output From other processes
    !
    CALL MPI_GATHERV(pSend(iy(1))   ,    iscnt, MPI_REAL,            &
                     pRecv,              ircnt, idisp, MPI_REAL,     &
                     0, MPI_COMM_WORLD, errCode)

    retCode = errCode

    !
    !--- set pixel size calculated in defined total pixel size
    !
    ix_first = 1            
    ix_last  = 1   ! ix(2)    

    !
    !
    !-- Write MPI Log
    !
    iy_size = iy(2) - iy(1)
    WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(F15.5))'

    IF (myrank==0) THEN
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(iy(1):iy(2)))
        ENDDO

        DO j=1, gci_nprocs

            n_mpi_tmp = j
            iy_size = iyy(n_mpi_tmp,2) - iyy(n_mpi_tmp,1)
            WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(F15.5))'
            DO i=ix_first, ix_last
                WRITE(unit_num, w_fmt),                                                 &
                        'myrank = [', n_mpi_tmp-1, '], irecv =',                        &
                        pRecv(iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2))
            ENDDO
        ENDDO
    ELSE
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(iy(1):iy(2)))
        ENDDO
    ENDIF

    DEALLOCATE(ircnt)
    DEALLOCATE(idisp)

    RETURN
END SUBROUTINE GEMS_Share_MOD_MPI_Process_1D_R4


!
!--- 1Dim 64bit 실수형 변수를 병렬처리
!
SUBROUTINE GEMS_Share_MOD_MPI_Process_1D_R8(myrank, pSend, pRecv, ix, iyy, retCode)

! argument parameter
    INTEGER(KIND=4),          INTENT(IN)    :: myrank
    REAL(KIND=8),    POINTER, INTENT(IN)    :: pSend(:)
    REAL(KIND=8),    POINTER, INTENT(IN)    :: pRecv(:)
    INTEGER(KIND=4),          INTENT(IN)    :: ix(:)
    INTEGER(KIND=4),          INTENT(IN)    :: iyy(:, :)
    INTEGER(KIND=4),          INTENT(OUT)   :: retCode

! local variables
    INTEGER(KIND=4)             :: errCode = 0
    CHARACTER(LEN=120)          :: w_fmt

    INTEGER(KIND=4)             :: iscnt
    INTEGER(KIND=4)             :: iscnt2
    INTEGER(KIND=4),ALLOCATABLE :: ircnt(:), idisp(:)

    INTEGER(KIND=4)             :: n_mpi
    INTEGER(KIND=4)             :: n_mpi_tmp
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

    WRITE(unit_num, '(/, A15, I7, A, /)'),                 &
            '###  myrank = [',  myrank,                    &
            '] => 1Dim 64bit 실수형 변수를 병렬처리'

    !
    !--- point the 1D pointer variable to 1D target variable
    !
    WRITE(unit_num, '(A15, I7, 2(A17,I7,A1))'),             &
            '###  myrank = [',  myrank,                     &
            '], shape(pSend)=[',   shape(pSend), ']',       &
            '   shape(pRecv)=[',   shape(pRecv), ']'

    ixTotSize = SIZE(pSend(:))
    iyTotSize = SIZE(pSend(:))

    WRITE(unit_num, '(A15, I7, 2(A32,I7,A1))'),                   &
            '###  myrank = [',  myrank,                           &
            '], SIZE(pSend(:))=ixTotSize=[',  ixTotSize, ']',     &  
            '], SIZE(pSend(:))=iyTotSize=[',  iyTotSize, ']'
    
    !
    !-- MPI variable 
    !
    ALLOCATE (ircnt(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Receive Count Allocation Error!!!'
            RETURN
         ENDIF
    ALLOCATE (idisp(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Display Count Allocation Error!!!'
            RETURN
         ENDIF

    !
    !-- Set send count, receive count, display position for MPI
    !
    iscnt   = SIZE(pSend(iy(1):iy(2)))        ! 마스터 프로세스(rank=0)로 전송할 데이터 개수
  
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

            n_mpi_tmp = i
            iscnt2   = iscnt2 + SIZE(pSend(iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))

            !
            !-- set receive count for MPI
            !
            ircnt(i) = SIZE(pSend(iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))
        ENDDO
    ENDIF

    !
    !-- Write MPI Log
    !
    WRITE(unit_num, '(A15, I7, A12, I7, A1, A15, I7, A1, I7, A2)'),     &
            '###  myrank = [', myrank, '], iscnt = [', iscnt, ']',      &
            ' => SIZE(pSend(', iy(1), ':', iy(2), '))'
    
    IF ( myrank == 0 ) THEN
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      ircnt  = ['
        DO i=1, SIZE(ircnt)
            WRITE(unit_num, '(A14, I9)'), ' ', ircnt(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      idisp  = ['
        DO i=1, SIZE(idisp)
            WRITE(unit_num, '(A14, I9)'), ' ', idisp(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
    ENDIF

    WRITE(unit_num, '(A3, I5)') "---iy(1)", iy(1)

    !
    !-- Gather the output From other processes
    !
    CALL MPI_GATHERV(pSend(iy(1))   ,    iscnt, MPI_DOUBLE,            &
                     pRecv,              ircnt, idisp, MPI_DOUBLE,     &
                     0, MPI_COMM_WORLD, errCode)

    retCode = errCode

    !
    !--- set pixel size calculated in defined total pixel size
    !
    ix_first = 1            
    ix_last  = 1   ! ix(2)    

    !
    !
    !-- Write MPI Log
    !
    iy_size = iy(2) - iy(1)
    WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(F15.5))'

    IF (myrank==0) THEN
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(iy(1):iy(2)))
        ENDDO

        DO j=1, gci_nprocs

            n_mpi_tmp = j
            iy_size = iyy(n_mpi_tmp,2) - iyy(n_mpi_tmp,1)
            WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(F15.5))'
            DO i=ix_first, ix_last
                WRITE(unit_num, w_fmt),                                                 &
                        'myrank = [', n_mpi_tmp-1, '], irecv =',                        &
                        pRecv(iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2))
            ENDDO
        ENDDO
    ELSE
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(iy(1):iy(2)))
        ENDDO
    ENDIF

    DEALLOCATE(ircnt)
    DEALLOCATE(idisp)

    RETURN
END SUBROUTINE GEMS_Share_MOD_MPI_Process_1D_R8


!
!--- 2Dim 8bit 정수형 변수를 병렬처리
!
SUBROUTINE GEMS_Share_MOD_MPI_Process_2D_I1(myrank, pSend, pRecv, ix, iyy, retCode)

! argument parameter
    INTEGER(KIND=4),          INTENT(IN)    :: myrank
    INTEGER(KIND=1), POINTER, INTENT(IN)    :: pSend(:, :)
    INTEGER(KIND=1), POINTER, INTENT(IN)    :: pRecv(:, :)
    INTEGER(KIND=4),          INTENT(IN)    :: ix(:)
    INTEGER(KIND=4),          INTENT(IN)    :: iyy(:, :)
    INTEGER(KIND=4),          INTENT(OUT)   :: retCode

! local variables
    INTEGER(KIND=4)             :: errCode = 0
    CHARACTER(LEN=120)          :: w_fmt

    INTEGER(KIND=4)             :: iscnt
    INTEGER(KIND=4)             :: iscnt2
    INTEGER(KIND=4),ALLOCATABLE :: ircnt(:), idisp(:)

    INTEGER(KIND=4)             :: n_mpi
    INTEGER(KIND=4)             :: n_mpi_tmp
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

    WRITE(unit_num, '(/, A15, I7, A, /)'),                 &
            '###  myrank = [',  myrank,                    &
            '] => 2Dim 8bit 정수형 변수를 병렬처리'

    !
    !--- point the 2D pointer variable to 2D target variable
    !
    WRITE(unit_num, '(A15, I7, 2(A17,2I7,A1))'),            &
            '###  myrank = [',  myrank,                     &
            '], shape(pSend)=[',   shape(pSend), ']',       &
            '   shape(pRecv)=[',   shape(pRecv), ']'

    ixTotSize = SIZE(pSend(:,1))
    iyTotSize = SIZE(pSend(1,:))

    WRITE(unit_num, '(A15, I7, 2(A32,I7,A1))'),                   &
            '###  myrank = [',  myrank,                           &
            '], SIZE(pSend(:,1))=ixTotSize=[',  ixTotSize, ']',   &  
            '], SIZE(pSend(1,:))=iyTotSize=[',  iyTotSize, ']'
    
    !
    !-- MPI variable 
    !
    ALLOCATE (ircnt(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Receive Count Allocation Error!!!'
            RETURN
         ENDIF
    ALLOCATE (idisp(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Display Count Allocation Error!!!'
            RETURN
         ENDIF

    !
    !-- Set send count, receive count, display position for MPI
    !
    iscnt   = SIZE(pSend(:, iy(1):iy(2)))        ! 마스터 프로세스(rank=0)로 전송할 데이터 개수
  
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

            n_mpi_tmp = i
            iscnt2   = iscnt2 + SIZE(pSend(:, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))

            !
            !-- set receive count for MPI
            !
            ircnt(i) = SIZE(pSend(:, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))
        ENDDO
    ENDIF

    !
    !-- Write MPI Log
    !
    WRITE(unit_num, '(A15, I7, A12, I7, A1, A18, I7, A1, I7, A2)'),     &
            '###  myrank = [', myrank, '], iscnt = [', iscnt, ']',      &
            ' => SIZE(pSend(:, ', iy(1), ':', iy(2), '))'
    
    IF ( myrank == 0 ) THEN
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      ircnt  = ['
        DO i=1, SIZE(ircnt)
            WRITE(unit_num, '(A14, I9)'), ' ', ircnt(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      idisp  = ['
        DO i=1, SIZE(idisp)
            WRITE(unit_num, '(A14, I9)'), ' ', idisp(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
    ENDIF

    WRITE(unit_num, '(A3, I5)') "---iy(1)", iy(1)

    !
    !-- Gather the output From other processes
    !
    CALL MPI_GATHERV(pSend(:, iy(1)),    iscnt, MPI_CHARACTER,            &
                     pRecv,              ircnt, idisp, MPI_CHARACTER,     &
                     0, MPI_COMM_WORLD, errCode)

    retCode = errCode

    !
    !--- set pixel size calculated in defined total pixel size
    !
    ix_first = 1            
    ix_last  = ix(2)    

    !
    !
    !-- Write MPI Log
    !
    iy_size = iy(2) - iy(1)
    WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(I10))'

    IF (myrank==0) THEN
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(i, iy(1):iy(2)))
        ENDDO

        DO j=1, gci_nprocs

            n_mpi_tmp = j
            iy_size = iyy(n_mpi_tmp,2) - iyy(n_mpi_tmp,1)
            WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(I10))'
            DO i=ix_first, ix_last
                WRITE(unit_num, w_fmt),                                                 &
                        'myrank = [', n_mpi_tmp-1, '], irecv =',                        &
                        pRecv(i, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2))
            ENDDO
        ENDDO
    ELSE
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(i, iy(1):iy(2)))
        ENDDO
    ENDIF

    DEALLOCATE(ircnt)
    DEALLOCATE(idisp)

    RETURN
END SUBROUTINE GEMS_Share_MOD_MPI_Process_2D_I1


!
!--- 2Dim 16bit 정수형 변수를 병렬처리
!
SUBROUTINE GEMS_Share_MOD_MPI_Process_2D_I2(myrank, pSend, pRecv, ix, iyy, retCode)

! argument parameter
    INTEGER(KIND=4),          INTENT(IN)    :: myrank
    INTEGER(KIND=2), POINTER, INTENT(IN)    :: pSend(:, :)
    INTEGER(KIND=2), POINTER, INTENT(IN)    :: pRecv(:, :)
    INTEGER(KIND=4),          INTENT(IN)    :: ix(:)
    INTEGER(KIND=4),          INTENT(IN)    :: iyy(:, :)
    INTEGER(KIND=4),          INTENT(OUT)   :: retCode

! local variables
    INTEGER(KIND=4)             :: errCode = 0
    CHARACTER(LEN=120)          :: w_fmt

    INTEGER(KIND=4)             :: iscnt
    INTEGER(KIND=4)             :: iscnt2
    INTEGER(KIND=4),ALLOCATABLE :: ircnt(:), idisp(:)

    INTEGER(KIND=4)             :: n_mpi
    INTEGER(KIND=4)             :: n_mpi_tmp
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

    WRITE(unit_num, '(/, A15, I7, A, /)'),                 &
            '###  myrank = [',  myrank,                    &
            '] => 2Dim 16bit 정수형 변수를 병렬처리'

    !
    !--- point the 2D pointer variable to 2D target variable
    !
    WRITE(unit_num, '(A15, I7, 2(A17,2I7,A1))'),            &
            '###  myrank = [',  myrank,                     &
            '], shape(pSend)=[',   shape(pSend), ']',       &
            '   shape(pRecv)=[',   shape(pRecv), ']'

    ixTotSize = SIZE(pSend(:,1))
    iyTotSize = SIZE(pSend(1,:))

    WRITE(unit_num, '(A15, I7, 2(A32,I7,A1))'),                   &
            '###  myrank = [',  myrank,                           &
            '], SIZE(pSend(:,1))=ixTotSize=[',  ixTotSize, ']',   &  
            '], SIZE(pSend(1,:))=iyTotSize=[',  iyTotSize, ']'
    
    !
    !-- MPI variable 
    !
    ALLOCATE (ircnt(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Receive Count Allocation Error!!!'
            RETURN
         ENDIF
    ALLOCATE (idisp(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Display Count Allocation Error!!!'
            RETURN
         ENDIF

    !
    !-- Set send count, receive count, display position for MPI
    !
    iscnt   = SIZE(pSend(:, iy(1):iy(2)))        ! 마스터 프로세스(rank=0)로 전송할 데이터 개수
  
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

            n_mpi_tmp = i
            iscnt2   = iscnt2 + SIZE(pSend(:, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))

            !
            !-- set receive count for MPI
            !
            ircnt(i) = SIZE(pSend(:, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))
        ENDDO
    ENDIF

    !
    !-- Write MPI Log
    !
    WRITE(unit_num, '(A15, I7, A12, I7, A1, A18, I7, A1, I7, A2)'),     &
            '###  myrank = [', myrank, '], iscnt = [', iscnt, ']',      &
            ' => SIZE(pSend(:, ', iy(1), ':', iy(2), '))'
    
    IF ( myrank == 0 ) THEN
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      ircnt  = ['
        DO i=1, SIZE(ircnt)
            WRITE(unit_num, '(A14, I9)'), ' ', ircnt(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      idisp  = ['
        DO i=1, SIZE(idisp)
            WRITE(unit_num, '(A14, I9)'), ' ', idisp(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
    ENDIF

    WRITE(unit_num, '(A3, I5)') "---iy(1)", iy(1)

    !
    !-- Gather the output From other processes
    !
    CALL MPI_GATHERV(pSend(:, iy(1)),    iscnt, MPI_SHORT,            &
                     pRecv,              ircnt, idisp, MPI_SHORT,     &
                     0, MPI_COMM_WORLD, errCode)

    retCode = errCode

    !
    !--- set pixel size calculated in defined total pixel size
    !
    ix_first = 1            
    ix_last  = ix(2)    

    !
    !
    !-- Write MPI Log
    !
    iy_size = iy(2) - iy(1)
    WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(I10))'

    IF (myrank==0) THEN
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(i, iy(1):iy(2)))
        ENDDO

        DO j=1, gci_nprocs

            n_mpi_tmp = j
            iy_size = iyy(n_mpi_tmp,2) - iyy(n_mpi_tmp,1)
            WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(I10))'
            DO i=ix_first, ix_last
                WRITE(unit_num, w_fmt),                                                 &
                        'myrank = [', n_mpi_tmp-1, '], irecv =',                        &
                        pRecv(i, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2))
            ENDDO
        ENDDO
    ELSE
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(i, iy(1):iy(2)))
        ENDDO
    ENDIF

    DEALLOCATE(ircnt)
    DEALLOCATE(idisp)

    RETURN
END SUBROUTINE GEMS_Share_MOD_MPI_Process_2D_I2


!
!--- 2Dim 32bit 정수형 변수를 병렬처리
!
SUBROUTINE GEMS_Share_MOD_MPI_Process_2D_I4(myrank, pSend, pRecv, ix, iyy, retCode)

! argument parameter
    INTEGER(KIND=4),          INTENT(IN)    :: myrank
    INTEGER(KIND=4), POINTER, INTENT(IN)    :: pSend(:, :)
    INTEGER(KIND=4), POINTER, INTENT(IN)    :: pRecv(:, :)
    INTEGER(KIND=4),          INTENT(IN)    :: ix(:)
    INTEGER(KIND=4),          INTENT(IN)    :: iyy(:, :)
    INTEGER(KIND=4),          INTENT(OUT)   :: retCode

! local variables
    INTEGER(KIND=4)             :: errCode = 0
    CHARACTER(LEN=120)          :: w_fmt

    INTEGER(KIND=4)             :: iscnt
    INTEGER(KIND=4)             :: iscnt2
    INTEGER(KIND=4),ALLOCATABLE :: ircnt(:), idisp(:)

    INTEGER(KIND=4)             :: n_mpi
    INTEGER(KIND=4)             :: n_mpi_tmp
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

    WRITE(unit_num, '(/, A15, I7, A, /)'),                 &
            '###  myrank = [',  myrank,                    &
            '] => 2Dim 32bit 정수형 변수를 병렬처리'

    !
    !--- point the 2D pointer variable to 2D target variable
    !
    WRITE(unit_num, '(A15, I7, 2(A17,2I7,A1))'),            &
            '###  myrank = [',  myrank,                     &
            '], shape(pSend)=[',   shape(pSend), ']',       &
            '   shape(pRecv)=[',   shape(pRecv), ']'

    ixTotSize = SIZE(pSend(:,1))
    iyTotSize = SIZE(pSend(1,:))

    WRITE(unit_num, '(A15, I7, 2(A32,I7,A1))'),                   &
            '###  myrank = [',  myrank,                           &
            '], SIZE(pSend(:,1))=ixTotSize=[',  ixTotSize, ']',   &  
            '], SIZE(pSend(1,:))=iyTotSize=[',  iyTotSize, ']'
    
    !
    !-- MPI variable 
    !
    ALLOCATE (ircnt(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Receive Count Allocation Error!!!'
            RETURN
         ENDIF
    ALLOCATE (idisp(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Display Count Allocation Error!!!'
            RETURN
         ENDIF

    !
    !-- Set send count, receive count, display position for MPI
    !
    iscnt   = SIZE(pSend(:, iy(1):iy(2)))        ! 마스터 프로세스(rank=0)로 전송할 데이터 개수
  
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

            n_mpi_tmp = i
            iscnt2   = iscnt2 + SIZE(pSend(:, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))

            !
            !-- set receive count for MPI
            !
            ircnt(i) = SIZE(pSend(:, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))
        ENDDO
    ENDIF

    !
    !-- Write MPI Log
    !
    WRITE(unit_num, '(A15, I7, A12, I7, A1, A18, I7, A1, I7, A2)'),     &
            '###  myrank = [', myrank, '], iscnt = [', iscnt, ']',      &
            ' => SIZE(pSend(:, ', iy(1), ':', iy(2), '))'
    
    IF ( myrank == 0 ) THEN
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      ircnt  = ['
        DO i=1, SIZE(ircnt)
            WRITE(unit_num, '(A14, I9)'), ' ', ircnt(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      idisp  = ['
        DO i=1, SIZE(idisp)
            WRITE(unit_num, '(A14, I9)'), ' ', idisp(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
    ENDIF

    WRITE(unit_num, '(A3, I5)') "---iy(1)", iy(1)

    !
    !-- Gather the output From other processes
    !
    CALL MPI_GATHERV(pSend(:, iy(1)),    iscnt, MPI_INTEGER,            &
                     pRecv,              ircnt, idisp, MPI_INTEGER,     &
                     0, MPI_COMM_WORLD, errCode)

    retCode = errCode

    !
    !--- set pixel size calculated in defined total pixel size
    !
    ix_first = 1            
    ix_last  = ix(2)    

    !
    !
    !-- Write MPI Log
    !
    iy_size = iy(2) - iy(1)
    WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(I10))'

    IF (myrank==0) THEN
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(i, iy(1):iy(2)))
        ENDDO

        DO j=1, gci_nprocs

            n_mpi_tmp = j
            iy_size = iyy(n_mpi_tmp,2) - iyy(n_mpi_tmp,1)
            WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(I10))'
            DO i=ix_first, ix_last
                WRITE(unit_num, w_fmt),                                                 &
                        'myrank = [', n_mpi_tmp-1, '], irecv =',                        &
                        pRecv(i, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2))
            ENDDO
        ENDDO
    ELSE
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(i, iy(1):iy(2)))
        ENDDO
    ENDIF

    DEALLOCATE(ircnt)
    DEALLOCATE(idisp)

    RETURN
END SUBROUTINE GEMS_Share_MOD_MPI_Process_2D_I4

!
!--- 2Dim 32bit 실수형 변수를 병렬처리
!
SUBROUTINE GEMS_Share_MOD_MPI_Process_2D_R4(myrank, pSend, pRecv, ix, iyy, retCode)

! argument parameter
    INTEGER(KIND=4),          INTENT(IN)    :: myrank
    REAL(KIND=4),    POINTER, INTENT(IN)    :: pSend(:, :)
    REAL(KIND=4),    POINTER, INTENT(IN)    :: pRecv(:, :)
    INTEGER(KIND=4),          INTENT(IN)    :: ix(:)
    INTEGER(KIND=4),          INTENT(IN)    :: iyy(:, :)
    INTEGER(KIND=4),          INTENT(OUT)   :: retCode

! local variables
    INTEGER(KIND=4)             :: errCode = 0
    CHARACTER(LEN=120)          :: w_fmt

    INTEGER(KIND=4)             :: iscnt
    INTEGER(KIND=4)             :: iscnt2
    INTEGER(KIND=4),ALLOCATABLE :: ircnt(:), idisp(:)

    INTEGER(KIND=4)             :: n_mpi
    INTEGER(KIND=4)             :: n_mpi_tmp
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

    WRITE(unit_num, '(/, A15, I7, A, /)'),                 &
            '###  myrank = [',  myrank,                    &
            '] => 2Dim 32bit 실수형 변수를 병렬처리'

    !
    !--- point the 2D pointer variable to 2D target variable
    !
    WRITE(unit_num, '(A15, I7, 2(A17,2I7,A1))'),            &
            '###  myrank = [',  myrank,                     &
            '], shape(pSend)=[',   shape(pSend), ']',       &
            '   shape(pRecv)=[',   shape(pRecv), ']'

    ixTotSize = SIZE(pSend(:,1))
    iyTotSize = SIZE(pSend(1,:))

    WRITE(unit_num, '(A15, I7, 2(A32,I7,A1))'),                   &
            '###  myrank = [',  myrank,                           &
            '], SIZE(pSend(:,1))=ixTotSize=[',  ixTotSize, ']',   &  
            '], SIZE(pSend(1,:))=iyTotSize=[',  iyTotSize, ']'
    
    !
    !-- MPI variable 
    !
    ALLOCATE (ircnt(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Receive Count Allocation Error!!!'
            RETURN
         ENDIF
    ALLOCATE (idisp(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Display Count Allocation Error!!!'
            RETURN
         ENDIF

    !
    !-- Set send count, receive count, display position for MPI
    !
    iscnt   = SIZE(pSend(:, iy(1):iy(2)))        ! 마스터 프로세스(rank=0)로 전송할 데이터 개수
  
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

            n_mpi_tmp = i
            iscnt2   = iscnt2 + SIZE(pSend(:, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))

            !
            !-- set receive count for MPI
            !
            ircnt(i) = SIZE(pSend(:, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))
        ENDDO
    ENDIF

    !
    !-- Write MPI Log
    !
    WRITE(unit_num, '(A15, I7, A12, I7, A1, A18, I7, A1, I7, A2)'),     &
            '###  myrank = [', myrank, '], iscnt = [', iscnt, ']',      &
            ' => SIZE(pSend(:, ', iy(1), ':', iy(2), '))'
    
    IF ( myrank == 0 ) THEN
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      ircnt  = ['
        DO i=1, SIZE(ircnt)
            WRITE(unit_num, '(A14, I9)'), ' ', ircnt(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      idisp  = ['
        DO i=1, SIZE(idisp)
            WRITE(unit_num, '(A14, I9)'), ' ', idisp(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
    ENDIF

    WRITE(unit_num, '(A3, I5)') "---iy(1)", iy(1)

    !
    !-- Gather the output From other processes
    !
    CALL MPI_GATHERV(pSend(:, iy(1)),    iscnt, MPI_REAL,            &
                     pRecv,              ircnt, idisp, MPI_REAL,     &
                     0, MPI_COMM_WORLD, errCode)

    retCode = errCode

    !
    !--- set pixel size calculated in defined total pixel size
    !
    ix_first = 1            
    ix_last  = ix(2)    

    !
    !
    !-- Write MPI Log
    !
    iy_size = iy(2) - iy(1)
    WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(F15.5))'

    IF (myrank==0) THEN
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(i, iy(1):iy(2)))
        ENDDO

        DO j=1, gci_nprocs

            n_mpi_tmp = j
            iy_size = iyy(n_mpi_tmp,2) - iyy(n_mpi_tmp,1)
            WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(F15.5))'
            DO i=ix_first, ix_last
                WRITE(unit_num, w_fmt),                                                 &
                        'myrank = [', n_mpi_tmp-1, '], irecv =',                        &
                        pRecv(i, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2))
            ENDDO
        ENDDO
    ELSE
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(i, iy(1):iy(2)))
        ENDDO
    ENDIF

    DEALLOCATE(ircnt)
    DEALLOCATE(idisp)

    RETURN
END SUBROUTINE GEMS_Share_MOD_MPI_Process_2D_R4


!
!--- 2Dim 64bit 실수형 변수를 병렬처리
!
SUBROUTINE GEMS_Share_MOD_MPI_Process_2D_R8(myrank, pSend, pRecv, ix, iyy, retCode)

! argument parameter
    INTEGER(KIND=4),          INTENT(IN)    :: myrank
    REAL(KIND=8),    POINTER, INTENT(IN)    :: pSend(:, :)
    REAL(KIND=8),    POINTER, INTENT(IN)    :: pRecv(:, :)
    INTEGER(KIND=4),          INTENT(IN)    :: ix(:)
    INTEGER(KIND=4),          INTENT(IN)    :: iyy(:, :)
    INTEGER(KIND=4),          INTENT(OUT)   :: retCode

! local variables
    INTEGER(KIND=4)             :: errCode = 0
    CHARACTER(LEN=120)          :: w_fmt

    INTEGER(KIND=4)             :: iscnt
    INTEGER(KIND=4)             :: iscnt2
    INTEGER(KIND=4),ALLOCATABLE :: ircnt(:), idisp(:)

    INTEGER(KIND=4)             :: n_mpi
    INTEGER(KIND=4)             :: n_mpi_tmp
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

    WRITE(unit_num, '(/, A15, I7, A, /)'),                 &
            '###  myrank = [',  myrank,                    &
            '] => 2Dim 64bit 실수형 변수를 병렬처리'

    !
    !--- point the 2D pointer variable to 2D target variable
    !
    WRITE(unit_num, '(A15, I7, 2(A17,2I7,A1))'),            &
            '###  myrank = [',  myrank,                     &
            '], shape(pSend)=[',   shape(pSend), ']',       &
            '   shape(pRecv)=[',   shape(pRecv), ']'

    ixTotSize = SIZE(pSend(:,1))
    iyTotSize = SIZE(pSend(1,:))

    WRITE(unit_num, '(A15, I7, 2(A32,I7,A1))'),                   &
            '###  myrank = [',  myrank,                           &
            '], SIZE(pSend(:,1))=ixTotSize=[',  ixTotSize, ']',   &  
            '], SIZE(pSend(1,:))=iyTotSize=[',  iyTotSize, ']'
    
    !
    !-- MPI variable 
    !
    ALLOCATE (ircnt(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Receive Count Allocation Error!!!'
            RETURN
         ENDIF
    ALLOCATE (idisp(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Display Count Allocation Error!!!'
            RETURN
         ENDIF

    !
    !-- Set send count, receive count, display position for MPI
    !
    iscnt   = SIZE(pSend(:, iy(1):iy(2)))        ! 마스터 프로세스(rank=0)로 전송할 데이터 개수
  
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

            n_mpi_tmp = i
            iscnt2   = iscnt2 + SIZE(pSend(:, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))

            !
            !-- set receive count for MPI
            !
            ircnt(i) = SIZE(pSend(:, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))
        ENDDO
    ENDIF

    !
    !-- Write MPI Log
    !
    WRITE(unit_num, '(A15, I7, A12, I7, A1, A18, I7, A1, I7, A2)'),     &
            '###  myrank = [', myrank, '], iscnt = [', iscnt, ']',      &
            ' => SIZE(pSend(:, ', iy(1), ':', iy(2), '))'
    
    IF ( myrank == 0 ) THEN
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      ircnt  = ['
        DO i=1, SIZE(ircnt)
            WRITE(unit_num, '(A14, I9)'), ' ', ircnt(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      idisp  = ['
        DO i=1, SIZE(idisp)
            WRITE(unit_num, '(A14, I9)'), ' ', idisp(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
    ENDIF

    WRITE(unit_num, '(A3, I5)') "---iy(1)", iy(1)

    !
    !-- Gather the output From other processes
    !
    CALL MPI_GATHERV(pSend(:, iy(1)),    iscnt, MPI_DOUBLE,            &
                     pRecv,              ircnt, idisp, MPI_DOUBLE,     &
                     0, MPI_COMM_WORLD, errCode)

    retCode = errCode

    !
    !--- set pixel size calculated in defined total pixel size
    !
    ix_first = 1            
    ix_last  = ix(2)    

    !
    !
    !-- Write MPI Log
    !
    iy_size = iy(2) - iy(1)
    WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(F15.5))'

    IF (myrank==0) THEN
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(i, iy(1):iy(2)))
        ENDDO

        DO j=1, gci_nprocs

            n_mpi_tmp = j
            iy_size = iyy(n_mpi_tmp,2) - iyy(n_mpi_tmp,1)
            WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(F15.5))'
            DO i=ix_first, ix_last
                WRITE(unit_num, w_fmt),                                                 &
                        'myrank = [', n_mpi_tmp-1, '], irecv =',                        &
                        pRecv(i, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2))
            ENDDO
        ENDDO
    ELSE
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(i, iy(1):iy(2)))
        ENDDO
    ENDIF

    DEALLOCATE(ircnt)
    DEALLOCATE(idisp)

    RETURN
END SUBROUTINE GEMS_Share_MOD_MPI_Process_2D_R8


!
!--- 3Dim 8bit 정수형 변수를 병렬처리
!
SUBROUTINE GEMS_Share_MOD_MPI_Process_3D_I1(myrank, pSend, pRecv, ix, iyy, retCode)

! argument parameter
    INTEGER(KIND=4),          INTENT(IN)    :: myrank
    INTEGER(KIND=1), POINTER, INTENT(IN)    :: pSend(:, :, :)
    INTEGER(KIND=1), POINTER, INTENT(IN)    :: pRecv(:, :, :)
    INTEGER(KIND=4),          INTENT(IN)    :: ix(:)
    INTEGER(KIND=4),          INTENT(IN)    :: iyy(:, :)
    INTEGER(KIND=4),          INTENT(OUT)   :: retCode

! local variables
    INTEGER(KIND=4)             :: errCode = 0
    CHARACTER(LEN=120)          :: w_fmt

    INTEGER(KIND=4)             :: iscnt
    INTEGER(KIND=4)             :: iscnt2
    INTEGER(KIND=4),ALLOCATABLE :: ircnt(:), idisp(:)

    INTEGER(KIND=4)             :: n_mpi
    INTEGER(KIND=4)             :: n_mpi_tmp
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

    WRITE(unit_num, '(/, A15, I7, A, /)'),                 &
            '###  myrank = [',  myrank,                    &
            '] => 3Dim 8bit 정수형 변수를 병렬처리'

    !
    !--- point the 3D pointer variable to 3D target variable
    !
    WRITE(unit_num, '(A15, I7, 2(A17,3I7,A1))'),            &
            '###  myrank = [',  myrank,                     &
            '], shape(pSend)=[',   shape(pSend), ']',       &
            '   shape(pRecv)=[',   shape(pRecv), ']'

    ixTotSize = SIZE(pSend(1,:,1))
    iyTotSize = SIZE(pSend(1,1,:))

    WRITE(unit_num, '(A15, I7, 2(A32,I7,A1))'),                    &
            '###  myrank = [',  myrank,                            &
            '], SIZE(pSend(:, :,1))=ixTotSize=[',  ixTotSize, ']', &  
            '], SIZE(pSend(:,1,:))=iyTotSize=[',  iyTotSize, ']'
    
    !
    !-- MPI variable 
    !
    ALLOCATE (ircnt(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Receive Count Allocation Error!!!'
            RETURN
         ENDIF
    ALLOCATE (idisp(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Display Count Allocation Error!!!'
            RETURN
         ENDIF

    !
    !-- Set send count, receive count, display position for MPI
    !
    iscnt   = SIZE(pSend(:, :, iy(1):iy(2)))        ! 마스터 프로세스(rank=0)로 전송할 데이터 개수
  
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

            n_mpi_tmp = i
            iscnt2   = iscnt2 + SIZE(pSend(:, :, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))

            !
            !-- set receive count for MPI
            !
            ircnt(i) = SIZE(pSend(:, :, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))
        ENDDO
    ENDIF

    !
    !-- Write MPI Log
    !
    WRITE(unit_num, '(A15, I7, A12, I7, A1, A21, I7, A1, I7, A2)'),     &
            '###  myrank = [', myrank, '], iscnt = [', iscnt, ']',      &
            ' => SIZE(pSend(:, :, ', iy(1), ':', iy(2), '))'
    
    IF ( myrank == 0 ) THEN
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      ircnt  = ['
        DO i=1, SIZE(ircnt)
            WRITE(unit_num, '(A14, I9)'), ' ', ircnt(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      idisp  = ['
        DO i=1, SIZE(idisp)
            WRITE(unit_num, '(A14, I9)'), ' ', idisp(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
    ENDIF

    WRITE(unit_num, '(A3, I5)') "---iy(1)", iy(1)

    !
    !-- Gather the output From other processes
    !
    CALL MPI_GATHERV(pSend(:, :, iy(1)), iscnt, MPI_CHARACTER,            &
                     pRecv,              ircnt, idisp, MPI_CHARACTER,     &
                     0, MPI_COMM_WORLD, errCode)

    retCode = errCode

    !
    !--- set pixel size calculated in defined total pixel size
    !
    ix_first = 1            
    ix_last  = ix(2)    

    !
    !
    !-- Write MPI Log
    !
    iy_size = iy(2) - iy(1)
    WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(I10))'

    IF (myrank==0) THEN
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(:, i, iy(1):iy(2)))
        ENDDO

        DO j=1, gci_nprocs

            n_mpi_tmp = j
            iy_size = iyy(n_mpi_tmp,2) - iyy(n_mpi_tmp,1)
            WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(I10))'
            DO i=ix_first, ix_last
                WRITE(unit_num, w_fmt),                                                 &
                        'myrank = [', n_mpi_tmp-1, '], irecv =',                        &
                        pRecv(:, i, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2))
            ENDDO
        ENDDO
    ELSE
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(:, i, iy(1):iy(2)))
        ENDDO
    ENDIF

    DEALLOCATE(ircnt)
    DEALLOCATE(idisp)

    RETURN
END SUBROUTINE GEMS_Share_MOD_MPI_Process_3D_I1


!
!--- 3Dim 16bit 정수형 변수를 병렬처리
!
SUBROUTINE GEMS_Share_MOD_MPI_Process_3D_I2(myrank, pSend, pRecv, ix, iyy, retCode)

! argument parameter
    INTEGER(KIND=4),          INTENT(IN)    :: myrank
    INTEGER(KIND=2), POINTER, INTENT(IN)    :: pSend(:, :, :)
    INTEGER(KIND=2), POINTER, INTENT(IN)    :: pRecv(:, :, :)
    INTEGER(KIND=4),          INTENT(IN)    :: ix(:)
    INTEGER(KIND=4),          INTENT(IN)    :: iyy(:, :)
    INTEGER(KIND=4),          INTENT(OUT)   :: retCode

! local variables
    INTEGER(KIND=4)             :: errCode = 0
    CHARACTER(LEN=120)          :: w_fmt

    INTEGER(KIND=4)             :: iscnt
    INTEGER(KIND=4)             :: iscnt2
    INTEGER(KIND=4),ALLOCATABLE :: ircnt(:), idisp(:)

    INTEGER(KIND=4)             :: n_mpi
    INTEGER(KIND=4)             :: n_mpi_tmp
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

    WRITE(unit_num, '(/, A15, I7, A, /)'),                 &
            '###  myrank = [',  myrank,                    &
            '] => 3Dim 16bit 정수형 변수를 병렬처리'

    !
    !--- point the 3D pointer variable to 3D target variable
    !
    WRITE(unit_num, '(A15, I7, 2(A17,3I7,A1))'),            &
            '###  myrank = [',  myrank,                     &
            '], shape(pSend)=[',   shape(pSend), ']',       &
            '   shape(pRecv)=[',   shape(pRecv), ']'

    ixTotSize = SIZE(pSend(1,:,1))
    iyTotSize = SIZE(pSend(1,1,:))

    WRITE(unit_num, '(A15, I7, 2(A32,I7,A1))'),                      &
            '###  myrank = [',  myrank,                              &
            '], SIZE(pSend(:, :,1))=ixTotSize=[',  ixTotSize, ']',   &  
            '], SIZE(pSend(:,1,:))=iyTotSize=[',  iyTotSize, ']'
    
    !
    !-- MPI variable 
    !
    ALLOCATE (ircnt(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Receive Count Allocation Error!!!'
            RETURN
         ENDIF
    ALLOCATE (idisp(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Display Count Allocation Error!!!'
            RETURN
         ENDIF

    !
    !-- Set send count, receive count, display position for MPI
    !
    iscnt   = SIZE(pSend(:, :, iy(1):iy(2)))        ! 마스터 프로세스(rank=0)로 전송할 데이터 개수
  
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

            n_mpi_tmp = i
            iscnt2   = iscnt2 + SIZE(pSend(:, :, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))

            !
            !-- set receive count for MPI
            !
            ircnt(i) = SIZE(pSend(:, :, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))
        ENDDO
    ENDIF

    !
    !-- Write MPI Log
    !
    WRITE(unit_num, '(A15, I7, A12, I7, A1, A21, I7, A1, I7, A2)'),     &
            '###  myrank = [', myrank, '], iscnt = [', iscnt, ']',      &
            ' => SIZE(pSend(:, :, ', iy(1), ':', iy(2), '))'

    IF ( myrank == 0 ) THEN
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      ircnt  = ['
        DO i=1, SIZE(ircnt)
            WRITE(unit_num, '(A14, I9)'), ' ', ircnt(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      idisp  = ['
        DO i=1, SIZE(idisp)
            WRITE(unit_num, '(A14, I9)'), ' ', idisp(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
    ENDIF

    WRITE(unit_num, '(A3, I5)') "---iy(1)", iy(1)

    !
    !-- Gather the output From other processes
    !
    CALL MPI_GATHERV(pSend(:, :, iy(1)), iscnt, MPI_SHORT,            &
                     pRecv,              ircnt, idisp, MPI_SHORT,     &
                     0, MPI_COMM_WORLD, errCode)

    retCode = errCode

    !
    !--- set pixel size calculated in defined total pixel size
    !
    ix_first = 1            
    ix_last  = ix(2)    

    !
    !
    !-- Write MPI Log
    !
    iy_size = iy(2) - iy(1)
    WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(I10))'

    IF (myrank==0) THEN
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(:, i, iy(1):iy(2)))
        ENDDO

        DO j=1, gci_nprocs

            n_mpi_tmp = j
            iy_size = iyy(n_mpi_tmp,2) - iyy(n_mpi_tmp,1)
            WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(I10))'
            DO i=ix_first, ix_last
                WRITE(unit_num, w_fmt),                                                 &
                        'myrank = [', n_mpi_tmp-1, '], irecv =',                        &
                        pRecv(:, i, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2))
            ENDDO
        ENDDO
    ELSE
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(:, i, iy(1):iy(2)))
        ENDDO
    ENDIF

    DEALLOCATE(ircnt)
    DEALLOCATE(idisp)

    RETURN
END SUBROUTINE GEMS_Share_MOD_MPI_Process_3D_I2

!
!--- 3Dim 32bit 정수형 변수를 병렬처리
!
SUBROUTINE GEMS_Share_MOD_MPI_Process_3D_I4(myrank, pSend, pRecv, ix, iyy, retCode)

! argument parameter
    INTEGER(KIND=4),          INTENT(IN)    :: myrank
    INTEGER(KIND=4), POINTER, INTENT(IN)    :: pSend(:, :, :)
    INTEGER(KIND=4), POINTER, INTENT(IN)    :: pRecv(:, :, :)
    INTEGER(KIND=4),          INTENT(IN)    :: ix(:)
    INTEGER(KIND=4),          INTENT(IN)    :: iyy(:, :)
    INTEGER(KIND=4),          INTENT(OUT)   :: retCode

! local variables
    INTEGER(KIND=4)             :: errCode = 0
    CHARACTER(LEN=120)          :: w_fmt

    INTEGER(KIND=4)             :: iscnt
    INTEGER(KIND=4)             :: iscnt2
    INTEGER(KIND=4),ALLOCATABLE :: ircnt(:), idisp(:)

    INTEGER(KIND=4)             :: n_mpi
    INTEGER(KIND=4)             :: n_mpi_tmp
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

    WRITE(unit_num, '(/, A15, I7, A, /)'),                 &
            '###  myrank = [',  myrank,                    &
            '] => 3Dim 32bit 정수형 변수를 병렬처리'

    !
    !--- point the 3D pointer variable to 3D target variable
    !
    WRITE(unit_num, '(A15, I7, 2(A17,3I7,A1))'),            &
            '###  myrank = [',  myrank,                     &
            '], shape(pSend)=[',   shape(pSend), ']',       &
            '   shape(pRecv)=[',   shape(pRecv), ']'

    ixTotSize = SIZE(pSend(1,:,1))
    iyTotSize = SIZE(pSend(1,1,:))

    WRITE(unit_num, '(A15, I7, 2(A32,I7,A1))'),                     &
            '###  myrank = [',  myrank,                             &
            '], SIZE(pSend(1,:,1))=ixTotSize=[',  ixTotSize, ']',   &  
            '], SIZE(pSend(1,1,:))=iyTotSize=[',  iyTotSize, ']'
    
    !
    !-- MPI variable 
    !
    ALLOCATE (ircnt(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Receive Count Allocation Error!!!'
            RETURN
         ENDIF
    ALLOCATE (idisp(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Display Count Allocation Error!!!'
            RETURN
         ENDIF

    !
    !-- Set send count, receive count, display position for MPI
    !
    iscnt   = SIZE(pSend(:, :, iy(1):iy(2)))        ! 마스터 프로세스(rank=0)로 전송할 데이터 개수
  
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

            n_mpi_tmp = i
            iscnt2   = iscnt2 + SIZE(pSend(:, :, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))

            !
            !-- set receive count for MPI
            !
            ircnt(i) = SIZE(pSend(:, :, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))
        ENDDO
    ENDIF

    !
    !-- Write MPI Log
    !
    WRITE(unit_num, '(A15, I7, A12, I7, A1, A21, I7, A1, I7, A2)'),     &
            '###  myrank = [', myrank, '], iscnt = [', iscnt, ']',      &
            ' => SIZE(pSend(:, :, ', iy(1), ':', iy(2), '))'
    
    IF ( myrank == 0 ) THEN
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      ircnt  = ['
        DO i=1, SIZE(ircnt)
            WRITE(unit_num, '(A14, I9)'), ' ', ircnt(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      idisp  = ['
        DO i=1, SIZE(idisp)
            WRITE(unit_num, '(A14, I9)'), ' ', idisp(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
    ENDIF

    WRITE(unit_num, '(A3, I5)') "---iy(1)", iy(1)

    !
    !-- Gather the output From other processes
    !
    CALL MPI_GATHERV(pSend(:, :, iy(1)), iscnt, MPI_INTEGER,            &
                     pRecv,              ircnt, idisp, MPI_INTEGER,     &
                     0, MPI_COMM_WORLD, errCode)

    retCode = errCode

    !
    !--- set pixel size calculated in defined total pixel size
    !
    ix_first = 1            
    ix_last  = ix(2)    

    !
    !
    !-- Write MPI Log
    !
    iy_size = iy(2) - iy(1)
    WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(I10))'

    IF (myrank==0) THEN
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(:, i, iy(1):iy(2)))
        ENDDO

        DO j=1, gci_nprocs

            n_mpi_tmp = j
            iy_size = iyy(n_mpi_tmp,2) - iyy(n_mpi_tmp,1)
            WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(I10))'
            DO i=ix_first, ix_last
                WRITE(unit_num, w_fmt),                                                 &
                        'myrank = [', n_mpi_tmp-1, '], irecv =',                        &
                        pRecv(:, i, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2))
            ENDDO
        ENDDO
    ELSE
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(:, i, iy(1):iy(2)))
        ENDDO
    ENDIF

    DEALLOCATE(ircnt)
    DEALLOCATE(idisp)

    RETURN
END SUBROUTINE GEMS_Share_MOD_MPI_Process_3D_I4

!
!--- 3Dim 32bit 실수형 변수를 병렬처리
!
SUBROUTINE GEMS_Share_MOD_MPI_Process_3D_R4(myrank, pSend, pRecv, ix, iyy, retCode)

! argument parameter
    INTEGER(KIND=4),          INTENT(IN)    :: myrank
    REAL(KIND=4),    POINTER, INTENT(IN)    :: pSend(:, :, :)
    REAL(KIND=4),    POINTER, INTENT(IN)    :: pRecv(:, :, :)
    INTEGER(KIND=4),          INTENT(IN)    :: ix(:)
    INTEGER(KIND=4),          INTENT(IN)    :: iyy(:, :)
    INTEGER(KIND=4),          INTENT(OUT)   :: retCode

! local variables
    INTEGER(KIND=4)             :: errCode = 0
    CHARACTER(LEN=120)          :: w_fmt

    INTEGER(KIND=4)             :: iscnt
    INTEGER(KIND=4)             :: iscnt2
    INTEGER(KIND=4),ALLOCATABLE :: ircnt(:), idisp(:)

    INTEGER(KIND=4)             :: n_mpi
    INTEGER(KIND=4)             :: n_mpi_tmp
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

    WRITE(unit_num, '(/, A15, I7, A, /)'),                    &
            '###  myrank = [',  myrank,                       &
            '] => 3Dim 32bit 실수형 변수를 병렬처리'

    !
    !--- point the 3D pointer variable to 3D target variable
    !
    WRITE(unit_num, '(A15, I7, 2(A17,3I7,A1))'),            &
            '###  myrank = [',  myrank,                     &
            '], shape(pSend)=[',   shape(pSend), ']',       &
            '   shape(pRecv)=[',   shape(pRecv), ']'

    ixTotSize = SIZE(pSend(1,:,1))
    iyTotSize = SIZE(pSend(1,1,:))

    WRITE(unit_num, '(A15, I7, 2(A32,I7,A1))'),                     &
            '###  myrank = [',  myrank,                             &
            '], SIZE(pSend(1,:,1))=ixTotSize=[',  ixTotSize, ']',   &  
            '], SIZE(pSend(1,1,:))=iyTotSize=[',  iyTotSize, ']'
    
    !
    !-- MPI variable 
    !
    ALLOCATE (ircnt(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Receive Count Allocation Error!!!'
            RETURN
         ENDIF
    ALLOCATE (idisp(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Display Count Allocation Error!!!'
            RETURN
         ENDIF

    !
    !-- Set send count, receive count, display position for MPI
    !
    iscnt   = SIZE(pSend(:, :, iy(1):iy(2)))        ! 마스터 프로세스(rank=0)로 전송할 데이터 개수
  
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

            n_mpi_tmp = i
            iscnt2   = iscnt2 + SIZE(pSend(:, :, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))

            !
            !-- set receive count for MPI
            !
            ircnt(i) = SIZE(pSend(:, :, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))
        ENDDO
    ENDIF

    !
    !-- Write MPI Log
    !
    WRITE(unit_num, '(A15, I7, A12, I7, A1, A21, I7, A1, I7, A2)'),     &
            '###  myrank = [', myrank, '], iscnt = [', iscnt, ']',      &
            ' => SIZE(pSend(:, :, ', iy(1), ':', iy(2), '))'
    
    IF ( myrank == 0 ) THEN
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      ircnt  = ['
        DO i=1, SIZE(ircnt)
            WRITE(unit_num, '(A14, I9)'), ' ', ircnt(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      idisp  = ['
        DO i=1, SIZE(idisp)
            WRITE(unit_num, '(A14, I9)'), ' ', idisp(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
    ENDIF

    WRITE(unit_num, '(A3, I5)') "---iy(1)", iy(1)

    !
    !-- Gather the output From other processes
    !
    CALL MPI_GATHERV(pSend(:, :, iy(1)), iscnt, MPI_REAL,            &
                     pRecv,              ircnt, idisp, MPI_REAL,     &
                     0, MPI_COMM_WORLD, errCode)

    retCode = errCode

    !
    !--- set pixel size calculated in defined total pixel size
    !
    ix_first = 1            
    ix_last  = ix(2)    

    !
    !
    !-- Write MPI Log
    !
    iy_size = iy(2) - iy(1)
    WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(F15.5))'

    IF (myrank==0) THEN
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(:, i, iy(1):iy(2)))
        ENDDO

        DO j=1, gci_nprocs

            n_mpi_tmp = j
            iy_size = iyy(n_mpi_tmp,2) - iyy(n_mpi_tmp,1)
            WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(F15.5))'
            DO i=ix_first, ix_last
                WRITE(unit_num, w_fmt),                                                 &
                        'myrank = [', n_mpi_tmp-1, '], irecv =',                        &
                        pRecv(:, i, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2))
            ENDDO
        ENDDO
    ELSE
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(:, i, iy(1):iy(2)))
        ENDDO
    ENDIF

    DEALLOCATE(ircnt)
    DEALLOCATE(idisp)

    RETURN
END SUBROUTINE GEMS_Share_MOD_MPI_Process_3D_R4


!
!--- 3Dim 64bit 실수형 변수를 병렬처리
!
SUBROUTINE GEMS_Share_MOD_MPI_Process_3D_R8(myrank, pSend, pRecv, ix, iyy, retCode)

! argument parameter
    INTEGER(KIND=4),          INTENT(IN)    :: myrank
    REAL(KIND=8),    POINTER, INTENT(IN)    :: pSend(:, :, :)
    REAL(KIND=8),    POINTER, INTENT(IN)    :: pRecv(:, :, :)
    INTEGER(KIND=4),          INTENT(IN)    :: ix(:)
    INTEGER(KIND=4),          INTENT(IN)    :: iyy(:, :)
    INTEGER(KIND=4),          INTENT(OUT)   :: retCode

! local variables
    INTEGER(KIND=4)             :: errCode = 0
    CHARACTER(LEN=120)          :: w_fmt

    INTEGER(KIND=4)             :: iscnt
    INTEGER(KIND=4)             :: iscnt2
    INTEGER(KIND=4),ALLOCATABLE :: ircnt(:), idisp(:)

    INTEGER(KIND=4)             :: n_mpi
    INTEGER(KIND=4)             :: n_mpi_tmp
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

    WRITE(unit_num, '(/, A15, I7, A, /)'),                    &
            '###  myrank = [',  myrank,                       &
            '] => 3Dim 64bit 실수형 변수를 병렬처리'
    !
    !--- point the 3D pointer variable to 3D target variable
    !
    WRITE(unit_num, '(A15, I7, 2(A17,3I7,A1))'),            &
            '###  myrank = [',  myrank,                     &
            '], shape(pSend)=[',   shape(pSend), ']',       &
            '   shape(pRecv)=[',   shape(pRecv), ']'

    ixTotSize = SIZE(pSend(1,:,1))
    iyTotSize = SIZE(pSend(1,1,:))

    WRITE(unit_num, '(A15, I7, 2(A32,I7,A1))'),                     &
            '###  myrank = [',  myrank,                             &
            '], SIZE(pSend(1,:,1))=ixTotSize=[',  ixTotSize, ']',   &  
            '], SIZE(pSend(1,1,:))=iyTotSize=[',  iyTotSize, ']'
    
    !
    !-- MPI variable 
    !
    ALLOCATE (ircnt(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Receive Count Allocation Error!!!'
            RETURN
         ENDIF
    ALLOCATE (idisp(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Display Count Allocation Error!!!'
            RETURN
         ENDIF

    !
    !-- Set send count, receive count, display position for MPI
    !
    iscnt   = SIZE(pSend(:, :, iy(1):iy(2)))        ! 마스터 프로세스(rank=0)로 전송할 데이터 개수
  
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

            n_mpi_tmp = i
            iscnt2   = iscnt2 + SIZE(pSend(:, :, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))

            !
            !-- set receive count for MPI
            !
            ircnt(i) = SIZE(pSend(:, :, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))
        ENDDO
    ENDIF

    !
    !-- Write MPI Log
    !
    WRITE(unit_num, '(A15, I7, A12, I7, A1, A21, I7, A1, I7, A2)'),     &
            '###  myrank = [', myrank, '], iscnt = [', iscnt, ']',      &
            ' => SIZE(pSend(:, :, ', iy(1), ':', iy(2), '))'
    
    IF ( myrank == 0 ) THEN
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      ircnt  = ['
        DO i=1, SIZE(ircnt)
            WRITE(unit_num, '(A14, I9)'), ' ', ircnt(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      idisp  = ['
        DO i=1, SIZE(idisp)
            WRITE(unit_num, '(A14, I9)'), ' ', idisp(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
    ENDIF

    WRITE(unit_num, '(A3, I5)') "---iy(1)", iy(1)

    !
    !-- Gather the output From other processes
    !
    CALL MPI_GATHERV(pSend(:, :, iy(1)), iscnt, MPI_DOUBLE,            &
                     pRecv,              ircnt, idisp, MPI_DOUBLE,     &
                     0, MPI_COMM_WORLD, errCode)

    retCode = errCode

    !
    !--- set pixel size calculated in defined total pixel size
    !
    ix_first = 1            
    ix_last  = ix(2)    

    !
    !
    !-- Write MPI Log
    !
    iy_size = iy(2) - iy(1)
    WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(F15.5))'

    IF (myrank==0) THEN
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(:, i, iy(1):iy(2)))
        ENDDO

        DO j=1, gci_nprocs

            n_mpi_tmp = j
            iy_size = iyy(n_mpi_tmp,2) - iyy(n_mpi_tmp,1)
            WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(F15.5))'
            DO i=ix_first, ix_last
                WRITE(unit_num, w_fmt),                                                 &
                        'myrank = [', n_mpi_tmp-1, '], irecv =',                        &
                        pRecv(:, i, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2))
            ENDDO
        ENDDO
    ELSE
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(:, i, iy(1):iy(2)))
        ENDDO
    ENDIF

    DEALLOCATE(ircnt)
    DEALLOCATE(idisp)

    RETURN
END SUBROUTINE GEMS_Share_MOD_MPI_Process_3D_R8


!
!--- 4Dim 8bit 정수형 변수를 병렬처리
!
SUBROUTINE GEMS_Share_MOD_MPI_Process_4D_I1(myrank, pSend, pRecv, ix, iyy, retCode)

! argument parameter
    INTEGER(KIND=4),          INTENT(IN)    :: myrank
    INTEGER(KIND=1), POINTER, INTENT(IN)    :: pSend(:, :, :, :)
    INTEGER(KIND=1), POINTER, INTENT(IN)    :: pRecv(:, :, :, :)
    INTEGER(KIND=4),          INTENT(IN)    :: ix(:)
    INTEGER(KIND=4),          INTENT(IN)    :: iyy(:, :)
    INTEGER(KIND=4),          INTENT(OUT)   :: retCode

! local variables
    INTEGER(KIND=4)             :: errCode = 0
    CHARACTER(LEN=120)          :: w_fmt

    INTEGER(KIND=4)             :: iscnt
    INTEGER(KIND=4)             :: iscnt2
    INTEGER(KIND=4),ALLOCATABLE :: ircnt(:), idisp(:)

    INTEGER(KIND=4)             :: n_mpi
    INTEGER(KIND=4)             :: n_mpi_tmp
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

    WRITE(unit_num, '(/, A15, I7, A, /)'),                 &
            '###  myrank = [',  myrank,                    &
            '] => 4Dim 8bit 정수형 변수를 병렬처리'

    !
    !--- point the 4D pointer variable to 4D target variable
    !
    WRITE(unit_num, '(A15, I7, 2(A17,4I7,A1))'),            &
            '###  myrank = [',  myrank,                     &
            '], shape(pSend)=[',   shape(pSend), ']',       &
            '   shape(pRecv)=[',   shape(pRecv), ']'

    ixTotSize = SIZE(pSend(1,1,:,1))
    iyTotSize = SIZE(pSend(1,1,1,:))

    WRITE(unit_num, '(A15, I7, 2(A32,I7,A1))'),                       &
            '###  myrank = [',  myrank,                               &
            '], SIZE(pSend(1,1,:,1))=ixTotSize=[',  ixTotSize, ']',   &  
            '], SIZE(pSend(1,1,1,:))=iyTotSize=[',  iyTotSize, ']'
    
    !
    !-- MPI variable 
    !
    ALLOCATE (ircnt(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Receive Count Allocation Error!!!'
            RETURN
         ENDIF
    ALLOCATE (idisp(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Display Count Allocation Error!!!'
            RETURN
         ENDIF

    !
    !-- Set send count, receive count, display position for MPI
    !
    iscnt   = SIZE(pSend(:, :, :, iy(1):iy(2)))        ! 마스터 프로세스(rank=0)로 전송할 데이터 개수
  
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

            n_mpi_tmp = i
            iscnt2   = iscnt2 + SIZE(pSend(:, :, :, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))

            !
            !-- set receive count for MPI
            !
            ircnt(i) = SIZE(pSend(:, :, :, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))
        ENDDO
    ENDIF

    !
    !-- Write MPI Log
    !
    WRITE(unit_num, '(A15, I7, A12, I7, A1, A24, I7, A1, I7, A2)'),     &
            '###  myrank = [', myrank, '], iscnt = [', iscnt, ']',      &
            ' => SIZE(pSend(:, :, :, ', iy(1), ':', iy(2), '))'
    
    IF ( myrank == 0 ) THEN
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      ircnt  = ['
        DO i=1, SIZE(ircnt)
            WRITE(unit_num, '(A14, I9)'), ' ', ircnt(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      idisp  = ['
        DO i=1, SIZE(idisp)
            WRITE(unit_num, '(A14, I9)'), ' ', idisp(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
    ENDIF

    WRITE(unit_num, '(A3, I5)') "---iy(1)", iy(1)

    !
    !-- Gather the output From other processes
    !
    CALL MPI_GATHERV(pSend(:,:,:,iy(1)), iscnt, MPI_CHARACTER,            &
                     pRecv,              ircnt, idisp, MPI_CHARACTER,     &
                     0, MPI_COMM_WORLD, errCode)

    retCode = errCode

    !
    !--- set pixel size calculated in defined total pixel size
    !
    ix_first = 1            
    ix_last  = ix(2)    

    !
    !
    !-- Write MPI Log
    !
    iy_size = iy(2) - iy(1)
    WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(I10))'

    IF (myrank==0) THEN
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(:, :, i, iy(1):iy(2)))
        ENDDO

        DO j=1, gci_nprocs

            n_mpi_tmp = j
            iy_size = iyy(n_mpi_tmp,2) - iyy(n_mpi_tmp,1)
            WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(I10))'
            DO i=ix_first, ix_last
                WRITE(unit_num, w_fmt),                                                 &
                        'myrank = [', n_mpi_tmp-1, '], irecv =',                        &
                        pRecv(:, :, i, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2))
            ENDDO
        ENDDO
    ELSE
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(:, :, i, iy(1):iy(2)))
        ENDDO
    ENDIF

    DEALLOCATE(ircnt)
    DEALLOCATE(idisp)

    RETURN
END SUBROUTINE GEMS_Share_MOD_MPI_Process_4D_I1


!
!--- 4Dim 16bit 정수형 변수를 병렬처리
!
SUBROUTINE GEMS_Share_MOD_MPI_Process_4D_I2(myrank, pSend, pRecv, ix, iyy, retCode)

! argument parameter
    INTEGER(KIND=4),          INTENT(IN)    :: myrank
    INTEGER(KIND=2), POINTER, INTENT(IN)    :: pSend(:, :, :, :)
    INTEGER(KIND=2), POINTER, INTENT(IN)    :: pRecv(:, :, :, :)
    INTEGER(KIND=4),          INTENT(IN)    :: ix(:)
    INTEGER(KIND=4),          INTENT(IN)    :: iyy(:, :)
    INTEGER(KIND=4),          INTENT(OUT)   :: retCode

! local variables
    INTEGER(KIND=4)             :: errCode = 0
    CHARACTER(LEN=120)          :: w_fmt

    INTEGER(KIND=4)             :: iscnt
    INTEGER(KIND=4)             :: iscnt2
    INTEGER(KIND=4),ALLOCATABLE :: ircnt(:), idisp(:)

    INTEGER(KIND=4)             :: n_mpi
    INTEGER(KIND=4)             :: n_mpi_tmp
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

    WRITE(unit_num, '(/, A15, I7, A, /)'),                 &
            '###  myrank = [',  myrank,                    &
            '] => 4Dim 16bit 정수형 변수를 병렬처리'

    !
    !--- point the 4D pointer variable to 4D target variable
    !
    WRITE(unit_num, '(A15, I7, 2(A17,4I7,A1))'),            &
            '###  myrank = [',  myrank,                     &
            '], shape(pSend)=[',   shape(pSend), ']',       &
            '   shape(pRecv)=[',   shape(pRecv), ']'

    ixTotSize = SIZE(pSend(1,1,:,1))
    iyTotSize = SIZE(pSend(1,1,1,:))

    WRITE(unit_num, '(A15, I7, 2(A32,I7,A1))'),                        &
            '###  myrank = [',  myrank,                                &
            '], SIZE(pSend(1,1,:,1))=ixTotSize=[',  ixTotSize, ']',     &  
            '], SIZE(pSend(1,1,1,:))=iyTotSize=[',  iyTotSize, ']'
    
    !
    !-- MPI variable 
    !
    ALLOCATE (ircnt(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Receive Count Allocation Error!!!'
            RETURN
         ENDIF
    ALLOCATE (idisp(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Display Count Allocation Error!!!'
            RETURN
         ENDIF

    !
    !-- Set send count, receive count, display position for MPI
    !
    iscnt   = SIZE(pSend(:, :, :, iy(1):iy(2)))        ! 마스터 프로세스(rank=0)로 전송할 데이터 개수
  
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

            n_mpi_tmp = i
            iscnt2   = iscnt2 + SIZE(pSend(:, :, :, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))

            !
            !-- set receive count for MPI
            !
            ircnt(i) = SIZE(pSend(:, :, :, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))
        ENDDO
    ENDIF

    !
    !-- Write MPI Log
    !
    WRITE(unit_num, '(A15, I7, A12, I7, A1, A24, I7, A1, I7, A2)'),     &
            '###  myrank = [', myrank, '], iscnt = [', iscnt, ']',      &
            ' => SIZE(pSend(:, :, :, ', iy(1), ':', iy(2), '))'

    IF ( myrank == 0 ) THEN
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      ircnt  = ['
        DO i=1, SIZE(ircnt)
            WRITE(unit_num, '(A14, I9)'), ' ', ircnt(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      idisp  = ['
        DO i=1, SIZE(idisp)
            WRITE(unit_num, '(A14, I9)'), ' ', idisp(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
    ENDIF

    WRITE(unit_num, '(A3, I5)') "---iy(1)", iy(1)

    !
    !-- Gather the output From other processes
    !
    CALL MPI_GATHERV(pSend(:,:,:,iy(1)), iscnt, MPI_SHORT,            &
                     pRecv,              ircnt, idisp, MPI_SHORT,     &
                     0, MPI_COMM_WORLD, errCode)

    retCode = errCode

    !
    !--- set pixel size calculated in defined total pixel size
    !
    ix_first = 1            
    ix_last  = ix(2)    

    !
    !
    !-- Write MPI Log
    !
    iy_size = iy(2) - iy(1)
    WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(I10))'

    IF (myrank==0) THEN
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(:, :, i, iy(1):iy(2)))
        ENDDO

        DO j=1, gci_nprocs

            n_mpi_tmp = j
            iy_size = iyy(n_mpi_tmp,2) - iyy(n_mpi_tmp,1)
            WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(I10))'
            DO i=ix_first, ix_last
                WRITE(unit_num, w_fmt),                                                 &
                        'myrank = [', n_mpi_tmp-1, '], irecv =',                        &
                        pRecv(:, :, i, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2))
            ENDDO
        ENDDO
    ELSE
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(:, :, i, iy(1):iy(2)))
        ENDDO
    ENDIF

    DEALLOCATE(ircnt)
    DEALLOCATE(idisp)

    RETURN
END SUBROUTINE GEMS_Share_MOD_MPI_Process_4D_I2

!
!--- 4Dim 32bit 정수형 변수를 병렬처리
!
SUBROUTINE GEMS_Share_MOD_MPI_Process_4D_I4(myrank, pSend, pRecv, ix, iyy, retCode)

! argument parameter
    INTEGER(KIND=4),          INTENT(IN)    :: myrank
    INTEGER(KIND=4), POINTER, INTENT(IN)    :: pSend(:, :, :, :)
    INTEGER(KIND=4), POINTER, INTENT(IN)    :: pRecv(:, :, :, :)
    INTEGER(KIND=4),          INTENT(IN)    :: ix(:)
    INTEGER(KIND=4),          INTENT(IN)    :: iyy(:, :)
    INTEGER(KIND=4),          INTENT(OUT)   :: retCode

! local variables
    INTEGER(KIND=4)             :: errCode = 0
    CHARACTER(LEN=120)          :: w_fmt

    INTEGER(KIND=4)             :: iscnt
    INTEGER(KIND=4)             :: iscnt2
    INTEGER(KIND=4),ALLOCATABLE :: ircnt(:), idisp(:)

    INTEGER(KIND=4)             :: n_mpi
    INTEGER(KIND=4)             :: n_mpi_tmp
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

    WRITE(unit_num, '(/, A15, I7, A, /)'),                 &
            '###  myrank = [',  myrank,                    &
            '] => 4Dim 32bit 정수형 변수를 병렬처리'

    !
    !--- point the 4D pointer variable to 4D target variable
    !
    WRITE(unit_num, '(A15, I7, 2(A17,4I7,A1))'),            &
            '###  myrank = [',  myrank,                     &
            '], shape(pSend)=[',   shape(pSend), ']',       &
            '   shape(pRecv)=[',   shape(pRecv), ']'

    ixTotSize = SIZE(pSend(1,1,:,1))
    iyTotSize = SIZE(pSend(1,1,1,:))

    WRITE(unit_num, '(A15, I7, 2(A32,I7,A1))'),                     &
            '###  myrank = [',  myrank,                             &
            '], SIZE(pSend(1,1,:,1))=ixTotSize=[',  ixTotSize, ']',   &  
            '], SIZE(pSend(1,1,1,:))=iyTotSize=[',  iyTotSize, ']'
    
    !
    !-- MPI variable 
    !
    ALLOCATE (ircnt(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Receive Count Allocation Error!!!'
            RETURN
         ENDIF
    ALLOCATE (idisp(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Display Count Allocation Error!!!'
            RETURN
         ENDIF

    !
    !-- Set send count, receive count, display position for MPI
    !
    iscnt   = SIZE(pSend(:, :, :, iy(1):iy(2)))        ! 마스터 프로세스(rank=0)로 전송할 데이터 개수
  
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

            n_mpi_tmp = i
            iscnt2   = iscnt2 + SIZE(pSend(:, :, :, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))

            !
            !-- set receive count for MPI
            !
            ircnt(i) = SIZE(pSend(:, :, :, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))
        ENDDO
    ENDIF

    !
    !-- Write MPI Log
    !
    WRITE(unit_num, '(A15, I7, A12, I7, A1, A24, I7, A1, I7, A2)'),     &
            '###  myrank = [', myrank, '], iscnt = [', iscnt, ']',      &
            ' => SIZE(pSend(:, :, :, ', iy(1), ':', iy(2), '))'
    
    IF ( myrank == 0 ) THEN
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      ircnt  = ['
        DO i=1, SIZE(ircnt)
            WRITE(unit_num, '(A14, I9)'), ' ', ircnt(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      idisp  = ['
        DO i=1, SIZE(idisp)
            WRITE(unit_num, '(A14, I9)'), ' ', idisp(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
    ENDIF

    WRITE(unit_num, '(A3, I5)') "---iy(1)", iy(1)

    !
    !-- Gather the output From other processes
    !
    CALL MPI_GATHERV(pSend(:,:,:,iy(1)), iscnt, MPI_INTEGER,            &
                     pRecv,              ircnt, idisp, MPI_INTEGER,     &
                     0, MPI_COMM_WORLD, errCode)

    retCode = errCode

    !
    !--- set pixel size calculated in defined total pixel size
    !
    ix_first = 1            
    ix_last  = ix(2)    

    !
    !
    !-- Write MPI Log
    !
    iy_size = iy(2) - iy(1)
    WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(I10))'

    IF (myrank==0) THEN
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(:, :, i, iy(1):iy(2)))
        ENDDO

        DO j=1, gci_nprocs

            n_mpi_tmp = j
            iy_size = iyy(n_mpi_tmp,2) - iyy(n_mpi_tmp,1)
            WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(I10))'
            DO i=ix_first, ix_last
                WRITE(unit_num, w_fmt),                                                 &
                        'myrank = [', n_mpi_tmp-1, '], irecv =',                        &
                        pRecv(:, :, i, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2))
            ENDDO
        ENDDO
    ELSE
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(:, :, i, iy(1):iy(2)))
        ENDDO
    ENDIF

    DEALLOCATE(ircnt)
    DEALLOCATE(idisp)

    RETURN
END SUBROUTINE GEMS_Share_MOD_MPI_Process_4D_I4

!
!--- 4Dim 32bit 실수형 변수를 병렬처리
!
SUBROUTINE GEMS_Share_MOD_MPI_Process_4D_R4(myrank, pSend, pRecv, ix, iyy, retCode)

! argument parameter
    INTEGER(KIND=4),          INTENT(IN)    :: myrank
    REAL(KIND=4),    POINTER, INTENT(IN)    :: pSend(:, :, :, :)
    REAL(KIND=4),    POINTER, INTENT(IN)    :: pRecv(:, :, :, :)
    INTEGER(KIND=4),          INTENT(IN)    :: ix(:)
    INTEGER(KIND=4),          INTENT(IN)    :: iyy(:, :)
    INTEGER(KIND=4),          INTENT(OUT)   :: retCode

! local variables
    INTEGER(KIND=4)             :: errCode = 0
    CHARACTER(LEN=120)          :: w_fmt

    INTEGER(KIND=4)             :: iscnt
    INTEGER(KIND=4)             :: iscnt2
    INTEGER(KIND=4),ALLOCATABLE :: ircnt(:), idisp(:)

    INTEGER(KIND=4)             :: n_mpi
    INTEGER(KIND=4)             :: n_mpi_tmp
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

    WRITE(unit_num, '(/, A15, I7, A, /)'),                    &
            '###  myrank = [',  myrank,                       &
            '] => 4Dim 32bit 실수형 변수를 병렬처리'

    !
    !--- point the 4D pointer variable to 4D target variable
    !
    WRITE(unit_num, '(A15, I7, 2(A17,4I7,A1))'),            &
            '###  myrank = [',  myrank,                     &
            '], shape(pSend)=[',   shape(pSend), ']',       &
            '   shape(pRecv)=[',   shape(pRecv), ']'

    ixTotSize = SIZE(pSend(1,1,:,1))
    iyTotSize = SIZE(pSend(1,1,1,:))

    WRITE(unit_num, '(A15, I7, 2(A32,I7,A1))'),                     &
            '###  myrank = [',  myrank,                             &
            '], SIZE(pSend(1,1,:,1))=ixTotSize=[',  ixTotSize, ']',   &  
            '], SIZE(pSend(1,1,1,:))=iyTotSize=[',  iyTotSize, ']'
    
    !
    !-- MPI variable 
    !
    ALLOCATE (ircnt(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Receive Count Allocation Error!!!'
            RETURN
         ENDIF
    ALLOCATE (idisp(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Display Count Allocation Error!!!'
            RETURN
         ENDIF

    !
    !-- Set send count, receive count, display position for MPI
    !
    iscnt   = SIZE(pSend(:, :, :, iy(1):iy(2)))        ! 마스터 프로세스(rank=0)로 전송할 데이터 개수
  
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

            n_mpi_tmp = i
            iscnt2   = iscnt2 + SIZE(pSend(:, :, :, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))

            !
            !-- set receive count for MPI
            !
            ircnt(i) = SIZE(pSend(:, :, :, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))
        ENDDO
    ENDIF

    !
    !-- Write MPI Log
    !
    WRITE(unit_num, '(A15, I7, A12, I7, A1, A24, I7, A1, I7, A2)'),     &
            '###  myrank = [', myrank, '], iscnt = [', iscnt, ']',      &
            ' => SIZE(pSend(:, :, :, ', iy(1), ':', iy(2), '))'
    
    IF ( myrank == 0 ) THEN
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      ircnt  = ['
        DO i=1, SIZE(ircnt)
            WRITE(unit_num, '(A14, I9)'), ' ', ircnt(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      idisp  = ['
        DO i=1, SIZE(idisp)
            WRITE(unit_num, '(A14, I9)'), ' ', idisp(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
    ENDIF

    WRITE(unit_num, '(A3, I5)') "---iy(1)", iy(1)

    !
    !-- Gather the output From other processes
    !
    CALL MPI_GATHERV(pSend(:,:,:,iy(1)), iscnt, MPI_REAL,            &
                     pRecv,              ircnt, idisp, MPI_REAL,     &
                     0, MPI_COMM_WORLD, errCode)

    retCode = errCode

    !
    !--- set pixel size calculated in defined total pixel size
    !
    ix_first = 1            
    ix_last  = ix(2)    

    !
    !
    !-- Write MPI Log
    !
    iy_size = iy(2) - iy(1)
    WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(F15.5))'

    IF (myrank==0) THEN
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(:, :, i, iy(1):iy(2)))
        ENDDO

        DO j=1, gci_nprocs

            n_mpi_tmp = j
            iy_size = iyy(n_mpi_tmp,2) - iyy(n_mpi_tmp,1)
            WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(F15.5))'
            DO i=ix_first, ix_last
                WRITE(unit_num, w_fmt),                                                 &
                        'myrank = [', n_mpi_tmp-1, '], irecv =',                        &
                        pRecv(:, :, i, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2))
            ENDDO
        ENDDO
    ELSE
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(:, :, i, iy(1):iy(2)))
        ENDDO
    ENDIF

    DEALLOCATE(ircnt)
    DEALLOCATE(idisp)

    RETURN
END SUBROUTINE GEMS_Share_MOD_MPI_Process_4D_R4


!
!--- 4Dim 64bit 실수형 변수를 병렬처리
!
SUBROUTINE GEMS_Share_MOD_MPI_Process_4D_R8(myrank, pSend, pRecv, ix, iyy, retCode)

! argument parameter
    INTEGER(KIND=4),          INTENT(IN)    :: myrank
    REAL(KIND=8),    POINTER, INTENT(IN)    :: pSend(:, :, :, :)
    REAL(KIND=8),    POINTER, INTENT(IN)    :: pRecv(:, :, :, :)
    INTEGER(KIND=4),          INTENT(IN)    :: ix(:)
    INTEGER(KIND=4),          INTENT(IN)    :: iyy(:, :)
    INTEGER(KIND=4),          INTENT(OUT)   :: retCode

! local variables
    INTEGER(KIND=4)             :: errCode = 0
    CHARACTER(LEN=120)          :: w_fmt

    INTEGER(KIND=4)             :: iscnt
    INTEGER(KIND=4)             :: iscnt2
    INTEGER(KIND=4),ALLOCATABLE :: ircnt(:), idisp(:)

    INTEGER(KIND=4)             :: n_mpi
    INTEGER(KIND=4)             :: n_mpi_tmp
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

    WRITE(unit_num, '(/, A15, I7, A, /)'),                    &
            '###  myrank = [',  myrank,                       &
            '] => 4Dim 64bit 실수형 변수를 병렬처리'
    !
    !--- point the 4D pointer variable to 4D target variable
    !
    WRITE(unit_num, '(A15, I7, 2(A17,4I7,A1))'),            &
            '###  myrank = [',  myrank,                     &
            '], shape(pSend)=[',   shape(pSend), ']',       &
            '   shape(pRecv)=[',   shape(pRecv), ']'

    ixTotSize = SIZE(pSend(1,1,:,1))
    iyTotSize = SIZE(pSend(1,1,1,:))

    WRITE(unit_num, '(A15, I7, 2(A32,I7,A1))'),                      &
            '###  myrank = [',  myrank,                              &
            '], SIZE(pSend(1,1,:,1))=ixTotSize=[',  ixTotSize, ']',  &  
            '], SIZE(pSend(1,1,1,:))=iyTotSize=[',  iyTotSize, ']'
    
    !
    !-- MPI variable 
    !
    ALLOCATE (ircnt(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Receive Count Allocation Error!!!'
            RETURN
         ENDIF
    ALLOCATE (idisp(gci_nprocs), STAT = errCode)    ! 프로세스 개수 만큼의 방이 필요
         IF(errCode .NE. 0) THEN
            PRINT*, 'MPI Display Count Allocation Error!!!'
            RETURN
         ENDIF

    !
    !-- Set send count, receive count, display position for MPI
    !
    iscnt   = SIZE(pSend(:, :, :, iy(1):iy(2)))        ! 마스터 프로세스(rank=0)로 전송할 데이터 개수
  
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

            n_mpi_tmp = i
            iscnt2   = iscnt2 + SIZE(pSend(:, :, :, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))

            !
            !-- set receive count for MPI
            !
            ircnt(i) = SIZE(pSend(:, :, :, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2)))
        ENDDO
    ENDIF

    !
    !-- Write MPI Log
    !
    WRITE(unit_num, '(A15, I7, A12, I7, A1, A24, I7, A1, I7, A2)'),     &
            '###  myrank = [', myrank, '], iscnt = [', iscnt, ']',      &
            ' => SIZE(pSend(:, :, :, ', iy(1), ':', iy(2), '))'
    
    IF ( myrank == 0 ) THEN
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      ircnt  = ['
        DO i=1, SIZE(ircnt)
            WRITE(unit_num, '(A14, I9)'), ' ', ircnt(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
        WRITE(unit_num, '(A15, I7, A1)'),                               &
                '###  myrank = [', myrank, ']'
        WRITE(unit_num, '(A15, I9, A1)'),                               &
                '      idisp  = ['
        DO i=1, SIZE(idisp)
            WRITE(unit_num, '(A14, I9)'), ' ', idisp(i)
        ENDDO
        WRITE(unit_num, '(A6)'), '     ]'
    ENDIF

    WRITE(unit_num, '(A3, I5)') "---iy(1)", iy(1)

    !
    !-- Gather the output From other processes
    !
    CALL MPI_GATHERV(pSend(:,:,:,iy(1)), iscnt, MPI_DOUBLE,            &
                     pRecv,              ircnt, idisp, MPI_DOUBLE,     &
                     0, MPI_COMM_WORLD, errCode)

    retCode = errCode

    !
    !--- set pixel size calculated in defined total pixel size
    !
    ix_first = 1            
    ix_last  = ix(2)    

    !
    !
    !-- Write MPI Log
    !
    iy_size = iy(2) - iy(1)
    WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(F15.5))'

    IF (myrank==0) THEN
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(:, :, i, iy(1):iy(2)))
        ENDDO

        DO j=1, gci_nprocs

            n_mpi_tmp = j
            iy_size = iyy(n_mpi_tmp,2) - iyy(n_mpi_tmp,1)
            WRITE(w_fmt, '(A14, I7, A8)'), '(A10, I7, A11,', (iy_size+1),'(F15.5))'
            DO i=ix_first, ix_last
                WRITE(unit_num, w_fmt),                                                 &
                        'myrank = [', n_mpi_tmp-1, '], irecv =',                        &
                        pRecv(:, :, i, iyy(n_mpi_tmp,1):iyy(n_mpi_tmp,2))
            ENDDO
        ENDDO
    ELSE
        DO i=ix_first, ix_last
            WRITE(unit_num, w_fmt),                                                     &
                        'myrank = [', myrank, '], isend =',                             &
                        (pSend(:, :, i, iy(1):iy(2)))
        ENDDO
    ENDIF

    DEALLOCATE(ircnt)
    DEALLOCATE(idisp)

    RETURN
END SUBROUTINE GEMS_Share_MOD_MPI_Process_4D_R8


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
