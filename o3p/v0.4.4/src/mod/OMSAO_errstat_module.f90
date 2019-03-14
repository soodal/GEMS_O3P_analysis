MODULE OMSAO_errstat_module

  ! =================================================================
  !
  ! This module defines variables associated with error handling. It
  ! also loads/includes all (SDPTK) files that define error messages
  ! and generally deal with error handling.
  !
  ! =================================================================

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: maxchlen

  IMPLICIT NONE

  ! Toolkit include files (sytem)
!!  INCLUDE 'PGS_SMF.f'
!!  INCLUDE 'PGS_PC.f'
!!  INCLUDE 'PGS_PC_9.f'
!!  INCLUDE 'PGS_IO.f'

  ! Toolkit include files (PGE specific)
  INCLUDE 'PGS_OMSAO_52500.f'

  ! OMI L1B reader include file
!  INCLUDE 'PGS_OMI_1900.f'

  ! ----------------------------------------------------
  ! String carrying a potential error or warning message
  ! ----------------------------------------------------
  CHARACTER (LEN=maxchlen) :: pge_errstat_msg

  ! ----------------
  ! Some error stati
  ! ----------------
  INTEGER, PARAMETER :: &
       pge_errstat_ok = 0, pge_errstat_warning  = 1, &
       pge_errstat_error = 2, pge_errstat_fatal = 3
  
  ! -----------------------------------------------------------------
  ! Status variables for READ process; used to identify status of
  ! current READ requests, like "O.K.", "FAILED", "END OF DATA", etc.
  ! -----------------------------------------------------------------
  INTEGER, PARAMETER :: &
       file_read_ok = 0, file_read_failed = 1, file_read_missing = 3, &
       file_read_eof = -1

  ! --------------------------------------------
  ! Status parameters for HE5 interface routines.
  ! --------------------------------------------
  INTEGER, PARAMETER :: he5_stat_ok = 0, he5_stat_fail = -1



END MODULE OMSAO_errstat_module
