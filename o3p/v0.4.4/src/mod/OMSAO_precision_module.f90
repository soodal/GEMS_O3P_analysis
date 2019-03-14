MODULE OMSAO_precision_module

  IMPLICIT NONE

  ! =====================================================
  ! Define KIND variables for single and double precision
  ! =====================================================

  INTEGER, PARAMETER :: i1 = 1
  INTEGER, PARAMETER :: i2 = 2
  INTEGER, PARAMETER :: i3 = 3
  INTEGER, PARAMETER :: i4 = 4
  INTEGER, PARAMETER :: i8 = 8

  INTEGER, PARAMETER :: sp = KIND(1.0)
  INTEGER, PARAMETER :: dp = KIND(1.0D0) !dp = KIND(1.0)
  INTEGER, PARAMETER :: r4 = KIND(1.0)
  INTEGER, PARAMETER :: r8 = KIND(1.0D0)

END MODULE OMSAO_precision_module

