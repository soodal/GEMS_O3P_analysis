!-------------------------------------------------------------------------------
!+Module to log the processing of GEMS Algorithm 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE Share_MOD_Log

!-------------------------------------------------------------------------------
!+Description: 
!      GEMS Algorithm Program Execution Logging.
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2014.09.12 Fisrt Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Modules:

!-------------------------------------------------------------------------------
!+INTERFACE:

  INTERFACE GEMS_Share_MOD_Log
  
      MODULE PROCEDURE GEMS_Share_MOD_Log_F32,          &
                       GEMS_Share_MOD_Log_F32Array1,    &
                       GEMS_Share_MOD_Log_F32Array2,    &
                       GEMS_Share_MOD_Log_F32Array3,    &
                       GEMS_Share_MOD_Log_F32Array3_1,  &
                       GEMS_Share_MOD_Log_F32Array3_2,  &
                       GEMS_Share_MOD_Log_F32Array4,    &
                       GEMS_Share_MOD_Log_F64,          &
                       GEMS_Share_MOD_Log_F64Array1,    &
                       GEMS_Share_MOD_Log_F64Array2,    &
                       GEMS_Share_MOD_Log_F64Array3,    &
                       GEMS_Share_MOD_Log_F64Array3_1,  &
                       GEMS_Share_MOD_Log_F64Array3_2,  &
                       GEMS_Share_MOD_Log_F64Array4,    &
                       GEMS_Share_MOD_Log_I8,           &
                       GEMS_Share_MOD_Log_I8Array1,     &
                       GEMS_Share_MOD_Log_I8Array2,     &
                       GEMS_Share_MOD_Log_I8Array3,     &
                       GEMS_Share_MOD_Log_I8Array3_1,   &
                       GEMS_Share_MOD_Log_I8Array4,     &
                       GEMS_Share_MOD_Log_I8Array4_1,   &
                       GEMS_Share_MOD_Log_I16,          &
                       GEMS_Share_MOD_Log_I16Array1,    &
                       GEMS_Share_MOD_Log_I16Array2,    &
                       GEMS_Share_MOD_Log_I16Array2_1,  &
                       GEMS_Share_MOD_Log_I16Array3,    &
                       GEMS_Share_MOD_Log_I16Array3_1,  &
                       GEMS_Share_MOD_Log_I16Array4,    &
                       GEMS_Share_MOD_Log_I16Array4_1,  &
                       GEMS_Share_MOD_Log_I32,          &
                       GEMS_Share_MOD_Log_I32Array1,    &
                       GEMS_Share_MOD_Log_I32Array2,    &
                       GEMS_Share_MOD_Log_I32Array3,    &
                       GEMS_Share_MOD_Log_I32Array3_1,  &
                       GEMS_Share_MOD_Log_I64,          &
                       GEMS_Share_MOD_Log_I64Array1,    &
                       GEMS_Share_MOD_Log_I64Array2,    &
                       GEMS_Share_MOD_Log_I64Array3,    &
                       GEMS_Share_MOD_Log_I64Array3_1,  &
                       GEMS_Share_MOD_Log_str,          &
                       GEMS_Share_MOD_Log_str2

  END INTERFACE

!-------------------------------------------------------------------------------
!+Global Variable:

    !
    ! External Function Declaration
    !
    INTEGER  lprint
    external lprint

!--------------------------------
    !
    ! External Function Declaration For Using GEMS HDF5 API
    !
!    INTEGER  GEMS_Share_Hdf5ReadData
!    INTEGER  GEMS_Share_Hdf5WriteData
!    INTEGER  GEMS_Share_Hdf5ReadAttr
!    INTEGER  GEMS_Share_Hdf5WriteAttr
!    INTEGER  GEMS_Share_Hdf5CreateL2File
!    INTEGER  GEMS_Share_Hdf5CreateL1BFile
!    external GEMS_Share_Hdf5ReadData
!    external GEMS_Share_Hdf5WriteData
!    external GEMS_Share_Hdf5ReadAttr
!    external GEMS_Share_Hdf5WriteAttr
!    external GEMS_Share_Hdf5CreateL2File
!    external GEMS_Share_Hdf5CreateL1BFile
!--------------------------------


  CONTAINS 

!-------------------------------------------------------------------------------
!+SUBROUTINE:

!
! ---------- GEMS_Share_MOD_Log_F32
!
    SUBROUTINE GEMS_Share_MOD_Log_F32(llvl, F32, str)
    integer          :: llvl
    real(KIND=4)     :: F32
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_f32(llvl, F32, str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_F32

!
! ---------- GEMS_Share_MOD_Log_F32Array1
!
    SUBROUTINE GEMS_Share_MOD_Log_F32Array1(llvl, F32_dim1, sz, str)
    integer          :: llvl
    real(KIND=4), DIMENSION(:)  :: F32_dim1
    integer          :: sz
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_f32array1(llvl, F32_dim1, sz, str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_F32Array1

!
! ---------- GEMS_Share_MOD_Log_F32Array2
!
    SUBROUTINE GEMS_Share_MOD_Log_F32Array2(llvl, F32_dim2, isz, jsz, str)
    integer          :: llvl
    real(KIND=4), DIMENSION(:, :)  :: F32_dim2
    integer          :: isz
    integer          :: jsz
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_f32array2(llvl, F32_dim2, jsz, isz, str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_F32Array2

!
! ---------- GEMS_Share_MOD_Log_F32Array3
!
    SUBROUTINE GEMS_Share_MOD_Log_F32Array3(llvl, F32_dim3, isz, jsz, ksz, str)
    integer          :: llvl
    real(KIND=4), DIMENSION(:, :, :)  :: F32_dim3
    integer          :: isz
    integer          :: jsz
    integer          :: ksz
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_f32array3(llvl, F32_dim3, ksz, jsz, isz, str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_F32Array3

!
! ---------- GEMS_Share_MOD_Log_F32Array3_1
!
    SUBROUTINE GEMS_Share_MOD_Log_F32Array3_1(llvl, F32_dim3, isz, jsz, ksz, po, str)
    integer          :: llvl
    real(KIND=4), DIMENSION(:, :, :)  :: F32_dim3
    integer          :: isz
    integer          :: jsz
    integer          :: ksz
    integer          :: po
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_f32array3_1(llvl, F32_dim3, ksz, jsz, isz, po, str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_F32Array3_1

!
! ---------- GEMS_Share_MOD_Log_F32Array3_2
!
    SUBROUTINE GEMS_Share_MOD_Log_F32Array3_2(llvl, F32_dim3, isz, jsz, ksz, fo, str)
    integer          :: llvl
    real(KIND=4), DIMENSION(:, :, :)  :: F32_dim3
    integer          :: isz
    integer          :: jsz
    integer          :: ksz
    character(LEN=*) :: fo
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_f32array3_2(llvl, F32_dim3, ksz, jsz, isz, fo, str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_F32Array3_2

!
! ---------- GEMS_Share_MOD_Log_F32Array4
!
    SUBROUTINE GEMS_Share_MOD_Log_F32Array4(llvl, F32_dim4, isz, jsz, ksz, lsz, str)
    integer          :: llvl
    real(KIND=4), DIMENSION(:, :, :,:)  :: F32_dim4
    integer          :: isz
    integer          :: jsz
    integer          :: ksz
    integer          :: lsz
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_f32array4(llvl, F32_dim4, lsz, ksz, jsz, isz, str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_F32Array4

!
!
! ---------- GEMS_Share_MOD_Log_F64
!
    SUBROUTINE GEMS_Share_MOD_Log_F64(llvl, F64, str)
    integer          :: llvl
    real(KIND=8)     :: F64
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_f64(llvl, F64, str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_F64

!
! ---------- GEMS_Share_MOD_Log_F64Array1
!
    SUBROUTINE GEMS_Share_MOD_Log_F64Array1(llvl, F64_dim1, sz, str)
    integer          :: llvl
    real(KIND=8), DIMENSION(:)  :: F64_dim1
    integer          :: sz
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_f64array1(llvl, F64_dim1, sz, str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_F64Array1

!
! ---------- GEMS_Share_MOD_Log_F64Array2
!
    SUBROUTINE GEMS_Share_MOD_Log_F64Array2(llvl, F64_dim2, isz, jsz, str)
    integer          :: llvl
    real(KIND=8), DIMENSION(:, :)  :: F64_dim2
    integer          :: isz
    integer          :: jsz
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_f64array2(llvl, F64_dim2, jsz, isz, str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_F64Array2

!
! ---------- GEMS_Share_MOD_Log_F64Array3
!
    SUBROUTINE GEMS_Share_MOD_Log_F64Array3(llvl, F64_dim3, isz, jsz, ksz, str)
    integer          :: llvl
    real(KIND=8), DIMENSION(:, :, :)  :: F64_dim3
    integer          :: isz
    integer          :: jsz
    integer          :: ksz
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_f64array3(llvl, F64_dim3, ksz, jsz, isz, str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_F64Array3

!
! ---------- GEMS_Share_MOD_Log_F64Array3_1
!
    SUBROUTINE GEMS_Share_MOD_Log_F64Array3_1(llvl, F64_dim3, isz, jsz, ksz, po, str)
    integer          :: llvl
    real(KIND=8), DIMENSION(:, :, :)  :: F64_dim3
    integer          :: isz
    integer          :: jsz
    integer          :: ksz
    integer          :: po
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_f64array3_1(llvl, F64_dim3, ksz, jsz, isz, po, str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_F64Array3_1

!
! ---------- GEMS_Share_MOD_Log_F64Array3_2
!
    SUBROUTINE GEMS_Share_MOD_Log_F64Array3_2(llvl, F64_dim3, isz, jsz, ksz, fo, str)
    integer          :: llvl
    real(KIND=8), DIMENSION(:, :, :)  :: F64_dim3
    integer          :: isz
    integer          :: jsz
    integer          :: ksz
    character(LEN=*) :: fo
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_f64array3_2(llvl, F64_dim3, ksz, jsz, isz, fo, str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_F64Array3_2

!
! ---------- GEMS_Share_MOD_Log_F64Array4
!
    SUBROUTINE GEMS_Share_MOD_Log_F64Array4(llvl, F64_dim4, isz, jsz, ksz, lsz, str)
    integer          :: llvl
    real(KIND=8), DIMENSION(:, :, :,:)  :: F64_dim4
    integer          :: isz
    integer          :: jsz
    integer          :: ksz
    integer          :: lsz
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_f64array4(llvl, F64_dim4, lsz, ksz, jsz, isz, str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_F64Array4

!
! ---------- GEMS_Share_MOD_Log_I8
!
    SUBROUTINE GEMS_Share_MOD_Log_I8(llvl, I8, str)
    integer          :: llvl
    integer(KIND=1)  :: I8
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_i8(llvl, I8, str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_I8

!
! ---------- GEMS_Share_MOD_Log_I8Array1
!
    SUBROUTINE GEMS_Share_MOD_Log_I8Array1(llvl, I8_dim1, sz, str)
    integer          :: llvl
    integer(KIND=1), DIMENSION(:)  :: I8_dim1
    integer          :: sz
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_i8array1(llvl, I8_dim1, sz, str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_I8Array1

!
! ---------- GEMS_Share_MOD_Log_I8Array2
!
    SUBROUTINE GEMS_Share_MOD_Log_I8Array2(llvl, I8_dim2, isz, jsz, str)
    integer          :: llvl
    integer(KIND=1), DIMENSION(:, :)  :: I8_dim2
    integer          :: isz
    integer          :: jsz
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_i8array2(llvl, I8_dim2, jsz, isz, str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_I8Array2

!
! ---------- GEMS_Share_MOD_Log_I8Array3
!
    SUBROUTINE GEMS_Share_MOD_Log_I8Array3(llvl, I8_dim3, isz, jsz, ksz, str)
    integer          :: llvl
    integer(KIND=1), DIMENSION(:, :, :)  :: I8_dim3
    integer          :: isz
    integer          :: jsz
    integer          :: ksz
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_i8array3(llvl, I8_dim3, ksz, jsz, isz, str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_I8Array3

!
! ---------- GEMS_Share_MOD_Log_I8Array3_1
!
    SUBROUTINE GEMS_Share_MOD_Log_I8Array3_1(llvl, I8_dim3, isz, jsz, ksz, po, str)
    integer          :: llvl
    integer(KIND=1), DIMENSION(:, :, :)  :: I8_dim3
    integer          :: isz
    integer          :: jsz
    integer          :: ksz
    integer          :: po      ! option to print, if po = -1, print all DATAs.
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_i8array3_1(llvl, I8_dim3, ksz, jsz, isz, po, str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_I8Array3_1

!
! ---------- GEMS_Share_MOD_Log_I8Array4
!
    SUBROUTINE GEMS_Share_MOD_Log_I8Array4(llvl, I8_dim4, isz, jsz, ksz, lsz, str)
    integer          :: llvl
    integer(KIND=1), DIMENSION(:, :, :, :)  :: I8_dim4
    integer          :: isz
    integer          :: jsz
    integer          :: ksz
    integer          :: lsz
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_i8array4(llvl, I8_dim4, lsz, ksz, jsz, isz, str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_I8Array4

!
! ---------- GEMS_Share_MOD_Log_I8Array4_1
!
    SUBROUTINE GEMS_Share_MOD_Log_I8Array4_1(llvl, I8_dim4, isz, jsz, ksz, lsz, po, str)
    integer          :: llvl
    integer(KIND=1), DIMENSION(:, :, :, :)  :: I8_dim4
    integer          :: isz
    integer          :: jsz
    integer          :: ksz
    integer          :: lsz
    integer          :: po      ! option to print, if po = -1, print all DATAs.
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_i8array4_1(llvl, I8_dim4, ksz, jsz, isz, po, str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_I8Array4_1

!
! ---------- GEMS_Share_MOD_Log_I16
!
    SUBROUTINE GEMS_Share_MOD_Log_I16(llvl, I16, str)
    integer          :: llvl
    integer(KIND=2)  :: I16
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_i16(llvl, I16, str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_I16

!
! ---------- GEMS_Share_MOD_Log_I16Array1
!
    SUBROUTINE GEMS_Share_MOD_Log_I16Array1(llvl, I16_dim1, sz, str)
    integer          :: llvl
    integer(KIND=2), DIMENSION(:)  :: I16_dim1
    integer          :: sz
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_i16array1(llvl, I16_dim1, sz, str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_I16Array1

!
! ---------- GEMS_Share_MOD_Log_I16Array2
!
    SUBROUTINE GEMS_Share_MOD_Log_I16Array2(llvl, I16_dim2, isz, jsz, str)
    integer          :: llvl
    integer(KIND=2), DIMENSION(:, :)  :: I16_dim2
    integer          :: isz
    integer          :: jsz
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_i16array2(llvl, I16_dim2, jsz, isz, str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_I16Array2

!
! ---------- GEMS_Share_MOD_Log_I16Array2_1
!
    SUBROUTINE GEMS_Share_MOD_Log_I16Array2_1(llvl, I16_dim2, isz, jsz, so, str)
    integer          :: llvl
    integer(KIND=2), DIMENSION(:, :)  :: I16_dim2
    integer          :: isz
    integer          :: jsz
    integer          :: so      ! if so = 1, this number type is signed number type.
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_i16array2_1(llvl, I16_dim2, jsz, isz, so, str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_I16Array2_1

!
! ---------- GEMS_Share_MOD_Log_I16Array3
!
    SUBROUTINE GEMS_Share_MOD_Log_I16Array3(llvl, I16_dim3, isz, jsz, ksz, str)
    integer          :: llvl
    integer(KIND=2), DIMENSION(:, :, :)  :: I16_dim3
    integer          :: isz
    integer          :: jsz
    integer          :: ksz
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_i16array3(llvl, I16_dim3, ksz, jsz, isz,  str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_I16Array3

!
! ---------- GEMS_Share_MOD_Log_I16Array3_1
!
    SUBROUTINE GEMS_Share_MOD_Log_I16Array3_1(llvl, I16_dim3, isz, jsz, ksz, po, str)
    integer          :: llvl
    integer(KIND=2), DIMENSION(:, :, :)  :: I16_dim3
    integer          :: isz
    integer          :: jsz
    integer          :: ksz
    integer          :: po      ! option to print, if po = -1, print all DATAs.
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_i16array3_1(llvl, I16_dim3, ksz, jsz, isz,  po, str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_I16Array3_1

!
! ---------- GEMS_Share_MOD_Log_I16Array4
!
    SUBROUTINE GEMS_Share_MOD_Log_I16Array4(llvl, I16_dim4, isz, jsz, ksz, lsz, str)
    integer          :: llvl
    integer(KIND=2), DIMENSION(:, :, :, :)  :: I16_dim4
    integer          :: isz
    integer          :: jsz
    integer          :: ksz
    integer          :: lsz
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_i16array4(llvl, I16_dim4, lsz, ksz, jsz, isz,  str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_I16Array4

!
! ---------- GEMS_Share_MOD_Log_I16Array4_1
!
    SUBROUTINE GEMS_Share_MOD_Log_I16Array4_1(llvl, I16_dim4, isz, jsz, ksz, lsz, po, str)
    integer          :: llvl
    integer(KIND=2), DIMENSION(:, :, :, :)  :: I16_dim4
    integer          :: isz
    integer          :: jsz
    integer          :: ksz
    integer          :: lsz
    integer          :: po      ! option to print, if po = -1, print all DATAs.
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_i16array4_1(llvl, I16_dim4, lsz, ksz, jsz, isz,  po, str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_I16Array4_1

!
! ---------- GEMS_Share_MOD_Log_I32
!
    SUBROUTINE GEMS_Share_MOD_Log_I32(llvl, I32, str)
    integer          :: llvl
    integer(KIND=4)  :: I32
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_i32(llvl, I32, str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_I32

!
! ---------- GEMS_Share_MOD_Log_I32Array1
!
    SUBROUTINE GEMS_Share_MOD_Log_I32Array1(llvl, I32_dim1, sz, str)
    integer          :: llvl
    integer(KIND=4), DIMENSION(:)  :: I32_dim1
    integer          :: sz
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_i32array1(llvl, I32_dim1, sz, str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_I32Array1

!
! ---------- GEMS_Share_MOD_Log_I32Array2
!
    SUBROUTINE GEMS_Share_MOD_Log_I32Array2(llvl, I32_dim2, isz, jsz, str)
    integer          :: llvl
    integer(KIND=4), DIMENSION(:, :)  :: I32_dim2
    integer          :: isz
    integer          :: jsz
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_i32array2(llvl, I32_dim2, jsz, isz, str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_I32Array2

!
! ---------- GEMS_Share_MOD_Log_I32Array3
!
    SUBROUTINE GEMS_Share_MOD_Log_I32Array3(llvl, I32_dim3, isz, jsz, ksz, str)
    integer          :: llvl
    integer(KIND=4), DIMENSION(:, :, :)  :: I32_dim3
    integer          :: isz
    integer          :: jsz
    integer          :: ksz
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_i32array3(llvl, I32_dim3, ksz, jsz, isz,  str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_I32Array3

!
! ---------- GEMS_Share_MOD_Log_I32Array3_1
!
    SUBROUTINE GEMS_Share_MOD_Log_I32Array3_1(llvl, I32_dim3, isz, jsz, ksz, po, str)
    integer          :: llvl
    integer(KIND=4), DIMENSION(:, :, :)  :: I32_dim3
    integer          :: isz
    integer          :: jsz
    integer          :: ksz
    integer          :: po      ! option to print, if po = -1, print all DATAs.
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_i32array3_1(llvl, I32_dim3, ksz, jsz, isz,  po, str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_I32Array3_1

!
! ---------- GEMS_Share_MOD_Log_I64
!
    SUBROUTINE GEMS_Share_MOD_Log_I64(llvl, I64, str)
    integer          :: llvl
    integer(KIND=8)  :: I64
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_i64(llvl, I64, str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_I64

!
! ---------- GEMS_Share_MOD_Log_I64Array1
!
    SUBROUTINE GEMS_Share_MOD_Log_I64Array1(llvl, I64_dim1, sz, str)
    integer          :: llvl
    integer(KIND=8), DIMENSION(:)  :: I64_dim1
    integer          :: sz
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_i64array1(llvl, I64_dim1, sz, str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_I64Array1

!
! ---------- GEMS_Share_MOD_Log_I64Array2
!
    SUBROUTINE GEMS_Share_MOD_Log_I64Array2(llvl, I64_dim2, isz, jsz, str)
    integer          :: llvl
    integer(KIND=8), DIMENSION(:, :)  :: I64_dim2
    integer          :: isz
    integer          :: jsz
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_i64array2(llvl, I64_dim2, jsz, isz, str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_I64Array2

!
! ---------- GEMS_Share_MOD_Log_I64Array3
!
    SUBROUTINE GEMS_Share_MOD_Log_I64Array3(llvl, I64_dim3, isz, jsz, ksz, str)
    integer          :: llvl
    integer(KIND=8), DIMENSION(:, :, :)  :: I64_dim3
    integer          :: isz
    integer          :: jsz
    integer          :: ksz
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_i64array3(llvl, I64_dim3, ksz, jsz, isz,  str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_I64Array3

!
! ---------- GEMS_Share_MOD_Log_I64Array3_1
!
    SUBROUTINE GEMS_Share_MOD_Log_I64Array3_1(llvl, I64_dim3, isz, jsz, ksz, po, str)
    integer          :: llvl
    integer(KIND=8), DIMENSION(:, :, :)  :: I64_dim3
    integer          :: isz
    integer          :: jsz
    integer          :: ksz
    integer          :: po      ! option to print, if po = -1, print all DATAs.
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_i64array3_1(llvl, I64_dim3, ksz, jsz, isz,  po, str)

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_I64Array3_1

!
! ---------- GEMS_Share_MOD_Log_str
!
    SUBROUTINE GEMS_Share_MOD_Log_str(llvl, str)
    integer          :: llvl
    character(LEN=*) :: str

    integer          :: rv

    rv = gems_share_log_str(llvl, str)
#ifdef PRNLOG
    print*, trim(str)
#endif

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_str

!
! ---------- GEMS_Share_MOD_Log_str2
!
    SUBROUTINE GEMS_Share_MOD_Log_str2(llvl, str, str2)
    integer          :: llvl
    character(LEN=*) :: str
    character(LEN=*) :: str2

    integer          :: rv

    rv = gems_share_log_str2(llvl, str, str2)
#ifdef PRNLOG
    print*, trim(str)
#endif

    RETURN
    END SUBROUTINE GEMS_Share_MOD_Log_str2

END MODULE Share_MOD_Log

!-------------------------------------------------------------------------------
