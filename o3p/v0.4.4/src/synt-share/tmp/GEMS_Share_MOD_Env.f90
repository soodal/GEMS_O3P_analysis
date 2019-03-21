!-------------------------------------------------------------------------------
!+Module to handle Environment Variables

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE Share_MOD_Env

!-------------------------------------------------------------------------------
!+Description: 
!     Module to handle Environment Variables
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2015.02.24 Fisrt Code (YuGeun Ki, SAEAsoft) 
!-------------------------------------------------------------------------------
      
! Declarations:

! Modules used:

!**********************************************************!
IMPLICIT NONE
    CHARACTER(LEN=256)  :: gvc_in_l1b_fpath
    CHARACTER(LEN=256)  :: gvc_in_l1b_fname
    CHARACTER(LEN=256)  :: gvc_nml_lv2_fpath
    CHARACTER(LEN=256)  :: gvc_nml_lv2_fname
    CHARACTER(LEN=256)  :: gvc_in_lv2_fpath
    CHARACTER(LEN=256)  :: gvc_in_lv2_fname
    CHARACTER(LEN=256)  :: gvc_out_lv2_fpath
    CHARACTER(LEN=256)  :: gvc_out_lv2_fname

    TYPE :: Lvl2Env
        CHARACTER(LEN=256)  :: nml_lv2_fpath
        CHARACTER(LEN=256)  :: nml_lv2_fname
        CHARACTER(LEN=256)  :: in_lv2_fpath
        CHARACTER(LEN=256)  :: in_lv2_fname
        CHARACTER(LEN=256)  :: out_lv2_fpath
        CHARACTER(LEN=256)  :: out_lv2_fname
        INTEGER(KIND=1)     :: setUp = 0
    END TYPE Lvl2Env

    TYPE(Lvl2Env)           :: gds_CldL2Env
    TYPE(Lvl2Env)           :: gds_HchoL2Env
    TYPE(Lvl2Env)           :: gds_AodL2Env
    TYPE(Lvl2Env)           :: gds_O3tL2Env
    TYPE(Lvl2Env)           :: gds_O3pL2Env
    TYPE(Lvl2Env)           :: gds_No2L2Env
    TYPE(Lvl2Env)           :: gds_So2L2Env
    TYPE(Lvl2Env)           :: gds_AlbdL2Env
    TYPE(Lvl2Env)           :: gds_AlhL2Env

!**********************************************************!
CONTAINS


!-------------------------------------------------------------------------------
!+Description: 
!      subroutine to get L1B and L2 file path infomation from namelist.
!
! Method:
!
! input  :
!       env_fpath = the environment file path
!
! output :
!
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1.0   2015.03.06 Fisrt Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------
SUBROUTINE GEMS_Share_GetEnv(env_fpath)
    CHARACTER(LEN=*),   INTENT(IN)  :: env_fpath

    Namelist /LEVEL1B_File_Info/gvc_in_l1b_fpath,   &
                                gvc_in_l1b_fname
    Namelist /LEVEL2_File_Info/gvc_nml_lv2_fpath,   &
                               gvc_nml_lv2_fname,   & 
                               gvc_in_lv2_fpath,    &
                               gvc_in_lv2_fname,    &
                               gvc_out_lv2_fpath,   &
                               gvc_out_lv2_fname

    !--------------------------------
    !---     Read Control file    ---
    !--------------------------------
    OPEN(10, file=trim(env_fpath), status='old')
    READ(10, LEVEL1B_File_Info)
    READ(10, LEVEL2_File_Info)
    CLOSE(10)

    RETURN
END SUBROUTINE GEMS_Share_GetEnv

!-------------------------------------------------------------------------------
!+Description: 
!      subroutine to get L1B file path infomation from namelist.
!
! Method:
!
! input  :
!       env_fpath = the environment file path
!
! output :
!       rtCode    = the result code of reading L1B environment
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1.0   2015.03.06 Fisrt Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------
SUBROUTINE GEMS_Share_GetL1BEnv(env_fpath, rtCode)
    CHARACTER(LEN=*),   INTENT(IN)  :: env_fpath
    INTEGER(KIND=4),    INTENT(OUT) :: rtCode

    Namelist /LEVEL1B_File_Info/gvc_in_l1b_fpath,   &
                                gvc_in_l1b_fname

    !--------------------------------
    !---     Read Control file    ---
    !--------------------------------
    OPEN(10, file=trim(env_fpath), status='old')
    READ(10, LEVEL1B_File_Info, iostat=rtCode)
    CLOSE(10)

    RETURN
END SUBROUTINE GEMS_Share_GetL1BEnv

!-------------------------------------------------------------------------------
!+Description: 
!      subroutine to get L2 file path infomation of L2 CLD from namelist.
!
! Method:
!
! input  :
!       env_fpath = the environment file path
!

!       rtCode    = the result code of reading L2 CLD environment
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1.0   2015.03.06 Fisrt Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------
SUBROUTINE GEMS_Share_GetL2CldEnv(env_fpath, rtCode)
    CHARACTER(LEN=*),   INTENT(IN)  :: env_fpath
    INTEGER(KIND=4),    INTENT(OUT) :: rtCode

    Namelist /LEVEL2_File_Info_Cld/gvc_nml_lv2_fpath,   &
                                   gvc_nml_lv2_fname,   & 
                                   gvc_in_lv2_fpath,    &
                                   gvc_in_lv2_fname,    &
                                   gvc_out_lv2_fpath,   &
                                   gvc_out_lv2_fname

    !--------------------------------
    !---     Read Control file    ---
    !--------------------------------
    OPEN(10, file=trim(env_fpath), status='old')
    READ(10, LEVEL2_File_Info_Cld, iostat=rtCode)

    IF ( rtCode < 0 ) THEN 
        CLOSE(10)
        RETURN
    ENDIF

    gds_CldL2Env%nml_lv2_fpath  = gvc_nml_lv2_fpath
    gds_CldL2Env%nml_lv2_fname  = gvc_nml_lv2_fname
    gds_CldL2Env%in_lv2_fpath   = gvc_in_lv2_fpath  
    gds_CldL2Env%in_lv2_fname   = gvc_in_lv2_fname  
    gds_CldL2Env%out_lv2_fpath  = gvc_out_lv2_fpath
    gds_CldL2Env%out_lv2_fname  = gvc_out_lv2_fname  
    gds_CldL2Env%setUp          = 1

    CLOSE(10)

    rtCode = 0
    RETURN
END SUBROUTINE GEMS_Share_GetL2CldEnv

!-------------------------------------------------------------------------------
!+Description: 
!      subroutine to get L2 file path infomation of L2 HCHO from namelist.
!
! Method:
!
! input  :
!       env_fpath = the environment file path
!
! output :
!       rtCode    = the result code of reading L2 HCHO environment
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1.0   2015.03.06 Fisrt Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------
SUBROUTINE GEMS_Share_GetL2HchoEnv(env_fpath, rtCode)
    CHARACTER(LEN=*),   INTENT(IN)  :: env_fpath
    INTEGER(KIND=4),    INTENT(OUT) :: rtCode

    Namelist /LEVEL2_File_Info_Hcho/gvc_nml_lv2_fpath,   &
                                    gvc_nml_lv2_fname,   & 
                                    gvc_in_lv2_fpath,    &
                                    gvc_in_lv2_fname,    &
                                    gvc_out_lv2_fpath,   &
                                    gvc_out_lv2_fname

    !--------------------------------
    !---     Read Control file    ---
    !--------------------------------
    OPEN(10, file=trim(env_fpath), status='old')
    READ(10, LEVEL2_File_Info_Hcho, iostat=rtCode)
    IF ( rtCode < 0 ) THEN 
        CLOSE(10)
        RETURN
    ENDIF

    gds_HchoL2Env%nml_lv2_fpath  = gvc_nml_lv2_fpath
    gds_HchoL2Env%nml_lv2_fname  = gvc_nml_lv2_fname
    gds_HchoL2Env%in_lv2_fpath   = gvc_in_lv2_fpath  
    gds_HchoL2Env%in_lv2_fname   = gvc_in_lv2_fname  
    gds_HchoL2Env%out_lv2_fpath  = gvc_out_lv2_fpath
    gds_HchoL2Env%out_lv2_fname  = gvc_out_lv2_fname  
    gds_HchoL2Env%setUp          = 1

    CLOSE(10)

    rtCode = 0
    RETURN
END SUBROUTINE GEMS_Share_GetL2HchoEnv

!-------------------------------------------------------------------------------
!+Description: 
!      subroutine to get L2 file path infomation of L2 AOD from namelist.
!
! Method:
!
! input  :
!       env_fpath = the environment file path
!
! output :
!       rtCode    = the result code of reading L2 AOD environment
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1.0   2015.04.03 Fisrt Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------
SUBROUTINE GEMS_Share_GetL2AodEnv(env_fpath, rtCode)
    CHARACTER(LEN=*),   INTENT(IN)  :: env_fpath
    INTEGER(KIND=4),    INTENT(OUT) :: rtCode

    Namelist /LEVEL2_File_Info_Aod/gvc_nml_lv2_fpath,   &
                                    gvc_nml_lv2_fname,   & 
                                    gvc_in_lv2_fpath,    &
                                    gvc_in_lv2_fname,    &
                                    gvc_out_lv2_fpath,   &
                                    gvc_out_lv2_fname

    !--------------------------------
    !---     Read Control file    ---
    !--------------------------------
    OPEN(10, file=trim(env_fpath), status='old')
    READ(10, LEVEL2_File_Info_Aod, iostat=rtCode)
    IF ( rtCode < 0 ) THEN 
        CLOSE(10)
        RETURN
    ENDIF

    gds_AodL2Env%nml_lv2_fpath  = gvc_nml_lv2_fpath
    gds_AodL2Env%nml_lv2_fname  = gvc_nml_lv2_fname
    gds_AodL2Env%in_lv2_fpath   = gvc_in_lv2_fpath  
    gds_AodL2Env%in_lv2_fname   = gvc_in_lv2_fname  
    gds_AodL2Env%out_lv2_fpath  = gvc_out_lv2_fpath
    gds_AodL2Env%out_lv2_fname  = gvc_out_lv2_fname  
    gds_AodL2Env%setUp          = 1

    CLOSE(10)

    rtCode = 0
    RETURN
END SUBROUTINE GEMS_Share_GetL2AodEnv


!-------------------------------------------------------------------------------
!+Description: 
!      subroutine to get L2 file path infomation of L2 O3T from namelist.
!
! Method:
!
! input  :
!       env_fpath = the environment file path
!
! output :
!       rtCode    = the result code of reading L2 O3T environment
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1.0   2015.04.03 Fisrt Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------
SUBROUTINE GEMS_Share_GetL2O3tEnv(env_fpath, rtCode)
    CHARACTER(LEN=*),   INTENT(IN)  :: env_fpath
    INTEGER(KIND=4),    INTENT(OUT) :: rtCode

    Namelist /LEVEL2_File_Info_O3t/gvc_nml_lv2_fpath,   &
                                   gvc_nml_lv2_fname,   & 
                                   gvc_in_lv2_fpath,    &
                                   gvc_in_lv2_fname,    &
                                   gvc_out_lv2_fpath,   &
                                   gvc_out_lv2_fname

    !--------------------------------
    !---     Read Control file    ---
    !--------------------------------
    OPEN(10, file=trim(env_fpath), status='old')
    READ(10, LEVEL2_File_Info_O3t, iostat=rtCode)
    IF ( rtCode < 0 ) THEN 
        CLOSE(10)
        RETURN
    ENDIF

    gds_O3tL2Env%nml_lv2_fpath  = gvc_nml_lv2_fpath
    gds_O3tL2Env%nml_lv2_fname  = gvc_nml_lv2_fname
    gds_O3tL2Env%in_lv2_fpath   = gvc_in_lv2_fpath  
    gds_O3tL2Env%in_lv2_fname   = gvc_in_lv2_fname  
    gds_O3tL2Env%out_lv2_fpath  = gvc_out_lv2_fpath
    gds_O3tL2Env%out_lv2_fname  = gvc_out_lv2_fname  
    gds_O3tL2Env%setUp          = 1

    CLOSE(10)

    rtCode = 0
    RETURN
END SUBROUTINE GEMS_Share_GetL2O3tEnv


!-------------------------------------------------------------------------------
!+Description: 
!      subroutine to get L2 file path infomation of L2 O3P from namelist.
!
! Method:
!
! input  :
!       env_fpath = the environment file path
!
! output :
!       rtCode    = the result code of reading L2 O3P environment
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1.0   2015.03.06 Fisrt Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------
SUBROUTINE GEMS_Share_GetL2O3pEnv(env_fpath, rtCode)
    CHARACTER(LEN=*),   INTENT(IN)  :: env_fpath
    INTEGER(KIND=4),    INTENT(OUT) :: rtCode

    Namelist /LEVEL2_File_Info_O3p/gvc_nml_lv2_fpath,   &
                                   gvc_nml_lv2_fname,   & 
                                   gvc_in_lv2_fpath,    &
                                   gvc_in_lv2_fname,    &
                                   gvc_out_lv2_fpath,   &
                                   gvc_out_lv2_fname

    !--------------------------------
    !---     Read Control file    ---
    !--------------------------------
    OPEN(10, file=trim(env_fpath), status='old')
    READ(10, LEVEL2_File_Info_O3p, iostat=rtCode)

    IF ( rtCode < 0 ) THEN 
        CLOSE(10)
        RETURN
    ENDIF

    gds_O3pL2Env%nml_lv2_fpath  = gvc_nml_lv2_fpath
    gds_O3pL2Env%nml_lv2_fname  = gvc_nml_lv2_fname
    gds_O3pL2Env%in_lv2_fpath   = gvc_in_lv2_fpath  
    gds_O3pL2Env%in_lv2_fname   = gvc_in_lv2_fname  
    gds_O3pL2Env%out_lv2_fpath  = gvc_out_lv2_fpath
    gds_O3pL2Env%out_lv2_fname  = gvc_out_lv2_fname  
    gds_O3pL2Env%setUp          = 1

    CLOSE(10)

    rtCode = 0
    RETURN
END SUBROUTINE GEMS_Share_GetL2O3pEnv

!-------------------------------------------------------------------------------
!+Description: 
!      subroutine to get L2 file path infomation of L2 NO2 from namelist.
!
! Method:
!
! input  :
!       env_fpath = the environment file path
!
! output :
!       rtCode    = the result code of reading L2 NO2 environment
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1.0   2015.03.06 Fisrt Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------
SUBROUTINE GEMS_Share_GetL2No2Env(env_fpath, rtCode)
    CHARACTER(LEN=*),   INTENT(IN)  :: env_fpath
    INTEGER(KIND=4),    INTENT(OUT) :: rtCode

    Namelist /LEVEL2_File_Info_No2/gvc_nml_lv2_fpath,   &
                                   gvc_nml_lv2_fname,   & 
                                   gvc_in_lv2_fpath,    &
                                   gvc_in_lv2_fname,    &
                                   gvc_out_lv2_fpath,   &
                                   gvc_out_lv2_fname

    !--------------------------------
    !---     Read Control file    ---
    !--------------------------------
    OPEN(10, file=trim(env_fpath), status='old')
    READ(10, LEVEL2_File_Info_No2, iostat=rtCode)
    IF ( rtCode < 0 ) THEN 
        CLOSE(10)
        RETURN
    ENDIF

    gds_No2L2Env%nml_lv2_fpath  = gvc_nml_lv2_fpath
    gds_No2L2Env%nml_lv2_fname  = gvc_nml_lv2_fname
    gds_No2L2Env%in_lv2_fpath   = gvc_in_lv2_fpath  
    gds_No2L2Env%in_lv2_fname   = gvc_in_lv2_fname  
    gds_No2L2Env%out_lv2_fpath  = gvc_out_lv2_fpath
    gds_No2L2Env%out_lv2_fname  = gvc_out_lv2_fname  
    gds_No2L2Env%setUp          = 1

    CLOSE(10)

    rtCode = 0
    RETURN
END SUBROUTINE GEMS_Share_GetL2No2Env

!-------------------------------------------------------------------------------
!+Description: 
!      subroutine to get L2 file path infomation of L2 SO2 from namelist.
!
! Method:
!
! input  :
!       env_fpath = the environment file path
!
! output :
!       rtCode    = the result code of reading L2 SO2 environment
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1.0   2015.04.03 Fisrt Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------
SUBROUTINE GEMS_Share_GetL2So2Env(env_fpath, rtCode)
    CHARACTER(LEN=*),   INTENT(IN)  :: env_fpath
    INTEGER(KIND=4),    INTENT(OUT) :: rtCode

    Namelist /LEVEL2_File_Info_So2/gvc_nml_lv2_fpath,   &
                                   gvc_nml_lv2_fname,   & 
                                   gvc_in_lv2_fpath,    &
                                   gvc_in_lv2_fname,    &
                                   gvc_out_lv2_fpath,   &
                                   gvc_out_lv2_fname

    !--------------------------------
    !---     Read Control file    ---
    !--------------------------------
    OPEN(10, file=trim(env_fpath), status='old')
    READ(10, LEVEL2_File_Info_So2, iostat=rtCode)
    IF ( rtCode < 0 ) THEN 
        CLOSE(10)
        RETURN
    ENDIF

    gds_So2L2Env%nml_lv2_fpath  = gvc_nml_lv2_fpath
    gds_So2L2Env%nml_lv2_fname  = gvc_nml_lv2_fname
    gds_So2L2Env%in_lv2_fpath   = gvc_in_lv2_fpath  
    gds_So2L2Env%in_lv2_fname   = gvc_in_lv2_fname  
    gds_So2L2Env%out_lv2_fpath  = gvc_out_lv2_fpath
    gds_So2L2Env%out_lv2_fname  = gvc_out_lv2_fname  
    gds_So2L2Env%setUp          = 1

    CLOSE(10)

    rtCode = 0
    RETURN
END SUBROUTINE GEMS_Share_GetL2So2Env


!-------------------------------------------------------------------------------
!+Description: 
!      subroutine to get L2 file path infomation of L2 ALBD from namelist.
!
! Method:
!
! input  :
!       env_fpath = the environment file path
!
! output :
!       rtCode    = the result code of reading L2 ALBD environment
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1.0   2015.04.03 Fisrt Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------
SUBROUTINE GEMS_Share_GetL2AlbdEnv(env_fpath, rtCode)
    CHARACTER(LEN=*),   INTENT(IN)  :: env_fpath
    INTEGER(KIND=4),    INTENT(OUT) :: rtCode

    Namelist /LEVEL2_File_Info_Albd/gvc_nml_lv2_fpath,   &
                                    gvc_nml_lv2_fname,   & 
                                    gvc_in_lv2_fpath,    &
                                    gvc_in_lv2_fname,    &
                                    gvc_out_lv2_fpath,   &
                                    gvc_out_lv2_fname

    !--------------------------------
    !---     Read Control file    ---
    !--------------------------------
    OPEN(10, file=trim(env_fpath), status='old')
    READ(10, LEVEL2_File_Info_Albd, iostat=rtCode)
    IF ( rtCode < 0 ) THEN 
        CLOSE(10)
        RETURN
    ENDIF

    gds_AlbdL2Env%nml_lv2_fpath  = gvc_nml_lv2_fpath
    gds_AlbdL2Env%nml_lv2_fname  = gvc_nml_lv2_fname
    gds_AlbdL2Env%in_lv2_fpath   = gvc_in_lv2_fpath  
    gds_AlbdL2Env%in_lv2_fname   = gvc_in_lv2_fname  
    gds_AlbdL2Env%out_lv2_fpath  = gvc_out_lv2_fpath
    gds_AlbdL2Env%out_lv2_fname  = gvc_out_lv2_fname  
    gds_AlbdL2Env%setUp          = 1

    CLOSE(10)

    rtCode = 0
    RETURN
END SUBROUTINE GEMS_Share_GetL2AlbdEnv

!-------------------------------------------------------------------------------
!+Description: 
!      subroutine to get L2 file path infomation of L2 Alh from namelist.
!
! Method:
!
! input  :
!       env_fpath = the environment file path
!
! output :
!       rtCode    = the result code of reading L2 Alh environment
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1.0   2016.12.12 Fisrt Code (YuGeun Ki) 
!-------------------------------------------------------------------------------
SUBROUTINE GEMS_Share_GetL2AlhEnv(env_fpath, rtCode)
    CHARACTER(LEN=*),   INTENT(IN)  :: env_fpath
    INTEGER(KIND=4),    INTENT(OUT) :: rtCode

    Namelist /LEVEL2_File_Info_Alh/gvc_nml_lv2_fpath,   &
                                    gvc_nml_lv2_fname,   & 
                                    gvc_in_lv2_fpath,    &
                                    gvc_in_lv2_fname,    &
                                    gvc_out_lv2_fpath,   &
                                    gvc_out_lv2_fname

    !--------------------------------
    !---     Read Control file    ---
    !--------------------------------
    OPEN(10, file=trim(env_fpath), status='old')
    READ(10, LEVEL2_File_Info_Alh, iostat=rtCode)
    IF ( rtCode < 0 ) THEN 
        CLOSE(10)
        RETURN
    ENDIF

    gds_AlhL2Env%nml_lv2_fpath  = gvc_nml_lv2_fpath
    gds_AlhL2Env%nml_lv2_fname  = gvc_nml_lv2_fname
    gds_AlhL2Env%in_lv2_fpath   = gvc_in_lv2_fpath  
    gds_AlhL2Env%in_lv2_fname   = gvc_in_lv2_fname  
    gds_AlhL2Env%out_lv2_fpath  = gvc_out_lv2_fpath
    gds_AlhL2Env%out_lv2_fname  = gvc_out_lv2_fname  
    gds_AlhL2Env%setUp          = 1

    CLOSE(10)

    rtCode = 0
    RETURN
END SUBROUTINE GEMS_Share_GetL2AlhEnv

END MODULE

!-------------------------------------------------------------------------------
