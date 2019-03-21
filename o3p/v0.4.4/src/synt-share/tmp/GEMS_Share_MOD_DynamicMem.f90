!-------------------------------------------------------------------------------
!+Module to process the dynamic memory

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE Share_MOD_DynamicMem

!-------------------------------------------------------------------------------
!+Description: 
!     Allocate or Deallocate the dynamic memory.
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2015.01.27 Fisrt Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Modules:

USE Share_MOD_Constants_Variables, ONLY: sz_kind,           &
                                         gvf_angles,        &
                                         gvs_size,          &
                                         gvf_latlon,        &
                                         gvf_euv1_rad,      &
                                         gvf_euv1_radwave,  &
                                         gvf_euv1_prec,     &
                                         gvf_euv2_rad,      &
                                         gvf_euv2_radwave,  &
                                         gvf_euv2_prec,     &
                                         gvf_evis_rad,      &
                                         gvf_evis_radwave,  &
                                         gvf_suv2_irrad,    &
                                         gvf_suv2_prec,     &
                                         gvf_suv2_irrwave,  &
                                         gvf_svis_irrad,    &
                                         gvf_svis_irrwave,  &
                                         gvf_albedo,        &
                                         gvf_elevation,     &
                                         gvf_amfgeo,        &
                                         gvi_flg_ls,        &
                                         gvi_flg_cm,        &
                                         gvi_flg_px
                                                        

IMPLICIT NONE
SAVE

CONTAINS

!-------------------------------------------------------------------------------
!+Description: 
!      전역변수 메모리를 할당해줌
!
!
! Method:
!
!
! Input files:
!
! Output files:
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2015.01.27 Fisrt Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------
SUBROUTINE GEMS_Share_AllocMem

    CALL GEMS_Share_AllocMem1st

END SUBROUTINE GEMS_Share_AllocMem

SUBROUTINE GEMS_Share_AllocMem2

    CALL GEMS_Share_AllocMem1st
    CALL GEMS_Share_AllocMem2nd

END SUBROUTINE GEMS_Share_AllocMem2

SUBROUTINE GEMS_Share_AllocMem1st

    ALLOCATE( gvf_angles%vza( gvs_size%angles%nx,      gvs_size%angles%ny) )
    ALLOCATE( gvf_angles%sza( gvs_size%angles%nx,      gvs_size%angles%ny) )
    ALLOCATE( gvf_angles%raa( gvs_size%angles%nx,      gvs_size%angles%ny) )

    ALLOCATE( gvf_latlon%lat( gvs_size%latlon%nx,      gvs_size%latlon%ny) )
    ALLOCATE( gvf_latlon%lon( gvs_size%latlon%nx,      gvs_size%latlon%ny) )

    ALLOCATE( gvf_albedo(     gvs_size%albedo%nx,      gvs_size%albedo%ny,  gvs_size%albedo%nz) )

    ALLOCATE( gvf_elevation(  gvs_size%elevation%nx,   gvs_size%elevation%ny) )
    ALLOCATE( gvf_amfgeo(     gvs_size%amfgeo%nx,      gvs_size%amfgeo%ny) )

    ALLOCATE( gvi_flg_ls(     gvs_size%flg_ls%nx,      gvs_size%flg_ls%ny) )
    ALLOCATE( gvi_flg_cm(     gvs_size%flg_cm%nx,      gvs_size%flg_cm%ny) )
    ALLOCATE( gvi_flg_px(     gvs_size%flg_px%nx,      gvs_size%flg_px%ny) )

END SUBROUTINE GEMS_Share_AllocMem1st

SUBROUTINE GEMS_Share_AllocMem2nd

    CALL GEMS_Share_AllocMem2nd_4L1BRUG_EUV1
    CALL GEMS_Share_AllocMem2nd_4L1BRUG_EUV2
    CALL GEMS_Share_AllocMem2nd_4LIBRVG     
    CALL GEMS_Share_AllocMem2nd_4LIBIRR_SUV2
    CALL GEMS_Share_AllocMem2nd_4LIBIRR_SVIS

END SUBROUTINE GEMS_Share_AllocMem2nd

SUBROUTINE GEMS_Share_AllocMem2nd_4L1BRUG_EUV1

    ALLOCATE( gvf_euv1_rad    ( gvs_size%euv1_rad%nx,    gvs_size%euv1_rad%ny,    gvs_size%euv1_rad%nz    ) )
    ALLOCATE( gvf_euv1_radwave( gvs_size%euv1_radwav%nx, gvs_size%euv1_radwav%ny, gvs_size%euv1_radwav%nz ) )
    ALLOCATE( gvf_euv1_prec   ( gvs_size%euv1_prec%nx,   gvs_size%euv1_prec%ny,   gvs_size%euv1_prec%nz   ) )

END SUBROUTINE GEMS_Share_AllocMem2nd_4L1BRUG_EUV1

SUBROUTINE GEMS_Share_AllocMem2nd_4L1BRUG_EUV2

    ALLOCATE( gvf_euv2_rad    ( gvs_size%euv2_rad%nx,    gvs_size%euv2_rad%ny,    gvs_size%euv2_rad%nz    ) )
    ALLOCATE( gvf_euv2_radwave( gvs_size%euv2_radwav%nx, gvs_size%euv2_radwav%ny, gvs_size%euv2_radwav%nz ) )
    ALLOCATE( gvf_euv2_prec   ( gvs_size%euv2_prec%nx,   gvs_size%euv2_prec%ny,   gvs_size%euv2_prec%nz   ) )

END SUBROUTINE GEMS_Share_AllocMem2nd_4L1BRUG_EUV2

SUBROUTINE GEMS_Share_AllocMem2nd_4LIBRVG

    ALLOCATE( gvf_evis_rad    ( gvs_size%evis_rad%nx,    gvs_size%evis_rad%ny,    gvs_size%evis_rad%nz    ) )
    ALLOCATE( gvf_evis_radwave( gvs_size%evis_radwav%nx, gvs_size%evis_radwav%ny, gvs_size%evis_radwav%nz ) )

END SUBROUTINE GEMS_Share_AllocMem2nd_4LIBRVG

SUBROUTINE GEMS_Share_AllocMem2nd_4LIBIRR_SUV2

    ALLOCATE( gvf_suv2_irrad  ( gvs_size%suv2_irrad%nx,  gvs_size%suv2_irrad%ny  ) )
    ALLOCATE( gvf_suv2_prec   ( gvs_size%suv2_prec%nx,   gvs_size%suv2_prec%ny   ) )
    ALLOCATE( gvf_suv2_irrwave( gvs_size%suv2_irrwav%nx, gvs_size%suv2_irrwav%ny ) )

END SUBROUTINE GEMS_Share_AllocMem2nd_4LIBIRR_SUV2

SUBROUTINE GEMS_Share_AllocMem2nd_4LIBIRR_SVIS

    ALLOCATE( gvf_svis_irrad  ( gvs_size%svis_irrad%nx,  gvs_size%svis_irrad%ny  ) )
    ALLOCATE( gvf_svis_irrwave( gvs_size%svis_irrwav%nx, gvs_size%svis_irrwav%ny ) )

END SUBROUTINE GEMS_Share_AllocMem2nd_4LIBIRR_SVIS

!-------------------------------------------------------------------------------
!+Description: 
!      전역변수 메모리를 해제함
!
!
! Method:
!
!
! Input files:
!
! Output files:
!
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.1     2015.01.27 Fisrt Code (YuGeun Ki, Seasoft) 
!-------------------------------------------------------------------------------
SUBROUTINE GEMS_Share_DeallocMem

    CALL GEMS_Share_DeallocMem1st

END SUBROUTINE GEMS_Share_DeallocMem

SUBROUTINE GEMS_Share_DeallocMem2

    CALL GEMS_Share_DeallocMem1st
    CALL GEMS_Share_DeallocMem2nd

END SUBROUTINE GEMS_Share_DeallocMem2

SUBROUTINE GEMS_Share_DeallocMem1st

    DEALLOCATE( gvf_angles%vza )
    DEALLOCATE( gvf_angles%sza )
    DEALLOCATE( gvf_angles%raa )

    DEALLOCATE( gvf_latlon%lat )
    DEALLOCATE( gvf_latlon%lon )

    DEALLOCATE( gvf_albedo     )

    DEALLOCATE( gvf_elevation  )
    DEALLOCATE( gvf_amfgeo     )

    DEALLOCATE( gvi_flg_ls     )
    DEALLOCATE( gvi_flg_cm     )
    DEALLOCATE( gvi_flg_px     )

END SUBROUTINE GEMS_Share_DeallocMem1st

SUBROUTINE GEMS_Share_DeallocMem2nd

    CALL GEMS_Share_DeallocMem2nd_4L1BRUG_EUV1
    CALL GEMS_Share_DeallocMem2nd_4L1BRUG_EUV2
    CALL GEMS_Share_DeallocMem2nd_4LIBRVG     
    CALL GEMS_Share_DeallocMem2nd_4LIBIRR_SUV2
    CALL GEMS_Share_DeallocMem2nd_4LIBIRR_SVIS

END SUBROUTINE GEMS_Share_DeallocMem2nd

SUBROUTINE GEMS_Share_DeallocMem2nd_4L1BRUG_EUV1

    IF(ASSOCIATED(gvf_euv1_rad    )) DEALLOCATE( gvf_euv1_rad    )
    IF(ASSOCIATED(gvf_euv1_radwave)) DEALLOCATE( gvf_euv1_radwave)
    IF(ASSOCIATED(gvf_euv1_prec   )) DEALLOCATE( gvf_euv1_prec    )

    NULLIFY( gvf_euv1_rad    )
    NULLIFY( gvf_euv1_radwave)
    NULLIFY( gvf_euv1_prec    )

    RETURN

END SUBROUTINE GEMS_Share_DeallocMem2nd_4L1BRUG_EUV1

SUBROUTINE GEMS_Share_DeallocMem2nd_4L1BRUG_EUV2

    IF(ASSOCIATED(gvf_euv2_rad    )) DEALLOCATE( gvf_euv2_rad      )
    IF(ASSOCIATED(gvf_euv2_radwave)) DEALLOCATE( gvf_euv2_radwave  )
    IF(ASSOCIATED(gvf_euv2_prec   )) DEALLOCATE( gvf_euv2_prec    )

    NULLIFY( gvf_euv2_rad      ) 
    NULLIFY( gvf_euv2_radwave  ) 
    NULLIFY( gvf_euv2_prec    )

    RETURN
END SUBROUTINE GEMS_Share_DeallocMem2nd_4L1BRUG_EUV2

SUBROUTINE GEMS_Share_DeallocMem2nd_4LIBRVG

    IF(ASSOCIATED(gvf_evis_rad    )) DEALLOCATE( gvf_evis_rad     )
    IF(ASSOCIATED(gvf_evis_radwave)) DEALLOCATE( gvf_evis_radwave )

    NULLIFY( gvf_evis_rad     )
    NULLIFY( gvf_evis_radwave )

    RETURN
END SUBROUTINE GEMS_Share_DeallocMem2nd_4LIBRVG

SUBROUTINE GEMS_Share_DeallocMem2nd_4LIBIRR_SUV2

    IF(ASSOCIATED(gvf_suv2_irrad  )) DEALLOCATE( gvf_suv2_irrad   )
    IF(ASSOCIATED(gvf_suv2_prec   )) DEALLOCATE( gvf_suv2_prec    )
    IF(ASSOCIATED(gvf_suv2_irrwave)) DEALLOCATE( gvf_suv2_irrwave )

    NULLIFY( gvf_suv2_irrad   )
    NULLIFY( gvf_suv2_prec    )
    NULLIFY( gvf_suv2_irrwave )

    RETURN
END SUBROUTINE GEMS_Share_DeallocMem2nd_4LIBIRR_SUV2

SUBROUTINE GEMS_Share_DeallocMem2nd_4LIBIRR_SVIS

    IF(ASSOCIATED(gvf_svis_irrad  )) DEALLOCATE( gvf_svis_irrad   )
    IF(ASSOCIATED(gvf_svis_irrwave)) DEALLOCATE( gvf_svis_irrwave )

    NULLIFY( gvf_svis_irrad   ) 
    NULLIFY( gvf_svis_irrwave ) 

    RETURN
END SUBROUTINE GEMS_Share_DeallocMem2nd_4LIBIRR_SVIS

END MODULE Share_MOD_DynamicMem

!-------------------------------------------------------------------------------
