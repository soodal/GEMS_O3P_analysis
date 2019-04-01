!-------------------------------------------------------------------------------
!+

MODULE O3P_MOD_Output
!-------------------------------------------------------------------------------
!+Description: 
!     .
! 
!+Version Date       Comment
! ------- ---------- -------
! 0.0     2016.05.12 Fisrt Code (R&D, SaeaSoft Co., Ltd.) 
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Modules:
 USE OMSAO_precision_module
 USE GEMS_O3P_gemsdata_module, ONLY: gems_ny
!**********************************************************!
 IMPLICIT NONE

!**********************************************************!


CONTAINS


SUBROUTINE GEMS_O3P_Write ( fitcol, dfitcol, exval)


USE Share_MOD_Env
USE Share_MOD_Constants_Variables
USE Share_l2_o3p_mod_write_read ! /share/src/mall/fsrc
USE O3P_MOD_L2wr
USE OMSAO_precision_module
USE GEMS_O3P_gemsdata_module, ONLY: nfxtrack, gems_time, gems_sza, gems_vza, gems_aza, &
                                     gems_lon, gems_lat,gems_clouds, nxbin, nybin, &
                                     first_line, last_line, first_pix, last_pix, offline , gems_ny, lineloc,pixloc

USE ozprof_data_module,      ONLY: nlay,nsfc,nalb,ntp, nlay_fit,num_iter,use_lograd, use_oe,&       
       atmosprof, ozprof_ap, ozprof_apstd,ozprof, ozprof_std,ozprof_nstd,strataod, stratsca, tropaod, tropsca, & 
       covar, contri, avg_kernel,favg_kernel,weight_function,ozdfs, ozinfo,  &  
       ozprof_start_index, ozprof_end_index, ozfit_start_index, ozfit_end_index, start_layer, end_layer,&
       which_alb,which_clima,which_cld,aerwavs, maxawin, actawin,eff_alb, eff_alb_init,  &      
       the_ai, the_cld_flg, the_ctp, the_cfrac, the_cod, lambcld_refl, do_lambcld, &       
       do_simu,ozwrtres, ozwrtavgk, ozwrtvar,ozwrtcorr, ozwrtcovar, ozwrtcontri, ozwrtfavgk, ozwrtsnr,ozwrtwf,radcalwrt,&
       fgasidxs, ngas, nfgas, gaswrt, tracegas,trace_prof, trace_avgk, trace_contri ,    &        
       saa_flag,nsaa_spike, the_snowice , the_landwater_flag, the_glint_flag, glintprob
USE OMSAO_variables_module, ONLY : currloop, currpix,&                    
                                   the_sza_atm, the_vza_atm, the_aza_atm, nloc,the_lons, the_lats,the_time,&
                                   n_fitvar_rad,mask_fitvar_rad, fitvar_rad_str, fitvar_rad,fitvar_rad_apriori, &
		                   fitvar_rad_std,fitvar_rad_nstd,  &
                                   n_rad_wvl, fitwavs,fitweights,fitspec_rad, fitres_rad, linenum_lim,   &
                                   simspec_rad, clmspec_rad, actspec_rad, numwin, nradpix,database, refidx  

IMPLICIT NONE

!----------
! IN/OUT variables
!----------
 INTEGER, INTENT(IN) :: exval
 REAL (KIND=dp), DIMENSION(3), INTENT (IN)    :: fitcol
 REAL (KIND=dp), DIMENSION(3, 2), INTENT (IN) :: dfitcol
!----------
! Local variables
!----------

 CHARACTER(LEN=128)  :: ctl_fpath
 INTEGER(KIND=4)     :: errCode = 0, the_x,the_y, &
                        i,j,the_i,the_j, fidx, lidx
 REAL (KIND=dp) :: avgres, rms
 REAL (KIND=dp), DIMENSION (numwin)       :: allrms, allavgres
 LOGICAL, SAVE :: first = .TRUE.

 
IF ( first) THEN 

    !------------------------------------------
    !---- 1. Getting Configuration ------------
    !------------------------------------------
    ctl_fpath = "../../../share/conf/gems.conf"
    CALL GEMS_Share_GetL2O3pEnv(ctl_fpath, errCode)
    IF ( errCode /= 0 ) GOTO 1999

    !------------------------------------------
    !---- 2. Initialize Global Constants ------
    !------------------------------------------
    CALL GEMS_Share_Init_GlobalConstants

    !------------------------------------------
    !---- 3. Initialize L2O3P Global Constants-
    !------------------------------------------
    CALL GEMS_Share_Init_L2O3P_GlobalConstants

    !------------------------------------------
    !---- 4. Memory Allocation of TypeVariable
    !------------------------------------------
    errCode = GEMS_Share_MOD_L2O3P_MemAlloc()
    IF ( errCode /= 0 ) GOTO 1999

    errCode = GEMS_Share_MOD_L2O3P_MemInit()
    IF ( errCode /= 0 ) GOTO 1999

    IF (.not.ASSOCIATED(gds_L2O3P%o3p_dfld%OzRet)) THEN
        errCode = -1
        WRITE(*, '( A53 )') "[gds_L2O3P%o3p_dfld%OzRet] Memory Allocation Failed !"
        
    ENDIF
    first = .FALSE.  
ENDIF

    the_x = pixloc(currpix) ; the_y = lineloc(currloop)
    
    !-------------------------------------------------------------------------
    !--------  계산된 결과를 L2 O3P 전역변수에 저장                   --------
    !-------------------------------------------------------------------------
  
     fidx = 1
     DO i = 1, numwin
        lidx = fidx + nradpix(i) - 1
        allrms(i) = SQRT(SUM((fitres_rad(fidx:lidx)/fitweights(fidx:lidx))**2) / nradpix(i))
        fidx = lidx + 1
     ENDDO
     rms = SQRT(SUM((ABS(fitres_rad(1:n_rad_wvl)) / &
     fitweights(1:n_rad_wvl))**2.0)/n_rad_wvl)

     simspec_rad(1:n_rad_wvl) = fitspec_rad(1:n_rad_wvl) - fitres_rad(1:n_rad_wvl)
     IF (use_lograd) THEN
        fitspec_rad(1:n_rad_wvl) = EXP(fitspec_rad(1:n_rad_wvl))
        simspec_rad(1:n_rad_wvl) = EXP(simspec_rad(1:n_rad_wvl))
        fitres_rad(1:n_rad_wvl) = fitspec_rad(1:n_rad_wvl) -  simspec_rad(1:n_rad_wvl)
     END IF
     avgres = SQRT(SUM((ABS(fitres_rad(1:n_rad_wvl)) / &
          fitspec_rad(1:n_rad_wvl))**2.0)/n_rad_wvl)*100.0   

     fidx = 1
     DO i = 1, numwin
        lidx = fidx + nradpix(i) - 1
        allavgres(i) = SQRT(SUM((ABS(fitres_rad(fidx:lidx)) / &
             fitspec_rad(fidx:lidx))**2.0)/nradpix(i))*100.0
        fidx = lidx + 1
     ENDDO

    !==================================      
    ! GEOlocation Data
    !================================== 
    gds_L2O3P%o3p_gfld%LAT(the_x,the_y)         =  gems_lat  (currpix, currloop)
    gds_L2O3P%o3p_gfld%LON(the_x,the_y)         =  gems_lon  (currpix, currloop)
    gds_L2O3P%o3p_gfld%SolZenAng(the_x,the_y)   =  gems_sza  (currpix, currloop)
    gds_L2O3P%o3p_gfld%ViewZenAng(the_x,the_y)  =  gems_vza  (currpix, currloop)
    gds_L2O3P%o3p_gfld%RelAziAng(the_x,the_y)   =  gems_aza  (currpix, currloop)
    gds_L2O3P%o3p_gfld%Time(the_y)              =  gems_time (currloop)
    gds_L2O3P%o3p_gfld%Pix_tmp(the_x,the_y)     =  the_x    ! geun
    gds_L2O3P%o3p_gfld%Line_tmp(the_x,the_y)    =  the_y    ! geun

    if (exval <= -9) RETURN

    gds_L2O3P%o3p_gfld%P(1:nlay+1,the_x,the_y)     =  atmosprof(1,0:nlay)
    gds_L2O3P%o3p_gfld%Alt(1:nlay+1,the_x,the_y)   =  atmosprof(2,0:nlay)
    gds_L2O3P%o3p_gfld%TEMP(1:nlay+1,the_x,the_y)  =  atmosprof(3,0:nlay) ! change dimension
    !gds_L2O3P%o3p_gfld%TerrP(the_x,the_y) =  atmosprof(1,nlay)
    !gds_L2O3P%o3p_gfld%ptrp(the_x,the_y)  =  atmosprof(1, ntp)
    !==================================      
    ! Retrieved Data
    !================================== 
    gds_L2O3P%o3p_dfld%AvgK(1:nlay,1:nlay,the_x,the_y) =  avg_kernel(ozfit_start_index:ozfit_end_index,ozfit_start_index:ozfit_end_index)
    gds_L2O3P%o3p_dfld%OzRet(1:nlay,the_x,the_y)   = ozprof(1:nlay)
    gds_L2O3P%o3p_dfld%OzNErr(1:nlay,the_x,the_y)  = ozprof_nstd(1:nlay)
    gds_L2O3P%o3p_dfld%OzSMErr(1:nlay,the_x,the_y) = ozprof_std (1:nlay)
    gds_L2O3P%o3p_dfld%OzAP(1:nlay,the_x,the_y)    = ozprof_ap(1:nlay)
    gds_L2O3P%o3p_dfld%OzApErr(1:nlay,the_x,the_y) = ozprof_apstd(1:nlay)
    gds_L2O3P%o3p_dfld%ColAmtOz(:,the_x,the_y)     = fitcol(1:3)
    gds_L2O3P%o3p_dfld%DFS(1,the_x,the_y)          = ozdfs

    gds_L2O3P%o3p_dfld%CldP(the_x,the_y)         =  the_ctp
    gds_L2O3P%o3p_dfld%EffCldFrac(the_x,the_y) = the_cfrac
    gds_L2O3P%o3p_dfld%SfcAlb(the_x,the_y)     = eff_alb(2)

    gds_L2O3P%o3p_dfld%PrsQFlag(the_x,the_y)       = exval
    gds_L2O3P%o3p_dfld%NumIter(the_x,the_y)        = num_iter





    rms    = SQRT(SUM((fitres_rad(1:n_rad_wvl)/fitweights(1:n_rad_wvl))**2) / n_rad_wvl)
    simspec_rad(1:n_rad_wvl) = fitspec_rad(1:n_rad_wvl) - fitres_rad(1:n_rad_wvl)
    IF (use_lograd) THEN
        fitspec_rad(1:n_rad_wvl) = EXP(fitspec_rad(1:n_rad_wvl))
        simspec_rad(1:n_rad_wvl) = EXP(simspec_rad(1:n_rad_wvl))
        fitres_rad(1:n_rad_wvl) = fitspec_rad(1:n_rad_wvl) -  simspec_rad(1:n_rad_wvl)
    END IF
    avgres = SQRT(SUM((ABS(fitres_rad(1:n_rad_wvl)) / &
             fitspec_rad(1:n_rad_wvl))**2.0)/n_rad_wvl)*100.0   

    gds_L2O3P%o3p_dfld%AvgRes(1,the_x,the_y) = avgres
    gds_L2O3P%o3p_dfld%RMS(the_x,the_y)      = rms

 

    
  RETURN 



1999  CONTINUE


2999  CONTINUE

    print*, "ERROR !"

    !------------------------------------------
    !---- 6. Memory DeAllocation of TypeVariable
    !------------------------------------------    

    errCode = GEMS_Share_MOD_L2O3P_MemDeAlloc()

    RETURN 


END SUBROUTINE GEMS_O3P_Write

SUBROUTINE GEMS_O3P_SUB5_Output(L2O3P_WR, retCode)
  USE Share_l2_o3p_mod_write_read ! define L2_O3P
  USE O3P_MOD_L2wr 

    IMPLICIT NONE

    TYPE(L2_o3p),    INTENT(INOUT)  :: L2O3P_WR
    INTEGER(KIND=4), INTENT(OUT)    :: retCode         ! return error code

    print*, "+++ 5. OUTPUT"
    print*, "Here is GEMS_O3P_SUB5_Output() subroutine. !!!"
    print*, " "
    print*, " "

    CALL GEMS_O3P_L2wr(L2O3P_wr, retCode)
    !write(*,*) L2O3P_wr%o3p_gfld%lat(2,:)
    IF ( retCode /= 0 ) THEN  
      WRITE(*, '( A12 I5)') "result code=", retCode
      RETURN
    ENDIF

    RETURN

END SUBROUTINE GEMS_O3P_SUB5_Output


END MODULE O3P_MOD_Output
