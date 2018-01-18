!--------------------------------------------------------------
!+Description :
!     Main Subroutine for preparing GEMS measurement and delivering results to specfit_ozprof
!     
! Method:
!     Algorithm description file - 13 page
!
! Input files:
!     
!
! Output files:
!    
!
!+Version  Date        Comment
! -------- ----------- -------
!  0.2.0   2015.07.01  Modified Code for GEMS (Bak Juseon, PNU UNIV.)
!--------------------------------------------------------------


SUBROUTINE GEMS_O3P_SUB4_fitting_process ( pge_error_status)


  USE OMSAO_precision_module
  USE OMSAO_variables_module,  ONLY: currpix,currloop,currtrack, currline,currtime, &
                                     the_lons, the_lats,ozabs_convl, so2crs_convl, wavcal, scnwrt, &
                                     npix_fitting, npix_fitted,n_fitvar_rad, fitvar_rad_saved,mask_fitvar_rad 
  USE ozprof_data_module, ONLY:num_iter
 !USE OMSAO_slitfunction_module
  USE GEMS_o3P_geo_module
  USE OMSAO_errstat_module
  USE GEMS_O3P_gemsdata_module, ONLY: ncoadd, nxbin, nybin, lineloc, pixloc, &
                                      first_pix, last_pix, first_line,last_line, offline, &
                                      gems_rad, gems_irrad, gems_exitval, gems_fitvar, gems_initval
  USE Share_l2_o3p_mod_write_read
  USE O3P_MOD_Output, ONLY:GEMS_O3P_Write
  IMPLICIT NONE

  ! -----------------------
  ! Input/Output variables
  ! ----------------------
  INTEGER, INTENT (OUT) :: pge_error_status

  ! -------------------------
  ! Local variables 
  ! -------------------------
  INTEGER :: initval, errstat, exval
  INTEGER :: curr_fitted_line, nxcoadd, iy, ix
  REAL (KIND=dp), DIMENSION(3)    :: fitcol
  REAL (KIND=dp), DIMENSION(3, 2) :: dfitcol
  REAL (KIND=dp)     :: fitcol_avg, rms_avg, dfitcol_avg, drel_fitcol_avg, rms
  LOGICAL            :: reduce_resolution_save

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=23), PARAMETER :: modulename = 'GEMS_O3P_fitting_process'

  pge_error_status = pge_errstat_ok
 
  nxcoadd    = ncoadd*nxbin
  
  !------------------------------
  ! Read GeoLocation data

  !------------------------------
  WRITE(*, '(A)') ' => Preparing Geolocation data'
  CALL gems_o3p_prep_geo (last_line, offline, pge_error_status)
  IF ( pge_error_status >= pge_errstat_error ) RETURN
  

  !------------------------------
  ! Read Cloud Pressure
  !------------------------------
  WRITE(*, '(A)') ' => Preparing Cloud data'
  CALL gems_o3p_prep_cld (last_line, offline, pge_error_status)
  IF ( pge_error_status >= pge_errstat_error ) RETURN


  !---------------------------
  ! Read radiance
  !------------------------------
  WRITE(*,'(A, i4,"-",3i4)') ' => Preparing Radiance data'
  CALL GEMS_O3P_read_l1b_rad ( nxcoadd, first_pix, last_pix,  last_line, offline, pge_error_status)
  If (pge_error_status >= pge_errstat_error) RETURN
  

  
  IF (wavcal) THEN
     !CALL gems_o3p_rad_cross_calibrate (first_pix, last_pix, last_line, offline, pge_error_status)
     IF ( pge_error_status >= pge_errstat_error ) RETURN
     IF (scnwrt) WRITE(*, '(A,I4,A,I4)') '  (6.4) Finishing calibrating radiances for lines: ',first_line+offline, offline + nybin*last_line
  ENDIf


  ! Initialize fitting statistics
  npix_fitting = 0       ! number of pixels (failure + success)
  npix_fitted  = 0       ! number of successfully fitted pixels   
  fitcol_avg   = 0.0;  rms_avg = 0.0; dfitcol_avg = 0.0; drel_fitcol_avg = 0.0
  curr_fitted_line =0
  gems_exitval(:, :) = -10
  gems_initval(:,:)  = 0

  WRITE(*, '(6A5,2A8,A6, a5,A8)') 'LINE','XPix','CLine','Cpix','Lloc','Ploc', 'Lon','Lat','Exval','RMS','NUM'

  GEMS_PIX: DO currpix = first_pix, last_pix   
     
     ! Load/adjust irradiances and slit calibration parameters
     CALL gems_o3p_adj_solar_data (pge_error_status)
     IF ( pge_error_status >= pge_errstat_error ) CYCLE
     ozabs_convl =.TRUE. ; so2crs_convl=.TRUE.


     GEMS_LINE : DO currloop = first_line, last_line

      currtrack = (currpix-1) * nxbin  + 1
      currline  = (currloop-1) * nybin + offline + 1
            
      IF (gems_rad%errstat(currpix, currloop) == pge_errstat_error .OR. &
          gems_irrad%errstat(currpix) == pge_errstat_error ) &
          gems_exitval(currpix, currloop) = -9
    
      ! Load/adjust radiances/geolocations fields for a particular pixel
      ! Prepare databases for the first pixel (ifitline == 1)   
      IF (gems_exitval(currpix, currloop) == -10) THEN
        ! print * , curr_fitted_line, currloop
         CALL gems_o3p_adj_earthshine_data (curr_fitted_line,last_line, pge_error_status)
	  curr_fitted_line = curr_fitted_line + 1
         IF ( pge_error_status >= pge_errstat_error ) gems_exitval(currpix, currloop) = -9
      ENDIF


     CALL timestamp(currtime)
      IF (gems_exitval(currpix, currloop) == -10) THEN
              
         initval = gems_initval(currpix, currloop)
         CALL specfit_ozprof (initval, fitcol, dfitcol, rms, exval)
         gems_exitval(currpix, currloop) = exval  ! Sotre exit status for current pixel
         gems_fitvar(currpix, currloop, 1:n_fitvar_rad) = fitvar_rad_saved(mask_fitvar_rad(1:n_fitvar_rad))
     ELSE
         exval = -9
     ENDIF
     
     WRITE(*, '(6I5,2f8.2,I6, f8.3,i5, A27)') currline, currtrack,currloop, currpix,lineloc(currloop), pixloc(currpix), &
                                              the_lons(5), the_lats(5) ,exval, rms, num_iter, currtime

     ! Write retrieval
      CALL gems_o3p_write (fitcol, dfitcol,exval) ! if( exval > -9) 
    
     IF ( exval >= 0 .AND. fitcol(1) > 0.0 .AND. dfitcol(1, 1) >= 0.0 ) THEN                     
       ! ----------------------------------------------------------------------
       ! Some general statistics on the average fitted column and uncertainty.
       ! Again, we make sure that only "good" fits are included in the average.
       ! ----------------------------------------------------------------------              
       fitcol_avg       = fitcol_avg + fitcol(1)
       rms_avg          = rms_avg + rms
       dfitcol_avg      = dfitcol_avg + dfitcol(1, 1)
       drel_fitcol_avg  = drel_fitcol_avg + dfitcol(1, 1) / fitcol(1)
       npix_fitted      = npix_fitted + 1
       curr_fitted_line = curr_fitted_line + 1
     ENDIF
     npix_fitting = npix_fitting + 1   
  ENDDO GEMS_LINE
  ENDDO GEMS_PIX
  
  IF ( npix_fitted == 0) npix_fitted = 1
  !IF (scnwrt) CALL write_final(fitcol_avg, rms_avg, dfitcol_avg, drel_fitcol_avg, npix_fitted)
  IF (scnwrt) WRITE(*, '(2(A,I5))') 'Number of pixels = ', &
       npix_fitting, '   Number of fitted pixels = ', npix_fitted

   
  RETURN
END SUBROUTINE GEMS_O3P_SUB4_fitting_process    




