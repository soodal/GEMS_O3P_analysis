!NOTE MOre check for coading process

MODULE GEMS_o3P_geo_module

  ! ============================================================================== !
  !
  ! 1. subset variables for all cross-track and defined along-track 
  !    NOTE : mpi calculation divide  along-track pixels amongh the three
  !    Example : gems_lat (1:10) = gems_geo%lat(11:20)
  ! 2. interpret binary flags
  ! 3. apply spatial co-adding
  !                                                   Juseon Bak @PNU
  ! ============================================================================== !

  USE OMSAO_precision_module
  USE OMSAO_errstat_module

  IMPLICIT NONE



  ! ----------------
  ! Local Parameters
  ! ---------------------------------------------------------------------
  ! * Values for Pi (rad, deg) and Conversions between Degree and Radians
  ! ---------------------------------------------------------------------
  REAL (KIND=r8), PARAMETER, PRIVATE :: pi         = 3.14159265358979_r8  ! 2*ASIN(1.0_r8)
  REAL (KIND=r8), PARAMETER, PRIVATE :: pihalf     = 0.5_r8  * pi
  REAL (KIND=r8), PARAMETER, PRIVATE :: twopi      = 2.0_r8  * pi
  REAL (KIND=r8), PARAMETER, PRIVATE :: pi_deg     = 180.0_r8
  REAL (KIND=r8), PARAMETER, PRIVATE :: pihalf_deg =  90.0_r8
  REAL (KIND=r8), PARAMETER, PRIVATE :: twopi_deg  = 360.0_r8
  REAL (KIND=r8), PARAMETER, PRIVATE :: deg2rad    = pi / 180.0_r8
  REAL (KIND=r8), PARAMETER, PRIVATE :: rad2deg    = 180.0_r8 / pi
  ! ---------------------------------------------------------------------
  ! * Precison for DEG <-> RAD conversion - anything less than EPS
  !   is effectively ZERO.
  ! ---------------------------------------------------------------------
  REAL (KIND=r8), PARAMETER, PRIVATE :: eps = 1.0E-10_r8
  ! ---------------------------------------------------------------------
  REAL (KIND=r4), PARAMETER, PRIVATE :: rearth0 = 6378  ! equatorial radius
  REAL (KIND=r4), PARAMETER, PRIVATE :: minza = 0.0, maxza=90.0, minaza = -360., maxaza = 360.0


CONTAINS

  SUBROUTINE GEMS_O3P_prep_geo ( ny, offline, pge_error_status )

  USE OMSAO_variables_module,   ONLY: band_selectors, numwin
  USE GEMS_O3P_gemsdata_module, ONLY: nchannel, nfxtrack,  ntimes, nxbin, nybin, &
                                      gemsl1btype, gems_uv1, gems_uv2, & 
                                      ! output
                                      gems_lat, gems_lon,  gems_elat, gems_elon, gems_clat, gems_clon, & 
                                      gems_sza, gems_vza,  gems_aza, gems_sca, &
                                      gems_time, gems_height, &
                                      gems_XTrackQFlg, gems_xflag, & 
                                      gems_groundQflg, gems_gflag
 
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER , INTENT (IN) :: ny, offline
    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i4), INTENT (OUT) :: pge_error_status

    ! ---------------
    ! Local variables
    ! ---------------
    TYPE (gemsl1btype) :: gems_geo    
    INTEGER :: i,j, ix,iix, iy, iiy,sline, eline, nline,nx, nbx,error, nbits, nt
    INTEGER :: ysidx, yeidx, ymidx, xsidx, xeidx, xmidx
    REAL (KIND=r8), DIMENSION (nfxtrack, ntimes)  :: lat, lon, sza, vza, saz, vaz
    INTEGER (KINd=i2), DIMENSION(nfxtrack, ntimes) :: hgt

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
     CHARACTER (LEN=23), PARAMETER :: modulename = 'GEMS_O3P_read_geofield'

    pge_error_status = pge_errstat_ok


   ! -----------------------------------------------------------
   !  Open data block (UV-1, if both are selected) 
   !  NOTE: will be removed for GEMS
   ! -----------------------------------------------------------

   IF (ANY(band_selectors(1:numwin) >= 1)) THEN 
      gems_geo = gems_uv1
   ELSE 
      gems_geo = gems_uv2
   ENDIF

   ! Initialize output varialbes

    gems_lon  (:,:) = 0.0  ;  gems_lat (:,:)  = 0.0
    gems_clon (:,:) = 0.0 ;  gems_clat (:,:) = 0.0
    gems_elon (:,:) = 0.0 ;  gems_elat (:,:) = 0.0
    gems_sza (:,:) = 0.0 ;   gems_vza (:,:)  = 0.0 
    gems_aza (:,:) = 0.0 ;   gems_sca(:,:)   = 0.0
    gems_time (:)  = 0.0 ;   gems_height(:,:)  = 0.0
    
   !===============================================
   ! angle & geolocation
   !===============================================

    sline = offline +1 ; eline  = offline + ny * nybin 
    IF (eline == sline) eline = sline + 1
    nline = eline - sline + 1
    nx = nfxtrack

    lat(1:nx, 1:nline) = gems_geo%lat(1:nx, sline:eline)*1.0_r8
    lon(1:nx, 1:nline) = gems_geo%lon(1:nx, sline:eline)*1.0_r8
    sza(1:nx, 1:nline) = gems_geo%sza(1:nx, sline:eline)*1.0_r8
    vza(1:nx, 1:nline) = gems_geo%vza(1:nx, sline:eline)*1.0_r8
    saz(1:nx, 1:nline) = gems_geo%saz(1:nx, sline:eline)*1.0_r8
    vaz(1:nx, 1:nline) = gems_geo%vaz(1:nx, sline:eline)*1.0_r8  

    CALL get_sphgeoview_corners (nx, nline, nxbin, nybin , &
         lon(1:nx, 1:nline), lat(1:nx, 1:nline),         &
         sza(1:nx, 1:nline), saz(1:nx, 1:nline),      &
         vza(1:nx, 1:nline), vaz(1:nx, 1:nline),      &
         gems_clon(0:nx, 1:nline+1), gems_clat(0:nx, 1:nline+1), &
         gems_elon(0:nx, 1:nline), gems_elat(0:nx, 1:nline),   &
         gems_sza (1:nx, 1:nline), gems_vza(1:nx, 1:nline),    &
         gems_aza (1:nx, 1:nline), gems_sca(1:nx, 1:nline))
  
   gems_lon (1:nx, 1:nline) = lon(1:nx, 1:nline)
   gems_lat (1:nx, 1:nline) = lat(1:nx, 1:nline)

   !================================================
   ! time, height, flags
   !================================================  
   sline = offline +1 ; eline  = offline + ny * nybin 
   nline = eline - sline + 1
   nx    = nfxtrack

   gems_time(1:nline)             = gems_geo%time(sline:eline) ! not averaged, just selected
   gems_height(1:nx, 1:nline)     = gems_geo%height(1:nx, sline:eline)
   gems_groundQflg(1:nx, 1:nline) = gems_geo%groundQflg(1:nx, sline:eline)

   !--------------------------------------
   ! For XtrackQuality flags, always use that from UV2 because UV1 and UV2 are different
   ! and more pixels are filtered in UV2
   !------------------------------------
   nx = gems_uv2%nxtrack
   gems_xtrackqflg (1:nx, 1:nline) = gems_uv2%xtrackqflg(1:nx, sline:eline)
  
   ! Coadd the flags
   IF (nchannel == 2) THEN
     ! correct the bug on 2016 08 10 
     i = 0 
     nbits = 8
     DO ix = 1, nfxtrack
       i = ix * 2 - 1
       CALL coadd_byte_qflgs(nbits, nline, gems_Xtrackqflg(i, 1:nline), gems_Xtrackqflg(i + 1, 1:nline)) 
       gems_XTrackQFlg(ix, 1:nline) = gems_Xtrackqflg(i, 1:nline)
    ENDDO
   ENDIF

   WHERE (gems_XTrackQFlg(1:nfxtrack, 1:nline) == -127)
        gems_XTrackQFlg(1:nfxtrack, 1:nline) = 0
   ENDWHERE


   IF (nybin > 1) THEN 
      hgt (1:nfxtrack, 1:nline) = gems_height (1:nfxtrack,1:nline)*1.0_r8
      gems_height (1:nfxtrack,1:nline) = 0
   ENDIF

   DO iy = 1, ny

      IF ( nybin > 1 ) THEN 
        ysidx = (iy - 1) * nybin +1
        yeidx = ysidx + nybin - 1
        ymidx = ysidx + nybin / 2

        gems_time(iy) = gems_time(ymidx)
        gems_XTrackQFlg(1:nfxtrack, iy)= gems_XTrackQFlg(1:nfxtrack, ymidx)
        gems_groundQflg(1:nfxtrack,iy) = gems_groundQflg(1:nfxtrack,ymidx) 
        
        DO iiy = 1, nybin 
             gems_height(1:nfxtrack,iy) = gems_height(1:nfxtrack,iy)  + hgt(1:nfxtrack,ysidx+iiy-1) 
        ENDDO
             gems_height(1:nfxtrack, iy) = NINT (1.0 * gems_height(1:nfxtrack, iy) / nybin)   
      ENDIF

     CALL convert_gpqualflag_info (nfxtrack, gems_groundQflg(1:nfxtrack, iy), &
          gems_gflag%land_water(1:nfxtrack,iy), &
          gems_gflag%glint(1:nfxtrack,iy), &
          gems_gflag%snow_ice(1:nfxtrack,iy))
  ENDDO

  ! Resample the water/land, sea glint, and snow/oce flags if xbin, surface altitude
  IF (nxbin > 1) THEN 
     DO ix = 1, nfxtrack / nxbin
        xsidx = (ix - 1) * nxbin + 1
        xeidx = xsidx + nxbin - 1
        xmidx = xsidx + nxbin / 2
        iix = (ix - 1) * nxbin + NINT(nxbin / 2.0)  ! 1 3 5 7
        gems_gflag%land_water(ix, 1:ny) = gems_gflag%land_water(iix, 1:ny)
        gems_gflag%snow_ice(ix, 1:ny)   = gems_gflag%snow_ice(iix, 1:ny)
	gems_XTrackQFlg(ix, 1:ny)       = gems_XTrackQFlg(xmidx, 1:ny)
        iix = (ix - 1) * nxbin ! 0, 2, 4
       DO i = 1, ny 
        gems_height(ix, i) = NINT(SUM(1.0 * gems_height(iix+1:iix+nxbin, i)) / nxbin)
        gems_gflag%glint(ix, i)  = NINT(SUM(1.0 * gems_gflag%glint(iix+1:iix+nxbin, i))  / nxbin)
       ENDDO
     ENDDO
  ENDIF
  ! Derive Row Anomaly Related Flags
  nx = nfxtrack/nxbin
  DO i = 1, ny 
    CALL convert_xtrackqfag_info ( nx, gems_xtrackqflg(1:nx, i), gems_xflag%rowanomaly(1:nx, i), &
    gems_xflag%waveshift(1:nx, i), gems_xflag%blockage(1:nx, i),gems_xflag%straysun(1:nx, i), gems_xflag%strayearth(1:nx, i) )    
  ENDDO

  RETURN
  END SUBROUTINE GEMS_O3P_prep_geo


  SUBROUTINE get_sphgeoview_corners (nxtrack, ntimes, nxbin, nybin, lon, lat, sza, saza, vza, vaza, &
         clon, clat, elon, elat, esza, evza, eaza, esca)

    IMPLICIT NONE
    
    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER,                                               INTENT(IN)    :: nxtrack, ntimes, nxbin, nybin
    REAL    (KIND=r8), DIMENSION (1:nxtrack, 0:ntimes-1), INTENT(INOUT) :: lon, lat, sza, saza, vza, vaza
    REAL    (KIND=r8), DIMENSION (0:nxtrack, 0:ntimes),   INTENT(OUT)   :: clon, clat
    REAL    (KIND=r8), DIMENSION (0:nxtrack, 0:ntimes-1), INTENT(OUT)   :: elon, elat
    REAL    (KIND=r8), DIMENSION (1:nxtrack, 0:ntimes-1), INTENT(OUT)   :: esza, evza, eaza, esca     

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                                   :: i, j, jj, ix, mpix, nx, ny
    INTEGER                                             :: errstat
    !REAL    (KIND=r8)                                  :: lat1, lat2, phi1, phi2, dis
    !REAL    (KIND=r8)                                  :: a0, b0, c0, gam0, alp0, a, gam
    REAL    (KIND=r8), DIMENSION (1:nxtrack,0:ntimes-1) :: omixsize
    REAL    (KIND=r8), DIMENSION (0:nxtrack, 0:ntimes-1):: edsza, edsazm, edvza, edvazm
    REAL    (KIND=r8), DIMENSION (1:nxtrack)            :: tmpdisx, xsize, tmpxmid
    REAL    (KIND=r8), DIMENSION (0:nxtrack)            :: tmpx, tmpsza, tmpvza, tmpsaza, tmpvaza
    REAL    (KIND=r8), DIMENSION (0:ntimes-1)           :: tmpdisy
    REAL    (KIND=dp), DIMENSION (3)                    :: zen0, zen, sazm, vazm, relaza

    ! ------------------------------------------------------------------
    ! Convert geolocation to radians; do everything in R8 rather than R4
    ! ------------------------------------------------------------------
    lon = lon * deg2rad; lat = lat * deg2rad

    ! -------------------------
    ! Initialize some variables
    ! -------------------------
    clon   = -999.9 ; clat = -999.9
    elon   = -999.9 ; elat = -999.9
    esza   = -999.9 ; evza = -999.9; eaza = -999.9; esca = -999.9

    ! Perform interpolation across the track
    DO i = 0, ntimes - 1

       ! Compute the distances between two pixels: (x1 + x2) / 2.
       DO ix = 1, nxtrack - 1 
          tmpdisx(ix) = circle_rdis(lat(ix, i), lon(ix, i), lat(ix+1, i), lon(ix+1, i))
       ENDDO

       ! Compute the pixel size across the track
       ! Assume the center two pixels have equal pixel size (which causes about < 0.1 km error for UV-2)
       mpix = nxtrack / 2
       xsize(mpix) = tmpdisx(mpix) / 2.0;  xsize(mpix + 1) = xsize(mpix)
       
       DO ix = mpix -1, 1, -1
          xsize(ix) = tmpdisx(ix) - xsize(ix + 1)
       ENDDO
       DO ix = mpix + 2, nxtrack
          xsize(ix) = tmpdisx(ix - 1) - xsize(ix - 1)
       ENDDO
       omixsize(:, i) = xsize

       !!  This is to test SUBROUTINE sphergeom_intermediate
       !!  Works for both interpolation and extrapolation (with certain limitation) 
       !lat1 = -88.0 * deg2rad; lat2 = -88.0 * deg2rad
       !phi1 = 0.0  * deg2rad; phi2 = 180.0 * deg2rad 
       !dis  = circle_rdis(lat1, phi1, lat2, phi2)     
       !print *, dis
       !CALL sphergeom_intermediate ( lat1, phi1, lat2, phi2, dis, dis*0.8, a, gam )
       !WRITE(*, '(6D14.6)')  lat1, phi1, a, gam, lat2, phi2
       !print *, circle_rdis(lat1, phi1, a, gam)/dis
       !print *, circle_rdis(a, gam, lat2, phi2)/dis
  
       ! Perform interpolation 
       DO ix = 1, nxtrack - 1          
          CALL sphergeom_intermediate(lat(ix, i), lon(ix, i), lat(ix+1, i), lon(ix+1, i), &
               tmpdisx(ix), xsize(ix), elat(ix, i), elon(ix, i))
       ENDDO
       ix = 1
       CALL sphergeom_intermediate(lat(ix, i), lon(ix, i), lat(ix+1, i), lon(ix+1, i),   &
            tmpdisx(ix), -xsize(ix), elat(ix-1, i), elon(ix-1, i))
       ix = nxtrack
       CALL sphergeom_intermediate(lat(ix, i), lon(ix, i), lat(ix-1, i), lon(ix-1, i),   &
            tmpdisx(ix-1), -xsize(ix), elat(ix, i), elon(ix, i))   

       ! Compute viewing geometry for west and east edge
       ! Performal interpolation/extrapolation (2 points) along the spherical lines (good enough)
       tmpx = 0.0
       DO ix = 1, nxtrack
          tmpx(ix) = tmpx(ix-1) + xsize(ix)
       ENDDO
       tmpxmid(1:nxtrack) = (tmpx(0:nxtrack-1) + tmpx(1:nxtrack)) / 2.0
       CALL interpol(tmpxmid(1:nxtrack), sza(1:nxtrack, i),  nxtrack, tmpx(0:nxtrack), &
            tmpsza(0:nxtrack),  nxtrack+1, errstat)      
       CALL interpol(tmpxmid(1:nxtrack), saza(1:nxtrack, i), nxtrack, tmpx(0:nxtrack), &
            tmpsaza(0:nxtrack), nxtrack+1, errstat)
       CALL interpol(tmpxmid(1:nxtrack), vza(1:nxtrack, i),  nxtrack, tmpx(0:nxtrack), &
            tmpvza(0:nxtrack),  nxtrack+1, errstat)
       CALL interpol(tmpxmid(1:nxtrack), vaza(1:nxtrack, i), nxtrack, tmpx(0:nxtrack), &
            tmpvaza(0:nxtrack), nxtrack+1, errstat)

       ! Check center pixel
       DO ix = mpix-1, mpix + 1
          IF (tmpvaza(ix) < 0) THEN
             tmpvaza(ix) = -tmpvaza(ix-1); EXIT
          ENDIF
       ENDDO

       edsza(0:nxtrack, i) = tmpsza(0:nxtrack); edsazm(0:nxtrack, i) = tmpsaza(0:nxtrack)
       edvza(0:nxtrack, i) = tmpvza(0:nxtrack); edvazm(0:nxtrack, i) = tmpvaza(0:nxtrack)
    ENDDO

    !PRINT *
    !WRITE(*, *) 'Solar Angle'
    !WRITE(*, '(10F9.2)') esza(1:nxtrack, 0)
    !WRITE(*, *) 'View Angle'
    !WRITE(*, '(10F9.2)') evza(1:nxtrack, 0)
    !WRITE(*, *) 'Relative Azmimuthal Angle'
    !WRITE(*, '(10F9.2)') eaza(1:nxtrack, 0)
    !WRITE(*, *) 'Scattering Angle'
    !WRITE(*, '(10F9.2)') esca(1:nxtrack, 0)
 
    ! Perform interpolation along the track with a simipler but simpler (center) approach
    ! since the pixel size along the track does not vary much
    IF (ntimes == 1) tmpdisy(0) = 0.00212031  ! ~ 13.5 km
    DO ix = 0, nxtrack
       DO i = 0, ntimes - 2 
          tmpdisy(i) = circle_rdis(elat(ix, i), elon(ix, i), elat(ix, i+1), elon(ix, i+1))
          CALL sphergeom_intermediate(elat(ix, i), elon(ix, i), elat(ix, i+1), &
               elon(ix, i+1), tmpdisy(i), tmpdisy(i)*0.5, clat(ix, i+1), clon(ix, i+1))  
          !print *, i+1, clat(ix, i+1), clon(ix, i+1)
       ENDDO

       i = 0
       CALL sphergeom_intermediate(elat(ix, i), elon(ix, i), elat(ix, i+1),    &
            elon(ix, i+1), tmpdisy(i), -tmpdisy(i)*0.5, clat(ix, i), clon(ix, i))
       !print *, i, clat(ix, i), clon(ix, i)

       i = ntimes - 1       
       CALL sphergeom_intermediate(elat(ix, i), elon(ix, i), elat(ix, i-1),    &
            elon(ix, i-1), tmpdisy(i-1), -tmpdisy(i-1)*0.5, clat(ix, i+1), clon(ix, i+1))
       !print *, i+1, clat(ix, i+1), clon(ix, i+1)
    ENDDO

    ! Perform coadding
    IF (nxbin > 1 .OR. nybin > 1) THEN
       nx = nxtrack / nxbin    
       ny = ntimes  / nybin

       ! cornor coordinates (only need sampling)
       j = 0
       DO ix = 0, nxtrack, nxbin
          clon(j, :) = clon(ix, :)
          clat(j, :) = clat(ix, :)
          j = j + 1
       ENDDO

       j = 0
       DO i = 0, ntimes, nybin
          clon(0:nx, j) = clon(0:nx, i)
          clat(0:nx, j) = clat(0:nx, i)
          j = j + 1
       ENDDO

       ! edge coordinates (easy to be re-computed from cornor coordinates)
       DO ix = 0, nx
          DO i = 0, ny - 1 
             tmpdisy(i) = circle_rdis(clat(ix, i), clon(ix, i), clat(ix, i+1), clon(ix, i+1))
             CALL sphergeom_intermediate(clat(ix, i), clon(ix, i), clat(ix, i+1), clon(ix, i+1), &
                  tmpdisy(i), tmpdisy(i)*0.5, elat(ix, i), elon(ix, i))  
             !print *, tmpdisy(i), elat(ix, i)*rad2deg, elon(ix, i) * rad2deg
             !print *, clat(ix, i)*rad2deg, clon(ix, i)*rad2deg, clat(ix, i+1)*rad2deg, clon(ix, i+1)*rad2deg
          ENDDO
       ENDDO

       ! Center coordinates (computed from edge coordinates)
       DO ix = 1, nx
          DO i = 0, ny - 1 
             tmpdisx(ix) = circle_rdis(elat(ix-1, i), elon(ix-1, i), elat(ix, i), elon(ix, i))
             CALL sphergeom_intermediate(elat(ix-1, i), elon(ix-1, i), elat(ix, i), elon(ix, i), &
                  tmpdisx(ix), tmpdisx(ix)*0.5, lat(ix, i), lon(ix, i))
             !print *, elat(ix-1, i)*rad2deg, elon(ix-1, i)*rad2deg, elat(ix, i)*rad2deg, elon(ix, i)*rad2deg
             !print *, lat(ix, i) *rad2deg, lon(ix, i)*rad2deg
          ENDDO
       ENDDO 

       ! Average edge viewing geometries along the track, sample along the track
       i = 0
       DO ix = 0, nxtrack, nxbin
          jj = 0
          DO j = 0, ntimes-1, nybin
             edsza (i, jj) = SUM(edsza (ix, j:j+nybin-1)) / nybin
             edsazm(i, jj) = SUM(edsazm(ix, j:j+nybin-1)) / nybin
             edvza (i, jj) = SUM(edvza (ix, j:j+nybin-1)) / nybin
             edvazm(i, jj) = SUM(edvazm(ix, j:j+nybin-1)) / nybin
             jj = jj + 1
          ENDDO
          i = i + 1
       ENDDO   

       ! Compute center viewing geometries (interpolate across the track)
       DO i = 0, ny-1
          tmpx = 0.0
          DO ix = 1, nx
             tmpx(ix) = tmpx(ix-1) + circle_rdis(elat(ix-1, i), elon(ix-1, i), elat(ix, i), elon(ix, i))
          ENDDO
          tmpxmid(1:nx) = (tmpx(0:nx-1) + tmpx(1:nx)) / 2.0

          CALL interpol(tmpx(0:nx), edsza (0:nx, i),  nx+1, tmpxmid(1:nx), sza (1:nx, i),  nx, errstat) 
          CALL interpol(tmpx(0:nx), edsazm(0:nx, i),  nx+1, tmpxmid(1:nx), saza(1:nx, i),  nx, errstat)  
     
          CALL interpol(tmpx(0:nx), edvza (0:nx, i),  nx+1, tmpxmid(1:nx), vza (1:nx, i),  nx, errstat) 
          CALL interpol(tmpx(0:nx), edvazm(0:nx, i),  nx+1, tmpxmid(1:nx), vaza(1:nx, i),  nx, errstat)      

          ! Check center pixel
          mpix = nx / 2
          DO ix = mpix - 1, mpix + 1
             IF (vaza(ix, i) < 0) THEN
                vaza(ix, i) = -vaza(ix-1, i); EXIT
             ENDIF
          ENDDO  
       ENDDO
    ELSE
       nx = nxtrack;       ny = ntimes
    ENDIF
    
    ! Now compute effective viewing geometry
    ! Compute effective viewing geometry for each pixel (at a certain atmosphere) 
    DO i = 0, ny-1
       DO ix = 1, nx
          IF ( sza(ix, i)  >= minza  .AND. sza(ix, i)  < maxza  .AND. &
               vza(ix, i)  >= minza  .AND. vza(ix, i)  < maxza  .AND. &
               saza(ix, i) >= minaza .AND. saza(ix, i) < maxaza .AND. &
               vaza(ix, i) >= minaza .AND. vaza(ix, i) < maxaza) THEN
             zen0(1) = edsza(ix-1, i) ; zen0(2) = sza(ix, i) ; zen0(3) = edsza(ix, i)
             zen(1)  = edvza(ix-1, i) ; zen(2)  = vza(ix, i) ; zen(3)  = edvza(ix, i)
             sazm(1) = edsazm(ix-1, i); sazm(2) = saza(ix, i); sazm(3) = edsazm(ix, i)
             vazm(1) = edvazm(ix-1, i); vazm(2) = vaza(ix, i); vazm(3) = edvazm(ix, i)
             relaza  = ABS(vazm - sazm)
             
             WHERE (vazm > 0) 
                zen = - zen
             ENDWHERE
             
             ! -180 < relaza < 180          
             WHERE (relaza > 180.0)
                relaza = 360.0 - relaza
             ENDWHERE
             
             CALL omi_angle_sat2toa (3, zen0, zen, relaza, esza(ix, i), evza(ix, i), eaza(ix, i), esca(ix, i) )
             !print *, zen0 * rad2deg
             !print *, zen  * rad2deg
             !print *, 180.0 - relaza * rad2deg
             !print *, esza(ix, i), evza(ix, i), eaza(ix, i), esca(ix, i)
             !STOP             
          ENDIF
       ENDDO
    ENDDO
    
    clon(0:nx, 0:ny)   = clon(0:nx, 0:ny)   * rad2deg    
    clat(0:nx, 0:ny)   = clat(0:nx, 0:ny)   * rad2deg
    lon(1:nx, 0:ny-1)  = lon(1:nx, 0:ny-1)  * rad2deg    
    lat(1:nx, 0:ny-1)  = lat(1:nx, 0:ny-1)  * rad2deg
    elon(0:nx, 0:ny-1) = elon(0:nx, 0:ny-1) * rad2deg
    elat(0:nx, 0:ny-1) = elat(0:nx, 0:ny-1) * rad2deg

    WHERE (clon(0:nx, 0:ny) > 180.0)
       clon(0:nx, 0:ny) = clon(0:nx, 0:ny) - 360.0
    ENDWHERE
    WHERE (clon(0:nx, 0:ny) < -180.0)
       clon(0:nx, 0:ny) = clon(0:nx, 0:ny) + 360.0
    ENDWHERE

    WHERE (lon(1:nx, 0:ny-1) > 180.0)
       lon(1:nx, 0:ny-1) = lon(1:nx, 0:ny-1) - 360.0
    ENDWHERE
    WHERE (lon(1:nx, 0:ny-1) < -180.0)
       lon(1:nx, 0:ny-1) = lon(1:nx, 0:ny-1) + 360.0
    ENDWHERE

    WHERE (elon(0:nx, 0:ny-1) > 180.0)
       elon(0:nx, 0:ny-1) = elon(0:nx, 0:ny-1) - 360.0
    ENDWHERE
    WHERE (elon(0:nx, 0:ny-1) < -180.0)
       elon(0:nx, 0:ny-1) = elon(0:nx, 0:ny-1) + 360.0
    ENDWHERE

    RETURN
  
  END SUBROUTINE get_sphgeoview_corners

  ! Compute spherical distance between two points
  ! Haversine Formula (from R.W. Sinnott, "virtue of the Haversine", 
  ! Sky and Telescope V68 (2), 1984, p159  
  ! lat1, lon1, lat2, lon2 in radians
  FUNCTION circle_dis(lat1, lon1, lat2, lon2) RESULT(dis)

    IMPLICIT NONE

    ! ----------------------
    ! Input/output variables
    ! -----------------------
    REAL (KIND=r8), INTENT (IN) :: lat1, lon1, lat2, lon2
    REAL (KIND=r8)              :: dis

    ! Local variable
    REAL (KIND=r8)              :: dlon, dlat, a, rdis, mlat
     
    dlat = lat2 - lat1;    dlon = lon2 - lon1; mlat = (lat1 + lat2) / 2.0

    a = MIN(1.0, SQRT( SIN(dlat/2.0)**2.0 + COS(lat1) * COS(lat2) * SIN(dlon/2.0)**2.0 )  )
    rdis = 2.0 * ASIN(a)                                    ! relative distaince in radiances
    dis = rdis * (rearth0 - 21.0 * SIN(mlat))               ! in the same unit of rearth0

    RETURN

  END FUNCTION circle_dis


  FUNCTION circle_rdis(lat1, lon1, lat2, lon2) RESULT(rdis)

    IMPLICIT NONE

    ! ----------------------
    ! Input/output variables
    ! -----------------------
    REAL (KIND=r8), INTENT (IN) :: lat1, lon1, lat2, lon2
    REAL (KIND=r8)              :: rdis

    ! Local variable
    REAL (KIND=r8)              :: dlon, dlat, a

    dlat = lat2 - lat1;    dlon = lon2 - lon1
    a = MIN(1.0, SQRT( SIN(dlat/2.0)**2.0 + COS(lat1) * COS(lat2) * SIN(dlon/2.0)**2.0 )  )
    rdis = 2.0 * ASIN(a)                                    ! relative distaince in radiances
 
    RETURN

  END FUNCTION circle_rdis

  ! This one works, gam is the longitude difference (-pi < gam0 < pi) between 2 and 1 (phi2 - phi1)
  SUBROUTINE sphergeom_intermediate ( lat1, lon1, lat2, lon2, c0, c, lat, lon )

    ! -----------------------------------------------------------------
    ! Finds the co-ordinates of C the baseline extended from two
    ! lon/lat points (A, B) on a sphere given the hypotenuse C_IN.
    ! ----------------------------------------------------------------
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    REAL (KIND=r8),    INTENT (IN) :: lat1, lat2, lon1, lon2, c0, c

    ! ----------------
    ! Output variables
    ! ----------------
    REAL (KIND=r8),  INTENT (OUT)  :: lat, lon

    ! ---------------
    ! Local variables
    ! ---------------
    REAL (KIND=r8)  :: x, y, z, tmp1, tmp2, frc, gamsign, theta, gam0, gam

    lat = 0.0_r8 ; lon = 0.0_r8    
    gam0 = angle_minus_twopi ( lon2 - lon1, pi )
    gamsign = ABS(gam0) / gam0;  gam0 = ABS(gam0)
    if (gam0 == 0 ) gamsign = 1
    ! Get straight line (AB) segment fraction frc intercepted by the line from center to C
    ! If frc < 0, extrapolation, but it is limited to |c| < (180-c0)/2.0
    tmp1 = SIN(c)
    frc = tmp1 / (SIN(c0 - c) + tmp1)
    ! Work in Cartesian Coordinate 
    tmp1 = frc * COS(lat2); tmp2 = 1.0 - frc
    x = tmp2 *   COS(lat1) + tmp1 * COS(gam0)
    y = tmp1 *   SIN(gam0)
    z = tmp2 *   SIN(lat1) + frc * SIN(lat2)

    gam = ATAN(y/x)                          ! -90 < gam < 90
    IF (frc >= 0) THEN
       IF (gam < 0) gam = gam + pi           ! 0 <= gam <= 180
    ELSE
       IF (gam > 0) gam = gam - pi           ! -180 <= gam <= 0
    ENDIF
    gam = gamsign * gam                      ! Get correct sign
    lon = gam + lon1

    theta = ATAN (SQRT(x**2 + y**2) / z)     ! -90 < theta < 90
    IF (theta < 0) theta = theta + pi        ! 0 <= theta <= 180
    lat = pihalf - theta                     ! Convert to latitude
    RETURN
  END SUBROUTINE sphergeom_intermediate

  SUBROUTINE sphergeom_baseline_comp ( a0, b0, gam0, c0 )
    ! -------------------------------------------------------
    ! Finds the lengh of the baseline of a spherical triangle
    ! -------------------------------------------------------
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    REAL (KIND=r8), INTENT (IN) :: a0, b0, gam0

    ! ---------------
    ! Output variable
    ! ---------------
    REAL (KIND=r8), INTENT (OUT) :: c0

    ! --------------
    ! Local variable
    ! --------------
    REAL (KIND=r8) :: tmp

    ! --------------------------------------------------
    ! Initialize output variables to keep compiler happy
    ! --------------------------------------------------
    c0 = 0.0_r8

    ! -------------------------------------------------
    ! Compute length of baseline between the two points
    ! -------------------------------------------------
    tmp = COS(a0) * COS(b0) + SIN(a0) * SIN(b0) * COS(gam0)
    IF ( ABS(tmp) < eps ) THEN
       c0 = pihalf
    ELSE
       c0 = ACOS(tmp)
    END IF

    RETURN
  END SUBROUTINE sphergeom_baseline_comp

  SUBROUTINE lonlat_to_pi ( lon, lat )

    IMPLICIT NONE

    ! ------------------
    ! Modified variables
    ! ------------------
    REAL (KIND=r8), INTENT (INOUT) :: lon, lat

    ! ------------------------------------
    ! Adjust longitude values to [-pi,+pi]
    ! ------------------------------------
    IF ( ABS(lon) > twopi ) lon = MOD(lon, twopi)
    IF ( lon >  pi ) lon = lon - twopi
    IF ( lon < -pi ) lon = lon + twopi
    ! ---------------------------------------
    ! Adjust latitude values to [-pi/2,+pi/2]
    ! ---------------------------------------
    IF ( ABS(lat) > pihalf ) lat = MOD(lat, pihalf)
    IF ( lat >  pihalf ) lat =   pi - lat
    IF ( lat < -pihalf ) lat = -(pi + lat)
    
    RETURN
  END SUBROUTINE lonlat_to_pi

  REAL (KIND=r8) FUNCTION angle_minus_twopi ( gamma0, pival ) RESULT ( gamma )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    REAL (KIND=r8), INTENT (IN) :: gamma0, pival

    IF ( gamma0 > pival ) THEN
       gamma = gamma0 - 2.0_r8 * pival !SIGN(2.0_r8*pival - gamma0, gamma0)
    ELSE IF ( gamma0 < -pival ) THEN
       gamma = gamma0 + 2.0_r8 * pival 
    ELSE
       gamma = gamma0
    END IF

    RETURN
  END FUNCTION angle_minus_twopi
END MODULE GEMS_o3P_geo_module
