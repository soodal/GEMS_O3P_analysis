MODULE SYNT_read_l1b
  USE SYNT_data_module
  USE NETCDF
  USE OMSAO_parameters_module,  ONLY: maxchlen,normweight
  USE GEMS_O3P_gemsdata_module, ONLY: gemsraddate

 
  CONTAINS
  SUBROUTINE get_synt_dims(fname, nx, ny, nw)
      IMPLICIT NONE
      !----------------------
      ! INPUT/OUTPUT variables
      !----------------------
      CHARACTER(maxchlen), INTENT(IN) :: fname
      INTEGER , INTENT(OUT) :: nx,  ny, nw
      nw = 1001
      nx = 30 
      ny = 899
    RETURN
  END SUBROUTINE get_synt_dims

  SUBROUTINE read_synt_l1b (rad_fname, pge_error_status)
    ! USE HDF5
    IMPLICIT NONE
    !------------------------
    ! INPUT/OUTPUT varialbes
    !------------------------
    CHARACTER(maxchlen), INTENT(IN) :: rad_fname 
    INTEGER,     INTENT(INOUT)      :: pge_error_status
    INTEGER :: ncid, varid    
    INTEGER :: i , j
    REAL, DIMENSION (nxtrack, nytrack) :: sza, vza, saa, vaa
    REAL, DIMENSION (nxtrack, nytrack) :: lon, lat, hter 
    REAL, DIMENSION (nxtrack, nytrack, nwavel) :: rad, snr
    REAL, DIMENSION (nytrack, nwavel) :: irrad, wave, wave_irr

    call check( nf90_open(adjustl(trim(rad_fname)), NF90_NOWRITE, ncid))

    !----------------------------
    ! Read Geolocation Field
    !---------------------------- 
    call check( nf90_inq_varid(ncid, "SZA", varid) ) ! nx, ny
    call check( nf90_get_var(ncid, varid, sza) )

    call check( nf90_inq_varid(ncid, "VZA", varid) ) ! nx, ny
    call check( nf90_get_var(ncid, varid, vza) )

    call check( nf90_inq_varid(ncid, "SAA", varid) ) ! nx, ny
    call check( nf90_get_var(ncid, varid, saa) )

    call check( nf90_inq_varid(ncid, "VAA", varid) ) ! nx, ny
    call check( nf90_get_var(ncid, varid, vaa) )

    call check( nf90_inq_varid(ncid, "Longitude", varid) ) ! nx
    call check( nf90_get_var(ncid, varid, lon) )

    call check( nf90_inq_varid(ncid, "Latitude", varid) ) ! nx, ny
    call check( nf90_get_var(ncid, varid, lat) )

    !call check( nf90_inq_varid(ncid, "TerrainHeight", varid) ) ! nx, ny
    !call check( nf90_get_var(ncid, varid, hter) )

    call check( nf90_inq_varid(ncid, "SNR", varid) ) ! nx, ny
    call check( nf90_get_var(ncid, varid, snr) )

    !----------------------------
    ! Read Radiance Field
    !----------------------------

    call check( nf90_inq_varid(ncid, "Wavelength", varid) ) !nx, ny, nw
    call check( nf90_get_var(ncid, varid, wave) )

    call check( nf90_inq_varid(ncid, "radiance", varid) ) ! nx, ny, nw
    call check( nf90_get_var(ncid, varid, rad) )

    !----------------------------
    ! Read Irradiance Field
    !----------------------------
    call check( nf90_inq_varid(ncid, "irradiance", varid) ) !ny
    call check( nf90_get_var(ncid, varid, irrad) )

    call check( nf90_inq_varid(ncid, "Wavelength_irr", varid) ) !nx, ny, nw
    call check( nf90_get_var(ncid, varid, wave_irr) )

    synt_szenith = sza
    synt_vzenith = vza
    synt_saa     = saa
    synt_vaa     = vaa
    synt_latitude  = lat    
    synt_Longitude = lon
    synt_hter      = hter

    DO i = 1, nxtrack
      synt_irrad_wavl(:,i) = wave_irr(100,:)
      synt_irrad_spec(:,i) = irrad(100,:)
      DO j = 1, nytrack
        synt_rad_wavl (:,i,j) = wave(j,:)
        synt_rad_snr (:,i,j) = snr(i,j,:)
        synt_rad_spec (:,i,j) = rad (i,j,:)
        !synt_rad_spec (:,i,j) = rad (i,j,:) * irrad(100,:)    ! synt data includes normalized radiance
      ENDDO
    ENDDO

    ! Read date from radiance file (used for correcting sun-earth distance when using backupirradiance)
    !i = INDEX(rad_fname, '.nc')-19   ! total data
    i = INDEX(rad_fname, '.nc')-23    ! sliced data
    gemsraddate = rad_fname(i : i + 3)//'m'//rad_fname(i+4 : i+7)
    RETURN
  END SUBROUTINE read_synt_l1b
  

  SUBROUTINE read_synt_alb (alb_fname, synt_surfalb, pge_error_status)
    IMPLICIT NONE
    !------------------------
    ! INPUT/OUTPUT varialbes
    !------------------------
    CHARACTER(maxchlen), INTENT(IN) :: alb_fname
    INTEGER,     INTENT(INOUT)      :: pge_error_status
    INTEGER :: ncid, varid    
    INTEGER :: i , j
    REAL, DIMENSION (nxtrack, nytrack) :: alb
    REAL, DIMENSION (nxtrack, nytrack, nwavel) :: synt_surfalb

    !----------------------------
    ! OPEN
    !---------------------------- 
    call check( nf90_open(adjustl(trim(alb_fname)), NF90_NOWRITE, ncid)) 

    !----------------------------
    ! Read Geolocation Field
    !---------------------------- 
    call check( nf90_inq_varid(ncid, "SurfaceAlbedo", varid) ) ! nx, ny
    call check( nf90_get_var(ncid, varid, alb) )

    DO j = 1, nytrack
    DO i = 1, nxtrack
      synt_surfalb(i,j,:)=alb(i,j)
    ENDDO
    ENDDO

    RETURN
  END SUBROUTINE read_synt_alb



  SUBROUTINE check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  END SUBROUTINE check  

END MODULE SYNT_read_l1b
