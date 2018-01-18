MODULE OMSAO_gome_data_module

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: maxchlen, max_fit_pts, maxwin
  IMPLICIT NONE


  ! ===================================
  ! GOME instrument and data parameters
  ! ===================================
  INTEGER, PARAMETER :: n_gome_pmd = 3, n_gome_pmd_pix = 16, n_gome_data_dim = 5
       
  ! ----------------------------------------------------------
  ! Index for the GOME East, Center, West, and Back-Scan Pixel
  ! ----------------------------------------------------------
  INTEGER, PARAMETER :: &
       gome_east_idx = 0, gome_cent_idx = 1, gome_west_idx = 2, gome_bcks_idx = 3

  ! =========================================================
  ! Landmarks to look for when reading extracted Level-1 file
  ! =========================================================
  ! --------------------------------------
  ! First some character length parameters
  ! --------------------------------------
  INTEGER, PARAMETER :: &
       gome_lm_chan_len =  9, gome_lm_gp_len  = 12, gome_lm_pmd_len = 3, &
       gome_lm_band_len =  7, gome_lm_ers_len = 15, gome_lm_esh_len = 19, &
       gome_lm_sols_len = 14
  ! ----------------------------------------------------------
  ! Now the character land marks (lm) in the extracted L1 file
  ! ----------------------------------------------------------
  CHARACTER (LEN = gome_lm_chan_len), PARAMETER :: &
       lm_gome_chan_1  = 'Channel 1', lm_gome_chan_2  = 'Channel 2', &
       lm_gome_chan_3  = 'Channel 3', lm_gome_chan_4  = 'Channel 4'   
  CHARACTER (LEN = gome_lm_band_len), PARAMETER :: &
       lm_gome_band_1  = 'Band 1 ', &
       lm_gome_band_2a = 'Band 2a', lm_gome_band_2b = 'Band 2b', &
       lm_gome_band_3  = 'Band 3 ', lm_gome_band_4  = 'Band 4 '
  CHARACTER (LEN = gome_lm_sols_len), PARAMETER :: lm_gome_solspec   = 'Solar Spectrum'
  CHARACTER (LEN = gome_lm_gp_len  ), PARAMETER :: lm_gome_groundpix = 'Ground Pixel'
  CHARACTER (LEN = gome_lm_pmd_len ), PARAMETER :: lm_gome_pmd       = 'PMD'
  CHARACTER (LEN = gome_lm_ers_len ), PARAMETER :: lm_gome_ersinfo   = 'ERS Information'
  CHARACTER (LEN = gome_lm_esh_len ), PARAMETER :: lm_gome_eshine    = 'Earthshine Spectrum'


  ! ====================================================
  ! Variables associated with GOME ground pixel numbers, 
  ! wavelength coverage, and other things
  ! ====================================================
  ! Maximum number of spectral points in a GOME channel
  INTEGER, PARAMETER :: n_gome_max_pts = 1024
  ! Flag for missing spectral channel data
  INTEGER, PARAMETER :: gome_spec_missing = 0
  ! Indices for Zenith and Azimuth angles
  INTEGER, PARAMETER :: zen0_idx = 1, azm0_idx = 2, zen_idx = 3, azm_idx = 4, los_idx = 2
  ! Indices for Latitudes and Longitudes
  INTEGER, PARAMETER :: lat_idx = 1, lon_idx = 2
  ! Number of geolocation and angle measurements per pixel
  INTEGER, PARAMETER :: n_gome_ang = 3, n_gome_geo = 5
  

  ! =================
  ! Variables section
  ! =================

  ! ------------------------------------------------------------
  ! GOME data filenames - temporary until OMI data are available
  ! ------------------------------------------------------------
  CHARACTER (LEN=maxchlen) :: gome_data_input_fname, gome_data_output_fname

  INTEGER :: obs_month
  INTEGER :: n_gome_solpts, gome_curpix, gome_curqual, gome_curscan, n_gome_radpts
  INTEGER :: gome_stpix, gome_endpix, gome_npix
  INTEGER :: orbnum                ! orbit number from the launch
  CHARACTER (LEN=24) :: gome_pixdate
  CHARACTER (LEN=3)  :: gome_orbc  ! one of the 14 orbits in a day
  REAL (KIND=dp)                                           :: ers2_alt, earth_curv, gome_integt
  REAL (KIND=dp), DIMENSION (azm_idx,n_gome_ang)           :: gome_angles_wrtn, gome_angles_wrts
  REAL (KIND=dp), DIMENSION (lon_idx,n_gome_geo)           :: gome_geoloc
  REAL (KIND=dp), DIMENSION (n_gome_data_dim, max_fit_pts) :: gome_radspec, gome_solspec

  ! Maximum ratio of # channel 2 pixels to # channel 1 ratios
  INTEGER, PARAMETER :: maxc1c2r = 8

  ! Index of the current pixels in those ch2 pixels
  INTEGER            :: thecurpix  

  CHARACTER (LEN=24), DIMENSION(maxc1c2r) :: allpixdate
  INTEGER           , DIMENSION(maxc1c2r) :: allcurscan, allcurpix
  REAL (KIND=dp), DIMENSION (maxc1c2r, lon_idx, n_gome_geo)  :: allgeoloc 
  REAL (KIND=dp), DIMENSION (maxc1c2r) :: allsza, allvza, allaza, allsca
  REAL (KIND=dp), DIMENSION (maxc1c2r, maxwin, n_gome_max_pts, n_gome_data_dim) :: allspec
  INTEGER, DIMENSION(maxc1c2r, maxwin) :: nwpos, stpos

END MODULE OMSAO_gome_data_module
