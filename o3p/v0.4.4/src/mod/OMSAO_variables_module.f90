 ! *********************** Modification History *******************
  ! xiong liu, July 2003
  ! 1. Add three variables: the_sza_atm, the_vza_atm, the_aza_atm 
  !    Effective siewing geometry at TOA averaged on A, B, C for
  !    inputting into LIDORT
  ! 2. Add the_lat, the_lon (used in preparing atmospheric profiles)
  ! 3. Add yn_varyslit, hwlarr, hwrarr, vglarr, vgrarr, slitwav, n_slit_pts,
  !    n_slit_interval, slit_fname, slit_redo, wavcal_redo, 
  !    wavcal_fname, shiarr, squarr, sswav, for implementing 
  !    variable slit width ( for voigt profile shape only)
  ! 5. Add varaible use_meas_sig (use measurement error, otherwise
  !    use normal weight = 1 for all pixels except edge pixels)
  ! 6. Add variable for using pixel_bin, else use interpoaltion
  ! 7. Add # of wavelengths in Channel 1 and 2 respectively
  ! ****************************************************************

MODULE OMSAO_variables_module

  USE OMSAO_precision_module
  USE OMSAO_indices_module, ONLY: &
       max_rs_idx, max_calfit_idx, n_max_fitpars, mxs_idx, sig_idx, icf_idx, &
       amf_idx, n_amftab_ang_max,  n_amftab_dim_max
  USE OMSAO_parameters_module,   ONLY: &
       maxchlen, max_spec_pts, max_fit_pts, n_sol_winwav, n_rad_winwav,      &
       max_mol_fit, maxwin, maxview, maxloc, max_ref_pts

  IMPLICIT NONE


  ! -------------------------------------
  ! GOME data fitting or OMI data fitting
  ! -------------------------------------
  INTEGER, PARAMETER :: omi_idx = 1, gome_idx = 2, scia_idx = 3, gome2_idx = 4,gems_idx = 5
  INTEGER, PARAMETER :: max_instrument_idx = gems_idx
  CHARACTER (LEN=4), DIMENSION (max_instrument_idx), PARAMETER :: &
       which_instrument = (/ 'OMI ', 'GOME', 'SCIA', 'GOM2','GEMS' /)
  INTEGER :: instrument_idx

  ! Time of a Granule
  REAL (KIND=8) :: TAI93At0ZOfGranule, TAI93StartOfGranule, GranuleSecond
  INTEGER       :: GranuleYear, GranuleMonth, GranuleDay, GranuleHour,       &
       GranuleMinute, GranuleJDay

  CHARACTER (LEN=24) :: currtime
  ! -----------------------------
  ! Variables read from PCF file
  ! -----------------------------
  ! * Current PGE name and index
  ! -----------------------------
  INTEGER           :: pge_idx
  CHARACTER (LEN=6) :: pge_name

  ! ---------------------
  ! * Verbosity threshold
  ! ---------------------
  CHARACTER (LEN=1) :: verb_thresh_char
  INTEGER           :: verb_thresh_lev
  ! -------------------------------
  ! * Swath name (really from PCF?)
  ! -------------------------------
  CHARACTER (LEN=maxchlen) :: swath_name
  ! --------------
  ! * Orbit number
  ! --------------
  INTEGER          :: pcf_orbit_number

  ! -------------------------------------------------
  ! Variables defined in preamble of original program
  ! -------------------------------------------------
  LOGICAL :: yn_smooth, yn_doas, yn_varyslit, use_meas_sig
  LOGICAL :: correct_merr              ! Correct OMI COL3 measurement error, xliu: 09/25/12
  LOGICAL :: xbin_decerr, ybin_decerr  ! Reduce meas error when coadding in x/y direction, xliu: 09/25/12
  LOGICAL :: weight_sun, weight_rad, renorm

  ! -------------------------------------------
  ! A special beast: The undersampling spectrum
  ! -------------------------------------------
  LOGICAL :: have_undersampling

  ! xliu, 01/03/2007, add variables for reduce spectral resolution
  LOGICAL                                :: reduce_resolution, use_redfixwav
  INTEGER                                :: reduce_slit, nredfixwav
  REAL (KIND=dp)                         :: redslw, redsampr, redlam
  REAL (KIND=dp), DIMENSION(max_fit_pts) :: redfixwav
  CHARACTER (LEN=maxchlen)               :: redfixwav_fname

  ! Number of pixels read and number of pixels with successful retrieval
  INTEGER :: npix_fitting, npix_fitted
    
  INTEGER                                       :: n_fitvar_rad, n_fitvar_sol
  INTEGER,           DIMENSION (max_calfit_idx) :: mask_fitvar_sol, rmask_fitvar_sol
  REAL (KIND=dp),    DIMENSION (max_calfit_idx) :: fitvar_sol, fitvar_sol_init,&
       fitvar_sol_saved, lo_sunbnd, up_sunbnd, lo_sunbnd_init, up_sunbnd_init
  CHARACTER (LEN=6), DIMENSION (max_calfit_idx) :: fitvar_sol_str
  
  INTEGER,           DIMENSION (n_max_fitpars)  :: mask_fitvar_rad, rmask_fitvar_rad, &
       database_indices, fothvarpos
  REAL (KIND=dp),    DIMENSION (n_max_fitpars)  :: fitvar_rad_init, fitvar_rad, &
       fitvar_rad_apriori, fitvar_rad_aperror, fitvar_rad_saved, fitvar_rad_std, &
       lo_radbnd, up_radbnd, lo_radbnd_init, up_radbnd_init, fitvar_rad_nstd, fitvar_rad_init_saved
  CHARACTER (LEN=6),  DIMENSION (n_max_fitpars)  :: fitvar_rad_str
  CHARACTER (LEN=15), DIMENSION (n_max_fitpars)  :: fitvar_rad_unit

  ! fitspec_rad: fitted spectrum    (I/F, after removing non-ozone and albedo compoments)
  ! simspec_rad: simulated spectrum (only ozone and albedo terms)
  ! fitres_rad : fitspec_rad - simspec_rad
  ! actspec_rad: I/F (without removing anything)
  ! clmspec_rad: (simulated spectrum, ozone and albedo terms only, but with a priori climatology)
  REAL (KIND=dp),    DIMENSION (max_fit_pts)    :: fitspec_rad, fitres_rad, actspec_rad, simspec_rad, clmspec_rad
  
  REAL (KIND=dp), DIMENSION (max_rs_idx, max_ref_pts) :: database, database_shiwf, database_save
  REAL (KIND=dp), DIMENSION (max_ref_pts)             :: slwf, dfdsl
  REAL (KIND=dp), DIMENSION (max_ref_pts)             :: refwvl, refwvl_sav
  INTEGER, DIMENSION (max_ref_pts)                    :: refidx, refsol_idx
  INTEGER, DIMENSION (max_ref_pts)                    :: refidx_sav

  INTEGER                                             :: n_refwvl, n_refwvl_sav, the_nspc
  REAL (KIND=dp), DIMENSION (3)     :: sza_atm, vza_atm, aza_atm
  REAL (KIND=dp)                    :: the_sza_atm, the_vza_atm, the_aza_atm, &
       the_sca_atm, the_lat, the_lon, the_surfalt, avgsza, avgvza, avgaza, avgsca
  INTEGER                           :: the_month, the_year, the_day 
  INTEGER                           :: nview, nloc
  REAL (KIND=dp)                    :: the_time
  REAL (KIND=dp), DIMENSION(maxloc) :: the_lons, the_lats
  REAL (KIND=dp), DIMENSION(2)      :: edgelons, edgelats 
  CHARACTER (LEN = 28)              :: the_utc

        
  ! --------------------------------------------------------
  ! Identifier string for irradiacne and radiance input file
  ! --------------------------------------------------------
  CHARACTER(LEN=6)         :: sol_identifier
  CHARACTER(LEN=6)         :: rad_identifier
  CHARACTER(LEN=maxchlen)  :: ch1_out_file
  CHARACTER(LEN=maxchlen)  :: outdir
  
  ! -------------------------------------
  ! Variables related to Air Mass Factors
  ! -------------------------------------
  INTEGER            :: n_amftab_dim, n_amftab_ang
  LOGICAL            :: have_amf, have_amftable
  REAL    (KIND=dp)  :: amfgeo, amf, sol_zen_eff
  REAL    (KIND=dp)  :: amf_esza_min, amf_esza_max
  REAL    (KIND=dp), DIMENSION (n_amftab_ang_max, n_amftab_dim_max) :: amf_table

  ! -----------------------------
  ! Previously IMPLICIT variables
  ! -----------------------------
  REAL (KIND=dp) :: phase, szamax, zatmos, chisq, sol_wav_avg, rad_wav_avg

  ! ------------------------------------------------
  ! File names for solar and other reference spectra
  ! ------------------------------------------------
  CHARACTER (LEN=maxchlen) :: fitctrl_fname

  ! ---------------------------------------------------------
  ! Directory for reference spectra and atmospheric databases
  ! ---------------------------------------------------------
  CHARACTER (LEN=maxchlen) :: atmdbdir, refdbdir

  ! -----------------------------------------------------------
  ! Variables related to reference spectra
  !  * number of reference spectra:       N_REfSPEC
  !  * indentification strings:           FITPAR_IDXNAME
  !  * file names with reference spectra: REFSPEC_FNAME
  !  * original (uniterpolated) data:     REFSPEC_ORIG_DATA
  !  * number of spectral points:         N_REFSPEC_PTS
  !  * first and last wavelenghts:        REFSPEC_FIRSTLAST_WAV
  ! -----------------------------------------------------------
  INTEGER                                                         :: n_refspec
  CHARACTER (LEN=6),        DIMENSION (max_rs_idx)                :: fitpar_idxname
  CHARACTER (LEN=maxchlen), DIMENSION (max_rs_idx)                :: refspec_fname
  REAL (KIND=dp),           DIMENSION (max_rs_idx,max_spec_pts,3) :: refspec_orig_data
  REAL (KIND=dp),           DIMENSION (max_spec_pts)              :: solar_refspec
  REAL (KIND=dp),           DIMENSION (max_rs_idx,2)              :: refspec_firstlast_wav
  REAL (KIND=dp),           DIMENSION (max_rs_idx)                :: refspec_norm
  INTEGER,                  DIMENSION (max_rs_idx)                :: n_refspec_pts

  ! Special for ozone, use Tdependent coefficients or at several T
  !INTEGER, PARAMETER :: maxozabs = 5
  !INTEGER            :: numozabs = 0, n_ozref_pts = 0
  !LOGICAL            :: oztdepend, ozconv
  !REAL (KIND=dp)            :: ozrefspec_norm
  !REAL (KIND=dp), DIMENSION (maxozabs)                  :: ozrefts
  !REAL (KIND=dp), DIMENSION (maxozabs, max_spec_pts)    :: ozrefspec
  !REAL (KIND=dp), DIMENSION (max_spec_pts) :: ozrefpos
  !REAL (KIND=dp), DIMENSION (maxozabs, max_fit_pts+4)   :: ozdb, ozdb_shiwf


  CHARACTER (LEN=maxchlen), DIMENSION (icf_idx:amf_idx) :: static_input_fnames

  REAL (KIND=dp), DIMENSION (max_spec_pts) :: cubic_x, cubic_y, cubic_w
  REAL (KIND=dp), DIMENSION (max_spec_pts) :: poly_x, poly_y, poly_w
  INTEGER                                  :: poly_order

  REAL (KIND=dp), DIMENSION (max_spec_pts) :: step2_y
  REAL (KIND=dp), DIMENSION (max_spec_pts, max_fit_pts) :: step2_dyda
  
  ! --------------------------------------
  ! Solar and Earth shine wavlength limits
  ! --------------------------------------
  REAL (KIND=dp), DIMENSION (n_sol_winwav) :: sol_winwav_lim
  REAL (KIND=dp), DIMENSION (n_rad_winwav) :: rad_winwav_lim
  REAL (KIND=dp)                           :: winwav_min, winwav_max

  ! ------------------------------------------------------------------
  ! Indices of fitting window defining wavelengths in current spectrum
  ! ------------------------------------------------------------------
  INTEGER, DIMENSION (n_rad_winwav)        :: rad_winwav_idx

  ! --------------------------------------------------------------------------
  ! The current solar and radiance spectrum, including wavelengths and weights
  ! --------------------------------------------------------------------------
  INTEGER  :: currloop, currpix, currline, currtrack
  CHARACTER (LEN=3) :: currpixchar
  REAL (KIND=dp), DIMENSION (sig_idx, max_fit_pts) :: curr_rad_spec
  REAL (KIND=dp), DIMENSION (sig_idx, max_fit_pts) :: curr_sol_spec
  REAL (KIND=dp), DIMENSION (max_fit_pts):: fitwavs, fitweights, currspec

  INTEGER                                       :: nsol_ring, sring_fidx, sring_lidx
  REAL (KIND=dp), DIMENSION (2, max_fit_pts)    :: sol_spec_ring

  ! ----------------------
  ! variable slit width
  ! ----------------------
  REAL (KIND=dp)                          :: slit_trunc_limit
  REAL (KIND=dp), DIMENSION (max_fit_pts) :: slitwav, slitwav_rad=0., &
       slitwav_sol=0.0, sswav_rad=0.0, sswav_sol=0.0, slitdis=0.0
  REAL (KIND=dp), DIMENSION (max_fit_pts, max_calfit_idx, 2) :: solslitfit=0.0, &
       radslitfit=0.0, solwavfit=0.0, radwavfit=0.0, slitfit=0.0
  INTEGER  :: slit_fit_pts, n_slit_step, nslit,nslit_rad, nslit_sol, &
       wavcal_fit_pts, n_wavcal_step, nwavcal_sol, nwavcal_rad
  CHARACTER (LEN=maxchlen)  :: slit_fname, rslit_fname, swavcal_fname, wavcal_fname
  LOGICAL                   :: slit_redo, wavcal_redo, wavcal_sol, wavcal, slit_rad
  LOGICAL                   :: fixslitcal, smooth_slit, slitcal
  INTEGER                   :: which_slit   ! 1. Gauss 2. Voigt 3. Triangle Other: Gauss


  ! hw1e, e_asym, shi, squ, hwl, hwr, vgl, vgr at each window 
  ! (value, standard deviation)
  REAL (KIND=dp), DIMENSION(maxwin,max_calfit_idx,2) :: solwinfit, radwinfit
  REAL (KIND=dp), DIMENSION(maxwin)                  :: wincal_wav 
  
  ! ----------------------
  ! calibration
  ! ----------------------
  INTEGER                                :: numwin
  REAL (KIND=dp), DIMENSION(maxwin, 2)   :: winlim
  INTEGER, DIMENSION(maxwin, 2)          :: winpix
  INTEGER, DIMENSION(maxwin)             :: nradpix, nsolpix, &
       n_band_avg, n_band_samp, nradpix_sav,  band_selectors

  ! radnhtrunc: number of unused radiance pixels at each end of a spectralregion
  ! refnhextra: number of extra pixels added to the reference spectra
  ! The main purpose is to avoid extrapolation

  INTEGER                                :: radnhtrunc, refnhextra
  LOGICAL                                :: do_bandavg
  INTEGER                                :: n_radwvl_sav
  REAL (KIND=dp), DIMENSION(max_fit_pts) :: radwvl_sav, i0sav
  

  ! Determine whether to coadd UV2 spectrum to the UV-1 resolution
  ! When both UV-1 and UV-2 exist
  LOGICAL                                :: coadd_uv2

  ! Use solar composite for destriping
  LOGICAL                                ::use_backup, use_solcomp, avg_solcomp, avgsol_allorb  

  ! Whether to perform wavelength calibration (due to spatial smile) before spatial coadding 
  LOGICAL                                :: wcal_bef_coadd
  
  ! Whehter to filter spectral pixels around 280 and 285 nm
  LOGICAL                                :: rm_mgline
  REAL (KIND=dp)                         :: dwavmax

  ! Need to convolve high-resolution ozone absorption cross section at the beginning of 
  ! processing each scan position of each block
  ! Once the xsection is convolved, it will be set to false in ROUTINE getabs_crs
  LOGICAL                                :: ozabs_convl, so2crs_convl

  ! Writing limited output to screen for debugging and testing
  LOGICAL                                :: scnwrt
  
  ! -------------------------------------------
  !band 1a and 1b boundary
  ! -------------------------------------------
  REAL (KIND=dp) :: b1ab_div_wav 

  ! ---------------------------------------------
  ! Ground Scan Lines and xtrack pixels Limits
  ! --------------------------------------------
  INTEGER, DIMENSION (2) :: linenum_lim          ! scan lines or pixel numbers
  INTEGER, DIMENSION (2) :: pixnum_lim           ! across the track pixels

  ! --------------------------------------------
  ! Frequency of radiance wavelength calibration
  ! --------------------------------------------
  INTEGER :: radwavcal_freq

  ! ------------------------------------------------------------------------
  ! Variables connected with ELSUNC numerical precision/convergence criteria
  ! ------------------------------------------------------------------------
  REAL (KIND=dp) :: tol,  epsrel,  epsabs,  epsx

  ! ----------------------------------------
  ! Variable for +1.0 or -1.0 multiplication
  ! ----------------------------------------
  REAL (KIND=dp) :: pm_one

  ! ----------------------------------------------------------------------
  ! Index for the fitting parameters carrying the fitted column value.
  !
  ! N_MOL_FIT:    Number of "molecules" that carry the final column; this
  !               can be one molecule at different temperatures.
  ! FITCOL_IDX:   The main molecule indices, corresponding to the list of
  !               reference spectra.
  ! FINCOL_IDX:   For the final summation of the fitted column: The total
  !               number is the number of different molecules times the
  !               allowed sub-indices. The second dimension is for the
  !               reference spectrum index - this eases the final sum over
  !               the fitted columns (see RADIANCE_FIT subroutine).
  !                includes
  !               the subindices, hence the dimension.
  ! N_FINCOL_IDX: Number of final column indices.
  ! ----------------------------------------------------------------------
  INTEGER                                    :: n_mol_fit, n_fincol_idx
  INTEGER, DIMENSION (max_mol_fit)           :: fitcol_idx
  INTEGER, DIMENSION (2,max_mol_fit*mxs_idx) :: fincol_idx

  ! ------------------------------------
  ! Maximum number of fitting iterations
  ! ------------------------------------
  INTEGER :: max_itnum_sol, max_itnum_rad

  ! ---------------------
  ! L1B and L2 file names
  ! ---------------------
  INTEGER :: l2_hdf_flag
  CHARACTER (LEN=maxchlen) :: l1b_rad_filename, l1b_irrad_filename, l2_filename, &
       l2_cld_filename, l2_swathname

  ! -----------------------------------------------------------------
  ! Generic dimension variables (initialized from either GOME or OMI)
  ! -----------------------------------------------------------------
  INTEGER :: n_irrad_wvl, n_rad_wvl

  ! -----------------------------------------------------------------
  !
  ! This module defines variables associated with the SAO PGEs.
  !
  ! -----------------------------------------------------------------

  ! debug variable
  LOGICAL :: debug_boreas

  ! Define SAA boundary 
  REAL (KIND=sp) :: saa_minlon = -75.0, saa_maxlon = 0.00, &
                    saa_minlat = -50.0, saa_maxlat = -5.00, &
                    saa_minlon1= 0.00 , saa_maxlon1= 30.0, &
                    saa_minlat1= -35.0, saa_maxlat1= -15.0

END MODULE OMSAO_variables_module
