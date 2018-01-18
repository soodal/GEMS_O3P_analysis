  ! *********************** Modification History ********
  ! 1. xliu: specify maximum of layers in ozone profile
  ! 2. xliu: specify the tolerance value used in FMIN
  ! *****************************************************

MODULE OMSAO_parameters_module

  USE OMSAO_precision_module
  IMPLICIT NONE

  ! ================================================
  ! Define array dimensions and other parameters for
  ! SAO Trace Gas Fitting Algorithms
  ! ================================================

  ! -----------------------------------
  ! Maximum length of CHARACTER strings
  ! -----------------------------------
  INTEGER, PARAMETER :: maxchlen = 256

  ! ----------------------------------
  ! Maximum iteration number for loops
  ! ----------------------------------
  INTEGER, PARAMETER :: forever = 9999, max_iter=50

  ! -------------------------------------------------------------------------
  ! Maximum numbers for fitting parameters, GOME pixels, spectral points, ...
  ! -------------------------------------------------------------------------
  INTEGER, PARAMETER :: max_spec_pts = 13505             ! Original reference spectrum
  INTEGER, PARAMETER :: max_fit_pts  = 800               ! fitted array
  INTEGER, PARAMETER :: max_ref_pts  = max_fit_pts + 20  ! convolved reference spectrum
  INTEGER, PARAMETER :: max_ring_pts = max_fit_pts + 50  ! solar spectrum for ring calculation
  INTEGER, PARAMETER :: maxloc       =  5  ! # of points to describe a pixel
  INTEGER, PARAMETER :: maxview      =  5  ! # of view angles to describe a pixel

  ! -----------------------------------------------------------------------------
  ! Maximum numbers for layers in retrievals and in radiative transfer calculation
  ! ------------------------------------------------------------------------------
  INTEGER, PARAMETER        :: maxlay = 50, mflay = 60 
  REAL (KIND=dp), PARAMETER :: lcurve_tol = 1.0D-8, smallval = 1.0E-10_dp
  
  ! -----------------------------------------------
  ! Number of other gases that is retrieved with o3
  ! useing assumed profile shape
  ! -----------------------------------------------
  INTEGER, PARAMETER :: maxog = 10
  
  ! -----------------------------------------------
  ! Number of windows to be used: 1, 2, 3, 4
  ! -----------------------------------------------
  INTEGER, PARAMETER :: maxband = 2   ! UV-1, UV-2
  INTEGER, PARAMETER :: maxwin = 5

  ! number of wavelengths around 370 nm used for calculating cloud fraction
  INTEGER, PARAMETER :: mrefl = 100
  
  ! ==================
  ! Physical constants
  ! ==================
  REAL (KIND=dp), PARAMETER :: pi      = 3.14159265358979_dp
  REAL (KIND=dp), PARAMETER :: deg2rad = pi / 180.0_dp, rad2deg = 180.0_dp / pi
  REAL (KIND=dp), PARAMETER :: du2mol  = 2.686763D16

  REAL (KIND=dp), PARAMETER :: p0      = 1013.25         ! Default surface pressure (mb)
  REAL (KIND=dp), PARAMETER :: boltz   = 1.3806505D-23   ! Boltzmann's constant (J/K)
  REAL (KIND=dp), PARAMETER :: xmair   = 28.97           ! Mean molecular mass of air (g/mol)
  REAL (KIND=dp), PARAMETER :: accgrav = 9.80665         ! Acceleration due to gravity (m/s/s)
  REAL (KIND=dp), PARAMETER :: ugc     = 8.314472        ! Universal gas constant (J/mol/K)
  REAL (KIND=dp), PARAMETER :: avo     = 6.0221415D23    ! Avogadro's number (molecules/mol)
  REAL (KIND=dp), PARAMETER :: rearth  = 6367.45         ! Mean radius of Earth (km)
  REAL (KIND=dp), PARAMETER :: zerok   = 273.15          ! Mean radius of Earth (km)
  REAL (KIND=dp), PARAMETER :: O2mix= 0.20949858, N2mix= 0.78079469, CO2mix=0.0003668  ! Air composition

  ! =============================================================
  ! Large weight for wavelengths to be excluded from the fitting;
  ! small weight for wavelengths to be included in the fitting
  ! =============================================================
  REAL (KIND=dp), PARAMETER :: downweight = 1.0E+10_dp, normweight = 1.0_dp


  ! =========================================
  ! Order of polynomial for DOAS baseline fit
  ! =========================================
  INTEGER, PARAMETER :: doas_npol = 4

  ! ===================================================
  ! Dimension parameters for ELSUNC auxiliary variables
  ! ===================================================
  INTEGER, PARAMETER :: elsunc_np = 11, elsunc_nw = 6

  ! ===============================================
  ! Number of Solar and Earthshine fitting windows.
  ! ===============================================
  INTEGER, PARAMETER :: n_sol_window = 1, n_rad_window = 1

  ! ===============================================================
  ! Number of wavelengths defining the Solar and Earthshine fitting
  ! windows. N windows require 2*N+2 bounding wavelenghts: Two each
  ! for each window, plus first and last wavlength
  ! ===============================================================
  INTEGER, PARAMETER :: n_sol_winwav = 2*n_sol_window+2
  INTEGER, PARAMETER :: n_rad_winwav = 2*n_rad_window+2

  ! ----------------
  ! ZERO_SPEC string
  ! ----------------
  CHARACTER (LEN=9), PARAMETER :: zerospec_string = 'Zero_Spec'

  ! -----------------------------------------------------------------
  ! Maximum number of molecules to fit. This number can be greater
  ! than ONE: We may include multiple reference spectra for the same
  ! molecule in the fit, e.g., at different temperatures. In this
  ! case we need to collect the fitted columns from multiple indices.
  ! -----------------------------------------------------------------
  INTEGER, PARAMETER :: max_mol_fit = 3

  ! =======================================================================
  ! Missing value for column and uncertainty in cases of non-converged fits
  ! =======================================================================
  REAL    (KIND=dp), PARAMETER :: missing_value_dp = -9.9999E+02_dp
  REAL    (KIND=sp), PARAMETER :: missing_value_sp = -9.9999E+02_sp
  INTEGER (KIND=i4), PARAMETER :: missing_value_i4 = -99999_i4

  ! -------------------------
  ! "start of table" landmark
  ! -------------------------
  CHARACTER (LEN=39), PARAMETER :: &
       lm_start_of_table = "start of table (don't delete this line)"


  ! -----------
  ! Fill values
  ! -----------
  !INTEGER (KIND=i1), PARAMETER :: int8_fill  = -127,        uint8_fill  = 255
  !INTEGER (KIND=i2), PARAMETER :: int16_fill = -32767,      uint16_fill = 65535
  !INTEGER (KIND=i4), PARAMETER :: int32_fill = -2147483647, uint32_fill = 4294967295

  REAL (KIND=r4),  PARAMETER :: float32_fill = -1.0E+30_r4
  REAL (KIND=r8),  PARAMETER :: float64_fill = -1.0E+30_r8

  ! --------------------------
  ! Verbosity Threshold Levels
  ! --------------------------
  ! vb_lev_default   0  default (one line only, at completion)
  ! vb_lev_omidebug  1  OMI support team debugging
  ! vb_lev_develop   3  PGE development
  ! vb_lev_lt1mb     4  LOG files of size <= 1MB
  ! vb_lev_gt1mb     5  LOG files of size >  1MB
  ! vb_lev_screen    6  Screen output (normally not desired)
  INTEGER, PARAMETER :: &
       vb_lev_default = 0, vb_lev_omidebug = 1, vb_lev_develop = 3, vb_lev_1mb = 4, &
       vb_lev_gt1mb = 5, vb_lev_screen = 6

  ! -------------------------------------------------------------------------------
  ! Missing values: We define two sets of the same value, but with different names.
  ! --------------- one set follows the naming convention of the OMI-GDPS-IODS
  ! document (Table 4-14), the other one follows more closely the KIND definition
  ! used in the present PGEs.
  !
  ! NOTE that the values for R4_MISSVAL and R8_MISSVAL are identical. This is for
  ! the case of quantities like the Effective Solar Zenith Angle, which are defined
  ! as R8 in the PGE but are written as R4 to the output file. If R8_MISSVAL was a
  ! truly R8 value (e.g., -1.0E+300), then the conversionto R4 and the HE5 write
  ! would fail.
  !
  ! -------------------------------------------------------------------------------
  CHARACTER (LEN=9),   PARAMETER :: str_missval     = "undefined"
  INTEGER   (KIND=i1), PARAMETER :: int8_missval    = -100         !-127
  INTEGER   (KIND=i2), PARAMETER :: int16_missval   = -30000       !-32767
  INTEGER   (KIND=i4), PARAMETER :: int32_missval   = -2000000000  !-2147483647
  REAL      (KIND=r4), PARAMETER :: float32_missval = -1.0E+30_r4  !-1.0_r4*(2.0_r4**100)
  REAL      (KIND=r8), PARAMETER :: float64_missval = -1.0E+30_r8  !-HUGE(1.0_r8) !-1.0_r8*(2.0_r8**100)

  INTEGER   (KIND=i1), PARAMETER :: i1_missval = int8_missval
  INTEGER   (KIND=i2), PARAMETER :: i2_missval = int16_missval
  INTEGER   (KIND=i4), PARAMETER :: i4_missval = int32_missval
  REAL      (KIND=r4), PARAMETER :: r4_missval = float32_missval
  REAL      (KIND=r8), PARAMETER :: r8_missval = float64_missval

  ! --------------------------------------------------
  ! Some generic Valid Ranges of the output quantities
  ! --------------------------------------------------
  REAL (KIND=r8), PARAMETER :: &
       valid_min_r8 = -9.9E+99_r8, valid_max_r8 = +9.9E+99_r8, zero_r8 = 0.0_r8, one_r8 = 1.0_r8, &
       valid_max_i2 = REAL(HUGE(1_i2), KIND=r8), valid_max_i4 = REAL(HUGE(1_i4), KIND=r8)
  
  ! ---------------------------------------
  ! ELSUNC Abnormal Termination Exit Values
  ! ---------------------------------------
  INTEGER (KIND=i2), PARAMETER :: &
       elsunc_nostart_eval =    -1, &  ! Wrong dimensions or wrong starting point
       elsunc_maxiter_eval =    -2, &  ! Interation exceeded maximum allowed number
       elsunc_hessian_eval =    -3, &  ! Interation exceeded maximum allowed number
       elsunc_no2deri_eval =    -4, &  ! No Second Derivative
       elsunc_newtstp_eval =    -5, &  ! Undamped Newton Step is a failure
       elsunc_nodecst_eval =    -6, &  ! Last step was not descending
       elsunc_onesolu_eval =    -7, &  ! Only one feasible point
       elsunc_parsoob_eval =   -11, &  ! User-defined: Fitting parameters out of bounds
       elsunc_infloop_eval =   -12, &  ! User-defined: Hit "infinite loop" snag
       elsunc_highest_eval = 12344, &  ! Largest possible exit value
       elsunc_usrstop_eval = elsunc_infloop_eval
  REAL    (KIND=r8), PARAMETER ::   &    ! R8 versions for Valid entries
       elsunc_usrstop_eval_r8 =   -10_r8, &
       elsunc_highest_eval_r8 = 12344_r8

  ! ----------------------------------
  ! Entries for main quality data flag
  ! ----------------------------------
  INTEGER (KIND=i2), PARAMETER :: &
       main_qa_missing = -1_i2, main_qa_good = 0_i2, main_qa_suspect = 1_i2, main_qa_bad = 2_i2
  INTEGER (KIND=i2), PARAMETER :: &
       main_qa_min_flag = main_qa_missing, main_qa_max_flag = main_qa_bad
  REAL    (KIND=r8), PARAMETER :: &
       main_qa_min_flag_r8 = REAL(main_qa_missing, KIND=r8), &
       main_qa_max_flag_r8 = REAL(main_qa_bad,     KIND=r8)

  ! -----------------------------------------------------------------
  ! Blank strings of various lengths.
  !
  ! Some compilers don't allow to define CHARACTER PARAMETER arrays
  ! which are initialized with field of unequal length. The following
  ! are a few padding strings that we attach to shorter entries. Not
  ! very stylish, but effective. Note that what matters is the LEN
  ! declaration - we don't need to initialize the strings with the
  ! appropriate number of blanks. Using "" is perfectly fine.
  ! -----------------------------------------------------------------
  CHARACTER (LEN=13), PARAMETER :: blank13 = ""
  CHARACTER (LEN=21), PARAMETER :: blank21 = ""
  CHARACTER (LEN=23), PARAMETER :: blank23 = ""
  CHARACTER (LEN=24), PARAMETER :: blank24 = ""
  CHARACTER (LEN=25), PARAMETER :: blank25 = ""
  CHARACTER (LEN=27), PARAMETER :: blank27 = ""
  CHARACTER (LEN=30), PARAMETER :: blank30 = ""
  

END MODULE OMSAO_parameters_module
