  ! *********************** Modification History ***************
  ! Xiong Liu, July 2003 (!xliu)
  ! 1. Increase the parameter: n_max_fitpars by maxlay + 3
  ! 2. Add four voigt parameters (left/right Gaussian width, and
  !    left/right ratios of Lorentz to Gaussian widths)
  ! *************************************************************


MODULE OMSAO_indices_module

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: maxlay, maxwin

  IMPLICIT NONE

  ! ================================================================
  ! In this MODULE we collect all indices that we have defined to
  ! generalize the match between reference spectra and their
  ! associated fitting parameters. Some of the indices overlap with
  ! those required for unique idenfication of I/O files for the OMI
  ! PGE PCF file.
  ! ================================================================

  ! --------------------------------------------------------------------
  ! The following list represents all the reference spectra we can
  ! encounter during OMI fitting of OClO, BrO, and HCHO (the "PGEs").
  ! Not all entries are required for all three SAO PGEs, but each of 
  ! them is present in at least one PGE. 
  !
  ! The entry for the fitting input control file (ICF) is not a 
  ! reference spectrum, and is therefore assigned entry "0".
  !
  ! The entry for the air mas factor (AMF) is also not a reference 
  ! spectrum, and is therefore assigned the last entry.
  !
  ! Each index (except ICF and AMF) has an identification string  
  ! associated with it, so that it can be used for the identification
  ! of reference spectra and fitting parameters in the ICF.
  !
  !     icf:    Input control file (not a reference spectrum)
  !
  !     solar:  Solar reference (usually Kitt Peak)
  !     ring:   Ring
  !     o3_t1:  O3, first  temperature
  !     o3_t2:  O3, second temperature
  !     o3_t3:  O3, third  temperature
  !     no2_t1: NO2, first  temperature
  !     no2_t2: NO2, second temperature
  !     o2o2:   O2-O2 collision complex (usually Greenblatt)
  !     so2:    SO2
  !     bro:    BrO
  !     oclo:   OClO
  !     hcho:   HCHO
  !     comm:   Common mode
  !     us1:    First undersampling spectrum
  !     us2:    Second undersampling spectrum
  !     bro_tc: BrO total column (from OMBrO   PGE fitting)
  !     o3_tc:  BrO total column (from OMHCHO  PGE fitting)
  !     pseudo: Pseudo absorber
  !     polcor: Polarization correction
  !     o2gam:  o2 gamma
  !     h2o:    h2o
  !     ring1:  second ring effect
  !     comod1: second common mode residual
  !     fsl   : irradiance stray light
  !     rsl   : radiance stray light    
  !     noname: Not yet determined, dummy placeholder
  !
  !     amf:    Air mass factor (not a reference spectrum)
  ! ---------------------------------------------------------------------
  INTEGER, PARAMETER :: &
       icf_idx    =  0, solar_idx  =  1, ring_idx   =  2, o3_t1_idx  =  3, &
       o3_t2_idx  =  4, o3_t3_idx  =  5, no2_t1_idx =  6, no2_t2_idx =  7, &
       o2o2_idx   =  8, so2_idx    =  9, bro_idx    = 10, oclo_idx   = 11, &
       hcho_idx   = 12, comm_idx   = 13, us1_idx    = 14, us2_idx    = 15, &
       bro_tc_idx = 16, o3_tc_idx  = 17, pabs_idx   = 18, polcor_idx = 19, &
       o2gam_idx  = 20, h2o_idx    = 21, ring1_idx  = 22, com1_idx   = 23, &
       so2v_idx   = 24, bro2_idx   = 25, glyox_idx  = 26, io_idx     = 27, &
       vraman_idx = 28, fsl_idx    = 29, rsl_idx    = 30, noname_idx = 31, &
       amf_idx    = 32

  INTEGER           :: comfidx, cm1fidx, comvidx, cm1vidx

  ! ----------------------------------------------------------
  ! The minimum and maximum indices of reference spectra (rs).
  ! ----------------------------------------------------------
  INTEGER, PARAMETER :: min_rs_idx = solar_idx, max_rs_idx = noname_idx

  ! --------------------------------------------------------------
  ! The identification strings associated with the fitting indices
  ! --------------------------------------------------------------
  CHARACTER (LEN=6), DIMENSION (min_rs_idx:max_rs_idx), PARAMETER :: &
       refspec_strings = (/ &
       'solar ', 'ring  ', 'o3_t1 ', 'o3_t2 ', 'o3_t3 ', 'no2_t1', 'no2_t2', &
       'o2o2  ', 'so2   ', 'bro   ', 'oclo  ', 'hcho  ', 'commod', 'usamp1', &
       'usamp2', 'bro_tc', 'o3_tc ', 'pseudo', 'polcor', 'o2gam ', 'h2o   ', &
       'ring1 ', 'comod1', 'so2v',   'bro2',   'glyox ', 'io    ', 'vraman', &
       'fsl',    'rsl',    'noname'  /)


  ! ==============================================
  ! Now we define the specific fitting parameters.
  ! ==============================================

  ! ----------------------------------------------------------------------------
  ! (1) Solar and Earthshine radiance fitting indices and identification strings
  ! ----------------------------------------------------------------------------
  INTEGER, PARAMETER :: solfit_idx = 1, radfit_idx = 2
  CHARACTER (LEN=12), DIMENSION (radfit_idx), PARAMETER :: &
       solradfit_str = (/ 'solar_fit   ', 'radiance_fit' /)

  ! ---------------------------------------------------------------------
  ! (2) Particular fitting parameters: Solar fit and radiance calibration
  ! ---------------------------------------------------------------------
  !     bl0: baseline, 0 order
  !     bl1: baseline, 1 order
  !     bl2: baseline, 2 order
  !     bl3: baseline, 3 order
  !     sc0: scaling,  0 order
  !     sc1: scaling,  1 order
  !     sc2: scaling,  2 order
  !     sc3: scaling,  3 order
  !     sin: solar intensity
  !     hwe: slit width at 1/e
  !     asy: slit function asymmetry
  !     shi: spectral shift
  !     squ: spectral squeeze
  !     vgl: Left ratio of Lorentz to Gaussian width
  !     vgr: Right ratio of Lorentz to Gaussian Width
  !     hwl: Left gaussian width at 1/e
  !     hwr: right gaussian width at 1/e
  ! --------------------------
  INTEGER, PARAMETER :: &
       bl0_idx =  1, bl1_idx =  2, bl2_idx =  3, bl3_idx =  4, sc0_idx =  5, sc1_idx =  6, &
       sc2_idx =  7, sc3_idx =  8, sin_idx =  9,  hwe_idx = 10, asy_idx = 11, &
       shi_idx = 12, squ_idx = 13, vgl_idx = 14, vgr_idx = 15, hwl_idx =  16, &
       hwr_idx = 17,  max_calfit_idx = hwr_idx

  CHARACTER (LEN=3), DIMENSION (max_calfit_idx), PARAMETER :: calfit_strings = (/ &
       'bl0', 'bl1', 'bl2', 'bl3', 'sc0', 'sc1', 'sc2', 'sc3', 'sin', 'hwe', 'asy', &
       'shi', 'squ', 'vgl', 'vgr', 'hwl', 'hwr'/)

  ! ------------------------------------------------------------
  ! (3) Particular fitting parameters: Radiance spectral fitting
  ! ------------------------------------------------------------
  !     ad1: first added contribution
  !     lbe: Lambert-Beer terms
  !     ad2: second added contribution
  !     mns: minimum radiance fitting sub-index
  !     mxs: maximum radiance fitting sub-index
  ! -------------------------------------------
  INTEGER, PARAMETER :: ad1_idx = 1, lbe_idx = 2, ad2_idx = 3
  INTEGER, PARAMETER :: mns_idx = ad1_idx, mxs_idx = ad2_idx
  CHARACTER (LEN=3), DIMENSION (mxs_idx), PARAMETER :: &
       radfit_strings = (/ 'ad1', 'lbe', 'ad2' /)

  ! -----------------------------------------------------------------
  ! Total number of (radiance) fitting parameters: The MAX_CALFIT_IDX
  ! calibration parameters, plus MAX_RS_IDX*AD2_IDX parameters
  ! associated with external or online computed reference spectra.
  ! -----------------------------------------------------------------
  ! xliu: add  max_rs_idx + maxlay*2 + maxother terms
  ! 1. shift parameter for each species with reference spectrum
  ! 2. ozone and temperature profiles
  ! 3. 6 terms for albedo (channel 1, 2, 3, 4) plus 1st and 2nd albedo
  ! 4. shift parameter for o3 cross section (4), sol/rad slit diff (4)  
  !    sol/rad pos diff (4)
  !    paramteres in 5 can be 5 x 5 terms
  ! 5. Ring effect, zero order (all windows)
  ! 6. Ring effect (1, 2, 3 for channel 2 only)
  ! 7. internal scattering for irradiance
  !    internal scattering for radiance
  !    degradation correction 1
  !    degradation corrrection 2
  !    degradation corrrection 3

  INTEGER, PARAMETER :: maxalb = 12, maxwfc = 12, maxoth = 4, maxgrp = 7, maxcldaer=8
  INTEGER, PARAMETER :: n_max_fitpars = max_calfit_idx + mxs_idx * max_rs_idx &
       + max_rs_idx + maxlay * 2  + maxalb + maxwfc + maxcldaer + (maxgrp * maxoth) * maxwin
  INTEGER, PARAMETER :: shift_offset = max_calfit_idx + max_rs_idx * mxs_idx

  
  ! --------------------------------
  ! ELSUNC fitting constraints index
  ! --------------------------------
  INTEGER, PARAMETER :: elsunc_unconstrained = 0, elsunc_same_lower = 1,&
       elsunc_userdef = 2
  
  ! ------------------------------------------------------
  ! Indices for fitting wavelengths, spectrum, and weights
  ! ------------------------------------------------------
  INTEGER, PARAMETER :: wvl_idx = 1, spc_idx = 2, sig_idx = 3

  ! -----------------------------------------
  ! Indices for angles and AMFs in AMF tables
  ! -----------------------------------------
  INTEGER, PARAMETER :: amftab_ang_idx = 1, amftab_amf_idx = 2
  INTEGER, PARAMETER :: n_amftab_ang_max = 16, n_amftab_dim_max = 2


  ! =================================================================
  !
  ! This module defines a range of indices for the OMI SAO PGE fitting
  ! processes. These are mostly associated with file unit numbers as
  ! they are required by the PGE Process Control File (PCF). But they
  ! also include the number of SAO PGEs, the official PGE numbers for
  ! the molecules, etc.
  !
  ! =================================================================

  ! ------------------------------------------------------------------------
  ! First define the number of SAO PGEs. We are fitting three molecules, so:
  ! ------------------------------------------------------------------------
  INTEGER (KIND=i4), PARAMETER :: n_sao_pge = 4
  
  ! -----------------------------------------------------------------
  ! The operational OMI environment assigns reference numbers to each
  ! PGE. Here we adopt the official setting for future indexing.
  ! -----------------------------------------------------------------
  INTEGER (KIND=i4), PARAMETER :: &
       pge_oclo_idx = 10, pge_bro_idx = 11, pge_hcho_idx = 12, pge_o3_idx = 13, &
       pge_no2_idx = 14, pge_so2_idx = 15, pge_gly_idx = 16, pge_no2d_idx = 17
  INTEGER (KIND=i4), PARAMETER :: sao_pge_min_idx = pge_oclo_idx, sao_pge_max_idx = pge_no2d_idx

  ! ------------------------------------------------
  ! Identifications strings and indices for SAO PGEs
  ! ------------------------------------------------
  CHARACTER (LEN=8), DIMENSION (sao_pge_min_idx:sao_pge_max_idx), PARAMETER :: &
       sao_pge_names      = (/ 'OMOCLO  ', 'OMBRO   ', 'OMHCHO  ', 'OMSAO3  ', &
       'OMSAONO2', 'OMSAOSO2', 'OMCHOCHO', 'OMNO2D  ' /)
  CHARACTER (LEN=6), DIMENSION (sao_pge_min_idx:sao_pge_max_idx), PARAMETER :: &
       sao_molecule_names = (/ 'OClO  ',   'BrO   ',   'HCHO  ',   'O3    ',   &
       'NO2   ',   'SO2   ',   'CHOCHO',   'NO2   '   /)

  ! ----------------------------------------------------------------------
  ! What follows are the input file unit numbers required for uniquely 
  ! identifying files in the PCF. Remember that each PGE has a separate
  ! OMI PGE number, and thus we have to define three plus one sets of 
  ! file units. The fourth set is that for O3, which requires a full PGE
  ! set of inputs without having the benefit of being an official SAO PGE.
  ! The easiest way to do this is to assign two full sets of input numbers
  ! for the HCHO PGE (#12).
  !
  ! Note on RESHAPE: 
  !   This array functions reshapes a 1-dim array into a multi-dimensional
  !   one. The first (new) array dimension varies fastest, followed by the
  !   second, and so on.
  ! ----------------------------------------------------------------------

  ! -----------------------------------------------
  ! LUNs. Some registered, some not, some temporary
  ! -----------------------------------------------
  INTEGER, PARAMETER :: pge_mol_id             = 700000   ! PGE molecule ID; registered
  INTEGER, PARAMETER :: pge_l1b_radiance_lun   = 749000   ! PGE L1B randiance  file; not yet registered
  INTEGER, PARAMETER :: pge_l1b_irradiance_lun = 749001   ! PGE L1B irradiance file; not yet registered
  INTEGER, PARAMETER :: verbosity_lun          = 200100   ! PGE verbosity threshold
  INTEGER, PARAMETER :: orbitnumber_lun        = 200200   ! LUN for orbit number in PCF
  INTEGER, PARAMETER :: granule_s_lun          =  10258   ! LUN for GranuleStartTime in PCF
  INTEGER, PARAMETER :: granule_e_lun          =  10259   ! LUN for GranuleStartTime in PCF
  INTEGER, DIMENSION (sao_pge_min_idx:sao_pge_max_idx), PARAMETER :: &
       mcf_lun  = (/ 710001, 711001, 712001, 712002, 712003, 712004, 712005, 712006 /)  ! MCF file; not yet registered


  ! -------------------------------------------------------
  ! GOME data LUNs - temporary until OMI data are available
  ! -------------------------------------------------------
  INTEGER, PARAMETER :: gome_data_input_lun    = 749003   ! GOME data input; temporary
  INTEGER, PARAMETER :: gome_data_output_lun   = 749004   ! GOME output;     temporary


  ! ----------------------------------------------------------------------
  ! * Input file LUNs for reference spectra and algorithm control file.
  !   First (ICF) and last (AMF) lun are for static, non-reference spectra
  !   input files.
  ! ----------------------------------------------------------------------
  INTEGER, DIMENSION (icf_idx:amf_idx, pge_oclo_idx:pge_hcho_idx+1), PARAMETER :: &
       pge_static_input_luns = RESHAPE ( (/ &
       710100,  710101,  710102,  710103,  710104,  710105,  710106,  710107, &
       710108,  710109,  710110,  710111,  710112,  710113,  710114,  710115, &
       710116,  710117,  710118,  710119,  710120,  710121,  710122,  710123, &
       710124,  710125,  710126,  710127,  710128,  710129,  710130,  710131, &
       710132,                                                                &
       711100,  711101,  711102,  711103,  711104,  711105,  711106,  711107, &
       711108,  711109,  711110,  711111,  711112,  711113,  711114,  711115, &
       711116,  711117,  711118,  711119,  711120,  711121,  711122,  711123, &
       711124,  711125,  711126,  711127,  711128,  711129,  711130,  711131, &
       711132,                                                                &
       712100,  712101,  712102,  712103,  712104,  712105,  712106,  712107, &
       712108,  712109,  712110,  712111,  712112,  712113,  712114,  712115, &
       712116,  712117,  712118,  712119,  712120,  712121,  712122,  712123, &
       712124,  712125,  712126,  712127,  712128,  712129,  712130,  712131, &
       712132,                                                                &
       712200,  712201,  712202,  712203,  712204,  712205,  712206,  712207, &
       712208,  712209,  712210,  712211,  712212,  712213,  712214,  712215, &
       712216,  712217,  712218,  712219,  712220,  712221,  712222,  712223, &
       712224,  712225,  712226,  712227,  712228,  712229,  712230,  712231, &
       712232/),&
        (/ amf_idx+1, n_sao_pge /) )

  ! ------------------------
  ! * LUNs for PGE L2 Output
  ! ------------------------
  INTEGER, DIMENSION (pge_oclo_idx:pge_hcho_idx), PARAMETER :: &
       pge_l2_output_luns = (/ 710999, 711999, 712999 /)

  ! -------------------------------------------------------
  ! Strings to search for in the fitting control input file
  ! -------------------------------------------------------
  CHARACTER (LEN=18), PARAMETER :: molline_str = 'Molecule(s) to fit'
  CHARACTER (LEN=26), PARAMETER :: genline_str = 'General fitting parameters'
  CHARACTER (LEN=28), PARAMETER :: socline_str = 'Solar calibration parameters'
  CHARACTER (LEN=31), PARAMETER :: racline_str = 'Radiance calibration parameters'
  CHARACTER (LEN=27), PARAMETER :: rafline_str = 'Radiance fitting parameters'
  CHARACTER (LEN=24), PARAMETER :: rspline_str = 'Input reference spectra'
  CHARACTER (LEN=24), PARAMETER :: iofline_str = 'Input/output data files'

  ! --------------------------------------------------------------
  ! End-Of-Input string; used to terminate a list of string inputs
  ! --------------------------------------------------------------
  CHARACTER (LEN=12), PARAMETER :: eoi_str = 'end_of_input'
  CHARACTER (LEN= 3), PARAMETER :: eoi3str = 'eoi'
  
  ! -----------------------------------------------------------------------
  ! Input/output units for reading/writing data, ancillary and output files
  ! -----------------------------------------------------------------------
  INTEGER, PARAMETER :: &
       l1b_input_unit = 11, l2_output_unit = 12, static_input_unit = 13

  ! ----------------
  ! Metadata indices
  ! ----------------
  INTEGER, PARAMETER :: md_inventory_idx = 2, md_archive_idx = 3

END MODULE OMSAO_indices_module


