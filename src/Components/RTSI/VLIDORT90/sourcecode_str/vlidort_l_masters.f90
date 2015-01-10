! ###############################################################
! #                                                             #
! #                    THE VLIDORT  MODEL                       #
! #                                                             #
! #  Vectorized LInearized Discrete Ordinate Radiative Transfer #
! #  -          --         -        -        -         -        #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Author :      Robert. J. D. Spurr                          #
! #                                                             #
! #  Address :     RT Solutions, inc.                           #
! #                9 Channing Street                            #
! #                Cambridge, MA 02138, USA                     #
! #                Tel: (617) 492 1183                          #
! #                                                             #
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  Versions     :   2.0, 2.2, 2.3, 2.4, 2.4R, 2.4RT, 2.4RTC,  #
! #                   2.5, 2.6                                  #
! #  Release Date :   December 2005  (2.0)                      #
! #  Release Date :   March 2007     (2.2)                      #
! #  Release Date :   October 2007   (2.3)                      #
! #  Release Date :   December 2008  (2.4)                      #
! #  Release Date :   April 2009     (2.4R)                     #
! #  Release Date :   July 2009      (2.4RT)                    #
! #  Release Date :   October 2010   (2.4RTC)                   #
! #  Release Date :   March 2011     (2.5)                      #
! #  Release Date :   May 2012       (2.6)                      #
! #                                                             #
! #       NEW: TOTAL COLUMN JACOBIANS         (2.4)             #
! #       NEW: BPDF Land-surface KERNELS      (2.4R)            #
! #       NEW: Thermal Emission Treatment     (2.4RT)           #
! #       Consolidated BRDF treatment         (2.4RTC)          #
! #       f77/f90 Release                     (2.5)             #
! #       External SS / New I/O Structures    (2.6)             #
! #                                                             #
! ###############################################################

!    #####################################################
!    #                                                   #
!    #   This Version of VLIDORT comes with a GNU-style  #
!    #   license. Please read the license carefully.     #
!    #                                                   #
!    #####################################################

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #            VLIDORT_L_MASTER (master)                        #
! #            VLIDORT_L_FOURIER (master)                       #
! #                                                             #
! ###############################################################


      MODULE vlidort_l_masters

      PRIVATE
      PUBLIC :: VLIDORT_L_MASTER

      CONTAINS

      SUBROUTINE VLIDORT_L_MASTER ( &
        VLIDORT_FixIn, &
        VLIDORT_ModIn, &
        VLIDORT_Sup, &
        VLIDORT_Out, &
        VLIDORT_LinFixIn, &
        VLIDORT_LinModIn, &
        VLIDORT_LinSup, &
        VLIDORT_LinOut )

      USE VLIDORT_PARS

      USE VLIDORT_Inputs_def
      USE VLIDORT_Sup_InOut_def
      USE VLIDORT_Outputs_def
      USE VLIDORT_LinInputs_def
      USE VLIDORT_LinSup_InOut_def
      USE VLIDORT_LinOutputs_def

      USE VLIDORT_INPUTS
      USE VLIDORT_MISCSETUPS_MODULE
      USE VLIDORT_CORRECTIONS
      USE VLIDORT_CONVERGE_MODULE
      USE VLIDORT_WRITEMODULES

      USE VLIDORT_L_INPUTS
      USE VLIDORT_L_MISCSETUPS_MODULE
      USE VLIDORT_L_CORRECTIONS
      USE VLIDORT_L_CONVERGE_MODULE
      USE VLIDORT_L_WRITEMODULES

      IMPLICIT NONE

!  VLIDORT input structures

      TYPE(VLIDORT_Fixed_Inputs), INTENT (IN)       :: &
        VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs), INTENT (INOUT) :: &
        VLIDORT_ModIn

!  VLIDORT supplements structure

      TYPE(VLIDORT_Sup_InOut), INTENT (INOUT)       :: &
        VLIDORT_Sup

!  VLIDORT output structure

      TYPE(VLIDORT_Outputs), INTENT (OUT)           :: &
        VLIDORT_Out

!  VLIDORT linearized input structures

      TYPE(VLIDORT_Fixed_LinInputs), INTENT (IN)       :: &
        VLIDORT_LinFixIn
      TYPE(VLIDORT_Modified_LinInputs), INTENT (INOUT) :: &
        VLIDORT_LinModIn

!  VLIDORT linearized supplements structure

      TYPE(VLIDORT_LinSup_InOut), INTENT (INOUT)       :: &
        VLIDORT_LinSup

!  VLIDORT linearized output structure

      TYPE(VLIDORT_LinOutputs), INTENT (OUT)           :: &
        VLIDORT_LinOut

!  -----------------------
!  Standard Inputs - Fixed
!  -----------------------

      LOGICAL ::            DO_FULLRAD_MODE
      LOGICAL ::            DO_SSCORR_TRUNCATION
      LOGICAL ::            DO_SSFULL
      LOGICAL ::            DO_THERMAL_EMISSION
      LOGICAL ::            DO_SURFACE_EMISSION
      LOGICAL ::            DO_PLANE_PARALLEL
      LOGICAL ::            DO_UPWELLING
      LOGICAL ::            DO_DNWELLING
      LOGICAL ::            DO_QUAD_OUTPUT
      LOGICAL ::            DO_TOA_CONTRIBS
      LOGICAL ::            DO_LAMBERTIAN_SURFACE
      LOGICAL ::            DO_SPECIALIST_OPTION_1
      LOGICAL ::            DO_SPECIALIST_OPTION_2
      LOGICAL ::            DO_SPECIALIST_OPTION_3

!  New 15 March 2012
      LOGICAL ::            DO_SS_EXTERNAL

!  New 17 May 2012
      LOGICAL ::            DO_SURFACE_LEAVING
      LOGICAL ::            DO_SL_ISOTROPIC

      INTEGER ::            NSTOKES
      INTEGER ::            NSTREAMS
      INTEGER ::            NLAYERS
      INTEGER ::            NFINELAYERS
      INTEGER ::            N_THERMAL_COEFFS
      DOUBLE PRECISION ::   VLIDORT_ACCURACY
      INTEGER ::            NLAYERS_NOMS
      INTEGER ::            NLAYERS_CUTOFF

      DOUBLE PRECISION ::   SZANGLES ( MAX_SZANGLES )

      INTEGER ::            N_USER_VZANGLES
      INTEGER ::            N_USER_LEVELS
      DOUBLE PRECISION ::   USER_LEVELS ( MAX_USER_LEVELS )

      DOUBLE PRECISION ::   HEIGHT_GRID ( 0:MAXLAYERS )
      DOUBLE PRECISION ::   PRESSURE_GRID ( 0:MAXLAYERS )
      DOUBLE PRECISION ::   TEMPERATURE_GRID ( 0:MAXLAYERS )
      INTEGER ::            FINEGRID ( MAXLAYERS )
      DOUBLE PRECISION ::   RFINDEX_PARAMETER

      DOUBLE PRECISION ::   DELTAU_VERT_INPUT ( MAXLAYERS )
      DOUBLE PRECISION ::   GREEKMAT_TOTAL_INPUT &
          ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION ::   THERMAL_BB_INPUT ( 0:MAXLAYERS )
      DOUBLE PRECISION ::   LAMBERTIAN_ALBEDO
      DOUBLE PRECISION ::   SURFBB
      DOUBLE PRECISION ::   LTE_DELTAU_VERT_INPUT ( 2, MAXLAYERS )
      DOUBLE PRECISION ::   LTE_THERMAL_BB_INPUT ( 0:MAXLAYERS )

      LOGICAL ::            DO_DEBUG_WRITE
      LOGICAL ::            DO_WRITE_INPUT
      LOGICAL ::            DO_WRITE_SCENARIO
      LOGICAL ::            DO_WRITE_FOURIER
      LOGICAL ::            DO_WRITE_RESULTS
      CHARACTER (LEN=60) :: INPUT_WRITE_FILENAME
      CHARACTER (LEN=60) :: SCENARIO_WRITE_FILENAME
      CHARACTER (LEN=60) :: FOURIER_WRITE_FILENAME
      CHARACTER (LEN=60) :: RESULTS_WRITE_FILENAME

!  --------------------------
!  Standard Inputs - Variable
!  --------------------------

      LOGICAL ::          DO_SSCORR_NADIR
      LOGICAL ::          DO_SSCORR_OUTGOING
      LOGICAL ::          DO_DOUBLE_CONVTEST
      LOGICAL ::          DO_SOLAR_SOURCES
      LOGICAL ::          DO_REFRACTIVE_GEOMETRY
      LOGICAL ::          DO_CHAPMAN_FUNCTION
      LOGICAL ::          DO_RAYLEIGH_ONLY
      LOGICAL ::          DO_DELTAM_SCALING
      LOGICAL ::          DO_SOLUTION_SAVING
      LOGICAL ::          DO_BVP_TELESCOPING
      LOGICAL ::          DO_USER_VZANGLES
      LOGICAL ::          DO_ADDITIONAL_MVOUT
      LOGICAL ::          DO_MVOUT_ONLY
      LOGICAL ::          DO_THERMAL_TRANSONLY

      INTEGER ::          NGREEK_MOMENTS_INPUT

      DOUBLE PRECISION :: FLUX_FACTOR
      INTEGER ::          N_SZANGLES

      INTEGER ::          N_USER_RELAZMS
      DOUBLE PRECISION :: USER_RELAZMS  ( MAX_USER_RELAZMS )
      DOUBLE PRECISION :: USER_VZANGLES ( MAX_USER_VZANGLES )
      DOUBLE PRECISION :: GEOMETRY_SPECHEIGHT

      DOUBLE PRECISION :: EARTH_RADIUS

      DOUBLE PRECISION :: OMEGA_TOTAL_INPUT ( MAXLAYERS )

!  -----------------------
!  Standard Supplement I/O
!  -----------------------

!  BRDF Inputs
!  -----------

      DOUBLE PRECISION ::   EXACTDB_BRDFUNC &
          ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )
      DOUBLE PRECISION ::   BRDF_F_0 &
          ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION ::   BRDF_F &
          ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION ::   USER_BRDF_F_0 &
          ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION ::   USER_BRDF_F &
          ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )
      DOUBLE PRECISION ::   EMISSIVITY ( MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION ::   USER_EMISSIVITY &
          ( MAXSTOKES, MAX_USER_STREAMS )

!  SS and DB I/O
!  -------------

!  SS Intensity Results at all angles and optical depths

      DOUBLE PRECISION ::  STOKES_SS &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )

!  DB Intensity Results at all angles and optical depths

      DOUBLE PRECISION ::  STOKES_DB &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  Contribution function (TOA Upwelling only)

!  SS component of Diffuse Field

!      DOUBLE PRECISION ::  SS_CONTRIBS &
!          ( MAX_GEOMETRIES, MAXSTOKES, MAXLAYERS )

!  New Surface-Leaving Inputs, 17 May 12
!  -------------------------------------

      DOUBLE PRECISION ::   SLTERM_ISOTROPIC ( MAXSTOKES )
      DOUBLE PRECISION ::   SLTERM_USERANGLES &
          ( MAXSTOKES, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS  )

      DOUBLE PRECISION ::   SLTERM_F_0 &
          ( 0:MAXMOMENTS, MAXSTOKES, MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION ::   USER_SLTERM_F_0 &
          ( 0:MAXMOMENTS, MAXSTOKES, MAX_USER_STREAMS, MAXBEAMS )

!  ----------------
!  Standard Outputs
!  ----------------

!  Fourier values

!      DOUBLE PRECISION :: STOKES_F &
!          ( MAX_USER_LEVELS, MAX_USER_VZANGLES, &
!            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  Fourier-summed values

      DOUBLE PRECISION :: STOKES &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )

!  Results for mean-value output
!  -----------------------------

!  Complete Actinic and Regular Fluxes (including Direct terms)

      DOUBLE PRECISION :: MEAN_STOKES &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: FLUX_STOKES &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  Direct Fluxes only

      DOUBLE PRECISION :: MEAN_DIRECT &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES )
      DOUBLE PRECISION :: FLUX_DIRECT &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES )

!  Contribution functions (TOA Upwelling only)
!  -------------------------------------------

!  Fourier component of Diffuse Field

!      DOUBLE PRECISION ::  MS_CONTRIBS_F &
!          ( MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAXLAYERS  )

!  Fourier-summed values

!      DOUBLE PRECISION ::  CONTRIBS &
!          ( MAX_GEOMETRIES, MAXSTOKES, MAXLAYERS )

!  Ancillary Output
!  ----------------

!  Fourier numbers used

      INTEGER ::          FOURIER_SAVED ( MAX_SZANGLES )

!  Number of geometries computed

      INTEGER ::          N_GEOMETRIES

!  Offsets for indexing geometries

      INTEGER ::          SZA_OFFSETS ( MAX_SZANGLES )
      INTEGER ::          VZA_OFFSETS ( MAX_SZANGLES, MAX_USER_VZANGLES )

!  -------------------------
!  Linearized Inputs - Fixed
!  -------------------------

      LOGICAL ::            DO_SIMULATION_ONLY
      LOGICAL ::            DO_LINEARIZATION

      LOGICAL ::            DO_LTE_LINEARIZATION
      LOGICAL ::            DO_SURFBB_LINEARIZATION

      LOGICAL ::            LAYER_VARY_FLAG ( MAXLAYERS )
      INTEGER ::            LAYER_VARY_NUMBER ( MAXLAYERS )

      INTEGER ::            N_TOTALCOLUMN_WFS
      INTEGER ::            N_TOTALPROFILE_WFS
      INTEGER ::            N_SURFACE_WFS

      CHARACTER (LEN=31) :: COLUMNWF_NAMES ( MAX_ATMOSWFS )
      CHARACTER (LEN=31) :: PROFILEWF_NAMES ( MAX_ATMOSWFS )

      DOUBLE PRECISION ::   L_DELTAU_VERT_INPUT &
          ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION ::   L_OMEGA_TOTAL_INPUT &
          ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION ::   L_GREEKMAT_TOTAL_INPUT &
          ( MAX_ATMOSWFS, 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )

!  ----------------------------
!  Linearized Inputs - Variable
!  ----------------------------

      LOGICAL ::         DO_ATMOS_LINEARIZATION
      LOGICAL ::         DO_COLUMN_LINEARIZATION
      LOGICAL ::         DO_PROFILE_LINEARIZATION
      LOGICAL ::         DO_SURFACE_LINEARIZATION

!  -------------------------
!  Linearized Supplement I/O
!  -------------------------

!  BRDF Inputs
!  -----------

      DOUBLE PRECISION ::   LS_EXACTDB_BRDFUNC &
          ( MAX_SURFACEWFS, MAXSTOKES_SQ, &
            MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )
      DOUBLE PRECISION ::   LS_BRDF_F_0 &
          ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTOKES_SQ, &
            MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION ::   LS_BRDF_F &
          ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTOKES_SQ, &
            MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION ::   LS_USER_BRDF_F_0 &
          ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTOKES_SQ, &
            MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION ::   LS_USER_BRDF_F &
          ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTOKES_SQ, &
            MAX_USER_STREAMS, MAXSTREAMS )

      DOUBLE PRECISION ::   LS_EMISSIVITY &
          ( MAX_SURFACEWFS, MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION ::   LS_USER_EMISSIVITY &
          ( MAX_SURFACEWFS, MAXSTOKES, MAX_USER_STREAMS )

!  SS and DB I/O
!  -------------

!  SS Profile Jacobian Results at all angles and optical depths

      DOUBLE PRECISION :: PROFILEWF_SS &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
            MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )

!  DB Profile Jacobian Results at all angles and optical depths

      DOUBLE PRECISION :: PROFILEWF_DB &
          ( MAX_ATMOSWFS, MAXLAYERS, &
            MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  SS Column Jacobian Results at all angles and optical depths

      DOUBLE PRECISION :: COLUMNWF_SS &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, &
            MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )

!  DB Column Jacobian Results at all angles and optical depths

      DOUBLE PRECISION :: COLUMNWF_DB &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  DB Surface Jacobian Results at all angles and optical depths

      DOUBLE PRECISION :: SURFACEWF_DB &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  New Surface-Leaving Inputs, 17 May 12
!  -------------------------------------

!  Empty slot

!  ------------------
!  Linearized Outputs
!  ------------------

! Fourier values

!      DOUBLE PRECISION :: ATMOSWF_F &
!         ( MAX_ATMOSWFS, 0:MAXLAYERS, MAX_USER_LEVELS, &
!           MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!      DOUBLE PRECISION :: SURFACEWF_F &
!         ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_VZANGLES, &
!           MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  Fourier-summed values

      DOUBLE PRECISION :: PROFILEWF &
         ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
           MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )

      DOUBLE PRECISION :: COLUMNWF &
         ( MAX_ATMOSWFS, MAX_USER_LEVELS, &
           MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )

      DOUBLE PRECISION :: SURFACEWF &
         ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
           MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )

!  mean intensity (actinic flux) weighting functions. Total Fluxes.
!   ( atmospheric WFs have separate Direct Beam WFs. Introduced 09.24.09

      DOUBLE PRECISION :: MINT_ATMOSWF &
          ( MAX_ATMOSWFS, 0:MAXLAYERS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

      DOUBLE PRECISION :: MINT_ATMOSWF_DIRECT &
          ( MAX_ATMOSWFS, 0:MAXLAYERS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES )

      DOUBLE PRECISION :: MINT_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  flux weighting functions. Total Fluxes.
!   ( atmospheric WFs have separate Direct Beam WFs. Introduced 09.24.09

      DOUBLE PRECISION :: FLUX_ATMOSWF &
          ( MAX_ATMOSWFS, 0:MAXLAYERS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

      DOUBLE PRECISION :: FLUX_ATMOSWF_DIRECT &
          ( MAX_ATMOSWFS, 0:MAXLAYERS, MAX_USER_LEVELS, &
           MAX_SZANGLES, MAXSTOKES )

      DOUBLE PRECISION :: FLUX_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  RT Solutions Inc. RJD Spurr.  11 September 2009.
!    LTE linearization: Introduction of T-Jacobians for BB functions
!   Only works with pure thermal emission (no scattering)
!  Introduced for the GEOCAPE study

      DOUBLE PRECISION :: LTE_ATMOSWF &
          ( 0:MAXLAYERS, MAX_USER_LEVELS, &
            MAX_USER_VZANGLES, MAX_DIRECTIONS )

!  Error handling

      INTEGER ::             STATUS_INPUTCHECK
      INTEGER ::             NCHECKMESSAGES
      CHARACTER (LEN=120) :: CHECKMESSAGES (0:MAX_MESSAGES)
      CHARACTER (LEN=120) :: ACTIONS (0:MAX_MESSAGES)
      INTEGER ::             STATUS_CALCULATION
      CHARACTER (LEN=120) :: MESSAGE, TRACE_1, TRACE_2, TRACE_3

!  Local variables
!  ---------------

      LOGICAL ::          LOCAL_DO_NO_AZIMUTH
      LOGICAL ::          SAVE_DO_NO_AZIMUTH
      LOGICAL ::          LOCAL_ITERATION
      LOGICAL ::          DO_INCLUDE_SURFACE
      LOGICAL ::          DO_INCLUDE_SURFEMISS
      LOGICAL ::          DO_INCLUDE_THERMEMISS

      INTEGER ::          FOURIER_COMPONENT
      INTEGER ::          N_FOURIER_COMPONENTS

      INTEGER ::          UM, UA, TESTCONV, L
      INTEGER ::          LOCAL_N_USERAZM, STATUS_SUB
      INTEGER ::          IUNIT, SUNIT, RUNIT

      DOUBLE PRECISION :: AZM_ARGUMENT, DFC
      DOUBLE PRECISION :: AZMFAC &
          ( MAX_USER_STREAMS, MAXBEAMS, MAX_USER_RELAZMS, MAXSTOKES )

      INTEGER ::          IBEAM_COUNT, IBEAM, IB, NSOURCES
      LOGICAL ::          BEAM_ITERATION ( MAXBEAMS )
      INTEGER ::          BEAM_TESTCONV  ( MAXBEAMS )

      LOGICAL ::          ADJUST_SURFACE

      DOUBLE PRECISION :: SS_FLUX_MULTIPLIER, MODIFIED_ERADIUS

!  Helper variables

      LOGICAL ::          DO_ALL_FOURIER
      LOGICAL ::          DO_DIRECT_BEAM
      LOGICAL ::          DO_CLASSICAL_SOLUTION
      LOGICAL ::          DO_DBCORRECTION
      DOUBLE PRECISION :: FLUXVEC ( MAXSTOKES )
      DOUBLE PRECISION :: COS_SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION :: SIN_SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION :: SUN_SZA_COSINES ( MAXLAYERS, MAX_SZANGLES )
      LOGICAL ::          DO_MULTIBEAM ( MAXBEAMS, 0:MAXFOURIER )
      DOUBLE PRECISION :: QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION :: QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION :: QUAD_STRMWTS ( MAXSTREAMS )
      DOUBLE PRECISION :: QUAD_HALFWTS ( MAXSTREAMS )
      DOUBLE PRECISION :: QUAD_SINES   ( MAXSTREAMS )
      DOUBLE PRECISION :: QUAD_ANGLES  ( MAXSTREAMS )
      LOGICAL ::          DO_MSMODE_VLIDORT
      LOGICAL ::          DO_MSMODE_THERMAL
      LOGICAL ::          DO_NO_AZIMUTH
      INTEGER ::          NMOMENTS
      INTEGER ::          NSTREAMS_2
      INTEGER ::          NTOTAL
      INTEGER ::          N_SUBDIAG
      INTEGER ::          N_SUPDIAG
      INTEGER ::          NSTOKES_SQ
      INTEGER ::          NSTKS_NSTRMS
      INTEGER ::          NSTKS_NSTRMS_2
      INTEGER ::          NBEAMS
      !INTEGER ::          N_USER_STREAMS
      INTEGER ::          NPARTICSOLS

      INTEGER ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: DMAT ( MAXSTOKES, MAXSTOKES )

      INTEGER ::          GREEKMAT_INDEX ( 6 )
      LOGICAL ::          DO_REAL_EIGENSOLVER ( 0:MAXMOMENTS, MAXLAYERS )
      LOGICAL ::          BVP_REGULAR_FLAG ( 0:MAXMOMENTS )
      INTEGER ::          LAYER_MAXMOMENTS ( MAXLAYERS )
      LOGICAL ::          DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )
      INTEGER ::          N_CONVTESTS
      INTEGER ::          N_CONV_STREAMS
      !LOGICAL ::          DO_USER_STREAMS
      INTEGER ::          LOCAL_UM_START
      INTEGER ::          N_DIRECTIONS
      INTEGER ::          WHICH_DIRECTIONS ( MAX_DIRECTIONS )
      DOUBLE PRECISION :: USER_VZANGLES_ADJUST ( MAX_USER_VZANGLES )
      DOUBLE PRECISION :: SZANGLES_ADJUST &
          ( MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS )
      DOUBLE PRECISION :: USER_RELAZMS_ADJUST &
          ( MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS )
      DOUBLE PRECISION :: USER_STREAMS ( MAX_USER_STREAMS )
      DOUBLE PRECISION :: USER_SINES   ( MAX_USER_STREAMS )
      DOUBLE PRECISION :: USER_SECANTS ( MAX_USER_STREAMS )
      INTEGER ::          N_OUT_STREAMS
      DOUBLE PRECISION :: OUT_ANGLES ( MAX_USER_STREAMS )
      LOGICAL ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER ::          UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER ::          UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )
      LOGICAL ::          DO_PARTLAYERS
      INTEGER ::          N_PARTLAYERS
      INTEGER ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      DOUBLE PRECISION :: PARTLAYERS_VALUES   ( MAX_PARTLAYERS )
      INTEGER ::          N_LAYERSOURCE_UP
      INTEGER ::          N_LAYERSOURCE_DN
      INTEGER ::          N_ALLLAYERS_UP
      INTEGER ::          N_ALLLAYERS_DN
      LOGICAL ::          STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL ::          STERM_LAYERMASK_DN ( MAXLAYERS )
      INTEGER ::          N_VIEWING
      DOUBLE PRECISION :: CHAPMAN_FACTORS &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: DFLUX ( MAXSTOKES )
      DOUBLE PRECISION :: TAUGRID_INPUT ( 0:MAXLAYERS )
      DOUBLE PRECISION :: TAUGRID ( 0:MAXLAYERS )
      DOUBLE PRECISION :: OMEGA_TOTAL ( MAXLAYERS )
      DOUBLE PRECISION :: GREEKMAT_TOTAL &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION :: SZA_LOCAL_INPUT &
          ( 0:MAXLAYERS, MAX_SZANGLES )

      DOUBLE PRECISION :: STOKES_F &
          ( MAX_USER_LEVELS, MAX_USER_VZANGLES, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

      DOUBLE PRECISION :: MS_CONTRIBS_F &
          ( MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAXLAYERS  )
      DOUBLE PRECISION :: CONTRIBS &
          ( MAX_GEOMETRIES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION :: SS_CONTRIBS &
          ( MAX_GEOMETRIES, MAXSTOKES, MAXLAYERS )

      DOUBLE PRECISION :: ATMOSWF_F &
          ( MAX_ATMOSWFS, 0:MAXLAYERS, MAX_USER_LEVELS, &
            MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: SURFACEWF_F &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_VZANGLES, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  Local error handling

      LOGICAL ::          FAIL

!  Test variables

      LOGICAL ::          DO_FDTEST=.FALSE.
      LOGICAL ::          DO_DEBUG_INPUT=.FALSE.

!  ====================================
!  BEGIN COPY INPUTS TO LOCAL VARIABLES
!  ====================================

!  Fixed Boolean inputs

      DO_FULLRAD_MODE        = VLIDORT_FixIn%Bool%TS_DO_FULLRAD_MODE
      DO_SSCORR_TRUNCATION   = VLIDORT_FixIn%Bool%TS_DO_SSCORR_TRUNCATION
      DO_SSFULL              = VLIDORT_FixIn%Bool%TS_DO_SSFULL
      DO_THERMAL_EMISSION    = VLIDORT_FixIn%Bool%TS_DO_THERMAL_EMISSION
      DO_SURFACE_EMISSION    = VLIDORT_FixIn%Bool%TS_DO_SURFACE_EMISSION
      DO_PLANE_PARALLEL      = VLIDORT_FixIn%Bool%TS_DO_PLANE_PARALLEL
      !DO_BRDF_SURFACE        = VLIDORT_FixIn%Bool%TS_DO_BRDF_SURFACE
      DO_UPWELLING           = VLIDORT_FixIn%Bool%TS_DO_UPWELLING
      DO_DNWELLING           = VLIDORT_FixIn%Bool%TS_DO_DNWELLING
      DO_QUAD_OUTPUT         = VLIDORT_FixIn%Bool%TS_DO_QUAD_OUTPUT
      DO_TOA_CONTRIBS        = VLIDORT_FixIn%Bool%TS_DO_TOA_CONTRIBS
      DO_LAMBERTIAN_SURFACE  = VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE
      DO_SPECIALIST_OPTION_1 = VLIDORT_FixIn%Bool%TS_DO_SPECIALIST_OPTION_1
      DO_SPECIALIST_OPTION_2 = VLIDORT_FixIn%Bool%TS_DO_SPECIALIST_OPTION_2
      DO_SPECIALIST_OPTION_3 = VLIDORT_FixIn%Bool%TS_DO_SPECIALIST_OPTION_3

!  New 15 March 2012
      DO_SS_EXTERNAL         = VLIDORT_FixIn%Bool%TS_DO_SS_EXTERNAL

!  New 17 May 2012
      DO_SURFACE_LEAVING     = VLIDORT_FixIn%Bool%TS_DO_SURFACE_LEAVING
      DO_SL_ISOTROPIC        = VLIDORT_FixIn%Bool%TS_DO_SL_ISOTROPIC

!  Modified Boolean inputs

      DO_SSCORR_NADIR        = VLIDORT_ModIn%MBool%TS_DO_SSCORR_NADIR
      DO_SSCORR_OUTGOING     = VLIDORT_ModIn%MBool%TS_DO_SSCORR_OUTGOING
      DO_DOUBLE_CONVTEST     = VLIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST
      DO_SOLAR_SOURCES       = VLIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES
      DO_REFRACTIVE_GEOMETRY = VLIDORT_ModIn%MBool%TS_DO_REFRACTIVE_GEOMETRY
      DO_CHAPMAN_FUNCTION    = VLIDORT_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION
      DO_RAYLEIGH_ONLY       = VLIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY
      !DO_ISOTROPIC_ONLY      = VLIDORT_ModIn%MBool%TS_DO_ISOTROPIC_ONLY
      !DO_NO_AZIMUTH          = VLIDORT_ModIn%MBool%TS_DO_NO_AZIMUTH
      !DO_ALL_FOURIER         = VLIDORT_ModIn%MBool%TS_DO_ALL_FOURIER
      DO_DELTAM_SCALING      = VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING
      DO_SOLUTION_SAVING     = VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING
      DO_BVP_TELESCOPING     = VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING
      !DO_USER_STREAMS        = VLIDORT_ModIn%MBool%TS_DO_USER_STREAMS
      DO_USER_VZANGLES       = VLIDORT_ModIn%MBool%TS_DO_USER_VZANGLES
      DO_ADDITIONAL_MVOUT    = VLIDORT_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT
      DO_MVOUT_ONLY          = VLIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY
      DO_THERMAL_TRANSONLY   = VLIDORT_ModIn%MBool%TS_DO_THERMAL_TRANSONLY

!  Fixed Control inputs

      NSTOKES          = VLIDORT_FixIn%Cont%TS_NSTOKES
      NSTREAMS         = VLIDORT_FixIn%Cont%TS_NSTREAMS
      NLAYERS          = VLIDORT_FixIn%Cont%TS_NLAYERS
      NFINELAYERS      = VLIDORT_FixIn%Cont%TS_NFINELAYERS
      N_THERMAL_COEFFS = VLIDORT_FixIn%Cont%TS_N_THERMAL_COEFFS
      VLIDORT_ACCURACY = VLIDORT_FixIn%Cont%TS_VLIDORT_ACCURACY

      NLAYERS_NOMS     = VLIDORT_FixIn%Cont%TS_NLAYERS_NOMS
      NLAYERS_CUTOFF   = VLIDORT_FixIn%Cont%TS_NLAYERS_CUTOFF

!  Modified Control inputs

      NGREEK_MOMENTS_INPUT = VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT

!  Fixed Beam inputs

      !BEAM_SZAS   = VLIDORT_FixIn%Sunrays%TS_BEAM_SZAS
      SZANGLES    = VLIDORT_FixIn%Sunrays%TS_SZANGLES

!  Modified Beam inputs

      FLUX_FACTOR = VLIDORT_ModIn%MSunrays%TS_FLUX_FACTOR
      !NBEAMS      = VLIDORT_ModIn%MSunrays%TS_NBEAMS
      N_SZANGLES  = VLIDORT_ModIn%MSunrays%TS_N_SZANGLES

!  Fixed User Value inputs

      !N_USER_STREAMS      = VLIDORT_FixIn%UserVal%TS_N_USER_STREAMS
      N_USER_VZANGLES     = VLIDORT_FixIn%UserVal%TS_N_USER_VZANGLES

      N_USER_LEVELS       = VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS
      USER_LEVELS         = VLIDORT_FixIn%UserVal%TS_USER_LEVELS

!  Modified User Value inputs

      N_USER_RELAZMS      = VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS
      USER_RELAZMS        = VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS

      !USER_ANGLES         = VLIDORT_ModIn%MUserVal%TS_USER_ANGLES_INPUT
      USER_VZANGLES       = VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT

      GEOMETRY_SPECHEIGHT = VLIDORT_ModIn%MUserVal%TS_GEOMETRY_SPECHEIGHT

!  Fixed Chapman Function inputs

      HEIGHT_GRID         = VLIDORT_FixIn%Chapman%TS_HEIGHT_GRID
      PRESSURE_GRID       = VLIDORT_FixIn%Chapman%TS_PRESSURE_GRID
      TEMPERATURE_GRID    = VLIDORT_FixIn%Chapman%TS_TEMPERATURE_GRID
      FINEGRID            = VLIDORT_FixIn%Chapman%TS_FINEGRID
      RFINDEX_PARAMETER   = VLIDORT_FixIn%Chapman%TS_RFINDEX_PARAMETER

!  Modified Chapman Function inputs

      EARTH_RADIUS        = VLIDORT_ModIn%MChapman%TS_EARTH_RADIUS

!  Fixed Optical inputs

      DELTAU_VERT_INPUT     = VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT
      GREEKMAT_TOTAL_INPUT  = VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT
      THERMAL_BB_INPUT      = VLIDORT_FixIn%Optical%TS_THERMAL_BB_INPUT

!  @@@ Rob Fix 1/31/11, TS_EMISSIVITY, TS_USER_EMISSIVITY are defined in
!                       in Structure VLIDORT_SupIn%BRDF, so do not need
!                       to be copied here.
!      EMISSIVITY            = VLIDORT_FixIn%Optical%TS_EMISSIVITY
!      USER_EMISSIVITY       = VLIDORT_FixIn%Optical%TS_USER_EMISSIVITY

      LAMBERTIAN_ALBEDO     = VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO
      SURFBB                = VLIDORT_FixIn%Optical%TS_SURFACE_BB_INPUT

      LTE_DELTAU_VERT_INPUT = VLIDORT_FixIn%Optical%TS_LTE_DELTAU_VERT_INPUT
      LTE_THERMAL_BB_INPUT  = VLIDORT_FixIn%Optical%TS_LTE_THERMAL_BB_INPUT

!  Modified Optical inputs

      OMEGA_TOTAL_INPUT     = VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT

!  Fixed Write inputs

      DO_DEBUG_WRITE          = VLIDORT_FixIn%Write%TS_DO_DEBUG_WRITE

      DO_WRITE_INPUT          = VLIDORT_FixIn%Write%TS_DO_WRITE_INPUT
      INPUT_WRITE_FILENAME    = VLIDORT_FixIn%Write%TS_INPUT_WRITE_FILENAME

      DO_WRITE_SCENARIO       = VLIDORT_FixIn%Write%TS_DO_WRITE_SCENARIO
      SCENARIO_WRITE_FILENAME = VLIDORT_FixIn%Write%TS_SCENARIO_WRITE_FILENAME

      DO_WRITE_FOURIER        = VLIDORT_FixIn%Write%TS_DO_WRITE_FOURIER
      FOURIER_WRITE_FILENAME  = VLIDORT_FixIn%Write%TS_FOURIER_WRITE_FILENAME

      DO_WRITE_RESULTS        = VLIDORT_FixIn%Write%TS_DO_WRITE_RESULTS
      RESULTS_WRITE_FILENAME  = VLIDORT_FixIn%Write%TS_RESULTS_WRITE_FILENAME

!  BRDF Supplement Inputs

      EXACTDB_BRDFUNC = VLIDORT_Sup%BRDF%TS_EXACTDB_BRDFUNC
      BRDF_F_0        = VLIDORT_Sup%BRDF%TS_BRDF_F_0
      BRDF_F          = VLIDORT_Sup%BRDF%TS_BRDF_F
      USER_BRDF_F_0   = VLIDORT_Sup%BRDF%TS_USER_BRDF_F_0
      USER_BRDF_F     = VLIDORT_Sup%BRDF%TS_USER_BRDF_F

!  Surface-Leaving Supplement inputs
!    (This code introduced 17 May 2012)

      SLTERM_ISOTROPIC  = VLIDORT_Sup%SLEAVE%TS_SLTERM_ISOTROPIC
      SLTERM_USERANGLES = VLIDORT_Sup%SLEAVE%TS_SLTERM_USERANGLES
      SLTERM_F_0        = VLIDORT_Sup%SLEAVE%TS_SLTERM_F_0
      USER_SLTERM_F_0   = VLIDORT_Sup%SLEAVE%TS_USER_SLTERM_F_0

!  @@@ Rob fix 1/31/11, Emissivities from earlier structure
!                       were not copied --> wrong answers for BRDF cases
!        Lambertian case was OK, as used internal definitions

      EMISSIVITY      = VLIDORT_Sup%BRDF%TS_EMISSIVITY
      USER_EMISSIVITY = VLIDORT_Sup%BRDF%TS_USER_EMISSIVITY

!  new 12 March 2012
!   IF SS results already available copy them !
!  SS Inputs

      IF ( DO_SS_EXTERNAL ) THEN
         STOKES_SS = VLIDORT_Sup%SS%TS_STOKES_SS
         STOKES_DB = VLIDORT_Sup%SS%TS_STOKES_DB
      ENDIF

!  Fixed Linearized Control inputs

      DO_SIMULATION_ONLY       = &
        VLIDORT_LinFixIn%Cont%TS_DO_SIMULATION_ONLY
      DO_LTE_LINEARIZATION     = &
        VLIDORT_LinFixIn%Cont%TS_DO_LTE_LINEARIZATION
      DO_SURFBB_LINEARIZATION  = &
        VLIDORT_LinFixIn%Cont%TS_DO_SURFBB_LINEARIZATION

      LAYER_VARY_FLAG          = VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG
      LAYER_VARY_NUMBER        = VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER

      N_TOTALCOLUMN_WFS        = VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS
      N_TOTALPROFILE_WFS       = VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS
      N_SURFACE_WFS            = VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS

      COLUMNWF_NAMES           = VLIDORT_LinFixIn%Cont%TS_COLUMNWF_NAMES
      PROFILEWF_NAMES          = VLIDORT_LinFixIn%Cont%TS_PROFILEWF_NAMES

!  Modified Linearized Control inputs

      DO_COLUMN_LINEARIZATION  = &
        VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION
      DO_PROFILE_LINEARIZATION = &
        VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION
      DO_ATMOS_LINEARIZATION   = &
        VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION

      DO_SURFACE_LINEARIZATION = &
        VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION
      DO_LINEARIZATION         = &
        VLIDORT_LinModIn%MCont%TS_DO_LINEARIZATION

!  Linearized Optical inputs

      L_DELTAU_VERT_INPUT    = &
        VLIDORT_LinFixIn%Optical%TS_L_DELTAU_VERT_INPUT
      L_OMEGA_TOTAL_INPUT    = &
        VLIDORT_LinFixIn%Optical%TS_L_OMEGA_TOTAL_INPUT
      L_GREEKMAT_TOTAL_INPUT = &
        VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT

!  @@@ Rob Fix 1/31/11, TS_LS_EMISSIVITY, TS_LS_USER_EMISSIVITY are defined in
!                       in Structure VLIDORT_LinSup%BRDF, so do not need
!                       to be copied here.
!      LS_EMISSIVITY          = VLIDORT_LinFixIn%Optical%TS_LS_EMISSIVITY
!      LS_USER_EMISSIVITY     = VLIDORT_LinFixIn%Optical%TS_LS_USER_EMISSIVITY

!  BRDF Linearized Supplement inputs

      LS_EXACTDB_BRDFUNC = VLIDORT_LinSup%BRDF%TS_LS_EXACTDB_BRDFUNC
      LS_BRDF_F_0        = VLIDORT_LinSup%BRDF%TS_LS_BRDF_F_0
      LS_BRDF_F          = VLIDORT_LinSup%BRDF%TS_LS_BRDF_F
      LS_USER_BRDF_F_0   = VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F_0
      LS_USER_BRDF_F     = VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F

!  @@@ Rob fix 1/31/11, Linearized Emissivities from earlier structure
!                       were not copied --> wrong answers for BRDF cases
!        Lambertian case was OK, as used internal definitions

      LS_EMISSIVITY      = VLIDORT_LinSup%BRDF%TS_LS_EMISSIVITY
      LS_USER_EMISSIVITY = VLIDORT_LinSup%BRDF%TS_LS_USER_EMISSIVITY

!  new 12 March 2012
!   IF SS results already available copy them !
!  Linearized SS Inputs

      IF ( DO_SS_EXTERNAL ) THEN
         !SS atmosphere weighting functions

         COLUMNWF_SS   = VLIDORT_LinSup%SS%Atmos%TS_COLUMNWF_SS
         COLUMNWF_DB   = VLIDORT_LinSup%SS%Atmos%TS_COLUMNWF_DB

         PROFILEWF_SS  = VLIDORT_LinSup%SS%Atmos%TS_PROFILEWF_SS
         PROFILEWF_DB  = VLIDORT_LinSup%SS%Atmos%TS_PROFILEWF_DB

         !SS surface weighting functions

         !SURFACEWF_SS  = VLIDORT_LinSup%SS%Surf%TS_SURFACEWF_SS
         SURFACEWF_DB  = VLIDORT_LinSup%SS%Surf%TS_SURFACEWF_DB
      ENDIF

!  ==================================
!  END COPY INPUTS TO LOCAL VARIABLES
!  ==================================

!  VLIDORT input debug

      IF (DO_DEBUG_INPUT) THEN
        CALL VLIDORT_DEBUG_INPUT_MASTER()
        CALL VLIDORT_DEBUG_LIN_INPUT_MASTER()
      END IF

!  initialize outputs
!  ------------------

!  Status

      STATUS_CALCULATION = VLIDORT_SUCCESS
      STATUS_INPUTCHECK  = VLIDORT_SUCCESS

!  Model calculation

      MESSAGE = ' '
      TRACE_1 = ' '
      TRACE_2 = ' '
      TRACE_3 = ' '

!mick fix
!  Main outputs

      !Radiances and fluxes
      STOKES = ZERO
      MEAN_STOKES = ZERO
      FLUX_STOKES = ZERO
      MEAN_DIRECT = ZERO
      FLUX_DIRECT = ZERO

!  New 15 March 2012
      IF ( .NOT. DO_SS_EXTERNAL ) THEN
         STOKES_SS = ZERO
         STOKES_DB = ZERO
      ENDIF

      !Misc
      FOURIER_SAVED = 0
      N_GEOMETRIES  = 0
      SZA_OFFSETS   = 0
      VZA_OFFSETS   = 0

      !Column weighting functions
      COLUMNWF = ZERO

      !Profile weighting functions
      PROFILEWF = ZERO

      !Combined mean weighting function array
      MINT_ATMOSWF = ZERO
      FLUX_ATMOSWF = ZERO

      !LTE atmos weighting functions
      LTE_ATMOSWF = ZERO

      !Combined mean weighting function array (direct beam)
      MINT_ATMOSWF_DIRECT = ZERO
      FLUX_ATMOSWF_DIRECT = ZERO

      !Surface weighting functions
      SURFACEWF = ZERO
      MINT_SURFACEWF = ZERO
      FLUX_SURFACEWF = ZERO

!  New 15 March 2012
      IF ( .NOT. DO_SS_EXTERNAL ) THEN
         !SS atmosphere weighting functions
         COLUMNWF_SS = ZERO
         COLUMNWF_DB = ZERO

         PROFILEWF_SS = ZERO
         PROFILEWF_DB = ZERO

         !SS surface weighting functions
         !SURFACEWF_SS = ZERO
         SURFACEWF_DB = ZERO
      ENDIF

!  set multiplier

      SS_FLUX_MULTIPLIER = FLUX_FACTOR / PI4

!  Check input
!  -----------

!  Standard input variables check

      CALL VLIDORT_CHECK_INPUT ( &
        DO_SSFULL, DO_PLANE_PARALLEL, &
        DO_UPWELLING, DO_DNWELLING, &
        DO_QUAD_OUTPUT, DO_ADDITIONAL_MVOUT, &
        DO_MVOUT_ONLY, NSTOKES, &
        NSTREAMS, NLAYERS, &
        NFINELAYERS, SZANGLES, &
        USER_RELAZMS, N_USER_VZANGLES, &
        N_USER_LEVELS, USER_LEVELS, &
        HEIGHT_GRID, OMEGA_TOTAL_INPUT, &
        GREEKMAT_TOTAL_INPUT, DO_LAMBERTIAN_SURFACE, &
        DO_THERMAL_EMISSION, &
        DO_SPECIALIST_OPTION_2, DO_SPECIALIST_OPTION_3, &
        NLAYERS_NOMS, NLAYERS_CUTOFF, &
        DO_TOA_CONTRIBS, DO_NO_AZIMUTH, DO_SS_EXTERNAL, &
        DO_SSCORR_NADIR, DO_SSCORR_OUTGOING, &
        DO_SOLAR_SOURCES, DO_REFRACTIVE_GEOMETRY, &
        DO_CHAPMAN_FUNCTION, DO_RAYLEIGH_ONLY, &
        DO_DELTAM_SCALING, DO_SOLUTION_SAVING, &
        DO_BVP_TELESCOPING, DO_USER_VZANGLES, &
        NGREEK_MOMENTS_INPUT, N_SZANGLES, &
        EARTH_RADIUS, GEOMETRY_SPECHEIGHT, &
        N_USER_RELAZMS, USER_VZANGLES, &
        DO_THERMAL_TRANSONLY, DO_ALL_FOURIER, &
        DO_DIRECT_BEAM, DO_CLASSICAL_SOLUTION, &
        N_OUT_STREAMS, OUT_ANGLES, &
        STATUS_SUB, NCHECKMESSAGES, CHECKMESSAGES, ACTIONS )

!  Exception handling

      IF ( STATUS_SUB .EQ. VLIDORT_SERIOUS ) THEN
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
        VLIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
        VLIDORT_Out%Status%TS_NCHECKMESSAGES = NCHECKMESSAGES
        VLIDORT_Out%Status%TS_CHECKMESSAGES  = CHECKMESSAGES
        VLIDORT_Out%Status%TS_ACTIONS        = ACTIONS
        RETURN
      ELSE IF ( STATUS_SUB .EQ. VLIDORT_WARNING ) THEN
        STATUS_INPUTCHECK = VLIDORT_WARNING
        VLIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
!  ########## Rob Change, 3/28/2011 ################
!  Program will execute - these outputs are set at the end.
!          VLIDORT_Out%Status%TS_NCHECKMESSAGES = NCHECKMESSAGES
!          VLIDORT_Out%Status%TS_CHECKMESSAGES  = CHECKMESSAGES
!          VLIDORT_Out%Status%TS_ACTIONS        = ACTIONS
!  ########## Rob Change, 3/28/2011 ################
      ENDIF

!  Extended input variables check

      IF ( DO_LINEARIZATION ) THEN

        CALL VLIDORT_L_CHECK_INPUT ( &
          DO_SIMULATION_ONLY, DO_SURFBB_LINEARIZATION, &
          N_SURFACE_WFS, &
          DO_PROFILE_LINEARIZATION, DO_COLUMN_LINEARIZATION, &
          DO_ATMOS_LINEARIZATION, DO_SURFACE_LINEARIZATION, &
          STATUS_SUB, NCHECKMESSAGES, CHECKMESSAGES, ACTIONS )

        IF ( STATUS_SUB .EQ. VLIDORT_SERIOUS ) THEN
          STATUS_INPUTCHECK = VLIDORT_SERIOUS
          VLIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
          VLIDORT_Out%Status%TS_NCHECKMESSAGES = NCHECKMESSAGES
          VLIDORT_Out%Status%TS_CHECKMESSAGES  = CHECKMESSAGES
          VLIDORT_Out%Status%TS_ACTIONS        = ACTIONS
          RETURN
        ELSE IF ( STATUS_SUB .EQ. VLIDORT_WARNING ) THEN
          STATUS_INPUTCHECK = VLIDORT_WARNING
          VLIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
!  ########## Rob Change, 3/28/2011 ################
!  Program will execute - these outputs are set at the end.
!          VLIDORT_Out%Status%TS_NCHECKMESSAGES = NCHECKMESSAGES
!          VLIDORT_Out%Status%TS_CHECKMESSAGES  = CHECKMESSAGES
!          VLIDORT_Out%Status%TS_ACTIONS        = ACTIONS
!  ########## Rob Change, 3/28/2011 ################
        ENDIF

      ENDIF

!  Geometry adjustment
!  -------------------

!  Adjust surface condition

      ADJUST_SURFACE = .FALSE.
      IF ( DO_SSCORR_OUTGOING ) THEN
        IF (HEIGHT_GRID(NLAYERS).GT.GEOMETRY_SPECHEIGHT ) THEN
         ADJUST_SURFACE = .TRUE.
        ENDIF
      ENDIF

!  Perform adjustment

      MODIFIED_ERADIUS = EARTH_RADIUS + GEOMETRY_SPECHEIGHT
      CALL MULTI_OUTGOING_ADJUSTGEOM &
         ( MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS, &
           N_USER_VZANGLES,   N_SZANGLES,   N_USER_RELAZMS, &
           HEIGHT_GRID(NLAYERS), MODIFIED_ERADIUS, ADJUST_SURFACE, &
           USER_VZANGLES,  SZANGLES, USER_RELAZMS, &
           USER_VZANGLES_ADJUST, SZANGLES_ADJUST, USER_RELAZMS_ADJUST, &
           FAIL, MESSAGE, TRACE_1 )

!  Update exception handling. October 2010 2p4RTC

      if ( fail ) then
        TRACE_2 = ' Failure in multi_outgoing_adjustgeom'
        TRACE_3 = ' ** VLIDORT_L_MASTER '
        STATUS_CALCULATION = VLIDORT_SERIOUS
        VLIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION
        VLIDORT_Out%Status%TS_MESSAGE = MESSAGE
        VLIDORT_Out%Status%TS_TRACE_1 = TRACE_1
        VLIDORT_Out%Status%TS_TRACE_2 = TRACE_2
        VLIDORT_Out%Status%TS_TRACE_3 = TRACE_3
        RETURN
      ENDIF

!  Chapman function calculation
!  ----------------------------

!mick fix - comment out 1st IF
      !IF ( DO_SOLAR_SOURCES ) THEN
        IF ( DO_CHAPMAN_FUNCTION ) THEN

          CALL VLIDORT_CHAPMAN ( &
            DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, &
            NLAYERS, N_SZANGLES, &
            SZANGLES, &
            EARTH_RADIUS, RFINDEX_PARAMETER, &
            HEIGHT_GRID, PRESSURE_GRID, &
            TEMPERATURE_GRID, FINEGRID, &
            SZA_LOCAL_INPUT, CHAPMAN_FACTORS, &
            FAIL, MESSAGE, TRACE_1 )

          IF (FAIL) THEN
            TRACE_2 = 'Direct call in VLIDORT_L_MASTER'
            TRACE_3 = ' ** VLIDORT_L_MASTER '
            STATUS_CALCULATION = VLIDORT_SERIOUS
            VLIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION
            VLIDORT_Out%Status%TS_MESSAGE = MESSAGE
            VLIDORT_Out%Status%TS_TRACE_1 = TRACE_1
            VLIDORT_Out%Status%TS_TRACE_2 = TRACE_2
            VLIDORT_Out%Status%TS_TRACE_3 = TRACE_3
            RETURN
          ENDIF

        ENDIF
      !ENDIF

!  write input variables
!  ---------------------

!  open file, call standard/linearized input writes, close file

      IF ( DO_WRITE_INPUT ) THEN
        IUNIT = VLIDORT_INUNIT
        OPEN(IUNIT,FILE=INPUT_WRITE_FILENAME,STATUS='UNKNOWN')

        CALL VLIDORT_WRITEINPUT ( &
          IUNIT, &
          DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, &
          DO_RAYLEIGH_ONLY, DO_QUAD_OUTPUT, &
          DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY, &
          NSTREAMS, NLAYERS, &
          NGREEK_MOMENTS_INPUT, VLIDORT_ACCURACY, &
          FLUX_FACTOR, N_USER_RELAZMS, &
          N_USER_LEVELS, &
          DO_LAMBERTIAN_SURFACE, &
          DO_THERMAL_EMISSION, &
          DO_SURFACE_EMISSION, &
          DO_DIRECT_BEAM, DO_CLASSICAL_SOLUTION, &
          DO_NO_AZIMUTH, N_SZANGLES, &
          N_USER_VZANGLES, DO_USER_VZANGLES )

        IF ( DO_LINEARIZATION ) THEN
          CALL VLIDORT_L_WRITEINPUT ( &
            IUNIT, &
            DO_SIMULATION_ONLY, DO_PROFILE_LINEARIZATION, &
            DO_COLUMN_LINEARIZATION, N_TOTALCOLUMN_WFS, &
            N_TOTALPROFILE_WFS, DO_SURFACE_LINEARIZATION, &
            DO_SURFBB_LINEARIZATION, PROFILEWF_NAMES, &
            COLUMNWF_NAMES )
        ENDIF
        CLOSE(IUNIT)
      ENDIF

!  Get derived inputs
!  ==================

!  Miscellaneous and layer input.

      CALL VLIDORT_DERIVE_INPUT ( &
        DO_FULLRAD_MODE, DO_SSCORR_NADIR, DO_SSCORR_OUTGOING, DO_SSFULL, &
        DO_SOLAR_SOURCES, DO_REFRACTIVE_GEOMETRY, DO_RAYLEIGH_ONLY, &
        DO_UPWELLING, DO_DNWELLING, DO_USER_VZANGLES, NSTOKES, &
        NSTREAMS, NLAYERS, NGREEK_MOMENTS_INPUT, &
        N_SZANGLES, SZANGLES, SZA_LOCAL_INPUT, N_USER_RELAZMS, &
        N_USER_VZANGLES, USER_VZANGLES, N_USER_LEVELS, USER_LEVELS, &
        DELTAU_VERT_INPUT, &
        GREEKMAT_TOTAL_INPUT, DO_SURFACE_EMISSION, SURFBB, & 
        DO_THERMAL_TRANSONLY, DO_SPECIALIST_OPTION_2, NLAYERS_NOMS, &
        USER_VZANGLES_ADJUST, N_OUT_STREAMS, DO_SOLUTION_SAVING, &
        DO_BVP_TELESCOPING, DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY, &
        FLUX_FACTOR, OMEGA_TOTAL_INPUT, DO_DOUBLE_CONVTEST, &
        DO_ALL_FOURIER, DO_DIRECT_BEAM, DO_CLASSICAL_SOLUTION, & 
        DO_DBCORRECTION, DO_NO_AZIMUTH, FLUXVEC, COS_SZANGLES, &
        SIN_SZANGLES, SUN_SZA_COSINES, QUAD_STREAMS, QUAD_WEIGHTS, &
        QUAD_STRMWTS, QUAD_HALFWTS, QUAD_SINES, QUAD_ANGLES, &
        DO_MSMODE_VLIDORT, NMOMENTS, NSTREAMS_2, &
        NTOTAL, N_SUBDIAG, N_SUPDIAG, NSTOKES_SQ, &
        NSTKS_NSTRMS, NSTKS_NSTRMS_2, NBEAMS, &
        NPARTICSOLS, MUELLER_INDEX, DMAT, GREEKMAT_INDEX, DO_REAL_EIGENSOLVER, &
        BVP_REGULAR_FLAG, LAYER_MAXMOMENTS, DO_LAYER_SCATTERING, N_CONVTESTS, &
        N_CONV_STREAMS, N_DIRECTIONS, WHICH_DIRECTIONS, &
        USER_STREAMS, USER_SINES, USER_SECANTS, PARTLAYERS_OUTFLAG, &
        PARTLAYERS_OUTINDEX, UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, & 
        DO_PARTLAYERS, N_PARTLAYERS, PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES, &
        N_LAYERSOURCE_UP, N_LAYERSOURCE_DN, N_ALLLAYERS_UP, &
        N_ALLLAYERS_DN, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, &
        DFLUX, TAUGRID_INPUT, &
        STATUS_SUB, MESSAGE )

      DO_MSMODE_THERMAL = (.not.DO_FULLRAD_MODE) .and. &
             ( DO_SURFACE_EMISSION .and.DO_THERMAL_EMISSION )

!  Exception handling

!  ########## Rob Change, 3/28/2011 ################
!    Output from DERIVE_INPUTS are checks, not execution failures
!  Old code---------
!      IF ( STATUS_SUB .EQ. VLIDORT_WARNING ) THEN
!        TRACE_1 = ''
!        TRACE_2 = 'Derive_Input Call in VLIDORT_L_MASTER'
!        TRACE_3 = ' ** VLIDORT_L_MASTER'
!        STATUS_INPUTCHECK = VLIDORT_WARNING
!        VLIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION
!        VLIDORT_Out%Status%TS_MESSAGE = MESSAGE
!        VLIDORT_Out%Status%TS_TRACE_1 = TRACE_1
!        VLIDORT_Out%Status%TS_TRACE_2 = TRACE_2
!        VLIDORT_Out%Status%TS_TRACE_3 = TRACE_3
!      ENDIF
!  New code-----------
      IF ( STATUS_SUB .EQ. VLIDORT_WARNING ) THEN
        STATUS_INPUTCHECK = VLIDORT_WARNING
        VLIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
        NCHECKMESSAGES = NCHECKMESSAGES + 1
        CHECKMESSAGES(NCHECKMESSAGES)  = TRIM(ADJUSTL(MESSAGE))
        ACTIONS(NCHECKMESSAGES) &
          = ' Action taken in VLIDORT_DERIVE_INPUT to set internal default'
      ENDIF
!  ########## Rob Change, 3/28/2011 ################

!  Initialise Fourier loop
!  =======================

!  Single scatter correction: flux multiplier

      IF ( DO_SSCORR_NADIR .OR. DO_SSCORR_OUTGOING ) THEN
        SS_FLUX_MULTIPLIER = FLUX_FACTOR / PI4
      ELSE
        SS_FLUX_MULTIPLIER = FLUX_FACTOR / PI4
      ENDIF

!  Set Number of Fourier terms (NMOMENTS = Maximum).
!    ( Starting from 0 = Fundamental )

      SAVE_DO_NO_AZIMUTH  = DO_NO_AZIMUTH
      LOCAL_DO_NO_AZIMUTH = DO_NO_AZIMUTH

!  No azimuth dependency for following cases
!   Other cases now disabled, 17 Janaury 2006

!      IF ( DO_TRANSMITTANCE_ONLY ) THEN
!        LOCAL_DO_NO_AZIMUTH = .TRUE.
!      ENDIF

!      IF ( DO_ISOTROPIC_ONLY  ) THEN
!        LOCAL_DO_NO_AZIMUTH = .TRUE.
!      ENDIF

      IF ( .NOT. DO_SOLAR_SOURCES  ) THEN
        LOCAL_DO_NO_AZIMUTH = .TRUE.
      ENDIF

      IF ( DO_MVOUT_ONLY  ) THEN
        LOCAL_DO_NO_AZIMUTH = .TRUE.
      ENDIF

!  set Fourier number (2 for Rayleigh only)

      IF ( LOCAL_DO_NO_AZIMUTH .OR. DO_SSFULL ) THEN
        N_FOURIER_COMPONENTS = 0
      ELSE
        IF ( DO_RAYLEIGH_ONLY  ) THEN
          N_FOURIER_COMPONENTS = 2
        ELSE
          N_FOURIER_COMPONENTS = NMOMENTS
        ENDIF
      ENDIF

!  re-set no-azimuth flag

      DO_NO_AZIMUTH = LOCAL_DO_NO_AZIMUTH

!  if there's no azimuth dependence, just do one value in azimuth loop

      IF ( DO_NO_AZIMUTH ) THEN
        LOCAL_N_USERAZM = 1
       ELSE
        LOCAL_N_USERAZM = N_USER_RELAZMS
      ENDIF

!  Number of sources

      IF ( DO_SOLAR_SOURCES ) THEN
        NSOURCES = NBEAMS
      ELSE
        NSOURCES = 1
        NBEAMS   = 1
      ENDIF

!  save some offsets for indexing geometries

      N_VIEWING    = N_USER_VZANGLES * LOCAL_N_USERAZM
      N_GEOMETRIES = NSOURCES * N_VIEWING
      DO IBEAM = 1, NBEAMS
        SZA_OFFSETS(IBEAM) = N_VIEWING * ( IBEAM - 1 )
        DO UM = 1, N_USER_VZANGLES
          VZA_OFFSETS(IBEAM,UM) = &
                 SZA_OFFSETS(IBEAM) + LOCAL_N_USERAZM * (UM - 1)
        END DO
      END DO

!  Fourier loop
!  ============

!  Initialize

      LOCAL_ITERATION   = .TRUE.
      FOURIER_COMPONENT = -1
      TESTCONV          = 0

!  set up solar beam flags. Required in all cases.
!   ---Even required for the thermal.....

      DO IBEAM = 1, NBEAMS
        BEAM_TESTCONV  ( IBEAM )  = 0
        BEAM_ITERATION ( IBEAM ) = .TRUE.
        DO L = 0, MAXFOURIER
          DO_MULTIBEAM   ( IBEAM, L ) = .TRUE.
        ENDDO
      ENDDO

!  start loop

      DO WHILE ( LOCAL_ITERATION .AND. &
                   FOURIER_COMPONENT.LT.N_FOURIER_COMPONENTS )

!  Fourier counter

        FOURIER_COMPONENT = FOURIER_COMPONENT + 1

!  Local start of user-defined streams. Should always be 1.
!    No zenith tolerance now.

        LOCAL_UM_START = 1

!  azimuth cosine/sine factors, using adjust geometries.

        IF ( FOURIER_COMPONENT .GT. 0 ) THEN
          DFC = DBLE(FOURIER_COMPONENT)
          DO UA = 1, LOCAL_N_USERAZM
            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_VZANGLES
                AZM_ARGUMENT = USER_RELAZMS_ADJUST(UM,IB,UA) * DFC
                AZMFAC(UM,IB,UA,1)   = DCOS(DEG_TO_RAD*AZM_ARGUMENT)
                AZMFAC(UM,IB,UA,2)   = AZMFAC(UM,IB,UA,1)
                AZMFAC(UM,IB,UA,3)   = DSIN(DEG_TO_RAD*AZM_ARGUMENT)
                AZMFAC(UM,IB,UA,4)   = AZMFAC(UM,IB,UA,3)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  Main call to VLidort Fourier module

!        write(*,*)' ..calculating fourier component',FOURIER_COMPONENT

        CALL VLIDORT_L_FOURIER ( &
          FOURIER_COMPONENT, SS_FLUX_MULTIPLIER, DO_SSCORR_NADIR, DO_SSCORR_OUTGOING, &
          DO_SSCORR_TRUNCATION, DO_SSFULL, DO_SOLAR_SOURCES, DO_PLANE_PARALLEL, &
          DO_REFRACTIVE_GEOMETRY, DO_DELTAM_SCALING, DO_SOLUTION_SAVING, DO_BVP_TELESCOPING, &
          DO_UPWELLING, DO_DNWELLING, DO_QUAD_OUTPUT, DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY, &
          DO_DEBUG_WRITE, DO_FDTEST, NSTOKES, NSTREAMS, NLAYERS, NGREEK_MOMENTS_INPUT, &
          FLUX_FACTOR, SZA_LOCAL_INPUT, N_USER_RELAZMS, USER_RELAZMS, N_USER_LEVELS, &
          NFINELAYERS, EARTH_RADIUS, HEIGHT_GRID, &
          OMEGA_TOTAL_INPUT, DELTAU_VERT_INPUT, GREEKMAT_TOTAL_INPUT, &
          DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
          EXACTDB_BRDFUNC, BRDF_F, BRDF_F_0, USER_BRDF_F, USER_BRDF_F_0, &
          EMISSIVITY, USER_EMISSIVITY, DO_THERMAL_EMISSION, N_THERMAL_COEFFS, THERMAL_BB_INPUT, &
          DO_SURFACE_EMISSION, SURFBB, DO_THERMAL_TRANSONLY, DO_SPECIALIST_OPTION_2, DO_SPECIALIST_OPTION_3, &
          DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, SLTERM_ISOTROPIC, SLTERM_USERANGLES, SLTERM_F_0, USER_SLTERM_F_0, &
          NLAYERS_CUTOFF, DO_TOA_CONTRIBS, DO_DIRECT_BEAM, &
          DO_CLASSICAL_SOLUTION, DO_DBCORRECTION, FLUXVEC, &
          COS_SZANGLES, SIN_SZANGLES, SUN_SZA_COSINES, DO_MULTIBEAM, &
          QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS, QUAD_HALFWTS, &
          DO_MSMODE_VLIDORT, DO_MSMODE_THERMAL, NMOMENTS, NSTREAMS_2, NTOTAL, &
          N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2, NBEAMS, N_USER_VZANGLES, &
          MUELLER_INDEX, DMAT, DO_REAL_EIGENSOLVER, BVP_REGULAR_FLAG, LAYER_MAXMOMENTS, &
          DO_LAYER_SCATTERING, DO_USER_VZANGLES, LOCAL_UM_START, N_DIRECTIONS, WHICH_DIRECTIONS, &
          USER_VZANGLES_ADJUST, SZANGLES_ADJUST, USER_RELAZMS_ADJUST, USER_STREAMS, USER_SECANTS, &
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, &
          DO_PARTLAYERS, N_PARTLAYERS, PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES, &
          N_ALLLAYERS_UP, N_ALLLAYERS_DN, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, &
          DFLUX, N_GEOMETRIES, VZA_OFFSETS, TAUGRID_INPUT, CHAPMAN_FACTORS, &
          NSTOKES_SQ, DO_SIMULATION_ONLY, DO_PROFILE_LINEARIZATION, DO_COLUMN_LINEARIZATION, &
          DO_ATMOS_LINEARIZATION, DO_LTE_LINEARIZATION, N_TOTALCOLUMN_WFS, &
          LAYER_VARY_FLAG, LAYER_VARY_NUMBER, L_OMEGA_TOTAL_INPUT, L_DELTAU_VERT_INPUT, &
          L_GREEKMAT_TOTAL_INPUT, LTE_DELTAU_VERT_INPUT, LTE_THERMAL_BB_INPUT, &
          DO_SURFACE_LINEARIZATION, N_SURFACE_WFS, LS_EXACTDB_BRDFUNC, LS_BRDF_F, LS_BRDF_F_0, &
          LS_USER_BRDF_F, LS_USER_BRDF_F_0, LS_USER_EMISSIVITY, LS_EMISSIVITY, &
          TAUGRID, OMEGA_TOTAL, GREEKMAT_TOTAL, &
          STOKES_SS, STOKES_DB, SS_CONTRIBS, MEAN_STOKES, FLUX_STOKES, MEAN_DIRECT, FLUX_DIRECT, &
          MS_CONTRIBS_F, STOKES_F, DO_INCLUDE_SURFACE, DO_INCLUDE_SURFEMISS, DO_INCLUDE_THERMEMISS, &
          PROFILEWF_SS, COLUMNWF_SS, PROFILEWF_DB, COLUMNWF_DB, SURFACEWF_DB, &
          ATMOSWF_F, MINT_ATMOSWF, MINT_ATMOSWF_DIRECT, FLUX_ATMOSWF, FLUX_ATMOSWF_DIRECT, &
          SURFACEWF_F, MINT_SURFACEWF, FLUX_SURFACEWF, LTE_ATMOSWF, &
          STATUS_SUB, MESSAGE, TRACE_1, TRACE_2 )

!  error handling

        IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
          TRACE_3 = ' Called by VLIDORT_L_MASTER '
          STATUS_CALCULATION = VLIDORT_SERIOUS
          VLIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION
          VLIDORT_Out%Status%TS_MESSAGE = MESSAGE
          VLIDORT_Out%Status%TS_TRACE_1 = TRACE_1
          VLIDORT_Out%Status%TS_TRACE_2 = TRACE_2
          VLIDORT_Out%Status%TS_TRACE_3 = TRACE_3
          RETURN
        ENDIF

!  Fourier summation and Convergence examination
!  ---------------------------------------------

!   -- only done for beams which are still not converged.
!      This is controlled by flag DO_MULTIBEAM

!   -- new criterion, SS is added for Fourier = 0, as this means that
!      higher-order terms will be relatively smaller, which implies
!      faster convergence in some circumstances (generally not often).

        IBEAM_COUNT = 0
        DO IBEAM = 1, NBEAMS
          IF ( DO_MULTIBEAM ( IBEAM, FOURIER_COMPONENT ) ) THEN

!  Convergence and radiance summation

            CALL VLIDORT_CONVERGE ( &
              AZMFAC, FOURIER_COMPONENT, LOCAL_N_USERAZM, &
              DO_SSCORR_NADIR, DO_SSCORR_OUTGOING, &
              DO_SS_EXTERNAL, &
              DO_SSFULL, DO_DOUBLE_CONVTEST, &
              DO_RAYLEIGH_ONLY, DO_UPWELLING, &
              NSTOKES, NSTREAMS, &
              NLAYERS, VLIDORT_ACCURACY, &
              N_USER_RELAZMS, N_USER_LEVELS, &
              DO_TOA_CONTRIBS, &
              DO_ALL_FOURIER, DO_DBCORRECTION, &
              DO_NO_AZIMUTH, N_CONVTESTS, &
              LOCAL_UM_START, N_DIRECTIONS, &
              WHICH_DIRECTIONS, &
              N_OUT_STREAMS, VZA_OFFSETS, &
              STOKES_SS, STOKES_DB, SS_CONTRIBS, &
              STOKES_F, MS_CONTRIBS_F, IBEAM, &
              STOKES, FOURIER_SAVED, CONTRIBS, &
              BEAM_TESTCONV(IBEAM), BEAM_ITERATION(IBEAM) )

!  Fourier summation of linearization quantities

            IF ( DO_LINEARIZATION ) THEN
              CALL VLIDORT_L_CONVERGE ( &
                FOURIER_COMPONENT, LOCAL_N_USERAZM, &
                IBEAM, AZMFAC, &
                DO_SSCORR_NADIR, DO_SSCORR_OUTGOING, &
                DO_SS_EXTERNAL, & ! New 15 March 2012
                DO_SSFULL, DO_UPWELLING, &
                NSTOKES, NLAYERS, &
                N_USER_LEVELS, &
                DO_LAMBERTIAN_SURFACE, &
                DO_DBCORRECTION, DO_NO_AZIMUTH, &
                N_USER_VZANGLES, LOCAL_UM_START, &
                N_DIRECTIONS, WHICH_DIRECTIONS, &
                N_OUT_STREAMS, &
                VZA_OFFSETS, &
                DO_PROFILE_LINEARIZATION, DO_COLUMN_LINEARIZATION, &
                N_TOTALCOLUMN_WFS, LAYER_VARY_FLAG, &
                LAYER_VARY_NUMBER, DO_SURFACE_LINEARIZATION, &
                N_SURFACE_WFS, &
                PROFILEWF_SS, COLUMNWF_SS, &
                PROFILEWF_DB, COLUMNWF_DB, SURFACEWF_DB, &
                SURFACEWF_F, ATMOSWF_F, &
                SURFACEWF, PROFILEWF, COLUMNWF )
            ENDIF

!  Check number of beams already converged

            IF ( BEAM_ITERATION(IBEAM) ) THEN
              IBEAM_COUNT = IBEAM_COUNT + 1
            ELSE
              DO L = FOURIER_COMPONENT+1,MAXFOURIER
                DO_MULTIBEAM (IBEAM,L) = .FALSE.
              ENDDO
            ENDIF

!  end beam count loop

          ENDIF
        END DO

!  If all beams have converged, stop iteration

        IF ( IBEAM_COUNT .EQ. 0 ) LOCAL_ITERATION = .FALSE.

!  Fourier output
!  --------------

!  Open file if Fourier = 0
!  Write Standard   Fourier output (radiances)
!  Write additional Fourier output (linearizations)
!  Close file if iteration has finished
!  New comment:
!    If the SS correction is set, Fourier=0 will include SS field

!        IF ( DO_WRITE_FOURIER ) THEN
!          FUNIT = VLIDORT_FUNIT
!          IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
!            OPEN(FUNIT,FILE=FOURIER_WRITE_FILENAME,STATUS='UNKNOWN')
!          ENDIF
!          CALL VLIDORT_WRITEFOURIER ( FUNIT, FOURIER_COMPONENT )
!          IF ( DO_LINEARIZATION ) THEN
!            CALL VLIDORT_L_WRITEFOURIER
!     I         ( DO_INCLUDE_SURFACE, FUNIT, FOURIER_COMPONENT )
!          ENDIF
!          IF ( .NOT.LOCAL_ITERATION ) CLOSE ( FUNIT )
!        ENDIF

!  end iteration loop

      ENDDO

!  restore no azimuth flag

      DO_NO_AZIMUTH = SAVE_DO_NO_AZIMUTH

!  Major result output
!  ===================

!  Standard + linearized output
!  ----------------------------

      IF ( DO_WRITE_RESULTS ) THEN
        RUNIT = VLIDORT_RESUNIT
        OPEN(RUNIT,FILE=RESULTS_WRITE_FILENAME,STATUS='UNKNOWN')

        CALL VLIDORT_WRITERESULTS ( &
          RUNIT, &
          DO_FULLRAD_MODE, DO_SSCORR_NADIR, &
          DO_SSCORR_OUTGOING, DO_DOUBLE_CONVTEST, &
          DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, &
          DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY, &
          NSTOKES, VLIDORT_ACCURACY, &
          SZANGLES, N_USER_RELAZMS, &
          USER_RELAZMS, N_USER_LEVELS, &
          USER_LEVELS, HEIGHT_GRID, &
          DELTAU_VERT_INPUT, &
          DO_CLASSICAL_SOLUTION, DO_NO_AZIMUTH, &
          NBEAMS, N_DIRECTIONS, &
          WHICH_DIRECTIONS, N_OUT_STREAMS, &
          OUT_ANGLES, PARTLAYERS_OUTFLAG, &
          VZA_OFFSETS, &
          TAUGRID_INPUT, DO_MULTIBEAM, &
          STOKES, MEAN_STOKES, &
          FLUX_STOKES, FOURIER_SAVED )

        IF ( DO_LINEARIZATION ) THEN
          CALL VLIDORT_L_WRITERESULTS ( &
            RUNIT, &
            DO_FULLRAD_MODE, DO_SSCORR_NADIR, &
            DO_SSCORR_OUTGOING, DO_DOUBLE_CONVTEST, &
            DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, &
            DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY, &
            NSTOKES, NLAYERS, &
            VLIDORT_ACCURACY, SZANGLES, &
            N_USER_RELAZMS, USER_RELAZMS, &
            N_USER_LEVELS, USER_LEVELS, &
            DO_LAMBERTIAN_SURFACE, &
            DO_NO_AZIMUTH, NBEAMS, &
            N_DIRECTIONS, WHICH_DIRECTIONS, &
            N_OUT_STREAMS, &
            OUT_ANGLES, VZA_OFFSETS, &
            DO_MULTIBEAM, FOURIER_SAVED, &
            DO_PROFILE_LINEARIZATION, DO_COLUMN_LINEARIZATION, &
            N_TOTALCOLUMN_WFS, LAYER_VARY_FLAG, &
            LAYER_VARY_NUMBER, DO_SURFACE_LINEARIZATION, &
            DO_SURFBB_LINEARIZATION, N_SURFACE_WFS, &
            PROFILEWF_NAMES, COLUMNWF_NAMES, &
            SURFACEWF, PROFILEWF, &
            COLUMNWF, MINT_SURFACEWF, &
            MINT_ATMOSWF, FLUX_SURFACEWF, &
            FLUX_ATMOSWF )
        ENDIF

        CLOSE(RUNIT)
      ENDIF

!  Geophysical input (scenario) write
!  ----------------------------------

      IF ( DO_WRITE_SCENARIO ) THEN
        SUNIT = VLIDORT_SCENUNIT
        OPEN(SUNIT,FILE=SCENARIO_WRITE_FILENAME,STATUS='UNKNOWN')

        CALL VLIDORT_WRITESCEN ( &
          SUNIT, &
          DO_DELTAM_SCALING, NSTREAMS, &
          NLAYERS, NGREEK_MOMENTS_INPUT, &
          N_SZANGLES, SZANGLES, &
          N_USER_RELAZMS, USER_RELAZMS, &
          N_USER_VZANGLES, USER_VZANGLES, &
          N_USER_LEVELS, USER_LEVELS, &
          OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT, &
          DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
          DO_THERMAL_EMISSION, N_THERMAL_COEFFS, &
          DO_SURFACE_EMISSION, SURFBB, &
          QUAD_STREAMS, QUAD_WEIGHTS, &
          QUAD_ANGLES, DO_NO_AZIMUTH, &
          NMOMENTS, GREEKMAT_INDEX, &
          TAUGRID_INPUT, OMEGA_TOTAL, &
          GREEKMAT_TOTAL, TAUGRID )

        IF ( DO_LINEARIZATION ) THEN
          CALL VLIDORT_L_WRITESCEN ( &
            SUNIT, NLAYERS, DO_ATMOS_LINEARIZATION, &
            LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
            L_OMEGA_TOTAL_INPUT, L_DELTAU_VERT_INPUT )
        ENDIF

        CLOSE(SUNIT)
      ENDIF

!  ========================================================
!  BEGIN COPY LOCAL VARIABLES TO OUTPUTS (IN/OUT variables)
!  ========================================================

!  Standard variables

      VLIDORT_ModIn%MBool%TS_DO_SSCORR_NADIR        = DO_SSCORR_NADIR
      VLIDORT_ModIn%MBool%TS_DO_SSCORR_OUTGOING     = DO_SSCORR_OUTGOING

      VLIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST     = DO_DOUBLE_CONVTEST
      VLIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES       = DO_SOLAR_SOURCES

      VLIDORT_ModIn%MBool%TS_DO_REFRACTIVE_GEOMETRY = DO_REFRACTIVE_GEOMETRY
      VLIDORT_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION    = DO_CHAPMAN_FUNCTION

      VLIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY       = DO_RAYLEIGH_ONLY
      !VLIDORT_ModIn%MBool%TS_DO_ISOTROPIC_ONLY      = DO_ISOTROPIC_ONLY
      !VLIDORT_ModIn%MBool%TS_DO_NO_AZIMUTH          = DO_NO_AZIMUTH
      !VLIDORT_ModIn%MBool%TS_DO_ALL_FOURIER         = DO_ALL_FOURIER

      VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING      = DO_DELTAM_SCALING

      VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING     = DO_SOLUTION_SAVING
      VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING     = DO_BVP_TELESCOPING

      VLIDORT_ModIn%MBool%TS_DO_USER_VZANGLES       = DO_USER_VZANGLES

      VLIDORT_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT    = DO_ADDITIONAL_MVOUT
      VLIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY          = DO_MVOUT_ONLY

      VLIDORT_ModIn%MBool%TS_DO_THERMAL_TRANSONLY   = DO_THERMAL_TRANSONLY

      VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT   = NGREEK_MOMENTS_INPUT

      VLIDORT_ModIn%MSunrays%TS_FLUX_FACTOR         = FLUX_FACTOR
      VLIDORT_ModIn%MSunrays%TS_N_SZANGLES          = N_SZANGLES

      VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS      = N_USER_RELAZMS
      VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS        = USER_RELAZMS
      VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT = USER_VZANGLES
      VLIDORT_ModIn%MUserVal%TS_GEOMETRY_SPECHEIGHT = GEOMETRY_SPECHEIGHT

      VLIDORT_ModIn%MChapman%TS_EARTH_RADIUS        = EARTH_RADIUS

      VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT   = OMEGA_TOTAL_INPUT

!  Linearized variables

      VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION  = &
        DO_COLUMN_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION = &
        DO_PROFILE_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION   = &
        DO_ATMOS_LINEARIZATION

      VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION = &
        DO_SURFACE_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_LINEARIZATION = &
        DO_LINEARIZATION

!  ==========================================================
!  BEGIN COPY LOCAL VARIABLES TO OUTPUTS (pure OUT variables)
!  ==========================================================

!  Radiances and fluxes
!    Output direct Mean/Flux added 17 May 2012

      VLIDORT_Out%Main%TS_STOKES      = STOKES
      VLIDORT_Out%Main%TS_MEAN_STOKES = MEAN_STOKES
      VLIDORT_Out%Main%TS_FLUX_STOKES = FLUX_STOKES

      VLIDORT_Out%Main%TS_MEAN_DIRECT = MEAN_DIRECT
      VLIDORT_Out%Main%TS_FLUX_DIRECT = FLUX_DIRECT

!  new 12 March 2012
!   IF SS results already available, no need to copy them !

      IF ( .NOT. DO_SS_EXTERNAL ) THEN
         VLIDORT_Sup%SS%TS_STOKES_SS = STOKES_SS
         VLIDORT_Sup%SS%TS_STOKES_DB = STOKES_DB
      ENDIF

!  Bookkeeping

      VLIDORT_Out%Main%TS_FOURIER_SAVED = FOURIER_SAVED
      VLIDORT_Out%Main%TS_N_GEOMETRIES  = N_GEOMETRIES
      VLIDORT_Out%Main%TS_SZA_OFFSETS   = SZA_OFFSETS

!  Atmosphere weighting functions
!    Output direct Mean/Flux added 17 May 2012

      VLIDORT_LinOut%Atmos%TS_COLUMNWF      = COLUMNWF
      VLIDORT_LinOut%Atmos%TS_PROFILEWF     = PROFILEWF
      VLIDORT_LinOut%Atmos%TS_MINT_ATMOSWF  = MINT_ATMOSWF
      VLIDORT_LinOut%Atmos%TS_FLUX_ATMOSWF  = FLUX_ATMOSWF
      VLIDORT_LinOut%Atmos%TS_LTE_ATMOSWF   = LTE_ATMOSWF

      VLIDORT_LinOut%Atmos%TS_MINT_ATMOSWF_DIRECT  = MINT_ATMOSWF_DIRECT
      VLIDORT_LinOut%Atmos%TS_FLUX_ATMOSWF_DIRECT  = FLUX_ATMOSWF_DIRECT

!  Surface weighting functions

      VLIDORT_LinOut%Surf%TS_SURFACEWF      = SURFACEWF
      VLIDORT_LinOut%Surf%TS_MINT_SURFACEWF = MINT_SURFACEWF
      VLIDORT_LinOut%Surf%TS_FLUX_SURFACEWF = FLUX_SURFACEWF

      IF ( .NOT. DO_SS_EXTERNAL ) THEN
         !SS atmosphere weighting functions

         VLIDORT_LinSup%SS%Atmos%TS_COLUMNWF_SS   = COLUMNWF_SS
         VLIDORT_LinSup%SS%Atmos%TS_COLUMNWF_DB   = COLUMNWF_DB

         VLIDORT_LinSup%SS%Atmos%TS_PROFILEWF_SS  = PROFILEWF_SS
         VLIDORT_LinSup%SS%Atmos%TS_PROFILEWF_DB  = PROFILEWF_DB

         !SS surface weighting functions

         !VLIDORT_LinSup%SS%Surf%TS_SURFACEWF_SS   = SURFACEWF_SS
         VLIDORT_LinSup%SS%Surf%TS_SURFACEWF_DB   = SURFACEWF_DB
      ENDIF

!  Exception handling

      VLIDORT_Out%Status%TS_STATUS_INPUTCHECK  = STATUS_INPUTCHECK
      VLIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION

      VLIDORT_Out%Status%TS_NCHECKMESSAGES = NCHECKMESSAGES
      VLIDORT_Out%Status%TS_CHECKMESSAGES  = CHECKMESSAGES
      VLIDORT_Out%Status%TS_ACTIONS        = ACTIONS

      VLIDORT_Out%Status%TS_MESSAGE = MESSAGE
      VLIDORT_Out%Status%TS_TRACE_1 = TRACE_1
      VLIDORT_Out%Status%TS_TRACE_2 = TRACE_2
      VLIDORT_Out%Status%TS_TRACE_3 = TRACE_3

!  ===================================
!  END COPY LOCAL VARIABLES TO OUTPUTS
!  ===================================

!  Finish

      RETURN

      CONTAINS

      SUBROUTINE VLIDORT_DEBUG_INPUT_MASTER()

      CALL VLIDORT_WRITE_STD_INPUT ( &
        DO_FULLRAD_MODE,DO_SSCORR_TRUNCATION,DO_SS_EXTERNAL,DO_SSFULL,&
        DO_THERMAL_EMISSION,DO_SURFACE_EMISSION,DO_PLANE_PARALLEL,&
        DO_UPWELLING,DO_DNWELLING,DO_QUAD_OUTPUT,&
        DO_TOA_CONTRIBS,DO_LAMBERTIAN_SURFACE,&
        DO_SPECIALIST_OPTION_1,DO_SPECIALIST_OPTION_2,DO_SPECIALIST_OPTION_3,&
        DO_SURFACE_LEAVING,DO_SL_ISOTROPIC,&
        NSTOKES,NSTREAMS,NLAYERS,&
        NFINELAYERS,N_THERMAL_COEFFS,VLIDORT_ACCURACY,&
        NLAYERS_NOMS,NLAYERS_CUTOFF,SZANGLES,&
        N_USER_VZANGLES,N_USER_LEVELS,USER_LEVELS,&
        HEIGHT_GRID,PRESSURE_GRID,TEMPERATURE_GRID,&
        FINEGRID,RFINDEX_PARAMETER,&
        DELTAU_VERT_INPUT,GREEKMAT_TOTAL_INPUT,THERMAL_BB_INPUT,&
        LAMBERTIAN_ALBEDO,SURFBB,LTE_DELTAU_VERT_INPUT,LTE_THERMAL_BB_INPUT,&
        DO_DEBUG_WRITE,DO_WRITE_INPUT,DO_WRITE_SCENARIO,&
        DO_WRITE_FOURIER,DO_WRITE_RESULTS,&
        INPUT_WRITE_FILENAME,SCENARIO_WRITE_FILENAME,&
        FOURIER_WRITE_FILENAME,RESULTS_WRITE_FILENAME,&
        DO_SSCORR_NADIR,DO_SSCORR_OUTGOING,DO_DOUBLE_CONVTEST,&
        DO_SOLAR_SOURCES,DO_REFRACTIVE_GEOMETRY,DO_CHAPMAN_FUNCTION,&
        DO_RAYLEIGH_ONLY,DO_DELTAM_SCALING,DO_SOLUTION_SAVING,&
        DO_BVP_TELESCOPING,DO_USER_VZANGLES,DO_ADDITIONAL_MVOUT,&
        DO_MVOUT_ONLY,DO_THERMAL_TRANSONLY,&
        NGREEK_MOMENTS_INPUT,FLUX_FACTOR,N_SZANGLES,&
        N_USER_RELAZMS,USER_RELAZMS,&
        USER_VZANGLES,GEOMETRY_SPECHEIGHT,&
        EARTH_RADIUS,OMEGA_TOTAL_INPUT)

      IF (.NOT. DO_LAMBERTIAN_SURFACE) THEN
        CALL VLIDORT_WRITE_SUP_BRDF_INPUT ( &
          EXACTDB_BRDFUNC,BRDF_F_0,BRDF_F,USER_BRDF_F_0,USER_BRDF_F,&
          EMISSIVITY,USER_EMISSIVITY)
      END IF

      IF (DO_SS_EXTERNAL) THEN
        CALL VLIDORT_WRITE_SUP_SS_INPUT ( &
          STOKES_SS,STOKES_DB)
      END IF

      IF (DO_SURFACE_LEAVING .OR. DO_SL_ISOTROPIC) THEN
        CALL VLIDORT_WRITE_SUP_SLEAVE_INPUT ( &
          SLTERM_ISOTROPIC,SLTERM_USERANGLES,SLTERM_F_0,USER_SLTERM_F_0)
      END IF

      END SUBROUTINE VLIDORT_DEBUG_INPUT_MASTER

      SUBROUTINE VLIDORT_DEBUG_LIN_INPUT_MASTER()

      CALL VLIDORT_WRITE_LIN_INPUT ( &
        DO_SIMULATION_ONLY,DO_SURFBB_LINEARIZATION,&
        LAYER_VARY_FLAG,LAYER_VARY_NUMBER,&
        N_TOTALCOLUMN_WFS,N_TOTALPROFILE_WFS,N_SURFACE_WFS,&
        COLUMNWF_NAMES,PROFILEWF_NAMES,DO_LTE_LINEARIZATION,&
        L_DELTAU_VERT_INPUT,L_OMEGA_TOTAL_INPUT,L_GREEKMAT_TOTAL_INPUT,&
        DO_COLUMN_LINEARIZATION,DO_PROFILE_LINEARIZATION,DO_ATMOS_LINEARIZATION,&
        DO_SURFACE_LINEARIZATION,DO_LINEARIZATION)

      IF (.NOT. DO_LAMBERTIAN_SURFACE) THEN
        CALL VLIDORT_WRITE_LIN_SUP_BRDF_INPUT ( &
          LS_EXACTDB_BRDFUNC,LS_BRDF_F_0,LS_BRDF_F,LS_USER_BRDF_F_0,LS_USER_BRDF_F,&
          LS_EMISSIVITY,LS_USER_EMISSIVITY)
      END IF

      IF (DO_SS_EXTERNAL) THEN
        CALL VLIDORT_WRITE_LIN_SUP_SS_INPUT ( &
          COLUMNWF_SS,COLUMNWF_DB,PROFILEWF_SS,PROFILEWF_DB,SURFACEWF_DB)
      END IF

      END SUBROUTINE VLIDORT_DEBUG_LIN_INPUT_MASTER

      END SUBROUTINE VLIDORT_L_MASTER

!

      SUBROUTINE VLIDORT_L_FOURIER ( &
        FOURIER_COMPONENT, SS_FLUX_MULTIPLIER, DO_SSCORR_NADIR, DO_SSCORR_OUTGOING, &
        DO_SSCORR_TRUNCATION, DO_SSFULL, DO_SOLAR_SOURCES, DO_PLANE_PARALLEL, &
        DO_REFRACTIVE_GEOMETRY, DO_DELTAM_SCALING, DO_SOLUTION_SAVING, DO_BVP_TELESCOPING, &
        DO_UPWELLING, DO_DNWELLING, DO_QUAD_OUTPUT, DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY, &
        DO_DEBUG_WRITE, DO_FDTEST, NSTOKES, NSTREAMS, NLAYERS, NGREEK_MOMENTS_INPUT, &
        FLUX_FACTOR, SZA_LOCAL_INPUT, N_USER_RELAZMS, USER_RELAZMS, N_USER_LEVELS, &
        NFINELAYERS, EARTH_RADIUS, HEIGHT_GRID, &
        OMEGA_TOTAL_INPUT, DELTAU_VERT_INPUT, GREEKMAT_TOTAL_INPUT, &
        DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
        EXACTDB_BRDFUNC, BRDF_F, BRDF_F_0, USER_BRDF_F, USER_BRDF_F_0, &
        EMISSIVITY, USER_EMISSIVITY, DO_THERMAL_EMISSION, N_THERMAL_COEFFS, THERMAL_BB_INPUT, &
        DO_SURFACE_EMISSION, SURFBB, DO_THERMAL_TRANSONLY, DO_SPECIALIST_OPTION_2, DO_SPECIALIST_OPTION_3, &
        DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, SLTERM_ISOTROPIC, SLTERM_USERANGLES, SLTERM_F_0, USER_SLTERM_F_0, &
        NLAYERS_CUTOFF, DO_TOA_CONTRIBS, DO_DIRECT_BEAM, &
        DO_CLASSICAL_SOLUTION, DO_DBCORRECTION, FLUXVEC, &
        COS_SZANGLES, SIN_SZANGLES, SUN_SZA_COSINES, DO_MULTIBEAM, &
        QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS, QUAD_HALFWTS, &
        DO_MSMODE_VLIDORT, DO_MSMODE_THERMAL, NMOMENTS, NSTREAMS_2, NTOTAL, &
        N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2, NBEAMS, N_USER_STREAMS, &
        MUELLER_INDEX, DMAT, DO_REAL_EIGENSOLVER, BVP_REGULAR_FLAG, LAYER_MAXMOMENTS, &
        DO_LAYER_SCATTERING, DO_USER_STREAMS, LOCAL_UM_START, N_DIRECTIONS, WHICH_DIRECTIONS, &
        USER_VZANGLES_ADJUST, SZANGLES_ADJUST, USER_RELAZMS_ADJUST, USER_STREAMS, USER_SECANTS, &
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, &
        DO_PARTLAYERS, N_PARTLAYERS, PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES, &
        N_ALLLAYERS_UP, N_ALLLAYERS_DN, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, &
        DFLUX, N_GEOMETRIES, VZA_OFFSETS, TAUGRID_INPUT, CHAPMAN_FACTORS, &
        NSTOKES_SQ, DO_SIMULATION_ONLY, DO_PROFILE_LINEARIZATION, DO_COLUMN_LINEARIZATION, &
        DO_ATMOS_LINEARIZATION, DO_LTE_LINEARIZATION, N_TOTALCOLUMN_WFS, &
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER, L_OMEGA_TOTAL_INPUT, L_DELTAU_VERT_INPUT, &
        L_GREEKMAT_TOTAL_INPUT, LTE_DELTAU_VERT_INPUT, LTE_THERMAL_BB_INPUT, &
        DO_SURFACE_LINEARIZATION, N_SURFACE_WFS, LS_EXACTDB_BRDFUNC, LS_BRDF_F, LS_BRDF_F_0, &
        LS_USER_BRDF_F, LS_USER_BRDF_F_0, LS_USER_EMISSIVITY, LS_EMISSIVITY, &
        TAUGRID, OMEGA_TOTAL, GREEKMAT_TOTAL, &
        STOKES_SS, STOKES_DB, SS_CONTRIBS, MEAN_STOKES, FLUX_STOKES, MEAN_DIRECT, FLUX_DIRECT, &
        MS_CONTRIBS_F, STOKES_F, DO_INCLUDE_SURFACE, DO_INCLUDE_SURFEMISS, DO_INCLUDE_THERMEMISS, &
        PROFILEWF_SS, COLUMNWF_SS, PROFILEWF_DB, COLUMNWF_DB, SURFACEWF_DB, &
        ATMOSWF_F, MINT_ATMOSWF, MINT_ATMOSWF_DIRECT, FLUX_ATMOSWF, FLUX_ATMOSWF_DIRECT, &
        SURFACEWF_F, MINT_SURFACEWF, FLUX_SURFACEWF, LTE_ATMOSWF, &
        STATUS, MESSAGE, TRACE_1, TRACE_2 )

!  Complete Fourier component calculation for the Extended Code
!    Stokes vector and Weighting function computations

      USE VLIDORT_PARS

      USE VLIDORT_MISCSETUPS_MODULE
      USE VLIDORT_THERMALSUP
      USE VLIDORT_MULTIPLIERS
      USE VLIDORT_CORRECTIONS
      USE VLIDORT_SOLUTIONS
      USE VLIDORT_BVPROBLEM
      USE VLIDORT_THERMALSUP
      USE VLIDORT_INTENSITY

      USE VLIDORT_L_MISCSETUPS_MODULE
      USE VLIDORT_L_THERMALSUP
      USE VLIDORT_L_MULTIPLIERS
      USE VLIDORT_L_CORRECTIONS
      USE VLIDORT_L_SOLUTIONS
      USE VLIDORT_L_BVPROBLEM
      USE VLIDORT_L_THERMALSUP
      USE VLIDORT_L_WFATMOS
      USE VLIDORT_L_WFSURFACE

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          FOURIER_COMPONENT
      DOUBLE PRECISION, INTENT (IN) :: SS_FLUX_MULTIPLIER
      LOGICAL, INTENT (IN) ::          DO_SSCORR_NADIR
      LOGICAL, INTENT (IN) ::          DO_SSCORR_OUTGOING
      LOGICAL, INTENT (IN) ::          DO_SSCORR_TRUNCATION
      LOGICAL, INTENT (IN) ::          DO_SSFULL
      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL
      LOGICAL, INTENT (IN) ::          DO_REFRACTIVE_GEOMETRY
      LOGICAL, INTENT (IN) ::          DO_DELTAM_SCALING
      LOGICAL, INTENT (IN) ::          DO_SOLUTION_SAVING
      LOGICAL, INTENT (IN) ::          DO_BVP_TELESCOPING
      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_DNWELLING
      LOGICAL, INTENT (IN) ::          DO_QUAD_OUTPUT
      LOGICAL, INTENT (IN) ::          DO_ADDITIONAL_MVOUT
      LOGICAL, INTENT (IN) ::          DO_MVOUT_ONLY
      LOGICAL, INTENT (IN) ::          DO_DEBUG_WRITE
      LOGICAL, INTENT (IN) ::          DO_FDTEST
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NGREEK_MOMENTS_INPUT
      DOUBLE PRECISION, INTENT (IN) :: FLUX_FACTOR
      DOUBLE PRECISION, INTENT (IN) :: SZA_LOCAL_INPUT &
          ( 0:MAXLAYERS, MAX_SZANGLES )
      INTEGER, INTENT (IN) ::          N_USER_RELAZMS
      DOUBLE PRECISION, INTENT (IN) :: USER_RELAZMS  ( MAX_USER_RELAZMS )
      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      INTEGER, INTENT (IN) ::          NFINELAYERS
      DOUBLE PRECISION, INTENT (IN) :: EARTH_RADIUS
      DOUBLE PRECISION, INTENT (IN) :: HEIGHT_GRID ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: OMEGA_TOTAL_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: GREEKMAT_TOTAL_INPUT &
          ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      DOUBLE PRECISION, INTENT (IN) :: LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: EXACTDB_BRDFUNC &
          ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: BRDF_F &
          ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: BRDF_F_0 &
          ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS   )
      DOUBLE PRECISION, INTENT (IN) :: USER_BRDF_F &
          ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: USER_BRDF_F_0 &
          ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS   )
      DOUBLE PRECISION, INTENT (IN) :: EMISSIVITY ( MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: USER_EMISSIVITY &
          ( MAXSTOKES, MAX_USER_STREAMS )
      LOGICAL, INTENT (IN) ::          DO_THERMAL_EMISSION
      INTEGER, INTENT (IN) ::          N_THERMAL_COEFFS
      DOUBLE PRECISION, INTENT (IN) :: THERMAL_BB_INPUT ( 0:MAXLAYERS )
      LOGICAL, INTENT (IN) ::          DO_SURFACE_EMISSION
      DOUBLE PRECISION, INTENT (IN) :: SURFBB
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_SPECIALIST_OPTION_2
      LOGICAL, INTENT (IN) ::          DO_SPECIALIST_OPTION_3

!  New Surface-Leaving stuff 17 May 2012
      LOGICAL, INTENT (IN) ::            DO_SURFACE_LEAVING
      LOGICAL, INTENT (IN) ::            DO_SL_ISOTROPIC
      DOUBLE PRECISION, INTENT (IN) ::   SLTERM_ISOTROPIC ( MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::   SLTERM_USERANGLES &
          ( MAXSTOKES, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS  )
      DOUBLE PRECISION, INTENT (IN) ::   SLTERM_F_0 &
          ( 0:MAXMOMENTS, MAXSTOKES, MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::   USER_SLTERM_F_0 &
          ( 0:MAXMOMENTS, MAXSTOKES, MAX_USER_STREAMS, MAXBEAMS )

      INTEGER, INTENT (IN) ::          NLAYERS_CUTOFF
      LOGICAL, INTENT (IN) ::          DO_TOA_CONTRIBS
      LOGICAL, INTENT (IN) ::          DO_DIRECT_BEAM
      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION
      LOGICAL, INTENT (IN) ::          DO_DBCORRECTION
      DOUBLE PRECISION, INTENT (IN) :: FLUXVEC ( MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: COS_SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) :: SIN_SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) :: SUN_SZA_COSINES &
          ( MAXLAYERS, MAX_SZANGLES )
      LOGICAL, INTENT (IN) ::          DO_MULTIBEAM ( MAXBEAMS, 0:MAXFOURIER )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_HALFWTS ( MAXSTREAMS )
      LOGICAL, INTENT (IN) ::          DO_MSMODE_VLIDORT
      LOGICAL, INTENT (IN) ::          DO_MSMODE_THERMAL
      INTEGER, INTENT (IN) ::          NMOMENTS
      INTEGER, INTENT (IN) ::          NSTREAMS_2
      INTEGER, INTENT (IN) ::          NTOTAL
      INTEGER, INTENT (IN) ::          N_SUBDIAG
      INTEGER, INTENT (IN) ::          N_SUPDIAG
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS_2
      INTEGER, INTENT (IN) ::          NBEAMS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: DMAT ( MAXSTOKES, MAXSTOKES )
      LOGICAL, INTENT (IN) ::          DO_REAL_EIGENSOLVER &
          ( 0:MAXMOMENTS, MAXLAYERS )
      LOGICAL, INTENT (IN) ::          BVP_REGULAR_FLAG ( 0:MAXMOMENTS )
      INTEGER, INTENT (IN) ::          LAYER_MAXMOMENTS ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::          DO_LAYER_SCATTERING &
          ( 0:MAXMOMENTS, MAXLAYERS )
      LOGICAL, INTENT (IN) ::          DO_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      INTEGER, INTENT (IN) ::          N_DIRECTIONS
      INTEGER, INTENT (IN) ::          WHICH_DIRECTIONS ( MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (IN) :: USER_VZANGLES_ADJUST &
          ( MAX_USER_VZANGLES )
      DOUBLE PRECISION, INTENT (IN) :: SZANGLES_ADJUST &
          ( MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS )
      DOUBLE PRECISION, INTENT (IN) :: USER_RELAZMS_ADJUST &
          ( MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS )
      DOUBLE PRECISION, INTENT (IN) :: USER_STREAMS  ( MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: USER_SECANTS  ( MAX_USER_STREAMS )
      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )
      LOGICAL, INTENT (IN) ::          DO_PARTLAYERS
      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PARTLAYERS_VALUES ( MAX_PARTLAYERS )
      INTEGER, INTENT (IN) ::          N_ALLLAYERS_UP
      INTEGER, INTENT (IN) ::          N_ALLLAYERS_DN
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_DN ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DFLUX ( MAXSTOKES )
      INTEGER, INTENT (IN) ::          N_GEOMETRIES
      INTEGER, INTENT (IN) ::          VZA_OFFSETS &
          ( MAX_SZANGLES, MAX_USER_VZANGLES )
      DOUBLE PRECISION, INTENT (IN) :: TAUGRID_INPUT ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: CHAPMAN_FACTORS &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

      INTEGER, INTENT (IN) ::          NSTOKES_SQ
      LOGICAL, INTENT (IN) ::          DO_SIMULATION_ONLY
      LOGICAL, INTENT (IN) ::          DO_PROFILE_LINEARIZATION
      LOGICAL, INTENT (IN) ::          DO_COLUMN_LINEARIZATION
      LOGICAL, INTENT (IN) ::          DO_ATMOS_LINEARIZATION
      LOGICAL, INTENT (IN) ::          DO_LTE_LINEARIZATION
      INTEGER, INTENT (IN) ::          N_TOTALCOLUMN_WFS
      LOGICAL, INTENT (IN) ::          LAYER_VARY_FLAG  ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          LAYER_VARY_NUMBER ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_TOTAL_INPUT &
          ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT_INPUT &
          ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_GREEKMAT_TOTAL_INPUT &
           ( MAX_ATMOSWFS, 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION, INTENT (IN) :: LTE_DELTAU_VERT_INPUT ( 2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LTE_THERMAL_BB_INPUT ( 0:MAXLAYERS )
      LOGICAL, INTENT (IN) ::          DO_SURFACE_LINEARIZATION
      INTEGER, INTENT (IN) ::          N_SURFACE_WFS
      DOUBLE PRECISION, INTENT (IN) :: LS_EXACTDB_BRDFUNC &
          ( MAX_SURFACEWFS, MAXSTOKES_SQ, &
            MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: LS_BRDF_F &
          ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTOKES_SQ, &
            MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: LS_BRDF_F_0 &
          ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTOKES_SQ, &
            MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: LS_USER_BRDF_F &
          ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTOKES_SQ, &
            MAX_USER_STREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: LS_USER_BRDF_F_0 &
          ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTOKES_SQ, &
            MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: LS_USER_EMISSIVITY &
          ( MAX_SURFACEWFS, MAXSTOKES,MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: LS_EMISSIVITY &
          ( MAX_SURFACEWFS, MAXSTOKES,MAXSTREAMS )

      DOUBLE PRECISION, INTENT (INOUT) ::  TAUGRID ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) ::  OMEGA_TOTAL ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) ::  GREEKMAT_TOTAL &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION, INTENT (INOUT) ::  STOKES_SS &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) ::  STOKES_DB &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) ::  SS_CONTRIBS &
          ( MAX_GEOMETRIES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) ::  MEAN_STOKES &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) ::  FLUX_STOKES &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES ,MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) ::  MEAN_DIRECT &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) ::  FLUX_DIRECT &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES )

      DOUBLE PRECISION, INTENT (OUT) ::  MS_CONTRIBS_F &
          ( MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAXLAYERS  )
      DOUBLE PRECISION, INTENT (OUT) ::  STOKES_F &
          ( MAX_USER_LEVELS, MAX_USER_VZANGLES, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      LOGICAL, INTENT (OUT) ::           DO_INCLUDE_SURFACE
      LOGICAL, INTENT (OUT) ::           DO_INCLUDE_SURFEMISS
      LOGICAL, INTENT (OUT) ::           DO_INCLUDE_THERMEMISS

      DOUBLE PRECISION, INTENT (INOUT) :: PROFILEWF_SS &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
            MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) :: COLUMNWF_SS &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, &
            MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) :: PROFILEWF_DB &
          ( MAX_ATMOSWFS, MAXLAYERS, &
            MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) :: COLUMNWF_DB &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) :: SURFACEWF_DB &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

      DOUBLE PRECISION, INTENT (OUT) ::  ATMOSWF_F &
          ( MAX_ATMOSWFS, 0:MAXLAYERS, MAX_USER_LEVELS, &
            MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

      DOUBLE PRECISION, INTENT (INOUT) ::  MINT_ATMOSWF &
          ( MAX_ATMOSWFS, 0:MAXLAYERS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) ::  MINT_ATMOSWF_DIRECT &
          ( MAX_ATMOSWFS, 0:MAXLAYERS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) ::  FLUX_ATMOSWF &
          ( MAX_ATMOSWFS, 0:MAXLAYERS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) ::  FLUX_ATMOSWF_DIRECT &
          ( MAX_ATMOSWFS, 0:MAXLAYERS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES )

      DOUBLE PRECISION, INTENT (OUT) ::  SURFACEWF_F &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_VZANGLES, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

      DOUBLE PRECISION, INTENT (INOUT) ::  MINT_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) ::  FLUX_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) ::  LTE_ATMOSWF &
          ( 0:MAXLAYERS, MAX_USER_LEVELS, &
            MAX_USER_VZANGLES, MAX_DIRECTIONS )

      INTEGER, INTENT (OUT) ::           STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE_1
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE_2

!  Local variables
!  ---------------

      INTEGER ::          LAYER, IBEAM, IPARTIC, I, O1, N

!  local inclusion flags

      LOGICAL ::          DO_INCLUDE_MVOUTPUT
      LOGICAL ::          DO_INCLUDE_DIRECTBEAM

!  Flux multiplier and Fourier component numbers

      DOUBLE PRECISION :: FLUX_MULTIPLIER
      DOUBLE PRECISION :: DELTA_FACTOR
      DOUBLE PRECISION :: SURFACE_FACTOR

!  weighting function indices

      INTEGER ::          N_LAYER_WFS, LAYER_TO_VARY, N_PARAMETERS
      LOGICAL ::          DO_RTSOL_VARY
      INTEGER ::          NPARAMS_VARY, VARIATION_INDEX

!  error tracing variables

      CHARACTER (LEN=2) :: CF
      INTEGER ::           STATUS_SUB
      LOGICAL ::           FAIL

!  progress

      LOGICAL, PARAMETER :: DO_WRITE_SCREEN = .FALSE.

!  helper

      DOUBLE PRECISION :: WUPPER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION :: WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION :: T_UTUP_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION :: T_UTDN_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )

      DOUBLE PRECISION :: PI_XQP &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PI_XQM &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PI_XUP &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PI_XUM &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PI_X0P &
          ( 0:MAXMOMENTS, MAXBEAMS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PI_XQM_POST &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PI_XQM_PRE &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PI_XQP_PRE &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PI_XUM_POST &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PI_XUP_PRE &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )

      DOUBLE PRECISION :: T_WUPPER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION :: T_WLOWER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION :: UT_T_PARTIC ( MAXSTREAMS_2, MAX_PARTLAYERS )
      DOUBLE PRECISION :: U_TPOS1 ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: U_TNEG1 ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: U_TPOS2 ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: U_TNEG2 ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: LAYER_TSUP_UP ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: LAYER_TSUP_DN ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: LAYER_TSUP_UTUP ( MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION :: LAYER_TSUP_UTDN ( MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION :: SAB &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION :: DAB &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION :: EIGENMAT_SAVE ( MAXEVALUES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION :: REAL_KSQ ( MAXSTRMSTKS )
      DOUBLE PRECISION :: IMAG_KSQ ( MAXSTRMSTKS )
      DOUBLE PRECISION :: LEFT_EVEC ( MAXSTRMSTKS, MAXSTRMSTKS )
      DOUBLE PRECISION :: RITE_EVEC ( MAXSTRMSTKS, MAXSTRMSTKS )
      LOGICAL ::          EIGENDEGEN ( MAXSTRMSTKS, MAXLAYERS )
      INTEGER ::          EIGENMASK_R ( MAXEVALUES )
      INTEGER ::          EIGENMASK_C ( MAXEVALUES )
      DOUBLE PRECISION :: FWD_SUMVEC &
          ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION :: FWD_DIFVEC &
          ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      INTEGER ::          K_REAL ( MAXLAYERS )
      INTEGER ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION :: KEIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION :: KEIGEN_CSQ ( MAXEVALUES )
      DOUBLE PRECISION :: SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION :: SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION :: QSUMVEC_SAVE ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION :: QDIFVEC_SAVE ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION :: QVEC_SAVE ( MAXSTRMSTKS )
      DOUBLE PRECISION :: QDIF_SAVE ( MAXSTRMSTKS )
      DOUBLE PRECISION :: QMAT_SAVE ( MAXSTRMSTKS, MAXSTRMSTKS )
      INTEGER ::          QPIVOT ( MAXSTRMSTKS )
      DOUBLE PRECISION :: BVEC ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION :: BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER ::          IPIVOT ( MAXTOTAL )
      DOUBLE PRECISION :: SMAT2 ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )
      INTEGER ::          SIPIVOT ( MAXSTRMSTKS_2 )
      DOUBLE PRECISION :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION :: MCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION :: BANDTELMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER ::          IPIVOTTEL ( MAXTOTAL )
      DOUBLE PRECISION :: UHOM_DNDN &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION :: UHOM_DNUP &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION :: UHOM_UPDN &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION :: UHOM_UPUP &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION :: BOA_THTONLY_SOURCE ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION :: DIRECT_BEAM &
          ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION :: USER_DIRECT_BEAM &
          ( MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )
      LOGICAL ::          HSINGO (MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: HMULT_1 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: HMULT_2 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: UT_HMULT_UU &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION :: UT_HMULT_UD &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION :: UT_HMULT_DU &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION :: UT_HMULT_DD &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION :: AXBID_F &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES_SQ )
      DOUBLE PRECISION :: HELPSTOKES &
          ( 0:MAXMOMENTS, MAXEVALUES, MAXSTOKES )
      DOUBLE PRECISION :: ZETA_M &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: ZETA_P &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      INTEGER ::          ACTIVE_LAYERS ( MAXLAYERS )
      DOUBLE PRECISION :: ATMOS_ATTN ( MAXBEAMS )
      INTEGER ::          N_BVTELMATRIX_SIZE
      INTEGER ::          N_BVTELMATRIX_SUPDIAG
      INTEGER ::          N_BVTELMATRIX_SUBDIAG
      INTEGER ::          NLAYERS_TEL
      DOUBLE PRECISION :: R2_HOMP ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION :: R2_HOMM ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION :: R2_BEAM ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION :: HELPSTOKES_BEAM ( 0:MAXMOMENTS, MAXSTOKES )
      DOUBLE PRECISION :: UPAR_DN_1 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION :: UPAR_DN_2 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION :: UPAR_UP_1 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION :: UPAR_UP_2 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION :: STOKES_DOWNSURF ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION :: CUMSOURCE_UP &
          ( MAX_USER_STREAMS, MAXSTOKES, 0:MAXLAYERS )
      DOUBLE PRECISION :: CUMSOURCE_DN &
          ( MAX_USER_STREAMS, MAXSTOKES, 0:MAXLAYERS )

      DOUBLE PRECISION :: L_T_WUPPER &
          ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_WLOWER &
          ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS)
      DOUBLE PRECISION :: L_UT_T_PARTIC &
          ( MAXSTREAMS_2, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_U_TPOS1 &
          ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_U_TNEG1 &
          ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_U_TPOS2 &
          ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_U_TNEG2 &
          ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_LAYER_TSUP_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_LAYER_TSUP_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_LAYER_TSUP_UTUP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_LAYER_TSUP_UTDN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_UHOM_DNDN ( MAX_USER_STREAMS, &
          MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_UHOM_DNUP ( MAX_USER_STREAMS, &
          MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_UHOM_UPDN ( MAX_USER_STREAMS, &
          MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_UHOM_UPUP ( MAX_USER_STREAMS, &
          MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_BOA_THTONLY_SOURCE (MAXSTREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_HMULT_1 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_HMULT_2 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_UT_HMULT_UU &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_UT_HMULT_UD &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_UT_HMULT_DU &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_UT_HMULT_DD &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_ZETA_M &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS,  MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_ZETA_P &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS,  MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_SAB &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, &
            MAXSTOKES,  MAXLAYERS,  MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_DAB &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, &
            MAXSTOKES,  MAXLAYERS,  MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_EIGENMAT &
          ( MAXEVALUES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_KEIGEN &
          ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_DELT_EIGEN &
          ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_UTUP_EIGEN &
          ( MAXEVALUES, MAX_USER_LEVELS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_UTDN_EIGEN &
          ( MAXEVALUES, MAX_USER_LEVELS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_BVEC &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, 0:MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_UPAR_DN_1 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_UPAR_UP_1 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_UPAR_DN_2 ( MAX_USER_STREAMS, &
          MAXSTOKES, MAXLAYERS, 0:MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_UPAR_UP_2 ( MAX_USER_STREAMS, &
          MAXSTOKES, MAXLAYERS, 0:MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_WUPPER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_WLOWER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: NCON &
          ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: PCON &
          ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )

!  Saved local variables
!  ---------------------

      LOGICAL, SAVE ::          DO_BVTEL_INITIAL

!FROM VLIDORT_MISCSETUPS:
      DOUBLE PRECISION, SAVE :: DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, SAVE :: PARTAU_VERT ( MAX_PARTLAYERS )
      DOUBLE PRECISION, SAVE :: DELTAU_SLANT &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, SAVE :: TRUNC_FACTOR ( MAXLAYERS )
      DOUBLE PRECISION, SAVE :: FAC1 ( MAXLAYERS )
      DOUBLE PRECISION, SAVE :: OMEGA_GREEK &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      INTEGER, SAVE ::          LAYER_PIS_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, SAVE :: SOLAR_BEAM_OPDEP ( MAXBEAMS )
      LOGICAL, SAVE ::          DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )
      DOUBLE PRECISION, SAVE :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, SAVE :: T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, SAVE :: T_DISORDS_UTDN ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, SAVE :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, SAVE :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, SAVE :: T_UTUP_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, SAVE :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, SAVE :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, SAVE :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, SAVE :: CUMTRANS ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, SAVE :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, SAVE :: ITRANS_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, SAVE :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, SAVE :: LOCAL_CSZA ( 0:MAXLAYERS, MAXBEAMS )
!FROM THERMAL_SETUP:
      DOUBLE PRECISION, SAVE :: THERMCOEFFS ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, SAVE :: DELTAU_POWER ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, SAVE :: XTAU_POWER ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, SAVE :: TCOM1 ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, SAVE :: T_DIRECT_UP ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, SAVE :: T_DIRECT_DN ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, SAVE :: T_UT_DIRECT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, SAVE :: T_UT_DIRECT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS )
!FROM EMULT_MASTER:
      LOGICAL, SAVE ::          EMULT_HOPRULE &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, SAVE :: SIGMA_M &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, SAVE :: SIGMA_P &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, SAVE :: EMULT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, SAVE :: EMULT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, SAVE :: UT_EMULT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, SAVE :: UT_EMULT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )

!FROM VLIDORT_L_MISCSETUPS:
      DOUBLE PRECISION, SAVE :: L_OMEGA_TOTAL &
          ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, SAVE :: L_DELTAU_VERT &
          ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, SAVE :: L_GREEKMAT_TOTAL &
          ( MAX_ATMOSWFS, 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION, SAVE :: L_DELTAU_SLANT &
          ( MAX_ATMOSWFS, MAXLAYERS, MAXLAYERS, MAXBEAMS )
      LOGICAL, SAVE ::          DO_SCATMAT_VARIATION &
          ( MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, SAVE :: L_TRUNC_FACTOR &
          ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, SAVE :: L_OMEGA_GREEK &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION, SAVE :: L_AVERAGE_SECANT &
          ( MAXLAYERS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, SAVE :: L_INITIAL_TRANS &
          ( MAXLAYERS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, SAVE :: L_T_DELT_DISORDS &
          ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, SAVE :: L_T_DISORDS_UTDN &
          ( MAXSTREAMS, MAX_USER_LEVELS, MAX_ATMOSWFS )
      DOUBLE PRECISION, SAVE :: L_T_DISORDS_UTUP &
          ( MAXSTREAMS, MAX_USER_LEVELS, MAX_ATMOSWFS )
      DOUBLE PRECISION, SAVE :: L_T_DELT_MUBAR &
          ( MAXLAYERS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, SAVE :: L_T_UTDN_MUBAR &
          ( MAX_USER_LEVELS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)
      DOUBLE PRECISION, SAVE :: L_T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, SAVE :: L_T_UTDN_USERM &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, SAVE :: L_T_UTUP_USERM &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS )
!FROM THERMAL_SETUP_PLUS:
      DOUBLE PRECISION, SAVE :: L_THERMCOEFFS &
          ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION, SAVE :: L_DELTAU_POWER &
          ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION, SAVE :: L_XTAU_POWER &
          ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION, SAVE :: L_TCOM1 &
          ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION, SAVE :: L_T_DIRECT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, SAVE :: L_T_DIRECT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, SAVE :: L_T_UT_DIRECT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS)
      DOUBLE PRECISION, SAVE :: L_T_UT_DIRECT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
!FROM L_EMULT_MASTER:
      DOUBLE PRECISION, SAVE :: L_EMULT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, 0:MAXLAYERS, &
            MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, SAVE :: L_EMULT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, 0:MAXLAYERS, &
            MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, SAVE :: L_UT_EMULT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, 0:MAXLAYERS, &
            MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, SAVE :: L_UT_EMULT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, 0:MAXLAYERS, &
            MAXBEAMS, MAX_ATMOSWFS )

!FROM VLIDORT_SSCORR_NADIR:
      DOUBLE PRECISION, SAVE :: ZMAT_UP &
          ( MAX_GEOMETRIES, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, SAVE :: ZMAT_DN &
          ( MAX_GEOMETRIES, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, SAVE :: TMS ( MAXLAYERS )
      DOUBLE PRECISION, SAVE :: SSFDEL ( MAXLAYERS )
      DOUBLE PRECISION, SAVE :: SS_CUMSOURCE_UP &
          ( MAX_GEOMETRIES, MAXSTOKES, 0:MAXLAYERS )
      DOUBLE PRECISION, SAVE :: SS_CUMSOURCE_DN &
          ( MAX_GEOMETRIES, MAXSTOKES, 0:MAXLAYERS )
!FROM VLIDORT_SSCORR_OUTGOING:
      DOUBLE PRECISION, SAVE :: UP_MULTIPLIERS &
          ( MAXLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION, SAVE :: DN_MULTIPLIERS &
          ( MAXLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION, SAVE :: UP_MULTIPLIERS_UT &
          ( MAX_PARTLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION, SAVE :: DN_MULTIPLIERS_UT &
          ( MAX_PARTLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION, SAVE :: UP_LOSTRANS &
          ( MAXLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION, SAVE :: DN_LOSTRANS &
          ( MAXLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION, SAVE :: UP_LOSTRANS_UT &
          ( MAX_PARTLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION, SAVE :: DN_LOSTRANS_UT &
          ( MAX_PARTLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION, SAVE :: BOA_ATTN ( MAX_GEOMETRIES )
!FROM VLIDORT_DBCORRECTION:
      DOUBLE PRECISION, SAVE :: DB_CUMSOURCE &
          ( MAX_GEOMETRIES, MAXSTOKES, 0:MAXLAYERS )
      DOUBLE PRECISION, SAVE :: EXACTDB_SOURCE &
          ( MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, SAVE :: ATTN_DB_SAVE ( MAX_GEOMETRIES )

!FROM VLIDORT_LP_SSCORR_NADIR:
!      DOUBLE PRECISION, SAVE :: L_ZMAT_UP &
!          ( MAX_ATMOSWFS, MAX_GEOMETRIES, MAXLAYERS, MAXSTOKES, MAXSTOKES )
!      DOUBLE PRECISION, SAVE :: L_ZMAT_DN &
!          ( MAX_ATMOSWFS, MAX_GEOMETRIES, MAXLAYERS, MAXSTOKES, MAXSTOKES )
!      DOUBLE PRECISION, SAVE :: L_SS_CUMSOURCE &
!          ( MAX_ATMOSWFS, MAX_GEOMETRIES, MAXSTOKES )
!      DOUBLE PRECISION, SAVE :: L_SSFDEL ( MAXLAYERS, MAX_ATMOSWFS )
!FROM VLIDORT_L_SSCORR_OUTGOING:
      DOUBLE PRECISION, SAVE :: L_UP_MULTIPLIERS &
          ( MAXLAYERS, 0:MAXLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )
      DOUBLE PRECISION, SAVE :: L_DN_MULTIPLIERS &
          ( MAXLAYERS, 0:MAXLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )
      DOUBLE PRECISION, SAVE :: L_UP_LOSTRANS &
          ( MAXLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )
      DOUBLE PRECISION, SAVE :: L_DN_LOSTRANS &
          ( MAXLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )
      DOUBLE PRECISION, SAVE :: L_UP_MULTIPLIERS_UT &
          ( MAX_PARTLAYERS, 0:MAXLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )
      DOUBLE PRECISION, SAVE :: L_DN_MULTIPLIERS_UT &
          ( MAX_PARTLAYERS, 0:MAXLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )
      DOUBLE PRECISION, SAVE :: L_UP_LOSTRANS_UT &
          ( MAX_PARTLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )
      DOUBLE PRECISION, SAVE :: L_DN_LOSTRANS_UT &
          ( MAX_PARTLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )
      DOUBLE PRECISION, SAVE :: L_BOA_ATTN &
          ( 0:MAXLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )
!FROM VLIDORT_LAP_DBCORRECTION:
      DOUBLE PRECISION, SAVE :: L_EXACTDB_SOURCE &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, SAVE :: L_DB_CUMSOURCE &
          ( MAX_GEOMETRIES, MAXSTOKES )
!FROM VLIDORT_LAC_DBCORRECTION:
      DOUBLE PRECISION, SAVE :: L_EXACTDBC_SOURCE &
          ( MAX_ATMOSWFS, MAX_GEOMETRIES, MAXSTOKES )
!FROM VLIDORT_LS_DBCORRECTION:
      DOUBLE PRECISION, SAVE :: LS_EXACTDB_SOURCE &
          ( MAX_SURFACEWFS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, SAVE :: LS_DB_CUMSOURCE &
          ( MAX_GEOMETRIES, MAXSTOKES )

!TESTING
      INTEGER :: Q,UM,LVARY,K
      LOGICAL :: DO_DEBUG

!  ##############
!  initialization
!  ##############

!  module status and message initialization

      STATUS  = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE_1 = ' '
      TRACE_2 = ' '

!  Set local flags
!  ---------------

!  inclusion of thermal surface emission term, only for Fourier = 0

      DO_INCLUDE_SURFEMISS = .FALSE.
      IF ( DO_SURFACE_EMISSION ) THEN
        IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
          DO_INCLUDE_SURFEMISS = .TRUE.
        ENDIF
      ENDIF

!  inclusion of thermal emission term, only for Fourier = 0

      DO_INCLUDE_THERMEMISS = .FALSE.
      IF ( DO_THERMAL_EMISSION ) THEN
        IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
          DO_INCLUDE_THERMEMISS = .TRUE.
        ENDIF
      ENDIF

!  Surface flag (for inclusion of some kind of reflecting boundary)

      DO_INCLUDE_SURFACE = .TRUE.
      IF ( DO_LAMBERTIAN_SURFACE ) THEN
        IF ( FOURIER_COMPONENT .NE. 0 ) THEN
            DO_INCLUDE_SURFACE = .FALSE.
!mick fix 2/15/12 - lower portion of IF block commented out.
!                   "DO_INCLUDE_SURFACE = .TRUE." needed when doing
!                   analytic surfacewf even when LAMBERTIAN_ALBEDO = ZERO
        !ELSE
        !  IF ( LAMBERTIAN_ALBEDO .EQ. ZERO ) THEN
        !    DO_INCLUDE_SURFACE = .FALSE.
        !  ENDIF
        ENDIF
      ENDIF

!  Direct beam flag (only if above albedo flag has been set)
!   (DO_REFLECTED_DIRECTBEAM is a stored variable)

      IF ( DO_DIRECT_BEAM ) THEN
        IF ( DO_INCLUDE_SURFACE ) THEN
          DO IBEAM = 1, NBEAMS
            DO_REFLECTED_DIRECTBEAM(IBEAM) = .TRUE.
          ENDDO
        ELSE
          DO IBEAM = 1, NBEAMS
            DO_REFLECTED_DIRECTBEAM(IBEAM) = .FALSE.
          ENDDO
        ENDIF
      ENDIF

!  surface reflectance factors

      IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
        SURFACE_FACTOR = TWO
        DELTA_FACTOR   = ONE
      ELSE
        SURFACE_FACTOR = ONE
        DELTA_FACTOR   = TWO
      ENDIF

!  Flux multipliers

!       = 1 / 4.pi with beam sources,  = 1 for Thermal alone.    ! OLD
!       = 1                                                      ! NEW

!      FLUX_MULTIPLIER   = SS_FLUX_MULTIPLIER * DELTA_FACTOR     ! OLD
      FLUX_MULTIPLIER = DELTA_FACTOR                             ! NEW

      IF ( DO_INCLUDE_THERMEMISS.AND..NOT.DO_SOLAR_SOURCES) THEN
        FLUX_MULTIPLIER = DELTA_FACTOR
      ENDIF

!  inclusion of mean value output

      DO_INCLUDE_MVOUTPUT = .FALSE.
      IF ( DO_ADDITIONAL_MVOUT .OR. DO_MVOUT_ONLY ) THEN
        IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
          DO_INCLUDE_MVOUTPUT = .TRUE.
        ENDIF
      ENDIF

!  Initialise BVP telescoping. Important.

      IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
        DO_BVTEL_INITIAL = DO_BVP_TELESCOPING
      ENDIF

!  #################
!  Set up operations
!  #################

!  Miscellaneous setups for Fourier = 0
!   MISCSETUPS   : Delta-M, average-secant formulation, transmittances
!   EMULT_MASTER : Beam source function multipliers. Not required for the
!                  Full SS calculation in outgoing mode

!  ############################ IMPORTANT NOTE, READ CAREFULLY #########
!
!  30 January 2008
!  Telescoping setup is done in DERIVE_INPUTS (unlike LIDORT scalar code)
!
!  ############################ IMPORTANT NOTE, READ CAREFULLY #########

      IF ( FOURIER_COMPONENT .EQ. 0 ) THEN

!  With linearization.

        IF ( DO_ATMOS_LINEARIZATION ) THEN

          CALL VLIDORT_MISCSETUPS ( &
            DO_DELTAM_SCALING, NSTOKES, &
            NLAYERS, N_USER_LEVELS, &
            OMEGA_TOTAL_INPUT, &
            DELTAU_VERT_INPUT, GREEKMAT_TOTAL_INPUT, &
            NMOMENTS, NBEAMS, &
            PARTLAYERS_OUTFLAG, DO_PARTLAYERS, &
            PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES, &
            TAUGRID_INPUT, CHAPMAN_FACTORS, &
            MUELLER_INDEX, &
            DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, &
            DO_SPECIALIST_OPTION_3, NLAYERS_CUTOFF, &
            COS_SZANGLES, SUN_SZA_COSINES, &
            DO_SOLUTION_SAVING, NSTREAMS, &
            DO_TOA_CONTRIBS, QUAD_STREAMS, &
            N_USER_STREAMS, DO_USER_STREAMS, &
            USER_SECANTS, N_PARTLAYERS, &
            STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, &
            OMEGA_TOTAL, DELTAU_VERT, &
            PARTAU_VERT, GREEKMAT_TOTAL, &
            TAUGRID, DELTAU_SLANT, &
            TRUNC_FACTOR, FAC1, &
            OMEGA_GREEK, LAYER_PIS_CUTOFF, &
            SOLAR_BEAM_OPDEP, DO_REFLECTED_DIRECTBEAM, &
            INITIAL_TRANS, AVERAGE_SECANT, LOCAL_CSZA, &
            T_DELT_DISORDS, T_DISORDS_UTUP, &
            T_DISORDS_UTDN, T_DELT_MUBAR, &
            T_UTDN_MUBAR, T_UTUP_MUBAR, &
            T_DELT_USERM, T_UTDN_USERM, &
            T_UTUP_USERM, CUMTRANS, ITRANS_USERM )

          CALL VLIDORT_L_MISCSETUPS ( &
            DO_DELTAM_SCALING, NSTOKES, &
            NLAYERS, OMEGA_TOTAL_INPUT, &
            DELTAU_VERT_INPUT, GREEKMAT_TOTAL_INPUT, &
            NMOMENTS, NSTOKES_SQ, NBEAMS, &
            DO_ATMOS_LINEARIZATION, LAYER_VARY_FLAG, &
            LAYER_VARY_NUMBER, L_OMEGA_TOTAL_INPUT, &
            L_DELTAU_VERT_INPUT, L_GREEKMAT_TOTAL_INPUT, &
            DELTAU_SLANT, TRUNC_FACTOR, FAC1, &
            MUELLER_INDEX, OMEGA_GREEK, &
            DO_PLANE_PARALLEL, DO_SOLUTION_SAVING, &
            NSTREAMS, LAYER_PIS_CUTOFF, QUAD_STREAMS, &
            N_USER_STREAMS, DO_USER_STREAMS, &
            USER_SECANTS, N_PARTLAYERS, &
            PARTLAYERS_LAYERIDX, &
            DO_PROFILE_LINEARIZATION, DO_COLUMN_LINEARIZATION, &
            DELTAU_VERT, PARTAU_VERT, T_DELT_DISORDS, &
            T_DISORDS_UTUP, T_DISORDS_UTDN, &
            T_DELT_MUBAR, T_UTDN_MUBAR, &
            T_DELT_USERM, T_UTDN_USERM, &
            T_UTUP_USERM, AVERAGE_SECANT, &
            L_OMEGA_TOTAL, L_DELTAU_VERT, &
            L_GREEKMAT_TOTAL, L_DELTAU_SLANT, &
            DO_SCATMAT_VARIATION, L_TRUNC_FACTOR, &
            L_OMEGA_GREEK, L_AVERAGE_SECANT, &
            L_INITIAL_TRANS, L_T_DELT_DISORDS, &
            L_T_DISORDS_UTDN, L_T_DISORDS_UTUP, &
            L_T_DELT_MUBAR, L_T_UTDN_MUBAR, &
            L_T_DELT_USERM, L_T_UTDN_USERM, &
            L_T_UTUP_USERM )

          IF ( DO_INCLUDE_THERMEMISS ) THEN
            CALL THERMAL_SETUP_PLUS ( &
              DO_UPWELLING, DO_DNWELLING, DO_MSMODE_THERMAL, &
              NLAYERS, N_THERMAL_COEFFS, &
              THERMAL_BB_INPUT, DO_THERMAL_TRANSONLY, &
              N_USER_STREAMS, DO_USER_STREAMS, &
              USER_STREAMS, DO_PARTLAYERS, &
              N_PARTLAYERS, PARTLAYERS_LAYERIDX, &
              STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, &
              OMEGA_TOTAL, DELTAU_VERT, &
              PARTAU_VERT, T_DELT_USERM, &
              T_UTDN_USERM, T_UTUP_USERM, &
              DO_ATMOS_LINEARIZATION, LAYER_VARY_FLAG, &
              LAYER_VARY_NUMBER, &
              L_OMEGA_TOTAL, L_DELTAU_VERT, &
              L_T_DELT_USERM, L_T_UTDN_USERM, &
              L_T_UTUP_USERM, &
              THERMCOEFFS, DELTAU_POWER, XTAU_POWER, &
              TCOM1, T_DIRECT_UP, &
              T_DIRECT_DN, T_UT_DIRECT_UP, &
              T_UT_DIRECT_DN, &
              L_THERMCOEFFS, L_DELTAU_POWER, &
              L_XTAU_POWER, L_TCOM1, &
              L_T_DIRECT_UP, L_T_DIRECT_DN, &
              L_T_UT_DIRECT_UP, L_T_UT_DIRECT_DN )
          ENDIF

          IF ( DO_SOLAR_SOURCES ) THEN
            IF (.NOT.DO_SSFULL .OR. (DO_SSFULL.AND.DO_SSCORR_NADIR)) THEN
              CALL EMULT_MASTER ( &
                DO_UPWELLING, DO_DNWELLING, &
                NLAYERS, &
                LAYER_PIS_CUTOFF, NBEAMS, &
                N_USER_STREAMS, &
                USER_SECANTS, N_PARTLAYERS, &
                PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, &
                STERM_LAYERMASK_DN, &
                DELTAU_VERT, PARTAU_VERT, &
                T_DELT_MUBAR, T_UTDN_MUBAR, &
                T_DELT_USERM, T_UTDN_USERM, &
                T_UTUP_USERM, ITRANS_USERM, &
                AVERAGE_SECANT, &
                EMULT_HOPRULE, &
                SIGMA_M, SIGMA_P, &
                EMULT_UP, EMULT_DN, &
                UT_EMULT_UP, UT_EMULT_DN )

              CALL L_EMULT_MASTER ( &
                DO_UPWELLING, DO_DNWELLING, &
                NLAYERS, N_USER_LEVELS, &
                PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
                PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, &
                STERM_LAYERMASK_DN, &
                EMULT_UP, EMULT_DN, &
                DO_PROFILE_LINEARIZATION, DO_COLUMN_LINEARIZATION, &
                N_TOTALCOLUMN_WFS, LAYER_VARY_FLAG, &
                LAYER_VARY_NUMBER, DO_PLANE_PARALLEL, &
                LAYER_PIS_CUTOFF, NBEAMS, &
                N_USER_STREAMS, &
                T_DELT_MUBAR, T_DELT_USERM, &
                ITRANS_USERM, SIGMA_P, &
                L_AVERAGE_SECANT, L_INITIAL_TRANS, &
                L_T_DELT_MUBAR, L_T_DELT_USERM, &
                USER_SECANTS, &
                DELTAU_VERT, INITIAL_TRANS, &
                EMULT_HOPRULE, SIGMA_M, &
                L_DELTAU_VERT, T_UTUP_USERM, &
                UT_EMULT_UP, &
                L_T_UTDN_MUBAR, L_T_UTUP_USERM, &
                PARTAU_VERT, UT_EMULT_DN, &
                L_T_UTDN_USERM, &
                L_EMULT_UP, L_EMULT_DN, &
                L_UT_EMULT_UP, L_UT_EMULT_DN )
            ENDIF
          ENDIF

!  No linearization

        ELSE

          CALL VLIDORT_MISCSETUPS ( &
            DO_DELTAM_SCALING, NSTOKES, &
            NLAYERS, N_USER_LEVELS, &
            OMEGA_TOTAL_INPUT, &
            DELTAU_VERT_INPUT, GREEKMAT_TOTAL_INPUT, &
            NMOMENTS, NBEAMS, &
            PARTLAYERS_OUTFLAG, DO_PARTLAYERS, &
            PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES, &
            TAUGRID_INPUT, CHAPMAN_FACTORS, &
            MUELLER_INDEX, &
            DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, &
            DO_SPECIALIST_OPTION_3, NLAYERS_CUTOFF, &
            COS_SZANGLES, SUN_SZA_COSINES, &
            DO_SOLUTION_SAVING, NSTREAMS, &
            DO_TOA_CONTRIBS, QUAD_STREAMS, &
            N_USER_STREAMS, DO_USER_STREAMS, &
            USER_SECANTS, N_PARTLAYERS, &
            STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, &
            OMEGA_TOTAL, DELTAU_VERT, &
            PARTAU_VERT, GREEKMAT_TOTAL, &
            TAUGRID, DELTAU_SLANT, &
            TRUNC_FACTOR, FAC1, &
            OMEGA_GREEK, LAYER_PIS_CUTOFF, &
            SOLAR_BEAM_OPDEP, DO_REFLECTED_DIRECTBEAM, &
            INITIAL_TRANS, AVERAGE_SECANT, LOCAL_CSZA, &
            T_DELT_DISORDS, T_DISORDS_UTUP, &
            T_DISORDS_UTDN, T_DELT_MUBAR, &
            T_UTDN_MUBAR, T_UTUP_MUBAR, &
            T_DELT_USERM, T_UTDN_USERM, &
            T_UTUP_USERM, CUMTRANS, ITRANS_USERM )

          IF ( DO_INCLUDE_THERMEMISS ) THEN
            CALL THERMAL_SETUP ( &
              DO_UPWELLING, DO_DNWELLING, DO_MSMODE_THERMAL, &
              NLAYERS, N_THERMAL_COEFFS, &
              THERMAL_BB_INPUT, DO_THERMAL_TRANSONLY, &
              N_USER_STREAMS, DO_USER_STREAMS, &
              USER_STREAMS, DO_PARTLAYERS, &
              N_PARTLAYERS, PARTLAYERS_LAYERIDX, &
              STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, &
              OMEGA_TOTAL, DELTAU_VERT, &
              PARTAU_VERT, T_DELT_USERM, &
              T_UTDN_USERM, T_UTUP_USERM, &
              THERMCOEFFS, DELTAU_POWER, &
              XTAU_POWER, &
              TCOM1, T_DIRECT_UP, &
              T_DIRECT_DN, T_UT_DIRECT_UP, &
              T_UT_DIRECT_DN )
          ENDIF

          IF ( DO_SOLAR_SOURCES ) THEN
            IF (.NOT.DO_SSFULL .OR. (DO_SSFULL.AND.DO_SSCORR_NADIR)) THEN
              CALL EMULT_MASTER ( &
                DO_UPWELLING, DO_DNWELLING, &
                NLAYERS, &
                LAYER_PIS_CUTOFF, NBEAMS, &
                N_USER_STREAMS, &
                USER_SECANTS, N_PARTLAYERS, &
                PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, &
                STERM_LAYERMASK_DN, &
                DELTAU_VERT, PARTAU_VERT, &
                T_DELT_MUBAR, T_UTDN_MUBAR, &
                T_DELT_USERM, T_UTDN_USERM, &
                T_UTUP_USERM, ITRANS_USERM, &
                AVERAGE_SECANT, &
                EMULT_HOPRULE, &
                SIGMA_M, SIGMA_P, &
                EMULT_UP, EMULT_DN, &
                UT_EMULT_UP, UT_EMULT_DN )
            ENDIF
          ENDIF

        ENDIF

!  End Fourier 0 clause

      ENDIF

!  #####################
!  Correction operations
!  #####################

!  Not required if no solar sources

      IF ( DO_SOLAR_SOURCES ) THEN

!  Single scatter correction (pre-calculation)
!  -------------------------------------------

!   Must be done after MISCSETUPS, as we need the multipliers and
!    transmittance factors for the SUN and LOS paths.
!   Standard Code added 6 May 2005. Replaces call in Master routine.
!   Linearized code 13 July 2005.  Replaces former code in L_Master module

!      Version 2.2. Added call to the new outgoing sphericity correction
!      Version 2.4. Added call bulk/column linearizations

       IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
        IF ( DO_USER_STREAMS ) THEN

!  Regular (nadir view) SS correction

          IF ( DO_SSCORR_NADIR ) THEN
            CALL VLIDORT_SSCORR_NADIR ( &
              SS_FLUX_MULTIPLIER, &
              DO_SSCORR_TRUNCATION, DO_REFRACTIVE_GEOMETRY, &
              DO_DELTAM_SCALING, DO_UPWELLING, &
              DO_DNWELLING, NSTOKES, &
              NLAYERS, NGREEK_MOMENTS_INPUT, &
              N_USER_RELAZMS, USER_RELAZMS, &
              N_USER_LEVELS, &
              OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT, &
              DO_TOA_CONTRIBS, &
              FLUXVEC, COS_SZANGLES, &
              SIN_SZANGLES, SUN_SZA_COSINES, &
              NBEAMS, N_USER_STREAMS, &
              LAYER_MAXMOMENTS, USER_STREAMS, &
              PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
              UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, &
              PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, &
              STERM_LAYERMASK_DN, N_GEOMETRIES, &
              VZA_OFFSETS, GREEKMAT_TOTAL, &
              TRUNC_FACTOR, T_DELT_USERM, &
              T_UTDN_USERM, T_UTUP_USERM, &
              CUMTRANS, &
              EMULT_UP, EMULT_DN, &
              UT_EMULT_UP, UT_EMULT_DN, &
              ZMAT_UP, ZMAT_DN, TMS, &
              SSFDEL, SS_CUMSOURCE_UP, SS_CUMSOURCE_DN, &
              STOKES_SS, SS_CONTRIBS )

            IF ( DO_PROFILE_LINEARIZATION ) THEN
              DO LAYER_TO_VARY = 1, NLAYERS

                IF ( LAYER_VARY_FLAG(LAYER_TO_VARY) ) THEN
                  N_PARAMETERS = LAYER_VARY_NUMBER(LAYER_TO_VARY)
                  CALL VLIDORT_LP_SSCORR_NADIR ( &
                    SS_FLUX_MULTIPLIER, LAYER_TO_VARY, N_PARAMETERS, &
                    DO_SSCORR_TRUNCATION, DO_REFRACTIVE_GEOMETRY, &
                    DO_DELTAM_SCALING, DO_UPWELLING, &
                    DO_DNWELLING, NSTOKES, &
                    NLAYERS, NGREEK_MOMENTS_INPUT, &
                    N_USER_RELAZMS, USER_RELAZMS, &
                    N_USER_LEVELS, GREEKMAT_TOTAL_INPUT, &
                    FLUXVEC, COS_SZANGLES, &
                    SIN_SZANGLES, SUN_SZA_COSINES, &
                    NMOMENTS, NBEAMS, &
                    N_USER_STREAMS, LAYER_MAXMOMENTS, &
                    USER_STREAMS, PARTLAYERS_OUTFLAG, &
                    PARTLAYERS_OUTINDEX, UTAU_LEVEL_MASK_UP, &
                    UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX, &
                    STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, &
                    N_GEOMETRIES, VZA_OFFSETS, &
                    L_OMEGA_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT, &
                    TRUNC_FACTOR, T_DELT_USERM, &
                    T_UTDN_USERM, T_UTUP_USERM, &
                    EMULT_UP, EMULT_DN, &
                    UT_EMULT_UP, UT_EMULT_DN, &
                    ZMAT_UP, ZMAT_DN, &
                    TMS, SSFDEL, &
                    SS_CUMSOURCE_UP, SS_CUMSOURCE_DN, &
                    DO_SCATMAT_VARIATION, L_T_DELT_USERM, &
                    L_T_UTDN_USERM, L_T_UTUP_USERM, &
                    L_EMULT_UP, L_EMULT_DN, &
                    L_UT_EMULT_UP, L_UT_EMULT_DN, &
                    PROFILEWF_SS )

                ENDIF
              ENDDO
            ENDIF

            IF ( DO_COLUMN_LINEARIZATION ) THEN
              CALL VLIDORT_LC_SSCORR_NADIR ( &
                SS_FLUX_MULTIPLIER, N_TOTALCOLUMN_WFS, &
                DO_SSCORR_TRUNCATION, DO_REFRACTIVE_GEOMETRY, &
                DO_DELTAM_SCALING, DO_UPWELLING, &
                DO_DNWELLING, NSTOKES, &
                NLAYERS, NGREEK_MOMENTS_INPUT, &
                N_USER_RELAZMS, USER_RELAZMS, &
                N_USER_LEVELS, GREEKMAT_TOTAL_INPUT, &
                FLUXVEC, COS_SZANGLES, &
                SIN_SZANGLES, SUN_SZA_COSINES, &
                NMOMENTS, NBEAMS, &
                N_USER_STREAMS, LAYER_MAXMOMENTS, &
                USER_STREAMS, PARTLAYERS_OUTFLAG, &
                PARTLAYERS_OUTINDEX, UTAU_LEVEL_MASK_UP, &
                UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX, &
                STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, &
                N_GEOMETRIES, VZA_OFFSETS, &
                L_OMEGA_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT, &
                TRUNC_FACTOR, T_DELT_USERM, &
                T_UTDN_USERM, T_UTUP_USERM, &
                EMULT_UP, EMULT_DN, &
                UT_EMULT_UP, UT_EMULT_DN, &
                ZMAT_UP, ZMAT_DN, &
                TMS, SSFDEL, &
                SS_CUMSOURCE_UP, SS_CUMSOURCE_DN, &
                DO_SCATMAT_VARIATION, L_T_DELT_USERM, &
                L_T_UTDN_USERM, L_T_UTUP_USERM, &
                L_EMULT_UP, L_EMULT_DN, &
                L_UT_EMULT_UP, L_UT_EMULT_DN, &
                COLUMNWF_SS )
            ENDIF
          ENDIF

!  Outgoing sphericity correction ------ Added 20 March 2007.
!     Calling not changed for column Jacobians (Version 2.4, December 20

          IF ( DO_SSCORR_OUTGOING ) THEN
           IF ( DO_ATMOS_LINEARIZATION ) THEN

            CALL VLIDORT_L_SSCORR_OUTGOING ( &
              SS_FLUX_MULTIPLIER, DO_SSCORR_TRUNCATION, DO_DELTAM_SCALING, &
              DO_UPWELLING, DO_DNWELLING, NSTOKES, NLAYERS, &
              NGREEK_MOMENTS_INPUT, N_USER_RELAZMS, N_USER_LEVELS, &
              EARTH_RADIUS, HEIGHT_GRID, &
              OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT, &
              DO_TOA_CONTRIBS, FLUXVEC, NMOMENTS, &
              NBEAMS, N_USER_STREAMS, &
              LAYER_MAXMOMENTS, USER_VZANGLES_ADJUST, &
              SZANGLES_ADJUST, USER_RELAZMS_ADJUST, &
              N_PARTLAYERS, NFINELAYERS, &
              PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
              UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, &
              PARTLAYERS_LAYERIDX, &
              STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, &
              N_GEOMETRIES, VZA_OFFSETS, &
              DELTAU_VERT, PARTAU_VERT, TRUNC_FACTOR, &
              DO_PROFILE_LINEARIZATION, DO_COLUMN_LINEARIZATION, &
              DO_ATMOS_LINEARIZATION, N_TOTALCOLUMN_WFS, &
              LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
              L_OMEGA_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT, &
              L_DELTAU_VERT, DO_SCATMAT_VARIATION, &
              UP_MULTIPLIERS, DN_MULTIPLIERS, UP_LOSTRANS, DN_LOSTRANS, &
              UP_MULTIPLIERS_UT, DN_MULTIPLIERS_UT, &
              UP_LOSTRANS_UT, DN_LOSTRANS_UT, &
              ZMAT_UP, ZMAT_DN, TMS, SSFDEL, &
              SS_CUMSOURCE_UP, SS_CUMSOURCE_DN, &
              BOA_ATTN, STOKES_SS, SS_CONTRIBS, &
              L_UP_MULTIPLIERS, L_DN_MULTIPLIERS, &
              L_UP_LOSTRANS, L_DN_LOSTRANS, &
              L_UP_MULTIPLIERS_UT, L_DN_MULTIPLIERS_UT, &
              L_UP_LOSTRANS_UT, L_DN_LOSTRANS_UT, &
              L_BOA_ATTN, &
              PROFILEWF_SS, COLUMNWF_SS, &
              FAIL, MESSAGE, TRACE_1 )

            IF ( FAIL ) THEN
              TRACE_2 = 'Linearized SS correction outgoing '// &
                           'failed in VLIDORT_L_FOURIER'
              STATUS  = VLIDORT_SERIOUS
              RETURN
            ENDIF

           ELSE

            CALL VLIDORT_SSCORR_OUTGOING ( &
              SS_FLUX_MULTIPLIER, &
              DO_SSCORR_TRUNCATION, DO_DELTAM_SCALING, &
              DO_UPWELLING, DO_DNWELLING, &
              NSTOKES, NLAYERS, &
              NGREEK_MOMENTS_INPUT, &
              N_USER_RELAZMS, N_USER_LEVELS, &
              EARTH_RADIUS, HEIGHT_GRID, &
              OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT, &
              DO_TOA_CONTRIBS, &
              FLUXVEC, NBEAMS, &
              N_USER_STREAMS, LAYER_MAXMOMENTS, &
              USER_VZANGLES_ADJUST, SZANGLES_ADJUST, &
              USER_RELAZMS_ADJUST, &
              UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, &
              STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, &
              N_GEOMETRIES, VZA_OFFSETS, &
              DELTAU_VERT, TRUNC_FACTOR, &
              N_PARTLAYERS, NFINELAYERS, &
              PARTLAYERS_LAYERIDX, PARTAU_VERT, &
              PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
              UP_MULTIPLIERS, DN_MULTIPLIERS, &
              UP_MULTIPLIERS_UT, DN_MULTIPLIERS_UT, &
              UP_LOSTRANS, DN_LOSTRANS, &
              UP_LOSTRANS_UT, DN_LOSTRANS_UT, &
              ZMAT_UP, ZMAT_DN, TMS, &
              SSFDEL, SS_CUMSOURCE_UP, SS_CUMSOURCE_DN, &
              BOA_ATTN, STOKES_SS, SS_CONTRIBS, &
              FAIL, MESSAGE, TRACE_1)

            IF ( FAIL ) THEN
              TRACE_2 = 'SS correction outgoing '// &
                         'failed in VLIDORT_L_FOURIER'
              STATUS  = VLIDORT_SERIOUS
              RETURN
            ENDIF

           ENDIF
          ENDIF

!  Exact direct beam Calculation
!  ------------------------------

!   Must be done after MISCSETUPS, as we need the SUN transmittance fact
!   Code added 6 May 2005. Standard case.
!   Profile/surface Linearization code added July 13 2005.
!     Column        Linearization code added December 2008

!  Reflected Direct beam attenuation + DB correction (if flagged)
!  Surface leaving option and variables added 17 May 2012
!  Full single scatter calculation, return after completion

          IF ( DO_SSFULL ) THEN
            CALL VLIDORT_DBCORRECTION ( &
              SS_FLUX_MULTIPLIER, &
              DO_SSCORR_OUTGOING, DO_REFRACTIVE_GEOMETRY, &
              DO_UPWELLING, NSTOKES, &
              NLAYERS, &
              SZA_LOCAL_INPUT, N_USER_RELAZMS, &
              N_USER_LEVELS, &
              DO_LAMBERTIAN_SURFACE, &
              LAMBERTIAN_ALBEDO, EXACTDB_BRDFUNC, &
              DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, &
              SLTERM_ISOTROPIC, SLTERM_USERANGLES, &
              FLUXVEC, COS_SZANGLES, &
              NBEAMS, N_USER_STREAMS, &
              MUELLER_INDEX, SZANGLES_ADJUST, &
              PARTLAYERS_OUTFLAG, &
              PARTLAYERS_OUTINDEX, UTAU_LEVEL_MASK_UP, &
              PARTLAYERS_LAYERIDX, N_GEOMETRIES, &
              VZA_OFFSETS, &
              SOLAR_BEAM_OPDEP, DO_REFLECTED_DIRECTBEAM, &
              T_DELT_USERM, T_UTUP_USERM, &
              UP_LOSTRANS, UP_LOSTRANS_UT, BOA_ATTN, &
              DB_CUMSOURCE, EXACTDB_SOURCE, &
              ATTN_DB_SAVE, STOKES_DB )

            IF ( DO_PROFILE_LINEARIZATION ) THEN
              CALL VLIDORT_LAP_DBCORRECTION ( &
                SS_FLUX_MULTIPLIER, &
                DO_SSCORR_OUTGOING, DO_REFRACTIVE_GEOMETRY, &
                DO_UPWELLING, NSTOKES, NLAYERS, &
                SZA_LOCAL_INPUT, N_USER_RELAZMS, &
                N_USER_LEVELS, &
                DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
                COS_SZANGLES, NBEAMS, &
                N_USER_STREAMS, SZANGLES_ADJUST, &
                PARTLAYERS_OUTFLAG, &
                PARTLAYERS_OUTINDEX, UTAU_LEVEL_MASK_UP, &
                PARTLAYERS_LAYERIDX, N_GEOMETRIES, &
                VZA_OFFSETS, SOLAR_BEAM_OPDEP, &
                DO_REFLECTED_DIRECTBEAM, T_DELT_USERM, &
                T_UTUP_USERM, &
                LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
                L_DELTAU_SLANT, L_T_DELT_USERM, &
                L_T_UTUP_USERM, &
                DB_CUMSOURCE, EXACTDB_SOURCE, &
                UP_LOSTRANS, UP_LOSTRANS_UT, &
                BOA_ATTN, ATTN_DB_SAVE, &
                L_UP_LOSTRANS, L_UP_LOSTRANS_UT, &
                L_BOA_ATTN, L_EXACTDB_SOURCE, &
                L_DB_CUMSOURCE, PROFILEWF_DB )
            ENDIF

            IF ( DO_COLUMN_LINEARIZATION ) THEN
              CALL VLIDORT_LAC_DBCORRECTION ( &
                SS_FLUX_MULTIPLIER, &
                DO_SSCORR_OUTGOING, DO_REFRACTIVE_GEOMETRY, &
                DO_UPWELLING, NSTOKES, NLAYERS, &
                SZA_LOCAL_INPUT, N_USER_RELAZMS, &
                N_USER_LEVELS, &
                DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
                COS_SZANGLES, NBEAMS, &
                N_USER_STREAMS, SZANGLES_ADJUST, &
                PARTLAYERS_OUTFLAG, &
                PARTLAYERS_OUTINDEX, UTAU_LEVEL_MASK_UP, &
                PARTLAYERS_LAYERIDX, N_GEOMETRIES, &
                VZA_OFFSETS, DELTAU_SLANT, &
                SOLAR_BEAM_OPDEP, DO_REFLECTED_DIRECTBEAM, &
                T_DELT_USERM, T_UTUP_USERM, &
                N_TOTALCOLUMN_WFS, &
                L_DELTAU_VERT, L_T_DELT_USERM, &
                L_T_UTUP_USERM, &
                DB_CUMSOURCE, EXACTDB_SOURCE, &
                UP_LOSTRANS, UP_LOSTRANS_UT, &
                BOA_ATTN, ATTN_DB_SAVE, &
                L_UP_LOSTRANS, L_UP_LOSTRANS_UT, &
                L_BOA_ATTN, L_EXACTDBC_SOURCE, &
                L_DB_CUMSOURCE, COLUMNWF_DB )
            ENDIF

            IF ( DO_SURFACE_LINEARIZATION ) THEN
              CALL VLIDORT_LS_DBCORRECTION ( &
                SS_FLUX_MULTIPLIER, &
                DO_SSCORR_OUTGOING, DO_UPWELLING, &
                NSTOKES, NLAYERS, &
                N_USER_RELAZMS, N_USER_LEVELS, &
                DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
                FLUXVEC, NBEAMS, &
                N_USER_STREAMS, MUELLER_INDEX, &
                PARTLAYERS_OUTFLAG, &
                PARTLAYERS_OUTINDEX, UTAU_LEVEL_MASK_UP, &
                PARTLAYERS_LAYERIDX, N_GEOMETRIES, &
                VZA_OFFSETS, &
                DO_REFLECTED_DIRECTBEAM, T_DELT_USERM, &
                T_UTUP_USERM, ATTN_DB_SAVE, &
                N_SURFACE_WFS, LS_EXACTDB_BRDFUNC, &
                UP_LOSTRANS, UP_LOSTRANS_UT, &
                LS_EXACTDB_SOURCE, LS_DB_CUMSOURCE, &
                SURFACEWF_DB )
            ENDIF

            RETURN
          ELSE
            IF ( DO_DBCORRECTION ) THEN
              CALL VLIDORT_DBCORRECTION ( &
                SS_FLUX_MULTIPLIER, &
                DO_SSCORR_OUTGOING, DO_REFRACTIVE_GEOMETRY, &
                DO_UPWELLING, NSTOKES, &
                NLAYERS, &
                SZA_LOCAL_INPUT, N_USER_RELAZMS, &
                N_USER_LEVELS, &
                DO_LAMBERTIAN_SURFACE, &
                LAMBERTIAN_ALBEDO, EXACTDB_BRDFUNC, &
                DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, &
                SLTERM_ISOTROPIC, SLTERM_USERANGLES, &
                FLUXVEC, COS_SZANGLES, &
                NBEAMS, N_USER_STREAMS, &
                MUELLER_INDEX, SZANGLES_ADJUST, &
                PARTLAYERS_OUTFLAG, &
                PARTLAYERS_OUTINDEX, UTAU_LEVEL_MASK_UP, &
                PARTLAYERS_LAYERIDX, N_GEOMETRIES, &
                VZA_OFFSETS, &
                SOLAR_BEAM_OPDEP, DO_REFLECTED_DIRECTBEAM, &
                T_DELT_USERM, T_UTUP_USERM, &
                UP_LOSTRANS, UP_LOSTRANS_UT, BOA_ATTN, &
                DB_CUMSOURCE, EXACTDB_SOURCE, &
                ATTN_DB_SAVE, STOKES_DB )

              IF ( DO_PROFILE_LINEARIZATION ) THEN
                CALL VLIDORT_LAP_DBCORRECTION ( &
                  SS_FLUX_MULTIPLIER, &
                  DO_SSCORR_OUTGOING, DO_REFRACTIVE_GEOMETRY, &
                  DO_UPWELLING, NSTOKES, NLAYERS, &
                  SZA_LOCAL_INPUT, N_USER_RELAZMS, &
                  N_USER_LEVELS, &
                  DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
                  COS_SZANGLES, NBEAMS, &
                  N_USER_STREAMS, SZANGLES_ADJUST, &
                  PARTLAYERS_OUTFLAG, &
                  PARTLAYERS_OUTINDEX, UTAU_LEVEL_MASK_UP, &
                  PARTLAYERS_LAYERIDX, N_GEOMETRIES, &
                  VZA_OFFSETS, SOLAR_BEAM_OPDEP, &
                  DO_REFLECTED_DIRECTBEAM, T_DELT_USERM, &
                  T_UTUP_USERM, &
                  LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
                  L_DELTAU_SLANT, L_T_DELT_USERM, &
                  L_T_UTUP_USERM, &
                  DB_CUMSOURCE, EXACTDB_SOURCE, &
                  UP_LOSTRANS, UP_LOSTRANS_UT, &
                  BOA_ATTN, ATTN_DB_SAVE, &
                  L_UP_LOSTRANS, L_UP_LOSTRANS_UT, &
                  L_BOA_ATTN, L_EXACTDB_SOURCE, &
                  L_DB_CUMSOURCE, PROFILEWF_DB )
              ENDIF

              IF ( DO_COLUMN_LINEARIZATION ) THEN
                CALL VLIDORT_LAC_DBCORRECTION ( &
                  SS_FLUX_MULTIPLIER, &
                  DO_SSCORR_OUTGOING, DO_REFRACTIVE_GEOMETRY, &
                  DO_UPWELLING, NSTOKES, NLAYERS, &
                  SZA_LOCAL_INPUT, N_USER_RELAZMS, &
                  N_USER_LEVELS, &
                  DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
                  COS_SZANGLES, NBEAMS, &
                  N_USER_STREAMS, SZANGLES_ADJUST, &
                  PARTLAYERS_OUTFLAG, &
                  PARTLAYERS_OUTINDEX, UTAU_LEVEL_MASK_UP, &
                  PARTLAYERS_LAYERIDX, N_GEOMETRIES, &
                  VZA_OFFSETS, DELTAU_SLANT, &
                  SOLAR_BEAM_OPDEP, DO_REFLECTED_DIRECTBEAM, &
                  T_DELT_USERM, T_UTUP_USERM, &
                  N_TOTALCOLUMN_WFS, &
                  L_DELTAU_VERT, L_T_DELT_USERM, &
                  L_T_UTUP_USERM, &
                  DB_CUMSOURCE, EXACTDB_SOURCE, &
                  UP_LOSTRANS, UP_LOSTRANS_UT, &
                  BOA_ATTN, ATTN_DB_SAVE, &
                  L_UP_LOSTRANS, L_UP_LOSTRANS_UT, &
                  L_BOA_ATTN, L_EXACTDBC_SOURCE, &
                  L_DB_CUMSOURCE, COLUMNWF_DB )
              ENDIF

              IF ( DO_SURFACE_LINEARIZATION ) THEN
                CALL VLIDORT_LS_DBCORRECTION ( &
                  SS_FLUX_MULTIPLIER, &
                  DO_SSCORR_OUTGOING, DO_UPWELLING, &
                  NSTOKES, NLAYERS, &
                  N_USER_RELAZMS, N_USER_LEVELS, &
                  DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
                  FLUXVEC, NBEAMS, &
                  N_USER_STREAMS, MUELLER_INDEX, &
                  PARTLAYERS_OUTFLAG, &
                  PARTLAYERS_OUTINDEX, UTAU_LEVEL_MASK_UP, &
                  PARTLAYERS_LAYERIDX, N_GEOMETRIES, &
                  VZA_OFFSETS, &
                  DO_REFLECTED_DIRECTBEAM, T_DELT_USERM, &
                  T_UTUP_USERM, ATTN_DB_SAVE, &
                  N_SURFACE_WFS, LS_EXACTDB_BRDFUNC, &
                  UP_LOSTRANS, UP_LOSTRANS_UT, &
                  LS_EXACTDB_SOURCE, LS_DB_CUMSOURCE, &
                  SURFACEWF_DB )
              ENDIF

            ENDIF
          ENDIF

!  End fourier = 0, user stream and solar source  clauses

        ENDIF
       ENDIF

!  Surface direct beam
!    New Surface-Leaving arguments added, 17 May 2012

       IF ( .not. DO_SSFULL ) THEN
         CALL VLIDORT_SURFACE_DIRECTBEAM ( &
           DELTA_FACTOR, FOURIER_COMPONENT, &
           DO_REFRACTIVE_GEOMETRY, NSTOKES, &
           NSTREAMS, NLAYERS, &
           FLUX_FACTOR, &
           SZA_LOCAL_INPUT, DO_LAMBERTIAN_SURFACE, &
           LAMBERTIAN_ALBEDO, BRDF_F, &
           BRDF_F_0, USER_BRDF_F, &
           USER_BRDF_F_0, &
           DO_SURFACE_LEAVING, DO_SL_ISOTROPIC,  &
           SLTERM_ISOTROPIC, SLTERM_F_0, USER_SLTERM_F_0, &
           FLUXVEC, COS_SZANGLES, &
           NBEAMS, N_USER_STREAMS, &
           MUELLER_INDEX, DO_USER_STREAMS, &
           SOLAR_BEAM_OPDEP, DO_REFLECTED_DIRECTBEAM, &
           ATMOS_ATTN, DIRECT_BEAM, &
           USER_DIRECT_BEAM )
       ENDIF

!  End solar sources only clause for corrections

      ENDIF

!  ###############################################################
!     END OF SSCORRECTION STUFF
!  ###############################################################

!  Get Pi Matrices for this Fourier component
!    Including all beam angles !!

      IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
        CALL VLIDORT_PIMATRIX_SETUP ( &
          FOURIER_COMPONENT, &
          DO_REFRACTIVE_GEOMETRY, NSTOKES, &
          NSTREAMS, NLAYERS, &
          COS_SZANGLES, SUN_SZA_COSINES, &
          QUAD_STREAMS, NMOMENTS, &
          NBEAMS, N_USER_STREAMS, &
          DO_USER_STREAMS, USER_STREAMS, &
          MUELLER_INDEX, DMAT, &
          PI_XQP, PI_XQM, &
          PI_XUP, PI_XUM, &
          PI_X0P, &
          PI_XQM_POST, PI_XQM_PRE, &
          PI_XQP_PRE, PI_XUM_POST, &
          PI_XUP_PRE )
      ENDIF

!  #########################################################
!  RT differential equation Eigensolutions + linearizations
!  #########################################################

!  Go to continuation point for thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) GO TO 8899

!  Start layer loop

      DO LAYER = 1, NLAYERS

!  Get Discrete ordinate solutions for this layer

        CALL VLIDORT_QHOM_SOLUTION ( &
          LAYER, FOURIER_COMPONENT, &
          DO_SOLUTION_SAVING, NSTOKES, &
          NSTREAMS, N_USER_LEVELS, &
          QUAD_STREAMS, QUAD_HALFWTS, &
          NMOMENTS, NSTKS_NSTRMS, &
          DO_REAL_EIGENSOLVER, DO_LAYER_SCATTERING, &
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
          PARTLAYERS_LAYERIDX, &
          DELTAU_VERT, PARTAU_VERT, &
          OMEGA_GREEK, T_DELT_DISORDS, &
          T_DISORDS_UTUP, T_DISORDS_UTDN, &
          PI_XQP, PI_XQM, &
          PI_XQM_PRE, &
          SAB, DAB, &
          EIGENMAT_SAVE, REAL_KSQ, &
          IMAG_KSQ, LEFT_EVEC, &
          RITE_EVEC, EIGENDEGEN, &
          EIGENMASK_R, EIGENMASK_C, &
          FWD_SUMVEC, FWD_DIFVEC, &
          T_DELT_EIGEN, T_UTUP_EIGEN, &
          T_UTDN_EIGEN, &
          K_REAL, K_COMPLEX, &
          KEIGEN, KEIGEN_CSQ, &
          SOLA_XPOS, SOLB_XNEG, &
          STATUS_SUB, MESSAGE, TRACE_1 )

!  .. error tracing

        IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
          write(CF,'(I2)')FOURIER_COMPONENT
          TRACE_2 = &
          'Called in VLIDORT_L_FOURIER, Fourier component '//CF
          STATUS = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  Get Post-processing ("user") solutions for this layer

        IF  ( STERM_LAYERMASK_UP(LAYER) .OR. &
              STERM_LAYERMASK_DN(LAYER) ) THEN
          IF ( DO_USER_STREAMS ) THEN
            CALL VLIDORT_UHOM_SOLUTION ( &
              LAYER, FOURIER_COMPONENT, &
              DO_UPWELLING, DO_DNWELLING, &
              DO_DEBUG_WRITE, DO_FDTEST, &
              NSTOKES, NSTREAMS, &
              QUAD_HALFWTS, NMOMENTS, &
              N_USER_STREAMS, DO_LAYER_SCATTERING, &
              LOCAL_UM_START, USER_STREAMS, &
              USER_SECANTS, &
              OMEGA_GREEK, PI_XQP, &
              PI_XUP, PI_XUM, PI_XQM_PRE, &
              PI_XUM_POST, PI_XUP_PRE, &
              K_REAL, K_COMPLEX, KEIGEN, &
              KEIGEN_CSQ, SOLA_XPOS, &
              SOLB_XNEG, &
              UHOM_DNDN, UHOM_DNUP, &
              UHOM_UPDN, UHOM_UPUP, &
              HSINGO, HELPSTOKES, &
              ZETA_M, ZETA_P )
          ENDIF
        ENDIF

!  Additional Solutions for Linearization
!  --------------------------------------

        IF ( DO_ATMOS_LINEARIZATION ) THEN

!  Parameter control

          IF ( DO_PROFILE_LINEARIZATION ) THEN
            DO_RTSOL_VARY = LAYER_VARY_FLAG(LAYER)
            NPARAMS_VARY  = LAYER_VARY_NUMBER(LAYER)
          ELSE IF ( DO_COLUMN_LINEARIZATION ) THEN
            DO_RTSOL_VARY = .TRUE.
            NPARAMS_VARY  = N_TOTALCOLUMN_WFS
          ENDIF

!  Linearizations of Discrete Ordinate solution

!         if ( do_write_screen) write(*,*)'l_homsolution',layer
          CALL VLIDORT_L_QHOM_SOLUTION ( &
            LAYER, FOURIER_COMPONENT, DO_RTSOL_VARY, NPARAMS_VARY, &
            DO_SOLUTION_SAVING, DO_DEBUG_WRITE, &
            NSTOKES, NSTREAMS, &
            N_USER_LEVELS, &
            QUAD_STREAMS, QUAD_HALFWTS, &
            NMOMENTS, NSTREAMS_2, &
            NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
            DO_REAL_EIGENSOLVER, DO_LAYER_SCATTERING, &
            PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
            PARTLAYERS_LAYERIDX, DMAT, &
            DELTAU_VERT, PARTAU_VERT, &
            T_DELT_EIGEN, T_UTUP_EIGEN, &
            T_UTDN_EIGEN, &
            PI_XQP, PI_XQM_PRE, &
            SAB, DAB, &
            EIGENMAT_SAVE, REAL_KSQ, &
            IMAG_KSQ, LEFT_EVEC, &
            RITE_EVEC, EIGENDEGEN, &
            EIGENMASK_R, EIGENMASK_C, &
            K_REAL, K_COMPLEX, &
            KEIGEN, KEIGEN_CSQ, &
            FWD_SUMVEC, FWD_DIFVEC, &
            L_DELTAU_VERT, L_OMEGA_GREEK, &
            L_T_DELT_DISORDS, &
            L_T_DISORDS_UTDN, L_T_DISORDS_UTUP, &
            L_SAB, L_DAB, &
            L_EIGENMAT, L_KEIGEN, &
            L_SOLA_XPOS, L_SOLB_XNEG, &
            L_T_DELT_EIGEN, L_T_UTUP_EIGEN, &
            L_T_UTDN_EIGEN, &
            STATUS_SUB, MESSAGE, TRACE_1 )

!  .. error tracing

          IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
            write(CF,'(I2)')FOURIER_COMPONENT
            TRACE_2 = &
             'Called in VLIDORT_L_FOURIER, Fourier component '//CF
            STATUS = VLIDORT_SERIOUS
            RETURN
          ENDIF

!  Get Linearizations of ("user") solutions for this layer

          CALL VLIDORT_L_UHOM_SOLUTION ( &
            LAYER, FOURIER_COMPONENT, DO_RTSOL_VARY, &
            NPARAMS_VARY, &
            DO_UPWELLING, DO_DNWELLING, &
            DO_DEBUG_WRITE, &
            NSTOKES, NSTREAMS, &
            QUAD_HALFWTS, NMOMENTS, &
            N_USER_STREAMS, DO_LAYER_SCATTERING, &
            LOCAL_UM_START, &
            USER_SECANTS, STERM_LAYERMASK_UP, &
            STERM_LAYERMASK_DN, &
            OMEGA_GREEK, PI_XQP, PI_XUP, &
            PI_XUM, PI_XQM_PRE, &
            PI_XUM_POST, PI_XUP_PRE, &
            K_REAL, K_COMPLEX, &
            KEIGEN, HELPSTOKES, &
            UHOM_DNDN, UHOM_DNUP, &
            UHOM_UPDN, UHOM_UPUP, &
            ZETA_M, ZETA_P, L_KEIGEN, &
            L_SOLA_XPOS, L_SOLB_XNEG, &
            L_OMEGA_GREEK, &
            L_UHOM_DNDN, L_UHOM_DNUP, &
            L_UHOM_UPDN, L_UHOM_UPUP, &
            L_ZETA_M, L_ZETA_P )

!  End linearization control

        ENDIF

!  end layer loop

      ENDDO

!  Prepare homogeneous solution multipliers

      CALL HMULT_MASTER ( &
        DO_UPWELLING, DO_DNWELLING, &
        N_USER_LEVELS, &
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
        PARTLAYERS_LAYERIDX, NLAYERS, &
        N_USER_STREAMS, LOCAL_UM_START, &
        USER_SECANTS, &
        STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, &
        DELTAU_VERT, T_DELT_EIGEN, &
        T_DELT_USERM, &
        K_REAL, K_COMPLEX, &
        HSINGO, ZETA_M, ZETA_P, &
        PARTAU_VERT, &
        T_UTUP_EIGEN, T_UTDN_EIGEN, &
        T_UTUP_USERM, T_UTDN_USERM, &
        HMULT_1, HMULT_2, &
        UT_HMULT_UU, UT_HMULT_UD, &
        UT_HMULT_DU, UT_HMULT_DD )

      IF ( DO_ATMOS_LINEARIZATION ) THEN
        CALL L_HMULT_MASTER ( &
          DO_UPWELLING, DO_DNWELLING, &
          NLAYERS, N_USER_LEVELS, &
          N_USER_STREAMS, LOCAL_UM_START, &
          USER_SECANTS, &
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
          PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, &
          STERM_LAYERMASK_DN, &
          DELTAU_VERT, T_DELT_EIGEN, &
          T_UTUP_EIGEN, T_UTDN_EIGEN, &
          T_DELT_USERM, T_UTDN_USERM, &
          T_UTUP_USERM, &
          K_REAL, K_COMPLEX, HSINGO, &
          HMULT_1, HMULT_2, ZETA_M, ZETA_P, &
          DO_PROFILE_LINEARIZATION, DO_COLUMN_LINEARIZATION, &
          N_TOTALCOLUMN_WFS, LAYER_VARY_FLAG, &
          LAYER_VARY_NUMBER, &
          L_DELTAU_VERT, L_T_DELT_EIGEN, &
          L_T_UTUP_EIGEN, L_T_UTDN_EIGEN, &
          L_T_DELT_USERM, L_T_UTDN_USERM, &
          L_T_UTUP_USERM, L_KEIGEN, &
          L_ZETA_M, L_ZETA_P, &
          L_HMULT_1, L_HMULT_2, &
          L_UT_HMULT_UU, L_UT_HMULT_UD, &
          L_UT_HMULT_DU, L_UT_HMULT_DD )
      ENDIF

!  ############################################
!   boundary value problem - MATRIX PREPARATION
!  ############################################

!  standard case using compression of band matrices, etc..

      IF ( BVP_REGULAR_FLAG(FOURIER_COMPONENT) ) THEN

        CALL BVP_MATRIXSETUP_MASTER ( &
          DO_INCLUDE_SURFACE, FOURIER_COMPONENT, &
          SURFACE_FACTOR, &
          NSTOKES, NSTREAMS, &
          NLAYERS, DO_LAMBERTIAN_SURFACE, &
          LAMBERTIAN_ALBEDO, BRDF_F, &
          QUAD_STRMWTS, NSTKS_NSTRMS, &
          MUELLER_INDEX, &
          K_REAL, K_COMPLEX, &
          SOLA_XPOS, SOLB_XNEG, &
          NSTREAMS_2, NTOTAL, &
          N_SUBDIAG, N_SUPDIAG, &
          NSTKS_NSTRMS_2, &
          T_DELT_EIGEN, &
          R2_HOMP, R2_HOMM, AXBID_F, &
          BANDMAT2, IPIVOT, &
          SMAT2, SIPIVOT, &
          STATUS_SUB, MESSAGE, TRACE_1 )

        IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
          write(CF,'(I2)')FOURIER_COMPONENT
          TRACE_2 = 'Error from BVP_MATRIXSETUP_MASTER, '// &
            'Called in VLIDORT_L_FOURIER, Fourier # '//CF
          STATUS = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  Telescoped case

      ELSE

        CALL BVPTEL_MATRIXSETUP_MASTER ( &
          DO_INCLUDE_SURFACE, FOURIER_COMPONENT, &
          SURFACE_FACTOR, &
          NLAYERS, DO_SPECIALIST_OPTION_2, &
          NSTKS_NSTRMS, DO_LAYER_SCATTERING, &
          NSTOKES, NSTREAMS, &
          DO_LAMBERTIAN_SURFACE, &
          LAMBERTIAN_ALBEDO, BRDF_F, &
          QUAD_STRMWTS, &
          MUELLER_INDEX, K_REAL, K_COMPLEX, &
          SOLA_XPOS, SOLB_XNEG, &
          NSTREAMS_2, NSTKS_NSTRMS_2, T_DELT_EIGEN, &
          DO_BVTEL_INITIAL, &
          N_BVTELMATRIX_SIZE, N_BVTELMATRIX_SUPDIAG, &
          N_BVTELMATRIX_SUBDIAG, NLAYERS_TEL, &
          ACTIVE_LAYERS, &
          BANDTELMAT2, IPIVOTTEL, &
          R2_HOMP, R2_HOMM, AXBID_F, &
          SMAT2, SIPIVOT, &
          STATUS_SUB, MESSAGE, TRACE_1 )

        IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
          write(CF,'(I2)')FOURIER_COMPONENT
          TRACE_2 = 'Error from BVPTEL_MATRIXSETUP_MASTER, '// &
            'Called in VLIDORT_L_FOURIER, Fourier # '//CF
          STATUS = VLIDORT_SERIOUS
          RETURN
        ENDIF

      ENDIF

!  ################
!  Thermal Solution (2 steps)
!  ################

!  Continuation point for avoiding the scattering calculations

 8899 continue

!  THERMAL SOLUTIONS
!  =================

!  Separate calls if linearization is required

!  1. Find the Particular solution (also for transmittance only)
!  2. Compute thermal layer source terms. (Upwelling and Downwelling)
!    These will be scaled up by factor 4.pi if solar beams as well

      IF ( DO_INCLUDE_THERMEMISS ) THEN

!  With linearization

        IF ( DO_ATMOS_LINEARIZATION ) THEN
          CALL THERMAL_CLSOLUTION_PLUS ( &
            DO_UPWELLING, DO_DNWELLING, &
            DO_QUAD_OUTPUT, DO_ADDITIONAL_MVOUT, &
            NSTREAMS, NLAYERS, &
            N_THERMAL_COEFFS, DO_THERMAL_TRANSONLY, &
            DO_ATMOS_LINEARIZATION, LAYER_VARY_FLAG, &
            LAYER_VARY_NUMBER, &
            QUAD_STREAMS, QUAD_HALFWTS, &
            NMOMENTS, NSTREAMS_2, &
            N_USER_STREAMS, LOCAL_UM_START, &
            N_PARTLAYERS, &
            PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, &
            STERM_LAYERMASK_DN, &
            OMEGA_GREEK, T_DELT_DISORDS, &
            T_DISORDS_UTUP, T_DISORDS_UTDN, &
            PI_XQP, PI_XUP, PI_XUM, PI_XQM_PRE, &
            L_OMEGA_GREEK, L_T_DELT_DISORDS, &
            L_T_DISORDS_UTDN, L_T_DISORDS_UTUP, &
            SAB, DAB, &
            L_SAB, L_DAB, &
            THERMCOEFFS, DELTAU_POWER, &
            XTAU_POWER, TCOM1, &
            L_THERMCOEFFS, L_DELTAU_POWER, &
            L_XTAU_POWER, L_TCOM1, &
            T_WUPPER, T_WLOWER, &
            UT_T_PARTIC, U_TPOS1, &
            U_TNEG1, U_TPOS2, &
            U_TNEG2, &
            L_T_WUPPER, L_T_WLOWER, &
            L_UT_T_PARTIC, L_U_TPOS1, &
            L_U_TNEG1, L_U_TPOS2, &
            L_U_TNEG2, &
            STATUS_SUB, MESSAGE, TRACE_1 )

          IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
            TRACE_2 = 'Called in VLIDORT_L_FOURIER, Fourier 0'
            STATUS = VLIDORT_SERIOUS
            RETURN
          ENDIF

          IF ( DO_UPWELLING ) THEN
            CALL THERMAL_STERMS_UP_PLUS ( &
              DO_SOLAR_SOURCES, NLAYERS, &
              N_THERMAL_COEFFS, DO_THERMAL_TRANSONLY, &
              DO_ATMOS_LINEARIZATION, LAYER_VARY_FLAG, &
              LAYER_VARY_NUMBER, &
              N_USER_STREAMS, LOCAL_UM_START, &
              USER_STREAMS, DO_PARTLAYERS, &
              N_PARTLAYERS, PARTLAYERS_LAYERIDX, &
              N_ALLLAYERS_UP, STERM_LAYERMASK_UP, &
              T_DELT_USERM, T_UTUP_USERM, &
              L_T_DELT_USERM, L_T_UTUP_USERM, &
              DELTAU_POWER, XTAU_POWER, &
              U_TPOS1, U_TPOS2, &
              T_DIRECT_UP, T_UT_DIRECT_UP, &
              L_DELTAU_POWER, L_XTAU_POWER, &
              L_U_TPOS1, L_U_TPOS2, &
              L_T_DIRECT_UP, L_T_UT_DIRECT_UP, &
              LAYER_TSUP_UP, LAYER_TSUP_UTUP, &
              L_LAYER_TSUP_UP, L_LAYER_TSUP_UTUP )
          ENDIF

          IF ( DO_DNWELLING ) THEN
            CALL THERMAL_STERMS_DN_PLUS ( &
              DO_SOLAR_SOURCES, NLAYERS, &
              N_THERMAL_COEFFS, DO_THERMAL_TRANSONLY, &
              DO_ATMOS_LINEARIZATION, LAYER_VARY_FLAG, &
              LAYER_VARY_NUMBER, &
              N_USER_STREAMS, LOCAL_UM_START, &
              USER_STREAMS, DO_PARTLAYERS, &
              N_PARTLAYERS, PARTLAYERS_LAYERIDX, &
              N_ALLLAYERS_DN, STERM_LAYERMASK_DN, &
              T_DELT_USERM, T_UTDN_USERM, &
              L_T_DELT_USERM, L_T_UTDN_USERM, &
              DELTAU_POWER, XTAU_POWER, &
              U_TNEG1, U_TNEG2, &
              T_DIRECT_DN, T_UT_DIRECT_DN, &
              L_DELTAU_POWER, L_XTAU_POWER, &
              L_U_TNEG1, L_U_TNEG2, &
              L_T_DIRECT_DN, L_T_UT_DIRECT_DN, &
              LAYER_TSUP_DN, LAYER_TSUP_UTDN, &
              L_LAYER_TSUP_DN, L_LAYER_TSUP_UTDN )
          ENDIF

!  No linearization
!  @@@ Rob fix 1/31/11, This line was wrong (see below)
!          PI_XQP, PI_XQM, PI_XUP, PI_XQM_PRE, &

        ELSE

          CALL THERMAL_CLSOLUTION ( &
            DO_UPWELLING, DO_DNWELLING, &
            DO_QUAD_OUTPUT, DO_ADDITIONAL_MVOUT, &
            NSTREAMS, NLAYERS, &
            N_THERMAL_COEFFS, DO_THERMAL_TRANSONLY, &
            QUAD_STREAMS, QUAD_HALFWTS, &
            NMOMENTS, NSTREAMS_2, &
            N_USER_STREAMS, LOCAL_UM_START, &
            N_PARTLAYERS, &
            PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, &
            STERM_LAYERMASK_DN, &
            OMEGA_GREEK, T_DELT_DISORDS, &
            T_DISORDS_UTUP, T_DISORDS_UTDN, &
            PI_XQP, PI_XUP, PI_XUM, PI_XQM_PRE, &     ! Corrected 1/31/11
            SAB, DAB, &
            THERMCOEFFS, DELTAU_POWER, &
            XTAU_POWER, TCOM1, &
            T_WUPPER, T_WLOWER, &
            UT_T_PARTIC, U_TPOS1, &
            U_TNEG1, U_TPOS2, &
            U_TNEG2, &
            STATUS_SUB, MESSAGE, TRACE_1 )

          IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
            TRACE_2 = 'Called in VLIDORT_L_FOURIER, Fourier 0'
            STATUS = VLIDORT_SERIOUS
            RETURN
          ENDIF

          IF ( DO_UPWELLING ) THEN
            CALL THERMAL_STERMS_UP ( &
              DO_SOLAR_SOURCES, NLAYERS, &
              N_THERMAL_COEFFS, DO_THERMAL_TRANSONLY, &
              N_USER_STREAMS, LOCAL_UM_START, &
              USER_STREAMS, DO_PARTLAYERS, &
              N_PARTLAYERS, PARTLAYERS_LAYERIDX, &
              N_ALLLAYERS_UP, STERM_LAYERMASK_UP, &
              T_DELT_USERM, T_UTUP_USERM, &
              DELTAU_POWER, XTAU_POWER, &
              U_TPOS1, U_TPOS2, &
              T_DIRECT_UP, T_UT_DIRECT_UP, &
              LAYER_TSUP_UP, LAYER_TSUP_UTUP )
          ENDIF

          IF ( DO_DNWELLING ) THEN
            CALL THERMAL_STERMS_DN ( &
              DO_SOLAR_SOURCES, NLAYERS, &
              N_THERMAL_COEFFS, DO_THERMAL_TRANSONLY, &
              N_USER_STREAMS, LOCAL_UM_START, &
              USER_STREAMS, DO_PARTLAYERS, &
              N_PARTLAYERS, PARTLAYERS_LAYERIDX, &
              N_ALLLAYERS_DN, STERM_LAYERMASK_DN, &
              T_DELT_USERM, T_UTDN_USERM, &
              DELTAU_POWER, XTAU_POWER, &
              U_TNEG1, U_TNEG2, &
              T_DIRECT_DN, T_UT_DIRECT_DN, &
              LAYER_TSUP_DN, LAYER_TSUP_UTDN )
          ENDIF

        ENDIF

!  End include thermal emission

      ENDIF

!  Skip the thermal-only section if there are solar sources

      IF ( DO_SOLAR_SOURCES ) GO TO 455

!  ####################################################
!  Complete Radiation Field with Thermal-only solutions
!  ####################################################

!    Only one solution, local direct_beam flag NOT set

      IPARTIC = 1
      DO_INCLUDE_DIRECTBEAM   = .FALSE.

!  Avoid the scattering solutions if not flagged

      IF ( DO_THERMAL_TRANSONLY ) GO TO 566

!  Thermal-only. Find the BVP solution and intensity field
!  -------------------------------------------------------

!  set the BVP PI solution at the lower/upper boundaries

!mick fix 2/17/11 - include remaining stokes components as needed
!      O1 = 1
!      DO LAYER = 1, NLAYERS
!        DO I = 1, NSTREAMS_2
!          WUPPER(I,O1,LAYER) = T_WUPPER(I,LAYER)
!          WLOWER(I,O1,LAYER) = T_WLOWER(I,LAYER)
!        ENDDO
!      ENDDO
      IF (NSTOKES .EQ. 1) THEN
        O1 = 1
        DO LAYER = 1, NLAYERS
          DO I = 1, NSTREAMS_2
            WUPPER(I,O1,LAYER) = T_WUPPER(I,LAYER)
            WLOWER(I,O1,LAYER) = T_WLOWER(I,LAYER)
          ENDDO
        ENDDO
      ELSE
        DO LAYER = 1, NLAYERS
          DO O1 = 1, NSTOKES
            DO I = 1, NSTREAMS_2
              WUPPER(I,O1,LAYER) = T_WUPPER(I,LAYER)
              WLOWER(I,O1,LAYER) = T_WLOWER(I,LAYER)
            ENDDO
          ENDDO
        ENDDO
      END IF

!  Solve the boundary value problem. No telescoping here

      CALL BVP_SOLUTION_MASTER ( &
        DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM, &
        DO_INCLUDE_SURFEMISS, FOURIER_COMPONENT, &
        IPARTIC, SURFACE_FACTOR, &
        NSTOKES, NSTREAMS, &
        NLAYERS, DO_LAMBERTIAN_SURFACE, &
        LAMBERTIAN_ALBEDO, &
        QUAD_STRMWTS, MUELLER_INDEX, &
        WLOWER, &
        EMISSIVITY, SURFBB, &
        NSTREAMS_2, NTOTAL, &
        NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
        DIRECT_BEAM, &
        WUPPER, &
        DO_DEBUG_WRITE, DO_FDTEST, &
        N_SUBDIAG, N_SUPDIAG, &
        K_REAL, K_COMPLEX, &
        BANDMAT2, IPIVOT, &
        SMAT2, SIPIVOT, &
        AXBID_F, &
        R2_BEAM, LCON, MCON, &
        STATUS_SUB, MESSAGE, TRACE_1 )

!  Exception handling

      IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
        write(CF,'(I2)')FOURIER_COMPONENT
        TRACE_2 = 'Error return from BVP_SOLUTION_MASTER, '// &
           'Called in VLIDORT_L_FOURIER, Fourier # '//CF
        STATUS = VLIDORT_SERIOUS
        RETURN
      ENDIF

!  Continuation point for avoiding thermal scattering

 566  CONTINUE

!  Post-processing - upwelling thermal-only field

      IF ( DO_UPWELLING ) THEN
        CALL VLIDORT_UPUSER_INTENSITY ( &
          DO_INCLUDE_SURFACE, DO_INCLUDE_SURFEMISS, &
          DO_INCLUDE_THERMEMISS, DO_INCLUDE_DIRECTBEAM, &
          DO_INCLUDE_MVOUTPUT, FOURIER_COMPONENT, IPARTIC, &
          FLUX_MULTIPLIER, SURFACE_FACTOR, DO_DEBUG_WRITE, DO_FDTEST, &
          NSTOKES, NLAYERS, N_USER_LEVELS, DO_TOA_CONTRIBS, &
          N_USER_STREAMS, DO_LAYER_SCATTERING, &
          DO_USER_STREAMS, LOCAL_UM_START, &
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX, &
          UTAU_LEVEL_MASK_UP, T_DELT_USERM, T_UTUP_USERM, &
          CUMTRANS, &
          DO_QUAD_OUTPUT, NSTREAMS, DO_LAMBERTIAN_SURFACE, &
          LAMBERTIAN_ALBEDO, USER_BRDF_F, EMISSIVITY, USER_EMISSIVITY, &
          SURFBB, DO_THERMAL_TRANSONLY, DO_DBCORRECTION, QUAD_WEIGHTS, &
          QUAD_STRMWTS, MUELLER_INDEX, K_REAL, K_COMPLEX, &
          SOLA_XPOS, SOLB_XNEG, WLOWER, LCON, MCON, &
          USER_DIRECT_BEAM, &
          T_DELT_DISORDS, T_DELT_EIGEN, T_WLOWER, &
          DO_SOLAR_SOURCES, DO_CLASSICAL_SOLUTION, &
          DO_MSMODE_VLIDORT, DO_MSMODE_THERMAL, &
          UHOM_UPDN, UHOM_UPUP, UPAR_UP_1, UPAR_UP_2, &
          HMULT_1, HMULT_2, EMULT_UP, LAYER_TSUP_UP, &
          UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP, LAYER_TSUP_UTUP, BRDF_F, &
          STOKES_DOWNSURF, &
          CUMSOURCE_UP, STOKES_F, MS_CONTRIBS_F, BOA_THTONLY_SOURCE )
      ENDIF

!  Post-processing - Downwelling thermal-only field

      IF ( DO_DNWELLING ) THEN
        CALL VLIDORT_DNUSER_INTENSITY ( &
          DO_INCLUDE_THERMEMISS, FOURIER_COMPONENT, &
          IPARTIC, FLUX_MULTIPLIER, &
          DO_DEBUG_WRITE, NSTOKES, &
          N_USER_LEVELS, &
          N_USER_STREAMS, DO_LAYER_SCATTERING, &
          DO_USER_STREAMS, LOCAL_UM_START, &
          PARTLAYERS_OUTFLAG, &
          PARTLAYERS_OUTINDEX, UTAU_LEVEL_MASK_DN, &
          PARTLAYERS_LAYERIDX, &
          T_DELT_USERM, T_UTDN_USERM, &
          DO_SOLAR_SOURCES, DO_THERMAL_TRANSONLY, &
          DO_CLASSICAL_SOLUTION, DO_MSMODE_VLIDORT, &
          K_REAL, K_COMPLEX, &
          LCON, MCON, &
          UHOM_DNDN, UHOM_DNUP, &
          UPAR_DN_1, UPAR_DN_2, &
          HMULT_1, HMULT_2, &
          EMULT_DN, LAYER_TSUP_DN, &
          UT_HMULT_DU, UT_HMULT_DD, &
          UT_EMULT_DN, LAYER_TSUP_UTDN, &
          CUMSOURCE_DN, STOKES_F )
      ENDIF

!  mean value (integrated) output for thermal-only field
!    ------- Can also use this to get debug Quadrature output
! @@@ Rob fix 1/31/11, - added FLUX_FACTOR argument (3rd line)

      IF ( DO_INCLUDE_MVOUTPUT.OR.DO_QUAD_OUTPUT ) THEN
        CALL VLIDORT_INTEGRATED_OUTPUT ( &
          DO_INCLUDE_MVOUTPUT, DO_INCLUDE_DIRECTBEAM, &
          DO_INCLUDE_THERMEMISS, &
          IPARTIC, FLUX_MULTIPLIER, FLUX_FACTOR, &    !  added here
          NSTOKES, NSTREAMS, &
          N_USER_LEVELS, &
          FLUXVEC, LAYER_PIS_CUTOFF, &
          QUAD_WEIGHTS, QUAD_STRMWTS, &
          N_DIRECTIONS, WHICH_DIRECTIONS, &
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
          UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, &
          PARTLAYERS_LAYERIDX, &
          T_DELT_MUBAR, T_UTDN_MUBAR, &
          INITIAL_TRANS, LOCAL_CSZA, &
          DO_SOLAR_SOURCES, NLAYERS, &
          DO_THERMAL_TRANSONLY, &
          DO_CLASSICAL_SOLUTION, QUAD_STREAMS, &
          T_DELT_DISORDS, T_DISORDS_UTUP, &
          T_UTUP_EIGEN, T_UTDN_EIGEN, &
          K_REAL, K_COMPLEX, &
          SOLA_XPOS, SOLB_XNEG, &
          WUPPER, LCON, MCON, &
          T_WUPPER, UT_T_PARTIC, &
          BOA_THTONLY_SOURCE, &
          T_DELT_EIGEN, WLOWER, &
          T_DISORDS_UTDN, T_WLOWER, &
          MEAN_STOKES, FLUX_STOKES, &
          MEAN_DIRECT, FLUX_DIRECT )
      ENDIF

!  LTE linearization (new section, 14 September 2009)
!  --------------------------------------------------

      IF ( DO_THERMAL_TRANSONLY .AND. DO_LTE_LINEARIZATION ) THEN
        CALL THERMAL_LTE_LINEARIZATION ( &
          DO_INCLUDE_SURFACE, SURFACE_FACTOR, &
          FLUX_MULTIPLIER, &
          DO_UPWELLING, DO_DNWELLING, &
          NSTREAMS, NLAYERS, N_USER_LEVELS, &
          DELTAU_VERT_INPUT, LAMBERTIAN_ALBEDO, &
          THERMAL_BB_INPUT, &
          QUAD_STREAMS, QUAD_STRMWTS, &
          N_USER_STREAMS, LOCAL_UM_START, USER_STREAMS, &
          UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, &
          STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, &
          T_DELT_DISORDS, T_DELT_USERM, &
          CUMSOURCE_UP, CUMSOURCE_DN, &
          THERMCOEFFS, &
          LTE_DELTAU_VERT_INPUT, LTE_THERMAL_BB_INPUT, &
          LTE_ATMOSWF )
      ENDIF

!  Thermal-only: Atmospheric profile linearization
!  -----------------------------------------------

      IF ( DO_PROFILE_LINEARIZATION ) THEN

!  Start over layers that will contain variations
!    - Set flag for variation of layer, and number of variations

        DO LAYER_TO_VARY = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(LAYER_TO_VARY) ) THEN
            N_LAYER_WFS = LAYER_VARY_NUMBER(LAYER_TO_VARY)

!  Avoidance of BVP problem for transmittance only

            IF ( DO_THERMAL_TRANSONLY ) GO TO 788

!  Solve the Regular BVP linearization

            IF ( BVP_REGULAR_FLAG(FOURIER_COMPONENT) ) THEN

!  Get the BVP solution (regular case) for this beam component

              CALL L_BVP_SOLUTION_MASTER ( &
                DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM, &
                DO_INCLUDE_THERMEMISS, &
                LAYER_TO_VARY, N_LAYER_WFS, &
                FOURIER_COMPONENT, IPARTIC, &
                SURFACE_FACTOR, DO_PROFILE_LINEARIZATION, &
                DO_COLUMN_LINEARIZATION, &
                DO_SOLAR_SOURCES, NSTOKES, &
                NSTREAMS, NLAYERS, NSTREAMS_2, &
                NTOTAL, NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
                DELTAU_SLANT, T_DELT_EIGEN, &
                K_REAL, K_COMPLEX, &
                SOLA_XPOS, SOLB_XNEG, &
                LCON, MCON, DIRECT_BEAM, &
                L_DELTAU_VERT, L_T_DELT_EIGEN, &
                L_SOLA_XPOS, L_SOLB_XNEG, &
                L_T_WUPPER, L_T_WLOWER, &
                DO_CLASSICAL_SOLUTION, LAYER_PIS_CUTOFF, &
                DO_LAYER_SCATTERING, &
                T_DELT_MUBAR, INITIAL_TRANS, BVEC, &
                L_INITIAL_TRANS, L_T_DELT_MUBAR, L_BVEC, &
                DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
                BRDF_F, QUAD_STRMWTS, MUELLER_INDEX, &
                LAYER_TO_VARY, N_LAYER_WFS, &
                DO_PLANE_PARALLEL, &
                N_SUBDIAG, N_SUPDIAG, &
                BANDMAT2, IPIVOT, &
                SMAT2, SIPIVOT, &
                L_WLOWER, L_WUPPER, NCON, PCON, &
                STATUS_SUB, MESSAGE, TRACE_1 )

!  Exception handling

              IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
                write(CF,'(I2)')FOURIER_COMPONENT
                TRACE_2 = 'Error return from L_BVP_SOLUTION_MASTER, '// &
        'Profile Jacobians, Called in VLIDORT_L_FOURIER, Fourier # '//CF
                STATUS = VLIDORT_SERIOUS
                RETURN
              ENDIF

            ENDIF

!  Continuation point for avoiding scattering solution

 788        CONTINUE

!  Post-processing for the weighting functions, Sequence:
!    Upwelling   Atmospheric weighting functions
!    Downwelling Atmospheric weighting functions
!    Mean-value  Atmospheric weighting functions


            IF ( DO_UPWELLING ) THEN
              CALL VLIDORT_UPUSER_ATMOSWF ( &
                DO_INCLUDE_SURFACE, DO_INCLUDE_THERMEMISS, &
                DO_INCLUDE_MVOUTPUT, DO_INCLUDE_DIRECTBEAM, &
                SURFACE_FACTOR, FLUX_MULTIPLIER, FOURIER_COMPONENT, &
                IPARTIC, LAYER_TO_VARY, N_LAYER_WFS, &
                NSTOKES, NLAYERS, N_USER_LEVELS, &
                N_USER_STREAMS, DO_LAYER_SCATTERING, &
                DO_USER_STREAMS, LOCAL_UM_START, &
                PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
                UTAU_LEVEL_MASK_UP, PARTLAYERS_LAYERIDX, &
                T_DELT_USERM, T_UTUP_USERM, CUMSOURCE_UP, &
                L_T_DELT_USERM, L_T_UTUP_USERM, &
                DO_SOLAR_SOURCES, DO_QUAD_OUTPUT, &
                NSTREAMS, DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
                BRDF_F, USER_BRDF_F, DO_THERMAL_TRANSONLY, DO_DBCORRECTION, &
                QUAD_WEIGHTS, QUAD_STRMWTS, MUELLER_INDEX, &
                DELTAU_SLANT, T_DELT_DISORDS, T_DELT_EIGEN, &
                K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, &
                LCON, MCON, T_WLOWER, USER_DIRECT_BEAM, &
                L_DELTAU_VERT, L_T_DELT_EIGEN, L_T_DELT_DISORDS, &
                L_SOLA_XPOS, L_SOLB_XNEG, &
                L_WLOWER, NCON, PCON, L_T_WLOWER, &
                DO_MSMODE_VLIDORT, &
                UHOM_UPDN, UHOM_UPUP, UPAR_UP_1, UPAR_UP_2, &
                HMULT_1, HMULT_2, EMULT_UP, &
                L_UHOM_UPDN, L_UHOM_UPUP, L_UPAR_UP_1, L_UPAR_UP_2, &
                L_HMULT_1, L_HMULT_2, L_EMULT_UP, &
                L_LAYER_TSUP_UP, &
                UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP, &
                L_UT_HMULT_UU, L_UT_HMULT_UD, L_UT_EMULT_UP, &
                L_LAYER_TSUP_UTUP, &
                L_BOA_THTONLY_SOURCE, ATMOSWF_F )
            ENDIF

            IF ( DO_DNWELLING ) THEN
              CALL VLIDORT_DNUSER_ATMOSWF ( &
                DO_INCLUDE_THERMEMISS, FLUX_MULTIPLIER, &
                FOURIER_COMPONENT, IPARTIC, &
                LAYER_TO_VARY, N_LAYER_WFS, &
                NSTOKES, N_USER_LEVELS, &
                N_USER_STREAMS, DO_LAYER_SCATTERING, &
                DO_USER_STREAMS, LOCAL_UM_START, &
                PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
                UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX, &
                T_DELT_USERM, T_UTDN_USERM, CUMSOURCE_DN, &
                L_T_DELT_USERM, L_T_UTDN_USERM, &
                DO_SOLAR_SOURCES, DO_THERMAL_TRANSONLY, &
                DO_MSMODE_VLIDORT, &
                K_REAL, K_COMPLEX, &
                LCON, MCON, &
                UHOM_DNDN, UHOM_DNUP, &
                UPAR_DN_1, UPAR_DN_2, &
                HMULT_1, HMULT_2, EMULT_DN, &
                NCON, PCON, &
                L_UHOM_DNDN, L_UHOM_DNUP, &
                L_UPAR_DN_1, L_UPAR_DN_2, &
                L_HMULT_1, L_HMULT_2, L_EMULT_DN, &
                L_LAYER_TSUP_DN, &
                UT_HMULT_DU, UT_HMULT_DD, UT_EMULT_DN, &
                L_UT_HMULT_DU, L_UT_HMULT_DD, &
                L_UT_EMULT_DN, &
                L_LAYER_TSUP_UTDN, &
                ATMOSWF_F )
            ENDIF

! @@@ Rob fix 1/31/11, - added FLUX_FACTOR argument

            IF ( DO_INCLUDE_MVOUTPUT .OR. DO_QUAD_OUTPUT ) THEN
              CALL VLIDORT_L_INTEGRATED_OUTPUT ( &
                DO_INCLUDE_MVOUTPUT, DO_INCLUDE_DIRECTBEAM, &
                DO_INCLUDE_THERMEMISS, FLUX_MULTIPLIER, FLUX_FACTOR, & ! added
                IPARTIC, LAYER_TO_VARY, N_LAYER_WFS, &
                NSTOKES, NSTREAMS, N_USER_LEVELS, &
                FLUXVEC, LAYER_PIS_CUTOFF, &
                QUAD_WEIGHTS, QUAD_STRMWTS, &
                N_DIRECTIONS, WHICH_DIRECTIONS, &
                PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
                UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, &
                PARTLAYERS_LAYERIDX, &
                T_DELT_MUBAR, T_UTDN_MUBAR, &
                INITIAL_TRANS, LOCAL_CSZA, &
                L_INITIAL_TRANS, L_T_DELT_MUBAR, &
                L_T_UTDN_MUBAR, &
                DO_SOLAR_SOURCES, &
                NLAYERS, DO_THERMAL_TRANSONLY, &
                DO_CLASSICAL_SOLUTION, QUAD_STREAMS, &
                T_DELT_DISORDS, T_DISORDS_UTUP, &
                T_UTUP_EIGEN, T_UTDN_EIGEN, &
                K_REAL, K_COMPLEX, &
                SOLA_XPOS, SOLB_XNEG, &
                WUPPER, LCON, MCON, T_WUPPER, &
                BOA_THTONLY_SOURCE, &
                L_T_UTUP_EIGEN, L_T_UTDN_EIGEN, &
                L_T_DELT_DISORDS, L_T_DISORDS_UTUP, &
                L_SOLA_XPOS, L_SOLB_XNEG, &
                L_WUPPER, NCON, PCON, &
                L_T_WUPPER, L_UT_T_PARTIC, &
                L_BOA_THTONLY_SOURCE, &
                T_DELT_EIGEN, L_T_DELT_EIGEN, &
                L_WLOWER, &
                T_DISORDS_UTDN, L_T_DISORDS_UTDN, &
                L_T_WLOWER, T_WLOWER, &
                MINT_ATMOSWF, MINT_ATMOSWF_DIRECT, &
                FLUX_ATMOSWF, FLUX_ATMOSWF_DIRECT )
            ENDIF

!  Finish loop over layers with variation

          ENDIF
        ENDDO

!  End atmospheric profile weighting functions

      ENDIF

!  Thermal Only: Atmospheric Bulk weighting functions
!  --------------------------------------------------

      IF ( DO_COLUMN_LINEARIZATION ) THEN

!   variation index = 0

        VARIATION_INDEX = 0

!  Avoidance of BVP problem for transmittance only

        IF ( DO_THERMAL_TRANSONLY ) GO TO 789

!  Get the linearized BVP solution for this beam component

        CALL L_BVP_SOLUTION_MASTER ( &
          DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM, &
          DO_INCLUDE_THERMEMISS, &
          VARIATION_INDEX, N_TOTALCOLUMN_WFS, &
          FOURIER_COMPONENT, IBEAM, &
          SURFACE_FACTOR, DO_PROFILE_LINEARIZATION, &
          DO_COLUMN_LINEARIZATION, &
          DO_SOLAR_SOURCES, NSTOKES, &
          NSTREAMS, NLAYERS, NSTREAMS_2, &
          NTOTAL, NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
          DELTAU_SLANT, T_DELT_EIGEN, &
          K_REAL, K_COMPLEX, &
          SOLA_XPOS, SOLB_XNEG, &
          LCON, MCON, DIRECT_BEAM, &
          L_DELTAU_VERT, L_T_DELT_EIGEN, &
          L_SOLA_XPOS, L_SOLB_XNEG, &
          L_T_WUPPER, L_T_WLOWER, &
          DO_CLASSICAL_SOLUTION, LAYER_PIS_CUTOFF, &
          DO_LAYER_SCATTERING, &
          T_DELT_MUBAR, INITIAL_TRANS, BVEC, &
          L_INITIAL_TRANS, L_T_DELT_MUBAR, L_BVEC, &
          DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
          BRDF_F, QUAD_STRMWTS, MUELLER_INDEX, &
          LAYER_TO_VARY, N_LAYER_WFS, &
          DO_PLANE_PARALLEL, &
          N_SUBDIAG, N_SUPDIAG, &
          BANDMAT2, IPIVOT, &
          SMAT2, SIPIVOT, &
          L_WLOWER, L_WUPPER, NCON, PCON, &
          STATUS_SUB, MESSAGE, TRACE_1 )

        IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
          write(CF,'(I2)')FOURIER_COMPONENT
          TRACE_2 = 'Error return from L_BVP_SOLUTION_MASTER, '// &
          'Column Jacobians, Called in VLIDORT_L_FOURIER, Fourier # '//CF
          STATUS = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  Continuation point for avoiding scattering solution

 789    CONTINUE

!  Post-processing for the weighting functions, Sequence:
!    Upwelling   Atmospheric weighting functions
!    Downwelling Atmospheric weighting functions
!    Mean-value  Atmospheric weighting functions

        IF ( DO_UPWELLING ) THEN
          CALL VLIDORT_UPUSER_ATMOSWF ( &
            DO_INCLUDE_SURFACE, DO_INCLUDE_THERMEMISS, &
            DO_INCLUDE_MVOUTPUT, DO_INCLUDE_DIRECTBEAM, &
            SURFACE_FACTOR, FLUX_MULTIPLIER, FOURIER_COMPONENT, &
            IPARTIC, VARIATION_INDEX, N_TOTALCOLUMN_WFS, &
            NSTOKES, NLAYERS, N_USER_LEVELS, &
            N_USER_STREAMS, DO_LAYER_SCATTERING, &
            DO_USER_STREAMS, LOCAL_UM_START, &
            PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
            UTAU_LEVEL_MASK_UP, PARTLAYERS_LAYERIDX, &
            T_DELT_USERM, T_UTUP_USERM, CUMSOURCE_UP, &
            L_T_DELT_USERM, L_T_UTUP_USERM, &
            DO_SOLAR_SOURCES, DO_QUAD_OUTPUT, &
            NSTREAMS, DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
            BRDF_F, USER_BRDF_F, DO_THERMAL_TRANSONLY, DO_DBCORRECTION, &
            QUAD_WEIGHTS, QUAD_STRMWTS, MUELLER_INDEX, &
            DELTAU_SLANT, T_DELT_DISORDS, T_DELT_EIGEN, &
            K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, &
            LCON, MCON, T_WLOWER, USER_DIRECT_BEAM, &
            L_DELTAU_VERT, L_T_DELT_EIGEN, L_T_DELT_DISORDS, &
            L_SOLA_XPOS, L_SOLB_XNEG, &
            L_WLOWER, NCON, PCON, L_T_WLOWER, &
            DO_MSMODE_VLIDORT, &
            UHOM_UPDN, UHOM_UPUP, UPAR_UP_1, UPAR_UP_2, &
            HMULT_1, HMULT_2, EMULT_UP, &
            L_UHOM_UPDN, L_UHOM_UPUP, L_UPAR_UP_1, L_UPAR_UP_2, &
            L_HMULT_1, L_HMULT_2, L_EMULT_UP, &
            L_LAYER_TSUP_UP, &
            UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP, &
            L_UT_HMULT_UU, L_UT_HMULT_UD, L_UT_EMULT_UP, &
            L_LAYER_TSUP_UTUP, &
            L_BOA_THTONLY_SOURCE, ATMOSWF_F )
        ENDIF

        IF ( DO_DNWELLING ) THEN
          CALL VLIDORT_DNUSER_ATMOSWF ( &
            DO_INCLUDE_THERMEMISS, FLUX_MULTIPLIER, &
            FOURIER_COMPONENT, IPARTIC, &
            VARIATION_INDEX, N_TOTALCOLUMN_WFS, &
            NSTOKES, N_USER_LEVELS, &
            N_USER_STREAMS, DO_LAYER_SCATTERING, &
            DO_USER_STREAMS, LOCAL_UM_START, &
            PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
            UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX, &
            T_DELT_USERM, T_UTDN_USERM, CUMSOURCE_DN, &
            L_T_DELT_USERM, L_T_UTDN_USERM, &
            DO_SOLAR_SOURCES, DO_THERMAL_TRANSONLY, &
            DO_MSMODE_VLIDORT, &
            K_REAL, K_COMPLEX, &
            LCON, MCON, &
            UHOM_DNDN, UHOM_DNUP, &
            UPAR_DN_1, UPAR_DN_2, &
            HMULT_1, HMULT_2, EMULT_DN, &
            NCON, PCON, &
            L_UHOM_DNDN, L_UHOM_DNUP, &
            L_UPAR_DN_1, L_UPAR_DN_2, &
            L_HMULT_1, L_HMULT_2, L_EMULT_DN, &
            L_LAYER_TSUP_DN, &
            UT_HMULT_DU, UT_HMULT_DD, UT_EMULT_DN, &
            L_UT_HMULT_DU, L_UT_HMULT_DD, &
            L_UT_EMULT_DN, &
            L_LAYER_TSUP_UTDN, &
            ATMOSWF_F )
        ENDIF

! @@@ Rob fix 1/31/11, - added FLUX_FACTOR argument

        IF ( DO_INCLUDE_MVOUTPUT .OR. DO_QUAD_OUTPUT ) THEN
          CALL VLIDORT_L_INTEGRATED_OUTPUT ( &
            DO_INCLUDE_MVOUTPUT, DO_INCLUDE_DIRECTBEAM, &
            DO_INCLUDE_THERMEMISS, FLUX_MULTIPLIER, FLUX_FACTOR, & ! added
            IPARTIC, VARIATION_INDEX, N_TOTALCOLUMN_WFS, &
            NSTOKES, NSTREAMS, N_USER_LEVELS, &
            FLUXVEC, LAYER_PIS_CUTOFF, &
            QUAD_WEIGHTS, QUAD_STRMWTS, &
            N_DIRECTIONS, WHICH_DIRECTIONS, &
            PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
            UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, &
            PARTLAYERS_LAYERIDX, &
            T_DELT_MUBAR, T_UTDN_MUBAR, &
            INITIAL_TRANS, LOCAL_CSZA, &
            L_INITIAL_TRANS, L_T_DELT_MUBAR, &
            L_T_UTDN_MUBAR, &
            DO_SOLAR_SOURCES, &
            NLAYERS, DO_THERMAL_TRANSONLY, &
            DO_CLASSICAL_SOLUTION, QUAD_STREAMS, &
            T_DELT_DISORDS, T_DISORDS_UTUP, &
            T_UTUP_EIGEN, T_UTDN_EIGEN, &
            K_REAL, K_COMPLEX, &
            SOLA_XPOS, SOLB_XNEG, &
            WUPPER, LCON, MCON, T_WUPPER, &
            BOA_THTONLY_SOURCE, &
            L_T_UTUP_EIGEN, L_T_UTDN_EIGEN, &
            L_T_DELT_DISORDS, L_T_DISORDS_UTUP, &
            L_SOLA_XPOS, L_SOLB_XNEG, &
            L_WUPPER, NCON, PCON, &
            L_T_WUPPER, L_UT_T_PARTIC, &
            L_BOA_THTONLY_SOURCE, &
            T_DELT_EIGEN, L_T_DELT_EIGEN, &
            L_WLOWER, &
            T_DISORDS_UTDN, L_T_DISORDS_UTDN, &
            L_T_WLOWER, T_WLOWER, &
            MINT_ATMOSWF, MINT_ATMOSWF_DIRECT, &
            FLUX_ATMOSWF, FLUX_ATMOSWF_DIRECT )
        ENDIF

!  End atmospheric column weighting functions

      ENDIF

!  Thermal-only: Surface Reflectance weighting functions
!  -----------------------------------------------------

      IF ( DO_SURFACE_LINEARIZATION ) THEN

!  Only do this if there is surface reflectance !

        IF ( DO_INCLUDE_SURFACE ) THEN

          CALL VLIDORT_SURFACE_WFS ( &
            DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFEMISS, &
            DO_MSMODE_THERMAL, &
            DO_INCLUDE_MVOUTPUT, N_SURFACE_WFS, &
            FOURIER_COMPONENT, IPARTIC, SURFACE_FACTOR, &
            FLUX_MULTIPLIER, DO_UPWELLING, DO_DNWELLING, &
            DO_QUAD_OUTPUT, NLAYERS, NTOTAL, &
            N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
            K_REAL, K_COMPLEX, &
            BANDMAT2, IPIVOT, SMAT2, SIPIVOT, &
            NSTOKES, NSTREAMS, &
            DO_LAMBERTIAN_SURFACE, SURFBB, LS_BRDF_F, &
            LS_BRDF_F_0, LS_EMISSIVITY, QUAD_STRMWTS, &
            MUELLER_INDEX, T_DELT_EIGEN, SOLA_XPOS, SOLB_XNEG, &
            WLOWER, LCON, MCON, ATMOS_ATTN, &
            R2_HOMP, R2_HOMM, R2_BEAM, DIRECT_BEAM, &
            LAMBERTIAN_ALBEDO, USER_BRDF_F, &
            DO_THERMAL_TRANSONLY, LS_USER_BRDF_F, LS_USER_BRDF_F_0, &
            LS_USER_EMISSIVITY, &
            DO_DBCORRECTION, FLUXVEC, &
            N_USER_STREAMS, DO_USER_STREAMS, LOCAL_UM_START, &
            STOKES_DOWNSURF, N_USER_LEVELS, &
            PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
            UTAU_LEVEL_MASK_UP, PARTLAYERS_LAYERIDX, &
            T_DELT_USERM, T_UTUP_USERM, &
            UHOM_UPDN, UHOM_UPUP, HMULT_1, HMULT_2, &
            UT_HMULT_UU, UT_HMULT_UD, &
            UTAU_LEVEL_MASK_DN, T_UTDN_USERM, &
            UHOM_DNDN, UHOM_DNUP, UT_HMULT_DU, UT_HMULT_DD, &
            QUAD_WEIGHTS, N_DIRECTIONS, WHICH_DIRECTIONS, &
            T_DELT_DISORDS, T_DISORDS_UTUP, &
            T_UTUP_EIGEN, T_UTDN_EIGEN, &
            SURFACEWF_F, MINT_SURFACEWF, FLUX_SURFACEWF, &
            STATUS_SUB, MESSAGE, TRACE_1 )

!  Exception handling

          IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
            write(CF,'(I2)')FOURIER_COMPONENT
            TRACE_2 = 'Error return from VLIDORT_SURFACE_WFS, '// &
            'Thermal-only Called in VLIDORT_L_FOURIER, Fourier # '//CF
            STATUS = VLIDORT_SERIOUS
            RETURN
          ENDIF

        ENDIF

!  end of surface weighting functions

      ENDIF

!  Finish Thermal only.

      RETURN

!  ##################################################
!  Complete Radiation Field with Solar Beam solutions
!  ##################################################

!  Continuation point

 455  CONTINUE

!  Start loop over desired beam solutions

      DO IBEAM = 1, NBEAMS

        IF ( DO_MULTIBEAM(IBEAM,FOURIER_COMPONENT) ) THEN

!  Step 1. Solar beam Particular solutions + linearizations
!  --------------------------------------------------------

          DO LAYER = 1, NLAYERS

!  Parameter control for the linearization

            IF ( DO_ATMOS_LINEARIZATION ) THEN
              IF ( DO_PROFILE_LINEARIZATION ) THEN
                DO_RTSOL_VARY = LAYER_VARY_FLAG(LAYER)
                NPARAMS_VARY  = LAYER_VARY_NUMBER(LAYER)
              ELSE IF ( DO_COLUMN_LINEARIZATION ) THEN
                DO_RTSOL_VARY = .TRUE.
                NPARAMS_VARY  = N_TOTALCOLUMN_WFS
              ENDIF
            ENDIF

!  A. the classical solution
!     ++++++++++++++++++++++

            IF ( DO_CLASSICAL_SOLUTION ) THEN

!  For the particular integral itself

              CALL VLIDORT_QBEAM_SOLUTION ( &
                LAYER, FOURIER_COMPONENT, IBEAM, &
                DO_DEBUG_WRITE, DO_FDTEST, &
                NSTOKES, NSTREAMS, FLUX_FACTOR, &
                LAYER_PIS_CUTOFF, QUAD_STREAMS, &
                NMOMENTS, NSTREAMS_2, &
                NSTKS_NSTRMS, DO_LAYER_SCATTERING, &
                DMAT, DFLUX, &
                OMEGA_GREEK, T_DELT_MUBAR, &
                INITIAL_TRANS, AVERAGE_SECANT, &
                PI_XQP, PI_XQM, &
                PI_X0P, PI_XQM_POST, &
                SAB, DAB, &
                EIGENMAT_SAVE, &
                QSUMVEC_SAVE, QDIFVEC_SAVE, &
                QVEC_SAVE, QDIF_SAVE, &
                QMAT_SAVE, QPIVOT, &
                BVEC, WUPPER, WLOWER, &
                STATUS_SUB, MESSAGE, TRACE_1 )

!  error flag

              IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
                write(CF,'(I2)')FOURIER_COMPONENT
                TRACE_2 = 'Error return from VLIDORT_QBEAM_SOLUTION, '// &
                  'Called in VLIDORT_L_FOURIER, Fourier # '//CF
                STATUS = VLIDORT_SERIOUS
                RETURN
              ENDIF

!  For the linearizations of the classical solution

              IF ( DO_ATMOS_LINEARIZATION ) THEN

                CALL VLIDORT_L_QBEAM_SOLUTION ( &
                  LAYER, FOURIER_COMPONENT, IBEAM, &
                  DO_RTSOL_VARY, NPARAMS_VARY, &
                  DO_PLANE_PARALLEL, NSTOKES, &
                  NSTREAMS, FLUX_FACTOR, &
                  LAYER_PIS_CUTOFF, QUAD_STREAMS, &
                  NMOMENTS, NSTREAMS_2, &
                  NSTKS_NSTRMS, DO_LAYER_SCATTERING, &
                  DMAT, DFLUX, &
                  AVERAGE_SECANT, &
                  PI_XQP, PI_X0P, PI_XQM_POST, &
                  SAB, DAB, &
                  QSUMVEC_SAVE, QDIFVEC_SAVE, &
                  QVEC_SAVE, QDIF_SAVE, &
                  QMAT_SAVE, QPIVOT, &
                  DO_PROFILE_LINEARIZATION, DO_COLUMN_LINEARIZATION, &
                  LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
                  L_SAB, L_DAB, L_EIGENMAT, &
                  L_OMEGA_GREEK, L_AVERAGE_SECANT, &
                  L_BVEC, &
                  STATUS_SUB, MESSAGE, TRACE_1 )

                IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
                  write(CF,'(I2)')FOURIER_COMPONENT
                  TRACE_2 = 'Error from VLIDORT_L_QBEAM_SOLUTION, '// &
                          'Called in VLIDORT_L_FOURIER, Fourier # '//CF
                  STATUS = VLIDORT_SERIOUS
                  RETURN
                ENDIF

              ENDIF

!  B. the Green's function solution
!  ++++++++++++++++++++++++++++++++

            ELSE

!  $$$$$$$$$$$$$$$$$$$$$$    P L A C E H O L D E R    $$$$$$$$$$$$$$$$$$
!            CALL V_GREENFUNC_SOLUTION ( FOURIER)
!  $$$$$$$$$$$$$$$$$$$$$$    P L A C E H O L D E R    $$$$$$$$$$$$$$$$$$

            ENDIF

!  C. user solutions
!  +++++++++++++++++

            IF  ( STERM_LAYERMASK_UP(LAYER) .OR. &
                  STERM_LAYERMASK_DN(LAYER) ) THEN
              IF ( DO_USER_STREAMS ) THEN

!  (classical)
                IF ( DO_CLASSICAL_SOLUTION ) THEN

                  CALL VLIDORT_UBEAM_SOLUTION ( &
                    LAYER, FOURIER_COMPONENT, IBEAM, &
                    DO_UPWELLING, DO_DNWELLING, &
                    NSTOKES, NSTREAMS, FLUX_FACTOR, &
                    LAYER_PIS_CUTOFF, QUAD_HALFWTS, &
                    NMOMENTS, N_USER_STREAMS, &
                    DO_LAYER_SCATTERING, LOCAL_UM_START, &
                    STERM_LAYERMASK_UP, &
                    STERM_LAYERMASK_DN, DFLUX, &
                    OMEGA_GREEK, PI_XQP, PI_XUP, &
                    PI_XUM, PI_X0P, PI_XQM_PRE, &
                    BVEC, HELPSTOKES_BEAM, &
                    UPAR_DN_1, UPAR_DN_2, &
                    UPAR_UP_1, UPAR_UP_2 )

                  IF ( DO_ATMOS_LINEARIZATION ) THEN
                    CALL VLIDORT_L_UBEAM_SOLUTION ( &
                      LAYER, FOURIER_COMPONENT, IBEAM, &
                      DO_RTSOL_VARY, NPARAMS_VARY, &
                      DO_PLANE_PARALLEL, DO_UPWELLING, &
                      DO_DNWELLING, NSTOKES, &
                      NSTREAMS, FLUX_FACTOR, &
                      LAYER_PIS_CUTOFF, QUAD_HALFWTS, &
                      NMOMENTS, N_USER_STREAMS, &
                      DO_LAYER_SCATTERING, LOCAL_UM_START, &
                      STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, &
                      DFLUX, OMEGA_GREEK, &
                      PI_XQP, PI_XUP, PI_XUM, PI_X0P, &
                      PI_XQM_PRE, HELPSTOKES_BEAM, &
                      DO_PROFILE_LINEARIZATION, DO_COLUMN_LINEARIZATION, &
                      LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
                      L_BVEC, L_OMEGA_GREEK, &
                      L_UPAR_DN_1, L_UPAR_UP_1, &
                      L_UPAR_DN_2, L_UPAR_UP_2 )
                  ENDIF

!  (Green's function)

                ELSE
!  $$$$$$$$$$$$$$$$$$$$$$    P L A C E H O L D E R    $$$$$$$$$$$$$$$$$$
!            CALL V_GREENFUNC_SOLUTION ( FOURIER)
!  $$$$$$$$$$$$$$$$$$$$$$    P L A C E H O L D E R    $$$$$$$$$$$$$$$$$$
                ENDIF

              END IF
            END IF

!  end layer loop

          END DO

!  Add thermal solutions if flagged
!    4.PI modulus on the thermal contribution.

          IF ( DO_INCLUDE_THERMEMISS ) THEN
            O1 = 1
            DO N = 1, NLAYERS
              DO I = 1, NSTREAMS_2
! 1/26/10       WUPPER(I,O1,N) = WUPPER(I,O1,N) + pi4*T_WUPPER(I,N)
! 1/26/10       WLOWER(I,O1,N) = WLOWER(I,O1,N) + pi4*T_WLOWER(I,N)
                WUPPER(I,O1,N) = WUPPER(I,O1,N) + T_WUPPER(I,N)
                WLOWER(I,O1,N) = WLOWER(I,O1,N) + T_WLOWER(I,N)
              ENDDO
            ENDDO
          ENDIF

!  Adding the linearized thermal solutions is done later.....

!  Step 2. Solve boundary value problem and get Stokes vector
!  ----------------------------------------------------------

!  2A. Boundary Value problem
!  ==========================

!  standard case using compressed-band matrices, etc..

!       if ( do_write_screen) write(*,*)'bvp solution',ibeam
          IF ( BVP_REGULAR_FLAG(FOURIER_COMPONENT) ) THEN

!  Get the BVP solution (regular case) for this beam component

            CALL BVP_SOLUTION_MASTER ( &
              DO_INCLUDE_SURFACE, &
              DO_REFLECTED_DIRECTBEAM(IBEAM), &
              DO_INCLUDE_SURFEMISS, FOURIER_COMPONENT, &
              IBEAM, SURFACE_FACTOR, &
              NSTOKES, NSTREAMS, &
              NLAYERS, DO_LAMBERTIAN_SURFACE, &
              LAMBERTIAN_ALBEDO, &
              QUAD_STRMWTS, MUELLER_INDEX, &
              WLOWER, &
              EMISSIVITY, SURFBB, &
              NSTREAMS_2, NTOTAL, &
              NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
              DIRECT_BEAM, &
              WUPPER, &
              DO_DEBUG_WRITE, DO_FDTEST, &
              N_SUBDIAG, N_SUPDIAG, &
              K_REAL, K_COMPLEX, &
              BANDMAT2, IPIVOT, &
              SMAT2, SIPIVOT, &
              AXBID_F, &
              R2_BEAM, LCON, MCON, &
              STATUS_SUB, MESSAGE, TRACE_1 )

            IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
              write(CF,'(I2)')FOURIER_COMPONENT
              TRACE_2 = 'Error return from BVP_SOLUTION_MASTER, '// &
                'Called in VLIDORT_L_FOURIER, Fourier # '//CF
              STATUS = VLIDORT_SERIOUS
              RETURN
            ENDIF

!  Telescoped case

          ELSE

            CALL BVPTEL_SOLUTION_MASTER ( &
              DO_INCLUDE_SURFACE, &
              DO_REFLECTED_DIRECTBEAM(IBEAM), &
              FOURIER_COMPONENT, IBEAM, &
              SURFACE_FACTOR, &
              NLAYERS, DO_SPECIALIST_OPTION_2, &
              NLAYERS_TEL, ACTIVE_LAYERS, &
              NSTOKES, NSTREAMS, &
              DO_LAMBERTIAN_SURFACE, &
              LAMBERTIAN_ALBEDO, &
              QUAD_STRMWTS, MUELLER_INDEX, WLOWER, &
              NSTREAMS_2, NSTKS_NSTRMS, &
              NSTKS_NSTRMS_2, &
              DIRECT_BEAM, WUPPER, &
              T_DELT_DISORDS, T_DELT_EIGEN, &
              K_REAL, K_COMPLEX, &
              SOLA_XPOS, SOLB_XNEG, &
              IPIVOT, SMAT2, SIPIVOT, &
              N_BVTELMATRIX_SIZE, &
              N_BVTELMATRIX_SUPDIAG, N_BVTELMATRIX_SUBDIAG, &
              BANDTELMAT2, IPIVOTTEL, &
              AXBID_F, R2_BEAM, &
              LCON, MCON, &
              STATUS_SUB, MESSAGE, TRACE_1 )

            IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
              write(CF,'(I2)')FOURIER_COMPONENT
              TRACE_2 = 'Error return from BVPTEL_SOLUTION_MASTER, '// &
                'Called in VLIDORT_L_FOURIER, Fourier # '//CF
              STATUS = VLIDORT_SERIOUS
              RETURN
            ENDIF

          ENDIF

!  2B. Stokes vector Post Processing
!  =================================

!  upwelling

          IF ( DO_UPWELLING ) THEN

! LIDORT COMMENTS:  Direct beam inclusion flag:
!   This now has the DBCORRECTION option: if the DBCORRECTION Flag
!   is set, then we will be doing exact calculations of the reflected
!   directbeam, so we do not need to include it in the Post-processing.
!   However, the direct beam will need to be included in the basic RT
!   solution (the BVP), and this is controlled separately by the
!   DO_REFLECTED_DIRECTBEAM(IBEAM) flags.
!     R. Spurr, RT Solutions, Inc., 19 August 2005.

!Rob fix 4/12/12 - added DO_MSMODE_VLIDORT
            DO_INCLUDE_DIRECTBEAM = ( DO_UPWELLING .AND. &
              (DO_REFLECTED_DIRECTBEAM(IBEAM).AND..NOT.DO_DBCORRECTION) ) &
              .AND. .NOT.DO_MSMODE_VLIDORT

            CALL VLIDORT_UPUSER_INTENSITY ( &
              DO_INCLUDE_SURFACE, DO_INCLUDE_SURFEMISS, &
              DO_INCLUDE_THERMEMISS, DO_INCLUDE_DIRECTBEAM, &
              DO_INCLUDE_MVOUTPUT, FOURIER_COMPONENT, IBEAM, &
              FLUX_MULTIPLIER, SURFACE_FACTOR, DO_DEBUG_WRITE, DO_FDTEST, &
              NSTOKES, NLAYERS, N_USER_LEVELS, DO_TOA_CONTRIBS, &
              N_USER_STREAMS, DO_LAYER_SCATTERING, &
              DO_USER_STREAMS, LOCAL_UM_START, &
              PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX, &
              UTAU_LEVEL_MASK_UP, T_DELT_USERM, T_UTUP_USERM, &
              CUMTRANS, &
              DO_QUAD_OUTPUT, NSTREAMS, DO_LAMBERTIAN_SURFACE, &
              LAMBERTIAN_ALBEDO, USER_BRDF_F, EMISSIVITY, USER_EMISSIVITY, &
              SURFBB, DO_THERMAL_TRANSONLY, DO_DBCORRECTION, QUAD_WEIGHTS, &
              QUAD_STRMWTS, MUELLER_INDEX, K_REAL, K_COMPLEX, &
              SOLA_XPOS, SOLB_XNEG, WLOWER, LCON, MCON, &
              USER_DIRECT_BEAM, &
              T_DELT_DISORDS, T_DELT_EIGEN, T_WLOWER, &
              DO_SOLAR_SOURCES, DO_CLASSICAL_SOLUTION, &
              DO_MSMODE_VLIDORT, DO_MSMODE_THERMAL, &
              UHOM_UPDN, UHOM_UPUP, UPAR_UP_1, UPAR_UP_2, &
              HMULT_1, HMULT_2, EMULT_UP, LAYER_TSUP_UP, &
              UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP, LAYER_TSUP_UTUP, BRDF_F, &
              STOKES_DOWNSURF, &
              CUMSOURCE_UP, STOKES_F, MS_CONTRIBS_F, BOA_THTONLY_SOURCE )
          ENDIF

!  Downwelling

          IF ( DO_DNWELLING ) THEN
            CALL VLIDORT_DNUSER_INTENSITY ( &
              DO_INCLUDE_THERMEMISS, FOURIER_COMPONENT, &
              IBEAM, FLUX_MULTIPLIER, &
              DO_DEBUG_WRITE, NSTOKES, &
              N_USER_LEVELS, &
              N_USER_STREAMS, DO_LAYER_SCATTERING, &
              DO_USER_STREAMS, LOCAL_UM_START, &
              PARTLAYERS_OUTFLAG, &
              PARTLAYERS_OUTINDEX, UTAU_LEVEL_MASK_DN, &
              PARTLAYERS_LAYERIDX, &
              T_DELT_USERM, T_UTDN_USERM, &
              DO_SOLAR_SOURCES, DO_THERMAL_TRANSONLY, &
              DO_CLASSICAL_SOLUTION, DO_MSMODE_VLIDORT, &
              K_REAL, K_COMPLEX, &
              LCON, MCON, &
              UHOM_DNDN, UHOM_DNUP, &
              UPAR_DN_1, UPAR_DN_2, &
              HMULT_1, HMULT_2, &
              EMULT_DN, LAYER_TSUP_DN, &
              UT_HMULT_DU, UT_HMULT_DD, &
              UT_EMULT_DN, LAYER_TSUP_UTDN, &
              CUMSOURCE_DN, STOKES_F )
          ENDIF

!  mean value (integrated) output.
!    ------- Can also use this to get debug Quadrature output
! @@@ Rob fix 1/31/11, - added FLUX_FACTOR argument

          IF ( DO_INCLUDE_MVOUTPUT.OR.DO_QUAD_OUTPUT ) THEN
            CALL VLIDORT_INTEGRATED_OUTPUT ( &
              DO_INCLUDE_MVOUTPUT, &
              DO_REFLECTED_DIRECTBEAM(IBEAM), &
              DO_INCLUDE_THERMEMISS, &
              IBEAM, FLUX_MULTIPLIER, FLUX_FACTOR, &    !  added here
              NSTOKES, NSTREAMS, &
              N_USER_LEVELS, &
              FLUXVEC, LAYER_PIS_CUTOFF, &
              QUAD_WEIGHTS, QUAD_STRMWTS, &
              N_DIRECTIONS, WHICH_DIRECTIONS, &
              PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
              UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, &
              PARTLAYERS_LAYERIDX, &
              T_DELT_MUBAR, T_UTDN_MUBAR, &
              INITIAL_TRANS, LOCAL_CSZA, &
              DO_SOLAR_SOURCES, NLAYERS, &
              DO_THERMAL_TRANSONLY, &
              DO_CLASSICAL_SOLUTION, QUAD_STREAMS, &
              T_DELT_DISORDS, T_DISORDS_UTUP, &
              T_UTUP_EIGEN, T_UTDN_EIGEN, &
              K_REAL, K_COMPLEX, &
              SOLA_XPOS, SOLB_XNEG, &
              WUPPER, LCON, MCON, &
              T_WUPPER, UT_T_PARTIC, &
              BOA_THTONLY_SOURCE, &
              T_DELT_EIGEN, WLOWER, &
              T_DISORDS_UTDN, T_WLOWER, &
              MEAN_STOKES, FLUX_STOKES, &
              MEAN_DIRECT, FLUX_DIRECT )
          ENDIF

!  Finished this Beam solution, if only Stokes vector is required

          IF ( DO_SIMULATION_ONLY ) GO TO 4000

!  Step 3. Atmospheric Profile weighting function loop
!  ---------------------------------------------------

          IF ( DO_PROFILE_LINEARIZATION ) THEN

!  Start over layers that will contain variations
!    - Set flag for variation of layer, and number of variations

           DO LAYER_TO_VARY = 1, NLAYERS
            IF ( LAYER_VARY_FLAG(LAYER_TO_VARY) ) THEN
             N_LAYER_WFS = LAYER_VARY_NUMBER(LAYER_TO_VARY)

!  3A. Solve the linearized BVP
!  ============================

!  Get the linearized BVP solution for this beam component

!    (a) Regular case (full matrix linearization)

             IF ( BVP_REGULAR_FLAG(FOURIER_COMPONENT) ) THEN

              CALL L_BVP_SOLUTION_MASTER ( &
                DO_INCLUDE_SURFACE, &
                DO_REFLECTED_DIRECTBEAM(IBEAM), &
                DO_INCLUDE_THERMEMISS, &
                LAYER_TO_VARY, N_LAYER_WFS, &
                FOURIER_COMPONENT, IBEAM, &
                SURFACE_FACTOR, DO_PROFILE_LINEARIZATION, &
                DO_COLUMN_LINEARIZATION, &
                DO_SOLAR_SOURCES, NSTOKES, &
                NSTREAMS, NLAYERS, NSTREAMS_2, &
                NTOTAL, NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
                DELTAU_SLANT, T_DELT_EIGEN, &
                K_REAL, K_COMPLEX, &
                SOLA_XPOS, SOLB_XNEG, &
                LCON, MCON, DIRECT_BEAM, &
                L_DELTAU_VERT, L_T_DELT_EIGEN, &
                L_SOLA_XPOS, L_SOLB_XNEG, &
                L_T_WUPPER, L_T_WLOWER, &
                DO_CLASSICAL_SOLUTION, LAYER_PIS_CUTOFF, &
                DO_LAYER_SCATTERING, &
                T_DELT_MUBAR, INITIAL_TRANS, BVEC, &
                L_INITIAL_TRANS, L_T_DELT_MUBAR, L_BVEC, &
                DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
                BRDF_F, QUAD_STRMWTS, MUELLER_INDEX, &
                LAYER_TO_VARY, N_LAYER_WFS, &
                DO_PLANE_PARALLEL, &
                N_SUBDIAG, N_SUPDIAG, &
                BANDMAT2, IPIVOT, &
                SMAT2, SIPIVOT, &
                L_WLOWER, L_WUPPER, NCON, PCON, &
                STATUS_SUB, MESSAGE, TRACE_1 )

              IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
                write(CF,'(I2)')FOURIER_COMPONENT
                TRACE_2 = 'Error return from L_BVP_SOLUTION_MASTER, '// &
        'Profile Jacobians, Called in VLIDORT_L_FOURIER, Fourier # '//CF
                STATUS = VLIDORT_SERIOUS
                RETURN
              ENDIF

!    (b) telescoped case

            ELSE

             CALL L_BVPTEL_SOLUTION_MASTER ( &
               DO_INCLUDE_SURFACE, &
               DO_REFLECTED_DIRECTBEAM(IBEAM), &
               LAYER_TO_VARY, N_LAYER_WFS, &
               FOURIER_COMPONENT, IBEAM, SURFACE_FACTOR, &
               NSTOKES, NSTREAMS, NLAYERS, &
               NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
               K_REAL, K_COMPLEX, &
               SOLA_XPOS, SOLB_XNEG, &
               LCON, MCON, N_BVTELMATRIX_SIZE, &
               N_BVTELMATRIX_SUPDIAG, N_BVTELMATRIX_SUBDIAG, &
               NLAYERS_TEL, ACTIVE_LAYERS, &
               T_DELT_DISORDS, T_DELT_EIGEN, &
               DO_PROFILE_LINEARIZATION, DO_COLUMN_LINEARIZATION, &
               L_T_DELT_EIGEN, L_T_DELT_DISORDS, &
               L_SOLA_XPOS, L_SOLB_XNEG, &
               DO_SPECIALIST_OPTION_2, NSTREAMS_2, &
               DELTAU_SLANT, &
               L_DELTAU_VERT, DIRECT_BEAM, &
               DO_CLASSICAL_SOLUTION, LAYER_PIS_CUTOFF, &
               DO_LAYER_SCATTERING, DO_PLANE_PARALLEL, &
               T_DELT_MUBAR, INITIAL_TRANS, BVEC, &
               L_INITIAL_TRANS, L_T_DELT_MUBAR, L_BVEC, &
               DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
               BRDF_F, QUAD_STRMWTS, MUELLER_INDEX, &
               BANDTELMAT2, IPIVOTTEL, &
               SMAT2, SIPIVOT, &
               L_WLOWER, L_WUPPER, NCON, PCON, &
               STATUS_SUB, MESSAGE, TRACE_1 )

              IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
                write(CF,'(I2)')FOURIER_COMPONENT
                TRACE_2 = 'Error from L_BVPTEL_SOLUTION_MASTER, '// &
                'Profile Jacobians, Called in VLIDORT_L_FOURIER, Fourier # '//CF
                STATUS = VLIDORT_SERIOUS
                RETURN
              ENDIF

            ENDIF

!  3B. Post-processing for the weighting functions
!  ===============================================

!  Sequence:
!    Upwelling   Atmospheric weighting functions
!    Downwelling Atmospheric weighting functions
!    Mean-value  Atmospheric weighting functions

             IF ( DO_UPWELLING ) THEN
!            if ( do_write_screen) write(*,*)'wf up',ibeam,layer_to_vary
              CALL VLIDORT_UPUSER_ATMOSWF ( &
                DO_INCLUDE_SURFACE, DO_INCLUDE_THERMEMISS, &
                DO_INCLUDE_MVOUTPUT, DO_INCLUDE_DIRECTBEAM, &
                SURFACE_FACTOR, FLUX_MULTIPLIER, FOURIER_COMPONENT, &
                IBEAM, LAYER_TO_VARY, N_LAYER_WFS, &
                NSTOKES, NLAYERS, N_USER_LEVELS, &
                N_USER_STREAMS, DO_LAYER_SCATTERING, &
                DO_USER_STREAMS, LOCAL_UM_START, &
                PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
                UTAU_LEVEL_MASK_UP, PARTLAYERS_LAYERIDX, &
                T_DELT_USERM, T_UTUP_USERM, CUMSOURCE_UP, &
                L_T_DELT_USERM, L_T_UTUP_USERM, &
                DO_SOLAR_SOURCES, DO_QUAD_OUTPUT, &
                NSTREAMS, DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
                BRDF_F, USER_BRDF_F, DO_THERMAL_TRANSONLY, DO_DBCORRECTION, &
                QUAD_WEIGHTS, QUAD_STRMWTS, MUELLER_INDEX, &
                DELTAU_SLANT, T_DELT_DISORDS, T_DELT_EIGEN, &
                K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, &
                LCON, MCON, T_WLOWER, USER_DIRECT_BEAM, &
                L_DELTAU_VERT, L_T_DELT_EIGEN, L_T_DELT_DISORDS, &
                L_SOLA_XPOS, L_SOLB_XNEG, &
                L_WLOWER, NCON, PCON, L_T_WLOWER, &
                DO_MSMODE_VLIDORT, &
                UHOM_UPDN, UHOM_UPUP, UPAR_UP_1, UPAR_UP_2, &
                HMULT_1, HMULT_2, EMULT_UP, &
                L_UHOM_UPDN, L_UHOM_UPUP, L_UPAR_UP_1, L_UPAR_UP_2, &
                L_HMULT_1, L_HMULT_2, L_EMULT_UP, &
                L_LAYER_TSUP_UP, &
                UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP, &
                L_UT_HMULT_UU, L_UT_HMULT_UD, L_UT_EMULT_UP, &
                L_LAYER_TSUP_UTUP, &
                L_BOA_THTONLY_SOURCE, ATMOSWF_F )
             ENDIF

             IF ( DO_DNWELLING ) THEN
!           if ( do_write_screen) write(*,*)'wf dn',ibeam,layer_to_vary
              CALL VLIDORT_DNUSER_ATMOSWF ( &
                DO_INCLUDE_THERMEMISS, FLUX_MULTIPLIER, &
                FOURIER_COMPONENT, IBEAM, &
                LAYER_TO_VARY, N_LAYER_WFS, &
                NSTOKES, N_USER_LEVELS, &
                N_USER_STREAMS, DO_LAYER_SCATTERING, &
                DO_USER_STREAMS, LOCAL_UM_START, &
                PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
                UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX, &
                T_DELT_USERM, T_UTDN_USERM, CUMSOURCE_DN, &
                L_T_DELT_USERM, L_T_UTDN_USERM, &
                DO_SOLAR_SOURCES, DO_THERMAL_TRANSONLY, &
                DO_MSMODE_VLIDORT, &
                K_REAL, K_COMPLEX, &
                LCON, MCON, &
                UHOM_DNDN, UHOM_DNUP, &
                UPAR_DN_1, UPAR_DN_2, &
                HMULT_1, HMULT_2, EMULT_DN, &
                NCON, PCON, &
                L_UHOM_DNDN, L_UHOM_DNUP, &
                L_UPAR_DN_1, L_UPAR_DN_2, &
                L_HMULT_1, L_HMULT_2, L_EMULT_DN, &
                L_LAYER_TSUP_DN, &
                UT_HMULT_DU, UT_HMULT_DD, UT_EMULT_DN, &
                L_UT_HMULT_DU, L_UT_HMULT_DD, &
                L_UT_EMULT_DN, &
                L_LAYER_TSUP_UTDN, &
                ATMOSWF_F )
             ENDIF

! @@@ Rob fix 1/31/11, - added FLUX_FACTOR argument

             IF ( DO_INCLUDE_MVOUTPUT .OR. DO_QUAD_OUTPUT ) THEN
              CALL VLIDORT_L_INTEGRATED_OUTPUT ( &
                DO_INCLUDE_MVOUTPUT, &
                DO_REFLECTED_DIRECTBEAM(IBEAM), &
                DO_INCLUDE_THERMEMISS, FLUX_MULTIPLIER, FLUX_FACTOR, & ! added
                IBEAM, LAYER_TO_VARY, N_LAYER_WFS, &
                NSTOKES, NSTREAMS, N_USER_LEVELS, &
                FLUXVEC, LAYER_PIS_CUTOFF, &
                QUAD_WEIGHTS, QUAD_STRMWTS, &
                N_DIRECTIONS, WHICH_DIRECTIONS, &
                PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
                UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, &
                PARTLAYERS_LAYERIDX, &
                T_DELT_MUBAR, T_UTDN_MUBAR, &
                INITIAL_TRANS, LOCAL_CSZA, &
                L_INITIAL_TRANS, L_T_DELT_MUBAR, &
                L_T_UTDN_MUBAR, &
                DO_SOLAR_SOURCES, &
                NLAYERS, DO_THERMAL_TRANSONLY, &
                DO_CLASSICAL_SOLUTION, QUAD_STREAMS, &
                T_DELT_DISORDS, T_DISORDS_UTUP, &
                T_UTUP_EIGEN, T_UTDN_EIGEN, &
                K_REAL, K_COMPLEX, &
                SOLA_XPOS, SOLB_XNEG, &
                WUPPER, LCON, MCON, T_WUPPER, &
                BOA_THTONLY_SOURCE, &
                L_T_UTUP_EIGEN, L_T_UTDN_EIGEN, &
                L_T_DELT_DISORDS, L_T_DISORDS_UTUP, &
                L_SOLA_XPOS, L_SOLB_XNEG, &
                L_WUPPER, NCON, PCON, &
                L_T_WUPPER, L_UT_T_PARTIC, &
                L_BOA_THTONLY_SOURCE, &
                T_DELT_EIGEN, L_T_DELT_EIGEN, &
                L_WLOWER, &
                T_DISORDS_UTDN, L_T_DISORDS_UTDN, &
                L_T_WLOWER, T_WLOWER, &
                MINT_ATMOSWF, MINT_ATMOSWF_DIRECT, &
                FLUX_ATMOSWF, FLUX_ATMOSWF_DIRECT )
             ENDIF

!  Finish loop over layers with variation

            ENDIF
           ENDDO

!  End profile atmospheric weighting functions

          ENDIF

!  Step 4. Atmospheric Bulk weighting functions
!  --------------------------------------------

          IF ( DO_COLUMN_LINEARIZATION ) THEN

!   variation index = 0

            VARIATION_INDEX = 0

!  4A. Solve the linearized BVP
!  ============================

!  (a) Regular BVP linearization

            IF ( BVP_REGULAR_FLAG(FOURIER_COMPONENT) ) THEN

!  Get the linearized BVP solution for this beam component

              CALL L_BVP_SOLUTION_MASTER ( &
                DO_INCLUDE_SURFACE, &
                DO_REFLECTED_DIRECTBEAM(IBEAM), &
                DO_INCLUDE_THERMEMISS, &
                VARIATION_INDEX, N_TOTALCOLUMN_WFS, &
                FOURIER_COMPONENT, IBEAM, &
                SURFACE_FACTOR, DO_PROFILE_LINEARIZATION, &
                DO_COLUMN_LINEARIZATION, &
                DO_SOLAR_SOURCES, NSTOKES, &
                NSTREAMS, NLAYERS, NSTREAMS_2, &
                NTOTAL, NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
                DELTAU_SLANT, T_DELT_EIGEN, &
                K_REAL, K_COMPLEX, &
                SOLA_XPOS, SOLB_XNEG, &
                LCON, MCON, DIRECT_BEAM, &
                L_DELTAU_VERT, L_T_DELT_EIGEN, &
                L_SOLA_XPOS, L_SOLB_XNEG, &
                L_T_WUPPER, L_T_WLOWER, &
                DO_CLASSICAL_SOLUTION, LAYER_PIS_CUTOFF, &
                DO_LAYER_SCATTERING, &
                T_DELT_MUBAR, INITIAL_TRANS, BVEC, &
                L_INITIAL_TRANS, L_T_DELT_MUBAR, L_BVEC, &
                DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
                BRDF_F, QUAD_STRMWTS, MUELLER_INDEX, &
                LAYER_TO_VARY, N_LAYER_WFS, &
                DO_PLANE_PARALLEL, &
                N_SUBDIAG, N_SUPDIAG, &
                BANDMAT2, IPIVOT, &
                SMAT2, SIPIVOT, &
                L_WLOWER, L_WUPPER, NCON, PCON, &
                STATUS_SUB, MESSAGE, TRACE_1 )

              IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
                write(CF,'(I2)')FOURIER_COMPONENT
                TRACE_2 = 'Error return from L_BVP_SOLUTION_MASTER, '// &
                'Column Jacobians, Called in VLIDORT_L_FOURIER, Fourier # '//CF
                STATUS = VLIDORT_SERIOUS
                RETURN
              ENDIF

!  (b) telescoped BVP linearization

            ELSE

              CALL L_BVPTEL_SOLUTION_MASTER ( &
                DO_INCLUDE_SURFACE, &
                DO_REFLECTED_DIRECTBEAM(IBEAM), &
                VARIATION_INDEX, N_TOTALCOLUMN_WFS, &
                FOURIER_COMPONENT, IBEAM, SURFACE_FACTOR, &
                NSTOKES, NSTREAMS, NLAYERS, &
                NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
                K_REAL, K_COMPLEX, &
                SOLA_XPOS, SOLB_XNEG, &
                LCON, MCON, N_BVTELMATRIX_SIZE, &
                N_BVTELMATRIX_SUPDIAG, N_BVTELMATRIX_SUBDIAG, &
                NLAYERS_TEL, ACTIVE_LAYERS, &
                T_DELT_DISORDS, T_DELT_EIGEN, &
                DO_PROFILE_LINEARIZATION, DO_COLUMN_LINEARIZATION, &
                L_T_DELT_EIGEN, L_T_DELT_DISORDS, &
                L_SOLA_XPOS, L_SOLB_XNEG, &
                DO_SPECIALIST_OPTION_2, NSTREAMS_2, &
                DELTAU_SLANT, &
                L_DELTAU_VERT, DIRECT_BEAM, &
                DO_CLASSICAL_SOLUTION, LAYER_PIS_CUTOFF, &
                DO_LAYER_SCATTERING, DO_PLANE_PARALLEL, &
                T_DELT_MUBAR, INITIAL_TRANS, BVEC, &
                L_INITIAL_TRANS, L_T_DELT_MUBAR, L_BVEC, &
                DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
                BRDF_F, QUAD_STRMWTS, MUELLER_INDEX, &
                BANDTELMAT2, IPIVOTTEL, &
                SMAT2, SIPIVOT, &
                L_WLOWER, L_WUPPER, NCON, PCON, &
                STATUS_SUB, MESSAGE, TRACE_1 )

              IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
                write(CF,'(I2)')FOURIER_COMPONENT
                TRACE_2 = 'Error from L_BVPTEL_SOLUTION_MASTER, '// &
                'Column Jacobians, Called in VLIDORT_L_FOURIER, Fourier # '//CF
                STATUS = VLIDORT_SERIOUS
                RETURN
              ENDIF

            ENDIF

!  4B. Post-processing for the weighting functions
!  ===============================================

!  Sequence:
!    Upwelling   Atmospheric weighting functions
!    Downwelling Atmospheric weighting functions
!    Mean-value  Atmospheric weighting functions

            IF ( DO_UPWELLING ) THEN
!            if ( do_write_screen) write(*,*)'wf up',ibeam,layer_to_vary
              CALL VLIDORT_UPUSER_ATMOSWF ( &
                DO_INCLUDE_SURFACE, DO_INCLUDE_THERMEMISS, &
                DO_INCLUDE_MVOUTPUT, DO_INCLUDE_DIRECTBEAM, &
                SURFACE_FACTOR, FLUX_MULTIPLIER, FOURIER_COMPONENT, &
                IBEAM, VARIATION_INDEX, N_TOTALCOLUMN_WFS, &
                NSTOKES, NLAYERS, N_USER_LEVELS, &
                N_USER_STREAMS, DO_LAYER_SCATTERING, &
                DO_USER_STREAMS, LOCAL_UM_START, &
                PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
                UTAU_LEVEL_MASK_UP, PARTLAYERS_LAYERIDX, &
                T_DELT_USERM, T_UTUP_USERM, CUMSOURCE_UP, &
                L_T_DELT_USERM, L_T_UTUP_USERM, &
                DO_SOLAR_SOURCES, DO_QUAD_OUTPUT, &
                NSTREAMS, DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
                BRDF_F, USER_BRDF_F, DO_THERMAL_TRANSONLY, DO_DBCORRECTION, &
                QUAD_WEIGHTS, QUAD_STRMWTS, MUELLER_INDEX, &
                DELTAU_SLANT, T_DELT_DISORDS, T_DELT_EIGEN, &
                K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, &
                LCON, MCON, T_WLOWER, USER_DIRECT_BEAM, &
                L_DELTAU_VERT, L_T_DELT_EIGEN, L_T_DELT_DISORDS, &
                L_SOLA_XPOS, L_SOLB_XNEG, &
                L_WLOWER, NCON, PCON, L_T_WLOWER, &
                DO_MSMODE_VLIDORT, &
                UHOM_UPDN, UHOM_UPUP, UPAR_UP_1, UPAR_UP_2, &
                HMULT_1, HMULT_2, EMULT_UP, &
                L_UHOM_UPDN, L_UHOM_UPUP, L_UPAR_UP_1, L_UPAR_UP_2, &
                L_HMULT_1, L_HMULT_2, L_EMULT_UP, &
                L_LAYER_TSUP_UP, &
                UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP, &
                L_UT_HMULT_UU, L_UT_HMULT_UD, L_UT_EMULT_UP, &
                L_LAYER_TSUP_UTUP, &
                L_BOA_THTONLY_SOURCE, ATMOSWF_F )
            ENDIF

            IF ( DO_DNWELLING ) THEN
!           if ( do_write_screen) write(*,*)'wf dn',ibeam,layer_to_vary
              CALL VLIDORT_DNUSER_ATMOSWF ( &
                DO_INCLUDE_THERMEMISS, FLUX_MULTIPLIER, &
                FOURIER_COMPONENT, IBEAM, &
                VARIATION_INDEX, N_TOTALCOLUMN_WFS, &
                NSTOKES, N_USER_LEVELS, &
                N_USER_STREAMS, DO_LAYER_SCATTERING, &
                DO_USER_STREAMS, LOCAL_UM_START, &
                PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
                UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX, &
                T_DELT_USERM, T_UTDN_USERM, CUMSOURCE_DN, &
                L_T_DELT_USERM, L_T_UTDN_USERM, &
                DO_SOLAR_SOURCES, DO_THERMAL_TRANSONLY, &
                DO_MSMODE_VLIDORT, &
                K_REAL, K_COMPLEX, &
                LCON, MCON, &
                UHOM_DNDN, UHOM_DNUP, &
                UPAR_DN_1, UPAR_DN_2, &
                HMULT_1, HMULT_2, EMULT_DN, &
                NCON, PCON, &
                L_UHOM_DNDN, L_UHOM_DNUP, &
                L_UPAR_DN_1, L_UPAR_DN_2, &
                L_HMULT_1, L_HMULT_2, L_EMULT_DN, &
                L_LAYER_TSUP_DN, &
                UT_HMULT_DU, UT_HMULT_DD, UT_EMULT_DN, &
                L_UT_HMULT_DU, L_UT_HMULT_DD, &
                L_UT_EMULT_DN, &
                L_LAYER_TSUP_UTDN, &
                ATMOSWF_F )
            ENDIF

! @@@ Rob fix 1/31/11, - added FLUX_FACTOR argument

            IF ( DO_INCLUDE_MVOUTPUT .OR. DO_QUAD_OUTPUT ) THEN
              CALL VLIDORT_L_INTEGRATED_OUTPUT ( &
                DO_INCLUDE_MVOUTPUT, &
                DO_REFLECTED_DIRECTBEAM(IBEAM), &
                DO_INCLUDE_THERMEMISS, FLUX_MULTIPLIER, FLUX_FACTOR, & ! added
                IBEAM, VARIATION_INDEX, N_TOTALCOLUMN_WFS, &
                NSTOKES, NSTREAMS, N_USER_LEVELS, &
                FLUXVEC, LAYER_PIS_CUTOFF, &
                QUAD_WEIGHTS, QUAD_STRMWTS, &
                N_DIRECTIONS, WHICH_DIRECTIONS, &
                PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
                UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, &
                PARTLAYERS_LAYERIDX, &
                T_DELT_MUBAR, T_UTDN_MUBAR, &
                INITIAL_TRANS, LOCAL_CSZA, &
                L_INITIAL_TRANS, L_T_DELT_MUBAR, &
                L_T_UTDN_MUBAR, &
                DO_SOLAR_SOURCES, &
                NLAYERS, DO_THERMAL_TRANSONLY, &
                DO_CLASSICAL_SOLUTION, QUAD_STREAMS, &
                T_DELT_DISORDS, T_DISORDS_UTUP, &
                T_UTUP_EIGEN, T_UTDN_EIGEN, &
                K_REAL, K_COMPLEX, &
                SOLA_XPOS, SOLB_XNEG, &
                WUPPER, LCON, MCON, T_WUPPER, &
                BOA_THTONLY_SOURCE, &
                L_T_UTUP_EIGEN, L_T_UTDN_EIGEN, &
                L_T_DELT_DISORDS, L_T_DISORDS_UTUP, &
                L_SOLA_XPOS, L_SOLB_XNEG, &
                L_WUPPER, NCON, PCON, &
                L_T_WUPPER, L_UT_T_PARTIC, &
                L_BOA_THTONLY_SOURCE, &
                T_DELT_EIGEN, L_T_DELT_EIGEN, &
                L_WLOWER, &
                T_DISORDS_UTDN, L_T_DISORDS_UTDN, &
                L_T_WLOWER, T_WLOWER, &
                MINT_ATMOSWF, MINT_ATMOSWF_DIRECT, &
                FLUX_ATMOSWF, FLUX_ATMOSWF_DIRECT )
            ENDIF

!  End atmospheric column weighting functions

          ENDIF

!  Step 5. Surface Reflectance weighting functions
!  -----------------------------------------------

          IF ( DO_SURFACE_LINEARIZATION ) THEN

!  Only do this if there is surface reflectance !

            IF ( DO_INCLUDE_SURFACE ) THEN
              CALL VLIDORT_SURFACE_WFS ( &
                DO_REFLECTED_DIRECTBEAM(IBEAM), DO_INCLUDE_SURFEMISS, &
                DO_MSMODE_THERMAL, &
                DO_INCLUDE_MVOUTPUT, N_SURFACE_WFS, &
                FOURIER_COMPONENT, IBEAM, SURFACE_FACTOR, &
                FLUX_MULTIPLIER, DO_UPWELLING, DO_DNWELLING, &
                DO_QUAD_OUTPUT, NLAYERS, NTOTAL, &
                N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
                K_REAL, K_COMPLEX, &
                BANDMAT2, IPIVOT, SMAT2, SIPIVOT, &
                NSTOKES, NSTREAMS, &
                DO_LAMBERTIAN_SURFACE, SURFBB, LS_BRDF_F, &
                LS_BRDF_F_0, LS_EMISSIVITY, QUAD_STRMWTS, &
                MUELLER_INDEX, T_DELT_EIGEN, SOLA_XPOS, SOLB_XNEG, &
                WLOWER, LCON, MCON, ATMOS_ATTN, &
                R2_HOMP, R2_HOMM, R2_BEAM, DIRECT_BEAM, &
                LAMBERTIAN_ALBEDO, USER_BRDF_F, &
                DO_THERMAL_TRANSONLY, LS_USER_BRDF_F, LS_USER_BRDF_F_0, &
                LS_USER_EMISSIVITY, &
                DO_DBCORRECTION, FLUXVEC, &
                N_USER_STREAMS, DO_USER_STREAMS, LOCAL_UM_START, &
                STOKES_DOWNSURF, N_USER_LEVELS, &
                PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
                UTAU_LEVEL_MASK_UP, PARTLAYERS_LAYERIDX, &
                T_DELT_USERM, T_UTUP_USERM, &
                UHOM_UPDN, UHOM_UPUP, HMULT_1, HMULT_2, &
                UT_HMULT_UU, UT_HMULT_UD, &
                UTAU_LEVEL_MASK_DN, T_UTDN_USERM, &
                UHOM_DNDN, UHOM_DNUP, UT_HMULT_DU, UT_HMULT_DD, &
                QUAD_WEIGHTS, N_DIRECTIONS, WHICH_DIRECTIONS, &
                T_DELT_DISORDS, T_DISORDS_UTUP, &
                T_UTUP_EIGEN, T_UTDN_EIGEN, &
                SURFACEWF_F, MINT_SURFACEWF, FLUX_SURFACEWF, &
                STATUS_SUB, MESSAGE, TRACE_1 )

!  Exception handling

              IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
                write(CF,'(I2)')FOURIER_COMPONENT
                TRACE_2 = 'Error return from VLIDORT_SURFACE_WFS, '// &
                 'Beam Called in VLIDORT_L_FOURIER, Fourier # '//CF
                STATUS = VLIDORT_SERIOUS
                RETURN
              ENDIF

            ENDIF

!  end of surface weighting functions

          ENDIF

!  Continuation point for avoiding weighting functions

 4000     CONTINUE

!  End loop over beam solutions

        END IF
      END DO

!  debug write

      IF ( DO_DEBUG_WRITE ) THEN
        write(54,'(a,I3,a,100I3)') &
           'Fourier ',FOURIER_COMPONENT, &
         ' : # complex evalues by layer: ', &
             (K_COMPLEX(LAYER),LAYER=1,NLAYERS)
      ENDIF

!  ######
!  finish
!  ######

      RETURN
      END SUBROUTINE VLIDORT_L_FOURIER

      END MODULE vlidort_l_masters

