! ###############################################################
! #                                                             #
! #                    THE VECTOR LIDORT MODEL                  #
! #                                                             #
! #  (Vector LInearized Discrete Ordinate Radiative Transfer)   #
! #   -      --         -        -        -         -           #
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
! #            VLIDORT_MASTER (master)                          #
! #            VLIDORT_FOURIER (master)                         #
! #                                                             #
! ###############################################################


      MODULE vlidort_masters

      PRIVATE
      PUBLIC :: VLIDORT_MASTER

      CONTAINS

      SUBROUTINE VLIDORT_MASTER ( &
        VLIDORT_FixIn, &
        VLIDORT_ModIn, &
        VLIDORT_Sup,   &
        VLIDORT_Out )

      USE VLIDORT_PARS

      USE VLIDORT_Inputs_def
      USE VLIDORT_Sup_InOut_def
      USE VLIDORT_Outputs_def

      USE VLIDORT_INPUTS
      USE VLIDORT_MISCSETUPS_MODULE
      USE VLIDORT_CORRECTIONS
      USE VLIDORT_CONVERGE_MODULE
      USE VLIDORT_WRITEMODULES

      IMPLICIT NONE

!  VLIDORT input structures

      TYPE(VLIDORT_Fixed_Inputs), INTENT (IN)       :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs), INTENT (INOUT) :: VLIDORT_ModIn

!  VLIDORT supplements structure

      TYPE(VLIDORT_Sup_InOut), INTENT (INOUT)       :: VLIDORT_Sup

!  VLIDORT output structure

      TYPE(VLIDORT_Outputs), INTENT (OUT)           :: VLIDORT_Out

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
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES ,MAX_DIRECTIONS )

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

!  Error handling

      INTEGER ::             STATUS_INPUTCHECK
      INTEGER ::             NCHECKMESSAGES
      CHARACTER (LEN=120) :: CHECKMESSAGES(0:MAX_MESSAGES)
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
!                       in Structure VLIDORT_Sup%BRDF, so do not need
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

!  BRDF Supplement inputs

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

!  ==================================
!  END COPY INPUTS TO LOCAL VARIABLES
!  ==================================

!  VLIDORT input debug

      IF (DO_DEBUG_INPUT) THEN
        CALL VLIDORT_DEBUG_INPUT_MASTER()
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

      STOKES      = ZERO
      MEAN_STOKES = ZERO
      FLUX_STOKES = ZERO
      MEAN_DIRECT = ZERO
      FLUX_DIRECT = ZERO

!  New 15 March 2012
      IF ( .NOT. DO_SS_EXTERNAL ) THEN
         STOKES_SS = ZERO
         STOKES_DB = ZERO
      ENDIF

      FOURIER_SAVED = 0
      N_GEOMETRIES  = 0
      SZA_OFFSETS   = 0
      VZA_OFFSETS   = 0

!  set multiplier

      SS_FLUX_MULTIPLIER = FLUX_FACTOR / PI4

!  Check input
!  -----------

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
           USER_VZANGLES, SZANGLES, USER_RELAZMS, &
           USER_VZANGLES_ADJUST, SZANGLES_ADJUST, USER_RELAZMS_ADJUST, &
           FAIL, MESSAGE, TRACE_1 )

!  Update exception handling. October 2010 2p4RTC

      IF ( FAIL ) THEN
        TRACE_2 = ' Failure in multi_outgoing_adjustgeom'
        TRACE_3 = ' ** VLIDORT_MASTER '
        STATUS_CALCULATION = VLIDORT_SERIOUS
        VLIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION
        VLIDORT_Out%Status%TS_MESSAGE = MESSAGE
        VLIDORT_Out%Status%TS_TRACE_1 = TRACE_1
        VLIDORT_Out%Status%TS_TRACE_2 = TRACE_2
        VLIDORT_Out%Status%TS_TRACE_3 = TRACE_3
        RETURN
      ENDIF

!  Chapman function calculation of TAUTHICK_INPUT
!  ----------------------------------------------

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
            TRACE_2 = 'Direct call in VLIDORT_MASTER'
            TRACE_3 = ' ** VLIDORT_MASTER '
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

!  open file, call standard input write, close file

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
             ( DO_SURFACE_EMISSION .and. DO_THERMAL_EMISSION )

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
!    Thermal-only treatment is isotropic, no azimuth

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
        NSOURCES = N_SZANGLES
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

      DO IBEAM = 1, NSOURCES
        BEAM_TESTCONV  ( IBEAM )  = 0
        BEAM_ITERATION ( IBEAM ) = .TRUE.
        DO L = 0, MAXFOURIER
          DO_MULTIBEAM   ( IBEAM, L ) = .TRUE.
        ENDDO
      ENDDO

!  start loop

      DO WHILE ( LOCAL_ITERATION .AND. &
                 FOURIER_COMPONENT .LT. N_FOURIER_COMPONENTS )

!  Fourier counter

        FOURIER_COMPONENT = FOURIER_COMPONENT + 1

!  Local start of user-defined streams. Should always be 1.
!    No zenith tolerance now.

        LOCAL_UM_START = 1

!  azimuth cosine/sine factors, using adjust geometries.

        IF ( FOURIER_COMPONENT .GT. 0 ) THEN
          DFC = DBLE(FOURIER_COMPONENT)
          DO UA = 1, LOCAL_N_USERAZM
            DO IB = 1, NSOURCES
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

        CALL VLIDORT_FOURIER ( &
          FOURIER_COMPONENT, SS_FLUX_MULTIPLIER, &
          DO_SSCORR_NADIR, DO_SSCORR_OUTGOING, &
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
          DO_SURFACE_EMISSION, SURFBB, DO_THERMAL_TRANSONLY, &
          DO_SPECIALIST_OPTION_2, DO_SPECIALIST_OPTION_3, &
          DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, SLTERM_ISOTROPIC, &
          SLTERM_USERANGLES, SLTERM_F_0, USER_SLTERM_F_0, &
          NLAYERS_CUTOFF, DO_TOA_CONTRIBS, DO_DIRECT_BEAM, &
          DO_CLASSICAL_SOLUTION, DO_DBCORRECTION, FLUXVEC, &
          COS_SZANGLES, SIN_SZANGLES, SUN_SZA_COSINES, DO_MULTIBEAM, &
          QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS, QUAD_HALFWTS, &
          DO_MSMODE_VLIDORT, DO_MSMODE_THERMAL, NMOMENTS, NSTREAMS_2, NTOTAL, &
          N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2, NBEAMS, N_USER_VZANGLES, &
          MUELLER_INDEX, DMAT, &
          DO_REAL_EIGENSOLVER, BVP_REGULAR_FLAG, LAYER_MAXMOMENTS, &
          DO_LAYER_SCATTERING, DO_USER_VZANGLES, LOCAL_UM_START, N_DIRECTIONS, WHICH_DIRECTIONS, &
          USER_VZANGLES_ADJUST, SZANGLES_ADJUST, USER_RELAZMS_ADJUST, USER_STREAMS, USER_SECANTS, &
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, &
          DO_PARTLAYERS, N_PARTLAYERS, PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES, &
          N_ALLLAYERS_UP, N_ALLLAYERS_DN, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, &
          DFLUX, N_GEOMETRIES, VZA_OFFSETS, &
          TAUGRID_INPUT, CHAPMAN_FACTORS, &
          TAUGRID, OMEGA_TOTAL, GREEKMAT_TOTAL, &
          STOKES_SS, STOKES_DB, SS_CONTRIBS, MEAN_STOKES, FLUX_STOKES, &
          MEAN_DIRECT, FLUX_DIRECT, MS_CONTRIBS_F, STOKES_F, &
          DO_INCLUDE_SURFACE, DO_INCLUDE_SURFEMISS, DO_INCLUDE_THERMEMISS, &
          STATUS_SUB, MESSAGE, TRACE_1, TRACE_2 )

!  error handling

        IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
          TRACE_3 = ' Called by VLIDORT_MASTER '
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

!   -- only done for beams which are still not converged
!      This is controlled by flag DO_MULTIBEAM

!   -- new criterion, SS is added for Fourier = 0, as this means that
!      higher-order terms will be relatively smaller, which implies
!      faster convergence in some circumstances (generally not often).

        IBEAM_COUNT = 0
        DO IBEAM = 1, NSOURCES
          IF ( DO_MULTIBEAM ( IBEAM, FOURIER_COMPONENT ) ) THEN

!  Convergence and radiance summation

            CALL VLIDORT_CONVERGE ( &
              AZMFAC, FOURIER_COMPONENT, LOCAL_N_USERAZM, &
              DO_SSCORR_NADIR, DO_SSCORR_OUTGOING, &
              DO_SS_EXTERNAL, & ! New 15 March 2012
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
!  Write Standard Fourier output
!  Close file if iteration has finished

!  New comment:
!    If the SS correction is set, Fourier=0 will include SS field

!        IF ( DO_WRITE_FOURIER ) THEN
!          FUNIT = VLIDORT_FUNIT
!          IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
!            OPEN(FUNIT,FILE=FOURIER_WRITE_FILENAME,STATUS='UNKNOWN')
!          ENDIF
!          CALL VLIDORT_WRITEFOURIER ( FUNIT, FOURIER_COMPONENT )
!          IF ( .NOT.LOCAL_ITERATION ) CLOSE ( FUNIT )
!        ENDIF

!  end iteration loop

      ENDDO

!  restore no azimuth flag

      DO_NO_AZIMUTH = SAVE_DO_NO_AZIMUTH

!  Major result output
!  ===================

!  Standard output
!  ---------------

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

        CLOSE(SUNIT)
      ENDIF

!  ========================================================
!  BEGIN COPY LOCAL VARIABLES TO OUTPUTS (IN/OUT variables)
!  ========================================================

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

      VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS       = N_USER_RELAZMS
      VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS         = USER_RELAZMS
      VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT  = USER_VZANGLES
      VLIDORT_ModIn%MUserVal%TS_GEOMETRY_SPECHEIGHT  = GEOMETRY_SPECHEIGHT

      VLIDORT_ModIn%MChapman%TS_EARTH_RADIUS         = EARTH_RADIUS

      VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT    = OMEGA_TOTAL_INPUT

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
      VLIDORT_Out%Main%TS_VZA_OFFSETS   = VZA_OFFSETS

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

      END SUBROUTINE VLIDORT_MASTER

!

      SUBROUTINE VLIDORT_FOURIER ( &
        FOURIER_COMPONENT, SS_FLUX_MULTIPLIER, &
        DO_SSCORR_NADIR, DO_SSCORR_OUTGOING, &
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
        DO_SURFACE_EMISSION, SURFBB, DO_THERMAL_TRANSONLY, &
        DO_SPECIALIST_OPTION_2, DO_SPECIALIST_OPTION_3, &
        DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, SLTERM_ISOTROPIC, &
        SLTERM_USERANGLES, SLTERM_F_0, USER_SLTERM_F_0, &
        NLAYERS_CUTOFF, DO_TOA_CONTRIBS, DO_DIRECT_BEAM, &
        DO_CLASSICAL_SOLUTION, DO_DBCORRECTION, FLUXVEC, &
        COS_SZANGLES, SIN_SZANGLES, SUN_SZA_COSINES, DO_MULTIBEAM, &
        QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS, QUAD_HALFWTS, &
        DO_MSMODE_VLIDORT, DO_MSMODE_THERMAL, NMOMENTS, NSTREAMS_2, NTOTAL, &
        N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2, NBEAMS, N_USER_STREAMS, &
        MUELLER_INDEX, DMAT, &
        DO_REAL_EIGENSOLVER, BVP_REGULAR_FLAG, LAYER_MAXMOMENTS, &
        DO_LAYER_SCATTERING, DO_USER_STREAMS, LOCAL_UM_START, N_DIRECTIONS, WHICH_DIRECTIONS, &
        USER_VZANGLES_ADJUST, SZANGLES_ADJUST, USER_RELAZMS_ADJUST, USER_STREAMS, USER_SECANTS, &
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, &
        DO_PARTLAYERS, N_PARTLAYERS, PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES, &
        N_ALLLAYERS_UP, N_ALLLAYERS_DN, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, &
        DFLUX, N_GEOMETRIES, VZA_OFFSETS, &
        TAUGRID_INPUT, CHAPMAN_FACTORS, &
        TAUGRID, OMEGA_TOTAL, GREEKMAT_TOTAL, &
        STOKES_SS, STOKES_DB, SS_CONTRIBS, MEAN_STOKES, FLUX_STOKES, &
        MEAN_DIRECT, FLUX_DIRECT, MS_CONTRIBS_F, STOKES_F, &
        DO_INCLUDE_SURFACE, DO_INCLUDE_SURFEMISS, DO_INCLUDE_THERMEMISS, &
        STATUS, MESSAGE, TRACE_1, TRACE_2 )

!  Complete Fourier component calculation for the Standard Code.
!    Consolidated V2p4RTC: Lambertian or BRDF Surface

      USE VLIDORT_PARS

      USE VLIDORT_MISCSETUPS_MODULE
      USE VLIDORT_THERMALSUP
      USE VLIDORT_MULTIPLIERS
      USE VLIDORT_CORRECTIONS
      USE VLIDORT_SOLUTIONS
      USE VLIDORT_BVPROBLEM
      USE VLIDORT_THERMALSUP
      USE VLIDORT_INTENSITY

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
          ( 0:MAXMOMENTS,MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: USER_BRDF_F_0 &
          ( 0:MAXMOMENTS,MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS   )
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

      DOUBLE PRECISION, INTENT (INOUT) :: TAUGRID ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: OMEGA_TOTAL ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: GREEKMAT_TOTAL &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION, INTENT (INOUT) :: STOKES_SS &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) :: STOKES_DB &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) :: SS_CONTRIBS &
          ( MAX_GEOMETRIES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: MEAN_STOKES &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) :: FLUX_STOKES &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES ,MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) :: MEAN_DIRECT &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) :: FLUX_DIRECT &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES )

      DOUBLE PRECISION, INTENT (OUT) ::  MS_CONTRIBS_F &
          ( MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) ::  STOKES_F &
          ( MAX_USER_LEVELS, MAX_USER_VZANGLES, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      LOGICAL, INTENT (OUT) ::           DO_INCLUDE_SURFACE
      LOGICAL, INTENT (OUT) ::           DO_INCLUDE_SURFEMISS
      LOGICAL, INTENT (OUT) ::           DO_INCLUDE_THERMEMISS

      INTEGER, INTENT (OUT) ::           STATUS
      CHARACTER (LEN=*), INTENT (OUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (OUT) :: TRACE_1
      CHARACTER (LEN=*), INTENT (OUT) :: TRACE_2

!  Local variables
!  ---------------

      INTEGER ::           LAYER, IBEAM, IPARTIC, I, O1, N

!  local inclusion flags

      LOGICAL ::           DO_INCLUDE_MVOUTPUT
      LOGICAL ::           DO_INCLUDE_DIRECTBEAM

!  Flux multiplier and Fourier component numbers

      DOUBLE PRECISION ::  FLUX_MULTIPLIER
      DOUBLE PRECISION ::  DELTA_FACTOR
      DOUBLE PRECISION ::  SURFACE_FACTOR

!  error tracing

      LOGICAL ::           FAIL
      INTEGER ::           STATUS_SUB
      CHARACTER (LEN=2) :: CF

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
      DOUBLE PRECISION :: STOKES_DOWNSURF &
          ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION :: CUMSOURCE_UP &
          ( MAX_USER_STREAMS, MAXSTOKES, 0:MAXLAYERS )
      DOUBLE PRECISION :: CUMSOURCE_DN &
          ( MAX_USER_STREAMS, MAXSTOKES, 0:MAXLAYERS )

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

!  Albedo flag (for inclusion of some kind of reflecting boundary)

      DO_INCLUDE_SURFACE = .TRUE.
      IF ( DO_LAMBERTIAN_SURFACE ) THEN
        IF ( FOURIER_COMPONENT .NE. 0 ) THEN
            DO_INCLUDE_SURFACE = .FALSE.
        ELSE
          IF ( LAMBERTIAN_ALBEDO .EQ. ZERO ) THEN
            DO_INCLUDE_SURFACE = .FALSE.
          ENDIF
        ENDIF
      ENDIF

!  Direct beam flag (only if above albedo flag has been set)
!   (DO_REFLECTED_DIRECTBEAM is a stored Bookkeeping variable)

      IF ( DO_DIRECT_BEAM .and. DO_SOLAR_SOURCES ) THEN
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

!  Flux multiplier

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

!  Setups for Fourier = 0

!   MISCSETUPS    : Delta-M, average-secant formulation, transmittances
!   THERMAL_SETUP : Coefficients, direct multipliers
!   EMULT_MASTER  : Beam source function multipliers. Not required for the
!                   Full SS calculation in outgoing mode

!  30 January 2008
!  Telescoping setup is done in DERIVE_INPUTS (unlike LIDORT scalar code)

!!!!!!  WARNING CHECK ARGUMENTS on VLIDORT_MISCSETUPS

      IF ( FOURIER_COMPONENT .EQ. 0 ) THEN

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

!        CALL VLIDORT_MISCSETUPS (STATUS_SUB,MAIL)
!        IF ( STATUS_SUB .EQ. VLIDORT_WARNING ) THEN
!          TRACE= 'BVP Telescoping Warning from VLIDORT_MISCSETUPS'
!          STATUS = VLIDORT_WARNING
!          CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
!        ENDIF

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
        END IF

        IF ( DO_SOLAR_SOURCES ) THEN
          IF (.NOT.DO_SSFULL.OR.(DO_SSFULL.AND.DO_SSCORR_NADIR)) THEN
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

!  #####################
!  Correction operations
!  #####################

!  Single scatter correction (pre-calculation)
!  -------------------------------------------

!  Not required if no solar sources

      IF ( DO_SOLAR_SOURCES ) THEN

!  SSCORR_nadir:
!    Must be done after MISCSETUPS, as we need the multipliers and
!    transmittance factors for the SUN and LOS paths (SSCORR_NADIR)
!    Code added 6 May 2005. Replaces call in Master routine.

!  SSCORR_outgoing
!      Version 2.1. Added call to the new outgoing sphericity correction

       IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
        IF ( DO_USER_STREAMS ) THEN

!  regular (nadir view) SS correction

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
          ENDIF

!  New outgoing sphericity correction

          IF ( DO_SSCORR_OUTGOING ) THEN
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
              TRACE_2 = 'SS correction outgoing failed, VLIDORT_FOURIER'
              STATUS  = VLIDORT_SERIOUS
              RETURN
            ENDIF
          ENDIF

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
            ENDIF
          ENDIF

!  End Fourier = 0 and User stream

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

!      write(*,*)stokes_ss(1,1,1,1),stokes_db(1,1,1)

!  ###################
!  Spherical functions
!  ###################

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

!  ########################################
!  RT differential equation Eigensolutions
!  ########################################

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
          'Called in VLIDORT_FOURIER, Fourier component '//CF
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
            'Called in VLIDORT_FOURIER, Fourier # '//CF
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
            'Called in VLIDORT_FOURIER_MASTER, Fourier # '//CF
          STATUS = VLIDORT_SERIOUS
          RETURN
        ENDIF

      ENDIF

!  ################
!  Thermal Solution
!  ################

!  Continuation point for avoiding the scattering calculations

 8899 continue

!  1. Find the Particular solution (also for transmittance only)
!  2. Compute thermal layer source terms. (Upwelling and Downwelling)
!    These will be scaled up by factor 4.pi if solar beams as well
!  REMARK. No Green's function treatment here

!  @@@ Rob fix 1/31/11, This line was wrong (see below)
!          PI_XQP, PI_XQM, PI_XUP, PI_XQM_PRE, &

      IF ( DO_INCLUDE_THERMEMISS ) THEN

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
          PI_XQP, PI_XUP, PI_XUM, PI_XQM_PRE, &    ! Corrected 1/31/11
          SAB, DAB, &
          THERMCOEFFS, DELTAU_POWER, &
          XTAU_POWER, TCOM1, &
          T_WUPPER, T_WLOWER, &
          UT_T_PARTIC, U_TPOS1, &
          U_TNEG1, U_TPOS2, &
          U_TNEG2, &
          STATUS_SUB, MESSAGE, TRACE_1 )

        IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
          TRACE_2 = 'Called in VLIDORT_FOURIER, Fourier 0'
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

!  Skip the thermal-only section if there are solar sources

      IF ( DO_SOLAR_SOURCES ) GO TO 455

!  ####################################################
!  Complete Radiation Field with Thermal-only solutions
!  ####################################################

!  Only one solution, local direct_beam flag NOT set

      IPARTIC = 1
      DO_INCLUDE_DIRECTBEAM = .FALSE.

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
           'Called in VLIDORT_FOURIER, Fourier # '//CF
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

!  Solar beam Particular solutions (classical method)
!  --------------------------------------------------

          IF ( DO_CLASSICAL_SOLUTION ) THEN

            DO LAYER = 1, NLAYERS

!  get the solution

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
                  'Called in VLIDORT_FOURIER, Fourier # '//CF
                STATUS = VLIDORT_SERIOUS
                RETURN
              ENDIF

!  user solutions

              IF ( STERM_LAYERMASK_UP(LAYER) .OR. &
                   STERM_LAYERMASK_DN(LAYER) ) THEN
                IF ( DO_USER_STREAMS ) THEN
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
                END IF
              END IF

!  end layer loop

            END DO

!  Solar beam Particular solutions (Green's function)
!  -------------------------------------------------

          ELSE
!  $$$$$$$$$$$$$$$$$$$$$$    P L A C E H O L D E R    $$$$$$$$$$$$$$$$$$
!            CALL V_GREENFUNC_SOLUTION ( FOURIER)
!  $$$$$$$$$$$$$$$$$$$$$$    P L A C E H O L D E R    $$$$$$$$$$$$$$$$$$
          ENDIF

!  Add thermal solutions if flagged
!    NO modulus on the thermal contribution.
!   Bug fixed 26 January 2010

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

!  Solve boundary value problem
!  ----------------------------

!  standard case using compressed-band matrices, etc..

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
                'Called in VLIDORT_FOURIER, Fourier # '//CF
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
                'Called in VLIDORT_FOURIER, Fourier # '//CF
              STATUS = VLIDORT_SERIOUS
              RETURN
            ENDIF

          ENDIF

!  Stokes Vector Post Processing
!  -----------------------------

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
              (DO_REFLECTED_DIRECTBEAM(IBEAM) .AND. .NOT.DO_DBCORRECTION) ) &
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
! @@@ Rob fix 1/31/11, - added FLUX_FACTOR argument (3rd line)

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

!  End loop over beam solutions

        END IF
      END DO

!  Debug write

      IF ( DO_DEBUG_WRITE ) THEN
        write(*,'(a,I3,a,100I3)') &
           'Fourier ',FOURIER_COMPONENT, &
         ' : # complex evalues by layer: ', &
             (K_COMPLEX(LAYER),LAYER=1,NLAYERS)
      ENDIF

!  ######
!  finish
!  ######

      RETURN
      END SUBROUTINE VLIDORT_FOURIER


      END MODULE vlidort_masters

