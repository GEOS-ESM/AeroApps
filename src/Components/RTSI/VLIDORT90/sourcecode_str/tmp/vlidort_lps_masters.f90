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
! #                   2.5, 2.6, 2.7                             #
! #  Release Date :   December 2005  (2.0)                      #
! #  Release Date :   March 2007     (2.2)                      #
! #  Release Date :   October 2007   (2.3)                      #
! #  Release Date :   December 2008  (2.4)                      #
! #  Release Date :   April 2009     (2.4R)                     #
! #  Release Date :   July 2009      (2.4RT)                    #
! #  Release Date :   October 2010   (2.4RTC)                   #
! #  Release Date :   March 2011     (2.5)                      #
! #  Release Date :   May 2012       (2.6)                      #
! #  Release Date :   August 2014    (2.7)                      #
! #                                                             #
! #       NEW: TOTAL COLUMN JACOBIANS         (2.4)             #
! #       NEW: BPDF Land-surface KERNELS      (2.4R)            #
! #       NEW: Thermal Emission Treatment     (2.4RT)           #
! #       Consolidated BRDF treatment         (2.4RTC)          #
! #       f77/f90 Release                     (2.5)             #
! #       External SS / New I/O Structures    (2.6)             #
! #                                                             #
! #       SURFACE-LEAVING / BRDF-SCALING      (2.7)             #
! #       TAYLOR Series / OMP THREADSAFE      (2.7)             #
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
! #            VLIDORT_LPS_MASTER (master)                      #
! #            VLIDORT_LPS_FOURIER (master)                     #
! #                                                             #
! ###############################################################


      MODULE vlidort_lps_masters

      PRIVATE
      PUBLIC :: VLIDORT_LPS_MASTER

      CONTAINS

      SUBROUTINE VLIDORT_LPS_MASTER ( &
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

      USE VLIDORT_Work_def
      USE VLIDORT_LinWork_def

!  @@@ Rob Fix 5/10/13. Taylor series routines

      USE vlidort_Taylor_m

      USE VLIDORT_INPUTS
      USE VLIDORT_GEOMETRY
      USE VLIDORT_MISCSETUPS_MODULE
      USE VLIDORT_THERMALSUP
      USE VLIDORT_MULTIPLIERS
      USE VLIDORT_CORRECTIONS
      USE VLIDORT_INTENSITY, ONLY: VLIDORT_CONVERGE, &
                                   VLIDORT_CONVERGE_OBSGEO
      USE VLIDORT_WRITEMODULES

      USE VLIDORT_L_INPUTS
      USE VLIDORT_LA_MISCSETUPS
      USE VLIDORT_LP_MISCSETUPS
      USE VLIDORT_LA_CORRECTIONS
      USE VLIDORT_LP_CORRECTIONS
      USE VLIDORT_LS_CORRECTIONS
      USE VLIDORT_L_THERMALSUP
      USE VLIDORT_LS_WFSLEAVE

      USE vlidort_lbbf_jacobians_m    ! Need this, New 28 March 2014

      USE VLIDORT_LP_WFATMOS, ONLY: VLIDORT_LPS_CONVERGE, &
                                    VLIDORT_LPS_CONVERGE_OBSGEO
      USE VLIDORT_L_WRITEMODULES

      USE VLIDORT_PACK
      USE VLIDORT_L_PACK
      USE VLIDORT_LP_PACK

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

!  New 15 March 2012
      LOGICAL ::            DO_SS_EXTERNAL

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

!  New 17 May 2012
      LOGICAL ::            DO_SURFACE_LEAVING
      LOGICAL ::            DO_SL_ISOTROPIC

!  Order of Taylor series (including terms up to EPS^n).
!        Introduced 2/19/14 for Version 2.7
      
      INTEGER ::            TAYLOR_ORDER

      INTEGER ::            NSTOKES
      INTEGER ::            NSTREAMS
      INTEGER ::            NLAYERS
      INTEGER ::            NFINELAYERS
      INTEGER ::            N_THERMAL_COEFFS
      DOUBLE PRECISION ::   VLIDORT_ACCURACY
      INTEGER ::            NLAYERS_NOMS
      INTEGER ::            NLAYERS_CUTOFF

      DOUBLE PRECISION ::   FLUX_FACTOR

      INTEGER ::            N_USER_LEVELS

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

!  Dont need these anymore, superseded in Version 2.7 by New Code.
!      DOUBLE PRECISION ::   LTE_DELTAU_VERT_INPUT ( 2, MAXLAYERS )
!      DOUBLE PRECISION ::   LTE_THERMAL_BB_INPUT ( 0:MAXLAYERS )

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
      LOGICAL ::          DO_FO_CALC  !New 02 Jul 2013
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
      LOGICAL ::          DO_OBSERVATION_GEOMETRY

      INTEGER ::          NGREEK_MOMENTS_INPUT

      INTEGER ::          N_SZANGLES
      DOUBLE PRECISION :: SZANGLES ( MAX_SZANGLES )

      INTEGER ::          N_USER_RELAZMS
      DOUBLE PRECISION :: USER_RELAZMS  ( MAX_USER_RELAZMS )
      INTEGER ::          N_USER_VZANGLES
      DOUBLE PRECISION :: USER_VZANGLES ( MAX_USER_VZANGLES )
      DOUBLE PRECISION :: USER_LEVELS ( MAX_USER_LEVELS )
      DOUBLE PRECISION :: GEOMETRY_SPECHEIGHT
      INTEGER ::          N_USER_OBSGEOMS
      DOUBLE PRECISION :: USER_OBSGEOMS ( MAX_USER_OBSGEOMS, 3 )

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

      DOUBLE PRECISION ::   SLTERM_ISOTROPIC ( MAXSTOKES, MAXBEAMS )
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

      LOGICAL ::            LAYER_VARY_FLAG ( MAXLAYERS )
      INTEGER ::            LAYER_VARY_NUMBER ( MAXLAYERS )

      INTEGER ::            N_TOTALCOLUMN_WFS
      INTEGER ::            N_TOTALPROFILE_WFS
      INTEGER ::            N_SURFACE_WFS
      INTEGER ::            N_SLEAVE_WFS

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

      LOGICAL ::         DO_COLUMN_LINEARIZATION
      LOGICAL ::         DO_PROFILE_LINEARIZATION
      LOGICAL ::         DO_ATMOS_LINEARIZATION
      LOGICAL ::         DO_SURFACE_LINEARIZATION
      LOGICAL ::         DO_LINEARIZATION
      LOGICAL ::         DO_SLEAVE_WFS

!  Control for  Blackbody Jacobians, New 28 March 2014
!   Replaces the two commented out flags

      LOGICAL :: DO_ATMOS_LBBF, DO_SURFACE_LBBF
!      LOGICAL ::            DO_LTE_LINEARIZATION
!      LOGICAL ::            DO_SURFBB_LINEARIZATION

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

!  SS Column Jacobian Results at all angles and optical depths

      DOUBLE PRECISION :: COLUMNWF_SS &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, &
            MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )

!  DB Column Jacobian Results at all angles and optical depths

      DOUBLE PRECISION :: COLUMNWF_DB &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  SS Profile Jacobian Results at all angles and optical depths

      DOUBLE PRECISION :: PROFILEWF_SS &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
            MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )

!  DB Profile Jacobian Results at all angles and optical depths

      DOUBLE PRECISION :: PROFILEWF_DB &
          ( MAX_ATMOSWFS, MAXLAYERS, &
            MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  DB Surface Jacobian Results at all angles and optical depths

      DOUBLE PRECISION :: SURFACEWF_DB &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  New Surface-Leaving Inputs, 22 Aug 12
!  -------------------------------------

!@@@ Addition of SLEAVE WF stuff, R. Spurr, 22 August 2012 @@@@@@@@@

      DOUBLE PRECISION ::   LSSL_SLTERM_ISOTROPIC &
          ( MAX_SLEAVEWFS, MAXSTOKES, MAXBEAMS )
      DOUBLE PRECISION ::   LSSL_SLTERM_USERANGLES &
          ( MAX_SLEAVEWFS, MAXSTOKES, MAX_USER_STREAMS, &
            MAX_USER_RELAZMS, MAXBEAMS  )
      DOUBLE PRECISION ::   LSSL_SLTERM_F_0 &
          ( MAX_SLEAVEWFS, 0:MAXMOMENTS, MAXSTOKES, &
            MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION ::   LSSL_USER_SLTERM_F_0 &
          ( MAX_SLEAVEWFS, 0:MAXMOMENTS, MAXSTOKES, &
            MAX_USER_STREAMS, MAXBEAMS )

!  ------------------
!  Linearized Outputs
!  ------------------

! Fourier values

!      DOUBLE PRECISION :: PROFILEWF_F &
!         ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
!           MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!      DOUBLE PRECISION :: SURFACEWF_F &
!         ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_VZANGLES, &
!           MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  Fourier-summed values

      DOUBLE PRECISION :: PROFILEWF &
         ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
           MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )

      DOUBLE PRECISION :: SURFACEWF &
         ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
           MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )

!  Mean intensity (actinic flux) and regular flux weighting functions.
!  Total.

      DOUBLE PRECISION :: MINT_PROFILEWF &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

      DOUBLE PRECISION :: FLUX_PROFILEWF &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  Mean intensity (actinic flux) and regular flux weighting functions.
!  Direct beam only.

      DOUBLE PRECISION :: MINT_PROFILEWF_DIRECT &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES )

      DOUBLE PRECISION :: FLUX_PROFILEWF_DIRECT &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
           MAX_SZANGLES, MAXSTOKES )

!  Surface weighting functions

      DOUBLE PRECISION :: MINT_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

      DOUBLE PRECISION :: FLUX_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  BLACKBODY Jacobians, New 18 March 2014. V ersion 2.7
!  ====================================================

!  Postprocessed Jacobians.

      DOUBLE PRECISION :: ABBWFS_JACOBIANS &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, 0:MAXLAYERS, MAXSTOKES, MAX_DIRECTIONS)
      DOUBLE PRECISION :: SBBWFS_JACOBIANS &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAXSTOKES, MAX_DIRECTIONS)

!  Flux Jacobians.

      DOUBLE PRECISION :: ABBWFS_FLUXES &
          ( MAX_USER_LEVELS, 2, 0:MAXLAYERS, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: SBBWFS_FLUXES &
          ( MAX_USER_LEVELS, 2, MAXSTOKES, MAX_DIRECTIONS)

!  RT Solutions Inc. RJD Spurr.  11 September 2009.
!    LTE linearization: Introduction of T-Jacobians for BB functions
!     Only works with pure thermal emission (no scattering)
!      Introduced for the GEOCAPE study. Version 2.6 Superseded
!      DOUBLE PRECISION :: LTE_ATMOSWF &
!          ( 0:MAXLAYERS, MAX_USER_LEVELS, &
!            MAX_USER_VZANGLES, MAX_DIRECTIONS )

!  Error handling

      INTEGER ::             STATUS_INPUTCHECK
      INTEGER ::             NCHECKMESSAGES
      CHARACTER (LEN=120) :: CHECKMESSAGES (0:MAX_MESSAGES)
      CHARACTER (LEN=120) :: ACTIONS (0:MAX_MESSAGES)
      INTEGER ::             STATUS_CALCULATION
      CHARACTER (LEN=120) :: MESSAGE, TRACE_1, TRACE_2, TRACE_3

! ######################################################################

!                       Local Arguments
!                       +++++++++++++++

! ######################################################################

!  Local bookkeeping
!  ----------------

!               Intent(In) To the Fourier  routine
!               Intent(In) To the Converge routine

!  Mode of operation

      LOGICAL         :: DO_MSMODE_VLIDORT

! @@@@@@@@@@@@@ Robfix, 13 January 2012.
!               Add MSMODE_THERMAL flag
      logical         :: DO_MSMODE_THERMAL
! @@@@@@@@@@@@@ End Robfix, 13 January 2012.

!  Actual number of moments used in calculations
!   ( Normally 2 x NSTREAMS - 1 )

      INTEGER         :: NMOMENTS

!  NSTREAMS_2 = 2*NSTREAMS
!  total number of layers and streams NTOTAL = NSTREAMS_2 x NLAYERS
!  Number of super and sub diagonals in Band Matrix storage

      INTEGER         :: NSTREAMS_2
      INTEGER         :: NTOTAL
      INTEGER         :: N_SUBDIAG, N_SUPDIAG

!  Quadrature weights and abscissae, and product

      DOUBLE PRECISION :: QUAD_STREAMS (MAXSTREAMS)
      DOUBLE PRECISION :: QUAD_WEIGHTS (MAXSTREAMS)
      DOUBLE PRECISION :: QUAD_STRMWTS (MAXSTREAMS)

!  Angles/Cosines/sines of user-defined (off-quadrature) stream angles

      !DOUBLE PRECISION :: USER_ANGLES   (MAX_USER_STREAMS)
      DOUBLE PRECISION :: USER_STREAMS  (MAX_USER_STREAMS)
      DOUBLE PRECISION :: USER_SINES    (MAX_USER_STREAMS)
      DOUBLE PRECISION :: USER_SECANTS  (MAX_USER_STREAMS)

!  Output optical depth masks and indices

      LOGICAL         :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER         :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER         :: UTAU_LEVEL_MASK_UP  (MAX_USER_LEVELS)
      INTEGER         :: UTAU_LEVEL_MASK_DN  (MAX_USER_LEVELS)

!  Off-grid optical depths (values, masks, indices)

      INTEGER         :: N_PARTLAYERS
      INTEGER         :: PARTLAYERS_LAYERIDX     (MAX_PARTLAYERS)
      DOUBLE PRECISION :: PARTLAYERS_VALUES       (MAX_PARTLAYERS)

!  Layer masks for doing integrated source terms

      LOGICAL         :: STERM_LAYERMASK_UP(MAXLAYERS)
      LOGICAL         :: STERM_LAYERMASK_DN(MAXLAYERS)

!  Indexing numbers

      INTEGER         :: N_VIEWING

!  Offsets for geometry indexing

      !INTEGER         :: IBOFF(MAXBEAMS)
      !INTEGER         :: UMOFF(MAXBEAMS,MAX_USER_STREAMS)

!  Local input solar zenith angles Cosines
!  ( Only required for refractive geometry attenuation of the solar beam)

      DOUBLE PRECISION :: SUN_SZA_COSINES(MAXLAYERS,MAXBEAMS)

!  Local solar zenith angles Cosines (regular case)

      !DOUBLE PRECISION :: BEAM_COSINES(MAXBEAMS)

!  Solar beam flags (always internal)

      LOGICAL         :: DO_MULTIBEAM (MAXBEAMS,0:MAXFOURIER)

!  Number of directions (1 or 2) and directional array

      INTEGER         :: N_DIRECTIONS
      INTEGER         :: WHICH_DIRECTIONS ( MAX_DIRECTIONS )

!  Number of convergence tests

      INTEGER         :: N_CONVTESTS

!  Adjusted geometries. New, 2007.
!  -------------------------------

!               Intent(Out) from the Adjust-geometry routine
!               Intent(In)  to   the Correction      routine

      DOUBLE PRECISION :: USER_VZANGLES_ADJUST ( MAX_USER_VZANGLES )
      DOUBLE PRECISION :: SZANGLES_ADJUST &
          ( MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS )
      DOUBLE PRECISION :: USER_RELAZMS_ADJUST &
          ( MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS )

!  Arrays for setups and Corrections
!  ---------------------------------

!               Intent(In) To the Fourier routine

!  Local flags for the solution saving option

      INTEGER         :: LAYER_MAXMOMENTS (MAXLAYERS)

!  Initial transmittances * (secants)

      DOUBLE PRECISION :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Saved arrays for truncation factor and Delta-M scaling

      DOUBLE PRECISION :: TRUNC_FACTOR(MAXLAYERS)
      DOUBLE PRECISION :: FAC1(MAXLAYERS)

!  Derived Slant optical thickness inputs

      !DOUBLE PRECISION :: TAUSLANT    ( 0:MAXLAYERS, MAXBEAMS )
      !DOUBLE PRECISION :: TAUGRID     ( 0:MAXLAYERS )
      DOUBLE PRECISION :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      !DOUBLE PRECISION :: DELTAU_SLANT_UNSCALED( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Scaled SSAs and phase function moments

      DOUBLE PRECISION :: OMEGA_TOTAL    ( MAXLAYERS )
      DOUBLE PRECISION :: GREEKMAT_TOTAL &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )

!  Linearized Input optical depths after delta-M scaling

      DOUBLE PRECISION :: L_OMEGA_TOTAL    ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION :: L_GREEKMAT_TOTAL &
          ( MAX_ATMOSWFS, 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )

!  Derived Slant optical thickness inputs

      DOUBLE PRECISION :: L_DELTAU_SLANT (MAX_ATMOSWFS,MAXLAYERS,MAXLAYERS,MAXBEAMS)

!  Linearized truncation factor

      DOUBLE PRECISION :: L_TRUNC_FACTOR(MAX_ATMOSWFS,MAXLAYERS)

!  L'Hopital's rule logical variables

      LOGICAL          :: EMULT_HOPRULE (MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  Coefficient functions for user-defined angles

      DOUBLE PRECISION :: SIGMA_M(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)
      DOUBLE PRECISION :: SIGMA_P(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  Fourier component output
!  ------------------------

!               Intent(Out) from the Fourier routine
!               Intent(in)  to   the Converge routine

!  Fourier comonents User-defined solutions

      DOUBLE PRECISION :: STOKES_F &
          ( MAX_USER_LEVELS, MAX_USER_VZANGLES, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  Fourier-component Profile weighting functions at user angles

      DOUBLE PRECISION :: PROFILEWF_F &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
            MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  Fourier-component surface weighting functions at user angles

      DOUBLE PRECISION :: SURFACEWF_F &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_VZANGLES, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  Help arrays from the SS/DB correction routines
!  ==============================================

!  Saved Legendre polynomials

      DOUBLE PRECISION :: SS_PLEG_UP(MAX_GEOMETRIES,MAXLAYERS,0:MAXMOMENTS_INPUT)
      DOUBLE PRECISION :: SS_PLEG_DN(MAX_GEOMETRIES,MAXLAYERS,0:MAXMOMENTS_INPUT)

!  Saved TMS (Nakajima-Tanaka) factor

      DOUBLE PRECISION :: TMS(MAXLAYERS)

!  Local truncation factors for additional DELTAM scaling

      DOUBLE PRECISION :: SSFDEL ( MAXLAYERS )

!  Exact Phase function calculations

      !DOUBLE PRECISION :: EXACTSCAT_UP(MAX_GEOMETRIES,MAXLAYERS)
      !DOUBLE PRECISION :: EXACTSCAT_DN(MAX_GEOMETRIES,MAXLAYERS)

!  Cumulative single scatter source terms

      DOUBLE PRECISION :: SS_CUMSOURCE_UP(MAX_GEOMETRIES,MAXSTOKES,0:MAXLAYERS)
      DOUBLE PRECISION :: SS_CUMSOURCE_DN(MAX_GEOMETRIES,MAXSTOKES,0:MAXLAYERS)

!  Atmospheric attenuation before reflection

      DOUBLE PRECISION :: ATTN_DB_SAVE(MAX_GEOMETRIES)

!  Exact direct beam source terms

      DOUBLE PRECISION :: EXACTDB_SOURCE(MAX_GEOMETRIES, MAXSTOKES)

!  Cumulative direct bounce source terms

      DOUBLE PRECISION :: DB_CUMSOURCE(MAX_GEOMETRIES,MAXSTOKES,0:MAXLAYERS)

!  Solar beam attenuation to BOA (required for exact DB calculation)

      DOUBLE PRECISION :: BOA_ATTN(MAX_GEOMETRIES)

!  Outgoing sphericity stuff
!  Whole and part-layer LOS transmittance factors

      DOUBLE PRECISION :: UP_LOSTRANS(MAXLAYERS,MAX_GEOMETRIES)
      DOUBLE PRECISION :: UP_LOSTRANS_UT(MAX_PARTLAYERS,MAX_GEOMETRIES)

!  Outgoing sphericity stuff
!    Linearized output. Output required for later on.

      DOUBLE PRECISION :: L_UP_LOSTRANS   (MAXLAYERS,     MAX_ATMOSWFS,MAX_GEOMETRIES)
      DOUBLE PRECISION :: L_UP_LOSTRANS_UT(MAX_PARTLAYERS,MAX_ATMOSWFS,MAX_GEOMETRIES)

!  Solar beam attenuation to BOA (required for exact DB calculation)

      DOUBLE PRECISION :: LP_BOA_ATTN(MAXLAYERS,MAX_ATMOSWFS,MAX_GEOMETRIES)

!  Arrays required at the Top level
!  ================================

!               Intent(In) To the Fourier routine

!  Input optical properties after delta-M scaling

      DOUBLE PRECISION :: DELTAU_VERT    ( MAXLAYERS )
      DOUBLE PRECISION :: PARTAU_VERT    ( MAX_PARTLAYERS )
      DOUBLE PRECISION :: OMEGA_GREEK &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )

!  Linearized input optical properties after delta-M scaling

      LOGICAL          :: DO_SCATMAT_VARIATION ( MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_DELTAU_VERT  ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION :: L_OMEGA_GREEK &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )

!  Local input solar zenith angles by levels
!  ( Only required for refractive geometry attenuation of the solar beam)
!  These will be set internally if the refraction flag is set.

      DOUBLE PRECISION :: SZA_LOCAL_INPUT ( 0:MAXLAYERS, MAX_SZANGLES )

!  Last layer to include Particular integral solution

      INTEGER         :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Average-secant and initial tramsittance factors for solar beams.

      DOUBLE PRECISION :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: LOCAL_CSZA     ( 0:MAXLAYERS, MAXBEAMS )

!  Solar beam attenuation

      DOUBLE PRECISION :: TRANS_SOLAR_BEAM ( MAXBEAMS )

!  Rob fix 11/17/14. Proxy for new output.
!     NOTE. SOLARBEAM_BOATRANS is computed with Unscaled optical depths, TRANS_SOLAR_BEAM with scaled ODs.
!           [ They are the same if there is no Deltam-scaling]

      DOUBLE PRECISION :: SOLARBEAM_BOATRANS ( MAXBEAMS )

!  Linearized Average-secant and initial tramsittance factors for solar beams.

      DOUBLE PRECISION :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )
      DOUBLE PRECISION :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )

!  Reflectance flags

      LOGICAL         :: DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )

!  Local flags for the solution saving option

      LOGICAL         :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Local flags,  BVP telescoping enhancement

      LOGICAL         :: BVP_REGULAR_FLAG (0:MAXMOMENTS)

!  Masking for regular case. Required again for linearization

      !? INTEGER         :: LCONMASK(MAXSTREAMS,MAXLAYERS)
      !? INTEGER         :: MCONMASK(MAXSTREAMS,MAXLAYERS)

!  Telescoping initial flag (modified argument), Layer bookkeeping
!  Number of telescoped layers, active layers,  Size of BVP matrix 

      LOGICAL         :: DO_BVTEL_INITIAL
      INTEGER         :: BVTEL_FOURIER_COMPONENT
      !INTEGER         :: NLAYERS_TEL
      !INTEGER         :: ACTIVE_LAYERS ( MAXLAYERS )
      !INTEGER         :: N_BVTELMATRIX_SIZE

!  Set up for band matrix compression

      !? INTEGER         :: BMAT_ROWMASK(MAXTOTAL,MAXTOTAL)
      !? INTEGER         :: BTELMAT_ROWMASK(MAXTOTAL,MAXTOTAL)

!  Transmittance Setups
!  --------------------

!               Intent(In) To the Fourier routine

!  Discrete ordinate factors (BVP telescoping, solutions saving)
!  Code added by R. Spurr, RT SOLUTIONS Inc., 30 August 2005.

      DOUBLE PRECISION :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      DOUBLE PRECISION :: T_DISORDS_UTUP(MAXSTREAMS,MAX_PARTLAYERS)
      DOUBLE PRECISION :: T_DISORDS_UTDN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage

      DOUBLE PRECISION :: T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      DOUBLE PRECISION :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: T_UTUP_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage 

      DOUBLE PRECISION :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      DOUBLE PRECISION :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Cumulative transmission

      DOUBLE PRECISION :: CUMTRANS ( MAXLAYERS, MAX_USER_STREAMS )

!  Multiplier arrays
!  -----------------

!               Intent(In) To the Fourier routine

!  Forcing term multipliers (saved for whole atmosphere)

      DOUBLE PRECISION :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)
      DOUBLE PRECISION :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Partial layer multipliers

      DOUBLE PRECISION :: UT_EMULT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)
      DOUBLE PRECISION :: UT_EMULT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)

!  LINEARIZED Arrays required at the Top level
!  ===========================================

!  Linearized Transmittance Setups
!  -------------------------------

!  Linearized discrete ordinate factors (BVP telescoping, solutions saving)
!  Code added by R. Spurr, RT SOLUTIONS Inc., 30 August 2005.

      DOUBLE PRECISION :: L_T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS,      MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_DISORDS_UTUP ( MAXSTREAMS, MAX_USER_LEVELS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_DISORDS_UTDN ( MAXSTREAMS, MAX_USER_LEVELS, MAX_ATMOSWFS )

!  Linearized Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      DOUBLE PRECISION :: L_T_DELT_USERM ( MAXLAYERS,       MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_UTDN_USERM ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_UTUP_USERM ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS )

!  Linearized transmittances, solar beam

      DOUBLE PRECISION :: LP_T_DELT_MUBAR ( MAXLAYERS,      MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LP_T_UTDN_MUBAR ( MAX_USER_LEVELS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Beam multipliers
!  ---------------------------

!  Linearized whole layer multipliers

      DOUBLE PRECISION :: LP_EMULT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LP_EMULT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized part layer multipliers

      DOUBLE PRECISION :: LP_UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LP_UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Thermal Setup outputs
!  ---------------------

!  Optical depth powers

      DOUBLE PRECISION :: DELTAU_POWER (MAXLAYERS,     MAX_THERMAL_COEFFS)
      DOUBLE PRECISION :: XTAU_POWER   (MAX_PARTLAYERS,MAX_THERMAL_COEFFS)

!

      DOUBLE PRECISION :: TCOM1 ( MAXLAYERS, MAX_THERMAL_COEFFS )

!  Thermal coefficients

      DOUBLE PRECISION :: THERMCOEFFS (MAXLAYERS,MAX_THERMAL_COEFFS)

!  Tranmsittance solutions

      DOUBLE PRECISION :: T_DIRECT_UP ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: T_DIRECT_DN ( MAX_USER_STREAMS, MAXLAYERS )

      DOUBLE PRECISION :: T_UT_DIRECT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS)
      DOUBLE PRECISION :: T_UT_DIRECT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Linearized optical depth powers

      DOUBLE PRECISION :: L_DELTAU_POWER &
         (MAXLAYERS,     MAX_THERMAL_COEFFS,MAX_ATMOSWFS)
      DOUBLE PRECISION :: L_XTAU_POWER &
         (MAX_PARTLAYERS,MAX_THERMAL_COEFFS,MAX_ATMOSWFS)

!

      DOUBLE PRECISION :: L_TCOM1 &
          ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )

!  Linearized Thermal coefficients

      DOUBLE PRECISION :: L_THERMCOEFFS &
         (MAXLAYERS,MAX_THERMAL_COEFFS,MAX_ATMOSWFS)

!  Linearized Tranmsittance solutions

      DOUBLE PRECISION :: L_T_DIRECT_UP &
         ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_DIRECT_DN  &
         ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION :: L_T_UT_DIRECT_UP &
         ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_UT_DIRECT_DN &
         ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

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

      INTEGER ::          UM, UA, TESTCONV, L, LUM, LUA
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
      !DOUBLE PRECISION :: SUN_SZA_COSINES ( MAXLAYERS, MAX_SZANGLES )
      !LOGICAL ::          DO_MULTIBEAM ( MAXBEAMS, 0:MAXFOURIER )
      !DOUBLE PRECISION :: QUAD_STREAMS ( MAXSTREAMS )
      !DOUBLE PRECISION :: QUAD_WEIGHTS ( MAXSTREAMS )
      !DOUBLE PRECISION :: QUAD_STRMWTS ( MAXSTREAMS )
      DOUBLE PRECISION :: QUAD_HALFWTS ( MAXSTREAMS )
      DOUBLE PRECISION :: QUAD_SINES   ( MAXSTREAMS )
      DOUBLE PRECISION :: QUAD_ANGLES  ( MAXSTREAMS )
      !LOGICAL ::          DO_MSMODE_VLIDORT
      !LOGICAL ::          DO_MSMODE_THERMAL
      LOGICAL ::          DO_NO_AZIMUTH
      !INTEGER ::          NMOMENTS
      !INTEGER ::          NSTREAMS_2
      !INTEGER ::          NTOTAL
      !INTEGER ::          N_SUBDIAG
      !INTEGER ::          N_SUPDIAG
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
      !LOGICAL ::          BVP_REGULAR_FLAG ( 0:MAXMOMENTS )
      !INTEGER ::          LAYER_MAXMOMENTS ( MAXLAYERS )
      !LOGICAL ::          DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )
      !INTEGER ::          N_CONVTESTS
      INTEGER ::          N_CONV_STREAMS
      !LOGICAL ::          DO_USER_STREAMS
      INTEGER ::          LOCAL_UM_START
      !INTEGER ::          N_DIRECTIONS
      !INTEGER ::          WHICH_DIRECTIONS ( MAX_DIRECTIONS )
      !DOUBLE PRECISION :: USER_VZANGLES_ADJUST ( MAX_USER_VZANGLES )
      !DOUBLE PRECISION :: SZANGLES_ADJUST &
      !    ( MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS )
      !DOUBLE PRECISION :: USER_RELAZMS_ADJUST &
      !    ( MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS )
      !DOUBLE PRECISION :: USER_STREAMS ( MAX_USER_STREAMS )
      !DOUBLE PRECISION :: USER_SINES   ( MAX_USER_STREAMS )
      !DOUBLE PRECISION :: USER_SECANTS ( MAX_USER_STREAMS )
      INTEGER ::          N_OUT_STREAMS
      DOUBLE PRECISION :: OUT_ANGLES ( MAX_USER_STREAMS )
      !LOGICAL ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      !INTEGER ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      !INTEGER ::          UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      !INTEGER ::          UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )
      LOGICAL ::          DO_PARTLAYERS
      !INTEGER ::          N_PARTLAYERS
      !INTEGER ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      !DOUBLE PRECISION :: PARTLAYERS_VALUES   ( MAX_PARTLAYERS )
      INTEGER ::          N_LAYERSOURCE_UP
      INTEGER ::          N_LAYERSOURCE_DN
      INTEGER ::          N_ALLLAYERS_UP
      INTEGER ::          N_ALLLAYERS_DN
      !LOGICAL ::          STERM_LAYERMASK_UP ( MAXLAYERS )
      !LOGICAL ::          STERM_LAYERMASK_DN ( MAXLAYERS )
      !INTEGER ::          N_VIEWING
      DOUBLE PRECISION :: CHAPMAN_FACTORS &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: DFLUX ( MAXSTOKES )
      DOUBLE PRECISION :: TAUGRID_INPUT ( 0:MAXLAYERS )
      DOUBLE PRECISION :: TAUGRID ( 0:MAXLAYERS )
      !DOUBLE PRECISION :: OMEGA_TOTAL ( MAXLAYERS )
      !DOUBLE PRECISION :: GREEKMAT_TOTAL &
      !    ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )
      !DOUBLE PRECISION :: SZA_LOCAL_INPUT &
      !    ( 0:MAXLAYERS, MAX_SZANGLES )

      INTEGER ::          N_TOTALSURFACE_WFS

      !DOUBLE PRECISION :: STOKES_F &
      !    ( MAX_USER_LEVELS, MAX_USER_VZANGLES, &
      !      MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

      DOUBLE PRECISION :: MS_CONTRIBS_F &
          ( MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAXLAYERS  )
      DOUBLE PRECISION :: CONTRIBS &
          ( MAX_GEOMETRIES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION :: SS_CONTRIBS &
          ( MAX_GEOMETRIES, MAXSTOKES, MAXLAYERS )

      !DOUBLE PRECISION :: PROFILEWF_F &
      !    ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
      !      MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      !DOUBLE PRECISION :: SURFACEWF_F &
      !    ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_VZANGLES, &
      !      MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!FROM VLIDORT_SSCORR_NADIR:
      DOUBLE PRECISION :: ZMAT_UP &
          ( MAX_GEOMETRIES, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: ZMAT_DN &
          ( MAX_GEOMETRIES, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      !DOUBLE PRECISION :: TMS ( MAXLAYERS )
      !DOUBLE PRECISION :: SSFDEL ( MAXLAYERS )
      !DOUBLE PRECISION :: SS_CUMSOURCE_UP &
      !    ( MAX_GEOMETRIES, MAXSTOKES, 0:MAXLAYERS )
      !DOUBLE PRECISION :: SS_CUMSOURCE_DN &
      !    ( MAX_GEOMETRIES, MAXSTOKES, 0:MAXLAYERS )

!FROM VLIDORT_SSCORR_OUTGOING:
      DOUBLE PRECISION :: UP_MULTIPLIERS &
          ( MAXLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION :: DN_MULTIPLIERS &
          ( MAXLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION :: UP_MULTIPLIERS_UT &
          ( MAX_PARTLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION :: DN_MULTIPLIERS_UT &
          ( MAX_PARTLAYERS, MAX_GEOMETRIES )
      !DOUBLE PRECISION :: UP_LOSTRANS &
      !    ( MAXLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION :: DN_LOSTRANS &
          ( MAXLAYERS, MAX_GEOMETRIES )
      !DOUBLE PRECISION :: UP_LOSTRANS_UT &
      !    ( MAX_PARTLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION :: DN_LOSTRANS_UT &
          ( MAX_PARTLAYERS, MAX_GEOMETRIES )
      !DOUBLE PRECISION :: BOA_ATTN ( MAX_GEOMETRIES )

!FROM VLIDORT_DBCORRECTION:
      !DOUBLE PRECISION :: DB_CUMSOURCE &
      !    ( MAX_GEOMETRIES, MAXSTOKES, 0:MAXLAYERS )
      !DOUBLE PRECISION :: EXACTDB_SOURCE &
      !    ( MAX_GEOMETRIES, MAXSTOKES )
      !DOUBLE PRECISION :: ATTN_DB_SAVE ( MAX_GEOMETRIES )

!FROM VLIDORT_LP_SSCORR_NADIR:
!      DOUBLE PRECISION :: L_ZMAT_UP &
!          ( MAX_ATMOSWFS, MAX_GEOMETRIES, MAXLAYERS, MAXSTOKES, MAXSTOKES )
!      DOUBLE PRECISION :: L_ZMAT_DN &
!          ( MAX_ATMOSWFS, MAX_GEOMETRIES, MAXLAYERS, MAXSTOKES, MAXSTOKES )
!      DOUBLE PRECISION :: L_SS_CUMSOURCE &
!          ( MAX_ATMOSWFS, MAX_GEOMETRIES, MAXSTOKES )
!      DOUBLE PRECISION :: L_SSFDEL ( MAXLAYERS, MAX_ATMOSWFS )

!FROM VLIDORT_LP_SSCORR_OUTGOING:
      DOUBLE PRECISION :: LP_UP_MULTIPLIERS &
          ( MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )
      DOUBLE PRECISION :: LP_DN_MULTIPLIERS &
          ( MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )
      !DOUBLE PRECISION :: L_UP_LOSTRANS &
      !    ( MAXLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )
      DOUBLE PRECISION :: L_DN_LOSTRANS &
          ( MAXLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )
      DOUBLE PRECISION :: LP_UP_MULTIPLIERS_UT &
          ( MAX_PARTLAYERS, MAXLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )
      DOUBLE PRECISION :: LP_DN_MULTIPLIERS_UT &
          ( MAX_PARTLAYERS, MAXLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )
      !DOUBLE PRECISION :: L_UP_LOSTRANS_UT &
      !    ( MAX_PARTLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )
      DOUBLE PRECISION :: L_DN_LOSTRANS_UT &
          ( MAX_PARTLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )
      !DOUBLE PRECISION :: LP_BOA_ATTN &
      !    ( MAXLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )

!FROM VLIDORT_LAP_DBCORRECTION:
      DOUBLE PRECISION :: L_EXACTDB_SOURCE &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION :: L_DB_CUMSOURCE &
          ( MAX_GEOMETRIES, MAXSTOKES )

!FROM VLIDORT_LS_DBCORRECTION:
      DOUBLE PRECISION :: LS_EXACTDB_SOURCE &
          ( MAX_SURFACEWFS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION :: LS_DB_CUMSOURCE &
          ( MAX_GEOMETRIES, MAXSTOKES )

!FROM VFO_LPS_MASTER_INTERFACE:
      DOUBLE PRECISION :: FO_STOKES_SS &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: FO_STOKES_DB &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION :: FO_STOKES_DTA &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: FO_STOKES_DTS &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION :: FO_STOKES &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )

      DOUBLE PRECISION :: FO_PROFILEWF_SS &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS,&
            MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: FO_PROFILEWF_DB &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION :: FO_SURFACEWF_DB &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION :: FO_PROFILEWF_DTA &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS,&
            MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: FO_PROFILEWF_DTS &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS,&
            MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION :: FO_SURFACEWF_DTS &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION :: FO_PROFILEWF &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS,&
            MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: FO_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  Work type structures

      TYPE(VLIDORT_Work_Miscellanous) :: Misc
      TYPE(VLIDORT_Work_Thermal)      :: Therm
      TYPE(VLIDORT_Work_Multiplier)   :: Mult
      !TYPE(VLIDORT_Work_Corrections)  :: Corr
      !TYPE(VLIDORT_Work_FirstOrder)   :: Fo

      TYPE(VLIDORT_LinWork_Miscellanous) :: LAP_Misc
      TYPE(VLIDORT_LinWork_Thermal)      :: L_Therm
      TYPE(VLIDORT_LinWork_Multiplier)   :: LP_Mult

!  Intermediate quantities

      DOUBLE PRECISION :: PIMM_11 ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION :: PIMM_KM ( MAX_ALLSTRMS_P1 )

!  Weighting function indices

      INTEGER ::          LAYER_TO_VARY, N_PARAMETERS

!  Local error handling

      LOGICAL ::          FAIL

!  Test variables

      LOGICAL ::          DO_FDTEST=.FALSE.

!      LOGICAL ::          DO_DEBUG_INPUT=.TRUE.
      LOGICAL ::          DO_DEBUG_INPUT=.FALSE.

!  Initialize some variables
!  -------------------------

!  Main status

      STATUS_CALCULATION = VLIDORT_SUCCESS
      STATUS_INPUTCHECK  = VLIDORT_SUCCESS

!mick fix 6/29/11 - initialize "Input checks"
!  Input checks

      NCHECKMESSAGES = 0
      CHECKMESSAGES  = ' '
      ACTIONS        = ' '

!  Model calculation

      MESSAGE = ' '
      TRACE_1 = ' '
      TRACE_2 = ' '
      TRACE_3 = ' '

!  Local user indices

      LUM = 1
      LUA = 1

!  Check input dimensions. New for Version 2.7
!  -------------------------------------------

!  regular

      CALL VLIDORT_CHECK_INPUT_DIMS &
      ( VLIDORT_FixIn, VLIDORT_ModIn, &
        STATUS_SUB, NCHECKMESSAGES, CHECKMESSAGES, ACTIONS )

!  Exception handling

      IF ( STATUS_SUB .EQ. VLIDORT_SERIOUS ) THEN
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
        VLIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
        VLIDORT_Out%Status%TS_NCHECKMESSAGES    = NCHECKMESSAGES
        VLIDORT_Out%Status%TS_CHECKMESSAGES     = CHECKMESSAGES
        VLIDORT_Out%Status%TS_ACTIONS           = ACTIONS
        RETURN
      ENDIF

!  Linearized

      !CALL VLIDORT_L_CHECK_INPUT_DIMS &
      !( DO_COLUMN_LINEARIZATION,  DO_PROFILE_LINEARIZATION, & ! Flags
      !  DO_SURFACE_LINEARIZATION, DO_SLEAVE_WFS,            & ! Flags
      !  N_TOTALCOLUMN_WFS, N_TOTALPROFILE_WFS,              & ! Dimension-checking
      !  N_SURFACE_WFS, N_SLEAVE_WFS,                        & ! Dimension-checking
      !  STATUS_SUB, NCHECKMESSAGES, CHECKMESSAGES, ACTIONS )

      CALL VLIDORT_L_CHECK_INPUT_DIMS &
      ( VLIDORT_LinFixIn, VLIDORT_LinModIn, &
        STATUS_SUB, NCHECKMESSAGES, CHECKMESSAGES, ACTIONS )

!  Exception handling

      IF ( STATUS_SUB .EQ. VLIDORT_SERIOUS ) THEN
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
        VLIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
        VLIDORT_Out%Status%TS_NCHECKMESSAGES    = NCHECKMESSAGES
        VLIDORT_Out%Status%TS_CHECKMESSAGES     = CHECKMESSAGES
        VLIDORT_Out%Status%TS_ACTIONS           = ACTIONS
        RETURN
      ENDIF

!  ====================================
!  BEGIN COPY INPUTS TO LOCAL VARIABLES
!  ====================================

!  Fixed Boolean inputs

      DO_FULLRAD_MODE        = VLIDORT_FixIn%Bool%TS_DO_FULLRAD_MODE
      DO_SSCORR_TRUNCATION   = VLIDORT_FixIn%Bool%TS_DO_SSCORR_TRUNCATION

!  New 15 March 2012
      DO_SS_EXTERNAL         = VLIDORT_FixIn%Bool%TS_DO_SS_EXTERNAL

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

!  New 17 May 2012
      DO_SURFACE_LEAVING     = VLIDORT_FixIn%Bool%TS_DO_SURFACE_LEAVING
      DO_SL_ISOTROPIC        = VLIDORT_FixIn%Bool%TS_DO_SL_ISOTROPIC

!  Fixed Control inputs, Taylor parameter new, Version 2p7

      TAYLOR_ORDER     = VLIDORT_FixIn%Cont%TS_TAYLOR_ORDER 

      NSTOKES          = VLIDORT_FixIn%Cont%TS_NSTOKES
      NSTREAMS         = VLIDORT_FixIn%Cont%TS_NSTREAMS
      NLAYERS          = VLIDORT_FixIn%Cont%TS_NLAYERS
      NFINELAYERS      = VLIDORT_FixIn%Cont%TS_NFINELAYERS
      N_THERMAL_COEFFS = VLIDORT_FixIn%Cont%TS_N_THERMAL_COEFFS
      VLIDORT_ACCURACY = VLIDORT_FixIn%Cont%TS_VLIDORT_ACCURACY

      NLAYERS_NOMS     = VLIDORT_FixIn%Cont%TS_NLAYERS_NOMS
      NLAYERS_CUTOFF   = VLIDORT_FixIn%Cont%TS_NLAYERS_CUTOFF

!  Fixed Beam inputs

      FLUX_FACTOR = VLIDORT_FixIn%Sunrays%TS_FLUX_FACTOR

!  Fixed User Value inputs

      N_USER_LEVELS       = VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS

!  Fixed Chapman Function inputs

      HEIGHT_GRID(0:NLAYERS)      = &
        VLIDORT_FixIn%Chapman%TS_HEIGHT_GRID(0:NLAYERS)
      PRESSURE_GRID(0:NLAYERS)    = &
        VLIDORT_FixIn%Chapman%TS_PRESSURE_GRID(0:NLAYERS)
      TEMPERATURE_GRID(0:NLAYERS) = &
        VLIDORT_FixIn%Chapman%TS_TEMPERATURE_GRID(0:NLAYERS)
      FINEGRID(1:NLAYERS)           = &
        VLIDORT_FixIn%Chapman%TS_FINEGRID(1:NLAYERS)

      RFINDEX_PARAMETER     = VLIDORT_FixIn%Chapman%TS_RFINDEX_PARAMETER

!  Fixed Optical inputs

      DELTAU_VERT_INPUT(1:NLAYERS) = &
        VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT(1:NLAYERS)
      GREEKMAT_TOTAL_INPUT&
          (0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT,1:NLAYERS,:) = &
        VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT&
          (0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT,1:NLAYERS,:)
      THERMAL_BB_INPUT(0:NLAYERS)  = &
        VLIDORT_FixIn%Optical%TS_THERMAL_BB_INPUT(0:NLAYERS)

!  @@@ Rob Fix 1/31/11, TS_EMISSIVITY, TS_USER_EMISSIVITY are defined in
!                       in Structure VLIDORT_Sup%BRDF, so do not need
!                       to be copied here.
!      EMISSIVITY            = VLIDORT_FixIn%Optical%TS_EMISSIVITY
!      USER_EMISSIVITY       = VLIDORT_FixIn%Optical%TS_USER_EMISSIVITY

      LAMBERTIAN_ALBEDO     = VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO
      SURFBB                = VLIDORT_FixIn%Optical%TS_SURFACE_BB_INPUT

!  Rob, 28 March 2014. Superseded by LBBF Jacobians

!      LTE_DELTAU_VERT_INPUT(1:2,1:NLAYERS) = &
!        VLIDORT_FixIn%Optical%TS_LTE_DELTAU_VERT_INPUT(1:2,1:NLAYERS)
!      LTE_THERMAL_BB_INPUT(0:NLAYERS)      = &
!        VLIDORT_FixIn%Optical%TS_LTE_THERMAL_BB_INPUT(0:NLAYERS)

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

!  Modified Boolean inputs

      DO_SSCORR_NADIR         = VLIDORT_ModIn%MBool%TS_DO_SSCORR_NADIR
      DO_SSCORR_OUTGOING      = VLIDORT_ModIn%MBool%TS_DO_SSCORR_OUTGOING
      DO_FO_CALC              = VLIDORT_ModIn%MBool%TS_DO_FO_CALC  !New 02 Jul 2013
      DO_DOUBLE_CONVTEST      = VLIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST
      DO_SOLAR_SOURCES        = VLIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES
      DO_REFRACTIVE_GEOMETRY  = VLIDORT_ModIn%MBool%TS_DO_REFRACTIVE_GEOMETRY
      DO_CHAPMAN_FUNCTION     = VLIDORT_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION
      DO_RAYLEIGH_ONLY        = VLIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY
      !DO_ISOTROPIC_ONLY       = VLIDORT_ModIn%MBool%TS_DO_ISOTROPIC_ONLY
      !DO_NO_AZIMUTH           = VLIDORT_ModIn%MBool%TS_DO_NO_AZIMUTH
      !DO_ALL_FOURIER          = VLIDORT_ModIn%MBool%TS_DO_ALL_FOURIER
      DO_DELTAM_SCALING       = VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING
      DO_SOLUTION_SAVING      = VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING
      DO_BVP_TELESCOPING      = VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING
      !DO_USER_STREAMS         = VLIDORT_ModIn%MBool%TS_DO_USER_STREAMS
      DO_USER_VZANGLES        = VLIDORT_ModIn%MBool%TS_DO_USER_VZANGLES
      DO_ADDITIONAL_MVOUT     = VLIDORT_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT
      DO_MVOUT_ONLY           = VLIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY
      DO_THERMAL_TRANSONLY    = VLIDORT_ModIn%MBool%TS_DO_THERMAL_TRANSONLY
      DO_OBSERVATION_GEOMETRY = VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY

!  Modified Control inputs

      NGREEK_MOMENTS_INPUT = VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT

!  Modified Beam inputs

      !NBEAMS      = VLIDORT_ModIn%MSunrays%TS_NBEAMS
      N_SZANGLES  = VLIDORT_ModIn%MSunrays%TS_N_SZANGLES
      !BEAM_SZAS(1:NBEAMS)    = VLIDORT_ModIn%MSunrays%TS_BEAM_SZAS(1:NBEAMS)
      SZANGLES(1:N_SZANGLES) = VLIDORT_ModIn%MSunrays%TS_SZANGLES(1:N_SZANGLES)

!  Modified User Value inputs

      N_USER_RELAZMS      = VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS
      USER_RELAZMS(1:N_USER_RELAZMS)        = &
        VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(1:N_USER_RELAZMS)

      !N_USER_STREAMS     = VLIDORT_ModIn%MUserVal%TS_N_USER_STREAMS
      N_USER_VZANGLES     = VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES
      !USER_ANGLES(1:N_USER_STREAMS)        = &
      !  VLIDORT_ModIn%MUserVal%TS_USER_ANGLES_INPUT(1:N_USER_STREAMS)
      USER_VZANGLES(1:N_USER_VZANGLES)      = &
        VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(1:N_USER_VZANGLES)

      USER_LEVELS(1:N_USER_LEVELS)          = &
        VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(1:N_USER_LEVELS)

      GEOMETRY_SPECHEIGHT = VLIDORT_ModIn%MUserVal%TS_GEOMETRY_SPECHEIGHT

      N_USER_OBSGEOMS     = VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS
      USER_OBSGEOMS(1:N_USER_OBSGEOMS,1:3)  = &
        VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1:N_USER_OBSGEOMS,1:3)

!  Modified Chapman Function inputs

      EARTH_RADIUS        = VLIDORT_ModIn%MChapman%TS_EARTH_RADIUS
      !CHAPMAN_FACTORS     = VLIDORT_ModIn%MChapman%TS_CHAPMAN_FACTORS

!  Modified Optical inputs

      OMEGA_TOTAL_INPUT(1:NLAYERS) = &
        VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(1:NLAYERS)

!  BRDF Supplement Inputs

      EXACTDB_BRDFUNC&
          (:,1:N_USER_VZANGLES,1:N_USER_RELAZMS,1:N_SZANGLES) = &
        VLIDORT_Sup%BRDF%TS_EXACTDB_BRDFUNC&
          (:,1:N_USER_VZANGLES,1:N_USER_RELAZMS,1:N_SZANGLES)

      BRDF_F_0(0:MAXMOMENTS,:,1:NSTREAMS,1:N_SZANGLES)       = &
        VLIDORT_Sup%BRDF%TS_BRDF_F_0&
          (0:MAXMOMENTS,:,1:NSTREAMS,1:N_SZANGLES)
      BRDF_F(0:MAXMOMENTS,:,1:NSTREAMS,1:NSTREAMS)           = &
        VLIDORT_Sup%BRDF%TS_BRDF_F&
          (0:MAXMOMENTS,:,1:NSTREAMS,1:NSTREAMS)

      USER_BRDF_F_0&
          (0:MAXMOMENTS,:,1:N_USER_VZANGLES,1:N_SZANGLES)     = &
        VLIDORT_Sup%BRDF%TS_USER_BRDF_F_0&
          (0:MAXMOMENTS,:,1:N_USER_VZANGLES,1:N_SZANGLES)
      USER_BRDF_F&
          (0:MAXMOMENTS,:,1:N_USER_VZANGLES,1:NSTREAMS)       = &
        VLIDORT_Sup%BRDF%TS_USER_BRDF_F&
          (0:MAXMOMENTS,:,1:N_USER_VZANGLES,1:NSTREAMS)

!  @@@ Rob fix 1/31/11, Emissivities from earlier structure
!                       were not copied --> wrong answers for BRDF cases
!        Lambertian case was OK, as used internal definitions

      EMISSIVITY(1:NSTOKES,1:NSTREAMS)            = &
        VLIDORT_Sup%BRDF%TS_EMISSIVITY(1:NSTOKES,1:NSTREAMS)
      USER_EMISSIVITY(1:NSTOKES,1:N_USER_VZANGLES) = &
        VLIDORT_Sup%BRDF%TS_USER_EMISSIVITY(1:NSTOKES,1:N_USER_VZANGLES)

!  Surface-Leaving Supplement inputs
!    (This code introduced 17 May 2012)

      SLTERM_ISOTROPIC(1:NSTOKES,1:N_SZANGLES)                       = &
        VLIDORT_Sup%SLEAVE%TS_SLTERM_ISOTROPIC(1:NSTOKES,1:N_SZANGLES)
      SLTERM_USERANGLES&
          (1:NSTOKES,1:N_USER_VZANGLES,1:N_USER_RELAZMS,1:N_SZANGLES) = &
        VLIDORT_Sup%SLEAVE%TS_SLTERM_USERANGLES&
          (1:NSTOKES,1:N_USER_VZANGLES,1:N_USER_RELAZMS,1:N_SZANGLES)
      SLTERM_F_0&
          (0:MAXMOMENTS,1:NSTOKES,1:NSTREAMS,1:N_SZANGLES)           = &
        VLIDORT_Sup%SLEAVE%TS_SLTERM_F_0&
          (0:MAXMOMENTS,1:NSTOKES,1:NSTREAMS,1:N_SZANGLES)
      USER_SLTERM_F_0&
          (0:MAXMOMENTS,1:NSTOKES,1:N_USER_VZANGLES,1:N_SZANGLES)     = &
        VLIDORT_Sup%SLEAVE%TS_USER_SLTERM_F_0&
          (0:MAXMOMENTS,1:NSTOKES,1:N_USER_VZANGLES,1:N_SZANGLES)

!  new 12 March 2012
!   IF SS results already available copy them !
!  SS Inputs

      IF ( DO_SS_EXTERNAL ) THEN
        STOKES_SS(1:N_USER_LEVELS,:,1:NSTOKES,:) = &
          VLIDORT_Sup%SS%TS_STOKES_SS
        STOKES_DB(1:N_USER_LEVELS,:,1:NSTOKES) = &
          VLIDORT_Sup%SS%TS_STOKES_DB
      ENDIF

!  Fixed Linearized Control inputs

      LAYER_VARY_FLAG(1:NLAYERS)   = &
        VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG(1:NLAYERS)
      LAYER_VARY_NUMBER(1:NLAYERS) = &
        VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER(1:NLAYERS)

      N_TOTALCOLUMN_WFS        = VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS
      N_TOTALPROFILE_WFS       = VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS
      N_SURFACE_WFS            = VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS
      N_SLEAVE_WFS             = VLIDORT_LinFixIn%Cont%TS_N_SLEAVE_WFS

      N_TOTALSURFACE_WFS = N_SURFACE_WFS + N_SLEAVE_WFS

      COLUMNWF_NAMES           = VLIDORT_LinFixIn%Cont%TS_COLUMNWF_NAMES
      PROFILEWF_NAMES          = VLIDORT_LinFixIn%Cont%TS_PROFILEWF_NAMES

!  Fixed Linearized Optical inputs

      L_DELTAU_VERT_INPUT(1:N_TOTALPROFILE_WFS,1:NLAYERS) = &
        VLIDORT_LinFixIn%Optical%TS_L_DELTAU_VERT_INPUT&
          (1:N_TOTALPROFILE_WFS,1:NLAYERS)
      L_OMEGA_TOTAL_INPUT(1:N_TOTALPROFILE_WFS,1:NLAYERS) = &
        VLIDORT_LinFixIn%Optical%TS_L_OMEGA_TOTAL_INPUT&
          (1:N_TOTALPROFILE_WFS,1:NLAYERS)
      L_GREEKMAT_TOTAL_INPUT&
          (1:N_TOTALPROFILE_WFS,0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT,1:NLAYERS,:) = &
        VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT&
          (1:N_TOTALPROFILE_WFS,0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT,1:NLAYERS,:)

!  Modified Linearized Control inputs
!   (First three were formerly in Fixed Lin Control)

      DO_SIMULATION_ONLY       = VLIDORT_LinModIn%MCont%TS_DO_SIMULATION_ONLY

!  Superseded in Version 2.7 by LBBF flags

!      DO_LTE_LINEARIZATION    = VLIDORT_LinFixIn%Cont%TS_DO_LTE_LINEARIZATION
!      DO_SURFBB_LINEARIZATION = VLIDORT_LinFixIn%Cont%TS_DO_SURFBB_LINEARIZATION
      DO_ATMOS_LBBF   = VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LBBF
      DO_SURFACE_LBBF = VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LBBF

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

      DO_SLEAVE_WFS            = &
        VLIDORT_LinModIn%MCont%TS_DO_SLEAVE_WFS

!  @@@ Rob Fix 1/31/11, TS_LS_EMISSIVITY, TS_LS_USER_EMISSIVITY are defined in
!                       in Structure VLIDORT_LinSup%BRDF, so do not need
!                       to be copied here.
!      LS_EMISSIVITY          = VLIDORT_LinFixIn%Optical%TS_LS_EMISSIVITY
!      LS_USER_EMISSIVITY     = VLIDORT_LinFixIn%Optical%TS_LS_USER_EMISSIVITY

!  BRDF Linearized Supplement inputs

      LS_EXACTDB_BRDFUNC&
          (1:N_SURFACE_WFS,:,1:N_USER_VZANGLES,1:N_USER_RELAZMS,1:N_SZANGLES) = &
        VLIDORT_LinSup%BRDF%TS_LS_EXACTDB_BRDFUNC&
          (1:N_SURFACE_WFS,:,1:N_USER_VZANGLES,1:N_USER_RELAZMS,1:N_SZANGLES)

      LS_BRDF_F_0&
          (1:N_SURFACE_WFS,0:MAXMOMENTS,:,1:NSTREAMS,1:N_SZANGLES) = &
        VLIDORT_LinSup%BRDF%TS_LS_BRDF_F_0&
          (1:N_SURFACE_WFS,0:MAXMOMENTS,:,1:NSTREAMS,1:N_SZANGLES)
      LS_BRDF_F&
          (1:N_SURFACE_WFS,0:MAXMOMENTS,:,1:NSTREAMS,1:NSTREAMS)   = &
        VLIDORT_LinSup%BRDF%TS_LS_BRDF_F&
          (1:N_SURFACE_WFS,0:MAXMOMENTS,:,1:NSTREAMS,1:NSTREAMS)

      LS_USER_BRDF_F_0&
          (1:N_SURFACE_WFS,0:MAXMOMENTS,:,1:N_USER_VZANGLES,1:N_SZANGLES) = &
        VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F_0&
          (1:N_SURFACE_WFS,0:MAXMOMENTS,:,1:N_USER_VZANGLES,1:N_SZANGLES)
      LS_USER_BRDF_F&
          (1:N_SURFACE_WFS,0:MAXMOMENTS,:,1:N_USER_VZANGLES,1:NSTREAMS)   = &
        VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F&
          (1:N_SURFACE_WFS,0:MAXMOMENTS,:,1:N_USER_VZANGLES,1:NSTREAMS)

!  @@@ Rob fix 1/31/11, Linearized Emissivities from earlier structure
!                       were not copied --> wrong answers for BRDF cases
!        Lambertian case was OK, as used internal definitions

      LS_EMISSIVITY(1:N_SURFACE_WFS,1:NSTOKES,1:NSTREAMS)            = &
        VLIDORT_LinSup%BRDF%TS_LS_EMISSIVITY&
          (1:N_SURFACE_WFS,1:NSTOKES,1:NSTREAMS)
      LS_USER_EMISSIVITY(1:N_SURFACE_WFS,1:NSTOKES,1:N_USER_VZANGLES) = &
        VLIDORT_LinSup%BRDF%TS_LS_USER_EMISSIVITY&
          (1:N_SURFACE_WFS,1:NSTOKES,1:N_USER_VZANGLES)

!  new 12 March 2012
!   IF SS results already available copy them !
!  Linearized SS Inputs

      IF ( DO_SS_EXTERNAL ) THEN
         !SS atmosphere weighting functions

         PROFILEWF_SS&
             (1:N_TOTALPROFILE_WFS,1:NLAYERS,1:N_USER_LEVELS,:,1:NSTOKES,:) = &
           VLIDORT_LinSup%SS%Prof%TS_PROFILEWF_SS&
             (1:N_TOTALPROFILE_WFS,1:NLAYERS,1:N_USER_LEVELS,:,1:NSTOKES,:)
         PROFILEWF_DB&
             (1:N_TOTALPROFILE_WFS,1:NLAYERS,1:N_USER_LEVELS,:,1:NSTOKES)                  = &
           VLIDORT_LinSup%SS%Prof%TS_PROFILEWF_DB&
             (1:N_TOTALPROFILE_WFS,1:NLAYERS,1:N_USER_LEVELS,:,1:NSTOKES)

         !SS surface weighting functions

         !SURFACEWF_SS&
         !    (1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,:,1:NSTOKES,:) = &
         !  VLIDORT_LinSup%SS%Surf%TS_SURFACEWF_SS
         !    (1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,:,1:NSTOKES,:)
         SURFACEWF_DB&
             (1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,:,1:NSTOKES)  = &
           VLIDORT_LinSup%SS%Surf%TS_SURFACEWF_DB&
             (1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,:,1:NSTOKES)

      ENDIF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@ Addition of SLEAVE WF stuff, R. Spurr, 22 August 2012 @@@@@@@@@

      LSSL_SLTERM_ISOTROPIC&
          (1:N_SLEAVE_WFS,1:NSTOKES,1:N_SZANGLES) = &
        VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_ISOTROPIC&
          (1:N_SLEAVE_WFS,1:NSTOKES,1:N_SZANGLES)
      LSSL_SLTERM_USERANGLES&
          (1:N_SLEAVE_WFS,1:NSTOKES,1:N_USER_VZANGLES,1:N_USER_RELAZMS,1:N_SZANGLES) = &
        VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_USERANGLES&
          (1:N_SLEAVE_WFS,1:NSTOKES,1:N_USER_VZANGLES,1:N_USER_RELAZMS,1:N_SZANGLES)
      LSSL_SLTERM_F_0&
          (1:N_SLEAVE_WFS,0:MAXMOMENTS,1:NSTOKES,1:NSTREAMS,1:N_SZANGLES) = &
        VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_F_0&
          (1:N_SLEAVE_WFS,0:MAXMOMENTS,1:NSTOKES,1:NSTREAMS,1:N_SZANGLES)
      LSSL_USER_SLTERM_F_0&
          (1:N_SLEAVE_WFS,0:MAXMOMENTS,1:NSTOKES,1:N_USER_VZANGLES,1:N_SZANGLES) = &
        VLIDORT_LinSup%SLEAVE%TS_LSSL_USER_SLTERM_F_0&
          (1:N_SLEAVE_WFS,0:MAXMOMENTS,1:NSTOKES,1:N_USER_VZANGLES,1:N_SZANGLES)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  ==================================
!  END COPY INPUTS TO LOCAL VARIABLES
!  ==================================

!  VLIDORT input debug

      IF (DO_DEBUG_INPUT) THEN
        CALL VLIDORT_DEBUG_INPUT_MASTER()
        IF (DO_PROFILE_LINEARIZATION .OR. &
            DO_SURFACE_LINEARIZATION) &
          CALL VLIDORT_DEBUG_LIN_INPUT_MASTER()
      END IF

!  initialize outputs
!  ------------------

!mick fix
!  Main outputs

      !Radiances and fluxes
      STOKES(1:N_USER_LEVELS,:,1:NSTOKES,:) = ZERO
      MEAN_STOKES(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:) = ZERO
      FLUX_STOKES(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:) = ZERO
      MEAN_DIRECT(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES)   = ZERO
      FLUX_DIRECT(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES)   = ZERO

!  New 15 March 2012
      IF ( .NOT. DO_SS_EXTERNAL ) THEN
         STOKES_SS(1:N_USER_LEVELS,:,1:NSTOKES,:) = ZERO
         STOKES_DB(1:N_USER_LEVELS,:,1:NSTOKES)   = ZERO
      ENDIF

      !Misc
      FOURIER_SAVED(1:N_SZANGLES) = 0
      N_GEOMETRIES = 0
      SZA_OFFSETS(1:N_SZANGLES)   = 0
      VZA_OFFSETS(1:N_SZANGLES,1:N_USER_VZANGLES) = 0

      !Profile weighting functions
      PROFILEWF&
        (1:N_TOTALPROFILE_WFS,1:NLAYERS,1:N_USER_LEVELS,&
         :,1:NSTOKES,:) = ZERO

      !Total mean intensity and flux weighting functions
      MINT_PROFILEWF&
        (1:N_TOTALPROFILE_WFS,1:NLAYERS,1:N_USER_LEVELS,&
         1:N_SZANGLES,1:NSTOKES,:) = ZERO
      FLUX_PROFILEWF&
        (1:N_TOTALPROFILE_WFS,1:NLAYERS,1:N_USER_LEVELS,&
         1:N_SZANGLES,1:NSTOKES,:) = ZERO

      !Mean intensity and flux weighting functions.  Direct beam only.
      MINT_PROFILEWF_DIRECT&
        (1:N_TOTALPROFILE_WFS,1:NLAYERS,1:N_USER_LEVELS,&
         1:N_SZANGLES,1:NSTOKES) = ZERO
      FLUX_PROFILEWF_DIRECT&
        (1:N_TOTALPROFILE_WFS,1:NLAYERS,1:N_USER_LEVELS,&
         1:N_SZANGLES,1:NSTOKES) = ZERO

!  Removed. 28 March 2014
!      !LTE profile weighting functions
!      LTE_ATMOSWF(0:NLAYERS,1:N_USER_LEVELS,1:N_USER_VZANGLES,:) = ZERO

      !Surface weighting functions
      SURFACEWF(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,:,1:NSTOKES,:) = ZERO
      MINT_SURFACEWF&
        (1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:) = ZERO
      FLUX_SURFACEWF&
        (1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:) = ZERO

!  New 15 March 2012
      IF ( .NOT. DO_SS_EXTERNAL ) THEN
         !SS atmosphere weighting functions
         PROFILEWF_SS&
           (1:N_TOTALPROFILE_WFS,1:NLAYERS,1:N_USER_LEVELS,&
           :,1:NSTOKES,:) = ZERO
         PROFILEWF_DB&
           (1:N_TOTALPROFILE_WFS,1:NLAYERS,1:N_USER_LEVELS,&
           :,1:NSTOKES)   = ZERO

         !SS surface weighting functions
         !SURFACEWF_SS = ZERO
         SURFACEWF_DB(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,:,1:NSTOKES) = ZERO
      ENDIF

!  New 28 March 2014. BLACKBODY Linearization, Version 2.7
!  -------------------------------------------------------

      ABBWFS_JACOBIANS(1:N_USER_LEVELS,:,:,:,:)  = ZERO
      ABBWFS_FLUXES   (1:N_USER_LEVELS,:,:,:,:)  = ZERO
      SBBWFS_JACOBIANS(1:N_USER_LEVELS,:,:,:)    = ZERO
      SBBWFS_FLUXES   (1:N_USER_LEVELS,:,:,:)    = ZERO

!  Single scatter correction: flux multiplier
!    Now always F / 4pi

      SS_FLUX_MULTIPLIER = FLUX_FACTOR / PI4

!  Check input
!  -----------

!  Standard input variables check

!    Major revision of I/O output list, 25 October 2012.
!     ---- Observational Geometry control, New, 25 October 2012
!     ---- Automatic setting of NBEAMS, N_USER_VZANGLES, N_USER_RELAZMS, DO_USER_STREAMS

!  %% DO_FULLRAD_MODE argument added (First line). R. Spurr, 05 March 2013
!  %%   Needed to ensure MS-only output in all cases when flagged

      CALL VLIDORT_CHECK_INPUT ( &
        DO_FULLRAD_MODE, DO_SSFULL, DO_PLANE_PARALLEL, &    ! %%
        DO_UPWELLING, DO_DNWELLING, &
        DO_QUAD_OUTPUT, DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY, &
        TAYLOR_ORDER, NSTOKES, NSTREAMS, NLAYERS, &   ! Version 2p7.
        NFINELAYERS, N_USER_LEVELS, &
        HEIGHT_GRID, OMEGA_TOTAL_INPUT, &
        GREEKMAT_TOTAL_INPUT, DO_LAMBERTIAN_SURFACE, &
        DO_THERMAL_EMISSION, &
        DO_SPECIALIST_OPTION_2, DO_SPECIALIST_OPTION_3, &
        NLAYERS_NOMS, NLAYERS_CUTOFF, &
        DO_TOA_CONTRIBS, DO_NO_AZIMUTH, DO_SS_EXTERNAL, &
        DO_SSCORR_NADIR, DO_SSCORR_OUTGOING, DO_FO_CALC, &
        DO_SOLAR_SOURCES, DO_REFRACTIVE_GEOMETRY, &
        DO_CHAPMAN_FUNCTION, DO_RAYLEIGH_ONLY, &
        DO_DELTAM_SCALING, DO_SOLUTION_SAVING, &
        DO_BVP_TELESCOPING, DO_USER_VZANGLES, &
        DO_OBSERVATION_GEOMETRY, &
        NGREEK_MOMENTS_INPUT, N_SZANGLES, SZANGLES, &
        EARTH_RADIUS, GEOMETRY_SPECHEIGHT, &
        N_USER_RELAZMS, USER_RELAZMS, &
        N_USER_VZANGLES, USER_VZANGLES, &
        USER_LEVELS, &
        N_USER_OBSGEOMS, USER_OBSGEOMS, &
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
          DO_SIMULATION_ONLY, N_SURFACE_WFS, &
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
!mick hold - 9/26/2012
!      IF ( DO_SOLAR_SOURCES ) THEN

      IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
        CALL MULTI_OUTGOING_ADJUSTGEOM                                &
         ( MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS,         & ! Input
           N_USER_VZANGLES,   N_SZANGLES,   N_USER_RELAZMS,           & ! Input
           HEIGHT_GRID(NLAYERS), MODIFIED_ERADIUS, ADJUST_SURFACE,    & ! Input
           USER_VZANGLES, SZANGLES, USER_RELAZMS,                     & ! Input
           USER_VZANGLES_ADJUST, SZANGLES_ADJUST, USER_RELAZMS_ADJUST,& ! Output
           FAIL, MESSAGE, TRACE_1 )                                     ! Output
      ELSE
        CALL OBSGEOM_OUTGOING_ADJUSTGEOM                               &
         ( MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS,          & ! Input
           N_USER_VZANGLES,   N_SZANGLES,   N_USER_RELAZMS,            & ! Input
           HEIGHT_GRID(NLAYERS), MODIFIED_ERADIUS, ADJUST_SURFACE,     & ! Input
           USER_VZANGLES, SZANGLES, USER_RELAZMS,                      & ! Input
           USER_VZANGLES_ADJUST, SZANGLES_ADJUST, USER_RELAZMS_ADJUST, & ! Output
           FAIL, MESSAGE, TRACE_1 )                                      ! Output
      ENDIF

!      ELSE
!        CALL LOSONLY_OUTGOING_ADJUSTGEOM                           &
!         ( MAX_USER_VZANGLES, N_USER_VZANGLES,                     & ! Input
!           HEIGHT_GRID(NLAYERS), MODIFIED_ERADIUS, ADJUST_SURFACE, & ! Input
!           USER_VZANGLES,                                          & ! Input
!           USER_VZANGLES_ADJUST,                                   & ! Output
!           FAIL, MESSAGE, TRACE_1 )                                  ! Output
!      ENDIF

!  Update exception handling. October 2010 2p4RTC

      if ( fail ) then
        TRACE_2 = ' Failure in multi_outgoing_adjustgeom'
        TRACE_3 = ' ** VLIDORT_LPS_MASTER '
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
            TRACE_2 = 'Direct call in VLIDORT_LPS_MASTER'
            TRACE_3 = ' ** VLIDORT_LPS_MASTER '
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
            DO_ATMOS_LBBF, DO_SURFACE_LBBF, PROFILEWF_NAMES, &
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
        DO_UPWELLING, DO_DNWELLING, DO_USER_VZANGLES, DO_OBSERVATION_GEOMETRY, &
        NSTOKES, NSTREAMS, NLAYERS, NGREEK_MOMENTS_INPUT, &
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
!        TRACE_2 = 'Derive_Input Call in VLIDORT_LPS_MASTER'
!        TRACE_3 = ' ** VLIDORT_LPS_MASTER'
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

!!  New 02 Jul 2013
      IF (DO_FO_CALC .AND. DO_PARTLAYERS) THEN
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
        VLIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
        NCHECKMESSAGES = NCHECKMESSAGES + 1
        CHECKMESSAGES(NCHECKMESSAGES) &
          = ' Internal SS calculation using FO code: Must only have whole layer (layer boundary) output'
        ACTIONS(NCHECKMESSAGES) &
          = ' Set output for only layer boundaries'
      ENDIF

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
!   This section revised for the Observational Geometry option

      IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
        N_VIEWING    = N_USER_VZANGLES * LOCAL_N_USERAZM
        N_GEOMETRIES = NSOURCES * N_VIEWING
        DO IBEAM = 1, NBEAMS
          SZA_OFFSETS(IBEAM) = N_VIEWING * ( IBEAM - 1 )
          DO UM = 1, N_USER_VZANGLES
            VZA_OFFSETS(IBEAM,UM) = &
                   SZA_OFFSETS(IBEAM) + LOCAL_N_USERAZM * (UM - 1)
          END DO
        END DO
      ELSE
        N_VIEWING    = N_USER_OBSGEOMS
        N_GEOMETRIES = N_USER_OBSGEOMS
        SZA_OFFSETS = 0 ; VZA_OFFSETS = 0
      ENDIF

!  #################
!  Set up operations
!  #################

!  Setups for Fourier = 0:
!  -----------------------

!  VLIDORT_MISCSETUPS    : Delta-M, average-secant formulation, transmittances
!  THERMAL_SETUP         : Coefficients, direct multipliers
!  EMULT_MASTER /        : Beam source function multipliers.  Not required for the
!    EMULT_MASTER_OBSGEO     Full SS calculation in outgoing mode

!  VLIDORT_LAP_MISCSETUPS
!  THERMAL_SETUP_PLUS
!  LP_EMULT_MASTER /
!    LP_EMULT_MASTER_OBSGEO

!  30 January 2008
!  Telescoping setup is done in DERIVE_INPUTS (unlike LIDORT scalar code)

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
          N_USER_VZANGLES, DO_USER_VZANGLES, &
          DO_OBSERVATION_GEOMETRY, &
          USER_SECANTS, N_PARTLAYERS, &
          STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, &
          OMEGA_TOTAL, DELTAU_VERT, &
          PARTAU_VERT, GREEKMAT_TOTAL, &
          TAUGRID, DELTAU_SLANT, SOLARBEAM_BOATRANS, &   ! Rob fix 11/17/14 added argument
          TRUNC_FACTOR, FAC1, &
          OMEGA_GREEK, LAYER_PIS_CUTOFF, &
          TRANS_SOLAR_BEAM, DO_REFLECTED_DIRECTBEAM, &
          INITIAL_TRANS, AVERAGE_SECANT, LOCAL_CSZA, &
          T_DELT_DISORDS, T_DISORDS_UTUP, &
          T_DISORDS_UTDN, T_DELT_MUBAR, &
          T_UTDN_MUBAR, T_UTUP_MUBAR, &
          T_DELT_USERM, T_UTDN_USERM, &
          T_UTUP_USERM, CUMTRANS, ITRANS_USERM )

        CALL VLIDORT_PACK_MISC ( &
          NSTOKES, NSTREAMS, NLAYERS, N_PARTLAYERS, N_USER_VZANGLES, N_SZANGLES,& ! Input 
          NMOMENTS,                                                             & ! Input
          DELTAU_VERT, PARTAU_VERT, DELTAU_SLANT, OMEGA_GREEK,                  & ! Input
          LAYER_PIS_CUTOFF, TRANS_SOLAR_BEAM, DO_REFLECTED_DIRECTBEAM,          & ! Input
          T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,                       & ! Input
          T_DELT_MUBAR, T_UTDN_MUBAR, T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM, & ! Input
          CUMTRANS, INITIAL_TRANS, AVERAGE_SECANT, LOCAL_CSZA,                  & ! Input
          Misc )                                                                  ! Output

        CALL VLIDORT_LAP_MISCSETUPS ( &
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
          N_USER_VZANGLES, DO_USER_VZANGLES, &
          USER_SECANTS, N_PARTLAYERS, &
          PARTLAYERS_LAYERIDX, &
          DELTAU_VERT, PARTAU_VERT, T_DELT_DISORDS, &
          T_DISORDS_UTUP, T_DISORDS_UTDN, &
          T_DELT_MUBAR, T_UTDN_MUBAR, &
          T_DELT_USERM, T_UTDN_USERM, &
          T_UTUP_USERM, AVERAGE_SECANT, &
          L_OMEGA_TOTAL, L_DELTAU_VERT, &
          L_GREEKMAT_TOTAL, L_DELTAU_SLANT, &
          DO_SCATMAT_VARIATION, L_TRUNC_FACTOR, &
          L_OMEGA_GREEK, LP_AVERAGE_SECANT, &
          LP_INITIAL_TRANS, L_T_DELT_DISORDS, &
          L_T_DISORDS_UTDN, L_T_DISORDS_UTUP, &
          LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR, &
          L_T_DELT_USERM, L_T_UTDN_USERM, &
          L_T_UTUP_USERM )

        CALL VLIDORT_PACK_LAP_MISC ( &
          NSTOKES, NSTOKES_SQ, NSTREAMS, NLAYERS, N_PARTLAYERS,      & ! Input
          N_SZANGLES, N_USER_VZANGLES, NMOMENTS, N_TOTALPROFILE_WFS, & ! Input
          L_DELTAU_VERT, L_OMEGA_GREEK,                              & ! Input
          LP_AVERAGE_SECANT, LP_INITIAL_TRANS,                       & ! Input
          L_T_DELT_DISORDS, L_T_DISORDS_UTDN, L_T_DISORDS_UTUP,      & ! Input
          LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR,                          & ! Input
          L_T_DELT_USERM, L_T_UTDN_USERM, L_T_UTUP_USERM,            & ! Input
          LAP_Misc )                                                   ! Output

        IF ( DO_THERMAL_EMISSION ) THEN
          CALL THERMAL_SETUP_PLUS ( &
            DO_UPWELLING, DO_DNWELLING, DO_MSMODE_THERMAL, &
            NLAYERS, N_THERMAL_COEFFS, &
            THERMAL_BB_INPUT, DO_THERMAL_TRANSONLY, &
            N_USER_VZANGLES, DO_USER_VZANGLES, &
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

          CALL VLIDORT_PACK_THERM ( &
            NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_VZANGLES, & ! Input
            THERMCOEFFS, DELTAU_POWER, XTAU_POWER, TCOM1,             & ! Input
            T_DIRECT_UP, T_DIRECT_DN, T_UT_DIRECT_UP, T_UT_DIRECT_DN, & ! Input
            Therm )                                                     ! Output

          CALL VLIDORT_PACK_L_THERM ( &
            NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_VZANGLES, N_TOTALPROFILE_WFS, & ! Input
            L_THERMCOEFFS, L_DELTAU_POWER, L_XTAU_POWER, L_TCOM1,                         & ! Input
            L_T_DIRECT_UP, L_T_DIRECT_DN, L_T_UT_DIRECT_UP, L_T_UT_DIRECT_DN,             & ! Input
            L_Therm )                                                                       ! Output
        ENDIF

        IF ( DO_SOLAR_SOURCES ) THEN
          IF (.NOT.DO_SSFULL .OR. (DO_SSFULL.AND.DO_SSCORR_NADIR)) THEN
            IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
              CALL EMULT_MASTER ( &
                DO_UPWELLING, DO_DNWELLING, TAYLOR_ORDER, & ! Add, Version 2.7
                NLAYERS, &
                LAYER_PIS_CUTOFF, NBEAMS, &
                N_USER_VZANGLES, &
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

              CALL VLIDORT_PACK_MULT ( &
                NLAYERS, N_PARTLAYERS, N_SZANGLES, N_USER_VZANGLES, & ! Input
                EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN,       & ! Input
                Mult )                                                ! Output

              CALL LP_EMULT_MASTER ( &
                DO_UPWELLING, DO_DNWELLING, TAYLOR_ORDER, & ! Add, Version 2.7
                NLAYERS, N_USER_LEVELS, &
                PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
                PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, &
                STERM_LAYERMASK_DN, &
                EMULT_UP, EMULT_DN, &
                LAYER_VARY_FLAG, &
                LAYER_VARY_NUMBER, DO_PLANE_PARALLEL, &
                LAYER_PIS_CUTOFF, NBEAMS, &
                N_USER_VZANGLES, &
                T_DELT_MUBAR, T_DELT_USERM, &
                ITRANS_USERM, SIGMA_P, &
                LP_AVERAGE_SECANT, LP_INITIAL_TRANS, &
                LP_T_DELT_MUBAR, L_T_DELT_USERM, &
                USER_SECANTS, &
                DELTAU_VERT, EMULT_HOPRULE, SIGMA_M, &
                L_DELTAU_VERT, T_UTUP_USERM, T_UTDN_USERM, & ! Rob Fix 5/10/13 added
                UT_EMULT_UP, &
                LP_T_UTDN_MUBAR, L_T_UTUP_USERM, &
                PARTAU_VERT, UT_EMULT_DN, &
                L_T_UTDN_USERM, &
                LP_EMULT_UP, LP_EMULT_DN, &
                LP_UT_EMULT_UP, LP_UT_EMULT_DN )

              CALL VLIDORT_PACK_LP_MULT ( &
                NLAYERS, N_PARTLAYERS, N_SZANGLES, N_USER_VZANGLES, N_TOTALPROFILE_WFS, & ! Input
                LP_EMULT_UP, LP_EMULT_DN, LP_UT_EMULT_UP, LP_UT_EMULT_DN,               & ! Input
                LP_Mult )                                                                 ! Output
            ELSE
              CALL EMULT_MASTER_OBSGEO ( &
                DO_UPWELLING, DO_DNWELLING, TAYLOR_ORDER, & ! Add, Version 2.7
                NLAYERS, &
                LAYER_PIS_CUTOFF, NBEAMS, &
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

              CALL VLIDORT_PACK_MULT ( &
                NLAYERS, N_PARTLAYERS, N_SZANGLES, N_USER_VZANGLES, & ! Input
                EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN,       & ! Input
                Mult )                                                ! Output

              CALL LP_EMULT_MASTER_OBSGEO ( &
                DO_UPWELLING, DO_DNWELLING, TAYLOR_ORDER, & ! Add, Version 2.7
                NLAYERS, N_USER_LEVELS, &
                PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
                PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, &
                STERM_LAYERMASK_DN, &
                EMULT_UP, EMULT_DN, &
                LAYER_VARY_FLAG, &
                LAYER_VARY_NUMBER, DO_PLANE_PARALLEL, &
                LAYER_PIS_CUTOFF, NBEAMS, &
                T_DELT_MUBAR, T_DELT_USERM, &
                ITRANS_USERM, SIGMA_P, &
                LP_AVERAGE_SECANT, LP_INITIAL_TRANS, &
                LP_T_DELT_MUBAR, L_T_DELT_USERM, &
                USER_SECANTS, &
                DELTAU_VERT, EMULT_HOPRULE, SIGMA_M, &
                L_DELTAU_VERT, T_UTUP_USERM, T_UTDN_USERM, &  ! Rob Fix 5/10/13 added
                UT_EMULT_UP, &
                LP_T_UTDN_MUBAR, L_T_UTUP_USERM, &
                PARTAU_VERT, UT_EMULT_DN, &
                L_T_UTDN_USERM, &
                LP_EMULT_UP, LP_EMULT_DN, &
                LP_UT_EMULT_UP, LP_UT_EMULT_DN )

              CALL VLIDORT_PACK_LP_MULT ( &
                NLAYERS, N_PARTLAYERS, N_SZANGLES, N_USER_VZANGLES,       & ! Input
                N_TOTALPROFILE_WFS,                                       & ! Input
                LP_EMULT_UP, LP_EMULT_DN, LP_UT_EMULT_UP, LP_UT_EMULT_DN, & ! Input
                LP_Mult )                                                   ! Output
            ENDIF
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
          N_USER_VZANGLES, DO_USER_VZANGLES, &
          DO_OBSERVATION_GEOMETRY, &
          USER_SECANTS, N_PARTLAYERS, &
          STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, &
          OMEGA_TOTAL, DELTAU_VERT, &
          PARTAU_VERT, GREEKMAT_TOTAL, &
          TAUGRID, DELTAU_SLANT, SOLARBEAM_BOATRANS, &   ! Rob fix 11/17/14 added argument
          TRUNC_FACTOR, FAC1, &
          OMEGA_GREEK, LAYER_PIS_CUTOFF, &
          TRANS_SOLAR_BEAM, DO_REFLECTED_DIRECTBEAM, &
          INITIAL_TRANS, AVERAGE_SECANT, LOCAL_CSZA, &
          T_DELT_DISORDS, T_DISORDS_UTUP, &
          T_DISORDS_UTDN, T_DELT_MUBAR, &
          T_UTDN_MUBAR, T_UTUP_MUBAR, &
          T_DELT_USERM, T_UTDN_USERM, &
          T_UTUP_USERM, CUMTRANS, ITRANS_USERM )

        CALL VLIDORT_PACK_MISC ( &
          NSTOKES, NSTREAMS, NLAYERS, N_PARTLAYERS, N_USER_VZANGLES, N_SZANGLES,& ! Input 
          NMOMENTS,                                                             & ! Input
          DELTAU_VERT, PARTAU_VERT, DELTAU_SLANT, OMEGA_GREEK,                  & ! Input
          LAYER_PIS_CUTOFF, TRANS_SOLAR_BEAM, DO_REFLECTED_DIRECTBEAM,          & ! Input
          T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,                       & ! Input
          T_DELT_MUBAR, T_UTDN_MUBAR, T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM, & ! Input
          CUMTRANS, INITIAL_TRANS, AVERAGE_SECANT, LOCAL_CSZA,                  & ! Input
          Misc )                                                                  ! Output

        IF ( DO_THERMAL_EMISSION ) THEN
          CALL THERMAL_SETUP ( &
            DO_UPWELLING, DO_DNWELLING, DO_MSMODE_THERMAL, &
            NLAYERS, N_THERMAL_COEFFS, &
            THERMAL_BB_INPUT, DO_THERMAL_TRANSONLY, &
            N_USER_VZANGLES, DO_USER_VZANGLES, &
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

          CALL VLIDORT_PACK_THERM ( &
            NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_VZANGLES, & ! Input
            THERMCOEFFS, DELTAU_POWER, XTAU_POWER, TCOM1,             & ! Input
            T_DIRECT_UP, T_DIRECT_DN, T_UT_DIRECT_UP, T_UT_DIRECT_DN, & ! Input
            Therm )                                                     ! Output
        ENDIF

        IF ( DO_SOLAR_SOURCES ) THEN
          IF (.NOT.DO_SSFULL .OR. (DO_SSFULL.AND.DO_SSCORR_NADIR)) THEN
            IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
              CALL EMULT_MASTER ( &
                DO_UPWELLING, DO_DNWELLING, TAYLOR_ORDER, & ! Add, Version 2.7
                NLAYERS, &
                LAYER_PIS_CUTOFF, NBEAMS, &
                N_USER_VZANGLES, &
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
            ELSE
              CALL EMULT_MASTER_OBSGEO ( &
                DO_UPWELLING, DO_DNWELLING, TAYLOR_ORDER, & ! Add, Version 2.7
                NLAYERS, &
                LAYER_PIS_CUTOFF, NBEAMS, &
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

            CALL VLIDORT_PACK_MULT ( &
              NLAYERS, N_PARTLAYERS, N_SZANGLES, N_USER_VZANGLES, & ! Input
              EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN,       & ! Input
              Mult )                                                ! Output
          ENDIF
        ENDIF

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

        IF ( DO_USER_VZANGLES ) THEN

!  Regular (nadir view) SS correction

          IF ( DO_SSCORR_NADIR .AND. .NOT.DO_FO_CALC ) THEN

            CALL VLIDORT_SSCORR_NADIR ( &
              SS_FLUX_MULTIPLIER, &
              DO_SSCORR_TRUNCATION, DO_REFRACTIVE_GEOMETRY, &
              DO_DELTAM_SCALING, DO_UPWELLING, &
              DO_DNWELLING, DO_OBSERVATION_GEOMETRY, NSTOKES, &
              NLAYERS, NGREEK_MOMENTS_INPUT, &
              N_USER_RELAZMS, USER_RELAZMS, &
              N_USER_LEVELS, &
              OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT, &
              DO_TOA_CONTRIBS, &
              FLUXVEC, COS_SZANGLES, &
              SIN_SZANGLES, SUN_SZA_COSINES, &
              NBEAMS, N_USER_VZANGLES, &
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
                    DO_DNWELLING, DO_OBSERVATION_GEOMETRY, &
                    NSTOKES, NLAYERS, NGREEK_MOMENTS_INPUT, &
                    N_USER_RELAZMS, USER_RELAZMS, &
                    N_USER_LEVELS, GREEKMAT_TOTAL_INPUT, &
                    FLUXVEC, COS_SZANGLES, &
                    SIN_SZANGLES, SUN_SZA_COSINES, &
                    NMOMENTS, NBEAMS, &
                    N_USER_VZANGLES, LAYER_MAXMOMENTS, &
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
                    LP_EMULT_UP, LP_EMULT_DN, &
                    LP_UT_EMULT_UP, LP_UT_EMULT_DN, &
                    PROFILEWF_SS )

                ENDIF
              ENDDO
            ENDIF
          ENDIF

!  Outgoing sphericity correction ------ Added 20 March 2007.
!     Calling not changed for column Jacobians (Version 2.4, December 20

          IF ( DO_SSCORR_OUTGOING .AND. .NOT.DO_FO_CALC ) THEN

           IF ( DO_ATMOS_LINEARIZATION ) THEN

            CALL VLIDORT_LP_SSCORR_OUTGOING ( &
              SS_FLUX_MULTIPLIER, DO_SSCORR_TRUNCATION, DO_DELTAM_SCALING, &
              DO_UPWELLING, DO_DNWELLING, DO_OBSERVATION_GEOMETRY, &
              NSTOKES, NLAYERS, &
              NGREEK_MOMENTS_INPUT, N_USER_RELAZMS, N_USER_LEVELS, &
              EARTH_RADIUS, HEIGHT_GRID, &
              OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT, &
              DO_TOA_CONTRIBS, FLUXVEC, NMOMENTS, &
              NBEAMS, N_USER_VZANGLES, &
              LAYER_MAXMOMENTS, USER_VZANGLES_ADJUST, &
              SZANGLES_ADJUST, USER_RELAZMS_ADJUST, &
              N_PARTLAYERS, NFINELAYERS, &
              PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
              UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, &
              PARTLAYERS_LAYERIDX, &
              STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, &
              N_GEOMETRIES, VZA_OFFSETS, &
              DELTAU_VERT, PARTAU_VERT, TRUNC_FACTOR, &
              DO_PROFILE_LINEARIZATION, &
              LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
              L_OMEGA_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT, &
              L_DELTAU_VERT, DO_SCATMAT_VARIATION, &
              UP_MULTIPLIERS, DN_MULTIPLIERS, UP_LOSTRANS, DN_LOSTRANS, &
              UP_MULTIPLIERS_UT, DN_MULTIPLIERS_UT, &
              UP_LOSTRANS_UT, DN_LOSTRANS_UT, &
              ZMAT_UP, ZMAT_DN, TMS, SSFDEL, &
              SS_CUMSOURCE_UP, SS_CUMSOURCE_DN, &
              BOA_ATTN, STOKES_SS, SS_CONTRIBS, &
              LP_UP_MULTIPLIERS, LP_DN_MULTIPLIERS, &
              L_UP_LOSTRANS, L_DN_LOSTRANS, &
              LP_UP_MULTIPLIERS_UT, LP_DN_MULTIPLIERS_UT, &
              L_UP_LOSTRANS_UT, L_DN_LOSTRANS_UT, &
              LP_BOA_ATTN, PROFILEWF_SS, &
              FAIL, MESSAGE, TRACE_1 )

            IF ( FAIL ) THEN
              TRACE_2 = 'Linearized SS correction outgoing '// &
                           'failed in VLIDORT_LPS_MASTER'
              STATUS_SUB  = VLIDORT_SERIOUS
              RETURN
            ENDIF

           ELSE

            CALL VLIDORT_SSCORR_OUTGOING ( &
              SS_FLUX_MULTIPLIER, &
              DO_SSCORR_TRUNCATION, DO_DELTAM_SCALING, &
              DO_UPWELLING, DO_DNWELLING, DO_OBSERVATION_GEOMETRY, &
              NSTOKES, NLAYERS, &
              NGREEK_MOMENTS_INPUT, &
              N_USER_RELAZMS, N_USER_LEVELS, &
              EARTH_RADIUS, HEIGHT_GRID, &
              OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT, &
              DO_TOA_CONTRIBS, &
              FLUXVEC, NBEAMS, &
              N_USER_VZANGLES, LAYER_MAXMOMENTS, &
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
                         'failed in VLIDORT_LPS_MASTER'
              STATUS_SUB  = VLIDORT_SERIOUS
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

          IF ( DO_SSFULL .AND. .NOT.DO_FO_CALC ) THEN
            CALL VLIDORT_DBCORRECTION ( &
              SS_FLUX_MULTIPLIER, &
              DO_SSCORR_OUTGOING, DO_REFRACTIVE_GEOMETRY, &
              DO_UPWELLING, DO_OBSERVATION_GEOMETRY, &
              NSTOKES, NLAYERS, &
              SZA_LOCAL_INPUT, N_USER_RELAZMS, &
              N_USER_LEVELS, &
              DO_LAMBERTIAN_SURFACE, &
              LAMBERTIAN_ALBEDO, EXACTDB_BRDFUNC, &
              DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, &
              SLTERM_ISOTROPIC, SLTERM_USERANGLES, &
              FLUXVEC, COS_SZANGLES, &
              NBEAMS, N_USER_VZANGLES, &
              MUELLER_INDEX, SZANGLES_ADJUST, &
              PARTLAYERS_OUTFLAG, &
              PARTLAYERS_OUTINDEX, UTAU_LEVEL_MASK_UP, &
              PARTLAYERS_LAYERIDX, N_GEOMETRIES, &
              VZA_OFFSETS, &
              TRANS_SOLAR_BEAM, DO_REFLECTED_DIRECTBEAM, &
              T_DELT_USERM, T_UTUP_USERM, &
              UP_LOSTRANS, UP_LOSTRANS_UT, BOA_ATTN, &
              DB_CUMSOURCE, EXACTDB_SOURCE, &
              ATTN_DB_SAVE, STOKES_DB )

            IF ( DO_PROFILE_LINEARIZATION ) THEN
              CALL VLIDORT_LAP_DBCORRECTION ( &
                SS_FLUX_MULTIPLIER, &
                DO_SSCORR_OUTGOING, DO_REFRACTIVE_GEOMETRY, &
                DO_UPWELLING, DO_OBSERVATION_GEOMETRY, &
                NSTOKES, NLAYERS, &
                SZA_LOCAL_INPUT, N_USER_RELAZMS, &
                N_USER_LEVELS, &
                DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
                FLUXVEC, COS_SZANGLES, NBEAMS, &
                N_USER_VZANGLES, SZANGLES_ADJUST, &
                PARTLAYERS_OUTFLAG, &
                PARTLAYERS_OUTINDEX, UTAU_LEVEL_MASK_UP, &
                PARTLAYERS_LAYERIDX, N_GEOMETRIES, &
                VZA_OFFSETS, TRANS_SOLAR_BEAM, &
                DO_REFLECTED_DIRECTBEAM, T_DELT_USERM, &
                T_UTUP_USERM, &
                LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
                L_DELTAU_SLANT, L_T_DELT_USERM, &
                L_T_UTUP_USERM, &
                DB_CUMSOURCE, EXACTDB_SOURCE, &
                UP_LOSTRANS, UP_LOSTRANS_UT, &
                BOA_ATTN, ATTN_DB_SAVE, &
                L_UP_LOSTRANS, L_UP_LOSTRANS_UT, &
                LP_BOA_ATTN, L_EXACTDB_SOURCE, &
                L_DB_CUMSOURCE, PROFILEWF_DB )
            ENDIF

            IF ( DO_SURFACE_LINEARIZATION ) THEN
              CALL VLIDORT_LS_DBCORRECTION ( &
                SS_FLUX_MULTIPLIER, &
                DO_SSCORR_OUTGOING, DO_UPWELLING, &
                DO_OBSERVATION_GEOMETRY, &
                NSTOKES, NLAYERS, &
                N_USER_RELAZMS, N_USER_LEVELS, &
                DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
                FLUXVEC, NBEAMS, &
                N_USER_VZANGLES, MUELLER_INDEX, &
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

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@ Addition of SLEAVE WFs to Corrected Directbeam @@@@@@@@@@
!    R. Spurr, 22 August 2012
!
            IF ( DO_SURFACE_LEAVING .and. DO_SLEAVE_WFS ) then
              CALL VLIDORT_LSSL_DBCORRECTION ( &
                DO_SSCORR_OUTGOING, DO_UPWELLING, DO_SL_ISOTROPIC,       &
                DO_OBSERVATION_GEOMETRY,                                 &
                NSTOKES, NLAYERS, NBEAMS, N_SLEAVE_WFS, N_SURFACE_WFS,   &
                N_USER_VZANGLES, N_USER_RELAZMS, N_USER_LEVELS,           &
                PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                 &
                PARTLAYERS_LAYERIDX, UTAU_LEVEL_MASK_UP,                 &
                N_GEOMETRIES, VZA_OFFSETS, DO_REFLECTED_DIRECTBEAM,      &
                SS_FLUX_MULTIPLIER, LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_USERANGLES, &
                T_DELT_USERM, T_UTUP_USERM, UP_LOSTRANS, UP_LOSTRANS_UT, &
                SURFACEWF_DB )
            ENDIF

!@@@@@@@@@@ END Addition of SLEAVE Corrected Directbeam @@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

          ELSE
            IF ( DO_DBCORRECTION .AND. .NOT.DO_FO_CALC ) THEN
              CALL VLIDORT_DBCORRECTION ( &
                SS_FLUX_MULTIPLIER, &
                DO_SSCORR_OUTGOING, DO_REFRACTIVE_GEOMETRY, &
                DO_UPWELLING, DO_OBSERVATION_GEOMETRY, &
                NSTOKES, NLAYERS, &
                SZA_LOCAL_INPUT, N_USER_RELAZMS, &
                N_USER_LEVELS, &
                DO_LAMBERTIAN_SURFACE, &
                LAMBERTIAN_ALBEDO, EXACTDB_BRDFUNC, &
                DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, &
                SLTERM_ISOTROPIC, SLTERM_USERANGLES, &
                FLUXVEC, COS_SZANGLES, &
                NBEAMS, N_USER_VZANGLES, &
                MUELLER_INDEX, SZANGLES_ADJUST, &
                PARTLAYERS_OUTFLAG, &
                PARTLAYERS_OUTINDEX, UTAU_LEVEL_MASK_UP, &
                PARTLAYERS_LAYERIDX, N_GEOMETRIES, &
                VZA_OFFSETS, &
                TRANS_SOLAR_BEAM, DO_REFLECTED_DIRECTBEAM, &
                T_DELT_USERM, T_UTUP_USERM, &
                UP_LOSTRANS, UP_LOSTRANS_UT, BOA_ATTN, &
                DB_CUMSOURCE, EXACTDB_SOURCE, &
                ATTN_DB_SAVE, STOKES_DB )

              IF ( DO_PROFILE_LINEARIZATION ) THEN
                CALL VLIDORT_LAP_DBCORRECTION ( &
                  SS_FLUX_MULTIPLIER, &
                  DO_SSCORR_OUTGOING, DO_REFRACTIVE_GEOMETRY, &
                  DO_UPWELLING, DO_OBSERVATION_GEOMETRY, &
                  NSTOKES, NLAYERS, &
                  SZA_LOCAL_INPUT, N_USER_RELAZMS, &
                  N_USER_LEVELS, &
                  DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
                  FLUXVEC, COS_SZANGLES, NBEAMS, &
                  N_USER_VZANGLES, SZANGLES_ADJUST, &
                  PARTLAYERS_OUTFLAG, &
                  PARTLAYERS_OUTINDEX, UTAU_LEVEL_MASK_UP, &
                  PARTLAYERS_LAYERIDX, N_GEOMETRIES, &
                  VZA_OFFSETS, TRANS_SOLAR_BEAM, &
                  DO_REFLECTED_DIRECTBEAM, T_DELT_USERM, &
                  T_UTUP_USERM, &
                  LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
                  L_DELTAU_SLANT, L_T_DELT_USERM, &
                  L_T_UTUP_USERM, &
                  DB_CUMSOURCE, EXACTDB_SOURCE, &
                  UP_LOSTRANS, UP_LOSTRANS_UT, &
                  BOA_ATTN, ATTN_DB_SAVE, &
                  L_UP_LOSTRANS, L_UP_LOSTRANS_UT, &
                  LP_BOA_ATTN, L_EXACTDB_SOURCE, &
                  L_DB_CUMSOURCE, PROFILEWF_DB )
              ENDIF

              IF ( DO_SURFACE_LINEARIZATION ) THEN
                CALL VLIDORT_LS_DBCORRECTION ( &
                  SS_FLUX_MULTIPLIER, &
                  DO_SSCORR_OUTGOING, DO_UPWELLING, &
                  DO_OBSERVATION_GEOMETRY, &
                  NSTOKES, NLAYERS, &
                  N_USER_RELAZMS, N_USER_LEVELS, &
                  DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
                  FLUXVEC, NBEAMS, &
                  N_USER_VZANGLES, MUELLER_INDEX, &
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

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@ Addition of SLEAVE WFs to Corrected Directbeam @@@@@@@@@@
!    R. Spurr, 22 August 2012
!
              IF ( DO_SURFACE_LEAVING .and. DO_SLEAVE_WFS ) then
                CALL VLIDORT_LSSL_DBCORRECTION ( &
                  DO_SSCORR_OUTGOING, DO_UPWELLING, DO_SL_ISOTROPIC,       &
                  DO_OBSERVATION_GEOMETRY,                                 &
                  NSTOKES, NLAYERS, NBEAMS, N_SLEAVE_WFS, N_SURFACE_WFS,   &
                  N_USER_VZANGLES, N_USER_RELAZMS, N_USER_LEVELS,           &
                  PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                 &
                  PARTLAYERS_LAYERIDX, UTAU_LEVEL_MASK_UP,                 &
                  N_GEOMETRIES, VZA_OFFSETS, DO_REFLECTED_DIRECTBEAM,      &
                  SS_FLUX_MULTIPLIER, LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_USERANGLES, &
                  T_DELT_USERM, T_UTUP_USERM, UP_LOSTRANS, UP_LOSTRANS_UT, &
                  SURFACEWF_DB )
              ENDIF

!@@@@@@@@@@ END Addition of SLEAVE Corrected Directbeam @@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

            ENDIF
          ENDIF

! ################################################################
!  New first-order code added 02 Jul 2013

          IF ( DO_FO_CALC ) THEN
            CALL VFO_LPS_MASTER_INTERFACE ( &
              DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,    & ! Input flags
              DO_PLANE_PARALLEL, DO_SSCORR_NADIR, DO_SSCORR_OUTGOING,        & ! Input flags
              DO_DELTAM_SCALING, DO_UPWELLING, DO_DNWELLING,                 & ! Input flags
              DO_LAMBERTIAN_SURFACE, DO_OBSERVATION_GEOMETRY,                & ! Input flags
              DO_PROFILE_LINEARIZATION, DO_SURFACE_LINEARIZATION,            & ! Input flags
              NSTOKES, NLAYERS, NFINELAYERS, NMOMENTS, NSTREAMS,             & ! Input numbers
              NGREEK_MOMENTS_INPUT, N_SURFACE_WFS,                           & ! Input numbers
              NBEAMS, SZANGLES, N_USER_VZANGLES, USER_VZANGLES,               & ! Input geometry
              N_USER_RELAZMS, USER_RELAZMS,                                  & ! Input geometry
              N_USER_LEVELS, USER_LEVELS, EARTH_RADIUS, HEIGHT_GRID,         & ! Input other
              SS_FLUX_MULTIPLIER, FLUXVEC,                                   & ! Inputs (Optical - Atmos)
              DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,    & ! Inputs (Optical - Atmos)
              DELTAU_VERT, OMEGA_TOTAL, GREEKMAT_TOTAL, TRUNC_FACTOR,        & ! Inputs (Optical - Atmos)
              THERMAL_BB_INPUT,                                              & ! Inputs (Optical - Atmos)
              LAMBERTIAN_ALBEDO, EXACTDB_BRDFUNC, SURFBB, USER_EMISSIVITY,   & ! Inputs (Optical - Surf)
              LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                            & ! Inputs (Lin control)
              L_DELTAU_VERT_INPUT, L_OMEGA_TOTAL_INPUT,                      & ! Inputs (Optical - Lin Atmos)
              L_GREEKMAT_TOTAL_INPUT, L_DELTAU_VERT, L_TRUNC_FACTOR,         & ! Inputs (Optical - Lin Atmos)
              LS_EXACTDB_BRDFUNC, LS_USER_EMISSIVITY,                        & ! Inputs (Optical - Lin Surf)
              FO_STOKES_SS, FO_STOKES_DB, FO_STOKES_DTA, FO_STOKES_DTS,      & ! Output - Std
              FO_STOKES,                                                     & ! Output - Std
              FO_PROFILEWF_SS,  FO_PROFILEWF_DB,  FO_SURFACEWF_DB,           & ! Output - Lin
              FO_PROFILEWF_DTA, FO_PROFILEWF_DTS, FO_SURFACEWF_DTS,          & ! Output - Lin
              FO_PROFILEWF, FO_SURFACEWF,                                    & ! Output - Lin
              FAIL, MESSAGE, TRACE_1 )                                         ! Output
            IF ( FAIL ) THEN
              TRACE_2 = 'VFO_LPS_MASTER_INTERFACE failed, VLIDORT_LPS_MASTER'
              STATUS_SUB  = VLIDORT_SERIOUS
              RETURN
            ENDIF

            !Note: set for ObsGeo mode only at present
            STOKES_SS(1:N_USER_LEVELS,1:NBEAMS,1:NSTOKES,:)   = &
              FO_STOKES_SS(1:N_USER_LEVELS,1:NBEAMS,1:NSTOKES,:)
            STOKES_DB(1:N_USER_LEVELS,1:NBEAMS,1:NSTOKES)     = &
              FO_STOKES_DB(1:N_USER_LEVELS,1:NBEAMS,1:NSTOKES)

            PROFILEWF_SS(1:N_TOTALPROFILE_WFS,1:NLAYERS,&
                         1:N_USER_LEVELS,1:NBEAMS,1:NSTOKES,:) = &
              FO_PROFILEWF_SS(1:N_TOTALPROFILE_WFS,1:NLAYERS,&
                         1:N_USER_LEVELS,1:NBEAMS,1:NSTOKES,:)
            PROFILEWF_DB(1:N_TOTALPROFILE_WFS,1:NLAYERS,&
                         1:N_USER_LEVELS,1:NBEAMS,1:NSTOKES) = &
              FO_PROFILEWF_DB(1:N_TOTALPROFILE_WFS,1:NLAYERS,&
                         1:N_USER_LEVELS,1:NBEAMS,1:NSTOKES)
            SURFACEWF_DB(1:N_SURFACE_WFS,1:N_USER_LEVELS,&
                         1:NBEAMS,1:NSTOKES) = &
              FO_SURFACEWF_DB(1:N_SURFACE_WFS,1:N_USER_LEVELS,&
                         1:NBEAMS,1:NSTOKES)

!IF ( DO_SSFULL ) write(*,*) 'DO_SSFULL = .TRUE. -> DONE ...'

          ENDIF
! ################################################################

!  End user stream if block

        ENDIF

!  End solar srcs if block

      ENDIF

!  ####################
!   MAIN FOURIER LOOP
!  ####################

!  Initialise Fourier loop
!  =======================

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

!  Initialise BVP telescoping. Important.

      DO_BVTEL_INITIAL = DO_BVP_TELESCOPING

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
          IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
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
          ELSE
            DO IB = 1, NSOURCES
              !AZM_ARGUMENT = USER_RELAZMS_ADJUST(LUM,IB,LUA) * DFC
              AZM_ARGUMENT = USER_RELAZMS(IB) * DFC
              AZMFAC(LUM,IB,LUA,1)   = DCOS(DEG_TO_RAD*AZM_ARGUMENT)
              AZMFAC(LUM,IB,LUA,2)   = AZMFAC(LUM,IB,LUA,1)
              AZMFAC(LUM,IB,LUA,3)   = DSIN(DEG_TO_RAD*AZM_ARGUMENT)
              AZMFAC(LUM,IB,LUA,4)   = AZMFAC(LUM,IB,LUA,3)
            ENDDO
          ENDIF
        ENDIF

!  Main call to VLidort Fourier module

!        write(*,*)' ..calculating fourier component',FOURIER_COMPONENT

        CALL VLIDORT_LPS_FOURIER ( FOURIER_COMPONENT,  &
          DO_OBSERVATION_GEOMETRY, DO_SSCORR_NADIR, DO_SSFULL, DO_SOLAR_SOURCES, DO_PLANE_PARALLEL,                       & !Input
          DO_REFRACTIVE_GEOMETRY, DO_SOLUTION_SAVING, DO_UPWELLING, DO_DNWELLING,                                         & !Input
          DO_QUAD_OUTPUT, DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY, DO_TOA_CONTRIBS, DO_DIRECT_BEAM, DO_CLASSICAL_SOLUTION,     & !Input
          DO_DBCORRECTION, DO_MULTIBEAM, DO_USER_VZANGLES, DO_PARTLAYERS, DO_MSMODE_VLIDORT, DO_THERMAL_TRANSONLY,        & !Input
          DO_MSMODE_THERMAL, DO_LAMBERTIAN_SURFACE, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION, DO_REAL_EIGENSOLVER,        & !Input
          DO_LAYER_SCATTERING, DO_SPECIALIST_OPTION_2, DO_DEBUG_WRITE, DO_FDTEST,                                         & !Input
          NSTOKES, NSTREAMS, NLAYERS, NBEAMS, N_USER_VZANGLES, N_USER_LEVELS, N_THERMAL_COEFFS, NMOMENTS, NSTREAMS_2,     & !Input
          NTOTAL, N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2, N_PARTLAYERS, N_DIRECTIONS,  TAYLOR_ORDER,          & !Input
          N_ALLLAYERS_UP, N_ALLLAYERS_DN, FLUX_FACTOR, FLUXVEC, COS_SZANGLES, SZA_LOCAL_INPUT, SUN_SZA_COSINES,           & !Input
          USER_STREAMS, USER_SECANTS, QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS, QUAD_HALFWTS,                             & !Input
          MUELLER_INDEX, DMAT, BVP_REGULAR_FLAG, LOCAL_UM_START, WHICH_DIRECTIONS,                                        & !Input
          UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,           & !Input
          STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, DFLUX, DELTAU_VERT_INPUT, THERMAL_BB_INPUT, SURFBB,                     & !Input
          DO_SIMULATION_ONLY, DO_ATMOS_LINEARIZATION, DO_PROFILE_LINEARIZATION, DO_SURFACE_LINEARIZATION,                 & !Input
          DO_ATMOS_LBBF, DO_SURFACE_LBBF, NSTOKES_SQ, N_TOTALPROFILE_WFS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER,             & !Input
          EMISSIVITY, USER_EMISSIVITY, LS_USER_EMISSIVITY, LS_EMISSIVITY,                                                 & !Input
          Misc, Therm, Mult, LAP_Misc, L_Therm, LP_Mult,                                                                  & !Input
          LAMBERTIAN_ALBEDO, BRDF_F, BRDF_F_0, USER_BRDF_F, USER_BRDF_F_0,                                                & !Input
          N_SURFACE_WFS, LS_BRDF_F, LS_BRDF_F_0, LS_USER_BRDF_F, LS_USER_BRDF_F_0,                                        & !Input
          DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, SLTERM_ISOTROPIC, SLTERM_F_0, USER_SLTERM_F_0,                             & !Input
          DO_SLEAVE_WFS, N_SLEAVE_WFS, LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0,LSSL_USER_SLTERM_F_0,                       & !Input
          PIMM_11, PIMM_KM, DO_BVTEL_INITIAL, BVTEL_FOURIER_COMPONENT,                                                    & !InOut
          MEAN_STOKES, FLUX_STOKES, MEAN_DIRECT, FLUX_DIRECT, MINT_PROFILEWF, FLUX_PROFILEWF,                             & !InOut
          MINT_PROFILEWF_DIRECT, FLUX_PROFILEWF_DIRECT, MINT_SURFACEWF, FLUX_SURFACEWF,                                   & !InOut
          ABBWFS_JACOBIANS, ABBWFS_FLUXES, SBBWFS_JACOBIANS, SBBWFS_FLUXES,                                               & !InOut
          MESSAGE, TRACE_1, TRACE_2,                                                                                      & !InOut
          DO_INCLUDE_THERMEMISS, DO_INCLUDE_SURFACE, DO_INCLUDE_SURFEMISS,                                                & !Output
          MS_CONTRIBS_F, STOKES_F, PROFILEWF_F, SURFACEWF_F, STATUS_SUB )                                                   !Output

!  error handling

        IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
          TRACE_3 = ' Called by VLIDORT_LPS_MASTER '
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

            IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
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

!  Fourier summation of linearization quantities

!mick fix 9/6/2012 - added N_SLEAVE_WFS to call
              IF ( DO_LINEARIZATION ) THEN
                CALL VLIDORT_LPS_CONVERGE ( &
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
                  DO_PROFILE_LINEARIZATION, &
                  LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
                  DO_SURFACE_LINEARIZATION, &
                  N_SURFACE_WFS, N_SLEAVE_WFS, &
                  PROFILEWF_SS, PROFILEWF_DB, PROFILEWF_F, &
                  SURFACEWF_DB, SURFACEWF_F, &
                  PROFILEWF, SURFACEWF )
              ENDIF

            ELSE

              CALL VLIDORT_CONVERGE_OBSGEO ( &
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
                WHICH_DIRECTIONS, VZA_OFFSETS, &
                STOKES_SS, STOKES_DB, SS_CONTRIBS, &
                STOKES_F, MS_CONTRIBS_F, IBEAM, &
                STOKES, FOURIER_SAVED, CONTRIBS, &
                BEAM_TESTCONV(IBEAM), BEAM_ITERATION(IBEAM) )

!  Fourier summation of linearization quantities

!mick fix 9/6/2012 - added N_SLEAVE_WFS to call
              IF ( DO_LINEARIZATION ) THEN
                CALL VLIDORT_LPS_CONVERGE_OBSGEO ( &
                  FOURIER_COMPONENT, LOCAL_N_USERAZM, &
                  IBEAM, AZMFAC, &
                  DO_SSCORR_NADIR, DO_SSCORR_OUTGOING, &
                  DO_SS_EXTERNAL, & ! New 15 March 2012
                  DO_SSFULL, DO_UPWELLING, &
                  NSTOKES, NLAYERS, &
                  N_USER_LEVELS, &
                  DO_LAMBERTIAN_SURFACE, &
                  DO_DBCORRECTION, DO_NO_AZIMUTH, &
                  N_DIRECTIONS, WHICH_DIRECTIONS, &
                  VZA_OFFSETS, &
                  DO_PROFILE_LINEARIZATION, &
                  LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
                  DO_SURFACE_LINEARIZATION, &
                  N_SURFACE_WFS, N_SLEAVE_WFS, &
                  PROFILEWF_SS, PROFILEWF_DB, PROFILEWF_F, &
                  SURFACEWF_DB, SURFACEWF_F, &
                  PROFILEWF, SURFACEWF )
              ENDIF

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

        !IF ( DO_LINEARIZATION ) THEN
        !  CALL VLIDORT_L_WRITERESULTS ( &
        !    RUNIT, &
        !    DO_FULLRAD_MODE, DO_SSCORR_NADIR, &
        !    DO_SSCORR_OUTGOING, DO_DOUBLE_CONVTEST, &
        !    DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, &
        !    DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY, &
        !    NSTOKES, NLAYERS, &
        !    VLIDORT_ACCURACY, SZANGLES, &
        !    N_USER_RELAZMS, USER_RELAZMS, &
        !    N_USER_LEVELS, USER_LEVELS, &
        !    DO_LAMBERTIAN_SURFACE, &
        !    DO_NO_AZIMUTH, NBEAMS, &
        !    N_DIRECTIONS, WHICH_DIRECTIONS, &
        !    N_OUT_STREAMS, &
        !    OUT_ANGLES, VZA_OFFSETS, &
        !    DO_MULTIBEAM, FOURIER_SAVED, &
        !    DO_PROFILE_LINEARIZATION, DO_COLUMN_LINEARIZATION, &
        !    N_TOTALCOLUMN_WFS, LAYER_VARY_FLAG, &
        !    LAYER_VARY_NUMBER, DO_SURFACE_LINEARIZATION, &
        !    DO_SURFBB_LINEARIZATION, N_SURFACE_WFS, &
        !    PROFILEWF_NAMES, COLUMNWF_NAMES, &
        !    SURFACEWF, PROFILEWF, &
        !    COLUMNWF, MINT_SURFACEWF, &
        !    MINT_PROFILEWF, FLUX_SURFACEWF, &
        !    FLUX_PROFILEWF )
        !ENDIF

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

!  Modified Boolean inputs

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

      VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY= DO_OBSERVATION_GEOMETRY

!  Modified control inputs

      VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT   = NGREEK_MOMENTS_INPUT

!  Modified beam inputs

      VLIDORT_ModIn%MSunrays%TS_N_SZANGLES             = N_SZANGLES
      VLIDORT_ModIn%MSunRays%TS_SZANGLES(1:N_SZANGLES) = SZANGLES(1:N_SZANGLES)

!  Modified user value inputs

      VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS      = N_USER_RELAZMS
      VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(1:N_USER_RELAZMS)           = &
        USER_RELAZMS(1:N_USER_RELAZMS)

      VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES     = N_USER_VZANGLES
      VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(1:N_USER_VZANGLES)   = &
        USER_VZANGLES(1:N_USER_VZANGLES)

      VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(1:N_USER_LEVELS)             = &
        USER_LEVELS(1:N_USER_LEVELS)

      VLIDORT_ModIn%MUserVal%TS_GEOMETRY_SPECHEIGHT = GEOMETRY_SPECHEIGHT

      VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS     = N_USER_OBSGEOMS
      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1:N_USER_OBSGEOMS,:) = &
        USER_OBSGEOMS(1:N_USER_OBSGEOMS,:)

!  Modified Chapman function inputs

      VLIDORT_ModIn%MChapman%TS_EARTH_RADIUS        = EARTH_RADIUS

!  Modified optical inputs

      VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(1:NLAYERS) = &
        OMEGA_TOTAL_INPUT(1:NLAYERS)

!  Modified linearized control variables

!      First three are additional variables, Version 2.7

      VLIDORT_LinModIn%MCont%TS_DO_SIMULATION_ONLY = DO_SIMULATION_ONLY
      VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LBBF      = DO_ATMOS_LBBF
      VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LBBF    = DO_SURFACE_LBBF

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

      VLIDORT_Out%Main%TS_STOKES(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,:)    = &
        STOKES(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,:)
      VLIDORT_Out%Main%TS_MEAN_STOKES(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:) = &
        MEAN_STOKES(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:)
      VLIDORT_Out%Main%TS_FLUX_STOKES(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:) = &
        FLUX_STOKES(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:)

      VLIDORT_Out%Main%TS_MEAN_DIRECT(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES) = &
        MEAN_DIRECT(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES)
      VLIDORT_Out%Main%TS_FLUX_DIRECT(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES) = &
        FLUX_DIRECT(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES)

!  new 12 March 2012
!   IF SS results already available, no need to copy them !

      IF ( .NOT. DO_SS_EXTERNAL ) THEN
         VLIDORT_Sup%SS%TS_STOKES_SS(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,:) = &
           STOKES_SS(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,:)
         VLIDORT_Sup%SS%TS_STOKES_DB(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES) = &
           STOKES_DB(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES)
      ENDIF

!  Bookkeeping

      VLIDORT_Out%Main%TS_FOURIER_SAVED(1:N_SZANGLES) = &
        FOURIER_SAVED(1:N_SZANGLES)
      VLIDORT_Out%Main%TS_N_GEOMETRIES = &
        N_GEOMETRIES
      VLIDORT_Out%Main%TS_SZA_OFFSETS(1:N_SZANGLES)   = &
        SZA_OFFSETS(1:N_SZANGLES)
      VLIDORT_Out%Main%TS_VZA_OFFSETS(1:N_SZANGLES,1:N_USER_VZANGLES) = &
        VZA_OFFSETS(1:N_SZANGLES,1:N_USER_VZANGLES)

!  Solar Beam Transmittance to BOA
!  rob fix 11/17/2014, for diagnostic use only

      VLIDORT_Out%Main%TS_SOLARBEAM_BOATRANS (1:NBEAMS) = SOLARBEAM_BOATRANS (1:NBEAMS)

!  Atmosphere weighting functions
!    Output direct Mean/Flux added 17 May 2012

      VLIDORT_LinOut%Prof%TS_PROFILEWF&
          (1:N_TOTALPROFILE_WFS,1:NLAYERS,1:N_USER_LEVELS,&
           1:N_GEOMETRIES,1:NSTOKES,:) = &
        PROFILEWF&
          (1:N_TOTALPROFILE_WFS,1:NLAYERS,1:N_USER_LEVELS,&
           1:N_GEOMETRIES,1:NSTOKES,:)

      VLIDORT_LinOut%Prof%TS_MINT_PROFILEWF&
          (1:N_TOTALPROFILE_WFS,1:NLAYERS,1:N_USER_LEVELS,&
           1:N_SZANGLES,1:NSTOKES,:) = &
        MINT_PROFILEWF&
          (1:N_TOTALPROFILE_WFS,1:NLAYERS,1:N_USER_LEVELS,&
           1:N_SZANGLES,1:NSTOKES,:)
      VLIDORT_LinOut%Prof%TS_FLUX_PROFILEWF&
          (1:N_TOTALPROFILE_WFS,1:NLAYERS,1:N_USER_LEVELS,&
           1:N_SZANGLES,1:NSTOKES,:) = &
        FLUX_PROFILEWF&
          (1:N_TOTALPROFILE_WFS,1:NLAYERS,1:N_USER_LEVELS,&
           1:N_SZANGLES,1:NSTOKES,:)

      VLIDORT_LinOut%Prof%TS_MINT_PROFILEWF_DIRECT&
          (1:N_TOTALPROFILE_WFS,1:NLAYERS,1:N_USER_LEVELS,&
           1:N_SZANGLES,1:NSTOKES) = &
        MINT_PROFILEWF_DIRECT&
          (1:N_TOTALPROFILE_WFS,1:NLAYERS,1:N_USER_LEVELS,&
           1:N_SZANGLES,1:NSTOKES)
      VLIDORT_LinOut%Prof%TS_FLUX_PROFILEWF_DIRECT&
          (1:N_TOTALPROFILE_WFS,1:NLAYERS,1:N_USER_LEVELS,&
           1:N_SZANGLES,1:NSTOKES) = &
        FLUX_PROFILEWF_DIRECT&
          (1:N_TOTALPROFILE_WFS,1:NLAYERS,1:N_USER_LEVELS,&
           1:N_SZANGLES,1:NSTOKES)

!  Superseded, Version 2.7
!      VLIDORT_LinOut%Atmos%TS_LTE_ATMOSWF&
!          (0:NLAYERS,1:N_USER_LEVELS,1:N_USER_VZANGLES,:) = &
!        LTE_ATMOSWF&
!          (0:NLAYERS,1:N_USER_LEVELS,1:N_USER_VZANGLES,:)

!  Surface weighting functions

      VLIDORT_LinOut%Surf%TS_SURFACEWF&
          (1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,:) = &
        SURFACEWF&
          (1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,:)
      VLIDORT_LinOut%Surf%TS_MINT_SURFACEWF&
          (1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:) = &
        MINT_SURFACEWF&
          (1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:)
      VLIDORT_LinOut%Surf%TS_FLUX_SURFACEWF&
          (1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:) = &
        FLUX_SURFACEWF&
          (1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:)

      IF ( .NOT. DO_SS_EXTERNAL ) THEN
         !SS atmosphere weighting functions

         VLIDORT_LinSup%SS%Prof%TS_PROFILEWF_SS&
             (1:N_TOTALPROFILE_WFS,1:NLAYERS,1:N_USER_LEVELS,&
              1:N_GEOMETRIES,1:NSTOKES,:) = &
           PROFILEWF_SS&
             (1:N_TOTALPROFILE_WFS,1:NLAYERS,1:N_USER_LEVELS,&
              1:N_GEOMETRIES,1:NSTOKES,:)
         VLIDORT_LinSup%SS%Prof%TS_PROFILEWF_DB&
             (1:N_TOTALPROFILE_WFS,1:NLAYERS,1:N_USER_LEVELS,&
              1:N_GEOMETRIES,1:NSTOKES) = &
           PROFILEWF_DB&
             (1:N_TOTALPROFILE_WFS,1:NLAYERS,1:N_USER_LEVELS,&
              1:N_GEOMETRIES,1:NSTOKES)

         !SS surface weighting functions

         !VLIDORT_LinSup%SS%Surf%TS_SURFACEWF_SS   = SURFACEWF_SS
         VLIDORT_LinSup%SS%Surf%TS_SURFACEWF_DB&
             (1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES) = &
           SURFACEWF_DB&
             (1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES)
      ENDIF

!  BLACKBODY JACOBIANS, New 28 March 2014. 

      VLIDORT_LinOut%Atmos%TS_ABBWFS_JACOBIANS(1:N_USER_LEVELS,1:N_USER_VZANGLES,0:NLAYERS,1:NSTOKES,1:2) = &
          ABBWFS_JACOBIANS(1:N_USER_LEVELS,1:N_USER_VZANGLES,0:NLAYERS,1:NSTOKES,1:2)
      VLIDORT_LinOut%Atmos%TS_ABBWFS_FLUXES(1:N_USER_LEVELS,1:2,0:NLAYERS,1:NSTOKES,1:2) = &
          ABBWFS_FLUXES(1:N_USER_LEVELS,1:2,0:NLAYERS,1:NSTOKES,1:2)

      VLIDORT_LinOut%Surf%TS_SBBWFS_JACOBIANS(1:N_USER_LEVELS,1:N_USER_VZANGLES,1:NSTOKES,1:2) = &
          SBBWFS_JACOBIANS(1:N_USER_LEVELS,1:N_USER_VZANGLES,1:NSTOKES,1:2)
      VLIDORT_LinOut%Surf%TS_SBBWFS_FLUXES(1:N_USER_LEVELS,1:2,1:NSTOKES,1:2) = &
          SBBWFS_FLUXES(1:N_USER_LEVELS,1:2,1:NSTOKES,1:2)

!  Exception handling

      VLIDORT_Out%Status%TS_STATUS_INPUTCHECK  = STATUS_INPUTCHECK
      VLIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION

      VLIDORT_Out%Status%TS_NCHECKMESSAGES = &
        NCHECKMESSAGES
      VLIDORT_Out%Status%TS_CHECKMESSAGES(0:NCHECKMESSAGES)  = &
        CHECKMESSAGES(0:NCHECKMESSAGES)
      VLIDORT_Out%Status%TS_ACTIONS(0:NCHECKMESSAGES)        = &
        ACTIONS(0:NCHECKMESSAGES)

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

!  3/28/14. Changes for Version 2.7. remove LTE linearization references

      CALL VLIDORT_WRITE_STD_INPUT ( &
        DO_FULLRAD_MODE,DO_SSCORR_TRUNCATION,DO_SS_EXTERNAL,DO_SSFULL,&
        DO_THERMAL_EMISSION,DO_SURFACE_EMISSION,DO_PLANE_PARALLEL,&
        DO_UPWELLING,DO_DNWELLING,DO_QUAD_OUTPUT,&
        DO_TOA_CONTRIBS,DO_LAMBERTIAN_SURFACE,&
        DO_SPECIALIST_OPTION_1,DO_SPECIALIST_OPTION_2,DO_SPECIALIST_OPTION_3,&
        DO_SURFACE_LEAVING,DO_SL_ISOTROPIC,&
        NSTOKES,NSTREAMS,NLAYERS,&
        NFINELAYERS,N_THERMAL_COEFFS,VLIDORT_ACCURACY,&
        NLAYERS_NOMS,NLAYERS_CUTOFF,FLUX_FACTOR,N_USER_LEVELS,&
        HEIGHT_GRID,PRESSURE_GRID,TEMPERATURE_GRID,&
        FINEGRID,RFINDEX_PARAMETER,&
        DELTAU_VERT_INPUT,GREEKMAT_TOTAL_INPUT,THERMAL_BB_INPUT,&
        LAMBERTIAN_ALBEDO,SURFBB,& ! *****
        DO_DEBUG_WRITE,DO_WRITE_INPUT,DO_WRITE_SCENARIO,&
        DO_WRITE_FOURIER,DO_WRITE_RESULTS,&
        INPUT_WRITE_FILENAME,SCENARIO_WRITE_FILENAME,&
        FOURIER_WRITE_FILENAME,RESULTS_WRITE_FILENAME,&
        DO_SSCORR_NADIR,DO_SSCORR_OUTGOING,DO_FO_CALC,DO_DOUBLE_CONVTEST,&
        DO_SOLAR_SOURCES,DO_REFRACTIVE_GEOMETRY,DO_CHAPMAN_FUNCTION,&
        DO_RAYLEIGH_ONLY,DO_DELTAM_SCALING,DO_SOLUTION_SAVING,&
        DO_BVP_TELESCOPING,DO_USER_VZANGLES,DO_ADDITIONAL_MVOUT,&
        DO_MVOUT_ONLY,DO_THERMAL_TRANSONLY,DO_OBSERVATION_GEOMETRY,&
        NGREEK_MOMENTS_INPUT,N_SZANGLES,SZANGLES,&
        N_USER_RELAZMS,USER_RELAZMS,N_USER_VZANGLES,USER_VZANGLES,&
        USER_LEVELS,GEOMETRY_SPECHEIGHT,N_USER_OBSGEOMS,USER_OBSGEOMS,&
        EARTH_RADIUS,OMEGA_TOTAL_INPUT)

      IF (.NOT. DO_LAMBERTIAN_SURFACE) THEN
        CALL VLIDORT_WRITE_SUP_BRDF_INPUT ( &
          NSTOKES,NSTREAMS,N_SZANGLES,N_USER_VZANGLES,N_USER_RELAZMS,&
          EXACTDB_BRDFUNC,BRDF_F_0,BRDF_F,USER_BRDF_F_0,USER_BRDF_F,&
          EMISSIVITY,USER_EMISSIVITY)
      END IF

      IF (DO_SS_EXTERNAL) THEN
        CALL VLIDORT_WRITE_SUP_SS_INPUT ( &
          NSTOKES,N_USER_LEVELS,&
          STOKES_SS,STOKES_DB)
      END IF

      IF (DO_SURFACE_LEAVING) THEN
        CALL VLIDORT_WRITE_SUP_SLEAVE_INPUT ( &
          NSTOKES,NSTREAMS,N_SZANGLES,N_USER_VZANGLES,N_USER_RELAZMS,&
          SLTERM_ISOTROPIC,SLTERM_USERANGLES,SLTERM_F_0,USER_SLTERM_F_0)
      END IF

      END SUBROUTINE VLIDORT_DEBUG_INPUT_MASTER

      SUBROUTINE VLIDORT_DEBUG_LIN_INPUT_MASTER()

!  3/28/14. Changes for Version 2.7. remove LTE linearization references. Add LBBF

      CALL VLIDORT_WRITE_LIN_INPUT ( &
        NSTOKES,NLAYERS,NGREEK_MOMENTS_INPUT,&
        DO_SIMULATION_ONLY,&
        LAYER_VARY_FLAG,LAYER_VARY_NUMBER,&
        N_TOTALCOLUMN_WFS,N_TOTALPROFILE_WFS,N_SURFACE_WFS,N_SLEAVE_WFS,&
        COLUMNWF_NAMES,PROFILEWF_NAMES,DO_ATMOS_LBBF,DO_SURFACE_LBBF,& !**********
        L_DELTAU_VERT_INPUT,L_OMEGA_TOTAL_INPUT,L_GREEKMAT_TOTAL_INPUT,&
        DO_COLUMN_LINEARIZATION,DO_PROFILE_LINEARIZATION,DO_ATMOS_LINEARIZATION,&
        DO_SURFACE_LINEARIZATION,DO_LINEARIZATION,DO_SLEAVE_WFS)

      IF (.NOT. DO_LAMBERTIAN_SURFACE .AND. DO_SURFACE_LINEARIZATION) THEN
        CALL VLIDORT_WRITE_LIN_SUP_BRDF_INPUT ( &
          NSTOKES,NSTREAMS,N_SZANGLES,N_USER_VZANGLES,&
          N_USER_RELAZMS,N_SURFACE_WFS,&
          LS_EXACTDB_BRDFUNC,LS_BRDF_F_0,LS_BRDF_F,&
          LS_USER_BRDF_F_0,LS_USER_BRDF_F,&
          LS_EMISSIVITY,LS_USER_EMISSIVITY)
      END IF

      IF (DO_SS_EXTERNAL) THEN
        COLUMNWF_SS(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,:,1:NSTOKES,:) = ZERO
        COLUMNWF_DB(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,:,1:NSTOKES)   = ZERO

        CALL VLIDORT_WRITE_LIN_SUP_SS_INPUT ( &
          DO_COLUMN_LINEARIZATION,DO_PROFILE_LINEARIZATION,&
          DO_SURFACE_LINEARIZATION,&
          NSTOKES,NLAYERS,N_USER_LEVELS,&
          N_TOTALCOLUMN_WFS,N_TOTALPROFILE_WFS,N_TOTALSURFACE_WFS,&
          COLUMNWF_SS,COLUMNWF_DB,PROFILEWF_SS,PROFILEWF_DB,SURFACEWF_DB)
      END IF

      IF (DO_SURFACE_LEAVING .AND. DO_SURFACE_LINEARIZATION &
          .AND. DO_SLEAVE_WFS) THEN
        CALL VLIDORT_WRITE_LIN_SUP_SLEAVE_INPUT ( &
          NSTOKES,NSTREAMS,N_SZANGLES,N_USER_VZANGLES,&
          N_USER_RELAZMS,N_SLEAVE_WFS,&
          LSSL_SLTERM_ISOTROPIC,LSSL_SLTERM_USERANGLES,&
          LSSL_SLTERM_F_0,LSSL_USER_SLTERM_F_0)
      END IF

      END SUBROUTINE VLIDORT_DEBUG_LIN_INPUT_MASTER

      END SUBROUTINE VLIDORT_LPS_MASTER

!

      SUBROUTINE VLIDORT_LPS_FOURIER ( FOURIER_COMPONENT,  &
        DO_OBSERVATION_GEOMETRY, DO_SSCORR_NADIR, DO_SSFULL, DO_SOLAR_SOURCES, DO_PLANE_PARALLEL,                       & !Input
        DO_REFRACTIVE_GEOMETRY, DO_SOLUTION_SAVING, DO_UPWELLING, DO_DNWELLING,                                         & !Input
        DO_QUAD_OUTPUT, DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY, DO_TOA_CONTRIBS, DO_DIRECT_BEAM, DO_CLASSICAL_SOLUTION,     & !Input
        DO_DBCORRECTION, DO_MULTIBEAM, DO_USER_STREAMS, DO_PARTLAYERS, DO_MSMODE_VLIDORT, DO_THERMAL_TRANSONLY,         & !Input
        DO_MSMODE_THERMAL, DO_LAMBERTIAN_SURFACE, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION, DO_REAL_EIGENSOLVER,        & !Input
        DO_LAYER_SCATTERING, DO_SPECIALIST_OPTION_2, DO_DEBUG_WRITE, DO_FDTEST,                                         & !Input
        NSTOKES, NSTREAMS, NLAYERS, NBEAMS, N_USER_STREAMS, N_USER_LEVELS, N_THERMAL_COEFFS, NMOMENTS, NSTREAMS_2,      & !Input
        NTOTAL, N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2, N_PARTLAYERS, N_DIRECTIONS, TAYLOR_ORDER,           & !Input
        N_ALLLAYERS_UP, N_ALLLAYERS_DN, FLUX_FACTOR, FLUXVEC, COS_SZANGLES, SZA_LOCAL_INPUT, SUN_SZA_COSINES,           & !Input
        USER_STREAMS, USER_SECANTS, QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS, QUAD_HALFWTS,                             & !Input
        MUELLER_INDEX, DMAT, BVP_REGULAR_FLAG, LOCAL_UM_START, WHICH_DIRECTIONS,                                        & !Input
        UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,           & !Input
        STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, DFLUX, DELTAU_VERT_INPUT, THERMAL_BB_INPUT, SURFBB,                     & !Input
        DO_SIMULATION_ONLY, DO_ATMOS_LINEARIZATION, DO_PROFILE_LINEARIZATION, DO_SURFACE_LINEARIZATION,                 & !Input
        DO_ATMOS_LBBF, DO_SURFACE_LBBF, NSTOKES_SQ, N_TOTALPROFILE_WFS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER,             & !Input
        EMISSIVITY, USER_EMISSIVITY, LS_USER_EMISSIVITY, LS_EMISSIVITY,                                                 & !Input
        Misc, Therm, Mult, LAP_Misc, L_Therm, LP_Mult,                                                                  & !Input
        LAMBERTIAN_ALBEDO, BRDF_F, BRDF_F_0, USER_BRDF_F, USER_BRDF_F_0,                                                & !Input
        N_SURFACE_WFS, LS_BRDF_F, LS_BRDF_F_0, LS_USER_BRDF_F, LS_USER_BRDF_F_0,                                        & !Input
        DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, SLTERM_ISOTROPIC, SLTERM_F_0, USER_SLTERM_F_0,                             & !Input
        DO_SLEAVE_WFS, N_SLEAVE_WFS, LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0,LSSL_USER_SLTERM_F_0,                       & !Input
        PIMM_11, PIMM_KM, DO_BVTEL_INITIAL, BVTEL_FOURIER_COMPONENT,                                                    & !InOut
        MEAN_STOKES, FLUX_STOKES, MEAN_DIRECT, FLUX_DIRECT, MINT_PROFILEWF, FLUX_PROFILEWF,                             & !InOut
        MINT_PROFILEWF_DIRECT, FLUX_PROFILEWF_DIRECT, MINT_SURFACEWF, FLUX_SURFACEWF,                                   & !InOut
        ABBWFS_JACOBIANS, ABBWFS_FLUXES, SBBWFS_JACOBIANS, SBBWFS_FLUXES,                                               & !InOut
        MESSAGE, TRACE_1, TRACE_2,                                                                                      & !InOut
        DO_INCLUDE_THERMEMISS, DO_INCLUDE_SURFACE, DO_INCLUDE_SURFEMISS,                                                & !Output
        MS_CONTRIBS_F, STOKES_F, PROFILEWF_F, SURFACEWF_F, STATUS )                                                       !Output

!  Complete Fourier component calculation for the Extended Code
!    Stokes vector and Weighting function computations

      USE VLIDORT_PARS

      USE VLIDORT_Work_def
      USE VLIDORT_LinWork_def

      USE VLIDORT_MISCSETUPS_MODULE
      USE VLIDORT_THERMALSUP
      USE VLIDORT_MULTIPLIERS
      USE VLIDORT_SOLUTIONS
      USE VLIDORT_BVPROBLEM
      USE VLIDORT_INTENSITY

      USE VLIDORT_LA_MISCSETUPS
      USE VLIDORT_LP_MISCSETUPS
      USE VLIDORT_L_THERMALSUP
!      USE VLIDORT_LA_CORRECTIONS
!      USE VLIDORT_LP_CORRECTIONS
!      USE VLIDORT_LS_CORRECTIONS
      USE VLIDORT_LPC_SOLUTIONS
      USE VLIDORT_LP_SOLUTIONS
      USE VLIDORT_LPC_BVPROBLEM
      USE VLIDORT_LP_BVPROBLEM
      USE VLIDORT_LP_WFATMOS
      USE VLIDORT_LS_WFSURFACE

      USE VLIDORT_UNPACK
      USE VLIDORT_L_UNPACK
      USE VLIDORT_LP_UNPACK

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@ Addition of SLEAVE WF stuff, R. Spurr, 22 August 2012 @@@@@@@@@
      USE VLIDORT_LS_WFSLEAVE
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Addition of LBBF Jacobians, Version 2.7, 3/28/14

      use vlidort_lbbf_jacobians_m

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          FOURIER_COMPONENT

      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY
      LOGICAL, INTENT (IN) ::          DO_SSCORR_NADIR
      !LOGICAL, INTENT (IN) ::          DO_SSCORR_OUTGOING
      !LOGICAL, INTENT (IN) ::          DO_FO_CALC !New 02 Jul 2013
      !LOGICAL, INTENT (IN) ::          DO_SSCORR_TRUNCATION
      LOGICAL, INTENT (IN) ::          DO_SSFULL
      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL
      LOGICAL, INTENT (IN) ::          DO_REFRACTIVE_GEOMETRY
      !LOGICAL, INTENT (IN) ::          DO_DELTAM_SCALING
      LOGICAL, INTENT (IN) ::          DO_SOLUTION_SAVING
      !LOGICAL, INTENT (IN) ::          DO_BVP_TELESCOPING
      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_DNWELLING
      LOGICAL, INTENT (IN) ::          DO_QUAD_OUTPUT
      LOGICAL, INTENT (IN) ::          DO_ADDITIONAL_MVOUT
      LOGICAL, INTENT (IN) ::          DO_MVOUT_ONLY
      LOGICAL, INTENT (IN) ::          DO_TOA_CONTRIBS
      LOGICAL, INTENT (IN) ::          DO_DIRECT_BEAM
      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION
      LOGICAL, INTENT (IN) ::          DO_DBCORRECTION
      LOGICAL, INTENT (IN) ::          DO_MULTIBEAM &
          ( MAXBEAMS, 0:MAXFOURIER )
      LOGICAL, INTENT (IN) ::          DO_USER_STREAMS
      LOGICAL, INTENT (IN) ::          DO_PARTLAYERS
      LOGICAL, INTENT (IN) ::          DO_MSMODE_VLIDORT
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_MSMODE_THERMAL
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT (IN) ::          DO_THERMAL_EMISSION
      LOGICAL, INTENT (IN) ::          DO_SURFACE_EMISSION
      LOGICAL, INTENT (IN) ::          DO_REAL_EIGENSOLVER &
          ( 0:MAXMOMENTS, MAXLAYERS )
      LOGICAL, INTENT (IN) ::          DO_LAYER_SCATTERING &
          ( 0:MAXMOMENTS, MAXLAYERS )
      LOGICAL, INTENT (IN) ::          DO_SPECIALIST_OPTION_2
      !LOGICAL, INTENT (IN) ::          DO_SPECIALIST_OPTION_3
      LOGICAL, INTENT (IN) ::          DO_DEBUG_WRITE
      LOGICAL, INTENT (IN) ::          DO_FDTEST

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTOKES_SQ
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      !INTEGER, INTENT (IN) ::          NGREEK_MOMENTS_INPUT
      INTEGER, INTENT (IN) ::          NBEAMS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      !INTEGER, INTENT (IN) ::          N_USER_RELAZMS
      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      INTEGER, INTENT (IN) ::          N_THERMAL_COEFFS
      !INTEGER, INTENT (IN) ::          NFINELAYERS
      INTEGER, INTENT (IN) ::          NMOMENTS
      INTEGER, INTENT (IN) ::          NSTREAMS_2
      INTEGER, INTENT (IN) ::          NTOTAL
      !INTEGER, INTENT (IN) ::          N_GEOMETRIES
      INTEGER, INTENT (IN) ::          N_SUBDIAG
      INTEGER, INTENT (IN) ::          N_SUPDIAG
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS_2
      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      INTEGER, INTENT (IN) ::          N_DIRECTIONS

!  Version 2p7 input, 2/19/14
      INTEGER, INTENT (IN) ::          TAYLOR_ORDER

      !INTEGER, INTENT (IN) ::          NLAYERS_CUTOFF
      INTEGER, INTENT (IN) ::          N_ALLLAYERS_UP
      INTEGER, INTENT (IN) ::          N_ALLLAYERS_DN

      DOUBLE PRECISION, INTENT (IN) :: FLUX_FACTOR
      DOUBLE PRECISION, INTENT (IN) :: FLUXVEC ( MAXSTOKES )
      !DOUBLE PRECISION, INTENT (IN) :: SS_FLUX_MULTIPLIER

      !DOUBLE PRECISION, INTENT (IN) :: SZANGLES ( MAX_SZANGLES )
      !DOUBLE PRECISION, INTENT (IN) :: SZANGLES_ADJUST &
      !    ( MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS )
      DOUBLE PRECISION, INTENT (IN) :: COS_SZANGLES ( MAX_SZANGLES )
      !DOUBLE PRECISION, INTENT (IN) :: SINBEAMS ( MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) :: SZA_LOCAL_INPUT &
          ( 0:MAXLAYERS, MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) :: SUN_SZA_COSINES &
          ( MAXLAYERS, MAX_SZANGLES )

      !DOUBLE PRECISION, INTENT (IN) :: USER_VZANGLES ( MAX_USER_VZANGLES )
      !DOUBLE PRECISION, INTENT (IN) :: USER_VZANGLES_ADJUST ( MAX_USER_VZANGLES )
      DOUBLE PRECISION, INTENT (IN) :: USER_STREAMS  ( MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: USER_SECANTS  ( MAX_USER_STREAMS )
      !INTEGER, INTENT (IN) ::          VZA_OFFSETS &
      !    ( MAX_SZANGLES, MAX_USER_VZANGLES )

      !DOUBLE PRECISION, INTENT (IN) :: USER_RELAZMS  ( MAX_USER_RELAZMS )
      !DOUBLE PRECISION, INTENT (IN) :: USER_RELAZMS_ADJUST &
      !    ( MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS )
      !DOUBLE PRECISION, INTENT (IN) :: USER_LEVELS   ( MAX_USER_LEVELS )

      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_HALFWTS ( MAXSTREAMS )

      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: DMAT ( MAXSTOKES, MAXSTOKES )
      LOGICAL, INTENT (IN) ::          BVP_REGULAR_FLAG ( 0:MAXMOMENTS )
      !INTEGER, INTENT (IN) ::          LAYER_MAXMOMENTS ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      INTEGER, INTENT (IN) ::          WHICH_DIRECTIONS ( MAX_DIRECTIONS )

      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )
      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      !DOUBLE PRECISION, INTENT (IN) :: PARTLAYERS_VALUES ( MAX_PARTLAYERS )

      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_DN ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DFLUX ( MAXSTOKES )

      !DOUBLE PRECISION, INTENT (IN) :: CHAPMAN_FACTORS &
      !    ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      !DOUBLE PRECISION, INTENT (IN) :: EARTH_RADIUS
      !DOUBLE PRECISION, INTENT (IN) :: HEIGHT_GRID ( 0:MAXLAYERS )

      !DOUBLE PRECISION, INTENT (IN) :: OMEGA_TOTAL_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT_INPUT ( MAXLAYERS )
      !DOUBLE PRECISION, INTENT (IN) :: GREEKMAT_TOTAL_INPUT &
      !    ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
      !DOUBLE PRECISION, INTENT (IN) :: TAUGRID_INPUT ( 0:MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: THERMAL_BB_INPUT ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SURFBB

      !DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_TOTAL_INPUT &
      !    ( MAX_ATMOSWFS, MAXLAYERS )
      !DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT_INPUT &
      !    ( MAX_ATMOSWFS, MAXLAYERS )
      !DOUBLE PRECISION, INTENT (IN) :: L_GREEKMAT_TOTAL_INPUT &
      !     ( MAX_ATMOSWFS, 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )

!  These have been removed for version 2.7
!      DOUBLE PRECISION, INTENT (IN) :: LTE_DELTAU_VERT_INPUT ( 2, MAXLAYERS )
!      DOUBLE PRECISION, INTENT (IN) :: LTE_THERMAL_BB_INPUT ( 0:MAXLAYERS )

      LOGICAL, INTENT (INOUT) ::          DO_SIMULATION_ONLY        ! Intent changed, Version 2.7
      LOGICAL, INTENT (IN) ::          DO_ATMOS_LINEARIZATION
      LOGICAL, INTENT (IN) ::          DO_PROFILE_LINEARIZATION
      LOGICAL, INTENT (IN) ::          DO_SURFACE_LINEARIZATION
      INTEGER, INTENT (IN) ::          N_TOTALPROFILE_WFS
      LOGICAL, INTENT (IN) ::          LAYER_VARY_FLAG  ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          LAYER_VARY_NUMBER ( MAXLAYERS )

!  Control for  Blackbody Jacobians, New 28 March 2014. Version 2.7
!   Replaces the two commented out flags

      LOGICAL, INTENT (INOUT)  :: DO_ATMOS_LBBF, DO_SURFACE_LBBF
!      LOGICAL ::            DO_LTE_LINEARIZATION
!      LOGICAL ::            DO_SURFBB_LINEARIZATION


!  Work type structures

      TYPE(VLIDORT_Work_Miscellanous), INTENT (IN) :: Misc
      TYPE(VLIDORT_Work_Thermal),      INTENT (IN) :: Therm
      TYPE(VLIDORT_Work_Multiplier),   INTENT (IN) :: Mult
      !TYPE(VLIDORT_Work_Corrections),  INTENT (IN) :: Corr
      !TYPE(VLIDORT_Work_FirstOrder),   INTENT (IN) :: Fo

      TYPE(VLIDORT_LinWork_Miscellanous), INTENT (IN) :: LAP_Misc
      TYPE(VLIDORT_LinWork_Thermal),      INTENT (IN) :: L_Therm
      TYPE(VLIDORT_LinWork_Multiplier),   INTENT (IN) :: LP_Mult
      !TYPE(VLIDORT_Work_Corrections),  INTENT (IN) :: Corr
      !TYPE(VLIDORT_Work_FirstOrder),   INTENT (IN) :: Fo

!  From BRDF supplement

      DOUBLE PRECISION, INTENT (IN) :: LAMBERTIAN_ALBEDO
      !DOUBLE PRECISION, INTENT (IN) :: EXACTDB_BRDFUNC &
      !    ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: BRDF_F &
          ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: BRDF_F_0 &
          ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS   )
      DOUBLE PRECISION, INTENT (IN) :: USER_BRDF_F &
          ( 0:MAXMOMENTS,MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: USER_BRDF_F_0 &
          ( 0:MAXMOMENTS,MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS   )

      INTEGER, INTENT (IN) ::          N_SURFACE_WFS
      !DOUBLE PRECISION, INTENT (IN) :: LS_EXACTDB_BRDFUNC &
      !    ( MAX_SURFACEWFS, MAXSTOKES_SQ, &
      !      MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )
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

      DOUBLE PRECISION, INTENT (IN) :: EMISSIVITY ( MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: USER_EMISSIVITY &
          ( MAXSTOKES, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: LS_USER_EMISSIVITY &
          ( MAX_SURFACEWFS, MAXSTOKES,MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: LS_EMISSIVITY &
          ( MAX_SURFACEWFS, MAXSTOKES,MAXSTREAMS )

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  New Surface-Leaving stuff 17 May 2012
      LOGICAL, INTENT (IN) ::            DO_SURFACE_LEAVING
      LOGICAL, INTENT (IN) ::            DO_SL_ISOTROPIC
      DOUBLE PRECISION, INTENT (IN) ::   SLTERM_ISOTROPIC &
          ( MAXSTOKES, MAXBEAMS )
      !DOUBLE PRECISION, INTENT (IN) ::   SLTERM_USERANGLES &
      !    ( MAXSTOKES, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS  )
      DOUBLE PRECISION, INTENT (IN) ::   SLTERM_F_0 &
          ( 0:MAXMOMENTS, MAXSTOKES, MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::   USER_SLTERM_F_0 &
          ( 0:MAXMOMENTS, MAXSTOKES, MAX_USER_STREAMS, MAXBEAMS )

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@ Addition of SLEAVE WF stuff, R. Spurr, 22 August 2012 @@@@@@@@@

      LOGICAL         , INTENT (IN) ::   DO_SLEAVE_WFS
      INTEGER         , INTENT (IN) ::   N_SLEAVE_WFS

      DOUBLE PRECISION, INTENT (IN) ::   LSSL_SLTERM_ISOTROPIC &
          ( MAX_SLEAVEWFS, MAXSTOKES, MAXBEAMS )
      !DOUBLE PRECISION, INTENT (IN) ::   LSSL_SLTERM_USERANGLES &
      !    ( MAX_SLEAVEWFS, MAXSTOKES, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS  )
      DOUBLE PRECISION, INTENT (IN) ::   LSSL_SLTERM_F_0 &
          ( MAX_SLEAVEWFS, 0:MAXMOMENTS, MAXSTOKES, MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::   LSSL_USER_SLTERM_F_0 &
          ( MAX_SLEAVEWFS, 0:MAXMOMENTS, MAXSTOKES, MAX_USER_STREAMS, MAXBEAMS )

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      DOUBLE PRECISION, INTENT (INOUT) :: PIMM_11 ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, INTENT (INOUT) :: PIMM_KM ( MAX_ALLSTRMS_P1 )

      LOGICAL, INTENT (INOUT) ::           DO_BVTEL_INITIAL
      INTEGER, INTENT (INOUT) ::           BVTEL_FOURIER_COMPONENT
      !DOUBLE PRECISION, INTENT (INOUT) ::  TAUGRID ( 0:MAXLAYERS )
      !DOUBLE PRECISION, INTENT (INOUT) ::  OMEGA_TOTAL ( MAXLAYERS )
      !DOUBLE PRECISION, INTENT (INOUT) ::  GREEKMAT_TOTAL &
      !    ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )

      !DOUBLE PRECISION, INTENT (INOUT) ::  STOKES_SS &
      !    ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      !DOUBLE PRECISION, INTENT (INOUT) ::  STOKES_DB &
      !    ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      !DOUBLE PRECISION, INTENT (INOUT) ::  SS_CONTRIBS &
      !    ( MAX_GEOMETRIES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) ::  MEAN_STOKES &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) ::  FLUX_STOKES &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES ,MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) ::  MEAN_DIRECT &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) ::  FLUX_DIRECT &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES )

      !DOUBLE PRECISION, INTENT (INOUT) :: PROFILEWF_SS &
      !    ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
      !      MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      !DOUBLE PRECISION, INTENT (INOUT) :: PROFILEWF_DB &
      !    ( MAX_ATMOSWFS, MAXLAYERS, &
      !      MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      !DOUBLE PRECISION, INTENT (INOUT) :: SURFACEWF_DB &
      !    ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

      DOUBLE PRECISION, INTENT (INOUT) ::  MINT_PROFILEWF &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) ::  MINT_PROFILEWF_DIRECT &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) ::  FLUX_PROFILEWF &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) ::  FLUX_PROFILEWF_DIRECT &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES )

      DOUBLE PRECISION, INTENT (INOUT) ::  MINT_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) ::  FLUX_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  New Code Version 2.7. LBBF outputs. Superceded LTE_ATMOSWF output (now dropped)
!     Outputs are all Pre-zeroed in the calling Masters
!     Postprocessed and Flux Jacobians.

      DOUBLE PRECISION, INTENT(INOUT) :: ABBWFS_JACOBIANS &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, 0:MAXLAYERS, MAXSTOKES, MAX_DIRECTIONS)
      DOUBLE PRECISION, INTENT(INOUT) :: SBBWFS_JACOBIANS &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAXSTOKES, MAX_DIRECTIONS)

      DOUBLE PRECISION, INTENT(INOUT) :: ABBWFS_FLUXES &
          ( MAX_USER_LEVELS, 2, 0:MAXLAYERS, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT(INOUT) :: SBBWFS_FLUXES &
          ( MAX_USER_LEVELS, 2, MAXSTOKES, MAX_DIRECTIONS)

!  Dropped code from Version 2p6
!      DOUBLE PRECISION, INTENT (INOUT) ::  LTE_ATMOSWF &
!          ( 0:MAXLAYERS, MAX_USER_LEVELS, &
!            MAX_USER_VZANGLES, MAX_DIRECTIONS )

      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE_1
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE_2

      LOGICAL, INTENT (OUT) ::           DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (OUT) ::           DO_INCLUDE_SURFACE
      LOGICAL, INTENT (OUT) ::           DO_INCLUDE_SURFEMISS

      DOUBLE PRECISION, INTENT (OUT) ::  MS_CONTRIBS_F &
          ( MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAXLAYERS  )
      DOUBLE PRECISION, INTENT (OUT) ::  STOKES_F &
          ( MAX_USER_LEVELS, MAX_USER_VZANGLES, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (OUT) ::  PROFILEWF_F &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
            MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (OUT) ::  SURFACEWF_F &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_VZANGLES, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      INTEGER, INTENT (OUT) ::           STATUS

!  Local variables
!  ---------------

!  FROM VLIDORT_MISCSETUPS

      DOUBLE PRECISION :: DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION :: PARTAU_VERT ( MAX_PARTLAYERS )
      DOUBLE PRECISION :: DELTAU_SLANT &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      !DOUBLE PRECISION :: TRUNC_FACTOR ( MAXLAYERS )
      !DOUBLE PRECISION :: FAC1 ( MAXLAYERS )
      DOUBLE PRECISION :: OMEGA_GREEK &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      INTEGER ::          LAYER_PIS_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION :: TRANS_SOLAR_BEAM ( MAXBEAMS )
      LOGICAL ::          DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )
      DOUBLE PRECISION :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION :: T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION :: T_DISORDS_UTDN ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )
      !DOUBLE PRECISION :: T_UTUP_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION :: CUMTRANS ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      !DOUBLE PRECISION :: ITRANS_USERM &
      !    ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: LOCAL_CSZA ( 0:MAXLAYERS, MAXBEAMS )

!  FROM THERMAL_SETUP

      DOUBLE PRECISION :: THERMCOEFFS ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION :: DELTAU_POWER ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION :: XTAU_POWER ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION :: TCOM1 ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION :: T_DIRECT_UP ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: T_DIRECT_DN ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: T_UT_DIRECT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION :: T_UT_DIRECT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS )

!  FROM EMULT_MASTER

      !LOGICAL ::          EMULT_HOPRULE &
      !    ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      !DOUBLE PRECISION :: SIGMA_M &
      !    ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      !DOUBLE PRECISION :: SIGMA_P &
      !    ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION :: EMULT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: EMULT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: UT_EMULT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: UT_EMULT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )

!  FROM VLIDORT_L_MISCSETUPS

      !DOUBLE PRECISION :: L_OMEGA_TOTAL &
      !    ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION :: L_DELTAU_VERT &
          ( MAX_ATMOSWFS, MAXLAYERS )
      !DOUBLE PRECISION :: L_GREEKMAT_TOTAL &
      !    ( MAX_ATMOSWFS, 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )
      !DOUBLE PRECISION :: L_DELTAU_SLANT &
      !    ( MAX_ATMOSWFS, MAXLAYERS, MAXLAYERS, MAXBEAMS )
      !LOGICAL ::          DO_SCATMAT_VARIATION &
      !    ( MAXLAYERS, MAX_ATMOSWFS )
      !DOUBLE PRECISION :: L_TRUNC_FACTOR &
      !    ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION :: L_OMEGA_GREEK &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LP_AVERAGE_SECANT &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LP_INITIAL_TRANS &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_DELT_DISORDS &
          ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_DISORDS_UTDN &
          ( MAXSTREAMS, MAX_USER_LEVELS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_DISORDS_UTUP &
          ( MAXSTREAMS, MAX_USER_LEVELS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LP_T_DELT_MUBAR &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LP_T_UTDN_MUBAR &
          ( MAX_USER_LEVELS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)
      DOUBLE PRECISION :: L_T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_UTDN_USERM &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_UTUP_USERM &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS )

!  FROM THERMAL_SETUP_PLUS

      DOUBLE PRECISION :: L_THERMCOEFFS &
          ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_DELTAU_POWER &
          ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_XTAU_POWER &
          ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_TCOM1 &
          ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_DIRECT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_DIRECT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_UT_DIRECT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS)
      DOUBLE PRECISION :: L_T_UT_DIRECT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  FROM LP_EMULT_MASTER

      DOUBLE PRECISION :: LP_EMULT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXLAYERS, &
            MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LP_EMULT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXLAYERS, &
            MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LP_UT_EMULT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, &
            MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LP_UT_EMULT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, &
            MAXBEAMS, MAX_ATMOSWFS )

!  Local inclusion flags

      LOGICAL ::          DO_INCLUDE_MVOUTPUT
      LOGICAL ::          DO_INCLUDE_DIRECTBEAM

!  Flux multiplier and Fourier component numbers

      DOUBLE PRECISION :: FLUX_MULTIPLIER
      DOUBLE PRECISION :: DELTA_FACTOR
      DOUBLE PRECISION :: SURFACE_FACTOR

!  Weighting function indices

      INTEGER ::          N_LAYER_WFS, LAYER_TO_VARY, N_PARAMETERS
      LOGICAL ::          DO_RTSOL_VARY
      INTEGER ::          NPARAMS_VARY, VARIATION_INDEX

!  Error tracing variables

      CHARACTER (LEN=2) :: CF
      INTEGER ::           STATUS_SUB
      LOGICAL ::           FAIL

!  Progress

      LOGICAL, PARAMETER :: DO_WRITE_SCREEN = .FALSE.

!  Helper

      INTEGER ::          LAYER, IBEAM, IPARTIC, I, O1, N

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

!mick fix 7/29/2014 - added to make VLIDORT threadsafe
      DOUBLE PRECISION :: COL2 ( MAXTOTAL, MAXBEAMS )
      DOUBLE PRECISION :: COLTEL2 ( MAXTOTAL, MAXBEAMS )
      DOUBLE PRECISION :: SCOL2 ( MAXSTRMSTKS_2, MAXBEAMS )

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
      DOUBLE PRECISION :: LP_BVEC &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_UPAR_DN_1 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_UPAR_UP_1 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LP_UPAR_DN_2 ( MAX_USER_STREAMS, &
          MAXSTOKES, MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LP_UPAR_UP_2 ( MAX_USER_STREAMS, &
          MAXSTOKES, MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_WUPPER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_WLOWER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: NCON &
          ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: PCON &
          ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@ Addition of SLEAVE WF stuff, R. Spurr, 22 August 2012 @@@@@@@@@

!FROM VLIDORT_LSDL_DBSETUPS
      DOUBLE PRECISION :: LSSL_DIRECT_BEAM &
          ( MAX_SLEAVEWFS, MAXSTREAMS, MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION :: LSSL_USER_DIRECT_BEAM &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      INTEGER :: LVARY

!TESTING
      INTEGER :: K,Q,UM,UTA,V
      LOGICAL :: DO_DEBUG=.FALSE.,&
                 DO_LIN_DEBUG=.TRUE.

!  ##############
!  initialization
!  ##############

!  module status and message initialization

      STATUS  = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE_1 = ' '
      TRACE_2 = ' '

!  Local simulation only
!   Definition extended, 28 March 2014

      DO_SIMULATION_ONLY = ( .NOT. DO_PROFILE_LINEARIZATION .AND. &
                             .NOT. DO_SURFACE_LINEARIZATION .and. &
                   .not.do_ATMOS_LBBF .and. .not. do_surface_LBBF )

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

!  ##################
!  mick chg 7/17/2014 - moved misc setup, thermal setup, emult setup, sscor codes (2),
!                       db code, and fo code up to VLIDORT_LPS_MASTER.  There, the needed
!                       vars are packed into type structures to bypass f90 continuation
!                       line limits and passed down.  Here, they are unpacked.

      CALL VLIDORT_UNPACK_MISC ( Misc,                                        & ! Input 
        NSTOKES, NSTREAMS, NLAYERS, N_PARTLAYERS, N_USER_STREAMS, NBEAMS,     & ! Input 
        NMOMENTS,                                                             & ! Input    
        DELTAU_VERT, PARTAU_VERT, DELTAU_SLANT, OMEGA_GREEK,                  & ! Output
        LAYER_PIS_CUTOFF, TRANS_SOLAR_BEAM, DO_REFLECTED_DIRECTBEAM,          & ! Output
        T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,                       & ! Output
        T_DELT_MUBAR, T_UTDN_MUBAR, T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM, & ! Output
        CUMTRANS, INITIAL_TRANS, AVERAGE_SECANT, LOCAL_CSZA )                   ! Output   

      IF ( DO_THERMAL_EMISSION ) THEN
        CALL VLIDORT_UNPACK_THERM ( Therm,                          & ! Input
          NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_STREAMS,  & ! Input
          THERMCOEFFS, DELTAU_POWER, XTAU_POWER, TCOM1,             & ! Output
          T_DIRECT_UP, T_DIRECT_DN, T_UT_DIRECT_UP, T_UT_DIRECT_DN )  ! Output
      END IF

      IF ( DO_SOLAR_SOURCES ) THEN
        IF (.NOT.DO_SSFULL.OR.(DO_SSFULL.AND.DO_SSCORR_NADIR)) THEN
          CALL VLIDORT_UNPACK_MULT ( Mult,                 & ! Input
            NLAYERS, N_PARTLAYERS, NBEAMS, N_USER_STREAMS, & ! Input
            EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN )   ! Output
        ENDIF
      ENDIF

      IF ( DO_ATMOS_LINEARIZATION ) THEN
        CALL VLIDORT_UNPACK_LAP_MISC ( LAP_Misc,                & ! Input
          NSTOKES, NSTOKES_SQ, NSTREAMS, NLAYERS, N_PARTLAYERS, & ! Input
          NBEAMS, N_USER_STREAMS, NMOMENTS, N_TOTALPROFILE_WFS, & ! Input
          L_DELTAU_VERT, L_OMEGA_GREEK,                         & ! Output
          LP_AVERAGE_SECANT, LP_INITIAL_TRANS,                  & ! Output
          L_T_DELT_DISORDS, L_T_DISORDS_UTDN, L_T_DISORDS_UTUP, & ! Output
          LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR,                     & ! Output
          L_T_DELT_USERM, L_T_UTDN_USERM, L_T_UTUP_USERM )        ! Output

        IF ( DO_THERMAL_EMISSION ) THEN
          CALL VLIDORT_UNPACK_L_THERM ( L_Therm,                                         & ! Input
            NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_STREAMS, N_TOTALPROFILE_WFS, & ! Input
            L_THERMCOEFFS, L_DELTAU_POWER, L_XTAU_POWER, L_TCOM1,                        & ! Output
            L_T_DIRECT_UP, L_T_DIRECT_DN, L_T_UT_DIRECT_UP, L_T_UT_DIRECT_DN )             ! Output
        END IF

        IF ( DO_SOLAR_SOURCES ) THEN
          IF (.NOT.DO_SSFULL.OR.(DO_SSFULL.AND.DO_SSCORR_NADIR)) THEN
            CALL VLIDORT_UNPACK_LP_MULT ( LP_Mult,                     & ! Input
              NLAYERS, N_PARTLAYERS, NBEAMS, N_USER_STREAMS,           & ! Input
              N_TOTALPROFILE_WFS,                                      & ! Input
              LP_EMULT_UP, LP_EMULT_DN, LP_UT_EMULT_UP, LP_UT_EMULT_DN ) ! Output
          ENDIF
        ENDIF
      ENDIF

!  ##################

!  mick fix 7/17/2014 - if block moved from above and modified to accomodate move of
!                       "Fourier = 0" subroutines
!  Direct beam flag (only if above albedo flag has been set)
!   (DO_REFLECTED_DIRECTBEAM is a stored Bookkeeping variable)

      IF ( FOURIER_COMPONENT .GT. 0 ) THEN
        IF ( DO_DIRECT_BEAM .AND. DO_SOLAR_SOURCES ) THEN
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
      ENDIF

!  Not required if no solar sources

      IF ( DO_SOLAR_SOURCES ) THEN

!  Surface direct beam
!    New Surface-Leaving arguments added, 17 May 2012

        IF ( .not. DO_SSFULL ) THEN
          CALL VLIDORT_DIRECTBEAM ( &
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
            DO_OBSERVATION_GEOMETRY, &
            TRANS_SOLAR_BEAM, DO_REFLECTED_DIRECTBEAM, &
            ATMOS_ATTN, DIRECT_BEAM, &
            USER_DIRECT_BEAM )
        ENDIF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@ Addition of SLEAVE weighting function setups @@@@@@@@@@@@
!    R. Spurr, 22 August 2012

       IF ( DO_SURFACE_LEAVING .and. DO_SLEAVE_WFS ) then
         CALL VLIDORT_LSSL_DBSETUPS ( &
           DO_OBSERVATION_GEOMETRY, DO_SL_ISOTROPIC, DO_USER_STREAMS,    &
           DO_REFLECTED_DIRECTBEAM, FOURIER_COMPONENT,                   &
           NSTOKES, NSTREAMS, NBEAMS, N_USER_STREAMS, N_SLEAVE_WFS,      &
           FLUX_FACTOR, DELTA_FACTOR,                                    &
           LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0, LSSL_USER_SLTERM_F_0, &
           LSSL_DIRECT_BEAM, LSSL_USER_DIRECT_BEAM )
       ENDIF

!@@@@@@@@@@ END Addition of SLEAVE weighting function setups @@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  End solar sources only clause for corrections

      ENDIF

!  ###################
!  Spherical functions
!  ###################

!  Get Pi Matrices for this Fourier component
!    Including all beam angles !!

      IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
        CALL VLIDORT_PIMATRIX_SETUP_OMP ( &
          FOURIER_COMPONENT, &
          DO_REFRACTIVE_GEOMETRY, NSTOKES, &
          NSTREAMS, NLAYERS, &
          COS_SZANGLES, SUN_SZA_COSINES, &
          QUAD_STREAMS, NMOMENTS, &
          NBEAMS, N_USER_STREAMS, &
          DO_USER_STREAMS, USER_STREAMS, &
          MUELLER_INDEX, DMAT, &
          PIMM_11, PIMM_KM, &
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
          'Called in VLIDORT_LPS_FOURIER, Fourier component '//CF
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
             'Called in VLIDORT_LPS_FOURIER, Fourier component '//CF
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
        DO_UPWELLING, DO_DNWELLING, TAYLOR_ORDER, & ! 2p7 Taylor-order
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

      IF ( DO_PROFILE_LINEARIZATION ) THEN
        CALL L_HMULT_MASTER ( &
          DO_UPWELLING, DO_DNWELLING, TAYLOR_ORDER, & ! 2p7 Taylor-order
          NLAYERS, N_USER_LEVELS, &
          N_USER_STREAMS, LOCAL_UM_START, &
          USER_SECANTS, &
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
          PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, &
          STERM_LAYERMASK_DN, &
          DELTAU_VERT, PARTAU_VERT, T_DELT_EIGEN, &  !@@@ Rob Fix 5/10/13 added
          T_UTUP_EIGEN, T_UTDN_EIGEN, &
          T_DELT_USERM, T_UTDN_USERM, &
          T_UTUP_USERM, &
          K_REAL, K_COMPLEX, HSINGO, &
          HMULT_1, HMULT_2, ZETA_M, ZETA_P, &
          LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
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
            'Called in VLIDORT_LPS_FOURIER, Fourier # '//CF
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
          DO_BVTEL_INITIAL, BVTEL_FOURIER_COMPONENT, &
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
            'Called in VLIDORT_LPS_FOURIER, Fourier # '//CF
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
            TRACE_2 = 'Called in VLIDORT_LPS_FOURIER, Fourier 0'
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
            TRACE_2 = 'Called in VLIDORT_LPS_FOURIER, Fourier 0'
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
        COL2, SCOL2, &
        R2_BEAM, LCON, MCON, &
        STATUS_SUB, MESSAGE, TRACE_1 )

!  Exception handling

      IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
        write(CF,'(I2)')FOURIER_COMPONENT
        TRACE_2 = 'Error return from BVP_SOLUTION_MASTER, '// &
           'Called in VLIDORT_LPS_FOURIER, Fourier # '//CF
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
          DO_USER_STREAMS, DO_OBSERVATION_GEOMETRY, LOCAL_UM_START, &
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
          DO_USER_STREAMS, DO_OBSERVATION_GEOMETRY, &
          LOCAL_UM_START, &
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

!  This routine is now defunct in Verison 2.7. Replaced by LBBF stuff. 3/28/14

!  LTE linearization (new section, 14 September 2009)
!  --------------------------------------------------

!      IF ( DO_THERMAL_TRANSONLY .AND. DO_LTE_LINEARIZATION ) THEN
!        CALL THERMAL_LTE_LINEARIZATION ( &
!          DO_INCLUDE_SURFACE, SURFACE_FACTOR, &
!          FLUX_MULTIPLIER, &
!          DO_UPWELLING, DO_DNWELLING, &
!          NSTREAMS, NLAYERS, N_USER_LEVELS, &
!          DELTAU_VERT_INPUT, LAMBERTIAN_ALBEDO, &
!          THERMAL_BB_INPUT, &
!          QUAD_STREAMS, QUAD_STRMWTS, &
!          N_USER_STREAMS, LOCAL_UM_START, USER_STREAMS, &
!          UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, &
!          STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, &
!          T_DELT_DISORDS, T_DELT_USERM, &
!          CUMSOURCE_UP, CUMSOURCE_DN, &
!          THERMCOEFFS, &
!          LTE_DELTAU_VERT_INPUT, LTE_THERMAL_BB_INPUT, &
!          LTE_ATMOSWF )
!      ENDIF

!  Thermal only. Avoid weighting functions all together
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      IF ( .not. DO_SIMULATION_ONLY ) THEN

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

              CALL LP_BVP_SOLUTION_MASTER ( &
                DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM, &
                DO_INCLUDE_THERMEMISS, &
                LAYER_TO_VARY, N_LAYER_WFS, &
                FOURIER_COMPONENT, IPARTIC, &
                SURFACE_FACTOR, DO_PROFILE_LINEARIZATION, &
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
                LP_INITIAL_TRANS, LP_T_DELT_MUBAR, LP_BVEC, &
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
                TRACE_2 = 'Error return from LP_BVP_SOLUTION_MASTER, '// &
        'Profile Jacobians, Called in VLIDORT_LPS_FOURIER, Fourier # '//CF
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
              CALL UPUSER_PROFILEWF ( &
                DO_OBSERVATION_GEOMETRY, &
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
                L_UHOM_UPDN, L_UHOM_UPUP, L_UPAR_UP_1, LP_UPAR_UP_2, &
                L_HMULT_1, L_HMULT_2, LP_EMULT_UP, &
                L_LAYER_TSUP_UP, &
                UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP, &
                L_UT_HMULT_UU, L_UT_HMULT_UD, LP_UT_EMULT_UP, &
                L_LAYER_TSUP_UTUP, &
                L_BOA_THTONLY_SOURCE, PROFILEWF_F )
            ENDIF

            IF ( DO_DNWELLING ) THEN
              CALL DNUSER_PROFILEWF ( &
                DO_OBSERVATION_GEOMETRY, &
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
                L_UPAR_DN_1, LP_UPAR_DN_2, &
                L_HMULT_1, L_HMULT_2, LP_EMULT_DN, &
                L_LAYER_TSUP_DN, &
                UT_HMULT_DU, UT_HMULT_DD, UT_EMULT_DN, &
                L_UT_HMULT_DU, L_UT_HMULT_DD, &
                LP_UT_EMULT_DN, &
                L_LAYER_TSUP_UTDN, &
                PROFILEWF_F )
            ENDIF

! @@@ Rob fix 1/31/11, - added FLUX_FACTOR argument

            IF ( DO_INCLUDE_MVOUTPUT .OR. DO_QUAD_OUTPUT ) THEN
              CALL MIFLUX_PROFILEWF ( &
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
                LP_INITIAL_TRANS, LP_T_DELT_MUBAR, &
                LP_T_UTDN_MUBAR, &
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
                MINT_PROFILEWF, MINT_PROFILEWF_DIRECT, &
                FLUX_PROFILEWF, FLUX_PROFILEWF_DIRECT )
            ENDIF

!  Finish loop over layers with variation

          ENDIF
        ENDDO

!  End atmospheric profile weighting functions

      ENDIF

!  Thermal-only: Surface Reflectance weighting functions
!  -----------------------------------------------------

!mick fix 9/6/2012 - added (N_SURFACE_WFS > 0) condition to IF
      IF ( DO_SURFACE_LINEARIZATION .AND. DO_INCLUDE_SURFACE .AND. &
           (N_SURFACE_WFS > 0) ) THEN

          CALL SURFACEWF_MASTER ( &
            DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFEMISS, &
            DO_MSMODE_THERMAL, DO_OBSERVATION_GEOMETRY, &
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
            TRACE_2 = 'Error return from SURFACEWF_MASTER, '// &
            'Thermal-only Called in VLIDORT_LPS_FOURIER, Fourier # '//CF
            STATUS = VLIDORT_SERIOUS
            RETURN
          ENDIF

!  end of surface weighting functions

      ENDIF

!  New, 28 March 2014. Linearization for BLACKBODY

!      IF ( ( DO_ATMOS_LBBF .or. DO_SURFACE_LBBF ) .and. Fourier_component.eq.0 ) THEN
!         if ( n_partlayers .gt. 0 ) then
!            call vlidort_lbbf_jacobians_wpartials &
!               ( DO_ATMOS_LBBF, DO_SURFACE_LBBF, DO_THERMAL_TRANSONLY,      & ! Input
!                 DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES,              & ! Input
!                 DO_MSMODE_THERMAL, DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT,   & ! input
!                 DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,                       & ! input
!                 NSTOKES, NLAYERS, NSTREAMS, N_USER_STREAMS, N_USER_LEVELS, & ! input
!                 NMOMENTS, NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2,        & ! Input
!                 NTOTAL, N_SUPDIAG, N_SUBDIAG,  MUELLER_INDEX,              & ! input
!                 N_PARTLAYERS, PARTLAYERS_LAYERIDX,                         & ! Input
!                 PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                   & ! Input
!                 UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,                    & ! Input
!                 USER_STREAMS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,      & ! Input
!                 QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS,                  & ! input
!                 SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,               & ! input
!                 EMISSIVITY, USER_EMISSIVITY, FLUX_MULTIPLIER,              & ! input
!                 OMEGA_TOTAL, DELTAU_VERT, PARTAU_VERT, OMEGA_GREEK,        & ! Input
!                 T_DELT_DISORDS, T_UTDN_DISORDS, T_UTUP_DISORDS,            & ! input
!                 PI_XQP, PI_XQP_PRE, PI_XUP, PI_XUM, SAB, DAB,              & ! input
!                 K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,                   & ! input
!                 T_DELT_EIGEN, T_UTDN_EIGEN, T_UTUP_EIGEN,                  & ! input
!                 T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM,                  & ! Input
!                 BANDMAT2, IPIVOT, SMAT2, SIPIVOT, HMULT_1, HMULT_2,        & ! Input
!                 UT_HMULT_UU, UT_HMULT_UD, UT_HMULT_DU, UT_HMULT_DD,        & ! Input
!                 UHOM_DNDN, UHOM_DNUP, UHOM_UPDN, UHOM_UPUP,                & ! input
!                 ABBWFS_JACOBIANS, ABBWFS_FLUXES,                           & ! Output
!                 SBBWFS_JACOBIANS, SBBWFS_FLUXES,                           & ! Output
!                 STATUS_SUB, MESSAGE, TRACE_1 )                               ! Output
!         else
!            call vlidort_lbbf_jacobians_whole &
!               ( DO_ATMOS_LBBF, DO_SURFACE_LBBF, DO_THERMAL_TRANSONLY,      & ! Input
!                 DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES,              & ! Input
!                 DO_MSMODE_THERMAL, DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT,   & ! input
!                 DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,                       & ! input
!                 NSTOKES, NLAYERS, NSTREAMS, N_USER_STREAMS, N_USER_LEVELS, & ! input
!                 NMOMENTS, NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2,        & ! Input
!                 NTOTAL, N_SUPDIAG, N_SUBDIAG,  MUELLER_INDEX,              & ! input
!                 UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,                    & ! Input
!                 USER_STREAMS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,      & ! Input
!                 QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS,                  & ! input
!                 SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,               & ! input
!                 EMISSIVITY, USER_EMISSIVITY,                               & ! input
!                 FLUX_MULTIPLIER, OMEGA_TOTAL, DELTAU_VERT, OMEGA_GREEK,    & ! Input
!                 T_DELT_DISORDS, T_DELT_USERM,                              & ! input
!                 PI_XQP, PI_XQP_PRE, PI_XUP, PI_XUM, SAB, DAB,              & ! input
!                 K_REAL, K_COMPLEX, T_DELT_EIGEN, SOLA_XPOS, SOLB_XNEG,     & ! input
!                 BANDMAT2, IPIVOT, SMAT2, SIPIVOT, HMULT_1, HMULT_2,        & ! Input
!                 UHOM_DNDN, UHOM_DNUP, UHOM_UPDN, UHOM_UPUP,                & ! input
!                 ABBWFS_JACOBIANS, ABBWFS_FLUXES,                           & ! Output
!                 SBBWFS_JACOBIANS, SBBWFS_FLUXES,                           & ! Output
!                 STATUS_SUB, MESSAGE, TRACE_1 )                               ! Output
!         endif

!  Error handling

!         IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
!            write(CF,'(I2)')FOURIER_COMPONENT
!            TRACE_2 = 'Error return from vlidort_lbbf_jacobians, thermal only'// &
!                        'Called in VLIDORT_LPS_FOURIER_MASTER, Fourier # '//CF
!            STATUS = VLIDORT_SERIOUS
!            RETURN
!         ENDIF

!  ENd of BB Jacobians

!        ENDIF

!  End simulation clause

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
                  'Called in VLIDORT_LPS_FOURIER, Fourier # '//CF
                STATUS = VLIDORT_SERIOUS
                RETURN
              ENDIF

!  For the linearizations of the classical solution

              IF ( DO_ATMOS_LINEARIZATION ) THEN

                CALL VLIDORT_LP_QBEAM_SOLUTION ( &
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
                  LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
                  L_SAB, L_DAB, L_EIGENMAT, &
                  L_OMEGA_GREEK, LP_AVERAGE_SECANT, &
                  LP_BVEC, &
                  STATUS_SUB, MESSAGE, TRACE_1 )

                IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
                  write(CF,'(I2)')FOURIER_COMPONENT
                  TRACE_2 = 'Error from VLIDORT_LP_QBEAM_SOLUTION, '// &
                          'Called in VLIDORT_LPS_FOURIER, Fourier # '//CF
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
                    DO_OBSERVATION_GEOMETRY, &
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
                    CALL VLIDORT_LP_UBEAM_SOLUTION ( &
                      LAYER, FOURIER_COMPONENT, IBEAM, &
                      DO_RTSOL_VARY, NPARAMS_VARY, &
                      DO_PLANE_PARALLEL, DO_UPWELLING, &
                      DO_DNWELLING, DO_OBSERVATION_GEOMETRY, &
                      NSTOKES, NSTREAMS, FLUX_FACTOR, &
                      LAYER_PIS_CUTOFF, QUAD_HALFWTS, &
                      NMOMENTS, N_USER_STREAMS, &
                      DO_LAYER_SCATTERING, LOCAL_UM_START, &
                      STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, &
                      DFLUX, OMEGA_GREEK, &
                      PI_XQP, PI_XUP, PI_XUM, PI_X0P, &
                      PI_XQM_PRE, HELPSTOKES_BEAM, &
                      LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
                      LP_BVEC, L_OMEGA_GREEK, &
                      L_UPAR_DN_1, L_UPAR_UP_1, &
                      LP_UPAR_DN_2, LP_UPAR_UP_2 )
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
              COL2, SCOL2, &
              R2_BEAM, LCON, MCON, &
              STATUS_SUB, MESSAGE, TRACE_1 )

            IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
              write(CF,'(I2)')FOURIER_COMPONENT
              TRACE_2 = 'Error return from BVP_SOLUTION_MASTER, '// &
                'Called in VLIDORT_LPS_FOURIER, Fourier # '//CF
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
              AXBID_F, &
              COLTEL2, SCOL2, &
              R2_BEAM, LCON, MCON, &
              STATUS_SUB, MESSAGE, TRACE_1 )

            IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
              write(CF,'(I2)')FOURIER_COMPONENT
              TRACE_2 = 'Error return from BVPTEL_SOLUTION_MASTER, '// &
                'Called in VLIDORT_LPS_FOURIER, Fourier # '//CF
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
              DO_USER_STREAMS, DO_OBSERVATION_GEOMETRY, LOCAL_UM_START, &
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
              DO_USER_STREAMS, DO_OBSERVATION_GEOMETRY, &
              LOCAL_UM_START, &
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

              CALL LP_BVP_SOLUTION_MASTER ( &
                DO_INCLUDE_SURFACE, &
                DO_REFLECTED_DIRECTBEAM(IBEAM), &
                DO_INCLUDE_THERMEMISS, &
                LAYER_TO_VARY, N_LAYER_WFS, &
                FOURIER_COMPONENT, IBEAM, &
                SURFACE_FACTOR, DO_PROFILE_LINEARIZATION, &
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
                LP_INITIAL_TRANS, LP_T_DELT_MUBAR, LP_BVEC, &
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
                TRACE_2 = 'Error return from LP_BVP_SOLUTION_MASTER, '// &
        'Profile Jacobians, Called in VLIDORT_LPS_FOURIER, Fourier # '//CF
                STATUS = VLIDORT_SERIOUS
                RETURN
              ENDIF

!    (b) telescoped case

            ELSE

             CALL LP_BVPTEL_SOLUTION_MASTER ( &
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
               DO_PROFILE_LINEARIZATION, &
               L_T_DELT_EIGEN, L_T_DELT_DISORDS, &
               L_SOLA_XPOS, L_SOLB_XNEG, &
               DO_SPECIALIST_OPTION_2, NSTREAMS_2, &
               DELTAU_SLANT, &
               L_DELTAU_VERT, DIRECT_BEAM, &
               DO_CLASSICAL_SOLUTION, LAYER_PIS_CUTOFF, &
               DO_LAYER_SCATTERING, DO_PLANE_PARALLEL, &
               T_DELT_MUBAR, INITIAL_TRANS, BVEC, &
               LP_INITIAL_TRANS, LP_T_DELT_MUBAR, LP_BVEC, &
               DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
               BRDF_F, QUAD_STRMWTS, MUELLER_INDEX, &
               BANDTELMAT2, IPIVOTTEL, &
               SMAT2, SIPIVOT, &
               L_WLOWER, L_WUPPER, NCON, PCON, &
               STATUS_SUB, MESSAGE, TRACE_1 )

              IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
                write(CF,'(I2)')FOURIER_COMPONENT
                TRACE_2 = 'Error from LP_BVPTEL_SOLUTION_MASTER, '// &
                'Profile Jacobians, Called in VLIDORT_LPS_FOURIER, Fourier # '//CF
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
              CALL UPUSER_PROFILEWF ( &
                DO_OBSERVATION_GEOMETRY, &
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
                L_UHOM_UPDN, L_UHOM_UPUP, L_UPAR_UP_1, LP_UPAR_UP_2, &
                L_HMULT_1, L_HMULT_2, LP_EMULT_UP, &
                L_LAYER_TSUP_UP, &
                UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP, &
                L_UT_HMULT_UU, L_UT_HMULT_UD, LP_UT_EMULT_UP, &
                L_LAYER_TSUP_UTUP, &
                L_BOA_THTONLY_SOURCE, PROFILEWF_F )
             ENDIF

             IF ( DO_DNWELLING ) THEN
!           if ( do_write_screen) write(*,*)'wf dn',ibeam,layer_to_vary
              CALL DNUSER_PROFILEWF ( &
                DO_OBSERVATION_GEOMETRY, &
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
                L_UPAR_DN_1, LP_UPAR_DN_2, &
                L_HMULT_1, L_HMULT_2, LP_EMULT_DN, &
                L_LAYER_TSUP_DN, &
                UT_HMULT_DU, UT_HMULT_DD, UT_EMULT_DN, &
                L_UT_HMULT_DU, L_UT_HMULT_DD, &
                LP_UT_EMULT_DN, &
                L_LAYER_TSUP_UTDN, &
                PROFILEWF_F )
             ENDIF

! @@@ Rob fix 1/31/11, - added FLUX_FACTOR argument

             IF ( DO_INCLUDE_MVOUTPUT .OR. DO_QUAD_OUTPUT ) THEN
              CALL MIFLUX_PROFILEWF ( &
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
                LP_INITIAL_TRANS, LP_T_DELT_MUBAR, &
                LP_T_UTDN_MUBAR, &
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
                MINT_PROFILEWF, MINT_PROFILEWF_DIRECT, &
                FLUX_PROFILEWF, FLUX_PROFILEWF_DIRECT )
             ENDIF

!  Finish loop over layers with variation

            ENDIF
           ENDDO

!  End profile atmospheric weighting functions

          ENDIF

!  Step 5. Surface Reflectance weighting functions
!  -----------------------------------------------

          IF ( DO_SURFACE_LINEARIZATION .AND. DO_INCLUDE_SURFACE ) THEN

!  Only do this if there is surface reflectance !

!mick fix 9/6/2012 - added (N_SURFACE_WFS > 0) IF condition
            IF ( N_SURFACE_WFS > 0 ) THEN

              CALL SURFACEWF_MASTER ( &
                DO_REFLECTED_DIRECTBEAM(IBEAM), DO_INCLUDE_SURFEMISS, &
                DO_MSMODE_THERMAL, DO_OBSERVATION_GEOMETRY, &
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
                TRACE_2 = 'Error return from SURFACEWF_MASTER, '// &
                 'Beam Called in VLIDORT_LPS_FOURIER, Fourier # '//CF
                STATUS = VLIDORT_SERIOUS
                RETURN
              ENDIF

!  End surface-property linearizations

            ENDIF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@ Addition of SLEAVE weighting function @@@@@@@@@@@@@@@@@@@
!    R. Spurr, 22 August 2012
!
            IF ( DO_SURFACE_LEAVING .and. DO_SLEAVE_WFS ) then

!  Direct-beam Flag (repeated here, just so we know!)

              DO_INCLUDE_DIRECTBEAM = ( DO_UPWELLING .AND. &
                (DO_REFLECTED_DIRECTBEAM(IBEAM).AND..NOT.DO_DBCORRECTION) ) &
                 .AND. .NOT.DO_MSMODE_VLIDORT

!  Rob Fix @@@ 11 Sep 12, Add New Line to the Argument list

              CALL VLIDORT_LSSL_WFS ( &
                DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_MVOUTPUT,                  &
                DO_QUAD_OUTPUT, DO_USER_STREAMS, DO_UPWELLING, DO_DNWELLING, &
                DO_SL_ISOTROPIC, DO_OBSERVATION_GEOMETRY,                    &
                N_SLEAVE_WFS, N_SURFACE_WFS, NSTOKES, NSTREAMS, NLAYERS,     &
                N_USER_STREAMS, LOCAL_UM_START, FOURIER_COMPONENT, IBEAM,    &
                N_USER_LEVELS, N_DIRECTIONS, WHICH_DIRECTIONS,               &
                NTOTAL, N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2,  &
                PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                     &
                UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX, &
                DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, USER_BRDF_F,       & ! New Line
                SURFACE_FACTOR, FLUX_MULTIPLIER, QUAD_WEIGHTS, QUAD_STRMWTS, &
                LSSL_DIRECT_BEAM, LSSL_USER_DIRECT_BEAM,            &
                K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,            &
                T_DELT_EIGEN, BANDMAT2, IPIVOT, SMAT2, SIPIVOT,     &
                T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,           &
                T_UTUP_EIGEN, T_UTDN_EIGEN,  HMULT_1, HMULT_2,      &
                UHOM_UPDN, UHOM_UPUP, UHOM_DNDN, UHOM_DNUP,         &
                UT_HMULT_UU, UT_HMULT_UD, UT_HMULT_DU, UT_HMULT_DD, &
                SURFACEWF_F, MINT_SURFACEWF, FLUX_SURFACEWF,        &
                STATUS_SUB, MESSAGE, TRACE_1 )

!  Exception handling

              IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
                write(CF,'(I2)')FOURIER_COMPONENT
                TRACE_2 = 'Error return from VLIDORT_LSSL_WFS, '// &
                 'Beam Called in VLIDORT_LPS_FOURIER, Fourier # '//CF
                STATUS = VLIDORT_SERIOUS
                RETURN
              ENDIF

!  End surface-leaving linearizations

            ENDIF

!@@@@@@@@@@ END Addition of SLEAVE weighting function @@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  End of surface weighting functions

          ENDIF

!  New, 28 March 2014. Linearization for BLACKBODY

!          IF ( ( DO_ATMOS_LBBF .or. DO_SURFACE_LBBF ) &
!                 .and. IBEAM.eq.1 .and. Fourier_component.eq.0 )  THEN
!           if ( n_partlayers .gt. 0 ) then
!             call vlidort_lbbf_jacobians_wpartials &
!               ( DO_ATMOS_LBBF, DO_SURFACE_LBBF, DO_THERMAL_TRANSONLY,      & ! Input
!                 DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES,              & ! Input
!                 DO_MSMODE_THERMAL, DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT,   & ! input
!                 DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,                       & ! input
!                 NSTOKES, NLAYERS, NSTREAMS, N_USER_STREAMS, N_USER_LEVELS, & ! input
!                 NMOMENTS, NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2,        & ! Input
!                 NTOTAL, N_SUPDIAG, N_SUBDIAG,  MUELLER_INDEX,              & ! input
!                 N_PARTLAYERS, PARTLAYERS_LAYERIDX,                         & ! Input
!                 PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                   & ! Input
!                 UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,                    & ! Input
!                 USER_STREAMS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,      & ! Input
!                 QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS,                  & ! input
!                 SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,               & ! input
!                 EMISSIVITY, USER_EMISSIVITY, FLUX_MULTIPLIER,              & ! input
!                 OMEGA_TOTAL, DELTAU_VERT, PARTAU_VERT, OMEGA_GREEK,        & ! Input
!                 T_DELT_DISORDS, T_UTDN_DISORDS, T_UTUP_DISORDS,            & ! input
!                 PI_XQP, PI_XQP_PRE, PI_XUP, PI_XUM, SAB, DAB,              & ! input
!                 K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,                   & ! input
!                 T_DELT_EIGEN, T_UTDN_EIGEN, T_UTUP_EIGEN,                  & ! input
!                 T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM,                  & ! Input
!                 BANDMAT2, IPIVOT, SMAT2, SIPIVOT, HMULT_1, HMULT_2,        & ! Input
!                 UT_HMULT_UU, UT_HMULT_UD, UT_HMULT_DU, UT_HMULT_DD,        & ! Input
!                 UHOM_DNDN, UHOM_DNUP, UHOM_UPDN, UHOM_UPUP,                & ! input
!                 ABBWFS_JACOBIANS, ABBWFS_FLUXES,                           & ! Output
!                 SBBWFS_JACOBIANS, SBBWFS_FLUXES,                           & ! Output
!                 STATUS_SUB, MESSAGE, TRACE_1 )                               ! Output
!         else
!            call vlidort_lbbf_jacobians_whole &
!               ( DO_ATMOS_LBBF, DO_SURFACE_LBBF, DO_THERMAL_TRANSONLY,      & ! Input
!                 DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES,              & ! Input
!                 DO_MSMODE_THERMAL, DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT,   & ! input
!                 DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,                       & ! input
!                 NSTOKES, NLAYERS, NSTREAMS, N_USER_STREAMS, N_USER_LEVELS, & ! input
!                 NMOMENTS, NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2,        & ! Input
!                 NTOTAL, N_SUPDIAG, N_SUBDIAG,  MUELLER_INDEX,              & ! input
!                 UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,                    & ! Input
!                 USER_STREAMS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,      & ! Input
!                 QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS,                  & ! input
!                 SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,               & ! input
!                 EMISSIVITY, USER_EMISSIVITY,                               & ! input
!                 FLUX_MULTIPLIER, OMEGA_TOTAL, DELTAU_VERT, OMEGA_GREEK,    & ! Input
!                 T_DELT_DISORDS, T_DELT_USERM,                              & ! input
!                 PI_XQP, PI_XQP_PRE, PI_XUP, PI_XUM, SAB, DAB,              & ! input
!                 K_REAL, K_COMPLEX, T_DELT_EIGEN, SOLA_XPOS, SOLB_XNEG,     & ! input
!                 BANDMAT2, IPIVOT, SMAT2, SIPIVOT, HMULT_1, HMULT_2,        & ! Input
!                 UHOM_DNDN, UHOM_DNUP, UHOM_UPDN, UHOM_UPUP,                & ! input
!                 ABBWFS_JACOBIANS, ABBWFS_FLUXES,                           & ! Output
!                 SBBWFS_JACOBIANS, SBBWFS_FLUXES,                           & ! Output
!                 STATUS_SUB, MESSAGE, TRACE_1 )                               ! Output
!         endif

!  Error handling

!         IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
!            write(CF,'(I2)')FOURIER_COMPONENT
!            TRACE_2 = 'Error return from vlidort_lbbf_jacobians, With Solar Term'// &
!                        'Called in VLIDORT_LPS_FOURIER_MASTER, Fourier # '//CF
!            STATUS = VLIDORT_SERIOUS
!            RETURN
!         ENDIF

!  ENd of BB Jacobians

!        ENDIF

!  Continuation point for avoiding weighting functions

 4000   CONTINUE

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
      END SUBROUTINE VLIDORT_LPS_FOURIER

      END MODULE vlidort_lps_masters

