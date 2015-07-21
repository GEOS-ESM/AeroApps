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
! #            VLIDORT_L_INPUT_MASTER (master)                  #
! #                                                             #
! #    Input-Master routine calls the following                 #
! #                                                             #
! #            (4 routines from vlidort_inputs)                 #
! #            VLIDORT_L_INIT_VARS                              #
! #            VLIDORT_L_INPUTREAD                              #
! #                                                             #
! #    These routines called by Main VLIDORT_LPS/LCS module     #
! #                                                             #
! #            VLIDORT_L_CHECK_INPUT_DIMS                       #
! #            VLIDORT_L_CHECK_INPUT                            #
! #                                                             #
! ###############################################################

      MODULE vlidort_l_inputs

      PRIVATE
      PUBLIC :: VLIDORT_L_INPUT_MASTER, &
                VLIDORT_L_CHECK_INPUT_DIMS, &
                VLIDORT_L_CHECK_INPUT

      CONTAINS

      SUBROUTINE VLIDORT_L_INPUT_MASTER ( &
        FILNAM,             & ! Input
        VLIDORT_FixIn,      & ! Outputs
        VLIDORT_ModIn,      & ! Outputs
        VLIDORT_LinFixIn,   & ! Outputs
        VLIDORT_LinModIn,   & ! Outputs
        VLIDORT_InputStatus ) ! Outputs

      USE VLIDORT_PARS
      USE VLIDORT_Inputs_def
      USE VLIDORT_LinInputs_def
      USE VLIDORT_Outputs_def
      USE VLIDORT_AUX
      USE VLIDORT_INPUTS

      IMPLICIT NONE

!  Input data filename

      CHARACTER (LEN=*), INTENT (IN) ::   FILNAM

!  Outputs

      TYPE(VLIDORT_Fixed_Inputs), INTENT (OUT)             :: &
        VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs), INTENT (OUT)          :: &
        VLIDORT_ModIn
      TYPE(VLIDORT_Fixed_LinInputs), INTENT (OUT)          :: &
        VLIDORT_LinFixIn
      TYPE(VLIDORT_Modified_LinInputs), INTENT (OUT)       :: &
        VLIDORT_LinModIn
      TYPE(VLIDORT_Input_Exception_Handling), INTENT (OUT) :: &
        VLIDORT_InputStatus

!  Local variables

      LOGICAL ::            DO_FULLRAD_MODE
      LOGICAL ::            DO_SSCORR_NADIR
      LOGICAL ::            DO_SSCORR_OUTGOING

!  New 02 Jul 2013
      LOGICAL ::            DO_FO_CALC

      LOGICAL ::            DO_SSCORR_TRUNCATION

!  New 15 March 2012
      LOGICAL ::            DO_SS_EXTERNAL

!  New 17 May 2012
      LOGICAL ::            DO_SURFACE_LEAVING
      LOGICAL ::            DO_SL_ISOTROPIC

      LOGICAL ::            DO_SSFULL
      LOGICAL ::            DO_DOUBLE_CONVTEST
      LOGICAL ::            DO_SOLAR_SOURCES
      LOGICAL ::            DO_PLANE_PARALLEL
      LOGICAL ::            DO_REFRACTIVE_GEOMETRY
      LOGICAL ::            DO_CHAPMAN_FUNCTION
      LOGICAL ::            DO_RAYLEIGH_ONLY
      LOGICAL ::            DO_DELTAM_SCALING
      LOGICAL ::            DO_SOLUTION_SAVING
      LOGICAL ::            DO_BVP_TELESCOPING
      LOGICAL ::            DO_UPWELLING
      LOGICAL ::            DO_DNWELLING
      LOGICAL ::            DO_QUAD_OUTPUT
      LOGICAL ::            DO_USER_VZANGLES
      LOGICAL ::            DO_ADDITIONAL_MVOUT
      LOGICAL ::            DO_MVOUT_ONLY
      LOGICAL ::            DO_DEBUG_WRITE
      LOGICAL ::            DO_WRITE_INPUT
      LOGICAL ::            DO_WRITE_SCENARIO
      LOGICAL ::            DO_WRITE_FOURIER
      LOGICAL ::            DO_WRITE_RESULTS
      CHARACTER (LEN=60) :: INPUT_WRITE_FILENAME
      CHARACTER (LEN=60) :: SCENARIO_WRITE_FILENAME
      CHARACTER (LEN=60) :: FOURIER_WRITE_FILENAME
      CHARACTER (LEN=60) :: RESULTS_WRITE_FILENAME

!  Order of Taylor series (including terms up to EPS^n).
!    Introduced 2/19/14 for Version 2.7
      
      INTEGER :: TAYLOR_ORDER

      INTEGER ::            NSTOKES
      INTEGER ::            NSTREAMS
      INTEGER ::            NLAYERS
      INTEGER ::            NFINELAYERS
      INTEGER ::            NGREEK_MOMENTS_INPUT
      DOUBLE PRECISION ::   VLIDORT_ACCURACY
      DOUBLE PRECISION ::   FLUX_FACTOR
      INTEGER ::            N_SZANGLES
      DOUBLE PRECISION ::   SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION ::   EARTH_RADIUS
      DOUBLE PRECISION ::   RFINDEX_PARAMETER
      DOUBLE PRECISION ::   GEOMETRY_SPECHEIGHT
      INTEGER ::            N_USER_RELAZMS
      DOUBLE PRECISION ::   USER_RELAZMS  ( MAX_USER_RELAZMS )
      INTEGER ::            N_USER_VZANGLES
      DOUBLE PRECISION ::   USER_VZANGLES ( MAX_USER_VZANGLES )
      INTEGER ::            N_USER_LEVELS
      DOUBLE PRECISION ::   USER_LEVELS ( MAX_USER_LEVELS )
      LOGICAL ::            DO_LAMBERTIAN_SURFACE
      DOUBLE PRECISION ::   LAMBERTIAN_ALBEDO
      LOGICAL ::            DO_OBSERVATION_GEOMETRY
      INTEGER ::            N_USER_OBSGEOMS
      DOUBLE PRECISION ::   USER_OBSGEOMS ( MAX_USER_OBSGEOMS, 3 )

      LOGICAL ::            DO_THERMAL_EMISSION
      INTEGER ::            N_THERMAL_COEFFS
      DOUBLE PRECISION ::   THERMAL_BB_INPUT ( 0:MAXLAYERS )
      LOGICAL ::            DO_SURFACE_EMISSION
      DOUBLE PRECISION ::   SURFBB
      LOGICAL ::            DO_THERMAL_TRANSONLY

      LOGICAL ::            DO_SPECIALIST_OPTION_1
      LOGICAL ::            DO_SPECIALIST_OPTION_2
      LOGICAL ::            DO_SPECIALIST_OPTION_3
      LOGICAL ::            DO_TOA_CONTRIBS

      LOGICAL ::            DO_SIMULATION_ONLY
      LOGICAL ::            DO_LINEARIZATION
      LOGICAL ::            DO_PROFILE_LINEARIZATION
      LOGICAL ::            DO_COLUMN_LINEARIZATION
      LOGICAL ::            DO_ATMOS_LINEARIZATION
!      LOGICAL ::            DO_LTE_LINEARIZATION  ! 2p6 variable replaced
      LOGICAL ::            DO_ATMOS_LBBF          ! 2p7 variable new

      INTEGER ::            N_TOTALPROFILE_WFS
      INTEGER ::            N_TOTALCOLUMN_WFS
      LOGICAL ::            LAYER_VARY_FLAG  ( MAXLAYERS )
      INTEGER ::            LAYER_VARY_NUMBER ( MAXLAYERS )
      LOGICAL ::            DO_SURFACE_LINEARIZATION
!      LOGICAL ::            DO_SURFBB_LINEARIZATION ! 2p6 variable replaced
      LOGICAL ::            DO_SURFACE_LBBF          ! 2p7 variable new
      LOGICAL ::            DO_SLEAVE_WFS
      INTEGER ::            N_SURFACE_WFS
      INTEGER ::            N_SLEAVE_WFS
      CHARACTER (LEN=31) :: PROFILEWF_NAMES ( MAX_ATMOSWFS )
      CHARACTER (LEN=31) :: COLUMNWF_NAMES  ( MAX_ATMOSWFS )

      INTEGER ::             STATUS
      INTEGER ::             NMESSAGES
      CHARACTER (LEN=120) :: MESSAGES (0:MAX_MESSAGES)
      CHARACTER (LEN=120) :: ACTIONS (0:MAX_MESSAGES)

      INTEGER ::             STATUS_SUB, FILUNIT

!  initialize status

      STATUS = VLIDORT_SUCCESS

      MESSAGES(1:MAX_MESSAGES) = ' '
      ACTIONS (1:MAX_MESSAGES) = ' '

      NMESSAGES       = 0
      MESSAGES(0)     = 'Successful Read of VLIDORT Input file'
      ACTIONS(0)      = 'No Action required for this Task'

!  Initialize standard variables
!    Surface-leaving terms added 17 May 2012

      CALL VLIDORT_INIT_CONTROL_VARS ( &
        DO_FULLRAD_MODE, DO_FO_CALC, DO_SSCORR_NADIR, &
        DO_SSCORR_OUTGOING, DO_SSCORR_TRUNCATION, &
        DO_SS_EXTERNAL, DO_SSFULL, DO_DOUBLE_CONVTEST, &
        DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, &
        DO_CHAPMAN_FUNCTION, DO_RAYLEIGH_ONLY, &
        DO_DELTAM_SCALING, DO_SOLUTION_SAVING, &
        DO_BVP_TELESCOPING, DO_UPWELLING, &
        DO_DNWELLING, DO_QUAD_OUTPUT, &
        DO_USER_VZANGLES, DO_ADDITIONAL_MVOUT, &
        DO_MVOUT_ONLY, DO_DEBUG_WRITE, &
        DO_WRITE_INPUT, DO_WRITE_SCENARIO, &
        DO_WRITE_FOURIER, DO_WRITE_RESULTS, &
        DO_LAMBERTIAN_SURFACE, DO_OBSERVATION_GEOMETRY, &
        DO_SPECIALIST_OPTION_1, DO_SPECIALIST_OPTION_2, &
        DO_SPECIALIST_OPTION_3, DO_TOA_CONTRIBS, &
        DO_SURFACE_LEAVING, DO_SL_ISOTROPIC )

      CALL VLIDORT_INIT_THERMAL_VARS ( &
        DO_THERMAL_EMISSION, N_THERMAL_COEFFS, &
        THERMAL_BB_INPUT, DO_SURFACE_EMISSION, &
        SURFBB, DO_THERMAL_TRANSONLY )

      CALL VLIDORT_INIT_MODEL_VARS ( &
        TAYLOR_ORDER, NSTOKES, NSTREAMS, &
        NLAYERS, NFINELAYERS, &
        NGREEK_MOMENTS_INPUT, VLIDORT_ACCURACY, &
        FLUX_FACTOR, N_SZANGLES, &
        SZANGLES, EARTH_RADIUS, &
        RFINDEX_PARAMETER, GEOMETRY_SPECHEIGHT, &
        N_USER_RELAZMS, USER_RELAZMS, &
        N_USER_VZANGLES, USER_VZANGLES, &
        N_USER_LEVELS, USER_LEVELS, &
        N_USER_OBSGEOMS, USER_OBSGEOMS, &
        LAMBERTIAN_ALBEDO )

!  initialize linearization variables

      CALL VLIDORT_L_INIT_VARS ( &
        DO_SIMULATION_ONLY, DO_LINEARIZATION, &
        DO_PROFILE_LINEARIZATION, DO_COLUMN_LINEARIZATION, &
        DO_ATMOS_LINEARIZATION, DO_ATMOS_LBBF, &
        N_TOTALPROFILE_WFS, N_TOTALCOLUMN_WFS, &
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
        DO_SURFACE_LINEARIZATION, DO_SURFACE_LBBF, &
        DO_SLEAVE_WFS, &
        N_SURFACE_WFS, N_SLEAVE_WFS, &
        PROFILEWF_NAMES, COLUMNWF_NAMES )

!  Open file

      FILUNIT = VLIDORT_INUNIT
      OPEN(VLIDORT_INUNIT,FILE=FILNAM,ERR=300,STATUS='OLD')

!  Read standard inputs
!    Surface-leaving terms added 17 May 2012

      CALL VLIDORT_READ_INPUTS ( &
        DO_FULLRAD_MODE, DO_SSCORR_NADIR, &
        DO_SSCORR_OUTGOING, DO_SSCORR_TRUNCATION, &
        DO_SS_EXTERNAL, DO_SSFULL, DO_DOUBLE_CONVTEST, &
        DO_SOLAR_SOURCES, DO_PLANE_PARALLEL, &
        DO_REFRACTIVE_GEOMETRY, DO_CHAPMAN_FUNCTION, &
        DO_RAYLEIGH_ONLY, DO_DELTAM_SCALING, &
        DO_SOLUTION_SAVING, DO_BVP_TELESCOPING, &
        DO_UPWELLING, DO_DNWELLING, &
        DO_QUAD_OUTPUT, DO_USER_VZANGLES, &
        DO_SURFACE_LEAVING, DO_SL_ISOTROPIC,&
        DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY, &
        DO_OBSERVATION_GEOMETRY, &
        DO_DEBUG_WRITE, DO_WRITE_INPUT, &
        DO_WRITE_SCENARIO, DO_WRITE_FOURIER, &
        DO_WRITE_RESULTS, INPUT_WRITE_FILENAME, &
        SCENARIO_WRITE_FILENAME, FOURIER_WRITE_FILENAME, &
        RESULTS_WRITE_FILENAME, TAYLOR_ORDER, NSTOKES, &   ! New 2p7
        NSTREAMS, NLAYERS, &
        NFINELAYERS, NGREEK_MOMENTS_INPUT, &
        VLIDORT_ACCURACY, FLUX_FACTOR, &
        N_SZANGLES, SZANGLES, &
        EARTH_RADIUS, RFINDEX_PARAMETER, &
        GEOMETRY_SPECHEIGHT, N_USER_RELAZMS, &
        USER_RELAZMS, N_USER_VZANGLES, &
        USER_VZANGLES, N_USER_LEVELS, &
        USER_LEVELS, N_USER_OBSGEOMS, &
        USER_OBSGEOMS, DO_LAMBERTIAN_SURFACE, &
        LAMBERTIAN_ALBEDO, DO_THERMAL_EMISSION, &
        N_THERMAL_COEFFS, DO_SURFACE_EMISSION, &
        DO_THERMAL_TRANSONLY, &
        STATUS_SUB, NMESSAGES, MESSAGES, ACTIONS )

      IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
        STATUS = VLIDORT_SERIOUS
        VLIDORT_InputStatus%TS_STATUS_INPUTREAD = STATUS
        VLIDORT_InputStatus%TS_NINPUTMESSAGES   = NMESSAGES
        VLIDORT_InputStatus%TS_INPUTMESSAGES    = MESSAGES
        VLIDORT_InputStatus%TS_INPUTACTIONS     = ACTIONS
        CLOSE(FILUNIT)
        RETURN
      ENDIF

!  read additional inputs for linearization

      CALL VLIDORT_L_INPUTREAD ( &
        NLAYERS, DO_THERMAL_TRANSONLY, &
        DO_SIMULATION_ONLY, DO_LINEARIZATION, &
        DO_PROFILE_LINEARIZATION, DO_COLUMN_LINEARIZATION, &
        DO_ATMOS_LINEARIZATION, DO_ATMOS_LBBF, &
        N_TOTALPROFILE_WFS, N_TOTALCOLUMN_WFS, &
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
        DO_SURFACE_LINEARIZATION, DO_SURFACE_LBBF, &
        PROFILEWF_NAMES, COLUMNWF_NAMES, &
        STATUS_SUB, NMESSAGES, MESSAGES, ACTIONS )

      IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
        STATUS = VLIDORT_SERIOUS
        VLIDORT_InputStatus%TS_STATUS_INPUTREAD = STATUS
        VLIDORT_InputStatus%TS_NINPUTMESSAGES   = NMESSAGES
        VLIDORT_InputStatus%TS_INPUTMESSAGES    = MESSAGES
        VLIDORT_InputStatus%TS_INPUTACTIONS     = ACTIONS
        CLOSE(FILUNIT)
        RETURN
      ENDIF

!  normal execution: Copy all variables and return
!  -----------------------------------------------

!  First close file

      CLOSE(FILUNIT)

!  Copy data to structure variables:

!  Fixed Boolean inputs

      VLIDORT_FixIn%Bool%TS_DO_FULLRAD_MODE        = DO_FULLRAD_MODE
      VLIDORT_FixIn%Bool%TS_DO_SSCORR_TRUNCATION   = DO_SSCORR_TRUNCATION
!  New 15 march 2012
      VLIDORT_FixIn%Bool%TS_DO_SS_EXTERNAL         = DO_SS_EXTERNAL
      VLIDORT_FixIn%Bool%TS_DO_SSFULL              = DO_SSFULL
      VLIDORT_FixIn%Bool%TS_DO_THERMAL_EMISSION    = DO_THERMAL_EMISSION
      VLIDORT_FixIn%Bool%TS_DO_SURFACE_EMISSION    = DO_SURFACE_EMISSION
      VLIDORT_FixIn%Bool%TS_DO_PLANE_PARALLEL      = DO_PLANE_PARALLEL
      !VLIDORT_FixIn%Bool%TS_DO_BRDF_SURFACE        = DO_BRDF_SURFACE
      VLIDORT_FixIn%Bool%TS_DO_UPWELLING           = DO_UPWELLING
      VLIDORT_FixIn%Bool%TS_DO_DNWELLING           = DO_DNWELLING
      VLIDORT_FixIn%Bool%TS_DO_QUAD_OUTPUT         = DO_QUAD_OUTPUT
      VLIDORT_FixIn%Bool%TS_DO_TOA_CONTRIBS        = DO_TOA_CONTRIBS
      VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE  = DO_LAMBERTIAN_SURFACE
      VLIDORT_FixIn%Bool%TS_DO_SPECIALIST_OPTION_1 = DO_SPECIALIST_OPTION_1
      VLIDORT_FixIn%Bool%TS_DO_SPECIALIST_OPTION_2 = DO_SPECIALIST_OPTION_2
      VLIDORT_FixIn%Bool%TS_DO_SPECIALIST_OPTION_3 = DO_SPECIALIST_OPTION_3

!  New 17 May 2012
      VLIDORT_FixIn%Bool%TS_DO_SURFACE_LEAVING     = DO_SURFACE_LEAVING
      VLIDORT_FixIn%Bool%TS_DO_SL_ISOTROPIC        = DO_SL_ISOTROPIC

!  Fixed Control inputs. New 2p7 Taylor-order.

      VLIDORT_FixIn%Cont%TS_TAYLOR_ORDER     = TAYLOR_ORDER 
      VLIDORT_FixIn%Cont%TS_NSTOKES          = NSTOKES
      VLIDORT_FixIn%Cont%TS_NSTREAMS         = NSTREAMS
      VLIDORT_FixIn%Cont%TS_NLAYERS          = NLAYERS
      VLIDORT_FixIn%Cont%TS_NFINELAYERS      = NFINELAYERS
      VLIDORT_FixIn%Cont%TS_N_THERMAL_COEFFS = N_THERMAL_COEFFS
      VLIDORT_FixIn%Cont%TS_VLIDORT_ACCURACY = VLIDORT_ACCURACY

!  Fixed Beam inputs

      VLIDORT_FixIn%SunRays%TS_FLUX_FACTOR = FLUX_FACTOR

!  Fixed User Value inputs

      VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS       = N_USER_LEVELS

!  Fixed Chapman Function inputs

      !VLIDORT_FixIn%Chapman%TS_HEIGHT_GRID       = HEIGHT_GRID
      !VLIDORT_FixIn%Chapman%TS_PRESSURE_GRID     = PRESSURE_GRID
      !VLIDORT_FixIn%Chapman%TS_TEMPERATURE_GRID  = TEMPERATURE_GRID
      !VLIDORT_FixIn%Chapman%TS_FINEGRID          = FINEGRID
      VLIDORT_FixIn%Chapman%TS_RFINDEX_PARAMETER = RFINDEX_PARAMETER

!  Fixed Optical inputs

      !VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT    = DELTAU_VERT_INPUT
      !VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT = GREEKMAT_TOTAL_INPUT
      VLIDORT_FixIn%Optical%TS_THERMAL_BB_INPUT     = THERMAL_BB_INPUT
      VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO    = LAMBERTIAN_ALBEDO
      VLIDORT_FixIn%Optical%TS_SURFACE_BB_INPUT     = SURFBB

!  Fixed Write inputs

      VLIDORT_FixIn%Write%TS_DO_DEBUG_WRITE           = DO_DEBUG_WRITE

      VLIDORT_FixIn%Write%TS_DO_WRITE_INPUT           = DO_WRITE_INPUT
      VLIDORT_FixIn%Write%TS_INPUT_WRITE_FILENAME     = INPUT_WRITE_FILENAME

      VLIDORT_FixIn%Write%TS_DO_WRITE_SCENARIO        = DO_WRITE_SCENARIO
      VLIDORT_FixIn%Write%TS_SCENARIO_WRITE_FILENAME  = SCENARIO_WRITE_FILENAME

      VLIDORT_FixIn%Write%TS_DO_WRITE_FOURIER         = DO_WRITE_FOURIER
      VLIDORT_FixIn%Write%TS_FOURIER_WRITE_FILENAME   = FOURIER_WRITE_FILENAME

      VLIDORT_FixIn%Write%TS_DO_WRITE_RESULTS         = DO_WRITE_RESULTS
      VLIDORT_FixIn%Write%TS_RESULTS_WRITE_FILENAME   = RESULTS_WRITE_FILENAME

!  Modified Boolean inputs

      VLIDORT_ModIn%MBool%TS_DO_SSCORR_NADIR         = DO_SSCORR_NADIR
      VLIDORT_ModIn%MBool%TS_DO_SSCORR_OUTGOING      = DO_SSCORR_OUTGOING
      VLIDORT_ModIn%MBool%TS_DO_FO_CALC              = DO_FO_CALC

      VLIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST      = DO_DOUBLE_CONVTEST
      VLIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES        = DO_SOLAR_SOURCES

      VLIDORT_ModIn%MBool%TS_DO_REFRACTIVE_GEOMETRY  = DO_REFRACTIVE_GEOMETRY
      VLIDORT_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION     = DO_CHAPMAN_FUNCTION

      VLIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY        = DO_RAYLEIGH_ONLY
      !VLIDORT_ModIn%MBool%TS_DO_ISOTROPIC_ONLY       = DO_ISOTROPIC_ONLY
      !VLIDORT_ModIn%MBool%TS_DO_NO_AZIMUTH           = DO_NO_AZIMUTH
      !VLIDORT_ModIn%MBool%TS_DO_ALL_FOURIER          = DO_ALL_FOURIER

      VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING       = DO_DELTAM_SCALING

      VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING      = DO_SOLUTION_SAVING
      VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING      = DO_BVP_TELESCOPING

      !VLIDORT_ModIn%MBool%TS_DO_USER_STREAMS         = DO_USER_STREAMS
      VLIDORT_ModIn%MBool%TS_DO_USER_VZANGLES        = DO_USER_VZANGLES

      VLIDORT_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT     = DO_ADDITIONAL_MVOUT
      VLIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY           = DO_MVOUT_ONLY

      VLIDORT_ModIn%MBool%TS_DO_THERMAL_TRANSONLY    = DO_THERMAL_TRANSONLY

      VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY = DO_OBSERVATION_GEOMETRY

!  Modified Control inputs

      VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT = NGREEK_MOMENTS_INPUT

!  Modified Beam inputs

      !VLIDORT_ModIn%MSunRays%TS_NBEAMS      = NBEAMS
      VLIDORT_ModIn%MSunRays%TS_N_SZANGLES  = N_SZANGLES
      !VLIDORT_ModIn%MSunRays%TS_BEAM_SZAS   = BEAM_SZAS
      VLIDORT_ModIn%MSunRays%TS_SZANGLES    = SZANGLES

!  Modified User Value inputs

      VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS      = N_USER_RELAZMS
      VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS        = USER_RELAZMS

      !VLIDORT_ModIn%MUserVal%TS_N_USER_STREAMS      = N_USER_STREAMS
      VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES     = N_USER_VZANGLES
      !VLIDORT_ModIn%MUserVal%TS_USER_ANGLES_INPUT   = USER_ANGLES
      VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT = USER_VZANGLES

      VLIDORT_ModIn%MUserVal%TS_USER_LEVELS         = USER_LEVELS

      VLIDORT_ModIn%MUserVal%TS_GEOMETRY_SPECHEIGHT = GEOMETRY_SPECHEIGHT

      VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS     = N_USER_OBSGEOMS
      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT = USER_OBSGEOMS

!  Modified Chapman Function inputs

      VLIDORT_ModIn%MChapman%TS_EARTH_RADIUS      = EARTH_RADIUS

!  Modified Optical inputs

      !VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT    = OMEGA_TOTAL_INPUT

!  Fixed Linearized Control inputs

      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG    = &
                               LAYER_VARY_FLAG
      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER  = &
                               LAYER_VARY_NUMBER

      VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS  = &
                               N_TOTALCOLUMN_WFS
      VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS = &
                               N_TOTALPROFILE_WFS
      VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS      = &
                               N_SURFACE_WFS
      VLIDORT_LinFixIn%Cont%TS_N_SLEAVE_WFS       = &
                               N_SLEAVE_WFS

      VLIDORT_LinFixIn%Cont%TS_COLUMNWF_NAMES  = &
                               COLUMNWF_NAMES
      VLIDORT_LinFixIn%Cont%TS_PROFILEWF_NAMES = &
                               PROFILEWF_NAMES

!  Modified linearized control inputs
!  (First three were formerely in Fixed Lin Cont)

      VLIDORT_LinModIn%MCont%TS_DO_SIMULATION_ONLY   = &
                                DO_SIMULATION_ONLY

!  replaced 2p6 variables
!      VLIDORT_LinFixIn%Cont%TS_DO_LTE_LINEARIZATION    = DO_LTE_LINEARIZATION
!      VLIDORT_LinFixIn%Cont%TS_DO_SURFBB_LINEARIZATION = DO_SURFBB_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LBBF   = DO_ATMOS_LBBF
      VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LBBF = DO_SURFACE_LBBF

      VLIDORT_LinModIn%MCont%TS_DO_SLEAVE_WFS   = DO_SLEAVE_WFS

      VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION  = &
                                DO_COLUMN_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION = &
                                DO_PROFILE_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION   = &
                                DO_ATMOS_LINEARIZATION

      VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION = &
                                DO_SURFACE_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_LINEARIZATION         = &
                                DO_LINEARIZATION

!  Exception handling

      VLIDORT_InputStatus%TS_STATUS_INPUTREAD = STATUS
      VLIDORT_InputStatus%TS_NINPUTMESSAGES   = NMESSAGES
      VLIDORT_InputStatus%TS_INPUTMESSAGES    = MESSAGES
      VLIDORT_InputStatus%TS_INPUTACTIONS     = ACTIONS

!  Return

      RETURN

!  Open file error
!  ---------------

300   CONTINUE
      STATUS = VLIDORT_SERIOUS
      NMESSAGES = NMESSAGES + 1
      MESSAGES(NMESSAGES) = &
             'openfile failure for '//FILNAM(1:LEN_STRING(FILNAM))
      ACTIONS(NMESSAGES)  = 'Find the Right File!!'
      CLOSE(FILUNIT)

      VLIDORT_InputStatus%TS_STATUS_INPUTREAD = STATUS
      VLIDORT_InputStatus%TS_NINPUTMESSAGES   = NMESSAGES
      VLIDORT_InputStatus%TS_INPUTMESSAGES    = MESSAGES
      VLIDORT_InputStatus%TS_INPUTACTIONS     = ACTIONS

!  Finish

      END SUBROUTINE VLIDORT_L_INPUT_MASTER

!

      SUBROUTINE VLIDORT_L_INIT_VARS ( &
        DO_SIMULATION_ONLY, DO_LINEARIZATION, &
        DO_PROFILE_LINEARIZATION, DO_COLUMN_LINEARIZATION, &
        DO_ATMOS_LINEARIZATION, DO_ATMOS_LBBF, &
        N_TOTALPROFILE_WFS, N_TOTALCOLUMN_WFS, &
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
        DO_SURFACE_LINEARIZATION, DO_SURFACE_LBBF, &
        DO_SLEAVE_WFS, &
        N_SURFACE_WFS, N_SLEAVE_WFS, &
        PROFILEWF_NAMES, COLUMNWF_NAMES )

      USE VLIDORT_PARS

      IMPLICIT NONE

!mick fix - added DO_LTE_LINEARIZATION, PROFILEWF_NAMES, COLUMNWF_NAMES
! Rob Fix Version 2.7, replaced DO_LTE_LINEARIZATION and DO_SURFBB_LINEARIZATION
!                      with DO_ATMOS_LBBF, DO_SURFACE_LBBF
                            
      LOGICAL, INTENT (OUT) :: DO_SIMULATION_ONLY
      LOGICAL, INTENT (OUT) :: DO_LINEARIZATION
      LOGICAL, INTENT (OUT) :: DO_PROFILE_LINEARIZATION
      LOGICAL, INTENT (OUT) :: DO_COLUMN_LINEARIZATION
      LOGICAL, INTENT (OUT) :: DO_ATMOS_LINEARIZATION
!      LOGICAL, INTENT (OUT) :: DO_LTE_LINEARIZATION  ! 2p6 variable, Replaced 
      LOGICAL, INTENT (OUT) :: DO_ATMOS_LBBF          ! 2p7 variable 
      INTEGER, INTENT (OUT) :: N_TOTALPROFILE_WFS
      INTEGER, INTENT (OUT) :: N_TOTALCOLUMN_WFS
      LOGICAL, INTENT (OUT) :: LAYER_VARY_FLAG  ( MAXLAYERS )
      INTEGER, INTENT (OUT) :: LAYER_VARY_NUMBER ( MAXLAYERS )
      LOGICAL, INTENT (OUT) :: DO_SURFACE_LINEARIZATION
!      LOGICAL, INTENT (OUT) :: DO_SURFBB_LINEARIZATION  ! 2p6 variable, Replaced 
      LOGICAL, INTENT (OUT) :: DO_SURFACE_LBBF           ! 2p7 variable
      LOGICAL, INTENT (OUT) :: DO_SLEAVE_WFS
      INTEGER, INTENT (OUT) :: N_SURFACE_WFS
      INTEGER, INTENT (OUT) :: N_SLEAVE_WFS
      CHARACTER (LEN=*), INTENT (OUT) :: PROFILEWF_NAMES ( MAX_ATMOSWFS )
      CHARACTER (LEN=*), INTENT (OUT) :: COLUMNWF_NAMES  ( MAX_ATMOSWFS )

!  local variables

      INTEGER :: N

!  initialise overall control inputs

      DO_SIMULATION_ONLY        = .FALSE.
      DO_LINEARIZATION          = .FALSE.

!  initialise atmospheric weighting function control inputs

      DO_PROFILE_LINEARIZATION  = .FALSE.
      DO_COLUMN_LINEARIZATION   = .FALSE.
      DO_ATMOS_LINEARIZATION    = .FALSE.
!      DO_LTE_LINEARIZATION      = .FALSE.
      DO_ATMOS_LBBF             = .FALSE.

      N_TOTALCOLUMN_WFS  = 0
      N_TOTALPROFILE_WFS = 0
      N_SURFACE_WFS      = 0
      N_SLEAVE_WFS       = 0
      DO N = 1, MAXLAYERS
        LAYER_VARY_FLAG(N)   = .FALSE.
        LAYER_VARY_NUMBER(N) = 0
      ENDDO

!  initialise surface weighting function control inputs

      DO_SURFACE_LINEARIZATION = .FALSE.
!      DO_SURFBB_LINEARIZATION  = .FALSE.
      DO_SURFACE_LBBF          = .FALSE.
      DO_SLEAVE_WFS            = .FALSE.

!  initialise atmospheric weighting function names

      PROFILEWF_NAMES = ' '
      COLUMNWF_NAMES  = ' '

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_L_INIT_VARS

!

      SUBROUTINE VLIDORT_L_INPUTREAD ( &
        NLAYERS, DO_THERMAL_TRANSONLY, &
        DO_SIMULATION_ONLY, DO_LINEARIZATION, &
        DO_PROFILE_LINEARIZATION, DO_COLUMN_LINEARIZATION, &
        DO_ATMOS_LINEARIZATION, DO_ATMOS_LBBF, &
        N_TOTALPROFILE_WFS, N_TOTALCOLUMN_WFS, &
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
        DO_SURFACE_LINEARIZATION, DO_SURFACE_LBBF, &
        PROFILEWF_NAMES, COLUMNWF_NAMES, &
        STATUS, NMESSAGES, MESSAGES, ACTIONS )

      USE VLIDORT_PARS

      USE VLIDORT_AUX

      INTEGER, INTENT (IN) ::             NLAYERS
      LOGICAL, INTENT (IN) ::             DO_THERMAL_TRANSONLY

      LOGICAL, INTENT (INOUT) ::            DO_SIMULATION_ONLY
      LOGICAL, INTENT (INOUT) ::            DO_LINEARIZATION
      LOGICAL, INTENT (INOUT) ::            DO_PROFILE_LINEARIZATION
      LOGICAL, INTENT (INOUT) ::            DO_COLUMN_LINEARIZATION
      LOGICAL, INTENT (INOUT) ::            DO_ATMOS_LINEARIZATION
!      LOGICAL, INTENT (INOUT) ::            DO_LTE_LINEARIZATION ! 2p6 variable 
      LOGICAL, INTENT (INOUT) ::            DO_ATMOS_LBBF         ! 2p7 variable 
      INTEGER, INTENT (INOUT) ::            N_TOTALPROFILE_WFS
      INTEGER, INTENT (INOUT) ::            N_TOTALCOLUMN_WFS
      LOGICAL, INTENT (INOUT) ::            LAYER_VARY_FLAG  ( MAXLAYERS )
      INTEGER, INTENT (INOUT) ::            LAYER_VARY_NUMBER ( MAXLAYERS )
      LOGICAL, INTENT (INOUT) ::            DO_SURFACE_LINEARIZATION
!      LOGICAL, INTENT (INOUT) ::            DO_SURFBB_LINEARIZATION ! 2p6 var 
      LOGICAL, INTENT (INOUT) ::            DO_SURFACE_LBBF          ! 2p7 var
      CHARACTER (LEN=*), INTENT (INOUT) ::  PROFILEWF_NAMES ( MAX_ATMOSWFS )
      CHARACTER (LEN=*), INTENT (INOUT) ::  COLUMNWF_NAMES  ( MAX_ATMOSWFS )

      INTEGER, INTENT (OUT) ::            STATUS

      INTEGER, INTENT (INOUT) ::            NMESSAGES
      CHARACTER (LEN=*), INTENT (INOUT) ::  MESSAGES(0:MAX_MESSAGES)
      CHARACTER (LEN=*), INTENT (INOUT) ::  ACTIONS (0:MAX_MESSAGES)

!  local variables

      CHARACTER (LEN=9), PARAMETER :: &
        PREFIX = 'VLIDORT -'

      LOGICAL ::            ERROR
      CHARACTER (LEN=80) :: PAR_STR
      INTEGER ::            FILUNIT, I, NM

!  initialize status

      STATUS = VLIDORT_SUCCESS
      ERROR  = .FALSE.
      NM     = 0

!  These are already initialized in calling routine
!      MESSAGES(1:MAX_MESSAGES) = ' '
!      ACTIONS (1:MAX_MESSAGES) = ' '
!      NMESSAGES       = 0
!      MESSAGES(0)     = 'Successful Read of VLIDORT Input file'
!      ACTIONS(0)      = 'No Action required for this Task'

!  file unit

      FILUNIT = VLIDORT_INUNIT

!  linearization control

      PAR_STR = 'Do simulation only?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_SIMULATION_ONLY
      CALL FINDPAR_ERROR &
        ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Do atmospheric profile weighting functions?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_PROFILE_LINEARIZATION
      CALL FINDPAR_ERROR &
        ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Do atmospheric column weighting functions?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_COLUMN_LINEARIZATION
      CALL FINDPAR_ERROR &
        ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Do surface property weighting functions?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_SURFACE_LINEARIZATION
      CALL FINDPAR_ERROR &
        ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  New section, LBBF linearization flags. Version 2.7
!  No restrictions on usage of these options !!!!!

      PAR_STR = 'Atmospheric BB emission weighting functions?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_ATMOS_LBBF
      CALL FINDPAR_ERROR &
        ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Surface BB emission weighting functions?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_SURFACE_LBBF
      CALL FINDPAR_ERROR &
        ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Options currently disabled. 07 August 2014

      IF ( DO_ATMOS_LBBF .or. DO_SURFACE_LBBF ) then
         NM = NM + 1
         STATUS       = VLIDORT_SERIOUS
         MESSAGES(NM) = 'DO_ATMOS_LBBF and DO_SURFACE_LBBF not enabled 8/7/14'
         ACTIONS(NM)  = 'Turn both flags off: DO_ATMOS_LBBF, DO_SURFACE_LBBF'
         NMESSAGES    = NM
         RETURN
      ENDIF

!  LTE linearization flag. 2.6 code Replaced in Version 2.7
!      PAR_STR = 'Do LTE temperature weighting functions?'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
!           READ (FILUNIT,*,ERR=998) DO_LTE_LINEARIZATION
!      CALL FINDPAR_ERROR &
!        ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Check LTE linearization. Version 2.6 code superceded.
!      IF ( DO_LTE_LINEARIZATION ) THEN
!        IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
!          NM = NM + 1
!          MESSAGES(NM) = 'Must have THERMAL_transonly for LTE linearization'
!          ACTIONS(NM)  = 'Re-set input value'
!          STATUS  = VLIDORT_SERIOUS
!          NMESSAGES = NM
!          RETURN
!        ELSE
!          IF( DO_PROFILE_LINEARIZATION .OR. DO_COLUMN_LINEARIZATION )THEN
!            NM = NM + 1
!            MESSAGES(NM) = 'No other kind of atmospheric linearization allowed'
!            ACTIONS(NM)  = 'Re-set input value'
!            STATUS  = VLIDORT_SERIOUS
!            NMESSAGES = NM
!            RETURN
!          ENDIF
!        ENDIF
!      ENDIF

!  disabled Version 2.6 code
!      PAR_STR = 'Surface emission weighting functions?'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR))
!     &     READ (FILUNIT,*,ERR=998) DO_SURFBB_LINEARIZATION
!      CALL FINDPAR_ERROR
!     &  ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Atmospheric profile weighting function input control
!  individual layer number does not exceed whole number

      IF ( DO_PROFILE_LINEARIZATION ) THEN

        PAR_STR='Number of atmospheric profile weighting functions (total)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
          READ (FILUNIT,*,ERR=998) N_TOTALPROFILE_WFS
        CALL FINDPAR_ERROR &
        ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  This code should be prepared, not taken from file-read
!        PAR_STR = 'Atmospheric weighting functions, layer control'
!        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
!          DO I = 1, NLAYERS
!            READ (FILUNIT,*,ERR=998) LAYER_VARY_FLAG(I),
!     &                               LAYER_VARY_NUMBER(I)
!          ENDDO
!        ENDIF
!        CALL FINDPAR_ERROR
!     &  ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      ENDIF

!  atmospheric Bulk/column weighting function input control

       IF ( DO_COLUMN_LINEARIZATION ) THEN

         PAR_STR='Number of atmospheric column weighting functions (total)'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
          READ (FILUNIT,*,ERR=998) N_TOTALCOLUMN_WFS
         CALL FINDPAR_ERROR &
        ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

         DO I = 1, NLAYERS
           LAYER_VARY_FLAG(I)   = .TRUE.
           LAYER_VARY_NUMBER(I) = N_TOTALCOLUMN_WFS
         ENDDO

       ENDIF

!  weighting function names

      IF ( DO_PROFILE_LINEARIZATION ) THEN
        PAR_STR='Atmospheric profile Jacobian names (character*31)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_TOTALPROFILE_WFS
            READ (FILUNIT,'(a31)',ERR=998) PROFILEWF_NAMES(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR &
        ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ENDIF

      IF ( DO_COLUMN_LINEARIZATION ) THEN
        PAR_STR='Atmospheric column Jacobian names (character*31)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_TOTALCOLUMN_WFS
            READ (FILUNIT,'(a31)',ERR=998) COLUMNWF_NAMES(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR &
        ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ENDIF

!  set overall linearization flags

      DO_ATMOS_LINEARIZATION = ( DO_PROFILE_LINEARIZATION.OR. &
                                 DO_COLUMN_LINEARIZATION .or. &
                                 DO_ATMOS_LBBF )

      DO_LINEARIZATION = ( DO_ATMOS_LINEARIZATION   .OR. &
                           DO_SURFACE_LINEARIZATION .or. &
                           DO_SURFACE_LBBF )

!mick fix
      NMESSAGES = NM

!  normal return

      RETURN

!  line read error - abort immediately

998   CONTINUE
      NM = NM + 1
      STATUS       = VLIDORT_SERIOUS
      MESSAGES(NM) = 'Read failure for entry below String: ' &
              //PAR_STR(1:LEN_STRING(PAR_STR))
      ACTIONS(NM)  = &
                 'Re-set value: Entry wrongly formatted in Input file'
      NMESSAGES    = NM

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_L_INPUTREAD

!

      SUBROUTINE VLIDORT_L_CHECK_INPUT_DIMS_HOLD &
      ( DO_COLUMN_LINEARIZATION,  DO_PROFILE_LINEARIZATION, & ! Flags
        DO_SURFACE_LINEARIZATION, DO_SLEAVE_WFS,            & ! Flags
        N_TOTALCOLUMN_WFS, N_TOTALPROFILE_WFS,              & ! Dimension-checking
        N_SURFACE_WFS, N_SLEAVE_WFS,                        & ! Dimension-checking
        STATUS, NMESSAGES, MESSAGES, ACTIONS )                ! Status and output

!  Check input dimensions

!  module, dimensions and numbers

      USE VLIDORT_pars, only : MAX_ATMOSWFS, MAX_SURFACEWFS, &
                               MAX_SLEAVEWFS, MAX_MESSAGES, &
                               VLIDORT_SUCCESS, VLIDORT_SERIOUS

      IMPLICIT NONE

!  Subroutine inputs
!  -----------------

!  Flags

      LOGICAL  , intent(in) :: DO_COLUMN_LINEARIZATION
      LOGICAL  , intent(in) :: DO_PROFILE_LINEARIZATION
      LOGICAL  , intent(in) :: DO_SURFACE_LINEARIZATION
      LOGICAL  , intent(in) :: DO_SLEAVE_WFS

!  Numbers to be checked for Dimensioning

      INTEGER, intent(inout)   :: N_TOTALCOLUMN_WFS, N_TOTALPROFILE_WFS
      INTEGER, intent(inout)   :: N_SURFACE_WFS,     N_SLEAVE_WFS

!  Exception handling.  Message Length should be at least 120 Characters

      INTEGER      , intent(out)   :: STATUS
      INTEGER      , intent(inout) :: NMESSAGES
      CHARACTER*(*), intent(inout) :: MESSAGES(0:MAX_MESSAGES)
      CHARACTER*(*), intent(inout) :: ACTIONS (0:MAX_MESSAGES)

!  Local variables

      INTEGER :: NM

!  Initialize Exception handling

      STATUS = VLIDORT_SUCCESS
      NM = NMESSAGES

!  Check LIDORT linearized input dimensions against maximum dimensions
!  ===================================================================

      IF ( DO_COLUMN_LINEARIZATION ) THEN
        IF ( N_TOTALCOLUMN_WFS .GT. MAX_ATMOSWFS ) THEN
         NM = NM + 1
         MESSAGES(NM) = 'Bad error: Insuffient dimensioning for column WFs'
         ACTIONS(NM)  = 'Action: Decrease N_TOTALCOLUMN_WFS or increase MAX_ATMOSWFS in VLIDORT.PARS'
         STATUS = VLIDORT_SERIOUS
        ENDIF
      ENDIF

      IF ( DO_PROFILE_LINEARIZATION ) THEN
        IF ( N_TOTALPROFILE_WFS .GT. MAX_ATMOSWFS ) THEN
         NM = NM + 1
         MESSAGES(NM) = 'Bad error: Insuffient dimensioning for profile WFs'
         ACTIONS(NM)  = 'Action: Decrease N_TOTALPROFILE_WFS or increase MAX_ATMOSWFS in VLIDORT.PARS'
         STATUS = VLIDORT_SERIOUS
        ENDIF
      ENDIF

      IF ( DO_SURFACE_LINEARIZATION ) THEN
        IF ( N_SURFACE_WFS .GT. MAX_SURFACEWFS ) THEN
         NM = NM + 1
         MESSAGES(NM) = 'Bad error: Insuffient dimensioning for Surface WFs'
         ACTIONS(NM)  = 'Action: Decrease N_SURFACE_WFS or increase MAX_SURFACE_WFS in VLIDORT.PARS'
         STATUS = VLIDORT_SERIOUS
        ENDIF
      ENDIF

      IF ( DO_SLEAVE_WFS ) THEN
        IF ( N_SLEAVE_WFS .GT. MAX_SLEAVEWFS ) THEN
         NM = NM + 1
         MESSAGES(NM) =  'Bad error: Insuffient dimensioning for water surface-leaving (SLEAVE) WFs'
         ACTIONS(NM)  =  'Action: Decrease N_SLEAVE_WFS or increase MAX_SLEAVEWFS in VLIDORT.PARS'
         STATUS = VLIDORT_SERIOUS
        ENDIF
      ENDIF

!  Update NMESSAGES

      NMESSAGES = NM

!  Finish

      END SUBROUTINE VLIDORT_L_CHECK_INPUT_DIMS_HOLD

!

      SUBROUTINE VLIDORT_L_CHECK_INPUT_DIMS &
        ( VLIDORT_LinFixIn, VLIDORT_LinModIn, &
          STATUS, NMESSAGES, MESSAGES, ACTIONS )

!  Check input dimensions

!  module, dimensions and numbers

      USE VLIDORT_pars, only : MAX_ATMOSWFS, MAX_SURFACEWFS, &
                               MAX_SLEAVEWFS, MAX_MESSAGES, &
                               VLIDORT_SUCCESS, VLIDORT_SERIOUS

      USE VLIDORT_LinInputs_def

      IMPLICIT NONE

!  Subroutine inputs
!  -----------------

!  VLIDORT input structures

      TYPE(VLIDORT_Fixed_LinInputs), INTENT (IN) :: &
        VLIDORT_LinFixIn
      TYPE(VLIDORT_Modified_LinInputs), INTENT (IN) :: &
        VLIDORT_LinModIn

!  Exception handling.  Message Length should be at least 120 Characters

      INTEGER      , intent(out)   :: STATUS
      INTEGER      , intent(inout) :: NMESSAGES
      CHARACTER*(*), intent(inout) :: MESSAGES(0:MAX_MESSAGES)
      CHARACTER*(*), intent(inout) :: ACTIONS (0:MAX_MESSAGES)

!  Local variables

      INTEGER :: NM

!  Initialize Exception handling

      STATUS = VLIDORT_SUCCESS
      NM = NMESSAGES

!  Check VLIDORT linearized input dimensions
!    against maximum dimensions
!  =========================================

      IF ( VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION ) THEN
        IF ( VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS .GT. MAX_ATMOSWFS ) THEN
         NM = NM + 1
         MESSAGES(NM) = &
             'Bad error: Insuffient dimensioning for column WFs'
         ACTIONS(NM)  = &
             'Action: Decrease N_TOTALCOLUMN_WFS or increase MAX_ATMOSWFS in VLIDORT.PARS'
         STATUS = VLIDORT_SERIOUS
        ENDIF
      ENDIF

      IF ( VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION ) THEN
        IF ( VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS .GT. MAX_ATMOSWFS ) THEN
         NM = NM + 1
         MESSAGES(NM) = &
             'Bad error: Insuffient dimensioning for profile WFs'
         ACTIONS(NM)  = &
             'Action: Decrease N_TOTALPROFILE_WFS or increase MAX_ATMOSWFS in VLIDORT.PARS'
         STATUS = VLIDORT_SERIOUS
        ENDIF
      ENDIF

      IF ( VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION ) THEN
        IF ( VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS .GT. MAX_SURFACEWFS ) THEN
         NM = NM + 1
         MESSAGES(NM) = &
             'Bad error: Insuffient dimensioning for Surface WFs'
         ACTIONS(NM)  = &
             'Action: Decrease N_SURFACE_WFS or increase MAX_SURFACE_WFS in VLIDORT.PARS'
         STATUS = VLIDORT_SERIOUS
        ENDIF
      ENDIF

      IF ( VLIDORT_LinModIn%MCont%TS_DO_SLEAVE_WFS ) THEN
        IF ( VLIDORT_LinFixIn%Cont%TS_N_SLEAVE_WFS .GT. MAX_SLEAVEWFS ) THEN
         NM = NM + 1
         MESSAGES(NM) = &
             'Bad error: Insuffient dimensioning for water surface-leaving (SLEAVE) WFs'
         ACTIONS(NM)  = &
             'Action: Decrease N_SLEAVE_WFS or increase MAX_SLEAVEWFS in VLIDORT.PARS'
         STATUS = VLIDORT_SERIOUS
        ENDIF
      ENDIF

!  Update NMESSAGES

      NMESSAGES = NM

!  Finish

      END SUBROUTINE VLIDORT_L_CHECK_INPUT_DIMS

!

      SUBROUTINE VLIDORT_L_CHECK_INPUT ( &
        DO_SIMULATION_ONLY, N_SURFACE_WFS, &
        DO_PROFILE_LINEARIZATION, DO_COLUMN_LINEARIZATION, &
        DO_ATMOS_LINEARIZATION, DO_SURFACE_LINEARIZATION, &
        STATUS, NMESSAGES, MESSAGES, ACTIONS )

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::            DO_SIMULATION_ONLY
      INTEGER, INTENT (IN) ::            N_SURFACE_WFS

      LOGICAL, INTENT (INOUT) ::         DO_PROFILE_LINEARIZATION
      LOGICAL, INTENT (INOUT) ::         DO_COLUMN_LINEARIZATION
      LOGICAL, INTENT (INOUT) ::         DO_ATMOS_LINEARIZATION
      LOGICAL, INTENT (INOUT) ::         DO_SURFACE_LINEARIZATION

      INTEGER, INTENT (OUT) ::           STATUS

!  %%%%%%%%%%%% Rob Change, 28 March 2011 %%%%%%%%%%%%%%%
      INTEGER, INTENT (INOUT) ::           NMESSAGES
!   Messages and actions should be INTENT(inout)
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGES(0:MAX_MESSAGES)
      CHARACTER (LEN=*), INTENT (INOUT) :: ACTIONS (0:MAX_MESSAGES)
!      INTEGER, INTENT (OUT) ::           NMESSAGES
!      CHARACTER (LEN=*), INTENT (OUT) :: MESSAGES(0:MAX_MESSAGES)
!      CHARACTER (LEN=*), INTENT (OUT) :: ACTIONS (0:MAX_MESSAGES)
!  %%%%%%%%%%%% Rob Change, 28 March 2011 %%%%%%%%%%%%%%%

!  local variables

      INTEGER :: NM

!  Initialize output status

      STATUS = VLIDORT_SUCCESS

!  %%%%%%%%%%%% Rob Change, 28 March 2011 %%%%%%%%%%%%%%%
!   Messages and actions should be INTENT(inout)
!      MESSAGES(1:MAX_MESSAGES) = ' '
!      ACTIONS (1:MAX_MESSAGES) = ' '
!      NMESSAGES       = 0
!      MESSAGES(0)     = 'Successful Check of VLIDORT Basic Input'
!      ACTIONS(0)      = 'No Action required for this Task'
!  %%%%%%%%%%%% Rob Change, 28 March 2011 %%%%%%%%%%%%%%%

!  %%% Initialize local message count to 0, 05 March 2013
!      NM     = NMESSAGES
      NM     = 0

!  Do some checking
!  ================

!  Check formatting; should not exceed 10 entries for the number of
!  layer weighting functions. This is needed to avoid errors with
!  compilers not using variable formatting.

!      IF ( MAX_ATMOSWFS .GT. 10 ) THEN
!        NM = NM + 1
!        MESSAGES(NM) =
!     &     'Max layer variations > 10, get formatting problems!'
!        ACTIONS(NM)  = 'Re-set number to avoid this'
!        STATUS = VLIDORT_SERIOUS
!      ENDIF

!  Check something is being varied

      IF ( .NOT. DO_SIMULATION_ONLY ) THEN
        IF ( .NOT. DO_ATMOS_LINEARIZATION .AND. &
               .NOT. DO_SURFACE_LINEARIZATION ) THEN
         NM = NM + 1
         MESSAGES(NM) = &
               'Bad input: Nothing being varied or calculated'
         ACTIONS(NM)  = 'Re-set input values '
         STATUS = VLIDORT_SERIOUS
!         NMESSAGES = NM             ! @@@ Robfix 01 aug 2012, Miscount NMESSAGES
        ENDIF
      ENDIF

!  Check and fix AF flags if you want simulation only

      IF ( DO_SIMULATION_ONLY ) THEN
        IF ( DO_ATMOS_LINEARIZATION ) THEN
          NM = NM + 1
          MESSAGES(NM) = &
             'Input warning; simulation only ---> NO Atmos WFs'
          ACTIONS(NM)  = &
            'Action: VLIDORT switched off WF flag for you!'
          STATUS = VLIDORT_WARNING
          DO_PROFILE_LINEARIZATION = .FALSE.
          DO_COLUMN_LINEARIZATION  = .FALSE.
          DO_ATMOS_LINEARIZATION   = .FALSE.
        ENDIF
        IF ( DO_SURFACE_LINEARIZATION ) THEN
          NM = NM + 1
          MESSAGES(NM) = &
             'Input warning; simulation only ---> NO Surface WFs'
          ACTIONS(NM)  = &
            'Action: VLIDORT switched off WF flag for you!'
          STATUS = VLIDORT_WARNING
          DO_SURFACE_LINEARIZATION = .FALSE.
        ENDIF
      ENDIF

!  %%%%%%%%%%%% Rob Change, 28 March 2011 %%%%%%%%%%%%%%%
!  Number of messages
      NMESSAGES = NMESSAGES + NM
!  %%%%%%%%%%%% Rob Change, 28 March 2011 %%%%%%%%%%%%%%%

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_L_CHECK_INPUT

      END MODULE vlidort_l_inputs

