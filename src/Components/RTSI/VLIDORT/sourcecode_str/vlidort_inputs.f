C ###############################################################
C #                                                             #
C #                    THE VLIDORT  MODEL                       #
C #                                                             #
C #  Vectorized LInearized Discrete Ordinate Radiative Transfer #
C #  -          --         -        -        -         -        #
C #                                                             #
C ###############################################################

C ###############################################################
C #                                                             #
C #  Author :      Robert. J. D. Spurr                          #
C #                                                             #
C #  Address :      RT Solutions, Inc.                          #
C #            9 Channing Street                                #
C #             Cambridge, MA 02138, USA                        #
C #            Tel: (617) 492 1183                              #
C #                                                             #
C #  Email :      rtsolutions@verizon.net                       #
C #                                                             #
C #  Versions     :   2.0, 2.2, 2.3, 2.4, 2.4R, 2.4RT           #
C #  Release Date :   December 2005  (2.0)                      #
C #  Release Date :   March 2007     (2.2)                      #
C #  Release Date :   October 2007   (2.3)                      #
C #  Release Date :   December 2008  (2.4)                      #
C #  Release Date :   April/May 2009 (2.4R)                     #
C #  Release Date :   July 2009      (2.4RT)                    #
C #                                                             #
C #       NEW: TOTAL COLUMN JACOBIANS         (2.4)             #
C #       NEW: BPDF Land-surface KERNELS      (2.4R)            #
C #       NEW: Thermal Emission Treatment     (2.4RT)           #
C #                                                             #
C ###############################################################

C    #####################################################
C    #                                                   #
C    #   This Version of VLIDORT comes with a GNU-style  #
C    #   license. Please read the license carefully.     #
C    #                                                   #
C    #####################################################

C ###############################################################
C #                                                             #
C # Subroutines in this Module                                  #
C #                                                             #
C #            VLIDORT_INPUT_MASTER (master)                    #
C #                                                             #
C #    Input-Master routine calls the following                 #
C #                                                             #
C #            VLIDORT_INIT_CONTROL_VARS                        #
C #            VLIDORT_INIT_SURFACE_VARS                        #
C #            VLIDORT_INIT_MODEL_VARS                          #
C #            VLIDORT_INIT_THERMAL_VARS                        #
C #            VLIDORT_INPUTREAD                                #
C #                                                             #
C #    These routines are called by the Main VLIDORT module     #
C #                                                             #
C #            VLIDORT_CHECK_INPUT                              #
C #            VLIDORT_DERIVE_INPUT                             #
C #                                                             #
C ###############################################################

      SUBROUTINE VLIDORT_INPUT_MASTER ( FILNAM, EFILNAM, STATUS )

C  Read all control inputs for VLIDORT
C  -----------------------------------

C  include file of constants

      INCLUDE '../includes/VLIDORT.PARS'

C  include file of input variables (CONTROL section)

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'

C  Module arguments (input filename, error filename,output status)

      CHARACTER*(*)  FILNAM
      CHARACTER*(*)  EFILNAM
      INTEGER        STATUS

C  local variables

      INTEGER        STATUS_SUB, FILUNIT, LEN_STRING
      CHARACTER*70   MAIL
      EXTERNAL       LEN_STRING

C  initialize status

      STATUS = VLIDORT_SUCCESS

C  initialize error file

      VLIDORT_ERROR_FILENAME = EFILNAM
      VLIDORT_ERROR_INIT     = .TRUE.

C  initialize variables

      CALL VLIDORT_INIT_CONTROL_VARS
      CALL VLIDORT_INIT_SURFACE_VARS
      CALL VLIDORT_INIT_THERMAL_VARS
      CALL VLIDORT_INIT_MODEL_VARS

C  Open file

      FILUNIT = VLIDORT_INUNIT
      OPEN(VLIDORT_INUNIT,FILE=FILNAM,ERR=300,STATUS='OLD')

C  read standard inputs

      CALL VLIDORT_INPUTREAD ( STATUS_SUB )
      IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
        STATUS = VLIDORT_SERIOUS
        MAIL = 'Error from VLIDORT_INPUTREAD'
        CALL VLIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
        CLOSE(FILUNIT)
        RETURN
      ENDIF

C  normal return

      CLOSE(FILUNIT)
      RETURN
      
C  Open file error

300   CONTINUE
      STATUS = VLIDORT_SERIOUS
      MAIL = 'openfile failure for '//FILNAM(1:LEN_STRING(FILNAM))
      CALL VLIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
      CLOSE(FILUNIT)

C  Finish

      RETURN
      END

C

      SUBROUTINE VLIDORT_INIT_CONTROL_VARS

C  Initialises all control inputs for VLIDORT
C  ------------------------------------------

C  include file of constants

      INCLUDE '../includes/VLIDORT.PARS'

C  include file of input variables (CONTROL section)

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'

C  mode calculations
C  -----------------

      DO_FULLRAD_MODE    = .FALSE.

C  Nadir single scatter correction (renamed).
C  Outgoing sphericity correction introduced, 31 January 2007.

      DO_SSCORR_NADIR    = .FALSE.
      DO_SSCORR_OUTGOING = .FALSE.

C  Additional deltam scaling for the single scatter corrections
C    Code added by R. Spurr, 07 September 2007

      DO_SSCORR_TRUNCATION = .FALSE.

C  direct beam correction only with BRDFs

      DO_DBCORRECTION    = .FALSE.

C  Full-up single scattering calculation

      DO_SSFULL = .FALSE.

C  Convergence testing twice

      DO_DOUBLE_CONVTEST = .FALSE.

C  solar beam options
C  ------------------

      DO_PLANE_PARALLEL      = .FALSE.
      DO_REFRACTIVE_GEOMETRY = .FALSE.
      DO_CHAPMAN_FUNCTION    = .FALSE.

C  special options
C  ---------------

      DO_RAYLEIGH_ONLY  = .FALSE.

C  do no azimuth relegated to a Bookkeeping variable
C      DO_NO_AZIMUTH      = .FALSE.
C  Isotropic only now been removed
C      DO_ISOTROPIC_ONLY = .FALSE.

C  Performance enhancements
C  ========================

C  delta-M scaling now active version 2.

      DO_DELTAM_SCALING  = .FALSE.

C  New flags for version 2. RTSolutions 4/11/05

      DO_SOLUTION_SAVING = .FALSE.
      DO_BVP_TELESCOPING = .FALSE.

C  write options
C  =============

      DO_DEBUG_WRITE       = .FALSE.
      DO_WRITE_INPUT       = .FALSE.
      DO_WRITE_RESULTS     = .FALSE.
      DO_WRITE_SCENARIO    = .FALSE.
      DO_WRITE_FOURIER     = .FALSE.

C  User-output options
C  -------------------

      DO_ADDITIONAL_MVOUT = .FALSE.
      DO_MVOUT_ONLY       = .FALSE.

      DO_USER_VZANGLES = .FALSE.

      DO_UPWELLING = .FALSE.
      DO_DNWELLING = .FALSE.

C  Automatic settings in the bookkeeping routine

c      DO_DIRECT_BEAM        = .FALSE.
c      DO_CLASSICAL_SOLUTION = .FALSE.
c      DO_ALL_FOURIER     = .FALSE.

C  Specialist options. Should always be initialized here

      DO_SPECIALIST_OPTION_1 = .FALSE.
      DO_SPECIALIST_OPTION_2 = .FALSE.
      DO_SPECIALIST_OPTION_3 = .FALSE.

C  Quadrature output is a debug flag only.

      DO_QUAD_OUTPUT  = .FALSE.

C  finish

      RETURN
      END

C

      SUBROUTINE VLIDORT_INIT_SURFACE_VARS

C  include files
C  -------------

C  include file of constants

      INCLUDE '../includes/VLIDORT.PARS'

C  include files to be initialised (surface variables)

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'

C  Initialises all surface control inputs for LIDORT
C  -------------------------------------------------

C  Local variables

      INTEGER             K, L

C  Initialize the Lambertian surface control

      DO_LAMBERTIAN_SURFACE       = .FALSE.
      DO K = 1, MAX_BRDF_KERNELS
        LAMBERTIAN_KERNEL_FLAG(K) = .FALSE.
      ENDDO

C  Initialize BRDF control
C  =======================

      NSTREAMS_BRDF = 0
      N_BRDF_KERNELS = 0
      N_BRDF_STOKESSQ  = 0
      DO_SHADOW_EFFECT = .FALSE.
      DO_COXMUNK_DBMS  = .FALSE.

      DO K = 1, MAX_BRDF_KERNELS
        BRDF_FACTORS(K) = ZERO
        DO L = 1, MAX_BRDF_PARAMETERS
          BRDF_PARAMETERS(K,L) = ZERO
        ENDDO
      ENDDO

C  finish

      RETURN
      END

C

      SUBROUTINE VLIDORT_INIT_THERMAL_VARS

C  include files
C  -------------

C  include file of constants

      INCLUDE '../includes/VLIDORT.PARS'

C  include file of input variables (Surface section)

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'

C  Initialises all thermal inputs
C  ------------------------------

C  Local variables

      INTEGER             N

C  initial emissivity control
C  ==========================

      DO_THERMAL_EMISSION = .FALSE.
      N_THERMAL_COEFFS = 0
      DO N = 0, MAXLAYERS
        THERMAL_BB_INPUT(N) = ZERO
      ENDDO

      DO_SURFACE_EMISSION = .FALSE.
      SURFBB    = ZERO

C  This flag introduced 31 July 2007 in LIDORT.

      DO_THERMAL_TRANSONLY = .FALSE.

C  finish

      RETURN
      END

C

      SUBROUTINE VLIDORT_INIT_MODEL_VARS

C  Initialises all file-read model inputs for VLIDORT
C  --------------------------------------------------

C  include file of constants

      INCLUDE '../includes/VLIDORT.PARS'

C  include files to be initialised (fileread variables)

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'

C  local variables

      INTEGER             I

C  basic integer inputs

      NSTOKES  = 0
      NSTREAMS = 0
      NLAYERS  = 0
      NFINELAYERS = 0

      NGREEK_MOMENTS_INPUT = 0

C  accuracy. No more "zenith_tolerance"!

      VLIDORT_ACCURACY  = ZERO

C  FLux factor

      FLUX_FACTOR = ZERO

C  solar beam (flux factor fixed to 1.0 in Bookkeeping, 17 January 2006)

      N_SZANGLES   = 0
      DO I = 1, MAX_SZANGLES
        SZANGLES(I) = ZERO
      ENDDO

C  Pseudo-spherical

      EARTH_RADIUS      = ZERO
      RFINDEX_PARAMETER = ZERO

C  Geometry specification height
C   (New, 06 August 2007)

      GEOMETRY_SPECHEIGHT = ZERO

C  user angles and levels

      N_USER_VZANGLES = 0
      N_USER_RELAZMS  = 0
      N_USER_LEVELS   = 0

      DO I = 1, MAX_USER_VZANGLES
        USER_VZANGLES(I) = ZERO
      ENDDO

      DO I = 1, MAX_USER_RELAZMS
        USER_RELAZMS(I) = ZERO
      ENDDO

      DO I = 1, MAX_USER_LEVELS
        USER_LEVELS(I) = ZERO
      ENDDO

C finish

      RETURN
      END

C

      SUBROUTINE VLIDORT_INPUTREAD ( STATUS )

C  Read all control inputs for VLIDORT
C  -----------------------------------

C  include file of constants

      INCLUDE '../includes/VLIDORT.PARS'

C  include file to be filled out with input data

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'

C  Module arguments (output status)

      INTEGER         STATUS

C  local variables

      CHARACTER*9     PREFIX
      PARAMETER     ( PREFIX = 'VLIDORT -' )
      LOGICAL         ERROR
      CHARACTER * 80  PAR_STR
      LOGICAL         GFINDPAR
      EXTERNAL        GFINDPAR
      INTEGER         I, K, FILUNIT, LEN_STRING
      CHARACTER*70    MAIL, ACTION
      EXTERNAL        LEN_STRING

C  BRDF Kernel names (check)

C      DATA BRDF_CHECK_NAMES /
C     &     'Lambertian',
C     &     'Ross-thin ',
C     &     'Ross-thick',
C     &     'Li-sparse ',
C     &     'Li-dense  ',
C     &     'Hapke-fn  ',
C     &     'Roujean   ',
C     &     'Cox-Munk  '/

C  initialize status

      STATUS = VLIDORT_SUCCESS
      ERROR  = .FALSE.

C  file unit

      FILUNIT = VLIDORT_INUNIT

C  1. read all CONTROL variables
C  =============================

C  operation modes
C  ---------------

C  full Stokes vector calculation

      PAR_STR = 'Do full Stokes vector calculation?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR))
     &     READ (FILUNIT,*,ERR=998) DO_FULLRAD_MODE
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  Single scatter Corrections, calculations
C  ----------------------------------------

C  Nadir single scatter correction

      PAR_STR = 'Do nadir single scatter correction?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR))
     &     READ (FILUNIT,*,ERR=998) DO_SSCORR_NADIR
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  Outgoing single scatter correction

      PAR_STR = 'Do outgoing single scatter correction?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR))
     &     READ (FILUNIT,*,ERR=998) DO_SSCORR_OUTGOING
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  Full-up single scatter calculation

      PAR_STR = 'Do Full-up single scatter calculation?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR))
     &     READ (FILUNIT,*,ERR=998) DO_SSFULL
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  Direct beam correction (BRDF options only)

      PAR_STR = 'Do direct beam correction?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR))
     &     READ (FILUNIT,*,ERR=998) DO_DBCORRECTION
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  multiple scatter source function output control
C    Removed. 30 March 2007

c      PAR_STR = 'Output multiple scatter layer source functions?'
c      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
c     &     READ (FILUNIT,*,ERR=998) SAVE_LAYER_MSST
c      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  double convergence test

      PAR_STR = 'Perform double convergence test?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_DOUBLE_CONVTEST
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  Solar beam control
C  ------------------

c      PAR_STR = 'Include direct beam?'
c      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR))
c     &     READ (FILUNIT,*,ERR=998) DO_DIRECT_BEAM
c      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  Solution method (Classical only). Variable is set in DERIVE_INPUTS.

C      DO_CLASSICAL_SOLUTION = .TRUE.
c      PAR_STR = 'Use classical beam solution?'
c      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR))
c     &     READ (FILUNIT,*,ERR=998) DO_CLASSICAL_SOLUTION
c      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  Basic control

      PAR_STR = 'Include solar sources?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR))
     &     READ (FILUNIT,*,ERR=998) DO_SOLAR_SOURCES
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  Other control options for the solar sources

      IF ( DO_SOLAR_SOURCES ) THEN

C  Pseudo-spherical control

        PAR_STR = 'Plane-parallel treatment of direct beam?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR))
     &     READ (FILUNIT,*,ERR=998) DO_PLANE_PARALLEL
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  internal Chapman function calculation

        PAR_STR = 'Perform internal Chapman function calculation?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR))
     &     READ (FILUNIT,*,ERR=998) DO_CHAPMAN_FUNCTION
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  refractive atmosphere

        PAR_STR = 'Beam path with refractive atmosphere?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &      READ (FILUNIT,*,ERR=998) DO_REFRACTIVE_GEOMETRY
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  End control

      ENDIF

C  numerical control (azimuth series)
C  ----------------------------------

C  scatterers and phase function control

      PAR_STR='Rayleigh atmosphere only?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_RAYLEIGH_ONLY
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  Removed isotropic-only option, 17 January 2006.
c      PAR_STR='Isotropic atmosphere only?'
c      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
c     &     READ (FILUNIT,*,ERR=998) DO_ISOTROPIC_ONLY
c      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  No azimuth dependence (TRUE means Fourier m = 0 only )
C    Removed. Now this is a bookkeeping variable.
C    17 January 2006.
c      PAR_STR = 'No azimuth dependence in the solution?'
c      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
c     &     READ (FILUNIT,*,ERR=998) DO_NO_AZIMUTH
c      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  All possible Fourier components (2N-1). Debug only
c      PAR_STR = 'Compute all Fourier components?'
c      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
c     &     READ (FILUNIT,*,ERR=998) DO_ALL_FOURIER
c      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  performance control
C  -------------------

C  Delta-M scaling

C    Should only be set for solar beam sources

      IF ( DO_SOLAR_SOURCES ) THEN
       PAR_STR = 'Include Delta-M scaling?'
       IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_DELTAM_SCALING
       CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )
      ENDIF

C  Additional deltam scaling for either single-scatter corrections

      IF ( DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING.OR.DO_SSFULL )  THEN
        PAR_STR = 'Additional Delta-M scaling for SS correction?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR))
     &     READ (FILUNIT,*,ERR=998) DO_SSCORR_TRUNCATION
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )
      ENDIF

C  Solution saving mode.
C    New code, RJDS, RT Solutions, Inc. 4/11/05.

      PAR_STR = 'Solution saving mode?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_SOLUTION_SAVING
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  Boundary value problem (BVP) telescope mode.
C    New code, RJDS, RT Solutions, Inc. 4/11/05.

      PAR_STR = 'Boundary value problem telescoping mode?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_BVP_TELESCOPING
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  User-defined output control
C  ---------------------------

C  directional output control

      PAR_STR = 'Upwelling output?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_UPWELLING
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      PAR_STR = 'Downwelling output?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_DNWELLING
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  Stream angle and optical depth output control

c      PAR_STR = 'Include quadrature angles in output?'
c      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
c     &     READ (FILUNIT,*,ERR=998) DO_QUAD_OUTPUT
c      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      PAR_STR = 'User-defined viewing zenith angles?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &   READ (FILUNIT,*,ERR=998) DO_USER_VZANGLES
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  Old code, version 2.3 and earlier............
c      PAR_STR = 'User-defined optical depths?'
c      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
c     &   READ (FILUNIT,*,ERR=998) DO_USER_TAUS
c      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

c      PAR_STR = 'Layer boundary optical depths?'
c      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
c     &   READ (FILUNIT,*,ERR=998) DO_LBOUND_TAUS
c      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  mean-value output control

      PAR_STR = 'Generate mean value output additionally?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_ADDITIONAL_MVOUT
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      PAR_STR = 'Generate only mean value output?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_MVOUT_ONLY
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  Write control
C  -------------

C  output write flags

      PAR_STR = 'Debug write?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_DEBUG_WRITE
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      PAR_STR = 'Input control write?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_WRITE_INPUT
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      PAR_STR = 'Input scenario write?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_WRITE_SCENARIO
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      PAR_STR = 'Fourier component output write?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_WRITE_FOURIER
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      PAR_STR = 'Results write?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_WRITE_RESULTS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  output filenames

      PAR_STR = 'filename for input write'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,'(a)',ERR=998) INPUT_WRITE_FILENAME
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      PAR_STR = 'filename for scenario write'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,'(a)',ERR=998) SCENARIO_WRITE_FILENAME
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      PAR_STR = 'Fourier output filename'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,'(a)',ERR=998) FOURIER_WRITE_FILENAME
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      PAR_STR = 'filename for main output'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,'(a)',ERR=998) RESULTS_WRITE_FILENAME
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  2. Read all model variables
C  ===========================

C  Stokes/streams/layers/moments (INTEGER input)

      PAR_STR = 'Number of Stokes vector components'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) NSTOKES
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      PAR_STR = 'Number of half-space streams'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) NSTREAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      PAR_STR = 'Number of atmospheric layers'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) NLAYERS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      IF ( DO_SSCORR_OUTGOING .OR. DO_SSFULL ) THEN
        PAR_STR =
     &   'Number of fine layers (outgoing sphericity correction only)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &       READ (FILUNIT,*,ERR=998) NFINELAYERS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )
      ENDIF

      PAR_STR = 'Number of scattering matrix expansion coefficients'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) NGREEK_MOMENTS_INPUT
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  All numbers are now checked against maximum diensions

      IF ( NSTOKES .GT. MAXSTOKES ) THEN
        MAIL = 'Number of Stokes parameters > maximum dimension'
        STATUS  = VLIDORT_SERIOUS
        CALL VLIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
      ENDIF
      IF ( NSTREAMS .GT. MAXSTREAMS ) THEN
        MAIL = 'Number of half-space streams > maximum dimension'
        STATUS  = VLIDORT_SERIOUS
        CALL VLIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
      ENDIF
      IF ( NLAYERS .GT. MAXLAYERS ) THEN
        MAIL = 'Number of layers > maximum dimension'
        STATUS  = VLIDORT_SERIOUS
        CALL VLIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
      ENDIF
      IF ( DO_SSCORR_OUTGOING ) THEN
        IF ( NFINELAYERS .GT. MAXFINELAYERS ) THEN
          MAIL   = 'Number of fine layers > maximum dimension'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
          RETURN
        ENDIF
      ENDIF
      IF ( NGREEK_MOMENTS_INPUT .GT. MAXMOMENTS_INPUT) THEN
        MAIL = 'Number of Expansion coefficients > maximum dimension'
        STATUS  = VLIDORT_SERIOUS
        CALL VLIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
      ENDIF

C  accuracy input
C  --------------

      PAR_STR = 'Fourier series convergence'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) VLIDORT_ACCURACY
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  Zenith tolerance level. Now removed. 17 January 2006.

C  Flux constant. Should be set to 1 if no solar sources.
C   Formerly a Book-keeping variable set to 1 in "derive inputs"
C   Now (July 2009) allowed to vary because of thermal emission
C   Must take physical values if using solar + thermal.

      if ( DO_SOLAR_SOURCES ) THEN
        PAR_STR = 'Solar flux constant'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) FLUX_FACTOR
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )
      ELSE
        FLUX_FACTOR = 1.0d0
      ENDIF

C  Note the following possibility (not yet allowed for)
C        PAR_STR = 'TOA flux vector'
C        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
C     &     READ (FILUNIT,*,ERR=998) (FLUXVEC(I),I=1,MAXSTOKES)
C      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  Solar zenith angles
C  -------------------

C  number of Solar zenith angles

      PAR_STR = 'Number of solar zenith angles'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) N_SZANGLES
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  check not exceeding dimensioned number

      IF ( N_SZANGLES .GT. MAX_SZANGLES ) THEN
        MAIL = 'Number of solar zenith angles > maximum dimension'
        STATUS  = VLIDORT_SERIOUS
        CALL VLIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
      ENDIF

C  BOA solar zenith angle inputs

      PAR_STR = 'Solar zenith angles (degrees)'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
        DO I = 1, N_SZANGLES
          READ (FILUNIT,*,ERR=998) SZANGLES(I)
        ENDDO
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  Pseudo-spherical inputs
C  -----------------------

C (only for Chapman function calculation)

      IF ( DO_CHAPMAN_FUNCTION ) THEN

        PAR_STR = 'Earth radius (km)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &       READ (FILUNIT,*,ERR=998) EARTH_RADIUS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
          PAR_STR = 'Refractive index parameter'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &         READ (FILUNIT,*,ERR=998) RFINDEX_PARAMETER
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )
        ENDIF

      ENDIF

C  User defined output control
C  ---------------------------

C  Input Geometry specfication height

      IF ( DO_SSCORR_OUTGOING ) THEN
        PAR_STR = 'Input geometry specification height [km]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) GEOMETRY_SPECHEIGHT
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )
      ENDIF

C  1.  Azimuthal input values

C  Number of azimuths

      PAR_STR = 'Number of user-defined relative azimuth angles'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) N_USER_RELAZMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  Check dimensioning for number of azimuths

      IF ( N_USER_RELAZMS .GT. MAX_USER_RELAZMS ) THEN
        MAIL = 'Number of relative azimuths > maximum dimension'
        STATUS  = VLIDORT_SERIOUS
        CALL VLIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
      ENDIF

C  Read in azimuths

      PAR_STR = 'User-defined relative azimuth angles'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
        DO I = 1, N_USER_RELAZMS
          READ (FILUNIT,*,ERR=998) USER_RELAZMS(I)
        ENDDO
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  2. User defined viewing zenith angles (should be positive)
C       ( Number, check dimension, values)

      IF ( DO_USER_VZANGLES ) THEN

        PAR_STR = 'Number of user-defined stream angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) N_USER_VZANGLES
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        IF ( N_USER_VZANGLES .GT. MAX_USER_VZANGLES ) THEN
          MAIL = 'Number of user streams > maximum dimension'
          STATUS  = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
        ENDIF

        PAR_STR = 'User-defined stream angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_USER_VZANGLES
            READ (FILUNIT,*,ERR=998) USER_VZANGLES(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      ENDIF

C  3. User defined level output   *** NEW SECTION ***
C  ============================

C  This is designed to ensure that output is not related to optical
C  depth (which is wavelength-dependent); we use a height-based system.

C  User defined boundary (whole layer) and off-boundary (partial layer)
C  output choices are specified as follows.
C       USER_LEVELS(1) = 0.0    Top of the first layer
C       USER_LEVELS(2) = 1.0    Bottom of the first layer
C       USER_LEVELS(3) = 17.49  Output is in Layer 18, at a distance of
C                               0.49 of the way down from top (in height)

C    -- Number is checked now, to see if exceeds dimensioning

      PAR_STR = 'Number of user-defined vertical output levels'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &        READ (FILUNIT,*,ERR=998) N_USER_LEVELS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      IF ( N_USER_LEVELS .GT. MAX_USER_LEVELS ) THEN
        MAIL=
     &    'Number of user vertical output levels > maximum dimension'
        ACTION = ' Re-set input value'
        STATUS = VLIDORT_SERIOUS
        CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        RETURN
      ENDIF

      PAR_STR = 'User-defined vertical output levels'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
        DO I = 1, N_USER_LEVELS
          READ (FILUNIT,*,ERR=998) USER_LEVELS(I)
        ENDDO
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  4. read all surface variables
C  =============================

C  Lambertian surface

      PAR_STR = 'Lambertian surface?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_LAMBERTIAN_SURFACE
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      IF ( DO_LAMBERTIAN_SURFACE ) THEN
        PAR_STR = 'Lambertian (isotropic) input albedo'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &       READ (FILUNIT,*,ERR=998) LAMBERTIAN_ALBEDO
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )
      ENDIF

C  BRDF input
C  ----------

      IF ( .NOT. DO_LAMBERTIAN_SURFACE ) THEN

C  number of kernels, check this value

        PAR_STR = 'Number of bidirectional reflectance kernels'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &       READ (FILUNIT,*,ERR=998) N_BRDF_KERNELS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        IF ( N_BRDF_KERNELS .GT. MAX_BRDF_KERNELS ) THEN
          MAIL = 'Number of BRDF kernels > maximum dimension'
          STATUS  = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
          RETURN
        ENDIF

C  number of BRDF azimuth streams, check this value

        PAR_STR = 'Number of bidirectional reflectance streams'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &       READ (FILUNIT,*,ERR=998) NSTREAMS_BRDF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        IF ( NSTREAMS_BRDF .GT. MAXSTREAMS_BRDF ) THEN
          MAIL = 'Number of BRDF streams > maximum dimension'
          STATUS  = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
          RETURN
        ENDIF

C  Main kernel input

        PAR_STR = 
     & 'Kernel names, indices, amplitudes, # parameters, parameters'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_BRDF_KERNELS
            READ (FILUNIT,56,ERR=998)
     &         BRDF_NAMES(I), WHICH_BRDF(I), BRDF_FACTORS(I),
     &        N_BRDF_PARAMETERS(I),(BRDF_PARAMETERS(I,K),K=1,3)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )
 56     FORMAT( A10, I2, F6.2, I2, 3F12.6 )

C  Shadowing input (for Cox-Munk types)

        DO I = 1, N_BRDF_KERNELS
          IF ( BRDF_NAMES(I) .EQ. 'Cox-Munk  ' .OR.
     &         BRDF_NAMES(I) .EQ. 'GissCoxMnk' ) THEN
           PAR_STR = 'Do shadow effect for Cox-Munk kernels?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
            READ (FILUNIT,*,ERR=998)DO_SHADOW_EFFECT
           ENDIF
          ENDIF
        ENDDO
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  Multiple reflectance DB correction (for Cox-Munk types)

        DO I = 1, N_BRDF_KERNELS
         IF ( BRDF_NAMES(I) .EQ. 'Cox-Munk  ' .OR.
     &        BRDF_NAMES(I) .EQ. 'GissCoxMnk' ) THEN
           PAR_STR = 'Do multiple reflectance for Cox-Munk kernels?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
            READ (FILUNIT,*,ERR=998)DO_COXMUNK_DBMS
           ENDIF
         ENDIF
        ENDDO
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      ENDIF

C  5. read all thermal emission input variables
C  ============================================

C  Thermal controls

      PAR_STR = 'Do thermal emission?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_THERMAL_EMISSION
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  Thermal control, transmittance only

      IF ( DO_THERMAL_EMISSION ) THEN
        PAR_STR = 'Do thermal emission, transmittance only?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_THERMAL_TRANSONLY
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )
      ENDIF

C  Number of coefficients (includes a dimensioning check)

      IF ( DO_THERMAL_EMISSION ) THEN

        PAR_STR = 'Number of thermal coefficients'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) N_THERMAL_COEFFS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        IF ( N_THERMAL_COEFFS .GT. MAX_THERMAL_COEFFS ) THEN
          MAIL= 'Number of thermal coefficients > maximum dimension'
          ACTION = ' Re-set input value'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          RETURN
        ENDIF

      ENDIF

C  surface emission control

      PAR_STR = 'Do Surface emission?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_SURFACE_EMISSION
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  normal return

      RETURN

C  line read error - abort immediately

998   CONTINUE
      STATUS = VLIDORT_SERIOUS
      MAIL = 'read failure for '//PAR_STR(1:LEN_STRING(PAR_STR))
      CALL VLIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )

C  Finish

      RETURN
      END

c^L

      SUBROUTINE VLIDORT_CHECK_INPUT (STATUS)

C  check inputs (both file-read and derived)

      INCLUDE '../includes/VLIDORT.PARS'

C  include file with input variables to be checked out

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'

C  Some bookkeeping is assigned nere

      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  Module output

      INTEGER          STATUS

C  local variables

      CHARACTER*70     MAIL
      CHARACTER*70     ACTION

      INTEGER          I, K, L, N, UTA, NSTART
      INTEGER          NALLSTREAMS
      INTEGER          INDEX_ANGLES ( MAX_USER_STREAMS )
      DOUBLE PRECISION XT, ALL_ANGLES ( MAX_USER_STREAMS )
      CHARACTER*3      C3
      CHARACTER*2      C2
      LOGICAL          LOOP

C  BRDF Kernel names (check).
C   Add 2 New Land BRDFs (Rondeaux-Herman and Breon), 01 July 2008
C   BPDF kernels added for Version 2.4R.

      CHARACTER*10 BRDF_CHECK_NAMES ( MAXBRDF_IDX )
      DATA BRDF_CHECK_NAMES /
     &     'Lambertian',
     &     'Ross-thin ',
     &     'Ross-thick',
     &     'Li-sparse ',
     &     'Li-dense  ',
     &     'Hapke     ',
     &     'Roujean   ',
     &     'Rahman    ',  
     &     'Cox-Munk  ',  
     &     'GissCoxMnk',  
     &     'RondHerman',  
     &     'Breon     ',  
     &     'BPDF08Veg ',  
     &     'BPDF08Soil',
     &     'BPDF2009  '/  

C  Initialize output status

      STATUS = VLIDORT_SUCCESS

C  Automatic input
C    Flux factor set to unity. (21 December 2005)

      DO_ALL_FOURIER        = .FALSE.
      DO_CLASSICAL_SOLUTION = .TRUE.
      DO_DIRECT_BEAM        = .TRUE.
!      FLUX_FACTOR           = ONE

C  Check top level options Solar/THermal, set warnings
C  ===================================================

C  Check thermal or Solar sources present

      IF ( .NOT.DO_SOLAR_SOURCES.AND..NOT.DO_THERMAL_EMISSION ) THEN
        MAIL = 'Bad input: No solar or thermal sources'
        ACTION = 'Abort: must set one of the source flags!'
        STATUS = VLIDORT_SERIOUS
        CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
      ENDIF

C  Switch off several flags with thermal-only option
C    Set default, regardless of whether solar sources are on.

      IF ( .NOT.DO_SOLAR_SOURCES.AND.DO_THERMAL_EMISSION ) THEN
       IF ( DO_SSCORR_NADIR ) THEN
         MAIL = 'Switch off SS correction, not needed for thermal-only'
         ACTION = 'Warning: SS correction flag turned off internally'
         STATUS = VLIDORT_WARNING
         DO_SSCORR_NADIR = .FALSE.
         CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
       IF ( DO_DBCORRECTION ) THEN
         MAIL = 'Switch off DB correction, not needed for thermal-only'
         ACTION = 'Warning: DB correction flag turned off internally'
         STATUS = VLIDORT_WARNING
         DO_DBCORRECTION = .FALSE.
         CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
      ENDIF

C  Set number of sources NBEAMS to 1 for the thermal-only default

      IF ( .NOT.DO_SOLAR_SOURCES.AND.DO_THERMAL_EMISSION ) THEN
        IF ( N_SZANGLES .NE. 1 ) THEN
         MAIL = 'Bad input: N_Szas set to 1 for thermal-only'
         ACTION = 'Warning: N_Szas set to 1 internally'
         STATUS = VLIDORT_WARNING
         N_SZANGLES = 1
         CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
      ENDIF

C  Set number of sources NBEAMS to 1 for the thermal-only default

      IF ( .NOT.DO_SOLAR_SOURCES.AND.DO_THERMAL_EMISSION ) THEN
        IF ( N_USER_RELAZMS .NE. 1 ) THEN
         MAIL = 'Bad input: N_AZIMUTHS set to 1 for thermal-only'
         ACTION = 'Warning: N_USER_RELAZMS set to 1 internally'
         STATUS = VLIDORT_WARNING
         N_USER_RELAZMS = 1
         CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
      ENDIF

C  Check dimensions (safety first). All errors are FATAL
C  -----------------------------------------------------

C  1. Brdf dimensioning

      IF ( .NOT. DO_LAMBERTIAN_SURFACE ) THEN

C  number of kernels, check this value

        IF ( N_BRDF_KERNELS .GT. MAX_BRDF_KERNELS ) THEN
          MAIL   = 'Number of BRDF kernels > maximum dimension'
          ACTION = 'Re-set input value'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          RETURN
        ENDIF

C  number of BRDF azimuth streams, check this value

        IF ( NSTREAMS_BRDF .GT. MAXSTREAMS_BRDF ) THEN
          MAIL   = 'Number of BRDF streams > maximum dimension'
          ACTION = 'Re-set input value'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          RETURN
        ENDIF

      ENDIF

C  2. Basic layer/streams/moments

      IF ( NSTREAMS .GT. MAXSTREAMS ) THEN
        MAIL   = 'Number of half-space streams > maximum dimension'
        ACTION = 'Re-set input value'
        STATUS = VLIDORT_SERIOUS
        CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        RETURN
      ENDIF

      IF ( NLAYERS .GT. MAXLAYERS ) THEN
        MAIL   = 'Number of layers > maximum dimension'
        ACTION = 'Re-set input value'
        STATUS = VLIDORT_SERIOUS
        CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        RETURN
      ENDIF

      IF ( DO_SSCORR_OUTGOING ) THEN
        IF ( NFINELAYERS .GT. MAXFINELAYERS ) THEN
          MAIL   = 'Number of fine layers > maximum dimension'
          ACTION = 'Re-set input value'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          RETURN
        ENDIF
      ENDIF

      IF ( NGREEK_MOMENTS_INPUT .GT. MAXMOMENTS_INPUT ) THEN
        MAIL   = 'Number of Expansion Coeffs. > maximum dimension'
        ACTION = 'Re-set input value'
        STATUS = VLIDORT_SERIOUS
        CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        RETURN
      ENDIF

C  3. Basic geometry

      IF ( N_SZANGLES .GT. MAX_SZANGLES ) THEN
        MAIL   = 'Number of solar zenith angles > maximum dimension'
        ACTION = 'Re-set input value'
        STATUS = VLIDORT_SERIOUS
        CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        RETURN
      ENDIF

      IF ( DO_USER_VZANGLES ) THEN
        IF ( N_USER_VZANGLES .GT. MAX_USER_VZANGLES ) THEN
          MAIL   = 
     &    'Number of user view zenith angles > maximum dimension'
          ACTION = 'Re-set input value'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          RETURN
        ENDIF
      ENDIF

      IF ( N_USER_RELAZMS .GT. MAX_USER_RELAZMS ) THEN
        MAIL   = 
     &    'Number of user relative azimuths > maximum dimension'
        ACTION = 'Re-set input value'
        STATUS = VLIDORT_SERIOUS
        CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        RETURN
      ENDIF

C  4. optical depth output

      IF ( N_USER_LEVELS .GT. MAX_USER_LEVELS ) THEN
        MAIL=
     &    'Number of user vertical output levels > maximum dimension'
        ACTION = ' Re-set input value'
        STATUS = VLIDORT_SERIOUS
        CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        RETURN
      ENDIF

C  check inputs (both file-read and derived)
C  -----------------------------------------

C  Check Chapman function options

      IF ( DO_SOLAR_SOURCES ) THEN
       IF ( .NOT. DO_CHAPMAN_FUNCTION ) THEN
        IF ( DO_PLANE_PARALLEL ) THEN
          MAIL   = 'Chapman Function not set, plane parallel'
          ACTION = 'Warning: Chapman function set internally'
          STATUS = VLIDORT_WARNING
          DO_CHAPMAN_FUNCTION = .TRUE.
          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ELSE
          MAIL   = 'Chapman Function not set, pseudo-spherical'
          ACTION = 'Have you set the CHAPMAN_FACTORS values?'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
       ENDIF
      ENDIF

C  Check sphericity corrections....Cannot both be turned on
C    --------- New code 31 January 2007

      IF ( DO_SOLAR_SOURCES ) THEN
       IF ( DO_SSCORR_NADIR .and. DO_SSCORR_OUTGOING ) THEN
        MAIL   = 'Cannot have both single scatter corrections on'
        ACTION = 'Turn off DO_SSCORR_NADIR and/or DO_SSCORR_OUTGOING'
        STATUS = VLIDORT_SERIOUS
        CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
       ENDIF
      ENDIF

C  This check no longer applies.....        
c      IF ( DO_SSCORR_OUTGOING ) THEN
c        IF ( DO_USER_TAUS ) THEN
c          MAIL   = 'SS correction outgoing only for LBOUND_TAUS'
c          ACTION = 'Turn off DO_USER_TAUS, use only DO_LBOUND_TAUS'
c          STATUS = VLIDORT_SERIOUS
c          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
c        ENDIF
c        IF ( DO_DNWELLING ) THEN
c          MAIL   = 'SS correction outgoing only for Upwelling'
c          ACTION = 'Turn off DO_DNWELLING'
c          STATUS = VLIDORT_SERIOUS
c          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
c        ENDIF
c      ENDIF
   
C  Check beam mode operation

      IF ( DO_SOLAR_SOURCES ) THEN
       IF ( DO_PLANE_PARALLEL ) THEN
        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
         MAIL =
     &  'Bad input: plane-parallel and refractive flags both set'
         ACTION = 'Warning: turn off Refraction internally'
         STATUS = VLIDORT_WARNING
         DO_REFRACTIVE_GEOMETRY = .FALSE.
         CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
       ENDIF
      ENDIF

C  check consistency of mean value input control
C  ---------------------------------------------

      IF ( DO_ADDITIONAL_MVOUT ) THEN
        IF ( DO_MVOUT_ONLY ) THEN
          MAIL   = 'Bad input: Cannot have both mean-value flags set'
          ACTION = 'Warning: disable DO_MVOUT_ONLY flag internally'
          STATUS = VLIDORT_WARNING
          DO_MVOUT_ONLY = .FALSE.
          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
      ENDIF

C  remove section on DO_NO_AZIMUTH. 17 January 2006

      IF ( .NOT.DO_ADDITIONAL_MVOUT ) THEN
        IF ( DO_MVOUT_ONLY ) THEN
          IF ( DO_USER_VZANGLES ) THEN
            MAIL   =
     &           'Bad input: Mean-value option needs quadratures only'
            ACTION = 
     &         'Warning: DO_USER_VZANGLES flag disabled internally'
            STATUS = VLIDORT_WARNING
            DO_USER_VZANGLES = .FALSE.
            CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          ENDIF
        ENDIF
      ENDIF

C  check consistency of multiple scatter source term output control
C   ---SPecialist options. Removed 30 March 2007.

c      IF ( SAVE_LAYER_MSST ) THEN
c        IF ( .NOT. DO_LBOUND_TAUS ) THEN
c          MAIL =
c     & 'Bad input: MSCAT. source term - layer boundary TAU flag not set'
c          ACTION = 'Check DO_LBOUND_TAUS and SAVE_LAYER_MSST flags'
c          STATUS = VLIDORT_SERIOUS
c          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
c        ENDIF
c        IF ( .NOT. DO_USER_VZANGLES ) THEN
c          MAIL   =
c     &   'Bad input: MSCAT. source term - USER_VZANGLES flag not set'
c          ACTION = 'Check DO_USER_VZANGLES and SAVE_LAYER_MSST flags'
c          STATUS = VLIDORT_SERIOUS
c          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
c        ENDIF
c      ENDIF

C  Check consistency of Quadrature output flag
C    - Quadrature output only with Uncorrected Full-radiance calculation.
C    - All warnings. DO_QUAD_OUTPUT disabled internally

c      IF ( DO_QUAD_OUTPUT ) THEN
c        IF ( .NOT. DO_FULLRAD_MODE ) THEN
c          MAIL =
c     &     'Bad input: Quadrature output not with Full radiance mode'
c          ACTION = 'Warning: turn off DO_QUAD_OUTPUT internally'
c          STATUS = VLIDORT_WARNING
c          DO_QUAD_OUTPUT = .FALSE.
c          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
c        ELSE
c         IF ( DO_SOLAR_SOURCES ) THEN
c          IF ( DO_SSCORR_NADIR .OR. DO_SSCORR_OUTGOING ) THEN
c            MAIL =
c     &       'Bad input: Quadrature output not with SS correction'
c            ACTION = 'Warning: turn off DO_QUAD_OUTPUT internally'
c            STATUS = VLIDORT_WARNING
c            DO_QUAD_OUTPUT = .FALSE.
c            CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
c          ENDIF
c          IF ( CO_DBCORRECTION ) THEN
c            MAIL =
c     &       'Bad input: Quadrature output not with DB correction'
c            ACTION = 'Warning: turn off DO_QUAD_OUTPUT internally'
c            STATUS = VLIDORT_WARNING
c            DO_QUAD_OUTPUT = .FALSE.
c            CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
c          ENDIF
c        ENDIF
c       ENDIF
c      ENDIF

C  Check consistency of BVP_TELESCOPING and SOLUTION_SAVING flags
C  ---Warning. Set solution-saving internally

      IF (DO_BVP_TELESCOPING.AND..NOT.DO_SOLUTION_SAVING) THEN
        MAIL   =
     &  'Bad input: BVP telescoping -> solution saving must be set'
        ACTION = 'Warning:  Solution saveing was set internally'
        STATUS = VLIDORT_WARNING
        DO_SOLUTION_SAVING = .TRUE.
        CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
      ENDIF

C  Check consistency of Rayleigh-only and Isotropic-only cases
C   Removed 17 January 2006, Isotropic option removed.
c      IF ( DO_RAYLEIGH_ONLY .AND. DO_ISOTROPIC_ONLY ) THEN
c        MAIL   =
c     &   'Bad input: Isotropic_only & Rayleigh-only flags both set'
c        ACTION =
c     &     'Check DO_RAYLEIGH_ONLY and DO_ISOTROPIC_ONLY flags'
c        STATUS = VLIDORT_SERIOUS
c        CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
c      ENDIF

C  -----------Note the following in the scalar code -------------------

C  no Delta-M scaling with Rayleigh only
C   ---Warning. Turn off delta-M scaling.

      IF ( DO_RAYLEIGH_ONLY ) THEN
        IF ( DO_DELTAM_SCALING ) THEN
          MAIL  ='Bad input: No delta-M scaling with Rayleigh-only'
          ACTION = 'Warning: DO_DELTAM_SCALING turned off internally'
          STATUS = VLIDORT_WARNING
          DO_DELTAM_SCALING = .FALSE.
          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
      ENDIF

C  no Delta-M scaling with Isotropic only
C   ---Warning. Turn off delta-M scaling.
C   Removed 17 January 2006, Isotropic option removed.
c      IF ( DO_ISOTROPIC_ONLY ) THEN
c        IF ( DO_DELTAM_SCALING ) THEN
c          MAIL  ='Bad input: No delta-M scaling with Isotropic-only'
c          ACTION = 'Warning: DO_DELTAM_SCALING turned off internally'
c          STATUS = VLIDORT_WARNING
c          DO_DELTAM_SCALING = .FALSE.
c          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
c        ENDIF
c      ENDIF

C---------------------------------------------------------------------

C  check optical depth flags

c      IF ( DO_USER_TAUS .AND. DO_LBOUND_TAUS ) THEN
c        MAIL   = 'Bad input: both optical depth flags cannot be set'
c        ACTION = 'Check DO_USER_TAUS and DO_LBOUND_TAUS flags'
c        STATUS = VLIDORT_SERIOUS
c        CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
c      ELSE IF ( .NOT. DO_USER_TAUS .AND. .NOT. DO_LBOUND_TAUS ) THEN
c        MAIL   = 'Bad input: Neither optical depth flags are set'
c        ACTION = 'Check DO_USER_TAUS and DO_LBOUND_TAUS flags'
c        STATUS = VLIDORT_SERIOUS
c        CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
c      ENDIF
        
C  check directional input

      IF ( .NOT.DO_UPWELLING .AND. .NOT. DO_DNWELLING ) THEN
        MAIL   = 'Bad input: no directional input is set'
        ACTION = 'Check DO_UPWELLING & DO_DNWELLING: one must be set!'
        STATUS = VLIDORT_SERIOUS
        CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
      ENDIF

C  check number of input expansion coefficient  moments (non-Rayleigh)
C  ====================================================

C   Isotropic part Removed 17 January 2006, Isotropic option removed.

      IF ( .NOT.DO_RAYLEIGH_ONLY .AND. .NOT. DO_SSFULL ) THEN
        IF ( DO_DELTAM_SCALING ) THEN
          IF ( NGREEK_MOMENTS_INPUT.LT.2*NSTREAMS ) THEN
            MAIL =
     &      'Bad input: Fewer than 2N expansion moments with delta-M'
            ACTION = 
     &       'Warning: Re-set NGREEK_MOMENTS_INPUT to 2N internally'
            STATUS = VLIDORT_WARNING
            NGREEK_MOMENTS_INPUT = 2*NSTREAMS
            CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          ENDIF
        ELSE
          IF ( DO_SSCORR_NADIR .OR. DO_SSCORR_OUTGOING ) THEN
            IF ( NGREEK_MOMENTS_INPUT.LT.2*NSTREAMS-1 ) THEN
              MAIL =
     &   'Bad input: Fewer than 2N-1 expansion moments without delta-M'
              ACTION = 
     &       'Warning: Re-set NGREEK_MOMENTS_INPUT to 2N-1 internally'
              STATUS = VLIDORT_WARNING
              NGREEK_MOMENTS_INPUT = 2*NSTREAMS - 1
              CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
            ENDIF
          ENDIF
        ENDIF

      ELSE

C  Checks for Rayleigh only option
C   All warnings.

        IF ( DO_RAYLEIGH_ONLY ) THEN

          IF ( NGREEK_MOMENTS_INPUT.NE.2 ) THEN
            MAIL = 'Bad input: Rayleigh-only, expansion momemts NOT = 2'
            ACTION = 'Warning: Set NGREEK_MOMENTS_INPUT = 2 internally'
            STATUS = VLIDORT_WARNING
            NGREEK_MOMENTS_INPUT = 2
            CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          ENDIF

          IF ( DO_BVP_TELESCOPING ) THEN
            MAIL   = 
     &   'Bad input: Bvp telescoping not possible, Rayleigh only'
            ACTION = 'Warning: Turn off BVP_TELESCOPING internally'
            STATUS = VLIDORT_WARNING
            CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          ENDIF

          IF ( DO_SOLUTION_SAVING ) THEN
            MAIL   =
     &   'Bad input: Solution saving not possible, Rayleigh only'
            ACTION = 'Warning: Turn off SOLUTION_SAVING internally'
            STATUS = VLIDORT_WARNING
            CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          ENDIF

        ENDIF

C  Checks for Isotropic only option. All removed, 17 January 2006.

      ENDIF

C  reset solution saving and bvp telescoping flags
C  Do not need the isotropic-only options

c      DO_SOLUTION_SAVING =  ( DO_SOLUTION_SAVING .AND.
c     &     ((.NOT.DO_RAYLEIGH_ONLY).OR.(.NOT.DO_ISOTROPIC_ONLY)) )
c      DO_BVP_TELESCOPING =  ( DO_BVP_TELESCOPING .AND.
c     &     ((.NOT.DO_RAYLEIGH_ONLY).OR.(.NOT.DO_ISOTROPIC_ONLY)) )

      DO_SOLUTION_SAVING =  ( DO_SOLUTION_SAVING .AND.
     &     (.NOT.DO_RAYLEIGH_ONLY) )
      DO_BVP_TELESCOPING =  ( DO_BVP_TELESCOPING .AND.
     &               (.NOT.DO_RAYLEIGH_ONLY) )

C  BVP telescoping doesn't work with non-Lambertian surfaces
C   Reason, not yet coded. Theory already worked out.

      IF (  DO_BVP_TELESCOPING ) THEN
        IF ( .NOT. DO_LAMBERTIAN_SURFACE ) THEN
          MAIL   = 'BVP telescoping disabled for non-Lambertian'
          ACTION = 'Turn off DO_BVP_TELESCOPING flag, done internally'
          STATUS = VLIDORT_WARNING
          DO_BVP_TELESCOPING = .FALSE.
          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
      ENDIF

C  BVP telescoping doesn't work with non-Lambertian surfaces
C   Reason, not yet coded. Theory already worked out.
C    ---WARNING. BVP telescoping Flag turned off

      IF (  DO_BVP_TELESCOPING ) THEN
        IF ( .NOT. DO_LAMBERTIAN_SURFACE ) THEN
          MAIL   = 'BVP telescoping must be disabled, non-Lambertian'
          ACTION = 'Warning: DO_BVP_TELESCOPING turned off internally'
          STATUS = VLIDORT_WARNING
          DO_BVP_TELESCOPING = .FALSE.
          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
      ENDIF
   
C  Check azimuth-only conditions
C  =============================

C  Check no-Azimuth flag. Now set internally
C    ---WARNING. Do-no-Azimuth Flag turned on
c      IF ( .NOT.DO_NO_AZIMUTH ) THEN
c        IF ( DO_USER_VZANGLES. AND. N_USER_VZANGLES.EQ.1 ) THEN
c          IF ( USER_ANGLES_INPUT(1) .EQ. ZERO ) THEN
c            MAIL   ='Bad input: zenith-sky output requires no azimuth'
c            ACTION ='Warning: DO_NO_AZIMUTH flag set true internally'
c            STATUS = VLIDORT_WARNING
c            CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
c          ENDIF
c        ENDIF
c      ENDIF

C  Checks for Isotropic only option
C    ---WARNING. Do-no-Azimuth Flag turned on
C   Code disabled, 17 January 2006. Isotropic-only option no longer.
c      IF ( DO_ISOTROPIC_ONLY ) THEN
c        IF ( .NOT.DO_NO_AZIMUTH ) THEN
c          MAIL   =
c     &       'Bad input: no azimuth dependence for isotropic_only'
c          ACTION = 'Warning: DO_NO_AZIMUTH turned on internally'
c          STATUS = VLIDORT_WARNING
c          DO_NO_AZIMUTH = .TRUE.
c          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
c        ENDIF
c      ENDIF

C  Check single scattering correction and Do no Azimuth
C    ---WARNING. Do-no-Azimuth Flag turned off
C   Now been removed to the bookkeeping section
c      IF ( DO_SSCORR_NADIR ) THEN
c        IF ( DO_NO_AZIMUTH ) THEN
c          MAIL   =
c     &       'Bad input: need azimuth dependence for SS correction'
c          ACTION = 'Warning: DO_NO_AZIMUTH turned off internally'
c          STATUS = VLIDORT_WARNING
c          DO_NO_AZIMUTH = .FALSE.
c          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
c        ENDIF
c      ENDIF

C  Check: OLD single scattering correction and Do Rayleigh
C    ---WARNING. SS Flag turned off
C  Check only required for the diffuse field calculations (Version 2.3)

      IF ( DO_SSCORR_NADIR ) THEN
        IF ( DO_RAYLEIGH_ONLY .AND. .NOT. DO_SSFULL ) THEN
          MAIL   = 'Bad input: No SS correction for Rayleigh only'
          ACTION = 'Warning: DO_SSCORR_NADIR turned off internally'
          STATUS = VLIDORT_WARNING
          DO_SSCORR_NADIR = .FALSE.
          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
      ENDIF

C  Full-up single scatter, enabled 25 September 2007.
C   Single scatter corrections must be turned on

      IF ( DO_SSFULL ) THEN
        IF ( .not.DO_SSCORR_NADIR.and..not.DO_SSCORR_OUTGOING ) THEN
          MAIL   = 'Bad input: Full SS, must have one SSCORR flag set'
          ACTION = 'Full SS: default to use outgoing SS correction'
          STATUS = VLIDORT_WARNING
          DO_SSCORR_NADIR    = .FALSE.
          DO_SSCORR_OUTGOING = .TRUE.
          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
      ENDIF

C  Full-up single scatter, enabled 25 September 2007.
C   Diffuse-field Delta-M scaling must be turned off

      IF ( DO_SSFULL ) THEN
        IF ( DO_DELTAM_SCALING ) THEN        
          MAIL   = 'Bad input: Full SS, diffuse-field delta-M on'
          ACTION = 'Full SS: default to deltam_scaling = false'
          STATUS = VLIDORT_WARNING
          DO_DELTAM_SCALING   = .FALSE.
          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
      ENDIF

C  Check DB correction only for non-Lambertian surface
C    ---WARNING. Flag turned off
C    SS FULL calculation includes automatic Lambertian DB correction
C   This has been disabled. 26 November 2008

c      IF ( DO_DBCORRECTION ) THEN
c        IF ( DO_LAMBERTIAN_SURFACE ) THEN
c          MAIL   = 'Bad input: No DB correction for Lambertian only'
c          ACTION = 'Warning: DO_DBCORRECTION turned off internally'
c          STATUS = VLIDORT_WARNING
c          DO_DBCORRECTION = .FALSE.
c          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
c        ENDIF
c      ENDIF

C  Check surface inputs
C  ====================

      IF ( .NOT. DO_LAMBERTIAN_SURFACE ) THEN
       DO K = 1, N_BRDF_KERNELS
        IF ( WHICH_BRDF(K).GT.MAXBRDF_IDX.OR.WHICH_BRDF(K).LE.0) THEN
          MAIL   = 'Bad input: BRDF Index not on list of indices'
          ACTION = 'look in VLIDORT.PARS for correct index'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ELSE
          IF ( BRDF_NAMES(K).NE.BRDF_CHECK_NAMES(WHICH_BRDF(K)) ) THEN
            MAIL   = 'Bad input: BRDF kernel name not corresponding'
            ACTION = 'look in VLIDORT.PARS for correct name'
            STATUS = VLIDORT_SERIOUS
            CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          ENDIF
        ENDIF
       ENDDO
      ENDIF

C  make sure the Lambertian surface is defined
C  set defaults for Lambertian case (Derived input)

      IF ( DO_LAMBERTIAN_SURFACE ) THEN
        N_BRDF_KERNELS = 1
        LAMBERTIAN_KERNEL_FLAG(1) = .TRUE.
        BRDF_FACTORS(1) = LAMBERTIAN_ALBEDO
      ENDIF

C  Check thermal inputs
C  ====================

C  If thermal transmittance only, check thermal flag

      IF ( .NOT. DO_THERMAL_EMISSION ) THEN
       IF ( DO_THERMAL_TRANSONLY ) THEN
         MAIL = 'Bad input: No thermal, must turn off transonly flag'
         ACTION = 'Warning: DO_THERMAL_TRANSONLY turned off internally'
         STATUS = VLIDORT_WARNING
         DO_THERMAL_TRANSONLY = .FALSE.
         CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
       ENDIF
      ENDIF

C  Switch off a bunch of flags

      IF ( DO_THERMAL_EMISSION ) THEN
       IF ( DO_THERMAL_TRANSONLY ) THEN
         DO_RAYLEIGH_ONLY  = .FALSE.
         DO_DELTAM_SCALING = .FALSE.
       ENDIF
      ENDIF

C  No solar sources for thermal transmittance

      IF ( DO_THERMAL_TRANSONLY ) THEN
        IF ( DO_SOLAR_SOURCES ) THEN
         MAIL = 'Bad input: thermal tranmsittance, must turn off solar'
         ACTION = 'Warning: DO_SOLAR_SOURCES turned off internally'
         STATUS = VLIDORT_WARNING
         DO_SOLAR_SOURCES = .FALSE.
         CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
       ENDIF
      ENDIF

C  check viewing geometry input
C  ============================

C  Check earth radius (Chapman function only)
C    ---WARNING. Default value of 6371.0 will be set

      IF ( DO_CHAPMAN_FUNCTION ) THEN
        IF ( .NOT. DO_PLANE_PARALLEL ) THEN
          IF ( EARTH_RADIUS.LT.6320.0D0 .OR.
     &         EARTH_RADIUS.GT.6420.0D0 ) THEN
            MAIL   = 'Bad input: Earth radius outside of [6320-6420]'
            ACTION = 'Warning: default value of 6371.0 was set'
            STATUS = VLIDORT_WARNING
            EARTH_RADIUS = 6371.0D0
            CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          ENDIF
        ENDIF
      ENDIF

C  Check dimensioning on Legendre numbers (refractive geometry only)

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN
        NALLSTREAMS = N_SZANGLES*NLAYERS + NSTREAMS + N_USER_VZANGLES
        IF ( NALLSTREAMS .GT. MAX_ALLSTRMS_P1 ) THEN
          MAIL   = 'Dimensioning error for refractive beam angles'
          ACTION = 'Increase dimension MAX_ALLSTRMS_P1'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
      ENDIF

C  Check GEOMETRY_SPECHEIGHT (only for outgoing sphericity correction)
C    GEOMETRY_SPECHEIGHT cannot be greater than HEIGHT_GRID(NLAYERS)

      IF ( DO_SSCORR_OUTGOING ) THEN
        IF ( GEOMETRY_SPECHEIGHT .GT. HEIGHT_GRID(NLAYERS) ) THEN
          MAIL   = 'GEOMETRY_SPECHEIGHT must be =< Input BOA-HEIGHT '
          ACTION = 'Warning: Internal Re-set of GEOMETRY_SPECHEIGHT '
          STATUS = VLIDORT_WARNING
          GEOMETRY_SPECHEIGHT  = HEIGHT_GRID(NLAYERS)
          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
      ENDIF
      
C  Check solar zenith angle input

      DO I = 1, N_SZANGLES
        IF ( SZANGLES(I) .LT. ZERO .OR.
     &       SZANGLES(I).GE.90.0D0 ) THEN
          WRITE(C2,'(I2)')I
          MAIL   = 'Bad input: out-of-range solar angle, no. '//C2
          ACTION = 'Look at SZANGLES input, should be < 90 & > 0'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
      ENDDO

C  Check zenith tolerance input. NO longer required, 01/17/06.
C    ---Formerly: WARNING. Default of 0.001 will be set

C  Check relative azimuths

      LOOP = .TRUE.
      I = 0
      DO WHILE (LOOP .AND. I.LT.N_USER_RELAZMS)
        I = I + 1
        IF ( USER_RELAZMS(I) .GT. 360.0D0   .OR.
     &       USER_RELAZMS(I) .LT. ZERO ) THEN
          WRITE(C2,'(I2)')I
          MAIL = 'Bad input: out-of-range azimuth angle, no. '//C2
          ACTION = 'Look at azimuth angle input, should be in [0,360]'
          LOOP = .FALSE.
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
      ENDDO

C  limits on user-defined options

c      IF ( .NOT. DO_USER_VZANGLES .AND..NOT.DO_QUAD_OUTPUT ) THEN
c        MAIL   = 'Bad input: No angular stream output is specified'
c        ACTION = 'Check DO_USER_VZANGLES and DO_QUAD_OUTPUT flags'
c        STATUS = VLIDORT_SERIOUS
c        CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
c      ENDIF

C  check user-defined stream angles (should always be [0,90])

      IF ( DO_USER_VZANGLES ) THEN
        LOOP = .TRUE.
        I = 0
        DO WHILE (LOOP .AND. I.LT.N_USER_VZANGLES)
          I = I + 1
          IF ( USER_VZANGLES(I) .GT. 90.0   .OR.
     &         USER_VZANGLES(I) .LT. ZERO ) THEN
            WRITE(C2,'(I2)')I
            MAIL   = 'Bad input: out-of-range user stream, no. '//C2
            ACTION = 'Look at user viewing zenith angle input'
            LOOP = .FALSE.
            STATUS = VLIDORT_SERIOUS
            CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          ENDIF
        ENDDO
      ENDIF

C  Re-order the input angles
c  -------------------------

C  This section has been truncated. 28 March 2007

      IF ( DO_USER_VZANGLES ) THEN
        N_OUT_STREAMS  = N_USER_VZANGLES
        IF ( N_OUT_STREAMS .EQ. 1 ) THEN
          OUT_ANGLES(1) =  USER_VZANGLES(1)
        ELSE
          DO I = 1, N_USER_VZANGLES
            ALL_ANGLES(I) = USER_VZANGLES(I)
          ENDDO
          CALL INDEXX
     &          ( N_OUT_STREAMS, USER_VZANGLES, INDEX_ANGLES )
          DO I = 1, N_OUT_STREAMS
            OUT_ANGLES(I) = ALL_ANGLES(INDEX_ANGLES(I))
            USER_VZANGLES(I) = OUT_ANGLES(I)
          ENDDO
        ENDIF
      ENDIF

C  Check height grid input (Chapman function only)

      IF ( DO_CHAPMAN_FUNCTION ) THEN
        LOOP = .TRUE.
        I = 0
        DO WHILE (LOOP .AND. I.LT.NLAYERS)
          I = I + 1
          IF ( HEIGHT_GRID(I-1).LE.HEIGHT_GRID(I) ) THEN
            WRITE(C2,'(I2)')I
            MAIL   =
     &'Bad input: Height-grid not monotonically decreasing; Layer '//C2
            ACTION = 'Look at Height-grid input'
            LOOP = .FALSE.
            STATUS = VLIDORT_SERIOUS
            CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          ENDIF
        ENDDO
      ENDIF

C  Check vertical outputs
C  ----------------------

C  check vertical output levels (should always be within atmosphere!)

      LOOP = .TRUE.
      I = 0
      DO WHILE (LOOP .AND. I.LT.N_USER_LEVELS)
        I = I + 1
        IF ( USER_LEVELS(I) .GT. DBLE(NLAYERS) .OR.
     &       USER_LEVELS(I) .LT. ZERO )  THEN
          WRITE(C2,'(I2)')I
          MAIL   = 'Bad input: Out of range for level choice # '//C2
          ACTION = 'Re-set level output '
          LOOP = .FALSE.
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
      ENDDO

C  check repetition of vertical output choices

      UTA = 0
      LOOP = .TRUE.
      DO WHILE ( LOOP .AND. UTA .LT. N_USER_LEVELS )
        UTA = UTA + 1
        XT = USER_LEVELS(UTA)
        NSTART = 0
        DO N = 1, N_USER_LEVELS
          IF ( XT .EQ. USER_LEVELS(N)) NSTART = NSTART + 1
        ENDDO
        IF ( NSTART .NE. 1 ) THEN
          LOOP = .FALSE.
          MAIL   = 'Bad input: repetition of vertical output choice'
          ACTION = 'Re-set level output '
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
      ENDDO


C  Check optical depth input  (OLD CODE)
C  =====================================

c      TAUATM     = ZERO
c      LOCALTAU(0) = TAUATM
c      DO N = 1, NLAYERS
c        TAUATM = TAUATM + DELTAU_VERT_INPUT(N)
c        LOCALTAU(N) = TAUATM
c      ENDDO
C  Check non-negative optical thickness values
c      DO L = 1, NLAYERS
c        IF ( DELTAU_VERT_INPUT(L).LE.ZERO ) THEN
c          WRITE(C3,'(I3)')L
c          MAIL = 'Bad input: optical thickness <= 0, layer '//C3
c          ACTION = 'Check optical thickness input'
c          STATUS = VLIDORT_SERIOUS
c          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
c        ENDIF
c      ENDDO
C  check user-defined optical depths (should always be within atmosphere!)
c      IF ( DO_USER_TAUS ) THEN
c        LOOP = .TRUE.
c        I = 0
c        DO WHILE (LOOP .AND. I.LT.N_OUT_USERTAUS)
c          I = I + 1
c          IF ( USER_TAUS_INPUT(I) .GT. TAUATM )  THEN
c            WRITE(C2,'(I2)')I
c            MAIL   = 'Bad input: user optical depth too big, '//C2
c            ACTION = 'Check optical depth values in range'
c            LOOP = .FALSE.
c            STATUS = VLIDORT_SERIOUS
c            CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
c          ENDIF
c        ENDDO
c      ELSE IF ( DO_LBOUND_TAUS ) THEN
c        IF ( N_OUT_USERTAUS .GT. NLAYERS + 1 ) THEN
c          MAIL   = 'Bad input: Too many Layer boundary optical depths'
c          ACTION = 'Check N_OUT_USERTAUS = NLAYERS + 1'
c          LOOP = .FALSE.
c          STATUS = VLIDORT_SERIOUS
c          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
c        ENDIF
c        DO I = 1, N_OUT_USERTAUS
c          IF ( LBOUND_TAUS_INPUT(I) .GT. NLAYERS ) THEN
c            WRITE(C2,'(I2)')LBOUND_TAUS_INPUT(I)
c            MAIL   = 'Bad input: Layer boundary index > NLAYERS: '//C2
c            ACTION = 'Check LBOUND_TAUS_INPUT layer boundary indices'
c            LOOP = .FALSE.
c            STATUS = VLIDORT_SERIOUS
c            CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
c          ENDIF
c        ENDDO
c      ENDIF
C  check Number of offgrid user-defined optical depths
c      IF ( DO_USER_TAUS ) THEN
c        NSTART = 0
c        DO UTA = 1, N_OUT_USERTAUS
c          TAU = USER_TAUS_INPUT(UTA)
c          LOOP = .TRUE.
c          N = 0
c          DO WHILE (LOOP.AND.N.LT.NLAYERS)
c            N = N + 1
c            IF ( DABS(TAU-LOCALTAU(N-1)).LT.SMALLNUM ) THEN
c              NSTART = NSTART + 1
c            ENDIF
c          ENDDO
c        ENDDO
c        N_OFFGRID_LOCAL = N_OUT_USERTAUS - NSTART
c        IF ( N_OFFGRID_LOCAL .GT. MAX_OFFGRID_USERTAUS ) THEN
c          WRITE(C2,'(I2)')N_OFFGRID_LOCAL
c          MAIL   = 'Bad input: too many offgrid optical depths : '//C2
c          ACTION =
c     &        'Check number of offgrid Tau does not exceed dimensioning'
c          LOOP = .FALSE.
c          STATUS = VLIDORT_SERIOUS
c          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
c        ENDIF
c      ENDIF
C  check repetition of user-defined optical depth input
c      IF ( DO_USER_TAUS ) THEN
c        UTA = 0
c        LOOP = .TRUE.
c        DO WHILE ( LOOP .AND. UTA .LT. N_OUT_USERTAUS )
c          UTA = UTA + 1
c          TAU = USER_TAUS_INPUT(UTA)
c          NSTART = 0
c          DO N = 1, N_OUT_USERTAUS
c            IF ( TAU .EQ. USER_TAUS_INPUT(N)) NSTART = NSTART + 1
c          ENDDO
c          IF ( NSTART .NE. 1 ) THEN
c            LOOP = .FALSE.
c            MAIL   = 'Bad input: repetition of optical depth input'
c            ACTION = 'Check user-defined optical depth input'
c            STATUS = VLIDORT_SERIOUS
c            CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
c          ENDIF
c        ENDDO
c      ENDIF
C  check repetition of layer boundary optical depth input
c      IF ( DO_LBOUND_TAUS ) THEN
c        UTA = 0
c        LOOP = .TRUE.
c        DO WHILE (LOOP .AND. UTA .LT. N_OUT_USERTAUS )
c          UTA = UTA + 1
c          N1 = LBOUND_TAUS_INPUT(UTA)
c          NSTART = 0
c          DO N = 1, N_OUT_USERTAUS
c            IF ( LBOUND_TAUS_INPUT(N) .EQ. N1 ) NSTART = NSTART + 1
c          ENDDO
c          IF ( NSTART .GT. 1 ) THEN
c            LOOP = .FALSE.
c            MAIL   = 'Bad input: repetition of optical depth input'
c            ACTION = 'Check optical depth layer boundary input'
c            STATUS = VLIDORT_SERIOUS
c            CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
c          ENDIF
c        ENDDO
c      ENDIF

C  check geophysical scattering inputs
C  -----------------------------------

C  check single scatter albedos

      IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
       DO L = 1, NLAYERS
        IF ( OMEGA_TOTAL_INPUT(L).GT.ONE-OMEGA_SMALLNUM ) THEN
          WRITE(C3,'(I3)')L
          MAIL = 'Bad input: SS-albedo too close to 1, layer '//C3
          ACTION = 'Check SS-albedo input'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ELSE IF ( OMEGA_TOTAL_INPUT(L).LT.OMEGA_SMALLNUM ) THEN
          WRITE(C3,'(I3)')L
          MAIL = 'Bad input: SS-albedo too close to 0, layer '//C3
          ACTION = 'Check SS-albedo input'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
       ENDDO
      ENDIF

C  solar beam, cannot be too small

      IF ( DO_SOLAR_SOURCES ) THEN
        DO L = 1, NLAYERS
         IF ( OMEGA_TOTAL_INPUT(L).LT.OMEGA_SMALLNUM ) THEN
          WRITE(C3,'(I3)')L
          MAIL = 'Bad input: SS-albedo too close to 0, layer '//C3
          ACTION = 'Check SS-albedo input'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
         ENDIF
        ENDDO
      ENDIF
      
C  Check first phase function moments

      IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
       DO L = 1, NLAYERS
        IF ( GREEKMAT_TOTAL_INPUT(0,L,1).NE.ONE ) THEN
          WRITE(C3,'(I3)')L
          MAIL ='First phase moment (GREEK_11) not 1 for layer '//C3
          ACTION = 'Check First phase function moment'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
       ENDDO
      ENDIF

C  Additional check on non-negativity of (1,1) moments.
C    Furnished 27 September 2007, as a result of input-checking.
c      DO N = 1, NLAYERS
c       DO L = 1, NGREEK_MOMENTS_INPUT
c        IF ( GREEKMAT_TOTAL_INPUT(L,N,1).LT.ZERO ) THEN
c         WRITE(C3,'(I3)')N
c         MAIL ='Some moments (GREEK_11) are NEGATIVE for layer '//C3
c         ACTION = 'Check Greek moments input - some bad values!'
c         STATUS = VLIDORT_SERIOUS
c         CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
c        ENDIF
c       ENDDO
c      ENDDO

C  Specialist options, Check the given input

      IF ( DO_SPECIALIST_OPTION_2 ) THEN
        IF ( NLAYERS_NOMS .GT. NLAYERS-1 ) THEN
          MAIL   = 'Bad input specialist option 2: NLAYERS_NOMS'
          ACTION = 'Check NLAYERS_NOMS must be less than NLAYERS'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
      ENDIF

      IF ( DO_SPECIALIST_OPTION_3 ) THEN
        IF ( NLAYERS_CUTOFF .lt. 2 ) THEN
          MAIL   = 'Bad input specialist option 3: NLAYERS_CUTOFF'
          ACTION = 'Check NLAYERS_CUTOFF must be greater than 1'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE VLIDORT_DERIVE_INPUT ( STATUS_SUB, MESSAGE )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include file with input variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'

C  include file with bookkeeping variables to be assigned

      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include file with setups (for taugrid_input)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  Status output

      INTEGER          STATUS_SUB
      CHARACTER*(*)    MESSAGE

C  local variables
C  ---------------

      DOUBLE PRECISION MU1, MU2, DT, RT
      INTEGER          I, UT, N, UTA, NSTART, M, NAP, NS
      INTEGER          O1, O2, O1_S, MS2, P, NA, Q, QC, L, EC
      LOGICAL          LOOP, LOCAL_NADIR_ONLY
      INTEGER          CONV_START_STREAMS

C  set additional numbers (derived input)
C  ======================

C  set status

      STATUS_SUB = VLIDORT_SUCCESS
      MESSAGE    = ' '

C  Automatic input
C    Flux factor set to unity. (21 December 2005)

      DO_ALL_FOURIER        = .FALSE.
      DO_CLASSICAL_SOLUTION = .TRUE.
      DO_DIRECT_BEAM        = .TRUE.
!      FLUX_FACTOR           = ONE

C  SSFULL flag cancels some other flags. Not the DB correction !!!!

      IF ( DO_SSFULL ) THEN
        DO_DOUBLE_CONVTEST = .FALSE.
        DO_SOLUTION_SAVING = .FALSE.
        DO_BVP_TELESCOPING = .FALSE.
        DO_ADDITIONAL_MVOUT = .FALSE.
        DO_MVOUT_ONLY       = .FALSE.
      ENDIF

C  Mode of operation
C   SS outgoing sphericity option, added 31 January 2007
C   SS full calculation option added 25 September 2007.

      DO_MSMODE_VLIDORT = .FALSE.
      IF ( DO_FULLRAD_MODE ) THEN
        IF ( .NOT. DO_SSFULL ) THEN 
          IF ( DO_SSCORR_NADIR .OR. DO_SSCORR_OUTGOING ) THEN
            DO_MSMODE_VLIDORT = .TRUE.
          ENDIF
        ENDIF
      ELSE
        IF ( .NOT. DO_SSFULL ) THEN
          DO_MSMODE_VLIDORT = .TRUE.
        ENDIF
      ENDIF

C  General flag

      DO_USER_STREAMS = DO_USER_VZANGLES

C  New section. Setting the DO_NO_AZIMUTH flag
C    Rt Solutions. 17 January 2006. R. Spurr and V. Natraj.

C  DO_NO_AZIMUTH should be set internally only when:
C     (a) NSTOKES = 1 and nadir view only, and no SS correction
C     (b) MIFLUX_ONLY flag is true.
C     (c) SSFULL calculation flag is set

      LOCAL_NADIR_ONLY = .FALSE.
      IF ( N_USER_VZANGLES.EQ.1.and.USER_VZANGLES(1).EQ.ZERO ) THEN
        LOCAL_NADIR_ONLY = .TRUE.
      ENDIF
      DO_NO_AZIMUTH = .FALSE.
      IF ( (NSTOKES.EQ.1.AND.LOCAL_NADIR_ONLY.AND.(.NOT.DO_SSFULL.OR.
     &      .NOT.DO_SSCORR_NADIR.OR..not.DO_SSCORR_OUTGOING))
     &       .OR.DO_MVOUT_ONLY  ) THEN
        DO_NO_AZIMUTH = .TRUE.
      ENDIF

C  Directional indices

      IF ( DO_UPWELLING .AND. DO_DNWELLING ) THEN
        N_DIRECTIONS = 2
        WHICH_DIRECTIONS(1) = UPIDX
        WHICH_DIRECTIONS(2) = DNIDX
      ELSE
        N_DIRECTIONS = 1
        WHICH_DIRECTIONS(2) = 0
        IF ( DO_UPWELLING ) THEN
          WHICH_DIRECTIONS(1) = UPIDX
        ELSE IF ( DO_DNWELLING) THEN
          WHICH_DIRECTIONS(1) = DNIDX
        ENDIF
      ENDIF

C  Flux vector set unity input.

      FLUXVEC(1) = ONE
      FLUXVEC(2) = ZERO
      FLUXVEC(3) = ZERO
      FLUXVEC(4) = ZERO

C  Convert Surface emission input.............. NOT USED....
C  Input values should be Watts/ sq m
C      IF ( DO_SURFACE_EMISSION ) THEN
C        FP_SURFBB = PI4 * SURFBB
C      ENDIF

C  make sure the Lambertian surface is defined
C  set defaults for Lambertian case (Derived input)

      IF ( DO_LAMBERTIAN_SURFACE ) THEN
        N_BRDF_KERNELS = 1
        LAMBERTIAN_KERNEL_FLAG(1) = .TRUE.
        BRDF_FACTORS(1) = LAMBERTIAN_ALBEDO
      ENDIF

C  make sure the Lambertian surface is defined if part of a Kernel BRDF
C   Code Introduced 28 August 2009

      IF ( .NOT.DO_LAMBERTIAN_SURFACE ) THEN
        DO I = 1, N_BRDF_KERNELS
          IF ( BRDF_NAMES(I) .EQ. 'Lambertian' ) THEN
             LAMBERTIAN_KERNEL_FLAG(I) = .TRUE.
          ELSE
             LAMBERTIAN_KERNEL_FLAG(I) = .FALSE.
          ENDIF
        ENDDO
      ENDIF

C  Bookkeeping for surface kernel Cox-Munk types
C   Only the Giss CoxMunk kernel is vectorized (as of 19 January 2006)

      N_BRDF_STOKESSQ  = 1
      IF ( .NOT. DO_LAMBERTIAN_SURFACE ) THEN
        DO I = 1, N_BRDF_KERNELS
          IF ( BRDF_NAMES(I) .EQ. 'Cox-Munk  ' .OR.
     &         BRDF_NAMES(I) .EQ. 'GissCoxMnk' ) THEN
            N_BRDF_PARAMETERS(I) = 3
            IF ( DO_SHADOW_EFFECT ) THEN
              BRDF_PARAMETERS(I,3) = ONE
            ELSE
              BRDF_PARAMETERS(I,3) = ZERO
            ENDIF
          ENDIF
          IF ( BRDF_NAMES(I) .EQ. 'GissCoxMnk' ) THEN
            N_BRDF_STOKESSQ = NSTOKES * NSTOKES
          ENDIF
       ENDDO
      ENDIF

C  set (tempo)  the number of BRDF entries. Not used Lambertian
!      N_EMISS_STOKESSQ = 1

C  Number of moments. Isotropic option removed 17 January 2006.

      IF ( DO_RAYLEIGH_ONLY ) THEN
        NMOMENTS = 2
      ENDIF
      IF ( .NOT.DO_RAYLEIGH_ONLY ) THEN
        NMOMENTS = MIN ( 2 * NSTREAMS - 1, NGREEK_MOMENTS_INPUT )
      ENDIF

C  total quadratures (up and down)

      NSTREAMS_2 = 2*NSTREAMS

C  Additional quantities (Stokes)

      NSTOKES_SQ     = NSTOKES * NSTOKES
      NSTKS_NSTRMS   = NSTOKES * NSTREAMS
      NSTKS_NSTRMS_2 = NSTOKES * NSTREAMS_2

C  Mueller index

      DO O1 = 1, MAXSTOKES
        O1_S = MAXSTOKES*(O1 - 1)
        DO O2 = 1, MAXSTOKES
          MUELLER_INDEX(O1,O2) = O1_S + O2
        ENDDO
      ENDDO

C  Greek matrix index

      GREEKMAT_INDEX(1) = 1
      GREEKMAT_INDEX(2) = 6
      GREEKMAT_INDEX(3) = 2
      GREEKMAT_INDEX(4) = 11
      GREEKMAT_INDEX(5) = 12
      GREEKMAT_INDEX(6) = 16

C  Check number of particular solution modes
C   Current default is NPARTICSOLS = 1

      IF ( FLUXVEC(3).EQ.ZERO. AND.FLUXVEC(4).EQ.ZERO ) THEN
        NPARTICSOLS = 1
      ELSE
        NPARTICSOLS = 2
      ENDIF

C  D matrices, and multiply by Fluxvector

      DO O1 = 1, MAXSTOKES
        DFLUX(O1) = ZERO
        DO O2 = 1, MAXSTOKES
          DMAT(O1,O2) = ZERO
          I44(O1,O2)  = ZERO
          DO P = 1, NPARTICSOLS
            DMAT_PSOLS(O1,O2,P) = ZERO
          ENDDO
        ENDDO
        I44(O1,O1)  = ONE
      ENDDO
      MS2 = MAXSTOKES/2
      DO O1 = 1, MS2
        O2 = O1 + MS2
        DANTE(O1) = ONE
        DANTE(O2) = -ONE
        DMAT(O1,O1) = ONE
        DMAT(O2,O2) = -ONE
        DMAT_PSOLS(O1,O1,1) = ONE
        DMAT_PSOLS_FLUXVEC(O1,1) = FLUXVEC(O1)
        DFLUX(O1) = FLUXVEC(O1)
      ENDDO
      IF ( NPARTICSOLS .EQ. 2 ) THEN
        DO O1 = 1, MS2
          O2 = O1 + MS2
          DMAT_PSOLS(O2,O2,2) = ONE
          DMAT_PSOLS_FLUXVEC(O2,1) = FLUXVEC(O2)
          DFLUX(O2) = FLUXVEC(O2)
        ENDDO
      ENDIF

C  Set Quadrature abscissae and weights

      CALL GAULEG(ZERO,ONE,QUAD_STREAMS,QUAD_WEIGHTS,NSTREAMS)

C  Following code disabled from earlier versions
C      IF ( DO_FULL_QUADRATURE ) THEN
C        CALL GAULEG(-ONE,ONE,X2,A2,NSTR2)
C        DO I = 1, NSTREAMS
C          I1 = I + NSTREAMS
C          X(I) = X2(I1)
C          A(I) = A2(I1)
C        ENDDO
C      ENDIF

C  set auxiliary quantities

      DO I = 1, NSTREAMS
        QUAD_STRMWTS(I) = QUAD_STREAMS(I)*QUAD_WEIGHTS(I)
        QUAD_HALFWTS(I) = HALF * QUAD_WEIGHTS(I)
        QUAD_ANGLES(I)  = DACOS(QUAD_STREAMS(I))/DEG_TO_RAD
        QUAD_SINES(I)   = DSQRT(ONE-QUAD_STREAMS(I)*QUAD_STREAMS(I))
      ENDDO

C  size of boundary value problem matrices and vectors

      NTOTAL = NLAYERS*NSTKS_NSTRMS_2

C  number of sub and super diagonals in band matrix (boundary value problem)

      IF ( NLAYERS .EQ. 1 ) THEN
        N_SUBDIAG = 2*NSTKS_NSTRMS - 1
        N_SUPDIAG = 2*NSTKS_NSTRMS - 1
      ELSE
        N_SUBDIAG = 3*NSTKS_NSTRMS - 1
        N_SUPDIAG = 3*NSTKS_NSTRMS - 1
      ENDIF

C  solar zenith angle cosines/sines

      NBEAMS = N_SZANGLES
      DO I = 1, N_SZANGLES
        COS_SZANGLES(I)  = DCOS ( SZANGLES(I) * DEG_TO_RAD )
        SIN_SZANGLES(I) = DSQRT(ONE-COS_SZANGLES(I)*COS_SZANGLES(I))
      ENDDO

C  Set average cosines in the refractive geometry case

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN
        DO I = 1, N_SZANGLES
          MU1 = DCOS(SZA_LOCAL_INPUT(0,I)*DEG_TO_RAD)
          DO N = 1, NLAYERS
            MU2 = DCOS(SZA_LOCAL_INPUT(N,I)*DEG_TO_RAD)
            SUN_SZA_COSINES(N,I) = HALF * ( MU1 + MU2 )
            MU1 = MU2
         ENDDO
        ENDDO
      ENDIF

C  Set performance flags
C  ---------------------

C  New section for Version 2.0 and higher.

C  Specialist Option 2: Set solution saving. BVP Telescoping
C   Set the layer-scattering flags to be false for all N up to NLAYERS_NOMS

      IF ( DO_SPECIALIST_OPTION_2 ) THEN
        DO_SOLUTION_SAVING = .TRUE.
        DO_BVP_TELESCOPING = .TRUE.
        DO M = 0, NMOMENTS
          BVP_REGULAR_FLAG(M) = .FALSE.
          DO N = 1, NLAYERS_NOMS
            DO_LAYER_SCATTERING(M,N) = .FALSE.
          ENDDO
          DO N = NLAYERS_NOMS + 1, NLAYERS
            DO_LAYER_SCATTERING(M,N) = .TRUE.
          ENDDO
        ENDDO
        GO TO 4545
      ENDIF

C  Set the layer scattering flags
C  set the BVP telescoping flags

C  initialise (M = Fourier index) to normal mode.

      DO M = 0, NMOMENTS
        BVP_REGULAR_FLAG(M) = .TRUE.
        DO N = 1, NLAYERS
          DO_LAYER_SCATTERING(M,N) = .TRUE.
        ENDDO
      ENDDO

C  Special clause for transmittance-only calculation
C   New flag, 31 July 2007. Move on to other bookkeeping

      IF ( DO_THERMAL_TRANSONLY ) THEN
        DO_SOLUTION_SAVING = .TRUE.
        DO N = 1, NLAYERS
          OMEGA_TOTAL_INPUT(N) = ZERO
          DO_LAYER_SCATTERING(0,N) = .FALSE.
        ENDDO
        GO TO 5557
      ENDIF

C  for M > 2 terms, if solution saving flag is set...
C  .... examine  Greek matrix entries (1,1) - no scattering if
C      they are all zero for L > M - 1.
C  [Equivalent to examining phase function moments]

C  Addition of Solution_saving flag, 22 November 2009.
C    Bug spotted by V. Natraj, 20 November 2009

      IF ( DO_SOLAR_SOURCES ) THEN
       IF ( DO_SOLUTION_SAVING ) THEN
        DO M = 3, NMOMENTS
          QC = NMOMENTS - M + 1
          DO N = 1, NLAYERS
            Q = 0
            DO L = M, NMOMENTS
              IF(GREEKMAT_TOTAL_INPUT(L,N,1).EQ.ZERO)Q=Q+1
            ENDDO
            DO_LAYER_SCATTERING(M,N) = (Q.LT.QC)
          ENDDO
        ENDDO
       ENDIF
      ENDIF

C  BVP telescoping (only if do_solution_saving is set)

      IF ( DO_SOLAR_SOURCES ) THEN
       IF ( DO_SOLUTION_SAVING ) THEN
        IF ( DO_BVP_TELESCOPING ) THEN
          DO M = 3, NMOMENTS
            Q = 0
            DO N = 1, NLAYERS
              IF (.NOT.DO_LAYER_SCATTERING(M,N))Q = Q + 1
            ENDDO
            IF ( Q.GT.1) BVP_REGULAR_FLAG(M) = .FALSE.
          ENDDO
        ENDIF
       ENDIF
      ENDIF

C    Set of Telescoped layers must be contiguous

      IF ( DO_SOLAR_SOURCES ) THEN
       IF ( DO_SOLUTION_SAVING ) THEN

        IF ( DO_BVP_TELESCOPING ) THEN
         DO M = 3, NMOMENTS
          NS = 0
          NAP = 0
          DO N = 1, NLAYERS
           IF ( DO_LAYER_SCATTERING(M,N) ) THEN
            NS = NS + 1
            NA = N
            IF ( NS.GT.1)  THEN
              IF ( NA.NE.NAP+1 ) THEN
                STATUS_SUB = VLIDORT_WARNING
                GO TO 4564
              ENDIF
            ENDIF
            NAP = NA
           ENDIF
          ENDDO
         ENDDO
        ENDIF

C  Collect warning and re-set default option

 4564   continue
        if ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
         MESSAGE = 'Telescoped layers not contiguous: turn off option'
         DO_BVP_TELESCOPING = .FALSE.
         DO M = 3, NMOMENTS
           BVP_REGULAR_FLAG(M) = .TRUE.
         ENDDO
        ENDIF

C  End clause

       ENDIF
      ENDIF
      
C  Continuation point for avoiding the default, if using Specialist # 2

 4545 CONTINUE

C  Save calculation time for single scattering.
C   Remove isotropic-only option. 17 January 2006.
C   Bug. 22 November 2006. L must be < Ngreek_moments_input in DO WHILE
C   Bug 31 January 2007. LAYER_MAXMOMENTS should be = nmoments_input

      IF ( DO_SOLAR_SOURCES ) THEN
       IF ( DO_SSCORR_NADIR .OR. DO_SSCORR_OUTGOING ) THEN
        IF ( DO_RAYLEIGH_ONLY ) THEN
          DO N = 1, NLAYERS
            LAYER_MAXMOMENTS(N) = 2
          ENDDO
        ELSE
         DO N = 1, NLAYERS
            L = 2
            LOOP = .TRUE.
            DO WHILE (LOOP.AND.L.LT.NGREEK_MOMENTS_INPUT)
              L = L + 1
              LOOP= (GREEKMAT_TOTAL_INPUT(L,N,1).NE.ZERO)
            ENDDO
c            LAYER_MAXMOMENTS(N) = L - 1
            LAYER_MAXMOMENTS(N) = L 
          ENDDO
        ENDIF
       ENDIF
      ENDIF

C  Set the Eigensolver Mask. New section 21 December 2005.
C   1. Always set the Real Eigensolver for Rayleigh-only or Scalar-only
C   2. Otherwise, examine the occurrence of the "epsilon" Greek constant

      IF ( DO_RAYLEIGH_ONLY .OR. NSTOKES.EQ.1 ) THEN
       DO M = 0, NMOMENTS
        DO N = 1, NLAYERS
         DO_REAL_EIGENSOLVER(M,N) = .TRUE.
        ENDDO
       ENDDO
      ELSE
       DO M = 0, NMOMENTS
        DO N = 1, NLAYERS
         IF ( DO_LAYER_SCATTERING(M,N) ) THEN
           EC = 0
           DO L = M, NMOMENTS
            IF (GREEKMAT_TOTAL_INPUT(L,N,12).NE.ZERO) EC = EC + 1
           ENDDO
           DO_REAL_EIGENSOLVER(M,N) = (EC.EQ.ZERO)
         ENDIF
        ENDDO
       ENDDO
      ENDIF

C  debug

c      DO M = 0, NMOMENTS
c       write(*,'(I4,L2,5x,7L2)')M,BVP_REGULAR_FLAG(M),
c     &             (DO_LAYER_SCATTERING(M,N),N=1,7)
c        ENDDO
c      PAUSE

C  COntinuation point

 5557 continue

C  User stream cosines and secants

      IF ( DO_USER_STREAMS ) THEN
        N_USER_STREAMS = N_USER_VZANGLES
        DO I = 1, N_USER_STREAMS
c          USER_STREAMS(I) = DCOS(DEG_TO_RAD*USER_VZANGLES(I))
          USER_STREAMS(I) = DCOS(DEG_TO_RAD*USER_VZANGLES_ADJUST(I))
          USER_SECANTS(I) = ONE / USER_STREAMS(I)
          USER_SINES(I) = DSQRT(ONE-USER_STREAMS(I)*USER_STREAMS(I))
        ENDDO
      ENDIF

C  number of tests to be applied for convergence
C  number of streams for convergence
C    - Zenith tolerance no longer applies for the vector code
C    - Set CONV_START_STREAMS = 1. Do not examine for "Zenith_tolerance"

      CONV_START_STREAMS = 1
      N_CONV_STREAMS = N_OUT_STREAMS - CONV_START_STREAMS + 1
      N_CONVTESTS = N_USER_RELAZMS * N_CONV_STREAMS * N_DIRECTIONS
      N_CONVTESTS = N_CONVTESTS * N_USER_LEVELS

C  Sort out User vertical level outputs
C  ------------------------------------

C  Sort in ascending order

      IF ( N_USER_LEVELS .GT. 1 ) THEN
        CALL HPSORT(N_USER_LEVELS,USER_LEVELS)
      ENDIF

C  mark all output levels not equal to layer boundary values

      NSTART = 0
      UT = 0
      DO UTA = 1, N_USER_LEVELS
        DT = USER_LEVELS(UTA)
        RT = DT - DBLE(INT(DT))
        N = INT(DT) + 1
        IF ( RT.GT.ZERO) THEN
          UT = UT + 1
          PARTLAYERS_OUTFLAG(UTA)  = .TRUE.
          PARTLAYERS_OUTINDEX(UTA) = UT
          PARTLAYERS_LAYERIDX(UT)  = N
          UTAU_LEVEL_MASK_UP(UTA) = N
          UTAU_LEVEL_MASK_DN(UTA) = N - 1
          PARTLAYERS_VALUES(UT)    = RT
        ELSE
          PARTLAYERS_OUTFLAG(UTA)  = .FALSE.
          PARTLAYERS_OUTINDEX(UTA) =   0
          UTAU_LEVEL_MASK_UP(UTA) = N - 1
          UTAU_LEVEL_MASK_DN(UTA) = N - 1
        ENDIF
      ENDDO
      N_PARTLAYERS = UT
      DO_PARTLAYERS = ( N_PARTLAYERS .NE. 0 )

C  Set masking and number of layer source terms
C  --------------------------------------------

C   .. for upwelling

      IF ( DO_UPWELLING ) THEN
        DO N = 1, NLAYERS
          STERM_LAYERMASK_UP(N) = .FALSE.
        ENDDO
        UTA = 1
        UT  = 1
        IF ( .NOT. PARTLAYERS_OUTFLAG(UTA) ) THEN
          N_LAYERSOURCE_UP = UTAU_LEVEL_MASK_UP(UTA) + 1
          N_ALLLAYERS_UP   = N_LAYERSOURCE_UP
        ELSE
          N_LAYERSOURCE_UP = PARTLAYERS_LAYERIDX(UT) + 1
          N_ALLLAYERS_UP   = N_LAYERSOURCE_UP - 1
        ENDIF
        DO N = NLAYERS, N_ALLLAYERS_UP, -1
          STERM_LAYERMASK_UP(N) = .TRUE.
        ENDDO
      ENDIF

C   .. for downwelling

      IF ( DO_DNWELLING ) THEN
        DO N = 1, NLAYERS
          STERM_LAYERMASK_DN(N) = .FALSE.
        ENDDO
        UTA = N_USER_LEVELS
        UT  = N_PARTLAYERS
        IF ( .NOT. PARTLAYERS_OUTFLAG(UTA) ) THEN
          N_LAYERSOURCE_DN = UTAU_LEVEL_MASK_DN(UTA)
          N_ALLLAYERS_DN   = N_LAYERSOURCE_DN
        ELSE
          N_LAYERSOURCE_DN = PARTLAYERS_LAYERIDX(UT)
          N_ALLLAYERS_DN   = N_LAYERSOURCE_DN
        ENDIF
        DO N = 1, N_ALLLAYERS_DN
          STERM_LAYERMASK_DN(N) = .TRUE.
        ENDDO
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE FINDPAR_ERROR
     &     (  ERROR, PAR_STR, STATUS )

C  include file

      INCLUDE '../includes/VLIDORT.PARS'

C  subroutine arguments

      LOGICAL       ERROR
      CHARACTER*(*) PAR_STR
      INTEGER       STATUS
      
C  local variables

      INTEGER       LEN_STRING
      CHARACTER*70  MAIL
      EXTERNAL      LEN_STRING

      IF ( ERROR ) THEN
        STATUS = VLIDORT_SERIOUS
        MAIL = 'Cannot find string: '//PAR_STR(1:LEN_STRING(PAR_STR))
        CALL VLIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
      ENDIF

C  finish

      RETURN
      END
