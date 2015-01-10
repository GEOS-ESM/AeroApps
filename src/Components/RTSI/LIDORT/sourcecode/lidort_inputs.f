C ###########################################################
C #                                                         #
C #                    THE LIDORT FAMILY                    #
C #                                                         #
C #      (LInearized Discrete Ordinate Radiative Transfer)  #
C #        --           -            -        -        -    #
C #                                                         #
C ###########################################################

C ###########################################################
C #                                                         #
C #  Author :      Robert. J. D. Spurr                      #
C #                                                         #
C #  Address :     RT Solutions, Inc.                       #
C #                9 Channing Street                        #
C #                Cambridge, MA 02138, USA                 #
C #                                                         #
C #  Tel:          (617) 492 1183                           #
C #  Email :        rtsolutions@verizon.net                 #
C #                                                         #
C #  This Version :   3.3                                   #
C #  Release Date :   September 2007                        #
C #                                                         #
C #       NEW: THERMAL SUPPLEMENT INCLUDED    (3.2)         #
C #       NEW: OUTGOING SPHERICITY CORRECTION (3.2)         #
C #       NEW: TOTAL COLUMN JACOBIANS         (3.3)         #
C #                                                         #
C ###########################################################

C    #####################################################
C    #                                                   #
C    #   This Version of LIDORT comes with a GNU-style   #
C    #   license. Please read the license carefully.     #
C    #                                                   #
C    #####################################################

C ###############################################################
C #                                                             #
C # Subroutines in this Module                                  #
C #                                                             #
C #            LIDORT_INPUT_MASTER (master), calling:           #
C #             LIDORT_INIT_CONTROL_VARS                        #
C #             LIDORT_INIT_MODEL_VARS                          #
C #             LIDORT_INIT_SURFACE_VARS                        #
C #             LIDORT_INIT_THERMAL_VARS                        #
C #             LIDORT_READ_INPUT                               #
C #                                                             #
C #            LIDORT_CHECK_INPUT                               #
C #            LIDORT_DERIVE_INPUT                              #
C #                                                             #
C ###############################################################

      SUBROUTINE LIDORT_INPUT_MASTER ( FILNAM, EFILNAM, STATUS )

C  Read all control inputs for LIDORT
C  ----------------------------------

C  include file of constants

      INCLUDE '../includes/LIDORT.PARS'

C  include file of input variables (CONTROL section)

      INCLUDE '../includes/LIDORT_INPUTS.VARS'

C  Module arguments (input filename, error filename,output status)

      CHARACTER*(*)      FILNAM
      CHARACTER*(*)      EFILNAM
      INTEGER            STATUS

C  local variables

      INTEGER            STATUS_SUB, FILUNIT, LEN_STRING
      CHARACTER*70       MAIL
      EXTERNAL           LEN_STRING

C  initialize status

      STATUS = LIDORT_SUCCESS

C  initialize error file

      LIDORT_ERROR_FILENAME = EFILNAM

      LIDORT_ERROR_INIT     = .TRUE.

C  initialize variables

      CALL LIDORT_INIT_CONTROL_VARS

      CALL LIDORT_INIT_MODEL_VARS

      CALL LIDORT_INIT_SURFACE_VARS

      CALL LIDORT_INIT_THERMAL_VARS

C  Open file

      FILUNIT = LIDORT_INUNIT
      OPEN(LIDORT_INUNIT,FILE=FILNAM,ERR=300,STATUS='OLD')

C  read standard inputs

      CALL LIDORT_READ_INPUT ( STATUS_SUB )
      IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
        STATUS = LIDORT_SERIOUS
        MAIL = 'Error from LIDORT_INPUTREAD'
        CALL LIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
        CLOSE(FILUNIT)
        RETURN
      ENDIF

C  normal return

      CLOSE(FILUNIT)
      RETURN
      
C  Open file error

300   CONTINUE
      STATUS = LIDORT_SERIOUS
      MAIL = 'openfile failure for '//FILNAM(1:LEN_STRING(FILNAM))
      CALL LIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
      CLOSE(FILUNIT)
      RETURN

C  Finish

      END


      SUBROUTINE LIDORT_INIT_CONTROL_VARS

C  Initialises all control inputs for LIDORT
C  ------------------------------------------

C  include file of constants

      INCLUDE '../includes/LIDORT.PARS'

C  include file of input variables (CONTROL section)

      INCLUDE '../includes/LIDORT_INPUTS.VARS'

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

C  single scatter correction now active version 2.0
C  direct beam correction now active version 3.0
C  Outgoing sphericity correction 23 January 2007, version 3.2

c      SAVE_LAYER_MSST    = .FALSE.

C  solar beam options
C  ------------------

      DO_SOLAR_SOURCES       = .FALSE.
      DO_CLASSICAL_SOLUTION  = .FALSE.
      DO_PLANE_PARALLEL      = .FALSE.
      DO_REFRACTIVE_GEOMETRY = .FALSE.
      DO_CHAPMAN_FUNCTION    = .FALSE.

C  special options for Azimuth control
C  -----------------------------------

      DO_RAYLEIGH_ONLY  = .FALSE.
      DO_ISOTROPIC_ONLY = .FALSE.

      DO_NO_AZIMUTH      = .FALSE.
      DO_ALL_FOURIER     = .FALSE.

C  Performance enhancements
C  ========================

C  delta-M scaling now active version 3.

      DO_DELTAM_SCALING  = .FALSE.

C  New flags for version 3. RTSolutions 4/11/05

      DO_SOLUTION_SAVING = .FALSE.
      DO_BVP_TELESCOPING = .FALSE.

C  write options
C  =============

      DO_INFO_WRITE        = .FALSE.
      DO_DEBUG_WRITE       = .FALSE.
      DO_WRITE_INPUT       = .FALSE.
      DO_WRITE_RESULTS     = .FALSE.
      DO_WRITE_SCENARIO    = .FALSE.
      DO_WRITE_FOURIER     = .FALSE.

C  User-output options
C  -------------------

      DO_ADDITIONAL_MVOUT = .FALSE.
      DO_MVOUT_ONLY       = .FALSE.

      DO_QUAD_OUTPUT  = .FALSE.
      DO_USER_STREAMS = .FALSE.
      DO_USER_TAUS    = .FALSE.

      DO_UPWELLING = .FALSE.
      DO_DNWELLING = .FALSE.

C  finish

      RETURN
      END

C

      SUBROUTINE LIDORT_INIT_SURFACE_VARS

C  include files
C  -------------

C  include file of constants

      INCLUDE '../includes/LIDORT.PARS'

C  include file of input variables (Surface section)

      INCLUDE '../includes/LIDORT_INPUTS.VARS'

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
      DO_SHADOW_EFFECT = .FALSE.
      DO_GLITTER_DBMS  = .FALSE.

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

      SUBROUTINE LIDORT_INIT_THERMAL_VARS

C  include files
C  -------------

C  include file of constants

      INCLUDE '../includes/LIDORT.PARS'

C  include file of input variables (Surface section)

      INCLUDE '../includes/LIDORT_INPUTS.VARS'

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

C  This flag introduced 31 July 2007

      DO_THERMAL_TRANSONLY = .FALSE.

C  finish

      RETURN
      END

C

      SUBROUTINE LIDORT_INIT_MODEL_VARS

C  Initialises all file-read model inputs for LIDORT
C  --------------------------------------------------

C  include file of constants

      INCLUDE '../includes/LIDORT.PARS'

C  include file of input variables (fileread section)

      INCLUDE '../includes/LIDORT_INPUTS.VARS'

C  local variables

      INTEGER             I

C  basic integer inputs

      NSTREAMS = 0
      NLAYERS  = 0
      NFINELAYERS    = 0
      NMOMENTS_INPUT = 0

C  accuracy

      LIDORT_ACCURACY  = ZERO
      ZENITH_TOLERANCE  = ZERO

C  solar beam

      FLUX_FACTOR      = ZERO
      NBEAMS   = 0
      DO I = 1, MAXBEAMS
        BEAM_SZAS(I) = ZERO
      ENDDO

C  Pseudo-spherical

      EARTH_RADIUS      = ZERO
      RFINDEX_PARAMETER = ZERO

C  Geometry specification height
C   (New, 06 August 2007)

      GEOMETRY_SPECHEIGHT = ZERO

C  user angles

      N_USER_STREAMS = 0
      N_USER_RELAZMS = 0
      N_OUT_USERTAUS = 0
      DO I = 1, MAX_USER_STREAMS
        USER_ANGLES_INPUT(I) = ZERO
      ENDDO
      DO I = 1, MAX_USER_RELAZMS
        USER_RELAZMS(I) = ZERO
      ENDDO
      DO I = 1, MAX_OUT_USERTAUS
        USER_TAUS_INPUT(I) = ZERO
      ENDDO

C finish

      RETURN
      END

C

      SUBROUTINE LIDORT_READ_INPUT ( STATUS )

C  Read all control inputs for LIDORT
C  ----------------------------------

C  include file of constants

      INCLUDE '../includes/LIDORT.PARS'

C  include file of input variables to be filled out with input data

      INCLUDE '../includes/LIDORT_INPUTS.VARS'

C  Module arguments (output status)

      INTEGER         STATUS

C  local variables

      CHARACTER*8     PREFIX
      PARAMETER     ( PREFIX = 'LIDORT -' )
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

      STATUS = LIDORT_SUCCESS
      ERROR  = .FALSE.

C  file unit

      FILUNIT = LIDORT_INUNIT

C  1. read all variables in include file LIDORT_CONTROL.VARS
C  =========================================================

C  operation modes
C  ---------------

C  full Radiance calculation

      PAR_STR = 'Do full Radiance calculation?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR))
     &     READ (FILUNIT,*,ERR=998) DO_FULLRAD_MODE
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  Corrections

      PAR_STR = 'Do nadir single scatter correction?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR))
     &     READ (FILUNIT,*,ERR=998) DO_SSCORR_NADIR
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      PAR_STR = 'Do outgoing single scatter correction?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR))
     &     READ (FILUNIT,*,ERR=998) DO_SSCORR_OUTGOING
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      PAR_STR = 'Do direct beam correction?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR))
     &     READ (FILUNIT,*,ERR=998) DO_DBCORRECTION
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  Full-up single scatter calculation

      PAR_STR = 'Do Full-up single scatter calculation?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR))
     &     READ (FILUNIT,*,ERR=998) DO_SSFULL
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  multiple scatter source function output control
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

C  Basic control

      PAR_STR = 'Include solar sources?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR))
     &     READ (FILUNIT,*,ERR=998) DO_SOLAR_SOURCES
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  Other control options for the solar sources

      IF ( DO_SOLAR_SOURCES ) THEN

C  Solution method

       PAR_STR = 'Use classical beam solution?'
       IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR))
     &     READ (FILUNIT,*,ERR=998) DO_CLASSICAL_SOLUTION
       CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

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
     &    READ (FILUNIT,*,ERR=998) DO_REFRACTIVE_GEOMETRY
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

      PAR_STR='Isotropic atmosphere only?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_ISOTROPIC_ONLY
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  No azimuth dependence (TRUE means Fourier m = 0 only )

      PAR_STR = 'No azimuth dependence in the solution?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_NO_AZIMUTH
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  All possible Fourier components (2N-1). Debug only

      PAR_STR = 'Compute all Fourier components?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_ALL_FOURIER
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

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

      PAR_STR = 'Include quadrature angles in output?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_QUAD_OUTPUT
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      PAR_STR = 'User-defined stream angles?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &   READ (FILUNIT,*,ERR=998) DO_USER_STREAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      PAR_STR = 'User-defined optical depths?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &   READ (FILUNIT,*,ERR=998) DO_USER_TAUS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      PAR_STR = 'Layer boundary optical depths?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &   READ (FILUNIT,*,ERR=998) DO_LBOUND_TAUS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

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

C  2.  read variables in LIDORT_INPUT.VARS, group B
C  ================================================

C  streams/layers/moments (INTEGER input)

      PAR_STR = 'Number of half-space streams'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) NSTREAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      PAR_STR = 'Number of atmospheric layers'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) NLAYERS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      IF ( DO_SSCORR_OUTGOING ) THEN
        PAR_STR =
     &   'Number of fine layers (outgoing sphericity correction only)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &       READ (FILUNIT,*,ERR=998) NFINELAYERS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )
      ENDIF

      PAR_STR = 'Number of input Legendre moments'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) NMOMENTS_INPUT
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  All numbers are now checked against maximum diensions

      IF ( NSTREAMS .GT. MAXSTREAMS ) THEN
        MAIL   = 'Number of half-space streams > maximum dimension'
        ACTION = ' Re-set input value'
        STATUS = LIDORT_SERIOUS
        CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        RETURN
      ENDIF
      IF ( NLAYERS .GT. MAXLAYERS ) THEN
        MAIL   = 'Number of layers > maximum dimension'
        ACTION = ' Re-set input value'
        STATUS = LIDORT_SERIOUS
        CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        RETURN
      ENDIF
      IF ( DO_SSCORR_OUTGOING ) THEN
        IF ( NFINELAYERS .GT. MAXFINELAYERS ) THEN
          MAIL   = 'Number of fine layers > maximum dimension'
          ACTION = 'Re-set input value'
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          RETURN
        ENDIF
      ENDIF
      IF ( NMOMENTS_INPUT .GT. MAXMOMENTS_INPUT ) THEN
        MAIL   = 'Number of Legendre moments > maximum dimension'
        ACTION = ' Re-set input value'
        STATUS = LIDORT_SERIOUS
        CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        RETURN
      ENDIF

C  accuracy input
C  --------------

      PAR_STR = 'Fourier series convergence'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) LIDORT_ACCURACY
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      PAR_STR = 'Zenith tolerance level'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) ZENITH_TOLERANCE
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  Solar Beams input
C  -----------------

C  Flux constant

      PAR_STR = 'Solar flux constant'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) FLUX_FACTOR
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  number of Solar zenith angles

      PAR_STR = 'Number of solar zenith angles'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) NBEAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  check not exceeding dimensioned number

      IF ( NBEAMS .GT. MAXBEAMS ) THEN
        MAIL   = 'Number of solar zenith angles > maximum dimension'
        ACTION = ' Re-set input value'
        STATUS = LIDORT_SERIOUS
        CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        RETURN
      ENDIF

C  TOA solar zenith angle inputs

      PAR_STR = 'TOA solar zenith angles (degrees)'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
        DO I = 1, NBEAMS
          READ (FILUNIT,*,ERR=998) BEAM_SZAS(I)
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

C  1. always need some azimuth angles if Azimuth flag set
C      --- Set number of azimuths locally to 1 if flag not set

      IF ( .NOT. DO_NO_AZIMUTH ) THEN

        PAR_STR = 'Number of user-defined relative azimuth angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) N_USER_RELAZMS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        IF ( N_USER_RELAZMS .GT. MAX_USER_RELAZMS ) THEN
          MAIL   = 'Number of relative azimuths > maximum dimension'
          ACTION = ' Re-set input value'
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          RETURN
        ENDIF

        PAR_STR = 'User-defined relative azimuth angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_USER_RELAZMS
            READ (FILUNIT,*,ERR=998) USER_RELAZMS(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      ELSE

        N_USER_RELAZMS = 1

      ENDIF

C  2. User defined stream angles (should be positive)

      IF ( DO_USER_STREAMS ) THEN

        PAR_STR = 'Number of user-defined stream angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) N_USER_STREAMS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        IF ( N_USER_STREAMS .GT. MAX_USER_STREAMS ) THEN
          MAIL   = 'Number of user streams > maximum dimension'
          ACTION = ' Re-set input value'
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          RETURN
        ENDIF

        PAR_STR = 'User-defined stream angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_USER_STREAMS
            READ (FILUNIT,*,ERR=998) USER_ANGLES_INPUT(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      ENDIF

C  3. User defined optical depths (should be positive)

      IF ( DO_USER_TAUS ) THEN

        PAR_STR = 'Number of user-defined optical depths'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &        READ (FILUNIT,*,ERR=998) N_OUT_USERTAUS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        IF ( N_OUT_USERTAUS .GT. MAX_OUT_USERTAUS ) THEN
          MAIL=
     &      'Number of user-defined optical depths > maximum dimension'
          ACTION = ' Re-set input value'
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          RETURN
        ENDIF

        PAR_STR = 'User-defined optical depths'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_OUT_USERTAUS
            READ (FILUNIT,*,ERR=998) USER_TAUS_INPUT(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C   ( override if layer boundaries are required )

      ELSE IF ( DO_LBOUND_TAUS ) THEN

        PAR_STR = 'Number of layer boundary optical depths'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &        READ (FILUNIT,*,ERR=998) N_OUT_USERTAUS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        IF ( N_OUT_USERTAUS .GT. MAX_OUT_USERTAUS ) THEN
          MAIL=
     &    'Number of layer boundary optical depths > maximum dimension'
          ACTION = ' Re-set input value'
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          RETURN
        ENDIF

        PAR_STR = 'Which layer boundaries for optical depth'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_OUT_USERTAUS
            READ (FILUNIT,*,ERR=998) LBOUND_TAUS_INPUT(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      ENDIF

C  3. read all variables in LIDORT_INPUT.VARS, group D
C  ===================================================

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

C  set defaults for Lambertian case (Derived input)

      IF ( DO_LAMBERTIAN_SURFACE ) THEN
        N_BRDF_KERNELS = 1
        LAMBERTIAN_KERNEL_FLAG(1) = .TRUE.
        BRDF_FACTORS(1) = LAMBERTIAN_ALBEDO
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
          ACTION = ' Re-set input value'
          STATUS  = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          RETURN
        ENDIF

C  number of BRDF azimuth streams, check this value

        PAR_STR = 'Number of bidirectional reflectance streams'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &       READ (FILUNIT,*,ERR=998) NSTREAMS_BRDF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        IF ( NSTREAMS_BRDF .GT. MAXSTREAMS_BRDF ) THEN
          MAIL = 'Number of BRDF streams > maximum dimension'
          ACTION = ' Re-set input value'
          STATUS  = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
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

C  Shadowing input (for Cox-Munk type)

        DO I = 1, N_BRDF_KERNELS
         IF ( BRDF_NAMES(I) .EQ. 'Cox-Munk  ' ) THEN
           PAR_STR = 'Do shadow effect for glitter kernels?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
            READ (FILUNIT,*,ERR=998)DO_SHADOW_EFFECT
           ENDIF
         ENDIF
        ENDDO
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  Multiple reflectance DB correction (for Cox-Munk type)

        DO I = 1, N_BRDF_KERNELS
         IF ( BRDF_NAMES(I) .EQ. 'Cox-Munk  ' ) THEN
           PAR_STR = 'Do multiple reflectance for glitter kernels?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
            READ (FILUNIT,*,ERR=998)DO_GLITTER_DBMS
           ENDIF
         ENDIF
        ENDDO
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      ENDIF

C  4. read all variables in LIDORT_INPUT.VARS, group E
C  ===================================================

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
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
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
      STATUS = LIDORT_SERIOUS
      MAIL = 'read failure for '//PAR_STR(1:LEN_STRING(PAR_STR))
      CALL LIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_CHECK_INPUT (STATUS)

C  check inputs (both file-read and derived)

      INCLUDE '../includes/LIDORT.PARS'

C  include file with input variables to be checked out

      INCLUDE '../includes/LIDORT_INPUTS.VARS'

C  Module output

      INTEGER          STATUS

C  local variables

      CHARACTER*70     MAIL
      CHARACTER*70     ACTION

      INTEGER          I, K, L, N, UTA, NSTART, N1
      INTEGER          NALLSTREAMS, N_OFFGRID_LOCAL
      DOUBLE PRECISION TAU, TAUATM, LOCALTAU(0:MAXLAYERS)
      CHARACTER*3      C3
      CHARACTER*2      C2
      LOGICAL          LOOP

C  BRDF Kernel names (check)

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
     &     'Gisssoil  ',  
     &     'Gisssnow  '/

C  Initialize output status

      STATUS = LIDORT_SUCCESS

C  Check top level options, set warnings
C  =====================================

C  Check thermal or Solar sources present

      IF ( .NOT.DO_SOLAR_SOURCES.AND..NOT.DO_THERMAL_EMISSION ) THEN
        MAIL = 'Bad input: No solar or thermal sources'
        ACTION = 'Abort: must set one of the source flags!'
        STATUS = LIDORT_SERIOUS
        CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
      ENDIF

C  Switch off several flags with thermal-only option
C    Set default, regardless of whether solar sources are on.

      IF ( .NOT.DO_SOLAR_SOURCES.AND.DO_THERMAL_EMISSION ) THEN
       IF ( DO_SSCORR_NADIR ) THEN
         MAIL = 'Switch off SS correction, not needed for thermal-only'
         ACTION = 'Warning: SS correction flag turned off internally'
         STATUS = LIDORT_WARNING
         DO_SSCORR_NADIR = .FALSE.
         CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
       IF ( DO_DBCORRECTION ) THEN
         MAIL = 'Switch off DB correction, not needed for thermal-only'
         ACTION = 'Warning: DB correction flag turned off internally'
         STATUS = LIDORT_WARNING
         DO_DBCORRECTION = .FALSE.
         CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
      ENDIF

C  Set number of sources NBEAMS to 1 for the thermal-only default

      IF ( .NOT.DO_SOLAR_SOURCES.AND.DO_THERMAL_EMISSION ) THEN
c        MAIL = 'Bad input: NBEAMS set to 1 for thermal-only'
c        ACTION = 'Warning: NBEAMS set to 1 internally'
c        STATUS = LIDORT_WARNING
        NBEAMS = 1
c        CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
      ENDIF

C  Check dimensions (safety first). All errors are FATAL
C  -----------------------------------------------------

C  1. Brdf dimensioning

      IF ( .NOT. DO_LAMBERTIAN_SURFACE ) THEN

C  number of kernels, check this value

        IF ( N_BRDF_KERNELS .GT. MAX_BRDF_KERNELS ) THEN
          MAIL   = 'Number of BRDF kernels > maximum dimension'
          ACTION = 'Re-set input value'
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          RETURN
        ENDIF

C  number of BRDF azimuth streams, check this value

        IF ( NSTREAMS_BRDF .GT. MAXSTREAMS_BRDF ) THEN
          MAIL   = 'Number of BRDF streams > maximum dimension'
          ACTION = 'Re-set input value'
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          RETURN
        ENDIF

      ENDIF

C  2. Basic layer/streams/moments

      IF ( NSTREAMS .GT. MAXSTREAMS ) THEN
        MAIL   = 'Number of half-space streams > maximum dimension'
        ACTION = 'Re-set input value'
        STATUS = LIDORT_SERIOUS
        CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        RETURN
      ENDIF

      IF ( NLAYERS .GT. MAXLAYERS ) THEN
        MAIL   = 'Number of layers > maximum dimension'
        ACTION = 'Re-set input value'
        STATUS = LIDORT_SERIOUS
        CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        RETURN
      ENDIF

      IF ( DO_SSCORR_OUTGOING ) THEN
        IF ( NFINELAYERS .GT. MAXFINELAYERS ) THEN
          MAIL   = 'Number of fine layers > maximum dimension'
          ACTION = 'Re-set input value'
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          RETURN
        ENDIF
      ENDIF

      IF ( NMOMENTS_INPUT .GT. MAXMOMENTS_INPUT ) THEN
        MAIL   = 'Number of Legendre moments > maximum dimension'
        ACTION = 'Re-set input value'
        STATUS = LIDORT_SERIOUS
        CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        RETURN
      ENDIF

C  3. Basic geometry

      IF ( NBEAMS .GT. MAXBEAMS ) THEN
        MAIL   = 'Number of solar zenith angles > maximum dimension'
        ACTION = 'Re-set input value'
        STATUS = LIDORT_SERIOUS
        CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        RETURN
      ENDIF

      IF ( DO_USER_STREAMS ) THEN
        IF ( N_USER_STREAMS .GT. MAX_USER_STREAMS ) THEN
          MAIL   = 'Number of user streams > maximum dimension'
          ACTION = 'Re-set input value'
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          RETURN
        ENDIF
      ENDIF

      IF ( N_USER_RELAZMS .GT. MAX_USER_RELAZMS ) THEN
        MAIL   = 'Number of relative azimuths > maximum dimension'
        ACTION = 'Re-set input value'
        STATUS = LIDORT_SERIOUS
        CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        RETURN
      ENDIF

C  4. optical depth output

      IF ( DO_LBOUND_TAUS ) THEN
        IF ( N_OUT_USERTAUS .GT. MAX_OUT_USERTAUS ) THEN
          MAIL=
     &    'Number of layer boundary optical depths > maximum dimension'
          ACTION = 'Re-set input value'
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          RETURN
        ENDIF
      ENDIF

      IF ( DO_USER_TAUS ) THEN
        IF ( N_OUT_USERTAUS .GT. MAX_OUT_USERTAUS ) THEN
          MAIL=
     &      'Number of user-defined optical depths > maximum dimension'
          ACTION = 'Re-set input value'
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          RETURN
        ENDIF
      ENDIF

C  check inputs (both file-read and derived)
C  -----------------------------------------

C  Check Chapman function options

      IF ( DO_SOLAR_SOURCES ) THEN
       IF ( .NOT. DO_CHAPMAN_FUNCTION ) THEN
        IF ( DO_PLANE_PARALLEL ) THEN
          MAIL   = 'Chapman Function not set, plane parallel'
          ACTION = 'Warning: Chapman function set internally'
          STATUS = LIDORT_WARNING
          DO_CHAPMAN_FUNCTION = .TRUE.
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ELSE
          MAIL   = 'Chapman Function not set, pseudo-spherical'
          ACTION = 'Have you set the DELTAU_SLANT_INPUT values?'
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
       ENDIF
      ENDIF

C  Check sphericity corrections....Cannot both be turned on
C    --------- New code 31 January 2007

      IF ( DO_SOLAR_SOURCES ) THEN
       IF ( DO_SSCORR_NADIR .and. DO_SSCORR_OUTGOING ) THEN
        MAIL   = 'Cannot have both single scatter corrections on'
        ACTION = 'Turn off DO_SSCORR_NADIR and/or DO_SSCORR_OUTGOING'
        STATUS = LIDORT_SERIOUS
        CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
       ENDIF     
c       IF ( DO_SSCORR_OUTGOING ) THEN
c        IF ( DO_USER_TAUS ) THEN
c          MAIL   = 'SS correction outgoing only for LBOUND_TAUS'
c          ACTION = 'Turn off DO_USER_TAUS, use only DO_LBOUND_TAUS'
c          STATUS = LIDORT_SERIOUS
c          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
c        ENDIF
c        IF ( DO_DNWELLING ) THEN
c          MAIL   = 'SS correction outgoing only for Upwelling'
c          ACTION = 'Turn off DO_DNWELLING'
c          STATUS = LIDORT_SERIOUS
c          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
c        ENDIF
c       ENDIF
      ENDIF

C  Check beam mode operation. Warning

      IF ( DO_SOLAR_SOURCES ) THEN
       IF ( DO_PLANE_PARALLEL ) THEN
        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
         MAIL =
     &  'Bad input: plane-parallel and refractive flags both set'
         ACTION = 'Warning: turn off Refraction internally'
         STATUS = LIDORT_WARNING
         DO_REFRACTIVE_GEOMETRY = .FALSE.
         CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
       ENDIF
      ENDIF

C  check consistency of mean value input control
C  ---------------------------------------------

      IF ( DO_ADDITIONAL_MVOUT ) THEN
        IF ( DO_MVOUT_ONLY ) THEN
          MAIL   = 'Bad input: Cannot have both mean-value flags set'
          ACTION = 'Warning: disable DO_MVOUT_ONLY flag internally'
          STATUS = LIDORT_WARNING
          DO_MVOUT_ONLY = .FALSE.
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
      ENDIF

      IF ( .NOT.DO_ADDITIONAL_MVOUT ) THEN
        IF ( DO_MVOUT_ONLY ) THEN
          IF ( .NOT. DO_NO_AZIMUTH ) THEN
            MAIL   =
     $          'Bad input: Mean-value option requires NO azimuth'
            ACTION = 'Warning: DO_NO_AZIMUTH flag was set internally'
            STATUS = LIDORT_WARNING
            DO_NO_AZIMUTH = .TRUE.
            CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          ENDIF
          IF ( DO_USER_STREAMS ) THEN
            MAIL   =
     &           'Bad input: Mean-value option needs quadratures only'
            ACTION = 'Warning: DO_USER_STREAMS flag disabled internally'
            STATUS = LIDORT_WARNING
            DO_USER_STREAMS = .FALSE.
            CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          ENDIF
        ENDIF
      ENDIF

C  check consistency of multiple scatter source term output control
C   ---SPecialist options. Commented out 2 April 2007
c      IF ( SAVE_LAYER_MSST ) THEN
c        IF ( .NOT. DO_LBOUND_TAUS ) THEN
c          MAIL =
c     & 'Bad input: MSCAT. source term - layer boundary TAU flag not set'
c          ACTION = 'Check DO_LBOUND_TAUS and SAVE_LAYER_MSST flags'
c          STATUS = LIDORT_SERIOUS
c          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
c        ENDIF
c        IF ( .NOT. DO_USER_STREAMS ) THEN
c          MAIL   =
c     &   'Bad input: MSCAT. source term - user_streams flag not set'
c          ACTION = 'Check DO_USER_STREAMS and SAVE_LAYER_MSST flags'
c          STATUS = LIDORT_SERIOUS
c          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
c        ENDIF
c      ENDIF

C  Check consistency of Quadrature output flag
C    - Quadrature output only with Uncorrected Full-radiance calculation.
C    - All warnings. DO_QUAD_OUTPUT disabled internally

      IF ( DO_QUAD_OUTPUT ) THEN
        IF ( .NOT. DO_FULLRAD_MODE ) THEN
          MAIL =
     &     'Bad input: Quadrature output not with Full radiance mode'
          ACTION = 'Warning: turn off DO_QUAD_OUTPUT internally'
          STATUS = LIDORT_WARNING
          DO_QUAD_OUTPUT = .FALSE.
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ELSE
         IF ( DO_SOLAR_SOURCES ) THEN
          IF ( DO_SSCORR_NADIR .OR. DO_SSCORR_OUTGOING ) THEN
            MAIL =
     &       'Bad input: Quadrature output not with SS corrections'
            ACTION = 'Warning: turn off DO_QUAD_OUTPUT internally'
            STATUS = LIDORT_WARNING
            DO_QUAD_OUTPUT = .FALSE.
            CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          ENDIF
          IF ( DO_DBCORRECTION ) THEN
            MAIL =
     &       'Bad input: Quadrature output not with DB correction'
            ACTION = 'Warning: turn off DO_QUAD_OUTPUT internally'
            STATUS = LIDORT_WARNING
            DO_QUAD_OUTPUT = .FALSE.
            CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          ENDIF
         ENDIF
        ENDIF
      ENDIF

C  Check consistency of BVP_TELESCOPING and SOLUTION_SAVING flags
C  ---Warning. Set solution-saving internally

      IF (DO_BVP_TELESCOPING.AND..NOT.DO_SOLUTION_SAVING) THEN
        MAIL   =
     &  'Bad input: BVP telescoping -> solution saving must be set'
        ACTION = 'Warning:  Solution saveing was set internally'
        STATUS = LIDORT_WARNING
        DO_SOLUTION_SAVING = .TRUE.
        CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
      ENDIF

C  Check consistency of Rayleigh-only and Isotropic-only cases

      IF ( DO_RAYLEIGH_ONLY .AND. DO_ISOTROPIC_ONLY ) THEN
        MAIL   =
     &   'Bad input: Isotropic_only & Rayleigh-only flags both set'
        ACTION =
     &     'Check DO_RAYLEIGH_ONLY and DO_ISOTROPIC_ONLY flags'
        STATUS = LIDORT_SERIOUS
        CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
      ENDIF

C  -----------Note the following in the scalar code -------------------

C  no Delta-M scaling with Rayleigh only
C   ---Warning. Turn off delta-M scaling.

      IF ( DO_RAYLEIGH_ONLY ) THEN
        IF ( DO_DELTAM_SCALING ) THEN
          MAIL  ='Bad input: No delta-M scaling with Rayleigh-only'
          ACTION = 'Warning: DO_DELTAM_SCALING turned off internally'
          STATUS = LIDORT_WARNING
          DO_DELTAM_SCALING = .FALSE.
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
      ENDIF

C  no Delta-M scaling with Isotropic only
C   ---Warning. Turn off delta-M scaling.

      IF ( DO_ISOTROPIC_ONLY ) THEN
        IF ( DO_DELTAM_SCALING ) THEN
          MAIL  ='Bad input: No delta-M scaling with Isotropic-only'
          ACTION = 'Warning: DO_DELTAM_SCALING turned off internally'
          STATUS = LIDORT_WARNING
          DO_DELTAM_SCALING = .FALSE.
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
      ENDIF

C---------------------------------------------------------------------

C  check optical depth flags

      IF ( DO_USER_TAUS .AND. DO_LBOUND_TAUS ) THEN
        MAIL   = 'Bad input: both optical depth flags cannot be set'
        ACTION = 'Check DO_USER_TAUS and DO_LBOUND_TAUS flags'
        STATUS = LIDORT_SERIOUS
        CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
      ELSE IF ( .NOT. DO_USER_TAUS .AND. .NOT. DO_LBOUND_TAUS ) THEN
        MAIL   = 'Bad input: Neither optical depth flags are set'
        ACTION = 'Check DO_USER_TAUS and DO_LBOUND_TAUS flags'
        STATUS = LIDORT_SERIOUS
        CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
      ENDIF
        
C  check directional input

      IF ( .NOT.DO_UPWELLING .AND. .NOT. DO_DNWELLING ) THEN
        MAIL   = 'Bad input: no directional input is set'
        ACTION = 'Check DO_UPWELLING & DO_DNWELLING: one must be set!'
        STATUS = LIDORT_SERIOUS
        CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
      ENDIF

C  check number of input Legendre moments
C  ======================================

C  (general scattering case)

      IF ( .NOT.DO_RAYLEIGH_ONLY.AND..NOT.DO_ISOTROPIC_ONLY
     &       .AND..NOT. DO_SSFULL ) THEN

        IF ( DO_DELTAM_SCALING ) THEN
          IF ( NMOMENTS_INPUT.LT.2*NSTREAMS ) THEN
            MAIL =
     &      'Bad input: Fewer than 2N moments with delta-M'
            ACTION = 
     &       'Warning: Re-set NMOMENTS_INPUT to 2N internally'
            STATUS = LIDORT_WARNING
            NMOMENTS_INPUT = 2*NSTREAMS
            CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          ENDIF
        ELSE
          IF ( DO_SSCORR_NADIR .OR. DO_SSCORR_OUTGOING ) THEN
            IF ( NMOMENTS_INPUT.LT.2*NSTREAMS-1 ) THEN
              MAIL =
     &      'Bad input: Fewer than 2N-1 moments without delta-M'
              ACTION = 
     &       'Warning: Re-set NMOMENTS_INPUT to 2N-1 internally'
              STATUS = LIDORT_WARNING
              NMOMENTS_INPUT = 2*NSTREAMS - 1
              CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
            ENDIF
          ENDIF
        ENDIF

      ELSE

C  Checks for Rayleigh only option
C   All warnings.

        IF ( DO_RAYLEIGH_ONLY ) THEN

          IF ( NMOMENTS_INPUT.NE.2 ) THEN
            MAIL   = 'Bad input: Rayleigh-only phase momemts not = 2'
            ACTION = 'Warning: Set NMOMENTS_INPUT = 2 internally'
            STATUS = LIDORT_WARNING
            NMOMENTS_INPUT = 2
            CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          ENDIF

          IF ( DO_BVP_TELESCOPING ) THEN
            MAIL   = 
     &   'Bad input: Bvp telescoping not possible, Rayleigh only'
            ACTION = 'Warning: Turn off BVP_TELESCOPING internally'
            STATUS = LIDORT_WARNING
            DO_BVP_TELESCOPING = .FALSE.
            CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          ENDIF

          IF ( DO_SOLUTION_SAVING ) THEN
            MAIL   =
     &   'Bad input: Solution saving not possible, Rayleigh only'
            ACTION = 'Warning: Turn off SOLUTION_SAVING internally'
            STATUS = LIDORT_WARNING
            DO_SOLUTION_SAVING = .FALSE.
            CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          ENDIF

        ENDIF

C  Checks for Isotropic only option
C   All warning messages

        IF ( DO_ISOTROPIC_ONLY ) THEN

          IF ( NMOMENTS_INPUT.NE.0 ) THEN
            MAIL   =
     &       'Bad input: No phase function moments for isotropic-only'
            ACTION = 'Warning: NMOMENTS_INPUT = 0 was set internally'
            STATUS = LIDORT_WARNING
            NMOMENTS_INPUT = 0
            CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          ENDIF

          IF ( DO_BVP_TELESCOPING ) THEN
            MAIL   = 
     &   'Bad input: Bvp telescoping not possible, Isotropic only'
            ACTION = 'Warning: Turn off BVP_TELESCOPING internally'
            STATUS = LIDORT_WARNING
            DO_BVP_TELESCOPING = .FALSE.
            CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          ENDIF

          IF ( DO_SOLUTION_SAVING ) THEN
            MAIL   = 
     &  'Bad input: Solution saving not possible, Isotropic only'
            ACTION = 'Warning" Turn off SOLUTION_SAVING internally'
            STATUS = LIDORT_WARNING
            DO_SOLUTION_SAVING = .FALSE.
            CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          ENDIF

        ENDIF

      ENDIF

C  reset solution saving and bvp telescoping flags

      DO_SOLUTION_SAVING =  ( DO_SOLUTION_SAVING .AND.
     &     ((.NOT.DO_RAYLEIGH_ONLY).OR.(.NOT.DO_ISOTROPIC_ONLY)) )

      DO_BVP_TELESCOPING =  ( DO_BVP_TELESCOPING .AND.
     &     ((.NOT.DO_RAYLEIGH_ONLY).OR.(.NOT.DO_ISOTROPIC_ONLY)) )

C  BVP telescoping doesn't work with non-Lambertian surfaces
C   Reason, not yet coded. Theory already worked out.
C    ---WARNING. BVP telescoping Flag turned off

      IF (  DO_BVP_TELESCOPING ) THEN
        IF ( .NOT. DO_LAMBERTIAN_SURFACE ) THEN
          MAIL   = 'BVP telescoping must be disabled, non-Lambertian'
          ACTION = 'Warning: DO_BVP_TELESCOPING turned off internally'
          STATUS = LIDORT_WARNING
          DO_BVP_TELESCOPING = .FALSE.
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
      ENDIF
   
C  Check azimuth-only conditions
C  =============================

C  Check no-Azimuth flag
C    ---WARNING. Do-no-Azimuth Flag turned on

      IF ( .NOT.DO_NO_AZIMUTH ) THEN
        IF ( DO_USER_STREAMS. AND. N_USER_STREAMS.EQ.1 ) THEN
          IF ( USER_ANGLES_INPUT(1) .EQ. ZERO ) THEN
            MAIL   ='Bad input: zenith-sky output requires no azimuth'
            ACTION ='Warning: DO_NO_AZIMUTH flag set true internally'
            STATUS = LIDORT_WARNING
            CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          ENDIF
        ENDIF
      ENDIF

C  Checks for Isotropic only option
C    ---WARNING. Do-no-Azimuth Flag turned on

      IF ( DO_ISOTROPIC_ONLY ) THEN
        IF ( .NOT.DO_NO_AZIMUTH ) THEN
          MAIL   =
     &       'Bad input: no azimuth dependence for isotropic_only'
          ACTION = 'Warning: DO_NO_AZIMUTH turned on internally'
          STATUS = LIDORT_WARNING
          DO_NO_AZIMUTH = .TRUE.
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
      ENDIF

C  Check single scattering correction and Do no Azimuth
C    ---WARNING. Do-no-Azimuth Flag turned off

      IF ( DO_SSCORR_NADIR .OR. DO_SSCORR_OUTGOING ) THEN
        IF ( DO_NO_AZIMUTH ) THEN
          MAIL   =
     &       'Bad input: need azimuth dependence for SS corrections'
          ACTION = 'Warning: DO_NO_AZIMUTH turned off internally'
          STATUS = LIDORT_WARNING
          DO_NO_AZIMUTH = .FALSE.
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
      ENDIF

C  Check: OLD single scattering correction and Do Rayleigh
C    ---WARNING. SS Flag turned off
C  Check only required for the diffuse field calculations (Version 2.3)

      IF ( DO_SSCORR_NADIR ) THEN
        IF ( DO_RAYLEIGH_ONLY .AND. .NOT. DO_SSFULL ) THEN
          MAIL   = 'Bad input: No SS correction for Rayleigh only'
          ACTION = 'Warning: DO_SSCORR_NADIR turned off internally'
          STATUS = LIDORT_WARNING
          DO_SSCORR_NADIR = .FALSE.
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
      ENDIF

C  Full-up single scatter, enabled 25 September 2007.
C   Single scatter corrections must be turned on

      IF ( DO_SSFULL ) THEN
        IF ( .not.DO_SSCORR_NADIR.and..not.DO_SSCORR_OUTGOING ) THEN        
          MAIL   = 'Bad input: Full SS, must have one SSCORR flag set'
          ACTION = 'Full SS: default to use outgoing SS correction'
          STATUS = LIDORT_WARNING
          DO_SSCORR_NADIR    = .FALSE.
          DO_SSCORR_OUTGOING = .TRUE.
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
      ENDIF

C  Full-up single scatter, enabled 25 September 2007.
C   Diffuse-field Delta-M scaling must be turned off

      IF ( DO_SSFULL ) THEN
        IF ( DO_DELTAM_SCALING ) THEN        
          MAIL   = 'Bad input: Full SS, diffuse-field delta-M on'
          ACTION = 'Full SS: default to deltam_scaling = false'
          STATUS = LIDORT_WARNING
          DO_DELTAM_SCALING   = .FALSE.
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
      ENDIF

C  Check DB correction only for non-Lambertian surface
C    ---WARNING. Flag turned off

c      IF ( DO_DBCORRECTION ) THEN
c        IF ( DO_LAMBERTIAN_SURFACE ) THEN
c          MAIL   = 'Bad input: No DB correction for Lambertian only'
c          ACTION = 'Warning: DO_DBCORRECTION turned off internally'
c          STATUS = LIDORT_WARNING
c          DO_DBCORRECTION = .FALSE.
c          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
c        ENDIF
c      ENDIF

C  Check surface inputs
C  ====================

C  Bookkeeping for surface kernel Cox-Munk types

      IF ( .NOT. DO_LAMBERTIAN_SURFACE ) THEN
        DO I = 1, N_BRDF_KERNELS
          IF ( BRDF_NAMES(I) .EQ. 'Cox-Munk  ' ) THEN
            N_BRDF_PARAMETERS(I) = 3
            IF ( DO_SHADOW_EFFECT ) THEN
              BRDF_PARAMETERS(I,3) = ONE
            ELSE
              BRDF_PARAMETERS(I,3) = ZERO
            ENDIF
          ENDIF
       ENDDO
      ENDIF

C  Check list of names (non-Lambertian)

      IF ( .NOT. DO_LAMBERTIAN_SURFACE ) THEN
       DO K = 1, N_BRDF_KERNELS
        IF ( WHICH_BRDF(K).GT.MAXBRDF_IDX.OR.WHICH_BRDF(K).LE.0) THEN
          MAIL   = 'Bad input: BRDF Index not on list of indices'
          ACTION = 'look in LIDORT.PARS for correct index'
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ELSE
          IF ( BRDF_NAMES(K).NE.BRDF_CHECK_NAMES(WHICH_BRDF(K)) ) THEN
            MAIL   = 'Bad input: BRDF kernel name not corresponding'
            ACTION = 'look in LIDORT.PARS for correct name'
            STATUS = LIDORT_SERIOUS
            CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
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
         STATUS = LIDORT_WARNING
         DO_THERMAL_TRANSONLY = .FALSE.
         CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
       ENDIF
      ENDIF

C  Check use of Green's function with thermal
C    Set default, regardless of whether solar sources are on.
C   Only require if using scattering with thermal input (7/31/07)

      IF ( DO_THERMAL_EMISSION ) THEN
       IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
        IF ( DO_CLASSICAL_SOLUTION ) THEN
         MAIL = 'Bad input: Must use Greens function for thermal'
         ACTION = 'Warning: Greens function turned on internally'
         STATUS = LIDORT_WARNING
         DO_CLASSICAL_SOLUTION = .FALSE.
         CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
         ENDIF
       ENDIF
      ENDIF

C  Switch off a bunch of flags

      IF ( DO_THERMAL_EMISSION ) THEN
       IF ( DO_THERMAL_TRANSONLY ) THEN
         DO_RAYLEIGH_ONLY  = .FALSE.
         DO_ISOTROPIC_ONLY = .FALSE.
         DO_DELTAM_SCALING = .FALSE.
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
            STATUS = LIDORT_WARNING
            EARTH_RADIUS = 6371.0D0
            CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          ENDIF
        ENDIF
      ENDIF

C  Check dimensioning on Legendre numbers (refractive geometry only)

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN
        NALLSTREAMS = NBEAMS*NLAYERS + NSTREAMS + N_USER_STREAMS
        IF ( NALLSTREAMS .GT. MAX_ALLSTRMS_P1 ) THEN
          MAIL   = 'Dimensioning error for refractive beam angles'
          ACTION = 'Increase dimension MAX_ALLSTRMS_P1'
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
      ENDIF

C  Check solar zenith angle input

      DO I = 1, NBEAMS
        IF ( BEAM_SZAS(I) .LT. ZERO .OR.
     &       BEAM_SZAS(I).GE.90.0D0 ) THEN
          WRITE(C2,'(I2)')I
          MAIL   = 'Bad input: out-of-range beam angle, no. '//C2
          ACTION = 'Look at BEAM_SZAS input, should be < 90 & > 0'
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
      ENDDO

C  Check zenith tolerance input
C    ---WARNING. Default of 0.001 will be set

      IF ( ZENITH_TOLERANCE.LE.ZERO .OR.
     &     ZENITH_TOLERANCE.GT.0.001 ) THEN
        MAIL   = 'Bad input: Zenith tolerance level out of bounds'
        ACTION = 'Warning: ZENITH_TOLERANCE set to 0.001 internally'
        STATUS = LIDORT_WARNING
        ZENITH_TOLERANCE = 0.001D0
        CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
      ENDIF

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
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
      ENDDO

C  limits on user-defined options

      IF ( .NOT.DO_USER_STREAMS .AND..NOT.DO_QUAD_OUTPUT ) THEN
        MAIL   = 'Bad input: No angular stream output is specified'
        ACTION = 'Check DO_USER_STREAMS and DO_QUAD_OUTPUT flags'
        STATUS = LIDORT_SERIOUS
        CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
      ENDIF

C  check user-defined stream angles (should always be [0,90])

      IF ( DO_USER_STREAMS ) THEN
        LOOP = .TRUE.
        I = 0
        DO WHILE (LOOP .AND. I.LT.N_USER_STREAMS)
          I = I + 1
          IF ( USER_ANGLES_INPUT(I) .GT. 90.0   .OR.
     &         USER_ANGLES_INPUT(I) .LT. ZERO ) THEN
            WRITE(C2,'(I2)')I
            MAIL   = 'Bad input: out-of-range user stream, no. '//C2
            ACTION = 'Look at user-defined angle input'
            LOOP = .FALSE.
            STATUS = LIDORT_SERIOUS
            CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          ENDIF
        ENDDO
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
            STATUS = LIDORT_SERIOUS
            CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          ENDIF
        ENDDO
      ENDIF

C  Check optical depth input
C  =========================

C  set derived Taugrid input

      TAUATM     = ZERO
      LOCALTAU(0) = TAUATM
      DO N = 1, NLAYERS
        TAUATM = TAUATM + DELTAU_VERT_INPUT(N)
        LOCALTAU(N) = TAUATM
      ENDDO

C  Check non-negative optical thickness values

      DO L = 1, NLAYERS
        IF ( DELTAU_VERT_INPUT(L).LE.ZERO ) THEN
          WRITE(C3,'(I3)')L
          MAIL = 'Bad input: optical thickness <= 0, layer '//C3
          ACTION = 'Check optical thickness input'
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
      ENDDO

C  check user-defined optical depths (should always be within atmosphere!)

      IF ( DO_USER_TAUS ) THEN
        LOOP = .TRUE.
        I = 0
        DO WHILE (LOOP .AND. I.LT.N_OUT_USERTAUS)
          I = I + 1
          IF ( USER_TAUS_INPUT(I) .GT. TAUATM )  THEN
            WRITE(C2,'(I2)')I
            MAIL   = 'Bad input: user optical depth too big, '//C2
            ACTION = 'Check optical depth values in range'
            LOOP = .FALSE.
            STATUS = LIDORT_SERIOUS
            CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          ENDIF
        ENDDO
      ELSE IF ( DO_LBOUND_TAUS ) THEN
        IF ( N_OUT_USERTAUS .GT. NLAYERS + 1 ) THEN
          MAIL   = 'Bad input: Too many Layer boundary optical depths'
          ACTION = 'Check N_OUT_USERTAUS = NLAYERS + 1'
          LOOP = .FALSE.
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
        DO I = 1, N_OUT_USERTAUS
          IF ( LBOUND_TAUS_INPUT(I) .GT. NLAYERS ) THEN
            WRITE(C2,'(I2)')LBOUND_TAUS_INPUT(I)
            MAIL   = 'Bad input: Layer boundary index > NLAYERS: '//C2
            ACTION = 'Check LBOUND_TAUS_INPUT layer boundary indices'
            LOOP = .FALSE.
            STATUS = LIDORT_SERIOUS
            CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          ENDIF
        ENDDO
      ENDIF

C  check Number of offgrid user-defined optical depths

      IF ( DO_USER_TAUS ) THEN
        NSTART = 0
        DO UTA = 1, N_OUT_USERTAUS
          TAU = USER_TAUS_INPUT(UTA)
          LOOP = .TRUE.
          N = 0
          DO WHILE (LOOP.AND.N.LT.NLAYERS)
            N = N + 1
            IF ( DABS(TAU-LOCALTAU(N-1)).LT.SMALLNUM ) THEN
              NSTART = NSTART + 1
            ENDIF
          ENDDO
        ENDDO
        N_OFFGRID_LOCAL = N_OUT_USERTAUS - NSTART
        IF ( N_OFFGRID_LOCAL .GT. MAX_OFFGRID_USERTAUS ) THEN
          WRITE(C2,'(I2)')N_OFFGRID_LOCAL
          MAIL   = 'Bad input: too many offgrid optical depths : '//C2
          ACTION =
     &        'Check number of offgrid Tau does not exceed dimensioning'
          LOOP = .FALSE.
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
      ENDIF

C  check repetition of user-defined optical depth input

      IF ( DO_USER_TAUS ) THEN
        UTA = 0
        LOOP = .TRUE.
        DO WHILE ( LOOP .AND. UTA .LT. N_OUT_USERTAUS )
          UTA = UTA + 1
          TAU = USER_TAUS_INPUT(UTA)
          NSTART = 0
          DO N = 1, N_OUT_USERTAUS
            IF ( TAU .EQ. USER_TAUS_INPUT(N)) NSTART = NSTART + 1
          ENDDO
          IF ( NSTART .NE. 1 ) THEN
            LOOP = .FALSE.
            MAIL   = 'Bad input: repetition of optical depth input'
            ACTION = 'Check user-defined optical depth input'
            STATUS = LIDORT_SERIOUS
            CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          ENDIF
        ENDDO
      ENDIF
    
C  check repetition of layer boundary optical depth input

      IF ( DO_LBOUND_TAUS ) THEN
        UTA = 0
        LOOP = .TRUE.
        DO WHILE (LOOP .AND. UTA .LT. N_OUT_USERTAUS )
          UTA = UTA + 1
          N1 = LBOUND_TAUS_INPUT(UTA)
          NSTART = 0
          DO N = 1, N_OUT_USERTAUS
            IF ( LBOUND_TAUS_INPUT(N) .EQ. N1 ) NSTART = NSTART + 1
          ENDDO
          IF ( NSTART .GT. 1 ) THEN
            LOOP = .FALSE.
            MAIL   = 'Bad input: repetition of optical depth input'
            ACTION = 'Check optical depth layer boundary input'
            STATUS = LIDORT_SERIOUS
            CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
          ENDIF
        ENDDO
      ENDIF

C  check geophysical scattering inputs
C  -----------------------------------

C  check single scatter albedos

      IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
       DO L = 1, NLAYERS
        IF ( OMEGA_TOTAL_INPUT(L).GT.ONE-OMEGA_SMALLNUM ) THEN
          WRITE(C3,'(I3)')L
          MAIL = 'Bad input: SS-albedo too close to 1, layer '//C3
          ACTION = 'Check SS-albedo input'
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
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
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
         ENDIF
        ENDDO
      ENDIF
      
C  Check first phase function moments

      IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
       DO L = 1, NLAYERS
        IF ( PHASMOMS_TOTAL_INPUT(0,L).NE.ONE ) THEN
          WRITE(C3,'(I3)')L
          MAIL   ='First phase moment not 1 for layer '//C3
          ACTION = 'Moment has been set internally to 1.0'
          STATUS = LIDORT_WARNING
          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
        ENDIF
       ENDDO
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_DERIVE_INPUT

C  Read, check and write to file of all control input

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  include files with input variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'

C  include file with bookkeeping variables (output from this module)

      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  setups

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
 
C  local variables
C  ---------------

      INTEGER          INDEX_ANGLES ( MAX_OUT_STREAMS )
      DOUBLE PRECISION ALL_ANGLES ( MAX_OUT_STREAMS ), TAU, MU1, MU2
      INTEGER          I, I1, UT, N, N1, UTA, NSTART, JINDEX, J
      LOGICAL          LOOP
      INTEGER          CONV_START_STREAMS

C  set additional numbers (derived input)
C  ======================

C  Mode of operation

      IF ( DO_FULLRAD_MODE ) THEN
        IF ( DO_SSCORR_NADIR .OR. DO_SSCORR_OUTGOING ) THEN
          DO_MSMODE_LIDORT = .TRUE.
        ELSE
          DO_MSMODE_LIDORT = .FALSE.
        ENDIF
      ELSE
        DO_MSMODE_LIDORT = .TRUE.
      ENDIF

C  SSFULL flag cancels some other flags

      IF ( DO_SSFULL ) THEN
c        DO_DBCORRECTION    = .FALSE.
        DO_DOUBLE_CONVTEST = .FALSE.
        DO_SOLUTION_SAVING = .FALSE.
        DO_BVP_TELESCOPING = .FALSE.
        DO_ADDITIONAL_MVOUT = .FALSE.
        DO_MVOUT_ONLY       = .FALSE.
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

C  make sure the Lambertian surface is defined
C  set defaults for Lambertian case (Derived input)

      IF ( DO_LAMBERTIAN_SURFACE ) THEN
        N_BRDF_KERNELS = 1
        LAMBERTIAN_KERNEL_FLAG(1) = .TRUE.
        BRDF_FACTORS(1) = LAMBERTIAN_ALBEDO
      ENDIF

C  Number of moments

      IF ( DO_RAYLEIGH_ONLY ) THEN
        NMOMENTS = 2
      ENDIF
      IF ( DO_ISOTROPIC_ONLY ) THEN
        NMOMENTS = 0
      ENDIF
      IF ( .NOT.DO_RAYLEIGH_ONLY.AND..NOT.DO_ISOTROPIC_ONLY ) THEN
        NMOMENTS = MIN ( 2 * NSTREAMS - 1, NMOMENTS_INPUT )
      ENDIF

C  total quadratures (up and down)

      NSTREAMS_2 = 2*NSTREAMS

C  Set Quadrature abscissae and weights
C    ------- 2-stream, use 1/rt(3)

      if ( nstreams .gt. 1 ) then
        CALL GAULEG(ZERO,ONE,X,A,NSTREAMS)
      else
        if ( do_full_quadrature ) then
          x(1) = 1.d0/dsqrt(3.d0)
          a(1) = 1.d0
        else
          x(1) = 0.5d0
          a(1) = 1.d0
        endif
      endif

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
        AX(I)   = X(I)*A(I)
        A5(I)   = HALF * A(I)
        XANG(I) = DACOS(X(I))/DEG_TO_RAD
        SX(I)   = DSQRT(ONE-X(I)*X(I))
      ENDDO

C  size of boundary value problem matrices and vectors

      NTOTAL = NLAYERS*2*NSTREAMS

C  number of sub and super diagonals in band matrix (boundary value problem)

      IF ( NLAYERS .EQ. 1 ) THEN
        N_SUBDIAG = 2*NSTREAMS - 1
        N_SUPDIAG = 2*NSTREAMS - 1
      ELSE
        N_SUBDIAG = 3*NSTREAMS - 1
        N_SUPDIAG = 3*NSTREAMS - 1
      ENDIF

C  solar zenith angle cosines/sines

      IF ( DO_SOLAR_SOURCES ) THEN
        DO I = 1, NBEAMS
          X0(I)  = DCOS ( BEAM_SZAS(I) * DEG_TO_RAD )
          SX0(I) = DSQRT(ONE-X0(I)*X0(I))
        ENDDO 
      ENDIF

C  Set average cosines in the refractive geometry case

      IF ( DO_SOLAR_SOURCES ) THEN
       IF ( DO_REFRACTIVE_GEOMETRY ) THEN
        DO I = 1, NBEAMS
          MU1 = DCOS(SZA_LOCAL_INPUT(0,I)*DEG_TO_RAD)
          DO N = 1, NLAYERS
            MU2 = DCOS(SZA_LOCAL_INPUT(N,I)*DEG_TO_RAD)
            SUN_SZA_COSINES(N,I) = HALF * ( MU1 + MU2 )
            MU1 = MU2
         ENDDO
        ENDDO
       ENDIF
      ENDIF

C  Set the angle masks
C  ===================

C  initialize
C  ----------

      DO I = 1, MAX_USER_STREAMS
        USEROUTPUT_INDEX(I) = 0
      ENDDO
      DO I = 1, MAXSTREAMS
        QUADOUTPUT_INDEX(I) = 0
      ENDDO

C  If quadrature output is required
C  --------------------------------

      IF ( DO_QUAD_OUTPUT ) THEN

C  .. For quadrature only, OUT_ANGLES are just the quadrature angles

        IF ( .NOT. DO_USER_STREAMS ) THEN

C  .. number of streams for convergence and output

          N_OUT_STREAMS  = NSTREAMS
          CONV_START_STREAMS = 1

          DO I = 1, N_OUT_STREAMS
            I1 = NSTREAMS + 1 - I
            QUADOUTPUT_INDEX(I) = I1
            OUT_ANGLES(I) = XANG(I1)
          ENDDO

C  For Quadrature AND user-defined stream angles

        ELSE

C  .. number of streams for output

          N_OUT_STREAMS  = NSTREAMS + N_USER_STREAMS

C   .. Set a dummy array ALL_ANGLES of all angles
C      ( order first with the user streams )

          DO I = 1, N_USER_STREAMS
            ALL_ANGLES(I) = USER_ANGLES_INPUT(I)
          ENDDO
          DO I = 1, NSTREAMS
            I1 = I + N_USER_STREAMS
            ALL_ANGLES(I1) = XANG(I)
          ENDDO

C   .. Use an indexing routine INDEXX to sort ALL_ANGLES

          CALL INDEXX(N_OUT_STREAMS, ALL_ANGLES, INDEX_ANGLES )

C   .. Set the output angles according to this indexing

          DO I = 1, N_OUT_STREAMS
            OUT_ANGLES(I) = ALL_ANGLES(INDEX_ANGLES(I))
          ENDDO

C   .. Set the masks for quadrature and user-defined streams

          DO I = 1, NSTREAMS
            I1 = I + N_USER_STREAMS
            LOOP = .TRUE.
            J = 0
            DO WHILE ( LOOP )
              J = J + 1
              JINDEX = INDEX_ANGLES(J)
              IF ( JINDEX .EQ. I1 ) LOOP = .FALSE.
            ENDDO
            QUADOUTPUT_INDEX(I) = J
          ENDDO

          DO I = 1, N_USER_STREAMS
            LOOP = .TRUE.
            J = 0
            DO WHILE ( LOOP )
              J = J + 1
              JINDEX = INDEX_ANGLES(J)
              IF ( JINDEX .EQ. I ) LOOP = .FALSE.
            ENDDO
            USEROUTPUT_INDEX(I) = J
          ENDDO

C  number of streams for convergence (uses zenith tolerance)

          J = 0
          DO I = 1, N_USER_STREAMS
            IF ( DABS(ONE-DCOS(OUT_ANGLES(I)*DEG_TO_RAD)).LT.
     &                     ZENITH_TOLERANCE ) THEN
              J = J + 1
            ENDIF
          ENDDO
          CONV_START_STREAMS = J + 1

        ENDIF

C  No Quadratures - only user-defined streams

      ELSE

        IF ( DO_USER_STREAMS ) THEN

C  ..  number of streams for output

          N_OUT_STREAMS  = N_USER_STREAMS

C   .. Index the user-defined angles and set the output anges & masks
C      ( Check already made that User-defined option is set)

          DO I = 1, N_USER_STREAMS
            ALL_ANGLES(I) = USER_ANGLES_INPUT(I)
          ENDDO

          IF ( N_OUT_STREAMS .EQ. 1 ) THEN
            OUT_ANGLES(1) = ALL_ANGLES(1)
            USEROUTPUT_INDEX(1) = 1
          ELSE
            CALL INDEXX
     &          ( N_OUT_STREAMS, USER_ANGLES_INPUT, INDEX_ANGLES )
            DO I = 1, N_OUT_STREAMS
              OUT_ANGLES(I) = ALL_ANGLES(INDEX_ANGLES(I))
              LOOP = .TRUE.
              J = 0
              DO WHILE ( LOOP )
                J = J + 1
                JINDEX = INDEX_ANGLES(J)
                IF ( JINDEX .EQ. I ) LOOP = .FALSE.
              ENDDO
              USEROUTPUT_INDEX(I) = J
            ENDDO
          ENDIF

C  number of streams for convergence (uses zenith tolerance)

          J = 0
          DO I = 1, N_USER_STREAMS
            IF ( DABS(ONE-DCOS(OUT_ANGLES(I)*DEG_TO_RAD)).LT.
     &                     ZENITH_TOLERANCE ) THEN
              J = J + 1
            ENDIF
          ENDDO
          CONV_START_STREAMS = J + 1

        ENDIF

      ENDIF

C  User stream cosines and secants

      IF ( DO_USER_STREAMS ) THEN
        DO I = 1, N_USER_STREAMS
          USER_STREAMS(I) = DCOS(DEG_TO_RAD*USER_ANGLES_INPUT(I))
          USER_SECANTS(I) = ONE / USER_STREAMS(I)
          USER_SINES(I) = DSQRT(ONE-USER_STREAMS(I)*USER_STREAMS(I))
        ENDDO
      ENDIF

C  number of tests to be applied for convergence

C  tempo

      CONV_START_STREAMS = 1

      N_CONV_STREAMS = N_OUT_STREAMS - CONV_START_STREAMS + 1
      N_CONVTESTS = N_USER_RELAZMS * N_CONV_STREAMS * N_DIRECTIONS
      N_CONVTESTS = N_CONVTESTS * N_OUT_USERTAUS 

C  Sort out User optical depths
C  ----------------------------

C  set derived Taugrid input

      TAUGRID_INPUT(0) = ZERO
      DO N = 1, NLAYERS
        TAUGRID_INPUT(N) = TAUGRID_INPUT(N-1)+ DELTAU_VERT_INPUT(N)
      ENDDO

C  ( These should be checked beforehand to lie within range of TAUGRID )

      IF ( DO_USER_TAUS ) THEN

C  Sort in ascending order

        IF ( N_OUT_USERTAUS .GT. 1 ) THEN
          CALL HPSORT(N_OUT_USERTAUS,USER_TAUS_INPUT)
        ENDIF

C  mark all optical depths not equal to layer boundary values

        NSTART = 0
        UT = 0
        DO UTA = 1, N_OUT_USERTAUS
          TAU = USER_TAUS_INPUT(UTA)
          OFFGRID_UTAU_OUTFLAG(UTA)  = .FALSE.
          OFFGRID_UTAU_OUTINDEX(UTA) =   0
          LOOP = .TRUE.
          N = NSTART
          DO WHILE (LOOP.AND.N.LT.NLAYERS)
            N = N + 1
            IF ( (TAU-TAUGRID_INPUT(N-1)).GT.SMALLNUM .AND.
     &             (TAUGRID_INPUT(N)-TAU).GT.SMALLNUM ) THEN
              UT = UT + 1
              OFFGRID_UTAU_OUTINDEX(UTA) =   UT
              OFFGRID_UTAU_OUTFLAG(UTA)  = .TRUE.
              OFFGRID_UTAU_LAYERIDX(UT) = N
              UTAU_LEVEL_MASK_UP(UTA) = N
              UTAU_LEVEL_MASK_DN(UTA) = N - 1
              LOOP = .FALSE.
              NSTART = N - 1
            ENDIF
          ENDDO
        ENDDO
        N_OFFGRID_USERTAUS = UT

C  set remaining optical depth mask (for values at layer boundaries)

        NSTART = 0
        DO UTA = 1, N_OUT_USERTAUS
          IF ( .NOT. OFFGRID_UTAU_OUTFLAG(UTA) ) THEN
            TAU = USER_TAUS_INPUT(UTA)
            LOOP = .TRUE.
            N = NSTART
            DO WHILE (LOOP.AND.N.LE.NLAYERS)
              N1 = N
              N = N + 1
              IF ( DABS(TAU-TAUGRID_INPUT(N1)).LT.SMALLNUM ) THEN
                LOOP = .FALSE.
                UTAU_LEVEL_MASK_UP(UTA) = N1
                UTAU_LEVEL_MASK_DN(UTA) = N1
                NSTART = N1
              ENDIF
            ENDDO
          ENDIF
        ENDDO

C  Otherwise if the optical depths are set to be at all layer boundaries

      ELSE IF ( DO_LBOUND_TAUS ) THEN

C  assign

        N_OFFGRID_USERTAUS = 0
        DO UTA = 1, N_OUT_USERTAUS
          USER_TAUS_INPUT(UTA) = TAUGRID_INPUT(LBOUND_TAUS_INPUT(UTA))
          OFFGRID_UTAU_OUTFLAG(UTA) = .FALSE.
        ENDDO

C  Sort in ascending order

        IF ( N_OUT_USERTAUS .GT. 1 ) THEN
          CALL HPSORT(N_OUT_USERTAUS,USER_TAUS_INPUT)
        ENDIF

C  set optical depth mask (for values at layer boundaries)

        NSTART = 0
        DO UTA = 1, N_OUT_USERTAUS
          TAU = USER_TAUS_INPUT(UTA)
          LOOP = .TRUE.
          N = NSTART
          DO WHILE (LOOP.AND.N.LE.NLAYERS)
            N1 = N
            N = N + 1
            IF ( DABS(TAU-TAUGRID_INPUT(N1)).LT.SMALLNUM ) THEN
              LOOP = .FALSE.
              UTAU_LEVEL_MASK_UP(UTA) = N1
              UTAU_LEVEL_MASK_DN(UTA) = N1
              NSTART = N1
            ENDIF
          ENDDO
        ENDDO

      ENDIF

C  Set masking and number of layer source terms
C  --------------------------------------------

C   .. for upwelling

      IF ( DO_UPWELLING ) THEN
        DO N = 1, NLAYERS
          STERM_LAYERMASK_UP(N) = .FALSE.
        ENDDO
        UTA = 1
        UT  = 1
        IF ( .NOT. OFFGRID_UTAU_OUTFLAG(UTA) ) THEN
          N_LAYERSOURCE_UP = UTAU_LEVEL_MASK_UP(UTA) + 1
          N_ALLLAYERS_UP   = N_LAYERSOURCE_UP
        ELSE
          N_LAYERSOURCE_UP = OFFGRID_UTAU_LAYERIDX(UT) + 1
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
        UTA = N_OUT_USERTAUS
        UT  = N_OFFGRID_USERTAUS
        IF ( .NOT. OFFGRID_UTAU_OUTFLAG(UTA) ) THEN
          N_LAYERSOURCE_DN = UTAU_LEVEL_MASK_DN(UTA)
          N_ALLLAYERS_DN   = N_LAYERSOURCE_DN
        ELSE
          N_LAYERSOURCE_DN = OFFGRID_UTAU_LAYERIDX(UT)
          N_ALLLAYERS_DN   = N_LAYERSOURCE_DN
        ENDIF
        DO N = 1, N_ALLLAYERS_DN
          STERM_LAYERMASK_DN(N) = .TRUE.
        ENDDO
      ENDIF

C  Finish

      RETURN
      END
