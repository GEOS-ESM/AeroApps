C ###########################################################
C #                                                         #
C #                    THE LIDORT FAMILY                    #
C #                                                         #
C #      (LInearized Discrete Ordinate Radiative Transfer)  #
C #       --         -        -        -         -          #
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
C #            LIDORT_L_INPUT_MASTER (master)                   #
C #            LIDORT_L_READ_INPUT                              #
C #            LIDORT_L_ZERO_INPUT                              #
C #            LIDORT_L_CHECK_INPUT                             #
C #                                                             #
C ###############################################################

      SUBROUTINE LIDORT_L_INPUT_MASTER ( FILNAM, EFILNAM, STATUS )

C  Read all control inputs for LIDORT
C  ----------------------------------

C  include file of constants and dimensions

      INCLUDE '../includes/LIDORT.PARS'

C  include file of input variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'

C  Module arguments (input filename, error filename,output status)

      CHARACTER*(*)      FILNAM
      CHARACTER*(*)      EFILNAM
      INTEGER            STATUS

C  local variables

      INTEGER        STATUS_SUB, FILUNIT, LEN_STRING
      CHARACTER*70   MAIL
      EXTERNAL       LEN_STRING

C  initialize status

      STATUS = LIDORT_SUCCESS

C  initialize error file

      LIDORT_ERROR_FILENAME = EFILNAM
      LIDORT_ERROR_INIT     = .TRUE.

C  initialize standard variables
 
      CALL LIDORT_INIT_CONTROL_VARS
      CALL LIDORT_INIT_SURFACE_VARS
      CALL LIDORT_INIT_MODEL_VARS

C  initialize linearization variables

      CALL LIDORT_INIT_L_INPUTS_VARS

C  Open file

      FILUNIT = LIDORT_INUNIT
      OPEN(LIDORT_INUNIT,FILE=FILNAM,ERR=300,STATUS='OLD')

C  read standard inputs

      CALL LIDORT_READ_INPUT ( STATUS_SUB )
      IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
        STATUS = LIDORT_SERIOUS
        MAIL = 'Error from LIDORT_READINPUT'
        CALL LIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
        CLOSE(FILUNIT)
        RETURN
      ENDIF

C  read additional inputs for linearization

      CALL LIDORT_L_READ_INPUT ( STATUS_SUB )
      IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
        STATUS = LIDORT_SERIOUS
        MAIL = 'Error from LIDORT_L_INPUTREAD'
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

C

      SUBROUTINE LIDORT_L_READ_INPUT ( STATUS )

C  Read all linearization control inputs for LIDORT
C  ------------------------------------------------

C  include file of constants and dimensions

      INCLUDE '../includes/LIDORT.PARS'

C  include file of regular inputs

      INCLUDE '../includes/LIDORT_INPUTS.VARS'

C  include file of linearization inputs to be filled out

      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'

C  output status

      INTEGER         STATUS

C  local variables

      CHARACTER*8     PREFIX
      PARAMETER       ( PREFIX = 'LIDORT -' )
      LOGICAL         ERROR
      CHARACTER * 80  PAR_STR
      LOGICAL         GFINDPAR
      EXTERNAL        GFINDPAR
      INTEGER         FILUNIT, LEN_STRING, I, J
      INTEGER         DUM_INDEX, DUM_NPARS, NC
      CHARACTER*10    DUM_NAME
      CHARACTER*70    MAIL
      EXTERNAL        LEN_STRING

C  initialize status

      STATUS = LIDORT_SUCCESS
      ERROR  = .FALSE.
      NC = 0

C  file unit

      FILUNIT = LIDORT_INUNIT

C  linearization control

      PAR_STR = 'Simulation only?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_SIMULATION_ONLY
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      PAR_STR = 'Atmospheric profile weighting functions?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_PROFILE_LINEARIZATION
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      PAR_STR = 'Atmospheric column weighting functions?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_COLUMN_LINEARIZATION
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      PAR_STR = 'Surface property weighting functions?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_SURFACE_LINEARIZATION
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

c      PAR_STR = 'Surface emission weighting functions?'
c      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
c     &     READ (FILUNIT,*,ERR=998) DO_SURFBB_LINEARIZATION
c      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  Atmospheric profile weighting function input control
C  individual layer number does not exceed whole number

      IF ( DO_PROFILE_LINEARIZATION ) THEN

        PAR_STR='Atmospheric profile weighting functions, total number'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR))
     &    READ (FILUNIT,*,ERR=998) N_TOTALPROFILE_WFS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

c        PAR_STR = 'Atmospheric weighting functions, layer control'
c        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
c          DO I = 1, NLAYERS
c            READ (FILUNIT,*,ERR=998) LAYER_VARY_FLAG(I),
c     &                               LAYER_VARY_NUMBER(I)
c            IF (LAYER_VARY_NUMBER(I) .GT. N_TOTALATMOS_WFS ) THEN
c              LAYER_VARY_NUMBER(I) = N_TOTALATMOS_WFS
c            ENDIF
c          ENDDO
c        ENDIF
c        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      ENDIF

C  atmospheric Bulk/column weighting function input control

       IF ( DO_COLUMN_LINEARIZATION ) THEN

         PAR_STR='Atmospheric column weighting functions, total number'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR))
     &    READ (FILUNIT,*,ERR=998) N_TOTALCOLUMN_WFS
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

         DO I = 1, NLAYERS
           LAYER_VARY_FLAG(I)   = .TRUE.
           LAYER_VARY_NUMBER(I) = N_TOTALCOLUMN_WFS
         ENDDO

       ENDIF

C  Surface weighting function input control

      IF ( DO_SURFACE_LINEARIZATION ) THEN

C  Lambertian WF flag

       IF ( DO_LAMBERTIAN_SURFACE ) THEN

        PAR_STR = 'Do Lambertian albedo weighting function?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          READ (FILUNIT,*,ERR=998)DO_KERNEL_FACTOR_WFS(1)
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )
        N_TOTALBRDF_WFS = 1

C  kernel input. Names should agree with earlier entries.

       ELSE
        PAR_STR = 'Kernels, indices, # pars, Jacobian flags'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_BRDF_KERNELS
            READ (FILUNIT,57,ERR=998)
     &       DUM_NAME, DUM_INDEX,DUM_NPARS,DO_KERNEL_FACTOR_WFS(I),
     &           (DO_KERNEL_PARAMS_WFS(I,J),J=1,3)

            IF ( DUM_NAME .NE. BRDF_NAMES(I) ) THEN
              MAIL = 'Input BRDF name not same as earlier list'
              STATUS  = LIDORT_SERIOUS
              CALL LIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
              RETURN
            ENDIF

            IF ( DUM_INDEX .NE. WHICH_BRDF(I) ) THEN
              MAIL = 'Input BRDF index not same as earlier list'
              STATUS  = LIDORT_SERIOUS
              CALL LIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
              RETURN
            ENDIF

            IF ( DUM_NPARS .NE. N_BRDF_PARAMETERS(I) ) THEN
              MAIL = 'Input BRDF parameter no. not same as earlier'
              STATUS  = LIDORT_SERIOUS
              CALL LIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
              RETURN
            ENDIF

C  Compute total number of pars

            IF ( DO_KERNEL_FACTOR_WFS(I) ) NC = NC + 1
            DO J = 1, N_BRDF_PARAMETERS(I)
              IF ( DO_KERNEL_PARAMS_WFS(I,J) ) NC = NC + 1
            ENDDO

          ENDDO
          N_TOTALBRDF_WFS = NC
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )
 57     FORMAT( A10, I3, I2, 1X, L2, 2X, 3L2 ) 

C  Check total number of BRDF weighting functions is not out of bounds

        IF ( N_TOTALBRDF_WFS .GT. MAX_SURFACEWFS ) THEN
         MAIL  = 'Bad error: Insuffient dimensioning for Surface WFs'
         STATUS = LIDORT_SERIOUS
         CALL LIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
         RETURN
        ENDIF
   
       ENDIF
      ENDIF

C  weighting function names

      IF ( DO_PROFILE_LINEARIZATION ) THEN
        PAR_STR='Atmospheric profile Jacobian names (character*31)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_TOTALPROFILE_WFS
            READ (FILUNIT,'(a31)',ERR=998) PROFILEWF_NAMES(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )
      ENDIF

      IF ( DO_COLUMN_LINEARIZATION ) THEN
        PAR_STR='Atmospheric column Jacobian names (character*31)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_TOTALCOLUMN_WFS
            READ (FILUNIT,'(a31)',ERR=998) COLUMNWF_NAMES(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )
      ENDIF

      IF ( DO_SURFACE_LINEARIZATION ) THEN
       IF ( .NOT. DO_LAMBERTIAN_SURFACE ) THEN
        PAR_STR = 'Surface Jacobian names (character*22)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_TOTALBRDF_WFS
            READ (FILUNIT,'(a22)',ERR=998) SURFACEWF_NAMES(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )
       ELSE
        SURFACEWF_NAMES(1) = 'Lambertian Amplitude  '
       ENDIF
      ENDIF

C  set overall linearization flags

      DO_ATMOS_LINEARIZATION = ( DO_PROFILE_LINEARIZATION.OR.
     &                           DO_COLUMN_LINEARIZATION )

      DO_LINEARIZATION = ( DO_ATMOS_LINEARIZATION   .OR.
     &                     DO_SURFACE_LINEARIZATION .OR.
     &                     DO_SURFBB_LINEARIZATION )

C  normal return

      RETURN

C  line read error - abort immediately

998   CONTINUE
      STATUS = LIDORT_SERIOUS
      MAIL = 'read failure for '//PAR_STR(1:LEN_STRING(PAR_STR))
      CALL LIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
      RETURN

C  Finish

      END

C

      SUBROUTINE LIDORT_INIT_L_INPUTS_VARS

C  Initialises all linearization control inputs for LIDORT
C  -------------------------------------------------------

C  include file of constants and dimensions

      INCLUDE '../includes/LIDORT.PARS'

C  include file to be initialised

      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'

C  local variables

      INTEGER             B, Q, N

C  initialise overall control inputs

      DO_SIMULATION_ONLY        = .FALSE.
      DO_LINEARIZATION          = .FALSE.

C  initialise atmospheric weighting function control inputs

      DO_PROFILE_LINEARIZATION  = .FALSE.
      DO_COLUMN_LINEARIZATION   = .FALSE.
      DO_ATMOS_LINEARIZATION    = .FALSE.
      N_TOTALCOLUMN_WFS  = 0
      N_TOTALPROFILE_WFS = 0
      DO N = 1, MAXLAYERS
        LAYER_VARY_NUMBER(N) = 0
        LAYER_VARY_FLAG(N)   = .FALSE.
      ENDDO

C  initialise surface weighting function control inputs

      DO_SURFACE_LINEARIZATION = .FALSE.
      DO_SURFBB_LINEARIZATION  = .FALSE.
      N_TOTALBRDF_WFS  = 0
      DO B = 1, MAX_BRDF_KERNELS
        DO_KERNEL_FACTOR_WFS(B) = .FALSE.
        DO Q = 1, MAX_BRDF_PARAMETERS
          DO_KERNEL_PARAMS_WFS(B,Q) = .FALSE.
        ENDDO
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_L_CHECK_INPUT (STATUS)

C  include files
C  -------------

C  include file of constants and dimensions

      INCLUDE '../includes/LIDORT.PARS'

C  include file of regular inputs

      INCLUDE '../includes/LIDORT_INPUTS.VARS'

C  include file of linearization inputs to be checked out

      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'

C  Module output
C  -------------

      INTEGER           STATUS

C  local variables
C  ---------------

      CHARACTER*70      MAIL, TRACE
      CHARACTER*1       NOTRACE
      PARAMETER         ( NOTRACE = ' ')
      
C  Initialize output status

      STATUS = LIDORT_SUCCESS

C  Do some checking
C  ================

C  Check formatting; should not exceed 10 entries for the number of
C  layer weighting functions. This is needed to avoid errors with
C  compilers not using variable formatting.

C      IF ( MAX_PARAMETERS .GT. 10 ) THEN
C        MAIL='Max layer variations > 10, get formatting problems!'
C        STATUS = LIDORT_SERIOUS
C        CALL LIDORT_ERROR_TRACE ( MAIL, NOTRACE, STATUS )
C      ENDIF

C  Check something is being varied

      IF ( .NOT. DO_SIMULATION_ONLY ) THEN
        IF ( .NOT. DO_ATMOS_LINEARIZATION .AND.
     &         .NOT. DO_SURFACE_LINEARIZATION. AND.
     &         .NOT. DO_SURFBB_LINEARIZATION  ) THEN
          MAIL='Bad input: Nothing being varied or calculated'
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, NOTRACE, STATUS )
        ENDIF
      ENDIF

C  Check and fix AF flags if you want simulation only

      IF ( DO_SIMULATION_ONLY ) THEN
        IF ( DO_ATMOS_LINEARIZATION ) THEN
          MAIL  = 'Input warning; simulation only ---> NO Atmos WFs'
          TRACE = 'Action: LIDORT switched off WF flag for you!'
          STATUS = LIDORT_WARNING
          CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          DO_PROFILE_LINEARIZATION = .FALSE.
          DO_COLUMN_LINEARIZATION  = .FALSE.
          DO_ATMOS_LINEARIZATION   = .FALSE.
        ENDIF
        IF ( DO_SURFACE_LINEARIZATION ) THEN
          MAIL  = 'Input warning; simulation only ---> NO Surface WFs'
          TRACE = 'Action: LIDORT switched off WF flag for you!'
          STATUS = LIDORT_WARNING
          CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          DO_SURFACE_LINEARIZATION = .FALSE.
        ENDIF
      ENDIF

C  Check total number of BRDF weighting functions is not out of bounds

      IF ( DO_SURFACE_LINEARIZATION ) THEN
       IF ( .NOT.DO_LAMBERTIAN_SURFACE ) THEN
        IF ( N_TOTALBRDF_WFS .GT. MAX_SURFACEWFS ) THEN
         MAIL  = 'Bad error: Insuffient dimensioning for Surface WFs'
         TRACE = 'Action: increase MAX_SURFACE_WFS in LIDORT.PARS!'
         STATUS = LIDORT_SERIOUS
         CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
        ENDIF
       ENDIF
      ENDIF

C  Finish

      RETURN
      END
