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
C #            LIDORT_L_FOURIER_MASTER (master)                 #
C #                                                             #
C ###############################################################

      SUBROUTINE LIDORT_L_FOURIER_MASTER
     I       ( FOURIER_COMPONENT,
     I         SS_FLUX_MULTIPLIER,
     O         DO_INCLUDE_SURFACE,
     O         DO_INCLUDE_SURFEMISS,
     O         DO_INCLUDE_THERMEMISS,
     O         STATUS )

C  Complete Fourier component calculation for the Extended Code
C    Radiance and Weighting function computations

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file of Linearization input variables

      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'

C  Include files of setup variables (output from SETUPS routine)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_THERMALSUP.VARS'

C  include this file for debug only

      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  module arguments
C  ----------------

C  Input Fourier component number

      INTEGER          FOURIER_COMPONENT

C  single scatter correction flux multiplier

      DOUBLE PRECISION SS_FLUX_MULTIPLIER

C  Output

      LOGICAL          DO_INCLUDE_SURFACE
      LOGICAL          DO_INCLUDE_SURFEMISS
      LOGICAL          DO_INCLUDE_THERMEMISS
      INTEGER          STATUS

C  Local variables
C  ---------------

      INTEGER          LAYER, IBEAM, I, IPARTIC

C  local inclusion flags

      LOGICAL          DO_INCLUDE_MVOUTPUT
      LOGICAL          DO_INCLUDE_DIRECTBEAM

C  Flux multiplier and Fourier component numbers

      DOUBLE PRECISION FLUX_MULTIPLIER
      DOUBLE PRECISION DELTA_FACTOR
      DOUBLE PRECISION SURFACE_FACTOR

C  weighting function indices

      INTEGER          N_LAYER_WFS, LAYER_TO_VARY, N_PARAMETERS
      INTEGER          BRDF_INDEX, PAR_INDEX, WF_INDEX
      LOGICAL          DO_RTSOL_VARY
      INTEGER          NPARAMS_VARY, VARIATION_INDEX

C  error tracing variables

      CHARACTER*(2)    C2
      CHARACTER*(70)   MAIL, TRACE
      INTEGER          STATUS_SUB
      LOGICAL          FAIL

C  progress

      logical          do_write_screen
      parameter       ( do_write_screen = .false. )

C  ##############
C  initialization
C  ##############

C  module status

      STATUS = LIDORT_SUCCESS

C  Set local flags
C  ---------------

C  inclusion of thermal surface emission term, only for Fourier = 0

      DO_INCLUDE_SURFEMISS = .FALSE.
      IF ( DO_SURFACE_EMISSION ) THEN
        IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
          DO_INCLUDE_SURFEMISS = .TRUE.
        ENDIF
      ENDIF

C  inclusion of thermal emission term, only for Fourier = 0

      DO_INCLUDE_THERMEMISS = .FALSE.
      IF ( DO_THERMAL_EMISSION ) THEN
        IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
          DO_INCLUDE_THERMEMISS = .TRUE.
        ENDIF
      ENDIF

C  Surface flag (for inclusion of some kind of reflecting boundary)

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

C  Direct beam flag (only if above albedo flag has been set)

      IF ( DO_SOLAR_SOURCES ) THEN
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

C  surface reflectance factors

      IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
        SURFACE_FACTOR = TWO
        DELTA_FACTOR   = ONE
      ELSE
        SURFACE_FACTOR = ONE
        DELTA_FACTOR   = TWO
      ENDIF

C  Flux multipliers
C   = 1 / 4.pi with beam sources,  = 1 for Thermal alone.

      FLUX_MULTIPLIER   = DELTA_FACTOR
      IF ( DO_INCLUDE_THERMEMISS.AND..NOT.DO_SOLAR_SOURCES) THEN
        FLUX_MULTIPLIER = DELTA_FACTOR
      ENDIF

C  inclusion of mean value output

      DO_INCLUDE_MVOUTPUT = .FALSE.
      IF ( DO_ADDITIONAL_MVOUT. OR. DO_MVOUT_ONLY ) THEN
        IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
          DO_INCLUDE_MVOUTPUT = .TRUE.
        ENDIF
      ENDIF

C  Initialise BVP telescoping. Important.

      IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
        DO_BVTEL_INITIAL = DO_BVP_TELESCOPING
      ENDIF

C  #################
C  Set up operations 
C  #################

C  Miscellaneous setups for Fourier = 0

C   MISCSETUPS    : Delta-M, average-secant formulation, transmittances
C   THERMAL_SETUP : Coefficients, direct multipliers
C   EMULT_MASTER  : Beam source function multipliers. Not required for the
C                  Full SS calculation in outgoing mode

      IF ( FOURIER_COMPONENT .EQ. 0 ) THEN

C  With linearization. Thermal setup + linearization in one routine.

        IF ( DO_ATMOS_LINEARIZATION ) THEN
          CALL LIDORT_MISCSETUPS
          CALL LIDORT_L_MISCSETUPS
          IF ( DO_INCLUDE_THERMEMISS ) CALL THERMAL_SETUP_PLUS
          IF ( DO_SOLAR_SOURCES ) THEN 
            IF (.NOT.DO_SSFULL.OR.(DO_SSFULL.AND.DO_SSCORR_NADIR)) THEN
              CALL EMULT_MASTER
              CALL L_EMULT_MASTER
            ENDIF
          ENDIF
        ENDIF

C  No linearization

        IF ( .NOT. DO_ATMOS_LINEARIZATION ) THEN
          CALL LIDORT_MISCSETUPS
          IF ( DO_INCLUDE_THERMEMISS ) CALL THERMAL_SETUP
          IF ( DO_SOLAR_SOURCES ) THEN 
            IF (.NOT.DO_SSFULL.OR.(DO_SSFULL.AND.DO_SSCORR_NADIR)) THEN
              CALL EMULT_MASTER
            ENDIF
          ENDIF
        ENDIF

      ENDIF

C  #####################
C  Correction operations 
C  #####################

C  Not required if no solar sources

      IF ( DO_SOLAR_SOURCES ) THEN

C  Single scatter correction (pre-calculation)
C   Must be done after MISCSETUPS, as we need the multipliers and
C    transmittance factors for the SUN and LOS paths.
C   Standard Code added 6 May 2005. Replaces call in Master routine.
C   Linearized code 23 May 2005. Replaces former code in L_Master module.

C      Version 3.2. Added call to the new outgoing sphericity correction
C      Version 3.3. Added call bulk/column linearizations

       IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
        IF ( DO_USER_STREAMS ) THEN

C  Regular (nadir-view) ss correction

          IF ( DO_SSCORR_NADIR ) THEN
            CALL LIDORT_SSCORR_NADIR ( SS_FLUX_MULTIPLIER )
            IF ( DO_PROFILE_LINEARIZATION ) THEN
              DO LAYER_TO_VARY = 1, NLAYERS
                IF ( LAYER_VARY_FLAG(LAYER_TO_VARY) ) THEN
                  N_PARAMETERS = LAYER_VARY_NUMBER(LAYER_TO_VARY)
                  CALL LIDORT_LP_SSCORR_NADIR
     &            ( SS_FLUX_MULTIPLIER, LAYER_TO_VARY, N_PARAMETERS)
                ENDIF
              ENDDO
            ENDIF
            IF ( DO_COLUMN_LINEARIZATION ) THEN
              CALL LIDORT_LC_SSCORR_NADIR
     &            ( SS_FLUX_MULTIPLIER, N_TOTALCOLUMN_WFS)
            ENDIF
          ENDIF

C  Outgoing sphericity correction ------ Added 20 March 2007.
C     Calling not changed for column Jacobians (23 May 2007)

          IF ( DO_SSCORR_OUTGOING ) THEN
            CALL LIDORT_L_SSCORR_OUTGOING
     &     ( SS_FLUX_MULTIPLIER, FAIL, MAIL )
            IF ( FAIL ) THEN
              TRACE= 'Linearized SS correction outgoing failed'
              STATUS = LIDORT_SERIOUS
              CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
              RETURN
            ENDIF
          ENDIF

C  End fourier = 0 and user stream clauses

        ENDIF
       ENDIF
      
C  Exact direct beam BRDF (if flagged)
C   Must be done after MISCSETUPS, as we need the SUN transmittance factors
C   BRDF Kernels already known (called before this module)
C   Code added 6 May 2005. Standard case
C   Profile/surface Linearization code added 23 May 2005.
C     Column        Linearization code added 14 May 2007.

       IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
        IF ( DO_USER_STREAMS ) THEN
          IF ( DO_DBCORRECTION ) THEN
            CALL LIDORT_DBCORRECTION ( SS_FLUX_MULTIPLIER )
            IF ( DO_PROFILE_LINEARIZATION ) THEN
              CALL LIDORT_LAP_DBCORRECTION ( SS_FLUX_MULTIPLIER )
            ENDIF
            IF ( DO_COLUMN_LINEARIZATION ) THEN
              CALL LIDORT_LAC_DBCORRECTION ( SS_FLUX_MULTIPLIER )
            ENDIF
            IF ( DO_SURFACE_LINEARIZATION ) THEN
              CALL LIDORT_LS_DBCORRECTION ( SS_FLUX_MULTIPLIER )
            ENDIF
          ENDIF
        ENDIF
       ENDIF

C  End correction loop (DO_SOLAR_SOURCES)

      ENDIF

C  BRDF Fourier components - coefficient setup
C   (also computes emissivity if required)
C   (DO_REFLECTED_DIRECTBEAM is a stored variable)
C   Should always call this, even with Lambertian surfaces!

      IF (.NOT. DO_SSFULL ) THEN
        CALL LIDORT_BRDF_FOURIER
     I    ( DO_INCLUDE_SURFACE,
     I      DO_INCLUDE_SURFEMISS,
     I      FOURIER_COMPONENT,
     I      DELTA_FACTOR )
      ENDIF

C  Linearization of BRDF Fourier components

      IF ( DO_SURFACE_LINEARIZATION ) THEN
       IF (.NOT. DO_SSFULL ) THEN
        CALL LIDORT_L_BRDF_FOURIER
     I    ( DO_INCLUDE_SURFACE,
     I      DO_INCLUDE_SURFEMISS,
     I      FOURIER_COMPONENT,
     I      DELTA_FACTOR )
       ENDIF
      ENDIF

C  Reflected Direct beam attenuation (BRDF). Diffuse field
C  Full single scatter calculation, return after completion.
C    Must add Lambertian DB term if flagged (BRDF DB term is done above)

      IF ( DO_SSFULL. AND. DO_SOLAR_SOURCES ) THEN
        IF ( DO_LAMBERTIAN_SURFACE ) THEN
          CALL LIDORT_LAMBERTIAN_DBCORRECTION (SS_FLUX_MULTIPLIER)
          IF ( DO_PROFILE_LINEARIZATION ) THEN
            CALL LIDORT_LAP_LAMBERTIAN_DBCORRECTION (SS_FLUX_MULTIPLIER)
          ENDIF
          IF ( DO_COLUMN_LINEARIZATION ) THEN
            CALL LIDORT_LAC_LAMBERTIAN_DBCORRECTION (SS_FLUX_MULTIPLIER)
          ENDIF
          IF ( DO_SURFACE_LINEARIZATION ) THEN
            CALL LIDORT_LS_LAMBERTIAN_DBCORRECTION
          ENDIF
        ENDIF
        RETURN
      ELSE
        CALL LIDORT_DIRECTBEAM
     I    ( FOURIER_COMPONENT, DELTA_FACTOR )
      ENDIF

C  Get Legendre polynomials for this Fourier component
C    Including all beam angles !!

      IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
       CALL LIDORT_LEGENDRE_SETUP ( FOURIER_COMPONENT )
       IF ( DO_USER_STREAMS ) THEN
         CALL LIDORT_USERLEGENDRE_SETUP ( FOURIER_COMPONENT )
       ENDIF
      ENDIF

C  ###################################################
C  RT differential equation solutions + linearizations
C  ###################################################

C  Go to continuation point for thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) GO TO 8899

C  Start layer loop

      DO LAYER = 1, NLAYERS

C  Standard Solutions
C  ------------------

C  Get Discrete ordinate solutions for this layer

c      if ( do_write_screen) write(*,*)'homsolution',layer
        CALL LIDORT_HOM_SOLUTION
     I    ( LAYER, FOURIER_COMPONENT,
     O      STATUS_SUB )

C  .. error tracing

        IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
          MAIL = 'Error return from module LIDORT_HOM_SOLUTION'
          TRACE= 'Call #1 in LIDORT_L_FOURIER_MASTER'
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ENDIF

C  Get Post-processing ("user") solutions for this layer

        IF  ( STERM_LAYERMASK_UP(LAYER) .OR.
     &        STERM_LAYERMASK_DN(LAYER) ) THEN
          IF ( DO_USER_STREAMS ) THEN
            CALL LIDORT_HOM_USERSOLUTION
     I           ( LAYER, FOURIER_COMPONENT )
          ENDIF
        ENDIF

C  Additional Solutions for Linearization
C  --------------------------------------

        IF ( DO_ATMOS_LINEARIZATION ) THEN

C  Parameter control

          IF ( DO_PROFILE_LINEARIZATION ) THEN
            DO_RTSOL_VARY = LAYER_VARY_FLAG(LAYER)
            NPARAMS_VARY  = LAYER_VARY_NUMBER(LAYER)
          ELSE IF ( DO_COLUMN_LINEARIZATION ) THEN
            DO_RTSOL_VARY = .TRUE.
            NPARAMS_VARY  = N_TOTALCOLUMN_WFS
          ENDIF

C  Linearizations of Discrete Ordinate solution

c      if ( do_write_screen) write(*,*)'l_homsolution',layer
          CALL LIDORT_L_HOM_SOLUTION
     I      ( LAYER, FOURIER_COMPONENT,
     I        DO_RTSOL_VARY,
     I        NPARAMS_VARY,
     O        STATUS_SUB )

C  .. error tracing

          IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
            MAIL = 'Error return from module LIDORT_L_HOM_SOLUTION'
            TRACE= 'Call #2 in LIDORT_L_FOURIER_MASTER'
            STATUS = LIDORT_SERIOUS
            CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
            RETURN
          ENDIF

C  Get Linearizations of ("user") solutions for this layer

          CALL LIDORT_L_HOM_USERSOLUTION
     I      ( LAYER, FOURIER_COMPONENT,
     I        DO_RTSOL_VARY,
     I        NPARAMS_VARY )

C  End linearization control

        ENDIF

C  end layer loop

      ENDDO

C  Prepare Eigenvalue transmittances and any linearization

c      if ( do_write_screen) write(*,*)'eigentrans'
      CALL LIDORT_HOM_EIGENTRANS ( FOURIER_COMPONENT )
      IF ( DO_ATMOS_LINEARIZATION ) THEN
        CALL LIDORT_L_HOM_EIGENTRANS ( FOURIER_COMPONENT )
      ENDIF

C  Prepare solution norms if Green's function is in operation

      IF ( .NOT. DO_CLASSICAL_SOLUTION ) THEN
        CALL LIDORT_HOM_NORMS
        IF ( DO_ATMOS_LINEARIZATION ) THEN
          CALL LIDORT_L_HOM_NORMS
        ENDIF
      ENDIF

C  Prepare homogeneous solution multipliers

c      if ( do_write_screen) write(*,*)'hmult'
      CALL HMULT_MASTER
      IF ( DO_ATMOS_LINEARIZATION ) THEN
        CALL L_HMULT_MASTER
      ENDIF

C  ############################################
C   boundary value problem - MATRIX PREPARATION
C  ############################################

c      if ( do_write_screen) write(*,*)'bvpmatrix'

C  standard case using compression of band matrices, etc..

      IF ( BVP_REGULAR_FLAG(FOURIER_COMPONENT) ) THEN

        CALL BVP_MATRIXSETUP_MASTER
     I     ( DO_INCLUDE_SURFACE,
     I       FOURIER_COMPONENT,
     I       SURFACE_FACTOR,
     O       STATUS_SUB )

        IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
          MAIL = 'Error return from module BVPMATRIX_SETUP_MASTER'
          TRACE= 'Call #3A in LIDORT_L_FOURIER_MASTER'
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ENDIF

C  Telescoped case

      ELSE

        CALL BVPTEL_MATRIXSETUP_MASTER
     I     ( DO_INCLUDE_SURFACE,
     I       FOURIER_COMPONENT, 
     O       STATUS_SUB )

        IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
          MAIL = 'Error return from module BVPTELMATRIX_SETUP_MASTER'
          TRACE= 'Call #3B in LIDORT_L_FOURIER_MASTER'
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ENDIF

      ENDIF

C  ################
C  Thermal Solution (2 steps)
C  ################

C  Continuation point for avoiding the scattering calculations

 8899 continue

C  Separate calls if linearization is required

C  1. Find the Greens function solution.
C  2. Compute thermal layer source terms. (Upwelling and Downwelling)

C    These will be scaled up by factor 4.pi if solar beams as well

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        IF ( DO_ATMOS_LINEARIZATION ) THEN 
          CALL THERMAL_GFSOLUTION_PLUS
          IF ( DO_UPWELLING ) CALL THERMAL_STERMS_UP_PLUS
          IF ( DO_DNWELLING ) CALL THERMAL_STERMS_DN_PLUS
        ELSE
          CALL THERMAL_GFSOLUTION
          IF ( DO_UPWELLING ) CALL THERMAL_STERMS_UP
          IF ( DO_DNWELLING ) CALL THERMAL_STERMS_DN
        ENDIF
      ENDIF

C  Skip the thermal-only section if there are solar sources

      IF ( DO_SOLAR_SOURCES ) GO TO 455

C  ####################################################
C  Complete Radiation Field with Thermal-only solutions
C  ####################################################

C    Only one solution, local direct_beam flag NOT set

      IPARTIC = 1
      DO_INCLUDE_DIRECTBEAM   = .FALSE.

C  Avoid the scattering solutions if not flagged

      IF ( DO_THERMAL_TRANSONLY ) GO TO 566

C  Thermal-only. Find the BVP solution and intensity field
C  -------------------------------------------------------

C  set the BVP PI solution at the lower/upper boundaries

      DO LAYER = 1, NLAYERS
        DO I = 1, NSTREAMS_2
          WUPPER(I,LAYER) = T_WUPPER(I,LAYER)
          WLOWER(I,LAYER) = T_WLOWER(I,LAYER)
        ENDDO
      ENDDO

C  Solve the boundary value problem

      CALL BVP_SOLUTION_MASTER
     I       ( DO_INCLUDE_SURFACE,
     I         DO_INCLUDE_DIRECTBEAM,
     I         DO_INCLUDE_SURFEMISS,
     I         FOURIER_COMPONENT, IPARTIC, SURFACE_FACTOR,
     O         STATUS_SUB )

C  Error handling

      IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
        MAIL = 'Error return from module BVP_SOLUTION_MASTER'
        TRACE= 'Call #5A in LIDORT_FOURIER_MASTER'
        STATUS = LIDORT_SERIOUS
        CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
        RETURN
      ENDIF

C  Continuation point for avoiding thermal scattering

 566  CONTINUE

C  Post-processing - upwelling thermal-only field

      IF ( DO_UPWELLING ) THEN
        CALL UPUSER_INTENSITY
     I        ( DO_INCLUDE_SURFACE,
     I          DO_INCLUDE_SURFEMISS,
     I          DO_INCLUDE_THERMEMISS,
     I          DO_INCLUDE_MVOUTPUT,
     I          DO_INCLUDE_DIRECTBEAM,
     I          FOURIER_COMPONENT,
     I          IPARTIC,
     I          SURFACE_FACTOR,
     I          FLUX_MULTIPLIER )
      ENDIF

C  Post-processing - Downwelling thermal-only field

      IF ( DO_DNWELLING ) THEN
        CALL DNUSER_INTENSITY
     I        ( DO_INCLUDE_THERMEMISS,
     I          DO_INCLUDE_MVOUTPUT,
     I          FLUX_MULTIPLIER, 
     I          FOURIER_COMPONENT,
     I          IPARTIC )
      ENDIF

C  mean value output for thermal-only field

      IF ( DO_INCLUDE_MVOUTPUT ) THEN
        CALL MIFLUX_INTENSITY
     I          ( DO_INCLUDE_DIRECTBEAM, IPARTIC )
      ENDIF

C  Thermal-only: Atmospheric profile linearization
C  -----------------------------------------------

      IF ( DO_PROFILE_LINEARIZATION ) THEN

C  Start over layers that will contain variations
C    - Set flag for variation of layer, and number of variations

        DO LAYER_TO_VARY = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(LAYER_TO_VARY) ) THEN
            N_LAYER_WFS = LAYER_VARY_NUMBER(LAYER_TO_VARY)

C  Avoidance of BVP problem for transmittance only

            IF ( DO_THERMAL_TRANSONLY ) GO TO 788

C  Solve the Regular BVP linearization

            IF ( BVP_REGULAR_FLAG(FOURIER_COMPONENT) ) THEN

              CALL L_BVP_SOLUTION_MASTER
     I       ( DO_INCLUDE_SURFACE,
     I         DO_INCLUDE_SURFEMISS,
     I         DO_INCLUDE_THERMEMISS,
     I         DO_INCLUDE_DIRECTBEAM,
     I         LAYER_TO_VARY, N_LAYER_WFS,
     I         FOURIER_COMPONENT, IPARTIC, SURFACE_FACTOR,
     O         STATUS_SUB )

              IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
               MAIL = 'Error in module L_BVP_SOLUTION_MASTER, Prof/Lin'
               TRACE= 'Call #8A in LIDORT_L_FOURIER_MASTER'
               STATUS = LIDORT_SERIOUS
               CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
               RETURN
              ENDIF

            ENDIF

C  Continuation point for avoiding scattering solution

 788        CONTINUE

C  Post-processing for the weighting functions

C   if ( do_write_screen) write(*,*)'wf up',ibeam,layer_to_vary

            IF ( DO_UPWELLING ) THEN
              CALL UPUSER_ATMOSWF
     I         ( DO_INCLUDE_SURFACE,
     I           DO_INCLUDE_SURFEMISS,
     I           DO_INCLUDE_THERMEMISS,
     I           DO_INCLUDE_MVOUTPUT, 
     I           DO_INCLUDE_DIRECTBEAM,
     I           SURFACE_FACTOR,     
     I           FLUX_MULTIPLIER,
     I           FOURIER_COMPONENT,   IPARTIC,
     I           LAYER_TO_VARY,       N_LAYER_WFS )
            ENDIF

C  if ( do_write_screen) write(*,*)'wf dn',ibeam,layer_to_vary

            IF ( DO_DNWELLING ) THEN
               CALL DNUSER_ATMOSWF
     I         ( DO_INCLUDE_THERMEMISS,
     I           DO_INCLUDE_MVOUTPUT, 
     I           FLUX_MULTIPLIER,
     I           FOURIER_COMPONENT,   IPARTIC,
     I           LAYER_TO_VARY,       N_LAYER_WFS )
            ENDIF

            IF ( DO_INCLUDE_MVOUTPUT ) THEN
              CALL MIFLUX_ATMOSWF
     I         ( DO_INCLUDE_DIRECTBEAM, 
     I           IPARTIC, LAYER_TO_VARY, N_LAYER_WFS )
            ENDIF

C  Finish loop over layers with variation

          ENDIF
        ENDDO

C  End atmospheric profile weighting functions

      ENDIF

C  Thermal Only: Atmospheric Bulk weighting functions
C  --------------------------------------------------

      IF ( DO_COLUMN_LINEARIZATION ) THEN

C   variation index = 0

        VARIATION_INDEX = 0

C  Avoidance of BVP problem for transmittance only

        IF ( DO_THERMAL_TRANSONLY ) GO TO 789

C  Solve the linearized BVP (No telescoped solution)

        CALL L_BVP_SOLUTION_MASTER
     I       ( DO_INCLUDE_SURFACE,
     I         DO_INCLUDE_SURFEMISS,
     I         DO_INCLUDE_THERMEMISS,
     I         DO_INCLUDE_DIRECTBEAM,
     I         VARIATION_INDEX, N_TOTALCOLUMN_WFS,
     I         FOURIER_COMPONENT, IPARTIC, SURFACE_FACTOR,
     O         STATUS_SUB )

        IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
          MAIL = 'Error in module L_BVP_SOLUTION_MASTER, Col/Lin.'
          TRACE= 'Call #9A in LIDORT_L_FOURIER_MASTER'
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ENDIF

C  Continuation point for avoiding scattering solution

 789    CONTINUE

C  Post-processing for the weighting functions
C  -------------------------------------------

C  Sequence:
C    Upwelling   Atmospheric weighting functions
C    Downwelling Atmospheric weighting functions
C    Mean-value  Atmospheric weighting functions

        IF ( DO_UPWELLING ) THEN
          CALL UPUSER_ATMOSWF
     I         ( DO_INCLUDE_SURFACE,
     I           DO_INCLUDE_SURFEMISS,
     I           DO_INCLUDE_THERMEMISS,
     I           DO_INCLUDE_MVOUTPUT, 
     I           DO_INCLUDE_DIRECTBEAM,
     I           SURFACE_FACTOR,     
     I           FLUX_MULTIPLIER,
     I           FOURIER_COMPONENT,   IPARTIC,
     I           VARIATION_INDEX, N_TOTALCOLUMN_WFS )
        ENDIF

        IF ( DO_DNWELLING ) THEN
              CALL DNUSER_ATMOSWF
     I         ( DO_INCLUDE_THERMEMISS,
     I           DO_INCLUDE_MVOUTPUT, 
     I           FLUX_MULTIPLIER,
     I           FOURIER_COMPONENT,   IPARTIC,
     I           VARIATION_INDEX, N_TOTALCOLUMN_WFS )
        ENDIF

        IF ( DO_INCLUDE_MVOUTPUT ) THEN
          CALL MIFLUX_ATMOSWF
     I     ( DO_INCLUDE_DIRECTBEAM, 
     I       IPARTIC, VARIATION_INDEX, N_TOTALCOLUMN_WFS )
        ENDIF

C  End atmospheric column weighting functions

      ENDIF

C  Thermal-only: Surface Reflectance weighting functions
C  -----------------------------------------------------

      IF ( DO_SURFACE_LINEARIZATION ) THEN

C  initialise weighting function count

        WF_INDEX = 0

C  Only do this if there is surface reflectance !

        IF ( DO_INCLUDE_SURFACE ) THEN

C  Start loop over BRDF Kernels

          DO BRDF_INDEX = 1, N_BRDF_KERNELS

C  Factor weighting functions

            IF ( DO_KERNEL_FACTOR_WFS(BRDF_INDEX) ) THEN

              WF_INDEX = WF_INDEX + 1
              CALL KFACTOR_BRDFWF
     I            ( DO_INCLUDE_SURFACE,
     I              DO_INCLUDE_SURFEMISS,
     I              DO_INCLUDE_MVOUTPUT,
     I              DO_INCLUDE_DIRECTBEAM,
     I              FOURIER_COMPONENT, IPARTIC,
     I              SURFACE_FACTOR,
     I              FLUX_MULTIPLIER,
     I              BRDF_INDEX, WF_INDEX,
     O              STATUS_SUB )

              IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
                MAIL = 'Error in module LIDORT_KFACTOR_BRDFWF'
                TRACE= 'Call #10 in LIDORT_L_FOURIER_MASTER'
                STATUS = LIDORT_SERIOUS
                CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
                RETURN
              ENDIF

            ENDIF

C  Kernel parameter weighting functions

            DO PAR_INDEX = 1, N_BRDF_PARAMETERS(BRDF_INDEX)
              IF (DO_KERNEL_PARAMS_WFS(BRDF_INDEX,PAR_INDEX)) THEN

                WF_INDEX = WF_INDEX + 1
                CALL KPARAMS_BRDFWF
     I               ( DO_INCLUDE_SURFACE,
     I                 DO_INCLUDE_SURFEMISS,
     I                 DO_INCLUDE_MVOUTPUT, 
     I                 DO_INCLUDE_DIRECTBEAM,
     I                 FOURIER_COMPONENT, IPARTIC,
     I                 SURFACE_FACTOR,
     I                 FLUX_MULTIPLIER,
     I                 BRDF_INDEX, PAR_INDEX, WF_INDEX,
     O                 STATUS_SUB )

                IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
                  MAIL = 'Error in module LIDORT_KPARAMS_BRDFWF'
                  TRACE= 'Call #11 in LIDORT_L_FOURIER_MASTER'
                  STATUS = LIDORT_SERIOUS
                  CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
                  RETURN
                ENDIF

              ENDIF
            ENDDO

C  end of albedo weighting functions

          ENDDO
        ENDIF
      ENDIF

C  Finish Thermal only.

      RETURN

C  ##################################################
C  Complete Radiation Field with Solar Beam solutions
C  ##################################################

C  Continuation point

 455  CONTINUE

C  Start beam loop

      DO IBEAM = 1, NBEAMS

        IF ( DO_MULTIBEAM(IBEAM,FOURIER_COMPONENT) ) THEN

C  Step 1. Solar beam Particular solutions + linearizations
C  --------------------------------------------------------

          DO LAYER = 1, NLAYERS

C  Parameter control for the linearization

            IF ( DO_ATMOS_LINEARIZATION ) THEN
              IF ( DO_PROFILE_LINEARIZATION ) THEN
                DO_RTSOL_VARY = LAYER_VARY_FLAG(LAYER)
                NPARAMS_VARY  = LAYER_VARY_NUMBER(LAYER)
              ELSE IF ( DO_COLUMN_LINEARIZATION ) THEN
                DO_RTSOL_VARY = .TRUE.
                NPARAMS_VARY  = N_TOTALCOLUMN_WFS
              ENDIF
            ENDIF

C  A. the classical solution
C     ++++++++++++++++++++++

            IF ( DO_CLASSICAL_SOLUTION ) THEN

C  For the particular integral itself

c              if ( do_write_screen) write(*,*)'qbeam',layer
              CALL  LIDORT_QBEAM_SOLUTION
     I         ( LAYER, FOURIER_COMPONENT, IBEAM, STATUS_SUB )

              IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
                WRITE(C2,'(I2)')LAYER
                MAIL = 'Error from LIDORT_QBEAM_SOLUTION, layer'//C2
                TRACE= 'Call #5 in LIDORT_L_FOURIER_MASTER'
                STATUS = LIDORT_SERIOUS
                CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
                RETURN
              ENDIF

C  For the linearizations of the classical solution

              IF ( DO_ATMOS_LINEARIZATION ) THEN

c              if ( do_write_screen) write(*,*)'L_qbeam',layer
                CALL LIDORT_L_QBEAM_SOLUTION
     I            ( LAYER, FOURIER_COMPONENT, IBEAM,
     I              DO_RTSOL_VARY, 
     I              NPARAMS_VARY,
     O              STATUS_SUB )

                IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
                  MAIL =
     &             'Error return from module LIDORT_L_QBEAM_SOLUTION'
                  TRACE= 'Call #6 in LIDORT_L_FOURIER_MASTER'
                  STATUS = LIDORT_SERIOUS
                  CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
                  RETURN
                ENDIF

              ENDIF

C  B. the Green's function solution
C  ++++++++++++++++++++++++++++++++

            ELSE

C  For the particular solution itself

              CALL LIDORT_GBEAM_SOLUTION
     I      ( LAYER, FOURIER_COMPONENT, IBEAM )

C  For the linearization of this particular solution

              IF ( DO_ATMOS_LINEARIZATION ) THEN

                CALL LIDORT_L_GBEAM_SOLUTION
     I            ( LAYER, FOURIER_COMPONENT, IBEAM,
     I              DO_RTSOL_VARY, 
     I              NPARAMS_VARY )

              ENDIF

            ENDIF

C  C. user solutions
C  +++++++++++++++++

            IF  ( STERM_LAYERMASK_UP(LAYER) .OR.
     &            STERM_LAYERMASK_DN(LAYER) ) THEN
              IF ( DO_USER_STREAMS ) THEN

C  (classical)
                IF ( DO_CLASSICAL_SOLUTION ) THEN

                  CALL LIDORT_QBEAM_USERSOLUTION
     I             ( LAYER, FOURIER_COMPONENT, IBEAM )
     
                  IF ( DO_ATMOS_LINEARIZATION ) THEN
                    CALL LIDORT_L_QBEAM_USERSOLUTION
     I               ( LAYER, FOURIER_COMPONENT, IBEAM,
     I                 DO_RTSOL_VARY, 
     I                 NPARAMS_VARY )
                  ENDIF

C  (Green's function)

                ELSE
                  CALL LIDORT_GBEAM_USERSOLUTION
     I             ( LAYER, FOURIER_COMPONENT, IBEAM )

                  IF ( DO_ATMOS_LINEARIZATION ) THEN
                    CALL LIDORT_L_GBEAM_USERSOLUTION
     I               ( LAYER, FOURIER_COMPONENT, IBEAM,
     I                 DO_RTSOL_VARY, 
     I                 NPARAMS_VARY )
                  ENDIF
                ENDIF

              END IF
            END IF

C  end layer loop

          END DO

C  Add thermal solutions if flagged
C    Notice the 4.PI modulus on the thermal contribution.

          IF ( DO_INCLUDE_THERMEMISS ) THEN
           DO LAYER = 1, NLAYERS
            DO I = 1, NSTREAMS_2
             WUPPER(I,LAYER) = WUPPER(I,LAYER) + T_WUPPER(I,LAYER)
             WLOWER(I,LAYER) = WLOWER(I,LAYER) + T_WLOWER(I,LAYER)
            ENDDO
           ENDDO
          ENDIF

C  Adding the linearized thermal solutions is done later.....

C  Solve boundary value problem
C  -----------------------------

C  Get the REGULAR BVP solution for this beam component

c       if ( do_write_screen) write(*,*)'bvp solution',ibeam

          IF ( BVP_REGULAR_FLAG(FOURIER_COMPONENT) ) THEN

           CALL BVP_SOLUTION_MASTER
     I       ( DO_INCLUDE_SURFACE,
     I         DO_REFLECTED_DIRECTBEAM(IBEAM),
     I         DO_INCLUDE_SURFEMISS,
     I         FOURIER_COMPONENT, IBEAM, SURFACE_FACTOR,
     O         STATUS_SUB )

           IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
            MAIL = 'Error return from module BVP_SOLUTION_MASTER'
            TRACE= 'Call #7A in LIDORT_L_FOURIER_MASTER'
            STATUS = LIDORT_SERIOUS
            CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
            RETURN
           ENDIF

C  Get the telescoped boundary value result

          ELSE

           CALL BVPTEL_SOLUTION_MASTER
     I       ( FOURIER_COMPONENT, IBEAM,
     O         STATUS_SUB )

           IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
            MAIL = 'Error return from module BVPTEL_SOLUTION_MASTER'
            TRACE= 'Call #7B in LIDORT_L_FOURIER_MASTER'
            STATUS = LIDORT_SERIOUS
            CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
            RETURN
           ENDIF

          ENDIF

C  Radiance Field Post Processing
C  ------------------------------

C  Direct beam inclusion flag:
C   This now has the DBCORRECTION option: if the DBCORRECTION Flag
C   is set, then we will be doing exact calculations of the reflected
C   directbeam, so we do not need to include it in the Post-processing.
C   However, the direct beam will need to be included in the basic RT
C   solution (the BVP), and this is controlled separately by the
C   DO_REFLECTED_DIRECTBEAM(IBEAM) flags.
C     R. Spurr, RT Solutions, Inc., 19 August 2005.

          DO_INCLUDE_DIRECTBEAM = ( DO_UPWELLING .AND.
     &     (DO_REFLECTED_DIRECTBEAM(IBEAM).AND..NOT.DO_DBCORRECTION))

C  former code
c          DO_INCLUDE_DIRECTBEAM = ( DO_UPWELLING .AND.
c     &     (DO_REFLECTED_DIRECTBEAM(IBEAM)))

C  upwelling

          IF ( DO_UPWELLING ) THEN
c            if ( do_write_screen) write(*,*)'intensity-up ',ibeam
            CALL UPUSER_INTENSITY
     I        ( DO_INCLUDE_SURFACE,
     I          DO_INCLUDE_SURFEMISS,
     I          DO_INCLUDE_THERMEMISS,
     I          DO_INCLUDE_MVOUTPUT,
     I          DO_INCLUDE_DIRECTBEAM,
     I          FOURIER_COMPONENT,
     I          IBEAM,
     I          SURFACE_FACTOR,
     I          FLUX_MULTIPLIER )
          ENDIF

C  Downwelling

          IF ( DO_DNWELLING ) THEN
c            if ( do_write_screen) write(*,*)'intensity-dn ',ibeam
            CALL DNUSER_INTENSITY
     I        ( DO_INCLUDE_THERMEMISS,
     I          DO_INCLUDE_MVOUTPUT,
     I          FLUX_MULTIPLIER, 
     I          FOURIER_COMPONENT,
     I          IBEAM )
          ENDIF

C  mean value output

          IF ( DO_INCLUDE_MVOUTPUT ) THEN
            CALL MIFLUX_INTENSITY
     I          ( DO_INCLUDE_DIRECTBEAM, IBEAM )
          ENDIF

C  Finished this Beam solution, if only intensity is required

          IF ( DO_SIMULATION_ONLY ) GO TO 4000

C  Step 3. Atmospheric Profile weighting function loop
C  ---------------------------------------------------

          IF ( DO_PROFILE_LINEARIZATION ) THEN

C  Start over layers that will contain variations
C    - Set flag for variation of layer, and number of variations

           DO LAYER_TO_VARY = 1, NLAYERS
            IF ( LAYER_VARY_FLAG(LAYER_TO_VARY) ) THEN
             N_LAYER_WFS = LAYER_VARY_NUMBER(LAYER_TO_VARY)

C  3A. Solve the linearized BVP
C  ============================

C  Get the linearized BVP solution for this beam component

C            if ( do_write_screen) write(*,*)'l_bvp',ibeam,layer_to_vary

C  Regular BVP linearization

             IF ( BVP_REGULAR_FLAG(FOURIER_COMPONENT) ) THEN

              CALL L_BVP_SOLUTION_MASTER
     I       ( DO_INCLUDE_SURFACE,
     I         DO_INCLUDE_SURFEMISS,
     I         DO_INCLUDE_THERMEMISS,
     I         DO_REFLECTED_DIRECTBEAM(IBEAM),
     I         LAYER_TO_VARY, N_LAYER_WFS,
     I         FOURIER_COMPONENT, IBEAM, SURFACE_FACTOR,
     O         STATUS_SUB )

              IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
               MAIL = 'Error in module L_BVP_SOLUTION_MASTER'
               TRACE= 'Call #8A in LIDORT_L_FOURIER_MASTER'
               STATUS = LIDORT_SERIOUS
               CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
               RETURN
              ENDIF

C  telescoped BVP linearization

             ELSE

              CALL L_BVPTEL_SOLUTION_MASTER
     I       ( LAYER_TO_VARY, N_LAYER_WFS,
     I         FOURIER_COMPONENT, IBEAM,
     O         STATUS_SUB )

              IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
               MAIL = 'Error in module L_BVPTEL_SOLUTION_MASTER'
               TRACE= 'Call #8B in LIDORT_L_FOURIER_MASTER'
               STATUS = LIDORT_SERIOUS
               CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
               RETURN
              ENDIF

             ENDIF

C  3B. Post-processing for the weighting functions
C  ===============================================

C  Sequence:
C    Upwelling   Atmospheric weighting functions
C    Downwelling Atmospheric weighting functions
C    Mean-value  Atmospheric weighting functions

             IF ( DO_UPWELLING ) THEN

C            if ( do_write_screen) write(*,*)'wf up',ibeam,layer_to_vary
 
              CALL UPUSER_ATMOSWF
     I         ( DO_INCLUDE_SURFACE,
     I           DO_INCLUDE_SURFEMISS,
     I           DO_INCLUDE_THERMEMISS,
     I           DO_INCLUDE_MVOUTPUT, 
     I           DO_INCLUDE_DIRECTBEAM,
     I           SURFACE_FACTOR,     
     I           FLUX_MULTIPLIER,
     I           FOURIER_COMPONENT,   IBEAM,
     I           LAYER_TO_VARY,       N_LAYER_WFS )
             ENDIF

             IF ( DO_DNWELLING ) THEN
 
C           if ( do_write_screen) write(*,*)'wf dn',ibeam,layer_to_vary

              CALL DNUSER_ATMOSWF
     I         ( DO_INCLUDE_THERMEMISS,
     I           DO_INCLUDE_MVOUTPUT, 
     I           FLUX_MULTIPLIER,
     I           FOURIER_COMPONENT,   IBEAM,
     I           LAYER_TO_VARY,       N_LAYER_WFS )
             ENDIF

             IF ( DO_INCLUDE_MVOUTPUT ) THEN
              CALL MIFLUX_ATMOSWF
     I         ( DO_INCLUDE_DIRECTBEAM, 
     I           IBEAM, LAYER_TO_VARY, N_LAYER_WFS )
             ENDIF

C  Finish loop over layers with variation

            ENDIF
           ENDDO

C  End atmospheric weighting functions

          ENDIF

C  Step 4. Atmospheric Bulk weighting functions
C  ============================================

          IF ( DO_COLUMN_LINEARIZATION ) THEN

C   variation index = 0

            VARIATION_INDEX = 0

C  4A. Solve the linearized BVP
C  ----------------------------

C  Regular BVP linearization

            IF ( BVP_REGULAR_FLAG(FOURIER_COMPONENT) ) THEN

C  Get the linearized BVP solution for this beam component

              CALL L_BVP_SOLUTION_MASTER
     I       ( DO_INCLUDE_SURFACE,
     I         DO_INCLUDE_SURFEMISS,
     I         DO_INCLUDE_THERMEMISS,
     I         DO_REFLECTED_DIRECTBEAM(IBEAM),
     I         VARIATION_INDEX, N_TOTALCOLUMN_WFS,
     I         FOURIER_COMPONENT, IBEAM, SURFACE_FACTOR,
     O         STATUS_SUB )

              IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
                MAIL = 'Error in module L_BVP_SOLUTION_MASTER, Col/Lin.'
                TRACE= 'Call #9A in LIDORT_L_FOURIER_MASTER'
                STATUS = LIDORT_SERIOUS
                CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
                RETURN
              ENDIF

C  telescoped BVP linearization

             ELSE

              CALL L_BVPTEL_SOLUTION_MASTER
     I       ( VARIATION_INDEX, N_TOTALCOLUMN_WFS,
     I         FOURIER_COMPONENT, IBEAM,
     O         STATUS_SUB )

              IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
               MAIL = 'Error in module L_BVPTEL_SOLUTION_MASTER'
               TRACE= 'Call #9B in LIDORT_L_FOURIER_MASTER'
               STATUS = LIDORT_SERIOUS
               CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
               RETURN
              ENDIF

             ENDIF

C  4B. Post-processing for the weighting functions
C  -----------------------------------------------

C  Sequence:
C    Upwelling   Atmospheric weighting functions
C    Downwelling Atmospheric weighting functions
C    Mean-value  Atmospheric weighting functions

            IF ( DO_UPWELLING ) THEN
C            if ( do_write_screen) write(*,*)'wf up',ibeam,layer_to_vary 
              CALL UPUSER_ATMOSWF
     I         ( DO_INCLUDE_SURFACE,
     I           DO_INCLUDE_SURFEMISS,
     I           DO_INCLUDE_THERMEMISS,
     I           DO_INCLUDE_MVOUTPUT,
     I           DO_INCLUDE_DIRECTBEAM,
     I           SURFACE_FACTOR,      FLUX_MULTIPLIER,
     I           FOURIER_COMPONENT,   IBEAM,
     I           VARIATION_INDEX, N_TOTALCOLUMN_WFS )
            ENDIF

            IF ( DO_DNWELLING ) THEN
C           if ( do_write_screen) write(*,*)'wf dn',ibeam,layer_to_vary
              CALL DNUSER_ATMOSWF
     I         ( DO_INCLUDE_THERMEMISS,
     I           DO_INCLUDE_MVOUTPUT, 
     I           FLUX_MULTIPLIER,
     I           FOURIER_COMPONENT,   IBEAM,
     I           VARIATION_INDEX, N_TOTALCOLUMN_WFS )
            ENDIF

            IF ( DO_INCLUDE_MVOUTPUT ) THEN
              CALL MIFLUX_ATMOSWF
     I         ( DO_INCLUDE_DIRECTBEAM,
     I           IBEAM, VARIATION_INDEX, N_TOTALCOLUMN_WFS )
            ENDIF

C  End atmospheric column weighting functions

          ENDIF

C  Step 5. Surface Reflectance weighting functions
C  ===============================================

          IF ( DO_SURFACE_LINEARIZATION ) THEN

C  initialise weighting function count

            WF_INDEX = 0

C  Only do this if there is surface reflectance !

            IF ( DO_INCLUDE_SURFACE ) THEN

C  Start loop over BRDF Kernels

              DO BRDF_INDEX = 1, N_BRDF_KERNELS

C  Factor weighting functions

                IF ( DO_KERNEL_FACTOR_WFS(BRDF_INDEX) ) THEN

                  WF_INDEX = WF_INDEX + 1
                  CALL KFACTOR_BRDFWF
     I            ( DO_INCLUDE_SURFACE,
     I              DO_INCLUDE_SURFEMISS,
     I              DO_INCLUDE_MVOUTPUT,
     I              DO_INCLUDE_DIRECTBEAM,
     I              FOURIER_COMPONENT, IBEAM,
     I              SURFACE_FACTOR,
     I              FLUX_MULTIPLIER,
     I              BRDF_INDEX, WF_INDEX,
     O              STATUS_SUB )

                  IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
                    MAIL = 'Error in module LIDORT_KFACTOR_BRDFWF'
                    TRACE= 'Call #10 in LIDORT_L_FOURIER_MASTER'
                    STATUS = LIDORT_SERIOUS
                    CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
                    RETURN
                  ENDIF

                ENDIF

C  Kernel parameter weighting functions

                DO PAR_INDEX = 1, N_BRDF_PARAMETERS(BRDF_INDEX)
  
                  IF (DO_KERNEL_PARAMS_WFS(BRDF_INDEX,PAR_INDEX)) THEN

                    WF_INDEX = WF_INDEX + 1
                    CALL KPARAMS_BRDFWF
     I               ( DO_INCLUDE_SURFACE,
     I                 DO_INCLUDE_SURFEMISS,
     I                 DO_INCLUDE_MVOUTPUT, 
     I                 DO_INCLUDE_DIRECTBEAM,
     I                 FOURIER_COMPONENT, IBEAM,
     I                 SURFACE_FACTOR,
     I                 FLUX_MULTIPLIER,
     I                 BRDF_INDEX, PAR_INDEX, WF_INDEX,
     O                 STATUS_SUB )

                    IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
                      MAIL = 'Error in module LIDORT_KPARAMS_BRDFWF'
                      TRACE= 'Call #11 in LIDORT_L_FOURIER_MASTER'
                      STATUS = LIDORT_SERIOUS
                      CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
                      RETURN
                    ENDIF

                  ENDIF
                ENDDO

C  end of albedo weighting functions

              ENDDO
            ENDIF
          ENDIF

C  Step 6. SurfBB weighting functions
C  ----------------------------------

c          IF ( DO_INCLUDE_SURFEMISS ) THEN
c            IF ( DO_SURFBB_LINEARIZATION ) THEN
c              CALL LIDORT_SURFBBWF
c     I            ( DO_INCLUDE_SURFACE,
c     I              DO_INCLUDE_MVOUTPUT,
c     I              FOURIER_COMPONENT, IBEAM,
c     I              SURFACE_FACTOR,
c     I              FLUX_MULTIPLIER,
c     O              STATUS_SUB )
c              IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
c                MAIL = 'Error in module LIDORT_SURFBBWF'
c                TRACE= 'Call #12 in LIDORT_L_FOURIER_MASTER'
c                STATUS = LIDORT_SERIOUS
c                CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
c                RETURN
c              ENDIF
c            ENDIF
c          ENDIF

C  Continuation point for avoiding weighting functions

 4000     CONTINUE

C  End loop over beam solutions

        END IF
      END DO

C  ######
C  finish
C  ######

      RETURN
      END
