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
C #            LIDORT_FOURIER_MASTER (master)                   #
C #                                                             #
C ###############################################################

      SUBROUTINE LIDORT_FOURIER_MASTER
     I       ( FOURIER_COMPONENT,
     I         SS_FLUX_MULTIPLIER,
     O         DO_INCLUDE_SURFACE,
     O         DO_INCLUDE_SURFEMISS,
     O         DO_INCLUDE_THERMEMISS,
     O         STATUS )

C  Complete Fourier component calculation for the Standard Code.

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include file of input variables
C  Include file of bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  Include files of setup variables (output from SETUPS routine)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_THERMALSUP.VARS'

C  include the solution file for debug only

      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  include this file for results (debug only)

      INCLUDE '../includes/LIDORT_RESULTS.VARS'

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

      INTEGER          LAYER, IBEAM, IPARTIC, UM, NLEVEL, I

C  local inclusion flags

      LOGICAL          DO_INCLUDE_MVOUTPUT
      LOGICAL          DO_INCLUDE_DIRECTBEAM

C  Flux multiplier and Fourier component numbers

      DOUBLE PRECISION FLUX_MULTIPLIER
      DOUBLE PRECISION DELTA_FACTOR
      DOUBLE PRECISION SURFACE_FACTOR

C  error tracing

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

C  inclusion of thermal surface emission term. only for Fourier = 0

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

C  Setups for Fourier = 0

C   MISCSETUPS    : Delta-M, average-secant formulation, transmittances
C   THERMAL_SETUP : Coefficients, direct multipliers
C   EMULT_MASTER  : Beam source function multipliers. Not required for the
C                  Full SS calculation in outgoing mode

      IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
        CALL LIDORT_MISCSETUPS
        IF ( DO_INCLUDE_THERMEMISS ) CALL THERMAL_SETUP
        IF ( DO_SOLAR_SOURCES      ) THEN
          IF (.NOT.DO_SSFULL.OR.(DO_SSFULL.AND.DO_SSCORR_NADIR)) THEN
            CALL EMULT_MASTER
          ENDIF
        ENDIF
      ENDIF

C  #####################
C  Correction operations 
C  #####################

C  Not required if no solar sources

      IF ( DO_SOLAR_SOURCES ) THEN

C  Single scatter correction (pre-calculation)
C     Must be done after MISCSETUPS and EMULT_MASTER, as we need
C      multipliers and transmittance factors for SUN and LOS paths.
C      Code added 6 May 2005. Replaces call in Master routine.
C      Version 3.1. Added call to the new outgoing sphericity correction

       IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
        IF ( DO_USER_STREAMS ) THEN

C  regular nadir-scattering SS correction

          IF ( DO_SSCORR_NADIR ) THEN
            CALL LIDORT_SSCORR_NADIR ( SS_FLUX_MULTIPLIER )
          ENDIF

C  New outgoing sphericity correction

          IF ( DO_SSCORR_OUTGOING ) THEN
            CALL LIDORT_SSCORR_OUTGOING
     &     ( SS_FLUX_MULTIPLIER, FAIL, MAIL )
            IF ( FAIL ) THEN
              TRACE= 'SS correction outgoing failed'
              STATUS = LIDORT_SERIOUS
              CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
              RETURN
            ENDIF 
          ENDIF

C  End Fourier = 0 and User streams clauses

        ENDIF
       ENDIF

C  Exact direct beam BRDF (if flagged)
C   Must be done after MISCSETUPS, as we need the SUN transmittance factors
C   BRDF Kernels already known (called before this module)
C   Code added 6 May 2005.

       IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
        IF ( DO_USER_STREAMS ) THEN
          IF ( DO_LAMBERTIAN_SURFACE.AND.DO_SSFULL ) THEN
            CALL LIDORT_LAMBERTIAN_DBCORRECTION (SS_FLUX_MULTIPLIER)
            RETURN
          ELSE IF ( DO_LAMBERTIAN_SURFACE.AND.DO_DBCORRECTION ) THEN
            CALL LIDORT_LAMBERTIAN_DBCORRECTION (SS_FLUX_MULTIPLIER)
          ELSE
            IF ( DO_DBCORRECTION .OR. DO_SSFULL ) THEN
              CALL LIDORT_DBCORRECTION ( SS_FLUX_MULTIPLIER )
            ENDIF
          ENDIF
          IF ( DO_SSFULL ) RETURN
        ENDIF
       ENDIF

C  End solar sources only clause for corrections

      ENDIF

C  BRDF Fourier components - coefficient setup
C   (also computes emissivity if required)
C   (DO_REFLECTED_DIRECTBEAM is a stored variable)
C   Should always call this, even with Lambertian surfaces!

      CALL LIDORT_BRDF_FOURIER
     I    ( DO_INCLUDE_SURFACE,
     I      DO_INCLUDE_SURFEMISS,
     I      FOURIER_COMPONENT,
     I      DELTA_FACTOR )
      
C  Reflected Direct beam attenuation
C   (DO_REFLECTED_DIRECTBEAM is a stored variable)

c      if ( do_write_screen) write(*,*)'directbeam'

      IF ( DO_SOLAR_SOURCES ) THEN
        CALL LIDORT_DIRECTBEAM
     I    ( FOURIER_COMPONENT,
     I      DELTA_FACTOR )
      ENDIF

C  Get Legendre polynomials for this Fourier component
C    Including all beam angles !!

      IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
       CALL LIDORT_LEGENDRE_SETUP ( FOURIER_COMPONENT )
       IF ( DO_USER_STREAMS ) THEN
         CALL LIDORT_USERLEGENDRE_SETUP ( FOURIER_COMPONENT )
       ENDIF
      ENDIF

C  ########################################
C  RT differential equation Eigensolutions
C  ########################################

C  Go to continuation point for thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) GO TO 8899

C  Start layer loop

      DO LAYER = 1, NLAYERS

C  Get Discrete ordinate solutions for this layer

c      if ( do_write_screen) write(*,*)'homsolution',layer
        CALL LIDORT_HOM_SOLUTION
     I    ( LAYER, FOURIER_COMPONENT,
     O      STATUS_SUB )

C  .. error tracing

        IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
          MAIL = 'Error return from module LIDORT_HOM_SOLUTION'
          TRACE= 'Call #1 in LIDORT_FOURIER_MASTER'
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

C  end layer loop

      ENDDO

C  prepare eigenstream tranmsittances

c      if ( do_write_screen) write(*,*)'eigentrans'
      CALL LIDORT_HOM_EIGENTRANS
     I   ( FOURIER_COMPONENT )

C  Prepare solution norms if Green's function is in operation

      IF ( .NOT. DO_CLASSICAL_SOLUTION ) THEN
        CALL LIDORT_HOM_NORMS
      ENDIF

C  Prepare homogeneous solution multipliers

c      if ( do_write_screen) write(*,*)'hmult'
      CALL HMULT_MASTER

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
          TRACE= 'Call #2A in LIDORT_FOURIER_MASTER'
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
          MAIL = 'Error return from module BVPTEL_MATRIXSETUP_MASTER'
          TRACE= 'Call #2B in LIDORT_FOURIER_MASTER'
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ENDIF

      ENDIF

C  ################
C  Thermal Solution
C  ################

C  Continuation point for avoiding the scattering calculations

 8899 continue

C  1. Find the Green's function solution (also transmittance only)
C  2. Compute thermal layer source terms. (Upwelling and Downwelling)

C    These will be scaled up by factor 4.pi if solar beams as well

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        CALL THERMAL_GFSOLUTION
        IF ( DO_UPWELLING ) CALL THERMAL_STERMS_UP
        IF ( DO_DNWELLING ) CALL THERMAL_STERMS_DN
      ENDIF

C  Skip the thermal-only section if there are solar sources

      IF ( DO_SOLAR_SOURCES ) GO TO 455

C  ####################################################
C  Complete Radiation Field with Thermal-only solutions
C  ####################################################

C  Only one solution, local direct_beam flag NOT set

      IPARTIC = 1
      DO_INCLUDE_DIRECTBEAM = .FALSE.

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

C  Finish Thermal only.

      RETURN

C  ##################################################
C  Complete Radiation Field with Solar Beam solutions
C  ##################################################

C  Continuation point

 455  CONTINUE

C  Start loop over varisou solar beams

      DO IBEAM = 1, NBEAMS

C  Only calculate if still not converged

        IF ( DO_MULTIBEAM(IBEAM,FOURIER_COMPONENT) ) THEN

C  Solar beam Particular solutions
C  -------------------------------

C  start layer loop

          DO LAYER = 1, NLAYERS

C  get the classical solution

            IF ( DO_CLASSICAL_SOLUTION ) THEN

              CALL  LIDORT_QBEAM_SOLUTION
     I         ( LAYER, FOURIER_COMPONENT, IBEAM, STATUS_SUB )

              IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
                WRITE(C2,'(I2)')LAYER
                MAIL = 'Error from LIDORT_QBEAM_SOLUTION, layer'//C2
                TRACE= 'Call #4 in LIDORT_FOURIER_MASTER'
                STATUS = LIDORT_SERIOUS
                CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
                RETURN
              ENDIF

C  or get the Green's function solution

            ELSE

              CALL LIDORT_GBEAM_SOLUTION
     I      ( LAYER, FOURIER_COMPONENT, IBEAM )

            ENDIF

C  user solutions

            IF  ( STERM_LAYERMASK_UP(LAYER) .OR.
     &            STERM_LAYERMASK_DN(LAYER) ) THEN
              IF ( DO_USER_STREAMS ) THEN
                IF ( DO_CLASSICAL_SOLUTION ) THEN
                  CALL LIDORT_QBEAM_USERSOLUTION
     I             ( LAYER, FOURIER_COMPONENT, IBEAM )
                ELSE
                  CALL LIDORT_GBEAM_USERSOLUTION
     I             ( LAYER, FOURIER_COMPONENT, IBEAM )
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

C  Solve boundary value problem
C  ----------------------------

c       if ( do_write_screen) write(*,*)'bvp solution',ibeam

C  Get the Regular BVP solution

          IF ( BVP_REGULAR_FLAG(FOURIER_COMPONENT) ) THEN

           CALL BVP_SOLUTION_MASTER
     I       ( DO_INCLUDE_SURFACE,
     I         DO_REFLECTED_DIRECTBEAM(IBEAM),
     I         DO_INCLUDE_SURFEMISS,
     I         FOURIER_COMPONENT, IBEAM, SURFACE_FACTOR,
     O         STATUS_SUB )

           IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
            MAIL = 'Error return from module BVP_SOLUTION_MASTER'
            TRACE= 'Call #5A in LIDORT_FOURIER_MASTER'
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
            TRACE= 'Call #5B in LIDORT_FOURIER_MASTER'
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

C  End loop over beam solutions

        END IF
      END DO

C  IF0 Rayleigh assignment
C  -----------------------

      IF ( DO_RAYLEIGH_ONLY .AND. DO_UPWELLING ) THEN
       IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
         NLEVEL = UTAU_LEVEL_MASK_UP(1)
         IF ( NLEVEL .EQ. 0 ) THEN
           DO IBEAM = 1, NBEAMS
             DO UM = 1, N_USER_STREAMS
c               write(*,*)INTENSITY_F(1,UM,IBEAM,UPIDX)
               IF0_RAYLEIGH(IBEAM,UM) = INTENSITY_F(1,UM,IBEAM,UPIDX)
             ENDDO
           ENDDO
         ENDIF
       ENDIF
      ENDIF

C  ######
C  finish
C  ######

      RETURN
      END

