C ###############################################################
C #                                                             #
C #                    THE VECTOR LIDORT MODEL                  #
C #                                                             #
C #  (Vector LInearized Discrete Ordinate Radiative Transfer)   #
C #   -      --         -        -        -         -           #
C #                                                             #
C ###############################################################

C ###############################################################
C #                                                             #
C #  Author :      Robert. J. D. Spurr                          #
C #                                                             #
C #  Address :      RT Solutions, inc.                          #
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
C #            VLIDORT_MASTER_LAMBERTIAN (master)               #
C #            VLIDORT_FOURIER_LAMBERTIAN (master)              #
C #                                                             #
C ###############################################################
 
      SUBROUTINE VLIDORT_MASTER_LAMBERTIAN
     &     ( STATUS_INPUTCHECK, STATUS_CALCULATION )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  subroutine arguments
C  --------------------

      INTEGER          STATUS_INPUTCHECK
      INTEGER          STATUS_CALCULATION

C  local variables
C  ---------------

      LOGICAL          LOCAL_DO_NO_AZIMUTH
      LOGICAL          SAVE_DO_NO_AZIMUTH
      LOGICAL          LOCAL_ITERATION
      LOGICAL          DO_INCLUDE_SURFACE
      LOGICAL          DO_INCLUDE_SURFEMISS
      LOGICAL          DO_INCLUDE_THERMEMISS

      INTEGER          FOURIER_COMPONENT
      INTEGER          N_FOURIER_COMPONENTS

      INTEGER          UM, UA, TESTCONV, L
      INTEGER          LOCAL_N_USERAZM, STATUS_SUB
      INTEGER          IUNIT, SUNIT, FUNIT, RUNIT

      DOUBLE PRECISION AZM_ARGUMENT, DFC
      DOUBLE PRECISION AZMFAC
     &     (MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS,MAXSTOKES)

      INTEGER          IBEAM_COUNT, IBEAM, IB, NSOURCES
      LOGICAL          BEAM_ITERATION ( MAXBEAMS )
      INTEGER          BEAM_TESTCONV  ( MAXBEAMS )

      LOGICAL          ADJUST_SURFACE

      DOUBLE PRECISION SS_FLUX_MULTIPLIER, modified_eradius

C  Local error handling

      LOGICAL          FAIL
      CHARACTER*70     MAIL, TRACE
      CHARACTER*3      CF

C  initialize output status
C  ------------------------

      STATUS_CALCULATION = VLIDORT_SUCCESS
      STATUS_INPUTCHECK  = VLIDORT_SUCCESS

C  set multiplier

      SS_FLUX_MULTIPLIER = FLUX_FACTOR / PI4

C  Check input
C  -----------

      CALL VLIDORT_CHECK_INPUT ( STATUS_SUB )

      IF ( STATUS_SUB .EQ. VLIDORT_SERIOUS ) THEN
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
        MAIL = ' VLIDORT_CHECK_INPUT failed'
        TRACE = ' Called in VLIDORT_MASTER_LAMBERTIAN '
        CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS_INPUTCHECK )
        CLOSE(VLIDORT_ERRUNIT)
        RETURN
      ELSE IF ( STATUS_SUB .EQ. VLIDORT_WARNING ) THEN
        STATUS_INPUTCHECK = VLIDORT_WARNING
        MAIL = ' VLIDORT_CHECK_INPUT warning message'
        TRACE = ' Called in VLIDORT_MASTER_LAMBERTIAN '
        CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS_INPUTCHECK )
      ENDIF

C  Geometry adjustment
C  -------------------

C  Adjust surface condition

      ADJUST_SURFACE = .FALSE.
      IF ( DO_SSCORR_OUTGOING ) THEN
        IF (HEIGHT_GRID(NLAYERS).GT.GEOMETRY_SPECHEIGHT ) THEN
         ADJUST_SURFACE = .TRUE.
        ENDIF
      ENDIF

C  Perform adjustment

      modified_eradius = earth_radius + GEOMETRY_SPECHEIGHT
      CALL multi_outgoing_adjustgeom
     i   ( MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS,
     i     N_USER_VZANGLES,   N_SZANGLES,   N_USER_RELAZMS,
     i     height_grid(nlayers), modified_eradius, adjust_surface,
     i     user_vzangles,  szangles, user_relazms,
     o     user_vzangles_adjust, szangles_adjust, user_relazms_adjust,
     o     fail, mail )
      if ( fail ) return       

C  Chapman function calculation of TAUTHICK_INPUT
C  ----------------------------------------------

      IF ( DO_SOLAR_SOURCES ) THEN
       IF ( DO_CHAPMAN_FUNCTION ) THEN
        CALL VLIDORT_CHAPMAN ( FAIL, MAIL, TRACE )
        IF (FAIL) THEN
          STATUS_CALCULATION = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS_CALCULATION )
          CLOSE(VLIDORT_ERRUNIT)
          RETURN
        ENDIF
       ENDIF
      ENDIF

C  write input variables
C  ---------------------

C  open file, call standard input write, close file

      IF ( DO_WRITE_INPUT ) THEN
        IUNIT = VLIDORT_INUNIT
        OPEN(IUNIT,FILE=INPUT_WRITE_FILENAME,STATUS='UNKNOWN')
        CALL VLIDORT_WRITEINPUT ( IUNIT )
        CLOSE(IUNIT)
      ENDIF

C  Get derived inputs
C  ==================

C  Miscellaneous and layer input.

      CALL VLIDORT_DERIVE_INPUT ( STATUS_SUB, MAIL )

      IF ( STATUS_SUB .EQ. VLIDORT_WARNING ) THEN
        STATUS_INPUTCHECK = VLIDORT_WARNING
        TRACE = 'Derive_Input Call in VLIDORT_MASTER_LAMBERTIAN'
        CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS_INPUTCHECK )
      ENDIF

C  BRDF inputs. Commented out in the Lambertian-only version.

c      IF ( .NOT. DO_LAMBERTIAN_SURFACE ) THEN
c        CALL VLIDORT_BRDF_CREATOR
c      ENDIF

C  Initialise Fourier loop
C  =======================

C  Single scatter correction: flux multiplier

      IF ( DO_SSCORR_NADIR .OR. DO_SSCORR_OUTGOING ) THEN
        SS_FLUX_MULTIPLIER = FLUX_FACTOR / PI4
      ELSE
        SS_FLUX_MULTIPLIER = FLUX_FACTOR / PI4
      ENDIF

C  Set Number of Fourier terms (NMOMENTS = Maximum).
C    ( Starting from 0 = Fundamental )

      SAVE_DO_NO_AZIMUTH  = DO_NO_AZIMUTH
      LOCAL_DO_NO_AZIMUTH = DO_NO_AZIMUTH

C  No azimuth dependency for following cases
C    Thermal-only treatment is isotropic, no azimuth

C      IF ( DO_TRANSMITTANCE_ONLY ) THEN
C        LOCAL_DO_NO_AZIMUTH = .TRUE.
C      ENDIF

      IF ( .NOT. DO_SOLAR_SOURCES  ) THEN
        LOCAL_DO_NO_AZIMUTH = .TRUE.
      ENDIF

      IF ( DO_MVOUT_ONLY  ) THEN
        LOCAL_DO_NO_AZIMUTH = .TRUE.
      ENDIF

C  set Fourier number (2 for Rayleigh only)

      IF ( LOCAL_DO_NO_AZIMUTH .OR. DO_SSFULL ) THEN
        N_FOURIER_COMPONENTS = 0
      ELSE
        IF ( DO_RAYLEIGH_ONLY  ) THEN
          N_FOURIER_COMPONENTS = 2
        ELSE
          N_FOURIER_COMPONENTS = NMOMENTS
        ENDIF
      ENDIF

C  re-set no-azimuth flag

      DO_NO_AZIMUTH = LOCAL_DO_NO_AZIMUTH

C  if there's no azimuth dependence, just do one value in azimuth loop

      IF ( DO_NO_AZIMUTH ) THEN
        LOCAL_N_USERAZM = 1
      ELSE
        LOCAL_N_USERAZM = N_USER_RELAZMS
      ENDIF

C  Number of sources

      IF ( DO_SOLAR_SOURCES ) THEN
        NSOURCES = N_SZANGLES
      ELSE
        NSOURCES = 1
        NBEAMS   = 1
      ENDIF

C  save some offsets for indexing geometries
C     These results are now stored in VLIDORT_BOOKKEEP.VARS

      N_VIEWING    = N_USER_VZANGLES * LOCAL_N_USERAZM
      N_GEOMETRIES = NSOURCES * N_VIEWING
      DO IBEAM = 1, NBEAMS
        SZA_OFFSETS(IBEAM) = N_VIEWING * ( IBEAM - 1 )
        DO UM = 1, N_USER_VZANGLES
          VZA_OFFSETS(IBEAM,UM) = 
     &           SZA_OFFSETS(IBEAM) + LOCAL_N_USERAZM * (UM - 1)
        END DO
      END DO

C  Fourier loop
C  ============

      LOCAL_ITERATION   = .TRUE.
      FOURIER_COMPONENT = -1
      TESTCONV          = 0

C  set up solar beam flags. Required in all cases.
C   ---Even required for the thermal.....

      DO IBEAM = 1, NSOURCES
        BEAM_TESTCONV  ( IBEAM )  = 0
        BEAM_ITERATION ( IBEAM ) = .TRUE.
        DO L = 0, MAXFOURIER
          DO_MULTIBEAM   ( IBEAM, L ) = .TRUE.
        ENDDO
      ENDDO

C  start loop

      DO WHILE ( LOCAL_ITERATION .AND.
     &              FOURIER_COMPONENT.LT.N_FOURIER_COMPONENTS )

C  Fourier counter

        FOURIER_COMPONENT = FOURIER_COMPONENT + 1

C  Local start of user-defined streams. Should always be 1.
C    No zenith tolerance now.

        LOCAL_UM_START = 1

C  azimuth cosine/sine factors, using adjust geometries.

        IF ( FOURIER_COMPONENT .GT. 0 ) THEN
          DFC = DBLE(FOURIER_COMPONENT)
          DO UA = 1, LOCAL_N_USERAZM
            DO IB = 1, NSOURCES
              DO UM = 1, N_USER_STREAMS
                AZM_ARGUMENT = USER_RELAZMS_ADJUST(UM,IB,UA) * DFC
                AZMFAC(UM,IB,UA,1)   = DCOS(DEG_TO_RAD*AZM_ARGUMENT)
                AZMFAC(UM,IB,UA,2)   = AZMFAC(UM,IB,UA,1)
                AZMFAC(UM,IB,UA,3)   = DSIN(DEG_TO_RAD*AZM_ARGUMENT)
                AZMFAC(UM,IB,UA,4)   = AZMFAC(UM,IB,UA,3)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

C  Main call to VLidort Fourier module

c        write(*,*)' ..calculating fourier component',FOURIER_COMPONENT

        CALL VLIDORT_FOURIER_LAMBERTIAN
     I        ( FOURIER_COMPONENT,
     I          SS_FLUX_MULTIPLIER,
     O          DO_INCLUDE_SURFACE,
     O          DO_INCLUDE_SURFEMISS,
     O          DO_INCLUDE_THERMEMISS,
     O          STATUS_SUB )

C  error handling

        IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
          STATUS_CALCULATION = VLIDORT_SERIOUS
          WRITE(CF,'(I3)')FOURIER_COMPONENT
          MAIL  = 
     &  'Error from VLIDORT_FOURIER_LAMBERTIAN, Fourier = '//CF
          TRACE = 'Called by VLIDORT_MASTER_LAMBERTIAN'
          CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS_CALCULATION )
          CLOSE(VLIDORT_ERRUNIT)
          RETURN
        ENDIF

C  Fourier summation and Convergence examination
C  ---------------------------------------------

C   -- only done for beams which are still not converged
C      This is controlled by flag DO_MULTIBEAM

C   -- new criterion, SS is added for Fourier = 0, as this means that
C      higher-order terms will be relatively smaller, which implies
C      faster convergence in some circumstances (generally not often).

        IBEAM_COUNT = 0
        DO IBEAM = 1, NSOURCES
          IF ( DO_MULTIBEAM ( IBEAM, FOURIER_COMPONENT ) ) THEN

C  Convergence and radiance summation

            CALL VLIDORT_CONVERGE
     I      ( AZMFAC, FOURIER_COMPONENT, LOCAL_N_USERAZM,
     O        BEAM_TESTCONV(IBEAM), IBEAM, BEAM_ITERATION(IBEAM) )

            IF ( BEAM_ITERATION(IBEAM) ) THEN
              IBEAM_COUNT = IBEAM_COUNT + 1
            ELSE
              DO L = FOURIER_COMPONENT+1,MAXFOURIER
                DO_MULTIBEAM (IBEAM,L) = .FALSE.
              ENDDO
            ENDIF

          ENDIF
        END DO

C  If all beams have converged, stop iteration

        IF ( IBEAM_COUNT .EQ. 0 ) LOCAL_ITERATION = .FALSE.

C  Fourier output
C  --------------

C  Open file if Fourier = 0
C  Write Standard Fourier output
C  Close file if iteration has finished

C  New comment:
C    If the SS correction is set, Fourier=0 will include SS field

        IF ( DO_WRITE_FOURIER ) THEN
          FUNIT = VLIDORT_FUNIT
          IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
            OPEN(FUNIT,FILE=FOURIER_WRITE_FILENAME,STATUS='UNKNOWN')
          ENDIF
          CALL VLIDORT_WRITEFOURIER ( FUNIT, FOURIER_COMPONENT )
          IF ( .NOT.LOCAL_ITERATION ) CLOSE ( FUNIT )
        ENDIF

C  end iteration loop

      ENDDO

C  restore no azimuth flag

      DO_NO_AZIMUTH = SAVE_DO_NO_AZIMUTH

C  Major result output
C  ===================

      IF ( DO_WRITE_RESULTS ) THEN
        RUNIT = VLIDORT_RESUNIT
        OPEN(RUNIT,FILE=RESULTS_WRITE_FILENAME,STATUS='UNKNOWN')
        CALL VLIDORT_WRITERESULTS ( RUNIT )
        CLOSE(RUNIT)
      ENDIF

C  Geophysical input (scenario) write
C  ----------------------------------

      IF ( DO_WRITE_SCENARIO ) THEN
        SUNIT = VLIDORT_SCENUNIT
        OPEN(SUNIT,FILE=SCENARIO_WRITE_FILENAME,STATUS='UNKNOWN')
        CALL VLIDORT_WRITESCEN ( SUNIT )
        CLOSE(SUNIT)
      ENDIF

C  close Error file if it was used

      IF ( .NOT. VLIDORT_ERROR_INIT ) CLOSE(VLIDORT_ERRUNIT)

C  Finish

      RETURN
      END

c

      SUBROUTINE VLIDORT_FOURIER_LAMBERTIAN
     I       ( FOURIER_COMPONENT,
     I         SS_FLUX_MULTIPLIER,
     O         DO_INCLUDE_SURFACE,
     O         DO_INCLUDE_SURFEMISS,
     O         DO_INCLUDE_THERMEMISS,
     O         STATUS )

C  Complete Fourier component calculation for the Standard Code.
C    Lambertian version: No BRDFs

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables 

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  Include file of setup variables (output from SETUPS routine)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_THERMALSUP.VARS'

C  include file for debug

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'

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

      INTEGER          LAYER, IBEAM, IPARTIC, I, O1, N

C  local inclusion flags

      LOGICAL          DO_INCLUDE_MVOUTPUT
      LOGICAL          DO_INCLUDE_DIRECTBEAM

C  Flux multiplier and Fourier component numbers

      DOUBLE PRECISION FLUX_MULTIPLIER
      DOUBLE PRECISION DELTA_FACTOR
      DOUBLE PRECISION SURFACE_FACTOR

C  Lambertian surface:  external function calls

      EXTERNAL         BVP_LAMBERTIAN_SURFACE_HOM
      EXTERNAL         BVP_LAMBERTIAN_SURFACE_BEAM
      EXTERNAL         BOA_LAMBERTIAN_SOURCE

C  error tracing

      CHARACTER*(2)    C2
      CHARACTER*(70)   MAIL, TRACE
      INTEGER          STATUS_SUB
      LOGICAL          FAIL

C  ##############
C  initialization
C  ##############

C  module status

      STATUS = VLIDORT_SUCCESS

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

C  Albedo flag (for inclusion of some kind of reflecting boundary)

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
C   (DO_REFLECTED_DIRECTBEAM is a stored Bookkeeping variable)

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

      FLUX_MULTIPLIER   = SS_FLUX_MULTIPLIER * DELTA_FACTOR
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

C  30 January 2008
C  Telescoping setup is done in DERIVE_INPUTS (unlike LIDORT scalar code)

!!!!!!  WARNING CHECK ARGUMENTS on VLIDORT_MISCSETUPS

      IF ( FOURIER_COMPONENT .EQ. 0 ) THEN

        CALL VLIDORT_MISCSETUPS
c        CALL VLIDORT_MISCSETUPS (STATUS_SUB,MAIL)
c        IF ( STATUS_SUB .EQ. VLIDORT_WARNING ) THEN
c          TRACE= 'BVP Telescoping Warning from VLIDORT_MISCSETUPS'
c          STATUS = VLIDORT_WARNING
c          CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
c        ENDIF 
 
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

C  Single scatter correction (pre-calculation)
C  -------------------------------------------

C  Not required if no solar sources

      IF ( DO_SOLAR_SOURCES ) THEN

C  SSCORR_nadir:
C    Must be done after MISCSETUPS, as we need the multipliers and
C    transmittance factors for the SUN and LOS paths (SSCORR_NADIR)
C    Code added 6 May 2005. Replaces call in Master routine.

C  SSCORR_outgoing
C      Version 2.1. Added call to the new outgoing sphericity correction

       IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
        IF ( DO_USER_STREAMS ) THEN

C  regular (nadir view) SS correction

          IF ( DO_SSCORR_NADIR ) THEN
            CALL VLIDORT_SSCORR_NADIR ( SS_FLUX_MULTIPLIER )
          ENDIF

C  New outgoing sphericity correction

          IF ( DO_SSCORR_OUTGOING ) THEN
            CALL VLIDORT_SSCORR_OUTGOING
     &     ( SS_FLUX_MULTIPLIER, FAIL, MAIL )
            IF ( FAIL ) THEN
              TRACE= 'SS correction outgoing failed'
              STATUS = VLIDORT_SERIOUS
              CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
              RETURN
            ENDIF 
          ENDIF

C  End Fourier = 0 and User stream


        ENDIF
       ENDIF

C--------------------------------------------------- COMMENTED OUT ------
C  Exact direct beam BRDF (if flagged)
C   Must be done after MISCSETUPS, as we need the SUN transmittance factors
C   BRDF Kernels already known (called before this module)
C   Code added 6 May 2005.
c      IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
c        IF ( DO_USER_STREAMS ) THEN
c          IF ( DO_DBCORRECTION ) THEN
c            CALL VLIDORT_DBCORRECTION ( SS_FLUX_MULTIPLIER )
c          ENDIF
c        ENDIF
c      ENDIF
C  BRDF Fourier components - coefficient setup
C   (also computes emissivity if required)
c      IF ( .NOT. DO_LAMBERTIAN_SURFACE ) THEN
c        CALL VLIDORT_BRDF_FRTERMS
c     I    ( DO_INCLUDE_SURFACE,
c     I      DO_INCLUDE_SURFEMISS,
c     I      FOURIER_COMPONENT,
c     I      DELTA_FACTOR )
c      ENDIF
C--------------------------------------------------- COMMENTED OUT ------

C  Reflected Direct beam attenuation (Lambertian) + DB correction (if flagged)
C  Full single scatter calculation, return after completion

       IF ( DO_SSFULL ) THEN
         CALL VLIDORT_LAMBERTIAN_DBCORR (SS_FLUX_MULTIPLIER)
         RETURN
       ELSE
         IF ( DO_DBCORRECTION ) THEN
           CALL VLIDORT_LAMBERTIAN_DBCORR (SS_FLUX_MULTIPLIER)
         ENDIF
         CALL LAMBERTIAN_SURFACE_DIRECTBEAM
     I    ( FOURIER_COMPONENT, DELTA_FACTOR )
       ENDIF

C  End solar sources only clause for corrections

      ENDIF

C  ###################
C  Spherical functions
C  ###################

C  Get Pi Matrices for this Fourier component
C    Including all beam angles !!

      IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
        CALL VLIDORT_PIMATRIX_SETUP ( FOURIER_COMPONENT )
      ENDIF
  
C  ########################################
C  RT differential equation Eigensolutions
C  ########################################

C  Go to continuation point for thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) GO TO 8899

C  Start layer loop

      DO LAYER = 1, NLAYERS

C  Get Discrete ordinate solutions for this layer

        CALL VLIDORT_QHOM_SOLUTION
     I    ( LAYER, FOURIER_COMPONENT,
     O      STATUS_SUB )

C  .. error tracing

        IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
          MAIL = 'Error return from module VLIDORT_QHOM_SOLUTION'
          TRACE= 'Call #1 in VLIDORT_FOURIER_LAMBERTIAN'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ENDIF

C  Get Post-processing ("user") solutions for this layer

        IF  ( STERM_LAYERMASK_UP(LAYER) .OR.
     &        STERM_LAYERMASK_DN(LAYER) ) THEN
          IF ( DO_USER_STREAMS ) THEN
            CALL VLIDORT_UHOM_SOLUTION
     I           ( LAYER, FOURIER_COMPONENT )
          ENDIF
        ENDIF

C  end layer loop

      ENDDO

C  Prepare homogeneous solution multipliers

      CALL HMULT_MASTER

C  ############################################
C   boundary value problem - MATRIX PREPARATION
C  ############################################

C  standard case using compression of band matrices, etc..

      IF ( BVP_REGULAR_FLAG(FOURIER_COMPONENT) ) THEN

        CALL BVP_MATRIXSETUP_MASTER
     I     ( DO_INCLUDE_SURFACE,
     I       FOURIER_COMPONENT,
     I       SURFACE_FACTOR,
     I       BVP_LAMBERTIAN_SURFACE_HOM,
     O       STATUS_SUB )

        IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
          MAIL = 'Error return from module BVP_MATRIXSETUP_MASTER'
          TRACE= 'Call #2 in VLIDORT_FOURIER_LAMBERTIAN'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ENDIF

C  Telescoped case

      ELSE

        CALL BVPTEL_MATRIXSETUP_MASTER
     I     ( DO_INCLUDE_SURFACE,
     I       FOURIER_COMPONENT, 
     I       SURFACE_FACTOR,
     I       BVP_LAMBERTIAN_SURFACE_HOM,
     O       STATUS_SUB )

        IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
          MAIL = 'Error return from module BVPTEL_MATRIXSETUP_MASTER'
          TRACE= 'Call #3 in VLIDORT_FOURIER_LAMBERTIAN'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ENDIF

      ENDIF

C  ################
C  Thermal Solution
C  ################

C  Continuation point for avoiding the scattering calculations

 8899 continue

C  1. Find the Particular solution (also for transmittance only)
C  2. Compute thermal layer source terms. (Upwelling and Downwelling)
C    These will be scaled up by factor 4.pi if solar beams as well
C  REMARK. No Green's function treatment here

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        CALL THERMAL_CLSOLUTION ( STATUS_SUB, MAIL )
        IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
          TRACE= 'Call #4 in VLIDORT_FOURIER_LAMBERTIAN'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ENDIF
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

      O1 = 1
      DO LAYER = 1, NLAYERS
        DO I = 1, NSTREAMS_2
          WUPPER(I,O1,LAYER) = T_WUPPER(I,LAYER)
          WLOWER(I,O1,LAYER) = T_WLOWER(I,LAYER)
       ENDDO
      ENDDO

C  Solve the boundary value problem. No telescoping here

      CALL BVP_SOLUTION_MASTER
     I       ( DO_INCLUDE_SURFACE,
     I         DO_INCLUDE_DIRECTBEAM,
     I         DO_INCLUDE_SURFEMISS,
     I         FOURIER_COMPONENT, IPARTIC,
     I         SURFACE_FACTOR,
     I         BVP_LAMBERTIAN_SURFACE_BEAM,
     O         STATUS_SUB )

C  Error handling

      IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
        MAIL = 'Error return from module BVP_SOLUTION_MASTER'
        TRACE= 'Call #5A in VLIDORT_FOURIER_LAMBERTIAN'
        STATUS = VLIDORT_SERIOUS
        CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
        RETURN
      ENDIF

C  Continuation point for avoiding thermal scattering

 566  CONTINUE

C  Post-processing - upwelling thermal-only field

      IF ( DO_UPWELLING ) THEN
        CALL VLIDORT_UPUSER_INTENSITY
     I        ( DO_INCLUDE_SURFACE,
     I          DO_INCLUDE_SURFEMISS,
     I          DO_INCLUDE_THERMEMISS,
     I          DO_INCLUDE_DIRECTBEAM,
     I          DO_INCLUDE_MVOUTPUT,
     I          FOURIER_COMPONENT,  IPARTIC,
     I          SURFACE_FACTOR,
     I          BOA_LAMBERTIAN_SOURCE,
     I          FLUX_MULTIPLIER )
      ENDIF

C  Post-processing - Downwelling thermal-only field

      IF ( DO_DNWELLING ) THEN
        CALL VLIDORT_DNUSER_INTENSITY
     I        ( DO_INCLUDE_THERMEMISS,
     I          FLUX_MULTIPLIER, 
     I          FOURIER_COMPONENT,
     I          IPARTIC )
      ENDIF

C  mean value (integrated) output for thermal-only field
C    ------- Can also use this to get debug Quadrature output

      IF ( DO_INCLUDE_MVOUTPUT.OR.DO_QUAD_OUTPUT ) THEN
         CALL VLIDORT_INTEGRATED_OUTPUT
     I      ( DO_INCLUDE_MVOUTPUT,
     I        DO_INCLUDE_DIRECTBEAM,
     I        DO_INCLUDE_THERMEMISS,
     I        FLUX_MULTIPLIER, IPARTIC )
      ENDIF

C  Finish Thermal only.

      RETURN

C  ##################################################
C  Complete Radiation Field with Solar Beam solutions
C  ##################################################

C  Continuation point

 455  CONTINUE

C  Start loop over desired beam solutions

      DO IBEAM = 1, NBEAMS

        IF ( DO_MULTIBEAM(IBEAM,FOURIER_COMPONENT) ) THEN

C  Solar beam Particular solutions (classical method)
C  --------------------------------------------------

          IF ( DO_CLASSICAL_SOLUTION ) THEN

            DO LAYER = 1, NLAYERS

C  get the solution

              CALL  VLIDORT_QBEAM_SOLUTION
     I           ( LAYER, FOURIER_COMPONENT, IBEAM, STATUS_SUB )

C  error flag

              IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
                WRITE(C2,'(I2)')LAYER
                MAIL = 'Error from VLIDORT_QBEAM_SOLUTION, layer'//C2
                TRACE= 'Call #4 in VLIDORT_FOURIER_LAMBERTIAN'
                STATUS = VLIDORT_SERIOUS
                CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
                RETURN
              ENDIF

C  user solutions

              IF  ( STERM_LAYERMASK_UP(LAYER) .OR.
     &              STERM_LAYERMASK_DN(LAYER) ) THEN
                IF ( DO_USER_STREAMS ) THEN
                  CALL VLIDORT_UBEAM_SOLUTION
     I             ( LAYER, FOURIER_COMPONENT, IBEAM )
                END IF
              END IF

C  end layer loop

            END DO

C  Solar beam Particular solutions (Green's function)
C  -------------------------------------------------

          ELSE
C  $$$$$$$$$$$$$$$$$$$$$$    P L A C E H O L D E R    $$$$$$$$$$$$$$$$$$$$
C            CALL V_GREENFUNC_SOLUTION ( FOURIER)
C  $$$$$$$$$$$$$$$$$$$$$$    P L A C E H O L D E R    $$$$$$$$$$$$$$$$$$$$
          ENDIF

C  Add thermal solutions if flagged
C    NO 4.PI modulus on the thermal contribution.

          IF ( DO_INCLUDE_THERMEMISS ) THEN
            O1 = 1
            DO N = 1, NLAYERS
              DO I = 1, NSTREAMS_2
                WUPPER(I,O1,N) = WUPPER(I,O1,N) + pi4*T_WUPPER(I,N)
                WLOWER(I,O1,N) = WLOWER(I,O1,N) + pi4*T_WLOWER(I,N)
              ENDDO
            ENDDO
          ENDIF

C  Solve boundary value problem
C  ----------------------------

C  standard case using compressed-band matrices, etc..

          IF ( BVP_REGULAR_FLAG(FOURIER_COMPONENT) ) THEN
            
C  Get the BVP solution (regular case) for this beam component

            CALL BVP_SOLUTION_MASTER
     I       ( DO_INCLUDE_SURFACE,
     I         DO_REFLECTED_DIRECTBEAM(IBEAM),
     I         DO_INCLUDE_SURFEMISS,
     I         FOURIER_COMPONENT, IBEAM,
     I         SURFACE_FACTOR,
     I         BVP_LAMBERTIAN_SURFACE_BEAM,
     O         STATUS_SUB )

            IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
              MAIL = 'Error return from module BVP_SOLUTION_MASTER'
              TRACE= 'Call #5 in VLIDORT_FOURIER_LAMBERTIAN'
              STATUS = VLIDORT_SERIOUS
              CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
              RETURN
            ENDIF

C  Telescoped case

          ELSE

            CALL BVPTEL_SOLUTION_MASTER
     I         ( DO_INCLUDE_SURFACE,
     I           DO_REFLECTED_DIRECTBEAM(IBEAM),
     I           FOURIER_COMPONENT, IBEAM,
     I           SURFACE_FACTOR,
     I           BVP_LAMBERTIAN_SURFACE_BEAM,
     O           STATUS_SUB )

            IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
              MAIL = 'Error return from module BVPTEL_SOLUTION_MASTER'
              TRACE= 'Call #6 in VLIDORT_FOURIER_LAMBERTIAN'
              STATUS = VLIDORT_SERIOUS
              CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
              RETURN
            ENDIF

          ENDIF

C  Stokes Vector Post Processing
C  -----------------------------

C  upwelling

          IF ( DO_UPWELLING ) THEN

C LIDORT COMMENTS:  Direct beam inclusion flag:
C   This now has the DBCORRECTION option: if the DBCORRECTION Flag
C   is set, then we will be doing exact calculations of the reflected
C   directbeam, so we do not need to include it in the Post-processing.
C   However, the direct beam will need to be included in the basic RT
C   solution (the BVP), and this is controlled separately by the
C   DO_REFLECTED_DIRECTBEAM(IBEAM) flags.
C     R. Spurr, RT Solutions, Inc., 19 August 2005.

            DO_INCLUDE_DIRECTBEAM = ( DO_UPWELLING .AND.
     &     (DO_REFLECTED_DIRECTBEAM(IBEAM).AND..NOT.DO_DBCORRECTION))

            CALL VLIDORT_UPUSER_INTENSITY
     I        ( DO_INCLUDE_SURFACE,
     I          DO_INCLUDE_SURFEMISS,
     I          DO_INCLUDE_THERMEMISS,
     I          DO_INCLUDE_DIRECTBEAM,
     I          DO_INCLUDE_MVOUTPUT,
     I          FOURIER_COMPONENT,  IBEAM,
     I          SURFACE_FACTOR,
     I          BOA_LAMBERTIAN_SOURCE,
     I          FLUX_MULTIPLIER )
          ENDIF

C  Downwelling

          IF ( DO_DNWELLING ) THEN
            CALL VLIDORT_DNUSER_INTENSITY
     I        ( DO_INCLUDE_THERMEMISS,
     I          FLUX_MULTIPLIER, 
     I          FOURIER_COMPONENT, IBEAM )
          ENDIF

C  mean value (integrated) output.
C    ------- Can also use this to get debug Quadrature output

          IF ( DO_INCLUDE_MVOUTPUT.OR.DO_QUAD_OUTPUT ) THEN
            CALL VLIDORT_INTEGRATED_OUTPUT
     I      ( DO_INCLUDE_MVOUTPUT,
     I        DO_REFLECTED_DIRECTBEAM(IBEAM),
     I        DO_INCLUDE_THERMEMISS,
     I        FLUX_MULTIPLIER, IBEAM )
          ENDIF

C  End loop over beam solutions

        END IF
      END DO

c  Debug write

      IF ( DO_DEBUG_WRITE ) THEN
        write(*,'(a,I3,a,100I3)')
     &     'Fourier ',FOURIER_COMPONENT,
     &   ' : # complex evalues by layer: ',
     &       (K_COMPLEX(LAYER),LAYER=1,NLAYERS)
      ENDIF

C  ######
C  finish
C  ######

      RETURN
      END

