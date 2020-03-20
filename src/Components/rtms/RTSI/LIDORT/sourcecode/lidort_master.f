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

      SUBROUTINE LIDORT_MASTER
     &     ( STATUS_INPUTCHECK, STATUS_CALCULATION )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file for results (outputs)

      INCLUDE '../includes/LIDORT_RESULTS.VARS'

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

      INTEGER          UA, UM, IB, TESTCONV, L
      INTEGER          LOCAL_N_USERAZM, STATUS_SUB
      INTEGER          IUNIT, SUNIT, FUNIT, RUNIT

      DOUBLE PRECISION AZM_ARGUMENT, DFC
      DOUBLE PRECISION AZMFAC
     &     (MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)

      INTEGER          IBEAM_COUNT, IBEAM, NSOURCES
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

      STATUS_CALCULATION = LIDORT_SUCCESS
      STATUS_INPUTCHECK  = LIDORT_SUCCESS

C  Single scatter correction: flux multiplier
C    Now always F / 4pi

      SS_FLUX_MULTIPLIER = FLUX_FACTOR / PI4
      
C  Check input
C  -----------

      CALL LIDORT_CHECK_INPUT ( STATUS_SUB )

      IF ( STATUS_SUB .EQ. LIDORT_SERIOUS ) THEN
        STATUS_INPUTCHECK = LIDORT_SERIOUS
        MAIL = ' LIDORT_CHECK_INPUT failed'
        TRACE = ' Called in LIDORT_MASTER '
        CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS_INPUTCHECK )
        CLOSE(LIDORT_ERRUNIT)
        RETURN
      ELSE IF ( STATUS_SUB .EQ. LIDORT_WARNING ) THEN
        STATUS_INPUTCHECK = LIDORT_WARNING
        MAIL = ' LIDORT_CHECK_INPUT warning message'
        TRACE = ' Called in LIDORT_MASTER '
        CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS_INPUTCHECK )
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
     i   ( MAX_USER_STREAMS, MAXBEAMS, MAX_USER_RELAZMS,
     i     N_USER_STREAMS,   NBEAMS,   N_USER_RELAZMS,
     i     height_grid(nlayers), modified_eradius, adjust_surface,
     i     user_angles_input,  beam_szas, user_relazms,
     o     user_angles_adjust, beam_szas_adjust, user_relazms_adjust,
     o     fail, mail )
      if ( fail ) return       

C  Chapman function calculation
C  ----------------------------

      IF ( DO_SOLAR_SOURCES ) THEN
       IF ( DO_CHAPMAN_FUNCTION ) THEN
        CALL LIDORT_CHAPMAN ( FAIL, MAIL, TRACE )
        IF (FAIL) THEN
          STATUS_CALCULATION = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS_CALCULATION )
          CLOSE(LIDORT_ERRUNIT)
          RETURN
        ENDIF
       ENDIF
      ENDIF

C  write input variables
C  ---------------------

C  open file, call standard input write, close file

      IF ( DO_WRITE_INPUT ) THEN
        IUNIT = LIDORT_INUNIT
        OPEN(IUNIT,FILE=INPUT_WRITE_FILENAME,STATUS='UNKNOWN')
        CALL LIDORT_WRITEINPUT ( IUNIT )
        CLOSE(IUNIT)
      ENDIF

C  Get derived inputs
C  ==================

C  Miscellaneous and layer input

      CALL LIDORT_DERIVE_INPUT

C  BRDF input functions (non-Lambertian)

      IF ( .NOT. DO_LAMBERTIAN_SURFACE ) THEN
        CALL LIDORT_BRDF_MASTER
      ENDIF

C  Initialise Fourier loop
C  =======================

C  Set Number of Fourier terms (NMOMENTS = Maximum).
C    ( Starting from 0 = Fundamental )

      SAVE_DO_NO_AZIMUTH  = DO_NO_AZIMUTH
      LOCAL_DO_NO_AZIMUTH = DO_NO_AZIMUTH

C  No azimuth dependency for following four cases

C      IF ( DO_TRANSMITTANCE_ONLY ) THEN
C        LOCAL_DO_NO_AZIMUTH = .TRUE.
C      ENDIF

      IF ( .NOT. DO_SOLAR_SOURCES  ) THEN
        LOCAL_DO_NO_AZIMUTH = .TRUE.
      ENDIF

      IF ( DO_ISOTROPIC_ONLY  ) THEN
        LOCAL_DO_NO_AZIMUTH = .TRUE.
      ENDIF

      IF ( DO_MVOUT_ONLY  ) THEN
        LOCAL_DO_NO_AZIMUTH = .TRUE.
      ENDIF

C  set Fourier number (2 for Rayleigh only, otherwise 2*Nstreams))

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
        NSOURCES = NBEAMS
      ELSE
        NSOURCES = 1
      ENDIF

C  save some offsets for indexing geometries
C     These results are now stored in LIDORT_BOOKKEEP.VARS

      N_VIEWING    = N_OUT_STREAMS * LOCAL_N_USERAZM
      N_GEOMETRIES = NSOURCES * N_VIEWING

      DO IBEAM = 1, NBEAMS
        IBOFF(IBEAM) = N_VIEWING * ( IBEAM - 1 )
        DO UM = 1, N_OUT_STREAMS
          UMOFF(IBEAM,UM) = IBOFF(IBEAM) + LOCAL_N_USERAZM * (UM - 1)
        END DO
      END DO

C  Fourier loop
C  ============

C  initialize

      LOCAL_ITERATION   = .TRUE.
      FOURIER_COMPONENT = -1
      TESTCONV          = 0

C  set up solar beam flags

      DO IBEAM = 1, NBEAMS
        BEAM_TESTCONV  ( IBEAM )  = 0
        BEAM_ITERATION ( IBEAM ) = .TRUE.
        DO L = 0, MAXFOURIER
          DO_MULTIBEAM   ( IBEAM, L ) = .TRUE.
        ENDDO
      ENDDO

C  start loop
C  ----------

      DO WHILE ( LOCAL_ITERATION .AND.
     &              FOURIER_COMPONENT.LT.N_FOURIER_COMPONENTS )

C  Fourier counter

        FOURIER_COMPONENT = FOURIER_COMPONENT + 1

C  Local start of user-defined streams
C  Now set = 1. Fudging in earlier versions caused Havoc.
C  Here is the older version fudge
C        IF ( FOURIER_COMPONENT. EQ. 0 ) THEN
C          LOCAL_UM_START = 1
C        ELSE
C          LOCAL_UM_START = N_OUT_STREAMS - N_CONV_STREAMS + 1
C        ENDIF

        LOCAL_UM_START = 1

C  azimuth cosine factor, using adjust geometries.

        IF ( FOURIER_COMPONENT .GT. 0 ) THEN
          DFC = DBLE(FOURIER_COMPONENT)
          DO UA = 1, LOCAL_N_USERAZM
            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                AZM_ARGUMENT = USER_RELAZMS_ADJUST(UM,IB,UA) * DFC
                AZMFAC(UM,IB,UA)   = DCOS(DEG_TO_RAD*AZM_ARGUMENT)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

C  Main call to Lidort Fourier module
C  ----------------------------------

C        write(*,*)' ..calculating fourier component',FOURIER_COMPONENT

        CALL LIDORT_FOURIER_MASTER
     I        ( FOURIER_COMPONENT,
     I          SS_FLUX_MULTIPLIER,
     O          DO_INCLUDE_SURFACE,
     O          DO_INCLUDE_SURFEMISS,
     O          DO_INCLUDE_THERMEMISS,
     O          STATUS_SUB )

C  error handling

        IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
          STATUS_CALCULATION = LIDORT_SERIOUS
          WRITE(CF,'(I3)')FOURIER_COMPONENT
          MAIL = 'Error from LIDORT_FOURIER_MASTER, Fourier = '//CF
          TRACE = 'Called by LIDORT_MASTER'
          CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS_CALCULATION )
          CLOSE(LIDORT_ERRUNIT)
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
        DO IBEAM = 1, NBEAMS
          IF ( DO_MULTIBEAM ( IBEAM, FOURIER_COMPONENT ) ) THEN

C  Convergence and radiance summation

            CALL LIDORT_CONVERGE
     I      ( AZMFAC, FOURIER_COMPONENT, LOCAL_N_USERAZM,
     O        BEAM_TESTCONV(IBEAM), IBEAM, BEAM_ITERATION(IBEAM) )

C  Check number of beams already converged

            IF ( BEAM_ITERATION(IBEAM) ) THEN
              IBEAM_COUNT = IBEAM_COUNT + 1
            ELSE
              DO L = FOURIER_COMPONENT+1,MAXFOURIER
                DO_MULTIBEAM (IBEAM,L) = .FALSE.
              ENDDO
            ENDIF

C  end beam count loop

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
          FUNIT = LIDORT_FUNIT
          IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
            OPEN(FUNIT,FILE=FOURIER_WRITE_FILENAME,STATUS='UNKNOWN')
          ENDIF
          CALL LIDORT_WRITEFOURIER ( FUNIT, FOURIER_COMPONENT )
          IF ( .NOT.LOCAL_ITERATION ) CLOSE ( FUNIT )
        ENDIF

C  end iteration loop

      ENDDO

C  restore no azimuth flag

      DO_NO_AZIMUTH = SAVE_DO_NO_AZIMUTH

C  Major result output
C  -------------------

      IF ( DO_WRITE_RESULTS ) THEN
        RUNIT = LIDORT_RESUNIT
        OPEN(RUNIT,FILE=RESULTS_WRITE_FILENAME,STATUS='UNKNOWN')
        CALL LIDORT_WRITERESULTS ( RUNIT )
        CLOSE(RUNIT)
      ENDIF

C  Geophysical input (scenario) write
C  ----------------------------------

      IF ( DO_WRITE_SCENARIO ) THEN
        SUNIT = LIDORT_SCENUNIT
        OPEN(SUNIT,FILE=SCENARIO_WRITE_FILENAME,STATUS='UNKNOWN')
        CALL LIDORT_WRITESCEN ( SUNIT )
        CLOSE(SUNIT)
      ENDIF

C  close Error file if it was used

      IF ( .NOT. LIDORT_ERROR_INIT ) CLOSE(LIDORT_ERRUNIT)

C  Finish

      RETURN
      END

