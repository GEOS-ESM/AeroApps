C ###########################################################
C #                                                         #
C #             THE TWOSTREAM LIDORT MODEL                  #
C #                                                         #
C #      (LInearized Discrete Ordinate Radiative Transfer)  #
C #       --         -        -        -         -          #
C #                                                         #
C ###########################################################

C ###########################################################
C #                                                         #
C #  Authors :      Robert. J. D. Spurr (1)                 #
C #                 Vijay Natraj        (2)                 #
C #                                                         #
C #  Address (1) :     RT Solutions, Inc.                   #
C #                    9 Channing Street                    #
C #                    Cambridge, MA 02138, USA             #
C #  Tel:             (617) 492 1183                        #
C #  Email :           rtsolutions@verizon.net              #
C #                                                         #
C #  Address (2) :     CalTech                              #
C #                    Department of Planetary Sciences     #
C #                    1200 East California Boulevard       #
C #                    Pasadena, CA 91125                   #
C #  Tel:             (626) 395 6962                        #
C #  Email :           vijay@gps.caltech.edu                #
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
C #            TWOSTREAM_CHECK_INPUTS_BASIC                     #
C #            TWOSTREAM_CHECK_INPUTS_THREAD                    #
C #                                                             #
C #            TWOSTREAM_ERROR_TRACE                            #
C #                                                             #
C ###############################################################

C

      SUBROUTINE TWOSTREAM_CHECK_INPUTS_BASIC
     I   ( DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL,
     I     NLAYERS, NFINELAYERS, NBEAMS, NMOMENTS_INPUT,
     I     N_USER_STREAMS, N_USER_RELAZMS,
     I     BEAM_SZAS, USER_ANGLES, USER_RELAZMS,
     I     EARTH_RADIUS, HEIGHT_GRID, INIT,
     O     STATUS )

C  Inputs
C  ------

C  Flags

      LOGICAL          DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL

C  Single scatter flags omitted from this streamlined version
C      LOGICAL          DO_SSCORR_OUTGOING, DO_SSCORR_NADIR, DO_SSFULL

C  Numbers

      INTEGER          NLAYERS, NFINELAYERS, NMOMENTS_INPUT
      INTEGER          NBEAMS, N_USER_STREAMS, N_USER_RELAZMS  

C  Geometry

      DOUBLE PRECISION BEAM_SZAS    ( NBEAMS )
      DOUBLE PRECISION USER_ANGLES  ( N_USER_STREAMS )
      DOUBLE PRECISION USER_RELAZMS ( N_USER_RELAZMS )

C  height and earth radius

      DOUBLE PRECISION EARTH_RADIUS
      DOUBLE PRECISION HEIGHT_GRID ( 0:NLAYERS )

C  Error log file initialization flag

      LOGICAL          INIT

C  Module output

      INTEGER          STATUS

C  local variables

      CHARACTER*70     MAIL
      CHARACTER*70     ACTION

      INTEGER          I
      CHARACTER*2      C2
      LOGICAL          LOOP

C  Initialize output status

      STATUS = 0

C  Single scattering stuff omitted
c      IF ( DO_SSCORR_OUTGOING ) THEN
c        IF ( NFINELAYERS .GT. 10 ) THEN
c          MAIL   = 'Number of fine layers > 10'
c          ACTION = 'Re-set input value'
c          STATUS = 1
c          CALL TWOSTREAM_ERROR_TRACE ( INIT, MAIL, ACTION, STATUS )
c          RETURN
c        ENDIF
c      ENDIF

      IF ( NMOMENTS_INPUT .GT. 300 ) THEN
        MAIL   = 'Number of Legendre moments > 300'
        ACTION = 'Re-set input value'
c        STATUS = 1
c        CALL TWOSTREAM_ERROR_TRACE ( INIT, MAIL, ACTION, STATUS )
        RETURN
      ENDIF

C  3. Basic geometry

      IF ( NBEAMS .GT. 10 ) THEN
        MAIL   = 'Number of solar zenith angles > 10'
        ACTION = 'Re-set input value'
        STATUS = 1
        CALL TWOSTREAM_ERROR_TRACE ( INIT, MAIL, ACTION, STATUS )
        RETURN
      ENDIF

      IF ( N_USER_STREAMS .GT. 10 ) THEN
        MAIL   = 'Number of user streams > 10'
        ACTION = 'Re-set input value'
        STATUS = 1
        CALL TWOSTREAM_ERROR_TRACE ( INIT, MAIL, ACTION, STATUS )
        RETURN
      ENDIF

      IF ( N_USER_RELAZMS .GT. 10 ) THEN
        MAIL   = 'Number of relative azimuths > 10'
        ACTION = 'Re-set input value'
        STATUS = 1
        CALL TWOSTREAM_ERROR_TRACE ( INIT, MAIL, ACTION, STATUS )
        RETURN
      ENDIF

C  check directional input

      IF ( .NOT.DO_UPWELLING .AND. .NOT. DO_DNWELLING ) THEN
        MAIL   = 'Bad input: no directional input is set'
        ACTION = 'Check DO_UPWELLING & DO_DNWELLING: one must be set!'
        STATUS = 1
        CALL TWOSTREAM_ERROR_TRACE ( INIT, MAIL, ACTION, STATUS )
      ENDIF

C  check SS inputs (both file-read and derived)
C  --------------------------------------------

C  This section commented out in the streamlined version

C  Check sphericity corrections....Cannot both be turned on
c      IF ( DO_SSCORR_NADIR .and. DO_SSCORR_OUTGOING ) THEN
c        MAIL   = 'Cannot have both single scatter corrections on'
c        ACTION = 'Turn off DO_SSCORR_NADIR and/or DO_SSCORR_OUTGOING'
c        STATUS = 1
c        CALL TWOSTREAM_ERROR_TRACE ( INIT, MAIL, ACTION, STATUS )
c      ENDIF     

C   Single scatter corrections must be turned on
c      IF ( DO_SSFULL ) THEN
c        IF ( .not.DO_SSCORR_NADIR.and..not.DO_SSCORR_OUTGOING ) THEN        
c          MAIL   = 'Bad input: Full SS, must have one SSCORR flag set'
c          ACTION = 'Full SS: default to use outgoing SS correction'
c          DO_SSCORR_NADIR    = .FALSE.
c          DO_SSCORR_OUTGOING = .TRUE.
c          CALL TWOSTREAM_ERROR_TRACE ( INIT, MAIL, ACTION, STATUS )
c        ENDIF
c      ENDIF

C  check viewing geometry input
C  ============================

C  Check earth radius (Chapman function only)
C    ---WARNING. Default value of 6371.0 will be set

      IF ( .NOT. DO_PLANE_PARALLEL ) THEN
        IF ( EARTH_RADIUS.LT.6320.0D0 .OR.
     &       EARTH_RADIUS.GT.6420.0D0 ) THEN
          MAIL   = 'Bad input: Earth radius outside of [6320-6420]'
          ACTION = 'Warning: default value of 6371.0 was set'
          EARTH_RADIUS = 6371.0D0
          CALL TWOSTREAM_ERROR_TRACE ( INIT, MAIL, ACTION, STATUS )
        ENDIF
      ENDIF

C  Check solar zenith angle input

      DO I = 1, NBEAMS
        IF ( BEAM_SZAS(I) .LT. 0.0d0 .OR.
     &       BEAM_SZAS(I).GE.90.0D0 ) THEN
          WRITE(C2,'(I2)')I
          MAIL   = 'Bad input: out-of-range beam angle, no. '//C2
          ACTION = 'Look at BEAM_SZAS input, should be < 90 & > 0'
          STATUS = 1
          CALL TWOSTREAM_ERROR_TRACE ( INIT, MAIL, ACTION, STATUS )
        ENDIF
      ENDDO

C  Check relative azimuths

      LOOP = .TRUE.
      I = 0
      DO WHILE (LOOP .AND. I.LT.N_USER_RELAZMS)
        I = I + 1
        IF ( USER_RELAZMS(I) .GT. 360.0D0   .OR.
     &       USER_RELAZMS(I) .LT. 0.0d0 ) THEN
          WRITE(C2,'(I2)')I
          MAIL = 'Bad input: out-of-range azimuth angle, no. '//C2
          ACTION = 'Look at azimuth angle input, should be in [0,360]'
          LOOP = .FALSE.
          STATUS = 1
          CALL TWOSTREAM_ERROR_TRACE ( INIT, MAIL, ACTION, STATUS )
        ENDIF
      ENDDO

C  check user-defined stream angles (should always be [0,90])

      LOOP = .TRUE.
      I = 0
      DO WHILE (LOOP .AND. I.LT.N_USER_STREAMS)
        I = I + 1
        IF ( USER_ANGLES(I) .GT. 90.0   .OR.
     &       USER_ANGLES(I) .LT. 0.0d0 ) THEN
          WRITE(C2,'(I2)')I
          MAIL   = 'Bad input: out-of-range user stream, no. '//C2
          ACTION = 'Look at user-defined angle input'
          LOOP = .FALSE.
          STATUS = 1
          CALL TWOSTREAM_ERROR_TRACE ( INIT, MAIL, ACTION, STATUS )
        ENDIF
      ENDDO

C  Check height grid input (Chapman function only)

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
          STATUS = 1
          CALL TWOSTREAM_ERROR_TRACE ( INIT, MAIL, ACTION, STATUS )
        ENDIF
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE TWOSTREAM_CHECK_INPUTS_THREAD
     I     ( INIT, NLAYERS, NTHREADS, THREAD, 
     I       SURFTYPE, LAMBERTIAN_ALBEDO,
     I       DELTAU_VERT, OMEGA_TOTAL, ASYMM_TOTAL,
     O       STATUS )

C  Module inputs
C  -------------

C  Error intial8ziation flag

      LOGICAL          INIT

C  Numbers

      INTEGER          NLAYERS, NTHREADS

C  Thread number

      INTEGER          THREAD

C  surface type

      INTEGER          SURFTYPE

C  Optical properties

      DOUBLE PRECISION DELTAU_VERT(NLAYERS, NTHREADS)
      DOUBLE PRECISION OMEGA_TOTAL(NLAYERS, NTHREADS)
      DOUBLE PRECISION ASYMM_TOTAL(NLAYERS, NTHREADS)
      DOUBLE PRECISION LAMBERTIAN_ALBEDO(NTHREADS)

C  Module output
C  -------------

      INTEGER          STATUS

C  local variables

      CHARACTER*70     MAIL
      CHARACTER*70     ACTION
      INTEGER          L
      CHARACTER*3      C3, WTHREAD

C  Initialize output status

      STATUS = 0

C  thread number

      WTHREAD = '000'
      IF (THREAD.LT.10)WRITE(WTHREAD(3:3),'(I1)')THREAD
      IF (THREAD.GT.99)WRITE(WTHREAD(1:3),'(I3)')THREAD
      IF (THREAD.GE.10.and.THREAD.LE.99)WRITE(WTHREAD(2:3),'(I2)')THREAD

C  check Thread-dependent optical property inputs
C  ----------------------------------------------

C  make sure the Lambertian surface is in range

      IF ( SURFTYPE .EQ. 1 ) THEN
       IF ( LAMBERTIAN_ALBEDO(THREAD) .LT.0.0d0 .OR.
     &     LAMBERTIAN_ALBEDO(THREAD) .GT.1.0d0 ) THEN
        MAIL = 'Bad input: Lambertian albedo not in range [0,1]'
        ACTION = 'Check albedo input, thread # '//wthread
        STATUS = 1
        CALL TWOSTREAM_ERROR_TRACE ( INIT, MAIL, ACTION, STATUS )
       ENDIF
      ENDIF

C  Check non-negative optical thickness values

      DO L = 1, NLAYERS
        IF ( DELTAU_VERT(L,THREAD).LE.0.0d0 ) THEN
          WRITE(C3,'(I3)')L
          MAIL = 'Bad input: optical thickness <= 0, layer '//C3
          ACTION = 'Check optical thickness input, thread # '//wthread
          STATUS = 1
          CALL TWOSTREAM_ERROR_TRACE ( INIT, MAIL, ACTION, STATUS )
        ENDIF
      ENDDO

C  check single scatter albedos

      DO L = 1, NLAYERS
        IF ( OMEGA_TOTAL(L,THREAD).GT.0.999999d0 ) THEN
          WRITE(C3,'(I3)')L
          MAIL = 'Bad input: SS-albedo too close to 1, layer '//C3
          ACTION = 'Check SS-albedo input, thread # '//wthread
          STATUS = 1
          CALL TWOSTREAM_ERROR_TRACE ( INIT, MAIL, ACTION, STATUS )
        ENDIF
      ENDDO

C  solar beam, cannot be too small

      DO L = 1, NLAYERS
        IF ( OMEGA_TOTAL(L,THREAD).LT.1.0d-06 ) THEN
          WRITE(C3,'(I3)')L
          MAIL = 'Bad input: SS-albedo too close to 0, layer '//C3
          ACTION = 'Check SS-albedo input, thread # '//wthread
          STATUS = 1
          CALL TWOSTREAM_ERROR_TRACE ( INIT, MAIL, ACTION, STATUS )
        ENDIF
      ENDDO
      
C  Finish

      RETURN
      END

C

      SUBROUTINE TWOSTREAM_ERROR_TRACE ( INIT, MESSAGE, TRACE, STATUS )

C  subroutine arguments

      LOGICAL             INIT
      INTEGER             STATUS
      CHARACTER*(*)       TRACE, MESSAGE

C  first call - open file

      IF ( INIT ) THEN
        OPEN ( 25, FILE = '2stream_errors.log', STATUS = 'UNKNOWN' )
      ENDIF

C  write messages

      IF ( STATUS .NE. 0 ) THEN
        IF ( TRACE .EQ. ' ' ) THEN
          WRITE(25,'(/a)')MESSAGE
        ELSE
          WRITE(25,'(/a)')MESSAGE
          WRITE(25,'(a)')TRACE
        ENDIF
      ENDIF

C  turn off error initialization flag

      INIT = .FALSE.

C  finish

      RETURN
      END
