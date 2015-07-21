! ###############################################################
! #                                                             #
! #                    THE VECTOR LIDORT MODEL                  #
! #                                                             #
! #  (Vector LInearized Discrete Ordinate Radiative Transfer)   #
! #   -      --         -        -        -         -           #
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
! #            VSLEAVE_LIN_INPUTMASTER                          #
! #            VSLEAVE_LIN_MAINMASTER (master)                  #
! #                                                             #
! ###############################################################

      MODULE vsleave_linsup_masters_m

      PRIVATE
      PUBLIC :: VSLEAVE_LIN_INPUTMASTER,&
                VSLEAVE_LIN_MAINMASTER

      CONTAINS

      SUBROUTINE VSLEAVE_LIN_INPUTMASTER ( &
        FILNAM, VSLEAVE_Sup_In, VSLEAVE_LinSup_In, &
        VSLEAVE_Sup_InputStatus )

!  Input routine for VSLEAVE program

!  Version 2.6 Notes
!  -----------------

!  INDWAT, MORCASIWAT Routines taken straight from Clark Weaver code
!      Compiled here by R. Spurr, 11 July 2012
!  get_fluorescence_755 Routines taken straight from Chris O'dell code
!      Compiled here by R. Spurr, 12 July 2012

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Upgrade for Version 2.7
!  -----------------------

!  New Water-leaving code developed/tested by R. Spurr, April 2014
!  Based in part on Modified-6S code by A. Sayer.
!  Validated against Modified-6 OCEABRDF.F code, 24  April 2014.

!  Water-leaving upgrade according to the modified 6S specification
!  New code calculates transmittances into and out of ocean, using
!  usual sun-glint rough surface approximations. In addition, the
!  water-leaving term itself is now SZA-dependent (A. Sayer), and
!  there is now a correction for Whitecaps (again, from 6S)

!  This upgrade gives the  water-leaving terms some directionality,
!  but they are still azimuth-independent

!  Earlier version inputs were just Wavelength/Salinity/PigmentConc
!  This was enough for the isotropic case (Fast Option) in Version 2.6.
!  For Version 2.7, we require additional inputs, including:
!    - Wind-speed and Direction (direction was not used in earlier version) 
!    - flags to control use sunglint shadowing and foam (whitecaps) correction.

!  This water-leaving option is designed to work alongside the "NewCM" glint
!  reflectance option in the VBRDF code. The glint and whitecap calculations
!  in the two supplements are the same.

!  You need to make sure that the wind input information is the same as that
!  used for the "NewCM" glint option in the VBRDF supplement. Also, the Foam
!  correction applied here in the surface-leaving code should also be 
!  applied for "NewCM" glint.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      USE VLIDORT_PARS
      USE VSLEAVE_FINDPAR_M

      USE vsleave_sup_inputs_def
      USE vsleave_linsup_inputs_def

      USE vsleave_sup_outputs_def

      IMPLICIT NONE

!  Arguments
!  ---------

      CHARACTER (LEN=*), INTENT(IN) :: FILNAM

      TYPE(VSLEAVE_Sup_inputs)   , INTENT(OUT) :: VSLEAVE_Sup_In
      TYPE(VSLEAVE_LinSup_inputs), INTENT(OUT) :: VSLEAVE_LinSup_In

      TYPE(VSLEAVE_Input_Exception_Handling), INTENT(OUT) :: &
        VSLEAVE_Sup_InputStatus

!  Local variables
!  ===============

!  Main Boolean flags
!  ------------------

!  Inclusion flag (not really necessary, Brian)

      LOGICAL :: DO_SLEAVING

!  Isotropic flag

      LOGICAL :: DO_ISOTROPIC

!  Flo flag

      LOGICAL :: DO_FLUORESCENCE

!   !@@ Solar sources + Observational Geometry flag !@@

      LOGICAL :: DO_SOLAR_SOURCES
      LOGICAL :: DO_USER_OBSGEOMS

!  Exact flag (!@@) and Exact only flag --> no Fourier term calculations

      LOGICAL :: DO_EXACT
      LOGICAL :: DO_EXACTONLY

!  Linearizations

      LOGICAL :: DO_SL_JACOBIANS
      LOGICAL :: DO_ISO_JACOBIANS
      INTEGER :: N_SLEAVE_WFS

!  Geometry and control
!  --------------------

!  Stream angle flag

      LOGICAL :: DO_USER_STREAMS

!  Number of Stokes components

      INTEGER :: NSTOKES

!  Number of discrete ordinate streams

      INTEGER :: NSTREAMS

!  Local angle control

      INTEGER :: NBEAMS
      INTEGER :: N_USER_STREAMS
      INTEGER :: N_USER_RELAZMS

!  Angles

      REAL(fpk) :: BEAM_SZAS   (MAXBEAMS)
      REAL(fpk) :: USER_RELAZMS(MAX_USER_RELAZMS)
      REAL(fpk) :: USER_ANGLES (MAX_USER_STREAMS)

!  !@@ Local Observational Geometry control and angles

      INTEGER    :: N_USER_OBSGEOMS
      REAL(fpk)  :: USER_OBSGEOMS (MAX_USER_OBSGEOMS,3)

!  Water-leaving variables
!  -----------------------

!  Input Salinity in [ppt]

      REAL(fpk) :: SALINITY

!  Input Chlorophyll concentration in [mg/M]

      REAL(fpk) :: CHLORCONC

!  Input wavelenth in [Microns]

      REAL(fpk) :: WAVELENGTH

!  Changed for Version 2.7
!     Input Wind speed in m/s, and azimuth directions relative to Sun positions

      REAL(fpk) :: WINDSPEED, WINDDIR ( MAXBEAMS )

!  Removed, Version 2.7 --> Quadrature is internal. 
!     Number of azimuth quadrature streams for reflectivity 
!        (only for non-isotropic water leaving)
!      INTEGER :: NSTREAMS_AZQUAD

!  New for Version 2.7.
!    Flags for glint shadowing, Foam Correction, facet Isotropy

      LOGICAL   :: DO_GlintShadow
      LOGICAL   :: DO_FoamOption
      LOGICAL   :: DO_FacetIsotropy

!  Linearization flags

      LOGICAL :: DO_CHLORCONC_WF
      LOGICAL :: DO_WINDSPEED_WF

!  Fluorescence variables
!  ----------------------

!  Input wavelength in [nm]

      REAL(fpk) :: FL_Wavelength

!  Input Latitude/Longitude in [degs]

      REAL(fpk) :: FL_Latitude, FL_Longitude 

!  Input Epoch

      INTEGER   :: FL_Epoch(6)

!  Input F755 Amplitude

      REAL(fpk)  :: FL_Amplitude755

!  Flag for using Data Gaussian parameters

      LOGICAL          :: FL_DO_DataGaussian

!  Linearization of Gaussians only when the data option is not set

      LOGICAL          :: FL_F755_JACOBIANS
      LOGICAL          :: FL_GAUSS_JACOBIANS(6)

!  Exception handling
!  ------------------

!     Message Length should be at least 120 Characters

      INTEGER ::             STATUS
      INTEGER ::             NMESSAGES
      CHARACTER (LEN=120) :: MESSAGES ( 0:MAX_MESSAGES )
      CHARACTER (LEN=120) :: ACTIONS ( 0:MAX_MESSAGES )

!  local variables
!  ===============

      CHARACTER (LEN=12), PARAMETER :: PREFIX = 'VSLEAVESUP -'

      LOGICAL ::            ERROR
      CHARACTER (LEN=80) :: PAR_STR
      INTEGER ::            I, FILUNIT, NM, M, QM

!  Initialize Exception handling

      STATUS = VLIDORT_SUCCESS

      MESSAGES(1:MAX_MESSAGES) = ' '
      ACTIONS (1:MAX_MESSAGES) = ' '

      NMESSAGES       = 0
      MESSAGES(0)     = 'Successful Read of VLIDORT Input file'
      ACTIONS(0)      = 'No Action required for this Task'

!  Local error handling initialization

      ERROR  = .FALSE.
      NM     = NMESSAGES

!  Open file

      FILUNIT = VLIDORT_INUNIT
      OPEN(VLIDORT_INUNIT,FILE=FILNAM,ERR=300,STATUS='OLD')

!  Initialize Angle control
!  ========================

      DO_SOLAR_SOURCES = .FALSE.  !@@ New line
      DO_USER_OBSGEOMS = .FALSE.  !@@ New line

      DO_USER_STREAMS = .FALSE.
      NSTREAMS = 0
      NSTOKES = 0

      NBEAMS   = 0
      DO I = 1, MAXBEAMS
        BEAM_SZAS(I) = ZERO
      ENDDO
      N_USER_STREAMS = 0
      DO I = 1, MAX_USER_STREAMS
        USER_ANGLES(I) = ZERO
      ENDDO
      N_USER_RELAZMS = 0
      DO I = 1, MAX_USER_RELAZMS
        USER_RELAZMS(I) = ZERO
      ENDDO

! !@@ Observational Geometry

      N_USER_OBSGEOMS = 0
      DO I = 1, MAX_USER_OBSGEOMS
        USER_OBSGEOMS(I,1:3) = ZERO
      ENDDO

!  Initialize Surface stuff
!  ========================

!  Control flags

      DO_EXACT        = .FALSE.       !@@  New line
      DO_EXACTONLY    = .FALSE.
      DO_ISOTROPIC    = .FALSE.
      DO_SLEAVING     = .FALSE.
      DO_FLUORESCENCE = .FALSE.

      DO_SL_JACOBIANS   = .false.
      DO_ISO_JACOBIANS  = .false.
      N_SLEAVE_WFS = 0

!  Fluorescence variables

      FL_LATITUDE   = ZERO
      FL_LONGITUDE  = ZERO
      FL_EPOCH      = 0
      FL_WAVELENGTH = ZERO
      FL_Amplitude755     = ZERO
      FL_DO_DataGaussian  = .false.

      FL_F755_JACOBIANS   = .false.
      FL_GAUSS_JACOBIANS  = .false.

!  Water-leaving variables

      SALINITY   = ZERO
      CHLORCONC  = ZERO
      WAVELENGTH = ZERO

      WINDSPEED  = ZERO
      WINDDIR    = ZERO

      DO_GlintShadow   = .FALSE.
      DO_FoamOption    = .FALSE.
      DO_FacetIsotropy = .FALSE.

      DO_CHLORCONC_WF  = .false.
      DO_WINDSPEED_WF   = .false.

!  Geometry and Input Control
!  ==========================

!  !@@ Solar sources is True, always

      DO_SOLAR_SOURCES = .TRUE.

!  User-defined Stream angle

      PAR_STR = 'Use user-defined viewing zenith angles?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
          READ (FILUNIT,*,ERR=998) DO_USER_STREAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  Number of Stokes components

      PAR_STR = 'Number of Stokes vector components'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
          READ (FILUNIT,*,ERR=998) NSTOKES
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  Discrete ordinates

      PAR_STR = 'Number of half-space streams'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
            READ (FILUNIT,*,ERR=998) NSTREAMS
      CALL FINDPAR_ERROR (ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                          ACTIONS )

!  All numbers are now checked against maximum dimensions

      IF ( NSTREAMS .GT. MAXSTREAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = &
        'Number of half-space streams > maximum dimension'
        ACTIONS(NM)  = &
         'Re-set input value or increase MAXSTREAMS dimension'
        STATUS = VLIDORT_SERIOUS
        NMESSAGES = NM
        GO TO 764
      ENDIF

!  Observational Geometry    !@@
!  ======================

!  !@@ New flag, Observational Geometry

      IF ( DO_SOLAR_SOURCES ) THEN
         PAR_STR = 'Do Observation Geometry?'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
              READ (FILUNIT,*,ERR=998) DO_USER_OBSGEOMS
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ENDIF

!  !@@ Observational Geometry control
!     ---- check not exceeding dimensioned number

      IF ( DO_USER_OBSGEOMS ) THEN
        PAR_STR = 'Number of Observation Geometry inputs'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) N_USER_OBSGEOMS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
        IF ( N_USER_OBSGEOMS .GT. MAX_USER_OBSGEOMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of ObsGeometry inputs > Maximum dimension'
          ACTIONS(NM)  = 'Re-set input or increase MAX_USER_OBSGEOMS dimension'
          STATUS       = VLIDORT_SERIOUS
          NMESSAGES    = NM
          GO TO 764
        ENDIF
      ENDIF

!  !@@ Observational Geometry control
!     Automatic setting of NBEAMS, N_USER_STREAMS, N_USER_RELAZMS

      IF ( DO_USER_OBSGEOMS ) THEN
         NBEAMS          = N_USER_OBSGEOMS
         N_USER_STREAMS  = N_USER_OBSGEOMS
         N_USER_RELAZMS  = N_USER_OBSGEOMS
         DO_USER_STREAMS = .TRUE.
      ENDIF

!  !@@ Observational Geometry control

      IF ( DO_USER_OBSGEOMS ) THEN
        PAR_STR = 'Observation Geometry inputs'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
           DO I = 1, N_USER_OBSGEOMS
             READ (FILUNIT,*,ERR=998)&
               USER_OBSGEOMS(I,1), USER_OBSGEOMS(I,2), USER_OBSGEOMS(I,3)
           ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ENDIF

!  Set angles

      IF ( DO_USER_OBSGEOMS ) THEN
         BEAM_SZAS   (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,1)
         USER_ANGLES (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,2)
         USER_RELAZMS(1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,3)
         GO TO 5667
      ENDIF

!  Solar beams
!  ===========

!  Number of Solar zenith angles

      PAR_STR = 'Number of solar zenith angles'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
            READ (FILUNIT,*,ERR=998) NBEAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  Check not exceeding dimensioned number

      IF ( NBEAMS .GT. MAXBEAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = &
        'Number of solar zenith angles > maximum dimension'
        ACTIONS(NM)  = &
        'Re-set input value or increase MAXBEAMS dimension'
        STATUS = VLIDORT_SERIOUS
        NMESSAGES = NM
        GO TO 764
      ENDIF

!  TOA solar zenith angle inputs

      PAR_STR = 'Solar zenith angles (degrees)'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
        DO I = 1, NBEAMS
          READ (FILUNIT,*,ERR=998) BEAM_SZAS(I)
        ENDDO
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  Azimuth angles
!  ==============

!  Number of Azimuth angles

      PAR_STR = 'Number of user-defined relative azimuth angles'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
            READ (FILUNIT,*,ERR=998) N_USER_RELAZMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  Check not exceeding dimensioned number

      IF ( N_USER_RELAZMS .GT. MAX_USER_RELAZMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = &
         'Number of relative azimuth angles > maximum dimension'
        ACTIONS(NM)  = &
         'Re-set input value or increase MAX_USER_RELAZMS dimension'
        STATUS       = VLIDORT_SERIOUS
        NMESSAGES    = NM
        GO TO 764
      ENDIF

! Azimuth angles

      PAR_STR = 'User-defined relative azimuth angles (degrees)'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
        DO I = 1, N_USER_RELAZMS
          READ (FILUNIT,*,ERR=998) USER_RELAZMS(I)
        ENDDO
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  User-defined viewing zenith angles  (should be positive)
!  ==================================

      IF ( DO_USER_STREAMS ) THEN

!  Number of User-defined Viewing zenith angles

        PAR_STR = 'Number of user-defined viewing zenith angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
            READ (FILUNIT,*,ERR=998) N_USER_STREAMS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                             ACTIONS )

!  Check dimension

        IF ( N_USER_STREAMS .GT. MAX_USER_STREAMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = &
          'Number of viewing zenith angles > maximum dimension'
          ACTIONS(NM)  = &
          'Re-set input value or increase MAX_USER_STREAMS dimension'
          STATUS = VLIDORT_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF

!  User-defined Viewing zenith angles

        PAR_STR = 'User-defined viewing zenith angles (degrees)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_USER_STREAMS
            READ (FILUNIT,*,ERR=998) USER_ANGLES(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                             ACTIONS )

      ENDIF

!  !@@ Continuation point for Skipping the Lattice-input angles

5667  continue

!  Surface stuff
!  =============

!  General SLEAVING input
!  ----------------------

!  Basic flag

      PAR_STR = 'Do surface-leaving Contributions?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_SLEAVING
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  Skip if not doing surface leaving
!   @@ Rob Fix 04 August 2014

      if ( .not.DO_SLEAVING ) GOTO 652

!  Isotropic flag

      PAR_STR = 'Do Isotropic surface-leaving?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
         READ (FILUNIT,*,ERR=998) DO_ISOTROPIC
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                         ACTIONS )

!  !@@ Overall-Exact flag

      PAR_STR = 'Do Overall-Exact surface-leaving?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
         READ (FILUNIT,*,ERR=998)DO_EXACT
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Exact only flag. Only if above is set (!@@)

      IF ( DO_EXACT ) THEN
        PAR_STR = 'Do Exact-only (no Fourier-term contributions)?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
           READ (FILUNIT,*,ERR=998)DO_EXACTONLY
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ENDIF

!  Basic source

      PAR_STR = 'Do surface-leaving Fluorescence?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_FLUORESCENCE
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  Basic linearization flags. Should be true in Isotropic cases

      PAR_STR = 'Do surface-leaving Jacobians?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_SL_JACOBIANS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

      IF ( DO_SL_JACOBIANS ) THEN
        PAR_STR = 'Do Isotropic surface-leaving Jacobians?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) DO_ISO_JACOBIANS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                             ACTIONS )
      ENDIF

!  Continuation point

652   continue

!  Local count for Jacobians

      QM = 0

!  Inputs for Water-leaving (Non-Fluorescence case)
!  ------------------------------------------------

      IF ( DO_SLEAVING.and..not.DO_FLUORESCENCE ) THEN

!  salinity, chlorophyll concentration, wavelength

        PAR_STR = 'Ocean water salinity [ppt]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) SALINITY
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        PAR_STR = 'Chlorophyll concentration in [mg/M]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) CHLORCONC
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        PAR_STR = 'Wavelength in [Microns]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) WAVELENGTH
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  New for Version 2.7.   Wind-speed and directions, used for the first time

        PAR_STR = 'Windspeed in [m/s]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) WINDSPEED
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        PAR_STR = 'Wind directions (degrees) relative to sun positions'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, NBEAMS
            READ (FILUNIT,*,ERR=998) WINDDIR(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  New for Version 2.7. Control Flag for including Whitecaps

      PAR_STR = 'Do whitecap (foam) calculation?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) Do_FoamOption
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  New for Version 2.7. Glint, Control Flags for Facet Isotropy and Shadowing
!     These are only needed for the non-Isotropic case

      if ( .not. do_Isotropic ) then
         PAR_STR = 'Do glint calculation with facet isotropy?'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) Do_FacetIsotropy
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
         PAR_STR = 'Do glint calculation with shadowing?'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_GlintShadow
!           READ (FILUNIT,*,ERR=998) Do_FacetIsotropy ! Bug Fix 9/27/14
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      endif

!  Removed for Version 2.7. --> Quadrature is now internal
!    Non-isotropic input = number of azimuth streams, check this value
!        IF ( .not. DO_ISOTROPIC ) THEN
!          PAR_STR = 'Number of azimuth quadrature streams'
!          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
!              READ (FILUNIT,*,ERR=998) NSTREAMS_AZQUAD
!          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
!          IF ( NSTREAMS_AZQUAD .GT. MAXSTREAMS_BRDF ) THEN
!            NM = NM + 1
!            MESSAGES(NM) = 'Number of AZQUAD streams > maximum dimension'
!            ACTIONS(NM)  = 'Re-set input value or increase MAXSTREAMS_BRDF dimension'
!            STATUS = VLIDORT_SERIOUS
!            NMESSAGES = NM
!            GO TO 764
!          ENDIF
!        ENDIF

!  Check Added 9/27/14. Multiple SZAS not allowed with Facet ANISOTROPY

        IF ( .not. Do_FacetIsotropy .and. NBEAMS.gt.1 ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Facet Anisotropy (wind direction) not Allowed for NBEAMS > 1'
          ACTIONS(NM)  = 'Either re-set NBEAMS = 1 or re-set Do_FacetIsotropy = .true.'
          STATUS = VLIDORT_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF

!  New for Version 2.7. Linearization flags

      IF ( DO_SL_JACOBIANS ) THEN
         PAR_STR = 'Do pigment concentration weighting function?'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) DO_CHLORCONC_WF
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
          IF ( DO_CHLORCONC_WF ) QM = QM + 1
         PAR_STR = 'Do wind speed weighting function?'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) DO_WINDSPEED_WF
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
          IF ( DO_WINDSPEED_WF ) QM = QM + 1
      ENDIF

!  Total number of Sleave weighting functions
!   Chack not out of bounds

        N_SLEAVE_WFS = QM
        if ( N_SLEAVE_WFS .gt. MAX_SLEAVEWFS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of Sleave Water-Leaving Jacobians > maximum dimension'
          ACTIONS(NM)  = 'Re-set input values or increase MAX_SLEAVEWFS dimension'
          STATUS = VLIDORT_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF

!  Inputs for Fluorescence Case
!  ----------------------------

      ELSE IF ( DO_SLEAVING.and.DO_FLUORESCENCE ) THEN

!  Temporary Check 

        IF ( .not. DO_ISOTROPIC ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'DO_ISOTROPIC was set to .FALSE. in fluorescence case'
          ACTIONS(NM)  = 'Tempo! Set DO_ISOTROPIC to .TRUE. if doing fluorescence'
          STATUS = VLIDORT_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF

!  Use of Data Gaussians (New, 8 August 2012)
!    IF NOT SET, YOU MUST USE YOUR OWN PARAMETERS

        PAR_STR = 'Do Data Gaussians in Fluorescence?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) FL_DO_DataGaussian
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  Amplitude for FS755 (Nominally, this is one)

        PAR_STR = 'Amplitude for Fluorescence model at 755 nm'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) FL_Amplitude755
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  Lat/Long, day-of-year, wavelength

        PAR_STR = 'Latitude for Fluorescence model [degs]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) FL_LATITUDE
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

        PAR_STR = 'Longitude for Fluorescence model [degs]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) FL_LONGITUDE
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

        PAR_STR = 'Epoch for Fluorescence model'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) FL_EPOCH(1:6)
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

        PAR_STR = 'Wavelength for Fluorescence model in [nm]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) FL_WAVELENGTH
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  Linearization inputs for Fluorescence

        IF ( DO_SL_JACOBIANS .and. DO_ISO_JACOBIANS ) THEN
           PAR_STR = 'Do Jacobians for F755 Fluorescence value?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) FL_F755_JACOBIANS
           CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                              ACTIONS )
        ENDIF

!  Gaussian linearizations
!     6 parameters: 1,2,3 for first Gaussian, 4,5,6 for  second Gaussian

        IF ( DO_SL_JACOBIANS .and. DO_ISO_JACOBIANS .and..not. FL_DO_DataGaussian ) THEN
           PAR_STR = 'Do Gaussian parameter Jacobians for Fluorescence?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
              DO M = 1, 6
                READ (FILUNIT,*,ERR=998) FL_GAUSS_JACOBIANS(m)
                IF ( FL_GAUSS_JACOBIANS(m) ) QM = QM + 1
              ENDDO
           ENDIF
           CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
        ENDIF

!  Total number of Sleave weighting functions
!   Chack not out of bounds

        N_SLEAVE_WFS = QM
        if ( N_SLEAVE_WFS .gt. MAX_SLEAVEWFS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of Sleave Fluorescence Jacobians > maximum dimension'
          ACTIONS(NM)  = 'Re-set input values or increase MAX_SLEAVEWFS dimension'
          STATUS = VLIDORT_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF

!  End fluorescence parameters

      ENDIF

!  Successful finish

      CLOSE(FILUNIT)

!mick fix
      NMESSAGES = NM

!  Copy General Control inputs

      VSLEAVE_Sup_In%SL_DO_USER_STREAMS  = DO_USER_STREAMS
      VSLEAVE_Sup_In%SL_DO_SLEAVING      = DO_SLEAVING
      VSLEAVE_Sup_In%SL_DO_FLUORESCENCE  = DO_FLUORESCENCE
      VSLEAVE_Sup_In%SL_DO_ISOTROPIC     = DO_ISOTROPIC
      VSLEAVE_Sup_In%SL_DO_EXACT         = DO_EXACT         !@@
      VSLEAVE_Sup_In%SL_DO_EXACTONLY     = DO_EXACTONLY
      VSLEAVE_Sup_In%SL_DO_SOLAR_SOURCES = DO_SOLAR_SOURCES   !@@
      VSLEAVE_Sup_In%SL_DO_USER_OBSGEOMS = DO_USER_OBSGEOMS   !@@

      VSLEAVE_LinSup_In%SL_DO_SL_JACOBIANS   = DO_SL_JACOBIANS
      VSLEAVE_LinSup_In%SL_DO_ISO_JACOBIANS  = DO_ISO_JACOBIANS

!  Copy Geometry results

      VSLEAVE_Sup_In%SL_NSTOKES           = NSTOKES
      VSLEAVE_Sup_In%SL_NSTREAMS          = NSTREAMS
      VSLEAVE_Sup_In%SL_NBEAMS            = NBEAMS
      VSLEAVE_Sup_In%SL_BEAM_SZAS         = BEAM_SZAS
      VSLEAVE_Sup_In%SL_N_USER_RELAZMS    = N_USER_RELAZMS
      VSLEAVE_Sup_In%SL_USER_RELAZMS      = USER_RELAZMS
      VSLEAVE_Sup_In%SL_N_USER_STREAMS    = N_USER_STREAMS
      VSLEAVE_Sup_In%SL_USER_ANGLES_INPUT = USER_ANGLES
      VSLEAVE_Sup_In%SL_N_USER_OBSGEOMS   = N_USER_OBSGEOMS !@@
      VSLEAVE_Sup_In%SL_USER_OBSGEOMS     = USER_OBSGEOMS   !@@

!  Copy Water-leaving inputs
!  -------------------------

!  Original

      VSLEAVE_Sup_In%SL_SALINITY         = SALINITY
      VSLEAVE_Sup_In%SL_CHLORCONC        = CHLORCONC
      VSLEAVE_Sup_In%SL_WAVELENGTH       = WAVELENGTH

!  Version 2.7 changes

!      VSLEAVE_Sup_In%SL_NSTREAMS_AZQUAD  = NSTREAMS_AZQUAD
      VSLEAVE_Sup_In%SL_WINDSPEED        = WINDSPEED
      VSLEAVE_Sup_In%SL_WINDDIR          = WINDDIR

      VSLEAVE_Sup_In%SL_DO_GlintShadow   = DO_GlintShadow
      VSLEAVE_Sup_In%SL_DO_FoamOption    = DO_FoamOption
      VSLEAVE_Sup_In%SL_DO_FacetIsotropy = DO_FacetIsotropy

      VSLEAVE_LinSup_In%SL_DO_CHLORCONC_WF = DO_CHLORCONC_WF
      VSLEAVE_LinSup_In%SL_DO_WINDSPEED_WF = DO_WINDSPEED_WF

!  Copy Fluorescence inputs
!  ------------------------

      VSLEAVE_Sup_In%SL_FL_LATITUDE        = FL_LATITUDE
      VSLEAVE_Sup_In%SL_FL_LONGITUDE       = FL_LONGITUDE
      VSLEAVE_Sup_In%SL_FL_WAVELENGTH      = FL_WAVELENGTH
      VSLEAVE_Sup_In%SL_FL_EPOCH           = FL_EPOCH
      VSLEAVE_Sup_In%SL_FL_Amplitude755    = FL_Amplitude755
      VSLEAVE_Sup_In%SL_FL_DO_DataGaussian = FL_DO_DataGaussian

      VSLEAVE_LinSup_In%SL_FL_F755_JACOBIANS  = FL_F755_JACOBIANS
      VSLEAVE_LinSup_In%SL_FL_GAUSS_JACOBIANS = FL_GAUSS_JACOBIANS

!  total number of weighting functions

      VSLEAVE_LinSup_In%SL_N_SLEAVE_WFS = N_SLEAVE_WFS

!  Exception handling

      VSLEAVE_Sup_InputStatus%SL_STATUS_INPUTREAD = STATUS
      VSLEAVE_Sup_InputStatus%SL_NINPUTMESSAGES   = NMESSAGES
      VSLEAVE_Sup_InputStatus%SL_INPUTMESSAGES    = MESSAGES
      VSLEAVE_Sup_InputStatus%SL_INPUTACTIONS     = ACTIONS

!  Normal return

      RETURN

!  Open file error

300   CONTINUE
      STATUS = VLIDORT_SERIOUS
      NMESSAGES = NMESSAGES + 1
      MESSAGES(NMESSAGES) = 'openfile failure for '//adjustl(trim(FILNAM))
      ACTIONS(NMESSAGES)  = 'Find the Right input file!!'
      CLOSE(FILUNIT)
      GO TO 764

!  line read error - abort immediately

998   CONTINUE
      STATUS = VLIDORT_SERIOUS
      NMESSAGES = NMESSAGES + 1
      MESSAGES(NMESSAGES) = 'read failure for '//adjustl(trim(FILNAM))
      ACTIONS(NMESSAGES)  = 'Re-set: Entry is incorrect in input file'
      CLOSE(FILUNIT)
      GO TO 764

!  Final error copying

764   CONTINUE

      VSLEAVE_Sup_InputStatus%SL_STATUS_INPUTREAD = STATUS
      VSLEAVE_Sup_InputStatus%SL_NINPUTMESSAGES   = NMESSAGES
      VSLEAVE_Sup_InputStatus%SL_INPUTMESSAGES    = MESSAGES
      VSLEAVE_Sup_InputStatus%SL_INPUTACTIONS     = ACTIONS

!  Finish

      RETURN
      END SUBROUTINE VSLEAVE_LIN_INPUTMASTER

!

      SUBROUTINE VSLEAVE_LIN_MAINMASTER ( &
        VSLEAVE_Sup_In,            & ! Inputs
        VSLEAVE_LinSup_In,         & ! Inputs
        VSLEAVE_Sup_Out,           & ! Outputs
        VSLEAVE_LinSup_Out )         ! Outputs

!  Prepares the Surface Leaving (and linearizations) necessary for VLIDORT.

!  Version 2.6 Notes
!  -----------------

!  INDWAT, MORCASIWAT Routines taken straight from Clark Weaver code
!      Compiled here by R. Spurr, 11 July 2012
!  get_fluorescence_755 Routines taken straight from Chris O'dell code
!      Compiled here by R. Spurr, 12 July 2012

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Upgrade for Version 2.7
!  -----------------------

!  New Water-leaving code developed/tested by R. Spurr, April 2014
!  Based in part on Modified-6S code by A. Sayer.
!  Validated against Modified-6 OCEABRDF.F code, 24  April 2014.

!  Water-leaving upgrade according to the modified 6S specification
!  New code calculates transmittances into and out of ocean, using
!  usual sun-glint rough surface approximations. In addition, the
!  water-leaving term itself is now SZA-dependent (A. Sayer), and
!  there is now a correction for Whitecaps (again, from 6S)

!  This upgrade gives the  water-leaving terms some directionality,
!  but they are still azimuth-independent

!  Earlier version inputs were just Wavelength/Salinity/PigmentConc
!  This was enough for the isotropic case (Fast Option) in Version 2.6.
!  For Version 2.7, we require additional inputs, including:
!    - Wind-speed and Direction (direction was not used in earlier version) 
!    - flags to control use sunglint shadowing and foam (whitecaps) correction.

!  This water-leaving option is designed to work alongside the "NewCM" glint
!  reflectance option in the VBRDF code. The glint and whitecap calculations
!  in the two supplements are the same.

!  You need to make sure that the wind input information is the same as that
!  used for the "NewCM" glint option in the VBRDF supplement. Also, the Foam
!  correction applied here in the surface-leaving code should also be 
!  applied for "NewCM" glint.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      USE VLIDORT_PARS

      USE vsleave_sup_inputs_def
      USE vsleave_linsup_inputs_def

      USE vsleave_sup_outputs_def
      USE vsleave_linsup_outputs_def

      USE vsleave_sup_aux_m          , only : VSleave_GAULEG
      USE vsleave_sup_routines_m     , only : WaterLeaving,              &
                                              get_fluorescence_755,      &
                                              solar_spec_irradiance
      USE vsleave_lin_sup_routines_m , only : Linearized_WaterLeaving

      IMPLICIT NONE

!  Input structures
!  ----------------

      TYPE(VSLEAVE_Sup_Inputs)   , INTENT(IN)   :: VSLEAVE_Sup_In
      TYPE(VSLEAVE_LinSup_Inputs), INTENT(IN)   :: VSLEAVE_LinSup_In

!  Output structures
!  -----------------

      TYPE(VSLEAVE_Sup_Outputs)   , INTENT(OUT) :: VSLEAVE_Sup_Out
      TYPE(VSLEAVE_LinSup_Outputs), INTENT(OUT) :: VSLEAVE_LinSup_Out

!  VLIDORT local variables
!  +++++++++++++++++++++++

!  Input arguments
!  ===============

!  Main Boolean flags
!  ------------------

!  Inclusion flag (not really necessary, Brian)

      LOGICAL :: DO_SLEAVING

!  Isotropic flag

      LOGICAL :: DO_ISOTROPIC

!  Flo flag

      LOGICAL :: DO_FLUORESCENCE

!  !@@ Solar sources + Observational Geometry flag !@@

      LOGICAL :: DO_SOLAR_SOURCES
      LOGICAL :: DO_USER_OBSGEOMS

!  Exact flag (!@@) and Exact only flag --> no Fourier term calculations

      LOGICAL :: DO_EXACT
      LOGICAL :: DO_EXACTONLY

!  Linearizations

      LOGICAL :: DO_SL_JACOBIANS
      LOGICAL :: DO_ISO_JACOBIANS

!  Geometry and control
!  --------------------

!  Stream angle flag

      LOGICAL ::   DO_USER_STREAMS

!  Number of Stokes components

      INTEGER ::   NSTOKES

!  Number of discrete ordinate streams, quadrature
!   Version 2.7, added the qudrature arrays

      INTEGER   :: NSTREAMS
      REAL(fpk) :: STREAMS (MAXSTREAMS)
      REAL(fpk) :: WEIGHTS (MAXSTREAMS)

!  Local angle control

      INTEGER   :: NBEAMS
      INTEGER   :: N_USER_STREAMS
      INTEGER   :: N_USER_RELAZMS

!  Angles

      REAL(fpk) :: BEAM_SZAS   (MAXBEAMS)
      REAL(fpk) :: USER_RELAZMS(MAX_USER_RELAZMS)
      REAL(fpk) :: USER_ANGLES (MAX_USER_STREAMS)

!  !@@ Local Observational Geometry control and angles

      INTEGER    :: N_USER_OBSGEOMS
      REAL(fpk)  :: USER_OBSGEOMS (MAX_USER_OBSGEOMS,3)

!  Water-leaving variables
!  -----------------------

!  Input Salinity in [ppt]

      REAL(fpk) :: SALINITY

!  Input Chlorophyll concentration in [mg/M]

      REAL(fpk) :: CHLORCONC

!  Input wavelenth in [Microns]

      REAL(fpk) :: WAVELENGTH

!  Changed for Version 2.7
!     Input Wind speed in m/s, and azimuth directions relative to Sun positions

      REAL(fpk) :: WINDSPEED, WINDDIR ( MAXBEAMS )

!  Removed, Version 2.7 --> Quadrature is internal. 
!     Number of azimuth quadrature streams for reflectivity 
!        (only for non-isotropic water leaving)
!      INTEGER :: NSTREAMS_AZQUAD

!  New for Version 2.7.
!    Flags for glint shadowing, Foam Correction, facet Isotropy

      LOGICAL   :: DO_GlintShadow
      LOGICAL   :: DO_FoamOption
      LOGICAL   :: DO_FacetIsotropy

!  Linearization flags

      LOGICAL   :: DO_CHLORCONC_WF
      LOGICAL   :: DO_WINDSPEED_WF

!  Fluorescence variables
!  ----------------------

!  Input wavelength in [nm]

      REAL(fpk) :: FL_Wavelength

!  Input Latitude/Longitude in [degs]

      REAL(fpk) :: FL_Latitude, FL_Longitude

!  Input Epoch

      INTEGER   :: FL_Epoch(6)

!  Input F755 Amplitude

      REAL(fpk)  :: FL_Amplitude755

!  Flag for using Data Gaussian parameters

      LOGICAL    :: FL_DO_DataGaussian

!  Linearization flags
!     Linearization of Gaussians only when the data option is not set

      LOGICAL          :: FL_F755_JACOBIANS, FL_GAUSS_JACOBIANS(6)

!  Local functions
!  ===============

!  Fluorescence Isotropic Surface leaving term

   REAL(fpk) :: Fluor_ISOTROPIC

!  Isotropic value. Fast calculation

   REAL(fpk)    :: WLeaving_ISO    ( MAXBEAMS )
   REAL(fpk)    :: dC_WLeaving_ISO ( MAXBEAMS )
   REAL(fpk)    :: dW_WLeaving_ISO ( MAXBEAMS )

!  Input solar, output stream angles

   REAL(fpk)    :: WLeaving_SD    ( MAXBEAMS, MAXSTREAMS )
   REAL(fpk)    :: dC_WLeaving_SD ( MAXBEAMS, MAXSTREAMS )
   REAL(fpk)    :: dW_WLeaving_SD ( MAXBEAMS, MAXSTREAMS )

!  input solar, output view angles

   REAL(fpk)    :: WLeaving_SV    ( MAXBEAMS, MAX_USER_STREAMS )
   REAL(fpk)    :: dC_WLeaving_SV ( MAXBEAMS, MAX_USER_STREAMS )
   REAL(fpk)    :: dW_WLeaving_SV ( MAXBEAMS, MAX_USER_STREAMS )

!  Exact Surface-Leaving term
!   REAL(fpk), dimension ( MAXSTOKES, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS ) :: SLTERM_USERANGLES
!  Fourier components of Surface-leaving terms:
!    Every solar direction, SL-transmitted quadrature streams
!    Every solar direction, SL-transmitted user streams
!   REAL(fpk), dimension ( 0:MAXMOMENTS, MAXSTOKES, MAXSTREAMS, MAXBEAMS )       :: SLTERM_F_0
!   REAL(fpk), dimension ( 0:MAXMOMENTS, MAXSTOKES, MAX_USER_STREAMS, MAXBEAMS ) :: USER_SLTERM_F_0

!  Other local variables
!  =====================

!  (Version 2.6 code). Water-leaving model
!      REAL :: WAV,CHL,RW,SAL,A,REFR,REFI,N12,RWB,TDS,TDV

!  Fluorescence Gaussian parameters (Data)
!     Parameters of the fluorescence Gaussian spectral shape model.
!           Gaussian    A (Wm−2 μm−1 sr−1) Lambda(nm) Sigma(nm)
!              1           1.445           736.8        21.2
!              2           0.868           685.2        9.55

      REAL(FPK) :: FL_DataGAUSSIANS(3,2), FL_GAUSSIANS(3,2)
      data FL_DataGAUSSIANS(1,1) / 1.445d0 /
      data FL_DataGAUSSIANS(2,1) / 736.8d0 /
      data FL_DataGAUSSIANS(3,1) / 21.2d0  /
      data FL_DataGAUSSIANS(1,2) / 0.868d0 /
      data FL_DataGAUSSIANS(2,2) / 685.2d0 /
      data FL_DataGAUSSIANS(3,2) / 9.55d0  /

!  Solar spectral radiance model

      REAL(FPK) :: ssr_wvl

!  Fluorescence model

      CHARACTER*60 :: Fluofile
      INTEGER   :: IB, K, nslpars, M, Q, UM
      REAL(FPK) :: Fs755(MAXBEAMS), FL_SunSpec, FsSum, d_Fssum(6)
      REAL(FPK) :: ampli, lamda, sigma, arg, gauss, var, ff, expl, help1, help2
      !REAL(FPK) :: solar_spec_irradiance

      INTEGER, PARAMETER :: LUM = 1   !@@
      INTEGER, PARAMETER :: LUA = 1   !@@

!  Copy from input structure
!  -------------------------

!  Copy Top-level general Control inputs

      DO_USER_STREAMS  = VSLEAVE_Sup_In%SL_DO_USER_STREAMS
      DO_SLEAVING      = VSLEAVE_Sup_In%SL_DO_SLEAVING

      DO_EXACT         = VSLEAVE_Sup_In%SL_DO_EXACT          !@@
      DO_EXACTONLY     = VSLEAVE_Sup_In%SL_DO_EXACTONLY
      DO_FLUORESCENCE  = VSLEAVE_Sup_In%SL_DO_FLUORESCENCE
      DO_ISOTROPIC     = VSLEAVE_Sup_In%SL_DO_ISOTROPIC

!  !@@ New lines

      DO_SOLAR_SOURCES = VSLEAVE_Sup_In%SL_DO_SOLAR_SOURCES
      DO_USER_OBSGEOMS = VSLEAVE_Sup_In%SL_DO_USER_OBSGEOMS

      DO_SL_JACOBIANS  = VSLEAVE_LinSup_In%SL_DO_SL_JACOBIANS
      DO_ISO_JACOBIANS = VSLEAVE_LinSup_In%SL_DO_ISO_JACOBIANS

!  Set number of stokes elements and streams
!   Stream Quadrature new for Version 2.7

      NSTOKES  = VSLEAVE_Sup_In%SL_NSTOKES
      NSTREAMS = VSLEAVE_Sup_In%SL_NSTREAMS
      CALL VSleave_Gauleg ( zero, one, STREAMS(1:nstreams), WEIGHTS(1:nstreams), NSTREAMS )

!   !@@ Observational Geometry + Solar sources Optionalities
!   !@@ Either set from User Observational Geometry
!          Or Copy from Usual lattice input

      IF ( DO_USER_OBSGEOMS ) THEN
        N_USER_OBSGEOMS = VSLEAVE_Sup_In%SL_N_USER_OBSGEOMS
        USER_OBSGEOMS   = VSLEAVE_Sup_In%SL_USER_OBSGEOMS
        IF ( DO_SOLAR_SOURCES ) THEN
          NBEAMS          = N_USER_OBSGEOMS
          N_USER_STREAMS  = N_USER_OBSGEOMS
          N_USER_RELAZMS  = N_USER_OBSGEOMS
          BEAM_SZAS   (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,1)
          USER_ANGLES (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,2)
          USER_RELAZMS(1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,3)
        ELSE
          NBEAMS         = 1 ; BEAM_SZAS      = ZERO
          N_USER_RELAZMS = 1 ; USER_RELAZMS   = ZERO
          N_USER_STREAMS = N_USER_OBSGEOMS
          USER_ANGLES(1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,2)
        ENDIF
      ELSE
        IF ( DO_SOLAR_SOURCES ) THEN
          NBEAMS         = VSLEAVE_Sup_In%SL_NBEAMS
          BEAM_SZAS      = VSLEAVE_Sup_In%SL_BEAM_SZAS
          N_USER_RELAZMS = VSLEAVE_Sup_In%SL_N_USER_RELAZMS
          USER_RELAZMS   = VSLEAVE_Sup_In%SL_USER_RELAZMS
          N_USER_STREAMS = VSLEAVE_Sup_In%SL_N_USER_STREAMS
          USER_ANGLES    = VSLEAVE_Sup_In%SL_USER_ANGLES_INPUT
        ELSE
          NBEAMS         = 1 ; BEAM_SZAS      = ZERO
          N_USER_RELAZMS = 1 ; USER_RELAZMS   = ZERO
          N_USER_STREAMS = VSLEAVE_Sup_In%SL_N_USER_STREAMS
          USER_ANGLES    = VSLEAVE_Sup_In%SL_USER_ANGLES_INPUT
        ENDIF
      ENDIF

!  Copy Water-leaving inputs
!  -------------------------

!  Original

      SALINITY        = VSLEAVE_Sup_In%SL_SALINITY
      CHLORCONC       = VSLEAVE_Sup_In%SL_CHLORCONC
      WAVELENGTH      = VSLEAVE_Sup_In%SL_WAVELENGTH

!  Version 2.7 changes

!      NSTREAMS_AZQUAD = VSLEAVE_Sup_In%SL_NSTREAMS_AZQUAD
      WINDSPEED       = VSLEAVE_Sup_In%SL_WINDSPEED
      WINDDIR         = VSLEAVE_Sup_In%SL_WINDDIR

      DO_GlintShadow   = VSLEAVE_Sup_In%SL_DO_GlintShadow
      DO_FoamOption    = VSLEAVE_Sup_In%SL_DO_FoamOption
      DO_FacetIsotropy = VSLEAVE_Sup_In%SL_DO_FacetIsotropy

      DO_ChlorConc_WF  = VSLEAVE_LinSup_In%SL_DO_ChlorConc_WF
      DO_WindSpeed_WF  = VSLEAVE_LinSup_In%SL_DO_WindSpeed_WF

!  Copy Fluorescence inputs
!  ------------------------

!  Main variables

      FL_Wavelength   = VSLEAVE_Sup_In%SL_FL_Wavelength
      FL_Latitude     = VSLEAVE_Sup_In%SL_FL_Latitude
      FL_Longitude    = VSLEAVE_Sup_In%SL_FL_Longitude
      FL_Epoch        = VSLEAVE_Sup_In%SL_FL_Epoch
      FL_Amplitude755 = VSLEAVE_Sup_In%SL_FL_Amplitude755
      FL_DO_DataGaussian = VSLEAVE_Sup_In%SL_FL_DO_DataGaussian

!mick fix 8/31/2012 - added outer if block
      if (DO_FLUORESCENCE) then
        if ( FL_DO_DataGaussian ) then
           FL_GAUSSIANS(1:3,1) = FL_DataGAUSSIANS(1:3,1) 
           FL_GAUSSIANS(1:3,2) = FL_DataGAUSSIANS(1:3,2)
        else
           FL_GAUSSIANS(1:3,1) = VSLEAVE_Sup_In%SL_FL_InputGAUSSIANS(1:3,1) 
           FL_GAUSSIANS(1:3,2) = VSLEAVE_Sup_In%SL_FL_InputGAUSSIANS(1:3,2)
        endif
      endif

      FL_F755_JACOBIANS   = VSLEAVE_LinSup_In%SL_FL_F755_JACOBIANS
      FL_GAUSS_JACOBIANS  = VSLEAVE_LinSup_In%SL_FL_GAUSS_JACOBIANS

!  Main code
!  ---------

!  Zero the output

      VSLEAVE_Sup_Out%SL_SLTERM_ISOTROPIC  = ZERO
      VSLEAVE_Sup_Out%SL_SLTERM_USERANGLES = ZERO
      VSLEAVE_Sup_Out%SL_SLTERM_F_0        = ZERO
      VSLEAVE_Sup_Out%SL_USER_SLTERM_F_0   = ZERO

!  Zero the linearized output (8/8/12, only isotropic defined so far)

      VSLEAVE_LinSup_Out%SL_N_SLEAVE_WFS         = 0
      VSLEAVE_LinSup_Out%SL_LS_SLTERM_ISOTROPIC  = ZERO
      VSLEAVE_LinSup_Out%SL_LS_SLTERM_USERANGLES = ZERO
      VSLEAVE_LinSup_Out%SL_LS_SLTERM_F_0        = ZERO
      VSLEAVE_LinSup_Out%SL_LS_USER_SLTERM_F_0   = ZERO

!  Return if not called
!   ROb Fix, 04 August 2014

      IF ( .not. DO_SLEAVING ) RETURN

!  Fluorescence
!  ============

      IF ( DO_FLUORESCENCE ) THEN

!  Temporary - Only Isotropic yet.

        IF ( .not.DO_ISOTROPIC ) &
          Stop'Non-isotropic not allowed yet if doing fluorescence'

!  F_755 data file

        Fluofile = 'vlidort_v_test/data/fluor_data_2009_fortran.dat'

!  Get solar spectral irradiance, in (W m−2 μm−1), to normalize data

        !FL_SunSpec = 1.0d0  ! Temporary

        ssr_wvl = FL_Wavelength*1.0d-3 !convert from nm to um
        FL_SunSpec = solar_spec_irradiance( ssr_wvl )

!  factor: After  some fiddling, this is 1.0 (July 30th, 2012)
!    4pi factor is required in DB correction ter,

!         FF = PI4
        FF = 1.0d0
!        FF = 0.0d0

!  Zero the running total of weighting functions

        NSLPARS = 0

!  For each Solar zenith angle

        DO IB = 1, NBEAMS

 !  Get the F_755 data from the subroutine

          CALL get_fluorescence_755 &
   ( FL_Latitude, FL_Longitude, FL_Epoch, BEAM_SZAS(IB), FluoFile, Fs755(IB) )

!  Compute double Gaussian sums

          FsSum = zero
          do k = 1, 2
             ampli = FL_Gaussians(1,k)
             lamda = FL_Gaussians(2,k)
             sigma = FL_Gaussians(3,k)
             var = 0.5d0/sigma/sigma
             arg = ( FL_Wavelength - lamda ) * ( FL_Wavelength - lamda ) * var
             Gauss = zero
             if ( arg.lt.88.0d0 ) gauss = ampli * dexp ( - arg )
             FsSum = FsSum + Gauss
          enddo

!  Assign output Fluorescence (Apply Amplitude)
!  multiply by Fs755, and normalize to solar spectrum
!   FF is the fudge factor

          help2 = FF * FsSum / FL_SunSpec
          Fluor_ISOTROPIC = help2 * Fs755(IB)
          VSLEAVE_Sup_Out%SL_SLTERM_ISOTROPIC(1,IB) = FL_Amplitude755 * Fluor_ISOTROPIC

!  Apply Main linearization (w.r.t Fs755 Amplitude)

          if ( DO_ISO_JACOBIANS .and. FL_F755_JACOBIANS ) then
             NSLPARS = 1 ; VSLEAVE_LinSup_Out%SL_LS_SLTERM_ISOTROPIC(1,1,IB) = Fluor_ISOTROPIC
          endif

!  Applyl Gaussian linearizations if flagged.

          if ( DO_ISO_JACOBIANS .and. .not. FL_DO_DataGaussian ) then
             help1 = FL_Amplitude755 * FF * Fs755(IB) / FL_SunSpec;  d_Fssum = zero
             do k = 1, 2
                ampli = FL_Gaussians(1,k)
                lamda = FL_Gaussians(2,k)
                sigma = FL_Gaussians(3,k)
                var = 0.5d0/sigma/sigma
                arg = ( FL_Wavelength - lamda ) * ( FL_Wavelength - lamda ) * var
                m = 3*(k-1);  Gauss = zero
                if ( arg.lt.88.0d0 ) then
                   expl  = dexp ( - arg )
                   gauss = ampli * expl
                   d_Fssum(m+1)  = expl
                   d_Fssum(m+2)  = + gauss * two * ( FL_Wavelength - lamda ) * var
                   d_Fssum(m+3)  = + gauss * two * arg / sigma
                endif
             enddo
             do m = 1, 6
                IF ( FL_GAUSS_JACOBIANS(m) ) then
                   NSLPARS = NSLPARS + 1
                   VSLEAVE_LinSup_Out%SL_LS_SLTERM_ISOTROPIC(NSLPARS,1,IB) = help1 * d_Fssum(m)
                ENDIF
             enddo
          endif

!  Total number of weighting functions

          VSLEAVE_LinSup_Out%SL_N_SLEAVE_WFS = NSLPARS

!        write(*,*)FL_Wavelength,FsSum, FL_SunSpec, SLTERM_ISOTROPIC

!  End Beam loop

        ENDDO

      ENDIF

!  WATER-LEAVING
!  =============

      IF ( .not. DO_FLUORESCENCE ) THEN

!  Call the routine for Version 2.7

         if ( .not. DO_CHLORCONC_WF .and. .not. DO_WINDSPEED_WF ) then
            Call WaterLeaving &
             ( Maxbeams, Max_User_streams, Maxstreams,                            &
               do_Isotropic, Do_FoamOption, Do_GlintShadow, Do_FacetIsotropy,     &
               Wavelength, Salinity, ChlorConc, Windspeed, WindDir,               &
               nbeams, n_user_streams, nstreams, beam_szas, user_angles, streams, &
               WLeaving_ISO, WLeaving_SD, WLeaving_SV )
         else
            Call Linearized_WaterLeaving &
             ( Maxbeams, Max_User_streams, Maxstreams,                            &
               do_Isotropic, Do_FoamOption, Do_GlintShadow, Do_FacetIsotropy,     &
               Wavelength, Salinity, ChlorConc, Windspeed, WindDir,               &
               nbeams, n_user_streams, nstreams, beam_szas, user_angles, streams, &
               WLeaving_ISO, WLeaving_SD, WLeaving_SV,                            &
               dC_WLeaving_ISO, dC_WLeaving_SD, dC_WLeaving_SV,                   &
               dW_WLeaving_ISO, dW_WLeaving_SD, dW_WLeaving_SV )
         endif

!  Copy to Type structure arrays
!  -----------------------------

!    If Isotropic, VLIDORT takes care of the rest

         VSLEAVE_Sup_Out%SL_SLTERM_ISOTROPIC(1,1:NBEAMS) = WLeaving_ISO(1:NBEAMS)
         if ( .not. do_Isotropic ) then
            do ib = 1, nbeams
               do um = 1, n_user_streams
                  VSLEAVE_Sup_Out%SL_SLTERM_USERANGLES(1,um,:1:n_user_relazms,ib) = WLeaving_SV(ib,um)
                  VSLEAVE_Sup_Out%SL_USER_SLTERM_F_0  (0,1,um,ib)                 = WLeaving_SV(ib,um)
               enddo
               VSLEAVE_Sup_Out%SL_SLTERM_F_0(0,1,1:nstreams,ib) = WLeaving_SD(ib,1:nstreams)
            enddo
         endif

!  Copy to linearized Type structure arrays
!  ----------------------------------------

         Q = 0 ! start Jacobian count

!   Pigment

         if ( DO_CHLORCONC_WF ) then
            Q = Q + 1 ; VSLEAVE_LinSup_Out%SL_LS_SLTERM_ISOTROPIC(Q,1,1:NBEAMS) = dC_WLeaving_ISO(1:NBEAMS)
            if ( .not. do_Isotropic ) then
               do ib = 1, nbeams
                  do um = 1, n_user_streams
                     VSLEAVE_LinSup_Out%SL_LS_SLTERM_USERANGLES(Q,1,um,1:n_user_relazms,ib) = dC_WLeaving_SV(ib,um)
                     VSLEAVE_LinSup_Out%SL_LS_USER_SLTERM_F_0  (Q,0,1,um,ib)                = dC_WLeaving_SV(ib,um)
                  enddo
                  VSLEAVE_LinSup_Out%SL_LS_SLTERM_F_0(Q,0,1,1:nstreams,ib) = dC_WLeaving_SD(ib,1:nstreams)
               enddo
            endif
         endif

!   Windspeed

         if ( DO_WINDSPEED_WF ) then
            Q = Q + 1 ; VSLEAVE_LinSup_Out%SL_LS_SLTERM_ISOTROPIC(Q,1,1:NBEAMS) = dW_WLeaving_ISO(1:NBEAMS)
            if ( .not. do_Isotropic ) then
               do ib = 1, nbeams
                  do um = 1, n_user_streams
                     VSLEAVE_LinSup_Out%SL_LS_SLTERM_USERANGLES(Q,1,um,1:n_user_relazms,ib) = dW_WLeaving_SV(ib,um)
                     VSLEAVE_LinSup_Out%SL_LS_USER_SLTERM_F_0  (Q,0,1,um,ib)                = dW_WLeaving_SV(ib,um)
                  enddo
                  VSLEAVE_LinSup_Out%SL_LS_SLTERM_F_0(Q,0,1,1:nstreams,ib) = dW_WLeaving_SD(ib,1:nstreams)
               enddo
            endif
         endif

!  Total number of WFs

         VSLEAVE_LinSup_Out%SL_N_SLEAVE_WFS = Q

!  Here is the compressed Version 2.6 code ---------------------------------------
!    INDWAT/MORCASIWAT calls . use single precision routine
!        SAL = REAL(SALINITY) ; WAV = REAL(WAVELENGTH)
!        CALL INDWAT(WAV,SAL,refr,refi)
!        CHL = REAL(CHLORCONC) ; WAV = REAL(WAVELENGTH)
!        CALL MORCASIWAT(WAV,CHL,RW,.false.)
!     Perfect Transmittance. Add change in solid angle through surface
!     that accounts for 1/(n12*n12) decrease in directional reflectivity
!        a   = 0.485 ; tds = 1.0 ; tdv = 1.0
!        n12 = refr*refr + refi*refi  ; n12 = sqrt(n12)
!        Rwb=(1.0/(n12*n12))*tds*tdv*Rw/(1-a*Rw)
!        VSLEAVE_Sup_Out%SL_SLTERM_ISOTROPIC(1,1:NBEAMS)= DBLE(Rwb)

!  end water leaving

      endif

!  Finish

      RETURN
      END SUBROUTINE VSLEAVE_LIN_MAINMASTER

      END MODULE vsleave_linsup_masters_m

