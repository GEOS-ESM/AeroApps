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
! #                   2.5, 2.6                                  #
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
! #            VLIDORT_VBRDF_INPUT_CHECKER                      #
! #            VLIDORT_VSLEAVE_INPUT_CHECKER                    #
! #            VBRDF_VSLEAVE_INPUT_CHECKER                      #
! #                                                             #
! ###############################################################


      MODULE vlidort_sup_accessories

      PRIVATE
      PUBLIC :: VLIDORT_VBRDF_INPUT_CHECKER,   &
                VLIDORT_VSLEAVE_INPUT_CHECKER, &
                VBRDF_VSLEAVE_INPUT_CHECKER

      CONTAINS

      SUBROUTINE VLIDORT_VBRDF_INPUT_CHECKER ( &
        VBRDF_Sup_In,             & ! Inputs
        VLIDORT_FixIn,            & ! Inputs
        VLIDORT_ModIn,            & ! Inputs
        VLIDORT_VBRDFCheck_Status ) ! Outputs

      USE VLIDORT_PARS

      USE VBRDF_Sup_Inputs_def
      USE VLIDORT_Inputs_def
      USE VLIDORT_Outputs_def

      IMPLICIT NONE

      TYPE(VBRDF_Sup_Inputs), INTENT(IN)             :: VBRDF_Sup_In

      TYPE(VLIDORT_Fixed_Inputs), INTENT (IN)        :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs), INTENT (IN)     :: VLIDORT_ModIn

      TYPE(VLIDORT_Exception_Handling), INTENT(OUT)  :: VLIDORT_VBRDFCheck_Status

!  ---------------
!  Local variables
!  ---------------

!  BRDF supplement inputs
!  ----------------------

!  User stream, BRDF surface and surface emission flags

      LOGICAL ::             BS_DO_BRDF_SURFACE
      LOGICAL ::             BS_DO_USER_STREAMS
      LOGICAL ::             BS_DO_SURFACE_EMISSION

!  Number of Stokes vector components

      INTEGER ::             BS_NSTOKES

!  Number of discrete ordinate streams

      INTEGER ::             BS_NSTREAMS

!  BOA solar zenith angles

      INTEGER ::             BS_NBEAMS
      DOUBLE PRECISION ::    BS_BEAM_SZAS ( MAXBEAMS )

!  User-defined relative azimuths (mandatory for Fourier > 0)

      INTEGER ::             BS_N_USER_RELAZMS
      DOUBLE PRECISION ::    BS_USER_RELAZMS (MAX_USER_RELAZMS)

!  User-defined zenith angle input

      INTEGER ::             BS_N_USER_STREAMS
      DOUBLE PRECISION ::    BS_USER_ANGLES_INPUT (MAX_USER_STREAMS)

!  VLIDORT Main inputs
!  -------------------

      LOGICAL ::             DO_USER_VZANGLES
      LOGICAL ::             DO_MVOUT_ONLY
      LOGICAL ::             DO_LAMBERTIAN_SURFACE
      LOGICAL ::             DO_SURFACE_EMISSION

      INTEGER ::             NSTOKES
      INTEGER ::             NSTREAMS

      INTEGER ::             N_SZANGLES
      DOUBLE PRECISION ::    SZANGLES ( MAX_SZANGLES )
      INTEGER ::             N_USER_RELAZMS
      DOUBLE PRECISION ::    USER_RELAZMS ( MAX_USER_RELAZMS )
      INTEGER ::             N_USER_VZANGLES
      DOUBLE PRECISION ::    USER_VZANGLES ( MAX_USER_VZANGLES )

!  Exception handling

      INTEGER ::             STATUS_INPUTCHECK
      INTEGER ::             NMESSAGES
      CHARACTER (LEN=120) :: MESSAGES ( 0:MAX_MESSAGES )
      CHARACTER (LEN=120) :: ACTIONS ( 0:MAX_MESSAGES )

!  Other

      INTEGER          :: NM, I
      CHARACTER(Len=2) :: C2

!  ====================================
!  BEGIN COPY INPUTS TO LOCAL VARIABLES
!  ====================================

!  BRDF Control inputs

      BS_DO_USER_STREAMS     = VBRDF_Sup_In%BS_DO_USER_STREAMS
      BS_DO_BRDF_SURFACE     = VBRDF_Sup_In%BS_DO_BRDF_SURFACE
      BS_DO_SURFACE_EMISSION = VBRDF_Sup_In%BS_DO_SURFACE_EMISSION

!  BRDF Geometry inputs

      BS_NSTOKES             = VBRDF_Sup_In%BS_NSTOKES
      BS_NSTREAMS            = VBRDF_Sup_In%BS_NSTREAMS
      BS_NBEAMS              = VBRDF_Sup_In%BS_NBEAMS
      BS_BEAM_SZAS           = VBRDF_Sup_In%BS_BEAM_SZAS
      BS_N_USER_RELAZMS      = VBRDF_Sup_In%BS_N_USER_RELAZMS
      BS_USER_RELAZMS        = VBRDF_Sup_In%BS_USER_RELAZMS
      BS_N_USER_STREAMS      = VBRDF_Sup_In%BS_N_USER_STREAMS
      BS_USER_ANGLES_INPUT   = VBRDF_Sup_In%BS_USER_ANGLES_INPUT

!  VLIDORT Fixed Boolean inputs

      DO_SURFACE_EMISSION   = VLIDORT_FixIn%Bool%TS_DO_SURFACE_EMISSION
      DO_LAMBERTIAN_SURFACE = VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE

!  VLIDORT Modified Boolean inputs

      !DO_USER_STREAMS  = VLIDORT_ModIn%MBool%TS_DO_USER_STREAMS
      DO_USER_VZANGLES = VLIDORT_ModIn%MBool%TS_DO_USER_VZANGLES
      DO_MVOUT_ONLY    = VLIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY

!  VLIDORT Fixed Control inputs

      NSTOKES  = VLIDORT_FixIn%Cont%TS_NSTOKES
      NSTREAMS = VLIDORT_FixIn%Cont%TS_NSTREAMS

!  VLIDORT Fixed Beam inputs

      !BEAM_SZAS  = VLIDORT_ModIn%MSunRays%TS_BEAM_SZAS
      SZANGLES   = VLIDORT_ModIn%MSunRays%TS_SZANGLES

!  VLIDORT Modified Beam inputs

      !NBEAMS     = VLIDORT_ModIn%MSunRays%TS_NBEAMS
      N_SZANGLES = VLIDORT_ModIn%MSunRays%TS_N_SZANGLES

!  VLIDORT Fixed User Value inputs

      !N_USER_STREAMS  = VLIDORT_FixIn%UserVal%TS_N_USER_STREAMS
      N_USER_VZANGLES = VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES

!  VLIDORT Modified User Value inputs

      N_USER_RELAZMS  = VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS
      USER_RELAZMS    = VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS

      !USER_ANGLES     = VLIDORT_ModIn%MUserVal%TS_USER_ANGLES_INPUT
      USER_VZANGLES   = VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT

!  ==================================
!  END COPY INPUTS TO LOCAL VARIABLES
!  ==================================

!  Initialize output status

      STATUS_INPUTCHECK = VLIDORT_SUCCESS
      MESSAGES(1:MAX_MESSAGES) = ' '
      ACTIONS (1:MAX_MESSAGES) = ' '

      NMESSAGES   = 0
      MESSAGES(0) = 'Successful Check of BRDF/MAIN compatibility'
      ACTIONS(0)  = 'No Action required for this Task'

      NM = NMESSAGES

!  Checks

      IF ( BS_DO_BRDF_SURFACE.neqv.(.not.DO_LAMBERTIAN_SURFACE) ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'BRDF surface not set for VLIDORT Main'
        ACTIONS(NM)  = 'DO_LAMBERTIAN_SURFACE should be False!'
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
      ENDIF

      IF ( BS_DO_SURFACE_EMISSION.neqv.DO_SURFACE_EMISSION ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Surface emission flags do not agree'
        ACTIONS(NM)  = 'Check flag compatibility!'
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
      ENDIF

      IF ( BS_DO_USER_STREAMS .neqv. DO_USER_VZANGLES) THEN
        NM = NM + 1
        MESSAGES(NM) = 'User Streams/VZangles flags do not agree'
        ACTIONS(NM)  = 'Check flag compatibility!'
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
      ENDIF

      IF ( BS_DO_USER_STREAMS .neqv. (.not.DO_MVOUT_ONLY) ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'User Streams and DO_MVOUT_ONLY not agree'
        ACTIONS(NM)  = 'Check flag compatibility!'
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
      ENDIF

      IF ( BS_NSTREAMS .ne. NSTREAMS) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of discrete ordinates does not agree'
        ACTIONS(NM)  = 'Check NSTREAMS input'
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
      ENDIF

      IF ( BS_NSTOKES .ne. NSTOKES) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of Stokes components does not agree'
        ACTIONS(NM)  = 'Check NSTOKES input'
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
      ENDIF

!  Angles

      IF ( BS_NBEAMS .ne. N_SZANGLES) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of Solar beams does not agree'
        ACTIONS(NM)  = 'Check BS_NBEAMS and N_SZANGLES input'
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
      ELSE
        DO I = 1, N_SZANGLES
          if ( BS_BEAM_SZAS(I) .ne. SZANGLES(I) ) THEN
            write(C2,'(I2)')I
            NM = NM + 1
            MESSAGES(NM) = 'Solar beam angle does not agree, # '//C2
            ACTIONS(NM)  = 'Check BS_BEAM_SZAS and SZANGLES input'
            STATUS_INPUTCHECK = VLIDORT_SERIOUS
          endif
        ENDDO
      ENDIF

      IF ( BS_N_USER_STREAMS .ne. N_USER_VZANGLES) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of viewing zenith angles does not agree'
        ACTIONS(NM)  = 'Check N_USER_STREAMS and N_USER_VZANGLES input'
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
      ELSE
        DO I = 1, N_USER_VZANGLES
          if ( BS_USER_ANGLES_INPUT(I) .ne. USER_VZANGLES(I) ) THEN
            write(C2,'(I2)')I
            NM = NM + 1
            MESSAGES(NM) = 'View zenith angle does not agree, # '//C2
            ACTIONS(NM)  = 'Check BS_USER_ANGLES & USER_VZANGLES input'
            STATUS_INPUTCHECK = VLIDORT_SERIOUS
          endif
        ENDDO
      ENDIF

      IF ( BS_N_USER_RELAZMS .ne. N_USER_RELAZMS) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of viewing zenith angles does not agree'
        ACTIONS(NM)  = 'Check BS_N_USER_RELAZMS & N_USER_RELAZMS input'
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
      ELSE
        DO I = 1, N_USER_RELAZMS
          if ( BS_USER_RELAZMS(I) .ne.USER_RELAZMS(I) ) THEN
            write(C2,'(I2)')I
            NM = NM + 1
            MESSAGES(NM) = 'Azimuth angle does not agree, # '//C2
            ACTIONS(NM)  = 'Check BS_USER_RELAZMS & USER_RELAZMS input'
            STATUS_INPUTCHECK = VLIDORT_SERIOUS
          endif
        ENDDO
      ENDIF

!  Tally up messages

      NMESSAGES = NM

!  Copy Exception handling output

      VLIDORT_VBRDFCheck_Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
      VLIDORT_VBRDFCheck_Status%TS_NCHECKMESSAGES    = NMESSAGES
      VLIDORT_VBRDFCheck_Status%TS_CHECKMESSAGES     = MESSAGES
      VLIDORT_VBRDFCheck_Status%TS_ACTIONS           = ACTIONS

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_VBRDF_INPUT_CHECKER

!

      SUBROUTINE VLIDORT_VSLEAVE_INPUT_CHECKER ( &
        VSLEAVE_Sup_In,             & ! Inputs
        VLIDORT_FixIn,              & ! Inputs
        VLIDORT_ModIn,              & ! Inputs
        VLIDORT_VSLEAVECheck_Status ) ! Outputs

      USE VLIDORT_PARS

      USE VSLEAVE_Sup_Inputs_def
      USE VLIDORT_Inputs_def
      USE VLIDORT_Outputs_def

      IMPLICIT NONE

      TYPE(VSLEAVE_Sup_Inputs), INTENT(IN)           :: VSLEAVE_Sup_In

      TYPE(VLIDORT_Fixed_Inputs), INTENT (IN)        :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs), INTENT (IN)     :: VLIDORT_ModIn

      TYPE(VLIDORT_Exception_Handling), INTENT(OUT)  :: VLIDORT_VSLEAVECheck_Status

!  ---------------
!  Local variables
!  ---------------

!  VSLEAVE supplement inputs
!  -------------------------

!  Surface-leaving control flags

      LOGICAL :: SL_DO_SLEAVING
      LOGICAL :: SL_DO_ISOTROPIC
      LOGICAL :: SL_DO_EXACTONLY
      LOGICAL :: SL_DO_USER_STREAMS

!  Number of Stokes vector components

      INTEGER ::             SL_NSTOKES

!  Number of discrete ordinate streams

      INTEGER ::             SL_NSTREAMS

!  BOA solar zenith angles

      INTEGER ::             SL_NBEAMS
      DOUBLE PRECISION ::    SL_BEAM_SZAS ( MAXBEAMS )

!  User-defined relative azimuths (mandatory for Fourier > 0)

      INTEGER ::             SL_N_USER_RELAZMS
      DOUBLE PRECISION ::    SL_USER_RELAZMS (MAX_USER_RELAZMS)

!  User-defined zenith angle input

      INTEGER ::             SL_N_USER_STREAMS
      DOUBLE PRECISION ::    SL_USER_ANGLES_INPUT (MAX_USER_STREAMS)

!  VLIDORT Main inputs
!  -------------------

      LOGICAL ::             DO_USER_VZANGLES
      LOGICAL ::             DO_MVOUT_ONLY
      LOGICAL ::             DO_SURFACE_LEAVING
      LOGICAL ::             DO_SL_ISOTROPIC

      INTEGER ::             NSTOKES
      INTEGER ::             NSTREAMS

      INTEGER ::             N_SZANGLES
      DOUBLE PRECISION ::    SZANGLES ( MAX_SZANGLES )
      INTEGER ::             N_USER_RELAZMS
      DOUBLE PRECISION ::    USER_RELAZMS ( MAX_USER_RELAZMS )
      INTEGER ::             N_USER_VZANGLES
      DOUBLE PRECISION ::    USER_VZANGLES ( MAX_USER_VZANGLES )

!  Exception handling

      INTEGER ::             STATUS_INPUTCHECK
      INTEGER ::             NMESSAGES
      CHARACTER (LEN=120) :: MESSAGES ( 0:MAX_MESSAGES )
      CHARACTER (LEN=120) :: ACTIONS ( 0:MAX_MESSAGES )

!  Other

      INTEGER          :: NM, I
      CHARACTER(Len=2) :: C2

!  ====================================
!  BEGIN COPY INPUTS TO LOCAL VARIABLES
!  ====================================

!  VSLEAVE Control inputs

      SL_DO_SLEAVING         = VSLEAVE_Sup_In%SL_DO_SLEAVING
      SL_DO_ISOTROPIC        = VSLEAVE_Sup_In%SL_DO_ISOTROPIC
      SL_DO_EXACTONLY        = VSLEAVE_Sup_In%SL_DO_EXACTONLY
      SL_DO_USER_STREAMS     = VSLEAVE_Sup_In%SL_DO_USER_STREAMS

!  VSLEAVE Geometry inputs

      SL_NSTOKES             = VSLEAVE_Sup_In%SL_NSTOKES
      SL_NSTREAMS            = VSLEAVE_Sup_In%SL_NSTREAMS
      SL_NBEAMS              = VSLEAVE_Sup_In%SL_NBEAMS
      SL_BEAM_SZAS           = VSLEAVE_Sup_In%SL_BEAM_SZAS
      SL_N_USER_RELAZMS      = VSLEAVE_Sup_In%SL_N_USER_RELAZMS
      SL_USER_RELAZMS        = VSLEAVE_Sup_In%SL_USER_RELAZMS
      SL_N_USER_STREAMS      = VSLEAVE_Sup_In%SL_N_USER_STREAMS
      SL_USER_ANGLES_INPUT   = VSLEAVE_Sup_In%SL_USER_ANGLES_INPUT

!  VLIDORT Fixed Boolean inputs

      DO_SURFACE_LEAVING = VLIDORT_FixIn%Bool%TS_DO_SURFACE_LEAVING
      DO_SL_ISOTROPIC    = VLIDORT_FixIn%Bool%TS_DO_SL_ISOTROPIC

!  VLIDORT Modified Boolean inputs

      !DO_USER_STREAMS  = VLIDORT_ModIn%MBool%TS_DO_USER_STREAMS
      DO_USER_VZANGLES = VLIDORT_ModIn%MBool%TS_DO_USER_VZANGLES
      DO_MVOUT_ONLY    = VLIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY

!  VLIDORT Fixed Control inputs

      NSTOKES  = VLIDORT_FixIn%Cont%TS_NSTOKES
      NSTREAMS = VLIDORT_FixIn%Cont%TS_NSTREAMS

!  VLIDORT Fixed Beam inputs

      !BEAM_SZAS  = VLIDORT_ModIn%MSunRays%TS_BEAM_SZAS
      SZANGLES   = VLIDORT_ModIn%MSunRays%TS_SZANGLES

!  VLIDORT Modified Beam inputs

      !NBEAMS     = VLIDORT_ModIn%MSunRays%TS_NBEAMS
      N_SZANGLES = VLIDORT_ModIn%MSunRays%TS_N_SZANGLES

!  VLIDORT Fixed User Value inputs

      !N_USER_STREAMS  = VLIDORT_FixIn%UserVal%TS_N_USER_STREAMS
      N_USER_VZANGLES = VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES

!  VLIDORT Modified User Value inputs

      N_USER_RELAZMS  = VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS
      USER_RELAZMS    = VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS

      !USER_ANGLES     = VLIDORT_ModIn%MUserVal%TS_USER_ANGLES_INPUT
      USER_VZANGLES   = VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT

!  ==================================
!  END COPY INPUTS TO LOCAL VARIABLES
!  ==================================

!  Initialize output status

      STATUS_INPUTCHECK = VLIDORT_SUCCESS
      MESSAGES(1:MAX_MESSAGES) = ' '
      ACTIONS (1:MAX_MESSAGES) = ' '

      NMESSAGES   = 0
      MESSAGES(0) = 'Successful Check of VSLEAVE/MAIN compatibility'
      ACTIONS(0)  = 'No Action required for this Task'

      NM = NMESSAGES

!  Checks

      IF ( SL_DO_SLEAVING.neqv.DO_SURFACE_LEAVING ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Surface-leaving control flags do not agree'
        ACTIONS(NM)  = 'Check flag compatibility!'
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
      ENDIF

      IF ( SL_DO_ISOTROPIC.neqv.DO_SL_ISOTROPIC ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Surface-leaving isotropic flags do not agree'
        ACTIONS(NM)  = 'Check flag compatibility!'
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
      ENDIF

      IF ( SL_DO_USER_STREAMS .neqv. DO_USER_VZANGLES) THEN
        NM = NM + 1
        MESSAGES(NM) = 'User Streams/VZangles flags do not agree'
        ACTIONS(NM)  = 'Check flag compatibility!'
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
      ENDIF

      IF ( SL_DO_USER_STREAMS .neqv. (.not.DO_MVOUT_ONLY) ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'User Streams and DO_MVOUT_ONLY not agree'
        ACTIONS(NM)  = 'Check flag compatibility!'
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
      ENDIF

      IF ( SL_NSTREAMS .ne. NSTREAMS) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of discrete ordinates does not agree'
        ACTIONS(NM)  = 'Check SL_NSTREAMS and NSTREAMS input'
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
      ENDIF

      IF ( SL_NSTOKES .ne. NSTOKES) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of Stokes components does not agree'
        ACTIONS(NM)  = 'Check SL_NSTOKES and NSTOKES input'
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
      ENDIF

!  Angles

      IF ( SL_NBEAMS .ne. N_SZANGLES) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of Solar beams does not agree'
        ACTIONS(NM)  = 'Check SL_NBEAMS and N_SZANGLES input'
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
      ELSE
        DO I = 1, N_SZANGLES
          if ( SL_BEAM_SZAS(I) .ne. SZANGLES(I) ) THEN
            write(C2,'(I2)')I
            NM = NM + 1
            MESSAGES(NM) = 'Solar beam angle does not agree, # '//C2
            ACTIONS(NM)  = 'Check SL_BEAM_SZAS and SZANGLES input'
            STATUS_INPUTCHECK = VLIDORT_SERIOUS
          endif
        ENDDO
      ENDIF

      IF ( SL_N_USER_STREAMS .ne. N_USER_VZANGLES) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of viewing zenith angles does not agree'
        ACTIONS(NM)  = 'Check SL_N_USER_STREAMS and N_USER_VZANGLES input'
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
      ELSE
        DO I = 1, N_USER_VZANGLES
          if ( SL_USER_ANGLES_INPUT(I) .ne. USER_VZANGLES(I) ) THEN
            write(C2,'(I2)')I
            NM = NM + 1
            MESSAGES(NM) = 'View zenith angle does not agree, # '//C2
            ACTIONS(NM)  = 'Check SL_USER_ANGLES_INPUT & USER_VZANGLES input'
            STATUS_INPUTCHECK = VLIDORT_SERIOUS
          endif
        ENDDO
      ENDIF

      IF ( SL_N_USER_RELAZMS .ne. N_USER_RELAZMS) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of viewing zenith angles does not agree'
        ACTIONS(NM)  = 'Check SL_N_USER_RELAZMS & N_USER_RELAZMS input'
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
      ELSE
        DO I = 1, N_USER_RELAZMS
          if ( SL_USER_RELAZMS(I) .ne.USER_RELAZMS(I) ) THEN
            write(C2,'(I2)')I
            NM = NM + 1
            MESSAGES(NM) = 'Azimuth angle does not agree, # '//C2
            ACTIONS(NM)  = 'Check SL_USER_RELAZMS & USER_RELAZMS input'
            STATUS_INPUTCHECK = VLIDORT_SERIOUS
          endif
        ENDDO
      ENDIF

!  Tally up messages

      NMESSAGES = NM

!  Copy Exception handling output

      VLIDORT_VSLEAVECheck_Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
      VLIDORT_VSLEAVECheck_Status%TS_NCHECKMESSAGES    = NMESSAGES
      VLIDORT_VSLEAVECheck_Status%TS_CHECKMESSAGES     = MESSAGES
      VLIDORT_VSLEAVECheck_Status%TS_ACTIONS           = ACTIONS

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_VSLEAVE_INPUT_CHECKER

!

      SUBROUTINE VBRDF_VSLEAVE_INPUT_CHECKER ( &
        VSLEAVE_Sup_In,             & ! Inputs
        VBRDF_Sup_In,               & ! Inputs
        VBRDF_VSLEAVECheck_Status )   ! Outputs

      USE VLIDORT_PARS

      USE VSLEAVE_Sup_Inputs_def
      USE VBRDF_Sup_Inputs_def
      USE VLIDORT_Outputs_def

      IMPLICIT NONE

      TYPE(VSLEAVE_Sup_Inputs), INTENT(IN)           :: VSLEAVE_Sup_In
      TYPE(VBRDF_Sup_Inputs), INTENT(IN)             :: VBRDF_Sup_In

      TYPE(VLIDORT_Exception_Handling), INTENT(OUT)  :: VBRDF_VSLEAVECheck_Status

!  ---------------
!  Local variables
!  ---------------

!  VSLEAVE supplement inputs
!  -------------------------

!  Surface-leaving control flags #1 (used for gatekeeping checks)

      LOGICAL :: SL_DO_ISOTROPIC
      LOGICAL :: SL_DO_FLUORESCENCE

!  Number of solar zenith angles

      INTEGER :: SL_NBEAMS

!  Surface-leaving control flags #2

      LOGICAL :: SL_DO_GlintShadow
      LOGICAL :: SL_DO_FoamOption
      LOGICAL :: SL_DO_FacetIsotropy

!  Salinity

      DOUBLE PRECISION :: SL_SALINITY

!  Wind-speed and directions

      DOUBLE PRECISION :: SL_WINDSPEED
      DOUBLE PRECISION :: SL_WINDDIR ( MAXBEAMS )

!  VBRDF supplement inputs
!  -----------------------

!  Number of solar zenith angles

      INTEGER :: BS_NBEAMS

!  Surface-leaving control flags

      LOGICAL :: BS_DO_GlintShadow
      LOGICAL :: BS_DO_FoamOption
      LOGICAL :: BS_DO_FacetIsotropy

!  Salinity

      DOUBLE PRECISION :: BS_SALINITY

!  Wind-speed and directions

      DOUBLE PRECISION :: BS_WINDSPEED
      DOUBLE PRECISION :: BS_WINDDIR ( MAXBEAMS )

!  Exception handling

      INTEGER ::             STATUS_INPUTCHECK
      INTEGER ::             NMESSAGES
      CHARACTER (LEN=120) :: MESSAGES ( 0:MAX_MESSAGES )
      CHARACTER (LEN=120) :: ACTIONS ( 0:MAX_MESSAGES )

!  Other

      INTEGER          :: NM, I
      CHARACTER(Len=2) :: C2

!  ====================================
!  BEGIN COPY INPUTS TO LOCAL VARIABLES
!  ====================================

!  VSLEAVE Control inputs #1 (used for gatekeeping checks)

      SL_DO_ISOTROPIC     = VSLEAVE_Sup_In%SL_DO_ISOTROPIC
      SL_DO_FLUORESCENCE  = VSLEAVE_Sup_In%SL_DO_FLUORESCENCE

!  VSLEAVE Geometry inputs

      SL_NBEAMS           = VSLEAVE_Sup_In%SL_NBEAMS

!  VSLEAVE Control inputs #2

      SL_DO_GlintShadow   = VSLEAVE_Sup_In%SL_DO_GlintShadow
      SL_DO_FoamOption    = VSLEAVE_Sup_In%SL_DO_FoamOption
      SL_DO_FacetIsotropy = VSLEAVE_Sup_In%SL_DO_FacetIsotropy

!  VSLEAVE Other inputs

      SL_SALINITY         = VSLEAVE_Sup_In%SL_SALINITY
      SL_WINDSPEED        = VSLEAVE_Sup_In%SL_WINDSPEED
      SL_WINDDIR          = VSLEAVE_Sup_In%SL_WINDDIR

!  VBRDF Geometry inputs

      BS_NBEAMS           = VBRDF_Sup_In%BS_NBEAMS

!  VBRDF Control inputs

      BS_DO_GlintShadow   = VBRDF_Sup_In%BS_DO_GlintShadow
      BS_DO_FoamOption    = VBRDF_Sup_In%BS_DO_FoamOption
      BS_DO_FacetIsotropy = VBRDF_Sup_In%BS_DO_FacetIsotropy

!  VBRDF Other inputs

      BS_SALINITY         = VBRDF_Sup_In%BS_SALINITY
      BS_WINDSPEED        = VBRDF_Sup_In%BS_WINDSPEED
      BS_WINDDIR          = VBRDF_Sup_In%BS_WINDDIR

!mick debug
!write(*,*)
!write(*,*) 'mick debug'
!write(*,*) 'SL_DO_GlintShadow   = ',SL_DO_GlintShadow
!write(*,*) 'SL_DO_FoamOption    = ',SL_DO_FoamOption
!write(*,*) 'SL_DO_FacetIsotropy = ',SL_DO_FacetIsotropy
!write(*,*) 'SL_SALINITY  = ',SL_SALINITY
!write(*,*) 'SL_WINDSPEED = ',SL_WINDSPEED
!write(*,*) 'SL_WINDDIR   = ',SL_WINDDIR
!write(*,*)
!write(*,*) 'BS_DO_GlintShadow   = ',BS_DO_GlintShadow
!write(*,*) 'BS_DO_FoamOption    = ',BS_DO_FoamOption
!write(*,*) 'BS_DO_FacetIsotropy = ',BS_DO_FacetIsotropy
!write(*,*) 'BS_SALINITY  = ',BS_SALINITY
!write(*,*) 'BS_WINDSPEED = ',BS_WINDSPEED
!write(*,*) 'BS_WINDDIR   = ',BS_WINDDIR

!  ==================================
!  END COPY INPUTS TO LOCAL VARIABLES
!  ==================================

!  Initialize output status

      STATUS_INPUTCHECK = VLIDORT_SUCCESS
      MESSAGES(1:MAX_MESSAGES) = ' '
      ACTIONS (1:MAX_MESSAGES) = ' '

      NMESSAGES   = 0
      MESSAGES(0) = 'Successful Check of VSLEAVE/VBRDF compatibility'
      ACTIONS(0)  = 'No Action required for this Task'

      NM = NMESSAGES

!  Checks
!  ------

!  NewCMGLINT/VSLEAVE gatekeeper checks
!  (if any fails, don't look at inputs further until correct)

      IF ( SL_DO_ISOTROPIC ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'New Cox-Munk glint BRDF is active and surface-leaving DO_ISOTROPIC is active'
        ACTIONS(NM)  = 'Deactivate surface-leaving isotropy!'
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
      ENDIF

      IF ( SL_DO_FLUORESCENCE ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'New Cox-Munk glint BRDF is active and surface-leaving DO_FLUORESCENCE is active'
        ACTIONS(NM)  = 'Deactivate surface-leaving fluorescence!'
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
      ENDIF

      IF ( STATUS_INPUTCHECK == VLIDORT_SERIOUS ) GOTO 500

!  Control flags

      IF ( SL_DO_GlintShadow.neqv.BS_DO_GlintShadow ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Surface-leaving glint shadow flag does not agree'
        ACTIONS(NM)  = 'Check flag compatibility!'
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
      ENDIF

      IF ( SL_DO_FoamOption.neqv.BS_DO_FoamOption ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Surface-leaving foam option flag does not agree'
        ACTIONS(NM)  = 'Check flag compatibility!'
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
      ENDIF

      IF ( SL_DO_FacetIsotropy.neqv.BS_DO_FacetIsotropy ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Surface-leaving facet isotropy flag does not agree'
        ACTIONS(NM)  = 'Check flag compatibility!'
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
      ENDIF

!  Salinity

      IF ( SL_SALINITY .ne. BS_SALINITY) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Surface-leaving salinity does not agree'
        ACTIONS(NM)  = 'Check SL_SALINITY and BS_SALINITY input'
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
      ENDIF

!  Wind-speed and directions

      IF ( SL_WINDSPEED .ne. BS_WINDSPEED) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Surface-leaving wind-speed does not agree'
        ACTIONS(NM)  = 'Check SL_WINDSPEED and BS_WINDSPEED input'
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
      ENDIF

      DO I = 1, BS_NBEAMS
        if ( SL_WINDDIR(I) .ne. BS_WINDDIR(I) ) THEN
          write(C2,'(I2)')I
          NM = NM + 1
          MESSAGES(NM) = 'Wind direction angle does not agree, # '//C2
          ACTIONS(NM)  = 'Check SL_WINDDIR and BS_WINDDIR input'
          STATUS_INPUTCHECK = VLIDORT_SERIOUS
        endif
      ENDDO

500   CONTINUE

!  Tally up messages

      NMESSAGES = NM

!  Copy Exception handling output

      VBRDF_VSLEAVECheck_Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
      VBRDF_VSLEAVECheck_Status%TS_NCHECKMESSAGES    = NMESSAGES
      VBRDF_VSLEAVECheck_Status%TS_CHECKMESSAGES     = MESSAGES
      VBRDF_VSLEAVECheck_Status%TS_ACTIONS           = ACTIONS

!  Finish

      RETURN
      END SUBROUTINE VBRDF_VSLEAVE_INPUT_CHECKER

!  Finish Module

      END MODULE vlidort_sup_accessories

