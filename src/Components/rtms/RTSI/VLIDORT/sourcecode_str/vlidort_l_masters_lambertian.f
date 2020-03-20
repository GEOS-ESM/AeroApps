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
C #            VLIDORT_L_MASTER_LAMBERTIAN (master)             #
C #            VLIDORT_L_FOURIER_LAMBERTIAN (master)            #
C #                                                             #
C ###############################################################

      SUBROUTINE VLIDORT_L_MASTER_LAMBERTIAN
     &     ( STATUS_INPUTCHECK, STATUS_CALCULATION )

C  Standard model include files
C  ----------------------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  Include file of linearization input variables

      INCLUDE '../includes/VLIDORT_L_INPUTS.VARS'

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

C  Standard input variables check

      CALL VLIDORT_CHECK_INPUT ( STATUS_SUB )

      IF ( STATUS_SUB .EQ. VLIDORT_SERIOUS ) THEN
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
        MAIL = ' VLIDORT_CHECK_INPUT failed'
        TRACE = ' Called in VLIDORT_L_MASTER '
        CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS_INPUTCHECK )
        CLOSE(VLIDORT_ERRUNIT)
        RETURN
      ELSE IF ( STATUS_SUB .EQ. VLIDORT_WARNING ) THEN
        STATUS_INPUTCHECK = VLIDORT_WARNING
        MAIL = ' VLIDORT_CHECK_INPUT warning message'
        TRACE = ' Called in VLIDORT_L_MASTER '
        CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS_INPUTCHECK )
      ENDIF

C  Extended input variables check

      IF ( DO_LINEARIZATION ) THEN

        CALL VLIDORT_L_CHECK_INPUT ( STATUS_SUB )

        IF ( STATUS_SUB .EQ. VLIDORT_SERIOUS ) THEN
          STATUS_INPUTCHECK = VLIDORT_SERIOUS
          MAIL = ' VLIDORT_L_CHECK_INPUT failed'
          TRACE = ' Called in VLIDORT_L_MASTER_LAMBERTIAN'
          CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS_INPUTCHECK )
          CLOSE(VLIDORT_ERRUNIT)
          RETURN
        ELSE IF ( STATUS_SUB .EQ. VLIDORT_WARNING ) THEN
          STATUS_INPUTCHECK = VLIDORT_WARNING
          MAIL = ' VLIDORT_L_CHECK_INPUT warning message'
          TRACE = ' Called in VLIDORT_L_MASTER_LAMBERTIAN'
          CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS_INPUTCHECK )
        ENDIF

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

      modified_eradius = earth_radius + geometry_specheight
      CALL multi_outgoing_adjustgeom
     i   ( MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS,
     i     N_USER_VZANGLES,   N_SZANGLES,   N_USER_RELAZMS,
     i     height_grid(nlayers), modified_eradius, adjust_surface,
     i     user_vzangles,  szangles, user_relazms,
     o     user_vzangles_adjust, szangles_adjust, user_relazms_adjust,
     o     fail, mail )
      if ( fail ) return       

C  Chapman function calculation
C  ----------------------------

      IF ( DO_CHAPMAN_FUNCTION ) THEN
        CALL VLIDORT_CHAPMAN ( FAIL, MAIL, TRACE )
        IF (FAIL) THEN
          STATUS_CALCULATION = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS_CALCULATION )
          CLOSE(VLIDORT_ERRUNIT)
          RETURN
        ENDIF
      ENDIF

C  write input variables
C  ---------------------

C  open file, call standard/linearized input writes, close file

      IF ( DO_WRITE_INPUT ) THEN
        IUNIT = VLIDORT_INUNIT
        OPEN(IUNIT,FILE=INPUT_WRITE_FILENAME,STATUS='UNKNOWN')
        CALL VLIDORT_WRITEINPUT ( IUNIT )
        IF ( DO_LINEARIZATION ) THEN
          CALL VLIDORT_L_WRITEINPUT ( IUNIT )
        ENDIF
        CLOSE(IUNIT)
      ENDIF

C  Get derived inputs
C  ==================

C  Miscellaneous and layer input.

      CALL VLIDORT_DERIVE_INPUT ( STATUS_SUB, MAIL )

      IF ( STATUS_SUB .EQ. VLIDORT_WARNING ) THEN
        STATUS_INPUTCHECK = VLIDORT_WARNING
        TRACE = 'Derive_Input Call in VLIDORT_L_MASTERS_LAMBERTIAN'
        CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS_INPUTCHECK )
      ENDIF

C  Derive number of Surface weighting functions
C  Automatically 1 in this version without the BRDFs

      N_TOTALBRDF_WFS = 1

c      IF ( .NOT. DO_LAMBERTIAN_SURFACE ) THEN
c       COUNT = 0
c       DO KL = 1, N_BRDF_KERNELS
c        IF ( DO_KERNEL_FACTOR_WFS(KL) ) COUNT = COUNT + 1
c        DO_KPARAMS_DERIVS(KL) = .FALSE.
c        DO KQ = 1, N_BRDF_PARAMETERS(KL)
c          IF ( DO_KERNEL_PARAMS_WFS(KL,KQ) ) THEN
c            COUNT = COUNT + 1
c            DO_KPARAMS_DERIVS(KL) = .TRUE.
c          ENDIF
c        ENDDO
c       ENDDO
c      ENDIF
C      N_TOTALBRDF_WFS = COUNT

C  standard + linearized BRDF inputs, both together

c      IF ( .NOT. DO_LAMBERTIAN_SURFACE ) THEN
c        CALL VLIDORT_L_BRDF_CREATOR
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
C   Other cases now disabled, 17 Janaury 2006

C      IF ( DO_TRANSMITTANCE_ONLY ) THEN
C        LOCAL_DO_NO_AZIMUTH = .TRUE.
C      ENDIF

c      IF ( DO_ISOTROPIC_ONLY  ) THEN
c        LOCAL_DO_NO_AZIMUTH = .TRUE.
c      ENDIF

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
        NSOURCES = NBEAMS
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

C  Initialize

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

      DO WHILE ( LOCAL_ITERATION .AND. 
     &             FOURIER_COMPONENT.LT.N_FOURIER_COMPONENTS )

C  Fourier counter

        FOURIER_COMPONENT = FOURIER_COMPONENT + 1

C  Local start of user-defined streams. Should always be 1.
C    No zenith tolerance now.

        LOCAL_UM_START = 1

C  azimuth cosine/sine factors, using adjust geometries.

        IF ( FOURIER_COMPONENT .GT. 0 ) THEN
          DFC = DBLE(FOURIER_COMPONENT)
          DO UA = 1, LOCAL_N_USERAZM
            DO IB = 1, NBEAMS
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

        CALL VLIDORT_L_FOURIER_LAMBERTIAN
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
          MAIL = 'Error from VLIDORT_L_FOURIER_LAMBERTIAN, Fourier '//CF
          TRACE = 'Called by VLIDORT_L_MASTER_LAMBERTIAN'
          CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS_CALCULATION )
          CLOSE(VLIDORT_ERRUNIT)
          RETURN
        ENDIF

C  Fourier summation and Convergence examination
C  ---------------------------------------------

C   -- only done for beams which are still not converged.
C      This is controlled by flag DO_MULTIBEAM

C   -- new criterion, SS is added for Fourier = 0, as this means that
C      higher-order terms will be relatively smaller, which implies
C      faster convergence in some circumstances (generally not often).

        IBEAM_COUNT = 0
        DO IBEAM = 1, NBEAMS
          IF ( DO_MULTIBEAM ( IBEAM, FOURIER_COMPONENT ) ) THEN

C  Convergence and radiance summation

            CALL VLIDORT_CONVERGE
     I      ( AZMFAC, FOURIER_COMPONENT, LOCAL_N_USERAZM,
     O        BEAM_TESTCONV(IBEAM), IBEAM, BEAM_ITERATION(IBEAM) )

C  Fourier summation of linearization quantities

            IF ( DO_LINEARIZATION ) THEN
              CALL VLIDORT_L_CONVERGE
     &      ( AZMFAC, FOURIER_COMPONENT, LOCAL_N_USERAZM, IBEAM )
            ENDIF

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
C  Write Standard   Fourier output (radiances)
C  Write additional Fourier output (linearizations)
C  Close file if iteration has finished
C  New comment:
C    If the SS correction is set, Fourier=0 will include SS field

        IF ( DO_WRITE_FOURIER ) THEN
          FUNIT = VLIDORT_FUNIT
          IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
            OPEN(FUNIT,FILE=FOURIER_WRITE_FILENAME,STATUS='UNKNOWN')
          ENDIF
          CALL VLIDORT_WRITEFOURIER ( FUNIT, FOURIER_COMPONENT )
          IF ( DO_LINEARIZATION ) THEN
            CALL VLIDORT_L_WRITEFOURIER
     I         ( DO_INCLUDE_SURFACE, FUNIT, FOURIER_COMPONENT )
          ENDIF
          IF ( .NOT.LOCAL_ITERATION ) CLOSE ( FUNIT )
        ENDIF

C  end iteration loop

      ENDDO

C  restore no azimuth flag

      DO_NO_AZIMUTH = SAVE_DO_NO_AZIMUTH

C  Major result output
C  ===================

C  Standard + linearized output
C  ----------------------------

      IF ( DO_WRITE_RESULTS ) THEN
        RUNIT = VLIDORT_RESUNIT
        OPEN(RUNIT,FILE=RESULTS_WRITE_FILENAME,STATUS='UNKNOWN')
        CALL VLIDORT_WRITERESULTS ( RUNIT )
        IF ( DO_LINEARIZATION ) THEN
          CALL VLIDORT_L_WRITERESULTS ( RUNIT )
        ENDIF
        CLOSE(RUNIT)
      ENDIF

C  Geophysical input (scenario) write
C  ----------------------------------

      IF ( DO_WRITE_SCENARIO ) THEN
        SUNIT = VLIDORT_SCENUNIT
        OPEN(SUNIT,FILE=SCENARIO_WRITE_FILENAME,STATUS='UNKNOWN')
        CALL VLIDORT_WRITESCEN ( SUNIT )
        IF ( DO_LINEARIZATION ) THEN
          CALL VLIDORT_L_WRITESCEN ( SUNIT )
        ENDIF
        CLOSE(SUNIT)
      ENDIF

C  Closing business
C  ================

C  close Error file if it was used

      IF ( .NOT. VLIDORT_ERROR_INIT ) CLOSE(VLIDORT_ERRUNIT)

C  Finish

      RETURN
      END

C

      SUBROUTINE VLIDORT_L_FOURIER_LAMBERTIAN
     I       ( FOURIER_COMPONENT,
     I         SS_FLUX_MULTIPLIER,
     O         DO_INCLUDE_SURFACE,
     O         DO_INCLUDE_SURFEMISS,
     O         DO_INCLUDE_THERMEMISS,
     O         STATUS )

C  Complete Fourier component calculation for the Extended Code
C    Stokes vector and Weighting function computations

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include file of input weighting function control variables

      INCLUDE '../includes/VLIDORT_L_INPUTS.VARS'

C  Include file of setup variables (output from SETUPS routine)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_THERMALSUP.VARS'

C  include this file file for debug only

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

      EXTERNAL         L_BVP_LAMBERTIAN_SURFACE
      EXTERNAL         L_BOA_LAMBERTIAN_SOURCE

      EXTERNAL         LS_BVP_LAMBERTIAN_SURFACE
      EXTERNAL         LS_BOA_LAMBERTIAN_SOURCE

C  weighting function indices

      INTEGER          N_LAYER_WFS, LAYER_TO_VARY, N_PARAMETERS
      INTEGER          BRDF_INDEX, WF_INDEX
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
C   (DO_REFLECTED_DIRECTBEAM is a stored variable)

      IF ( DO_DIRECT_BEAM ) THEN
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

      FLUX_MULTIPLIER    = SS_FLUX_MULTIPLIER * DELTA_FACTOR
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
C   MISCSETUPS   : Delta-M, average-secant formulation, transmittances
C   EMULT_MASTER : Beam source function multipliers. Not required for the
C                  Full SS calculation in outgoing mode

C  ############################ IMPORTANT NOTE, READ CAREFULLY #########
C
C  30 January 2008
C  Telescoping setup is done in DERIVE_INPUTS (unlike LIDORT scalar code)
C 
C  ############################ IMPORTANT NOTE, READ CAREFULLY #########

      IF ( FOURIER_COMPONENT .EQ. 0 ) THEN

C  With linearization.

        IF ( DO_ATMOS_LINEARIZATION ) THEN

          CALL VLIDORT_MISCSETUPS 
          CALL VLIDORT_L_MISCSETUPS
          IF ( DO_INCLUDE_THERMEMISS ) THEN
            CALL THERMAL_SETUP_PLUS
          ENDIF
          IF ( DO_SOLAR_SOURCES      ) THEN
            IF (.NOT.DO_SSFULL.OR.(DO_SSFULL.AND.DO_SSCORR_NADIR)) THEN
              CALL EMULT_MASTER
              CALL L_EMULT_MASTER
            ENDIF
          ENDIF

C  No linearization

        ELSE

          CALL VLIDORT_MISCSETUPS
          IF ( DO_INCLUDE_THERMEMISS ) THEN
            CALL THERMAL_SETUP
          ENDIF
          IF ( DO_SOLAR_SOURCES      ) THEN
            IF (.NOT.DO_SSFULL.OR.(DO_SSFULL.AND.DO_SSCORR_NADIR)) THEN
              CALL EMULT_MASTER
            ENDIF
          ENDIF

        ENDIF

C  End Fourier 0 clause

      ENDIF

C  #####################
C  Correction operations 
C  #####################

C  Not required if no solar sources

      IF ( DO_SOLAR_SOURCES ) THEN

C  Single scatter correction (pre-calculation)
C  -------------------------------------------

C   Must be done after MISCSETUPS, as we need the multipliers and
C    transmittance factors for the SUN and LOS paths.
C   Standard Code added 6 May 2005. Replaces call in Master routine.
C   Linearized code 13 July 2005.  Replaces former code in L_Master module.

C      Version 2.2. Added call to the new outgoing sphericity correction
C      Version 2.4. Added call bulk/column linearizations

       IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
        IF ( DO_USER_STREAMS ) THEN

C  Regular (nadir view) SS correction

          IF ( DO_SSCORR_NADIR ) THEN
            CALL VLIDORT_SSCORR_NADIR ( SS_FLUX_MULTIPLIER )
            IF ( DO_PROFILE_LINEARIZATION ) THEN
              DO LAYER_TO_VARY = 1, NLAYERS
                IF ( LAYER_VARY_FLAG(LAYER_TO_VARY) ) THEN
                  N_PARAMETERS = LAYER_VARY_NUMBER(LAYER_TO_VARY)
                  CALL VLIDORT_LP_SSCORR_NADIR
     &            ( SS_FLUX_MULTIPLIER, LAYER_TO_VARY, N_PARAMETERS)
                ENDIF
              ENDDO
            ENDIF
            IF ( DO_COLUMN_LINEARIZATION ) THEN
              CALL VLIDORT_LC_SSCORR_NADIR
     &            ( SS_FLUX_MULTIPLIER, N_TOTALCOLUMN_WFS)
            ENDIF
          ENDIF

C  Outgoing sphericity correction ------ Added 20 March 2007.
C     Calling not changed for column Jacobians (Version 2.4, December 2008)

          IF ( DO_SSCORR_OUTGOING ) THEN
           IF ( DO_ATMOS_LINEARIZATION ) THEN
            CALL VLIDORT_L_SSCORR_OUTGOING
     &     ( SS_FLUX_MULTIPLIER, FAIL, MAIL )
            IF ( FAIL ) THEN
              TRACE= 'Linearized SS correction outgoing failed'
              STATUS = VLIDORT_SERIOUS
              CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
              RETURN
            ENDIF
           ELSE
            CALL VLIDORT_SSCORR_OUTGOING
     &     ( SS_FLUX_MULTIPLIER, FAIL, MAIL )
            IF ( FAIL ) THEN
              TRACE= 'SS correction outgoing failed'
              STATUS = VLIDORT_SERIOUS
              CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
              RETURN
            ENDIF
           ENDIF
          ENDIF

C  End fourier = 0, user stream and solar source  clauses

        ENDIF
       ENDIF
      ENDIF

C  Exact direct beam BRDF (if flagged)
C  -----------------------------------

C  THIS IS THE LAMBERTIAN MASTER, SO CODE NOT REQUIRED

C   Must be done after MISCSETUPS, as we need the SUN transmittance factors
C   BRDF Kernels already known (called before this module)
C   Code added 6 May 2005. Standard case.
C   Profile/surface Linearization code added July 13 2005.
C     Column        Linearization code added December 2008

c      IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
c        IF ( DO_USER_STREAMS ) THEN
c          IF ( DO_DBCORRECTION ) THEN
c            CALL VLIDORT_DBCORRECTION ( SS_FLUX_MULTIPLIER )
c            IF ( DO_PROFILE_LINEARIZATION ) THEN
c              DO LAYER_TO_VARY = 1, NLAYERS
c                IF ( LAYER_VARY_FLAG(LAYER_TO_VARY) ) THEN
c                  N_PARAMETERS = LAYER_VARY_NUMBER(LAYER_TO_VARY)
c                  CALL VLIDORT_LAP_DBCORRECTION
c     &            ( SS_FLUX_MULTIPLIER, LAYER_TO_VARY, N_PARAMETERS)
c                ENDIF
c              ENDDO
c            ENDIF
c            IF ( DO_COLUMN_LINEARIZATION ) THEN
c              CALL VLIDORT_LAC_DBCORRECTION
c     &           ( SS_FLUX_MULTIPLIER, N_TOTALCOLUMN_WFS )
c            ENDIF
c            IF ( DO_SURFACE_LINEARIZATION ) THEN
c              CALL VLIDORT_LS_DBCORRECTION
c            ENDIF
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

C  Linearization of BRDF Fourier components

c      IF ( .NOT. DO_LAMBERTIAN_SURFACE ) THEN
c        CALL VLIDORT_L_BRDF_FRTERMS
c     I    ( DO_INCLUDE_SURFACE,
c     I      DO_INCLUDE_SURFEMISS,
c     I      DELTA_FACTOR )
c      ENDIF

C  Do the DB-corrected Reflected Direct beam attenuation (Lambertian)
C         --- when DO_SSFULL (only single scatter result)
C         --- when Full Radinace + DB correction is flagged
C  Otherwise, direct beam term is set by internal post-processing

      IF ( DO_SOLAR_SOURCES ) THEN
       IF ( DO_SSFULL.OR.(.NOT.DO_SSFULL.AND.DO_DBCORRECTION) ) THEN
        CALL VLIDORT_LAMBERTIAN_DBCORR (SS_FLUX_MULTIPLIER)
        IF ( DO_PROFILE_LINEARIZATION ) THEN
          CALL VLIDORT_LAP_LAMBERTIAN_DBCORR (SS_FLUX_MULTIPLIER)
        ENDIF
        IF ( DO_COLUMN_LINEARIZATION ) THEN
          CALL VLIDORT_LAC_LAMBERTIAN_DBCORR (SS_FLUX_MULTIPLIER)
        ENDIF
        IF ( DO_SURFACE_LINEARIZATION ) THEN
          CALL VLIDORT_LS_LAMBERTIAN_DBCORR
        ENDIF
       ELSE
        CALL LAMBERTIAN_SURFACE_DIRECTBEAM
     I    ( FOURIER_COMPONENT, DELTA_FACTOR )
       ENDIF
      ENDIF

C  Full single scatter calculation, return after completion

      IF ( DO_SSFULL ) RETURN

C  ###############################################################
C     END OF SSCORRECTION STUFF
C  ###############################################################

C  Get Pi Matrices for this Fourier component
C    Including all beam angles !!

      IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
        CALL VLIDORT_PIMATRIX_SETUP ( FOURIER_COMPONENT )
      ENDIF

C  #########################################################
C  RT differential equation Eigensolutions + linearizations
C  #########################################################

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
          TRACE= 'Call #1 in VLIDORT_L_FOURIER_LAMBERTIAN'
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

c         if ( do_write_screen) write(*,*)'l_homsolution',layer
          CALL VLIDORT_L_QHOM_SOLUTION
     I      ( LAYER, FOURIER_COMPONENT,
     I        DO_RTSOL_VARY,
     I        NPARAMS_VARY,
     O        STATUS_SUB )

C  .. error tracing

          IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
            MAIL = 'Error return from module VLIDORT_L_QHOM_SOLUTION'
            TRACE= 'Call #2 in VLIDORT_L_FOURIER_LAMBERTIAN'
            STATUS = VLIDORT_SERIOUS
            CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
            RETURN
          ENDIF

C  Get Linearizations of ("user") solutions for this layer

          CALL VLIDORT_L_UHOM_SOLUTION
     I      ( LAYER, FOURIER_COMPONENT,
     I        DO_RTSOL_VARY,
     I        NPARAMS_VARY )

C  End linearization control

        ENDIF

C  end layer loop

      ENDDO

C  Prepare homogeneous solution multipliers

      CALL HMULT_MASTER
      IF ( DO_ATMOS_LINEARIZATION ) THEN
        CALL L_HMULT_MASTER
      ENDIF

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
          MAIL = 'Error return from module BVPMATRIX_SETUP_MASTER'
          TRACE= 'Call #3 in VLIDORT_L_FOURIER_LAMBERTIAN'
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
          TRACE= 'Call #4 in VLIDORT_L_FOURIER_LAMBERTIAN'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ENDIF

      ENDIF

C  ################
C  Thermal Solution (2 steps)
C  ################

C  Continuation point for avoiding the scattering calculations

 8899 continue

C  THERMAL SOLUTIONS
C  =================

C  Separate calls if linearization is required

C  1. Find the Particular solution (also for transmittance only)
C  2. Compute thermal layer source terms. (Upwelling and Downwelling)
C    These will be scaled up by factor 4.pi if solar beams as well

      IF ( DO_INCLUDE_THERMEMISS ) THEN 

C  With linearization

       IF ( DO_ATMOS_LINEARIZATION ) THEN
         CALL THERMAL_CLSOLUTION_PLUS ( STATUS_SUB, MAIL )
         IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
           TRACE  = 'Call #4A in VLIDORT_L_FOURIER_LAMBERTIAN'
           STATUS = VLIDORT_SERIOUS
           CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
           RETURN
         ENDIF
         IF ( DO_UPWELLING ) CALL THERMAL_STERMS_UP_PLUS
         IF ( DO_DNWELLING ) CALL THERMAL_STERMS_DN_PLUS

C  No linearization

       ELSE
         CALL THERMAL_CLSOLUTION ( STATUS_SUB, MAIL )
         IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
           TRACE  = 'Call #4B in VLIDORT_L_FOURIER_LAMBERTIAN'
           STATUS = VLIDORT_SERIOUS
           CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
           RETURN
         ENDIF
         IF ( DO_UPWELLING ) CALL THERMAL_STERMS_UP
         IF ( DO_DNWELLING ) CALL THERMAL_STERMS_DN
       ENDIF

!  End include thermal emission

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
        TRACE= 'Call #5A in VLIDORT_L_FOURIER_LAMBERTIAN'
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

C  LTE linearization (new section, 14 September 2009)
C  --------------------------------------------------

      if ( DO_THERMAL_TRANSONLY .and. do_LTE_linearization ) then
         CALL THERMAL_LTE_LINEARIZATION
     &   ( DO_INCLUDE_SURFACE, SURFACE_FACTOR, FLUX_MULTIPLIER )
      endif

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

C  Get the BVP solution (regular case) for this beam component

              CALL L_BVP_SOLUTION_MASTER
     I       ( DO_INCLUDE_SURFACE,
     I         DO_INCLUDE_DIRECTBEAM,
     I         DO_INCLUDE_THERMEMISS,
     I         LAYER_TO_VARY, N_LAYER_WFS,
     I         FOURIER_COMPONENT, IPARTIC, SURFACE_FACTOR,
     I         L_BVP_LAMBERTIAN_SURFACE,
     O         STATUS_SUB )

              IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
               MAIL = 'Error in L_BVP_SOLUTION_MASTER, Prof/Lin'
               TRACE= 'Call #5B in VLIDORT_L_FOURIER_LAMBERTIAN'
               STATUS = VLIDORT_SERIOUS
               CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
               RETURN
              ENDIF

            ENDIF

C  Continuation point for avoiding scattering solution

 788        CONTINUE

C  Post-processing for the weighting functions, Sequence:
C    Upwelling   Atmospheric weighting functions
C    Downwelling Atmospheric weighting functions
C    Mean-value  Atmospheric weighting functions

            IF ( DO_UPWELLING ) THEN
              CALL VLIDORT_UPUSER_ATMOSWF
     I         ( DO_INCLUDE_SURFACE,
     I           DO_INCLUDE_THERMEMISS,
     I           DO_INCLUDE_MVOUTPUT,
     I           DO_INCLUDE_DIRECTBEAM,
     I           SURFACE_FACTOR,      FLUX_MULTIPLIER,
     I           FOURIER_COMPONENT,   IPARTIC,
     I           LAYER_TO_VARY,       N_LAYER_WFS,
     I           L_BOA_LAMBERTIAN_SOURCE )
            ENDIF

            IF ( DO_DNWELLING ) THEN
              CALL VLIDORT_DNUSER_ATMOSWF
     I         ( DO_INCLUDE_THERMEMISS,
     I           FLUX_MULTIPLIER,
     I           FOURIER_COMPONENT,   IPARTIC,
     I           LAYER_TO_VARY,       N_LAYER_WFS )
            ENDIF

            IF ( DO_INCLUDE_MVOUTPUT .OR. DO_QUAD_OUTPUT ) THEN
              CALL VLIDORT_L_INTEGRATED_OUTPUT
     I            ( DO_INCLUDE_MVOUTPUT,
     I              DO_INCLUDE_DIRECTBEAM,
     I              DO_INCLUDE_THERMEMISS,
     I              FLUX_MULTIPLIER,
     I              IPARTIC, LAYER_TO_VARY, N_LAYER_WFS )
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

C  Get the linearized BVP solution for this beam component

        CALL L_BVP_SOLUTION_MASTER
     I       ( DO_INCLUDE_SURFACE, 
     I         DO_INCLUDE_DIRECTBEAM,
     I         DO_INCLUDE_THERMEMISS,
     I         VARIATION_INDEX, N_TOTALCOLUMN_WFS,
     I         FOURIER_COMPONENT, IBEAM, SURFACE_FACTOR,
     I         L_BVP_LAMBERTIAN_SURFACE,
     O         STATUS_SUB )

        IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
          MAIL = 'Error in module L_BVP_SOLUTION_MASTER, Col/lin'
          TRACE= 'Call #5C in VLIDORT_L_FOURIER_LAMBERTIAN'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ENDIF

C  Continuation point for avoiding scattering solution

 789    CONTINUE

C  Post-processing for the weighting functions, Sequence:
C    Upwelling   Atmospheric weighting functions
C    Downwelling Atmospheric weighting functions
C    Mean-value  Atmospheric weighting functions

        IF ( DO_UPWELLING ) THEN
          CALL VLIDORT_UPUSER_ATMOSWF
     I         ( DO_INCLUDE_SURFACE,
     I           DO_INCLUDE_THERMEMISS,
     I           DO_INCLUDE_MVOUTPUT,
     I           DO_INCLUDE_DIRECTBEAM,
     I           SURFACE_FACTOR,      FLUX_MULTIPLIER,
     I           FOURIER_COMPONENT,   IPARTIC,
     I           VARIATION_INDEX,     N_TOTALCOLUMN_WFS,
     I           L_BOA_LAMBERTIAN_SOURCE )
        ENDIF

        IF ( DO_DNWELLING ) THEN
          CALL VLIDORT_DNUSER_ATMOSWF
     I         ( DO_INCLUDE_THERMEMISS,
     I           FLUX_MULTIPLIER,
     I           FOURIER_COMPONENT,   IPARTIC,
     I           VARIATION_INDEX,     N_TOTALCOLUMN_WFS )
        ENDIF

        IF ( DO_INCLUDE_MVOUTPUT .OR. DO_QUAD_OUTPUT ) THEN
          CALL VLIDORT_L_INTEGRATED_OUTPUT
     I            ( DO_INCLUDE_MVOUTPUT,
     I              DO_INCLUDE_DIRECTBEAM,
     I              DO_INCLUDE_THERMEMISS,
     I              FLUX_MULTIPLIER,
     I              IPARTIC, VARIATION_INDEX, N_TOTALCOLUMN_WFS )
        ENDIF

C  End atmospheric column weighting functions

      ENDIF

C  Thermal-only: Surface Reflectance weighting functions
C  -----------------------------------------------------

      IF ( DO_SURFACE_LINEARIZATION ) THEN

C  initialise weighting function count

        WF_INDEX = 1
        BRDF_INDEX = 1

C  Only do this if there is surface reflectance !

        IF ( DO_INCLUDE_SURFACE ) THEN

C  Lambertian albedo weighting function

          CALL VLIDORT_ALBEDO_WFS
     I           (  DO_REFLECTED_DIRECTBEAM,
     I              DO_INCLUDE_SURFEMISS,
     I              DO_INCLUDE_MVOUTPUT,
     I              BVP_REGULAR_FLAG(FOURIER_COMPONENT),
     I              FOURIER_COMPONENT, IPARTIC,
     I              SURFACE_FACTOR,
     I              FLUX_MULTIPLIER,
     I              BRDF_INDEX, WF_INDEX,
     I              LS_BVP_LAMBERTIAN_SURFACE,
     I              LS_BOA_LAMBERTIAN_SOURCE,
     O              STATUS_SUB )

          IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
            MAIL = 'Error in module VLIDORT_ALBEDO_WFS'
            TRACE= 'Call #5D in VLIDORT_L_FOURIER_LAMBERTIAN'
            STATUS = VLIDORT_SERIOUS
            CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
            RETURN
          ENDIF

        ENDIF

C  end of albedo weighting functions

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

              CALL  VLIDORT_QBEAM_SOLUTION
     I         ( LAYER, FOURIER_COMPONENT, IBEAM, STATUS_SUB )

              IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
                WRITE(C2,'(I2)')LAYER
                MAIL = 'Error from VLIDORT_QBEAM_SOLUTION, layer'//C2
                TRACE= 'Call #6A in VLIDORT_L_FOURIER_LAMBERTIAN'
                STATUS = VLIDORT_SERIOUS
                CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
                RETURN
              ENDIF

C  For the linearizations of the classical solution

              IF ( DO_ATMOS_LINEARIZATION ) THEN

                CALL VLIDORT_L_QBEAM_SOLUTION
     I            ( LAYER, FOURIER_COMPONENT, IBEAM,
     I              DO_RTSOL_VARY, 
     I              NPARAMS_VARY,
     O              STATUS_SUB )

                IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
                  MAIL =
     &             'Error return from module VLIDORT_L_QBEAM_SOLUTION'
                  TRACE= 'Call #6B in VLIDORT_L_FOURIER_LAMBERTIAN'
                  STATUS = VLIDORT_SERIOUS
                  CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
                  RETURN
                ENDIF

              ENDIF

C  B. the Green's function solution
C  ++++++++++++++++++++++++++++++++

            ELSE

C  $$$$$$$$$$$$$$$$$$$$$$    P L A C E H O L D E R    $$$$$$$$$$$$$$$$$$$$
C            CALL V_GREENFUNC_SOLUTION ( FOURIER)
C  $$$$$$$$$$$$$$$$$$$$$$    P L A C E H O L D E R    $$$$$$$$$$$$$$$$$$$$

            ENDIF

C  C. user solutions
C  +++++++++++++++++

            IF  ( STERM_LAYERMASK_UP(LAYER) .OR.
     &            STERM_LAYERMASK_DN(LAYER) ) THEN
              IF ( DO_USER_STREAMS ) THEN

C  (classical)
                IF ( DO_CLASSICAL_SOLUTION ) THEN

                  CALL VLIDORT_UBEAM_SOLUTION
     I             ( LAYER, FOURIER_COMPONENT, IBEAM )
                  IF ( DO_ATMOS_LINEARIZATION ) THEN
                    CALL VLIDORT_L_UBEAM_SOLUTION
     I               ( LAYER, FOURIER_COMPONENT, IBEAM,
     I                 DO_RTSOL_VARY, 
     I                 NPARAMS_VARY )
                  ENDIF

C  (Green's function)

                ELSE
C  $$$$$$$$$$$$$$$$$$$$$$    P L A C E H O L D E R    $$$$$$$$$$$$$$$$$$$$
C            CALL V_GREENFUNC_SOLUTION ( FOURIER)
C  $$$$$$$$$$$$$$$$$$$$$$    P L A C E H O L D E R    $$$$$$$$$$$$$$$$$$$$
                ENDIF

              END IF
            END IF

C  end layer loop

          END DO

C  Add thermal solutions if flagged
C    4.PI modulus on the thermal contribution.

          IF ( DO_INCLUDE_THERMEMISS ) THEN
            O1 = 1
            DO N = 1, NLAYERS
              DO I = 1, NSTREAMS_2
                WUPPER(I,O1,N) = WUPPER(I,O1,N) + PI4 * T_WUPPER(I,N)
                WLOWER(I,O1,N) = WLOWER(I,O1,N) + PI4 * T_WLOWER(I,N)
              ENDDO
            ENDDO
          ENDIF

C  Adding the linearized thermal solutions is done later.....

C  Step 2. Solve boundary value problem and get Stokes vector
C  ----------------------------------------------------------

C  2A. Boundary Value problem
C  ==========================

C  standard case using compressed-band matrices, etc..

c       if ( do_write_screen) write(*,*)'bvp solution',ibeam
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
              TRACE= 'Call #7A in VLIDORT_L_FOURIER_LAMBERTIAN'
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
              TRACE= 'Call #7B in VLIDORT_L_FOURIER_LAMBERTIAN'
              STATUS = VLIDORT_SERIOUS
              CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
              RETURN
            ENDIF

          ENDIF

C  2B. Stokes vector Post Processing
C  =================================

C  upwelling

          IF ( DO_UPWELLING ) THEN

            DO_INCLUDE_DIRECTBEAM = 
     &      (DO_REFLECTED_DIRECTBEAM(IBEAM).AND..NOT.DO_DBCORRECTION)

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
     &        FLUX_MULTIPLIER, IBEAM )
          ENDIF

C  Finished this Beam solution, if only Stokes vector is required

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

C    (a) Regular case (full matrix linearization)

             IF ( BVP_REGULAR_FLAG(FOURIER_COMPONENT) ) THEN

              CALL L_BVP_SOLUTION_MASTER
     I       ( DO_INCLUDE_SURFACE,
     I         DO_REFLECTED_DIRECTBEAM(IBEAM),
     I         DO_INCLUDE_THERMEMISS,
     I         LAYER_TO_VARY, N_LAYER_WFS,
     I         FOURIER_COMPONENT, IBEAM, SURFACE_FACTOR,
     I         L_BVP_LAMBERTIAN_SURFACE,
     O         STATUS_SUB )

              IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
               MAIL = 'Error in module L_BVP_SOLUTION_MASTER'
               TRACE= 'Call #8A in VLIDORT_L_FOURIER_LAMBERTIAN'
               STATUS = VLIDORT_SERIOUS
               CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
               RETURN
              ENDIF

C    (b) telescoped case 

            ELSE

             CALL L_BVPTEL_SOLUTION_MASTER
     I       ( DO_INCLUDE_SURFACE,
     I         DO_REFLECTED_DIRECTBEAM(IBEAM),
     I         LAYER_TO_VARY, N_LAYER_WFS,
     I         FOURIER_COMPONENT, IBEAM, SURFACE_FACTOR,
     I         L_BVP_LAMBERTIAN_SURFACE,
     O         STATUS_SUB )

              IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
               MAIL = 'Error in module L_BVPTEL_SOLUTION_MASTER'
               TRACE= 'Call #8B in VLIDORT_L_FOURIER_LAMBERTIAN'
               STATUS = VLIDORT_SERIOUS
               CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
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
               CALL VLIDORT_UPUSER_ATMOSWF
     I         ( DO_INCLUDE_SURFACE,
     I           DO_INCLUDE_THERMEMISS,
     I           DO_INCLUDE_MVOUTPUT,
     I           DO_REFLECTED_DIRECTBEAM(IBEAM),
     I           SURFACE_FACTOR,      FLUX_MULTIPLIER,
     I           FOURIER_COMPONENT,   IBEAM,
     I           LAYER_TO_VARY,       N_LAYER_WFS,
     I           L_BOA_LAMBERTIAN_SOURCE )
             ENDIF

             IF ( DO_DNWELLING ) THEN
C           if ( do_write_screen) write(*,*)'wf dn',ibeam,layer_to_vary
              CALL VLIDORT_DNUSER_ATMOSWF
     I         ( DO_INCLUDE_THERMEMISS,
     I           FLUX_MULTIPLIER,
     I           FOURIER_COMPONENT,   IBEAM,
     I           LAYER_TO_VARY,       N_LAYER_WFS )
             ENDIF

             IF ( DO_INCLUDE_MVOUTPUT .OR. DO_QUAD_OUTPUT ) THEN
              CALL VLIDORT_L_INTEGRATED_OUTPUT
     I            ( DO_INCLUDE_MVOUTPUT,
     I              DO_INCLUDE_DIRECTBEAM,
     I              DO_INCLUDE_THERMEMISS,
     I              FLUX_MULTIPLIER,
     I              IBEAM, LAYER_TO_VARY, N_LAYER_WFS )
             ENDIF

C  Finish loop over layers with variation

            ENDIF
           ENDDO

C  End profile atmospheric weighting functions

          ENDIF

C  Step 4. Atmospheric Bulk weighting functions
C  --------------------------------------------

          IF ( DO_COLUMN_LINEARIZATION ) THEN

C   variation index = 0

            VARIATION_INDEX = 0

C  4A. Solve the linearized BVP
C  ============================

C  (a) Regular BVP linearization

            IF ( BVP_REGULAR_FLAG(FOURIER_COMPONENT) ) THEN

C  Get the linearized BVP solution for this beam component

              CALL L_BVP_SOLUTION_MASTER
     I       ( DO_INCLUDE_SURFACE,
     I         DO_REFLECTED_DIRECTBEAM(IBEAM),
     I         DO_INCLUDE_THERMEMISS,
     I         VARIATION_INDEX, N_TOTALCOLUMN_WFS,
     I         FOURIER_COMPONENT, IBEAM, SURFACE_FACTOR,
     I         L_BVP_LAMBERTIAN_SURFACE,
     O         STATUS_SUB )

              IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
               MAIL = 'Error in module L_BVP_SOLUTION_MASTER'
               TRACE= 'Call #9A in VLIDORT_L_FOURIER_LAMBERTIAN'
               STATUS = VLIDORT_SERIOUS
               CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
               RETURN
              ENDIF

C  (b) telescoped BVP linearization

            ELSE

              CALL L_BVPTEL_SOLUTION_MASTER
     I       ( DO_INCLUDE_SURFACE,
     I         DO_REFLECTED_DIRECTBEAM(IBEAM),
     I         VARIATION_INDEX, N_TOTALCOLUMN_WFS,
     I         FOURIER_COMPONENT, IBEAM, SURFACE_FACTOR,
     I         L_BVP_LAMBERTIAN_SURFACE,
     O         STATUS_SUB )

              IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
               MAIL = 'Error in module L_BVPTEL_SOLUTION_MASTER'
               TRACE= 'Call #9B in VLIDORT_L_FOURIER_LAMBERTIAN'
               STATUS = VLIDORT_SERIOUS
               CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
               RETURN
              ENDIF

            ENDIF

C  4B. Post-processing for the weighting functions
C  ===============================================

C  Sequence:
C    Upwelling   Atmospheric weighting functions
C    Downwelling Atmospheric weighting functions
C    Mean-value  Atmospheric weighting functions

            IF ( DO_UPWELLING ) THEN
c            if ( do_write_screen) write(*,*)'wf up',ibeam,layer_to_vary
               CALL VLIDORT_UPUSER_ATMOSWF
     I         ( DO_INCLUDE_SURFACE,
     I           DO_INCLUDE_THERMEMISS,
     I           DO_INCLUDE_MVOUTPUT,
     I           DO_REFLECTED_DIRECTBEAM(IBEAM),
     I           SURFACE_FACTOR,      FLUX_MULTIPLIER,
     I           FOURIER_COMPONENT,   IBEAM,
     I           VARIATION_INDEX,     N_TOTALCOLUMN_WFS,
     I           L_BOA_LAMBERTIAN_SOURCE )
            ENDIF

            IF ( DO_DNWELLING ) THEN
c           if ( do_write_screen) write(*,*)'wf dn',ibeam,layer_to_vary
              CALL VLIDORT_DNUSER_ATMOSWF
     I         ( DO_INCLUDE_THERMEMISS,
     I           FLUX_MULTIPLIER,
     I           FOURIER_COMPONENT,   IBEAM,
     I           VARIATION_INDEX,     N_TOTALCOLUMN_WFS )
            ENDIF

            IF ( DO_INCLUDE_MVOUTPUT .OR. DO_QUAD_OUTPUT ) THEN
              CALL VLIDORT_L_INTEGRATED_OUTPUT
     I            ( DO_INCLUDE_MVOUTPUT,
     I              DO_REFLECTED_DIRECTBEAM(IBEAM),
     I              DO_INCLUDE_THERMEMISS,
     I              FLUX_MULTIPLIER,
     I              IBEAM, VARIATION_INDEX, N_TOTALCOLUMN_WFS )
            ENDIF

C  End atmospheric column weighting functions

          ENDIF

C  Step 5. Surface Reflectance weighting functions
C  -----------------------------------------------

          IF ( DO_SURFACE_LINEARIZATION ) THEN

C  initialise weighting function count

            WF_INDEX   = 1
            BRDF_INDEX = 1

C  Only do this if there is surface reflectance !

            IF ( DO_INCLUDE_SURFACE ) THEN

C  Lambertian albedo weighting function

              CALL VLIDORT_ALBEDO_WFS
     I           (  DO_REFLECTED_DIRECTBEAM(IBEAM),
     I              DO_INCLUDE_SURFEMISS,
     I              DO_INCLUDE_MVOUTPUT,
     I              BVP_REGULAR_FLAG(FOURIER_COMPONENT),
     I              FOURIER_COMPONENT, IBEAM,
     I              SURFACE_FACTOR,
     I              FLUX_MULTIPLIER,
     I              BRDF_INDEX, WF_INDEX,
     I              LS_BVP_LAMBERTIAN_SURFACE,
     I              LS_BOA_LAMBERTIAN_SOURCE,
     O              STATUS_SUB )

              IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
                MAIL = 'Error in module VLIDORT_ALBEDO_WFS'
                TRACE= 'Call #10 in VLIDORT_L_FOURIER_LAMBERTIAN'
                STATUS = VLIDORT_SERIOUS
                CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
                RETURN
              ENDIF

            ENDIF

C  end of albedo weighting functions

          ENDIF

C  Continuation point for avoiding weighting functions

 4000     CONTINUE

C  End loop over beam solutions

        END IF
      END DO

C  debug write

      IF ( DO_DEBUG_WRITE ) THEN
        write(54,'(a,I3,a,100I3)')
     &     'Fourier ',FOURIER_COMPONENT,
     &   ' : # complex evalues by layer: ',
     &       (K_COMPLEX(LAYER),LAYER=1,NLAYERS)
      ENDIF

C  ######
C  finish
C  ######

      RETURN
      END
