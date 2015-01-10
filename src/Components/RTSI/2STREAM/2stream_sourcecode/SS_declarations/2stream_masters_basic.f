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
C #            TWOSTREAM_MASTER (top-level master)              #
C #            TWOSTREAM_FOURIER_MASTER                         #
C #            TWOSTREAM_CONVERGE (master)                      #
C #                                                             #
C ###############################################################

      SUBROUTINE TWOSTREAM_MASTER
     I   ( DO_SSCORR_OUTGOING, DO_SSCORR_NADIR, DO_SSFULL,
     I     DO_UPWELLING, DO_DNWELLING, 
     I     DO_PLANE_PARALLEL, DO_DB_CORRECTION,
     I     NLAYERS, NTOTAL, NFINELAYERS, NTHREADS, NMOMENTS_INPUT,
     I     NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, N_GEOMETRIES,
     I     STREAM_VALUE, FLUX_FACTOR, GEOMETRY_SPECHEIGHT,
     I     BEAM_SZAS, USER_ANGLES, USER_RELAZMS,
     I     SURFTYPE, EARTH_RADIUS, HEIGHT_GRID, THREAD,
     I     LAMBERTIAN_ALBEDO, DELTAU_VERT, OMEGA_TOTAL, ASYMM_TOTAL,
     O     INTENSITY_TOA, INTENSITY_BOA,
     O     STATUS_INPUTCHECK, STATUS_CALCULATION )

C  subroutine input arguments
C  --------------------------

C  Flags

      LOGICAL          DO_SSCORR_OUTGOING, DO_SSCORR_NADIR, DO_SSFULL
      LOGICAL          DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL
      LOGICAL          DO_DB_CORRECTION

C  Numbers

      INTEGER          NLAYERS, NTOTAL, NFINELAYERS
      INTEGER          NMOMENTS_INPUT, NTHREADS
      INTEGER          NBEAMS, N_USER_STREAMS, N_USER_RELAZMS  
      INTEGER          N_GEOMETRIES

C  Geometry

      DOUBLE PRECISION BEAM_SZAS    ( NBEAMS )
      DOUBLE PRECISION USER_ANGLES  ( N_USER_STREAMS )
      DOUBLE PRECISION USER_RELAZMS ( N_USER_RELAZMS )

C  height and earth radius

      DOUBLE PRECISION EARTH_RADIUS
      DOUBLE PRECISION HEIGHT_GRID ( 0:NLAYERS )

C  Stream value

      DOUBLE PRECISION STREAM_VALUE

C  Flux factor

      DOUBLE PRECISION FLUX_FACTOR

C  Geometry specification height

      DOUBLE PRECISION GEOMETRY_SPECHEIGHT

C  surface type

      INTEGER          SURFTYPE

C  Optical properties

      DOUBLE PRECISION DELTAU_VERT(NLAYERS, NTHREADS)
      DOUBLE PRECISION OMEGA_TOTAL(NLAYERS, NTHREADS)
      DOUBLE PRECISION ASYMM_TOTAL(NLAYERS, NTHREADS)
      DOUBLE PRECISION LAMBERTIAN_ALBEDO(NTHREADS)

C  Input thread

      INTEGER          THREAD

C  output
C  ------

C  Results

      DOUBLE PRECISION INTENSITY_TOA(N_GEOMETRIES,NTHREADS)
      DOUBLE PRECISION INTENSITY_BOA(N_GEOMETRIES,NTHREADS)

C  output status

      INTEGER          STATUS_INPUTCHECK
      INTEGER          STATUS_CALCULATION

C  Local definitions
C  -----------------

C  MS-only flag

      LOGICAL          DO_MSMODE_TWOSTREAM

C  Cosines

      DOUBLE PRECISION X0 ( NBEAMS )
      DOUBLE PRECISION USER_STREAMS ( N_USER_STREAMS )

C  Chapman factors (from pseudo-spherical geometry)

      DOUBLE PRECISION CHAPMAN_FACTORS
     &              ( NLAYERS, NLAYERS, NBEAMS )

C     Last layer to include Particular integral solution
C     Average-secant and initial tramsittance factors for solar beams.
C     Solar beam attenuation

      INTEGER          LAYER_PIS_CUTOFF(NBEAMS)
      DOUBLE PRECISION
     &     INITIAL_TRANS  ( NLAYERS, NBEAMS ),
     &     AVERAGE_SECANT ( NLAYERS, NBEAMS ),
     &     LOCAL_SZA      ( NLAYERS, NBEAMS )
      DOUBLE PRECISION SOLAR_BEAM_OPDEP ( NBEAMS )

C  Derived optical thickness inputs

      DOUBLE PRECISION DELTAU_SLANT ( NLAYERS, NLAYERS, NBEAMS )
      DOUBLE PRECISION TAUSLANT    ( 0:NLAYERS, NBEAMS )

C  reflectance flags

      LOGICAL          DO_REFLECTED_DIRECTBEAM ( NBEAMS )

C  Transmittance Setups
C  --------------------

C  Transmittance factors for average secant stream
C    Computed in the initial setup stage for Fourier m = 0

      DOUBLE PRECISION
     &     T_DELT_MUBAR ( NLAYERS,      NBEAMS )

C  Transmittance factors for user-defined stream angles
C    Computed in the initial setup stage for Fourier m = 0

      DOUBLE PRECISION
     &     T_DELT_USERM ( NLAYERS,      N_USER_STREAMS )

C  Multiplier arrays
C  -----------------

C  forcing term multipliers (saved for whole atmosphere)

      DOUBLE PRECISION EMULT_UP
     &       (N_USER_STREAMS,NLAYERS,NBEAMS)

      DOUBLE PRECISION EMULT_DN
     &       (N_USER_STREAMS,NLAYERS,NBEAMS)

C  Band compression index

      INTEGER          BMAT_ROWMASK(NTOTAL,NTOTAL)

C  User-defined solutions

      DOUBLE PRECISION INTENSITY_F_UP
     D   (N_USER_STREAMS,NBEAMS)
      DOUBLE PRECISION INTENSITY_F_DN
     D   (N_USER_STREAMS,NBEAMS)

C  Single scatter solutions

      DOUBLE PRECISION INTENSITY_SS_UP(N_GEOMETRIES)
      DOUBLE PRECISION INTENSITY_SS_DN(N_GEOMETRIES)

C  Other arrays
C  ------------

C  Local error handling

      LOGICAL          INIT
      CHARACTER*70     MAIL, TRACE
      CHARACTER*3      CF, WTHREAD

C  help variables

      LOGICAL          ADJUST_SURFACE

      INTEGER          FOURIER
      INTEGER          N_FOURIERS

      INTEGER          UA, UM, IB, N_VIEWING, IBEAM, I
      INTEGER          STATUS_SUB

      DOUBLE PRECISION SS_FLUX_MULTIPLIER, modified_eradius
      DOUBLE PRECISION AZM_ARGUMENT, DFC

C  Local azimuth factors

      INTEGER          IBOFF ( NBEAMS )
      INTEGER          UMOFF ( NBEAMS, N_USER_STREAMS )
      DOUBLE PRECISION AZMFAC
     &     (N_USER_STREAMS,NBEAMS,N_USER_RELAZMS)


      DOUBLE PRECISION DEG_TO_RAD, PI4

C  initialize output status
C  ------------------------

      STATUS_CALCULATION = 0
      STATUS_INPUTCHECK  = 0

C  initialise error flag

      INIT = .TRUE.

C  constants

      DEG_TO_RAD = DACOS(-1.0d0)/180.0d0
      PI4 = DEG_TO_RAD * 720.0d0

C  thread number

      WTHREAD = '000'
      IF (THREAD.LT.10)WRITE(WTHREAD(3:3),'(I1)')THREAD
      IF (THREAD.GT.99)WRITE(WTHREAD(1:3),'(I3)')THREAD
      IF (THREAD.GE.10.and.THREAD.LE.99)WRITE(WTHREAD(2:3),'(I2)')THREAD

C  Single scatter correction: flux multiplier
C    Now always F / 4pi

      SS_FLUX_MULTIPLIER = FLUX_FACTOR / PI4

C  Check input Basic. This could be put outside the thread loop.

      CALL  TWOSTREAM_CHECK_INPUTS_BASIC
     I   ( DO_SSCORR_OUTGOING, DO_SSCORR_NADIR, DO_SSFULL,
     I     DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL,
     I     NLAYERS, NFINELAYERS, NBEAMS, NMOMENTS_INPUT,
     I     N_USER_STREAMS, N_USER_RELAZMS,
     I     BEAM_SZAS, USER_ANGLES, USER_RELAZMS,
     I     EARTH_RADIUS, HEIGHT_GRID, INIT,
     O     STATUS_SUB )

      IF ( STATUS_SUB .EQ. 1 ) THEN
        STATUS_INPUTCHECK = 1
        MAIL = ' TWOSTREAM_CHECK_INPUT_BASIC failed'      
        TRACE = ' Called in TWOSTREAM_MASTER, initial check'
        CALL TWOSTREAM_ERROR_TRACE
     &      ( INIT, MAIL, TRACE, STATUS_INPUTCHECK )
        RETURN
      ENDIF

C  Check input threaded values (IOPs, albedo)

      CALL TWOSTREAM_CHECK_INPUTS_THREAD
     I     ( INIT, NLAYERS, NTHREADS, THREAD, 
     I       SURFTYPE, LAMBERTIAN_ALBEDO,
     I       DELTAU_VERT, OMEGA_TOTAL, ASYMM_TOTAL,
     O       STATUS_SUB )
      IF ( STATUS_SUB .EQ. 1 ) THEN
        STATUS_INPUTCHECK = 1
        MAIL = ' TWOSTREAM_CHECK_INPUT_THREAD failed'      
        TRACE = ' Called in TWOSTREAM_MASTER , THREAD # '//wthread
        CALL TWOSTREAM_ERROR_TRACE
     &      ( INIT, MAIL, TRACE, STATUS_INPUTCHECK )
        RETURN
      ENDIF

C  save some offsets for indexing geometries

      N_VIEWING    = N_USER_STREAMS * N_USER_RELAZMS
c      N_GEOMETRIES = NBEAMS * N_VIEWING

      DO IBEAM = 1, NBEAMS
        IBOFF(IBEAM) = N_VIEWING * ( IBEAM - 1 )
        DO UM = 1, N_USER_STREAMS
          UMOFF(IBEAM,UM) = IBOFF(IBEAM) +  N_USER_RELAZMS * (UM - 1)
        END DO
      END DO

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
c      CALL multi_outgoing_adjustgeom
c     i   ( N_USER_STREAMS, NBEAMS, N_USER_RELAZMS,
c     i     N_USER_STREAMS,   NBEAMS,   N_USER_RELAZMS,
c     i     height_grid(nlayers), modified_eradius, adjust_surface,
c     i     user_angles,  beam_szas, user_relazms,
c     o     user_angles_adjust, beam_szas_adjust, user_relazms_adjust,
c     o     fail, mail )
c      if ( fail ) return
     
C  Chapman function calculation
C  ----------------------------

C  Chapman function calculation
C  ----------------------------

C  start beam loop

      DO IB = 1, NBEAMS

C  Get the factors

        CALL TWOSTREAM_BEAM_GEOMETRY_PREPARE
     I     ( DO_PLANE_PARALLEL, NBEAMS, NLAYERS, IB,
     I       BEAM_SZAS(IB), EARTH_RADIUS, HEIGHT_GRID,
     O       CHAPMAN_FACTORS, LOCAL_SZA(1,IB) )
C  End Beam loop

      ENDDO

C  Get derived inputs
C  ==================

C  Mode of operation

      IF ( DO_SSCORR_NADIR .OR. DO_SSCORR_OUTGOING ) THEN
        DO_MSMODE_TWOSTREAM = .TRUE.
      ELSE
        DO_MSMODE_TWOSTREAM = .FALSE.
      ENDIF

C  solar zenith angle cosines/sines

      DO IB = 1, NBEAMS
        X0(IB)  = DCOS ( BEAM_SZAS(IB) * DEG_TO_RAD )
      ENDDO 

C  User stream cosines and secants

      DO I = 1, N_USER_STREAMS
        USER_STREAMS(I) = DCOS(DEG_TO_RAD*USER_ANGLES(I))
      ENDDO

C  Initialise Fourier loop
C  =======================

C  set Fourier number (2 for Rayleigh only, otherwise 2*Nstreams))

      IF ( DO_SSFULL ) THEN
        N_FOURIERS = 0
      ELSE
        N_FOURIERS = 1
      ENDIF

C  Fourier loop
C  ============

      DO  FOURIER = 0, N_FOURIERS

C  azimuth cosine factor, using adjust geometries.

        IF ( FOURIER .GT. 0 ) THEN
          DFC = DBLE(FOURIER)
          DO UA = 1, N_USER_RELAZMS
            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
c                AZM_ARGUMENT = USER_RELAZMS_ADJUST(UM,IB,UA) * DFC
                AZM_ARGUMENT = USER_RELAZMS(UA) * DFC
                AZMFAC(UM,IB,UA)   = DCOS(DEG_TO_RAD*AZM_ARGUMENT)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

C  Main call to Lidort Fourier module
C  ----------------------------------

        CALL TWOSTREAM_FOURIER_MASTER
     I   ( DO_SSCORR_OUTGOING, DO_SSCORR_NADIR, DO_SSFULL,
     I     DO_UPWELLING, DO_DNWELLING, DO_MSMODE_TWOSTREAM,
     I     DO_PLANE_PARALLEL, DO_DB_CORRECTION,
     I     NLAYERS, NTOTAL, NFINELAYERS, NTHREADS, NMOMENTS_INPUT,
     I     NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, N_GEOMETRIES,
     I     THREAD, FOURIER, SURFTYPE, INIT,
     I     STREAM_VALUE, FLUX_FACTOR, SS_FLUX_MULTIPLIER, 
     I     CHAPMAN_FACTORS, LAMBERTIAN_ALBEDO,
     I     X0, USER_STREAMS, DELTAU_VERT, OMEGA_TOTAL, ASYMM_TOTAL,
     H     INITIAL_TRANS, AVERAGE_SECANT, LOCAL_SZA, LAYER_PIS_CUTOFF,
     H     DELTAU_SLANT, TAUSLANT, SOLAR_BEAM_OPDEP,
     H     DO_REFLECTED_DIRECTBEAM, BMAT_ROWMASK,
     H     T_DELT_MUBAR, T_DELT_USERM, EMULT_UP, EMULT_DN, 
     O     INTENSITY_F_UP, INTENSITY_SS_UP,
     O     INTENSITY_F_DN, INTENSITY_SS_DN,
     O     STATUS_SUB )

C  error handling

        IF ( STATUS_SUB .NE. 0 ) THEN
          STATUS_CALCULATION = 1
          WRITE(CF,'(I3)')FOURIER
          MAIL = 'Error from TWOSTREAM_FOURIER_MASTER, Fourier = '//CF
          TRACE = ' Called in TWOSTREAM_MASTER , THREAD # '//wthread
          CALL TWOSTREAM_ERROR_TRACE
     &        ( INIT, MAIL, TRACE, STATUS_CALCULATION )
          RETURN
        ENDIF

C  Fourier summation and Convergence examination
C  ---------------------------------------------

C   -- only done for beams which are still not converged
C      This is controlled by flag DO_MULTIBEAM

C   -- new criterion, SS is added for Fourier = 0, as this means that
C      higher-order terms will be relatively smaller, which implies
C      faster convergence in some circumstances (generally not often).

        DO IBEAM = 1, NBEAMS
          CALL TWOSTREAM_CONVERGE
     I    ( DO_UPWELLING, DO_DNWELLING,
     I      DO_SSFULL, DO_SSCORR_OUTGOING, DO_SSCORR_NADIR,
     I      N_GEOMETRIES, NBEAMS, NTHREADS,
     I      N_USER_STREAMS, N_USER_RELAZMS, AZMFAC, UMOFF,
     I      THREAD, IBEAM, FOURIER,
     I      INTENSITY_F_UP,  INTENSITY_F_DN, 
     I      INTENSITY_SS_UP, INTENSITY_SS_DN,
     O      INTENSITY_TOA, INTENSITY_BOA )
        END DO

C  end iteration loop

      ENDDO

C  close Error file if it was used

      IF ( .NOT. INIT ) CLOSE(25)

C  Finish

      RETURN
      END

C

      SUBROUTINE TWOSTREAM_FOURIER_MASTER
     I   ( DO_SSCORR_OUTGOING, DO_SSCORR_NADIR, DO_SSFULL,
     I     DO_UPWELLING, DO_DNWELLING, DO_MSMODE_TWOSTREAM,
     I     DO_PLANE_PARALLEL, DO_DB_CORRECTION,
     I     NLAYERS, NTOTAL, NFINELAYERS, NTHREADS, NMOMENTS_INPUT,
     I     NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, N_GEOMETRIES,
     I     THREAD, FOURIER, SURFTYPE, INIT,
     I     STREAM_VALUE, FLUX_FACTOR, SS_FLUX_MULTIPLIER,
     I     CHAPMAN_FACTORS, LAMBERTIAN_ALBEDO,
     I     X0, USER_STREAMS, DELTAU_VERT, OMEGA_TOTAL, ASYMM_TOTAL,
     H     INITIAL_TRANS, AVERAGE_SECANT, LOCAL_SZA, LAYER_PIS_CUTOFF,
     H     DELTAU_SLANT, TAUSLANT, SOLAR_BEAM_OPDEP,
     H     DO_REFLECTED_DIRECTBEAM, BMAT_ROWMASK,
     H     T_DELT_MUBAR, T_DELT_USERM, EMULT_UP, EMULT_DN, 
     O     INTENSITY_F_UP, INTENSITY_SS_UP,
     O     INTENSITY_F_DN, INTENSITY_SS_DN,
     O     STATUS )

C  input
C  -----

C  Flags

      LOGICAL          DO_SSCORR_OUTGOING, DO_SSCORR_NADIR, DO_SSFULL
      LOGICAL          DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL
      LOGICAL          DO_DB_CORRECTION, DO_MSMODE_TWOSTREAM

C  Numbers

      INTEGER          NLAYERS, NTOTAL, NFINELAYERS
      INTEGER          NMOMENTS_INPUT, NTHREADS
      INTEGER          NBEAMS, N_USER_STREAMS, N_USER_RELAZMS  
      INTEGER          N_GEOMETRIES

C  Geometry

      DOUBLE PRECISION X0           ( NBEAMS )
      DOUBLE PRECISION USER_STREAMS ( N_USER_STREAMS )

C  Optical properties

      DOUBLE PRECISION DELTAU_VERT(NLAYERS, NTHREADS)
      DOUBLE PRECISION OMEGA_TOTAL(NLAYERS, NTHREADS)
      DOUBLE PRECISION ASYMM_TOTAL(NLAYERS, NTHREADS)
      DOUBLE PRECISION LAMBERTIAN_ALBEDO(NTHREADS)

C  Thread number

      INTEGER          THREAD

C  Input Fourier component number

      INTEGER          FOURIER

C  surface type

      INTEGER          SURFTYPE

C  Stream value

      DOUBLE PRECISION STREAM_VALUE

C  Flux factor

      DOUBLE PRECISION FLUX_FACTOR

C  error initialize flag

      LOGICAL          INIT

C  single scatter correction flux multiplier

      DOUBLE PRECISION SS_FLUX_MULTIPLIER

C  Output
C  ------

C  User-defined solutions

      DOUBLE PRECISION INTENSITY_F_UP
     D   (N_USER_STREAMS,NBEAMS)
      DOUBLE PRECISION INTENSITY_F_DN
     D   (N_USER_STREAMS,NBEAMS)

C  Single scatter solutions

      DOUBLE PRECISION INTENSITY_SS_UP (N_GEOMETRIES)
      DOUBLE PRECISION INTENSITY_SS_DN (N_GEOMETRIES)

C  status

      INTEGER          STATUS

C  Arrays required at the Top level
C  ================================

C  Solar beam pseudo-spherical setup
C  ---------------------------------

C  Chapman factors (from pseudo-spherical geometry)

      DOUBLE PRECISION CHAPMAN_FACTORS
     &              ( NLAYERS, NLAYERS, NBEAMS )

C     Last layer to include Particular integral solution
C     Average-secant and initial tramsittance factors for solar beams.
C     Solar beam attenuation

      INTEGER          LAYER_PIS_CUTOFF(NBEAMS)
      DOUBLE PRECISION
     &     INITIAL_TRANS  ( NLAYERS, NBEAMS ),
     &     AVERAGE_SECANT ( NLAYERS, NBEAMS ),
     &     LOCAL_SZA      ( NLAYERS, NBEAMS )
      DOUBLE PRECISION SOLAR_BEAM_OPDEP ( NBEAMS )

C  Derived optical thickness inputs

      DOUBLE PRECISION DELTAU_SLANT ( NLAYERS, NLAYERS, NBEAMS )
      DOUBLE PRECISION TAUSLANT    ( 0:NLAYERS, NBEAMS )

C  reflectance flags

      LOGICAL          DO_REFLECTED_DIRECTBEAM ( NBEAMS )

C  Transmittance Setups
C  --------------------

C  Transmittance factors for average secant stream
C    Computed in the initial setup stage for Fourier m = 0

      DOUBLE PRECISION
     &     T_DELT_MUBAR ( NLAYERS, NBEAMS )

C  Transmittance factors for user-defined stream angles
C    Computed in the initial setup stage for Fourier m = 0

      DOUBLE PRECISION
     &     T_DELT_USERM ( NLAYERS, N_USER_STREAMS )

C  Multiplier arrays
C  -----------------

C  forcing term multipliers (saved for whole atmosphere)

      DOUBLE PRECISION EMULT_UP
     &       (N_USER_STREAMS,NLAYERS,NBEAMS)

      DOUBLE PRECISION EMULT_DN
     &       (N_USER_STREAMS,NLAYERS,NBEAMS)

C  Local Arrays for argument passing
C  =================================

C  Geometry arrays
C  ---------------

C  These just save some Polynomial expansions

      DOUBLE PRECISION ULP ( N_USER_STREAMS )
      DOUBLE PRECISION POX  ( NBEAMS )
      DOUBLE PRECISION PX0X ( NBEAMS )
      DOUBLE PRECISION PX11, PXSQ

C  Solar beam pseudo-spherical setup
C  ---------------------------------

C  Atmospheric attenuation

      DOUBLE PRECISION ATMOS_ATTN ( NBEAMS )

C  Direct beam solutions

      DOUBLE PRECISION
     &        DIRECT_BEAM      ( NBEAMS ),
     &        USER_DIRECT_BEAM ( N_USER_STREAMS, NBEAMS )

C  Transmittance factor

      DOUBLE PRECISION
     &     ITRANS_USERM   ( NLAYERS, N_USER_STREAMS, NBEAMS )

C  Multiplier arrays
C  -----------------

C  coefficient functions for user-defined angles

      DOUBLE PRECISION 
     &      SIGMA_M(NLAYERS,N_USER_STREAMS,NBEAMS),
     &      SIGMA_P(NLAYERS,N_USER_STREAMS,NBEAMS)
      
C  L'Hopital's rule logical variables

      LOGICAL          EMULT_HOPRULE
     &       (NLAYERS,N_USER_STREAMS,NBEAMS)

C  coefficient functions for user-defined angles

      DOUBLE PRECISION 
     &      ZETA_M(N_USER_STREAMS,NLAYERS),
     &      ZETA_P(N_USER_STREAMS,NLAYERS)

C  Integrated homogeneous solution multipliers, whole layer

      DOUBLE PRECISION 
     &      HMULT_1(N_USER_STREAMS,NLAYERS),
     &      HMULT_2(N_USER_STREAMS,NLAYERS)

C  Solutions to the homogeneous RT equations 
C  -----------------------------------------

C  local matrices for eigenvalue computation

      DOUBLE PRECISION SAB(NLAYERS), DAB(NLAYERS)

C  Eigensolutions

      DOUBLE PRECISION EIGENVALUE(NLAYERS)
      DOUBLE PRECISION EIGENTRANS(NLAYERS)

C  Eigenvector solutions

      DOUBLE PRECISION XPOS(2,NLAYERS)
      DOUBLE PRECISION XNEG(2,NLAYERS)

C  Saved help variables

      DOUBLE PRECISION U_HELP_P(0:1)
      DOUBLE PRECISION U_HELP_M(0:1)

C  Eigenvectors defined at user-defined stream angles
C     EP for the positive KEIGEN values, EM for -ve KEIGEN

      DOUBLE PRECISION
     U        U_XPOS(N_USER_STREAMS,NLAYERS),
     U        U_XNEG(N_USER_STREAMS,NLAYERS)

C  Reflected solutions

      DOUBLE PRECISION R2_HOMP
      DOUBLE PRECISION R2_HOMM

C  Boundary Value Problem
C  ----------------------

C  Single Matrix, Band-matrices

      DOUBLE PRECISION SMAT2      (2,2)
      DOUBLE PRECISION BANDMAT2   (7,NTOTAL)
      INTEGER          BMAT_ROWMASK(NTOTAL,NTOTAL)
      INTEGER          IPIVOT     (NTOTAL)

C  particular integrals
C  --------------------

C  Beam solution

      DOUBLE PRECISION WVEC(2,NLAYERS)

C  Solutions at layer boundaries

      DOUBLE PRECISION WUPPER(2,NLAYERS)
      DOUBLE PRECISION WLOWER(2,NLAYERS)

C  Auxiliary vectors

      DOUBLE PRECISION QDIFVEC(NLAYERS)
      DOUBLE PRECISION QSUMVEC(NLAYERS)
      DOUBLE PRECISION QVEC   (NLAYERS)

C  reflected solution

      DOUBLE PRECISION R2_PARTIC

C  Saved help variables

      DOUBLE PRECISION W_HELP(0:1)

C  Particular beam solutions at user-defined stream angles

      DOUBLE PRECISION
     U        U_WPOS1(N_USER_STREAMS,NLAYERS),
     U        U_WNEG1(N_USER_STREAMS,NLAYERS),
     U        U_WPOS2(N_USER_STREAMS,NLAYERS),
     U        U_WNEG2(N_USER_STREAMS,NLAYERS)

C  output
C  ------

C  Column vectors for solving BCs

      DOUBLE PRECISION COL2    (NTOTAL,NBEAMS)
      DOUBLE PRECISION SCOL2   (2,NBEAMS)

C  Solution constants of integration, and related quantities

      DOUBLE PRECISION LCON(NLAYERS)
      DOUBLE PRECISION MCON(NLAYERS)
      DOUBLE PRECISION LCON_XVEC(2,NLAYERS)
      DOUBLE PRECISION MCON_XVEC(2,NLAYERS)

C  Post-processing variables
C  -------------------------

C  BOA source terms

      DOUBLE PRECISION BOA_SOURCE        ( N_USER_STREAMS )
      DOUBLE PRECISION DIRECT_BOA_SOURCE ( N_USER_STREAMS )

C  Reflectance integrand  a(j).x(j).I(-j)

      DOUBLE PRECISION IDOWNSURF

C  Cumulative source terms

      DOUBLE PRECISION CUMSOURCE_UP(N_USER_STREAMS,0:NLAYERS)
      DOUBLE PRECISION CUMSOURCE_DN(N_USER_STREAMS,0:NLAYERS)

C  Local help variables
C  --------------------

      INTEGER          N, IBEAM

C  local inclusion flags

      LOGICAL          DO_INCLUDE_DIRECTBEAM
      LOGICAL          DO_INCLUDE_SURFACE

C  Flux multiplier and Fourier component numbers

      DOUBLE PRECISION FLUX_MULTIPLIER
      DOUBLE PRECISION DELTA_FACTOR
      DOUBLE PRECISION SURFACE_FACTOR, ALBEDO

C  error tracing

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

      STATUS = 0

C  Set local flags
C  ---------------

C  Surface flag (for inclusion of some kind of reflecting boundary)
C     Lambertian ony - Placeholder for the rest

      IF ( SURFTYPE .EQ. 1 ) THEN
        ALBEDO = LAMBERTIAN_ALBEDO(THREAD)
        DO_INCLUDE_SURFACE = .TRUE.
        IF ( FOURIER .NE. 0 ) THEN
          DO_INCLUDE_SURFACE = .FALSE.
        ELSE
          IF ( ALBEDO .EQ. 0.0d0 ) THEN
            DO_INCLUDE_SURFACE = .FALSE.
          ENDIF
        ENDIF
      ENDIF

C  Direct beam flag (only if above albedo flag has been set)

      IF ( DO_INCLUDE_SURFACE ) THEN
        DO IBEAM = 1, NBEAMS
          DO_REFLECTED_DIRECTBEAM(IBEAM) = .TRUE.
        ENDDO
      ELSE
        DO IBEAM = 1, NBEAMS
          DO_REFLECTED_DIRECTBEAM(IBEAM) = .FALSE.
        ENDDO
      ENDIF

C  surface reflectance factors

      IF ( FOURIER .EQ. 0 ) THEN
        SURFACE_FACTOR = 2.0d0
        DELTA_FACTOR   = 1.0d0
      ELSE
        SURFACE_FACTOR = 1.0d0
        DELTA_FACTOR   = 2.0d0
      ENDIF

C  Flux multipliers
C   = 1 / 4.pi with beam sources,  = 1 for Thermal alone.

      FLUX_MULTIPLIER   = DELTA_FACTOR
      
C  ###################################
C  Set up operations (for Fourier = 0)
C  ###################################

      IF ( FOURIER .EQ. 0 ) THEN

C   MISCSETUPS (4 subroutines)  :
C       Performance Setup,
C       Delta-M scaling,
C       average-secant formulation,
C       transmittance setup

C  Prepare quasi-spherical attenuation

        CALL TWOSTREAM_QSPREP
     I  ( NLAYERS, NBEAMS, DO_PLANE_PARALLEL, 
     I    DELTAU_VERT, CHAPMAN_FACTORS, X0,
     O    INITIAL_TRANS, AVERAGE_SECANT,
     O    LOCAL_SZA, LAYER_PIS_CUTOFF,
     O    DELTAU_SLANT, TAUSLANT, 
     O    SOLAR_BEAM_OPDEP, DO_REFLECTED_DIRECTBEAM )

C  Transmittances and Transmittance factors

        CALL TWOSTREAM_PREPTRANS
     I   ( NLAYERS, N_USER_STREAMS, NBEAMS, DELTAU_VERT, 
     I     STREAM_VALUE, USER_STREAMS,
     I     INITIAL_TRANS, AVERAGE_SECANT, LAYER_PIS_CUTOFF,
     O     T_DELT_MUBAR, T_DELT_USERM,
     O     ITRANS_USERM  ) 

C   EMULT_MASTER  : Beam source function multipliers. Not required for the
C                  Full SS calculation in outgoing mode

        IF (.NOT.DO_SSFULL.OR.(DO_SSFULL.AND.DO_SSCORR_NADIR)) THEN
          CALL TWOSTREAM_EMULTMASTER
     I     ( DO_UPWELLING, DO_DNWELLING,
     I       NLAYERS, NBEAMS, N_USER_STREAMS, DELTAU_VERT, 
     I       USER_STREAMS, T_DELT_MUBAR, T_DELT_USERM, 
     I       ITRANS_USERM, AVERAGE_SECANT, LAYER_PIS_CUTOFF,
     O       SIGMA_M, SIGMA_P, EMULT_HOPRULE,
     O       EMULT_UP, EMULT_DN )
        ENDIF

C  End setups operation

      ENDIF

C  #####################
C  Correction operations 
C  #####################

C  Single scatter correction (pre-calculation)
C     Must be done after MISCSETUPS and EMULT_MASTER, as we need
C      multipliers and transmittance factors for SUN and LOS paths.
C      Code added 6 May 2005. Replaces call in Master routine.
C      Version 3.1. Added call to the new outgoing sphericity correction

      IF ( FOURIER .EQ. 0 ) THEN

C  regular nadir-scattering SS correction

        IF ( DO_SSCORR_NADIR ) THEN
c          CALL TWOSTREAM_SSCORR_NADIR ( SS_FLUX_MULTIPLIER )
        ENDIF

C  New outgoing sphericity correction

        IF ( DO_SSCORR_OUTGOING ) THEN
c          CALL TWOSTREAM_SSCORR_OUTGOING
c     &     ( SS_FLUX_MULTIPLIER, FAIL, MAIL )
         IF ( FAIL ) THEN
            TRACE= 'SS correction outgoing failed'
            STATUS = 1
            CALL TWOSTREAM_ERROR_TRACE ( INIT, MAIL, TRACE, STATUS )
            RETURN
          ENDIF 
        ENDIF

      ENDIF

C  Do Lambertian term if single scatter is being done

      IF ( FOURIER .EQ. 0 ) THEN
        IF ( DO_SSFULL ) THEN
c            CALL TWOSTREAM_LAMBERTIAN_DBCORRECTION (SS_FLUX_MULTIPLIER)
          RETURN
        ENDIF
      ENDIF

C  Reflected Direct beam attenuation

      CALL TWOSTREAM_DIRECTBEAM
     I     ( NBEAMS, N_USER_STREAMS, DO_INCLUDE_SURFACE,
     I       DO_DB_CORRECTION, DELTA_FACTOR, ALBEDO,
     I       SURFTYPE, FLUX_FACTOR, X0,
     I       SOLAR_BEAM_OPDEP, DO_REFLECTED_DIRECTBEAM,
     O       ATMOS_ATTN, DIRECT_BEAM, USER_DIRECT_BEAM )

C  Auxiliary Geometry

      CALL TWOSTREAM_AUXGEOM
     I ( N_USER_STREAMS, NBEAMS, FOURIER, 
     I   X0, USER_STREAMS, STREAM_VALUE,
     O   PX11, PXSQ, POX, PX0X, ULP )

C  ########################################
C  RT differential equation Eigensolutions
C  ########################################

C  Start layer loop

      DO N = 1, NLAYERS

C  Get Discrete ordinate solutions for this layer

        CALL TWOSTREAM_HOM_SOLUTION
     I    ( NLAYERS, N, FOURIER, STREAM_VALUE, PXSQ,
     I      OMEGA_TOTAL(1,THREAD),
     I      ASYMM_TOTAL(1,THREAD),
     I      DELTAU_VERT(1,THREAD), 
     O      SAB, DAB, EIGENVALUE, EIGENTRANS,
     O      XPOS, XNEG )

C  Get Post-processing ("user") solutions for this layer

        CALL TWOSTREAM_HOM_USERSOLUTION
     I    ( NLAYERS, N_USER_STREAMS, N, FOURIER, 
     I      STREAM_VALUE, PX11, USER_STREAMS, ULP, XPOS, XNEG,
     I      OMEGA_TOTAL(1,THREAD),
     I      ASYMM_TOTAL(1,THREAD),
     O      U_XPOS, U_XNEG, U_HELP_P, U_HELP_M )

C  end layer loop

      ENDDO

C  Prepare homogeneous solution multipliers

      CALL TWOSTREAM_HMULT_MASTER
     I     ( NLAYERS, N_USER_STREAMS, USER_STREAMS,
     I       EIGENVALUE, EIGENTRANS, T_DELT_USERM,  
     O       ZETA_M, ZETA_P, HMULT_1, HMULT_2 )

C  ############################################
C   boundary value problem - MATRIX PREPARATION
C  ############################################

C  standard case using compression of band matrices, etc..

      CALL TWOSTREAM_BVP_MATRIXSETUP_MASTER
     I     ( DO_INCLUDE_SURFACE, FOURIER,
     I       SURFACE_FACTOR, ALBEDO, SURFTYPE, NLAYERS, NTOTAL, 
     I       XPOS, XNEG, EIGENTRANS, STREAM_VALUE,
     O       R2_HOMP, R2_HOMM,
     O       BMAT_ROWMASK, BANDMAT2, SMAT2, IPIVOT,
     O       STATUS_SUB )

      IF ( STATUS_SUB .NE. 0 ) THEN
        MAIL = 'Error return from module BVPMATRIX_SETUP_MASTER'
        TRACE= 'Call #2A in TWOSTREAM_FOURIER_MASTER'
        STATUS = 1
        CALL TWOSTREAM_ERROR_TRACE ( INIT, MAIL, TRACE, STATUS )
        RETURN
      ENDIF

C  Start loop over various solar beams

      DO IBEAM = 1, NBEAMS

C  Solar beam Particular solution
C  ------------------------------

C  start layer loop

        DO N = 1, NLAYERS

C  stream solution

          CALL TWOSTREAM_BEAM_SOLUTION
     I    ( NLAYERS, NBEAMS, N, FOURIER, IBEAM,
     I      FLUX_FACTOR, LAYER_PIS_CUTOFF, STREAM_VALUE, X0, PX0X,
     I      AVERAGE_SECANT, INITIAL_TRANS, T_DELT_MUBAR,
     I      OMEGA_TOTAL(1,THREAD),
     I      ASYMM_TOTAL(1,THREAD),
     I      DELTAU_VERT(1,THREAD), 
     I      SAB, DAB, EIGENVALUE,
     O      QSUMVEC, QDIFVEC, QVEC,
     O      WVEC, WUPPER, WLOWER )

C  user solutions

          CALL TWOSTREAM_BEAM_USERSOLUTION
     I    ( DO_UPWELLING, DO_DNWELLING,
     I      NLAYERS, NBEAMS, N_USER_STREAMS, N, FOURIER, IBEAM,
     I      FLUX_FACTOR, LAYER_PIS_CUTOFF, STREAM_VALUE, PX11, X0, POX,
     I      OMEGA_TOTAL(1,THREAD),
     I      ASYMM_TOTAL(1,THREAD),
     I      USER_STREAMS, ULP, WVEC,
     O      W_HELP, U_WPOS1, U_WPOS2, U_WNEG1, U_WNEG2 )

C  end layer loop

        END DO

C  Solve boundary value problem
C  ----------------------------

        CALL TWOSTREAM_BVP_SOLUTION_MASTER
     I       ( DO_INCLUDE_SURFACE, DO_REFLECTED_DIRECTBEAM(IBEAM),
     I         FOURIER, IBEAM, NBEAMS, NLAYERS, NTOTAL,
     I         SURFTYPE, SURFACE_FACTOR, ALBEDO, DIRECT_BEAM,
     I         XPOS, XNEG, WUPPER, WLOWER, STREAM_VALUE,
     I         BANDMAT2, SMAT2, IPIVOT,
     O         R2_PARTIC, COL2, SCOL2, LCON, MCON, 
     O         LCON_XVEC,  MCON_XVEC, 
     O         STATUS_SUB )

        IF ( STATUS_SUB .NE. 0 ) THEN
          MAIL = 'Error return from module BVP_SOLUTION_MASTER'
          TRACE= 'Call #5A in TWOSTREAM_FOURIER_MASTER'
          STATUS = 1
          CALL TWOSTREAM_ERROR_TRACE ( INIT, MAIL, TRACE, STATUS )
          RETURN
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

c          DO_INCLUDE_DIRECTBEAM = ( DO_UPWELLING .AND.
c     &     (DO_REFLECTED_DIRECTBEAM(IBEAM).AND..NOT.DO_DBCORRECTION))

C  former code

        DO_INCLUDE_DIRECTBEAM = ( DO_UPWELLING .AND.
     &     (DO_REFLECTED_DIRECTBEAM(IBEAM)))

C  upwelling

        IF ( DO_UPWELLING ) THEN
          CALL TWOSTREAM_UPUSER_INTENSITY
     I      ( DO_INCLUDE_SURFACE, DO_MSMODE_TWOSTREAM,
     I        DO_INCLUDE_DIRECTBEAM, DO_DB_CORRECTION,
     I        FOURIER, NLAYERS, NBEAMS, N_USER_STREAMS,
     I        IBEAM,  FLUX_MULTIPLIER,
     I        SURFACE_FACTOR, SURFTYPE, ALBEDO, USER_DIRECT_BEAM,
     I        EIGENTRANS, T_DELT_USERM, STREAM_VALUE,
     I        XPOS, XNEG, WUPPER, WLOWER,
     I        LCON, LCON_XVEC, MCON, MCON_XVEC,
     I        U_XPOS, U_XNEG, U_WPOS1, U_WPOS2,
     I        HMULT_1, HMULT_2, EMULT_UP,
     I        BOA_SOURCE, DIRECT_BOA_SOURCE, IDOWNSURF,
     O        INTENSITY_F_UP, CUMSOURCE_UP )
        ENDIF

C  Downwelling

        IF ( DO_DNWELLING ) THEN
          CALL TWOSTREAM_DNUSER_INTENSITY
     I      ( FOURIER, DO_MSMODE_TWOSTREAM,
     I        IBEAM, FLUX_MULTIPLIER,
     I        NLAYERS, NBEAMS, N_USER_STREAMS,
     I        EIGENTRANS, T_DELT_USERM, XPOS, XNEG,
     I        LCON, LCON_XVEC, MCON, MCON_XVEC,
     I        U_XPOS, U_XNEG, U_WNEG1, U_WNEG2,
     I        HMULT_1, HMULT_2, EMULT_DN,
     O        INTENSITY_F_DN,  CUMSOURCE_DN )
        ENDIF

C  End loop over beam solutions

      END DO

C  ######
C  finish
C  ######

      RETURN
      END

C

      SUBROUTINE TWOSTREAM_CONVERGE
     I    ( DO_UPWELLING, DO_DNWELLING, 
     I      DO_SSFULL, DO_SSCORR_OUTGOING, DO_SSCORR_NADIR,
     I      N_GEOMETRIES, NBEAMS, NTHREADS,
     I      N_USER_STREAMS, N_USER_RELAZMS, AZMFAC, UMOFF,
     I      THREAD, IBEAM, FOURIER_COMPONENT,
     I      INTENSITY_F_UP,  INTENSITY_F_DN, 
     I      INTENSITY_SS_UP, INTENSITY_SS_DN,
     O      INTENSITY_TOA, INTENSITY_BOA )

C  input variables
C  ---------------

C  Control

      LOGICAL          DO_UPWELLING, DO_DNWELLING
      LOGICAL          DO_SSFULL, DO_SSCORR_OUTGOING, DO_SSCORR_NADIR

C  Numbers

      INTEGER          N_GEOMETRIES, NBEAMS, NTHREADS
      INTEGER          N_USER_STREAMS, N_USER_RELAZMS

C  FOurier component and thread, beam

      INTEGER          FOURIER_COMPONENT, THREAD, IBEAM

C  Local  azimuth factors

      INTEGER          UMOFF ( NBEAMS, N_USER_STREAMS )
      DOUBLE PRECISION AZMFAC
     &     (N_USER_STREAMS,NBEAMS,N_USER_RELAZMS)

C  User-defined solutions

      DOUBLE PRECISION INTENSITY_F_UP
     D   (N_USER_STREAMS,NBEAMS)
      DOUBLE PRECISION INTENSITY_F_DN
     D   (N_USER_STREAMS,NBEAMS)

C  Single scatter solutions

      DOUBLE PRECISION INTENSITY_SS_UP(N_GEOMETRIES)
      DOUBLE PRECISION INTENSITY_SS_DN(N_GEOMETRIES)

C  Output
C  ------

      DOUBLE PRECISION INTENSITY_TOA(N_GEOMETRIES,NTHREADS)
      DOUBLE PRECISION INTENSITY_BOA(N_GEOMETRIES,NTHREADS)

C  local variables
C  ---------------

      INTEGER          I, UA, V
      DOUBLE PRECISION TOLD, TAZM

C  ###################
C  Fourier 0 component
C  ###################

      IF ( FOURIER_COMPONENT.EQ.0 ) THEN

C  Copy DIFFUSE Fourier component at all output angles and optical depths
C    If no SSCORR and no DBCORR, then two options apply:
C     (a) Convergence on RADIANCE = DIFFUSE + SSTRUNCATED + DBTRUNCATED
C              (full radiance, no SS correction, no DB correction)
C     (b) Convergence on RADIANCE = DIFFUSE alone (MS only mode)
C              (SSTRUNCATED + DBTRUNCATED do not get calculated)

        IF ( .not. DO_SSFULL ) THEN
          DO I = 1, N_USER_STREAMS
            DO UA = 1, N_USER_RELAZMS
              V = UMOFF(IBEAM,I) + UA
              IF ( DO_UPWELLING ) THEN
                INTENSITY_TOA(V,THREAD) = INTENSITY_F_UP(I,IBEAM)
              ENDIF
              IF ( DO_DNWELLING ) THEN
                INTENSITY_BOA(V,THREAD) = INTENSITY_F_DN(I,IBEAM)
              ENDIF
            ENDDO
          ENDDO
        ELSE
           DO I = 1, N_USER_STREAMS
            DO UA = 1, N_USER_RELAZMS
              V = UMOFF(IBEAM,I) + UA
              IF ( DO_UPWELLING ) THEN
                INTENSITY_TOA(V,THREAD) = 0.0d0
              ENDIF
              IF ( DO_DNWELLING ) THEN
                INTENSITY_BOA(V,THREAD) = 0.0d0
              ENDIF
            ENDDO
          ENDDO
        ENDIF

C    Add the single scatter component if flagged

        IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
           DO I = 1, N_USER_STREAMS
            DO UA = 1, N_USER_RELAZMS
              V = UMOFF(IBEAM,I) + UA
              IF ( DO_UPWELLING ) THEN
                INTENSITY_TOA(V,THREAD) = 
     &             INTENSITY_TOA(V,THREAD) + INTENSITY_SS_UP(V)
              ENDIF
              IF ( DO_DNWELLING ) THEN
                INTENSITY_BOA(V,THREAD) = 
     &             INTENSITY_BOA(V,THREAD) + INTENSITY_SS_DN(V)
              ENDIF
            ENDDO
          ENDDO
       ENDIF

C  ######################
C  Fourier component = 1
C  ######################

      ELSE

C  No examination of convergence
C  -----------------------------

        DO UA = 1, N_USER_RELAZMS
          DO I = 1, N_USER_STREAMS
            V = UMOFF(IBEAM,I) + UA
            IF ( DO_UPWELLING ) THEN
              TOLD = INTENSITY_TOA(V,THREAD)
              TAZM = AZMFAC(I,IBEAM,UA)*INTENSITY_F_UP(I,IBEAM)
              INTENSITY_TOA(V,THREAD) = TOLD + TAZM
            ENDIF
            IF ( DO_DNWELLING ) THEN
              TOLD = INTENSITY_BOA(V,THREAD)
              TAZM = AZMFAC(I,IBEAM,UA)*INTENSITY_F_DN(I,IBEAM)
              INTENSITY_BOA(V,THREAD) = TOLD + TAZM
            ENDIF
          ENDDO
        ENDDO

C  Finish Fourier
   
      ENDIF

C  Finish

      RETURN
      END
