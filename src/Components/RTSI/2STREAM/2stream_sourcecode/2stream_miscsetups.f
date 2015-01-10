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
C #              TWOSTREAM_AUXGEOM                              #
C #              TWOSTREAM_QSPREP                               #
C #              TWOSTREAM_PREPTRANS                            #
C #              TWOSTREAM_DIRECTBEAM                           #
C #              TWOSTREAM_EMULTMASTER                          #
C #                                                             #
C #              TWOSTREAM_BEAM_GEOMETRY_PREPARE                #
C #                                                             #
C ###############################################################

      SUBROUTINE TWOSTREAM_AUXGEOM
     I ( N_USER_STREAMS, NBEAMS, FOURIER,
     I   X0, USER_STREAMS, STREAM_VALUE,
     O   PX11, PXSQ, POX, PX0X, ULP )

C  Input arguments
C  ---------------

C  Numbers

      INTEGER          N_USER_STREAMS, NBEAMS

C  Fourier

      INTEGER          FOURIER

C  stream directions

      DOUBLE PRECISION X0 ( NBEAMS )
      DOUBLE PRECISION USER_STREAMS ( N_USER_STREAMS )
      DOUBLE PRECISION STREAM_VALUE

C  Output
C  ------

      DOUBLE PRECISION ULP ( N_USER_STREAMS )
      DOUBLE PRECISION POX  ( NBEAMS )
      DOUBLE PRECISION PX0X ( NBEAMS )
      DOUBLE PRECISION PXSQ, PX11

C  Local variables

      INTEGER          UM, IBEAM
      DOUBLE PRECISION MU, MU0

C  Saved quantities

      DO UM = 1, N_USER_STREAMS
        MU = USER_STREAMS(UM)
        ULP(UM) =  -DSQRT(0.5d0*(1.0d0-MU*MU))
      ENDDO

      if ( fourier.eq.0) then
        PXSQ = STREAM_VALUE * STREAM_VALUE
      else if ( fourier.eq.1 ) then
        PX11 = DSQRT(0.5d0*(1.0d0-STREAM_VALUE*STREAM_VALUE))
        PXSQ = PX11 * PX11
      endif

      DO IBEAM = 1, NBEAMS
        MU0 = x0(ibeam)
        POX(IBEAM) = DSQRT(0.5d0*(1.0d0-MU0*MU0))
        if ( fourier.eq.0) then
          PX0X(IBEAM) = X0(IBEAM)  * STREAM_VALUE
        else
          PX0X(IBEAM) = POX(IBEAM) * PX11
        endif
      ENDDO

C  FInish

      RETURN
      END

C

      SUBROUTINE TWOSTREAM_QSPREP
     I  ( NLAYERS, NBEAMS, DO_PLANE_PARALLEL, 
     I    DELTAU_VERT, CHAPMAN_FACTORS, X0,
     O    INITIAL_TRANS, AVERAGE_SECANT,
     O    LOCAL_SZA, LAYER_PIS_CUTOFF,
     O    DELTAU_SLANT, TAUSLANT, 
     O    SOLAR_BEAM_OPDEP, DO_REFLECTED_DIRECTBEAM )

C  Inputs
C  ------

C  Control

      INTEGER          NLAYERS, NBEAMS

C  optical thickness input

      DOUBLE PRECISION DELTAU_VERT ( NLAYERS )

C  Path segments distances (km)

      DOUBLE PRECISION CHAPMAN_FACTORS(NLAYERS,NLAYERS,NBEAMS)

C  Plane parallel control

      LOGICAL          DO_PLANE_PARALLEL

C  Beam SZA cosines

      DOUBLE PRECISION X0(NBEAMS)

C  Output
C  ------

C  Last layer to include Particular integral solution

      INTEGER          LAYER_PIS_CUTOFF(NBEAMS)

C  Average-secant and initial tramsittance factors for solar beams.

      DOUBLE PRECISION
     &     INITIAL_TRANS  ( NLAYERS, NBEAMS ),
     &     AVERAGE_SECANT ( NLAYERS, NBEAMS ),
     &     LOCAL_SZA      ( NLAYERS, NBEAMS )

C  Derived optical thickness inputs

      DOUBLE PRECISION TAUSLANT     ( 0:NLAYERS, NBEAMS )
      DOUBLE PRECISION DELTAU_SLANT ( NLAYERS, NLAYERS, NBEAMS )

C  Solar beam attenuations and reflectance flags

      DOUBLE PRECISION SOLAR_BEAM_OPDEP        ( NBEAMS )
      LOGICAL          DO_REFLECTED_DIRECTBEAM ( NBEAMS )

C  Local variables
C  ---------------

      INTEGER          N, K, IB
      DOUBLE PRECISION S_T_0, S_T_1, SEC0, TAU, DELS
      DOUBLE PRECISION TAUGRID (0:NLAYERS), TAU_SOLAR(NBEAMS)
      DOUBLE PRECISION MAX_TAU_PATH
      PARAMETER        ( MAX_TAU_PATH = 88.0d0 )

C  Complete grid

      TAUGRID(0) = 0.0d0
      DO N = 1, NLAYERS
        TAUGRID(N) = TAUGRID(N-1) + DELTAU_VERT(N)
      ENDDO

C  slant optical thickness values

      DO IB = 1, NBEAMS
        DO N = 1, NLAYERS
          DO K = 1, N
            DELS = CHAPMAN_FACTORS(N,K,IB)
            DELTAU_SLANT(N,K,IB) = DELTAU_VERT(K) * DELS
          ENDDO
        ENDDO
      ENDDO

C  plane-parallel case
C  -------------------

      IF ( DO_PLANE_PARALLEL ) THEN

       S_T_0 = 1.0d0
       DO IB = 1, NBEAMS
        SEC0 = 1.0d0 / X0(IB)
        LAYER_PIS_CUTOFF(IB) = NLAYERS
        DO N = 1, NLAYERS
          TAUSLANT(N,IB) = TAUGRID(N) * SEC0
          IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
            IF ( TAUSLANT(N,IB) .GT. MAX_TAU_PATH ) THEN
              LAYER_PIS_CUTOFF(IB) = N
            ENDIF
            AVERAGE_SECANT(N,IB) = SEC0
            INITIAL_TRANS(N,IB)  = DEXP ( - TAUGRID(N-1) * SEC0 )
            LOCAL_SZA(N,IB)      = X0(IB)
          ELSE
            AVERAGE_SECANT(N,IB) = 0.0d0
            INITIAL_TRANS(N,IB)  = 0.0d0
            LOCAL_SZA(N,IB)      = 0.0d0
          ENDIF
        ENDDO
        TAU_SOLAR(IB) = TAUSLANT(NLAYERS,IB)
       ENDDO

      ELSE

C  pseudo-spherical case
C  ---------------------

       DO IB = 1, NBEAMS

C  Get the total spherical attenuation from layer thickness sums

        TAUSLANT(0,IB) = 0.0d0
        DO N = 1, NLAYERS
          TAU = 0.0d0
          DO K = 1, N
            TAU = TAU + DELTAU_SLANT(N,K,IB)
          ENDDO
          TAUSLANT(N,IB) = TAU
        ENDDO
        TAU_SOLAR(IB) = TAUSLANT(NLAYERS,IB)

C  set up the average secant formulation

        S_T_0 = 1.0d0
        S_T_1 = 0.0d0
        LAYER_PIS_CUTOFF(IB) = NLAYERS
        DO N = 1, NLAYERS
          IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
            IF ( TAUSLANT(N,IB) .GT. MAX_TAU_PATH ) THEN
              LAYER_PIS_CUTOFF(IB) = N
            ELSE
              S_T_1 = DEXP ( - TAUSLANT(N,IB) )
            ENDIF
            AVERAGE_SECANT(N,IB) = (TAUSLANT(N,IB)-TAUSLANT(N-1,IB)) 
     &                               / DELTAU_VERT(N)
            INITIAL_TRANS(N,IB)  = S_T_0
            LOCAL_SZA(N,IB)      = X0(IB)
            S_T_0             = S_T_1
          ELSE
            AVERAGE_SECANT(N,IB) = 0.0d0
            INITIAL_TRANS(N,IB)  = 0.0d0
          ENDIF
        ENDDO

C  Set the Local solar zenith angles

        DO N = 1, NLAYERS
          IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
            LOCAL_SZA(N,IB) = X0(IB)
          ELSE
            LOCAL_SZA(N,IB) = 0.0d0
          ENDIF
        ENDDO

       ENDDO
      ENDIF

C  debug
c      do N = 1, nlayers
c      write(*,'(i3,1p2e20.10)')n,initial_trans(n,1),average_secant(n,1)
c      enddo

C  Set Direct Beam Flag and solar beam total attenuation to surface

      DO IB = 1, NBEAMS
        IF ( TAU_SOLAR(IB) .GT. MAX_TAU_PATH ) THEN
          SOLAR_BEAM_OPDEP(IB) = 0.0d0
          DO_REFLECTED_DIRECTBEAM(IB) = .FALSE.
        ELSE
          SOLAR_BEAM_OPDEP(IB) = DEXP( - TAU_SOLAR(IB) )
          DO_REFLECTED_DIRECTBEAM(IB) = .TRUE.
        ENDIF
      ENDDO

C  finish

      RETURN
      END

C

      SUBROUTINE TWOSTREAM_PREPTRANS
     I   ( NLAYERS, N_USER_STREAMS, NBEAMS, DELTAU_VERT, 
     I     STREAM_VALUE, USER_STREAMS,
     I     INITIAL_TRANS, AVERAGE_SECANT, LAYER_PIS_CUTOFF,
     O     T_DELT_MUBAR, T_DELT_USERM,
     O     ITRANS_USERM  )     

C  Prepare transmittances and transmittance factors

C  Inputs
C  ------

C  Control

      INTEGER          NLAYERS, N_USER_STREAMS, NBEAMS

C  Input optical depths after delta-M scaling and Chapman function

      DOUBLE PRECISION DELTAU_VERT    ( NLAYERS )

C  Stream value

      DOUBLE PRECISION STREAM_VALUE

C  User streams

      DOUBLE PRECISION USER_STREAMS ( N_USER_STREAMS )

C  Last layer to include Particular integral solution

      INTEGER          LAYER_PIS_CUTOFF(NBEAMS)

C  Average-secant and initial tramsittance factors for solar beams.

      DOUBLE PRECISION
     &     INITIAL_TRANS  ( NLAYERS, NBEAMS ),
     &     AVERAGE_SECANT ( NLAYERS, NBEAMS )

C  Outputs
C  -------

C  Transmittance factors for average secant stream

      DOUBLE PRECISION T_DELT_MUBAR ( NLAYERS, NBEAMS )

C  Transmittance factors for user-defined stream angles

      DOUBLE PRECISION T_DELT_USERM ( NLAYERS, N_USER_STREAMS )
      DOUBLE PRECISION
     &     ITRANS_USERM ( NLAYERS, N_USER_STREAMS, NBEAMS )

C  local variables
C  ---------------

      INTEGER          N, UM, IB
      DOUBLE PRECISION SPHER
      DOUBLE PRECISION MAX_TAU_PATH
      PARAMETER        ( MAX_TAU_PATH = 88.0d0 )

C  Transmittance factors for average secant streams
C  ================================================

      DO IB = 1, NBEAMS
        DO N = 1, NLAYERS
          IF ( N. GT. LAYER_PIS_CUTOFF(IB) ) THEN
            T_DELT_MUBAR(N,IB) = 0.0d0
          ELSE
            SPHER = DELTAU_VERT(N) * AVERAGE_SECANT(N,IB)
            IF ( SPHER .GT. MAX_TAU_PATH ) THEN
              T_DELT_MUBAR(N,IB) = 0.0d0
            ELSE
              T_DELT_MUBAR(N,IB) = DEXP ( - SPHER )
            ENDIF
          ENDIF
        ENDDO
      ENDDO

C  Transmittances for User Streams
C  ===============================

C  Initial transmittances divided by user streams

      DO IB = 1, NBEAMS
        DO N = 1, NLAYERS
         DO UM = 1, N_USER_STREAMS
          ITRANS_USERM(N,UM,IB) = INITIAL_TRANS(N,IB)/USER_STREAMS(UM)
         ENDDO
        ENDDO
      ENDDO

C  Whole Layer transmittances

      DO N = 1, NLAYERS
        DO UM = 1, N_USER_STREAMS
          SPHER = DELTAU_VERT(N) / USER_STREAMS(UM)
          IF ( SPHER.GT.MAX_TAU_PATH ) THEN
            T_DELT_USERM(N,UM) = 0.0d0
          ELSE
            T_DELT_USERM(N,UM) = DEXP ( - SPHER )
          ENDIF
        ENDDO
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE TWOSTREAM_DIRECTBEAM
     I     ( NBEAMS, N_USER_STREAMS, DO_INCLUDE_SURFACE,
     I       DELTA_FACTOR, ALBEDO,
     I       SURFTYPE, FLUX_FACTOR, X0,
     I       SOLAR_BEAM_OPDEP, DO_REFLECTED_DIRECTBEAM,
     O       ATMOS_ATTN, DIRECT_BEAM, USER_DIRECT_BEAM )

C  input arguments
C  ---------------

C  Numbers

      INTEGER          NBEAMS, N_USER_STREAMS

C  DB Control, not required in this version
c      LOGICAL          DO_DB_CORRECTION

C  Control

      LOGICAL          DO_INCLUDE_SURFACE

C  Surface inputs

      DOUBLE PRECISION ALBEDO
      DOUBLE PRECISION DELTA_FACTOR
      INTEGER          SURFTYPE

C  FLux factor

      DOUBLE PRECISION FLUX_FACTOR

C  Solar beams, cosines

      DOUBLE PRECISIOn X0(NBEAMS)

C  Solar beam attenuations and reflectance flags

      DOUBLE PRECISION SOLAR_BEAM_OPDEP        ( NBEAMS )
      LOGICAL          DO_REFLECTED_DIRECTBEAM ( NBEAMS )

C  output arguments
C  ----------------

C  Atmospheric attenuation

      DOUBLE PRECISION ATMOS_ATTN ( NBEAMS )

C  Direct beam solutions

      DOUBLE PRECISION
     &        DIRECT_BEAM      ( NBEAMS ),
     &        USER_DIRECT_BEAM ( N_USER_STREAMS, NBEAMS )

C  Local variables
C  ---------------

      DOUBLE PRECISION PI4, X0_FLUX, ATTN, REFL_ATTN
      INTEGER          UI, IB

C  Initialize
C  ----------

C   Safety first!  Return if there is no reflection.

      DO IB = 1, NBEAMS
        DIRECT_BEAM(IB) = 0.0d0
        DO UI = 1, N_USER_STREAMS
          USER_DIRECT_BEAM(UI,IB) = 0.0d0
        ENDDO
      ENDDO

C  Return if no surface

      IF ( .NOT.DO_INCLUDE_SURFACE ) RETURN

C  Attenuation of solar beam
C  -------------------------

      PI4 = DACOS(-1.0d0) * 4.0d0
      DO IB = 1, NBEAMS
       IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN

C  There should be no flux factor here.
C    Bug fixed 18 November 2005. Earlier Italian Job!!
C    Flux Factor put back, 1 March 2007. Using 1 / pi4

        X0_FLUX        = 4.0d0 * X0(IB) / DELTA_FACTOR
        X0_FLUX        = FLUX_FACTOR * X0_FLUX / PI4
        ATTN           = X0_FLUX * SOLAR_BEAM_OPDEP(IB)
        ATMOS_ATTN(IB) = ATTN

C  Total contributions, Lambertian case

        IF ( SURFTYPE .EQ. 1 ) THEN
          REFL_ATTN      = ATTN * ALBEDO
          DIRECT_BEAM(IB) = REFL_ATTN
c          IF ( .NOT.DO_DB_CORRECTION ) THEN
            DO UI = 1, N_USER_STREAMS
              USER_DIRECT_BEAM(UI,IB) = REFL_ATTN
            ENDDO
c          ENDIF
        ENDIF

C  Total contributions, Surface Type BRDF

        IF ( SURFTYPE .NE. 1 ) THEN
C  Placeholder
        ENDIF

C  end direct beam calculation

       ENDIF
      ENDDO

C  finish

      RETURN
      END

C

      SUBROUTINE TWOSTREAM_EMULTMASTER
     I     ( DO_UPWELLING, DO_DNWELLING,
     I       NLAYERS, NBEAMS, N_USER_STREAMS, DELTAU_VERT, 
     I       USER_STREAMS, T_DELT_MUBAR, T_DELT_USERM, 
     I       ITRANS_USERM, AVERAGE_SECANT, LAYER_PIS_CUTOFF,
     O       SIGMA_M, SIGMA_P, EMULT_HOPRULE,
     O       EMULT_UP, EMULT_DN )

C  Prepare multipliers for the Beam source terms

C  Input arguments
C  ===============

C  Numbers

      INTEGER          NLAYERS, NBEAMS, N_USER_STREAMS

C  Control

      LOGICAL          DO_UPWELLING, DO_DNWELLING

C  Layer optical thickness

      DOUBLE PRECISION DELTAU_VERT ( NLAYERS )

C  User streams

      DOUBLE PRECISION USER_STREAMS ( N_USER_STREAMS )

C  Transmittance factors for user-defined stream angles
C    Computed in the initial setup stage for Fourier m = 0

      DOUBLE PRECISION
     &     T_DELT_USERM ( NLAYERS, N_USER_STREAMS )

C  Transmittance factors for average secant stream

      DOUBLE PRECISION
     &     T_DELT_MUBAR ( NLAYERS, NBEAMS )

C  Average-secant and initial tramsittance factors for solar beams.

      DOUBLE PRECISION
     &     ITRANS_USERM   ( NLAYERS, N_USER_STREAMS, NBEAMS ),
     &     AVERAGE_SECANT ( NLAYERS, NBEAMS )

C  Last layer to include Particular integral solution

      INTEGER          LAYER_PIS_CUTOFF(NBEAMS)

C  Output = Global multipliers
C  ===========================

C  coefficient functions for user-defined angles

      DOUBLE PRECISION 
     &      SIGMA_M(NLAYERS,N_USER_STREAMS,NBEAMS),
     &      SIGMA_P(NLAYERS,N_USER_STREAMS,NBEAMS)
      
C  forcing term multipliers (saved for whole atmosphere)

      DOUBLE PRECISION EMULT_UP
     &       (N_USER_STREAMS,NLAYERS,NBEAMS)

      DOUBLE PRECISION EMULT_DN
     &       (N_USER_STREAMS,NLAYERS,NBEAMS)

C  L'Hopital's rule logical variables

      LOGICAL          EMULT_HOPRULE
     &       (NLAYERS,N_USER_STREAMS,NBEAMS)

C  local variables
C  ---------------

      INTEGER          N, UM, IB
      DOUBLE PRECISION WDEL, WUDEL, UDEL
      DOUBLE PRECISION DIFF, SB, SU, SD, SM
      DOUBLE PRECISION HOPITAL_TOLERANCE
      PARAMETER        ( HOPITAL_TOLERANCE = 0.001d0 )

C  L'Hopital's Rule flags for Downwelling EMULT
C  --------------------------------------------

      IF ( DO_DNWELLING ) THEN
       DO N = 1, NLAYERS
         DO IB = 1, NBEAMS
          SB = AVERAGE_SECANT(N,IB)
          DO UM = 1, N_USER_STREAMS
            SM = 1.0d0 / USER_STREAMS(UM)
            DIFF = DABS ( SM - SB )
            IF ( DIFF .LT. HOPITAL_TOLERANCE ) THEN
              EMULT_HOPRULE(N,UM,IB) = .TRUE.
            ELSE
              EMULT_HOPRULE(N,UM,IB) = .FALSE.
            ENDIF
          ENDDO
         ENDDO
       ENDDO
      ENDIF

C  sigma functions (all layers)
C  ----------------------------

      DO N = 1, NLAYERS
       DO IB = 1, NBEAMS
        SB = AVERAGE_SECANT(N,IB)
        DO UM = 1, N_USER_STREAMS
          SM = 1.0d0 / USER_STREAMS(UM)
          SIGMA_P(N,UM,IB) = SB + SM
          SIGMA_M(N,UM,IB) = SB - SM
        ENDDO
       ENDDO
      ENDDO

C  upwelling External source function multipliers
C  ----------------------------------------------

      IF ( DO_UPWELLING ) THEN
        DO N = 1, NLAYERS
          DO IB = 1, NBEAMS
            IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
              DO UM = 1, N_USER_STREAMS
                EMULT_UP(UM,N,IB) = 0.0d0
              ENDDO
            ELSE
              WDEL = T_DELT_MUBAR(N,IB)
              DO UM = 1, N_USER_STREAMS
                WUDEL = WDEL * T_DELT_USERM(N,UM)
                SU = ( 1.0d0 - WUDEL ) / SIGMA_P(N,UM,IB)
                EMULT_UP(UM,N,IB) = ITRANS_USERM(N,UM,IB) * SU
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDIF

C  debug
c      do n = 1, nlayers
c        write(57,'(i3,1pe24.12)')n,EMULT_UP(1,N,1)
c        write(57,'(i3,1pe24.12)')n,T_DELT_MUBAR(N,1)
c      enddo

C  downwelling External source function multipliers
C  ------------------------------------------------

C    .. Note use of L'Hopitals Rule

      IF ( DO_DNWELLING ) THEN
        DO N = 1, NLAYERS
          DO IB = 1, NBEAMS
            IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
              DO UM = 1, N_USER_STREAMS
                EMULT_DN(UM,N,IB) = 0.0d0
              ENDDO
            ELSE
              WDEL = T_DELT_MUBAR(N,IB)
              DO UM = 1, N_USER_STREAMS
                UDEL = T_DELT_USERM(N,UM)
                IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
                  SD = DELTAU_VERT(N) * UDEL
                ELSE
                  SD = ( UDEL - WDEL ) / SIGMA_M(N,UM,IB)
                ENDIF 
                EMULT_DN(UM,N,IB) = ITRANS_USERM(N,UM,IB) * SD
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDIF

C  debug
c      do n = 1, nlayers
c        write(57,'(i3,1pe24.12)')n,EMULT_DN(1,N,1)
c      enddo

C  Finish

      RETURN
      END


      SUBROUTINE TWOSTREAM_BEAM_GEOMETRY_PREPARE
     I     ( DO_PLANE_PARALLEL, NBEAMS, NLAYERS, IBEAM,
     I       SZA_GEOM_TRUE, REARTH, HEIGHTS,
     O       CHAPMAN_FACTORS, SZA_LEVEL_OUTPUT )
        
C  Generate path CHAPMAN_FACTORS and SZA angles SZA_LEVEL_OUTPUT
C  for a curved ray-traced beam through a multilayer atmosphere.


C  Input arguments
C  ===============

C  flag for plane parallel case

      LOGICAL          DO_PLANE_PARALLEL
        
C  Beam index

      INTEGER          IBEAM

C  number of beams and layers

      INTEGER          NBEAMS, NLAYERS

C  True solar zenith angle (degrees)

      DOUBLE PRECISION SZA_GEOM_TRUE

C  Earth radius (km)

      DOUBLE PRECISION REARTH
        
C  Coarse grids of heights, pressures and temperatures

      DOUBLE PRECISION HEIGHTS     (0:NLAYERS)

C  Output arguments
C  ================

C  Path segments distances (km)

      DOUBLE PRECISION CHAPMAN_FACTORS(NLAYERS,NLAYERS,NBEAMS)

C  solar zenith angles at nadir

      DOUBLE PRECISION SZA_LEVEL_OUTPUT(0:NLAYERS,NBEAMS)

C  Local variables
C  ===============

C  local height arrays

      DOUBLE PRECISION H(0:NLAYERS)
      DOUBLE PRECISION DELZ(NLAYERS)

C  help variables

      INTEGER          N, K, IB
      DOUBLE PRECISION GM_TOA, TH_TOA, MU_TOA
      DOUBLE PRECISION DEG_TO_RAD
      DOUBLE PRECISION 
     &      STH1, SINTH1, STH2, SINTH2, PHI, SINPHI, 
     &      RE_LOWER, RE_UPPER, DIST, STH2D

C  Some setup operations
C  =====================

C  initialise output

      IB = IBEAM
      SZA_LEVEL_OUTPUT(0,IB) = 0.0D0
      DO N = 1, NLAYERS
        SZA_LEVEL_OUTPUT(N,IB) = 0.0D0
        DO K = 1, NLAYERS
          CHAPMAN_FACTORS(N,K,IB) = 0.0D0
        ENDDO
      ENDDO
        
C  earth radii and heights differences

      DO N = 0, NLAYERS
        H(N) = HEIGHTS(N) + REARTH
      ENDDO

      DO N = 1, NLAYERS
        DELZ(N) = HEIGHTS(N-1)-HEIGHTS(N)
      ENDDO

C  TOA values

      SZA_LEVEL_OUTPUT(0,IB) = SZA_GEOM_TRUE
      DEG_TO_RAD = DATAN(1.0D0) / 45.0D0
      TH_TOA = SZA_GEOM_TRUE * DEG_TO_RAD
      MU_TOA = DCOS(TH_TOA)
      GM_TOA = DSQRT ( 1.0D0 - MU_TOA * MU_TOA )

C  initialize

      STH2D  = 0.0D0
        
C  plane-parallel case
C  ===================

      IF ( DO_PLANE_PARALLEL ) THEN
        DO N = 1, NLAYERS
          SZA_LEVEL_OUTPUT(N,IB) = SZA_GEOM_TRUE
          DO K = 1, N
            CHAPMAN_FACTORS(N,K,IB) = 1.0D0 / MU_TOA
          ENDDO
        ENDDO
        RETURN
      ENDIF

C  Straight line geometry
C  ======================

      DO N = 1, NLAYERS

C  start values

        SINTH1 = GM_TOA * H(N) / H(0)
        STH1   = DASIN(SINTH1)
        RE_UPPER = H(0)

C  solar zenith angles are all the same = input value

        SZA_LEVEL_OUTPUT(N,IB) = SZA_GEOM_TRUE

C loop over layers K from 1 to layer N

        DO K = 1, N

C  sine-rule; PHI = earth-centered angle

          RE_LOWER = RE_UPPER - DELZ(K)
          SINTH2 = RE_UPPER * SINTH1 / RE_LOWER
          STH2   = DASIN(SINTH2)
          PHI    = STH2 - STH1
          SINPHI = DSIN(PHI)
          DIST = RE_UPPER * SINPHI / SINTH2
          CHAPMAN_FACTORS(N,K,IB) = DIST / DELZ(K)

C  re-set

          RE_UPPER = RE_LOWER
          SINTH1 = SINTH2
          STH1   = STH2

        ENDDO

C  finish main layer loop

      ENDDO

C  end of routine

      RETURN
      END
