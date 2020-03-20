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
C #              TWOSTREAM_L_QSPREP                             #
C #              TWOSTREAM_L_EMULTMASTER                        #
C #                TWOSTREAM_L_EMULT_UP                         #
C #                TWOSTREAM_L_EMULT_DN                         #
C #                                                             #
C ###############################################################

      SUBROUTINE TWOSTREAM_L_QSPREP
     I  ( NLAYERS, NBEAMS, NPARS, N_USER_STREAMS,
     I    DO_PLANE_PARALLEL, DO_COLUMN_WFS, DO_PROFILE_WFS,
     I    LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_COLUMN_WFS,
     I    DELTAU_VERT, L_DELTAU_VERT, CHAPMAN_FACTORS, 
     I    X0, USER_STREAMS, T_DELT_USERM, LAYER_PIS_CUTOFF,
     I    INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,
     O    L_INITIAL_TRANS, L_AVERAGE_SECANT, 
     O    L_T_DELT_MUBAR, L_T_DELT_USERM )

C  Inputs
C  ------

C  Numbers

      INTEGER          NLAYERS, NBEAMS, NPARS, N_USER_STREAMS

C  optical thickness inputs

      DOUBLE PRECISION DELTAU_VERT   ( NLAYERS )
      DOUBLE PRECISION L_DELTAU_VERT ( NLAYERS, NPARS )

C  Chapman factors

      DOUBLE PRECISION CHAPMAN_FACTORS ( NLAYERS, NLAYERS, NBEAMS )

C  Plane parallel control

      LOGICAL          DO_PLANE_PARALLEL

C  LInearization flags

      LOGICAL          DO_COLUMN_WFS
      LOGICAL          DO_PROFILE_WFS

C  Linearization control

      LOGICAL          LAYER_VARY_FLAG(NLAYERS)
      INTEGER          LAYER_VARY_NUMBER(NLAYERS)
      INTEGER          N_COLUMN_WFS

C  Beam SZA cosines

      DOUBLE PRECISION X0(NBEAMS)

C  Last layer to include Particular integral solution

      INTEGER          LAYER_PIS_CUTOFF(NBEAMS)

C  Average-secant and initial tramsittance factors for solar beams.

      DOUBLE PRECISION
     &     INITIAL_TRANS  ( NLAYERS, NBEAMS ),
     &     AVERAGE_SECANT ( NLAYERS, NBEAMS )

C  User streams

      DOUBLE PRECISION USER_STREAMS ( N_USER_STREAMS )

C  Transmittance factors for average secant stream

      DOUBLE PRECISION T_DELT_MUBAR ( NLAYERS, NBEAMS )

C  Transmittance factors for user-defined stream angles

      DOUBLE PRECISION T_DELT_USERM ( NLAYERS, N_USER_STREAMS )

C  Output
C  ------

C  Linearized Average-secant and initial tramsittance factors

      DOUBLE PRECISION
     &     L_INITIAL_TRANS  ( NLAYERS, NBEAMS, 0:NLAYERS, NPARS ),
     &     L_AVERAGE_SECANT ( NLAYERS, NBEAMS, 0:NLAYERS, NPARS )

C  Transmittance factors for average secant stream

      DOUBLE PRECISION L_T_DELT_MUBAR
     *        ( NLAYERS, NBEAMS, 0:NLAYERS, NPARS )

C  Transmittance factors for user-defined stream angles

      DOUBLE PRECISION L_T_DELT_USERM 
     *        ( NLAYERS, N_USER_STREAMS, NPARS )

C  Local variables
C  ---------------

      INTEGER          N, K, IB, Q, UM
      DOUBLE PRECISION MAX_TAU_PATH, SUM, DELT, LAMDA, FAC
      DOUBLE PRECISION VAR, WDEL, RHO, TRANS, CF
      PARAMETER        ( MAX_TAU_PATH = 88.0d0 )

C  linearization of Initial transmittances
C  =======================================

C         Use Logarithmic derivative !!!!
C         Reason: avoids exceptions if INITIAL_TRANS underflows

C  Profile linearization

      IF ( DO_PROFILE_WFS ) THEN
       DO IB = 1, NBEAMS
        DO N = 1, NLAYERS
          IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_INITIAL_TRANS(N,IB,N,Q) = 0.0d0
            ENDDO
            IF ( N .GT. 1 ) THEN
              DO K = 1, N-1
                DO Q = 1, LAYER_VARY_NUMBER(K)
                  CF = CHAPMAN_FACTORS(N-1,K,IB)
                  FAC = - L_DELTAU_VERT(K,Q) * CF
!                  FAC = FAC * INITIAL_TRANS(N,IB)   ! Non-log Derivative
                  L_INITIAL_TRANS(N,IB,K,Q) = FAC
                ENDDO
              ENDDO
            ENDIF
          ELSE
            DO K = 1, N
              DO Q = 1, LAYER_VARY_NUMBER(K)
                L_INITIAL_TRANS(N,IB,K,Q) = 0.0d0
              ENDDO
            ENDDO
          ENDIF
        ENDDO
       ENDDO
      ENDIF

C  Column weighting functions

      IF ( DO_COLUMN_WFS ) THEN
        DO IB = 1, NBEAMS
          N = 1
          DO Q = 1, N_COLUMN_WFS
            L_INITIAL_TRANS(N,IB,0,Q) = 0.0d0
          ENDDO
          DO N = 2, NLAYERS
            IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
              DO Q = 1, N_COLUMN_WFS
                SUM = 0.0d0
                DO K = 1, N-1
                  CF = CHAPMAN_FACTORS(N-1,K,IB)
                  SUM = SUM + L_DELTAU_VERT(K,Q) * CF
                ENDDO
!                SUM = SUM * INITIAL_TRANS(N,IB)   ! Non-log derivative
                L_INITIAL_TRANS(N,IB,0,Q) = - SUM
              ENDDO
            ELSE
              DO Q = 1, N_COLUMN_WFS
                L_INITIAL_TRANS(N,IB,0,Q) = 0.0d0
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDIF

C  linearization of average secants for pseudo-spherical case
C  ==========================================================

C   (average secant = 1/mu-0 = constant for plane parallel)

      IF ( .NOT. DO_PLANE_PARALLEL ) THEN

C  Profile linearization

        IF( DO_PROFILE_WFS ) THEN
         DO IB = 1, NBEAMS
          DO N = 1, NLAYERS
            IF ( N .EQ. 1 ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                L_AVERAGE_SECANT(N,IB,N,Q) = 0.0D0
              ENDDO
            ELSE
              IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
                DELT  = DELTAU_VERT(N)
                LAMDA = AVERAGE_SECANT(N,IB)
                FAC   = ( CHAPMAN_FACTORS(N,N,IB) - LAMDA ) / DELT
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  L_AVERAGE_SECANT(N,IB,N,Q) = L_DELTAU_VERT(N,Q) * FAC
                ENDDO
                DO K = 1, N-1
                  FAC = ( CHAPMAN_FACTORS(N,K,IB) -
     &                    CHAPMAN_FACTORS(N-1,K,IB) ) / DELT
                  DO Q = 1, LAYER_VARY_NUMBER(K)
                    L_AVERAGE_SECANT(N,IB,K,Q) =
     &                  L_DELTAU_VERT(K,Q) * FAC
                  ENDDO
                ENDDO
              ELSE
                DO K = 1, N
                  DO Q = 1, LAYER_VARY_NUMBER(K)
                    L_AVERAGE_SECANT(N,K,IB,Q) = 0.0D0
                  ENDDO
                ENDDO
              ENDIF
            ENDIF
          ENDDO
         ENDDO
        ENDIF

C  Column linearization

        IF( DO_COLUMN_WFS ) THEN
          DO IB = 1, NBEAMS
            N = 1
            DO Q = 1, N_COLUMN_WFS
              L_AVERAGE_SECANT(N,IB,0,Q) = 0.0D0
            ENDDO
            DO N = 2, NLAYERS
              IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
                DELT  = DELTAU_VERT(N)
                LAMDA = AVERAGE_SECANT(N,IB)
                FAC   = ( CHAPMAN_FACTORS(N,N,IB) - LAMDA ) / DELT
                DO Q = 1, N_COLUMN_WFS
                  L_AVERAGE_SECANT(N,IB,0,Q) = L_DELTAU_VERT(N,Q) * FAC
                ENDDO
                DO K = 1, N-1
                  FAC = ( CHAPMAN_FACTORS(N,K,IB) -
     &                     CHAPMAN_FACTORS(N-1,K,IB) ) / DELT
                  DO Q = 1, N_COLUMN_WFS
                    L_AVERAGE_SECANT(N,IB,0,Q) =
     &              L_AVERAGE_SECANT(N,IB,0,Q) + L_DELTAU_VERT(K,Q)*FAC
                  ENDDO
                ENDDO
              ELSE
                DO Q = 1, N_COLUMN_WFS
                  L_AVERAGE_SECANT(N,IB,0,Q) = 0.0D0
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDIF

C  End pseudo-spherical clause
C  Plane parallel (safety first)

      ELSE

        DO IB = 1, NBEAMS
          DO N = 1, NLAYERS
            DO K = 0, NLAYERS
              DO Q = 1, NPARS
                L_AVERAGE_SECANT(N,IB,K,Q) = 0.0D0
              ENDDO
            ENDDO
          ENDDO
        ENDDO

      ENDIF

C  debug
c      do N = 1, nlayers
c       write(*,'(2i3,1p3e20.10)')n,1,(l_average_secant(n,1,k,1),k=1,n)
c       write(*,'(2i3,1p3e20.10)')n,2,(l_average_secant(n,1,k,2),k=1,n)
c      enddo
c      pause
C  debug
c      do N = 1, nlayers
c       write(*,'(a,i3,1p2e20.10)')'lin',n,l_initial_trans(n,1,0,1),
c     *             l_average_secant(n,1,0,1)
c      enddo

C  Linearization of Whole layer Transmittance factors
C  ==================================================

C  profile linearization
C  ---------------------

      IF ( DO_PROFILE_WFS ) THEN
       DO IB = 1, NBEAMS
        DO N = 1, NLAYERS

         WDEL  = T_DELT_MUBAR(N,IB)
         LAMDA = AVERAGE_SECANT(N,IB)
         VAR   = - DELTAU_VERT(N) * WDEL
         FAC   = - WDEL * AVERAGE_SECANT(N,IB)
 
C  Pseudo-spherical

         IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          IF ( N .EQ. 1 ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_T_DELT_MUBAR(N,IB,N,Q) = FAC * L_DELTAU_VERT(N,Q)
            ENDDO
          ELSE
            IF  ( N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                RHO = L_AVERAGE_SECANT(N,IB,N,Q)
                L_T_DELT_MUBAR(N,IB,N,Q) = L_DELTAU_VERT(N,Q) * FAC
     &                                     + VAR * RHO
              ENDDO
              DO K = 1, N-1
                DO Q = 1, LAYER_VARY_NUMBER(K)
                  RHO = L_AVERAGE_SECANT(N,IB,K,Q)
                  L_T_DELT_MUBAR(N,IB,K,Q) = VAR * RHO
                ENDDO
              ENDDO
            ELSE
              DO K = 1, N
                DO Q = 1, LAYER_VARY_NUMBER(K)
                  L_T_DELT_MUBAR(N,IB,K,Q) = 0.0D0
                ENDDO
              ENDDO
            ENDIF
          ENDIF

C  Plane-parallel

         ELSE IF ( DO_PLANE_PARALLEL ) THEN

          IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_T_DELT_MUBAR(N,IB,N,Q) = FAC * L_DELTAU_VERT(N,Q)
            ENDDO
            DO K = 1, N-1
              DO Q = 1, LAYER_VARY_NUMBER(K)
                L_T_DELT_MUBAR(N,IB,K,Q) = 0.0D0
              ENDDO
            ENDDO
          ELSE
            DO K = 1, N
              DO Q = 1, LAYER_VARY_NUMBER(K)
                L_T_DELT_MUBAR(N,IB,K,Q) = 0.0D0
              ENDDO
            ENDDO
          ENDIF

         ENDIF

C  end layer and beam loops

        ENDDO
       ENDDO
      ENDIF

C  Column linearization
C  --------------------

      IF ( DO_COLUMN_WFS ) THEN
       DO IB = 1, NBEAMS
        DO N = 1, NLAYERS

         WDEL  = T_DELT_MUBAR(N,IB)
         LAMDA = AVERAGE_SECANT(N,IB)
         VAR   = - DELTAU_VERT(N) * WDEL
         FAC   = - WDEL * AVERAGE_SECANT(N,IB)

C  Pseudo-spherical

         IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          IF ( N .EQ. 1 ) THEN
            DO Q = 1, N_COLUMN_WFS
              L_T_DELT_MUBAR(N,IB,0,Q) = FAC * L_DELTAU_VERT(N,Q)
            ENDDO
          ELSE
            IF  ( N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
              DO Q = 1, N_COLUMN_WFS
                RHO = L_AVERAGE_SECANT(N,IB,0,Q)
                L_T_DELT_MUBAR(N,IB,0,Q) = L_DELTAU_VERT(N,Q) * FAC
     &                                     + VAR * RHO
              ENDDO
            ENDIF
          ENDIF

C  Plane-parallel

         ELSE IF ( DO_PLANE_PARALLEL ) THEN

          IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
            DO Q = 1, N_COLUMN_WFS
              L_T_DELT_MUBAR(N,IB,0,Q) = FAC * L_DELTAU_VERT(N,Q)
            ENDDO
          ENDIF

         ENDIF

C  end layer and beam loops

        ENDDO
       ENDDO
      ENDIF

C  Linearization of Transmittance factors for User Streams
C  -------------------------------------------------------

      IF ( DO_PROFILE_WFS ) THEN
       DO N = 1, NLAYERS
        IF ( LAYER_VARY_FLAG(N) ) THEN
          DO UM = 1, N_USER_STREAMS
            TRANS = T_DELT_USERM(N,UM) / USER_STREAMS(UM)
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_T_DELT_USERM(N,UM,Q) = - TRANS * L_DELTAU_VERT(N,Q)
            ENDDO
          ENDDO
        ENDIF
       ENDDO
      ENDIF

      IF ( DO_COLUMN_WFS ) THEN
        DO N = 1, NLAYERS
          DO UM = 1, N_USER_STREAMS
            TRANS = T_DELT_USERM(N,UM) / USER_STREAMS(UM)
            DO Q = 1, N_COLUMN_WFS
              L_T_DELT_USERM(N,UM,Q) = - TRANS * L_DELTAU_VERT(N,Q)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

C  FInish

      RETURN
      END

C

      SUBROUTINE TWOSTREAM_L_EMULTMASTER
     I     ( DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL,
     I       DO_PROFILE_WFS, DO_COLUMN_WFS,
     I       NLAYERS, NBEAMS, N_USER_STREAMS, NPARS,
     I       LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_COLUMN_WFS,
     I       DELTAU_VERT, L_DELTAU_VERT,
     I       USER_STREAMS, T_DELT_MUBAR, T_DELT_USERM, 
     I       ITRANS_USERM, AVERAGE_SECANT, LAYER_PIS_CUTOFF,
     I       INITIAL_TRANS, L_INITIAL_TRANS, L_AVERAGE_SECANT,
     I       L_T_DELT_MUBAR,  L_T_DELT_USERM,
     I       SIGMA_M, SIGMA_P, EMULT_HOPRULE, EMULT_UP, EMULT_DN,
     O       L_EMULT_UP, L_EMULT_DN )

C  Prepare multipliers for the Beam source terms

C  Input arguments
C  ===============

C  Numbers

      INTEGER          NLAYERS, NBEAMS, N_USER_STREAMS, NPARS

C  Control

      LOGICAL          DO_UPWELLING, DO_DNWELLING

C  Linearization flags

      LOGICAL          DO_PROFILE_WFS
      LOGICAL          DO_COLUMN_WFS

C  BEAM control

      LOGICAL          DO_PLANE_PARALLEL

C  Layer optical thickness and linearization

      DOUBLE PRECISION DELTAU_VERT   ( NLAYERS )
      DOUBLE PRECISION L_DELTAU_VERT ( NLAYERS, NPARS )

C  Linearization control

      LOGICAL          LAYER_VARY_FLAG   ( NLAYERS )
      INTEGER          LAYER_VARY_NUMBER ( NLAYERS )
      INTEGER          N_COLUMN_WFS

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
     &     AVERAGE_SECANT ( NLAYERS, NBEAMS ),
     &     INITIAL_TRANS  ( NLAYERS, NBEAMS )

C  Last layer to include Particular integral solution

      INTEGER          LAYER_PIS_CUTOFF(NBEAMS)

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

C  Linearizations of Transmittance factors 

      DOUBLE PRECISION
     &     L_T_DELT_USERM ( NLAYERS, N_USER_STREAMS, NPARS )
      DOUBLE PRECISION
     &     L_T_DELT_MUBAR   ( NLAYERS, NBEAMS, 0:NLAYERS, NPARS )
      DOUBLE PRECISION
     &     L_AVERAGE_SECANT ( NLAYERS, NBEAMS, 0:NLAYERS, NPARS )
      DOUBLE PRECISION
     &     L_INITIAL_TRANS  ( NLAYERS, NBEAMS, 0:NLAYERS, NPARS )

C  Output
C  ======

C  Linearized forcing term multipliers (saved for whole atmosphere)

      DOUBLE PRECISION L_EMULT_UP
     &       (N_USER_STREAMS,NLAYERS,NBEAMS,0:NLAYERS,NPARS)
      DOUBLE PRECISION L_EMULT_DN
     &       (N_USER_STREAMS,NLAYERS,NBEAMS,0:NLAYERS,NPARS)

C  local variables
C  ---------------

      INTEGER          N, K

C  Upwelling
C  =========

      IF ( DO_UPWELLING ) THEN

C    Profiles:  loop over all varying layers K such that K </= N 

        IF ( DO_PROFILE_WFS ) THEN  
          DO N = 1, NLAYERS
            DO K = 1, N
	      CALL TWOSTREAM_L_EMULT_UP
     I          ( NLAYERS, NBEAMS, N_USER_STREAMS, NPARS,
     I            LAYER_VARY_FLAG(K), N, K, LAYER_VARY_NUMBER(K),
     I            DO_PLANE_PARALLEL, LAYER_PIS_CUTOFF,
     I            ITRANS_USERM, AVERAGE_SECANT,
     I            T_DELT_MUBAR, SIGMA_P, SIGMA_M,
     I            T_DELT_USERM, USER_STREAMS, EMULT_UP,
     I            L_INITIAL_TRANS, L_AVERAGE_SECANT,
     I            L_T_DELT_MUBAR,  L_T_DELT_USERM,
     O            L_EMULT_UP )       
            ENDDO
          ENDDO
        ENDIF

C  debug
c          DO N = 1, NLAYERS
c            write(56,'(i3,1p3e24.12)')n,(L_EMULT_UP(1,n,1,k,1),k=1,n)
c            write(56,'(i3,1p3e24.12)')n,(L_EMULT_UP(1,n,1,k,2),k=1,n)
c            write(56,'(i3,1p3e24.12)')n,(L_T_DELT_MUBAR(n,1,k,1),k=1,n)
c            write(56,'(i3,1p3e24.12)')n,(L_T_DELT_MUBAR(n,1,k,2),k=1,n)
c          enddo

C  Column linearization

        IF ( DO_COLUMN_WFS ) THEN
          DO N = 1, NLAYERS
           CALL TWOSTREAM_L_EMULT_UP
     I          ( NLAYERS, NBEAMS, N_USER_STREAMS, NPARS,
     I            .TRUE., N, 0, N_COLUMN_WFS,
     I            DO_PLANE_PARALLEL, LAYER_PIS_CUTOFF,
     I            ITRANS_USERM, AVERAGE_SECANT,
     I            T_DELT_MUBAR, SIGMA_P, SIGMA_M,
     I            T_DELT_USERM, USER_STREAMS, EMULT_UP,
     I            L_INITIAL_TRANS, L_AVERAGE_SECANT,
     I            L_T_DELT_MUBAR,  L_T_DELT_USERM,
     O            L_EMULT_UP )       
          ENDDO
        ENDIF

      ENDIF

C  Downwelling
C  ===========

      IF ( DO_DNWELLING ) THEN

C    Profiles:  loop over all varying layers K such that K </= N 

        IF ( DO_PROFILE_WFS ) THEN  
          DO N = 1, NLAYERS
            DO K = 1, N
	      CALL TWOSTREAM_L_EMULT_DN
     I          ( NLAYERS, NBEAMS, N_USER_STREAMS, NPARS,
     I            LAYER_VARY_FLAG(K), N, K, LAYER_VARY_NUMBER(K),
     I            DELTAU_VERT, L_DELTAU_VERT,
     I            DO_PLANE_PARALLEL, LAYER_PIS_CUTOFF,
     I            ITRANS_USERM, AVERAGE_SECANT, INITIAL_TRANS,
     I            T_DELT_MUBAR, SIGMA_P, SIGMA_M, EMULT_HOPRULE,
     I            T_DELT_USERM, USER_STREAMS, EMULT_DN,
     I            L_INITIAL_TRANS, L_AVERAGE_SECANT,
     I            L_T_DELT_MUBAR,  L_T_DELT_USERM,
     O            L_EMULT_DN )       
            ENDDO
          ENDDO
        ENDIF

C  debug
c          DO N = 1, NLAYERS
c            write(56,'(i3,1p3e24.12)')n,(L_EMULT_DN(1,n,1,k,1),k=1,n)
c            write(56,'(i3,1p3e24.12)')n,(L_EMULT_DN(1,n,1,k,2),k=1,n)
c          enddo

C  Column linearization

        IF ( DO_COLUMN_WFS ) THEN
          DO N = 1, NLAYERS
           CALL TWOSTREAM_L_EMULT_DN
     I          ( NLAYERS, NBEAMS, N_USER_STREAMS, NPARS,
     I            .TRUE., N, 0, N_COLUMN_WFS,
     I            DELTAU_VERT, L_DELTAU_VERT,
     I            DO_PLANE_PARALLEL, LAYER_PIS_CUTOFF,
     I            ITRANS_USERM, AVERAGE_SECANT, INITIAL_TRANS,
     I            T_DELT_MUBAR, SIGMA_P, SIGMA_M, EMULT_HOPRULE,
     I            T_DELT_USERM, USER_STREAMS, EMULT_DN,
     I            L_INITIAL_TRANS, L_AVERAGE_SECANT,
     I            L_T_DELT_MUBAR,  L_T_DELT_USERM,
     O            L_EMULT_DN )       
          ENDDO
        ENDIF

      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE TWOSTREAM_L_EMULT_UP
     I  ( NLAYERS, NBEAMS, N_USER_STREAMS, NPARS,
     I    DOVARY, N, K, K_PARAMETERS,
     I    DO_PLANE_PARALLEL, LAYER_PIS_CUTOFF,
     I    ITRANS_USERM, AVERAGE_SECANT,
     I    T_DELT_MUBAR, SIGMA_P, SIGMA_M,
     I    T_DELT_USERM, USER_STREAMS, EMULT_UP,
     I    L_INITIAL_TRANS, L_AVERAGE_SECANT,
     I    L_T_DELT_MUBAR,  L_T_DELT_USERM,
     O    L_EMULT_UP )       

C  Numbers

      INTEGER          NLAYERS, NBEAMS, N_USER_STREAMS, NPARS

C  Linearization control

      LOGICAL          DOVARY
      INTEGER          N, K, K_PARAMETERS

C  BEAM control

      LOGICAL          DO_PLANE_PARALLEL

C  Last layer to include Particular integral solution

      INTEGER          LAYER_PIS_CUTOFF(NBEAMS)

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

C  coefficient functions for user-defined angles

      DOUBLE PRECISION 
     &      SIGMA_M(NLAYERS,N_USER_STREAMS,NBEAMS),
     &      SIGMA_P(NLAYERS,N_USER_STREAMS,NBEAMS)
      
C  forcing term multipliers (saved for whole atmosphere)

      DOUBLE PRECISION EMULT_UP
     &       (N_USER_STREAMS,NLAYERS,NBEAMS)

C  Linearizations of Transmittance factors 

      DOUBLE PRECISION
     &     L_T_DELT_USERM ( NLAYERS, N_USER_STREAMS, NPARS )
      DOUBLE PRECISION
     &     L_T_DELT_MUBAR   ( NLAYERS, NBEAMS, 0:NLAYERS, NPARS )
      DOUBLE PRECISION
     &     L_AVERAGE_SECANT ( NLAYERS, NBEAMS, 0:NLAYERS, NPARS )
      DOUBLE PRECISION
     &     L_INITIAL_TRANS  ( NLAYERS, NBEAMS, 0:NLAYERS, NPARS )

C  Output
C  ======

C  Linearized forcing term multipliers (saved for whole atmosphere)

      DOUBLE PRECISION L_EMULT_UP
     &       (N_USER_STREAMS,NLAYERS,NBEAMS,0:NLAYERS,NPARS)

C  Local variables
C  ===============

      DOUBLE PRECISION SU, V1, V2, WDEL, UDEL
      INTEGER          UM, Q, IB

C  Start Beam loop

      DO IB = 1, NBEAMS

C  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( .NOT.DOVARY .OR. (N.GT.LAYER_PIS_CUTOFF(IB)) ) THEN
         DO UM = 1, N_USER_STREAMS
           DO Q = 1, K_PARAMETERS
             L_EMULT_UP(UM,N,IB,K,Q) = 0.0d0
           ENDDO
         ENDDO
         GO TO 5678
       ENDIF

C  Profile linearizations: Two cases --------
C  (a) If N = K, multiplier for due to variations in the layer N
C  (b) If N > K, multiplier due to variations in a higher layer K
C  Column linearizations: One case ----------
C  (a) If K = 0, Multiplier for bulk (column) variations
  
C  transmittance factor

       WDEL = T_DELT_MUBAR(N,IB)

C  For the pseudo-spherical case
C  -----------------------------

       IF ( .NOT. DO_PLANE_PARALLEL ) THEN

C  Case(a)

        IF ( K.EQ.N .OR. K.EQ.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            UDEL = T_DELT_USERM(N,UM)
            SU   = - ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
            DO Q = 1, K_PARAMETERS
              V1 = - L_AVERAGE_SECANT(N,IB,K,Q) / SIGMA_P(N,UM,IB)
              IF ( K.EQ.0 ) V1 = V1 + L_INITIAL_TRANS (N,IB,K,Q)
              V2 = WDEL * L_T_DELT_USERM(N,UM,Q) +
     *             UDEL * L_T_DELT_MUBAR(N,IB,K,Q)
              L_EMULT_UP(UM,N,IB,K,Q) = EMULT_UP(UM,N,IB) * V1 + SU * V2
            ENDDO
          ENDDO
        ENDIF

C  Case (b)

        IF ( N.GT.K .AND. K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            UDEL = T_DELT_USERM(N,UM)
            SU = - ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
            DO Q = 1, K_PARAMETERS
              V1 = L_INITIAL_TRANS (N,IB,K,Q) -
     &           ( L_AVERAGE_SECANT(N,IB,K,Q) / SIGMA_P(N,UM,IB) )
              V2 =  UDEL * L_T_DELT_MUBAR(N,IB,K,Q)
              L_EMULT_UP(UM,N,IB,K,Q) = EMULT_UP(UM,N,IB) * V1 + SU * V2
            ENDDO
          ENDDO
        ENDIF

C  For the plane-parallel case
C  ---------------------------

      ELSE

C  Case (a)

        IF ( K.EQ.N .OR. K.EQ.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            UDEL = T_DELT_USERM(N,UM)
            SU = - ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
            DO Q = 1, K_PARAMETERS
              V1 = 0.0d0
              IF ( K.EQ.0 ) V1 = L_INITIAL_TRANS (N,IB,K,Q)
              V2 = WDEL * L_T_DELT_USERM(N,UM,Q) +
     &             UDEL * L_T_DELT_MUBAR(N,IB,K,Q)
              L_EMULT_UP(UM,N,IB,K,Q) = EMULT_UP(UM,N,IB)*V1 + SU * V2
            ENDDO
          ENDDO
        ENDIF

C  Case (b)

        IF ( N.GT.K .AND. K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              V1 = L_INITIAL_TRANS(N,IB,K,Q)
              L_EMULT_UP(UM,N,IB,K,Q) = EMULT_UP(UM,N,IB) * V1
            ENDDO
          ENDDO
        ENDIF

C  End clause pseudo-spherical versus plane-parallel

       ENDIF

C  continuation point for next beam

 5678  CONTINUE

      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE TWOSTREAM_L_EMULT_DN
     I  ( NLAYERS, NBEAMS, N_USER_STREAMS, NPARS,
     I    DOVARY, N, K, K_PARAMETERS,
     I    DELTAU_VERT, L_DELTAU_VERT,
     I    DO_PLANE_PARALLEL, LAYER_PIS_CUTOFF,
     I    ITRANS_USERM, AVERAGE_SECANT, INITIAL_TRANS,
     I    T_DELT_MUBAR, SIGMA_P, SIGMA_M, EMULT_HOPRULE,
     I    T_DELT_USERM, USER_STREAMS, EMULT_DN,
     I    L_INITIAL_TRANS, L_AVERAGE_SECANT,
     I    L_T_DELT_MUBAR,  L_T_DELT_USERM,
     O    L_EMULT_DN )       

C  Numbers

      INTEGER          NLAYERS, NBEAMS, N_USER_STREAMS, NPARS

C  Linearization control

      LOGICAL          DOVARY
      INTEGER          N, K, K_PARAMETERS

C  BEAM control

      LOGICAL          DO_PLANE_PARALLEL

C  Last layer to include Particular integral solution

      INTEGER          LAYER_PIS_CUTOFF(NBEAMS)

C  Layer optical thickness and linearization

      DOUBLE PRECISION DELTAU_VERT   ( NLAYERS )
      DOUBLE PRECISION L_DELTAU_VERT ( NLAYERS, NPARS )

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
     &     AVERAGE_SECANT ( NLAYERS, NBEAMS ),
     &     INITIAL_TRANS  ( NLAYERS, NBEAMS )

C  coefficient functions for user-defined angles

      DOUBLE PRECISION 
     &      SIGMA_M(NLAYERS,N_USER_STREAMS,NBEAMS),
     &      SIGMA_P(NLAYERS,N_USER_STREAMS,NBEAMS)
      
C  L'Hopital's rule logical variables

      LOGICAL          EMULT_HOPRULE
     &       (NLAYERS,N_USER_STREAMS,NBEAMS)

C  forcing term multipliers (saved for whole atmosphere)

      DOUBLE PRECISION EMULT_DN
     &       (N_USER_STREAMS,NLAYERS,NBEAMS)

C  Linearizations of Transmittance factors 

      DOUBLE PRECISION
     &     L_T_DELT_USERM ( NLAYERS, N_USER_STREAMS, NPARS )
      DOUBLE PRECISION
     &     L_T_DELT_MUBAR   ( NLAYERS, NBEAMS, 0:NLAYERS, NPARS )
      DOUBLE PRECISION
     &     L_AVERAGE_SECANT ( NLAYERS, NBEAMS, 0:NLAYERS, NPARS )
      DOUBLE PRECISION
     &     L_INITIAL_TRANS  ( NLAYERS, NBEAMS, 0:NLAYERS, NPARS )

C  Output
C  ======

C  Linearized forcing term multipliers (saved for whole atmosphere)

      DOUBLE PRECISION L_EMULT_DN
     &       (N_USER_STREAMS,NLAYERS,NBEAMS,0:NLAYERS,NPARS)

C  Local variables
C  ===============

      DOUBLE PRECISION SD, V1, V2, V3
      INTEGER          UM, Q, IB

C  Start Beam loop

      DO IB = 1, NBEAMS

C  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
         DO UM = 1, N_USER_STREAMS
           DO Q = 1, K_PARAMETERS
             L_EMULT_DN(UM,N,IB,K,Q) = 0.0d0
           ENDDO
         ENDDO
         GO TO 5678
       ENDIF

C  Profile linearizations: Two cases --------
C  (a) If N = K, multiplier for due to variations in the layer N
C  (b) If N > K, multiplier due to variations in a higher layer K
C  Column linearizations: One case ----------
C  (a) If K = 0, Multiplier for bulk (column) variations
  
C  NOTE - use of L'Hopital's Rule is present in this module

C  For the pseudo-spherical case
C  -----------------------------

       IF ( .NOT. DO_PLANE_PARALLEL ) THEN

C  Case(a). Note the use of L'Hopital's Rule flag.

        IF ( K.EQ.N .OR. K.EQ.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
              V1 = 1.0d0 - DELTAU_VERT(N) / USER_STREAMS(UM)
              V2 = - 0.5d0 * DELTAU_VERT(N)
              DO Q = 1, K_PARAMETERS
                SD = V1 * L_DELTAU_VERT(N,Q) +
     &               V2 * L_AVERAGE_SECANT(N,IB,K,Q)
                V3 = 0.0d0
                IF ( K.EQ.0) V3 = V3 + L_INITIAL_TRANS (N,IB,K,Q)
                L_EMULT_DN(UM,N,IB,K,Q) = EMULT_DN(UM,N,IB) * (SD+V3)
              ENDDO
            ELSE
              SD = ITRANS_USERM(N,UM,IB) / SIGMA_M(N,UM,IB)
              DO Q = 1, K_PARAMETERS
                V1 = - L_AVERAGE_SECANT(N,IB,K,Q) / SIGMA_M(N,UM,IB)
                IF ( K.EQ.0 ) V1 = V1 + L_INITIAL_TRANS (N,IB,K,Q)
                V2 = L_T_DELT_USERM(N,UM,Q) - L_T_DELT_MUBAR(N,IB,K,Q)
                L_EMULT_DN(UM,N,IB,K,Q) = EMULT_DN(UM,N,IB)*V1 + SD*V2
              ENDDO
            ENDIF
          ENDDO
        ENDIF

C  Case (b)

        IF ( N.GT.K .AND. K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
              V2 = - 0.5d0 * DELTAU_VERT(N)
              DO Q = 1, K_PARAMETERS
                SD =        L_INITIAL_TRANS (N,IB,K,Q) +
     &                 V2 * L_AVERAGE_SECANT(N,IB,K,Q)
                L_EMULT_DN(UM,N,IB,K,Q) = EMULT_DN(UM,N,IB) * SD
              ENDDO
            ELSE
              SD = ITRANS_USERM(N,UM,IB) / SIGMA_M(N,UM,IB)
              DO Q = 1, K_PARAMETERS
                V1 =   L_INITIAL_TRANS(N,IB,K,Q) -
     &               ( L_AVERAGE_SECANT(N,IB,K,Q) / SIGMA_M(N,UM,IB) )
                V2 = - L_T_DELT_MUBAR(N,IB,K,Q)
                L_EMULT_DN(UM,N,IB,K,Q) = EMULT_DN(UM,N,IB)*V1 + SD*V2
              ENDDO
            ENDIF
          ENDDO
        ENDIF

C  For the plane-parallel case
C  ---------------------------

       ELSE

C  Case (a)

        IF ( K.EQ.N .OR. K.EQ.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
              V1 = 1.0d0 - DELTAU_VERT(N) / USER_STREAMS(UM)
              DO Q = 1, K_PARAMETERS
                SD = V1 * L_DELTAU_VERT(N,Q)
                V2 = 0.0d0
                IF ( K.EQ.0.AND.INITIAL_TRANS(N,IB).NE.0.0d0)
     &              V2 = V2 + L_INITIAL_TRANS (N,IB,K,Q)
                L_EMULT_DN(UM,N,IB,K,Q) = EMULT_DN(UM,N,IB) * (SD+V2)
              ENDDO
            ELSE
              SD = ITRANS_USERM(N,UM,IB) / SIGMA_M(N,UM,IB)
              DO Q = 1, K_PARAMETERS
                V1 = 0.0d0
                IF ( K.EQ.0 ) V1 = L_INITIAL_TRANS (N,IB,K,Q)
                V2 = L_T_DELT_USERM(N,UM,Q) - L_T_DELT_MUBAR(N,IB,K,Q)
                L_EMULT_DN(UM,N,IB,K,Q) = EMULT_DN(UM,N,IB)*V1 + SD*V2
              ENDDO
            ENDIF
          ENDDO
        ENDIF

C  Case (b)

        IF ( N.GT.K .AND. K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              V1 = L_INITIAL_TRANS(N,IB,K,Q)
              L_EMULT_DN(UM,N,IB,K,Q) = EMULT_DN(UM,N,IB) * V1
            ENDDO
          ENDDO
        ENDIF

C  End clause pseudo-spherical versus plane-parallel

       ENDIF

C  continuation point for next beam

 5678  CONTINUE

      ENDDO

C  Finish

      RETURN
      END

