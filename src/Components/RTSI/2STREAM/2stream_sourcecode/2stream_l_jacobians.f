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

C ###########################################################
C #                                                         #
C #   Contains the following Master subroutines             #
C #                                                         #
C #          TWOSTREAM_UPUSER_ATMOSWF   (master)            #
C #          TWOSTREAM_DNUSER_ATMOSWF   (master)            #
C #          TWOSTREAM_STYPE1_SURFACEWF (master)            #
C #          TWOSTREAM_L_CONVERGE (master)                  #
C #                                                         #
C ###########################################################

      SUBROUTINE TWOSTREAM_UPUSER_ATMOSWF
     I      ( DO_INCLUDE_SURFACE, DO_MSMODE_LIDORT,
     I        DO_INCLUDE_DIRECTBEAM,
     I        FOURIER_COMPONENT, IPARTIC, FLUX_MULTIPLIER,
     I        NLAYERS, NBEAMS, N_USER_STREAMS, NPARS,
     I        VARIATION_INDEX, K_PARAMETERS, L_DELTAU_VERT,
     I        SURFACE_FACTOR, SURFTYPE, ALBEDO,
     I        USER_DIRECT_BEAM, CHAPMAN_FACTORS, 
     I        T_DELT_EIGEN, T_DELT_USERM, STREAM_VALUE,
     I        L_T_DELT_EIGEN, L_T_DELT_USERM,
     I        XPOS, XNEG, L_XPOS, L_XNEG, L_WLOWER,
     I        LCON, MCON, NCON, PCON,
     I        LCON_XVEC, MCON_XVEC, NCON_XVEC, PCON_XVEC,
     I        U_XPOS, U_XNEG, U_WPOS1, U_WPOS2,
     I        L_U_XPOS, L_U_XNEG, L_U_WPOS1, L_U_WPOS2,
     I        HMULT_1, HMULT_2, EMULT_UP, CUMSOURCE_UP,
     I        L_HMULT_1, L_HMULT_2, L_EMULT_UP,
     O        ATMOSWF_F_UP )

C  Subroutine input arguments
C  --------------------------

C  local control flags

      LOGICAL          DO_MSMODE_LIDORT
      LOGICAL          DO_INCLUDE_SURFACE
      LOGICAL          DO_INCLUDE_DIRECTBEAM

C  DB correction flag not required, streamlined initial version
C      LOGICAL          DO_DBCORRECTION

C  Fourier component, beam index

      INTEGER          FOURIER_COMPONENT
      INTEGER          IPARTIC

C  Numbers

      INTEGER          NLAYERS, NBEAMS, N_USER_STREAMS, NPARS

C  multipliers

      DOUBLE PRECISION FLUX_MULTIPLIER

C  Linearization control

      INTEGER          VARIATION_INDEX, K_PARAMETERS

C  Linearized optical properties

      DOUBLE PRECISION L_DELTAU_VERT(NLAYERS,NPARS)

C  Surface stuff

      INTEGER          SURFTYPE
      DOUBLE PRECISION SURFACE_FACTOR, ALBEDO
      DOUBLE PRECISION USER_DIRECT_BEAM(N_USER_STREAMS,NBEAMS)
      DOUBLE PRECISION CHAPMAN_FACTORS ( NLAYERS, NLAYERS, NBEAMS )

C  transmittance factors for +/- eigenvalues

      DOUBLE PRECISION T_DELT_EIGEN  (NLAYERS)
      DOUBLE PRECISION L_T_DELT_EIGEN(NLAYERS,NPARS)

C  Transmittance factors for user-defined stream angles

      DOUBLE PRECISION T_DELT_USERM   (NLAYERS,N_USER_STREAMS)
      DOUBLE PRECISION L_T_DELT_USERM (NLAYERS,N_USER_STREAMS,NPARS)

C  Stream value

      DOUBLE PRECISION STREAM_VALUE

C  Eigenvector solutions

      DOUBLE PRECISION XPOS(2,NLAYERS)
      DOUBLE PRECISION XNEG(2,NLAYERS)
      DOUBLE PRECISION L_XPOS(2,NLAYERS,NPARS)
      DOUBLE PRECISION L_XNEG(2,NLAYERS,NPARS)

C  Solution constants of integration

      DOUBLE PRECISION LCON(NLAYERS)
      DOUBLE PRECISION MCON(NLAYERS)
      DOUBLE PRECISION NCON(NLAYERS,NPARS)
      DOUBLE PRECISION PCON(NLAYERS,NPARS)

C  Solution constants of integration multiplied by homogeneous solutions

      DOUBLE PRECISION LCON_XVEC(2,NLAYERS)
      DOUBLE PRECISION MCON_XVEC(2,NLAYERS)
      DOUBLE PRECISION NCON_XVEC(2,NLAYERS,NPARS)
      DOUBLE PRECISION PCON_XVEC(2,NLAYERS,NPARS)

C  General beam solutions at the Upper/Lower boundary

      DOUBLE PRECISION L_WLOWER(2,NLAYERS,NPARS)

C  Eigenvectors defined at user-defined stream angles

      DOUBLE PRECISION U_XPOS(N_USER_STREAMS,NLAYERS)
      DOUBLE PRECISION U_XNEG(N_USER_STREAMS,NLAYERS)
      DOUBLE PRECISION L_U_XPOS(N_USER_STREAMS,NLAYERS,NPARS)
      DOUBLE PRECISION L_U_XNEG(N_USER_STREAMS,NLAYERS,NPARS)

C  Particular beam solutions at user-defined stream angles

      DOUBLE PRECISION U_WPOS1(N_USER_STREAMS,NLAYERS)
      DOUBLE PRECISION U_WPOS2(N_USER_STREAMS,NLAYERS)
      DOUBLE PRECISION L_U_WPOS1(N_USER_STREAMS,NLAYERS,NPARS)
      DOUBLE PRECISION L_U_WPOS2(N_USER_STREAMS,NLAYERS,0:NLAYERS,NPARS)

C  solution multipliers 

      DOUBLE PRECISION HMULT_1(N_USER_STREAMS,NLAYERS)
      DOUBLE PRECISION HMULT_2(N_USER_STREAMS,NLAYERS)
      DOUBLE PRECISION EMULT_UP(N_USER_STREAMS,NLAYERS,NBEAMS)

      DOUBLE PRECISION L_HMULT_1(N_USER_STREAMS,NLAYERS,NPARS)
      DOUBLE PRECISION L_HMULT_2(N_USER_STREAMS,NLAYERS,NPARS)
      DOUBLE PRECISION L_EMULT_UP
     &        (N_USER_STREAMS,NLAYERS,NBEAMS,0:NLAYERS,NPARS)

C  Cumulative source terms

      DOUBLE PRECISION
     U    CUMSOURCE_UP(N_USER_STREAMS,0:NLAYERS)

C  Outputs
C  -------

C  User-defined Jacobians, Fourier component

      DOUBLE PRECISION ATMOSWF_F_UP
     &         (N_USER_STREAMS,NBEAMS,0:NLAYERS,NPARS)

C  local variables
C  ---------------

C  BOA source terms

      DOUBLE PRECISION L_BOA_MSSOURCE ( N_USER_STREAMS, NPARS )
      DOUBLE PRECISION L_BOA_DBSOURCE ( N_USER_STREAMS, NPARS )

C  Reflectance integrand  a(j).x(j).I(-j)

      DOUBLE PRECISION L_IDOWN(NPARS)

C  Local layer and cumulative source terms

      DOUBLE PRECISION L_LAYERSOURCE ( N_USER_STREAMS, NPARS )
      DOUBLE PRECISION L_CUMULSOURCE ( N_USER_STREAMS, NPARS )

C  help variables

      INTEGER          UM, N, NC, Q, K, K1, IB, M
      DOUBLE PRECISION SFOR1, SFOR2, KMULT, H1, H2, H3, H4, H5, H6, LB
      DOUBLE PRECISION LCON_UXVEC, MCON_UXVEC, FAC
      DOUBLE PRECISION NCON_UXVEC, PCON_UXVEC

C  index

      K  = VARIATION_INDEX
      IB = IPARTIC
      M  = FOURIER_COMPONENT

C  Zero all Fourier components - New rule, better for safety
C    Only did this for components close to zenith (formerly)

      DO UM = 1, N_USER_STREAMS
        DO Q = 1, K_PARAMETERS
          ATMOSWF_F_UP(UM,IPARTIC,K,Q) = 0.0d0
        ENDDO
      ENDDO

C  BOA source terms
C  ----------------

C  initialise boa source terms

      DO UM = 1, N_USER_STREAMS
        DO Q = 1, K_PARAMETERS
          L_BOA_MSSOURCE(UM,Q) = 0.0d0
          L_BOA_DBSOURCE(UM,Q) = 0.0d0
        ENDDO
      ENDDO

C  Surface factor

      IF ( SURFTYPE .EQ. 1 ) THEN
        KMULT = SURFACE_FACTOR * ALBEDO * STREAM_VALUE
      ENDIF

C  PLACEHOLDER SURFTYPE 2

C  Surface Calculation, multiple scatter intensity at user defined-angles

      IF ( DO_INCLUDE_SURFACE )THEN
        N = NLAYERS
        DO Q = 1, K_PARAMETERS
          L_IDOWN(Q) = L_WLOWER(1,N,Q)
          IF ( K.EQ.N .OR. K.EQ.0 ) THEN
            H1 = NCON_XVEC(1,N,Q) * T_DELT_EIGEN(N)
            H2 = LCON_XVEC(1,N) * L_T_DELT_EIGEN(N,Q)
            H3 = LCON(N)*T_DELT_EIGEN(N)*L_XPOS(1,N,Q)
            H4 = PCON_XVEC(1,N,Q)
            H5 = MCON(N) * L_XNEG(1,N,Q)
            L_IDOWN(Q) = L_IDOWN(Q) + H1 + H2 + H3 + H4 + H5
          ELSE IF (K.LT.N.AND.K.NE.0) THEN
            H1 = NCON_XVEC(1,N,Q) * T_DELT_EIGEN(N)
            H2 = PCON_XVEC(1,N,Q) 
            L_IDOWN(Q) = L_IDOWN(Q) + H1 + H2
          ENDIF
          DO UM = 1, N_USER_STREAMS
            L_BOA_MSSOURCE(UM,Q) = KMULT * L_IDOWN(Q)
          ENDDO
        ENDDO

C  Add direct beam if flagged
C    For K > 0, profile weighting functions
C    For K = 0, Bulk (column) weighting functions

C        IF ( DO_INCLUDE_DIRECTBEAM .AND. .NOT. DO_DBCORRECTION ) THEN

        IF ( DO_INCLUDE_DIRECTBEAM ) THEN
          IF ( K .NE. 0 ) THEN
            DO UM = 1, N_USER_STREAMS
              FAC = - USER_DIRECT_BEAM(UM,IB) * CHAPMAN_FACTORS(N,K,IB) 
              DO Q = 1, K_PARAMETERS
                L_BOA_DBSOURCE(UM,Q) = L_DELTAU_VERT(K,Q) * FAC
              ENDDO
            ENDDO
          ELSE
            DO UM = 1, N_USER_STREAMS
              FAC = - USER_DIRECT_BEAM(UM,IB) 
              DO Q = 1, K_PARAMETERS
                LB = 0.0d0
                DO K1 = 1, N
                  LB = LB + L_DELTAU_VERT(K1,Q)*CHAPMAN_FACTORS(N,K1,IB)
                ENDDO
                L_BOA_DBSOURCE(UM,Q) = LB * FAC
              ENDDO
            ENDDO
          ENDIF
        ENDIF

      ENDIF

C  Initialize post-processing recursion
C  ====================================

C  Set the cumulative source term equal to BOA values

      NC = 0
      DO UM = 1, N_USER_STREAMS
        DO Q = 1, K_PARAMETERS
          L_CUMULSOURCE(UM,Q) = L_BOA_MSSOURCE(UM,Q) + 
     &                          L_BOA_DBSOURCE(UM,Q)
        ENDDO
      ENDDO

C  Recursion Loop in Source function integration
C  =============================================

      DO N = NLAYERS, 1, -1
        NC = NLAYERS + 1 - N

C  Homogeneous: Special case when N = K, or K = 0 (bulk)
C  Other cases when N not equal to K (only variation of Integ-Cons)

        IF ( N.EQ.K .OR. K.EQ.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            LCON_UXVEC = LCON(N) * U_XPOS(UM,N)
            MCON_UXVEC = MCON(N) * U_XNEG(UM,N)
            DO Q = 1, K_PARAMETERS
              NCON_UXVEC = NCON(N,Q) * U_XPOS(UM,N)
              PCON_UXVEC = PCON(N,Q) * U_XNEG(UM,N)
              H1 = LCON_UXVEC * L_HMULT_2(UM,N,Q)
              H2 = NCON_UXVEC *   HMULT_2(UM,N)
              H3 = LCON(N)*L_U_XPOS(UM,N,Q)*HMULT_2(UM,N)
              H4 = MCON_UXVEC * L_HMULT_1(UM,N,Q)
              H5 = PCON_UXVEC *   HMULT_1(UM,N)
              H6 = MCON(N)*L_U_XNEG(UM,N,Q)*HMULT_1(UM,N)
              L_LAYERSOURCE(UM,Q) = H1 + H2 + H3 + H4 + H5 + H6
            ENDDO
          ENDDO
        ELSE IF ( N.NE.K .AND. K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              NCON_UXVEC = NCON(N,Q) * U_XPOS(UM,N)
              PCON_UXVEC = PCON(N,Q) * U_XNEG(UM,N)
              H2 = NCON_UXVEC * HMULT_2(UM,N)
              H5 = PCON_UXVEC * HMULT_1(UM,N)
              L_LAYERSOURCE(UM,Q) = H2 + H5
            ENDDO
          ENDDO
        ENDIF

C  Particular for  N = K, or K = 0 (bulk)

        IF ( N.EQ.K .OR. K.EQ.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              SFOR2 = L_EMULT_UP(UM,N,IB,K,Q) *   U_WPOS2(UM,N)
     &                + EMULT_UP(UM,N,IB)     * L_U_WPOS2(UM,N,K,Q)
              L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR2 
            ENDDO
          ENDDO
          IF ( .NOT. DO_MSMODE_LIDORT ) THEN
            DO UM = 1, N_USER_STREAMS
              DO Q = 1, K_PARAMETERS
                SFOR1 = L_U_WPOS1(UM,N,Q) *   EMULT_UP(UM,N,IB) +
     &                    U_WPOS1(UM,N)   * L_EMULT_UP(UM,N,IB,K,Q)
                L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR1
              ENDDO
            ENDDO
          ENDIF

C  Particular for N > K (profile only)

        ELSE IF ( N.GT.K .AND. K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              SFOR2 = L_EMULT_UP(UM,N,IB,K,Q) *   U_WPOS2(UM,N)
     &                + EMULT_UP(UM,N,IB)     * L_U_WPOS2(UM,N,K,Q)
              L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR2 
            ENDDO
          ENDDO
          IF ( .NOT. DO_MSMODE_LIDORT ) THEN
            DO UM = 1, N_USER_STREAMS
              DO Q = 1, K_PARAMETERS
                SFOR1 = U_WPOS1(UM,N)   * L_EMULT_UP(UM,N,IB,K,Q)
                L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR1
              ENDDO
            ENDDO
          ENDIF
        ENDIF

C  Add to Linearized cumulative source sterm

        IF ( N.EQ.K .OR. K.EQ.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              L_CUMULSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q)  +
     &               T_DELT_USERM(N,UM)*L_CUMULSOURCE(UM,Q) +
     &             L_T_DELT_USERM(N,UM,Q)*CUMSOURCE_UP(UM,NC-1)
            ENDDO
          ENDDO
        ELSE IF ( N.NE.K.AND.K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              L_CUMULSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q)  +
     &             T_DELT_USERM(N,UM)*L_CUMULSOURCE(UM,Q)
            ENDDO
          ENDDO
        ENDIF

C  End layer loop

      ENDDO

C  User-defined stream output, just set to the cumulative source term

      DO UM = 1, N_USER_STREAMS
        DO Q = 1, K_PARAMETERS
          ATMOSWF_F_UP(UM,IPARTIC,K,Q) =
     &             FLUX_MULTIPLIER * L_CUMULSOURCE(UM,Q)
        ENDDO
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE TWOSTREAM_DNUSER_ATMOSWF
     I      ( DO_MSMODE_LIDORT, 
     I        FOURIER_COMPONENT, IPARTIC, FLUX_MULTIPLIER,
     I        NLAYERS, NBEAMS, N_USER_STREAMS,
     I        NPARS, VARIATION_INDEX, K_PARAMETERS,
     I        T_DELT_EIGEN, T_DELT_USERM,
     I        L_T_DELT_EIGEN, L_T_DELT_USERM,
     I        LCON, MCON, NCON, PCON,
     I        U_XPOS, U_XNEG, U_WNEG1, U_WNEG2,
     I        L_U_XPOS, L_U_XNEG, L_U_WNEG1, L_U_WNEG2,
     I        HMULT_1, HMULT_2, EMULT_DN, CUMSOURCE_DN,
     I        L_HMULT_1, L_HMULT_2, L_EMULT_DN,
     O        ATMOSWF_F_DN )

C  Subroutine input arguments
C  --------------------------

C  multipliers

      DOUBLE PRECISION FLUX_MULTIPLIER

C  Fourier component

      INTEGER          FOURIER_COMPONENT

C  Control

      LOGICAL          DO_MSMODE_LIDORT

C  beam index

      INTEGER          IPARTIC

C  Numbers

      INTEGER          NLAYERS, NBEAMS, N_USER_STREAMS, NPARS

C  Linearization control

      INTEGER          VARIATION_INDEX, K_PARAMETERS

C  transmittance factors for +/- eigenvalues

      DOUBLE PRECISION T_DELT_EIGEN(NLAYERS)
      DOUBLE PRECISION L_T_DELT_EIGEN(NLAYERS,NPARS)

C  Transmittance factors for user-defined stream angles

      DOUBLE PRECISION T_DELT_USERM   (NLAYERS,N_USER_STREAMS)
      DOUBLE PRECISION L_T_DELT_USERM (NLAYERS,N_USER_STREAMS,NPARS)

C  Solution constants of integration

      DOUBLE PRECISION LCON(NLAYERS)
      DOUBLE PRECISION MCON(NLAYERS)
      DOUBLE PRECISION NCON(NLAYERS,NPARS)
      DOUBLE PRECISION PCON(NLAYERS,NPARS)

C  Eigenvectors defined at user-defined stream angles

      DOUBLE PRECISION U_XPOS(N_USER_STREAMS,NLAYERS)
      DOUBLE PRECISION U_XNEG(N_USER_STREAMS,NLAYERS)
      DOUBLE PRECISION L_U_XPOS(N_USER_STREAMS,NLAYERS,NPARS)
      DOUBLE PRECISION L_U_XNEG(N_USER_STREAMS,NLAYERS,NPARS)

C  Particular beam solutions at user-defined stream angles

      DOUBLE PRECISION U_WNEG1(N_USER_STREAMS,NLAYERS)
      DOUBLE PRECISION U_WNEG2(N_USER_STREAMS,NLAYERS)
      DOUBLE PRECISION L_U_WNEG1(N_USER_STREAMS,NLAYERS,NPARS)
      DOUBLE PRECISION L_U_WNEG2(N_USER_STREAMS,NLAYERS,0:NLAYERS,NPARS)

C  solution multipliers 

      DOUBLE PRECISION HMULT_1(N_USER_STREAMS,NLAYERS)
      DOUBLE PRECISION HMULT_2(N_USER_STREAMS,NLAYERS)
      DOUBLE PRECISION EMULT_DN(N_USER_STREAMS,NLAYERS,NBEAMS)

      DOUBLE PRECISION L_HMULT_1(N_USER_STREAMS,NLAYERS,NPARS)
      DOUBLE PRECISION L_HMULT_2(N_USER_STREAMS,NLAYERS,NPARS)
      DOUBLE PRECISION L_EMULT_DN
     &        (N_USER_STREAMS,NLAYERS,NBEAMS,0:NLAYERS,NPARS)

C  Cumulative source terms

      DOUBLE PRECISION
     U    CUMSOURCE_DN(N_USER_STREAMS,0:NLAYERS)

C  Outputs
C  -------

C  User-defined Jacobians, Fourier component

      DOUBLE PRECISION ATMOSWF_F_DN
     &         (N_USER_STREAMS,NBEAMS,0:NLAYERS,NPARS)

C  local variables
C  ---------------

C  Local layer and cumulative source terms

      DOUBLE PRECISION L_LAYERSOURCE ( N_USER_STREAMS, NPARS )
      DOUBLE PRECISION L_CUMULSOURCE ( N_USER_STREAMS, NPARS )

C  Help variables

      INTEGER          UM, N, NC, Q, K, IB, M
      DOUBLE PRECISION SFOR1, SFOR2, H1, H2, H3, H4, H5, H6
      DOUBLE PRECISION LCON_UXVEC, MCON_UXVEC, NCON_UXVEC, PCON_UXVEC

C  index

      K  = VARIATION_INDEX
      IB = IPARTIC
      M  = FOURIER_COMPONENT

C  Zero all Fourier components - New rule, better for safety
C    Only did this for components close to zenith (formerly)

      DO UM = 1, N_USER_STREAMS
        DO Q = 1, K_PARAMETERS
          ATMOSWF_F_DN(UM,IPARTIC,K,Q) = 0.0d0
        ENDDO
      ENDDO

C  Set the cumulative source term equal to TOA values

      NC = 0
      DO UM = 1, N_USER_STREAMS
        DO Q = 1, K_PARAMETERS
          L_CUMULSOURCE(UM,Q) = 0.0d0
        ENDDO
      ENDDO

C  Cumulative source terms to layer NUT (user-defined stream angles only)
C    1. Get layer source terms
C    2. Find cumulative source term

      DO N = 1, NLAYERS
        NC = N

C  Homogeneous: Special case when N = K, or K = 0 (bulk)
C  Other cases when N not equal to K (only variation of Integ-Cons)

        IF ( N.EQ.K .OR. K.EQ.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            LCON_UXVEC = LCON(N) * U_XNEG(UM,N)
            MCON_UXVEC = MCON(N) * U_XPOS(UM,N)
            DO Q = 1, K_PARAMETERS
              NCON_UXVEC = NCON(N,Q) * U_XNEG(UM,N)
              PCON_UXVEC = PCON(N,Q) * U_XPOS(UM,N)
              H1 = LCON_UXVEC * L_HMULT_1(UM,N,Q)
              H2 = NCON_UXVEC *   HMULT_1(UM,N)
              H3 = LCON(N)*L_U_XNEG(UM,N,Q)*HMULT_1(UM,N)
              H4 = MCON_UXVEC * L_HMULT_2(UM,N,Q)
              H5 = PCON_UXVEC *   HMULT_2(UM,N)
              H6 = MCON(N)*L_U_XPOS(UM,N,Q)*HMULT_2(UM,N)
              L_LAYERSOURCE(UM,Q) = H1 + H2 + H3 + H4 + H5 + H6
            ENDDO
          ENDDO
        ELSE IF ( N.NE.K .AND. K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              NCON_UXVEC = NCON(N,Q) * U_XNEG(UM,N)
              PCON_UXVEC = PCON(N,Q) * U_XPOS(UM,N)
              H2 = NCON_UXVEC * HMULT_1(UM,N)
              H5 = PCON_UXVEC * HMULT_2(UM,N)
              L_LAYERSOURCE(UM,Q) = H2 + H5
            ENDDO
          ENDDO
        ENDIF

C  Particular for  N = K, or K = 0 (bulk)

        IF ( N.EQ.K .OR. K.EQ.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              SFOR2 = L_EMULT_DN(UM,N,IB,K,Q) *   U_WNEG2(UM,N)
     &                + EMULT_DN(UM,N,IB)     * L_U_WNEG2(UM,N,K,Q)
              L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR2 
            ENDDO
          ENDDO
          IF ( .NOT. DO_MSMODE_LIDORT ) THEN
            DO UM = 1, N_USER_STREAMS
              DO Q = 1, K_PARAMETERS
                SFOR1 = L_U_WNEG1(UM,N,Q) *   EMULT_DN(UM,N,IB) +
     &                    U_WNEG1(UM,N)   * L_EMULT_DN(UM,N,IB,K,Q)
                L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR1
              ENDDO
            ENDDO
          ENDIF

C  Particular for N > K (profile only)

        ELSE IF ( N.GT.K .AND. K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              SFOR2 = L_EMULT_DN(UM,N,IB,K,Q) *   U_WNEG2(UM,N)
     &                + EMULT_DN(UM,N,IB)     * L_U_WNEG2(UM,N,K,Q)
              L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR2 
            ENDDO
          ENDDO
          IF ( .NOT. DO_MSMODE_LIDORT ) THEN
            DO UM = 1, N_USER_STREAMS
              DO Q = 1, K_PARAMETERS
                SFOR1 = U_WNEG1(UM,N)   * L_EMULT_DN(UM,N,IB,K,Q)
                L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR1
              ENDDO
            ENDDO
          ENDIF
        ENDIF

C  Add to Linearized cumulative source sterm

        IF ( N.EQ.K .OR. K.EQ.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              L_CUMULSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q)  +
     &               T_DELT_USERM(N,UM)*L_CUMULSOURCE(UM,Q) +
     &             L_T_DELT_USERM(N,UM,Q)*CUMSOURCE_DN(UM,NC-1)
            ENDDO
          ENDDO
        ELSE IF ( N.NE.K.AND.K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              L_CUMULSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q)  +
     &             T_DELT_USERM(N,UM)*L_CUMULSOURCE(UM,Q)
            ENDDO
          ENDDO
        ENDIF

C  End layer loop

      ENDDO

C  User-defined stream output, just set to the cumulative source term

      DO UM = 1, N_USER_STREAMS
        DO Q = 1, K_PARAMETERS
          ATMOSWF_F_DN(UM,IPARTIC,K,Q) =
     &             FLUX_MULTIPLIER * L_CUMULSOURCE(UM,Q)
        ENDDO
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE TWOSTREAM_STYPE1_SURFACEWF
     I      ( DO_UPWELLING, DO_DNWELLING,
     I        DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM, 
     I        FOURIER_COMPONENT, IBEAM, FLUX_MULT,
     I        NLAYERS, NBEAMS, N_USER_STREAMS, NSPARS, ALBEDO,
     I        SURFACE_FACTOR, USER_DIRECT_BEAM, BOA_MSSOURCE,
     I        T_DELT_EIGEN, T_DELT_USERM, STREAM_VALUE,
     I        XPOS, XNEG, NCONALB, PCONALB,
     I        U_XPOS, U_XNEG, HMULT_1, HMULT_2,
     O        SURFACEWF_F_UP, SURFACEWF_F_DN )

C  Subroutine input arguments
C  --------------------------

C  Direction flags

      LOGICAL          DO_UPWELLING
      LOGICAL          DO_DNWELLING

C  local control flags

      LOGICAL          DO_INCLUDE_SURFACE
      LOGICAL          DO_INCLUDE_DIRECTBEAM

C  multipliers

      DOUBLE PRECISION FLUX_MULT

C  Fourier component

      INTEGER          FOURIER_COMPONENT

C  beam index

      INTEGER          IBEAM

C  Numbers

      INTEGER          NLAYERS, NBEAMS, N_USER_STREAMS, NSPARS

C  Surface stuff

      DOUBLE PRECISION ALBEDO
      DOUBLE PRECISION SURFACE_FACTOR
      DOUBLE PRECISION USER_DIRECT_BEAM(N_USER_STREAMS,NBEAMS)
      DOUBLE PRECISION BOA_MSSOURCE(N_USER_STREAMS)

C  Stream value

      DOUBLE PRECISION STREAM_VALUE

C  transmittance factors for +/- eigenvalues

      DOUBLE PRECISION T_DELT_EIGEN(NLAYERS)

C  Transmittance factors for user-defined stream angles

      DOUBLE PRECISION T_DELT_USERM   (NLAYERS,N_USER_STREAMS)

C  Solution constants of integration

      DOUBLE PRECISION NCONALB(NLAYERS)
      DOUBLE PRECISION PCONALB(NLAYERS)

C  Eigensolutions defined at quandrature stream angles

      DOUBLE PRECISION XPOS(2,NLAYERS)
      DOUBLE PRECISION XNEG(2,NLAYERS)

C  Eigensolutions defined at user-defined stream angles

      DOUBLE PRECISION U_XPOS(N_USER_STREAMS,NLAYERS)
      DOUBLE PRECISION U_XNEG(N_USER_STREAMS,NLAYERS)

C  solution multipliers 

      DOUBLE PRECISION HMULT_1(N_USER_STREAMS,NLAYERS)
      DOUBLE PRECISION HMULT_2(N_USER_STREAMS,NLAYERS)

C  Outputs
C  -------

C  User-defined Jacobians, Fourier component

      DOUBLE PRECISION SURFACEWF_F_UP(N_USER_STREAMS,NBEAMS,NSPARS)
      DOUBLE PRECISION SURFACEWF_F_DN(N_USER_STREAMS,NBEAMS,NSPARS)

C  local variables
C  ---------------

C  Local layer and cumulative source terms

      DOUBLE PRECISION L_LAYERSOURCE ( N_USER_STREAMS )
      DOUBLE PRECISION L_CUMULSOURCE ( N_USER_STREAMS )
      DOUBLE PRECISION L_BOA_SOURCE  ( N_USER_STREAMS )

C  Help variables

      INTEGER          UM, IB, N
      DOUBLE PRECISION H1, H2, INTEGRAND, REFLEC
      DOUBLE PRECISION NCON_HELP, PCON_HELP

C  Get the weighting functions
C  ---------------------------

      IF ( DO_UPWELLING ) THEN

C  Get the surface term (L_BOA_SOURCE)

        N  = NLAYERS
        IB = IBEAM

C  initialise Derivative of BOA source function

        DO UM = 1, N_USER_STREAMS
          L_BOA_SOURCE(UM) = 0.0d0
        ENDDO

C  Move on to continuation point if no surface

        IF ( .NOT. DO_INCLUDE_SURFACE ) GO TO 4556

C  Contribution due to derivatives of BV constants
C  -----------------------------------------------

C  First compute derivative of downward intensity Integrand at stream angles 
C        .. reflectance integrand  = a(j).x(j).I_DOWN(-j)

        H1  = NCONALB(N)*XPOS(1,N)*T_DELT_EIGEN(N)
        H2  = PCONALB(N)*XNEG(1,N) 
        INTEGRAND = ALBEDO * STREAM_VALUE  * ( H1 + H2 )
        REFLEC = SURFACE_FACTOR * INTEGRAND 
        DO UM = 1, N_USER_STREAMS
          L_BOA_SOURCE(UM) = REFLEC
        ENDDO

C   Add linearization due to kernel factor variation, diffuse term

        DO UM = 1, N_USER_STREAMS
          L_BOA_SOURCE(UM) = L_BOA_SOURCE(UM) + BOA_MSSOURCE(UM)
        ENDDO

C  Due to factor variation of direct beam (only if flagged)

        IF ( DO_INCLUDE_DIRECTBEAM ) THEN
          DO UM = 1, N_USER_STREAMS
            L_BOA_SOURCE(UM) =
     &              L_BOA_SOURCE(UM) + USER_DIRECT_BEAM(UM,IB)
          ENDDO
        ENDIF

C  Continuation point

 4556   continue

C  Upwelling Albedo weighting functions
C  ------------------------------------

C  Zero all Fourier component output

        DO UM = 1, N_USER_STREAMS
          SURFACEWF_F_UP(UM,IB,1) = 0.0d0
        ENDDO

C  Initialize recursion for user-defined stream angles only

        DO UM = 1, N_USER_STREAMS
          L_CUMULSOURCE(UM) = L_BOA_SOURCE(UM)
        ENDDO

C  Cumulative source terms to layer NUT (user-defined stream angles only)

        DO N = NLAYERS, 1, -1
          DO UM = 1, N_USER_STREAMS
            NCON_HELP =  NCONALB(N) * U_XPOS(UM,N)
            PCON_HELP =  PCONALB(N) * U_XNEG(UM,N)
            H1 = NCON_HELP * HMULT_2(UM,N)
            H2 = PCON_HELP * HMULT_1(UM,N)
            L_LAYERSOURCE(UM) = H1 + H2
            L_CUMULSOURCE(UM) = L_LAYERSOURCE(UM)  +
     &                 T_DELT_USERM(N,UM) * L_CUMULSOURCE(UM)
          ENDDO
        ENDDO

C  User-defined stream output, just set to the cumulative source term

        DO UM = 1, N_USER_STREAMS
          SURFACEWF_F_UP(UM,IB,1) = FLUX_MULT * L_CUMULSOURCE(UM)
        ENDDO

C  Finish upwelling

      ENDIF

C  Downwelling Albedo weighting functions
C  --------------------------------------

      IF ( DO_DNWELLING ) THEN

C  Zero all Fourier component output

        DO UM = 1, N_USER_STREAMS
          SURFACEWF_F_DN(UM,IB,1) = 0.0d0
        ENDDO

C  Initialize recursion for user-defined stream angles only

        DO UM = 1, N_USER_STREAMS
          L_CUMULSOURCE(UM) = 0.0d0
        ENDDO

C  Cumulative source terms to layer NUT (user-defined stream angles only)

        DO N = 1, NLAYERS
          DO UM = 1, N_USER_STREAMS
            NCON_HELP =  NCONALB(N) * U_XNEG(UM,N)
            PCON_HELP =  PCONALB(N) * U_XPOS(UM,N)
            H1 = NCON_HELP * HMULT_1(UM,N)
            H2 = PCON_HELP * HMULT_2(UM,N)
            L_LAYERSOURCE(UM) = H1 + H2
            L_CUMULSOURCE(UM) = L_LAYERSOURCE(UM)  +
     &                 T_DELT_USERM(N,UM) * L_CUMULSOURCE(UM)
          ENDDO
        ENDDO

C  User-defined stream output, just set to the cumulative source term

        DO UM = 1, N_USER_STREAMS
          SURFACEWF_F_DN(UM,IB,1) = FLUX_MULT * L_CUMULSOURCE(UM)
        ENDDO

C  Finish downwelling

      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE TWOSTREAM_L_CONVERGE
     I    ( DO_UPWELLING, DO_DNWELLING, 
     I      DO_PROFILE_WFS, DO_COLUMN_WFS, DO_SURFACE_WFS,
     I      N_GEOMS, NBEAMS, NTHREADS, NLAYERS, NPARS, NSPARS,
     I      N_USER_STREAMS, N_USER_RELAZMS, AZMFAC, UMOFF,
     I      THREAD, IBEAM, FOURIER_COMPONENT,
     I      LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_COLUMN_WFS,
     I      ATMOSWF_F_UP,    ATMOSWF_F_DN, 
     I      SURFACEWF_F_UP,  SURFACEWF_F_DN, 
     O      PROFILEWF_TOA,   PROFILEWF_BOA,
     O      COLUMNWF_TOA,    COLUMNWF_BOA,
     O      SURFACEWF_TOA,   SURFACEWF_BOA )

C  input variables
C  ---------------

C  Control

      LOGICAL          DO_UPWELLING, DO_DNWELLING

C  SS control, not required in this streamlined version
C      LOGICAL          DO_SSFULL, DO_SSCORR_OUTGOING, DO_SSCORR_NADIR

C  Linearization control

      LOGICAL          DO_PROFILE_WFS
      LOGICAL          DO_COLUMN_WFS
      LOGICAL          DO_SURFACE_WFS

C  Numbers

      INTEGER          NTHREADS, NBEAMS, NLAYERS, NPARS, NSPARS
      INTEGER          N_GEOMS, N_USER_STREAMS, N_USER_RELAZMS

C  Fourier component and thread, beam

      INTEGER          FOURIER_COMPONENT, THREAD, IBEAM

C  Local linearization control

      LOGICAL          LAYER_VARY_FLAG   ( NLAYERS )
      INTEGER          LAYER_VARY_NUMBER ( NLAYERS )
      INTEGER          N_COLUMN_WFS

C  Local  azimuth factors

      INTEGER          UMOFF ( NBEAMS, N_USER_STREAMS )
      DOUBLE PRECISION AZMFAC
     &     (N_USER_STREAMS,NBEAMS,N_USER_RELAZMS)

C  User-defined solutions

      DOUBLE PRECISION ATMOSWF_F_UP
     D   (N_USER_STREAMS,NBEAMS,0:NLAYERS,NPARS)
      DOUBLE PRECISION ATMOSWF_F_DN
     D   (N_USER_STREAMS,NBEAMS,0:NLAYERS,NPARS)

      DOUBLE PRECISION SURFACEWF_F_UP(N_USER_STREAMS,NBEAMS,NSPARS)
      DOUBLE PRECISION SURFACEWF_F_DN(N_USER_STREAMS,NBEAMS,NSPARS)

C  Single scatter solutions
C    Commented out in this streamlined version
c      DOUBLE PRECISION PROFILEWF_SS_UP(N_GEOMS,NLAYERS,NPARS)
c      DOUBLE PRECISION PROFILEWF_SS_DN(N_GEOMS,NLAYERS,NPARS)
c      DOUBLE PRECISION COLUMNWF_SS_UP(N_GEOMS,NPARS)
c      DOUBLE PRECISION COLUMNWF_SS_DN(N_GEOMS,NPARS)
c      DOUBLE PRECISION SURFACEWF_DB(N_GEOMS,NSPARS)

C  Output
C  ------

      DOUBLE PRECISION PROFILEWF_TOA(N_GEOMS,NLAYERS,NPARS,NTHREADS)
      DOUBLE PRECISION PROFILEWF_BOA(N_GEOMS,NLAYERS,NPARS,NTHREADS)
      DOUBLE PRECISION COLUMNWF_TOA(N_GEOMS,NPARS,NTHREADS)
      DOUBLE PRECISION COLUMNWF_BOA(N_GEOMS,NPARS,NTHREADS)

      DOUBLE PRECISION SURFACEWF_TOA(N_GEOMS,NSPARS,NTHREADS)
      DOUBLE PRECISION SURFACEWF_BOA(N_GEOMS,NSPARS,NTHREADS)

C  local variables
C  ---------------

      INTEGER          I, UA, V, Q, Z, N
      DOUBLE PRECISION TOLD, TAZM

C  ###################
C  Fourier 0 component
C  ###################

      IF ( FOURIER_COMPONENT.EQ.0 ) THEN

C  Profile atmospheric weighting functions
C  ---------------------------------------

       IF ( DO_PROFILE_WFS ) THEN

C  SS stuff commented out..................................

C  Copy DIFFUSE Fourier component at all output angles and optical depths
C    If no SSCORR and no DBCORR, then two options apply:
C     (a) Convergence on RADIANCE = DIFFUSE + SSTRUNCATED + DBTRUNCATED
C              (full radiance, no SS correction, no DB correction)
C     (b) Convergence on RADIANCE = DIFFUSE alone (MS only mode)
C              (SSTRUNCATED + DBTRUNCATED do not get calculated)

c        IF ( .not. DO_SSFULL ) THEN
         DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
           DO Q = 1, LAYER_VARY_NUMBER(N)
            DO I = 1, N_USER_STREAMS
             DO UA = 1, N_USER_RELAZMS
              V = UMOFF(IBEAM,I) + UA
              IF ( DO_UPWELLING ) THEN
               PROFILEWF_TOA(V,N,Q,THREAD) = ATMOSWF_F_UP(I,IBEAM,N,Q)
              ENDIF
              IF ( DO_DNWELLING ) THEN
               PROFILEWF_BOA(V,N,Q,THREAD) = ATMOSWF_F_DN(I,IBEAM,N,Q)
              ENDIF
             ENDDO
            ENDDO
           ENDDO
          ENDIF
         ENDDO
        
c        IF ( DO_SSFULL ) THEN
c         DO N = 1, NLAYERS
c          IF ( LAYER_VARY_FLAG(N) ) THEN
c           DO Q = 1, LAYER_VARY_NUMBER(N)
c            DO I = 1, N_USER_STREAMS
c             DO UA = 1, N_USER_RELAZMS
c              V = UMOFF(IBEAM,I) + UA
c              IF ( DO_UPWELLING ) THEN
c               PROFILEWF_TOA(V,N,Q,THREAD) = 0.0d0
c              ENDIF
c              IF ( DO_DNWELLING ) THEN
c               PROFILEWF_BOA(V,N,Q,THREAD) = 0.0d0
c              ENDIF
c             ENDDO
c            ENDDO
c           ENDDO
c          ENDIF
c         ENDDO
c        ENDIF

C    Add the single scatter component if flagged
c        IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
c         DO N = 1, NLAYERS
c          IF ( LAYER_VARY_FLAG(N) ) THEN
c           DO Q = 1, LAYER_VARY_NUMBER(N)
c            DO I = 1, N_USER_STREAMS
c             DO UA = 1, N_USER_RELAZMS
c              V = UMOFF(IBEAM,I) + UA
c              IF ( DO_UPWELLING ) THEN
c               PROFILEWF_TOA(V,N,Q,THREAD) =
c     &           PROFILEWF_TOA(V,N,Q,THREAD) + PROFILEWF_SS_UP(V,N,Q)
c              ENDIF
c              IF ( DO_DNWELLING ) THEN
c               PROFILEWF_BOA(V,N,Q,THREAD) =
c     &           PROFILEWF_BOA(V,N,Q,THREAD) + PROFILEWF_SS_DN(V,N,Q)
c              ENDIF
c             ENDDO
c            ENDDO
c           ENDDO
c          ENDIF
c         ENDDO
c        ENDIF

C  end profile atmospheric WF clause

       ENDIF

C  Bulk/column atmospheric weighting functions (Version 3.3)
C  ---------------------------------------------------------

C  This section newly written, 26 September 2006, installed 14 May 2007.

       IF ( DO_COLUMN_WFS ) THEN

C  Diffuse field at all output angles

c        IF ( .NOT. DO_SSFULL ) THEN
          DO Q = 1, N_COLUMN_WFS
            DO I = 1, N_USER_STREAMS
              DO UA = 1, N_USER_RELAZMS
              V = UMOFF(IBEAM,I) + UA
              IF ( DO_UPWELLING )
     &          COLUMNWF_TOA(V,Q,THREAD) = ATMOSWF_F_UP(I,IBEAM,0,Q)
              IF ( DO_DNWELLING )
     &          COLUMNWF_BOA(V,Q,THREAD) = ATMOSWF_F_DN(I,IBEAM,0,Q)
              ENDDO
            ENDDO
          ENDDO

C  SS stuff commented out...........................
c        IF ( DO_SSFULL ) THEN
c          DO Q = 1, N_COLUMN_WFS
c            DO I = 1, N_USER_STREAMS
c              DO UA = 1, N_USER_RELAZMS
c              V = UMOFF(IBEAM,I) + UA
c              IF ( DO_UPWELLING ) COLUMNWF_TOA(V,Q,THREAD) = 0.0d0
c              IF ( DO_DNWELLING ) COLUMNWF_BOA(V,Q,THREAD) = 0.0d0
c              ENDDO
c            ENDDO
c          ENDDO
c        ENDIF

C    Add the single scatter component if flagged
C  SS stuff commented out...........................
c        IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
c          DO Q = 1, N_COLUMN_WFS
c            DO I = 1, N_USER_STREAMS
c              DO UA = 1, N_USER_RELAZMS
c              V = UMOFF(IBEAM,I) + UA
c              IF ( DO_UPWELLING ) COLUMNWF_TOA(V,Q,THREAD) =
c     &          COLUMNWF_TOA(V,Q,THREAD) + COLUMNWF_SS_UP(V,Q)
c              IF ( DO_DNWELLING ) COLUMNWF_BOA(V,Q,THREAD) =
c     &          COLUMNWF_BOA(V,Q,THREAD) + COLUMNWF_SS_DN(V,Q)
c              ENDDO
c            ENDDO
c          ENDDO
c        ENDIF

C  end bulk/column atmospheric WF clause

       ENDIF

C  Reflectance weighting functions
C  -------------------------------

       IF ( DO_SURFACE_WFS ) THEN

C  DIffuse field at all output angles
C  Alternative - zero the output (single scatter case)

c        IF ( .not. DO_SSFULL ) THEN
          DO Z = 1, NSPARS
            DO I = 1, N_USER_STREAMS
              DO UA = 1, N_USER_RELAZMS
                V = UMOFF(IBEAM,I) + UA
                IF ( DO_UPWELLING )
     &            SURFACEWF_TOA(V,Z,THREAD) = SURFACEWF_F_UP(I,IBEAM,Z)
                IF ( DO_DNWELLING )
     &            SURFACEWF_BOA(V,Z,THREAD) = SURFACEWF_F_DN(I,IBEAM,Z)
              ENDDO
            ENDDO
          ENDDO

C  SS stuff commented out...........................
c        IF ( DO_SSFULL ) THEN
c          DO Z = 1, NSPARS
c            DO I = 1, N_USER_STREAMS
c              DO UA = 1, N_USER_RELAZMS
c                V = UMOFF(IBEAM,I) + UA
c                IF ( DO_UPWELLING ) SURFACEWF_TOA(V,Z,THREAD) = 0.0d0
c                IF ( DO_DNWELLING ) SURFACEWF_BOA(V,Z,THREAD) = 0.0d0
c              ENDDO
c            ENDDO
c          ENDDO
c        ENDIF

C    Add the  component if flagged (upwelling only)
C  SS stuff commented out...........................
c        IF ( DO_SSFULL.AND.DO_UPWELLING ) THEN
c          DO Z = 1, NSPARS
c            DO I = 1, N_USER_STREAMS
c              DO UA = 1, N_USER_RELAZMS
c                V = UMOFF(IBEAM,I) + UA
c                SURFACEWF_TOA(V,Z,THREAD) =
c     &             SURFACEWF_TOA(V,Z,THREAD) + SURFACEWF_DB(V,Z)
c              ENDDO
c            ENDDO
c          ENDDO
c        ENDIF

C  ENd surface WFS

       ENDIF

C  ######################
C  Fourier component = 1
C  ######################

      ELSE

C  Profile atmospheric weighting functions
C  ---------------------------------------

       IF ( DO_PROFILE_WFS ) THEN
        DO N = 1, NLAYERS
         IF ( LAYER_VARY_FLAG(N) ) THEN
          DO Q = 1, LAYER_VARY_NUMBER(N)
           DO UA = 1, N_USER_RELAZMS
            DO I = 1, N_USER_STREAMS
             V = UMOFF(IBEAM,I) + UA
             IF ( DO_UPWELLING ) THEN
              TOLD = PROFILEWF_TOA(V,N,Q,THREAD)
              TAZM = AZMFAC(I,IBEAM,UA)*ATMOSWF_F_UP(I,IBEAM,N,Q)
              PROFILEWF_TOA(V,N,Q,THREAD) = TOLD + TAZM
             ENDIF
             IF ( DO_DNWELLING ) THEN
              TOLD = PROFILEWF_BOA(V,N,Q,THREAD)
              TAZM = AZMFAC(I,IBEAM,UA)*ATMOSWF_F_DN(I,IBEAM,N,Q)
              PROFILEWF_BOA(V,N,Q,THREAD) = TOLD + TAZM
             ENDIF
            ENDDO
           ENDDO
          ENDDO
         ENDIF
        ENDDO
       ENDIF

C  Column atmospheric weighting functions
C  --------------------------------------

       IF ( DO_COLUMN_WFS ) THEN
        DO Q = 1, N_COLUMN_WFS
         DO UA = 1, N_USER_RELAZMS
          DO I = 1, N_USER_STREAMS
           V = UMOFF(IBEAM,I) + UA
           IF ( DO_UPWELLING ) THEN
            TOLD = COLUMNWF_TOA(V,Q,THREAD)
            TAZM = AZMFAC(I,IBEAM,UA)*ATMOSWF_F_UP(I,IBEAM,0,Q)
            COLUMNWF_TOA(V,Q,THREAD) = TOLD + TAZM
           ENDIF
           IF ( DO_DNWELLING ) THEN
            TOLD = COLUMNWF_BOA(V,Q,THREAD)
            TAZM = AZMFAC(I,IBEAM,UA)*ATMOSWF_F_DN(I,IBEAM,0,Q)
            COLUMNWF_BOA(V,Q,THREAD) = TOLD + TAZM
           ENDIF
          ENDDO
         ENDDO
        ENDDO
       ENDIF

C  Surface atmospheric weighting functions
C  ---------------------------------------

C  PLACEHOLDER (only for NON-LAMBERTIAN)

C  Finish Fourier
   
      ENDIF

C  Finish

      RETURN
      END

