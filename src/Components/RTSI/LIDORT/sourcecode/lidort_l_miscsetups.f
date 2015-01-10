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

C ###############################################################
C #                                                             #
C # Subroutines in this Module                                  #
C #                                                             #
C #            LIDORT_L_MISCSETUPS (master)                     #
C #            LIDORT_L_DELTAMSCALE_TOTAL                       #
C #            LIDORT_L_SSALBINIT_TOTAL                         #
C #            LIDORT_L_PREPTRANS                               #
C #                                                             #
C ###############################################################

C  Extension to Version 3.3.
C  ------------------------

C  RT Solutions Inc. RJD Spurr.  14 May 2007. 
C   Introduction of Bulk-property linearization.
C    Based on old LIDORT Version 2.4 and 2.5 code.

C  Final output arrays are now defined separately for the PROFILE and
C  COLUMN weighting functions

      SUBROUTINE LIDORT_L_MISCSETUPS

C  miscellaneous setup operations for linearized quantities

C  Deltam scaling of variational quantities

      CALL LIDORT_L_DELTAMSCALE

C  initialise single scatter albedo variational quantities

      CALL LIDORT_L_SSALBINIT

C  Transmittance factors variational quantities

      CALL LIDORT_L_PREPTRANS

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_L_DELTAMSCALE

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file of linearized input variables

      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'

C  include file of setup variables (input)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  include file of linearized setup variables (output to this module)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'

C  local variables
C  ---------------

      DOUBLE PRECISION BLD, DL, OF1, F, F1, UQ, EQ
      DOUBLE PRECISION FZM, T1, L_SLANT, L_FAC1, DELS
      DOUBLE PRECISION ZMQ, ZLQ, UZQ_SCALE, ZQ_SCALE, UQ_SCALE
      INTEGER          N, Q, NM1, L, K, IB
      LOGICAL          LOOP

C  slant optical thickness values

      DO IB = 1, NBEAMS
        DO N = 1, NLAYERS
          DO K = 1, N
            DELS = CHAPMAN_FACTORS(N,K,IB)
            IF ( LAYER_VARY_FLAG(K) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(K)
                L_DELTAU_SLANT_INPUT(Q,N,K,IB) =
     &             L_DELTAU_VERT_INPUT(Q,K) * DELS
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO

C  Set phase function linearization flag
C    Dimensioning bug, 21 March 2007. Care with using NM1

      IF ( DO_DELTAM_SCALING ) THEN
        NM1 = NMOMENTS+1
      ELSE
        NM1 = NMOMENTS
      ENDIF

      IF ( DO_ATMOS_LINEARIZATION.AND..NOT.DO_THERMAL_TRANSONLY ) THEN
        DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              LOOP = .TRUE.
              L = 0
              DO WHILE (LOOP.AND.L.LT.NM1)
                L = L + 1
                LOOP = ( DABS(L_PHASMOMS_TOTAL_INPUT(Q,L,N))
     &                    .LT.1000.0*SMALLNUM )
              ENDDO
              DO_PHASFUNC_VARIATION(Q,N) = .NOT.LOOP
            ENDDO
          ENDIF
        ENDDO
      ENDIF

C  DELTAM SCALING
C  ==============

      IF ( DO_DELTAM_SCALING ) THEN

        IF ( DO_ATMOS_LINEARIZATION ) THEN

          NM1 = NMOMENTS+1
          DO N = 1, NLAYERS
            IF ( LAYER_VARY_FLAG(N) ) THEN
              OF1 = ( ONE - FAC1(N) ) / FAC1(N)
              F   = TRUNC_FACTOR(N)
              F1  = ONE - F
              DO Q = 1, LAYER_VARY_NUMBER(N)

C  scale phase function linearization additionally

                IF ( DO_PHASFUNC_VARIATION(Q,N) ) THEN

                  UQ  = L_OMEGA_TOTAL_INPUT(Q,N)
                  EQ  = L_DELTAU_VERT_INPUT(Q,N)
                  ZMQ = L_PHASMOMS_TOTAL_INPUT(Q,NM1,N)
                  FZM = F * ZMQ
                  L_TRUNC_FACTOR(Q,N) = FZM
                  UZQ_SCALE = ( UQ + ZMQ ) * OF1
                  ZQ_SCALE = FZM / F1
                  L_OMEGA_TOTAL(Q,N)      = UQ + UZQ_SCALE - ZQ_SCALE
                  L_DELTAU_VERT(Q,N)      = EQ - UZQ_SCALE
                  L_PHASMOMS_TOTAL(Q,0,N) = ZERO
                  DO L = 1, NMOMENTS
                    ZLQ = L_PHASMOMS_TOTAL_INPUT(Q,L,N)
                    DL  = DBLE(2*L+1)
                    BLD = PHASMOMS_TOTAL_INPUT(L,N) / DL
                    T1  = ( BLD*ZLQ - FZM ) / ( BLD - F )
                    L_PHASMOMS_TOTAL(Q,L,N) = T1 + ZQ_SCALE
                  ENDDO

C  No phase function linearization
C   Zero all linearized phase function quantities now;

                ELSE

                  UQ = L_OMEGA_TOTAL_INPUT(Q,N)
                  EQ = L_DELTAU_VERT_INPUT(Q,N)
                  L_TRUNC_FACTOR(Q,N) = ZERO
                  UQ_SCALE = UQ * OF1
                  L_OMEGA_TOTAL(Q,N) = UQ + UQ_SCALE
                  L_DELTAU_VERT(Q,N) = EQ - UQ_SCALE
                  DO L = 0, NMOMENTS
                    L_PHASMOMS_TOTAL(Q,L,N) = ZERO
                  ENDDO

                ENDIF

C  End layer loop

              ENDDO
            ENDIF
          ENDDO

C  Scale layer path thickness values (linearized)

          IF ( DO_SOLAR_SOURCES ) THEN
           DO N = 1, NLAYERS
            DO K = 1, N
             IF ( LAYER_VARY_FLAG(K) ) THEN
               DO Q = 1, LAYER_VARY_NUMBER(K)
                 L_FAC1 = L_TRUNC_FACTOR(Q,K) *   OMEGA_TOTAL_INPUT(K)
     &                    + TRUNC_FACTOR(K)   * L_OMEGA_TOTAL_INPUT(Q,K)
                 DO IB = 1, NBEAMS
                   L_SLANT =  -L_FAC1 *   DELTAU_SLANT_INPUT(N,K,IB) + 
     &                      FAC1(K)  * L_DELTAU_SLANT_INPUT(Q,N,K,IB)
                   L_DELTAU_SLANT(Q,N,K,IB) = L_SLANT
                 ENDDO
               ENDDO
             ENDIF
           ENDDO
          ENDDO
         ENDIF

C  End linearization clause

        ENDIF

C  NO DELTAM SCALING
C  =================

C  move input geophysical variables to Workspace quantities

      ELSE

        IF ( DO_ATMOS_LINEARIZATION ) THEN

C  Optical thickness

          DO N = 1, NLAYERS
            IF ( LAYER_VARY_FLAG(N) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                L_DELTAU_VERT(Q,N)  = L_DELTAU_VERT_INPUT(Q,N)
              ENDDO
            ENDIF
          ENDDO

C  Scattering variables just copied

          IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
            DO N = 1, NLAYERS
              IF ( LAYER_VARY_FLAG(N) ) THEN
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  L_TRUNC_FACTOR(Q,N) = ZERO
                  L_OMEGA_TOTAL(Q,N)  = L_OMEGA_TOTAL_INPUT(Q,N)
                  IF ( DO_PHASFUNC_VARIATION(Q,N) ) THEN
                    DO L = 0, MAXMOMENTS
                      L_PHASMOMS_TOTAL(Q,L,N) = 
     &                           L_PHASMOMS_TOTAL_INPUT(Q,L,N)
                    ENDDO
                  ELSE
                    DO L = 0, MAXMOMENTS
                      L_PHASMOMS_TOTAL(Q,L,N) = ZERO
                    ENDDO
                  ENDIF
                ENDDO
              ENDIF
            ENDDO
          ENDIF

C  Scale layer path thickness values (linearized)

          IF ( DO_SOLAR_SOURCES ) THEN
            DO N = 1, NLAYERS
              DO K = 1, N
                IF ( LAYER_VARY_FLAG(K) ) THEN
                  DO Q = 1, LAYER_VARY_NUMBER(K)
                    DO IB = 1, NBEAMS
                      L_SLANT =  L_DELTAU_SLANT_INPUT(Q,N,K,IB)
                      L_DELTAU_SLANT(Q,N,K,IB) = L_SLANT
                    ENDDO
                  ENDDO
                ENDIF
              ENDDO
            ENDDO
          ENDIF

C  End linearization clause

        ENDIF

C  End delta-m clause

      ENDIF

C  Finish module

      RETURN
      END

C
    
      SUBROUTINE LIDORT_L_SSALBINIT

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file of linearized input variables

      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'

C  include file of setup variables (input)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  include file of linearized setup variables (output to this module)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'

C  local variables
C  ---------------

      INTEGER          N, L, Q
      DOUBLE PRECISION VAR_L

C  phase moment-weighted OMEGA and linearizations
C  Including phase function linearization

      IF ( DO_ATMOS_LINEARIZATION ) THEN
        DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              DO L = 0, NMOMENTS
                VAR_L = L_OMEGA_TOTAL   (Q,N)   +
     &                  L_PHASMOMS_TOTAL(Q,L,N)
                L_OMEGA_MOMS(Q,N,L) = OMEGA_MOMS(N,L) * VAR_L
               ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_L_PREPTRANS

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file of linearized input variables

      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'

C  include files of setup variables (input)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  include files of linearized setup variables (output to this module)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'

C  local variables
C  ---------------

      INTEGER          N, Q, UT, UM, K, IB, I
      DOUBLE PRECISION VD, VU, TRANS, UX, TRANS_D, TRANS_U
      DOUBLE PRECISION WDEL, VAR, RHO, FAC, DELT, LAMDA
      DOUBLE PRECISION CONST, L_IT, L_WDEL, L_TAU, L_TDEL
      DOUBLE PRECISION L_TD, L_TU, LDN, LUP, XT, SUM

C  Linearization of discrete ordinate transmittances
C  Only required for the solution saving option
C    (automatic for BVP telescoping)
C  Completed by R. Spurr, RTSOLUTIONS Inc., 30 August 2005

      IF ( DO_SOLUTION_SAVING ) THEN

C  whole layers

        DO N = 1, NLAYERS
          DO Q = 1, LAYER_VARY_NUMBER(N)
            L_TAU = L_DELTAU_VERT(Q,N) * DELTAU_VERT(N)
            DO I = 1, NSTREAMS 
              L_TDEL = - L_TAU / X(I)
              L_T_DELT_DISORDS(I,N,Q) = T_DELT_DISORDS(I,N) * L_TDEL
            ENDDO
          ENDDO
        ENDDO

C  Partial layers

        DO UT = 1, N_OFFGRID_USERTAUS
          XT  = OFFGRID_UTAU_VALUES(UT)
          N   = OFFGRID_UTAU_LAYERIDX(UT)
          DO Q = 1, LAYER_VARY_NUMBER(N)
            L_TD = - L_DELTAU_VERT(Q,N) * XT
            L_TU = - L_DELTAU_VERT(Q,N) * ( DELTAU_VERT(N) - XT )
            DO I = 1, NSTREAMS
              LDN = L_TD / X(I)
              LUP = L_TU / X(I)
              L_T_DISORDS_UTDN(I,UT,Q) = T_DISORDS_UTDN(I,UT) * LDN
              L_T_DISORDS_UTUP(I,UT,Q) = T_DISORDS_UTUP(I,UT) * LUP
            ENDDO
          ENDDO
        ENDDO

C  end solution saving option

      ENDIF

C  linearization of Initial transmittances
C  =======================================

C   Bug fixed, 12 August 2005 for linearization of INITIAL_TRANS
C         Use Logarithmic derivative !!!!
C         Reason: avoids exceptions if INITIAL_TRANS underflows

C  Profile linearization

      IF ( DO_PROFILE_LINEARIZATION ) THEN
       DO IB = 1, NBEAMS
        DO N = 1, NLAYERS
          IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_INITIAL_TRANS(N,N,IB,Q) = ZERO
            ENDDO
            IF ( N .GT. 1 ) THEN
              DO K = 1, N-1
                DO Q = 1, LAYER_VARY_NUMBER(K)
                  L_INITIAL_TRANS(N,K,IB,Q) =
     &           - L_DELTAU_VERT(Q,K) * DELTAU_SLANT(N-1,K,IB)
                ENDDO
              ENDDO
            ENDIF
          ELSE
            DO K = 1, N
              DO Q = 1, LAYER_VARY_NUMBER(K)
                L_INITIAL_TRANS(N,K,IB,Q) = ZERO
              ENDDO
            ENDDO
          ENDIF
        ENDDO
       ENDDO
      ENDIF

C  Column weighting functions

      IF ( DO_COLUMN_LINEARIZATION ) THEN
        DO IB = 1, NBEAMS
          N = 1
          DO Q = 1, LAYER_VARY_NUMBER(N)
            L_INITIAL_TRANS(N,0,IB,Q) = ZERO
          ENDDO
          DO N = 2, NLAYERS
            IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                SUM = ZERO  
                DO K = 1, N-1
                  SUM = SUM + L_DELTAU_VERT(Q,K)*DELTAU_SLANT(N-1,K,IB)
                ENDDO
                L_INITIAL_TRANS(N,0,IB,Q) = - SUM
              ENDDO
            ELSE
              DO Q = 1, LAYER_VARY_NUMBER(N)
                L_INITIAL_TRANS(N,0,IB,Q) = ZERO
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

        IF( DO_PROFILE_LINEARIZATION ) THEN
         DO IB = 1, NBEAMS
          DO N = 1, NLAYERS
            IF ( N .EQ. 1 ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                L_AVERAGE_SECANT(N,N,IB,Q) = ZERO
              ENDDO
            ELSE
              IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
                DELT  = DELTAU_VERT(N)
                LAMDA = AVERAGE_SECANT(N,IB)
                FAC   = ( DELTAU_SLANT(N,N,IB) / DELT ) - LAMDA
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  L_AVERAGE_SECANT(N,N,IB,Q) = L_DELTAU_VERT(Q,N) * FAC
                ENDDO
                DO K = 1, N-1
                  FAC = ( DELTAU_SLANT(N,K,IB)   -
     &                    DELTAU_SLANT(N-1,K,IB) ) / DELT
                  DO Q = 1, LAYER_VARY_NUMBER(K)
                    L_AVERAGE_SECANT(N,K,IB,Q) =
     &                  L_DELTAU_VERT(Q,K) * FAC
                  ENDDO
                ENDDO
              ELSE
                DO K = 1, N
                  DO Q = 1, LAYER_VARY_NUMBER(K)
                    L_AVERAGE_SECANT(N,K,IB,Q) = ZERO
                  ENDDO
                ENDDO
              ENDIF
            ENDIF
          ENDDO
         ENDDO
        ENDIF

C  Column linearization

        IF( DO_COLUMN_LINEARIZATION ) THEN
          DO IB = 1, NBEAMS
            N = 1
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_AVERAGE_SECANT(N,0,IB,Q) = ZERO
            ENDDO
            DO N = 2, NLAYERS
              IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
                DELT  = DELTAU_VERT(N)
                LAMDA = AVERAGE_SECANT(N,IB)
                FAC   = ( DELTAU_SLANT(N,N,IB) / DELT ) - LAMDA
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  L_AVERAGE_SECANT(N,0,IB,Q) = L_DELTAU_VERT(Q,N) * FAC
                ENDDO
                DO K = 1, N-1
                  FAC = ( DELTAU_SLANT(N,K,IB)   -
     &                    DELTAU_SLANT(N-1,K,IB) ) / DELT
                  DO Q = 1, LAYER_VARY_NUMBER(K)
                    L_AVERAGE_SECANT(N,0,IB,Q) =
     &              L_AVERAGE_SECANT(N,0,IB,Q) + L_DELTAU_VERT(Q,K)*FAC
                  ENDDO
                ENDDO
              ELSE
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  L_AVERAGE_SECANT(N,0,IB,Q) = ZERO
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDIF

C  End pseudo-spherical clause

      ENDIF

C  debug
c      do N = 1, nlayers
c       write(*,'(2i3,1p3e20.10)')n,1,(l_average_secant(n,k,1,1),k=1,n)
c       write(*,'(2i3,1p3e20.10)')n,2,(l_average_secant(n,k,1,2),k=1,n)
c      enddo
c      pause

C  Linearization of Whole layer Transmittance factors
C  ==================================================

C  profile linearization
C  ---------------------

      IF ( DO_PROFILE_LINEARIZATION ) THEN
       DO IB = 1, NBEAMS
        DO N = 1, NLAYERS

         WDEL  = T_DELT_MUBAR(N,IB)
         VAR   = - DELTAU_VERT(N) * WDEL
         LAMDA = AVERAGE_SECANT(N,IB)
         FAC   = VAR * AVERAGE_SECANT(N,IB)
 
C  Pseudo-spherical

         IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          IF ( N .EQ. 1 ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_T_DELT_MUBAR(N,N,IB,Q) = FAC * L_DELTAU_VERT(Q,N)
            ENDDO
          ELSE
            IF  ( N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                RHO = L_AVERAGE_SECANT(N,N,IB,Q)
                L_T_DELT_MUBAR(N,N,IB,Q) = L_DELTAU_VERT(Q,N) * FAC
     &                                     + VAR * RHO
              ENDDO
              DO K = 1, N-1
                DO Q = 1, LAYER_VARY_NUMBER(K)
                  RHO = L_AVERAGE_SECANT(N,K,IB,Q)
                  L_T_DELT_MUBAR(N,K,IB,Q) = VAR * RHO
                ENDDO
              ENDDO
            ELSE
              DO K = 1, N
                DO Q = 1, LAYER_VARY_NUMBER(K)
                  L_T_DELT_MUBAR(N,K,IB,Q) = ZERO
                ENDDO
              ENDDO
            ENDIF
          ENDIF

C  Plane-parallel

         ELSE IF ( DO_PLANE_PARALLEL ) THEN

          IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_T_DELT_MUBAR(N,N,IB,Q) = FAC * L_DELTAU_VERT(Q,N)
            ENDDO
            DO K = 1, N-1
              DO Q = 1, LAYER_VARY_NUMBER(K)
                L_T_DELT_MUBAR(N,K,IB,Q) = ZERO
              ENDDO
            ENDDO
          ELSE
            DO K = 1, N
              DO Q = 1, LAYER_VARY_NUMBER(K)
                L_T_DELT_MUBAR(N,K,IB,Q) = ZERO
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

      IF ( DO_COLUMN_LINEARIZATION ) THEN
       DO IB = 1, NBEAMS
        DO N = 1, NLAYERS

         WDEL  = T_DELT_MUBAR(N,IB)
         VAR   = - DELTAU_VERT(N) * WDEL
         LAMDA = AVERAGE_SECANT(N,IB)
         FAC   = VAR * AVERAGE_SECANT(N,IB)

C  Pseudo-spherical

         IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          IF ( N .EQ. 1 ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_T_DELT_MUBAR(N,0,IB,Q) = FAC * L_DELTAU_VERT(Q,N)
            ENDDO
          ELSE
            IF  ( N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                RHO = L_AVERAGE_SECANT(N,0,IB,Q)
                L_T_DELT_MUBAR(N,0,IB,Q) = L_DELTAU_VERT(Q,N) * FAC
     &                                     + VAR * RHO
              ENDDO
            ENDIF
          ENDIF

C  Plane-parallel

         ELSE IF ( DO_PLANE_PARALLEL ) THEN

          IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_T_DELT_MUBAR(N,0,IB,Q) = FAC * L_DELTAU_VERT(Q,N)
            ENDDO
          ENDIF

         ENDIF

C  end layer and beam loops

        ENDDO
       ENDDO
      ENDIF

C  Partial layer transmittance factors (for off-grid optical depths)
C  =================================================================

C  Profile linearization
C  ---------------------

      IF ( DO_PROFILE_LINEARIZATION ) THEN

       DO IB = 1, NBEAMS
        DO UT = 1, N_OFFGRID_USERTAUS

         N   = OFFGRID_UTAU_LAYERIDX(UT)
         VAR = - OFFGRID_UTAU_VALUES(UT) * T_UTDN_MUBAR(UT,IB)
         FAC = VAR * AVERAGE_SECANT(N,IB)

C  Pseudo-spherical

         IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          IF ( N .EQ. 1 ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_T_UTDN_MUBAR(UT,N,IB,Q) = FAC *  L_DELTAU_VERT(Q,N)
            ENDDO
          ELSE
            IF ( N. LE. LAYER_PIS_CUTOFF(IB) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                RHO = L_AVERAGE_SECANT(N,N,IB,Q)
                L_T_UTDN_MUBAR(UT,N,IB,Q) =  L_DELTAU_VERT(Q,N)* FAC 
     &                                       + VAR * RHO
              ENDDO
              DO K = 1, N-1
                DO Q = 1, LAYER_VARY_NUMBER(K)
                  RHO = L_AVERAGE_SECANT(N,K,IB,Q)
                  L_T_UTDN_MUBAR(UT,K,IB,Q) = VAR * RHO
                ENDDO
              ENDDO
            ELSE
              DO K = 1, N
                DO Q = 1, LAYER_VARY_NUMBER(K)
                  L_T_UTDN_MUBAR(UT,K,IB,Q) = ZERO
                ENDDO
              ENDDO
            ENDIF
          ENDIF

C  Plane-parallel

         ELSE IF ( DO_PLANE_PARALLEL ) THEN

          IF ( N. LE. LAYER_PIS_CUTOFF(IB) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_T_UTDN_MUBAR(UT,N,IB,Q) = FAC * L_DELTAU_VERT(Q,N)
            ENDDO
            DO K = 1, N-1
              DO Q = 1, LAYER_VARY_NUMBER(K)
                L_T_UTDN_MUBAR(UT,K,IB,Q) = ZERO
              ENDDO
            ENDDO
          ELSE
            DO K = 1, N
              DO Q = 1, LAYER_VARY_NUMBER(K)
                L_T_UTDN_MUBAR(UT,K,IB,Q) = ZERO
              ENDDO
            ENDDO
          ENDIF

         ENDIF

C  End optical depth and beam loops

        ENDDO
       ENDDO
      ENDIF

C  Column linearization
C  --------------------

      IF ( DO_COLUMN_LINEARIZATION ) THEN

       DO IB = 1, NBEAMS

C  zero it

	DO UT = 1, N_OFFGRID_USERTAUS
         N   = OFFGRID_UTAU_LAYERIDX(UT)
	  DO Q = 1, LAYER_VARY_NUMBER(N)
	    L_T_UTDN_MUBAR(UT,0,IB,Q) = ZERO
	  ENDDO
	ENDDO

        DO UT = 1, N_OFFGRID_USERTAUS

         N   = OFFGRID_UTAU_LAYERIDX(UT)
         VAR = - OFFGRID_UTAU_VALUES(UT) * T_UTDN_MUBAR(UT,IB)
         FAC = VAR * AVERAGE_SECANT(N,IB)

C  Pseudo-spherical

         IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          IF ( N .EQ. 1 ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_T_UTDN_MUBAR(UT,0,IB,Q) = FAC *  L_DELTAU_VERT(Q,N)
            ENDDO
          ELSE
            IF ( N. LE. LAYER_PIS_CUTOFF(IB) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                RHO = L_AVERAGE_SECANT(N,0,IB,Q)
                L_T_UTDN_MUBAR(UT,0,IB,Q) =  L_DELTAU_VERT(Q,N)* FAC 
     &                                       + VAR * RHO
              ENDDO
            ENDIF
          ENDIF

C  Plane-parallel

         ELSE IF ( DO_PLANE_PARALLEL ) THEN

          IF ( N. LE. LAYER_PIS_CUTOFF(IB) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_T_UTDN_MUBAR(UT,0,IB,Q) = FAC * L_DELTAU_VERT(Q,N)
            ENDDO
          ENDIF

         ENDIF

C  End optical depth and beam loops

        ENDDO
       ENDDO
      ENDIF

C  Linearization of Transmittance factors for User Streams
C  =======================================================

C  If no user streams, then return

      IF ( .NOT. DO_USER_STREAMS  ) RETURN

C  help arrays for Green's function linearization

      IF ( .NOT. DO_CLASSICAL_SOLUTION ) THEN
       IF ( .NOT. DO_PLANE_PARALLEL ) THEN

C  profile linearization

        IF ( DO_PROFILE_LINEARIZATION ) THEN
         DO IB = 1, NBEAMS
          DO N = 1, NLAYERS
            IF ( N .GT. 1 ) THEN
              WDEL  = T_DELT_MUBAR (N,IB)
              CONST = INITIAL_TRANS(N,IB)
              DO K = 1, N-1
                DO Q = 1, LAYER_VARY_NUMBER(K)
                  L_WDEL = L_T_DELT_MUBAR (N,K,IB,Q)
                  L_IT   = L_INITIAL_TRANS(N,K,IB,Q)
                  HELP_AQ(N,K,IB,Q) = CONST * L_IT
                  HELP_BQ(N,K,IB,Q) = - CONST * ( L_IT*WDEL + L_WDEL )
                ENDDO
              ENDDO
            ENDIF
          ENDDO
         ENDDO
        ENDIF

C  Column linearization
C    -----------------------------Do we need this????  PROBABLY NOT

        IF ( DO_COLUMN_LINEARIZATION ) THEN
         DO IB = 1, NBEAMS
          DO N = 1, NLAYERS
            IF ( N .GT. 1 ) THEN
              WDEL  = T_DELT_MUBAR (N,IB)
              CONST = INITIAL_TRANS(N,IB)
              DO Q = 1, LAYER_VARY_NUMBER(N)
                L_WDEL = L_T_DELT_MUBAR (N,0,IB,Q)
                L_IT   = L_INITIAL_TRANS(N,0,IB,Q)
                HELP_AQ(N,0,IB,Q) = CONST * L_IT
                HELP_BQ(N,0,IB,Q) = - CONST * ( L_IT*WDEL + L_WDEL )
              ENDDO
            ENDIF
           ENDDO
         ENDDO
        ENDIF

       ENDIF
      ENDIF

C  Whole Layer transmittance factors
C  ---------------------------------

      DO N = 1, NLAYERS
        IF ( LAYER_VARY_FLAG(N) ) THEN
          DO UM = 1, N_USER_STREAMS
            TRANS = T_DELT_USERM(N,UM) * 
     &                USER_SECANTS(UM) * DELTAU_VERT(N)
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_T_DELT_USERM(N,UM,Q) = - TRANS * L_DELTAU_VERT(Q,N)
            ENDDO
          ENDDO
        ENDIF
      ENDDO

C  Partial Layer transmittance factors for off-grid optical depths
C  ---------------------------------------------------------------

      DO UT = 1, N_OFFGRID_USERTAUS
        N  = OFFGRID_UTAU_LAYERIDX(UT)
        UX = OFFGRID_UTAU_VALUES(UT)
        IF ( LAYER_VARY_FLAG(N) ) THEN
          DO UM = 1, N_USER_STREAMS
            TRANS_D = T_UTDN_USERM(UT,UM) * USER_SECANTS(UM)
            TRANS_U = T_UTUP_USERM(UT,UM) * USER_SECANTS(UM)
            DO Q = 1, LAYER_VARY_NUMBER(N)
              VD = L_DELTAU_VERT(Q,N) * UX
              VU = L_DELTAU_VERT(Q,N) * ( DELTAU_VERT(N) - UX )
              L_T_UTDN_USERM(UT,UM,Q) = - TRANS_D * VD
              L_T_UTUP_USERM(UT,UM,Q) = - TRANS_U * VU
            ENDDO
          ENDDO
        ENDIF
      ENDDO

C  debug

c      k = 96
c      do n = 1, nlayers
c       write(k,'(1p9e15.7)')(l_average_secant(n,12,ib,1),ib=1,nbeams)
c       write(k,'(1p9e15.7)')(l_initial_trans(n,12,ib,1),ib=1,nbeams)
c      enddo
c      if ( do_fdtest ) pause

C  Finish

      RETURN
      END
