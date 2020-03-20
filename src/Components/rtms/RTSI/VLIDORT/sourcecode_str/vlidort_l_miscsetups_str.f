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
C #            VLIDORT_L_MISCSETUPS (master)                    #
C #            VLIDORT_L_DELTAMSCALE_TOTAL                      #
C #            VLIDORT_L_SSALBINIT_TOTAL                        #
C #            VLIDORT_L_PREPTRANS                              #
C #                                                             #
C ###############################################################

      SUBROUTINE VLIDORT_L_MISCSETUPS

C  miscellaneous setup operations for linearized quantities

C  Deltam scaling of variational quantities

      CALL VLIDORT_L_DELTAMSCALE

C  initialise single scatter albedo variational quantities

      CALL VLIDORT_L_SSALBINIT

C  Transmittance factors variational quantities

      CALL VLIDORT_L_PREPTRANS

C  Finish

      RETURN
      END

C

      SUBROUTINE VLIDORT_L_DELTAMSCALE

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include file of VLIDORT linearized input variables

      INCLUDE '../includes/VLIDORT_L_INPUTS.VARS'

C  include files of setup variables (input)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  include files of linearized setup variables (output to this module)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'

C  local variables

      DOUBLE PRECISION BLD, DL, OF1, F, F1, UQ, EQ
      DOUBLE PRECISION FZM, T1, L_FAC1, DELS
      DOUBLE PRECISION ZMQ, ZLQ, UZQ_SCALE, ZQ_SCALE, UQ_SCALE
      INTEGER          N, Q, NM1, L, K, IB, K1, K2, O1
      LOGICAL          LOOP

C  Indexing

      INTEGER          KTYPE1(4), KTYPE2(4)
      DATA             KTYPE1 / 1, 6, 11, 16 /
      DATA             KTYPE2 / 2, 5, 12, 15 /

C  slant optical thickness values
C    (Revert to the input values of deltau_slant)
C   Commented out line is wrong !!!!!!!!!!!!!!!!!!!!!!

      DO IB = 1, NBEAMS
        DO N = 1, NLAYERS
          DO K = 1, N
c            DELS = DELTAU_SLANT(N,K,IB)/DELTAU_VERT_INPUT(K)/FAC1(K)
            DELS = DELTAU_SLANT(N,K,IB)/FAC1(K)
            IF ( LAYER_VARY_FLAG(K) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(K)
                L_DELTAU_SLANT(Q,N,K,IB) =
     &             L_DELTAU_VERT_INPUT(Q,K) * DELS
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO

C  Set Scattering matrix linearization flag
C   Only examine the (1,1) entry (the phase function coefficient)
C    Dimensioning bug, 21 March 2007. Care with using NM1

      K1 = 1
      IF ( DO_DELTAM_SCALING ) THEN
        NM1 = NMOMENTS+1
      ELSE
        NM1 = NMOMENTS
      ENDIF

      IF ( DO_ATMOS_LINEARIZATION ) THEN
        DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              LOOP = .TRUE.
              L = 0
              DO WHILE (LOOP.AND.L.LT.NM1)
                L = L + 1
                LOOP = ( DABS(L_GREEKMAT_TOTAL_INPUT(Q,L,N,K1))
     &                    .LT.1000.0*SMALLNUM )
              ENDDO
              DO_SCATMAT_VARIATION(N,Q) = .NOT.LOOP
            ENDDO
          ENDIF
        ENDDO
      ENDIF

C  DELTAM SCALING
C  ==============

C  New section added 21 December 2005 by R. Spurr
C   Linearization of the scaled delta-M inputs

      IF ( DO_DELTAM_SCALING ) THEN

        IF ( DO_ATMOS_LINEARIZATION ) THEN

          NM1 = NMOMENTS+1

C  start loop over layers

          DO N = 1, NLAYERS

C  If the layer is varying

            IF ( LAYER_VARY_FLAG(N) ) THEN

              OF1 = ( ONE - FAC1(N) ) / FAC1(N)
              F   = TRUNC_FACTOR(N)
              F1  = ONE - F

C  Start loop over varying parameters in this layer

              DO Q = 1, LAYER_VARY_NUMBER(N)

C  scale scattering matrix linearization additionally

                IF ( DO_SCATMAT_VARIATION(N,Q) ) THEN

C  set bulk property values

                  UQ  = L_OMEGA_TOTAL_INPUT(Q,N)
                  EQ  = L_DELTAU_VERT_INPUT(Q,N)
                  ZMQ = L_GREEKMAT_TOTAL_INPUT(Q,NM1,N,1)
                  FZM = F * ZMQ
                  UZQ_SCALE = ( UQ + ZMQ ) * OF1
                  ZQ_SCALE  = FZM / F1
                  L_TRUNC_FACTOR(Q,N) = FZM
                  L_OMEGA_TOTAL(Q,N)      = UQ + UZQ_SCALE - ZQ_SCALE
                  L_DELTAU_VERT(Q,N)      = EQ - UZQ_SCALE

C  do the phase function first

                  K1 = 1
                  L_GREEKMAT_TOTAL(Q,0,N,K1) = ZERO
                  FZM = F * L_GREEKMAT_TOTAL_INPUT(Q,NM1,N,1)
                  DO L = 1, NMOMENTS
                    DL  = DFLOAT(2*L+1)
                    ZLQ = L_GREEKMAT_TOTAL_INPUT(Q,L,N,K1)
                    BLD = GREEKMAT_TOTAL_INPUT(L,N,K1) / DL
                    T1  = ( BLD*ZLQ - FZM ) / ( BLD - F )
                    L_GREEKMAT_TOTAL(Q,L,N,K1) = T1 + ZQ_SCALE
                  ENDDO

C  If there is polarization.........

                  IF ( NSTOKES .GT. 1 ) THEN

C  Now do the other Type 1 elements

                    DO L = 0, NMOMENTS
                      DL  = DFLOAT(2*L+1)
                      DO K = 2, 4
                        K1 = KTYPE1(K)
                        ZLQ = L_GREEKMAT_TOTAL_INPUT(Q,L,N,K1)
                        BLD = GREEKMAT_TOTAL_INPUT(L,N,K1) / DL
                        T1  = ( BLD*ZLQ - FZM ) / ( BLD - F )
                        L_GREEKMAT_TOTAL(Q,L,N,K1) = T1 + ZQ_SCALE
                      ENDDO
                    ENDDO

C  Now do the Type 2 elements
C    Bug discovered 24 September 2007
C    Formerly, we were overwriting the Type 1 elements

                    DO L = 0, NMOMENTS
                      DL  = DFLOAT(2*L+1)
                      DO K = 1, 4
                        K2 = KTYPE2(K)
                        ZLQ = L_GREEKMAT_TOTAL_INPUT(Q,L,N,K2)
c  !!!! WRONG           L_GREEKMAT_TOTAL(Q,L,N,K1) = ZLQ + ZQ_SCALE
                        L_GREEKMAT_TOTAL(Q,L,N,K2) = ZLQ + ZQ_SCALE
                      ENDDO
                    ENDDO

C  Finish scaled Greek matrix linearization

                  ENDIF

C  No scattering matrix linearization

                ELSE

C  Bulk property linearization

                  UQ = L_OMEGA_TOTAL_INPUT(Q,N)
                  EQ = L_DELTAU_VERT_INPUT(Q,N)
                  L_TRUNC_FACTOR(Q,N) = ZERO
                  UQ_SCALE = UQ * OF1
                  L_OMEGA_TOTAL(Q,N) = UQ + UQ_SCALE
                  L_DELTAU_VERT(Q,N) = EQ - UQ_SCALE

C   Zero all linearized scattering matrix quantities now;

                  DO L = 0, NMOMENTS
                   DO O1 = 1, NSTOKES_SQ
                    L_GREEKMAT_TOTAL(Q,L,N,O1) = ZERO
                   ENDDO
                  ENDDO

                ENDIF

              ENDDO
            ENDIF
          ENDDO

C  Scale linearized slant path optical thickness values

          DO N = 1, NLAYERS
           DO K = 1, N

C  Buggy code discovered, 30 October 2007
c            IF ( LAYER_VARY_FLAG(K) ) THEN
c              DO Q = 1, LAYER_VARY_NUMBER(K)
c                L_FAC1 = L_TRUNC_FACTOR(Q,K) *   OMEGA_TOTAL_INPUT(K)
c     &                   + TRUNC_FACTOR(K)   * L_OMEGA_TOTAL_INPUT(Q,K)
c                DO IB = 1, NBEAMS
c                  DELS = DELTAU_SLANT(N,K,IB)/FAC1(K)
c                  L_DELTAU_SLANT(Q,N,K,IB) = - L_FAC1 * DELS + 
c     &                           FAC1(K) * L_DELTAU_SLANT(Q,N,K,IB)  
c                ENDDO
c              ENDDO
c            ENDIF

            IF ( LAYER_VARY_FLAG(K) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(K)
                L_FAC1 = L_TRUNC_FACTOR(Q,K) 
     &                   + TRUNC_FACTOR(K)   * L_OMEGA_TOTAL_INPUT(Q,K)
                L_FAC1 = L_FAC1 *  OMEGA_TOTAL_INPUT(K)
                DO IB = 1, NBEAMS
                  DELS = DELTAU_SLANT(N,K,IB)/FAC1(K)
                  L_DELTAU_SLANT(Q,N,K,IB) = - L_FAC1 * DELS + 
     &                           FAC1(K) * L_DELTAU_SLANT(Q,N,K,IB)  
                ENDDO
              ENDDO
            ENDIF

          ENDDO
         ENDDO

        ENDIF

C  NO DELTAM SCALING
C  =================

C  move input geophysical variables to Workspace quantities

      ELSE

        IF ( DO_ATMOS_LINEARIZATION ) THEN

          DO N = 1, NLAYERS
            IF ( LAYER_VARY_FLAG(N) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                L_TRUNC_FACTOR(Q,N) = ZERO
                L_DELTAU_VERT(Q,N)  = L_DELTAU_VERT_INPUT(Q,N)
                L_OMEGA_TOTAL(Q,N)  = L_OMEGA_TOTAL_INPUT(Q,N)
                IF ( DO_SCATMAT_VARIATION(N,Q) ) THEN
                  DO L = 0, MAXMOMENTS
                    DO O1 = 1, NSTOKES_SQ
                      L_GREEKMAT_TOTAL(Q,L,N,O1) = 
     &                         L_GREEKMAT_TOTAL_INPUT(Q,L,N,O1)
                    ENDDO
                  ENDDO
                ELSE
                  DO L = 0, MAXMOMENTS
                    DO O1 = 1, NSTOKES_SQ
                      L_GREEKMAT_TOTAL(Q,L,N,O1) = ZERO
                    ENDDO
                  ENDDO
                ENDIF
              ENDDO
            ENDIF
          ENDDO

        ENDIF

      ENDIF

C  Finish module

      RETURN
      END

C
    
      SUBROUTINE VLIDORT_L_SSALBINIT

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include file of LIDORT linearized input variables

      INCLUDE '../includes/VLIDORT_L_INPUTS.VARS'

C  include files of setup variables (input)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  include files of linearized setup variables (output to this module)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'

C  local variables
C  ---------------

      INTEGER          N, L, Q, O1, O2, OM
      DOUBLE PRECISION VL

C  phase moment-weighted OMEGA and linearizations
C  Including phase function linearization

      IF ( DO_ATMOS_LINEARIZATION ) THEN
        DO N = 1, NLAYERS
         IF ( LAYER_VARY_FLAG(N) ) THEN
          DO Q = 1, LAYER_VARY_NUMBER(N)
           DO L = 0, NMOMENTS
            DO O1 = 1, NSTOKES
             DO O2 = 1, NSTOKES
              OM = MUELLER_INDEX(O1,O2)
              VL = L_OMEGA_TOTAL   (Q,N)   +
     &             L_GREEKMAT_TOTAL(Q,L,N,OM)
              L_OMEGA_GREEK(L,N,O1,O2,Q) = OMEGA_GREEK(L,N,O1,O2) * VL
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDIF
        ENDDO
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE VLIDORT_L_PREPTRANS

C  Linearization of various transmittance factors
C  This routine is the same as the scalar version

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include file of VLIDORT linearized input variables

      INCLUDE '../includes/VLIDORT_L_INPUTS.VARS'

C  include files of setup variables (input)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  include files of linearized setup variables (output to this module)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'

C  local variables
C  ---------------

      INTEGER          N, Q, UT, UM, K, IB, I
      DOUBLE PRECISION VD, VU, TRANS, UX, TRANS_D, TRANS_U
      DOUBLE PRECISION WDEL, VAR, RHO, FAC, DELT, LAMDA, XT
      DOUBLE PRECISION L_TAU, L_TDEL, L_TD, L_TU, LDN, LUP, SUM

C  Linearization of discrete ordinate transmittances
C  Required for the solution saving option

      IF ( DO_SOLUTION_SAVING ) THEN

C  Whole layers

        DO N = 1, NLAYERS
          DO Q = 1, LAYER_VARY_NUMBER(N)
            L_TAU = L_DELTAU_VERT(Q,N) * DELTAU_VERT(N)
            DO I = 1, NSTREAMS 
              L_TDEL = - L_TAU / QUAD_STREAMS(I)
              L_T_DELT_DISORDS(I,N,Q) = T_DELT_DISORDS(I,N) * L_TDEL
            ENDDO
          ENDDO
        ENDDO

C  Partial layers

        DO UT = 1, N_PARTLAYERS
          XT  = PARTAU_VERT(UT)
          N   = PARTLAYERS_LAYERIDX(UT)
          DO Q = 1, LAYER_VARY_NUMBER(N)
            L_TD = - L_DELTAU_VERT(Q,N) * XT
            L_TU = - L_DELTAU_VERT(Q,N) * ( DELTAU_VERT(N) - XT )
            DO I = 1, NSTREAMS
              LDN = L_TD / QUAD_STREAMS(I)
              LUP = L_TU / QUAD_STREAMS(I)
              L_T_DISORDS_UTDN(I,UT,Q) = T_DISORDS_UTDN(I,UT) * LDN
              L_T_DISORDS_UTUP(I,UT,Q) = T_DISORDS_UTUP(I,UT) * LUP
            ENDDO
          ENDDO
        ENDDO

C  end solution saving option

      ENDIF

C  linearization of Initial transmittances
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
     &         - L_DELTAU_VERT(Q,K) * DELTAU_SLANT(N-1,K,IB)
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
        DO UT = 1, N_PARTLAYERS

         N   = PARTLAYERS_LAYERIDX(UT)
         VAR = - PARTAU_VERT(UT) * T_UTDN_MUBAR(UT,IB)
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

        DO UT = 1, N_PARTLAYERS
         N   = PARTLAYERS_LAYERIDX(UT)
	  DO Q = 1, LAYER_VARY_NUMBER(N)
	    L_T_UTDN_MUBAR(UT,0,IB,Q) = ZERO
	  ENDDO
	ENDDO

        DO UT = 1, N_PARTLAYERS

         N   = PARTLAYERS_LAYERIDX(UT)
         VAR = - PARTAU_VERT(UT) * T_UTDN_MUBAR(UT,IB)
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

C  Help arrays for Green's function removed (see scalar routine)

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

      DO UT = 1, N_PARTLAYERS
        N = PARTLAYERS_LAYERIDX(UT) 
        UX = PARTAU_VERT(UT)
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

C  Finish

      RETURN
      END
