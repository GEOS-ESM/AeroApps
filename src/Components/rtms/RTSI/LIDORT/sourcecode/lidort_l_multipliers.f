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
C #            L_HMULT_MASTER (master)                          #
C #                                                             #
C #            L_EMULT_MASTER (master)                          #
C #                L_WHOLELAYER_EMULT_UP                        #
C #                L_WHOLELAYER_EMULT_DN                        #
C #                L_PARTLAYER_EMULT_UP                         #
C #                L_PARTLAYER_EMULT_DN                         #
C #                                                             #
C #            (Green's function linearizations)                #
C #                L_QUAD_GFUNCMULT                             #
C #                L_WHOLELAYER_GMULT_UP                        #
C #                L_WHOLELAYER_GMULT_DN                        #
C #                L_PARTLAYER_GMULT_UP                         #
C #                L_PARTLAYER_GMULT_DN                         #
C #                                                             #
C ###############################################################

      SUBROUTINE L_HMULT_MASTER

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  standard files
C  --------------

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setups and multipliers (inputs to this module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  Linearization files
C  -------------------

C  include file of input control variables for the linearization

      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'

C  include files of linearized setup and solution variables (input)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'

C  linearized multiplier (output goes in this file)

      INCLUDE '../includes/LIDORT_L_MULTIPLIERS.VARS'

C  Local variables
C  ---------------

      LOGICAL          DO_RTSOL_VARY ( MAXLAYERS )
      INTEGER          NPARAMS_VARY  ( MAXLAYERS )
      INTEGER          N, UT, UTA, UM, AA, Q
      DOUBLE PRECISION UDEL, ZDEL, L_UDEL, L_ZDEL
      DOUBLE PRECISION UX_UP, L_UX_UP, UX_DN, L_UX_DN
      DOUBLE PRECISION L_T2, L_T1, H1, H2, HOM1, HOM2, SM, FA

C  Local control
C  -------------

      IF ( DO_PROFILE_LINEARIZATION  ) THEN
        DO N = 1, NLAYERS
          DO_RTSOL_VARY(N) = LAYER_VARY_FLAG(N)
          NPARAMS_VARY(N)  = LAYER_VARY_NUMBER(N)
        ENDDO
      ELSE IF ( DO_COLUMN_LINEARIZATION ) THEN
        DO N = 1, NLAYERS
          DO_RTSOL_VARY(N) = .TRUE.
          NPARAMS_VARY(N)  = N_TOTALCOLUMN_WFS
        ENDDO
      ENDIF

C  whole layer multipliers
C  -----------------------

C    Only done if layers are flagged

      DO N = 1, NLAYERS
       IF ( STERM_LAYERMASK_UP(N).OR.STERM_LAYERMASK_DN(N) ) THEN
        IF ( DO_RTSOL_VARY(N) ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            UDEL = T_DELT_USERM(N,UM)
            SM = USER_SECANTS(UM)
            DO AA = 1, NSTREAMS
              ZDEL  = T_DELT_EIGEN(AA,N)
              DO Q = 1, NPARAMS_VARY(N)
                L_ZDEL = L_T_DELT_EIGEN(AA,N,Q)
                L_UDEL = L_T_DELT_USERM(N,UM,Q)
                L_T2 = - ZDEL * L_UDEL - L_ZDEL * UDEL
                L_T1 = L_ZDEL - L_UDEL
                HOM1 =   L_KEIGEN(AA,N,Q)*HMULT_1(AA,UM,N) + SM*L_T1
                HOM2 = - L_KEIGEN(AA,N,Q)*HMULT_2(AA,UM,N) + SM*L_T2
                L_HMULT_1(AA,UM,N,Q) = ZETA_M(AA,UM,N) * HOM1
                L_HMULT_2(AA,UM,N,Q) = ZETA_P(AA,UM,N) * HOM2
              ENDDO
            ENDDO
          ENDDO
        ENDIF 
       ENDIF
      ENDDO

C  partial layer multipliers
C  -------------------------

C  start loop over partial layers

      DO UTA = 1, N_OUT_USERTAUS
        IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN
          UT = OFFGRID_UTAU_OUTINDEX(UTA)
          N  = OFFGRID_UTAU_LAYERIDX(UT)

C  upwelling

          IF ( DO_UPWELLING. AND. STERM_LAYERMASK_UP(N) ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              UX_UP = T_UTUP_USERM(UT,UM)
              SM = USER_SECANTS(UM)
              DO AA = 1, NSTREAMS
                ZDEL  = T_DELT_EIGEN(AA,N)
                DO Q = 1, NPARAMS_VARY(N)
                  FA = L_KEIGEN(AA,N,Q)
                  L_ZDEL  = L_T_DELT_EIGEN(AA,N,Q)
                  L_UX_UP = L_T_UTUP_USERM(UT,UM,Q)
                  L_T1 = - ZDEL * L_UX_UP - L_ZDEL * UX_UP
                  L_T1 = L_T_UTDN_EIGEN(AA,UT,Q) + L_T1
                  L_T2 = L_T_UTUP_EIGEN(AA,UT,Q) - L_UX_UP
                  H1 = - FA * UT_HMULT_UD(AA,UM,UT) + SM*L_T1
                  H2 =   FA * UT_HMULT_UU(AA,UM,UT) + SM*L_T2
                  L_UT_HMULT_UU(AA,UM,UT,Q) = ZETA_M(AA,UM,N) * H2
                  L_UT_HMULT_UD(AA,UM,UT,Q) = ZETA_P(AA,UM,N) * H1
                ENDDO
              ENDDO
            ENDDO
          ENDIF

C  Downwelling

          IF ( DO_DNWELLING .AND. STERM_LAYERMASK_DN(N) ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              UX_DN = T_UTDN_USERM(UT,UM)
              SM = USER_SECANTS(UM)
              DO AA = 1, NSTREAMS
                ZDEL  = T_DELT_EIGEN(AA,N)
                DO Q = 1, NPARAMS_VARY(N)
                  FA = L_KEIGEN(AA,N,Q)
                  L_ZDEL  = L_T_DELT_EIGEN(AA,N,Q)
                  L_UX_DN = L_T_UTDN_USERM(UT,UM,Q)
                  L_T2 = - ZDEL * L_UX_DN - L_ZDEL * UX_DN
                  L_T2 = L_T_UTUP_EIGEN(AA,UT,Q) + L_T2
                  L_T1 = L_T_UTDN_EIGEN(AA,UT,Q) - L_UX_DN
                  H1 =   FA * UT_HMULT_DD(AA,UM,UT) + SM*L_T1
                  H2 = - FA * UT_HMULT_DU(AA,UM,UT) + SM*L_T2
                  L_UT_HMULT_DU(AA,UM,UT,Q) = ZETA_P(AA,UM,N) * H2
                  L_UT_HMULT_DD(AA,UM,UT,Q) = ZETA_M(AA,UM,N) * H1
                ENDDO
              ENDDO
            ENDDO
          ENDIF

C  End loop over partial layers

        ENDIF
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE L_EMULT_MASTER

C  Linearized multipliers for the Beam source terms

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  standard files
C  --------------

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setups and multipliers (inputs to this module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  Linearization files
C  -------------------

C  include file of input control variables for the linearization

      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'

C  include files of linearized setup and solution variables (input)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'

C  linearized multiplier (output goes in this file)

      INCLUDE '../includes/LIDORT_L_MULTIPLIERS.VARS'

C  Local variables
C  ---------------

      INTEGER          N, UT, UTA, K, K_PARAMETERS

C  Upwelling
C  =========

      IF ( DO_UPWELLING ) THEN

C  Whole layer upwelling
C  ---------------------

C  Loop over all  model  layers N
C    Profiles:  loop over all varying layers K such that K </= N 
C    Columns :  K = 0

       DO N = 1, NLAYERS
        IF ( STERM_LAYERMASK_UP(N) ) THEN
         IF ( DO_PROFILE_LINEARIZATION ) THEN  
          DO K = 1, NLAYERS
	   IF ( N.GE.K ) THEN
            IF ( LAYER_VARY_FLAG(K) ) THEN
             K_PARAMETERS = LAYER_VARY_NUMBER(K)
	     CALL L_WHOLELAYER_EMULT_UP ( N, K, K_PARAMETERS )
            ENDIF
           ENDIF
          ENDDO
         ELSE IF ( DO_COLUMN_LINEARIZATION ) THEN
          K = 0
          CALL L_WHOLELAYER_EMULT_UP ( N, K, N_TOTALCOLUMN_WFS )
         ENDIF
        ENDIF
       ENDDO

C  Partial layer upwelling
C  -----------------------

C  Start loop over all partial output UT occuring in layers N
C    Profiles:  loop over all varying layers K such that K </= N 
C    Columns :  K = 0

       DO UTA = 1, N_OUT_USERTAUS
        IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN
         UT = OFFGRID_UTAU_OUTINDEX(UTA)
         N  = OFFGRID_UTAU_LAYERIDX(UT)
         IF ( STERM_LAYERMASK_UP(N) ) THEN
          IF ( DO_PROFILE_LINEARIZATION ) THEN  
           DO K = 1, NLAYERS
	    IF ( N.GE.K ) THEN
             IF ( LAYER_VARY_FLAG(K) ) THEN
              K_PARAMETERS = LAYER_VARY_NUMBER(K)
	      CALL L_PARTLAYER_EMULT_UP ( N, UT, K, K_PARAMETERS )
             ENDIF
            ENDIF
           ENDDO
          ELSE IF ( DO_COLUMN_LINEARIZATION ) THEN
           K = 0
           CALL L_PARTLAYER_EMULT_UP ( N, UT, K, N_TOTALCOLUMN_WFS )
          ENDIF
         ENDIF
        ENDIF
       ENDDO

C  end upwelling

      ENDIF

C  Downwelling
C  ===========

      IF ( DO_DNWELLING ) THEN

C  Whole layer downwelling
C  -----------------------

C  Start loop over all  model  layers N
C  Start loop over all varying layers K such that K </= N 

       DO N = 1, NLAYERS
        IF ( STERM_LAYERMASK_DN(N) ) THEN
         IF ( DO_PROFILE_LINEARIZATION ) THEN  
          DO K = 1, NLAYERS
	   IF ( N.GE.K ) THEN
            IF ( LAYER_VARY_FLAG(K) ) THEN
             K_PARAMETERS = LAYER_VARY_NUMBER(K)
	     CALL L_WHOLELAYER_EMULT_DN ( N, K, K_PARAMETERS )
            ENDIF
           ENDIF
          ENDDO
         ELSE IF  ( DO_COLUMN_LINEARIZATION ) THEN
          K = 0
          CALL L_WHOLELAYER_EMULT_DN ( N, K, N_TOTALCOLUMN_WFS )
         ENDIF 
        ENDIF
       ENDDO

C  Partial layer downwelling
C  -------------------------

C  Start loop over all partial output UT occuring in layers N
C  Start loop over all varying layers K such that K </= N 

       DO UTA = 1, N_OUT_USERTAUS
        IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN
         UT = OFFGRID_UTAU_OUTINDEX(UTA)
         N  = OFFGRID_UTAU_LAYERIDX(UT)
         IF ( STERM_LAYERMASK_DN(N) ) THEN
          IF ( DO_PROFILE_LINEARIZATION ) THEN  
           DO K = 1, NLAYERS
            IF ( N.GE.K ) THEN
             IF ( LAYER_VARY_FLAG(K) ) THEN
              K_PARAMETERS = LAYER_VARY_NUMBER(K)
	      CALL L_PARTLAYER_EMULT_DN ( N, UT, K, K_PARAMETERS )
             ENDIF
            ENDIF
           ENDDO
          ELSE IF ( DO_COLUMN_LINEARIZATION ) THEN
           K = 0
           CALL L_PARTLAYER_EMULT_DN ( N, UT, K, N_TOTALCOLUMN_WFS )
          ENDIF 
         ENDIF
        ENDIF
       ENDDO

C  end downwelling

      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE L_WHOLELAYER_EMULT_UP ( N, K, K_PARAMETERS )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setups and multiplier variables (input)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  include file of linearized setup variables (input)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'

C  linearized multiplier output goes in this file

      INCLUDE '../includes/LIDORT_L_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index, varying layer index and number of parameters

      INTEGER          N, K, K_PARAMETERS

C  local variables
C  ---------------

      DOUBLE PRECISION SU, V1, V2, WDEL, UDEL
      INTEGER          UM, Q, IB

C  Start Beam loop
C  ===============

      DO IB = 1, NBEAMS

C  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
         DO UM = 1, N_USER_STREAMS
           DO Q = 1, K_PARAMETERS
             L_EMULT_UP(UM,N,K,IB,Q) = ZERO
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
            SU = - ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
            DO Q = 1, K_PARAMETERS
              V1 = -L_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_P(N,UM,IB)
              IF ( K.EQ.0 ) V1 = V1 + L_INITIAL_TRANS (N,K,IB,Q)
              V2 = WDEL * L_T_DELT_USERM(N,UM,Q) +
     *             UDEL * L_T_DELT_MUBAR(N,K,IB,Q)
              L_EMULT_UP(UM,N,K,IB,Q) = EMULT_UP(UM,N,IB) * V1 + SU * V2
            ENDDO
          ENDDO
        ENDIF

C  Case (b)

        IF ( N.GT.K .AND. K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            UDEL = T_DELT_USERM(N,UM)
            SU = - ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
            DO Q = 1, K_PARAMETERS
              V1 = L_INITIAL_TRANS (N,K,IB,Q) -
     &           ( L_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_P(N,UM,IB) )
              V2 =  UDEL * L_T_DELT_MUBAR(N,K,IB,Q)
              L_EMULT_UP(UM,N,K,IB,Q) = EMULT_UP(UM,N,IB) * V1 + SU * V2
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
              V1 = ZERO
              IF ( K.EQ.0 ) V1 = L_INITIAL_TRANS (N,K,IB,Q)
              V2 = WDEL * L_T_DELT_USERM(N,UM,Q) +
     &             UDEL * L_T_DELT_MUBAR(N,K,IB,Q)
c              L_EMULT_UP(UM,N,K,IB,Q) =  SU * V2
              L_EMULT_UP(UM,N,K,IB,Q) = EMULT_UP(UM,N,IB)*V1 + SU * V2
            ENDDO
          ENDDO
        ENDIF

C  Case (b)

        IF ( N.GT.K .AND. K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              V1 = L_INITIAL_TRANS(N,K,IB,Q)
              L_EMULT_UP(UM,N,K,IB,Q) = EMULT_UP(UM,N,IB) * V1
            ENDDO
          ENDDO
        ENDIF

C  End clause pseudo-spherical versus plaen-parallel

       ENDIF

C  continuation point for next beam

 5678  CONTINUE

      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE L_WHOLELAYER_EMULT_DN ( N, K, K_PARAMETERS )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setups and multiplier variables (input)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  include file of linearized setup variables (input)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'

C  linearized multiplier output goes in this file

      INCLUDE '../includes/LIDORT_L_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index, varying layer index and number of parameters

      INTEGER          N, K, K_PARAMETERS

C  local variables
C  ---------------

      DOUBLE PRECISION SD, V1, V2, V3
      INTEGER          UM, Q, IB

C  Start Beam loop
C  ===============

      DO IB = 1, NBEAMS

C  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
         DO UM = 1, N_USER_STREAMS
           DO Q = 1, K_PARAMETERS
             L_EMULT_DN(UM,N,K,IB,Q) = ZERO
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
              V1 = ONE - DELTAU_VERT(N) * USER_SECANTS(UM)
              V2 = - HALF * DELTAU_VERT(N)
              DO Q = 1, K_PARAMETERS
                SD = V1 * L_DELTAU_VERT(Q,N) +
     &               V2 * L_AVERAGE_SECANT(N,K,IB,Q)
                V3 = ZERO
                IF ( K.EQ.0) V3 = V3 + L_INITIAL_TRANS (N,K,IB,Q)
                L_EMULT_DN(UM,N,K,IB,Q) = EMULT_DN(UM,N,IB) * (SD+V3)
              ENDDO
            ELSE
              SD = ITRANS_USERM(N,UM,IB) / SIGMA_M(N,UM,IB)
              DO Q = 1, K_PARAMETERS
                V1 = - L_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_M(N,UM,IB)
                IF ( K.EQ.0 ) V1 = V1 + L_INITIAL_TRANS (N,K,IB,Q)
                V2 = L_T_DELT_USERM(N,UM,Q) - L_T_DELT_MUBAR(N,K,IB,Q)
                L_EMULT_DN(UM,N,K,IB,Q) = EMULT_DN(UM,N,IB)*V1 + SD*V2
              ENDDO
            ENDIF
          ENDDO
        ENDIF

C  Case (b)

        IF ( N.GT.K .AND. K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
              V2 = - HALF * DELTAU_VERT(N)
              DO Q = 1, K_PARAMETERS
                SD =        L_INITIAL_TRANS (N,K,IB,Q) +
     &                 V2 * L_AVERAGE_SECANT(N,K,IB,Q)
                L_EMULT_DN(UM,N,K,IB,Q) = EMULT_DN(UM,N,IB) * SD
              ENDDO
            ELSE
              SD = ITRANS_USERM(N,UM,IB) / SIGMA_M(N,UM,IB)
              DO Q = 1, K_PARAMETERS
                V1 =   L_INITIAL_TRANS(N,K,IB,Q) -
     &               ( L_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_M(N,UM,IB) )
                V2 = - L_T_DELT_MUBAR(N,K,IB,Q)
                L_EMULT_DN(UM,N,K,IB,Q) = EMULT_DN(UM,N,IB)*V1 + SD*V2
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
              V1 = ONE - DELTAU_VERT(N) * USER_SECANTS(UM)
              DO Q = 1, K_PARAMETERS
                SD = V1 * L_DELTAU_VERT(Q,N)
                V2 = ZERO
                IF ( K.EQ.0.AND.INITIAL_TRANS(N,IB).NE.ZERO)
     &              V2 = V2 + L_INITIAL_TRANS (N,K,IB,Q)
                L_EMULT_DN(UM,N,K,IB,Q) = EMULT_DN(UM,N,IB) * (SD+V2)
              ENDDO
            ELSE
              SD = ITRANS_USERM(N,UM,IB) / SIGMA_M(N,UM,IB)
              DO Q = 1, K_PARAMETERS
                V1 = ZERO
                IF ( K.EQ.0 ) V1 = L_INITIAL_TRANS (N,K,IB,Q)
                V2 = L_T_DELT_USERM(N,UM,Q) - L_T_DELT_MUBAR(N,K,IB,Q)
                L_EMULT_DN(UM,N,K,IB,Q) = EMULT_DN(UM,N,IB)*V1 + SD*V2
c               L_EMULT_DN(UM,N,K,IB,Q) =  SD * V2
              ENDDO
            ENDIF
          ENDDO
        ENDIF

C  Case (b)

        IF ( N.GT.K .AND. K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              V1 = L_INITIAL_TRANS(N,K,IB,Q)
              L_EMULT_DN(UM,N,K,IB,Q) = EMULT_DN(UM,N,IB) * V1
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

      SUBROUTINE L_PARTLAYER_EMULT_UP ( N, UT, K, K_PARAMETERS )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setups and multiplier variables (input)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  include file of linearized setup variables (input)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'

C  linearized multiplier output goes in this file

      INCLUDE '../includes/LIDORT_L_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Given offgrid indices, varying layer index and number of parameters

      INTEGER          N, UT, K, K_PARAMETERS

C  local variables
C  ---------------

      DOUBLE PRECISION SU, V1, V2, WDEL, UX_UP
      INTEGER          UM, Q, IB

C  Start Beam loop
C  ===============

      DO IB = 1, NBEAMS

C  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
         DO UM = 1, N_USER_STREAMS
           DO Q = 1, K_PARAMETERS
             L_UT_EMULT_UP(UM,UT,K,IB,Q) = ZERO
           ENDDO
         ENDDO
         GO TO 5678
       ENDIF

C  transmittance factor

       WDEL = T_DELT_MUBAR(N,IB)

C  Profile linearizations: Two cases --------
C  (a) If N = K, multiplier for due to variations in the layer N
C  (b) If N > K, multiplier due to variations in a higher layer K
C  Column linearizations: One case ----------
C  (a) If K = 0, Multiplier for bulk (column) variations
  
C  For the pseudo-spherical case
C  -----------------------------

       IF ( .NOT. DO_PLANE_PARALLEL ) THEN

C  Case(a)

        IF ( K.EQ.N .OR. K.EQ.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            UX_UP = T_UTUP_USERM(UT,UM)
            SU = ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
            DO Q = 1, K_PARAMETERS
              V1 = -L_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_P(N,UM,IB)
              IF ( K.EQ.0 ) V1 = V1 + L_INITIAL_TRANS (N,K,IB,Q)
              V2 =           L_T_UTDN_MUBAR(UT,K,IB,Q) -
     *               UX_UP * L_T_DELT_MUBAR(N, K,IB,Q) -
     *               WDEL  * L_T_UTUP_USERM(UT,UM,Q)
              L_UT_EMULT_UP(UM,UT,K,IB,Q) =       SU * V2 +
     &                           UT_EMULT_UP(UM,UT,IB) * V1
            ENDDO
          ENDDO
        ENDIF

C  ..(b)

        IF ( N.GT.K .AND. K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            UX_UP = T_UTUP_USERM(UT,UM)
            SU = ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
            DO Q = 1, K_PARAMETERS
              V1 = L_INITIAL_TRANS(N,K,IB,Q) -
     &             ( L_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_P(N,UM,IB) )
              V2 =           L_T_UTDN_MUBAR(UT,K,IB,Q) -
     *               UX_UP * L_T_DELT_MUBAR( N,K,IB,Q)
              L_UT_EMULT_UP(UM,UT,K,IB,Q) =       SU * V2 +
     &                           UT_EMULT_UP(UM,UT,IB) * V1
            ENDDO
          ENDDO
        ENDIF

C  For the plane-parallel case
C  ---------------------------

       ELSE

C  Case (a)

        IF ( K.EQ.N .OR. K.EQ.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            UX_UP = T_UTUP_USERM(UT,UM)
            SU = ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
            DO Q = 1, K_PARAMETERS
              V1 = ZERO
              IF ( K.EQ.0 ) V1 = L_INITIAL_TRANS (N,K,IB,Q)
              V2 =           L_T_UTDN_MUBAR(UT,K,IB,Q) -
     *               UX_UP * L_T_DELT_MUBAR( N,K,IB,Q) -
     *               WDEL  * L_T_UTUP_USERM(UT,UM,Q)
              L_UT_EMULT_UP(UM,UT,K,IB,Q) =       SU * V2 +
     &                           UT_EMULT_UP(UM,UT,IB) * V1
c              L_UT_EMULT_UP(UM,UT,K,IB,Q) =     SU * V2
            ENDDO
          ENDDO
        ENDIF

C  Case (b)

        IF ( N.GT.K .AND. K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              V1 = L_INITIAL_TRANS(N,K,IB,Q)
              L_UT_EMULT_UP(UM,UT,K,IB,Q) = UT_EMULT_UP(UM,UT,IB) * V1
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

      SUBROUTINE L_PARTLAYER_EMULT_DN ( N, UT, K, K_PARAMETERS )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setups and multiplier variables (input)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  include file of linearized setup variables (input)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'

C  linearized multiplier output goes in this file

      INCLUDE '../includes/LIDORT_L_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Given offgrid indices, varying layer index and number of parameters

      INTEGER          N, UT, K, K_PARAMETERS

C  local variables
C  ---------------

      DOUBLE PRECISION SD, V1, V2, V3
      INTEGER          UM, Q, IB

C  Start Beam loop
C  ===============

      DO IB = 1, NBEAMS

C  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
         DO UM = 1, N_USER_STREAMS
           DO Q = 1, K_PARAMETERS
             L_UT_EMULT_DN(UM,UT,K,IB,Q) = ZERO
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

C  Case(a)

        IF ( K.EQ.N .OR. K.EQ.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
              V1 = ONE - OFFGRID_UTAU_VALUES(UT) * USER_SECANTS(UM)
              V2 = - HALF * OFFGRID_UTAU_VALUES(UT)
              DO Q = 1, K_PARAMETERS
                SD = V1 * L_DELTAU_VERT(Q,N) +
     &               V2 * L_AVERAGE_SECANT(N,K,IB,Q)
                V3 = ZERO
                IF ( K.EQ.0 ) V3 = L_INITIAL_TRANS (N,K,IB,Q)                
                L_UT_EMULT_DN(UM,UT,K,IB,Q) = 
     &                   UT_EMULT_DN(UM,UT,IB) * ( SD + V3 )
              ENDDO
            ELSE
              SD = ITRANS_USERM(N,UM,IB) / SIGMA_M(N,UM,IB)
              DO Q = 1, K_PARAMETERS
                V1 = - L_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_M(N,UM,IB)
                IF ( K.EQ.0 ) V1 = V1 + L_INITIAL_TRANS (N,K,IB,Q)
                V2 = L_T_UTDN_USERM(UT,UM,Q) - L_T_UTDN_MUBAR(UT,K,IB,Q)
                L_UT_EMULT_DN(UM,UT,K,IB,Q) =       SD * V2 +
     &                           UT_EMULT_DN(UM,UT,IB) * V1
              ENDDO
            ENDIF
          ENDDO
        ENDIF

C  Case (b), profile only

        IF ( N.GT.K .AND. K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
              V2 = - HALF * OFFGRID_UTAU_VALUES(UT)
              DO Q = 1, K_PARAMETERS
                SD =        L_INITIAL_TRANS (N,K,IB,Q) +
     &                 V2 * L_AVERAGE_SECANT(N,K,IB,Q)
                L_UT_EMULT_DN(UM,UT,K,IB,Q) = UT_EMULT_DN(UM,UT,IB)*SD
              ENDDO
            ELSE
              SD = ITRANS_USERM(N,UM,IB) / SIGMA_M(N,UM,IB)
              DO Q = 1, K_PARAMETERS
                V1 = L_INITIAL_TRANS(N,K,IB,Q) -
     &               ( L_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_M(N,UM,IB) )
                V2 = - L_T_UTDN_MUBAR(UT,K,IB,Q)
                L_UT_EMULT_DN(UM,UT,K,IB,Q) =       SD * V2 +
     &                           UT_EMULT_DN(UM,UT,IB) * V1
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
              V1 = ONE - OFFGRID_UTAU_VALUES(UT) * USER_SECANTS(UM)
              DO Q = 1, K_PARAMETERS
                V3 = ZERO
                IF ( K.EQ.0 ) V3 = L_INITIAL_TRANS (N,K,IB,Q)
                SD = V1 * L_DELTAU_VERT(Q,N)
                L_UT_EMULT_DN(UM,UT,K,IB,Q) =
     &                UT_EMULT_DN(UM,UT,IB)  * ( SD + V3 )
              ENDDO
             ELSE
              SD = ITRANS_USERM(N,UM,IB) / SIGMA_M(N,UM,IB)
              DO Q = 1, K_PARAMETERS
                V1 = ZERO
                IF ( K.EQ.0 ) V1 = L_INITIAL_TRANS (N,K,IB,Q)
                V2 = L_T_UTDN_USERM(UT,UM,Q) - L_T_UTDN_MUBAR(UT,K,IB,Q)
c                L_UT_EMULT_DN(UM,UT,K,IB,Q) =  SD * V2 
                L_UT_EMULT_DN(UM,UT,K,IB,Q) =       SD * V2 +
     &                           UT_EMULT_DN(UM,UT,IB) * V1
              ENDDO
            ENDIF
          ENDDO
        ENDIF

C   Case (b), profile only

        IF ( N.GT.K .AND. K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              V1 = L_INITIAL_TRANS(N,K,IB,Q)
              L_UT_EMULT_DN(UM,UT,K,IB,Q) = UT_EMULT_DN(UM,UT,IB) * V1
            ENDDO
          ENDDO
        ENDIF

C  End clause pseudo-spherical versus plaen-parallel

       ENDIF

C  continuation point for next beam

 5678  CONTINUE

      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE L_QUAD_GFUNCMULT ( IB, UT, N, K, K_PARAMETERS )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setups and multiplier variables (input)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  include files of linearized setup and solution variables (input)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'

C  linearized multiplier output goes in this file

      INCLUDE '../includes/LIDORT_L_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Beam index,  offgrid indices
C  varying layer index and number of parameters

      INTEGER          IB, N, UT, K, K_PARAMETERS

C  local variables
C  ---------------

      INTEGER          Q, AA
      DOUBLE PRECISION SD, SU, T0, TD, TU 
      DOUBLE PRECISION ZX_DN, ZX_UP, ZW, WX, WDEL
      DOUBLE PRECISION L_ZX_DN, L_ZX_UP, L_ZW, L_WX, L_WDEL

C  No particular solution beyond the cutoff layer.
C    [ Zero the multiplier values and exit )

      IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
        DO AA = 1, NSTREAMS
          DO Q = 1, K_PARAMETERS
            L_UT_GMULT_UP(AA,UT,K,Q) = ZERO
            L_UT_GMULT_DN(AA,UT,K,Q) = ZERO
           ENDDO
        ENDDO
        RETURN
      ENDIF

C  Layer constant terms

      WX    = T_UTDN_MUBAR(UT,IB)
      WDEL  = T_DELT_MUBAR(N,IB)

C  For the Pseudo-spherical (average secant) multipliers
C  =====================================================

      IF ( .NOT. DO_PLANE_PARALLEL ) THEN

C  If the layer containing UT is the same as the varying layer
C   Or if we are doing column linearization

        IF ( N.EQ.K .OR. K.EQ.0 ) THEN

          DO AA = 1, NSTREAMS

            ZX_DN = T_UTDN_EIGEN(AA,UT)
            ZX_UP = T_UTUP_EIGEN(AA,UT)
            ZW    = WDEL * ZX_UP

            DO Q = 1, K_PARAMETERS

              L_WDEL = L_T_DELT_MUBAR(N,K,IB,Q)
              L_WX   = L_T_UTDN_MUBAR(UT,K,IB,Q)
              L_ZX_DN = L_T_UTDN_EIGEN(AA,UT,Q)
              L_ZX_UP = L_T_UTUP_EIGEN(AA,UT,Q)
              L_ZW    = WDEL * L_ZX_UP + L_WDEL * ZX_UP
            
              SD  = ( L_ZX_DN - L_WX ) / ( ZX_DN - WX )
              TD = L_ATERM_SAVE(AA,N,Q) + L_GAMMA_M(AA,N,N,Q) + SD
              L_UT_GMULT_DN(AA,UT,K,Q) = TD * UT_GMULT_DN(AA,UT)

              SU  = ( L_WX    - L_ZW ) / ( WX    - ZW )
              TU = L_BTERM_SAVE(AA,N,Q) + L_GAMMA_P(AA,N,N,Q) + SU
              L_UT_GMULT_UP(AA,UT,K,Q) = TU * UT_GMULT_UP(AA,UT)

            ENDDO

          ENDDO

C  If the varying layer is above layer N, profile only

        ELSE IF ( N.GT.K .AND. K.NE.0 ) THEN

          DO Q = 1, K_PARAMETERS

            L_WDEL = L_T_DELT_MUBAR(N,K,IB,Q)
            L_WX   = L_T_UTDN_MUBAR(UT,K,IB,Q)

            DO AA = 1, NSTREAMS

              ZX_DN = T_UTDN_EIGEN(AA,UT)
              ZX_UP = T_UTUP_EIGEN(AA,UT)
              ZW    = WDEL * ZX_UP
              L_ZW  = L_WDEL * ZX_UP
            
              SD  = - L_WX / ( ZX_DN - WX )
              TD = L_GAMMA_M(AA,N,K,Q) + SD + L_INITIAL_TRANS(N,K,IB,Q)
              L_UT_GMULT_DN(AA,UT,K,Q) = TD * UT_GMULT_DN(AA,UT)

              SU  = ( L_WX  - L_ZW ) / ( WX - ZW )
              TU = L_GAMMA_P(AA,N,K,Q) + SU + L_INITIAL_TRANS(N,K,IB,Q)
              L_UT_GMULT_UP(AA,UT,K,Q) = TU * UT_GMULT_UP(AA,UT)

            ENDDO

          ENDDO

        ENDIF

C  Plane parallel case
C  ===================

      ELSE IF ( DO_PLANE_PARALLEL ) THEN

C  If the layer containing UT is the same as the varying layer
C   Or if this is the column linearization

        IF ( N.EQ.K .OR. K.EQ.0 ) THEN

          DO AA = 1, NSTREAMS

            ZX_DN = T_UTDN_EIGEN(AA,UT)
            ZX_UP = T_UTUP_EIGEN(AA,UT)
            ZW    = WDEL * ZX_UP

            DO Q = 1, K_PARAMETERS

              L_WDEL = L_T_DELT_MUBAR(N,N,IB,Q)
              L_WX   = L_T_UTDN_MUBAR(UT,N,IB,Q)
              L_ZX_DN = L_T_UTDN_EIGEN(AA,UT,Q)
              L_ZX_UP = L_T_UTUP_EIGEN(AA,UT,Q)
              L_ZW    = WDEL * L_ZX_UP + L_WDEL * ZX_UP
            
              SD  = ( L_ZX_DN - L_WX ) / ( ZX_DN - WX )
              TD = L_ATERM_SAVE(AA,N,Q) + L_GAMMA_M(AA,N,N,Q) + SD
              L_UT_GMULT_DN(AA,UT,K,Q) = TD * UT_GMULT_DN(AA,UT)

              SU  = ( L_WX - L_ZW ) / ( WX - ZW )
              TU = L_BTERM_SAVE(AA,N,Q) + L_GAMMA_P(AA,N,N,Q) + SU
              L_UT_GMULT_UP(AA,UT,K,Q) = TU * UT_GMULT_UP(AA,UT)

            ENDDO

          ENDDO

C  If the varying layer is above layer N
C  Case where N > K, profile only (K not 0)

        ELSE IF ( N.GT.K .AND. K.NE.0 ) THEN

          DO Q = 1, K_PARAMETERS
            T0 = L_INITIAL_TRANS(N,K,IB,Q)
            DO AA = 1, NSTREAMS            
              L_UT_GMULT_DN(AA,UT,K,Q) = T0 * UT_GMULT_DN(AA,UT)
              L_UT_GMULT_UP(AA,UT,K,Q) = T0 * UT_GMULT_UP(AA,UT)
            ENDDO
          ENDDO

        ENDIF

      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE L_WHOLELAYER_GMULT_UP ( IB, N, K, K_PARAMETERS )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setups, solution and multiplier variables (input)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  include files of linearized setup and solution variables (input)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'

C  linearized multiplier output goes in this file

      INCLUDE '../includes/LIDORT_L_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Beam index, layer index
C  Varying layer and number of parameters

      INTEGER          IB, N, K, K_PARAMETERS

C  Local variables
C  ---------------

      INTEGER          UM, AA, Q
      DOUBLE PRECISION WDEL, L_WDEL, CONST, T1, L_SD, L_SU, LTI

C  No particular solution beyond the cutoff layer.
C    [ Zero the multiplier values and exit )

      IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO AA = 1, NSTREAMS
            DO Q = 1, K_PARAMETERS
              L_SGMULT_UP(AA,UM,Q) = ZERO
              L_SGMULT_DN(AA,UM,Q) = ZERO
            ENDDO
          ENDDO
        ENDDO
        RETURN
      ENDIF

C  transmittance

      WDEL  = T_DELT_MUBAR(N,IB)
      CONST = INITIAL_TRANS(N,IB)

C  For N = K, or K = 0, Quasi-spherical and plane parallel the same.

      IF ( N.EQ.K .OR. K.EQ.0 ) THEN

        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO AA = 1, NSTREAMS
            DO Q = 1, K_PARAMETERS
              L_WDEL = L_T_DELT_MUBAR(N,K,IB,Q)
              L_SD = CONST * L_HMULT_2(AA,UM,N,Q)
              L_SD = L_SD - L_EMULT_UP(UM,N,K,IB,Q)
              T1 = L_ATERM_SAVE(AA,N,Q) + L_GAMMA_M(AA,N,K,Q)
              T1 = T1 * SGMULT_UD(AA,UM,N)
              L_SGMULT_DN(AA,UM,Q) = AGM(AA,N) * L_SD + T1
              IF ( K.EQ.0 ) THEN
                LTI = L_INITIAL_TRANS(N,K,IB,Q) * CONST
                L_SGMULT_DN(AA,UM,Q) = L_SGMULT_DN(AA,UM,Q) + 
     &             AGM(AA,N)*HMULT_2(AA,UM,N) * LTI
              ENDIF
              L_SU = - CONST * ( L_HMULT_1(AA,UM,N,Q) * WDEL +
     &                             HMULT_1(AA,UM,N)   * L_WDEL )
              L_SU =  L_SU + L_EMULT_UP(UM,N,K,IB,Q)
              T1 = L_BTERM_SAVE(AA,N,Q) + L_GAMMA_P(AA,N,K,Q)
              T1 = T1 * SGMULT_UU(AA,UM,N)
              L_SGMULT_UP(AA,UM,Q) = BGP(AA,N) * L_SU + T1
              IF ( K.EQ.0 ) THEN
                LTI = L_INITIAL_TRANS(N,K,IB,Q) * CONST
                L_SGMULT_UP(AA,UM,Q) = L_SGMULT_UP(AA,UM,Q) - 
     &               WDEL*BGP(AA,N)*HMULT_1(AA,UM,N) * LTI
              ENDIF
            ENDDO
          ENDDO
        ENDDO

C  Case where N > K, profile only (K not 0)

      ELSE IF ( N.GT.K .AND. K.NE.0 ) THEN

C  Pseudo-spherical

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO AA = 1, NSTREAMS
              DO Q = 1, K_PARAMETERS
                L_SD = HELP_AQ(N,K,IB,Q) * HMULT_2(AA,UM,N) -
     &                            L_EMULT_UP(UM,N,K,IB,Q)
                T1 = L_GAMMA_M(AA,N,K,Q) * SGMULT_UD(AA,UM,N)
                L_SGMULT_DN(AA,UM,Q) = AGM(AA,N) * L_SD + T1
                L_SU = HELP_BQ(N,K,IB,Q) * HMULT_1(AA,UM,N) +
     &                             L_EMULT_UP(UM,N,K,IB,Q)
                T1 = L_GAMMA_P(AA,N,K,Q) * SGMULT_UU(AA,UM,N)
                L_SGMULT_UP(AA,UM,Q) = BGP(AA,N) * L_SU + T1
              ENDDO
            ENDDO
          ENDDO
        ENDIF

C  Plane parallel case

        IF ( DO_PLANE_PARALLEL ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO AA = 1, NSTREAMS
              DO Q = 1, K_PARAMETERS
                T1 = L_INITIAL_TRANS(N,K,IB,Q)
                L_SGMULT_DN(AA,UM,Q) = T1 * SGMULT_UD(AA,UM,N)
                L_SGMULT_UP(AA,UM,Q) = T1 * SGMULT_UU(AA,UM,N)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE L_WHOLELAYER_GMULT_DN ( IB, N, K, K_PARAMETERS )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setups, solution and multiplier variables (input)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  include files of linearized setup and solution variables (input)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'

C  linearized multiplier output goes in this file

      INCLUDE '../includes/LIDORT_L_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Beam index, layer index
C  Varying layer and number of parameters

      INTEGER          IB, N, K, K_PARAMETERS

C  Local variables
C  ---------------

      INTEGER          UM, AA, Q
      DOUBLE PRECISION WDEL, L_WDEL, CONST, T1, L_SD, L_SU, LTI

C  No particular solution beyond the cutoff layer.
C    [ Zero the multiplier values and exit )

      IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO AA = 1, NSTREAMS
            DO Q = 1, K_PARAMETERS
              L_SGMULT_UP(AA,UM,Q) = ZERO
              L_SGMULT_DN(AA,UM,Q) = ZERO
            ENDDO
          ENDDO
        ENDDO
        RETURN
      ENDIF

C  transmittance

      WDEL  = T_DELT_MUBAR(N,IB)
      CONST = INITIAL_TRANS(N,IB)

C  For N = K, or K = 0, Quasi-spherical and plane parallel the same.

      IF ( N.EQ.K .OR. K.EQ.0 ) THEN

        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO AA = 1, NSTREAMS
            DO Q = 1, K_PARAMETERS
              L_WDEL = L_T_DELT_MUBAR(N,K,IB,Q)
              L_SD = CONST * L_HMULT_1(AA,UM,N,Q)
              L_SD = L_SD - L_EMULT_DN(UM,N,K,IB,Q)
              T1 = L_ATERM_SAVE(AA,N,Q) + L_GAMMA_M(AA,N,K,Q)
              T1 = T1 * SGMULT_DD(AA,UM,N)
              L_SGMULT_DN(AA,UM,Q) = AGM(AA,N) * L_SD + T1
              IF ( K.EQ.0 ) THEN
                LTI = L_INITIAL_TRANS(N,K,IB,Q) * CONST
                L_SGMULT_DN(AA,UM,Q) = L_SGMULT_DN(AA,UM,Q) + 
     &            AGM(AA,N)*HMULT_1(AA,UM,N) * LTI
              ENDIF
              L_SU = - CONST * ( L_HMULT_2(AA,UM,N,Q) * WDEL +
     &                             HMULT_2(AA,UM,N)   * L_WDEL )
              L_SU =  L_SU + L_EMULT_DN(UM,N,K,IB,Q)
              T1 = L_BTERM_SAVE(AA,N,Q) + L_GAMMA_P(AA,N,K,Q)
              T1 = T1 * SGMULT_DU(AA,UM,N)
              L_SGMULT_UP(AA,UM,Q) = BGP(AA,N) * L_SU + T1
              IF ( K.EQ.0 ) THEN
                LTI = L_INITIAL_TRANS(N,K,IB,Q) * CONST
                L_SGMULT_UP(AA,UM,Q) = L_SGMULT_UP(AA,UM,Q) - 
     &            WDEL*BGP(AA,N)*HMULT_2(AA,UM,N) * LTI
              ENDIF
            ENDDO
          ENDDO
        ENDDO

C  Case where N > K, profile only (K not 0)

      ELSE IF ( N.GT.K .AND. K.NE.0 ) THEN

C  Pseudo-spherical

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO AA = 1, NSTREAMS
              DO Q = 1, K_PARAMETERS
                L_SD = HELP_AQ(N,K,IB,Q) * HMULT_1(AA,UM,N) -
     &                            L_EMULT_DN(UM,N,K,IB,Q)
                T1 = L_GAMMA_M(AA,N,K,Q) * SGMULT_DD(AA,UM,N)
                L_SGMULT_DN(AA,UM,Q) = AGM(AA,N) * L_SD + T1
                L_SU = HELP_BQ(N,K,IB,Q) * HMULT_2(AA,UM,N) +
     &                             L_EMULT_DN(UM,N,K,IB,Q)
                T1 = L_GAMMA_P(AA,N,K,Q) * SGMULT_DU(AA,UM,N)
                L_SGMULT_UP(AA,UM,Q) = BGP(AA,N) * L_SU + T1
              ENDDO
            ENDDO
          ENDDO
        ENDIF

C  Plane parallel case

        IF ( DO_PLANE_PARALLEL ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO AA = 1, NSTREAMS
              DO Q = 1, K_PARAMETERS
                T1 = L_INITIAL_TRANS(N,K,IB,Q)
                L_SGMULT_DN(AA,UM,Q) = T1 * SGMULT_DD(AA,UM,N)
                L_SGMULT_UP(AA,UM,Q) = T1 * SGMULT_DU(AA,UM,N)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE L_PARTLAYER_GMULT_UP ( IB, UT, N, K, K_PARAMETERS )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setups, solution and multiplier variables (input)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  include files of linearized setup and solution variables (input)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'

C  linearized multiplier output goes in this file

      INCLUDE '../includes/LIDORT_L_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Beam index, offgrid indices
C  Varying layer and number of parameters

      INTEGER          IB, N, UT, K, K_PARAMETERS

C  Local variables
C  ---------------

      INTEGER          UM, AA, Q
      DOUBLE PRECISION WDEL, L_WDEL, CONST, T1, L_SD, L_SU, LTI

C  No particular solution beyond the cutoff layer.
C    [ Zero the multiplier values and exit )

      IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO AA = 1, NSTREAMS
            DO Q = 1, K_PARAMETERS
              L_SGMULT_UP(AA,UM,Q) = ZERO
              L_SGMULT_DN(AA,UM,Q) = ZERO
            ENDDO
          ENDDO
        ENDDO
        RETURN
      ENDIF

C  transmittance

      WDEL  = T_DELT_MUBAR(N,IB)
      CONST = INITIAL_TRANS(N,IB)

C  For N = K, or K = 0, Pseudo-spherical and plane parallel the same.

      IF ( N.EQ.K .OR. K.EQ.0 ) THEN

        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO AA = 1, NSTREAMS
            DO Q = 1, K_PARAMETERS
              L_WDEL = L_T_DELT_MUBAR(N,K,IB,Q)
              L_SD = CONST * L_UT_HMULT_UD(AA,UM,UT,Q)
              L_SD = L_SD -  L_UT_EMULT_UP(UM,UT,K,IB,Q)
              T1 = L_ATERM_SAVE(AA,N,Q) + L_GAMMA_M(AA,N,K,Q)
              T1 = T1 * UT_SGMULT_UD(AA,UM,UT)
              L_SGMULT_DN(AA,UM,Q) = AGM(AA,N) * L_SD + T1
              IF ( K.EQ.0 ) THEN
                LTI = L_INITIAL_TRANS(N,K,IB,Q) * CONST
                L_SGMULT_DN(AA,UM,Q) = L_SGMULT_DN(AA,UM,Q) + 
     &            AGM(AA,N)*UT_HMULT_UD(AA,UM,UT) * LTI
              ENDIF
              L_SU = - CONST * ( L_UT_HMULT_UU(AA,UM,UT,Q) * WDEL +
     &                             UT_HMULT_UU(AA,UM,UT)   * L_WDEL )
              L_SU =  L_SU + L_UT_EMULT_UP(UM,UT,K,IB,Q)
              T1 = L_BTERM_SAVE(AA,N,Q) + L_GAMMA_P(AA,N,K,Q)
              T1 = T1 * UT_SGMULT_UU(AA,UM,UT)
              L_SGMULT_UP(AA,UM,Q) = BGP(AA,N) * L_SU + T1
              IF ( K.EQ.0 ) THEN
                LTI = L_INITIAL_TRANS(N,K,IB,Q) * CONST
                L_SGMULT_UP(AA,UM,Q) = L_SGMULT_UP(AA,UM,Q) - 
     &            WDEL*BGP(AA,N)*UT_HMULT_UU(AA,UM,UT) * LTI
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF

C  Case where N > K, profile only (K not 0)

      IF ( N.GT.K .AND. K.NE.0 ) THEN

C  Pseudo-spherical

        IF ( .NOT.DO_PLANE_PARALLEL ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO AA = 1, NSTREAMS
              DO Q = 1, K_PARAMETERS
                L_SD = HELP_AQ(N,K,IB,Q) * UT_HMULT_UD(AA,UM,UT) - 
     &                            L_UT_EMULT_UP(UM,UT,K,IB,Q)
                T1 = L_GAMMA_M(AA,N,K,Q) * UT_SGMULT_UD(AA,UM,UT)
                L_SGMULT_DN(AA,UM,Q) = AGM(AA,N) * L_SD + T1
                L_SU = HELP_BQ(N,K,IB,Q) * UT_HMULT_UU(AA,UM,UT) +
     &                            L_UT_EMULT_UP(UM,UT,K,IB,Q)
                T1 = L_GAMMA_P(AA,N,K,Q) * UT_SGMULT_UU(AA,UM,UT)
                L_SGMULT_UP(AA,UM,Q) = BGP(AA,N) * L_SU + T1
              ENDDO
            ENDDO
          ENDDO
        ENDIF

C  Plane parallel case

        IF ( DO_PLANE_PARALLEL ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO AA = 1, NSTREAMS
              DO Q = 1, K_PARAMETERS
                T1 = L_INITIAL_TRANS(N,K,IB,Q)
                L_SGMULT_DN(AA,UM,Q) = T1 * UT_SGMULT_UD(AA,UM,UT)
                L_SGMULT_UP(AA,UM,Q) = T1 * UT_SGMULT_UU(AA,UM,UT)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE L_PARTLAYER_GMULT_DN ( IB, UT, N, K, K_PARAMETERS )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setups, solution and multiplier variables (input)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  include files of linearized setup and solution variables (input)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'

C  linearized multiplier output goes in this file

      INCLUDE '../includes/LIDORT_L_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Beam index, offgrid indices
C  Varying layer and number of parameters

      INTEGER          IB, N, UT, K, K_PARAMETERS

C  Local variables
C  ---------------

      INTEGER          UM, AA, Q
      DOUBLE PRECISION WDEL, L_WDEL, CONST, T1, L_SD, L_SU, LTI

C  No particular solution beyond the cutoff layer.
C    [ Zero the multiplier values and exit )

      IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO AA = 1, NSTREAMS
            DO Q = 1, K_PARAMETERS
              L_SGMULT_UP(AA,UM,Q) = ZERO
              L_SGMULT_DN(AA,UM,Q) = ZERO
            ENDDO
          ENDDO
        ENDDO
        RETURN
      ENDIF

C  transmittance

      WDEL  = T_DELT_MUBAR(N,IB)
      CONST = INITIAL_TRANS(N,IB)

C  For N = K, or K = 0, Pseudo-spherical and plane parallel the same.

      IF ( N.EQ.K .OR. K.EQ.0 ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO AA = 1, NSTREAMS
            DO Q = 1, K_PARAMETERS
              L_WDEL = L_T_DELT_MUBAR(N,K,IB,Q)
              L_SD = CONST * L_UT_HMULT_DD(AA,UM,UT,Q)
              L_SD = L_SD -  L_UT_EMULT_DN(UM,UT,K,IB,Q)
              T1 = L_ATERM_SAVE(AA,N,Q) + L_GAMMA_M(AA,N,K,Q)
              T1 = T1 * UT_SGMULT_DD(AA,UM,UT)
              L_SGMULT_DN(AA,UM,Q) = AGM(AA,N) * L_SD + T1
              IF ( K.EQ.0 ) THEN
                LTI = L_INITIAL_TRANS(N,K,IB,Q) * CONST
                L_SGMULT_DN(AA,UM,Q) = L_SGMULT_DN(AA,UM,Q) + 
     &            AGM(AA,N)*UT_HMULT_DD(AA,UM,UT) * LTI
              ENDIF
              L_SU = - CONST * ( L_UT_HMULT_DU(AA,UM,UT,Q) * WDEL +
     &                             UT_HMULT_DU(AA,UM,UT)   * L_WDEL )
              L_SU =  L_SU + L_UT_EMULT_DN(UM,UT,K,IB,Q)
              T1 = L_BTERM_SAVE(AA,N,Q) + L_GAMMA_P(AA,N,K,Q)
              T1 = T1 * UT_SGMULT_DU(AA,UM,UT)
              L_SGMULT_UP(AA,UM,Q) = BGP(AA,N) * L_SU + T1
              IF ( K.EQ.0 ) THEN
                LTI = L_INITIAL_TRANS(N,K,IB,Q) * CONST
                L_SGMULT_UP(AA,UM,Q) = L_SGMULT_UP(AA,UM,Q) - 
     &            WDEL*BGP(AA,N)*UT_HMULT_DU(AA,UM,UT) * LTI
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF

C  Case where N > K, profile only (K not 0)

      IF ( N.GT.K .AND. K.NE.0 ) THEN

C  Pseudo-spherical

        IF ( .NOT.DO_PLANE_PARALLEL ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO AA = 1, NSTREAMS
              DO Q = 1, K_PARAMETERS
                L_SD = HELP_AQ(N,K,IB,Q) * UT_HMULT_DD(AA,UM,UT) - 
     &                            L_UT_EMULT_DN(UM,UT,K,IB,Q)
                T1 = L_GAMMA_M(AA,N,K,Q) * UT_SGMULT_DD(AA,UM,UT)
                L_SGMULT_DN(AA,UM,Q) = AGM(AA,N) * L_SD + T1
                L_SU = HELP_BQ(N,K,IB,Q) * UT_HMULT_DU(AA,UM,UT) +
     &                            L_UT_EMULT_DN(UM,UT,K,IB,Q)
                T1 = L_GAMMA_P(AA,N,K,Q) * UT_SGMULT_DU(AA,UM,UT)
                L_SGMULT_UP(AA,UM,Q) = BGP(AA,N) * L_SU + T1
              ENDDO
            ENDDO
          ENDDO
        ENDIF

C  Plane parallel case

        IF ( DO_PLANE_PARALLEL ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO AA = 1, NSTREAMS
              DO Q = 1, K_PARAMETERS
                T1 = L_INITIAL_TRANS(N,K,IB,Q)
                L_SGMULT_DN(AA,UM,Q) = T1 * UT_SGMULT_DD(AA,UM,UT)
                L_SGMULT_UP(AA,UM,Q) = T1 * UT_SGMULT_DU(AA,UM,UT)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

      ENDIF

C  Finish

      RETURN
      END
