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
C # Multiplier Subroutines in this Module                       #
C #                                                             #
C #   1. For the homogeneous solutions ------                   #
C #              HMULT_MASTER (master)                          #
C #                                                             #
C #   2. For the solar beam source  ------                      #
C #              EMULT_MASTER (master)                          #
C #                                                             #
C #   3. For the  Green's function ------                       #
C #              QUAD_GFUNCMULT                                 #
C #              WHOLELAYER_GMULT_UP                            #
C #              WHOLELAYER_GMULT_DN                            #
C #              PARTLAYER_GMULT_UP                             #
C #              PARTLAYER_GMULT_DN                             #
C #                                                             #
C ###############################################################

      SUBROUTINE HMULT_MASTER

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup/solution stuff (inputs to this module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  include file of Multiplier variables (output to this module)

      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  Local variables
C  ---------------

      INTEGER          UM, K, N, UT, UTA
      DOUBLE PRECISION UDEL, SM
      DOUBLE PRECISION ZDEL, ZUDEL, THETA_1, THETA_2
      DOUBLE PRECISION UX_UP, UX_DN, ZX_UP, ZX_DN
      DOUBLE PRECISION THETA_DN, THETA_UP
      LOGICAL          STERM_EXIST(MAXLAYERS)

C  Existence flags

      DO N = 1, NLAYERS
        STERM_EXIST(N) = 
     &   ( STERM_LAYERMASK_UP(N).OR.STERM_LAYERMASK_DN(N) )
      ENDDO

C  whole layer multipliers
C  -----------------------

C  Start loops over layers and user-streams
C    Only done if layers are flagged

      DO N = 1, NLAYERS
        IF ( STERM_EXIST(N) ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            UDEL = T_DELT_USERM(N,UM)
            SM   = USER_SECANTS(UM)
            DO K = 1, NSTREAMS
              ZDEL    = T_DELT_EIGEN(K,N)
              ZUDEL   = ZDEL * UDEL
              THETA_2 = ONE    - ZUDEL
              THETA_1 = ZDEL - UDEL
              HMULT_1(K,UM,N) = SM * THETA_1 * ZETA_M(K,UM,N)
              HMULT_2(K,UM,N) = SM * THETA_2 * ZETA_P(K,UM,N)
            ENDDO
          ENDDO
        ENDIF
      ENDDO

C  partial layer multipliers
C  -------------------------

      DO UTA = 1, N_OUT_USERTAUS
        IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN
          UT = OFFGRID_UTAU_OUTINDEX(UTA)
          N  = OFFGRID_UTAU_LAYERIDX(UT)
          IF ( STERM_EXIST(N) ) THEN

C  Upwelling

            IF ( DO_UPWELLING ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                UX_UP = T_UTUP_USERM(UT,UM)
                SM    = USER_SECANTS(UM)
                DO K = 1, NSTREAMS
                  ZDEL     = T_DELT_EIGEN(K,N)
                  ZX_UP    = T_UTUP_EIGEN(K,UT)
                  ZX_DN    = T_UTDN_EIGEN(K,UT)
                  THETA_DN = ZX_DN - ZDEL * UX_UP
                  THETA_UP = ZX_UP - UX_UP
                  UT_HMULT_UD(K,UM,UT) = SM * THETA_DN * ZETA_P(K,UM,N)
                  UT_HMULT_UU(K,UM,UT) = SM * THETA_UP * ZETA_M(K,UM,N)
                ENDDO
              ENDDO
            ENDIF

C  Downwelling

            IF ( DO_DNWELLING ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                UX_DN = T_UTDN_USERM(UT,UM)
                SM    = USER_SECANTS(UM)
                DO K = 1, NSTREAMS
                  ZDEL     = T_DELT_EIGEN(K,N)
                  ZX_UP    = T_UTUP_EIGEN(K,UT)
                  ZX_DN    = T_UTDN_EIGEN(K,UT)
                  THETA_DN = ZX_DN - UX_DN
                  THETA_UP = ZX_UP - ZDEL * UX_DN
                  UT_HMULT_DD(K,UM,UT) = SM * THETA_DN * ZETA_M(K,UM,N)
                  UT_HMULT_DU(K,UM,UT) = SM * THETA_UP * ZETA_P(K,UM,N)
                ENDDO
              ENDDO
            ENDIF

C  end loop over partial layers

          ENDIF
        ENDIF
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE EMULT_MASTER

C  Prepare multipliers for the Beam source terms

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup and solution variables

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  include file of Multiplier variables (output stored here)

      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  local variables
C  ---------------

      INTEGER          N, UT, UM, IB
      DOUBLE PRECISION WDEL, WX, WUDEL, UDEL
      DOUBLE PRECISION DIFF, SB, SECMUM, SU, SD
      DOUBLE PRECISION UX_DN, UX_UP, WDEL_UXUP

C  L'Hopital's Rule flags for Downwelling EMULT
C  --------------------------------------------

      IF ( DO_DNWELLING ) THEN
       DO N = 1, NLAYERS
        IF ( STERM_LAYERMASK_DN(N) ) THEN
         DO IB = 1, NBEAMS
          SB = AVERAGE_SECANT(N,IB)
          DO UM = 1, N_USER_STREAMS
            DIFF = DABS ( USER_SECANTS(UM) - SB )
            IF ( DIFF .LT. HOPITAL_TOLERANCE ) THEN
              EMULT_HOPRULE(N,UM,IB) = .TRUE.
            ELSE
              EMULT_HOPRULE(N,UM,IB) = .FALSE.
            ENDIF
          ENDDO
         ENDDO
        ENDIF
       ENDDO
      ENDIF

C  sigma functions (all layers)
C  ----------------------------

      DO N = 1, NLAYERS
       DO IB = 1, NBEAMS
        SB = AVERAGE_SECANT(N,IB)
        DO UM = 1, N_USER_STREAMS
          SECMUM = USER_SECANTS(UM)
          SIGMA_P(N,UM,IB) = SB + SECMUM
          SIGMA_M(N,UM,IB) = SB - SECMUM
        ENDDO
       ENDDO
      ENDDO

C  upwelling External source function multipliers
C  ----------------------------------------------

      IF ( DO_UPWELLING ) THEN

C  whole layer

        DO N = 1, NLAYERS
         IF ( STERM_LAYERMASK_UP(N) ) THEN
          DO IB = 1, NBEAMS
            IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
              DO UM = 1, N_USER_STREAMS
                EMULT_UP(UM,N,IB) = ZERO
              ENDDO
            ELSE
              WDEL = T_DELT_MUBAR(N,IB)
              DO UM = 1, N_USER_STREAMS
                WUDEL = WDEL * T_DELT_USERM(N,UM)
                SU = ( ONE - WUDEL ) / SIGMA_P(N,UM,IB)
                EMULT_UP(UM,N,IB) = ITRANS_USERM(N,UM,IB) * SU
              ENDDO
            ENDIF
          ENDDO
         ENDIF
        ENDDO

C  Partial layer

        DO UT = 1, N_OFFGRID_USERTAUS
         N  = OFFGRID_UTAU_LAYERIDX(UT)
         DO IB = 1, NBEAMS
          IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
            DO UM = 1, N_USER_STREAMS
              UT_EMULT_UP(UM,UT,IB) = ZERO
            ENDDO
          ELSE
            WX   = T_UTDN_MUBAR(UT,IB)
            WDEL = T_DELT_MUBAR(N,IB)
            DO UM = 1, N_USER_STREAMS
              UX_UP = T_UTUP_USERM(UT,UM)
              WDEL_UXUP = UX_UP * WDEL
              SU = ( WX - WDEL_UXUP ) / SIGMA_P(N,UM,IB)
              UT_EMULT_UP(UM,UT,IB) = ITRANS_USERM(N,UM,IB) * SU
            ENDDO
          ENDIF
         ENDDO
        ENDDO

      ENDIF

C  downwelling External source function multipliers
C  ------------------------------------------------

C    .. Note use of L'Hopitals Rule

      IF ( DO_DNWELLING ) THEN

C  whole layer

        DO N = 1, NLAYERS
         IF ( STERM_LAYERMASK_DN(N) ) THEN
          DO IB = 1, NBEAMS
            IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
              DO UM = 1, N_USER_STREAMS
                EMULT_DN(UM,N,IB) = ZERO
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
         ENDIF
        ENDDO

C  Partial layer

        DO UT = 1, N_OFFGRID_USERTAUS
         N  = OFFGRID_UTAU_LAYERIDX(UT)
         DO IB = 1, NBEAMS
          IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
            DO UM = 1, N_USER_STREAMS
              UT_EMULT_DN(UM,UT,IB) = ZERO
            ENDDO
          ELSE
            WX   = T_UTDN_MUBAR(UT,IB)
            DO UM = 1, N_USER_STREAMS
              UX_DN = T_UTDN_USERM(UT,UM)
              IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
                SD = OFFGRID_UTAU_VALUES(UT) * UX_DN
              ELSE
                SD = ( UX_DN - WX ) / SIGMA_M(N,UM,IB)
              ENDIF
              UT_EMULT_DN(UM,UT,IB) = ITRANS_USERM(N,UM,IB) * SD
            ENDDO
          ENDIF
         ENDDO
        ENDDO

      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE QUAD_GFUNCMULT ( IB, UT, N )

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup and solution variables (input to this module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  include file of Multiplier variables (output to this module)

      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

      INTEGER          N, UT, IB

C  Local variables
C  ---------------

      INTEGER          AA
      DOUBLE PRECISION SD, SU, WDEL 
      DOUBLE PRECISION ZX_DN, ZX_UP, ZW, WX, CONST

C  No particular solution beyond the cutoff layer.
C    [ Zero the multiplier values and exit )

      IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
        DO AA = 1, NSTREAMS
          UT_GMULT_DN(AA,UT) = ZERO
          UT_GMULT_UP(AA,UT) = ZERO
        ENDDO
        RETURN
      ENDIF
        
C  Layer constant terms

      WX    = T_UTDN_MUBAR(UT,IB)
      WDEL  = T_DELT_MUBAR(N,IB)
      CONST = INITIAL_TRANS(N,IB)

C  Tau integration without coefficients (average secant approximation)

      DO AA = 1, NSTREAMS
        ZX_DN = T_UTDN_EIGEN(AA,UT)
        ZX_UP = T_UTUP_EIGEN(AA,UT)
        ZW    = WDEL * ZX_UP
        SD =  ( ZX_DN - WX ) * GAMMA_M(AA,N)
        SU =  ( WX    - ZW ) * GAMMA_P(AA,N)
        UT_GMULT_DN(AA,UT) = SD * ATERM_SAVE(AA,N) * CONST
        UT_GMULT_UP(AA,UT) = SU * BTERM_SAVE(AA,N) * CONST
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE WHOLELAYER_GMULT_UP ( N, IB )

C  INTEGRATED GREEN FUNCTION MULTIPLIER (Upwelling Whole layer)

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup and solution stuff (input to this module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  include file of Multiplier variables (output to this module)

      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer and beam index

      INTEGER          N, IB

C  Local variables
C  ---------------

      INTEGER          UM, AA
      DOUBLE PRECISION WDEL, CONS, CONSWDEL, SD, SU, TD, TU

C  No particular solution beyond the cutoff layer.
C    [ Zero the multiplier values and exit )

      IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO AA = 1, NSTREAMS
            SGMULT_UD(AA,UM,N) = ZERO
            SGMULT_UU(AA,UM,N) = ZERO
          ENDDO
        ENDDO
        RETURN
      ENDIF

C  Layer quantities

      WDEL = T_DELT_MUBAR(N,IB)
      CONS = INITIAL_TRANS(N,IB)
      CONSWDEL = - CONS * WDEL

C  Get the multipliers

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO AA = 1, NSTREAMS
          TD = CONS     * HMULT_2(AA,UM,N)
          TU = CONSWDEL * HMULT_1(AA,UM,N)
          SD  = GAMMA_M(AA,N) * ( TD - EMULT_UP(UM,N,IB) )
          SU  = GAMMA_P(AA,N) * ( TU + EMULT_UP(UM,N,IB) )
          SGMULT_UD(AA,UM,N) = SD * ATERM_SAVE(AA,N)
          SGMULT_UU(AA,UM,N) = SU * BTERM_SAVE(AA,N)
        ENDDO
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE WHOLELAYER_GMULT_DN ( N, IB )

C  INTEGRATED GREEN FUNCTION MULTIPLIER (Downwelling Whole layer)

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup and solution stuff (input to this module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  include file of Multiplier variables (output to this module)

      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer and beam index

      INTEGER          N, IB

C  Local variables
C  ---------------

      INTEGER          UM, AA
      DOUBLE PRECISION WDEL, CONS, CONSWDEL, SD, SU, TD, TU

C  No particular solution beyond the cutoff layer.
C    [ Zero the multiplier values and exit )

      IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO AA = 1, NSTREAMS
            SGMULT_DD(AA,UM,N) = ZERO
            SGMULT_DU(AA,UM,N) = ZERO
          ENDDO
        ENDDO
        RETURN
      ENDIF

C  Layer quantities

      WDEL = T_DELT_MUBAR(N,IB)
      CONS = INITIAL_TRANS(N,IB)
      CONSWDEL = - CONS * WDEL

C  Get the basic multipliers

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO AA = 1, NSTREAMS
          TD = CONS     * HMULT_1(AA,UM,N)
          TU = CONSWDEL * HMULT_2(AA,UM,N)
          SD  = GAMMA_M(AA,N) * ( TD - EMULT_DN(UM,N,IB) )
          SU  = GAMMA_P(AA,N) * ( TU + EMULT_DN(UM,N,IB) )
          SGMULT_DD(AA,UM,N) = SD * ATERM_SAVE(AA,N)
          SGMULT_DU(AA,UM,N) = SU * BTERM_SAVE(AA,N)
        ENDDO
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE PARTLAYER_GMULT_UP ( N, UT, IB )

C  INTEGRATED GREEN FUNCTION MULTIPLIER (Upwelling partial layer)

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup and solution stuff (input to this module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  include file of Multiplier variables (output to this module)

      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index, offgrid index, beam index

      INTEGER          N, UT, IB

C  Local variables
C  ---------------

      INTEGER          UM, AA
      DOUBLE PRECISION WDEL, CONS, CONSWDEL, SD, SU, TD, TU

C  No particular solution beyond the cutoff layer.
C    [ Zero the multiplier values and exit )

      IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO AA = 1, NSTREAMS
            UT_SGMULT_UD(AA,UM,UT) = ZERO
            UT_SGMULT_UU(AA,UM,UT) = ZERO
          ENDDO
        ENDDO
        RETURN
      ENDIF

C  Layer quantities

      WDEL = T_DELT_MUBAR(N,IB)
      CONS = INITIAL_TRANS(N,IB)
      CONSWDEL = - CONS * WDEL

C  Get the multipliers

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO AA = 1, NSTREAMS
          TD =   CONS     * UT_HMULT_UD(AA,UM,UT)
          TU =   CONSWDEL * UT_HMULT_UU(AA,UM,UT)
          SD  = GAMMA_M(AA,N) * ( TD - UT_EMULT_UP(UM,UT,IB) )
          SU  = GAMMA_P(AA,N) * ( TU + UT_EMULT_UP(UM,UT,IB) )
          UT_SGMULT_UD(AA,UM,UT) = SD * ATERM_SAVE(AA,N)
          UT_SGMULT_UU(AA,UM,UT) = SU * BTERM_SAVE(AA,N)
        ENDDO
       ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE PARTLAYER_GMULT_DN ( N, UT, IB )

C  INTEGRATED GREEN FUNCTION MULTIPLIER (Downwelling partial layer)

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup and solution stuff (input to this module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  include file of Multiplier variables (output to this module)

      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index, offgrid index, beam index

      INTEGER          N, UT, IB

C  Local variables
C  ---------------

      INTEGER          UM, AA
      DOUBLE PRECISION WDEL, CONS, CONSWDEL, SD, SU, TD, TU

C  No particular solution beyond the cutoff layer.
C    [ Zero the multiplier values and exit )

      IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO AA = 1, NSTREAMS
            UT_SGMULT_DD(AA,UM,UT) = ZERO
            UT_SGMULT_DU(AA,UM,UT) = ZERO
          ENDDO
        ENDDO
        RETURN
      ENDIF

C  Layer quantities

      WDEL = T_DELT_MUBAR(N,IB)
      CONS = INITIAL_TRANS(N,IB)
      CONSWDEL = - CONS * WDEL

C  Get the multipliers

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO AA = 1, NSTREAMS
          TD =   CONS     * UT_HMULT_DD(AA,UM,UT)
          TU =   CONSWDEL * UT_HMULT_DU(AA,UM,UT)
          SD  = GAMMA_M(AA,N) * ( TD - UT_EMULT_DN(UM,UT,IB) )
          SU  = GAMMA_P(AA,N) * ( TU + UT_EMULT_DN(UM,UT,IB) )
          UT_SGMULT_DD(AA,UM,UT) = SD * ATERM_SAVE(AA,N)
          UT_SGMULT_DU(AA,UM,UT) = SU * BTERM_SAVE(AA,N)
        ENDDO
      ENDDO

C  Finish

      RETURN
      END

