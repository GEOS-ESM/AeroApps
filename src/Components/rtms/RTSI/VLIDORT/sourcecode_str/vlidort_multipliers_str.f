C ###############################################################
C #                                                             #
C #                    THE VECTOR LIDORT MODEL                  #
C #                                                             #
C #  (Vector LInearized Discrete Ordinate Radiative Transfer)   #
C #   -      --         -        -        -         -           #
C #                                                             #
C ###############################################################

C ###############################################################
C #                                                             #
C #  Author :      Robert. J. D. Spurr                          #
C #                                                             #
C #  Address :      RT Solutions, inc.                          #
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
C #              HMULT_MASTER (master)                          #
C #                WHOLELAYER_HMULT                             #
C #                PARTLAYER_HMULT_UP                           #
C #                PARTLAYER_HMULT_DN                           #
C #                                                             #
C #              EMULT_MASTER (master)                          #
C #                                                             #
C ###############################################################

      SUBROUTINE HMULT_MASTER

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup/solution stuff (inputs to this module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'

C  include files of Multiplier variables (output to this module)

      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  Local variables
C  ---------------

      INTEGER          N, UT, UTA

C  whole layer multipliers

      CALL WHOLELAYER_HMULT

C  partial layer multipliers

      DO UTA = 1, N_USER_LEVELS
        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)
          IF ( DO_UPWELLING ) CALL PARTLAYER_HMULT_UP ( N, UT )
          IF ( DO_DNWELLING ) CALL PARTLAYER_HMULT_DN ( N, UT )
        ENDIF
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE WHOLELAYER_HMULT

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup/solution stuff (inputs to this module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'

C  include files of Multiplier variables (output to this module)

      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  Local variables
C  ---------------

      INTEGER          UM, K, N, KO1, K0, K1, K2
      DOUBLE PRECISION UDEL, SM
      DOUBLE PRECISION ZDEL_R, ZUDEL_R, THETA_1_R, THETA_2_R
      DOUBLE PRECISION ZDEL_CR,    ZDEL_CI
      DOUBLE PRECISION THETA_1_CR, THETA_1_CI
      DOUBLE PRECISION THETA_2_CR, THETA_2_CI

C  Start loops over layers and user-streams
C   Only done if layers are flagged

C  Small numbers analysis, added 30 October 2007

      DO N = 1, NLAYERS
        IF ( STERM_LAYERMASK_UP(N).OR.STERM_LAYERMASK_DN(N) ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS

C  setup

            UDEL = T_DELT_USERM(N,UM)
            SM = USER_SECANTS(UM)

C  Get the real multipliers

            DO K = 1, K_REAL(N)
              ZDEL_R    = T_DELT_EIGEN(K,N)
              ZUDEL_R   = ZDEL_R * UDEL
              THETA_2_R = ONE    - ZUDEL_R
              THETA_1_R = ZDEL_R - UDEL
              IF ( HSINGO(K,UM,N) )THEN
                HMULT_1(K,UM,N) = SM * UDEL * DELTAU_VERT(N)
              ELSE
                HMULT_1(K,UM,N) = SM * THETA_1_R * ZETA_M(K,UM,N)
              ENDIF
              HMULT_2(K,UM,N) = SM * THETA_2_R * ZETA_P(K,UM,N)
            ENDDO

C  Get the complex multipliers

            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              ZDEL_CR    = T_DELT_EIGEN(K1,N)
              ZDEL_CI    = T_DELT_EIGEN(K2,N)
              THETA_2_CR = ONE - ZDEL_CR * UDEL
              THETA_2_CI =     - ZDEL_CI * UDEL
              THETA_1_CR = ZDEL_CR - UDEL
              THETA_1_CI = ZDEL_CI
              IF ( HSINGO(K1,UM,N) )THEN
                HMULT_1(K1,UM,N) = SM * UDEL * DELTAU_VERT(N)
                HMULT_1(K2,UM,N) = ZERO
              ELSE
                HMULT_1(K1,UM,N) = SM * 
     &             ( THETA_1_CR * ZETA_M(K1,UM,N) -
     &               THETA_1_CI * ZETA_M(K2,UM,N) )
                HMULT_1(K2,UM,N) = SM * 
     &             ( THETA_1_CR * ZETA_M(K2,UM,N) +
     &               THETA_1_CI * ZETA_M(K1,UM,N) )
              ENDIF
              HMULT_2(K1,UM,N) = SM * 
     &             ( THETA_2_CR * ZETA_P(K1,UM,N) -
     &               THETA_2_CI * ZETA_P(K2,UM,N) )
              HMULT_2(K2,UM,N) = SM * 
     &             ( THETA_2_CR * ZETA_P(K2,UM,N) +
     &               THETA_2_CI * ZETA_P(K1,UM,N) )
           ENDDO

C  Start loops over layers and user-streams

          ENDDO
        ENDIF
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE PARTLAYER_HMULT_UP ( N, UT )

C  Partial layer INTEGRATED homogeneous solution MULTIPLIERS (Upwelling)

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup/solution stuff (inputs to this module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'

C  include files of Multiplier variables (output to this module)

      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer indices

      INTEGER          N, UT

C  Local variables
C  ---------------

      INTEGER          UM, K, KO1, K0, K1, K2
      DOUBLE PRECISION UX_UP, SM, DX

C  real multipliers

      DOUBLE PRECISION ZDEL_R, ZX_UP_R, ZX_DN_R
      DOUBLE PRECISION THETA_DN_R, THETA_UP_R

C  complex multipliers

      DOUBLE PRECISION THETA_DN_CR, THETA_UP_CR
      DOUBLE PRECISION THETA_DN_CI, THETA_UP_CI
      DOUBLE PRECISION ZDEL_CR, ZX_UP_CR, ZX_DN_CR
      DOUBLE PRECISION ZDEL_CI, ZX_UP_CI, ZX_DN_CI

C  Partial layer multipliers

      DO UM = LOCAL_UM_START, N_USER_STREAMS

C  set-up

        UX_UP = T_UTUP_USERM(UT,UM)
        SM    = USER_SECANTS(UM)

C  real multipliers

        DO K = 1, K_REAL(N)
          ZDEL_R     = T_DELT_EIGEN(K,N)
          ZX_UP_R    = T_UTUP_EIGEN(K,UT)
          ZX_DN_R    = T_UTDN_EIGEN(K,UT)
          THETA_DN_R = ZX_DN_R - ZDEL_R * UX_UP
          THETA_UP_R = ZX_UP_R - UX_UP
          UT_HMULT_UD(K,UM,UT) = SM * THETA_DN_R * ZETA_P(K,UM,N)
          IF ( HSINGO(K,UM,N) )THEN
            DX = DELTAU_VERT(N) - PARTAU_VERT(UT)
            UT_HMULT_UU(K,UM,UT) = SM * UX_UP * DX
          ELSE
            UT_HMULT_UU(K,UM,UT) = SM * THETA_UP_R * ZETA_M(K,UM,N)
          ENDIF
        ENDDO

C  Complex multipliers

        KO1 = K_REAL(N) + 1

        DO K = 1, K_COMPLEX(N)

          K0 = 2*K - 2
          K1 = KO1 + K0
          K2 = K1  + 1

          ZDEL_CR     = T_DELT_EIGEN(K1,N)
          ZDEL_CI     = T_DELT_EIGEN(K2,N)
          ZX_UP_CR    = T_UTUP_EIGEN(K1,UT)
          ZX_UP_CI    = T_UTUP_EIGEN(K2,UT)
          ZX_DN_CR    = T_UTDN_EIGEN(K1,UT)
          ZX_DN_CI    = T_UTDN_EIGEN(K2,UT)

          THETA_DN_CR = ZX_DN_CR - ZDEL_CR * UX_UP
          THETA_DN_CI = ZX_DN_CI - ZDEL_CI * UX_UP
          THETA_UP_CR = ZX_UP_CR - UX_UP
          THETA_UP_CI = ZX_UP_CI

          UT_HMULT_UD(K1,UM,UT) = SM * 
     &             ( THETA_DN_CR * ZETA_P(K1,UM,N) -
     &               THETA_DN_CI * ZETA_P(K2,UM,N) )
          UT_HMULT_UD(K2,UM,UT) = SM * 
     &             ( THETA_DN_CR * ZETA_P(K2,UM,N) +
     &               THETA_DN_CI * ZETA_P(K1,UM,N) )
          IF ( HSINGO(K1,UM,N) )THEN
            DX = DELTAU_VERT(N) - PARTAU_VERT(UT)
            UT_HMULT_UU(K1,UM,UT) = SM * UX_UP * DX
            UT_HMULT_UU(K2,UM,UT) = ZERO
          ELSE
            UT_HMULT_UU(K1,UM,UT) = SM * 
     &             ( THETA_UP_CR * ZETA_M(K1,UM,N) - 
     &               THETA_UP_CI * ZETA_M(K2,UM,N) )
            UT_HMULT_UU(K2,UM,UT) = SM * 
     &             ( THETA_UP_CR * ZETA_M(K2,UM,N) + 
     &               THETA_UP_CI * ZETA_M(K1,UM,N) )
          ENDIF

        ENDDO
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE PARTLAYER_HMULT_DN ( N, UT )

C  Partial layer INTEGRATED homogeneous solution MULTIPLIERS (Downwelling)

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup/solution stuff (inputs to this module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'

C  include files of Multiplier variables (output to this module)

      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer indices

      INTEGER          N, UT

C  Local variables
C  ---------------

      INTEGER          UM, K, KO1, K0, K1, K2
      DOUBLE PRECISION UX_DN, SM, DX

C  real multipliers

      DOUBLE PRECISION ZDEL_R, ZX_UP_R, ZX_DN_R
      DOUBLE PRECISION THETA_DN_R, THETA_UP_R

C  complex multipliers

      DOUBLE PRECISION THETA_DN_CR, THETA_UP_CR
      DOUBLE PRECISION THETA_DN_CI, THETA_UP_CI
      DOUBLE PRECISION ZDEL_CR, ZX_UP_CR, ZX_DN_CR
      DOUBLE PRECISION ZDEL_CI, ZX_UP_CI, ZX_DN_CI

C  Partial layer multipliers

      DO UM = LOCAL_UM_START, N_USER_STREAMS

C  set-up

        UX_DN = T_UTDN_USERM(UT,UM)
        SM    = USER_SECANTS(UM)

C  real multipliers

        DO K = 1, K_REAL(N)
          ZDEL_R     = T_DELT_EIGEN(K,N)
          ZX_UP_R    = T_UTUP_EIGEN(K,UT)
          ZX_DN_R    = T_UTDN_EIGEN(K,UT)
          THETA_DN_R = ZX_DN_R - UX_DN
          THETA_UP_R = ZX_UP_R - ZDEL_R * UX_DN
          IF ( HSINGO(K,UM,N) ) THEN
            DX = DELTAU_VERT(N) - PARTAU_VERT(UT)
            UT_HMULT_DD(K,UM,UT) = SM * UX_DN * DX
          ELSE
            UT_HMULT_DD(K,UM,UT) = SM * THETA_DN_R * ZETA_M(K,UM,N)
          ENDIF
          UT_HMULT_DU(K,UM,UT) = SM * THETA_UP_R * ZETA_P(K,UM,N)
        ENDDO

C  Complex multipliers

        KO1 = K_REAL(N) + 1

        DO K = 1, K_COMPLEX(N)

          K0 = 2*K - 2
          K1 = KO1 + K0
          K2 = K1  + 1

          ZDEL_CR     = T_DELT_EIGEN(K1,N)
          ZDEL_CI     = T_DELT_EIGEN(K2,N)
          ZX_UP_CR    = T_UTUP_EIGEN(K1,UT)
          ZX_UP_CI    = T_UTUP_EIGEN(K2,UT)
          ZX_DN_CR    = T_UTDN_EIGEN(K1,UT)
          ZX_DN_CI    = T_UTDN_EIGEN(K2,UT)

          THETA_DN_CR = ZX_DN_CR - UX_DN
          THETA_DN_CI = ZX_DN_CI
          THETA_UP_CR = ZX_UP_CR - ZDEL_CR * UX_DN
          THETA_UP_CI = ZX_UP_CI - ZDEL_CI * UX_DN

          IF ( HSINGO(K1,UM,N) )THEN
            DX = DELTAU_VERT(N) - PARTAU_VERT(UT)
            UT_HMULT_DD(K1,UM,UT) = SM * UX_DN * DX
            UT_HMULT_DD(K2,UM,UT) = ZERO
          ELSE
            UT_HMULT_DD(K1,UM,UT) = SM * 
     &             ( THETA_DN_CR * ZETA_M(K1,UM,N) -
     &               THETA_DN_CI * ZETA_M(K2,UM,N) )
            UT_HMULT_DD(K2,UM,UT) = SM * 
     &             ( THETA_DN_CR * ZETA_M(K2,UM,N) +
     &               THETA_DN_CI * ZETA_M(K1,UM,N) )
          ENDIF
          UT_HMULT_DU(K1,UM,UT) = SM * 
     &             ( THETA_UP_CR * ZETA_P(K1,UM,N) - 
     &               THETA_UP_CI * ZETA_P(K2,UM,N) )
          UT_HMULT_DU(K2,UM,UT) = SM * 
     &             ( THETA_UP_CR * ZETA_P(K2,UM,N) + 
     &               THETA_UP_CI * ZETA_P(K1,UM,N) )

        ENDDO
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

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup and solution variables

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'

C  include file of Multiplier variables (output stored here)

      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

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

        DO UT = 1, N_PARTLAYERS
         N  = PARTLAYERS_LAYERIDX(UT)
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
C       Retaining only the first order term

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
                  SD = DELTAU_VERT(N) * UDEL             ! First order
! Second order    SD = SD*(ONE-HALF*DELTAU_VERT(N)*SIGMA_M(N,UM,IB))
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

        DO UT = 1, N_PARTLAYERS
         N  = PARTLAYERS_LAYERIDX(UT)
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
                SD = PARTAU_VERT(UT) * UX_DN
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
