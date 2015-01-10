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
C #            L_EMULT_MASTER (master), calling:                #
C #                L_WHOLELAYER_EMULT_UP                        #
C #                L_WHOLELAYER_EMULT_DN                        #
C #                L_PARTLAYER_EMULT_UP                         #
C #                L_PARTLAYER_EMULT_DN                         #
C #                                                             #
C ###############################################################

      SUBROUTINE L_HMULT_MASTER

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include file of setup/solution stuff (inputs to this module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'

C  include files of Multiplier variables (output to this module)

      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  Linearization files
C  -------------------

C  include file of input control variables for the linearization

      INCLUDE '../includes/VLIDORT_L_INPUTS.VARS'

C  include files of linearized setup and solution variables (input)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'

C  linearized multiplier (output goes in this file)

      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'

C  Local variables
C  ---------------

C  Control

      LOGICAL          DO_RTSOL_VARY ( MAXLAYERS )
      INTEGER          NPARAMS_VARY  ( MAXLAYERS )

C  integers

      INTEGER          N, UT, UTA, UM, K, Q, KO1, K0, K1, K2

C  real variables

      DOUBLE PRECISION UDEL, ZDEL_R, L_UDEL, L_ZDEL_R
      DOUBLE PRECISION L_T2, L_T1, HOM1, HOM2, SM
      DOUBLE PRECISION UX_UP, L_UX_UP(MAX_ATMOSWFS)
      DOUBLE PRECISION UX_DN, L_UX_DN(MAX_ATMOSWFS)
      DOUBLE PRECISION ZX_UP_R, ZX_DN_R
      DOUBLE PRECISION THETA_DN_R, THETA_UP_R
      DOUBLE PRECISION L_ZX_UP_R, L_ZX_DN_R
      DOUBLE PRECISION L_THETA_DN_R, L_THETA_UP_R

C  variables for the complex multiplier linearization

      DOUBLE PRECISION ZDEL_CR, ZDEL_CI
      DOUBLE PRECISION THETA_1_CR, L_THETA_1_CR
      DOUBLE PRECISION THETA_2_CR, L_THETA_2_CR
      DOUBLE PRECISION THETA_1_CI, L_THETA_1_CI
      DOUBLE PRECISION THETA_2_CI, L_THETA_2_CI
      DOUBLE PRECISION THETA_DN_CR, THETA_UP_CR
      DOUBLE PRECISION THETA_DN_CI, THETA_UP_CI
      DOUBLE PRECISION ZX_UP_CR, ZX_DN_CR
      DOUBLE PRECISION ZX_UP_CI, ZX_DN_CI
      DOUBLE PRECISION L_THETA_DN_CR, L_THETA_UP_CR
      DOUBLE PRECISION L_THETA_DN_CI, L_THETA_UP_CI
      DOUBLE PRECISION L_ZDEL_CR, L_ZX_UP_CR, L_ZX_DN_CR
      DOUBLE PRECISION L_ZDEL_CI, L_ZX_UP_CI, L_ZX_DN_CI

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

C  Start loops over layers and user-streams
C   Only done if layers are flagged

      DO N = 1, NLAYERS
       IF ( STERM_LAYERMASK_UP(N).OR.STERM_LAYERMASK_DN(N) ) THEN
        IF ( DO_RTSOL_VARY(N) ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS

C  setup

            UDEL = T_DELT_USERM(N,UM)
            SM   = USER_SECANTS(UM)

C  linearized Real multipliers
C    Chain Rule: similar to scalar code.

            DO K = 1, K_REAL(N)
              ZDEL_R    = T_DELT_EIGEN(K,N)
              DO Q = 1, NPARAMS_VARY(N)
                L_ZDEL_R = L_T_DELT_EIGEN(K,N,Q)
                L_UDEL   = L_T_DELT_USERM(N,UM,Q)
                IF ( HSINGO(K,UM,N) ) THEN
                  HOM1 = L_DELTAU_VERT(Q,N)*(ONE-SM) - 
     &                   HALF * L_KEIGEN(K,N,Q) * DELTAU_VERT(N)
                  L_HMULT_1(K,UM,N,Q) = HMULT_1(K,UM,N) * HOM1
                ELSE
                  L_T1 = L_ZDEL_R - L_UDEL
                  HOM1 =   L_KEIGEN(K,N,Q)*HMULT_1(K,UM,N) + SM*L_T1
                  L_HMULT_1(K,UM,N,Q) = ZETA_M(K,UM,N) * HOM1
                ENDIF
                L_T2 = - ZDEL_R * L_UDEL - L_ZDEL_R * UDEL
                HOM2 = - L_KEIGEN(K,N,Q)*HMULT_2(K,UM,N) + SM*L_T2
                L_HMULT_2(K,UM,N,Q) = ZETA_P(K,UM,N) * HOM2
              ENDDO
            ENDDO

C  Linearized Complex multipliers
C    Chain rule using real and complex parts.

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
              DO Q = 1, NPARAMS_VARY(N)
                L_ZDEL_CR = L_T_DELT_EIGEN(K1,N,Q)
                L_ZDEL_CI = L_T_DELT_EIGEN(K2,N,Q)
                L_UDEL    = L_T_DELT_USERM(N,UM,Q)
                L_THETA_2_CR = - L_ZDEL_CR * UDEL - ZDEL_CR * L_UDEL
                L_THETA_2_CI = - L_ZDEL_CI * UDEL - ZDEL_CI * L_UDEL
                L_THETA_1_CR = L_ZDEL_CR - L_UDEL
                L_THETA_1_CI = L_ZDEL_CI

                IF ( HSINGO(K1,UM,N) ) THEN
                  HOM1 = L_DELTAU_VERT(Q,N)*(ONE-SM) - 
     &                   HALF * L_KEIGEN(K1,N,Q) * DELTAU_VERT(N)
                  L_HMULT_1(K1,UM,N,Q) = HMULT_1(K1,UM,N) * HOM1
                  L_HMULT_1(K2,UM,N,Q) = ZERO
                ELSE
                 L_HMULT_1(K1,UM,N,Q) = SM * 
     &             (   THETA_1_CR * L_ZETA_M(K1,UM,N,Q) +
     &               L_THETA_1_CR *   ZETA_M(K1,UM,N)   -
     &                 THETA_1_CI * L_ZETA_M(K2,UM,N,Q) -
     &               L_THETA_1_CI *   ZETA_M(K2,UM,N)   )
                 L_HMULT_1(K2,UM,N,Q) = SM * 
     &             (   THETA_1_CR * L_ZETA_M(K2,UM,N,Q) +
     &               L_THETA_1_CR *   ZETA_M(K2,UM,N)   +
     &                 THETA_1_CI * L_ZETA_M(K1,UM,N,Q) +
     &               L_THETA_1_CI *   ZETA_M(K1,UM,N)   )
                ENDIF

                L_HMULT_2(K1,UM,N,Q) = SM * 
     &             (   THETA_2_CR * L_ZETA_P(K1,UM,N,Q) +
     &               L_THETA_2_CR *   ZETA_P(K1,UM,N)   -
     &                 THETA_2_CI * L_ZETA_P(K2,UM,N,Q) -
     &               L_THETA_2_CI *   ZETA_P(K2,UM,N)   )
                L_HMULT_2(K2,UM,N,Q) = SM * 
     &             (   THETA_2_CR * L_ZETA_P(K2,UM,N,Q) +
     &               L_THETA_2_CR *   ZETA_P(K2,UM,N)   +
     &                 THETA_2_CI * L_ZETA_P(K1,UM,N,Q) +
     &               L_THETA_2_CI *   ZETA_P(K1,UM,N)   )

              ENDDO
            ENDDO

C  End loops over user angles and layers

          ENDDO
        ENDIF 
       ENDIF
      ENDDO

C  partial layer multipliers
C  -------------------------

      DO UTA = 1, N_USER_LEVELS
       IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
        UT = PARTLAYERS_OUTINDEX(UTA)
        N  = PARTLAYERS_LAYERIDX(UT)

C  UPWELLING
C  ---------

        IF ( DO_UPWELLING.AND.DO_RTSOL_VARY(N) ) THEN

C  start code with loop over user angles

         DO UM = LOCAL_UM_START, N_USER_STREAMS

C  set-up

          UX_UP = T_UTUP_USERM(UT,UM)
          SM    = USER_SECANTS(UM)
          DO Q = 1, NPARAMS_VARY(N)
           L_UX_UP(Q) = L_T_UTUP_USERM(UT,UM,Q)
          ENDDO

C  Linearization of Real multipliers
C    Use Chain rule. Similar to scalar case.

          DO K = 1, K_REAL(N)
           ZDEL_R     = T_DELT_EIGEN(K,N)
           ZX_UP_R    = T_UTUP_EIGEN(K,UT)
           ZX_DN_R    = T_UTDN_EIGEN(K,UT)
           THETA_DN_R = ZX_DN_R - ZDEL_R * UX_UP
           THETA_UP_R = ZX_UP_R - UX_UP
           DO Q = 1, NPARAMS_VARY(N)
            L_ZDEL_R     = L_T_DELT_EIGEN(K,N,Q)
            L_ZX_UP_R    = L_T_UTUP_EIGEN(K,UT,Q)
            L_ZX_DN_R    = L_T_UTDN_EIGEN(K,UT,Q)
            L_THETA_DN_R = L_ZX_DN_R
     &             - L_ZDEL_R * UX_UP - ZDEL_R * L_UX_UP(Q)
            L_THETA_UP_R = L_ZX_UP_R - L_UX_UP(Q)
            L_UT_HMULT_UD(K,UM,UT,Q) = SM * 
     &       ( L_THETA_DN_R *   ZETA_P(K,UM,N) +
     &           THETA_DN_R * L_ZETA_P(K,UM,N,Q) )
            L_UT_HMULT_UU(K,UM,UT,Q) = SM * 
     &       ( L_THETA_UP_R *   ZETA_M(K,UM,N) +
     &           THETA_UP_R * L_ZETA_M(K,UM,N,Q) )
           ENDDO
          ENDDO

C  Linearization of Complex multipliers
C    Use Chain rule

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

           DO Q = 1, NPARAMS_VARY(N)
            L_ZDEL_CR     = L_T_DELT_EIGEN(K1,N,Q)
            L_ZDEL_CI     = L_T_DELT_EIGEN(K2,N,Q)
            L_ZX_UP_CR    = L_T_UTUP_EIGEN(K1,UT,Q)
            L_ZX_UP_CI    = L_T_UTUP_EIGEN(K2,UT,Q)
            L_ZX_DN_CR    = L_T_UTDN_EIGEN(K1,UT,Q)
            L_ZX_DN_CI    = L_T_UTDN_EIGEN(K2,UT,Q)

            L_THETA_DN_CR = L_ZX_DN_CR
     &             - L_ZDEL_CR * UX_UP - ZDEL_CR * L_UX_UP(Q)
            L_THETA_DN_CI = L_ZX_DN_CI
     &             - L_ZDEL_CI * UX_UP - ZDEL_CI * L_UX_UP(Q)
            L_THETA_UP_CR = L_ZX_UP_CR - L_UX_UP(Q)
            L_THETA_UP_CI = L_ZX_UP_CI
  
            L_UT_HMULT_UD(K1,UM,UT,Q) = SM * 
     &             ( L_THETA_DN_CR *   ZETA_P(K1,UM,N)   +
     &                 THETA_DN_CR * L_ZETA_P(K1,UM,N,Q) -
     &               L_THETA_DN_CI *   ZETA_P(K2,UM,N)   -
     &                 THETA_DN_CI * L_ZETA_P(K2,UM,N,Q) )

            L_UT_HMULT_UD(K2,UM,UT,Q) = SM * 
     &             ( L_THETA_DN_CR *   ZETA_P(K2,UM,N)   +
     &                 THETA_DN_CR * L_ZETA_P(K2,UM,N,Q) +
     &               L_THETA_DN_CI *   ZETA_P(K1,UM,N)   +
     &                 THETA_DN_CI * L_ZETA_P(K1,UM,N,Q) )

            L_UT_HMULT_UU(K1,UM,UT,Q) = SM * 
     &             ( L_THETA_UP_CR *   ZETA_M(K1,UM,N)   +
     &                 THETA_UP_CR * L_ZETA_M(K1,UM,N,Q) -
     &               L_THETA_UP_CI *   ZETA_M(K2,UM,N)   -
     &                 THETA_UP_CI * L_ZETA_M(K2,UM,N,Q) )

            L_UT_HMULT_UU(K2,UM,UT,Q) = SM * 
     &             ( L_THETA_UP_CR *   ZETA_M(K2,UM,N)   +
     &                 THETA_UP_CR * L_ZETA_M(K2,UM,N,Q) +
     &               L_THETA_UP_CI *   ZETA_M(K1,UM,N)   +
     &                 THETA_UP_CI * L_ZETA_M(K1,UM,N,Q) )

           ENDDO
          ENDDO

C  Finish loop over user angles

         ENDDO

C  End upwelling clause

        ENDIF

C  DOWNWELLING
C  -----------

        IF ( DO_DNWELLING.AND.DO_RTSOL_VARY(N) ) THEN

C  start code with loop over user angles

         DO UM = LOCAL_UM_START, N_USER_STREAMS

C  set-up

          UX_DN = T_UTDN_USERM(UT,UM)
          SM    = USER_SECANTS(UM)
          DO Q = 1, NPARAMS_VARY(N)
           L_UX_DN(Q) = L_T_UTDN_USERM(UT,UM,Q)
          ENDDO

C  Linearization of Real multipliers
C    Use Chain rule. Similar to scalar case.

          DO K = 1, K_REAL(N)
           ZDEL_R     = T_DELT_EIGEN(K,N)
           ZX_UP_R    = T_UTUP_EIGEN(K,UT)
           ZX_DN_R    = T_UTDN_EIGEN(K,UT)
           THETA_DN_R = ZX_DN_R - UX_DN
           THETA_UP_R = ZX_UP_R - ZDEL_R * UX_DN
           DO Q = 1, NPARAMS_VARY(N)
            L_ZDEL_R     = L_T_DELT_EIGEN(K,N,Q)
            L_ZX_UP_R    = L_T_UTUP_EIGEN(K,UT,Q)
            L_ZX_DN_R    = L_T_UTDN_EIGEN(K,UT,Q)
            L_THETA_UP_R = L_ZX_UP_R
     &             - L_ZDEL_R * UX_DN - ZDEL_R * L_UX_DN(Q)
            L_THETA_DN_R = L_ZX_DN_R - L_UX_DN(Q)
            L_UT_HMULT_DD(K,UM,UT,Q) = SM * 
     &       ( L_THETA_DN_R *   ZETA_M(K,UM,N) +
     &           THETA_DN_R * L_ZETA_M(K,UM,N,Q) )
            L_UT_HMULT_DU(K,UM,UT,Q) = SM * 
     &       ( L_THETA_UP_R *   ZETA_P(K,UM,N) +
     &           THETA_UP_R * L_ZETA_P(K,UM,N,Q) )
           ENDDO
          ENDDO

C  Linearization of Complex multipliers
C    Use Chain rule

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

           THETA_UP_CR = ZX_UP_CR - ZDEL_CR * UX_DN
           THETA_UP_CI = ZX_UP_CI - ZDEL_CI * UX_DN
           THETA_DN_CR = ZX_DN_CR - UX_DN
           THETA_DN_CI = ZX_DN_CI

           DO Q = 1, NPARAMS_VARY(N)
            L_ZDEL_CR     = L_T_DELT_EIGEN(K1,N,Q)
            L_ZDEL_CI     = L_T_DELT_EIGEN(K2,N,Q)
            L_ZX_UP_CR    = L_T_UTUP_EIGEN(K1,UT,Q)
            L_ZX_UP_CI    = L_T_UTUP_EIGEN(K2,UT,Q)
            L_ZX_DN_CR    = L_T_UTDN_EIGEN(K1,UT,Q)
            L_ZX_DN_CI    = L_T_UTDN_EIGEN(K2,UT,Q)

            L_THETA_UP_CR = L_ZX_UP_CR
     &             - L_ZDEL_CR * UX_DN - ZDEL_CR * L_UX_DN(Q)
            L_THETA_UP_CI = L_ZX_UP_CI
     &             - L_ZDEL_CI * UX_DN - ZDEL_CI * L_UX_DN(Q)
            L_THETA_DN_CR = L_ZX_DN_CR - L_UX_DN(Q)
            L_THETA_DN_CI = L_ZX_DN_CI
  
            L_UT_HMULT_DD(K1,UM,UT,Q) = SM * 
     &             ( L_THETA_DN_CR *   ZETA_M(K1,UM,N)   +
     &                 THETA_DN_CR * L_ZETA_M(K1,UM,N,Q) -
     &               L_THETA_DN_CI *   ZETA_M(K2,UM,N)   -
     &                 THETA_DN_CI * L_ZETA_M(K2,UM,N,Q) )

            L_UT_HMULT_DD(K2,UM,UT,Q) = SM * 
     &             ( L_THETA_DN_CR *   ZETA_M(K2,UM,N)   +
     &                 THETA_DN_CR * L_ZETA_M(K2,UM,N,Q) +
     &               L_THETA_DN_CI *   ZETA_M(K1,UM,N)   +
     &                 THETA_DN_CI * L_ZETA_M(K1,UM,N,Q) )

            L_UT_HMULT_DU(K1,UM,UT,Q) = SM * 
     &             ( L_THETA_UP_CR *   ZETA_P(K1,UM,N)   +
     &                 THETA_UP_CR * L_ZETA_P(K1,UM,N,Q) -
     &               L_THETA_UP_CI *   ZETA_P(K2,UM,N)   -
     &                 THETA_UP_CI * L_ZETA_P(K2,UM,N,Q) )

            L_UT_HMULT_DU(K2,UM,UT,Q) = SM * 
     &             ( L_THETA_UP_CR *   ZETA_P(K2,UM,N)   +
     &                 THETA_UP_CR * L_ZETA_P(K2,UM,N,Q) +
     &               L_THETA_UP_CI *   ZETA_P(K1,UM,N)   +
     &                 THETA_UP_CI * L_ZETA_P(K1,UM,N,Q) )

           ENDDO
          ENDDO

C  Finish loop over user angles

         ENDDO

C  end of downwelling clause

        ENDIF

C  End off-grid loop

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

      INCLUDE '../includes/VLIDORT.PARS'

C  standard files
C  --------------

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setups and multipliers (inputs to this module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  Linearization files
C  -------------------

C  include file of input control variables for the linearization

      INCLUDE '../includes/VLIDORT_L_INPUTS.VARS'

C  include files of linearized setup and solution variables (input)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'

C  linearized multiplier (output goes in this file)

      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'

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
C  Start loop over all varying layers K such that K </= N 

       DO UTA = 1, N_USER_LEVELS
        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
         UT = PARTLAYERS_OUTINDEX(UTA)
         N  = PARTLAYERS_LAYERIDX(UT)
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

       DO UTA = 1, N_USER_LEVELS
        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
         UT = PARTLAYERS_OUTINDEX(UTA)
         N  = PARTLAYERS_LAYERIDX(UT)
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

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setups and multiplier variables (input)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  include file of linearized setup variables (input)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'

C  linearized multiplier output goes in this file

      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'

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

C  End clause pseudo-spherical versus plane-parallel

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

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setups and multiplier variables (input)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  include file of linearized setup variables (input)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'

C  linearized multiplier output goes in this file

      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'

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

C  End clause pseudo-spherical versus plaen-parallel

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

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setups and multiplier variables (input)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  include file of linearized setup variables (input)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'

C  linearized multiplier output goes in this file

      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'

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

C  End clause pseudo-spherical versus plaen-parallel

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

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setups and multiplier variables (input)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  include file of linearized setup variables (input)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'

C  linearized multiplier output goes in this file

      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'

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
              V1 = ONE - PARTAU_VERT(UT) * USER_SECANTS(UM)
              V2 = - HALF * PARTAU_VERT(UT)
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
              V2 = - HALF * PARTAU_VERT(UT)
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
              V1 = ONE - PARTAU_VERT(UT) * USER_SECANTS(UM)
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

