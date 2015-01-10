! ###############################################################
! #                                                             #
! #                    THE VLIDORT  MODEL                       #
! #                                                             #
! #  Vectorized LInearized Discrete Ordinate Radiative Transfer #
! #  -          --         -        -        -         -        #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Author :      Robert. J. D. Spurr                          #
! #                                                             #
! #  Address :     RT Solutions, inc.                           #
! #                9 Channing Street                            #
! #                Cambridge, MA 02138, USA                     #
! #                Tel: (617) 492 1183                          #
! #                                                             #
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  Versions     :   2.0, 2.2, 2.3, 2.4, 2.4R, 2.4RT, 2.4RTC,  #
! #                   2.5, 2.6                                  #
! #  Release Date :   December 2005  (2.0)                      #
! #  Release Date :   March 2007     (2.2)                      #
! #  Release Date :   October 2007   (2.3)                      #
! #  Release Date :   December 2008  (2.4)                      #
! #  Release Date :   April 2009     (2.4R)                     #
! #  Release Date :   July 2009      (2.4RT)                    #
! #  Release Date :   October 2010   (2.4RTC)                   #
! #  Release Date :   March 2011     (2.5)                      #
! #  Release Date :   May 2012       (2.6)                      #
! #                                                             #
! #       NEW: TOTAL COLUMN JACOBIANS         (2.4)             #
! #       NEW: BPDF Land-surface KERNELS      (2.4R)            #
! #       NEW: Thermal Emission Treatment     (2.4RT)           #
! #       Consolidated BRDF treatment         (2.4RTC)          #
! #       f77/f90 Release                     (2.5)             #
! #       External SS / New I/O Structures    (2.6)             #
! #                                                             #
! ###############################################################

!    #####################################################
!    #                                                   #
!    #   This Version of LIDORT comes with a GNU-style   #
!    #   license. Please read the license carefully.     #
!    #                                                   #
!    #####################################################

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #            L_HMULT_MASTER (master)                          #
! #                                                             #
! #            L_EMULT_MASTER (master), calling:                #
! #                L_WHOLELAYER_EMULT_UP                        #
! #                L_WHOLELAYER_EMULT_DN                        #
! #                L_PARTLAYER_EMULT_UP                         #
! #                L_PARTLAYER_EMULT_DN                         #
! #                                                             #
! ###############################################################


      MODULE vlidort_l_multipliers

      PRIVATE
      PUBLIC :: L_HMULT_MASTER, &
                L_EMULT_MASTER

      CONTAINS

      SUBROUTINE L_HMULT_MASTER ( &
        DO_UPWELLING, DO_DNWELLING, &
        NLAYERS, N_USER_LEVELS, &
        N_USER_STREAMS, LOCAL_UM_START, &
        USER_SECANTS, &
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
        PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, &
        STERM_LAYERMASK_DN, &
        DELTAU_VERT, T_DELT_EIGEN, &
        T_UTUP_EIGEN, T_UTDN_EIGEN, &
        T_DELT_USERM, T_UTDN_USERM, &
        T_UTUP_USERM, &
        K_REAL, K_COMPLEX, HSINGO, &
        HMULT_1, HMULT_2, ZETA_M, ZETA_P, &
        DO_PROFILE_LINEARIZATION, DO_COLUMN_LINEARIZATION, &
        N_TOTALCOLUMN_WFS, LAYER_VARY_FLAG, &
        LAYER_VARY_NUMBER, &
        L_DELTAU_VERT, L_T_DELT_EIGEN, &
        L_T_UTUP_EIGEN, L_T_UTDN_EIGEN, &
        L_T_DELT_USERM, L_T_UTDN_USERM, &
        L_T_UTUP_USERM, L_KEIGEN, &
        L_ZETA_M, L_ZETA_P, &
        L_HMULT_1, L_HMULT_2, &
        L_UT_HMULT_UU, L_UT_HMULT_UD, &
        L_UT_HMULT_DU, L_UT_HMULT_DD )

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_DNWELLING
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      DOUBLE PRECISION, INTENT (IN) :: USER_SECANTS  ( MAX_USER_STREAMS )
      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_DN ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN &
          ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN &
          ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::          HSINGO &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_1 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_2 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: ZETA_M &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: ZETA_P &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      LOGICAL, INTENT (IN) ::          DO_PROFILE_LINEARIZATION
      LOGICAL, INTENT (IN) ::          DO_COLUMN_LINEARIZATION
      INTEGER, INTENT (IN) ::          N_TOTALCOLUMN_WFS
      LOGICAL, INTENT (IN) ::          LAYER_VARY_FLAG  ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          LAYER_VARY_NUMBER ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_EIGEN &
          ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTUP_EIGEN &
          ( MAXEVALUES, MAX_USER_LEVELS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_EIGEN &
          ( MAXEVALUES, MAX_USER_LEVELS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_USERM &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTUP_USERM &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_KEIGEN &
          ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_ZETA_M ( MAXEVALUES, &
            MAX_USER_STREAMS, MAXLAYERS,  MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_ZETA_P ( MAXEVALUES, &
            MAX_USER_STREAMS, MAXLAYERS,  MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (OUT) :: L_HMULT_1 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_HMULT_2 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_UT_HMULT_UU ( MAXEVALUES, &
            MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_UT_HMULT_UD ( MAXEVALUES, &
            MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_UT_HMULT_DU ( MAXEVALUES, &
            MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_UT_HMULT_DD ( MAXEVALUES, &
            MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Local variables
!  ---------------

!  Control

      LOGICAL ::          DO_RTSOL_VARY ( MAXLAYERS )
      INTEGER ::          NPARAMS_VARY  ( MAXLAYERS )

!  integers

      INTEGER ::          N, UT, UTA, UM, K, Q, KO1, K0, K1, K2

!  real variables

      DOUBLE PRECISION :: UDEL, ZDEL_R, L_UDEL, L_ZDEL_R
      DOUBLE PRECISION :: L_T2, L_T1, HOM1, HOM2, SM
      DOUBLE PRECISION :: UX_UP, L_UX_UP(MAX_ATMOSWFS)
      DOUBLE PRECISION :: UX_DN, L_UX_DN(MAX_ATMOSWFS)
      DOUBLE PRECISION :: ZX_UP_R, ZX_DN_R
      DOUBLE PRECISION :: THETA_DN_R, THETA_UP_R
      DOUBLE PRECISION :: L_ZX_UP_R, L_ZX_DN_R
      DOUBLE PRECISION :: L_THETA_DN_R, L_THETA_UP_R

!  variables for the complex multiplier linearization

      DOUBLE PRECISION :: ZDEL_CR, ZDEL_CI
      DOUBLE PRECISION :: THETA_1_CR, L_THETA_1_CR
      DOUBLE PRECISION :: THETA_2_CR, L_THETA_2_CR
      DOUBLE PRECISION :: THETA_1_CI, L_THETA_1_CI
      DOUBLE PRECISION :: THETA_2_CI, L_THETA_2_CI
      DOUBLE PRECISION :: THETA_DN_CR, THETA_UP_CR
      DOUBLE PRECISION :: THETA_DN_CI, THETA_UP_CI
      DOUBLE PRECISION :: ZX_UP_CR, ZX_DN_CR
      DOUBLE PRECISION :: ZX_UP_CI, ZX_DN_CI
      DOUBLE PRECISION :: L_THETA_DN_CR, L_THETA_UP_CR
      DOUBLE PRECISION :: L_THETA_DN_CI, L_THETA_UP_CI
      DOUBLE PRECISION :: L_ZDEL_CR, L_ZX_UP_CR, L_ZX_DN_CR
      DOUBLE PRECISION :: L_ZDEL_CI, L_ZX_UP_CI, L_ZX_DN_CI

!  Local control
!  -------------

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

!  whole layer multipliers
!  -----------------------

!  Start loops over layers and user-streams
!   Only done if layers are flagged

      DO N = 1, NLAYERS
       IF ( STERM_LAYERMASK_UP(N).OR.STERM_LAYERMASK_DN(N) ) THEN
        IF ( DO_RTSOL_VARY(N) ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS

!  setup

            UDEL = T_DELT_USERM(N,UM)
            SM   = USER_SECANTS(UM)

!  linearized Real multipliers
!    Chain Rule: similar to scalar code.

            DO K = 1, K_REAL(N)
              ZDEL_R    = T_DELT_EIGEN(K,N)
              DO Q = 1, NPARAMS_VARY(N)
                L_ZDEL_R = L_T_DELT_EIGEN(K,N,Q)
                L_UDEL   = L_T_DELT_USERM(N,UM,Q)
                IF ( HSINGO(K,UM,N) ) THEN
                  HOM1 = L_DELTAU_VERT(Q,N)*(ONE-SM) - &
                         HALF * L_KEIGEN(K,N,Q) * DELTAU_VERT(N)
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

!  Linearized Complex multipliers
!    Chain rule using real and complex parts.

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
                  HOM1 = L_DELTAU_VERT(Q,N)*(ONE-SM) - &
                         HALF * L_KEIGEN(K1,N,Q) * DELTAU_VERT(N)
                  L_HMULT_1(K1,UM,N,Q) = HMULT_1(K1,UM,N) * HOM1
                  L_HMULT_1(K2,UM,N,Q) = ZERO
                ELSE
                 L_HMULT_1(K1,UM,N,Q) = SM * &
                   (   THETA_1_CR * L_ZETA_M(K1,UM,N,Q) + &
                     L_THETA_1_CR *   ZETA_M(K1,UM,N)   - &
                       THETA_1_CI * L_ZETA_M(K2,UM,N,Q) - &
                     L_THETA_1_CI *   ZETA_M(K2,UM,N)   )
                 L_HMULT_1(K2,UM,N,Q) = SM * &
                   (   THETA_1_CR * L_ZETA_M(K2,UM,N,Q) + &
                     L_THETA_1_CR *   ZETA_M(K2,UM,N)   + &
                       THETA_1_CI * L_ZETA_M(K1,UM,N,Q) + &
                     L_THETA_1_CI *   ZETA_M(K1,UM,N)   )
                ENDIF

                L_HMULT_2(K1,UM,N,Q) = SM * &
                   (   THETA_2_CR * L_ZETA_P(K1,UM,N,Q) + &
                     L_THETA_2_CR *   ZETA_P(K1,UM,N)   - &
                       THETA_2_CI * L_ZETA_P(K2,UM,N,Q) - &
                     L_THETA_2_CI *   ZETA_P(K2,UM,N)   )
                L_HMULT_2(K2,UM,N,Q) = SM * &
                   (   THETA_2_CR * L_ZETA_P(K2,UM,N,Q) + &
                     L_THETA_2_CR *   ZETA_P(K2,UM,N)   + &
                       THETA_2_CI * L_ZETA_P(K1,UM,N,Q) + &
                     L_THETA_2_CI *   ZETA_P(K1,UM,N)   )

              ENDDO
            ENDDO

!  End loops over user angles and layers

          ENDDO
        ENDIF
       ENDIF
      ENDDO

!  partial layer multipliers
!  -------------------------

      DO UTA = 1, N_USER_LEVELS
       IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
        UT = PARTLAYERS_OUTINDEX(UTA)
        N  = PARTLAYERS_LAYERIDX(UT)

!  UPWELLING
!  ---------

        IF ( DO_UPWELLING.AND.DO_RTSOL_VARY(N) ) THEN

!  start code with loop over user angles

         DO UM = LOCAL_UM_START, N_USER_STREAMS

!  set-up

          UX_UP = T_UTUP_USERM(UT,UM)
          SM    = USER_SECANTS(UM)
          DO Q = 1, NPARAMS_VARY(N)
           L_UX_UP(Q) = L_T_UTUP_USERM(UT,UM,Q)
          ENDDO

!  Linearization of Real multipliers
!    Use Chain rule. Similar to scalar case.

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
            L_THETA_DN_R = L_ZX_DN_R &
                   - L_ZDEL_R * UX_UP - ZDEL_R * L_UX_UP(Q)
            L_THETA_UP_R = L_ZX_UP_R - L_UX_UP(Q)
            L_UT_HMULT_UD(K,UM,UT,Q) = SM * &
             ( L_THETA_DN_R *   ZETA_P(K,UM,N) + &
                 THETA_DN_R * L_ZETA_P(K,UM,N,Q) )
            L_UT_HMULT_UU(K,UM,UT,Q) = SM * &
             ( L_THETA_UP_R *   ZETA_M(K,UM,N) + &
                 THETA_UP_R * L_ZETA_M(K,UM,N,Q) )
           ENDDO
          ENDDO

!  Linearization of Complex multipliers
!    Use Chain rule

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

            L_THETA_DN_CR = L_ZX_DN_CR &
                   - L_ZDEL_CR * UX_UP - ZDEL_CR * L_UX_UP(Q)
            L_THETA_DN_CI = L_ZX_DN_CI &
                   - L_ZDEL_CI * UX_UP - ZDEL_CI * L_UX_UP(Q)
            L_THETA_UP_CR = L_ZX_UP_CR - L_UX_UP(Q)
            L_THETA_UP_CI = L_ZX_UP_CI

            L_UT_HMULT_UD(K1,UM,UT,Q) = SM * &
                   ( L_THETA_DN_CR *   ZETA_P(K1,UM,N)   + &
                       THETA_DN_CR * L_ZETA_P(K1,UM,N,Q) - &
                     L_THETA_DN_CI *   ZETA_P(K2,UM,N)   - &
                       THETA_DN_CI * L_ZETA_P(K2,UM,N,Q) )

            L_UT_HMULT_UD(K2,UM,UT,Q) = SM * &
                   ( L_THETA_DN_CR *   ZETA_P(K2,UM,N)   + &
                       THETA_DN_CR * L_ZETA_P(K2,UM,N,Q) + &
                     L_THETA_DN_CI *   ZETA_P(K1,UM,N)   + &
                       THETA_DN_CI * L_ZETA_P(K1,UM,N,Q) )

            L_UT_HMULT_UU(K1,UM,UT,Q) = SM * &
                   ( L_THETA_UP_CR *   ZETA_M(K1,UM,N)   + &
                       THETA_UP_CR * L_ZETA_M(K1,UM,N,Q) - &
                     L_THETA_UP_CI *   ZETA_M(K2,UM,N)   - &
                       THETA_UP_CI * L_ZETA_M(K2,UM,N,Q) )

            L_UT_HMULT_UU(K2,UM,UT,Q) = SM * &
                   ( L_THETA_UP_CR *   ZETA_M(K2,UM,N)   + &
                       THETA_UP_CR * L_ZETA_M(K2,UM,N,Q) + &
                     L_THETA_UP_CI *   ZETA_M(K1,UM,N)   + &
                       THETA_UP_CI * L_ZETA_M(K1,UM,N,Q) )

           ENDDO
          ENDDO

!  Finish loop over user angles

         ENDDO

!  End upwelling clause

        ENDIF

!  DOWNWELLING
!  -----------

        IF ( DO_DNWELLING.AND.DO_RTSOL_VARY(N) ) THEN

!  start code with loop over user angles

         DO UM = LOCAL_UM_START, N_USER_STREAMS

!  set-up

          UX_DN = T_UTDN_USERM(UT,UM)
          SM    = USER_SECANTS(UM)
          DO Q = 1, NPARAMS_VARY(N)
           L_UX_DN(Q) = L_T_UTDN_USERM(UT,UM,Q)
          ENDDO

!  Linearization of Real multipliers
!    Use Chain rule. Similar to scalar case.

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
            L_THETA_UP_R = L_ZX_UP_R &
                   - L_ZDEL_R * UX_DN - ZDEL_R * L_UX_DN(Q)
            L_THETA_DN_R = L_ZX_DN_R - L_UX_DN(Q)
            L_UT_HMULT_DD(K,UM,UT,Q) = SM * &
             ( L_THETA_DN_R *   ZETA_M(K,UM,N) + &
                 THETA_DN_R * L_ZETA_M(K,UM,N,Q) )
            L_UT_HMULT_DU(K,UM,UT,Q) = SM * &
             ( L_THETA_UP_R *   ZETA_P(K,UM,N) + &
                 THETA_UP_R * L_ZETA_P(K,UM,N,Q) )
           ENDDO
          ENDDO

!  Linearization of Complex multipliers
!    Use Chain rule

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

            L_THETA_UP_CR = L_ZX_UP_CR &
                   - L_ZDEL_CR * UX_DN - ZDEL_CR * L_UX_DN(Q)
            L_THETA_UP_CI = L_ZX_UP_CI &
                   - L_ZDEL_CI * UX_DN - ZDEL_CI * L_UX_DN(Q)
            L_THETA_DN_CR = L_ZX_DN_CR - L_UX_DN(Q)
            L_THETA_DN_CI = L_ZX_DN_CI

            L_UT_HMULT_DD(K1,UM,UT,Q) = SM * &
                   ( L_THETA_DN_CR *   ZETA_M(K1,UM,N)   + &
                       THETA_DN_CR * L_ZETA_M(K1,UM,N,Q) - &
                     L_THETA_DN_CI *   ZETA_M(K2,UM,N)   - &
                       THETA_DN_CI * L_ZETA_M(K2,UM,N,Q) )

            L_UT_HMULT_DD(K2,UM,UT,Q) = SM * &
                   ( L_THETA_DN_CR *   ZETA_M(K2,UM,N)   + &
                       THETA_DN_CR * L_ZETA_M(K2,UM,N,Q) + &
                     L_THETA_DN_CI *   ZETA_M(K1,UM,N)   + &
                       THETA_DN_CI * L_ZETA_M(K1,UM,N,Q) )

            L_UT_HMULT_DU(K1,UM,UT,Q) = SM * &
                   ( L_THETA_UP_CR *   ZETA_P(K1,UM,N)   + &
                       THETA_UP_CR * L_ZETA_P(K1,UM,N,Q) - &
                     L_THETA_UP_CI *   ZETA_P(K2,UM,N)   - &
                       THETA_UP_CI * L_ZETA_P(K2,UM,N,Q) )

            L_UT_HMULT_DU(K2,UM,UT,Q) = SM * &
                   ( L_THETA_UP_CR *   ZETA_P(K2,UM,N)   + &
                       THETA_UP_CR * L_ZETA_P(K2,UM,N,Q) + &
                     L_THETA_UP_CI *   ZETA_P(K1,UM,N)   + &
                       THETA_UP_CI * L_ZETA_P(K1,UM,N,Q) )

           ENDDO
          ENDDO

!  Finish loop over user angles

         ENDDO

!  end of downwelling clause

        ENDIF

!  End off-grid loop

       ENDIF
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE L_HMULT_MASTER

!

      SUBROUTINE L_EMULT_MASTER ( &
        DO_UPWELLING, DO_DNWELLING, &
        NLAYERS, N_USER_LEVELS, &
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
        PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, &
        STERM_LAYERMASK_DN, &
        EMULT_UP, EMULT_DN, &
        DO_PROFILE_LINEARIZATION, DO_COLUMN_LINEARIZATION, &
        N_TOTALCOLUMN_WFS, LAYER_VARY_FLAG, &
        LAYER_VARY_NUMBER, DO_PLANE_PARALLEL, &
        LAYER_PIS_CUTOFF, NBEAMS, &
        N_USER_STREAMS, &
        T_DELT_MUBAR, T_DELT_USERM, &
        ITRANS_USERM, SIGMA_P, &
        L_AVERAGE_SECANT, L_INITIAL_TRANS, &
        L_T_DELT_MUBAR, L_T_DELT_USERM, &
        USER_SECANTS, &
        DELTAU_VERT, INITIAL_TRANS, &
        EMULT_HOPRULE, SIGMA_M, &
        L_DELTAU_VERT, T_UTUP_USERM, &
        UT_EMULT_UP, &
        L_T_UTDN_MUBAR, L_T_UTUP_USERM, &
        PARTAU_VERT, UT_EMULT_DN, &
        L_T_UTDN_USERM, &
        L_EMULT_UP, L_EMULT_DN, &
        L_UT_EMULT_UP, L_UT_EMULT_DN )

!  Linearized multipliers for the Beam source terms

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_DNWELLING
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_DN ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: EMULT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: EMULT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      LOGICAL, INTENT (IN) ::          DO_PROFILE_LINEARIZATION
      LOGICAL, INTENT (IN) ::          DO_COLUMN_LINEARIZATION
      INTEGER, INTENT (IN) ::          N_TOTALCOLUMN_WFS
      LOGICAL, INTENT (IN) ::          LAYER_VARY_FLAG  ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          LAYER_VARY_NUMBER ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL
      INTEGER, INTENT (IN) ::          LAYER_PIS_CUTOFF ( MAXBEAMS )
      INTEGER, INTENT (IN) ::          NBEAMS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: ITRANS_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: SIGMA_P &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_AVERAGE_SECANT &
          ( MAXLAYERS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_INITIAL_TRANS &
          ( MAXLAYERS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_MUBAR &
          ( MAXLAYERS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: USER_SECANTS ( MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      LOGICAL, INTENT (IN) ::          EMULT_HOPRULE &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: SIGMA_M &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: UT_EMULT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_MUBAR &
          ( MAX_USER_LEVELS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTUP_USERM &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PARTAU_VERT ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_EMULT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_USERM &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (OUT) :: L_EMULT_UP ( MAX_USER_STREAMS, &
            MAXLAYERS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_EMULT_DN ( MAX_USER_STREAMS, &
                 MAXLAYERS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_UT_EMULT_UP ( MAX_USER_STREAMS, &
            MAX_PARTLAYERS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_UT_EMULT_DN ( MAX_USER_STREAMS, &
            MAX_PARTLAYERS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Local variables
!  ---------------

      INTEGER :: N, UT, UTA, K, K_PARAMETERS

!  Upwelling
!  =========

      IF ( DO_UPWELLING ) THEN

!  Whole layer upwelling
!  ---------------------

!  Loop over all  model  layers N
!    Profiles:  loop over all varying layers K such that K </= N
!    Columns :  K = 0

       DO N = 1, NLAYERS
        IF ( STERM_LAYERMASK_UP(N) ) THEN
         IF ( DO_PROFILE_LINEARIZATION ) THEN
          DO K = 1, NLAYERS
           IF ( N.GE.K ) THEN
            IF ( LAYER_VARY_FLAG(K) ) THEN
             K_PARAMETERS = LAYER_VARY_NUMBER(K)
             CALL L_WHOLELAYER_EMULT_UP ( &
               N, K, K_PARAMETERS, &
               DO_PLANE_PARALLEL, &
               LAYER_PIS_CUTOFF, NBEAMS, &
               N_USER_STREAMS, &
               T_DELT_MUBAR, T_DELT_USERM, &
               ITRANS_USERM, &
               EMULT_UP, SIGMA_P, &
               L_AVERAGE_SECANT, L_INITIAL_TRANS, &
               L_T_DELT_MUBAR, L_T_DELT_USERM, &
               L_EMULT_UP )

            ENDIF
           ENDIF
          ENDDO
         ELSE IF ( DO_COLUMN_LINEARIZATION ) THEN
          K = 0
          CALL L_WHOLELAYER_EMULT_UP ( &
            N, K, N_TOTALCOLUMN_WFS, &
            DO_PLANE_PARALLEL, &
            LAYER_PIS_CUTOFF, NBEAMS, &
            N_USER_STREAMS, &
            T_DELT_MUBAR, T_DELT_USERM, &
            ITRANS_USERM, &
            EMULT_UP, SIGMA_P, &
            L_AVERAGE_SECANT, L_INITIAL_TRANS, &
            L_T_DELT_MUBAR, L_T_DELT_USERM, &
            L_EMULT_UP )
         ENDIF
        ENDIF
       ENDDO

!  Partial layer upwelling
!  -----------------------

!  Start loop over all partial output UT occuring in layers N
!  Start loop over all varying layers K such that K </= N

       DO UTA = 1, N_USER_LEVELS
        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
         UT = PARTLAYERS_OUTINDEX(UTA)
         N  = PARTLAYERS_LAYERIDX(UT)
         IF ( STERM_LAYERMASK_UP(N) ) THEN
          IF ( DO_PROFILE_LINEARIZATION ) THEN
           DO K = 1, NLAYERS
            !L_UT_EMULT_UP(:,UT,K,:,:) = 0.0d0   ! ROBZERO 12/21
            IF ( N.GE.K ) THEN
             IF ( LAYER_VARY_FLAG(K) ) THEN
              K_PARAMETERS = LAYER_VARY_NUMBER(K)
              CALL L_PARTLAYER_EMULT_UP ( &
                N, UT, K, K_PARAMETERS, &
                DO_PLANE_PARALLEL, &
                LAYER_PIS_CUTOFF, NBEAMS, &
                N_USER_STREAMS, T_DELT_MUBAR, &
                T_UTUP_USERM, ITRANS_USERM, &
                UT_EMULT_UP, SIGMA_P, &
                L_AVERAGE_SECANT, L_INITIAL_TRANS, &
                L_T_DELT_MUBAR, L_T_UTDN_MUBAR, &
                L_T_UTUP_USERM, &
                L_UT_EMULT_UP )
             ENDIF
            ENDIF
           ENDDO
          ELSE IF ( DO_COLUMN_LINEARIZATION ) THEN
           K = 0
           CALL L_PARTLAYER_EMULT_UP ( &
             N, UT, K, N_TOTALCOLUMN_WFS, &
             DO_PLANE_PARALLEL, &
             LAYER_PIS_CUTOFF, NBEAMS, &
             N_USER_STREAMS, T_DELT_MUBAR, &
             T_UTUP_USERM, ITRANS_USERM, &
             UT_EMULT_UP, SIGMA_P, &
             L_AVERAGE_SECANT, L_INITIAL_TRANS, &
             L_T_DELT_MUBAR, L_T_UTDN_MUBAR, &
             L_T_UTUP_USERM, &
             L_UT_EMULT_UP )
          ENDIF
         ENDIF
        ENDIF
       ENDDO

!  end upwelling

      ENDIF

!  Downwelling
!  ===========

      IF ( DO_DNWELLING ) THEN

!  Whole layer downwelling
!  -----------------------

!  Start loop over all  model  layers N
!  Start loop over all varying layers K such that K </= N

       DO N = 1, NLAYERS
        IF ( STERM_LAYERMASK_DN(N) ) THEN
         IF ( DO_PROFILE_LINEARIZATION ) THEN
          DO K = 1, NLAYERS
           IF ( N.GE.K ) THEN
            IF ( LAYER_VARY_FLAG(K) ) THEN
             K_PARAMETERS = LAYER_VARY_NUMBER(K)
             CALL L_WHOLELAYER_EMULT_DN ( &
               N, K, K_PARAMETERS, &
               DO_PLANE_PARALLEL, &
               LAYER_PIS_CUTOFF, NBEAMS, &
               N_USER_STREAMS, USER_SECANTS, &
               DELTAU_VERT, INITIAL_TRANS, &
               ITRANS_USERM, &
               EMULT_DN, EMULT_HOPRULE, &
               SIGMA_M, &
               L_DELTAU_VERT, L_AVERAGE_SECANT, &
               L_INITIAL_TRANS, L_T_DELT_MUBAR, &
               L_T_DELT_USERM, &
               L_EMULT_DN )
            ENDIF
           ENDIF
          ENDDO
         ELSE IF  ( DO_COLUMN_LINEARIZATION ) THEN
          K = 0
          CALL L_WHOLELAYER_EMULT_DN ( &
            N, K, N_TOTALCOLUMN_WFS, &
            DO_PLANE_PARALLEL, &
            LAYER_PIS_CUTOFF, NBEAMS, &
            N_USER_STREAMS, USER_SECANTS, &
            DELTAU_VERT, INITIAL_TRANS, &
            ITRANS_USERM, &
            EMULT_DN, EMULT_HOPRULE, &
            SIGMA_M, &
            L_DELTAU_VERT, L_AVERAGE_SECANT, &
            L_INITIAL_TRANS, L_T_DELT_MUBAR, &
            L_T_DELT_USERM, &
            L_EMULT_DN )
         ENDIF
        ENDIF
       ENDDO

!  Partial layer downwelling
!  -------------------------

!  Start loop over all partial output UT occuring in layers N
!  Start loop over all varying layers K such that K </= N

       DO UTA = 1, N_USER_LEVELS
        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
         UT = PARTLAYERS_OUTINDEX(UTA)
         N  = PARTLAYERS_LAYERIDX(UT)
         IF ( STERM_LAYERMASK_DN(N) ) THEN
          IF ( DO_PROFILE_LINEARIZATION ) THEN
           DO K = 1, NLAYERS
            !L_UT_EMULT_DN(:,UT,K,:,:) = 0.0d0   ! ROBZERO 12/21
            IF ( N.GE.K ) THEN
             IF ( LAYER_VARY_FLAG(K) ) THEN
              K_PARAMETERS = LAYER_VARY_NUMBER(K)
              CALL L_PARTLAYER_EMULT_DN ( &
                N, UT, K, K_PARAMETERS, &
                DO_PLANE_PARALLEL, &
                LAYER_PIS_CUTOFF, NBEAMS, &
                N_USER_STREAMS, USER_SECANTS, &
                PARTAU_VERT, ITRANS_USERM, &
                UT_EMULT_DN, &
                EMULT_HOPRULE, SIGMA_M, &
                L_DELTAU_VERT, L_AVERAGE_SECANT, &
                L_INITIAL_TRANS, L_T_UTDN_MUBAR, &
                L_T_UTDN_USERM, &
                L_UT_EMULT_DN )
             ENDIF
            ENDIF
           ENDDO
          ELSE IF ( DO_COLUMN_LINEARIZATION ) THEN
           K = 0
           CALL L_PARTLAYER_EMULT_DN ( &
             N, UT, K, N_TOTALCOLUMN_WFS, &
             DO_PLANE_PARALLEL, &
             LAYER_PIS_CUTOFF, NBEAMS, &
             N_USER_STREAMS, USER_SECANTS, &
             PARTAU_VERT, ITRANS_USERM, &
             UT_EMULT_DN, &
             EMULT_HOPRULE, SIGMA_M, &
             L_DELTAU_VERT, L_AVERAGE_SECANT, &
             L_INITIAL_TRANS, L_T_UTDN_MUBAR, &
             L_T_UTDN_USERM, &
             L_UT_EMULT_DN )
          ENDIF
         ENDIF
        ENDIF
       ENDDO

!  end downwelling

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE L_EMULT_MASTER

!

      SUBROUTINE L_WHOLELAYER_EMULT_UP ( &
        N, K, K_PARAMETERS, &
        DO_PLANE_PARALLEL, &
        LAYER_PIS_CUTOFF, NBEAMS, &
        N_USER_STREAMS, &
        T_DELT_MUBAR, T_DELT_USERM, &
        ITRANS_USERM, &
        EMULT_UP, SIGMA_P, &
        L_AVERAGE_SECANT, L_INITIAL_TRANS, &
        L_T_DELT_MUBAR, L_T_DELT_USERM, &
        L_EMULT_UP )

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          N
      INTEGER, INTENT (IN) ::          K
      INTEGER, INTENT (IN) ::          K_PARAMETERS
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL
      INTEGER, INTENT (IN) ::          LAYER_PIS_CUTOFF ( MAXBEAMS )
      INTEGER, INTENT (IN) ::          NBEAMS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: ITRANS_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: EMULT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: SIGMA_P &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_AVERAGE_SECANT &
          ( MAXLAYERS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_INITIAL_TRANS &
          ( MAXLAYERS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_MUBAR &
          ( MAXLAYERS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (INOUT) :: L_EMULT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, 0:MAXLAYERS, &
            MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      DOUBLE PRECISION :: SU, V1, V2, WDEL, UDEL
      INTEGER          :: UM, Q, IB

!  Start Beam loop
!  ===============

      DO IB = 1, NBEAMS

!  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
         DO UM = 1, N_USER_STREAMS
           DO Q = 1, K_PARAMETERS
             L_EMULT_UP(UM,N,K,IB,Q) = ZERO
           ENDDO
         ENDDO
         GO TO 5678
       ENDIF

!  Profile linearizations: Two cases --------
!  (a) If N = K, multiplier for due to variations in the layer N
!  (b) If N > K, multiplier due to variations in a higher layer K
!  Column linearizations: One case ----------
!  (a) If K = 0, Multiplier for bulk (column) variations

!  transmittance factor

       WDEL = T_DELT_MUBAR(N,IB)

!  For the pseudo-spherical case
!  -----------------------------

       IF ( .NOT. DO_PLANE_PARALLEL ) THEN

!  Case(a)

        IF ( K.EQ.N .OR. K.EQ.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            UDEL = T_DELT_USERM(N,UM)
            SU = - ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
            DO Q = 1, K_PARAMETERS
              V1 = -L_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_P(N,UM,IB)
              IF ( K.EQ.0 ) V1 = V1 + L_INITIAL_TRANS (N,K,IB,Q)
              V2 = WDEL * L_T_DELT_USERM(N,UM,Q) + &
                   UDEL * L_T_DELT_MUBAR(N,K,IB,Q)
              L_EMULT_UP(UM,N,K,IB,Q) = EMULT_UP(UM,N,IB) * V1 + SU * V2
           ENDDO
          ENDDO
        ENDIF

!  Case (b)

        IF ( N.GT.K .AND. K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            UDEL = T_DELT_USERM(N,UM)
            SU = - ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
            DO Q = 1, K_PARAMETERS
              V1 = L_INITIAL_TRANS (N,K,IB,Q) - &
                 ( L_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_P(N,UM,IB) )
              V2 =  UDEL * L_T_DELT_MUBAR(N,K,IB,Q)
              L_EMULT_UP(UM,N,K,IB,Q) = EMULT_UP(UM,N,IB) * V1 + SU * V2
            ENDDO
          ENDDO
        ENDIF

!  For the plane-parallel case
!  ---------------------------

      ELSE

!  Case (a)

        IF ( K.EQ.N .OR. K.EQ.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            UDEL = T_DELT_USERM(N,UM)
            SU = - ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
            DO Q = 1, K_PARAMETERS
              V1 = ZERO
              IF ( K.EQ.0 ) V1 = L_INITIAL_TRANS (N,K,IB,Q)
              V2 = WDEL * L_T_DELT_USERM(N,UM,Q) + &
                   UDEL * L_T_DELT_MUBAR(N,K,IB,Q)
              L_EMULT_UP(UM,N,K,IB,Q) = EMULT_UP(UM,N,IB)*V1 + SU * V2
            ENDDO
          ENDDO
        ENDIF

!  Case (b)

        IF ( N.GT.K .AND. K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              V1 = L_INITIAL_TRANS(N,K,IB,Q)
              L_EMULT_UP(UM,N,K,IB,Q) = EMULT_UP(UM,N,IB) * V1
            ENDDO
          ENDDO
        ENDIF

!  End clause pseudo-spherical versus plane-parallel

       ENDIF

!  continuation point for next beam

 5678  CONTINUE

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE L_WHOLELAYER_EMULT_UP

!

      SUBROUTINE L_WHOLELAYER_EMULT_DN ( &
        N, K, K_PARAMETERS, &
        DO_PLANE_PARALLEL, &
        LAYER_PIS_CUTOFF, NBEAMS, &
        N_USER_STREAMS, USER_SECANTS, &
        DELTAU_VERT, INITIAL_TRANS, &
        ITRANS_USERM, &
        EMULT_DN, EMULT_HOPRULE, &
        SIGMA_M, &
        L_DELTAU_VERT, L_AVERAGE_SECANT, &
        L_INITIAL_TRANS, L_T_DELT_MUBAR, &
        L_T_DELT_USERM, &
        L_EMULT_DN )

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          N
      INTEGER, INTENT (IN) ::          K
      INTEGER, INTENT (IN) ::          K_PARAMETERS
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL
      INTEGER, INTENT (IN) ::          LAYER_PIS_CUTOFF ( MAXBEAMS )
      INTEGER, INTENT (IN) ::          NBEAMS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      DOUBLE PRECISION, INTENT (IN) :: USER_SECANTS  ( MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: ITRANS_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: EMULT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      LOGICAL, INTENT (IN) ::          EMULT_HOPRULE &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: SIGMA_M &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_AVERAGE_SECANT &
          ( MAXLAYERS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_INITIAL_TRANS &
          ( MAXLAYERS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_MUBAR &
          ( MAXLAYERS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (INOUT) :: L_EMULT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, 0:MAXLAYERS, &
            MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      DOUBLE PRECISION :: SD, V1, V2, V3
      INTEGER          :: UM, Q, IB

!  Start Beam loop
!  ===============

      DO IB = 1, NBEAMS

!  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
         DO UM = 1, N_USER_STREAMS
           DO Q = 1, K_PARAMETERS
             L_EMULT_DN(UM,N,K,IB,Q) = ZERO
           ENDDO
         ENDDO
         GO TO 5678
       ENDIF

!  Profile linearizations: Two cases --------
!  (a) If N = K, multiplier for due to variations in the layer N
!  (b) If N > K, multiplier due to variations in a higher layer K
!  Column linearizations: One case ----------
!  (a) If K = 0, Multiplier for bulk (column) variations

!  NOTE - use of L'Hopital's Rule is present in this module

!  For the pseudo-spherical case
!  -----------------------------

       IF ( .NOT. DO_PLANE_PARALLEL ) THEN

!  Case(a). Note the use of L'Hopital's Rule flag.

        IF ( K.EQ.N .OR. K.EQ.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
              V1 = ONE - DELTAU_VERT(N) * USER_SECANTS(UM)
              V2 = - HALF * DELTAU_VERT(N)
              DO Q = 1, K_PARAMETERS
                SD = V1 * L_DELTAU_VERT(Q,N) + &
                     V2 * L_AVERAGE_SECANT(N,K,IB,Q)
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

!  Case (b)

        IF ( N.GT.K .AND. K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
              V2 = - HALF * DELTAU_VERT(N)
              DO Q = 1, K_PARAMETERS
                SD =        L_INITIAL_TRANS (N,K,IB,Q) + &
                       V2 * L_AVERAGE_SECANT(N,K,IB,Q)
                L_EMULT_DN(UM,N,K,IB,Q) = EMULT_DN(UM,N,IB) * SD
              ENDDO
            ELSE
              SD = ITRANS_USERM(N,UM,IB) / SIGMA_M(N,UM,IB)
              DO Q = 1, K_PARAMETERS
                V1 =   L_INITIAL_TRANS(N,K,IB,Q) - &
                     ( L_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_M(N,UM,IB) )
                V2 = - L_T_DELT_MUBAR(N,K,IB,Q)
                L_EMULT_DN(UM,N,K,IB,Q) = EMULT_DN(UM,N,IB)*V1 + SD*V2
              ENDDO
            ENDIF
          ENDDO
        ENDIF

!  For the plane-parallel case
!  ---------------------------

       ELSE

!  Case (a)

        IF ( K.EQ.N .OR. K.EQ.0 ) THEN
          DO UM = 1, N_USER_STREAMS

            IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
              V1 = ONE - DELTAU_VERT(N) * USER_SECANTS(UM)
              DO Q = 1, K_PARAMETERS
                SD = V1 * L_DELTAU_VERT(Q,N)
                V2 = ZERO
                IF ( K.EQ.0.AND.INITIAL_TRANS(N,IB).NE.ZERO) &
                    V2 = V2 + L_INITIAL_TRANS (N,K,IB,Q)
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

!  Case (b)

        IF ( N.GT.K .AND. K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              V1 = L_INITIAL_TRANS(N,K,IB,Q)
              L_EMULT_DN(UM,N,K,IB,Q) = EMULT_DN(UM,N,IB) * V1
            ENDDO
          ENDDO
        ENDIF

!  End clause pseudo-spherical versus plaen-parallel

       ENDIF

!  continuation point for next beam

 5678  CONTINUE

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE L_WHOLELAYER_EMULT_DN

!

      SUBROUTINE L_PARTLAYER_EMULT_UP ( &
        N, UT, K, K_PARAMETERS, &
        DO_PLANE_PARALLEL, &
        LAYER_PIS_CUTOFF, NBEAMS, &
        N_USER_STREAMS, T_DELT_MUBAR, &
        T_UTUP_USERM, ITRANS_USERM, &
        UT_EMULT_UP, SIGMA_P, &
        L_AVERAGE_SECANT, L_INITIAL_TRANS, &
        L_T_DELT_MUBAR, L_T_UTDN_MUBAR, &
        L_T_UTUP_USERM, &
        L_UT_EMULT_UP )

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          N
      INTEGER, INTENT (IN) ::          UT
      INTEGER, INTENT (IN) ::          K
      INTEGER, INTENT (IN) ::          K_PARAMETERS
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL
      INTEGER, INTENT (IN) ::          LAYER_PIS_CUTOFF ( MAXBEAMS )
      INTEGER, INTENT (IN) ::          NBEAMS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: ITRANS_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: UT_EMULT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: SIGMA_P &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_AVERAGE_SECANT &
          ( MAXLAYERS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_INITIAL_TRANS &
          ( MAXLAYERS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_MUBAR &
          ( MAXLAYERS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_MUBAR &
          ( MAX_USER_LEVELS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTUP_USERM &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (INOUT) :: L_UT_EMULT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, 0:MAXLAYERS, &
            MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      DOUBLE PRECISION :: SU, V1, V2, WDEL, UX_UP
      INTEGER          :: UM, Q, IB

!  Start Beam loop
!  ===============

      DO IB = 1, NBEAMS

!  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
         DO UM = 1, N_USER_STREAMS
           DO Q = 1, K_PARAMETERS
             L_UT_EMULT_UP(UM,UT,K,IB,Q) = ZERO
           ENDDO
         ENDDO
         GO TO 5678
       ENDIF

!  transmittance factor

       WDEL = T_DELT_MUBAR(N,IB)

!  Profile linearizations: Two cases --------
!  (a) If N = K, multiplier for due to variations in the layer N
!  (b) If N > K, multiplier due to variations in a higher layer K
!  Column linearizations: One case ----------
!  (a) If K = 0, Multiplier for bulk (column) variations

!  For the pseudo-spherical case
!  -----------------------------

       IF ( .NOT. DO_PLANE_PARALLEL ) THEN

!  Case(a)

        IF ( K.EQ.N .OR. K.EQ.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            UX_UP = T_UTUP_USERM(UT,UM)
            SU = ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
            DO Q = 1, K_PARAMETERS
              V1 = -L_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_P(N,UM,IB)
              IF ( K.EQ.0 ) V1 = V1 + L_INITIAL_TRANS (N,K,IB,Q)
              V2 =           L_T_UTDN_MUBAR(UT,K,IB,Q) - &
                     UX_UP * L_T_DELT_MUBAR(N, K,IB,Q) - &
                     WDEL  * L_T_UTUP_USERM(UT,UM,Q)
              L_UT_EMULT_UP(UM,UT,K,IB,Q) =       SU * V2 + &
                                 UT_EMULT_UP(UM,UT,IB) * V1
            ENDDO
          ENDDO
        ENDIF

!  ..(b)

        IF ( N.GT.K .AND. K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            UX_UP = T_UTUP_USERM(UT,UM)
            SU = ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
            DO Q = 1, K_PARAMETERS
              V1 = L_INITIAL_TRANS(N,K,IB,Q) - &
                   ( L_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_P(N,UM,IB) )
              V2 =           L_T_UTDN_MUBAR(UT,K,IB,Q) - &
                     UX_UP * L_T_DELT_MUBAR( N,K,IB,Q)
              L_UT_EMULT_UP(UM,UT,K,IB,Q) =       SU * V2 + &
                                 UT_EMULT_UP(UM,UT,IB) * V1
            ENDDO
          ENDDO
        ENDIF

!  For the plane-parallel case
!  ---------------------------

       ELSE

!  Case (a)

        IF ( K.EQ.N .OR. K.EQ.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            UX_UP = T_UTUP_USERM(UT,UM)
            SU = ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
            DO Q = 1, K_PARAMETERS
              V1 = ZERO
              IF ( K.EQ.0 ) V1 = L_INITIAL_TRANS (N,K,IB,Q)
              V2 =           L_T_UTDN_MUBAR(UT,K,IB,Q) - &
                     UX_UP * L_T_DELT_MUBAR( N,K,IB,Q) - &
                     WDEL  * L_T_UTUP_USERM(UT,UM,Q)
              L_UT_EMULT_UP(UM,UT,K,IB,Q) =       SU * V2 + &
                                 UT_EMULT_UP(UM,UT,IB) * V1
            ENDDO
          ENDDO
        ENDIF

!  Case (b)

        IF ( N.GT.K .AND. K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              V1 = L_INITIAL_TRANS(N,K,IB,Q)
              L_UT_EMULT_UP(UM,UT,K,IB,Q) = UT_EMULT_UP(UM,UT,IB) * V1
            ENDDO
          ENDDO
        ENDIF

!  End clause pseudo-spherical versus plaen-parallel

       ENDIF

!  continuation point for next beam

 5678  CONTINUE

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE L_PARTLAYER_EMULT_UP

!

      SUBROUTINE L_PARTLAYER_EMULT_DN ( &
        N, UT, K, K_PARAMETERS, &
        DO_PLANE_PARALLEL, &
        LAYER_PIS_CUTOFF, NBEAMS, &
        N_USER_STREAMS, USER_SECANTS, &
        PARTAU_VERT, ITRANS_USERM, &
        UT_EMULT_DN, &
        EMULT_HOPRULE, SIGMA_M, &
        L_DELTAU_VERT, L_AVERAGE_SECANT, &
        L_INITIAL_TRANS, L_T_UTDN_MUBAR, &
        L_T_UTDN_USERM, &
        L_UT_EMULT_DN )

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          N
      INTEGER, INTENT (IN) ::          UT
      INTEGER, INTENT (IN) ::          K
      INTEGER, INTENT (IN) ::          K_PARAMETERS
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL
      INTEGER, INTENT (IN) ::          LAYER_PIS_CUTOFF ( MAXBEAMS )
      INTEGER, INTENT (IN) ::          NBEAMS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      DOUBLE PRECISION, INTENT (IN) :: USER_SECANTS  ( MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: PARTAU_VERT ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: ITRANS_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: UT_EMULT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      LOGICAL, INTENT (IN) ::          EMULT_HOPRULE &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: SIGMA_M &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_AVERAGE_SECANT &
          ( MAXLAYERS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_INITIAL_TRANS &
          ( MAXLAYERS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_MUBAR &
          ( MAX_USER_LEVELS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_USERM &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (INOUT) :: L_UT_EMULT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, 0:MAXLAYERS, &
            MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      DOUBLE PRECISION :: SD, V1, V2, V3
      INTEGER          :: UM, Q, IB

!  Start Beam loop
!  ===============

      DO IB = 1, NBEAMS

!  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
         DO UM = 1, N_USER_STREAMS
           DO Q = 1, K_PARAMETERS
             L_UT_EMULT_DN(UM,UT,K,IB,Q) = ZERO
           ENDDO
         ENDDO
         GO TO 5678
       ENDIF

!  Profile linearizations: Two cases --------
!  (a) If N = K, multiplier for due to variations in the layer N
!  (b) If N > K, multiplier due to variations in a higher layer K
!  Column linearizations: One case ----------
!  (a) If K = 0, Multiplier for bulk (column) variations

!  NOTE - use of L'Hopital's Rule is present in this module

!  For the pseudo-spherical case
!  -----------------------------

       IF ( .NOT. DO_PLANE_PARALLEL ) THEN

!  Case(a)

        IF ( K.EQ.N .OR. K.EQ.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
              V1 = ONE - PARTAU_VERT(UT) * USER_SECANTS(UM)
              V2 = - HALF * PARTAU_VERT(UT)
              DO Q = 1, K_PARAMETERS
                SD = V1 * L_DELTAU_VERT(Q,N) + &
                     V2 * L_AVERAGE_SECANT(N,K,IB,Q)
                V3 = ZERO
                IF ( K.EQ.0 ) V3 = L_INITIAL_TRANS (N,K,IB,Q)
                L_UT_EMULT_DN(UM,UT,K,IB,Q) = &
                         UT_EMULT_DN(UM,UT,IB) * ( SD + V3 )
              ENDDO
            ELSE
              SD = ITRANS_USERM(N,UM,IB) / SIGMA_M(N,UM,IB)
              DO Q = 1, K_PARAMETERS
                V1 = - L_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_M(N,UM,IB)
                IF ( K.EQ.0 ) V1 = V1 + L_INITIAL_TRANS (N,K,IB,Q)
                V2 = L_T_UTDN_USERM(UT,UM,Q) - L_T_UTDN_MUBAR(UT,K,IB,Q)
                L_UT_EMULT_DN(UM,UT,K,IB,Q) =       SD * V2 + &
                                 UT_EMULT_DN(UM,UT,IB) * V1
              ENDDO
            ENDIF
          ENDDO
        ENDIF

!  Case (b), profile only

        IF ( N.GT.K .AND. K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
              V2 = - HALF * PARTAU_VERT(UT)
              DO Q = 1, K_PARAMETERS
                SD =        L_INITIAL_TRANS (N,K,IB,Q) + &
                       V2 * L_AVERAGE_SECANT(N,K,IB,Q)
                L_UT_EMULT_DN(UM,UT,K,IB,Q) = UT_EMULT_DN(UM,UT,IB)*SD
              ENDDO
            ELSE
              SD = ITRANS_USERM(N,UM,IB) / SIGMA_M(N,UM,IB)
              DO Q = 1, K_PARAMETERS
                V1 = L_INITIAL_TRANS(N,K,IB,Q) - &
                     ( L_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_M(N,UM,IB) )
                V2 = - L_T_UTDN_MUBAR(UT,K,IB,Q)
                L_UT_EMULT_DN(UM,UT,K,IB,Q) =       SD * V2 + &
                                 UT_EMULT_DN(UM,UT,IB) * V1
              ENDDO
            ENDIF
          ENDDO
        ENDIF

!  For the plane-parallel case
!  ---------------------------

      ELSE

!  Case (a)

        IF ( K.EQ.N .OR. K.EQ.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
              V1 = ONE - PARTAU_VERT(UT) * USER_SECANTS(UM)
              DO Q = 1, K_PARAMETERS
                V3 = ZERO
                IF ( K.EQ.0 ) V3 = L_INITIAL_TRANS (N,K,IB,Q)
                SD = V1 * L_DELTAU_VERT(Q,N)
                L_UT_EMULT_DN(UM,UT,K,IB,Q) = &
                      UT_EMULT_DN(UM,UT,IB)  * ( SD + V3 )
              ENDDO
            ELSE
              SD = ITRANS_USERM(N,UM,IB) / SIGMA_M(N,UM,IB)
              DO Q = 1, K_PARAMETERS
                V1 = ZERO
                IF ( K.EQ.0 ) V1 = L_INITIAL_TRANS (N,K,IB,Q)
                V2 = L_T_UTDN_USERM(UT,UM,Q) - L_T_UTDN_MUBAR(UT,K,IB,Q)
                L_UT_EMULT_DN(UM,UT,K,IB,Q) =       SD * V2 + &
                                 UT_EMULT_DN(UM,UT,IB) * V1
              ENDDO
            ENDIF
          ENDDO
        ENDIF

!   Case (b), profile only

        IF ( N.GT.K .AND. K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              V1 = L_INITIAL_TRANS(N,K,IB,Q)
              L_UT_EMULT_DN(UM,UT,K,IB,Q) = UT_EMULT_DN(UM,UT,IB) * V1
            ENDDO
          ENDDO
        ENDIF

!  End clause pseudo-spherical versus plaen-parallel

       ENDIF

!  continuation point for next beam

 5678  CONTINUE

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE L_PARTLAYER_EMULT_DN


      END MODULE vlidort_l_multipliers

