! ###############################################################
! #                                                             #
! #                    THE VECTOR LIDORT MODEL                  #
! #                                                             #
! #  (Vector LInearized Discrete Ordinate Radiative Transfer)   #
! #   -      --         -        -        -         -           #
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
!    #   This Version of VLIDORT comes with a GNU-style  #
!    #   license. Please read the license carefully.     #
!    #                                                   #
!    #####################################################

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #   setups                                                    #
! #          THERMAL_SETUP_PLUS                                 #
! #                                                             #
! #   discrete ordinate particular integrals                    #
! #          THERMAL_CLSOLUTION_PLUS                            #
! #                                                             #
! #   postprocessing source terms                               #
! #          THERMAL_STERMS_UP_PLUS                             #
! #          THERMAL_STERMS_DN_PLUS                             #
! #                                                             #
! #   LTE Linearization   (new 9/11/09)                         #
! #          THERMAL_LTE_LINEARIZATION                          #
! #                                                             #
! ###############################################################


      MODULE vlidort_l_thermalsup

      PRIVATE
      PUBLIC :: THERMAL_SETUP_PLUS, &
                THERMAL_CLSOLUTION_PLUS, &
                THERMAL_STERMS_UP_PLUS, &
                THERMAL_STERMS_DN_PLUS, &
                THERMAL_LTE_LINEARIZATION

      CONTAINS

      SUBROUTINE THERMAL_SETUP_PLUS ( &
        DO_UPWELLING, DO_DNWELLING, DO_MSMODE_THERMAL, &
        NLAYERS, N_THERMAL_COEFFS, &
        THERMAL_BB_INPUT, DO_THERMAL_TRANSONLY, &
        N_USER_STREAMS, DO_USER_STREAMS, &
        USER_STREAMS, DO_PARTLAYERS, &
        N_PARTLAYERS, PARTLAYERS_LAYERIDX, &
        STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, &
        OMEGA_TOTAL, DELTAU_VERT, &
        PARTAU_VERT, T_DELT_USERM, &
        T_UTDN_USERM, T_UTUP_USERM, &
        DO_ATMOS_LINEARIZATION, LAYER_VARY_FLAG, &
        LAYER_VARY_NUMBER, &
        L_OMEGA_TOTAL, L_DELTAU_VERT, &
        L_T_DELT_USERM, L_T_UTDN_USERM, &
        L_T_UTUP_USERM, &
        THERMCOEFFS, DELTAU_POWER, XTAU_POWER, &
        TCOM1, T_DIRECT_UP, &
        T_DIRECT_DN, T_UT_DIRECT_UP, &
        T_UT_DIRECT_DN, &
        L_THERMCOEFFS, L_DELTAU_POWER, &
        L_XTAU_POWER, L_TCOM1, &
        L_T_DIRECT_UP, L_T_DIRECT_DN, &
        L_T_UT_DIRECT_UP, L_T_UT_DIRECT_DN )

!  SET-UP OF THERMAL EXPANSION COEFFICIENTS, ALWAYS DONE AFTER DELTA-M.
!  SET-UP OF LINEARIZATION OF THERMAL EXPANSION COEFFICIENTS.

!  LINEARIZATION W.R.T PROFILE VARIABLES, NOT BB INPUT.

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_DNWELLING
      LOGICAL, INTENT (IN) ::          DO_MSMODE_THERMAL
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_THERMAL_COEFFS
      DOUBLE PRECISION, INTENT (IN) :: THERMAL_BB_INPUT ( 0:MAXLAYERS )
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      LOGICAL, INTENT (IN) ::          DO_USER_STREAMS
      DOUBLE PRECISION, INTENT (IN) :: USER_STREAMS  ( MAX_USER_STREAMS )
      LOGICAL, INTENT (IN) ::          DO_PARTLAYERS
      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_DN ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: OMEGA_TOTAL ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PARTAU_VERT ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      LOGICAL, INTENT (IN) ::          DO_ATMOS_LINEARIZATION
      LOGICAL, INTENT (IN) ::          LAYER_VARY_FLAG  ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          LAYER_VARY_NUMBER ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_TOTAL ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_USERM &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTUP_USERM &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (OUT) :: THERMCOEFFS &
          ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (OUT) :: DELTAU_POWER &
          ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (OUT) :: XTAU_POWER &
          ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (OUT) :: TCOM1 ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (OUT) :: T_DIRECT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: T_DIRECT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: T_UT_DIRECT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: T_UT_DIRECT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS )

      DOUBLE PRECISION, INTENT (OUT) :: L_THERMCOEFFS &
          ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_DELTAU_POWER &
          ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_XTAU_POWER &
          ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_TCOM1 &
          ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_T_DIRECT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_T_DIRECT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_T_UT_DIRECT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS)
      DOUBLE PRECISION, INTENT (OUT) :: L_T_UT_DIRECT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  LOCAL VARIABLES
!  ---------------

      INTEGER          :: N, S, UT, UM, Q, NPARS, NT
      DOUBLE PRECISION :: XTAU, SUM, COSMUM, L_D, L_X
      DOUBLE PRECISION :: OMEGAS1 ( MAXLAYERS )
      DOUBLE PRECISION :: HELP(MAXLAYERS)
      DOUBLE PRECISION :: T_MULT_UP(MAXLAYERS,0:MAX_THERMAL_COEFFS)
      DOUBLE PRECISION :: T_MULT_DN(MAXLAYERS,0:MAX_THERMAL_COEFFS)
      DOUBLE PRECISION :: L_T_MULT_UP &
          (MAXLAYERS,0:MAX_THERMAL_COEFFS,MAX_ATMOSWFS)
      DOUBLE PRECISION :: L_T_MULT_DN &
          (MAXLAYERS,0:MAX_THERMAL_COEFFS,MAX_ATMOSWFS)

!  POWERS OF OPTICAL THICKNESS
!  ---------------------------

!  WHOLE LAYER

      DO N = 1, NLAYERS
       DELTAU_POWER(N,1) = ONE
       DO S = 2, N_THERMAL_COEFFS
        DELTAU_POWER(N,S) = DELTAU_VERT(N) * DELTAU_POWER(N,S-1)
       END DO
      ENDDO

      IF ( DO_ATMOS_LINEARIZATION ) THEN
       DO N = 1, NLAYERS
        DO Q = 1, LAYER_VARY_NUMBER(N)
         L_DELTAU_POWER(N,1,Q) = ZERO
         L_D = L_DELTAU_VERT(Q,N) * DELTAU_VERT(N)
         DO S = 2, N_THERMAL_COEFFS
          L_DELTAU_POWER(N,S,Q) = DBLE(S-1) * L_D * DELTAU_POWER(N,S-1)
         ENDDO
        ENDDO
       ENDDO
      ENDIF

!  PARTIAL LAYER

      IF ( DO_PARTLAYERS ) THEN
       DO UT = 1, N_PARTLAYERS

        XTAU = PARTAU_VERT(UT)
        XTAU_POWER(UT,1) = ONE
        DO S = 2, N_THERMAL_COEFFS
         XTAU_POWER(UT,S) = XTAU * XTAU_POWER(UT,S-1)
        END DO

        IF ( DO_ATMOS_LINEARIZATION ) THEN
         N = PARTLAYERS_LAYERIDX(UT)
         DO Q = 1, LAYER_VARY_NUMBER(N)
          L_X = L_DELTAU_VERT(Q,N) *  PARTAU_VERT(UT)
          L_XTAU_POWER(UT,1,Q) = ZERO
          DO S = 2, N_THERMAL_COEFFS
           L_XTAU_POWER(UT,S,Q) = DBLE(S-1) * L_X * XTAU_POWER(UT,S-1)
          ENDDO
         ENDDO
        ENDIF

       ENDDO
      ENDIF

!  INITIAL SET OF COEFFICIENTS

      DO N = 1, NLAYERS
       THERMCOEFFS(N,1) = THERMAL_BB_INPUT(N-1)
       HELP(N) = &
          (THERMAL_BB_INPUT(N)-THERMAL_BB_INPUT(N-1)) / DELTAU_VERT(N)
      END DO

!  PIECEWISE CONTINUOUS FOR LINEAR REGIME

      IF ( N_THERMAL_COEFFS == 2 ) THEN
       DO N = 1, NLAYERS
        THERMCOEFFS(N,2) = HELP(N)
       END DO
      END IF

      IF ( N_THERMAL_COEFFS .NE. 2 ) STOP' NUMBER OF COEFFICIENTS NOT 2'

!  AUXILIARY QUANTITIES

      IF ( DO_THERMAL_TRANSONLY ) THEN
        DO N = 1, NLAYERS
          DO S = 1, N_THERMAL_COEFFS
           TCOM1(N,S) = THERMCOEFFS(N,S)
          END DO
        END DO
      ELSE
        DO N = 1, NLAYERS
          OMEGAS1(N) = ONE - OMEGA_TOTAL(N)
          DO S = 1, N_THERMAL_COEFFS
           TCOM1(N,S) = THERMCOEFFS(N,S) * OMEGAS1(N)
          END DO
        END DO
      ENDIF

! LINEARIZED COEFFICIENTS AND AUXILIARY QUANTITIES

      IF ( DO_THERMAL_TRANSONLY ) THEN
       DO N = 1, NLAYERS
        IF ( LAYER_VARY_FLAG(N) ) THEN
         DO Q = 1, LAYER_VARY_NUMBER(N)
          L_THERMCOEFFS(N,1,Q) = ZERO
          L_THERMCOEFFS(N,2,Q) = - L_DELTAU_VERT(Q,N)*THERMCOEFFS(N,2)
          DO S = 1, N_THERMAL_COEFFS
           L_TCOM1(N,S,Q) =  L_THERMCOEFFS(N,S,Q)
          ENDDO
         ENDDO
        ELSE
         DO S = 1, N_THERMAL_COEFFS
           L_TCOM1(N,S,1) = ZERO
         ENDDO
        ENDIF
       ENDDO
      ELSE
       DO N = 1, NLAYERS
        IF ( LAYER_VARY_FLAG(N) ) THEN
         DO Q = 1, LAYER_VARY_NUMBER(N)
          L_THERMCOEFFS(N,1,Q) = ZERO
          L_THERMCOEFFS(N,2,Q) = - L_DELTAU_VERT(Q,N)*THERMCOEFFS(N,2)
          DO S = 1, N_THERMAL_COEFFS
           L_TCOM1(N,S,Q) =  OMEGAS1(N) * L_THERMCOEFFS(N,S,Q) &
             - THERMCOEFFS(N,S) * L_OMEGA_TOTAL(Q,N) * OMEGA_TOTAL(N)
          ENDDO
         ENDDO
        ELSE
         DO S = 1, N_THERMAL_COEFFS
           L_TCOM1(N,S,1) = ZERO
         ENDDO
        ENDIF
       ENDDO
      ENDIF

! RETURN IF POST-PROCESSING NOT FLAGGED

      IF ( .NOT. DO_USER_STREAMS ) RETURN

!  ZERO DIRECT SOLUTIONS if working in MSMODE only, then return

      IF ( DO_MSMODE_THERMAL ) THEN
         T_DIRECT_UP = 0.0d0   ; T_UT_DIRECT_UP = 0.0d0
         T_DIRECT_DN = 0.0d0   ; T_UT_DIRECT_DN = 0.0d0
         L_T_DIRECT_UP = 0.0d0 ; L_T_UT_DIRECT_UP = 0.0d0
         L_T_DIRECT_DN = 0.0d0 ; L_T_UT_DIRECT_DN = 0.0d0
         RETURN
      ENDIF

!  SHORT HAND

      NT = N_THERMAL_COEFFS

!  UPWELLING DIRECT SOLUTION SOURCE TERMS
!  --------------------------------------

      IF ( DO_UPWELLING ) THEN

!  START USER STREAM LOOP

       DO UM = 1, N_USER_STREAMS
        COSMUM = USER_STREAMS(UM)

!  DIRECT SOLUTION: WHOLE LAYER SOURCE TERMS
!   NOTE: T_DELT_USERM(N,UM) WAS INDEXED OPPOSITELY

!  START LAYER LOOP

        DO N = 1, NLAYERS

         IF ( STERM_LAYERMASK_UP(N) ) THEN
          T_MULT_UP(N,NT) = TCOM1(N,NT)
          DO S = NT - 1, 1, -1
           T_MULT_UP(N,S) = TCOM1(N,S) + S * COSMUM * T_MULT_UP(N,S+1)
          ENDDO
          SUM = T_MULT_UP(N,1)
          DO S = 2, NT
           SUM = SUM + T_MULT_UP(N,S) * DELTAU_POWER(N,S)
          ENDDO
          T_MULT_UP(N,0) = - SUM
          T_DIRECT_UP(UM,N) = T_MULT_UP(N,0) * T_DELT_USERM(N,UM) &
                            + T_MULT_UP(N,1)
         ENDIF

!  LINEARIZATION OF DIRECT SOLUTION: WHOLE LAYER SOURCE TERMS

         IF ( DO_ATMOS_LINEARIZATION ) THEN
          IF ( STERM_LAYERMASK_UP(N).AND.LAYER_VARY_FLAG(N) ) THEN
           NPARS = LAYER_VARY_NUMBER(N)
           DO Q = 1, NPARS
            L_T_MULT_UP(N,NT,Q) = L_TCOM1(N,NT,Q)
            DO S = N_THERMAL_COEFFS - 1, 1, -1
             L_T_MULT_UP(N,S,Q) = L_TCOM1(N,S,Q) &
                        + S * COSMUM * L_T_MULT_UP(N,S+1,Q)
            END DO
            SUM = L_T_MULT_UP(N,1,Q)
            DO S = 2, N_THERMAL_COEFFS
             SUM = SUM + L_T_MULT_UP(N,S,Q) *   DELTAU_POWER(N,S) &
                       +   T_MULT_UP(N,S)   * L_DELTAU_POWER(N,S,Q)
            END DO
            L_T_MULT_UP(N,0,Q) = - SUM
            L_T_DIRECT_UP(UM,N,Q) = L_T_MULT_UP(N,1,Q) &
                  + L_T_MULT_UP(N,0,Q) *   T_DELT_USERM(N,UM) &
                  +   T_MULT_UP(N,0)   * L_T_DELT_USERM(N,UM,Q)
           END DO
          ENDIF
         ENDIF

!  END WHOLE LAYER LOOP

        ENDDO

!  DIRECT SOLUTION: PARTIAL LAYER SOURCE TERMS

        IF ( DO_PARTLAYERS ) THEN
         DO UT = 1, N_PARTLAYERS
          N  = PARTLAYERS_LAYERIDX(UT)
          T_UT_DIRECT_UP(UM,UT) = T_MULT_UP(N,0) * T_UTUP_USERM(UT,UM)
          SUM = ZERO
          DO S = 1, NT
           SUM = SUM + T_MULT_UP(N,S) * XTAU_POWER(UT,S)
          END DO
          T_UT_DIRECT_UP(UM,UT) = T_UT_DIRECT_UP(UM,UT) + SUM
         END DO
        ENDIF

!  DIRECT SOLUTION: LINEARIZED PARTIAL LAYER SOURCE TERMS

        IF ( DO_ATMOS_LINEARIZATION .AND. DO_PARTLAYERS ) THEN
         DO UT = 1, N_PARTLAYERS
          N  = PARTLAYERS_LAYERIDX(UT)
          IF ( LAYER_VARY_FLAG(N) ) THEN
           NPARS = LAYER_VARY_NUMBER(N)
           DO Q = 1, NPARS
            L_T_UT_DIRECT_UP(UM,UT,Q) = &
                  L_T_MULT_UP(N,0,Q) *   T_UTUP_USERM(UT,UM) &
                 +  T_MULT_UP(N,0)   * L_T_UTUP_USERM(UT,UM,Q)
            SUM = ZERO
            DO S = 1, N_THERMAL_COEFFS
             SUM = SUM + L_T_MULT_UP(N,S,Q) *   XTAU_POWER(UT,S) &
                       +   T_MULT_UP(N,S)   * L_XTAU_POWER(UT,S,Q)
            END DO
            L_T_UT_DIRECT_UP(UM,UT,Q) = L_T_UT_DIRECT_UP(UM,UT,Q) + SUM
           ENDDO
          END IF
         END DO
        ENDIF

!  END USER STREAM LOOP, AND UPWELLING

       ENDDO
      ENDIF

!  DOWNWELLING DIRECT SOLUTION SOURCE TERMS
!  ----------------------------------------

      IF ( DO_DNWELLING ) THEN

!  START USER STREAM LOOP

       DO UM = 1, N_USER_STREAMS
        COSMUM = USER_STREAMS(UM)

!  DIRECT SOLUTION: WHOLE LAYER SOURCE TERMS
!   NOTE: T_DELT_USERM(N,UM) WAS INDEXED OPPOSITELY

!  START WHOLE LAYER LOOP

        DO N = 1, NLAYERS

         IF ( STERM_LAYERMASK_DN(N) ) THEN
          T_MULT_DN(N,NT) = TCOM1(N,NT)
          DO S = NT - 1, 1, -1
           T_MULT_DN(N,S) = TCOM1(N,S) - S * COSMUM * T_MULT_DN(N,S+1)
          END DO
          T_MULT_DN(N,0) = - T_MULT_DN(N,1)
          T_DIRECT_DN(UM,N) = T_MULT_DN(N,0) * T_DELT_USERM(N,UM)
          SUM = ZERO
          DO S = 1, NT
           SUM = SUM + T_MULT_DN(N,S) * DELTAU_POWER(N,S)
          END DO
          T_DIRECT_DN(UM,N) = T_DIRECT_DN(UM,N) + SUM
         ENDIF

!  LINEARIZATION OF DIRECT SOLUTION: WHOLE LAYER SOURCE TERMS

         IF ( DO_ATMOS_LINEARIZATION ) THEN
          IF ( STERM_LAYERMASK_DN(N).AND.LAYER_VARY_FLAG(N) ) THEN
           NPARS = LAYER_VARY_NUMBER(N)
           DO Q = 1, NPARS
            L_T_MULT_DN(N,NT,Q) = L_TCOM1(N,NT,Q)
            DO S = NT - 1, 1, -1
             L_T_MULT_DN(N,S,Q) = L_TCOM1(N,S,Q) - &
                        S * COSMUM * L_T_MULT_DN(N,S+1,Q)
            ENDDO
            L_T_MULT_DN(N,0,Q) = - L_T_MULT_DN(N,1,Q)
            L_T_DIRECT_DN(UM,N,Q) = &
                 L_T_MULT_DN(N,0,Q) *   T_DELT_USERM(N,UM) &
                 + T_MULT_DN(N,0)   * L_T_DELT_USERM(N,UM,Q)
            SUM = ZERO
            DO S = 1, N_THERMAL_COEFFS
             SUM = SUM + L_T_MULT_DN(N,S,Q) *   DELTAU_POWER(N,S) &
                       +   T_MULT_DN(N,S)   * L_DELTAU_POWER(N,S,Q)
            END DO
            L_T_DIRECT_DN(UM,N,Q) = L_T_DIRECT_DN(UM,N,Q) + SUM
           ENDDO
          ENDIF
         ENDIF

!  END WHOLE LAYER LOOP

        ENDDO

!  DIRECT SOLUTION: PARTIAL LAYER SOURCE TERMS

        IF ( DO_PARTLAYERS ) THEN
         DO UT = 1, N_PARTLAYERS
          N  = PARTLAYERS_LAYERIDX(UT)
          T_UT_DIRECT_DN(UM,UT) = T_MULT_DN(N,0) * T_UTDN_USERM(UT,UM)
          SUM = ZERO
          DO S = 1, NT
           SUM = SUM + T_MULT_DN(N,S) * XTAU_POWER(UT,S)
          END DO
          T_UT_DIRECT_DN(UM,UT) = T_UT_DIRECT_DN(UM,UT) + SUM
         END DO
        ENDIF

!  DIRECT SOLUTION: LINEARIZED PARTIAL LAYER SOURCE TERMS

        IF ( DO_ATMOS_LINEARIZATION .AND. DO_PARTLAYERS ) THEN
         DO UT = 1, N_PARTLAYERS
          N  = PARTLAYERS_LAYERIDX(UT)
          IF ( LAYER_VARY_FLAG(N) ) THEN
           NPARS = LAYER_VARY_NUMBER(N)
           DO Q = 1, NPARS
            L_T_UT_DIRECT_DN(UM,UT,Q) = &
                   L_T_MULT_DN(N,0,Q) *   T_UTDN_USERM(UT,UM) &
                   + T_MULT_DN(N,0)   * L_T_UTDN_USERM(UT,UM,Q)
            SUM = ZERO
            DO S = 1, N_THERMAL_COEFFS
             SUM = SUM + L_T_MULT_DN(N,S,Q) *   XTAU_POWER(UT,S) &
                       +   T_MULT_DN(N,S)   * L_XTAU_POWER(UT,S,Q)
            END DO
            L_T_UT_DIRECT_DN(UM,UT,Q) = L_T_UT_DIRECT_DN(UM,UT,Q) + SUM
           ENDDO
          ENDIF
         ENDDO
        ENDIF

!  END USER STREAM LOOP, AND DOWNWELLING

       ENDDO
      ENDIF

!  FINISH

      RETURN
      END SUBROUTINE THERMAL_SETUP_PLUS

!

      SUBROUTINE THERMAL_CLSOLUTION_PLUS ( &
        DO_UPWELLING, DO_DNWELLING, &
        DO_QUAD_OUTPUT, DO_ADDITIONAL_MVOUT, &
        NSTREAMS, NLAYERS, &
        N_THERMAL_COEFFS, DO_THERMAL_TRANSONLY, &
        DO_ATMOS_LINEARIZATION, LAYER_VARY_FLAG, &
        LAYER_VARY_NUMBER, &
        QUAD_STREAMS, QUAD_HALFWTS, &
        NMOMENTS, NSTREAMS_2, &
        N_USER_STREAMS, LOCAL_UM_START, &
        N_PARTLAYERS, &
        PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, &
        STERM_LAYERMASK_DN, &
        OMEGA_GREEK, T_DELT_DISORDS, &
        T_DISORDS_UTUP, T_DISORDS_UTDN, &
        PI_XQP, PI_XUP, PI_XUM, PI_XQM_PRE, &
        L_OMEGA_GREEK, L_T_DELT_DISORDS, &
        L_T_DISORDS_UTDN, L_T_DISORDS_UTUP, &
        SAB, DAB, &
        L_SAB, L_DAB, &
        THERMCOEFFS, DELTAU_POWER, &
        XTAU_POWER, TCOM1, &
        L_THERMCOEFFS, L_DELTAU_POWER, &
        L_XTAU_POWER, L_TCOM1, &
        T_WUPPER, T_WLOWER, &
        UT_T_PARTIC, U_TPOS1, &
        U_TNEG1, U_TPOS2, &
        U_TNEG2, &
        L_T_WUPPER, L_T_WLOWER, &
        L_UT_T_PARTIC, L_U_TPOS1, &
        L_U_TNEG1, L_U_TPOS2, &
        L_U_TNEG2, &
        STATUS, MESSAGE, TRACE )

!  LAYER LINEARIZATIONS OF CLASSICAL THERMAL PARTICULAR INTEGRALS

      USE VLIDORT_PARS
      USE LAPACK_TOOLS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_DNWELLING
      LOGICAL, INTENT (IN) ::          DO_QUAD_OUTPUT
      LOGICAL, INTENT (IN) ::          DO_ADDITIONAL_MVOUT
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_THERMAL_COEFFS
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_ATMOS_LINEARIZATION
      LOGICAL, INTENT (IN) ::          LAYER_VARY_FLAG  ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          LAYER_VARY_NUMBER ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_HALFWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::          NMOMENTS
      INTEGER, INTENT (IN) ::          NSTREAMS_2
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_DN ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: OMEGA_GREEK &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTUP &
          ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTDN &
          ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PI_XQP &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_XUP &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_XUM &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_XQM_PRE &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_GREEK &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_DISORDS &
          ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DISORDS_UTDN &
          ( MAXSTREAMS, MAX_USER_LEVELS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DISORDS_UTUP &
          ( MAXSTREAMS, MAX_USER_LEVELS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: SAB &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DAB &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_SAB &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, &
            MAXSTOKES,  MAXLAYERS,  MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_DAB &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, &
            MAXSTOKES,  MAXLAYERS,  MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: THERMCOEFFS &
          ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_POWER &
          ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (IN) :: XTAU_POWER &
          ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (IN) :: TCOM1 ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (IN) :: L_THERMCOEFFS &
          ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_POWER &
          ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_XTAU_POWER &
          ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_TCOM1 &
          ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (OUT) ::  T_WUPPER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) ::  T_WLOWER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) ::  UT_T_PARTIC &
          ( MAXSTREAMS_2, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (OUT) ::  U_TPOS1 ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) ::  U_TNEG1 ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) ::  U_TPOS2 ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) ::  U_TNEG2 ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) ::  L_T_WUPPER &
          ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) ::  L_T_WLOWER &
          ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS)
      DOUBLE PRECISION, INTENT (OUT) ::  L_UT_T_PARTIC &
          ( MAXSTREAMS_2, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) ::  L_U_TPOS1 &
          ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) ::  L_U_TNEG1 &
          ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) ::  L_U_TPOS2 &
          ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) ::  L_U_TNEG2 &
          ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )

      INTEGER, INTENT (OUT) ::           STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

! LOCAL VARIABLES
! ---------------

      LOGICAL ::           DO_VARY
      INTEGER ::           AA, AA1, I, I1, J, J1, L, M, O1
      INTEGER ::           S, N, UT, NT, INFO, UM, Q, NVARY
      DOUBLE PRECISION ::  SUM_M, SUM_P, TK, K1, SD, SU, SUM
      DOUBLE PRECISION ::  POS1, POS2, NEG1, NEG2, SUM1, SUM2, A5
      DOUBLE PRECISION ::  TERM1, L_TERM1(MAX_ATMOSWFS), L_SUM1, L_SUM2
      DOUBLE PRECISION ::  TERM2, L_TERM2(MAX_ATMOSWFS)
      DOUBLE PRECISION ::  LSU, LSD, L_M, L_P, L_TK
      CHARACTER (LEN=3) :: CI, C3

!  LOCAL MATRICES AND VECTORS

      DOUBLE PRECISION :: HVEC1(MAXSTREAMS)
      DOUBLE PRECISION :: HVEC2(MAXSTREAMS)
      DOUBLE PRECISION :: JVEC1(MAXSTREAMS)
      DOUBLE PRECISION :: TVEC1(MAXSTREAMS_2,MAXLAYERS)
      DOUBLE PRECISION :: TVEC2(MAXSTREAMS_2,MAXLAYERS)
      DOUBLE PRECISION :: TMAT(MAXSTREAMS,MAXSTREAMS)
      INTEGER          :: TPIVOT(MAXSTREAMS)

      DOUBLE PRECISION :: T_HELP1(0:MAXMOMENTS)
      DOUBLE PRECISION :: T_HELP2(0:MAXMOMENTS)
      DOUBLE PRECISION :: L_T_HELP1(0:MAXMOMENTS,MAX_ATMOSWFS)
      DOUBLE PRECISION :: L_T_HELP2(0:MAXMOMENTS,MAX_ATMOSWFS)

      DOUBLE PRECISION :: LHVEC1(MAXSTREAMS,MAX_ATMOSWFS)
      DOUBLE PRECISION :: LHVEC2(MAXSTREAMS,MAX_ATMOSWFS)
      DOUBLE PRECISION :: LJVEC1(MAXSTREAMS,MAX_ATMOSWFS)
      DOUBLE PRECISION :: L_TVEC1(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      DOUBLE PRECISION :: L_TVEC2(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

      DOUBLE PRECISION :: &
         T_C_MINUS(MAXSTREAMS,MAXLAYERS,0:MAX_THERMAL_COEFFS), &
         T_C_PLUS (MAXSTREAMS,MAXLAYERS,0:MAX_THERMAL_COEFFS)

      DOUBLE PRECISION :: L_T_C_MINUS &
        ( MAXSTREAMS, MAXLAYERS, 0:MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_C_PLUS &
        ( MAXSTREAMS, MAXLAYERS, 0:MAX_THERMAL_COEFFS, MAX_ATMOSWFS )

! -------------------------------
!  ZERO THE BOUNDARY LAYER VALUES
! -------------------------------

      DO I = 1, NSTREAMS_2
       DO N = 1, NLAYERS
        T_WUPPER(I,N) = ZERO
        T_WLOWER(I,N) = ZERO
       ENDDO
      ENDDO

      IF ( DO_ATMOS_LINEARIZATION ) THEN
       DO I = 1, NSTREAMS_2
        DO N = 1, NLAYERS
         DO Q = 1, LAYER_VARY_NUMBER(N)
           L_T_WUPPER(I,N,Q) = ZERO
           L_T_WLOWER(I,N,Q) = ZERO
         ENDDO
        ENDDO
       ENDDO
      ENDIF

!  (1,1) COMPONENT ONLY

      O1 = 1

!  INITIALIZE EXCEPTION HANDLING

      STATUS  = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  SHORTHAND

      NT = N_THERMAL_COEFFS

!  FOURIER COMPONENT = 0 (ALWAYS)

      M = 0

!  THERMAL TRANSMITTANCE ONLY, QUADRATURE SOLUTIONS
!  ================================================

      IF ( DO_THERMAL_TRANSONLY ) THEN

!  WHOLE LAYER SOLUTIONS

       DO N = 1, NLAYERS
        DO AA = 1, NSTREAMS
         AA1 = AA + NSTREAMS
         K1 = QUAD_STREAMS(AA)
         TK = T_DELT_DISORDS(AA,N)
         T_C_MINUS(AA,N,NT)  = K1 * THERMCOEFFS(N,NT)
         T_C_PLUS(AA,N,NT)   = K1 * THERMCOEFFS(N,NT)
         DO S = NT - 1, 1, -1
          T_C_MINUS(AA,N,S)= K1*(THERMCOEFFS(N,S)-S*T_C_MINUS(AA,N,S+1))
          T_C_PLUS(AA,N,S) = K1*(THERMCOEFFS(N,S)+S*T_C_PLUS (AA,N,S+1))
         END DO
         SUM_P = T_C_PLUS (AA,N,1)
         SUM_M = T_C_MINUS(AA,N,1)
         DO S = 2, NT
          SUM_M = SUM_M + T_C_MINUS(AA,N,S) * DELTAU_POWER(N,S)
          SUM_P = SUM_P + T_C_PLUS(AA,N,S)  * DELTAU_POWER(N,S)
         END DO
         T_C_MINUS(AA,N,0) = - T_C_MINUS(AA,N,1)
         T_C_PLUS(AA,N,0)  = - SUM_P
         T_WLOWER(AA,N)  = TK * T_C_MINUS(AA,N,0) + SUM_M
         T_WUPPER(AA1,N) = TK * T_C_PLUS(AA,N,0)  + T_C_PLUS(AA,N,1)
        END DO

!  LINEARIZED TRANSMITTANCE-ONLY THERMAL SOLUTIONS

        IF ( DO_ATMOS_LINEARIZATION ) THEN
         DO Q = 1, LAYER_VARY_NUMBER(N)
          DO AA = 1, NSTREAMS
           AA1 = AA + NSTREAMS
           K1 = QUAD_STREAMS(AA)
           TK = T_DELT_DISORDS(AA,N)
           L_TK = L_T_DELT_DISORDS(AA,N,Q)
           L_T_C_MINUS(AA,N,NT,Q)  =  K1 * L_THERMCOEFFS(N,NT,Q)
           L_T_C_PLUS(AA,N,NT,Q)   = L_T_C_MINUS(AA,N,NT,Q)
           DO S = NT - 1, 1, -1
            L_T_C_MINUS(AA,N,S,Q) = &
             + K1 * (L_THERMCOEFFS(N,S,Q) - S * L_T_C_MINUS(AA,N,S+1,Q))
            L_T_C_PLUS(AA,N,S,Q)  = &
             + K1 * (L_THERMCOEFFS(N,S,Q) + S * L_T_C_PLUS(AA,N,S+1,Q))
           END DO
           L_P   = L_T_C_PLUS (AA,N,1,Q)
           L_M   = L_T_C_MINUS(AA,N,1,Q)
           DO S = 2, NT
            L_M = L_M + L_T_C_MINUS(AA,N,S,Q) *   DELTAU_POWER(N,S) &
                      +   T_C_MINUS(AA,N,S)   * L_DELTAU_POWER(N,S,Q)
            L_P = L_P + L_T_C_PLUS(AA,N,S,Q)  *   DELTAU_POWER(N,S) &
                      +   T_C_PLUS(AA,N,S)    * L_DELTAU_POWER(N,S,Q)
           END DO
           L_T_C_MINUS(AA,N,0,Q) = - L_T_C_MINUS(AA,N,1,Q)
           L_T_C_PLUS(AA,N,0,Q)  = - L_P
           L_T_WLOWER(AA,N,Q) = &
              (L_TK*T_C_MINUS(AA,N,0)+TK*L_T_C_MINUS(AA,N,0,Q)+L_M)
           L_T_WUPPER(AA1,N,Q) = &
                ( L_TK*T_C_PLUS(AA,N,0) + TK*L_T_C_PLUS(AA,N,0,Q) &
                                + L_T_C_PLUS(AA,N,1,Q) )
          END DO
         ENDDO
        ENDIF

!  END LAYER LOOP

       END DO

!  OFFGRID: ONLY FOR QUADRATURE OR MEAN-VALUE OUTPUT

       IF ( DO_QUAD_OUTPUT .OR. DO_ADDITIONAL_MVOUT ) THEN
        IF ( N_PARTLAYERS .GT. 0 ) THEN
         DO UT = 1, N_PARTLAYERS
          N  = PARTLAYERS_LAYERIDX(UT)

!  REGULAR

          DO AA = 1, NSTREAMS
           AA1 = AA + NSTREAMS
           SD = T_C_MINUS(AA,N,0) * T_DISORDS_UTDN(AA,UT)
           SU = T_C_PLUS(AA,N,0)  * T_DISORDS_UTUP(AA,UT)
           DO S = 1, NT
            SD = SD + T_C_MINUS(AA,N,S) * XTAU_POWER(UT,S)
            SU = SU + T_C_PLUS(AA,N,S)  * XTAU_POWER(UT,S)
           END DO
           UT_T_PARTIC(AA,UT)  = SD
           UT_T_PARTIC(AA1,UT) = SU
          END DO

!  LINEARIZED

          IF ( DO_ATMOS_LINEARIZATION ) THEN
           IF ( LAYER_VARY_FLAG(N) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
             DO AA = 1, NSTREAMS
              AA1 = AA + NSTREAMS
              LSD = L_T_C_MINUS(AA,N,0,Q) *   T_DISORDS_UTDN(AA,UT) &
                    + T_C_MINUS(AA,N,0)   * L_T_DISORDS_UTDN(AA,UT,Q)
              LSU = L_T_C_PLUS(AA,N,0,Q)  *   T_DISORDS_UTUP(AA,UT) &
                    + T_C_PLUS(AA,N,0)    * L_T_DISORDS_UTUP(AA,UT,Q)
              DO S = 1, NT
               LSD = LSD + L_T_C_MINUS(AA,N,S,Q) *   XTAU_POWER(UT,S) &
                         +   T_C_MINUS(AA,N,S)   * L_XTAU_POWER(UT,S,Q)
               LSU = LSU + L_T_C_PLUS(AA,N,S,Q)  * XTAU_POWER(UT,S) &
                         +   T_C_PLUS(AA,N,S)    * L_XTAU_POWER(UT,S,Q)
              END DO
              L_UT_T_PARTIC(AA,UT,Q)  = LSD
              L_UT_T_PARTIC(AA1,UT,Q) = LSU
             END DO
            ENDDO
           ENDIF
          ENDIF

!  END OFFGRID OUTPUT

         ENDDO
        ENDIF
       END IF

!  RETURN THERMAL TRANSMITTANCE-ONLY

       RETURN

!  END THERMAL TRANSMITTANCE-ONLY CLAUSE

      ENDIF

! ---------------------------------------
! CLASSICAL SOLUTIONS FOR ALL LAYERS
! ---------------------------------------

      DO N = 1, NLAYERS

!  VARIATION CONDITIONS

        DO_VARY = DO_ATMOS_LINEARIZATION.AND.LAYER_VARY_FLAG(N)
        NVARY   = LAYER_VARY_NUMBER(N)

!  SOURCE CONSTANTS

        TERM1 = TWO * TCOM1(N,1)
        TERM2 = TWO * TCOM1(N,2)
        IF ( DO_VARY ) THEN
          DO Q = 1, NVARY
            L_TERM1(Q) = TWO * L_TCOM1(N,1,Q)
            L_TERM2(Q) = TWO * L_TCOM1(N,2,Q)
          ENDDO
        ENDIF

!  H SOLUTIONS
!  -----------

!  SOLUTION MATRIX FOR THE REDUCED PROBLEM
!  ( MATRIX SHOULD BE SAVED IN THE LU DECOMPOSITION FORM)

        DO I = 1, NSTREAMS
          DO J = 1, NSTREAMS
            TMAT(I,J) = SAB(I,J,O1,O1,N) * QUAD_STREAMS(I)
          ENDDO
        ENDDO

!  L-U DECOMPOSITION OF THE SOLUTION MATRIX

        CALL DGETRF(NSTREAMS,NSTREAMS,TMAT,MAXSTREAMS,TPIVOT,INFO)

        IF ( INFO .GT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(C3, '(I3)' ) N
          MESSAGE = 'ARGUMENT I ILLEGAL VALUE, FOR I = '//CI
          TRACE   = &
            'DGETRF CALL FOR H, THERMAL_CLSOLUTION_PLUS, LAYER '//C3
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  H VECTOR_1 AND SOLUTION BY BACK-SUBSTITUTION

        DO I = 1, NSTREAMS
          HVEC1(I) = TERM1
        ENDDO
        CALL DGETRS ('N',NSTREAMS,1,TMAT,MAXSTREAMS,TPIVOT, &
           HVEC1,MAXSTREAMS,INFO)

        IF ( INFO .NE. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(C3, '(I3)' ) N
          MESSAGE = 'ARGUMENT I ILLEGAL VALUE, FOR I = '//CI
          TRACE   = &
            'DGETRS CALL FOR H_1, THERMAL_CLSOLUTION_PLUS, LAYER '//C3
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  LINEARIZATION H VECTOR_1

        IF ( DO_VARY ) THEN
          DO Q = 1, NVARY
            DO I = 1, NSTREAMS
              SUM = ZERO
              DO J = 1, NSTREAMS
                SUM = SUM - &
                     L_SAB(I,J,O1,O1,N,Q) * QUAD_STREAMS(I) * HVEC1(J)
              ENDDO
              LHVEC1(I,Q) = L_TERM1(Q) + SUM
            ENDDO
          ENDDO
          CALL DGETRS ('N',NSTREAMS,NVARY,TMAT,MAXSTREAMS,TPIVOT, &
                        LHVEC1,MAXSTREAMS,INFO)
          IF ( INFO .NE. 0 ) THEN
            WRITE(CI, '(I3)' ) INFO
            WRITE(C3, '(I3)' ) N
            MESSAGE = 'ARGUMENT I ILLEGAL VALUE, FOR I = '//CI
            TRACE   = &
              'DGETRS CALL LH_1, THERMAL_CLSOLUTION_PLUS, LAYER '//C3
            STATUS  = VLIDORT_SERIOUS
            RETURN
          ENDIF

        ENDIF

!  H VECTOR_2 AND SOLUTION BY BACK-SUBSTITUTION

        DO I = 1, NSTREAMS
          HVEC2(I) = TERM2
        ENDDO

        CALL DGETRS &
          ('N',NSTREAMS,1,TMAT,MAXSTREAMS,TPIVOT, &
           HVEC2,MAXSTREAMS,INFO)

        IF ( INFO .NE. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(C3, '(I3)' ) N
          MESSAGE = 'ARGUMENT I ILLEGAL VALUE, FOR I = '//CI
          TRACE   = &
            'DGETRS CALL FOR H_2, THERMAL_CLSOLUTION_PLUS, LAYER '//C3
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  LINEARIZATION H VECTOR_2

        IF ( DO_VARY ) THEN
          DO Q = 1, NVARY
            DO I = 1, NSTREAMS
              SUM = ZERO
              DO J = 1, NSTREAMS
                SUM = SUM - &
                     L_SAB(I,J,O1,O1,N,Q) * QUAD_STREAMS(I) * HVEC2(J)
              ENDDO
              LHVEC2(I,Q) = L_TERM2(Q) + SUM
            ENDDO
          ENDDO
          CALL DGETRS &
          ('N',NSTREAMS,NVARY,TMAT,MAXSTREAMS,TPIVOT, &
           LHVEC2,MAXSTREAMS,INFO)

          IF ( INFO .NE. 0 ) THEN
            WRITE(CI, '(I3)' ) INFO
            WRITE(C3, '(I3)' ) N
            MESSAGE = 'ARGUMENT I ILLEGAL VALUE, FOR I = '//CI
            TRACE   = &
              'DGETRS CALL LIN H_2, THERMAL_CLSOLUTION_PLUS, LAYER '//C3
            STATUS  = VLIDORT_SERIOUS
            RETURN
          ENDIF

        ENDIF

!  J SOLUTION
!  ----------

!  SOLUTION MATRIX AND VECTOR FOR THE REDUCED PROBLEM
!  ( MATRIX SHOULD BE SAVED IN THE LU DECOMPOSITION FORM)

        DO I = 1, NSTREAMS
          DO J = 1, NSTREAMS
            TMAT(I,J) = - DAB(I,J,O1,O1,N)
          ENDDO
          JVEC1(I) = HVEC2(I)
        ENDDO

!  L-U DECOMPOSITION OF THE SOLUTION MATRIX

        CALL DGETRF(NSTREAMS,NSTREAMS,TMAT,MAXSTREAMS,TPIVOT,INFO)

        IF ( INFO .GT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(C3, '(I3)' ) N
          MESSAGE = 'ARGUMENT I ILLEGAL VALUE, FOR I = '//CI
          TRACE   = &
             'DGETRF CALL FOR J, THERMAL_CLSOLUTION_PLUS, LAYER '//C3
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF


!  J VECTOR_1 SOLUTION BY BACK-SUBSTITUTION

        CALL DGETRS &
          ('N',NSTREAMS,1,TMAT,MAXSTREAMS,TPIVOT, &
           JVEC1,MAXSTREAMS,INFO)

        IF ( INFO .NE. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(C3, '(I3)' ) N
          MESSAGE = 'ARGUMENT I ILLEGAL VALUE, FOR I = '//CI
          TRACE   = &
             'DGETRS CALL FOR J, THERMAL_CLSOLUTION_PLUS, LAYER '//C3
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  LINEARIZATION J VECTOR_1

        IF ( DO_VARY ) THEN
          DO Q = 1, NVARY
            DO I = 1, NSTREAMS
              SUM = ZERO
              DO J = 1, NSTREAMS
                SUM = SUM + L_DAB(I,J,O1,O1,N,Q) * JVEC1(J)
              ENDDO
              LJVEC1(I,Q) = LHVEC2(I,Q) + SUM
            ENDDO
          ENDDO

          CALL DGETRS &
          ('N',NSTREAMS,NVARY,TMAT,MAXSTREAMS,TPIVOT, &
           LJVEC1,MAXSTREAMS,INFO)

          IF ( INFO .NE. 0 ) THEN
            WRITE(CI, '(I3)' ) INFO
            WRITE(C3, '(I3)' ) N
            MESSAGE = 'ARGUMENT I ILLEGAL VALUE, FOR I = '//CI
            TRACE   = &
          'DGETRS CALL FOR LIN J, THERMAL_CLSOLUTION_PLUS, LAYER '//C3
            STATUS  = VLIDORT_SERIOUS
            RETURN
          ENDIF

        ENDIF

!  SET SOLUTION
!  ------------

!  EXPANSION

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          TVEC1(I,N)  = HALF * (HVEC1(I) + JVEC1(I))
          TVEC1(I1,N) = HALF * (HVEC1(I) - JVEC1(I))
          TVEC2(I,N)  = HALF * HVEC2(I)
          TVEC2(I1,N) = TVEC2(I,N)
        ENDDO

!  LINEARIZED EXPANSION

        IF ( DO_VARY ) THEN
          DO Q = 1, NVARY
            DO I = 1, NSTREAMS
              I1 = I + NSTREAMS
              L_TVEC1(I,N,Q)  = HALF * (LHVEC1(I,Q) + LJVEC1(I,Q))
              L_TVEC1(I1,N,Q) = HALF * (LHVEC1(I,Q) - LJVEC1(I,Q))
              L_TVEC2(I,N,Q)  = HALF * LHVEC2(I,Q)
              L_TVEC2(I1,N,Q) = L_TVEC2(I,N,Q)
            ENDDO
          ENDDO
        ENDIF

!  VALUES AT THE LAYER BOUNDARIES

        DO I = 1, NSTREAMS_2
          T_WUPPER(I,N) = TVEC1(I,N)
          T_WLOWER(I,N) = TVEC1(I,N) + TVEC2(I,N) * DELTAU_POWER(N,2)
        ENDDO

!  LINEARIZED VALUES AT THE LAYER BOUNDARIES

        IF ( DO_VARY ) THEN
          DO Q = 1, NVARY
            DO I = 1, NSTREAMS_2
              L_T_WUPPER(I,N,Q) = L_TVEC1(I,N,Q)
              L_T_WLOWER(I,N,Q) = L_TVEC1(I,N,Q) &
                   + L_TVEC2(I,N,Q) *   DELTAU_POWER(N,2) &
                   +   TVEC2(I,N)   * L_DELTAU_POWER(N,2,Q)
            ENDDO
          ENDDO
        ENDIF

!  USER SOLUTIONS - HELP VECTORS

        DO L = M, NMOMENTS

          SUM1 = ZERO
          SUM2 = ZERO
          DO  J = 1, NSTREAMS
            A5 = QUAD_HALFWTS(J)
            J1 = J + NSTREAMS
            POS1 = TVEC1(J1,N) * A5 * PI_XQP    (L,J,O1,O1)
            NEG1 = TVEC1(J,N)  * A5 * PI_XQM_PRE(L,J,O1,O1)
            SUM1 = SUM1 + POS1 + NEG1
            POS2 = TVEC2(J,N)  * A5 * PI_XQP    (L,J,O1,O1)
            NEG2 = TVEC2(J1,N) * A5 * PI_XQM_PRE(L,J,O1,O1)
            SUM2 = SUM2 + POS2 + NEG2
          ENDDO
          T_HELP1(L) = SUM1 * OMEGA_GREEK(L,N,O1,O1)
          T_HELP2(L) = SUM2 * OMEGA_GREEK(L,N,O1,O1)

!  LINEARIZED USER SOLUTIONS, HELP VECTORS

          IF ( DO_VARY ) THEN
            DO Q = 1, NVARY
              L_SUM1 = ZERO
              L_SUM2 = ZERO
              DO  J = 1, NSTREAMS
                J1 = J + NSTREAMS
                A5 = QUAD_HALFWTS(J)
                POS1 = L_TVEC1(J1,N,Q) * A5 * PI_XQP    (L,J,O1,O1)
                NEG1 = L_TVEC1(J,N,Q)  * A5 * PI_XQM_PRE(L,J,O1,O1)
                L_SUM1 = L_SUM1 + POS1 + NEG1
                POS2 = L_TVEC2(J,N,Q)  * A5 * PI_XQP    (L,J,O1,O1)
                NEG2 = L_TVEC2(J1,N,Q) * A5 * PI_XQM_PRE(L,J,O1,O1)
                L_SUM2 = L_SUM2 + POS2 + NEG2
              ENDDO
              L_T_HELP1(L,Q) = L_SUM1 *   OMEGA_GREEK(L,N,O1,O1) &
                               + SUM1 * L_OMEGA_GREEK(L,N,O1,O1,Q)
              L_T_HELP2(L,Q) = L_SUM2 *   OMEGA_GREEK(L,N,O1,O1) &
                               + SUM2 * L_OMEGA_GREEK(L,N,O1,O1,Q)
            ENDDO
          ENDIF
        ENDDO

!  UPWELLING: SUM OVER ALL HARMONIC CONTRIBUTIONS, EACH USER STREAM

        IF ( DO_UPWELLING.AND.STERM_LAYERMASK_UP(N) ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            POS1 = ZERO
            POS2 = ZERO
            DO L = M, NMOMENTS
              POS1 = POS1 + T_HELP1(L) * PI_XUP(L,UM,O1,O1)
              POS2 = POS2 + T_HELP2(L) * PI_XUP(L,UM,O1,O1)
            ENDDO
            U_TPOS1(UM,N) = POS1
            U_TPOS2(UM,N) = POS2
            IF ( DO_VARY ) THEN
              DO Q = 1, NVARY
                POS1 = ZERO
                POS2 = ZERO
                DO L = M, NMOMENTS
                  POS1 = POS1 + L_T_HELP1(L,Q) * PI_XUP(L,UM,O1,O1)
                  POS2 = POS2 + L_T_HELP2(L,Q) * PI_XUP(L,UM,O1,O1)
                ENDDO
                L_U_TPOS1(UM,N,Q) = POS1
                L_U_TPOS2(UM,N,Q) = POS2
              ENDDO
            ENDIF
          ENDDO
        ENDIF

!  DNWELLING: SUM OVER ALL HARMONIC CONTRIBUTIONS, EACH USER STREAM

        IF ( DO_DNWELLING.AND.STERM_LAYERMASK_DN(N) ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            NEG1 = ZERO
            NEG2 = ZERO
            DO L = M, NMOMENTS
              NEG1 = NEG1 + T_HELP1(L) * PI_XUM(L,UM,O1,O1)
              NEG2 = NEG2 + T_HELP2(L) * PI_XUM(L,UM,O1,O1)
            ENDDO
            U_TNEG1(UM,N) = NEG1
            U_TNEG2(UM,N) = NEG2
            IF ( DO_VARY ) THEN
              DO Q = 1, NVARY
                NEG1 = ZERO
                NEG2 = ZERO
                DO L = M, NMOMENTS
                  NEG1 = NEG1 + L_T_HELP1(L,Q) * PI_XUM(L,UM,O1,O1)
                  NEG2 = NEG2 + L_T_HELP2(L,Q) * PI_XUM(L,UM,O1,O1)
                ENDDO
                L_U_TNEG1(UM,N,Q) = NEG1
                L_U_TNEG2(UM,N,Q) = NEG2
              ENDDO
            ENDIF
          ENDDO
        ENDIF

!  END LAYER LOOP

      ENDDO

!  OFFGRID: ONLY FOR QUADRATURE OR MEAN-VALUE OUTPUT
!  =================================================

      IF ( DO_QUAD_OUTPUT .OR. DO_ADDITIONAL_MVOUT ) THEN
        IF ( N_PARTLAYERS .GT. 0 ) THEN
          DO UT = 1, N_PARTLAYERS
            N  = PARTLAYERS_LAYERIDX(UT)
            DO I = 1, NSTREAMS_2
              UT_T_PARTIC(I,UT) = &
                  TVEC1(I,N) + TVEC2(I,N) * XTAU_POWER(UT,2)
            ENDDO
            IF ( LAYER_VARY_FLAG(N) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                DO I = 1, NSTREAMS_2
                  L_UT_T_PARTIC(I,UT,Q) = L_TVEC1(I,N,Q) &
                    + L_TVEC2(I,N,Q) *   XTAU_POWER(UT,2) &
                    +   TVEC2(I,N)   * L_XTAU_POWER(UT,2,Q)
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDIF
      END IF



!  FINISH

      RETURN
      END SUBROUTINE THERMAL_CLSOLUTION_PLUS

!

      SUBROUTINE THERMAL_STERMS_UP_PLUS ( &
        DO_SOLAR_SOURCES, NLAYERS, &
        N_THERMAL_COEFFS, DO_THERMAL_TRANSONLY, &
        DO_ATMOS_LINEARIZATION, LAYER_VARY_FLAG, &
        LAYER_VARY_NUMBER, &
        N_USER_STREAMS, LOCAL_UM_START, &
        USER_STREAMS, DO_PARTLAYERS, &
        N_PARTLAYERS, PARTLAYERS_LAYERIDX, &
        N_ALLLAYERS_UP, STERM_LAYERMASK_UP, &
        T_DELT_USERM, T_UTUP_USERM, &
        L_T_DELT_USERM, L_T_UTUP_USERM, &
        DELTAU_POWER, XTAU_POWER, &
        U_TPOS1, U_TPOS2, &
        T_DIRECT_UP, T_UT_DIRECT_UP, &
        L_DELTAU_POWER, L_XTAU_POWER, &
        L_U_TPOS1, L_U_TPOS2, &
        L_T_DIRECT_UP, L_T_UT_DIRECT_UP, &
        LAYER_TSUP_UP, LAYER_TSUP_UTUP, &
        L_LAYER_TSUP_UP, L_LAYER_TSUP_UTUP )

!  THERMAL CONTRIBUTIONS TO LAYER SOURCE TERMS (UPWELLING)
!  LINEARIZED THERMAL CONTRIBUTIONS TO LAYER SOURCE TERMS (UPWELLING)

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_THERMAL_COEFFS
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_ATMOS_LINEARIZATION
      LOGICAL, INTENT (IN) ::          LAYER_VARY_FLAG  ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          LAYER_VARY_NUMBER ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      DOUBLE PRECISION, INTENT (IN) :: USER_STREAMS  ( MAX_USER_STREAMS )
      LOGICAL, INTENT (IN) ::          DO_PARTLAYERS
      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      INTEGER, INTENT (IN) ::          N_ALLLAYERS_UP
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_UP ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTUP_USERM &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_POWER &
          ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (IN) :: XTAU_POWER &
          ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (IN) :: U_TPOS1 ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: U_TPOS2 ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DIRECT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UT_DIRECT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_POWER &
          ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_XTAU_POWER &
          ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_U_TPOS1 &
          ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_U_TPOS2 &
          ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DIRECT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UT_DIRECT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS)

      DOUBLE PRECISION, INTENT (OUT) :: LAYER_TSUP_UP &
          ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: LAYER_TSUP_UTUP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: L_LAYER_TSUP_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_LAYER_TSUP_UTUP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  LOCAL VARIABLES

      INTEGER ::          UM, N, UT, S, NT, Q
      DOUBLE PRECISION :: SPAR, COSMUM, FAC, SUM

!  MULTIPLIERS (WHOLE LAYER)

      DOUBLE PRECISION :: T_MULT_UP(MAXLAYERS,0:MAX_THERMAL_COEFFS)
      DOUBLE PRECISION :: L_T_MULT_UP &
          (MAXLAYERS,0:MAX_THERMAL_COEFFS,MAX_ATMOSWFS)

!  PARTICULAR SOLUTION LAYER SOURCE TERMS ( GREEN'S FUNCTION SOLUTION )
!  --------------------------------------------------------------------

!  INITIAL MODULUS = 4.PI IF SOLAR SOURCES ARE INCLUDED

      FAC   = ONE
      IF ( DO_SOLAR_SOURCES ) FAC = PI4

!  SHORTHAND

      NT = N_THERMAL_COEFFS

!  START USER ANGLE LOOP

      DO UM = LOCAL_UM_START, N_USER_STREAMS

!  LOCAL COSINE

       COSMUM = USER_STREAMS(UM)

!  DIRECT TERMS TO START

       DO N = N_ALLLAYERS_UP, NLAYERS
         LAYER_TSUP_UP(UM,N) = FAC * T_DIRECT_UP(UM,N)
       ENDDO
       IF ( DO_PARTLAYERS ) THEN
        DO UT = 1, N_PARTLAYERS
         LAYER_TSUP_UTUP(UM,UT) = FAC*T_UT_DIRECT_UP(UM,UT)
        ENDDO
       ENDIF

!  LINEARIZED DIRECT TERMS TO START

       IF ( DO_ATMOS_LINEARIZATION ) THEN
        DO N = N_ALLLAYERS_UP, NLAYERS
         DO Q = 1, LAYER_VARY_NUMBER(N)
          L_LAYER_TSUP_UP(UM,N,Q) = FAC * L_T_DIRECT_UP(UM,N,Q)
         ENDDO
        ENDDO
        IF ( DO_PARTLAYERS ) THEN
         DO UT = 1, N_PARTLAYERS
          N = PARTLAYERS_LAYERIDX(UT)
          DO Q = 1, LAYER_VARY_NUMBER(N)
           L_LAYER_TSUP_UTUP(UM,UT,Q) = FAC*L_T_UT_DIRECT_UP(UM,UT,Q)
          ENDDO
         ENDDO
        ENDIF
       ENDIF

!  FINISH IF TRANSMITTANCE ONLY

       IF ( DO_THERMAL_TRANSONLY ) GO TO 678

!  CLASSICAL SECTION PARTICULAR INTEGRAL
!  -------------------------------------

!  LAYER LOOP

       DO N = 1, NLAYERS
         IF ( STERM_LAYERMASK_UP(N) ) THEN

!  WHOLE LAYER SOURCE TERMS
!   NOTE: T_DELT_USERM(N,UM) WAS INDEXED OPPOSITELY

          T_MULT_UP(N,2) = U_TPOS2(UM,N)
          T_MULT_UP(N,1) = U_TPOS1(UM,N) + COSMUM * T_MULT_UP(N,2)
          SUM = T_MULT_UP(N,1)
          SUM = SUM + T_MULT_UP(N,2) * DELTAU_POWER(N,2)
          T_MULT_UP(N,0)   = - SUM
          SPAR = T_MULT_UP(N,0) * T_DELT_USERM(N,UM) &
               + T_MULT_UP(N,1)
          LAYER_TSUP_UP(UM,N) = LAYER_TSUP_UP(UM,N) + SPAR * FAC

!  LINEARIZED WHOLE LAYER SOURCE TERMS

          IF ( LAYER_VARY_FLAG(N) ) THEN
           DO Q = 1, LAYER_VARY_NUMBER(N)
            L_T_MULT_UP(N,2,Q) = L_U_TPOS2(UM,N,Q)
            L_T_MULT_UP(N,1,Q) = L_U_TPOS1(UM,N,Q) &
                   + COSMUM * L_T_MULT_UP(N,2,Q)
            SUM = L_T_MULT_UP(N,1,Q)
            SUM = SUM + L_T_MULT_UP(N,2,Q) *   DELTAU_POWER(N,2) &
                      +   T_MULT_UP(N,2)   * L_DELTAU_POWER(N,2,Q)
            L_T_MULT_UP(N,0,Q)   = - SUM
            SPAR = T_MULT_UP(N,0)   * L_T_DELT_USERM(N,UM,Q) &
               + L_T_MULT_UP(N,0,Q) *   T_DELT_USERM(N,UM) &
               + L_T_MULT_UP(N,1,Q)
            L_LAYER_TSUP_UP(UM,N,Q) = &
                L_LAYER_TSUP_UP(UM,N,Q) + SPAR * FAC
           ENDDO
          ENDIF
         ENDIF
       ENDDO

!  PARTIAL LAYER SOURCES - ADD THERMAL DIFFUSE TERM

       IF ( DO_PARTLAYERS ) THEN
         DO UT = 1, N_PARTLAYERS
           N = PARTLAYERS_LAYERIDX(UT)
           SUM = T_MULT_UP(N,0) * T_UTUP_USERM(UT,UM)
           DO S = 1, NT
             SUM = SUM + T_MULT_UP(N,S) * XTAU_POWER(UT,S)
           END DO
           LAYER_TSUP_UTUP(UM,UT) = &
                 LAYER_TSUP_UTUP(UM,UT) + SUM*FAC
           IF ( LAYER_VARY_FLAG(N) ) THEN
             DO Q = 1, LAYER_VARY_NUMBER(N)
              SUM = L_T_MULT_UP(N,0,Q) *   T_UTUP_USERM(UT,UM) &
                    + T_MULT_UP(N,0)   * L_T_UTUP_USERM(UT,UM,Q)
              DO S = 1, NT
                SUM = SUM + L_T_MULT_UP(N,S,Q) *   XTAU_POWER(UT,S) &
                          +   T_MULT_UP(N,S)   * L_XTAU_POWER(UT,S,Q)
              END DO
              L_LAYER_TSUP_UTUP(UM,UT,Q) = &
                 L_LAYER_TSUP_UTUP(UM,UT,Q) + SUM*FAC
             ENDDO
           ENDIF
         ENDDO
       ENDIF

!  CONTINUATION POINT FOR AVOIDING SCATTERING CALCULATIONS

 678   CONTINUE

!  END USER-STREAM LOOP

      END DO

!  FINISH

      RETURN
      END SUBROUTINE THERMAL_STERMS_UP_PLUS

!

      SUBROUTINE THERMAL_STERMS_DN_PLUS ( &
        DO_SOLAR_SOURCES, NLAYERS, &
        N_THERMAL_COEFFS, DO_THERMAL_TRANSONLY, &
        DO_ATMOS_LINEARIZATION, LAYER_VARY_FLAG, &
        LAYER_VARY_NUMBER, &
        N_USER_STREAMS, LOCAL_UM_START, &
        USER_STREAMS, DO_PARTLAYERS, &
        N_PARTLAYERS, PARTLAYERS_LAYERIDX, &
        N_ALLLAYERS_DN, STERM_LAYERMASK_DN, &
        T_DELT_USERM, T_UTDN_USERM, &
        L_T_DELT_USERM, L_T_UTDN_USERM, &
        DELTAU_POWER, XTAU_POWER, &
        U_TNEG1, U_TNEG2, &
        T_DIRECT_DN, T_UT_DIRECT_DN, &
        L_DELTAU_POWER, L_XTAU_POWER, &
        L_U_TNEG1, L_U_TNEG2, &
        L_T_DIRECT_DN, L_T_UT_DIRECT_DN, &
        LAYER_TSUP_DN, LAYER_TSUP_UTDN, &
        L_LAYER_TSUP_DN, L_LAYER_TSUP_UTDN )

!  THERMAL CONTRIBUTIONS TO LAYER SOURCE TERMS (DOWNWELLING)
!  LINEARIZED THERMAL CONTRIBUTIONS TO LAYER SOURCE TERMS (DOWNWELLING)

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_THERMAL_COEFFS
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_ATMOS_LINEARIZATION
      LOGICAL, INTENT (IN) ::          LAYER_VARY_FLAG  ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          LAYER_VARY_NUMBER ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      DOUBLE PRECISION, INTENT (IN) :: USER_STREAMS  ( MAX_USER_STREAMS )
      LOGICAL, INTENT (IN) ::          DO_PARTLAYERS
      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      INTEGER, INTENT (IN) ::          N_ALLLAYERS_DN
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_DN ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_USERM &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_POWER &
          ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (IN) :: XTAU_POWER &
          ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (IN) :: U_TNEG1 ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: U_TNEG2 ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DIRECT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UT_DIRECT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_POWER &
          ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_XTAU_POWER &
          ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_U_TNEG1 &
          ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_U_TNEG2 &
          ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DIRECT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UT_DIRECT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (OUT) :: LAYER_TSUP_DN &
          ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: LAYER_TSUP_UTDN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: L_LAYER_TSUP_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_LAYER_TSUP_UTDN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  LOCAL VARIABLES

      INTEGER ::          UM, N, UT, S, NT, Q
      DOUBLE PRECISION :: SPAR, COSMUM, FAC, SUM

!  MULTIPLIERS (WHOLE LAYER)

      DOUBLE PRECISION :: T_MULT_DN(MAXLAYERS,0:MAX_THERMAL_COEFFS)
      DOUBLE PRECISION :: L_T_MULT_DN &
          (MAXLAYERS,0:MAX_THERMAL_COEFFS,MAX_ATMOSWFS)

!  PARTICULAR SOLUTION LAYER SOURCE TERMS ( GREEN'S FUNCTION SOLUTION )
!  --------------------------------------------------------------------

!  INITIAL MODULUS = 4.PI IF SOLAR SOURCES ARE INCLUDED

      FAC = ONE
      IF ( DO_SOLAR_SOURCES ) FAC = PI4

!  SHORTHAND

      NT = N_THERMAL_COEFFS

!  START USER ANGLE LOOP

      DO UM = LOCAL_UM_START, N_USER_STREAMS

!  LOCAL COSINE

       COSMUM = USER_STREAMS(UM)

!  DIRECT TERMS TO START

       DO N = 1, N_ALLLAYERS_DN
         LAYER_TSUP_DN(UM,N) = FAC * T_DIRECT_DN(UM,N)
       ENDDO
       IF ( DO_PARTLAYERS ) THEN
        DO UT = 1, N_PARTLAYERS
         LAYER_TSUP_UTDN(UM,UT) = FAC*T_UT_DIRECT_DN(UM,UT)
        ENDDO
       ENDIF

!  LINEARIZED DIRECT TERMS TO START

       IF ( DO_ATMOS_LINEARIZATION ) THEN
        DO N = 1, N_ALLLAYERS_DN
         DO Q = 1, LAYER_VARY_NUMBER(N)
          L_LAYER_TSUP_DN(UM,N,Q) = FAC * L_T_DIRECT_DN(UM,N,Q)
         ENDDO
        ENDDO
        IF ( DO_PARTLAYERS ) THEN
         DO UT = 1, N_PARTLAYERS
          N = PARTLAYERS_LAYERIDX(UT)
          DO Q = 1, LAYER_VARY_NUMBER(N)
           L_LAYER_TSUP_UTDN(UM,UT,Q) = FAC*L_T_UT_DIRECT_DN(UM,UT,Q)
          ENDDO
         ENDDO
        ENDIF
       ENDIF

!  FINISH IF TRANSMITTANCE ONLY

       IF ( DO_THERMAL_TRANSONLY ) GO TO 678

!  CLASSICAL SECTION PARTICULAR INTEGRAL
!  =====================================

       DO N = 1, NLAYERS
         IF ( STERM_LAYERMASK_DN(N) ) THEN

!  WHOLE LAYER SOURCE TERMS
!   NOTE: T_DELT_USERM(N,UM) WAS INDEXED OPPOSITELY

          T_MULT_DN(N,2) = U_TNEG2(UM,N)
          T_MULT_DN(N,1) = U_TNEG1(UM,N) - COSMUM * T_MULT_DN(N,2)
          T_MULT_DN(N,0) = - T_MULT_DN(N,1)
          SPAR = T_MULT_DN(N,0) * T_DELT_USERM(N,UM)
          SPAR = SPAR + T_MULT_DN(N,1)
          SPAR = SPAR + T_MULT_DN(N,2) * DELTAU_POWER(N,2)
          LAYER_TSUP_DN(UM,N) = LAYER_TSUP_DN(UM,N) + SPAR * FAC

          IF ( LAYER_VARY_FLAG(N) ) THEN
           DO Q = 1, LAYER_VARY_NUMBER(N)
            L_T_MULT_DN(N,2,Q) = L_U_TNEG2(UM,N,Q)
            L_T_MULT_DN(N,1,Q) = L_U_TNEG1(UM,N,Q) &
                   - COSMUM * L_T_MULT_DN(N,2,Q)
            L_T_MULT_DN(N,0,Q) = - L_T_MULT_DN(N,1,Q)
            SPAR =     T_MULT_DN(N,0)   * L_T_DELT_USERM(N,UM,Q) &
                   + L_T_MULT_DN(N,0,Q) *   T_DELT_USERM(N,UM)
            SPAR = SPAR + L_T_MULT_DN(N,1,Q)
            SPAR = SPAR + L_T_MULT_DN(N,2,Q) *   DELTAU_POWER(N,2) &
                        +   T_MULT_DN(N,2)   * L_DELTAU_POWER(N,2,Q)
            L_LAYER_TSUP_DN(UM,N,Q) = &
                  L_LAYER_TSUP_DN(UM,N,Q) + SPAR * FAC
           ENDDO
          ENDIF
         ENDIF
       ENDDO

!  PARTIAL LAYER SOURCE TERMS, ADD THERMAL DIFFUSE PART

       IF ( DO_PARTLAYERS ) THEN
         DO UT = 1, N_PARTLAYERS
           N  = PARTLAYERS_LAYERIDX(UT)
           SUM = T_MULT_DN(N,0) * T_UTDN_USERM(UT,UM)
           DO S = 1, NT
             SUM = SUM + T_MULT_DN(N,S) * XTAU_POWER(UT,S)
           END DO
           LAYER_TSUP_UTDN(UM,UT) = &
                 LAYER_TSUP_UTDN(UM,UT) + SUM*FAC
           IF ( LAYER_VARY_FLAG(N) ) THEN
             DO Q = 1, LAYER_VARY_NUMBER(N)
              SUM = L_T_MULT_DN(N,0,Q) *   T_UTDN_USERM(UT,UM) &
                    + T_MULT_DN(N,0)   * L_T_UTDN_USERM(UT,UM,Q)
              DO S = 1, NT
                SUM = SUM + L_T_MULT_DN(N,S,Q) *   XTAU_POWER(UT,S) &
                          +   T_MULT_DN(N,S)   * L_XTAU_POWER(UT,S,Q)
              END DO
              L_LAYER_TSUP_UTDN(UM,UT,Q) = &
                   L_LAYER_TSUP_UTDN(UM,UT,Q) + SUM*FAC
             ENDDO
           ENDIF
         END DO
       END IF

!  CONTINUATION POINT FOR AVOIDING SCATTERING CALCULATIONS

 678   CONTINUE

!  END USER-STREAM LOOP

      ENDDO

!  FINISH

      RETURN
      END SUBROUTINE THERMAL_STERMS_DN_PLUS

!

      SUBROUTINE THERMAL_LTE_LINEARIZATION ( &
        DO_INCLUDE_SURFACE, SURFACE_FACTOR, &
        FLUX_MULTIPLIER, &
        DO_UPWELLING, DO_DNWELLING, &
        NSTREAMS, NLAYERS, N_USER_LEVELS, &
        DELTAU_VERT_INPUT, LAMBERTIAN_ALBEDO, &
        THERMAL_BB_INPUT, &
        QUAD_STREAMS, QUAD_STRMWTS, &
        N_USER_STREAMS, LOCAL_UM_START, USER_STREAMS, &
        UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, &
        STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, &
        T_DELT_DISORDS, T_DELT_USERM, &
        CUMSOURCE_UP, CUMSOURCE_DN, &
        THERMCOEFFS, &
        LTE_DELTAU_VERT_INPUT, LTE_THERMAL_BB_INPUT, &
        LTE_ATMOSWF )

!  SET-UP OF LTE LINEARIZATION OF THERMAL EXPANSION COEFFICIENTS.

!  LINEARIZATION W.R.T PROFILE OPTICAL DEPTH AND THE BB INPUT.

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_SURFACE
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER
      DOUBLE PRECISION, INTENT (IN) :: SURFACE_FACTOR
      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_DNWELLING
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: THERMAL_BB_INPUT ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      DOUBLE PRECISION, INTENT (IN) :: USER_STREAMS  ( MAX_USER_STREAMS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_DN ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: CUMSOURCE_UP &
          ( MAX_USER_STREAMS, MAXSTOKES, 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: CUMSOURCE_DN &
          ( MAX_USER_STREAMS, MAXSTOKES, 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: THERMCOEFFS &
          ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (IN) :: LTE_DELTAU_VERT_INPUT ( 2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LTE_THERMAL_BB_INPUT ( 0:MAXLAYERS )

      DOUBLE PRECISION, INTENT (OUT) :: LTE_ATMOSWF &
          ( 0:MAXLAYERS, MAX_USER_LEVELS, &
            MAX_USER_VZANGLES, MAX_DIRECTIONS )

!  LOCAL VARIABLES
!  ---------------

      INTEGER ::          K, I, UI, Q, UTA
      INTEGER ::          N, NC, NUT, NSTART, NUT_PREV, NLEVEL
      DOUBLE PRECISION :: MU, XI, TERM1, TERM2, KT, KT1, AMB, BD, KD
      DOUBLE PRECISION :: LTE_KT, LTE_AMB, LTE_BD, FINAL_SOURCE
      DOUBLE PRECISION :: KMULT, REFLEC

!  LOCAL CUMULATIVE SOURCE

      DOUBLE PRECISION :: LTE_CUMSOURCE ( MAX_USER_STREAMS )

!  LINEARIZED COEFFICIENTS

      DOUBLE PRECISION :: LTE_TCOEFFS(2,MAXLAYERS,MAX_THERMAL_COEFFS)

!  LINEARIZED TRANSMITTANCES

      DOUBLE PRECISION :: LTE_TRANS_DISORDS(2,MAXLAYERS,MAXSTREAMS)
      DOUBLE PRECISION :: LTE_TRANS_USERM  (2,MAXLAYERS,MAX_USER_STREAMS)

!  LINEARIZED SOURCE TERMS (USER STREAMS)

      DOUBLE PRECISION :: LTE_T_DIRECT_UP(2,MAXLAYERS,MAX_USER_STREAMS)
      DOUBLE PRECISION :: LTE_T_DIRECT_DN(2,MAXLAYERS,MAX_USER_STREAMS)

!  SOURCE TERMS (AND LINEARIZED) FOR DISCRETE ORDINATES

      DOUBLE PRECISION :: T_DISORDS_DN (MAXLAYERS,MAXSTREAMS)
      DOUBLE PRECISION :: LTE_T_DISORDS_DN(2,MAXLAYERS,MAXSTREAMS)

!  SURFACE REFLECTION STUFF

      DOUBLE PRECISION :: DOWNCUMS (0:MAXLAYERS,MAXSTREAMS)
      DOUBLE PRECISION :: LTE_DOWNSURF (MAXSTREAMS)
      DOUBLE PRECISION :: LTE_BOA_SOURCES (0:MAXLAYERS,MAXSTREAMS)

!  LTE-LINEARIZED COEFFICIENTS
!  ---------------------------

      DO N = 1, NLAYERS
         LTE_TCOEFFS(1,N,1) = LTE_THERMAL_BB_INPUT(N-1)
         LTE_TCOEFFS(2,N,1) = ZERO
         LTE_TCOEFFS(1,N,2) = - LTE_THERMAL_BB_INPUT(N-1) &
                    - THERMCOEFFS(N,2) * LTE_DELTAU_VERT_INPUT(1,N)
         LTE_TCOEFFS(2,N,2) = + LTE_THERMAL_BB_INPUT(N) &
                    - THERMCOEFFS(N,2) * LTE_DELTAU_VERT_INPUT(2,N)
         LTE_TCOEFFS(1,N,2) = LTE_TCOEFFS(1,N,2)/DELTAU_VERT_INPUT(N)
         LTE_TCOEFFS(2,N,2) = LTE_TCOEFFS(2,N,2)/DELTAU_VERT_INPUT(N)
      ENDDO

!  LTE-LINEARIZATION OF TRANSMITTANCES
!  -----------------------------------

!  USER STREAMS, ALWAYS NEED THIS

      DO UI = 1, N_USER_STREAMS
        MU = USER_STREAMS(UI)
        DO N = 1, NLAYERS
          KT  = T_DELT_USERM(N,UI)
          KD = - KT / MU
          DO Q = 1, 2
            LTE_TRANS_USERM(Q,N,UI) = KD * LTE_DELTAU_VERT_INPUT(Q,N)
          ENDDO
        ENDDO
      ENDDO

!  DISCRETE ORDINATE STSREAMS, ONLY FOR SURFACE CONTRIBUTION FOR UPWELLING

      IF ( DO_UPWELLING .AND. DO_INCLUDE_SURFACE ) THEN
        DO I = 1, NSTREAMS
          XI = QUAD_STREAMS(I)
          DO N = 1, NLAYERS
            KT  = T_DELT_DISORDS(I,N)
            KD = - KT / XI
            DO Q = 1, 2
              LTE_TRANS_DISORDS(Q,N,I) = KD * LTE_DELTAU_VERT_INPUT(Q,N)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  UPWELLING DIRECT SOLUTION SOURCE TERMS
!  --------------------------------------

      IF ( DO_UPWELLING ) THEN
        DO UI = 1, N_USER_STREAMS
          MU = USER_STREAMS(UI)
          DO N = 1, NLAYERS
            IF ( STERM_LAYERMASK_UP(N) ) THEN
              AMB = THERMCOEFFS(N,1) + MU * THERMCOEFFS(N,2)
              BD  = THERMCOEFFS(N,2) * DELTAU_VERT_INPUT(N)
              KT  = T_DELT_USERM(N,UI)
              KT1 = ONE - KT
              DO Q = 1, 2
                LTE_KT  = LTE_TRANS_USERM(Q,N,UI)
                LTE_AMB = LTE_TCOEFFS(Q,N,1) + MU * LTE_TCOEFFS(Q,N,2)
                LTE_BD  = THERMCOEFFS(N,2) * LTE_DELTAU_VERT_INPUT(Q,N) &
                      + LTE_TCOEFFS(Q,N,2) * DELTAU_VERT_INPUT(N)
                TERM1 = - KT  * LTE_BD  - LTE_KT * BD
                TERM2 = + KT1 * LTE_AMB - LTE_KT * AMB
                LTE_T_DIRECT_UP(Q,N,UI) = TERM1 + TERM2
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDIF

!  DOWNWELLING DISCRETE ORDINATE SOURCE TERMS
!  ------------------------------------------

!  THESE ARE REQUIRED IF THERE IS SURFACE REFLECTANCE

      IF ( DO_UPWELLING .AND. DO_INCLUDE_SURFACE ) THEN
        DO I = 1, NSTREAMS
          XI = QUAD_STREAMS(I)
          DO N = 1, NLAYERS
            AMB = THERMCOEFFS(N,1) - XI * THERMCOEFFS(N,2)
            BD  = THERMCOEFFS(N,2) * DELTAU_VERT_INPUT(N)
            KT  = T_DELT_DISORDS(I,N)
            KT1 = ONE - KT
            DO Q = 1, 2
              LTE_KT  = LTE_TRANS_DISORDS(Q,N,I)
              LTE_AMB = LTE_TCOEFFS(Q,N,1) - XI * LTE_TCOEFFS(Q,N,2)
              LTE_BD  = THERMCOEFFS(N,2)   * LTE_DELTAU_VERT_INPUT(Q,N) &
                      + LTE_TCOEFFS(Q,N,2) * DELTAU_VERT_INPUT(N)
              TERM1 = + LTE_BD
              TERM2 = + KT1 * LTE_AMB - LTE_KT * AMB
              LTE_T_DISORDS_DN(Q,N,I) = TERM1 + TERM2
            ENDDO
            T_DISORDS_DN(N,I) = AMB * KT1 + BD
          ENDDO
        ENDDO
      ENDIF

!  DOWNWELLING DIRECT SOLUTION SOURCE TERMS
!  ----------------------------------------

      IF ( DO_DNWELLING ) THEN
        DO UI = 1, N_USER_STREAMS
          MU = USER_STREAMS(UI)
          DO N = 1, NLAYERS
            IF ( STERM_LAYERMASK_DN(N) ) THEN
              AMB = THERMCOEFFS(N,1) - MU * THERMCOEFFS(N,2)
              BD  = THERMCOEFFS(N,2) * DELTAU_VERT_INPUT(N)
              KT  = T_DELT_USERM(N,UI)
              KT1 = ONE - KT
              DO Q = 1, 2
                LTE_KT  = LTE_TRANS_USERM(Q,N,UI)
                LTE_AMB = LTE_TCOEFFS(Q,N,1) - MU * LTE_TCOEFFS(Q,N,2)
                LTE_BD  = THERMCOEFFS(N,2)  * LTE_DELTAU_VERT_INPUT(Q,N) &
                       + LTE_TCOEFFS(Q,N,2) * DELTAU_VERT_INPUT(N)
                TERM1 = + LTE_BD
                TERM2 = + KT1 * LTE_AMB - LTE_KT * AMB
                LTE_T_DIRECT_DN(Q,N,UI) = TERM1 + TERM2
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDIF

!  ==============================
!  NOW DO THE WEIGHTING FUNCTIONS
!  ==============================

!  FOR UPWELLING, HAVE TO DO LINEARIZATION OF STOKES_DOWNSURF FIRST

      IF ( DO_UPWELLING ) THEN
        IF ( DO_INCLUDE_SURFACE ) THEN

! DEVELOP THE CUMULATIVE SOURCES

          DO I = 1, NSTREAMS
            DOWNCUMS(0,I) = ZERO
          ENDDO
          DO N = 1, NLAYERS
            DO I = 1, NSTREAMS
              DOWNCUMS(N,I) = T_DISORDS_DN(N,I) + &
                         T_DELT_DISORDS(I,N) * DOWNCUMS(N-1,I)
            ENDDO
          ENDDO

!  FOR EACH LTE WEIGHTING FUNCTION

          DO K = 0, NLAYERS

!  LINEARIZE THE DOWNWELLING SOLUTION

            DO I = 1, NSTREAMS
              LTE_DOWNSURF(I) = ZERO
            ENDDO
            DO N = 1, NLAYERS
              IF ( K .EQ. N ) THEN
                DO I = 1, NSTREAMS
                  LTE_DOWNSURF(I) = LTE_T_DISORDS_DN(2,N,I) + &
                        T_DELT_DISORDS(I,N)  * LTE_DOWNSURF(I) + &
                    LTE_TRANS_DISORDS(2,N,I) * DOWNCUMS(N-1,I)
                ENDDO
              ELSE IF ( K .EQ. N-1 ) THEN
                DO I = 1, NSTREAMS
                  LTE_DOWNSURF(I) = LTE_T_DISORDS_DN(1,N,I) + &
                        T_DELT_DISORDS(I,N)  * LTE_DOWNSURF(I) + &
                    LTE_TRANS_DISORDS(1,N,I) * DOWNCUMS(N-1,I)
                ENDDO
              ELSE
                DO I = 1, NSTREAMS
                  LTE_DOWNSURF(I) = &
                        T_DELT_DISORDS(I,N)   * LTE_DOWNSURF(I)
                ENDDO
              ENDIF
            ENDDO

!  REFLECT OFF THE SURFACE (LAMBERTIAN ONLY)

            KMULT = SURFACE_FACTOR * LAMBERTIAN_ALBEDO
            REFLEC = 0.0D0
            DO I = 1, NSTREAMS
              REFLEC = REFLEC + QUAD_STRMWTS(I) * LTE_DOWNSURF(I)
            ENDDO
            REFLEC = REFLEC * KMULT

!  SET THE BOA SOURCE TERMS (LAMBERTIAN ONLY)

            DO UI = 1, N_USER_STREAMS
              LTE_BOA_SOURCES(K,UI) = REFLEC
            ENDDO

!  END LTE WEIGHTING FUNCTION LOOP

          ENDDO

!  END SURFACE REFLECTION CLAUSE

        ENDIF
      ENDIF

!  ZERO THIS CONTRIBUTION IF NO SURFACE

      IF ( .NOT. DO_INCLUDE_SURFACE ) THEN
        DO K = 0, NLAYERS
          DO UI = 1, N_USER_STREAMS
            LTE_BOA_SOURCES(K,UI) = 0.0D0
          ENDDO
        ENDDO
      ENDIF

!  START UPWELLING LTE WEIGHTING FUNCTIONS

      IF ( DO_UPWELLING ) THEN

        DO K = 0, NLAYERS

!  INITIALISE CUMULATIVE SOURCE TERM LOOP

          NC  = 0
          NUT = 0
          NSTART   = NLAYERS
          NUT_PREV = NSTART + 1

!  INITIALISE SOURCE

          DO UI = LOCAL_UM_START, N_USER_STREAMS
            LTE_CUMSOURCE(UI) = LTE_BOA_SOURCES(K,UI)
          ENDDO

!  LOOP OVER ALL OUTPUT OPTICAL DEPTHS
!  -----------------------------------

          DO UTA = N_USER_LEVELS, 1, -1
            NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
            NUT = NLEVEL + 1

!  CUMULATIVE LAYER PROCESSING

            DO N = NSTART, NUT, -1
              NC = NLAYERS + 1 - N
              IF ( K .EQ. N ) THEN
                DO UI = LOCAL_UM_START, N_USER_STREAMS
                  LTE_CUMSOURCE(UI) = LTE_T_DIRECT_UP(2,N,UI) + &
                      T_DELT_USERM(N,UI)  * LTE_CUMSOURCE(UI) + &
                  LTE_TRANS_USERM(2,N,UI) * CUMSOURCE_UP(UI,1,NC-1)
                ENDDO
              ELSE IF ( K .EQ. N -1 ) THEN
                DO UI = LOCAL_UM_START, N_USER_STREAMS
                  LTE_CUMSOURCE(UI) = LTE_T_DIRECT_UP(1,N,UI)  + &
                       T_DELT_USERM(N,UI)  * LTE_CUMSOURCE(UI) + &
                   LTE_TRANS_USERM(1,N,UI) * CUMSOURCE_UP(UI,1,NC-1)
                ENDDO
              ELSE
                DO UI = LOCAL_UM_START, N_USER_STREAMS
                  LTE_CUMSOURCE(UI) = &
                     T_DELT_USERM(N,UI) * LTE_CUMSOURCE(UI)
                ENDDO
              ENDIF
            ENDDO

!  USER-DEFINED STREAM OUTPUT, JUST SET TO THE CUMULATIVE SOURCE TERM

            DO UI = LOCAL_UM_START, N_USER_STREAMS
              FINAL_SOURCE = FLUX_MULTIPLIER * LTE_CUMSOURCE(UI)
              LTE_ATMOSWF(K,UTA,UI,UPIDX) = FINAL_SOURCE
            ENDDO

!  CHECK FOR UPDATING THE RECURSION

            IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
            NUT_PREV = NUT

!  END LOOP OVER OPTICAL DEPTH

          ENDDO

!  END UPWELLING LTE WEIGHTING FUNCTIONS

        ENDDO
      ENDIF

!  START DOWNWELLING LTE WEIGHTING FUNCTIONS
!  -----------------------------------------

      IF ( DO_DNWELLING ) THEN

        DO K = 0, NLAYERS

!  INITIALISE CUMULATIVE SOURCE TERM LOOP

          NC  = 0
          NUT = 0
          NSTART   = 1
          NUT_PREV = NSTART - 1

!  INITIALISE SOURCE

          DO UI = LOCAL_UM_START, N_USER_STREAMS
            LTE_CUMSOURCE(UI) = ZERO
          ENDDO

!  LOOP OVER ALL OUTPUT OPTICAL DEPTHS

          DO UTA = 1, N_USER_LEVELS
            NLEVEL = UTAU_LEVEL_MASK_DN(UTA)
            NUT = NLEVEL

!  CUMULATIVE LAYER PROCESSING

            DO N = NSTART, NUT
              NC = N
              IF ( K .EQ. N ) THEN
                DO UI = LOCAL_UM_START, N_USER_STREAMS
                  LTE_CUMSOURCE(UI) = LTE_T_DIRECT_DN(2,N,UI)  + &
                      T_DELT_USERM(N,UI)  * LTE_CUMSOURCE(UI) + &
                  LTE_TRANS_USERM(2,N,UI) * CUMSOURCE_DN(UI,1,NC-1)
                ENDDO
              ELSE IF ( K .EQ. N -1 ) THEN
                DO UI = LOCAL_UM_START, N_USER_STREAMS
                  LTE_CUMSOURCE(UI) = LTE_T_DIRECT_DN(1,N,UI) + &
                      T_DELT_USERM(N,UI)  * LTE_CUMSOURCE(UI) + &
                  LTE_TRANS_USERM(1,N,UI) * CUMSOURCE_DN(UI,1,NC-1)
                ENDDO
              ELSE
                DO UI = LOCAL_UM_START, N_USER_STREAMS
                  LTE_CUMSOURCE(UI) = &
                      T_DELT_USERM(N,UI) * LTE_CUMSOURCE(UI)
                ENDDO
              ENDIF
            ENDDO

!  USER-DEFINED STREAM OUTPUT, JUST SET TO THE CUMULATIVE SOURCE TERM

            DO UI = LOCAL_UM_START, N_USER_STREAMS
              FINAL_SOURCE = FLUX_MULTIPLIER * LTE_CUMSOURCE(UI)
              LTE_ATMOSWF(K,UTA,UI,DNIDX) = FINAL_SOURCE
            ENDDO

!  CHECK FOR UPDATING THE RECURSION

            IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
            NUT_PREV = NUT

!  END LOOP OVER OPTICAL DEPTH

          ENDDO

!  END DOWNWELLING LTE WEIGHTING FUNCTIONS

        ENDDO
      ENDIF

!  FINISH

      RETURN
      END SUBROUTINE THERMAL_LTE_LINEARIZATION


      END MODULE vlidort_l_thermalsup

