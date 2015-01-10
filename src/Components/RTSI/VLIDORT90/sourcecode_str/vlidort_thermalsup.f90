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
! #          THERMAL_SETUP                                      #
! #                                                             #
! #   discrete ordinate particular integral                     #
! #          THERMAL_CLSOLUTION                                 #
! #                                                             #
! #   postprocessing source terms                               #
! #          THERMAL_STERMS_UP                                  #
! #          THERMAL_STERMS_DN                                  #
! #                                                             #
! ###############################################################


      MODULE vlidort_thermalsup

      PRIVATE
      PUBLIC :: THERMAL_SETUP, &
                THERMAL_CLSOLUTION, &
                THERMAL_STERMS_UP, &
                THERMAL_STERMS_DN

      CONTAINS

      SUBROUTINE THERMAL_SETUP ( &
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
        THERMCOEFFS, DELTAU_POWER, &
        XTAU_POWER, &
        TCOM1, T_DIRECT_UP, &
        T_DIRECT_DN, T_UT_DIRECT_UP, &
        T_UT_DIRECT_DN )

!  SET-UP OF THERMAL EXPANSION COEFFICIENTS, ALWAYS DONE AFTER DELTA-M.

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

!  LOCAL VARIABLES
!  ---------------

      INTEGER ::          N, S, UT, UM, NT
      DOUBLE PRECISION :: HELP ( MAXLAYERS )
      DOUBLE PRECISION :: XTAU, SUM, COSMUM
      DOUBLE PRECISION :: OMEGAS1 ( MAXLAYERS )

!  MULTIPLIERS (WHOLE LAYER)

      DOUBLE PRECISION :: T_MULT_UP ( MAXLAYERS, 0:MAX_THERMAL_COEFFS )
      DOUBLE PRECISION :: T_MULT_DN ( MAXLAYERS, 0:MAX_THERMAL_COEFFS )

!  POWERS OF OPTICAL THICKNESS
!  ---------------------------

!  WHOLE LAYER

      DO N = 1, NLAYERS
       DELTAU_POWER(N,1) = ONE
       DO S = 2, N_THERMAL_COEFFS
        DELTAU_POWER(N,S) = DELTAU_VERT(N) * DELTAU_POWER(N,S-1)
       END DO
      END DO

!  PARTIAL LAYER

      DO UT = 1, N_PARTLAYERS
       XTAU = PARTAU_VERT(UT)
       XTAU_POWER(UT,1) = ONE
       DO S = 2, N_THERMAL_COEFFS
        XTAU_POWER(UT,S) = XTAU * XTAU_POWER(UT,S-1)
       END DO
      END DO

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

!  DERIVATIVE CONTINUITY FOR QUADRATIC REGIME
!    ( FIRST LAYER IS LINEAR; BETTER THAN USING A FREE DERIVATIVE AT TOA
!  IF ( N_THERMAL_COEFFS == 3 ) THEN
!     THERMCOEFFS(1,3) = ZERO
!     THERMCOEFFS(1,2) = HELP(1)
!     DO N = 1, N_COMP_LAYERS - 1
!        N1 = N + 1
!        THERMCOEFFS(N1,2) = THERMCOEFFS(N,2) + TWO * DELTAUS(N) * THERM
!        THERMCOEFFS(N1,3) = ( HELP(N1) - THERMCOEFFS(N1,2) ) / DELTAUS(
!     END DO
!  END IF

!  ALTERNATIVE SCHEME: BACKWARD LOOKING QUADRATICS
!    ( FIRST LAYER IS LINEAR; BETTER THAN USING A FREE DERIVATIVE AT TOA

!      IF ( N_THERMAL_COEFFS == 3 ) THEN
!       THERMCOEFFS(1,3) = ZERO
!       THERMCOEFFS(1,2) = HELP(1)
!       DO N = 1, NLAYERS - 1
!        N1 = N + 1
!        SUM = ( (THERMCOEFFS(N,1)-THERMCOEFFS(N1,1)) / DELTAU_VERT(N) )
!     &         + HELP(N1)
!       THERMCOEFFS(N1,3) = SUM / ( DELTAU_VERT(N1) + DELTAU_VERT(N) )
!        THERMCOEFFS(N1,2) = HELP(N1) - DELTAU_VERT(N1)*THERMCOEFFS(N1,3
!       END DO
!      END IF

!  DEBUG CHECK

!      XTAU = ZERO
!      DO N = 1, NLAYERS
!       SUM = DELTAU_VERT(N) / 10.0
!       DO S = 1, 10
!        XTAU = XTAU + SUM
!        WRITE(34,'(1P2E15.6)')
!     &      XTAU,THERMCOEFFS(N,1)+THERMCOEFFS(N,2)*S*SUM+
!     &         THERMCOEFFS(N,3)*S*SUM*S*SUM
!       END DO
!      END DO

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

! RETURN IF POST-PROCESSING NOT FLAGGED

      IF ( .NOT. DO_USER_STREAMS ) RETURN

!  ZERO DIRECT SOLUTIONS if working in MSMODE only, then return

      IF ( DO_MSMODE_THERMAL ) THEN
         T_DIRECT_UP = 0.0d0 ; T_UT_DIRECT_UP = 0.0d0
         T_DIRECT_DN = 0.0d0 ; T_UT_DIRECT_DN = 0.0d0
         RETURN
      ENDIF

!  SHORT HAND

      NT = N_THERMAL_COEFFS

!  UPWELLING DIRECT SOLUTION SOURCE TERMS
! ---------------------------------------

      IF ( DO_UPWELLING ) THEN

!  START USER STREAM LOOP

       DO UM = 1, N_USER_STREAMS
        COSMUM = USER_STREAMS(UM)

!  DIRECT SOLUTION: WHOLE LAYER SOURCE TERMS
!   NOTE: T_DELT_USERM(N,UM) WAS INDEXED OPPOSITELY

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

!  END USER STREAM LOOP, AND UPWELLING

       ENDDO
      ENDIF

!  DOWNWELLING DIRECT SOLUTION SOURCE TERMS
! -----------------------------------------

      IF ( DO_DNWELLING ) THEN

!  START USER STREAM LOOP

       DO UM = 1, N_USER_STREAMS
        COSMUM = USER_STREAMS(UM)

!  DIRECT SOLUTION: WHOLE LAYER SOURCE TERMS
!   NOTE: T_DELT_USERM(N,UM) WAS INDEXED OPPOSITELY

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
         END IF
        END DO

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
        END IF

!  END USER STREAM LOOP, AND DOWNWELLING

       ENDDO
      ENDIF

!  FINISH

      RETURN
      END SUBROUTINE THERMAL_SETUP

!

      SUBROUTINE THERMAL_CLSOLUTION ( &
        DO_UPWELLING, DO_DNWELLING, &
        DO_QUAD_OUTPUT, DO_ADDITIONAL_MVOUT, &
        NSTREAMS, NLAYERS, &
        N_THERMAL_COEFFS, DO_THERMAL_TRANSONLY, &
        QUAD_STREAMS, QUAD_HALFWTS, &
        NMOMENTS, NSTREAMS_2, &
        N_USER_STREAMS, LOCAL_UM_START, &
        N_PARTLAYERS, &
        PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, &
        STERM_LAYERMASK_DN, &
        OMEGA_GREEK, T_DELT_DISORDS, &
        T_DISORDS_UTUP, T_DISORDS_UTDN, &
        PI_XQP, PI_XUP, PI_XUM, PI_XQM_PRE, &
        SAB, DAB, &
        THERMCOEFFS, DELTAU_POWER, &
        XTAU_POWER, TCOM1, &
        T_WUPPER, T_WLOWER, &
        UT_T_PARTIC, U_TPOS1, &
        U_TNEG1, U_TPOS2, &
        U_TNEG2, &
        STATUS, MESSAGE, TRACE )

!  CLASSICAL FUNCTION THERMAL PARTICULAR INTEGRAL, ALL LAYERS.
!  USES COEFFICIENT EXPANSION OF ATTENUATION.

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
      DOUBLE PRECISION, INTENT (IN) :: SAB &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DAB &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: THERMCOEFFS &
          ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_POWER &
          ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (IN) :: XTAU_POWER &
          ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (IN) :: TCOM1 ( MAXLAYERS, MAX_THERMAL_COEFFS )

      DOUBLE PRECISION, INTENT (OUT) ::  T_WUPPER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) ::  T_WLOWER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) ::  UT_T_PARTIC &
          ( MAXSTREAMS_2, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (OUT) ::  U_TPOS1 ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) ::  U_TNEG1 ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) ::  U_TPOS2 ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) ::  U_TNEG2 ( MAX_USER_STREAMS, MAXLAYERS )

      INTEGER, INTENT (OUT) ::           STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

! LOCAL VARIABLES
! ---------------

      INTEGER ::           AA, AA1, I, I1, J, J1, L, M, O1
      INTEGER ::           S, N, NT, INFO, UM, UT
      DOUBLE PRECISION  :: SUM_M, SUM_P, TK, K1, SD, SU, A5
      DOUBLE PRECISION  :: TERM1, TERM2, SUM1, SUM2
      DOUBLE PRECISION  :: POS1, POS2, NEG1, NEG2
      CHARACTER (LEN=3) :: CI, C3

!  LOCAL MATRICES AND VECTORS

      DOUBLE PRECISION :: &
         T_C_MINUS ( MAXSTREAMS, MAXLAYERS, 0:MAX_THERMAL_COEFFS ), &
         T_C_PLUS  ( MAXSTREAMS, MAXLAYERS, 0:MAX_THERMAL_COEFFS )

      DOUBLE PRECISION :: HVEC1 ( MAXSTREAMS )
      DOUBLE PRECISION :: HVEC2 ( MAXSTREAMS )
      DOUBLE PRECISION :: JVEC1 ( MAXSTREAMS )
      DOUBLE PRECISION :: TVEC1 ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION :: TVEC2 ( MAXSTREAMS_2, MAXLAYERS )

      DOUBLE PRECISION :: T_HELP1 ( 0:MAXMOMENTS )
      DOUBLE PRECISION :: T_HELP2 ( 0:MAXMOMENTS )

      DOUBLE PRECISION :: TMAT ( MAXSTREAMS, MAXSTREAMS )
      INTEGER          :: TPIVOT ( MAXSTREAMS )

! -------------------------------
!  ZERO THE BOUNDARY LAYER VALUES
! -------------------------------

      DO I = 1, NSTREAMS_2
       DO N = 1, NLAYERS
        T_WUPPER(I,N) = ZERO
        T_WLOWER(I,N) = ZERO
       ENDDO
      ENDDO

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
       END DO

!  OFFGRID: ONLY FOR QUADRATURE OR MEAN-VALUE OUTPUT

       IF ( DO_QUAD_OUTPUT .OR. DO_ADDITIONAL_MVOUT ) THEN
        IF ( N_PARTLAYERS .GT. 0 ) THEN
         DO UT = 1, N_PARTLAYERS

!  REGULAR OFF-GRID SOLUTION

          N  = PARTLAYERS_LAYERIDX(UT)
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

!  SOURCE CONSTANTS

        TERM1 = TWO * TCOM1(N,1)
        TERM2 = TWO * TCOM1(N,2)

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
          TRACE = 'DGETRF CALL FOR H, THERMAL_CLSOLUTION, LAYER '//C3
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  H VECTOR_1 AND SOLUTION BY BACK-SUBSTITUTION

        DO I = 1, NSTREAMS
          HVEC1(I) = TERM1
        ENDDO

        CALL DGETRS  ('N',NSTREAMS,1,TMAT,MAXSTREAMS,TPIVOT, &
           HVEC1,MAXSTREAMS,INFO)

        IF ( INFO .NE. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(C3, '(I3)' ) N
          MESSAGE = 'ARGUMENT I ILLEGAL VALUE, FOR I = '//CI
          TRACE = 'DGETRS CALL FOR H_1, THERMAL_CLSOLUTION, LAYER '//C3
          STATUS  = VLIDORT_SERIOUS
          RETURN
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
          TRACE = 'DGETRS CALL FOR H_2, THERMAL_CLSOLUTION, LAYER '//C3
          STATUS  = VLIDORT_SERIOUS
          RETURN
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
          TRACE = 'DGETRF CALL FOR J, THERMAL_CLSOLUTION, LAYER '//C3
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
          TRACE = 'DGETRS CALL FOR J, THERMAL_CLSOLUTION, LAYER '//C3
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  SET SOLUTION
!  ============

!  EXPANSION

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          TVEC1(I,N)  = HALF * (HVEC1(I) + JVEC1(I))
          TVEC1(I1,N) = HALF * (HVEC1(I) - JVEC1(I))
          TVEC2(I,N)  = HALF * HVEC2(I)
          TVEC2(I1,N) = TVEC2(I,N)
        ENDDO

!  VALUES AT THE LAYER BOUNDARIES

        DO I = 1, NSTREAMS_2
          T_WUPPER(I,N) = TVEC1(I,N)
          T_WLOWER(I,N) = TVEC1(I,N) + TVEC2(I,N) * DELTAU_POWER(N,2)
!          WRITE(16,'(2I4,1P2E10.10)')I,N,TVEC1(I,N),TVEC2(I,N)
!         write(*,*)i,n,T_WUPPER(I,N),T_WLOWER(I,N)
        ENDDO
!        if ( n.eq.nlayers)pause'GRONK'

!  USER SOLUTIONS

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
          ENDDO
        ENDIF

!  DNWELLING: SUM OVER ALL HARMONIC CONTRIBUTIONS, EACH USER STREAM

        IF ( DO_DNWELLING.AND.STERM_LAYERMASK_DN(N) ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            NEG1 = ZERO
            NEG2 = ZERO
            DO L = M, NMOMENTS
              NEG1 = NEG1 + T_HELP1(L)*PI_XUM(L,UM,O1,O1)
              NEG2 = NEG2 + T_HELP2(L)*PI_XUM(L,UM,O1,O1)
            ENDDO
            U_TNEG1(UM,N) = NEG1
            U_TNEG2(UM,N) = NEG2
          ENDDO
        ENDIF

!  END LAYER LOOP

      ENDDO

!  OFFGRID: ONLY FOR QUADRATURE OR MEAN-VALUE OUTPUT

      IF ( DO_QUAD_OUTPUT .OR. DO_ADDITIONAL_MVOUT ) THEN
        IF ( N_PARTLAYERS .GT. 0 ) THEN
          DO UT = 1, N_PARTLAYERS
            N  = PARTLAYERS_LAYERIDX(UT)
            DO I = 1, NSTREAMS_2
              UT_T_PARTIC(I,UT) = &
                  TVEC1(I,N) + TVEC2(I,N) * XTAU_POWER(UT,2)
            ENDDO
          ENDDO
        ENDIF
      END IF

!  FINISH

      RETURN
      END SUBROUTINE THERMAL_CLSOLUTION

!

      SUBROUTINE THERMAL_STERMS_UP ( &
        DO_SOLAR_SOURCES, NLAYERS, &
        N_THERMAL_COEFFS, DO_THERMAL_TRANSONLY, &
        N_USER_STREAMS, LOCAL_UM_START, &
        USER_STREAMS, DO_PARTLAYERS, &
        N_PARTLAYERS, PARTLAYERS_LAYERIDX, &
        N_ALLLAYERS_UP, STERM_LAYERMASK_UP, &
        T_DELT_USERM, T_UTUP_USERM, &
        DELTAU_POWER, XTAU_POWER, &
        U_TPOS1, U_TPOS2, &
        T_DIRECT_UP, T_UT_DIRECT_UP, &
        LAYER_TSUP_UP, LAYER_TSUP_UTUP )

!  THERMAL CONTRIBUTIONS TO LAYER SOURCE TERMS (UPWELLING)

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_THERMAL_COEFFS
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
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

      DOUBLE PRECISION, INTENT (OUT) :: LAYER_TSUP_UP &
          ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: LAYER_TSUP_UTUP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS )

!  LOCAL VARIABLES

      INTEGER          :: UM, N, UT, S, NT
      DOUBLE PRECISION :: SPAR, COSMUM, FAC, SUM

!  MULTIPLIERS (WHOLE LAYER)

      DOUBLE PRECISION :: T_MULT_UP ( MAXLAYERS, 0:MAX_THERMAL_COEFFS )

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

       DO N = N_ALLLAYERS_UP, NLAYERS
         LAYER_TSUP_UP(UM,N) = FAC * T_DIRECT_UP(UM,N)
       ENDDO
       IF ( DO_PARTLAYERS ) THEN
        DO UT = 1, N_PARTLAYERS
         LAYER_TSUP_UTUP(UM,UT) = FAC*T_UT_DIRECT_UP(UM,UT)
        ENDDO
       ENDIF

!  FINISH IF TRANSMITTANCE ONLY

       IF ( DO_THERMAL_TRANSONLY ) GO TO 678

!  CLASSICAL SECTION PARTICULAR INTEGRAL
!  -------------------------------------

!  WHOLE LAYER SOURCE TERMS
!   NOTE: T_DELT_USERM(N,UM) WAS INDEXED OPPOSITELY

       DO N = 1, NLAYERS
         IF ( STERM_LAYERMASK_UP(N) ) THEN
          T_MULT_UP(N,2) = U_TPOS2(UM,N)
          T_MULT_UP(N,1) = U_TPOS1(UM,N) + COSMUM * T_MULT_UP(N,2)
          SUM = T_MULT_UP(N,1)
          DO S = 2, NT
           SUM = SUM + T_MULT_UP(N,S) * DELTAU_POWER(N,S)
          ENDDO
          T_MULT_UP(N,0)   = - SUM
          SPAR = T_MULT_UP(N,0) * T_DELT_USERM(N,UM) &
               + T_MULT_UP(N,1)
          LAYER_TSUP_UP(UM,N) = LAYER_TSUP_UP(UM,N) + SPAR * FAC
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
         ENDDO
       ENDIF

!  CONTINUATION POINT FOR AVOIDING SCATTERING CALCULATIONS

 678   CONTINUE

!  END USER-STREAM LOOP

      ENDDO

!  FINISH

      RETURN
      END SUBROUTINE THERMAL_STERMS_UP

!

      SUBROUTINE THERMAL_STERMS_DN ( &
        DO_SOLAR_SOURCES, NLAYERS, &
        N_THERMAL_COEFFS, DO_THERMAL_TRANSONLY, &
        N_USER_STREAMS, LOCAL_UM_START, &
        USER_STREAMS, DO_PARTLAYERS, &
        N_PARTLAYERS, PARTLAYERS_LAYERIDX, &
        N_ALLLAYERS_DN, STERM_LAYERMASK_DN, &
        T_DELT_USERM, T_UTDN_USERM, &
        DELTAU_POWER, XTAU_POWER, &
        U_TNEG1, U_TNEG2, &
        T_DIRECT_DN, T_UT_DIRECT_DN, &
        LAYER_TSUP_DN, LAYER_TSUP_UTDN )

!  THERMAL CONTRIBUTIONS TO LAYER SOURCE TERMS (DOWNWELLING)

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_THERMAL_COEFFS
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      DOUBLE PRECISION, INTENT (IN) :: USER_STREAMS ( MAX_USER_STREAMS )
      LOGICAL, INTENT (IN) ::          DO_PARTLAYERS
      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      INTEGER, INTENT (IN) ::          N_ALLLAYERS_DN
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_DN ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
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

      DOUBLE PRECISION, INTENT (OUT) :: LAYER_TSUP_DN &
          ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: LAYER_TSUP_UTDN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS )

!  LOCAL VARIABLES

      INTEGER          :: UM, N, UT, S, NT
      DOUBLE PRECISION :: SPAR, COSMUM, FAC, SUM

!  MULTIPLIERS (WHOLE LAYER)

      DOUBLE PRECISION :: T_MULT_DN ( MAXLAYERS, 0:MAX_THERMAL_COEFFS )

!  PARTICULAR SOLUTION LAYER SOURCE TERMS ( GREEN'S FUNCTION SOLUTION )
!  --------------------------------------------------------------------

!  INITIAL MODULUS = 4.PI IF SOLAR SOURCES ARE INCLUDED

      FAC = ONE
      IF ( DO_SOLAR_SOURCES ) FAC = PI4

!  SHORTHAND

      NT = N_THERMAL_COEFFS

!  START USER ANGLE LOOP

      DO UM = LOCAL_UM_START, N_USER_STREAMS

!  DIRECT TERMS TO START

       DO N = 1, N_ALLLAYERS_DN
         LAYER_TSUP_DN(UM,N) = FAC * T_DIRECT_DN(UM,N)
       ENDDO
       IF ( DO_PARTLAYERS ) THEN
        DO UT = 1, N_PARTLAYERS
         LAYER_TSUP_UTDN(UM,UT) = FAC*T_UT_DIRECT_DN(UM,UT)
        ENDDO
       ENDIF

!  FINISH IF TRANSMITTANCE ONLY

       IF ( DO_THERMAL_TRANSONLY ) GO TO 678

!  LOCAL COSINE

       COSMUM = USER_STREAMS(UM)

!  CLASSICAL SECTION PARTICULAR INTEGRAL

!  WHOLE LAYER SOURCE TERMS
!   NOTE: T_DELT_USERM(N,UM) WAS INDEXED OPPOSITELY

       DO N = 1, NLAYERS
         IF ( STERM_LAYERMASK_DN(N) ) THEN
          T_MULT_DN(N,2) = U_TNEG2(UM,N)
          T_MULT_DN(N,1) = U_TNEG1(UM,N) - COSMUM * T_MULT_DN(N,2)
          T_MULT_DN(N,0) = - T_MULT_DN(N,1)
          SPAR = T_MULT_DN(N,0) * T_DELT_USERM(N,UM)
          SPAR = SPAR + T_MULT_DN(N,1)
          SPAR = SPAR + T_MULT_DN(N,2) * DELTAU_POWER(N,2)
          LAYER_TSUP_DN(UM,N) = LAYER_TSUP_DN(UM,N) + SPAR * FAC
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
         END DO
       END IF

!  CONTINUATION POINT FOR AVOIDING SCATTERING CALCULATIONS

 678   CONTINUE

!  END USER-STREAM LOOP

      ENDDO

!  FINISH

      RETURN
      END SUBROUTINE THERMAL_STERMS_DN

      END MODULE vlidort_thermalsup

