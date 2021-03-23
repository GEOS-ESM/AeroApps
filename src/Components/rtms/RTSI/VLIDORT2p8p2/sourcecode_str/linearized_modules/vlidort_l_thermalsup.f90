
! ###############################################################
! #                                                             #
! #                       VLIDORT_2p8p2                         #
! #                                                             #
! #  Vectorized LInearized Discrete Ordinate Radiative Transfer #
! #  -          --         -        -        -         -        #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Authors :     Robert. J. D. Spurr (1)                      #
! #                Matt Christi                                 #
! #                                                             #
! #  Address (1) : RT Solutions, inc.                           #
! #                9 Channing Street                            #
! #                Cambridge, MA 02138, USA                     #
! #                                                             #
! #  Tel:          (617) 492 1183                               #
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  This Version :   VLIDORT_2p8p2                             #
! #  Release Date :   15 April 2020                             #
! #                                                             #
! #  Previous VLIDORT Versions under Standard GPL 3.0:          #
! #  ------------------------------------------------           #
! #                                                             #
! #      2.7   F90, released August 2014                        #
! #      2.8   F90, released May    2017                        #
! #      2.8.1 F90, released August 2019                        # 
! #                                                             #
! #  Features Summary of Recent VLIDORT Versions:               #
! #  -------------------------------------------                #
! #                                                             #
! #      NEW: TOTAL COLUMN JACOBIANS         (2.4)              #
! #      NEW: BPDF Land-surface KERNELS      (2.4R)             #
! #      NEW: Thermal Emission Treatment     (2.4RT)            #
! #      Consolidated BRDF treatment         (2.4RTC)           #
! #      f77/f90 Release                     (2.5)              #
! #      External SS / New I/O Structures    (2.6)              #
! #                                                             #
! #      SURFACE-LEAVING / BRDF-SCALING      (2.7)              #
! #      TAYLOR Series / OMP THREADSAFE      (2.7)              #
! #      New Water-Leaving Treatment         (2.8)              #
! #      LBBF & BRDF-Telescoping, enabled    (2.8)              #
! #      Several Performance Enhancements    (2.8)              #
! #      Water-leaving coupled code          (2.8.1)            #
! #      Planetary problem, media properties (2.8.1)            #
! #                                                             #
! #  Features Summary of This VLIDORT Version                   #
! #  ----------------------------------------                   #
! #                                                             #
! #   2.8.2, released 15 April 2020.                            #
! #     ==> Geometry (FO/MS), check/derive separation           #
! #     ==> New setup_master for Geometry/Check/Derive          #
! #     ==> Reduction of zeroing, some dynamic memory           #
! #     ==> Use of F-matrixes only in FO code                   #
! #     ==> Use I/O type structures directly                    #
! #     ==> Doublet geometry post-processing option             #
! #                                                             #
! ###############################################################

! ###################################################################
! #                                                                 #
! # This is Version 2.8.2 of the VLIDORT_2p8 software library.      #
! # This library comes with the Standard GNU General Public License,#
! # Version 3.0, 29 June 2007. Please read this license carefully.  #
! #                                                                 #
! #      VLIDORT Copyright (c) 2003-2020.                           #
! #          Robert Spurr, RT Solutions, Inc.                       #
! #          9 Channing Street, Cambridge, MA 02138, USA.           #
! #                                                                 #
! # This file is part of VLIDORT_2p8p2 ( Version 2.8.2 )            #
! #                                                                 #
! # VLIDORT_2p8p2 is free software: you can redistribute it         #
! # and/or modify it under the terms of the Standard GNU GPL        #
! # (General Public License) as published by the Free Software      #
! # Foundation, either version 3.0 of the License, or any           #
! # later version.                                                  #
! #                                                                 #
! # VLIDORT_2p8p2 is distributed in the hope that it will be        #
! # useful, but WITHOUT ANY WARRANTY; without even the implied      #
! # warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR         #
! # PURPOSE. See the Standard GNU General Public License (GPL)      #
! # for more details.                                               #
! #                                                                 #
! # You should have received a copy of the Standard GNU General     #
! # Public License (GPL) Version 3.0, along with the VLIDORT_2p8p2  #
! # code package. If not, see <http://www.gnu.org/licenses/>.       #
! #                                                                 #
! ###################################################################

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
! #                                                             #
! ###############################################################

!  4/15/20. Version 2.8.2. No Changes.

!  Note that the routine THERMAL_LTE_LINEARIZATION, which was 
!   introduced 9/11/09 into Version 2.5, and disabled in Version
!   2.7 by the "LBBF" code, has now been removed altogether from
!   its place in the present module "vlidort_l_thermalsup_m".

!  Note also that the LBBF code introduced into Version 2.7 was
!  never functioning, but has now been debugged and enabled for Version 2.8.

      MODULE vlidort_l_thermalsup_m

      PUBLIC :: THERMAL_SETUP_PLUS,       &
                THERMAL_CLSOLUTION_PLUS,  &
                THERMAL_STERMS_UP_PLUS,   &
                THERMAL_STERMS_DN_PLUS

!                THERMAL_LTE_LINEARIZATION

      CONTAINS

      SUBROUTINE THERMAL_SETUP_PLUS ( &
        DO_USER_STREAMS, DO_UPWELLING, DO_DNWELLING,                         & ! Flags
        DO_PARTLAYERS, DO_THERMAL_TRANSONLY, DO_MSMODE_THERMAL,              & ! Flags
        DO_ATMOS_LINEARIZATION, LAYER_VARY_FLAG, LAYER_VARY_NUMBER,          & ! Linearization control
        NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_STREAMS,             & ! Numbers basic
        THERMAL_BB_INPUT, USER_STREAMS,                                      & ! thermal input, streams
        PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,         & ! Level control
        OMEGA_TOTAL, DELTAU_VERT, PARTAU_VERT, L_OMEGA_TOTAL, L_DELTAU_VERT, & ! Input optical+linearized
        T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM,                            & ! Input transmittances
        L_T_DELT_USERM, L_T_UTDN_USERM, L_T_UTUP_USERM,                      & ! Input linearized transmittances
        THERMCOEFFS, DELTAU_POWER, XTAU_POWER, TCOM1,                        & ! output thermal setups
        T_DIRECT_UP, T_DIRECT_DN, T_UT_DIRECT_UP, T_UT_DIRECT_DN,            & ! output thermal direct solutions
        L_THERMCOEFFS, L_DELTAU_POWER, L_XTAU_POWER, L_TCOM1,                & ! output Linearized thermal setups
        L_T_DIRECT_UP, L_T_DIRECT_DN, L_T_UT_DIRECT_UP, L_T_UT_DIRECT_DN )     ! output Linearized thermal direct solutions

!  SET-UP OF THERMAL EXPANSION COEFFICIENTS, ALWAYS DONE AFTER DELTA-M.
!  SET-UP OF LINEARIZATION OF THERMAL EXPANSION COEFFICIENTS.

!  LINEARIZATION W.R.T PROFILE VARIABLES, NOT BB INPUT.

      USE VLIDORT_PARS_m, Only : MAXLAYERS, MAX_PARTLAYERS, MAX_THERMAL_COEFFS, &
                                 MAX_USER_STREAMS, MAX_USER_LEVELS, MAX_ATMOSWFS, ZERO, ONE

      IMPLICIT NONE

!  flags

      LOGICAL, INTENT (IN) ::          DO_USER_STREAMS
      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_DNWELLING

      LOGICAL, INTENT (IN) ::          DO_PARTLAYERS
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_MSMODE_THERMAL

!  Linearization control

      LOGICAL, INTENT (IN) ::          DO_ATMOS_LINEARIZATION
      LOGICAL, INTENT (IN) ::          LAYER_VARY_FLAG   ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          LAYER_VARY_NUMBER ( MAXLAYERS )

!  Numbers

      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_THERMAL_COEFFS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          N_PARTLAYERS

!  thermal input and streams

      DOUBLE PRECISION, INTENT (IN) :: THERMAL_BB_INPUT ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: USER_STREAMS  ( MAX_USER_STREAMS )

!  level output control

      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_DN ( MAXLAYERS )

!  Optical

      DOUBLE PRECISION, INTENT (IN) :: OMEGA_TOTAL ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PARTAU_VERT ( MAX_PARTLAYERS )

!  Linearized optical

      DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_TOTAL ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  linearized Transmittances

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

!  outputs
!  =======

!  Auxiliary setups

      DOUBLE PRECISION, INTENT (OUT) :: THERMCOEFFS  ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (OUT) :: DELTAU_POWER ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (OUT) :: XTAU_POWER   ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (OUT) :: TCOM1 ( MAXLAYERS, MAX_THERMAL_COEFFS )

!  thermal direct solutions

      DOUBLE PRECISION, INTENT (OUT) :: T_DIRECT_UP    ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: T_DIRECT_DN    ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: T_UT_DIRECT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: T_UT_DIRECT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS )

!  Linearized Auxiliary setups

      DOUBLE PRECISION, INTENT (OUT) :: L_THERMCOEFFS  ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_DELTAU_POWER ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_XTAU_POWER   ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_TCOM1        ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )

!  Linearized thermal direct solutions

      DOUBLE PRECISION, INTENT (OUT) :: L_T_DIRECT_UP    ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_T_DIRECT_DN    ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_T_UT_DIRECT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS)
      DOUBLE PRECISION, INTENT (OUT) :: L_T_UT_DIRECT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  LOCAL VARIABLES
!  ---------------

      INTEGER          :: N, S, UT, UM, Q, NPARS, NT
      DOUBLE PRECISION :: XTAU, SUM, COSMUM, L_D, L_X
      DOUBLE PRECISION :: OMEGAS1 ( MAXLAYERS )
      DOUBLE PRECISION :: HELP(MAXLAYERS)

!  MULTIPLIERS (WHOLE LAYER), and their linearizations

      DOUBLE PRECISION :: T_MULT_UP   (MAXLAYERS,0:MAX_THERMAL_COEFFS)
      DOUBLE PRECISION :: T_MULT_DN   (MAXLAYERS,0:MAX_THERMAL_COEFFS)
      DOUBLE PRECISION :: L_T_MULT_UP (MAXLAYERS,0:MAX_THERMAL_COEFFS,MAX_ATMOSWFS)
      DOUBLE PRECISION :: L_T_MULT_DN (MAXLAYERS,0:MAX_THERMAL_COEFFS,MAX_ATMOSWFS)

!Rob fix 7/8/2016 - initi alized for packing
      T_DIRECT_UP   = ZERO ; T_UT_DIRECT_UP   = ZERO
      T_DIRECT_DN   = ZERO ; T_UT_DIRECT_DN   = ZERO
      L_T_DIRECT_UP = ZERO ; L_T_UT_DIRECT_UP = ZERO
      L_T_DIRECT_DN = ZERO ; L_T_UT_DIRECT_DN = ZERO

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
        END DO
        IF ( DO_ATMOS_LINEARIZATION ) THEN
          DO UT = 1, N_PARTLAYERS
            N = PARTLAYERS_LAYERIDX(UT)
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_X = L_DELTAU_VERT(Q,N) *  PARTAU_VERT(UT)
              L_XTAU_POWER(UT,1,Q) = ZERO
              DO S = 2, N_THERMAL_COEFFS
                L_XTAU_POWER(UT,S,Q) = DBLE(S-1) * L_X * XTAU_POWER(UT,S-1)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  INITIAL SET OF COEFFICIENTS

      DO N = 1, NLAYERS
        THERMCOEFFS(N,1) = THERMAL_BB_INPUT(N-1)
        HELP(N) = (THERMAL_BB_INPUT(N)-THERMAL_BB_INPUT(N-1)) / DELTAU_VERT(N)
      END DO

!  PIECEWISE CONTINUOUS FOR LINEAR REGIME

      IF ( N_THERMAL_COEFFS == 2 ) THEN
        DO N = 1, NLAYERS
          THERMCOEFFS(N,2) = HELP(N)
        END DO
      END IF

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

!  LINEARIZED COEFFICIENTS AND AUXILIARY QUANTITIES
!  ------------------------------------------------

!  Transmittance only

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
      ENDIF

!  General thermal solution

      IF ( .not.DO_THERMAL_TRANSONLY ) THEN
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
!   (Zeroing is done at the beginning of the routine now)

      IF ( DO_MSMODE_THERMAL ) RETURN

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
          T_DIRECT_UP(UM,N) = T_MULT_UP(N,0) * T_DELT_USERM(N,UM) + T_MULT_UP(N,1)
         ENDIF

!  LINEARIZATION OF DIRECT SOLUTION: WHOLE LAYER SOURCE TERMS

         IF ( DO_ATMOS_LINEARIZATION ) THEN
          IF ( STERM_LAYERMASK_UP(N).AND.LAYER_VARY_FLAG(N) ) THEN
           NPARS = LAYER_VARY_NUMBER(N)
           DO Q = 1, NPARS
            L_T_MULT_UP(N,NT,Q) = L_TCOM1(N,NT,Q)
            DO S = N_THERMAL_COEFFS - 1, 1, -1
             L_T_MULT_UP(N,S,Q) = L_TCOM1(N,S,Q) + S * COSMUM * L_T_MULT_UP(N,S+1,Q)
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
        DO_UPWELLING, DO_DNWELLING, DO_USER_STREAMS,              & ! input flags
        DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT, DO_THERMAL_TRANSONLY, & ! Input flags
        DO_ATMOS_WFS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER,     & ! Input linearization control
        NSTREAMS, NLAYERS, N_THERMAL_COEFFS, N_USER_STREAMS,  & ! Input basic numbers
        NMOMENTS, NSTREAMS_2, LOCAL_UM_START, N_PARTLAYERS,   & ! Input numbers
        PARTLAYERS_LAYERIDX, STERM_LMASK_UP, STERM_LMASK_DN,  & ! Input Level control
        QUAD_STREAMS, QUAD_HALFWTS, OMEGA_GREEK, SAB, DAB,    & ! Input optical and SAB/DAB
        T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,       & ! Input Discrete Ord. Trans.
        PI_XQP, PI_XUP, PI_XUM, PI_XQM_PRE,                   & ! Input PI matrices, 
        THERMCOEFFS, DELTAU_POWER, XTAU_POWER, TCOM1,         & ! Input thermal setups
        L_OMEGA_GREEK, L_SAB, L_DAB,                          & ! Input Linearized optical and SAB/DAB
        L_T_DELT_DISORDS, L_T_DISORDS_UTDN, L_T_DISORDS_UTUP, & ! Input Linearized Discrete Ord. Trans.
        L_THERMCOEFFS, L_DELTAU_POWER, L_XTAU_POWER, L_TCOM1, & ! Input Linearized thermal setups
        T_WUPPER, T_WLOWER, UT_T_PARTIC,                      & ! Output thermal solutions
        U_TPOS1, U_TNEG1, U_TPOS2, U_TNEG2,                   & ! Output User thermal solutions
        L_T_WUPPER, L_T_WLOWER, L_UT_T_PARTIC,                & ! Output Linearized thermal solutions
        L_U_TPOS1, L_U_TNEG1, L_U_TPOS2, L_U_TNEG2,           & ! Output Linearized User thermal solutions
        STATUS, MESSAGE, TRACE )                                ! Exception handling

!  CLASSICAL THERMAL PARTICULAR INTEGRAL, ALL LAYERS + LINEARIZATIONS
!  USES COEFFICIENT EXPANSION OF ATTENUATION.

      USE VLIDORT_PARS_m, Only : MAXSTOKES, MAXSTREAMS, MAXLAYERS, MAX_PARTLAYERS, &
                                 MAXMOMENTS, MAX_USER_STREAMS, MAX_USER_LEVELS,    &
                                 MAX_THERMAL_COEFFS, MAXSTREAMS_2, MAX_ATMOSWFS,   &
                                 VLIDORT_SUCCESS, VLIDORT_SERIOUS, ZERO, HALF, TWO

      USE LAPACK_TOOLS_m, Only : DGETRF, DGETRS

      IMPLICIT NONE


!  Flags 
!mick fix 3/30/2015 - added DO_MVOUT_ONLY to input
!mick fix 9/19/2017 - added DO_USER_STREAMS to input

      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_DNWELLING
      LOGICAL, INTENT (IN) ::          DO_USER_STREAMS
      LOGICAL, INTENT (IN) ::          DO_MVOUT_ONLY 

!  removed 7/6/16 version 2.8
!      LOGICAL, INTENT (IN) ::          DO_QUAD_OUTPUT

      LOGICAL, INTENT (IN) ::          DO_ADDITIONAL_MVOUT
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY

!  Linearization control

      LOGICAL, INTENT (IN) ::          DO_ATMOS_WFS
      LOGICAL, INTENT (IN) ::          LAYER_VARY_FLAG   ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          LAYER_VARY_NUMBER ( MAXLAYERS )

!  basic numbers

      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_THERMAL_COEFFS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS

!  other numbers


      INTEGER, INTENT (IN) ::          NMOMENTS
      INTEGER, INTENT (IN) ::          NSTREAMS_2
      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      INTEGER, INTENT (IN) ::          N_PARTLAYERS

!  Level control

      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      LOGICAL, INTENT (IN) ::          STERM_LMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::          STERM_LMASK_DN ( MAXLAYERS )

!  Quadrature and Optical, SAB, DAB

      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_HALFWTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: SAB ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DAB ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )

!  discrete ordinate transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTDN ( MAXSTREAMS, MAX_PARTLAYERS )

!  PI matrices

      DOUBLE PRECISION, INTENT (IN) :: PI_XQP     ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_XUP     ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_XUM     ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_XQM_PRE ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )

!  thermal setups

      DOUBLE PRECISION, INTENT (IN) :: THERMCOEFFS  ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_POWER ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (IN) :: XTAU_POWER   ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (IN) :: TCOM1        ( MAXLAYERS, MAX_THERMAL_COEFFS )

!  Linearized Optical, SAB, DAB

      DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SAB &
                    ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES,  MAXLAYERS,  MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_DAB &
                    ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES,  MAXLAYERS,  MAX_ATMOSWFS )

!  Linearized thermal setups

      DOUBLE PRECISION, INTENT (IN) :: L_THERMCOEFFS  ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_POWER ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_XTAU_POWER   ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_TCOM1        ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )

!  Linearized discrete ordinate transmittances

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DISORDS_UTDN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  OUTPUT
!  ======

!  discrete ordinate thermal solutions

      DOUBLE PRECISION, INTENT (OUT) ::  T_WUPPER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) ::  T_WLOWER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) ::  UT_T_PARTIC( MAXSTREAMS_2, MAX_PARTLAYERS )

!  user-stream thermal solutions

      DOUBLE PRECISION, INTENT (OUT) ::  U_TPOS1 ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) ::  U_TNEG1 ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) ::  U_TPOS2 ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) ::  U_TNEG2 ( MAX_USER_STREAMS, MAXLAYERS )

!  Linearized discrete ordinate thermal solutions

      DOUBLE PRECISION, INTENT (OUT) ::  L_T_WUPPER     ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) ::  L_T_WLOWER     ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS)
      DOUBLE PRECISION, INTENT (OUT) ::  L_UT_T_PARTIC  ( MAXSTREAMS_2, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Linearized user-stream thermal solutions

      DOUBLE PRECISION, INTENT (OUT) ::  L_U_TPOS1 ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) ::  L_U_TNEG1 ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) ::  L_U_TPOS2 ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) ::  L_U_TNEG2 ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )

!  exception handling

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

      DOUBLE PRECISION :: T_C_MINUS ( MAXSTREAMS, MAXLAYERS, 0:MAX_THERMAL_COEFFS )
      DOUBLE PRECISION :: T_C_PLUS  ( MAXSTREAMS, MAXLAYERS, 0:MAX_THERMAL_COEFFS )

      DOUBLE PRECISION :: L_T_C_MINUS ( MAXSTREAMS, MAXLAYERS, 0:MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_C_PLUS  ( MAXSTREAMS, MAXLAYERS, 0:MAX_THERMAL_COEFFS, MAX_ATMOSWFS )

! -------------------------------
!  ZERO THE BOUNDARY LAYER VALUES
! -------------------------------

      DO I = 1, NSTREAMS_2
       DO N = 1, NLAYERS
        T_WUPPER(I,N) = ZERO
        T_WLOWER(I,N) = ZERO
       ENDDO
      ENDDO

      IF ( DO_ATMOS_WFS ) THEN
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

        IF ( DO_ATMOS_WFS ) THEN
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

! removed Quad_output flag, Version 2.8, 7/8/16
!  OFFGRID: ONLY FOR QUADRATURE OR MEAN-VALUE OUTPUT
!       IF ( DO_QUAD_OUTPUT .OR. DO_MVOUT_ONLY .OR. DO_ADDITIONAL_MVOUT ) THEN

       IF ( DO_MVOUT_ONLY .OR. DO_ADDITIONAL_MVOUT ) THEN
        IF ( N_PARTLAYERS .GT. 0 ) THEN
         DO UT = 1, N_PARTLAYERS
          N  = PARTLAYERS_LAYERIDX(UT)

!  REGULAR OFF-GRID SOLUTION

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

!  LINEARIZED OFF-GRID SOLUTION

          IF ( DO_ATMOS_WFS ) THEN
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

        DO_VARY = DO_ATMOS_WFS.AND.LAYER_VARY_FLAG(N)
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
          TRACE   = 'DGETRF CALL FOR H, THERMAL_CLSOLUTION_PLUS, LAYER '//C3
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  H VECTOR_1 AND SOLUTION BY BACK-SUBSTITUTION

        DO I = 1, NSTREAMS
          HVEC1(I) = TERM1
        ENDDO
        CALL DGETRS ('N',NSTREAMS,1,TMAT,MAXSTREAMS,TPIVOT, HVEC1,MAXSTREAMS,INFO)

        IF ( INFO .NE. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(C3, '(I3)' ) N
          MESSAGE = 'ARGUMENT I ILLEGAL VALUE, FOR I = '//CI
          TRACE   = 'DGETRS CALL FOR H_1, THERMAL_CLSOLUTION_PLUS, LAYER '//C3
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  LINEARIZATION H VECTOR_1

        IF ( DO_VARY ) THEN
          DO Q = 1, NVARY
            DO I = 1, NSTREAMS
              SUM = ZERO
              DO J = 1, NSTREAMS
                SUM = SUM - L_SAB(I,J,O1,O1,N,Q) * QUAD_STREAMS(I) * HVEC1(J)
              ENDDO
              LHVEC1(I,Q) = L_TERM1(Q) + SUM
            ENDDO
          ENDDO
          CALL DGETRS ('N',NSTREAMS,NVARY,TMAT,MAXSTREAMS,TPIVOT, LHVEC1,MAXSTREAMS,INFO)
          IF ( INFO .NE. 0 ) THEN
            WRITE(CI, '(I3)' ) INFO
            WRITE(C3, '(I3)' ) N
            MESSAGE = 'ARGUMENT I ILLEGAL VALUE, FOR I = '//CI
            TRACE   = 'DGETRS CALL LH_1, THERMAL_CLSOLUTION_PLUS, LAYER '//C3
            STATUS  = VLIDORT_SERIOUS
            RETURN
          ENDIF

        ENDIF

!  H VECTOR_2 AND SOLUTION BY BACK-SUBSTITUTION

        DO I = 1, NSTREAMS
          HVEC2(I) = TERM2
        ENDDO
        CALL DGETRS ('N',NSTREAMS,1,TMAT,MAXSTREAMS,TPIVOT, HVEC2,MAXSTREAMS,INFO)

        IF ( INFO .NE. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(C3, '(I3)' ) N
          MESSAGE = 'ARGUMENT I ILLEGAL VALUE, FOR I = '//CI
          TRACE   = 'DGETRS CALL FOR H_2, THERMAL_CLSOLUTION_PLUS, LAYER '//C3
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  LINEARIZATION H VECTOR_2

        IF ( DO_VARY ) THEN
          DO Q = 1, NVARY
            DO I = 1, NSTREAMS
              SUM = ZERO
              DO J = 1, NSTREAMS
                SUM = SUM - L_SAB(I,J,O1,O1,N,Q) * QUAD_STREAMS(I) * HVEC2(J)
              ENDDO
              LHVEC2(I,Q) = L_TERM2(Q) + SUM
            ENDDO
          ENDDO
          CALL DGETRS ('N',NSTREAMS,NVARY,TMAT,MAXSTREAMS,TPIVOT, LHVEC2,MAXSTREAMS,INFO)

          IF ( INFO .NE. 0 ) THEN
            WRITE(CI, '(I3)' ) INFO
            WRITE(C3, '(I3)' ) N
            MESSAGE = 'ARGUMENT I ILLEGAL VALUE, FOR I = '//CI
            TRACE   = 'DGETRS CALL LIN H_2, THERMAL_CLSOLUTION_PLUS, LAYER '//C3
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
          TRACE   = 'DGETRF CALL FOR J, THERMAL_CLSOLUTION_PLUS, LAYER '//C3
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  J VECTOR_1 SOLUTION BY BACK-SUBSTITUTION

        CALL DGETRS ('N',NSTREAMS,1,TMAT,MAXSTREAMS,TPIVOT,JVEC1,MAXSTREAMS,INFO)

        IF ( INFO .NE. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(C3, '(I3)' ) N
          MESSAGE = 'ARGUMENT I ILLEGAL VALUE, FOR I = '//CI
          TRACE   = 'DGETRS CALL FOR J, THERMAL_CLSOLUTION_PLUS, LAYER '//C3
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

          CALL DGETRS ('N',NSTREAMS,NVARY,TMAT,MAXSTREAMS,TPIVOT,LJVEC1,MAXSTREAMS,INFO)

          IF ( INFO .NE. 0 ) THEN
            WRITE(CI, '(I3)' ) INFO
            WRITE(C3, '(I3)' ) N
            MESSAGE = 'ARGUMENT I ILLEGAL VALUE, FOR I = '//CI
            TRACE   = 'DGETRS CALL FOR LIN J, THERMAL_CLSOLUTION_PLUS, LAYER '//C3
            STATUS  = VLIDORT_SERIOUS
            RETURN
          ENDIF

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

!mick fix 9/19/2017 - added IF condition
        IF ( DO_USER_STREAMS ) THEN

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
                L_T_HELP1(L,Q) = L_SUM1 *   OMEGA_GREEK(L,N,O1,O1) + SUM1 * L_OMEGA_GREEK(L,N,O1,O1,Q)
                L_T_HELP2(L,Q) = L_SUM2 *   OMEGA_GREEK(L,N,O1,O1) + SUM2 * L_OMEGA_GREEK(L,N,O1,O1,Q)
              ENDDO
            ENDIF
          ENDDO

!  UPWELLING: SUM OVER ALL HARMONIC CONTRIBUTIONS, EACH USER STREAM

          IF ( DO_UPWELLING.AND.STERM_LMASK_UP(N) ) THEN
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

          IF ( DO_DNWELLING.AND.STERM_LMASK_DN(N) ) THEN
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

        ENDIF

!  END LAYER LOOP

      ENDDO

!  OFFGRID: ONLY FOR QUADRATURE OR MEAN-VALUE OUTPUT
!  =================================================

! removed Quad_output flag, Version 2.8, 7/6/16
!       IF ( DO_QUAD_OUTPUT .OR. DO_MVOUT_ONLY .OR. DO_ADDITIONAL_MVOUT ) THEN

      IF ( DO_MVOUT_ONLY .OR. DO_ADDITIONAL_MVOUT ) THEN
        IF ( N_PARTLAYERS .GT. 0 ) THEN
          DO UT = 1, N_PARTLAYERS
            N  = PARTLAYERS_LAYERIDX(UT)
            DO I = 1, NSTREAMS_2
              UT_T_PARTIC(I,UT) = TVEC1(I,N) + TVEC2(I,N) * XTAU_POWER(UT,2)
            ENDDO
            IF ( LAYER_VARY_FLAG(N) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                DO I = 1, NSTREAMS_2
                  L_UT_T_PARTIC(I,UT,Q) = L_TVEC1(I,N,Q) + L_TVEC2(I,N,Q) *   XTAU_POWER(UT,2) &
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
        DO_SOLAR_SOURCES, DO_THERMAL_TRANSONLY, DO_PARTLAYERS,               & ! Input flags
        DO_ATMOS_LINEARIZATION, LAYER_VARY_FLAG, LAYER_VARY_NUMBER,          & ! Input Linearization control
        NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_STREAMS,             & ! Input numbers
        N_ALLLAYERS_UP, PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP,             & ! Input level control
        LOCAL_UM_START, USER_STREAMS,                                        & ! Input User streams
        T_DELT_USERM, T_UTUP_USERM, L_T_DELT_USERM, L_T_UTUP_USERM,          & ! Input transmittances
        DELTAU_POWER, XTAU_POWER, L_DELTAU_POWER, L_XTAU_POWER,              & ! Input Thermal setups
        U_TPOS1, U_TPOS2, L_U_TPOS1, L_U_TPOS2,                              & ! Input Thermal solutions               
        T_DIRECT_UP, T_UT_DIRECT_UP, L_T_DIRECT_UP, L_T_UT_DIRECT_UP,        & ! Input thermal direct solutions
        LAYER_TSUP_UP, LAYER_TSUP_UTUP, L_LAYER_TSUP_UP, L_LAYER_TSUP_UTUP )   ! Output user RTE thermal

!  THERMAL CONTRIBUTIONS TO LAYER SOURCE TERMS (UPWELLING)
!  LINEARIZED THERMAL CONTRIBUTIONS TO LAYER SOURCE TERMS (UPWELLING)

      USE VLIDORT_PARS_m, Only : MAXLAYERS, MAX_PARTLAYERS, MAX_USER_STREAMS, &
                                 MAX_USER_LEVELS, MAX_ATMOSWFS, MAX_THERMAL_COEFFS, ONE, PI4

      IMPLICIT NONE

!  Inputs
!  ======

!  flags

      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_PARTLAYERS

!  Linearization control

      LOGICAL, INTENT (IN) ::          DO_ATMOS_LINEARIZATION
      LOGICAL, INTENT (IN) ::          LAYER_VARY_FLAG   ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          LAYER_VARY_NUMBER ( MAXLAYERS )

!  Numbers

      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      INTEGER, INTENT (IN) ::          N_THERMAL_COEFFS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS

!  Level control

      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      INTEGER, INTENT (IN) ::          N_ALLLAYERS_UP
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_UP ( MAXLAYERS )

!  User Streams/Transm.

      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      DOUBLE PRECISION, INTENT (IN) :: USER_STREAMS ( MAX_USER_STREAMS )

!  User Transmittances and Linearizations

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM   ( MAXLAYERS,      MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_USERM   ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

!  Thermal setups and Linearizations

      DOUBLE PRECISION, INTENT (IN) :: DELTAU_POWER   ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_POWER ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (IN) :: XTAU_POWER   ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (IN) :: L_XTAU_POWER ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )

!  User-stream thermal solutions and Linearizations

      DOUBLE PRECISION, INTENT (IN) :: U_TPOS1   ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: U_TPOS2   ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_U_TPOS1 ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_U_TPOS2 ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )

!  Thermal direct solutions and Linearizations

      DOUBLE PRECISION, INTENT (IN) :: T_DIRECT_UP    ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DIRECT_UP  ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (IN) :: T_UT_DIRECT_UP   ( MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UT_DIRECT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS)

!  OUTPUT
!  ======

!  Postprocessed thermal solutions and Linearizations

      DOUBLE PRECISION, INTENT (OUT) :: LAYER_TSUP_UP     ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: L_LAYER_TSUP_UP   ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (OUT) :: LAYER_TSUP_UTUP   ( MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: L_LAYER_TSUP_UTUP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  LOCAL VARIABLES
!  ---------------

      INTEGER ::          UM, N, UT, S, NT, Q
      DOUBLE PRECISION :: SPAR, COSMUM, FAC, SUM

!  MULTIPLIERS (WHOLE LAYER)

      DOUBLE PRECISION :: T_MULT_UP   ( MAXLAYERS, 0:MAX_THERMAL_COEFFS )
      DOUBLE PRECISION :: L_T_MULT_UP ( MAXLAYERS, 0:MAX_THERMAL_COEFFS, MAX_ATMOSWFS )

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

!  ONLY DO THE NEXT SECTION FOR SCATTERING SOLUTIONS
!   -- Version 2.8, 7/7/16 remove goto 
!       IF ( DO_THERMAL_TRANSONLY ) GO TO 678

       IF ( .not. DO_THERMAL_TRANSONLY ) THEN

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

!  End scattering clause

       ENDIF

!  CONTINUATION POINT FOR AVOIDING SCATTERING CALCULATIONS
! 678   CONTINUE

!  END USER-STREAM LOOP

      END DO

!  FINISH

      RETURN
      END SUBROUTINE THERMAL_STERMS_UP_PLUS

!

      SUBROUTINE THERMAL_STERMS_DN_PLUS ( &
        DO_SOLAR_SOURCES, DO_THERMAL_TRANSONLY, DO_PARTLAYERS,               & ! Input flags
        DO_ATMOS_LINEARIZATION, LAYER_VARY_FLAG, LAYER_VARY_NUMBER,          & ! Input Linearization control
        NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_STREAMS,             & ! Input numbers
        N_ALLLAYERS_DN, PARTLAYERS_LAYERIDX, STERM_LAYERMASK_DN,             & ! Input level control
        LOCAL_UM_START, USER_STREAMS,                                        & ! Input User streams
        T_DELT_USERM, T_UTDN_USERM, L_T_DELT_USERM, L_T_UTDN_USERM,          & ! Input transmittances
        DELTAU_POWER, XTAU_POWER, L_DELTAU_POWER, L_XTAU_POWER,              & ! Input Thermal setups
        U_TNEG1, U_TNEG2, L_U_TNEG1, L_U_TNEG2,                              & ! Input Thermal solutions               
        T_DIRECT_DN, T_UT_DIRECT_DN, L_T_DIRECT_DN, L_T_UT_DIRECT_DN,        & ! Input thermal direct solutions
        LAYER_TSUP_DN, LAYER_TSUP_UTDN, L_LAYER_TSUP_DN, L_LAYER_TSUP_UTDN )   ! Output user RTE thermal

!  THERMAL CONTRIBUTIONS TO LAYER SOURCE TERMS (DOWNWELLING)
!  LINEARIZED THERMAL CONTRIBUTIONS TO LAYER SOURCE TERMS (DOWNWELLING)

      USE VLIDORT_PARS_m, Only : MAXLAYERS, MAX_PARTLAYERS, MAX_USER_STREAMS, &
                                 MAX_USER_LEVELS, MAX_ATMOSWFS, MAX_THERMAL_COEFFS, ONE, PI4

      IMPLICIT NONE

!  Inputs
!  ======

!  flags

      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_PARTLAYERS

!  Linearization control

      LOGICAL, INTENT (IN) ::          DO_ATMOS_LINEARIZATION
      LOGICAL, INTENT (IN) ::          LAYER_VARY_FLAG   ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          LAYER_VARY_NUMBER ( MAXLAYERS )

!  Numbers

      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      INTEGER, INTENT (IN) ::          N_THERMAL_COEFFS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS

!  Level control

      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      INTEGER, INTENT (IN) ::          N_ALLLAYERS_DN
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_DN ( MAXLAYERS )

!  User Streams/Transm.

      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      DOUBLE PRECISION, INTENT (IN) :: USER_STREAMS ( MAX_USER_STREAMS )

!  User Transmittances and Linearizations

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM   ( MAXLAYERS,      MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_USERM   ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

!  Thermal setups and Linearizations

      DOUBLE PRECISION, INTENT (IN) :: DELTAU_POWER   ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_POWER ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (IN) :: XTAU_POWER   ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION, INTENT (IN) :: L_XTAU_POWER ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )

!  User-stream thermal solutions and Linearizations

      DOUBLE PRECISION, INTENT (IN) :: U_TNEG1   ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: U_TNEG2   ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_U_TNEG1 ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_U_TNEG2 ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )

!  Thermal direct solutions and Linearizations

      DOUBLE PRECISION, INTENT (IN) :: T_DIRECT_DN    ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DIRECT_DN  ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (IN) :: T_UT_DIRECT_DN   ( MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UT_DIRECT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS)

!  OUTPUT
!  ======

!  Postprocessed thermal solutions and Linearizations

      DOUBLE PRECISION, INTENT (OUT) :: LAYER_TSUP_DN     ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: L_LAYER_TSUP_DN   ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (OUT) :: LAYER_TSUP_UTDN   ( MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: L_LAYER_TSUP_UTDN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  LOCAL VARIABLES
!  ---------------

      INTEGER ::          UM, N, UT, S, NT, Q
      DOUBLE PRECISION :: SPAR, COSMUM, FAC, SUM

!  MULTIPLIERS (WHOLE LAYER)

      DOUBLE PRECISION :: T_MULT_DN   ( MAXLAYERS, 0:MAX_THERMAL_COEFFS )
      DOUBLE PRECISION :: L_T_MULT_DN ( MAXLAYERS, 0:MAX_THERMAL_COEFFS, MAX_ATMOSWFS )

!  PARTICULAR SOLUTION LAYER SOURCE TERMS
!  --------------------------------------

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

!  ONLY DO THE NEXT SECTION FOR SCATTERING SOLUTIONS
!   -- Version 2.8, 7/7/16 remove goto 
!       IF ( DO_THERMAL_TRANSONLY ) GO TO 678

       IF ( .not. DO_THERMAL_TRANSONLY ) THEN

!  CLASSICAL SECTION PARTICULAR INTEGRAL
!  =====================================

!  LOCAL COSINE

        COSMUM = USER_STREAMS(UM)

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

!  End scattering clause

       ENDIF

!  CONTINUATION POINT FOR AVOIDING SCATTERING CALCULATIONS
! 678 CONTINUE

!  END USER-STREAM LOOP

      ENDDO

!  FINISH

      RETURN
      END SUBROUTINE THERMAL_STERMS_DN_PLUS

!  Finish Module

      END MODULE vlidort_l_thermalsup_m

