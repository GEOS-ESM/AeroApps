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
! #                   2.5, 2.6, 2.7                             #
! #  Release Date :   December 2005  (2.0)                      #
! #  Release Date :   March 2007     (2.2)                      #
! #  Release Date :   October 2007   (2.3)                      #
! #  Release Date :   December 2008  (2.4)                      #
! #  Release Date :   April 2009     (2.4R)                     #
! #  Release Date :   July 2009      (2.4RT)                    #
! #  Release Date :   October 2010   (2.4RTC)                   #
! #  Release Date :   March 2011     (2.5)                      #
! #  Release Date :   May 2012       (2.6)                      #
! #  Release Date :   August 2014    (2.7)                      #
! #                                                             #
! #       NEW: TOTAL COLUMN JACOBIANS         (2.4)             #
! #       NEW: BPDF Land-surface KERNELS      (2.4R)            #
! #       NEW: Thermal Emission Treatment     (2.4RT)           #
! #       Consolidated BRDF treatment         (2.4RTC)          #
! #       f77/f90 Release                     (2.5)             #
! #       External SS / New I/O Structures    (2.6)             #
! #                                                             #
! #       SURFACE-LEAVING / BRDF-SCALING      (2.7)             #
! #       TAYLOR Series / OMP THREADSAFE      (2.7)             #
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
! #            VLIDORT_LAC_MISCSETUPS (master)                  #
! #              VLIDORT_LC_PREPTRANS                           #
! #                                                             #
! #            LC_EMULT_MASTER (master), calling:               #
! #                LC_WHOLELAYER_EMULT_UP                       #
! #                LC_WHOLELAYER_EMULT_DN                       #
! #                LC_PARTLAYER_EMULT_UP                        #
! #                LC_PARTLAYER_EMULT_DN                        #
! #                                                             #
! #            LC_EMULT_MASTER_OBSGEO (master), calling:        #
! #                LC_WHOLELAYER_EMULT_OG_UP                    #
! #                LC_WHOLELAYER_EMULT_OG_DN                    #
! #                LC_PARTLAYER_EMULT_OG_UP                     #
! #                LC_PARTLAYER_EMULT_OG_DN                     #
! #                                                             #
! ###############################################################

!  Notes for Version 2.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob  Fix 05/10/13  - Introduce TAYLOR_LIMIT parameters in module VLIDORT_PARS, instead of "HOPITAL_TOLERANCE"
!     Rob  Fix 05/10/13  - L'Hopitals Rule replaced by Taylor series (original calculation was first term in series!)
!     Rob  Fix 02/19/14  - Final  Taylor series stuff.

      MODULE vlidort_lc_miscsetups

!  @@@ Rob Fix 10 May 13 - need the Taylor series routines

      USE vlidort_Taylor_m, ONLY : TAYLOR_SERIES_L_1

      PRIVATE
      PUBLIC :: VLIDORT_LAC_MISCSETUPS,&
                LC_EMULT_MASTER,&
                LC_EMULT_MASTER_OBSGEO

      CONTAINS

      SUBROUTINE VLIDORT_LAC_MISCSETUPS ( &
        DO_DELTAM_SCALING, NSTOKES, &
        NLAYERS, OMEGA_TOTAL_INPUT, &
        DELTAU_VERT_INPUT, GREEKMAT_TOTAL_INPUT, &
        NMOMENTS, NSTOKES_SQ, NBEAMS, &
        DO_ATMOS_LINEARIZATION, LAYER_VARY_FLAG, &
        LAYER_VARY_NUMBER, L_OMEGA_TOTAL_INPUT, &
        L_DELTAU_VERT_INPUT, L_GREEKMAT_TOTAL_INPUT, &
        DELTAU_SLANT, TRUNC_FACTOR, FAC1, &
        MUELLER_INDEX, OMEGA_GREEK, &
        DO_PLANE_PARALLEL, DO_SOLUTION_SAVING, &
        NSTREAMS, LAYER_PIS_CUTOFF, QUAD_STREAMS, &
        N_USER_STREAMS, DO_USER_STREAMS, &
        USER_SECANTS, N_PARTLAYERS, &
        PARTLAYERS_LAYERIDX, &
        DELTAU_VERT, PARTAU_VERT, T_DELT_DISORDS, &
        T_DISORDS_UTUP, T_DISORDS_UTDN, &
        T_DELT_MUBAR, T_UTDN_MUBAR, &
        T_DELT_USERM, T_UTDN_USERM, &
        T_UTUP_USERM, AVERAGE_SECANT, &
        L_OMEGA_TOTAL, L_DELTAU_VERT, &
        L_GREEKMAT_TOTAL, L_DELTAU_SLANT, &
        DO_SCATMAT_VARIATION, L_TRUNC_FACTOR, &
        L_OMEGA_GREEK, LC_AVERAGE_SECANT, &
        LC_INITIAL_TRANS, L_T_DELT_DISORDS, &
        L_T_DISORDS_UTDN, L_T_DISORDS_UTUP, &
        LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR, &
        L_T_DELT_USERM, L_T_UTDN_USERM, &
        L_T_UTUP_USERM )

      USE VLIDORT_PARS
      USE VLIDORT_LA_MISCSETUPS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_DELTAM_SCALING
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NLAYERS
      DOUBLE PRECISION, INTENT (IN) :: OMEGA_TOTAL_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: GREEKMAT_TOTAL_INPUT &
          ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
      INTEGER, INTENT (IN) ::          NMOMENTS
      INTEGER, INTENT (IN) ::          NSTOKES_SQ
      INTEGER, INTENT (IN) ::          NBEAMS
      LOGICAL, INTENT (IN) ::          DO_ATMOS_LINEARIZATION
      LOGICAL, INTENT (IN) ::          LAYER_VARY_FLAG  ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          LAYER_VARY_NUMBER ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_TOTAL_INPUT &
          ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT_INPUT &
          ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_GREEKMAT_TOTAL_INPUT &
           ( MAX_ATMOSWFS, 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_SLANT &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: TRUNC_FACTOR ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: FAC1 ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: OMEGA_GREEK &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL
      LOGICAL, INTENT (IN) ::          DO_SOLUTION_SAVING
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          LAYER_PIS_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      LOGICAL, INTENT (IN) ::          DO_USER_STREAMS
      DOUBLE PRECISION, INTENT (IN) :: USER_SECANTS  ( MAX_USER_STREAMS )
      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PARTAU_VERT ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTUP &
          ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTDN &
          ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )

      DOUBLE PRECISION, INTENT (OUT) :: L_OMEGA_TOTAL &
          ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: L_DELTAU_VERT &
          ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: L_GREEKMAT_TOTAL &
          ( MAX_ATMOSWFS, 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION, INTENT (OUT) :: L_DELTAU_SLANT &
          ( MAX_ATMOSWFS, MAXLAYERS, MAXLAYERS, MAXBEAMS )
      LOGICAL, INTENT (OUT) ::          DO_SCATMAT_VARIATION &
          ( MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_TRUNC_FACTOR &
          ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: L_OMEGA_GREEK &
         ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LC_AVERAGE_SECANT &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LC_INITIAL_TRANS &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_T_DELT_DISORDS &
          ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_T_DISORDS_UTDN &
          ( MAXSTREAMS, MAX_USER_LEVELS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_T_DISORDS_UTUP &
          ( MAXSTREAMS, MAX_USER_LEVELS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LC_T_DELT_MUBAR &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LC_T_UTDN_MUBAR &
          ( MAX_USER_LEVELS, MAXBEAMS, MAX_ATMOSWFS)
      DOUBLE PRECISION, INTENT (OUT) :: L_T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_T_UTDN_USERM &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_T_UTUP_USERM &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS )

!  Miscellaneous setup operations for linearized quantities

!  Deltam scaling of variational quantities

      CALL VLIDORT_LA_DELTAMSCALE ( &
        DO_DELTAM_SCALING, NSTOKES, &
        NLAYERS, OMEGA_TOTAL_INPUT, &
        DELTAU_VERT_INPUT, GREEKMAT_TOTAL_INPUT, &
        NMOMENTS, NSTOKES_SQ, &
        NBEAMS, &
        DO_ATMOS_LINEARIZATION, LAYER_VARY_FLAG, &
        LAYER_VARY_NUMBER, L_OMEGA_TOTAL_INPUT, &
        L_DELTAU_VERT_INPUT, L_GREEKMAT_TOTAL_INPUT, &
        DELTAU_SLANT, TRUNC_FACTOR, FAC1, &
        L_OMEGA_TOTAL, L_DELTAU_VERT, &
        L_GREEKMAT_TOTAL, L_DELTAU_SLANT, &
        DO_SCATMAT_VARIATION, L_TRUNC_FACTOR )

!  Initialise single scatter albedo variational quantities

      CALL VLIDORT_LA_SSALBINIT ( &
        NSTOKES, NLAYERS, &
        NMOMENTS, MUELLER_INDEX, &
        DO_ATMOS_LINEARIZATION, LAYER_VARY_FLAG, &
        LAYER_VARY_NUMBER, OMEGA_GREEK, &
        L_OMEGA_TOTAL, L_GREEKMAT_TOTAL, &
        L_OMEGA_GREEK )

!  Linearization of transmittances

      CALL VLIDORT_LA_PREPTRANS ( &
        DO_SOLUTION_SAVING, DO_USER_STREAMS, &
        NLAYERS, NSTREAMS, QUAD_STREAMS, &
        N_USER_STREAMS, USER_SECANTS, &
        N_PARTLAYERS, PARTLAYERS_LAYERIDX, &
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
        DELTAU_VERT, PARTAU_VERT, &
        T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN, &
        T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM, &
        L_DELTAU_VERT, &
        L_T_DELT_DISORDS, L_T_DISORDS_UTDN, L_T_DISORDS_UTUP, &
        L_T_DELT_USERM, L_T_UTDN_USERM, L_T_UTUP_USERM )

!  Linearization of pseudo-spherical setup

      CALL VLIDORT_LC_PREPTRANS ( &
        DO_PLANE_PARALLEL, NLAYERS, NBEAMS, N_PARTLAYERS,      & ! Input
        PARTLAYERS_LAYERIDX, LAYER_VARY_NUMBER,                & ! Input
        DELTAU_VERT, PARTAU_VERT, DELTAU_SLANT, L_DELTAU_VERT, & ! Input
        AVERAGE_SECANT, LAYER_PIS_CUTOFF,                      & ! Input
        T_DELT_MUBAR,   T_UTDN_MUBAR,                          & ! Input
        LC_T_DELT_MUBAR,  LC_T_UTDN_MUBAR,                     & ! Output
        LC_INITIAL_TRANS, LC_AVERAGE_SECANT )                    ! Output

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_LAC_MISCSETUPS

!

      SUBROUTINE VLIDORT_LC_PREPTRANS ( &
        DO_PLANE_PARALLEL, NLAYERS, NBEAMS, N_PARTLAYERS,        & ! Input
        PARTLAYERS_LAYERIDX, LAYER_VARY_NUMBER,                  & ! Input
        DELTAU_VERT, PARTAU_VERT, DELTAU_SLANT, L_DELTAU_VERT,   & ! Input
        AVERAGE_SECANT, LAYER_PIS_CUTOFF,                        & ! Input
        T_DELT_MUBAR,   T_UTDN_MUBAR,                            & ! Input
        LC_T_DELT_MUBAR,  LC_T_UTDN_MUBAR,                       & ! Output
        LC_INITIAL_TRANS, LC_AVERAGE_SECANT )                      ! Output

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NBEAMS
      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      INTEGER, INTENT (IN) ::          LAYER_VARY_NUMBER ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PARTAU_VERT ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_SLANT &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      INTEGER, INTENT (IN) ::          LAYER_PIS_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

      DOUBLE PRECISION, INTENT (OUT) :: LC_T_DELT_MUBAR &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LC_T_UTDN_MUBAR &
          ( MAX_USER_LEVELS, MAXBEAMS, MAX_ATMOSWFS)
      DOUBLE PRECISION, INTENT (OUT) :: LC_INITIAL_TRANS &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LC_AVERAGE_SECANT &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      INTEGER ::          N, Q, UT, K, IB
      DOUBLE PRECISION :: WDEL, VAR, RHO, FAC, DELT, LAMDA, SUM

!  linearization of Initial transmittances
!  =======================================

!   Bug fixed, 12 August 2005 for linearization of INITIAL_TRANS
!         Use Logarithmic derivative !!!!
!         Reason: avoids exceptions if INITIAL_TRANS underflows

      DO IB = 1, NBEAMS
        N = 1
        DO Q = 1, LAYER_VARY_NUMBER(N)
          LC_INITIAL_TRANS(N,IB,Q) = ZERO
        ENDDO
        DO N = 2, NLAYERS
          IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              SUM = ZERO  
              DO K = 1, N-1
                SUM = SUM + L_DELTAU_VERT(Q,K)*DELTAU_SLANT(N-1,K,IB)
              ENDDO
              LC_INITIAL_TRANS(N,IB,Q) = - SUM
            ENDDO
          ELSE
            DO Q = 1, LAYER_VARY_NUMBER(N)
              LC_INITIAL_TRANS(N,IB,Q) = ZERO
            ENDDO
          ENDIF
        ENDDO
      ENDDO

!  linearization of average secants for pseudo-spherical case
!  ==========================================================

!   (average secant = 1/mu-0 = constant for plane parallel)

      IF ( .NOT. DO_PLANE_PARALLEL ) THEN
        DO IB = 1, NBEAMS
          N = 1
          DO Q = 1, LAYER_VARY_NUMBER(N)
            LC_AVERAGE_SECANT(N,IB,Q) = ZERO
          ENDDO
          DO N = 2, NLAYERS
            IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
              DELT  = DELTAU_VERT(N)
              LAMDA = AVERAGE_SECANT(N,IB)
              FAC   = ( DELTAU_SLANT(N,N,IB) / DELT ) - LAMDA
              DO Q = 1, LAYER_VARY_NUMBER(N)
                LC_AVERAGE_SECANT(N,IB,Q) = L_DELTAU_VERT(Q,N) * FAC
              ENDDO
              DO K = 1, N-1
                FAC = ( DELTAU_SLANT(N,K,IB) - DELTAU_SLANT(N-1,K,IB) ) / DELT
                DO Q = 1, LAYER_VARY_NUMBER(K)
                  LC_AVERAGE_SECANT(N,IB,Q) = &
                     LC_AVERAGE_SECANT(N,IB,Q) + L_DELTAU_VERT(Q,K)*FAC
                ENDDO
              ENDDO
            ELSE
              DO Q = 1, LAYER_VARY_NUMBER(N)
                LC_AVERAGE_SECANT(N,IB,Q) = ZERO
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDIF

!  Linearization of Whole layer Transmittance factors
!  ==================================================

      DO IB = 1, NBEAMS
        DO N = 1, NLAYERS

         WDEL  = T_DELT_MUBAR(N,IB)
         VAR   = - DELTAU_VERT(N) * WDEL
         LAMDA = AVERAGE_SECANT(N,IB)
         FAC   = VAR * AVERAGE_SECANT(N,IB)

!  Pseudo-spherical

         IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          IF ( N .EQ. 1 ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              LC_T_DELT_MUBAR(N,IB,Q) = FAC * L_DELTAU_VERT(Q,N)
            ENDDO
          ELSE
            IF  ( N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                RHO = LC_AVERAGE_SECANT(N,IB,Q)
                LC_T_DELT_MUBAR(N,IB,Q) = L_DELTAU_VERT(Q,N) * FAC + VAR * RHO
              ENDDO
            ENDIF
          ENDIF

!  Plane-parallel

         ELSE IF ( DO_PLANE_PARALLEL ) THEN

          IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              LC_T_DELT_MUBAR(N,IB,Q) = FAC * L_DELTAU_VERT(Q,N)
            ENDDO
          ENDIF

         ENDIF

!  end layer and beam loops

        ENDDO
      ENDDO

!  Partial layer transmittance factors (for off-grid optical depths)
!  =================================================================

      DO IB = 1, NBEAMS

!  zero it

        DO UT = 1, N_PARTLAYERS
          N = PARTLAYERS_LAYERIDX(UT)
          DO Q = 1, LAYER_VARY_NUMBER(N)
            LC_T_UTDN_MUBAR(UT,IB,Q) = ZERO
          ENDDO
        ENDDO

        DO UT = 1, N_PARTLAYERS
         N = PARTLAYERS_LAYERIDX(UT)
         VAR = - PARTAU_VERT(UT) * T_UTDN_MUBAR(UT,IB)
         FAC = VAR * AVERAGE_SECANT(N,IB)

!  Pseudo-spherical

         IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          IF ( N .EQ. 1 ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              LC_T_UTDN_MUBAR(UT,IB,Q) = FAC *  L_DELTAU_VERT(Q,N)
            ENDDO
          ELSE
            IF ( N .LE. LAYER_PIS_CUTOFF(IB) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                RHO = LC_AVERAGE_SECANT(N,IB,Q)
                LC_T_UTDN_MUBAR(UT,IB,Q) = L_DELTAU_VERT(Q,N)* FAC + VAR * RHO
              ENDDO
            ENDIF
          ENDIF

!  Plane-parallel

         ELSE IF ( DO_PLANE_PARALLEL ) THEN

          IF ( N .LE. LAYER_PIS_CUTOFF(IB) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              LC_T_UTDN_MUBAR(UT,IB,Q) = FAC * L_DELTAU_VERT(Q,N)
            ENDDO
          ENDIF

         ENDIF

!  End optical depth and beam loops

        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_LC_PREPTRANS

!

      SUBROUTINE LC_EMULT_MASTER ( &
        DO_UPWELLING, DO_DNWELLING, TAYLOR_ORDER, & ! New for 2p7
        NLAYERS, N_USER_LEVELS, &
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
        PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, &
        STERM_LAYERMASK_DN, &
        EMULT_UP, EMULT_DN, &
        N_TOTALCOLUMN_WFS, DO_PLANE_PARALLEL, &
        LAYER_PIS_CUTOFF, NBEAMS, &
        N_USER_STREAMS, &
        T_DELT_MUBAR, T_DELT_USERM, &
        ITRANS_USERM, SIGMA_P, &
        LC_AVERAGE_SECANT, LC_INITIAL_TRANS, &
        LC_T_DELT_MUBAR, L_T_DELT_USERM, &
        USER_SECANTS, &
        DELTAU_VERT, EMULT_HOPRULE, SIGMA_M, &
        L_DELTAU_VERT, T_UTUP_USERM, T_UTDN_USERM, & ! @@@ 5/10/13
        UT_EMULT_UP, &
        LC_T_UTDN_MUBAR, L_T_UTUP_USERM, &
        PARTAU_VERT, UT_EMULT_DN, &
        L_T_UTDN_USERM, &
        LC_EMULT_UP, LC_EMULT_DN, &
        LC_UT_EMULT_UP, LC_UT_EMULT_DN )

!  Linearized multipliers for the Beam source terms

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_DNWELLING
      INTEGER, intent(in)  ::          TAYLOR_ORDER ! 2p7:  Order of Taylor series (including terms up to EPS^n)
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
      INTEGER, INTENT (IN) ::          N_TOTALCOLUMN_WFS
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
      DOUBLE PRECISION, INTENT (IN) :: LC_AVERAGE_SECANT &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_INITIAL_TRANS &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_DELT_MUBAR &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: USER_SECANTS ( MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::          EMULT_HOPRULE &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: SIGMA_M &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_USERM &          ! @@@ 5/10/13
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: UT_EMULT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_UTDN_MUBAR &
          ( MAX_USER_LEVELS, MAXBEAMS, MAX_ATMOSWFS)
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTUP_USERM &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PARTAU_VERT ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_EMULT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_USERM &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (OUT) :: LC_EMULT_UP ( MAX_USER_STREAMS, &
            MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LC_EMULT_DN ( MAX_USER_STREAMS, &
            MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LC_UT_EMULT_UP ( MAX_USER_STREAMS, &
            MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LC_UT_EMULT_DN ( MAX_USER_STREAMS, &
            MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Local variables
!  ---------------

      INTEGER :: N, UT, UTA

!mick fix 7/23/2014 - initialized for packing
      LC_EMULT_UP = ZERO
      LC_EMULT_DN = ZERO
      LC_UT_EMULT_UP = ZERO
      LC_UT_EMULT_DN = ZERO

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
          CALL LC_WHOLELAYER_EMULT_UP ( &
            N, N_TOTALCOLUMN_WFS, &
            DO_PLANE_PARALLEL, &
            LAYER_PIS_CUTOFF, NBEAMS, &
            N_USER_STREAMS, &
            T_DELT_MUBAR, T_DELT_USERM, &
            ITRANS_USERM, &
            EMULT_UP, SIGMA_P, &
            LC_AVERAGE_SECANT, LC_INITIAL_TRANS, &
            LC_T_DELT_MUBAR, L_T_DELT_USERM, &
            LC_EMULT_UP )
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
           CALL LC_PARTLAYER_EMULT_UP ( &
             N, UT, N_TOTALCOLUMN_WFS, &
             DO_PLANE_PARALLEL, &
             LAYER_PIS_CUTOFF, NBEAMS, &
             N_USER_STREAMS, T_DELT_MUBAR, &
             T_UTUP_USERM, ITRANS_USERM, &
             UT_EMULT_UP, SIGMA_P, &
             LC_AVERAGE_SECANT, LC_INITIAL_TRANS, &
             LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR, &
             L_T_UTUP_USERM, &
             LC_UT_EMULT_UP )
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

       DO N = 1, NLAYERS
        IF ( STERM_LAYERMASK_DN(N) ) THEN
          CALL LC_WHOLELAYER_EMULT_DN ( &
            N, N_TOTALCOLUMN_WFS, TAYLOR_ORDER, & ! New for 2p7
            DO_PLANE_PARALLEL, &
            LAYER_PIS_CUTOFF, NBEAMS, &
            N_USER_STREAMS, USER_SECANTS, &
            DELTAU_VERT, T_DELT_USERM, ITRANS_USERM, &
            EMULT_DN, EMULT_HOPRULE, SIGMA_M, &
            L_DELTAU_VERT, LC_AVERAGE_SECANT, &
            LC_INITIAL_TRANS, LC_T_DELT_MUBAR, &
            L_T_DELT_USERM, &
            LC_EMULT_DN )
        ENDIF
       ENDDO

!  Partial layer downwelling
!  -------------------------

!  Start loop over all partial output UT occuring in layers N

       DO UTA = 1, N_USER_LEVELS
        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
         UT = PARTLAYERS_OUTINDEX(UTA)
         N  = PARTLAYERS_LAYERIDX(UT)
         IF ( STERM_LAYERMASK_DN(N) ) THEN
           CALL LC_PARTLAYER_EMULT_DN ( &
             N, UT, N_TOTALCOLUMN_WFS, TAYLOR_ORDER, & ! New for 2p7
             DO_PLANE_PARALLEL, &
             LAYER_PIS_CUTOFF, NBEAMS, &
             N_USER_STREAMS, USER_SECANTS, &
             PARTAU_VERT, T_UTDN_USERM,  ITRANS_USERM, &
             UT_EMULT_DN, EMULT_HOPRULE, SIGMA_M, &
             L_DELTAU_VERT, LC_AVERAGE_SECANT, &
             LC_INITIAL_TRANS, LC_T_UTDN_MUBAR, &
             L_T_UTDN_USERM, &
             LC_UT_EMULT_DN )
         ENDIF
        ENDIF
       ENDDO

!  end downwelling

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LC_EMULT_MASTER

!

      SUBROUTINE LC_WHOLELAYER_EMULT_UP ( &
        N, K_PARAMETERS, &
        DO_PLANE_PARALLEL, &
        LAYER_PIS_CUTOFF, NBEAMS, &
        N_USER_STREAMS, &
        T_DELT_MUBAR, T_DELT_USERM, &
        ITRANS_USERM, &
        EMULT_UP, SIGMA_P, &
        LC_AVERAGE_SECANT, LC_INITIAL_TRANS, &
        LC_T_DELT_MUBAR, L_T_DELT_USERM, &
        LC_EMULT_UP )

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          N
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
      DOUBLE PRECISION, INTENT (IN) :: LC_AVERAGE_SECANT &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_INITIAL_TRANS &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_DELT_MUBAR &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (INOUT) :: LC_EMULT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

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
             LC_EMULT_UP(UM,N,IB,Q) = ZERO
           ENDDO
         ENDDO

       ELSE

!  transmittance factor

        WDEL = T_DELT_MUBAR(N,IB)

!  For the pseudo-spherical case
!  -----------------------------

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          DO UM = 1, N_USER_STREAMS
            UDEL = T_DELT_USERM(N,UM)
            SU = - ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
            DO Q = 1, K_PARAMETERS
              V1 = -LC_AVERAGE_SECANT(N,IB,Q) / SIGMA_P(N,UM,IB)
              V1 = V1 + LC_INITIAL_TRANS (N,IB,Q)
              V2 = WDEL * L_T_DELT_USERM(N,UM,Q) + &
                   UDEL * LC_T_DELT_MUBAR(N,IB,Q)
              LC_EMULT_UP(UM,N,IB,Q) = EMULT_UP(UM,N,IB) * V1 + SU * V2
            ENDDO
          ENDDO

!  For the plane-parallel case
!  ---------------------------

        ELSE

          DO UM = 1, N_USER_STREAMS
            UDEL = T_DELT_USERM(N,UM)
            SU = - ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
            DO Q = 1, K_PARAMETERS
              V1 = LC_INITIAL_TRANS (N,IB,Q)
              V2 = WDEL * L_T_DELT_USERM(N,UM,Q) + &
                   UDEL * LC_T_DELT_MUBAR(N,IB,Q)
              LC_EMULT_UP(UM,N,IB,Q) = EMULT_UP(UM,N,IB)*V1 + SU * V2
            ENDDO
          ENDDO

!  End clause pseudo-spherical versus plane-parallel

        ENDIF

!  continuation point for next beam

       ENDIF

!  End beam loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LC_WHOLELAYER_EMULT_UP

!

      SUBROUTINE LC_WHOLELAYER_EMULT_DN ( &
        N, K_PARAMETERS, TAYLOR_ORDER, & ! New for 2p7
        DO_PLANE_PARALLEL, &
        LAYER_PIS_CUTOFF, NBEAMS, &
        N_USER_STREAMS, USER_SECANTS, &
        DELTAU_VERT, T_DELT_USERM,  ITRANS_USERM, &
        EMULT_DN, EMULT_HOPRULE, SIGMA_M, &
        L_DELTAU_VERT, LC_AVERAGE_SECANT, &
        LC_INITIAL_TRANS, LC_T_DELT_MUBAR, &
        L_T_DELT_USERM, &
        LC_EMULT_DN )

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          N
      INTEGER, INTENT (IN) ::          K_PARAMETERS
      INTEGER, intent(in)  ::          TAYLOR_ORDER ! 2p7:  Order of Taylor series (including terms up to EPS^n)
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL
      INTEGER, INTENT (IN) ::          LAYER_PIS_CUTOFF ( MAXBEAMS )
      INTEGER, INTENT (IN) ::          NBEAMS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      DOUBLE PRECISION, INTENT (IN) :: USER_SECANTS  ( MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT ( MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )

      DOUBLE PRECISION, INTENT (IN) :: ITRANS_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: EMULT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      LOGICAL, INTENT (IN) ::          EMULT_HOPRULE &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: SIGMA_M &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LC_AVERAGE_SECANT &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_INITIAL_TRANS &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_DELT_MUBAR &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (INOUT) :: LC_EMULT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      DOUBLE PRECISION :: SD, V1, V2
      DOUBLE PRECISION :: EPS, SM, DELTA, MULT, TMEW, L_DELTA, UDEL, L_LAM, L_MULT
      INTEGER          :: UM, Q, IB

!  Start Beam loop
!  ===============

      DO IB = 1, NBEAMS

!  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN

         DO UM = 1, N_USER_STREAMS
           DO Q = 1, K_PARAMETERS
             LC_EMULT_DN(UM,N,IB,Q) = ZERO
           ENDDO
         ENDDO

       ELSE

!  NOTE - use of L'Hopital's Rule is present in this module

!  For the pseudo-spherical case
!  -----------------------------

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN

!  Note the use of L'Hopital's Rule flag.

         DO UM = 1, N_USER_STREAMS
           SM   = USER_SECANTS(UM) ; MULT = EMULT_DN(UM,N,IB) ; TMEW  = ITRANS_USERM(N,UM,IB) 
           IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
             UDEL = T_DELT_USERM(N,UM) ; EPS  = - SIGMA_M(N,UM,IB) ; DELTA = DELTAU_VERT(N)
             DO Q = 1, K_PARAMETERS
               L_LAM   = LC_AVERAGE_SECANT(N,IB,Q)
               L_DELTA = L_DELTAU_VERT(Q,N) * DELTA ! Input is single normalized
               CALL TAYLOR_SERIES_L_1 &
                     ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, L_LAM, ZERO, UDEL, SM, L_mult )
               LC_EMULT_DN(UM,N,IB,Q) = LC_INITIAL_TRANS(N,IB,Q) * MULT + TMEW * L_MULT
             ENDDO
           ELSE
             SD = TMEW / SIGMA_M(N,UM,IB)
             DO Q = 1, K_PARAMETERS
               V1 = - LC_AVERAGE_SECANT(N,IB,Q) / SIGMA_M(N,UM,IB)
               V1 = V1 + LC_INITIAL_TRANS (N,IB,Q)
               V2 = L_T_DELT_USERM(N,UM,Q) - LC_T_DELT_MUBAR(N,IB,Q)
               LC_EMULT_DN(UM,N,IB,Q) = EMULT_DN(UM,N,IB)*V1 + SD*V2
             ENDDO
           ENDIF
         ENDDO

!  For the plane-parallel case
!  ---------------------------

        ELSE

         DO UM = 1, N_USER_STREAMS
           SM   = USER_SECANTS(UM)   ; MULT = EMULT_DN(UM,N,IB)  ; TMEW  = ITRANS_USERM(N,UM,IB) 
           IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
             UDEL = T_DELT_USERM(N,UM) ; EPS  = - SIGMA_M(N,UM,IB) ; DELTA = DELTAU_VERT(N)
             DO Q = 1, K_PARAMETERS
               L_DELTA = L_DELTAU_VERT(Q,N) * DELTA ! Input is single normalized
               CALL TAYLOR_SERIES_L_1 &
                     ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, ZERO, ZERO, UDEL, SM, L_mult )
               LC_EMULT_DN(UM,N,IB,Q) = LC_INITIAL_TRANS(N,IB,Q) * MULT + TMEW * L_MULT
             ENDDO
           ELSE
             SD = TMEW / SIGMA_M(N,UM,IB)
             DO Q = 1, K_PARAMETERS
               V1 = LC_INITIAL_TRANS (N,IB,Q)
               V2 = L_T_DELT_USERM(N,UM,Q) - LC_T_DELT_MUBAR(N,IB,Q)
               LC_EMULT_DN(UM,N,IB,Q) = EMULT_DN(UM,N,IB)*V1 + SD*V2
             ENDDO
           ENDIF
         ENDDO

!  End clause pseudo-spherical versus plaen-parallel

        ENDIF

!  continuation point for next beam

       ENDIF

!  End beam loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LC_WHOLELAYER_EMULT_DN

!

      SUBROUTINE LC_PARTLAYER_EMULT_UP ( &
        N, UT, K_PARAMETERS, &
        DO_PLANE_PARALLEL, &
        LAYER_PIS_CUTOFF, NBEAMS, &
        N_USER_STREAMS, T_DELT_MUBAR, &
        T_UTUP_USERM, ITRANS_USERM, &
        UT_EMULT_UP, SIGMA_P, &
        LC_AVERAGE_SECANT, LC_INITIAL_TRANS, &
        LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR, &
        L_T_UTUP_USERM, &
        LC_UT_EMULT_UP )

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          N
      INTEGER, INTENT (IN) ::          UT
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
      DOUBLE PRECISION, INTENT (IN) :: LC_AVERAGE_SECANT &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_INITIAL_TRANS &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_DELT_MUBAR &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_UTDN_MUBAR &
          ( MAX_USER_LEVELS, MAXBEAMS, MAX_ATMOSWFS)
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTUP_USERM &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (INOUT) :: LC_UT_EMULT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

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
             LC_UT_EMULT_UP(UM,UT,IB,Q) = ZERO
           ENDDO
         ENDDO

       ELSE

!  transmittance factor

        WDEL = T_DELT_MUBAR(N,IB)

!  For the pseudo-spherical case
!  -----------------------------

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          DO UM = 1, N_USER_STREAMS
            UX_UP = T_UTUP_USERM(UT,UM)
            SU = ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
            DO Q = 1, K_PARAMETERS
              V1 = -LC_AVERAGE_SECANT(N,IB,Q) / SIGMA_P(N,UM,IB)
              V1 = V1 + LC_INITIAL_TRANS (N,IB,Q)
              V2 =           LC_T_UTDN_MUBAR(UT,IB,Q) - &
                     UX_UP * LC_T_DELT_MUBAR(N, IB,Q) - &
                     WDEL  * L_T_UTUP_USERM(UT,UM,Q)
              LC_UT_EMULT_UP(UM,UT,IB,Q) =       SU * V2 + &
                                 UT_EMULT_UP(UM,UT,IB) * V1
            ENDDO
          ENDDO

!  For the plane-parallel case
!  ---------------------------

        ELSE

          DO UM = 1, N_USER_STREAMS
            UX_UP = T_UTUP_USERM(UT,UM)
            SU = ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
            DO Q = 1, K_PARAMETERS
              V1 = LC_INITIAL_TRANS (N,IB,Q)
              V2 =           LC_T_UTDN_MUBAR(UT,IB,Q) - &
                     UX_UP * LC_T_DELT_MUBAR( N,IB,Q) - &
                     WDEL  * L_T_UTUP_USERM(UT,UM,Q)
              LC_UT_EMULT_UP(UM,UT,IB,Q) =       SU * V2 + &
                                 UT_EMULT_UP(UM,UT,IB) * V1
            ENDDO
          ENDDO

!  End clause pseudo-spherical versus plaen-parallel

        ENDIF

!  continuation point for next beam

       ENDIF

!  End beam loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LC_PARTLAYER_EMULT_UP

!

      SUBROUTINE LC_PARTLAYER_EMULT_DN ( &
        N, UT, K_PARAMETERS, TAYLOR_ORDER, & ! New for 2p7
        DO_PLANE_PARALLEL, &
        LAYER_PIS_CUTOFF, NBEAMS, &
        N_USER_STREAMS, USER_SECANTS, &
        PARTAU_VERT, T_UTDN_USERM, ITRANS_USERM, & ! @@@ 5/10/13
        UT_EMULT_DN, EMULT_HOPRULE, SIGMA_M, &
        L_DELTAU_VERT, LC_AVERAGE_SECANT, &
        LC_INITIAL_TRANS, LC_T_UTDN_MUBAR, &
        L_T_UTDN_USERM, &
        LC_UT_EMULT_DN )

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          N
      INTEGER, INTENT (IN) ::          UT
      INTEGER, INTENT (IN) ::          K_PARAMETERS
      INTEGER, intent(in)  ::          TAYLOR_ORDER ! 2p7:  Order of Taylor series (including terms up to EPS^n)
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL
      INTEGER, INTENT (IN) ::          LAYER_PIS_CUTOFF ( MAXBEAMS )
      INTEGER, INTENT (IN) ::          NBEAMS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      DOUBLE PRECISION, INTENT (IN) :: USER_SECANTS  ( MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: PARTAU_VERT ( MAX_PARTLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS ) ! @@@ 5/10/13

      DOUBLE PRECISION, INTENT (IN) :: ITRANS_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: UT_EMULT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      LOGICAL, INTENT (IN) ::          EMULT_HOPRULE &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: SIGMA_M &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LC_AVERAGE_SECANT &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_INITIAL_TRANS &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_UTDN_MUBAR &
          ( MAX_USER_LEVELS, MAXBEAMS, MAX_ATMOSWFS)
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_USERM &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (INOUT) :: LC_UT_EMULT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      DOUBLE PRECISION :: SD, V1, V2
      DOUBLE PRECISION :: EPS, SM, MULT, TMEW, DELTA, L_DELTA, UXDN, L_LAM, L_MULT
      INTEGER          :: UM, Q, IB

!  Start Beam loop
!  ===============

      DO IB = 1, NBEAMS

!  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN

         DO UM = 1, N_USER_STREAMS
           DO Q = 1, K_PARAMETERS
             LC_UT_EMULT_DN(UM,UT,IB,Q) = ZERO
           ENDDO
         ENDDO

       ELSE

!  NOTE - use of L'Hopital's Rule is present in this module

!  For the pseudo-spherical case
!  -----------------------------

       IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          DO UM = 1, N_USER_STREAMS
            SM    = USER_SECANTS(UM) ; MULT = UT_EMULT_DN(UM,UT,IB) ; TMEW  = ITRANS_USERM(N,UM,IB) 
            IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
              UXDN = T_UTDN_USERM(UT,UM) ; EPS = - SIGMA_M(N,UM,IB)  ; DELTA = PARTAU_VERT(UT)
              DO Q = 1, K_PARAMETERS
                L_LAM   = LC_AVERAGE_SECANT(N,IB,Q)
                L_DELTA = L_DELTAU_VERT(Q,N) * DELTA ! Input is single normalized
                CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, L_LAM, ZERO, UXDN, SM, L_mult )
                LC_UT_EMULT_DN(UM,UT,IB,Q) = LC_INITIAL_TRANS(N,IB,Q) * MULT + TMEW * L_MULT
              ENDDO
            ELSE
              SD = TMEW / SIGMA_M(N,UM,IB)
              DO Q = 1, K_PARAMETERS
                V1 = - LC_AVERAGE_SECANT(N,IB,Q) / SIGMA_M(N,UM,IB)
                V1 = V1 + LC_INITIAL_TRANS (N,IB,Q)
                V2 = L_T_UTDN_USERM(UT,UM,Q) - LC_T_UTDN_MUBAR(UT,IB,Q)
                LC_UT_EMULT_DN(UM,UT,IB,Q) =       SD * V2 + &
                                 UT_EMULT_DN(UM,UT,IB) * V1
              ENDDO
            ENDIF
          ENDDO

!  For the plane-parallel case
!  ---------------------------

        ELSE

          DO UM = 1, N_USER_STREAMS
            SM    = USER_SECANTS(UM) ; MULT = UT_EMULT_DN(UM,UT,IB) ; TMEW  = ITRANS_USERM(N,UM,IB) 
            IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
              UXDN = T_UTDN_USERM(UT,UM) ; EPS = - SIGMA_M(N,UM,IB)  ; DELTA = PARTAU_VERT(UT)
              DO Q = 1, K_PARAMETERS
                L_DELTA = L_DELTAU_VERT(Q,N) * DELTA ! Input is single normalized
                CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, ZERO, ZERO, UXDN, SM, L_mult )
                LC_UT_EMULT_DN(UM,UT,IB,Q) = LC_INITIAL_TRANS(N,IB,Q) * MULT + TMEW * L_MULT
              ENDDO
            ELSE
              SD = TMEW / SIGMA_M(N,UM,IB)
              DO Q = 1, K_PARAMETERS
                V1 = LC_INITIAL_TRANS (N,IB,Q)
                V2 = L_T_UTDN_USERM(UT,UM,Q) - LC_T_UTDN_MUBAR(UT,IB,Q)
                LC_UT_EMULT_DN(UM,UT,IB,Q) =       SD * V2 + &
                                 UT_EMULT_DN(UM,UT,IB) * V1
              ENDDO
            ENDIF
          ENDDO

!  End clause pseudo-spherical versus plane-parallel

        ENDIF

!  continuation point for next beam

       ENDIF

!  End beam loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LC_PARTLAYER_EMULT_DN

!

      SUBROUTINE LC_EMULT_MASTER_OBSGEO ( &
        DO_UPWELLING, DO_DNWELLING, TAYLOR_ORDER, & ! New for 2p7
        NLAYERS, N_USER_LEVELS, &
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
        PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, &
        STERM_LAYERMASK_DN, &
        EMULT_UP, EMULT_DN, &
        N_TOTALCOLUMN_WFS, DO_PLANE_PARALLEL, &
        LAYER_PIS_CUTOFF, NBEAMS, &
        T_DELT_MUBAR, T_DELT_USERM, &
        ITRANS_USERM, SIGMA_P, &
        LC_AVERAGE_SECANT, LC_INITIAL_TRANS, &
        LC_T_DELT_MUBAR, L_T_DELT_USERM, &
        USER_SECANTS, &
        DELTAU_VERT, EMULT_HOPRULE, SIGMA_M, &
        L_DELTAU_VERT, T_UTUP_USERM, T_UTDN_USERM, &
        UT_EMULT_UP, &
        LC_T_UTDN_MUBAR, L_T_UTUP_USERM, &
        PARTAU_VERT, UT_EMULT_DN, &
        L_T_UTDN_USERM, &
        LC_EMULT_UP, LC_EMULT_DN, &
        LC_UT_EMULT_UP, LC_UT_EMULT_DN )

!  Linearized multipliers for the Beam source terms

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_DNWELLING
      INTEGER, intent(in)  ::          TAYLOR_ORDER ! 2p7:  Order of Taylor series (including terms up to EPS^n)
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
      INTEGER, INTENT (IN) ::          N_TOTALCOLUMN_WFS
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL
      INTEGER, INTENT (IN) ::          LAYER_PIS_CUTOFF ( MAXBEAMS )
      INTEGER, INTENT (IN) ::          NBEAMS
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: ITRANS_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: SIGMA_P &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: LC_AVERAGE_SECANT &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_INITIAL_TRANS &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_DELT_MUBAR &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: USER_SECANTS ( MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::          EMULT_HOPRULE &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: SIGMA_M &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: UT_EMULT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_UTDN_MUBAR &
          ( MAX_USER_LEVELS, MAXBEAMS, MAX_ATMOSWFS)
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTUP_USERM &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PARTAU_VERT ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_EMULT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_USERM &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (OUT) :: LC_EMULT_UP ( MAX_USER_STREAMS, &
            MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LC_EMULT_DN ( MAX_USER_STREAMS, &
            MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LC_UT_EMULT_UP ( MAX_USER_STREAMS, &
            MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LC_UT_EMULT_DN ( MAX_USER_STREAMS, &
            MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Local variables
!  ---------------

      INTEGER :: N, UT, UTA

!mick fix 7/23/2014 - initialized for packing
      LC_EMULT_UP = ZERO
      LC_EMULT_DN = ZERO
      LC_UT_EMULT_UP = ZERO
      LC_UT_EMULT_DN = ZERO

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
          CALL LC_WHOLELAYER_EMULT_OG_UP ( &
            N, N_TOTALCOLUMN_WFS, &
            DO_PLANE_PARALLEL, &
            LAYER_PIS_CUTOFF, NBEAMS, &
            T_DELT_MUBAR, T_DELT_USERM, &
            ITRANS_USERM, &
            EMULT_UP, SIGMA_P, &
            LC_AVERAGE_SECANT, LC_INITIAL_TRANS, &
            LC_T_DELT_MUBAR, L_T_DELT_USERM, &
            LC_EMULT_UP )
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
           CALL LC_PARTLAYER_EMULT_OG_UP ( &
             N, UT, N_TOTALCOLUMN_WFS, &
             DO_PLANE_PARALLEL, &
             LAYER_PIS_CUTOFF, NBEAMS, &
             T_DELT_MUBAR, &
             T_UTUP_USERM, ITRANS_USERM, &
             UT_EMULT_UP, SIGMA_P, &
             LC_AVERAGE_SECANT, LC_INITIAL_TRANS, &
             LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR, &
             L_T_UTUP_USERM, &
             LC_UT_EMULT_UP )
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
          CALL LC_WHOLELAYER_EMULT_OG_DN ( &
            N, N_TOTALCOLUMN_WFS, TAYLOR_ORDER, & ! New for 2p7
            DO_PLANE_PARALLEL, &
            LAYER_PIS_CUTOFF, NBEAMS, &
            USER_SECANTS, &
            DELTAU_VERT, T_DELT_USERM, ITRANS_USERM, &
            EMULT_DN, EMULT_HOPRULE, SIGMA_M, &
            L_DELTAU_VERT, LC_AVERAGE_SECANT, &
            LC_INITIAL_TRANS, LC_T_DELT_MUBAR, &
            L_T_DELT_USERM, &
            LC_EMULT_DN )
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
           CALL LC_PARTLAYER_EMULT_OG_DN ( &
             N, UT, N_TOTALCOLUMN_WFS, TAYLOR_ORDER, & ! New for 2p7
             DO_PLANE_PARALLEL, &
             LAYER_PIS_CUTOFF, NBEAMS, &
             USER_SECANTS, &
             PARTAU_VERT, T_UTDN_USERM, ITRANS_USERM, &
             UT_EMULT_DN, EMULT_HOPRULE, SIGMA_M, &
             L_DELTAU_VERT, LC_AVERAGE_SECANT, &
             LC_INITIAL_TRANS, LC_T_UTDN_MUBAR, &
             L_T_UTDN_USERM, &
             LC_UT_EMULT_DN )
         ENDIF
        ENDIF
       ENDDO

!  end downwelling

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LC_EMULT_MASTER_OBSGEO

!

      SUBROUTINE LC_WHOLELAYER_EMULT_OG_UP ( &
        N, K_PARAMETERS, &
        DO_PLANE_PARALLEL, &
        LAYER_PIS_CUTOFF, NBEAMS, &
        T_DELT_MUBAR, T_DELT_USERM, &
        ITRANS_USERM, &
        EMULT_UP, SIGMA_P, &
        LC_AVERAGE_SECANT, LC_INITIAL_TRANS, &
        LC_T_DELT_MUBAR, L_T_DELT_USERM, &
        LC_EMULT_UP )

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          N
      INTEGER, INTENT (IN) ::          K_PARAMETERS
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL
      INTEGER, INTENT (IN) ::          LAYER_PIS_CUTOFF ( MAXBEAMS )
      INTEGER, INTENT (IN) ::          NBEAMS
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: ITRANS_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: EMULT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: SIGMA_P &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: LC_AVERAGE_SECANT &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_INITIAL_TRANS &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_DELT_MUBAR &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (INOUT) :: LC_EMULT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      DOUBLE PRECISION :: SU, V1, V2, WDEL, UDEL
      INTEGER          :: Q, IB, LUM

!  Start Beam loop
!  ===============

      LUM = 1

      DO IB = 1, NBEAMS

!  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN

         DO Q = 1, K_PARAMETERS
           LC_EMULT_UP(LUM,N,IB,Q) = ZERO
         ENDDO

       ELSE

!  transmittance factor

        WDEL = T_DELT_MUBAR(N,IB)

!  For the pseudo-spherical case
!  -----------------------------

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          UDEL = T_DELT_USERM(N,IB)
          SU = - ITRANS_USERM(N,LUM,IB) / SIGMA_P(N,LUM,IB)
          DO Q = 1, K_PARAMETERS
            V1 = -LC_AVERAGE_SECANT(N,IB,Q) / SIGMA_P(N,LUM,IB)
            V1 = V1 + LC_INITIAL_TRANS(N,IB,Q)
            V2 = WDEL * L_T_DELT_USERM(N,IB,Q) + UDEL * LC_T_DELT_MUBAR(N,IB,Q)
            LC_EMULT_UP(LUM,N,IB,Q) = EMULT_UP(LUM,N,IB) * V1 + SU * V2
          ENDDO

!  For the plane-parallel case
!  ---------------------------

        ELSE

          UDEL = T_DELT_USERM(N,IB)
          SU = - ITRANS_USERM(N,LUM,IB) / SIGMA_P(N,LUM,IB)
          DO Q = 1, K_PARAMETERS
            V1 = LC_INITIAL_TRANS (N,IB,Q)
            V2 = WDEL * L_T_DELT_USERM(N,IB,Q) + UDEL * LC_T_DELT_MUBAR(N,IB,Q)
            LC_EMULT_UP(LUM,N,IB,Q) = EMULT_UP(LUM,N,IB)*V1 + SU * V2
          ENDDO

!  End clause pseudo-spherical versus plane-parallel

        ENDIF

!  continuation point for next beam

       ENDIF

!  End beam loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LC_WHOLELAYER_EMULT_OG_UP

!

      SUBROUTINE LC_WHOLELAYER_EMULT_OG_DN ( &
        N, K_PARAMETERS, TAYLOR_ORDER, & ! New for 2p7
        DO_PLANE_PARALLEL, &
        LAYER_PIS_CUTOFF, NBEAMS, &
        USER_SECANTS, &
        DELTAU_VERT, T_DELT_USERM, ITRANS_USERM, &
        EMULT_DN, EMULT_HOPRULE, SIGMA_M, &
        L_DELTAU_VERT, LC_AVERAGE_SECANT,  &
        LC_INITIAL_TRANS, LC_T_DELT_MUBAR, &
        L_T_DELT_USERM, &
        LC_EMULT_DN )

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          N
      INTEGER, INTENT (IN) ::          K_PARAMETERS
      INTEGER, intent(in)  ::          TAYLOR_ORDER ! 2p7:  Order of Taylor series (including terms up to EPS^n)
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL
      INTEGER, INTENT (IN) ::          LAYER_PIS_CUTOFF ( MAXBEAMS )
      INTEGER, INTENT (IN) ::          NBEAMS
      DOUBLE PRECISION, INTENT (IN) :: USER_SECANTS  ( MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT ( MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )

      DOUBLE PRECISION, INTENT (IN) :: ITRANS_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: EMULT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      LOGICAL, INTENT (IN) ::          EMULT_HOPRULE &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: SIGMA_M &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LC_AVERAGE_SECANT &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_INITIAL_TRANS &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_DELT_MUBAR &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (INOUT) :: LC_EMULT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      DOUBLE PRECISION :: SD, V1, V2
      DOUBLE PRECISION :: EPS, SM, MULT, TMEW, DELTA, L_DELTA, UDEL, L_LAM, L_MULT
      INTEGER          :: Q, IB, LUM

!  Start Beam loop
!  ===============

      LUM = 1

      DO IB = 1, NBEAMS

!  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN

         DO Q = 1, K_PARAMETERS
           LC_EMULT_DN(LUM,N,IB,Q) = ZERO
         ENDDO

       ELSE

!  NOTE - use of L'Hopital's Rule is present in this module

!  For the pseudo-spherical case
!  -----------------------------

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          SM    = USER_SECANTS(IB) ; MULT = EMULT_DN(LUM,N,IB)   ; TMEW  = ITRANS_USERM(N,LUM,IB) 
          IF ( EMULT_HOPRULE(N,LUM,IB) ) THEN
            UDEL = T_DELT_USERM(N,IB) ; EPS  = - SIGMA_M(N,LUM,IB) ; DELTA = DELTAU_VERT(N)
            DO Q = 1, K_PARAMETERS
              L_LAM   = LC_AVERAGE_SECANT(N,IB,Q)
              L_DELTA = L_DELTAU_VERT(Q,N) * DELTA ! Input is single normalized
              CALL TAYLOR_SERIES_L_1 &
                     ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, L_LAM, ZERO, UDEL, SM, L_mult )
              LC_EMULT_DN(LUM,N,IB,Q) = LC_INITIAL_TRANS(N,IB,Q) * MULT + TMEW * L_MULT
            ENDDO
          ELSE
            SD = TMEW / SIGMA_M(N,LUM,IB)
            DO Q = 1, K_PARAMETERS
              V1 = - LC_AVERAGE_SECANT(N,IB,Q) / SIGMA_M(N,LUM,IB)
              V1 = V1 + LC_INITIAL_TRANS(N,IB,Q)
              V2 = L_T_DELT_USERM(N,IB,Q) - LC_T_DELT_MUBAR(N,IB,Q)
              LC_EMULT_DN(LUM,N,IB,Q) = EMULT_DN(LUM,N,IB)*V1 + SD*V2
            ENDDO
          ENDIF

!  For the plane-parallel case
!  ---------------------------

        ELSE

          SM    = USER_SECANTS(IB) ; MULT = EMULT_DN(LUM,N,IB)   ; TMEW  = ITRANS_USERM(N,LUM,IB) 
          IF ( EMULT_HOPRULE(N,LUM,IB) ) THEN
            UDEL = T_DELT_USERM(N,IB) ; EPS  = - SIGMA_M(N,LUM,IB) ; DELTA = DELTAU_VERT(N)
            DO Q = 1, K_PARAMETERS
              L_DELTA = L_DELTAU_VERT(Q,N) * DELTA ! Input is single normalized
              CALL TAYLOR_SERIES_L_1 &
                     ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, ZERO, ZERO, UDEL, SM, L_mult )
              LC_EMULT_DN(LUM,N,IB,Q) = LC_INITIAL_TRANS(N,IB,Q) * MULT + TMEW * L_MULT
            ENDDO
          ELSE
            SD = TMEW / SIGMA_M(N,LUM,IB)
            DO Q = 1, K_PARAMETERS
              V1 = ZERO
              V1 = LC_INITIAL_TRANS(N,IB,Q)
              V2 = L_T_DELT_USERM(N,IB,Q) - LC_T_DELT_MUBAR(N,IB,Q)
              LC_EMULT_DN(LUM,N,IB,Q) = EMULT_DN(LUM,N,IB)*V1 + SD*V2
            ENDDO
          ENDIF

!  End clause pseudo-spherical versus plane-parallel

        ENDIF

!  continuation point for next beam

       ENDIF

!  End beam loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LC_WHOLELAYER_EMULT_OG_DN

!

      SUBROUTINE LC_PARTLAYER_EMULT_OG_UP ( &
        N, UT, K_PARAMETERS, &
        DO_PLANE_PARALLEL, &
        LAYER_PIS_CUTOFF, NBEAMS, &
        T_DELT_MUBAR, &
        T_UTUP_USERM, ITRANS_USERM, &
        UT_EMULT_UP, SIGMA_P, &
        LC_AVERAGE_SECANT, LC_INITIAL_TRANS, &
        LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR, &
        L_T_UTUP_USERM, &
        LC_UT_EMULT_UP )

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          N
      INTEGER, INTENT (IN) ::          UT
      INTEGER, INTENT (IN) ::          K_PARAMETERS
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL
      INTEGER, INTENT (IN) ::          LAYER_PIS_CUTOFF ( MAXBEAMS )
      INTEGER, INTENT (IN) ::          NBEAMS
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: ITRANS_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: UT_EMULT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: SIGMA_P &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: LC_AVERAGE_SECANT &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_INITIAL_TRANS &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_DELT_MUBAR &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_UTDN_MUBAR &
          ( MAX_USER_LEVELS, MAXBEAMS, MAX_ATMOSWFS)
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTUP_USERM &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (INOUT) :: LC_UT_EMULT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      DOUBLE PRECISION :: SU, V1, V2, WDEL, UX_UP
      INTEGER          :: Q, IB, LUM

!  Start Beam loop
!  ===============

      LUM = 1

      DO IB = 1, NBEAMS

!  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN

         DO Q = 1, K_PARAMETERS
           LC_UT_EMULT_UP(LUM,UT,IB,Q) = ZERO
         ENDDO

       ELSE

!  transmittance factor

        WDEL = T_DELT_MUBAR(N,IB)

!  For the pseudo-spherical case
!  -----------------------------

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          UX_UP = T_UTUP_USERM(UT,IB)
          SU = ITRANS_USERM(N,LUM,IB) / SIGMA_P(N,LUM,IB)
          DO Q = 1, K_PARAMETERS
            V1 = - LC_AVERAGE_SECANT(N,IB,Q) / SIGMA_P(N,LUM,IB)
            V1 = V1 + LC_INITIAL_TRANS(N,IB,Q)
            V2 = LC_T_UTDN_MUBAR(UT,IB,Q) - &
                  UX_UP * LC_T_DELT_MUBAR(N,IB,Q) - WDEL * L_T_UTUP_USERM(UT,IB,Q)
            LC_UT_EMULT_UP(LUM,UT,IB,Q) = SU * V2 + UT_EMULT_UP(LUM,UT,IB) * V1
          ENDDO

!  For the plane-parallel case
!  ---------------------------

        ELSE

          UX_UP = T_UTUP_USERM(UT,IB)
          SU = ITRANS_USERM(N,LUM,IB) / SIGMA_P(N,LUM,IB)
          DO Q = 1, K_PARAMETERS
            V1 = LC_INITIAL_TRANS(N,IB,Q)
            V2 = LC_T_UTDN_MUBAR(UT,IB,Q) - &
                  UX_UP * LC_T_DELT_MUBAR(N,IB,Q) - WDEL * L_T_UTUP_USERM(UT,IB,Q)
            LC_UT_EMULT_UP(LUM,UT,IB,Q) =  SU * V2 + UT_EMULT_UP(LUM,UT,IB) * V1
          ENDDO

        ENDIF

!  continuation point for next beam

       ENDIF

!  End beam loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LC_PARTLAYER_EMULT_OG_UP

!

      SUBROUTINE LC_PARTLAYER_EMULT_OG_DN ( &
        N, UT, K_PARAMETERS, TAYLOR_ORDER, & ! New for 2p7
        DO_PLANE_PARALLEL, &
        LAYER_PIS_CUTOFF, NBEAMS, &
        USER_SECANTS, &
        PARTAU_VERT, T_UTDN_USERM,  ITRANS_USERM, &
        UT_EMULT_DN, EMULT_HOPRULE, SIGMA_M, &
        L_DELTAU_VERT, LC_AVERAGE_SECANT, &
        LC_INITIAL_TRANS, LC_T_UTDN_MUBAR, &
        L_T_UTDN_USERM, &
        LC_UT_EMULT_DN )

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          N
      INTEGER, INTENT (IN) ::          UT
      INTEGER, INTENT (IN) ::          K_PARAMETERS
      INTEGER, intent(in)  ::          TAYLOR_ORDER ! 2p7:  Order of Taylor series (including terms up to EPS^n)
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL
      INTEGER, INTENT (IN) ::          LAYER_PIS_CUTOFF ( MAXBEAMS )
      INTEGER, INTENT (IN) ::          NBEAMS
      DOUBLE PRECISION, INTENT (IN) :: USER_SECANTS  ( MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: PARTAU_VERT ( MAX_PARTLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

      DOUBLE PRECISION, INTENT (IN) :: ITRANS_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: UT_EMULT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      LOGICAL, INTENT (IN) ::          EMULT_HOPRULE &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: SIGMA_M &
          ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LC_AVERAGE_SECANT &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_INITIAL_TRANS &
          ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_UTDN_MUBAR &
          ( MAX_USER_LEVELS, MAXBEAMS, MAX_ATMOSWFS)
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_USERM &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (INOUT) :: LC_UT_EMULT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      DOUBLE PRECISION :: SD, V1, V2
      DOUBLE PRECISION :: EPS, SM, MULT, TMEW, DELTA, L_DELTA, UXDN, L_LAM, L_MULT
      INTEGER          :: Q, IB, LUM

!  Start Beam loop
!  ===============

      LUM = 1

      DO IB = 1, NBEAMS

!  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN

         DO Q = 1, K_PARAMETERS
           LC_UT_EMULT_DN(LUM,UT,IB,Q) = ZERO
         ENDDO

       ELSE

!  NOTE - use of L'Hopital's Rule is present in this module

!  For the pseudo-spherical case
!  -----------------------------

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          SM  = USER_SECANTS(IB) ;  MULT = UT_EMULT_DN(LUM,UT,IB) ; TMEW  = ITRANS_USERM(N,LUM,IB) 
          IF ( EMULT_HOPRULE(N,LUM,IB) ) THEN
            UXDN = T_UTDN_USERM(UT,IB) ; EPS = - SIGMA_M(N,LUM,IB)  ; DELTA = PARTAU_VERT(UT)
            DO Q = 1, K_PARAMETERS
              L_LAM   = LC_AVERAGE_SECANT(N,IB,Q)
              L_DELTA = L_DELTAU_VERT(Q,N) * DELTA ! Input is single normalized
              CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, L_LAM, ZERO, UXDN, SM, L_mult )
              LC_UT_EMULT_DN(LUM,UT,IB,Q) = LC_INITIAL_TRANS(N,IB,Q) * MULT + TMEW * L_MULT
            ENDDO
          ELSE
            SD = TMEW / SIGMA_M(N,LUM,IB)
            DO Q = 1, K_PARAMETERS
              V1 = - LC_AVERAGE_SECANT(N,IB,Q) / SIGMA_M(N,LUM,IB)
              V1 = V1 + LC_INITIAL_TRANS(N,IB,Q)
              V2 = L_T_UTDN_USERM(UT,IB,Q) - LC_T_UTDN_MUBAR(UT,IB,Q)
              LC_UT_EMULT_DN(LUM,UT,IB,Q) = SD * V2 + UT_EMULT_DN(LUM,UT,IB) * V1
            ENDDO
          ENDIF

!  For the plane-parallel case
!  ---------------------------

        ELSE

          SM  = USER_SECANTS(IB) ;  MULT = UT_EMULT_DN(LUM,UT,IB) ; TMEW  = ITRANS_USERM(N,LUM,IB) 
          IF ( EMULT_HOPRULE(N,LUM,IB) ) THEN
            UXDN = T_UTDN_USERM(UT,IB) ; EPS = - SIGMA_M(N,LUM,IB)  ; DELTA = PARTAU_VERT(UT)
            DO Q = 1, K_PARAMETERS
              L_DELTA = L_DELTAU_VERT(Q,N) * DELTA ! Input is single normalized
              CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, ZERO, ZERO, UXDN, SM, L_mult )
              LC_UT_EMULT_DN(LUM,UT,IB,Q) = LC_INITIAL_TRANS(N,IB,Q) * MULT + TMEW * L_MULT
            ENDDO
          ELSE
            DO Q = 1, K_PARAMETERS
              V1 = LC_INITIAL_TRANS(N,IB,Q)
              V2 = L_T_UTDN_USERM(UT,IB,Q) - LC_T_UTDN_MUBAR(UT,IB,Q)
              LC_UT_EMULT_DN(LUM,UT,IB,Q) = SD * V2 + UT_EMULT_DN(LUM,UT,IB) * V1
            ENDDO
          ENDIF

        ENDIF

!  continuation point for next beam

       ENDIF

!  End beam loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LC_PARTLAYER_EMULT_OG_DN

      END MODULE vlidort_lc_miscsetups

