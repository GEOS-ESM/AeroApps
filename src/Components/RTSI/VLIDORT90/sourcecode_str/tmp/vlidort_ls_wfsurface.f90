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

! ##########################################################
! #                                                        #
! # Subroutines in this Module                             #
! #                                                        #
! #     Top level routines--------------                   #
! #            SURFACEWF_MASTER  (master)                  #
! #                                                        #
! #     Linearized BVP Column, surface WFs ---------       #
! #            SURFACEWF_COLSETUP                          #
! #                                                        #
! #     BOA surface source terms ---------                 #
! #            BOA_SURFACEWF                               #
! #                                                        #
! #     Recursion relations ---------                      #
! #            UPUSER_SURFACEWF                            #
! #            DNUSER_SURFACEWF                            #
! #                                                        #
! #     Post-processing at user angles --------            #
! #            LS_WHOLELAYER_STERM_UP                      #
! #            LS_WHOLELAYER_STERM_DN                      #
! #            LS_PARTLAYER_STERM_UP                       #
! #            LS_PARTLAYER_STERM_DN                       #
! #                                                        #
! ##########################################################
! #                                                        #
! #     High level Jacobian routine   ---------            #
! #            VLIDORT_LS_INTEGRATED_OUTPUT (master)       #
! #                                                        #
! #            QUADSURFACEWF_LEVEL_UP                      #
! #            QUADSURFACEWF_LEVEL_DN                      #
! #            QUADSURFACEWF_OFFGRID_UP                    #
! #            QUADSURFACEWF_OFFGRID_DN                    #
! #                                                        #
! ##########################################################

      MODULE vlidort_ls_wfsurface

      PRIVATE
      PUBLIC :: SURFACEWF_MASTER

      CONTAINS

      SUBROUTINE SURFACEWF_MASTER ( &
        DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFEMISS, &
        DO_MSMODE_THERMAL, DO_OBSERVATION_GEOMETRY, &
        DO_INCLUDE_MVOUTPUT, N_SURFACE_WFS, &
        FOURIER_COMPONENT, IBEAM, SURFACE_FACTOR, &
        FLUX_MULTIPLIER, DO_UPWELLING, DO_DNWELLING, &
        DO_QUAD_OUTPUT, NLAYERS, NTOTAL, &
        N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
        K_REAL, K_COMPLEX, &
        BANDMAT2, IPIVOT, SMAT2, SIPIVOT, &
        NSTOKES, NSTREAMS, &
        DO_LAMBERTIAN_SURFACE, SURFBB, LS_BRDF_F, &
        LS_BRDF_F_0, LS_EMISSIVITY, QUAD_STRMWTS, &
        MUELLER_INDEX, T_DELT_EIGEN, SOLA_XPOS, SOLB_XNEG, &
        WLOWER, LCON, MCON, ATMOS_ATTN, &
        R2_HOMP, R2_HOMM, R2_BEAM, DIRECT_BEAM, &
        LAMBERTIAN_ALBEDO, USER_BRDF_F, &
        DO_THERMAL_TRANSONLY, LS_USER_BRDF_F, LS_USER_BRDF_F_0, &
        LS_USER_EMISSIVITY, &
        DO_DBCORRECTION, FLUXVEC, &
        N_USER_STREAMS, DO_USER_STREAMS, LOCAL_UM_START, &
        STOKES_DOWNSURF, N_USER_LEVELS, &
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
        UTAU_LEVEL_MASK_UP, PARTLAYERS_LAYERIDX, &
        T_DELT_USERM, T_UTUP_USERM, &
        UHOM_UPDN, UHOM_UPUP, HMULT_1, HMULT_2, &
        UT_HMULT_UU, UT_HMULT_UD, &
        UTAU_LEVEL_MASK_DN, T_UTDN_USERM, &
        UHOM_DNDN, UHOM_DNUP, UT_HMULT_DU, UT_HMULT_DD, &
        QUAD_WEIGHTS, N_DIRECTIONS, WHICH_DIRECTIONS, &
        T_DELT_DISORDS, T_DISORDS_UTUP, &
        T_UTUP_EIGEN, T_UTDN_EIGEN, &
        SURFACEWF_F, MINT_SURFACEWF, FLUX_SURFACEWF, &
        STATUS, MESSAGE, TRACE )

      USE VLIDORT_PARS
      USE LAPACK_TOOLS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_SURFEMISS
      LOGICAL, INTENT (IN) ::          DO_MSMODE_THERMAL
      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_MVOUTPUT
      INTEGER, INTENT (IN) ::          N_SURFACE_WFS
      DOUBLE PRECISION, INTENT (IN) :: SURFACE_FACTOR
      INTEGER, INTENT (IN) ::          FOURIER_COMPONENT, IBEAM
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER
      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_DNWELLING
      LOGICAL, INTENT (IN) ::          DO_QUAD_OUTPUT
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NTOTAL
      INTEGER, INTENT (IN) ::          N_SUBDIAG
      INTEGER, INTENT (IN) ::          N_SUPDIAG
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS_2
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER, INTENT (IN) ::          IPIVOT ( MAXTOTAL )
      DOUBLE PRECISION, INTENT (IN) :: SMAT2 ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )
      INTEGER, INTENT (IN) ::          SIPIVOT ( MAXSTRMSTKS_2 )
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      DOUBLE PRECISION, INTENT (IN) :: SURFBB
      DOUBLE PRECISION, INTENT (IN) :: LS_BRDF_F &
          ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTOKES_SQ, &
            MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: LS_BRDF_F_0 &
          ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTOKES_SQ, &
            MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: LS_EMISSIVITY &
          ( MAX_SURFACEWFS, MAXSTOKES,MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: WLOWER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: ATMOS_ATTN ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: R2_HOMP &
          ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION, INTENT (IN) :: R2_HOMM &
          ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION, INTENT (IN) :: R2_BEAM ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: DIRECT_BEAM &
          ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: USER_BRDF_F &
          ( 0:MAXMOMENTS,MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      DOUBLE PRECISION, INTENT (IN) :: LS_USER_BRDF_F &
          ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTOKES_SQ, &
            MAX_USER_STREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: LS_USER_BRDF_F_0 &
          ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTOKES_SQ, &
            MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: LS_USER_EMISSIVITY &
           ( MAX_SURFACEWFS, MAXSTOKES, MAX_USER_STREAMS )
      LOGICAL, INTENT (IN) ::          DO_DBCORRECTION
      DOUBLE PRECISION, INTENT (IN) :: FLUXVEC ( MAXSTOKES )
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      LOGICAL, INTENT (IN) ::          DO_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      DOUBLE PRECISION, INTENT (IN) :: STOKES_DOWNSURF ( MAXSTREAMS, MAXSTOKES )
      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPDN &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPUP &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_1 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_2 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_UU &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_UD &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNDN &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNUP &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_DU &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_DD &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_WEIGHTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::          N_DIRECTIONS
      INTEGER, INTENT (IN) ::          WHICH_DIRECTIONS ( MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTUP &
          ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN &
          ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN &
          ( MAXEVALUES, MAX_PARTLAYERS )

!  Linearized surface output

      DOUBLE PRECISION, INTENT (INOUT) ::  SURFACEWF_F &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_VZANGLES, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) ::  MINT_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) ::  FLUX_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  Exception handling

      INTEGER, INTENT (OUT) ::             STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  Local variables
!  ---------------

!  error tracing variables

      INTEGER ::           INFO
      CHARACTER (LEN=3) :: CI

!  Linearized BOA terms

      DOUBLE PRECISION :: LS_BOA_SOURCE &
          ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION :: LS_BOA_THTONLY_SOURCE &
          ( MAX_SURFACEWFS, MAXSTREAMS, MAXSTOKES )

      DOUBLE PRECISION :: COL2_WFALB ( MAXTOTAL, MAX_SURFACEWFS )
      DOUBLE PRECISION :: SCOL2_WFALB ( MAXSTRMSTKS_2, MAX_SURFACEWFS )

      DOUBLE PRECISION :: NCON_ALB &
          ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION :: PCON_ALB &
          ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )

!  Other local variables

      INTEGER ::          N, Q, K, K0, K1, K2, KO1, C0, LUM
      INTEGER ::          IROW, IROW1, IROW_S, IROW1_S

!  Initialise status

      STATUS = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Local user index

      LUM = 1

!  Regular BVP Solution
!  ====================

!   NO TELESCOPING HERE

!  BV solution for perturbed integration constants
!  -----------------------------------------------

!  Compute the main column B' where AX = B'

!       write(*,*)DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFEMISS,
!     I      FOURIER_COMPONENT, N_SURFACE_WFS

      CALL SURFACEWF_COLSETUP ( &
        DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFEMISS, &
        SURFACE_FACTOR, FOURIER_COMPONENT, &
        IBEAM, NSTOKES, NSTREAMS, &
        NLAYERS, DO_LAMBERTIAN_SURFACE, &
        SURFBB, N_SURFACE_WFS, LS_BRDF_F, &
        LS_BRDF_F_0, LS_EMISSIVITY, &
        QUAD_STRMWTS, NTOTAL, &
        NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
        MUELLER_INDEX, T_DELT_EIGEN, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        WLOWER, LCON, MCON, &
        ATMOS_ATTN, R2_HOMP, &
        R2_HOMM, R2_BEAM, DIRECT_BEAM, &
        COL2_WFALB, SCOL2_WFALB )

!  BVP back-substitution: With compression (multilayers)
!  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  LAPACK substitution (DGBTRS) using RHS column vector COL2_WF
!  BV solution for perturbed integration constants
!    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS &
           ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, N_SURFACE_WFS, &
              BANDMAT2, MAXBANDTOTAL, IPIVOT, &
              COL2_WFALB, MAXTOTAL, INFO )

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRS call (multilayer) in VLIDORT_SURFACE WFS'
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  Set Linearized integration constants NCON_ALB and PCON_ALB, all layer

        DO Q = 1, N_SURFACE_WFS
          DO N = 1, NLAYERS
            C0 = (N-1)*NSTKS_NSTRMS_2
            DO K = 1, K_REAL(N)
              IROW = K
              IROW1 = IROW + NSTKS_NSTRMS
              NCON_ALB(Q,K,N) = COL2_WFALB(C0+IROW,Q)
              PCON_ALB(Q,K,N) = COL2_WFALB(C0+IROW1,Q)
            ENDDO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2 * K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              IROW    = K + K_REAL(N)
              IROW1   = IROW + NSTKS_NSTRMS
              IROW_S  = IROW + K_COMPLEX(N)
              IROW1_S = IROW_S + NSTKS_NSTRMS
              NCON_ALB(Q,K1,N) = COL2_WFALB(C0+IROW,   Q)
              NCON_ALB(Q,K2,N) = COL2_WFALB(C0+IROW_S, Q)
              PCON_ALB(Q,K1,N) = COL2_WFALB(C0+IROW1,  Q)
              PCON_ALB(Q,K2,N) = COL2_WFALB(C0+IROW1_S,Q)
            ENDDO
          ENDDO
        ENDDO

!  Solve the boundary problem: No compression, Single Layer only
!  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2_WFALB

        CALL DGETRS &
           ( 'N', NTOTAL, N_SURFACE_WFS, SMAT2, MAXSTRMSTKS_2, &
              SIPIVOT, SCOL2_WFALB, MAXSTRMSTKS_2, INFO )

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRS call (Reg. 1 layer) in SURFACEWF_MASTER'
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  Set Linearized integration constants NCON_ALB and PCON_ALB, 1 layer

        DO Q = 1, N_SURFACE_WFS
          N = 1
          DO K = 1, K_REAL(N)
            IROW = K
            IROW1 = IROW + NSTKS_NSTRMS
            NCON_ALB(Q,K,N) = SCOL2_WFALB(IROW,Q)
            PCON_ALB(Q,K,N) = SCOL2_WFALB(IROW1,Q)
          ENDDO
          KO1 = K_REAL(N) + 1
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            IROW    = K + K_REAL(N)
            IROW1   = IROW + NSTKS_NSTRMS
            IROW_S  = IROW + K_COMPLEX(N)
            IROW1_S = IROW_S + NSTKS_NSTRMS
            NCON_ALB(Q,K1,N) = SCOL2_WFALB(IROW,   Q)
            NCON_ALB(Q,K2,N) = SCOL2_WFALB(IROW_S, Q)
            PCON_ALB(Q,K1,N) = SCOL2_WFALB(IROW1,  Q)
            PCON_ALB(Q,K2,N) = SCOL2_WFALB(IROW1_S,Q)
          ENDDO
        ENDDO

!  end clause

      ENDIF

!  debug------------------------------------------
!        if ( do_debug_write.and.fourier_component.eq.0 ) then
!         DO N = 1, NLAYERS
!          DO K = 1, K_REAL(N)
!           write(86,'(3i2,1p6e13.5)')FOURIER_COMPONENT,N,K,
!     &                LCON(K,N), MCON(K,N),
!     &                NCON_ALB(1,K,N),PCON_ALB(1,K,N),
!     &                NCON_ALB(1,K,N),PCON_ALB(1,K,N)
!          ENDDO
!         ENDDO
!        ENDIF


!  Get the Post-processed weighting functions
!  ==========================================

!  Upwelling weighting functions
!  -----------------------------

      IF ( DO_UPWELLING ) THEN

!  Get the surface term (L_BOA_SOURCE). External Function

        CALL BOA_SURFACEWF ( &
          DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFEMISS, &
          DO_MSMODE_THERMAL, DO_OBSERVATION_GEOMETRY, &
          DO_INCLUDE_MVOUTPUT, SURFACE_FACTOR, &
          FOURIER_COMPONENT, IBEAM, &
          DO_QUAD_OUTPUT, NSTOKES, &
          NSTREAMS, NLAYERS, &
          DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
          USER_BRDF_F, &
          SURFBB, DO_THERMAL_TRANSONLY, &
          N_SURFACE_WFS, LS_BRDF_F, &
          LS_USER_BRDF_F, LS_USER_BRDF_F_0, &
          LS_USER_EMISSIVITY, LS_EMISSIVITY, &
          DO_DBCORRECTION, FLUXVEC, &
          QUAD_STRMWTS, N_USER_STREAMS, &
          MUELLER_INDEX, DO_USER_STREAMS, &
          LOCAL_UM_START, T_DELT_EIGEN, &
          K_REAL, K_COMPLEX, &
          SOLA_XPOS, SOLB_XNEG, &
          ATMOS_ATTN, STOKES_DOWNSURF, &
          NCON_ALB, PCON_ALB, &
          LS_BOA_SOURCE, LS_BOA_THTONLY_SOURCE )

        CALL UPUSER_SURFACEWF ( &
          DO_OBSERVATION_GEOMETRY, &
          FLUX_MULTIPLIER, IBEAM, &
          N_SURFACE_WFS, LS_BOA_SOURCE, &
          NSTOKES, NLAYERS, &
          N_USER_LEVELS, N_USER_STREAMS, &
          DO_USER_STREAMS, LOCAL_UM_START, &
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
          UTAU_LEVEL_MASK_UP, PARTLAYERS_LAYERIDX, &
          T_DELT_USERM, T_UTUP_USERM, &
          DO_THERMAL_TRANSONLY, &
          K_REAL, K_COMPLEX, &
          UHOM_UPDN, UHOM_UPUP, &
          HMULT_1, HMULT_2, &
          NCON_ALB, PCON_ALB, &
          UT_HMULT_UU, UT_HMULT_UD, &
          SURFACEWF_F )
!mick temp fix 9/7/2012 - ELSE added
      ELSE
        IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
          SURFACEWF_F(1:N_SURFACE_WFS,1:N_USER_LEVELS,&
                      1:N_USER_STREAMS,IBEAM,1:NSTOKES,UPIDX) = ZERO
        ELSE
          SURFACEWF_F(1:N_SURFACE_WFS,1:N_USER_LEVELS,&
                      LUM,IBEAM,1:NSTOKES,UPIDX) = ZERO
        ENDIF
      ENDIF

!  Downwelling Albedo weighting functions
!  --------------------------------------

      IF ( DO_DNWELLING ) THEN
        CALL DNUSER_SURFACEWF ( &
          DO_OBSERVATION_GEOMETRY, &
          FLUX_MULTIPLIER, IBEAM, &
          N_SURFACE_WFS, NSTOKES, N_USER_LEVELS, &
          N_USER_STREAMS, DO_USER_STREAMS, &
          LOCAL_UM_START, &
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
          UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX, &
          T_DELT_USERM, T_UTDN_USERM, &
          DO_THERMAL_TRANSONLY, &
          K_REAL, K_COMPLEX, &
          UHOM_DNDN, UHOM_DNUP, &
          HMULT_1, HMULT_2, &
          NCON_ALB, PCON_ALB, &
          UT_HMULT_DU, UT_HMULT_DD, &
          SURFACEWF_F )
!mick temp fix 9/7/2012 - ELSE added
      ELSE
        IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
          SURFACEWF_F(1:N_SURFACE_WFS,1:N_USER_LEVELS,&
                      1:N_USER_STREAMS,IBEAM,1:NSTOKES,DNIDX) = ZERO
        ELSE
          SURFACEWF_F(1:N_SURFACE_WFS,1:N_USER_LEVELS,&
                      LUM,IBEAM,1:NSTOKES,DNIDX) = ZERO
        ENDIF
      ENDIF

!  mean value output
!  -----------------

      IF ( DO_INCLUDE_MVOUTPUT.OR.DO_QUAD_OUTPUT ) THEN
        CALL VLIDORT_LS_INTEGRATED_OUTPUT ( &
          DO_INCLUDE_MVOUTPUT, FLUX_MULTIPLIER, &
          IBEAM, N_SURFACE_WFS, LS_BOA_THTONLY_SOURCE, &
          NSTOKES, NSTREAMS, N_USER_LEVELS, &
          QUAD_WEIGHTS, QUAD_STRMWTS, &
          N_DIRECTIONS, WHICH_DIRECTIONS, &
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
          UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, &
          PARTLAYERS_LAYERIDX, &
          NLAYERS, DO_THERMAL_TRANSONLY, &
          T_DELT_DISORDS, T_DISORDS_UTUP, &
          T_UTUP_EIGEN, T_UTDN_EIGEN, &
          K_REAL, K_COMPLEX, &
          SOLA_XPOS, SOLB_XNEG, &
          NCON_ALB, PCON_ALB, &
          T_DELT_EIGEN, &
          MINT_SURFACEWF, FLUX_SURFACEWF )
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE SURFACEWF_MASTER

!

      SUBROUTINE SURFACEWF_COLSETUP ( &
        DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFEMISS, &
        SURFACE_FACTOR, FOURIER_COMPONENT, &
        IBEAM_INDEX, NSTOKES, NSTREAMS, &
        NLAYERS, DO_LAMBERTIAN_SURFACE, &
        SURFBB, N_SURFACE_WFS, LS_BRDF_F, &
        LS_BRDF_F_0, LS_EMISSIVITY, &
        QUAD_STRMWTS, NTOTAL, &
        NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
        MUELLER_INDEX, T_DELT_EIGEN, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        WLOWER, LCON, MCON, &
        ATMOS_ATTN, R2_HOMP, &
        R2_HOMM, R2_BEAM, DIRECT_BEAM, &
        COL2_WFALB, SCOL2_WFALB )

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_SURFEMISS
      DOUBLE PRECISION, INTENT (IN) :: SURFACE_FACTOR
      INTEGER, INTENT (IN) ::          FOURIER_COMPONENT
      INTEGER, INTENT (IN) ::          IBEAM_INDEX
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      DOUBLE PRECISION, INTENT (IN) :: SURFBB
      INTEGER, INTENT (IN) ::          N_SURFACE_WFS
      DOUBLE PRECISION, INTENT (IN) :: LS_BRDF_F &
          ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTOKES_SQ, &
            MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: LS_BRDF_F_0 &
          ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTOKES_SQ, &
            MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: LS_EMISSIVITY &
          ( MAX_SURFACEWFS, MAXSTOKES,MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::          NTOTAL
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS_2
      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: WLOWER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: ATMOS_ATTN ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: R2_HOMP &
          ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION, INTENT (IN) :: R2_HOMM &
          ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION, INTENT (IN) :: R2_BEAM ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: DIRECT_BEAM &
          ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )

      DOUBLE PRECISION, INTENT (OUT) :: COL2_WFALB &
          ( MAXTOTAL, MAX_SURFACEWFS )
      DOUBLE PRECISION, INTENT (OUT) :: SCOL2_WFALB &
          ( MAXSTRMSTKS_2, MAX_SURFACEWFS )

!  local variables
!  ---------------

      INTEGER ::          N, IB, Q, M, NSTOKES_ACTUAL
      INTEGER ::          K, KO1, K0, K1, K2
      INTEGER ::          I, J, O1, O2, OM, IR, IROW, C0, CM
      DOUBLE PRECISION :: L_BEAM, L_HOM_R, L_HOM_CR, H1R, H1I, REFL_B
      DOUBLE PRECISION :: FACTOR_BRDF, REFL_ATTN, AWF_DIRECT, AWF_EMISS
      DOUBLE PRECISION :: H_1,  H_2,  H_1_S,  H_2_S,  H1,  H2,  DBR
      DOUBLE PRECISION :: H_1_CR,   H_2_CR,   H_1_CI,   H_2_CI, EMISS_VAR
      DOUBLE PRECISION :: H_1_S_CR, H_2_S_CR, H_1_S_CI, H_2_S_CI

!  Help arrays

      DOUBLE PRECISION :: PV_W ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION :: HV_P ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION :: HV_M ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )

!  Ground level boundary condition
!  -------------------------------

!  Initialise

      IB = IBEAM_INDEX
      M = FOURIER_COMPONENT
      FACTOR_BRDF = SURFACE_FACTOR

!  initialise. Vitally necessary
!    We found that when commented out, the whole thing didn't work.
!    R. Spurr and V. Natraj, 20 january 2006

      DO Q = 1, N_SURFACE_WFS
        DO I = 1, NTOTAL
          COL2_WFALB(I,Q) = ZERO
        ENDDO
      ENDDO

!  If this is the surface

      N  = NLAYERS
      C0 = N*NSTKS_NSTRMS_2 - NSTKS_NSTRMS

!  Save some quantities
!  --------------------

!  Only want 1 component for Lambertian

      NSTOKES_ACTUAL = NSTOKES
      IF ( DO_LAMBERTIAN_SURFACE ) NSTOKES_ACTUAL = 1

!  (This is a repetition of earlier code and could be stored)

!  start loops

      DO J = 1, NSTREAMS
        DO O1 = 1, NSTOKES_ACTUAL

!  Beam

          PV_W(J,O1) = WLOWER(J,O1,N) * QUAD_STRMWTS(J)

!  real homogeneous solution contributions

          DO K = 1, K_REAL(N)
            H1 = SOLA_XPOS(J,O1,K,N)
            H2 = SOLB_XNEG(J,O1,K,N)
            HV_P(J,O1,K) = QUAD_STRMWTS(J)*H1
            HV_M(J,O1,K) = QUAD_STRMWTS(J)*H2
          ENDDO

!  Complex homogeneous solution contributions

          KO1 = K_REAL(N) + 1
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1 + 1
            HV_P(J,O1,K1) = QUAD_STRMWTS(J)* SOLA_XPOS(J,O1,K1,N)
            HV_P(J,O1,K2) = QUAD_STRMWTS(J)* SOLA_XPOS(J,O1,K2,N)
            HV_M(J,O1,K1) = QUAD_STRMWTS(J)* SOLB_XNEG(J,O1,K1,N)
            HV_M(J,O1,K2) = QUAD_STRMWTS(J)* SOLB_XNEG(J,O1,K2,N)
          ENDDO

!  End loops

        ENDDO
      ENDDO

!  Lambertian case
!  ===============

!  Skip if not flagged

      IF ( .not. DO_LAMBERTIAN_SURFACE ) goto 998

!  Only 1 weighting function. Q = 1.

      DO Q = 1, N_SURFACE_WFS
        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM   = C0 + IROW

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  OLD CODE, 10 January 2011
!  Beam contribution for this Kernel
!
!            L_BEAM = R2_BEAM(I,O1)
!
!  Real homogeneous solutions for this Kernel
!
!            L_HOM_R = ZERO
!            DO K = 1, K_REAL(N)
!              L_HOM_R = L_HOM_R &
!                + LCON(K,N) * R2_HOMP(I,O1,K) * T_DELT_EIGEN(K,N) &
!                + MCON(K,N) * R2_HOMM(I,O1,K)
!            ENDDO
!  Complex homogeneous solutions for this Kernel
!            L_HOM_CR = ZERO
!            KO1 = K_REAL(N) + 1
!            DO K = 1, K_COMPLEX(N)
!              K0 = 2*K-2
!              K1 = KO1 + K0
!              K2 = K1 + 1
!              H1R =   R2_HOMP(I,O1,K1) * T_DELT_EIGEN(K1,N) &
!                    - R2_HOMP(I,O1,K2) * T_DELT_EIGEN(K2,N)
!              H1I =   R2_HOMP(I,O1,K1) * T_DELT_EIGEN(K2,N) &
!                    + R2_HOMP(I,O1,K2) * T_DELT_EIGEN(K1,N)
!              L_HOM_CR = L_HOM_CR &
!                        + LCON(K1,N) * H1R  - LCON(K2,N) * H1I &
!                        + MCON(K1,N) * R2_HOMM(I,O1,K1) &
!                        - MCON(K2,N) * R2_HOMM(I,O1,K2)
!            ENDDO
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  @@@@@@@@@@@@@@@    START New Section   @@@@@@@@@@@@@@@@@@@@@@@@@

!  Beam contribution for this Kernel, only for NStokes = 1

            L_BEAM = ZERO
            IF ( O1 .eq. 1 ) THEN
              REFL_B = ZERO
              DO J = 1, NSTREAMS
                REFL_B = REFL_B + PV_W(J,O1)
              ENDDO
              L_BEAM = REFL_B * SURFACE_FACTOR
            ENDIF
!  @@@ Rob fix 17jan11, this L_BEAM statement was outside the IF loop
!            L_BEAM = REFL_B * SURFACE_FACTOR

!  Real homogeneous solutions contribution to this Kernel

            L_HOM_R = ZERO
            IF ( O1 .eq. 1 ) THEN
              DO K = 1, K_REAL(N)
                H_1 = ZERO
                H_2 = ZERO
                DO J = 1, NSTREAMS
                  H_1 = H_1 + HV_P(J,O1,K) * SURFACE_FACTOR
                  H_2 = H_2 + HV_M(J,O1,K) * SURFACE_FACTOR
                ENDDO
                L_HOM_R = L_HOM_R + LCON(K,N) * H_1 * T_DELT_EIGEN(K,N) &
                                  + MCON(K,N) * H_2
              ENDDO
            ENDIF

!  Complex homogeneous solutions for this Kernel

            L_HOM_CR = ZERO
            KO1 = K_REAL(N) + 1
            IF ( O1 .eq. 1 ) THEN
              DO K = 1, K_COMPLEX(N)
                K0 = 2 * K - 2
                K1 = KO1 + K0
                K2 = K1 + 1
                H_1_CR = ZERO
                H_2_CR = ZERO
                H_1_CI = ZERO
                H_2_CI = ZERO
                DO J = 1, NSTREAMS
                  H_1_CR = H_1_CR + HV_P(J,O1,K1) * SURFACE_FACTOR
                  H_2_CR = H_2_CR + HV_M(J,O1,K1) * SURFACE_FACTOR
                  H_1_CI = H_1_CI + HV_P(J,O1,K2) * SURFACE_FACTOR
                  H_2_CI = H_2_CI + HV_M(J,O1,K2) * SURFACE_FACTOR
                ENDDO
                H1R =   H_1_CR * T_DELT_EIGEN(K1,N) &
                      - H_1_CI * T_DELT_EIGEN(K2,N)
                H1I =   H_1_CR * T_DELT_EIGEN(K2,N) &
                      + H_1_CI * T_DELT_EIGEN(K1,N)
                L_HOM_CR = L_HOM_CR &
                  + LCON(K1,N) *   H1R  - LCON(K2,N) *   H1I  &
                  + MCON(K1,N) * H_2_CR - MCON(K2,N) * H_2_CI
              ENDDO
            ENDIF
!  @@@@@@@@@@@@@@@    End of New Section   @@@@@@@@@@@@@@@@@@@@@@@@@

!  Final contribution

            COL2_WFALB(CM,Q) = L_BEAM + L_HOM_R + L_HOM_CR

!  End Streams and Stokes loops

          ENDDO
        ENDDO

!  Add direct beam variation of albedo

        IF ( DO_INCLUDE_DIRECTBEAM ) THEN
          O1 = 1
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            IROW = IR + O1
            CM   = C0 + IROW
            COL2_WFALB(CM,Q) = COL2_WFALB(CM,Q) + ATMOS_ATTN(IB)
          ENDDO
        ENDIF

!  @@ Rob Fix: Use ATMOS_ATTN instead of Direct Beam (Unnormalized WF!!)
!              COL2_WFALB(CM,Q) = COL2_WFALB(CM,Q) + DIRECT_BEAM(I,IB,O1)
!              OOnly need O1 = 1 component (Lambertian)
!        IF ( DO_INCLUDE_DIRECTBEAM ) THEN
!          DO I = 1, NSTREAMS
!            IR = NSTOKES*(I-1)
!            DO O1 = 1, NSTOKES
!              IROW = IR + O1
!              CM   = C0 + IROW
!              COL2_WFALB(CM,Q) =  COL2_WFALB(CM,Q) + DIRECT_BEAM(I,IB,O1)
!            ENDDO
!          ENDDO
!        ENDIF

!  If surface emission, include emissivity variation
!    This code added for Version 2.4RT

        IF ( DO_INCLUDE_SURFEMISS ) THEN
          O1   = 1
          EMISS_VAR = SURFBB
          DO I = 1, NSTREAMS
            IR   = NSTOKES*(I-1)
            IROW = IR + O1
            CM   = C0 + IROW
            COL2_WFALB(CM,Q) = COL2_WFALB(CM,Q) - EMISS_VAR
          ENDDO
        ENDIF

!  @@@ Rob Fix FORMER CODE in F77 Version
!             Changed after input from V. Natraj 11/14/2010.
!  @@@ Rob Fix, exclude this line
!          IF ( DO_SOLAR_SOURCES ) EMISS_VAR = SURFBB * PI4
!  @@@ Rob Fix, do not multiply EMISS_VAR by albedo
!          EMISS_VAR = LAMBERTIAN_ALBEDO * EMISS_VAR
!        IF ( DO_INCLUDE_SURFEMISS ) THEN
!          EMISS_VAR = SURFBB
!          IF ( DO_SOLAR_SOURCES ) EMISS_VAR = SURFBB * PI4
!          EMISS_VAR = LAMBERTIAN_ALBEDO * EMISS_VAR
!          O1   = 1
!          DO I = 1, NSTREAMS
!            IR   = NSTOKES*(I-1)
!            IROW = IR + O1
!            CM   = C0 + IROW
!            COL2_WFALB(CM,Q) = COL2_WFALB(CM,Q) - SURFBB
!          ENDDO
!        ENDIF

!  Copy for the single layer case

        IF ( NLAYERS .EQ. 1 ) THEN
          DO N = 1, NTOTAL
            SCOL2_WFALB(N,Q) = COL2_WFALB(N,Q)
          ENDDO
        ENDIF

!  End parameter loop

      ENDDO

!  Return after finishing Lambertian case

      RETURN

!  BRDF boundary conditions
!  ========================

!  Continuation point

998   continue

!  Diffuse scatter contributions
!  -----------------------------

!  Start weighting function loop

      DO Q = 1, N_SURFACE_WFS

!  start loops

        DO I = 1, NSTREAMS
         IR = NSTOKES*(I-1)
         DO O1 = 1, NSTOKES
          IROW = IR + O1
          CM   = C0 + IROW

!  Beam contribution for this Kernel
!  @@@@@@@@@@@ Rob Fix, 2/9/11, (I,J) instead of (J,I)

          REFL_B = ZERO
          DO J = 1, NSTREAMS
            DO O2 = 1, NSTOKES
              OM = MUELLER_INDEX(O1,O2)
!              REFL_B = REFL_B + PV_W(J,O2) * LS_BRDF_F(Q,M,OM,J,I)
              REFL_B = REFL_B + PV_W(J,O2) * LS_BRDF_F(Q,M,OM,I,J)
            ENDDO
          ENDDO
          L_BEAM = REFL_B * FACTOR_BRDF

!  Real homogeneous solutions contribution to this Kernel
!  @@@@@@@@@@@ Rob Fix, 2/9/11, (I,J) instead of (J,I)

          L_HOM_R = ZERO
          DO K = 1, K_REAL(N)
            H_1 = ZERO
            H_2 = ZERO
            DO J = 1, NSTREAMS
              H_1_S = ZERO
              H_2_S = ZERO
              DO O2 = 1, NSTOKES
                OM  = MUELLER_INDEX(O1,O2)
!                H_1_S = H_1_S + HV_P(J,O2,K) * LS_BRDF_F(Q,M,OM,J,I)
!                H_2_S = H_2_S + HV_M(J,O2,K) * LS_BRDF_F(Q,M,OM,J,I)
                H_1_S = H_1_S + HV_P(J,O2,K) * LS_BRDF_F(Q,M,OM,I,J)
                H_2_S = H_2_S + HV_M(J,O2,K) * LS_BRDF_F(Q,M,OM,I,J)
              ENDDO
              H_1 = H_1 + H_1_S
              H_2 = H_2 + H_2_S
            ENDDO
            H_1 = FACTOR_BRDF * H_1
            H_2 = FACTOR_BRDF * H_2
            L_HOM_R = L_HOM_R + LCON(K,N) * H_1 * T_DELT_EIGEN(K,N) &
                              + MCON(K,N) * H_2
          ENDDO

!  homogeneous complex solutions

          L_HOM_CR = ZERO
          KO1 = K_REAL(N) + 1
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1 + 1
            H_1_CR = ZERO
            H_2_CR = ZERO
            H_1_CI = ZERO
            H_2_CI = ZERO
            DO J = 1, NSTREAMS
              H_1_S_CR = ZERO
              H_2_S_CR = ZERO
              H_1_S_CI = ZERO
              H_2_S_CI = ZERO
              DO O2 = 1, NSTOKES
                OM = MUELLER_INDEX(O1,O2)
                DBR = LS_BRDF_F(Q,M,OM,J,I)
                H_1_S_CR = H_1_S_CR + HV_P(J,O2,K1) * DBR
                H_2_S_CR = H_2_S_CR + HV_M(J,O2,K1) * DBR
                H_1_S_CI = H_1_S_CI + HV_P(J,O2,K2) * DBR
                H_2_S_CI = H_2_S_CI + HV_M(J,O2,K2) * DBR
              ENDDO
              H_1_CR = H_1_CR + H_1_S_CR
              H_2_CR = H_2_CR + H_2_S_CR
              H_1_CI = H_1_CI + H_1_S_CI
              H_2_CI = H_2_CI + H_2_S_CI
            ENDDO
            H_1_CR = FACTOR_BRDF  * H_1_CR
            H_1_CI = FACTOR_BRDF  * H_1_CI
            H_2_CR = FACTOR_BRDF  * H_2_CR
            H_2_CI = FACTOR_BRDF  * H_2_CI
            H1R =   H_1_CR * T_DELT_EIGEN(K1,N) &
                  - H_1_CI * T_DELT_EIGEN(K2,N)
            H1I =   H_1_CR * T_DELT_EIGEN(K2,N) &
                  + H_1_CI * T_DELT_EIGEN(K1,N)
            L_HOM_CR = L_HOM_CR &
              + LCON(K1,N) *   H1R  - LCON(K2,N) *   H1I &
              + MCON(K1,N) * H_2_CR - MCON(K2,N) * H_2_CI
          ENDDO

!  Final contribution

          COL2_WFALB(CM,Q) = L_BEAM + L_HOM_R + L_HOM_CR

!  End loops

         ENDDO
        ENDDO

!  Direct beam reflection

        IF ( DO_INCLUDE_DIRECTBEAM ) THEN
          REFL_ATTN = ATMOS_ATTN(IB)
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM   = C0 + IROW
              OM = MUELLER_INDEX(O1,1)
              AWF_DIRECT = REFL_ATTN * LS_BRDF_F_0(Q,M,OM,I,IB)
              COL2_WFALB(CM,Q) = COL2_WFALB(CM,Q) + AWF_DIRECT
            ENDDO
          ENDDO
        ENDIF

!  If surface emission, include emissivity variation, BRDF surface

        IF ( DO_INCLUDE_SURFEMISS ) THEN
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM   = C0 + IROW
              AWF_EMISS = SURFBB * LS_EMISSIVITY(Q,O1,I)
              COL2_WFALB(CM,Q) = COL2_WFALB(CM,Q) + AWF_EMISS
            ENDDO
          ENDDO
        ENDIF

!  Copy for the single layer case

        IF ( NLAYERS .EQ. 1 ) THEN
           DO N = 1, NTOTAL
            SCOL2_WFALB(N,Q) = COL2_WFALB(N,Q)
          ENDDO
        ENDIF

!  End parameter loop

      ENDDO

!  debug

!      if ( do_debug_write ) then
!        DO N = 1, NTOTAL
!          write(85,'(2i4,1p4e17.9)')IBEAM_INDEX, N, COL2_WFALB(N,1)
!        ENDDO
!        pause
!      ENDIF

!  Finish

      RETURN
      END SUBROUTINE SURFACEWF_COLSETUP

!

      SUBROUTINE BOA_SURFACEWF ( &
        DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFEMISS, &
        DO_MSMODE_THERMAL, DO_OBSERVATION_GEOMETRY, &
        DO_INCLUDE_MVOUTPUT, SURFACE_FACTOR, &
        FOURIER_COMPONENT, IBEAM, &
        DO_QUAD_OUTPUT, NSTOKES, &
        NSTREAMS, NLAYERS, &
        DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
        USER_BRDF_F, &
        SURFBB, DO_THERMAL_TRANSONLY, &
        N_SURFACE_WFS, LS_BRDF_F, &
        LS_USER_BRDF_F, LS_USER_BRDF_F_0, &
        LS_USER_EMISSIVITY, LS_EMISSIVITY, &
        DO_DBCORRECTION, FLUXVEC, &
        QUAD_STRMWTS, N_USER_STREAMS, &
        MUELLER_INDEX, DO_USER_STREAMS, &
        LOCAL_UM_START, T_DELT_EIGEN, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        ATMOS_ATTN, STOKES_DOWNSURF, &
        NCON_ALB, PCON_ALB, &
        LS_BOA_SOURCE, LS_BOA_THTONLY_SOURCE )

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_SURFEMISS
      LOGICAL, INTENT (IN) ::          DO_MSMODE_THERMAL
      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_MVOUTPUT
      DOUBLE PRECISION, INTENT (IN) :: SURFACE_FACTOR
      INTEGER, INTENT (IN) ::          FOURIER_COMPONENT
      INTEGER, INTENT (IN) ::          IBEAM
      LOGICAL, INTENT (IN) ::          DO_QUAD_OUTPUT
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      DOUBLE PRECISION, INTENT (IN) :: LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: USER_BRDF_F &
          ( 0:MAXMOMENTS,MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: SURFBB
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      INTEGER, INTENT (IN) ::          N_SURFACE_WFS
      DOUBLE PRECISION, INTENT (IN) :: LS_BRDF_F &
          ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTOKES_SQ, &
            MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: LS_USER_BRDF_F &
          ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTOKES_SQ, &
            MAX_USER_STREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: LS_USER_BRDF_F_0 &
          ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTOKES_SQ, &
            MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: LS_USER_EMISSIVITY &
           ( MAX_SURFACEWFS, MAXSTOKES, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: LS_EMISSIVITY &
           ( MAX_SURFACEWFS, MAXSTOKES, MAXSTREAMS )
      LOGICAL, INTENT (IN) ::          DO_DBCORRECTION
      DOUBLE PRECISION, INTENT (IN) :: FLUXVEC ( MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      LOGICAL, INTENT (IN) ::          DO_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: ATMOS_ATTN ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: STOKES_DOWNSURF ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: NCON_ALB &
          ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_ALB &
          ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )

      DOUBLE PRECISION, INTENT (OUT) :: LS_BOA_SOURCE &
          ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: LS_BOA_THTONLY_SOURCE &
          ( MAX_SURFACEWFS, MAXSTREAMS, MAXSTOKES )

!  Local variables
!  ---------------

      LOGICAL ::          DO_QTHTONLY
      INTEGER ::          UM, I, J, N, IB, O1, O2, M, OM, Q, LUM
      INTEGER ::          K, KO1, K0, K1, K2
      DOUBLE PRECISION :: INTEGRAND ( MAX_SURFACEWFS, MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION :: SUM_R, SUM_CR, REFLEC, S_REFLEC, REFL_ATTN
      DOUBLE PRECISION :: H1, H2, NXR, PXR, NXR1, NXR2, PXR1

!  Initialise
!  ----------

!  Shorthand

      N   = NLAYERS
      KO1 = K_REAL(N) + 1
      IB  = IBEAM
      M   = FOURIER_COMPONENT
      LUM = 1

!  Special flag

      DO_QTHTONLY = ( DO_THERMAL_TRANSONLY ) .AND. &
            ( DO_QUAD_OUTPUT .OR. DO_INCLUDE_MVOUTPUT )

!  initialise derivative of BOA source function

      IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
        DO Q = 1, N_SURFACE_WFS
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              LS_BOA_SOURCE(Q,UM,O1) = ZERO
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO Q = 1, N_SURFACE_WFS
          DO O1 = 1, NSTOKES
            LS_BOA_SOURCE(Q,IB,O1) = ZERO
          ENDDO
        ENDDO
      ENDIF

!  Thermal tranmsittance only, special term

      IF ( DO_QTHTONLY ) THEN
        DO Q = 1, N_SURFACE_WFS
          DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES
              LS_BOA_THTONLY_SOURCE(Q,I,O1) = ZERO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Skip diffuse-field variation for thermal transmittance-only

      IF ( DO_THERMAL_TRANSONLY .OR. .NOT. DO_USER_STREAMS ) GO TO 599

!  Contribution due to derivatives of BV constants
!  -----------------------------------------------

!  First compute derivative of downward intensity Integrand at stream an
!        .. reflectance integrand  = a(j).x(j).I_DOWN(-j)

!  start loops

      DO Q = 1, N_SURFACE_WFS
       DO I = 1, NSTREAMS
        DO O1 = 1, NSTOKES

!  Real homogeneous solutions

         SUM_R = ZERO
         DO K = 1, K_REAL(N)
          NXR = NCON_ALB(Q,K,N) * SOLA_XPOS(I,O1,K,N)
          PXR = PCON_ALB(Q,K,N) * SOLB_XNEG(I,O1,K,N)
          SUM_R = SUM_R + NXR*T_DELT_EIGEN(K,N) + PXR
         ENDDO

!  Complex solutions

         SUM_CR = ZERO
         DO K = 1, K_COMPLEX(N)
          K0 = 2 * K - 2
          K1 = KO1 + K0
          K2 = K1  + 1
          NXR1 =   NCON_ALB(Q,K1,N) * SOLA_XPOS(I,O1,K1,N) &
                 - NCON_ALB(Q,K2,N) * SOLA_XPOS(I,O1,K2,N)
          NXR2 =   NCON_ALB(Q,K1,N) * SOLA_XPOS(I,O1,K2,N) &
                 + NCON_ALB(Q,K2,N) * SOLA_XPOS(I,O1,K1,N)
          PXR1 =   PCON_ALB(Q,K1,N) * SOLB_XNEG(I,O1,K1,N) &
                 - PCON_ALB(Q,K2,N) * SOLB_XNEG(I,O1,K2,N)
          H1 =  NXR1 * T_DELT_EIGEN(K1,N) &
               -NXR2 * T_DELT_EIGEN(K2,N)
          H2 =  PXR1
          SUM_CR = SUM_CR + H1 + H2
         ENDDO

!  Final result

         INTEGRAND(Q,I,O1) = QUAD_STRMWTS(I) * ( SUM_R + SUM_CR )

!  end loops

        ENDDO
       ENDDO
      ENDDO

!  integrated reflectance term
!  ---------------------------

!  integrate Lambertian case, same for all user-streams

      IF ( DO_LAMBERTIAN_SURFACE ) THEN
       IF ( FOURIER_COMPONENT.EQ.0 ) THEN
        O1 = 1
        Q = 1
        REFLEC = ZERO
        DO J = 1, NSTREAMS
           REFLEC = REFLEC + INTEGRAND(Q,J,O1)
        ENDDO
        REFLEC = SURFACE_FACTOR * REFLEC * LAMBERTIAN_ALBEDO
        IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            LS_BOA_SOURCE(Q,UM,O1) = REFLEC
          ENDDO
        ELSE
          LS_BOA_SOURCE(Q,IB,O1) = REFLEC
        ENDIF
       ENDIF
      ENDIF

!  BRDF case

      IF ( .not. DO_LAMBERTIAN_SURFACE ) THEN
        IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
          DO Q = 1, N_SURFACE_WFS
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                REFLEC = ZERO
                DO J = 1, NSTREAMS
                  S_REFLEC = ZERO
                  DO O2 = 1, NSTOKES
                    OM = MUELLER_INDEX(O1,O2)
                    S_REFLEC = S_REFLEC + INTEGRAND(Q,J,O2) * &
                                   USER_BRDF_F(M,OM,UM,J)
                  ENDDO
                  REFLEC = REFLEC + S_REFLEC
                ENDDO
                LS_BOA_SOURCE(Q,UM,O1) = REFLEC * SURFACE_FACTOR
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO Q = 1, N_SURFACE_WFS
            DO O1 = 1, NSTOKES
              REFLEC = ZERO
              DO J = 1, NSTREAMS
                S_REFLEC = ZERO
                DO O2 = 1, NSTOKES
                  OM = MUELLER_INDEX(O1,O2)
                  S_REFLEC = S_REFLEC + INTEGRAND(Q,J,O2) * &
                                 USER_BRDF_F(M,OM,IB,J)
                ENDDO
                REFLEC = REFLEC + S_REFLEC
              ENDDO
              LS_BOA_SOURCE(Q,IB,O1) = REFLEC * SURFACE_FACTOR
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  Contributions due to direct variation of kernel parameter
!  ---------------------------------------------------------

!  Lambertian (this is the albedo)

      IF ( DO_LAMBERTIAN_SURFACE ) THEN
        IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
          O1 = 1
          Q  = 1
          IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              REFLEC = ZERO
              DO J = 1, NSTREAMS
                REFLEC = REFLEC + STOKES_DOWNSURF(J,O1)
              ENDDO
!  @@@ Rob Fix. Leave out Albedo ( Un-normalzied WFs )
!             REFLEC = REFLEC * LAMBERTIAN_ALBEDO * SURFACE_FACTOR
              REFLEC = REFLEC * SURFACE_FACTOR
              LS_BOA_SOURCE(Q,UM,O1) = LS_BOA_SOURCE(Q,UM,O1) + REFLEC
            ENDDO
          ELSE
            REFLEC = ZERO
            DO J = 1, NSTREAMS
              REFLEC = REFLEC + STOKES_DOWNSURF(J,O1)
            ENDDO
!  @@@ Rob Fix. Leave out Albedo ( Un-normalzied WFs )
!           REFLEC = REFLEC * LAMBERTIAN_ALBEDO * SURFACE_FACTOR
            REFLEC = REFLEC * SURFACE_FACTOR
            LS_BOA_SOURCE(Q,IB,O1) = LS_BOA_SOURCE(Q,IB,O1) + REFLEC
          ENDIF
        ENDIF
      ENDIF

!   Non-Lambertian (need derivative of BRDF Fourier term)

      IF ( .not.DO_LAMBERTIAN_SURFACE ) THEN
        IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
          DO Q = 1, N_SURFACE_WFS
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                REFLEC = ZERO
                DO J = 1, NSTREAMS
                  S_REFLEC = ZERO
                  DO O2 = 1, NSTOKES
                    OM = MUELLER_INDEX(O1,O2)
                    S_REFLEC = S_REFLEC + STOKES_DOWNSURF(J,O2) * &
                                      LS_USER_BRDF_F(Q,M,OM,UM,J)
                  ENDDO
                  REFLEC = REFLEC + S_REFLEC
                ENDDO
                REFLEC = REFLEC * SURFACE_FACTOR
                LS_BOA_SOURCE(Q,UM,O1)=LS_BOA_SOURCE(Q,UM,O1) + REFLEC
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO Q = 1, N_SURFACE_WFS
            DO O1 = 1, NSTOKES
              REFLEC = ZERO
              DO J = 1, NSTREAMS
                S_REFLEC = ZERO
                DO O2 = 1, NSTOKES
                  OM = MUELLER_INDEX(O1,O2)
                  S_REFLEC = S_REFLEC + STOKES_DOWNSURF(J,O2) * &
                                    LS_USER_BRDF_F(Q,M,OM,IB,J)
                ENDDO
                REFLEC = REFLEC + S_REFLEC
              ENDDO
              REFLEC = REFLEC * SURFACE_FACTOR
              LS_BOA_SOURCE(Q,IB,O1)=LS_BOA_SOURCE(Q,IB,O1) + REFLEC
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  Continuation point for avoiding diffuse field computation

 599  continue

!  Thermal tranmsittance and Integrated output (quadrature terms)
!     This is the reflectance term

      IF ( DO_QTHTONLY ) THEN
        IF ( DO_LAMBERTIAN_SURFACE ) THEN
          O1 = 1
          Q = 1
          REFLEC = ZERO
          DO J = 1, NSTREAMS
            REFLEC = REFLEC + STOKES_DOWNSURF(J,O1)
          ENDDO
          REFLEC = SURFACE_FACTOR * REFLEC
          DO I = 1, NSTREAMS
            LS_BOA_THTONLY_SOURCE(Q,I,O1) = &
                LS_BOA_THTONLY_SOURCE(Q,I,O1)  + REFLEC
          ENDDO
        ELSE
          DO Q = 1, N_SURFACE_WFS
            DO I = 1, NSTREAMS
              DO O1 = 1, NSTOKES
                REFLEC = ZERO
                DO J = 1, NSTREAMS
                  S_REFLEC = ZERO
                  DO O2 = 1, NSTOKES
                    OM = MUELLER_INDEX(O1,O2)
                    S_REFLEC = S_REFLEC + STOKES_DOWNSURF(J,O2) * &
                                    LS_BRDF_F(Q,M,OM,I,J)
                  ENDDO
                  REFLEC = REFLEC + S_REFLEC
                ENDDO
                REFLEC = REFLEC * SURFACE_FACTOR
                LS_BOA_THTONLY_SOURCE(Q,I,O1) = &
                  LS_BOA_THTONLY_SOURCE(Q,I,O1)  + REFLEC
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  Add linearization term for variation of direct beam reflectance

      go to 567

!  @@@ Rob Fix, add DO_USER_STREAMS to If Block
      IF ( DO_USER_STREAMS ) THEN
       IF ( DO_INCLUDE_DIRECTBEAM.and..not.DO_DBCORRECTION ) THEN
        REFL_ATTN = ATMOS_ATTN(IB)
        IF ( DO_LAMBERTIAN_SURFACE.AND.M.EQ.0) THEN
          O1 = 1
          Q  = 1
!  @@@ Rob Fix. Leave out Albedo ( Un-normalzied WFs )
!          REFLEC = REFL_ATTN * LAMBERTIAN_ALBEDO
          REFLEC = REFL_ATTN
          IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              LS_BOA_SOURCE(Q,UM,O1) = LS_BOA_SOURCE(Q,UM,O1) + REFLEC
            ENDDO
          ELSE
            LS_BOA_SOURCE(Q,IB,O1) = LS_BOA_SOURCE(Q,IB,O1) + REFLEC
          ENDIF
        ELSE IF ( .not. DO_LAMBERTIAN_SURFACE ) THEN
          IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
            DO Q = 1, N_SURFACE_WFS
!mick fix 1/3/2013 - added UM loop
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO O1 = 1, NSTOKES
                  REFLEC = ZERO
                  DO O2 = 1, NSTOKES
                    OM = MUELLER_INDEX(O1,O2)
                    REFLEC = REFLEC + &
                       LS_USER_BRDF_F_0(Q,M,OM,UM,IB) * FLUXVEC(O2)
                  ENDDO
                  REFLEC = REFLEC * REFL_ATTN
                  LS_BOA_SOURCE(Q,UM,O1) = LS_BOA_SOURCE(Q,UM,O1) + REFLEC
                ENDDO
              ENDDO
            ENDDO
          ELSE
            DO Q = 1, N_SURFACE_WFS
              DO O1 = 1, NSTOKES
                REFLEC = ZERO
                DO O2 = 1, NSTOKES
                  OM = MUELLER_INDEX(O1,O2)
                  REFLEC = REFLEC + &
                     LS_USER_BRDF_F_0(Q,M,OM,LUM,IB) * FLUXVEC(O2)
                ENDDO
                REFLEC = REFLEC * REFL_ATTN
                LS_BOA_SOURCE(Q,IB,O1) = LS_BOA_SOURCE(Q,IB,O1) + REFLEC
              ENDDO
            ENDDO
          ENDIF
        ENDIF
       ENDIF
      ENDIF    ! @@@ If (DO_USER_STREAM) clause, Added January 2011

 567 continue

!  Add emissivity variation at user defined angles.
!    Apparenly only present for Fourier zero
!  (expression for emissivity variation follows from Kirchhoff's law)

!  @@@ Rob Fix, add DO_USER_STREAMS to If Block
      IF ( DO_USER_STREAMS ) THEN
       IF ( DO_INCLUDE_SURFEMISS .and. M.eq.0 .and..not.DO_MSMODE_THERMAL ) THEN
        IF ( DO_LAMBERTIAN_SURFACE ) THEN
          O1 = 1
          Q  = 1
          IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              LS_BOA_SOURCE(Q,UM,O1) = LS_BOA_SOURCE(Q,UM,O1) - SURFBB
            ENDDO
          ELSE
            LS_BOA_SOURCE(Q,IB,O1) = LS_BOA_SOURCE(Q,IB,O1) - SURFBB
          ENDIF
        ELSE
          IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
            DO Q = 1, N_SURFACE_WFS
              DO O1 = 1, NSTOKES
                DO UM = LOCAL_UM_START, N_USER_STREAMS
                  LS_BOA_SOURCE(Q,UM,O1) = LS_BOA_SOURCE(Q,UM,O1) &
                              - SURFBB * LS_USER_EMISSIVITY(Q,O1,UM)
                ENDDO
              ENDDO
            ENDDO
          ELSE
            DO Q = 1, N_SURFACE_WFS
              DO O1 = 1, NSTOKES
                LS_BOA_SOURCE(Q,IB,O1) = LS_BOA_SOURCE(Q,IB,O1) &
                            - SURFBB * LS_USER_EMISSIVITY(Q,O1,IB)
              ENDDO
            ENDDO
          ENDIF
        ENDIF
       ENDIF
      ENDIF    ! @@@ If (DO_USER_STREAM) clause, Added January 2011

!  Thermal Transmittance only (For Integrated product)
!    Surface emissivity term

      IF ( DO_INCLUDE_SURFEMISS .and. DO_QTHTONLY ) THEN
        IF ( DO_LAMBERTIAN_SURFACE ) THEN
          O1 = 1
          Q = 1
          DO I = 1, NSTREAMS
            LS_BOA_THTONLY_SOURCE(Q,I,O1) = &
             LS_BOA_THTONLY_SOURCE(Q,I,O1) -  SURFBB
          ENDDO
        ELSE
          DO Q = 1, N_SURFACE_WFS
            DO O1 = 1, NSTOKES
              DO I = 1, NSTREAMS
                LS_BOA_THTONLY_SOURCE(Q,I,O1) = &
                      LS_BOA_THTONLY_SOURCE(Q,I,O1) &
                         - SURFBB * LS_EMISSIVITY(Q,O1,I)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BOA_SURFACEWF

!

      SUBROUTINE UPUSER_SURFACEWF ( &
        DO_OBSERVATION_GEOMETRY, &
        FLUX_MULTIPLIER, IBEAM, &
        N_SURFACE_WFS, LS_BOA_SOURCE, &
        NSTOKES, NLAYERS, &
        N_USER_LEVELS, N_USER_STREAMS, &
        DO_USER_STREAMS, LOCAL_UM_START, &
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
        UTAU_LEVEL_MASK_UP, PARTLAYERS_LAYERIDX, &
        T_DELT_USERM, T_UTUP_USERM, &
        DO_THERMAL_TRANSONLY, &
        K_REAL, K_COMPLEX, &
        UHOM_UPDN, UHOM_UPUP, &
        HMULT_1, HMULT_2, &
        NCON_ALB, PCON_ALB, &
        UT_HMULT_UU, UT_HMULT_UD, &
        SURFACEWF_F )

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER
      INTEGER, INTENT (IN) ::          IBEAM
      INTEGER, INTENT (IN) ::          N_SURFACE_WFS
      DOUBLE PRECISION, INTENT (IN) :: LS_BOA_SOURCE &
          ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      LOGICAL, INTENT (IN) ::          DO_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPDN &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPUP &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_1 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_2 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: NCON_ALB &
          ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_ALB &
          ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_UU &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_UD &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )

      DOUBLE PRECISION, INTENT (INOUT) :: SURFACEWF_F &
         ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_VZANGLES, &
           MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  local variables
!  ---------------

      INTEGER ::          N, NUT, NSTART, NUT_PREV, NLEVEL, O1
      INTEGER ::          UTA, UM, Q, UT, IB, LUM

      DOUBLE PRECISION :: LS_CUMUL_SOURCE &
          ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION :: LS_LAYER_SOURCE &
          ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION :: LS_FINAL_SOURCE

!  Local indices

      IB  = IBEAM
      LUM = 1

!  Zero all Fourier components - New rule, better for safety
!    Only did this for components close to zenith (formerly)

      IF ( DO_USER_STREAMS ) THEN
        DO UTA = 1, N_USER_LEVELS
          IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
            DO UM = 1, N_USER_STREAMS
              DO Q = 1, N_SURFACE_WFS
                DO O1 = 1, NSTOKES
                  SURFACEWF_F(Q,UTA,UM,IB,O1,UPIDX) = ZERO
                ENDDO
              ENDDO
            ENDDO
          ELSE
            DO Q = 1, N_SURFACE_WFS
              DO O1 = 1, NSTOKES
                SURFACEWF_F(Q,UTA,LUM,IB,O1,UPIDX) = ZERO
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDIF

!  Initialize post-processing recursion
!  ====================================

      IF ( DO_USER_STREAMS ) THEN

!  Set the cumulative source term equal to the BOA sum

        IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
           DO Q = 1, N_SURFACE_WFS
            DO O1 = 1, NSTOKES
             LS_CUMUL_SOURCE(Q,UM,O1) = LS_BOA_SOURCE(Q,UM,O1)
            ENDDO
           ENDDO
          ENDDO
        ELSE
          DO Q = 1, N_SURFACE_WFS
           DO O1 = 1, NSTOKES
            LS_CUMUL_SOURCE(Q,IB,O1) = LS_BOA_SOURCE(Q,IB,O1)
           ENDDO
          ENDDO
        ENDIF

      ENDIF

!  Recursion Loop for linearized Post-processing
!  =============================================

!  initialise cumulative source term loop

      NUT = 0
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = N_USER_LEVELS, 1, -1

!  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_UP(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only
!    1. Get layer source terms
!    2. Find cumulative source term
!    3. Set multiple scatter source term (MSST) output if flagged

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL + 1
          DO N = NSTART, NUT, -1

            CALL LS_WHOLELAYER_STERM_UP ( &
              DO_OBSERVATION_GEOMETRY, &
              N_SURFACE_WFS, IB, N, &
              NSTOKES, DO_THERMAL_TRANSONLY, &
              N_USER_STREAMS, LOCAL_UM_START, &
              K_REAL, K_COMPLEX, &
              UHOM_UPDN, UHOM_UPUP, &
              HMULT_1, HMULT_2, &
              NCON_ALB, PCON_ALB, &
              LS_LAYER_SOURCE )

            IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO Q = 1, N_SURFACE_WFS
                  DO O1 = 1, NSTOKES
                    LS_CUMUL_SOURCE(Q,UM,O1) = LS_LAYER_SOURCE(Q,UM,O1) &
                       + T_DELT_USERM(N,UM) * LS_CUMUL_SOURCE(Q,UM,O1)
                  ENDDO
               ENDDO
              ENDDO
            ELSE
              DO Q = 1, N_SURFACE_WFS
                DO O1 = 1, NSTOKES
                  LS_CUMUL_SOURCE(Q,IB,O1) = LS_LAYER_SOURCE(Q,IB,O1) &
                     + T_DELT_USERM(N,IB) * LS_CUMUL_SOURCE(Q,IB,O1)
                ENDDO
              ENDDO
            ENDIF

          ENDDO
        ENDIF

!  Offgrid output
!  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)

!  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN
            CALL LS_PARTLAYER_STERM_UP ( &
              DO_OBSERVATION_GEOMETRY, &
              N_SURFACE_WFS, IB, UT, N, &
              NSTOKES, DO_THERMAL_TRANSONLY, &
              N_USER_STREAMS, LOCAL_UM_START, &
              K_REAL, K_COMPLEX, &
              UHOM_UPDN, UHOM_UPUP, &
              UT_HMULT_UU, UT_HMULT_UD, &
              NCON_ALB, PCON_ALB, &
              LS_LAYER_SOURCE )

            IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO Q = 1, N_SURFACE_WFS
                  DO O1 = 1, NSTOKES
                    LS_FINAL_SOURCE = LS_LAYER_SOURCE(Q,UM,O1) &
                  + T_UTUP_USERM(UT,UM) * LS_CUMUL_SOURCE(Q,UM,O1)
                    SURFACEWF_F(Q,UTA,UM,IB,O1,UPIDX) = &
                          FLUX_MULTIPLIER * LS_FINAL_SOURCE
                  ENDDO
                ENDDO
              ENDDO
            ELSE
              DO Q = 1, N_SURFACE_WFS
                DO O1 = 1, NSTOKES
                  LS_FINAL_SOURCE = LS_LAYER_SOURCE(Q,IB,O1) &
                + T_UTUP_USERM(UT,IB) * LS_CUMUL_SOURCE(Q,IB,O1)
                  SURFACEWF_F(Q,UTA,LUM,IB,O1,UPIDX) = &
                         FLUX_MULTIPLIER * LS_FINAL_SOURCE
                ENDDO
              ENDDO
            ENDIF

          ENDIF

!  Ongrid output
!  -------------

        ELSE

!  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
            IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO Q = 1, N_SURFACE_WFS
                  DO O1 = 1, NSTOKES
                    SURFACEWF_F(Q,UTA,UM,IB,O1,UPIDX) = &
                             FLUX_MULTIPLIER * LS_CUMUL_SOURCE(Q,UM,O1)
                  ENDDO
                ENDDO
              ENDDO
            ELSE
              DO Q = 1, N_SURFACE_WFS
                DO O1 = 1, NSTOKES
                  SURFACEWF_F(Q,UTA,LUM,IB,O1,UPIDX) = &
                           FLUX_MULTIPLIER * LS_CUMUL_SOURCE(Q,IB,O1)
                ENDDO
              ENDDO
            ENDIF
          ENDIF

        ENDIF

!  Check for updating the recursion

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT
        ENDIF

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE UPUSER_SURFACEWF

!

      SUBROUTINE DNUSER_SURFACEWF ( &
        DO_OBSERVATION_GEOMETRY, &
        FLUX_MULTIPLIER, IBEAM, &
        N_SURFACE_WFS, NSTOKES, N_USER_LEVELS, &
        N_USER_STREAMS, DO_USER_STREAMS, &
        LOCAL_UM_START, &
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
        UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX, &
        T_DELT_USERM, T_UTDN_USERM, &
        DO_THERMAL_TRANSONLY, &
        K_REAL, K_COMPLEX, &
        UHOM_DNDN, UHOM_DNUP, &
        HMULT_1, HMULT_2, &
        NCON_ALB, PCON_ALB, &
        UT_HMULT_DU, UT_HMULT_DD, &
        SURFACEWF_F )

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY
      INTEGER, INTENT (IN) ::          IBEAM
      INTEGER, INTENT (IN) ::          N_SURFACE_WFS
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      LOGICAL, INTENT (IN) ::          DO_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNDN &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNUP &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_1 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_2 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: NCON_ALB &
          ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_ALB &
          ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_DU &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_DD &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )

      DOUBLE PRECISION, INTENT (INOUT) :: SURFACEWF_F &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_VZANGLES, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  local variables
!  ---------------

      INTEGER ::          N, NUT, NSTART, NUT_PREV, NLEVEL, O1
      INTEGER ::          UTA, UM, Q, UT, IB, LUM

      DOUBLE PRECISION :: LS_CUMUL_SOURCE &
          ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION :: LS_LAYER_SOURCE &
          ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION :: LS_TOA_SOURCE &
          ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION :: LS_FINAL_SOURCE

!  Initialise

      IB  = IBEAM
      LUM = 1

!  Zero all Fourier component output

      IF ( DO_USER_STREAMS ) THEN
        DO UTA = 1, N_USER_LEVELS
          IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
!mick fix 1/8/2013 - changed limit of do loop
            !DO UM = 1, LOCAL_UM_START
            DO UM = 1, N_USER_STREAMS
              DO Q = 1, N_SURFACE_WFS
                DO O1 = 1, NSTOKES
                  SURFACEWF_F(Q,UTA,UM,IB,O1,DNIDX) = ZERO
                ENDDO
              ENDDO
            ENDDO
          ELSE
            DO Q = 1, N_SURFACE_WFS
              DO O1 = 1, NSTOKES
                SURFACEWF_F(Q,UTA,LUM,IB,O1,DNIDX) = ZERO
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDIF

!  Initialize post-processing recursion
!  ====================================

!  Get the linearized TOA source terms

      IF ( DO_USER_STREAMS ) THEN
        IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO Q = 1, N_SURFACE_WFS
              DO O1 = 1, NSTOKES
                LS_TOA_SOURCE(Q,UM,O1)   = ZERO
                LS_CUMUL_SOURCE(Q,UM,O1) = LS_TOA_SOURCE(Q,UM,O1)
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO Q = 1, N_SURFACE_WFS
            DO O1 = 1, NSTOKES
              LS_TOA_SOURCE(Q,IB,O1)   = ZERO
              LS_CUMUL_SOURCE(Q,IB,O1) = LS_TOA_SOURCE(Q,IB,O1)
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  Recursion Loop for linearized Post-processing
!  =============================================

!  initialise cumulative source term loop

      NUT = 0
      NSTART = 1
      NUT_PREV = NSTART - 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = 1, N_USER_LEVELS

!  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_DN(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only
!    1. Get layer source terms
!    2. Find cumulative source term
!    3. Set multiple scatter source term output if flagged

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL
          DO N = NSTART, NUT
            CALL LS_WHOLELAYER_STERM_DN ( &
              DO_OBSERVATION_GEOMETRY, &
              N_SURFACE_WFS, IB, N, &
              NSTOKES, DO_THERMAL_TRANSONLY, &
              N_USER_STREAMS, LOCAL_UM_START, &
              K_REAL, K_COMPLEX, &
              UHOM_DNDN, UHOM_DNUP, &
              HMULT_1, HMULT_2, &
              NCON_ALB, PCON_ALB, &
              LS_LAYER_SOURCE )

            IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO Q = 1, N_SURFACE_WFS
                  DO O1 = 1, NSTOKES
                    LS_CUMUL_SOURCE(Q,UM,O1) = LS_LAYER_SOURCE(Q,UM,O1) &
                     + T_DELT_USERM(N,UM) * LS_CUMUL_SOURCE(Q,UM,O1)
                  ENDDO
                ENDDO
              ENDDO
            ELSE
              DO Q = 1, N_SURFACE_WFS
                DO O1 = 1, NSTOKES
                  LS_CUMUL_SOURCE(Q,IB,O1) = LS_LAYER_SOURCE(Q,IB,O1) &
                       + T_DELT_USERM(N,IB) * LS_CUMUL_SOURCE(Q,IB,O1)
                ENDDO
              ENDDO
            ENDIF

          ENDDO
        ENDIF

!  Offgrid output
!  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)

!  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN
            CALL LS_PARTLAYER_STERM_DN ( &
              DO_OBSERVATION_GEOMETRY, &
              N_SURFACE_WFS, IB, UT, N, &
              NSTOKES, DO_THERMAL_TRANSONLY, &
              N_USER_STREAMS, LOCAL_UM_START, &
              K_REAL, K_COMPLEX, &
              UHOM_DNDN, UHOM_DNUP, &
              UT_HMULT_DU, UT_HMULT_DD, &
              NCON_ALB, PCON_ALB, &
              LS_LAYER_SOURCE )

            IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO Q = 1, N_SURFACE_WFS
                  DO O1 = 1, NSTOKES
                    LS_FINAL_SOURCE = LS_LAYER_SOURCE(Q,UM,O1) &
                     + T_UTDN_USERM(UT,UM) * LS_CUMUL_SOURCE(Q,UM,O1)
                    SURFACEWF_F(Q,UTA,UM,IB,O1,DNIDX) = &
                       FLUX_MULTIPLIER * LS_FINAL_SOURCE
                 ENDDO
               ENDDO
              ENDDO
            ELSE
              DO Q = 1, N_SURFACE_WFS
                DO O1 = 1, NSTOKES
                  LS_FINAL_SOURCE = LS_LAYER_SOURCE(Q,IB,O1) &
                   + T_UTDN_USERM(UT,IB) * LS_CUMUL_SOURCE(Q,IB,O1)
                  SURFACEWF_F(Q,UTA,LUM,IB,O1,DNIDX) = &
                        FLUX_MULTIPLIER * LS_FINAL_SOURCE
                ENDDO
              ENDDO
            ENDIF
          ENDIF

!  Ongrid output
!  -------------

        ELSE

!  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
            DO Q = 1, N_SURFACE_WFS
              IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
                DO UM = LOCAL_UM_START, N_USER_STREAMS
                  DO O1 = 1, NSTOKES
                    SURFACEWF_F(Q,UTA,UM,IB,O1,DNIDX) = &
                            FLUX_MULTIPLIER * LS_CUMUL_SOURCE(Q,UM,O1)
                  ENDDO
                ENDDO
              ELSE
                DO O1 = 1, NSTOKES
                  SURFACEWF_F(Q,UTA,LUM,IB,O1,DNIDX) = &
                          FLUX_MULTIPLIER * LS_CUMUL_SOURCE(Q,IB,O1)
                ENDDO
              ENDIF
            ENDDO
          ENDIF

        ENDIF

!  Check for updating the recursion

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT.NE. NUT_PREV ) NSTART = NUT + 1
          NUT_PREV = NUT
        ENDIF

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE DNUSER_SURFACEWF

!

      SUBROUTINE LS_WHOLELAYER_STERM_UP ( &
        DO_OBSERVATION_GEOMETRY, &
        N_SURFACE_WFS, IBEAM, GIVEN_LAYER, &
        NSTOKES, DO_THERMAL_TRANSONLY, &
        N_USER_STREAMS, LOCAL_UM_START, &
        K_REAL, K_COMPLEX, &
        UHOM_UPDN, UHOM_UPUP, &
        HMULT_1, HMULT_2, &
        NCON_ALB, PCON_ALB, &
        LS_LAYERSOURCE )

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY
      INTEGER, INTENT (IN) ::          N_SURFACE_WFS
      INTEGER, INTENT (IN) ::          IBEAM, GIVEN_LAYER
      INTEGER, INTENT (IN) ::          NSTOKES
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPDN &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPUP &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_1 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_2 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: NCON_ALB &
          ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_ALB &
          ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )

      DOUBLE PRECISION, INTENT (OUT) :: LS_LAYERSOURCE &
          ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          N, UM, O1, IB, K, KO1, K0, K1, K2, Q
      INTEGER ::          UM_START, UM_END
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, H1, H2
      DOUBLE PRECISION :: NUXR, PUXR, NUXR1, NUXR2, PUXR1, PUXR2

!  local indices

      N   = GIVEN_LAYER
      KO1 = K_REAL(N) + 1
      IB  = IBEAM

      IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
        UM_START = LOCAL_UM_START
        UM_END   = N_USER_STREAMS
      ELSE
        UM_START = IB
        UM_END   = IB
      ENDIF

!  Thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) THEN
        DO Q = 1, N_SURFACE_WFS
          DO O1 = 1, NSTOKES
            DO UM = UM_START,UM_END
              LS_LAYERSOURCE(Q,UM,O1) = ZERO
            ENDDO
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Homogeneous solutions
!  =====================

!  Loops over user angles, weighting functions and Stokes

      DO UM = UM_START,UM_END
       DO Q = 1, N_SURFACE_WFS
        DO O1 = 1, NSTOKES

!  Real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, K_REAL(N)
            NUXR = NCON_ALB(Q,K,N)*UHOM_UPDN(UM,O1,K,N)
            PUXR = PCON_ALB(Q,K,N)*UHOM_UPUP(UM,O1,K,N)
            H1 =  NUXR * HMULT_2(K,UM,N)
            H2 =  PUXR * HMULT_1(K,UM,N)
            SHOM_R = SHOM_R + H1 + H2
          ENDDO

!  Complex homogeneous solutions

          SHOM_CR = ZERO
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            NUXR1 =   NCON_ALB(Q,K1,N) * UHOM_UPDN(UM,O1,K1,N) &
                    - NCON_ALB(Q,K2,N) * UHOM_UPDN(UM,O1,K2,N)
            NUXR2 =   NCON_ALB(Q,K1,N) * UHOM_UPDN(UM,O1,K2,N) &
                    + NCON_ALB(Q,K2,N) * UHOM_UPDN(UM,O1,K1,N)
            PUXR1 =   PCON_ALB(Q,K1,N) * UHOM_UPUP(UM,O1,K1,N) &
                    - PCON_ALB(Q,K2,N) * UHOM_UPUP(UM,O1,K2,N)
            PUXR2 =   PCON_ALB(Q,K1,N) * UHOM_UPUP(UM,O1,K2,N) &
                    + PCON_ALB(Q,K2,N) * UHOM_UPUP(UM,O1,K1,N)
            H1 =   NUXR1 * HMULT_2(K1,UM,N) &
                 - NUXR2 * HMULT_2(K2,UM,N)
            H2 =   PUXR1 * HMULT_1(K1,UM,N) &
                 - PUXR2 * HMULT_1(K2,UM,N)
            SHOM_CR = SHOM_CR + H1 + H2
          ENDDO

!  homogeneous contribution

          LS_LAYERSOURCE(Q,UM,O1) = SHOM_R + SHOM_CR

!  End loops over Q, O1 and UM

        ENDDO
       ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LS_WHOLELAYER_STERM_UP

!

      SUBROUTINE LS_WHOLELAYER_STERM_DN ( &
        DO_OBSERVATION_GEOMETRY, &
        N_SURFACE_WFS, IBEAM, GIVEN_LAYER, &
        NSTOKES, DO_THERMAL_TRANSONLY, &
        N_USER_STREAMS, LOCAL_UM_START, &
        K_REAL, K_COMPLEX, &
        UHOM_DNDN, UHOM_DNUP, &
        HMULT_1, HMULT_2, &
        NCON_ALB, PCON_ALB, &
        LS_LAYERSOURCE )

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY
      INTEGER, INTENT (IN) ::          N_SURFACE_WFS
      INTEGER, INTENT (IN) ::          IBEAM, GIVEN_LAYER
      INTEGER, INTENT (IN) ::          NSTOKES
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNDN &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNUP &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_1 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_2 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: NCON_ALB &
          ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_ALB &
          ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )

      DOUBLE PRECISION, INTENT (OUT) :: LS_LAYERSOURCE &
          ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          N, UM, O1, IB, K, KO1, K0, K1, K2, Q
      INTEGER ::          UM_START, UM_END
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, H1, H2
      DOUBLE PRECISION :: NUXR, PUXR, NUXR1, NUXR2, PUXR1, PUXR2

!  Local indices

      N   = GIVEN_LAYER
      KO1 = K_REAL(N) + 1
      IB  = IBEAM

      IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
        UM_START = LOCAL_UM_START
        UM_END   = N_USER_STREAMS
      ELSE
        UM_START = IB
        UM_END   = IB
      ENDIF

!  Thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) THEN
        DO Q = 1, N_SURFACE_WFS
          DO O1 = 1, NSTOKES
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              LS_LAYERSOURCE(Q,UM,O1) = ZERO
            ENDDO
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Homogeneous solutions
!  =====================

!  Loop over user angles, weighting functions and STokes

      DO UM = LOCAL_UM_START, N_USER_STREAMS
       DO Q = 1, N_SURFACE_WFS
        DO O1 = 1, NSTOKES

!  Real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, K_REAL(N)
            NUXR = NCON_ALB(Q,K,N)*UHOM_DNDN(UM,O1,K,N)
            PUXR = PCON_ALB(Q,K,N)*UHOM_DNUP(UM,O1,K,N)
            H1 =  NUXR * HMULT_1(K,UM,N)
            H2 =  PUXR * HMULT_2(K,UM,N)
            SHOM_R = SHOM_R + H1 + H2
          ENDDO

!  Complex homogeneous solutions

          SHOM_CR = ZERO
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            NUXR1 =   NCON_ALB(Q,K1,N) * UHOM_DNDN(UM,O1,K1,N) &
                    - NCON_ALB(Q,K2,N) * UHOM_DNDN(UM,O1,K2,N)
            NUXR2 =   NCON_ALB(Q,K1,N) * UHOM_DNDN(UM,O1,K2,N) &
                    + NCON_ALB(Q,K2,N) * UHOM_DNDN(UM,O1,K1,N)
            PUXR1 =   PCON_ALB(Q,K1,N) * UHOM_DNUP(UM,O1,K1,N) &
                    - PCON_ALB(Q,K2,N) * UHOM_DNUP(UM,O1,K2,N)
            PUXR2 =   PCON_ALB(Q,K1,N) * UHOM_DNUP(UM,O1,K2,N) &
                    + PCON_ALB(Q,K2,N) * UHOM_DNUP(UM,O1,K1,N)
            H1 =   NUXR1 * HMULT_1(K1,UM,N) &
                 - NUXR2 * HMULT_1(K2,UM,N)
            H2 =   PUXR1 * HMULT_2(K1,UM,N) &
                 - PUXR2 * HMULT_2(K2,UM,N)
            SHOM_CR = SHOM_CR + H1 + H2
          ENDDO

!  homogeneous contribution

          LS_LAYERSOURCE(Q,UM,O1) = SHOM_R + SHOM_CR

!  End loops over Q, O1 and UM

        ENDDO
       ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LS_WHOLELAYER_STERM_DN

!

      SUBROUTINE LS_PARTLAYER_STERM_UP ( &
        DO_OBSERVATION_GEOMETRY, &
        N_SURFACE_WFS, IBEAM, &
        OFFGRID_INDEX, GIVEN_LAYER, &
        NSTOKES, DO_THERMAL_TRANSONLY, &
        N_USER_STREAMS, LOCAL_UM_START, &
        K_REAL, K_COMPLEX, &
        UHOM_UPDN, UHOM_UPUP, &
        UT_HMULT_UU, UT_HMULT_UD, &
        NCON_ALB, PCON_ALB, &
        LS_LAYERSOURCE )

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY
      INTEGER, INTENT (IN) ::          N_SURFACE_WFS
      INTEGER, INTENT (IN) ::          GIVEN_LAYER, IBEAM
      INTEGER, INTENT (IN) ::          OFFGRID_INDEX
      INTEGER, INTENT (IN) ::          NSTOKES
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPDN &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPUP &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_UU &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_UD &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: NCON_ALB &
          ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_ALB &
          ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )

      DOUBLE PRECISION, INTENT (OUT) :: LS_LAYERSOURCE &
          ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          N, UM, O1, IB, UT, K, KO1, K0, K1, K2, Q
      INTEGER ::          UM_START, UM_END
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, H1, H2
      DOUBLE PRECISION :: NUXR, PUXR, NUXR1, NUXR2, PUXR1, PUXR2

!  Local indices

      N   = GIVEN_LAYER
      KO1 = K_REAL(N) + 1
      UT  = OFFGRID_INDEX
      IB  = IBEAM

      IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
        UM_START = LOCAL_UM_START
        UM_END   = N_USER_STREAMS
      ELSE
        UM_START = IB
        UM_END   = IB
      ENDIF

!  Thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) THEN
        DO Q = 1, N_SURFACE_WFS
          DO O1 = 1, NSTOKES
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              LS_LAYERSOURCE(Q,UM,O1) = ZERO
            ENDDO
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Partial layer source function ( Homogeneous/constants variation )
!  =================================================================

!  Loop over user angles, weighting functions and STokes

      DO UM = LOCAL_UM_START, N_USER_STREAMS
       DO Q = 1, N_SURFACE_WFS
        DO O1 = 1, NSTOKES

!  Real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, K_REAL(N)
            NUXR = NCON_ALB(Q,K,N) * UHOM_UPDN(UM,O1,K,N)
            PUXR = PCON_ALB(Q,K,N) * UHOM_UPUP(UM,O1,K,N)
            H1 =  NUXR * UT_HMULT_UD(K,UM,UT)
            H2 =  PUXR * UT_HMULT_UU(K,UM,UT)
            SHOM_R = SHOM_R + H1 + H2
          ENDDO

!  Complex homogeneous solutions
!     Rob, 12/21/10, Changed N -> UT in UT_HMULT_DD, UT_HMULT_DU

          SHOM_CR = ZERO
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            NUXR1 =   NCON_ALB(Q,K1,N) * UHOM_UPDN(UM,O1,K1,N) &
                    - NCON_ALB(Q,K2,N) * UHOM_UPDN(UM,O1,K2,N)
            NUXR2 =   NCON_ALB(Q,K1,N) * UHOM_UPDN(UM,O1,K2,N) &
                    + NCON_ALB(Q,K2,N) * UHOM_UPDN(UM,O1,K1,N)
            PUXR1 =   PCON_ALB(Q,K1,N) * UHOM_UPUP(UM,O1,K1,N) &
                    - PCON_ALB(Q,K2,N) * UHOM_UPUP(UM,O1,K2,N)
            PUXR2 =   PCON_ALB(Q,K1,N) * UHOM_UPUP(UM,O1,K2,N) &
                    + PCON_ALB(Q,K2,N) * UHOM_UPUP(UM,O1,K1,N)
            H1 =   NUXR1 * UT_HMULT_UD(K1,UM,UT) &
                 - NUXR2 * UT_HMULT_UD(K2,UM,UT)
            H2 =   PUXR1 * UT_HMULT_UU(K1,UM,UT) &
                 - PUXR2 * UT_HMULT_UU(K2,UM,UT)
            SHOM_CR = SHOM_CR + H1 + H2
          ENDDO

!  homogeneous contribution

          LS_LAYERSOURCE(Q,UM,O1) = SHOM_R + SHOM_CR

!  End loops over Q, O1 and UM

        ENDDO
       ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LS_PARTLAYER_STERM_UP

!

      SUBROUTINE LS_PARTLAYER_STERM_DN ( &
        DO_OBSERVATION_GEOMETRY, &
        N_SURFACE_WFS, IBEAM, &
        OFFGRID_INDEX, GIVEN_LAYER, &
        NSTOKES, DO_THERMAL_TRANSONLY, &
        N_USER_STREAMS, LOCAL_UM_START, &
        K_REAL, K_COMPLEX, &
        UHOM_DNDN, UHOM_DNUP, &
        UT_HMULT_DU, UT_HMULT_DD, &
        NCON_ALB, PCON_ALB, &
        LS_LAYERSOURCE )

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY
      INTEGER, INTENT (IN) ::          N_SURFACE_WFS
      INTEGER, INTENT (IN) ::          GIVEN_LAYER, IBEAM
      INTEGER, INTENT (IN) ::          OFFGRID_INDEX
      INTEGER, INTENT (IN) ::          NSTOKES
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNDN &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNUP &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_DU &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_DD &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: NCON_ALB &
          ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_ALB &
          ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )

      DOUBLE PRECISION, INTENT (OUT) :: LS_LAYERSOURCE &
          ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          N, UM, O1, IB, UT, K, KO1, K0, K1, K2, Q
      INTEGER ::          UM_START, UM_END
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, H1, H2
      DOUBLE PRECISION :: NUXR, PUXR, NUXR1, NUXR2, PUXR1, PUXR2

!  Local indices

      N   = GIVEN_LAYER
      KO1 = K_REAL(N) + 1
      UT  = OFFGRID_INDEX
      IB  = IBEAM

      IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
        UM_START = LOCAL_UM_START
        UM_END   = N_USER_STREAMS
      ELSE
        UM_START = IB
        UM_END   = IB
      ENDIF

!  Thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) THEN
        DO Q = 1, N_SURFACE_WFS
          DO O1 = 1, NSTOKES
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              LS_LAYERSOURCE(Q,UM,O1) = ZERO
            ENDDO
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Partial layer source function ( Homogeneous/constants variation )
!  =================================================================

!  Loop over user angles, weighting functions and STokes

      DO UM = LOCAL_UM_START, N_USER_STREAMS
       DO Q = 1, N_SURFACE_WFS
        DO O1 = 1, NSTOKES

!  Real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, K_REAL(N)
            NUXR = NCON_ALB(Q,K,N) * UHOM_DNDN(UM,O1,K,N)
            PUXR = PCON_ALB(Q,K,N) * UHOM_DNUP(UM,O1,K,N)
            H1 =  NUXR * UT_HMULT_DD(K,UM,UT)
            H2 =  PUXR * UT_HMULT_DU(K,UM,UT)
            SHOM_R = SHOM_R + H1 + H2
          ENDDO

!  Complex homogeneous solutions
!     Rob, 12/21/10, Changed N -> UT in UT_HMULT_DD, UT_HMULT_DU

          SHOM_CR = ZERO
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            NUXR1 =   NCON_ALB(Q,K1,N) * UHOM_DNDN(UM,O1,K1,N) &
                    - NCON_ALB(Q,K2,N) * UHOM_DNDN(UM,O1,K2,N)
            NUXR2 =   NCON_ALB(Q,K1,N) * UHOM_DNDN(UM,O1,K2,N) &
                    + NCON_ALB(Q,K2,N) * UHOM_DNDN(UM,O1,K1,N)
            PUXR1 =   PCON_ALB(Q,K1,N) * UHOM_DNUP(UM,O1,K1,N) &
                    - PCON_ALB(Q,K2,N) * UHOM_DNUP(UM,O1,K2,N)
            PUXR2 =   PCON_ALB(Q,K1,N) * UHOM_DNUP(UM,O1,K2,N) &
                    + PCON_ALB(Q,K2,N) * UHOM_DNUP(UM,O1,K1,N)
            H1 =   NUXR1 * UT_HMULT_DD(K1,UM,UT) &
                 - NUXR2 * UT_HMULT_DD(K2,UM,UT)
            H2 =   PUXR1 * UT_HMULT_DU(K1,UM,UT) &
                 - PUXR2 * UT_HMULT_DU(K2,UM,UT)
            SHOM_CR = SHOM_CR + H1 + H2
          ENDDO

!  homogeneous contribution

          LS_LAYERSOURCE(Q,UM,O1) = SHOM_R + SHOM_CR

!  End loops over Q, O1 and UM

        ENDDO
       ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LS_PARTLAYER_STERM_DN

!

      SUBROUTINE VLIDORT_LS_INTEGRATED_OUTPUT ( &
        DO_INCLUDE_MVOUTPUT, FLUX_MULTIPLIER, &
        IBEAM, N_SURFACE_WFS, LS_BOA_THTONLY_SOURCE, &
        NSTOKES, NSTREAMS, N_USER_LEVELS, &
        QUAD_WEIGHTS, QUAD_STRMWTS, &
        N_DIRECTIONS, WHICH_DIRECTIONS, &
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
        UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, &
        PARTLAYERS_LAYERIDX, &
        NLAYERS, DO_THERMAL_TRANSONLY, &
        T_DELT_DISORDS, T_DISORDS_UTUP, &
        T_UTUP_EIGEN, T_UTDN_EIGEN, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        NCON_ALB, PCON_ALB, &
        T_DELT_EIGEN, &
        MINT_SURFACEWF, FLUX_SURFACEWF )

!  Quadrature output at offgrid or ongrid optical depths
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_MVOUTPUT
      INTEGER, INTENT (IN) ::          IBEAM
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER
      INTEGER, INTENT (IN) ::          N_SURFACE_WFS
      DOUBLE PRECISION, INTENT (IN) :: LS_BOA_THTONLY_SOURCE &
          ( MAX_SURFACEWFS, MAXSTREAMS, MAXSTOKES )
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      DOUBLE PRECISION, INTENT (IN) :: QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::          N_DIRECTIONS
      INTEGER, INTENT (IN) ::          WHICH_DIRECTIONS ( MAX_DIRECTIONS )
      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      INTEGER, INTENT (IN) ::          NLAYERS
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTUP &
          ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN &
          ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN &
          ( MAXEVALUES, MAX_PARTLAYERS )
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: NCON_ALB &
          ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_ALB &
          ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (INOUT) :: MINT_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) :: FLUX_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  local variables
!  ---------------

!  Local quadrature output (for debug)

      DOUBLE PRECISION :: QSURFACEWF_F &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  Help variables

      INTEGER ::          I, IDIR, WDIR, UTA, UT, N, NLEVEL, O1, Q
      DOUBLE PRECISION :: SMI, SFX

!  direction loop

      DO IDIR = 1, N_DIRECTIONS
        WDIR = WHICH_DIRECTIONS(IDIR)

!  Upwelling Jacobian output at Quadrature angles

        IF ( WDIR .EQ. UPIDX ) THEN
          DO UTA = 1, N_USER_LEVELS
            NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
              UT = PARTLAYERS_OUTINDEX(UTA)
              N  = PARTLAYERS_LAYERIDX(UT)

              CALL QUADSURFACEWF_OFFGRID_UP ( &
                UTA, UT, N, N_SURFACE_WFS, &
                FLUX_MULTIPLIER, LS_BOA_THTONLY_SOURCE, &
                NSTOKES, NSTREAMS, NLAYERS, &
                DO_THERMAL_TRANSONLY, &
                T_DELT_DISORDS, T_DISORDS_UTUP, &
                T_UTUP_EIGEN, T_UTDN_EIGEN, &
                K_REAL, K_COMPLEX, &
                SOLA_XPOS, SOLB_XNEG, &
                NCON_ALB, PCON_ALB, &
                QSURFACEWF_F )

            ELSE
              CALL QUADSURFACEWF_LEVEL_UP ( &
                UTA, NLEVEL, N_SURFACE_WFS, &
                FLUX_MULTIPLIER, LS_BOA_THTONLY_SOURCE, &
                NSTOKES, NSTREAMS, NLAYERS, &
                DO_THERMAL_TRANSONLY, &
                T_DELT_DISORDS, T_DELT_EIGEN, &
                K_REAL, K_COMPLEX, &
                SOLA_XPOS, SOLB_XNEG, &
                NCON_ALB, PCON_ALB, &
                QSURFACEWF_F )

            ENDIF
          ENDDO
        ENDIF

!  Downwelling Jacobian output at Quadrature angles

        IF ( WDIR .EQ. DNIDX ) THEN
          DO UTA = 1, N_USER_LEVELS
            NLEVEL = UTAU_LEVEL_MASK_DN(UTA)
            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
              UT = PARTLAYERS_OUTINDEX(UTA)
              N  = PARTLAYERS_LAYERIDX(UT)

              CALL QUADSURFACEWF_OFFGRID_DN ( &
                UTA, UT, N, N_SURFACE_WFS, &
                FLUX_MULTIPLIER, &
                NSTOKES, NSTREAMS, &
                DO_THERMAL_TRANSONLY, &
                T_UTUP_EIGEN, T_UTDN_EIGEN, &
                K_REAL, K_COMPLEX, &
                SOLA_XPOS, SOLB_XNEG, &
                NCON_ALB, PCON_ALB, &
                QSURFACEWF_F )

            ELSE
              CALL QUADSURFACEWF_LEVEL_DN ( &
                UTA, NLEVEL, N_SURFACE_WFS, &
                FLUX_MULTIPLIER, &
                NSTOKES, NSTREAMS, &
                DO_THERMAL_TRANSONLY, &
                T_DELT_EIGEN, &
                K_REAL, K_COMPLEX, &
                SOLA_XPOS, SOLB_XNEG, &
                NCON_ALB, PCON_ALB, &
                QSURFACEWF_F )

            ENDIF
          ENDDO
        ENDIF

!  Mean Intensity and Flux  output (diffuse term only)

        IF ( DO_INCLUDE_MVOUTPUT ) THEN
          DO UTA = 1, N_USER_LEVELS
           DO Q = 1, N_SURFACE_WFS
            DO O1 = 1, NSTOKES
              SMI = ZERO
              SFX = ZERO
              DO I = 1, NSTREAMS
                SMI = SMI + QUAD_WEIGHTS(I) * QSURFACEWF_F(Q,UTA,I,O1)
                SFX = SFX + QUAD_STRMWTS(I) * QSURFACEWF_F(Q,UTA,I,O1)
              ENDDO
              MINT_SURFACEWF(Q,UTA,IBEAM,O1,WDIR) = SMI * HALF
              FLUX_SURFACEWF(Q,UTA,IBEAM,O1,WDIR) = SFX * PI2
            ENDDO
           ENDDO
          ENDDO
        ENDIF

!  End directions loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_LS_INTEGRATED_OUTPUT

!

      SUBROUTINE QUADSURFACEWF_LEVEL_UP ( &
        UTA, NL, N_SURFACE_WFS, &
        FLUX_MULTIPLIER, LS_BOA_THTONLY_SOURCE, &
        NSTOKES, NSTREAMS, NLAYERS, &
        DO_THERMAL_TRANSONLY, &
        T_DELT_DISORDS, T_DELT_EIGEN, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        NCON_ALB, PCON_ALB, &
        QSURFACEWF_F )

!  Upwelling weighting function Fourier components at level boundary NL
!  Quadrature angles only

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          UTA, NL
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER
      INTEGER, INTENT (IN) ::          N_SURFACE_WFS
      DOUBLE PRECISION, INTENT (IN) :: LS_BOA_THTONLY_SOURCE &
          ( MAX_SURFACEWFS, MAXSTREAMS, MAXSTOKES)

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: NCON_ALB &
          ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_ALB &
          ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )

      DOUBLE PRECISION, INTENT (INOUT) :: QSURFACEWF_F &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          N, I, I1, O1, K, KO1, K0, K1, K2, LAY, Q
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, H1, H2, THELP
      DOUBLE PRECISION :: NXR, NXR1, NXR2, PXR, PXR1, PXR2

!  This depends on the level mask - if this is 0 to NLAYERS - 1, then we
!  looking at the perturbation field at the top of these layers. The
!  case where the level mask = NLAYERS is the upwelling perturbed fields
!  at the bottom of the atmosphere (treated separately).

      N = NL + 1

!  For the lowest level
!  ====================

!  Thermal transmittance-only solution

      IF ( NL.EQ.NLAYERS .and. DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO Q = 1, N_SURFACE_WFS
          DO I = 1, NSTREAMS
            THELP = LS_BOA_THTONLY_SOURCE(Q,I,O1)
            QSURFACEWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER * THELP
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Scattering solutions
!  --------------------

      IF ( NL .EQ. NLAYERS ) THEN

!  Offset

       KO1 = K_REAL(NL) + 1

!  parameter loop

       DO Q = 1, N_SURFACE_WFS

!  Stokes and streams loops

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES

!  Real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(NL)
              NXR = NCON_ALB(Q,K,NL) * SOLA_XPOS(I1,O1,K,NL)
              PXR = PCON_ALB(Q,K,NL) * SOLB_XNEG(I1,O1,K,NL)
              SHOM_R = SHOM_R + NXR*T_DELT_EIGEN(K,NL) + PXR
            ENDDO

!  Complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(NL)
             K0 = 2 * K - 2
             K1 = KO1 + K0
             K2 = K1  + 1
             NXR1 =   NCON_ALB(Q,K1,NL) * SOLA_XPOS(I1,O1,K1,NL) &
                    - NCON_ALB(Q,K2,NL) * SOLA_XPOS(I1,O1,K2,NL)
             NXR2 =   NCON_ALB(Q,K1,NL) * SOLA_XPOS(I1,O1,K2,NL) &
                    + NCON_ALB(Q,K2,NL) * SOLA_XPOS(I1,O1,K1,NL)
             PXR1 =   PCON_ALB(Q,K1,NL) * SOLB_XNEG(I1,O1,K1,NL) &
                    - PCON_ALB(Q,K2,NL) * SOLB_XNEG(I1,O1,K2,NL)
             H1 =   NXR1 * T_DELT_EIGEN(K1,NL) &
                  - NXR2 * T_DELT_EIGEN(K2,NL)
             H2 =  PXR1
             SHOM_CR = SHOM_CR + H1 + H2
            ENDDO

!  collect solution

            QSURFACEWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER*(SHOM_R+SHOM_CR)

!  Finish loops

          ENDDO
        ENDDO
       ENDDO

      ENDIF

!  For other levels in the atmosphere
!  ==================================

!  Thermal transmittance-only solution

      IF ( NL.NE.NLAYERS .and. DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO Q = 1, N_SURFACE_WFS
          DO I = 1, NSTREAMS
            THELP = LS_BOA_THTONLY_SOURCE(Q,I,O1)
            DO LAY = NLAYERS, N, -1
              THELP = THELP*T_DELT_DISORDS(I,LAY)
            ENDDO
            QSURFACEWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER * THELP
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Scattering solutions
!  --------------------

      IF ( NL .NE. NLAYERS ) THEN

!  Offset

       KO1 = K_REAL(N) + 1

!  parameter loop

       DO Q = 1, N_SURFACE_WFS

!  Stokes and streams loops

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES

!  Real homogeneous solutions

             SHOM_R = ZERO
            DO K = 1, K_REAL(N)
              NXR = NCON_ALB(Q,K,N) * SOLA_XPOS(I1,O1,K,N)
              PXR = PCON_ALB(Q,K,N) * SOLB_XNEG(I1,O1,K,N)
              SHOM_R = SHOM_R + PXR*T_DELT_EIGEN(K,N) + NXR
            ENDDO

!  Complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(N)
             K0 = 2 * K - 2
             K1 = KO1 + K0
             K2 = K1  + 1
             NXR1 =   NCON_ALB(Q,K1,N) * SOLA_XPOS(I1,O1,K1,N) &
                    - NCON_ALB(Q,K2,N) * SOLA_XPOS(I1,O1,K2,N)
             PXR1 =   PCON_ALB(Q,K1,N) * SOLB_XNEG(I1,O1,K1,N) &
                    - PCON_ALB(Q,K2,N) * SOLB_XNEG(I1,O1,K2,N)
             PXR2 =   PCON_ALB(Q,K1,N) * SOLB_XNEG(I1,O1,K2,N) &
                    + PCON_ALB(Q,K2,N) * SOLB_XNEG(I1,O1,K1,N)
             H1 =   PXR1 * T_DELT_EIGEN(K1,N) &
                  - PXR2 * T_DELT_EIGEN(K2,N)
             H2 =  NXR1
             SHOM_CR = SHOM_CR + H1 + H2
            ENDDO

!  collect solution

           QSURFACEWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER*(SHOM_R+SHOM_CR)

!  Finish loops

          ENDDO
        ENDDO
       ENDDO

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE QUADSURFACEWF_LEVEL_UP

!

      SUBROUTINE QUADSURFACEWF_LEVEL_DN ( &
        UTA, NL, N_SURFACE_WFS, &
        FLUX_MULTIPLIER, &
        NSTOKES, NSTREAMS, &
        DO_THERMAL_TRANSONLY, &
        T_DELT_EIGEN, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        NCON_ALB, PCON_ALB, &
        QSURFACEWF_F )

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          UTA, NL
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER
      INTEGER, INTENT (IN) ::          N_SURFACE_WFS
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: NCON_ALB &
          ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_ALB &
          ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )

      DOUBLE PRECISION, INTENT (INOUT) :: QSURFACEWF_F &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          N, I, O1, K, KO1, K0, K1, K2, Q
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, H1, H2
      DOUBLE PRECISION :: NXR, NXR1, NXR2, PXR, PXR1

      N = NL

!  Zero level = TOA no solutions
!  =============================

!  Downwelling weighting function at TOA ( or N = 0 ) is zero

      IF ( NL .EQ. 0 ) THEN
        DO Q = 1, N_SURFACE_WFS
          DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES
              QSURFACEWF_F(Q,UTA,I,O1) = ZERO
            ENDDO
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  For other levels in the atmosphere
!  ==================================

      IF ( NL .NE. 0 .and. DO_THERMAL_TRANSONLY) THEN
        DO Q = 1, N_SURFACE_WFS
          DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES
              QSURFACEWF_F(Q,UTA,I,O1) = ZERO
            ENDDO
          ENDDO
        ENDDO
        RETURN
      endif

!  Scattering solutions
!  --------------------

      IF ( NL .NE. 0 ) THEN

!  Offset

       KO1 = K_REAL(N) + 1

!  parameter loop

       DO Q = 1, N_SURFACE_WFS

!  Stokes and streams loops

        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES

!  Real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(N)
              NXR = NCON_ALB(Q,K,N) * SOLA_XPOS(I,O1,K,N)
              PXR = PCON_ALB(Q,K,N) * SOLB_XNEG(I,O1,K,N)
              SHOM_R = SHOM_R + NXR*T_DELT_EIGEN(K,N) + PXR
            ENDDO

!  Complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(N)
             K0 = 2 * K - 2
             K1 = KO1 + K0
             K2 = K1  + 1
             NXR1 =   NCON_ALB(Q,K1,N) * SOLA_XPOS(I,O1,K1,N) &
                    - NCON_ALB(Q,K2,N) * SOLA_XPOS(I,O1,K2,N)
             NXR2 =   NCON_ALB(Q,K1,N) * SOLA_XPOS(I,O1,K2,N) &
                    + NCON_ALB(Q,K2,N) * SOLA_XPOS(I,O1,K1,N)
             PXR1 =   PCON_ALB(Q,K1,N) * SOLB_XNEG(I,O1,K1,N) &
                    - PCON_ALB(Q,K2,N) * SOLB_XNEG(I,O1,K2,N)
             H1 =   NXR1 * T_DELT_EIGEN(K1,N) &
                  - NXR2 * T_DELT_EIGEN(K2,N)
             H2 =  PXR1
             SHOM_CR = SHOM_CR + H1 + H2
            ENDDO

!  collect solution

            QSURFACEWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER*(SHOM_R+SHOM_CR)

!  Finish loops

          ENDDO
        ENDDO
       ENDDO

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE QUADSURFACEWF_LEVEL_DN

!

      SUBROUTINE QUADSURFACEWF_OFFGRID_UP ( &
        UTA, UT, N, N_SURFACE_WFS, &
        FLUX_MULTIPLIER, LS_BOA_THTONLY_SOURCE, &
        NSTOKES, NSTREAMS, NLAYERS, &
        DO_THERMAL_TRANSONLY, &
        T_DELT_DISORDS, T_DISORDS_UTUP, &
        T_UTUP_EIGEN, T_UTDN_EIGEN, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        NCON_ALB, PCON_ALB, &
        QSURFACEWF_F )

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          UTA, UT, N
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER
      INTEGER, INTENT (IN) ::          N_SURFACE_WFS
      DOUBLE PRECISION, INTENT (IN) :: LS_BOA_THTONLY_SOURCE &
          ( MAX_SURFACEWFS, MAXSTREAMS, MAXSTOKES )
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTUP &
          ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN &
          ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN &
          ( MAXEVALUES, MAX_PARTLAYERS )
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: NCON_ALB &
          ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_ALB &
          ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )

      DOUBLE PRECISION, INTENT (INOUT) :: QSURFACEWF_F &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          I, I1, O1, K, KO1, K0, K1, K2, LAY, Q
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, H1, H2, THELP
      DOUBLE PRECISION :: NXR, NXR1, NXR2, PXR, PXR1, PXR2

!  For thermal transmittance-only
!  ------------------------------

      IF ( DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO Q = 1, N_SURFACE_WFS
          DO I = 1, NSTREAMS
            THELP = LS_BOA_THTONLY_SOURCE(Q,I,O1)
            DO LAY = NLAYERS, N+1, -1
              THELP = THELP*T_DELT_DISORDS(I,LAY)
            ENDDO
            THELP = THELP*T_DISORDS_UTUP(I,UT)
            QSURFACEWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER * THELP
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  For those optical depths at off-grid levels
!  -------------------------------------------

!  Offset

      KO1 = K_REAL(N) + 1

!  parameter loop

      DO Q = 1, N_SURFACE_WFS

!  stream directions

       DO I = 1, NSTREAMS
        I1 = I + NSTREAMS
        DO O1 = 1, NSTOKES

!  real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, K_REAL(N)
            NXR = NCON_ALB(Q,K,N) * SOLA_XPOS(I1,O1,K,N)
            PXR = PCON_ALB(Q,K,N) * SOLB_XNEG(I1,O1,K,N)
            SHOM_R = SHOM_R + NXR * T_UTDN_EIGEN(K,UT) &
                            + PXR * T_UTUP_EIGEN(K,UT)
!mick fix
            !SHOM_R = SHOM_R + H1 + H2
          ENDDO

!  complex homogeneous solutions
!     Rob, 12/21/10, Changed N -> UT in T_UTDN_EIGEN, T_UTUP_EIGEN

          SHOM_CR = ZERO
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            NXR1 =   NCON_ALB(Q,K1,N) * SOLA_XPOS(I1,O1,K1,N) &
                   - NCON_ALB(Q,K2,N) * SOLA_XPOS(I1,O1,K2,N)
            NXR2 =   NCON_ALB(Q,K1,N) * SOLA_XPOS(I1,O1,K2,N) &
                   + NCON_ALB(Q,K2,N) * SOLA_XPOS(I1,O1,K1,N)
            PXR1 =   PCON_ALB(Q,K1,N) * SOLB_XNEG(I1,O1,K1,N) &
                   - PCON_ALB(Q,K2,N) * SOLB_XNEG(I1,O1,K2,N)
            PXR2 =   PCON_ALB(Q,K1,N) * SOLB_XNEG(I1,O1,K2,N) &
                   + PCON_ALB(Q,K2,N) * SOLB_XNEG(I1,O1,K1,N)
            H1 =   NXR1 * T_UTDN_EIGEN(K1,UT) &
                 - NXR2 * T_UTDN_EIGEN(K2,UT)
            H2 =   PXR1 * T_UTUP_EIGEN(K1,UT) &
                 - PXR2 * T_UTUP_EIGEN(K2,UT)
            SHOM_CR = SHOM_CR + H1 + H2
          ENDDO

!  set result

          QSURFACEWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER*(SHOM_R+SHOM_CR)

!  end O1 and I loops, end Q loop

        ENDDO
       ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE QUADSURFACEWF_OFFGRID_UP

!

      SUBROUTINE QUADSURFACEWF_OFFGRID_DN ( &
        UTA, UT, N, N_SURFACE_WFS, &
        FLUX_MULTIPLIER, &
        NSTOKES, NSTREAMS, &
        DO_THERMAL_TRANSONLY, &
        T_UTUP_EIGEN, T_UTDN_EIGEN, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        NCON_ALB, PCON_ALB, &
        QSURFACEWF_F )

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          UTA, UT, N
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER
      INTEGER, INTENT (IN) ::          N_SURFACE_WFS
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN &
          ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN &
          ( MAXEVALUES, MAX_PARTLAYERS )
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: NCON_ALB &
          ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_ALB &
          ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )

      DOUBLE PRECISION, INTENT (INOUT) :: QSURFACEWF_F &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          I, O1, K, KO1, K0, K1, K2, Q
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, H1, H2
      DOUBLE PRECISION :: NXR, NXR1, NXR2, PXR, PXR1, PXR2

!  For thermal transmittance-only, no contribution
!  -----------------------------------------------

      IF ( DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO Q = 1, N_SURFACE_WFS
          DO I = 1, NSTREAMS
            QSURFACEWF_F(Q,UTA,I,O1) = ZERO
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  For those optical depths at off-grid levels
!  -------------------------------------------

!  Offset

      KO1 = K_REAL(N) + 1

!  parameter loop

      DO Q = 1, N_SURFACE_WFS

!  stream directions

       DO I = 1, NSTREAMS
        DO O1 = 1, NSTOKES

!  real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, K_REAL(N)
            NXR = NCON_ALB(Q,K,N) * SOLA_XPOS(I,O1,K,N)
            PXR = PCON_ALB(Q,K,N) * SOLB_XNEG(I,O1,K,N)
            SHOM_R = SHOM_R + NXR * T_UTDN_EIGEN(K,UT) &
                            + PXR * T_UTUP_EIGEN(K,UT)
!mick fix
            !SHOM_R = SHOM_R + H1 + H2
          ENDDO

!  complex homogeneous solutions
!     Rob, 12/21/10, Changed N -> UT in T_UTDN_EIGEN, T_UTUP_EIGEN

          SHOM_CR = ZERO
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            NXR1 =   NCON_ALB(Q,K1,N) * SOLA_XPOS(I,O1,K1,N) &
                   - NCON_ALB(Q,K2,N) * SOLA_XPOS(I,O1,K2,N)
            NXR2 =   NCON_ALB(Q,K1,N) * SOLA_XPOS(I,O1,K2,N) &
                   + NCON_ALB(Q,K2,N) * SOLA_XPOS(I,O1,K1,N)
            PXR1 =   PCON_ALB(Q,K1,N) * SOLB_XNEG(I,O1,K1,N) &
                   - PCON_ALB(Q,K2,N) * SOLB_XNEG(I,O1,K2,N)
            PXR2 =   PCON_ALB(Q,K1,N) * SOLB_XNEG(I,O1,K2,N) &
                   + PCON_ALB(Q,K2,N) * SOLB_XNEG(I,O1,K1,N)
            H1 =   NXR1 * T_UTDN_EIGEN(K1,UT) &
                 - NXR2 * T_UTDN_EIGEN(K2,UT)
            H2 =   PXR1 * T_UTUP_EIGEN(K1,UT) &
                 - PXR2 * T_UTUP_EIGEN(K2,UT)
            SHOM_CR = SHOM_CR + H1 + H2
          ENDDO

!  set result

          QSURFACEWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER*(SHOM_R+SHOM_CR)

!  end O1 and I loops, end Q loop

        ENDDO
       ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE QUADSURFACEWF_OFFGRID_DN

      END MODULE vlidort_ls_wfsurface

