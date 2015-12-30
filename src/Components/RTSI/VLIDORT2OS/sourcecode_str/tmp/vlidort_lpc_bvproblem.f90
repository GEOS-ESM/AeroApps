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
! #            L_BVP_BACKSUB                                    #
! #            L_BVP_SURFACE_SETUP                              #
! #                                                             #
! ###############################################################


      MODULE vlidort_lpc_bvproblem

      PRIVATE
      PUBLIC :: L_BVP_BACKSUB, &
                L_BVP_SURFACE_SETUP

      CONTAINS

      SUBROUTINE L_BVP_BACKSUB ( &
        LAYER_TO_VARY, N_LAYER_WFS, &
        NLAYERS, NTOTAL, N_SUBDIAG, &
        N_SUPDIAG, NSTKS_NSTRMS, &
        NSTKS_NSTRMS_2, &
        K_REAL, K_COMPLEX, &
        BANDMAT2, IPIVOT, &
        SMAT2, SIPIVOT, &
        LCON, MCON, &
        COL2_WF, SCOL2_WF, &
        NCON, PCON, &
        STATUS, MESSAGE, TRACE )

!  Solves the linearized boundary value problem.

      USE VLIDORT_PARS
      USE LAPACK_TOOLS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          LAYER_TO_VARY
      INTEGER, INTENT (IN) ::          N_LAYER_WFS
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
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )

      DOUBLE PRECISION, INTENT (INOUT) :: COL2_WF ( MAXTOTAL, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: SCOL2_WF &
          ( MAXSTRMSTKS_2, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (OUT) ::  NCON &
          ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) ::  PCON &
          ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )

      INTEGER, INTENT (OUT) ::             STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  local variables
!  ---------------

      INTEGER ::           C0, LAY, K, Q, KO1, K0, K1, K2
      INTEGER ::           IROW, IROW1, IROW_S, IROW1_S
      INTEGER ::           INFO
      CHARACTER (LEN=3) :: CI, CN

!  Intialize status

      STATUS = VLIDORT_SUCCESS

!  BVP back-substitution: With compression (multilayers)
!  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  LAPACK substitution (DGBTRS) using RHS column vector COL2_WF
!  BV solution for perturbed integration constants
!    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS &
           ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, N_LAYER_WFS, &
              BANDMAT2, MAXBANDTOTAL, IPIVOT, &
              COL2_WF, MAXTOTAL, INFO )

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(CN, '(I3)' ) LAYER_TO_VARY
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'Atmos_Wfs for layer '//CN// &
                        '; DGBTRS call in L_BVP_BACKSUB'
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  Set Linearized integration constants NCON and PCON, all layers

        DO LAY = 1, NLAYERS
          C0 = (LAY-1)*NSTKS_NSTRMS_2
          KO1 = K_REAL(LAY) + 1
          DO K = 1, K_REAL(LAY)
            IROW = K
            IROW1 = IROW + NSTKS_NSTRMS
            DO Q = 1, N_LAYER_WFS
              NCON(K,LAY,Q) = COL2_WF(C0+IROW,Q)
              PCON(K,LAY,Q) = COL2_WF(C0+IROW1,Q)
            ENDDO
          ENDDO
          DO K = 1, K_COMPLEX(LAY)
            K0 = 2*K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            IROW    = K + K_REAL(LAY)
            IROW1   = IROW + NSTKS_NSTRMS
            IROW_S  = IROW + K_COMPLEX(LAY)
            IROW1_S = IROW_S + NSTKS_NSTRMS
            DO Q = 1, N_LAYER_WFS
              NCON(K1,LAY,Q) = COL2_WF(C0+IROW,   Q)
              NCON(K2,LAY,Q) = COL2_WF(C0+IROW_S, Q)
              PCON(K1,LAY,Q) = COL2_WF(C0+IROW1,  Q)
              PCON(K2,LAY,Q) = COL2_WF(C0+IROW1_S,Q)
            ENDDO
          ENDDO
        ENDDO

!  Solve the boundary problem: No compression, Single Layer only
!  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2

        CALL DGETRS &
           ( 'N', NTOTAL, N_LAYER_WFS, SMAT2, MAXSTRMSTKS_2, SIPIVOT, &
              SCOL2_WF, MAXSTRMSTKS_2, INFO )

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE = 'Atmos_Wfs for 1-layer: DGETRS call in L_BVP_BACKSUB'
          STATUS = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  Set linearized integration constants NCON and PCON, 1 layer

        LAY = 1
        KO1 = K_REAL(LAY) + 1
        DO K = 1, K_REAL(LAY)
          IROW = K
          IROW1 = IROW + NSTKS_NSTRMS
          DO Q = 1, N_LAYER_WFS
            NCON(K,LAY,Q) = SCOL2_WF(IROW, Q)
            PCON(K,LAY,Q) = SCOL2_WF(IROW1,Q)
          ENDDO
        ENDDO
        DO K = 1, K_COMPLEX(LAY)
          K0 = 2*K - 2
          K1 = KO1 + K0
          K2 = K1  + 1
          IROW    = K + K_REAL(LAY)
          IROW1   = IROW + NSTKS_NSTRMS
          IROW_S  = K + K_REAL(LAY) + K_COMPLEX(LAY)
          IROW1_S = IROW_S + NSTKS_NSTRMS
          DO Q = 1, N_LAYER_WFS
            NCON(K1,LAY,Q) = SCOL2_WF(IROW,    Q)
            NCON(K2,LAY,Q) = SCOL2_WF(IROW_S,  Q)
            PCON(K1,LAY,Q) = SCOL2_WF(IROW1,   Q)
            PCON(K2,LAY,Q) = SCOL2_WF(IROW1_S, Q)
          ENDDO
        ENDDO

      ENDIF

!  debug------------------------------------------

!       if ( fourier .EQ.3 .and.ibeam.eq.1 ) then
!        do n = 1, nlayers
!          write(76,'(a,2i3,1p2e18.10)')
!     &     'hey',n,variation_index,ncon(3,n,1),pcon(3,n,1)
!        enddo
!       endif

!        if ( do_debug_write ) then
!        if ( fourier.eq.1 ) then
!         Q = 1
!         DO LAY = 1, NLAYERS
!          KO1 = K_REAL(LAY) + 1
!          DO K = 1, K_REAL(LAY)
!           write(51,'(4i3,1p4e20.10)')FOURIER,IBEAM,K,LAY,
!     &                LCON(K,LAY),  MCON(K,LAY),
!     &                NCON(K,LAY,Q),PCON(K,LAY,Q)
!          ENDDO
!          DO K = 1, K_COMPLEX(LAY)
!           K0 = 2*K - 2
!           K1 = KO1 + K0
!           K2 = K1  + 1
!           write(51,'(4i3,1p8e20.10)')FOURIER,IBEAM,K,LAY,
!     &                LCON(K1,LAY), MCON(K1,LAY),
!     &                LCON(K2,LAY), MCON(K2,LAY),
!     &                NCON(K1,LAY,Q),PCON(K1,LAY,Q),
!     &                NCON(K2,LAY,Q),PCON(K2,LAY,Q)
!          ENDDO
!         ENDDO
!        endif
!        ENDIF

!  finish

      RETURN
      END SUBROUTINE L_BVP_BACKSUB

!

      SUBROUTINE L_BVP_SURFACE_SETUP ( &
        DO_INCLUDE_SURFACE, MODIFIED_BCL4, &
        IBEAM, FOURIER_COMPONENT, &
        SURFACE_FACTOR, N_LAYER_WFS, &
        NSTOKES, NSTREAMS, &
        NLAYERS, DO_LAMBERTIAN_SURFACE, &
        LAMBERTIAN_ALBEDO, BRDF_F, &
        QUAD_STRMWTS, NSTKS_NSTRMS, &
        MUELLER_INDEX, T_DELT_EIGEN, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, L_T_DELT_EIGEN, &
        L_SOLA_XPOS, L_SOLB_XNEG, &
        L_WLOWER, &
        R2_L_BEAM, R2_L_HOMP, R2_L_HOMM )

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::          MODIFIED_BCL4
      INTEGER, INTENT (IN) ::          IBEAM
      INTEGER, INTENT (IN) ::          FOURIER_COMPONENT
      DOUBLE PRECISION, INTENT (IN) :: SURFACE_FACTOR
      INTEGER, INTENT (IN) ::          N_LAYER_WFS
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      DOUBLE PRECISION, INTENT (IN) :: LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: BRDF_F &
          ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_EIGEN &
          ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_WLOWER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (OUT) :: R2_L_BEAM &
          ( MAXSTREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: R2_L_HOMP &
          ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: R2_L_HOMM &
          ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )

!  Local variables
!  ---------------

      DOUBLE PRECISION :: PV_W &
          ( MAXSTREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: HV_P &
          ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: HV_M &
          ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )

      INTEGER ::          I, J, O1, O2, OM, IB, Q, N, M, NELEMENTS
      INTEGER ::          K, KO1, K0, K1, K2

      DOUBLE PRECISION :: H1R, H1I, REFL_P, REFL_B, REFL_M, KMULT
      DOUBLE PRECISION :: H1, H2, HP, HM, BIREF

      DOUBLE PRECISION :: H1_CR,   H2_CR,   H1_CI,   H2_CI
      DOUBLE PRECISION :: H1_S_CR, H2_S_CR, H1_S_CI, H2_S_CI

!  Initial section
!  ---------------

!  Beam index

      IB = IBEAM

!  Fourier component

      M = FOURIER_COMPONENT

!  Always zero the result to start

      DO I = 1, NSTREAMS
       DO O1 = 1, NSTOKES
        DO Q = 1, N_LAYER_WFS
          R2_L_BEAM(I,O1,Q)    = ZERO
          DO K = 1, NSTKS_NSTRMS
            R2_L_HOMP(I,O1,K,Q) = ZERO
            R2_L_HOMM(I,O1,K,Q) = ZERO
          ENDDO
        ENDDO
       ENDDO
      ENDDO

!  Return if no albedo

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  Return if Fourier component > 0 (Lambertian)

      IF ( DO_LAMBERTIAN_SURFACE .and. M.gt.0 ) RETURN

!  Set up Auxiliary arrays
!  -----------------------

!  Only require (1,1) component

      O1 = 1

!  Last layer

      N = NLAYERS

!  number of elements

      IF ( DO_LAMBERTIAN_SURFACE ) THEN
        NELEMENTS = 1
      ELSE
        NELEMENTS = NSTOKES
      ENDIF

!  Need Beam arrays regardless
!   Always something from layers above the PBL and also from PBL

      DO J = 1, NSTREAMS
        DO O1 = 1, NSTOKES
          DO Q = 1, N_LAYER_WFS
            PV_W(J,O1,Q) = L_WLOWER(J,O1,N,Q) * QUAD_STRMWTS(J)
          ENDDO
        ENDDO
      ENDDO

!   Modified boundary condition
!     This applies when PBL is varying layer.
!       Require homogeneous solutions linearizations

      IF ( MODIFIED_BCL4 ) THEN

!  start loops

        DO J = 1, NSTREAMS
         DO O1 = 1, NSTOKES
          DO Q = 1, N_LAYER_WFS

!  real homogeneous solution contributions

           DO K = 1, K_REAL(N)
            H1 = L_SOLA_XPOS(J,O1,K,N,Q) *   T_DELT_EIGEN(K,N) + &
                   SOLA_XPOS(J,O1,K,N)   * L_T_DELT_EIGEN(K,N,Q)
            H2 = L_SOLB_XNEG(J,O1,K,N,Q)
            HV_P(J,O1,K,Q) = QUAD_STRMWTS(J)*H1
            HV_M(J,O1,K,Q) = QUAD_STRMWTS(J)*H2
           ENDDO

!  Complex homogeneous solution contributions

           KO1 = K_REAL(N) + 1
           DO K = 1, K_COMPLEX(N)
            K0 = 2*K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            H1R = L_SOLA_XPOS(J,O1,K1,N,Q) *   T_DELT_EIGEN(K1,N) &
                - L_SOLA_XPOS(J,O1,K2,N,Q) *   T_DELT_EIGEN(K2,N) &
                +   SOLA_XPOS(J,O1,K1,N)   * L_T_DELT_EIGEN(K1,N,Q) &
                -   SOLA_XPOS(J,O1,K2,N)   * L_T_DELT_EIGEN(K2,N,Q)
            H1I = L_SOLA_XPOS(J,O1,K1,N,Q) *   T_DELT_EIGEN(K2,N) &
                + L_SOLA_XPOS(J,O1,K2,N,Q) *   T_DELT_EIGEN(K1,N) &
                +   SOLA_XPOS(J,O1,K1,N)   * L_T_DELT_EIGEN(K2,N,Q) &
                +   SOLA_XPOS(J,O1,K2,N)   * L_T_DELT_EIGEN(K1,N,Q)
            HV_P(J,O1,K1,Q) = QUAD_STRMWTS(J)* H1R
            HV_P(J,O1,K2,Q) = QUAD_STRMWTS(J)* H1I
            HV_M(J,O1,K1,Q) = QUAD_STRMWTS(J)* L_SOLB_XNEG(J,O1,K1,N,Q)
            HV_M(J,O1,K2,Q) = QUAD_STRMWTS(J)* L_SOLB_XNEG(J,O1,K2,N,Q)
           ENDDO

!  End loops

          ENDDO
         ENDDO
        ENDDO

!  End modified BCL4 condition

      ENDIF

!  Lambertian condition
!  ====================

!  Skip if BRDF surface

      if ( .not. DO_LAMBERTIAN_SURFACE ) go to 988

!  reflection

      KMULT = SURFACE_FACTOR * LAMBERTIAN_ALBEDO

!  only 1 component

      O1 = 1

!  Integrated Downward reflection (Calculation)
!  --------------------------------------------

      DO Q = 1, N_LAYER_WFS

!  Particular solution (only for the first Stokes component)

       REFL_B = ZERO
       DO J = 1, NSTREAMS
         REFL_B = REFL_B + PV_W(J,O1,Q)
       ENDDO
       REFL_B = REFL_B * KMULT
       DO I = 1, NSTREAMS
         R2_L_BEAM(I,O1,Q) = REFL_B
       ENDDO

!  Homogeneous solutions for the modified condition

       IF ( MODIFIED_BCL4 ) THEN

!  Homogeneous real solutions

         DO K = 1, K_REAL(NLAYERS)
          REFL_P = ZERO
          REFL_M = ZERO
          DO J = 1, NSTREAMS
           REFL_P = REFL_P + HV_P(J,O1,K,Q)
           REFL_M = REFL_M + HV_M(J,O1,K,Q)
          ENDDO
          REFL_P = REFL_P * KMULT
          REFL_M = REFL_M * KMULT
          DO I = 1, NSTREAMS
           R2_L_HOMP(I,O1,K,Q) = REFL_P
           R2_L_HOMM(I,O1,K,Q) = REFL_M
          ENDDO
         ENDDO

!  Homogeneous complex solutions

         KO1 = K_REAL(NLAYERS) + 1
         DO K = 1, K_COMPLEX(NLAYERS)
          K0 = 2*K - 2
          K1 = KO1 + K0
          K2 = K1  + 1
          H1_CR = ZERO
          H1_CI = ZERO
          H2_CR = ZERO
          H2_CI = ZERO
          DO J = 1, NSTREAMS
           H1_CR = H1_CR + HV_P(J,O1,K1,Q)
           H1_CI = H1_CI + HV_P(J,O1,K2,Q)
           H2_CR = H2_CR + HV_M(J,O1,K1,Q)
           H2_CI = H2_CI + HV_M(J,O1,K2,Q)
          ENDDO
          H1_CR = H1_CR * KMULT
          H1_CI = H1_CI * KMULT
          H2_CR = H2_CR * KMULT
          H2_CI = H2_CI * KMULT
          DO I = 1, NSTREAMS
            R2_L_HOMP(I,O1,K1,Q) = H1_CR
            R2_L_HOMP(I,O1,K2,Q) = H1_CI
            R2_L_HOMM(I,O1,K1,Q) = H2_CR
            R2_L_HOMM(I,O1,K2,Q) = H2_CI
          ENDDO
         ENDDO

!  End modified boundary condition clause

        ENDIF

!  end parameter loop

      ENDDO

!  Complete Lambertian case

      RETURN

!  BRDF surface condition
!  ======================

!  Continuation point

988   continue

!  start loops

      DO I = 1, NSTREAMS
        DO O1 = 1, NSTOKES
          DO Q = 1, N_LAYER_WFS

!  Particular solution

            REFL_B = ZERO
            DO J = 1, NSTREAMS
              H1 = ZERO
              DO O2 = 1, NSTOKES
               OM = MUELLER_INDEX(O1,O2)
               H1 = H1 + PV_W(J,O2,Q)*BRDF_F(M,OM,J,I)
              ENDDO
              REFL_B = REFL_B + H1
            ENDDO
            REFL_B = REFL_B * SURFACE_FACTOR
            R2_L_BEAM(I,O1,Q) = R2_L_BEAM(I,O1,Q) + REFL_B

!  Homogeneous solutions for the modified condition

            IF ( MODIFIED_BCL4 ) THEN

!  Homogeneous real solutions

              DO K = 1, K_REAL(NLAYERS)
               REFL_P = ZERO
               REFL_M = ZERO
               DO J = 1, NSTREAMS
                HP = ZERO
                HM = ZERO
                DO O2 = 1, NSTOKES
                 OM = MUELLER_INDEX(O1,O2)
                 HP = HP + HV_P(J,O2,K,Q)*BRDF_F(M,OM,J,I)
                 HM = HM + HV_M(J,O2,K,Q)*BRDF_F(M,OM,J,I)
                ENDDO
                REFL_P = REFL_P + HP
                REFL_M = REFL_M + HM
               ENDDO
               R2_L_HOMP(I,O1,K,Q) =  REFL_P * SURFACE_FACTOR
               R2_L_HOMM(I,O1,K,Q) =  REFL_M * SURFACE_FACTOR
              ENDDO

!  homogeneous complex solutions

              KO1 = K_REAL(NLAYERS) + 1
              DO K = 1, K_COMPLEX(NLAYERS)
               K0 = 2*K - 2
               K1 = KO1 + K0
               K2 = K1  + 1
               H1_CR = ZERO
               H2_CR = ZERO
               H1_CI = ZERO
               H2_CI = ZERO
               DO J = 1, NSTREAMS
                H1_S_CR = ZERO
                H2_S_CR = ZERO
                H1_S_CI = ZERO
                H2_S_CI = ZERO
                DO O2 = 1, NSTOKES
                 OM    = MUELLER_INDEX(O1,O2)
                 BIREF = BRDF_F(M,OM,J,I)
                 H1_S_CR = H1_S_CR + HV_P(J,O2,K1,Q)*BIREF
                 H2_S_CR = H2_S_CR + HV_M(J,O2,K1,Q)*BIREF
                 H1_S_CI = H1_S_CI + HV_P(J,O2,K2,Q)*BIREF
                 H2_S_CI = H2_S_CI + HV_M(J,O2,K2,Q)*BIREF
                ENDDO
                H1_CR = H1_CR + H1_S_CR
                H2_CR = H2_CR + H2_S_CR
                H1_CI = H1_CI + H1_S_CI
                H2_CI = H2_CI + H2_S_CI
               ENDDO
               R2_L_HOMP(I,O1,K1,Q) = H1_CR * SURFACE_FACTOR
               R2_L_HOMM(I,O1,K1,Q) = H2_CR * SURFACE_FACTOR
               R2_L_HOMP(I,O1,K2,Q) = H1_CI * SURFACE_FACTOR
               R2_L_HOMM(I,O1,K2,Q) = H2_CI * SURFACE_FACTOR
              ENDDO

!  End modified condition

            ENDIF

!  end parameter and stream and Stokes loops

          ENDDO
        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE L_BVP_SURFACE_SETUP

      END MODULE vlidort_lpc_bvproblem
