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
! #              VLIDORT_L_QHOM_SOLUTION                        #
! #              VLIDORT_L_UHOM_SOLUTION                        #
! #                                                             #
! #              L_HMULT_MASTER                                 #
! #                                                             #
! ###############################################################

!  Notes for Version 2.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob Fix 05/10/13  - HSINGO defined using TAYLOR_SMALL
!     Rob Fix 02/19/14  - ZETA_M = RHO_M for the degenerate case....
!                         otherwise ZETA_M = ONE / RHO_M
!                         COMPLEX HSINGO removed.

      MODULE vlidort_lpc_solutions

      USE vlidort_Taylor_m, ONLY : TAYLOR_SERIES_L_1

      PRIVATE
      PUBLIC :: VLIDORT_L_QHOM_SOLUTION, &
                VLIDORT_L_UHOM_SOLUTION, &
                L_HMULT_MASTER

      CONTAINS

      SUBROUTINE VLIDORT_L_QHOM_SOLUTION ( &
        GIVEN_LAYER, FOURIER, DO_VARY, N_PARAMETERS, &
        DO_SOLUTION_SAVING, DO_DEBUG_WRITE, &
        NSTOKES, NSTREAMS, &
        N_USER_LEVELS, &
        QUAD_STREAMS, QUAD_HALFWTS, &
        NMOMENTS, NSTREAMS_2, &
        NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
        DO_REAL_EIGENSOLVER, DO_LAYER_SCATTERING, &
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
        PARTLAYERS_LAYERIDX, DMAT, &
        DELTAU_VERT, PARTAU_VERT, &
        T_DELT_EIGEN, T_UTUP_EIGEN, &
        T_UTDN_EIGEN, &
        PI_XQP, PI_XQM_PRE, &
        SAB, DAB, &
        EIGENMAT_SAVE, REAL_KSQ, &
        IMAG_KSQ, LEFT_EVEC, &
        RITE_EVEC, EIGENDEGEN, &
        EIGENMASK_R, EIGENMASK_C, &
        K_REAL, K_COMPLEX, &
        KEIGEN, KEIGEN_CSQ, &
        FWD_SUMVEC, FWD_DIFVEC, &
        L_DELTAU_VERT, L_OMEGA_GREEK, &
        L_T_DELT_DISORDS, &
        L_T_DISORDS_UTDN, L_T_DISORDS_UTUP, &
        L_SAB, L_DAB, &
        L_EIGENMAT, L_KEIGEN, &
        L_SOLA_XPOS, L_SOLB_XNEG, &
        L_T_DELT_EIGEN, L_T_UTUP_EIGEN, &
        L_T_UTDN_EIGEN, &
        STATUS, MESSAGE, TRACE )

!  Linearization of the homogeneous solutions.

      USE VLIDORT_PARS
      USE LAPACK_TOOLS

      IMPLICIT NONE

!  Strictly in (non-linearized)

      INTEGER, INTENT (IN) ::          GIVEN_LAYER
      INTEGER, INTENT (IN) ::          FOURIER
      LOGICAL, INTENT (IN) ::          DO_VARY
      INTEGER, INTENT (IN) ::          N_PARAMETERS
      LOGICAL, INTENT (IN) ::          DO_SOLUTION_SAVING
      LOGICAL, INTENT (IN) ::          DO_DEBUG_WRITE
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_HALFWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::          NMOMENTS
      INTEGER, INTENT (IN) ::          NSTREAMS_2
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS_2
      LOGICAL, INTENT (IN) ::          DO_REAL_EIGENSOLVER &
          ( 0:MAXMOMENTS, MAXLAYERS )
      LOGICAL, INTENT (IN) ::          DO_LAYER_SCATTERING &
          ( 0:MAXMOMENTS, MAXLAYERS )
      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DMAT ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PARTAU_VERT ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN &
          ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN &
          ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PI_XQP &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_XQM_PRE &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )

!  Inout (Ostensibly in) ( non-linearized )

      DOUBLE PRECISION, INTENT (INOUT) :: SAB &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: DAB &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: EIGENMAT_SAVE &
          ( MAXEVALUES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: REAL_KSQ ( MAXSTRMSTKS )
      DOUBLE PRECISION, INTENT (INOUT) :: IMAG_KSQ ( MAXSTRMSTKS )
      DOUBLE PRECISION, INTENT (INOUT) :: LEFT_EVEC ( MAXSTRMSTKS, MAXSTRMSTKS )
      DOUBLE PRECISION, INTENT (INOUT) :: RITE_EVEC ( MAXSTRMSTKS, MAXSTRMSTKS )
      LOGICAL, INTENT (INOUT) ::          EIGENDEGEN ( MAXSTRMSTKS, MAXLAYERS )
      INTEGER, INTENT (INOUT) ::          EIGENMASK_R ( MAXEVALUES )
      INTEGER, INTENT (INOUT) ::          EIGENMASK_C ( MAXEVALUES )
      INTEGER, INTENT (INOUT) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (INOUT) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: KEIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: KEIGEN_CSQ ( MAXEVALUES )
      DOUBLE PRECISION, INTENT (INOUT) :: FWD_SUMVEC &
          ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION, INTENT (INOUT) :: FWD_DIFVEC &
          ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )

!  LInearized, Strictly in

      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_GREEK &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_DISORDS &
          ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DISORDS_UTDN &
          ( MAXSTREAMS, MAX_USER_LEVELS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DISORDS_UTUP &
          ( MAXSTREAMS, MAX_USER_LEVELS, MAX_ATMOSWFS )

!  Linearized, inout (Ostensibly out)

      DOUBLE PRECISION, INTENT (INOUT) ::  L_SAB &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, &
            MAXSTOKES,  MAXLAYERS,  MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) ::  L_DAB &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, &
            MAXSTOKES,  MAXLAYERS,  MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) ::  L_EIGENMAT &
          ( MAXEVALUES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) ::  L_KEIGEN &
          ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) ::  L_SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) ::  L_SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) ::  L_T_DELT_EIGEN &
          ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) ::  L_T_UTUP_EIGEN &
          ( MAXEVALUES, MAX_USER_LEVELS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) ::  L_T_UTDN_EIGEN &
          ( MAXEVALUES, MAX_USER_LEVELS, MAX_ATMOSWFS )

!  Exception (Strictly out)

      INTEGER, INTENT (OUT) ::           STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  Local variables
!  ---------------

!  local matrices for eigenvalue computation
!  -----------------------------------------

!  For the real calculation based on DGEEV

      DOUBLE PRECISION :: &
            HMAT1(MAXSTRMSTKS,MAXSTRMSTKS), &
            HVEC1(MAXSTRMSTKS,MAX_ATMOSWFS), &
            HVEC_AUX(MAXSTRMSTKS,MAX_ATMOSWFS)

!  For the complex calculations based on DGEEV

      DOUBLE PRECISION :: &
            HMAT2(MAXSTRMSTKS_2,MAXSTRMSTKS_2), &
            HVEC2(MAXSTRMSTKS_2,MAX_ATMOSWFS)
      DOUBLE PRECISION :: &
            HMATRED(MAXSTRMSTKS_21,MAXSTRMSTKS_21), &
            HVECRED(MAXSTRMSTKS_21,MAX_ATMOSWFS)
      DOUBLE PRECISION :: &
            HVEC_CR(MAXSTRMSTKS,MAX_ATMOSWFS), &
            HVEC_CI(MAXSTRMSTKS,MAX_ATMOSWFS)
      DOUBLE PRECISION :: &
            ADJN_CR(MAXSTRMSTKS), &
            ADJN_CI(MAXSTRMSTKS)

!  For the real calculation based on ASYMTX

      DOUBLE PRECISION :: &
            HMAT1A(MAXSTRMSTKS_P1,MAXSTRMSTKS_P1), &
            HVEC1A(MAXSTRMSTKS_P1,MAX_ATMOSWFS)

!  Local linearized vectors

      DOUBLE PRECISION :: &
            L_FWD_DIFVEC_R(MAXSTREAMS,MAXSTOKES), &
            L_FWD_SUMVEC_R(MAXSTREAMS,MAXSTOKES,MAX_ATMOSWFS)
      DOUBLE PRECISION :: &
            L_FWD_DIFVEC_CR(MAXSTREAMS,MAXSTOKES), &
            L_FWD_SUMVEC_CR(MAXSTREAMS,MAXSTOKES)
      DOUBLE PRECISION :: &
            L_FWD_DIFVEC_CI(MAXSTREAMS,MAXSTOKES), &
            L_FWD_SUMVEC_CI(MAXSTREAMS,MAXSTOKES)

!  Pivoting arrays for LAPACK linear-algebra solvers

      INTEGER ::          LAPACK_IPIV1(MAXSTRMSTKS)
      INTEGER ::          LAPACK_IPIV21(MAXSTRMSTKS_21)
      INTEGER ::          LAPACK_INFO, NSTKS_NSTRMS_21, NSTMSTKS_P1
      INTEGER ::          LAPACK_IPIV1A(MAXSTRMSTKS_P1)

!  Miscellaneous local variables

      INTEGER ::          I, J, I1, J1, JC, IR, JCOL, IROW, IROW1, ILAST
      INTEGER ::          L, N, M, Q, O1, O2, O3, O4, NA, UTA, UT
      INTEGER ::          K, KO1, K0, K1, K2
      DOUBLE PRECISION :: DP, DM, KSQD, KVAL, KVAL2, XINV, FAC
      DOUBLE PRECISION :: H1P, H1M, H2P, H2M, H1, H2, L_OG
      DOUBLE PRECISION :: TBAS, TBAS_CR, TBAS_CI, LTQ, LCR, LCI
      DOUBLE PRECISION :: TAU_UP, TAU_DN, H_DN, H_UP, VAR
      DOUBLE PRECISION :: TBASUP_CR, TBASUP_CI, TBASDN_CR, TBASDN_CI

!  additional variables

      INTEGER           :: AA, AA1, TC, TR, I_NONTRIVIAL
      DOUBLE PRECISION  :: SR, SCR, SCI, HR, HCR, HCI, MODULUS
      DOUBLE PRECISION  :: KVAL_CR, KVAL_CI, KVAL2_CR, KVAL2_CI
      DOUBLE PRECISION  :: MOD_KVAL, SEPCON_CR, SEPCON_CI
      DOUBLE PRECISION  :: DENOM_CR, DENOM_CI, LFPD, LFND
      DOUBLE PRECISION  :: LFPD1, LFPD2, LFND1, LFND2
      CHARACTER (LEN=3) :: CN, CI

!      added by xliu for debugging
!      INTEGER count
!      data count /0/
!      SAVE count

!  Start of code
!  -------------
!      count = count + 1

!  initialise status

      STATUS  = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Layer and Fourier number

      N = GIVEN_LAYER
      M = FOURIER

!  nothing to do if not varying

      IF ( .NOT. DO_VARY ) RETURN

!  For solution saving option with no scattering,
!  linearized solution vectors are all zero
!  Exit the routine by going to continuation point.

      IF ( DO_SOLUTION_SAVING ) THEN
        IF ( .NOT. DO_LAYER_SCATTERING(M,N) ) THEN
          DO Q = 1, N_PARAMETERS
            DO I = 1, NSTREAMS_2
              DO O1 = 1, NSTOKES
                DO K = 1, K_REAL(N)
                  L_SOLA_XPOS(I,O1,K,N,Q)  = ZERO
                  L_SOLB_XNEG(I,O1,K,N,Q)  = ZERO
                ENDDO
              ENDDO
            ENDDO
            DO K = 1, K_REAL(N)
              L_KEIGEN(K,N,Q)  = ZERO
            ENDDO
          ENDDO
          GO TO 3456
        ENDIF
      ENDIF

!  Linearize the Eigenmatrix
!  =========================

!  set up linearizations of SAB and DAB. All real variables.
!   Thoroughly debugged, November 2005

      DO I = 1, NSTREAMS
       XINV = - ONE / QUAD_STREAMS(I)
       DO J = 1, NSTREAMS
        FAC =  QUAD_HALFWTS(J) * XINV
        DO Q = 1, N_PARAMETERS
         DO O1 = 1, NSTOKES
          DO O2 = 1, NSTOKES
           DP = ZERO
           DM = ZERO
           DO L = M, NMOMENTS
            H1P = ZERO
            H1M = ZERO
            DO O3 = 1, NSTOKES
             H2P = ZERO
             H2M = ZERO
             DO O4 = 1, NSTOKES
              L_OG = L_OMEGA_GREEK(L,N,O3,O4,Q)
              H2P = H2P + L_OG * PI_XQP(L,J,O4,O2)
              H2M = H2M + L_OG * PI_XQM_PRE(L,J,O4,O2)
             ENDDO
             H1P = H1P + PI_XQP(L,I,O1,O3)*H2P
             H1M = H1M + PI_XQP(L,I,O1,O3)*H2M
            ENDDO
            DP = DP + H1P
            DM = DM + H1M
           ENDDO
           L_SAB(I,J,O1,O2,N,Q) = FAC * ( DP + DM )
           L_DAB(I,J,O1,O2,N,Q) = FAC * ( DP - DM )
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

!  set up linearized eigenmatrices, saved for all layers. Real variables
!   Thoroughly debugged, November 2005

      DO Q = 1, N_PARAMETERS
       DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        DO J = 1, NSTREAMS
         JC = NSTOKES*(J-1)
         DO O1 = 1, NSTOKES
          IROW = IR + O1
          DO O2 = 1, NSTOKES
           JCOL = JC + O2
           H1 = ZERO
           DO K = 1, NSTREAMS
            H2 = ZERO
            DO O3 = 1, NSTOKES
             H2 = H2 + L_DAB(I,K,O1,O3,N,Q) *   SAB(K,J,O3,O2,N) &
                     +   DAB(I,K,O1,O3,N)   * L_SAB(K,J,O3,O2,N,Q)
            ENDDO
            H1 = H1 + H2
           ENDDO
           L_EIGENMAT(IROW,JCOL,N,Q) = H1
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     Linearization of Real Eigensolutions
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  start eigenvalue loop

      DO K = 1, K_REAL(N)

!  some basic quantities

        AA    = EIGENMASK_R(K)
        KVAL  = KEIGEN(K,N)
        KVAL2 = 2.0D0 * KVAL
        KSQD  = KVAL * KVAL

!   REAL SOLUTIONS ONLY (ASYMTX linearization)
!   ==========================================

       IF ( DO_REAL_EIGENSOLVER(M,N) ) THEN

!  Count the number of entries in the eigenvector ( = TC)
!    If the eigenvector is trivial (one entry, TC = 1).....
!    ......Zero the results and skip linearization for this eigenvalue

        TC = 0
        DO JCOL = 1,  NSTKS_NSTRMS
          if ( RITE_EVEC(JCOL,AA).NE.ZERO ) TC = TC + 1
        ENDDO
!Rob fix 4/12/12
        IF ( TC.EQ.1 .and. NSTKS_NSTRMS.gt.1 ) THEN
          DO Q = 1, N_PARAMETERS
            HVEC1A(1,Q) = ZERO
            DO I = 1, NSTKS_NSTRMS
              HVEC1A(I+1,Q) = ZERO
            ENDDO
          ENDDO
          GO TO 4555
        ENDIF

!  Nothing to do also if there is degeneracy
!    NEW CODE by R. Spurr 22 August 2008.

        IF ( EIGENDEGEN(K,N) ) THEN
          DO Q = 1, N_PARAMETERS
            HVEC1A(1,Q) = ZERO
            DO I = 1, NSTKS_NSTRMS
              HVEC1A(I+1,Q) = ZERO
            ENDDO
          ENDDO
          GO TO 4555
        ENDIF

!  initialise solution matrix HMAT1A (important to do this!)
!    HMAT1A is destroyed upon passage through the LAPACK linear modules

        DO I = 1, MAXSTRMSTKS_P1
!          LAPACK_IPIV1A(I) = 0
          DO J = 1, MAXSTRMSTKS_P1
            HMAT1A(I,J) = ZERO
          ENDDO
        ENDDO

!  Determine solution matrix HMAT (A in the matrix equation A.x = B)

        NSTMSTKS_P1 = NSTKS_NSTRMS + 1
        DO I = 1, NSTKS_NSTRMS
          DO J = 1, NSTKS_NSTRMS
            HMAT1A(I,J+1) = - EIGENMAT_SAVE(I,J,N)
          ENDDO
          HMAT1A(I,I+1) = HMAT1A(I,I+1) + KSQD
          HMAT1A(I,1)   = TWO * KVAL * RITE_EVEC(I,AA)
          HMAT1A(NSTMSTKS_P1,I+1) = RITE_EVEC(I,AA)
        ENDDO
        HMAT1A(NSTMSTKS_P1,1) = ZERO

!        HMAT1A(NSTMSTKS_P1,1) = 1.879d-15   ! shoudl technically be zer

!  solution column vectors (B in the matrix equation A.x = B).
!    (One for each parameter to be varied)

        DO Q = 1, N_PARAMETERS
          DO I = 1, NSTKS_NSTRMS
            H1 = ZERO
            DO J = 1, NSTKS_NSTRMS
              H1 = H1 + L_EIGENMAT(I,J,N,Q) * RITE_EVEC(J,AA)
            ENDDO
            HVEC1A(I,Q) = H1
          ENDDO
          HVEC1A(NSTMSTKS_P1,Q) = ZERO
        ENDDO

!  Find any degenerate columns J, and replace the (J,J) entry with
!   a small number. This is a quick fix that ensures that the LAPACK
!   routine does not fall over because of degeneracy.
!    --This problem is only present in the vector case.

!  Fiddle code

!        DO J = 1, NSTKS_NSTRMS
!          NC = 0
!          DO I = 1,  NSTKS_NSTRMS
!           IF ( HMAT1A(I,J+1).EQ.ZERO) NC = NC + 1
!          ENDDO
!          IF ( NC .EQ. NSTKS_NSTRMS ) THEN
!            do I = 1, NSTKS_NSTRMS
!               HMAT1A(J,I+1) = 1.237d-15*DBLE(i)
!            enddo
!            if ( HMAT1A(J+1,J+1).EQ.ZERO)
!     &    HMAT1A(J+1,J+1) = 1.056d-15
!          ENDIF
!        ENDDO

!  ultimate fiddle

!        DO J = 1, NSTMSTKS_P1
!           Do I = 1, NSTMSTKS_P1
!             IF (HMAT1A(I,J).EQ.ZERO )HMAT1A(I,J) = 1.0d-15*(I+J)
!           ENDDO
!        ENDDO

!  solve System for the linearized sum-vector + eigenvalue, together
!     Solve matrix system  using LAPACK modules DGETRF and DEGTRS

!        IF (m.eq.2.and.n.eq.18) THEN
!           WRITE(45, *) m, k, N, NSTMSTKS_P1, LAPACK_INFO
!           WRITE(45, '(17D24.12)')
!     &          ((HMAT1A(j, i), i = 1, NSTMSTKS_P1), j = 1, NSTMSTKS_P1
!c           if(lapack_info.ne.0)STOP'2,18,9'
!        ENDIF

!#######################################################################

!  Resume normal code

        CALL DGETRF  ( NSTMSTKS_P1, NSTMSTKS_P1, HMAT1A, &
                       MAXSTRMSTKS_P1,  LAPACK_IPIV1A, LAPACK_INFO )

!  Exception handling 1

        IF ( LAPACK_INFO .GT. 0 ) THEN
          WRITE(CI, '(I3)' ) LAPACK_INFO
          WRITE(CN, '(I3)' ) N
          MESSAGE = 'Singular matrix, u(i,i)=0, for i = '//CI
          TRACE   = 'DGETRF call # 1 in VLIDORT_L_QHOM_SOLUTION,'// &
                      ' layer # '//CN
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ELSE IF ( LAPACK_INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) LAPACK_INFO
          WRITE(CN, '(I3)' ) N
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGETRF call # 1 in VLIDORT_L_QHOM_SOLUTION,'// &
                      ' layer # '//CN
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

        CALL DGETRS ('N', NSTMSTKS_P1, N_PARAMETERS, HMAT1A, &
                     MAXSTRMSTKS_P1, LAPACK_IPIV1A,  HVEC1A, &
                     MAXSTRMSTKS_P1, LAPACK_INFO )

!  Exception handling 1

        IF ( LAPACK_INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) LAPACK_INFO
          WRITE(CN, '(I3)' ) N
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGETRS call # 1 in VLIDORT_L_QHOM_SOLUTION,'// &
                      ' layer # '//CN
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  Continuation point if linearized skipped

 4555   CONTINUE

!  Start loop over varying parameters, to assign answers

        DO Q = 1, N_PARAMETERS

          L_KEIGEN(K,N,Q) = HVEC1A(1,Q)
          DO I = 1, NSTREAMS
           IR = NSTOKES*(I-1)
           DO O1 = 1, NSTOKES
            IROW = IR + O1
            L_FWD_SUMVEC_R(I,O1,Q) = HVEC1A(IROW+1,Q)
           ENDDO
          ENDDO

!  Debug (important)
!   Tested 28 December 2005.

!         IF ( DO_DEBUG_WRITE ) THEN
!          IF (Q.EQ.1.and.m.eq.2.and.n.eq.1) THEN
!           write(37,'(2i3,L3,1p2e20.12)') &
!                   M,K,EIGENDEGEN(K,N),KEIGEN(K,N),L_KEIGEN(K,N,Q)
!           DO I = 1, NSTREAMS
!            DO O1 = 1, NSTOKES
!                 write(37,'(4i3,1p2e20.12)') &
!             m,k,i,o1,RITE_EVEC(I,k ),L_FWD_SUMVEC_R(I,O1,Q)
!            ENDDO
!           ENDDO
!          ENDIF
!         ENDIF

        ENDDO

!   REAL SOLUTIONS (DGEEV linearization)
!   ====================================

       ELSE

!  initialise

       I_NONTRIVIAL = 0

!  Count the number of entries in the eigenvector ( = TC)
!    If the eigenvector is trivial (one entry, TC = 1).....
!    ......Zero the results and skip linearization for this eigenvalue

        TC = 0
        DO JCOL = 1,  NSTKS_NSTRMS
          if ( RITE_EVEC(JCOL,AA).NE.ZERO ) TC = TC + 1
        ENDDO
        IF ( TC .EQ. 1 ) THEN
          DO Q = 1, N_PARAMETERS
            L_KEIGEN(K,N,Q)  = ZERO
            DO I = 1, NSTKS_NSTRMS
              HVEC1(I,Q) = ZERO
            ENDDO
          ENDDO
          GO TO 4556
        ENDIF

!  Establish the first non-trivial row
!      This is necessary to avoid underflow problems

        TR = 0
        DO JCOL = 1,  NSTKS_NSTRMS
          if ( DABS(RITE_EVEC(JCOL,AA)).GT.1.0D-10 ) THEN
            TR = JCOL
            GO TO 3442
          ENDIF
        ENDDO
 3442   CONTINUE

!  Use the adjoint method to get L(k) = <X^,L(G)X> / 2k <X^,X>

        SR  = ZERO
        DO JCOL = 1,  NSTKS_NSTRMS
          SR = SR + LEFT_EVEC(JCOL,AA) * RITE_EVEC(JCOL,AA)
        ENDDO
        MODULUS = ONE / KVAL2 / SR
        DO Q = 1, N_PARAMETERS
          DO IROW = 1, NSTKS_NSTRMS
            HR = ZERO
            DO JCOL = 1, NSTKS_NSTRMS
              HR = HR + L_EIGENMAT(IROW,JCOL,N,Q) * RITE_EVEC(JCOL,AA)
            ENDDO
            HVEC_AUX(IROW,Q) = HR
          ENDDO
          HR = ZERO
          DO JCOL = 1,  NSTKS_NSTRMS
            HR = HR + LEFT_EVEC(JCOL,AA) * HVEC_AUX(JCOL,Q)
          ENDDO
          L_KEIGEN(K,N,Q) = HR * MODULUS
        ENDDO

!  Now solve for the linearized eigenvector L(X) from A.L(X) = B
!     Determine matrix HMAT (=A) and vectors HVEC (B)

!  Initialise first, as HMAT is destroyed after LAPACK routines

        DO IROW = 1, NSTKS_NSTRMS
          DO JCOL = 1, NSTKS_NSTRMS
            HMAT1(IROW,JCOL) = ZERO
          ENDDO
        ENDDO

!  set up the matrix HMAT
!    Use the linearization of the Eigen-equation for all but one
!    of the rows. For the FIRST row, use the linearization of the
!    normalization condition on the eigenvectors.

        DO IROW = 1, NSTKS_NSTRMS
          IF ( IROW.NE.TR ) THEN
           DO JCOL = 1, NSTKS_NSTRMS
            HMAT1(IROW,JCOL) = - EIGENMAT_SAVE(IROW,JCOL,N)
           ENDDO
           HMAT1(IROW,IROW) = HMAT1(IROW,IROW) + KSQD
          ELSE
           DO JCOL = 1, NSTKS_NSTRMS
            HMAT1(IROW,JCOL) = RITE_EVEC(JCOL,AA)
           ENDDO
          ENDIF
        ENDDO

!  set up the HVEC:
!    Use the linearization of the Eigen-equation for all but one
!    of the rows. For the FIRST row, use the linearization of the
!    normalization condition on the eigenvectors.

        DO Q = 1, N_PARAMETERS
          I_NONTRIVIAL = 0
          DO IROW = 1, NSTKS_NSTRMS
           IF ( IROW.NE.TR ) THEN
            HVEC1(IROW,Q) = HVEC_AUX(IROW,Q) - &
               KVAL2 * L_KEIGEN(K,N,Q) * RITE_EVEC(IROW,AA)
            IF ( DABS(HVEC1(IROW,Q)).LT.1.0D-12)THEN
              HVEC1(IROW,Q) = ZERO
            ELSE
              I_NONTRIVIAL = I_NONTRIVIAL + 1
            ENDIF
           ELSE
            HVEC1(IROW,Q) = ZERO
           ENDIF
          ENDDO
        ENDDO

!  Solve matrix system  using LAPACK modules DGETRF and DEGTRS
!  Only do if vector has non-trivial entries

        IF ( I_NONTRIVIAL.GT.0) THEN

          CALL DGETRF  ( NSTKS_NSTRMS, NSTKS_NSTRMS, HMAT1, &
                       MAXSTRMSTKS,  LAPACK_IPIV1, LAPACK_INFO )

!  Exception handling 2

          IF ( LAPACK_INFO .GT. 0 ) THEN
            WRITE(CI, '(I3)' ) LAPACK_INFO
            WRITE(CN, '(I3)' ) N
            MESSAGE = 'Singular matrix, u(i,i)=0, for i = '//CI
            TRACE   = 'DGETRF call # 2 in VLIDORT_L_QHOM_SOLUTION,'// &
                      ' layer # '//CN
            STATUS  = VLIDORT_SERIOUS
            RETURN
          ELSE IF ( LAPACK_INFO .LT. 0 ) THEN
            WRITE(CI, '(I3)' ) LAPACK_INFO
            WRITE(CN, '(I3)' ) N
            MESSAGE = 'argument i illegal value, for i = '//CI
            TRACE   = 'DGETRF call # 2 in VLIDORT_L_QHOM_SOLUTION,'// &
                      ' layer # '//CN
            STATUS  = VLIDORT_SERIOUS
            RETURN
          ENDIF

          CALL DGETRS ('N', NSTKS_NSTRMS, N_PARAMETERS, HMAT1, &
                     MAXSTRMSTKS, LAPACK_IPIV1,       HVEC1, &
                     MAXSTRMSTKS, LAPACK_INFO )

!  Exception handling 2

          IF ( LAPACK_INFO .LT. 0 ) THEN
            WRITE(CI, '(I3)' ) LAPACK_INFO
            WRITE(CN, '(I3)' ) N
            MESSAGE = 'argument i illegal value, for i = '//CI
            TRACE   = 'DGETRS call # 2 in VLIDORT_L_QHOM_SOLUTION,'// &
                      ' layer # '//CN
            STATUS  = VLIDORT_SERIOUS
            RETURN
          ENDIF

        ENDIF

!  Continuation point if linearization is skipped

 4556   CONTINUE

!  Start loop over varying parameters, and assign
!    linearized sum vector = linearized eigensolution

        DO Q = 1, N_PARAMETERS

         DO I = 1, NSTREAMS
           IR = NSTOKES*(I-1)
           DO O1 = 1, NSTOKES
            IROW = IR + O1
            L_FWD_SUMVEC_R(I,O1,Q) = HVEC1(IROW,Q)
           ENDDO
         ENDDO

!  Debug (important)
!   Tested 28 December 2005.

!         IF ( DO_DEBUG_WRITE ) THEN
!          IF (Q.EQ.1.and.m.eq.0.and.n.eq.1) THEN
!           write(37,'(2i3,1p2e20.12)')
!     &              M,K,KEIGEN(K,N),L_KEIGEN(K,N,Q)
!           DO I = 1, NSTREAMS
!            DO O1 = 1, NSTOKES
!                 write(37,'(4i3,1p2e20.12)')
!     &        m,k,i,o1,FWD_SUMVEC(I,O1,K),L_FWD_SUMVEC_R(I,O1,Q)
!            ENDDO
!           ENDDO
!          ENDIF
!         ENDIF

        ENDDO

!  End Clause ASYMTX vs. DGEEV

       ENDIF

!  Resume General Code
!  ===================

!  Start loop over varying parameters

       DO Q = 1, N_PARAMETERS

!  linearized difference vector

         DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
           IROW = IR + O1
           SR = ZERO
           DO J = 1, NSTREAMS
            JC = NSTOKES*(J-1)
            DO O2 = 1, NSTOKES
             JCOL = JC + O2
             SR = SR - L_SAB(I,J,O1,O2,N,Q) *   FWD_SUMVEC(J,O2,K) &
                     -   SAB(I,J,O1,O2,N)   * L_FWD_SUMVEC_R(J,O2,Q)
            ENDDO
           ENDDO
           HR = ( - SR - L_KEIGEN(K,N,Q) * FWD_DIFVEC(I,O1,K) )
           L_FWD_DIFVEC_R(I,O1) = HR / KVAL
          ENDDO
         ENDDO

!  assign Forward solutions PHI (Siewert's notation)
!    --->  first N are "DOWN", last N are "UP" (streams)
!    --->  Use symmetry properties to set -ve eigensolutions

         DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES
           LFPD = HALF * ( L_FWD_SUMVEC_R(I,O1,Q)+L_FWD_DIFVEC_R(I,O1) )
           LFND = HALF * ( L_FWD_SUMVEC_R(I,O1,Q)-L_FWD_DIFVEC_R(I,O1) )
           L_SOLA_XPOS(I,O1,K,N,Q)  = LFPD
           L_SOLB_XNEG(I,O1,K,N,Q)  = LFND
           L_SOLA_XPOS(I1,O1,K,N,Q)  = LFND
           L_SOLB_XNEG(I1,O1,K,N,Q)  = LFPD
           IF ( O1 .GT. 2 ) THEN
            L_SOLA_XPOS(I1,O1,K,N,Q)  = - LFND
            L_SOLB_XNEG(I1,O1,K,N,Q)  = - LFPD
           ENDIF
          ENDDO
         ENDDO

!  Older code; Uses additional wasteful holding arrays
!  solution assignation (only for BVP)

!         DO I = 1, NSTREAMS
!          I1 = I + NSTREAMS
!          DO O1 = 1, NSTOKES
!           L_FWD_XPOS_R(I,O1,K,N,Q)  = HALF *
!     &          ( L_FWD_SUMVEC_R(I,O1,Q) + L_FWD_DIFVEC_R(I,O1) )
!           L_FWD_XPOS_R(I1,O1,K,N,Q) = HALF *
!     &          ( L_FWD_SUMVEC_R(I,O1,Q) - L_FWD_DIFVEC_R(I,O1) )
!           L_FWD_XNEG_R(I1,O1,K,N,Q) =  L_FWD_XPOS_R(I,O1,K,N,Q)
!           L_FWD_XNEG_R(I,O1,K,N,Q)  =  L_FWD_XPOS_R(I1,O1,K,N,Q)
!          ENDDO
!         ENDDO
!         DO I = 1, NSTREAMS
!          I1 = I + NSTREAMS
!          DO O1 = 1, NSTOKES
!           L_SOLA_XPOS_R(I,O1,K,N,Q)  = L_FWD_XPOS_R(I,O1,K,N,Q)
!           L_SOLB_XNEG_R(I,O1,K,N,Q)  = L_FWD_XNEG_R(I,O1,K,N,Q)
!           H1 = ZERO
!           H2 = ZERO
!           DO O2 = 1, NSTOKES
!            H1 = H1 + DMAT(O1,O2)*L_FWD_XNEG_R(I,O2,K,N,Q)
!            H2 = H2 + DMAT(O1,O2)*L_FWD_XPOS_R(I,O2,K,N,Q)
!           ENDDO
!           L_SOLA_XPOS_R(I1,O1,K,N,Q)  = H1
!           L_SOLB_XNEG_R(I1,O1,K,N,Q)  = H2
!          ENDDO
!         ENDDO

!  End parameter loop

       ENDDO

!  End real eigenvalue loop

      ENDDO

!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     Linearization of COMPLEX Eigensolutions
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Only if flagged

      IF ( .NOT. DO_REAL_EIGENSOLVER(M,N) ) THEN

!  Offset

       KO1 = K_REAL(N) + 1

!  start loop over eigensolutions

       DO K = 1, K_COMPLEX(N)

!  Bookkeeping indices

        K0 = 2 * ( K - 1 )
        K1 = KO1 + K0
        K2 = K1  + 1

!  Basic quantities

        AA  = EIGENMASK_C(K)
        AA1 = AA + 1
        KVAL_CR  = KEIGEN(K1,N)
        KVAL_CI  = KEIGEN(K2,N)
        KVAL2_CR = TWO * KVAL_CR
        KVAL2_CI = TWO * KVAL_CI
        MOD_KVAL = KEIGEN_CSQ(K)
        SEPCON_CR =   KEIGEN(K1,N) / MOD_KVAL
        SEPCON_CI = - KEIGEN(K2,N) / MOD_KVAL

!  Establish the zero imaginary part of one eigenvector component

        TR = 0
        DO JCOL = 1,  NSTKS_NSTRMS
          if ( RITE_EVEC(JCOL,AA1).EQ.ZERO ) THEN
            TR = JCOL + NSTKS_NSTRMS
            GO TO 3443
          ENDIF
        ENDDO
 3443   CONTINUE

!  Use the adjoint theory to get L(k) from 2kL(k) = <X^,L(G)X> / <X^,X>
!      MUST USE THE COMPLEX CONJUGATE OF THE LEFT EIGENVECTOR.

!  --------------- Adjoint  norm

        SCR = ZERO
        SCI = ZERO
        DO JCOL = 1,  NSTKS_NSTRMS
          HCR = + LEFT_EVEC(JCOL,AA)  * RITE_EVEC(JCOL,AA) &
                + LEFT_EVEC(JCOL,AA1) * RITE_EVEC(JCOL,AA1)
          SCR = SCR + HCR
          HCI = - LEFT_EVEC(JCOL,AA1) * RITE_EVEC(JCOL,AA) &
                + LEFT_EVEC(JCOL,AA)  * RITE_EVEC(JCOL,AA1)
          SCI = SCI + HCI
          ADJN_CR(JCOL) = + KVAL2_CR * RITE_EVEC(JCOL,AA) &
                          - KVAL2_CI * RITE_EVEC(JCOL,AA1)
          ADJN_CI(JCOL) = + KVAL2_CI * RITE_EVEC(JCOL,AA) &
                          + KVAL2_CR * RITE_EVEC(JCOL,AA1)
        ENDDO
        DENOM_CR = KVAL2_CR * SCR - KVAL2_CI * SCI
        DENOM_CI = KVAL2_CR * SCI + KVAL2_CI * SCR
        MODULUS  = ONE / ( DENOM_CR*DENOM_CR + DENOM_CI * DENOM_CI )

!  ---------------- Adjoint numerator, and final computation

        DO Q = 1, N_PARAMETERS
         DO IROW = 1, NSTKS_NSTRMS
          HCR = ZERO
          HCI = ZERO
          DO JCOL = 1, NSTKS_NSTRMS
           HCR = HCR + L_EIGENMAT(IROW,JCOL,N,Q) * RITE_EVEC(JCOL,AA)
           HCI = HCI + L_EIGENMAT(IROW,JCOL,N,Q) * RITE_EVEC(JCOL,AA1)
          ENDDO
          HVEC_CR(IROW,Q) = HCR
          HVEC_CI(IROW,Q) = HCI
         ENDDO
         HCR = ZERO
         HCI = ZERO
         DO JCOL = 1,  NSTKS_NSTRMS
          HCR = HCR + LEFT_EVEC(JCOL,AA)  * HVEC_CR(JCOL,Q) &
                    + LEFT_EVEC(JCOL,AA1) * HVEC_CI(JCOL,Q)
          HCI = HCI - LEFT_EVEC(JCOL,AA1) * HVEC_CR(JCOL,Q) &
                    + LEFT_EVEC(JCOL,AA)  * HVEC_CI(JCOL,Q)
         ENDDO
         L_KEIGEN(K1,N,Q) = ( HCR*DENOM_CR + HCI*DENOM_CI )*MODULUS
         L_KEIGEN(K2,N,Q) = ( HCI*DENOM_CR - HCR*DENOM_CI )*MODULUS
        ENDDO

!  Now solve for the linearized eigenvector L(X) from A.L(X) = B
!     Determine matrix HMAT (=A) and vectors HVEC (B)

!  ---------- Initialise first, as HMAT is destroyed after LAPACK routin

        DO IROW = 1, NSTKS_NSTRMS_2
          DO JCOL = 1, NSTKS_NSTRMS_2
            HMAT2(IROW,JCOL) = ZERO
          ENDDO
        ENDDO

!  ---------- set up matrix for Full solution

        DO I = 1, NSTKS_NSTRMS
          I1 = I + NSTKS_NSTRMS
          DO J = 1, NSTKS_NSTRMS
            J1 = J + NSTKS_NSTRMS
            HMAT2(I,J)   = - EIGENMAT_SAVE(I,J,N)
            HMAT2(I1,J1) = - EIGENMAT_SAVE(I,J,N)
          ENDDO
          HMAT2(I,I)   = HMAT2(I,I)   + REAL_KSQ(AA)
          HMAT2(I1,I1) = HMAT2(I1,I1) + REAL_KSQ(AA)
          HMAT2(I,I1)  = - IMAG_KSQ(AA)
          HMAT2(I1,I)  = + IMAG_KSQ(AA)
        ENDDO

!  ---------- set up the vectors for full solution

        DO Q = 1, N_PARAMETERS
          DO I = 1, NSTKS_NSTRMS
           I1 = I + NSTKS_NSTRMS
           HCR = + L_KEIGEN(K1,N,Q) * ADJN_CR(I) &
                 - L_KEIGEN(K2,N,Q) * ADJN_CI(I)
           HCI = + L_KEIGEN(K1,N,Q) * ADJN_CI(I) &
                 + L_KEIGEN(K2,N,Q) * ADJN_CR(I)
           HVEC2(I,Q)  = HVEC_CR(I,Q) - HCR
           HVEC2(I1,Q) = HVEC_CI(I,Q) - HCI
          ENDDO
        ENDDO

!  ---------- Strike out TR row/column to get reduced matrix and vectors

!mick mod
        !NSTKS_NSTRMS_2  = 2 * NSTKS_NSTRMS
        NSTKS_NSTRMS_21 = 2 * NSTKS_NSTRMS - 1
        DO I = 1, TR - 1
          DO J = 1, TR - 1
            HMATRED(I,J) = HMAT2(I,J)
          ENDDO
          DO J = TR + 1, NSTKS_NSTRMS_2
            HMATRED(I,J-1) = HMAT2(I,J)
          ENDDO
          DO Q = 1, N_PARAMETERS
            HVECRED(I,Q) = HVEC2(I,Q)
          ENDDO
        ENDDO
        DO I = TR + 1, NSTKS_NSTRMS_2
          DO J = 1, TR - 1
            HMATRED(I-1,J) = HMAT2(I,J)
          ENDDO
          DO J = TR + 1, NSTKS_NSTRMS_2
            HMATRED(I-1,J-1) = HMAT2(I,J)
          ENDDO
          DO Q = 1, N_PARAMETERS
            HVECRED(I-1,Q) = HVEC2(I,Q)
          ENDDO
        ENDDO

!  ----------  Replace the last row by linearized normalization constrai
!                ( IF (TR.LT.NSTKS_NSTRMS_2) THEN  possible exception...
!                Should not be necessary to flag this case )
!
        ILAST =  NSTKS_NSTRMS_21
        DO J = 1, NSTKS_NSTRMS
          HMATRED(ILAST,J) = RITE_EVEC(J,AA)
        ENDDO
        DO J = NSTKS_NSTRMS+1, TR - 1
          HMATRED(ILAST,J) = RITE_EVEC(J-NSTKS_NSTRMS,AA1)
        ENDDO
        DO J = TR + 1, NSTKS_NSTRMS_2
          HMATRED(ILAST,J-1) = RITE_EVEC(J-NSTKS_NSTRMS,AA1)
        ENDDO
        DO Q = 1, N_PARAMETERS
          HVECRED(ILAST,Q) = ZERO
        ENDDO

!  ---------- Solve matrix system using LAPACK modules DGETRF and DEGTRS

        CALL DGETRF  ( NSTKS_NSTRMS_21, NSTKS_NSTRMS_21, HMATRED, &
                       MAXSTRMSTKS_21,  LAPACK_IPIV21,   LAPACK_INFO )

!  Exception handling 3

        IF ( LAPACK_INFO .GT. 0 ) THEN
          WRITE(CI, '(I3)' ) LAPACK_INFO
          WRITE(CN, '(I3)' ) N
          MESSAGE = 'Singular matrix, u(i,i)=0, for i = '//CI
          TRACE   = 'DGETRF call # 3 in VLIDORT_L_QHOM_SOLUTION,'// &
                      ' layer # '//CN
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ELSE IF ( LAPACK_INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) LAPACK_INFO
          WRITE(CN, '(I3)' ) N
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGETRF call # 3 in VLIDORT_L_QHOM_SOLUTION,'// &
                      ' layer # '//CN
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

        CALL DGETRS ('N', NSTKS_NSTRMS_21, N_PARAMETERS, HMATRED, &
                     MAXSTRMSTKS_21, LAPACK_IPIV21,      HVECRED, &
                     MAXSTRMSTKS_21, LAPACK_INFO )

!  Exception handling 3

        IF ( LAPACK_INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) LAPACK_INFO
          WRITE(CN, '(I3)' ) N
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGETRS call # 3 in VLIDORT_L_QHOM_SOLUTION,'// &
                      ' layer # '//CN
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  Start loop over varying parameters

        DO Q = 1, N_PARAMETERS

!  linearized sum vector = linearized eigensolution

         DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
           IROW  = IR + O1
           IROW1 = IROW + NSTKS_NSTRMS
           L_FWD_SUMVEC_CR(I,O1) = HVECRED(IROW,Q)
           IF ( IROW1.LT.TR)THEN
             L_FWD_SUMVEC_CI(I,O1) = HVECRED(IROW1,Q)
           ELSE IF (IROW1.EQ.TR) THEN
             L_FWD_SUMVEC_CI(I,O1) = ZERO
           ELSE
             L_FWD_SUMVEC_CI(I,O1) = HVECRED(IROW1-1,Q)
           ENDIF
          ENDDO
         ENDDO

!  Debug code (important)

!         IF (DO_DEBUG_WRITE) THEN
!          IF ( Q.EQ.1.AND.M.EQ.0.AND.N.EQ.1) THEN
!           WRITE(91,'(2I3,1P4E20.10)')
!     &            M,K,KEIGEN(K1,N),KEIGEN(K2,N),
!     &            L_KEIGEN(K1,N,Q),L_KEIGEN(K2,N,Q)
!           DO I = 1, NSTREAMS
!            DO O1 = 1, NSTOKES
!              WRITE(91,'(4I3,1P4E20.10)')
!     &        M,K,I,O1,FWD_SUMVEC(I,O1,K1),FWD_SUMVEC(I,O1,K2),
!     &         L_FWD_SUMVEC_CR(I,O1),L_FWD_SUMVEC_CI(I,O1)
!            ENDDO
!           ENDDO
!          ENDIF
!         ENDIF

!  Determine the linearized difference vector

         DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
           IROW = IR + O1
           SCR = ZERO
           SCI = ZERO
           DO J = 1, NSTREAMS
            JC = NSTOKES*(J-1)
            DO O2 = 1, NSTOKES
             JCOL = JC + O2
             SCR = SCR - L_SAB(I,J,O1,O2,N,Q) *   FWD_SUMVEC(J,O2,K1) &
                       -   SAB(I,J,O1,O2,N)   * L_FWD_SUMVEC_CR(J,O2)
             SCI = SCI - L_SAB(I,J,O1,O2,N,Q) *   FWD_SUMVEC(J,O2,K2) &
                       -   SAB(I,J,O1,O2,N)   * L_FWD_SUMVEC_CI(J,O2)
            ENDDO
           ENDDO
           HCR = - SCR - L_KEIGEN(K1,N,Q) * FWD_DIFVEC(I,O1,K1) &
                       + L_KEIGEN(K2,N,Q) * FWD_DIFVEC(I,O1,K2)
           HCI = - SCI - L_KEIGEN(K2,N,Q) * FWD_DIFVEC(I,O1,K1) &
                       - L_KEIGEN(K1,N,Q) * FWD_DIFVEC(I,O1,K2)
           L_FWD_DIFVEC_CR(I,O1) = HCR*SEPCON_CR - HCI*SEPCON_CI
           L_FWD_DIFVEC_CI(I,O1) = HCI*SEPCON_CR + HCR*SEPCON_CI
          ENDDO
         ENDDO

!  assign Forward solutions PHI (Siewert's notation)
!    --->  first N are "DOWN", last N are "UP" (streams)
!    --->  Use symmetry properties to set -ve eigensolutions

         DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES
           LFND1 = HALF * (L_FWD_SUMVEC_CR(I,O1)-L_FWD_DIFVEC_CR(I,O1))
           LFPD1 = HALF * (L_FWD_SUMVEC_CR(I,O1)+L_FWD_DIFVEC_CR(I,O1))
           LFND2 = HALF * (L_FWD_SUMVEC_CI(I,O1)-L_FWD_DIFVEC_CI(I,O1))
           LFPD2 = HALF * (L_FWD_SUMVEC_CI(I,O1)+L_FWD_DIFVEC_CI(I,O1))
           L_SOLA_XPOS(I,O1,K1,N,Q)  = LFPD1
           L_SOLB_XNEG(I,O1,K1,N,Q)  = LFND1
           L_SOLA_XPOS(I,O1,K2,N,Q)  = LFPD2
           L_SOLB_XNEG(I,O1,K2,N,Q)  = LFND2
           L_SOLA_XPOS(I1,O1,K1,N,Q)  = LFND1
           L_SOLB_XNEG(I1,O1,K1,N,Q)  = LFPD1
           L_SOLA_XPOS(I1,O1,K2,N,Q)  = LFND2
           L_SOLB_XNEG(I1,O1,K2,N,Q)  = LFPD2
           IF ( O1 .GT. 2 ) THEN
            L_SOLA_XPOS(I1,O1,K1,N,Q)  = - LFND1
            L_SOLB_XNEG(I1,O1,K1,N,Q)  = - LFPD1
            L_SOLA_XPOS(I1,O1,K2,N,Q)  = - LFND2
            L_SOLB_XNEG(I1,O1,K2,N,Q)  = - LFPD2
           ENDIF
          ENDDO
         ENDDO

!  Older code: Uses additional wasteful holding arrays
!  solution assignation (only for BVP)

!         DO I = 1, NSTREAMS
!          I1 = I + NSTREAMS
!          DO O1 = 1, NSTOKES
!           L_FWD_XPOS_CR(I1,O1,K,N,Q) = HALF *
!     &          ( L_FWD_SUMVEC_CR(I,O1) - L_FWD_DIFVEC_CR(I,O1) )
!           L_FWD_XPOS_CR(I,O1,K,N,Q)  = HALF *
!     &          ( L_FWD_SUMVEC_CR(I,O1) + L_FWD_DIFVEC_CR(I,O1) )
!           L_FWD_XPOS_CI(I1,O1,K,N,Q) = HALF *
!     &          ( L_FWD_SUMVEC_CI(I,O1) - L_FWD_DIFVEC_CI(I,O1) )
!           L_FWD_XPOS_CI(I,O1,K,N,Q)  = HALF *
!     &          ( L_FWD_SUMVEC_CI(I,O1) + L_FWD_DIFVEC_CI(I,O1) )
!           L_FWD_XNEG_CR(I1,O1,K,N,Q) =  L_FWD_XPOS_CR(I,O1,K,N,Q)
!           L_FWD_XNEG_CR(I,O1,K,N,Q)  =  L_FWD_XPOS_CR(I1,O1,K,N,Q)
!           L_FWD_XNEG_CI(I1,O1,K,N,Q) =  L_FWD_XPOS_CI(I,O1,K,N,Q)
!           L_FWD_XNEG_CI(I,O1,K,N,Q)  =  L_FWD_XPOS_CI(I1,O1,K,N,Q)
!          ENDDO
!         ENDDO
!         DO I = 1, NSTREAMS
!          I1 = I + NSTREAMS
!          DO O1 = 1, NSTOKES
!           L_SOLA_XPOS_CR(I,O1,K,N,Q)  = L_FWD_XPOS_CR(I,O1,K,N,Q)
!           L_SOLB_XNEG_CR(I,O1,K,N,Q)  = L_FWD_XNEG_CR(I,O1,K,N,Q)
!           L_SOLA_XPOS_CI(I,O1,K,N,Q)  = L_FWD_XPOS_CI(I,O1,K,N,Q)
!           L_SOLB_XNEG_CI(I,O1,K,N,Q)  = L_FWD_XNEG_CI(I,O1,K,N,Q)
!           H1 = ZERO
!           H2 = ZERO
!           DO O2 = 1, NSTOKES
!            H1 = H1 + DMAT(O1,O2)*L_FWD_XNEG_CR(I,O2,K,N,Q)
!            H2 = H2 + DMAT(O1,O2)*L_FWD_XPOS_CR(I,O2,K,N,Q)
!           ENDDO
!           L_SOLA_XPOS_CR(I1,O1,K,N,Q)  = H1
!           L_SOLB_XNEG_CR(I1,O1,K,N,Q)  = H2
!           H1 = ZERO
!           H2 = ZERO
!           DO O2 = 1, NSTOKES
!            H1 = H1 + DMAT(O1,O2)*L_FWD_XNEG_CI(I,O2,K,N,Q)
!            H2 = H2 + DMAT(O1,O2)*L_FWD_XPOS_CI(I,O2,K,N,Q)
!           ENDDO
!           L_SOLA_XPOS_CI(I1,O1,K,N,Q)  = H1
!           L_SOLB_XNEG_CI(I1,O1,K,N,Q)  = H2
!          ENDDO
!         ENDDO

!  End parameter loop

        ENDDO

!  End complex eigenvalue loop

       ENDDO

!  End clause for not using the real eigensolver

      ENDIF

!  control point for avoiding eigensolver
!  ======================================

 3456 CONTINUE

!  Linearized Eigenstream transmittance factors for whole layer
!  ------------------------------------------------------------

!  When the solution saving option is set, then if there is no
!  scattering in this layer, then linearized transmittances are
!  linearized discrete ordinate transmittances.

      IF ( DO_SOLUTION_SAVING .AND. &
                     .NOT.DO_LAYER_SCATTERING(M,N) ) THEN

        DO Q = 1, N_PARAMETERS
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              L_T_DELT_EIGEN(IROW,N,Q) = L_T_DELT_DISORDS(I,N,Q)
            ENDDO
          ENDDO
        ENDDO

!  Otherwise compute them as normal

      ELSE

!  Real

        DO K = 1, K_REAL(N)
          TBAS = - T_DELT_EIGEN(K,N) * DELTAU_VERT(N)
          DO Q = 1, N_PARAMETERS
            LTQ = L_KEIGEN(K,N,Q) + KEIGEN(K,N)*L_DELTAU_VERT(Q,N)
            L_T_DELT_EIGEN(K,N,Q) = TBAS * LTQ
          ENDDO
        ENDDO

!  Complex

        DO K = 1, K_COMPLEX(N)
          K0 = 2 * K - 2
          K1 = KO1 + K0
          K2 = K1  + 1
          TBAS_CR = - T_DELT_EIGEN(K1,N) *  DELTAU_VERT(N)
          TBAS_CI = - T_DELT_EIGEN(K2,N) *  DELTAU_VERT(N)
          DO Q = 1, N_PARAMETERS
            LCR = L_KEIGEN(K1,N,Q) + KEIGEN(K1,N) * L_DELTAU_VERT(Q,N)
            LCI = L_KEIGEN(K2,N,Q) + KEIGEN(K2,N) * L_DELTAU_VERT(Q,N)
            L_T_DELT_EIGEN(K1,N,Q) = TBAS_CR * LCR - TBAS_CI * LCI
            L_T_DELT_EIGEN(K2,N,Q) = TBAS_CI * LCR + TBAS_CR * LCI
          ENDDO
        ENDDO

      ENDIF

!  debug
!      if ( m.eq.0.and.n.eq.20 ) then
!       DO K = 1, min(6,K_REAL(N))
!         WRITE(*,'(2i4,1p6e24.12)')m,k,L_T_DELT_EIGEN(K,N,1)
!       ENDDO
!      ENDIF

!  Eigenstream transmittance factors for partial layers
!  ----------------------------------------------------

!  Done irrespective of solution saving
!    Loop over user optical depths

      DO UTA = 1, N_USER_LEVELS
        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
         UT = PARTLAYERS_OUTINDEX(UTA)
         NA  = PARTLAYERS_LAYERIDX(UT)

!  If the layer of occurrence is the given layer

         IF ( NA .EQ. N ) THEN

!  When the solution saving option is set, then if there is no
!  scattering in this layer, then transmittances are just the
!  discrete ordinate transmittances.

          IF ( DO_SOLUTION_SAVING .AND. &
               .NOT.DO_LAYER_SCATTERING(M,N) ) THEN

           DO Q = 1, N_PARAMETERS
            DO I = 1, NSTREAMS
             IR = NSTOKES*(I-1)
             DO O1 = 1, NSTOKES
              IROW = IR + O1
              L_T_UTDN_EIGEN(IROW,UT,Q)=L_T_DISORDS_UTDN(I,UT,Q)
              L_T_UTUP_EIGEN(IROW,UT,Q)=L_T_DISORDS_UTUP(I,UT,Q)
             ENDDO
            ENDDO
           ENDDO

!  Otherwise, Compute the Eigenstream transmittance factors

          ELSE

!  Partial layer optical depths

           TAU_DN = PARTAU_VERT(UT)
           TAU_UP = DELTAU_VERT(N) - TAU_DN

!  Real eigenvalues

           DO K = 1, K_REAL(N)
            H_DN = - TAU_DN * T_UTDN_EIGEN(K,UT)
            H_UP = - TAU_UP * T_UTUP_EIGEN(K,UT)
            DO Q = 1, N_PARAMETERS
              VAR = L_KEIGEN(K,N,Q) + KEIGEN(K,N) * L_DELTAU_VERT(Q,N)
              L_T_UTDN_EIGEN(K,UT,Q) = H_DN * VAR
              L_T_UTUP_EIGEN(K,UT,Q) = H_UP * VAR
            ENDDO
           ENDDO

!  Complex eigenvalues
!    Bug fixed 16 December 2005

           DO K = 1, K_COMPLEX(N)
            K0 = 2 * ( K - 1 )
            K1 = KO1 + K0
            K2 = K1  + 1
            TBASUP_CR = - TAU_UP * T_UTUP_EIGEN(K1,UT)
            TBASUP_CI = - TAU_UP * T_UTUP_EIGEN(K2,UT)
            TBASDN_CR = - TAU_DN * T_UTDN_EIGEN(K1,UT)
            TBASDN_CI = - TAU_DN * T_UTDN_EIGEN(K2,UT)
            DO Q = 1, N_PARAMETERS
             LCR = L_KEIGEN(K1,N,Q) + KEIGEN(K1,N) * L_DELTAU_VERT(Q,N)
             LCI = L_KEIGEN(K2,N,Q) + KEIGEN(K2,N) * L_DELTAU_VERT(Q,N)
             L_T_UTDN_EIGEN(K1,UT,Q) = TBASDN_CR*LCR - TBASDN_CI*LCI
             L_T_UTDN_EIGEN(K2,UT,Q) = TBASDN_CI*LCR + TBASDN_CR*LCI
             L_T_UTUP_EIGEN(K1,UT,Q) = TBASUP_CR*LCR - TBASUP_CI*LCI
             L_T_UTUP_EIGEN(K2,UT,Q) = TBASUP_CI*LCR + TBASUP_CR*LCI
            ENDDO
           ENDDO

!  End clause for using Discrete ordinate valuaes or not

          ENDIF

!  End loop over partial layer depths

         ENDIF
        ENDIF
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_L_QHOM_SOLUTION

!

      SUBROUTINE VLIDORT_L_UHOM_SOLUTION ( &
        GIVEN_LAYER, FOURIER, DO_VARY, &
        N_PARAMETERS, &
        DO_UPWELLING, DO_DNWELLING, &
        DO_DEBUG_WRITE, &
        NSTOKES, NSTREAMS, &
        QUAD_HALFWTS, NMOMENTS, &
        N_USER_STREAMS, DO_LAYER_SCATTERING, &
        LOCAL_UM_START, &
        USER_SECANTS, STERM_LAYERMASK_UP, &
        STERM_LAYERMASK_DN, &
        OMEGA_GREEK, PI_XQP, PI_XUP, &
        PI_XUM, PI_XQM_PRE, &
        PI_XUM_POST, PI_XUP_PRE, &
        K_REAL, K_COMPLEX, &
        KEIGEN, HELPSTOKES, &
        UHOM_DNDN, UHOM_DNUP, &
        UHOM_UPDN, UHOM_UPUP, &
        ZETA_M, ZETA_P, L_KEIGEN, &
        L_SOLA_XPOS, L_SOLB_XNEG, &
        L_OMEGA_GREEK, &
        L_UHOM_DNDN, L_UHOM_DNUP, &
        L_UHOM_UPDN, L_UHOM_UPUP, &
        L_ZETA_M, L_ZETA_P )

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Regular inputs (strictly in)

      INTEGER, INTENT (IN) ::          GIVEN_LAYER
      INTEGER, INTENT (IN) ::          FOURIER
      LOGICAL, INTENT (IN) ::          DO_VARY
      INTEGER, INTENT (IN) ::          N_PARAMETERS
      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_DNWELLING
      LOGICAL, INTENT (IN) ::          DO_DEBUG_WRITE
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      DOUBLE PRECISION, INTENT (IN) :: QUAD_HALFWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::          NMOMENTS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      LOGICAL, INTENT (IN) ::          DO_LAYER_SCATTERING &
          ( 0:MAXMOMENTS, MAXLAYERS )
      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      DOUBLE PRECISION, INTENT (IN) :: USER_SECANTS  ( MAX_USER_STREAMS )
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_DN ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: OMEGA_GREEK &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_XQP &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_XUP &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_XUM &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_XQM_PRE &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_XUM_POST &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_XUP_PRE &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )

!  Regular values (Ostensibly in)

      INTEGER, INTENT (INOUT) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (INOUT) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: KEIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: HELPSTOKES &
          ( 0:MAXMOMENTS, MAXEVALUES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) :: UHOM_DNDN &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: UHOM_DNUP &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: UHOM_UPDN &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: UHOM_UPUP &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: ZETA_M &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: ZETA_P &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )

! Linearized values (Ostensibly in)
 
      DOUBLE PRECISION, INTENT (INOUT) :: L_KEIGEN &
          ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: L_SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: L_SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized values (Strictly in)

      DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_GREEK &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )

! Linearized values (Ostensibly Out)

      DOUBLE PRECISION, INTENT (INOUT) :: L_UHOM_DNDN ( MAX_USER_STREAMS, &
          MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: L_UHOM_DNUP ( MAX_USER_STREAMS, &
          MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: L_UHOM_UPDN ( MAX_USER_STREAMS, &
          MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: L_UHOM_UPUP ( MAX_USER_STREAMS, &
          MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: L_ZETA_M ( MAXEVALUES, &
          MAX_USER_STREAMS, MAXLAYERS,  MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: L_ZETA_P ( MAXEVALUES, &
          MAX_USER_STREAMS, MAXLAYERS,  MAX_ATMOSWFS )

!  Local variables
!  ---------------

      DOUBLE PRECISION :: L_HELPSTOKES_R(MAXSTOKES)
      DOUBLE PRECISION :: L_HELPSTOKES_CR(MAXSTOKES)
      DOUBLE PRECISION :: L_HELPSTOKES_CI(MAXSTOKES)
      DOUBLE PRECISION :: L_GAUX_R(0:MAXMOMENTS,MAXSTOKES)
      DOUBLE PRECISION :: L_GAUX_CR(0:MAXMOMENTS,MAXSTOKES)
      DOUBLE PRECISION :: L_GAUX_CI(0:MAXMOMENTS,MAXSTOKES)

      LOGICAL ::          DOLAYER_PP
      INTEGER ::          UM, J, L, N, M, Q, O1, O2
      INTEGER ::          K, KO1, K0, K1, K2
      DOUBLE PRECISION :: SPS, H, SP, SM, SG, SN, SNS, HPC, HMC
      DOUBLE PRECISION :: SPSCR, SPSCI, SPCR, SPCI
      DOUBLE PRECISION :: SNSCR, SNSCI, SNCR, SNCI
      DOUBLE PRECISION :: SR, SI, SMCR, SMCI, ZEPSQ, ZEMSQ
      DOUBLE PRECISION :: KISQ, RHO_P_CR, RHO_M_CR, MP, MM
      DOUBLE PRECISION :: L_KI, L_KR, L_KISQ, LMP, LMM

!  Layer and Fourier

      N = GIVEN_LAYER
      M = FOURIER

!  general flag

      DOLAYER_PP = ( STERM_LAYERMASK_UP(N) .OR. &
                     STERM_LAYERMASK_DN(N) )

!  Linearized Zeta constants
!  =========================

      IF ( DOLAYER_PP ) THEN

       DO UM = LOCAL_UM_START, N_USER_STREAMS
        SM = USER_SECANTS(UM)

!  Linearized Real Zeta constants (always required)

        DO K = 1, K_REAL(N)
          ZEPSQ = ZETA_P(K,UM,N) * ZETA_P(K,UM,N)
          ZEMSQ = ZETA_M(K,UM,N) * ZETA_M(K,UM,N)
          DO Q = 1, N_PARAMETERS
            L_ZETA_P(K,UM,N,Q) = - L_KEIGEN(K,N,Q) * ZEPSQ
            L_ZETA_M(K,UM,N,Q) = + L_KEIGEN(K,N,Q) * ZEMSQ
          ENDDO
        ENDDO

!  linearization of Complex Zeta constants

        IF ( K_COMPLEX(N) .GT. 0 ) THEN
          KO1 = K_REAL(N) + 1
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            KISQ = KEIGEN(K2,N)*KEIGEN(K2,N)
            RHO_P_CR = SM + KEIGEN(K1,N)
            RHO_M_CR = SM - KEIGEN(K1,N)
            MP = ONE / ( RHO_P_CR*RHO_P_CR + KISQ )
            MM = ONE / ( RHO_M_CR*RHO_M_CR + KISQ )
            DO Q = 1, N_PARAMETERS
             L_KISQ = L_KEIGEN(K2,N,Q)*KEIGEN(K2,N)
             L_KI   = L_KEIGEN(K2,N,Q)
             L_KR   = L_KEIGEN(K1,N,Q)
             LMP    = TWO * ( L_KISQ + RHO_P_CR * L_KR )
             LMM    = TWO * ( L_KISQ - RHO_M_CR * L_KR )
             L_ZETA_P(K1,UM,N,Q) = - MP * ( ZETA_P(K1,UM,N)*LMP - L_KR )
             L_ZETA_P(K2,UM,N,Q) = - MP * ( ZETA_P(K2,UM,N)*LMP + L_KI )
             L_ZETA_M(K1,UM,N,Q) = - MM * ( ZETA_M(K1,UM,N)*LMM + L_KR )
             L_ZETA_M(K2,UM,N,Q) = - MM * ( ZETA_M(K2,UM,N)*LMM - L_KI )
            ENDDO
          ENDDO
        ENDIF

!  End zeta linearization

       ENDDO
      ENDIF

!  Eigenvector interpolation to user-defined angles
!  ================================================

!  Only do this if both flags are set,
!    Or if there is a solution for scattering in this layer

      IF ( .NOT.DOLAYER_PP .OR. .NOT. DO_VARY .OR. &
           .NOT. DO_LAYER_SCATTERING(M,N) ) RETURN

!  For each parameter

      DO Q = 1, N_PARAMETERS

!  User defined solutions (Real)
!  -----------------------------

       DO K = 1, K_REAL(N)

!  For each moment, do inner sum over computational angles
!  for the positive and negative linearized eigenvectors.
!  Linearized the product with OMEGA_GREEK ---> L_GAUX

        DO L = M, NMOMENTS
         DO O1 = 1, NSTOKES
          SPS = ZERO
          DO J = 1, NSTREAMS
           H = QUAD_HALFWTS(J)
           SP = ZERO
           SM = ZERO
           DO O2 = 1, NSTOKES
! older     SP = SP + H*PI_XQP    (L,J,O1,O2)*L_FWD_XPOS_R(J,O2,K,N,Q)
! older     SM = SM + H*PI_XQM_PRE(L,J,O1,O2)*L_FWD_XNEG_R(J,O2,K,N,Q)
            SP = SP + H*PI_XQP    (L,J,O1,O2)*L_SOLA_XPOS(J,O2,K,N,Q)
            SM = SM + H*PI_XQM_PRE(L,J,O1,O2)*L_SOLB_XNEG(J,O2,K,N,Q)
           ENDDO
           SPS = SPS + SP + SM
          ENDDO
          L_HELPSTOKES_R(O1) = SPS
         ENDDO
         DO O1 = 1, NSTOKES
          SG = ZERO
          DO O2 = 1, NSTOKES
           SG = SG + OMEGA_GREEK(L,N,O1,O2)     * L_HELPSTOKES_R(O2) &
                   + L_OMEGA_GREEK(L,N,O1,O2,Q) *   HELPSTOKES(L,K,O2)
          ENDDO
          L_GAUX_R(L,O1) = SG
         ENDDO
        ENDDO

!  Now sum over all harmonic contributions (downwelling)

        IF ( DO_DNWELLING .AND. STERM_LAYERMASK_DN(N)) THEN
         DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
           SP = ZERO
           SN = ZERO
           DO L = M, NMOMENTS
            SPS = ZERO
            SNS = ZERO
            DO O2 = 1, NSTOKES
             SPS = SPS + PI_XUP(L,UM,O1,O2)      * L_GAUX_R(L,O2)
             SNS = SNS + PI_XUM_POST(L,UM,O1,O2) * L_GAUX_R(L,O2)
            ENDDO
            SP = SP + SPS
            SN = SN + SNS
           ENDDO
           L_UHOM_DNDN(UM,O1,K,N,Q) = SP
           L_UHOM_DNUP(UM,O1,K,N,Q) = SN
          ENDDO
         ENDDO
        ENDIF

!  Now sum over all harmonic contributions (Upwelling)

        IF ( DO_UPWELLING ) THEN
         DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
           SP = ZERO
           SN = ZERO
           DO L = M, NMOMENTS
            SPS = ZERO
            SNS = ZERO
            DO O2 = 1, NSTOKES
             SPS = SPS + PI_XUP_PRE(L,UM,O1,O2)  * L_GAUX_R(L,O2)
             SNS = SNS + PI_XUM(L,UM,O1,O2)      * L_GAUX_R(L,O2)
            ENDDO
            SP = SP + SPS
            SN = SN + SNS
           ENDDO
           L_UHOM_UPDN(UM,O1,K,N,Q) = SN
           L_UHOM_UPUP(UM,O1,K,N,Q) = SP
          ENDDO
         ENDDO
        ENDIF

!  Debug (important)

        IF ( DO_DEBUG_WRITE.AND.Q.EQ.2) THEN
         DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            WRITE(95,'(4I3,1P8E20.10)')M,K,UM,O1, &
               L_UHOM_UPUP(UM,O1,K,N,Q), &
               L_UHOM_UPDN(UM,O1,K,N,Q), &
               L_UHOM_DNUP(UM,O1,K,N,Q), &
               L_UHOM_DNDN(UM,O1,K,N,Q), &
               UHOM_UPUP(UM,O1,K,N), &
               UHOM_UPDN(UM,O1,K,N), &
               UHOM_DNUP(UM,O1,K,N), &
               UHOM_DNDN(UM,O1,K,N)
           ENDDO
         ENDDO
        ENDIF

!  end loop over real eigensolutions

       ENDDO

!  User defined solutions (Complex)
!  -------------------------------

!  offset

       KO1 = K_REAL(N) + 1

       DO K = 1, K_COMPLEX(N)

        K0 = 2 * K - 2
        K1 = KO1 + K0
        K2 = K1  + 1

!  For each moment, do inner sum over computational angles
!  for the positive and negative eigenvectors.
!    Gvec = STOKES. Eq 34. GAUx = Greekmat.Gvec

        DO L = M, NMOMENTS
         DO O1 = 1, NSTOKES
          SPSCR = ZERO
          SPSCI = ZERO
          DO J = 1, NSTREAMS
           SPCR = ZERO
           SPCI = ZERO
           SMCR = ZERO
           SMCI = ZERO
           DO O2 = 1, NSTOKES
            HPC  =  PI_XQP(L,J,O1,O2)*QUAD_HALFWTS(J)
! older     SPCR = SPCR + HPC * L_FWD_XPOS_CR(J,O2,K,N,Q)
! older     SPCI = SPCI + HPC * L_FWD_XPOS_CI(J,O2,K,N,Q)
            SPCR = SPCR + HPC * L_SOLA_XPOS(J,O2,K1,N,Q)
            SPCI = SPCI + HPC * L_SOLA_XPOS(J,O2,K2,N,Q)
            HMC  =  PI_XQM_PRE(L,J,O1,O2)*QUAD_HALFWTS(J)
! older     SMCR = SMCR + HMC * L_FWD_XNEG_CR(J,O2,K,N,Q)
! older     SMCI = SMCI + HMC * L_FWD_XNEG_CI(J,O2,K,N,Q)
            SMCR = SMCR + HMC * L_SOLB_XNEG(J,O2,K1,N,Q)
            SMCI = SMCI + HMC * L_SOLB_XNEG(J,O2,K2,N,Q)
           ENDDO
           SPSCR = SPSCR + SPCR + SMCR
           SPSCI = SPSCI + SPCI + SMCI
          ENDDO
          L_HELPSTOKES_CR(O1) = SPSCR
          L_HELPSTOKES_CI(O1) = SPSCI
         ENDDO
         DO O1 = 1, NSTOKES
          SR = ZERO
          SI = ZERO
          DO O2 = 1, NSTOKES
           SR = SR + OMEGA_GREEK(L,N,O1,O2)     *L_HELPSTOKES_CR(O2) &
                   + L_OMEGA_GREEK(L,N,O1,O2,Q) *  HELPSTOKES(L,K1,O2)
           SI = SI + OMEGA_GREEK(L,N,O1,O2)     *L_HELPSTOKES_CI(O2) &
                   + L_OMEGA_GREEK(L,N,O1,O2,Q) *  HELPSTOKES(L,K2,O2)
          ENDDO
          L_GAUX_CR(L,O1) = SR
          L_GAUX_CI(L,O1) = SI
         ENDDO
        ENDDO

!  Now sum over all harmonic contributions (downwelling)

        IF ( DO_DNWELLING ) THEN
         DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
           SPCR = ZERO
           SNCR = ZERO
           SPCI = ZERO
           SNCI = ZERO
           DO L = M, NMOMENTS
            SPSCR = ZERO
            SNSCR = ZERO
            SPSCI = ZERO
            SNSCI = ZERO
            DO O2 = 1, NSTOKES
             HPC = PI_XUP(L,UM,O1,O2)
             SPSCR = SPSCR + HPC * L_GAUX_CR(L,O2)
             SPSCI = SPSCI + HPC * L_GAUX_CI(L,O2)
             HMC = PI_XUM_POST(L,UM,O1,O2)
             SNSCR = SNSCR + HMC * L_GAUX_CR(L,O2)
             SNSCI = SNSCI + HMC * L_GAUX_CI(L,O2)
            ENDDO
            SPCR = SPCR + SPSCR
            SPCI = SPCI + SPSCI
            SNCR = SNCR + SNSCR
            SNCI = SNCI + SNSCI
           ENDDO
           L_UHOM_DNDN(UM,O1,K1,N,Q) = SPCR
           L_UHOM_DNUP(UM,O1,K1,N,Q) = SNCR
           L_UHOM_DNDN(UM,O1,K2,N,Q) = SPCI
           L_UHOM_DNUP(UM,O1,K2,N,Q) = SNCI
          ENDDO
         ENDDO
        ENDIF

!  Now sum over all harmonic contributions (upwelling)

        IF ( DO_UPWELLING ) THEN
         DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
           SPCR = ZERO
           SNCR = ZERO
           SPCI = ZERO
           SNCI = ZERO
           DO L = M, NMOMENTS
            SPSCR = ZERO
            SNSCR = ZERO
            SPSCI = ZERO
            SNSCI = ZERO
            DO O2 = 1, NSTOKES
             HPC = PI_XUP_PRE(L,UM,O1,O2)
             SPSCR = SPSCR + HPC * L_GAUX_CR(L,O2)
             SPSCI = SPSCI + HPC * L_GAUX_CI(L,O2)
             HMC = PI_XUM(L,UM,O1,O2)
             SNSCR = SNSCR + HMC * L_GAUX_CR(L,O2)
             SNSCI = SNSCI + HMC * L_GAUX_CI(L,O2)
            ENDDO
            SPCR = SPCR + SPSCR
            SPCI = SPCI + SPSCI
            SNCR = SNCR + SNSCR
            SNCI = SNCI + SNSCI
           ENDDO
           L_UHOM_UPDN(UM,O1,K1,N,Q) = SNCR
           L_UHOM_UPUP(UM,O1,K1,N,Q) = SPCR
           L_UHOM_UPDN(UM,O1,K2,N,Q) = SNCI
           L_UHOM_UPUP(UM,O1,K2,N,Q) = SPCI
          ENDDO
         ENDDO
        ENDIF

!  Debug (important)

        IF ( DO_DEBUG_WRITE.AND.Q.EQ.2) THEN
         DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            WRITE(97,'(4I3,1P8E20.10)')M,K,UM,O1, &
               L_UHOM_UPUP(UM,O1,K2,N,Q), &
               L_UHOM_UPDN(UM,O1,K2,N,Q), &
               L_UHOM_DNUP(UM,O1,K2,N,Q), &
               L_UHOM_DNDN(UM,O1,K2,N,Q), &
               UHOM_UPUP(UM,O1,K2,N), &
               UHOM_UPDN(UM,O1,K2,N), &
               UHOM_DNUP(UM,O1,K2,N), &
               UHOM_DNDN(UM,O1,K2,N)
           ENDDO
         ENDDO
        ENDIF

!  end loop over complex eigensolutions

       ENDDO

!  End loop over parameters

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_L_UHOM_SOLUTION

!

      SUBROUTINE L_HMULT_MASTER ( &
        DO_UPWELLING, DO_DNWELLING, TAYLOR_ORDER, & ! 2p7 Taylor-order argument added
        NLAYERS, N_USER_LEVELS, &
        N_USER_STREAMS, LOCAL_UM_START, &
        USER_SECANTS, &
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
        PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, &
        STERM_LAYERMASK_DN, &
        DELTAU_VERT, PARTAU_VERT, T_DELT_EIGEN, & ! @@@ Rob Fix 5/10/13 added PARTAU_VERT
        T_UTUP_EIGEN, T_UTDN_EIGEN, &
        T_DELT_USERM, T_UTDN_USERM, &
        T_UTUP_USERM, &
        K_REAL, K_COMPLEX, HSINGO, &
        HMULT_1, HMULT_2, ZETA_M, ZETA_P, &
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
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
      INTEGER, INTENT (IN) ::          TAYLOR_ORDER  ! 2p7 Taylor-order argument added
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
      DOUBLE PRECISION, INTENT (IN) :: PARTAU_VERT ( MAX_PARTLAYERS )  ! @@@ Rob Fix 5/10/13 added
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
      DOUBLE PRECISION :: L_T2, L_T1, HOM1, HOM2, SM, EPS, DELTA, L_DELTA, L_MULT
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

      DO N = 1, NLAYERS
        DO_RTSOL_VARY(N) = LAYER_VARY_FLAG(N)
        NPARAMS_VARY(N)  = LAYER_VARY_NUMBER(N)
      ENDDO

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
                  EPS      = ZETA_M(K,UM,N) ; DELTA   = DELTAU_VERT(N)
                  L_DELTA  = L_DELTAU_VERT(Q,N) * DELTA ! Input is single normalized
                  CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, &
                                         L_KEIGEN(K,N,Q), ZERO, UDEL, SM, L_MULT )
                  L_HMULT_1(K,UM,N,Q) = SM * L_MULT
!  10 May 13       CALL TAYLOR_SERIES_L_1 ( EPS, DELTAU_VERT(N), SM, ONE, UDEL, &
!  10 May 13                                L_KEIGEN(K,N,Q), L_DEL, L_HMULT_1(K,UM,N,Q) )
!   Old            HOM1 = L_DELTAU_VERT(Q,N)*(ONE-SM) - &
!   Old                  HALF * L_KEIGEN(K,N,Q) * DELTAU_VERT(N)
!   Old            L_HMULT_1(K,UM,N,Q) = HMULT_1(K,UM,N) * HOM1
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

            IF ( HSINGO(K,UM,N) )THEN
              EPS   = ZETA_M(K,UM,N) ; DELTA   = DELTAU_VERT(N) - PARTAU_VERT(UT)
              L_DELTA  = L_DELTAU_VERT(Q,N) * DELTA
              CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, &
                                       L_KEIGEN(K,N,Q), ZERO, UX_UP, SM, L_MULT )
              L_UT_HMULT_UU(K,UM,UT,Q)= SM * L_MULT
            ELSE
              L_UT_HMULT_UU(K,UM,UT,Q) = SM * &
               ( L_THETA_UP_R *   ZETA_M(K,UM,N) + &
                   THETA_UP_R * L_ZETA_M(K,UM,N,Q) )
            ENDIF

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

            IF ( HSINGO(K,UM,N) )THEN
              EPS      = ZETA_M(K,UM,N) ; DELTA   = PARTAU_VERT(UT)
              L_DELTA  = L_DELTAU_VERT(Q,N) * DELTA  ! Input is Single normalized (partau)
              CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, &
                                       L_KEIGEN(K,N,Q), ZERO, UX_DN, SM, L_MULT )
              L_UT_HMULT_DD(K,UM,UT,Q) = SM * L_MULT
            ELSE
              L_UT_HMULT_DD(K,UM,UT,Q) = SM * &
               ( L_THETA_DN_R *   ZETA_M(K,UM,N) + &
                   THETA_DN_R * L_ZETA_M(K,UM,N,Q) )
            ENDIF

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

! @@@ Rob Fix 5/10/13 - This is a good place to pause for debug
!      pause'End L_HMULT_MASTER'

      RETURN
      END SUBROUTINE L_HMULT_MASTER

      END MODULE vlidort_lpc_solutions

