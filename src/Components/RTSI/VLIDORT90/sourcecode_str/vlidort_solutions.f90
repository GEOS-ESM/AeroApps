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
! #              VLIDORT_QHOM_SOLUTION                          #
! #              VLIDORT_UHOM_SOLUTION                          #
! #              VLIDORT_QBEAM_SOLUTION                         #
! #              VLIDORT_UBEAM_SOLUTION                         #
! #              VLIDORT_PIMATRIX_SETUP                         #
! #                                                             #
! ###############################################################


      MODULE vlidort_solutions

      PRIVATE
      PUBLIC :: VLIDORT_QHOM_SOLUTION, &
                VLIDORT_UHOM_SOLUTION, &
                VLIDORT_QBEAM_SOLUTION, &
                VLIDORT_UBEAM_SOLUTION, &
                VLIDORT_PIMATRIX_SETUP

      CONTAINS

      SUBROUTINE VLIDORT_QHOM_SOLUTION ( &
        GIVEN_LAYER, FOURIER, &
        DO_SOLUTION_SAVING, NSTOKES, &
        NSTREAMS, N_USER_LEVELS, &
        QUAD_STREAMS, QUAD_HALFWTS, &
        NMOMENTS, NSTKS_NSTRMS, &
        DO_REAL_EIGENSOLVER, DO_LAYER_SCATTERING, &
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
        PARTLAYERS_LAYERIDX, &
        DELTAU_VERT, PARTAU_VERT, &
        OMEGA_GREEK, T_DELT_DISORDS, &
        T_DISORDS_UTUP, T_DISORDS_UTDN, &
        PI_XQP, PI_XQM, &
        PI_XQM_PRE, &
        SAB, DAB, &
        EIGENMAT_SAVE, REAL_KSQ, &
        IMAG_KSQ, LEFT_EVEC, &
        RITE_EVEC, EIGENDEGEN, &
        EIGENMASK_R, EIGENMASK_C, &
        FWD_SUMVEC, FWD_DIFVEC, &
        T_DELT_EIGEN, T_UTUP_EIGEN, &
        T_UTDN_EIGEN, &
        K_REAL, K_COMPLEX, &
        KEIGEN, KEIGEN_CSQ, &
        SOLA_XPOS, SOLB_XNEG, &
        STATUS, MESSAGE, TRACE )

!  Numerical solution of Eigenproblem.

      USE VLIDORT_PARS
      USE VLIDORT_AUX
      USE LAPACK_TOOLS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::           GIVEN_LAYER
      INTEGER, INTENT (IN) ::           FOURIER
      LOGICAL, INTENT (IN) ::           DO_SOLUTION_SAVING
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           N_USER_LEVELS
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_HALFWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::           NMOMENTS
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS
      LOGICAL, INTENT (IN) ::           DO_REAL_EIGENSOLVER &
          ( 0:MAXMOMENTS, MAXLAYERS )
      LOGICAL, INTENT (IN) ::           DO_LAYER_SCATTERING &
          ( 0:MAXMOMENTS, MAXLAYERS )
      LOGICAL, INTENT (IN) ::           PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  PARTAU_VERT ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  OMEGA_GREEK &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DISORDS_UTUP &
          ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DISORDS_UTDN &
          ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XQP &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XQM &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XQM_PRE &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )

      DOUBLE PRECISION, INTENT (INOUT) ::  SAB &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) ::  DAB &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) ::  EIGENMAT_SAVE &
          ( MAXEVALUES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) ::  REAL_KSQ ( MAXSTRMSTKS )
      DOUBLE PRECISION, INTENT (INOUT) ::  IMAG_KSQ ( MAXSTRMSTKS )
      DOUBLE PRECISION, INTENT (INOUT) ::  LEFT_EVEC &
          ( MAXSTRMSTKS, MAXSTRMSTKS )
      DOUBLE PRECISION, INTENT (INOUT) ::  RITE_EVEC &
          ( MAXSTRMSTKS, MAXSTRMSTKS )
      LOGICAL, INTENT (INOUT) ::           EIGENDEGEN ( MAXSTRMSTKS, MAXLAYERS )
      INTEGER, INTENT (INOUT) ::           EIGENMASK_R ( MAXEVALUES )
      INTEGER, INTENT (INOUT) ::           EIGENMASK_C ( MAXEVALUES )
      DOUBLE PRECISION, INTENT (INOUT) ::  FWD_SUMVEC &
          ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION, INTENT (INOUT) ::  FWD_DIFVEC &
          ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION, INTENT (INOUT) ::  T_DELT_EIGEN &
          ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) ::  T_UTUP_EIGEN &
          ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) ::  T_UTDN_EIGEN &
          ( MAXEVALUES, MAX_PARTLAYERS )
      INTEGER, INTENT (INOUT) ::           K_REAL ( MAXLAYERS )
      INTEGER, INTENT (INOUT) ::           K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) ::  KEIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) ::  KEIGEN_CSQ ( MAXEVALUES )
      DOUBLE PRECISION, INTENT (INOUT) ::  SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) ::  SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

      INTEGER, INTENT (OUT) ::           STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  Local variables
!  ---------------

!  local matrix for eigenvalue computation

      DOUBLE PRECISION :: EIGENMAT ( MAXSTRMSTKS, MAXSTRMSTKS )

!  output from Eigenpackage module ASYMTX
!    Following output now stored in commons:
!      DOUBLE PRECISION REAL_KSQ(MAXSTRMSTKS)
!      DOUBLE PRECISION RITE_EVEC(MAXSTRMSTKS,MAXSTRMSTKS)

!  Additional diagnostic output from ASYMTX, and Work space.
!    Dimension of scratch array WK increased to 4 times the basic dimens
!    R. Spurr, RTS, December 2005: More work space needed to avoid seg-f
!    Internal dimensioning in ASYMTX is the same.

!     DOUBLE PRECISION ::  WK(16*MAXSTRMSTKS)
      DOUBLE PRECISION ::  WK(4*MAXSTRMSTKS)
      INTEGER ::           IER
      LOGICAL ::           ASYMTX_FAILURE
      CHARACTER (LEN=3) :: CN, CI

!  output for Eigenpackage DGEEV ( LAPACK)
!    Following now stored in commons...............
!      DOUBLE PRECISION REAL_KSQ(MAXSTRMSTKS),IMAG_KSQ(MAXSTRMSTKS)
!      DOUBLE PRECISION LEFT_EVEC(MAXSTRMSTKS,MAXSTRMSTKS)
!      DOUBLE PRECISION RITE_EVEC(MAXSTRMSTKS,MAXSTRMSTKS)

!  Additional diagnostic output and Work space from DGEEV.
!    Dimension LWORK increased from 4 to 64 times the basic dimension
!    J. Kujanpaa, FMI, August 2005: More work space for DGEEV needed.

      INTEGER, PARAMETER :: LWORK = 64*MAXSTRMSTKS
      DOUBLE PRECISION ::   WORK ( LWORK )
      INTEGER ::            DGEEV_INFO

!  optimization of Sum and Difference matrices SAB and DAB
!  ...requires intermediate storage of these two arrays.
!   Timing tests by J. Kujanpaa, FMI, August 2005.

      DOUBLE PRECISION :: &
          H2PARR ( MAXSTOKES, 0:MAXMOMENTS, MAXSTOKES, MAXSTREAMS ), &
          H2MARR ( MAXSTOKES, 0:MAXMOMENTS, MAXSTOKES, MAXSTREAMS )

!  Miscellaneous local variables

      INTEGER ::          I, J, I1, L, N, NA, UTA, UT, NMINST
      INTEGER ::          M, AA, AA1, K, K_R, K_C, K_AA, KO1, K0, K1, K2
      INTEGER ::          O1, O2, O3, O4, IR, IROW, JC, JCOL
      DOUBLE PRECISION :: DP, DM, FAC, SC_R, XINV, HELP
      DOUBLE PRECISION :: H1P, H1M, NORM_R, TAU_DN, TAU_UP
      DOUBLE PRECISION :: FWD_H1, FWD_H2
      DOUBLE PRECISION :: XSQ_R, YSQ_R, RKSQ, IKSQ, KMUT, KVAL
      DOUBLE PRECISION :: CARG, SARG, ARG, MOD_KVAL, H1, H2
      DOUBLE PRECISION :: FWD_R, FWD_I, H2P, H2M, FPD, FND

!  original complex variables were:
!      COMPLEX*16      FWD_H1_C, FWD_H2_C, ADJ_H1_C, ADJ_H2_C
!      COMPLEX*16      HELP_C, SC_C, H1_C, H2_C

!  replacement for complex variables

      DOUBLE PRECISION :: HELP_CR, HELP_CI
      DOUBLE PRECISION :: FWD_H1_CR, FWD_H1_CI
      DOUBLE PRECISION :: FWD_H2_CR, FWD_H2_CI
      DOUBLE PRECISION :: SC_CR, SC_CI, MODKRT
      DOUBLE PRECISION :: FPD1, FPD2, FND1, FND2

!  Code start
!  ----------

!  initialise exception handling

      STATUS  = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Initialization

      KO1 = 0

!  initialise status

      STATUS = VLIDORT_SUCCESS

!  Layer and Fourier

      N = GIVEN_LAYER
      M = FOURIER

!  Flipper number

      NMINST = MIN(2,NSTOKES)

!  If there is no scattering in this layer, we can set
!  the solution vectors directly and move on.
!    Just a transmitttance case.

      IF ( DO_SOLUTION_SAVING ) THEN
       IF ( .NOT. DO_LAYER_SCATTERING(M,N) ) THEN
        K_COMPLEX(N) = 0
        K_REAL(N)    = NSTKS_NSTRMS
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          XINV = ONE / QUAD_STREAMS(I)
          DO O1 = 1, NSTOKES
            DO K = 1, K_REAL(N)
              SOLA_XPOS(I,O1,K,N)  = ZERO
              SOLA_XPOS(I1,O1,K,N) = ZERO
              SOLB_XNEG(I1,O1,K,N) = ZERO
              SOLB_XNEG(I,O1,K,N)  = ZERO
            ENDDO
          ENDDO
          DO O1 = 1, NMINST
            K = NSTOKES*(I-1) + O1
            KEIGEN(K,N) = XINV
            SOLA_XPOS(I,O1,K,N)  = ONE
            SOLB_XNEG(I1,O1,K,N) = ONE
          ENDDO
          DO O1 = 3, NSTOKES
            K = NSTOKES*(I-1) + O1
            KEIGEN(K,N) = XINV
            SOLA_XPOS(I,O1,K,N)  = MINUS_ONE
            SOLB_XNEG(I1,O1,K,N) = MINUS_ONE
          ENDDO
        ENDDO
        GO TO 3456
       ENDIF
      ENDIF

!  Scattering solutions
!  ====================

!  Construct Eigenmatrix
!  ---------------------

!  zero the eigensolution matrix
!   Not needed. Redundancy test by J. Kujanpaa, FMI, August 2005.
!      DO I = 1, NSTKS_NSTRMS
!        DO J = 1, NSTKS_NSTRMS
!          EIGENMAT(I,J) = ZERO
!        ENDDO
!      ENDDO

!  Develop Sum and Difference matrices, part I
!   Introduction of holding arrays by J. Kujanpaa, FMI, August 2005

      DO J = 1, NSTREAMS
       DO O2 = 1, NSTOKES
        DP = ZERO
        DM = ZERO
         DO L = M, NMOMENTS
          DO O3 = 1, NSTOKES
           H2P = ZERO
           H2M = ZERO
           DO O4 = 1, NSTOKES
            H2P = H2P + OMEGA_GREEK(L,N,O3,O4)*PI_XQP(L,J,O4,O2)
            H2M = H2M + OMEGA_GREEK(L,N,O3,O4)*PI_XQM_PRE(L,J,O4,O2)
           ENDDO
           H2PARR(O3,L,O2,J) = H2P
           H2MARR(O3,L,O2,J) = H2M
          ENDDO
         ENDDO
        ENDDO
      ENDDO

!  debug code for hparr, hmarr

!      DO J = 1, NSTREAMS
!       DO L = M, NMOMENTS
!         DO O2 = 1, min(NSTOKES,3)
!          DO O3 = 1, min(NSTOKES,3)
!           if (nstokes.eq.3.and.n.eq.6.and.m.eq.0)
!     *         write(87,'(4i3,1p2e17.8)')
!     *        J,L,O2,O3,H2PARR(O3,L,O2,J),H2MARR(O3,L,O2,J)
!           if (nstokes.eq.4.and.n.eq.6.and.m.eq.0)
!     *         write(88,'(4i3,1p2e17.8)')
!     *        J,L,O2,O3,H2PARR(O3,L,O2,J),H2MARR(O3,L,O2,J)
!          ENDDO
!         ENDDO
!        ENDDO
!      ENDDO
!      if (n.eq.6.and.m.eq.1)pause 'hp 6'

!  Develop Sum and Difference matrices, part II of calculation
!   Use of holding arrays by J. Kujanpaa, FMI, August 2005

      DO I = 1, NSTREAMS
       XINV = - ONE / QUAD_STREAMS(I)
       DO J = 1, NSTREAMS
        FAC = XINV * QUAD_HALFWTS(J)
        DO O1 = 1, NSTOKES
         DO O2 = 1, NSTOKES
          DP = ZERO
          DM = ZERO
          DO L = M, NMOMENTS
           H1P = ZERO
           H1M = ZERO
           DO O3 = 1, NSTOKES
            H1P = H1P + PI_XQP(L,I,O1,O3)*H2PARR(O3,L,O2,J)
            H1M = H1M + PI_XQP(L,I,O1,O3)*H2MARR(O3,L,O2,J)
           ENDDO
           DP = DP + H1P
           DM = DM + H1M
          ENDDO
          SAB(I,J,O1,O2,N) = FAC * ( DP + DM )
          DAB(I,J,O1,O2,N) = FAC * ( DP - DM )
         ENDDO
        ENDDO
       ENDDO
       DO O1 = 1, NSTOKES
        SAB(I,I,O1,O1,N) = SAB(I,I,O1,O1,N) - XINV
        DAB(I,I,O1,O1,N) = DAB(I,I,O1,O1,N) - XINV
       ENDDO
      ENDDO

!  Note. In some circumstances, H2PARR = H2MARR
!        Then DP = Dm in the above code, and DAB is Diagonal.
!        Then there are pair-degenerate eigenvalues, equal to the
!        discrete ordinate directions.

!  Compute Eigenmatrix

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
            H2 = H2 + DAB(I,K,O1,O3,N) * SAB(K,J,O3,O2,N)
           ENDDO
           H1 = H1 + H2
          ENDDO
          EIGENMAT(IROW,JCOL) = H1
         ENDDO
        ENDDO
       ENDDO
      ENDDO

!  save Eigenmatrix (original is destroyed, want to use later)
!    Looping order changed by J. Kujanpaa, FMI, August 2005

      DO J = 1, NSTKS_NSTRMS
        DO I = 1, NSTKS_NSTRMS
          EIGENMAT_SAVE(I,J,N) = EIGENMAT(I,J)
        ENDDO
      ENDDO

!      if (n.eq.24) then
!      do j = 1, 9
!      if(nstokes.eq.3)write(33,'(i4,1p12e11.3)')J,(EIGENMAT(i,J),i=1,9)
!      if(nstokes.eq.4)write(34,'(i4,1p12e11.3)')J,(EIGENMAT(i,J),i=1,12
!      enddo
!      pause
!      endif

!  Debug for testing new eigensolver
!  7 August 2006. Result --> NEW LAPACK ROUTINES NOT FASTER

!      IF ( N.EQ.1) then
!        OPEN(88,file='esolver.inp',status='unknown')
!      ENDIF
!      write(88,'(2i5,l2)')n,NSTKS_NSTRMS,DO_REAL_EIGENSOLVER(M,N)
!      DO I = 1, NSTKS_NSTRMS
!        write(88,'(1p32e20.10)')(eigenmat(i,j),j=1,NSTKS_NSTRMS)
!      enddo
!      if ( n.eq.nlayers)then
!        close(88)
!        pause'end debug'
!      ENDIF

!  Eigensolver package
!  -------------------

      IF ( DO_REAL_EIGENSOLVER(M,N) ) THEN

!  Here is the DISORT ASYMTX package, as used in LIDORT scalar codes
!    This is sufficient for Real Symmetric Eigenmatrices.

       CALL  ASYMTX &
            ( EIGENMAT, NSTKS_NSTRMS, MAXSTRMSTKS, MAXSTRMSTKS, &
              RITE_EVEC, REAL_KSQ, IER, WK, &
              MESSAGE, ASYMTX_FAILURE )

!  Exception handling 1

       IF ( ASYMTX_FAILURE  ) THEN
         WRITE(CN,'(I3)')N
         TRACE   ='ASYMTX error in VLIDORT_QHOM_SOLUTION, Layer='//CN
         STATUS  = VLIDORT_SERIOUS
         RETURN
       ENDIF

!  Exception handling 2

       IF ( IER.GT.0 ) THEN
         WRITE(CI,'(I3)')IER
         WRITE(CN,'(I3)')N
         MESSAGE = 'eigenvalue '//CI//' has not converged'
         TRACE   = 'ASYMTX error in VLIDORT_QHOM_SOLUTION, Layer ='//CN
         STATUS  = VLIDORT_SERIOUS
         RETURN
       ENDIF

!  renormalize - this is vital for ASMTX output

       DO AA = 1, NSTKS_NSTRMS
        NORM_R = ZERO
        DO I = 1, NSTKS_NSTRMS
          NORM_R = NORM_R + RITE_EVEC(I,AA)*RITE_EVEC(I,AA)
        ENDDO
        NORM_R = DSQRT(NORM_R)
        DO I = 1, NSTKS_NSTRMS
          RITE_EVEC(I,AA) = RITE_EVEC(I,AA) / NORM_R
        ENDDO
       ENDDO

      ELSE

!  Complex roots : Must use DGEEV (LAPACK module)
!  DGEEV is 1.7-2.0 times slower than ASYMTX because look for complex ro
!  Also first call to DGEEV is 20 times slower.

       CALL DGEEV &
            ( 'V', 'V', NSTKS_NSTRMS, EIGENMAT, &
              MAXSTRMSTKS, REAL_KSQ, IMAG_KSQ, &
              LEFT_EVEC, MAXSTRMSTKS, RITE_EVEC, MAXSTRMSTKS, &
              WORK, LWORK, DGEEV_INFO )

!  Exception handling

       IF ( DGEEV_INFO .NE. 0 ) THEN
        WRITE(CN,'(I3)')N
        TRACE = 'DGEEV call, layer '//CN//' in VLIDORT_QHOM_SOLUTION'
        IF ( DGEEV_INFO .LT. 0 ) THEN
         WRITE(CI,'(I3)')IABS(DGEEV_INFO)
         MESSAGE='Argument # '//CI//' had an illegal value'
        ELSE
         MESSAGE='QR algorithm in DGEEV failed to compute all e-values'
        ENDIF
        STATUS = VLIDORT_SERIOUS
        RETURN
       ENDIF

      ENDIF

!  Sort real and complex eigenvalues; get masks
!  --------------------------------------------

!  Initialise counts

       K_R  = 0
       K_C  = 0
       K_AA = 0

!  all eigenvalues are real if ASYMTX is used

      IF ( DO_REAL_EIGENSOLVER(M,N) ) THEN
       DO AA = 1, NSTKS_NSTRMS
        K_R = K_R + 1
        EIGENMASK_R(K_R) = AA
       ENDDO
      ENDIF

!  Degeneracy mask

      IF ( DO_REAL_EIGENSOLVER(M,N) ) THEN
       DO AA = 1, NSTKS_NSTRMS
         KVAL = DSQRT(REAL_KSQ(AA))
         EIGENDEGEN(AA,N) = .FALSE.
         DO K = 1, NSTKS_NSTRMS
           IF (K.NE.AA.AND..NOT.EIGENDEGEN(AA,N)) THEN
             KMUT= DSQRT(REAL_KSQ(K))
             EIGENDEGEN(AA,N) = (DABS(KMUT-KVAL).LT.1.0D-04)
           ENDIF
         ENDDO
       ENDDO
      ENDIF

!  For the DGEEV eigensolver--
!    real eigenvalues if IMAG_KSQ = 0, otherwise complex in pairs

      IF ( .NOT.DO_REAL_EIGENSOLVER(M,N) ) THEN
       DO AA = 1, NSTKS_NSTRMS
        IF ( IMAG_KSQ(AA) .EQ. ZERO ) THEN
          K_R = K_R + 1
          EIGENMASK_R(K_R) = AA
        ELSE
          K_AA = K_AA + 1
          IF ( K_AA .EQ. 1 ) THEN
            K_C = K_C + 1
            EIGENMASK_C(K_C) = AA
          ELSE IF ( K_AA .EQ. 2 ) THEN
            K_AA = 0
          ENDIF
        ENDIF
       ENDDO
      ENDIF

!  save number of eigenvalues

      K_REAL(N)    = K_R
      K_COMPLEX(N) = K_C

!  Real Solutions
!  --------------

!  start loop over real eigenvalues

      DO K = 1, K_REAL(N)

!  get eigenvector entry

        AA = EIGENMASK_R(K)

!  eigenvalue and separation constants (save the eigenvalues)

        KEIGEN(K,N) = DSQRT(REAL_KSQ(AA))
        SC_R    = ONE / KEIGEN(K,N)

!  set Right eigenvector norms
!    This piece of code should not be necessary as
!    (a) DGEEV norms are all = 1 upon output
!    (b) ASMTYX output has been renormalized to unity

        NORM_R = ZERO
        DO I = 1, NSTKS_NSTRMS
          NORM_R = NORM_R + RITE_EVEC(I,AA)*RITE_EVEC(I,AA)
        ENDDO
        NORM_R = DSQRT(NORM_R)

!  Forward Sum-vector = normalized Right eigenvector

        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            FWD_SUMVEC(I,O1,K) = RITE_EVEC(IROW,AA)/NORM_R
          ENDDO
        ENDDO

!  Debug (important)

!        if ( do_debug_write ) then
!         if (m.eq.1.and.n.gt.4) then
!          write(27,'(2I3,1pe20.12)') N,K,KEIGEN(K,N)
!          DO I = 1, NSTREAMS
!           DO O1 = 1, NSTOKES
!            write(27,'(4i3,1pe20.12)')N,k,i,o1,FWD_SUMVEC(I,O1,K)
!           ENDDO
!          ENDDO
!         ENDIF
!        endif

!  Find Forward difference vectors (Siewert's notation)

        DO I = 1, NSTREAMS
         DO O1 = 1, NSTOKES
          FWD_H1 = ZERO
          DO J = 1, NSTREAMS
           FWD_H2 = ZERO
           DO O2 = 1, NSTOKES
            FWD_H2 = FWD_H2 + SAB (I,J,O1,O2,N)*FWD_SUMVEC(J,O2,K)
           ENDDO
           FWD_H1 = FWD_H1 + FWD_H2
          ENDDO
          FWD_DIFVEC(I,O1,K) = FWD_H1 * SC_R
         ENDDO
        ENDDO

!  assign Forward solutions PHI (Siewert's notation)
!    --->  first N are "DOWN", last N are "UP" (streams)
!    --->  Use symmetry properties to set -ve eigensolutions

        DO I = 1, NSTREAMS
          XINV = HALF
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES
            FPD = XINV * ( FWD_SUMVEC(I,O1,K) + FWD_DIFVEC(I,O1,K) )
            FND = XINV * ( FWD_SUMVEC(I,O1,K) - FWD_DIFVEC(I,O1,K) )
            SOLA_XPOS(I,O1,K,N)  = FPD
            SOLB_XNEG(I,O1,K,N)  = FND
            SOLA_XPOS(I1,O1,K,N) = FND
            SOLB_XNEG(I1,O1,K,N) = FPD
            IF ( O1 .GT. 2 ) THEN
              SOLA_XPOS(I1,O1,K,N)  = - FND
              SOLB_XNEG(I1,O1,K,N)  = - FPD
            ENDIF
          ENDDO
        ENDDO

!  Older code: Use wasteful additional arrays...
!    solution assignation (only for BVP)

!        DO I = 1, NSTREAMS
!          XINV = HALF
!          I1 = I + NSTREAMS
!          FWD_XPOS(I,O1,K,N)  =
!     &      XINV * ( FWD_SUMVEC(I,O1,K) + FWD_DIFVEC(I,O1,K) )
!          FWD_XPOS(I1,O1,K,N)  =
!     &      XINV * ( FWD_SUMVEC(I,O1,K) - FWD_DIFVEC(I,O1,K) )
!          FWD_XNEG(I1,O1,K,N) = FWD_XPOS(I,O1,K,N)
!          FWD_XNEG(I,O1,K,N)  = FWD_XPOS(I1,O1,K,N)
!         ENDDO
!        ENDDO
!        DO I = 1, NSTREAMS
!          I1 = I + NSTREAMS
!          DO O1 = 1, NSTOKES
!            SOLA_XPOS(I,O1,K,N)  = FWD_XPOS(I,O1,K,N)
!            SOLB_XNEG(I,O1,K,N)  = FWD_XNEG(I,O1,K,N)
!          ENDDO
!          DO O1 = 1, NMINST
!            SOLA_XPOS(I1,O1,K,N)  = FWD_XNEG(I,O1,K,N)
!            SOLB_XNEG(I1,O1,K,N)  = FWD_XPOS(I,O1,K,N)
!          ENDDO
!          DO O1 = 3, NSTOKES
!            SOLA_XPOS(I1,O1,K,N)  = - FWD_XNEG(I,O1,K,N)
!            SOLB_XNEG(I1,O1,K,N)  = - FWD_XPOS(I,O1,K,N)
!          ENDDO
!        ENDDO

!  end loop over real eigenvalues

      ENDDO

!  Complex solutions
!  -----------------

!  Offsets

      KO1 = K_REAL(N) + 1

!  start loop over complex eigenvalues

      DO K = 1, K_COMPLEX(N)

!  Bookkeeping indices

        K0 = 2 * ( K - 1 )
        K1 = KO1 + K0
        K2 = K1  + 1

!  get eigenvector entry

        AA  = EIGENMASK_C(K)
        AA1 = AA + 1

!  set eigenvalue and separation constants for one of pair of Conjugates
!   ---> (the complex conjugate follows from this information)

        RKSQ = REAL_KSQ(AA) * REAL_KSQ(AA)
        IKSQ = IMAG_KSQ(AA) * IMAG_KSQ(AA)
        MOD_KVAL = DSQRT(RKSQ + IKSQ)
        ARG = HALF * DATAN ( IMAG_KSQ(AA) / REAL_KSQ(AA) )
        CARG = DCOS(ARG)
        SARG = DSIN(ARG)

!  original complex --------------------------------------------------
!        KEIGEN_C(K,N) = DSQRT(MOD_KVAL) * CMPLX ( CARG, SARG )
!        SC_C  = ONE_C / KEIGEN_C(K,N)
!  -------------------------------------------------------------------

!  replacement complex

        KEIGEN_CSQ(K) = MOD_KVAL
        MODKRT        = DSQRT(MOD_KVAL)
        KEIGEN(K1,N)  = MODKRT * CARG
        KEIGEN(K2,N)  = MODKRT * SARG
        SC_CR =   KEIGEN(K1,N) / MOD_KVAL
        SC_CI = - KEIGEN(K2,N) / MOD_KVAL

!  set Right eigenvector norms
!   ---> Follows convention of way that eigenvectors are fomatted
!        upon output from the LAPACK module DGEEV.

        NORM_R = ZERO
        DO I = 1, NSTKS_NSTRMS
          XSQ_R = RITE_EVEC(I,AA)  * RITE_EVEC(I,AA)
          YSQ_R = RITE_EVEC(I,AA1) * RITE_EVEC(I,AA1)
          NORM_R = NORM_R + XSQ_R + YSQ_R
        ENDDO
        NORM_R = DSQRT(NORM_R)

!  Forward Sum-vector = normalized Right eigenvector

        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            FWD_R = RITE_EVEC(IROW,AA) /NORM_R
            FWD_I = RITE_EVEC(IROW,AA1)/NORM_R
            FWD_SUMVEC(I,O1,K1) = FWD_R
            FWD_SUMVEC(I,O1,K2) = FWD_I
          ENDDO
        ENDDO

!  Debug (important)

!        if ( do_debug_write ) then
!         if (m.eq.0.and.n.eq.6)then
!          write(90,'(2i3,1p2e20.10)')
!     &            m,k,KEIGEN(K1,N),KEIGEN(K2,N)
!          DO I = 1, NSTREAMS
!           DO O1 = 1, NSTOKES
!             write(90,'(4i3,1p4e20.10)')
!     &      m,k,i,o1,FWD_SUMVEC(I,O1,K1),FWD_SUMVEC(I,O1,K2)
!           ENDDO
!          ENDDO
!         ENDIF
!        endif

!  Find Forward difference vectors (Siewert's notation)

        DO I = 1, NSTREAMS
         DO O1 = 1, NSTOKES
          FWD_H1_CR = ZERO
          FWD_H1_CI = ZERO
          DO J = 1, NSTREAMS
           FWD_H2_CR = ZERO
           FWD_H2_CI = ZERO
           DO O2 = 1, NSTOKES
            FWD_H2_CR = FWD_H2_CR + &
                 SAB(I,J,O1,O2,N) * FWD_SUMVEC(J,O2,K1)
            FWD_H2_CI = FWD_H2_CI + &
                 SAB(I,J,O1,O2,N) * FWD_SUMVEC(J,O2,K2)
           ENDDO
           FWD_H1_CR = FWD_H1_CR + FWD_H2_CR
           FWD_H1_CI = FWD_H1_CI + FWD_H2_CI
          ENDDO
          FWD_DIFVEC(I,O1,K1) = FWD_H1_CR*SC_CR - FWD_H1_CI*SC_CI
          FWD_DIFVEC(I,O1,K2) = FWD_H1_CR*SC_CI + FWD_H1_CI*SC_CR
         ENDDO
        ENDDO

!  assign Forward solutions PHI (Siewert's notation)
!    --->  first N are "DOWN", last N are "UP" (streams)
!    --->  Use symmetry properties to set -ve eigensolutions

        DO I = 1, NSTREAMS
          XINV = HALF
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES
            FPD1 = XINV * ( FWD_SUMVEC(I,O1,K1) + FWD_DIFVEC(I,O1,K1) )
            FPD2 = XINV * ( FWD_SUMVEC(I,O1,K2) + FWD_DIFVEC(I,O1,K2) )
            FND1 = XINV * ( FWD_SUMVEC(I,O1,K1) - FWD_DIFVEC(I,O1,K1) )
            FND2 = XINV * ( FWD_SUMVEC(I,O1,K2) - FWD_DIFVEC(I,O1,K2) )
            SOLA_XPOS(I,O1,K1,N)   = FPD1
            SOLB_XNEG(I,O1,K1,N)   = FND1
            SOLA_XPOS(I,O1,K2,N)   = FPD2
            SOLB_XNEG(I,O1,K2,N)   = FND2
            SOLA_XPOS(I1,O1,K1,N)  = FND1
            SOLB_XNEG(I1,O1,K1,N)  = FPD1
            SOLA_XPOS(I1,O1,K2,N)  = FND2
            SOLB_XNEG(I1,O1,K2,N)  = FPD2
            IF ( O1 .GT. 2 ) THEN
              SOLA_XPOS(I1,O1,K1,N)  = - FND1
              SOLB_XNEG(I1,O1,K1,N)  = - FPD1
              SOLA_XPOS(I1,O1,K2,N)  = - FND2
              SOLB_XNEG(I1,O1,K2,N)  = - FPD2
            ENDIF
          ENDDO
        ENDDO

!  Older code: Used additional wasteful holding arrays
!   solution assignation (only for BVP)

!        DO I = 1, NSTREAMS
!          XINV = HALF
!          I1 = I + NSTREAMS
!          DO O1 = 1, NSTOKES
!            FWD_XPOS(I,O1,K1,N)  = XINV *
!     &         ( FWD_SUMVEC(I,O1,K1) + FWD_DIFVEC(I,O1,K1) )
!            FWD_XPOS(I,O1,K2,N)  = XINV *
!     &         ( FWD_SUMVEC(I,O1,K2) + FWD_DIFVEC(I,O1,K2) )
!            FWD_XPOS(I1,O1,K1,N) = XINV *
!     &         ( FWD_SUMVEC(I,O1,K1) - FWD_DIFVEC(I,O1,K1) )
!            FWD_XPOS(I1,O1,K2,N) = XINV *
!     &         ( FWD_SUMVEC(I,O1,K2) - FWD_DIFVEC(I,O1,K2) )
!            FWD_XNEG(I1,O1,K1,N) = FWD_XPOS(I,O1,K1,N)
!            FWD_XNEG(I,O1,K1,N)  = FWD_XPOS(I1,O1,K1,N)
!            FWD_XNEG(I1,O1,K2,N) = FWD_XPOS(I,O1,K2,N)
!            FWD_XNEG(I,O1,K2,N)  = FWD_XPOS(I1,O1,K2,N)
!          ENDDO
!        ENDDO
!        DO I = 1, NSTREAMS
!          I1 = I + NSTREAMS
!          DO O1 = 1, NSTOKES
!            SOLA_XPOS(I,O1,K1,N)  = FWD_XPOS(I,O1,K1,N)
!            SOLB_XNEG(I,O1,K1,N)  = FWD_XNEG(I,O1,K1,N)
!            SOLA_XPOS(I,O1,K2,N)  = FWD_XPOS(I,O1,K2,N)
!            SOLB_XNEG(I,O1,K2,N)  = FWD_XNEG(I,O1,K2,N)
!          ENDDO
!          DO O1 = 1, NMINST
!            SOLA_XPOS(I1,O1,K1,N)  = FWD_XNEG(I,O1,K1,N)
!            SOLB_XNEG(I1,O1,K1,N)  = FWD_XPOS(I,O1,K1,N)
!            SOLA_XPOS(I1,O1,K2,N)  = FWD_XNEG(I,O1,K2,N)
!            SOLB_XNEG(I1,O1,K2,N)  = FWD_XPOS(I,O1,K2,N)
!          ENDDO
!          DO O1 = 3, NSTOKES
!            SOLA_XPOS(I1,O1,K1,N)  = -FWD_XNEG(I,O1,K1,N)
!            SOLB_XNEG(I1,O1,K1,N)  = -FWD_XPOS(I,O1,K1,N)
!            SOLA_XPOS(I1,O1,K2,N)  = -FWD_XNEG(I,O1,K2,N)
!            SOLB_XNEG(I1,O1,K2,N)  = -FWD_XPOS(I,O1,K2,N)
!          ENDDO
!        ENDDO

!  Try to see if the solutions are the same
!    Multilayer case - checks out 17 December.
!      IF ( M.le.1 ) THEN
!      IF ( N.EQ.2 ) THEN
!       DO I = 1, 2*NSTREAMS
!        WRITE(M+1,'(2i3,1p4e16.6)')K,I,(SOLA_XPOS(I,O1,K1,N),O1=1,4)
!        WRITE(M+1,'(2i3,1p4e16.6)')k,I,(SOLA_XPOS(I,O1,K2,N),O1=1,4)
!        WRITE(M+1,'(2i3,1p4e16.6)')k,I,(SOLB_XNEG(I,O1,K1,N),O1=1,4)
!        WRITE(M+1,'(2i3,1p4e16.6)')k,I,(SOLB_XNEG(I,O1,K2,N),O1=1,4)
!       ENDDO
!      ENDIF
!      ENDIF

!  end eigenvalue loop

      ENDDO

!  control point for avoiding eigensolver
!  ======================================

 3456 CONTINUE

!  Eigenstream transmittance factors for whole layer
!  -------------------------------------------------

!  When the solution saving option is set, then if there is no
!  scattering in this layer, then linearized transmittances are
!  linearized discrete ordinate transmittances.

      IF ( DO_SOLUTION_SAVING .AND. &
                   .NOT.DO_LAYER_SCATTERING(M,N) ) THEN

        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            T_DELT_EIGEN(IROW,N) = T_DELT_DISORDS(I,N)
          ENDDO
        ENDDO

!  Otherwise compute them as normal

      ELSE

        DO K = 1, K_REAL(N)
          HELP = KEIGEN(K,N) * DELTAU_VERT(N)
          IF ( HELP .GT. MAX_TAU_QPATH ) THEN
            T_DELT_EIGEN(K,N) = ZERO
          ELSE
            T_DELT_EIGEN(K,N) = DEXP(-HELP)
          ENDIF
        ENDDO

        DO K = 1, K_COMPLEX(N)
          K0 = 2 * K - 2
          K1 = KO1 + K0
          K2 = K1  + 1
          HELP = DELTAU_VERT(N) * KEIGEN(K1,N)
          IF ( HELP .GT. MAX_TAU_QPATH ) THEN
            T_DELT_EIGEN(K1,N) = ZERO
            T_DELT_EIGEN(K2,N) = ZERO
          ELSE
            HELP_CR = DEXP ( - HELP)
            HELP_CI = DELTAU_VERT(N) * KEIGEN(K2,N)
            T_DELT_EIGEN(K1,N) =   HELP_CR * DCOS ( HELP_CI )
            T_DELT_EIGEN(K2,N) = - HELP_CR * DSIN ( HELP_CI )
          ENDIF
        ENDDO

      ENDIF

!      if ( m.eq.0.and.n.eq.20 ) then
!       DO K = 1, min(6,K_REAL(N))
!         WRITE(*,'(2i4,1p6e24.12)')m,k,T_DELT_EIGEN(K,N)
!       ENDDO
!      ENDIF

!  Eigenstream transmittance factors for partial layers
!  ----------------------------------------------------

! Loop over user optical depths

      DO UTA = 1, N_USER_LEVELS
       IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
        UT = PARTLAYERS_OUTINDEX(UTA)
        NA = PARTLAYERS_LAYERIDX(UT)

!  If the layer of occurence = given layer then do the calculation

         IF ( NA .EQ. N ) THEN

!  When the solution saving option is set, then if there is no
!  scattering in this layer, then transmittances are just the
!  discrete ordinate transmittances.

          IF ( DO_SOLUTION_SAVING .AND. &
               .NOT.DO_LAYER_SCATTERING(M,N) ) THEN

           DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
             IROW = IR + O1
             T_UTDN_EIGEN(IROW,UT) = T_DISORDS_UTDN(I,UT)
             T_UTUP_EIGEN(IROW,UT) = T_DISORDS_UTUP(I,UT)
            ENDDO
           ENDDO

!  Otherwise compute them as Eigenstream transmittance factors

          ELSE

!  partial layer optical depths

           TAU_DN = PARTAU_VERT(UT)
           TAU_UP = DELTAU_VERT(N) - TAU_DN

!  Real eigenvalues

           DO K = 1, K_REAL(N)
            HELP = KEIGEN(K,N) * TAU_DN
            IF ( HELP .GT. MAX_TAU_QPATH ) THEN
              T_UTDN_EIGEN(K,UT) = ZERO
            ELSE
              T_UTDN_EIGEN(K,UT) = DEXP(-HELP)
            ENDIF
            HELP = KEIGEN(K,N) * TAU_UP
            IF ( HELP .GT. MAX_TAU_QPATH ) THEN
              T_UTUP_EIGEN(K,UT) = ZERO
            ELSE
              T_UTUP_EIGEN(K,UT) = DEXP(-HELP)
            ENDIF
           ENDDO

!  Complex eigenvalues

           DO K = 1, K_COMPLEX(N)
            K0 = 2 * ( K - 1 )
            K1 = KO1 + K0
            K2 = K1  + 1
            HELP = TAU_DN * KEIGEN(K1,N)
            IF ( HELP .GT. MAX_TAU_QPATH ) THEN
              T_UTDN_EIGEN(K1,UT) = ZERO
              T_UTDN_EIGEN(K2,UT) = ZERO
            ELSE
              HELP_CR = DEXP ( - HELP )
              HELP_CI = TAU_DN * KEIGEN(K2,N)
              T_UTDN_EIGEN(K1,UT) =   HELP_CR * DCOS ( HELP_CI )
              T_UTDN_EIGEN(K2,UT) = - HELP_CR * DSIN ( HELP_CI )
            ENDIF
            HELP = TAU_UP * KEIGEN(K1,N)
            IF ( HELP .GT. MAX_TAU_QPATH ) THEN
              T_UTUP_EIGEN(K1,UT) = ZERO
              T_UTUP_EIGEN(K2,UT) = ZERO
            ELSE
              HELP_CR = DEXP ( - HELP )
              HELP_CI = TAU_UP * KEIGEN(K2,N)
              T_UTUP_EIGEN(K1,UT) =   HELP_CR * DCOS ( HELP_CI )
              T_UTUP_EIGEN(K2,UT) = - HELP_CR * DSIN ( HELP_CI )
            ENDIF
           ENDDO

!  End clause to use Discrete ordinates or not

          ENDIF

!  end loop over off-grid optical depths

         ENDIF
        ENDIF
      ENDDO

!  Debug

!        if ( m.lt.3.and.n.gt.0) then
!        write(97,'(2i5)')M,N
!        DO K = 1, NSTREAMS
!          WRITE(97,'(I3,1p6e17.9)')K,KEIGEN(K,N),
!     *          (SOLA_XPOS(I,1,K,N),I=1,NSTREAMS,2)
!        ENDDO
!        endif
!c        if (n.eq.26)pause
!      if (n.eq.6) pause'end hom'

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_QHOM_SOLUTION

!

      SUBROUTINE VLIDORT_UHOM_SOLUTION ( &
        GIVEN_LAYER, FOURIER, &
        DO_UPWELLING, DO_DNWELLING, &
        DO_DEBUG_WRITE, DO_FDTEST, &
        NSTOKES, NSTREAMS, &
        QUAD_HALFWTS, NMOMENTS, &
        N_USER_STREAMS, DO_LAYER_SCATTERING, &
        LOCAL_UM_START, USER_STREAMS, &
        USER_SECANTS, &
        OMEGA_GREEK, PI_XQP, &
        PI_XUP, PI_XUM, PI_XQM_PRE, &
        PI_XUM_POST, PI_XUP_PRE, &
        K_REAL, K_COMPLEX, KEIGEN, &
        KEIGEN_CSQ, SOLA_XPOS, &
        SOLB_XNEG, &
        UHOM_DNDN, UHOM_DNUP, &
        UHOM_UPDN, UHOM_UPUP, &
        HSINGO, HELPSTOKES, &
        ZETA_M, ZETA_P )

!  Small-number Taylor series expansions, added 30 October 2007

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::           GIVEN_LAYER
      INTEGER, INTENT (IN) ::           FOURIER

      LOGICAL, INTENT (IN) ::           DO_UPWELLING
      LOGICAL, INTENT (IN) ::           DO_DNWELLING
      LOGICAL, INTENT (IN) ::           DO_DEBUG_WRITE
      LOGICAL, INTENT (IN) ::           DO_FDTEST
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_HALFWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::           NMOMENTS
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      LOGICAL, INTENT (IN) ::           DO_LAYER_SCATTERING &
          ( 0:MAXMOMENTS, MAXLAYERS )
      INTEGER, INTENT (IN) ::           LOCAL_UM_START
      DOUBLE PRECISION, INTENT (IN) ::  USER_STREAMS  ( MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  USER_SECANTS  ( MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  OMEGA_GREEK &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XQP &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XUP &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XUM &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XQM_PRE &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XUM_POST &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XUP_PRE &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      INTEGER, INTENT (IN) ::           K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  KEIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  KEIGEN_CSQ ( MAXEVALUES )
      DOUBLE PRECISION, INTENT (IN) ::  SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (INOUT) :: UHOM_DNDN &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: UHOM_DNUP &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: UHOM_UPDN &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: UHOM_UPUP &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      LOGICAL, INTENT (INOUT) ::          HSINGO &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: HELPSTOKES &
          ( 0:MAXMOMENTS, MAXEVALUES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) :: ZETA_M &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: ZETA_P &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )

!  Local variables
!  ---------------

      INTEGER ::          UM, J, L, N, M, K, O1, O2, KO1, K0, K1, K2

      DOUBLE PRECISION :: GAUX_R  ( 0:MAXMOMENTS, MAXSTOKES )
      DOUBLE PRECISION :: GAUX_CR ( 0:MAXMOMENTS, MAXSTOKES )
      DOUBLE PRECISION :: GAUX_CI ( 0:MAXMOMENTS, MAXSTOKES )

!  for real solutions

      DOUBLE PRECISION :: SM, SN, SP, SNS, SPS, A5
      DOUBLE PRECISION :: RHO_P, RHO_M, SECMUI, SUMG

!  for complex solutions

      DOUBLE PRECISION :: RHO_P_CR, RHO_M_CR, MODULUS_P, MODULUS_M
      DOUBLE PRECISION :: A_USER, A_HELP, B_USER, B_HELP, HPC, HMC
      DOUBLE PRECISION :: SMCR, SNCR, SPCR, SNSCR, SPSCR, SGCR
      DOUBLE PRECISION :: SMCI, SNCI, SPCI, SNSCI, SPSCI, SGCI

!  Layer and Fourier

      N = GIVEN_LAYER
      M = FOURIER

!  Initialization

      K   = 0

!  Zeta constants
!  ==============

!  Real Zeta constants (always required)

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        SECMUI = USER_SECANTS(UM)
        DO K = 1, K_REAL(N)
          RHO_P = SECMUI + KEIGEN(K,N)
          RHO_M = SECMUI - KEIGEN(K,N)
          ZETA_P(K,UM,N) = ONE / RHO_P
          IF ( DABS(RHO_M) .LT. SMALLNUM ) THEN
            HSINGO(K,UM,N) = .TRUE.
            ZETA_M(K,UM,N) = ONE
          ELSE
            HSINGO(K,UM,N) = .FALSE.
            ZETA_M(K,UM,N) = ONE / RHO_M
          ENDIF
        ENDDO
      ENDDO

!  Complex Zeta constants

      IF ( K_COMPLEX(N) .GT. 0 ) THEN
        KO1 = K_REAL(N) + 1
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          A_USER = USER_SECANTS(UM) * USER_SECANTS(UM)
          B_USER = TWO * USER_SECANTS(UM)
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            RHO_P_CR = USER_SECANTS(UM) + KEIGEN(K1,N)
            RHO_M_CR = USER_SECANTS(UM) - KEIGEN(K1,N)
            A_HELP = A_USER + KEIGEN_CSQ(K)
            B_HELP = B_USER * KEIGEN(K1,N)
            MODULUS_P = A_HELP + B_HELP
            MODULUS_M = A_HELP - B_HELP
            ZETA_P(K1,UM,N) = RHO_P_CR / MODULUS_P
            ZETA_P(K2,UM,N) = -KEIGEN(K2,N) / MODULUS_P
            IF ( DABS(MODULUS_M) .LT. SMALLNUM ) THEN
              HSINGO(K1,UM,N) = .TRUE.
              HSINGO(K2,UM,N) = .TRUE.
              ZETA_M(K1,UM,N) = ONE
              ZETA_M(K2,UM,N) = ZERO
            ELSE
              HSINGO(K1,UM,N) = .FALSE.
              HSINGO(K2,UM,N) = .FALSE.
              ZETA_M(K1,UM,N) = RHO_M_CR / MODULUS_M
              ZETA_M(K2,UM,N) = KEIGEN(K2,N) / MODULUS_M
            ENDIF

          ENDDO
        ENDDO
      ENDIF

!  If there is no scattering avoid these solutions
!  Warning - Solutions are zeroed in the scalar code before returning
!  @@@ Rob Fix. Zero them here too.

!mick fix 2/1/2011 - modified IF to initialize solutions
!      IF ( .NOT. DO_LAYER_SCATTERING(M,N) ) RETURN

      IF ( .NOT. DO_LAYER_SCATTERING(M,N) ) THEN
       IF ( DO_DNWELLING ) THEN
          UHOM_DNDN(:,:,:,N) = ZERO
          UHOM_DNUP(:,:,:,N) = ZERO
       ENDIF
       IF ( DO_UPWELLING ) THEN
          UHOM_UPDN(:,:,:,N) = ZERO
          UHOM_UPUP(:,:,:,N) = ZERO
       ENDIF
       HELPSTOKES(M,:,:)  = ZERO
       RETURN
      ENDIF

!  User defined solutions (Real)
!  -----------------------------

      DO K = 1, K_REAL(N)

!  For each moment, do inner sum over computational angles
!  for the positive and negative eigenvectors.
!    Gvec = STOKES. Eq 34. GAUx = Greekmat.Gvec

       DO L = M, NMOMENTS
        DO O1 = 1, NSTOKES
         SPS = ZERO
         DO J = 1, NSTREAMS
          SP = ZERO
          SM = ZERO
          A5 = QUAD_HALFWTS(J)
          DO O2 = 1, NSTOKES
!  older   SP = SP + PI_XQP    (L,J,O1,O2)*A5(J)*FWD_XPOS(J,O2,K,N)
!  older   SM = SM + PI_XQM_PRE(L,J,O1,O2)*A5(J)*FWD_XNEG(J,O2,K,N)
           SP = SP + PI_XQP    (L,J,O1,O2)*A5*SOLA_XPOS(J,O2,K,N)
           SM = SM + PI_XQM_PRE(L,J,O1,O2)*A5*SOLB_XNEG(J,O2,K,N)
          ENDDO
          SPS = SPS + SP + SM
         ENDDO
         HELPSTOKES(L,K,O1) = SPS
        ENDDO
        DO O1 = 1, NSTOKES
         SUMG = ZERO
         DO O2 = 1, NSTOKES
          SUMG = SUMG + OMEGA_GREEK(L,N,O1,O2) * HELPSTOKES(L,K,O2)
         ENDDO
         GAUX_R(L,O1) = SUMG
        ENDDO
       ENDDO

!  Now sum over all harmonic contributions (downwelling)

       IF ( DO_DNWELLING ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
         DO O1 = 1, NSTOKES
          SP = ZERO
          SN = ZERO
          DO L = M, NMOMENTS
           SPS = ZERO
           SNS = ZERO
           DO O2 = 1, NSTOKES
            SPS = SPS + PI_XUP(L,UM,O1,O2)      * GAUX_R(L,O2)
            SNS = SNS + PI_XUM_POST(L,UM,O1,O2) * GAUX_R(L,O2)
           ENDDO
           SP = SP + SPS
           SN = SN + SNS
          ENDDO
          UHOM_DNDN(UM,O1,K,N) = SP
          UHOM_DNUP(UM,O1,K,N) = SN
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
            SPS = SPS + PI_XUP_PRE(L,UM,O1,O2)  * GAUX_R(L,O2)
            SNS = SNS + PI_XUM(L,UM,O1,O2)      * GAUX_R(L,O2)
           ENDDO
           SP = SP + SPS
           SN = SN + SNS
          ENDDO
          UHOM_UPDN(UM,O1,K,N) = SN
          UHOM_UPUP(UM,O1,K,N) = SP
         ENDDO
        ENDDO
       ENDIF

!  Debug (important)

!        IF ( DO_DEBUG_WRITE .AND. DO_FDTEST) THEN
!         DO UM = LOCAL_UM_START, N_USER_STREAMS
!          DO O1 = 1, NSTOKES
!            WRITE(94,'(4I3,1P4E20.10)')M,K,UM,O1, &
!               UHOM_UPUP(UM,O1,K,N), &
!               UHOM_UPDN(UM,O1,K,N), &
!               UHOM_DNUP(UM,O1,K,N), &
!               UHOM_DNDN(UM,O1,K,N)
!           ENDDO
!         ENDDO
!        ENDIF

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
!           SPCR = SPCR + HPC * FWD_XPOS(J,O2,K1,N)         ! Older
!           SPCI = SPCI + HPC * FWD_XPOS(J,O2,K2,N)         ! Older
           SPCR = SPCR + HPC * SOLA_XPOS(J,O2,K1,N)
           SPCI = SPCI + HPC * SOLA_XPOS(J,O2,K2,N)
           HMC  =  PI_XQM_PRE(L,J,O1,O2)*QUAD_HALFWTS(J)
!           SMCR = SMCR + HMC * FWD_XNEG(J,O2,K1,N)         ! Older
!           SMCI = SMCI + HMC * FWD_XNEG(J,O2,K2,N)         ! Older
           SMCR = SMCR + HMC * SOLB_XNEG(J,O2,K1,N)
           SMCI = SMCI + HMC * SOLB_XNEG(J,O2,K2,N)
          ENDDO
          SPSCR = SPSCR + SPCR + SMCR
          SPSCI = SPSCI + SPCI + SMCI
         ENDDO
         HELPSTOKES(L,K1,O1) = SPSCR
         HELPSTOKES(L,K2,O1) = SPSCI
        ENDDO
        DO O1 = 1, NSTOKES
         SGCR = ZERO
         SGCI = ZERO
         DO O2 = 1, NSTOKES
          SGCR = SGCR + OMEGA_GREEK(L,N,O1,O2) * HELPSTOKES(L,K1,O2)
          SGCI = SGCI + OMEGA_GREEK(L,N,O1,O2) * HELPSTOKES(L,K2,O2)
         ENDDO
         GAUX_CR(L,O1) = SGCR
         GAUX_CI(L,O1) = SGCI
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
            SPSCR = SPSCR + HPC * GAUX_CR(L,O2)
            SPSCI = SPSCI + HPC * GAUX_CI(L,O2)
            HMC = PI_XUM_POST(L,UM,O1,O2)
            SNSCR = SNSCR + HMC * GAUX_CR(L,O2)
            SNSCI = SNSCI + HMC * GAUX_CI(L,O2)
           ENDDO
           SPCR = SPCR + SPSCR
           SPCI = SPCI + SPSCI
           SNCR = SNCR + SNSCR
           SNCI = SNCI + SNSCI
          ENDDO
          UHOM_DNDN(UM,O1,K1,N) = SPCR
          UHOM_DNUP(UM,O1,K1,N) = SNCR
          UHOM_DNDN(UM,O1,K2,N) = SPCI
          UHOM_DNUP(UM,O1,K2,N) = SNCI
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
            SPSCR = SPSCR + HPC * GAUX_CR(L,O2)
            SPSCI = SPSCI + HPC * GAUX_CI(L,O2)
            HMC = PI_XUM(L,UM,O1,O2)
            SNSCR = SNSCR + HMC * GAUX_CR(L,O2)
            SNSCI = SNSCI + HMC * GAUX_CI(L,O2)
           ENDDO
           SPCR = SPCR + SPSCR
           SPCI = SPCI + SPSCI
           SNCR = SNCR + SNSCR
           SNCI = SNCI + SNSCI
          ENDDO
          UHOM_UPDN(UM,O1,K1,N) = SNCR
          UHOM_UPUP(UM,O1,K1,N) = SPCR
          UHOM_UPDN(UM,O1,K2,N) = SNCI
          UHOM_UPUP(UM,O1,K2,N) = SPCI
         ENDDO
        ENDDO
       ENDIF

!  Debug (important)
!        if ( do_debug_write.and.do_fdtest) then
!         DO UM = LOCAL_UM_START, N_USER_STREAMS
!          DO O1 = 1, NSTOKES
!            write(96,'(4i3,1p4e20.10)')m,k,um,o1,
!     &         UHOM_UPUP(UM,O1,K2,N),
!     &         UHOM_UPDN(UM,O1,K2,N),
!     &         UHOM_DNUP(UM,O1,K2,N),
!     &         UHOM_DNDN(UM,O1,K2,N)
!           enddo
!         enddo
!        endif

!  Try to see if the solutions are the same
!    Multilayer case - checks out 17 December.
!        IF ( M.eq.1 ) THEN
!        IF ( N.EQ.2 ) THEN
!         DO UM = 1, N_USER_STREAMS
!          WRITE(N,'(2i3,1p4e16.6)')K,UM,(UHOM_UPUP(UM,O1,K1,N),O1=1,4)
!          WRITE(N,'(2i3,1p4e16.6)')K,UM,(UHOM_UPUP(UM,O1,K2,N),O1=1,4)
!          WRITE(N,'(2i3,1p4e16.6)')K,UM,(UHOM_UPDN(UM,O1,K1,N),O1=1,4)
!          WRITE(N,'(2i3,1p4e16.6)')K,UM,(UHOM_UPDN(UM,O1,K2,N),O1=1,4)
!          WRITE(N,'(2i3,1p4e16.6)')K,UM,(UHOM_DNUP(UM,O1,K1,N),O1=1,4)
!          WRITE(N,'(2i3,1p4e16.6)')K,UM,(UHOM_DNUP(UM,O1,K2,N),O1=1,4)
!          WRITE(N,'(2i3,1p4e16.6)')K,UM,(UHOM_DNDN(UM,O1,K1,N),O1=1,4)
!          WRITE(N,'(2i3,1p4e16.6)')K,UM,(UHOM_DNDN(UM,O1,K2,N),O1=1,4)
!         ENDDO
!         ENDIF
!        ENDIF

!  end loop over complex eigensolutions

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_UHOM_SOLUTION

!

      SUBROUTINE VLIDORT_QBEAM_SOLUTION ( &
        GIVEN_LAYER, FOURIER, IBEAM, &
        DO_DEBUG_WRITE, DO_FDTEST, &
        NSTOKES, NSTREAMS, FLUX_FACTOR, &
        LAYER_PIS_CUTOFF, QUAD_STREAMS, &
        NMOMENTS, NSTREAMS_2, &
        NSTKS_NSTRMS, DO_LAYER_SCATTERING, &
        DMAT, DFLUX, &
        OMEGA_GREEK, T_DELT_MUBAR, &
        INITIAL_TRANS, AVERAGE_SECANT, &
        PI_XQP, PI_XQM, &
        PI_X0P, PI_XQM_POST, &
        SAB, DAB, &
        EIGENMAT_SAVE, &
        QSUMVEC_SAVE, QDIFVEC_SAVE, &
        QVEC_SAVE, QDIF_SAVE, &
        QMAT_SAVE, QPIVOT, &
        BVEC, WUPPER, WLOWER, &
        STATUS, MESSAGE, TRACE )

!  This is the classical Chandrasekhar beam solution.
!  ( plane parallel or average secant only)
!  Linear Matrix algebra.

      USE VLIDORT_PARS
      USE LAPACK_TOOLS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::           GIVEN_LAYER
      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           IBEAM
      LOGICAL, INTENT (IN) ::           DO_DEBUG_WRITE
      LOGICAL, INTENT (IN) ::           DO_FDTEST
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      DOUBLE PRECISION, INTENT (IN) ::  FLUX_FACTOR
      INTEGER, INTENT (IN) ::           LAYER_PIS_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_STREAMS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::           NMOMENTS
      INTEGER, INTENT (IN) ::           NSTREAMS_2
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS
      LOGICAL, INTENT (IN) ::           DO_LAYER_SCATTERING &
          ( 0:MAXMOMENTS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  DMAT ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  DFLUX ( MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  OMEGA_GREEK &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XQP &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XQM &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_X0P &
          ( 0:MAXMOMENTS, MAXBEAMS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XQM_POST &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  SAB &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  DAB &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  EIGENMAT_SAVE &
          ( MAXEVALUES, MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (INOUT) ::  QSUMVEC_SAVE &
          ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) ::  QDIFVEC_SAVE &
          ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) ::  QVEC_SAVE ( MAXSTRMSTKS )
      DOUBLE PRECISION, INTENT (INOUT) ::  QDIF_SAVE ( MAXSTRMSTKS )
      DOUBLE PRECISION, INTENT (INOUT) ::  QMAT_SAVE &
          ( MAXSTRMSTKS, MAXSTRMSTKS )
      INTEGER, INTENT (INOUT) ::           QPIVOT ( MAXSTRMSTKS )
      DOUBLE PRECISION, INTENT (INOUT) ::  BVEC &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) ::  WUPPER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) ::  WLOWER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )

      INTEGER, INTENT (OUT) ::           STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  Local variables
!  ---------------

!  help variables

      INTEGER ::          I, J, I1, L, N, M
      INTEGER ::          O1, O2, O3, IR, JR, IROW, JROW
      DOUBLE PRECISION :: SUM, GSUM, INV_X0SQ, SECBAR, F1
      DOUBLE PRECISION :: WSUM, WDIF, TRANS1, TRANS2, WVEC_D, WVEC_U
      DOUBLE PRECISION :: TPA, TMA, S_TPA, S_TMA, S_HELP, HELP

      DOUBLE PRECISION :: HELP_QFUNC ( MAXLAYERS, 0:MAXMOMENTS, MAXSTOKES )

!  Local error handling

      INTEGER ::           INFO
      CHARACTER (LEN=3) :: CI, C3
      CHARACTER (LEN=2) :: CB

!  Initial section
!  ---------------

!  Initialize status

      STATUS  = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  initialise indices

      N = GIVEN_LAYER
      M = FOURIER
      F1 = FLUX_FACTOR / PI4

!  No particular solution beyond the cutoff layer
!  OR No scattering in this layer then:
!    ---> Zero the boundary layer values and exit

      IF ( N .GT. LAYER_PIS_CUTOFF(IBEAM) .OR. &
             .NOT. DO_LAYER_SCATTERING(M,N) ) THEN
        DO I = 1, NSTREAMS_2
          DO O1 = 1, NSTOKES
            WUPPER(I,O1,N) = ZERO
            WLOWER(I,O1,N) = ZERO
          ENDDO
        ENDDO

!mick fix
        QSUMVEC_SAVE = ZERO
        QDIFVEC_SAVE = ZERO
        QVEC_SAVE = ZERO
        QDIF_SAVE = ZERO
        QMAT_SAVE = ZERO
        QPIVOT = ZERO
        BVEC(1:NSTREAMS_2,1:NSTOKES,N) = ZERO

        RETURN
      ENDIF

!  set local values

      SECBAR   = AVERAGE_SECANT(N,IBEAM)
      INV_X0SQ = SECBAR * SECBAR

!  Initialise matrix and column vector for solution
!  Changed at FMI 28 september 2004

      DO J = 1, NSTKS_NSTRMS
        DO I = 1, NSTKS_NSTRMS
          QMAT_SAVE(I,J) = ZERO
        ENDDO
        QVEC_SAVE(J) = ZERO
      ENDDO

!  Set up sum and difference vectors for Beam source terms
!  ( sum vector may be required again in linearization )
!  Auxiliary matrix for Q functions

      DO L = M, NMOMENTS
        DO O1 = 1, NSTOKES
          SUM = ZERO
          DO O2 = 1, NSTOKES
            GSUM = ZERO
            DO O3 = 1, NSTOKES
              GSUM = GSUM + PI_X0P(L,IBEAM,N,O2,O3) * DFLUX(O3)
            ENDDO
            SUM = SUM + OMEGA_GREEK(L,N,O1,O2) * GSUM
          ENDDO
          HELP_QFUNC(N,L,O1) = SUM
        ENDDO
      ENDDO

      DO I = 1, NSTREAMS
        DO O1 = 1, NSTOKES
          S_TPA = ZERO
          S_TMA = ZERO
          DO L = M, NMOMENTS
            TPA = ZERO
            TMA = ZERO
            DO O2 = 1, NSTOKES
              TPA = TPA + PI_XQP(L,I,O1,O2)      * HELP_QFUNC(N,L,O2)
              TMA = TMA + PI_XQM_POST(L,I,O1,O2) * HELP_QFUNC(N,L,O2)
            ENDDO
            S_TPA = S_TPA + TPA
            S_TMA = S_TMA + TMA
          ENDDO
          QSUMVEC_SAVE(I,O1) = F1 * ( S_TPA + S_TMA ) / QUAD_STREAMS(I)
          QDIFVEC_SAVE(I,O1) = F1 * ( S_TPA - S_TMA ) / QUAD_STREAMS(I)
        ENDDO
      ENDDO

!  solution matrix for the reduced problem
!  ( matrix should be saved in the LU decomposition form)
!  Changed at FMI 28 september

      DO J = 1, NSTKS_NSTRMS
        DO I = 1, NSTKS_NSTRMS
          QMAT_SAVE(I,J) = EIGENMAT_SAVE(I,J,N)
        ENDDO
      ENDDO

      DO I = 1, NSTKS_NSTRMS
        QMAT_SAVE(I,I) = QMAT_SAVE(I,I) - INV_X0SQ
      ENDDO

!  RHS vector for the reduced problem
!  ( this vector will be the answer after the linear algebra solution,
!    and may be needed again if there is linearization )

      DO I = 1, NSTREAMS
        IR = ( I - 1 ) * NSTOKES
        DO O1 = 1, NSTOKES
          IROW = IR + O1
          HELP = ZERO
          DO J = 1, NSTREAMS
            S_HELP = ZERO
            DO O2 = 1, NSTOKES
              S_HELP = S_HELP + DAB(I,J,O1,O2,N)*QSUMVEC_SAVE(J,O2)
            ENDDO
            HELP = HELP + S_HELP
          ENDDO
          QVEC_SAVE(IROW) = HELP + QDIFVEC_SAVE(I,O1) * SECBAR
        ENDDO
      ENDDO

!  debug

!      write(*,'(i3)')N
!      DO I = 1, NSTREAMS
!       write(*,'(i3,1p5e11.3)')I,QVEC_SAVE(I),(QMAT_SAVE(I,J),J=1,4)
!       ENDDO
!      if (n.eq.2)pause'qsetup'

!  L-U decomposition of the solution matrix

      CALL DGETRF &
            ( NSTKS_NSTRMS, NSTKS_NSTRMS, &
              QMAT_SAVE, MAXSTRMSTKS, QPIVOT, INFO )

      IF ( INFO .GT. 0 ) THEN
        WRITE(CI, '(I3)' ) INFO
        WRITE(CB, '(I2)' ) IBEAM
        WRITE(C3, '(I3)' ) N
        MESSAGE = 'argument i illegal value, for i = '//CI
        TRACE   = 'DGETRF call in VLIDORT_QBEAM_SOLUTION, Beam/layer'// &
            ' numbers = '//CB//', '//C3
        STATUS  = VLIDORT_SERIOUS
        RETURN
      ENDIF

!  Solution of reduced problem by back-substitution

      CALL DGETRS &
          ( 'N', NSTKS_NSTRMS, 1, &
             QMAT_SAVE, MAXSTRMSTKS, QPIVOT, &
             QVEC_SAVE, MAXSTRMSTKS, INFO )

      IF ( INFO .LT. 0 ) THEN
        WRITE(CI, '(I3)' ) INFO
        WRITE(CB, '(I2)' ) IBEAM
        WRITE(C3, '(I3)' ) N
        MESSAGE = 'argument i illegal value, for i = '//CI
        TRACE   = 'DGETRS call in VLIDORT_QBEAM_SOLUTION, Beam/layer'// &
            ' numbers = '//CB//', '//C3
        STATUS  = VLIDORT_SERIOUS
        RETURN
      ENDIF

!  Assigning beam particular integral solution vector

      DO I = 1, NSTREAMS
        I1 = I + NSTREAMS
        IR = ( I - 1 ) * NSTOKES
        DO O1 = 1, NSTOKES
          IROW = IR + O1
          HELP = ZERO
          DO J = 1, NSTREAMS
            S_HELP = ZERO
            JR = ( J - 1 ) * NSTOKES
            DO O2 = 1, NSTOKES
              JROW = JR + O2
              S_HELP = S_HELP + SAB(I,J,O1,O2,N)*QVEC_SAVE(JROW)
            ENDDO
            HELP = HELP + S_HELP
          ENDDO
          QDIF_SAVE(IROW) = ( HELP - QSUMVEC_SAVE(I,O1) ) / SECBAR
          WSUM = QVEC_SAVE(IROW)
          WDIF = QDIF_SAVE(IROW)
          WVEC_D = HALF * ( WSUM + WDIF )
          WVEC_U = HALF * ( WSUM - WDIF )
          BVEC(I,O1,N)  = WVEC_D
          BVEC(I1,O1,N) = WVEC_U
          IF ( O1 .GT. 2 ) THEN
            BVEC(I1,O1,N) = - WVEC_U
          ENDIF
!          if(i.eq.4.and.n.eq.24)write(38,'(4i3,1p2e20.10)')
!     &           m,ibeam,i,O1,BVEC(I,O1,N)
        ENDDO
      ENDDO

!  Old code: solution assignation

!      DO I = 1, NSTREAMS
!        I1 = I + NSTREAMS
!        DO O1 = 1, NSTOKES
!          BVEC(I,O1,N) = WVEC(I,O1,N)
!          HELP = ZERO
!          DO O2 = 1, NSTOKES
!            HELP = HELP + DMAT(O1,O2) * WVEC(I1,O2,N)
!          ENDDO
!          BVEC(I1,O1,N) = HELP
!        ENDDO
!      ENDDO

!  Values at the layer boundaries
!  (transmittance factors have been determined in SETUPS module)

      TRANS1 = INITIAL_TRANS(N,IBEAM)
      TRANS2 = T_DELT_MUBAR(N,IBEAM) * TRANS1

      DO I = 1, NSTREAMS_2
        DO O1 = 1, NSTOKES
          WUPPER(I,O1,N) = BVEC(I,O1,N)*TRANS1
          WLOWER(I,O1,N) = BVEC(I,O1,N)*TRANS2
        ENDDO
      ENDDO

!  debug

!      IF ( DO_DEBUG_WRITE ) THEN
!       J = 97
!       IF ( DO_FDTEST ) J = 98
!       IF ( N.EQ.14.and.m.eq.1) THEN
!        DO I = 1, NSTREAMS*2
!         write(*,'(4i4,1p4e17.9)')IBEAM,M,N,I,BVEC(I,1,N),
!     &      INITIAL_TRANS(N,IBEAM), T_DELT_MUBAR(N,IBEAM)
!        ENDDO
!       ENDIF
!       IF ( DO_FDTEST. and. M.EQ.1) PAUSE
!      ENDIF

!  debug check solution Forward problem - Works
!    removed 2 July. See Older version, labeled _debug

!      IF ( N.le.6.and.m.eq.0) THEN
!       DO I = 1,40
!        write(91,'(3i3,1p4e17.7)')M,N,I,(BVEC(I,o1,N),o1=1,4)
!       ENDDO
!      ENDIF

!      if (n.eq.6) pause 'end beam'

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_QBEAM_SOLUTION

!

      SUBROUTINE VLIDORT_UBEAM_SOLUTION ( &
        GIVEN_LAYER, FOURIER, IBEAM, &
        DO_UPWELLING, DO_DNWELLING, &
        NSTOKES, NSTREAMS, FLUX_FACTOR, &
        LAYER_PIS_CUTOFF, QUAD_HALFWTS, &
        NMOMENTS, N_USER_STREAMS, &
        DO_LAYER_SCATTERING, LOCAL_UM_START, &
        STERM_LAYERMASK_UP, &
        STERM_LAYERMASK_DN, DFLUX, &
        OMEGA_GREEK, PI_XQP, PI_XUP, &
        PI_XUM, PI_X0P, PI_XQM_PRE, &
        BVEC, HELPSTOKES_BEAM, &
        UPAR_DN_1, UPAR_DN_2, &
        UPAR_UP_1, UPAR_UP_2 )

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::           GIVEN_LAYER
      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           IBEAM
      LOGICAL, INTENT (IN) ::           DO_UPWELLING
      LOGICAL, INTENT (IN) ::           DO_DNWELLING
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      DOUBLE PRECISION, INTENT (IN) ::  FLUX_FACTOR
      INTEGER, INTENT (IN) ::           LAYER_PIS_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_HALFWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::           NMOMENTS
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      LOGICAL, INTENT (IN) ::           DO_LAYER_SCATTERING &
          ( 0:MAXMOMENTS, MAXLAYERS )
      INTEGER, INTENT (IN) ::           LOCAL_UM_START
      LOGICAL, INTENT (IN) ::           STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::           STERM_LAYERMASK_DN ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  DFLUX ( MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  OMEGA_GREEK &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XQP &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XUP &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XUM &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_X0P &
          ( 0:MAXMOMENTS, MAXBEAMS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XQM_PRE &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  BVEC &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (INOUT) :: HELPSTOKES_BEAM &
          ( 0:MAXMOMENTS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) :: UPAR_DN_1 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: UPAR_DN_2 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: UPAR_UP_1 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: UPAR_UP_2 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )

!  Local variables
!  ---------------

      LOGICAL ::          DO_LAYER
      INTEGER ::          UM, L, N, M, O1, O2, O3, J, J1
      DOUBLE PRECISION :: SUM1, SUM2, T1, T2, A5, F1
      DOUBLE PRECISION :: SP, SM, SPS, SPT

      DOUBLE PRECISION :: HELP_Q1 ( 0:MAXMOMENTS, MAXSTOKES )
      DOUBLE PRECISION :: HELP_Q2 ( 0:MAXMOMENTS, MAXSTOKES )

      DOUBLE PRECISION, DIMENSION(MAXSTOKES)  :: &
        SGW =  (/ 1.0d0, 1.0d0, -1.0d0, -1.0d0 /)

!  Layer and Fourier

      N = GIVEN_LAYER
      M = FOURIER
      F1 = FLUX_FACTOR / PI4

!  No particular solution beyond the cutoff layer.
!  OR No scattering in this layer
!    --> Zero the boundary layer values and exit

      IF ( N .GT. LAYER_PIS_CUTOFF(IBEAM) .OR. &
             .NOT. DO_LAYER_SCATTERING(M,N) ) THEN
        IF ( DO_UPWELLING ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              UPAR_UP_1(UM,O1,N) = ZERO
              UPAR_UP_2(UM,O1,N) = ZERO
            ENDDO
          ENDDO
        ENDIF
        IF ( DO_DNWELLING ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              UPAR_DN_1(UM,O1,N) = ZERO
              UPAR_DN_2(UM,O1,N) = ZERO
            ENDDO
          ENDDO
        ENDIF
        RETURN
      ENDIF

!  existence flag

      DO_LAYER = ( STERM_LAYERMASK_UP(N) .OR. &
                   STERM_LAYERMASK_DN(N) )

!  first function

      IF ( DO_LAYER ) THEN
       DO L = M, NMOMENTS
        DO O1 = 1, NSTOKES
          SPS = ZERO
          DO O2 = 1, NSTOKES
            SP = ZERO
            DO O3 = 1, NSTOKES
              SP = SP + PI_X0P(L,IBEAM,N,O2,O3) * DFLUX(O3)
            ENDDO
            SPS = SPS + OMEGA_GREEK(L,N,O1,O2) * SP
          ENDDO
          HELP_Q1(L,O1) = F1 * SPS
        ENDDO
       ENDDO
      ENDIF

!  For each moment, do inner sum over computational angles
!    Equivalent to Gvec = GVECTOR. Eq 34. GAUx = Greekmat.Gvec

      IF ( DO_LAYER ) THEN
        DO L = M, NMOMENTS
         DO O1 = 1, NSTOKES
          SPS = ZERO
          DO J = 1, NSTREAMS
           J1 = J + NSTREAMS
           SP = ZERO
           SM = ZERO
           A5 = QUAD_HALFWTS(J)
           DO O2 = 1, NSTOKES
            SP = SP +         PI_XQP    (L,J,O1,O2)*A5*BVEC(J,O2,N)
            SM = SM + SGW(O2)*PI_XQM_PRE(L,J,O1,O2)*A5*BVEC(J1,O2,N)
           ENDDO
           SPS = SPS + SP + SM
          ENDDO
          HELPSTOKES_BEAM(L,O1) = SPS
         ENDDO
         DO O1 = 1, NSTOKES
          SPT = ZERO
          DO O2 = 1, NSTOKES
           SPT = SPT + OMEGA_GREEK(L,N,O1,O2) * HELPSTOKES_BEAM(L,O2)
          ENDDO
          HELP_Q2(L,O1) = SPT
         ENDDO
        ENDDO
      ENDIF

!  Now sum over all harmonic contributions (Upwelling)
!    Direct and integarated contributions

      IF ( DO_UPWELLING ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
         DO O1 = 1, NSTOKES
          T1 = ZERO
          T2 = ZERO
          DO L = M, NMOMENTS
           SUM1 = ZERO
           SUM2 = ZERO
           DO O2 = 1, NSTOKES
            SUM1 = SUM1 + HELP_Q1(L,O2)*PI_XUM(L,UM,O1,O2)
            SUM2 = SUM2 + HELP_Q2(L,O2)*PI_XUM(L,UM,O1,O2)
           ENDDO
           T1 = T1 + SUM1
           T2 = T2 + SUM2
          ENDDO
          UPAR_UP_1(UM,O1,N) = T1
          UPAR_UP_2(UM,O1,N) = T2
         ENDDO
        ENDDO
      ENDIF

!  Now sum over all harmonic contributions (Downwelling)
!    Direct and integarated contributions

      IF ( DO_DNWELLING ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
         DO O1 = 1, NSTOKES
          T1 = ZERO
          T2 = ZERO
          DO L = M, NMOMENTS
           SUM1 = ZERO
           SUM2 = ZERO
           DO O2 = 1, NSTOKES
            SUM1 = SUM1 + HELP_Q1(L,O2)*PI_XUP(L,UM,O1,O2)
            SUM2 = SUM2 + HELP_Q2(L,O2)*PI_XUP(L,UM,O1,O2)
           ENDDO
           T1 = T1 + SUM1
           T2 = T2 + SUM2
          ENDDO
          UPAR_DN_1(UM,O1,N) = T1
          UPAR_DN_2(UM,O1,N) = T2
         ENDDO
        ENDDO
      ENDIF

!  debug stokes = 1

!      DO UM = LOCAL_UM_START, N_USER_STREAMS
!        write(*,*)n,um,UPAR_UP_1(UM,1,N),UPAR_UP_2(UM,1,N)
!      ENDDO
!      DO UM = LOCAL_UM_START, N_USER_STREAMS
!        write(*,*)n,um,UPAR_DN_1(UM,1,N),UPAR_DN_2(UM,1,N)
!      ENDDO
!      PAUSE

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_UBEAM_SOLUTION

!

      SUBROUTINE VLIDORT_PIMATRIX_SETUP ( &
        FOURIER, &
        DO_REFRACTIVE_GEOMETRY, NSTOKES, &
        NSTREAMS, NLAYERS, &
        COS_SZANGLES, SUN_SZA_COSINES, &
        QUAD_STREAMS, NMOMENTS, &
        NBEAMS, N_USER_STREAMS, &
        DO_USER_STREAMS, USER_STREAMS, &
        MUELLER_INDEX, DMAT, &
        PI_XQP, PI_XQM, &
        PI_XUP, PI_XUM, &
        PI_X0P, &
        PI_XQM_POST, PI_XQM_PRE, &
        PI_XQP_PRE, PI_XUM_POST, &
        PI_XUP_PRE )

!  Notes
!  -----

!  This is equivalent to the Legendre setup modules in LIDORT

!---------old comments ---------------------------------
!  This needs modification for the case of Refractive atmosphere,
!  because the solar zenith angle is not constant.
!-------------------------------------------------------

!  PI-matrix setup, following recipe of Siewert (1982).
!  Single Fourier component only.

!  Tested against benchmark results in Vestrucci & Siewert (1984), Problem

!  original coding, September 2002, R. Spurr SAO
!  Multibeam SZA coding, July 2004, R. Spurr SAO
!  Coding for refractive geometry case, R. Spurr, RT Solutions, May 2005

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::           FOURIER
      LOGICAL, INTENT (IN) ::           DO_REFRACTIVE_GEOMETRY
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS
      DOUBLE PRECISION, INTENT (IN) ::  COS_SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) ::  SUN_SZA_COSINES &
          ( MAXLAYERS, MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_STREAMS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::           NMOMENTS
      INTEGER, INTENT (IN) ::           NBEAMS
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      LOGICAL, INTENT (IN) ::           DO_USER_STREAMS
      DOUBLE PRECISION, INTENT (IN) ::  USER_STREAMS  ( MAX_USER_STREAMS )

!  Indices and Dmatrix from VLIDORT_DERIVE_INPUTS

      INTEGER         , INTENT (IN)  :: MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN)  :: DMAT ( MAXSTOKES, MAXSTOKES )

!  Output (generalized spherical functions)

      DOUBLE PRECISION, INTENT (OUT) :: PI_XQP &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XQM &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XUP &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XUM &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_X0P &
          ( 0:MAXMOMENTS, MAXBEAMS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XQM_POST &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XQM_PRE &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XQP_PRE &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XUM_POST &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: PI_XUP_PRE &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )

!  Local matrices

      DOUBLE PRECISION :: XLM_DIAG ( 0:MAXMOMENTS, MAXSTOKES )
      DOUBLE PRECISION :: YLM_DIAG ( 0:MAXMOMENTS, MAXSTOKES )
      DOUBLE PRECISION :: ZLM ( 0:MAXMOMENTS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PI &
          ( 0:MAXMOMENTS, MAX_ALLSTRMS_P1, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PHI ( 0:MAXMOMENTS, 0:MAXMOMENTS )
      DOUBLE PRECISION :: PINORM ( 0:MAXMOMENTS, 0:MAXMOMENTS )

      DOUBLE PRECISION, SAVE :: DF_L ( 0:MAXMOMENTS )
      DOUBLE PRECISION, SAVE :: DF_LP1 ( 0:MAXMOMENTS )
      DOUBLE PRECISION, SAVE :: DF_2LP1 ( 0:MAXMOMENTS )
      DOUBLE PRECISION, SAVE :: DF_LSQM4 ( 0:MAXMOMENTS )
      DOUBLE PRECISION, SAVE :: DF_RT_LP3XLM1 ( 0:MAXMOMENTS )
      DOUBLE PRECISION, SAVE :: ZHELP ( 0:MAXMOMENTS )
      DOUBLE PRECISION, SAVE :: UMXSQ ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: RT_UMXSQ ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: UPXSQ_D_UMXSQ ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: PLEG20 ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: M2X_D_UMXSQ ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: XUMXSQ ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: X_RT_UMXSQ ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: PIMM_11 ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: PIMM_KM ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: X2LP1 ( 0:MAXMOMENTS, MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, SAVE :: PISIGN ( 0:MAXMOMENTS, 0:MAXMOMENTS )

!  local variables

!     INTEGER          :: LOCAL_NSTOKES
      DOUBLE PRECISION :: FAC, XSQ, DF_LPM, DF_LP1MM, ZLM_VALUE
      DOUBLE PRECISION :: HX, HY, HZ, H1, H2, RCONST_20, RCONST_21
      DOUBLE PRECISION :: PL20, RL20, PL21, RL21, TL21, XM(MAX_ALLSTRMS_P1)
      INTEGER          :: M, I, I1, L, SK1, SK2, SK3, SK4, N
      INTEGER          :: NS, NSA, M1, IB, IP

!  Set integer M = Fourier number

      M = FOURIER
      NSA = 0

!  Local NSTOKES

!      LOCAL_NSTOKES = NSTOKES
!      IF ( NSTOKES .GT. 1 ) LOCAL_NSTOKES = 4

!  total number of angles

      IF ( DO_USER_STREAMS ) THEN
        NS = NSTREAMS + N_USER_STREAMS
      ELSE
        NS = NSTREAMS
      ENDIF

!  add solar streams, depends on the use of refractive geometry

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN
!        NSA = NS + NBEAMS*NLAYERS
      ELSE
        NSA = NS + NBEAMS
      ENDIF

!  constants

      RCONST_20 = 0.25D0 * DSQRT(6.0D0)
      RCONST_21 = 2.0D0 * RCONST_20

!  Coefficient matrix PINORM (for normalization) = [ (L-m)!/(L+m)! ] ^1/
!  ---------------------------------------------------------------------

!  .. first entry = 1 / (2m)!    ---- be careful with overflow

      PHI(M,M) = ONE
      FAC = ONE
      DO L = 1, 2*M
        FAC = FAC * DBLE(L)
      ENDDO
      PHI(M,M) = PHI(M,M) / FAC

!  .. Other entries by recurrence

      DO L = M + 1, NMOMENTS
        FAC = DBLE(L-M)/DBLE(L+M)
        PHI(L,M) = PHI(L-1,M)*FAC
      ENDDO

!  .. Square root

      DO L = M , NMOMENTS
        PINORM(L,M) = DSQRT(PHI(L,M))
      ENDDO

!  Additional saved quantities (to commons in VLIDORT_PISETUP.VARS)
!  ----------------------------------------------------------------

!  Only for the first fundamental harmonic

      IF ( M .EQ. 0 ) THEN

!  Sign matrix

        DO M1 = 0, NMOMENTS
          DO L = M1, NMOMENTS
            IF (MOD((L-M1),2).EQ.0) THEN
              PISIGN(L,M1) = ONE
            ELSE
              PISIGN(L,M1) = -ONE
            ENDIF
          ENDDO
        ENDDO

!  floating point integer values

        DO L = 0, NMOMENTS
          DF_L(L)    = DBLE(L)
          DF_LP1(L)  = DBLE(L+1)
          DF_2LP1(L) = DBLE(2*L+1)
        ENDDO

        DF_LSQM4(2) = ZERO
        DO L = 3, NMOMENTS
          DF_LSQM4(L) = DSQRT(DBLE(L*L-4))
        ENDDO

        DO L = 2, NMOMENTS
          DF_RT_LP3XLM1(L) = DSQRT(DBLE((L+3)*(L-1)))
          ZHELP(L) = TWO*DF_2LP1(L)/DF_LP1(L)/DF_L(L)
        ENDDO

      ENDIF

!  local array of all streams (every Fourier component)

        DO I = 1, NSTREAMS
          XM(I) = QUAD_STREAMS(I)
        ENDDO
        IF ( DO_USER_STREAMS ) THEN
          DO I = 1, N_USER_STREAMS
            XM(I+NSTREAMS) = USER_STREAMS(I)
          ENDDO
        ENDIF

        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
!          DO N = 1, NLAYERS
!            DO IB = 1, NBEAMS
!              IP = NS + NBEAMS * (N-1) + IB
!              XM(IP) = SUN_SZA_COSINES(N,IB)
!            ENDDO
!          ENDDO
        ELSE
          DO IB = 1, NBEAMS
            IP = IB + NS
            XM(IP) = COS_SZANGLES(IB)
          ENDDO
        ENDIF

!  factors associated with stream values
!   Special case when XM(I) = ONE, factors are zeroed, and
!    those that give singularities are avoided later (see below)
!   R. Spurr and V. Natraj, 16 january 2006

      IF ( M .EQ. 0 )  THEN
        DO I = 1, NSA
          XSQ = XM(I) * XM(I)
          PLEG20(I)        = HALF*(THREE*XSQ-ONE)
          UMXSQ(I)         = ONE - XSQ
          XUMXSQ(I)        = UMXSQ(I) * XM(I)
          IF ( XM(I) .EQ. ONE ) THEN
           RT_UMXSQ(I)      = ZERO
           X_RT_UMXSQ(I)    = ZERO
           UPXSQ_D_UMXSQ(I) = ZERO
           M2X_D_UMXSQ(I)   = ZERO
          ELSE
           RT_UMXSQ(I)      = DSQRT(UMXSQ(I))
           X_RT_UMXSQ(I)    = RT_UMXSQ(I) * XM(I)
           UPXSQ_D_UMXSQ(I) = (ONE + XSQ) / UMXSQ(I)
           M2X_D_UMXSQ(I)   = - TWO * XM(I) / UMXSQ(I)
          ENDIF
          DO L = 0, NMOMENTS
            X2LP1(L,I) = XM(I) * DF_2LP1(L)
          ENDDO
        ENDDO

!  D-matrix and Mueller indices. NOW DONE in VLIDORT_DERIVE_INPUTS
!        DO SK1 = 1, MAXSTOKES
!          DO SK2 = 1, MAXSTOKES
!            MUELLER_INDEX(SK1,SK2) = MAXSTOKES*(SK1-1) + SK2
!            DMAT(SK1,SK2) = ZERO
!          ENDDO
!          IF ( SK1.GT.2) THEN
!            DMAT(SK1,SK1) = -ONE
!          ELSE
!            DMAT(SK1,SK1) = ONE
!          ENDIF
!          MUELLER_DIAGONAL_INDEX(SK1) = MUELLER_INDEX(SK1,SK1)
!        ENDDO

      ENDIF

!  XYZ matrices
!  ------------

!  Inverse XLM diagonal matrices (Siewert (1982), Eq. 35a)
!  YLM diagonal matrices (Siewert (1982), Eq. 35b)
!  ZLM matrices. (Siewert (1982), Eq. 35c)

!  .. for the azimuth-independent harmonic

      IF ( M .EQ. 0 ) THEN

        DO L = 2, NMOMENTS
          XLM_DIAG(L,1) = DF_LP1(L)
          XLM_DIAG(L,2) = DF_RT_LP3XLM1(L)
          XLM_DIAG(L,3) = XLM_DIAG(L,2)
          XLM_DIAG(L,4) = XLM_DIAG(L,1)
          DO SK1 = 1, NSTOKES
            XLM_DIAG(L,SK1) = ONE/XLM_DIAG(L,SK1)
          ENDDO
        ENDDO

        DO L = 2, NMOMENTS
          YLM_DIAG(L,1) = DF_L(L)
          YLM_DIAG(L,2) = DF_LSQM4(L)
          YLM_DIAG(L,3) = YLM_DIAG(L,2)
          YLM_DIAG(L,4) = YLM_DIAG(L,1)
        ENDDO

!  .. for the other harmonics

      ELSE

        DO L = M, NMOMENTS
          DF_LP1MM = DBLE ( L + 1 - M )
          XLM_DIAG(L,1) = DF_LP1MM
          XLM_DIAG(L,1) = ONE/XLM_DIAG(L,1)
          IF ( L .EQ. 1 ) THEN
            XLM_DIAG(L,2) = ZERO
          ELSE
            XLM_DIAG(L,2) = DF_RT_LP3XLM1(L) * DF_LP1MM / DF_LP1(L)
            XLM_DIAG(L,2) = ONE / XLM_DIAG(L,2)
          ENDIF
          XLM_DIAG(L,3) = XLM_DIAG(L,2)
          XLM_DIAG(L,4) = XLM_DIAG(L,1)
        ENDDO

        DO SK1 = 1, NSTOKES
          YLM_DIAG(M,SK1) = ZERO
        ENDDO
        DO L = M + 1, NMOMENTS
          DF_LPM = DBLE ( L + M )
          YLM_DIAG(L,1) = DF_LPM
          YLM_DIAG(L,2) = DF_LPM * DF_LSQM4(L) / DF_L(L)
          YLM_DIAG(L,3) = YLM_DIAG(L,2)
          YLM_DIAG(L,4) = YLM_DIAG(L,1)
        ENDDO

        DO L = 2, NMOMENTS
          ZLM_VALUE = DBLE(M) * ZHELP(L)
          DO SK1 = 1, NSTOKES
            DO SK2 = 1, NSTOKES
              ZLM(L,SK1,SK2) = ZERO
            ENDDO
          ENDDO
          ZLM(L,2,3) = ZLM_VALUE
          ZLM(L,3,2) = ZLM_VALUE
        ENDDO

      ENDIF

!  PI calculation
!  --------------

!  Initialise all the PI matrices for given harmonic

      DO I = 1, NSA
        DO L = M, NMOMENTS
          DO SK1 = 1, NSTOKES
            DO SK2 = 1, NSTOKES
              PI(L,I,SK1,SK2) = ZERO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!mick - alternate PI initialization
      PI = ZERO

!  M = 0 component. [Siewert (1982), Eqs (28) to (31)]

      IF ( M .EQ. 0 ) THEN

!  .. L = 0
        DO I = 1, NSA
          PI(0,I,1,1) = PI(0,I,1,1) + ONE
          PI(0,I,4,4) = PI(0,I,4,4) + ONE
        ENDDO
!  .. L = 1
        DO I = 1, NSA
          PI(1,I,1,1) = PI(1,I,1,1) + XM(I)
          PI(1,I,4,4) = PI(1,I,4,4) + XM(I)
        ENDDO
!  .. L = 2
        DO I = 1, NSA
          PL20 = PLEG20(I)
          RL20 = RCONST_20*UMXSQ(I)
          PI(2,I,1,1) = PI(2,I,1,1) + PL20
          PI(2,I,2,2) = PI(2,I,2,2) + RL20
          PI(2,I,3,3) = PI(2,I,2,2)
          PI(2,I,4,4) = PI(2,I,1,1)
        ENDDO
!  .. L > 2
        DO L = 2, NMOMENTS - 1
          DO I = 1, NSA
            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                HX = X2LP1(L,I)*PI(L,I,SK1,SK2)
                HY = YLM_DIAG(L,SK1)*PI(L-1,I,SK1,SK2)
                PI(L+1,I,SK1,SK2) = ( HX - HY ) * XLM_DIAG(L,SK1)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

!  M = 1 component. [Siewert (1982), Eqs (32) to (34), plus (35)]

      ELSE IF ( M .EQ. 1 ) THEN

!  .. L = 1
        DO I = 1, NSA
          PI(1,I,1,1) = PI(1,I,1,1) + RT_UMXSQ(I)
          PI(1,I,4,4) = PI(1,I,4,4) + RT_UMXSQ(I)
        ENDDO
!  .. L = 2
        DO I = 1, NSA
          PL21 =   THREE     * X_RT_UMXSQ(I)
          RL21 = - RCONST_21 * X_RT_UMXSQ(I)
          TL21 = + RCONST_21 * RT_UMXSQ(I)
          PI(2,I,1,1) = PI(2,I,1,1) + PL21
          PI(2,I,2,2) = PI(2,I,2,2) + RL21
          PI(2,I,2,3) = PI(2,I,2,3) + TL21
          PI(2,I,3,3) = PI(2,I,2,2)
          PI(2,I,4,4) = PI(2,I,1,1)
          PI(2,I,3,2) = PI(2,I,2,3)
        ENDDO
!  .. L > 2
        DO L = 2, NMOMENTS - 1
          DO I = 1, NSA
            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                HX = X2LP1(L,I)*PI(L,I,SK1,SK2)
                HY = YLM_DIAG(L,SK1)*PI(L-1,I,SK1,SK2)
                HZ = ZERO
                DO SK3 = 1, NSTOKES
                  HZ = HZ + ZLM(L,SK1,SK3)*PI(L,I,SK3,SK2)
                ENDDO
                PI(L+1,I,SK1,SK2) = ( HX - HY + HZ ) * XLM_DIAG(L,SK1)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

!  M > 1 components. [Siewert (1982), Eqs (36) to (38), plus (35)]
!  Limiting case of XM(I) = ONE requires special treatment to avoid NaN
!   R. Spurr and V. Natraj, 16 january 2006

      ELSE

!  .. L = M
        IF ( M .EQ. 2 ) THEN
          DO I = 1, NSA
            PIMM_11(I) =   THREE     * UMXSQ(I)
            IF ( XM(I).EQ.ONE) THEN
              PIMM_KM(I) =   RCONST_21
            ELSE
              PIMM_KM(I) =   RCONST_21 * UMXSQ(I)
            ENDIF
          ENDDO
        ELSE
          H1 = DF_2LP1(M-1)
          H2 = DF_LP1(M-1)*H1/DF_RT_LP3XLM1(M-1)
          DO I = 1, NSA
            PIMM_11(I) = H1 * RT_UMXSQ(I) * PIMM_11(I)
            IF ( XM(I).EQ.ONE) THEN
              PIMM_KM(I) = ZERO
            ELSE
              PIMM_KM(I) = H2 * RT_UMXSQ(I) * PIMM_KM(I)
            ENDIF
          ENDDO
        ENDIF
        DO I = 1, NSA
          PI(M,I,1,1) = PI(M,I,1,1) + PIMM_11(I)
          IF ( XM(I).EQ.ONE ) THEN
            PI(M,I,2,2) = PI(M,I,2,2) + PIMM_KM(I) * TWO
            PI(M,I,2,3) = PI(M,I,2,3) - PIMM_KM(I) * TWO
          ELSE
            PI(M,I,2,2) = PI(M,I,2,2) + PIMM_KM(I) * UPXSQ_D_UMXSQ(I)
            PI(M,I,2,3) = PI(M,I,2,3) + PIMM_KM(I) * M2X_D_UMXSQ(I)
          ENDIF
          PI(M,I,3,3) = PI(M,I,2,2)
          PI(M,I,4,4) = PI(M,I,1,1)
          PI(M,I,3,2) = PI(M,I,2,3)
        ENDDO
!  .. L = M + 1
        IF ( M .LT. NMOMENTS ) THEN
          DO I = 1, NSA
            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                HX = X2LP1(M,I)*PI(M,I,SK1,SK2)
                HZ = ZERO
                DO SK3 = 1, NSTOKES
                  HZ = HZ + ZLM(M,SK1,SK3)*PI(M,I,SK3,SK2)
                ENDDO
                PI(M+1,I,SK1,SK2) = ( HX + HZ ) * XLM_DIAG(M,SK1)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
!  .. L > M + 1
        DO L = M + 1, NMOMENTS - 1
          DO I = 1, NSA
            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                HX = X2LP1(L,I)*PI(L,I,SK1,SK2)
                HY = YLM_DIAG(L,SK1)*PI(L-1,I,SK1,SK2)
                HZ = ZERO
                DO SK3 = 1, NSTOKES
                  HZ = HZ + ZLM(L,SK1,SK3)*PI(L,I,SK3,SK2)
                ENDDO
                PI(L+1,I,SK1,SK2) = ( HX - HY + HZ ) * XLM_DIAG(L,SK1)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

      ENDIF

!  Normalized output.
!  ------------------

      DO L = M, NMOMENTS

!  .. at quadrature streams

        DO I = 1, NSTREAMS

          DO SK1 = 1, NSTOKES
            DO SK2 = 1, NSTOKES
              PI_XQP(L,I,SK1,SK2) = PI(L,I,SK1,SK2) * PINORM(L,M)
            ENDDO
          ENDDO

          DO SK1 = 1, NSTOKES
            DO SK2 = 1, NSTOKES

              H1 = ZERO
              H2 = ZERO
              DO SK3 = 1, NSTOKES
                H1 = H1 + DMAT(SK1,SK3) * PI_XQP(L,I,SK3,SK2)
                H2 = H2 + PI_XQP(L,I,SK1,SK3) * DMAT(SK3,SK2)
              ENDDO

              PI_XQP_PRE(L,I,SK1,SK2)  = H1
              PI_XQM_PRE(L,I,SK1,SK2)  = H1 * PISIGN(L,M)
              PI_XQM_POST(L,I,SK1,SK2) = H2 * PISIGN(L,M)

              H1 = ZERO
              DO SK3 = 1, NSTOKES
                H2 = ZERO
                DO SK4 = 1, NSTOKES
                  H2 = H2 + PI_XQP(L,I,SK3,SK4)*DMAT(SK4,SK2)
                ENDDO
                H1 = H1 + DMAT(SK1,SK3)*H2
              ENDDO
              PI_XQM(L,I,SK1,SK2)  = H1 * PISIGN(L,M)

            ENDDO
          ENDDO

        ENDDO

!  .. at positive user_defined angles

        IF ( DO_USER_STREAMS ) THEN

          DO I = 1, N_USER_STREAMS
            I1 = I + NSTREAMS

            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                PI_XUP(L,I,SK1,SK2) = PI(L,I1,SK1,SK2) * PINORM(L,M)
              ENDDO
            ENDDO

            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES

                H2 = ZERO
                H1 = ZERO
                DO SK3 = 1, NSTOKES
                  H1 = H1 + DMAT(SK1,SK3) * PI_XUP(L,I,SK3,SK2)
                  H2 = H2 + DMAT(SK3,SK2) * PI_XUP(L,I,SK1,SK3)
                ENDDO
                PI_XUP_PRE(L,I,SK1,SK2)  = H1
                PI_XUM_POST(L,I,SK1,SK2) = H2 * PISIGN(L,M)

                H1 = ZERO
                DO SK3 = 1, NSTOKES
                  H2 = ZERO
                  DO SK4 = 1, NSTOKES
                    H2 = H2 + PI_XUP(L,I,SK3,SK4)*DMAT(SK4,SK2)
                  ENDDO
                  H1 = H1 + DMAT(SK1,SK3)*H2
                ENDDO
                PI_XUM(L,I,SK1,SK2) = H1 * PISIGN(L,M)

              ENDDO
            ENDDO

          ENDDO

        ENDIF

!  .. at solar zenith angles

!  depends on use of refractive geometry

        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
!          DO N = 1, NLAYERS
!            DO IB = 1, NBEAMS
!              IP = NS + NBEAMS * (N-1) + IB
!              DO SK1 = 1, NSTOKES
!                DO SK2 = 1, NSTOKES
!                  PI_X0P(L,IB,N,SK1,SK2) =
!     &                 PI(L,IP,SK1,SK2) * PINORM(L,M)
!                ENDDO
!              ENDDO
!            ENDDO
!          ENDDO
        ELSE
          DO IB = 1, NBEAMS
            IP = IB + NS
            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                PI_X0P(L,IB,1,SK1,SK2) = &
                     PI(L,IP,SK1,SK2) * PINORM(L,M)
                DO N = 2, NLAYERS
                  PI_X0P(L,IB,N,SK1,SK2)= PI_X0P(L,IB,1,SK1,SK2)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  end loop moments

      ENDDO

!  debug

!      write(77,*)'Fourier',M
!      DO L = M, NMOMENTS
!       write(77,*)'Moment ',L
!       DO SK1 = 1, 4
!        WRITE(77,'(I3,1p4e15.6)')SK1,(PI_XUP(L,4,SK1,SK2),SK2=1,4)
!       ENDDO
!      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_PIMATRIX_SETUP

      END MODULE vlidort_solutions

