C ###############################################################
C #                                                             #
C #                    THE VECTOR LIDORT MODEL                  #
C #                                                             #
C #  (Vector LInearized Discrete Ordinate Radiative Transfer)   #
C #   -      --         -        -        -         -           #
C #                                                             #
C ###############################################################

C ###############################################################
C #                                                             #
C #  Author :      Robert. J. D. Spurr                          #
C #                                                             #
C #  Address :      RT Solutions, inc.                          #
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
C    #   This Version of VLIDORT comes with a GNU-style  #
C    #   license. Please read the license carefully.     #
C    #                                                   #
C    #####################################################

C ###############################################################
C #                                                             #
C # Subroutines in this Module                                  #
C #                                                             #
C #              VLIDORT_QHOM_SOLUTION                          #
C #              VLIDORT_UHOM_SOLUTION                          #
C #              VLIDORT_QBEAM_SOLUTION                         #
C #              VLIDORT_UBEAM_SOLUTION                         #
C #              VLIDORT_PIMATRIX_SETUP                         #
C #                                                             #
C ###############################################################

      SUBROUTINE VLIDORT_QHOM_SOLUTION
     I    ( GIVEN_LAYER, FOURIER, 
     O      STATUS )

C  Numerical solution of Eigenproblem.

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include file of Set up stuff (input/output to this module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  include file of main solution variables (output to this module)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index and Fourier number (inputs)

      INTEGER             GIVEN_LAYER
      INTEGER             FOURIER

C  output status

      INTEGER             STATUS

C  Local variables
C  ---------------

C  local matrix for eigenvalue computation

      DOUBLE PRECISION
     &      EIGENMAT(MAXSTRMSTKS,MAXSTRMSTKS)

C  output from Eigenpackage module ASYMTX
C    Following output now stored in commons:
C      DOUBLE PRECISION REAL_KSQ(MAXSTRMSTKS)
C      DOUBLE PRECISION RITE_EVEC(MAXSTRMSTKS,MAXSTRMSTKS)

C  Additional diagnostic output from ASYMTX, and Work space.
C    Dimension of scratch array WK increased to 4 times the basic dimension
C    R. Spurr, RTS, December 2005: More work space needed to avoid seg-faults.
C    Internal dimensioning in ASYMTX is the same.

c      DOUBLE PRECISION   WK(16*MAXSTRMSTKS)
      DOUBLE PRECISION   WK(4*MAXSTRMSTKS)
      INTEGER            IER
      LOGICAL            ASYMTX_FAILURE
      CHARACTER*3        CN, CI

C  output for Eigenpackage DGEEV ( LAPACK)
C    Following now stored in commons...............
C      DOUBLE PRECISION REAL_KSQ(MAXSTRMSTKS),IMAG_KSQ(MAXSTRMSTKS)
C      DOUBLE PRECISION LEFT_EVEC(MAXSTRMSTKS,MAXSTRMSTKS)
C      DOUBLE PRECISION RITE_EVEC(MAXSTRMSTKS,MAXSTRMSTKS)

C  Additional diagnostic output and Work space from DGEEV.
C    Dimension LWORK increased from 4 to 64 times the basic dimension
C    J. Kujanpaa, FMI, August 2005: More work space for DGEEV needed.

      INTEGER          LWORK
c      PARAMETER       (LWORK = 4*MAXSTRMSTKS)
      PARAMETER       (LWORK = 64*MAXSTRMSTKS)
      DOUBLE PRECISION WORK(LWORK)
      INTEGER          DGEEV_INFO

C  Exception handling

      CHARACTER*70     MESSAGE, TRACE

C  optimization of Sum and Difference matrices SAB and DAB
C  ...requires intermediate storage of these two arrays.
C   Timing tests by J. Kujanpaa, FMI, August 2005.

      DOUBLE PRECISION
     &     H2PARR(MAXSTOKES,0:MAXMOMENTS,MAXSTOKES,MAXSTREAMS),
     &     H2MARR(MAXSTOKES,0:MAXMOMENTS,MAXSTOKES,MAXSTREAMS)

C  Miscellaneous local variables

      INTEGER          I, J, I1, L, N, NA, UTA, UT, NMINST
      INTEGER          M, AA, AA1, K, K_R, K_C, K_AA, KO1, K0, K1, K2
      INTEGER          O1, O2, O3, O4, IR, IROW, JC, JCOL
      DOUBLE PRECISION DP, DM, FAC, SC_R, XINV, HELP
      DOUBLE PRECISION H1P, H1M, NORM_R, TAU_DN, TAU_UP
      DOUBLE PRECISION FWD_H1, FWD_H2
      DOUBLE PRECISION XSQ_R, YSQ_R, RKSQ, IKSQ, KMUT, KVAL
      DOUBLE PRECISION CARG, SARG, ARG, MOD_KVAL, H1, H2
      DOUBLE PRECISION FWD_R, FWD_I, H2P, H2M, FPD, FND

C  original complex variables were:
C      COMPLEX*16      FWD_H1_C, FWD_H2_C, ADJ_H1_C, ADJ_H2_C
C      COMPLEX*16      HELP_C, SC_C, H1_C, H2_C

C  replacement for complex variables

      DOUBLE PRECISION HELP_CR, HELP_CI
      DOUBLE PRECISION FWD_H1_CR, FWD_H1_CI
      DOUBLE PRECISION FWD_H2_CR, FWD_H2_CI
      DOUBLE PRECISION SC_CR, SC_CI, MODKRT
      DOUBLE PRECISION FPD1, FPD2, FND1, FND2

C  Code start
C  ----------

C  Initialization

      KO1 = 0

C  initialise status

      STATUS = VLIDORT_SUCCESS

C  Layer and Fourier

      N = GIVEN_LAYER
      M = FOURIER

C  Flipper number

      NMINST = MIN(2,NSTOKES)

C  If there is no scattering in this layer, we can set
C  the solution vectors directly and move on.
C    Just a transmitttance case.

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

C  Scattering solutions
C  ====================

C  Construct Eigenmatrix
C  ---------------------

C  zero the eigensolution matrix
C   Not needed. Redundancy test by J. Kujanpaa, FMI, August 2005.
c      DO I = 1, NSTKS_NSTRMS
c        DO J = 1, NSTKS_NSTRMS
c          EIGENMAT(I,J) = ZERO
c        ENDDO
c      ENDDO

C  Develop Sum and Difference matrices, part I
C   Introduction of holding arrays by J. Kujanpaa, FMI, August 2005

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

C  debug code for hparr, hmarr

c      DO J = 1, NSTREAMS
c       DO L = M, NMOMENTS
c         DO O2 = 1, min(NSTOKES,3)
c          DO O3 = 1, min(NSTOKES,3)
c           if (nstokes.eq.3.and.n.eq.6.and.m.eq.0)
c     *         write(87,'(4i3,1p2e17.8)')
c     *        J,L,O2,O3,H2PARR(O3,L,O2,J),H2MARR(O3,L,O2,J)
c           if (nstokes.eq.4.and.n.eq.6.and.m.eq.0)
c     *         write(88,'(4i3,1p2e17.8)')
c     *        J,L,O2,O3,H2PARR(O3,L,O2,J),H2MARR(O3,L,O2,J)
c          ENDDO
c         ENDDO
c        ENDDO
c      ENDDO
c      if (n.eq.6.and.m.eq.1)pause 'hp 6'

C  Develop Sum and Difference matrices, part II of calculation
C   Use of holding arrays by J. Kujanpaa, FMI, August 2005

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

C  Note. In some circumstances, H2PARR = H2MARR
C        Then DP = Dm in the above code, and DAB is Diagonal.
C        Then there are pair-degenerate eigenvalues, equal to the
C        discrete ordinate directions.

C  Compute Eigenmatrix 

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

C  save Eigenmatrix (original is destroyed, want to use later)
C    Looping order changed by J. Kujanpaa, FMI, August 2005

      DO J = 1, NSTKS_NSTRMS
        DO I = 1, NSTKS_NSTRMS
          EIGENMAT_SAVE(I,J,N) = EIGENMAT(I,J)
        ENDDO
      ENDDO

c      if (n.eq.24) then
c      do j = 1, 9       
c      if(nstokes.eq.3)write(33,'(i4,1p12e11.3)')J,(EIGENMAT(i,J),i=1,9)
c      if(nstokes.eq.4)write(34,'(i4,1p12e11.3)')J,(EIGENMAT(i,J),i=1,12)
c      enddo
c      pause
c      endif

C  Debug for testing new eigensolver
C  7 August 2006. Result --> NEW LAPACK ROUTINES NOT FASTER
      
c      IF ( N.EQ.1) then
c        OPEN(88,file='esolver.inp',status='unknown')
c      ENDIF      
c      write(88,'(2i5,l2)')n,NSTKS_NSTRMS,DO_REAL_EIGENSOLVER(M,N)
c      DO I = 1, NSTKS_NSTRMS
c        write(88,'(1p32e20.10)')(eigenmat(i,j),j=1,NSTKS_NSTRMS)
c      enddo
c      if ( n.eq.nlayers)then
c        close(88)
c        pause'end debug'
c      ENDIF

C  Eigensolver package
C  -------------------

      IF ( DO_REAL_EIGENSOLVER(M,N) ) THEN

C  Here is the DISORT ASYMTX package, as used in LIDORT scalar codes
C    This is sufficient for Real Symmetric Eigenmatrices.

       CALL  ASYMTX
     I      ( EIGENMAT, NSTKS_NSTRMS, MAXSTRMSTKS, MAXSTRMSTKS,
     O        RITE_EVEC, REAL_KSQ, IER, WK,
     O        MESSAGE, ASYMTX_FAILURE )

C  Exception handling

       IF ( ASYMTX_FAILURE .OR. IER.GT.0 ) THEN
        WRITE(CN,'(I3)')N
        TRACE = 'ASYMTX call, layer '//CN//' in VLIDORT_QHOM_SOLUTION'
        STATUS = VLIDORT_SERIOUS
        IF ( IER .GT. 0 ) THEN
         WRITE(CI,'(I3)')IER
         MESSAGE  = 'eigenvalue '//CI//' has not converged'
        ENDIF
        CALL VLIDORT_ERROR_TRACE ( MESSAGE, TRACE, STATUS )
        RETURN
       ENDIF

C  renormalize - this is vital for ASMTX output

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

C  Complex roots : Must use DGEEV (LAPACK module)
C  DGEEV is 1.7-2.0 times slower than ASYMTX because look for complex roots.
C  Also first call to DGEEV is 20 times slower.

       CALL DGEEV
     I      ( 'V', 'V', NSTKS_NSTRMS, EIGENMAT,
     O        MAXSTRMSTKS, REAL_KSQ, IMAG_KSQ,
     O        LEFT_EVEC, MAXSTRMSTKS, RITE_EVEC, MAXSTRMSTKS,
     W        WORK, LWORK, DGEEV_INFO )

C  Exception handling

       IF ( DGEEV_INFO .NE. 0 ) THEN
        WRITE(CN,'(I3)')N
        TRACE = 'DGEEV call, layer '//CN//' in VLIDORT_QHOM_SOLUTION'
        STATUS = VLIDORT_SERIOUS
        IF ( DGEEV_INFO .LT. 0 ) THEN
         WRITE(CI,'(I3)')IABS(DGEEV_INFO)
         MESSAGE='Argument # '//CI//' had an illegal value'
        ELSE
         MESSAGE='QR algorithm in DGEEV failed to compute all e-values'
        ENDIF
        CALL VLIDORT_ERROR_TRACE ( MESSAGE, TRACE, STATUS )
        RETURN
       ENDIF

      ENDIF

C  Sort real and complex eigenvalues; get masks
C  --------------------------------------------

C  Initialise counts

       K_R  = 0
       K_C  = 0
       K_AA = 0

C  all eigenvalues are real if ASYMTX is used

      IF ( DO_REAL_EIGENSOLVER(M,N) ) THEN
       DO AA = 1, NSTKS_NSTRMS
        K_R = K_R + 1
        EIGENMASK_R(K_R) = AA
       ENDDO
      ENDIF

C  Degeneracy mask

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

C  For the DGEEV eigensolver--
C    real eigenvalues if IMAG_KSQ = 0, otherwise complex in pairs

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

C  save number of eigenvalues

      K_REAL(N)    = K_R
      K_COMPLEX(N) = K_C

C  Real Solutions
C  --------------

C  start loop over real eigenvalues

      DO K = 1, K_REAL(N)

C  get eigenvector entry

        AA = EIGENMASK_R(K)

C  eigenvalue and separation constants (save the eigenvalues)

        KEIGEN(K,N) = DSQRT(REAL_KSQ(AA))
        SC_R    = ONE / KEIGEN(K,N)

C  set Right eigenvector norms
C    This piece of code should not be necessary as 
C    (a) DGEEV norms are all = 1 upon output
C    (b) ASMTYX output has been renormalized to unity

        NORM_R = ZERO
        DO I = 1, NSTKS_NSTRMS
          NORM_R = NORM_R + RITE_EVEC(I,AA)*RITE_EVEC(I,AA)
        ENDDO
        NORM_R = DSQRT(NORM_R)

C  Forward Sum-vector = normalized Right eigenvector

        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            FWD_SUMVEC(I,O1,K) = RITE_EVEC(IROW,AA)/NORM_R
          ENDDO
        ENDDO

C  Debug (important)

c        if ( do_debug_write ) then
c         if (m.eq.1.and.n.gt.4) then
c          write(27,'(2I3,1pe20.12)') N,K,KEIGEN(K,N)
c          DO I = 1, NSTREAMS
c           DO O1 = 1, NSTOKES
c            write(27,'(4i3,1pe20.12)')N,k,i,o1,FWD_SUMVEC(I,O1,K)
c           ENDDO
c          ENDDO
c         ENDIF
c        endif

C  Find Forward difference vectors (Siewert's notation)

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

C  assign Forward solutions PHI (Siewert's notation)
C    --->  first N are "DOWN", last N are "UP" (streams)
C    --->  Use symmetry properties to set -ve eigensolutions

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

C  Older code: Use wasteful additional arrays...
C    solution assignation (only for BVP)

c        DO I = 1, NSTREAMS
c          XINV = HALF
c          I1 = I + NSTREAMS
c          FWD_XPOS(I,O1,K,N)  =
c     &      XINV * ( FWD_SUMVEC(I,O1,K) + FWD_DIFVEC(I,O1,K) )
c          FWD_XPOS(I1,O1,K,N)  =
c     &      XINV * ( FWD_SUMVEC(I,O1,K) - FWD_DIFVEC(I,O1,K) )
c          FWD_XNEG(I1,O1,K,N) = FWD_XPOS(I,O1,K,N)
c          FWD_XNEG(I,O1,K,N)  = FWD_XPOS(I1,O1,K,N)
c         ENDDO
c        ENDDO
c        DO I = 1, NSTREAMS
c          I1 = I + NSTREAMS
c          DO O1 = 1, NSTOKES
c            SOLA_XPOS(I,O1,K,N)  = FWD_XPOS(I,O1,K,N)
c            SOLB_XNEG(I,O1,K,N)  = FWD_XNEG(I,O1,K,N)
c          ENDDO
c          DO O1 = 1, NMINST
c            SOLA_XPOS(I1,O1,K,N)  = FWD_XNEG(I,O1,K,N)
c            SOLB_XNEG(I1,O1,K,N)  = FWD_XPOS(I,O1,K,N)
c          ENDDO
c          DO O1 = 3, NSTOKES
c            SOLA_XPOS(I1,O1,K,N)  = - FWD_XNEG(I,O1,K,N)
c            SOLB_XNEG(I1,O1,K,N)  = - FWD_XPOS(I,O1,K,N)
c          ENDDO
c        ENDDO

C  end loop over real eigenvalues

      ENDDO

C  Complex solutions
C  -----------------

C  Offsets

      KO1 = K_REAL(N) + 1

C  start loop over complex eigenvalues

      DO K = 1, K_COMPLEX(N)

C  Bookkeeping indices

        K0 = 2 * ( K - 1 )
        K1 = KO1 + K0
        K2 = K1  + 1

C  get eigenvector entry

        AA  = EIGENMASK_C(K)
        AA1 = AA + 1

C  set eigenvalue and separation constants for one of pair of Conjugates
C   ---> (the complex conjugate follows from this information)

        RKSQ = REAL_KSQ(AA) * REAL_KSQ(AA)
        IKSQ = IMAG_KSQ(AA) * IMAG_KSQ(AA)
        MOD_KVAL = DSQRT(RKSQ + IKSQ)
        ARG = HALF * DATAN ( IMAG_KSQ(AA) / REAL_KSQ(AA) )
        CARG = DCOS(ARG)
        SARG = DSIN(ARG)

C  original complex --------------------------------------------------
C        KEIGEN_C(K,N) = DSQRT(MOD_KVAL) * CMPLX ( CARG, SARG )
C        SC_C  = ONE_C / KEIGEN_C(K,N)
C  -------------------------------------------------------------------

C  replacement complex

        KEIGEN_CSQ(K) = MOD_KVAL
        MODKRT        = DSQRT(MOD_KVAL)
        KEIGEN(K1,N)  = MODKRT * CARG
        KEIGEN(K2,N)  = MODKRT * SARG
        SC_CR =   KEIGEN(K1,N) / MOD_KVAL
        SC_CI = - KEIGEN(K2,N) / MOD_KVAL

C  set Right eigenvector norms
C   ---> Follows convention of way that eigenvectors are fomatted
C        upon output from the LAPACK module DGEEV.

        NORM_R = ZERO
        DO I = 1, NSTKS_NSTRMS
          XSQ_R = RITE_EVEC(I,AA)  * RITE_EVEC(I,AA)
          YSQ_R = RITE_EVEC(I,AA1) * RITE_EVEC(I,AA1)
          NORM_R = NORM_R + XSQ_R + YSQ_R
        ENDDO
        NORM_R = DSQRT(NORM_R)

C  Forward Sum-vector = normalized Right eigenvector

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

C  Debug (important)

c        if ( do_debug_write ) then
c         if (m.eq.0.and.n.eq.6)then
c          write(90,'(2i3,1p2e20.10)')
c     &            m,k,KEIGEN(K1,N),KEIGEN(K2,N)
c          DO I = 1, NSTREAMS
c           DO O1 = 1, NSTOKES
c             write(90,'(4i3,1p4e20.10)')
c     &      m,k,i,o1,FWD_SUMVEC(I,O1,K1),FWD_SUMVEC(I,O1,K2)
c           ENDDO
c          ENDDO
c         ENDIF
c        endif

C  Find Forward difference vectors (Siewert's notation)

        DO I = 1, NSTREAMS
         DO O1 = 1, NSTOKES
          FWD_H1_CR = ZERO
          FWD_H1_CI = ZERO
          DO J = 1, NSTREAMS
           FWD_H2_CR = ZERO
           FWD_H2_CI = ZERO
           DO O2 = 1, NSTOKES
            FWD_H2_CR = FWD_H2_CR +
     &           SAB(I,J,O1,O2,N) * FWD_SUMVEC(J,O2,K1)
            FWD_H2_CI = FWD_H2_CI +
     &           SAB(I,J,O1,O2,N) * FWD_SUMVEC(J,O2,K2)
           ENDDO
           FWD_H1_CR = FWD_H1_CR + FWD_H2_CR
           FWD_H1_CI = FWD_H1_CI + FWD_H2_CI
          ENDDO
          FWD_DIFVEC(I,O1,K1) = FWD_H1_CR*SC_CR - FWD_H1_CI*SC_CI
          FWD_DIFVEC(I,O1,K2) = FWD_H1_CR*SC_CI + FWD_H1_CI*SC_CR
         ENDDO
        ENDDO

C  assign Forward solutions PHI (Siewert's notation)
C    --->  first N are "DOWN", last N are "UP" (streams)
C    --->  Use symmetry properties to set -ve eigensolutions

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

C  Older code: Used additional wasteful holding arrays
C   solution assignation (only for BVP)

c        DO I = 1, NSTREAMS
c          XINV = HALF
c          I1 = I + NSTREAMS
c          DO O1 = 1, NSTOKES
c            FWD_XPOS(I,O1,K1,N)  = XINV * 
c     &         ( FWD_SUMVEC(I,O1,K1) + FWD_DIFVEC(I,O1,K1) )
c            FWD_XPOS(I,O1,K2,N)  = XINV * 
c     &         ( FWD_SUMVEC(I,O1,K2) + FWD_DIFVEC(I,O1,K2) )
c            FWD_XPOS(I1,O1,K1,N) = XINV * 
c     &         ( FWD_SUMVEC(I,O1,K1) - FWD_DIFVEC(I,O1,K1) )
c            FWD_XPOS(I1,O1,K2,N) = XINV * 
c     &         ( FWD_SUMVEC(I,O1,K2) - FWD_DIFVEC(I,O1,K2) )
c            FWD_XNEG(I1,O1,K1,N) = FWD_XPOS(I,O1,K1,N)
c            FWD_XNEG(I,O1,K1,N)  = FWD_XPOS(I1,O1,K1,N)
c            FWD_XNEG(I1,O1,K2,N) = FWD_XPOS(I,O1,K2,N)
c            FWD_XNEG(I,O1,K2,N)  = FWD_XPOS(I1,O1,K2,N)
c          ENDDO
c        ENDDO
c        DO I = 1, NSTREAMS
c          I1 = I + NSTREAMS
c          DO O1 = 1, NSTOKES
c            SOLA_XPOS(I,O1,K1,N)  = FWD_XPOS(I,O1,K1,N)
c            SOLB_XNEG(I,O1,K1,N)  = FWD_XNEG(I,O1,K1,N)
c            SOLA_XPOS(I,O1,K2,N)  = FWD_XPOS(I,O1,K2,N)
c            SOLB_XNEG(I,O1,K2,N)  = FWD_XNEG(I,O1,K2,N)
c          ENDDO
c          DO O1 = 1, NMINST
c            SOLA_XPOS(I1,O1,K1,N)  = FWD_XNEG(I,O1,K1,N)
c            SOLB_XNEG(I1,O1,K1,N)  = FWD_XPOS(I,O1,K1,N)
c            SOLA_XPOS(I1,O1,K2,N)  = FWD_XNEG(I,O1,K2,N)
c            SOLB_XNEG(I1,O1,K2,N)  = FWD_XPOS(I,O1,K2,N)
c          ENDDO
c          DO O1 = 3, NSTOKES
c            SOLA_XPOS(I1,O1,K1,N)  = -FWD_XNEG(I,O1,K1,N)
c            SOLB_XNEG(I1,O1,K1,N)  = -FWD_XPOS(I,O1,K1,N)
c            SOLA_XPOS(I1,O1,K2,N)  = -FWD_XNEG(I,O1,K2,N)
c            SOLB_XNEG(I1,O1,K2,N)  = -FWD_XPOS(I,O1,K2,N)
c          ENDDO
c        ENDDO

C  Try to see if the solutions are the same
C    Multilayer case - checks out 17 December.
c      IF ( M.le.1 ) THEN
c      IF ( N.EQ.2 ) THEN
c       DO I = 1, 2*NSTREAMS
c        WRITE(M+1,'(2i3,1p4e16.6)')K,I,(SOLA_XPOS(I,O1,K1,N),O1=1,4)
c        WRITE(M+1,'(2i3,1p4e16.6)')k,I,(SOLA_XPOS(I,O1,K2,N),O1=1,4)
c        WRITE(M+1,'(2i3,1p4e16.6)')k,I,(SOLB_XNEG(I,O1,K1,N),O1=1,4)
c        WRITE(M+1,'(2i3,1p4e16.6)')k,I,(SOLB_XNEG(I,O1,K2,N),O1=1,4)
c       ENDDO
c      ENDIF
c      ENDIF

C  end eigenvalue loop

      ENDDO

C  control point for avoiding eigensolver
C  ======================================

 3456 CONTINUE

C  Eigenstream transmittance factors for whole layer
C  -------------------------------------------------

C  When the solution saving option is set, then if there is no
C  scattering in this layer, then linearized transmittances are
C  linearized discrete ordinate transmittances.

      IF ( DO_SOLUTION_SAVING .AND.
     &             .NOT.DO_LAYER_SCATTERING(M,N) ) THEN

        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            T_DELT_EIGEN(IROW,N) = T_DELT_DISORDS(I,N)
          ENDDO
        ENDDO

C  Otherwise compute them as normal
 
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

c      if ( m.eq.0.and.n.eq.20 ) then
c       DO K = 1, min(6,K_REAL(N))
c         WRITE(*,'(2i4,1p6e24.12)')m,k,T_DELT_EIGEN(K,N)
c       ENDDO
c      ENDIF

C  Eigenstream transmittance factors for partial layers
C  ----------------------------------------------------

C Loop over user optical depths

      DO UTA = 1, N_USER_LEVELS
       IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
        UT = PARTLAYERS_OUTINDEX(UTA)
        NA = PARTLAYERS_LAYERIDX(UT)

C  If the layer of occurence = given layer then do the calculation

         IF ( NA .EQ. N ) THEN

C  When the solution saving option is set, then if there is no
C  scattering in this layer, then transmittances are just the
C  discrete ordinate transmittances.

          IF ( DO_SOLUTION_SAVING .AND.
     &         .NOT.DO_LAYER_SCATTERING(M,N) ) THEN

           DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
             IROW = IR + O1
             T_UTDN_EIGEN(IROW,UT) = T_DISORDS_UTDN(I,UT)
             T_UTUP_EIGEN(IROW,UT) = T_DISORDS_UTUP(I,UT)
            ENDDO
           ENDDO

C  Otherwise compute them as Eigenstream transmittance factors

          ELSE

C  partial layer optical depths

           TAU_DN = PARTAU_VERT(UT)
           TAU_UP = DELTAU_VERT(N) - TAU_DN

C  Real eigenvalues

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

C  Complex eigenvalues

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

C  End clause to use Discrete ordinates or not

          ENDIF

C  end loop over off-grid optical depths

         ENDIF
        ENDIF
      ENDDO

C  Debug

c        if ( m.lt.3.and.n.gt.0) then
c        write(97,'(2i5)')M,N
c        DO K = 1, NSTREAMS
c          WRITE(97,'(I3,1p6e17.9)')K,KEIGEN(K,N),
c     *          (SOLA_XPOS(I,1,K,N),I=1,NSTREAMS,2)
c        ENDDO
c        endif
cc        if (n.eq.26)pause
c      if (n.eq.6) pause'end hom'

C  Finish

      RETURN
      END

C

      SUBROUTINE VLIDORT_UHOM_SOLUTION
     I    ( GIVEN_LAYER, FOURIER )

C  Small-number Taylor series expansions, added 30 October 2007

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include file of setup stuff (input to this module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  include file of solution variables (output to this module)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'

C  include file of Multiplier coefficients

      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index and Fourier number (inputs)

      INTEGER          GIVEN_LAYER
      INTEGER          FOURIER

C  Local variables
C  ---------------

      INTEGER          UM, J, L, N, M, K, O1, O2, KO1, K0, K1, K2

C  for real solutions

      DOUBLE PRECISION SM, SN, SP, SNS, SPS, A5
      DOUBLE PRECISION RHO_P, RHO_M, SECMUI, SUMG

C  for complex solutions

      DOUBLE PRECISION RHO_P_CR, RHO_M_CR, MODULUS_P, MODULUS_M
      DOUBLE PRECISION A_USER, A_HELP, B_USER, B_HELP, HPC, HMC
      DOUBLE PRECISION SMCR, SNCR, SPCR, SNSCR, SPSCR, SGCR
      DOUBLE PRECISION SMCI, SNCI, SPCI, SNSCI, SPSCI, SGCI

C  Layer and Fourier

      N = GIVEN_LAYER
      M = FOURIER

C  Initialization

      K   = 0

C  Zeta constants
C  ==============

C  Real Zeta constants (always required)

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

C  Complex Zeta constants

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

C  If there is no scattering avoid these solutions
C  Warning - Solutions are zeroed in the scalar code before returning

      IF ( .NOT. DO_LAYER_SCATTERING(M,N) ) RETURN

C  User defined solutions (Real)
C  -----------------------------

      DO K = 1, K_REAL(N)

C  For each moment, do inner sum over computational angles
C  for the positive and negative eigenvectors.
C    Gvec = STOKES. Eq 34. GAUx = Greekmat.Gvec

       DO L = M, NMOMENTS
        DO O1 = 1, NSTOKES
         SPS = ZERO
         DO J = 1, NSTREAMS
          SP = ZERO
          SM = ZERO
          A5 = QUAD_HALFWTS(J)
          DO O2 = 1, NSTOKES
c  older   SP = SP + PI_XQP    (L,J,O1,O2)*A5(J)*FWD_XPOS(J,O2,K,N)
c  older   SM = SM + PI_XQM_PRE(L,J,O1,O2)*A5(J)*FWD_XNEG(J,O2,K,N)
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

C  Now sum over all harmonic contributions (downwelling)

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

C  Now sum over all harmonic contributions (Upwelling)

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

C  Debug (important)

        if ( do_debug_write.and.do_fdtest) then
         DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            write(94,'(4i3,1p4e20.10)')m,k,um,o1,
     &         UHOM_UPUP(UM,O1,K,N),
     &         UHOM_UPDN(UM,O1,K,N),
     &         UHOM_DNUP(UM,O1,K,N),
     &         UHOM_DNDN(UM,O1,K,N)
           enddo
         enddo
        endif

C  end loop over real eigensolutions

      ENDDO

C  User defined solutions (Complex)
C  -------------------------------

C  offset

      KO1 = K_REAL(N) + 1

      DO K = 1, K_COMPLEX(N)

       K0 = 2 * K - 2
       K1 = KO1 + K0
       K2 = K1  + 1

C  For each moment, do inner sum over computational angles
C  for the positive and negative eigenvectors.
C    Gvec = STOKES. Eq 34. GAUx = Greekmat.Gvec

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
c           SPCR = SPCR + HPC * FWD_XPOS(J,O2,K1,N)         ! Older
c           SPCI = SPCI + HPC * FWD_XPOS(J,O2,K2,N)         ! Older
           SPCR = SPCR + HPC * SOLA_XPOS(J,O2,K1,N)
           SPCI = SPCI + HPC * SOLA_XPOS(J,O2,K2,N)
           HMC  =  PI_XQM_PRE(L,J,O1,O2)*QUAD_HALFWTS(J)
c           SMCR = SMCR + HMC * FWD_XNEG(J,O2,K1,N)         ! Older
c           SMCI = SMCI + HMC * FWD_XNEG(J,O2,K2,N)         ! Older
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

C  Now sum over all harmonic contributions (downwelling)

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

C  Now sum over all harmonic contributions (upwelling)

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

C  Debug (important)
c        if ( do_debug_write.and.do_fdtest) then
c         DO UM = LOCAL_UM_START, N_USER_STREAMS
c          DO O1 = 1, NSTOKES
c            write(96,'(4i3,1p4e20.10)')m,k,um,o1,
c     &         UHOM_UPUP(UM,O1,K2,N),
c     &         UHOM_UPDN(UM,O1,K2,N),
c     &         UHOM_DNUP(UM,O1,K2,N),
c     &         UHOM_DNDN(UM,O1,K2,N)
c           enddo
c         enddo
c        endif

C  Try to see if the solutions are the same
C    Multilayer case - checks out 17 December.
c        IF ( M.eq.1 ) THEN
c        IF ( N.EQ.2 ) THEN
c         DO UM = 1, N_USER_STREAMS
c          WRITE(N,'(2i3,1p4e16.6)')K,UM,(UHOM_UPUP(UM,O1,K1,N),O1=1,4)
c          WRITE(N,'(2i3,1p4e16.6)')K,UM,(UHOM_UPUP(UM,O1,K2,N),O1=1,4)
c          WRITE(N,'(2i3,1p4e16.6)')K,UM,(UHOM_UPDN(UM,O1,K1,N),O1=1,4)
c          WRITE(N,'(2i3,1p4e16.6)')K,UM,(UHOM_UPDN(UM,O1,K2,N),O1=1,4)
c          WRITE(N,'(2i3,1p4e16.6)')K,UM,(UHOM_DNUP(UM,O1,K1,N),O1=1,4)
c          WRITE(N,'(2i3,1p4e16.6)')K,UM,(UHOM_DNUP(UM,O1,K2,N),O1=1,4)
c          WRITE(N,'(2i3,1p4e16.6)')K,UM,(UHOM_DNDN(UM,O1,K1,N),O1=1,4)
c          WRITE(N,'(2i3,1p4e16.6)')K,UM,(UHOM_DNDN(UM,O1,K2,N),O1=1,4)
c         ENDDO
c         ENDIF
c        ENDIF

C  end loop over complex eigensolutions

      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE VLIDORT_QBEAM_SOLUTION
     I    ( GIVEN_LAYER,
     I      FOURIER,
     I      IBEAM,
     O      STATUS )

C  This is the classical Chandrasekhar beam solution.
C  ( plane parallel or average secant only)
C  Linear Matrix algebra.

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include file of Set up stuff (input/output to this module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  include file of main solution variables (output to this module)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index and Fourier number (inputs)

      INTEGER             GIVEN_LAYER
      INTEGER             FOURIER
      INTEGER             IBEAM

C  output status

      INTEGER             STATUS

C  Local variables
C  ---------------

C  help variables

      CHARACTER*70     MAIL
      INTEGER          I, J, I1, L, N, M, INFO
      INTEGER          O1, O2, O3, IR, JR, IROW, JROW
      DOUBLE PRECISION SUM, GSUM, INV_X0SQ, SECBAR
      DOUBLE PRECISION WSUM, WDIF, TRANS1, TRANS2, WVEC_D, WVEC_U
      DOUBLE PRECISION TPA, TMA, S_TPA, S_TMA, S_HELP, HELP

C  debug local variables
C    removed 2 July 2004. See Older version, labeled _debug

C  Initial section
C  ---------------

C  initialise indices

      N = GIVEN_LAYER
      M = FOURIER

C  No particular solution beyond the cutoff layer
C  OR No scattering in this layer then:
C    ---> Zero the boundary layer values and exit

      IF ( N .GT. LAYER_PIS_CUTOFF(IBEAM) .OR.
     &       .NOT. DO_LAYER_SCATTERING(M,N) ) THEN
        DO I = 1, NSTREAMS_2
          DO O1 = 1, NSTOKES
            WUPPER(I,O1,N) = ZERO
            WLOWER(I,O1,N) = ZERO
          ENDDO
        ENDDO
        RETURN
      ENDIF

C  set local values

      SECBAR   = AVERAGE_SECANT(N,IBEAM)
      INV_X0SQ = SECBAR * SECBAR

C  Initialise matrix and column vector for solution
C  Changed at FMI 28 september 2004

      DO J = 1, NSTKS_NSTRMS
        DO I = 1, NSTKS_NSTRMS
          QMAT_SAVE(I,J) = ZERO
        ENDDO
        QVEC_SAVE(J) = ZERO
      ENDDO

C  Set up sum and difference vectors for Beam source terms
C  ( sum vector may be required again in linearization )
C  Auxiliary matrix for Q functions

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
          QSUMVEC_SAVE(I,O1) =  ( S_TPA + S_TMA ) / QUAD_STREAMS(I)
          QDIFVEC_SAVE(I,O1) =  ( S_TPA - S_TMA ) / QUAD_STREAMS(I)
        ENDDO
      ENDDO

C  solution matrix for the reduced problem
C  ( matrix should be saved in the LU decomposition form)
C  Changed at FMI 28 september

      DO J = 1, NSTKS_NSTRMS
        DO I = 1, NSTKS_NSTRMS
          QMAT_SAVE(I,J) = EIGENMAT_SAVE(I,J,N)
        ENDDO
      ENDDO

      DO I = 1, NSTKS_NSTRMS
        QMAT_SAVE(I,I) = QMAT_SAVE(I,I) - INV_X0SQ
      ENDDO

C  RHS vector for the reduced problem
C  ( this vector will be the answer after the linear algebra solution,
C    and may be needed again if there is linearization )

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

C  debug

C      write(*,'(i3)')N
c      DO I = 1, NSTREAMS
c       write(*,'(i3,1p5e11.3)')I,QVEC_SAVE(I),(QMAT_SAVE(I,J),J=1,4)
c       ENDDO
c      if (n.eq.2)pause'qsetup'

C  L-U decomposition of the solution matrix

      CALL DGETRF
     1      ( NSTKS_NSTRMS, NSTKS_NSTRMS,
     2        QMAT_SAVE, MAXSTRMSTKS, QPIVOT, INFO )

      IF ( INFO .NE. 0 ) THEN
        MAIL='Classical Beam solution LU decomposition (DGETRF)'
        CALL VLIDORT_LAPACK_ERROR ( INFO, N, MAIL, STATUS )
        IF ( STATUS .EQ. VLIDORT_SERIOUS ) RETURN
      ENDIF

C  Solution of reduced problem by back-substitution

      CALL DGETRS
     1    ( 'N', NSTKS_NSTRMS, 1,
     2       QMAT_SAVE, MAXSTRMSTKS, QPIVOT,
     3       QVEC_SAVE, MAXSTRMSTKS, INFO )

      IF ( INFO .NE. 0 ) THEN
        MAIL='Classical Beam solution back-substitution (DGETRS)'
        CALL VLIDORT_LAPACK_ERROR ( INFO, N, MAIL, STATUS )
        IF ( STATUS .EQ. VLIDORT_SERIOUS ) RETURN
      ENDIF

C  Assigning beam particular integral solution vector

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
c          if(i.eq.4.and.n.eq.24)write(38,'(4i3,1p2e20.10)')
c     &           m,ibeam,i,O1,BVEC(I,O1,N)
        ENDDO
      ENDDO

C  Old code: solution assignation

c      DO I = 1, NSTREAMS
c        I1 = I + NSTREAMS
c        DO O1 = 1, NSTOKES
c          BVEC(I,O1,N) = WVEC(I,O1,N)
c          HELP = ZERO
c          DO O2 = 1, NSTOKES
c            HELP = HELP + DMAT(O1,O2) * WVEC(I1,O2,N)
c          ENDDO
c          BVEC(I1,O1,N) = HELP
c        ENDDO
c      ENDDO

C  Values at the layer boundaries
C  (transmittance factors have been determined in SETUPS module)

      TRANS1 = INITIAL_TRANS(N,IBEAM)
      TRANS2 = T_DELT_MUBAR(N,IBEAM) * TRANS1

      DO I = 1, NSTREAMS_2
        DO O1 = 1, NSTOKES
          WUPPER(I,O1,N) = BVEC(I,O1,N)*TRANS1
          WLOWER(I,O1,N) = BVEC(I,O1,N)*TRANS2
        ENDDO
      ENDDO

C  debug

c      IF ( DO_DEBUG_WRITE ) THEN
c       J = 97
c       IF ( DO_FDTEST ) J = 98
c       IF ( N.EQ.14.and.m.eq.1) THEN
c        DO I = 1, NSTREAMS*2
c         write(*,'(4i4,1p4e17.9)')IBEAM,M,N,I,BVEC(I,1,N),
c     &      INITIAL_TRANS(N,IBEAM), T_DELT_MUBAR(N,IBEAM)
c        ENDDO
c       ENDIF
c       IF ( DO_FDTEST. and. M.EQ.1) PAUSE
c      ENDIF

C  debug check solution Forward problem - Works
C    removed 2 July. See Older version, labeled _debug

c      IF ( N.le.6.and.m.eq.0) THEN
c       DO I = 1,40
c        write(91,'(3i3,1p4e17.7)')M,N,I,(BVEC(I,o1,N),o1=1,4)
c       ENDDO
c      ENDIF

c      if (n.eq.6) pause 'end beam'

C  Finish

      RETURN
      END

C

      SUBROUTINE VLIDORT_UBEAM_SOLUTION
     I    ( GIVEN_LAYER, FOURIER, IBEAM )

C  Include files
C  =============

C  input
C  -----

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include file of setup stuff (input to this module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  include file of solution variables (output to this module)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index and Fourier number (inputs)

      INTEGER          GIVEN_LAYER
      INTEGER          FOURIER
      INTEGER          IBEAM

C  Local variables
C  ---------------

      LOGICAL          DO_LAYER
      INTEGER          UM, L, N, M, O1, O2, O3, J, J1
      DOUBLE PRECISION SUM1, SUM2, T1, T2, A5
      DOUBLE PRECISION SP, SM, SPS, SPT, SGW(MAXSTOKES)
      DOUBLE PRECISION HELP_Q1 ( 0:MAXMOMENTS, MAXSTOKES )
      DOUBLE PRECISION HELP_Q2 ( 0:MAXMOMENTS, MAXSTOKES )
      DATA             SGW / 1.0d0, 1.0d0, -1.0d0, -1.0d0 /

C  Layer and Fourier

      N = GIVEN_LAYER
      M = FOURIER

C  No particular solution beyond the cutoff layer.
C  OR No scattering in this layer
C    --> Zero the boundary layer values and exit

      IF ( N .GT. LAYER_PIS_CUTOFF(IBEAM) .OR.
     &       .NOT. DO_LAYER_SCATTERING(M,N) ) THEN
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

C  existence flag

      DO_LAYER = ( STERM_LAYERMASK_UP(N) .OR.
     &             STERM_LAYERMASK_DN(N) )

C  first function

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
          HELP_Q1(L,O1) = SPS
        ENDDO
       ENDDO
      ENDIF

C  For each moment, do inner sum over computational angles
C    Equivalent to Gvec = GVECTOR. Eq 34. GAUx = Greekmat.Gvec

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

C  Now sum over all harmonic contributions (Upwelling)
C    Direct and integarated contributions

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

C  Now sum over all harmonic contributions (Downwelling)
C    Direct and integarated contributions

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

C  debug stokes = 1

C      DO UM = LOCAL_UM_START, N_USER_STREAMS
c        write(*,*)n,um,UPAR_UP_1(UM,1,N),UPAR_UP_2(UM,1,N)
c      ENDDO
c      DO UM = LOCAL_UM_START, N_USER_STREAMS
c        write(*,*)n,um,UPAR_DN_1(UM,1,N),UPAR_DN_2(UM,1,N)
c      ENDDO
c      PAUSE

C  Finish

      RETURN
      END

C

      SUBROUTINE VLIDORT_PIMATRIX_SETUP (FOURIER)

C  Notes
C  -----

C  This is equivalent to the Legendre setup modules in LIDORT

C---------old comments ---------------------------------
C  This needs modification for the case of Refractive atmosphere,
C  because the solar zenith angle is not constant.
C-------------------------------------------------------

C  PI-matrix setup, following recipe of Siewert (1982).
C  Single Fourier component only.

C  Tested against benchmark results in Vestrucci & Siewert (1984), Problem I.

C  original coding, September 2002, R. Spurr SAO
C  Multibeam SZA coding, July 2004, R. Spurr SAO
C  Coding for refractive geometry case, R. Spurr, RT Solutions, May 2005

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include file of setup quantities (output)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  Subroutine arguments

      INTEGER          FOURIER

C  Local matrices

      DOUBLE PRECISION XLM_DIAG(0:MAXMOMENTS,MAXSTOKES)
      DOUBLE PRECISION YLM_DIAG(0:MAXMOMENTS,MAXSTOKES)
      DOUBLE PRECISION ZLM(0:MAXMOMENTS,MAXSTOKES,MAXSTOKES)
      DOUBLE PRECISION 
     &      PI(0:MAXMOMENTS,MAX_ALLSTRMS_P1,MAXSTOKES,MAXSTOKES)

C  local variables

c      INTEGER          LOCAL_NSTOKES
      DOUBLE PRECISION FAC, XSQ, DF_LPM, DF_LP1MM, ZLM_VALUE
      DOUBLE PRECISION HX, HY, HZ, H1, H2, RCONST_20, RCONST_21
      DOUBLE PRECISION PL20, RL20, PL21, RL21, TL21, XM(MAX_ALLSTRMS_P1)
      INTEGER          M, I, I1, L, SK1, SK2, SK3, SK4, N
      INTEGER          NS, NSA, M1, IB, IP

C  Set integer M = Fourier number

      M = FOURIER
      NSA = 0

C  Local NSTOKES

c      LOCAL_NSTOKES = NSTOKES
c      IF ( NSTOKES .GT. 1 ) LOCAL_NSTOKES = 4

C  total number of angles

      IF ( DO_USER_STREAMS ) THEN
        NS = NSTREAMS + N_USER_STREAMS
      ELSE
        NS = NSTREAMS
      ENDIF

C  add solar streams, depends on the use of refractive geometry

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN
c        NSA = NS + NBEAMS*NLAYERS
      ELSE
        NSA = NS + NBEAMS
      ENDIF

C  constants

      RCONST_20 = 0.25D0 * DSQRT(6.0D0)
      RCONST_21 = 2.0D0 * RCONST_20

C  Coefficient matrix PINORM (for normalization) = [ (L-m)!/(L+m)! ] ^1/2
C  ----------------------------------------------------------------------

C  .. first entry = 1 / (2m)!    ---- be careful with overflow

      PHI(M,M) = ONE
      FAC = ONE
      DO L = 1, 2*M
        FAC = FAC * DFLOAT(L)
      ENDDO
      PHI(M,M) = PHI(M,M) / FAC

C  .. Other entries by recurrence

      DO L = M + 1, NMOMENTS
        FAC = DFLOAT(L-M)/DFLOAT(L+M)
        PHI(L,M) = PHI(L-1,M)*FAC
      ENDDO

C  .. Square root

      DO L = M , NMOMENTS
        PINORM(L,M) = DSQRT(PHI(L,M))
      ENDDO

C  Additional saved quantities (to commons in VLIDORT_PISETUP.VARS)
C  ----------------------------------------------------------------

C  Only for the first fundamental harmonic

      IF ( M .EQ. 0 ) THEN

C  Sign matrix

        DO M1 = 0, NMOMENTS
          DO L = M1, NMOMENTS
            IF (MOD((L-M1),2).EQ.0) THEN
              PISIGN(L,M1) = ONE
            ELSE
              PISIGN(L,M1) = -ONE
            ENDIF
          ENDDO
        ENDDO

C  floating point integer values

        DO L = 0, NMOMENTS
          DF_L(L)    = DFLOAT(L)
          DF_LP1(L)  = DFLOAT(L+1)
          DF_2LP1(L) = DFLOAT(2*L+1)
        ENDDO

        DF_LSQM4(2) = ZERO
        DO L = 3, NMOMENTS
          DF_LSQM4(L) = DSQRT(DFLOAT(L*L-4))
        ENDDO

        DO L = 2, NMOMENTS
          DF_RT_LP3XLM1(L) = DSQRT(DFLOAT((L+3)*(L-1)))
          ZHELP(L) = TWO*DF_2LP1(L)/DF_LP1(L)/DF_L(L)
        ENDDO

      ENDIF

C  local array of all streams (every Fourier component)

        DO I = 1, NSTREAMS
          XM(I) = QUAD_STREAMS(I)
        ENDDO
        IF ( DO_USER_STREAMS ) THEN
          DO I = 1, N_USER_STREAMS
            XM(I+NSTREAMS) = USER_STREAMS(I)
          ENDDO
        ENDIF

        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
c          DO N = 1, NLAYERS
c            DO IB = 1, NBEAMS
c              IP = NS + NBEAMS * (N-1) + IB
c              XM(IP) = SUN_SZA_COSINES(N,IB)
c            ENDDO
c          ENDDO
        ELSE
          DO IB = 1, NBEAMS
            IP = IB + NS
            XM(IP) = COS_SZANGLES(IB)
          ENDDO
        ENDIF

C  factors associated with stream values
C   Special case when XM(I) = ONE, factors are zeroed, and
C    those that give singularities are avoided later (see below)
C   R. Spurr and V. Natraj, 16 january 2006

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

C  D-matrix and Mueller indices

        DO SK1 = 1, MAXSTOKES
          DO SK2 = 1, MAXSTOKES
            MUELLER_INDEX(SK1,SK2) = MAXSTOKES*(SK1-1) + SK2
            DMAT(SK1,SK2) = ZERO
          ENDDO
          IF ( SK1.GT.2) THEN
            DMAT(SK1,SK1) = -ONE
          ELSE
            DMAT(SK1,SK1) = ONE
          ENDIF
          MUELLER_DIAGONAL_INDEX(SK1) = MUELLER_INDEX(SK1,SK1)
        ENDDO

      ENDIF

C  XYZ matrices
C  ------------

C  Inverse XLM diagonal matrices (Siewert (1982), Eq. 35a)
C  YLM diagonal matrices (Siewert (1982), Eq. 35b)
C  ZLM matrices. (Siewert (1982), Eq. 35c)

C  .. for the azimuth-independent harmonic

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

C  .. for the other harmonics

      ELSE

        DO L = M, NMOMENTS
          DF_LP1MM = DFLOAT ( L + 1 - M )
          XLM_DIAG(L,1) = DF_LP1MM
          XLM_DIAG(L,1) = ONE/XLM_DIAG(L,1)          
          IF ( L. EQ. 1 ) THEN
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
          DF_LPM = DFLOAT ( L + M )
          YLM_DIAG(L,1) = DF_LPM
          YLM_DIAG(L,2) = DF_LPM * DF_LSQM4(L) / DF_L(L)
          YLM_DIAG(L,3) = YLM_DIAG(L,2)
          YLM_DIAG(L,4) = YLM_DIAG(L,1)
        ENDDO

        DO L = 2, NMOMENTS
          ZLM_VALUE = DFLOAT(M) * ZHELP(L)
          DO SK1 = 1, NSTOKES
            DO SK2 = 1, NSTOKES
              ZLM(L,SK1,SK2) = ZERO
            ENDDO
          ENDDO
          ZLM(L,2,3) = ZLM_VALUE
          ZLM(L,3,2) = ZLM_VALUE
        ENDDO

      ENDIF

C  PI calculation
C  --------------

C  Initialise all the PI matrices for given harmonic

      DO I = 1, NSA
        DO L = M, NMOMENTS
          DO SK1 = 1, NSTOKES
            DO SK2 = 1, NSTOKES
              PI(L,I,SK1,SK2) = ZERO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

C  M = 0 component. [Siewert (1982), Eqs (28) to (31)]

      IF ( M .EQ. 0 ) THEN

C  .. L = 0
        DO I = 1, NSA
          PI(0,I,1,1) = PI(0,I,1,1) + ONE
          PI(0,I,4,4) = PI(0,I,4,4) + ONE
        ENDDO
C  .. L = 1
        DO I = 1, NSA
          PI(1,I,1,1) = PI(1,I,1,1) + XM(I)
          PI(1,I,4,4) = PI(1,I,4,4) + XM(I)
        ENDDO
C  .. L = 2
        DO I = 1, NSA
          PL20 = PLEG20(I)
          RL20 = RCONST_20*UMXSQ(I)
          PI(2,I,1,1) = PI(2,I,1,1) + PL20
          PI(2,I,2,2) = PI(2,I,2,2) + RL20
          PI(2,I,3,3) = PI(2,I,2,2)
          PI(2,I,4,4) = PI(2,I,1,1)
        ENDDO
C  .. L > 2
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

C  M = 1 component. [Siewert (1982), Eqs (32) to (34), plus (35)]

      ELSE IF ( M .EQ. 1 ) THEN

C  .. L = 1
        DO I = 1, NSA
          PI(1,I,1,1) = PI(1,I,1,1) + RT_UMXSQ(I)
          PI(1,I,4,4) = PI(1,I,4,4) + RT_UMXSQ(I)
        ENDDO
C  .. L = 2
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
C  .. L > 2
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

C  M > 1 components. [Siewert (1982), Eqs (36) to (38), plus (35)]
C  Limiting case of XM(I) = ONE requires special treatment to avoid NaN
C   R. Spurr and V. Natraj, 16 january 2006

      ELSE

C  .. L = M
        IF ( M. EQ. 2 ) THEN
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
C  .. L = M + 1
        IF ( M. LT. NMOMENTS ) THEN
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
C  .. L > M + 1
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

C  Normalized output.
C  ------------------

      DO L = M, NMOMENTS

C  .. at quadrature streams

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

C  .. at positive user_defined angles

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

C  .. at solar zenith angles

C  depends on use of refractive geometry

        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
c          DO N = 1, NLAYERS
c            DO IB = 1, NBEAMS
c              IP = NS + NBEAMS * (N-1) + IB
c              DO SK1 = 1, NSTOKES
c                DO SK2 = 1, NSTOKES
c                  PI_X0P(L,IB,N,SK1,SK2) = 
c     &                 PI(L,IP,SK1,SK2) * PINORM(L,M)
c                ENDDO
c              ENDDO
c            ENDDO
c          ENDDO
        ELSE
          DO IB = 1, NBEAMS
            IP = IB + NS
            DO SK1 = 1, NSTOKES
              DO SK2 = 1, NSTOKES
                PI_X0P(L,IB,1,SK1,SK2) = 
     &               PI(L,IP,SK1,SK2) * PINORM(L,M)
                DO N = 2, NLAYERS
                  PI_X0P(L,IB,N,SK1,SK2)= PI_X0P(L,IB,1,SK1,SK2)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

C  end loop moments

      ENDDO

C  debug

C      write(77,*)'Fourier',M
C      DO L = M, NMOMENTS
C       write(77,*)'Moment ',L
C       DO SK1 = 1, 4
C        WRITE(77,'(I3,1p4e15.6)')SK1,(PI_XUP(L,4,SK1,SK2),SK2=1,4)
C       ENDDO
C      ENDDO
      
C  Finish

      RETURN
      END
    
