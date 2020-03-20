C ###############################################################
C #                                                             #
C #                    THE VLIDORT  MODEL                       #
C #                                                             #
C #  Vectorized LInearized Discrete Ordinate Radiative Transfer #
C #  -          --         -        -        -         -        #
C #                                                             #
C ###############################################################

C ###############################################################
C #                                                             #
C #  Author :      Robert. J. D. Spurr                          #
C #                                                             #
C #  Address :      RT Solutions, Inc.                          #
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
C #              VLIDORT_L_QHOM_SOLUTION                        #
C #              VLIDORT_L_UHOM_SOLUTION                        #
C #              VLIDORT_L_QBEAM_SOLUTION                       #
C #              VLIDORT_L_UBEAM_SOLUTION                       #
C #                                                             #
C ###############################################################

      SUBROUTINE VLIDORT_L_QHOM_SOLUTION
     I    ( GIVEN_LAYER, FOURIER, DO_VARY, N_PARAMETERS,
     O      STATUS )

C  Linearization of the homogeneous solutions.

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup and solution stuff (inputs to module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'

C  include file of main solution linearization variables (output)

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'

C  include file of setup linearization variables (output)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index and Fourier number (inputs)

      INTEGER         GIVEN_LAYER
      INTEGER         FOURIER

C  Variation flag for this layer

      LOGICAL         DO_VARY

C  number of varying parameters (input)

      INTEGER         N_PARAMETERS

C  output status

      INTEGER         STATUS

C  Local variables
C  ---------------

C  local matrices for eigenvalue computation
C  -----------------------------------------

C  For the real calculation based on DGEEV

      DOUBLE PRECISION
     &      HMAT1(MAXSTRMSTKS,MAXSTRMSTKS),
     &      HVEC1(MAXSTRMSTKS,MAX_ATMOSWFS),
     &      HVEC_AUX(MAXSTRMSTKS,MAX_ATMOSWFS)

C  For the complex calculations based on DGEEV

      DOUBLE PRECISION
     &      HMAT2(MAXSTRMSTKS_2,MAXSTRMSTKS_2),
     &      HVEC2(MAXSTRMSTKS_2,MAX_ATMOSWFS)
      DOUBLE PRECISION
     &      HMATRED(MAXSTRMSTKS_21,MAXSTRMSTKS_21),
     &      HVECRED(MAXSTRMSTKS_21,MAX_ATMOSWFS)
      DOUBLE PRECISION
     &      HVEC_CR(MAXSTRMSTKS,MAX_ATMOSWFS),
     &      HVEC_CI(MAXSTRMSTKS,MAX_ATMOSWFS)
      DOUBLE PRECISION
     &      ADJN_CR(MAXSTRMSTKS),
     &      ADJN_CI(MAXSTRMSTKS)

C  For the real calculation based on ASYMTX

      DOUBLE PRECISION
     &      HMAT1A(MAXSTRMSTKS_P1,MAXSTRMSTKS_P1),
     &      HVEC1A(MAXSTRMSTKS_P1,MAX_ATMOSWFS)

C  Local linearized vectors

      DOUBLE PRECISION 
     &      L_FWD_DIFVEC_R(MAXSTREAMS,MAXSTOKES),
     &      L_FWD_SUMVEC_R(MAXSTREAMS,MAXSTOKES,MAX_ATMOSWFS)
      DOUBLE PRECISION 
     &      L_FWD_DIFVEC_CR(MAXSTREAMS,MAXSTOKES),
     &      L_FWD_SUMVEC_CR(MAXSTREAMS,MAXSTOKES)
      DOUBLE PRECISION 
     &      L_FWD_DIFVEC_CI(MAXSTREAMS,MAXSTOKES),
     &      L_FWD_SUMVEC_CI(MAXSTREAMS,MAXSTOKES)

C  Pivoting arrays for LAPACK linear-algebra solvers

      INTEGER          LAPACK_IPIV1(MAXSTRMSTKS)
      INTEGER          LAPACK_IPIV21(MAXSTRMSTKS_21)
      INTEGER          LAPACK_INFO, NSTKS_NSTRMS_21, NSTMSTKS_P1
      INTEGER          LAPACK_IPIV1A(MAXSTRMSTKS_P1)

C  Miscellaneous local variables

      INTEGER          I, J, I1, J1, JC, IR, JCOL, IROW, IROW1, ILAST
      INTEGER          L, N, M, Q, O1, O2, O3, O4, NA, UTA, UT
      INTEGER          K, KO1, K0, K1, K2
      DOUBLE PRECISION DP, DM, KSQD, KVAL, KVAL2, XINV, FAC
      DOUBLE PRECISION H1P, H1M, H2P, H2M, H1, H2, L_OG
      DOUBLE PRECISION TBAS, TBAS_CR, TBAS_CI, LTQ, LCR, LCI
      DOUBLE PRECISION TAU_UP, TAU_DN, H_DN, H_UP, VAR
      DOUBLE PRECISION TBASUP_CR, TBASUP_CI, TBASDN_CR, TBASDN_CI

C  additional variables

      INTEGER          AA, AA1, TC, TR, I_NONTRIVIAL
      DOUBLE PRECISION SR, SCR, SCI, HR, HCR, HCI, MODULUS
      DOUBLE PRECISION KVAL_CR, KVAL_CI, KVAL2_CR, KVAL2_CI
      DOUBLE PRECISION MOD_KVAL, SEPCON_CR, SEPCON_CI
      DOUBLE PRECISION DENOM_CR, DENOM_CI, LFPD, LFND
      DOUBLE PRECISION LFPD1, LFPD2, LFND1, LFND2

c      added by xliu for debugging
c      INTEGER count
c      data count /0/
c      SAVE count

C  Start of code
C  -------------
C      count = count + 1

C  initialise status

      STATUS = VLIDORT_SUCCESS

C  Layer and Fourier number

      N = GIVEN_LAYER
      M = FOURIER

C  nothing to do if not varying

      IF ( .NOT. DO_VARY ) RETURN

C  For solution saving option with no scattering,
C  linearized solution vectors are all zero
C  Exit the routine by going to continuation point.

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

C  Linearize the Eigenmatrix
C  =========================

C  set up linearizations of SAB and DAB. All real variables.
C   Thoroughly debugged, November 2005

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

C  set up linearized eigenmatrices, saved for all layers. Real variables.
C   Thoroughly debugged, November 2005

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
             H2 = H2 + L_DAB(I,K,O1,O3,N,Q) *   SAB(K,J,O3,O2,N)
     &               +   DAB(I,K,O1,O3,N)   * L_SAB(K,J,O3,O2,N,Q)
            ENDDO
            H1 = H1 + H2
           ENDDO
           L_EIGENMAT(IROW,JCOL,N,Q) = H1
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

C  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C     Linearization of Real Eigensolutions
C  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

C  start eigenvalue loop

      DO K = 1, K_REAL(N)

C  some basic quantities

        AA    = EIGENMASK_R(K)
        KVAL  = KEIGEN(K,N)
        KVAL2 = 2.0D0 * KVAL
        KSQD  = KVAL * KVAL

C   REAL SOLUTIONS ONLY (ASYMTX linearization)
C   ==========================================

       IF ( DO_REAL_EIGENSOLVER(M,N) ) THEN

C  Count the number of entries in the eigenvector ( = TC)
C    If the eigenvector is trivial (one entry, TC = 1).....
C    ......Zero the results and skip linearization for this eigenvalue

        TC = 0
        DO JCOL = 1,  NSTKS_NSTRMS
          if ( RITE_EVEC(JCOL,AA).NE.ZERO ) TC = TC + 1
        ENDDO
        IF ( TC .EQ. 1 ) THEN
          DO Q = 1, N_PARAMETERS
            HVEC1A(1,Q) = ZERO
            DO I = 1, NSTKS_NSTRMS
              HVEC1A(I+1,Q) = ZERO
            ENDDO
          ENDDO
          GO TO 4555
        ENDIF

C  Nothing to do also if there is degeneracy
C    NEW CODE by R. Spurr 22 August 2008.

        IF ( EIGENDEGEN(K,N) ) THEN
          DO Q = 1, N_PARAMETERS
            HVEC1A(1,Q) = ZERO
            DO I = 1, NSTKS_NSTRMS
              HVEC1A(I+1,Q) = ZERO
            ENDDO
          ENDDO
          GO TO 4555
        ENDIF

C  initialise solution matrix HMAT1A (important to do this!)
C    HMAT1A is destroyed upon passage through the LAPACK linear modules

        DO I = 1, MAXSTRMSTKS_P1
c          LAPACK_IPIV1A(I) = 0
          DO J = 1, MAXSTRMSTKS_P1
            HMAT1A(I,J) = ZERO
          ENDDO
        ENDDO

C  Determine solution matrix HMAT (A in the matrix equation A.x = B)

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

c        HMAT1A(NSTMSTKS_P1,1) = 1.879d-15   ! shoudl technically be zero

C  solution column vectors (B in the matrix equation A.x = B).
C    (One for each parameter to be varied)

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

C  Find any degenerate columns J, and replace the (J,J) entry with
C   a small number. This is a quick fix that ensures that the LAPACK
C   routine does not fall over because of degeneracy.
C    --This problem is only present in the vector case.

C  Fiddle code

c        DO J = 1, NSTKS_NSTRMS
c          NC = 0
c          DO I = 1,  NSTKS_NSTRMS
c           IF ( HMAT1A(I,J+1).EQ.ZERO) NC = NC + 1
c          ENDDO
c          IF ( NC .EQ. NSTKS_NSTRMS ) THEN
c            do I = 1, NSTKS_NSTRMS
c               HMAT1A(J,I+1) = 1.237d-15*DBLE(i)
c            enddo
c            if ( HMAT1A(J+1,J+1).EQ.ZERO)
c     &    HMAT1A(J+1,J+1) = 1.056d-15
c          ENDIF
c        ENDDO

C  ultimate fiddle

c        DO J = 1, NSTMSTKS_P1
c           Do I = 1, NSTMSTKS_P1 
c             IF (HMAT1A(I,J).EQ.ZERO )HMAT1A(I,J) = 1.0d-15*(I+J)
c           ENDDO
c        ENDDO

C  solve System for the linearized sum-vector + eigenvalue, together
C     Solve matrix system  using LAPACK modules DGETRF and DEGTRS

!        IF (m.eq.2.and.n.eq.18) THEN
!           WRITE(45, *) m, k, N, NSTMSTKS_P1, LAPACK_INFO
!           WRITE(45, '(17D24.12)') 
!     &          ((HMAT1A(j, i), i = 1, NSTMSTKS_P1), j = 1, NSTMSTKS_P1)
!c           if(lapack_info.ne.0)STOP'2,18,9'
!        ENDIF

C#########################################################################

C  Resume normal code

        CALL DGETRF  ( NSTMSTKS_P1, NSTMSTKS_P1, HMAT1A,
     &                 MAXSTRMSTKS_P1,  LAPACK_IPIV1A, LAPACK_INFO )

        IF ( LAPACK_INFO .NE. 0 ) THEN
          CALL VLIDORT_LAPACK_ERROR( LAPACK_INFO, N, 
     I       ' Homogeneous Linearization DGETRF (1)', STATUS )
          IF ( STATUS .EQ. VLIDORT_SERIOUS ) RETURN
        ENDIF

        CALL DGETRS ('N', NSTMSTKS_P1, N_PARAMETERS, HMAT1A,
     &               MAXSTRMSTKS_P1, LAPACK_IPIV1A,  HVEC1A,
     &               MAXSTRMSTKS_P1, LAPACK_INFO )

        IF ( LAPACK_INFO .NE. 0 ) THEN
          CALL VLIDORT_LAPACK_ERROR( LAPACK_INFO, N, 
     I       ' Homogeneous Linearization DGETRS (1)', STATUS )
          IF ( STATUS .EQ. VLIDORT_SERIOUS ) RETURN
        ENDIF

C  Continuation point if linearized skipped

 4555   CONTINUE

C  Start loop over varying parameters, to assign answers

        DO Q = 1, N_PARAMETERS

          L_KEIGEN(K,N,Q) = HVEC1A(1,Q)
          DO I = 1, NSTREAMS
           IR = NSTOKES*(I-1)
           DO O1 = 1, NSTOKES
            IROW = IR + O1
            L_FWD_SUMVEC_R(I,O1,Q) = HVEC1A(IROW+1,Q)
           ENDDO
          ENDDO

C  Debug (important)
C   Tested 28 December 2005.

c         IF ( DO_DEBUG_WRITE ) THEN
c          IF (Q.EQ.1.and.m.eq.2.and.n.eq.1) THEN
c           write(37,'(2i3,L3,1p2e20.12)')
c     &              M,K,EIGENDEGEN(K,N),KEIGEN(K,N),L_KEIGEN(K,N,Q)
c           DO I = 1, NSTREAMS
c            DO O1 = 1, NSTOKES
c                 write(37,'(4i3,1p2e20.12)')
c     &        m,k,i,o1,FWD_SUMVEC(I,O1,K),L_FWD_SUMVEC_R(I,O1,Q)
c            ENDDO
c           ENDDO
c          ENDIF
c         ENDIF

        ENDDO

C   REAL SOLUTIONS (DGEEV linearization)
C   ====================================

       ELSE

C  initialise

       I_NONTRIVIAL = 0

C  Count the number of entries in the eigenvector ( = TC)
C    If the eigenvector is trivial (one entry, TC = 1).....
C    ......Zero the results and skip linearization for this eigenvalue

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

C  Establish the first non-trivial row
C      This is necessary to avoid underflow problems

        TR = 0
        DO JCOL = 1,  NSTKS_NSTRMS
          if ( DABS(RITE_EVEC(JCOL,AA)).GT.1.0D-10 ) THEN
            TR = JCOL
            GO TO 3442
          ENDIF
        ENDDO
 3442   CONTINUE

C  Use the adjoint method to get L(k) = <X^,L(G)X> / 2k <X^,X>

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

C  Now solve for the linearized eigenvector L(X) from A.L(X) = B
C     Determine matrix HMAT (=A) and vectors HVEC (B)

C  Initialise first, as HMAT is destroyed after LAPACK routines

        DO IROW = 1, NSTKS_NSTRMS
          DO JCOL = 1, NSTKS_NSTRMS
            HMAT1(IROW,JCOL) = ZERO
          ENDDO
        ENDDO

C  set up the matrix HMAT
C    Use the linearization of the Eigen-equation for all but one
C    of the rows. For the FIRST row, use the linearization of the
C    normalization condition on the eigenvectors.

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

C  set up the HVEC:
C    Use the linearization of the Eigen-equation for all but one
C    of the rows. For the FIRST row, use the linearization of the
C    normalization condition on the eigenvectors.

        DO Q = 1, N_PARAMETERS
          I_NONTRIVIAL = 0
          DO IROW = 1, NSTKS_NSTRMS
           IF ( IROW.NE.TR ) THEN
            HVEC1(IROW,Q) = HVEC_AUX(IROW,Q) -
     &         KVAL2 * L_KEIGEN(K,N,Q) * RITE_EVEC(IROW,AA)
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

C  Solve matrix system  using LAPACK modules DGETRF and DEGTRS
C  Only do if vector has non-trivial entries 

        IF ( I_NONTRIVIAL.GT.0) THEN

          CALL DGETRF  ( NSTKS_NSTRMS, NSTKS_NSTRMS, HMAT1,
     &                 MAXSTRMSTKS,  LAPACK_IPIV1, LAPACK_INFO )

          IF ( LAPACK_INFO .NE. 0 ) THEN
           CALL VLIDORT_LAPACK_ERROR( LAPACK_INFO, N, 
     I       ' Homogeneous Linearization DGETRF (2)', STATUS )
           IF ( STATUS .EQ. VLIDORT_SERIOUS ) RETURN
          ENDIF

          CALL DGETRS ('N', NSTKS_NSTRMS, N_PARAMETERS, HMAT1,
     &               MAXSTRMSTKS, LAPACK_IPIV1,       HVEC1,
     &               MAXSTRMSTKS, LAPACK_INFO )

          IF ( LAPACK_INFO .NE. 0 ) THEN
           CALL VLIDORT_LAPACK_ERROR( LAPACK_INFO, N, 
     I       ' Homogeneous Linearization DGETRS (2)', STATUS )
           IF ( STATUS .EQ. VLIDORT_SERIOUS ) RETURN
          ENDIF

        ENDIF

C  Continuation point if linearization is skipped

 4556   CONTINUE

C  Start loop over varying parameters, and assign
C    linearized sum vector = linearized eigensolution

        DO Q = 1, N_PARAMETERS

         DO I = 1, NSTREAMS
           IR = NSTOKES*(I-1)
           DO O1 = 1, NSTOKES
            IROW = IR + O1
            L_FWD_SUMVEC_R(I,O1,Q) = HVEC1(IROW,Q)
           ENDDO
         ENDDO

C  Debug (important)
C   Tested 28 December 2005.

c         IF ( DO_DEBUG_WRITE ) THEN
c          IF (Q.EQ.1.and.m.eq.0.and.n.eq.1) THEN
c           write(37,'(2i3,1p2e20.12)')
c     &              M,K,KEIGEN(K,N),L_KEIGEN(K,N,Q)
c           DO I = 1, NSTREAMS
c            DO O1 = 1, NSTOKES
c                 write(37,'(4i3,1p2e20.12)')
c     &        m,k,i,o1,FWD_SUMVEC(I,O1,K),L_FWD_SUMVEC_R(I,O1,Q)
c            ENDDO
c           ENDDO
c          ENDIF
c         ENDIF

        ENDDO
 
C  End Clause ASYMTX vs. DGEEV

       ENDIF

C  Resume General Code
C  ===================

C  Start loop over varying parameters

       DO Q = 1, N_PARAMETERS

C  linearized difference vector

         DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
           IROW = IR + O1
           SR = ZERO
           DO J = 1, NSTREAMS
            JC = NSTOKES*(J-1)
            DO O2 = 1, NSTOKES
             JCOL = JC + O2
             SR = SR - L_SAB(I,J,O1,O2,N,Q) *   FWD_SUMVEC(J,O2,K)
     &               -   SAB(I,J,O1,O2,N)   * L_FWD_SUMVEC_R(J,O2,Q)
            ENDDO
           ENDDO
           HR = ( - SR - L_KEIGEN(K,N,Q) * FWD_DIFVEC(I,O1,K) )
           L_FWD_DIFVEC_R(I,O1) = HR / KVAL
          ENDDO
         ENDDO

C  assign Forward solutions PHI (Siewert's notation)
C    --->  first N are "DOWN", last N are "UP" (streams)
C    --->  Use symmetry properties to set -ve eigensolutions

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

C  Older code; Uses additional wasteful holding arrays
C  solution assignation (only for BVP)

c         DO I = 1, NSTREAMS
c          I1 = I + NSTREAMS
c          DO O1 = 1, NSTOKES
c           L_FWD_XPOS_R(I,O1,K,N,Q)  = HALF * 
c     &          ( L_FWD_SUMVEC_R(I,O1,Q) + L_FWD_DIFVEC_R(I,O1) )
c           L_FWD_XPOS_R(I1,O1,K,N,Q) = HALF * 
c     &          ( L_FWD_SUMVEC_R(I,O1,Q) - L_FWD_DIFVEC_R(I,O1) )
c           L_FWD_XNEG_R(I1,O1,K,N,Q) =  L_FWD_XPOS_R(I,O1,K,N,Q)
c           L_FWD_XNEG_R(I,O1,K,N,Q)  =  L_FWD_XPOS_R(I1,O1,K,N,Q)
c          ENDDO
c         ENDDO
c         DO I = 1, NSTREAMS
c          I1 = I + NSTREAMS
c          DO O1 = 1, NSTOKES
c           L_SOLA_XPOS_R(I,O1,K,N,Q)  = L_FWD_XPOS_R(I,O1,K,N,Q)
c           L_SOLB_XNEG_R(I,O1,K,N,Q)  = L_FWD_XNEG_R(I,O1,K,N,Q)
c           H1 = ZERO
c           H2 = ZERO
c           DO O2 = 1, NSTOKES
c            H1 = H1 + DMAT(O1,O2)*L_FWD_XNEG_R(I,O2,K,N,Q)
c            H2 = H2 + DMAT(O1,O2)*L_FWD_XPOS_R(I,O2,K,N,Q)
c           ENDDO
c           L_SOLA_XPOS_R(I1,O1,K,N,Q)  = H1
c           L_SOLB_XNEG_R(I1,O1,K,N,Q)  = H2
c          ENDDO
c         ENDDO

C  End parameter loop

       ENDDO

C  End real eigenvalue loop

      ENDDO

C  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C     Linearization of COMPLEX Eigensolutions
C  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

C  Only if flagged

      IF ( .NOT. DO_REAL_EIGENSOLVER(M,N) ) THEN

C  Offset

       KO1 = K_REAL(N) + 1

C  start loop over eigensolutions

       DO K = 1, K_COMPLEX(N)

C  Bookkeeping indices

        K0 = 2 * ( K - 1 )
        K1 = KO1 + K0
        K2 = K1  + 1

C  Basic quantities

        AA  = EIGENMASK_C(K)
        AA1 = AA + 1
        KVAL_CR  = KEIGEN(K1,N)
        KVAL_CI  = KEIGEN(K2,N)
        KVAL2_CR = TWO * KVAL_CR
        KVAL2_CI = TWO * KVAL_CI
        MOD_KVAL = KEIGEN_CSQ(K)
        SEPCON_CR =   KEIGEN(K1,N) / MOD_KVAL
        SEPCON_CI = - KEIGEN(K2,N) / MOD_KVAL

C  Establish the zero imaginary part of one eigenvector component

        TR = 0
        DO JCOL = 1,  NSTKS_NSTRMS
          if ( RITE_EVEC(JCOL,AA1).EQ.ZERO ) THEN
            TR = JCOL + NSTKS_NSTRMS
            GO TO 3443
          ENDIF
        ENDDO
 3443   CONTINUE

C  Use the adjoint theory to get L(k) from 2kL(k) = <X^,L(G)X> / <X^,X>
C      MUST USE THE COMPLEX CONJUGATE OF THE LEFT EIGENVECTOR.
  
C  --------------- Adjoint  norm

        SCR = ZERO
        SCI = ZERO
        DO JCOL = 1,  NSTKS_NSTRMS
          HCR = + LEFT_EVEC(JCOL,AA)  * RITE_EVEC(JCOL,AA)
     &          + LEFT_EVEC(JCOL,AA1) * RITE_EVEC(JCOL,AA1)
          SCR = SCR + HCR
          HCI = - LEFT_EVEC(JCOL,AA1) * RITE_EVEC(JCOL,AA)
     &          + LEFT_EVEC(JCOL,AA)  * RITE_EVEC(JCOL,AA1)
          SCI = SCI + HCI
          ADJN_CR(JCOL) = + KVAL2_CR * RITE_EVEC(JCOL,AA)
     &                    - KVAL2_CI * RITE_EVEC(JCOL,AA1)
          ADJN_CI(JCOL) = + KVAL2_CI * RITE_EVEC(JCOL,AA)
     &                    + KVAL2_CR * RITE_EVEC(JCOL,AA1)
        ENDDO
        DENOM_CR = KVAL2_CR * SCR - KVAL2_CI * SCI
        DENOM_CI = KVAL2_CR * SCI + KVAL2_CI * SCR
        MODULUS  = ONE / ( DENOM_CR*DENOM_CR + DENOM_CI * DENOM_CI )

C  ---------------- Adjoint numerator, and final computation

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
          HCR = HCR + LEFT_EVEC(JCOL,AA)  * HVEC_CR(JCOL,Q)
     &              + LEFT_EVEC(JCOL,AA1) * HVEC_CI(JCOL,Q)
          HCI = HCI - LEFT_EVEC(JCOL,AA1) * HVEC_CR(JCOL,Q)
     &              + LEFT_EVEC(JCOL,AA)  * HVEC_CI(JCOL,Q)
         ENDDO
         L_KEIGEN(K1,N,Q) = ( HCR*DENOM_CR + HCI*DENOM_CI )*MODULUS
         L_KEIGEN(K2,N,Q) = ( HCI*DENOM_CR - HCR*DENOM_CI )*MODULUS
        ENDDO

C  Now solve for the linearized eigenvector L(X) from A.L(X) = B
C     Determine matrix HMAT (=A) and vectors HVEC (B)

C  ---------- Initialise first, as HMAT is destroyed after LAPACK routines

        DO IROW = 1, NSTKS_NSTRMS_2
          DO JCOL = 1, NSTKS_NSTRMS_2
            HMAT2(IROW,JCOL) = ZERO
          ENDDO
        ENDDO

C  ---------- set up matrix for Full solution

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

C  ---------- set up the vectors for full solution

        DO Q = 1, N_PARAMETERS
          DO I = 1, NSTKS_NSTRMS
           I1 = I + NSTKS_NSTRMS
           HCR = + L_KEIGEN(K1,N,Q) * ADJN_CR(I) 
     &           - L_KEIGEN(K2,N,Q) * ADJN_CI(I)
           HCI = + L_KEIGEN(K1,N,Q) * ADJN_CI(I) 
     &           + L_KEIGEN(K2,N,Q) * ADJN_CR(I)
           HVEC2(I,Q)  = HVEC_CR(I,Q) - HCR
           HVEC2(I1,Q) = HVEC_CI(I,Q) - HCI
          ENDDO
        ENDDO

C  ---------- Strike out TR row/column to get reduced matrix and vectors

        NSTKS_NSTRMS_2  = 2 * NSTKS_NSTRMS
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

C  ----------  Replace the last row by linearized normalization constraint
c                ( IF (TR.LT.NSTKS_NSTRMS_2) THEN  possible exception......
C                Should not be necessary to flag this case )
C   
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

C  ---------- Solve matrix system using LAPACK modules DGETRF and DEGTRS

        CALL DGETRF  ( NSTKS_NSTRMS_21, NSTKS_NSTRMS_21, HMATRED,
     &                 MAXSTRMSTKS_21,  LAPACK_IPIV21,   LAPACK_INFO )
        IF ( LAPACK_INFO .NE. 0 ) THEN
          CALL VLIDORT_LAPACK_ERROR( LAPACK_INFO, N, 
     I       ' Homogeneous Linearization DGETRF (2)', STATUS )
          IF ( STATUS .EQ. VLIDORT_SERIOUS ) RETURN
        ENDIF

        CALL DGETRS ('N', NSTKS_NSTRMS_21, N_PARAMETERS, HMATRED,
     &               MAXSTRMSTKS_21, LAPACK_IPIV21,      HVECRED,
     &               MAXSTRMSTKS_21, LAPACK_INFO )
        IF ( LAPACK_INFO .NE. 0 ) THEN
          CALL VLIDORT_LAPACK_ERROR( LAPACK_INFO, N, 
     I       ' Homogeneous Linearization DGETRS (2)', STATUS )
          IF ( STATUS .EQ. VLIDORT_SERIOUS ) RETURN
        ENDIF

C  Start loop over varying parameters

        DO Q = 1, N_PARAMETERS

C  linearized sum vector = linearized eigensolution

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

C  Debug code (important)

c         IF (do_debug_write) THEN
c          if ( Q.EQ.1.and.m.eq.0.and.n.eq.1) then
c           write(91,'(2i3,1p4e20.10)')
c     &            m,k,KEIGEN(K1,N),KEIGEN(K2,N),
c     &            L_KEIGEN(K1,N,Q),L_KEIGEN(K2,N,Q)
c           DO I = 1, NSTREAMS
c            DO O1 = 1, NSTOKES
c              write(91,'(4i3,1p4e20.10)')
c     &        m,k,i,o1,FWD_SUMVEC(I,O1,K1),FWD_SUMVEC(I,O1,K2),
c     &         L_FWD_SUMVEC_CR(I,O1),L_FWD_SUMVEC_CI(I,O1)
c            ENDDO
c           ENDDO
c          ENDIF
c         ENDIF

C  Determine the linearized difference vector

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
             SCR = SCR - L_SAB(I,J,O1,O2,N,Q) *   FWD_SUMVEC(J,O2,K1)
     &                 -   SAB(I,J,O1,O2,N)   * L_FWD_SUMVEC_CR(J,O2)
             SCI = SCI - L_SAB(I,J,O1,O2,N,Q) *   FWD_SUMVEC(J,O2,K2)
     &                 -   SAB(I,J,O1,O2,N)   * L_FWD_SUMVEC_CI(J,O2)
            ENDDO
           ENDDO
           HCR = - SCR - L_KEIGEN(K1,N,Q) * FWD_DIFVEC(I,O1,K1)
     &                 + L_KEIGEN(K2,N,Q) * FWD_DIFVEC(I,O1,K2)
           HCI = - SCI - L_KEIGEN(K2,N,Q) * FWD_DIFVEC(I,O1,K1)
     &                 - L_KEIGEN(K1,N,Q) * FWD_DIFVEC(I,O1,K2)
           L_FWD_DIFVEC_CR(I,O1) = HCR*SEPCON_CR - HCI*SEPCON_CI
           L_FWD_DIFVEC_CI(I,O1) = HCI*SEPCON_CR + HCR*SEPCON_CI
          ENDDO
         ENDDO

C  assign Forward solutions PHI (Siewert's notation)
C    --->  first N are "DOWN", last N are "UP" (streams)
C    --->  Use symmetry properties to set -ve eigensolutions

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

C  Older code: Uses additional wasteful holding arrays
C  solution assignation (only for BVP)

c         DO I = 1, NSTREAMS
c          I1 = I + NSTREAMS
c          DO O1 = 1, NSTOKES
c           L_FWD_XPOS_CR(I1,O1,K,N,Q) = HALF * 
c     &          ( L_FWD_SUMVEC_CR(I,O1) - L_FWD_DIFVEC_CR(I,O1) )
c           L_FWD_XPOS_CR(I,O1,K,N,Q)  = HALF * 
c     &          ( L_FWD_SUMVEC_CR(I,O1) + L_FWD_DIFVEC_CR(I,O1) )
c           L_FWD_XPOS_CI(I1,O1,K,N,Q) = HALF * 
c     &          ( L_FWD_SUMVEC_CI(I,O1) - L_FWD_DIFVEC_CI(I,O1) )
c           L_FWD_XPOS_CI(I,O1,K,N,Q)  = HALF * 
c     &          ( L_FWD_SUMVEC_CI(I,O1) + L_FWD_DIFVEC_CI(I,O1) )
c           L_FWD_XNEG_CR(I1,O1,K,N,Q) =  L_FWD_XPOS_CR(I,O1,K,N,Q)
c           L_FWD_XNEG_CR(I,O1,K,N,Q)  =  L_FWD_XPOS_CR(I1,O1,K,N,Q)
c           L_FWD_XNEG_CI(I1,O1,K,N,Q) =  L_FWD_XPOS_CI(I,O1,K,N,Q)
c           L_FWD_XNEG_CI(I,O1,K,N,Q)  =  L_FWD_XPOS_CI(I1,O1,K,N,Q)
c          ENDDO
c         ENDDO
c         DO I = 1, NSTREAMS
c          I1 = I + NSTREAMS
c          DO O1 = 1, NSTOKES
c           L_SOLA_XPOS_CR(I,O1,K,N,Q)  = L_FWD_XPOS_CR(I,O1,K,N,Q)
c           L_SOLB_XNEG_CR(I,O1,K,N,Q)  = L_FWD_XNEG_CR(I,O1,K,N,Q)
c           L_SOLA_XPOS_CI(I,O1,K,N,Q)  = L_FWD_XPOS_CI(I,O1,K,N,Q)
c           L_SOLB_XNEG_CI(I,O1,K,N,Q)  = L_FWD_XNEG_CI(I,O1,K,N,Q)
c           H1 = ZERO
c           H2 = ZERO
c           DO O2 = 1, NSTOKES
c            H1 = H1 + DMAT(O1,O2)*L_FWD_XNEG_CR(I,O2,K,N,Q)
c            H2 = H2 + DMAT(O1,O2)*L_FWD_XPOS_CR(I,O2,K,N,Q)
c           ENDDO
c           L_SOLA_XPOS_CR(I1,O1,K,N,Q)  = H1
c           L_SOLB_XNEG_CR(I1,O1,K,N,Q)  = H2
c           H1 = ZERO
c           H2 = ZERO
c           DO O2 = 1, NSTOKES
c            H1 = H1 + DMAT(O1,O2)*L_FWD_XNEG_CI(I,O2,K,N,Q)
c            H2 = H2 + DMAT(O1,O2)*L_FWD_XPOS_CI(I,O2,K,N,Q)
c           ENDDO
c           L_SOLA_XPOS_CI(I1,O1,K,N,Q)  = H1
c           L_SOLB_XNEG_CI(I1,O1,K,N,Q)  = H2
c          ENDDO
c         ENDDO

C  End parameter loop

        ENDDO

C  End complex eigenvalue loop

       ENDDO

C  End clause for not using the real eigensolver

      ENDIF

C  control point for avoiding eigensolver
C  ======================================

 3456 CONTINUE

C  Linearized Eigenstream transmittance factors for whole layer
C  ------------------------------------------------------------

C  When the solution saving option is set, then if there is no
C  scattering in this layer, then linearized transmittances are
C  linearized discrete ordinate transmittances.

      IF ( DO_SOLUTION_SAVING .AND.
     &               .NOT.DO_LAYER_SCATTERING(M,N) ) THEN

        DO Q = 1, N_PARAMETERS
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              L_T_DELT_EIGEN(IROW,N,Q) = L_T_DELT_DISORDS(I,N,Q)
            ENDDO
          ENDDO
        ENDDO

C  Otherwise compute them as normal
 
      ELSE

C  Real

        DO K = 1, K_REAL(N)
          TBAS = - T_DELT_EIGEN(K,N) * DELTAU_VERT(N)
          DO Q = 1, N_PARAMETERS
            LTQ = L_KEIGEN(K,N,Q) + KEIGEN(K,N)*L_DELTAU_VERT(Q,N)
            L_T_DELT_EIGEN(K,N,Q) = TBAS * LTQ
          ENDDO
        ENDDO

c  Complex

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

C  debug 
c      if ( m.eq.0.and.n.eq.20 ) then
c       DO K = 1, min(6,K_REAL(N))
c         WRITE(*,'(2i4,1p6e24.12)')m,k,L_T_DELT_EIGEN(K,N,1)
c       ENDDO
c      ENDIF

C  Eigenstream transmittance factors for partial layers
C  ----------------------------------------------------

C  Done irrespective of solution saving
C    Loop over user optical depths

      DO UTA = 1, N_USER_LEVELS
        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
         UT = PARTLAYERS_OUTINDEX(UTA)
         NA  = PARTLAYERS_LAYERIDX(UT)

C  If the layer of occurrence is the given layer

         IF ( NA .EQ. N ) THEN

C  When the solution saving option is set, then if there is no
C  scattering in this layer, then transmittances are just the
C  discrete ordinate transmittances.

          IF ( DO_SOLUTION_SAVING .AND.
     &         .NOT.DO_LAYER_SCATTERING(M,N) ) THEN

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

C  Otherwise, Compute the Eigenstream transmittance factors

          ELSE

C  Partial layer optical depths

           TAU_DN = PARTAU_VERT(UT)
           TAU_UP = DELTAU_VERT(N) - TAU_DN

C  Real eigenvalues

           DO K = 1, K_REAL(N)
            H_DN = - TAU_DN * T_UTDN_EIGEN(K,UT)
            H_UP = - TAU_UP * T_UTUP_EIGEN(K,UT)
            DO Q = 1, N_PARAMETERS
              VAR = L_KEIGEN(K,N,Q) + KEIGEN(K,N) * L_DELTAU_VERT(Q,N)
              L_T_UTDN_EIGEN(K,UT,Q) = H_DN * VAR
              L_T_UTUP_EIGEN(K,UT,Q) = H_UP * VAR
            ENDDO
           ENDDO

C  Complex eigenvalues
C    Bug fixed 16 December 2005

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

C  End clause for using Discrete ordinate valuaes or not

          ENDIF

C  End loop over partial layer depths

         ENDIF
        ENDIF
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE VLIDORT_L_UHOM_SOLUTION
     I    ( GIVEN_LAYER, FOURIER, DO_VARY, N_PARAMETERS )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include file of setups/solutions/multipliers stuff (inputs to module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  include file of main solution linearization variables (output)

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'

C  include file of setup linearization variables (output)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'

C  include file of multiplier linearization variables (output)

      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index, Fourier index, parameter number (inputs)

      INTEGER          GIVEN_LAYER
      INTEGER          FOURIER
      LOGICAL          DO_VARY
      INTEGER          N_PARAMETERS

C  Local variables
C  ---------------

      DOUBLE PRECISION L_HELPSTOKES_R(MAXSTOKES)
      DOUBLE PRECISION L_HELPSTOKES_CR(MAXSTOKES)
      DOUBLE PRECISION L_HELPSTOKES_CI(MAXSTOKES)
      DOUBLE PRECISION L_GAUX_R(0:MAXMOMENTS,MAXSTOKES)
      DOUBLE PRECISION L_GAUX_CR(0:MAXMOMENTS,MAXSTOKES)
      DOUBLE PRECISION L_GAUX_CI(0:MAXMOMENTS,MAXSTOKES)

      LOGICAL          DOLAYER_PP
      INTEGER          UM, J, L, N, M, Q, O1, O2
      INTEGER          K, KO1, K0, K1, K2
      DOUBLE PRECISION SPS, H, SP, SM, SG, SN, SNS, HPC, HMC
      DOUBLE PRECISION SPSCR, SPSCI, SPCR, SPCI
      DOUBLE PRECISION SNSCR, SNSCI, SNCR, SNCI
      DOUBLE PRECISION SR, SI, SMCR, SMCI, ZEPSQ, ZEMSQ
      DOUBLE PRECISION KISQ, RHO_P_CR, RHO_M_CR, MP, MM
      DOUBLE PRECISION L_KI, L_KR, L_KISQ, LMP, LMM

C  Layer and Fourier

      N = GIVEN_LAYER
      M = FOURIER

C  general flag

      DOLAYER_PP = ( STERM_LAYERMASK_UP(N) .OR.
     &               STERM_LAYERMASK_DN(N) )

C  Linearized Zeta constants
C  =========================

      IF ( DOLAYER_PP ) THEN

       DO UM = LOCAL_UM_START, N_USER_STREAMS
        SM = USER_SECANTS(UM)

C  Linearized Real Zeta constants (always required)

        DO K = 1, K_REAL(N)
          ZEPSQ = ZETA_P(K,UM,N) * ZETA_P(K,UM,N)
          ZEMSQ = ZETA_M(K,UM,N) * ZETA_M(K,UM,N)
          DO Q = 1, N_PARAMETERS
            L_ZETA_P(K,UM,N,Q) = - L_KEIGEN(K,N,Q) * ZEPSQ
            L_ZETA_M(K,UM,N,Q) = + L_KEIGEN(K,N,Q) * ZEMSQ
          ENDDO
        ENDDO

C  linearization of Complex Zeta constants

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

C  End zeta linearization

       ENDDO
      ENDIF

C  Eigenvector interpolation to user-defined angles
C  ================================================

C  Only do this if both flags are set,
C    Or if there is a solution for scattering in this layer

      IF ( .NOT.DOLAYER_PP .OR. .NOT. DO_VARY .OR.
     &     .NOT. DO_LAYER_SCATTERING(M,N) ) RETURN

C  For each parameter

      DO Q = 1, N_PARAMETERS

C  User defined solutions (Real)
C  -----------------------------

       DO K = 1, K_REAL(N)

C  For each moment, do inner sum over computational angles
C  for the positive and negative linearized eigenvectors.
C  Linearized the product with OMEGA_GREEK ---> L_GAUX

        DO L = M, NMOMENTS
         DO O1 = 1, NSTOKES
          SPS = ZERO
          DO J = 1, NSTREAMS
           H = QUAD_HALFWTS(J)
           SP = ZERO
           SM = ZERO
           DO O2 = 1, NSTOKES
c older     SP = SP + H*PI_XQP    (L,J,O1,O2)*L_FWD_XPOS_R(J,O2,K,N,Q)
c older     SM = SM + H*PI_XQM_PRE(L,J,O1,O2)*L_FWD_XNEG_R(J,O2,K,N,Q)
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
           SG = SG + OMEGA_GREEK(L,N,O1,O2)     * L_HELPSTOKES_R(O2)
     &             + L_OMEGA_GREEK(L,N,O1,O2,Q) *   HELPSTOKES(L,K,O2)
          ENDDO
          L_GAUX_R(L,O1) = SG
         ENDDO
        ENDDO

C  Now sum over all harmonic contributions (downwelling)

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

C  Debug (important)

        if ( do_debug_write.and.q.eq.2) then
         DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            write(95,'(4i3,1p8e20.10)')m,k,um,o1,
     &         L_UHOM_UPUP(UM,O1,K,N,Q),
     &         L_UHOM_UPDN(UM,O1,K,N,Q),
     &         L_UHOM_DNUP(UM,O1,K,N,Q),
     &         L_UHOM_DNDN(UM,O1,K,N,Q),
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
c older     SPCR = SPCR + HPC * L_FWD_XPOS_CR(J,O2,K,N,Q)
c older     SPCI = SPCI + HPC * L_FWD_XPOS_CI(J,O2,K,N,Q)
            SPCR = SPCR + HPC * L_SOLA_XPOS(J,O2,K1,N,Q)
            SPCI = SPCI + HPC * L_SOLA_XPOS(J,O2,K2,N,Q)
            HMC  =  PI_XQM_PRE(L,J,O1,O2)*QUAD_HALFWTS(J)
c older     SMCR = SMCR + HMC * L_FWD_XNEG_CR(J,O2,K,N,Q)
c older     SMCI = SMCI + HMC * L_FWD_XNEG_CI(J,O2,K,N,Q)
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
           SR = SR + OMEGA_GREEK(L,N,O1,O2)     *L_HELPSTOKES_CR(O2)
     &             + L_OMEGA_GREEK(L,N,O1,O2,Q) *  HELPSTOKES(L,K1,O2)
           SI = SI + OMEGA_GREEK(L,N,O1,O2)     *L_HELPSTOKES_CI(O2)
     &             + L_OMEGA_GREEK(L,N,O1,O2,Q) *  HELPSTOKES(L,K2,O2)
          ENDDO
          L_GAUX_CR(L,O1) = SR
          L_GAUX_CI(L,O1) = SI
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

C  Debug (important)

        if ( do_debug_write.and.q.eq.2) then
         DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            write(97,'(4i3,1p8e20.10)')m,k,um,o1,
     &         L_UHOM_UPUP(UM,O1,K2,N,Q),
     &         L_UHOM_UPDN(UM,O1,K2,N,Q),
     &         L_UHOM_DNUP(UM,O1,K2,N,Q),
     &         L_UHOM_DNDN(UM,O1,K2,N,Q),
     &         UHOM_UPUP(UM,O1,K2,N),
     &         UHOM_UPDN(UM,O1,K2,N),
     &         UHOM_DNUP(UM,O1,K2,N),
     &         UHOM_DNDN(UM,O1,K2,N)
           enddo
         enddo
        endif

C  end loop over complex eigensolutions

       ENDDO

C  End loop over parameters

      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE VLIDORT_L_QBEAM_SOLUTION
     I    ( GIVEN_LAYER, FOURIER, IBEAM, DOVARY, N_PARAMETERS,
     O      STATUS )

C  linearized values of the classical particular solution.

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include file of setup and solution stuff (inputs to module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'

C  include file of linearization inputs

      INCLUDE '../includes/VLIDORT_L_INPUTS.VARS'

C  include file of main solution linearization variables (output)

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'

C  include file of setup linearization variables (output)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index and Fourier number, beam index (inputs)

      INTEGER          GIVEN_LAYER
      INTEGER          FOURIER
      INTEGER          IBEAM

C  variation flag

      LOGICAL          DOVARY

C  number of varying parameters (input)

      INTEGER          N_PARAMETERS

C  output status

      INTEGER          STATUS

C  Local variables
C  ---------------

C  Saved vector

      DOUBLE PRECISION QSAVE(MAXSTRMSTKS)

C  linearization arrays

      DOUBLE PRECISION L_QVEC(MAXSTRMSTKS,MAX_ATMOSWFS)

      DOUBLE PRECISION L_QSUMVEC(MAXSTREAMS,MAXSTOKES,MAX_ATMOSWFS)
      DOUBLE PRECISION L_QDIFVEC(MAXSTREAMS,MAXSTOKES,MAX_ATMOSWFS)
      DOUBLE PRECISION HELP_TEMPO(0:MAXMOMENTS,MAXSTOKES)

C  help variables

      DOUBLE PRECISION HELP, WSUM, WDIF, SECBAR
      DOUBLE PRECISION TM1, TM2, TM3, L_WVEC_D, L_WVEC_U
      DOUBLE PRECISION SUM, GSUM, S_TPA, S_TMA, TPA, TMA

      INTEGER          I, J, I1, L, N, M, INFO, IB, LVARY
      INTEGER          O1, O2, O3, IR, JC, IROW, JCOL
      INTEGER          Q, K, K_PARAMETERS
      LOGICAL          DO_FIRST

C  Start of code
C  -------------

C  Set layer, fourier, beam index

      N  = GIVEN_LAYER
      M  = FOURIER
      IB = IBEAM

C  safety

      SECBAR = ZERO

C  set up driving vector (using saved results)
C  This must be done regardless of whether layer N is varying or not.

      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        DO O1 = 1, NSTOKES
          IROW = IR + O1
          QSAVE(IROW) = QDIFVEC_SAVE(I,O1) +
     &                 TWO * AVERAGE_SECANT(N,IB) * QVEC_SAVE(IROW)
        ENDDO
      ENDDO

C  Linearization for layer N is in two parts:
C    1A. Linearization due to variations in Layer N itself (Profiles)
C    1B. Linearization due to columns    
C    2. Linearization due to variations in layers K < N (Profiles)

C  Part 1.
C  =======

C  Set layer to vary

      IF ( DO_PROFILE_LINEARIZATION ) THEN
        LVARY = N
      ELSE IF ( DO_COLUMN_LINEARIZATION ) THEN
        LVARY = 0
      ENDIF

C  Only a solution if the layer is active and not below Cutoff.

      DO_FIRST = ( N .LE. LAYER_PIS_CUTOFF(IB) ) .AND.
     &                DO_LAYER_SCATTERING(M,N)

C  No particular solution beyond the cutoff layer.
C  Or no scattering in this layer
C    [ Zero the boundary layer values and start Part 2 )

      IF (.NOT. DOVARY .OR. .NOT. DO_FIRST  ) THEN
        DO I = 1, NSTREAMS_2
         DO O1 = 1, NSTOKES
          DO Q = 1, N_PARAMETERS
            L_BVEC(I,O1,N,LVARY,Q) = ZERO
          ENDDO
         ENDDO
        ENDDO
        GO TO 2222
      ENDIF

C  solar zenith cosine for this layer

      SECBAR = AVERAGE_SECANT(N,IB)

C  For each varying parameter

      DO Q = 1, N_PARAMETERS

C  Auxiliary matrix for Q functions

        DO L = M, NMOMENTS
          DO O1 = 1, NSTOKES
            SUM = ZERO
            DO O2 = 1, NSTOKES
              GSUM = ZERO
              DO O3 = 1, NSTOKES
                GSUM = GSUM + PI_X0P(L,IBEAM,N,O2,O3) * DFLUX(O3)
              ENDDO
              SUM = SUM + L_OMEGA_GREEK(L,N,O1,O2,Q) * GSUM
            ENDDO
            HELP_TEMPO(L,O1) = SUM
          ENDDO
        ENDDO

C  Set up linearized sum and difference vectors for Beam source terms

        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES
            S_TPA = ZERO
            S_TMA = ZERO
            DO L = M, NMOMENTS
              TPA = ZERO
              TMA = ZERO
              DO O2 = 1, NSTOKES
                TPA = TPA + PI_XQP(L,I,O1,O2)      * HELP_TEMPO(L,O2)
                TMA = TMA + PI_XQM_POST(L,I,O1,O2) * HELP_TEMPO(L,O2)
              ENDDO
              S_TPA = S_TPA + TPA
              S_TMA = S_TMA + TMA
            ENDDO
            L_QSUMVEC(I,O1,Q) =  ( S_TPA + S_TMA ) / QUAD_STREAMS(I)
            L_QDIFVEC(I,O1,Q) =  ( S_TPA - S_TMA ) / QUAD_STREAMS(I)
          ENDDO
        ENDDO

C   setup linearized RHS vector
C  ( use results from the original solution )
C    WARNING. BE CAREFUL OF SIGNS on TM1, TM2, TM3
C     slgihtly different from the scalar model case.

        DO I = 1, NSTREAMS
         IR = NSTOKES*(I-1)
         DO O1 = 1, NSTOKES
          IROW = IR + O1
          HELP = ZERO
          DO J = 1, NSTREAMS
           JC = NSTOKES*(J-1)
           DO O2 = 1, NSTOKES
            JCOL = JC + O2
            TM1 = L_EIGENMAT(IROW,JCOL,N,Q) * QVEC_SAVE(JCOL)
            TM2 =   DAB(I,J,O1,O2,N)   * L_QSUMVEC(J,O2,Q)
            TM3 = L_DAB(I,J,O1,O2,N,Q) *   QSUMVEC_SAVE(J,O2)
            HELP = HELP - TM1 + TM2 + TM3
c            HELP = HELP - TM1 - TM2 - TM3
           ENDDO
          ENDDO
          L_QVEC(IROW,Q)  = HELP + L_QDIFVEC(I,O1,Q) * SECBAR
         ENDDO
        ENDDO

C  end parameter loop

      ENDDO

C  additional terms for the quasi-spherical case
C  ( layers greater than one )

      IF ( .NOT. DO_PLANE_PARALLEL ) THEN
        IF ( N.GT.1 ) THEN
          DO Q = 1, N_PARAMETERS
            DO I = 1, NSTREAMS
              IR = NSTOKES*(I-1)
              DO O1 = 1, NSTOKES
                IROW = IR + O1
                HELP = QSAVE(IROW) * L_AVERAGE_SECANT(N,LVARY,IB,Q)
                L_QVEC(IROW,Q) = L_QVEC(IROW,Q) + HELP
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF

C  Solve problem by back substitution for all of the RHS vectors
C  ( uses L_U decomposition of QMAT_SAVE from original solution )

      CALL DGETRS
     &        ('N',NSTKS_NSTRMS,N_PARAMETERS,QMAT_SAVE,
     &          MAXSTRMSTKS,QPIVOT,L_QVEC,MAXSTRMSTKS,INFO)

      IF ( INFO .NE. 0 ) THEN
        CALL VLIDORT_LAPACK_ERROR
     I ( INFO, N, ' Beam P.I. linearization (N,N) (DGETRS)', STATUS )
        IF ( STATUS .EQ. VLIDORT_SERIOUS ) RETURN
      ENDIF

C  assign solutions for the quasi-spherical case, N > 1
C    Note change of sign with HELP, this is because SAB as
C    defined in VLIDORT is the negative of SAB in the scalar model.

      IF ( .NOT. DO_PLANE_PARALLEL .AND. N.GT.1 ) THEN

        DO I = 1, NSTREAMS
         IR = NSTOKES*(I-1)
         I1 = I + NSTREAMS
         DO O1 = 1, NSTOKES
          IROW = IR + O1
          TM3 = - QDIF_SAVE(IROW) / SECBAR
          DO Q = 1, N_PARAMETERS
           HELP = ZERO
           DO J = 1, NSTREAMS
            JC = NSTOKES*(J-1)
            DO O2 = 1, NSTOKES
             JCOL = JC + O2
             TM1 =   SAB(I,J,O1,O2,N)   * L_QVEC(JCOL,Q)
             TM2 = L_SAB(I,J,O1,O2,N,Q) *   QVEC_SAVE(JCOL)
             HELP = HELP + TM1 + TM2
            ENDDO
           ENDDO
           WSUM = L_QVEC(IROW,Q)
           TM2 = ( HELP - L_QSUMVEC(I,O1,Q) ) / SECBAR
           WDIF = L_AVERAGE_SECANT(N,LVARY,IB,Q) * TM3 + TM2
           L_WVEC_D = HALF * ( WSUM + WDIF )
           L_WVEC_U = HALF * ( WSUM - WDIF )
           L_BVEC(I,O1,N,LVARY,Q)  = L_WVEC_D
           L_BVEC(I1,O1,N,LVARY,Q) = L_WVEC_U
           IF ( O1 .GT. 2 ) L_BVEC(I1,O1,N,LVARY,Q) = - L_WVEC_U
C  Use this debug code
c         if(q.eq.1.and.i.eq.4.and.n.eq.24) write(39,'(4i3,1p2e20.10)')
c     &           m,ibeam,i,O1,BVEC(I,O1,N),L_BVEC(I,O1,N,LVARY,Q)
          ENDDO
         ENDDO
        ENDDO

C  assign solutions for plane/parallel & quasi-spherical case N = 1

      ELSE

        DO I = 1, NSTREAMS
         IR = NSTOKES*(I-1)
         I1 = I + NSTREAMS
         DO O1 = 1, NSTOKES
          IROW = IR + O1
          DO Q = 1, N_PARAMETERS
           HELP = ZERO
           DO J = 1, NSTREAMS
            JC = NSTOKES*(J-1)
            DO O2 = 1, NSTOKES
             JCOL = JC + O2
             TM1 =   SAB(I,J,O1,O2,N)   * L_QVEC(JCOL,Q)
             TM2 = L_SAB(I,J,O1,O2,N,Q) *   QVEC_SAVE(JCOL)
             HELP = HELP - TM1 - TM2
            ENDDO
           ENDDO
           WSUM = L_QVEC(IROW,Q)
           WDIF = ( - HELP - L_QSUMVEC(I,O1,Q) ) / SECBAR
           L_WVEC_D = HALF * ( WSUM + WDIF )
           L_WVEC_U = HALF * ( WSUM - WDIF )
           L_BVEC(I,O1,N,LVARY,Q)  = L_WVEC_D
           L_BVEC(I1,O1,N,LVARY,Q) = L_WVEC_U
           IF ( O1 .GT. 2 ) L_BVEC(I1,O1,N,LVARY,Q) = - L_WVEC_U
          ENDDO
         ENDDO
        ENDDO
        
      ENDIF


C  Old code: Toggle to L_BVEC

c      DO Q = 1, N_PARAMETERS
c       DO I = 1, NSTREAMS
c        I1 = I + NSTREAMS
c        DO O1 = 1, NSTOKES
c         L_BVEC(I,O1,N,LVARY,Q) = L_WVEC(I,O1,N,LVARY,Q)
c         HELP = ZERO
c         DO O2 = 1, NSTOKES
c          HELP = HELP + DMAT(O1,O2) * L_WVEC(I1,O2,N,LVARY,Q)
c         ENDDO
c         L_BVEC(I1,O1,N,LVARY,Q) = HELP
c        ENDDO
c       ENDDO
c      ENDDO


C  Part 2.
C  =======

C  Continuation point

 2222 CONTINUE

C  Only for the pseudo-spherical case, profile weighting functions
C  Also not required for the column weighting functions

      IF ( DO_PLANE_PARALLEL )       RETURN
      IF ( DO_COLUMN_LINEARIZATION ) RETURN

C  No particular solution beyond the cutoff layer.
C    [ Zero the boundary layer values and exit )

      IF ( .NOT. DO_FIRST ) THEN
        DO K = 1, N - 1
          K_PARAMETERS = LAYER_VARY_NUMBER(K)
          DO I = 1, NSTREAMS_2
            DO O1 = 1, NSTOKES
              DO Q = 1, K_PARAMETERS
                L_BVEC(I,O1,N,K,Q) = ZERO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        RETURN
      ENDIF

C  Loop over layers K above N

      DO K = 1, N - 1

C  If there is a varying layer

        IF ( LAYER_VARY_FLAG(K) ) THEN

C  number of varying parameters for this layer

          K_PARAMETERS = LAYER_VARY_NUMBER(K)

C  Set up vector

          DO I = 1, NSTREAMS
           IR = NSTOKES*(I-1)
           DO O1 = 1, NSTOKES
            IROW = IR + O1
            DO Q = 1, K_PARAMETERS
             L_QVEC(IROW,Q) = QSAVE(IROW) * L_AVERAGE_SECANT(N,K,IB,Q)
            ENDDO
           ENDDO
          ENDDO

C  Solve problem by back substitution for all of the RHS vectors
C  ( uses L_U decomposition of QMAT_SAVE from original solution )

          CALL DGETRS
     &        ('N',NSTKS_NSTRMS,K_PARAMETERS,QMAT_SAVE,
     &          MAXSTRMSTKS,QPIVOT,L_QVEC,MAXSTRMSTKS,INFO)

          IF ( INFO .NE. 0 ) THEN
            CALL VLIDORT_LAPACK_ERROR
     I ( INFO, N, ' Beam P.I. linearization (N,K) (DGETRS)', STATUS )
            IF ( STATUS .EQ. VLIDORT_SERIOUS ) RETURN
          ENDIF

C  assign linearized solutions for layer N due to variations in layer K
C    Note change of sign with HELP, this is because SAB as
C    defined in VLIDORT is the negative of SAB in the scalar model.

          DO I = 1, NSTREAMS
           IR = NSTOKES*(I-1)
           I1 = I + NSTREAMS
           DO O1 = 1, NSTOKES
            IROW = IR + O1
            TM3 = - QDIF_SAVE(IROW) / SECBAR
            DO Q = 1, K_PARAMETERS
             HELP = ZERO
             DO J = 1, NSTREAMS
              JC = NSTOKES*(J-1)
              DO O2 = 1, NSTOKES
               JCOL = JC + O2
               HELP = HELP + SAB(I,J,O1,O2,N) * L_QVEC(JCOL,Q)
              ENDDO
             ENDDO
             WSUM = L_QVEC(IROW,Q)
             WDIF = L_AVERAGE_SECANT(N,K,IB,Q) * TM3 + (HELP/SECBAR)
             L_WVEC_D = HALF * ( WSUM + WDIF )
             L_WVEC_U = HALF * ( WSUM - WDIF )
             L_BVEC(I,O1,N,K,Q)  = L_WVEC_D
             L_BVEC(I1,O1,N,K,Q) = L_WVEC_U
             IF ( O1 .GT. 2 ) L_BVEC(I1,O1,N,K,Q) = - L_WVEC_U
            ENDDO
           ENDDO
          ENDDO

C  Old code.Toggle to L_BVEC

c          DO Q = 1, K_PARAMETERS
c           DO I = 1, NSTREAMS
c            I1 = I + NSTREAMS
c            DO O1 = 1, NSTOKES
c             L_BVEC(I,O1,N,K,Q) = L_WVEC(I,O1,N,K,Q)
c             HELP = ZERO
c             DO O2 = 1, NSTOKES
c              HELP = HELP + DMAT(O1,O2) * L_WVEC(I1,O2,N,K,Q)
c             ENDDO
c             L_BVEC(I1,O1,N,K,Q) = HELP
c            ENDDO
c           ENDDO
c          ENDDO

C  end K-layer loop

        ENDIF
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE VLIDORT_L_UBEAM_SOLUTION
     I    ( GIVEN_LAYER, FOURIER, IBEAM, DO_VARY, N_PARAMETERS )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include file of setup and solution stuff (inputs to module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'

C  include file of linearization inputs

      INCLUDE '../includes/VLIDORT_L_INPUTS.VARS'
     
C  include file of main solution linearization variables (output)

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'

C  include file of setup linearization variables (output)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index, Fourier index, parameter number (inputs)

      INTEGER          GIVEN_LAYER
      INTEGER          FOURIER
      INTEGER          IBEAM

C  variation flag

      LOGICAL          DO_VARY
      INTEGER          N_PARAMETERS

C  Local variables
C  ---------------

      LOGICAL          DO_LAYER, DO_FIRST
      INTEGER          IB, UM, L, N, M, O1, O2, O3, J, J1
      INTEGER          Q, K, K_PARAMETERS, LVARY

      DOUBLE PRECISION SUM1, SUM2, T1, T2, H, HH
      DOUBLE PRECISION SP, SM, SPS, SPT, SGW(MAXSTOKES)

      DOUBLE PRECISION L_HELP_Q1(0:MAXMOMENTS,MAXSTOKES)
      DOUBLE PRECISION L_HELP_Q2(0:MAXMOMENTS,MAXSTOKES)
      DOUBLE PRECISION L_HELPSTOKES_BEAM(MAXSTOKES)
      DATA             SGW / 1.0d0, 1.0d0, -1.0d0, -1.0d0 /

C  Layer and Fourier

      N  = GIVEN_LAYER
      M  = FOURIER
      IB = IBEAM
      DO_LAYER = .FALSE.
      LVARY = 0

C  Part 1. Linearizations for quantities varying in Layer N
C  ========================================================

C  Profiles and Columns

      IF ( DO_PROFILE_LINEARIZATION ) THEN
        LVARY = N
      ELSE IF ( DO_COLUMN_LINEARIZATION ) THEN
        LVARY = 0
      ENDIF

C  Check existence
C  ---------------

C  Only a solution if the layer is active and not below Cutoff.

      DO_FIRST = ( N .LE. LAYER_PIS_CUTOFF(IB) ) .AND.
     &                DO_LAYER_SCATTERING(M,N)

C  If no solution or no variation, zero output and go to Part 2

      IF ( .NOT. DO_VARY .OR. .NOT. DO_FIRST ) THEN
        DO Q = 1, N_PARAMETERS
         IF ( DO_UPWELLING ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              L_UPAR_UP_1(UM,O1,N,Q)       = ZERO
              L_UPAR_UP_2(UM,O1,N,LVARY,Q) = ZERO
            ENDDO
          ENDDO
         ENDIF
         IF ( DO_DNWELLING ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              L_UPAR_DN_1(UM,O1,N,Q)       = ZERO
              L_UPAR_DN_2(UM,O1,N,LVARY,Q) = ZERO
            ENDDO
          ENDDO
         ENDIF
        ENDDO
        GO TO 2222
      ENDIF

C  existence flag

      DO_LAYER = ( STERM_LAYERMASK_UP(N) .OR.
     &             STERM_LAYERMASK_DN(N) )

C  start parameter loop

      DO Q = 1, N_PARAMETERS

C  first function
C    (linearization of primary scattering of beam)

       IF ( DO_LAYER ) THEN
        DO L = M, NMOMENTS
         DO O1 = 1, NSTOKES
          SPS = ZERO
          DO O2 = 1, NSTOKES
            SP = ZERO
            DO O3 = 1, NSTOKES
              SP = SP + PI_X0P(L,IBEAM,N,O2,O3) * DFLUX(O3)
            ENDDO
            SPS = SPS + L_OMEGA_GREEK(L,N,O1,O2,Q) * SP
          ENDDO
          L_HELP_Q1(L,O1) = SPS
         ENDDO
        ENDDO
       ENDIF

C  For each moment, do inner sum over computational angles
C   (Linearization of Diffuse scattering of beam)

       IF ( DO_LAYER ) THEN
        DO L = M, NMOMENTS
         DO O1 = 1, NSTOKES
          SPS = ZERO
          DO J = 1, NSTREAMS
           J1 = J + NSTREAMS
           SP = ZERO
           SM = ZERO
           H = QUAD_HALFWTS(J)
           DO O2 = 1, NSTOKES
            HH = H * SGW(O2)
            SP = SP + H  * PI_XQP    (L,J,O1,O2)*L_BVEC(J,O2,N,LVARY,Q)
            SM = SM + HH * PI_XQM_PRE(L,J,O1,O2)*L_BVEC(J1,O2,N,LVARY,Q)
           ENDDO
           SPS = SPS + SP + SM
          ENDDO
          L_HELPSTOKES_BEAM(O1) = SPS
         ENDDO
         DO O1 = 1, NSTOKES
          SPT = ZERO
          DO O2 = 1, NSTOKES
           SPT = SPT + OMEGA_GREEK(L,N,O1,O2)   * L_HELPSTOKES_BEAM(O2)
     &             + L_OMEGA_GREEK(L,N,O1,O2,Q) * HELPSTOKES_BEAM(L,O2)
          ENDDO
          L_HELP_Q2(L,O1) = SPT
         ENDDO
        ENDDO
       ENDIF

C  Now sum over all harmonic contributions (Upwelling)
C    Direct and integrated contributions

       IF ( DO_UPWELLING .AND. STERM_LAYERMASK_UP(N) ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
         DO O1 = 1, NSTOKES
          T1 = ZERO
          T2 = ZERO
          DO L = M, NMOMENTS
           SUM1 = ZERO
           SUM2 = ZERO
           DO O2 = 1, NSTOKES
            SUM1 = SUM1 + L_HELP_Q1(L,O2)*PI_XUM(L,UM,O1,O2)
            SUM2 = SUM2 + L_HELP_Q2(L,O2)*PI_XUM(L,UM,O1,O2)
           ENDDO
           T1 = T1 + SUM1
           T2 = T2 + SUM2
          ENDDO
          L_UPAR_UP_1(UM,O1,N,Q)       = T1
          L_UPAR_UP_2(UM,O1,N,LVARY,Q) = T2
         ENDDO
        ENDDO
       ENDIF

C  Now sum over all harmonic contributions (Downwelling)
C    Direct and integrated contributions

       IF ( DO_DNWELLING .AND. STERM_LAYERMASK_DN(N) ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
         DO O1 = 1, NSTOKES
          T1 = ZERO
          T2 = ZERO
          DO L = M, NMOMENTS
           SUM1 = ZERO
           SUM2 = ZERO
           DO O2 = 1, NSTOKES
            SUM1 = SUM1 + L_HELP_Q1(L,O2)*PI_XUP(L,UM,O1,O2)
            SUM2 = SUM2 + L_HELP_Q2(L,O2)*PI_XUP(L,UM,O1,O2)
           ENDDO
           T1 = T1 + SUM1
           T2 = T2 + SUM2
          ENDDO
          L_UPAR_DN_1(UM,O1,N,Q)       = T1
          L_UPAR_DN_2(UM,O1,N,LVARY,Q) = T2
         ENDDO
        ENDDO
       ENDIF

C  end parameter loop

      ENDDO

C  Part 2.
C  =======

C  Continuation point

 2222 CONTINUE

C  Only these extra variations when the following conditions
C  are not satisfied (pseudo-spherical, layer > 1)
C  Profiles only

      IF ( N .EQ. 1 )                RETURN
      IF ( DO_PLANE_PARALLEL )       RETURN
      IF ( DO_COLUMN_LINEARIZATION ) RETURN

C  No particular solution beyond the cutoff layer.
C   Or no scattering in this layer
C    [ Zero the boundary layer values and exit )

      IF ( .NOT. DO_FIRST ) THEN
       DO K = 1, N - 1
        K_PARAMETERS = LAYER_VARY_NUMBER(K)
        DO Q = 1, K_PARAMETERS
         IF ( DO_UPWELLING ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              L_UPAR_UP_2(UM,O1,N,K,Q) = ZERO
            ENDDO
          ENDDO
         ENDIF
         IF ( DO_DNWELLING ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              L_UPAR_DN_2(UM,O1,N,K,Q) = ZERO
            ENDDO
          ENDDO
         ENDIF
        ENDDO
       ENDDO
       RETURN
      ENDIF

C  start loop over all layers above N

      DO K = 1, N - 1

C  Start parameter loop -- only do if layer K has some variation

       IF ( LAYER_VARY_FLAG(K) ) THEN
        K_PARAMETERS = LAYER_VARY_NUMBER(K)
        DO Q = 1, K_PARAMETERS

C  For each moment, do inner sum over computational angles
C   (Linearization of Diffuse scattering of beam)

         IF ( DO_LAYER ) THEN
          DO L = M, NMOMENTS
           DO O1 = 1, NSTOKES
            SPS = ZERO
            DO J = 1, NSTREAMS
             J1 = J + NSTREAMS
             SP = ZERO
             SM = ZERO
             H = QUAD_HALFWTS(J)
             DO O2 = 1, NSTOKES
              HH = H * SGW(O2)
              SP = SP + H  * PI_XQP    (L,J,O1,O2)*L_BVEC(J,O2,N,K,Q)
              SM = SM + HH * PI_XQM_PRE(L,J,O1,O2)*L_BVEC(J1,O2,N,K,Q)
             ENDDO
             SPS = SPS + SP + SM
            ENDDO
            L_HELPSTOKES_BEAM(O1) = SPS
           ENDDO
           DO O1 = 1, NSTOKES
            SPT = ZERO
            DO O2 = 1, NSTOKES
             SPT = SPT + OMEGA_GREEK(L,N,O1,O2) * L_HELPSTOKES_BEAM(O2)
            ENDDO
            L_HELP_Q2(L,O1) = SPT
           ENDDO
          ENDDO
         ENDIF

C  Now sum over all harmonic contributions (Upwelling)
C    Diffuse (integrated) contributions only

         IF ( DO_UPWELLING .AND. STERM_LAYERMASK_UP(N) ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
           DO O1 = 1, NSTOKES
            T2 = ZERO
            DO L = M, NMOMENTS
             SUM2 = ZERO
             DO O2 = 1, NSTOKES
              SUM2 = SUM2 + L_HELP_Q2(L,O2)*PI_XUM(L,UM,O1,O2)
             ENDDO
             T2 = T2 + SUM2
            ENDDO
            L_UPAR_UP_2(UM,O1,N,K,Q) = T2
           ENDDO
          ENDDO
         ENDIF

C  Now sum over all harmonic contributions (Downwelling)
C    Diffuse (integrated) contributions only

         IF ( DO_DNWELLING .AND. STERM_LAYERMASK_DN(N) ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
           DO O1 = 1, NSTOKES
            T2 = ZERO
            DO L = M, NMOMENTS
             SUM2 = ZERO
             DO O2 = 1, NSTOKES
              SUM2 = SUM2 + L_HELP_Q2(L,O2)*PI_XUP(L,UM,O1,O2)
             ENDDO
             T2 = T2 + SUM2
            ENDDO
            L_UPAR_DN_2(UM,O1,N,K,Q) = T2
           ENDDO
          ENDDO
         ENDIF

C  end parameter and K-layer loop

        ENDDO
       ENDIF
      ENDDO

C  Finish

      RETURN
      END

