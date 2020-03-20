C ###########################################################
C #                                                         #
C #                    THE LIDORT FAMILY                    #
C #                                                         #
C #      (LInearized Discrete Ordinate Radiative Transfer)  #
C #       --         -        -        -         -          #
C #                                                         #
C ###########################################################

C ###########################################################
C #                                                         #
C #  Author :      Robert. J. D. Spurr                      #
C #                                                         #
C #  Address :     RT Solutions, Inc.                       #
C #                9 Channing Street                        #
C #                Cambridge, MA 02138, USA                 #
C #                                                         #
C #  Tel:          (617) 492 1183                           #
C #  Email :        rtsolutions@verizon.net                 #
C #                                                         #
C #  This Version :   3.3                                   #
C #  Release Date :   September 2007                        #
C #                                                         #
C #       NEW: THERMAL SUPPLEMENT INCLUDED    (3.2)         #
C #       NEW: OUTGOING SPHERICITY CORRECTION (3.2)         #
C #       NEW: TOTAL COLUMN JACOBIANS         (3.3)         #
C #                                                         #
C ###########################################################

C    #####################################################
C    #                                                   #
C    #   This Version of LIDORT comes with a GNU-style   #
C    #   license. Please read the license carefully.     #
C    #                                                   #
C    #####################################################

C ###############################################################
C #                                                             #
C # Subroutines in this Module                                  #
C #                                                             #
C #              LIDORT_HOM_SOLUTION                            #
C #              LIDORT_HOM_EIGENTRANS                          #
C #              LIDORT_HOM_USERSOLUTION                        #
C #              LIDORT_QBEAM_SOLUTION                          #
C #              LIDORT_QBEAM_USERSOLUTION                      #
C #              LIDORT_GBEAM_SOLUTION                          #
C #              LIDORT_GBEAM_USERSOLUTION                      #
C #              LIDORT_LEGENDRE_SETUP                          #
C #              LIDORT_USERLEGENDRE_SETUP                      #
C #                                                             #
C ###############################################################

      SUBROUTINE LIDORT_HOM_SOLUTION
     I    ( GIVEN_LAYER, FOURIER, 
     O      STATUS )

C  Numerical solution of Eigenproblem.

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file of Set up stuff (input/output to this module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  include file of main solution variables (output to this module)

      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index and Fourier number (inputs)

      INTEGER          GIVEN_LAYER
      INTEGER          FOURIER

C  output status

      INTEGER          STATUS

C  Local variables
C  ---------------

C  local matrices for eigenvalue computation

      DOUBLE PRECISION
     &      EIGENMAT(MAXSTREAMS,MAXSTREAMS)

C  (output from Eigenpackage module ASYMTX)

      DOUBLE PRECISION KSQ(MAXSTREAMS), WK(MAXSTREAMS_2)
      DOUBLE PRECISION EVEC(MAXSTREAMS,MAXSTREAMS)
      INTEGER          IER
      LOGICAL          ASYMTX_FAILURE
      CHARACTER*70     MESSAGE

C  Miscellaneous local variables

      INTEGER          I, J, I1, L, N, M, AA, K
      DOUBLE PRECISION DP, DM, SUM, FAC, KVAL, NORM, XINV
    
C  local variables for error tracing

      CHARACTER*(70)   MAIL, TRACE
      CHARACTER*3      CN, CI

C  initialise status

      STATUS = LIDORT_SUCCESS

C  Layer and Fourier

      N = GIVEN_LAYER
      M = FOURIER

C  When the solution saving option is set, then if there is no
C  scattering in this layer, then the eigenvectors are just the
C  discrete ordinates, the eigenvectors are unit vectors in the
C  discrete ordinate directions.

      IF ( DO_SOLUTION_SAVING ) THEN
        IF ( .NOT. DO_LAYER_SCATTERING(M,N) ) THEN
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            XINV = ONE / X(I)
            DO AA = 1, NSTREAMS
              XPOS(I,AA,N)  = ZERO
              XPOS(I1,AA,N) = ZERO
              XNEG(I1,AA,N) = ZERO
              XNEG(I,AA,N)  = ZERO
            ENDDO
            KEIGEN(I,N)  = XINV
            XPOS(I,I,N)  = ONE
            XNEG(I1,I,N) = ONE
          ENDDO
          GO TO 3456
        ENDIF
      ENDIF

C  Scattering solutions
C  ====================

C  Construct Eigenmatrix
C  ---------------------

C  zero the Eigenmatrix

      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          EIGENMAT(I,J) = ZERO
        ENDDO
      ENDDO

C  Develop Sum and Difference matrices

      DO I = 1, NSTREAMS
        XINV = ONE/X(I)
        DO J = 1, NSTREAMS
          FAC = XINV * HALF * A(J)
          DP = ZERO
          DM = ZERO
          DO L = M, NMOMENTS
            DP = DP + PLMI_PLMJ_P(I,J,L) * OMEGA_MOMS(N,L)
            DM = DM + PLMI_PLMJ_M(I,J,L) * OMEGA_MOMS(N,L)
          ENDDO
          SAB(I,J,N) = FAC * ( DP + DM ) 
          DAB(I,J,N) = FAC * ( DP - DM )
        ENDDO
        SAB(I,I,N) = SAB(I,I,N) - XINV
        DAB(I,I,N) = DAB(I,I,N) - XINV
      ENDDO

C  Compute Eigenmatrix

      DO I = 1, NSTREAMS 
        DO J = 1, NSTREAMS
          SUM = ZERO
          DO K = 1, NSTREAMS
            SUM = SUM + DAB(I,K,N) * SAB(K,J,N)
          ENDDO
          EIGENMAT(I,J) = SUM
        ENDDO
      ENDDO

C  save Eigenmatrix (original is destroyed by ASMTYX)

      DO I = 1, NSTREAMS 
        DO J = 1, NSTREAMS
          EIGENMAT_SAVE(I,J,N) = EIGENMAT(I,J)
        ENDDO
c        if (n.eq.14)write(*,'(2i3,1p4e15.7)')m,i,(eigenmat(i,j),j=1,4)
      ENDDO

C  Eigensolution package
C  ---------------------

C  Let's see how we get on with the DISORT package
C  tested 19 May 1999. Gives same results as Analytical approach.

C  second test using DGEEV (LAPACK module) - Worked against ASYMTX.
C  However, DGEEV is 1.7 - 2.0 times slower because it must look for
C  complex roots. [ASYMTX was specially written to avoid this search].
C  Also first call to DGEEV is 20 times slower. Tested June 24, 1999.
C  Conclusion stick with ASYMTX for now.

      CALL  ASYMTX
     I         ( EIGENMAT, NSTREAMS, MAXSTREAMS, MAXSTREAMS,
     O           EVEC, KSQ, IER, WK,
     O           MESSAGE, ASYMTX_FAILURE )

C  error tracing

      IF ( ASYMTX_FAILURE .OR. IER.GT.0 ) THEN
        WRITE(CI,'(I3)')IER
        WRITE(CN,'(I3)')N
        MAIL  = 'eigenvalue '//CI//' has not converged'
        IF ( ASYMTX_FAILURE ) MAIL = MESSAGE
        TRACE ='ASYMTX error, HOMOGENEOUS_SOLUTION, Layer='//CN
        STATUS = LIDORT_SERIOUS
        CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
        RETURN
      ENDIF

C  second test using DGEEV (LAPACK module). Here is what to use.
C        CALL DGEEV
C     I      ( 'N', 'V', NSTREAMS, EIGENMAT,
C     O        MAXSTREAMS, REAL_KSQ, IMAG_KSQ,
C     O        LEFT_EVEC, MAXSTREAMS, EVEC, MAXSTREAMS,
C     W        WORK, LWORK, LAPACK_INFO )

C  Find solution eigenvectors XPOS, XNEG for all eigenvalues
C  ---------------------------------------------------------

      DO AA = 1, NSTREAMS

C  Store positive values in output array for each layer

        KVAL = DSQRT(KSQ(AA))
        KEIGEN(AA,N) = KVAL

C  Normalize eigenvectors to 1

        NORM = ZERO
        DO I = 1, NSTREAMS
          NORM = NORM + EVEC(I,AA)*EVEC(I,AA)
        ENDDO
        NORM = DSQRT(NORM)

C  Find normalized eigenvector EIGENVEC_SAVE

        DO I = 1, NSTREAMS
          EIGENVEC_SAVE(I,AA) = EVEC(I,AA)/NORM
        ENDDO

C  Find difference eigenvector DIFVEC (Siewert's notation)

        DO I = 1, NSTREAMS
          SUM = ZERO
          DO K = 1, NSTREAMS
            SUM = SUM - SAB(I,K,N) * EIGENVEC_SAVE(K,AA)
          ENDDO
          DIFVEC_SAVE(I,AA) = SUM / KVAL
          ENDDO

C  assign original evectors; first N are "DOWN", last N are "UP" (streams)

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          XPOS(I,AA,N)  =
     &            HALF * ( EIGENVEC_SAVE(I,AA) + DIFVEC_SAVE(I,AA) )
          XPOS(I1,AA,N) =
     &            HALF * ( EIGENVEC_SAVE(I,AA) - DIFVEC_SAVE(I,AA) )
        ENDDO

C  debug
C   --linearization check
c        if ( do_debug_write ) then
c          write(88,'(3i4,1p20e15.7)')M,N,AA,(XPOS(I,AA,N),I=1,NSTREAMS)
c        endif

C  Use symmetry properties to set -ve eigenvectors

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          XNEG(I1,AA,N) = XPOS(I,AA,N)
          XNEG(I,AA,N)  = XPOS(I1,AA,N)
        ENDDO

C  End eigenstream loop

      ENDDO

C  control point (for debug purposed only)

 3456 CONTINUE

C  debug

c       k = 97
c       if (do_fdtest)k=98
c        if ( m.lt.3.and.n.gt.0) then
c        write(k,'(2i5)')M,N
c        DO AA = 1, NSTREAMS
c          WRITE(k,'(A3,I5,1p12e20.8)')'***',AA,
c     *          (XPOS(I,AA,N),I=1,NSTREAMS/2)
c        ENDDO
c        endif

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_HOM_EIGENTRANS
     I    ( FOURIER )

C  Eigenproblem, transmittance matrices

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file of main solution variables (input)

      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  include file of Setup stuff (input/output to this module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  subroutine arguments
C  --------------------

C  Fourier number (input)

      INTEGER          FOURIER

C  Local variables
C  ---------------

C  Miscellaneous local variables

      INTEGER          N, M, AA, UT, UTA
      DOUBLE PRECISION HELP, TAU_UP, TAU_DN
    
C  Fourier

      M = FOURIER

C  Layer transmittances for the eigenvalues
C  ========================================

C  start the layer loop

      DO N = 1, NLAYERS

C  When the solution saving option is set, then if there is no
C  scattering in this layer, then transmittances are just the
C  discrete ordinate transmittances.

        IF ( DO_SOLUTION_SAVING .AND.
     &         .NOT.DO_LAYER_SCATTERING(M,N) ) THEN

          DO AA = 1, NSTREAMS
            T_DELT_EIGEN(AA,N) = T_DELT_DISORDS(AA,N)
          ENDDO

C  Otherwise, get the full set of Eigenstream transmittance factors

        ELSE
          DO AA = 1, NSTREAMS
            HELP = KEIGEN(AA,N)*DELTAU_VERT(N)
            IF ( HELP .GT. MAX_TAU_QPATH ) THEN
              T_DELT_EIGEN(AA,N) = ZERO
            ELSE
              T_DELT_EIGEN(AA,N) = DEXP(-HELP)
            ENDIF
          ENDDO
        ENDIF

C  end layer loop

      ENDDO

C  Eigenstream transmittance factors for partial layers
C  ----------------------------------------------------

C  Code completed by R. Spurr, RT SOLUTIONS Inc. 30 August 2005

      DO UTA = 1, N_OUT_USERTAUS
       IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN
        UT = OFFGRID_UTAU_OUTINDEX(UTA)
        N  = OFFGRID_UTAU_LAYERIDX(UT)

C  When the solution saving option is set, then if there is no
C  scattering in this layer, then transmittances are just the
C  discrete ordinate transmittances.

        IF ( DO_SOLUTION_SAVING .AND.
     &         .NOT.DO_LAYER_SCATTERING(M,N) ) THEN

          DO AA = 1, NSTREAMS
            T_UTDN_EIGEN(AA,UT) = T_DISORDS_UTDN(AA,UT)
            T_UTUP_EIGEN(AA,UT) = T_DISORDS_UTUP(AA,UT)
          ENDDO

C  Otherwise, Compute the Eigenstream transmittance factors

        ELSE

          TAU_DN = OFFGRID_UTAU_VALUES(UT)
          TAU_UP = DELTAU_VERT(N) - TAU_DN
          DO AA = 1, NSTREAMS
            HELP = KEIGEN(AA,N) * TAU_DN
            IF ( HELP .GT. MAX_TAU_QPATH ) THEN
              T_UTDN_EIGEN(AA,UT) = ZERO
            ELSE
              T_UTDN_EIGEN(AA,UT) = DEXP(-HELP)
            ENDIF
            HELP = KEIGEN(AA,N) * TAU_UP
            IF ( HELP .GT. MAX_TAU_QPATH ) THEN
              T_UTUP_EIGEN(AA,UT) = ZERO
            ELSE
              T_UTUP_EIGEN(AA,UT) = DEXP(-HELP)
            ENDIF
          ENDDO

        ENDIF

C  end loop over off-grid optical depths

       ENDIF
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_HOM_NORMS

C  Eigenproblem, solution norms for Green's function

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file of main solution variables (input/output)

      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  Local variables
C  ---------------

C  Miscellaneous local variables

      INTEGER          N, AA, J, J1
      DOUBLE PRECISION T1, T2, NORM
    
C  For all layers, save the norms

      DO N = 1, NLAYERS
        DO AA = 1, NSTREAMS
          NORM = ZERO
          DO J = 1, NSTREAMS
            J1 = J + NSTREAMS
            T1 = XPOS(J,AA,N)  * XPOS(J,AA,N)
            T2 = XPOS(J1,AA,N) * XPOS(J1,AA,N)
            NORM = NORM + AX(J)*(T1-T2)
          ENDDO
          NORM_SAVED(N,AA) = NORM
        ENDDO
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_QBEAM_SOLUTION
     I    ( GIVEN_LAYER, FOURIER, IBEAM,
     O      STATUS )

C  This is the classical Chandrasekhar beam solution.
C  ( plane parallel or average secant only)
C  Linear Matrix algebra.

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file of Set up stuff (input/output to this module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  include file of main solution variables (output to this module)

      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index and Fourier number (inputs)

      INTEGER          GIVEN_LAYER
      INTEGER          FOURIER
      INTEGER          IBEAM

C  output status

      INTEGER          STATUS

C  Local variables
C  ---------------

C  help variables

      CHARACTER*70     MAIL
      INTEGER          I, J, I1, L, N, M, INFO
      DOUBLE PRECISION TP, TM, INV_X0SQ, SECBAR, XINV, F1,
     &                 HELP, WSUM, WDIF, TRANS1, TRANS2

C  Layer

      M = FOURIER
      N = GIVEN_LAYER
      F1 = FLUX_FACTOR / PI4

C  No particular solution beyond the cutoff layer
C  Or no scattering in this layer...
C  ... Zero the boundary layer values and exit

      IF ( N .GT. LAYER_PIS_CUTOFF(IBEAM) .OR.
     &         .NOT. DO_LAYER_SCATTERING(M,N) ) THEN
        DO I = 1, NSTREAMS_2
          WUPPER(I,N) = ZERO
          WLOWER(I,N) = ZERO
        ENDDO
        RETURN
      ENDIF

C  set local values

      SECBAR   = AVERAGE_SECANT(N,IBEAM)
      INV_X0SQ = SECBAR * SECBAR

C  Initialise matrix and column vector for solution

      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          QMAT_SAVE(I,J) = ZERO
        ENDDO
        QVEC_SAVE(I) = ZERO
      ENDDO

C  Set up sum and difference vectors for Beam source terms
C  ( sum vector may be required again in linearization )

      DO I = 1, NSTREAMS
        XINV = ONE / X(I)
        TP = ZERO
        TM = ZERO
        DO L = M, NMOMENTS
          TP = TP + PLMI_X0_P(I,L,N,IBEAM) * OMEGA_MOMS(N,L)
          TM = TM + PLMI_X0_M(I,L,N,IBEAM) * OMEGA_MOMS(N,L)
        ENDDO
        QSUMVEC_SAVE(I) =  F1 * ( TP + TM ) * XINV
        QDIFVEC_SAVE(I) =  F1 * ( TP - TM ) * XINV
      ENDDO

C  solution matrix for the reduced problem
C  ( matrix should be saved in the LU decomposition form)

      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          QMAT_SAVE(I,J) = EIGENMAT_SAVE(I,J,N)
        ENDDO
        QMAT_SAVE(I,I) = QMAT_SAVE(I,I) - INV_X0SQ
      ENDDO

C  RHS vector for the reduced problem
C  ( this vector will be the answer after the linear algebra solution,
C    and may be needed again if there is linearization )

      DO I = 1, NSTREAMS
        HELP = ZERO
        DO J = 1, NSTREAMS
          HELP = HELP - DAB(I,J,N)*QSUMVEC_SAVE(J)
        ENDDO
        QVEC_SAVE(I) = HELP + QDIFVEC_SAVE(I) * SECBAR
      ENDDO

C  L-U decomposition of the solution matrix

      CALL DGETRF(NSTREAMS,NSTREAMS,QMAT_SAVE,MAXSTREAMS,QPIVOT,INFO)

      IF ( INFO .NE. 0 ) THEN
        MAIL='Classical Beam solution LU decomposition (DGETRF)'
        CALL LIDORT_LAPACK_ERROR ( INFO, N, MAIL, STATUS )
        IF ( STATUS .EQ. LIDORT_SERIOUS ) RETURN
      ENDIF

C  Solution of reduced problem by back-substitution

      CALL DGETRS
     &    ('N',NSTREAMS,1,QMAT_SAVE,MAXSTREAMS,QPIVOT,
     &     QVEC_SAVE,MAXSTREAMS,INFO)

      IF ( INFO .NE. 0 ) THEN
        MAIL='Classical Beam solution back-substitution (DGETRS)'
        CALL LIDORT_LAPACK_ERROR ( INFO, N, MAIL, STATUS )
        IF ( STATUS .EQ. LIDORT_SERIOUS ) RETURN
      ENDIF

C  Assigning beam particular integral solution vector

      DO I = 1, NSTREAMS
        I1 = I + NSTREAMS
        HELP = ZERO
        DO J = 1, NSTREAMS
          HELP = HELP - SAB(I,J,N)*QVEC_SAVE(J)
        ENDDO
        QDIF_SAVE(I) = ( HELP - QSUMVEC_SAVE(I) ) / SECBAR
        WSUM = QVEC_SAVE(I)
        WDIF = QDIF_SAVE(I)
        WVEC(I,N)  = HALF * ( WSUM + WDIF )
        WVEC(I1,N) = HALF * ( WSUM - WDIF )
      ENDDO

C  Values at the layer boundaries
C  (transmittance factors have been determined in SETUPS module)

      TRANS1 = INITIAL_TRANS(N,IBEAM)
      TRANS2 = T_DELT_MUBAR(N,IBEAM) * TRANS1
      DO I = 1, NSTREAMS_2
        WUPPER(I,N) = WVEC(I,N)*TRANS1
        WLOWER(I,N) = WVEC(I,N)*TRANS2
      ENDDO

C  debug

c      IF ( DO_DEBUG_WRITE ) THEN
c       J = 97
c       IF ( DO_FDTEST ) J = 98
c       IF ( N.EQ.14.and.m.eq.1) THEN
c        DO I = 1, NSTREAMS*2
c         write(*,'(4i4,1p4e17.9)')IBEAM,M,N,I,WVEC(I,N),
c     &      INITIAL_TRANS(N,IBEAM), T_DELT_MUBAR(N,IBEAM)
c        ENDDO
c       ENDIF
c       IF ( DO_FDTEST. and. M.EQ.1) PAUSE
c      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_GBEAM_SOLUTION
     I    ( GIVEN_LAYER, FOURIER, IBEAM )

C  Green's function beam particular integral, one layer only.
C  Uses coefficient expansion of attenuation.

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file of setup stuff (input/output to this module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  include file of solution variables (input/output to this module)

      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  include file of Multiplier coefficients (input/output to this module)

      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  subroutine arguments
C  ====================

C  Given layer index and Fourier number (inputs)

      INTEGER          GIVEN_LAYER
      INTEGER          FOURIER
      INTEGER          IBEAM

C  local variables
C  ===============

      INTEGER          AA, L, I, I1, M, N
      DOUBLE PRECISION TPA, TMA, SUM_LA, SUM_LB
      DOUBLE PRECISION S_P_U, S_P_L, S_M_U, S_M_L

      DOUBLE PRECISION RHO_M, RHO_P, F1
      DOUBLE PRECISION CONST, SECBAR, WDEL, ZDEL, ZWDEL

C  initialise indices

      N = GIVEN_LAYER
      M = FOURIER

C  No particular solution beyond the cutoff layer
C  Or no scattering in this layer...
C  ... Zero the boundary layer values and exit

      IF ( N .GT. LAYER_PIS_CUTOFF(IBEAM) .OR.
     &         .NOT. DO_LAYER_SCATTERING(M,N) ) THEN
        DO I = 1, NSTREAMS_2
          WUPPER(I,N) = ZERO
          WLOWER(I,N) = ZERO
        ENDDO
        RETURN
      ENDIF

C  constants for the layer

      SECBAR = AVERAGE_SECANT(N,IBEAM)
      CONST  = INITIAL_TRANS(N,IBEAM)

C  Gamma constants

      DO AA = 1, NSTREAMS
        RHO_M = SECBAR - KEIGEN(AA,N)
        RHO_P = SECBAR + KEIGEN(AA,N)
        GAMMA_P(AA,N) = ONE / RHO_P
        GAMMA_M(AA,N) = ONE / RHO_M
      ENDDO

C  3. Optical depth integrations for the discrete ordinate solution
C     =============================================================

      WDEL    = T_DELT_MUBAR(N,IBEAM)
      DO AA = 1, NSTREAMS
        ZDEL  = T_DELT_EIGEN(AA,N)
        ZWDEL = ZDEL * WDEL
        CFUNC(AA,N)  = ( ZDEL - WDEL ) * GAMMA_M(AA,N)
        DFUNC(AA,N)  = ( ONE - ZWDEL ) * GAMMA_P(AA,N)
      ENDDO

C  4. Form quantities independent of optical depth
C     ============================================

C  set up help arrays (independent of eigenvector)

      F1 = FLUX_FACTOR / PI4
      DO I = 1, NSTREAMS
        TPA = ZERO
        TMA = ZERO
        DO L = M, NMOMENTS
          TPA = TPA + PLMI_X0_P(I,L,N,IBEAM) * OMEGA_MOMS(N,L)
          TMA = TMA + PLMI_X0_M(I,L,N,IBEAM) * OMEGA_MOMS(N,L)
        ENDDO
        DPI(I) = TPA * F1
        DMI(I) = TMA * F1
      ENDDO

C  For each eigenstream, get the terms ATERM_SAVE and BTERM_SAVE

      DO AA = 1, NSTREAMS
        SUM_LA = ZERO
        SUM_LB = ZERO
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          TPA = A(I)*(DPI(I)*XPOS(I,AA,N)+DMI(I)*XPOS(I1,AA,N))
          TMA = A(I)*(DMI(I)*XPOS(I,AA,N)+DPI(I)*XPOS(I1,AA,N))
          SUM_LA  = SUM_LA + TPA
          SUM_LB  = SUM_LB + TMA
        ENDDO
        ATERM_SAVE(AA,N) = SUM_LA / NORM_SAVED(N,AA)
        BTERM_SAVE(AA,N) = SUM_LB / NORM_SAVED(N,AA)
        AGM(AA,N)  = ATERM_SAVE(AA,N) * GAMMA_M(AA,N)
        BGP(AA,N)  = BTERM_SAVE(AA,N) * GAMMA_P(AA,N)
      ENDDO

C  5. Green function multipliers
C     ==========================

C  For each eigenstream

      DO AA = 1, NSTREAMS
        GFUNC_DN(AA,N) = CFUNC(AA,N) * ATERM_SAVE(AA,N) * CONST
        GFUNC_UP(AA,N) = DFUNC(AA,N) * BTERM_SAVE(AA,N) * CONST
      ENDDO

C  6. Set particular integral from Green function expansion
C     =====================================================

C  particular integrals at lower and upper boundaries

      DO I = 1, NSTREAMS
        I1 = I + NSTREAMS
        S_P_U = ZERO
        S_P_L = ZERO
        S_M_U = ZERO
        S_M_L = ZERO
        DO AA = 1, NSTREAMS
          S_P_U = S_P_U + GFUNC_UP(AA,N)*XPOS(I1,AA,N)
          S_M_U = S_M_U + GFUNC_UP(AA,N)*XPOS(I,AA,N)
          S_P_L = S_P_L + GFUNC_DN(AA,N)*XPOS(I,AA,N)
          S_M_L = S_M_L + GFUNC_DN(AA,N)*XPOS(I1,AA,N)
        ENDDO
        WUPPER(I,N)  = S_P_U
        WUPPER(I1,N) = S_M_U
        WLOWER(I1,N) = S_M_L
        WLOWER(I,N)  = S_P_L
      ENDDO

C  Finish

      END

C

      SUBROUTINE LIDORT_HOM_USERSOLUTION
     I    ( GIVEN_LAYER, FOURIER )

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file of setup stuff (input to this module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  include file of solution variables (input/output to this module)

      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  include file of Multiplier coefficients (output to this module)

      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index and Fourier number (inputs)

      INTEGER          GIVEN_LAYER
      INTEGER          FOURIER

C  Local variables
C  ---------------

      INTEGER          UM, J, J1, L, N, M, AA
      DOUBLE PRECISION SUM_NEG, SUM_POS, POS1, POS2, NEG1, NEG2
      DOUBLE PRECISION RHO_P, RHO_M, SECMUI, ULP

C  Layer and Fourier

      N = GIVEN_LAYER
      M = FOURIER

C  Zeta constants (always required)

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        SECMUI = USER_SECANTS(UM)
        DO AA = 1, NSTREAMS
          RHO_P = SECMUI + KEIGEN(AA,N)
          RHO_M = SECMUI - KEIGEN(AA,N)
          ZETA_P(AA,UM,N) = ONE / RHO_P
          ZETA_M(AA,UM,N) = ONE / RHO_M
        ENDDO
      ENDDO

C  If there is no scattering, zero the user solutions and exit

      IF ( .NOT. DO_LAYER_SCATTERING(M,N) ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO AA = 1, NSTREAMS
            U_XPOS(UM,AA,N) = ZERO
            U_XNEG(UM,AA,N) = ZERO
          ENDDO
        ENDDO
        RETURN
      ENDIF

C  Eigenvector interpolation to user-defined angles
C  ------------------------------------------------

C  For each eigenvector

      DO AA = 1, NSTREAMS

C  For each moment, do inner sum over computational angles
C  for the positive and negative eigenvectors

        DO L = M, NMOMENTS
          SUM_POS = ZERO
          SUM_NEG = ZERO
          DO  J = 1, NSTREAMS
            J1 = J + NSTREAMS
            POS1 = XPOS(J1,AA,N) * WT_LEGP(J,L)
            POS2 = XPOS(J,AA,N)  * WT_LEGM(J,L)
            NEG1 = XNEG(J1,AA,N) * WT_LEGP(J,L)
            NEG2 = XNEG(J,AA,N)  * WT_LEGM(J,L)
            SUM_POS = SUM_POS + POS1 + POS2
            SUM_NEG = SUM_NEG + NEG1 + NEG2
          ENDDO
          U_HELP_P(AA,L) = SUM_POS
          U_HELP_M(AA,L) = SUM_NEG
        ENDDO

C  Now sum over all harmonic contributions at each user-defined stream

        DO UM = LOCAL_UM_START, N_USER_STREAMS
          SUM_POS = ZERO
          SUM_NEG = ZERO
          DO L = M, NMOMENTS
            ULP = U_LEG_P(UM,L) * OMEGA_MOMS(N,L)
            SUM_POS = SUM_POS + U_HELP_P(AA,L) * ULP
            SUM_NEG = SUM_NEG + U_HELP_M(AA,L) * ULP
          ENDDO
          U_XPOS(UM,AA,N) = SUM_POS
          U_XNEG(UM,AA,N) = SUM_NEG
        ENDDO

C  end eigenvector loop

      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_QBEAM_USERSOLUTION
     I    ( GIVEN_LAYER, FOURIER, IBEAM )

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file of setup stuff (input to this module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  include file of solution variables (output to this module)

      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index and Fourier number (inputs)

      INTEGER          GIVEN_LAYER
      INTEGER          FOURIER
      INTEGER          IBEAM

C  Local variables
C  ---------------

      INTEGER          UM, J, J1, L, N, M
      DOUBLE PRECISION SUM, POS1, POS2, F1
      DOUBLE PRECISION HELP1(0:MAXMOMENTS)
      DOUBLE PRECISION HELP2(0:MAXMOMENTS)

C  Layer and Fourier

      N = GIVEN_LAYER
      M = FOURIER

C  No particular solution beyond the cutoff layer
C  Or no scattering in this layer...
C  ... Zero the user solutions and exit

      IF ( N .GT. LAYER_PIS_CUTOFF(IBEAM) .OR.
     &         .NOT. DO_LAYER_SCATTERING(M,N) ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          IF ( DO_UPWELLING ) THEN
            U_WPOS1(UM,N) = ZERO
            U_WPOS2(UM,N) = ZERO
          ENDIF
          IF ( DO_DNWELLING ) THEN
            U_WNEG1(UM,N) = ZERO
            U_WNEG2(UM,N) = ZERO
          ENDIF
        ENDDO
        RETURN
      ENDIF

C  Scattering solutions
C  ====================

C  For each moment & particulate, do inner sum over computational angles
C  (required for both directions)

      F1 = FLUX_FACTOR / PI4
      DO L = M, NMOMENTS
        SUM = ZERO
        DO  J = 1, NSTREAMS
          J1 = J + NSTREAMS
          POS1 = WVEC(J1,N) * WT_LEGP(J,L)
          POS2 = WVEC(J,N)  * WT_LEGM(J,L)
          SUM = SUM + POS1 + POS2
        ENDDO
        W_HELP(L) = SUM
        HELP1(L)  = LEG0_M(L,N,IBEAM) * OMEGA_MOMS(N,L) * F1
        HELP2(L)  = W_HELP(L)         * OMEGA_MOMS(N,L)
      ENDDO

C  Now sum over all harmonic contributions at each user-defined stream
C  Distinguish between upwelling and downwelling

      IF ( DO_UPWELLING ) THEN
        IF ( STERM_LAYERMASK_UP(N) ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            POS1 = ZERO
            DO L = M, NMOMENTS
              POS1 = POS1 + HELP1(L)*U_LEG_P(UM,L)
            ENDDO
            U_WPOS1(UM,N) = POS1
            POS1 = ZERO
            DO L = M, NMOMENTS
              POS1 = POS1 + HELP2(L)*U_LEG_P(UM,L)
            ENDDO
            U_WPOS2(UM,N) = POS1
          ENDDO
        ENDIF
      ENDIF

      IF ( DO_DNWELLING ) THEN
        IF ( STERM_LAYERMASK_DN(N) ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            POS2 = ZERO
            DO L = M, NMOMENTS
              POS2 = POS2 + HELP1(L)*U_LEG_M(UM,L)
            ENDDO
            U_WNEG1(UM,N) = POS2
            POS2 = ZERO
            DO L = M, NMOMENTS
              POS2 = POS2 + HELP2(L)*U_LEG_M(UM,L)
            ENDDO
            U_WNEG2(UM,N) = POS2
          ENDDO
        ENDIF
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_GBEAM_USERSOLUTION
     I    ( GIVEN_LAYER, FOURIER, IBEAM )

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file of setup stuff (input to this module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  include file of solution variables (output to this module)

      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index and Fourier number (inputs)

      INTEGER          GIVEN_LAYER
      INTEGER          FOURIER
      INTEGER          IBEAM

C  Local variables
C  ---------------

      INTEGER          UM, L, N, M
      DOUBLE PRECISION POS1, POS2, F1

C  Layer and Fourier

      N = GIVEN_LAYER
      M = FOURIER

C  No particular solution beyond the cutoff layer
C  Or no scattering in this layer...
C  ... Zero the user solutions and exit

      IF ( N .GT. LAYER_PIS_CUTOFF(IBEAM) .OR.
     &         .NOT. DO_LAYER_SCATTERING(M,N) ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          IF ( DO_UPWELLING ) THEN
            U_WPOS1(UM,N) = ZERO
            U_WPOS2(UM,N) = ZERO
          ENDIF
          IF ( DO_DNWELLING ) THEN
            U_WNEG1(UM,N) = ZERO
            U_WNEG2(UM,N) = ZERO
          ENDIF
        ENDDO
        RETURN
      ENDIF

C  Scattering solutions
C  ====================

C  For each moment do inner sum over computational angles

      F1 = FLUX_FACTOR / PI4
      DO L = M, NMOMENTS
        W_HELP(L) = LEG0_M(L,N,IBEAM)*OMEGA_MOMS(N,L)*F1
      ENDDO

C  Now sum over all harmonic contributions at each user-defined stream

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        POS1 = ZERO
        POS2 = ZERO
        DO L = M, NMOMENTS
          POS1 = POS1 + W_HELP(L)*U_LEG_P(UM,L)
          POS2 = POS2 + W_HELP(L)*U_LEG_M(UM,L)
        ENDDO
        U_WPOS1(UM,N) = POS1
        U_WNEG1(UM,N) = POS2
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_LEGENDRE_SETUP ( FOURIER )

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file of setup stuff (Legendre output to this module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  input

      INTEGER          FOURIER

C  local variables

      DOUBLE PRECISION CFPLG(0:MAXMOMENTS), WT
      INTEGER          M, I, J, L, LPM, N, IB

C  Set integer M = Fourier number

      M = FOURIER

C  Legendre polynomials
C  --------------------

C  .. positive computational angle streams

      DO I = 1, NSTREAMS
        CALL CFPLGARR (MAXMOMENTS, NMOMENTS, M, X(I), CFPLG)
        DO L = M, NMOMENTS
          LEG_P(I,L) = CFPLG(L)
        ENDDO
      ENDDO

C  .. negative streams by symmetry relations

      DO L = M, NMOMENTS
        LPM = L + M
        IF (MOD(LPM,2).EQ.0) THEN
          DO I = 1, NSTREAMS
            LEG_M(I,L) = LEG_P(I,L)
          ENDDO
        ELSE
          DO I = 1, NSTREAMS
            LEG_M(I,L) = - LEG_P(I,L)
          ENDDO
        ENDIF
      ENDDO

C  ..solar zenith angle values

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN
       DO IB = 1, NBEAMS
        DO N = 1, NLAYERS
          CALL CFPLGARR
     &        (MAXMOMENTS, NMOMENTS, M,  SUN_SZA_COSINES(N,IB), LEG0_P)
          DO L = M, NMOMENTS
            LPM = L + M
            IF (MOD(LPM,2).EQ.0) THEN
              LEG0_M(L,N,IB) = LEG0_P(L)
            ELSE
              LEG0_M(L,N,IB) = - LEG0_P(L)
            ENDIF    
            DO I = 1, NSTREAMS
              PLMI_X0_P(I,L,N,IB) = LEG0_P(L)*LEG_P(I,L)
              PLMI_X0_M(I,L,N,IB) = LEG0_P(L)*LEG_M(I,L)
            ENDDO
          ENDDO
        ENDDO
       ENDDO
      ELSE
       DO IB = 1, NBEAMS
        CALL CFPLGARR (MAXMOMENTS, NMOMENTS, M, X0(IB), LEG0_P)
        DO L = M, NMOMENTS
          LPM = L + M
          IF (MOD(LPM,2).EQ.0) THEN
            LEG0_M(L,1,IB) = LEG0_P(L)
          ELSE
            LEG0_M(L,1,IB) = - LEG0_P(L)
          ENDIF
          DO I = 1, NSTREAMS
            PLMI_X0_P(I,L,1,IB) = LEG0_P(L)*LEG_P(I,L)
            PLMI_X0_M(I,L,1,IB) = LEG0_P(L)*LEG_M(I,L)
          ENDDO
          DO N = 2, NLAYERS
            LEG0_M(L,N,IB) = LEG0_M(L,1,IB)
            DO I = 1, NSTREAMS
              PLMI_X0_P(I,L,N,IB) = PLMI_X0_P(I,L,1,IB)
              PLMI_X0_M(I,L,N,IB) = PLMI_X0_M(I,L,1,IB)
            ENDDO
          ENDDO
        ENDDO
       ENDDO
      ENDIF

C  set up products of Associated Legendre polynomials
C  --------------------------------------------------

C  set up PLM(x).PMM(x)

      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          DO L = M, NMOMENTS
            PLMI_PLMJ_P(I,J,L) = LEG_P(I,L) * LEG_P(J,L)
            PLMI_PLMJ_M(I,J,L) = LEG_P(I,L) * LEG_M(J,L)
          ENDDO
        ENDDO
      ENDDO

C  associated products

      DO  J = 1, NSTREAMS
        WT = HALF * A(J)
        DO L = M, NMOMENTS
          WT_LEGP(J,L) = LEG_P(J,L) * WT
          WT_LEGM(J,L) = LEG_M(J,L) * WT
        ENDDO
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_USERLEGENDRE_SETUP ( FOURIER )

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file of setup stuff (Legendre output to this module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  input

      INTEGER          FOURIER

C  local variables

      DOUBLE PRECISION CFPLG(0:MAXMOMENTS)
      INTEGER          M, L, UM, LPM

C  Set integer M = Fourier number

      M = FOURIER

C  Legendre polynomials defined on user stream angles
C  --------------------------------------------------

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        CALL CFPLGARR
     &       (  MAXMOMENTS, NMOMENTS, M, USER_STREAMS(UM), CFPLG)
        DO L = M, NMOMENTS
          U_LEG_P(UM,L) = CFPLG(L)
          LPM = L + M
          IF (MOD(LPM,2).EQ.0) THEN
            U_LEG_M(UM,L) = U_LEG_P(UM,L)
          ELSE
            U_LEG_M(UM,L) = - U_LEG_P(UM,L)
          ENDIF  
        ENDDO    
      ENDDO

C  Finish

      RETURN
      END
