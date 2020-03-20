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
C #              LIDORT_L_HOM_SOLUTION                          #
C #              LIDORT_L_HOM_EIGENTRANS                        #
C #              LIDORT_L_HOM_NORMS                             #
C #              LIDORT_L_HOM_USERSOLUTION                      #
C #              LIDORT_L_QBEAM_SOLUTION                        #
C #              LIDORT_L_QBEAM_USERSOLUTION                    #
C #              LIDORT_L_GBEAM_SOLUTION                        #
C #              LIDORT_L_GBEAM_USERSOLUTION                    #
C #                                                             #
C ###############################################################

      SUBROUTINE LIDORT_L_HOM_SOLUTION
     I    ( GIVEN_LAYER, FOURIER, DOVARY, N_PARAMETERS,
     O      STATUS )

C  Linearization of the homogeneous solutions.

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup and solution stuff (inputs to module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  include file of main solution linearization variables (output)

      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'

C  include file of setup linearization variables (output)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index and Fourier number (inputs)

      INTEGER         GIVEN_LAYER
      INTEGER         FOURIER

C  Variation flag for this layer

      LOGICAL         DOVARY

C  number of varying parameters (input)

      INTEGER         N_PARAMETERS

C  output status

      INTEGER         STATUS

C  Local variables
C  ---------------

C  local matrices for eigenvalue computation

      DOUBLE PRECISION
     &      HMAT(MAXSTREAMS_P1,MAXSTREAMS_P1),
     &      HVEC(MAXSTREAMS_P1,MAX_ATMOSWFS),
     &      L_DIFVEC(MAXSTREAMS), L_EIGENVEC(MAXSTREAMS)

C  extra output/input for linear-algebra solver (LAPACK)

      INTEGER         IPIV(MAXSTREAMS_P1), INFO

C  Miscellaneous local variables

      INTEGER          I, J, I1, L, N, M, AA, K, Q, NSTRM_P1
      DOUBLE PRECISION DP, DM, SUM, L_KVAL, KSQD,
     &                 KVAL, HELP, XINV, FAC

C  Start of code
C  -------------

C  initialise status

      STATUS = LIDORT_SUCCESS

C  Layer and Fourier number

      N = GIVEN_LAYER
      M = FOURIER
      NSTRM_P1 = NSTREAMS + 1

C  nothing to do if not varying

      IF ( .NOT. DOVARY ) RETURN

C  For solution saving option with no scattering,
C  linearized solution vectors are all zero
C  Exit the routine by going to continuation point.

      IF ( DO_SOLUTION_SAVING ) THEN
        IF ( .NOT. DO_LAYER_SCATTERING(M,N) ) THEN
          DO Q = 1, N_PARAMETERS
            DO I = 1, NSTREAMS_2
              DO AA = 1, NSTREAMS
                L_XPOS(I,AA,N,Q)  = ZERO
                L_XNEG(I,AA,N,Q)  = ZERO
              ENDDO
            ENDDO
            DO AA = 1, NSTREAMS
              L_KEIGEN(AA,N,Q)  = ZERO
            ENDDO
          ENDDO
          GO TO 3456
        ENDIF
      ENDIF

C  Linearize the Eigenmatrix
C  -------------------------

C  set up linearizations of SAB and DAB

      DO I = 1, NSTREAMS
        XINV = ONE/X(I)
        DO J = 1, NSTREAMS
          FAC = HALF * A(J) * XINV
          DO Q = 1, N_PARAMETERS
            DP = ZERO
            DM = ZERO
            DO L = M, NMOMENTS
              DP = DP + PLMI_PLMJ_P(I,J,L) * L_OMEGA_MOMS(Q,N,L)
              DM = DM + PLMI_PLMJ_M(I,J,L) * L_OMEGA_MOMS(Q,N,L)
            ENDDO
            L_SAB(I,J,N,Q) = FAC * ( DP + DM )
            L_DAB(I,J,N,Q) = FAC * ( DP - DM )
          ENDDO
        ENDDO
      ENDDO

C  set up linearized eigenmatrices
C   Need to save these for all layers

      DO Q = 1, N_PARAMETERS
        DO I = 1, NSTREAMS
          DO J = 1, NSTREAMS
            SUM = ZERO
            DO K = 1, NSTREAMS
              SUM = SUM + L_DAB(I,K,N,Q) *   SAB(K,J,N)
     &                  +   DAB(I,K,N)   * L_SAB(K,J,N,Q)
            ENDDO
            L_EIGENMAT(I,J,N,Q) = SUM
          ENDDO
        ENDDO
      ENDDO

C  Matrix and column setup
C  -----------------------

C  Do this for each eigenvactor

      DO AA = 1, NSTREAMS

        KVAL = KEIGEN(AA,N)
        KSQD = KVAL * KVAL

C  initialise solution matrix HMAT (important to do this!)

        DO I = 1, MAXSTREAMS_P1
          DO J = 1, MAXSTREAMS_P1
            HMAT(I,J) = ZERO
          ENDDO
        ENDDO

C  Determine solution matrix HMAT (A in the matrix equation A.x = B)

        DO I = 1, NSTREAMS
          DO J = 1, NSTREAMS
            HMAT(I,J+1) = - EIGENMAT_SAVE(I,J,N)
          ENDDO
          HMAT(I,I+1) = HMAT(I,I+1) + KSQD
          HMAT(I,1)   = TWO * KVAL * EIGENVEC_SAVE(I,AA)
          HMAT(NSTRM_P1,I+1) = EIGENVEC_SAVE(I,AA)
        ENDDO
        HMAT(NSTRM_P1,1) = ZERO

C  solution column vectors (B in the matrix equation A.x = B).
C    (One for each parameter to be varied)

        DO Q = 1, N_PARAMETERS
          DO I = 1, NSTREAMS
            HELP = ZERO
            DO J = 1, NSTREAMS
              HELP = HELP + L_EIGENMAT(I,J,N,Q) * EIGENVEC_SAVE(J,AA)
            ENDDO
            HVEC(I,Q) = HELP
          ENDDO
          HVEC(NSTRM_P1,Q) = ZERO
        ENDDO

C  Solve matrix system for linearization factors
C  ---------------------------------------------

C  solve for the linearized sum-vector + eigenvalue
C  July 1 1999, test using LAPACK modules DGETRF and DEGTRS successful

C   .. LU-decomposition of the matrix HMAT using DGETRF

        CALL DGETRF (NSTRM_P1,NSTRM_P1,HMAT,MAXSTREAMS_P1,IPIV,INFO)

        IF ( INFO .NE. 0 ) THEN
          CALL LIDORT_LAPACK_ERROR
     I      ( INFO, N, ' Homogeneous Linearization DGETRF', STATUS )
          IF ( STATUS .EQ. LIDORT_SERIOUS ) RETURN
        ENDIF

C   .. Back substitution for column vectors (one for each parameter varying)

        CALL DGETRS ('N',NSTRM_P1,N_PARAMETERS,HMAT,
     &                   MAXSTREAMS_P1,IPIV,HVEC,MAXSTREAMS_P1,INFO)

        IF ( INFO .NE. 0 ) THEN
          CALL LIDORT_LAPACK_ERROR
     I      ( INFO, N, ' Homogeneous Linearization DGETRS', STATUS )
          IF ( STATUS .EQ. LIDORT_SERIOUS ) RETURN
        ENDIF

C  Assign linearization factors for each scatterer
C  -----------------------------------------------

C  Start loop over varying parameters

        DO Q = 1, N_PARAMETERS

C  assign linearization for eigenvalue K

          L_KVAL = HVEC(1,Q)
          L_KEIGEN(AA,N,Q) = L_KVAL

C  linearization of the actual eigenvector

          DO I = 1, NSTREAMS
            L_EIGENVEC(I) = HVEC(I+1,Q)
          ENDDO

C  linearized difference vector

          DO I = 1, NSTREAMS
            SUM = ZERO
            DO K = 1, NSTREAMS
              SUM = SUM - L_SAB(I,K,N,Q) *   EIGENVEC_SAVE(K,AA)
     &                  -   SAB(I,K,N)   * L_EIGENVEC(K)
            ENDDO
            HELP = ( SUM - L_KVAL * DIFVEC_SAVE(I,AA) ) / KVAL
            L_DIFVEC(I) = HELP
          ENDDO

C  assign linearization for 'positive' homogeneous solution vectors

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            L_XPOS(I1,AA,N,Q) = HALF*(L_EIGENVEC(I)-L_DIFVEC(I))
            L_XPOS(I,AA,N,Q)  = HALF*(L_EIGENVEC(I)+L_DIFVEC(I))
          ENDDO

C  symmetry for linearized 'negative' homogeneous solution vectors

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            L_XNEG(I1,AA,N,Q) = L_XPOS(I,AA,N,Q)
            L_XNEG(I,AA,N,Q)  = L_XPOS(I1,AA,N,Q)
          ENDDO

C  End parameter loop

        ENDDO

C  debug
C   --linearization check
c        if ( do_debug_write ) THEN
c          write(*,*)m,aa,L_KEIGEN(aa,N,2),(L_XPOS(I,AA,N,2),I=1,3)
c        endif

C  End eigenvalue loop

      ENDDO

C  Continuation point

 3456 CONTINUE

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_L_HOM_EIGENTRANS
     I    ( FOURIER )

C  Linearization of the Eigenfunction transmittances.

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup and solution stuff (inputs to module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  include file of linearization control variables (input)

      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'

C  include file of linearization solution variables (input)

      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'

C  include file of setup linearization variables (output)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'

C  subroutine arguments
C  --------------------

C  Fourier number (input)

      INTEGER         FOURIER

C  Local variables
C  ---------------

C  Miscellaneous local variables

      INTEGER          N, M, AA, Q, UT, UTA
      DOUBLE PRECISION TAU_UP, H_UP, TAU_DN, H_DN, VAR, LTQ, TBAS

C  Fourier number

      M = FOURIER

C  Linearized layer transmittances for the eigenvalues
C  ===================================================

C  start the layer loop

      DO N = 1, NLAYERS

C  When the solution saving option is set, then if there is no
C  scattering in this layer, then linearized transmittances are
C  linearized discrete ordinate transmittances.

        IF ( DO_SOLUTION_SAVING .AND.
     &           .NOT.DO_LAYER_SCATTERING(M,N) ) THEN

          IF ( LAYER_VARY_FLAG(N) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              DO AA = 1, NSTREAMS
                L_T_DELT_EIGEN(AA,N,Q) = L_T_DELT_DISORDS(AA,N,Q)
              ENDDO
            ENDDO
          ENDIF

C  Otherwise, get linearized Eigenstream transmittance factors

        ELSE

          IF ( LAYER_VARY_FLAG(N) ) THEN
            DO AA = 1, NSTREAMS
             TBAS = - T_DELT_EIGEN(AA,N) * DELTAU_VERT(N)
             DO Q = 1, LAYER_VARY_NUMBER(N)
              LTQ = L_KEIGEN(AA,N,Q) + KEIGEN(AA,N)*L_DELTAU_VERT(Q,N)
              L_T_DELT_EIGEN(AA,N,Q) = TBAS * LTQ
             ENDDO
            ENDDO
          ENDIF

        ENDIF

C  debug
C   Linearization check
c       IF ( DO_DEBUG_WRITE ) THEN
c         WRITE(*,'(3i4,1p6e16.8)')m,N,1,(L_T_DELT_EIGEN(AA,N,1),AA=1,6)
c         WRITE(*,'(3i4,1p6e16.8)')m,N,2,(L_T_DELT_EIGEN(AA,N,2),AA=1,6)
c       ENDIF

C  end layer loop

      ENDDO

C  Eigenstream transmittance factors for partial layers
C  ----------------------------------------------------

C  Code completed by R. Spurr, RTSOLUTIONS inc., 30 August 2005

      DO UTA = 1, N_OUT_USERTAUS
        IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN

          UT = OFFGRID_UTAU_OUTINDEX(UTA)
          N  = OFFGRID_UTAU_LAYERIDX(UT)

C  When the solution saving option is set, then if there is no
C  scattering in this layer, then transmittances are just the
C  discrete ordinate transmittances.

          IF ( DO_SOLUTION_SAVING .AND.
     &         .NOT.DO_LAYER_SCATTERING(M,N) ) THEN

            IF ( LAYER_VARY_FLAG(N) ) THEN
             DO AA = 1, NSTREAMS
              DO Q = 1, LAYER_VARY_NUMBER(N)
               L_T_UTDN_EIGEN(AA,UT,Q)=L_T_DISORDS_UTDN(AA,UT,Q)
               L_T_UTUP_EIGEN(AA,UT,Q)=L_T_DISORDS_UTUP(AA,UT,Q)
              ENDDO
             ENDDO
            ENDIF

C  Otherwise, Compute the Eigenstream transmittance factors

          ELSE

            IF ( LAYER_VARY_FLAG(N) ) THEN
             TAU_DN = OFFGRID_UTAU_VALUES(UT)
             TAU_UP = DELTAU_VERT(N) - TAU_DN
             DO AA = 1, NSTREAMS
              H_DN = TAU_DN * T_UTDN_EIGEN(AA,UT)
              H_UP = TAU_UP * T_UTUP_EIGEN(AA,UT)
              DO Q = 1, LAYER_VARY_NUMBER(N)
               VAR = - L_KEIGEN(AA,N,Q) - 
     &              KEIGEN(AA,N) * L_DELTAU_VERT(Q,N)
               L_T_UTDN_EIGEN(AA,UT,Q) = H_DN * VAR
               L_T_UTUP_EIGEN(AA,UT,Q) = H_UP * VAR
              ENDDO
             ENDDO
            ENDIF

          ENDIF

C  Finish off-grid optical depths loop

        ENDIF
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_L_HOM_NORMS

C  Eigenproblem, linearized solution norms for Green's function

C  Include files
C  =============

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file of main solution variables (input)

      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  include file of linearized solution variables (output)

      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'

C  Local variables
C  ---------------

C  Miscellaneous local variables

      INTEGER          N, AA, J, J1, Q
      DOUBLE PRECISION T1, T2, NORM
    
C  For all layers, save the norms

      DO N = 1, NLAYERS
        DO AA = 1, NSTREAMS
          DO Q = 1, LAYER_VARY_NUMBER(N)
            NORM = ZERO
            DO J = 1, NSTREAMS
              J1 = J + NSTREAMS
              T1 = XPOS(J,AA,N)  * L_XPOS(J,AA,N,Q)
              T2 = XPOS(J1,AA,N) * L_XPOS(J1,AA,N,Q)
              NORM = NORM + AX(J)*(T1-T2)
            ENDDO
            L_NORM_SAVED(AA,N,Q) = TWO * NORM
          ENDDO
        ENDDO
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_L_QBEAM_SOLUTION
     I    ( GIVEN_LAYER, FOURIER, IBEAM, DOVARY, N_PARAMETERS,
     O      STATUS )

C  linearized values of the classical particular solution.

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  Include file of linearization input variables

      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'

C  include files of setup and solution stuff (inputs to module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  include file of main solution linearization variables (output)

      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'

C  include file of setup linearization variables (output)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'

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

      DOUBLE PRECISION QSAVE(MAXSTREAMS)

C  linearization arrays

      DOUBLE PRECISION L_QVEC(MAXSTREAMS,MAX_ATMOSWFS)
      DOUBLE PRECISION L_QSUMVEC(MAXSTREAMS,MAX_ATMOSWFS)
      DOUBLE PRECISION L_QDIFVEC(MAXSTREAMS,MAX_ATMOSWFS)

C  help variables

      INTEGER          I, J, I1, L, N, M, Q, IB, INFO
      INTEGER          LVARY, K, K_PARAMETERS
      DOUBLE PRECISION HELP, WSUM, WDIF, SECBAR, F1
      DOUBLE PRECISION TP, TM, TM1, TM2, TM3, XINV
      LOGICAL          DO_FIRST

C  Start of code
C  -------------

C  Set layer, fourier, beam index

      N  = GIVEN_LAYER
      M  = FOURIER
      IB = IBEAM
      F1 = FLUX_FACTOR / PI4
      LVARY = 0

C  set up driving vector (using saved results)
C  This must be done regardless of whether layer N is varying or not.

      DO I = 1, NSTREAMS
        QSAVE(I) = QDIFVEC_SAVE(I) +
     &                 TWO * AVERAGE_SECANT(N,IB) * QVEC_SAVE(I)
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

C  solar zenith cosine for this layer

      SECBAR = AVERAGE_SECANT(N,IB)

C  Skip this section if there is nothing varying
C    [ Zero the boundary layer values and start Part 2 ]

      IF ( .NOT. DOVARY .OR. .NOT. DO_FIRST ) THEN
        DO I = 1, NSTREAMS_2
          DO Q = 1, N_PARAMETERS
            L_WVEC(I,N,LVARY,Q) = ZERO
          ENDDO
        ENDDO
        GO TO 2222
      ENDIF

C  For each varying parameter

      DO Q = 1, N_PARAMETERS

C  Set up linearized sum and difference vectors for Beam source terms

        DO I = 1, NSTREAMS
          XINV = ONE / X(I)
          TP = ZERO
          TM = ZERO
          DO L = M, NMOMENTS
            TP = TP + PLMI_X0_P(I,L,N,IB) * L_OMEGA_MOMS(Q,N,L)
            TM = TM + PLMI_X0_M(I,L,N,IB) * L_OMEGA_MOMS(Q,N,L)
          ENDDO
          L_QSUMVEC(I,Q) =  F1 * ( TP + TM ) * XINV
          L_QDIFVEC(I,Q) =  F1 * ( TP - TM ) * XINV
        ENDDO

C   setup linearized RHS vector
C  ( use results from the original solution )

        DO I = 1, NSTREAMS
          HELP = ZERO
          DO J = 1, NSTREAMS
            TM1 = L_EIGENMAT(I,J,N,Q) * QVEC_SAVE(J)
            TM2 =   DAB(I,J,N)   * L_QSUMVEC(J,Q)
            TM3 = L_DAB(I,J,N,Q) *   QSUMVEC_SAVE(J)
            HELP = HELP - TM1 - TM2 - TM3
          ENDDO
          L_QVEC(I,Q)  = HELP + L_QDIFVEC(I,Q) * SECBAR
        ENDDO

      ENDDO

C  additional terms for the quasi-spherical case
C  ( layers greater than one )

      IF ( .NOT. DO_PLANE_PARALLEL ) THEN
        IF ( N.GT.1 ) THEN
          DO Q = 1, N_PARAMETERS
            DO I = 1, NSTREAMS
              HELP = QSAVE(I) * L_AVERAGE_SECANT(N,LVARY,IB,Q)
              L_QVEC(I,Q) = L_QVEC(I,Q) + HELP
            ENDDO
          ENDDO
        ENDIF
      ENDIF

C  Solve problem by back substitution for all of the RHS vectors
C  ( uses L_U decomposition of QMAT_SAVE from original solution )

      CALL DGETRS
     &        ('N',NSTREAMS,N_PARAMETERS,QMAT_SAVE,
     &          MAXSTREAMS,QPIVOT,L_QVEC,MAXSTREAMS,INFO)

      IF ( INFO .NE. 0 ) THEN
        CALL LIDORT_LAPACK_ERROR
     I ( INFO, N, ' Beam P.I. linearization (N,N) (DGETRS)', STATUS )
        IF ( STATUS .EQ. LIDORT_SERIOUS ) RETURN
      ENDIF

C  assign solutions for the quasi-spherical case, N > 1

      IF ( .NOT. DO_PLANE_PARALLEL .AND. N.GT.1 ) THEN

        DO I = 1, NSTREAMS
          TM3 = - QDIF_SAVE(I) / SECBAR
          DO Q = 1, N_PARAMETERS
            I1 = I + NSTREAMS
            HELP = ZERO
            DO J = 1, NSTREAMS
              TM1 =   SAB(I,J,N)   * L_QVEC(J,Q)
              TM2 = L_SAB(I,J,N,Q) *   QVEC_SAVE(J)
              HELP = HELP - TM1 - TM2
            ENDDO
            WSUM = L_QVEC(I,Q)
            TM2 = ( HELP - L_QSUMVEC(I,Q) ) / SECBAR
            WDIF = L_AVERAGE_SECANT(N,LVARY,IB,Q) * TM3 + TM2
            L_WVEC(I,N,LVARY,Q)  = HALF * ( WSUM + WDIF )
            L_WVEC(I1,N,LVARY,Q) = HALF * ( WSUM - WDIF )
          ENDDO
        ENDDO

C  assign solutions for plane/parallel & quasi-spherical case N = 1

      ELSE

        DO I = 1, NSTREAMS
          DO Q = 1, N_PARAMETERS
            I1 = I + NSTREAMS
            HELP = ZERO
            DO J = 1, NSTREAMS
              TM1 =   SAB(I,J,N)   * L_QVEC(J,Q)
              TM2 = L_SAB(I,J,N,Q) *   QVEC_SAVE(J)
              HELP = HELP - TM1 - TM2
            ENDDO
            WSUM = L_QVEC(I,Q)
            WDIF = ( HELP - L_QSUMVEC(I,Q) ) / SECBAR
            L_WVEC(I,N,LVARY,Q)  = HALF * ( WSUM + WDIF )
            L_WVEC(I1,N,LVARY,Q) = HALF * ( WSUM - WDIF )
          ENDDO
        ENDDO

      ENDIF

C  Part 2.
C  =======

C  Continuation point

 2222 CONTINUE

C  Only for the pseudo-spherical case, profile weighting functions
C  Also not required for the column weighting functions

      IF ( DO_PLANE_PARALLEL )       RETURN
      IF ( DO_COLUMN_LINEARIZATION ) RETURN

C  No particular solution beyond the cutoff layer.
C  No solution if layer is inactive
C    [ Zero the boundary layer values and exit ]

      IF ( .NOT. DO_FIRST ) THEN
        DO K = 1, N - 1
          K_PARAMETERS = LAYER_VARY_NUMBER(K)
          DO I = 1, NSTREAMS_2
            DO Q = 1, K_PARAMETERS
              L_WVEC(I,N,K,Q) = ZERO
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
            DO Q = 1, K_PARAMETERS
              L_QVEC(I,Q) = QSAVE(I) * L_AVERAGE_SECANT(N,K,IB,Q)
            ENDDO
          ENDDO

C  Solve problem by back substitution for all of the RHS vectors
C  ( uses L_U decomposition of QMAT_SAVE from original solution )

          CALL DGETRS
     &        ('N',NSTREAMS,K_PARAMETERS,QMAT_SAVE,
     &          MAXSTREAMS,QPIVOT,L_QVEC,MAXSTREAMS,INFO)

          IF ( INFO .NE. 0 ) THEN
            CALL LIDORT_LAPACK_ERROR
     I ( INFO, N, ' Beam P.I. linearization (N,K) (DGETRS)', STATUS )
            IF ( STATUS .EQ. LIDORT_SERIOUS ) RETURN
          ENDIF

C  assign linearized solutions for layer N due to variations in layer K

          DO I = 1, NSTREAMS
            TM3 = - QDIF_SAVE(I) / SECBAR
            DO Q = 1, K_PARAMETERS
              I1 = I + NSTREAMS
              HELP = ZERO
              DO J = 1, NSTREAMS
                HELP = HELP - SAB(I,J,N) * L_QVEC(J,Q)
              ENDDO
              WSUM = L_QVEC(I,Q)
              WDIF = L_AVERAGE_SECANT(N,K,IB,Q) * TM3 + (HELP/SECBAR)
              L_WVEC(I,N,K,Q)  = HALF * ( WSUM + WDIF )
              L_WVEC(I1,N,K,Q) = HALF * ( WSUM - WDIF )
            ENDDO
          ENDDO

C  end K-layer loop

        ENDIF
      ENDDO

C  debug
c        if (m.eq.0)then
c       write(*,'(2i4,1p3e24.12)')n,1,(L_WVEC(2,N,K,1),k=1,N)
c       write(*,'(2i4,1p3e24.12)')n,2,(L_WVEC(2,N,K,2),k=1,N)
c        endif
c        if (n.eq.3)pause

C  debug
c      k = 96
c      IF ( N.GT.10.and.m.eq.0) THEN
c      DO I = 1, NSTREAMS*2
c        write(K,'(4i4,1p4e17.9)')IBEAM,M,N,I,L_WVEC(I,N,12,1),
c     &      L_INITIAL_TRANS(N,12,IBEAM,1), L_T_DELT_MUBAR(N,12,IBEAM,1)
c      ENDDO
c      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_L_GBEAM_SOLUTION
     I    ( GIVEN_LAYER, FOURIER, IBEAM, DOVARY, N_PARAMETERS )

C  Green's function beam particular integral, one layer only.

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup and solution stuff (inputs to module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  include file of main solution linearization variables (output)

      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'

C  include file of setup linearization variables (output)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index and Fourier number, beam index (inputs)

      INTEGER          GIVEN_LAYER
      INTEGER          FOURIER
      INTEGER          IBEAM

C  Variation flag

      LOGICAL          DOVARY

C  number of varying parameters (input)

      INTEGER          N_PARAMETERS

C  Local variables
C  ---------------

C  help variables

      LOGICAL          DO_FIRST
      INTEGER          AA, J, J1, L, I, I1, M, N, Q, IB
      DOUBLE PRECISION NORM, T1, T2, F1
      DOUBLE PRECISION L_SUM_LA, L_SUM_LB, L_NORM, L_ATERM, L_BTERM
      DOUBLE PRECISION TA1, TA2, TB1, TB2, L_TPA, L_TMA

C  linearizations of component variables

      DOUBLE PRECISION
     &     L_DMI(MAXSTREAMS,MAX_ATMOSWFS),
     &     L_DPI(MAXSTREAMS,MAX_ATMOSWFS)

C  initialise indices

      N  = GIVEN_LAYER
      M  = FOURIER
      IB = IBEAM
      F1 = FLUX_FACTOR / PI4

C  Check existence
C  ===============

C  This section added by R. Spurr, RTSOLUTIONS Inc., 3/4/06.

C  Only a solution if the layer is active and not below Cutoff.

      DO_FIRST = ( N .LE. LAYER_PIS_CUTOFF(IB) ) .AND.
     &                DO_LAYER_SCATTERING(M,N)

C  If there is nothing varying or if there is no solution then
C        Zero the output values and exit 

      IF ( .NOT. DOVARY .OR. .NOT. DO_FIRST ) THEN
        DO AA = 1, NSTREAMS
          DO Q = 1, N_PARAMETERS
            L_ATERM_SAVE(AA,N,Q) = ZERO
            L_BTERM_SAVE(AA,N,Q) = ZERO
          ENDDO
        ENDDO
        RETURN
      ENDIF

C  Form quantities independent of optical depth
C  ============================================

C  set up linearizations of help arrays (independent of eigenvector)

      DO Q = 1, N_PARAMETERS
        DO I = 1, NSTREAMS
          L_TPA = ZERO
          L_TMA = ZERO
          DO L = M, NMOMENTS
            L_TPA = L_TPA + PLMI_X0_P(I,L,N,IB) * L_OMEGA_MOMS(Q,N,L)
            L_TMA = L_TMA + PLMI_X0_M(I,L,N,IB) * L_OMEGA_MOMS(Q,N,L)
          ENDDO
          L_DPI(I,Q) = L_TPA * F1
          L_DMI(I,Q) = L_TMA * F1
        ENDDO
      ENDDO

C  Set up linearized Norm
C    Must be done for each Beam angle (Bug, 16 August 2005)

      DO Q = 1, N_PARAMETERS
        DO AA = 1, NSTREAMS
          NORM = ZERO
          DO J = 1, NSTREAMS
            J1 = J + NSTREAMS
            T1 = XPOS(J,AA,N)  * L_XPOS(J,AA,N,Q)
            T2 = XPOS(J1,AA,N) * L_XPOS(J1,AA,N,Q)
            NORM = NORM + AX(J)*(T1-T2)
          ENDDO
          L_NORM_SAVED(AA,N,Q) = TWO * NORM
        ENDDO
      ENDDO

C  linearize quantities independent of TAU (L_ATERM_SAVE, L_BTERM_SAVE)

      DO AA = 1, NSTREAMS
        DO Q = 1, N_PARAMETERS
          L_SUM_LA = ZERO
          L_SUM_LB = ZERO
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            TA1 = L_DMI(I,Q)*XPOS(I1,AA,N) + DMI(I)*L_XPOS(I1,AA,N,Q)
            TA2 = L_DPI(I,Q)*XPOS(I,AA,N)  + DPI(I)*L_XPOS(I,AA,N,Q)
            L_SUM_LA  = L_SUM_LA + A(I) * ( TA1 + TA2 )
            TB1 = L_DMI(I,Q)*XPOS(I,AA,N)  + DMI(I)*L_XPOS(I,AA,N,Q)
            TB2 = L_DPI(I,Q)*XPOS(I1,AA,N) + DPI(I)*L_XPOS(I1,AA,N,Q)
            L_SUM_LB  = L_SUM_LB + A(I) * ( TB1 + TB2 )
          ENDDO
          L_NORM = L_NORM_SAVED(AA,N,Q)
          L_ATERM = ( L_SUM_LA / ATERM_SAVE(AA,N) ) - L_NORM 
          L_BTERM = ( L_SUM_LB / BTERM_SAVE(AA,N) ) - L_NORM 
          L_ATERM_SAVE(AA,N,Q) = L_ATERM / NORM_SAVED(N,AA)
          L_BTERM_SAVE(AA,N,Q) = L_BTERM / NORM_SAVED(N,AA)
        ENDDO
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_L_HOM_USERSOLUTION
     I    ( GIVEN_LAYER, FOURIER, DOVARY, N_PARAMETERS )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup and solution stuff (inputs to module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  include file of main solution linearization variables (output)

      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'

C  include file of setup linearization variables (output)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index, Fourier index, parameter number (inputs)

      INTEGER          GIVEN_LAYER
      INTEGER          FOURIER
      LOGICAL          DOVARY
      INTEGER          N_PARAMETERS

C  Local variables
C  ---------------

      INTEGER          UM, J, J1, L, N, M, AA, Q
      DOUBLE PRECISION L_U_HELP_P(MAXSTREAMS,0:MAXMOMENTS)
      DOUBLE PRECISION L_U_HELP_M(MAXSTREAMS,0:MAXMOMENTS)
      DOUBLE PRECISION SUM_NEG, SUM_POS, POS1, POS2, NEG1, NEG2
      DOUBLE PRECISION ULP, L_ULP
      LOGICAL          DOLAYER_PP

C  Layer and Fourier

      N = GIVEN_LAYER
      M = FOURIER

C  general flag for existence of solutions

      DOLAYER_PP =
     &  ( ( STERM_LAYERMASK_UP(N).OR.STERM_LAYERMASK_DN(N) )
     &        .AND. DO_LAYER_SCATTERING(M,N) )       

C Zero output and return if no solutions

      IF ( .NOT.DOLAYER_PP .OR. .NOT.DOVARY ) THEN
        DO Q = 1, N_PARAMETERS
          DO AA = 1, NSTREAMS
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              L_U_XPOS(UM,AA,N,Q) = ZERO
              L_U_XNEG(UM,AA,N,Q) = ZERO
            ENDDO
          ENDDO
        ENDDO 
        RETURN
      ENDIF

C  Eigenvector interpolation to user-defined angles
C  ------------------------------------------------

C  For each eigenvector and parameter

      DO Q = 1, N_PARAMETERS
        DO AA = 1, NSTREAMS

C  For each moment, do inner sum over computational angles
C  for the positive and negative linearized eigenvectors

          DO L = M, NMOMENTS
            SUM_POS = ZERO
            SUM_NEG = ZERO
            DO  J = 1, NSTREAMS
              J1 = J + NSTREAMS
              POS1 = L_XPOS(J1,AA,N,Q) * WT_LEGP(J,L)
              POS2 = L_XPOS(J,AA,N,Q)  * WT_LEGM(J,L)
              NEG1 = L_XNEG(J1,AA,N,Q) * WT_LEGP(J,L)
              NEG2 = L_XNEG(J,AA,N,Q)  * WT_LEGM(J,L)
              SUM_POS = SUM_POS + POS1 + POS2
              SUM_NEG = SUM_NEG + NEG1 + NEG2
            ENDDO
            L_U_HELP_P(AA,L) = SUM_POS
            L_U_HELP_M(AA,L) = SUM_NEG
          ENDDO

C  Now sum over all harmonic contributions at each user-defined stream

          DO UM = LOCAL_UM_START, N_USER_STREAMS
            SUM_POS = ZERO
            SUM_NEG = ZERO
            DO L = M, NMOMENTS
              ULP   = U_LEG_P(UM,L) * OMEGA_MOMS(N,L)
              L_ULP = U_LEG_P(UM,L) * L_OMEGA_MOMS(Q,N,L)
              SUM_POS = SUM_POS + L_U_HELP_P(AA,L) *   ULP +
     &                              U_HELP_P(AA,L) * L_ULP
              SUM_NEG = SUM_NEG + L_U_HELP_M(AA,L) *   ULP +
     &                              U_HELP_M(AA,L) * L_ULP
            ENDDO
            L_U_XPOS(UM,AA,N,Q) = SUM_POS
            L_U_XNEG(UM,AA,N,Q) = SUM_NEG
          ENDDO

C  end eigenvector and parameter loop

        ENDDO
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_L_QBEAM_USERSOLUTION
     I    ( GIVEN_LAYER, FOURIER, IBEAM, DOVARY, N_PARAMETERS )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  Include file of  Linearization input variables

      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'

C  include files of setup and solution stuff (inputs to module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  include file of main solution linearization variables (output)

      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'

C  include file of setup linearization variables (output)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index, Fourier index, parameter number (inputs)

      INTEGER          GIVEN_LAYER
      INTEGER          FOURIER
      INTEGER          IBEAM

C  variation flag

      LOGICAL          DOVARY
      INTEGER          N_PARAMETERS

C  Local variables
C  ---------------

      DOUBLE PRECISION HELP(0:MAXMOMENTS), SUM, POS1, POS2
      DOUBLE PRECISION HELP_FOR(0:MAXMOMENTS), F1
      INTEGER          IB, UM, J, J1, L, N, M, Q
      INTEGER          LVARY, K, K_PARAMETERS
      LOGICAL          DO_FIRST

C  Layer and Fourier

      N = GIVEN_LAYER
      M = FOURIER
      IB = IBEAM
      LVARY = 0
      F1 = FLUX_FACTOR / PI4

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

      IF ( .NOT. DOVARY .OR. .NOT. DO_FIRST ) THEN
        DO Q = 1, N_PARAMETERS
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            L_U_WFORPOS(UM,N,Q) = ZERO
            L_U_WFORNEG(UM,N,Q) = ZERO
            L_U_WPOS(UM,N,LVARY,Q) = ZERO
            L_U_WNEG(UM,N,LVARY,Q) = ZERO
          ENDDO
        ENDDO
        GO TO 2222
      ENDIF

C  start parameter loop

      DO Q = 1, N_PARAMETERS

C  For each moment do inner sum over computational angles
C  taking care to use chain rule on the linearization

        DO L = M, NMOMENTS
          SUM = ZERO
          DO  J = 1, NSTREAMS
            J1 = J + NSTREAMS
            POS1 = L_WVEC(J1,N,LVARY,Q) * WT_LEGP(J,L)
            POS2 = L_WVEC(J,N,LVARY,Q)  * WT_LEGM(J,L)
            SUM = SUM + POS1 + POS2
          ENDDO
          HELP_FOR(L) = LEG0_M(L,N,IB) * L_OMEGA_MOMS(Q,N,L) * F1
          HELP(L)     =       SUM *   OMEGA_MOMS(N,L) +
     &                  W_HELP(L) * L_OMEGA_MOMS(Q,N,L)
        ENDDO

C  Now sum over all harmonic contributions at each user-defined stream
C  Distinguish between upwelling and downwelling

        IF ( DO_UPWELLING ) THEN
          IF ( STERM_LAYERMASK_UP(N) ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              POS1 = ZERO
              DO L = M, NMOMENTS
                POS1 = POS1 + HELP_FOR(L)*U_LEG_P(UM,L)
              ENDDO
              L_U_WFORPOS(UM,N,Q) = POS1
              POS1 = ZERO
              DO L = M, NMOMENTS
                POS1 = POS1 + HELP(L)*U_LEG_P(UM,L)
              ENDDO
              L_U_WPOS(UM,N,LVARY,Q) = POS1
            ENDDO
          ENDIF
        ENDIF

        IF ( DO_DNWELLING ) THEN
          IF ( STERM_LAYERMASK_DN(N) ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              POS2 = ZERO
              DO L = M, NMOMENTS
                POS2 = POS2 + HELP_FOR(L)*U_LEG_M(UM,L)
              ENDDO
              L_U_WFORNEG(UM,N,Q) = POS2
              POS2 = ZERO
              DO L = M, NMOMENTS
                POS2 = POS2 + HELP(L)*U_LEG_M(UM,L)
              ENDDO
              L_U_WNEG(UM,N,LVARY,Q) =  POS2
            ENDDO
          ENDIF
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

C  If no solution, zero output and exit

      IF ( .NOT. DO_FIRST ) THEN
        DO K = 1, N - 1
          K_PARAMETERS = LAYER_VARY_NUMBER(K)
          DO Q = 1, K_PARAMETERS
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              L_U_WPOS(UM,N,K,Q) = ZERO
              L_U_WNEG(UM,N,K,Q) = ZERO
            ENDDO
          ENDDO
        ENDDO
        RETURN
      ENDIF

C  start loop over all layers above N

      DO K = 1, N - 1

C  only do if layer K has some variation

        IF ( LAYER_VARY_FLAG(K) ) THEN

          K_PARAMETERS = LAYER_VARY_NUMBER(K)

C  start parameter loop

          DO Q = 1, K_PARAMETERS

C  For each moment do inner sum over computational angles
C  taking care to use chain rule on the linearization

            DO L = M, NMOMENTS
              SUM = ZERO
              DO  J = 1, NSTREAMS
                J1 = J + NSTREAMS
                POS1 = L_WVEC(J1,N,K,Q) * WT_LEGP(J,L)
                POS2 = L_WVEC(J,N,K,Q)  * WT_LEGM(J,L)
                SUM = SUM + POS1 + POS2
              ENDDO
              HELP(L) =  SUM * OMEGA_MOMS(N,L)
            ENDDO

C  Now sum over all harmonic contributions at each user-defined stream
C  Distinguish between upwelling and downwelling

            IF ( DO_UPWELLING ) THEN
              IF ( STERM_LAYERMASK_UP(N) ) THEN
                DO UM = LOCAL_UM_START, N_USER_STREAMS
                  POS1 = ZERO
                  DO L = M, NMOMENTS
                    POS1 = POS1 + HELP(L)*U_LEG_P(UM,L)
                  ENDDO
                  L_U_WPOS(UM,N,K,Q) = POS1
                ENDDO
              ENDIF
            ENDIF

            IF ( DO_DNWELLING ) THEN
              IF ( STERM_LAYERMASK_DN(N) ) THEN
                DO UM = LOCAL_UM_START, N_USER_STREAMS
                  POS2 = ZERO
                  DO L = M, NMOMENTS
                    POS2 = POS2 + HELP(L)*U_LEG_M(UM,L)
                  ENDDO
                  L_U_WNEG(UM,N,K,Q) = POS2
                ENDDO
              ENDIF
            ENDIF

C  end parameter loop

          ENDDO

C  end K-layer loop

        ENDIF
      ENDDO

C  debug
c        if (fourier.eq.1)then
c       write(*,'(i4,1p3e24.12)')n,(L_U_WPOS(1,N,K,1),K=1,N)
c       write(*,'(i4,1p3e24.12)')n,(L_U_WPOS(1,N,K,2),K=1,N)
c        if ( n.eq.3)pause
c        endif

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_L_GBEAM_USERSOLUTION
     I    ( GIVEN_LAYER, FOURIER, IBEAM, DOVARY, N_PARAMETERS )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup and solution stuff (inputs to module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  include file of main solution linearization variables (output)

      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'

C  include file of setup linearization variables (output)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index, Fourier index, parameter number (inputs)

      INTEGER          GIVEN_LAYER
      INTEGER          FOURIER
      INTEGER          IBEAM

C  Variation flag

      LOGICAL          DOVARY
      INTEGER          N_PARAMETERS

C  Local variables
C  ---------------

      LOGICAL          DO_FIRST
      INTEGER          UM, L, N, M, Q, IB
      DOUBLE PRECISION HELP(0:MAXMOMENTS), POS1, POS2, F1

C  Layer and Fourier

      N = GIVEN_LAYER
      M = FOURIER
      IB = IBEAM
      F1 = FLUX_FACTOR / PI4

C  Check existence
C  ---------------

C  Only a solution if the layer is active and not below Cutoff.

      DO_FIRST = ( N .LE. LAYER_PIS_CUTOFF(IB) ) .AND.
     &                DO_LAYER_SCATTERING(M,N)

C  If no solution or no variation, zero output and exit

      IF ( .NOT. DOVARY .OR. .NOT. DO_FIRST ) THEN
        DO Q = 1, N_PARAMETERS
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            L_U_WFORPOS(UM,N,Q) = ZERO
            L_U_WFORNEG(UM,N,Q) = ZERO
          ENDDO
        ENDDO
        RETURN
      ENDIF

C  start parameter loop

      DO Q = 1, N_PARAMETERS

C  For each moment do inner sum over computational angles

        DO L = M, NMOMENTS
          HELP(L) = LEG0_M(L,N,IB)*L_OMEGA_MOMS(Q,N,L) * F1
        ENDDO

C  Now sum over all harmonic contributions at each user-defined stream

        DO UM = LOCAL_UM_START, N_USER_STREAMS
          POS1 = ZERO
          POS2 = ZERO
          DO L = M, NMOMENTS
            POS1 = POS1 + HELP(L)*U_LEG_P(UM,L)
            POS2 = POS2 + HELP(L)*U_LEG_M(UM,L)
          ENDDO
          L_U_WFORPOS(UM,N,Q) = POS1
          L_U_WFORNEG(UM,N,Q) = POS2
        ENDDO

C  end parameter loop

      ENDDO

C  Finish

      RETURN
      END
