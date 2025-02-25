! ###########################################################
! #                                                         #
! #                    THE LIDORT FAMILY                    #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #       --         -        -        -         -          #
! #                                                         #
! ###########################################################

! ###########################################################
! #                                                         #
! #  Author :      Robert. J. D. Spurr                      #
! #                                                         #
! #  Address :     RT Solutions, Inc.                       #
! #                9 Channing Street                        #
! #                Cambridge, MA 02138, USA                 #
! #                                                         #
! #  Tel:          (617) 492 1183                           #
! #  Email :        rtsolutions@verizon.net                 #
! #                                                         #
! #  This Version :   3.7 F90                               #
! #  Release Date :   June 2014                             #
! #                                                         #
! #       NEW: THERMAL SUPPLEMENT INCLUDED    (3.2)         #
! #       NEW: OUTGOING SPHERICITY CORRECTION (3.2)         #
! #       NEW: TOTAL COLUMN JACOBIANS         (3.3)         #
! #       VLIDORT COMPATIBILITY               (3.4)         #
! #       THREADED/OPTIMIZED F90 code         (3.5)         #
! #       EXTERNAL SS / NEW I/O STRUCTURES    (3.6)         #
! #                                                         #
! #       Surface-leaving, BRDF Albedo-scaling     (3.7)    # 
! #       Taylor series, BBF Jacobians, ThreadSafe (3.7)    #
! #                                                         #
! ###########################################################

!    #####################################################
!    #                                                   #
!    #   This Version of LIDORT comes with a GNU-style   #
!    #   license. Please read the license carefully.     #
!    #                                                   #
!    #####################################################

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #     Homogeneous solution                                    #
! #                                                             #
! #              LIDORT_HOM_SOLUTION                            #
! #              LIDORT_HOM_EIGENTRANS                          #
! #              LIDORT_HOM_NORMS                               #
! #              LIDORT_HOM_USERSOLUTION                        #
! #              HMULT_MASTER (master)                          #
! #                                                             #
! #     Green's function solution                               #
! #                                                             #
! #              LIDORT_GBEAM_SOLUTION                          #
! #              LIDORT_GBEAM_USERSOLUTION                      #
! #                                                             #
! ###############################################################

!  Notes for Version 3.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob Fix  05/06/13  - Introduce TAYLOR_LIMIT parameters in module LIDORT_PARS
!     Rob  Fix 05/06/13  - Moved Zeta calculations here, Introduce limiting case scenarios (Taylor series)
!     Mick Fix 09/11/13  - Redefined ZETAs to straight sum and difference
!     Rob  Fix 10/09/13  - Small numbers analysis finalized using Taylor_Order parameter

module lidort_solutions

!  Parameter types

   USE LIDORT_PARS, only : fpk

!  Dependencies

   USE lidort_aux     , only : ASYMTX
   USE lidort_Taylor_m, only : TAYLOR_SERIES_1


!private
public

contains

SUBROUTINE LIDORT_HOM_SOLUTION                                  &
         ( DO_SOLUTION_SAVING, NSTREAMS, NMOMENTS,              & ! Inputs
           GIVEN_LAYER, FOURIER, DO_LAYER_SCATTERING,           & ! Inputs
           OMEGA_MOMS, QUAD_STREAMS, QUAD_WEIGHTS,              & ! Inputs
           PLMI_PLMJ_P, PLMI_PLMJ_M,                            & ! Inputs
           SAB, DAB, EIGENMAT_SAVE, EIGENVEC_SAVE, DIFVEC_SAVE, & ! output
           KEIGEN, XPOS, XNEG,                                  & ! output
           STATUS, MESSAGE, TRACE )                               ! output

!  Numerical solution of Eigenproblem.

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS,  &
                              ZERO, ONE, HALF,                                  &
                              LIDORT_SUCCESS, LIDORT_SERIOUS

!  Implicit none

      IMPLICIT NONE

!  subroutine arguments
!  --------------------

!  Solution saving flag

      LOGICAL  , intent(in) :: DO_SOLUTION_SAVING

!  Number of streams and moments

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: NMOMENTS

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STREAMS( MAXSTREAMS )
      REAL(fpk), intent(in)  :: QUAD_WEIGHTS( MAXSTREAMS )

!  Given layer index and Fourier number (inputs)

      INTEGER  , intent(in)  :: GIVEN_LAYER
      INTEGER  , intent(in)  :: FOURIER

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Saved array involving product of OMEGA and phase function moments

      REAL(fpk), intent(in)  :: OMEGA_MOMS ( MAXLAYERS, 0:MAXMOMENTS )

!  Legendre polynomial products

      REAL(fpk), intent(in)  :: PLMI_PLMJ_P(MAXSTREAMS,MAXSTREAMS,0:MAXMOMENTS)
      REAL(fpk), intent(in)  :: PLMI_PLMJ_M(MAXSTREAMS,MAXSTREAMS,0:MAXMOMENTS)

!  Solutions to the homogeneous RT equations
!  -----------------------------------------

!  local matrices for eigenvalue computation

      REAL(fpk), intent(out) :: SAB(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk), intent(out) :: DAB(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk), intent(out) :: EIGENMAT_SAVE(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk), intent(out) :: EIGENVEC_SAVE(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk), intent(out) :: DIFVEC_SAVE  (MAXSTREAMS,MAXSTREAMS)

!mick fix 6/29/11 - changed next 3 outputs from "out" to "inout"

!  (Positive) Eigenvalues

      REAL(fpk), intent(inout) :: KEIGEN(MAXSTREAMS,MAXLAYERS)

!  Eigenvector solutions

      REAL(fpk), intent(inout) :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(inout) :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Exception handling. Updated 18 May 2010.

      INTEGER      , intent(out) :: STATUS
      CHARACTER*(*), intent(out) :: MESSAGE, TRACE

!  Local variables
!  ---------------

!  local matrices for eigenvalue computation

      REAL(fpk)    :: EIGENMAT(MAXSTREAMS,MAXSTREAMS)

!  (output from Eigenpackage module ASYMTX)

      REAL(fpk)    :: KSQ(MAXSTREAMS), WK(MAXSTREAMS_2)
      REAL(fpk)    :: EVEC(MAXSTREAMS,MAXSTREAMS)
      INTEGER      :: IER
      LOGICAL      :: ASYMTX_FAILURE

!  Miscellaneous local variables

      INTEGER      :: I, J, I1, L, N, M, AA, K
      REAL(fpk)    :: DP, DM, SUM, FAC, KVAL, NORM, XINV
      CHARACTER*3  :: CN, CI

!  initialise exception handling

      STATUS = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Layer and Fourier

      N = GIVEN_LAYER
      M = FOURIER

!  When the solution saving option is set, then if there is no
!  scattering in this layer, then the eigenvectors are just the
!  discrete ordinates, the eigenvectors are unit vectors in the
!  discrete ordinate directions.

      IF ( DO_SOLUTION_SAVING ) THEN
        IF ( .NOT. DO_LAYER_SCATTERING(M,N) ) THEN
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            XINV = ONE / QUAD_STREAMS(I)
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
          RETURN
        ENDIF
      ENDIF

!  Scattering solutions
!  ====================

!  Construct Eigenmatrix
!  ---------------------

!  zero the Eigenmatrix

      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          EIGENMAT(I,J) = ZERO
        ENDDO
      ENDDO

!  Develop Sum and Difference matrices

      DO I = 1, NSTREAMS
        XINV = ONE/QUAD_STREAMS(I)
        DO J = 1, NSTREAMS
          FAC = XINV * HALF * QUAD_WEIGHTS(J)
          DP = ZERO
          DM = ZERO
          DO L = M, NMOMENTS
            DP = DP + PLMI_PLMJ_P(I,J,L) * OMEGA_MOMS(N,L)
            DM = DM + PLMI_PLMJ_M(I,J,L) * OMEGA_MOMS(N,L)
          ENDDO
          SAB(I,J) = FAC * ( DP + DM ) 
          DAB(I,J) = FAC * ( DP - DM )
        ENDDO
        SAB(I,I) = SAB(I,I) - XINV
        DAB(I,I) = DAB(I,I) - XINV
      ENDDO

!  Compute Eigenmatrix

      DO I = 1, NSTREAMS 
        DO J = 1, NSTREAMS
          SUM = ZERO
          DO K = 1, NSTREAMS
            SUM = SUM + DAB(I,K) * SAB(K,J)
          ENDDO
          EIGENMAT(I,J) = SUM
        ENDDO
      ENDDO

!  save Eigenmatrix (original is destroyed by ASMTYX)

      DO I = 1, NSTREAMS 
        DO J = 1, NSTREAMS
          EIGENMAT_SAVE(I,J) = EIGENMAT(I,J)
        ENDDO
      ENDDO

!  Eigensolution package
!  ---------------------

!  Let's see how we get on with the DISORT package
!  tested 19 May 1999. Gives same results as Analytical approach.

!  second test using DGEEV (LAPACK module) - Worked against ASYMTX.
!  However, DGEEV is 1.7 - 2.0 times slower because it must look for
!  complex roots. [ASYMTX was specially written to avoid this search].
!  Also first call to DGEEV is 20 times slower. Tested June 24, 1999.
!  Conclusion stick with ASYMTX for now.

      CALL  ASYMTX ( EIGENMAT, NSTREAMS, MAXSTREAMS, MAXSTREAMS, &
                     EVEC, KSQ, IER, WK, MESSAGE, ASYMTX_FAILURE )

!  Exception handling 1

      IF ( ASYMTX_FAILURE  ) THEN
        WRITE(CN,'(I3)')N
        TRACE   ='ASYMTX error in LIDORT_HOM_SOLUTION, Layer='//CN
        STATUS  = LIDORT_SERIOUS
        RETURN
      ENDIF

! Exception handling 2

      IF ( IER.GT.0 ) THEN
        WRITE(CI,'(I3)')IER
        WRITE(CN,'(I3)')N
        MESSAGE = 'eigenvalue '//CI//' has not converged'
        TRACE   = 'ASYMTX error in LIDORT_HOM_SOLUTION, Layer='//CN
        STATUS = LIDORT_SERIOUS
        RETURN
      ENDIF

!  second test using DGEEV (LAPACK module). Here is what to use.
!        CALL DGEEV   ( 'N', 'V', NSTREAMS, EIGENMAT,
!     O        MAXSTREAMS, REAL_KSQ, IMAG_KSQ,
!     O        LEFT_EVEC, MAXSTREAMS, EVEC, MAXSTREAMS,
!     W        WORK, LWORK, LAPACK_INFO )

!  Find solution eigenvectors XPOS, XNEG for all eigenvalues
!  ---------------------------------------------------------

      DO AA = 1, NSTREAMS

!  Store positive values in output array for each layer

        KVAL = DSQRT(KSQ(AA))
        KEIGEN(AA,N) = KVAL

!  Normalize eigenvectors to 1

        NORM = ZERO
        DO I = 1, NSTREAMS
          NORM = NORM + EVEC(I,AA)*EVEC(I,AA)
        ENDDO
        NORM = DSQRT(NORM)

!  Find normalized eigenvector EIGENVEC_SAVE

        DO I = 1, NSTREAMS
          EIGENVEC_SAVE(I,AA) = EVEC(I,AA)/NORM
        ENDDO

!  Find difference eigenvector DIFVE! (Siewert's notation)

        DO I = 1, NSTREAMS
          SUM = ZERO
          DO K = 1, NSTREAMS
            SUM = SUM - SAB(I,K) * EIGENVEC_SAVE(K,AA)
          ENDDO
          DIFVEC_SAVE(I,AA) = SUM / KVAL
          ENDDO

!  assign original evectors; first N are "DOWN", last N are "UP" (streams)

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          XPOS(I,AA,N)  = HALF * ( EIGENVEC_SAVE(I,AA) + DIFVEC_SAVE(I,AA) )
          XPOS(I1,AA,N) = HALF * ( EIGENVEC_SAVE(I,AA) - DIFVEC_SAVE(I,AA) )
        ENDDO

!  Use symmetry properties to set -ve eigenvectors

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          XNEG(I1,AA,N) = XPOS(I,AA,N)
          XNEG(I,AA,N)  = XPOS(I1,AA,N)
        ENDDO

!  End eigenstream loop

      ENDDO

!  Finish

      RETURN
END SUBROUTINE LIDORT_HOM_SOLUTION

!

SUBROUTINE LIDORT_HOM_EIGENTRANS                                      &
    ( DO_SOLUTION_SAVING, NSTREAMS, NLAYERS, N_USER_LEVELS,           & ! input
      PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,   & ! input
      FOURIER, DO_LAYER_SCATTERING, DELTAU_VERT, PARTAU_VERT, KEIGEN, & ! input
      T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,                 & ! Output
      T_DELT_EIGEN,   T_UTUP_EIGEN,   T_UTDN_EIGEN )                    ! Output

!  Eigenproblem, transmittance matrices

!  module,  dimensions and numbers

      USE LIDORT_pars, only : MAXMOMENTS, MAXSTREAMS, MAX_USER_LEVELS, MAXLAYERS, MAX_PARTLAYERS, &
                              MAX_TAU_QPATH, ZERO

      IMPLICIT NONE

!  subroutine arguments
!  --------------------

!  Solution saving flag

      LOGICAL  , intent(in)  :: DO_SOLUTION_SAVING

!  Number of streams and layers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: NLAYERS

!  Number of output levels

      INTEGER  , intent(in)  :: N_USER_LEVELS

!  output optical depth masks and indices
!    off-grid optical depth mask

      LOGICAL  , intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  Fourier number (input)

      INTEGER  , intent(in)  :: FOURIER

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Optical depths

      REAL(fpk), intent(in)  :: DELTAU_VERT ( MAXLAYERS )
      REAL(fpk), intent(in)  :: PARTAU_VERT ( MAX_PARTLAYERS )

!  (Positive) Eigenvalues

      REAL(fpk), intent(in)  :: KEIGEN(MAXSTREAMS,MAXLAYERS)

!  discrete ordinate factors (BVP telescoping, solutions saving)
!  Code added by R. Spurr, RT SOLUTIONS Inc., 30 August 2005.

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_DISORDS_UTUP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_DISORDS_UTDN(MAXSTREAMS,MAX_PARTLAYERS)

!  Subroutine output arguments
!  ---------------------------

!  transmittance factors for +/- eigenvalues
!     Whole layer (DELTA), User optical depths (UTUP and UTDN)
!     These depend on eigensolutions and will change for each Fourier

      REAL(fpk), intent(out) :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(out) :: T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(out) :: T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Local variables
!  ---------------

!  Miscellaneous local variables

      INTEGER      :: N, M, AA, UT, UTA
      REAL(fpk)    :: HELP, TAU_UP, TAU_DN
    
!  Fourier

      M = FOURIER

!  Layer transmittances for the eigenvalues
!  ========================================

!  start the layer loop

      DO N = 1, NLAYERS

!  When the solution saving option is set, then if there is no
!  scattering in this layer, then transmittances are just the
!  discrete ordinate transmittances.

        IF ( DO_SOLUTION_SAVING .AND. &
              .NOT.DO_LAYER_SCATTERING(M,N) ) THEN

          DO AA = 1, NSTREAMS
            T_DELT_EIGEN(AA,N) = T_DELT_DISORDS(AA,N)
          ENDDO

!  Otherwise, get the full set of Eigenstream transmittance factors

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

!  end layer loop

      ENDDO

!  Eigenstream transmittance factors for partial layers
!  ----------------------------------------------------

!  Code completed by R. Spurr, RT SOLUTIONS Inc. 30 August 2005

      DO UTA = 1, N_USER_LEVELS
       IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
        UT = PARTLAYERS_OUTINDEX(UTA)
        N  = PARTLAYERS_LAYERIDX(UT)

!  When the solution saving option is set, then if there is no
!  scattering in this layer, then transmittances are just the
!  discrete ordinate transmittances.

        IF ( DO_SOLUTION_SAVING .AND. &
              .NOT.DO_LAYER_SCATTERING(M,N) ) THEN

          DO AA = 1, NSTREAMS
            T_UTDN_EIGEN(AA,UT) = T_DISORDS_UTDN(AA,UT)
            T_UTUP_EIGEN(AA,UT) = T_DISORDS_UTUP(AA,UT)
          ENDDO

!  Otherwise, Compute the Eigenstream transmittance factors

        ELSE

          TAU_DN = PARTAU_VERT(UT)
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

!  end loop over off-grid optical depths

       ENDIF
      ENDDO

!  Finish

      RETURN
END SUBROUTINE LIDORT_HOM_EIGENTRANS

!

SUBROUTINE LIDORT_HOM_NORMS                       &
         ( NSTREAMS, NLAYERS, QUAD_STRMWTS, XPOS, & ! Input
           NORM_SAVED )                             ! Output

!  Eigenproblem, solution norms for Green's function

!  module,  dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAXSTREAMS_2, MAXLAYERS, ZERO

      IMPLICIT NONE

!  subroutine arguments
!  --------------------

!  Number of streams and layers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: NLAYERS

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STRMWTS( MAXSTREAMS )

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  outputs
!  -------

!  Saved quantities for the Green function solution

      REAL(fpk), intent(out) :: NORM_SAVED(MAXLAYERS,MAXSTREAMS)

!  Local variables
!  ---------------

!  Miscellaneous local variables

      INTEGER      :: N, AA, J, J1
      REAL(fpk)    :: T1, T2, NORM
    
!  For all layers, save the norms

      DO N = 1, NLAYERS
        DO AA = 1, NSTREAMS
          NORM = ZERO
          DO J = 1, NSTREAMS
            J1 = J + NSTREAMS
            T1 = XPOS(J,AA,N)  * XPOS(J,AA,N)
            T2 = XPOS(J1,AA,N) * XPOS(J1,AA,N)
            NORM = NORM + QUAD_STRMWTS(J)*(T1-T2)
          ENDDO
          NORM_SAVED(N,AA) = NORM
        ENDDO
      ENDDO

!  Finish

      RETURN
END SUBROUTINE LIDORT_HOM_NORMS

!

SUBROUTINE LIDORT_HOM_USERSOLUTION                              &
         ( NSTREAMS, N_USER_STREAMS, NMOMENTS,                  & ! input
           GIVEN_LAYER, FOURIER, DO_LAYER_SCATTERING,           & ! input
           XPOS, XNEG, OMEGA_MOMS, WT_LEGP, WT_LEGM, U_LEG_P,   & ! input
           U_XPOS, U_XNEG, U_HELP_P, U_HELP_M )                   ! output

!   Rob fix 5/6/13 - Zeta_M and Zeta_P calculations moved.

!  module,  dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS, MAX_USER_STREAMS, &
                              ZERO, ONE

      IMPLICIT NONE

!  subroutine input arguments
!  --------------------------

!  Number of streams and moments

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_USER_STREAMS
      INTEGER  , intent(in)  :: NMOMENTS

!  Given layer index and Fourier number (inputs)

      INTEGER  , intent(in)  :: GIVEN_LAYER
      INTEGER  , intent(in)  :: FOURIER

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!Rob fix 5/6/13 - This argument removed
!      REAL(fpk), intent(in)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Saved array involving product of OMEGA and phase function moments

      REAL(fpk), intent(in)  :: OMEGA_MOMS ( MAXLAYERS, 0:MAXMOMENTS )

!Rob fix 5/6/13 - This argument removed
!      REAL(fpk), intent(in)  :: KEIGEN(MAXSTREAMS,MAXLAYERS)

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Polynomial-weight Legendre products

      REAL(fpk), intent(in)  :: WT_LEGP(MAXSTREAMS,0:MAXMOMENTS)
      REAL(fpk), intent(in)  :: WT_LEGM(MAXSTREAMS,0:MAXMOMENTS)

!  Legendre functions on User defined polar angles

      REAL(fpk), intent(in)  :: U_LEG_P(MAX_USER_STREAMS,0:MAXMOMENTS)

!  Subroutine output arguments
!  ---------------------------

!  Saved help variables

      REAL(fpk), intent(out) :: U_HELP_P(MAXSTREAMS,0:MAXMOMENTS)
      REAL(fpk), intent(out) :: U_HELP_M(MAXSTREAMS,0:MAXMOMENTS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(fpk), intent(inout) :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(inout) :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Local variables
!  ---------------

      INTEGER       :: UM, J, J1, L, N, M, AA
      REAL(fpk)    :: SUM_NEG, SUM_POS, POS1, POS2, NEG1, NEG2, ULP

!  Layer and Fourier

      N = GIVEN_LAYER
      M = FOURIER

!  If there is no scattering, zero the user solutions and exit

      IF ( .NOT. DO_LAYER_SCATTERING(M,N) ) THEN
        DO UM = 1, N_USER_STREAMS
          DO AA = 1, NSTREAMS
            U_XPOS(UM,AA,N) = ZERO
            U_XNEG(UM,AA,N) = ZERO
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Eigenvector interpolation to user-defined angles
!  ------------------------------------------------

!  For each eigenvector

      DO AA = 1, NSTREAMS

!  For each moment, do inner sum over computational angles
!  for the positive and negative eigenvectors

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

!  Now sum over all harmonic contributions at each user-defined stream

        DO UM = 1, N_USER_STREAMS
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

!  end eigenvector loop

      ENDDO

!  Finish

      RETURN
END SUBROUTINE LIDORT_HOM_USERSOLUTION

!

SUBROUTINE HMULT_MASTER                                                &
       ( M, DO_UPWELLING, DO_DNWELLING, TAYLOR_ORDER,                  & ! Input
         NSTREAMS, N_USER_STREAMS, NLAYERS, N_USER_LEVELS,             & ! Input
         PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX, & ! Input
         USER_SECANTS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,         & ! Input
         KEIGEN, T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,             & ! Input (Rob Fix 5/6/13, New argument)
         T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM, DELTAUS, PARTAUS,   & ! Input (Rob Fix 5/6/13, New argument)
         ZETA_M, ZETA_P, HMULT_1, HMULT_2,                             & ! Output
         UT_HMULT_UU, UT_HMULT_UD, UT_HMULT_DU, UT_HMULT_DD )            ! Output

!  Homogeneous solution multipliers

!  Notes for Version 3.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob Fix  05/06/13  - Introduce TAYLOR_LIMIT parameters in module LIDORT_PARS
!     Rob  Fix 05/06/13  - Moved Zeta calculations here, Introduce limiting case scenarios (Taylor series)
!     Mick Fix 09/11/13  - Redefined ZETAs to straight sum and difference
!     Rob  Fix 10/09/13  - Small numbers analysis finalized using Taylor_Order parameter

      USE LIDORT_pars, only : MAXSTREAMS, MAX_USER_STREAMS, MAX_USER_LEVELS, &
                              MAXLAYERS, MAX_PARTLAYERS, ONE, TAYLOR_SMALL

      IMPLICIT NONE

!  subroutine arguments
!  --------------------

!  FOURIER COMPONENT (DEBUG)

      INTEGER  , intent(in)  :: M

!  Direction flags

      LOGICAL  , intent(in)  :: DO_UPWELLING
      LOGICAL  , intent(in)  :: DO_DNWELLING

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Number of streams and layers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_USER_STREAMS
      INTEGER  , intent(in)  :: NLAYERS

!  Number of output levels

      INTEGER  , intent(in)  :: N_USER_LEVELS

!  output optical depth masks and indices
!    off-grid optical depth mask

      LOGICAL  , intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX     (MAX_PARTLAYERS)

!  User stream cosines

      REAL(fpk), intent(in)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Layer masks for doing integrated source terms

      LOGICAL  , intent(in)  :: STERM_LAYERMASK_UP(MAXLAYERS)
      LOGICAL  , intent(in)  :: STERM_LAYERMASK_DN(MAXLAYERS)

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  (Positive) Eigenvalues. Argument required for Zeta calculations

      REAL(fpk), intent(in)  :: KEIGEN(MAXSTREAMS,MAXLAYERS)

!  transmittance factors for +/- eigenvalues
!     Whole layer (DELTA), User optical depths (UTUP and UTDN)
!     These depend on eigensolutions and will change for each Fourier

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Input Optical depths required for Taylor-series limiting cases

      REAL(fpk), intent(in)  :: DELTAUS(MAXLAYERS)
      REAL(fpk), intent(in)  :: PARTAUS(MAX_PARTLAYERS)

!  Output = Global multipliers
!  ===========================

!  coefficient functions for user-defined angles. Intent "inout" (basically, Zetas are now output)

      REAL(fpk), intent(inout)  :: ZETA_M(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(inout)  :: ZETA_P(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  Integrated homogeneous solution multipliers, whole layer

      REAL(fpk), intent(out) :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(out) :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  Integrated homogeneous solution multipliers, partial layer

      REAL(fpk), intent(out) :: UT_HMULT_UU(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(out) :: UT_HMULT_UD(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(out) :: UT_HMULT_DU(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(out) :: UT_HMULT_DD(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Local variables
!  ---------------

      INTEGER      :: UM, K, N, UT, UTA
      REAL(fpk)    :: UDEL, SM, ZDEL, ZUDEL,  UX_UP, UX_DN, ZX_UP, ZX_DN, DELT, EPS
      LOGICAL      :: STERM_EXIST(MAXLAYERS)

!  Existence flags

      DO N = 1, NLAYERS
        STERM_EXIST(N) = ( STERM_LAYERMASK_UP(N).OR.STERM_LAYERMASK_DN(N) )
      ENDDO

!  Whole layer multipliers
!  -----------------------

!  Start loops over layers and user-streams
!    Only done if layers are flagged

      DO N = 1, NLAYERS
        IF ( STERM_EXIST(N) ) THEN
          DO UM = 1, N_USER_STREAMS
            UDEL = T_DELT_USERM(N,UM)
            SM   = USER_SECANTS(UM)
            DO K = 1, NSTREAMS
              ZETA_P(K,UM,N) = SM + KEIGEN(K,N)
              ZETA_M(K,UM,N) = SM - KEIGEN(K,N)
              ZDEL    = T_DELT_EIGEN(K,N)
              ZUDEL   = ZDEL * UDEL
              HMULT_2(K,UM,N) = SM * ( ONE - ZUDEL ) / ZETA_P(K,UM,N)
              IF ( ABS(ZETA_M(K,UM,N)) .LT. TAYLOR_SMALL ) THEN
                 EPS = ZETA_M(K,UM,N)
                 CALL TAYLOR_SERIES_1 ( TAYLOR_ORDER, EPS, DELTAUS(N), UDEL, SM, HMULT_1(K,UM,N) )
              ELSE
                 HMULT_1(K,UM,N) = SM * ( ZDEL - UDEL ) / ZETA_M(K,UM,N)
              ENDIF
!              if (m.eq.0.and.n.eq.23.and.k.eq.6)write(67,*)um,n,k,ABS(ZETA_M(K,UM,N)),&
!                         HMULT_1(K,UM,N),SM * ( ZDEL - UDEL ) / ZETA_M(K,UM,N),&
!                         1.0d0-(HMULT_1(K,UM,N)/(SM * ( ZDEL - UDEL ) / ZETA_M(K,UM,N)))
            ENDDO
          ENDDO
        ENDIF
      ENDDO

!  Partial layer multipliers
!  -------------------------

      DO UTA = 1, N_USER_LEVELS
        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)
          IF ( STERM_EXIST(N) ) THEN

!  Upwelling

            IF ( DO_UPWELLING ) THEN
              DO UM = 1, N_USER_STREAMS
                UX_UP = T_UTUP_USERM(UT,UM)
                SM    = USER_SECANTS(UM)
                DO K = 1, NSTREAMS
                  ZDEL     = T_DELT_EIGEN(K,N)
                  ZX_UP    = T_UTUP_EIGEN(K,UT)
                  ZX_DN    = T_UTDN_EIGEN(K,UT)
                  IF ( ABS(ZETA_M(K,UM,N)) .LT. TAYLOR_SMALL ) THEN
                     EPS = ZETA_M(K,UM,N) ; DELT = DELTAUS(N) - PARTAUS(UT)
                     CALL TAYLOR_SERIES_1 ( TAYLOR_ORDER, EPS, DELT, UX_UP, SM, UT_HMULT_UU(K,UM,UT) )
                  ELSE
                     UT_HMULT_UU(K,UM,UT) = SM * ( ZX_UP - UX_UP ) / ZETA_M(K,UM,N)
                  ENDIF
                  UT_HMULT_UD(K,UM,UT) = SM * ( ZX_DN - ZDEL * UX_UP)  / ZETA_P(K,UM,N)
                ENDDO
              ENDDO
            ENDIF

!  Downwelling

            IF ( DO_DNWELLING ) THEN
              DO UM = 1, N_USER_STREAMS
                UX_DN = T_UTDN_USERM(UT,UM)
                SM    = USER_SECANTS(UM)
                DO K = 1, NSTREAMS
                  ZDEL     = T_DELT_EIGEN(K,N)
                  ZX_UP    = T_UTUP_EIGEN(K,UT)
                  ZX_DN    = T_UTDN_EIGEN(K,UT)
                  IF ( ABS(ZETA_M(K,UM,N)) .LT. TAYLOR_SMALL ) THEN
                     EPS = ZETA_M(K,UM,N) ; DELT = PARTAUS(UT)
                     CALL TAYLOR_SERIES_1 ( TAYLOR_ORDER, EPS, DELT, UX_DN, SM, UT_HMULT_DD(K,UM,UT) )
                  ELSE
                     UT_HMULT_DD(K,UM,UT) = SM * ( ZX_DN - UX_DN ) / ZETA_M(K,UM,N)
                  ENDIF
!               if (m.eq.0.and.n.eq.23.and.k.eq.6)write(69,*)um,ut,k,ABS(ZETA_M(K,UM,N)),&
!                         UT_HMULT_DD(K,UM,UT),SM * ( ZX_DN - UX_DN ) / ZETA_M(K,UM,N),&
!                         1.0d0-(UT_HMULT_DD(K,UM,UT)/(SM * ( ZX_DN - UX_DN ) / ZETA_M(K,UM,N)))
                 UT_HMULT_DU(K,UM,UT) = SM *  ( ZX_UP - ZDEL * UX_DN ) / ZETA_P(K,UM,N)
                ENDDO
              ENDDO
            ENDIF

!  end loop over partial layers

          ENDIF
        ENDIF
      ENDDO

!  Finish

      RETURN
END SUBROUTINE HMULT_MASTER

!

SUBROUTINE LIDORT_GBEAM_SOLUTION ( TAYLOR_ORDER,                      & ! Input
           NSTREAMS, NSTREAMS_2, NMOMENTS, GIVEN_LAYER, FOURIER,      & ! input
           FLUX_FACTOR, IBEAM, DO_LAYER_SCATTERING, LAYER_PIS_CUTOFF, & ! input
           INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR, DELTAUS,      & ! input
           QUAD_WTS, OMEGA_MOMS, KEIGEN, XPOS, T_DELT_EIGEN,          & ! input
           NORM_SAVED, PLMI_X0_P, PLMI_X0_M,                          & ! input
           GAMMA_M, GAMMA_P, DMI, DPI, ATERM_SAVE, BTERM_SAVE,        & ! Output
           CFUNC, DFUNC, AGM, BGP, GFUNC_UP, GFUNC_DN,                & ! Output
           WUPPER, WLOWER )                                             ! Output

!  Green's function beam particular integral, one layer only.
!  Uses coefficient expansion of attenuation.

!  Notes for Version 3.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob Fix  05/06/13  - Introduce TAYLOR_LIMIT parameters in module LIDORT_PARS
!     Rob  Fix 05/06/13  - Move Gamma calculations here, Introduce limiting case scenarios (Taylor series)
!     Mick Fix 09/11/13  - Redefined GAMMAs to straight sum and difference
!     Rob  Fix 10/09/13  - Small numbers analysis finalized using Taylor_Order parameter

      USE LIDORT_pars, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, &
                              MAXLAYERS, MAXBEAMS,                  &
                              ZERO, ONE, PI4, TAYLOR_SMALL

      IMPLICIT NONE

!  subroutine input arguments
!  ==========================

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Number of streams and moments

      INTEGER  , intent(in)  :: NSTREAMS, NSTREAMS_2
      INTEGER  , intent(in)  :: NMOMENTS

!  Given layer index and Fourier number (inputs)

      INTEGER  , intent(in)  :: GIVEN_LAYER
      INTEGER  , intent(in)  :: FOURIER

!  Solar Flux

      REAL(fpk), intent(in)  :: FLUX_FACTOR

!  Beam index

      INTEGER  , intent(in)  :: IBEAM

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_WTS( MAXSTREAMS )

!  Saved array involving product of OMEGA and phase function moments

      REAL(fpk), intent(in)  :: OMEGA_MOMS ( MAXLAYERS, 0:MAXMOMENTS )

!  initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Transmittance factors for average secant stream

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  Average-secant for solar beams

      REAL(fpk), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )

!Rob Fix 5/6/13 - Input Optical depths required for Taylor-series limiting cases

      REAL(fpk), intent(in)  :: DELTAUS(MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  Legendre polynomial products

      REAL(fpk), intent(in)  :: PLMI_X0_P(MAXSTREAMS,0:MAXMOMENTS,MAXLAYERS,MAXBEAMS)
      REAL(fpk), intent(in)  :: PLMI_X0_M(MAXSTREAMS,0:MAXMOMENTS,MAXLAYERS,MAXBEAMS)

!  (Positive) Eigenvalues

      REAL(fpk), intent(in)  :: KEIGEN(MAXSTREAMS,MAXLAYERS)

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: NORM_SAVED(MAXLAYERS,MAXSTREAMS)

!  subroutine output arguments
!  ===========================

!  Green function solution
!  ***********************

!mick fix 6/29/11 - change most outputs from "out" to "inout"

!  Saved quantities for the Green function solution

      REAL(fpk), intent(inout) :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(inout) :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

      REAL(fpk), intent(out) :: DMI(MAXSTREAMS), DPI(MAXSTREAMS)

      REAL(fpk), intent(inout) :: AGM(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(inout) :: BGP(MAXSTREAMS,MAXLAYERS)

!  Layer C and D functions

      REAL(fpk), intent(inout) :: CFUNC(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(inout) :: DFUNC(MAXSTREAMS,MAXLAYERS)

!  Green function Multipliers for solution
!         ( GFUNC_DN = CFUNC * ATERM_SAVE )
!         ( GFUNC_UP = DFUNC * BTERM_SAVE )

      REAL(fpk), intent(inout) :: GFUNC_UP(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(inout) :: GFUNC_DN(MAXSTREAMS,MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(inout) :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(inout) :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  General beam solutions at the boundaries

      REAL(fpk), intent(inout) :: WUPPER(MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(inout) :: WLOWER(MAXSTREAMS_2,MAXLAYERS)

!  local variables
!  ===============

      INTEGER      :: AA, L, I, I1, M, N
      REAL(fpk)    :: TPA, TMA, SUM_LA, SUM_LB
      REAL(fpk)    :: S_P_U, S_P_L, S_M_U, S_M_L
      REAL(fpk)    :: EPS, CONST, SECBAR, WDEL, ZDEL, ZWDEL, F1

!  initialise indices

      N = GIVEN_LAYER
      M = FOURIER

!  No particular solution beyond the cutoff layer, Or no scattering in this layer...
!  ... Zero the boundary layer values and exit

      IF ( N .GT. LAYER_PIS_CUTOFF(IBEAM) .OR. .NOT. DO_LAYER_SCATTERING(M,N) ) THEN
        DO I = 1, NSTREAMS_2
          WUPPER(I,N) = ZERO
          WLOWER(I,N) = ZERO
        ENDDO
        RETURN
      ENDIF

!  constants for the layer

      SECBAR = AVERAGE_SECANT(N,IBEAM)
      CONST  = INITIAL_TRANS (N,IBEAM)
      WDEL   = T_DELT_MUBAR  (N,IBEAM)

!  3. Optical depth integrations for the discrete ordinate solution
!     =============================================================


      DO AA = 1, NSTREAMS
        GAMMA_P(AA,N) = SECBAR + KEIGEN(AA,N)
        GAMMA_M(AA,N) = SECBAR - KEIGEN(AA,N)
        ZDEL  = T_DELT_EIGEN(AA,N)
        ZWDEL = ZDEL * WDEL
        IF ( ABS(GAMMA_M(AA,N)) .LT. TAYLOR_SMALL ) THEN
           EPS = GAMMA_M(AA,N)
           CALL TAYLOR_SERIES_1 ( TAYLOR_ORDER, EPS, DELTAUS(N), WDEL, ONE, CFUNC(AA,N) )
        ELSE
           CFUNC(AA,N) =  ( ZDEL - WDEL ) / GAMMA_M(AA,N)
        ENDIF
        DFUNC(AA,N)  = ( ONE - ZWDEL ) / GAMMA_P(AA,N)
      ENDDO

!  4. Form quantities independent of optical depth
!     ============================================

!  set up help arrays (independent of eigenvector)

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

!  For each eigenstream, get the terms ATERM_SAVE and BTERM_SAVE

      DO AA = 1, NSTREAMS
        SUM_LA = ZERO
        SUM_LB = ZERO
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          TPA = QUAD_WTS(I)*(DPI(I)*XPOS(I,AA,N)+DMI(I)*XPOS(I1,AA,N))
          TMA = QUAD_WTS(I)*(DMI(I)*XPOS(I,AA,N)+DPI(I)*XPOS(I1,AA,N))
          SUM_LA  = SUM_LA + TPA
          SUM_LB  = SUM_LB + TMA
        ENDDO
        ATERM_SAVE(AA,N) = SUM_LA / NORM_SAVED(N,AA)
        BTERM_SAVE(AA,N) = SUM_LB / NORM_SAVED(N,AA)
        AGM(AA,N)  = ATERM_SAVE(AA,N) / GAMMA_M(AA,N) ! Watch for smallness
        BGP(AA,N)  = BTERM_SAVE(AA,N) / GAMMA_P(AA,N)
      ENDDO

!  5. Green function multipliers
!     ==========================

!  For each eigenstream

      DO AA = 1, NSTREAMS
        GFUNC_DN(AA,N) = CFUNC(AA,N) * ATERM_SAVE(AA,N) * CONST
        GFUNC_UP(AA,N) = DFUNC(AA,N) * BTERM_SAVE(AA,N) * CONST
      ENDDO

!  6. Set particular integral from Green function expansion
!     =====================================================

!  particular integrals at lower and upper boundaries

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

!  Finish

      RETURN
END SUBROUTINE LIDORT_GBEAM_SOLUTION

!

SUBROUTINE LIDORT_GBEAM_USERSOLUTION                             &
         ( DO_UPWELLING, DO_DNWELLING, DO_OBSERVATION_GEOMETRY,  & ! input
           N_USER_STREAMS, NMOMENTS,                             & ! input
           GIVEN_LAYER, FOURIER, IBEAM, FLUX_FACTOR,             & ! input
           DO_LAYER_SCATTERING, LAYER_PIS_CUTOFF,                & ! input
           OMEGA_MOMS, U_LEG_M, U_LEG_P, LEG0_M,                 & ! input
           U_WPOS1, U_WNEG1, W_HELP )                              ! Output

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAX_USER_STREAMS, MAXMOMENTS, MAXLAYERS, MAXBEAMS, &
                              ZERO, PI4

      IMPLICIT NONE

!  subroutine input arguments
!  --------------------------

!  Obervational Geometry flag

      LOGICAL  , intent(in)  :: DO_OBSERVATION_GEOMETRY

!  Direction flags

      LOGICAL  , intent(in)  :: DO_UPWELLING
      LOGICAL  , intent(in)  :: DO_DNWELLING

!  Number of streams

      INTEGER  , intent(in)  :: N_USER_STREAMS

!  Number of moments

      INTEGER  , intent(in)  :: NMOMENTS

!  Given layer index and Fourier number (inputs)

      INTEGER  , intent(in)  :: GIVEN_LAYER
      INTEGER  , intent(in)  :: FOURIER
      INTEGER  , intent(in)  :: IBEAM

!  Solar Flux

      REAL(fpk), intent(in)  :: FLUX_FACTOR

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Saved array involving product of OMEGA and phase function moments

      REAL(fpk), intent(in)  :: OMEGA_MOMS ( MAXLAYERS, 0:MAXMOMENTS )

!  Legendre functions on User defined polar angles

      REAL(fpk), intent(in)  :: U_LEG_P(MAX_USER_STREAMS,0:MAXMOMENTS)
      REAL(fpk), intent(in)  :: U_LEG_M(MAX_USER_STREAMS,0:MAXMOMENTS)

!  Legendre polynomials, LEG0_M holds stored quantities.

      REAL(fpk), intent(in)  :: LEG0_M(0:MAXMOMENTS,MAXLAYERS,MAXBEAMS)

!  Subroutine output arguments
!  ---------------------------

!  Saved help variables

      REAL(fpk), intent(out) :: W_HELP(0:MAXMOMENTS)

!mick fix 6/29/11 - change two outputs from "out" to "inout"

!  Particular beam solutions at user-defined stream angles

      REAL(fpk), intent(inout) :: U_WPOS1(MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(inout) :: U_WNEG1(MAX_USER_STREAMS,MAXLAYERS)

!  Local variables
!  ---------------

      INTEGER      :: UM, L, N, M
      REAL(fpk)    :: POS1, POS2, F1

!  Layer and Fourier

      N = GIVEN_LAYER
      M = FOURIER

!  No particular solution beyond the cutoff layer
!  Or no scattering in this layer...
!  ... Zero the user solutions and exit

      IF ( N .GT. LAYER_PIS_CUTOFF(IBEAM) .OR. &
              .NOT. DO_LAYER_SCATTERING(M,N) ) THEN
        DO UM = 1, N_USER_STREAMS
          IF ( DO_UPWELLING ) THEN
            U_WPOS1(UM,N) = ZERO
          ENDIF
          IF ( DO_DNWELLING ) THEN
            U_WNEG1(UM,N) = ZERO
          ENDIF
        ENDDO
        RETURN
      ENDIF

!  Scattering solutions
!  ====================

!  For each moment do inner sum over computational angles

      F1 = FLUX_FACTOR / PI4
      DO L = M, NMOMENTS
        W_HELP(L) = LEG0_M(L,N,IBEAM)*OMEGA_MOMS(N,L)*F1
      ENDDO

!  Now sum over all harmonic contributions at each user-defined stream
!   New: Observational Geometry, tie to the IB result.

      IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
        DO UM = 1, N_USER_STREAMS
          POS1 = ZERO
          POS2 = ZERO
          DO L = M, NMOMENTS
            POS1 = POS1 + W_HELP(L)*U_LEG_P(UM,L)
            POS2 = POS2 + W_HELP(L)*U_LEG_M(UM,L)
          ENDDO
          U_WPOS1(UM,N) = POS1
          U_WNEG1(UM,N) = POS2
        ENDDO
      ELSE
        POS1 = ZERO ; POS2 = ZERO
        DO L = M, NMOMENTS
          POS1 = POS1 + W_HELP(L)*U_LEG_P(IBEAM,L)
          POS2 = POS2 + W_HELP(L)*U_LEG_M(IBEAM,L)
        ENDDO
        U_WPOS1(IBEAM,N) = POS1
        U_WNEG1(IBEAM,N) = POS2
      ENDIF

!  Finish

      RETURN
END SUBROUTINE LIDORT_GBEAM_USERSOLUTION

!  End Module

end module lidort_solutions

