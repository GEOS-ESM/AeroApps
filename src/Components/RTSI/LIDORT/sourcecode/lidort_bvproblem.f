C ###########################################################
C #                                                         #
C #                    THE LIDORT FAMILY                    #
C #                                                         #
C #      (LInearized Discrete Ordinate Radiative Transfer)  #
C #       --         -        -        -          -         #
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
C # Regular BVP: Subroutines in this Module                     #
C #                                                             #
C #            BVP_MATRIXSETUP_MASTER      (master)             #
C #             BVP_SURFACE_SETUP_HOM                           #
C #             BVP_MATRIX_INIT                                 #
C #             BVP_MATRIX_SETUP                                #
C #             BVP_MATRIX_SVD                                  #
C #                                                             #
C #            BVP_SOLUTION_MASTER      (master)                #
C #             BVP_SURFACE_SETUP_SRC                           #
C #             BVP_COLUMN_SETUP                                #
C #             BVP_BACKSUB                                     #
C #                                                             #
C # Telescped BVP: Subroutines in this Module                   #
C #                                                             #
C #            BVPTEL_MATRIXSETUP_MASTER      (master)          #
C #             BVPTEL_MATRIX_INIT                              #
C #             BVPTEL_MATRIX_SETUP                             #
C #             BVPTEL_MATRIX_SVD                               #
C #                                                             #
C #            BVPTEL_SOLUTION_MASTER      (master)             #
C #             BVPTEL_COLUMN_SETUP                             #
C #             BVPTEL_BACKSUB                                  #
C #                                                             #
C ###############################################################

      SUBROUTINE BVP_MATRIXSETUP_MASTER
     I     ( DO_INCLUDE_SURFACE,
     I       FOURIER_COMPONENT,
     I       SURFACE_FACTOR,
     O       STATUS )

C  Include files
C  -------------

C  Additional sums for the final albedo-reflecting layer

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include file of input variables (control/model/surface/geophys)

      INCLUDE '../includes/LIDORT_INPUTS.VARS'

C  Include file of book-keeping inputs

      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file of solution stuff (i/o to this module)

      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  control input
C  -------------

      LOGICAL          DO_INCLUDE_SURFACE
      INTEGER          FOURIER_COMPONENT
      DOUBLE PRECISION SURFACE_FACTOR
      INTEGER          STATUS

C  local variables
C  ---------------

      CHARACTER*(70)   MAIL, TRACE
      INTEGER          STATUS_SUB
      INTEGER          I, I1, C0, LAY

C  Vector matrix masks

      IF (FOURIER_COMPONENT .EQ. 0 ) THEN
        DO LAY = 1, NLAYERS
          C0 = (LAY-1)*NSTREAMS_2
          DO I = 1, NSTREAMS
            I1 = I+NSTREAMS
            LCONMASK(I,LAY) = C0+I
            MCONMASK(I,LAY) = C0+I1
          ENDDO
        ENDDO
      ENDIF

C  Additional setups for the albedo layer

      CALL BVP_SURFACE_SETUP_HOM
     I  ( DO_INCLUDE_SURFACE,
     I    FOURIER_COMPONENT,
     I    SURFACE_FACTOR )

C  initialize compression matrix (Do this for every Fourier component)

      CALL BVP_MATRIX_INIT ( FOURIER_COMPONENT )

C  set up boundary values matrix in compressed form (the "A" as in AX=B)

      CALL BVP_MATRIX_SETUP
     I     ( DO_INCLUDE_SURFACE, FOURIER_COMPONENT )

C  SVD decomposition of compressed boundary values matrix

      CALL BVP_MATRIX_SVD ( STATUS_SUB )

C  error tracing

      IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
        MAIL = 'Error return from module BVP_MATRIX_SVD'
        TRACE= 'Call in BVP_MATRIXSETUP_MASTER'
        STATUS = LIDORT_SERIOUS
        CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
        RETURN
      ENDIF

C  finish

      RETURN
      END

C

      SUBROUTINE BVP_SURFACE_SETUP_HOM
     I     ( DO_INCLUDE_SURFACE,
     I       FOURIER_COMPONENT,
     I       SURFACE_FACTOR )

C  Additional sums for the final surface-reflecting layer

C  Include files
C  -------------

C  include file of dimensiopns and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include file of input variables (control/model/surface/geophys)
C  Include file of book-keeping inputs

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file of solution stuff (i/o to this module)

      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  control input
C  -------------

      LOGICAL          DO_INCLUDE_SURFACE
      INTEGER          FOURIER_COMPONENT
      DOUBLE PRECISION SURFACE_FACTOR

C  local variables

      DOUBLE PRECISION REFL_P, REFL_M, FACTOR_KERNEL
      INTEGER          I, AA, J, K

C  Initialization
C  ==============

C  Zero total reflected contributions

      DO I = 1, NSTREAMS
        DO AA = 1, NSTREAMS
          R2_HOMP(I,AA) = ZERO
          R2_HOMM(I,AA) = ZERO
        ENDDO
      ENDDO

C  Return with Zeroed values if albedo flag not set

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

C  ==================================================
C  ==================================================
C  Integrate Downward streams of particular solutions
C  ==================================================
C  ==================================================

      DO K = 1, N_BRDF_KERNELS

        FACTOR_KERNEL = SURFACE_FACTOR * BRDF_FACTORS(K)

C  For Lambertian reflectance, all streams are the same
C  ----------------------------------------------------

        IF ( LAMBERTIAN_KERNEL_FLAG(K) ) THEN

          IF ( FOURIER_COMPONENT .EQ. 0 ) THEN

            DO AA = 1, NSTREAMS
              REFL_P = ZERO
              REFL_M = ZERO
              DO J = 1, NSTREAMS
                REFL_P = REFL_P + XPOS(J,AA,NLAYERS) * AX(J)
                REFL_M = REFL_M + XNEG(J,AA,NLAYERS) * AX(J)
              ENDDO
              REFL_P = REFL_P * FACTOR_KERNEL
              REFL_M = REFL_M * FACTOR_KERNEL
              DO I = 1, NSTREAMS
                RA_HOMP(I,AA,K) = REFL_P
                RA_HOMM(I,AA,K) = REFL_M
              ENDDO
            ENDDO

C  zero these terms otherwise

          ELSE
            DO AA = 1, NSTREAMS
              DO I = 1, NSTREAMS
                RA_HOMP(I,AA,K) = ZERO
                RA_HOMM(I,AA,K) = ZERO
              ENDDO
            ENDDO
          ENDIF

C  For bidirectional reflecting surface
C  ------------------------------------

        ELSE

C  Save some quantities

          DO I = 1, NSTREAMS
            DO J = 1, NSTREAMS
              AXBID(I,J,K) = AX(J)*BIREFLEC(K,J,I)
            ENDDO
          ENDDO

C  get the reflected quantities

          DO AA = 1, NSTREAMS
            DO I = 1, NSTREAMS
              REFL_P = ZERO
              REFL_M = ZERO
              DO J = 1, NSTREAMS
                REFL_P = REFL_P + XPOS(J,AA,NLAYERS) * AXBID(I,J,K)
                REFL_M = REFL_M + XNEG(J,AA,NLAYERS) * AXBID(I,J,K)
              ENDDO
              RA_HOMP(I,AA,K) = REFL_P * FACTOR_KERNEL
              RA_HOMM(I,AA,K) = REFL_M * FACTOR_KERNEL
            ENDDO
          ENDDO

        ENDIF

C  End loop over albedo kernels

      ENDDO  

C  Total terms
C  -----------

      DO K = 1, N_BRDF_KERNELS
        DO I = 1, NSTREAMS
          DO AA = 1, NSTREAMS
            R2_HOMP(I,AA) = R2_HOMP(I,AA) + RA_HOMP(I,AA,K)
            R2_HOMM(I,AA) = R2_HOMM(I,AA) + RA_HOMM(I,AA,K)
          ENDDO
        ENDDO
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE BVP_MATRIX_INIT ( FOURIER_COMPONENT )

C  Initialize the compressed matrix

C  include files
C  -------------

C  include file of dimensiopns and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include file of input variables (control/model/surface/geophys)
C  Include file of book-keeping inputs

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of main model variables (local to this module)

      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  Input
C  -----

      INTEGER             FOURIER_COMPONENT

C  Local variables
C  ---------------

      INTEGER  NMIN(MAXTOTAL), NMAX(MAXTOTAL)
      INTEGER  I, J, N3, JS, JF, IS, LF, L, KALL, I1

C  special case

      IF ( NLAYERS .EQ. 1 ) THEN
        DO I = 1, NTOTAL
          DO J = 1, NTOTAL
            SMAT2(I,J)     = ZERO
          ENDDO
        ENDDO
        RETURN
      ENDIF

c  compression row indices

      DO J = 1, N_SUPDIAG + 1
        NMIN(J) = 1
      ENDDO
      DO J = N_SUPDIAG + 2, NTOTAL
        NMIN(J) = J - N_SUPDIAG
      ENDDO
      DO J = 1, NTOTAL - N_SUBDIAG
        NMAX(J) = J + N_SUBDIAG
      ENDDO
      DO J = NTOTAL - N_SUBDIAG + 1, NTOTAL
        NMAX(J) = NTOTAL
      ENDDO

      KALL = N_SUBDIAG + N_SUPDIAG + 1
      DO I = 1, NTOTAL
        DO J = 1, NTOTAL
          IF ( (I.GE.NMIN(J)) .AND. (I.LE.NMAX(J)) ) THEN
            BMAT_ROWMASK(I,J) = KALL + I - J
          ENDIF
        ENDDO
      ENDDO

C  compression matrix zeroing

      N3 = NSTREAMS_2 + NSTREAMS
      LF = NLAYERS - 2

C  upper band top

      JS = NSTREAMS_2 + 1
      JF = N3 - 1
      DO I = 1, NSTREAMS
        DO J = JS, JF + I
          BANDMAT2(BMAT_ROWMASK(I,J),J) = ZERO
        ENDDO
      ENDDO

C  upper band

      DO L = 1, LF
        IS = L*NSTREAMS_2 - NSTREAMS + 1
        JS = IS + N3
        JF = JS - 1
        DO I = 1, NSTREAMS_2-1
          I1 = I + IS
          DO J = JS, JF + I
            BANDMAT2(BMAT_ROWMASK(I1,J),J) = ZERO
          ENDDO
        ENDDO
      ENDDO

C  lower band

      DO L = 1, LF
        IS = L*NSTREAMS_2 + NSTREAMS
        JS = IS - N3 + 1
        JF = IS - NSTREAMS
        DO I = 1, NSTREAMS_2-1
          I1 = I + IS
          DO J = JS + I, JF
            BANDMAT2(BMAT_ROWMASK(I1,J),J) = ZERO
          ENDDO
        ENDDO
      ENDDO

C  lower band bottom

      JS = LF * NSTREAMS_2 + 1
      IS = JS + N3 - 1
      JF = IS - NSTREAMS
      DO I = 1, NSTREAMS
        I1 = I + IS
        DO J = JS + I, JF
          BANDMAT2(BMAT_ROWMASK(I1,J),J) = ZERO
        ENDDO
      ENDDO

C  finish

      RETURN
      END

C

      SUBROUTINE BVP_MATRIX_SETUP
     I     (  DO_INCLUDE_SURFACE, FOURIER_COMPONENT )

C  Fills up the compressed matrix directly

C  include files
C  -------------

C  include file of dimensiopns and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include file of input variables (control/model/surface/geophys)
C  Include file of book-keeping inputs

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of main model variables (local to this module)

      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  input arguments
C  ---------------

      LOGICAL             DO_INCLUDE_SURFACE
      INTEGER             FOURIER_COMPONENT

C  Local variables
C  ---------------

      INTEGER             I,EP,EM,N,N1,I1,AA
      INTEGER             C0,CE_OFFSET,CEM,CEP,CEM1,CEP1,CM,CP
      DOUBLE PRECISION    XPNET, XMNET

C  If Nlayers = 1, go to special case

      IF ( NLAYERS .EQ. 1 ) GO TO 345

C  top BC for layer 1: no downward diffuse radiation

      N = 1
      DO I = 1, NSTREAMS
        DO EP = 1, NSTREAMS
          EM = EP + NSTREAMS
          BANDMAT2(BMAT_ROWMASK(I,EP),EP)  =
     &                   XPOS(I,EP,N)
          BANDMAT2(BMAT_ROWMASK(I,EM),EM)  =
     &                   XNEG(I,EP,N)*T_DELT_EIGEN(EP,N)
        ENDDO
      ENDDO

C  intermediate layer boundaries (will not be done if NLAYERS = 1 )

      C0 = - NSTREAMS
      DO N = 2, NLAYERS
        N1 = N - 1
        C0   = C0 + NSTREAMS_2
        CE_OFFSET = C0 - NSTREAMS
        DO I = 1, NSTREAMS_2
          CM = C0 + I
          DO EP = 1, NSTREAMS
            CEP = CE_OFFSET + EP
            CEM = CEP + NSTREAMS
            CEP1 = CEP + NSTREAMS_2
            CEM1 = CEM + NSTREAMS_2
            BANDMAT2(BMAT_ROWMASK(CM,CEP),CEP)   =
     &                   T_DELT_EIGEN(EP,N1)*XPOS(I,EP,N1)
            BANDMAT2(BMAT_ROWMASK(CM,CEM),CEM)   =
     &                   XNEG(I,EP,N1)
            BANDMAT2(BMAT_ROWMASK(CM,CEP1),CEP1) =
     &                   -XPOS(I,EP,N)
            BANDMAT2(BMAT_ROWMASK(CM,CEM1),CEM1) =
     &                   -T_DELT_EIGEN(EP,N)*XNEG(I,EP,N)
          ENDDO
        ENDDO
      ENDDO

C  bottom BC (with albedo additions if flagged)

      N = NLAYERS
      C0 = C0 + NSTREAMS_2
      CE_OFFSET = C0 - NSTREAMS

      IF ( DO_INCLUDE_SURFACE ) THEN
        DO I = 1, NSTREAMS
          CP = C0 + I
          I1 = I + NSTREAMS
          DO AA = 1, NSTREAMS
            CEP = CE_OFFSET + AA
            CEM = CEP + NSTREAMS
            XPNET = XPOS(I1,AA,N) - R2_HOMP(I,AA)
            XMNET = XNEG(I1,AA,N) - R2_HOMM(I,AA)
            BANDMAT2(BMAT_ROWMASK(CP,CEP),CEP) =
     &                   T_DELT_EIGEN(AA,N) * XPNET
            BANDMAT2(BMAT_ROWMASK(CP,CEM),CEM) =
     &                   XMNET
          ENDDO
        ENDDO
      ELSE
        DO I = 1, NSTREAMS
          CP = C0 + I
          I1 = I + NSTREAMS
          DO AA = 1, NSTREAMS
            CEP = CE_OFFSET + AA
            CEM = CEP + NSTREAMS
            BANDMAT2(BMAT_ROWMASK(CP,CEP),CEP) =
     &                  T_DELT_EIGEN(AA,N) * XPOS(I1,AA,N)
            BANDMAT2(BMAT_ROWMASK(CP,CEM),CEM) =
     &                  XNEG(I1,AA,N)
          ENDDO
        ENDDO
      ENDIF

C  normal completion

      RETURN

C  special case

345   CONTINUE

C  top BC for layer 1: no downward diffuse radiation
C  -------------------------------------------------

      N = 1
      DO I = 1, NSTREAMS
        DO EP = 1, NSTREAMS
          EM = EP + NSTREAMS
          SMAT2(I,EP) = XPOS(I,EP,N)
          SMAT2(I,EM) = XNEG(I,EP,N)*T_DELT_EIGEN(EP,N)
        ENDDO
      ENDDO

C  bottom BC (with albedo additions)
C  ---------------------------------

      IF ( DO_INCLUDE_SURFACE ) THEN
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          CP = NSTREAMS + I
          DO EP = 1, NSTREAMS
            CEP = EP
            CEM = CEP + NSTREAMS
            XPNET = XPOS(I1,EP,N) - R2_HOMP(I,EP)
            XMNET = XNEG(I1,EP,N) - R2_HOMM(I,EP)
            SMAT2(CP,CEP) = T_DELT_EIGEN(EP,N) * XPNET
            SMAT2(CP,CEM) = XMNET
          ENDDO
        ENDDO
      ELSE
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          CP = NSTREAMS + I
          DO EP = 1, NSTREAMS
            CEP = EP
            CEM = CEP + NSTREAMS
            XPNET = XPOS(I1,EP,N)
            XMNET = XNEG(I1,EP,N)
            SMAT2(CP,CEP) = T_DELT_EIGEN(EP,N) * XPNET
            SMAT2(CP,CEM) = XMNET
          ENDDO
        ENDDO
      ENDIF

C  normal return and finish

      RETURN
      END

C

      SUBROUTINE BVP_MATRIX_SVD ( STATUS )

C  Solves the boundary value problem.

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include file of input variables (control/model/surface/geophys)
C  Include file of book-keeping inputs

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of main model variable (I/O to this module)

      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  subroutine arguments
C  --------------------

C  status flag

      INTEGER          STATUS

C  local variables
C  ---------------

      CHARACTER*(70)   MAIL, TRACE
      INTEGER          INFO
      CHARACTER*3      CI

C  Intialize status

      STATUS = LIDORT_SUCCESS

C  SVD the BVP matrix: With compression (multilayers)
C  --------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

C  LAPACK LU-decomposition for band matrix

        CALL DGBTRF
     &     ( NTOTAL, NTOTAL, N_SUBDIAG, N_SUPDIAG,
     &       BANDMAT2, MAXBANDTOTAL, IPIVOT, INFO )

C  (Error tracing)

        IF ( INFO .GT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MAIL  = 'Singular matrix, u(i,i)=0, for i = '//CI
          TRACE = 'DGBTRF call in BVP_MATRIX_SVD'
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ELSE IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MAIL  = 'argument i illegal value, for i = '//CI
          TRACE = 'DGBTRF call in BVP_MATRIX_SVD'
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ENDIF

C  SVD the BVP matrix: No compression, Single Layer only
C  -----------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

C  LAPACK LU-decomposition for  matrix

        CALL DGETRF
     &     ( NTOTAL, NTOTAL, SMAT2, MAXSTREAMS_2, SIPIVOT, INFO )

C  (Error tracing)

        IF ( INFO .GT. 0 ) THEN
          MAIL='BVP LU decomposition (nlayers=1) (DGETRF)'
          CALL LIDORT_LAPACK_ERROR ( INFO, 1, MAIL, STATUS )
          IF ( STATUS .EQ. LIDORT_SERIOUS ) RETURN
        ENDIF

      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE BVP_SOLUTION_MASTER
     I       ( DO_INCLUDE_SURFACE,
     I         DO_INCLUDE_DIRECTBEAM,
     I         DO_INCLUDE_SURFEMISS,
     I         FOURIER_COMPONENT, IPARTIC,
     I         SURFACE_FACTOR,
     O         STATUS )

C  include files
C  -------------

C  include file of dimensiopns and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include file of input variables (control/model/surface/geophys)
C  Include file of book-keeping inputs

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of solution variables (local to this module)

      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  input arguments
C  ---------------

C  inclusion flags

      LOGICAL          DO_INCLUDE_DIRECTBEAM
      LOGICAL          DO_INCLUDE_SURFACE
      LOGICAL          DO_INCLUDE_SURFEMISS

C  Fourier component and beam number

      INTEGER          FOURIER_COMPONENT, IPARTIC

C  surface factor = 1+delta(m,0)

      DOUBLE PRECISION SURFACE_FACTOR

C  output status
C  -------------

      INTEGER          STATUS

C  Local variables
C  ---------------

      CHARACTER*(70)   MAIL, TRACE
      INTEGER          STATUS_SUB

C  Regular BVP using compressed-band matrices, etc..
C  ==================================================

C  --Additional setups for the albedo layer

      CALL BVP_SURFACE_SETUP_SRC
     I       ( DO_INCLUDE_SURFACE,
     I         FOURIER_COMPONENT, IPARTIC,
     I         SURFACE_FACTOR )

C  --set up Column for solution vector (the "B" as in AX=B)

      CALL BVP_COLUMN_SETUP
     I     ( DO_INCLUDE_SURFACE,
     I       DO_INCLUDE_DIRECTBEAM,
     I       DO_INCLUDE_SURFEMISS, IPARTIC )

C  --Solve the boundary problem for this Fourier component (back substitution)

      CALL BVP_BACKSUB
     I   ( FOURIER_COMPONENT, IPARTIC, 
     O     STATUS_SUB )

C  error tracing

      IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
        MAIL  = 'Error return from module BVP_BACKSUB'
        TRACE = 'Call #1 in BVP_SOLUTION_MASTER'
        STATUS = LIDORT_SERIOUS
        CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE BVP_SURFACE_SETUP_SRC
     I     ( DO_INCLUDE_SURFACE,
     I       FOURIER_COMPONENT, IPARTIC,
     I       SURFACE_FACTOR )

C  Additional sums for the final surface-reflecting layer

C  Include files
C  -------------

C  include file of dimensiopns and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include file of input variables (control/model/surface/geophys)
C  Include file of book-keeping inputs

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file of solution stuff (i/o to this module)

      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  control input
C  -------------

      LOGICAL          DO_INCLUDE_SURFACE
      INTEGER          FOURIER_COMPONENT, IPARTIC
      DOUBLE PRECISION SURFACE_FACTOR

C  local variables

      DOUBLE PRECISION REFL_B, FACTOR_KERNEL
      INTEGER          I, J, K

C  Initialization
C  ==============

C  Zero total reflected contributions

      DO I = 1, NSTREAMS
        R2_PARTIC(I)    = ZERO
      ENDDO

C  Return with Zeroed values if albedo flag not set

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

C  ==================================================
C  ==================================================
C  Integrate Downward streams of particular solutions
C  ==================================================
C  ==================================================

      DO K = 1, N_BRDF_KERNELS

        FACTOR_KERNEL = SURFACE_FACTOR * BRDF_FACTORS(K)

C  For Lambertian reflectance, all streams are the same
C  ----------------------------------------------------

        IF ( LAMBERTIAN_KERNEL_FLAG(K) ) THEN

          IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
            REFL_B = ZERO
            DO J = 1, NSTREAMS
              REFL_B = REFL_B + AX(J)*WLOWER(J,NLAYERS)
            ENDDO
            REFL_B = REFL_B * FACTOR_KERNEL
            DO I = 1, NSTREAMS
              RA_PARTIC(I,K) = REFL_B
            ENDDO
          ELSE
            DO I = 1, NSTREAMS
              RA_PARTIC(I,K) = ZERO
            ENDDO
          ENDIF

C  For bidirectional reflecting surface
C  ------------------------------------

        ELSE

          DO I = 1, NSTREAMS
            REFL_B = ZERO
            DO J = 1, NSTREAMS
              REFL_B = REFL_B + WLOWER(J,NLAYERS) * AXBID(I,J,K)
            ENDDO
            RA_PARTIC(I,K) = FACTOR_KERNEL * REFL_B
          ENDDO

        ENDIF

C  End loop over albedo kernels

      ENDDO  

C  Total terms
C  -----------

      DO K = 1, N_BRDF_KERNELS
        DO I = 1, NSTREAMS
          R2_PARTIC(I) = R2_PARTIC(I) + RA_PARTIC(I,K)
         ENDDO
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE BVP_COLUMN_SETUP
     I     ( DO_INCLUDE_SURFACE,
     I       DO_INCLUDE_DIRECTBEAM,
     I       DO_INCLUDE_SURFEMISS, IPARTIC )

C  include files
C  -------------

C  include file of dimensiopns and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include file of input variables (control/model/surface/geophys)
C  Include file of book-keeping inputs

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of solution variables (local to this module)

      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  input arguments
C  ---------------

      LOGICAL          DO_INCLUDE_DIRECTBEAM
      LOGICAL          DO_INCLUDE_SURFACE
      LOGICAL          DO_INCLUDE_SURFEMISS
      INTEGER          IPARTIC

C  Local variables
C  ---------------

      INTEGER          I, I1, LAY, LAY1, C0, CM
      DOUBLE PRECISION FPSURF

C  If Nlayers = 1, go to special case

      IF ( NLAYERS .EQ. 1 ) GO TO 345

C  zero column vector

      DO I = 1, NTOTAL
        COL2(I,IPARTIC) = ZERO
      ENDDO

C  Upper boundary for layer 1: no downward diffuse radiation
C  ---------------------------------------------------------

      LAY = 1
      DO I = 1, NSTREAMS
        COL2(I,IPARTIC)   = - WUPPER(I,LAY)
      ENDDO

C  intermediate layer boundaries
C  -----------------------------

      DO LAY = 2, NLAYERS
        LAY1 = LAY - 1
          C0 = LAY1*NSTREAMS_2 - NSTREAMS
        DO I = 1, NSTREAMS_2
          CM = C0 + I
          COL2(CM,IPARTIC) = WUPPER(I,LAY) - WLOWER(I,LAY1)
        ENDDO
      ENDDO

C  lowest (surface) boundary with albedo (diffuse radiation terms only)
C  -------------------------------------

      LAY = NLAYERS
      C0 = (LAY-1)*NSTREAMS_2 + NSTREAMS

C  with non-zero surface term, include integrated downward reflectances

      IF ( DO_INCLUDE_SURFACE ) THEN

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          CM = C0 + I
          COL2(CM,IPARTIC) = - WLOWER(I1,LAY) + R2_PARTIC(I)
        ENDDO

C  no surface term, similar code excluding integrated reflectance

      ELSE

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          CM = C0 + I
          COL2(CM,IPARTIC) = - WLOWER(I1,LAY)
        ENDDO

      ENDIF

C  Add direct beam solution (only to final level)
C  ----------------------------------------------

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        IF ( DO_INCLUDE_SURFACE ) THEN
          DO I = 1, NSTREAMS
            CM = C0 + I
            COL2(CM,IPARTIC) = COL2(CM,IPARTIC) + DIRECT_BEAM(I,IPARTIC)
          ENDDO
        ENDIF
      ENDIF

C  Add thermal emission of ground surface (only to final level)
C  ------------------------------------------------------------

      IF ( DO_INCLUDE_SURFEMISS ) THEN
        FPSURF = SURFBB
        DO I = 1, NSTREAMS
          CM = C0 + I
          COL2(CM,IPARTIC) = COL2(CM,IPARTIC) + FPSURF * EMISSIVITY(I)
        ENDDO
      ENDIF

C  debug

c      do i = 1,ntotal
c        write(*,*)i,COL2(i,IPARTIC)
c      enddo
c      pause

C  normal completion

      RETURN

C  special case - Only one layer
C  =============================

345   CONTINUE

C  zero column vector

      DO I = 1, NTOTAL
        SCOL2(I,IPARTIC) = ZERO
      ENDDO

C  Upper boundary for layer 1: no downward diffuse radiation

      LAY = 1
      DO I = 1, NSTREAMS
        SCOL2(I,IPARTIC)   = - WUPPER(I,LAY)
      ENDDO

C  lowest (surface) boundary with albedo (diffuse radiation terms only)
C  with non-zero albedo, include integrated downward reflectances
C  no albedo, similar code excluding integrated reflectance

      C0 = NSTREAMS
      IF ( DO_INCLUDE_SURFACE ) THEN
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          CM = C0 + I
          SCOL2(CM,IPARTIC) = - WLOWER(I1,LAY) + R2_PARTIC(I)
        ENDDO
      ELSE
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          CM = C0 + I
          SCOL2(CM,IPARTIC) = - WLOWER(I1,LAY)
        ENDDO
      ENDIF

C  Add direct beam solution (only to final level)

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        IF ( DO_INCLUDE_SURFACE ) THEN
          DO I = 1, NSTREAMS
            CM = C0 + I
            SCOL2(CM,IPARTIC) = SCOL2(CM,IPARTIC)+DIRECT_BEAM(I,IPARTIC)
          ENDDO
        ENDIF
      ENDIF

C  Add thermal emission of ground surface (only to final level)

      IF ( DO_INCLUDE_SURFEMISS ) THEN
        FPSURF = SURFBB
        DO I = 1, NSTREAMS
          CM = C0 + I
          SCOL2(CM,IPARTIC) = SCOL2(CM,IPARTIC) + FPSURF*EMISSIVITY(I)
        ENDDO
      ENDIF

C  finish

      RETURN
      END

C

      SUBROUTINE BVP_BACKSUB
     I     ( FOURIER, IPARTIC,
     O       STATUS )

C  Solves the boundary value problem.

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include file of input variables (control/model/surface/geophys)
C  Include file of book-keeping inputs

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of main model variable (I/O to this module)

      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  subroutine arguments
C  --------------------

C  FOURIER and Beam index

      INTEGER          FOURIER, IPARTIC

C  status flag

      INTEGER          STATUS

C  local variables
C  ---------------

      INTEGER          I, C0, LAY, K, N, K1
      CHARACTER*(70)   MAIL, TRACE
      INTEGER          INFO
      CHARACTER*3      CI

C  Intialize status

      STATUS = LIDORT_SUCCESS

C  BVP back-substitution: With compression (multilayers)
C  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

C  LAPACK substitution (DGBTRS) using RHS column vector COL2

        CALL DGBTRS
     &     ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, IPARTIC,
     &        BANDMAT2, MAXBANDTOTAL, IPIVOT,
     &        COL2, MAXTOTAL, INFO )

C  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MAIL  = 'argument i illegal value, for i = '//CI
          TRACE = 'DGBTRS call in BVP_BACKSUB'
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
         ENDIF

C  Set integration constants LCON and MCON for -/+ eigensolutions, all layers

        DO LAY = 1, NLAYERS
          C0 = (LAY-1)*NSTREAMS_2
          DO K = 1, NSTREAMS
            K1 = K + NSTREAMS
            LCON(K,LAY) = COL2(C0+K, IPARTIC)
            MCON(K,LAY) = COL2(C0+K1,IPARTIC)
          ENDDO
        ENDDO

C  debug

c      IF ( DO_DEBUG_WRITE ) THEN
c       K= 97
c       i = 5
c       IF ( DO_FDTEST ) K= 98
c       if ( fourier.eq.0.and.ipartic.eq.1) then
c         DO N = 1, NLAYERS
c           DO I = 1, NSTREAMS
c           write(98,'(3i4,1p4e17.9)')IPARTIC,N,I, LCON(I,N), MCON(I,N)
c         ENDDO
c         ENDDO
c       ENDIF
c      IF ( DO_FDTEST .and. fourier.EQ.0) PAUSE
c      ENDIF

C  Debug -TOA upwelling for 10 streams. VERY USEFUL INFORMATION
c      do i = 11, 20
c       sum = wupper(i,1)
c       do aa = 1, nstreams
c        sum = sum + lcon(aa,1)*xpos(i,aa,1) +
c     &         mcon(aa,1)*xneg(i,aa,1)*t_delt_eigen(aa,1)
c       enddo
c        if ( fourier.eq.0)write(*,*)i,sum/pi4
c        if ( fourier.gt.0)write(*,*)i,sum/pi2
c      enddo
c      if (fourier.eq.2)pause

C  Solve the boundary problem: No compression, Single Layer only
C  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

C  LAPACK substitution (DGETRS) using RHS column vector SCOL2

        CALL DGETRS
     &     ( 'N', NTOTAL, IPARTIC, SMAT2, MAXSTREAMS_2, SIPIVOT,
     &        SCOL2, MAXSTREAMS_2, INFO )

C  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          MAIL='BVP back-substitution (nlayers=1) (DGETRS)'
          CALL LIDORT_LAPACK_ERROR ( INFO, 1, MAIL, STATUS )
          IF ( STATUS .EQ. LIDORT_SERIOUS ) RETURN
        ENDIF

C  Set integration constants LCON and MCON for -/+ eigensolutions, all layers

        LAY = 1
        DO K = 1, NSTREAMS
          K1 = K + NSTREAMS
          LCON(K,LAY) = SCOL2(K, IPARTIC)
          MCON(K,LAY) = SCOL2(K1,IPARTIC)
        ENDDO

C  debug

       IF ( DO_DEBUG_WRITE ) THEN
        K= 97
        i = 5
        IF ( DO_FDTEST ) K= 98
        if ( fourier.eq.0.and.ipartic.eq.1) then
          DO N = 1, NLAYERS
            DO I = 1, NSTREAMS
           write(K,'(3i4,1p4e17.9)')IPARTIC,N,I, LCON(I,N), MCON(I,N)
          ENDDO
          ENDDO
        ENDIF
c      IF ( DO_FDTEST .and. fourier.EQ.0) PAUSE
       ENDIF

      ENDIF

C  Associated quantities
C  ---------------------

      DO N = 1, NLAYERS
        DO I = 1, NSTREAMS_2
          DO K = 1, NSTREAMS
            LCON_XVEC(I,K,N) = LCON(K,N)*XPOS(I,K,N)
            MCON_XVEC(I,K,N) = MCON(K,N)*XNEG(I,K,N)
          ENDDO
        ENDDO
      ENDDO

C  Finish

      RETURN
      END

c

      SUBROUTINE BVPTEL_MATRIXSETUP_MASTER
     I     ( DO_INCLUDE_SURFACE,
     I       FOURIER_COMPONENT,
     O       STATUS )

C  Sets up the telescoped boundary value problem.
C    Standard case: Fourier > 0. With surface reflection term

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  control input
C  -------------

      LOGICAL          DO_INCLUDE_SURFACE
      INTEGER          FOURIER_COMPONENT
      INTEGER          STATUS

C  local variables
C  ---------------

      CHARACTER*70     MAIL, TRACE
      INTEGER          STATUS_SUB

C  start code
C  ----------

C  Intialize status

      STATUS = LIDORT_SUCCESS

C  initialize compression matrix (Do this for every Fourier component)

      CALL BVPTEL_MATRIX_INIT ( FOURIER_COMPONENT )

C  Not necessary to have reflected solutions in the lower boundary
C  layer - BVP Telescoping only for Lambertian albedo case.

C  set up boundary values matrix in compressed form (the "A" as in AX=B)

      CALL BVPTEL_MATRIX_SETUP
     I     ( DO_INCLUDE_SURFACE, FOURIER_COMPONENT )

C  SVD decomposition of compressed boundary values matrix

      CALL BVPTEL_MATRIX_SVD ( STATUS_SUB )

C  error tracing

      IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
        MAIL = 'Error return from module BVPTEL_MATRIX_SVD'
        TRACE= 'Call in BVPTEL_MATRIXSETUP_MASTER'
        STATUS = LIDORT_SERIOUS
        CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
        RETURN
      ENDIF

C  return

      RETURN
      END

C

      SUBROUTINE BVPTEL_MATRIX_INIT
     I    ( FOURIER_COMPONENT )

C  Initialise the compressed matrix

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include file of input variables
C  include file of bookkeeping inputs

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file of setup variables

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  include file of solution stuff (i/o to this module)

      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  Subroutine argument
C  -------------------

C  Fourier component

      INTEGER             FOURIER_COMPONENT

C  Local variables
C  ---------------

      INTEGER             I, J, KALL, NS, N, N3

C  set up
C  ------

      IF ( DO_BVTEL_INITIAL ) THEN

C  Determine active layers in atmosphere

        NS = 0
        DO N = 1, NLAYERS
         IF ( DO_LAYER_SCATTERING(FOURIER_COMPONENT,N) ) THEN
          NS = NS + 1
          ACTIVE_LAYERS(NS) = N
         ENDIF
        ENDDO
        NLAYERS_TEL = NS

C  size of Reduced BVTEL matrix

        N_BVTELMATRIX_SIZE = 2 * NSTREAMS * NLAYERS_TEL

C  Exit if only one active layer

        IF ( NLAYERS_TEL .EQ. 1 ) GO TO 345

C  Number of sub and super diagonals
  
        N_BVTELMATRIX_SUPDIAG = 3 * NSTREAMS - 1
        N_BVTELMATRIX_SUBDIAG = 3 * NSTREAMS - 1

c  compression row indices

        DO J = 1, N_BVTELMATRIX_SUPDIAG + 1
          NMINTEL(J) = 1
        ENDDO
        DO J = N_BVTELMATRIX_SUPDIAG + 2, N_BVTELMATRIX_SIZE
          NMINTEL(J) = J - N_BVTELMATRIX_SUPDIAG
        ENDDO
        DO J = 1, N_BVTELMATRIX_SIZE - N_BVTELMATRIX_SUBDIAG
          NMAXTEL(J) = J + N_BVTELMATRIX_SUBDIAG
        ENDDO
        N3 = N_BVTELMATRIX_SIZE - N_BVTELMATRIX_SUBDIAG + 1
        DO J = N3, N_BVTELMATRIX_SIZE
          NMAXTEL(J) = N_BVTELMATRIX_SIZE
        ENDDO

        KALL = N_BVTELMATRIX_SUBDIAG + N_BVTELMATRIX_SUPDIAG + 1
        DO I = 1, N_BVTELMATRIX_SIZE
          DO J = 1, N_BVTELMATRIX_SIZE
            IF ( (I.GE.NMINTEL(J)) .AND. (I.LE.NMAXTEL(J)) ) THEN
              BTELMAT_ROWMASK(I,J) = KALL + I - J
            ENDIF
          ENDDO
        ENDDO

C  Avoid fancy zeroing - adopt kludge
C   Potential Danger point

        DO I = 1, MAXBANDTOTAL
          DO J = 1, N_BVTELMATRIX_SIZE
            BANDTELMAT2(I,J) = ZERO
          ENDDO
        ENDDO

C  control point

 345    CONTINUE

C  reset

        DO_BVTEL_INITIAL = .FALSE.

C  end initialization

      ENDIF

C  finish

      RETURN
      END

C

      SUBROUTINE BVPTEL_MATRIX_SETUP
     I     (  DO_INCLUDE_SURFACE, FOURIER_COMPONENT )

C  Fills up the matrix directly (compressed or 1-layer)

C  include files
C  -------------

C  include file of dimensiopns and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include file of input variables
C  include file of bookkeeping inputs

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  Include files of setup variables

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  include file of solution BVP variables (output to this module)

      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  input arguments
C  ---------------

      LOGICAL             DO_INCLUDE_SURFACE
      INTEGER             FOURIER_COMPONENT

C  Local variables
C  ---------------

      INTEGER             I,J,EP,EM,N,N1,I1,AA, NS
      INTEGER             C0,CE_OFFSET,CEM,CEP,CEM1,CEP1,CM,CP
      DOUBLE PRECISION    XPNET, XMNET

C  If Nlayers = 1, go to special case

      IF ( NLAYERS_TEL .EQ. 1 ) GO TO 345

C  top BC for first active layer 1: no downward diffuse radiation

      NS = 1
      N = ACTIVE_LAYERS(NS)
      DO I = 1, NSTREAMS
        DO EP = 1, NSTREAMS
          EM = EP + NSTREAMS
          BANDTELMAT2(BTELMAT_ROWMASK(I,EP),EP)  =
     &                   XPOS(I,EP,N)
          BANDTELMAT2(BTELMAT_ROWMASK(I,EM),EM)  =
     &                   XNEG(I,EP,N)*T_DELT_EIGEN(EP,N)
        ENDDO
      ENDDO

C  intermediate layer boundaries

      C0 = - NSTREAMS
      DO NS = 2, NLAYERS_TEL
        N = ACTIVE_LAYERS(NS)
        N1 = N - 1
        C0   = C0 + NSTREAMS_2
        CE_OFFSET = C0 - NSTREAMS
        DO I = 1, NSTREAMS_2
          CM = C0 + I
          DO EP = 1, NSTREAMS
            CEP = CE_OFFSET + EP
            CEM = CEP + NSTREAMS
            CEP1 = CEP + NSTREAMS_2
            CEM1 = CEM + NSTREAMS_2
            BANDTELMAT2(BTELMAT_ROWMASK(CM,CEP),CEP)   =
     &                   T_DELT_EIGEN(EP,N1)*XPOS(I,EP,N1)
            BANDTELMAT2(BTELMAT_ROWMASK(CM,CEM),CEM)   =
     &                   XNEG(I,EP,N1)
            BANDTELMAT2(BTELMAT_ROWMASK(CM,CEP1),CEP1) =
     &                   -XPOS(I,EP,N)
            BANDTELMAT2(BTELMAT_ROWMASK(CM,CEM1),CEM1) =
     &                   -T_DELT_EIGEN(EP,N)*XNEG(I,EP,N)
          ENDDO
        ENDDO
      ENDDO

C  bottom BC (No albedo additions). Lowest active layer

      N = ACTIVE_LAYERS(NLAYERS_TEL)
      C0 = C0 + NSTREAMS_2
      CE_OFFSET = C0 - NSTREAMS

      IF ( DO_INCLUDE_SURFACE ) THEN
C  Place holder
      ELSE
        DO I = 1, NSTREAMS
          CP = C0 + I
          I1 = I + NSTREAMS
          DO AA = 1, NSTREAMS
            CEP = CE_OFFSET + AA
            CEM = CEP + NSTREAMS
            BANDTELMAT2(BTELMAT_ROWMASK(CP,CEP),CEP) =
     &                  T_DELT_EIGEN(AA,N) * XPOS(I1,AA,N)
            BANDTELMAT2(BTELMAT_ROWMASK(CP,CEM),CEM) =
     &                  XNEG(I1,AA,N)
          ENDDO
        ENDDO
      ENDIF

C  normal completion

      RETURN

C  special case. Only 1 active layer

345   CONTINUE

C  Set up BVP matrix for the active layer
C  ======================================

C  initialize using the SMAT2 matrix

      DO I = 1, NSTREAMS_2
        DO J = 1, NSTREAMS_2
          SMAT2(I,J) = ZERO
        ENDDO
      ENDDO

C  top BC for layer: no downward diffuse radiation
C  -----------------------------------------------

      N = ACTIVE_LAYERS(1)
      DO I = 1, NSTREAMS
        DO EP = 1, NSTREAMS
          EM = EP + NSTREAMS
          SMAT2(I,EP) = XPOS(I,EP,N)
          SMAT2(I,EM) = XNEG(I,EP,N)*T_DELT_EIGEN(EP,N)
        ENDDO
      ENDDO

C  bottom BC (No albedo additions)
C  -------------------------------

      IF ( DO_INCLUDE_SURFACE ) THEN
C  Place holder
      ELSE
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          CP = NSTREAMS + I
          DO EP = 1, NSTREAMS
            CEP = EP
            CEM = CEP + NSTREAMS
            XPNET = XPOS(I1,EP,N)
            XMNET = XNEG(I1,EP,N)
            SMAT2(CP,CEP) = T_DELT_EIGEN(EP,N) * XPNET
            SMAT2(CP,CEM) = XMNET
          ENDDO
        ENDDO
      ENDIF

C  normal return and finish

      RETURN
      END

C

      SUBROUTINE BVPTEL_MATRIX_SVD ( STATUS )

C TELESCOPED boundary value problem SVD decomposition.

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include file of input variables
C  include file of bookkeeping inputs

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file of solution variables (I/O to this module)

      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  subroutine arguments
C  --------------------

C  status flag

      INTEGER          STATUS

C  local variables
C  ---------------

      CHARACTER*(70)   MAIL, TRACE
      INTEGER          INFO
      CHARACTER*3      CI

C  Intialize status

      STATUS = LIDORT_SUCCESS

C  SVD the BVPTEL matrix: With compression (multilayers)
C  ----------------------------------------------------

      IF ( NLAYERS_TEL .GT. 1 ) THEN

C  LAPACK LU-decomposition for band matrix

        CALL DGBTRF
     &     ( N_BVTELMATRIX_SIZE, N_BVTELMATRIX_SIZE,
     &       N_BVTELMATRIX_SUBDIAG, N_BVTELMATRIX_SUPDIAG,
     &       BANDTELMAT2, MAXBANDTOTAL, IPIVOTTEL, INFO )

C  (Error tracing)

        IF ( INFO .GT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MAIL  = 'Singular matrix, u(i,i)=0, for i = '//CI
          TRACE = 'DGBTRF call in BVPTEL_MATRIX_SVD'
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ELSE IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MAIL  = 'argument i illegal value, for i = '//CI
          TRACE = 'DGBTRF call in BVPTEL_MATRIX_SVD'
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ENDIF

C  SVD the BVP matrix: No compression, Single Layer only
C  -----------------------------------------------------

      ELSE IF ( NLAYERS_TEL .EQ. 1 ) THEN

C  LAPACK LU-decomposition for single layer matrix

        CALL DGETRF
     &     (  NSTREAMS_2, NSTREAMS_2,
     &        SMAT2, MAXSTREAMS_2, SIPIVOT, INFO )

C  (Error tracing)

        IF ( INFO .GT. 0 ) THEN
          MAIL='Failure: BVP LU decomposition (telescoped) (DGETRF)'
          CALL LIDORT_LAPACK_ERROR ( INFO, 1, MAIL, STATUS )
          IF ( STATUS .EQ. LIDORT_SERIOUS ) RETURN
        ENDIF

      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE BVPTEL_SOLUTION_MASTER
     I       ( FOURIER_COMPONENT, IBEAM,
     O         STATUS )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  input arguments
C  ---------------

C  Fourier component and beam number

      INTEGER          FOURIER_COMPONENT, IBEAM

C  output status
C  -------------

      INTEGER          STATUS

C  Local variables
C  ---------------

      CHARACTER*(70)   MAIL, TRACE
      INTEGER          STATUS_SUB

C  This is suitable for Lambertian surfaces only.

C  --set up Column for solution vector (the "B" as in AX=B)

      CALL BVPTEL_COLUMN_SETUP  ( IBEAM, FOURIER_COMPONENT )

C  --Solve the boundary problem for this Fourier component (back substitution)

      CALL BVPTEL_BACKSUB
     I   ( FOURIER_COMPONENT, IBEAM, 
     O     STATUS_SUB )

C  error tracing

      IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
        MAIL  = 'Error return from module BVP_BACKSUB'
        TRACE = 'Call #1 in BVPTEL_BEAMSOLUTION_MASTER'
        STATUS = LIDORT_SERIOUS
        CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
      ENDIF

C  return

      RETURN
      END

C

      SUBROUTINE BVPTEL_COLUMN_SETUP ( IBEAM, FOURIER )

C  Sets up the telescoped boundary value problem, RHS vector
C    Standard case: Fourier > 0. No surface reflection term
C         Suitable for Lambertian surfaces.

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include file of input variables
C  include file of bookkeeping inputs

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of main model variable (I/O to this module)

      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  subroutine arguments
C  --------------------

C  beam and fourier index

      INTEGER          IBEAM, FOURIER

C  local variables
C  ---------------

      INTEGER          I, I1, N, N1, NS, C0, CM   

C  Go to special case for only 1 active layer

      IF ( NLAYERS_TEL .EQ. 1 ) GO TO 345

C  zero column vector

      DO I = 1, N_BVTELMATRIX_SIZE
        COLTEL2(I,1) = ZERO
      ENDDO

C  Upper boundary for first active layer: no downward diffuse radiation

      NS = 1
      N = ACTIVE_LAYERS(NS)
      DO I = 1, NSTREAMS
        COLTEL2(I,IBEAM)  =  - WUPPER(I,N)
      ENDDO

C  intermediate layer boundaries

      C0 = - NSTREAMS
      DO NS = 2, NLAYERS_TEL
        N = ACTIVE_LAYERS(NS)
        N1 = N - 1
        C0  = C0 + NSTREAMS_2
        DO I = 1, NSTREAMS_2
          CM = C0 + I
          COLTEL2(CM,IBEAM) = WUPPER(I,N) - WLOWER(I,N1)
        ENDDO
      ENDDO

C  Lower boundary for last active layer NLAYERS_TEL:
C  No albedo, as this is FOURIER > 0

      C0 = C0 + NSTREAMS_2
      NS = NLAYERS_TEL
      N = ACTIVE_LAYERS(NS)
      DO I = 1, NSTREAMS
        I1 = I + NSTREAMS
        CM = C0 + I
        COLTEL2(CM,IBEAM) = - WLOWER(I1,N)
      ENDDO

C  Normal return

      RETURN

C  Continuity point for single layer setup

 345  CONTINUE

C  Set up BVP column for the Single active layer
C  =============================================

C  initialize using the SCOL2 matrix

      DO I = 1, NSTREAMS_2
        SCOL2(I,IBEAM) = ZERO
      ENDDO

C  active layer for telescoped BVP

      N = ACTIVE_LAYERS(1)

C  Upper boundary for layer (downwelling only)
C  -------------------------------------------

      DO I = 1, NSTREAMS
        SCOL2(I,IBEAM)   = - WUPPER(I,N)
      ENDDO

C  lower boundary for layer (upwelling only)
C  -----------------------------------------

      C0 = NSTREAMS
      DO I = 1, NSTREAMS
        I1 = I + NSTREAMS
        CM = C0 + I
        SCOL2(CM,IBEAM) = - WLOWER(I1,N)
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE BVPTEL_BACKSUB
     I         ( FOURIER, IBEAM,
     O           STATUS )

C  Solves the telescoped boundary value problem.
C    Standard case: Fourier > 0. No surface reflection term
C         Suitable for Lambertian surfaces.

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include file of input variables
C  include file of bookkeeping inputs

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file of setup variables

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  include files of main model variable (I/O to this module)

      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  subroutine arguments
C  --------------------

C  Fourier component and beam index

      INTEGER          FOURIER, IBEAM

C  status flag

      INTEGER          STATUS

C  local variables
C  ---------------

      CHARACTER*70     MAIL, TRACE
      INTEGER          I, I1, N, N1, NAF, NAL, NS
      INTEGER          INFO, K, K1, C0, CM, CP    
      DOUBLE PRECISION SHOM, HOM1, HOM2
      CHARACTER*3      CI

C  start code
C  ----------

C  Intialize status

      STATUS = LIDORT_SUCCESS

C  Back-substitution for multi-layer BVP TEL
C  =========================================

      IF ( NLAYERS_TEL .GT. 1 ) THEN

C  LAPACK substitution using RHS column vector COL2

        CALL DGBTRS
     &     ( 'n', N_BVTELMATRIX_SIZE,
     &        N_BVTELMATRIX_SUBDIAG, N_BVTELMATRIX_SUPDIAG, IBEAM,
     &        BANDTELMAT2, MAXBANDTOTAL, IPIVOTTEL,
     &        COLTEL2, MAXTOTAL, INFO )

C  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MAIL  = 'argument i illegal value, for i = '//CI
          TRACE = 'DGBTRS call in BVPTEL_BACKSUB (telescoping)'
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ENDIF

C  set integration constants for active layers

        C0 = -NSTREAMS_2
        DO NS = 1, NLAYERS_TEL
          N = ACTIVE_LAYERS(NS)
          C0 = C0 + NSTREAMS_2
          DO I = 1, NSTREAMS
            CM = C0 + I
            CP = CM + NSTREAMS
            LCON(I,N) = COLTEL2(CM,IBEAM)
            MCON(I,N) = COLTEL2(CP,IBEAM)
          ENDDO
        ENDDO

C  Solve the boundary problem: Single Layer only
C  =============================================

      ELSE IF ( NLAYERS_TEL .EQ. 1 ) THEN

C  LAPACK substitution (DGETRS) using RHS column vector SCOL2

        CALL DGETRS
     &     ( 'N', NSTREAMS_2, IBEAM, 
     &        SMAT2, MAXSTREAMS_2, SIPIVOT,
     &        SCOL2, MAXSTREAMS_2, INFO )

C  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          MAIL='BVP back-substitution (telescoping) (DGETRS)'
          CALL LIDORT_LAPACK_ERROR ( INFO, 1, MAIL, STATUS )
          IF ( STATUS .EQ. LIDORT_SERIOUS ) RETURN
        ENDIF

C  Set integration constants LCON and MCON for active layer

        N = ACTIVE_LAYERS(1)
        DO K = 1, NSTREAMS
          K1 = K + NSTREAMS 
          LCON(K,N) = SCOL2(K, IBEAM)
          MCON(K,N) = SCOL2(K1,IBEAM)
        ENDDO

C  end clause for backsubstitution

      ENDIF

C  Associated quantities for active layers
C  ---------------------------------------

      DO NS = 1, NLAYERS_TEL
        N = ACTIVE_LAYERS(NS)
        DO I = 1, NSTREAMS_2
          DO K = 1, NSTREAMS
            LCON_XVEC(I,K,N) = LCON(K,N) * XPOS(I,K,N)
            MCON_XVEC(I,K,N) = MCON(K,N) * XNEG(I,K,N)
          ENDDO
        ENDDO
      ENDDO

C  Set integration constants for non-active layers
C  ===============================================

C  Transmittance layers ABOVE active layer(s)
C  -----------------------------------------

C   -- LCON values are zero (no downwelling radiation)
C   -- MCON values propagated upwards from top of first active layer

C  layer immediately above first active layer

      NAF = ACTIVE_LAYERS(1)
      IF ( NAF .GT. 1 ) THEN
        N1 = NAF - 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          SHOM = ZERO
          DO K = 1, NSTREAMS
            HOM1 = LCON_XVEC(I1,K,NAF)
            HOM2 = MCON_XVEC(I1,K,NAF) * T_DELT_EIGEN(K,NAF)
            SHOM = SHOM + HOM1 + HOM2
          ENDDO
          MCON(I,N1) = WUPPER(I1,NAF) + SHOM
          LCON(I,N1) = ZERO
        ENDDO
      ENDIF

C  other layers to top: propagate usng discrete ordinate tranmsittance

      DO N = NAF - 2, 1, -1
        N1 = N + 1
        DO I = 1, NSTREAMS
          LCON(I,N) = ZERO
          MCON(I,N) = T_DELT_DISORDS(I,N1) * MCON(I,N1)
        ENDDO
      ENDDO     

C  Transmittance layers below active layer
C  ---------------------------------------

C   -- MCON values are zero (no upwelling radiation)
C   -- LCON values propagated downwards from bottom of last active layer

C  layer immediately below last active layer

      NAL = ACTIVE_LAYERS(NLAYERS_TEL)
      IF ( NAL .LT. NLAYERS ) THEN
        N1 = NAL + 1
        DO I = 1, NSTREAMS
          SHOM = ZERO
          DO K = 1, NSTREAMS
            HOM1 = LCON_XVEC(I,K,NAL) * T_DELT_EIGEN(K,NAL)
            HOM2 = MCON_XVEC(I,K,NAL)
            SHOM = SHOM + HOM1 + HOM2
          ENDDO
          LCON(I,N1) = WLOWER(I,NAL) + SHOM
          MCON(I,N1) = ZERO
        ENDDO
      ENDIF

C  other layers to bottom

      DO N = NAL + 2, NLAYERS
        N1 = N - 1
        DO I = 1, NSTREAMS
          MCON(I,N) = ZERO
          LCON(I,N) = T_DELT_DISORDS(I,N1) * LCON(I,N1)
        ENDDO
      ENDDO

C  Associated quantities for non-active layers
C  -------------------------------------------

C  not efficient coding. Solutions are trivial.

      DO N = 1, NLAYERS
        IF ( N .LT. NAF .OR. N.GT.NAL ) THEN
          DO I = 1, NSTREAMS_2
            DO K = 1, NSTREAMS
              LCON_XVEC(I,K,N) = LCON(K,N) * XPOS(I,K,N)
              MCON_XVEC(I,K,N) = MCON(K,N) * XNEG(I,K,N)
            ENDDO
          ENDDO
        ENDIF
      ENDDO

C  debug

c      k = 97
c      if ( do_fdtest ) k = 98
c      if ( fourier. EQ.3 .and.ibeam.eq.1 ) then
c        do n = 1, nlayers
c          write(k,'(i3,1p2e18.10)')n,LCON(3,n), MCON(3,n)
c        enddo
c      endif

c      K= 98
c      i = 5
c      if ( fourier.eq.3.and.ibeam.eq.1) then
c        DO N = 1, NLAYERS
c          DO I = 1, NSTREAMS
c          write(K,'(3i4,1p4e17.9)')IBEAM,N,I, LCON(I,N), MCON(I,N)
c        ENDDO
c        ENDDO
c      ENDIF

C  Finish

      RETURN
      END

