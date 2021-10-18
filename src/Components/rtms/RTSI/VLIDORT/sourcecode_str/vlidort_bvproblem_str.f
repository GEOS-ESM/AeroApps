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
C #            BVP_MATRIXSETUP_MASTER      (master)             #
C #             BVP_LAMBERTIAN_SURFACE_HOM                      #
C #             BVP_MATRIX_INIT                                 #
C #             BVP_MATRIX_SETUP                                #
C #             BVP_MATRIX_SVD                                  #
C #                                                             #
C #            BVP_SOLUTION_MASTER      (master)                #
C #             BVP_LAMBERTIAN_SURFACE_BEAM                     #
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
C #            BMAT_ROWMASK    (integer function)               #
C #            BTELMAT_ROWMASK (integer function)               #
C #                                                             #
C ###############################################################

      SUBROUTINE BVP_MATRIXSETUP_MASTER
     I     ( DO_INCLUDE_SURFACE,
     I       FOURIER_COMPONENT,
     I       SURFACE_FACTOR,
     I       BVP_SURFACE_SETUP_HOM,
     O       STATUS )

C  Include files
C  -------------

C  Additional sums for the final albedo-reflecting layer

C  include file of dimensiopns and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of solution stuff (i/o to this module)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'

C  Arguments
C  ---------

C  control input

      LOGICAL          DO_INCLUDE_SURFACE
      INTEGER          FOURIER_COMPONENT
      DOUBLE PRECISION SURFACE_FACTOR

C  Surface reflectance function

      EXTERNAL         BVP_SURFACE_SETUP_HOM

C  output status

      INTEGER          STATUS

C  local variables
C  ---------------

      CHARACTER*(70)   MAIL, TRACE
      INTEGER          STATUS_SUB

C  Additional setups for the albedo layer

      CALL BVP_SURFACE_SETUP_HOM
     I  ( DO_INCLUDE_SURFACE,
     I    FOURIER_COMPONENT,
     I    SURFACE_FACTOR )

C  initialize compression matrix (Do this for every Fourier component)

      CALL BVP_MATRIX_INIT

C  set up boundary values matrix in compressed form (the "A" as in AX=B)

      CALL BVP_MATRIX_SETUP
     I     ( DO_INCLUDE_SURFACE, FOURIER_COMPONENT )

C  SVD decomposition of compressed boundary values matrix

      CALL BVP_MATRIX_SVD ( STATUS_SUB )

C  error tracing

      IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
        MAIL = 'Error return from module BVP_MATRIX_SVD'
        TRACE= 'Call in BVP_MATRIXSETUP_MASTER'
        STATUS = VLIDORT_SERIOUS
        CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
        RETURN
      ENDIF

C  finish

      RETURN
      END

C

      SUBROUTINE BVP_LAMBERTIAN_SURFACE_HOM
     I     ( DO_INCLUDE_SURFACE,
     I       FOURIER_COMPONENT,
     I       SURFACE_FACTOR )

C  This is the Lambertian surface routine
C     Reflected homogeneous solutions

C  Include files
C  -------------

C  Additional sums for the final albedo-reflecting layer

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of solution stuff (i/o to this module)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'

C  control input
C  -------------

      LOGICAL          DO_INCLUDE_SURFACE
      INTEGER          FOURIER_COMPONENT
      DOUBLE PRECISION SURFACE_FACTOR

C  local variables
C  ---------------

      DOUBLE PRECISION H_1, H_2, H_1_CR, H_2_CR, H_1_CI, H_2_CI
      DOUBLE PRECISION AXJ, KMULT
      INTEGER          I, J, O1, K, KO1, K0, K1, K2 

C  Initialization
C  ==============

C  Zero total reflected contributions

      DO I = 1, NSTREAMS
        DO K = 1, NSTKS_NSTRMS
          DO O1 = 1, NSTOKES
            R2_HOMP(I,O1,K) = ZERO
            R2_HOMM(I,O1,K) = ZERO
          ENDDO
        ENDDO
      ENDDO

C  Return with Zeroed values if albedo flag not set

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

C  Return if Fourier component not zero

      IF ( FOURIER_COMPONENT .GT. 0 ) RETURN

C  For Lambertian reflectance, all streams are the same
C  ----------------------------------------------------

C  Integrate Downward streams of particular solutions

      KMULT = SURFACE_FACTOR * LAMBERTIAN_ALBEDO

C  Homogeneous real solutions

      DO K = 1, K_REAL(NLAYERS)
        DO I = 1, NSTREAMS
          H_1 = ZERO
          H_2 = ZERO
          DO J = 1, NSTREAMS
            AXJ = QUAD_STRMWTS(J)
            H_1 = H_1 + AXJ*SOLA_XPOS(J,1,K,NLAYERS)
            H_2 = H_2 + AXJ*SOLB_XNEG(J,1,K,NLAYERS)
          ENDDO
          R2_HOMP(I,1,K) = H_1 * KMULT
          R2_HOMM(I,1,K) = H_2 * KMULT
        ENDDO
      ENDDO

C  Homogeneous complex solutions

      KO1 = K_REAL(NLAYERS) + 1
      DO K = 1, K_COMPLEX(NLAYERS)
        K0 = 2*K - 2
        K1 = KO1 + K0
        K2 = K1  + 1
        DO I = 1, NSTREAMS
          H_1_CR = ZERO
          H_2_CR = ZERO
          H_1_CI = ZERO
          H_2_CI = ZERO
          DO J = 1, NSTREAMS
            AXJ = QUAD_STRMWTS(J)
            H_1_CR = H_1_CR + AXJ*SOLA_XPOS(J,1,K1,NLAYERS)
            H_2_CR = H_2_CR + AXJ*SOLB_XNEG(J,1,K1,NLAYERS)
            H_1_CI = H_1_CI + AXJ*SOLA_XPOS(J,1,K2,NLAYERS)
            H_2_CI = H_2_CI + AXJ*SOLB_XNEG(J,1,K2,NLAYERS)
          ENDDO
          R2_HOMP(I,1,K1) = H_1_CR * KMULT
          R2_HOMM(I,1,K1) = H_2_CR * KMULT
          R2_HOMP(I,1,K2) = H_1_CI * KMULT
          R2_HOMM(I,1,K2) = H_2_CI * KMULT
        ENDDO
      ENDDO

C  debug

c      DO I = 1, NSTREAMS
c        DO O1 = 1, NSTOKES
c          DO K = 1, NSTKS_NSTRMS
c           write(*,'(3i3,1p2e14.6)')
c     &             I,O1,K,R2_HOMP(I,O1,K),R2_HOMM(I,O1,K)
c          ENDDO
c        ENDDO
c      ENDDO
c      pause'r2homp'
      
C  Finish

      RETURN
      END

C

      SUBROUTINE BVP_MATRIX_INIT

C  Initialisee the compressed matrix

C  include files
C  -------------

C  include file of dimensiopns and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of main model variables (local to this module)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'

C  Local variables
C  ---------------

C  compression function

      INTEGER  BMAT_ROWMASK
      EXTERNAL BMAT_ROWMASK

C  help variables

      INTEGER  I, J, N3, JS, JF, IS, LF, L, I1

C  special case

      IF ( NLAYERS .EQ. 1 ) THEN
        DO I = 1, NTOTAL
          DO J = 1, NTOTAL
            SMAT2(I,J)     = ZERO
          ENDDO
        ENDDO
        RETURN
      ENDIF

c  compression row indices, NMIN, NMAX and KALL
C    Results are now stored in commons

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

C  Former code to assign Array - now replaced by external function
C    Jukaa Kujanpaa, FMI, August 2005
c      DO I = 1, NTOTAL
c        DO J = 1, NTOTAL
c          IF ( (I.GE.NMIN(J)) .AND. (I.LE.NMAX(J)) ) THEN
c            BMAT_ROWMASK(I,J) = KALL + I - J
c          ENDIF
c        ENDDO
c      ENDDO

C  compression matrix zeroing

      N3 = NSTKS_NSTRMS_2 + NSTKS_NSTRMS
      LF = NLAYERS - 2

C  upper band top

      JS = NSTKS_NSTRMS_2 + 1
      JF = N3 - 1
      DO I = 1, NSTKS_NSTRMS
        DO J = JS, JF + I
          BANDMAT2(BMAT_ROWMASK(I,J),J) = ZERO
        ENDDO
      ENDDO

C  upper band

      DO L = 1, LF
        IS = L*NSTKS_NSTRMS_2 - NSTKS_NSTRMS + 1
        JS = IS + N3
        JF = JS - 1
        DO I = 1, NSTKS_NSTRMS_2-1
          I1 = I + IS
          DO J = JS, JF + I
            BANDMAT2(BMAT_ROWMASK(I1,J),J) = ZERO
          ENDDO
        ENDDO
      ENDDO

C  lower band

      DO L = 1, LF
        IS = L*NSTKS_NSTRMS_2 + NSTKS_NSTRMS
        JS = IS - N3 + 1
        JF = IS - NSTKS_NSTRMS
        DO I = 1, NSTREAMS_2-1
          I1 = I + IS
          DO J = JS + I, JF
            BANDMAT2(BMAT_ROWMASK(I1,J),J) = ZERO
          ENDDO
        ENDDO
      ENDDO

C  lower band bottom

      JS = LF * NSTKS_NSTRMS_2 + 1
      IS = JS + N3 - 1
      JF = IS - NSTKS_NSTRMS
      DO I = 1, NSTKS_NSTRMS
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

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of main model variables (local to this module)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'

C  input arguments
C  ---------------

      LOGICAL          DO_INCLUDE_SURFACE
      INTEGER          FOURIER_COMPONENT

C  compression function
C  --------------------

      INTEGER  BMAT_ROWMASK
      EXTERNAL BMAT_ROWMASK

C  Local variables
C  ---------------

      INTEGER          I, I1, N, N1, C0, CE_OFFSET, O1, M
      INTEGER          EP, EM, EPC, EPS, EMS
      INTEGER          CP, CM, CEP, CEM, CEM1, CEP1, CEPS, CEMS
      INTEGER          IR, IROW, KO11, KO1, K0, K1, K2
      DOUBLE PRECISION XPNET, XMNET
      DOUBLE PRECISION X1R, X2R, X1I, X2I, X1R_H, X1I_H

C  Initialization

      EP = 0

C  Fourier component

      M = FOURIER_COMPONENT

C  If Nlayers = 1, go to special case

      IF ( NLAYERS .EQ. 1 ) GO TO 345

C  Initialise - makes no difference
c      DO EP = 1, NTOTAL
c        DO EM = 1, MAXBANDTOTAL
c          BANDMAT2(EM,EP) = ZERO
c        ENDDO
c      ENDDO

C  top BC for layer 1: no downward diffuse radiation
C  -------------------------------------------------

      N = 1
      KO1 = K_REAL(N) + 1
      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        DO O1 = 1, NSTOKES
          IROW = IR + O1

C  real solutions

          DO EP = 1, K_REAL(N)
            EM = EP + NSTKS_NSTRMS
            BANDMAT2(BMAT_ROWMASK(IROW,EP),EP)  =
     &                   SOLA_XPOS(I,O1,EP,N)
            BANDMAT2(BMAT_ROWMASK(IROW,EM),EM)  =
     &                   SOLB_XNEG(I,O1,EP,N)*T_DELT_EIGEN(EP,N)
          ENDDO

C  complex solutions (REWORKED)

          DO EPC = 1, K_COMPLEX(N)
            K0 = 2*EPC - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            EP  = K_REAL(N) + EPC
            EPS = K_REAL(N) + K_COMPLEX(N) + EPC
            EM  = EP  + NSTKS_NSTRMS
            EMS = EPS + NSTKS_NSTRMS
            X1R = SOLA_XPOS(I,O1,K1,N)
            X1I = SOLA_XPOS(I,O1,K2,N)
            X2R = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K1,N)
     &           -SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K2,N)
            X2I = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K2,N)
     &           +SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K1,N)
            BANDMAT2(BMAT_ROWMASK(IROW,EP),EP)   =   X1R
            BANDMAT2(BMAT_ROWMASK(IROW,EPS),EPS) = - X1I
            BANDMAT2(BMAT_ROWMASK(IROW,EM),EM)   =   X2R
            BANDMAT2(BMAT_ROWMASK(IROW,EMS),EMS) = - X2I
          ENDDO

        ENDDO
      ENDDO

C  intermediate layer boundaries
C  -----------------------------

      C0 = - NSTKS_NSTRMS
      DO N = 2, NLAYERS        
        N1 = N - 1
        KO1  = K_REAL(N) + 1
        KO11 = K_REAL(N1) + 1
        C0   = C0 + NSTKS_NSTRMS_2
        CE_OFFSET = C0 - NSTKS_NSTRMS
        DO I = 1, NSTREAMS_2
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW

C  real solutions for layer above the boundary

            DO EP = 1, K_REAL(N1)
              CEP = CE_OFFSET + EP
              CEM = CEP + NSTKS_NSTRMS
              BANDMAT2(BMAT_ROWMASK(CM,CEP),CEP)   =
     &                 T_DELT_EIGEN(EP,N1)*SOLA_XPOS(I,O1,EP,N1)
              BANDMAT2(BMAT_ROWMASK(CM,CEM),CEM)   =
     &                   SOLB_XNEG(I,O1,EP,N1)
            ENDDO

C  real solutions for layer below the boundary
C   ( Note the change of sign !!! )

            DO EP = 1, K_REAL(N)
              CEP = CE_OFFSET + EP
              CEM = CEP + NSTKS_NSTRMS
              CEP1 = CEP + NSTKS_NSTRMS_2
              CEM1 = CEM + NSTKS_NSTRMS_2
              BANDMAT2(BMAT_ROWMASK(CM,CEP1),CEP1) =
     &                   - SOLA_XPOS(I,O1,EP,N)
              BANDMAT2(BMAT_ROWMASK(CM,CEM1),CEM1) =
     &                 - T_DELT_EIGEN(EP,N)*SOLB_XNEG(I,O1,EP,N)
            ENDDO

C  complex solutions for layer above boundary

            DO EPC = 1, K_COMPLEX(N1)
              K0 = 2*EPC - 2
              K1 = KO11 + K0
              K2 = K1  + 1
              EP  = K_REAL(N1) + EPC
              EPS = K_REAL(N1) + K_COMPLEX(N1) + EPC
              CEP  = CE_OFFSET + EP
              CEM  = CEP + NSTKS_NSTRMS
              CEPS = CE_OFFSET + EPS
              CEMS = CEPS + NSTKS_NSTRMS
              X1R = SOLA_XPOS(I,O1,K1,N1)*T_DELT_EIGEN(K1,N1)
     &             -SOLA_XPOS(I,O1,K2,N1)*T_DELT_EIGEN(K2,N1)
              X1I = SOLA_XPOS(I,O1,K1,N1)*T_DELT_EIGEN(K2,N1)
     &             +SOLA_XPOS(I,O1,K2,N1)*T_DELT_EIGEN(K1,N1)
              X2R = SOLB_XNEG(I,O1,K1,N1)
              X2I = SOLB_XNEG(I,O1,K2,N1)
              BANDMAT2(BMAT_ROWMASK(CM,CEP),CEP)   =   X1R
              BANDMAT2(BMAT_ROWMASK(CM,CEPS),CEPS) = - X1I
              BANDMAT2(BMAT_ROWMASK(CM,CEM),CEM)   =   X2R
              BANDMAT2(BMAT_ROWMASK(CM,CEMS),CEMS) = - X2I
            ENDDO

C  complex solutions for layer below boundary
C   ( Note the change of sign !!! )

            DO EPC = 1, K_COMPLEX(N)
              K0 = 2*EPC - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              EP  = K_REAL(N) + EPC
              EPS = K_REAL(N) + K_COMPLEX(N) + EPC
              CEP  = CE_OFFSET + EP + NSTKS_NSTRMS_2
              CEM  = CEP + NSTKS_NSTRMS
              CEPS = CE_OFFSET + EPS + NSTKS_NSTRMS_2
              CEMS = CEPS + NSTKS_NSTRMS
              X1R = SOLA_XPOS(I,O1,K1,N)
              X1I = SOLA_XPOS(I,O1,K2,N)
              X2R = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K1,N)
     &             -SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K2,N)
              X2I = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K2,N)
     &             +SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K1,N)
              BANDMAT2(BMAT_ROWMASK(CM,CEP),CEP)   = - X1R
              BANDMAT2(BMAT_ROWMASK(CM,CEPS),CEPS) =   X1I
              BANDMAT2(BMAT_ROWMASK(CM,CEM),CEM)   = - X2R
              BANDMAT2(BMAT_ROWMASK(CM,CEMS),CEMS) =   X2I
            ENDDO

          ENDDO
        ENDDO
      ENDDO

C  bottom BC (with albedo additions if flagged)

      N   = NLAYERS
      KO1 = K_REAL(N) + 1
      C0  = C0 + NSTKS_NSTRMS_2
      CE_OFFSET = C0 - NSTKS_NSTRMS

C  With albedo
C  -----------

      IF ( DO_INCLUDE_SURFACE ) THEN
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CP = C0 + IROW

C  real solutions

            DO EP = 1, K_REAL(N)
              CEP = CE_OFFSET + EP
              CEM = CEP + NSTKS_NSTRMS
              XPNET = SOLA_XPOS(I1,O1,EP,N) - R2_HOMP(I,O1,EP)
              XMNET = SOLB_XNEG(I1,O1,EP,N) - R2_HOMM(I,O1,EP)
              BANDMAT2(BMAT_ROWMASK(CP,CEP),CEP) =
     &                   T_DELT_EIGEN(EP,N) * XPNET
              BANDMAT2(BMAT_ROWMASK(CP,CEM),CEM) =
     &                   XMNET
            ENDDO

C  Complex solutions

            DO EPC = 1, K_COMPLEX(N)
              K0 = 2*EPC - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              EP  = K_REAL(N) + EPC
              EPS = K_REAL(N) + K_COMPLEX(N) + EPC
              CEP = CE_OFFSET + EP
              CEM = CEP + NSTKS_NSTRMS
              CEPS = CE_OFFSET + EPS
              CEMS = CEPS + NSTKS_NSTRMS
              X1R_H = SOLA_XPOS(I1,O1,K1,N) - R2_HOMP(I,O1,K1)
              X1I_H = SOLA_XPOS(I1,O1,K2,N) - R2_HOMP(I,O1,K2)
              X2R = SOLB_XNEG(I1,O1,K1,N) - R2_HOMM(I,O1,K1)
              X2I = SOLB_XNEG(I1,O1,K2,N) - R2_HOMM(I,O1,K2)
              X1R = X1R_H * T_DELT_EIGEN(K1,N) -
     &              X1I_H * T_DELT_EIGEN(K2,N)
              X1I = X1R_H * T_DELT_EIGEN(K2,N) +
     &              X1I_H * T_DELT_EIGEN(K1,N)
              BANDMAT2(BMAT_ROWMASK(CP,CEP),CEP)   =   X1R
              BANDMAT2(BMAT_ROWMASK(CP,CEPS),CEPS) = - X1I
              BANDMAT2(BMAT_ROWMASK(CP,CEM),CEM)   =   X2R
              BANDMAT2(BMAT_ROWMASK(CP,CEMS),CEMS) = - X2I

c         iunit = 20
c        IF (do_fdtest)iunit=19
c         write(iunit,4567)I,O1,EPC,
c     &         R2_HOMP(I,O1,EPC)*T_DELT_EIGEN(EPC,N)+
c     &         R2_HOMP(I,O1,EPC)*T_DELT_EIGEN(EPC,N),
c     &                   R2_HOMM(I,O1,EPC)
c 4567  format   (3i4,1p2e20.10)  

            ENDDO

          ENDDO
        ENDDO

C  No albedo

      ELSE

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CP = C0 + IROW

C  real solutions

            DO EP = 1, K_REAL(N)
              CEP = CE_OFFSET + EP
              CEM = CEP + NSTKS_NSTRMS
              XPNET = SOLA_XPOS(I1,O1,EP,N)
              XMNET = SOLB_XNEG(I1,O1,EP,N)
              BANDMAT2(BMAT_ROWMASK(CP,CEP),CEP) =
     &                   T_DELT_EIGEN(EP,N) * XPNET
              BANDMAT2(BMAT_ROWMASK(CP,CEM),CEM) =
     &                   XMNET
            ENDDO

C  Complex solutions
C    ----- Deep bug found for X1I, 17 December 2005.

            DO EPC = 1, K_COMPLEX(N)
              K0 = 2*EPC - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              EP  = K_REAL(N) + EPC
              EPS = K_REAL(N) + K_COMPLEX(N) + EPC
              CEP  = CE_OFFSET + EP
              CEM  = CEP + NSTKS_NSTRMS
              CEPS = CE_OFFSET + EPS
              CEMS = CEPS + NSTKS_NSTRMS
              X1R = SOLA_XPOS(I1,O1,K1,N) * T_DELT_EIGEN(K1,N)
     &             -SOLA_XPOS(I1,O1,K2,N) * T_DELT_EIGEN(K2,N)
              X1I = SOLA_XPOS(I1,O1,K1,N) * T_DELT_EIGEN(K2,N)
     &             +SOLA_XPOS(I1,O1,K2,N) * T_DELT_EIGEN(K1,N)
              X2R = SOLB_XNEG(I1,O1,K1,N)
              X2I = SOLB_XNEG(I1,O1,K2,N)
              BANDMAT2(BMAT_ROWMASK(CP,CEP),CEP)   =   X1R
              BANDMAT2(BMAT_ROWMASK(CP,CEPS),CEPS) = - X1I
              BANDMAT2(BMAT_ROWMASK(CP,CEM),CEM)   =   X2R
              BANDMAT2(BMAT_ROWMASK(CP,CEMS),CEMS) = - X2I
           ENDDO

          ENDDO
        ENDDO

      ENDIF

C  normal completion

      RETURN

C  special case for 1 layer only
C  =============================

345   CONTINUE

C  top BC for layer 1: no downward diffuse radiation
C  -------------------------------------------------

      N = 1
      KO1 = K_REAL(N) + 1
      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        DO O1 = 1, NSTOKES
          IROW = IR + O1

C  real solutions

          DO EP = 1, K_REAL(N)
            EM = EP + NSTKS_NSTRMS
            SMAT2(IROW,EP) = SOLA_XPOS(I,O1,EP,N)
            SMAT2(IROW,EM) = SOLB_XNEG(I,O1,EP,N)*T_DELT_EIGEN(EP,N)
          ENDDO

C  REWORKED complex solutions

          DO EPC = 1, K_COMPLEX(N)
            K0 = 2*EPC - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            EP  = K_REAL(N) + EPC
            EPS = K_REAL(N) + K_COMPLEX(N) + EPC
            EM  = EP  + NSTKS_NSTRMS
            EMS = EPS + NSTKS_NSTRMS
            X1R = SOLA_XPOS(I,O1,K1,N)
            X1I = SOLA_XPOS(I,O1,K2,N)
            X2R = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K1,N)
     &           -SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K2,N)
            X2I = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K2,N)
     &           +SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K1,N)
            SMAT2(IROW,EP)  =   X1R
            SMAT2(IROW,EPS) = - X1I
            SMAT2(IROW,EM)  =   X2R
            SMAT2(IROW,EMS) = - X2I
          ENDDO

        ENDDO
      ENDDO

C  bottom BC (with albedo additions)
C  ---------------------------------

      IF ( DO_INCLUDE_SURFACE ) THEN
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CP = NSTKS_NSTRMS + IROW

C  real solutions

            DO EP = 1, K_REAL(N)
              CEP = EP
              CEM = CEP + NSTKS_NSTRMS
              XPNET = SOLA_XPOS(I1,O1,EP,N) - R2_HOMP(I,O1,EP)
              XMNET = SOLB_XNEG(I1,O1,EP,N) - R2_HOMM(I,O1,EP)
              SMAT2(CP,CEP) = T_DELT_EIGEN(EP,N) * XPNET
              SMAT2(CP,CEM) = XMNET
            ENDDO

C  Complex solutions

            DO EPC = 1, K_COMPLEX(N)
              K0 = 2*EPC - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              EP  = K_REAL(N) + EPC
              EPS = K_REAL(N) + K_COMPLEX(N) + EPC
              EM  = EP  + NSTKS_NSTRMS
              CEP = EP
              CEM = EM
              CEPS = EPS
              CEMS = CEPS + NSTKS_NSTRMS
              X1R_H = SOLA_XPOS(I1,O1,K1,N) - R2_HOMP(I,O1,K1)
              X1I_H = SOLA_XPOS(I1,O1,K2,N) - R2_HOMP(I,O1,K2)
              X2R = SOLB_XNEG(I1,O1,K1,N) - R2_HOMM(I,O1,K1)
              X2I = SOLB_XNEG(I1,O1,K2,N) - R2_HOMM(I,O1,K2)
              X1R = X1R_H * T_DELT_EIGEN(K1,N) -
     &              X1I_H * T_DELT_EIGEN(K2,N) 
              X1I = X1R_H * T_DELT_EIGEN(K2,N) +
     &              X1I_H * T_DELT_EIGEN(K1,N) 
              SMAT2(CP,CEP)  =   X1R
              SMAT2(CP,CEPS) = - X1I
              SMAT2(CP,CEM)  =   X2R
              SMAT2(CP,CEMS) = - X2I
            ENDDO

          ENDDO
        ENDDO

C  bottom BC (with No albedo)
C  --------------------------

      ELSE

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CP = NSTKS_NSTRMS + IROW

C  real solutions

            DO EP = 1, K_REAL(N)
              CEP = EP
              CEM = CEP + NSTKS_NSTRMS
              XPNET = SOLA_XPOS(I1,O1,EP,N)
              XMNET = SOLB_XNEG(I1,O1,EP,N)
              SMAT2(CP,CEP) = T_DELT_EIGEN(EP,N) * XPNET
              SMAT2(CP,CEM) = XMNET
            ENDDO

C  REWORKED Complex solutions

            DO EPC = 1, K_COMPLEX(N)
              K0 = 2*EPC - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              EP  = K_REAL(N) + EPC
              EPS = K_REAL(N) + K_COMPLEX(N) + EPC
              EM  = EP  + NSTKS_NSTRMS
              CEP = EP
              CEM = EM
              CEPS = EPS
              CEMS = CEPS + NSTKS_NSTRMS
              X1R = SOLA_XPOS(I1,O1,K1,N)*T_DELT_EIGEN(K1,N)
     &             -SOLA_XPOS(I1,O1,K2,N)*T_DELT_EIGEN(K2,N)
              X1I = SOLA_XPOS(I1,O1,K1,N)*T_DELT_EIGEN(K2,N)
     &             +SOLA_XPOS(I1,O1,K2,N)*T_DELT_EIGEN(K1,N)
              X2R = SOLB_XNEG(I1,O1,K1,N)
              X2I = SOLB_XNEG(I1,O1,K2,N)
              SMAT2(CP,CEP)  =   X1R
              SMAT2(CP,CEPS) = - X1I
              SMAT2(CP,CEM)  =   X2R
              SMAT2(CP,CEMS) = - X2I
            ENDDO

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

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of main model variable (I/O to this module)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'

C  subroutine arguments
C  --------------------

C  status flag

      INTEGER          STATUS

C  local variables
C  ---------------

      CHARACTER*(70)   MAIL, TRACE
      INTEGER          INFO
      CHARACTER*3      CI

C  debug
C      DOUBLE PRECISION SUM, SUM_R, SUM_C
C      COMPLEX*16       X1_C, X2_C

C  Intialize status

      STATUS = VLIDORT_SUCCESS

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
          TRACE = 'DGBTRF call in VLIDORT_BVPROBLEM'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ELSE IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MAIL  = 'argument i illegal value, for i = '//CI
          TRACE = 'DGBTRF call in VLIDORT_BVPROBLEM'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ENDIF

C  SVD the BVP matrix: No compression, Single Layer only
C  -----------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

C  LAPACK LU-decomposition for  matrix

        CALL DGETRF
     &     ( NTOTAL, NTOTAL, SMAT2, MAXSTRMSTKS_2, SIPIVOT, INFO )

C  (Error tracing)

        IF ( INFO .GT. 0 ) THEN
          MAIL='BVP LU decomposition (nlayers=1) (DGETRF)'
          CALL VLIDORT_LAPACK_ERROR ( INFO, 1, MAIL, STATUS )
          IF ( STATUS .EQ. VLIDORT_SERIOUS ) RETURN
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
     I         FOURIER_COMPONENT, IBEAM,
     I         SURFACE_FACTOR,
     I         BVP_SURFACE_SETUP_BEAM,
     O         STATUS )

C  include files
C  -------------

C  include file of dimensiopns and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of solution variables (local to this module)

      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'

C  arguments
C  ---------

C  input control

      LOGICAL          DO_INCLUDE_DIRECTBEAM
      LOGICAL          DO_INCLUDE_SURFACE
      LOGICAL          DO_INCLUDE_SURFEMISS
      INTEGER          FOURIER_COMPONENT, IBEAM
      DOUBLE PRECISION SURFACE_FACTOR

C  Surface reflectance function

      EXTERNAL         BVP_SURFACE_SETUP_BEAM

C  output status

      INTEGER          STATUS

C  Local variables
C  ---------------

      CHARACTER*(70)   MAIL, TRACE
      INTEGER          STATUS_SUB

C  --Additional setups for the albedo layer

      CALL BVP_SURFACE_SETUP_BEAM
     I       ( DO_INCLUDE_SURFACE,
     I         FOURIER_COMPONENT,
     I         SURFACE_FACTOR )

C  --set up Column for solution vector (the "B" as in AX=B)

      CALL BVP_COLUMN_SETUP
     I     ( DO_INCLUDE_SURFACE,
     I       DO_INCLUDE_DIRECTBEAM,
     I       DO_INCLUDE_SURFEMISS, IBEAM )

C  --Solve the boundary problem for this Fourier component (back substitution)

      CALL BVP_BACKSUB
     &        ( FOURIER_COMPONENT, IBEAM, STATUS_SUB )

C  error tracing

      IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
        MAIL  = 'Error return from module BVP_BACKSUB'
        TRACE = 'Call in BVP_SOLUTION_MASTER'
        STATUS = VLIDORT_SERIOUS
        CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
        RETURN
      ENDIF

C  finish

      RETURN
      END

C

      SUBROUTINE BVP_LAMBERTIAN_SURFACE_BEAM
     I     ( DO_INCLUDE_SURFACE,
     I       FOURIER_COMPONENT,
     I       SURFACE_FACTOR )

C  This is the Lambertian surface routine
C     Reflected beam solution

C  Include files
C  -------------

C  Additional sums for the final albedo-reflecting layer

C  include file of dimensiopns and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of solution stuff (i/o to this module)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'

C  control input
C  -------------

      LOGICAL          DO_INCLUDE_SURFACE
      INTEGER          FOURIER_COMPONENT
      DOUBLE PRECISION SURFACE_FACTOR

C  local variables
C  ---------------

      DOUBLE PRECISION REFL_B, AXJ
      INTEGER          I, O1, J

C  Initialization
C  ==============

C  Zero total reflected contributions

      DO I = 1, NSTREAMS
        DO O1 = 1, NSTOKES
          R2_BEAM(I,O1)    = ZERO
        ENDDO
      ENDDO

C  Return with Zeroed values if albedo flag not set

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

C  Return if not the Fourier m = 0 component

      IF ( FOURIER_COMPONENT .GT. 0 ) RETURN

C  Integrate Downward streams of particular solutions

      REFL_B = ZERO
      DO J = 1, NSTREAMS
        AXJ = QUAD_STRMWTS(J)
        REFL_B = REFL_B + AXJ*WLOWER(J,1,NLAYERS)
      ENDDO

C  set reflectance

      REFL_B = REFL_B * SURFACE_FACTOR * LAMBERTIAN_ALBEDO
      DO I = 1, NSTREAMS
        R2_BEAM(I,1) = REFL_B
      ENDDO

C  debug

C      DO I = 1, NSTREAMS
C        DO O1 = 1, NSTOKES
C          write(*,'(2i3,1p2e14.6)')I,O1,R2_BEAM(I,O1)
C        ENDDO
C      ENDDO
C      pause'r2beam'
      
C  Finish

      RETURN
      END

C

      SUBROUTINE BVP_COLUMN_SETUP
     I     ( DO_INCLUDE_SURFACE,
     I       DO_INCLUDE_DIRECTBEAM,
     I       DO_INCLUDE_SURFEMISS, IBEAM )

C  include files
C  -------------

C  include file of dimensiopns and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of solution variables (local to this module)

      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  input arguments
C  ---------------

      LOGICAL         DO_INCLUDE_DIRECTBEAM
      LOGICAL         DO_INCLUDE_SURFACE
      LOGICAL         DO_INCLUDE_SURFEMISS
      INTEGER         IBEAM

C  Local variables
C  ---------------

      INTEGER          I,I1,LAY,LAY1,C0,CM,O1,IR,IROW,NSTOKES_FIRST
      double precision local_emiss, FP_SBB

C  If Nlayers = 1, go to special case

      IF ( NLAYERS .EQ. 1 ) GO TO 345

C  zero column vector

      DO I = 1, NTOTAL
        COL2(I,IBEAM) = ZERO
      ENDDO

C  Upper boundary for layer 1: no downward diffuse radiation
C  ---------------------------------------------------------

      LAY = 1
      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        DO O1 = 1, NSTOKES
          IROW = IR + O1
          COL2(IROW,IBEAM)   = - WUPPER(I,O1,LAY)
        ENDDO
      ENDDO

C  intermediate layer boundaries (will not be done if NLAYERS = 1 )
C  -----------------------------

      DO LAY = 2, NLAYERS

        LAY1 = LAY - 1
        C0 = LAY1*NSTKS_NSTRMS_2 - NSTKS_NSTRMS
        DO I = 1, NSTREAMS_2
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            COL2(CM,IBEAM) = WUPPER(I,O1,LAY) - WLOWER(I,O1,LAY1)
          ENDDO
        ENDDO

      ENDDO

C  lowest (surface) boundary with albedo (diffuse radiation terms only)
C  -------------------------------------

      LAY = NLAYERS
      C0 = (LAY-1)*NSTKS_NSTRMS_2 + NSTKS_NSTRMS

C  with non-zero albedo, include integrated downward reflectances

      IF ( DO_INCLUDE_SURFACE ) THEN

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            COL2(CM,IBEAM) = - WLOWER(I1,O1,LAY) + R2_BEAM(I,O1)
          ENDDO
        ENDDO

C  no albedo, similar code excluding integrated reflectance

      ELSE

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            COL2(CM,IBEAM) = - WLOWER(I1,O1,LAY)
          ENDDO
        ENDDO

      ENDIF

C  Add direct beam solution (only to final level)
C  ----------------------------------------------

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        IF ( DO_INCLUDE_SURFACE ) THEN
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              COL2(CM,IBEAM) = COL2(CM,IBEAM) + DIRECT_BEAM(I,IBEAM,O1)
            ENDDO
          ENDDO
        ENDIF
      ENDIF

C  Add thermal emission of ground surface (only to final level)
C  ------------------------------------------------------------

      IF ( DO_INCLUDE_SURFEMISS ) THEN
        O1 = 1
        local_emiss = one - lambertian_albedo
        FP_SBB = SURFBB
        IF ( DO_SOLAR_SOURCES ) FP_SBB = PI4 * SURFBB
        IF ( do_lambertian_surface ) THEN
          DO I = 1, NSTREAMS
            CM = C0 + NSTOKES*(I-1) + O1
            COL2(CM,IBEAM) = COL2(CM,IBEAM) + FP_SBB * local_emiss
          ENDDO
        ELSE
          DO I = 1, NSTREAMS
            CM = C0 + NSTOKES*(I-1) + O1
            COL2(CM,IBEAM) = COL2(CM,IBEAM) + FP_SBB * EMISSIVITY(I,O1)
          ENDDO
        ENDIF
      ENDIF

C  debug

c      do i = 1, ntotal
c        write(77,*)i,COL2(i,IBEAM)
c      enddo
c      pause'hello'

C  normal completion

      RETURN

C  special case

345   CONTINUE

C  zero column vector

      DO I = 1, NTOTAL
        SCOL2(I,IBEAM) = ZERO
      ENDDO

C  Upper boundary for layer 1: no downward diffuse radiation
C  ---------------------------------------------------------

      LAY = 1
      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        DO O1 = 1, NSTOKES
          IROW = IR + O1
          SCOL2(IROW,IBEAM)   = - WUPPER(I,O1,LAY)
        ENDDO
      ENDDO

C  lowest (surface) boundary with albedo (diffuse radiation terms only)
C  -------------------------------------------------------------------

C  with non-zero albedo, include integrated downward reflectances
C  no albedo, similar code excluding integrated reflectance

      C0 = NSTKS_NSTRMS
      IF ( DO_INCLUDE_SURFACE ) THEN
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            SCOL2(CM,IBEAM) = - WLOWER(I1,O1,LAY) + R2_BEAM(I,O1)
          ENDDO
        ENDDO
      ELSE
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            SCOL2(CM,IBEAM) = - WLOWER(I1,O1,LAY)
          ENDDO
        ENDDO
      ENDIF

C  Add direct beam solution (only to final level)
C  ----------------------------------------------

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        IF ( DO_INCLUDE_SURFACE ) THEN
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              SCOL2(CM,IBEAM) = SCOL2(CM,IBEAM)+DIRECT_BEAM(I,IBEAM,O1)
            ENDDO
          ENDDO
        ENDIF
      ENDIF

C  Add thermal emission of ground surface (only to final level)
C  ------------------------------------------------------------

      IF ( DO_INCLUDE_SURFEMISS ) THEN
        NSTOKES_FIRST = 1
        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES_FIRST
            IROW = IR + O1
            CM = C0 + IROW
            SCOL2(CM,IBEAM) = SCOL2(CM,IBEAM) + SURFBB*EMISSIVITY(I,O1)
          ENDDO
        ENDDO
      ENDIF

C
C  finish

      RETURN
      END

C

      SUBROUTINE BVP_BACKSUB
     I     ( FOURIER, IBEAM,
     O       STATUS )

C  Solves the boundary value problem.

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of main model variable (I/O to this module)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'

C  subroutine arguments
C  --------------------

C  FOURIER and Beam index

      INTEGER          FOURIER, IBEAM

C  status flag

      INTEGER          STATUS

C  local variables
C  ---------------

      INTEGER          C0, LAY, K, KO1, K0, K1, K2
      INTEGER          IROW, IROW1, IROW_S, IROW1_S
      CHARACTER*(70)   MAIL, TRACE
      INTEGER          INFO
      CHARACTER*3      CI

C  Intialize status

      STATUS = VLIDORT_SUCCESS

C  BVP back-substitution: With compression (multilayers)
C  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

C  LAPACK substitution (DGBTRS) using RHS column vector COL2

        CALL DGBTRS
     &     ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, IBEAM,
     &        BANDMAT2, MAXBANDTOTAL, IPIVOT,
     &        COL2, MAXTOTAL, INFO )

C  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MAIL  = 'argument i illegal value, for i = '//CI
          TRACE = 'DGBTRS call in BVP_BACKSUB'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
         ENDIF

C  Set integration constants LCON and MCON for -/+ eigensolutions, all layers

        DO LAY = 1, NLAYERS
          C0 = (LAY-1)*NSTKS_NSTRMS_2
          KO1 = K_REAL(LAY) + 1
          DO K = 1, K_REAL(LAY)
            IROW = K
            IROW1 = IROW + NSTKS_NSTRMS
            LCON(K,LAY) = COL2(C0+IROW,IBEAM)
            MCON(K,LAY) = COL2(C0+IROW1,IBEAM)
          ENDDO
          DO K = 1, K_COMPLEX(LAY)
            K0 = 2*K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            IROW    = K + K_REAL(LAY)
            IROW1   = IROW + NSTKS_NSTRMS
            IROW_S  = IROW + K_COMPLEX(LAY)
            IROW1_S = IROW_S + NSTKS_NSTRMS
            LCON(K1,LAY) = COL2(C0+IROW,   IBEAM)
            LCON(K2,LAY) = COL2(C0+IROW_S, IBEAM)
            MCON(K1,LAY) = COL2(C0+IROW1,  IBEAM)
            MCON(K2,LAY) = COL2(C0+IROW1_S,IBEAM)
          ENDDO
        ENDDO


C  debug-----------------------------------------------

c        if ( do_debug_write) then
c         if (fourier.eq.1 ) then
c         DO LAY = 1, NLAYERS
c          KO1 = K_REAL(LAY) + 1
c          DO K = 1, K_REAL(LAY)
c           write(50,'(4i3,1p2e20.10)')FOURIER,IBEAM,K,LAY,
c     &                LCON(K,LAY), MCON(K,LAY)
c          ENDDO
c          DO K = 1, K_COMPLEX(LAY)
c           K0 = 2*K - 2
c           K1 = KO1 + K0
c           K2 = K1  + 1
c           write(50,'(4i3,1p4e20.10)')FOURIER,IBEAM,K,LAY,
c     &                LCON(K1,LAY), MCON(K1,LAY),
c     &                LCON(K2,LAY), MCON(K2,LAY)
c          ENDDO
c         ENDDO
c        endif
c        ENDIF

C        DO LAY = 20, NLAYERS
C          DO K = 1, K_REAL(LAY)
C             write(97,'(3i5,1p2e15.7)')FOURIER,LAY,K,
C     &                LCON(K,LAY),MCON(K,LAY)
C            ENDDO
C           ENDDO
C          IF ( FOURIER.EQ.0)STOP

c        write(*,*)'reg',fourier,ibeam,lcon(1,7)

C  Solve the boundary problem: No compression, Single Layer only
C  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

C  LAPACK substitution (DGETRS) using RHS column vector SCOL2

        CALL DGETRS
     &     ( 'N', NTOTAL, IBEAM, SMAT2, MAXSTRMSTKS_2, SIPIVOT,
     &        SCOL2, MAXSTRMSTKS_2, INFO )

C  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          MAIL='BVP back-substitution (nlayers=1) (DGETRS)'
          CALL VLIDORT_LAPACK_ERROR ( INFO, 1, MAIL, STATUS )
          IF ( STATUS .EQ. VLIDORT_SERIOUS ) RETURN
        ENDIF

C  Set integration constants LCON and MCON for -/+ eigensolutions, all layers

        LAY = 1
        KO1 = K_REAL(LAY) + 1
        DO K = 1, K_REAL(LAY)
          IROW = K
          IROW1 = IROW + NSTKS_NSTRMS
          LCON(K,LAY) = SCOL2(IROW,IBEAM)
          MCON(K,LAY) = SCOL2(IROW1,IBEAM)
        ENDDO
        DO K = 1, K_COMPLEX(LAY)
          K0 = 2*K - 2
          K1 = KO1 + K0
          K2 = K1  + 1
          IROW = K + K_REAL(LAY)
          IROW1 = IROW + NSTKS_NSTRMS
          IROW_S = K + K_REAL(LAY) + K_COMPLEX(LAY)
          IROW1_S = IROW_S + NSTKS_NSTRMS
          LCON(K1,LAY) = SCOL2(IROW,    IBEAM)
          LCON(K2,LAY) = SCOL2(IROW_S,  IBEAM)
          MCON(K1,LAY) = SCOL2(IROW1,   IBEAM)
          MCON(K2,LAY) = SCOL2(IROW1_S, IBEAM)
        ENDDO

C  debug-----------------------------------------------

        if ( do_debug_write.and.do_fdtest ) then
         DO LAY = 1, NLAYERS
          KO1 = K_REAL(LAY) + 1
          DO K = 1, K_REAL(LAY)
           write(40,'(4i3,1p2e20.10)')FOURIER,IBEAM,K,LAY,
     &                LCON(K,LAY), MCON(K,LAY)
          ENDDO
          DO K = 1, K_COMPLEX(LAY)
           K0 = 2*K - 2
           K1 = KO1 + K0
           K2 = K1  + 1
           write(50,'(4i3,1p4e20.10)')FOURIER,IBEAM,K,LAY,
     &                LCON(K1,LAY), MCON(K1,LAY),
     &                LCON(K2,LAY), MCON(K2,LAY)
          ENDDO
         ENDDO
        ENDIF

c      IF ( DO_DEBUG_WRITE ) THEN
c       IF ( DO_FDTEST ) THEN
c        DO N = 1, NLAYERS
c          DO I = 1, NSTKS_NSTRMS
c           write(30,'(3i3,1p2e20.10)')
c     &         FOURIER,IBEAM,I,LCON(I,N),MCON(I,N)
c          ENDDO
c        ENDDO
c       ENDIF
c      ENDIF

      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE BVPTEL_MATRIXSETUP_MASTER
     I     ( DO_INCLUDE_SURFACE,
     I       FOURIER_COMPONENT,
     I       SURFACE_FACTOR,
     I       BVP_SURFACE_SETUP_HOM,
     O       STATUS )

C  Sets up the telescoped boundary value problem.
C    Standard case: Fourier > 0. With surface reflection term

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and solution variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'

C  control input
C  -------------

      LOGICAL          DO_INCLUDE_SURFACE
      INTEGER          FOURIER_COMPONENT
      DOUBLE PRECISION SURFACE_FACTOR

C  Surface reflectance function

      EXTERNAL         BVP_SURFACE_SETUP_HOM

C  status output

      INTEGER          STATUS

C  local variables
C  ---------------

      CHARACTER*70     MAIL, TRACE
      INTEGER          STATUS_SUB

C  start code
C  ----------

C  Intialize status

      STATUS = VLIDORT_SUCCESS

C  initialize compression matrix (Do this for every Fourier component)

      CALL BVPTEL_MATRIX_INIT ( FOURIER_COMPONENT )

C  Do the surface setup if required

C  Necessary to have reflected solutions in the lower boundary
C  layer, even when BVP Telescoping only for Lambertian albedo case.

C  Specialist option 2 only, and only if the lowest active layer
C    is the lowest atmospheric layer. Lambertian only.

      IF ( DO_SPECIALIST_OPTION_2 ) THEN
        IF ( ACTIVE_LAYERS(NLAYERS_TEL).EQ.NLAYERS ) THEN
          CALL BVP_SURFACE_SETUP_HOM
     I    ( DO_INCLUDE_SURFACE,
     I      FOURIER_COMPONENT,
     I      SURFACE_FACTOR )
        ENDIF
      ENDIF

C  set up boundary values matrix in compressed form (the "A" as in AX=B)

      CALL BVPTEL_MATRIX_SETUP
     I     ( DO_INCLUDE_SURFACE )

C  SVD decomposition of compressed boundary values matrix

      CALL BVPTEL_MATRIX_SVD ( STATUS_SUB )

C  error tracing

      IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
        MAIL = 'Error return from module BVPTEL_MATRIX_SVD'
        TRACE= 'Call in BVPTEL_MATRIXSETUP_MASTER'
        STATUS = VLIDORT_SERIOUS
        CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
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

      INCLUDE '../includes/VLIDORT.PARS'

C  Include file of input variables
C  include file of bookkeeping inputs

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include file of setup variables

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  include file of solution stuff (i/o to this module)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'

C  Subroutine argument
C  -------------------

C  Fourier component

      INTEGER             FOURIER_COMPONENT

C  Local variables
C  ---------------

      INTEGER             I, J, NS, N, N3

C  compression function

      INTEGER  BTELMAT_ROWMASK
      EXTERNAL BTELMAT_ROWMASK

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

        N_BVTELMATRIX_SIZE = 2 * NSTKS_NSTRMS * NLAYERS_TEL

C  Exit if only one active layer

        IF ( NLAYERS_TEL .EQ. 1 ) GO TO 345

C  Number of sub and super diagonals
  
        N_BVTELMATRIX_SUPDIAG = 3 * NSTKS_NSTRMS - 1
        N_BVTELMATRIX_SUBDIAG = 3 * NSTKS_NSTRMS - 1

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

        KALLTEL = N_BVTELMATRIX_SUBDIAG + N_BVTELMATRIX_SUPDIAG + 1

C  Avoid fancy zeroing - adopt kludge
C   Potential Danger point

        DO I = 1, MAXBANDTOTAL
          DO J = 1, N_BVTELMATRIX_SIZE
            BANDTELMAT2(I,J) = ZERO
          ENDDO
        ENDDO

C  Control

        GO TO 346

C  control point

 345    CONTINUE

C  single layer setting

        N_BVTELMATRIX_SIZE = 2 * NSTKS_NSTRMS

 346    continue

C  reset

        DO_BVTEL_INITIAL = .FALSE.

C  end initialization

      ENDIF

C  finish

      RETURN
      END

C

      SUBROUTINE BVPTEL_MATRIX_SETUP
     I     (  DO_INCLUDE_SURFACE )

C  Fills up the matrix directly (compressed or 1-layer)

C  include files
C  -------------

C  include file of dimensiopns and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include file of input variables
C  include file of bookkeeping inputs

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  Include files of setup variables

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  include file of reflectance variables

      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'

C  include file of solution BVP variables (output to this module)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'

C  input arguments
C  ---------------

      LOGICAL             DO_INCLUDE_SURFACE

C  Local variables
C  ---------------

      INTEGER          I, J, I1, N, N1, O1, IR, IROW, NS
      INTEGER          EP, EM, EPC, EPS, EMS
      INTEGER          CP, CM, CEP, CEM, CEPS, CEMS, CEP1, CEM1
      DOUBLE PRECISION XPNET, XMNET, X1R, X2R, X1I, X2I, X1R_H, X1I_H
      INTEGER          KO1, K0, K1, K2, KO11
      INTEGER          CE_OFFSET, C0

C  compression function

      INTEGER  BTELMAT_ROWMASK
      EXTERNAL BTELMAT_ROWMASK

C  Initialization

      EP = 0

C  If Nlayers = 1, go to special case

      IF ( NLAYERS_TEL .EQ. 1 ) GO TO 345

C  top BC for first active layer 1: no downward diffuse radiation

      NS = 1
      N = ACTIVE_LAYERS(NS)
      KO1 = K_REAL(N) + 1
      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        DO O1 = 1, NSTOKES
          IROW = IR + O1

C  real solutions

          DO EP = 1, K_REAL(N)
            EM = EP + NSTKS_NSTRMS
            BANDTELMAT2(BTELMAT_ROWMASK(IROW,EP),EP)  =
     &                   SOLA_XPOS(I,O1,EP,N)
            BANDTELMAT2(BTELMAT_ROWMASK(IROW,EM),EM)  =
     &                   SOLB_XNEG(I,O1,EP,N)*T_DELT_EIGEN(EP,N)
          ENDDO

C  complex solutions (REWORKED)

          DO EPC = 1, K_COMPLEX(N)
            K0 = 2*EPC - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            EP  = K_REAL(N) + EPC
            EPS = K_REAL(N) + K_COMPLEX(N) + EPC
            EM  = EP  + NSTKS_NSTRMS
            EMS = EPS + NSTKS_NSTRMS
            X1R = SOLA_XPOS(I,O1,K1,N)
            X1I = SOLA_XPOS(I,O1,K2,N)
            X2R = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K1,N)
     &           -SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K2,N)
            X2I = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K2,N)
     &           +SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K1,N)
            BANDTELMAT2(BTELMAT_ROWMASK(IROW,EP),EP)   =   X1R
            BANDTELMAT2(BTELMAT_ROWMASK(IROW,EPS),EPS) = - X1I
            BANDTELMAT2(BTELMAT_ROWMASK(IROW,EM),EM)   =   X2R
            BANDTELMAT2(BTELMAT_ROWMASK(IROW,EMS),EMS) = - X2I
          ENDDO

        ENDDO
      ENDDO

C  intermediate layer boundaries

      C0 = - NSTKS_NSTRMS
      DO NS = 2, NLAYERS_TEL
        N = ACTIVE_LAYERS(NS)
        N1 = N - 1
        KO1  = K_REAL(N) + 1
        KO11 = K_REAL(N1) + 1
        C0   = C0 + NSTKS_NSTRMS_2
        CE_OFFSET = C0 - NSTKS_NSTRMS
        DO I = 1, NSTREAMS_2
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW

C  real solutions for layer above the boundary

            DO EP = 1, K_REAL(N1)
              CEP = CE_OFFSET + EP
              CEM = CEP + NSTKS_NSTRMS
              BANDTELMAT2(BTELMAT_ROWMASK(CM,CEP),CEP)   =
     &                 T_DELT_EIGEN(EP,N1)*SOLA_XPOS(I,O1,EP,N1)
              BANDTELMAT2(BTELMAT_ROWMASK(CM,CEM),CEM)   =
     &                   SOLB_XNEG(I,O1,EP,N1)
            ENDDO

C  real solutions for layer below the boundary
C   ( Note the change of sign !!! )

            DO EP = 1, K_REAL(N)
              CEP = CE_OFFSET + EP
              CEM = CEP + NSTKS_NSTRMS
              CEP1 = CEP + NSTKS_NSTRMS_2
              CEM1 = CEM + NSTKS_NSTRMS_2
              BANDTELMAT2(BTELMAT_ROWMASK(CM,CEP1),CEP1) =
     &                   - SOLA_XPOS(I,O1,EP,N)
              BANDTELMAT2(BTELMAT_ROWMASK(CM,CEM1),CEM1) =
     &                 - T_DELT_EIGEN(EP,N)*SOLB_XNEG(I,O1,EP,N)
            ENDDO

C  complex solutions for layer above boundary

            DO EPC = 1, K_COMPLEX(N1)
              K0 = 2*EPC - 2
              K1 = KO11 + K0
              K2 = K1  + 1
              EP  = K_REAL(N1) + EPC
              EPS = K_REAL(N1) + K_COMPLEX(N1) + EPC
              CEP  = CE_OFFSET + EP
              CEM  = CEP + NSTKS_NSTRMS
              CEPS = CE_OFFSET + EPS
              CEMS = CEPS + NSTKS_NSTRMS
              X1R = SOLA_XPOS(I,O1,K1,N1)*T_DELT_EIGEN(K1,N1)
     &             -SOLA_XPOS(I,O1,K2,N1)*T_DELT_EIGEN(K2,N1)
              X1I = SOLA_XPOS(I,O1,K1,N1)*T_DELT_EIGEN(K2,N1)
     &             +SOLA_XPOS(I,O1,K2,N1)*T_DELT_EIGEN(K1,N1)
              X2R = SOLB_XNEG(I,O1,K1,N1)
              X2I = SOLB_XNEG(I,O1,K2,N1)
              BANDTELMAT2(BTELMAT_ROWMASK(CM,CEP),CEP)   =   X1R
              BANDTELMAT2(BTELMAT_ROWMASK(CM,CEPS),CEPS) = - X1I
              BANDTELMAT2(BTELMAT_ROWMASK(CM,CEM),CEM)   =   X2R
              BANDTELMAT2(BTELMAT_ROWMASK(CM,CEMS),CEMS) = - X2I
            ENDDO

C  complex solutions for layer below boundary
C   ( Note the change of sign !!! )

            DO EPC = 1, K_COMPLEX(N)
              K0 = 2*EPC - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              EP  = K_REAL(N) + EPC
              EPS = K_REAL(N) + K_COMPLEX(N) + EPC
              CEP  = CE_OFFSET + EP + NSTKS_NSTRMS_2
              CEM  = CEP + NSTKS_NSTRMS
              CEPS = CE_OFFSET + EPS + NSTKS_NSTRMS_2
              CEMS = CEPS + NSTKS_NSTRMS
              X1R = SOLA_XPOS(I,O1,K1,N)
              X1I = SOLA_XPOS(I,O1,K2,N)
              X2R = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K1,N)
     &             -SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K2,N)
              X2I = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K2,N)
     &             +SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K1,N)
              BANDTELMAT2(BTELMAT_ROWMASK(CM,CEP),CEP)   = - X1R
              BANDTELMAT2(BTELMAT_ROWMASK(CM,CEPS),CEPS) =   X1I
              BANDTELMAT2(BTELMAT_ROWMASK(CM,CEM),CEM)   = - X2R
              BANDTELMAT2(BTELMAT_ROWMASK(CM,CEMS),CEMS) =   X2I
            ENDDO

          ENDDO
        ENDDO
      ENDDO

C  bottom BC for Lowest active layer
C   Normally no surface additions, except Specialist # 2

      N = ACTIVE_LAYERS(NLAYERS_TEL)
      KO1 = K_REAL(N) + 1
      C0  = C0 + NSTKS_NSTRMS_2
      CE_OFFSET = C0 - NSTKS_NSTRMS

      IF ( DO_INCLUDE_SURFACE.AND.DO_SPECIALIST_OPTION_2 ) THEN

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CP = C0 + IROW

C  real solutions

            DO EP = 1, K_REAL(N)
              CEP = CE_OFFSET + EP
              CEM = CEP + NSTKS_NSTRMS
              XPNET = SOLA_XPOS(I1,O1,EP,N) - R2_HOMP(I,O1,EP)
              XMNET = SOLB_XNEG(I1,O1,EP,N) - R2_HOMM(I,O1,EP)
              BANDTELMAT2(BTELMAT_ROWMASK(CP,CEP),CEP) =
     &                   T_DELT_EIGEN(EP,N) * XPNET
              BANDTELMAT2(BTELMAT_ROWMASK(CP,CEM),CEM) =
     &                   XMNET
            ENDDO

C  Complex solutions

            DO EPC = 1, K_COMPLEX(N)
              K0 = 2*EPC - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              EP  = K_REAL(N) + EPC
              EPS = K_REAL(N) + K_COMPLEX(N) + EPC
              CEP = CE_OFFSET + EP
              CEM = CEP + NSTKS_NSTRMS
              CEPS = CE_OFFSET + EPS
              CEMS = CEPS + NSTKS_NSTRMS
              X1R_H = SOLA_XPOS(I1,O1,K1,N) - R2_HOMP(I,O1,K1)
              X1I_H = SOLA_XPOS(I1,O1,K2,N) - R2_HOMP(I,O1,K2)
              X2R = SOLB_XNEG(I1,O1,K1,N) - R2_HOMM(I,O1,K1)
              X2I = SOLB_XNEG(I1,O1,K2,N) - R2_HOMM(I,O1,K2)
              X1R = X1R_H * T_DELT_EIGEN(K1,N) -
     &              X1I_H * T_DELT_EIGEN(K2,N)
              X1I = X1R_H * T_DELT_EIGEN(K2,N) +
     &              X1I_H * T_DELT_EIGEN(K1,N)
              BANDTELMAT2(BTELMAT_ROWMASK(CP,CEP),CEP)   =   X1R
              BANDTELMAT2(BTELMAT_ROWMASK(CP,CEPS),CEPS) = - X1I
              BANDTELMAT2(BTELMAT_ROWMASK(CP,CEM),CEM)   =   X2R
              BANDTELMAT2(BTELMAT_ROWMASK(CP,CEMS),CEMS) = - X2I
            ENDDO

C  End stokes and streams loops

          ENDDO
        ENDDO

C  No surface albedo reflections

      ELSE

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CP = C0 + IROW

C  real solutions

            DO EP = 1, K_REAL(N)
              CEP = CE_OFFSET + EP
              CEM = CEP + NSTKS_NSTRMS
              XPNET = SOLA_XPOS(I1,O1,EP,N)
              XMNET = SOLB_XNEG(I1,O1,EP,N)
              BANDTELMAT2(BTELMAT_ROWMASK(CP,CEP),CEP) =
     &                   T_DELT_EIGEN(EP,N) * XPNET
              BANDTELMAT2(BTELMAT_ROWMASK(CP,CEM),CEM) =
     &                   XMNET
            ENDDO

C  Complex solutions

            DO EPC = 1, K_COMPLEX(N)
              K0 = 2*EPC - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              EP  = K_REAL(N) + EPC
              EPS = K_REAL(N) + K_COMPLEX(N) + EPC
              CEP  = CE_OFFSET + EP
              CEM  = CEP + NSTKS_NSTRMS
              CEPS = CE_OFFSET + EPS
              CEMS = CEPS + NSTKS_NSTRMS
              X1R = SOLA_XPOS(I1,O1,K1,N) * T_DELT_EIGEN(K1,N)
     &             -SOLA_XPOS(I1,O1,K2,N) * T_DELT_EIGEN(K2,N)
              X1I = SOLA_XPOS(I1,O1,K1,N) * T_DELT_EIGEN(K2,N)
     &             +SOLA_XPOS(I1,O1,K2,N) * T_DELT_EIGEN(K1,N)
              X2R = SOLB_XNEG(I1,O1,K1,N)
              X2I = SOLB_XNEG(I1,O1,K2,N)
              BANDTELMAT2(BTELMAT_ROWMASK(CP,CEP),CEP)   =   X1R
              BANDTELMAT2(BTELMAT_ROWMASK(CP,CEPS),CEPS) = - X1I
              BANDTELMAT2(BTELMAT_ROWMASK(CP,CEM),CEM)   =   X2R
              BANDTELMAT2(BTELMAT_ROWMASK(CP,CEMS),CEMS) = - X2I
           ENDDO

          ENDDO
        ENDDO

      ENDIF

C  normal completion

      RETURN

C  special case. Only 1 active layer

345   CONTINUE

C  Set up BVP matrix for the active layer
C  ======================================

C  active layer for telescoped BVP

      N = ACTIVE_LAYERS(1)
      KO1 = K_REAL(N) + 1

C  initialize using the SMAT2 matrix

      DO I = 1, NSTKS_NSTRMS_2
        DO J = 1, NSTKS_NSTRMS_2
          SMAT2(I,J)     = ZERO
        ENDDO
      ENDDO

C  top of the layer (downwelling only)
C  -----------------------------------

      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        DO O1 = 1, NSTOKES
          IROW = IR + O1

C  real solutions

          DO EP = 1, K_REAL(N)
            EM = EP + NSTKS_NSTRMS
            SMAT2(IROW,EP) = SOLA_XPOS(I,O1,EP,N)
            SMAT2(IROW,EM) = SOLB_XNEG(I,O1,EP,N)*T_DELT_EIGEN(EP,N)
          ENDDO

C  REWORKED complex solutions

          DO EPC = 1, K_COMPLEX(N)
            K0 = 2*EPC - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            EP  = K_REAL(N) + EPC
            EPS = K_REAL(N) + K_COMPLEX(N) + EPC
            EM  = EP  + NSTKS_NSTRMS
            EMS = EPS + NSTKS_NSTRMS
            X1R = SOLA_XPOS(I,O1,K1,N)
            X1I = SOLA_XPOS(I,O1,K2,N)
            X2R = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K1,N)
     &           -SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K2,N)
            X2I = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K2,N)
     &           +SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K1,N)
            SMAT2(IROW,EP)  =   X1R
            SMAT2(IROW,EPS) = - X1I
            SMAT2(IROW,EM)  =   X2R
            SMAT2(IROW,EMS) = - X2I
          ENDDO

        ENDDO
      ENDDO

C  bottom BC. Only albedo additions for specialist option 2
C  --------------------------------------------------------

      IF ( DO_INCLUDE_SURFACE.AND.DO_SPECIALIST_OPTION_2 ) THEN
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CP = NSTKS_NSTRMS + IROW

C  real solutions

            DO EP = 1, K_REAL(N)
              CEP = EP
              CEM = CEP + NSTKS_NSTRMS
              XPNET = SOLA_XPOS(I1,O1,EP,N) - R2_HOMP(I,O1,EP)
              XMNET = SOLB_XNEG(I1,O1,EP,N) - R2_HOMM(I,O1,EP)
              SMAT2(CP,CEP) = T_DELT_EIGEN(EP,N) * XPNET
              SMAT2(CP,CEM) = XMNET
            ENDDO

C  Complex solutions

            DO EPC = 1, K_COMPLEX(N)
              K0 = 2*EPC - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              EP  = K_REAL(N) + EPC
              EPS = K_REAL(N) + K_COMPLEX(N) + EPC
              EM  = EP  + NSTKS_NSTRMS
              CEP = EP
              CEM = EM
              CEPS = EPS
              CEMS = CEPS + NSTKS_NSTRMS
              X1R_H = SOLA_XPOS(I1,O1,K1,N) - R2_HOMP(I,O1,K1)
              X1I_H = SOLA_XPOS(I1,O1,K2,N) - R2_HOMP(I,O1,K2)
              X2R = SOLB_XNEG(I1,O1,K1,N) - R2_HOMM(I,O1,K1)
              X2I = SOLB_XNEG(I1,O1,K2,N) - R2_HOMM(I,O1,K2)
              X1R = X1R_H * T_DELT_EIGEN(K1,N) -
     &              X1I_H * T_DELT_EIGEN(K2,N) 
              X1I = X1R_H * T_DELT_EIGEN(K2,N) +
     &              X1I_H * T_DELT_EIGEN(K1,N) 
              SMAT2(CP,CEP)  =   X1R
              SMAT2(CP,CEPS) = - X1I
              SMAT2(CP,CEM)  =   X2R
              SMAT2(CP,CEMS) = - X2I
            ENDDO

          ENDDO
        ENDDO

C  bottom BC (with No albedo). This is the usual situation
C  -------------------------------------------------------

      ELSE

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CP = NSTKS_NSTRMS + IROW

C  real solutions

            DO EP = 1, K_REAL(N)
              CEP = EP
              CEM = CEP + NSTKS_NSTRMS
              XPNET = SOLA_XPOS(I1,O1,EP,N)
              XMNET = SOLB_XNEG(I1,O1,EP,N)
              SMAT2(CP,CEP) = T_DELT_EIGEN(EP,N) * XPNET
              SMAT2(CP,CEM) = XMNET
            ENDDO

C  REWORKED Complex solutions
C    Note to self. The second set of solutions seems correct.

            DO EPC = 1, K_COMPLEX(N)
              K0 = 2*EPC - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              EP  = K_REAL(N) + EPC
              EPS = K_REAL(N) + K_COMPLEX(N) + EPC
              EM  = EP  + NSTKS_NSTRMS
              CEP = EP
              CEM = EM
              CEPS = EPS
              CEMS = CEPS + NSTKS_NSTRMS
              X1R = SOLA_XPOS(I1,O1,K1,N)*T_DELT_EIGEN(K1,N)
     &             -SOLA_XPOS(I1,O1,K2,N)*T_DELT_EIGEN(K2,N)
              X1I = SOLA_XPOS(I1,O1,K1,N)*T_DELT_EIGEN(K2,N)
     &             +SOLA_XPOS(I1,O1,K2,N)*T_DELT_EIGEN(K1,N)
              X2R = SOLB_XNEG(I1,O1,K1,N)
              X2I = SOLB_XNEG(I1,O1,K2,N)
              IF ( N.EQ.NLAYERS ) THEN
                SMAT2(CP,CEP)  =   X1R   ! apparently wrong
                SMAT2(CP,CEPS) = - X1I   ! apparently wrong
                SMAT2(CP,CEM)  =   X2R   ! apparently wrong
                SMAT2(CP,CEMS) = - X2I   ! apparently wrong
              ELSE
                SMAT2(CP,CEP)  =   X1R  ! tested 29 December 2005
                SMAT2(CP,CEPS) = - X1I  ! tested 29 December 2005
                SMAT2(CP,CEM)  =   X2R  ! tested 29 December 2005
                SMAT2(CP,CEMS) = - X2I  ! tested 29 December 2005
              ENDIF
            ENDDO

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

      INCLUDE '../includes/VLIDORT.PARS'

C  Include file of input variables
C  include file of bookkeeping inputs

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include file of solution variables (I/O to this module)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'

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

      STATUS = VLIDORT_SUCCESS

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
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ELSE IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MAIL  = 'argument i illegal value, for i = '//CI
          TRACE = 'DGBTRF call in BVPTEL_MATRIX_SVD'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ENDIF

C  SVD the BVP matrix: No compression, Single Layer only
C  -----------------------------------------------------

      ELSE IF ( NLAYERS_TEL .EQ. 1 ) THEN

C  LAPACK LU-decomposition for single layer matrix

        CALL DGETRF
     &     (  NSTKS_NSTRMS_2, NSTKS_NSTRMS_2,
     &        SMAT2, MAXSTRMSTKS_2, SIPIVOT, INFO )

C  (Error tracing)

        IF ( INFO .GT. 0 ) THEN
          MAIL='Failure: BVP LU decomposition (telescoped) (DGETRF)'
          CALL VLIDORT_LAPACK_ERROR ( INFO, 1, MAIL, STATUS )
          IF ( STATUS .EQ. VLIDORT_SERIOUS ) RETURN
        ENDIF

      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE BVPTEL_SOLUTION_MASTER
     I       ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM,
     I         FOURIER_COMPONENT, IBEAM,
     I         SURFACE_FACTOR,
     I         BVP_SURFACE_SETUP_BEAM,
     O         STATUS )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and solution variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'

C  input arguments
C  ---------------

C  inclusion of surface

      LOGICAL          DO_INCLUDE_SURFACE
      LOGICAL          DO_INCLUDE_DIRECTBEAM

C  Fourier component and beam number

      INTEGER          FOURIER_COMPONENT, IBEAM

C  surface factor

      DOUBLE PRECISION SURFACE_FACTOR

C  Surface reflectance function

      EXTERNAL         BVP_SURFACE_SETUP_BEAM

C  output status
C  -------------

      INTEGER          STATUS

C  Local variables
C  ---------------

      CHARACTER*(70)   MAIL, TRACE
      INTEGER          STATUS_SUB

C  This is suitable for Lambertian surfaces only.
C  --Additional setups for the albedo layer.
C    Only required for the Specialise Option 2 case

      IF ( DO_SPECIALIST_OPTION_2 ) THEN
        IF ( ACTIVE_LAYERS(NLAYERS_TEL).EQ.NLAYERS ) THEN
          CALL BVP_SURFACE_SETUP_BEAM
     I       ( DO_INCLUDE_SURFACE,
     I         FOURIER_COMPONENT, IBEAM,
     I         SURFACE_FACTOR )
        ENDIF
      ENDIF

C  --set up Column for solution vector (the "B" as in AX=B)

      CALL BVPTEL_COLUMN_SETUP 
     &  ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM,
     &    IBEAM )

C  --Solve the boundary problem for this Fourier component (back substitution)

      CALL BVPTEL_BACKSUB
     I   ( IBEAM, 
     O     STATUS_SUB )

C  error tracing

      IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
        MAIL  = 'Error return from module BVP_BACKSUB'
        TRACE = 'Call #1 in BVPTEL_SOLUTION_MASTER'
        STATUS = VLIDORT_SERIOUS
        CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
      ENDIF

C  return

      RETURN
      END

C

      SUBROUTINE BVPTEL_COLUMN_SETUP
     &       ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM,
     &         IBEAM )

C  Sets up the telescoped boundary value problem, RHS vector
C    Standard case: Fourier > 0.
C         Suitable for Lambertian surfaces.
C    Surface reflection only appears in the Specialist Option 2 case.

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include file of input variables
C  include file of bookkeeping inputs

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include file of reflectance variables

      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'

C  include files of main model variable (I/O to this module)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'

C  subroutine arguments
C  --------------------

C  inclusion of surface

      LOGICAL          DO_INCLUDE_SURFACE
      LOGICAL          DO_INCLUDE_DIRECTBEAM

C  beam and fourier index

      INTEGER          IBEAM

C  local variables
C  ---------------

      INTEGER          I, I1, N, N1, NS, C0, CM
      INTEGER          IR, O1, IROW

C  Go to special case for only 1 active layer

      IF ( NLAYERS_TEL .EQ. 1 ) GO TO 345

C  zero column vector

      DO I = 1, N_BVTELMATRIX_SIZE
        COLTEL2(I,IBEAM) = ZERO
      ENDDO

C  Upper boundary for first active layer: no downward diffuse radiation

      NS = 1
      N = ACTIVE_LAYERS(NS)
      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        DO O1 = 1, NSTOKES
          IROW = IR + O1
          COLTEL2(IROW,IBEAM)   = - WUPPER(I,O1,N)
        ENDDO
      ENDDO

C  intermediate layer boundaries

      C0 = - NSTKS_NSTRMS
      DO NS = 2, NLAYERS_TEL
        N = ACTIVE_LAYERS(NS)
        N1 = N - 1
        C0 = C0 + NSTKS_NSTRMS_2
        DO I = 1, NSTREAMS_2
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            COLTEL2(CM,IBEAM) = WUPPER(I,O1,N) - WLOWER(I,O1,N1)
          ENDDO
        ENDDO
      ENDDO

C  Lower boundary for last active layer NLAYERS_TEL:
C  No albedo, as this is FOURIER > 0

      NS = NLAYERS_TEL
      N  = ACTIVE_LAYERS(NS)
      C0 = C0 + NSTKS_NSTRMS_2
      IF ( DO_INCLUDE_SURFACE. AND. DO_SPECIALIST_OPTION_2 ) THEN
        IF ( N.EQ.NLAYERS ) THEN
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              COLTEL2(CM,IBEAM) = - WLOWER(I1,O1,N) + R2_BEAM(I,O1)
            ENDDO
          ENDDO
          IF ( DO_INCLUDE_DIRECTBEAM ) THEN
            DO I = 1, NSTREAMS
              IR = NSTOKES*(I-1)
              DO O1 = 1, NSTOKES
                IROW = IR + O1
                CM = C0 + IROW
                COLTEL2(CM,IBEAM) = COLTEL2(CM,IBEAM)
     &                  + DIRECT_BEAM(I,IBEAM,O1)
              ENDDO
            ENDDO
          ENDIF
        ENDIF
      ELSE
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            COLTEL2(CM,IBEAM) = - WLOWER(I1,O1,N)
          ENDDO
        ENDDO
      ENDIF

C  Normal return

      RETURN

C  Continuity point for single layer setup

 345  CONTINUE

C  active layer for telescoped BVP

      N = ACTIVE_LAYERS(1)

C  initialize using the SCOL2 matrix

      DO I = 1, NSTKS_NSTRMS_2
        SCOL2(I,IBEAM) = ZERO
      ENDDO

C  Set up BVP column for the active layer
C  ======================================

C  Upper boundary for layer (downwelling only)
C  -------------------------------------------

      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        DO O1 = 1, NSTOKES
          IROW = IR + O1
          SCOL2(IROW,IBEAM)   = - WUPPER(I,O1,N)
        ENDDO
      ENDDO

C  lower boundary for layer.
C    Surface only required for the Specialist  2 option

      C0 = NSTKS_NSTRMS

      IF ( DO_INCLUDE_SURFACE. AND. DO_SPECIALIST_OPTION_2 ) THEN
        IF ( N.EQ.NLAYERS ) THEN
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              SCOL2(CM,IBEAM) = - WLOWER(I1,O1,N) + R2_BEAM(I,O1)
            ENDDO
          ENDDO
          IF ( DO_INCLUDE_DIRECTBEAM ) THEN
            DO I = 1, NSTREAMS
              IR = NSTOKES*(I-1)
              DO O1 = 1, NSTOKES
                IROW = IR + O1
                CM = C0 + IROW
                SCOL2(CM,IBEAM) = SCOL2(CM,IBEAM)
     &                  + DIRECT_BEAM(I,IBEAM,O1)
              ENDDO
            ENDDO
          ENDIF
        ENDIF
      ELSE
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            SCOL2(CM,IBEAM) = - WLOWER(I1,O1,N)
          ENDDO
        ENDDO
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE BVPTEL_BACKSUB
     I         ( IBEAM,
     O           STATUS )

C  Solves the telescoped boundary value problem.
C    Standard case: Fourier > 0. No surface reflection term
C         Suitable for Lambertian surfaces.

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include file of input variables
C  include file of bookkeeping inputs

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include file of setup variables

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  include files of main model variable (I/O to this module)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'

C  subroutine arguments
C  --------------------

C  Fourier component and beam index

      INTEGER          IBEAM

C  status flag

      INTEGER          STATUS

C  local variables
C  ---------------

      CHARACTER*70     MAIL, TRACE
      INTEGER          K, KO1, K0, K1, K2
      INTEGER          I, I1, N, NS, N1, NAC, O1, IC, ICOW, INFO
      INTEGER          C0, IR, IROW, IROW1, IROW_S, IROW1_S        
      DOUBLE PRECISION SPAR, SHOM, HOM1, HOM2, SHOM_R
      DOUBLE PRECISION SHOM_CR, HOM1CR, HOM2CR
      DOUBLE PRECISION LXR, MXR, LXR_CR, LXR_CI, MXR_CR, MXR_CI
      CHARACTER*3      CI

C  start code
C  ----------

C  Intialize status

      STATUS = VLIDORT_SUCCESS

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
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ENDIF

C  set integration constants for active layers

        C0 = - NSTKS_NSTRMS_2
        DO NS = 1, NLAYERS_TEL
          N = ACTIVE_LAYERS(NS)
          C0 = C0 + NSTKS_NSTRMS_2
          KO1 = K_REAL(N) + 1
          DO K = 1, K_REAL(N)
            IROW = K
            IROW1 = IROW + NSTKS_NSTRMS
            LCON(K,N) = COLTEL2(C0+IROW,IBEAM)
            MCON(K,N) = COLTEL2(C0+IROW1,IBEAM)
          ENDDO
          DO K = 1, K_COMPLEX(N)
            K0 = 2*K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            IROW    = K + K_REAL(N)
            IROW1   = IROW + NSTKS_NSTRMS
            IROW_S  = IROW + K_COMPLEX(N)
            IROW1_S = IROW_S + NSTKS_NSTRMS
            LCON(K1,N) = COLTEL2(C0+IROW,   IBEAM)
            LCON(K2,N) = COLTEL2(C0+IROW_S, IBEAM)
            MCON(K1,N) = COLTEL2(C0+IROW1,  IBEAM)
            MCON(K2,N) = COLTEL2(C0+IROW1_S,IBEAM)
          ENDDO
        ENDDO

c        write(*,*)'tel',fourier,ibeam,lcon(1,7)

C  Solve the boundary problem: Single Layer only
C  =============================================

      ELSE IF ( NLAYERS_TEL .EQ. 1 ) THEN

C  Active layer

        N = ACTIVE_LAYERS(1)
        KO1 = K_REAL(N) + 1

C  LAPACK substitution (DGETRS) using RHS column vector SCOL2

        CALL DGETRS
     &     ( 'N', NSTKS_NSTRMS_2, IBEAM, 
     &        SMAT2, MAXSTRMSTKS_2, SIPIVOT,
     &        SCOL2, MAXSTRMSTKS_2, INFO )

C  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          MAIL='BVP back-substitution (telescoping) (DGETRS)'
          CALL VLIDORT_LAPACK_ERROR ( INFO, 1, MAIL, STATUS )
          IF ( STATUS .EQ. VLIDORT_SERIOUS ) RETURN
        ENDIF

C  set real constants from the solution vector

        DO K = 1, K_REAL(N)
          IROW  = K
          IROW1 = IROW + NSTKS_NSTRMS
          LCON(K,N) = SCOL2(IROW,IBEAM)
          MCON(K,N) = SCOL2(IROW1,IBEAM)
        ENDDO

C  set complex constants from the solution vector

        DO K = 1, K_COMPLEX(N)
          K0 = 2*K - 2
          K1 = KO1 + K0
          K2 = K1  + 1
          IROW    = K + K_REAL(N)
          IROW1   = IROW + NSTKS_NSTRMS
          IROW_S  = K + K_REAL(N) + K_COMPLEX(N)
          IROW1_S = IROW_S + NSTKS_NSTRMS
          LCON(K1,N) = SCOL2(IROW,    IBEAM)
          LCON(K2,N) = SCOL2(IROW_S,  IBEAM)
          MCON(K1,N) = SCOL2(IROW1,   IBEAM)
          MCON(K2,N) = SCOL2(IROW1_S, IBEAM)
        ENDDO

C  end clause for backsubstitution

      ENDIF

C  Set integration constants for non-active layers
C  ===============================================

C  Transmittance layers above active layer
C  ---------------------------------------

C   -- LCON values are zero (no downwelling radiation)
C   -- MCON values propagated upwards from top of first active layer

C  layer immediately above first active layer
C    .... Require solutions at top of active layer

      NAC = ACTIVE_LAYERS(1)
      KO1 = K_REAL(NAC) + 1

      IF ( NAC .GT. 1 ) THEN

        N1 = NAC - 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = ( I1 - 1 ) * NSTOKES
          IC = ( I - 1  ) * NSTOKES
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            ICOW = IC + O1

C  real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(NAC)
             LXR = LCON(K,NAC) * SOLA_XPOS(I1,O1,K,NAC)
             MXR = MCON(K,NAC) * SOLB_XNEG(I1,O1,K,NAC)
             HOM1 = LXR
             HOM2 = MXR * T_DELT_EIGEN(K,NAC)
             SHOM_R = SHOM_R + HOM1 + HOM2
            ENDDO

C  complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(NAC)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              LXR_CR = LCON(K1,NAC) * SOLA_XPOS(I1,O1,K1,NAC) -
     &                 LCON(K2,NAC) * SOLA_XPOS(I1,O1,K2,NAC)
              MXR_CR = MCON(K1,NAC) * SOLB_XNEG(I1,O1,K1,NAC) -
     &                 MCON(K2,NAC) * SOLB_XNEG(I1,O1,K2,NAC)
              MXR_CI = MCON(K1,NAC) * SOLB_XNEG(I1,O1,K2,NAC) +
     &                 MCON(K2,NAC) * SOLB_XNEG(I1,O1,K1,NAC)
              HOM1CR = LXR_CR
              HOM2CR = MXR_CR*T_DELT_EIGEN(K1,NAC)
     &               - MXR_CI*T_DELT_EIGEN(K2,NAC)
              SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
            ENDDO

C  real part and add particular solution
C    ---Sets Real integration constants (no complex ones)

            SHOM = SHOM_R + SHOM_CR
            SPAR = WUPPER(I1,O1,NAC)
            MCON(ICOW,N1) = SPAR + SHOM
            LCON(ICOW,N1) = ZERO

          ENDDO
        ENDDO

      ENDIF

C  other layers to top, just propagate

      DO N = NAC - 2, 1, -1
        N1 = N + 1
        DO I = 1, NSTREAMS
          IR = ( I - 1 ) * NSTOKES
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            LCON(IROW,N) = ZERO
            MCON(IROW,N) = T_DELT_DISORDS(I,N1) * MCON(IROW,N1)
          ENDDO
        ENDDO
      ENDDO

C  Transmittance layers below active layer
C  ---------------------------------------

C  layer immediately below active layer
C    .... Require solutions at bottom of active layer

      NAC = ACTIVE_LAYERS(NLAYERS_TEL)
      KO1 = K_REAL(NAC) + 1
      IF ( NAC .LT. NLAYERS ) THEN
        N1 = NAC + 1
        DO I = 1, NSTREAMS
          IR = ( I - 1 ) * NSTOKES
          DO O1 = 1, NSTOKES
            IROW = IR + O1

C  real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(NAC)
              LXR = LCON(K,NAC) * SOLA_XPOS(I,O1,K,NAC)
              MXR = MCON(K,NAC) * SOLB_XNEG(I,O1,K,NAC)
              HOM1 = LXR * T_DELT_EIGEN(K,NAC)
              HOM2 = MXR
              SHOM_R = SHOM_R + HOM1 + HOM2
            ENDDO

C  complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(NAC)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              LXR_CR = LCON(K1,NAC) * SOLA_XPOS(I,O1,K1,NAC) -
     &                 LCON(K2,NAC) * SOLA_XPOS(I,O1,K2,NAC)
              LXR_CI = LCON(K1,NAC) * SOLA_XPOS(I,O1,K2,NAC) +
     &                 LCON(K2,NAC) * SOLA_XPOS(I,O1,K1,NAC)
              MXR_CR = MCON(K1,NAC) * SOLB_XNEG(I,O1,K1,NAC) -
     &                 MCON(K2,NAC) * SOLB_XNEG(I,O1,K2,NAC)
              HOM1CR = LXR_CR*T_DELT_EIGEN(K1,NAC)
     &                -LXR_CI*T_DELT_EIGEN(K2,NAC)
              HOM2CR = MXR_CR
              SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
            ENDDO

C  real part and add particular solution
C    Note to self. The  MINUS sign appears to be correct !  WHY?

C            SHOM = SHOM_R + SHOM_CR   ! apparently wrong
            SHOM = SHOM_R - SHOM_CR
            SPAR = WLOWER(I,O1,NAC)
            LCON(IROW,N1) = SPAR + SHOM
            MCON(IROW,N1) = ZERO

          ENDDO
        ENDDO

      ENDIF

C  other layers to bottom

      DO N = NAC + 2, NLAYERS
        N1 = N - 1
        DO I = 1, NSTREAMS
          IR = ( I - 1 ) * NSTOKES
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            MCON(IROW,N) = ZERO
            LCON(IROW,N) = T_DELT_DISORDS(I,N1) * LCON(IROW,N1)
          ENDDO
        ENDDO
      ENDDO

C  debug

C        IF ( FOURIER .EQ. 3 ) THEN
C        DO N = 1, NLAYERS
C          DO K = 1, K_REAL(N)
C             WRITE(98,'(3i5,1p2e15.7)')FOURIER,N,K,
C     &        LCON(K,N), MCON(K,N)
C            ENDDO
C        ENDDO
C        ENDIF

C  Finish

      RETURN
      END

C

      INTEGER FUNCTION BMAT_ROWMASK(I,J)

C  Integer function for Band-matrix compression
C   Replaces Array BMAT_ROWMASK(I,J) which was too large.
C   Replacement by J. Kujanpaa, FMI, August 2005.

C  include files

      INCLUDE '../includes/VLIDORT.PARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'

C  Arguments

      INTEGER I, J

C  Assign the function (Keep the hard-wired Stop for now)
C    KALL and NMIN,NMAX are assigned in BVP_MATRIX_INIT

      IF ( (I.GE.NMIN(J)) .AND. (I.LE.NMAX(J)) ) THEN
         BMAT_ROWMASK = KALL + I - J
      ELSE
         BMAT_ROWMASK = 0
         WRITE(*,*) 'ERROR=out of bounds, BMAT_ROWMASK:',I,J
         STOP
      ENDIF

C  Finish

      RETURN
      END

C

      INTEGER FUNCTION BTELMAT_ROWMASK(I,J)

C  Integer function for Band-matrix compression
C   Replaces Array BTELMAT_ROWMASK(I,J) which was too large.
C   Replacement by J. Kujanpaa, FMI, August 2005.

C  include files

      INCLUDE '../includes/VLIDORT.PARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'

C  Arguments

      INTEGER I, J

C  Assign the function (Keep the hard-wired Stop for now)
C    KALL and NMIN,NMAX are assigned in BVP_MATRIX_INIT

      IF ( (I.GE.NMINTEL(J)) .AND. (I.LE.NMAXTEL(J)) ) THEN
         BTELMAT_ROWMASK = KALLTEL + I - J
      ELSE
         BTELMAT_ROWMASK = 0
         WRITE(*,*) 'ERROR=out of bounds, BTELMAT_ROWMASK:',I,J
         STOP
      ENDIF

C  Finish

      RETURN
      END
