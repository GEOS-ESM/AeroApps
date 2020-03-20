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
C #  For linearizations involving atmospheric parameters        #
C #                                                             #
C #            L_BVP_SOLUTION_MASTER                            #
C #            L_BVP_P_COLUMN_SETUP (re-named, Version 2.4)     #
C #            L_BVP_C_COLUMN_SETUP (new for version 2.4)       #
C #            L_BVP_LAMBERTIAN_SURFACE                         #
C #            L_BEAMSOLUTION_P_NEQK (re-named)                 #
C #            L_BEAMSOLUTION_P_NNEK (re-named)                 #
C #            L_BEAMSOLUTION_C_NEQK (new for version 2.4)      #
C #                                                             #
C #  For linearizations with telescoped boundary value problem  #
C #                                                             #
C #            L_BVPTEL_SOLUTION_MASTER                         #
C #            L_BVPTEL_COLUMN_SETUP                            #
C #                                                             #
C #  For linearizations involving surface parameters            #
C #                                                             #
C #            LS_BVP_LAMBERTIAN_SURFACE                        #
C #                                                             #
C ###############################################################

      SUBROUTINE L_BVP_SOLUTION_MASTER
     I       ( DO_INCLUDE_SURFACE,
     I         DO_INCLUDE_DIRECTBEAM,
     I         DO_INCLUDE_THERMEMISS,
     I         VARIATION_INDEX, N_WEIGHTFUNCS,
     I         FOURIER_COMPONENT, IBEAM,
     I         SURFACE_FACTOR,
     I         L_BVP_SURFACE_SETUP,
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

C  include files of linearized input/setup/solution/multiplier variables

      INCLUDE '../includes/VLIDORT_L_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'

C  input arguments
C  ---------------

C  inclusion flags

      LOGICAL          DO_INCLUDE_DIRECTBEAM
      LOGICAL          DO_INCLUDE_SURFACE
      LOGICAL          DO_INCLUDE_THERMEMISS

C  Fourier component and beam number

      INTEGER          FOURIER_COMPONENT, IBEAM

C  surface factor = 1+delta(m,0)

      DOUBLE PRECISION SURFACE_FACTOR

C  Linearization control

      INTEGER          VARIATION_INDEX
      INTEGER          N_WEIGHTFUNCS

C  External function

      EXTERNAL         L_BVP_SURFACE_SETUP

C  output status
C  -------------

      INTEGER          STATUS

C  Local variables
C  ---------------

C  boundary condition flags

      LOGICAL          MODIFIED_BCL3, MODIFIED_BCL4

C  error tracing variables

      INTEGER          STATUS_SUB
      CHARACTER*70     MAIL, TRACE

C  Initialise

      STATUS = VLIDORT_SUCCESS
      
C  Linearization of the regular BVP case
C  =====================================

C  Set up the column vectors for Bulk/profile linearizations
C  ---------------------------------------------------------

C  Bulk: Compute the main column B' where AX = B'

      IF ( DO_COLUMN_LINEARIZATION ) THEN

        CALL L_BVP_C_COLUMN_SETUP
     I     ( DO_INCLUDE_SURFACE,
     I       DO_INCLUDE_DIRECTBEAM,
     I       DO_INCLUDE_THERMEMISS,
     I       N_WEIGHTFUNCS,
     I       FOURIER_COMPONENT,  IBEAM,
     I       SURFACE_FACTOR,
     I       L_BVP_SURFACE_SETUP )

C  Profile: Boundary condition flags for special cases
C  Profile: Compute the main column B' where AX = B'

      ELSE IF ( DO_PROFILE_LINEARIZATION ) THEN

C  Boundary condition flags for special cases

        MODIFIED_BCL3 = ( VARIATION_INDEX .EQ. 1 )
        MODIFIED_BCL4 = ( VARIATION_INDEX .EQ. NLAYERS )

        CALL L_BVP_P_COLUMN_SETUP
     I     ( DO_INCLUDE_SURFACE, 
     I       DO_INCLUDE_DIRECTBEAM,
     I       DO_INCLUDE_THERMEMISS,
     I       MODIFIED_BCL3,      MODIFIED_BCL4,
     I       VARIATION_INDEX,    N_WEIGHTFUNCS,
     I       FOURIER_COMPONENT,  IBEAM,
     I       SURFACE_FACTOR,
     I       L_BVP_SURFACE_SETUP )

      ENDIF

C  Back-substitution

      CALL L_BVP_BACKSUB
     &    ( VARIATION_INDEX,   N_WEIGHTFUNCS, STATUS_SUB )

C  error tracing

      IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
        MAIL  = 'Error return from module L_BVP_BACKSUB'
        TRACE = 'Call in L_BVP_SOLUTION_MASTER'
        STATUS = VLIDORT_SERIOUS
        CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
        RETURN
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE L_BVP_BACKSUB
     I     ( LAYER_TO_VARY, N_LAYER_WFS,
     O       STATUS )

C  Solves the linearized boundary value problem.

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of main model variable (INPUT to this module)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  Include file of linearized solution variables (OUTPUT)

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'

C  subroutine arguments
C  --------------------

C  linearization control

      INTEGER          LAYER_TO_VARY, N_LAYER_WFS

C  status flag

      INTEGER          STATUS

C  local variables
C  ---------------

      INTEGER          C0, LAY, K, Q, KO1, K0, K1, K2
      INTEGER          IROW, IROW1, IROW_S, IROW1_S
      CHARACTER*(70)   MAIL, TRACE
      INTEGER          INFO
      CHARACTER*3      CI, CN

C  Intialize status

      STATUS = VLIDORT_SUCCESS

C  BVP back-substitution: With compression (multilayers)
C  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

C  LAPACK substitution (DGBTRS) using RHS column vector COL2_WF
C  BV solution for perturbed integration constants
C    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS
     &     ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, N_LAYER_WFS,
     &        BANDMAT2, MAXBANDTOTAL, IPIVOT,
     &        COL2_WF, MAXTOTAL, INFO )

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(CN, '(I3)' ) LAYER_TO_VARY
          MAIL  = 'argument i illegal value, for i = '//CI
          TRACE = 'Atmos_Wfs for layer '//CN//
     *                  '; DGBTRS call in L_BVP_BACKSUB'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ENDIF

C  Set Linearized integration constants NCON and PCON, all layers

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

C  Solve the boundary problem: No compression, Single Layer only
C  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

C  LAPACK substitution (DGETRS) using RHS column vector SCOL2

        CALL DGETRS
     &     ( 'N', NTOTAL, N_LAYER_WFS, SMAT2, MAXSTRMSTKS_2, SIPIVOT,
     &        SCOL2_WF, MAXSTRMSTKS_2, INFO )

C  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MAIL  = 'argument i illegal value, for i = '//CI
          TRACE = 'Atmos_Wfs for 1-layer: DGETRS call in L_BVP_BACKSUB'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ENDIF

C  Set linearized integration constants NCON and PCON, 1 layer

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

C  debug------------------------------------------

c       if ( fourier .EQ.3 .and.ibeam.eq.1 ) then
c        do n = 1, nlayers
c          write(76,'(a,2i3,1p2e18.10)')
c     &     'hey',n,variation_index,ncon(3,n,1),pcon(3,n,1)
c        enddo
c       endif

c        if ( do_debug_write ) then
c        if ( fourier.eq.1 ) then
c         Q = 1
c         DO LAY = 1, NLAYERS
c          KO1 = K_REAL(LAY) + 1
c          DO K = 1, K_REAL(LAY)
c           write(51,'(4i3,1p4e20.10)')FOURIER,IBEAM,K,LAY,
c     &                LCON(K,LAY),  MCON(K,LAY),
c     &                NCON(K,LAY,Q),PCON(K,LAY,Q)
c          ENDDO
c          DO K = 1, K_COMPLEX(LAY)
c           K0 = 2*K - 2
c           K1 = KO1 + K0
c           K2 = K1  + 1
c           write(51,'(4i3,1p8e20.10)')FOURIER,IBEAM,K,LAY,
c     &                LCON(K1,LAY), MCON(K1,LAY),
c     &                LCON(K2,LAY), MCON(K2,LAY),
c     &                NCON(K1,LAY,Q),PCON(K1,LAY,Q),
c     &                NCON(K2,LAY,Q),PCON(K2,LAY,Q)
c          ENDDO
c         ENDDO
c        endif
c        ENDIF

C  Associated quantities
C  ---------------------

C  Linearized Real variables

c      DO N = 1, NLAYERS
c       DO I = 1, NSTREAMS_2
c        IR = ( I - 1 ) * NSTOKES
c        DO O1 = 1, NSTOKES
c         IROW = IR + O1
c         DO K = 1, K_REAL(N)
c          DO Q = 1, N_LAYER_WFS
c           NCON_XVEC_R(IROW,K,N,Q) = NCON_R(K,N,Q)*SOLA_XPOS_R(I,O1,K,N)
c           PCON_XVEC_R(IROW,K,N,Q) = PCON_R(K,N,Q)*SOLB_XNEG_R(I,O1,K,N)
c          ENDDO
c         ENDDO
c        ENDDO
c       ENDDO
c      ENDDO

C  Linearized complex variables

c      DO N = 1, NLAYERS
c        DO I = 1, NSTREAMS_2
c          IR = ( I - 1 ) * NSTOKES
c          DO O1 = 1, NSTOKES
c            IROW = IR + O1
c            DO K = 1, K_COMPLEX(N)
c              DO Q = 1, N_LAYER_WFS
c                NCON_XVEC_CR(IROW,K,N,Q) = 
c     &                  NCON_CR(K,N,Q) * SOLA_XPOS_CR(I,O1,K,N) -
c     &                  NCON(K,N,Q) * SOLA_XPOS(I,O1,K,N)
c                NCON_XVEC(IROW,K,N,Q) = 
c     &                  NCON_CR(K,N,Q) * SOLA_XPOS(I,O1,K,N) +
c     &                  NCON(K,N,Q) * SOLA_XPOS_CR(I,O1,K,N)
c                PCON_XVEC_CR(IROW,K,N,Q) =
c     &                  PCON_CR(K,N,Q) * SOLB_XNEG_CR(I,O1,K,N) -
c     &                  PCON(K,N,Q) * SOLB_XNEG(I,O1,K,N)
c                PCON_XVEC(IROW,K,N,Q) =
c     &                  PCON_CR(K,N,Q) * SOLB_XNEG(I,O1,K,N) +
c     &                  PCON(K,N,Q) * SOLB_XNEG_CR(I,O1,K,N)
c              ENDDO
c            ENDDO
c          ENDDO
c        ENDDO
c      ENDDO

C  debug--------------------------------------------
C        IF ( FOURIER .EQ. 3 ) THEN
C        DO N = 1, NLAYERS
C          DO K = 1, K_REAL(N)
C             WRITE(99,'(3i5,1p2e15.7)')FOURIER,N,K,
C     &        NCON_R(K,N,1), PCON_R(K,N,1)
C            ENDDO
C        ENDDO
C        ENDIF

C  finish

      RETURN
      END

C

      SUBROUTINE L_BVP_P_COLUMN_SETUP
     I     ( DO_INCLUDE_SURFACE,
     I       DO_INCLUDE_DIRECTBEAM,
     I       DO_INCLUDE_THERMEMISS,
     I       MODIFIED_BCL3,      MODIFIED_BCL4,
     I       LAYER_TO_VARY,      N_LAYER_WFS,
     I       FOURIER_COMPONENT,  IBEAM,
     I       SURFACE_FACTOR,
     I       L_BVP_SURFACE_SETUP )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and reflectance variables

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'

C  include file of linearized setup variables (input)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'

C  include file of linearized solution variables (output stored here)

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_THERMALSUP.VARS'

C  input arguments
C  ---------------

C  inclusion flags

      LOGICAL          DO_INCLUDE_DIRECTBEAM
      LOGICAL          DO_INCLUDE_SURFACE
      LOGICAL          DO_INCLUDE_THERMEMISS

C  boundary condition flags

      LOGICAL          MODIFIED_BCL3, MODIFIED_BCL4

C  Fourier component and beam number

      INTEGER          FOURIER_COMPONENT, IBEAM

C  surface factor = 1+delta(m,0)

      DOUBLE PRECISION SURFACE_FACTOR

C  External function

      EXTERNAL         L_BVP_SURFACE_SETUP

C  Linearization control

      INTEGER          LAYER_TO_VARY
      INTEGER          N_LAYER_WFS

C  local variables
C  ---------------

      INTEGER          Q, N, N1, I, I1, IR, IROW, CM, C0, O1, NV
      INTEGER          K, KO1, K0, K1, K2
      DOUBLE PRECISION CPOS, CNEG, L_HOM_R, L_HOM_CR, L_BEAM, FAC
      DOUBLE PRECISION T1, T2, T1R, T1I, T2R, T2I, TM

      DOUBLE PRECISION R2_L_BEAM
     &        ( MAXSTREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION R2_L_HOMP
     &        ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )
      DOUBLE PRECISION R2_L_HOMM
     &        ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )

      LOGICAL          REGULAR_BCL3, REGULAR_BCL4

C  initialise
C  ----------

C  zero the results vectors

      DO I = 1, NTOTAL
        DO Q = 1, MAX_ATMOSWFS
          COL2_WF(I,Q) = ZERO
        ENDDO
      ENDDO

C  Layer to vary

      NV = LAYER_TO_VARY

C  Copy already existing thermal linearizations
C    This is a very important zeroing.................!!!!!

      DO I = 1, NSTREAMS_2
        DO O1 = 1, NSTOKES
          DO Q = 1, N_LAYER_WFS
            DO N = 1, NLAYERS
              L_WUPPER(I,O1,N,Q) = ZERO
              L_WLOWER(I,O1,N,Q) = ZERO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      IF ( DO_INCLUDE_THERMEMISS ) THEN
        O1 = 1
        TM = ONE
        IF ( DO_SOLAR_SOURCES ) TM = PI4
        DO I = 1, NSTREAMS_2
          DO Q = 1, N_LAYER_WFS
            L_WUPPER(I,O1,NV,Q) = TM * L_T_WUPPER(I,NV,Q)
            L_WLOWER(I,O1,NV,Q) = TM * L_T_WLOWER(I,NV,Q)
          ENDDO
        ENDDO
      ENDIF
      
C  Get the linearized beam solution for the varying layer

      IF ( DO_SOLAR_SOURCES ) THEN
        CALL L_BEAMSOLUTION_P_NEQK
     I    ( FOURIER_COMPONENT, IBEAM, LAYER_TO_VARY, N_LAYER_WFS )
      ENDIF

C  complete boundary condition flags

      REGULAR_BCL3 = .NOT.MODIFIED_BCL3
      REGULAR_BCL4 = .NOT.MODIFIED_BCL4

C  BCL1 or BCL3M - top of first layer (TOA), UPPER boundary condition
C  ------------------------------------------------------------------

      N = 1

C    If this layer is the one that is varied, use MODIFIED_BCL3 (BCL3M)

      IF ( MODIFIED_BCL3 ) THEN

C  Start  loops

        DO I = 1, NSTREAMS
         IR = NSTOKES*(I-1)
         DO O1 = 1, NSTOKES
          IROW = IR + O1
          DO Q = 1, N_LAYER_WFS

C  beam solution linearization at top of layer

            L_BEAM = - L_WUPPER(I,O1,N,Q)

C  Linearized Real homogeneous solution contributions

            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N)
              CPOS = L_SOLA_XPOS(I,O1,K,N,Q)
              CNEG = T_DELT_EIGEN(K,N)   * L_SOLB_XNEG(I,O1,K,N,Q) + 
     &             L_T_DELT_EIGEN(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

C  Linearized Complex homogeneous solution contributions

            L_HOM_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T1 = L_SOLA_XPOS(I,O1,K1,N,Q) * LCON(K1,N) -
     &             L_SOLA_XPOS(I,O1,K2,N,Q) * LCON(K2,N)
              T2R =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K1,N,Q)
     &             - T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K2,N,Q)
     &             + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N)
     &             - L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)
              T2I =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K2,N,Q)
     &             + T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K1,N,Q)
     &             + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K2,N)
     &             + L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K1,N)
              T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

C  Final contribution

            COL2_WF(IROW,Q) = L_BEAM - L_HOM_R - L_HOM_CR

C  end loops

          ENDDO
         ENDDO
        ENDDO

C  No variation case (BCL1)

      ELSE

        DO I = 1, NSTREAMS
         IR = NSTOKES*(I-1)
         DO O1 = 1, NSTOKES
          IROW = IR + O1
          DO Q = 1, N_LAYER_WFS
            COL2_WF(IROW,Q) = ZERO
          ENDDO
         ENDDO
        ENDDO

      ENDIF

C  BCL2 Intermediate levels between top layer and varying layer
C  ------------------------------------------------------------

C  [not required if top layer is varying, case MODIFIED_BCL3 above]

      IF ( REGULAR_BCL3 ) THEN

C  .. nothing varying in these layers

        DO N = 2, LAYER_TO_VARY - 1
          N1 = N - 1
          C0 = N1*NSTKS_NSTRMS_2 - NSTKS_NSTRMS
          DO I = 1, NSTREAMS_2
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              DO Q = 1, N_LAYER_WFS
                COL2_WF(CM,Q) = ZERO
              ENDDO
            ENDDO
          ENDDO
        ENDDO

      ENDIF

C  BCL3 - regular upper boundary condition for layer that is varying
C  -----------------------------------------------------------------

      IF ( REGULAR_BCL3 ) THEN

C  offsets

        N = LAYER_TO_VARY
        N1  = N - 1
        C0  = N1*NSTKS_NSTRMS_2 - NSTKS_NSTRMS

C  start loops

        DO I = 1, NSTREAMS_2
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            DO Q = 1, N_LAYER_WFS

C  Beam contribution

            L_BEAM  = + L_WUPPER(I,O1,N,Q)

C  Linearized Real homogeneous solution contributions

            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N)
              CPOS = L_SOLA_XPOS(I,O1,K,N,Q)
              CNEG = T_DELT_EIGEN(K,N)   * L_SOLB_XNEG(I,O1,K,N,Q) + 
     &             L_T_DELT_EIGEN(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

C  Linearized Complex homogeneous solution contributions

            L_HOM_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T1 = L_SOLA_XPOS(I,O1,K1,N,Q) * LCON(K1,N) -
     &             L_SOLA_XPOS(I,O1,K2,N,Q) * LCON(K2,N)
              T2R =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K1,N,Q)
     &             - T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K2,N,Q)
     &             + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N)
     &             - L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)
              T2I =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K2,N,Q)
     &             + T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K1,N,Q)
     &             + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K2,N)
     &             + L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K1,N)
              T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

C  Final contribution

            COL2_WF(CM,Q) = L_BEAM + L_HOM_R + L_HOM_CR

C  end loops

          ENDDO
         ENDDO
        ENDDO

C  End BCL3 condition

      ENDIF
 
C  BCL4 - LOWER boundary condition for varying layer
C  -------------------------------------------------

C   special case when layer-to-vary = last (albedo) layer is treated
C   separately below under MODIFIED BCL4.

      IF ( REGULAR_BCL4 ) THEN

C  offsets

        N = LAYER_TO_VARY
        N1  = N + 1
        C0  = N*NSTKS_NSTRMS_2 - NSTKS_NSTRMS

C  Get the linearized beam solution for the next layer

        IF ( DO_SOLAR_SOURCES ) THEN
          CALL L_BEAMSOLUTION_P_NNEK
     I  ( FOURIER_COMPONENT, IBEAM, N1, LAYER_TO_VARY, N_LAYER_WFS )
        ENDIF

C  start loops

        DO I = 1, NSTREAMS_2
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            DO Q = 1, N_LAYER_WFS

C  Beam contributions

            L_BEAM = L_WUPPER(I,O1,N1,Q) - L_WLOWER(I,O1,N,Q)

C  Linearized Real homogeneous solution contributions

            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N)
              CNEG = L_SOLB_XNEG(I,O1,K,N,Q)
              CPOS = T_DELT_EIGEN(K,N)   * L_SOLA_XPOS(I,O1,K,N,Q) + 
     &             L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I,O1,K,N)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

C  Linearized Complex homogeneous solution contributions

            L_HOM_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T2 = L_SOLB_XNEG(I,O1,K1,N,Q) * MCON(K1,N) -
     &             L_SOLB_XNEG(I,O1,K2,N,Q) * MCON(K2,N)
              T1R =  T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I,O1,K1,N,Q)
     &             - T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I,O1,K2,N,Q)
     &             + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I,O1,K1,N)
     &             - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I,O1,K2,N)
              T1I =  T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I,O1,K2,N,Q)
     &             + T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I,O1,K1,N,Q)
     &             + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I,O1,K2,N)
     &             + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I,O1,K1,N)
              T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

C  Final contribution

            COL2_WF(CM,Q) = L_BEAM - L_HOM_R - L_HOM_CR

C  End loops

          ENDDO
         ENDDO
        ENDDO

C  End BCL4 condition

      ENDIF

C  BCL5 - Intermediate boundary conditions between varying layer & final layer
C  ---------------------------------------------------------------------------

      IF ( REGULAR_BCL4 ) THEN

        DO N = LAYER_TO_VARY + 1, NLAYERS - 1

C  offsets

          N1  = N + 1
          C0  = N*NSTKS_NSTRMS_2 - NSTKS_NSTRMS

C  Get the linearized beam solution for the next layer

          IF ( DO_SOLAR_SOURCES ) THEN
            CALL L_BEAMSOLUTION_P_NNEK
     I    ( FOURIER_COMPONENT, IBEAM, N1, LAYER_TO_VARY, N_LAYER_WFS )
          ENDIF

C  .. contributions from beam solution (direct assign). No homog. variation

          DO I = 1, NSTREAMS_2
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              DO Q = 1, N_LAYER_WFS
                L_BEAM = L_WUPPER(I,O1,N1,Q) - L_WLOWER(I,O1,N,Q)
                COL2_WF(CM,Q) = L_BEAM
              ENDDO
            ENDDO
          ENDDO

C  end layer loop

        ENDDO

C  end BCL5 boundary conditions

      ENDIF

C  Final layer - use BCL6 or BCL4M (last layer is varying)
C  -------------------------------------------------------

      N = NLAYERS

C  Modified BCL4M Component loop

      IF ( MODIFIED_BCL4 ) THEN

C  get the linearized downward-reflected term

        CALL L_BVP_SURFACE_SETUP
     I     ( DO_INCLUDE_SURFACE, MODIFIED_BCL4, IBEAM,
     I       FOURIER_COMPONENT, SURFACE_FACTOR, N_LAYER_WFS,
     O       R2_L_BEAM, R2_L_HOMP, R2_L_HOMM )

C  offsets

        C0  = (N-1)*NSTKS_NSTRMS_2 + NSTKS_NSTRMS

C  start loops

        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            DO Q = 1, N_LAYER_WFS

C  Beam contributions

            L_BEAM = L_WLOWER(I1,O1,N,Q) - R2_L_BEAM(I,O1,Q)

C  Linearized Real homogeneous solution contributions

            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N)
              CPOS =  T_DELT_EIGEN(K,N)   * L_SOLA_XPOS(I1,O1,K,N,Q) 
     &            + L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I1,O1,K,N)
              CPOS = CPOS - R2_L_HOMP(I,O1,K,Q)
              CNEG = L_SOLB_XNEG(I1,O1,K,N,Q)
              CNEG = CNEG - R2_L_HOMM(I,O1,K,Q)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

C  Linearized Complex homogeneous solution contributions
C    Bug Fixed 16 December 2005.

            L_HOM_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T1R = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q)
     &            - T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q)
     &            + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
     &            - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
              T1I = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q)
     &            + T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q)
     &            + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
     &            + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
              T1R = T1R - R2_L_HOMP(I,O1,K1,Q)
              T1I = T1I - R2_L_HOMP(I,O1,K2,Q)
              T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
              T2R = L_SOLB_XNEG(I1,O1,K1,N,Q) - R2_L_HOMM(I,O1,K1,Q)
              T2I = L_SOLB_XNEG(I1,O1,K2,N,Q) - R2_L_HOMM(I,O1,K2,Q)
              T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

C  Final contributions

            COL2_WF(CM,Q) = - L_BEAM - L_HOM_R - L_HOM_CR

C  End loops

          ENDDO
         ENDDO
        ENDDO

C  ordinary BCL6 Component loop
 
      ELSE

C  get the linearized downward-reflected term

        CALL L_BVP_SURFACE_SETUP
     I     ( DO_INCLUDE_SURFACE, MODIFIED_BCL4, IBEAM,
     I       FOURIER_COMPONENT, SURFACE_FACTOR, N_LAYER_WFS,
     O       R2_L_BEAM, R2_L_HOMP, R2_L_HOMM )

C  offsets

        C0  = (N-1)*NSTKS_NSTRMS_2 + NSTKS_NSTRMS

C  start loops : beam contributions only

        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            DO Q = 1, N_LAYER_WFS
              L_BEAM = L_WLOWER(I1,O1,N,Q) - R2_L_BEAM(I,O1,Q)
              COL2_WF(CM,Q) = - L_BEAM
            ENDDO
          ENDDO
        ENDDO

C  End lowest layer boundary value linearization

      ENDIF

C  Add direct beam variation to Final boundary
C  -------------------------------------------

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        IF ( DO_INCLUDE_SURFACE ) THEN
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              FAC = - DIRECT_BEAM(I,IBEAM,O1) * 
     &                DELTAU_SLANT(N,LAYER_TO_VARY,IBEAM) 
              DO Q = 1, N_LAYER_WFS
                L_BEAM = L_DELTAU_VERT(Q,LAYER_TO_VARY) * FAC
                COL2_WF(CM,Q) = COL2_WF(CM,Q) + L_BEAM
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF

C  copy the single layer vector

      IF ( NLAYERS .EQ. 1 ) THEN
        DO I = 1, NTOTAL
          DO Q = 1, N_LAYER_WFS
            SCOL2_WF(I,Q) = COL2_WF(I,Q)
          ENDDO
        ENDDO
      ENDIF

C  debug

c      if ( do_debug_write ) then
c        DO N = 1, NTOTAL
c          write(95,'(4i4,1p4e17.9)')
c     &      FOURIER_COMPONENT,IBEAM,LAYER_TO_VARY,N,
c     &                 COL2_WF(N,1),COL2_WF(N,2)
c        ENDDO
c      ENDIF
c       pause

C  finish

      RETURN
      END

C

      SUBROUTINE L_BVP_C_COLUMN_SETUP
     I     ( DO_INCLUDE_SURFACE,
     I       DO_INCLUDE_DIRECTBEAM,
     I       DO_INCLUDE_THERMEMISS,
     I       N_WEIGHTFUNCS,
     I       FOURIER_COMPONENT,  IBEAM,
     I       SURFACE_FACTOR,
     I       L_BVP_SURFACE_SETUP )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and reflectance variables

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'

C  include file of linearized setup variables (input)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'

C  include file of linearized solution variables (output stored here)

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_THERMALSUP.VARS'

C  input arguments
C  ---------------

C  inclusion flags

      LOGICAL          DO_INCLUDE_THERMEMISS
      LOGICAL          DO_INCLUDE_DIRECTBEAM
      LOGICAL          DO_INCLUDE_SURFACE

C  Fourier component and beam number

      INTEGER          FOURIER_COMPONENT, IBEAM

C  surface factor = 1+delta(m,0)

      DOUBLE PRECISION SURFACE_FACTOR

C  External function

      EXTERNAL         L_BVP_SURFACE_SETUP

C  Linearization control

      INTEGER          N_WEIGHTFUNCS

C  local variables
C  ---------------

      INTEGER          Q, N, N1, I, I1, IR, IROW, CM, C0, O1, NV
      INTEGER          K, KO1, K0, K1, K2
      DOUBLE PRECISION CPOS, CNEG, L_HOM_R, L_HOM_CR, L_BEAM, FAC
      DOUBLE PRECISION T1, T2, T1R, T1I, T2R, T2I, FAC3, TM
      DOUBLE PRECISION L_HOM_U_R  , L_HOM_D_R
      DOUBLE PRECISION L_HOM_U_CR , L_HOM_D_CR
      LOGICAL          MODIFIED_BOUNDARY

      DOUBLE PRECISION R2_L_BEAM
     &        ( MAXSTREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION R2_L_HOMP
     &        ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )
      DOUBLE PRECISION R2_L_HOMM
     &        ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )

C  initialise
C  ----------

C  zero the results vectors

      DO I = 1, NTOTAL
        DO Q = 1, MAX_ATMOSWFS
          COL2_WF(I,Q) = ZERO
        ENDDO
      ENDDO

C  Copy already existing thermal linearizations
C    This is a very important zeroing.................!!!!!

      DO NV = 1, NLAYERS
        IF ( DO_INCLUDE_THERMEMISS ) THEN
          TM = ONE
          IF ( DO_SOLAR_SOURCES ) TM = PI4
          DO I = 1, NSTREAMS_2
            DO Q = 1, N_WEIGHTFUNCS
              L_WUPPER(I,1,NV,Q) = TM * L_T_WUPPER(I,NV,Q)
              L_WLOWER(I,1,NV,Q) = TM * L_T_WLOWER(I,NV,Q)
              DO O1 = 2, NSTOKES
                L_WUPPER(I,O1,N,Q) = ZERO
                L_WLOWER(I,O1,N,Q) = ZERO
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO I = 1, NSTREAMS_2
            DO O1 = 1, NSTOKES
              DO Q = 1, N_WEIGHTFUNCS
                L_WUPPER(I,O1,NV,Q) = ZERO
                L_WLOWER(I,O1,NV,Q) = ZERO
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDDO

C  Top of first layer (TOA), UPPER boundary condition
C  --------------------------------------------------

      N = 1

C  Get the linearized beam solution for the first layer

      IF ( DO_SOLAR_SOURCES ) THEN
        CALL L_BEAMSOLUTION_C_NEQK
     I     ( FOURIER_COMPONENT, IBEAM, N, N_WEIGHTFUNCS )
      ENDIF

C  .. contribution WVAR from beam solution variations
C  .. contribution HVAR homogeneous (eigenvalue) solution variations

      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        DO O1 = 1, NSTOKES
          IROW = IR + O1
          DO Q = 1, N_WEIGHTFUNCS

C  beam solution linearization at top of layer

            L_BEAM = - L_WUPPER(I,O1,N,Q)

C  Linearized Real homogeneous solution contributions

            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N)
              CPOS = L_SOLA_XPOS(I,O1,K,N,Q)
              CNEG = T_DELT_EIGEN(K,N)   * L_SOLB_XNEG(I,O1,K,N,Q) + 
     &             L_T_DELT_EIGEN(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

C  Linearized Complex homogeneous solution contributions

            L_HOM_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T1 = L_SOLA_XPOS(I,O1,K1,N,Q) * LCON(K1,N) -
     &             L_SOLA_XPOS(I,O1,K2,N,Q) * LCON(K2,N)
              T2R =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K1,N,Q)
     &             - T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K2,N,Q)
     &             + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N)
     &             - L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)
              T2I =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K2,N,Q)
     &             + T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K1,N,Q)
     &             + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K2,N)
     &             + L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K1,N)
              T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

C  Final contribution

            COL2_WF(IROW,Q) = L_BEAM - L_HOM_R - L_HOM_CR
          ENDDO
        ENDDO
      ENDDO

C  Intermediate boundary conditions
C  --------------------------------

      DO N = 1, NLAYERS - 1

C  N1 is the layer below, C0 is the offset

        N1 = N + 1
        C0 = N*NSTKS_NSTRMS_2 - NSTKS_NSTRMS

C  Get the linearized beam solution for the next layer

        IF ( DO_SOLAR_SOURCES ) THEN
          CALL L_BEAMSOLUTION_C_NEQK
     I   ( FOURIER_COMPONENT, IBEAM, N1, N_WEIGHTFUNCS )
        ENDIF

C  .. 2 contributions to L_BEAM, from variations L_WUPPER L_WLOWER 
C  .. 2 contributions to L_HOM,  from variations above and below

        DO I = 1, NSTREAMS_2
         IR = NSTOKES*(I-1)
         DO O1 = 1, NSTOKES
          IROW = IR + O1
          CM = C0 + IROW
          DO Q = 1, N_WEIGHTFUNCS

C  Beam contributions

            L_BEAM  = + L_WUPPER(I,O1,N1,Q) - L_WLOWER(I,O1,N,Q)  

C  Linearized Real homogeneous solution contributions above

            L_HOM_U_R = ZERO
            DO K = 1, K_REAL(N1)
              CPOS = L_SOLA_XPOS(I,O1,K,N1,Q)
              CNEG = T_DELT_EIGEN(K,N1)   * L_SOLB_XNEG(I,O1,K,N1,Q) + 
     &             L_T_DELT_EIGEN(K,N1,Q) *   SOLB_XNEG(I,O1,K,N1)
              T1 = LCON(K,N1) * CPOS
              T2 = MCON(K,N1) * CNEG
              L_HOM_U_R = L_HOM_U_R + T1 + T2
            ENDDO

C  Linearized Real homogeneous solution contributions below

            L_HOM_D_R = ZERO
            DO K = 1, K_REAL(N)
              CNEG = L_SOLB_XNEG(I,O1,K,N,Q)
              CPOS = T_DELT_EIGEN(K,N)   * L_SOLA_XPOS(I,O1,K,N,Q) + 
     &             L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I,O1,K,N)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_D_R = L_HOM_D_R + T1 + T2
            ENDDO

C  Linearized Complex homogeneous solution contributions above

            L_HOM_U_CR  = ZERO
            KO1 = K_REAL(N1) + 1
            DO K = 1, K_COMPLEX(N1)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T1 = L_SOLA_XPOS(I,O1,K1,N1,Q) * LCON(K1,N1) -
     &             L_SOLA_XPOS(I,O1,K2,N1,Q) * LCON(K2,N1)
              T2R =  T_DELT_EIGEN(K1,N1)   * L_SOLB_XNEG(I,O1,K1,N1,Q)
     &             - T_DELT_EIGEN(K2,N1)   * L_SOLB_XNEG(I,O1,K2,N1,Q)
     &             + L_T_DELT_EIGEN(K1,N1,Q) *   SOLB_XNEG(I,O1,K1,N1)
     &             - L_T_DELT_EIGEN(K2,N1,Q) *   SOLB_XNEG(I,O1,K2,N1)
              T2I =  T_DELT_EIGEN(K1,N1)   * L_SOLB_XNEG(I,O1,K2,N1,Q)
     &             + T_DELT_EIGEN(K2,N1)   * L_SOLB_XNEG(I,O1,K1,N1,Q)
     &             + L_T_DELT_EIGEN(K1,N1,Q) *   SOLB_XNEG(I,O1,K2,N1)
     &             + L_T_DELT_EIGEN(K2,N1,Q) *   SOLB_XNEG(I,O1,K1,N1)
              T2 =  T2R * MCON(K1,N1) - T2I * MCON(K2,N1)
              L_HOM_U_CR = L_HOM_U_CR + T1 + T2
            ENDDO

C  Linearized Complex homogeneous solution contributions below

            L_HOM_D_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T2 = L_SOLB_XNEG(I,O1,K1,N,Q) * MCON(K1,N) -
     &             L_SOLB_XNEG(I,O1,K2,N,Q) * MCON(K2,N)
              T1R =  T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I,O1,K1,N,Q)
     &             - T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I,O1,K2,N,Q)
     &             + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I,O1,K1,N)
     &             - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I,O1,K2,N)
              T1I =  T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I,O1,K2,N,Q)
     &             + T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I,O1,K1,N,Q)
     &             + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I,O1,K2,N)
     &             + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I,O1,K1,N)
              T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
              L_HOM_D_CR = L_HOM_D_CR + T1 + T2
            ENDDO

C  Final contribution

            L_HOM_R  =  L_HOM_U_R  - L_HOM_D_R
            L_HOM_CR =  L_HOM_U_CR - L_HOM_D_CR
            COL2_WF(CM,Q) = L_BEAM + L_HOM_R + L_HOM_CR

c            if (fourier_component.lt.5.and.n.eq.10)
c     &          write(40,'(3i5,1p2e24.12)')fourier_component,o1,irow,
c     &               L_WUPPER(I,O1,N,Q),WUPPER(I,O1,N)

          ENDDO
         ENDDO
        ENDDO

C  End layer

      ENDDO

C  LOWER layer 
C  -----------

      N = NLAYERS
      MODIFIED_BOUNDARY = .TRUE.

C  get the linearized downward-reflected term

      CALL L_BVP_SURFACE_SETUP
     I     ( DO_INCLUDE_SURFACE, MODIFIED_BOUNDARY, IBEAM,
     I       FOURIER_COMPONENT, SURFACE_FACTOR, N_WEIGHTFUNCS,
     O       R2_L_BEAM, R2_L_HOMP, R2_L_HOMM )

C  offsets

      C0  = (N-1)*NSTKS_NSTRMS_2 + NSTKS_NSTRMS

C  start loops

      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        I1 = I + NSTREAMS
        DO O1 = 1, NSTOKES
          IROW = IR + O1
          CM = C0 + IROW
          DO Q = 1, N_WEIGHTFUNCS

C  Beam contributions

            L_BEAM = L_WLOWER(I1,O1,N,Q) - R2_L_BEAM(I,O1,Q)

C  Linearized Real homogeneous solution contributions

            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N)
              CPOS =  T_DELT_EIGEN(K,N)   * L_SOLA_XPOS(I1,O1,K,N,Q) 
     &            + L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I1,O1,K,N)
              CPOS = CPOS - R2_L_HOMP(I,O1,K,Q)
              CNEG = L_SOLB_XNEG(I1,O1,K,N,Q)
              CNEG = CNEG - R2_L_HOMM(I,O1,K,Q)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

C  Linearized Complex homogeneous solution contributions
C    Bug Fixed 16 December 2005.

            L_HOM_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T1R = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q)
     &            - T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q)
     &            + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
     &            - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
              T1I = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q)
     &            + T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q)
     &            + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
     &            + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
              T1R = T1R - R2_L_HOMP(I,O1,K1,Q)
              T1I = T1I - R2_L_HOMP(I,O1,K2,Q)
              T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
              T2R = L_SOLB_XNEG(I1,O1,K1,N,Q) - R2_L_HOMM(I,O1,K1,Q)
              T2I = L_SOLB_XNEG(I1,O1,K2,N,Q) - R2_L_HOMM(I,O1,K2,Q)
              T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

C  Final contributions

            COL2_WF(CM,Q) = - L_BEAM - L_HOM_R - L_HOM_CR

C  End loops

          ENDDO
        ENDDO
      ENDDO

C  Add direct beam variation to Final boundary
C  -------------------------------------------

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        IF ( DO_INCLUDE_SURFACE ) THEN
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              FAC = - DIRECT_BEAM(I,IBEAM,O1) 
              DO Q = 1, N_WEIGHTFUNCS
                L_BEAM = ZERO
                DO K = 1, NLAYERS
                  FAC3 = FAC * DELTAU_SLANT(N,K,IBEAM) 
                  L_BEAM = L_BEAM + L_DELTAU_VERT(Q,K) * FAC3
                ENDDO
                COL2_WF(CM,Q) = COL2_WF(CM,Q) + L_BEAM
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF

C  copy the single layer vector

      IF ( NLAYERS .EQ. 1 ) THEN
        DO I = 1, NTOTAL
          DO Q = 1, N_WEIGHTFUNCS
            SCOL2_WF(I,Q) = COL2_WF(I,Q)
          ENDDO
        ENDDO
      ENDIF

C  debug

c      if ( do_debug_write ) then
c        DO N = 1, NTOTAL
c          write(95,'(3i4,1p4e17.9)')
c     &      FOURIER_COMPONENT,IBEAM,N,
c     &                 COL2_WF(N,1),COL2_WF(N,2)
c        ENDDO
c      ENDIF
c       pause'f95'

C  finish

      RETURN
      END

C

      SUBROUTINE L_BVP_LAMBERTIAN_SURFACE
     I     ( DO_INCLUDE_SURFACE, MODIFIED_BCL4, IBEAM,
     I       FOURIER_COMPONENT, SURFACE_FACTOR, N_LAYER_WFS,
     O       R2_L_BEAM, R2_L_HOMP, R2_L_HOMM )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and reflectance variables

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'

C  include file of linearized setup variables (input)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'

C  include file of linearized solution variables (output stored here)

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'

C  input arguments
C  ---------------

C  Beam number

      INTEGER          IBEAM

C  Number of weighting functions

      INTEGER          N_LAYER_WFS

C  Fourier component and beam index

      INTEGER          FOURIER_COMPONENT

C  Flag for type of boundary condition

      LOGICAL          MODIFIED_BCL4

C  overall surface flag and surface factor

      LOGICAL          DO_INCLUDE_SURFACE
      DOUBLE PRECISION SURFACE_FACTOR

C  Output arguments
C  ----------------

      DOUBLE PRECISION R2_L_BEAM
     &        ( MAXSTREAMS, MAXSTOKES, MAX_ATMOSWFS )

      DOUBLE PRECISION R2_L_HOMP
     &        ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )
      DOUBLE PRECISION R2_L_HOMM
     &        ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )

C  Local variables
C  ---------------

      DOUBLE PRECISION PV_W ( MAXSTREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION HV_P ( MAXSTREAMS, MAXEVALUES, MAX_ATMOSWFS )
      DOUBLE PRECISION HV_M ( MAXSTREAMS, MAXEVALUES, MAX_ATMOSWFS )

      INTEGER          I, J, O1, Q, N, IB
      INTEGER          K, KO1, K0, K1, K2
      DOUBLE PRECISION H1R, H1I, REFL_P, REFL_B, REFL_M, KMULT
      DOUBLE PRECISION H1, H2, H1_CR, H2_CR, H1_CI, H2_CI

C  Initial section
C  ---------------

      IB = IBEAM

C  Always zero the result to start

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

C  Return if no albedo

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

C  Return if Fourier component > 0

      IF ( FOURIER_COMPONENT .GT. 0 ) RETURN

C  Set up Auxiliary arrays
C  -----------------------

C  Only require (1,1) component

      O1 = 1

C  Last layer

      N = NLAYERS

C  Need Beam arrays regardless
C   Always something from layers above the PBL and also from PBL

      DO J = 1, NSTREAMS
        DO Q = 1, N_LAYER_WFS
          PV_W(J,Q) = L_WLOWER(J,O1,N,Q) * QUAD_STRMWTS(J)
        ENDDO
      ENDDO

C   Modified boundary condition
C     This applies when PBL is varying layer.
C       Require homogeneous solutions linearizations

      IF ( MODIFIED_BCL4 ) THEN

C  start loops

        DO J = 1, NSTREAMS
          DO Q = 1, N_LAYER_WFS

C  real homogeneous solution contributions

           DO K = 1, K_REAL(N)
            H1 = L_SOLA_XPOS(J,O1,K,N,Q) *   T_DELT_EIGEN(K,N) +
     &             SOLA_XPOS(J,O1,K,N)   * L_T_DELT_EIGEN(K,N,Q)
            H2 = L_SOLB_XNEG(J,O1,K,N,Q)
            HV_P(J,K,Q) = QUAD_STRMWTS(J)*H1
            HV_M(J,K,Q) = QUAD_STRMWTS(J)*H2
           ENDDO

C  Complex homogeneous solution contributions

           KO1 = K_REAL(N) + 1
           DO K = 1, K_COMPLEX(N)
            K0 = 2*K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            H1R = L_SOLA_XPOS(J,O1,K1,N,Q) *   T_DELT_EIGEN(K1,N)
     &          - L_SOLA_XPOS(J,O1,K2,N,Q) *   T_DELT_EIGEN(K2,N)
     &          +   SOLA_XPOS(J,O1,K1,N)   * L_T_DELT_EIGEN(K1,N,Q)
     &          -   SOLA_XPOS(J,O1,K2,N)   * L_T_DELT_EIGEN(K2,N,Q)
            H1I = L_SOLA_XPOS(J,O1,K1,N,Q) *   T_DELT_EIGEN(K2,N)
     &          + L_SOLA_XPOS(J,O1,K2,N,Q) *   T_DELT_EIGEN(K1,N)
     &          +   SOLA_XPOS(J,O1,K1,N)   * L_T_DELT_EIGEN(K2,N,Q)
     &          +   SOLA_XPOS(J,O1,K2,N)   * L_T_DELT_EIGEN(K1,N,Q)
            HV_P(J,K1,Q) = QUAD_STRMWTS(J)* H1R
            HV_P(J,K2,Q) = QUAD_STRMWTS(J)* H1I
            HV_M(J,K1,Q) = QUAD_STRMWTS(J)* L_SOLB_XNEG(J,O1,K1,N,Q)
            HV_M(J,K2,Q) = QUAD_STRMWTS(J)* L_SOLB_XNEG(J,O1,K2,N,Q)
           ENDDO

C  End loops

          ENDDO
        ENDDO

C  End modified BCL4 condition

      ENDIF

C  reflection

      KMULT = SURFACE_FACTOR * LAMBERTIAN_ALBEDO

C  Integrated Downward reflection (Calculation)
C  --------------------------------------------

      DO Q = 1, N_LAYER_WFS

C  Particular solution (only for the first Stokes component)

       REFL_B = ZERO
       DO J = 1, NSTREAMS
         REFL_B = REFL_B + PV_W(J,Q)
       ENDDO
       REFL_B = REFL_B * KMULT
       DO I = 1, NSTREAMS
         R2_L_BEAM(I,O1,Q) = REFL_B
       ENDDO

C  Homogeneous solutions for the modified condition

       IF ( MODIFIED_BCL4 ) THEN

C  Homogeneous real solutions

         DO K = 1, K_REAL(NLAYERS)
          REFL_P = ZERO
          REFL_M = ZERO
          DO J = 1, NSTREAMS
           REFL_P = REFL_P + HV_P(J,K,Q)
           REFL_M = REFL_M + HV_M(J,K,Q)
          ENDDO
          REFL_P = REFL_P * KMULT
          REFL_M = REFL_M * KMULT
          DO I = 1, NSTREAMS
           R2_L_HOMP(I,O1,K,Q) = REFL_P
           R2_L_HOMM(I,O1,K,Q) = REFL_M
          ENDDO
         ENDDO

C  Homogeneous complex solutions

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
           H1_CR = H1_CR + HV_P(J,K1,Q)
           H1_CI = H1_CI + HV_P(J,K2,Q)
           H2_CR = H2_CR + HV_M(J,K1,Q)
           H2_CI = H2_CI + HV_M(J,K2,Q)
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

C  End modified boundary condition clause

        ENDIF

C  end parameter loop

      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE L_BEAMSOLUTION_P_NEQK
     I           ( FOURIER, IB, N, N_LAYER_WFS )

C  Linearization of beam particular integral in one layer.

C  In this module, this is the Layer that contains the variation.
C                   N = LAYER_TO_VARY

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and multiplier variables

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  include file of linearized setup variables (input)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'

C  include file of linearized solution variables (output stored here)

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'

C  include file of linearized Multiplier coefficients (output)

      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Beam index input

      INTEGER          IB

C  Given layer index (input)
C    This is also LAYER_TO_VARY

      INTEGER          N

C  number of varying parameters (input)

      INTEGER          N_LAYER_WFS

C  Fourier component

      INTEGER          FOURIER

C  Local variables
C  ---------------

      INTEGER          I, O1, Q
      DOUBLE PRECISION CONST, WDEL, VAR1
      DOUBLE PRECISION LBSOL(MAXSTREAMS_2,MAXSTOKES,MAX_ATMOSWFS,2)

C  No linearized particular solution beyond the cutoff layer. ALSO--
C  Nothing if the solution saving mode is on and layer is inactive
C    Exit (Solutions have already been zeroed).

      IF (.NOT.DO_LAYER_SCATTERING(FOURIER,N)
     &       .OR. (N .GT.LAYER_PIS_CUTOFF(IB))) THEN
        RETURN
      ENDIF

C  Classical solution
C  ==================

C  Very simple, same code for all situations

      IF ( DO_CLASSICAL_SOLUTION ) THEN
        CONST   = INITIAL_TRANS(N,IB)
        WDEL    = T_DELT_MUBAR(N,IB)
        DO Q = 1, N_LAYER_WFS
          VAR1 = L_T_DELT_MUBAR(N,N,IB,Q) * CONST
          DO I = 1, NSTREAMS_2
            DO O1 = 1, NSTOKES
              LBSOL(I,O1,Q,1) = CONST * L_BVEC(I,O1,N,N,Q)
              LBSOL(I,O1,Q,2) = WDEL  * LBSOL(I,O1,Q,1)
     &                              + VAR1 * BVEC(I,O1,N)
            ENDDO
          ENDDO
        ENDDO
      ELSE
C  GREENS FUNCTION PLACEHOLDER
      ENDIF

C Add to existing solution

      DO Q = 1, N_LAYER_WFS
        DO I = 1, NSTREAMS_2
          DO O1 = 1, NSTOKES
            L_WUPPER(I,O1,N,Q) = L_WUPPER(I,O1,N,Q) + LBSOL(I,O1,Q,1)
            L_WLOWER(I,O1,N,Q) = L_WLOWER(I,O1,N,Q) + LBSOL(I,O1,Q,2)
          ENDDO
        ENDDO
      ENDDO
    
C  Finish

      RETURN
      END

C

      SUBROUTINE L_BEAMSOLUTION_P_NNEK
     I     ( FOURIER, IB, N, NV, NV_PARAMETERS )

C  Linearization of beam particular integral in one layer only.
C  This is for a layer N not equal to varying layer NV.

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and multiplier variables

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  include file of linearized setup variables (input)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'

C  include file of linearized solution variables (output stored here)

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'

C  include file of linearized Multiplier coefficients (output)

      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Beam index

      INTEGER          IB

C  Given layer indices (inputs)

      INTEGER          N, NV

C  number of varying parameters (input)

      INTEGER          NV_PARAMETERS

C  Fourier component

      INTEGER          FOURIER

C  Local variables
C  ---------------

      INTEGER          I, O1, Q
      DOUBLE PRECISION CONST, WDEL, TRANS2, VAR1, VAR2, VAR_U

C  No linearized particular solution beyond the cutoff layer. ALSO--
C  Nothing if the solution saving mode is on and layer is inactive

      IF (.NOT.DO_LAYER_SCATTERING(FOURIER,N)
     &       .OR. (N .GT.LAYER_PIS_CUTOFF(IB))) THEN
        RETURN
      ENDIF

C  Classical solution
C  ==================

      IF ( DO_CLASSICAL_SOLUTION ) THEN

C  initialise

        CONST   = INITIAL_TRANS(N,IB)
        WDEL    = T_DELT_MUBAR(N,IB)
        TRANS2  = CONST * WDEL

C  Distinguish two cases
C  ..(a) quasi-spherical for n > 1 (only gets done for this case anyway)
C  ..(b) plane-parallel and qs for n = 1

        IF ( .NOT. DO_PLANE_PARALLEL  ) THEN
          DO Q = 1, NV_PARAMETERS
            VAR1 = L_T_DELT_MUBAR(N,NV,IB,Q) * CONST
            VAR2 = L_INITIAL_TRANS(N,NV,IB,Q)
            DO I = 1, NSTREAMS_2
              DO O1 = 1, NSTOKES
                VAR_U = VAR2 * BVEC(I,O1,N) + L_BVEC(I,O1,N,NV,Q)
                L_WUPPER(I,O1,N,Q) = CONST  * VAR_U
                L_WLOWER(I,O1,N,Q) = TRANS2 * VAR_U
     &                                 + VAR1 * BVEC(I,O1,N)
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO Q = 1, NV_PARAMETERS
            VAR1 = L_INITIAL_TRANS(N,NV,IB,Q) * CONST
            DO I = 1, NSTREAMS_2
              DO O1 = 1, NSTOKES
                L_WUPPER(I,O1,N,Q) = VAR1 * BVEC(I,O1,N)
                L_WLOWER(I,O1,N,Q) = L_WUPPER(I,O1,N,Q) * WDEL
              ENDDO
            ENDDO
          ENDDO
        ENDIF

      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE L_BEAMSOLUTION_C_NEQK
     I           ( FOURIER, IB, N, N_LAYER_WFS )

C  Linearization of beam particular integral in layer N, due to column.
C   This is the bulk property linearization

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and multiplier variables

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  include file of linearized setup variables (input)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'

C  include file of linearized solution variables (output stored here)

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'

C  include file of linearized Multiplier coefficients (output)

      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Beam index input

      INTEGER          IB

C  Given layer index (input)
C    This is also LAYER_TO_VARY

      INTEGER          N

C  number of varying parameters (input)

      INTEGER          N_LAYER_WFS

C  Fourier component

      INTEGER          FOURIER

C  Local variables
C  ---------------

      INTEGER          I, O1, Q, LVARY
      DOUBLE PRECISION CONST, WDEL, VAR1, VAR2, VAR_U, TRANS2
      DOUBLE PRECISION LBSOL(MAXSTREAMS_2,MAXSTOKES,MAX_ATMOSWFS,2)

C  No linearized particular solution beyond the cutoff layer. ALSO--
C  Nothing if the solution saving mode is on and layer is inactive
C   ALSO - Nothing if layer has no scattering (Do not need solution saving)

      IF (.NOT.DO_LAYER_SCATTERING(FOURIER,N)
     &       .OR. (N .GT.LAYER_PIS_CUTOFF(IB))) THEN
        RETURN
      ENDIF

C  Varying index is 0 in the holding arrays

      LVARY = 0

C  Classical solution
C  ==================

C  Very simple, same code for all situations

      IF ( DO_CLASSICAL_SOLUTION ) THEN
        CONST   = INITIAL_TRANS(N,IB)
        WDEL    = T_DELT_MUBAR(N,IB)
        TRANS2  = CONST * WDEL
        DO Q = 1, N_LAYER_WFS
          VAR1 = L_T_DELT_MUBAR (N,LVARY,IB,Q) * CONST
          VAR2 = L_INITIAL_TRANS(N,LVARY,IB,Q)
          DO I = 1, NSTREAMS_2
            DO O1 = 1, NSTOKES
	      VAR_U = VAR2 * BVEC(I,O1,N) + L_BVEC(I,O1,N,LVARY,Q)
              LBSOL(I,O1,Q,1) = CONST  * VAR_U
              LBSOL(I,O1,Q,2) = TRANS2 * VAR_U + VAR1 * BVEC(I,O1,N)
            ENDDO
          ENDDO
        ENDDO
      ELSE
C  Greens function placeholder
      ENDIF

C Add to existing solution

      DO Q = 1, N_LAYER_WFS
        DO I = 1, NSTREAMS_2
          DO O1 = 1, NSTOKES
            L_WUPPER(I,O1,N,Q) = L_WUPPER(I,O1,N,Q) + LBSOL(I,O1,Q,1)
            L_WLOWER(I,O1,N,Q) = L_WLOWER(I,O1,N,Q) + LBSOL(I,O1,Q,2)
          ENDDO
        ENDDO
      ENDDO
    
C  Finish

      RETURN
      END

C

      SUBROUTINE L_BVPTEL_SOLUTION_MASTER
     I       ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM,
     I         LAYER_TO_VARY, N_LAYER_WFS,
     I         FOURIER_COMPONENT, IBEAM, SURFACE_FACTOR,
     I         L_BVP_SURFACE_SETUP,
     O         STATUS )

C  include files
C  -------------

C  include file of dimensiopns and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of solution and setup variables (inputs to module)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  include files of linearized input/setup variables (inputs to module)

      INCLUDE '../includes/VLIDORT_L_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'

C  include file of linearized solution variables (input and output)

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'

C  input arguments
C  ---------------

C  inclusion flags

      LOGICAL          DO_INCLUDE_SURFACE
      LOGICAL          DO_INCLUDE_DIRECTBEAM

C  Fourier component and beam number

      INTEGER          FOURIER_COMPONENT, IBEAM

C  surface factor = 1+delta(m,0)

      DOUBLE PRECISION SURFACE_FACTOR

C  Linearization control

      INTEGER          LAYER_TO_VARY
      INTEGER          N_LAYER_WFS

C  External function

      EXTERNAL         L_BVP_SURFACE_SETUP

C  output status
C  -------------

      INTEGER          STATUS

C  Local variables
C  ---------------

C  error tracing variables

      INTEGER          INFO
      CHARACTER*70     MAIL, TRACE
      CHARACTER*3      CI

C  Other local help variables 

      INTEGER          I, I1, N, N1, NAC, O1, IC, ICOW, Q
      INTEGER          K, KO1, K0, K1, K2, C0, NS
      INTEGER          IR, IROW, IROW1, IROW_S, IROW1_S        
      DOUBLE PRECISION SPAR, SHOM, L_HOM1, L_HOM2, SHOM_R
      DOUBLE PRECISION SHOM_CR, L_HOM1CR, L_HOM2CR
      DOUBLE PRECISION LXR, MXR, NXR, PXR, LLXR, MLXR
      DOUBLE PRECISION LXR1, MXR1, NXR1, PXR1, LLXR1, MLXR1
      DOUBLE PRECISION LXR2, MXR2, NXR2, PXR2, LLXR2, MLXR2

C  Initialise

      STATUS = VLIDORT_SUCCESS
      Q = 0

C  Bulk: Compute the main column B' where AX = B'

      IF ( DO_COLUMN_LINEARIZATION ) THEN

        CALL L_BVPTEL_C_COLUMN_SETUP
     I     ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM,
     I       N_LAYER_WFS,
     I       FOURIER_COMPONENT,  IBEAM,
     I       SURFACE_FACTOR,
     I       L_BVP_SURFACE_SETUP )

C  Profile: Boundary condition flags for special cases
C  Profile: Compute the main column B' where AX = B'

      ELSE IF ( DO_PROFILE_LINEARIZATION ) THEN

        CALL L_BVPTEL_P_COLUMN_SETUP
     I     ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM,
     I       LAYER_TO_VARY,      N_LAYER_WFS,
     I       FOURIER_COMPONENT,  IBEAM,
     I       SURFACE_FACTOR,
     I       L_BVP_SURFACE_SETUP )

      ENDIF

C  Solve linearized BVP: Several Active layers
C  ===========================================

      IF ( NLAYERS_TEL .GT. 1 ) THEN

C  BV solution for linearized integration constants
C    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS
     &     ( 'n', N_BVTELMATRIX_SIZE, N_BVTELMATRIX_SUBDIAG, 
     &        N_BVTELMATRIX_SUPDIAG, N_LAYER_WFS,
     &        BANDTELMAT2, MAXBANDTOTAL, IPIVOTTEL,
     &        COLTEL2_WF, MAXTOTAL, INFO )

C  (error tracing)

        IF ( INFO .LT. 0 ) THEN
         WRITE(CI, '(I3)' ) INFO
         MAIL  = 'argument i illegal value, for i = '//CI
         TRACE = 'DGBTRS call in L_BVPTEL_SOLUTION_MASTER'
         STATUS = VLIDORT_SERIOUS
         CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
         RETURN
        ENDIF

C  Set linearized integration constants, active layers

        C0 = - NSTKS_NSTRMS_2
        DO NS = 1, NLAYERS_TEL
         N = ACTIVE_LAYERS(NS)
         C0 = C0 + NSTKS_NSTRMS_2

C  set real constants from the solution vector

         DO K = 1, K_REAL(N)
          IROW  = K
          IROW1 = IROW + NSTKS_NSTRMS
          DO Q = 1, N_LAYER_WFS
            NCON(K,N,Q) = COLTEL2_WF(C0+IROW,Q)
            PCON(K,N,Q) = COLTEL2_WF(C0+IROW1,Q)
          ENDDO
         ENDDO

C  set complex constants from the solution vector

         KO1 = K_REAL(N) + 1
         DO K = 1, K_COMPLEX(N)
          K0 = 2*K - 2
          K1 = KO1 + K0
          K2 = K1  + 1
          IROW    = K + K_REAL(N)
          IROW1   = IROW + NSTKS_NSTRMS
          IROW_S  = K + K_REAL(N) + K_COMPLEX(N)
          IROW1_S = IROW_S + NSTKS_NSTRMS
          DO Q = 1, N_LAYER_WFS
            NCON(K1,N,Q) = COLTEL2_WF(C0+IROW,    Q)
            NCON(K2,N,Q) = COLTEL2_WF(C0+IROW_S,  Q)
            PCON(K1,N,Q) = COLTEL2_WF(C0+IROW1,   Q)
            PCON(K2,N,Q) = COLTEL2_WF(C0+IROW1_S, Q)
          ENDDO
         ENDDO

C  End number of telescoped layers

        ENDDO

C  Solve linearized BVP: Single Layer only
C  =======================================

      ELSE IF ( NLAYERS_TEL .EQ. 1 ) THEN

        NAC = ACTIVE_LAYERS(1)

C  LAPACK substitution (DGETRS) using RHS column vector SCOL2_WF

        CALL DGETRS
     &     ( 'N', NSTKS_NSTRMS_2, N_LAYER_WFS, 
     &        SMAT2, MAXSTRMSTKS_2, SIPIVOT,
     &        SCOL2_WF, MAXSTRMSTKS_2, INFO )

C  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          MAIL='BVP back-substitution (telescoping) (DGETRS)'
          CALL VLIDORT_LAPACK_ERROR ( INFO, 1, MAIL, STATUS )
          IF ( STATUS .EQ. VLIDORT_SERIOUS ) RETURN
        ENDIF

C  set real constants from the solution vector

        DO K = 1, K_REAL(NAC)
          IROW  = K
          IROW1 = IROW + NSTKS_NSTRMS
          DO Q = 1, N_LAYER_WFS
            NCON(K,NAC,Q) = SCOL2_WF(IROW,Q)
            PCON(K,NAC,Q) = SCOL2_WF(IROW1,Q)
          ENDDO
        ENDDO

C  set complex constants from the solution vector

        KO1 = K_REAL(NAC) + 1
        DO K = 1, K_COMPLEX(NAC)
          K0 = 2*K - 2
          K1 = KO1 + K0
          K2 = K1  + 1
          IROW    = K + K_REAL(NAC)
          IROW1   = IROW + NSTKS_NSTRMS
          IROW_S  = IROW + K_COMPLEX(NAC)
          IROW1_S = IROW_S + NSTKS_NSTRMS
          DO Q = 1, N_LAYER_WFS
            NCON(K1,NAC,Q) = SCOL2_WF(IROW,    Q)
            NCON(K2,NAC,Q) = SCOL2_WF(IROW_S,  Q)
            PCON(K1,NAC,Q) = SCOL2_WF(IROW1,   Q)
            PCON(K2,NAC,Q) = SCOL2_WF(IROW1_S, Q)
          ENDDO
        ENDDO

      ENDIF

C  Set linearized integration constants for non-active layers
C  ==========================================================

C  Now we propagate the results upwards and downwards through the
C  appropriate non-active layers where there is no scattering.

C  Column linearization case is treated separately
C    New code by R. Spurr, 7 December 2008

      IF ( DO_COLUMN_LINEARIZATION ) GO TO 4554

C  Situation 1. Profile linearization
C  **********************************

C  Transmittance layers ABOVE active layer(s)
C  -----------------------------------------

C   --NCON values are zero (no downwelling radiation)
C   --PCON values propagated upwards from top of first active layer

C  layer immediately above first active layer
C   --- Require linearized solutions at top of first active layer
C   --- Additional linearizations required if the first active
C       layer is the varying layer

      NAC = ACTIVE_LAYERS(1)
      IF ( NAC .GT. 1 ) THEN

        N1 = NAC - 1

C  Case 1. If active layer is also the varying layer
C  -------------------------------------------------

        IF ( LAYER_TO_VARY.EQ.NAC ) THEN

C  start stream, stokes loops

         DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = ( I1 - 1 ) * NSTOKES
          IC = ( I - 1  ) * NSTOKES
          DO O1 = 1, NSTOKES
           IROW = IR + O1
           ICOW = IC + O1

C  start parameter loop

           DO Q = 1, N_LAYER_WFS

C  real homogeneous solutions
C    Needs checking------------------------ 22 March 2007

             SHOM_R = ZERO
             DO K = 1, K_REAL(NAC)
              NXR  = NCON(K,NAC,Q) *   SOLA_XPOS(I1,O1,K,NAC)
              PXR  = PCON(K,NAC,Q) *   SOLB_XNEG(I1,O1,K,NAC)
              MXR  = MCON(K,NAC)   *   SOLB_XNEG(I1,O1,K,NAC)
              LLXR = LCON(K,NAC)   * L_SOLA_XPOS(I1,O1,K,NAC,Q)
              MLXR = MCON(K,NAC)   * L_SOLB_XNEG(I1,O1,K,NAC,Q)             
              L_HOM1 = NXR + LLXR
              L_HOM2 =   T_DELT_EIGEN(K,NAC)   * ( PXR + MLXR ) +
     &                 L_T_DELT_EIGEN(K,NAC,Q) *   MXR
              SHOM_R = SHOM_R + L_HOM1 + L_HOM2
             ENDDO

C  complex homogeneous solutions

             SHOM_CR = ZERO
             KO1 = K_REAL(NAC) + 1
             DO K = 1, K_COMPLEX(NAC)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1

              NXR1  =   NCON(K1,NAC,Q) *   SOLA_XPOS(I1,O1,K1,NAC)
     &                - NCON(K2,NAC,Q) *   SOLA_XPOS(I1,O1,K2,NAC)
              PXR1  =   PCON(K1,NAC,Q) *   SOLB_XNEG(I1,O1,K1,NAC)
     &                - PCON(K2,NAC,Q) *   SOLB_XNEG(I1,O1,K2,NAC)
              PXR2  =   PCON(K1,NAC,Q) *   SOLB_XNEG(I1,O1,K2,NAC)
     &                + PCON(K2,NAC,Q) *   SOLB_XNEG(I1,O1,K1,NAC)

              MXR1  =   MCON(K1,NAC) *   SOLB_XNEG(I1,O1,K1,NAC)
     &                - MCON(K2,NAC) *   SOLB_XNEG(I1,O1,K2,NAC)
              MXR2  =   MCON(K1,NAC) *   SOLB_XNEG(I1,O1,K2,NAC)
     &                + MCON(K2,NAC) *   SOLB_XNEG(I1,O1,K1,NAC)

              LLXR1  =   LCON(K1,NAC) * L_SOLA_XPOS(I1,O1,K1,NAC,Q)
     &                 - LCON(K2,NAC) * L_SOLA_XPOS(I1,O1,K2,NAC,Q)
              MLXR1  =   MCON(K1,NAC) * L_SOLB_XNEG(I1,O1,K1,NAC,Q)
     &                 - MCON(K2,NAC) * L_SOLB_XNEG(I1,O1,K2,NAC,Q)
              MLXR2  =   MCON(K1,NAC) * L_SOLB_XNEG(I1,O1,K2,NAC,Q)
     &                 + MCON(K2,NAC) * L_SOLB_XNEG(I1,O1,K1,NAC,Q)

              L_HOM1CR = NXR1 + LLXR1
              L_HOM2CR = + T_DELT_EIGEN(K1,NAC)     * ( PXR1 + MLXR1 ) 
     &                   - T_DELT_EIGEN(K2,NAC)     * ( PXR2 + MLXR2 )
     &                   + L_T_DELT_EIGEN(K1,NAC,Q) *   MXR1
     &                   - L_T_DELT_EIGEN(K2,NAC,Q) *   MXR2

              SHOM_CR = SHOM_CR + L_HOM1CR + L_HOM2CR

             ENDDO

C  real part and add particular solution
C    ---Sets Real integration constants (no complex ones)

             SHOM = SHOM_R + SHOM_CR
             SPAR = L_WUPPER(I1,O1,NAC,Q)
             PCON(ICOW,N1,Q) = SPAR + SHOM
             NCON(ICOW,N1,Q) = ZERO
 
C  End loops

           ENDDO
          ENDDO
         ENDDO

C  End case 1

        ENDIF

C  Case 2. Active layer is below varying layer
C  -------------------------------------------

        IF ( LAYER_TO_VARY.LT.NAC) THEN

C  start stream, stokes loops

         DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = ( I1 - 1 ) * NSTOKES
          IC = ( I - 1  ) * NSTOKES
          DO O1 = 1, NSTOKES
           IROW = IR + O1
           ICOW = IC + O1

C  start parameter loop

           DO Q = 1, N_LAYER_WFS

C  real homogeneous solutions
C      Needs Checking ---------------------- 22 March 2007

             SHOM_R = ZERO
             DO K = 1, K_REAL(NAC)
              NXR = NCON(K,NAC,Q)*SOLA_XPOS(I1,O1,K,NAC)
              PXR = PCON(K,NAC,Q)*SOLB_XNEG(I1,O1,K,NAC)
              L_HOM1 = NXR
              L_HOM2 = PXR * T_DELT_EIGEN(K,NAC)
              SHOM_R = SHOM_R + L_HOM1 + L_HOM2
             ENDDO

C  complex homogeneous solutions

             SHOM_CR = ZERO
             KO1 = K_REAL(NAC) + 1
             DO K = 1, K_COMPLEX(NAC)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1

              NXR1  =   NCON(K1,NAC,Q) *   SOLA_XPOS(I1,O1,K1,NAC)
     &                - NCON(K2,NAC,Q) *   SOLA_XPOS(I1,O1,K2,NAC)
              PXR1  =   PCON(K1,NAC,Q) *   SOLB_XNEG(I1,O1,K1,NAC)
     &                - PCON(K2,NAC,Q) *   SOLB_XNEG(I1,O1,K2,NAC)
              PXR2  =   PCON(K1,NAC,Q) *   SOLB_XNEG(I1,O1,K2,NAC)
     &                + PCON(K2,NAC,Q) *   SOLB_XNEG(I1,O1,K1,NAC)

              L_HOM1CR = NXR1
              L_HOM2CR =  + T_DELT_EIGEN(K1,NAC) * PXR1
     &                    - T_DELT_EIGEN(K2,NAC) * PXR2
              SHOM_CR = SHOM_CR + L_HOM1CR + L_HOM2CR
             ENDDO

C  real part and add particular solution
C    ---Sets Real integration constants (no complex ones)

             SHOM = SHOM_R + SHOM_CR
             SPAR = L_WUPPER(I1,O1,NAC,Q)
             PCON(ICOW,N1,Q) = SPAR + SHOM
             NCON(ICOW,N1,Q) = ZERO
 
C  End loops

           ENDDO
          ENDDO
         ENDDO

C  End case 2

        ENDIF

C  Case 3. Active layer is above varying layer
C  -------------------------------------------

        IF ( LAYER_TO_VARY.GT.NAC) THEN

C  No solutions

         DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = ( I1 - 1 ) * NSTOKES
          IC = ( I - 1  ) * NSTOKES
          DO O1 = 1, NSTOKES
           IROW = IR + O1
           ICOW = IC + O1
           DO Q = 1, N_LAYER_WFS
             PCON(ICOW,N1,Q) = ZERO
             NCON(ICOW,N1,Q) = ZERO
           ENDDO
          ENDDO
         ENDDO

C  End case 3

        ENDIF

      ENDIF

C  For remaining non-active atmospheric layers to TOA, propagate upwards.
C   Additional linearizations if you are passing through the varying layer.

      DO N = NAC - 2, 1, -1
        N1 = N + 1
        DO I = 1, NSTREAMS
          IR = ( I - 1 ) * NSTOKES
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            DO Q = 1, N_LAYER_WFS
              NCON(IROW,N,Q) = ZERO
              PCON(IROW,N,Q) =
     &            T_DELT_DISORDS(I,N1) * PCON(IROW,N1,Q)
            ENDDO
          ENDDO
        ENDDO
        IF ( N1 .EQ. LAYER_TO_VARY ) THEN
         DO I = 1, NSTREAMS
          IR = ( I - 1 ) * NSTOKES
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            DO Q = 1, N_LAYER_WFS
              PCON(IROW,N,Q) = PCON(IROW,N,Q)
     &            + L_T_DELT_DISORDS(I,N1,Q) * MCON(IROW,N1)
            ENDDO
          ENDDO
         ENDDO
        ENDIF
      ENDDO

C  Transmittance layers below active layer(s)
C  -----------------------------------------

C   -- PCON values are zero (no upwelling radiation)
C   -- NCON values propagated downwards from bottom of last active layer

C  layer immediately below Last active layer
C    .... Require linearized solutions at bottom of last active layer
C    .... Additional linearizations is the last active layer is also
C         the varying layer.

      NAC = ACTIVE_LAYERS (NLAYERS_TEL)
      IF ( NAC .LT. NLAYERS ) THEN
        N1 = NAC + 1

C  Case 1. If active layer equals the varying layer
C  ------------------------------------------------

        IF ( LAYER_TO_VARY .EQ. NAC ) THEN

C  start stream, stokes loops

         DO I = 1, NSTREAMS
          IR = ( I - 1 ) * NSTOKES
          DO O1 = 1, NSTOKES
           IROW = IR + O1

C  start parameter loop

           DO Q = 1, N_LAYER_WFS

C  real homogeneous solutions
C   Needs checking------------------------- 22 March 2007

             SHOM_R = ZERO
             DO K = 1, K_REAL(NAC)
              NXR  = NCON(K,NAC,Q) *   SOLA_XPOS(I,O1,K,NAC)
              PXR  = PCON(K,NAC,Q) *   SOLB_XNEG(I,O1,K,NAC)
              LXR  = LCON(K,NAC)   *   SOLA_XPOS(I,O1,K,NAC)
              LLXR = LCON(K,NAC)   * L_SOLA_XPOS(I,O1,K,NAC,Q)
              MLXR = MCON(K,NAC)   * L_SOLB_XNEG(I,O1,K,NAC,Q)             
              L_HOM2 = PXR + MLXR
              L_HOM1 =   T_DELT_EIGEN(K,NAC)   * ( NXR + LLXR ) +
     &                 L_T_DELT_EIGEN(K,NAC,Q) *   LXR
              SHOM_R = SHOM_R + L_HOM1 + L_HOM2
             ENDDO

C  complex homogeneous solutions

             SHOM_CR = ZERO
             KO1 = K_REAL(NAC) + 1
             DO K = 1, K_COMPLEX(NAC)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1

              NXR1  =   NCON(K1,NAC,Q) *   SOLA_XPOS(I,O1,K1,NAC)
     &                - NCON(K2,NAC,Q) *   SOLA_XPOS(I,O1,K2,NAC)
              NXR2  =   NCON(K1,NAC,Q) *   SOLA_XPOS(I,O1,K2,NAC)
     &                + NCON(K2,NAC,Q) *   SOLA_XPOS(I,O1,K1,NAC)
              PXR1  =   PCON(K1,NAC,Q) *   SOLB_XNEG(I,O1,K1,NAC)
     &                - PCON(K2,NAC,Q) *   SOLB_XNEG(I,O1,K2,NAC)

              LXR1  =   LCON(K1,NAC) *   SOLA_XPOS(I,O1,K1,NAC)
     &                - LCON(K2,NAC) *   SOLA_XPOS(I,O1,K2,NAC)
              LXR2  =   LCON(K1,NAC) *   SOLA_XPOS(I,O1,K2,NAC)
     &                + LCON(K2,NAC) *   SOLA_XPOS(I,O1,K1,NAC)

              LLXR1  =   LCON(K1,NAC) * L_SOLA_XPOS(I,O1,K1,NAC,Q)
     &                 - LCON(K2,NAC) * L_SOLA_XPOS(I,O1,K2,NAC,Q)
              LLXR2  =   LCON(K1,NAC) * L_SOLA_XPOS(I,O1,K2,NAC,Q)
     &                 + LCON(K2,NAC) * L_SOLA_XPOS(I,O1,K1,NAC,Q)
              MLXR1  =   MCON(K1,NAC) * L_SOLB_XNEG(I,O1,K1,NAC,Q)
     &                 - MCON(K2,NAC) * L_SOLB_XNEG(I,O1,K2,NAC,Q)

              L_HOM2CR = PXR1 + MLXR1
              L_HOM1CR = + T_DELT_EIGEN(K1,NAC)     * ( NXR1 + LLXR1 ) 
     &                   - T_DELT_EIGEN(K2,NAC)     * ( NXR2 + LLXR2 )
     &                   + L_T_DELT_EIGEN(K1,NAC,Q) *   LXR1
     &                   - L_T_DELT_EIGEN(K2,NAC,Q) *   LXR2

              SHOM_CR = SHOM_CR + L_HOM1CR + L_HOM2CR

             ENDDO

C  real part and add particular solution
C    ---Sets Real integration constants (no complex ones)

             SHOM = SHOM_R + SHOM_CR
             SPAR = L_WLOWER(I,O1,NAC,Q)
c             PCON(IROW,N1,Q) = SPAR + SHOM
c             NCON(IROW,N1,Q) = ZERO
             NCON(IROW,N1,Q) = SPAR + SHOM
             PCON(IROW,N1,Q) = ZERO
 
C  End loops

           ENDDO
          ENDDO
         ENDDO

C  End case 1

        ENDIF

C  Case 2. Active layer is below varying layer
C  -------------------------------------------

        IF ( LAYER_TO_VARY.LT.NAC) THEN

C  start stream, stokes loops

         DO I = 1, NSTREAMS
          IR = ( I - 1 ) * NSTOKES
          DO O1 = 1, NSTOKES
           IROW = IR + O1

C  start parameter loop

           DO Q = 1, N_LAYER_WFS

C  real homogeneous solutions
C   Needs checking------------------------- 22 March 2007

             SHOM_R = ZERO
             DO K = 1, K_REAL(NAC)
              NXR = NCON(K,NAC,Q) * SOLA_XPOS(I,O1,K,NAC)
              PXR = PCON(K,NAC,Q) * SOLB_XNEG(I,O1,K,NAC)
              L_HOM2 = PXR
              L_HOM1 = NXR * T_DELT_EIGEN(K,NAC)
              SHOM_R = SHOM_R + L_HOM1 + L_HOM2
             ENDDO

C  complex homogeneous solutions

             SHOM_CR = ZERO
             KO1 = K_REAL(NAC) + 1
             DO K = 1, K_COMPLEX(NAC)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1

              NXR1  =   NCON(K1,NAC,Q) *   SOLA_XPOS(I,O1,K1,NAC)
     &                - NCON(K2,NAC,Q) *   SOLA_XPOS(I,O1,K2,NAC)
              NXR2  =   NCON(K1,NAC,Q) *   SOLA_XPOS(I,O1,K2,NAC)
     &                + NCON(K2,NAC,Q) *   SOLA_XPOS(I,O1,K1,NAC)
              PXR1  =   PCON(K1,NAC,Q) *   SOLB_XNEG(I,O1,K1,NAC)
     &                - PCON(K2,NAC,Q) *   SOLB_XNEG(I,O1,K2,NAC)

              L_HOM2CR = PXR1
              L_HOM1CR =  + T_DELT_EIGEN(K1,NAC) * NXR1
     &                    - T_DELT_EIGEN(K2,NAC) * NXR2
              SHOM_CR = SHOM_CR + L_HOM1CR + L_HOM2CR
             ENDDO

C  real part and add particular solution
C    ---Sets Real integration constants (no complex ones)

             SHOM = SHOM_R + SHOM_CR
             SPAR = L_WLOWER(I,O1,NAC,Q)
             NCON(IROW,N1,Q) = SPAR + SHOM
             PCON(IROW,N1,Q) = ZERO

C  End loops

           ENDDO
          ENDDO
         ENDDO

C  End case 2

        ENDIF

C  Case 3. Active layer is above varying layer
C  -------------------------------------------

        IF ( LAYER_TO_VARY.GT.NAC) THEN

C  No solutions

         DO I = 1, NSTREAMS
          IR = ( I - 1 ) * NSTOKES
          DO O1 = 1, NSTOKES
           IROW = IR + O1
           DO Q = 1, N_LAYER_WFS
             PCON(IROW,N1,Q) = ZERO
             NCON(IROW,N1,Q) = ZERO
           ENDDO
          ENDDO
         ENDDO

C  End case 3

        ENDIF

      ENDIF

C  other layers to bottom of medium: propagate downwards.
C   Additional variation if you are passing through the varying layer.

      DO N = NAC + 2, NLAYERS
        N1 = N - 1
        DO I = 1, NSTREAMS
          IR = ( I - 1 ) * NSTOKES
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            DO Q = 1, N_LAYER_WFS
              PCON(IROW,N,Q) = ZERO
              NCON(IROW,N,Q) =
     &            T_DELT_DISORDS(I,N1) * NCON(IROW,N1,Q)
            ENDDO
          ENDDO
        ENDDO
        IF ( N1 .EQ. LAYER_TO_VARY ) THEN
         DO I = 1, NSTREAMS
          IR = ( I - 1 ) * NSTOKES
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            DO Q = 1, N_LAYER_WFS
              NCON(IROW,N,Q) = NCON(IROW,N,Q)
     &            + L_T_DELT_DISORDS(I,N1,Q) * LCON(IROW,N1)
            ENDDO
          ENDDO
         ENDDO
        ENDIF
      ENDDO

C  debug

c      if ( layer_to_vary .eq. 1 ) then
c       if ( fourier_component .EQ.3 .and.ibeam.eq.1 ) then
c        do n = 1, nlayers
c          write(99,'(a,i3,1p2e18.10)')'hey',n,ncon(3,n,1),pcon(3,n,1)
c        enddo
c       endif
c      endif

C  debug

c      q = 96
c      if ( fourier_component.eq.0.and.ibeam.eq.1
c     &          .and.layer_to_vary.eq.12) then
c        DO N = 1, NLAYERS
c          DO I = 1, NSTREAMS*2
c          write(q,'(3i4,1p4e17.9)')IBEAM,N,I, NCON(I,N,1), PCON(I,N,1)
c         ENDDO
c        ENDDO
c      ENDIF

C  Go to final stuff

      GO TO 4555

C  Situation 2. Column linearization
C  *********************************

C Continuation point for doing this

 4554 continue

C   --NCON values are zero (no downwelling radiation)
C   --PCON values propagated upwards from top of first active layer

C  layer immediately above first active layer
C   --- Require linearized solutions at top of first active layer
C   --- Additional linearizations required, first layer is always active

      NAC = ACTIVE_LAYERS(1)
      IF ( NAC .GT. 1 ) THEN

       N1 = NAC - 1

C  start stream, stokes loops

       DO I = 1, NSTREAMS
         I1 = I + NSTREAMS
         IR = ( I1 - 1 ) * NSTOKES
         IC = ( I - 1  ) * NSTOKES
         DO O1 = 1, NSTOKES
           IROW = IR + O1
           ICOW = IC + O1

C  start parameter loop

           DO Q = 1, N_LAYER_WFS

C  real homogeneous solutions
C    Needs checking------------------------ 22 March 2007

             SHOM_R = ZERO
             DO K = 1, K_REAL(NAC)
              NXR  = NCON(K,NAC,Q) *   SOLA_XPOS(I1,O1,K,NAC)
              PXR  = PCON(K,NAC,Q) *   SOLB_XNEG(I1,O1,K,NAC)
              MXR  = MCON(K,NAC)   *   SOLB_XNEG(I1,O1,K,NAC)
              LLXR = LCON(K,NAC)   * L_SOLA_XPOS(I1,O1,K,NAC,Q)
              MLXR = MCON(K,NAC)   * L_SOLB_XNEG(I1,O1,K,NAC,Q)             
              L_HOM1 = NXR + LLXR
              L_HOM2 =   T_DELT_EIGEN(K,NAC)   * ( PXR + MLXR ) +
     &                 L_T_DELT_EIGEN(K,NAC,Q) *   MXR
              SHOM_R = SHOM_R + L_HOM1 + L_HOM2
             ENDDO

C  complex homogeneous solutions

             SHOM_CR = ZERO
             KO1 = K_REAL(NAC) + 1
             DO K = 1, K_COMPLEX(NAC)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1

              NXR1  =   NCON(K1,NAC,Q) *   SOLA_XPOS(I1,O1,K1,NAC)
     &                - NCON(K2,NAC,Q) *   SOLA_XPOS(I1,O1,K2,NAC)
              PXR1  =   PCON(K1,NAC,Q) *   SOLB_XNEG(I1,O1,K1,NAC)
     &                - PCON(K2,NAC,Q) *   SOLB_XNEG(I1,O1,K2,NAC)
              PXR2  =   PCON(K1,NAC,Q) *   SOLB_XNEG(I1,O1,K2,NAC)
     &                + PCON(K2,NAC,Q) *   SOLB_XNEG(I1,O1,K1,NAC)

              MXR1  =   MCON(K1,NAC) *   SOLB_XNEG(I1,O1,K1,NAC)
     &                - MCON(K2,NAC) *   SOLB_XNEG(I1,O1,K2,NAC)
              MXR2  =   MCON(K1,NAC) *   SOLB_XNEG(I1,O1,K2,NAC)
     &                + MCON(K2,NAC) *   SOLB_XNEG(I1,O1,K1,NAC)

              LLXR1  =   LCON(K1,NAC) * L_SOLA_XPOS(I1,O1,K1,NAC,Q)
     &                 - LCON(K2,NAC) * L_SOLA_XPOS(I1,O1,K2,NAC,Q)
              MLXR1  =   MCON(K1,NAC) * L_SOLB_XNEG(I1,O1,K1,NAC,Q)
     &                 - MCON(K2,NAC) * L_SOLB_XNEG(I1,O1,K2,NAC,Q)
              MLXR2  =   MCON(K1,NAC) * L_SOLB_XNEG(I1,O1,K2,NAC,Q)
     &                 + MCON(K2,NAC) * L_SOLB_XNEG(I1,O1,K1,NAC,Q)

              L_HOM1CR = NXR1 + LLXR1
              L_HOM2CR = + T_DELT_EIGEN(K1,NAC)     * ( PXR1 + MLXR1 ) 
     &                   - T_DELT_EIGEN(K2,NAC)     * ( PXR2 + MLXR2 )
     &                   + L_T_DELT_EIGEN(K1,NAC,Q) *   MXR1
     &                   - L_T_DELT_EIGEN(K2,NAC,Q) *   MXR2

              SHOM_CR = SHOM_CR + L_HOM1CR + L_HOM2CR

             ENDDO

C  real part and add particular solution
C    ---Sets Real integration constants (no complex ones)

             SHOM = SHOM_R + SHOM_CR
             SPAR = L_WUPPER(I1,O1,NAC,Q)
             PCON(ICOW,N1,Q) = SPAR + SHOM
             NCON(ICOW,N1,Q) = ZERO
 
C  End loops

           ENDDO
         ENDDO
       ENDDO

C  End active layer varying

      ENDIF

C  For remaining non-active atmospheric layers to TOA, propagate upwards.
C   Additional linearizations if you are passing through the varying layer.

      DO N = NAC - 2, 1, -1
        N1 = N + 1
        DO I = 1, NSTREAMS
          IR = ( I - 1 ) * NSTOKES
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            DO Q = 1, N_LAYER_WFS
              NCON(IROW,N,Q) = ZERO
              PCON(IROW,N,Q) =
     &                T_DELT_DISORDS(I,N1)   * PCON(IROW,N1,Q)
     &            + L_T_DELT_DISORDS(I,N1,Q) * MCON(IROW,N1)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

C  Transmittance layers below active layer(s)
C  -----------------------------------------

C   -- PCON values are zero (no upwelling radiation)
C   -- NCON values propagated downwards from bottom of last active layer

C  layer immediately below Last active layer
C    .... Require linearized solutions at bottom of last active layer
C    .... Additional linearizations in the last active layer, always active

      NAC = ACTIVE_LAYERS (NLAYERS_TEL)
      IF ( NAC .LT. NLAYERS ) THEN
       N1 = NAC + 1

C  start stream, stokes loops

       DO I = 1, NSTREAMS
         IR = ( I - 1 ) * NSTOKES
         DO O1 = 1, NSTOKES
           IROW = IR + O1

C  start parameter loop

           DO Q = 1, N_LAYER_WFS

C  real homogeneous solutions
C   Needs checking------------------------- 22 March 2007

             SHOM_R = ZERO
             DO K = 1, K_REAL(NAC)
              NXR  = NCON(K,NAC,Q) *   SOLA_XPOS(I,O1,K,NAC)
              PXR  = PCON(K,NAC,Q) *   SOLB_XNEG(I,O1,K,NAC)
              LXR  = LCON(K,NAC)   *   SOLA_XPOS(I,O1,K,NAC)
              LLXR = LCON(K,NAC)   * L_SOLA_XPOS(I,O1,K,NAC,Q)
              MLXR = MCON(K,NAC)   * L_SOLB_XNEG(I,O1,K,NAC,Q)             
              L_HOM2 = PXR + MLXR
              L_HOM1 =   T_DELT_EIGEN(K,NAC)   * ( NXR + LLXR ) +
     &                 L_T_DELT_EIGEN(K,NAC,Q) *   LXR
              SHOM_R = SHOM_R + L_HOM1 + L_HOM2
             ENDDO

C  complex homogeneous solutions

             SHOM_CR = ZERO
             KO1 = K_REAL(NAC) + 1
             DO K = 1, K_COMPLEX(NAC)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1

              NXR1  =   NCON(K1,NAC,Q) *   SOLA_XPOS(I,O1,K1,NAC)
     &                - NCON(K2,NAC,Q) *   SOLA_XPOS(I,O1,K2,NAC)
              NXR2  =   NCON(K1,NAC,Q) *   SOLA_XPOS(I,O1,K2,NAC)
     &                + NCON(K2,NAC,Q) *   SOLA_XPOS(I,O1,K1,NAC)
              PXR1  =   PCON(K1,NAC,Q) *   SOLB_XNEG(I,O1,K1,NAC)
     &                - PCON(K2,NAC,Q) *   SOLB_XNEG(I,O1,K2,NAC)

              LXR1  =   LCON(K1,NAC) *   SOLA_XPOS(I,O1,K1,NAC)
     &                - LCON(K2,NAC) *   SOLA_XPOS(I,O1,K2,NAC)
              LXR2  =   LCON(K1,NAC) *   SOLA_XPOS(I,O1,K2,NAC)
     &                + LCON(K2,NAC) *   SOLA_XPOS(I,O1,K1,NAC)

              LLXR1  =   LCON(K1,NAC) * L_SOLA_XPOS(I,O1,K1,NAC,Q)
     &                 - LCON(K2,NAC) * L_SOLA_XPOS(I,O1,K2,NAC,Q)
              LLXR2  =   LCON(K1,NAC) * L_SOLA_XPOS(I,O1,K2,NAC,Q)
     &                 + LCON(K2,NAC) * L_SOLA_XPOS(I,O1,K1,NAC,Q)
              MLXR1  =   MCON(K1,NAC) * L_SOLB_XNEG(I,O1,K1,NAC,Q)
     &                 - MCON(K2,NAC) * L_SOLB_XNEG(I,O1,K2,NAC,Q)

              L_HOM2CR = PXR1 + MLXR1
              L_HOM1CR = + T_DELT_EIGEN(K1,NAC)     * ( NXR1 + LLXR1 ) 
     &                   - T_DELT_EIGEN(K2,NAC)     * ( NXR2 + LLXR2 )
     &                   + L_T_DELT_EIGEN(K1,NAC,Q) *   LXR1
     &                   - L_T_DELT_EIGEN(K2,NAC,Q) *   LXR2

              SHOM_CR = SHOM_CR + L_HOM1CR + L_HOM2CR

             ENDDO

C  real part and add particular solution
C    ---Sets Real integration constants (no complex ones)

             SHOM = SHOM_R + SHOM_CR
             SPAR = L_WLOWER(I,O1,NAC,Q)
c             PCON(IROW,N1,Q) = SPAR + SHOM
c             NCON(IROW,N1,Q) = ZERO
             NCON(IROW,N1,Q) = SPAR + SHOM
             PCON(IROW,N1,Q) = ZERO
 
C  End loops

           ENDDO
         ENDDO
       ENDDO

C  End active layer

      ENDIF

C  other layers to bottom of medium: propagate downwards.
C   Additional variation if you are passing through the varying layer.

C  other layers to bottom of medium: propagate downwards.
C   Additional variation if you are passing through the varying layer.

      DO N = NAC + 2, NLAYERS
        N1 = N - 1
        DO I = 1, NSTREAMS
          IR = ( I - 1 ) * NSTOKES
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            DO Q = 1, N_LAYER_WFS
              PCON(IROW,N,Q) = ZERO
              NCON(IROW,N,Q) =
     &                T_DELT_DISORDS(I,N1)   * NCON(IROW,N1,Q)
     &            + L_T_DELT_DISORDS(I,N1,Q) * LCON(IROW,N1)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

C  Continuation point for final stuff

 4555 continue

C  finish

      RETURN
      END

C

      SUBROUTINE L_BVPTEL_P_COLUMN_SETUP
     I     ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM,
     I       LAYER_TO_VARY,      N_LAYER_WFS,
     I       FOURIER_COMPONENT,  IBEAM,
     I       SURFACE_FACTOR,
     I       L_BVP_SURFACE_SETUP )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup, solution variables

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'

C  include files of linearized setup variables (input)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'

C  include files of linearized solution variables (output stored here)

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'

C  Include file of reflectance variables

      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'

C  input arguments
C  ---------------

C  inclusion flags

      LOGICAL          DO_INCLUDE_SURFACE
      LOGICAL          DO_INCLUDE_DIRECTBEAM

C  Fourier component and beam number

      INTEGER          FOURIER_COMPONENT, IBEAM

C  surface factor = 1+delta(m,0)

      DOUBLE PRECISION SURFACE_FACTOR

C  Linearization control

      INTEGER          LAYER_TO_VARY
      INTEGER          N_LAYER_WFS

C  External function

      EXTERNAL         L_BVP_SURFACE_SETUP

C  output status
C  -------------

      INTEGER          STATUS

C  local variables
C  ---------------

      INTEGER          Q, N, I, I1, IR, CM, C0, IROW, O1, N1, NS
      INTEGER          K, KO1, K0, K1, K2
      DOUBLE PRECISION CPOS,CNEG,L_HOM_R,L_HOM_CR,L_BEAM
      DOUBLE PRECISION T1,T2,T1R,T1I,T2R,T2I,FAC

      DOUBLE PRECISION R2_L_BEAM
     &        ( MAXSTREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION R2_L_HOMP
     &        ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )
      DOUBLE PRECISION R2_L_HOMM
     &        ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )

C  Try this safety-first zeroing

      DO I = 1, NSTREAMS_2
        DO Q = 1, N_LAYER_WFS
          DO O1 = 1, NSTOKES
            L_WUPPER(I,O1,LAYER_TO_VARY,Q) = ZERO
            L_WLOWER(I,O1,LAYER_TO_VARY,Q) = ZERO
          ENDDO
        ENDDO
      ENDDO

C  status

      status = 0

C  Get the linearized solutions for the layer that is varying
C    Always need this, regardless of number of active layers

      DO NS = 1, NLAYERS_TEL
        N = ACTIVE_LAYERS(NS)
        IF ( N.EQ.LAYER_TO_VARY ) THEN
          CALL L_BEAMSOLUTION_P_NEQK
     I      ( FOURIER_COMPONENT, IBEAM, LAYER_TO_VARY, N_LAYER_WFS )
        ENDIF
      ENDDO

C  Go to special case for only 1 active layer

      IF ( NLAYERS_TEL .EQ. 1 ) GO TO 3456

C  zero column vector

      DO I = 1, N_BVTELMATRIX_SIZE
        DO Q = 1, MAX_ATMOSWFS
          COLTEL2_WF(I,Q) = ZERO
        ENDDO
      ENDDO

C  top of first active layer, first boundary condition
C  ---------------------------------------------------

      NS = 1
      N = ACTIVE_LAYERS(NS)
      C0 = 0

C  If this active layer = layer that is varying,
C       then require homogeneous and beam solution linearizations

      IF ( LAYER_TO_VARY .EQ. N ) THEN

C  start loops

        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            DO Q = 1, N_LAYER_WFS

C  Beam contribution

            L_BEAM  = - L_WUPPER(I,O1,N,Q)

C  Linearized Real homogeneous solution contributions

            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N)
              CPOS = L_SOLA_XPOS(I,O1,K,N,Q)
              CNEG = T_DELT_EIGEN(K,N)   * L_SOLB_XNEG(I,O1,K,N,Q) + 
     &             L_T_DELT_EIGEN(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

C  Linearized Complex homogeneous solution contributions

            L_HOM_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T1 = L_SOLA_XPOS(I,O1,K1,N,Q) * LCON(K1,N) -
     &             L_SOLA_XPOS(I,O1,K2,N,Q) * LCON(K2,N)
              T2R =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K1,N,Q)
     &             - T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K2,N,Q)
     &             + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N)
     &             - L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)
              T2I =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K2,N,Q)
     &             + T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K1,N,Q)
     &             + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K2,N)
     &             + L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K1,N)
              T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

C  Final contribution

            COLTEL2_WF(CM,Q) = L_BEAM - L_HOM_R - L_HOM_CR

C  end loops

          ENDDO
         ENDDO
        ENDDO

C  otherwise if varying layer is above first active layer, there are beam
C  solution contributions propagated downwards - find these by calling
C  the appropriate solution module = L_BEAMSOLUTION_NNEK

      ELSE IF ( LAYER_TO_VARY .LT. N ) THEN

        CALL L_BEAMSOLUTION_P_NNEK
     &    ( FOURIER_COMPONENT, IBEAM, N, LAYER_TO_VARY, N_LAYER_WFS )

        DO I = 1, NSTREAMS
         IR = NSTOKES*(I-1)
         DO O1 = 1, NSTOKES
          IROW = IR + O1
          CM = C0 + IROW
          DO Q = 1, N_LAYER_WFS
            COLTEL2_WF(CM,Q) = - L_WUPPER(I,O1,N,Q)
          ENDDO
         ENDDO
        ENDDO

      ENDIF

C  Intermediate boundaries between active layers
C  ---------------------------------------------

      DO NS = 1, NLAYERS_TEL - 1

C  offsets

       N  = ACTIVE_LAYERS(NS)
       N1 = N + 1
       C0 = NS*NSTKS_NSTRMS_2 - NSTKS_NSTRMS

C  if N is the varying layer, immediately above boundary

       IF ( N .EQ. LAYER_TO_VARY ) THEN

C  Get the linearized beam solution for the next layer N1

        CALL L_BEAMSOLUTION_P_NNEK
     I  ( FOURIER_COMPONENT, IBEAM, N1, LAYER_TO_VARY, N_LAYER_WFS )

C  start loops

        DO I = 1, NSTREAMS_2
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            DO Q = 1, N_LAYER_WFS

C  Beam contributions

            L_BEAM = L_WUPPER(I,O1,N1,Q) - L_WLOWER(I,O1,N,Q)

C  Linearized Real homogeneous solution contributions

            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N)
              CNEG = L_SOLB_XNEG(I,O1,K,N,Q)
              CPOS = T_DELT_EIGEN(K,N)   * L_SOLA_XPOS(I,O1,K,N,Q) + 
     &             L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I,O1,K,N)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

C  Linearized Complex homogeneous solution contributions

            L_HOM_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T2 = L_SOLB_XNEG(I,O1,K1,N,Q) * MCON(K1,N) -
     &             L_SOLB_XNEG(I,O1,K2,N,Q) * MCON(K2,N)
              T1R =  T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I,O1,K1,N,Q)
     &             - T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I,O1,K2,N,Q)
     &             + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I,O1,K1,N)
     &             - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I,O1,K2,N)
              T1I =  T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I,O1,K2,N,Q)
     &             + T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I,O1,K1,N,Q)
     &             + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I,O1,K2,N)
     &             + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I,O1,K1,N)
              T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

C  Final contribution

            COLTEL2_WF(CM,Q) = L_BEAM - L_HOM_R - L_HOM_CR

C  End loops

          ENDDO
         ENDDO
        ENDDO

C  If N1 is the varying layer, immediately below boundary
C    Only require contributions from this layer

       ELSE IF ( N1 .EQ. LAYER_TO_VARY ) THEN

C  start loops

        DO I = 1, NSTREAMS_2
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            DO Q = 1, N_LAYER_WFS

C  Beam contribution

            L_BEAM  = + L_WUPPER(I,O1,N1,Q)

C  Linearized Real homogeneous solution contributions

            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N1)
              CPOS = L_SOLA_XPOS(I,O1,K,N1,Q)
              CNEG = T_DELT_EIGEN(K,N1)   * L_SOLB_XNEG(I,O1,K,N1,Q) + 
     &             L_T_DELT_EIGEN(K,N1,Q) *   SOLB_XNEG(I,O1,K,N1)
              T1 = LCON(K,N1) * CPOS
              T2 = MCON(K,N1) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

C  Linearized Complex homogeneous solution contributions

            L_HOM_CR  = ZERO
            KO1 = K_REAL(N1) + 1
            DO K = 1, K_COMPLEX(N1)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T1 = L_SOLA_XPOS(I,O1,K1,N1,Q) * LCON(K1,N1) -
     &             L_SOLA_XPOS(I,O1,K2,N1,Q) * LCON(K2,N1)
              T2R =  T_DELT_EIGEN(K1,N1)   * L_SOLB_XNEG(I,O1,K1,N1,Q)
     &             - T_DELT_EIGEN(K2,N1)   * L_SOLB_XNEG(I,O1,K2,N1,Q)
     &             + L_T_DELT_EIGEN(K1,N1,Q) *   SOLB_XNEG(I,O1,K1,N1)
     &             - L_T_DELT_EIGEN(K2,N1,Q) *   SOLB_XNEG(I,O1,K2,N1)
              T2I =  T_DELT_EIGEN(K1,N1)   * L_SOLB_XNEG(I,O1,K2,N1,Q)
     &             + T_DELT_EIGEN(K2,N1)   * L_SOLB_XNEG(I,O1,K1,N1,Q)
     &             + L_T_DELT_EIGEN(K1,N1,Q) *   SOLB_XNEG(I,O1,K2,N1)
     &             + L_T_DELT_EIGEN(K2,N1,Q) *   SOLB_XNEG(I,O1,K1,N1)
              T2 =  T2R * MCON(K1,N1) - T2I * MCON(K2,N1)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

C  Final contribution

            COLTEL2_WF(CM,Q) = L_BEAM + L_HOM_R + L_HOM_CR

C  end loops

          ENDDO
         ENDDO
        ENDDO

C  non-zero variations if LAYER_TO_VARY is an active layer above N
C    Get the linearized beam solution for the next layer
C  .. contributions from beam solutions on both sides.

       ELSE IF ( LAYER_TO_VARY .LT. N ) THEN

C  Get the linearized beam solution for the next layer

        CALL L_BEAMSOLUTION_P_NNEK
     I    ( FOURIER_COMPONENT, IBEAM, N1, LAYER_TO_VARY, N_LAYER_WFS )

C  .. contributions from beam solution (direct assign). No homog. variation

        DO I = 1, NSTREAMS_2
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            DO Q = 1, N_LAYER_WFS
              L_BEAM = L_WUPPER(I,O1,N1,Q) - L_WLOWER(I,O1,N,Q)
              COLTEL2_WF(CM,Q) = L_BEAM
            ENDDO
          ENDDO
        ENDDO

C  Finish different types of boundary condition linearizations

       ENDIF

C  End loop over intermediate active layer boundaries

      ENDDO

C  Final boundary, bottom of lowest active layer
C  ---------------------------------------------

      NS = NLAYERS_TEL
      N  = ACTIVE_LAYERS(NS)      
      C0 = (NS-1)*NSTKS_NSTRMS_2 + NSTKS_NSTRMS

C  If this is the surface and Specialist option #2 is in place

      if ( DO_INCLUDE_SURFACE.AND.DO_SPECIALIST_OPTION_2
     &       .AND.  N.EQ.NLAYERS ) THEN

C  Modified BCL4M Component loop

        IF ( LAYER_TO_VARY .EQ. N ) THEN

C  get the linearized downward-reflected term

         CALL L_BVP_SURFACE_SETUP
     I     ( DO_INCLUDE_SURFACE, .TRUE., IBEAM,
     I       FOURIER_COMPONENT, SURFACE_FACTOR, N_LAYER_WFS,
     O       R2_L_BEAM, R2_L_HOMP, R2_L_HOMM )

C  start loops

         DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES
           IROW = IR + O1
           CM = C0 + IROW
           DO Q = 1, N_LAYER_WFS

C  Beam contributions

            L_BEAM = L_WLOWER(I1,O1,N,Q) - R2_L_BEAM(I,O1,Q)

C  Linearized Real homogeneous solution contributions

            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N)
              CPOS =  T_DELT_EIGEN(K,N)   * L_SOLA_XPOS(I1,O1,K,N,Q) 
     &            + L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I1,O1,K,N)
              CPOS = CPOS - R2_L_HOMP(I,O1,K,Q)
              CNEG = L_SOLB_XNEG(I1,O1,K,N,Q)
              CNEG = CNEG - R2_L_HOMM(I,O1,K,Q)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

C  Linearized Complex homogeneous solution contributions
C    Bug Fixed 16 December 2005.

            L_HOM_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T1R = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q)
     &            - T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q)
     &            + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
     &            - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
              T1I = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q)
     &            + T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q)
     &            + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
     &            + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
              T1R = T1R - R2_L_HOMP(I,O1,K1,Q)
              T1I = T1I - R2_L_HOMP(I,O1,K2,Q)
              T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
              T2R = L_SOLB_XNEG(I1,O1,K1,N,Q) - R2_L_HOMM(I,O1,K1,Q)
              T2I = L_SOLB_XNEG(I1,O1,K2,N,Q) - R2_L_HOMM(I,O1,K2,Q)
              T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

C  Final contributions

            COLTEL2_WF(CM,Q) = - L_BEAM - L_HOM_R - L_HOM_CR

C  End loops

           ENDDO
          ENDDO
         ENDDO

C  If the varying layer is above the surface layer then

        ELSE

C  get the linearized downward-reflected term

         CALL L_BVP_SURFACE_SETUP
     I     ( DO_INCLUDE_SURFACE, .FALSE., IBEAM,
     I       FOURIER_COMPONENT, SURFACE_FACTOR, N_LAYER_WFS,
     O       R2_L_BEAM, R2_L_HOMP, R2_L_HOMM )

C  start loops : beam contributions only

         DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES
           IROW = IR + O1
           CM = C0 + IROW
           DO Q = 1, N_LAYER_WFS
             L_BEAM = L_WLOWER(I1,O1,N,Q) - R2_L_BEAM(I,O1,Q)
             COLTEL2_WF(CM,Q) = - L_BEAM
           ENDDO
          ENDDO
         ENDDO

        ENDIF

C  Otherwise, there's no surface to consider.............

      ELSE

C  If active layer is varying layer, need full calculation.

       IF ( N .EQ. LAYER_TO_VARY ) THEN

C  start loops

        DO I = 1, NSTREAMS
         I1 = I + NSTREAMS
         IR = NSTOKES*(I-1)
         DO O1 = 1, NSTOKES
          IROW = IR + O1
          CM = C0 + IROW
          DO Q = 1, N_LAYER_WFS

C  Beam contributions

           L_BEAM = - L_WLOWER(I1,O1,N,Q)

C  Linearized Real homogeneous solution contributions

           L_HOM_R  = ZERO
           DO K = 1, K_REAL(N)
             CNEG = L_SOLB_XNEG(I1,O1,K,N,Q)
             CPOS = T_DELT_EIGEN(K,N) * L_SOLA_XPOS(I1,O1,K,N,Q) + 
     &            L_T_DELT_EIGEN(K,N,Q) * SOLA_XPOS(I1,O1,K,N)
             T1 = LCON(K,N) * CPOS
             T2 = MCON(K,N) * CNEG
             L_HOM_R = L_HOM_R + T1 + T2
           ENDDO

C  Linearized Complex homogeneous solution contributions

           L_HOM_CR  = ZERO
           KO1 = K_REAL(N) + 1
           DO K = 1, K_COMPLEX(N)
             K0 = 2*K - 2
             K1 = KO1 + K0
             K2 = K1  + 1
             T2 = L_SOLB_XNEG(I1,O1,K1,N,Q) * MCON(K1,N) -
     &            L_SOLB_XNEG(I1,O1,K2,N,Q) * MCON(K2,N)
             T1R =  T_DELT_EIGEN(K1,N)  * L_SOLA_XPOS(I1,O1,K1,N,Q)
     &            - T_DELT_EIGEN(K2,N)  * L_SOLA_XPOS(I1,O1,K2,N,Q)
     &            + L_T_DELT_EIGEN(K1,N,Q) *  SOLA_XPOS(I1,O1,K1,N)
     &            - L_T_DELT_EIGEN(K2,N,Q) *  SOLA_XPOS(I1,O1,K2,N)
             T1I =  T_DELT_EIGEN(K1,N)  * L_SOLA_XPOS(I1,O1,K2,N,Q)
     &            + T_DELT_EIGEN(K2,N)  * L_SOLA_XPOS(I1,O1,K1,N,Q)
     &            + L_T_DELT_EIGEN(K1,N,Q) *  SOLA_XPOS(I1,O1,K2,N)
     &            + L_T_DELT_EIGEN(K2,N,Q) *  SOLA_XPOS(I1,O1,K1,N)
             T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
             L_HOM_CR = L_HOM_CR + T1 + T2
           ENDDO

C  Final contribution

           COLTEL2_WF(CM,Q) = L_BEAM - L_HOM_R - L_HOM_CR

C  End loops

          ENDDO
         ENDDO
        ENDDO

C  otherwise use beam solution linearizations propagated downwards
C  from the layer that is varying (already computed for this layer)

       ELSE IF ( LAYER_TO_VARY .LT. N ) THEN

C  start loops

        DO I = 1, NSTREAMS
         I1 = I + NSTREAMS
         IR = NSTOKES*(I-1)
         DO O1 = 1, NSTOKES
          IROW = IR + O1
          CM = C0 + IROW
          DO Q = 1, N_LAYER_WFS
            COLTEL2_WF(CM,Q) = - L_WLOWER(I1,O1,N,Q)
          ENDDO
         ENDDO
        ENDDO

       ENDIF

C  Finally end the surface condition

      ENDIF

C  Add direct beam variation to Final boundary
C  -------------------------------------------

C  Only for the Specialist option # 2, and only if the surface
C   is included (Fourier = 0) and the lowest active layer is
C   the surface layer

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        IF ( DO_INCLUDE_SURFACE.AND.DO_SPECIALIST_OPTION_2 ) THEN
         IF ( ACTIVE_LAYERS(NLAYERS_TEL).EQ.NLAYERS ) THEN
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              FAC = - DIRECT_BEAM(I,IBEAM,O1) * 
     &                DELTAU_SLANT(N,LAYER_TO_VARY,IBEAM) 
              DO Q = 1, N_LAYER_WFS
                L_BEAM = L_DELTAU_VERT(Q,LAYER_TO_VARY) * FAC
                COLTEL2_WF(CM,Q) = COLTEL2_WF(CM,Q) + L_BEAM
              ENDDO
            ENDDO
          ENDDO
         ENDIF
        ENDIF
      ENDIF

C  Continuation point for single layer stuffs
C  ==========================================

 3456 CONTINUE

C  zero column vector. Extremely important.

      DO I = 1, NSTKS_NSTRMS_2
        DO Q = 1, MAX_ATMOSWFS
          SCOL2_WF(I,Q) = ZERO
        ENDDO
      ENDDO

C  top of active layer
C  -------------------

C  If active layer = layer that is varying,
C       then require homogeneous and beam solution linearizations

      IF ( LAYER_TO_VARY .EQ. N ) THEN

C  Start  loops

        DO I = 1, NSTREAMS
         IR = NSTOKES*(I-1)
         DO O1 = 1, NSTOKES
          IROW = IR + O1
          DO Q = 1, N_LAYER_WFS

C  beam solution linearization at top of layer

            L_BEAM = - L_WUPPER(I,O1,N,Q)

C  Linearized Real homogeneous solution contributions

            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N)
              CPOS = L_SOLA_XPOS(I,O1,K,N,Q)
              CNEG = T_DELT_EIGEN(K,N)   * L_SOLB_XNEG(I,O1,K,N,Q) + 
     &             L_T_DELT_EIGEN(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

C  Linearized Complex homogeneous solution contributions

            L_HOM_CR  = ZERO
             KO1 = K_REAL(N) + 1
             DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T1 = L_SOLA_XPOS(I,O1,K1,N,Q) * LCON(K1,N) -
     &             L_SOLA_XPOS(I,O1,K2,N,Q) * LCON(K2,N)
              T2R =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K1,N,Q)
     &             - T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K2,N,Q)
     &             + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N)
     &             - L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)
              T2I =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K2,N,Q)
     &             + T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K1,N,Q)
     &             + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K2,N)
     &             + L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K1,N)
              T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

C  Final contribution

            SCOL2_WF(IROW,Q) = L_BEAM - L_HOM_R - L_HOM_CR

C  end loops

          ENDDO
         ENDDO
        ENDDO

C  otherwise if varying layer is above active layer, there are beam
C  solution contributions propagated downwards - find these by calling
C  the appropriate solution module = L_BEAMSOLUTION_NNEK

      ELSE IF ( LAYER_TO_VARY .LT. N ) THEN

        CALL L_BEAMSOLUTION_P_NNEK
     &    ( FOURIER_COMPONENT, IBEAM, N, LAYER_TO_VARY, N_LAYER_WFS )

        DO I = 1, NSTREAMS
         IR = NSTOKES*(I-1)
         DO O1 = 1, NSTOKES
          IROW = IR + O1
          DO Q = 1, N_LAYER_WFS
            SCOL2_WF(IROW,Q) = - L_WUPPER(I,O1,N,Q)
          ENDDO
         ENDDO
        ENDDO

      ENDIF

C  Bottom of active layer
C  ----------------------

      C0 = NSTKS_NSTRMS

C  If this is the  surface and Specialist option #2 is in place

      if ( DO_INCLUDE_SURFACE.AND.DO_SPECIALIST_OPTION_2
     &       .AND.  N.EQ.NLAYERS ) THEN

C  Modified BCL4M Component loop

        IF ( LAYER_TO_VARY .EQ. N ) THEN

C  get the linearized downward-reflected term

         CALL L_BVP_SURFACE_SETUP
     I     ( DO_INCLUDE_SURFACE, .TRUE., IBEAM,
     I       FOURIER_COMPONENT, SURFACE_FACTOR, N_LAYER_WFS,
     O       R2_L_BEAM, R2_L_HOMP, R2_L_HOMM )

C  start loops

         DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES
           IROW = IR + O1
           CM = C0 + IROW
           DO Q = 1, N_LAYER_WFS

C  Beam contributions

            L_BEAM = L_WLOWER(I1,O1,N,Q) - R2_L_BEAM(I,O1,Q)

C  Linearized Real homogeneous solution contributions

            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N)
              CPOS =  T_DELT_EIGEN(K,N)   * L_SOLA_XPOS(I1,O1,K,N,Q) 
     &            + L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I1,O1,K,N)
              CPOS = CPOS - R2_L_HOMP(I,O1,K,Q)
              CNEG = L_SOLB_XNEG(I1,O1,K,N,Q)
              CNEG = CNEG - R2_L_HOMM(I,O1,K,Q)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

C  Linearized Complex homogeneous solution contributions
C    Bug Fixed 16 December 2005.

            L_HOM_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T1R = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q)
     &            - T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q)
     &            + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
     &            - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
              T1I = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q)
     &            + T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q)
     &            + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
     &            + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
              T1R = T1R - R2_L_HOMP(I,O1,K1,Q)
              T1I = T1I - R2_L_HOMP(I,O1,K2,Q)
              T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
              T2R = L_SOLB_XNEG(I1,O1,K1,N,Q) - R2_L_HOMM(I,O1,K1,Q)
              T2I = L_SOLB_XNEG(I1,O1,K2,N,Q) - R2_L_HOMM(I,O1,K2,Q)
              T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

C  Final contributions

            SCOL2_WF(CM,Q) = - L_BEAM - L_HOM_R - L_HOM_CR

C  End loops

           ENDDO
          ENDDO
         ENDDO

C  If the varying layer is above the surface layer then

        ELSE

C  get the linearized downward-reflected term

         CALL L_BVP_SURFACE_SETUP
     I     ( DO_INCLUDE_SURFACE, .FALSE., IBEAM,
     I       FOURIER_COMPONENT, SURFACE_FACTOR, N_LAYER_WFS,
     O       R2_L_BEAM, R2_L_HOMP, R2_L_HOMM )

C  start loops : beam contributions only

         DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES
           IROW = IR + O1
           CM = C0 + IROW
           DO Q = 1, N_LAYER_WFS
             L_BEAM = L_WLOWER(I1,O1,N,Q) - R2_L_BEAM(I,O1,Q)
             SCOL2_WF(CM,Q) = - L_BEAM
           ENDDO
          ENDDO
         ENDDO

        ENDIF

C  Otherwise, there's no surface to consider.............

      ELSE

C  If active layer is varying layer, need full calculation.

        IF ( N .EQ. LAYER_TO_VARY ) THEN
 
C  start loops

         DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
           IROW = IR + O1
           CM = C0 + IROW
           DO Q = 1, N_LAYER_WFS

C  Beam contributions

            L_BEAM = - L_WLOWER(I1,O1,N,Q)

C  Linearized Real homogeneous solution contributions

            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N)
              CNEG = L_SOLB_XNEG(I1,O1,K,N,Q)
              CPOS = T_DELT_EIGEN(K,N) * L_SOLA_XPOS(I1,O1,K,N,Q) + 
     &             L_T_DELT_EIGEN(K,N,Q) * SOLA_XPOS(I1,O1,K,N)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

C  Linearized Complex homogeneous solution contributions

            L_HOM_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T2 = L_SOLB_XNEG(I1,O1,K1,N,Q) * MCON(K1,N) -
     &             L_SOLB_XNEG(I1,O1,K2,N,Q) * MCON(K2,N)
              T1R =  T_DELT_EIGEN(K1,N)  * L_SOLA_XPOS(I1,O1,K1,N,Q)
     &             - T_DELT_EIGEN(K2,N)  * L_SOLA_XPOS(I1,O1,K2,N,Q)
     &             + L_T_DELT_EIGEN(K1,N,Q) *  SOLA_XPOS(I1,O1,K1,N)
     &             - L_T_DELT_EIGEN(K2,N,Q) *  SOLA_XPOS(I1,O1,K2,N)
              T1I =  T_DELT_EIGEN(K1,N)  * L_SOLA_XPOS(I1,O1,K2,N,Q)
     &             + T_DELT_EIGEN(K2,N)  * L_SOLA_XPOS(I1,O1,K1,N,Q)
     &             + L_T_DELT_EIGEN(K1,N,Q) *  SOLA_XPOS(I1,O1,K2,N)
     &             + L_T_DELT_EIGEN(K2,N,Q) *  SOLA_XPOS(I1,O1,K1,N)
              T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

C  Final contribution

            SCOL2_WF(CM,Q) = L_BEAM - L_HOM_R - L_HOM_CR

C  End loops

           ENDDO
          ENDDO
         ENDDO

C  otherwise use beam solution linearizations propagated downwards
C  from the layer that is varying (already computed for this layer)

        ELSE IF ( LAYER_TO_VARY .LT. N ) THEN

C  start loops

         DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
           IROW = IR + O1
           CM = C0 + IROW
           DO Q = 1, N_LAYER_WFS
             SCOL2_WF(CM,Q) = - L_WLOWER(I1,O1,N,Q)
           ENDDO
          ENDDO
         ENDDO

        ENDIF

C  End the surface condition

      ENDIF

C  Add direct beam variation to Final boundary
C  -------------------------------------------

C  Only for the Specialist option # 2, and only if the surface
C   is included (Fourier = 0) and the lowest active layer is
C   the surface layer

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        IF ( DO_INCLUDE_SURFACE.AND.DO_SPECIALIST_OPTION_2 ) THEN
         IF ( ACTIVE_LAYERS(NLAYERS_TEL).EQ.NLAYERS ) THEN
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              FAC = - DIRECT_BEAM(I,IBEAM,O1) * 
     &                DELTAU_SLANT(N,LAYER_TO_VARY,IBEAM) 
              DO Q = 1, N_LAYER_WFS
                L_BEAM = L_DELTAU_VERT(Q,LAYER_TO_VARY) * FAC
                SCOL2_WF(CM,Q) = SCOL2_WF(CM,Q) + L_BEAM
              ENDDO
            ENDDO
          ENDDO
         ENDIF
        ENDIF
      ENDIF

C  finish

      RETURN
      END

C

      SUBROUTINE L_BVPTEL_C_COLUMN_SETUP
     I     ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM,
     I       N_LAYER_WFS,
     I       FOURIER_COMPONENT,  IBEAM,
     I       SURFACE_FACTOR,
     I       L_BVP_SURFACE_SETUP )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup, solution variables

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'

C  include files of linearized setup variables (input)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'

C  include files of linearized solution variables (output stored here)

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'

C  Include file of reflectance variables

      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'

C  input arguments
C  ---------------

C  inclusion flags

      LOGICAL          DO_INCLUDE_SURFACE
      LOGICAL          DO_INCLUDE_DIRECTBEAM

C  Fourier component and beam number

      INTEGER          FOURIER_COMPONENT, IBEAM

C  surface factor = 1+delta(m,0)

      DOUBLE PRECISION SURFACE_FACTOR

C  Linearization control

      INTEGER          N_LAYER_WFS

C  External function

      EXTERNAL         L_BVP_SURFACE_SETUP

C  output status
C  -------------

      INTEGER          STATUS

C  local variables
C  ---------------

      INTEGER          Q, N, I, I1, IR, CM, C0, IROW, O1, N1, NS
      INTEGER          K, KO1, K0, K1, K2
      DOUBLE PRECISION CPOS,CNEG,L_HOM_R,L_HOM_CR,L_BEAM
      DOUBLE PRECISION L_HOM_U_R, L_HOM_U_CR, L_HOM_D_R, L_HOM_D_CR
      DOUBLE PRECISION T1,T2,T1R,T1I,T2R,T2I,FAC,FAC3

      DOUBLE PRECISION R2_L_BEAM
     &        ( MAXSTREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION R2_L_HOMP
     &        ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )
      DOUBLE PRECISION R2_L_HOMM
     &        ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )

C  status

      status = 0

C  Try this safety-first zeroing, copied from the scalar LIDORT code.
C    Very important to do this.   Bug solved, July 14th 2009

      DO NS = 1, NLAYERS_TEL
        N = ACTIVE_LAYERS(NS)
        DO I = 1, NSTREAMS_2
          DO Q = 1, N_LAYER_WFS
            DO O1 = 1, NSTOKES
              L_WUPPER(I,O1,N,Q) = ZERO
              L_WLOWER(I,O1,N,Q) = ZERO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

C  Get the linearized solutions for the layer that is varying
C    Always need this, regardless of number of active layers

      DO NS = 1, NLAYERS_TEL
        N = ACTIVE_LAYERS(NS)
        CALL L_BEAMSOLUTION_C_NEQK
     I      ( FOURIER_COMPONENT, IBEAM, N, N_LAYER_WFS )
      ENDDO

C  Go to special case for only 1 active layer

      IF ( NLAYERS_TEL .EQ. 1 ) GO TO 3456

C  zero column vector

      DO I = 1, N_BVTELMATRIX_SIZE
        DO Q = 1, MAX_ATMOSWFS
          COLTEL2_WF(I,Q) = ZERO
        ENDDO
      ENDDO

C  top of first active layer, first boundary condition
C  ---------------------------------------------------

      NS = 1
      N = ACTIVE_LAYERS(NS)
      C0 = 0

C  require homogeneous and beam solution linearizations


      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        DO O1 = 1, NSTOKES
          IROW = IR + O1
          CM = C0 + IROW
          DO Q = 1, N_LAYER_WFS

C  Beam contribution

            L_BEAM  = - L_WUPPER(I,O1,N,Q)

C  Linearized Real homogeneous solution contributions

            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N)
              CPOS = L_SOLA_XPOS(I,O1,K,N,Q)
              CNEG = T_DELT_EIGEN(K,N)   * L_SOLB_XNEG(I,O1,K,N,Q) + 
     &             L_T_DELT_EIGEN(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

C  Linearized Complex homogeneous solution contributions

            L_HOM_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T1 = L_SOLA_XPOS(I,O1,K1,N,Q) * LCON(K1,N) -
     &             L_SOLA_XPOS(I,O1,K2,N,Q) * LCON(K2,N)
              T2R =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K1,N,Q)
     &             - T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K2,N,Q)
     &             + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N)
     &             - L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)
              T2I =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K2,N,Q)
     &             + T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K1,N,Q)
     &             + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K2,N)
     &             + L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K1,N)
              T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

C  Final contribution

            COLTEL2_WF(CM,Q) = L_BEAM - L_HOM_R - L_HOM_CR

C  end loops

          ENDDO
        ENDDO
      ENDDO

C  Intermediate boundaries between active layers
C  ---------------------------------------------

      DO NS = 1, NLAYERS_TEL - 1

C  offsets

       N  = ACTIVE_LAYERS(NS)
       N1 = N + 1
       C0 = NS*NSTKS_NSTRMS_2 - NSTKS_NSTRMS

C  Get the linearized beam solution for the next layer N1

       DO I = 1, NSTREAMS_2
         IR = NSTOKES*(I-1)
         DO O1 = 1, NSTOKES
           IROW = IR + O1
           CM = C0 + IROW
           DO Q = 1, N_LAYER_WFS

C  Beam contributions

            L_BEAM = L_WUPPER(I,O1,N1,Q) - L_WLOWER(I,O1,N,Q)

C  Linearized Real homogeneous solution contributions

c            L_HOM_R  = ZERO
c            DO K = 1, K_REAL(N)
c              CNEG = L_SOLB_XNEG(I,O1,K,N,Q)
c              CPOS = T_DELT_EIGEN(K,N)   * L_SOLA_XPOS(I,O1,K,N,Q) + 
c     &             L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I,O1,K,N)
c              T1 = LCON(K,N) * CPOS
c              T2 = MCON(K,N) * CNEG
c              L_HOM_R = L_HOM_R + T1 + T2
c            ENDDO

C  Linearized Real homogeneous solution contributions above

            L_HOM_U_R = ZERO
            DO K = 1, K_REAL(N1)
              CPOS = L_SOLA_XPOS(I,O1,K,N1,Q)
              CNEG = T_DELT_EIGEN(K,N1)   * L_SOLB_XNEG(I,O1,K,N1,Q) + 
     &             L_T_DELT_EIGEN(K,N1,Q) *   SOLB_XNEG(I,O1,K,N1)
              T1 = LCON(K,N1) * CPOS
              T2 = MCON(K,N1) * CNEG
              L_HOM_U_R = L_HOM_U_R + T1 + T2
            ENDDO

C  Linearized Real homogeneous solution contributions below

            L_HOM_D_R = ZERO
            DO K = 1, K_REAL(N)
              CNEG = L_SOLB_XNEG(I,O1,K,N,Q)
              CPOS = T_DELT_EIGEN(K,N)   * L_SOLA_XPOS(I,O1,K,N,Q) + 
     &             L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I,O1,K,N)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_D_R = L_HOM_D_R + T1 + T2
            ENDDO

C  Linearized Complex homogeneous solution contributions above

            L_HOM_U_CR  = ZERO
            KO1 = K_REAL(N1) + 1
            DO K = 1, K_COMPLEX(N1)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T1 = L_SOLA_XPOS(I,O1,K1,N1,Q) * LCON(K1,N1) -
     &             L_SOLA_XPOS(I,O1,K2,N1,Q) * LCON(K2,N1)
              T2R =  T_DELT_EIGEN(K1,N1)   * L_SOLB_XNEG(I,O1,K1,N1,Q)
     &             - T_DELT_EIGEN(K2,N1)   * L_SOLB_XNEG(I,O1,K2,N1,Q)
     &             + L_T_DELT_EIGEN(K1,N1,Q) *   SOLB_XNEG(I,O1,K1,N1)
     &             - L_T_DELT_EIGEN(K2,N1,Q) *   SOLB_XNEG(I,O1,K2,N1)
              T2I =  T_DELT_EIGEN(K1,N1)   * L_SOLB_XNEG(I,O1,K2,N1,Q)
     &             + T_DELT_EIGEN(K2,N1)   * L_SOLB_XNEG(I,O1,K1,N1,Q)
     &             + L_T_DELT_EIGEN(K1,N1,Q) *   SOLB_XNEG(I,O1,K2,N1)
     &             + L_T_DELT_EIGEN(K2,N1,Q) *   SOLB_XNEG(I,O1,K1,N1)
              T2 =  T2R * MCON(K1,N1) - T2I * MCON(K2,N1)
              L_HOM_U_CR = L_HOM_U_CR + T1 + T2
            ENDDO

C  Linearized Complex homogeneous solution contributions below

            L_HOM_D_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T2 = L_SOLB_XNEG(I,O1,K1,N,Q) * MCON(K1,N) -
     &             L_SOLB_XNEG(I,O1,K2,N,Q) * MCON(K2,N)
              T1R =  T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I,O1,K1,N,Q)
     &             - T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I,O1,K2,N,Q)
     &             + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I,O1,K1,N)
     &             - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I,O1,K2,N)
              T1I =  T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I,O1,K2,N,Q)
     &             + T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I,O1,K1,N,Q)
     &             + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I,O1,K2,N)
     &             + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I,O1,K1,N)
              T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
              L_HOM_D_CR = L_HOM_D_CR + T1 + T2
            ENDDO

C  Linearized Complex homogeneous solution contributions

c            L_HOM_CR  = ZERO
c            KO1 = K_REAL(N) + 1
c            DO K = 1, K_COMPLEX(N)
c              K0 = 2*K - 2
c              K1 = KO1 + K0
c              K2 = K1  + 1
c              T2 = L_SOLB_XNEG(I,O1,K1,N,Q) * MCON(K1,N) -
c     &             L_SOLB_XNEG(I,O1,K2,N,Q) * MCON(K2,N)
c              T1R =  T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I,O1,K1,N,Q)
c     &             - T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I,O1,K2,N,Q)
c     &             + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I,O1,K1,N)
c     &             - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I,O1,K2,N)
c              T1I =  T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I,O1,K2,N,Q)
c     &             + T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I,O1,K1,N,Q)
c     &             + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I,O1,K2,N)
c     &             + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I,O1,K1,N)
c              T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
c              L_HOM_CR = L_HOM_CR + T1 + T2
c            ENDDO

C  Final contribution

            L_HOM_R  =  L_HOM_U_R  - L_HOM_D_R
            L_HOM_CR =  L_HOM_U_CR - L_HOM_D_CR
            COLTEL2_WF(CM,Q) = L_BEAM + L_HOM_R + L_HOM_CR

c            COLTEL2_WF(CM,Q) = L_BEAM - L_HOM_R - L_HOM_CR

C  End loops

           ENDDO
         ENDDO
       ENDDO

C  End loop over intermediate active layer boundaries

      ENDDO

C  Final boundary, bottom of lowest active layer
C  ---------------------------------------------

      NS = NLAYERS_TEL
      N  = ACTIVE_LAYERS(NS)      
      C0 = (NS-1)*NSTKS_NSTRMS_2 + NSTKS_NSTRMS

C  If this is the surface and Specialist option #2 is in place

      if ( DO_INCLUDE_SURFACE.AND.DO_SPECIALIST_OPTION_2
     &       .AND.  N.EQ.NLAYERS ) THEN

C  get the linearized downward-reflected term

        CALL L_BVP_SURFACE_SETUP
     I     ( DO_INCLUDE_SURFACE, .TRUE., IBEAM,
     I       FOURIER_COMPONENT, SURFACE_FACTOR, N_LAYER_WFS,
     O       R2_L_BEAM, R2_L_HOMP, R2_L_HOMM )

C  start loops

        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            DO Q = 1, N_LAYER_WFS

C  Beam contributions

             L_BEAM = L_WLOWER(I1,O1,N,Q) - R2_L_BEAM(I,O1,Q)

C  Linearized Real homogeneous solution contributions

             L_HOM_R  = ZERO
             DO K = 1, K_REAL(N)
              CPOS =  T_DELT_EIGEN(K,N)   * L_SOLA_XPOS(I1,O1,K,N,Q) 
     &            + L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I1,O1,K,N)
              CPOS = CPOS - R2_L_HOMP(I,O1,K,Q)
              CNEG = L_SOLB_XNEG(I1,O1,K,N,Q)
              CNEG = CNEG - R2_L_HOMM(I,O1,K,Q)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
             ENDDO

C  Linearized Complex homogeneous solution contributions
C    Bug Fixed 16 December 2005.

             L_HOM_CR  = ZERO
             KO1 = K_REAL(N) + 1
             DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T1R = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q)
     &            - T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q)
     &            + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
     &            - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
              T1I = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q)
     &            + T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q)
     &            + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
     &            + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
              T1R = T1R - R2_L_HOMP(I,O1,K1,Q)
              T1I = T1I - R2_L_HOMP(I,O1,K2,Q)
              T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
              T2R = L_SOLB_XNEG(I1,O1,K1,N,Q) - R2_L_HOMM(I,O1,K1,Q)
              T2I = L_SOLB_XNEG(I1,O1,K2,N,Q) - R2_L_HOMM(I,O1,K2,Q)
              T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
              L_HOM_CR = L_HOM_CR + T1 + T2
             ENDDO

C  Final contributions

             COLTEL2_WF(CM,Q) = - L_BEAM - L_HOM_R - L_HOM_CR

C  End loops

           ENDDO
         ENDDO
       ENDDO

C  Otherwise, there's no surface to consider.............

      ELSE

C  start loops

       DO I = 1, NSTREAMS
         I1 = I + NSTREAMS
         IR = NSTOKES*(I-1)
         DO O1 = 1, NSTOKES
           IROW = IR + O1
           CM = C0 + IROW
           DO Q = 1, N_LAYER_WFS

C  Beam contributions

            L_BEAM = - L_WLOWER(I1,O1,N,Q)

C  Linearized Real homogeneous solution contributions

            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N)
             CNEG = L_SOLB_XNEG(I1,O1,K,N,Q)
             CPOS = T_DELT_EIGEN(K,N) * L_SOLA_XPOS(I1,O1,K,N,Q) + 
     &            L_T_DELT_EIGEN(K,N,Q) * SOLA_XPOS(I1,O1,K,N)
             T1 = LCON(K,N) * CPOS
             T2 = MCON(K,N) * CNEG
             L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

C  Linearized Complex homogeneous solution contributions

            L_HOM_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
             K0 = 2*K - 2
             K1 = KO1 + K0
             K2 = K1  + 1
             T2 = L_SOLB_XNEG(I1,O1,K1,N,Q) * MCON(K1,N) -
     &            L_SOLB_XNEG(I1,O1,K2,N,Q) * MCON(K2,N)
             T1R =  T_DELT_EIGEN(K1,N)  * L_SOLA_XPOS(I1,O1,K1,N,Q)
     &            - T_DELT_EIGEN(K2,N)  * L_SOLA_XPOS(I1,O1,K2,N,Q)
     &            + L_T_DELT_EIGEN(K1,N,Q) *  SOLA_XPOS(I1,O1,K1,N)
     &            - L_T_DELT_EIGEN(K2,N,Q) *  SOLA_XPOS(I1,O1,K2,N)
             T1I =  T_DELT_EIGEN(K1,N)  * L_SOLA_XPOS(I1,O1,K2,N,Q)
     &            + T_DELT_EIGEN(K2,N)  * L_SOLA_XPOS(I1,O1,K1,N,Q)
     &            + L_T_DELT_EIGEN(K1,N,Q) *  SOLA_XPOS(I1,O1,K2,N)
     &            + L_T_DELT_EIGEN(K2,N,Q) *  SOLA_XPOS(I1,O1,K1,N)
             T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
             L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

C  Final contribution

            COLTEL2_WF(CM,Q) = L_BEAM - L_HOM_R - L_HOM_CR

C  End loops

           ENDDO
         ENDDO
       ENDDO

C  Finally end the surface condition

      ENDIF

C  Add direct beam variation to Final boundary
C  -------------------------------------------

C  Only for the Specialist option # 2, and only if the surface
C   is included (Fourier = 0) and the lowest active layer is
C   the surface layer

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        IF ( DO_INCLUDE_SURFACE.AND.DO_SPECIALIST_OPTION_2 ) THEN
         IF ( ACTIVE_LAYERS(NLAYERS_TEL).EQ.NLAYERS ) THEN
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              FAC = - DIRECT_BEAM(I,IBEAM,O1)
              DO Q = 1, N_LAYER_WFS
                L_BEAM = ZERO
                DO K = 1, NLAYERS
                  FAC3 = FAC * DELTAU_SLANT(N,K,IBEAM)
                  L_BEAM = L_BEAM + L_DELTAU_VERT(Q,K) * FAC3
                ENDDO
                COLTEL2_WF(CM,Q) = COLTEL2_WF(CM,Q) + L_BEAM
              ENDDO
            ENDDO
          ENDDO
         ENDIF
        ENDIF
      ENDIF

C  DEbug

c      do cm = 1, NSTREAMS_2 * NLAYERS_TEL
c        write(*,*)COLTEL2_WF(CM,1)
c      enddo
c      pause

C  Can return now

      RETURN

C  Continuation point for single layer stuffs
C  ==========================================

 3456 CONTINUE

C  zero column vector. Extremely important.

      DO I = 1, NSTKS_NSTRMS_2
        DO Q = 1, MAX_ATMOSWFS
          SCOL2_WF(I,Q) = ZERO
        ENDDO
      ENDDO

C  top of active layer
C  -------------------

      NS = 1
      N = ACTIVE_LAYERS(NS)

C  layer that is varying,
C       then require homogeneous and beam solution linearizations

      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        DO O1 = 1, NSTOKES
          IROW = IR + O1
          DO Q = 1, N_LAYER_WFS

C  beam solution linearization at top of layer

            L_BEAM = - L_WUPPER(I,O1,N,Q)

C  Linearized Real homogeneous solution contributions

            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N)
              CPOS = L_SOLA_XPOS(I,O1,K,N,Q)
              CNEG = T_DELT_EIGEN(K,N)   * L_SOLB_XNEG(I,O1,K,N,Q) + 
     &             L_T_DELT_EIGEN(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

C  Linearized Complex homogeneous solution contributions

            L_HOM_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T1 = L_SOLA_XPOS(I,O1,K1,N,Q) * LCON(K1,N) -
     &             L_SOLA_XPOS(I,O1,K2,N,Q) * LCON(K2,N)
              T2R =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K1,N,Q)
     &             - T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K2,N,Q)
     &             + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N)
     &             - L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)
              T2I =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K2,N,Q)
     &             + T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K1,N,Q)
     &             + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K2,N)
     &             + L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K1,N)
              T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

C  Final contribution

            SCOL2_WF(IROW,Q) = L_BEAM - L_HOM_R - L_HOM_CR

C  end loops

          ENDDO
        ENDDO
      ENDDO

C  Bottom of active layer
C  ----------------------

      C0 = NSTKS_NSTRMS

C  active layer is varying layer

C  If this is the  surface and Specialist option #2 is in place

      if ( DO_INCLUDE_SURFACE.AND.DO_SPECIALIST_OPTION_2
     &       .AND.  N.EQ.NLAYERS ) THEN

C  get the linearized downward-reflected term

       CALL L_BVP_SURFACE_SETUP
     I     ( DO_INCLUDE_SURFACE, .TRUE., IBEAM,
     I       FOURIER_COMPONENT, SURFACE_FACTOR, N_LAYER_WFS,
     O       R2_L_BEAM, R2_L_HOMP, R2_L_HOMM )

       DO I = 1, NSTREAMS
         IR = NSTOKES*(I-1)
         I1 = I + NSTREAMS
         DO O1 = 1, NSTOKES
           IROW = IR + O1
           CM = C0 + IROW
           DO Q = 1, N_LAYER_WFS

C  Beam contributions

            L_BEAM = L_WLOWER(I1,O1,N,Q) - R2_L_BEAM(I,O1,Q)

C  Linearized Real homogeneous solution contributions

            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N)
              CPOS =  T_DELT_EIGEN(K,N)   * L_SOLA_XPOS(I1,O1,K,N,Q) 
     &            + L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I1,O1,K,N)
              CPOS = CPOS - R2_L_HOMP(I,O1,K,Q)
              CNEG = L_SOLB_XNEG(I1,O1,K,N,Q)
              CNEG = CNEG - R2_L_HOMM(I,O1,K,Q)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

C  Linearized Complex homogeneous solution contributions
C    Bug Fixed 16 December 2005.

            L_HOM_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T1R = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q)
     &            - T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q)
     &            + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
     &            - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
              T1I = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q)
     &            + T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q)
     &            + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
     &            + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
              T1R = T1R - R2_L_HOMP(I,O1,K1,Q)
              T1I = T1I - R2_L_HOMP(I,O1,K2,Q)
              T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
              T2R = L_SOLB_XNEG(I1,O1,K1,N,Q) - R2_L_HOMM(I,O1,K1,Q)
              T2I = L_SOLB_XNEG(I1,O1,K2,N,Q) - R2_L_HOMM(I,O1,K2,Q)
              T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

C  Final contributions

            SCOL2_WF(CM,Q) = - L_BEAM - L_HOM_R - L_HOM_CR

C  End loops

           ENDDO
         ENDDO
       ENDDO

C  Otherwise, there's no surface to consider.............

      ELSE

       DO I = 1, NSTREAMS
         I1 = I + NSTREAMS
         IR = NSTOKES*(I-1)
         DO O1 = 1, NSTOKES
           IROW = IR + O1
           CM = C0 + IROW
           DO Q = 1, N_LAYER_WFS

C  Beam contributions

            L_BEAM = - L_WLOWER(I1,O1,N,Q)

C  Linearized Real homogeneous solution contributions

            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N)
              CNEG = L_SOLB_XNEG(I1,O1,K,N,Q)
              CPOS = T_DELT_EIGEN(K,N) * L_SOLA_XPOS(I1,O1,K,N,Q) + 
     &             L_T_DELT_EIGEN(K,N,Q) * SOLA_XPOS(I1,O1,K,N)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

C  Linearized Complex homogeneous solution contributions

            L_HOM_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T2 = L_SOLB_XNEG(I1,O1,K1,N,Q) * MCON(K1,N) -
     &             L_SOLB_XNEG(I1,O1,K2,N,Q) * MCON(K2,N)
              T1R =  T_DELT_EIGEN(K1,N)  * L_SOLA_XPOS(I1,O1,K1,N,Q)
     &             - T_DELT_EIGEN(K2,N)  * L_SOLA_XPOS(I1,O1,K2,N,Q)
     &             + L_T_DELT_EIGEN(K1,N,Q) *  SOLA_XPOS(I1,O1,K1,N)
     &             - L_T_DELT_EIGEN(K2,N,Q) *  SOLA_XPOS(I1,O1,K2,N)
              T1I =  T_DELT_EIGEN(K1,N)  * L_SOLA_XPOS(I1,O1,K2,N,Q)
     &             + T_DELT_EIGEN(K2,N)  * L_SOLA_XPOS(I1,O1,K1,N,Q)
     &             + L_T_DELT_EIGEN(K1,N,Q) *  SOLA_XPOS(I1,O1,K2,N)
     &             + L_T_DELT_EIGEN(K2,N,Q) *  SOLA_XPOS(I1,O1,K1,N)
              T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

C  Final contribution

            SCOL2_WF(CM,Q) = L_BEAM - L_HOM_R - L_HOM_CR

C  End loops

           ENDDO
         ENDDO
       ENDDO

C  End the surface condition

      ENDIF

C  Add direct beam variation to Final boundary
C  -------------------------------------------

C  Only for the Specialist option # 2, and only if the surface
C   is included (Fourier = 0) and the lowest active layer is
C   the surface layer

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        IF ( DO_INCLUDE_SURFACE.AND.DO_SPECIALIST_OPTION_2 ) THEN
         IF ( ACTIVE_LAYERS(NLAYERS_TEL).EQ.NLAYERS ) THEN
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              FAC = - DIRECT_BEAM(I,IBEAM,O1)
              DO Q = 1, N_LAYER_WFS
                L_BEAM = ZERO
                DO K = 1, NLAYERS
                  FAC3 = FAC * DELTAU_SLANT(N,K,IBEAM)
                  L_BEAM = L_BEAM + L_DELTAU_VERT(Q,K) * FAC3
                ENDDO
                SCOL2_WF(CM,Q) = SCOL2_WF(CM,Q) + L_BEAM
              ENDDO
            ENDDO
          ENDDO
         ENDIF
        ENDIF
      ENDIF

C  finish

      RETURN
      END

C

      SUBROUTINE LS_BVP_LAMBERTIAN_SURFACE
     I   (  DO_REGULAR_BVP, DO_INCLUDE_DIRECTBEAM,
     I      DO_INCLUDE_SURFEMISS,
     I      IBEAM_INDEX, BRDF_INDEX )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and reflectance variables

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'

C  include file of linearized solution variables (output stored here)

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'

C  inputs
C  ------

C  BVP control input

      LOGICAL          DO_REGULAR_BVP

C  .. direct beam inputs

      LOGICAL          DO_INCLUDE_DIRECTBEAM

C  .. surface emission input

      LOGICAL          DO_INCLUDE_SURFEMISS

C  Beam index

      INTEGER          IBEAM_INDEX

C  BRDF index (not used here)

      INTEGER          BRDF_INDEX

C  local variables
C  ---------------

      INTEGER          N, I, C0, CM, IB, IROW, IR, O1, NS, BR
      INTEGER          K, KO1, K0, K1, K2
      DOUBLE PRECISION L_BEAM, L_HOM_R, L_HOM_CR, H1R, H1I, EMISS_VAR

C  Ground level boundary condition
C  -------------------------------

C  Initialise

      IB = IBEAM_INDEX
      BR = BRDF_INDEX

      IF ( DO_REGULAR_BVP ) THEN
        N  = NLAYERS
        C0 = (N-1)*NSTKS_NSTRMS_2 + NSTKS_NSTRMS
      ELSE
        NS = NLAYERS_TEL
        N  = ACTIVE_LAYERS(NS)      
        C0 = (NS-1)*NSTKS_NSTRMS_2 + NSTKS_NSTRMS
      ENDIF

C  Zero the entries down to the last layer

      DO I = 1, C0
        COL2_WFALB(I,1) = ZERO
      ENDDO

C  Nothing to do if N is not the bottom layer
C    ------- IS THIS CORRECT ????????????????????????????????????

      IF ( N.NE.NLAYERS ) THEN
       DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        DO O1 = 1, NSTOKES
          IROW = IR + O1
          CM   = C0 + IROW
          COL2_WFALB(CM,1) = ZERO
        ENDDO
       ENDDO
       GO TO 457
      ENDIF

C  If this is the surface and Specialist option #2 is in place

      if ( N.EQ.NLAYERS ) THEN

C  Diffuse scatter contributions

        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM   = C0 + IROW

C  Beam contribution for this Kernel

            L_BEAM = R2_BEAM(I,O1)

C  Real homogeneous solutions for this Kernel

            L_HOM_R = ZERO
            DO K = 1, K_REAL(N)
              L_HOM_R = L_HOM_R
     &          + LCON(K,N) * R2_HOMP(I,O1,K) * T_DELT_EIGEN(K,N) 
     &          + MCON(K,N) * R2_HOMM(I,O1,K)
            ENDDO

C  Complex homogeneous solutions for this Kernel

            L_HOM_CR = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K-2
              K1 = KO1 + K0
              K2 = K1 + 1
              H1R =   R2_HOMP(I,O1,K1) * T_DELT_EIGEN(K1,N)
     &              - R2_HOMP(I,O1,K2) * T_DELT_EIGEN(K2,N)
              H1I =   R2_HOMP(I,O1,K1) * T_DELT_EIGEN(K2,N)
     &              + R2_HOMP(I,O1,K2) * T_DELT_EIGEN(K1,N)
              L_HOM_CR = L_HOM_CR
     &                  + LCON(K1,N) * H1R  - LCON(K2,N) * H1I 
     &                  + MCON(K1,N) * R2_HOMM(I,O1,K1)
     &                  - MCON(K2,N) * R2_HOMM(I,O1,K2)
            ENDDO

C  Final contribution

            COL2_WFALB(CM,1) = L_BEAM + L_HOM_R + L_HOM_CR

C  End Streams and Stokes loops

          ENDDO
        ENDDO

C  Add direct beam variation of albedo

        IF ( DO_INCLUDE_DIRECTBEAM ) THEN
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM   = C0 + IROW
              COL2_WFALB(CM,1) = 
     &         COL2_WFALB(CM,1) + DIRECT_BEAM(I,IB,O1)
            ENDDO
          ENDDO
        ENDIF

C  If surface emission, include emissivity variation
C    This code added for Version 2.4RT

        IF ( DO_INCLUDE_SURFEMISS ) THEN
          O1   = 1
          EMISS_VAR = SURFBB
          IF ( DO_SOLAR_SOURCES ) EMISS_VAR = SURFBB * PI4
          EMISS_VAR = LAMBERTIAN_ALBEDO * EMISS_VAR
          DO I = 1, NSTREAMS
            IR   = NSTOKES*(I-1)
            IROW = IR + O1
            CM   = C0 + IROW
            COL2_WFALB(CM,1) = COL2_WFALB(CM,1) - EMISS_VAR
          ENDDO
        ENDIF

C  End clause

      ENDIF

C  Continuation point 

 457  continue

C  Copy for the single layer case

      IF ( NLAYERS .EQ. 1 ) THEN
         DO N = 1, NTOTAL
          SCOL2_WFALB(N,1) = COL2_WFALB(N,1)
        ENDDO
      ENDIF

C  debug

c      if ( do_debug_write ) then
c        DO N = 1, NTOTAL
c          write(85,'(2i4,1p4e17.9)')IBEAM_INDEX, N, COL2_WFALB(N,1)
c        ENDDO
c        pause
c      ENDIF
           
C  Finish

      RETURN
      END
