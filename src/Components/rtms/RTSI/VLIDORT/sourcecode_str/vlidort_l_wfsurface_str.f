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

C ##########################################################
C #                                                        #
C # Subroutines in this Module                             #
C #                                                        #
C #     Top level routines--------------                   #
C #            VLIDORT_ALBEDO_WFS  (master)                #
C #                                                        #
C #     BOA surface source terms ---------                 #
C #            LS_BOA_LAMBERTIAN_SOURCE                    #
C #                                                        #
C #     Recursion relations ---------                      #
C #            UPUSER_SURFACEWF                            #
C #            DNUSER_SURFACEWF                            #
C #                                                        #
C #     Post-processing at user angles --------            #
C #            LS_WHOLELAYER_STERM_UP                      #
C #            LS_WHOLELAYER_STERM_DN                      #
C #            LS_PARTLAYER_STERM_UP                       #
C #            LS_PARTLAYER_STERM_DN                       #
C #                                                        #
C ##########################################################

      SUBROUTINE VLIDORT_ALBEDO_WFS
     I       (  DO_INCLUDE_DIRECTBEAM,
     I          DO_INCLUDE_SURFEMISS,
     I          DO_INCLUDE_MVOUTPUT,
     I          DO_REGULAR_BVP,
     I          FOURIER_COMPONENT, IBEAM,
     I          SURFACE_FACTOR,
     I          FLUX_MULTIPLIER,
     I          BRDF_INDEX, WF_INDEX,
     I          LS_BVP_SURFACE_SETUP,
     I          LS_BOA_SURFACE_SOURCE,
     O          STATUS )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include file of solution variables

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'

C  include file of Linearized solution variables

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'

C  Module arguments
C  ----------------

C  local control flags

      LOGICAL          DO_INCLUDE_DIRECTBEAM
      LOGICAL          DO_INCLUDE_SURFEMISS
      LOGICAL          DO_INCLUDE_MVOUTPUT

C  Flag for regular BVP calculation

      LOGICAL          DO_REGULAR_BVP

C  Fourier surface factor, Fourier number, Beam index, Flux multiplier

      DOUBLE PRECISION SURFACE_FACTOR
      INTEGER          FOURIER_COMPONENT, IBEAM
      DOUBLE PRECISION FLUX_MULTIPLIER

C  Weighting function index number

      INTEGER          WF_INDEX

C  BRDF index (indicates which BRDF kernel is active)

      INTEGER          BRDF_INDEX

C  External Functions

      EXTERNAL         LS_BVP_SURFACE_SETUP
      EXTERNAL         LS_BOA_SURFACE_SOURCE

C  Output status

      INTEGER          STATUS

C  Local variables
C  ---------------

C  error tracing variables

      INTEGER          INFO
      CHARACTER*3      CI
      CHARACTER*70     MAIL, TRACE

C  Linearized BOA terms

      DOUBLE PRECISION LS_BOA_SOURCE(MAX_USER_STREAMS, MAXSTOKES)
      DOUBLE PRECISION LS_BOA_THTONLY_SOURCE(MAXSTREAMS)

C  Other local variables

      INTEGER          I, I1, N, N1, NAC, O1, IC, ICOW
      INTEGER          K, KO1, K0, K1, K2, C0, NS
      INTEGER          IR, IROW, IROW1, IROW_S, IROW1_S        
      DOUBLE PRECISION L_HOM1, L_HOM2, SHOM_R
      DOUBLE PRECISION SHOM_CR, L_HOM1CR, L_HOM2CR
      DOUBLE PRECISION NXR, PXR, NXR1, PXR1, NXR2, PXR2

C  Initialise status

       STATUS = VLIDORT_SUCCESS

C  Go to telescoping solution if flagged
C  -------------------------------------

       IF ( .NOT. DO_REGULAR_BVP ) GO TO 987

C  Regular BVP Solution
C  ====================

C  BV solution for perturbed integration constants
C  -----------------------------------------------

C  Compute the main column B' where AX = B'

      CALL LS_BVP_SURFACE_SETUP
     I       (  DO_REGULAR_BVP,
     I          DO_INCLUDE_DIRECTBEAM,
     I          DO_INCLUDE_SURFEMISS,
     I          IBEAM, BRDF_INDEX )

C  BVP back-substitution: With compression (multilayers)
C  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

C  LAPACK substitution (DGBTRS) using RHS column vector COL2_WF
C  BV solution for perturbed integration constants
C    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS
     &     ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, 1,
     &        BANDMAT2, MAXBANDTOTAL, IPIVOT,
     &        COL2_WFALB, MAXTOTAL, INFO )

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MAIL  = 'argument i illegal value, for i = '//CI
          TRACE = 'DGBTRS call (Reg. multilayer) in VLIDORT_ALBEDO WFS'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ENDIF

C  Set Linearized integration constants NCON_ALB and PCON_ALB, all layers

        DO N = 1, NLAYERS
          C0 = (N-1)*NSTKS_NSTRMS_2
          DO K = 1, K_REAL(N)
            IROW = K
            IROW1 = IROW + NSTKS_NSTRMS
            NCON_ALB(K,N) = COL2_WFALB(C0+IROW,1)
            PCON_ALB(K,N) = COL2_WFALB(C0+IROW1,1)
          ENDDO
          KO1 = K_REAL(N) + 1
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            IROW    = K + K_REAL(N)
            IROW1   = IROW + NSTKS_NSTRMS
            IROW_S  = IROW + K_COMPLEX(N)
            IROW1_S = IROW_S + NSTKS_NSTRMS
            NCON_ALB(K1,N) = COL2_WFALB(C0+IROW,   1)
            NCON_ALB(K2,N) = COL2_WFALB(C0+IROW_S, 1)
            PCON_ALB(K1,N) = COL2_WFALB(C0+IROW1,  1)
            PCON_ALB(K2,N) = COL2_WFALB(C0+IROW1_S,1)
          ENDDO
        ENDDO

C  Solve the boundary problem: No compression, Single Layer only
C  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

C  LAPACK substitution (DGETRS) using RHS column vector SCOL2_WFALB

        CALL DGETRS
     &     ( 'N', NTOTAL, 1, SMAT2, MAXSTRMSTKS_2, SIPIVOT,
     &        SCOL2_WFALB, MAXSTRMSTKS_2, INFO )

C  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MAIL  = 'argument i illegal value, for i = '//CI
          TRACE = 'DGBTRS call (Reg. 1 layer) in VLIDORT_ALBEDO_WFS'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ENDIF

C  Set Linearized integration constants NCON_ALB and PCON_ALB, 1 layer

        N = 1
        DO K = 1, K_REAL(N)
          IROW = K
          IROW1 = IROW + NSTKS_NSTRMS
          NCON_ALB(K,N) = SCOL2_WFALB(IROW,1)
          PCON_ALB(K,N) = SCOL2_WFALB(IROW1,1)
        ENDDO
        KO1 = K_REAL(N) + 1
        DO K = 1, K_COMPLEX(N)
          K0 = 2 * K - 2
          K1 = KO1 + K0
          K2 = K1  + 1
          IROW    = K + K_REAL(N)
          IROW1   = IROW + NSTKS_NSTRMS
          IROW_S  = IROW + K_COMPLEX(N)
          IROW1_S = IROW_S + NSTKS_NSTRMS
          NCON_ALB(K1,N) = SCOL2_WFALB(IROW,   1)
          NCON_ALB(K2,N) = SCOL2_WFALB(IROW_S, 1)
          PCON_ALB(K1,N) = SCOL2_WFALB(IROW1,  1)
          PCON_ALB(K2,N) = SCOL2_WFALB(IROW1_S,1)
        ENDDO

C  end clause
 
      ENDIF

C  debug------------------------------------------
        if ( do_debug_write.and.fourier_component.eq.0 ) then
         DO N = 1, NLAYERS
          DO K = 1, K_REAL(N)
           write(86,'(3i2,1p6e13.5)')FOURIER_COMPONENT,N,K,
     &                LCON(K,N), MCON(K,N),
     &                NCON_ALB(K,N),PCON_ALB(K,N),
     &                NCON_ALB(K,N),PCON_ALB(K,N)
          ENDDO
         ENDDO
        ENDIF
c         pause
C  Resume calculation for post processing

        GO TO 6789

C  Continuation point for telescoped solution

 987    continue

C  Telescoped BVP solution
C  =======================

C  BV solution for perturbed integration constants
C  -----------------------------------------------

C  Compute the main column B' where AX = B'

      CALL LS_BVP_SURFACE_SETUP
     I       (  DO_REGULAR_BVP,
     I          DO_INCLUDE_DIRECTBEAM,
     I          DO_INCLUDE_SURFEMISS,
     I          SURFACE_FACTOR,
     I          IBEAM, BRDF_INDEX )

C  BVP back-substitution: With compression (multilayers)
C  -----------------------------------------------------

      IF ( NLAYERS_TEL .GT. 1 ) THEN

C  LAPACK substitution (DGBTRS) using RHS column vector COL2_WF
C  BV solution for perturbed integration constants
C    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS
     &     ( 'n', N_BVTELMATRIX_SIZE, N_BVTELMATRIX_SUBDIAG, 
     &        N_BVTELMATRIX_SUPDIAG, 1,
     &        BANDTELMAT2, MAXBANDTOTAL, IPIVOTTEL,
     &        COL2_WFALB, MAXTOTAL, INFO )

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MAIL  = 'argument i illegal value, for i = '//CI
          TRACE = 'DGBTRS call (Tel. multilayer) in VLIDORT_ALBEDO WFS'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ENDIF

C  Set linearized integration constants, active layers

        C0 = - NSTKS_NSTRMS_2
        DO NS = 1, NLAYERS_TEL
          N = ACTIVE_LAYERS(NS)
          C0 = C0 + NSTKS_NSTRMS_2
          DO K = 1, K_REAL(N)
            IROW  = K
            IROW1 = IROW + NSTKS_NSTRMS
            NCON_ALB(K,N) = COL2_WFALB(C0+IROW,1)
            PCON_ALB(K,N) = COL2_WFALB(C0+IROW1,1)
          ENDDO
          KO1 = K_REAL(N) + 1
          DO K = 1, K_COMPLEX(N)
            K0 = 2*K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            IROW    = K + K_REAL(N)
            IROW1   = IROW + NSTKS_NSTRMS
            IROW_S  = K + K_REAL(N) + K_COMPLEX(N)
            IROW1_S = IROW_S + NSTKS_NSTRMS
            NCON_ALB(K1,N) = COL2_WFALB(C0+IROW,1)
            NCON_ALB(K2,N) = COL2_WFALB(C0+IROW_S,1)
            PCON_ALB(K1,N) = COL2_WFALB(C0+IROW1,1)
            PCON_ALB(K2,N) = COL2_WFALB(C0+IROW1_S,1)
          ENDDO
        ENDDO

C  Solve the boundary problem: No compression, Single Layer only
C  -------------------------------------------------------------

      ELSE IF ( NLAYERS_TEL .EQ. 1 ) THEN

        NAC = ACTIVE_LAYERS(1)

C  LAPACK substitution (DGETRS) using RHS column vector SCOL2_WFALB

        CALL DGETRS
     &     ( 'N', NSTKS_NSTRMS_2, 1, 
     &        SMAT2, MAXSTRMSTKS_2, SIPIVOT,
     &        SCOL2_WFALB, MAXSTRMSTKS_2, INFO )

C  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MAIL  = 'argument i illegal value, for i = '//CI
          TRACE = 'DGBTRS call (Tel. 1 layer) in VLIDORT_ALBEDO_WFS'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ENDIF

C  Set Linearized integration constants NCON_ALB and PCON_ALB, 1 layer

        DO K = 1, K_REAL(NAC)
          IROW = K
          IROW1 = IROW + NSTKS_NSTRMS
          NCON_ALB(K,NAC) = SCOL2_WFALB(IROW,1)
          PCON_ALB(K,NAC) = SCOL2_WFALB(IROW1,1)
        ENDDO
        KO1 = K_REAL(NAC) + 1
        DO K = 1, K_COMPLEX(NAC)
          K0 = 2 * K - 2
          K1 = KO1 + K0
          K2 = K1  + 1
          IROW    = K + K_REAL(NAC)
          IROW1   = IROW + NSTKS_NSTRMS
          IROW_S  = IROW + K_COMPLEX(NAC)
          IROW1_S = IROW_S + NSTKS_NSTRMS
          NCON_ALB(K1,NAC) = SCOL2_WFALB(IROW,   1)
          NCON_ALB(K2,NAC) = SCOL2_WFALB(IROW_S, 1)
          PCON_ALB(K1,NAC) = SCOL2_WFALB(IROW1,  1)
          PCON_ALB(K2,NAC) = SCOL2_WFALB(IROW1_S,1)
        ENDDO

C  end clause
 
      ENDIF

C  Set linearized integration constants for non-active layers
C  ==========================================================

C  Now we propagate the results upwards and downwards through the
C  appropriate non-active layers where there is no scattering.

C  Transmittance layers ABOVE active layer(s)
C  -----------------------------------------

C   --NCON values are zero (no downwelling radiation)
C   --PCON values propagated upwards from top of first active layer

C  layer immediately above first active layer
C   --- Require linearized solutions at top of first active layer

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

C  real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(NAC)
              NXR = NCON_ALB(K,NAC)*SOLA_XPOS(I1,O1,K,NAC)
              PXR = PCON_ALB(K,NAC)*SOLB_XNEG(I1,O1,K,NAC)
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
              NXR1  =   NCON_ALB(K1,NAC) *   SOLA_XPOS(I1,O1,K1,NAC)
     &                - NCON_ALB(K2,NAC) *   SOLA_XPOS(I1,O1,K2,NAC)
              PXR1  =   PCON_ALB(K1,NAC) *   SOLB_XNEG(I1,O1,K1,NAC)
     &                - PCON_ALB(K2,NAC) *   SOLB_XNEG(I1,O1,K2,NAC)
              PXR2  =   PCON_ALB(K1,NAC) *   SOLB_XNEG(I1,O1,K2,NAC)
     &                + PCON_ALB(K2,NAC) *   SOLB_XNEG(I1,O1,K1,NAC)
              L_HOM1CR = NXR1
              L_HOM2CR =  + T_DELT_EIGEN(K1,NAC) * PXR1
     &                    - T_DELT_EIGEN(K2,NAC) * PXR2
              SHOM_CR = SHOM_CR + L_HOM1CR + L_HOM2CR
            ENDDO

C    ---Sets Real integration constants (no complex ones)

            PCON_ALB(ICOW,N1) = SHOM_R + SHOM_CR
            NCON_ALB(ICOW,N1) = ZERO
 
C  End loops

          ENDDO
        ENDDO

C  End first active layer above

      ENDIF

C  For remaining non-active atmospheric layers to TOA, propagate upwards.

      DO N = NAC - 2, 1, -1
        N1 = N + 1
        DO I = 1, NSTREAMS
          IR = ( I - 1 ) * NSTOKES
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            NCON_ALB(IROW,N) = ZERO
            PCON_ALB(IROW,N) =
     &            T_DELT_DISORDS(I,N1) * PCON_ALB(IROW,N1)
          ENDDO
        ENDDO
      ENDDO

C  Transmittance layers below active layer(s)
C  -----------------------------------------

C   -- PCON values are zero (no upwelling radiation)
C   -- NCON values propagated downwards from bottom of last active layer

C  layer immediately below Last active layer
C    .... Require linearized solutions at bottom of last active layer

      NAC = ACTIVE_LAYERS (NLAYERS_TEL)
      IF ( NAC .LT. NLAYERS ) THEN

        N1 = NAC + 1

C  start stream, stokes loops

        DO I = 1, NSTREAMS
          IR = ( I - 1 ) * NSTOKES
          DO O1 = 1, NSTOKES
            IROW = IR + O1
 
C  real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(NAC)
              NXR = NCON_ALB(K,NAC) * SOLA_XPOS(I,O1,K,NAC)
              PXR = PCON_ALB(K,NAC) * SOLB_XNEG(I,O1,K,NAC)
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
              NXR1  =   NCON_ALB(K1,NAC) *   SOLA_XPOS(I,O1,K1,NAC)
     &                - NCON_ALB(K2,NAC) *   SOLA_XPOS(I,O1,K2,NAC)
              NXR2  =   NCON_ALB(K1,NAC) *   SOLA_XPOS(I,O1,K2,NAC)
     &                + NCON_ALB(K2,NAC) *   SOLA_XPOS(I,O1,K1,NAC)
              PXR1  =   PCON_ALB(K1,NAC) *   SOLB_XNEG(I,O1,K1,NAC)
     &                - PCON_ALB(K2,NAC) *   SOLB_XNEG(I,O1,K2,NAC)
              L_HOM2CR = PXR1
              L_HOM1CR =  + T_DELT_EIGEN(K1,NAC) * NXR1
     &                    - T_DELT_EIGEN(K2,NAC) * NXR2
              SHOM_CR = SHOM_CR + L_HOM1CR + L_HOM2CR
            ENDDO

C    ---Sets Real integration constants (no complex ones)

            NCON_ALB(IROW,N1) = SHOM_R + SHOM_CR
            PCON_ALB(IROW,N1) = ZERO

C  End loops

          ENDDO
        ENDDO

C  End first active layer beneath

      ENDIF

C  other layers to bottom of medium: propagate downwards.
C   Additional variation if you are passing through the varying layer.

      DO N = NAC + 2, NLAYERS
        N1 = N - 1
        DO I = 1, NSTREAMS
          IR = ( I - 1 ) * NSTOKES
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            PCON_ALB(IROW,N) = ZERO
            NCON_ALB(IROW,N) =
     &            T_DELT_DISORDS(I,N1) * NCON_ALB(IROW,N1)
          ENDDO
        ENDDO
      ENDDO

C  Get the Post-processed weighting functions
C  ==========================================

C Continuation point

 6789 continue

C  Upwelling weighting functions
C  -----------------------------

      IF ( DO_UPWELLING ) THEN

C  Get the surface term (L_BOA_SOURCE). External Function

        CALL LS_BOA_SURFACE_SOURCE
     I        ( DO_INCLUDE_DIRECTBEAM,
     I          DO_INCLUDE_SURFEMISS,
     I          DO_INCLUDE_MVOUTPUT,
     I          SURFACE_FACTOR, 
     I          FOURIER_COMPONENT, IBEAM, 
     O          LS_BOA_SOURCE,
     O          LS_BOA_THTONLY_SOURCE )

C  Upwelling Albedo weighting functions

        CALL UPUSER_SURFACEWF
     I      ( FLUX_MULTIPLIER, IBEAM, WF_INDEX,
     I        LS_BOA_SOURCE )

      ENDIF

C  Downwelling Albedo weighting functions
C  --------------------------------------

      IF ( DO_DNWELLING ) THEN
        CALL DNUSER_SURFACEWF
     I       ( FLUX_MULTIPLIER, IBEAM,
     I         WF_INDEX )
      ENDIF

C  mean value output
C  -----------------

      IF ( DO_INCLUDE_MVOUTPUT.OR.DO_QUAD_OUTPUT ) THEN
        CALL VLIDORT_LS_INTEGRATED_OUTPUT
     I  ( DO_INCLUDE_MVOUTPUT, 
     I    FLUX_MULTIPLIER, IBEAM, WF_INDEX,
     I    LS_BOA_THTONLY_SOURCE )
      ENDIF

C  Finish

      RETURN
      END

c

      SUBROUTINE LS_BOA_LAMBERTIAN_SOURCE
     I        ( DO_INCLUDE_DIRECTBEAM,
     I          DO_INCLUDE_SURFEMISS,
     I          DO_INCLUDE_MVOUTPUT,
     I          SURFACE_FACTOR, 
     I          FOURIER_COMPONENT, 
     I          IBEAM, 
     O          LS_BOA_SOURCE,
     O          LS_BOA_THTONLY_SOURCE )

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

C  include file of Linearized solution variables

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_THERMALSUP.VARS'

C  Module arguments
C  ----------------

C  directbeam inclusion flag

      LOGICAL          DO_INCLUDE_DIRECTBEAM

C  surface emissivity inclusion flag

      LOGICAL          DO_INCLUDE_SURFEMISS

C  MV output inclusion flag

      LOGICAL          DO_INCLUDE_MVOUTPUT

C  Fourier surface factor

      DOUBLE PRECISION SURFACE_FACTOR

C  Fourier number

      INTEGER          FOURIER_COMPONENT

C  Beam index

      INTEGER          IBEAM

C  output

      DOUBLE PRECISION LS_BOA_SOURCE(MAX_USER_STREAMS, MAXSTOKES)
      DOUBLE PRECISION LS_BOA_THTONLY_SOURCE(MAXSTREAMS)

C  Local variables
C  ---------------

      LOGICAL          DO_QTHTONLY
      INTEGER          UM, I, J, N, IB, O1
      INTEGER          K, KO1, K0, K1, K2      
      DOUBLE PRECISION INTEGRAND(MAXSTREAMS,MAXSTOKES)
      DOUBLE PRECISION SUM_R, SUM_CR, REFLEC, FACTOR, EMISSVAR
      DOUBLE PRECISION H1, H2, NXR, PXR, NXR1, NXR2, PXR1

C  Initialise
C  ----------

C  Special flag

      DO_QTHTONLY = ( DO_THERMAL_TRANSONLY ) .AND.
     &      ( DO_QUAD_OUTPUT .OR. DO_INCLUDE_MVOUTPUT )

C  Shorthand

      N   = NLAYERS
      KO1 = K_REAL(N) + 1
      IB  = IBEAM

C  initialise Derivative of BOA source function

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO O1 = 1, NSTOKES
          LS_BOA_SOURCE(UM,O1) = ZERO
        ENDDO
      ENDDO

C  Thermal tranmsittance only, special term

      IF ( DO_QTHTONLY ) THEN
        DO I = 1, NSTREAMS
          LS_BOA_THTONLY_SOURCE(I) = ZERO
        ENDDO
      ENDIF

C  Skip diffuse-field variation for thermal transmittance-only

      IF ( DO_THERMAL_TRANSONLY .OR. .NOT. DO_USER_STREAMS ) GO TO 599

C  Contribution due to derivatives of BV constants
C  -----------------------------------------------

C  First compute derivative of downward intensity Integrand at stream angles 
C        .. reflectance integrand  = a(j).x(j).I_DOWN(-j)

C  start loops

      DO I = 1, NSTREAMS
       DO O1 = 1, NSTOKES

C  Real homogeneous solutions

        SUM_R = ZERO
        DO K = 1, K_REAL(N)
         NXR = NCON_ALB(K,N) * SOLA_XPOS(I,O1,K,N)
         PXR = PCON_ALB(K,N) * SOLB_XNEG(I,O1,K,N)
         SUM_R = SUM_R + NXR*T_DELT_EIGEN(K,N) + PXR
         ENDDO

C  Complex solutions

        SUM_CR = ZERO
        DO K = 1, K_COMPLEX(N)
         K0 = 2 * K - 2
         K1 = KO1 + K0
         K2 = K1  + 1
         NXR1 =   NCON_ALB(K1,N) * SOLA_XPOS(I,O1,K1,N)
     &          - NCON_ALB(K2,N) * SOLA_XPOS(I,O1,K2,N)
         NXR2 =   NCON_ALB(K1,N) * SOLA_XPOS(I,O1,K2,N)
     &          + NCON_ALB(K2,N) * SOLA_XPOS(I,O1,K1,N)
         PXR1 =   PCON_ALB(K1,N) * SOLB_XNEG(I,O1,K1,N)
     &          - PCON_ALB(K2,N) * SOLB_XNEG(I,O1,K2,N)
         H1 =  NXR1 * T_DELT_EIGEN(K1,N)
     &        -NXR2 * T_DELT_EIGEN(K2,N)
         H2 =  PXR1
         SUM_CR = SUM_CR + H1 + H2
        ENDDO

C  Final result

        INTEGRAND(I,O1) = QUAD_STRMWTS(I) * ( SUM_R + SUM_CR )

C  end loops

       ENDDO
      ENDDO

C  integrated reflectance term
C  ---------------------------

C  integrate Lambertian case, same for all user-streams

      IF ( FOURIER_COMPONENT.EQ.0 ) THEN
        O1 = 1
        REFLEC = ZERO
        DO J = 1, NSTREAMS
           REFLEC = REFLEC + INTEGRAND(J,O1)
        ENDDO
        FACTOR = SURFACE_FACTOR * LAMBERTIAN_ALBEDO
        REFLEC = FACTOR * REFLEC
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          LS_BOA_SOURCE(UM,O1) = REFLEC
        ENDDO
      ENDIF

C  Continuation point for avoiding diffuse field computation

 599  continue

C  Contributions due to direct variation of kernel factor
C  ------------------------------------------------------

C  Add linearization term due to kernel factor variation 

      IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
        EMISSVAR = SURFBB
        IF ( DO_SOLAR_SOURCES ) EMISSVAR = PI4 * SURFBB
        o1 = 1
        EMISSVAR = EMISSVAR * ( one - lambertian_albedo )
        DO UM = LOCAL_UM_START, N_USER_STREAMS
            LS_BOA_SOURCE(UM,O1) = LS_BOA_SOURCE(UM,O1)
     &                + BOA_DIFFUSE_SOURCE (UM,O1) - EMISSVAR
        ENDDO
      ENDIF

C  Thermal tranmsittance only (additional quadrature terms if flagged)

      IF ( DO_QTHTONLY ) THEN
        DO I = 1, NSTREAMS
           LS_BOA_THTONLY_SOURCE(I) = BOA_THTONLY_SOURCE(I)
        ENDDO
      ENDIF

C  Add linearization term due to kernel factor variation of Direct Beam
C    This term only applies if there is no DB_CORRECTION

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        IF ( FOURIER_COMPONENT.EQ.0.AND..NOT.DO_DBCORRECTION ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              LS_BOA_SOURCE(UM,O1) = LS_BOA_SOURCE(UM,O1)        
     &              + USER_DIRECT_BEAM ( UM, IB, O1 )
            ENDDO
          ENDDO
        ENDIF
      ENDIF

C  Add emissivity variation at user defined angles
C  (expression for emissivity variation follows from Kirchhoff's law)

      IF ( DO_INCLUDE_SURFEMISS ) THEN
        O1 = 1
        EMISSVAR = SURFBB
        IF ( DO_SOLAR_SOURCES ) EMISSVAR = PI4 * SURFBB
        EMISSVAR = EMISSVAR * LAMBERTIAN_ALBEDO
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          LS_BOA_SOURCE(UM,O1) = LS_BOA_SOURCE(UM,O1) - EMISSVAR
        ENDDO
        IF ( DO_QTHTONLY ) THEN
          DO I = 1, NSTREAMS
            LS_BOA_THTONLY_SOURCE(I) =
     &           LS_BOA_THTONLY_SOURCE(I) - EMISSVAR
          ENDDO
        ENDIF
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE UPUSER_SURFACEWF
     I   ( FLUX_MULTIPLIER, IBEAM, WF_INDEX,
     I     LS_BOA_SOURCE )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and reflectance variables (input)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'

C  include files of linearized setup and solution variables (input)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'

C  include files of result variables (module output stored here)

      INCLUDE '../includes/VLIDORT_L_RESULTS.VARS'

C  Subroutine input arguments
C  --------------------------

C  Flux multiplier

      DOUBLE PRECISION FLUX_MULTIPLIER

C  Beam index

      INTEGER          IBEAM

C  weighting function index

      INTEGER          WF_INDEX

C  derivatives of reflected surface upwelling intensity

      DOUBLE PRECISION LS_BOA_SOURCE(MAX_USER_STREAMS, MAXSTOKES)

C  local variables
C  ---------------

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL, O1
      INTEGER          UTA, UM, QW, UT, IB

      DOUBLE PRECISION LS_CUMUL_SOURCE(MAX_USER_STREAMS,MAXSTOKES)
      DOUBLE PRECISION LS_LAYER_SOURCE(MAX_USER_STREAMS,MAXSTOKES)
      DOUBLE PRECISION LS_FINAL_SOURCE

C  index

      QW = WF_INDEX
      IB = IBEAM

C  Zero all Fourier components - New rule, better for safety
C    Only did this for components close to zenith (formerly)

      IF ( DO_USER_STREAMS ) THEN
        DO UTA = 1, N_USER_LEVELS
          DO UM = 1, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              SURFACEWF_F(QW,UTA,UM,IB,O1,UPIDX) = ZERO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

C  Initialize post-processing recursion
C  ====================================

      IF ( DO_USER_STREAMS ) THEN

C  Set the cumuluative source term equal to the BOA sum

        DO UM = LOCAL_UM_START, N_USER_STREAMS
         DO O1 = 1, NSTOKES
           LS_CUMUL_SOURCE(UM,O1) = LS_BOA_SOURCE(UM,O1)
         ENDDO
        ENDDO

      ENDIF

C  Recursion Loop for linearized Post-processing
C  =============================================

C  initialise cumulative source term loop

      NUT = 0
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

C  loop over all output optical depths
C  -----------------------------------

      DO UTA = N_USER_LEVELS, 1, -1

C  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_UP(UTA)

C  Cumulative source terms to layer NUT (user-defined stream angles only)
C    1. Get layer source terms
C    2. Find cumulative source term
C    3. Set multiple scatter source term (MSST) output if flagged

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL + 1
          DO N = NSTART, NUT, -1
            CALL LS_WHOLELAYER_STERM_UP
     I       ( IB, N, LS_LAYER_SOURCE )
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                LS_CUMUL_SOURCE(UM,O1) = LS_LAYER_SOURCE(UM,O1)
     &               + T_DELT_USERM(N,UM) * LS_CUMUL_SOURCE(UM,O1)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

C  Offgrid output
C  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)

C  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN
            CALL LS_PARTLAYER_STERM_UP
     I        ( IB, UT, N, LS_LAYER_SOURCE )
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                LS_FINAL_SOURCE = LS_LAYER_SOURCE(UM,O1)
     &                + T_UTUP_USERM(UT,UM) * LS_CUMUL_SOURCE(UM,O1)
                SURFACEWF_F(QW,UTA,UM,IB,O1,UPIDX) =
     &                   FLUX_MULTIPLIER * LS_FINAL_SOURCE
              ENDDO
            ENDDO
          ENDIF

C  Ongrid output
C  -------------

        ELSE

C  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                SURFACEWF_F(QW,UTA,UM,IB,O1,UPIDX) =
     &                     FLUX_MULTIPLIER * LS_CUMUL_SOURCE(UM,O1)
              ENDDO
            ENDDO
          ENDIF

        ENDIF

C  Check for updating the recursion

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT. NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT
        ENDIF

C  end loop over optical depth

      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE DNUSER_SURFACEWF
     I    ( FLUX_MULTIPLIER,
     I      IBEAM, WF_INDEX )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and reflectance variables (input)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'

C  include files of linearized setup and solution variables (input)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'

C  include file of result variables (module output stored here)

      INCLUDE '../includes/VLIDORT_L_RESULTS.VARS'

C  Subroutine input arguments
C  --------------------------

C  Beam index

      INTEGER          IBEAM

C  weighting function index

      INTEGER          WF_INDEX

C  Flux multiplier = F/4.pi

      DOUBLE PRECISION FLUX_MULTIPLIER

C  local variables
C  ---------------

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL, O1
      INTEGER          UTA, UM, QW, UT, IB

      DOUBLE PRECISION LS_CUMUL_SOURCE ( MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION LS_TOA_SOURCE   ( MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION LS_LAYER_SOURCE ( MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION LS_FINAL_SOURCE

C  Initialise

      QW = WF_INDEX
      IB = IBEAM

C  Zero all Fourier component output

      IF ( DO_USER_STREAMS ) THEN
        DO UTA = 1, N_USER_LEVELS
          DO UM = 1, LOCAL_UM_START 
            DO O1 = 1, NSTOKES
              SURFACEWF_F(QW,UTA,UM,IB,O1,DNIDX) = ZERO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

C  Initialize post-processing recursion
C  ====================================

C  Get the linearized TOA source terms

      IF ( DO_USER_STREAMS ) THEN

        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            LS_TOA_SOURCE(UM,O1) = ZERO
          ENDDO
        ENDDO

        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            LS_CUMUL_SOURCE(UM,O1) = LS_TOA_SOURCE(UM,O1)
          ENDDO
        ENDDO

      ENDIF

C  Recursion Loop for linearized Post-processing
C  =============================================

C  initialise cumulative source term loop

      NUT = 0
      NSTART = 1
      NUT_PREV = NSTART - 1

C  loop over all output optical depths
C  -----------------------------------

      DO UTA = 1, N_USER_LEVELS

C  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_DN(UTA)

C  Cumulative source terms to layer NUT (user-defined stream angles only)
C    1. Get layer source terms
C    2. Find cumulative source term
C    3. Set multiple scatter source term output if flagged

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL
          DO N = NSTART, NUT
            CALL LS_WHOLELAYER_STERM_DN
     &       ( IB, N, LS_LAYER_SOURCE )
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                LS_CUMUL_SOURCE(UM,O1) = LS_LAYER_SOURCE(UM,O1)
     &               + T_DELT_USERM(N,UM) * LS_CUMUL_SOURCE(UM,O1)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

C  Offgrid output
C  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)

C  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN
            CALL LS_PARTLAYER_STERM_DN
     I        ( IB, UT, N, LS_LAYER_SOURCE )
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                LS_FINAL_SOURCE = LS_LAYER_SOURCE(UM,O1)
     &              + T_UTDN_USERM(UT,UM) * LS_CUMUL_SOURCE(UM,O1)
                SURFACEWF_F(QW,UTA,UM,IB,O1,DNIDX) =
     &                FLUX_MULTIPLIER * LS_FINAL_SOURCE
              ENDDO
            ENDDO
          ENDIF

C  Ongrid output
C  -------------

        ELSE

C  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                SURFACEWF_F(QW,UTA,UM,IB,O1,DNIDX) =
     &                    FLUX_MULTIPLIER * LS_CUMUL_SOURCE(UM,O1)
              ENDDO
            ENDDO
          ENDIF

        ENDIF

C  Check for updating the recursion

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT. NE. NUT_PREV ) NSTART = NUT + 1
          NUT_PREV = NUT
        ENDIF

C  end loop over optical depth

      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE LS_WHOLELAYER_STERM_UP
     I       ( IBEAM, GIVEN_LAYER,
     O         LS_LAYERSOURCE )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of ssolution and multiplier variables (input)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  include files of linearized multiplier and solution variables (input)

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'

C  Subroutine input arguments
C  --------------------------

C  Indices

      INTEGER          IBEAM, GIVEN_LAYER

C  Subroutine output arguments
C  ---------------------------

      DOUBLE PRECISION LS_LAYERSOURCE ( MAX_USER_STREAMS, MAXSTOKES )

C  local variables
C  ---------------

      INTEGER          N, UM, O1, IB, K, KO1, K0, K1, K2
      DOUBLE PRECISION SHOM_R, SHOM_CR, H1, H2
      DOUBLE PRECISION NUXR, PUXR, NUXR1, NUXR2, PUXR1, PUXR2

C  Thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) THEN
        DO O1 = 1, NSTOKES
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            LS_LAYERSOURCE(UM,O1) = ZERO
          ENDDO
        ENDDO 
        RETURN
      ENDIF

C  local indices

      N   = GIVEN_LAYER
      KO1 = K_REAL(N) + 1
      IB  = IBEAM

C  Homogeneous solutions
C  =====================

C  Loop over user angles and Stokes

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO O1 = 1, NSTOKES

C  Real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, K_REAL(N)
            NUXR = NCON_ALB(K,N)*UHOM_UPDN(UM,O1,K,N)
            PUXR = PCON_ALB(K,N)*UHOM_UPUP(UM,O1,K,N)
            H1 =  NUXR * HMULT_2(K,UM,N)
            H2 =  PUXR * HMULT_1(K,UM,N)
            SHOM_R = SHOM_R + H1 + H2
          ENDDO

C  Complex homogeneous solutions

          SHOM_CR = ZERO
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            NUXR1 =   NCON_ALB(K1,N) * UHOM_UPDN(UM,O1,K1,N)
     &              - NCON_ALB(K2,N) * UHOM_UPDN(UM,O1,K2,N)
            NUXR2 =   NCON_ALB(K1,N) * UHOM_UPDN(UM,O1,K2,N)
     &              + NCON_ALB(K2,N) * UHOM_UPDN(UM,O1,K1,N)
            PUXR1 =   PCON_ALB(K1,N) * UHOM_UPUP(UM,O1,K1,N)
     &              - PCON_ALB(K2,N) * UHOM_UPUP(UM,O1,K2,N)
            PUXR2 =   PCON_ALB(K1,N) * UHOM_UPUP(UM,O1,K2,N)
     &              + PCON_ALB(K2,N) * UHOM_UPUP(UM,O1,K1,N)
            H1 =   NUXR1 * HMULT_2(K1,UM,N)
     &           - NUXR2 * HMULT_2(K2,UM,N)
            H2 =   PUXR1 * HMULT_1(K1,UM,N)
     &           - PUXR2 * HMULT_1(K2,UM,N)
            SHOM_CR = SHOM_CR + H1 + H2
          ENDDO

C  homogeneous contribution

          LS_LAYERSOURCE(UM,O1) = SHOM_R + SHOM_CR

C  End loops over O1 and UM

        ENDDO
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE LS_WHOLELAYER_STERM_DN
     I       ( IBEAM, GIVEN_LAYER,
     O         LS_LAYERSOURCE )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of ssolution and multiplier variables (input)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  include files of linearized multiplier and solution variables (input)

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'

C  Subroutine input arguments
C  --------------------------

C  Indices

      INTEGER          IBEAM, GIVEN_LAYER

C  Subroutine output arguments
C  ---------------------------

      DOUBLE PRECISION LS_LAYERSOURCE ( MAX_USER_STREAMS, MAXSTOKES )

C  local variables
C  ---------------

      INTEGER          N, UM, O1, IB, K, KO1, K0, K1, K2
      DOUBLE PRECISION SHOM_R, SHOM_CR, H1, H2
      DOUBLE PRECISION NUXR, PUXR, NUXR1, NUXR2, PUXR1, PUXR2

C  Thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) THEN
        DO O1 = 1, NSTOKES
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            LS_LAYERSOURCE(UM,O1) = ZERO
          ENDDO
        ENDDO 
        RETURN
      ENDIF

C  local indices

      N   = GIVEN_LAYER
      KO1 = K_REAL(N) + 1
      IB  = IBEAM

C  Homogeneous solutions
C  =====================

C  Loop over user angles and STokes

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO O1 = 1, NSTOKES

C  Real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, K_REAL(N)
            NUXR = NCON_ALB(K,N)*UHOM_DNDN(UM,O1,K,N)
            PUXR = PCON_ALB(K,N)*UHOM_DNUP(UM,O1,K,N)
            H1 =  NUXR * HMULT_1(K,UM,N)
            H2 =  PUXR * HMULT_2(K,UM,N)
            SHOM_R = SHOM_R + H1 + H2
          ENDDO

C  Complex homogeneous solutions

          SHOM_CR = ZERO
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            NUXR1 =   NCON_ALB(K1,N) * UHOM_DNDN(UM,O1,K1,N)
     &              - NCON_ALB(K2,N) * UHOM_DNDN(UM,O1,K2,N)
            NUXR2 =   NCON_ALB(K1,N) * UHOM_DNDN(UM,O1,K2,N)
     &              + NCON_ALB(K2,N) * UHOM_DNDN(UM,O1,K1,N)
            PUXR1 =   PCON_ALB(K1,N) * UHOM_DNUP(UM,O1,K1,N)
     &              - PCON_ALB(K2,N) * UHOM_DNUP(UM,O1,K2,N)
            PUXR2 =   PCON_ALB(K1,N) * UHOM_DNUP(UM,O1,K2,N)
     &              + PCON_ALB(K2,N) * UHOM_DNUP(UM,O1,K1,N)
            H1 =   NUXR1 * HMULT_1(K1,UM,N)
     &           - NUXR2 * HMULT_1(K2,UM,N)
            H2 =   PUXR1 * HMULT_2(K1,UM,N)
     &           - PUXR2 * HMULT_2(K2,UM,N)
            SHOM_CR = SHOM_CR + H1 + H2
          ENDDO

C  homogeneous contribution

          LS_LAYERSOURCE(UM,O1) = SHOM_R + SHOM_CR

C  End loops over O1 and UM

        ENDDO
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE LS_PARTLAYER_STERM_UP
     I       ( IBEAM, OFFGRID_INDEX, GIVEN_LAYER,
     O         LS_LAYERSOURCE )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of ssolution and multiplier variables (input)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  include files of linearized multiplier and solution variables (input)

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'

C  Subroutine input arguments
C  --------------------------

C  layer and beam index

      INTEGER          GIVEN_LAYER, IBEAM

C  offgrid optical depth index

      INTEGER          OFFGRID_INDEX

C  Subroutine output arguments
C  ---------------------------

      DOUBLE PRECISION LS_LAYERSOURCE ( MAX_USER_STREAMS, MAXSTOKES )

C  local variables
C  ---------------

      INTEGER          N, UM, O1, IB, UT, K, KO1, K0, K1, K2
      DOUBLE PRECISION SHOM_R, SHOM_CR, H1, H2
      DOUBLE PRECISION NUXR, PUXR, NUXR1, NUXR2, PUXR1, PUXR2

C  Thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) THEN
        DO O1 = 1, NSTOKES
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            LS_LAYERSOURCE(UM,O1) = ZERO
          ENDDO
        ENDDO 
        RETURN
      ENDIF

C  local indices

      N   = GIVEN_LAYER
      KO1 = K_REAL(N) + 1
      UT  = OFFGRID_INDEX
      IB  = IBEAM

C  Partial layer source function ( Homogeneous/constants variation )
C  =================================================================

C  Loop over user angles and STokes

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO O1 = 1, NSTOKES

C  Real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, K_REAL(N)
            NUXR = NCON_ALB(K,N) * UHOM_UPDN(UM,O1,K,N)
            PUXR = PCON_ALB(K,N) * UHOM_UPUP(UM,O1,K,N)
            H1 =  NUXR * UT_HMULT_UD(K,UM,UT)
            H2 =  PUXR * UT_HMULT_UU(K,UM,UT)
            SHOM_R = SHOM_R + H1 + H2
          ENDDO

C  Complex homogeneous solutions

          SHOM_CR = ZERO
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            NUXR1 =   NCON_ALB(K1,N) * UHOM_UPDN(UM,O1,K1,N)
     &              - NCON_ALB(K2,N) * UHOM_UPDN(UM,O1,K2,N)
            NUXR2 =   NCON_ALB(K1,N) * UHOM_UPDN(UM,O1,K2,N)
     &              + NCON_ALB(K2,N) * UHOM_UPDN(UM,O1,K1,N)
            PUXR1 =   PCON_ALB(K1,N) * UHOM_UPUP(UM,O1,K1,N)
     &              - PCON_ALB(K2,N) * UHOM_UPUP(UM,O1,K2,N)
            PUXR2 =   PCON_ALB(K1,N) * UHOM_UPUP(UM,O1,K2,N)
     &              + PCON_ALB(K2,N) * UHOM_UPUP(UM,O1,K1,N)
            H1 =   NUXR1 * UT_HMULT_UD(K1,UM,N)
     &           - NUXR2 * UT_HMULT_UD(K2,UM,N)
            H2 =   PUXR1 * UT_HMULT_UU(K1,UM,N)
     &           - PUXR2 * UT_HMULT_UU(K2,UM,N)
            SHOM_CR = SHOM_CR + H1 + H2
          ENDDO

C  homogeneous contribution

          LS_LAYERSOURCE(UM,O1) = SHOM_R + SHOM_CR

C  End loops over O1 and UM

        ENDDO
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE LS_PARTLAYER_STERM_DN
     I       ( IBEAM, OFFGRID_INDEX, GIVEN_LAYER, 
     O         LS_LAYERSOURCE )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of ssolution and multiplier variables (input)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  include files of linearized multiplier and solution variables (input)

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'

C  Subroutine input arguments
C  --------------------------

C  layer and beam indices

      INTEGER          GIVEN_LAYER, IBEAM

C  offgrid optical depth index

      INTEGER          OFFGRID_INDEX

C  Subroutine output arguments
C  ---------------------------

      DOUBLE PRECISION LS_LAYERSOURCE ( MAX_USER_STREAMS, MAXSTOKES )

C  local variables
C  ---------------

      INTEGER          N, UM, O1, IB, UT, K, KO1, K0, K1, K2
      DOUBLE PRECISION SHOM_R, SHOM_CR, H1, H2
      DOUBLE PRECISION NUXR, PUXR, NUXR1, NUXR2, PUXR1, PUXR2

C  Thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) THEN
        DO O1 = 1, NSTOKES
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            LS_LAYERSOURCE(UM,O1) = ZERO
          ENDDO
        ENDDO 
        RETURN
      ENDIF

C  local indices

      N   = GIVEN_LAYER
      KO1 = K_REAL(N) + 1
      UT  = OFFGRID_INDEX
      IB  = IBEAM

C  Partial layer source function ( Homogeneous/constants variation )
C  =================================================================

C  Loop over user angles and STokes

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO O1 = 1, NSTOKES

C  Real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, K_REAL(N)
            NUXR = NCON_ALB(K,N) * UHOM_DNDN(UM,O1,K,N)
            PUXR = PCON_ALB(K,N) * UHOM_DNUP(UM,O1,K,N)
            H1 =  NUXR * UT_HMULT_DD(K,UM,UT)
            H2 =  PUXR * UT_HMULT_DU(K,UM,UT)
            SHOM_R = SHOM_R + H1 + H2
          ENDDO

C  Complex homogeneous solutions

          SHOM_CR = ZERO
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            NUXR1 =   NCON_ALB(K1,N) * UHOM_DNDN(UM,O1,K1,N)
     &              - NCON_ALB(K2,N) * UHOM_DNDN(UM,O1,K2,N)
            NUXR2 =   NCON_ALB(K1,N) * UHOM_DNDN(UM,O1,K2,N)
     &              + NCON_ALB(K2,N) * UHOM_DNDN(UM,O1,K1,N)
            PUXR1 =   PCON_ALB(K1,N) * UHOM_DNUP(UM,O1,K1,N)
     &              - PCON_ALB(K2,N) * UHOM_DNUP(UM,O1,K2,N)
            PUXR2 =   PCON_ALB(K1,N) * UHOM_DNUP(UM,O1,K2,N)
     &              + PCON_ALB(K2,N) * UHOM_DNUP(UM,O1,K1,N)
            H1 =   NUXR1 * UT_HMULT_DD(K1,UM,N)
     &           - NUXR2 * UT_HMULT_DD(K2,UM,N)
            H2 =   PUXR1 * UT_HMULT_DU(K1,UM,N)
     &           - PUXR2 * UT_HMULT_DU(K2,UM,N)
            SHOM_CR = SHOM_CR + H1 + H2
          ENDDO

C  homogeneous contribution

          LS_LAYERSOURCE(UM,O1) = SHOM_R + SHOM_CR

C  End loops over O1 and UM

        ENDDO
      ENDDO

C  Finish

      RETURN
      END

