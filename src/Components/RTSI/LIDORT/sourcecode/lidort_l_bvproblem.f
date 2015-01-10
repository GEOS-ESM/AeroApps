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
C #  For linearizations involving atmospheric parameters        #
C #            L_BVP_SOLUTION_MASTER                            #
C #            L_BVP_P_COLUMN_SETUP (re-named, Version 3.3)     #
C #            L_BVP_C_COLUMN_SETUP (new for version 3.3)       #
C #            L_BVP_SURFACE_SETUP                              #
C #            L_BEAMSOLUTION_P_NEQK (re-named)                 #
C #            L_BEAMSOLUTION_P_NNEK (re-named)                 #
C #            L_BEAMSOLUTION_C_NEQK (new for version 3.3)      #
C #                                                             #
C #  For linearizations involving surface parameters            #
C #            KFACTOR_WF_COLSETUP                              #
C #            KPARAMS_WF_COLSETUP                              #
C #                                                             #
C #  For linearizations (atmospheric) using telescoped BVP      #
C #            L_BVPTEL_SOLUTION_MASTER                         #
C #            L_BVPTEL_COLUMN_SETUP                            #
C #    Placeholder, Version 3.3, modify telecoped problem       #
C #                                                             #
C ###############################################################

      SUBROUTINE L_BVP_SOLUTION_MASTER
     I       ( DO_INCLUDE_SURFACE,
     I         DO_INCLUDE_SURFEMISS,
     I         DO_INCLUDE_THERMEMISS,
     I         DO_INCLUDE_DIRECTBEAM,
     I         VARIATION_INDEX, N_WEIGHTFUNCS,
     I         FOURIER_COMPONENT, IBEAM,
     I         SURFACE_FACTOR,
     O         STATUS )

C  include files
C  -------------

C  include file of dimensiopns and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of solution variables (local to this module)

      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  include files of linearized setup/solution/multiplier variables (input)

      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_L_MULTIPLIERS.VARS'

C  input arguments
C  ---------------

C  inclusion flags

      LOGICAL          DO_INCLUDE_DIRECTBEAM
      LOGICAL          DO_INCLUDE_SURFACE
      LOGICAL          DO_INCLUDE_SURFEMISS
      LOGICAL          DO_INCLUDE_THERMEMISS

C  Fourier component and beam number

      INTEGER          FOURIER_COMPONENT, IBEAM

C  surface factor = 1+delta(m,0)

      DOUBLE PRECISION SURFACE_FACTOR

C  Linearization control

      INTEGER          VARIATION_INDEX
      INTEGER          N_WEIGHTFUNCS

C  output status
C  -------------

      INTEGER          STATUS

C  Local variables
C  ---------------

C  boundary condition flags

      LOGICAL          MODIFIED_BCL3, MODIFIED_BCL4

C  error tracing variables

      INTEGER          INFO
      CHARACTER*3      CN, CI
      CHARACTER*70     MAIL, TRACE

C  Other local help variables 

      INTEGER          I, I1, Q, N, AA

C  Initialise

      STATUS = LIDORT_SUCCESS
      
C  Linearization of the regular BVP case
C  =====================================

C  Set up the column vectors for Bulk/profile linearizations
C  ---------------------------------------------------------

C  Bulk: Compute the main column B' where AX = B'

      IF ( DO_COLUMN_LINEARIZATION ) THEN

        CALL L_BVP_C_COLUMN_SETUP
     I     ( DO_INCLUDE_SURFACE,   DO_INCLUDE_DIRECTBEAM,
     I       DO_INCLUDE_SURFEMISS, DO_INCLUDE_THERMEMISS,
     I       N_WEIGHTFUNCS,
     I       FOURIER_COMPONENT,  IBEAM,
     I       SURFACE_FACTOR )

C  Profile: Boundary condition flags for special cases
C  Profile: Compute the main column B' where AX = B'

      ELSE IF ( DO_PROFILE_LINEARIZATION ) THEN

C  Boundary condition flags for special cases

        MODIFIED_BCL3 = ( VARIATION_INDEX .EQ. 1 )
        MODIFIED_BCL4 = ( VARIATION_INDEX .EQ. NLAYERS )

        CALL L_BVP_P_COLUMN_SETUP
     I     ( DO_INCLUDE_SURFACE,   DO_INCLUDE_DIRECTBEAM,
     I       DO_INCLUDE_SURFEMISS, DO_INCLUDE_THERMEMISS,
     I       MODIFIED_BCL3,      MODIFIED_BCL4,
     I       VARIATION_INDEX,    N_WEIGHTFUNCS,
     I       FOURIER_COMPONENT,  IBEAM,
     I       SURFACE_FACTOR )

      ENDIF

C  BVP back-substitution: With compression (multilayers)
C  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

C  LAPACK substitution (DGBTRS) using RHS column vector COL2_WF
C  BV solution for perturbed integration constants
C    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS
     &     ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, N_WEIGHTFUNCS,
     &        BANDMAT2, MAXBANDTOTAL, IPIVOT,
     &        COL2_WF, MAXTOTAL, INFO )

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(CN, '(I3)' ) VARIATION_INDEX
          MAIL  = 'argument i illegal value, for i = '//CI
          TRACE = 'Atmos_Wfs for layer '//CN//
     *         '; (multilayer) DGBTRS call in L_BVP_SOLUTION_MASTER'
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ENDIF

C  Set integration constants NCON and PCON for +/- eigensolutions

        DO N = 1, NLAYERS
          DO I = 1, NSTREAMS
            DO Q = 1, N_WEIGHTFUNCS
              NCON(I,N,Q) = COL2_WF(LCONMASK(I,N),Q)
              PCON(I,N,Q) = COL2_WF(MCONMASK(I,N),Q)
            ENDDO
          ENDDO
        ENDDO

C  Solve the boundary problem: No compression, Single Layer only
C  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

C  LAPACK substitution (DGETRS) using RHS column vector SCOL2_WF

        CALL DGETRS
     &     ( 'N', NTOTAL, N_WEIGHTFUNCS, SMAT2, MAXSTREAMS_2, SIPIVOT,
     &        SCOL2_WF, MAXSTREAMS_2, INFO )

C  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MAIL  = 'argument i illegal value, for i = '//CI
          TRACE = 
     &    'Atmos_Wfs for 1-layer: DGETRS call in L_BVP_SOLUTION_MASTER'
          STATUS = LIDORT_SERIOUS
          CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ENDIF

C  Set integration constants NCON and PCON for +/- eigensolutions

        N = 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO Q = 1, N_WEIGHTFUNCS
            NCON(I,N,Q) = SCOL2_WF(I,Q)
            PCON(I,N,Q) = SCOL2_WF(I1,Q)
          ENDDO
        ENDDO

      ENDIF

C  linearized BVP results
C  ======================

C  Associated quantities

      DO N = 1, NLAYERS
        DO I = 1, NSTREAMS_2
          DO AA = 1, NSTREAMS
             DO Q = 1, N_WEIGHTFUNCS
              NCON_XVEC(I,AA,N,Q) = NCON(AA,N,Q) * XPOS(I,AA,N)
              PCON_XVEC(I,AA,N,Q) = PCON(AA,N,Q) * XNEG(I,AA,N)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

C  debug

c      if ( variation_index .eq. 11 ) then
c       if ( fourier_component .EQ.3 .and.ibeam.eq.1 ) then
c        do n = 1, nlayers
c          write(98,'(a,2i3,1p2e18.10)')
c     &     'hey',n,variation_index,ncon(3,n,1),pcon(3,n,1)
c        enddo
c       endif
c      endif

C  debug------------------------------------------
c        if ( do_debug_write ) then
c         IF (ibeam .EQ. 1.and.fourier_component.eq.0 ) THEN
c          write(96,*)'newbie'
c          DO N = 1, NLAYERS
c           DO AA = 1, NSTREAMS
c            write(96,'(3i2,1p6e17.9)')FOURIER_COMPONENT,N,AA,
c     &                LCON(AA,N),MCON(AA,N),
c     &                NCON(AA,N,1),PCON(AA,N,1)
c           ENDDO
c          ENDDO
c         ENDIF
c        endif

C  debug

c      q = 96
c      if ( fourier_component.le.3.and.ibeam.eq.1
c     &          .and.variation_index.eq.4) then
c        write(96,*)fourier_component
c        DO N = 1, NLAYERS
c          write(q,'(2i4,1p4e17.9)')IBEAM,N, NCON(5,N,1), PCON(5,N,1)
c        ENDDO
c      ENDIF

C  finish

      RETURN
      END

C

      SUBROUTINE L_BVP_P_COLUMN_SETUP
     I     ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM,
     I       DO_INCLUDE_SURFEMISS, DO_INCLUDE_THERMEMISS,
     I       MODIFIED_BCL3,      MODIFIED_BCL4,
     I       LAYER_TO_VARY,      N_LAYER_WFS,
     I       FOURIER_COMPONENT,  IPARTIC,
     I       SURFACE_FACTOR )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and reflectance variables

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  include files of linearized setup variables (input)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'

C  include files of linearized solution variables (output stored here)

      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_L_THERMALSUP.VARS'

C  input arguments
C  ---------------

C  inclusion flags

      LOGICAL          DO_INCLUDE_DIRECTBEAM
      LOGICAL          DO_INCLUDE_SURFACE
      LOGICAL          DO_INCLUDE_SURFEMISS
      LOGICAL          DO_INCLUDE_THERMEMISS

C  boundary condition flags

      LOGICAL          MODIFIED_BCL3, MODIFIED_BCL4

C  Fourier component and beam number

      INTEGER          FOURIER_COMPONENT, IPARTIC

C  surface factor = 1+delta(m,0)

      DOUBLE PRECISION SURFACE_FACTOR

C  Linearization control

      INTEGER          LAYER_TO_VARY
      INTEGER          N_LAYER_WFS

C  local variables
C  ---------------

      INTEGER          Q,AA,N,N1,I,I1,CM,C0,K
      DOUBLE PRECISION CPOS, CNEG, L_HOM, L_BEAM, FAC, L_PARTIC

      DOUBLE PRECISION R2_L_PARTIC(MAXSTREAMS,MAX_ATMOSWFS)
      DOUBLE PRECISION R2_L_HOMP(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      DOUBLE PRECISION R2_L_HOMM(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)

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

      K = LAYER_TO_VARY

C  Copy already existing thermal linearizations

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        DO I = 1, NSTREAMS_2
          DO Q = 1, N_LAYER_WFS
            L_WUPPER(I,K,Q) = L_T_WUPPER(I,K,Q)
            L_WLOWER(I,K,Q) = L_T_WLOWER(I,K,Q)
          ENDDO
        ENDDO
      ELSE
        DO I = 1, NSTREAMS_2
          DO Q = 1, N_LAYER_WFS
            L_WUPPER(I,K,Q) = ZERO
            L_WLOWER(I,K,Q) = ZERO
          ENDDO
        ENDDO
      ENDIF
      
C  Get the linearized beam solution for the varying layer

      IF ( DO_SOLAR_SOURCES ) THEN
        CALL L_BEAMSOLUTION_P_NEQK
     I   ( FOURIER_COMPONENT, IPARTIC, LAYER_TO_VARY, N_LAYER_WFS )
      ENDIF
      
C  complete boundary condition flags

      REGULAR_BCL3 = .NOT.MODIFIED_BCL3
      REGULAR_BCL4 = .NOT.MODIFIED_BCL4

C  BCL1 or BCL3M - top of first layer (TOA), UPPER boundary condition
C  ------------------------------------------------------------------

      N = 1

C    If this layer is the one that is varied, use MODIFIED_BCL3 (BCL3M)

      IF ( MODIFIED_BCL3 ) THEN

C  .. contribution WVAR from beam solution variations
C  .. contribution HVAR homogeneous (eigenvalue) solution variations

        DO I = 1, NSTREAMS
          DO Q = 1, N_LAYER_WFS
            L_PARTIC = - L_WUPPER(I,N,Q)
            L_HOM  = ZERO
            DO AA = 1, NSTREAMS
              CPOS = L_XPOS(I,AA,N,Q)
              CNEG = T_DELT_EIGEN(AA,N)   * L_XNEG(I,AA,N,Q) + 
     &             L_T_DELT_EIGEN(AA,N,Q) *   XNEG(I,AA,N)
              L_HOM = L_HOM + LCON(AA,N) * CPOS + MCON(AA,N) * CNEG
            ENDDO
            COL2_WF(I,Q) = L_PARTIC - L_HOM
          ENDDO
        ENDDO

C  No variation case (BCL1)

      ELSE

        DO I = 1, NSTREAMS
          DO Q = 1, N_LAYER_WFS
            COL2_WF(I,Q) = ZERO
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
          C0  = N1*NSTREAMS_2 - NSTREAMS
          DO I = 1, NSTREAMS_2
            CM = C0 + I
            DO Q = 1, N_LAYER_WFS
              COL2_WF(CM,Q) = ZERO
            ENDDO
          ENDDO
        ENDDO

      ENDIF

C  BCL3 - regular upper boundary condition for layer that is varying
C  -----------------------------------------------------------------

      IF ( REGULAR_BCL3 ) THEN

        N = LAYER_TO_VARY
        N1  = N - 1
        C0  = N1*NSTREAMS_2 - NSTREAMS

C  .. contribution WVAR from beam solution variations
C  .. contribution HVAR homogeneous (eigenvalue) solution variations

        DO I = 1, NSTREAMS_2
          CM = C0 + I
          DO Q = 1, N_LAYER_WFS
            L_PARTIC  = + L_WUPPER(I,N,Q)
            L_HOM = ZERO
            DO AA = 1, NSTREAMS
              CPOS = L_XPOS(I,AA,N,Q)
              CNEG = T_DELT_EIGEN(AA,N)   * L_XNEG(I,AA,N,Q) + 
     &             L_T_DELT_EIGEN(AA,N,Q) *   XNEG(I,AA,N)
              L_HOM = L_HOM + LCON(AA,N) * CPOS + MCON(AA,N) * CNEG
            ENDDO
            COL2_WF(CM,Q) = L_PARTIC + L_HOM
          ENDDO
        ENDDO

      ENDIF
 
C  BCL4 - LOWER boundary condition for varying layer
C  -------------------------------------------------

C   special case when layer-to-vary = last (albedo) layer is treated
C   separately below under MODIFIED BCL4.

      IF ( REGULAR_BCL4 ) THEN

        N  = LAYER_TO_VARY
        N1 = N + 1
        C0 = N*NSTREAMS_2 - NSTREAMS

C  Get the linearized beam solution for the next layer

        IF ( DO_SOLAR_SOURCES ) THEN
          CALL L_BEAMSOLUTION_P_NNEK
     I  ( FOURIER_COMPONENT, IPARTIC, N1, LAYER_TO_VARY, N_LAYER_WFS )
        ENDIF
   
C  .. 2 contributions to WVAR from beam solution variations BEAM_V and BEAM_U 
C  .. contribution HVAR homogeneous (eigenvalue) solution variations

        DO I = 1, NSTREAMS_2
          CM = C0 + I
          DO Q = 1, N_LAYER_WFS
            IF ( DO_SOLAR_SOURCES ) THEN
              L_PARTIC = L_WUPPER(I,N1,Q) - L_WLOWER(I,N,Q)
            ELSE
              L_PARTIC = - L_WLOWER(I,N,Q)
            ENDIF
            L_HOM  = ZERO
            DO AA = 1, NSTREAMS
              CNEG = L_XNEG(I,AA,N,Q)
              CPOS = T_DELT_EIGEN(AA,N)   * L_XPOS(I,AA,N,Q) + 
     &             L_T_DELT_EIGEN(AA,N,Q) *   XPOS(I,AA,N)
              L_HOM = L_HOM + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
            ENDDO
            COL2_WF(CM,Q) = L_PARTIC - L_HOM
          ENDDO
        ENDDO

      ENDIF

C  BCL5 - Intermediate boundary conditions between varying layer & final layer
C  ---------------------------------------------------------------------------

      IF ( REGULAR_BCL4 ) THEN

        DO N = LAYER_TO_VARY + 1, NLAYERS - 1

          N1 = N + 1
          C0  = N*NSTREAMS_2 - NSTREAMS

C  Get the linearized beam solution for the next layer

          IF ( DO_SOLAR_SOURCES ) THEN
            CALL L_BEAMSOLUTION_P_NNEK
     I  ( FOURIER_COMPONENT, IPARTIC, N1, LAYER_TO_VARY, N_LAYER_WFS )
          ENDIF
          
C  .. contributions from beam solution (direct assign). No homog. variation

          DO I = 1, NSTREAMS_2
            CM = C0 + I
            DO Q = 1, N_LAYER_WFS
            IF ( DO_SOLAR_SOURCES ) THEN
              COL2_WF(CM,Q) = L_WUPPER(I,N1,Q) - L_WLOWER(I,N,Q)
            ELSE
              COL2_WF(CM,Q) = - L_WLOWER(I,N,Q)
            ENDIF
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
     I     ( DO_INCLUDE_SURFACE,
     I       MODIFIED_BCL4, IPARTIC,
     I       FOURIER_COMPONENT, SURFACE_FACTOR, N_LAYER_WFS,
     O       R2_L_PARTIC, R2_L_HOMP, R2_L_HOMM )

C  Compute the solution

        C0 = (N-1)*NSTREAMS_2 + NSTREAMS
        DO I = 1, NSTREAMS
          CM = C0 + I
          I1 = I + NSTREAMS
          DO Q = 1, N_LAYER_WFS
            L_PARTIC = L_WLOWER(I1,N,Q) - R2_L_PARTIC(I,Q)
            L_HOM  = ZERO
            DO AA = 1, NSTREAMS
              CPOS = T_DELT_EIGEN(AA,N)   * L_XPOS(I1,AA,N,Q) + 
     &             L_T_DELT_EIGEN(AA,N,Q) *   XPOS(I1,AA,N)
              CPOS =        CPOS       - R2_L_HOMP(I,AA,Q)
              CNEG = L_XNEG(I1,AA,N,Q) - R2_L_HOMM(I,AA,Q)
              L_HOM = L_HOM + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
            ENDDO
            COL2_WF(CM,Q) = - L_PARTIC - L_HOM
         ENDDO
        ENDDO

C  ordinary BCL6 Component loop
 
      ELSE

C  get the linearized downward-reflected term
C    ------- Only for the solar functions

        IF ( DO_SOLAR_SOURCES ) THEN
          CALL L_BVP_SURFACE_SETUP
     I     ( DO_INCLUDE_SURFACE,
     I       MODIFIED_BCL4, IPARTIC,
     I       FOURIER_COMPONENT, SURFACE_FACTOR, N_LAYER_WFS,
     O       R2_L_PARTIC, R2_L_HOMP, R2_L_HOMM )
        ENDIF
        
C  Compute the solution. Only present for the solar terms

        C0 = (N-1)*NSTREAMS_2 + NSTREAMS
        IF ( DO_SOLAR_SOURCES ) THEN
          DO I = 1, NSTREAMS
            CM = C0 + I
            I1 = I + NSTREAMS
            DO Q = 1, N_LAYER_WFS
              L_PARTIC = L_WLOWER(I1,N,Q) - R2_L_PARTIC(I,Q)              
              COL2_WF(CM,Q) = - L_PARTIC
            ENDDO
          ENDDO
        ENDIF
        
      ENDIF

C  Add direct beam variation to Final boundary
C  -------------------------------------------

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        DO I = 1, NSTREAMS
          CM = C0 + I
          FAC = - DIRECT_BEAM(I,IPARTIC) * 
     &                DELTAU_SLANT(N,LAYER_TO_VARY,IPARTIC) 
          DO Q = 1, N_LAYER_WFS
            L_BEAM = L_DELTAU_VERT(Q,LAYER_TO_VARY) * FAC
            COL2_WF(CM,Q) = COL2_WF(CM,Q) + L_BEAM
          ENDDO
        ENDDO
      ENDIF

C  Copy for the one-layer case

      IF ( NLAYERS .EQ. 1 ) THEN
        DO I = 1, NTOTAL
          DO Q = 1, N_LAYER_WFS
            SCOL2_WF(I,Q) = COL2_WF(I,Q)
          ENDDO
        ENDDO
      ENDIF
        
C  debug

c      if ( do_debug_write ) then
c       if ( layer_to_vary.eq.4.and.fourier_component.eq.3) then
c        DO N = 1, NTOTAL
c          write(95,'(3i4,1p4e17.9)')FOURIER_COMPONENT,IPARTIC,N,
c     &                            COL2_WF(N,1)
c        ENDDO
c       pause
c       endif
c      ENDIF

C  finish

      RETURN
      END

C

      SUBROUTINE L_BVP_C_COLUMN_SETUP
     I     ( DO_INCLUDE_SURFACE,
     I       DO_INCLUDE_DIRECTBEAM,
     I       DO_INCLUDE_SURFEMISS,
     I       DO_INCLUDE_THERMEMISS,
     I       N_WEIGHTFUNCS,
     I       FOURIER_COMPONENT, IPARTIC,
     I       SURFACE_FACTOR )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and reflectance variables

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  include files of linearized setup variables (input)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'

C  include files of linearized solution variables (output stored here)

      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_L_THERMALSUP.VARS'

C  input arguments
C  ---------------

      LOGICAL          DO_INCLUDE_SURFACE
      LOGICAL          DO_INCLUDE_DIRECTBEAM
      LOGICAL          DO_INCLUDE_SURFEMISS
      LOGICAL          DO_INCLUDE_THERMEMISS

      INTEGER          FOURIER_COMPONENT
      DOUBLE PRECISION SURFACE_FACTOR
      INTEGER          N_WEIGHTFUNCS
      INTEGER          IPARTIC

C  local variables
C  ---------------

      INTEGER          Q, AA, N, N1, I, I1, CM, C0, K
      DOUBLE PRECISION CPOS, CNEG, L_HOM, L_PARTIC
      DOUBLE PRECISION L_BEAM, L_HOMD, L_HOMU
      LOGICAL          MODIFIED_BOUNDARY

      DOUBLE PRECISION R2_L_PARTIC(MAXSTREAMS,MAX_ATMOSWFS), FAC, FAC3
      DOUBLE PRECISION R2_L_HOMP(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      DOUBLE PRECISION R2_L_HOMM(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)

C  initialise
C  ----------

C  zero the results vectors

      DO I = 1, NTOTAL
        DO Q = 1, MAX_ATMOSWFS
          COL2_WF(I,Q) = ZERO
        ENDDO
      ENDDO

C  Copy already existing thermal linearizations

      DO K = 1, NLAYERS
       IF ( DO_INCLUDE_THERMEMISS ) THEN
        DO I = 1, NSTREAMS_2
          DO Q = 1, N_WEIGHTFUNCS
            L_WUPPER(I,K,Q) = L_T_WUPPER(I,K,Q)
            L_WLOWER(I,K,Q) = L_T_WLOWER(I,K,Q)
          ENDDO
        ENDDO
       ELSE
        DO I = 1, NSTREAMS_2
          DO Q = 1, N_WEIGHTFUNCS
            L_WUPPER(I,K,Q) = ZERO
            L_WLOWER(I,K,Q) = ZERO
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
     I     ( FOURIER_COMPONENT, IPARTIC, N, N_WEIGHTFUNCS )
      ENDIF

C  .. contribution WVAR from beam solution variations
C  .. contribution HVAR homogeneous (eigenvalue) solution variations

      DO I = 1, NSTREAMS
        DO Q = 1, N_WEIGHTFUNCS
          L_PARTIC = - L_WUPPER(I,N,Q)
          L_HOM    = ZERO
          DO AA = 1, NSTREAMS
            CPOS = L_XPOS(I,AA,N,Q)
            CNEG =   T_DELT_EIGEN(AA,N)   * L_XNEG(I,AA,N,Q) + 
     &             L_T_DELT_EIGEN(AA,N,Q) *   XNEG(I,AA,N)
            L_HOM = L_HOM + LCON(AA,N) * CPOS + MCON(AA,N) * CNEG
          ENDDO
          COL2_WF(I,Q) = L_PARTIC - L_HOM
        ENDDO
      ENDDO
 
C  Intermediate boundary conditions
C  --------------------------------

      DO N = 1, NLAYERS - 1

C  N1 is the layer below, C0 is the offset

        N1 = N + 1
        C0 = N*NSTREAMS_2 - NSTREAMS

C  Get the linearized beam solution for the next layer

        IF ( DO_SOLAR_SOURCES ) THEN
          CALL L_BEAMSOLUTION_C_NEQK
     I   ( FOURIER_COMPONENT, IPARTIC, N1, N_WEIGHTFUNCS )
        ENDIF

C  .. 2 contributions to L_BEAM, from variations L_WUPPER L_WLOWER 
C  .. 2 contributions to L_HOM,  from variations above and below

        DO I = 1, NSTREAMS_2
          CM = C0 + I
          DO Q = 1, N_WEIGHTFUNCS
            L_HOMD = ZERO
            L_HOMU = ZERO
            DO AA = 1, NSTREAMS
              CPOS = L_XPOS(I,AA,N1,Q)
              CNEG =   T_DELT_EIGEN(AA,N1)   * L_XNEG(I,AA,N1,Q) + 
     &               L_T_DELT_EIGEN(AA,N1,Q) *   XNEG(I,AA,N1)
              L_HOMU = L_HOMU + LCON(AA,N1) * CPOS + MCON(AA,N1) * CNEG
              CNEG = L_XNEG(I,AA,N,Q)
              CPOS =   T_DELT_EIGEN(AA,N)   * L_XPOS(I,AA,N,Q) + 
     &               L_T_DELT_EIGEN(AA,N,Q) *   XPOS(I,AA,N)
              L_HOMD = L_HOMD + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
            ENDDO
            L_PARTIC      = L_WUPPER(I,N1,Q) - L_WLOWER(I,N,Q)      
            L_HOM         = L_HOMU - L_HOMD
            COL2_WF(CM,Q) = L_PARTIC + L_HOM
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
     I     ( DO_INCLUDE_SURFACE, MODIFIED_BOUNDARY, IPARTIC,
     I       FOURIER_COMPONENT, SURFACE_FACTOR, N_WEIGHTFUNCS,
     O       R2_L_PARTIC, R2_L_HOMP, R2_L_HOMM )

C  Compute the solution

      C0 = (N-1)*NSTREAMS_2 + NSTREAMS
      DO I = 1, NSTREAMS
        CM = C0 + I
        I1 = I + NSTREAMS
        DO Q = 1, N_WEIGHTFUNCS
          L_PARTIC = L_WLOWER(I1,N,Q) - R2_L_PARTIC(I,Q)
          L_HOM    = ZERO
          DO AA = 1, NSTREAMS
            CPOS =   T_DELT_EIGEN(AA,N)   * L_XPOS(I1,AA,N,Q) + 
     &             L_T_DELT_EIGEN(AA,N,Q) *   XPOS(I1,AA,N)
            CPOS =        CPOS       - R2_L_HOMP(I,AA,Q)
            CNEG = L_XNEG(I1,AA,N,Q) - R2_L_HOMM(I,AA,Q)
            L_HOM = L_HOM + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
          ENDDO
          COL2_WF(CM,Q) = - L_PARTIC - L_HOM
        ENDDO
      ENDDO

C  Add direct beam variation to Final boundary
C  -------------------------------------------

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        DO I = 1, NSTREAMS
          CM = C0 + I
          FAC = - DIRECT_BEAM(I,IPARTIC)
          DO Q = 1, N_WEIGHTFUNCS
            L_BEAM = ZERO
            DO K = 1, NLAYERS
              FAC3 = FAC * DELTAU_SLANT(N,K,IPARTIC)
              L_BEAM = L_BEAM + L_DELTAU_VERT(Q,K) * FAC3
            ENDDO
            COL2_WF(CM,Q) = COL2_WF(CM,Q) + L_BEAM
          ENDDO
        ENDDO
      ENDIF

C  finish

      RETURN
      END

C

      SUBROUTINE L_BVP_SURFACE_SETUP
     I     ( DO_INCLUDE_SURFACE, MODIFIED_BCL4, IPARTIC,
     I       FOURIER_COMPONENT, SURFACE_FACTOR, N_LAYER_WFS,
     O       R2_L_PARTIC, R2_L_HOMP, R2_L_HOMM )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and reflectance variables

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  include files of linearized setup variables (input)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'

C  include files of linearized solution variables (output stored here)

      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'

C  input arguments
C  ---------------

C  Number of weighting functions

      INTEGER          N_LAYER_WFS

C  Fourier component and PIS index (no strictly needed)

      INTEGER          FOURIER_COMPONENT, IPARTIC

C  Flag for type of boundary condition

      LOGICAL          MODIFIED_BCL4

C  overall surface flag and surface factor

      LOGICAL          DO_INCLUDE_SURFACE
      DOUBLE PRECISION SURFACE_FACTOR

C  Output arguments
C  ----------------

      DOUBLE PRECISION R2_L_PARTIC(MAXSTREAMS,MAX_ATMOSWFS)
      DOUBLE PRECISION R2_L_HOMP(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      DOUBLE PRECISION R2_L_HOMM(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)

C  Local variables
C  ---------------

      DOUBLE PRECISION PV_W(MAXSTREAMS,MAX_ATMOSWFS)
      DOUBLE PRECISION HV_P(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      DOUBLE PRECISION HV_M(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      DOUBLE PRECISION HSP_U, HSM_U, REFL_P, REFL_B, REFL_M
      DOUBLE PRECISION FACTOR_KERNEL
      INTEGER          AA, J, Q, N, I, K

C  Initial section
C  ---------------

C  Always zero the result to start

      DO I = 1, NSTREAMS
        DO Q = 1, N_LAYER_WFS
          R2_L_PARTIC(I,Q)    = ZERO
          DO AA = 1, NSTREAMS
            R2_L_HOMP(I,AA,Q) = ZERO
            R2_L_HOMM(I,AA,Q) = ZERO
          ENDDO
        ENDDO
      ENDDO

C  Return if no albedo

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

C  Set up Auxiliary arrays
C  -----------------------

      N = NLAYERS

C  Particular integral parts

      DO J = 1, NSTREAMS
        DO Q = 1, N_LAYER_WFS
          PV_W(J,Q) = L_WLOWER(J,N,Q) * AX(J)
        ENDDO
      ENDDO

C    Modified boundary condition: homogeneous parts

      IF ( MODIFIED_BCL4 ) THEN
        DO J = 1, NSTREAMS
          DO AA = 1, NSTREAMS
            DO Q = 1, N_LAYER_WFS
              HSP_U = L_XPOS(J,AA,N,Q) *   T_DELT_EIGEN(AA,N) +
     &                  XPOS(J,AA,N)   * L_T_DELT_EIGEN(AA,N,Q)
              HSM_U = L_XNEG(J,AA,N,Q)
              HV_P(J,AA,Q) = AX(J)*HSP_U
              HV_M(J,AA,Q) = AX(J)*HSM_U
            ENDDO
          ENDDO
        ENDDO
      ENDIF

C  Integrated Downward reflection (Calculation)
C  --------------------------------------------

C  Start loop over BRDF kernels

      DO K = 1, N_BRDF_KERNELS

C  Kernel amplitude

        FACTOR_KERNEL = SURFACE_FACTOR * BRDF_FACTORS(K)

C  Lambertian Reflection
C  =====================

        IF ( LAMBERTIAN_KERNEL_FLAG(K) ) THEN
          IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
            DO Q = 1, N_LAYER_WFS

C  Particular solution

              REFL_B = ZERO
              DO J = 1, NSTREAMS
                REFL_B = REFL_B + PV_W(J,Q)
              ENDDO
              REFL_B = REFL_B * FACTOR_KERNEL
              DO I = 1, NSTREAMS
                R2_L_PARTIC(I,Q) = R2_L_PARTIC(I,Q) + REFL_B
              ENDDO

C  Homogeneous solutions (only for modified BC)

              IF ( MODIFIED_BCL4 ) THEN
               DO AA = 1, NSTREAMS
                REFL_P = ZERO
                REFL_M = ZERO
                DO J = 1, NSTREAMS
                  REFL_P = REFL_P + HV_P(J,AA,Q)
                  REFL_M = REFL_M + HV_M(J,AA,Q)
                ENDDO
                REFL_P = REFL_P * FACTOR_KERNEL
                REFL_M = REFL_M * FACTOR_KERNEL
                DO I = 1, NSTREAMS
                  R2_L_HOMP(I,AA,Q) = R2_L_HOMP(I,AA,Q) + REFL_P
                  R2_L_HOMM(I,AA,Q) = R2_L_HOMM(I,AA,Q) + REFL_M
                ENDDO
               ENDDO
              ENDIF
              
C  end parameter loop

            ENDDO
          ENDIF

C  Bidirectional reflecting kernel
C  ===============================

        ELSE

          DO I = 1, NSTREAMS
            DO Q = 1, N_LAYER_WFS

C  particular solutions

              REFL_B = ZERO
              DO J = 1, NSTREAMS
                REFL_B = REFL_B + PV_W(J,Q) * BIREFLEC(K,J,I)
              ENDDO
              REFL_B = REFL_B * FACTOR_KERNEL
              R2_L_PARTIC(I,Q) = R2_L_PARTIC(I,Q) + REFL_B

C  homogeneous solutions

              IF ( MODIFIED_BCL4 ) THEN
               DO AA = 1, NSTREAMS
                REFL_P = ZERO
                REFL_M = ZERO
                DO J = 1, NSTREAMS
                  REFL_P = REFL_P + HV_P(J,AA,Q)*BIREFLEC(K,J,I)
                  REFL_M = REFL_M + HV_M(J,AA,Q)*BIREFLEC(K,J,I)
                ENDDO
                REFL_P = REFL_P * FACTOR_KERNEL
                REFL_M = REFL_M * FACTOR_KERNEL
                R2_L_HOMP(I,AA,Q) = R2_L_HOMP(I,AA,Q) + REFL_P
                R2_L_HOMM(I,AA,Q) = R2_L_HOMM(I,AA,Q) + REFL_M
               ENDDO
              ENDIF
              
C  end parameter and stream loops

            ENDDO
          ENDDO

        ENDIF

C  End kernels loop

      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE L_BEAMSOLUTION_P_NEQK ( M, IB, N, N_LAYER_WFS )

C  Linearization of beam particular integral in one layer.

C  In this module, this is the Layer that contains the variation.
C                   N = LAYER_TO_VARY

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and multiplier variables

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  include files of linearized setup variables (input)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'

C  include files of linearized solution variables (output stored here)

      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'

C  include file of linearized Multiplier coefficients (output)

      INCLUDE '../includes/LIDORT_L_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Fourier number

      INTEGER          M

C  Beam index input

      INTEGER          IB

C  Given layer index (input)
C    This is also LAYER_TO_VARY

      INTEGER          N

C  number of varying parameters (input)

      INTEGER          N_LAYER_WFS

C  Local variables
C  ---------------

      INTEGER          AA, I, I1, Q
      DOUBLE PRECISION S_P_U, S_P_L, S_M_U, S_M_L, L_SD, L_SU 
      DOUBLE PRECISION CONST, WDEL, ZDEL, ZWDEL, T1, PQ, FQ
      DOUBLE PRECISION L_WDEL, L_ZDEL, L_ZWDEL, VAR1
      DOUBLE PRECISION LBSOL(MAXSTREAMS_2,MAX_ATMOSWFS,2)

C  No linearized particular solution beyond the cutoff layer. ALSO -
C  Nothing if the solution saving mode is on and layer is inactive
C    Exit (Solutions have already been zeroed).

      IF ((DO_SOLUTION_SAVING.AND..NOT.DO_LAYER_SCATTERING(M,N))
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
            LBSOL(I,Q,1) = CONST * L_WVEC(I,N,N,Q)
            LBSOL(I,Q,2) = WDEL  * LBSOL(I,Q,1) + VAR1*WVEC(I,N)
          ENDDO
        ENDDO
      ENDIF

C  Green's function solution
C  =========================

C  More complicated!

      IF ( .NOT. DO_CLASSICAL_SOLUTION ) THEN

C  Set up linearizations of GAMMA constants
C  ----------------------------------------

C  Distinguish two cases:
C  ..(a) quasi-spherical for n > 1
C  ..(b) plane-parallel or QS for n=1

        IF ( .NOT. DO_PLANE_PARALLEL .AND. N.GT.1 ) THEN
          DO AA = 1, NSTREAMS
            DO Q = 1, N_LAYER_WFS
              FQ = L_KEIGEN(AA,N,Q)
              PQ = L_AVERAGE_SECANT(N,N,IB,Q)
              L_GAMMA_P(AA,N,N,Q) = - GAMMA_P(AA,N) * ( FQ + PQ )
              L_GAMMA_M(AA,N,N,Q) =   GAMMA_M(AA,N) * ( FQ - PQ )
            ENDDO
          ENDDO
        ELSE
          DO AA = 1, NSTREAMS
            DO Q = 1, N_LAYER_WFS
              FQ = L_KEIGEN(AA,N,Q)
              L_GAMMA_P(AA,N,N,Q) = - GAMMA_P(AA,N) * FQ 
              L_GAMMA_M(AA,N,N,Q) =   GAMMA_M(AA,N) * FQ
            ENDDO
          ENDDO
        ENDIF

C  Linearizations of optical depth integrations
C  Linearized Green function multipliers

         WDEL    = T_DELT_MUBAR(N,IB)
         DO AA = 1, NSTREAMS
          ZDEL  =   T_DELT_EIGEN(AA,N)
          ZWDEL = ZDEL * WDEL
          DO Q = 1, N_LAYER_WFS
            L_WDEL  = L_T_DELT_MUBAR(N,N,IB,Q)
            L_ZDEL  = L_T_DELT_EIGEN(AA,N,Q)
            L_ZWDEL = ZDEL * L_WDEL + L_ZDEL * WDEL
            L_SU  =        - L_ZWDEL     / ( ONE - ZWDEL )
            IF ( GFUNC_DN(AA,N).EQ.ZERO ) THEN
              L_GFUNC_DN(AA,Q) = ZERO
            ELSE
              L_SD  =  ( L_ZDEL - L_WDEL ) / ( ZDEL - WDEL )
              T1 = L_ATERM_SAVE(AA,N,Q) + L_GAMMA_M(AA,N,N,Q)
              L_GFUNC_DN(AA,Q) = GFUNC_DN(AA,N) * ( T1 + L_SD )
            ENDIF
            T1 = L_BTERM_SAVE(AA,N,Q) + L_GAMMA_P(AA,N,N,Q)
            L_GFUNC_UP(AA,Q) = GFUNC_UP(AA,N) * ( T1 + L_SU )
          ENDDO
        ENDDO

C  Set linearized form of particular integral at boundaries

        DO Q = 1, N_LAYER_WFS
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            S_P_U = ZERO
            S_P_L = ZERO
            S_M_U = ZERO
            S_M_L = ZERO
            DO AA = 1, NSTREAMS
              S_P_U = S_P_U + L_GFUNC_UP(AA,Q) *   XPOS(I1,AA,N) +
     &                          GFUNC_UP(AA,N) * L_XPOS(I1,AA,N,Q)
              S_M_U = S_M_U + L_GFUNC_UP(AA,Q) *   XPOS(I,AA,N) +
     &                          GFUNC_UP(AA,N)  * L_XPOS(I,AA,N,Q)
              S_P_L = S_P_L + L_GFUNC_DN(AA,Q) *   XPOS(I,AA,N) +
     &                          GFUNC_DN(AA,N) * L_XPOS(I,AA,N,Q)
              S_M_L = S_M_L + L_GFUNC_DN(AA,Q) *   XPOS(I1,AA,N) +
     &                          GFUNC_DN(AA,N) * L_XPOS(I1,AA,N,Q)
            ENDDO
            LBSOL(I,Q,1)  = S_P_U
            LBSOL(I1,Q,1) = S_M_U
            LBSOL(I1,Q,2) = S_M_L
            LBSOL(I,Q,2)  = S_P_L
          ENDDO
        ENDDO

C  end clause for choice of solution method

      ENDIF

C Add to existing solution

      DO Q = 1, N_LAYER_WFS
        DO I = 1, NSTREAMS_2
          L_WUPPER(I,N,Q) = L_WUPPER(I,N,Q) + LBSOL(I,Q,1)
          L_WLOWER(I,N,Q) = L_WLOWER(I,N,Q) + LBSOL(I,Q,2)
        ENDDO
      ENDDO
      
C  Finish

      RETURN
      END

C

      SUBROUTINE L_BEAMSOLUTION_P_NNEK ( M, IB, N, K, K_PARAMETERS )

C  Linearization of beam particular integral in one layer only.
C  This is for a layer N not equal to varying layer K.

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and multiplier variables

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  include files of linearized setup variables (input)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'

C  include files of linearized solution variables (output stored here)

      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'

C  include file of linearized Multiplier coefficients (output)

      INCLUDE '../includes/LIDORT_L_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Fourier index

      INTEGER          M

C  Beam index

      INTEGER          IB

C  Given layer indices (inputs)

      INTEGER          N, K

C  number of varying parameters (input)

      INTEGER          K_PARAMETERS

C  Local variables
C  ---------------

      INTEGER          AA, I, I1, Q
      DOUBLE PRECISION S_P_U, S_P_L, S_M_U, S_M_L, L_SD, L_SU 
      DOUBLE PRECISION CONST, WDEL, ZDEL, ZWDEL, T1
      DOUBLE PRECISION PQ, L_WDEL, L_ZWDEL
      DOUBLE PRECISION TRANS2, VAR1, VAR2, VAR_U

C  No linearized particular solution beyond the cutoff layer. ALSO--
C  Nothing if the solution saving mode is on and layer is inactive
C    Zero the solutions and exit.

      IF ((DO_SOLUTION_SAVING.AND..NOT.DO_LAYER_SCATTERING(M,N))
     &       .OR. (N .GT.LAYER_PIS_CUTOFF(IB))) THEN
        DO Q = 1, K_PARAMETERS
          DO I = 1, NSTREAMS_2
            L_WUPPER(I,N,Q) = ZERO
            L_WLOWER(I,N,Q) = ZERO
          ENDDO
        ENDDO
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

C  Please note the following:
C   Bug fixed 15 August 2005 --->
C     Consistent logarithmic linearization of INITIAL_TRANS

        IF ( .NOT. DO_PLANE_PARALLEL  ) THEN
          DO Q = 1, K_PARAMETERS
            VAR1 = L_T_DELT_MUBAR(N,K,IB,Q) * CONST
            VAR2 = L_INITIAL_TRANS(N,K,IB,Q)
            DO I = 1, NSTREAMS_2
              VAR_U = VAR2 * WVEC(I,N) + L_WVEC(I,N,K,Q)
              L_WUPPER(I,N,Q) = CONST  * VAR_U
              L_WLOWER(I,N,Q) = TRANS2 * VAR_U + VAR1 * WVEC(I,N)
            ENDDO
          ENDDO
        ELSE
          DO Q = 1, K_PARAMETERS
            VAR1 = L_INITIAL_TRANS(N,K,IB,Q) * CONST
            DO I = 1, NSTREAMS_2
              L_WUPPER(I,N,Q) = VAR1 * WVEC(I,N)
              L_WLOWER(I,N,Q) = L_WUPPER(I,N,Q) * WDEL
            ENDDO
          ENDDO
        ENDIF

C  Green's function solution
C  =========================

      ELSE

C  ..(a) quasi-spherical for n > 1 (only gets done for this case anyway)

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN

C  Set up linearizations of GAMMA constants

          DO AA = 1, NSTREAMS
            DO Q = 1, K_PARAMETERS
              PQ = L_AVERAGE_SECANT(N,K,IB,Q)
              L_GAMMA_P(AA,N,K,Q) = - GAMMA_P(AA,N) * PQ
              L_GAMMA_M(AA,N,K,Q) = - GAMMA_M(AA,N) * PQ
            ENDDO
          ENDDO

C  Linearizations of optical depth integrations
C  Linearized Green function multipliers

          WDEL    = T_DELT_MUBAR(N,IB)
          DO AA = 1, NSTREAMS
            ZDEL  =   T_DELT_EIGEN(AA,N)
            ZWDEL = ZDEL * WDEL
            DO Q = 1, K_PARAMETERS
              L_WDEL  = L_T_DELT_MUBAR(N,K,IB,Q)
              L_ZWDEL = ZDEL * L_WDEL
              L_SU  = - L_ZWDEL / ( ONE - ZWDEL )
              IF ( GFUNC_DN(AA,N).EQ.ZERO ) THEN
                L_GFUNC_DN(AA,Q) = ZERO
              ELSE
                L_SD  = - L_WDEL  / ( ZDEL - WDEL )
                T1 = L_INITIAL_TRANS(N,K,IB,Q) + L_GAMMA_M(AA,N,K,Q)
                L_GFUNC_DN(AA,Q) = GFUNC_DN(AA,N) * ( T1 + L_SD )
              ENDIF
              T1 = L_INITIAL_TRANS(N,K,IB,Q) + L_GAMMA_P(AA,N,K,Q)
              L_GFUNC_UP(AA,Q) = GFUNC_UP(AA,N) * ( T1 + L_SU )
            ENDDO
          ENDDO

C  ..(b) plane-parallel and qs for n = 1

        ELSE

          DO AA = 1, NSTREAMS
            DO Q = 1, K_PARAMETERS
              T1 = L_INITIAL_TRANS(N,K,IB,Q)
              L_GFUNC_DN(AA,Q) = GFUNC_DN(AA,N) * T1
              L_GFUNC_UP(AA,Q) = GFUNC_UP(AA,N) * T1
            ENDDO
          ENDDO

        ENDIF

C  Set linearized form of particular integral

        DO Q = 1, K_PARAMETERS
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            S_P_U = ZERO
            S_P_L = ZERO
            S_M_U = ZERO
            S_M_L = ZERO
            DO AA = 1, NSTREAMS
              S_P_U = S_P_U + L_GFUNC_UP(AA,Q) *   XPOS(I1,AA,N)
              S_M_U = S_M_U + L_GFUNC_UP(AA,Q) *   XPOS(I,AA,N)
              S_P_L = S_P_L + L_GFUNC_DN(AA,Q) *   XPOS(I,AA,N)
              S_M_L = S_M_L + L_GFUNC_DN(AA,Q) *   XPOS(I1,AA,N)
            ENDDO
            L_WUPPER(I,N,Q) = S_P_U
            L_WUPPER(I1,N,Q) = S_M_U
            L_WLOWER(I1,N,Q) = S_M_L
            L_WLOWER(I,N,Q) = S_P_L
          ENDDO
        ENDDO

C  End clause for choice of method

      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE L_BEAMSOLUTION_C_NEQK ( M, IB, N, N_WEIGHTFUNCS )

C  Linearization of beam particular integral in layer N, due to column.
C   This is the bulk property linearization

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and multiplier variables

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  include files of linearized setup variables (input)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'

C  include files of linearized solution variables (output stored here)

      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'

C  include file of linearized Multiplier coefficients (output)

      INCLUDE '../includes/LIDORT_L_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Fourier number

      INTEGER          M

C  Beam index input

      INTEGER          IB

C  Given layer index (input)
C    This is also LAYER_TO_VARY

      INTEGER          N

C  number of varying parameters (input)

      INTEGER          N_WEIGHTFUNCS

C  Local variables
C  ---------------

      INTEGER          AA, I, I1, Q, LVARY
      DOUBLE PRECISION S_P_U, S_P_L, S_M_U, S_M_L, L_SD, L_SU 
      DOUBLE PRECISION CONST, WDEL, ZDEL, ZWDEL, T1, PQ, FQ, LTI
      DOUBLE PRECISION L_WDEL, L_ZDEL, L_ZWDEL, VAR1, VAR2, VAR_U
      DOUBLE PRECISION LBSOL(MAXSTREAMS_2,MAX_ATMOSWFS,2)

C  No linearized particular solution beyond the cutoff layer. ALSO -
C  Nothing if the solution saving mode is on and layer is inactive
C    Exit (Solutions have not necessarily already been zeroed).

      IF ((DO_SOLUTION_SAVING.AND..NOT.DO_LAYER_SCATTERING(M,N))
     &       .OR. (N .GT.LAYER_PIS_CUTOFF(IB))) THEN
        DO I = 1, NSTREAMS_2
          DO Q = 1, N_WEIGHTFUNCS
            L_WUPPER(I,N,Q) = ZERO
            L_WLOWER(I,N,Q) = ZERO
          ENDDO
        ENDDO
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
        DO Q = 1, N_WEIGHTFUNCS
          VAR1 = L_T_DELT_MUBAR(N,LVARY,IB,Q) * CONST
          VAR2 = L_INITIAL_TRANS(N,LVARY,IB,Q)
          DO I = 1, NSTREAMS_2
	    VAR_U = VAR2 * WVEC(I,N) + L_WVEC(I,N,LVARY,Q)
            LBSOL(I,Q,1) = CONST * VAR_U
            LBSOL(I,Q,2) = WDEL * LBSOL(I,Q,1) + VAR1 * WVEC(I,N)
          ENDDO
        ENDDO
      ENDIF

C  Green's function solution
C  =========================

C  More complicated!

      IF ( .NOT. DO_CLASSICAL_SOLUTION ) THEN

C  Set up linearizations of GAMMA constants
C  ----------------------------------------

C  Distinguish two cases:
C  ..(a) quasi-spherical for n > 1
C  ..(b) plane-parallel or QS for n=1

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN
          DO Q = 1, N_WEIGHTFUNCS
            DO AA = 1, NSTREAMS
              FQ = L_KEIGEN(AA,N,Q)
              PQ = ZERO
              IF ( LVARY.NE.0.AND.N.GT.1.OR.LVARY.EQ.0 )
     &              PQ = L_AVERAGE_SECANT(N,LVARY,IB,Q)
              L_GAMMA_P(AA,N,LVARY,Q) = - GAMMA_P(AA,N) * ( FQ + PQ )
              L_GAMMA_M(AA,N,LVARY,Q) =   GAMMA_M(AA,N) * ( FQ - PQ )
            ENDDO
          ENDDO
        ELSE
          DO AA = 1, NSTREAMS
            DO Q = 1, N_WEIGHTFUNCS
              FQ = L_KEIGEN(AA,N,Q)
              L_GAMMA_P(AA,N,LVARY,Q) = - GAMMA_P(AA,N) * FQ 
              L_GAMMA_M(AA,N,LVARY,Q) =   GAMMA_M(AA,N) * FQ
            ENDDO
          ENDDO
        ENDIF

C  Linearizations of optical depth integrations
C  Linearized Green function multipliers

         CONST   = INITIAL_TRANS(N,IB)
         WDEL    = T_DELT_MUBAR(N,IB)
         DO AA = 1, NSTREAMS
          ZDEL  =   T_DELT_EIGEN(AA,N)
          ZWDEL = ZDEL * WDEL
          DO Q = 1, N_WEIGHTFUNCS
            L_WDEL  = L_T_DELT_MUBAR(N,LVARY,IB,Q)
            L_ZDEL  = L_T_DELT_EIGEN(AA,N,Q)
            L_ZWDEL = ZDEL * L_WDEL + L_ZDEL * WDEL
            L_SU  =        - L_ZWDEL     / ( ONE - ZWDEL )
            IF ( GFUNC_DN(AA,N).EQ.ZERO ) THEN
              L_GFUNC_DN(AA,Q) = ZERO
            ELSE
              L_SD  =  ( L_ZDEL - L_WDEL ) / ( ZDEL - WDEL )
              T1 = L_ATERM_SAVE(AA,N,Q) + L_GAMMA_M(AA,N,LVARY,Q)
              L_GFUNC_DN(AA,Q) = GFUNC_DN(AA,N) * ( T1 + L_SD )
            ENDIF
            T1 = L_BTERM_SAVE(AA,N,Q) + L_GAMMA_P(AA,N,LVARY,Q)
            L_GFUNC_UP(AA,Q) = GFUNC_UP(AA,N) * ( T1 + L_SU )
            IF ( CONST .NE. ZERO ) THEN
              LTI = L_INITIAL_TRANS(N,LVARY,IB,Q)
	      L_GFUNC_DN(AA,Q) = L_GFUNC_DN(AA,Q)+GFUNC_DN(AA,N)*LTI
	      L_GFUNC_UP(AA,Q) = L_GFUNC_UP(AA,Q)+GFUNC_UP(AA,N)*LTI
            ENDIF
          ENDDO
        ENDDO

C  Set linearized form of particular integral at boundaries

        DO Q = 1, N_WEIGHTFUNCS
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            S_P_U = ZERO
            S_P_L = ZERO
            S_M_U = ZERO
            S_M_L = ZERO
            DO AA = 1, NSTREAMS
              S_P_U = S_P_U + L_GFUNC_UP(AA,Q) *   XPOS(I1,AA,N) +
     &                          GFUNC_UP(AA,N) * L_XPOS(I1,AA,N,Q)
              S_M_U = S_M_U + L_GFUNC_UP(AA,Q) *   XPOS(I,AA,N) +
     &                          GFUNC_UP(AA,N)  * L_XPOS(I,AA,N,Q)
              S_P_L = S_P_L + L_GFUNC_DN(AA,Q) *   XPOS(I,AA,N) +
     &                          GFUNC_DN(AA,N) * L_XPOS(I,AA,N,Q)
              S_M_L = S_M_L + L_GFUNC_DN(AA,Q) *   XPOS(I1,AA,N) +
     &                          GFUNC_DN(AA,N) * L_XPOS(I1,AA,N,Q)
            ENDDO
            LBSOL(I,Q,1)  = S_P_U
            LBSOL(I1,Q,1) = S_M_U
            LBSOL(I1,Q,2) = S_M_L
            LBSOL(I,Q,2)  = S_P_L
          ENDDO
        ENDDO

C  end clause for choice of solution method

      ENDIF

C Add to existing solution

      DO Q = 1, N_WEIGHTFUNCS
        DO I = 1, NSTREAMS_2
          L_WUPPER(I,N,Q) = L_WUPPER(I,N,Q) + LBSOL(I,Q,1)
          L_WLOWER(I,N,Q) = L_WLOWER(I,N,Q) + LBSOL(I,Q,2)
        ENDDO
      ENDDO
     
C  Finish

      RETURN
      END

C

      SUBROUTINE KFACTOR_WF_COLSETUP
     I   (  DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFEMISS,
     I      SURFACE_FACTOR, IBEAM_INDEX, BRDF_INDEX )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and reflectance variables

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  include files of linearized solution variables (output stored here)

      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'

C  inputs
C  ------

C  .. direct beam inputs

      LOGICAL          DO_INCLUDE_DIRECTBEAM

C  .. surface emission inputs

      LOGICAL          DO_INCLUDE_SURFEMISS

C  surface factor = 1 + delta(m,0)

      DOUBLE PRECISION SURFACE_FACTOR

C  Beam index

      INTEGER          IBEAM_INDEX
  
C  BRDF index (indicates which albedo kernel)

      INTEGER          BRDF_INDEX


C  local variables
C  ---------------

      INTEGER          N, I, C0, CM, AA, N1, IB, B
      DOUBLE PRECISION HELP, EMISS_VAR

C  initialise

      DO I = 1, NTOTAL
        COL2_WFALB(I,1) = ZERO
      ENDDO

C  boundary conditions not changed for first layer upper (TOA)

      N = 1
      DO I = 1, NSTREAMS
        COL2_WFALB(I,1) = ZERO
      ENDDO

C  boundary conditions not changed for all intermediate levels

      DO N = 2, NLAYERS - 1
        N1 = N - 1
        C0  = N1*NSTREAMS_2 - NSTREAMS
        DO I = 1, NSTREAMS_2
          CM = C0 + I
          COL2_WFALB(CM,1) = ZERO
        ENDDO
      ENDDO

C  Ground level boundary condition
C  -------------------------------

C  Initialise

      N  = NLAYERS
      C0 = (N-1)*NSTREAMS_2 + NSTREAMS
      B  = BRDF_INDEX
      IB = IBEAM_INDEX

C  Diffuse scatter contributions

      DO I = 1, NSTREAMS
        CM = C0 + I
        HELP = RA_PARTIC(I,B)
        DO AA = 1, NSTREAMS
          HELP = HELP
     &        + LCON(AA,N)*T_DELT_EIGEN(AA,N)*RA_HOMP(I,AA,B)
     &        + MCON(AA,N)*RA_HOMM(I,AA,B)
        ENDDO
        COL2_WFALB(CM,1) = HELP
      ENDDO

C  Add direct beam variation of albedo

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        DO I = 1, NSTREAMS
          CM = C0 + I
          COL2_WFALB(CM,1) = COL2_WFALB(CM,1) + A_DIRECT_BEAM(I,IB,B)
        ENDDO
      ENDIF

C  If surface emission, include emissivity variation

      IF ( DO_INCLUDE_SURFEMISS ) THEN
        DO I = 1, NSTREAMS
          CM = C0 + I
          EMISS_VAR = ALBEDO_EMISSIVITY(I,B) * SURFBB
          COL2_WFALB(CM,1) = COL2_WFALB(CM,1) - EMISS_VAR
        ENDDO
      ENDIF

C  Copy for the single layer case

      IF ( NLAYERS .EQ. 1 ) THEN
        DO N = 1, NTOTAL
          SCOL2_WFALB(N,1) = COL2_WFALB(N,1)
        ENDDO
      ENDIF

C  debug

      if ( do_debug_write ) then
        DO N = 1, NTOTAL
          write(85,'(2i4,1p4e17.9)')IBEAM_INDEX,N,COL2_WFALB(N,1)
        ENDDO
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE KPARAMS_WF_COLSETUP
     I   (  DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFEMISS,
     I      SURFACE_FACTOR, IBEAM_INDEX, BRDF_INDEX, PAR_INDEX )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and reflectance variables

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  include files of linearized solution variables (output stored here)

      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'

C  include files of linearized reflectance variables (output stored here)

      INCLUDE '../includes/LIDORT_L_REFLECTANCE.VARS'

C  inputs
C  ------

C  .. direct beam inputs

      LOGICAL          DO_INCLUDE_DIRECTBEAM

C  .. surface emission inputs

      LOGICAL          DO_INCLUDE_SURFEMISS

C  Beam index

      INTEGER          IBEAM_INDEX
 
C  .. BRDF index (indicates which albedo kernel)

      INTEGER          BRDF_INDEX

C  -- parameter index

      INTEGER          PAR_INDEX

C  factor

      DOUBLE PRECISION SURFACE_FACTOR

C  local variables
C  ---------------

      INTEGER          N, I, C0, CM, AA, N1, J, B, IB, Q
      DOUBLE PRECISION HELP, REFL_B, REFL_P, REFL_M
      DOUBLE PRECISION AWF_BEAM, AWF_HOMP, AWF_HOMM
      DOUBLE PRECISION AWF_DIRECT, AWF_EMISS, PARVALUE
      DOUBLE PRECISION FACTOR_BRDF, REFL_ATTN

      DOUBLE PRECISION H_XPOS(MAXSTREAMS,MAXSTREAMS)
      DOUBLE PRECISION H_XNEG(MAXSTREAMS,MAXSTREAMS)
      DOUBLE PRECISION H_WLOWER(MAXSTREAMS)

C  initialise

      DO I = 1, NTOTAL
        COL2_WFALB(I,1) = ZERO
      ENDDO

C  boundary conditions not changed for first layer upper (TOA)

      N = 1
      DO I = 1, NSTREAMS
        COL2_WFALB(I,1) = ZERO
      ENDDO

C  boundary conditions not changed for all intermediate levels

      DO N = 2, NLAYERS - 1
        N1 = N - 1
        C0  = N1*NSTREAMS_2 - NSTREAMS
        DO I = 1, NSTREAMS_2
          CM = C0 + I
          COL2_WFALB(CM,1) = ZERO
        ENDDO
      ENDDO

C  Ground level boundary condition
C  -------------------------------

C  Initialise

      N  = NLAYERS
      C0 = (N-1)*NSTREAMS_2 + NSTREAMS
      Q  = PAR_INDEX
      B  = BRDF_INDEX
      IB = IBEAM_INDEX

      PARVALUE    = BRDF_PARAMETERS(B,Q)
      FACTOR_BRDF = SURFACE_FACTOR * BRDF_FACTORS(B) * PARVALUE

C  Save some quantities
C  ( This is a repetition of earlier code and could be stored)

      DO J = 1, NSTREAMS
        H_WLOWER(J) = AX(J)*WLOWER(J,NLAYERS)
      ENDDO
      DO AA = 1, NSTREAMS
        DO J = 1, NSTREAMS
          H_XPOS(J,AA) = AX(J)*XPOS(J,AA,NLAYERS)
          H_XNEG(J,AA) = AX(J)*XNEG(J,AA,NLAYERS)
        ENDDO
      ENDDO

C  Diffuse scatter contributions

      DO I = 1, NSTREAMS
        CM = C0 + I
        REFL_B = ZERO
        DO J = 1, NSTREAMS
          REFL_B = REFL_B + H_WLOWER(J) * D_BIREFLEC(B,Q,J,I)
        ENDDO
        AWF_BEAM = FACTOR_BRDF * REFL_B
        HELP = AWF_BEAM
        DO AA = 1, NSTREAMS
          REFL_P = ZERO
          REFL_M = ZERO
          DO J = 1, NSTREAMS
            REFL_P = REFL_P + H_XPOS(J,AA) * D_BIREFLEC(B,Q,J,I)
            REFL_M = REFL_M + H_XNEG(J,AA) * D_BIREFLEC(B,Q,J,I)
          ENDDO
          AWF_HOMP = REFL_P * FACTOR_BRDF
          AWF_HOMM = REFL_M * FACTOR_BRDF
          HELP = HELP + LCON(AA,N) * T_DELT_EIGEN(AA,N) * AWF_HOMP
     &                + MCON(AA,N) * AWF_HOMM
        ENDDO
        COL2_WFALB(CM,1) = HELP
      ENDDO

C  Direct beam reflection

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        REFL_ATTN = BRDF_FACTORS(B) * ATMOS_ATTN(IB)
        HELP      = PARVALUE        * REFL_ATTN
        DO I = 1, NSTREAMS
          CM = C0 + I
          AWF_DIRECT = HELP * D_BIREFLEC_0(B,Q,I,IB)
          COL2_WFALB(CM,1) = COL2_WFALB(CM,1) + AWF_DIRECT
        ENDDO
      ENDIF

C  If surface emission, emissivity has a derivative
C  Emissivity derivative is straight-up (not normalized) --> mult by PARVALUE
C  Emissivity derivative has correct sign already worked in
C   Linearization code tested with Thermal supplement, 02 March 2007.

      IF ( DO_INCLUDE_SURFEMISS ) THEN
        DO I = 1, NSTREAMS
          CM = C0 + I
          AWF_EMISS = SURFBB * PARVALUE * D_EMISSIVITY(B,Q,I)
          COL2_WFALB(CM,1) = COL2_WFALB(CM,1) + AWF_EMISS
        ENDDO
      ENDIF

C  Copy for the single layer case

      IF ( NLAYERS .EQ. 1 ) THEN
        DO N = 1, NTOTAL
          SCOL2_WFALB(N,1) = COL2_WFALB(N,1)
        ENDDO
      ENDIF

C  debug

      if ( do_debug_write ) then
        DO N = 1, NTOTAL
          write(85,'(2i4,1p4e17.9)')IBEAM_INDEX,N,COL2_WFALB(N,1)
        ENDDO
      ENDIF
           
C  Finish

      RETURN
      END

c

      SUBROUTINE L_BVPTEL_SOLUTION_MASTER
     I       ( LAYER_TO_VARY,     N_LAYER_WFS,
     I         FOURIER_COMPONENT, IBEAM,
     O         STATUS )

C  include files
C  -------------

C  include file of dimensiopns and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of solution and setup variables (inputs to module)

      INCLUDE '../includes/LIDORT_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  include file of linearized setup variables (input)

      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'

C  include file of linearized solution variables (input and output)

      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'

C  input arguments
C  ---------------

C  Fourier component and beam number

      INTEGER          FOURIER_COMPONENT, IBEAM

C  Linearization control

      INTEGER          LAYER_TO_VARY
      INTEGER          N_LAYER_WFS

C  output status
C  -------------

      INTEGER          STATUS

C  Local variables
C  ---------------

C  error tracing variables

      INTEGER          INFO
      CHARACTER*3      CI
      CHARACTER*70     MAIL, TRACE

C  Other local help variables 

      INTEGER          I, I1, K, K1, Q
      INTEGER          NS, N, N1, NAF, NAL, AA, C0, CM, CP
      DOUBLE PRECISION SHOM, L_HOM1, L_HOM2

C  Initialise

      STATUS = LIDORT_SUCCESS
      Q = 0

C  Set up linearized BVP, column vector  B' where AX = B'
C  ======================================================

C  Bulk: Compute the main column B' where AX = B'

      IF ( DO_COLUMN_LINEARIZATION ) THEN

        CALL L_BVPTEL_C_COLUMN_SETUP
     I     ( N_LAYER_WFS,
     I       FOURIER_COMPONENT,  IBEAM )

C  Profile: Boundary condition flags for special cases
C  Profile: Compute the main column B' where AX = B'

      ELSE IF ( DO_PROFILE_LINEARIZATION ) THEN

        CALL L_BVPTEL_P_COLUMN_SETUP
     I     ( LAYER_TO_VARY,      N_LAYER_WFS,
     I       FOURIER_COMPONENT,  IBEAM )

      ENDIF

C  Solve linearized BVP, several active layers
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
         TRACE = 'DGBTRS call in L_BVPTEL_BEAMSOLUTION_MASTER'
         STATUS = LIDORT_SERIOUS
         CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
         RETURN
        ENDIF

C  Set linearized integration constants, active layers

        C0 = - NSTREAMS_2
        DO NS = 1, NLAYERS_TEL
         N = ACTIVE_LAYERS(NS)
         C0 = C0 + NSTREAMS_2
         DO I = 1, NSTREAMS
          CM = C0 + I
          CP = CM + NSTREAMS
          DO Q = 1, N_LAYER_WFS
           NCON(I,N,Q) = COLTEL2_WF(CM,Q)
           PCON(I,N,Q) = COLTEL2_WF(CP,Q)
          ENDDO
         ENDDO
        ENDDO

C  Solve linearized BVP: Single Layer only
C  =======================================

      ELSE IF ( NLAYERS_TEL .EQ. 1 ) THEN

C  LAPACK substitution (DGETRS) using RHS column vector SCOL2_WF

        CALL DGETRS
     &     ( 'N', NSTREAMS_2, N_LAYER_WFS, 
     &        SMAT2, MAXSTREAMS_2, SIPIVOT,
     &        SCOL2_WF, MAXSTREAMS_2, INFO )

C  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          MAIL='BVP back-substitution (telescoping) (DGETRS)'
          CALL LIDORT_LAPACK_ERROR ( INFO, 1, MAIL, STATUS )
          IF ( STATUS .EQ. LIDORT_SERIOUS ) RETURN
        ENDIF

C  Set linearized integration constants for active layer

        N = ACTIVE_LAYERS(1)
        DO K = 1, NSTREAMS
         K1 = K + NSTREAMS 
         DO Q = 1, N_LAYER_WFS
           NCON(K,N,Q) = SCOL2_WF(K,Q)
           PCON(K,N,Q) = SCOL2_WF(K1,Q)
         ENDDO
        ENDDO

C  end clause for backsubstitution

      ENDIF

C  Associated quantities for active layers
C  ---------------------------------------

      DO NS = 1, NLAYERS_TEL
        N = ACTIVE_LAYERS(NS)
        DO I = 1, NSTREAMS_2
          DO K = 1, NSTREAMS
            DO Q = 1, N_LAYER_WFS
              NCON_XVEC(I,K,N,Q) = NCON(K,N,Q) * XPOS(I,K,N)
              PCON_XVEC(I,K,N,Q) = PCON(K,N,Q) * XNEG(I,K,N)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

C  Set linearized integration constants for non-active layers
C  ==========================================================

C  Now we propagate the results upwards and downwards through the
C  appropriate non-active layers where there is no scattering.

C  Column linearization case is treated separately
C    New code by R. Spurr, 31 October 2007

       IF ( DO_COLUMN_LINEARIZATION ) GO TO 4554

C  Case 1. Profile linearization
C  *****************************

C  Transmittance layers ABOVE active layer(s)
C  -----------------------------------------

C   --NCON values are zero (no downwelling radiation)
C   --PCON values propagated upwards from top of first active layer

C  layer immediately above first active layer
C   --- Require linearized solutions at top of first active layer
C   --- Additional linearizations required if the first active
C       layer is the varying layer

      NAF = ACTIVE_LAYERS(1)
      IF ( NAF .GT. 1 ) THEN
        N1 = NAF - 1
        IF ( LAYER_TO_VARY.EQ.NAF ) THEN
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO Q = 1, N_LAYER_WFS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                L_HOM1 = NCON_XVEC(I1,AA,NAF,Q) + 
     &             LCON(AA,NAF)*L_XPOS(I1,AA,NAF,Q)
                L_HOM2 = T_DELT_EIGEN(AA,NAF) * 
     &             ( PCON_XVEC(I1,AA,NAF,Q)             + 
     &               MCON(AA,NAF)*L_XNEG(I1,AA,NAF,Q) ) +
     &          L_T_DELT_EIGEN(AA,NAF,Q) * MCON_XVEC(I1,AA,NAF)
                SHOM = SHOM + L_HOM1 + L_HOM2
              ENDDO
              PCON(I,N1,Q) = L_WUPPER(I1,NAF,Q) + SHOM
              NCON(I,N1,Q) = ZERO
            ENDDO
          ENDDO
        ELSE IF ( LAYER_TO_VARY.LT.NAF) THEN
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO Q = 1, N_LAYER_WFS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                L_HOM1 = NCON_XVEC(I1,AA,NAF,Q)
                L_HOM2 = T_DELT_EIGEN(AA,NAF) * PCON_XVEC(I1,AA,NAF,Q)
                SHOM = SHOM + L_HOM1 + L_HOM2
              ENDDO
              PCON(I,N1,Q) = L_WUPPER(I1,NAF,Q) + SHOM
              NCON(I,N1,Q) = ZERO
            ENDDO
          ENDDO
        ELSE
          DO I = 1, NSTREAMS
            DO Q = 1, N_LAYER_WFS
              PCON(I,N1,Q) = ZERO
              NCON(I,N1,Q) = ZERO
            ENDDO
          ENDDO
        ENDIF
      ENDIF

C  For remaining non-active atmospheric layers to TOA, propagate upwards.
C   Additional linearizations if you are passing through the varying layer.

      DO N = NAF - 2, 1, -1
        N1 = N + 1
        DO I = 1, NSTREAMS
          DO Q = 1, N_LAYER_WFS
            NCON(I,N,Q) = ZERO
            PCON(I,N,Q) = T_DELT_DISORDS(I,N1) * PCON(I,N1,Q)
          ENDDO
        ENDDO
        IF ( N1 .EQ. LAYER_TO_VARY ) THEN
          DO I = 1, NSTREAMS
            DO Q = 1, N_LAYER_WFS
              PCON(I,N,Q) = PCON(I,N,Q)
     &            + L_T_DELT_DISORDS(I,N1,Q) * MCON(I,N1)
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

      NAL = ACTIVE_LAYERS(NLAYERS_TEL)
      IF ( NAL .LT. NLAYERS ) THEN
        N1 = NAL + 1
        IF ( LAYER_TO_VARY .EQ. NAL ) THEN
          DO I = 1, NSTREAMS
            DO Q = 1, N_LAYER_WFS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                L_HOM2 = PCON_XVEC(I,AA,NAL,Q) + 
     &             MCON(AA,NAL)*L_XNEG(I,AA,NAL,Q)
                L_HOM1 = T_DELT_EIGEN(AA,NAL) * 
     &             ( NCON_XVEC(I,AA,NAL,Q)             + 
     &               LCON(AA,NAL)*L_XPOS(I,AA,NAL,Q) ) +
     &          L_T_DELT_EIGEN(AA,NAL,Q) * LCON_XVEC(I,AA,NAL)
                SHOM = SHOM + L_HOM1 + L_HOM2
              ENDDO
              NCON(I,N1,Q) = L_WLOWER(I,NAL,Q) + SHOM
              PCON(I,N1,Q) = ZERO
            ENDDO
          ENDDO
        ELSE IF ( LAYER_TO_VARY. LT. NAL ) THEN
          DO I = 1, NSTREAMS
            DO Q = 1, N_LAYER_WFS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                L_HOM2 = PCON_XVEC(I,AA,NAL,Q) 
                L_HOM1 = T_DELT_EIGEN(AA,NAL) * NCON_XVEC(I,AA,NAL,Q)
                SHOM = SHOM + L_HOM1 + L_HOM2
              ENDDO
              NCON(I,N1,Q) = L_WLOWER(I,NAL,Q) + SHOM
              PCON(I,N1,Q) = ZERO
            ENDDO
          ENDDO
        ELSE
          DO I = 1, NSTREAMS
            DO Q = 1, N_LAYER_WFS
              PCON(I,N1,Q) = ZERO
              NCON(I,N1,Q) = ZERO
            ENDDO
          ENDDO
        ENDIF
      ENDIF

C  other layers to bottom of medium: propagate downwards.
C   Additional variation if you are passing through the varying layer.

      DO N = NAL + 2, NLAYERS
        N1 = N - 1
        DO I = 1, NSTREAMS
          DO Q = 1, N_LAYER_WFS
            PCON(I,N,Q) = ZERO
            NCON(I,N,Q) = T_DELT_DISORDS(I,N1) * NCON(I,N1,Q)
          ENDDO
        ENDDO
        IF ( N1 .EQ. LAYER_TO_VARY ) THEN
          DO I = 1, NSTREAMS
            DO Q = 1, N_LAYER_WFS
              NCON(I,N,Q) = NCON(I,N,Q)
     &              + L_T_DELT_DISORDS(I,N1,Q) * LCON(I,N1)
            ENDDO
          ENDDO
        ENDIF
      ENDDO

C  Go to final stuff

      GO TO 4555

C  Case 2. Column linearization
C  ****************************

C Continuation point for doing this

 4554 continue

C  Transmittance layers ABOVE active layer(s)
C  -----------------------------------------

C   --NCON values are zero (no downwelling radiation)
C   --PCON values propagated upwards from top of first active layer

C  layer immediately above first active layer
C   --- Require linearized solutions at top of first active layer
C   --- Additional linearizations required, first layer is always active

      NAF = ACTIVE_LAYERS(1)
      IF ( NAF .GT. 1 ) THEN
        N1 = NAF - 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO Q = 1, N_LAYER_WFS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              L_HOM1 = NCON_XVEC(I1,AA,NAF,Q) + 
     &           LCON(AA,NAF)*L_XPOS(I1,AA,NAF,Q)
              L_HOM2 = T_DELT_EIGEN(AA,NAF) * 
     &           ( PCON_XVEC(I1,AA,NAF,Q)             + 
     &             MCON(AA,NAF)*L_XNEG(I1,AA,NAF,Q) ) +
     &        L_T_DELT_EIGEN(AA,NAF,Q) * MCON_XVEC(I1,AA,NAF)
              SHOM = SHOM + L_HOM1 + L_HOM2
            ENDDO
            PCON(I,N1,Q) = L_WUPPER(I1,NAF,Q) + SHOM
            NCON(I,N1,Q) = ZERO
          ENDDO
        ENDDO
      ENDIF

C  For remaining non-active atmospheric layers to TOA, propagate upwards.
C   Additional linearizations if you are passing through the varying layer.

      DO N = NAF - 2, 1, -1
        N1 = N + 1
        DO I = 1, NSTREAMS
          DO Q = 1, N_LAYER_WFS
            NCON(I,N,Q) = ZERO
            PCON(I,N,Q) = T_DELT_DISORDS(I,N1) * PCON(I,N1,Q)
     &            + L_T_DELT_DISORDS(I,N1,Q) * MCON(I,N1)
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

      NAL = ACTIVE_LAYERS(NLAYERS_TEL)
      IF ( NAL .LT. NLAYERS ) THEN
        N1 = NAL + 1
        DO I = 1, NSTREAMS
          DO Q = 1, N_LAYER_WFS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              L_HOM2 = PCON_XVEC(I,AA,NAL,Q) + 
     &           MCON(AA,NAL)*L_XNEG(I,AA,NAL,Q)
              L_HOM1 = T_DELT_EIGEN(AA,NAL) * 
     &           ( NCON_XVEC(I,AA,NAL,Q)             + 
     &             LCON(AA,NAL)*L_XPOS(I,AA,NAL,Q) ) +
     &        L_T_DELT_EIGEN(AA,NAL,Q) * LCON_XVEC(I,AA,NAL)
              SHOM = SHOM + L_HOM1 + L_HOM2
            ENDDO
            NCON(I,N1,Q) = L_WLOWER(I,NAL,Q) + SHOM
            PCON(I,N1,Q) = ZERO
          ENDDO
        ENDDO
      ENDIF

C  other layers to bottom of medium: propagate downwards.
C   Additional variation if you are passing through the varying layer.

      DO N = NAL + 2, NLAYERS
        N1 = N - 1
        DO I = 1, NSTREAMS
          DO Q = 1, N_LAYER_WFS
            PCON(I,N,Q) = ZERO
            NCON(I,N,Q) = T_DELT_DISORDS(I,N1) * NCON(I,N1,Q)
     &              + L_T_DELT_DISORDS(I,N1,Q) * LCON(I,N1)
          ENDDO
        ENDDO
      ENDDO

C  Continuation point for final stuff

 4555 continue

C  Associated quantities for inactive layers
C  -----------------------------------------

C  atmosphere layers with no scattering

      DO N = 1, NLAYERS
        IF ( N .LT. NAF .OR. N.GT.NAL ) THEN
          DO I = 1, NSTREAMS_2
            DO AA = 1, NSTREAMS
              NCON_XVEC(I,AA,N,Q) = NCON(AA,N,Q) * XPOS(I,AA,N)
              PCON_XVEC(I,AA,N,Q) = PCON(AA,N,Q) * XNEG(I,AA,N)
            ENDDO
          ENDDO
        ENDIF
      ENDDO

C  debug

c      if ( layer_to_vary .eq. 11 ) then
c       if ( fourier_component .EQ.3 .and.ibeam.eq.1 ) then
c        do n = 1, nlayers
c          write(99,'(a,2i3,1p2e18.10)')
c     &         'hey',n,layer_to_vary,ncon(3,n,1),pcon(3,n,1)
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

C  finish

      RETURN
      END

C

      SUBROUTINE L_BVPTEL_P_COLUMN_SETUP
     I     ( LAYER_TO_VARY,      N_LAYER_WFS,
     I       FOURIER_COMPONENT,  IBEAM )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup, solution variables

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  include files of linearized setup variables (input)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'

C  include files of linearized solution variables (output stored here)

      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'

C  input arguments
C  ---------------

C  Fourier component and beam number

      INTEGER          FOURIER_COMPONENT, IBEAM

C  Linearization control

      INTEGER          LAYER_TO_VARY
      INTEGER          N_LAYER_WFS

C  local variables
C  ---------------

      INTEGER          Q,AA,N,N1,NS,I,I1,CM,C0,NAF
      DOUBLE PRECISION CPOS, CNEG, L_HOM, L_BEAM

C  Try this safety-first zeroing

      DO I = 1, NSTREAMS_2
        DO Q = 1, N_LAYER_WFS
          L_WUPPER(I,LAYER_TO_VARY,Q) = ZERO
          L_WLOWER(I,LAYER_TO_VARY,Q) = ZERO
        ENDDO
      ENDDO

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
          COLTEL2_WF(I,1) = ZERO
        ENDDO
      ENDDO

C  top of first active layer, first boundary condition
C  ---------------------------------------------------

      NS = 1
      N = ACTIVE_LAYERS(NS)
      NAF = N

C  If this active layer = layer that is varying,
C       then require homogeneous and beam solution linearizations

      IF ( LAYER_TO_VARY .EQ. N ) THEN

        DO I = 1, NSTREAMS
          DO Q = 1, N_LAYER_WFS
            L_BEAM = - L_WUPPER(I,N,Q)
            L_HOM  = ZERO
            DO AA = 1, NSTREAMS
              CPOS = L_XPOS(I,AA,N,Q)
              CNEG = T_DELT_EIGEN(AA,N)   * L_XNEG(I,AA,N,Q) + 
     &             L_T_DELT_EIGEN(AA,N,Q) *   XNEG(I,AA,N)
              L_HOM = L_HOM + LCON(AA,N) * CPOS + MCON(AA,N) * CNEG
            ENDDO
            COLTEL2_WF(I,Q) = L_BEAM - L_HOM
          ENDDO
        ENDDO

C  otherwise if varying layer is above first active layer, there are beam
C  solution contributions propagated downwards - find these by calling
C  the appropriate solution module = L_BEAMSOLUTION_NNEK

      ELSE IF ( LAYER_TO_VARY .LT. N ) THEN

        CALL L_BEAMSOLUTION_P_NNEK
     &    ( FOURIER_COMPONENT, IBEAM, N, LAYER_TO_VARY, N_LAYER_WFS )
        DO I = 1, NSTREAMS
          DO Q = 1, N_LAYER_WFS
             COLTEL2_WF(I,Q) = - L_WUPPER(I,N,Q)
          ENDDO
        ENDDO

      ENDIF

C  Intermediate boundaries between active layers
C  ---------------------------------------------

      DO NS = 1, NLAYERS_TEL - 1

C  offsets

       N  = ACTIVE_LAYERS(NS)
       N1 = N + 1
       C0 = NS*NSTREAMS_2 - NSTREAMS

C  if N is the varying layer, immediately above boundary
C  Get the linearized beam solution for the next layer N1

       IF ( N .EQ. LAYER_TO_VARY ) THEN

        CALL L_BEAMSOLUTION_P_NNEK
     &    ( FOURIER_COMPONENT, IBEAM, N1, LAYER_TO_VARY, N_LAYER_WFS )

        DO I = 1, NSTREAMS_2
          CM = C0 + I
          DO Q = 1, N_LAYER_WFS
            L_BEAM = L_WUPPER(I,N1,Q) - L_WLOWER(I,N,Q)      
            L_HOM  = ZERO
            DO AA = 1, NSTREAMS
              CNEG = L_XNEG(I,AA,N,Q)
              CPOS = T_DELT_EIGEN(AA,N)   * L_XPOS(I,AA,N,Q) + 
     &             L_T_DELT_EIGEN(AA,N,Q) *   XPOS(I,AA,N)
              L_HOM = L_HOM + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
            ENDDO
            COLTEL2_WF(CM,Q) = L_BEAM - L_HOM
          ENDDO
        ENDDO

C  If N1 is the varying layer, immediately below boundary
C    Only require contributions from this layer

       ELSE IF ( N1 .EQ. LAYER_TO_VARY ) THEN

        DO I = 1, NSTREAMS_2
          CM = C0 + I
          DO Q = 1, N_LAYER_WFS
            L_BEAM  = + L_WUPPER(I,N1,Q)
            L_HOM = ZERO
            DO AA = 1, NSTREAMS
              CPOS = L_XPOS(I,AA,N1,Q)
              CNEG = T_DELT_EIGEN(AA,N1)   * L_XNEG(I,AA,N1,Q) + 
     &             L_T_DELT_EIGEN(AA,N1,Q) *   XNEG(I,AA,N1)
              L_HOM = L_HOM + LCON(AA,N1) * CPOS + MCON(AA,N1) * CNEG
            ENDDO
            COLTEL2_WF(CM,Q) = L_BEAM + L_HOM
          ENDDO
        ENDDO

C  non-zero variations if LAYER_TO_VARY is an active layer above N
C    Get the linearized beam solution for the next layer
C  .. contributions from beam solutions on both sides.

       ELSE IF ( LAYER_TO_VARY .LT. N ) THEN

        CALL L_BEAMSOLUTION_P_NNEK
     &    ( FOURIER_COMPONENT, IBEAM, N1, LAYER_TO_VARY, N_LAYER_WFS )

        DO I = 1, NSTREAMS_2
          CM = C0 + I
          DO Q = 1, N_LAYER_WFS
            COLTEL2_WF(CM,Q) = L_WUPPER(I,N1,Q) - L_WLOWER(I,N,Q)
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
      C0 = (NS-1)*NSTREAMS_2 + NSTREAMS

C  If last active layer is varying

      IF ( N .EQ. LAYER_TO_VARY ) THEN

        DO I = 1, NSTREAMS
          CM = C0 + I
          I1 = I + NSTREAMS
          DO Q = 1, N_LAYER_WFS
            L_BEAM = L_WLOWER(I1,N,Q)
            L_HOM  = ZERO
            DO AA = 1, NSTREAMS
              CPOS = T_DELT_EIGEN(AA,N)   * L_XPOS(I1,AA,N,Q) + 
     &             L_T_DELT_EIGEN(AA,N,Q) *   XPOS(I1,AA,N)
              CNEG = L_XNEG(I1,AA,N,Q)
              L_HOM = L_HOM + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
            ENDDO
            COLTEL2_WF(CM,Q) = - L_BEAM - L_HOM
          ENDDO
        ENDDO

C  If varying layer is above layer and is active, beam contribution

      ELSE IF ( LAYER_TO_VARY .LT. N ) THEN

        DO I = 1, NSTREAMS
          CM = C0 + I
          I1 = I + NSTREAMS
          DO Q = 1, N_LAYER_WFS
            L_BEAM = L_WLOWER(I1,N,Q)
            COLTEL2_WF(CM,Q) = - L_BEAM
          ENDDO
        ENDDO

      ENDIF

C  normal return

      RETURN

C  Continuation point for the single-active-layer case

 3456 CONTINUE

C  zero column vector
C   Probably not needed

      DO I = 1, NSTREAMS_2
        DO Q = 1, MAX_ATMOSWFS
          SCOL2_WF(I,Q) = ZERO
        ENDDO
      ENDDO

C  top of active layer
C  -------------------

      NS = 1
      N = ACTIVE_LAYERS(NS)

C  If active layer = layer that is varying,
C       then require homogeneous and beam solution linearizations

      IF ( LAYER_TO_VARY .EQ. N ) THEN

        DO I = 1, NSTREAMS
          DO Q = 1, N_LAYER_WFS
            L_BEAM = - L_WUPPER(I,N,Q)
            L_HOM  = ZERO
            DO AA = 1, NSTREAMS
              CPOS = L_XPOS(I,AA,N,Q)
              CNEG = T_DELT_EIGEN(AA,N)   * L_XNEG(I,AA,N,Q) + 
     &             L_T_DELT_EIGEN(AA,N,Q) *   XNEG(I,AA,N)
              L_HOM = L_HOM + LCON(AA,N) * CPOS + MCON(AA,N) * CNEG
            ENDDO
            SCOL2_WF(I,Q) = L_BEAM - L_HOM
          ENDDO
        ENDDO

C  otherwise if varying layer is above active layer, there are beam
C  solution contributions propagated downwards - find these by calling
C  the appropriate solution module = L_BEAMSOLUTION_NNEK

      ELSE IF ( LAYER_TO_VARY .LT. N ) THEN

        CALL L_BEAMSOLUTION_P_NNEK
     &    ( FOURIER_COMPONENT, IBEAM, N, LAYER_TO_VARY, N_LAYER_WFS )
        DO I = 1, NSTREAMS
          DO Q = 1, N_LAYER_WFS
            SCOL2_WF(I,Q) = - L_WUPPER(I,N,Q)
          ENDDO
        ENDDO

      ENDIF

C  Bottom of active layer
C  ----------------------

C  If active layer is varying layer, need full calculation.

      IF ( N .EQ. LAYER_TO_VARY ) THEN

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO Q = 1, N_LAYER_WFS
            L_BEAM = L_WLOWER(I1,N,Q)
            L_HOM  = ZERO
            DO AA = 1, NSTREAMS
              CPOS = T_DELT_EIGEN(AA,N)   * L_XPOS(I1,AA,N,Q) + 
     &             L_T_DELT_EIGEN(AA,N,Q) *   XPOS(I1,AA,N)
              CNEG = L_XNEG(I1,AA,N,Q)
              L_HOM = L_HOM + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
            ENDDO
            SCOL2_WF(I1,Q) = - L_BEAM - L_HOM
          ENDDO
        ENDDO

C  otherwise use beam solution linearizations propagated downwards
C  from the layer that is varying (already computed for this layer)

      ELSE IF ( LAYER_TO_VARY .LT. N ) THEN

        DO I = 1, NSTREAMS
          DO Q = 1, N_LAYER_WFS
            I1 = I + NSTREAMS
            SCOL2_WF(I1,Q) = - L_WLOWER(I1,N,Q)
          ENDDO
        ENDDO

      ENDIF

C  finish

      RETURN
      END

C

      SUBROUTINE L_BVPTEL_C_COLUMN_SETUP
     I     ( N_LAYER_WFS,
     I       FOURIER_COMPONENT,  IBEAM )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup, solution variables

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SOLUTION.VARS'

C  include files of linearized setup variables (input)

      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'

C  include files of linearized solution variables (output stored here)

      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'

C  input arguments
C  ---------------

C  Fourier component and beam number

      INTEGER          FOURIER_COMPONENT, IBEAM

C  Linearization control

      INTEGER          N_LAYER_WFS

C  local variables
C  ---------------

      INTEGER          Q,AA,N,N1,NS,I,I1,CM,C0,NAF
      DOUBLE PRECISION CPOS, CNEG, L_HOM, L_BEAM

C  Get the linearized solutions for all active layers
C    Always need this, regardless of number of active layers

      DO NS = 1, NLAYERS_TEL
        N = ACTIVE_LAYERS(NS)
        CALL L_BEAMSOLUTION_C_NEQK
     I    ( FOURIER_COMPONENT, IBEAM, N, N_LAYER_WFS )
      ENDDO

C  Go to special case for only 1 active layer

      IF ( NLAYERS_TEL .EQ. 1 ) GO TO 3456

C  zero column vector

      DO I = 1, N_BVTELMATRIX_SIZE
        DO Q = 1, MAX_ATMOSWFS
          COLTEL2_WF(I,1) = ZERO
        ENDDO
      ENDDO

C  top of first active layer, first boundary condition
C  ---------------------------------------------------

      NS = 1
      N = ACTIVE_LAYERS(NS)
      NAF = N

C  require homogeneous and beam solution linearizations

      DO I = 1, NSTREAMS
        DO Q = 1, N_LAYER_WFS
          L_BEAM = - L_WUPPER(I,N,Q)
          L_HOM  = ZERO
           DO AA = 1, NSTREAMS
            CPOS = L_XPOS(I,AA,N,Q)
            CNEG = T_DELT_EIGEN(AA,N)   * L_XNEG(I,AA,N,Q) + 
     &           L_T_DELT_EIGEN(AA,N,Q) *   XNEG(I,AA,N)
            L_HOM = L_HOM + LCON(AA,N) * CPOS + MCON(AA,N) * CNEG
          ENDDO
          COLTEL2_WF(I,Q) = L_BEAM - L_HOM
        ENDDO
      ENDDO

C  Intermediate boundaries between active layers
C  ---------------------------------------------

      DO NS = 1, NLAYERS_TEL - 1

C  offsets

        N  = ACTIVE_LAYERS(NS)
        N1 = N + 1
        C0 = NS*NSTREAMS_2 - NSTREAMS

C  Get the linearized beam solution for the next layer N1

        DO I = 1, NSTREAMS_2
          CM = C0 + I
          DO Q = 1, N_LAYER_WFS
            L_BEAM = L_WUPPER(I,N1,Q) - L_WLOWER(I,N,Q)      
            L_HOM  = ZERO
            DO AA = 1, NSTREAMS
              CNEG = L_XNEG(I,AA,N,Q)
              CPOS = T_DELT_EIGEN(AA,N)   * L_XPOS(I,AA,N,Q) + 
     &             L_T_DELT_EIGEN(AA,N,Q) *   XPOS(I,AA,N)
              L_HOM = L_HOM + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
            ENDDO
            COLTEL2_WF(CM,Q) = L_BEAM - L_HOM
          ENDDO
        ENDDO

C  End loop over intermediate active layer boundaries

      ENDDO

C  Final boundary, bottom of lowest active layer
C  ---------------------------------------------

      NS = NLAYERS_TEL
      N  = ACTIVE_LAYERS(NS)      
      C0 = (NS-1)*NSTREAMS_2 + NSTREAMS

C  last active layer is varying

      DO I = 1, NSTREAMS
        CM = C0 + I
        I1 = I + NSTREAMS
        DO Q = 1, N_LAYER_WFS
          L_BEAM = L_WLOWER(I1,N,Q)
          L_HOM  = ZERO
          DO AA = 1, NSTREAMS
            CPOS = T_DELT_EIGEN(AA,N)   * L_XPOS(I1,AA,N,Q) + 
     &           L_T_DELT_EIGEN(AA,N,Q) *   XPOS(I1,AA,N)
            CNEG = L_XNEG(I1,AA,N,Q)
            L_HOM = L_HOM + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
          ENDDO
          COLTEL2_WF(CM,Q) = - L_BEAM - L_HOM
        ENDDO
      ENDDO

C  Continuation point for the single-active-layer case

 3456 CONTINUE

C  zero column vector
C   Probably not needed

      DO I = 1, NSTREAMS_2
        DO Q = 1, MAX_ATMOSWFS
          SCOL2_WF(I,Q) = ZERO
        ENDDO
      ENDDO

C  top of active layer
C  -------------------

      NS = 1
      N = ACTIVE_LAYERS(NS)

C  layer that is varying,
C    require homogeneous and beam solution linearizations

      DO I = 1, NSTREAMS
        DO Q = 1, N_LAYER_WFS
          L_BEAM = - L_WUPPER(I,N,Q)
          L_HOM  = ZERO
          DO AA = 1, NSTREAMS
           CPOS = L_XPOS(I,AA,N,Q)
            CNEG = T_DELT_EIGEN(AA,N)   * L_XNEG(I,AA,N,Q) + 
     &           L_T_DELT_EIGEN(AA,N,Q) *   XNEG(I,AA,N)
            L_HOM = L_HOM + LCON(AA,N) * CPOS + MCON(AA,N) * CNEG
          ENDDO
          SCOL2_WF(I,Q) = L_BEAM - L_HOM
        ENDDO
      ENDDO


C  Bottom of active layer
C  ----------------------

C  active layer is varying layer

      DO I = 1, NSTREAMS
        I1 = I + NSTREAMS
        DO Q = 1, N_LAYER_WFS
          L_BEAM = L_WLOWER(I1,N,Q)
          L_HOM  = ZERO
          DO AA = 1, NSTREAMS
            CPOS = T_DELT_EIGEN(AA,N)   * L_XPOS(I1,AA,N,Q) + 
     &           L_T_DELT_EIGEN(AA,N,Q) *   XPOS(I1,AA,N)
            CNEG = L_XNEG(I1,AA,N,Q)
            L_HOM = L_HOM + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
          ENDDO
          SCOL2_WF(I1,Q) = - L_BEAM - L_HOM
        ENDDO
      ENDDO

C  finish

      RETURN
      END

