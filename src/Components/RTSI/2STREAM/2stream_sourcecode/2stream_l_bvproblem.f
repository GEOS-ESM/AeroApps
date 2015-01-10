C ###########################################################
C #                                                         #
C #             THE TWOSTREAM LIDORT MODEL                  #
C #                                                         #
C #      (LInearized Discrete Ordinate Radiative Transfer)  #
C #       --         -        -        -         -          #
C #                                                         #
C ###########################################################

C ###########################################################
C #                                                         #
C #  Authors :      Robert. J. D. Spurr (1)                 #
C #                 Vijay Natraj        (2)                 #
C #                                                         #
C #  Address (1) :     RT Solutions, Inc.                   #
C #                    9 Channing Street                    #
C #                    Cambridge, MA 02138, USA             #
C #  Tel:             (617) 492 1183                        #
C #  Email :           rtsolutions@verizon.net              #
C #                                                         #
C #  Address (2) :     CalTech                              #
C #                    Department of Planetary Sciences     #
C #                    1200 East California Boulevard       #
C #                    Pasadena, CA 91125                   #
C #  Tel:             (626) 395 6962                        #
C #  Email :           vijay@gps.caltech.edu                #
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
C #            TWOSTREAM_BVP_L_SOLUTION_MASTER      (master)    #
C #            TWOSTREAM_L_BVP_P_COLUMN_SETUP                   #
C #            TWOSTREAM_L_BVP_C_COLUMN_SETUP                   #
C #                                                             #
C #            TWOSTREAM_BVP_LS_SOLUTION_MASTER      (master)   #
C #            TWOSTREAM_LS_BVP_COLUMN_SETUP                    #
C #                                                             #
C ###############################################################

      SUBROUTINE TWOSTREAM_BVP_L_SOLUTION_MASTER
     I       ( DO_INCLUDE_DIRECTBEAM, DO_PLANE_PARALLEL,
     I         DO_PROFILE_LINEARIZATION, DO_COLUMN_LINEARIZATION,
     I         FOURIER_COMPONENT, IPARTIC,
     I         NBEAMS, NLAYERS, NTOTAL, NPARS,
     I         VARIATION_INDEX, N_WEIGHTFUNCS,
     I         DO_INCLUDE_SURFACE, SURFTYPE, SURFACE_FACTOR,
     I         ALBEDO, STREAM_VALUE, DIRECT_BEAM, CHAPMAN_FACTORS, 
     I         INITIAL_TRANS, T_DELT_MUBAR, WVEC, WUPPER, WLOWER,
     I         EIGENTRANS, XPOS, XNEG, BANDMAT2, SMAT2, IPIVOT,
     I         LCON, MCON, L_DELTAU_VERT, 
     I         L_INITIAL_TRANS, L_T_DELT_MUBAR,
     I         L_EIGENTRANS, L_XPOS, L_XNEG, L_WVEC,
     O         L_WUPPER, L_WLOWER, COL2_WF, SCOL2_WF,
     O         NCON, PCON, NCON_XVEC, PCON_XVEC, 
     O         STATUS )

C  input arguments
C  ---------------

C  inclusion flags

      LOGICAL          DO_INCLUDE_DIRECTBEAM
      LOGICAL          DO_PLANE_PARALLEL

C  LInearization flags

      LOGICAL          DO_COLUMN_LINEARIZATION
      LOGICAL          DO_PROFILE_LINEARIZATION

C  Surface control

      INTEGER          SURFTYPE
      LOGICAL          DO_INCLUDE_SURFACE
      DOUBLE PRECISION SURFACE_FACTOR
      DOUBLE PRECISION ALBEDO

C  Numbers

      INTEGER          NBEAMS, NLAYERS, NTOTAL, NPARS

C  Stream

      DOUBLE PRECISION STREAM_VALUE

C  Fourier component and beam number

      INTEGER          FOURIER_COMPONENT, IPARTIC

C  Linearization control

      INTEGER          VARIATION_INDEX
      INTEGER          N_WEIGHTFUNCS

C  tramsittance factors for solar beams.

      DOUBLE PRECISION
     &     INITIAL_TRANS  ( NLAYERS, NBEAMS ),
     &     T_DELT_MUBAR   ( NLAYERS, NBEAMS )

C  Eigenvector solutions

      DOUBLE PRECISION EIGENTRANS(NLAYERS)
      DOUBLE PRECISION XPOS(2,NLAYERS)
      DOUBLE PRECISION XNEG(2,NLAYERS)

C  particular solutions

      DOUBLE PRECISION WVEC   ( 2, NLAYERS )
      DOUBLE PRECISION WLOWER ( 2, NLAYERS )
      DOUBLE PRECISION WUPPER ( 2, NLAYERS )

C  Direct beam

      DOUBLE PRECISION DIRECT_BEAM ( NBEAMS )

C  Matrix, Band-matrix

      DOUBLE PRECISION SMAT2   (2,2)
      DOUBLE PRECISION BANDMAT2(7,NTOTAL)

C  Pivot matrices

      INTEGER          IPIVOT  (NTOTAL)

C  Solution constants of integration, and related quantities

      DOUBLE PRECISION LCON(NLAYERS)
      DOUBLE PRECISION MCON(NLAYERS)

C  Linearized tramsittance factors for solar beams.

      DOUBLE PRECISION
     &     L_INITIAL_TRANS  ( NLAYERS, NBEAMS,0:NLAYERS,NPARS ),
     &     L_T_DELT_MUBAR   ( NLAYERS, NBEAMS,0:NLAYERS,NPARS )

C  Linearized Beam solutions

      DOUBLE PRECISION L_WVEC(2,NLAYERS,0:NLAYERS,NPARS)

C  Linearized up and down solutions to the homogeneous RT equations

      DOUBLE PRECISION L_EIGENTRANS(NLAYERS,NPARS)
      DOUBLE PRECISION L_XPOS(2,NLAYERS,NPARS)
      DOUBLE PRECISION L_XNEG(2,NLAYERS,NPARS)

C  Direct beam linearizations

      DOUBLE PRECISION L_DELTAU_VERT   ( NLAYERS, NPARS )
      DOUBLE PRECISION CHAPMAN_FACTORS ( NLAYERS, NLAYERS, NBEAMS )

C  output
C  ------

C  Linearized boundary conditions

      DOUBLE PRECISION L_WLOWER ( 2, NLAYERS, NPARS )
      DOUBLE PRECISION L_WUPPER ( 2, NLAYERS, NPARS )

C  Weighting function column matrices

      DOUBLE PRECISION COL2_WF  ( NTOTAL,NPARS)
      DOUBLE PRECISION SCOL2_WF ( 2,     NPARS)

C  Linearized Solution constants of integration, and related quantities

      DOUBLE PRECISION NCON(NLAYERS,NPARS)
      DOUBLE PRECISION PCON(NLAYERS,NPARS)

      DOUBLE PRECISION NCON_XVEC(2,NLAYERS,NPARS)
      DOUBLE PRECISION PCON_XVEC(2,NLAYERS,NPARS)

C  status

      INTEGER          STATUS

C  Local variables
C  ---------------

      INTEGER          N, I, C0, Q
      DOUBLE PRECISION A, B

C  boundary condition flags

      LOGICAL          MODIFIED_BCL3, MODIFIED_BCL4

C  error tracing variables

      INTEGER          INFO
      CHARACTER*3      CN, CI
      CHARACTER*70     MAIL, TRACE

C  Linearization of the regular BVP case
C  =====================================

C  set status

      STATUS = 0

C  Set up the column vectors for Bulk/profile linearizations
C  ---------------------------------------------------------

C  Bulk: Compute the main column B' where AX = B'

      IF ( DO_COLUMN_LINEARIZATION ) THEN

        CALL TWOSTREAM_L_BVP_C_COLUMN_SETUP
     I       ( DO_INCLUDE_DIRECTBEAM, DO_PLANE_PARALLEL,
     I         FOURIER_COMPONENT, IPARTIC, N_WEIGHTFUNCS,
     I         NBEAMS, NLAYERS, NTOTAL, NPARS,
     I         DO_INCLUDE_SURFACE, SURFTYPE, SURFACE_FACTOR,
     I         ALBEDO, STREAM_VALUE, DIRECT_BEAM,
     I         INITIAL_TRANS, T_DELT_MUBAR, WVEC, WUPPER, WLOWER,
     I         EIGENTRANS, XPOS, XNEG, 
     I         BANDMAT2, SMAT2, IPIVOT, LCON, MCON,
     I         L_INITIAL_TRANS, L_T_DELT_MUBAR, L_WVEC,
     I         L_EIGENTRANS, L_XPOS, L_XNEG, 
     I         L_DELTAU_VERT, CHAPMAN_FACTORS,
     O         L_WUPPER, L_WLOWER, COL2_WF, SCOL2_WF )

C  Profile: Boundary condition flags for special cases
C  Profile: Compute the main column B' where AX = B'

      ELSE IF ( DO_PROFILE_LINEARIZATION ) THEN

C  Boundary condition flags for special cases

        MODIFIED_BCL3 = ( VARIATION_INDEX .EQ. 1 )
        MODIFIED_BCL4 = ( VARIATION_INDEX .EQ. NLAYERS )

        CALL TWOSTREAM_L_BVP_P_COLUMN_SETUP
     I       ( DO_INCLUDE_DIRECTBEAM, DO_PLANE_PARALLEL,
     I         MODIFIED_BCL3,      MODIFIED_BCL4,
     I         VARIATION_INDEX,    N_WEIGHTFUNCS,
     I         FOURIER_COMPONENT, IPARTIC, 
     I         NBEAMS, NLAYERS, NTOTAL, NPARS,
     I         DO_INCLUDE_SURFACE, SURFTYPE, SURFACE_FACTOR,
     I         ALBEDO, STREAM_VALUE, DIRECT_BEAM,
     I         INITIAL_TRANS, T_DELT_MUBAR, WVEC, WUPPER, WLOWER,
     I         EIGENTRANS, XPOS, XNEG, 
     I         BANDMAT2, SMAT2, IPIVOT, LCON, MCON,
     I         L_INITIAL_TRANS, L_T_DELT_MUBAR, L_WVEC,
     I         L_EIGENTRANS, L_XPOS, L_XNEG, 
     I         L_DELTAU_VERT, CHAPMAN_FACTORS,
     O         L_WUPPER, L_WLOWER, COL2_WF, SCOL2_WF )

      ENDIF

C  BVP back-substitution: With compression (multilayers)
C  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

C  LAPACK substitution (DGBTRS) using RHS column vector COL2_WF
C  BV solution for perturbed integration constants
C    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS
     &     ( 'n', NTOTAL, 2, 2, N_WEIGHTFUNCS,
     &        BANDMAT2, 7, IPIVOT,
     &        COL2_WF, NTOTAL, INFO )

C  Error tracing

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(CN, '(I3)' ) VARIATION_INDEX
          MAIL  = 'argument i illegal value, for i = '//CI
          TRACE = 'Atmos_Wfs for layer '//CN//
     *         '; (multilayer) DGBTRS call in L_BVP_SOLUTION_MASTER'
          STATUS = 1
          CALL TWOSTREAM_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
         ENDIF

C  Set integration constants NCON and PCON for +/- eigensolutions

        DO N = 1, NLAYERS
          C0 = (N-1)*2
          DO Q = 1, N_WEIGHTFUNCS
            NCON(N,Q) = COL2_WF(C0+1,Q)
            PCON(N,Q) = COL2_WF(C0+2,Q)
          ENDDO
        ENDDO

C  debug
c      if ( fourier_component.eq.0 ) then
c        do n = 1, nlayers
c          write(*,'(i3,1p2e24.12)')n,NCON(N,1),PCON(N,1)
c        enddo
c      endif

C  Solve the boundary problem: No compression, Single Layer only
C  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN
        DO Q = 1, N_WEIGHTFUNCS
          A = SCOL2_WF(1,Q)
          B = SCOL2_WF(2,Q)
          SCOL2_WF(1,Q) = SMAT2(1,1) * A + smat2(1,2) * B
          SCOL2_WF(2,Q) = SMAT2(2,1) * A + smat2(2,2) * B
          NCON(1,Q) = SCOL2_WF(1,Q)
          PCON(1,Q) = SCOL2_WF(2,Q)
        ENDDO
      ENDIF

C  Associated quantities
C  ---------------------

      DO N = 1, NLAYERS
        DO Q = 1, N_WEIGHTFUNCS
          DO I = 1, 2
            NCON_XVEC(I,N,Q) = NCON(N,Q)*XPOS(I,N)
            PCON_XVEC(I,N,Q) = PCON(N,Q)*XNEG(I,N)
          ENDDO
        ENDDO
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE TWOSTREAM_L_BVP_P_COLUMN_SETUP
     I       ( DO_INCLUDE_DIRECTBEAM, DO_PLANE_PARALLEL,
     I         MODIFIED_BCL3,      MODIFIED_BCL4,
     I         LAYER_TO_VARY,      N_LAYER_WFS,
     I         FOURIER_COMPONENT,  IPARTIC, 
     I         NBEAMS, NLAYERS, NTOTAL, NPARS,
     I         DO_INCLUDE_SURFACE, SURFTYPE, SURFACE_FACTOR,
     I         ALBEDO, STREAM_VALUE, DIRECT_BEAM,
     I         INITIAL_TRANS, T_DELT_MUBAR, WVEC, WUPPER, WLOWER,
     I         EIGENTRANS, XPOS, XNEG, 
     I         BANDMAT2, SMAT2, IPIVOT, LCON, MCON,
     I         L_INITIAL_TRANS, L_T_DELT_MUBAR, L_WVEC,
     I         L_EIGENTRANS, L_XPOS, L_XNEG, 
     I         L_DELTAU_VERT, CHAPMAN_FACTORS,
     O         L_WUPPER, L_WLOWER, COL2_WF, SCOL2_WF )

C  input arguments
C  ---------------

C  inclusion flags

      LOGICAL          DO_INCLUDE_DIRECTBEAM
      LOGICAL          DO_PLANE_PARALLEL

C  boundary condition flags

      LOGICAL          MODIFIED_BCL3, MODIFIED_BCL4

C  Fourier component and beam number

      INTEGER          FOURIER_COMPONENT, IPARTIC

C  Linearization control

      INTEGER          LAYER_TO_VARY
      INTEGER          N_LAYER_WFS

C  Surface control, surface factor = 1+delta(m,0)

      INTEGER          SURFTYPE
      LOGICAL          DO_INCLUDE_SURFACE
      DOUBLE PRECISION SURFACE_FACTOR
      DOUBLE PRECISION ALBEDO

C  Numbers

      INTEGER          NBEAMS, NLAYERS, NTOTAL, NPARS

C  Stream

      DOUBLE PRECISION STREAM_VALUE

C  tramsittance factors for solar beams.

      DOUBLE PRECISION
     &     INITIAL_TRANS  ( NLAYERS, NBEAMS ),
     &     T_DELT_MUBAR   ( NLAYERS, NBEAMS )

C  Eigenvector solutions

      DOUBLE PRECISION EIGENTRANS(NLAYERS)
      DOUBLE PRECISION XPOS(2,NLAYERS)
      DOUBLE PRECISION XNEG(2,NLAYERS)

C  particular solutions

      DOUBLE PRECISION WVEC   ( 2, NLAYERS )
      DOUBLE PRECISION WLOWER ( 2, NLAYERS )
      DOUBLE PRECISION WUPPER ( 2, NLAYERS )

C  Direct beam

      DOUBLE PRECISION DIRECT_BEAM ( NBEAMS )

C  Matrix, Band-matrix

      DOUBLE PRECISION SMAT2   (2,2)
      DOUBLE PRECISION BANDMAT2(7,NTOTAL)

C  Pivot matrices

      INTEGER          IPIVOT  (NTOTAL)

C  Solution constants of integration, and related quantities

      DOUBLE PRECISION LCON(NLAYERS)
      DOUBLE PRECISION MCON(NLAYERS)

C  Linearized tramsittance factors for solar beams.

      DOUBLE PRECISION
     &     L_INITIAL_TRANS  ( NLAYERS, NBEAMS,0:NLAYERS,NPARS ),
     &     L_T_DELT_MUBAR   ( NLAYERS, NBEAMS,0:NLAYERS,NPARS )

C  Linearized Beam solutions

      DOUBLE PRECISION L_WVEC(2,NLAYERS,0:NLAYERS,NPARS)

C  Linearized up and down solutions to the homogeneous RT equations

      DOUBLE PRECISION L_EIGENTRANS(NLAYERS,NPARS)
      DOUBLE PRECISION L_XPOS(2,NLAYERS,NPARS)
      DOUBLE PRECISION L_XNEG(2,NLAYERS,NPARS)

C  Direct beam linearizations

      DOUBLE PRECISION L_DELTAU_VERT ( NLAYERS, NPARS )
      DOUBLE PRECISION CHAPMAN_FACTORS ( NLAYERS, NLAYERS, NBEAMS )

C  Outputs
C  -------

C  Linearized Beam solutions at boundaries

      DOUBLE PRECISION L_WUPPER(2,NLAYERS,NPARS)
      DOUBLE PRECISION L_WLOWER(2,NLAYERS,NPARS)

C  Weighting function column matrices

      DOUBLE PRECISION COL2_WF  ( NTOTAL,NPARS)
      DOUBLE PRECISION SCOL2_WF ( 2,     NPARS)

C  local variables
C  ---------------

      INTEGER          Q, N, N1, I, CM, C0, K
      DOUBLE PRECISION CPOS, CNEG, L_HOM, L_BEAM, FAC, L_PAR
      DOUBLE PRECISION HSP_U, HSM_U, CONST, WDEL
      DOUBLE PRECISION VAR1, VAR2, VAR_U, TRANS2, FACTOR
      LOGICAL          REGULAR_BCL3, REGULAR_BCL4

C  initialise
C  ----------

C  zero the results vectors

      DO I = 1, NTOTAL
        DO Q = 1, NPARS
          COL2_WF(I,Q) = 0.0d0
        ENDDO
      ENDDO

C  Layer to vary

      K = LAYER_TO_VARY

C    This is a very important zeroing.................!!!!!

      DO I = 1, 2
        DO Q = 1, N_LAYER_WFS
          DO N = 1, NLAYERS
            L_WUPPER(I,N,Q) = 0.0d0
            L_WLOWER(I,N,Q) = 0.0d0
          ENDDO
        ENDDO
      ENDDO

C  Get the linearized beam solution for the varying layer

      N = K
      CONST   = INITIAL_TRANS(N,IPARTIC)
      WDEL    = T_DELT_MUBAR(N,IPARTIC)
      DO Q = 1, N_LAYER_WFS
        VAR1 = L_T_DELT_MUBAR(N,IPARTIC,N,Q) * CONST
        DO I = 1, 2
          L_WUPPER(I,N,Q) = CONST * L_WVEC(I,N,N,Q)
          L_WLOWER(I,N,Q) = WDEL  * L_WUPPER(I,N,Q) + VAR1*WVEC(I,N)
        ENDDO
      ENDDO
      
C  complete boundary condition flags

      REGULAR_BCL3 = .NOT.MODIFIED_BCL3
      REGULAR_BCL4 = .NOT.MODIFIED_BCL4

C  BCL1 or BCL3M - top of first layer (TOA), UPPER boundary condition
C  ------------------------------------------------------------------

      N = 1

C    If this layer is the one that is varied, use MODIFIED_BCL3 (BCL3M)
C       .. contribution L_PARTIC from beam  solution variations
C       .. contribution L_HOM    from eigen solution variations
C    Otherwise, entry is zero

      IF ( MODIFIED_BCL3 ) THEN
       DO Q = 1, N_LAYER_WFS
         L_PAR = - L_WUPPER(1,N,Q)
         CPOS  = L_XPOS(1,N,Q)
         CNEG  = EIGENTRANS(N)   * L_XNEG(1,N,Q) + 
     &         L_EIGENTRANS(N,Q) *   XNEG(1,N)
         L_HOM = LCON(N) * CPOS + MCON(N) * CNEG
         COL2_WF(1,Q) = L_PAR - L_HOM
        ENDDO
      ELSE
        DO Q = 1, N_LAYER_WFS
          COL2_WF(1,Q) = 0.0d0
        ENDDO
      ENDIF

C  BCL2 Intermediate levels between top layer and varying layer
C  ------------------------------------------------------------

C  [not required if top layer is varying, case MODIFIED_BCL3 above]

      IF ( REGULAR_BCL3 ) THEN
        DO N = 2, LAYER_TO_VARY - 1
          N1 = N - 1
          C0  = 2*N1 - 1
          DO I = 1, 2
            CM = C0 + I
            DO Q = 1, N_LAYER_WFS
              COL2_WF(CM,Q) = 0.0d0
            ENDDO
          ENDDO
        ENDDO
      ENDIF

C  BCL3 - regular upper boundary condition for layer that is varying
C  -----------------------------------------------------------------

      IF ( REGULAR_BCL3 ) THEN
        N = LAYER_TO_VARY
        N1  = N - 1
        C0  = N1*2 - 1
        DO I = 1, 2
          CM = C0 + I
          DO Q = 1, N_LAYER_WFS
            L_PAR = + L_WUPPER(I,N,Q)
            CPOS  = L_XPOS(I,N,Q)
            CNEG  = EIGENTRANS(N)   * L_XNEG(I,N,Q) + 
     &            L_EIGENTRANS(N,Q) *   XNEG(I,N)
            L_HOM = LCON(N) * CPOS + MCON(N) * CNEG
            COL2_WF(CM,Q) = L_PAR + L_HOM
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
        C0 = N*2 - 1

C  Get the linearized beam solution for the next layer
C  Distinguish two cases
C  ..(a) quasi-spherical for n > 1 (only gets done for this case anyway)
C  ..(b) plane-parallel and qs for n = 1
C     Consistent logarithmic linearization of INITIAL_TRANS

        CONST   = INITIAL_TRANS(N1,IPARTIC)
        WDEL    = T_DELT_MUBAR(N1,IPARTIC)
        TRANS2  = CONST * WDEL
        IF ( .NOT. DO_PLANE_PARALLEL  ) THEN
          DO Q = 1, N_LAYER_WFS
            VAR1 = L_T_DELT_MUBAR(N1,IPARTIC,K,Q) * CONST
            VAR2 = L_INITIAL_TRANS(N1,IPARTIC,K,Q)
            DO I = 1, 2
              VAR_U = VAR2 * WVEC(I,N1) + L_WVEC(I,N1,K,Q)
              L_WUPPER(I,N1,Q) = CONST  * VAR_U
              L_WLOWER(I,N1,Q) = TRANS2 * VAR_U + VAR1 * WVEC(I,N1)
            ENDDO
          ENDDO
        ELSE
          DO Q = 1, N_LAYER_WFS
            VAR1 = L_INITIAL_TRANS(N1,IPARTIC,K,Q) * CONST
            DO I = 1, 2
              L_WUPPER(I,N1,Q) = VAR1 * WVEC(I,N1)
              L_WLOWER(I,N1,Q) = L_WUPPER(I,N1,Q) * WDEL
            ENDDO
          ENDDO
        ENDIF

C  .. 2 contributions to WVAR from beam solution variations BEAM_V and BEAM_U 
C  .. contribution HVAR homogeneous (eigenvalue) solution variations

        DO I = 1, 2
          CM = C0 + I
          DO Q = 1, N_LAYER_WFS
            L_PAR = L_WUPPER(I,N1,Q) - L_WLOWER(I,N,Q)
            CNEG  = L_XNEG(I,N,Q)
            CPOS  = EIGENTRANS(N)   * L_XPOS(I,N,Q) + 
     &            L_EIGENTRANS(N,Q) *   XPOS(I,N)
            L_HOM = LCON(N)*CPOS + MCON(N)*CNEG
            COL2_WF(CM,Q) = L_PAR - L_HOM
          ENDDO
        ENDDO

C  End Modified BCL4

      ENDIF

C  BCL5 - Intermediate boundary conditions between varying layer & final layer
C  ---------------------------------------------------------------------------

      IF ( REGULAR_BCL4 ) THEN

        DO N = LAYER_TO_VARY + 1, NLAYERS - 1

          N1 = N + 1
          C0  = N*2 - 1

C  Get the linearized beam solution for the next layer
C  Distinguish two cases
C  ..(a) quasi-spherical for n > 1 (only gets done for this case anyway)
C  ..(b) plane-parallel and qs for n1 = 1
C     Consistent logarithmic linearization of INITIAL_TRANS

          CONST   = INITIAL_TRANS(N1,IPARTIC)
          WDEL    = T_DELT_MUBAR(N1,IPARTIC)
          TRANS2  = CONST * WDEL
          IF ( .NOT. DO_PLANE_PARALLEL  ) THEN
            DO Q = 1, N_LAYER_WFS
              VAR1 = L_T_DELT_MUBAR(N1,IPARTIC,K,Q) * CONST
              VAR2 = L_INITIAL_TRANS(N1,IPARTIC,K,Q)
              DO I = 1, 2
                VAR_U = VAR2 * WVEC(I,N1) + L_WVEC(I,N1,K,Q)
                L_WUPPER(I,N1,Q) = CONST  * VAR_U
                L_WLOWER(I,N1,Q) = TRANS2 * VAR_U + VAR1 * WVEC(I,N1)
              ENDDO
            ENDDO
          ELSE
            DO Q = 1, N_LAYER_WFS
              VAR1 = L_INITIAL_TRANS(N1,IPARTIC,K,Q) * CONST
              DO I = 1, 2
                L_WUPPER(I,N1,Q) = VAR1 * WVEC(I,N1)
                L_WLOWER(I,N1,Q) = L_WUPPER(I,N1,Q) * WDEL
              ENDDO
            ENDDO
          ENDIF

C  .. contributions from beam solution (direct assign). No homog. variation

          DO I = 1, 2
            CM = C0 + I
            DO Q = 1, N_LAYER_WFS
              COL2_WF(CM,Q) = L_WUPPER(I,N1,Q) - L_WLOWER(I,N,Q)
            ENDDO
          ENDDO

C  end layer loop

        ENDDO

C  end BCL5 boundary conditions

      ENDIF

C  Final layer - use BCL6 or BCL4M (last layer is varying)
C  -------------------------------------------------------

      N = NLAYERS
      C0 = (N-1)*2 + 1
      CM = C0 + 1
      FACTOR = SURFACE_FACTOR * ALBEDO * STREAM_VALUE

C  Modified BCL4M Component loop

      IF ( MODIFIED_BCL4 ) THEN
       IF ( DO_INCLUDE_SURFACE ) THEN
        DO Q = 1, N_LAYER_WFS
          HSP_U = L_XPOS(1,N,Q) *   EIGENTRANS(N) +
     &              XPOS(1,N)   * L_EIGENTRANS(N,Q)
          HSM_U = L_XNEG(1,N,Q)
          L_PAR = L_WLOWER(2,N,Q) - L_WLOWER(1,N,Q)* FACTOR
          CPOS  =   EIGENTRANS(N)   * L_XPOS(2,N,Q) + 
     &            L_EIGENTRANS(N,Q) *   XPOS(2,N)
          CPOS  = CPOS          - HSP_U * FACTOR
          CNEG  = L_XNEG(2,N,Q) - HSM_U * FACTOR
          L_HOM = LCON(N)*CPOS + MCON(N)*CNEG
          COL2_WF(CM,Q) = - L_PAR - L_HOM
        ENDDO
       ELSE
        DO Q = 1, N_LAYER_WFS
          L_PAR = L_WLOWER(2,N,Q)
          CPOS  =   EIGENTRANS(N)   * L_XPOS(2,N,Q) + 
     &            L_EIGENTRANS(N,Q) *   XPOS(2,N)
          CNEG  = L_XNEG(2,N,Q)
          L_HOM = LCON(N)*CPOS + MCON(N)*CNEG
          COL2_WF(CM,Q) = - L_PAR - L_HOM
        ENDDO
       ENDIF
      ENDIF

C  ordinary BCL6 Component loop
 
      IF ( .NOT. MODIFIED_BCL4 ) THEN
       IF ( DO_INCLUDE_SURFACE ) THEN
        DO Q = 1, N_LAYER_WFS
          L_PAR = L_WLOWER(2,N,Q) - L_WLOWER(1,N,Q)* FACTOR
          COL2_WF(CM,Q) = - L_PAR
        ENDDO
       ELSE
        DO Q = 1, N_LAYER_WFS
          COL2_WF(CM,Q) = - L_WLOWER(2,N,Q)
        ENDDO
       ENDIF
      ENDIF

C  Add direct beam variation to Final boundary
C  -------------------------------------------

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        CM = C0 + 1
        FAC = - DIRECT_BEAM(IPARTIC) * 
     &              CHAPMAN_FACTORS(N,LAYER_TO_VARY,IPARTIC)
        DO Q = 1, N_LAYER_WFS
          L_BEAM = L_DELTAU_VERT(LAYER_TO_VARY,Q) * FAC
          COL2_WF(CM,Q) = COL2_WF(CM,Q) + L_BEAM
        ENDDO
      ENDIF

C  Copy for the one-layer case

      IF ( NLAYERS .EQ. 1 ) THEN
        DO I = 1, 1
          DO Q = 1, N_LAYER_WFS
            SCOL2_WF(I,Q) = COL2_WF(I,Q)
          ENDDO
        ENDDO
      ENDIF

C  finish

      RETURN
      END

C

      SUBROUTINE TWOSTREAM_L_BVP_C_COLUMN_SETUP
     I       ( DO_INCLUDE_DIRECTBEAM, DO_PLANE_PARALLEL,
     I         FOURIER_COMPONENT,  IPARTIC, N_WEIGHTFUNCS,
     I         NBEAMS, NLAYERS, NTOTAL, NPARS,
     I         DO_INCLUDE_SURFACE, SURFTYPE, SURFACE_FACTOR,
     I         ALBEDO, STREAM_VALUE, DIRECT_BEAM,
     I         INITIAL_TRANS, T_DELT_MUBAR, WVEC, WUPPER, WLOWER,
     I         EIGENTRANS, XPOS, XNEG, 
     I         BANDMAT2, SMAT2, IPIVOT, LCON, MCON,
     I         L_INITIAL_TRANS, L_T_DELT_MUBAR, L_WVEC,
     I         L_EIGENTRANS, L_XPOS, L_XNEG, 
     I         L_DELTAU_VERT, CHAPMAN_FACTORS,
     O         L_WUPPER, L_WLOWER, COL2_WF, SCOL2_WF )

C  input arguments
C  ---------------

C  inclusion flags

      LOGICAL          DO_INCLUDE_DIRECTBEAM
      LOGICAL          DO_PLANE_PARALLEL

C  Fourier component and beam number

      INTEGER          FOURIER_COMPONENT, IPARTIC

C  Linearization control

      INTEGER          N_WEIGHTFUNCS

C  Surface control

      INTEGER          SURFTYPE
      LOGICAL          DO_INCLUDE_SURFACE
      DOUBLE PRECISION SURFACE_FACTOR
      DOUBLE PRECISION ALBEDO

C  Numbers

      INTEGER          NBEAMS, NLAYERS, NTOTAL, NPARS

C  Stream

      DOUBLE PRECISION STREAM_VALUE

C  tramsittance factors for solar beams.

      DOUBLE PRECISION
     &     INITIAL_TRANS  ( NLAYERS, NBEAMS ),
     &     T_DELT_MUBAR   ( NLAYERS, NBEAMS )

C  Eigenvector solutions

      DOUBLE PRECISION EIGENTRANS(NLAYERS)
      DOUBLE PRECISION XPOS(2,NLAYERS)
      DOUBLE PRECISION XNEG(2,NLAYERS)

C  particular solutions

      DOUBLE PRECISION WVEC   ( 2, NLAYERS )
      DOUBLE PRECISION WLOWER ( 2, NLAYERS )
      DOUBLE PRECISION WUPPER ( 2, NLAYERS )

C  Direct beam

      DOUBLE PRECISION DIRECT_BEAM ( NBEAMS )

C  Matrix, Band-matrix

      DOUBLE PRECISION SMAT2   (2,2)
      DOUBLE PRECISION BANDMAT2(7,NTOTAL)

C  Pivot matrices

      INTEGER          IPIVOT  (NTOTAL)

C  Solution constants of integration, and related quantities

      DOUBLE PRECISION LCON(NLAYERS)
      DOUBLE PRECISION MCON(NLAYERS)

C  Linearized tramsittance factors for solar beams.

      DOUBLE PRECISION
     &     L_INITIAL_TRANS  ( NLAYERS, NBEAMS,0:NLAYERS,NPARS ),
     &     L_T_DELT_MUBAR   ( NLAYERS, NBEAMS,0:NLAYERS,NPARS )

C  Linearized Beam solutions

      DOUBLE PRECISION L_WVEC(2,NLAYERS,0:NLAYERS,NPARS)

C  Linearized up and down solutions to the homogeneous RT equations

      DOUBLE PRECISION L_EIGENTRANS(NLAYERS,NPARS)
      DOUBLE PRECISION L_XPOS(2,NLAYERS,NPARS)
      DOUBLE PRECISION L_XNEG(2,NLAYERS,NPARS)

C  Direct beam linearizations

c      DOUBLE PRECISION DELTAU_VERT   ( NLAYERS )
      DOUBLE PRECISION L_DELTAU_VERT ( NLAYERS, NPARS )
      DOUBLE PRECISION CHAPMAN_FACTORS  ( NLAYERS, NLAYERS, NBEAMS )

C  Outputs
C  -------

C  Linearized Beam solutions at boundaries

      DOUBLE PRECISION L_WUPPER(2,NLAYERS,NPARS)
      DOUBLE PRECISION L_WLOWER(2,NLAYERS,NPARS)

C  Weighting function column matrices

      DOUBLE PRECISION COL2_WF  ( NTOTAL,NPARS)
      DOUBLE PRECISION SCOL2_WF ( 2,     NPARS)

C  local variables
C  ---------------

      INTEGER          Q, N, N1, I, CM, C0, K
      DOUBLE PRECISION CPOS, CNEG, L_HOM, L_BEAM, FAC, FAC3
      DOUBLE PRECISION HSP_U, HSM_U, CONST, WDEL, FACTOR, L_PAR
      DOUBLE PRECISION VAR1, VAR2, VAR_U, TRANS2, L_HOMU, L_HOMD

C  initialise
C  ----------

C  zero the results vectors

      DO I = 1, NTOTAL
        DO Q = 1, N_WEIGHTFUNCS
          COL2_WF(I,Q) = 0.0d0
        ENDDO
      ENDDO

C  Very important zeroing

      DO K = 1, NLAYERS
        DO I = 1, 2
          DO Q = 1, N_WEIGHTFUNCS
            L_WUPPER(I,K,Q) = 0.0d0
            L_WLOWER(I,K,Q) = 0.0d0
          ENDDO
        ENDDO
      ENDDO
   
C  Top of first layer (TOA), UPPER boundary condition
C  --------------------------------------------------

      N = 1

C  Get the linearized beam solution for the first layer

      CONST   = INITIAL_TRANS(N,IPARTIC)
      WDEL    = T_DELT_MUBAR(N,IPARTIC)
      TRANS2  = CONST * WDEL
      DO Q = 1, N_WEIGHTFUNCS
        VAR1 = L_T_DELT_MUBAR(N,IPARTIC,0,Q) * CONST
        VAR2 = L_INITIAL_TRANS(N,IPARTIC,0,Q)
        DO I = 1, 2
          VAR_U = VAR2 * WVEC(I,N) + L_WVEC(I,N,0,Q)
          L_WUPPER(I,N,Q) = CONST  * VAR_U
          L_WLOWER(I,N,Q) = TRANS2 * VAR_U + VAR1 * WVEC(I,N)
        ENDDO
      ENDDO

C  COmpute the column contributions

      DO Q = 1, N_WEIGHTFUNCS
        L_PAR = - L_WUPPER(1,N,Q)
        CPOS  = L_XPOS(1,N,Q)
        CNEG  =   EIGENTRANS(N)   * L_XNEG(1,N,Q) + 
     &          L_EIGENTRANS(N,Q) *   XNEG(1,N)
        L_HOM = LCON(N) * CPOS + MCON(N) * CNEG
        COL2_WF(1,Q) = L_PAR - L_HOM
      ENDDO

C  Intermediate boundary conditions
C  --------------------------------

      DO N = 1, NLAYERS - 1

C  N1 is the layer below, C0 is the offset

        N1 = N + 1
        C0 = N*2 - 1

C  Get the linearized beam solution for the next layer

        CONST   = INITIAL_TRANS(N1,IPARTIC)
        WDEL    = T_DELT_MUBAR(N1,IPARTIC)
        TRANS2  = CONST * WDEL
        DO Q = 1, N_WEIGHTFUNCS
          VAR1 = L_T_DELT_MUBAR(N1,IPARTIC,0,Q) * CONST
          VAR2 = L_INITIAL_TRANS(N1,IPARTIC,0,Q)
          DO I = 1, 2
            VAR_U = VAR2 * WVEC(I,N1) + L_WVEC(I,N1,0,Q)
            L_WUPPER(I,N1,Q) = CONST  * VAR_U
            L_WLOWER(I,N1,Q) = TRANS2 * VAR_U + VAR1 * WVEC(I,N1)
          ENDDO
        ENDDO

C  .. 2 contributions to L_BEAM, from variations L_WUPPER L_WLOWER 
C  .. 2 contributions to L_HOM,  from variations above and below

        DO I = 1, 2
          CM = C0 + I
          DO Q = 1, N_WEIGHTFUNCS
            CPOS = L_XPOS(I,N1,Q)
            CNEG =   EIGENTRANS(N1)   * L_XNEG(I,N1,Q) + 
     &             L_EIGENTRANS(N1,Q) *   XNEG(I,N1)
            L_HOMU = LCON(N1) * CPOS + MCON(N1) * CNEG
            CNEG = L_XNEG(I,N,Q)
            CPOS =   EIGENTRANS(N)   * L_XPOS(I,N,Q) + 
     &             L_EIGENTRANS(N,Q) *   XPOS(I,N)
            L_HOMD = LCON(N)*CPOS + MCON(N)*CNEG
            L_PAR  = L_WUPPER(I,N1,Q) - L_WLOWER(I,N,Q)      
            L_HOM  = L_HOMU - L_HOMD
            COL2_WF(CM,Q) = L_PAR + L_HOM
          ENDDO
        ENDDO

C  End layer

      ENDDO

C  LOWEST layer 
C  ------------

      N = NLAYERS
      C0 = (N-1)*2 + 1
      CM = C0 + 1
      FACTOR = SURFACE_FACTOR * ALBEDO * STREAM_VALUE

      IF ( DO_INCLUDE_SURFACE ) THEN
        DO Q = 1, N_WEIGHTFUNCS
          HSP_U = L_XPOS(1,N,Q) *   EIGENTRANS(N) +
     &              XPOS(1,N)   * L_EIGENTRANS(N,Q)
          HSM_U = L_XNEG(1,N,Q)
          L_PAR = L_WLOWER(2,N,Q) - L_WLOWER(1,N,Q) * FACTOR
          CPOS  =   EIGENTRANS(N)   * L_XPOS(2,N,Q) + 
     &            L_EIGENTRANS(N,Q) *   XPOS(2,N)
          CPOS  = CPOS          - HSP_U * FACTOR
          CNEG  = L_XNEG(2,N,Q) - HSM_U * FACTOR
          L_HOM = LCON(N)*CPOS + MCON(N)*CNEG
          COL2_WF(CM,Q) = - L_PAR - L_HOM
        ENDDO
      ELSE
        DO Q = 1, N_WEIGHTFUNCS
          L_PAR = L_WLOWER(2,N,Q)
          CPOS  =   EIGENTRANS(N)   * L_XPOS(2,N,Q) + 
     &            L_EIGENTRANS(N,Q) *   XPOS(2,N)
          CNEG  = L_XNEG(2,N,Q)
          L_HOM = LCON(N)*CPOS + MCON(N)*CNEG
          COL2_WF(CM,Q) = - L_PAR - L_HOM
        ENDDO
      ENDIF

C  Add direct beam variation to Final boundary
C  -------------------------------------------

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        FAC = - DIRECT_BEAM(IPARTIC)
        DO Q = 1, N_WEIGHTFUNCS
          L_BEAM = 0.0d0
          DO K = 1, NLAYERS
            FAC3 = FAC * CHAPMAN_FACTORS(N,K,IPARTIC)
            L_BEAM = L_BEAM + L_DELTAU_VERT(K,Q) * FAC3
          ENDDO
          COL2_WF(CM,Q) = COL2_WF(CM,Q) + L_BEAM
        ENDDO
      ENDIF

C  Copy for the one-layer case

      IF ( NLAYERS .EQ. 1 ) THEN
        DO I = 1, 2
          DO Q = 1, N_WEIGHTFUNCS
            SCOL2_WF(I,Q) = COL2_WF(I,Q)
          ENDDO
        ENDDO
      ENDIF

C  finish

      RETURN
      END

C

      SUBROUTINE TWOSTREAM_BVP_LS_SOLUTION_MASTER
     I       ( DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFACE,
     I         FOURIER_COMPONENT, IBEAM, SURFTYPE,
     I         NBEAMS, NLAYERS, NTOTAL, DIRECT_BEAM,
     I         BANDMAT2, SMAT2, IPIVOT,
     I         LCON, MCON, EIGENTRANS, 
     I         R2_HOMP, R2_HOMM, R2_PARTIC,
     O         COL2_WFALB, SCOL2_WFALB,
     O         NCONALB, PCONALB,
     O         STATUS )

C  WARNING. Only valid for

C  input arguments
C  ---------------

C  inclusion flags

      LOGICAL          DO_INCLUDE_DIRECTBEAM

C  Surface control

      INTEGER          SURFTYPE
      LOGICAL          DO_INCLUDE_SURFACE

C  Numbers

      INTEGER          NBEAMS, NLAYERS, NTOTAL

C  Fourier component and beam number

      INTEGER          FOURIER_COMPONENT, IBEAM

C  ginesolution transmittances

      DOUBLE PRECISION EIGENTRANS(NLAYERS)

C  Diffuse solution at surface (stream value)

      DOUBLE PRECISION R2_HOMP, R2_HOMM, R2_PARTIC

C  Solution constants of integration, and related quantities

      DOUBLE PRECISION LCON(NLAYERS)
      DOUBLE PRECISION MCON(NLAYERS)

C  Reflected Direct beam at surface

      DOUBLE PRECISION DIRECT_BEAM ( NBEAMS )

C  Matrix, Band-matrix

      DOUBLE PRECISION SMAT2   (2,2)
      DOUBLE PRECISION BANDMAT2(7,NTOTAL)

C  Pivot matrices

      INTEGER          IPIVOT  (NTOTAL)

C  output
C  ------

C  Weighting function column matrices

      DOUBLE PRECISION COL2_WFALB  ( NTOTAL,1)
      DOUBLE PRECISION SCOL2_WFALB ( 2,     1)

C  Linearized Solution constants of integration, and related quantities

      DOUBLE PRECISION NCONALB(NLAYERS)
      DOUBLE PRECISION PCONALB(NLAYERS)

C  status

      INTEGER          STATUS

C  Local variables
C  ---------------

      INTEGER          N, I, C0, CM
      DOUBLE PRECISION A, B, IDOWNSURF_Q

C  error tracing variables

      INTEGER          INFO
      CHARACTER*3      CI
      CHARACTER*70     MAIL, TRACE

C  Linearization of the regular BVP case
C  =====================================

C  set status

      STATUS = 0

C  Exit if Fourier component is not 0 and Surface type  not 1
C    This is the Lambertian condition

      IF ( .NOT. DO_INCLUDE_SURFACE .OR.
     &    (FOURIER_COMPONENT.NE.0 .and. SURFTYPE.EQ.1) ) RETURN

C  Set up the column vectors for Surface linearizations
C  ----------------------------------------------------

C  initialise; Only contribution is from lowest layer
C  boundary conditions not changed for all higher levels

      DO I = 1, NTOTAL
        COL2_WFALB(I,1) = 0.0d0
      ENDDO

      N  = NLAYERS
      CM = (N-1) * 2 + 2
 
      IDOWNSURF_Q = R2_PARTIC + MCON(N) * R2_HOMM + 
     &                          LCON(N) * R2_HOMP * EIGENTRANS(N)
      COL2_WFALB(CM,1) = IDOWNSURF_Q
      IF ( DO_INCLUDE_DIRECTBEAM ) COL2_WFALB(CM,1) =
     &       COL2_WFALB(CM,1) + DIRECT_BEAM(IBEAM)

C  Copy for the single layer case

      IF ( NLAYERS .EQ. 1 ) THEN
        DO N = 1, NTOTAL
          SCOL2_WFALB(N,1) = COL2_WFALB(N,1)
        ENDDO
      ENDIF

C  BVP back-substitution: With compression (multilayers)
C  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

C  LAPACK substitution (DGBTRS) using RHS column vector COL2_WF
C  BV solution for perturbed integration constants
C    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS
     &     ( 'n', NTOTAL, 2, 2, 1,
     &        BANDMAT2, 7, IPIVOT,
     &        COL2_WFALB, NTOTAL, INFO )

C  Error tracing

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MAIL  = 'argument i illegal value, for i = '//CI
          TRACE = 'Surface_Wfs for Lambertian weighting function '//
     *         '; (multilayer) DGBTRS call in LS_BVP_SOLUTION_MASTER'
          STATUS = 1
          CALL TWOSTREAM_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
         ENDIF

C  Set integration constants NCON and PCON for +/- eigensolutions

        DO N = 1, NLAYERS
          C0 = (N-1)*2
          NCONALB(N) = COL2_WFALB(C0+1,1)
          PCONALB(N) = COL2_WFALB(C0+2,1)
        ENDDO

C  Solve the boundary problem: No compression, Single Layer only
C  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN
        A = SCOL2_WFALB(1,1)
        B = SCOL2_WFALB(2,1)
        SCOL2_WFALB(1,1) = SMAT2(1,1) * A + smat2(1,2) * B
        SCOL2_WFALB(2,1) = SMAT2(2,1) * A + smat2(2,2) * B
        NCONALB(1) = SCOL2_WFALB(1,1)
        PCONALB(1) = SCOL2_WFALB(2,1)
      ENDIF

C  debug
c      do n = 1, nlayers
c        write(*,*)n,NCONALB(n),pconalb(n)
c      enddo

C  Finish

      RETURN
      END

