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
C #            TWOSTREAM_BVP_MATRIXSETUP_MASTER   (master)      #
C #            TWOSTREAM_BVP_SOLUTION_MASTER      (master)      #
C #                                                             #
C ###############################################################

      SUBROUTINE TWOSTREAM_BVP_MATRIXSETUP_MASTER
     I     ( DO_INCLUDE_SURFACE, FOURIER_COMPONENT,
     I       SURFACE_FACTOR, ALBEDO, SURFTYPE, NLAYERS, NTOTAL, 
     I       XPOS, XNEG, T_DELT_EIGEN, STREAM_VALUE,
     O       R2_HOMP, R2_HOMM,
     O       BMAT_ROWMASK, BANDMAT2, SMAT2, IPIVOT,
     O       STATUS )

C  input
C  -----

C  control

      INTEGER          SURFTYPE
      LOGICAL          DO_INCLUDE_SURFACE
      INTEGER          FOURIER_COMPONENT
      DOUBLE PRECISION SURFACE_FACTOR
      DOUBLE PRECISION ALBEDO
      INTEGER          NLAYERS, NTOTAL

C  Eigenvector solutions

      DOUBLE PRECISION XPOS(2,NLAYERS)
      DOUBLE PRECISION XNEG(2,NLAYERS)

C  transmittance factors for +/- eigenvalues

      DOUBLE PRECISION T_DELT_EIGEN(NLAYERS)

C  Stream

      DOUBLE PRECISION STREAM_VALUE

C  Output
C  ------

C  Initialization of BVP matrix

      INTEGER          BMAT_ROWMASK(NTOTAL,NTOTAL)

C  Reflected solutions

      DOUBLE PRECISION R2_HOMP
      DOUBLE PRECISION R2_HOMM

C  Matrix, Band-matrix for solving BCs

      DOUBLE PRECISION BANDMAT2(7,NTOTAL)

C  square matrix for the single layer case

      DOUBLE PRECISION SMAT2   (2,2)

C  Pivot matrices

      INTEGER          IPIVOT  (NTOTAL)

C  status
 
      INTEGER          STATUS

C  local variables
C  ---------------

      CHARACTER*(70)   MAIL, TRACE

      INTEGER          J, C0, I, N, N1, INFO
      INTEGER          CEM,CEP,CEM1,CEP1,CM,CP
      INTEGER          NMIN(NTOTAL),NMAX(NTOTAL)

      DOUBLE PRECISION XPNET, XMNET, FACTOR, A, B, C, D, DET
      CHARACTER*3      CI

C  Status

      status = 0

C  Additional setups for the lowest layer
C  For Lambertian reflectance, all streams are the same
C   Placeholder for the other surface types

      R2_HOMP = 0.0d0
      R2_HOMM = 0.0d0
      IF ( DO_INCLUDE_SURFACE ) THEN
        IF ( SURFTYPE .EQ. 1 ) THEN
          FACTOR = SURFACE_FACTOR * ALBEDO
          R2_HOMP = FACTOR * XPOS(1,NLAYERS) * STREAM_VALUE
          R2_HOMM = FACTOR * XNEG(1,NLAYERS) * STREAM_VALUE
        ENDIF
      ENDIF

C  initialize compression matrix (Do this for every Fourier component)

      IF ( NLAYERS .gt. 1 ) THEN
        DO N = 1, NTOTAL
          DO J = 1, 7
            BANDMAT2(J,N) = 0.0d0
          ENDDO
        ENDDO
      ENDIF

C  For row indices (FOURIER = 0 only)

      IF ( FOURIER_COMPONENT .EQ. 0 ) THEN

c  compression row indices

        DO J = 1, 3
          NMIN(J) = 1
        ENDDO
        DO J = 4, NTOTAL
          NMIN(J) = J - 2
        ENDDO
        DO J = 1, NTOTAL - 2
          NMAX(J) = J + 2
        ENDDO
        DO J = NTOTAL - 1, NTOTAL
          NMAX(J) = NTOTAL
        ENDDO

C  Compression algorithm

        DO I = 1, NTOTAL
          DO J = 1, NTOTAL
            IF ( (I.GE.NMIN(J)) .AND. (I.LE.NMAX(J)) ) THEN
              BMAT_ROWMASK(I,J) = 5 + I - J
            ENDIF
          ENDDO
        ENDDO
       
      ENDIF

C  set up boundary values matrix in compressed form (the "A" as in AX=B)
C  ---------------------------------------------------------------------

C  If Nlayers = 1, go to special case

      IF ( NLAYERS .EQ. 1 ) GO TO 345

C  top BC for layer 1: no downward diffuse radiation

      N = 1
      BANDMAT2(BMAT_ROWMASK(1,1),1)  = XPOS(1,N)
      BANDMAT2(BMAT_ROWMASK(1,2),2)  = XNEG(1,N)*T_DELT_EIGEN(N)

C  intermediate layer boundaries (will not be done if NLAYERS = 1 )

      C0 = - 1
      DO N = 2, NLAYERS
        N1 = N - 1
        C0   = C0 + 2
        DO I = 1, 2
          CM = C0 + I
          CEP = C0
          CEM = CEP + 1
          CEP1 = CEP + 2
          CEM1 = CEM + 2
          BANDMAT2(BMAT_ROWMASK(CM,CEP),CEP)   =
     &                   T_DELT_EIGEN(N1)*XPOS(I,N1)
          BANDMAT2(BMAT_ROWMASK(CM,CEM),CEM)   =
     &                   XNEG(I,N1)
          BANDMAT2(BMAT_ROWMASK(CM,CEP1),CEP1) =
     &                   -XPOS(I,N)
          BANDMAT2(BMAT_ROWMASK(CM,CEM1),CEM1) =
     &                   -T_DELT_EIGEN(N)*XNEG(I,N)
        ENDDO
      ENDDO

C  bottom BC (with albedo additions if flagged)

      N = NLAYERS
      C0 = C0 + 2
      CP = C0 + 1
      CEP = C0
      CEM = CEP + 1
      IF ( DO_INCLUDE_SURFACE ) THEN
        XPNET = XPOS(2,N) - R2_HOMP
        XMNET = XNEG(2,N) - R2_HOMM
      ELSE
        XPNET = XPOS(2,N)
        XMNET = XNEG(2,N)
      ENDIF        
      BANDMAT2(BMAT_ROWMASK(CP,CEP),CEP) = T_DELT_EIGEN(N) * XPNET
      BANDMAT2(BMAT_ROWMASK(CP,CEM),CEM) = XMNET

C  normal completion

      GO TO 456

C  special case for 1 layer

345   CONTINUE

C  top BC : no downward diffuse radiation

      N = 1
      SMAT2(1,1) = XPOS(1,N)
      SMAT2(1,2) = XNEG(1,N)*T_DELT_EIGEN(N)

C  bottom BC (with surface reflection

      IF ( DO_INCLUDE_SURFACE ) THEN
        XPNET = XPOS(2,N) - R2_HOMP
        XMNET = XNEG(2,N) - R2_HOMM
      ELSE
        XPNET = XPOS(2,N)
        XMNET = XNEG(2,N)
      ENDIF        
      SMAT2(2,1) = T_DELT_EIGEN(N) * XPNET
      SMAT2(2,2) = XMNET

C  SVD decomposition of compressed boundary values matrix
C  ------------------------------------------------------

 456  continue

      IF ( NLAYERS .GT. 1 ) THEN

C  LAPACK LU-decomposition for band matrix

        CALL DGBTRF
     &     ( NTOTAL, NTOTAL, 2, 2,
     &       BANDMAT2, 7, IPIVOT, INFO )

C  (Error tracing)

        IF ( INFO .GT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MAIL  = 'Singular matrix, u(i,i)=0, for i = '//CI
          TRACE = 'DGBTRF call in BVP_MATRIXSETUP_MASTER'
          STATUS = 1
          CALL TWOSTREAM_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ELSE IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MAIL  = 'argument i illegal value, for i = '//CI
          TRACE = 'DGBTRF call in BVP_MATRIXSETUP_MASTER'
          STATUS = 1
          CALL TWOSTREAM_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ENDIF

C  Analytic solution Single Layer only
C  -----------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN
        A = SMAT2(1,1)
        B = SMAT2(1,2)
        C = SMAT2(2,1)
        D = SMAT2(2,2)
        DET = (A*D - B*C)
        IF ( DABS(DET) .LT. 1.0d-15 ) THEN
          MAIL=' Zero determinant (nlayers=1)'
          STATUS = 1
          CALL TWOSTREAM_ERROR_TRACE ( MAIL, ' ',STATUS )
        ELSE
          DET = 1.0d0 / DET
          SMAT2(1,1) =   D * DET 
          SMAT2(1,2) = - B * DET 
          SMAT2(2,1) = - C * DET 
          SMAT2(2,2) =   A * DET
        ENDIF
      ENDIF

C  finish

      RETURN
      END


C

      SUBROUTINE TWOSTREAM_BVP_SOLUTION_MASTER
     I       ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM,
     I         FOURIER_COMPONENT, IPARTIC, NBEAMS, NLAYERS, NTOTAL,
     I         SURFTYPE, SURFACE_FACTOR, ALBEDO, DIRECT_BEAM,
     I         XPOS, XNEG, WUPPER, WLOWER, STREAM_VALUE,
     I         BANDMAT2, SMAT2, IPIVOT,
     O         R2_PARTIC, COL2, SCOL2, LCON, MCON, 
     O         LCON_XVEC,  MCON_XVEC, 
     O         STATUS )

C  input arguments
C  ---------------

C  inclusion flags

      LOGICAL          DO_INCLUDE_DIRECTBEAM
      INTEGER          SURFTYPE
      LOGICAL          DO_INCLUDE_SURFACE
      DOUBLE PRECISION SURFACE_FACTOR
      DOUBLE PRECISION ALBEDO

      INTEGER          NBEAMS, NLAYERS, NTOTAL

C  Eigenvector solutions

      DOUBLE PRECISION XPOS(2,NLAYERS)
      DOUBLE PRECISION XNEG(2,NLAYERS)

C  Stream

      DOUBLE PRECISION STREAM_VALUE

C  Fourier component and beam number

      INTEGER          FOURIER_COMPONENT, IPARTIC

C  particular solutions

      DOUBLE PRECISION WLOWER ( 2, NLAYERS )
      DOUBLE PRECISION WUPPER ( 2, NLAYERS )

C  Direct beam

      DOUBLE PRECISION DIRECT_BEAM ( NBEAMS )

C  Matrix, Band-matrix

      DOUBLE PRECISION SMAT2   (2,2)
      DOUBLE PRECISION BANDMAT2(7,NTOTAL)

C  Pivot matrices

      INTEGER          IPIVOT  (NTOTAL)

C  output
C  ------

C  reflected solution

      DOUBLE PRECISION R2_PARTIC

C  Column vectors for solving BCs

      DOUBLE PRECISION COL2    (NTOTAL,NBEAMS)
      DOUBLE PRECISION SCOL2   (2,NBEAMS)

C  Solution constants of integration, and related quantities

      DOUBLE PRECISION LCON(NLAYERS)
      DOUBLE PRECISION MCON(NLAYERS)

      DOUBLE PRECISION LCON_XVEC(2,NLAYERS)
      DOUBLE PRECISION MCON_XVEC(2,NLAYERS)

C  status

      INTEGER          STATUS

C  Local variables
C  ---------------

      CHARACTER*(70)   MAIL, TRACE
      INTEGER          N, N1, I, INFO, C0, CM
      DOUBLE PRECISION FACTOR, A, B
      CHARACTER*3      CI

C  Regular BVP using compressed-band matrices, etc..
C  ==================================================

C  --Additional setups for the bottom layer
C  ----------------------------------------

C  Zero total reflected contributions
C  For Lambertian reflectance, all streams are the same

      R2_PARTIC = 0.0d0
      IF ( DO_INCLUDE_SURFACE ) THEN
        FACTOR = SURFACE_FACTOR * ALBEDO
        R2_PARTIC = FACTOR * WLOWER(1,NLAYERS) * STREAM_VALUE
      ENDIF

C  --set up Column for solution vector (the "B" as in AX=B)
C  --------------------------------------------------------

C  If Nlayers = 1, go to special case

      IF ( NLAYERS .EQ. 1 ) GO TO 345

C  zero column vector

      DO I = 1, NTOTAL
        COL2(I,IPARTIC) = 0.0d0
      ENDDO

C  Upper boundary for layer 1: no downward diffuse radiation

      COL2(1,IPARTIC)   = - WUPPER(1,1)

C  intermediate layer boundaries

      DO N = 2, NLAYERS
        N1 = N - 1
        C0 = N1 * 2 - 1
        DO I = 1, 2
          CM = C0 + I
          COL2(CM,IPARTIC) = WUPPER(I,N) - WLOWER(I,N1)
        ENDDO
      ENDDO

C  lowest (surface) boundary with albedo (diffuse radiation terms only)

      N = NLAYERS
      C0 = (N-1)*2 + 1

C  with non-zero surface term, include integrated downward reflectances

      CM = C0 + 1
      IF ( DO_INCLUDE_SURFACE ) THEN
        COL2(CM,IPARTIC) = - WLOWER(2,NLAYERS) + R2_PARTIC
      ELSE
        COL2(CM,IPARTIC) = - WLOWER(2,NLAYERS)
      ENDIF

C  Add direct beam solution (only to final level)

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        IF ( DO_INCLUDE_SURFACE ) THEN
          COL2(CM,IPARTIC) = COL2(CM,IPARTIC) + DIRECT_BEAM(IPARTIC)
        ENDIF
      ENDIF

C  special case - Only one layer

345   CONTINUE

C  zero column vector

      DO I = 1, 2
        SCOL2(I,IPARTIC) = 0.0d0
      ENDDO

C  Upper boundary for layer 1: no downward diffuse radiation

      SCOL2(1,IPARTIC)   = - WUPPER(1,1)

C  lowest (surface) boundary with albedo (diffuse radiation terms only)
C  with non-zero albedo, include integrated downward reflectances
C  no albedo, similar code excluding integrated reflectance

      IF ( DO_INCLUDE_SURFACE ) THEN
        SCOL2(2,IPARTIC) = - WLOWER(2,1) + R2_PARTIC
      ELSE
        SCOL2(2,IPARTIC) = - WLOWER(2,1)
      ENDIF

C  Add direct beam solution (only to final level)

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        IF ( DO_INCLUDE_SURFACE ) THEN
          SCOL2(2,IPARTIC) = SCOL2(2,IPARTIC)+DIRECT_BEAM(IPARTIC)
        ENDIF
      ENDIF

C  --Solve the boundary problem for this Fourier component (back substitution)
C ---------------------------------------------------------------------------

C  BVP back-substitution: With compression (multilayers)

      IF ( NLAYERS .GT. 1 ) THEN

C  LAPACK substitution (DGBTRS) using RHS column vector COL2

        CALL DGBTRS
     &     ( 'n', NTOTAL, 2, 2, IPARTIC,
     &        BANDMAT2, 7, IPIVOT,
     &        COL2, NTOTAL, INFO )

C  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MAIL  = 'argument i illegal value, for i = '//CI
          TRACE = 'DGBTRS call in BVP_SOLUTION_MASTER'
          STATUS = 1
          CALL TWOSTREAM_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
         ENDIF

C  Set integration constants LCON and MCON for -/+ eigensolutions, all layers

        DO N = 1, NLAYERS
         C0 = (N-1)*2
         LCON(N) = COL2(C0+1,IPARTIC)
         MCON(N) = COL2(C0+2,IPARTIC)
        ENDDO

C  Solve the boundary problem: No compression, Single Layer only

      ELSE IF ( NLAYERS .EQ. 1 ) THEN
        A = SCOL2(1,IPARTIC)
        B = SCOL2(2,IPARTIC)
        SCOL2(1,IPARTIC) = SMAT2(1,1) * A + smat2(1,2) * B
        SCOL2(2,IPARTIC) = SMAT2(2,1) * A + smat2(2,2) * B
        LCON(1) = SCOL2(1,IPARTIC)
        MCON(1) = SCOL2(2,IPARTIC)
      ENDIF

C  debug
c      if ( fourier_component.eq.0 ) then
c        do n = 1, nlayers
c          write(*,'(i3,1p2e24.12)')n,LCON(N),MCON(N)
c        enddo
c      endif

C  Associated quantities
C  ---------------------

      DO N = 1, NLAYERS
        DO I = 1, 2
          LCON_XVEC(I,N) = LCON(N)*XPOS(I,N)
          MCON_XVEC(I,N) = MCON(N)*XNEG(I,N)
        ENDDO
      ENDDO

C  Finish

      RETURN
      END


