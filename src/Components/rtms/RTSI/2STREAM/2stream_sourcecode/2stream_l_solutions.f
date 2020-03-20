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
C # Subroutines in this Module                                  #
C #                                                             #
C #     Homogeneous solution                                    #
C #                                                             #
C #              TWOSTREAM_L_HOM_SOLUTION                       #
C #              TWOSTREAM_L_HOM_USERSOLUTION                   #
C #              TWOSTREAM_L_HMULT_MASTER                       #
C #                                                             #
C #     Particular integrals                                    #
C #                                                             #
C #              TWOSTREAM_L_BEAM_SOLUTION                      #
C #              TWOSTREAM_L_BEAM_USERSOLUTION                  #
C #                                                             #
C ###############################################################

      SUBROUTINE TWOSTREAM_L_HOM_SOLUTION
     I    ( NLAYERS, NPARS, N, FOURIER, DOVARY, NVARY,
     I      STREAM_VALUE, PXSQ, OMEGA, ASYMM, DELTAU_VERT, 
     I      L_OMEGA, L_ASYMM, L_DELTAU_VERT, 
     I      SAB, DAB, EIGENVALUE, EIGENTRANS,
     O      L_EIGENVALUE, L_EIGENTRANS,
     O      L_SAB, L_DAB, L_XPOS, L_XNEG )

C  subroutine input arguments
C  --------------------------

C  Numbers

      INTEGER          NLAYERS
      INTEGER          NPARS

C  Given layer index and Fourier number (inputs)

      INTEGER          N
      INTEGER          FOURIER

C  Linearization control

      LOGICAL          DOVARY
      INTEGER          NVARY

C  Stream value

      DOUBLE PRECISION STREAM_VALUE

C  Legendre input

      DOUBLE PRECISION PXSQ

C  OMEGA and ASYMM and linearizations

      DOUBLE PRECISION OMEGA   ( NLAYERS )
      DOUBLE PRECISION ASYMM   ( NLAYERS )
      DOUBLE PRECISION L_OMEGA ( NLAYERS, NPARS  )
      DOUBLE PRECISION L_ASYMM ( NLAYERS, NPARS  )

C  optical thickness and its linearizations

      DOUBLE PRECISION DELTAU_VERT   ( NLAYERS )
      DOUBLE PRECISION L_DELTAU_VERT ( NLAYERS, NPARS )

C  local matrices for eigenvalue computation

      DOUBLE PRECISION SAB(NLAYERS), DAB(NLAYERS)

C  Eigensolutions

      DOUBLE PRECISION EIGENVALUE(NLAYERS)
      DOUBLE PRECISION EIGENTRANS(NLAYERS)

C  Output arguments 
C  ----------------

C  Eigensolutions

      DOUBLE PRECISION L_EIGENVALUE(NLAYERS,NPARS)
      DOUBLE PRECISION L_EIGENTRANS(NLAYERS,NPARS)

C  Linearized up and down solutions to the homogeneous RT equations

      DOUBLE PRECISION L_XPOS(2,NLAYERS,NPARS)
      DOUBLE PRECISION L_XNEG(2,NLAYERS,NPARS)

C  Saved Linearized sum and difference terms

      DOUBLE PRECISION L_SAB(NLAYERS,NPARS), L_DAB(NLAYERS,NPARS)

C  Local variables
C  ---------------

      INTEGER          Q
      DOUBLE PRECISION XINV, DIFVEC, LSD, LARG
      DOUBLE PRECISION L_OMEGA_ASYMM_3, L_DIFVEC, L_DP, L_DM

C  Initial values

      XINV = 1.0d0 / STREAM_VALUE
      DIFVEC = - SAB(N) / EIGENVALUE(N)

C  start parameter loop

      DO Q = 1, NPARS

C  Develop Linearization of Sum and Difference matrices

        L_OMEGA_ASYMM_3 = 3.0d0 * 
     &     ( L_OMEGA(N,Q)*ASYMM(N) + OMEGA(N)*L_ASYMM(N,Q) )
        if ( fourier.eq.0) then
          L_DP = L_OMEGA(N,Q) + PXSQ *  L_OMEGA_ASYMM_3
          L_DM = L_OMEGA(N,Q) - PXSQ *  L_OMEGA_ASYMM_3
        Else if ( fourier .eq. 1 ) then
          L_DP = L_OMEGA_ASYMM_3 * PXSQ
          L_DM = L_OMEGA_ASYMM_3 * PXSQ
        ENDIF
        L_DAB(N,Q) = ( L_DP - L_DM ) * 0.5d0 * XINV
        L_SAB(N,Q) = ( L_DP + L_DM ) * 0.5d0 * XINV
        LSD = L_DAB(N,Q) * SAB(N) + L_SAB(N,Q) * DAB(N)    

C   Use definitions to find linearizations of eigenproblem

        L_EIGENVALUE(N,Q) = 0.5d0 * LSD / EIGENVALUE(N)
        LARG =  L_EIGENVALUE(N,Q) *   DELTAU_VERT(N) +
     &            EIGENVALUE(N)   * L_DELTAU_VERT(N,Q)
        L_EIGENTRANS(N,Q) = - LARG * EIGENTRANS(N)

C  Auxiliary equation to get up and down solutions
C   Develop linearized solutions from auxiliary Eqn. definition

        L_DIFVEC = -(DIFVEC*L_EIGENVALUE(N,Q)+L_SAB(N,Q))/EIGENVALUE(N) 
        L_XPOS(1,N,Q) = 0.5d0 * L_DIFVEC
        L_XPOS(2,N,Q) = - L_XPOS(1,N,Q) 

C  Symmetry

        L_XNEG(1,N,Q) = L_XPOS(2,N,Q)
        L_XNEG(2,N,Q) = L_XPOS(1,N,Q)

C  debug
c        if (fourier.eq.1)
c     &  write(56,'(2i4,1p3e24.12)')n,q,
c     &         L_EIGENTRANS(N,Q),EIGENTRANS(N)

C finish parameter loop

      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE TWOSTREAM_L_BEAM_SOLUTION
     I    ( NLAYERS, NBEAMS, NPARS, N, FOURIER, IBEAM,
     I      DO_PLANE_PARALLEL, FLUX_FACTOR, 
     I      LAYER_PIS_CUTOFF, STREAM_VALUE, X0, PX0X,
     I      DO_COLUMN_LINEARIZATION, DO_PROFILE_LINEARIZATION,
     I      LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_COLUMN_WFS,
     I      OMEGA, ASYMM, DELTAU_VERT,
     I      L_OMEGA, L_ASYMM, L_DELTAU_VERT,
     I      SAB, DAB, EIGENVALUE, AVERAGE_SECANT, 
     I      QSUMVEC, QDIFVEC, QVEC,
     I      L_AVERAGE_SECANT, L_SAB, L_DAB, L_EIGENVALUE, 
     O      L_WVEC )

C  subroutine arguments
C  --------------------

C  Numbers

      INTEGER          NLAYERS, NBEAMS, NPARS

C  Given layer index and Fourier number, Beam number (inputs)

      INTEGER          N
      INTEGER          FOURIER
      INTEGER          IBEAM

C  Plane Parallel flag

      LOGICAL          DO_PLANE_PARALLEL

C  Flux factor

      DOUBLE PRECISION FLUX_FACTOR

C  Last layer to include Particular integral solution

      INTEGER          LAYER_PIS_CUTOFF(NBEAMS)

C  Stream value

      DOUBLE PRECISION STREAM_VALUE

C  Beam SZA cosines

      DOUBLE PRECISION X0  (NBEAMS)
      DOUBLE PRECISION PX0X(NBEAMS)

C  LInearization flags

      LOGICAL          DO_COLUMN_LINEARIZATION
      LOGICAL          DO_PROFILE_LINEARIZATION

C  Linearization control

      LOGICAL          LAYER_VARY_FLAG(NLAYERS)
      INTEGER          LAYER_VARY_NUMBER(NLAYERS)
      INTEGER          N_COLUMN_WFS

C  Average-secants

      DOUBLE PRECISION AVERAGE_SECANT ( NLAYERS, NBEAMS )

C  OMEGA and ASYMM

      DOUBLE PRECISION OMEGA ( NLAYERS )
      DOUBLE PRECISION ASYMM ( NLAYERS )

C  optical thickness

      DOUBLE PRECISION DELTAU_VERT(NLAYERS)

C  local matrices from eigenvalue computation

      DOUBLE PRECISION SAB(NLAYERS), DAB(NLAYERS)

C  Eigenvalues

      DOUBLE PRECISION EIGENVALUE(NLAYERS)

C  Auxiliary vectors

      DOUBLE PRECISION QSUMVEC(NLAYERS)
      DOUBLE PRECISION QDIFVEC(NLAYERS)
      DOUBLE PRECISION QVEC   (NLAYERS)

C  Linearizations

      DOUBLE PRECISION L_OMEGA       ( NLAYERS, NPARS )
      DOUBLE PRECISION L_ASYMM       ( NLAYERS, NPARS )
      DOUBLE PRECISION L_DELTAU_VERT ( NLAYERS, NPARS )
      DOUBLE PRECISION L_SAB         ( NLAYERS, NPARS )
      DOUBLE PRECISION L_DAB         ( NLAYERS, NPARS )
      DOUBLE PRECISION L_EIGENVALUE  ( NLAYERS, NPARS )

      DOUBLE PRECISION L_AVERAGE_SECANT
     &          ( NLAYERS, NBEAMS, 0:NLAYERS, NPARS )

C  Output variables
C  ----------------

C  Linearized Beam solution

      DOUBLE PRECISION L_WVEC(2,NLAYERS,0:NLAYERS,NPARS)

C  help variables
C  --------------

      LOGICAL          DO_FIRST, DO_PSVAR
      INTEGER          NV, NP, K, Q
      DOUBLE PRECISION TP, TM, SECBAR, XINV, F1
      DOUBLE PRECISION QMAT, QFIN, QAUX, PI4
      DOUBLE PRECISION L_SECBAR, L_INV_X0SQUARE
      DOUBLE PRECISION L_OMEGA_ASYMM, L_QMAT, L_EIGEN_SQUARE
      DOUBLE PRECISION L_QSUMVEC, L_QDIFVEC, L_HELP, L_QAUX, L_QVEC

C  Flux factor and other constants

      PI4 = 4.0d0 * dacos(-1.0d0)
      F1  = FLUX_FACTOR / PI4
      NV  = 0
      XINV = 1.0d0 / STREAM_VALUE

C  Only a solution if the layer is not below Cutoff.

      DO_FIRST = ( N .LE. LAYER_PIS_CUTOFF(IBEAM) )

C  Pseudo-spherical variation condition

      DO_PSVAR = ( .NOT. DO_PLANE_PARALLEL .AND. N.GT.1 )

C  solar zenith cosine for this layer

      SECBAR   = AVERAGE_SECANT(N,IBEAM)

C  set up driving vector (using saved results)
C  This must be done regardless of whether layer N is varying or not.

      QMAT = EIGENVALUE(N) * EIGENVALUE(N) - SECBAR * SECBAR
      QFIN = - SAB(N) * QVEC(N)
      QAUX = ( QFIN - QSUMVEC(N) ) / SECBAR

C  Linearization for layer N is in two parts:
C    1A. Linearization due to variations in Layer N itself (Profiles)
C    1B. Linearization due to columns    
C    2. Linearization due to variations in layers K < N (Profiles)

C  Part 1.
C  =======

C  Set layer to vary

      IF ( DO_PROFILE_LINEARIZATION ) THEN
        NV = N
        NP = LAYER_VARY_NUMBER(N)
      ELSE IF ( DO_COLUMN_LINEARIZATION ) THEN
        NV = 0
        NP = N_COLUMN_WFS
      ENDIF

C  Skip this section if there is nothing varying
C    [ Zero the boundary layer values and start Part 2 ]

      IF ( .NOT. LAYER_VARY_FLAG(N) .OR. .NOT. DO_FIRST) THEN
        DO Q = 1, NP
          L_WVEC(1,N,NV,Q) = 0.0d0
          L_WVEC(2,N,NV,Q) = 0.0d0
        ENDDO
        GO TO 2222
      ENDIF

C  For each varying parameter

      DO Q = 1, NP

C  linearizations

        L_SECBAR       = L_AVERAGE_SECANT(N,IBEAM,NV,Q)
        L_OMEGA_ASYMM  = 3.0d0 * 
     &       ( L_OMEGA(N,Q)*ASYMM(N) + OMEGA(N)*L_ASYMM(N,Q) )
        L_EIGEN_SQUARE = 2.0d0 * L_EIGENVALUE(N,Q) * EIGENVALUE(N)
        L_INV_X0SQUARE = 2.0d0 * L_SECBAR * SECBAR

C  Set up sum and differences for Beam source terms

        if ( fourier.eq.0) then
          TP = L_OMEGA(N,Q) + PX0X(IBEAM) * L_OMEGA_ASYMM
          TM = L_OMEGA(N,Q) - PX0X(IBEAM) * L_OMEGA_ASYMM
        Else if ( fourier .eq. 1 ) then
          TP = PX0X(IBEAM) * L_OMEGA_ASYMM
          TM = PX0X(IBEAM) * L_OMEGA_ASYMM
        ENDIF
        L_QSUMVEC = F1 * ( TP + TM ) * XINV
        L_QDIFVEC = F1 * ( TP - TM ) * XINV

C  Linearize the reduced problem

        L_HELP = - L_DAB(N,Q) * QSUMVEC(N) - DAB(N) * L_QSUMVEC
        L_QMAT = L_EIGEN_SQUARE - L_INV_X0SQUARE
        L_QVEC = L_HELP + L_QDIFVEC * SECBAR 
        IF ( DO_PSVAR ) THEN
          L_QVEC = L_QVEC + QDIFVEC(N) * L_SECBAR
        ENDIF
        L_QVEC = ( L_QVEC - QVEC(N) * L_QMAT ) / QMAT

C  Restore up and down solutions

        L_HELP = - L_SAB(N,Q) * QVEC(N) - SAB(N) * L_QVEC
        L_QAUX = ( L_HELP - L_QSUMVEC ) / SECBAR
        IF ( DO_PSVAR ) THEN
          L_QAUX = L_QAUX - QAUX * L_SECBAR / SECBAR
        ENDIF
        L_WVEC(1,N,NV,Q) = 0.5d0 * ( L_QVEC + L_QAUX )
        L_WVEC(2,N,NV,Q) = 0.5d0 * ( L_QVEC - L_QAUX )

C  End parameter loop

      ENDDO

C  debug
c        if (fourier.eq.0)write(*,'(2i4,1p3e24.12)')n,nv,
c     &         L_WVEC(1,N,NV,1),L_WVEC(1,N,NV,2)

C  Part 2.
C  =======

C  Continuation point

 2222 CONTINUE

C  Only for the pseudo-spherical case, profile weighting functions

C  Also not required for the column weighting functions

      IF ( DO_PLANE_PARALLEL )       RETURN
      IF ( DO_COLUMN_LINEARIZATION ) RETURN

C  Loop over layers K above N

      DO K = 1, N - 1

C  If there is a varying layer

        IF ( LAYER_VARY_FLAG(K) ) THEN

C  Start loop over parameters

          DO Q = 1, LAYER_VARY_NUMBER(K)

C  linearizations

            L_SECBAR       = L_AVERAGE_SECANT(N,IBEAM,K,Q)
            L_INV_X0SQUARE = 2.0d0 * L_SECBAR * SECBAR

C  Linearize the reduced problem

            L_QMAT = - L_INV_X0SQUARE
            L_QVEC = QDIFVEC(N) * L_SECBAR
            L_QVEC = ( L_QVEC - QVEC(N) * L_QMAT ) / QMAT

C  Restore up and down solutions

            L_QAUX =  - ( SAB(N) * L_QVEC + QAUX * L_SECBAR ) / SECBAR
            L_WVEC(1,N,K,Q) = 0.5d0 * ( L_QVEC + L_QAUX )
            L_WVEC(2,N,K,Q) = 0.5d0 * ( L_QVEC - L_QAUX )
   
C  End loops

          ENDDO
        ENDIF
      ENDDO

C  debug
c        if (fourier.eq.0)then
c       write(*,'(2i4,1p3e24.12)')n,1,(L_WVEC(2,N,K,1),k=1,N)
c       write(*,'(2i4,1p3e24.12)')n,2,(L_WVEC(2,N,K,2),k=1,N)
c        endif
c        if (n.eq.3)pause

C  Finish

      RETURN
      END

C

      SUBROUTINE TWOSTREAM_L_HMULT_MASTER
     I     ( NLAYERS, N_USER_STREAMS, NPARS,
     I       DO_PROFILE_WFS, DO_COLUMN_WFS, 
     I       LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_COLUMN_WFS,
     I       EIGENVALUE, EIGENTRANS, USER_STREAMS, T_DELT_USERM,  
     I       ZETA_M, ZETA_P, HMULT_1, HMULT_2,
     I       L_EIGENVALUE, L_EIGENTRANS, L_T_DELT_USERM,  
     O       L_HMULT_1, L_HMULT_2 )

C  Input arguments
C  ===============

C  Numbers

      INTEGER          NLAYERS, N_USER_STREAMS, NPARS

C  Flags

      LOGICAL          DO_PROFILE_WFS, DO_COLUMN_WFS

C  Linearization control

      LOGICAL          LAYER_VARY_FLAG   ( NLAYERS )
      INTEGER          LAYER_VARY_NUMBER ( NLAYERS )
      INTEGER          N_COLUMN_WFS

C  User streams

      DOUBLE PRECISION USER_STREAMS ( N_USER_STREAMS )

C  Transmittance factors for user-defined stream angles

      DOUBLE PRECISION
     &     T_DELT_USERM ( NLAYERS, N_USER_STREAMS )

C  Eigensolutions

      DOUBLE PRECISION EIGENVALUE(NLAYERS)
      DOUBLE PRECISION EIGENTRANS(NLAYERS)

C  coefficient functions for user-defined angles

      DOUBLE PRECISION 
     &      ZETA_M(N_USER_STREAMS,NLAYERS),
     &      ZETA_P(N_USER_STREAMS,NLAYERS)

C  Integrated homogeneous solution multipliers, whole layer

      DOUBLE PRECISION 
     &      HMULT_1(N_USER_STREAMS,NLAYERS),
     &      HMULT_2(N_USER_STREAMS,NLAYERS)

C  Linearized Transmittance factors for user-defined stream angles

      DOUBLE PRECISION
     &     L_T_DELT_USERM ( NLAYERS, N_USER_STREAMS, NPARS )

C  Linearized Eigensolutions

      DOUBLE PRECISION L_EIGENVALUE(NLAYERS,NPARS)
      DOUBLE PRECISION L_EIGENTRANS(NLAYERS,NPARS)

C  Output
C  ======

C  LInearized homogeneous solution multipliers

      DOUBLE PRECISION 
     &      L_HMULT_1(N_USER_STREAMS,NLAYERS,NPARS),
     &      L_HMULT_2(N_USER_STREAMS,NLAYERS,NPARS)

C  Local variables
C  ---------------

      INTEGER          Q, UM, N
      DOUBLE PRECISION SM, L_T_1, L_T_2

C  whole layer multipliers
C  -----------------------

C  Start loops over layers and user-streams
C    Only done if layers are flagged

      IF ( DO_PROFILE_WFS ) THEN
       DO UM = 1, N_USER_STREAMS
        SM   = 1.0d0 / USER_STREAMS(UM)
        DO N = 1, NLAYERS
         IF ( LAYER_VARY_FLAG(N) ) THEN
          DO Q = 1, LAYER_VARY_NUMBER(N)
           L_T_2 =   EIGENTRANS(N)   * L_T_DELT_USERM(N,UM,Q) +
     &             L_EIGENTRANS(N,Q) *   T_DELT_USERM(N,UM)
           L_T_1 = L_EIGENTRANS(N,Q) - L_T_DELT_USERM(N,UM,Q)
           L_HMULT_1(UM,N,Q) = ZETA_M(UM,N) * 
     &           ( SM*L_T_1 + L_EIGENVALUE(N,Q)*HMULT_1(UM,N) )
           L_HMULT_2(UM,N,Q) = ZETA_P(UM,N) * 
     &           ( - SM*L_T_2 - L_EIGENVALUE(N,Q)*HMULT_2(UM,N) )
          ENDDO
         ENDIF
        ENDDO
       ENDDO
      ENDIF

C  Column WFS:  All layers are flagged

      IF ( DO_COLUMN_WFS ) THEN
       DO UM = 1, N_USER_STREAMS
        SM   = 1.0d0 / USER_STREAMS(UM)
        DO N = 1, NLAYERS
         DO Q = 1, N_COLUMN_WFS
          L_T_2 =   EIGENTRANS(N)   * L_T_DELT_USERM(N,UM,Q) +
     &            L_EIGENTRANS(N,Q) *   T_DELT_USERM(N,UM)
          L_T_1 = L_EIGENTRANS(N,Q) - L_T_DELT_USERM(N,UM,Q)
          L_HMULT_1(UM,N,Q) = ZETA_M(UM,N) * 
     &          ( SM*L_T_1 + L_EIGENVALUE(N,Q)*HMULT_1(UM,N) )
          L_HMULT_2(UM,N,Q) = ZETA_P(UM,N) * 
     &          ( - SM*L_T_2 - L_EIGENVALUE(N,Q)*HMULT_2(UM,N) )
         ENDDO
        ENDDO
       ENDDO
      ENDIF

C  debug
c      do n = 1, 3
c        write(56,*)L_HMULT_1(1,N,1),L_HMULT_2(1,N,1)
c      enddo

C  Finish

      RETURN
      END

C

      SUBROUTINE TWOSTREAM_L_HOM_USERSOLUTION
     I    ( NLAYERS, N_USER_STREAMS, NPARS, 
     I      N, FOURIER, DOVARY, NVARY, STREAM_VALUE, PX11,
     I      USER_STREAMS, ULP, L_XPOS, L_XNEG, U_P, U_M,
     I      OMEGA, ASYMM, L_OMEGA, L_ASYMM,
     O      L_U_XPOS, L_U_XNEG )

C  subroutine input arguments
C  --------------------------

C  Numbers

      INTEGER          NLAYERS, N_USER_STREAMS, NPARS

C  Given layer index and Fourier number (inputs)

      INTEGER          N
      INTEGER          FOURIER

C  Linearization control

      LOGICAL          DOVARY
      INTEGER          NVARY

C  Stream value and polynomial

      DOUBLE PRECISION STREAM_VALUE, PX11

C  User-defined post-processing stream directions

      DOUBLE PRECISION USER_STREAMS ( N_USER_STREAMS )
      DOUBLE PRECISION ULP          ( N_USER_STREAMS )

C  OMEGA and ASYMM, + linearizations

      DOUBLE PRECISION OMEGA ( NLAYERS )
      DOUBLE PRECISION ASYMM ( NLAYERS )
      DOUBLE PRECISION L_OMEGA ( NLAYERS,NPARS )
      DOUBLE PRECISION L_ASYMM ( NLAYERS,NPARS )

C  Linearized UP and down eigensolutions

      DOUBLE PRECISION L_XPOS(2,NLAYERS,NPARS)
      DOUBLE PRECISION L_XNEG(2,NLAYERS,NPARS)

C  Saved help variables from the original routine for User solutions

      DOUBLE PRECISION U_P(0:1)
      DOUBLE PRECISION U_M(0:1)

C  Subroutine output arguments
C  ---------------------------

C  Linearized Eigensolutions defined at user-defined stream angles

      DOUBLE PRECISION
     U        L_U_XPOS(N_USER_STREAMS,NLAYERS,NPARS),
     U        L_U_XNEG(N_USER_STREAMS,NLAYERS,NPARS)

C  Local variables
C  ---------------

      INTEGER          UM, Q
      DOUBLE PRECISION SUM_NEG, SUM_POS
      DOUBLE PRECISION OMEGA_MOM, HMU_STREAM
      DOUBLE PRECISION L_OMEGA_MOM, lu_p(0:1), lu_m(0:1)
      DOUBLE PRECISION OM_MU, OM_ULP, L_OM_MU, L_OM_ULP
 
C  Zero output and return if no solutions

      IF ( .NOT. DOVARY ) THEN
        DO Q = 1, NVARY
          DO UM = 1, N_USER_STREAMS
            L_U_XPOS(UM,N,Q) = 0.0d0
            L_U_XNEG(UM,N,Q) = 0.0d0
          ENDDO
        ENDDO 
        RETURN
      ENDIF

C  Save some useful things

      HMU_STREAM = STREAM_VALUE *  0.5d0
      OMEGA_MOM  = 3.0d0 * OMEGA(N) * ASYMM(N)

C  Start parameter loop

      DO Q = 1, NVARY

C  basic optical property variation

        L_OMEGA_MOM = L_OMEGA(N,Q)*ASYMM(N)+OMEGA(N)*L_ASYMM(N,Q)
        L_OMEGA_MOM = 3.0d0 * L_OMEGA_MOM

C  Eigenvector interpolation to user-defined angles
C  ------------------------------------------------

C  For each moment, do inner sum over computational angles
C  for the positive and negative eigenvectors

        if ( fourier.eq.0) then
          lu_p(0) = ( L_XPOS(2,N,Q) + L_XPOS(1,N,Q) ) * 0.5d0
          lu_p(1) = ( L_XPOS(2,N,Q) - L_XPOS(1,N,Q) ) * HMU_STREAM
          lu_M(0) =   lu_p(0)
          lu_M(1) = - lu_p(1)
        else
          lu_p(1) = - ( L_XPOS(2,N,Q) + L_XPOS(1,N,Q) ) * PX11 * 0.5d0
          lu_M(1) = lu_p(1)
        endif

C  Now sum over all harmonic contributions at each user-defined stream

        DO UM = 1, N_USER_STREAMS
          if (fourier.eq.0 ) then
            OM_MU    =   OMEGA_MOM * USER_STREAMS(um)
            L_OM_MU  = L_OMEGA_MOM * USER_STREAMS(um)
            sum_pos = lu_p(0) *   omega(N)   + lu_p(1) *   OM_MU
     &               + u_p(0) * l_omega(N,Q) +  u_p(1) * L_OM_MU 
            sum_neg = lu_m(0) *   omega(N)   + lu_m(1) *   OM_MU
     &               + u_m(0) * l_omega(N,Q) +  u_m(1) * L_OM_MU 
          else
            OM_ULP   =   OMEGA_MOM * ulp(um)
            L_OM_ULP = L_OMEGA_MOM * ulp(um)
            sum_pos = lu_p(1) * OM_ULP + u_p(1) * L_OM_ULP
            sum_neg = lu_m(1) * OM_ULP + u_m(1) * L_OM_ULP
          endif
          L_U_XPOS(UM,N,Q) = SUM_POS
          L_U_XNEG(UM,N,Q) = SUM_NEG
        ENDDO

C  end parameter loop

      enddo

C  debug
c        if (fourier.eq.1)then
c       write(56,'(2i4,1p3e24.12)')n,1,L_U_XPOS(1,N,1),L_U_XNEG(1,N,1)
c       write(56,'(2i4,1p3e24.12)')n,2,L_U_XPOS(1,N,2),L_U_XNEG(1,N,2)
c        endif

C  Finish

      RETURN
      END

C 


      SUBROUTINE TWOSTREAM_L_BEAM_USERSOLUTION
     I    ( DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, 
     I      NLAYERS, NBEAMS, N_USER_STREAMS, NPARS, N, FOURIER, IBEAM,
     I      DO_COLUMN_LINEARIZATION, DO_PROFILE_LINEARIZATION,
     I      LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_COLUMN_WFS,
     I      FLUX_FACTOR, LAYER_PIS_CUTOFF,
     I      USER_STREAMS, STREAM_VALUE, PX11, ULP, X0, POX,
     I      OMEGA, ASYMM, L_OMEGA, L_ASYMM, W_HELP, L_WVEC,
     O      L_U_WPOS1, L_U_WPOS2, L_U_WNEG1, L_U_WNEG2 )

C  subroutine arguments
C  --------------------

C  Flags

      LOGICAL          DO_UPWELLING, DO_DNWELLING

C  Plane Parallel flag

      LOGICAL          DO_PLANE_PARALLEL

C  Numbers

      INTEGER          NLAYERS, NBEAMS, N_USER_STREAMS, NPARS

C  Given layer index and Fourier number, Beam number (inputs)

      INTEGER          N
      INTEGER          FOURIER
      INTEGER          IBEAM

C  LInearization flags

      LOGICAL          DO_COLUMN_LINEARIZATION
      LOGICAL          DO_PROFILE_LINEARIZATION

C  Linearization control

      LOGICAL          LAYER_VARY_FLAG(NLAYERS)
      INTEGER          LAYER_VARY_NUMBER(NLAYERS)
      INTEGER          N_COLUMN_WFS

C  Flux factor

      DOUBLE PRECISION FLUX_FACTOR

C  Last layer to include Particular integral solution

      INTEGER          LAYER_PIS_CUTOFF(NBEAMS)

C  Stream value and polynomial

      DOUBLE PRECISION STREAM_VALUE, PX11

C  Beam SZA cosines

      DOUBLE PRECISION X0(NBEAMS)
      DOUBLE PRECISION POX(NBEAMS)

C  OMEGA and ASYMM, + linearizations

      DOUBLE PRECISION OMEGA ( NLAYERS )
      DOUBLE PRECISION ASYMM ( NLAYERS )
      DOUBLE PRECISION L_OMEGA ( NLAYERS,NPARS )
      DOUBLE PRECISION L_ASYMM ( NLAYERS,NPARS )

C  User streams

      DOUBLE PRECISION USER_STREAMS ( N_USER_STREAMS )
      DOUBLE PRECISION ULP          ( N_USER_STREAMS )

C  Saved help variables

      DOUBLE PRECISION W_HELP(0:1)

C  Linearized Beam solutions

      DOUBLE PRECISION L_WVEC(2,NLAYERS,0:NLAYERS,NPARS)

C  Subroutine output arguments
C  ---------------------------

C  Linearized Particular beam solutions at user-defined stream angles
C       1 = Direct  term contribution
C       2 = Diffuse term contribution

      DOUBLE PRECISION
     U        L_U_WPOS1(N_USER_STREAMS,NLAYERS,NPARS),
     U        L_U_WPOS2(N_USER_STREAMS,NLAYERS,0:NLAYERS,NPARS)

      DOUBLE PRECISION
     U        L_U_WNEG1(N_USER_STREAMS,NLAYERS,NPARS),
     U        L_U_WNEG2(N_USER_STREAMS,NLAYERS,0:NLAYERS,NPARS)

C  Local variables
C  ---------------

      LOGICAL          DO_FIRST
      INTEGER          UM, NV, K, Q
      DOUBLE PRECISION POS1, POS2, F1, OMEGA_MOM, PI4, HMUS, HPX11
      DOUBLE PRECISION L_OMEGA_MOM, l_h1(0:1), l_h2(0:1), l_w(0:1)

C  Profiles and Columns

      IF ( DO_PROFILE_LINEARIZATION ) THEN
        NV = N
      ELSE IF ( DO_COLUMN_LINEARIZATION ) THEN
        NV = 0
      ENDIF

C  No particular solution beyond the cutoff layer
C  ... Zero the user solutions and exit

C  Only a solution if the layer is not below Cutoff.

      DO_FIRST = ( N .LE. LAYER_PIS_CUTOFF(IBEAM) )

C  If no solution or no variation, zero output and go to Part 2

      IF ( .NOT. LAYER_VARY_FLAG(N) .OR. .NOT. DO_FIRST ) THEN
        DO UM = 1, N_USER_STREAMS
          IF ( DO_UPWELLING ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_U_WPOS1(UM,N,Q)    = 0.0d0
              L_U_WPOS2(UM,N,NV,Q) = 0.0d0
            ENDDO
          ENDIF
          IF ( DO_DNWELLING ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_U_WNEG1(UM,N,Q)    = 0.0d0
              L_U_WNEG2(UM,N,NV,Q) = 0.0d0
            ENDDO
          ENDIF
        ENDDO
        GO TO 2222
      ENDIF

C  set up numbers

      PI4 = 4.0d0 * dacos(-1.0d0)
      F1 = FLUX_FACTOR / PI4
      OMEGA_MOM  = 3.0d0 * OMEGA(N) * ASYMM(N)
      HMUS       = STREAM_VALUE * 0.5d0
      IF ( FOURIER.EQ.1 ) HPX11 = PX11 * 0.5d0
      
C  Start parameter loop

      DO Q = 1, LAYER_VARY_NUMBER(N)

C  basic optical property variation

        L_OMEGA_MOM = L_OMEGA(N,Q)*ASYMM(N)+OMEGA(N)*L_ASYMM(N,Q)
        L_OMEGA_MOM = 3.0d0 * L_OMEGA_MOM

C  Eigenvector interpolation to user-defined angles
C  ------------------------------------------------

C  For each moment, do inner sum over computational angles
C  for the positive and negative eigenvectors

        if ( fourier.eq.0) then
          l_w(0)   = ( L_WVEC(2,N,NV,Q) + L_WVEC(1,N,NV,Q) ) * 0.5d0
          l_w(1)   = ( L_WVEC(2,N,NV,Q) - L_WVEC(1,N,NV,Q) ) * HMUS
          l_h1(0)  =   l_omega(n,q) * F1
          l_h1(1)  = - l_omega_mom  * F1 * x0(ibeam)
          l_h2(0)  = l_w(0) * omega(n)  + w_help(0) * l_omega(n,q)
          l_h2(1)  = l_w(1) * omega_mom + w_help(1) * l_omega_mom
        else
          l_w(1)   = - ( L_WVEC(2,N,NV,Q) + L_WVEC(1,N,NV,Q) )* HPX11
          l_h1(1)  = - l_omega_mom * F1 * POX(ibeam)
          l_h2(1)  = l_w(1) * omega_mom + w_help(1) * l_omega_mom
        endif

C  Now sum over all harmonic contributions at each user-defined stream
C  Distinguish between upwelling and downwelling

        IF ( DO_UPWELLING ) THEN
          DO UM = 1, N_USER_STREAMS
            if (fourier.eq.0 ) then
              pos1 = l_h1(0) + l_h1(1) * user_streams(UM)
              pos2 = l_h2(0) + l_h2(1) * user_streams(UM)
            else
              pos1 = l_h1(1) * ULP(UM)
              pos2 = l_h2(1) * ULP(UM)
            endif
            l_U_WPOS1(UM,N,Q)    = POS1
            l_U_WPOS2(UM,N,NV,Q) = POS2
          ENDDO
        ENDIF

        IF ( DO_DNWELLING ) THEN
          DO UM = 1, N_USER_STREAMS
            if (fourier.eq.0 ) then
              pos1 = l_h1(0) - l_h1(1) * user_streams(UM)
              pos2 = l_h2(0) - l_h2(1) * user_streams(UM)
            else
              pos1 = l_h1(1) * ULP(UM)
              pos2 = l_h2(1) * ULP(UM)
            endif
            l_U_WNEG1(UM,N,Q)    = POS1
            l_U_WNEG2(UM,N,NV,Q) = POS2
          ENDDO
        ENDIF

C  End parameter loop

      enddo

C  debug
c        if (fourier.eq.0)then
c       write(*,'(i4,1p3e24.12)')n,L_U_WPOS1(1,N,1),L_U_WPOS2(1,N,N,1)
c       write(*,'(i4,1p3e24.12)')n,L_U_WPOS1(1,N,2),L_U_WPOS2(1,N,N,2)
c        if ( n.eq.3)pause
c        endif

C  Part 2.
C  =======

C  Continuation point

 2222 continue

C  Only these extra variations when the following conditions
C  are not satisfied (pseudo-spherical, layer > 1)
C  Profiles only

      IF ( N .EQ. 1 )                RETURN
      IF ( DO_PLANE_PARALLEL )       RETURN
      IF ( DO_COLUMN_LINEARIZATION ) RETURN

C  If no solution, zero output and exit

      IF ( .NOT. DO_FIRST ) THEN
        DO K = 1, N - 1
          DO Q = 1, LAYER_VARY_NUMBER(K)
            DO UM = 1, N_USER_STREAMS
              L_U_WPOS2(UM,N,K,Q) = 0.0d0
              L_U_WNEG2(UM,N,K,Q) = 0.0d0
            ENDDO
          ENDDO
        ENDDO
        RETURN
      ENDIF

C  start loop over all layers above N

      DO K = 1, N - 1

C  only do if layer K has some variation

        IF ( LAYER_VARY_FLAG(K) ) THEN

C  start parameter loop

          DO Q = 1, LAYER_VARY_NUMBER(K)

C  For each moment, do inner sum over computational angles
C  for the positive and negative eigenvectors

            if ( fourier.eq.0) then
              l_w(0)   = ( L_WVEC(2,N,K,Q) + L_WVEC(1,N,K,Q) ) * 0.5d0
              l_w(1)   = ( L_WVEC(2,N,K,Q) - L_WVEC(1,N,K,Q) ) * HMUS
              l_h2(0)  = l_w(0) * omega(n) 
              l_h2(1)  = l_w(1) * omega_mom
            else
              l_w(1)   = - ( L_WVEC(2,N,K,Q) + L_WVEC(1,N,K,Q) ) * HPX11
              l_h2(1)  = l_w(1) * omega_mom 
            endif

C  Now sum over all harmonic contributions at each user-defined stream
C  Distinguish between upwelling and downwelling

            IF ( DO_UPWELLING ) THEN
              DO UM = 1, N_USER_STREAMS
                if (fourier.eq.0 ) then
                  pos2 = l_h2(0) + l_h2(1) * user_streams(UM)
                else
                  pos2 = l_h2(1) * ULP(UM)
                endif
                l_U_WPOS2(UM,N,K,Q) = POS2
              ENDDO
            ENDIF

            IF ( DO_DNWELLING ) THEN
              DO UM = 1, N_USER_STREAMS
                if (fourier.eq.0 ) then
                  pos2 = l_h2(0) - l_h2(1) * user_streams(UM)
                else
                  pos2 = l_h2(1) * ULP(UM)
                endif
                l_U_WNEG2(UM,N,K,Q) = POS2
              ENDDO
            ENDIF

C  end parameter loop

          ENDDO

C  end K-layer loop

        ENDIF
      ENDDO

C  debug
c        if (fourier.eq.1)then
c       write(*,'(i4,1p3e24.12)')n,(L_U_WPOS2(1,N,K,1),K=1,N)
c       write(*,'(i4,1p3e24.12)')n,(L_U_WPOS2(1,N,K,2),K=1,N)
c        if ( n.eq.3) pause
c        endif

C  Finish

      RETURN
      END

