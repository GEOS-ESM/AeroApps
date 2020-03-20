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
C #              TWOSTREAM_HOM_SOLUTION                         #
C #              TWOSTREAM_HOM_USERSOLUTION                     #
C #              TWOSTREAM_HMULT_MASTER                         #
C #                                                             #
C #     PArticular integrals                                    #
C #                                                             #
C #              TWOSTREAM_BEAM_SOLUTION                        #
C #              TWOSTREAM_BEAM_USERSOLUTION                    #
C #                                                             #
C ###############################################################

      SUBROUTINE TWOSTREAM_HOM_SOLUTION
     I    ( NLAYERS, N, FOURIER, STREAM_VALUE, PXSQ,
     I      OMEGA, ASYMM, DELTAU_VERT, 
     O      SAB, DAB, EIGENVALUE, EIGENTRANS,
     O      XPOS, XNEG )

C  subroutine arguments
C  --------------------

C  Given layer index and Fourier number (inputs)

      INTEGER          NLAYERS, N
      INTEGER          FOURIER

C  Stream value

      DOUBLE PRECISION STREAM_VALUE

C  Polynomials

      DOUBLE PRECISION PXSQ
      
C  OMEGA and ASYMM

      DOUBLE PRECISION OMEGA ( NLAYERS )
      DOUBLE PRECISION ASYMM ( NLAYERS )

C  optical thickness

      DOUBLE PRECISION DELTAU_VERT(NLAYERS)

C  Solutions to the homogeneous RT equations 
C  -----------------------------------------

C  local matrices for eigenvalue computation

      DOUBLE PRECISION SAB(NLAYERS), DAB(NLAYERS)

C  Eigensolutions

      DOUBLE PRECISION EIGENVALUE(NLAYERS)
      DOUBLE PRECISION EIGENTRANS(NLAYERS)

C  UP and down solutions

      DOUBLE PRECISION XPOS(2,NLAYERS)
      DOUBLE PRECISION XNEG(2,NLAYERS)

C  Local variables
C  ---------------

      DOUBLE PRECISION DP, DM, XINV, OMEGA_ASYMM_3, DIFVEC
    
C  Develop Sum and Difference matrices

      XINV = 1.0d0 / STREAM_VALUE
      OMEGA_ASYMM_3 = 3.0d0 * OMEGA(N) * ASYMM(N)
      if ( fourier.eq.0) then
        DP = OMEGA(N) + PXSQ * OMEGA_ASYMM_3
        DM = OMEGA(N) - PXSQ * OMEGA_ASYMM_3
      Else if ( fourier .eq. 1 ) then
        DP = OMEGA_ASYMM_3 * PXSQ
        DM = OMEGA_ASYMM_3 * PXSQ
      ENDIF
      SAB(N) = XINV * ( ( DP + DM ) * 0.5d0 - 1.0d0 )
      DAB(N) = XINV * ( ( DP - DM ) * 0.5d0 - 1.0d0 )
      EIGENVALUE(N) = DSQRT(SAB(N)*DAB(N))
      EIGENTRANS(N) = DEXP ( - EIGENVALUE(N) * DELTAU_VERT(N) )

C  Auxiliary equation to get up and down solutions

      DIFVEC = - SAB(N) / EIGENVALUE(N)
      XPOS(1,N) = 0.5d0 * ( 1.0d0 + DIFVEC )
      XPOS(2,N) = 0.5d0 * ( 1.0d0 - DIFVEC )

C  Symmetry

      XNEG(1,N) = XPOS(2,N)
      XNEG(2,N) = XPOS(1,N)

C  debug
c      if (fourier.eq.1)write(57,'(i4,1pe24.12)')n,EIGENTRANS(N)

C  Finish

      RETURN
      END

C

      SUBROUTINE TWOSTREAM_BEAM_SOLUTION
     I    ( NLAYERS, NBEAMS, N, FOURIER, IBEAM,
     I      FLUX_FACTOR, LAYER_PIS_CUTOFF, STREAM_VALUE, X0, PX0X,
     I      AVERAGE_SECANT, INITIAL_TRANS, T_DELT_MUBAR,
     I      OMEGA, ASYMM, DELTAU_VERT, SAB, DAB, EIGENVALUE,
     O      QSUMVEC, QDIFVEC, QVEC,
     O      WVEC, WUPPER, WLOWER )

C  subroutine arguments
C  --------------------

C  Numbers

      INTEGER          NLAYERS, NBEAMS

C  Given layer index and Fourier number, Beam number (inputs)

      INTEGER          N
      INTEGER          FOURIER
      INTEGER          IBEAM

C  Flux factor

      DOUBLE PRECISION FLUX_FACTOR

C  Last layer to include Particular integral solution

      INTEGER          LAYER_PIS_CUTOFF(NBEAMS)

C  Stream value

      DOUBLE PRECISION STREAM_VALUE

C  Beam SZA cosines

      DOUBLE PRECISION X0  (NBEAMS)
      DOUBLE PRECISION PX0X(NBEAMS)

C  Average-secant and initial tramsittance factors for solar beams.

      DOUBLE PRECISION
     &     INITIAL_TRANS  ( NLAYERS, NBEAMS ),
     &     AVERAGE_SECANT ( NLAYERS, NBEAMS ),
     &     T_DELT_MUBAR   ( NLAYERS, NBEAMS )

C  OMEGA and ASYMM

      DOUBLE PRECISION OMEGA ( NLAYERS )
      DOUBLE PRECISION ASYMM ( NLAYERS )

C  optical thickness

      DOUBLE PRECISION DELTAU_VERT(NLAYERS)

C  local matrices from eigenvalue computation

      DOUBLE PRECISION SAB(NLAYERS), DAB(NLAYERS)

C  Eigenvalues

      DOUBLE PRECISION EIGENVALUE(NLAYERS)

C  Output variables
C  ----------------

C  Beam solution

      DOUBLE PRECISION WVEC(2,NLAYERS)

C  Solutions at layer boundaries

      DOUBLE PRECISION WUPPER(2,NLAYERS)
      DOUBLE PRECISION WLOWER(2,NLAYERS)

C  Auxiliary vectors

      DOUBLE PRECISION QDIFVEC(NLAYERS)
      DOUBLE PRECISION QSUMVEC(NLAYERS)
      DOUBLE PRECISION QVEC   (NLAYERS)

C  help variables
C  --------------

      INTEGER          I
      DOUBLE PRECISION TP, TM, INV_X0SQ, SECBAR, XINV, F1
      DOUBLE PRECISION HELP, TRANS1, TRANS2
      DOUBLE PRECISION QMAT, QDIF, OMEGA_ASYMM, PI4

C  Flux factor

      PI4 = 4.0d0 * dacos(-1.0d0)
      F1 = FLUX_FACTOR / PI4

C  No particular solution beyond the cutoff layer
C  Or no scattering in this layer...
C  ... Zero the boundary layer values and exit

      IF ( N .GT. LAYER_PIS_CUTOFF(IBEAM) ) THEN
        DO I = 1, 2
          WUPPER(I,N) = 0.0d0
          WLOWER(I,N) = 0.0d0
        ENDDO
        RETURN
      ENDIF

C  set local values

      SECBAR   = AVERAGE_SECANT(N,IBEAM)
      INV_X0SQ = SECBAR * SECBAR

C  Set up sum and differences for Beam source terms
C  ( sum may be required again in linearization )

      XINV = 1.0d0 / STREAM_VALUE
      OMEGA_ASYMM = OMEGA(N) * ASYMM(N) * 3.0d0
      if ( fourier.eq.0) then
        TP = OMEGA(N) + PX0X(IBEAM) * OMEGA_ASYMM
        TM = OMEGA(N) - PX0X(IBEAM) * OMEGA_ASYMM
      Else if ( fourier .eq. 1 ) then
        TP = PX0X(IBEAM) * OMEGA_ASYMM
        TM = PX0X(IBEAM) * OMEGA_ASYMM
      ENDIF
      QSUMVEC(N) =  F1 * ( TP + TM ) * XINV
      QDIFVEC(N) =  F1 * ( TP - TM ) * XINV

C  the reduced problem: QMAT. W = QVEC (Overwrite QVEC)

      QMAT = EIGENVALUE(N) * EIGENVALUE(N) - INV_X0SQ
      HELP = - DAB(N) * QSUMVEC(N)
      QVEC(N) = HELP + QDIFVEC(N) * SECBAR
      QVEC(N) = QVEC(N) / QMAT

C  Restore up and down solutions

      HELP = - SAB(N) * QVEC(N)
      QDIF = ( HELP - QSUMVEC(N) ) / SECBAR
      WVEC(1,N) = 0.5d0 * ( QVEC(N) + QDIF )
      WVEC(2,N) = 0.5d0 * ( QVEC(N) - QDIF )

C  Values at the layer boundaries
C  (transmittance factors have been determined in SETUPS module)

      TRANS1 = INITIAL_TRANS(N,IBEAM)
      TRANS2 = T_DELT_MUBAR(N,IBEAM) * TRANS1
      DO I = 1, 2
        WUPPER(I,N) = WVEC(I,N)*TRANS1
        WLOWER(I,N) = WVEC(I,N)*TRANS2
      ENDDO

C  debug
c      if (fourier.eq.0)write(*,'(i4,1pe24.12)')n,WVEC(1,N)

C  Finish

      RETURN
      END

C

      SUBROUTINE TWOSTREAM_HMULT_MASTER
     I     ( NLAYERS, N_USER_STREAMS, USER_STREAMS,
     I       EIGENVALUE, EIGENTRANS, T_DELT_USERM,  
     O       ZETA_M, ZETA_P, HMULT_1, HMULT_2 )

C  Input arguments
C  ===============

C  Numbers

      INTEGER          NLAYERS, N_USER_STREAMS

C  User streams

      DOUBLE PRECISION USER_STREAMS ( N_USER_STREAMS )

C  Transmittance factors for user-defined stream angles

      DOUBLE PRECISION
     &     T_DELT_USERM ( NLAYERS, N_USER_STREAMS )

C  Eigensolutions

      DOUBLE PRECISION EIGENVALUE(NLAYERS)
      DOUBLE PRECISION EIGENTRANS(NLAYERS)

C  Output = Global multipliers
C  ===========================

C  coefficient functions for user-defined angles

      DOUBLE PRECISION 
     &      ZETA_M(N_USER_STREAMS,NLAYERS),
     &      ZETA_P(N_USER_STREAMS,NLAYERS)

C  Integrated homogeneous solution multipliers, whole layer

      DOUBLE PRECISION 
     &      HMULT_1(N_USER_STREAMS,NLAYERS),
     &      HMULT_2(N_USER_STREAMS,NLAYERS)

C  Local variables
C  ---------------

      INTEGER          UM, N
      DOUBLE PRECISION UDEL, SM, RHO_M, RHO_P
      DOUBLE PRECISION ZDEL, ZUDEL, THETA_1, THETA_2

C  whole layer multipliers
C  -----------------------

C  Start loops over layers and user-streams
C    Only done if layers are flagged

      DO N = 1, NLAYERS
        DO UM = 1, N_USER_STREAMS
          UDEL = T_DELT_USERM(N,UM)
          SM   = 1.0d0 / USER_STREAMS(UM)
          RHO_P = SM + EIGENVALUE(N)
          RHO_M = SM - EIGENVALUE(N)
          ZETA_P(UM,N) = 1.0d0 / RHO_P
          ZETA_M(UM,N) = 1.0d0 / RHO_M
          ZDEL    = EIGENTRANS(N)
          ZUDEL   = ZDEL * UDEL
          THETA_2 = 1.0d0 - ZUDEL
          THETA_1 = ZDEL  - UDEL
          HMULT_1(UM,N) = SM * THETA_1 * ZETA_M(UM,N)
          HMULT_2(UM,N) = SM * THETA_2 * ZETA_P(UM,N)
        ENDDO
      ENDDO

C  debug
c      do n = 1, 3
c        write(57,*)HMULT_1(1,N),HMULT_2(1,N)
c      enddo

C  Finish

      RETURN
      END


C

      SUBROUTINE TWOSTREAM_HOM_USERSOLUTION
     I    ( NLAYERS, N_USER_STREAMS, N, FOURIER, STREAM_VALUE, PX11,
     I      USER_STREAMS, ULP, XPOS, XNEG, OMEGA, ASYMM,
     O      U_XPOS, U_XNEG, U_HELP_P, U_HELP_M )

C  subroutine input arguments
C  --------------------------

C  Numbers

      INTEGER          NLAYERS, N_USER_STREAMS

C  Given layer index and Fourier number (inputs)

      INTEGER          N
      INTEGER          FOURIER

C  Stream value and polynomial

      DOUBLE PRECISION STREAM_VALUE, PX11

C  User-defined post-processing stream directions

      DOUBLE PRECISION USER_STREAMS ( N_USER_STREAMS )
      DOUBLE PRECISION ULP          ( N_USER_STREAMS )

C  OMEGA and ASYMM

      DOUBLE PRECISION OMEGA ( NLAYERS )
      DOUBLE PRECISION ASYMM ( NLAYERS )

C  UP and down solutions

      DOUBLE PRECISION XPOS(2,NLAYERS)
      DOUBLE PRECISION XNEG(2,NLAYERS)

C  Subroutine output arguments
C  ---------------------------

C  Saved help variables

      DOUBLE PRECISION U_HELP_P(0:1)
      DOUBLE PRECISION U_HELP_M(0:1)

C  Eigenvectors defined at user-defined stream angles
C     EP for the positive KEIGEN values, EM for -ve KEIGEN

      DOUBLE PRECISION
     U        U_XPOS(N_USER_STREAMS,NLAYERS),
     U        U_XNEG(N_USER_STREAMS,NLAYERS)

C  Local variables
C  ---------------

      INTEGER          UM
      DOUBLE PRECISION SUM_NEG, SUM_POS
      DOUBLE PRECISION OMEGA_MOM, HMU_STREAM

C  zero the user solutions

      DO UM = 1, N_USER_STREAMS
        U_XPOS(UM,N) = 0.0d0
        U_XNEG(UM,N) = 0.0d0
      ENDDO

C  Eigenvector interpolation to user-defined angles
C  ------------------------------------------------

C  For each moment, do inner sum over computational angles
C  for the positive and negative eigenvectors

      HMU_STREAM = 0.5d0 * STREAM_VALUE     
      if ( fourier.eq.0) then
        u_help_p(0) = ( XPOS(2,N) + XPOS(1,N) ) * 0.5d0
        u_help_p(1) = ( XPOS(2,N) - XPOS(1,N) ) * HMU_STREAM
        u_help_M(0) =   u_help_p(0)
        u_help_M(1) = - u_help_p(1)
      else
        u_help_p(1) = - ( XPOS(2,N) + XPOS(1,N) ) * PX11 * 0.5d0
        u_help_M(1) = u_help_p(1)
      endif

C  Now sum over all harmonic contributions at each user-defined stream

      OMEGA_MOM = 3.0d0 * OMEGA(N) * ASYMM(N)
      DO UM = 1, N_USER_STREAMS
        if (fourier.eq.0 ) then
          sum_pos = u_help_p(0) * omega(N)
     &           +  u_help_p(1) * omega_mom * user_streams(um)
          sum_neg = u_help_m(0) * omega(N)
     &           +  u_help_m(1) * omega_mom * user_streams(um)
        else
          sum_pos = u_help_p(1) * omega_mom * ulp(um)
          sum_neg = u_help_m(1) * omega_mom * ulp(um)
        endif
        U_XPOS(UM,N) = SUM_POS
        U_XNEG(UM,N) = SUM_NEG
      ENDDO

C  debug
c      if (fourier.eq.1)
c     &   write(57,'(i4,1p2e24.12)')n,u_xpos(1,N),u_xneg(1,N)

C  Finish

      RETURN
      END

C

      SUBROUTINE TWOSTREAM_BEAM_USERSOLUTION
     I    ( DO_UPWELLING, DO_DNWELLING,
     I      NLAYERS, NBEAMS, N_USER_STREAMS, N, FOURIER, IBEAM,
     I      FLUX_FACTOR, LAYER_PIS_CUTOFF, STREAM_VALUE, PX11, X0, POX,
     I      OMEGA, ASYMM, USER_STREAMS, ULP, WVEC,
     O      W_HELP, U_WPOS1, U_WPOS2, U_WNEG1, U_WNEG2 )

C  subroutine arguments
C  --------------------

C  Flags

      LOGICAL          DO_UPWELLING, DO_DNWELLING

C  Numbers

      INTEGER          NLAYERS, NBEAMS, N_USER_STREAMS

C  Given layer index and Fourier number, Beam number (inputs)

      INTEGER          N
      INTEGER          FOURIER
      INTEGER          IBEAM

C  Flux factor

      DOUBLE PRECISION FLUX_FACTOR

C  Last layer to include Particular integral solution

      INTEGER          LAYER_PIS_CUTOFF(NBEAMS)

C  Stream value and polynomial

      DOUBLE PRECISION STREAM_VALUE, PX11

C  Beam SZA cosines

      DOUBLE PRECISION X0(NBEAMS)
      DOUBLE PRECISION POX(NBEAMS)

C  OMEGA and ASYMM

      DOUBLE PRECISION OMEGA ( NLAYERS )
      DOUBLE PRECISION ASYMM ( NLAYERS )

C  User streams

      DOUBLE PRECISION USER_STREAMS ( N_USER_STREAMS )
      DOUBLE PRECISION ULP          ( N_USER_STREAMS )

C  Beam solution

      DOUBLE PRECISION WVEC(2,NLAYERS)

C  Subroutine output arguments
C  ---------------------------

C  Saved help variables

      DOUBLE PRECISION W_HELP(0:1)

C  Particular beam solutions at user-defined stream angles

      DOUBLE PRECISION
     U        U_WPOS1(N_USER_STREAMS,NLAYERS),
     U        U_WNEG1(N_USER_STREAMS,NLAYERS)

      DOUBLE PRECISION
     U        U_WPOS2(N_USER_STREAMS,NLAYERS),
     U        U_WNEG2(N_USER_STREAMS,NLAYERS)

C  Local variables
C  ---------------

      INTEGER          UM
      DOUBLE PRECISION POS1, POS2, F1, OMEGA_MOM, PI4
      DOUBLE PRECISION HELP1(0:1), HELP2(0:1), HMU_STREAM

C  No particular solution beyond the cutoff layer
C  ... Zero the user solutions and exit

      IF ( N .GT. LAYER_PIS_CUTOFF(IBEAM) ) THEN
        DO UM = 1, N_USER_STREAMS
          IF ( DO_UPWELLING ) THEN
            U_WPOS1(UM,N) = 0.0d0
            U_WPOS2(UM,N) = 0.0d0
          ENDIF
          IF ( DO_DNWELLING ) THEN
            U_WNEG1(UM,N) = 0.0d0
            U_WNEG2(UM,N) = 0.0d0
          ENDIF
        ENDDO
        RETURN
      ENDIF

C  Scattering solutions
C  ====================

C  Starter quantities

      PI4 = 4.0d0 * dacos(-1.0d0)
      F1 = FLUX_FACTOR / PI4
      OMEGA_MOM = 3.0d0 * OMEGA(N) * ASYMM(N)
      HMU_STREAM = STREAM_VALUE * 0.5d0

C  For each moment do inner sum over computational angles

      if ( fourier.eq.0) then
        w_help(0) = ( WVEC(2,N) + WVEC(1,N) ) * 0.5d0
        w_help(1) = ( WVEC(2,N) - WVEC(1,N) ) * HMU_STREAM
        help1(0)  =   omega(n)  * F1
        help1(1)  = - omega_mom * F1 * x0(ibeam)
        help2(0)  = w_help(0) * omega(n)
        help2(1)  = w_help(1) * omega_mom
      else
        w_help(1) = - ( WVEC(2,N) + WVEC(1,N) ) * PX11 * 0.5d0
        help1(1)  = - omega_mom * F1 * POX(IBEAM)
        help2(1)  = w_help(1) * omega_mom
      endif

C  Now sum over all harmonic contributions at each user-defined stream
C  Distinguish between upwelling and downwelling

      IF ( DO_UPWELLING ) THEN
        DO UM = 1, N_USER_STREAMS
          if (fourier.eq.0 ) then
            pos1 = help1(0) + help1(1)* user_streams(um)
            pos2 = help2(0) + help2(1)* user_streams(um)
          else
            pos1 = help1(1)* ULP(UM)
            pos2 = help2(1)* ULP(UM)
          endif
          U_WPOS1(UM,N) = POS1
          U_WPOS2(UM,N) = POS2
        ENDDO
      ENDIF

      IF ( DO_DNWELLING ) THEN
        DO UM = 1, N_USER_STREAMS
          if (fourier.eq.0 ) then
            pos1 = help1(0) - help1(1)* user_streams(um)
            pos2 = help2(0) - help2(1)* user_streams(um)
          else
            pos1 = help1(1) * ulp(UM)
            pos2 = help2(1) * ulp(UM)
          endif
          U_WNEG1(UM,N) = POS1
          U_WNEG2(UM,N) = POS2
        ENDDO
      ENDIF

C  debug
c      if (fourier.eq.0)
c     &   write(57,'(i4,1p2e24.12)')n,u_wpos2(1,N)

C  Finish

      RETURN
      END

