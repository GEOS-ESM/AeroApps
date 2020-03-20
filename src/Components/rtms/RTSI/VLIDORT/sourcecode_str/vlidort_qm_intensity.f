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
C #       VLIDORT_INTEGRATED_OUTPUT (master)                    #
C #                                                             #
C #              QUADINTENS_LEVEL_UP                            #
C #              QUADINTENS_LEVEL_DN                            #
C #              QUADINTENS_OFFGRID_UP                          #
C #              QUADINTENS_OFFGRID_DN                          #
C #                                                             #
C ###############################################################

      SUBROUTINE VLIDORT_INTEGRATED_OUTPUT
     I  ( DO_INCLUDE_MVOUTPUT,
     I    DO_INCLUDE_DIRECTBEAM,
     I    DO_INCLUDE_THERMEMISS,
     I    FLUX_MULTIPLIER, IBEAM )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of main model variables (input to this module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  include file of result variables (module output stored here)

      INCLUDE '../includes/VLIDORT_RESULTS.VARS'

C  Subroutine input arguments
C  --------------------------

C  Control inclusion

      LOGICAL          DO_INCLUDE_MVOUTPUT
      LOGICAL          DO_INCLUDE_DIRECTBEAM
      LOGICAL          DO_INCLUDE_THERMEMISS

C  PI number and flux multiplier

      INTEGER          IBEAM
      DOUBLE PRECISION FLUX_MULTIPLIER

C  local variables
C  ---------------

C  Local qudrature output (for debug)

      DOUBLE PRECISION QSTOKES_F 
     &       ( MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )  

C  help variables

      INTEGER          I, IDIR, WDIR, UTA, UT, N, NLEVEL, O1
      DOUBLE PRECISION SUM_MI, SUM_FX, FTRANS
      DOUBLE PRECISION DIRECT_TRANS, DIRECT_FLUX, DIRECT_MEANI

C  mean intensity and flux
C  -----------------------

C  direction loop

      DO IDIR = 1, N_DIRECTIONS
        WDIR = WHICH_DIRECTIONS(IDIR)

C  Upwelling Stokes output at Quadrature angles

        IF ( WDIR .EQ. UPIDX ) THEN
          DO UTA = 1, N_USER_LEVELS
            NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
              UT = PARTLAYERS_OUTINDEX(UTA)
              N  = PARTLAYERS_LAYERIDX(UT)
              CALL QUADINTENS_OFFGRID_UP ( N, UTA, UT, IBEAM,
     &            FLUX_MULTIPLIER, DO_INCLUDE_THERMEMISS, QSTOKES_F )
            ELSE
              CALL QUADINTENS_LEVEL_UP ( NLEVEL, UTA,
     &            FLUX_MULTIPLIER, QSTOKES_F )
            ENDIF
          ENDDO
        ENDIF

C  Downwelling Stokes output at Quadrature angles

        IF ( WDIR .EQ. DNIDX ) THEN
          DO UTA = 1, N_USER_LEVELS
            NLEVEL = UTAU_LEVEL_MASK_DN(UTA)
            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
              UT = PARTLAYERS_OUTINDEX(UTA)
              N  = PARTLAYERS_LAYERIDX(UT)
              CALL QUADINTENS_OFFGRID_DN ( N, UTA, UT, IBEAM,
     &            FLUX_MULTIPLIER, DO_INCLUDE_THERMEMISS, QSTOKES_F )
            ELSE
              CALL QUADINTENS_LEVEL_DN ( NLEVEL, UTA,
     &            FLUX_MULTIPLIER, QSTOKES_F )
            ENDIF
          ENDDO
        ENDIF

C  Mean Intensity and Flux  output
C  -------------------------------

        IF ( DO_INCLUDE_MVOUTPUT ) THEN

C  Diffuse term output

          DO UTA = 1, N_USER_LEVELS
            DO O1 = 1, NSTOKES
              SUM_MI = ZERO
              SUM_FX = ZERO
              DO I = 1, NSTREAMS
                SUM_MI = SUM_MI + QUAD_WEIGHTS(I) * QSTOKES_F(UTA,I,O1)
                SUM_FX = SUM_FX + QUAD_STRMWTS(I) * QSTOKES_F(UTA,I,O1)
              ENDDO
              MEAN_STOKES(UTA,IBEAM,O1,WDIR) = SUM_MI * HALF
              FLUX_STOKES(UTA,IBEAM,O1,WDIR) = SUM_FX * PI2
           ENDDO
          ENDDO

C  nothing to do if no solar sources

          IF ( .NOT. DO_INCLUDE_DIRECTBEAM ) GO TO 455

C  For the downward direction, add the direct beam contributions

          IF ( WDIR .EQ. DNIDX ) THEN

C  loop over all the output optical depths

           DO UTA = 1, N_USER_LEVELS

C  For the offgrid values.......
C     .....Only contributions for layers above the PI cutoff

            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
              UT = PARTLAYERS_OUTINDEX(UTA)
              N  = PARTLAYERS_LAYERIDX(UT)
              IF ( N .LE. LAYER_PIS_CUTOFF(IBEAM) ) THEN
                DIRECT_TRANS =
     &             INITIAL_TRANS(N,IBEAM) * T_UTDN_MUBAR(UT,IBEAM)
                DO O1 = 1, NSTOKES
                  FTRANS = FLUXVEC(O1) * DIRECT_TRANS
                  DIRECT_FLUX  = FTRANS * LOCAL_CSZA(N,IBEAM)
                  DIRECT_MEANI = FTRANS / PI4
                  MEAN_DIRECT(UTA,IBEAM,O1) = DIRECT_MEANI
                  FLUX_DIRECT(UTA,IBEAM,O1) = DIRECT_FLUX
                  MEAN_STOKES(UTA,IBEAM,O1,WDIR) =
     *                  MEAN_STOKES(UTA,IBEAM,O1,WDIR) + DIRECT_MEANI
                  FLUX_STOKES(UTA,IBEAM,O1,WDIR) =
     *                  FLUX_STOKES(UTA,IBEAM,O1,WDIR) + DIRECT_FLUX
                ENDDO
              ENDIF

C  For the on-grid values

             ELSE
              N = UTAU_LEVEL_MASK_DN(UTA)
              IF ( N .LE. LAYER_PIS_CUTOFF(IBEAM) ) THEN
                IF ( N .EQ. 0 ) THEN
                  DIRECT_TRANS = ONE
                ELSE
                  DIRECT_TRANS = 
     &               INITIAL_TRANS(N,IBEAM)*T_DELT_MUBAR(N,IBEAM)
                ENDIF
                DO O1 = 1, NSTOKES
                  FTRANS = FLUXVEC(O1) * DIRECT_TRANS
                  DIRECT_FLUX  = FTRANS * LOCAL_CSZA(N,IBEAM)
                  DIRECT_MEANI = FTRANS / PI4
                  MEAN_DIRECT(UTA,IBEAM,O1) = DIRECT_MEANI
                  FLUX_DIRECT(UTA,IBEAM,O1) = DIRECT_FLUX
                  MEAN_STOKES(UTA,IBEAM,O1,WDIR) =
     *                  MEAN_STOKES(UTA,IBEAM,O1,WDIR) + DIRECT_MEANI
                  FLUX_STOKES(UTA,IBEAM,O1,WDIR) =
     *                  FLUX_STOKES(UTA,IBEAM,O1,WDIR) + DIRECT_FLUX
                ENDDO
              ENDIF
             ENDIF
           ENDDO

C  Finish downwelling direct contribution

          ENDIF

C  Finish MV output

        ENDIF   

C  Continuation point for avoiding direct beam calculation

 455    CONTINUE

C  Finish directional loop

      ENDDO

C  Finish

      RETURN
      END
 
C

      SUBROUTINE QUADINTENS_LEVEL_UP
     I    ( NLEVEL, UTA, FLUX_MULTIPLIER,
     O      QSTOKES_F )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of main model variables (input to this module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_THERMALSUP.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'

C  Subroutine arguments
C  --------------------

C  inputs

      INTEGER          NLEVEL, UTA
      DOUBLE PRECISION FLUX_MULTIPLIER

C  output

      DOUBLE PRECISION QSTOKES_F 
     &       ( MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )  

C  local variables
C  ---------------

      INTEGER          N, NL, I, I1, O1, K, KO1, K0, K1, K2
      DOUBLE PRECISION SPAR, SHOM, HOM1, HOM2, SHOM_R
      DOUBLE PRECISION SHOM_CR, HOM1CR, HOM2CR, TPROP, THELP
      DOUBLE PRECISION LXR, MXR, LXR_CR, LXR_CI, MXR_CR, MXR_CI

C  For those optical depths at layer boundaries
C  --------------------------------------------

C  This depends on the level mask - if this is 0 to NLAYERS - 1, then we are
C  looking at the intensity at the top of these layers. The
C  case where the level mask = NLAYERS is the upwelling intensity
C  at the bottom of the atmosphere (treated separately).

      NL = NLEVEL
      N  = NL + 1

C  Lowest level contributions
C  ==========================

C  Lowest level, thermal transmittance only

      IF ( NL .EQ. NLAYERS .and. DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          THELP = BOA_THTONLY_SOURCE(I)
          QSTOKES_F(UTA,I,O1) = FLUX_MULTIPLIER * THELP
        ENDDO
        RETURN
      ENDIF

C  For the lowest level, scattering solution
C  -----------------------------------------

      IF ( NL .EQ. NLAYERS ) THEN

C  Start main loops

        KO1 = K_REAL(NL) + 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES

C  real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(NL)
              LXR = LCON(K,NL)*SOLA_XPOS(I1,O1,K,NL)
              MXR = MCON(K,NL)*SOLB_XNEG(I1,O1,K,NL)
              HOM1 = LXR * T_DELT_EIGEN(K,NL)
              HOM2 = MXR
              SHOM_R = SHOM_R + HOM1 + HOM2
            ENDDO

C  complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(NL)
              K0 = 2*K-2
              K1 = KO1 + K0
              K2 = K1  + 1
              LXR_CR =  LCON(K1,NL) * SOLA_XPOS(I1,O1,K1,NL) -
     &                  LCON(K2,NL) * SOLA_XPOS(I1,O1,K2,NL)
              LXR_CI =  LCON(K1,NL) * SOLA_XPOS(I1,O1,K2,NL) +
     &                  LCON(K2,NL) * SOLA_XPOS(I1,O1,K1,NL)
              MXR_CR =  MCON(K1,NL) * SOLB_XNEG(I1,O1,K1,NL) -
     &                  MCON(K2,NL) * SOLB_XNEG(I1,O1,K2,NL)
              HOM1CR =  LXR_CR*T_DELT_EIGEN(K1,NL)
     &                - LXR_CI*T_DELT_EIGEN(K2,NL)
              HOM2CR = MXR_CR
              SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
            ENDDO

C  real part, add particular solution, complete result

            SHOM = SHOM_R + SHOM_CR
            SPAR = WLOWER(I1,O1,NL)
            QSTOKES_F(UTA,I,O1) = FLUX_MULTIPLIER * ( SPAR + SHOM )

C  Finish streams/stokes loops

          ENDDO
        ENDDO

C  End lowest level clause

      ENDIF

C  For other levels in the atmosphere
C  ==================================

C  For other levels, thermal transmittance only

      IF ( NL .NE. NLAYERS .and. DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          THELP = BOA_THTONLY_SOURCE(I)
          DO K = NLAYERS, N, -1
            TPROP = T_WUPPER(I1,K)/QUAD_STREAMS(I)
            THELP = THELP * T_DELT_DISORDS(I,K) + TPROP
          ENDDO
          QSTOKES_F(UTA,I,O1) = FLUX_MULTIPLIER * THELP
        ENDDO
        RETURN
      ENDIF

C  For other levels, scattering solution
C  -------------------------------------

      IF ( NL .NE. NLAYERS ) THEN

C  start main loops

        KO1 = K_REAL(N) + 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES

C  real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(N)
              LXR = LCON(K,N)*SOLA_XPOS(I1,O1,K,N)
              MXR = MCON(K,N)*SOLB_XNEG(I1,O1,K,N)
              HOM1 = LXR
              HOM2 = MXR * T_DELT_EIGEN(K,N)
              SHOM_R = SHOM_R + HOM1 + HOM2
            ENDDO

C  complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K-2
              K1 = KO1 + K0
              K2 = K1  + 1
              LXR_CR =  LCON(K1,N) * SOLA_XPOS(I1,O1,K1,N) -
     &                  LCON(K2,N) * SOLA_XPOS(I1,O1,K2,N)
              MXR_CR =  MCON(K1,N) * SOLB_XNEG(I1,O1,K1,N) -
     &                  MCON(K2,N) * SOLB_XNEG(I1,O1,K2,N)
              MXR_CI =  MCON(K1,N) * SOLB_XNEG(I1,O1,K2,N) +
     &                  MCON(K2,N) * SOLB_XNEG(I1,O1,K1,N)
              HOM1CR = LXR_CR
              HOM2CR =  MXR_CR * T_DELT_EIGEN(K1,N)
     &                - MXR_CI * T_DELT_EIGEN(K2,N)
              SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
            ENDDO

C  real part, add particular solution, complete result

            SHOM = SHOM_R + SHOM_CR
            SPAR = WUPPER(I1,O1,N)
            QSTOKES_F(UTA,I,O1) = FLUX_MULTIPLIER * ( SPAR + SHOM )

C  Finish streams/stokes loops

          ENDDO
        ENDDO
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE QUADINTENS_LEVEL_DN
     I    ( NLEVEL, UTA, FLUX_MULTIPLIER,
     O      QSTOKES_F )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of main model variables (input to this module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_THERMALSUP.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'

C  Subroutine arguments
C  --------------------

C  inputs

      INTEGER          NLEVEL, UTA
      DOUBLE PRECISION FLUX_MULTIPLIER

C  output

      DOUBLE PRECISION QSTOKES_F 
     &       ( MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )  

C  local variables
C  ---------------

      INTEGER          N, I, O1, K, KO1, K0, K1, K2
      DOUBLE PRECISION SPAR, SHOM, HOM1, HOM2, SHOM_R
      DOUBLE PRECISION SHOM_CR, HOM1CR, HOM2CR, TPROP, THELP
      DOUBLE PRECISION LXR, MXR, LXR_CR, LXR_CI, MXR_CR

C  For those optical depths at layer boundaries
C  --------------------------------------------

      N = NLEVEL

C  Downwelling radiation at TOA ( or N = 0 ) is zero
C      A check has been performed on this.

      IF ( NLEVEL .EQ. 0 ) THEN
        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES
            QSTOKES_F(UTA,I,O1) = ZERO
          ENDDO
        ENDDO
        RETURN
      ENDIF

C  Other levels, Thermal transmittance-only solution

      IF ( NLEVEL .NE. 0 .and. DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          THELP = ZERO
          DO K = 1, N
            TPROP = T_WLOWER(I,K)/QUAD_STREAMS(I)
            THELP = THELP * T_DELT_DISORDS(I,K) + TPROP
          ENDDO
          QSTOKES_F(UTA,I,O1) = FLUX_MULTIPLIER * THELP
        ENDDO
        RETURN
      ENDIF

C  Other levels, scattering solution
C  ---------------------------------

      IF ( NLEVEL .NE. 0 ) THEN

C  start main loops

        KO1 = K_REAL(N) + 1
        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES

C  real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(N)
              LXR  = LCON(K,N)*SOLA_XPOS(I,O1,K,N)
              MXR  = MCON(K,N)*SOLB_XNEG(I,O1,K,N)
              HOM1 = LXR * T_DELT_EIGEN(K,N)
              HOM2 = MXR
              SHOM_R = SHOM_R + HOM1 + HOM2
            ENDDO

C  complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K-2
              K1 = KO1 + K0
              K2 = K1  + 1
              LXR_CR =  LCON(K1,N) * SOLA_XPOS(I,O1,K1,N) -
     &                  LCON(K2,N) * SOLA_XPOS(I,O1,K2,N)
              LXR_CI =  LCON(K1,N) * SOLA_XPOS(I,O1,K2,N) +
     &                  LCON(K2,N) * SOLA_XPOS(I,O1,K1,N)
              MXR_CR =  MCON(K1,N) * SOLB_XNEG(I,O1,K1,N) -
     &                  MCON(K2,N) * SOLB_XNEG(I,O1,K2,N)
              HOM1CR = LXR_CR*T_DELT_EIGEN(K1,N)
     &                -LXR_CI*T_DELT_EIGEN(K2,N)
              HOM2CR = MXR_CR
              SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
            ENDDO

C  real part, add particular solution, complete result

            SHOM = SHOM_R + SHOM_CR
            SPAR = WLOWER(I,O1,N)
            QSTOKES_F(UTA,I,O1) = FLUX_MULTIPLIER * ( SPAR + SHOM )

C  Finish streams/stokes loops

          ENDDO
        ENDDO

      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE QUADINTENS_OFFGRID_UP
     I  ( N, UTA, UT, IBEAM, FLUX_MULTIPLIER, DO_INCLUDE_THERMEMISS,
     O    QSTOKES_F )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of main model variables (input to this module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_THERMALSUP.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'

C  Subroutine arguments
C  --------------------

C  inputs

      INTEGER          N, UTA, UT, IBEAM
      DOUBLE PRECISION FLUX_MULTIPLIER
      LOGICAL          DO_INCLUDE_THERMEMISS

C  output

      DOUBLE PRECISION QSTOKES_F 
     &       ( MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )  

C  local variables
C  ---------------

      INTEGER          I, I1, O1, K, KO1, K0, K1, K2
      DOUBLE PRECISION SPAR, SHOM, HOM1, HOM2, SHOM_R
      DOUBLE PRECISION SHOM_CR, HOM1CR, HOM2CR, TPROP, THELP
      DOUBLE PRECISION LXR, MXR, LXR_CR, LXR_CI, MXR_CR, MXR_CI

C  Thermal Transmittance only
c  --------------------------

      IF ( DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          THELP = BOA_THTONLY_SOURCE(I)
          DO K = NLAYERS, N + 1, -1
            TPROP = T_WUPPER(I1,K) / QUAD_STREAMS(I)
            THELP = THELP * T_DELT_DISORDS(I,K) + TPROP
          ENDDO
          TPROP =  UT_T_PARTIC(I1,UT) / QUAD_STREAMS(I)
          THELP = THELP * T_DISORDS_UTUP(I,UT) + TPROP
          QSTOKES_F(UTA,I,O1) = FLUX_MULTIPLIER * THELP
        ENDDO
        RETURN
      ENDIF

C  For those optical depths at off-grid levels
C  -------------------------------------------

C  Homogeneous solution

      KO1 = K_REAL(N) + 1
      DO I = 1, NSTREAMS
        I1 = I + NSTREAMS
        DO O1 = 1, NSTOKES

C  real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, K_REAL(N)
            LXR  = LCON(K,N)*SOLA_XPOS(I1,O1,K,N)
            MXR  = MCON(K,N)*SOLB_XNEG(I1,O1,K,N)
            HOM1 = LXR * T_UTDN_EIGEN(K,UT)
            HOM2 = MXR * T_UTUP_EIGEN(K,UT)
            SHOM_R = SHOM_R + HOM1 + HOM2
          ENDDO

C  complex homogeneous solutions

          SHOM_CR = ZERO
          DO K = 1, K_COMPLEX(N)
            K0 = 2*K-2
            K1 = KO1 + K0
            K2 = K1  + 1
            LXR_CR =  LCON(K1,N) * SOLA_XPOS(I1,O1,K1,N) -
     &                LCON(K2,N) * SOLA_XPOS(I1,O1,K2,N)
            LXR_CI =  LCON(K1,N) * SOLA_XPOS(I1,O1,K2,N) +
     &                LCON(K2,N) * SOLA_XPOS(I1,O1,K1,N)
            MXR_CR =  MCON(K1,N) * SOLB_XNEG(I1,O1,K1,N) -
     &                MCON(K2,N) * SOLB_XNEG(I1,O1,K2,N)
            MXR_CI =  MCON(K1,N) * SOLB_XNEG(I1,O1,K2,N) +
     &                MCON(K2,N) * SOLB_XNEG(I1,O1,K1,N)
            HOM1CR =    LXR_CR * T_UTDN_EIGEN(K1,UT)
     &                - LXR_CI * T_UTDN_EIGEN(K2,UT)
            HOM2CR =    MXR_CR * T_UTUP_EIGEN(K1,UT)
     &                - MXR_CI * T_UTUP_EIGEN(K2,UT)
            SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
          ENDDO

C  real part

          SHOM = SHOM_R + SHOM_CR
          QSTOKES_F(UTA,I,O1) = FLUX_MULTIPLIER * SHOM

C  Finish streams/stokes loops

        ENDDO      
      ENDDO

C  Add the thermal solution  (if flagged)

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          SPAR = UT_T_PARTIC(I1, UT)
           QSTOKES_F(UTA,I,O1) = 
     &       QSTOKES_F(UTA,I,O1) + FLUX_MULTIPLIER * SPAR
        ENDDO
      ENDIF

C  Finished if no solar terms

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

C  Add the solar particular solutions
C     Green function solution is still a Placeholder....!

      IF ( DO_CLASSICAL_SOLUTION ) THEN
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES
            SPAR = WUPPER(I1,O1,N) * T_UTDN_MUBAR(UT,IBEAM)
            QSTOKES_F(UTA,I,O1) = 
     &        QSTOKES_F(UTA,I,O1) + FLUX_MULTIPLIER * SPAR
          ENDDO
        ENDDO
      ELSE
C                P L A C E H O L D E R
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE QUADINTENS_OFFGRID_DN
     I  ( N, UTA, UT, IBEAM, FLUX_MULTIPLIER, DO_INCLUDE_THERMEMISS,
     O    QSTOKES_F )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of main model variables (input to this module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_THERMALSUP.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  Subroutine arguments
C  --------------------

C  inputs

      INTEGER          N, UTA, UT, IBEAM
      DOUBLE PRECISION FLUX_MULTIPLIER
      LOGICAL          DO_INCLUDE_THERMEMISS

C  output

      DOUBLE PRECISION QSTOKES_F 
     &       ( MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )  

C  local variables
C  ---------------

      INTEGER          I, O1, K, KO1, K0, K1, K2
      DOUBLE PRECISION SPAR, SHOM, HOM1, HOM2, SHOM_R
      DOUBLE PRECISION SHOM_CR, HOM1CR, HOM2CR, TPROP, THELP
      DOUBLE PRECISION LXR, MXR, LXR_CR, LXR_CI, MXR_CR, MXR_CI

C  Thermal Transmittance only
c  --------------------------

      IF ( DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          THELP = ZERO
          DO K = 1, N - 1
            TPROP = T_WLOWER(I,K) / QUAD_STREAMS(I)
            THELP = THELP*T_DELT_DISORDS(I,K) + TPROP
          ENDDO
          TPROP = UT_T_PARTIC(I,UT) / QUAD_STREAMS(I)
          THELP = THELP * T_DISORDS_UTDN(I,UT) + TPROP
          QSTOKES_F(UTA,I,O1) = FLUX_MULTIPLIER * THELP
        ENDDO
        RETURN
      ENDIF

C  For those optical depths at off-grid levels
C  -------------------------------------------

C  Homogeneous solution

      KO1 = K_REAL(N) + 1
      DO I = 1, NSTREAMS
        DO O1 = 1, NSTOKES

C  real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, K_REAL(N)
            LXR  = LCON(K,N)*SOLA_XPOS(I,O1,K,N)
            MXR  = MCON(K,N)*SOLB_XNEG(I,O1,K,N)
            HOM1 = LXR * T_UTDN_EIGEN(K,UT)
            HOM2 = MXR * T_UTUP_EIGEN(K,UT)
            SHOM_R = SHOM_R + HOM1 + HOM2
          ENDDO

C  complex homogeneous solutions

          SHOM_CR = ZERO
          DO K = 1, K_COMPLEX(N)
            K0 = 2*K-2
            K1 = KO1 + K0
            K2 = K1  + 1
            LXR_CR =  LCON(K1,N) * SOLA_XPOS(I,O1,K1,N) -
     &                LCON(K2,N) * SOLA_XPOS(I,O1,K2,N)
            LXR_CI =  LCON(K1,N) * SOLA_XPOS(I,O1,K2,N) +
     &                LCON(K2,N) * SOLA_XPOS(I,O1,K1,N)
            MXR_CR =  MCON(K1,N) * SOLB_XNEG(I,O1,K1,N) -
     &                MCON(K2,N) * SOLB_XNEG(I,O1,K2,N)
            MXR_CI =  MCON(K1,N) * SOLB_XNEG(I,O1,K2,N) +
     &                MCON(K2,N) * SOLB_XNEG(I,O1,K1,N)
            HOM1CR =    LXR_CR * T_UTDN_EIGEN(K1,UT)
     &                - LXR_CI * T_UTDN_EIGEN(K2,UT)
            HOM2CR =    MXR_CR * T_UTUP_EIGEN(K1,UT)
     &                - MXR_CI * T_UTUP_EIGEN(K2,UT)
            SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
          ENDDO

C  real part

          SHOM = SHOM_R + SHOM_CR
          QSTOKES_F(UTA,I,O1) = FLUX_MULTIPLIER * SHOM

C  Finish streams/stokes loops

        ENDDO      
      ENDDO

C  Add the thermal solution  (if flagged)

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          SPAR = UT_T_PARTIC(I,UT)
           QSTOKES_F(UTA,I,O1) = 
     &       QSTOKES_F(UTA,I,O1) + FLUX_MULTIPLIER * SPAR
        ENDDO
      ENDIF

C  Finished if no solar terms

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

C  Add the solar particular solutions
C     Green function solution is still a Placeholder....!

      IF ( DO_CLASSICAL_SOLUTION ) THEN
        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES
            SPAR = WUPPER(I,O1,N) * T_UTDN_MUBAR(UT,IBEAM)
            QSTOKES_F(UTA,I,O1) = 
     &        QSTOKES_F(UTA,I,O1) + FLUX_MULTIPLIER * SPAR
          ENDDO
        ENDDO
      ELSE
C                P L A C E H O L D E R
      ENDIF

C  Finish

      RETURN
      END

