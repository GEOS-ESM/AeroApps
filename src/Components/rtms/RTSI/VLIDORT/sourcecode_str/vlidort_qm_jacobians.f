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
C #       VLIDORT_L_INTEGRATED_OUTPUT (master)                  #
C #                                                             #
C #            QUADATMOSWF_LEVEL_UP                             #
C #            QUADATMOSWF_LEVEL_DN                             #
C #            QUADATMOSWF_OFFGRID_UP                           #
C #            QUADATMOSWF_OFFGRID_DN                           #
C #                                                             #
C #       VLIDORT_LS_INTEGRATED_OUTPUT (master)                 #
C #                                                             #
C #            QUADSURFACEWF_LEVEL_UP                           #
C #            QUADSURFACEWF_LEVEL_DN                           #
C #            QUADSURFACEWF_OFFGRID_UP                         #
C #            QUADSURFACEWF_OFFGRID_DN                         #
C #                                                             #
C ###############################################################

      SUBROUTINE VLIDORT_L_INTEGRATED_OUTPUT
     &  ( DO_INCLUDE_MVOUTPUT,
     I    DO_INCLUDE_DIRECTBEAM,
     I    DO_INCLUDE_THERMEMISS,
     I    FLUX_MULTIPLIER,
     I    IBEAM, NV, NV_PARAMETERS )

C  Quadrature output at offgrid or ongrid optical depths
C  ( Required if mean-value calculations are to be done)
C    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include file of setup variables (input)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  include file of linearized setup variables (input)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'

C  include file of result variables (module output stored here)

      INCLUDE '../includes/VLIDORT_L_RESULTS.VARS'

C  Subroutine arguments
C  --------------------

C  Flags

      LOGICAL          DO_INCLUDE_THERMEMISS
      LOGICAL          DO_INCLUDE_DIRECTBEAM
      LOGICAL          DO_INCLUDE_MVOUTPUT

C  Linearization control

      INTEGER          NV, NV_PARAMETERS

C  Beam control

      INTEGER          IBEAM
      DOUBLE PRECISION FLUX_MULTIPLIER

C  local variables
C  ---------------

C  Local quadrature output (for debug)

      DOUBLE PRECISION QATMOSWF_F 
     &     ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

C  Help variables

      INTEGER          I, IDIR, WDIR, UTA, UT, Q, N, O1, NLEVEL
      DOUBLE PRECISION SMI, SFX, FTRANS
      DOUBLE PRECISION L_TRANS, L_DIRECT_FLUX, L_DIRECT_MEANI

C  direction loop

      DO IDIR = 1, N_DIRECTIONS
        WDIR = WHICH_DIRECTIONS(IDIR)

C  Upwelling Jacobian output at Quadrature angles

        IF ( WDIR .EQ. UPIDX ) THEN
          DO UTA = 1, N_USER_LEVELS
            NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
              UT = PARTLAYERS_OUTINDEX(UTA)
              N  = PARTLAYERS_LAYERIDX(UT)
              CALL QUADATMOSWF_OFFGRID_UP
     &          ( IBEAM, UTA, UT, N, NV, NV_PARAMETERS,
     &            DO_INCLUDE_THERMEMISS, FLUX_MULTIPLIER, QATMOSWF_F )
            ELSE
              CALL QUADATMOSWF_LEVEL_UP
     &          ( UTA, NLEVEL, NV, NV_PARAMETERS,
     &            FLUX_MULTIPLIER, QATMOSWF_F )
            ENDIF
          ENDDO
        ENDIF

C  Downwelling Jacobian output at Quadrature angles

        IF ( WDIR .EQ. DNIDX ) THEN
          DO UTA = 1, N_USER_LEVELS
            NLEVEL = UTAU_LEVEL_MASK_DN(UTA)
            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
              UT = PARTLAYERS_OUTINDEX(UTA)
              N  = PARTLAYERS_LAYERIDX(UT)
              CALL QUADATMOSWF_OFFGRID_DN
     &          ( IBEAM, UTA, UT, N, NV, NV_PARAMETERS,
     &            DO_INCLUDE_THERMEMISS, FLUX_MULTIPLIER, QATMOSWF_F )
            ELSE
              CALL QUADATMOSWF_LEVEL_DN
     &          ( UTA, NLEVEL, NV, NV_PARAMETERS,
     &            FLUX_MULTIPLIER, QATMOSWF_F )
            ENDIF
          ENDDO
        ENDIF

C  Mean Intensity and Flux  output
C  -------------------------------

        IF ( DO_INCLUDE_MVOUTPUT ) THEN

C  Diffuse term integrated output

          DO Q = 1, NV_PARAMETERS
           DO UTA = 1, N_USER_LEVELS
            DO O1 = 1, NSTOKES
              SMI = ZERO
              SFX = ZERO
              DO I = 1, NSTREAMS
                SMI = SMI + QUAD_WEIGHTS(I) * QATMOSWF_F(Q,UTA,I,O1)
                SFX = SFX + QUAD_STRMWTS(I) * QATMOSWF_F(Q,UTA,I,O1)
              ENDDO
              MINT_ATMOSWF(Q,NV,UTA,IBEAM,O1,WDIR) = SMI * HALF
              FLUX_ATMOSWF(Q,NV,UTA,IBEAM,O1,WDIR) = SFX * PI2
            ENDDO
           ENDDO
          ENDDO

C  nothing further to do if no solar sources

          IF ( .NOT. DO_INCLUDE_DIRECTBEAM ) GO TO 455

C  For the downward direction, add the direct beam contributions

          IF ( WDIR .EQ. DNIDX ) THEN

C  loop over all the output optical depths

           DO UTA = 1, N_USER_LEVELS

C  For the offgrid values

            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
              UT = PARTLAYERS_OUTINDEX(UTA)
              N  = PARTLAYERS_LAYERIDX(UT)

C  For the offgrid values.......
C     .....Only contributions for layers above the PI cutoff
C    L_INITIAL_TRANS is a logarithmic derivative

              IF ( N .LE. LAYER_PIS_CUTOFF(IBEAM) ) THEN
                IF ( NV.LE.N .OR.NV.EQ.0 ) THEN
                 DO Q = 1, NV_PARAMETERS
                  L_TRANS = L_T_UTDN_MUBAR(UT,NV,IBEAM,Q) +
     &           L_INITIAL_TRANS(N,NV,IBEAM,Q) * T_UTDN_MUBAR(UT,IBEAM)
                  L_TRANS = L_TRANS * INITIAL_TRANS(N,IBEAM)
                  DO O1 = 1, NSTOKES
                   FTRANS = FLUXVEC(O1) * L_TRANS
                   L_DIRECT_MEANI = FTRANS / PI4
                   L_DIRECT_FLUX  = FTRANS * LOCAL_CSZA(N,IBEAM)
                   MINT_ATMOSWF_DIRECT(Q,NV,UTA,IBEAM,O1)=L_DIRECT_MEANI
                   FLUX_ATMOSWF_DIRECT(Q,NV,UTA,IBEAM,O1)=L_DIRECT_FLUX
                   MINT_ATMOSWF(Q,NV,UTA,IBEAM,O1,WDIR) =
     *             MINT_ATMOSWF(Q,NV,UTA,IBEAM,O1,WDIR) + L_DIRECT_MEANI
                   FLUX_ATMOSWF(Q,NV,UTA,IBEAM,O1,WDIR) =
     *             FLUX_ATMOSWF(Q,NV,UTA,IBEAM,O1,WDIR) + L_DIRECT_FLUX
                  ENDDO
                 ENDDO
                ENDIF
              ENDIF

C  For the on-grid values
C    L_INITIAL_TRANS is a logarithmic derivative

             ELSE
              N = UTAU_LEVEL_MASK_DN(UTA)
              IF ( N .LE. LAYER_PIS_CUTOFF(IBEAM) ) THEN
               IF ( N.GT.0 ) THEN
                IF ( NV.LE.N .OR.NV.EQ.0 ) THEN
                 DO Q = 1, NV_PARAMETERS
                  L_TRANS = L_T_DELT_MUBAR(N,NV,IBEAM,Q) +
     &             L_INITIAL_TRANS(N,NV,IBEAM,Q) * T_DELT_MUBAR(N,IBEAM)
                  L_TRANS = L_TRANS * INITIAL_TRANS(N,IBEAM)
                  DO O1 = 1, NSTOKES
                   FTRANS = FLUXVEC(O1) * L_TRANS
                   L_DIRECT_MEANI = FTRANS / PI4
                   L_DIRECT_FLUX  = FTRANS * LOCAL_CSZA(N,IBEAM)
                   MINT_ATMOSWF_DIRECT(Q,NV,UTA,IBEAM,O1)=L_DIRECT_MEANI
                   FLUX_ATMOSWF_DIRECT(Q,NV,UTA,IBEAM,O1)=L_DIRECT_FLUX
                   MINT_ATMOSWF(Q,NV,UTA,IBEAM,O1,WDIR) =
     *             MINT_ATMOSWF(Q,NV,UTA,IBEAM,O1,WDIR) + L_DIRECT_MEANI
                   FLUX_ATMOSWF(Q,NV,UTA,IBEAM,O1,WDIR) =
     *             FLUX_ATMOSWF(Q,NV,UTA,IBEAM,O1,WDIR) + L_DIRECT_FLUX
                  ENDDO
                 ENDDO
                ENDIF
               ENDIF
              ENDIF
             ENDIF

C  End UTA loop

           ENDDO

C  Finish downwelling direct contribution

          ENDIF

C  Continuation point for avoiding direct beam calculation

 455      CONTINUE

C  Finish MV output

        ENDIF   

C  end direction loop

      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE QUADATMOSWF_LEVEL_UP
     &      ( UTA, NL, NV, NV_PARAMETERS,
     &        FLUX_MULTIPLIER, QATMOSWF_F )

C  Upwelling weighting function Fourier components at level boundary NL
C  Quadrature angles only

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup, solution variables (input)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_THERMALSUP.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'

C  include files of linearized setup and solution variables (input)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_THERMALSUP.VARS'

C  Subroutine arguments
C  --------------------

C  Inputs

      INTEGER          NV, NV_PARAMETERS
      INTEGER          UTA, NL
      DOUBLE PRECISION FLUX_MULTIPLIER

C  output

      DOUBLE PRECISION QATMOSWF_F 
     &     ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

C  local variables
C  ---------------

      INTEGER          N, I, I1, Q, O1, K, KO1, K0, K1, K2, LAY
      DOUBLE PRECISION SPAR, SHOM_R, SHOM_CR, SHOM
      DOUBLE PRECISION HOM1, HOM2, HOM1CR, HOM2CR, HOM3CR
      DOUBLE PRECISION NXR, PXR, NXR1, NXR2, PXR1, PXR2
      DOUBLE PRECISION LXR, MXR, LXR1, MXR1, LXR2, MXR2
      DOUBLE PRECISION LLXR, MLXR, LLXR1, MLXR1, LLXR2, MLXR2
      DOUBLE PRECISION TPROP, THELP, L_TPROP, L_THELP, FM

C  homogeneous and particular solution contributions SHOM and SPAR

C  This depends on the level mask - if this is 0 to NLAYERS - 1, then we are
C  looking at the perturbation field at the top of these layers. The
C  case where the level mask = NLAYERS is the upwelling perturbed fields
C  at the bottom of the atmosphere (treated separately).

      N  = NL + 1
      FM = FLUX_MULTIPLIER

C  For the lowest level
C  ====================

C  Thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY .and. NL.EQ.NLAYERS  ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          DO Q = 1, NV_PARAMETERS
            L_THELP = L_BOA_THTONLY_SOURCE(I,Q)
            QATMOSWF_F(Q,UTA,I,O1) = FM * L_THELP
          ENDDO
        ENDDO
        RETURN
      ENDIF

C  For the lowest level, scattering solution
C  -----------------------------------------

      IF ( NL .EQ. NLAYERS ) THEN

C  If this is also the layer that is varying, extra contributions

        IF ( NV .EQ. NL .OR. NV. EQ. 0 ) THEN

C  Stokes and streams loops

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES

C  Parameter loop

              DO Q = 1, NV_PARAMETERS

C  Real homogeneous solutions

               SHOM_R = ZERO
               DO K = 1, K_REAL(NL)

                LXR  = LCON(K,NL)   *   SOLA_XPOS(I1,O1,K,NL)
                LLXR = LCON(K,NL)   * L_SOLA_XPOS(I1,O1,K,NL,Q)
                MLXR = MCON(K,NL)   * L_SOLB_XNEG(I1,O1,K,NL,Q)
                NXR  = NCON(K,NL,Q) *   SOLA_XPOS(I1,O1,K,NL)
                PXR  = PCON(K,NL,Q) *   SOLB_XNEG(I1,O1,K,NL)

                HOM1 = ( LLXR + NXR )  *   T_DELT_EIGEN(K,NL)
     &                       +  LXR    * L_T_DELT_EIGEN(K,NL,Q)
                HOM2 = PXR + MLXR
                SHOM_R = SHOM_R + HOM1 + HOM2 
               ENDDO

C  Complex homogeneous solutions

               SHOM_CR = ZERO
               KO1 = K_REAL(NL) + 1
               DO K = 1, K_COMPLEX(NL)
                K0 = 2 * K - 2
                K1 = KO1 + K0
                K2 = K1  + 1

                NXR1  =   NCON(K1,NL,Q) *   SOLA_XPOS(I1,O1,K1,NL)
     &                  - NCON(K2,NL,Q) *   SOLA_XPOS(I1,O1,K2,NL)
                NXR2  =   NCON(K1,NL,Q) *   SOLA_XPOS(I1,O1,K2,NL)
     &                  + NCON(K2,NL,Q) *   SOLA_XPOS(I1,O1,K1,NL)
                PXR1  =   PCON(K1,NL,Q) *   SOLB_XNEG(I1,O1,K1,NL)
     &                  - PCON(K2,NL,Q) *   SOLB_XNEG(I1,O1,K2,NL)

                LXR1  =   LCON(K1,NL) *   SOLA_XPOS(I1,O1,K1,NL)
     &                  - LCON(K2,NL) *   SOLA_XPOS(I1,O1,K2,NL)
                LXR2  =   LCON(K1,NL) *   SOLA_XPOS(I1,O1,K2,NL)
     &                  + LCON(K2,NL) *   SOLA_XPOS(I1,O1,K1,NL)

                LLXR1  =   LCON(K1,NL) * L_SOLA_XPOS(I1,O1,K1,NL,Q)
     &                   - LCON(K2,NL) * L_SOLA_XPOS(I1,O1,K2,NL,Q)
                LLXR2  =   LCON(K1,NL) * L_SOLA_XPOS(I1,O1,K2,NL,Q)
     &                   + LCON(K2,NL) * L_SOLA_XPOS(I1,O1,K1,NL,Q)
                MLXR1  =   MCON(K1,NL) * L_SOLB_XNEG(I1,O1,K1,NL,Q)
     &                   - MCON(K2,NL) * L_SOLB_XNEG(I1,O1,K2,NL,Q)

                HOM1CR =   ( NXR1 + LLXR1 ) *   T_DELT_EIGEN(K1,NL)
     &                   - ( NXR2 + LLXR2 ) *   T_DELT_EIGEN(K2,NL)
                HOM2CR =             LXR1   * L_T_DELT_EIGEN(K1,NL,Q)
     &                             - LXR2   * L_T_DELT_EIGEN(K2,NL,Q)
                HOM3CR = PXR1 + MLXR1
                SHOM_CR = SHOM_CR + HOM1CR + HOM2CR + HOM3CR 
               ENDDO

C  real part, add particular solution, complete result

               SHOM = SHOM_R + SHOM_CR
               SPAR = L_WLOWER(I1,O1,NL,Q)
               QATMOSWF_F(Q,UTA,I,O1) = FM * ( SPAR + SHOM )

C  Finish Q, I and O1 loops

              ENDDO
            ENDDO
          ENDDO

C  non-varying lowest layer

        ELSE IF ( NV.LT.NL .AND. NV.NE.0 ) THEN

C  Stokes and streams loops

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES

C  Parameter loop

              DO Q = 1, NV_PARAMETERS

C  Real homogeneous solutions

               SHOM_R = ZERO
               DO K = 1, K_REAL(NL)
                NXR = NCON(K,NL,Q) *   SOLA_XPOS(I1,O1,K,NL)
                PXR = PCON(K,NL,Q) *   SOLB_XNEG(I1,O1,K,NL)
                HOM1 = NXR * T_DELT_EIGEN(K,NL)
                HOM2 = PXR
                SHOM_R = SHOM_R + HOM1 + HOM2 
               ENDDO

C  Complex homogeneous solutions

               SHOM_CR = ZERO
               KO1 = K_REAL(NL) + 1
               DO K = 1, K_COMPLEX(NL)
                K0 = 2 * K - 2
                K1 = KO1 + K0
                K2 = K1  + 1
                NXR1  =   NCON(K1,NL,Q) *   SOLA_XPOS(I1,O1,K1,NL)
     &                  - NCON(K2,NL,Q) *   SOLA_XPOS(I1,O1,K2,NL)
                NXR2  =   NCON(K1,NL,Q) *   SOLA_XPOS(I1,O1,K2,NL)
     &                  + NCON(K2,NL,Q) *   SOLA_XPOS(I1,O1,K1,NL)
                PXR1  =   PCON(K1,NL,Q) *   SOLB_XNEG(I1,O1,K1,NL)
     &                  - PCON(K2,NL,Q) *   SOLB_XNEG(I1,O1,K2,NL)
                HOM1CR =   NXR1 * T_DELT_EIGEN(K1,NL)
     &                   - NXR2 * T_DELT_EIGEN(K2,NL)
                HOM2CR = PXR1
                SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
               ENDDO

C  real part, add particular solution, complete result

               SHOM = SHOM_R + SHOM_CR
               SPAR = L_WLOWER(I1,O1,NL,Q)
               QATMOSWF_F(Q,UTA,I,O1) = FM * ( SPAR + SHOM )

C  Finish Q, I and O1 loops

              ENDDO
            ENDDO
          ENDDO

        ENDIF

C  End lowest level clause

      ENDIF

C  For other levels in the atmosphere
C  ==================================

C  For other levels, thermal transmittance only

      IF ( NL .NE. NLAYERS .and. DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO Q = 1, NV_PARAMETERS
            L_THELP = L_BOA_THTONLY_SOURCE(I,Q)
            THELP   =   BOA_THTONLY_SOURCE(I)
            DO LAY = NLAYERS, N, -1
              L_THELP = L_THELP * T_DELT_DISORDS(I,LAY) 
              IF ( LAY.EQ.NV .OR. NV.EQ.0 ) THEN
                L_TPROP = L_T_WUPPER(I1,LAY,Q) / QUAD_STREAMS(I)
                L_THELP = L_THELP + L_TPROP
     &                    + THELP * L_T_DELT_DISORDS(I,LAY,Q) 
              ENDIF
              TPROP = T_WUPPER(I1,LAY) / QUAD_STREAMS(I)
              THELP = THELP * T_DELT_DISORDS(I,LAY) + TPROP
            ENDDO
            QATMOSWF_F(Q,UTA,I,O1) = FM * L_THELP
          ENDDO
        ENDDO
        RETURN
      ENDIF

C  For other levels, scattering solution
C  -------------------------------------

      IF ( NL .NE. NLAYERS ) THEN

C  If this is also the layer that is varying, extra contributions

        IF ( NV .EQ. N .OR. NV. EQ. 0 ) THEN

C  Stokes and streams loops

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES

C  Parameter loop

              DO Q = 1, NV_PARAMETERS

C  Real homogeneous solutions

               SHOM_R = ZERO
               DO K = 1, K_REAL(N)
                LXR  = LCON(K,N)   *   SOLA_XPOS(I1,O1,K,N)
                MXR  = MCON(K,N)   *   SOLB_XNEG(I1,O1,K,N)
                LLXR = LCON(K,N)   * L_SOLA_XPOS(I1,O1,K,N,Q)
                MLXR = MCON(K,N)   * L_SOLB_XNEG(I1,O1,K,N,Q)
                NXR  = NCON(K,N,Q) *   SOLA_XPOS(I1,O1,K,N)
                PXR  = PCON(K,N,Q) *   SOLB_XNEG(I1,O1,K,N)
                HOM1 =   NXR + LLXR
                HOM2 = ( PXR + MLXR ) *   T_DELT_EIGEN(K,N)
     &                       + MXR    * L_T_DELT_EIGEN(K,N,Q)
                SHOM_R = SHOM_R + HOM1 + HOM2
               ENDDO

C  Complex homogeneous solutions

               SHOM_CR = ZERO
               KO1 = K_REAL(N) + 1
               DO K = 1, K_COMPLEX(N)
                K0 = 2 * K - 2
                K1 = KO1 + K0
                K2 = K1  + 1

                NXR1  =   NCON(K1,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
     &                  - NCON(K2,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
                PXR1  =   PCON(K1,N,Q) *   SOLB_XNEG(I1,O1,K1,N)
     &                  - PCON(K2,N,Q) *   SOLB_XNEG(I1,O1,K2,N)
                PXR2  =   PCON(K1,N,Q) *   SOLB_XNEG(I1,O1,K2,N)
     &                  + PCON(K2,N,Q) *   SOLB_XNEG(I1,O1,K1,N)

                MXR1  =   MCON(K1,N) *   SOLB_XNEG(I1,O1,K1,N)
     &                  - MCON(K2,N) *   SOLB_XNEG(I1,O1,K2,N)
                MXR2  =   MCON(K1,N) *   SOLB_XNEG(I1,O1,K2,N)
     &                  + MCON(K2,N) *   SOLB_XNEG(I1,O1,K1,N)

                LLXR1  =   LCON(K1,N) * L_SOLA_XPOS(I1,O1,K1,N,Q)
     &                   - LCON(K2,N) * L_SOLA_XPOS(I1,O1,K2,N,Q)
                MLXR1  =   MCON(K1,N) * L_SOLB_XNEG(I1,O1,K1,N,Q)
     &                   - MCON(K2,N) * L_SOLB_XNEG(I1,O1,K2,N,Q)
                MLXR2  =   MCON(K1,N) * L_SOLB_XNEG(I1,O1,K2,N,Q)
     &                   + MCON(K2,N) * L_SOLB_XNEG(I1,O1,K1,N,Q)

                HOM1CR =   ( PXR1 + MLXR1 ) *   T_DELT_EIGEN(K1,N)
     &                   - ( PXR2 + MLXR2 ) *   T_DELT_EIGEN(K2,N)
                HOM2CR =             MXR1   * L_T_DELT_EIGEN(K1,N,Q)
     &                             - MXR2   * L_T_DELT_EIGEN(K2,N,Q)
                HOM3CR = NXR1 + LLXR1
                SHOM_CR = SHOM_CR + HOM1CR + HOM2CR + HOM3CR 

               ENDDO

C  real part, add particular solution, complete result

               SHOM = SHOM_R + SHOM_CR
               SPAR = L_WUPPER(I1,O1,N,Q)
               QATMOSWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER * (SPAR+SHOM)

C  Finish Q, I and O1 loops

              ENDDO
            ENDDO
          ENDDO

C  For layers N beneath or above the varying layer NV

        ELSE IF ( NV.NE.N .AND. NV.NE.0 ) THEN

C  Stokes and streams loops

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES

C  Parameter loop

              DO Q = 1, NV_PARAMETERS

C  Real homogeneous solutions

               SHOM_R = ZERO
               DO K = 1, K_REAL(N)
                NXR  = NCON(K,N,Q) *   SOLA_XPOS(I1,O1,K,N)
                PXR  = PCON(K,N,Q) *   SOLB_XNEG(I1,O1,K,N)
                HOM1 = NXR
                HOM2 = PXR * T_DELT_EIGEN(K,N)
                SHOM_R = SHOM_R + HOM1 + HOM2 
               ENDDO

C  Complex homogeneous solutions

               SHOM_CR = ZERO
               KO1 = K_REAL(N) + 1
               DO K = 1, K_COMPLEX(N)
                K0 = 2 * K - 2
                K1 = KO1 + K0
                K2 = K1  + 1
                NXR1  =   NCON(K1,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
     &                  - NCON(K2,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
                PXR1  =   PCON(K1,N,Q) *   SOLB_XNEG(I1,O1,K1,N)
     &                  - PCON(K2,N,Q) *   SOLB_XNEG(I1,O1,K2,N)
                PXR2  =   PCON(K1,N,Q) *   SOLB_XNEG(I1,O1,K2,N)
     &                  + PCON(K2,N,Q) *   SOLB_XNEG(I1,O1,K1,N)
                HOM1CR =   PXR1 * T_DELT_EIGEN(K1,N)
     &                   - PXR2 * T_DELT_EIGEN(K2,N)
                HOM2CR = NXR1
                SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
               ENDDO

C  real part, add particular solution (only if N>NV), complete result

               SHOM = SHOM_R + SHOM_CR
               SPAR = ZERO
               IF (NV.LT.N)SPAR = L_WUPPER(I1,O1,N,Q)
               QATMOSWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER * (SPAR+SHOM)

C  Finish Q, I and O1 loops

              ENDDO
            ENDDO
          ENDDO

C  End variability clauses

        ENDIF
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE QUADATMOSWF_LEVEL_DN
     &      ( UTA, NL, NV, NV_PARAMETERS,
     &        FLUX_MULTIPLIER, QATMOSWF_F )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup, solution variables (input)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_THERMALSUP.VARS'

C  include files of linearized setup and solution variables (input)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_THERMALSUP.VARS'

C  Subroutine arguments
C  --------------------

C  Inputs

      INTEGER          NV, NV_PARAMETERS
      INTEGER          UTA, NL
      DOUBLE PRECISION FLUX_MULTIPLIER

C  output

      DOUBLE PRECISION QATMOSWF_F 
     &     ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

C  local variables
C  ---------------

      INTEGER          N, I, Q, O1, K, KO1, K0, K1, K2, LAY
      DOUBLE PRECISION SPAR, SHOM_R, SHOM_CR, SHOM
      DOUBLE PRECISION HOM1, HOM2, HOM1CR, HOM2CR, HOM3CR
      DOUBLE PRECISION NXR, PXR, NXR1, NXR2, PXR1
      DOUBLE PRECISION LXR, LXR1, LXR2
      DOUBLE PRECISION LLXR, MLXR, LLXR1, MLXR1, LLXR2
      DOUBLE PRECISION TPROP, THELP, L_TPROP, L_THELP, FM

C  Downwelling weighting function at TOA ( or N = 0 ) is zero
C    Zero and return

      IF ( NL .EQ. 0 ) THEN
        DO I = 1, NSTREAMS
          DO Q = 1, NV_PARAMETERS
            DO O1 = 1, NSTOKES
              QATMOSWF_F(Q,UTA,I,O1) = ZERO
            ENDDO
          ENDDO
        ENDDO
        RETURN
      ENDIF

C  For other levels in the atmosphere
C  ----------------------------------

C  Other levels, Thermal transmittance-only solution
C  Thermal transmittance solution, build from TOA downwards
C  Scattering solution, use the Discrete Ordinate solution

      IF ( NL.NE.0 .and. DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          DO Q = 1, NV_PARAMETERS
            L_THELP = ZERO
            THELP   = ZERO
            DO LAY = 1, N
              L_THELP = L_THELP * T_DELT_DISORDS(I,LAY)
              IF ( LAY.EQ.NV .OR. NV.EQ.0 ) THEN
                L_TPROP = L_T_WLOWER(I,LAY,Q) / QUAD_STREAMS(I)
                L_THELP = L_THELP +  L_TPROP
     &                    + THELP * L_T_DELT_DISORDS(I,LAY,Q) 
              ENDIF
              TPROP = T_WLOWER(I,LAY) / QUAD_STREAMS(I)
              THELP = THELP * T_DELT_DISORDS(I,LAY) + TPROP
            ENDDO
            QATMOSWF_F(Q,UTA,I,O1) = FM * L_THELP
          ENDDO
        ENDDO
        RETURN
      ENDIF

C  Scattering solutions
C  --------------------

      IF ( NL.NE.0 ) THEN

C  Shorthand

        N = NL

C  If this is also the layer that is varying, extra contributions

        IF ( NV .EQ. N .OR. NV.EQ.0 ) THEN

C  Stokes and streams loops

          DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES

C  Parameter loop

              DO Q = 1, NV_PARAMETERS

C  Real homogeneous solutions

               SHOM_R = ZERO
               DO K = 1, K_REAL(N)
                LXR  = LCON(K,N)   *   SOLA_XPOS(I,O1,K,N)
                LLXR = LCON(K,N)   * L_SOLA_XPOS(I,O1,K,N,Q)
                MLXR = MCON(K,N)   * L_SOLB_XNEG(I,O1,K,N,Q)
                NXR  = NCON(K,N,Q) *   SOLA_XPOS(I,O1,K,N)
                PXR  = PCON(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
                HOM1 = ( NXR + LLXR ) *   T_DELT_EIGEN(K,N)
     &                       +  LXR   * L_T_DELT_EIGEN(K,N,Q)
                HOM2 = PXR + MLXR
                SHOM_R = SHOM_R + HOM1 + HOM2 
               ENDDO

C  Complex homogeneous solutions

               SHOM_CR = ZERO
               KO1 = K_REAL(N) + 1
               DO K = 1, K_COMPLEX(N)
                K0 = 2 * K - 2
                K1 = KO1 + K0
                K2 = K1  + 1

                NXR1  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K1,N)
     &                  - NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K2,N)
                NXR2  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K2,N)
     &                  + NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K1,N)
                PXR1  =   PCON(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N)
     &                  - PCON(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)

                LXR1  =   LCON(K1,N) *   SOLA_XPOS(I,O1,K1,N)
     &                  - LCON(K2,N) *   SOLA_XPOS(I,O1,K2,N)
                LXR2  =   LCON(K1,N) *   SOLA_XPOS(I,O1,K2,N)
     &                  + LCON(K2,N) *   SOLA_XPOS(I,O1,K1,N)

                LLXR1  =   LCON(K1,N) * L_SOLA_XPOS(I,O1,K1,N,Q)
     &                   - LCON(K2,N) * L_SOLA_XPOS(I,O1,K2,N,Q)
                LLXR2  =   LCON(K1,N) * L_SOLA_XPOS(I,O1,K2,N,Q)
     &                   + LCON(K2,N) * L_SOLA_XPOS(I,O1,K1,N,Q)
                MLXR1  =   MCON(K1,N) * L_SOLB_XNEG(I,O1,K1,N,Q)
     &                   - MCON(K2,N) * L_SOLB_XNEG(I,O1,K2,N,Q)

                HOM1CR =   ( NXR1 + LLXR1 ) *   T_DELT_EIGEN(K1,N)
     &                   - ( NXR2 + LLXR2 ) *   T_DELT_EIGEN(K2,N)
                HOM2CR =             LXR1   * L_T_DELT_EIGEN(K1,N,Q)
     &                             - LXR2   * L_T_DELT_EIGEN(K2,N,Q)
                HOM3CR = PXR1 + MLXR1
                SHOM_CR = SHOM_CR + HOM1CR + HOM2CR + HOM3CR 

               ENDDO

C  real part, add particular solution, complete result

               SHOM = SHOM_R + SHOM_CR
               SPAR = L_WLOWER(I,O1,N,Q)
               QATMOSWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER * (SPAR+SHOM)

C  Finish Q, I and O1 loops

              ENDDO
            ENDDO
          ENDDO

C  varying layer above or below active layer

        ELSE IF ( NV.NE.N .AND. NV.NE.0 ) THEN

C  Stokes and streams loops

          DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES

C  Parameter loop

              DO Q = 1, NV_PARAMETERS

C  Real homogeneous solutions

               SHOM_R = ZERO
               DO K = 1, K_REAL(N)
                NXR = NCON(K,N,Q) *   SOLA_XPOS(I,O1,K,N)
                PXR = PCON(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
                HOM1 = NXR * T_DELT_EIGEN(K,N)
                HOM2 = PXR
                SHOM_R = SHOM_R + HOM1 + HOM2 
               ENDDO

C  Complex homogeneous solutions

               SHOM_CR = ZERO
               KO1 = K_REAL(N) + 1
               DO K = 1, K_COMPLEX(N)
                K0 = 2 * K - 2
                K1 = KO1 + K0
                K2 = K1  + 1
                NXR1  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K1,N)
     &                  - NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K2,N)
                NXR2  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K2,N)
     &                  + NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K1,N)
                PXR1  =   PCON(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N)
     &                  - PCON(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)
                HOM1CR =   NXR1 * T_DELT_EIGEN(K1,N)
     &                   - NXR2 * T_DELT_EIGEN(K2,N)
                HOM2CR = PXR1
                SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
               ENDDO

C  real part, add particular solution (only if N > NV), complete result

               SHOM = SHOM_R + SHOM_CR
               SPAR = ZERO
               IF ( N .GT. NV ) SPAR = L_WLOWER(I,O1,N,Q)
               QATMOSWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER * (SPAR+SHOM)

C  Finish Q, I and O1 loops

              ENDDO
            ENDDO
          ENDDO

C  Finish variability clauses

        ENDIF
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE QUADATMOSWF_OFFGRID_UP
     &     ( IB, UTA, UT, N, NV, NV_PARAMETERS,
     &       DO_INCLUDE_THERMEMISS, FLUX_MULTIPLIER, QATMOSWF_F )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and multiplier variables (input)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'
      INCLUDE '../includes/VLIDORT_THERMALSUP.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'

C  include files of linearized setup/solution/multiplier variables (input)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'
      INCLUDE '../includes/VLIDORT_L_THERMALSUP.VARS'

C  Subroutine arguments
C  --------------------

C  Inputs

      INTEGER          NV, NV_PARAMETERS
      INTEGER          IB, UTA, UT, N

      DOUBLE PRECISION FLUX_MULTIPLIER
      LOGICAL          DO_INCLUDE_THERMEMISS

C  output

      DOUBLE PRECISION QATMOSWF_F 
     &     ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

C  local variables
C  ---------------

      INTEGER          I, I1, Q, O1, K, KO1, K0, K1, K2, LAY
      DOUBLE PRECISION SPAR, SHOM_R, SHOM_CR, SHOM
      DOUBLE PRECISION HOM1, HOM2, HOM1CR, HOM2CR, HOM3CR, HOM4CR
      DOUBLE PRECISION NXR, PXR, NXR1, NXR2, PXR1, PXR2
      DOUBLE PRECISION LXR, MXR, LXR1, MXR1, LXR2, MXR2
      DOUBLE PRECISION LLXR, MLXR, LLXR1, MLXR1, LLXR2, MLXR2
      DOUBLE PRECISION TPROP, THELP, L_TPROP, L_THELP, FMULT

C  short hand

      FMULT = FLUX_MULTIPLIER

C  Thermal Transmittance only
c  --------------------------

      IF ( DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO Q = 1, NV_PARAMETERS
            THELP     =   BOA_THTONLY_SOURCE(I)
            L_THELP   = L_BOA_THTONLY_SOURCE(I,Q)
            DO LAY = NLAYERS, N+1, -1
              L_THELP = L_THELP *   T_DELT_DISORDS(I,LAY)
              IF ( LAY.EQ.NV .OR. NV.EQ.0 ) THEN
                L_TPROP = L_T_WUPPER(I1,LAY,Q) / QUAD_STREAMS(I)
                L_THELP = L_THELP + L_TPROP
     &                    + THELP * L_T_DELT_DISORDS(I,LAY,Q) 
              ENDIF
              TPROP = T_WUPPER(I1,Q) / QUAD_STREAMS(I)
              THELP = THELP * T_DELT_DISORDS(I,LAY) + TPROP
            ENDDO
            L_THELP = L_THELP * T_DISORDS_UTUP(I,UT)
            IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN
              L_TPROP = L_UT_T_PARTIC(I1,UT,Q) / QUAD_STREAMS(I)
              L_THELP = L_THELP + L_TPROP
     &                  + THELP * L_T_DISORDS_UTUP(I,UT,Q) 
            ENDIF
            QATMOSWF_F(Q,UTA,I,O1) = FMULT * L_THELP
          ENDDO
        ENDDO
        RETURN
      ENDIF

C  For those optical depths at off-grid levels
C  ###########################################
      
C  Homogeneous
C  -----------

C  solution for N being the varying layer NV

      IF ( N. EQ. NV .OR. NV.EQ.0 ) THEN

C  stream directions

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES

C  Parameter loop

            DO Q = 1, NV_PARAMETERS

C  real homogeneous solutions

              SHOM_R = ZERO
              DO K = 1, K_REAL(N)
                LXR  = LCON(K,N)   *   SOLA_XPOS(I1,O1,K,N)
                MXR  = MCON(K,N)   *   SOLB_XNEG(I1,O1,K,N)
                LLXR = LCON(K,N)   * L_SOLA_XPOS(I1,O1,K,N,Q)
                MLXR = MCON(K,N)   * L_SOLB_XNEG(I1,O1,K,N,Q)
                NXR  = NCON(K,N,Q) *   SOLA_XPOS(I1,O1,K,N)
                PXR  = PCON(K,N,Q) *   SOLB_XNEG(I1,O1,K,N)
                HOM1 = ( NXR + LLXR ) *   T_UTDN_EIGEN(K,UT)
     &                       +  LXR   * L_T_UTDN_EIGEN(K,UT,Q)
                HOM2 = ( PXR + MLXR ) *   T_UTUP_EIGEN(K,UT)
     &                       +  MXR   * L_T_UTUP_EIGEN(K,UT,Q)
                SHOM_R = SHOM_R + HOM1 + HOM2 
              ENDDO
  
C  complex homogeneous solutions

              SHOM_CR = ZERO
              KO1 = K_REAL(N) + 1
              DO K = 1, K_COMPLEX(N)
                K0 = 2 * K - 2
                K1 = KO1 + K0
                K2 = K1  + 1

                NXR1  =   NCON(K1,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
     &                  - NCON(K2,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
                NXR2  =   NCON(K1,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
     &                  + NCON(K2,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
                PXR1  =   PCON(K1,N,Q) *   SOLB_XNEG(I1,O1,K1,N)
     &                  - PCON(K2,N,Q) *   SOLB_XNEG(I1,O1,K2,N)
                PXR2  =   PCON(K1,N,Q) *   SOLB_XNEG(I1,O1,K2,N)
     &                  + PCON(K2,N,Q) *   SOLB_XNEG(I1,O1,K1,N)

                LXR1  =   LCON(K1,N) *   SOLA_XPOS(I1,O1,K1,N)
     &                  - LCON(K2,N) *   SOLA_XPOS(I1,O1,K2,N)
                LXR2  =   LCON(K1,N) *   SOLA_XPOS(I1,O1,K2,N)
     &                  + LCON(K2,N) *   SOLA_XPOS(I1,O1,K1,N)
                MXR1  =   MCON(K1,N) *   SOLB_XNEG(I1,O1,K1,N)
     &                  - MCON(K2,N) *   SOLB_XNEG(I1,O1,K2,N)
                MXR2  =   MCON(K1,N) *   SOLB_XNEG(I1,O1,K2,N)
     &                  + MCON(K2,N) *   SOLB_XNEG(I1,O1,K1,N)

                LLXR1  =   LCON(K1,N) * L_SOLA_XPOS(I1,O1,K1,N,Q)
     &                   - LCON(K2,N) * L_SOLA_XPOS(I1,O1,K2,N,Q)
                LLXR2  =   LCON(K1,N) * L_SOLA_XPOS(I1,O1,K2,N,Q)
     &                   + LCON(K2,N) * L_SOLA_XPOS(I1,O1,K1,N,Q)
                MLXR1  =   MCON(K1,N) * L_SOLB_XNEG(I1,O1,K1,N,Q)
     &                   - MCON(K2,N) * L_SOLB_XNEG(I1,O1,K2,N,Q)
                MLXR2  =   MCON(K1,N) * L_SOLB_XNEG(I1,O1,K2,N,Q)
     &                   + MCON(K2,N) * L_SOLB_XNEG(I1,O1,K1,N,Q)

                HOM1CR =   ( NXR1 + LLXR1 ) *   T_UTDN_EIGEN(K1,UT)
     &                   - ( NXR2 + LLXR2 ) *   T_UTDN_EIGEN(K2,UT)
                HOM2CR =             LXR1   * L_T_UTDN_EIGEN(K1,UT,Q)
     &                             - LXR2   * L_T_UTDN_EIGEN(K2,UT,Q)
                HOM3CR =   ( PXR1 + MLXR1 ) *   T_UTUP_EIGEN(K1,UT)
     &                   - ( PXR2 + MLXR2 ) *   T_UTUP_EIGEN(K2,UT)
                HOM4CR =             MXR1   * L_T_UTUP_EIGEN(K1,UT,Q)
     &                             - MXR2   * L_T_UTUP_EIGEN(K2,UT,Q)

                SHOM_CR = SHOM_CR + HOM1CR + HOM2CR + HOM3CR + HOM4CR

              ENDDO

C  real part

              SHOM = SHOM_R + SHOM_CR
              QATMOSWF_F(Q,UTA,I,O1) = FMULT * SHOM

C  Finish Q, I and O1 loops

            ENDDO
          ENDDO
        ENDDO

C  Solution for N  below/above the layer NV that varies

      ELSE IF ( NV.NE.N .AND. NV.NE.0 ) THEN

C  stream/stokes directions

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES

C  Parameter loop

            DO Q = 1, NV_PARAMETERS

C  real homogeneous solutions

              SHOM_R = ZERO
              DO K = 1, K_REAL(N)
                NXR  = NCON(K,N,Q) *   SOLA_XPOS(I1,O1,K,N)
                PXR  = PCON(K,N,Q) *   SOLB_XNEG(I1,O1,K,N)
                HOM1 = NXR * T_UTDN_EIGEN(K,UT)
                HOM2 = PXR * T_UTUP_EIGEN(K,UT)
                SHOM_R = SHOM_R + HOM1 + HOM2 
              ENDDO
  
C  complex homogeneous solutions

              SHOM_CR = ZERO
              KO1 = K_REAL(N) + 1
              DO K = 1, K_COMPLEX(N)
                K0 = 2 * K - 2
                K1 = KO1 + K0
                K2 = K1  + 1
                NXR1  =   NCON(K1,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
     &                  - NCON(K2,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
                NXR2  =   NCON(K1,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
     &                  + NCON(K2,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
                PXR1  =   PCON(K1,N,Q) *   SOLB_XNEG(I1,O1,K1,N)
     &                  - PCON(K2,N,Q) *   SOLB_XNEG(I1,O1,K2,N)
                PXR2  =   PCON(K1,N,Q) *   SOLB_XNEG(I1,O1,K2,N)
     &                  + PCON(K2,N,Q) *   SOLB_XNEG(I1,O1,K1,N)
                HOM1CR =   NXR1  * T_UTDN_EIGEN(K1,UT)
     &                   - NXR2  * T_UTDN_EIGEN(K2,UT)
                HOM2CR =   PXR1  * T_UTUP_EIGEN(K1,UT)
     &                   - PXR2  * T_UTUP_EIGEN(K2,UT)
                SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
              ENDDO

C  real part

              SHOM = SHOM_R + SHOM_CR
              QATMOSWF_F(Q,UTA,I,O1) = FMULT * SHOM

C  Finish Q, I and O1 loops

            ENDDO
          ENDDO
        ENDDO

C  end variability clause

      ENDIF

C  Add the linearized thermal solution  (if flagged)
C    ---Only present if N = NV, or NV = 0 (column linearization)
C   THIS is the solution with scattering

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        O1 = 1
        IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO Q = 1, NV_PARAMETERS
              SPAR = L_UT_T_PARTIC(I1,UT,Q)
              QATMOSWF_F(Q,UTA,I,O1) = 
     &            QATMOSWF_F(Q,UTA,I,O1) + FMULT * SPAR
            ENDDO
          ENDDO
        ENDIF
      ENDIF

C  Finished if no solar terms

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

C  Add Linearized Classical beam solution
C  Green function solution is still a Placeholder....!
C    solutions exist only for N = NV or N > NV or NV = 0

      IF ( DO_CLASSICAL_SOLUTION ) THEN
        IF ( N.GE.NV .OR. NV.EQ.0 ) THEN
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES
              DO Q = 1, NV_PARAMETERS
                SPAR = L_WUPPER(I1,O1,N,Q) *   T_UTDN_MUBAR(UT,IB) +
     &                   WUPPER(I1,O1,N)   * L_T_UTDN_MUBAR(UT,NV,IB,Q)
                QATMOSWF_F(Q,UTA,I,O1) = 
     &            QATMOSWF_F(Q,UTA,I,O1) + FMULT * SPAR
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ELSE
C                P L A C E H O L D E R
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE QUADATMOSWF_OFFGRID_DN
     &     ( IB, UTA, UT, N, NV, NV_PARAMETERS,
     &       DO_INCLUDE_THERMEMISS, FLUX_MULTIPLIER, QATMOSWF_F )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and multiplier variables (input)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'
      INCLUDE '../includes/VLIDORT_THERMALSUP.VARS'

C  include files of linearized setup/solution/multiplier variables (input)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'
      INCLUDE '../includes/VLIDORT_L_THERMALSUP.VARS'

C  Subroutine arguments
C  --------------------

C  Inputs

      INTEGER          NV, NV_PARAMETERS
      INTEGER          IB, UTA, UT, N

      DOUBLE PRECISION FLUX_MULTIPLIER
      LOGICAL          DO_INCLUDE_THERMEMISS

C  output

      DOUBLE PRECISION QATMOSWF_F 
     &     ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

C  local variables
C  ---------------

      INTEGER          I, Q, O1, K, KO1, K0, K1, K2, LAY
      DOUBLE PRECISION SPAR, SHOM_R, SHOM_CR, SHOM
      DOUBLE PRECISION HOM1, HOM2, HOM1CR, HOM2CR, HOM3CR, HOM4CR
      DOUBLE PRECISION NXR, PXR, NXR1, NXR2, PXR1, PXR2
      DOUBLE PRECISION LXR, MXR, LXR1, MXR1, LXR2, MXR2
      DOUBLE PRECISION LLXR, MLXR, LLXR1, MLXR1, LLXR2, MLXR2
      DOUBLE PRECISION TPROP, THELP, L_TPROP, L_THELP, FMULT

C  Short hand

      FMULT = FLUX_MULTIPLIER

C  Thermal Transmittance only
c  --------------------------

      IF ( DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          DO Q = 1, NV_PARAMETERS
            L_THELP = ZERO
            THELP   = ZERO
            DO LAY = 1, N-1
              L_THELP = L_THELP *   T_DELT_DISORDS(I,LAY)
              IF ( LAY.EQ.NV .OR. NV.EQ.0 ) THEN
                L_TPROP = L_T_WLOWER(I,LAY,Q) / QUAD_STREAMS(I)
                L_THELP = L_THELP + L_TPROP
     &                    + THELP * L_T_DELT_DISORDS(I,LAY,Q) 
              ENDIF
              TPROP = T_WLOWER(I,LAY) / QUAD_STREAMS(I)
              THELP = THELP * T_DELT_DISORDS(I,LAY) + TPROP
            ENDDO
            L_THELP = L_THELP * T_DISORDS_UTDN(I,UT)
            IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN
              L_TPROP = L_UT_T_PARTIC(I,UT,Q) / QUAD_STREAMS(I)
              L_THELP = L_THELP + L_TPROP
     &                  + THELP * L_T_DISORDS_UTDN(I,UT,Q) 
            ENDIF
            QATMOSWF_F(Q,UTA,I,O1) = FMULT * L_THELP
          ENDDO
        ENDDO
        RETURN
      ENDIF

C  For those optical depths at off-grid levels
C  ###########################################
      
C  Homogeneous
C  -----------

C  solution for N being the varying layer NV

      IF ( N .EQ. NV .OR. NV.EQ.0 ) THEN

C  stream and Stokes directions

        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES

C  Parameter loop

            DO Q = 1, NV_PARAMETERS

C  real homogeneous solutions

              SHOM_R = ZERO
              DO K = 1, K_REAL(N)
                LXR  = LCON(K,N)   *   SOLA_XPOS(I,O1,K,N)
                MXR  = MCON(K,N)   *   SOLB_XNEG(I,O1,K,N)
                LLXR = LCON(K,N)   * L_SOLA_XPOS(I,O1,K,N,Q)
                MLXR = MCON(K,N)   * L_SOLB_XNEG(I,O1,K,N,Q)
                NXR  = NCON(K,N,Q) *   SOLA_XPOS(I,O1,K,N)
                PXR  = PCON(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
                HOM1 = ( NXR + LLXR ) *   T_UTDN_EIGEN(K,UT)
     &                       +  LXR   * L_T_UTDN_EIGEN(K,UT,Q)
                HOM2 = ( PXR + MLXR ) *   T_UTUP_EIGEN(K,UT)
     &                       +  MXR   * L_T_UTUP_EIGEN(K,UT,Q)
                SHOM_R = SHOM_R + HOM1 + HOM2
              ENDDO
  
C  complex homogeneous solutions

              SHOM_CR = ZERO
              KO1 = K_REAL(N) + 1
              DO K = 1, K_COMPLEX(N)
                K0 = 2 * K - 2
                K1 = KO1 + K0
                K2 = K1  + 1

                NXR1  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K1,N)
     &                  - NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K2,N)
                NXR2  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K2,N)
     &                  + NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K1,N)
                PXR1  =   PCON(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N)
     &                  - PCON(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)
                PXR2  =   PCON(K1,N,Q) *   SOLB_XNEG(I,O1,K2,N)
     &                  + PCON(K2,N,Q) *   SOLB_XNEG(I,O1,K1,N)

                LXR1  =   LCON(K1,N) *   SOLA_XPOS(I,O1,K1,N)
     &                  - LCON(K2,N) *   SOLA_XPOS(I,O1,K2,N)
                LXR2  =   LCON(K1,N) *   SOLA_XPOS(I,O1,K2,N)
     &                  + LCON(K2,N) *   SOLA_XPOS(I,O1,K1,N)
                MXR1  =   MCON(K1,N) *   SOLB_XNEG(I,O1,K1,N)
     &                  - MCON(K2,N) *   SOLB_XNEG(I,O1,K2,N)
                MXR2  =   MCON(K1,N) *   SOLB_XNEG(I,O1,K2,N)
     &                  + MCON(K2,N) *   SOLB_XNEG(I,O1,K1,N)

                LLXR1  =   LCON(K1,N) * L_SOLA_XPOS(I,O1,K1,N,Q)
     &                   - LCON(K2,N) * L_SOLA_XPOS(I,O1,K2,N,Q)
                LLXR2  =   LCON(K1,N) * L_SOLA_XPOS(I,O1,K2,N,Q)
     &                   + LCON(K2,N) * L_SOLA_XPOS(I,O1,K1,N,Q)
                MLXR1  =   MCON(K1,N) * L_SOLB_XNEG(I,O1,K1,N,Q)
     &                   - MCON(K2,N) * L_SOLB_XNEG(I,O1,K2,N,Q)
                MLXR2  =   MCON(K1,N) * L_SOLB_XNEG(I,O1,K2,N,Q)
     &                   + MCON(K2,N) * L_SOLB_XNEG(I,O1,K1,N,Q)

                HOM1CR =   ( NXR1 + LLXR1 ) *   T_UTDN_EIGEN(K1,UT)
     &                   - ( NXR2 + LLXR2 ) *   T_UTDN_EIGEN(K2,UT)
                HOM2CR =             LXR1   * L_T_UTDN_EIGEN(K1,UT,Q)
     &                             - LXR2   * L_T_UTDN_EIGEN(K2,UT,Q)
                HOM3CR =   ( PXR1 + MLXR1 ) *   T_UTUP_EIGEN(K1,UT)
     &                   - ( PXR2 + MLXR2 ) *   T_UTUP_EIGEN(K2,UT)
                HOM4CR =             MXR1   * L_T_UTUP_EIGEN(K1,UT,Q)
     &                             - MXR2   * L_T_UTUP_EIGEN(K2,UT,Q)

                SHOM_CR = SHOM_CR + HOM1CR + HOM2CR + HOM3CR + HOM4CR
              ENDDO

C  real part

              SHOM = SHOM_R + SHOM_CR
              QATMOSWF_F(Q,UTA,I,O1) = FMULT * SHOM

C  Finish Q, I and O1 loops

            ENDDO
          ENDDO
        ENDDO

C  Solution for N  below/above the layer NV that varies

      ELSE IF ( NV.NE.N .AND. NV.NE.0 ) THEN

C  stream/stokes directions

        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES

C  Parameter loop

            DO Q = 1, NV_PARAMETERS

C  real homogeneous solutions

              SHOM_R = ZERO
              DO K = 1, K_REAL(N)
                NXR  = NCON(K,N,Q) *   SOLA_XPOS(I,O1,K,N)
                PXR  = PCON(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
                HOM1 = NXR * T_UTDN_EIGEN(K,UT)
                HOM2 = PXR * T_UTUP_EIGEN(K,UT)
                SHOM_R = SHOM_R + HOM1 + HOM2 
              ENDDO
  
C  complex homogeneous solutions

              SHOM_CR = ZERO
              KO1 = K_REAL(N) + 1
              DO K = 1, K_COMPLEX(N)
                K0 = 2 * K - 2
                K1 = KO1 + K0
                K2 = K1  + 1
                NXR1  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K1,N)
     &                  - NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K2,N)
                NXR2  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K2,N)
     &                  + NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K1,N)
                PXR1  =   PCON(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N)
     &                  - PCON(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)
                PXR2  =   PCON(K1,N,Q) *   SOLB_XNEG(I,O1,K2,N)
     &                  + PCON(K2,N,Q) *   SOLB_XNEG(I,O1,K1,N)
                HOM1CR =   NXR1  * T_UTDN_EIGEN(K1,UT)
     &                   - NXR2  * T_UTDN_EIGEN(K2,UT)
                HOM2CR =   PXR1  * T_UTUP_EIGEN(K1,UT)
     &                   - PXR2  * T_UTUP_EIGEN(K2,UT)
                SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
              ENDDO

C  real part

              SHOM = SHOM_R + SHOM_CR
              QATMOSWF_F(Q,UTA,I,O1) = FMULT * SHOM

C  Finish Q, I and O1 loops

            ENDDO
          ENDDO
        ENDDO

C  end variability clause

      ENDIF

C  Add the linearized thermal solution  (if flagged)
C    ---Only present if N = NV, or NV = 0
C   THIS IS THE SOLUTION with scattering

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        O1 = 1
        IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN
          DO I = 1, NSTREAMS
            DO Q = 1, NV_PARAMETERS
              SPAR = L_UT_T_PARTIC(I,UT,Q)
              QATMOSWF_F(Q,UTA,I,O1) =  
     &           QATMOSWF_F(Q,UTA,I,O1) + FMULT * SPAR
            ENDDO
          ENDDO
        ENDIF
      ENDIF

C  Finished if no solar terms

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

C  Add Linearized Classical beam solution
C  Green function solution is still a Placeholder....!
C    solutions exist only for N = NV or N > NV or NV = 0

      IF ( DO_CLASSICAL_SOLUTION ) THEN
        IF ( N.GE.NV .OR. NV.EQ.0 ) THEN
          DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES
              DO Q = 1, NV_PARAMETERS
                SPAR = L_WUPPER(I,O1,N,Q) *   T_UTDN_MUBAR(UT,IB) +
     &                   WUPPER(I,O1,N)   * L_T_UTDN_MUBAR(UT,NV,IB,Q)
                QATMOSWF_F(Q,UTA,I,O1) = 
     &            QATMOSWF_F(Q,UTA,I,O1) + FMULT * SPAR
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ELSE
C                P L A C E H O L D E R
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE VLIDORT_LS_INTEGRATED_OUTPUT
     I  ( DO_INCLUDE_MVOUTPUT, 
     I    FLUX_MULTIPLIER, IBEAM, WF_INDEX,
     I    LS_BOA_THTONLY_SOURCE )

C  Quadrature output at offgrid or ongrid optical depths
C  ( Required if mean-value calculations are to be done)
C    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include file of setup variables (input)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  include file of linearized setupr variables (input)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'

C  include file of result variables (module output stored here)

      INCLUDE '../includes/VLIDORT_L_RESULTS.VARS'

C  Subroutine arguments
C  --------------------

      LOGICAL          DO_INCLUDE_MVOUTPUT
      INTEGER          IBEAM, WF_INDEX
      DOUBLE PRECISION FLUX_MULTIPLIER

      DOUBLE PRECISION LS_BOA_THTONLY_SOURCE(MAXSTREAMS)

C  local variables
C  ---------------

C  Local qudrature output (for debug)

      DOUBLE PRECISION QSURFACEWF_F 
     &     ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

C  Help variables

      INTEGER          I, IDIR, WDIR, UTA, UT, N, NLEVEL, O1, QW
      DOUBLE PRECISION SMI, SFX

C  Index

      QW = WF_INDEX

C  direction loop

      DO IDIR = 1, N_DIRECTIONS
        WDIR = WHICH_DIRECTIONS(IDIR)

C  Upwelling Jacobian output at Quadrature angles

        IF ( WDIR .EQ. UPIDX ) THEN
          DO UTA = 1, N_USER_LEVELS
            NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
              UT = PARTLAYERS_OUTINDEX(UTA)
              N  = PARTLAYERS_LAYERIDX(UT)
              CALL QUADSURFACEWF_OFFGRID_UP
     &          ( UTA, UT, N, QW, FLUX_MULTIPLIER, 
     &            LS_BOA_THTONLY_SOURCE, QSURFACEWF_F )
            ELSE
              CALL QUADSURFACEWF_LEVEL_UP
     &          ( UTA, NLEVEL, QW, FLUX_MULTIPLIER, 
     &            LS_BOA_THTONLY_SOURCE, QSURFACEWF_F )
            ENDIF
          ENDDO
        ENDIF

C  Downwelling Jacobian output at Quadrature angles

        IF ( WDIR .EQ. DNIDX ) THEN
          DO UTA = 1, N_USER_LEVELS
            NLEVEL = UTAU_LEVEL_MASK_DN(UTA)
            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
              UT = PARTLAYERS_OUTINDEX(UTA)
              N  = PARTLAYERS_LAYERIDX(UT)
              CALL QUADSURFACEWF_OFFGRID_DN
     &          ( UTA, UT, N, QW, FLUX_MULTIPLIER, 
     &            QSURFACEWF_F )
            ELSE
              CALL QUADSURFACEWF_LEVEL_DN
     &          ( UTA, NLEVEL, QW, FLUX_MULTIPLIER, 
     &            QSURFACEWF_F )
            ENDIF
          ENDDO
        ENDIF

C  Mean Intensity and Flux  output (diffuse term only)

        IF ( DO_INCLUDE_MVOUTPUT ) THEN
          DO UTA = 1, N_USER_LEVELS
            DO O1 = 1, NSTOKES
              SMI = ZERO
              SFX = ZERO
              DO I = 1, NSTREAMS
                SMI = SMI + QUAD_WEIGHTS(I) * QSURFACEWF_F(QW,UTA,I,O1)
                SFX = SFX + QUAD_STRMWTS(I) * QSURFACEWF_F(QW,UTA,I,O1)
              ENDDO
              MINT_SURFACEWF(QW,UTA,IBEAM,O1,WDIR) = SMI * HALF
              FLUX_SURFACEWF(QW,UTA,IBEAM,O1,WDIR) = SFX * PI2
            ENDDO
          ENDDO
        ENDIF

C  End directions loop

      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE QUADSURFACEWF_LEVEL_UP
     &      ( UTA, NL, QW, FLUX_MULTIPLIER,
     &        LS_BOA_THTONLY_SOURCE, QSURFACEWF_F )

C  Upwelling weighting function Fourier components at level boundary NL
C  Quadrature angles only

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup, solution variables (input)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'

C  include files of linearized setup and solution variables (input)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'

C  Subroutine arguments
C  --------------------

C  Inputs

      INTEGER          UTA, NL, QW
      DOUBLE PRECISION FLUX_MULTIPLIER
      DOUBLE PRECISION LS_BOA_THTONLY_SOURCE(MAXSTREAMS)

C  output

      DOUBLE PRECISION QSURFACEWF_F 
     &     ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

C  local variables
C  ---------------

      INTEGER          N, I, I1, O1, K, KO1, K0, K1, K2, LAY
      DOUBLE PRECISION SHOM_R, SHOM_CR, H1, H2, THELP
      DOUBLE PRECISION NXR, NXR1, NXR2, PXR, PXR1, PXR2

C  This depends on the level mask - if this is 0 to NLAYERS - 1, then we are
C  looking at the perturbation field at the top of these layers. The
C  case where the level mask = NLAYERS is the upwelling perturbed fields
C  at the bottom of the atmosphere (treated separately).

      N = NL + 1

C  For the lowest level
!  ====================

C  Thermal transmittance-only solution

      IF ( NL.EQ.NLAYERS .and. DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          THELP = LS_BOA_THTONLY_SOURCE(I)
          QSURFACEWF_F(QW,UTA,I,O1) = FLUX_MULTIPLIER * THELP
        ENDDO
        RETURN
      ENDIF

C  Scattering solutions
C  --------------------

      IF ( NL .EQ. NLAYERS ) THEN

c  Offset

        KO1 = K_REAL(NL) + 1

C  Stokes and streams loops

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES

C  Real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(NL)
              NXR = NCON_ALB(K,NL) * SOLA_XPOS(I1,O1,K,NL)
              PXR = PCON_ALB(K,NL) * SOLB_XNEG(I1,O1,K,NL)
              SHOM_R = SHOM_R + NXR*T_DELT_EIGEN(K,NL) + PXR
            ENDDO

C  Complex nomoegeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(NL)
             K0 = 2 * K - 2
             K1 = KO1 + K0
             K2 = K1  + 1
             NXR1 =   NCON_ALB(K1,NL) * SOLA_XPOS(I1,O1,K1,NL)
     &              - NCON_ALB(K2,NL) * SOLA_XPOS(I1,O1,K2,NL)
             NXR2 =   NCON_ALB(K1,NL) * SOLA_XPOS(I1,O1,K2,NL)
     &              + NCON_ALB(K2,NL) * SOLA_XPOS(I1,O1,K1,NL)
             PXR1 =   PCON_ALB(K1,NL) * SOLB_XNEG(I1,O1,K1,NL)
     &              - PCON_ALB(K2,NL) * SOLB_XNEG(I1,O1,K2,NL)
             H1 =   NXR1 *  * T_DELT_EIGEN(K1,NL)
     &            - NXR2 *  * T_DELT_EIGEN(K2,NL)
             H2 =  PXR1
             SHOM_CR = SHOM_CR + H1 + H2
            ENDDO

C  collect solution

            QSURFACEWF_F(QW,UTA,I,O1) = FLUX_MULTIPLIER*(SHOM_R+SHOM_CR)

C  Finish loops

          ENDDO
        ENDDO

      ENDIF

C  For other levels in the atmosphere
C  ==================================

C  Thermal transmittance-only solution

      IF ( NL.NE.NLAYERS .and. DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          THELP = LS_BOA_THTONLY_SOURCE(I)
          DO LAY = NLAYERS, N, -1
            THELP = THELP*T_DELT_DISORDS(I,LAY)
          ENDDO
          QSURFACEWF_F(QW,UTA,I,O1) = FLUX_MULTIPLIER * THELP
        ENDDO
        RETURN
      ENDIF

C  Scattering solutions
C  --------------------

      IF ( NL .NE. NLAYERS ) THEN

c  Offset

        KO1 = K_REAL(N) + 1

C  Stokes and streams loops

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES

C  Real homogeneous solutions

             SHOM_R = ZERO
            DO K = 1, K_REAL(N)
              NXR = NCON_ALB(K,N) * SOLA_XPOS(I1,O1,K,N)
              PXR = PCON_ALB(K,N) * SOLB_XNEG(I1,O1,K,N)
              SHOM_R = SHOM_R + PXR*T_DELT_EIGEN(K,N) + NXR
            ENDDO

C  Complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(N)
             K0 = 2 * K - 2
             K1 = KO1 + K0
             K2 = K1  + 1
             NXR1 =   NCON_ALB(K1,N) * SOLA_XPOS(I1,O1,K1,N)
     &              - NCON_ALB(K2,N) * SOLA_XPOS(I1,O1,K2,N)
             PXR1 =   PCON_ALB(K1,N) * SOLB_XNEG(I1,O1,K1,N)
     &              - PCON_ALB(K2,N) * SOLB_XNEG(I1,O1,K2,N)
             PXR2 =   PCON_ALB(K1,N) * SOLB_XNEG(I1,O1,K2,N)
     &              + PCON_ALB(K2,N) * SOLB_XNEG(I1,O1,K1,N)
             H1 =   PXR1 *  * T_DELT_EIGEN(K1,N)
     &            - PXR2 *  * T_DELT_EIGEN(K2,N)
             H2 =  NXR1
             SHOM_CR = SHOM_CR + H1 + H2
            ENDDO

C  collect solution

           QSURFACEWF_F(QW,UTA,I,O1) = FLUX_MULTIPLIER*(SHOM_R+SHOM_CR)

C  Finish loops

          ENDDO
        ENDDO

      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE QUADSURFACEWF_LEVEL_DN
     &      ( UTA, NL, QW, FLUX_MULTIPLIER,
     &        QSURFACEWF_F )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup, solution variables (input)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'

C  include files of linearized setup and solution variables (input)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'

C  Subroutine arguments
C  --------------------

C  Inputs

      INTEGER          UTA, NL, QW
      DOUBLE PRECISION FLUX_MULTIPLIER

C  output

      DOUBLE PRECISION QSURFACEWF_F 
     &     ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

C  local variables
C  ---------------

      INTEGER          N, I, O1, K, KO1, K0, K1, K2
      DOUBLE PRECISION SHOM_R, SHOM_CR, H1, H2
      DOUBLE PRECISION NXR, NXR1, NXR2, PXR, PXR1

      N = NL

C  Zero level = TOA no solutions
C  =============================

C  Downwelling weighting function at TOA ( or N = 0 ) is zero

      IF ( NL .EQ. 0 ) THEN
        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES
            QSURFACEWF_F(QW,UTA,I,O1) = ZERO
          ENDDO
        ENDDO
        RETURN
      ENDIF

C  For other levels in the atmosphere
C  ==================================

      IF ( NL .NE. 0 .and. DO_THERMAL_TRANSONLY) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          QSURFACEWF_F(QW,UTA,I,O1) = ZERO
        ENDDO
        RETURN
      endif

C  Scattering solutions
C  --------------------

      IF ( NL .NE. 0 ) THEN

c  Offset

        KO1 = K_REAL(N) + 1

C  Stokes and streams loops

        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES

C  Real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(N)
              NXR = NCON_ALB(K,N) * SOLA_XPOS(I,O1,K,N)
              PXR = PCON_ALB(K,N) * SOLB_XNEG(I,O1,K,N)
              SHOM_R = SHOM_R + NXR*T_DELT_EIGEN(K,N) + PXR
            ENDDO

C  Complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(N)
             K0 = 2 * K - 2
             K1 = KO1 + K0
             K2 = K1  + 1
             NXR1 =   NCON_ALB(K1,N) * SOLA_XPOS(I,O1,K1,N)
     &              - NCON_ALB(K2,N) * SOLA_XPOS(I,O1,K2,N)
             NXR2 =   NCON_ALB(K1,N) * SOLA_XPOS(I,O1,K2,N)
     &              + NCON_ALB(K2,N) * SOLA_XPOS(I,O1,K1,N)
             PXR1 =   PCON_ALB(K1,N) * SOLB_XNEG(I,O1,K1,N)
     &              - PCON_ALB(K2,N) * SOLB_XNEG(I,O1,K2,N)
             H1 =   NXR1 *  * T_DELT_EIGEN(K1,N)
     &            - NXR2 *  * T_DELT_EIGEN(K2,N)
             H2 =  PXR1
             SHOM_CR = SHOM_CR + H1 + H2
            ENDDO

C  collect solution

            QSURFACEWF_F(QW,UTA,I,O1) = FLUX_MULTIPLIER*(SHOM_R+SHOM_CR)

C  Finish loops

          ENDDO
        ENDDO
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE QUADSURFACEWF_OFFGRID_UP
     &      ( UTA, UT, N, QW, FLUX_MULTIPLIER,
     &        LS_BOA_THTONLY_SOURCE, QSURFACEWF_F )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and multiplier variables (input)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  include files of linearized setup/solution/multiplier variables (input)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'

C  Subroutine arguments
C  --------------------

C  Inputs

      INTEGER          UTA, UT, N, QW
      DOUBLE PRECISION FLUX_MULTIPLIER
      DOUBLE PRECISION LS_BOA_THTONLY_SOURCE(MAXSTREAMS)

C  output

      DOUBLE PRECISION QSURFACEWF_F 
     &     ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

C  local variables
C  ---------------

      INTEGER          I, I1, O1, K, KO1, K0, K1, K2, LAY
      DOUBLE PRECISION SHOM_R, SHOM_CR, H1, H2, THELP
      DOUBLE PRECISION NXR, NXR1, NXR2, PXR, PXR1, PXR2

C  For thermal transmittance-only
C  ------------------------------

      IF ( DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          THELP = LS_BOA_THTONLY_SOURCE(I)
          DO LAY = NLAYERS, N+1, -1
            THELP = THELP*T_DELT_DISORDS(I,LAY)
          ENDDO
          THELP = THELP*T_DISORDS_UTUP(I,UT)
          QSURFACEWF_F(QW,UTA,I,O1) = FLUX_MULTIPLIER * THELP
        ENDDO
        RETURN
      ENDIF

C  For those optical depths at off-grid levels
C  -------------------------------------------

C  Offset

      KO1 = K_REAL(N) + 1

C  stream directions

      DO I = 1, NSTREAMS
        I1 = I + NSTREAMS
        DO O1 = 1, NSTOKES

C  real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, K_REAL(N)
            NXR = NCON_ALB(K,N) * SOLA_XPOS(I1,O1,K,N)
            PXR = PCON_ALB(K,N) * SOLB_XNEG(I1,O1,K,N)
            SHOM_R = SHOM_R + NXR * T_UTDN_EIGEN(K,UT)
     &                      + PXR * T_UTUP_EIGEN(K,UT)
            SHOM_R = SHOM_R + H1 + H2
          ENDDO
  
C  complex homogeneous solutions

          SHOM_CR = ZERO
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            NXR1 =   NCON_ALB(K1,N) * SOLA_XPOS(I1,O1,K1,N)
     &             - NCON_ALB(K2,N) * SOLA_XPOS(I1,O1,K2,N)
            NXR2 =   NCON_ALB(K1,N) * SOLA_XPOS(I1,O1,K2,N)
     &             + NCON_ALB(K2,N) * SOLA_XPOS(I1,O1,K1,N)
            PXR1 =   PCON_ALB(K1,N) * SOLB_XNEG(I1,O1,K1,N)
     &             - PCON_ALB(K2,N) * SOLB_XNEG(I1,O1,K2,N)
            PXR2 =   PCON_ALB(K1,N) * SOLB_XNEG(I1,O1,K2,N)
     &             + PCON_ALB(K2,N) * SOLB_XNEG(I1,O1,K1,N)
            H1 =   NXR1 *  * T_UTDN_EIGEN(K1,N)
     &           - NXR2 *  * T_UTDN_EIGEN(K2,N)
            H2 =   PXR1 *  * T_UTUP_EIGEN(K1,N)
     &           - PXR2 *  * T_UTUP_EIGEN(K2,N)
            SHOM_CR = SHOM_CR + H1 + H2
          ENDDO

C  set result

          QSURFACEWF_F(QW,UTA,I,O1) = FLUX_MULTIPLIER*(SHOM_R+SHOM_CR)

C  end O1 and I loops

        ENDDO
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE QUADSURFACEWF_OFFGRID_DN
     &      ( UTA, UT, N, QW, FLUX_MULTIPLIER,
     &        QSURFACEWF_F )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and multiplier variables (input)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  include files of linearized setup/solution/multiplier variables (input)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'

C  Subroutine arguments
C  --------------------

C  Inputs

      INTEGER          UTA, UT, N, QW
      DOUBLE PRECISION FLUX_MULTIPLIER

C  output

      DOUBLE PRECISION QSURFACEWF_F 
     &     ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

C  local variables
C  ---------------

      INTEGER          I, O1, K, KO1, K0, K1, K2
      DOUBLE PRECISION SHOM_R, SHOM_CR, H1, H2
      DOUBLE PRECISION NXR, NXR1, NXR2, PXR, PXR1, PXR2

C  For thermal transmittance-only, no contribution
C  -----------------------------------------------

      IF ( DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          QSURFACEWF_F(QW,UTA,I,O1) = ZERO
        ENDDO
        RETURN
      ENDIF
C  For those optical depths at off-grid levels
C  -------------------------------------------

C  Offset

      KO1 = K_REAL(N) + 1

C  stream directions

      DO I = 1, NSTREAMS
        DO O1 = 1, NSTOKES

C  real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, K_REAL(N)
            NXR = NCON_ALB(K,N) * SOLA_XPOS(I,O1,K,N)
            PXR = PCON_ALB(K,N) * SOLB_XNEG(I,O1,K,N)
            SHOM_R = SHOM_R + NXR * T_UTDN_EIGEN(K,UT)
     &                      + PXR * T_UTUP_EIGEN(K,UT)
            SHOM_R = SHOM_R + H1 + H2
          ENDDO
  
C  complex homogeneous solutions

          SHOM_CR = ZERO
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            NXR1 =   NCON_ALB(K1,N) * SOLA_XPOS(I,O1,K1,N)
     &             - NCON_ALB(K2,N) * SOLA_XPOS(I,O1,K2,N)
            NXR2 =   NCON_ALB(K1,N) * SOLA_XPOS(I,O1,K2,N)
     &             + NCON_ALB(K2,N) * SOLA_XPOS(I,O1,K1,N)
            PXR1 =   PCON_ALB(K1,N) * SOLB_XNEG(I,O1,K1,N)
     &             - PCON_ALB(K2,N) * SOLB_XNEG(I,O1,K2,N)
            PXR2 =   PCON_ALB(K1,N) * SOLB_XNEG(I,O1,K2,N)
     &             + PCON_ALB(K2,N) * SOLB_XNEG(I,O1,K1,N)
            H1 =   NXR1 *  * T_UTDN_EIGEN(K1,N)
     &           - NXR2 *  * T_UTDN_EIGEN(K2,N)
            H2 =   PXR1 *  * T_UTUP_EIGEN(K1,N)
     &           - PXR2 *  * T_UTUP_EIGEN(K2,N)
            SHOM_CR = SHOM_CR + H1 + H2
          ENDDO

C  set result

          QSURFACEWF_F(QW,UTA,I,O1) = FLUX_MULTIPLIER*(SHOM_R+SHOM_CR)

C  end O1 and I loops

        ENDDO
      ENDDO

C  Finish

      RETURN
      END
