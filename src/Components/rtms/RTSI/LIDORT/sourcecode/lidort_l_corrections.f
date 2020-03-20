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
C # Nadir single scatter corrections:                           #
C #                                                             #
C #          LIDORT_LP_SSCORR_NADIR(master) (renamed, V3.3)     #
C #          LIDORT_LC_SSCORR_NADIR(master) (new, Version 3.3)  #
C #                                                             #
C # Direct Beam corrections:                                    #
C #                                                             #
C #          LIDORT_LAP_DBCORRECTION(master) (renamed)          #
C #          LIDORT_LAC_DBCORRECTION(master) (new, Version 3.3) #
C #          LIDORT_LS_DBCORRECTION(master)                     #
C #                                                             #
C # Versions 3.2, 3.3: Outgoing sphericity correction:          #
C #               Profile weighting functions (3.2)             #
C #               Column weighting functions  (3.3)             #
C #                                                             #
C #            LIDORT_L_SSCORR_OUTGOING (master)                #
C #                 L_outgoing_integration_up                   #
C #                 L_outgoing_integration_dn                   #
C #                                                             #
C #      Version 3.3 outgoing Lambertian DB corrections         #
C #                                                             #
C #            LIDORT_LAP_LAMBERTIAN_DBCORRECTION               #
C #            LIDORT_LAC_LAMBERTIAN_DBCORRECTION               #
C #            LIDORT_LS_LAMBERTIAN_DBCORRECTION                #
C #                                                             #
C ###############################################################

      SUBROUTINE LIDORT_LP_SSCORR_NADIR
     &      ( SSFLUX, LAYER_TO_VARY, K_PARAMETERS)

C  Single scatter exact calculation
C   Programmed by R. Spurr, RT Solutions Inc.
C    Second Draft, April 14th 2005.
C    Third Draft,  May    6th 2005.
C       - additional code to deal with refraction.
C    Fourth Draft, May 2007, Names changed for output arrays.

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file of linearized input variables

      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'

C  include files of setup and multiplier variables (input)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  Include file of single scatter result variables

      INCLUDE '../includes/LIDORT_SINGSCAT.VARS'

C  include file with offset geometrical variables

      INCLUDE '../includes/LIDORT_RESULTS.VARS'

C  include files of linearized setup, solution and multiplier variables (input)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_L_MULTIPLIERS.VARS'

C  Include file of results (linearized single scatter weighting funcs)

      INCLUDE '../includes/LIDORT_L_SINGSCAT.VARS'

C  Input arguments
C  ---------------

C  Flux factor

      DOUBLE PRECISION SSFLUX

C  Linearization control

      INTEGER          LAYER_TO_VARY, K_PARAMETERS

C  local variables
C  ---------------

C  indices

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL, V
      INTEGER          UT, L, UTA, UM, IB, IA, NC, Q, NM1, K

C  other help variables

      DOUBLE PRECISION HELP, TTT, UVAR, AVAR, FT1, FT2, DNM1
      DOUBLE PRECISION L_FINAL_SOURCE, L_SS_LAYERSOURCE
      DOUBLE PRECISION L_SSCORRECTION, L_TR_CUMSOURCE
      DOUBLE PRECISION VAR_TMS(MAX_ATMOSWFS), LEGPOLY, L_PHASMOM
      DOUBLE PRECISION HELP1, HELP2, FDNL1, DNL1, FACT, GK11

C  Set up operations
C  -----------------

      K = LAYER_TO_VARY

C  Create Linearized TMS factor for Layer K
C   ( Use UNSCALED linearized inputs - the Nakajima-Tanaka way)

      IF ( DO_SSCORR_TRUNCATION ) THEN
        TTT = TMS(K) /  ( ONE - SSFDEL(K) )
      ELSE
        TTT = TMS(K)
      ENDIF
 
      IF ( DO_DELTAM_SCALING ) THEN
        NM1 = NMOMENTS+1
        FT1 = TRUNC_FACTOR(K) * TTT
        FT2 = ONE + FT1
        DO Q = 1, K_PARAMETERS
          IF ( DO_PHASFUNC_VARIATION(Q,K) ) THEN
            UVAR = L_OMEGA_TOTAL_INPUT(Q,K)
            AVAR = UVAR + L_PHASMOMS_TOTAL_INPUT(Q,NM1,K)
            VAR_TMS(Q) = UVAR + FT1 * AVAR
          ELSE
            UVAR = L_OMEGA_TOTAL_INPUT(Q,K)
            VAR_TMS(Q) = UVAR * FT2
          ENDIF
        ENDDO
      ELSE
        DO Q = 1, K_PARAMETERS
          UVAR = L_OMEGA_TOTAL_INPUT(Q,K)
          VAR_TMS(Q) = UVAR
        ENDDO
      ENDIF

C  Additional Delta-M scaling
C  --------------------------

C  New section. R. Spurr, 07 September 2007.

      IF ( DO_SSCORR_TRUNCATION ) THEN
        NM1  = NMOMENTS_INPUT
        DNM1 = DBLE(2*NM1+1)
        DO Q = 1, K_PARAMETERS
          IF ( DO_PHASFUNC_VARIATION(Q,K) ) THEN
            FT1 = ONE / ( ONE- SSFDEL(K) )
            L_SSFDEL(K,Q) = PHASMOMS_TOTAL_INPUT(NM1,K) *
     *                  L_PHASMOMS_TOTAL_INPUT(Q,NM1,K) / DNM1
            VAR_TMS(Q) = VAR_TMS(Q) - FT1 * L_SSFDEL(K,Q)
          ENDIF
        ENDDO
      ENDIF

C  ####################
C  #    UPWELLING     #
C  ####################

      IF ( DO_UPWELLING ) THEN

C  ===================================
C  Total phase function linearization (upwelling)
C  ===================================

        IF ( STERM_LAYERMASK_UP(K)) THEN

C  Loop over all geometries V 
C  Loop over varying parameters Q for layer K

          DO V = 1, N_GEOMETRIES
            DO Q = 1, K_PARAMETERS

C  Phase function moment variations
C    add TMS correction factor linearization

              IF ( DO_PHASFUNC_VARIATION(Q,K) ) THEN
                HELP = ZERO
                DO L = 0, NMOMENTS_INPUT
                  IF ( DO_SSCORR_TRUNCATION ) THEN
                    DNL1  = DBLE(2*L + 1 )
                    FDNL1 = SSFDEL(K) * DNL1
                    FACT  = ONE - SSFDEL(K)
                    GK11 = ( PHASMOMS_TOTAL_INPUT(L,K) - FDNL1 ) / FACT
                    HELP1 = L_PHASMOMS_TOTAL_INPUT(Q,L,K) * 
     &                        PHASMOMS_TOTAL_INPUT(L,K)
                    HELP2 = ( GK11 - DNL1 ) * L_SSFDEL(K,Q)
                    L_PHASMOM = ( HELP1 + HELP2 ) / FACT
                  ELSE
                    L_PHASMOM = L_PHASMOMS_TOTAL_INPUT(Q,L,K) *
     &                            PHASMOMS_TOTAL_INPUT(L,K)
                  ENDIF 
                  LEGPOLY = SS_PLEG_UP(V,K,L)
                  HELP = HELP + L_PHASMOM*LEGPOLY
                ENDDO
                L_EXACTSCAT_UP(V,K,Q) = HELP * TMS(K) +
     &                       EXACTSCAT_UP(V,K) * VAR_TMS(Q)
              ELSE
                L_EXACTSCAT_UP(V,K,Q) = EXACTSCAT_UP(V,K) * VAR_TMS(Q)
              ENDIF

C  end parameter loop and viewing directions loop

            ENDDO
          ENDDO

C  Only if layer K exists

        ENDIF

C  ===================================
C  Upwelling single scatter recurrence
C  ===================================

C  initialize cumulative source term

        NC = 0
        DO V = 1, N_GEOMETRIES
          DO Q = 1, K_PARAMETERS
            L_SS_CUMSOURCE(Q,V) = ZERO
          ENDDO
        ENDDO

C  initialise optical depth loop

        NSTART = NLAYERS
        NUT_PREV = NSTART + 1

C  Main loop over all output optical depths
C  ========================================

        DO UTA = N_OUT_USERTAUS, 1, -1

C  Layer index for given optical depth

          NLEVEL = UTAU_LEVEL_MASK_UP(UTA)

C  finishing layer

          NUT = NLEVEL + 1

C  Cumulative single scatter source terms to layer NUT
C  ---------------------------------------------------

C    1. Get layer source terms = Exact scattering * Multiplier
C    2. Loop over layers working upwards to NUT

          DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N

C  If N = K (the layer that is varying)

            IF ( N .EQ. K ) THEN

              DO IB = 1, NBEAMS
                DO UM = 1, N_USER_STREAMS
                  DO IA = 1, N_USER_RELAZMS
                    V = UMOFF(IB,UM) + IA
                    DO Q = 1, K_PARAMETERS
                      L_SS_LAYERSOURCE = 
     &                    EXACTSCAT_UP(V,N)   * L_EMULT_UP(UM,N,K,IB,Q)
     &                + L_EXACTSCAT_UP(V,N,Q) *   EMULT_UP(UM,N,IB)
                      L_SS_LAYERSOURCE = FLUX_FACTOR * L_SS_LAYERSOURCE
                      L_SS_CUMSOURCE(Q,V) = L_SS_LAYERSOURCE
     &            +   T_DELT_USERM(N,UM)   * L_SS_CUMSOURCE(Q,V)
     &            + L_T_DELT_USERM(N,UM,Q) *   SS_CUMSOURCE_UP(V,NC-1)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO

C  Variations when N > K

            ELSE IF ( N .GT. K ) THEN

              DO IB = 1, NBEAMS
                DO UM = 1, N_USER_STREAMS
                  DO IA = 1, N_USER_RELAZMS
                    V = UMOFF(IB,UM) + IA
                    DO Q = 1, K_PARAMETERS
                      L_SS_LAYERSOURCE = 
     &                  EXACTSCAT_UP(V,N) * L_EMULT_UP(UM,N,K,IB,Q)
                      L_SS_LAYERSOURCE = FLUX_FACTOR * L_SS_LAYERSOURCE
                      L_SS_CUMSOURCE(Q,V) = L_SS_LAYERSOURCE
     &                  + T_DELT_USERM(N,UM) * L_SS_CUMSOURCE(Q,V)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO

C  Transmittance for N < K

            ELSE

              DO IB = 1, NBEAMS
                DO UM = 1, N_USER_STREAMS
                  DO IA = 1, N_USER_RELAZMS
                    V = UMOFF(IB,UM) + IA
                    DO Q = 1, K_PARAMETERS
                      L_SS_CUMSOURCE(Q,V) =
     &                  T_DELT_USERM(N,UM) * L_SS_CUMSOURCE(Q,V)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO

            ENDIF

          ENDDO

C  Offgrid output
C  --------------

C  Set final cumulative source and Single scatter Weighting function

          IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN

            UT = OFFGRID_UTAU_OUTINDEX(UTA)
            N  = OFFGRID_UTAU_LAYERIDX(UT)

C  If N = K (the layer that is varying)
C    add Linearization of additional partial layer source term =
C        L_Exact_Scat(n) * Multiplier  +  Exact_Scat * L_Multiplier(n)

            IF ( N .EQ. K ) THEN

              DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = UMOFF(IB,UM) + IA
                 DO Q = 1, K_PARAMETERS
                  L_SS_LAYERSOURCE = 
     &               EXACTSCAT_UP(V,N)   * L_UT_EMULT_UP(UM,UT,N,IB,Q)
     &           + L_EXACTSCAT_UP(V,N,Q) *   UT_EMULT_UP(UM,UT,IB)
                  L_SS_LAYERSOURCE = FLUX_FACTOR * L_SS_LAYERSOURCE
                  L_TR_CUMSOURCE = 
     &           + L_T_UTUP_USERM(UT,UM,Q) *   SS_CUMSOURCE_UP(V,NC)
     &           +   T_UTUP_USERM(UT,UM)   * L_SS_CUMSOURCE(Q,V)
                  L_FINAL_SOURCE   = L_TR_CUMSOURCE + L_SS_LAYERSOURCE
                  L_SSCORRECTION   = SSFLUX * L_FINAL_SOURCE
                  PROFILEWF_SS(Q,K,UTA,V,UPIDX) = L_SSCORRECTION
                 ENDDO
                ENDDO
               ENDDO
              ENDDO

C  Variations when N > K
C    add Linearization of additional partial layer source term =
C         Exact_Scat * L_Multiplier(k)

            ELSE IF ( N .GT. K ) THEN

              DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = UMOFF(IB,UM) + IA
                 DO Q = 1, K_PARAMETERS
                  L_SS_LAYERSOURCE = 
     &               EXACTSCAT_UP(V,N) * L_UT_EMULT_UP(UM,UT,K,IB,Q)
                  L_SS_LAYERSOURCE = FLUX_FACTOR * L_SS_LAYERSOURCE
                  L_TR_CUMSOURCE = 
     &               T_UTUP_USERM(UT,UM)   * L_SS_CUMSOURCE(Q,V)
                  L_FINAL_SOURCE   = L_TR_CUMSOURCE + L_SS_LAYERSOURCE
                  L_SSCORRECTION   = SSFLUX * L_FINAL_SOURCE
                  PROFILEWF_SS(Q,K,UTA,V,UPIDX) = L_SSCORRECTION
                 ENDDO
                ENDDO
               ENDDO
              ENDDO

C  Transmittance for the other layers

            ELSE

              DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = UMOFF(IB,UM) + IA
                 DO Q = 1, K_PARAMETERS
                  L_FINAL_SOURCE = 
     &               T_UTUP_USERM(UT,UM)   * L_SS_CUMSOURCE(Q,V)
                  L_SSCORRECTION   = SSFLUX * L_FINAL_SOURCE
                  PROFILEWF_SS(Q,K,UTA,V,UPIDX) = L_SSCORRECTION
                 ENDDO
                ENDDO
               ENDDO
              ENDDO

            ENDIF

C  Ongrid output
C  -------------

C  just set to the cumulative source term 

          ELSE

            DO V = 1, N_GEOMETRIES
              DO Q = 1, K_PARAMETERS
                L_FINAL_SOURCE = L_SS_CUMSOURCE(Q,V)
                L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                PROFILEWF_SS(Q,K,UTA,V,UPIDX) = L_SSCORRECTION
              ENDDO
            ENDDO

          ENDIF

C  Check for updating the recursion

          IF ( NUT. NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT

C  end loop over optical depth

        ENDDO

C  end Upwelling clause

      ENDIF

C  ####################
C  #   DOWNWELLING    #
C  ####################

      IF ( DO_DNWELLING ) THEN

C  ===================================
C  Total phase function linearization (upwelling)
C  ===================================

        IF ( STERM_LAYERMASK_DN(K)) THEN

C  Loop over all geometries V 
C  Loop over varying parameters Q for layer K

          DO V = 1, N_GEOMETRIES
            DO Q = 1, K_PARAMETERS

C  Phase function moment variations
C    add TMS correction factor linearization

              IF ( DO_PHASFUNC_VARIATION(Q,K) ) THEN
                HELP = ZERO
                DO L = 0, NMOMENTS_INPUT
                  IF ( DO_SSCORR_TRUNCATION ) THEN
                    DNL1  = DBLE(2*L + 1 )
                    FDNL1 = SSFDEL(K) * DNL1
                    FACT  = ONE - SSFDEL(K)
                    GK11 = ( PHASMOMS_TOTAL_INPUT(L,K) - FDNL1 ) / FACT
                    HELP1 = L_PHASMOMS_TOTAL_INPUT(Q,L,K) * 
     &                        PHASMOMS_TOTAL_INPUT(L,K)
                    HELP2 = ( GK11 - DNL1 ) * L_SSFDEL(K,Q)
                    L_PHASMOM = ( HELP1 + HELP2 ) / FACT
                  ELSE
                    L_PHASMOM = L_PHASMOMS_TOTAL_INPUT(Q,L,K) *
     &                            PHASMOMS_TOTAL_INPUT(L,K)
                  ENDIF 
                  LEGPOLY = SS_PLEG_DN(V,K,L)
                  HELP = HELP + L_PHASMOM*LEGPOLY
                ENDDO
                L_EXACTSCAT_DN(V,K,Q) = HELP * TMS(K) +
     &                       EXACTSCAT_DN(V,K) * VAR_TMS(Q)
              ELSE
                L_EXACTSCAT_DN(V,K,Q) = EXACTSCAT_DN(V,K) * VAR_TMS(Q)
              ENDIF

C  end parameter loop and viewing directions loop

            ENDDO
          ENDDO

C  Only if layer K exists

        ENDIF

C  =====================================
C  Downwelling single scatter recurrence
C  =====================================

C  initialize cumulative source term

        NC = 0
        DO V = 1, N_GEOMETRIES
          DO Q = 1, K_PARAMETERS
            L_SS_CUMSOURCE(Q,V) = ZERO
          ENDDO
        ENDDO

C  initialise optical depth loop

        NSTART = 1
        NUT_PREV = NSTART - 1

C  Main loop over all output optical depths
C  ========================================

        DO UTA = 1, N_OUT_USERTAUS

C  Layer index for given optical depth

          NLEVEL = UTAU_LEVEL_MASK_DN(UTA)

C  finishing layer

          NUT = NLEVEL

C  Cumulative single scatter source terms to layer NUT
C  ---------------------------------------------------

C    1. Get layer source terms = Exact scattering * Multiplier
C    2. Loop over layers working upwards to NUT

          DO N = NSTART, NUT
            NC = N

C  If N = K (the layer that is varying)

            IF ( N .EQ. K ) THEN

              DO IB = 1, NBEAMS
                DO UM = 1, N_USER_STREAMS
                  DO IA = 1, N_USER_RELAZMS
                    V = UMOFF(IB,UM) + IA
                    DO Q = 1, K_PARAMETERS
                      L_SS_LAYERSOURCE = 
     &                    EXACTSCAT_DN(V,N)   * L_EMULT_DN(UM,N,K,IB,Q)
     &                + L_EXACTSCAT_DN(V,N,Q) *   EMULT_DN(UM,N,IB)
                      L_SS_LAYERSOURCE = FLUX_FACTOR * L_SS_LAYERSOURCE
                      L_SS_CUMSOURCE(Q,V) = L_SS_LAYERSOURCE
     &            +   T_DELT_USERM(N,UM)   * L_SS_CUMSOURCE(Q,V)
     &            + L_T_DELT_USERM(N,UM,Q) *   SS_CUMSOURCE_DN(V,NC-1)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO

C  Variations when N > K

           ELSE IF ( N .GT. K ) THEN

              DO IB = 1, NBEAMS
                DO UM = 1, N_USER_STREAMS
                  DO IA = 1, N_USER_RELAZMS
                    V = UMOFF(IB,UM) + IA
                    DO Q = 1, K_PARAMETERS
                      L_SS_LAYERSOURCE = 
     &                  EXACTSCAT_DN(V,N) * L_EMULT_DN(UM,N,K,IB,Q)
                      L_SS_LAYERSOURCE = FLUX_FACTOR * L_SS_LAYERSOURCE
                      L_SS_CUMSOURCE(Q,V) = L_SS_LAYERSOURCE
     &                  + T_DELT_USERM(N,UM) * L_SS_CUMSOURCE(Q,V)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO

C  Transmittance for N < K

            ELSE

              DO IB = 1, NBEAMS
                DO UM = 1, N_USER_STREAMS
                  DO IA = 1, N_USER_RELAZMS
                    V = UMOFF(IB,UM) + IA
                    DO Q = 1, K_PARAMETERS
                      L_SS_CUMSOURCE(Q,V) =
     &                  T_DELT_USERM(N,UM) * L_SS_CUMSOURCE(Q,V)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO

            ENDIF

          ENDDO

C  Offgrid output
C  --------------

C  add additional partial layer source term = Exact Scat * Multiplier
C  Set final cumulative source and Correct the intensity

          IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN

            UT = OFFGRID_UTAU_OUTINDEX(UTA)
            N  = OFFGRID_UTAU_LAYERIDX(UT)

C  If N = K (the layer that is varying)
C    add Linearization of additional partial layer source term =
C        L_Exact_Scat(n) * Multiplier  +  Exact_Scat * L_Multiplier(n)

            IF ( N .EQ. K ) THEN

              DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = UMOFF(IB,UM) + IA
                 DO Q = 1, K_PARAMETERS
                  L_SS_LAYERSOURCE = 
     &               EXACTSCAT_DN(V,N)   * L_UT_EMULT_DN(UM,UT,N,IB,Q)
     &           + L_EXACTSCAT_DN(V,N,Q) *   UT_EMULT_DN(UM,UT,IB)
                  L_SS_LAYERSOURCE = FLUX_FACTOR * L_SS_LAYERSOURCE
                  L_TR_CUMSOURCE = 
     &           + L_T_UTDN_USERM(UT,UM,Q) *   SS_CUMSOURCE_DN(V,NC)
     &           +   T_UTDN_USERM(UT,UM)   * L_SS_CUMSOURCE(Q,V)
                  L_FINAL_SOURCE   = L_TR_CUMSOURCE + L_SS_LAYERSOURCE
                  L_SSCORRECTION   = SSFLUX * L_FINAL_SOURCE
                  PROFILEWF_SS(Q,K,UTA,V,DNIDX) = L_SSCORRECTION
                 ENDDO
                ENDDO
               ENDDO
              ENDDO

C  Variations when N > K
C    add Linearization of additional partial layer source term =
C         Exact_Scat * L_Multiplier(k)

           ELSE IF ( N .GT. K ) THEN

              DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = UMOFF(IB,UM) + IA
                 DO Q = 1, K_PARAMETERS
                  L_SS_LAYERSOURCE = 
     &               EXACTSCAT_DN(V,N) * L_UT_EMULT_DN(UM,UT,K,IB,Q)
                  L_SS_LAYERSOURCE = FLUX_FACTOR * L_SS_LAYERSOURCE
                  L_TR_CUMSOURCE = 
     &               T_UTDN_USERM(UT,UM)   * L_SS_CUMSOURCE(Q,V)
                  L_FINAL_SOURCE   = L_TR_CUMSOURCE + L_SS_LAYERSOURCE
                  L_SSCORRECTION   = SSFLUX * L_FINAL_SOURCE
                  PROFILEWF_SS(Q,K,UTA,V,DNIDX) = L_SSCORRECTION
                 ENDDO
                ENDDO
               ENDDO
              ENDDO

C  Transmittance for the other layers

            ELSE

              DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = UMOFF(IB,UM) + IA
                 DO Q = 1, K_PARAMETERS
                  L_FINAL_SOURCE = 
     &               T_UTDN_USERM(UT,UM)   * L_SS_CUMSOURCE(Q,V)
                  L_SSCORRECTION   = SSFLUX * L_FINAL_SOURCE
                  PROFILEWF_SS(Q,K,UTA,V,DNIDX) = L_SSCORRECTION
                 ENDDO
                ENDDO
               ENDDO
              ENDDO

            ENDIF

C  Ongrid output
C  -------------

C  just set to the cumulative source term 

          ELSE

            DO V = 1, N_GEOMETRIES
              DO Q = 1, K_PARAMETERS
                L_FINAL_SOURCE = L_SS_CUMSOURCE(Q,V)
                L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                PROFILEWF_SS(Q,K,UTA,V,DNIDX) = L_SSCORRECTION
              ENDDO
            ENDDO

          ENDIF

C  Check for updating the recursion

          IF ( NUT. NE. NUT_PREV ) NSTART = NUT + 1
          NUT_PREV = NUT

C  end loop over optical depth

        ENDDO

C  end Downwelling clause

      ENDIF

C  debug

c      if ( k.eq.9) then
c      do v = 1, 8
c      write(99,*)v,PROFILEWF_SS(1,k,2,V,upIDX),PROFILEWF_SS(1,k,2,V,dnIDX) 
c      enddo
c      endif
  
C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_LC_SSCORR_NADIR
     &      ( SSFLUX, K_PARAMETERS)

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file of linearized input variables

      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'

C  include files of setup and multiplier variables (input)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  Include file of single scatter result variables

      INCLUDE '../includes/LIDORT_SINGSCAT.VARS'

C  include file with offset geometrical variables

      INCLUDE '../includes/LIDORT_RESULTS.VARS'

C  include files of linearized setup, solution and multiplier variables (input)

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_L_MULTIPLIERS.VARS'

C  Include file of results (linearized single scatter weighting funcs)

      INCLUDE '../includes/LIDORT_L_SINGSCAT.VARS'

C  Input arguments
C  ---------------

C  Flux factor

      DOUBLE PRECISION SSFLUX

C  Linearization control

      INTEGER          K_PARAMETERS

C  local variables
C  ---------------

C  indices

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL, V
      INTEGER          UT, L, UTA, UM, IB, IA, NC, Q, NM1, K

C  other help variables

      DOUBLE PRECISION HELP, HELP3, TTT, UVAR, AVAR, FT1, FT2
      DOUBLE PRECISION L_FINAL_SOURCE, L_SS_LAYERSOURCE
      DOUBLE PRECISION L_SSCORRECTION, L_TR_CUMSOURCE
      DOUBLE PRECISION VAR_TMS(MAX_ATMOSWFS,MAXLAYERS)
      DOUBLE PRECISION LEGPOLY, L_PHASMOM, DNM1
      DOUBLE PRECISION HELP1, HELP2, FDNL1, DNL1, FACT, GK11

C  Set up operations
C  -----------------

      DO K = 1, NLAYERS

C  Create Linearized TMS factors for all layers
C   ( Use UNSCALED linearized inputs - the Nakajima-Tanaka way)

        IF ( DO_SSCORR_TRUNCATION ) THEN
          TTT = TMS(K) /  ( ONE - SSFDEL(K) )
        ELSE
          TTT = TMS(K)
        ENDIF
 
        IF ( DO_DELTAM_SCALING ) THEN
          NM1 = NMOMENTS+1
          FT1 = TRUNC_FACTOR(K) * TTT
          FT2 = ONE + FT1
          DO Q = 1, K_PARAMETERS
            IF ( DO_PHASFUNC_VARIATION(Q,K) ) THEN
              UVAR = L_OMEGA_TOTAL_INPUT(Q,K)
              AVAR = UVAR + L_PHASMOMS_TOTAL_INPUT(Q,NM1,K)
              VAR_TMS(Q,K) = UVAR + FT1 * AVAR
            ELSE
              UVAR = L_OMEGA_TOTAL_INPUT(Q,K)
              VAR_TMS(Q,K) = UVAR * FT2
            ENDIF
          ENDDO
        ELSE
          DO Q = 1, K_PARAMETERS
            UVAR = L_OMEGA_TOTAL_INPUT(Q,K)
            VAR_TMS(Q,K) = UVAR
          ENDDO
        ENDIF

C  Additional Delta-M scaling
C  --------------------------

C  New section. R. Spurr, 07 September 2007.

        IF ( DO_SSCORR_TRUNCATION ) THEN
          NM1  = NMOMENTS_INPUT
          DNM1 = DBLE(2*NM1+1)
          DO Q = 1, K_PARAMETERS
            IF ( DO_PHASFUNC_VARIATION(Q,K) ) THEN
              FT1 = ONE / ( ONE- SSFDEL(K) )
              L_SSFDEL(K,Q) = PHASMOMS_TOTAL_INPUT(NM1,K) *
     *                    L_PHASMOMS_TOTAL_INPUT(Q,NM1,K) / DNM1
              VAR_TMS(Q,K) = VAR_TMS(Q,K) - FT1 * L_SSFDEL(K,Q)
            ENDIF
          ENDDO
        ENDIF

      ENDDO

C  ####################
C  #    UPWELLING     #
C  ####################

      IF ( DO_UPWELLING ) THEN

C  ===================================
C  Total phase function linearization (upwelling)
C  ===================================

C  All layers

        DO K = 1, NLAYERS
          IF ( STERM_LAYERMASK_UP(K)) THEN

C  Loop over all geometries V 
C  Loop over varying parameters Q for layer K

            DO V = 1, N_GEOMETRIES
              DO Q = 1, K_PARAMETERS

C  Phase function moment variations
C    add TMS correction factor linearization
C    No phase function variations - just add TMS linearization

                HELP3 = EXACTSCAT_UP(V,K) * VAR_TMS(Q,K)
                IF ( DO_PHASFUNC_VARIATION(Q,K) ) THEN
                  HELP = ZERO
                  DO L = 0, NMOMENTS_INPUT
                   IF ( DO_SSCORR_TRUNCATION ) THEN
                    DNL1  = DBLE(2*L + 1 )
                    FDNL1 = SSFDEL(K) * DNL1
                    FACT  = ONE - SSFDEL(K)
                    GK11 = ( PHASMOMS_TOTAL_INPUT(L,K) - FDNL1 ) / FACT
                    HELP1 = L_PHASMOMS_TOTAL_INPUT(Q,L,K) * 
     &                        PHASMOMS_TOTAL_INPUT(L,K)
                    HELP2 = ( GK11 - DNL1 ) * L_SSFDEL(K,Q)
                    L_PHASMOM = ( HELP1 + HELP2 ) / FACT
                   ELSE
                    L_PHASMOM = L_PHASMOMS_TOTAL_INPUT(Q,L,K) *
     &                            PHASMOMS_TOTAL_INPUT(L,K)
                   ENDIF 
                   LEGPOLY = SS_PLEG_UP(V,K,L)
                   HELP = HELP + L_PHASMOM*LEGPOLY
                  ENDDO
                  L_EXACTSCAT_UP(V,K,Q) = HELP * TMS(K) + HELP3
                ELSE
                  L_EXACTSCAT_UP(V,K,Q) = HELP3
                ENDIF

C  end parameter loop and viewing directions loop

              ENDDO
            ENDDO

C  Only if layer K exists

          ENDIF
        ENDDO

C  ===================================
C  Upwelling single scatter recurrence
C  ===================================

C  initialize cumulative source term

        NC = 0
        DO V = 1, N_GEOMETRIES
          DO Q = 1, K_PARAMETERS
            L_SS_CUMSOURCE(Q,V) = ZERO
          ENDDO
        ENDDO

C  initialise optical depth loop

        NSTART = NLAYERS
        NUT_PREV = NSTART + 1

C  Main loop over all output optical depths
C  ========================================

        DO UTA = N_OUT_USERTAUS, 1, -1

C  Layer index for given optical depth

          NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
          NUT = NLEVEL + 1

C  Cumulative single scatter source terms to layer NUT
C  ---------------------------------------------------

C    1. Get layer source terms = Exact scattering * Multiplier
C    2. Loop over layers working upwards to NUT

          DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N

C  ..All layers will have some variation

            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                  V = UMOFF(IB,UM) + IA
                  DO Q = 1, K_PARAMETERS
                    L_SS_LAYERSOURCE = 
     &                  EXACTSCAT_UP(V,N)   * L_EMULT_UP(UM,N,0,IB,Q)
     &              + L_EXACTSCAT_UP(V,N,Q) *   EMULT_UP(UM,N,IB)
                    L_SS_CUMSOURCE(Q,V) = L_SS_LAYERSOURCE
     &          +   T_DELT_USERM(N,UM)   * L_SS_CUMSOURCE(Q,V)
     &          + L_T_DELT_USERM(N,UM,Q) *   SS_CUMSOURCE_UP(V,NC-1)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO

C  Offgrid output
C  --------------

C  Set final cumulative source and Correct the Weighting function

          IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN

            UT = OFFGRID_UTAU_OUTINDEX(UTA)
            N  = OFFGRID_UTAU_LAYERIDX(UT)

C  ..(a) Multiplier for variations in the layer N
C    add Linearization of additional partial layer source term =
C        L_Exact_Scat(n) * Multiplier  +  Exact_Scat * L_Multiplier(n)

            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                  V = UMOFF(IB,UM) + IA
                  DO Q = 1, K_PARAMETERS
                    L_SS_LAYERSOURCE = 
     &               EXACTSCAT_UP(V,N)   * L_UT_EMULT_UP(UM,UT,0,IB,Q)
     &           + L_EXACTSCAT_UP(V,N,Q) *   UT_EMULT_UP(UM,UT,IB)
                    L_TR_CUMSOURCE = 
     &           + L_T_UTUP_USERM(UT,UM,Q) *   SS_CUMSOURCE_UP(V,NC)
     &           +   T_UTUP_USERM(UT,UM)   * L_SS_CUMSOURCE(Q,V)
                    L_FINAL_SOURCE   = L_TR_CUMSOURCE + L_SS_LAYERSOURCE
                    L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                    COLUMNWF_SS(Q,UTA,V,UPIDX) = L_SSCORRECTION
                  ENDDO
                ENDDO
              ENDDO
            ENDDO

C  Ongrid output
C  -------------

C  just set to the cumulative source term and Correct the Weighting Function

          ELSE

            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                  V = UMOFF(IB,UM) + IA
                  DO Q = 1, K_PARAMETERS
                    L_FINAL_SOURCE = L_SS_CUMSOURCE(Q,V)
                    L_SSCORRECTION = SSFLUX * L_SS_CUMSOURCE(Q,V)
                    COLUMNWF_SS(Q,UTA,V,UPIDX) = L_SSCORRECTION
                  ENDDO
                ENDDO
              ENDDO
            ENDDO

          ENDIF

C  Check for updating the recursion

          IF ( NUT. NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT

C  end loop over optical depth

        ENDDO

C  end Upwelling clause

      ENDIF

C  ####################
C  #   DOWNWELLING    #
C  ####################

      IF ( DO_DNWELLING ) THEN

C  ===================================
C  Total phase function linearization (downwelling)
C  ===================================

        DO K = 1, NLAYERS
          IF ( STERM_LAYERMASK_DN(K)) THEN

C  Loop over all geometries V 
C  Loop over varying parameters Q for layer K

            DO V = 1, N_GEOMETRIES
              DO Q = 1, K_PARAMETERS

C  Phase function moment variations
C    add TMS correction factor linearization

                HELP3 = EXACTSCAT_DN(V,K) * VAR_TMS(Q,K)
                IF ( DO_PHASFUNC_VARIATION(Q,K) ) THEN
                  HELP = ZERO
                  DO L = 0, NMOMENTS_INPUT
                   IF ( DO_SSCORR_TRUNCATION ) THEN
                    DNL1  = DBLE(2*L + 1 )
                    FDNL1 = SSFDEL(K) * DNL1
                    FACT  = ONE - SSFDEL(K)
                    GK11 = ( PHASMOMS_TOTAL_INPUT(L,K) - FDNL1 ) / FACT
                    HELP1 = L_PHASMOMS_TOTAL_INPUT(Q,L,K) * 
     &                        PHASMOMS_TOTAL_INPUT(L,K)
                    HELP2 = ( GK11 - DNL1 ) * L_SSFDEL(K,Q)
                    L_PHASMOM = ( HELP1 + HELP2 ) / FACT
                   ELSE
                    L_PHASMOM = L_PHASMOMS_TOTAL_INPUT(Q,L,K) *
     &                            PHASMOMS_TOTAL_INPUT(L,K)
                   ENDIF 
                   LEGPOLY = SS_PLEG_DN(V,K,L)
                   HELP = HELP + L_PHASMOM*LEGPOLY
                  ENDDO
                  L_EXACTSCAT_DN(V,K,Q) = HELP * TMS(K) + HELP3
                ELSE
                  L_EXACTSCAT_DN(V,K,Q) = HELP3
                ENDIF

C  end parameter loop and viewing directions loop

              ENDDO
            ENDDO

C  end layers 

          ENDIF
        ENDDO

C  =====================================
C  Downwelling single scatter recurrence
C  =====================================

C  initialize cumulative source term

        NC = 0
        DO V = 1, N_GEOMETRIES
          DO Q = 1, K_PARAMETERS
            L_SS_CUMSOURCE(Q,V) = ZERO
          ENDDO
        ENDDO

C  initialise optical depth loop

        NSTART = 1
        NUT_PREV = NSTART - 1

C  Main loop over all output optical depths
C  ========================================

        DO UTA = 1, N_OUT_USERTAUS

C  Layer index for given optical depth

          NLEVEL = UTAU_LEVEL_MASK_DN(UTA)

C  finishing layer

          NUT = NLEVEL

C  Cumulative single scatter source terms to layer NUT
C  ---------------------------------------------------

C    1. Get layer source terms = Exact scattering * Multiplier
C    2. Loop over layers working upwards to NUT

          DO N = NSTART, NUT
            NC = N
            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                  V = UMOFF(IB,UM) + IA
                  DO Q = 1, K_PARAMETERS
                    L_SS_LAYERSOURCE = 
     &                  EXACTSCAT_DN(V,N)   * L_EMULT_DN(UM,N,0,IB,Q)
     &              + L_EXACTSCAT_DN(V,N,Q) *   EMULT_DN(UM,N,IB)
                    L_SS_CUMSOURCE(Q,V) = L_SS_LAYERSOURCE
     &          +   T_DELT_USERM(N,UM)   * L_SS_CUMSOURCE(Q,V)
     &          + L_T_DELT_USERM(N,UM,Q) *   SS_CUMSOURCE_DN(V,NC-1)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO

C  Offgrid output
C  --------------

C  add additional partial layer source term = Exact Scat * Multiplier
C  Set final cumulative source and Correct the intensity

          IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN

            UT = OFFGRID_UTAU_OUTINDEX(UTA)
            N  = OFFGRID_UTAU_LAYERIDX(UT)

C    add Linearization of additional partial layer source term =
C        L_Exact_Scat(n) * Multiplier  +  Exact_Scat * L_Multiplier(n)

            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                  V = UMOFF(IB,UM) + IA
                  DO Q = 1, K_PARAMETERS
                    L_SS_LAYERSOURCE = 
     &               EXACTSCAT_DN(V,N)   * L_UT_EMULT_DN(UM,UT,0,IB,Q)
     &           + L_EXACTSCAT_DN(V,N,Q) *   UT_EMULT_DN(UM,UT,IB)
                    L_TR_CUMSOURCE = 
     &           + L_T_UTDN_USERM(UT,UM,Q) *   SS_CUMSOURCE_DN(V,NC)
     &           +   T_UTDN_USERM(UT,UM)   * L_SS_CUMSOURCE(Q,V)
                    L_FINAL_SOURCE   = L_TR_CUMSOURCE + L_SS_LAYERSOURCE
                    L_SSCORRECTION   = SSFLUX * L_FINAL_SOURCE
                    COLUMNWF_SS(Q,UTA,V,DNIDX) = L_SSCORRECTION
                  ENDDO
                ENDDO
              ENDDO
            ENDDO

C  Ongrid output
C  -------------

C  just set to the cumulative source term and Correct the Weighting Function

          ELSE

            DO V = 1, N_GEOMETRIES
              DO Q = 1, K_PARAMETERS
                L_FINAL_SOURCE = L_SS_CUMSOURCE(Q,V)
                L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                COLUMNWF_SS(Q,UTA,V,DNIDX) = L_SSCORRECTION
              ENDDO
            ENDDO

          ENDIF

C  Check for updating the recursion

          IF ( NUT. NE. NUT_PREV ) NSTART = NUT + 1
          NUT_PREV = NUT

C  end loop over optical depth

        ENDDO

C  end Downwelling clause

      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_LAP_DBCORRECTION(FLUXMULT)

C  Prepares Linearization of Exact Direct Beam reflection (BRDF case)
C   Linearization with respect to atmospheric profile variables

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file of linearized input variables

      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'

C  Include file of setup variables (Input to the present module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  Include files of Linearized setup variables

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'

C  Include files of Singlescatter and surface reflectance variables (Input)

      INCLUDE '../includes/LIDORT_SINGSCAT.VARS'
      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  Surface source output is stored in here

      INCLUDE '../includes/LIDORT_L_REFLECTANCE.VARS'

C  include file with offset geometrical variables

      INCLUDE '../includes/LIDORT_RESULTS.VARS'

C  Direct beam output stored in linearized correction file

      INCLUDE '../includes/LIDORT_L_SINGSCAT.VARS'

C  input argument

      DOUBLE PRECISION FLUXMULT

C  Local variables
C  ---------------

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER          UT, UTA, UM, UA, NC, IB, V, K, Q, KL
      DOUBLE PRECISION L_A_EXACTDB_SOURCE, L_ATTN, SUM, TR, LTR

C  first stage
C  -----------

C  return if no upwelling

      IF ( .NOT.DO_UPWELLING ) RETURN

C  Start weighting function loop

      DO K = 1, NLAYERS
        IF ( LAYER_VARY_FLAG(K) ) THEN
          DO Q = 1, LAYER_VARY_NUMBER(K)

C  Initialize the output results

            DO V = 1, N_GEOMETRIES
              L_EXACTDB_SOURCE(Q,K,V) = ZERO
              DO UTA = 1, N_OUT_USERTAUS
                PROFILEWF_DB(Q,K,UTA,V) = ZERO
              ENDDO
            ENDDO

C  For each solar beam
C    (a) Linearization factor L_ATTN for solar attenuation ATMOS_ATTN 
C    (b) for all geometreis, inner um over BRDF kernels
C   R. Spurr, 24 May 2005. RT Solutions Inc.
c            DO IB = 1, NBEAMS
c             IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
c              L_ATTN = - L_DELTAU_VERT(Q,K) * DELTAU_SLANT(NLAYERS,K,IB)   
c              DO UM = 1, N_USER_STREAMS
c               DO UA = 1, N_USER_RELAZMS
c                V = UMOFF(IB,UM) + UA
c                SUM = ZERO
c                DO KL = 1, N_BRDF_KERNELS
c                 L_A_EXACTDB_SOURCE = A_EXACTDB_SOURCE(V,KL)* L_ATTN
c                 SUM = SUM + L_A_EXACTDB_SOURCE
c                ENDDO
c                L_EXACTDB_SOURCE(Q,K,V) = SUM * FLUXMULT
c               ENDDO
c              ENDDO
c             ENDIF
c            ENDDO

C  New Code  R. Spurr, 6 August 2007. RT Solutions Inc.
C  ====================================================

C  Start geometry loops

            DO UM = 1, N_USER_STREAMS
             DO IB = 1, NBEAMS
              IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
               DO UA = 1, N_USER_RELAZMS
                V = UMOFF(IB,UM) + UA

c  Beam attenuation

                IF ( DO_SSCORR_OUTGOING ) THEN
                  L_ATTN = ZERO
                  IF (BOA_ATTN(V).GT.ZERO ) THEN
                    L_ATTN = L_BOA_ATTN(K,Q,V)/BOA_ATTN(V)
                  ENDIF
                ELSE
                  L_ATTN = - L_DELTAU_SLANT(Q,NLAYERS,K,IB)   
                ENDIF

C  Loop over albedo kernels, assign reflection

                SUM = ZERO
                DO KL = 1, N_BRDF_KERNELS
                  L_A_EXACTDB_SOURCE = A_EXACTDB_SOURCE(V,KL)* L_ATTN
                  SUM = SUM + L_A_EXACTDB_SOURCE
                ENDDO
                L_EXACTDB_SOURCE(Q,K,V) = SUM * FLUXMULT

C  Finish loops over geometries

               ENDDO
              ENDIF
             ENDDO
            ENDDO
       
C  Upwelling recurrence: transmittance of exact source term
C  --------------------------------------------------------

C  initialize cumulative source term
C    (Already flux-multiplied, because EXACTBD is a final result)

            NC =  0
            DO V = 1, N_GEOMETRIES
              L_DB_CUMSOURCE(V) = L_EXACTDB_SOURCE(Q,K,V)
            ENDDO

C  initialize optical depth loop

            NSTART = NLAYERS
            NUT_PREV = NSTART + 1

C  Main loop over all output optical depths

            DO UTA = N_OUT_USERTAUS, 1, -1

C  Layer index for given optical depth

              NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
              NUT    = NLEVEL + 1

C  Cumulative layer transmittance :
C    loop over layers working upwards to level NUT

              DO N = NSTART, NUT, -1
                NC = NLAYERS + 1 - N

C  If N  = K (varying layer), addition linearization of transmittance
C  If N ne K just transmitaance of the linearized source

                IF ( N.EQ.K) THEN
                 DO IB = 1, NBEAMS
                  DO UM = 1, N_USER_STREAMS
                   DO UA = 1, N_USER_RELAZMS
                    V = UMOFF(IB,UM) + UA
                    IF ( DO_SSCORR_OUTGOING ) THEN
                      TR  = UP_LOSTRANS(N,V)
                      LTR = L_UP_LOSTRANS(N,Q,V)
                    ELSE
                      TR  = T_DELT_USERM(N,UM)
                      LTR = L_T_DELT_USERM(N,UM,Q)
                    ENDIF
                    L_DB_CUMSOURCE(V) = TR * L_DB_CUMSOURCE(V)
     &                    + LTR * DB_CUMSOURCE(V,NC-1) 
                   ENDDO
                  ENDDO
                 ENDDO
                ELSE
                 DO IB = 1, NBEAMS
                  DO UM = 1, N_USER_STREAMS
                   DO UA = 1, N_USER_RELAZMS
                    V = UMOFF(IB,UM) + UA
                    IF ( DO_SSCORR_OUTGOING ) THEN
                      TR  = UP_LOSTRANS(N,V)
                    ELSE
                      TR  = T_DELT_USERM(N,UM)
                    ENDIF
                    L_DB_CUMSOURCE(V) = TR * L_DB_CUMSOURCE(V)
                   ENDDO
                  ENDDO
                 ENDDO
                ENDIF

C  end layer loop

              ENDDO

C  Offgrid output
C  ------------- 

C  Require partial layer transmittance
C  If N = K (varying layer), require additional linearization

              IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN

               UT = OFFGRID_UTAU_OUTINDEX(UTA)
               N  = OFFGRID_UTAU_LAYERIDX(UT)

               IF ( N.EQ.K) THEN
                DO IB = 1, NBEAMS
                 DO UM = 1, N_USER_STREAMS
                  DO UA = 1, N_USER_RELAZMS
                   V = UMOFF(IB,UM) + UA
                   IF ( DO_SSCORR_OUTGOING ) THEN
                     TR  = UP_LOSTRANS_UT(UT,V)
                     LTR = L_UP_LOSTRANS_UT(UT,Q,V)
                   ELSE
                     TR  = T_UTUP_USERM(UT,UM)
                     LTR = L_T_UTUP_USERM(UT,UM,Q)
                   ENDIF
                   PROFILEWF_DB(Q,K,UTA,V) = TR*L_DB_CUMSOURCE(V)
     &                                     + LTR*DB_CUMSOURCE(V,NC)
                  ENDDO
                 ENDDO
                ENDDO
               ELSE IF ( N.NE.K) THEN
                DO IB = 1, NBEAMS
                 DO UM = 1, N_USER_STREAMS
                  DO UA = 1, N_USER_RELAZMS
                   V = UMOFF(IB,UM) + UA
                   IF ( DO_SSCORR_OUTGOING ) THEN
                     TR  = UP_LOSTRANS_UT(UT,V)
                   ELSE
                     TR  = T_UTUP_USERM(UT,UM)
                   ENDIF
                   PROFILEWF_DB(Q,K,UTA,V) = TR*L_DB_CUMSOURCE(V)
                  ENDDO
                 ENDDO
                ENDDO
               ENDIF

C  Ongrid output : Set final cumulative source directly

              ELSE

               DO IB = 1, NBEAMS
                DO UM = 1, N_USER_STREAMS
                 DO UA = 1, N_USER_RELAZMS
                  V = UMOFF(IB,UM) + UA
                  PROFILEWF_DB(Q,K,UTA,V) = L_DB_CUMSOURCE(V)
                 ENDDO
                ENDDO
               ENDDO

              ENDIF

C  Check for updating the recursion 

             IF ( NUT. NE. NUT_PREV ) NSTART = NUT - 1
             NUT_PREV = NUT

C  end optical depth loop

            ENDDO

C  End weighting function loop

          ENDDO
        ENDIF
      ENDDO

C  debug

c      um = 99
c       do v = 1, 8
c        write(um,*)v, profilewf_db(1,12,2,v)
c       enddo

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_LAC_DBCORRECTION(FLUXMULT)

C  Prepares Linearization of Exact Direct Beam reflection (BRDF case)
C   Linearization with respect to bulk (total column) variables

C  New routine. 26 September 2006. R. Spurr. Integrated 14 May 2007.

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file of linearized input variables

      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'

C  Include file of setup variables (Input to the present module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  Include files of Linearized setup variables

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'

C  Include files of Singlescatter and surface reflectance variables (Input)

      INCLUDE '../includes/LIDORT_SINGSCAT.VARS'
      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  Surface source output is stored in here

      INCLUDE '../includes/LIDORT_L_REFLECTANCE.VARS'

C  include file with offset geometrical variables

      INCLUDE '../includes/LIDORT_RESULTS.VARS'

C  Direct beam output stored in linearized correction file

      INCLUDE '../includes/LIDORT_L_SINGSCAT.VARS'

C  input argument

      DOUBLE PRECISION FLUXMULT

C  Local variables
C  ---------------

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER          UT, UTA, UM, UA, NC, IB, V, K, Q, KL
      DOUBLE PRECISION L_A_EXACTDB_SOURCE, L_ATTN, SUM, TR, LTR

C  first stage
C  -----------

C  return if no upwelling

      IF ( .NOT.DO_UPWELLING ) RETURN

C  Start weighting function loop

      DO Q = 1, N_TOTALCOLUMN_WFS

C  Initialize the output results

        DO V = 1, N_GEOMETRIES
          L_EXACTDBC_SOURCE(Q,V) = ZERO
          DO UTA = 1, N_OUT_USERTAUS
            COLUMNWF_DB(Q,UTA,V)      = ZERO
          ENDDO
        ENDDO

C  For each solar beam
C    (a) Linearization factor L_ATTN for solar attenuation ATMOS_ATTN 
C    (b) for all geometreis, inner um over BRDF kernels
C   R. Spurr, 24 May 2005. RT Solutions Inc.
c        DO IB = 1, NBEAMS
c          IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
c            L_ATTN = ZERO
c            DO K = 1, NLAYERS
c              L_ATTN = L_ATTN - L_DELTAU_VERT(Q,K) * 
c     &                          DELTAU_SLANT(NLAYERS,K,IB)
c            ENDDO   
c            DO UM = 1, N_USER_STREAMS
c              DO UA = 1, N_USER_RELAZMS
c                V = UMOFF(IB,UM) + UA
c                SUM = ZERO
c                DO KL = 1, N_BRDF_KERNELS
c                  L_A_EXACTDB_SOURCE = A_EXACTDB_SOURCE(V,KL)*L_ATTN
c                  SUM = SUM + L_A_EXACTDB_SOURCE
c                ENDDO
c                L_EXACTDBC_SOURCE(Q,V) = SUM * FLUXMULT
c              ENDDO
c            ENDDO
c          ENDIF
c        ENDDO

C  New Code  R. Spurr, 6 August 2007. RT Solutions Inc.
C  ====================================================

C  Start geometry loops

            DO UM = 1, N_USER_STREAMS
             DO IB = 1, NBEAMS
              IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
               DO UA = 1, N_USER_RELAZMS
                V = UMOFF(IB,UM) + UA

c  Beam attenuation

                IF ( DO_SSCORR_OUTGOING ) THEN
                  L_ATTN = ZERO
                  IF (BOA_ATTN(V).GT.ZERO ) THEN
                    L_ATTN = L_BOA_ATTN(0,Q,V)/BOA_ATTN(V)
                  ENDIF
                ELSE
                  L_ATTN = ZERO
                  DO K = 1, NLAYERS
                    L_ATTN = L_ATTN - L_DELTAU_VERT(Q,K) * 
     &                                  DELTAU_SLANT(NLAYERS,K,IB)
                  ENDDO   
                ENDIF

C  Loop over albedo kernels, assign reflection

                SUM = ZERO
                DO KL = 1, N_BRDF_KERNELS
                  L_A_EXACTDB_SOURCE = A_EXACTDB_SOURCE(V,KL)* L_ATTN
                  SUM = SUM + L_A_EXACTDB_SOURCE
                ENDDO
                L_EXACTDBC_SOURCE(Q,V) = SUM * FLUXMULT

C  Finish loops over geometries

               ENDDO
              ENDIF
             ENDDO
            ENDDO
       
C  Upwelling recurrence: transmittance of exact source term
C  --------------------------------------------------------

C  initialize cumulative source term
C    (Already flux-multiplied, because EXACTBD is a final result)

        NC =  0
        DO V = 1, N_GEOMETRIES
          L_DB_CUMSOURCE(V) = L_EXACTDBC_SOURCE(Q,V)
        ENDDO

C  initialize optical depth loop

        NSTART = NLAYERS
        NUT_PREV = NSTART + 1

C  Main loop over all output optical depths

        DO UTA = N_OUT_USERTAUS, 1, -1

C  Layer index for given optical depth

          NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
          NUT    = NLEVEL + 1

C  Cumulative layer transmittance :
C    loop over layers working upwards to level NUT

          DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                DO UA = 1, N_USER_RELAZMS
                  V = UMOFF(IB,UM) + UA
                  IF ( DO_SSCORR_OUTGOING ) THEN
                    TR  = UP_LOSTRANS(N,V)
                    LTR = L_UP_LOSTRANS(N,Q,V)
                  ELSE
                    TR  = T_DELT_USERM(N,UM)
                    LTR = L_T_DELT_USERM(N,UM,Q)
                  ENDIF
                  L_DB_CUMSOURCE(V) = TR * L_DB_CUMSOURCE(V)
     &                  + LTR * DB_CUMSOURCE(V,NC-1) 
                ENDDO
              ENDDO
            ENDDO
          ENDDO

C  Offgrid output, Require partial layer transmittance

          IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN

            UT = OFFGRID_UTAU_OUTINDEX(UTA)
            N  = OFFGRID_UTAU_LAYERIDX(UT)
            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                DO UA = 1, N_USER_RELAZMS
                  V = UMOFF(IB,UM) + UA
                  IF ( DO_SSCORR_OUTGOING ) THEN
                    TR  = UP_LOSTRANS_UT(UT,V)
                    LTR = L_UP_LOSTRANS_UT(UT,Q,V)
                  ELSE
                    TR  = T_UTUP_USERM(UT,UM)
                    LTR = L_T_UTUP_USERM(UT,UM,Q)
                  ENDIF
                  COLUMNWF_DB(Q,UTA,V) = TR*L_DB_CUMSOURCE(V)
     &                                   + LTR*DB_CUMSOURCE(V,NC)
                ENDDO
              ENDDO
            ENDDO

C  Ongrid output : Set final cumulative source directly

          ELSE

            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                DO UA = 1, N_USER_RELAZMS
                  V = UMOFF(IB,UM) + UA
                  COLUMNWF_DB(Q,UTA,V) = L_DB_CUMSOURCE(V)
                ENDDO
              ENDDO
            ENDDO

          ENDIF

C  Check for updating the recursion 

          IF ( NUT. NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT

C  end optical depth loop

        ENDDO

C  End weighting function loop

      ENDDO

C  debug

c      um = 99
c       do v = 1, 8
c        write(99,*)v, columnwf_db(1,2,v)
c       enddo

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_LS_DBCORRECTION(FLUXMULT)

C  Prepares Linearization of Exact Direct Beam reflection (BRDF case)
C   Linearization with respect to kernel amplitude and parameter variables

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include file of linearized input variables

      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'

C  Include file of setup variables (Input to the present module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  Include files of Linearized setup variables

      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'

C  Include files of Singlescatter and surface reflectance variables (Input)

      INCLUDE '../includes/LIDORT_SINGSCAT.VARS'
      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  Surface source output is stored in here

      INCLUDE '../includes/LIDORT_L_REFLECTANCE.VARS'

C  Direct beam output stored in linearized correction file

      INCLUDE '../includes/LIDORT_L_SINGSCAT.VARS'

C  input argument

      DOUBLE PRECISION FLUXMULT

C  Local variables
C  ---------------

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER          UT, UTA, UM, UA, NC, IB, V, WF
      INTEGER          PAR_INDEX, BRDF_INDEX
      DOUBLE PRECISION REFL_ATTN, SUM, L_KER, PAR, TR

C  first stage
C  -----------

C  initialise output

      DO WF = 1, N_TOTALBRDF_WFS
        DO V = 1, N_GEOMETRIES
          LS_EXACTDB_SOURCE(WF,V) = ZERO
          DO UTA = 1, N_OUT_USERTAUS
            SURFACEWF_DB(WF,UTA,V)      = ZERO
          ENDDO
        ENDDO
      ENDDO

C  return if no upwelling

      IF ( .NOT.DO_UPWELLING ) RETURN

C  =========================================
C  1. Get the source function linearizations
C  =========================================

C  Initialise WF count

      WF = 0

C  Start loop over BRDF Kernels

      DO BRDF_INDEX = 1, N_BRDF_KERNELS

C  Amplitude weighting functions
C  -----------------------------

C  For each geometry, make Linearization w.r.t surface amplitude
C   R. Spurr, 24 May 2005. RT Solutions Inc.

        IF ( DO_KERNEL_FACTOR_WFS(BRDF_INDEX) ) THEN
          WF = WF + 1
          DO IB = 1, NBEAMS
           IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
            DO UM = 1, N_USER_STREAMS
             DO UA = 1, N_USER_RELAZMS
              V = UMOFF(IB,UM) + UA
              SUM= A_EXACTDB_SOURCE(V,BRDF_INDEX)
              LS_EXACTDB_SOURCE(WF,V) = SUM * FLUXMULT
             ENDDO
            ENDDO
           ENDIF
          ENDDO
        ENDIF

C  Parameter weighting functions
C  -----------------------------

C  For each geometry, make Linearization w.r.t nonlinear parameter
C   R. Spurr, 24 May 2005. RT Solutions Inc.

        DO PAR_INDEX = 1, N_BRDF_PARAMETERS(BRDF_INDEX)
         IF (DO_KERNEL_PARAMS_WFS(BRDF_INDEX,PAR_INDEX)) THEN
          WF = WF + 1
          PAR = BRDF_PARAMETERS(BRDF_INDEX,PAR_INDEX)
          DO UM = 1, N_USER_STREAMS
           DO IB = 1, NBEAMS
            IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
             DO UA = 1, N_USER_RELAZMS
              V = UMOFF(IB,UM) + UA
              REFL_ATTN =  BRDF_FACTORS(BRDF_INDEX) * ATTN_DB_SAVE(V)
              L_KER = D_EXACTDB_BRDFUNC(BRDF_INDEX,PAR_INDEX,UM,UA,IB)
              SUM = FLUXMULT * L_KER * PAR
              LS_EXACTDB_SOURCE(WF,V) = REFL_ATTN * SUM
             ENDDO
            ENDIF
           ENDDO
          ENDDO
         ENDIF
        ENDDO

C  End kernel loop

      ENDDO

C  ====================================================
C  2. Transmittance of source term: upwelling recursion
C  ====================================================

      DO WF = 1, N_TOTALBRDF_WFS

C  initialize cumulative source term
C    (Already flux-multiplied, because EXACTBD is a final result)

        NC =  0
        DO V = 1, N_GEOMETRIES
          LS_DB_CUMSOURCE(V) = LS_EXACTDB_SOURCE(WF,V)
        ENDDO

C  initialize optical depth loop

        NSTART = NLAYERS
        NUT_PREV = NSTART + 1

C  Main loop over all output optical depths

        DO UTA = N_OUT_USERTAUS, 1, -1

C  Layer index for given optical depth

          NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
          NUT    = NLEVEL + 1

C  Cumulative layer transmittance :
C    loop over layers working upwards to level NUT

          DO N = NSTART, NUT, -1
           NC = NLAYERS + 1 - N
           DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
             DO UA = 1, N_USER_RELAZMS
              V = UMOFF(IB,UM) + UA
              IF ( DO_SSCORR_OUTGOING ) THEN
                TR  = UP_LOSTRANS(N,V)
              ELSE
                TR  = T_DELT_USERM(N,UM)
              ENDIF
              LS_DB_CUMSOURCE(V) = TR * LS_DB_CUMSOURCE(V)
             ENDDO
            ENDDO
           ENDDO
          ENDDO

C  Offgrid output
C  -------------- 

C  Require partial layer transmittance

          IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN
           UT = OFFGRID_UTAU_OUTINDEX(UTA)
           N  = OFFGRID_UTAU_LAYERIDX(UT)
           DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
             DO UA = 1, N_USER_RELAZMS
              V = UMOFF(IB,UM) + UA
              IF ( DO_SSCORR_OUTGOING ) THEN
                TR  = UP_LOSTRANS_UT(UT,V)
              ELSE
                TR  = T_UTUP_USERM(UT,UM)
              ENDIF
              SURFACEWF_DB(WF,UTA,V) = TR * LS_DB_CUMSOURCE(V)
             ENDDO
            ENDDO
           ENDDO

C  Ongrid output : Set final cumulative source directly

          ELSE
           DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
             DO UA = 1, N_USER_RELAZMS
              V = UMOFF(IB,UM) + UA
              SURFACEWF_DB(WF,UTA,V) = LS_DB_CUMSOURCE(V)
             ENDDO
            ENDDO
           ENDDO
          ENDIF

C  Check for updating the recursion 

          IF ( NUT. NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT

C  end optical depth loop

        ENDDO

C  End Loop over all Surface  weighting functions

      ENDDO

C  debug

c      um = 99
c      do v = 1, 8
c       write(99,*)v, SURFACEWF_DB(3,2,v) 
c      enddo

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_L_SSCORR_OUTGOING
     &        ( SSFLUX, FAIL, MESSAGE)

C  Single scatter exact calculation for the outgoing LOS
C   This one with optional linearizations
C         - NEW for version 3.2

C   Programmed by R. Spurr, RT Solutions Inc.
C    First Draft, January 23rd 2007
C   Validated against TOMRAD, 29 March 2007.
C   Partial layer output added September 2007

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include file of input variables
C  Include file of bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  include files of setup and multiplier variables (input)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_MULTIPLIERS.VARS'

C  Include file of single scatter result variables

      INCLUDE '../includes/LIDORT_SINGSCAT.VARS'

C  include files of linearized inputs/setups/solution/multiplier variables

      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/LIDORT_L_MULTIPLIERS.VARS'

C  Include file of results (linearized single scatter weighting funcs)

      INCLUDE '../includes/LIDORT_L_SINGSCAT.VARS'

C  Arguments
C  ---------

      DOUBLE PRECISION SSFLUX
      logical          fail
      character*(*)    message

C  Geometry routine inputs and outputs
C  -----------------------------------

      logical          do_fine
      LOGICAL          do_partials
      double precision alpha_boa, theta_boa, phi_boa

C  main outputs (geometry)

      integer          ntraverse(0:maxlayers)
      double precision sunpaths(0:maxlayers,maxlayers)
      double precision radii   (0:maxlayers)
      double precision alpha_all  (0:maxlayers)

C  Fine level output (geometry)

      integer          ntraverse_fine(maxlayers,maxfinelayers)
      double precision sunpaths_fine (maxlayers,maxlayers,maxfinelayers)
      double precision radii_fine    (maxlayers,maxfinelayers)
      double precision alpha_fine    (maxlayers,maxfinelayers)

C  Partial layer output (geometry)

      integer          ntraverse_ut(max_offgrid_usertaus)
      double precision sunpaths_ut (maxlayers,max_offgrid_usertaus)
      double precision radii_ut    (max_offgrid_usertaus)
      double precision alpha_ut    (max_offgrid_usertaus)

C  Other (incidental) geometrical output

      double precision lospaths(maxlayers)
      double precision lospaths_ut_up(max_offgrid_usertaus)
      double precision lospaths_ut_dn(max_offgrid_usertaus)
      double precision theta_all  (0:maxlayers)
      double precision phi_all    (0:maxlayers)
      double precision cosscat_up (0:maxlayers)
      double precision cosscat_dn (0:maxlayers)

C  Extinction and var tms

      double precision extinction   (maxlayers)
      double precision L_extinction (maxlayers,max_atmoswfs)
      double precision var_tms      (maxlayers,max_atmoswfs)

C  Partial layer heights

      double precision height_grid_ut(max_offgrid_usertaus)

C  Bookkeeping

      integer          NTYPE_VARY
      logical          DO_RTSOL_VARY(maxlayers)
      integer          NPARAMS_VARY(maxlayers) 
      integer          KINDEX(maxlayers)

C  local variables
C  ---------------

C  Indices

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL, L
      INTEGER          UT, UTA, UM, IA, NC, IB, V, NM1
      integer          k, k_parameters, ks, q

C  help variables (double precision)

      DOUBLE PRECISION FINAL_SOURCE, SS_LAYERSOURCE, SS_CUMSOURCE
      DOUBLE PRECISION SSCORRECTION, COSSCAT, LEGPOLY, VAR1, TRANS
      DOUBLE PRECISION DF1(0:MAXMOMENTS_INPUT)
      DOUBLE PRECISION DF2(0:MAXMOMENTS_INPUT)
      DOUBLE PRECISION UVAR, AVAR, FT1, FT2, LSS, L_PHASMOM
      DOUBLE PRECISION L_FINAL_SOURCE, L_SSCORRECTION, DNM1, XT
      DOUBLE PRECISION FACT, DNL1, FDNL1, GK11, HELP, HELP1, HELP2

C  Set up operations
C  -----------------

C  Set up partials flag

      do_fine = .true.
      do_partials = ( n_offgrid_usertaus .gt. 0 )

C  Floating point numbers for Legendre polynomials

      DO L = 2, NMOMENTS_INPUT
        HELP = DBLE(L)
        DF1(L) = DBLE(2*L-1)/HELP
        DF2(L) = DBLE(L-1)/HELP
      ENDDO

C  Create TMS factors, these get stored
C    Delta-M Scaling introduced April 2005.

      DO N = 1, NLAYERS
        IF ( DO_DELTAM_SCALING ) THEN
          HELP   = ONE - TRUNC_FACTOR(N) * OMEGA_TOTAL_INPUT(N)
          TMS(N) = OMEGA_TOTAL_INPUT(N) / HELP
        ELSE
          TMS(N) = OMEGA_TOTAL_INPUT(N)
        ENDIF
      ENDDO

C  linearizations of the TMS factors
C  ---------------------------------

C  Create Linearized TMS factor for each layer
C   ( Use UNSCALED linearized inputs - the Nakajima-Tanaka way)
C    Distinguish between Deltam case or not.
C  Only if required

      IF ( DO_ATMOS_LINEARIZATION ) THEN

C  layer loop

       DO K = 1, NLAYERS
        IF ( LAYER_VARY_FLAG(K) ) THEN
         K_PARAMETERS = LAYER_VARY_NUMBER(K)

C  Deltam_scaling. Extra contributions if phase function moments are varying.
    
         IF ( DO_DELTAM_SCALING ) THEN
          NM1 = NMOMENTS+1
          FT1 = TRUNC_FACTOR(K) * TMS(K)
          FT2 = ONE + FT1
          DO Q = 1, K_PARAMETERS
           IF ( DO_PHASFUNC_VARIATION(Q,K) ) THEN
            UVAR = L_OMEGA_TOTAL_INPUT(Q,K)
            AVAR = UVAR + L_PHASMOMS_TOTAL_INPUT(Q,NM1,K)
            VAR_TMS(K,Q) = UVAR + FT1 * AVAR
           ELSE
            UVAR = L_OMEGA_TOTAL_INPUT(Q,K)
            VAR_TMS(K,Q) = UVAR * FT2
           ENDIF
          ENDDO
         ENDIF

C  No delta-M scaling, just copy

         IF ( .NOT. DO_DELTAM_SCALING ) THEN
          DO Q = 1, K_PARAMETERS
           UVAR = L_OMEGA_TOTAL_INPUT(Q,K)
           VAR_TMS(K,Q) = UVAR
          ENDDO
         ENDIF

C  End layer loop

        ENDIF
       ENDDO

C  End linearization clause

      ENDIF

C  Additional Delta-M scaling
C  --------------------------

C  New section. R. Spurr, 07 September 2007.
C   TMS gets modified by (1-F). Save the truncation factor.
C   Phase function moments are modified later on.

      IF ( DO_SSCORR_TRUNCATION ) THEN

C  Basic scaling

        NM1  = NMOMENTS_INPUT
        DNM1 = DBLE(2*NM1+1)
        DO N = 1, NLAYERS
          SSFDEL(N) = PHASMOMS_TOTAL_INPUT(NM1,N) / DNM1
          TMS(N) = TMS(N) * ( ONE - SSFDEL(N) )
        ENDDO

C  Linearization, Only change if SCATMAT variation is set

        IF ( DO_ATMOS_LINEARIZATION ) THEN
          DO K = 1, NLAYERS
            IF ( LAYER_VARY_FLAG(K) ) THEN
              K_PARAMETERS = LAYER_VARY_NUMBER(K)
              DO Q = 1, K_PARAMETERS
                IF ( DO_PHASFUNC_VARIATION(Q,K) ) THEN
                  FT1 = ONE / ( ONE- SSFDEL(K) )
                  L_SSFDEL(K,Q) = PHASMOMS_TOTAL_INPUT(NM1,K) *
     *                    L_PHASMOMS_TOTAL_INPUT(Q,NM1,K) / DNM1
                  VAR_TMS(K,Q) = VAR_TMS(K,Q) - FT1 * L_SSFDEL(K,Q)
                ENDIF
              ENDDO
            ENDIF
          ENDDO
        ENDIF

      ENDIF

C  Parameter control
C  -----------------

      IF ( DO_PROFILE_LINEARIZATION ) THEN
        NTYPE_VARY = NLAYERS
        DO N = 1, NLAYERS
         DO_RTSOL_VARY(n) = LAYER_VARY_FLAG(n)
         NPARAMS_VARY(n)  = LAYER_VARY_NUMBER(n)
         KINDEX(N) = N
        ENDDO
      ELSE IF ( DO_COLUMN_LINEARIZATION ) THEN
        NTYPE_VARY = 1
        DO_RTSOL_VARY(1) = .TRUE.
        NPARAMS_VARY(1)  = N_TOTALCOLUMN_WFS
        KINDEX(1) = 0
      ENDIF

C  Create extinctions
C  ------------------

C  Use basic definitions

      DO N = 1, NLAYERS
        HELP = HEIGHT_GRID(N-1) - HEIGHT_GRID(N)
        EXTINCTION(N) = DELTAU_VERT(N) / HELP
      ENDDO

C  Linearized extinctions

      IF ( DO_ATMOS_LINEARIZATION ) THEN
       DO K = 1, NLAYERS
        IF ( LAYER_VARY_FLAG(K) ) THEN
         K_PARAMETERS = LAYER_VARY_NUMBER(K)
         DO Q = 1, K_PARAMETERS
           L_EXTINCTION(K,Q) = L_DELTAU_VERT(Q,K) * EXTINCTION(K)
         ENDDO
        ENDIF
       ENDDO
      ENDIF

C  Create the heights of partial layers

      IF ( DO_PARTIALS ) THEN
        do ut = 1, n_offgrid_usertaus
          n = offgrid_utau_layeridx(ut)
          xt = deltau_vert(n) - offgrid_utau_values(ut)
          height_grid_ut(ut) = height_grid(n) + xt / extinction(n)
        enddo
      endif

C  Source function contributions
C  =============================

C  Start the main loop over all solar and viewing geometries
C    old code, not adjusted
c      DO IB = 1, NBEAMS
c        THETA_BOA = BEAM_SZAS(IB)
c        DO UM = 1, N_USER_STREAMS
c          ALPHA_BOA = USER_ANGLES_INPUT(UM)
c          DO IA = 1, N_USER_RELAZMS
c            PHI_BOA = USER_RELAZMS(IA)
c            V = UMOFF(IB,UM) + IA

C  Start the main loop over all solar and viewing geometries

      DO UM = 1, N_USER_STREAMS
        ALPHA_BOA = USER_ANGLES_ADJUST(UM)
        DO IB = 1, NBEAMS
          DO IA = 1, N_USER_RELAZMS
            THETA_BOA = BEAM_SZAS_ADJUST(UM,IB,IA)
            PHI_BOA   = USER_RELAZMS_ADJUST(UM,IB,IA)
            V = UMOFF(IB,UM) + IA
    
C  Call to geometry

            call outgoing_sphergeom_fine
     i  ( maxlayers, maxfinelayers, max_offgrid_usertaus,
     i    do_fine, do_partials, nlayers, nfinelayers,
     i    n_offgrid_usertaus, offgrid_utau_layeridx,
     i    height_grid, height_grid_ut, earth_radius,
     i    alpha_boa, theta_boa, phi_boa,
     o    sunpaths,      radii,      ntraverse,      alpha_all, 
     o    sunpaths_fine, radii_fine, ntraverse_fine, alpha_fine,
     o    sunpaths_ut,   radii_ut,   ntraverse_ut,   alpha_ut,
     o    offgrid_utau_layerfineidx, lospaths,
     o    lospaths_ut_up, lospaths_ut_dn, 
     o    theta_all, phi_all, cosscat_up, cosscat_dn,
     o    fail, message )

C  debug geometry
c            do n = 0, nlayers
c              write(45,'(i4,101f10.5)')n,(sunpaths(n,v),v=1,nlayers)
c            enddo
c            pause

C  Upwelling calculation
C  ---------------------

            IF ( DO_UPWELLING ) THEN

C  Multipliers and transmittances + linearizations

              call l_outgoing_integration_up
     i           ( nlayers, nfinelayers,
     i             do_fine, do_partials,
     i     do_profile_linearization, layer_vary_flag, layer_vary_number,
     i     do_column_linearization, n_totalcolumn_wfs,
     i             extinction, l_extinction,
     i             n_offgrid_usertaus, offgrid_utau_layeridx,
     i             offgrid_utau_layerfineidx,
     i             sunpaths, radii, ntraverse, alpha_all,
     i             sunpaths_ut,   ntraverse_ut,   alpha_ut,
     i             sunpaths_fine, ntraverse_fine, alpha_fine,
     o             up_multipliers(1,v), up_lostrans(1,v), boa_attn(v),
     o             up_multipliers_ut(1,v), up_lostrans_ut(1,v),
     o             l_up_multipliers(1,0,1,v),
     o             l_up_lostrans(1,1,v),
     o             l_boa_attn(0,1,v),
     o             l_up_multipliers_ut(1,0,1,v),
     o             l_up_lostrans_ut(1,1,v) )

C  DEbug
c              write(35,*)v
c              do n = 1, nlayers
c                write(35,'(1p2e18.10)')
c     &        up_lostrans(n,v),up_multipliers(n,v)
c              enddo
c              if ( v.eq.8)pause

c              if ( do_column_linearization ) then
c              do n = 1, nlayers
c          write(36,'(2i4,1p20e18.10)')v,n,
c     &       up_lostrans(n,v),up_multipliers(n,v),
c     &     l_up_lostrans(n,2,v),(l_up_multipliers(n,0,2,v))
c              enddo
c              else
c             do n = 1, nlayers
c               write(35,'(2i4,1p2e18.10)')v,n,
c     &       up_lostrans(n,v),up_multipliers(n,v)
c                enddo
c               endif

C  legendre polynomials

              COSSCAT = COSSCAT_UP(NLAYERS)
	      SS_PLEG_UP(V,1,0) = ONE
	      SS_PLEG_UP(V,1,1) = COSSCAT
	      DO L = 2, NMOMENTS_INPUT
	        SS_PLEG_UP(V,1,L) =
     &           DF1(L) * SS_PLEG_UP(V,1,L-1) * COSSCAT  -
     &           DF2(L) * SS_PLEG_UP(V,1,L-2)
              ENDDO

C  Phase functions (multiplied by TMS factor). Save them.

              DO N = 1, NLAYERS
                IF ( STERM_LAYERMASK_UP(N) ) THEN
                  HELP = ZERO
                  DO L = 0, NMOMENTS_INPUT
                   IF ( DO_SSCORR_TRUNCATION ) THEN
                    DNL1  = DBLE(2*L + 1 )
                    FDNL1 = SSFDEL(N) * DNL1
                    FACT  = ONE - SSFDEL(N)
                    GK11 = ( PHASMOMS_TOTAL_INPUT(L,N) - FDNL1 ) / FACT
                   ELSE
                    GK11 = PHASMOMS_TOTAL_INPUT(L,N)
                   ENDIF 
	           LEGPOLY = SS_PLEG_UP(V,1,L)
	           HELP = HELP + GK11 * LEGPOLY
	          ENDDO
	          EXACTSCAT_UP(V,N) = HELP * TMS(N)
                ENDIF
c                write(88,*)v,n,exactscat_up(v,n)
              ENDDO
c              if (v.eq.8)pause

C  Linearized phase functions
C   --extra terms for Phase function moment variations
C   -- must add TMS correction factor linearization

              IF ( DO_ATMOS_LINEARIZATION ) THEN
               DO K = 1, NLAYERS
                IF ( STERM_LAYERMASK_UP(K)) THEN
                 IF ( LAYER_VARY_FLAG(K) ) THEN
                  K_PARAMETERS = LAYER_VARY_NUMBER(K)
                  DO Q = 1, K_PARAMETERS
                   VAR1 = EXACTSCAT_UP(V,K) * VAR_TMS(K,Q)
                   IF ( DO_PHASFUNC_VARIATION(Q,K) ) THEN
                    HELP = ZERO
                    DO L = 0, NMOMENTS_INPUT
                     IF ( DO_SSCORR_TRUNCATION ) THEN
                      DNL1  = DBLE(2*L + 1 )
                      FDNL1 = SSFDEL(K) * DNL1
                      FACT  = ONE - SSFDEL(K)
                      GK11 = (PHASMOMS_TOTAL_INPUT(L,K)-FDNL1)/FACT
                      HELP1 = L_PHASMOMS_TOTAL_INPUT(Q,L,K) * 
     &                          PHASMOMS_TOTAL_INPUT(L,K)
                      HELP2 = ( GK11 - DNL1 ) * L_SSFDEL(K,Q)
                      L_PHASMOM = ( HELP1 + HELP2 ) / FACT
                     ELSE
                      L_PHASMOM = L_PHASMOMS_TOTAL_INPUT(Q,L,K) *
     &                              PHASMOMS_TOTAL_INPUT(L,K)
                     ENDIF 
                     LEGPOLY = SS_PLEG_UP(V,1,L)
                     HELP = HELP + L_PHASMOM*LEGPOLY
                    ENDDO
                    L_EXACTSCAT_UP(V,K,Q) = HELP * TMS(K) + VAR1             
                   ELSE
                    L_EXACTSCAT_UP(V,K,Q) = VAR1
                   ENDIF
                  ENDDO
                 ENDIF
                ENDIF 
               ENDDO
              ENDIF

C  End upwelling clause

            ENDIF

C  Downwelling calculation
C  -----------------------

            IF ( DO_DNWELLING ) THEN

C  Multipliers, transmittances + linearizations

              call l_outgoing_integration_dn
     i           ( nlayers, nfinelayers,
     i             do_fine, do_partials,
     i     do_profile_linearization, layer_vary_flag, layer_vary_number,
     i     do_column_linearization, n_totalcolumn_wfs,
     i             extinction, l_extinction,
     i             n_offgrid_usertaus, offgrid_utau_layeridx,
     i             offgrid_utau_layerfineidx,
     i             sunpaths, radii, ntraverse, alpha_all,
     i             sunpaths_ut,   ntraverse_ut,   alpha_ut,
     i             sunpaths_fine, ntraverse_fine, alpha_fine,
     o             dn_multipliers(1,v), dn_lostrans(1,v),
     o             dn_multipliers_ut(1,v), dn_lostrans_ut(1,v),
     o             l_dn_multipliers(1,0,1,v),
     o             l_dn_lostrans(1,1,v),
     o             l_dn_multipliers_ut(1,0,1,v),
     o             l_dn_lostrans_ut(1,1,v) )

C  Debug
c              do n = 1, nlayers
c                write(88,*)n,dn_lostrans(n,1),dn_multipliers(n,1)
c              enddo

C  Legendre polynomials

              COSSCAT = COSSCAT_DN(NLAYERS)
	      SS_PLEG_DN(V,1,0) = ONE
	      SS_PLEG_DN(V,1,1) = COSSCAT
	      DO L = 2, NMOMENTS_INPUT
	        SS_PLEG_DN(V,1,L) =
     &           DF1(L) * SS_PLEG_DN(V,1,L-1) * COSSCAT  -
     &           DF2(L) * SS_PLEG_DN(V,1,L-2)
              ENDDO

C  Phase functions (multiplied by TMS factor). Save them.

              DO N = 1, NLAYERS
                IF ( STERM_LAYERMASK_DN(N) ) THEN
                  HELP = ZERO
                  DO L = 0, NMOMENTS_INPUT
                   IF ( DO_SSCORR_TRUNCATION ) THEN
                    DNL1  = DBLE(2*L + 1 )
                    FDNL1 = SSFDEL(N) * DNL1
                    FACT  = ONE - SSFDEL(N)
                    GK11 = ( PHASMOMS_TOTAL_INPUT(L,N) - FDNL1 ) / FACT
                   ELSE
                    GK11 = PHASMOMS_TOTAL_INPUT(L,N)
                   ENDIF 
	           LEGPOLY = SS_PLEG_DN(V,1,L)
	           HELP = HELP + GK11 * LEGPOLY
	          ENDDO
	          EXACTSCAT_DN(V,N) = HELP * TMS(N)
                ENDIF
              ENDDO

C  Linearized phase functions
C   --extra terms for Phase function moment variations
C   -- must add TMS correction factor linearization

              IF ( DO_ATMOS_LINEARIZATION ) THEN
               DO K = 1, NLAYERS
                IF ( STERM_LAYERMASK_DN(K)) THEN
                 IF ( LAYER_VARY_FLAG(K) ) THEN
                  K_PARAMETERS = LAYER_VARY_NUMBER(K)
                  DO Q = 1, K_PARAMETERS
                   VAR1 = EXACTSCAT_DN(V,K) * VAR_TMS(K,Q)
                   IF ( DO_PHASFUNC_VARIATION(Q,K) ) THEN
                    HELP = ZERO
                    DO L = 0, NMOMENTS_INPUT
                     IF ( DO_SSCORR_TRUNCATION ) THEN
                      DNL1  = DBLE(2*L + 1 )
                      FDNL1 = SSFDEL(K) * DNL1
                      FACT  = ONE - SSFDEL(K)
                      GK11 = (PHASMOMS_TOTAL_INPUT(L,K)-FDNL1)/FACT
                      HELP1 = L_PHASMOMS_TOTAL_INPUT(Q,L,K) * 
     &                          PHASMOMS_TOTAL_INPUT(L,K)
                      HELP2 = ( GK11 - DNL1 ) * L_SSFDEL(K,Q)
                      L_PHASMOM = ( HELP1 + HELP2 ) / FACT
                     ELSE
                      L_PHASMOM = L_PHASMOMS_TOTAL_INPUT(Q,L,K) *
     &                              PHASMOMS_TOTAL_INPUT(L,K)
                     ENDIF 
                     LEGPOLY = SS_PLEG_DN(V,1,L)
                     HELP = HELP + L_PHASMOM*LEGPOLY
                    ENDDO
                    L_EXACTSCAT_DN(V,K,Q) = HELP * TMS(K) + VAR1             
                   ELSE
                    L_EXACTSCAT_DN(V,K,Q) = VAR1
                   ENDIF
                  ENDDO
                 ENDIF
                ENDIF 
               ENDDO
              ENDIF

C  End Downwelling clause

            ENDIF

C   Finish geometry loops

          ENDDO
        ENDDO
      ENDDO

C  Recurrence relation for the UPWELLING intensity
C  ===============================================

      IF ( DO_UPWELLING ) THEN

C  initialize cumulative source term, and optical depth loop

        NC =  0
        DO V = 1, N_GEOMETRIES
          SS_CUMSOURCE_UP(V,NC) = ZERO
        ENDDO

C  initialise optical depth loop

        NSTART = NLAYERS
        NUT_PREV = NSTART + 1

C  Main loop over all output optical depths

        DO UTA = N_OUT_USERTAUS, 1, -1

C  Layer index for given optical depth

          NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
          NUT    = NLEVEL + 1

C  Cumulative single scatter source terms :
C      For loop over layers working upwards to level NUT,
C      Get layer source terms = Exact Z-matrix * Multiplier
C  Multiplier using new integration scheme

          DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            DO V = 1, N_GEOMETRIES
              HELP = EXACTSCAT_UP(V,N) * FLUX_FACTOR
              SS_LAYERSOURCE = HELP * UP_MULTIPLIERS(N,V)
              SS_CUMSOURCE_UP(V,NC) = SS_LAYERSOURCE +
     &               UP_LOSTRANS(N,V) * SS_CUMSOURCE_UP(V,NC-1)
            ENDDO
          ENDDO

C  Offgrid output-------
C    Add additional partial layer source term = Exact Phase Func * Multiplier
C    Get final cumulative source and set the Single scatter results
C  Ongrid output--------
C     Set final cumulative source and single scatter intensity

          IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN
            UT = OFFGRID_UTAU_OUTINDEX(UTA)
            N  = OFFGRID_UTAU_LAYERIDX(UT)
            DO V = 1, N_GEOMETRIES
              HELP           = EXACTSCAT_UP(V,N) * FLUX_FACTOR
              SS_LAYERSOURCE = HELP * UP_MULTIPLIERS_UT(UT,V)
              SS_CUMSOURCE   = SS_CUMSOURCE_UP(V,NC)
              TRANS          = UP_LOSTRANS_UT(UT,V)
              FINAL_SOURCE   = TRANS*SS_CUMSOURCE + SS_LAYERSOURCE
              SSCORRECTION   = SSFLUX * FINAL_SOURCE
              INTENSITY_SS(UTA,V,UPIDX) = SSCORRECTION
            ENDDO
          ELSE
            DO V = 1, N_GEOMETRIES
              FINAL_SOURCE = SS_CUMSOURCE_UP(V,NC)
              SSCORRECTION = SSFLUX * FINAL_SOURCE
              INTENSITY_SS(UTA,V,UPIDX) = SSCORRECTION
            ENDDO
          ENDIF

C  Check for updating the recursion 

          IF ( NUT. NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT

C  end optical depth loop and Upwelling clause

        ENDDO
      ENDIF

C  Recurrence relation for the DOWNWELLING intensity
C  =================================================

      IF ( DO_DNWELLING ) THEN

C  initialize cumulative source term, and optical depth loop

        NC =  0
        DO V = 1, N_GEOMETRIES
          SS_CUMSOURCE_DN(V,NC) = ZERO
        ENDDO

C  initialise optical depth loop

        NSTART = 1
        NUT_PREV = NSTART - 1

C  Main loop over all output optical depths

        DO UTA = 1, N_OUT_USERTAUS

C  Layer index for given optical depth

          NLEVEL = UTAU_LEVEL_MASK_DN(UTA)
          NUT = NLEVEL

C  Cumulative single scatter source terms :
C      For loop over layers working downwards to NUT,
C      Get layer source terms = Exact Z-matrix * Multiplier
C      Multiplier by new integration method

          DO N = NSTART, NUT
            NC = N
            DO V = 1, N_GEOMETRIES
              HELP =  EXACTSCAT_DN(V,N) * FLUX_FACTOR
              SS_LAYERSOURCE = HELP * DN_MULTIPLIERS(N,V)
              SS_CUMSOURCE_DN(V,NC) = SS_LAYERSOURCE +
     &                DN_LOSTRANS(N,V)*SS_CUMSOURCE_DN(V,NC-1)
            ENDDO
          ENDDO

C  Offgrid output :
C    add additional partial layer source term = Exact Z-matrix * Multiplier
C    Set final cumulative source and Correct the intensity
C  Ongrid output :
C     Set final cumulative source and correct Stokes vector

          IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN
            UT = OFFGRID_UTAU_OUTINDEX(UTA)
            N  = OFFGRID_UTAU_LAYERIDX(UT)
            DO V = 1, N_GEOMETRIES
              HELP = EXACTSCAT_DN(V,N) * FLUX_FACTOR
              SS_LAYERSOURCE = HELP * DN_MULTIPLIERS_UT(UT,V)
              SS_CUMSOURCE   = SS_CUMSOURCE_DN(V,NC)
              TRANS          = DN_LOSTRANS_UT(UT,V)
              FINAL_SOURCE   = TRANS*SS_CUMSOURCE + SS_LAYERSOURCE
              SSCORRECTION   = SSFLUX * FINAL_SOURCE
              INTENSITY_SS(UTA,V,DNIDX) = SSCORRECTION
            ENDDO
          ELSE
            DO V = 1, N_GEOMETRIES
              FINAL_SOURCE = SS_CUMSOURCE_DN(V,NC)
              SSCORRECTION = SSFLUX * FINAL_SOURCE
              INTENSITY_SS(UTA,V,DNIDX) = SSCORRECTION
            ENDDO
          ENDIF

C  Check for updating the recursion 

          IF ( NUT. NE. NUT_PREV ) NSTART = NUT + 1
          NUT_PREV = NUT

C  end optical depth loop and Downwelling intensity clause

        ENDDO
      ENDIF

C  Recurrence relation for the UPWELLING Jacobians
C  ===============================================

      IF ( DO_UPWELLING .AND. DO_ATMOS_LINEARIZATION ) THEN

C  Start the main layer variation loop
C  -----------------------------------
   
        DO K = 1, NTYPE_VARY
         IF ( DO_RTSOL_VARY(K) ) THEN
          K_PARAMETERS = NPARAMS_VARY(K) 
          KS = KINDEX(K)

C  initialize cumulative source term

          NC = 0
          DO V = 1, N_GEOMETRIES
           DO Q = 1, K_PARAMETERS
            L_SS_CUMSOURCE(Q,V) = ZERO
           ENDDO
          ENDDO

C  initialise optical depth loop

          NSTART = NLAYERS
          NUT_PREV = NSTART + 1

C  Main loop over all output optical depths
C  ----------------------------------------

          DO UTA = N_OUT_USERTAUS, 1, -1

C  Layer index for given optical depth

           NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
           NUT = NLEVEL + 1

C  Cumulative single scatter source terms to layer NUT
C  ---------------------------------------------------

C    1. Get layer source terms = Exact scattering * Multiplier
C    2. Loop over layers working upwards to NUT

            DO N = NSTART, NUT, -1
             NC = NLAYERS + 1 - N

C  If N = K (profile) or KS = 0 (column), there are extra linearizations.
C  If N > K, N < K, transmittance + some multiplication linearizations

             IF ( N.EQ.K .OR. KS.EQ.0 ) THEN
              DO V = 1, N_GEOMETRIES
               DO Q = 1, K_PARAMETERS
                LSS = EXACTSCAT_UP(V,N)   * L_UP_MULTIPLIERS(N,KS,Q,V)
     &            + L_EXACTSCAT_UP(V,N,Q) *   UP_MULTIPLIERS(N,V)
                L_SS_CUMSOURCE(Q,V) = FLUX_FACTOR * LSS +
     &            +  UP_LOSTRANS(N,V)    * L_SS_CUMSOURCE(Q,V)
     &            + L_UP_LOSTRANS(N,Q,V) *   SS_CUMSOURCE_UP(V,NC-1)
               ENDDO
              ENDDO
             ELSE
              DO V = 1, N_GEOMETRIES
               DO Q = 1, K_PARAMETERS
                LSS = EXACTSCAT_UP(V,N)   * L_UP_MULTIPLIERS(N,KS,Q,V)
                L_SS_CUMSOURCE(Q,V) = FLUX_FACTOR * LSS +
     &               +  UP_LOSTRANS(N,V) * L_SS_CUMSOURCE(Q,V)
               ENDDO
              ENDDO
             ENDIF

C  End layer loop

            ENDDO

C  Offgrid output----------
C  Set final cumulative source and Single scatter Weighting function

           IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN

            UT = OFFGRID_UTAU_OUTINDEX(UTA)
            N  = OFFGRID_UTAU_LAYERIDX(UT)

C  If N = K (the layer that is varying)
C    add Linearization of additional partial layer source term =
C        L_Exact_Scat(n) * Multiplier  +  Exact_Scat * L_Multiplier(n)
C  Variations when N > K
C    add Linearization of additional partial layer source term =
C         Exact_Scat * L_Multiplier(k)

C  If N = K (profile) or KS = 0 (column), there are extra linearizations.
C  If N > K, N < K, transmittance + some multiplication linearizations

            IF ( N.EQ.K .OR. KS.EQ.0 ) THEN
              DO V = 1, N_GEOMETRIES
                DO Q = 1, K_PARAMETERS
                 LSS = EXACTSCAT_UP(V,N)*L_UP_MULTIPLIERS_UT(UT,KS,Q,V)
     &           + L_EXACTSCAT_UP(V,N,Q)* UP_MULTIPLIERS_UT(UT,V)
                 L_FINAL_SOURCE = FLUX_FACTOR * LSS +
     &          +  UP_LOSTRANS_UT(UT,V)    * L_SS_CUMSOURCE(Q,V)
     &         + L_UP_LOSTRANS_UT(UT,Q,V) *   SS_CUMSOURCE_UP(V,NC)
                 L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                 IF ( DO_PROFILE_LINEARIZATION ) THEN
                   PROFILEWF_SS(Q,K,UTA,V,UPIDX) = L_SSCORRECTION
                 ELSE IF ( DO_COLUMN_LINEARIZATION ) THEN
                   COLUMNWF_SS(Q,UTA,V,UPIDX) = L_SSCORRECTION
                 ENDIF
                ENDDO
              ENDDO
            ELSE
              DO V = 1, N_GEOMETRIES
                DO Q = 1, K_PARAMETERS
                 LSS = EXACTSCAT_UP(V,N)*L_UP_MULTIPLIERS_UT(UT,KS,Q,V)
                 L_FINAL_SOURCE = FLUX_FACTOR * LSS +
     &          +  UP_LOSTRANS_UT(UT,V)    * L_SS_CUMSOURCE(Q,V)
                 L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                 IF ( DO_PROFILE_LINEARIZATION ) THEN
                   PROFILEWF_SS(Q,K,UTA,V,UPIDX) = L_SSCORRECTION
                 ELSE IF ( DO_COLUMN_LINEARIZATION ) THEN
                   COLUMNWF_SS(Q,UTA,V,UPIDX) = L_SSCORRECTION
                 ENDIF
                ENDDO
              ENDDO
            ENDIF

C  Ongrid output---------
C  just set to the cumulative source term 

           ELSE

            IF ( DO_PROFILE_LINEARIZATION ) THEN
             DO V = 1, N_GEOMETRIES
              DO Q = 1, K_PARAMETERS
                L_FINAL_SOURCE = L_SS_CUMSOURCE(Q,V)
                L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                PROFILEWF_SS(Q,K,UTA,V,UPIDX) = L_SSCORRECTION
              ENDDO
             ENDDO
            ELSE IF ( DO_COLUMN_LINEARIZATION ) THEN
             DO V = 1, N_GEOMETRIES
              DO Q = 1, K_PARAMETERS
                L_FINAL_SOURCE = L_SS_CUMSOURCE(Q,V)
                L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                COLUMNWF_SS(Q,UTA,V,UPIDX) = L_SSCORRECTION
              ENDDO
             ENDDO
            ENDIF

           ENDIF

C  Check for updating the recursion

           IF ( NUT. NE. NUT_PREV ) NSTART = NUT - 1
           NUT_PREV = NUT

C  end loop over output optical depths

          ENDDO

C  end loop over varying layers

         ENDIF
        ENDDO

C  end Upwelling Jacobian clause

      ENDIF

C  Recurrence relation for the DOWNWELLING Jacobians
C  =================================================

      IF ( DO_DNWELLING .AND. DO_ATMOS_LINEARIZATION ) THEN

C  Start the main layer variation loop
C  -----------------------------------

        DO K = 1, NTYPE_VARY
         IF ( DO_RTSOL_VARY(K) ) THEN
          K_PARAMETERS = NPARAMS_VARY(K) 
          KS = KINDEX(K)

C  initialize cumulative source term

          NC = 0
          DO V = 1, N_GEOMETRIES
           DO Q = 1, K_PARAMETERS
            L_SS_CUMSOURCE(Q,V) = ZERO
           ENDDO
          ENDDO

C  initialise optical depth loop

          NSTART = 1
          NUT_PREV = NSTART - 1

C  Main loop over all output optical depths

          DO UTA = 1, N_OUT_USERTAUS

C  Layer index for given optical depth

            NLEVEL = UTAU_LEVEL_MASK_DN(UTA)
            NUT = NLEVEL

C  Cumulative single scatter source terms to layer NUT
C  ---------------------------------------------------

C    1. Get layer source terms = Exact scattering * Multiplier
C    2. Loop over layers working downwards to NUT

            DO N = NSTART, NUT
             NC = N

C  If N = K (profile) or KS = 0 (column), there are extra linearizations.
C  If N > K, N < K, transmittance + some multiplication linearizations

             IF ( N.EQ.K .OR. KS.EQ.0 ) THEN
              DO V = 1, N_GEOMETRIES
               DO Q = 1, K_PARAMETERS
                LSS = EXACTSCAT_DN(V,N)   * L_DN_MULTIPLIERS(N,KS,Q,V)
     &            + L_EXACTSCAT_DN(V,N,Q) *   DN_MULTIPLIERS(N,V)
                L_SS_CUMSOURCE(Q,V) = FLUX_FACTOR * LSS +
     &            +  DN_LOSTRANS(N,V)    * L_SS_CUMSOURCE(Q,V)
     &            + L_DN_LOSTRANS(N,Q,V) *   SS_CUMSOURCE_DN(V,NC-1)
               ENDDO
              ENDDO
             ELSE
              DO V = 1, N_GEOMETRIES
               DO Q = 1, K_PARAMETERS
                LSS = EXACTSCAT_DN(V,N)   * L_DN_MULTIPLIERS(N,KS,Q,V)
                L_SS_CUMSOURCE(Q,V) = FLUX_FACTOR * LSS +
     &            +  DN_LOSTRANS(N,V) * L_SS_CUMSOURCE(Q,V)
               ENDDO
              ENDDO
             ENDIF

C  End layer loop

            ENDDO

C  Offgrid output----------
C  Set final cumulative source and Single scatter Weighting function

           IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN

            UT = OFFGRID_UTAU_OUTINDEX(UTA)
            N  = OFFGRID_UTAU_LAYERIDX(UT)

C  If N = K (the layer that is varying)
C    add Linearization of additional partial layer source term =
C        L_Exact_Scat(n) * Multiplier  +  Exact_Scat * L_Multiplier(n)
C  Variations when N > K
C    add Linearization of additional partial layer source term =
C         Exact_Scat * L_Multiplier(k)

C  If N = K (profile) or KS = 0 (column), there are extra linearizations.
C  If N > K, N < K, transmittance + some multiplication linearizations

            IF ( N.EQ.K .OR. KS.EQ.0 ) THEN
              DO V = 1, N_GEOMETRIES
                DO Q = 1, K_PARAMETERS
                 LSS = EXACTSCAT_Dn(V,N)*L_DN_MULTIPLIERS_UT(UT,KS,Q,V)
     &           + L_EXACTSCAT_DN(V,N,Q)* DN_MULTIPLIERS_UT(UT,V)
                 L_FINAL_SOURCE = FLUX_FACTOR * LSS +
     &          +  DN_LOSTRANS_UT(UT,V)    * L_SS_CUMSOURCE(Q,V)
     &         + L_DN_LOSTRANS_UT(UT,Q,V) *   SS_CUMSOURCE_DN(V,NC)
                 L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                 IF ( DO_PROFILE_LINEARIZATION ) THEN
                   PROFILEWF_SS(Q,K,UTA,V,DNIDX) = L_SSCORRECTION
                 ELSE IF ( DO_COLUMN_LINEARIZATION ) THEN
                   COLUMNWF_SS(Q,UTA,V,DNIDX) = L_SSCORRECTION
                 ENDIF
                ENDDO
              ENDDO
            ELSE
              DO V = 1, N_GEOMETRIES
                DO Q = 1, K_PARAMETERS
                 LSS = EXACTSCAT_DN(V,N)*L_DN_MULTIPLIERS_UT(UT,KS,Q,V)
                 L_FINAL_SOURCE = FLUX_FACTOR * LSS +
     &          +  DN_LOSTRANS_UT(UT,V)    * L_SS_CUMSOURCE(Q,V)
                 L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                 IF ( DO_PROFILE_LINEARIZATION ) THEN
                   PROFILEWF_SS(Q,K,UTA,V,DNIDX) = L_SSCORRECTION
                 ELSE IF ( DO_COLUMN_LINEARIZATION ) THEN
                   COLUMNWF_SS(Q,UTA,V,DNIDX) = L_SSCORRECTION
                 ENDIF
                ENDDO
              ENDDO
            ENDIF

C  Ongrid output---------
C  just set to the cumulative source term 

           ELSE

            IF ( DO_PROFILE_LINEARIZATION ) THEN
             DO V = 1, N_GEOMETRIES
              DO Q = 1, K_PARAMETERS
                L_FINAL_SOURCE = L_SS_CUMSOURCE(Q,V)
                L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                PROFILEWF_SS(Q,K,UTA,V,DNIDX) = L_SSCORRECTION
              ENDDO
             ENDDO
            ELSE IF ( DO_COLUMN_LINEARIZATION ) THEN
             DO V = 1, N_GEOMETRIES
              DO Q = 1, K_PARAMETERS
                L_FINAL_SOURCE = L_SS_CUMSOURCE(Q,V)
                L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                COLUMNWF_SS(Q,UTA,V,DNIDX) = L_SSCORRECTION
              ENDDO
             ENDDO
            ENDIF

           ENDIF

C  Check for updating the recursion

           IF ( NUT. NE. NUT_PREV ) NSTART = NUT + 1
           NUT_PREV = NUT

C  end loop over output optical depths

          ENDDO

C  end loop over varying layers

         ENDIF
        ENDDO

C  end Downwelling Jacobian clause

      ENDIF

C  Finish

      RETURN
      END

C

      subroutine l_outgoing_integration_up
     i   ( nlayers, nfinelayers,
     i     do_fine, do_partials,
     i     do_profile_linearization, layer_vary_flag, layer_vary_number,
     i     do_column_linearization, n_totalcolumn_wfs,
     i     extinction, l_extinction,
     i     n_partials, partials_idx, partials_fineidx,
     i     sunpaths, radii, ntraverse, alpha_all,
     i     sunpaths_p,    ntraverse_p,    alpha_p,
     i     sunpaths_fine, ntraverse_fine, alpha_fine,
     o       multipliers,     lostrans,     boa_attn,
     o       multipliers_p,   lostrans_p,
     o     l_multipliers,   l_lostrans,   l_boa_attn,
     o     l_multipliers_p, l_lostrans_p )

C  Does the optical depth integration over layers.
C  Partial layer integration added September 2007.

C  include files

      include '../includes/LIDORT.PARS'

C  control

      logical          do_fine, do_partials
      integer          nfinelayers, nlayers
      integer          n_partials
      integer          partials_idx    (max_offgrid_usertaus)
      integer          partials_fineidx(max_offgrid_usertaus)

C  Profile linearization control inputs

      logical          do_profile_linearization
      logical          layer_vary_flag   (maxlayers)
      integer          layer_vary_number (maxlayers)

C  Column linearization control inputs

      logical          do_column_linearization
      integer          n_totalcolumn_wfs

C  Whole layers

      integer          ntraverse(0:maxlayers)
      double precision sunpaths(0:maxlayers,maxlayers)
      double precision radii   (0:maxlayers)
      double precision alpha_all  (0:maxlayers)

C  Fine level

      integer          ntraverse_fine(maxlayers,maxfinelayers)
      double precision sunpaths_fine (maxlayers,maxlayers,maxfinelayers)
      double precision alpha_fine    (maxlayers,maxfinelayers)

C  Partial layers

      integer          ntraverse_p(max_offgrid_usertaus)
      double precision sunpaths_p (max_offgrid_usertaus,maxlayers)
      double precision alpha_p    (max_offgrid_usertaus)

C  Extinction inputs

      double precision extinction   (maxlayers)
      double precision L_extinction (maxlayers,max_atmoswfs)

C  outputs
C  -------

      double precision multipliers  (maxlayers)
      double precision lostrans     (maxlayers)
      double precision L_multipliers(maxlayers,0:maxlayers,max_atmoswfs)
      double precision L_lostrans   (maxlayers,max_atmoswfs)

      double precision multipliers_p   (max_offgrid_usertaus)
      double precision lostrans_p      (max_offgrid_usertaus)
      double precision L_multipliers_p 
     &         (max_offgrid_usertaus,0:maxlayers,max_atmoswfs)
      double precision L_lostrans_p
     &         (max_offgrid_usertaus,max_atmoswfs)

      double precision boa_attn
      double precision L_boa_attn(0:maxlayers,max_atmoswfs)

C  local arrays
C  ------------

C  Local geoemetry arrays

      double precision csq_fine ( maxfinelayers)
      double precision cot_fine ( maxfinelayers)

C  Local attenuation factors

      double precision attn ( 0:maxlayers )
      double precision attn_fine ( maxlayers, maxfinelayers )
      double precision attn_p    ( max_offgrid_usertaus )

C  Local attenuation factors, linearizations

      double precision l_attn ( 0:maxlayers, 0:maxlayers, max_atmoswfs )
      double precision l_attn_fine ( maxlayers, maxfinelayers,
     &                               0:maxlayers, max_atmoswfs )
      double precision l_attn_p
     &       ( max_offgrid_usertaus, 0:maxlayers, max_atmoswfs  )

C  help variables
C  --------------

      integer          n, j, k, k0, q, ut, np, nfine_p

      double precision tau, sum, salpha, calpha, dfine1, raycon, kn
      double precision csq_1, csq_2, cot_1, cot_2, argm_1, argj
      double precision func_1, func_2, tran_1, tran(maxfinelayers)
      double precision term1, term2, dfine_p, step_f, step_p
      double precision l_term1, l_term2

      double precision l_func_1, l_func_2, l_tran_1, l_tran
      double precision l_sum, l_func, l_kn, l_iop, step, func

C  local optical thickness cutoff
C      (should be same as MAX_TAU_SPATH in LIDORT)

      DOUBLE PRECISION LOCAL_CUTOFF
      PARAMETER       ( LOCAL_CUTOFF = 32.0D0 )

C  initialise output
C  -----------------

C  Whole layers

      do n = 1, nlayers
        multipliers(n) = 0.0d0
        lostrans(n)    = 0.0d0
      enddo

C  Partial layers
 
      if ( do_partials ) then
        do ut = 1, n_partials
          multipliers_p(ut) = 0.0d0
          lostrans_p(ut)    = 0.0d0
        enddo
      endif

C  Whole layer linearizations

      if ( do_profile_linearization ) then
        do n = 1, nlayers
          do k = 1, nlayers
            if ( layer_vary_flag(k) ) then
              do q = 1, layer_vary_number(k)
                L_multipliers(n,k,q) = zero
              enddo
            endif
          enddo
          if ( layer_vary_flag(n) ) then
            do q = 1, layer_vary_number(n)
              L_lostrans(n,q) = zero
            enddo
          endif
        enddo
      endif

      if ( do_column_linearization ) then
        k0 = 0
        do n = 1, nlayers
          do q = 1, n_totalcolumn_wfs
            L_multipliers(n,k0,q) = zero
            L_lostrans(n,q)       = zero
          enddo
        enddo
      endif

C  Partial layer linearizations

      if ( do_partials ) then
       if ( do_profile_linearization ) then
        do ut = 1, n_partials
          n = partials_idx(ut)
          do k = 1, nlayers
            if ( layer_vary_flag(k) ) then
              do q = 1, layer_vary_number(k)
                L_multipliers_p(ut,k,q) = zero
              enddo
            endif
          enddo
          if ( layer_vary_flag(n) ) then
            do q = 1, layer_vary_number(n)
              L_lostrans_p(ut,q) = zero
            enddo
          endif
        enddo
       endif
       if ( do_column_linearization ) then
        k0 = 0
        do ut = 1, n_partials
          do q = 1, n_totalcolumn_wfs
            L_multipliers_p(ut,k0,q) = zero
            L_lostrans_p(ut,q)       = zero
          enddo
        enddo
       endif
      endif

C  Create slant path optical attenuations, solar beams
C  ---------------------------------------------------

C  Ray constant at TOA

      salpha = dsin(alpha_all(0))
      raycon  = radii(0) * salpha
      dfine1 = dble(nfinelayers+1)

C  attenuation functions, whole layers

      do n = 0, nlayers
       tau     = 0.0d0
       attn(n) = 0.0d0
       do k = 1, ntraverse(n)
         tau = tau + sunpaths(n,k) * extinction(k)
       enddo
       if ( tau .le. local_cutoff ) attn(n) = dexp(-tau)
      enddo
      boa_attn = attn(nlayers)

C  Attenuations for partials

      if ( do_partials ) then
        do ut = 1, n_partials
          tau = 0.0d0
          attn_p(ut) = 0.0d0
          do k = 1, ntraverse_p(ut)
            tau = tau + sunpaths_p(ut,k) * extinction(k)
          enddo
          if ( tau.le.local_cutoff) attn_p(ut) = dexp(-tau)
        enddo
      endif

C  Attenuations for fine layer stuff

      do n = 1, nlayers
       do j = 1, nfinelayers
        tau            = 0.0d0
        attn_fine(n,j) = 0.0d0
        do k = 1, ntraverse_fine(n,j)
          tau = tau + sunpaths_fine(n,k,j) * extinction(k)
        enddo
        if ( tau .le. local_cutoff ) attn_fine(n,j) = dexp(-tau)
       enddo
      enddo

C  Linearized attenuation factors. Profile linearization.

      if ( do_profile_linearization ) then
        do n = 0, nlayers
          do k = 1, ntraverse(n)
            if ( layer_vary_flag(k) ) then
              do q = 1, layer_vary_number(k)
                l_iop = - l_extinction(k,q)
                l_attn(n,k,q) = sunpaths(n,k) * attn(n) * l_iop
                if (n.eq.nlayers)l_boa_attn(k,q)=l_attn(n,k,q)
              enddo
            endif
          enddo
        enddo
        do n = 1, nlayers
         do j = 1, nfinelayers
          do k = 1, ntraverse_fine(n,j)
            if ( layer_vary_flag(k) ) then
              do q = 1, layer_vary_number(k)
                l_iop = - l_extinction(k,q) * attn_fine(n,j)
                l_attn_fine(n,j,k,q) = sunpaths_fine(n,k,j) * l_iop
              enddo
            endif
          enddo
         enddo
        enddo
        if ( do_partials ) then
         do ut = 1, n_partials
          do k = 1, ntraverse_p(ut)
            if ( layer_vary_flag(k) ) then
              do q = 1, layer_vary_number(k)
                l_iop = - l_extinction(k,q) * attn_p(ut)
                l_attn_p(ut,k,q) = sunpaths_p(ut,k) * l_iop
              enddo
            endif
          enddo
         enddo
        endif
      endif

C  Linearized attenuation factors. Column linearization

      if ( do_column_linearization ) then
        do n = 0, nlayers
          do q = 1, n_totalcolumn_wfs
            sum = zero
            do k = 1, ntraverse(n)
              l_iop = - l_extinction(k,q)
              sum = sum + sunpaths(n,k) * l_iop
            enddo
            l_attn(n,k0,q) = sum * attn(n) 
            if (n.eq.nlayers)l_boa_attn(k0,q)=l_attn(n,k0,q)
          enddo
        enddo
        do n = 1, nlayers
         do j = 1, nfinelayers
          do q = 1, n_totalcolumn_wfs
            sum = zero
            do k = 1, ntraverse_fine(n,j)
              l_iop = - l_extinction(k,q)
              sum = sum + sunpaths_fine(n,k,j) * l_iop
            enddo
            l_attn_fine(n,j,k0,q) = sum * attn_fine(n,j)
          enddo
         enddo
        enddo
        if ( do_partials ) then
         do ut = 1, n_partials
          do q = 1, n_totalcolumn_wfs
            sum = zero
            do k = 1, ntraverse_p(ut)
              l_iop = - l_extinction(k,q)
              sum = sum + sunpaths_p(ut,k) * l_iop
            enddo
            l_attn_p(ut,k0,q) = sum * attn_p(ut)
          enddo
         enddo
        endif
      endif

C  Work up from the bottom of the atmosphere
C  =========================================

C  initialise

      n = nlayers
      salpha = dsin(alpha_all(n))
      calpha = dcos(alpha_all(n))
      csq_1 = 1.0d0 / salpha / salpha
      cot_1 = calpha / salpha

C  Start layer loop

      do n = nlayers, 1, -1

C  Save some quantities

        kn = raycon * extinction(n)
        salpha = dsin(alpha_all(n-1))
        calpha = dcos(alpha_all(n-1))
        csq_2 = 1.0d0 / salpha / salpha
        cot_2 = calpha / salpha
        step    = (alpha_all(n) - alpha_all(n-1))/dfine1

        do j = 1, nfinelayers
          calpha = dcos(alpha_fine(n,j))
          salpha = dsin(alpha_fine(n,j))
          cot_fine(j) = calpha / salpha
          csq_fine(j) = 1.0d0 / salpha / salpha
        enddo

C  integrated source term multiplier + transmittance
C   ---- Trapezium Rule Integration

C   Surely this is wrong.................???? 01 October 2007
c        func_1 = attn(n-1) * csq_1
c        argm_2 = cot_2 - cot_1 
c        tran_2 = dexp ( - kn * argm_2 )
c        func_2 = attn(n) * csq_2 * tran_2

        func_2 = attn(n-1)   * csq_2 
        argm_1 = cot_2 - cot_1 
        tran_1 = dexp ( - kn * argm_1 )
        func_1 = attn(n) * csq_1 * tran_1
        sum = 0.5d0 * ( func_1 + func_2 )
        do j = 1, nfinelayers
          tran(j) = dexp ( -kn * ( cot_2 - cot_fine(j) ) )
          func = attn_fine(n,j) * tran(j) * csq_fine(j)
          sum = sum + func
        enddo
        lostrans(n)    = tran_1
        multipliers(n) = sum * step * kn

C  Linearization.
C  Checked 30 January 2007

        if ( do_profile_linearization ) then
         do k = 1, nlayers
          if ( layer_vary_flag(k) ) then
           do q = 1, layer_vary_number(k)
            if ( k.eq.n) then
             l_kn = raycon * l_extinction(n,q)
             l_func_2 = l_attn(n-1,n,q) * csq_2
             l_tran_1 = -l_kn * argm_1 * tran_1
             l_func_1 = csq_1 * ( l_attn(n,n,q) *   tran_1 
     &                            + attn(n)       * l_tran_1 )
             l_sum = 0.5d0 * ( l_func_1 + l_func_2 )
             do j = 1, nfinelayers
              argj = cot_2 - cot_fine(j)
              l_tran = - l_kn * argj * tran(j)
              l_func = csq_fine(j) * ( l_attn_fine(n,j,n,q) * tran(j)
     &                               +   attn_fine(n,j)     * l_tran )
              l_sum = l_sum + l_func
             enddo
             l_lostrans(n,q)      = l_tran_1
             l_multipliers(n,n,q) = step * ( l_kn * sum + kn * l_sum )
            else
             l_func_2 = l_attn(n-1,k,q) * csq_2
             l_func_1 = csq_1 * l_attn(n,k,q) * tran_1 
             l_sum = 0.5d0 * ( l_func_1 + l_func_2 )
             do j = 1, nfinelayers
              l_func = csq_fine(j) * l_attn_fine(n,j,k,q) * tran(j)
              l_sum = l_sum + l_func
             enddo
             l_multipliers(n,k,q) = step * kn * l_sum
            endif
           enddo
          endif
         enddo
        endif

C  Column Linearization.

        if ( do_column_linearization ) then
          do q = 1, n_totalcolumn_wfs
            l_kn = raycon * l_extinction(n,q)
            l_func_2 = l_attn(n-1,k0,q) * csq_2
            l_tran_1 = -l_kn * argm_1 * tran_1
            l_func_1 = csq_1 * ( l_attn(n,k0,q) *   tran_1 
     &                           + attn(n)       * l_tran_1 )
            l_sum = 0.5d0 * ( l_func_1 + l_func_2 )
            do j = 1, nfinelayers
              argj = cot_2 - cot_fine(j)
              l_tran = - l_kn * argj * tran(j)
              l_func = csq_fine(j) * ( l_attn_fine(n,j,k0,q) * tran(j)
     &                               +   attn_fine(n,j)      * l_tran )
              l_sum = l_sum + l_func
            enddo
            l_lostrans(n,q)      = l_tran_1
            l_multipliers(n,k0,q) = step * ( l_kn * sum + kn * l_sum )
          enddo
        endif

C  update the geometry

        cot_1 = cot_2
        csq_1 = csq_2

C  Finish layer

      ENDDO

C  Finish if no partials

      if ( .not. do_partials ) return
      
C  Partial layer functions
C  =======================

C  start partials loop

      do ut = 1, n_partials
      
C  layer of occurrence and fine-layer cutoff

        np     = partials_idx(ut)
        nfine_p = partials_fineidx(ut)
        if ( nfine_p .gt. 0 ) then
          dfine_p = dble(nfine_p)
          step_f = (alpha_all(np) - alpha_fine(np,nfine_p))/dfine_p
          step_p = (alpha_fine(np,nfine_p)-alpha_p(ut))
        else
          step_p = (alpha_all(np)-alpha_p(ut))
        endif
        
C  Bottom of layer of occurrence

        salpha = dsin(alpha_all(np))
        calpha = dcos(alpha_all(np))
        csq_1 = 1.0d0 / salpha / salpha
        cot_1 = calpha / salpha
        kn = raycon * extinction(np)

C  Top of integration = off-grid value

        salpha = dsin(alpha_p(ut))
        calpha = dcos(alpha_p(ut))
        csq_2 = 1.0d0 / salpha / salpha
        cot_2 = calpha / salpha

C  saved fine-grid quantities as far as cutoff

        do j = 1, nfine_p
          calpha = dcos(alpha_fine(np,j))
          salpha = dsin(alpha_fine(np,j))
          cot_fine(j) = calpha / salpha
          csq_fine(j) = 1.0d0 / salpha / salpha
        enddo

C  integrated source term multiplier + transmittance
C   ---- Trapezium Rule Integration
C  Careful with the last step of the integration 

        func_2 = attn_p(ut)   * csq_2 
        argm_1 = cot_2 - cot_1 
        tran_1 = dexp ( - kn * argm_1  )
        func_1 = attn(np) * csq_1 * tran_1

        if ( nfine_p.gt.0) then
          sum = 0.5d0 * func_1
          do j = 1, nfine_p
            tran(j) = dexp ( -kn * ( cot_2 - cot_fine(j) ) )
            func = attn_fine(np,j) * tran(j) * csq_fine(j)
            if ( j.lt.nfine_p ) sum = sum + func
            if ( j.eq.nfine_p ) sum = sum + 0.5d0*func
          enddo
          term1 = sum * step_f 
          term2 = (func + func_2) * 0.5d0 * step_p
        else
          term1 = 0.5d0 * func_1 * step_p
          term2 = 0.5d0 * func_2 * step_p
        endif
        lostrans_p(ut)    = tran_1
        multipliers_p(ut) = ( term1 + term2 ) * kn

C  Linearization.
C    Coded up 26 September 2007. 

        if ( do_profile_linearization ) then

         do k = 1, nlayers
          if ( layer_vary_flag(k) ) then
           do q = 1, layer_vary_number(k)
            if ( k.eq.np) then
             l_kn = raycon * l_extinction(np,q)
             l_func_2 = l_attn_p(ut,np,q) * csq_2
             l_tran_1 = -l_kn * argm_1 * tran_1
             l_func_1 = csq_1 * ( l_attn(np,np,q) *   tran_1 
     &                            + attn(np)      * l_tran_1 )
             if ( nfine_p.gt.0) then
              l_sum = 0.5d0 * l_func_1
              do j = 1, nfine_p
                argj = cot_2 - cot_fine(j)
                l_tran = - l_kn * argj * tran(j)
                l_func = csq_fine(j) * ( l_attn_fine(np,j,np,q)* tran(j)
     &                               +   attn_fine(np,j)     * l_tran )
                if ( j.lt.nfine_p ) l_sum = l_sum + l_func
                if ( j.eq.nfine_p ) l_sum = l_sum + 0.5d0*l_func
              enddo
              l_term1 = l_sum * step_f 
              l_term2 = ( l_func + l_func_2 ) * 0.5d0 * step_p
             else
              l_term1 = 0.5d0 * l_func_1 * step_p
              l_term2 = 0.5d0 * l_func_2 * step_p
             endif
             l_lostrans_p(ut,q)       = l_tran_1
             l_multipliers_p(ut,np,q) =  ( term1 +   term2 ) * l_kn +
     &                                 ( l_term1 + l_term2 ) * kn
            else
             l_func_2 = l_attn_p(ut,k,q) * csq_2
             l_func_1 = csq_1 * l_attn(np,k,q) * tran_1 
             if ( nfine_p.gt.0) then
              l_sum = 0.5d0 * l_func_1
              do j = 1, nfine_p
               l_func = csq_fine(j) * l_attn_fine(np,j,k,q) * tran(j)
               if ( j.lt.nfine_p ) l_sum = l_sum + l_func
               if ( j.eq.nfine_p ) l_sum = l_sum + 0.5d0*l_func
              enddo
              l_term1 = l_sum * step_f 
              l_term2 = ( l_func + l_func_2 ) * 0.5d0 * step_p
             else
              l_term1 = 0.5d0 * l_func_1 * step_p
              l_term2 = 0.5d0 * l_func_2 * step_p
             endif
             l_multipliers_p(ut,k,q) =  ( l_term1 + l_term2 ) * kn
            endif
           enddo
          endif
         enddo

        endif

C  Column linearization

        if ( do_column_linearization ) then
          do q = 1, n_totalcolumn_wfs
            l_kn = raycon * l_extinction(np,q)
            l_func_2 = l_attn_p(ut,k0,q) * csq_2
            l_tran_1 = -l_kn * argm_1 * tran_1
            l_func_1 = csq_1 * ( l_attn(np,k0,q) *   tran_1 
     &                           + attn(np)      * l_tran_1 )
            if ( nfine_p.gt.0) then
              l_sum = 0.5d0 * l_func_1
              do j = 1, nfine_p
                argj = cot_2 - cot_fine(j)
                l_tran = - l_kn * argj * tran(j)
                l_func = csq_fine(j) * ( l_attn_fine(np,j,k0,q)* tran(j)
     &                               +   attn_fine(np,j)     * l_tran )
                if ( j.lt.nfine_p ) l_sum = l_sum + l_func
                if ( j.eq.nfine_p ) l_sum = l_sum + 0.5d0*l_func
              enddo
              l_term1 = l_sum * step_f 
              l_term2 = ( l_func + l_func_2 ) * 0.5d0 * step_p
            else
              l_term1 = 0.5d0 * l_func_1 * step_p
              l_term2 = 0.5d0 * l_func_2 * step_p
            endif
            l_lostrans_p(ut,q)       = l_tran_1
            l_multipliers_p(ut,k0,q) =  ( term1 +   term2 ) * l_kn +
     &                                 ( l_term1 + l_term2 ) * kn
          enddo
        endif

C  Finish partials loop

      ENDDO

C  finish

      RETURN
      END

c

      subroutine l_outgoing_integration_dn
     i   ( nlayers, nfinelayers,
     i     do_fine, do_partials,
     i     do_profile_linearization, layer_vary_flag, layer_vary_number,
     i     do_column_linearization, n_totalcolumn_wfs,
     i     extinction, l_extinction,
     i     n_partials, partials_idx, partials_fineidx,
     i     sunpaths, radii, ntraverse, alpha_all,
     i     sunpaths_p,    ntraverse_p,    alpha_p,
     i     sunpaths_fine, ntraverse_fine, alpha_fine,
     o       multipliers,     lostrans,
     o       multipliers_p,   lostrans_p,
     o     l_multipliers,   l_lostrans,
     o     l_multipliers_p, l_lostrans_p )

C  Does the optical depth integration over layers.
C  Partial layer integration added September 2007.

C  include files

      include '../includes/LIDORT.PARS'

C  control

      logical          do_fine, do_partials
      integer          nfinelayers, nlayers
      integer          n_partials
      integer          partials_idx    (max_offgrid_usertaus)
      integer          partials_fineidx(max_offgrid_usertaus)

C  Profile linearization control inputs

      logical          do_profile_linearization
      logical          layer_vary_flag   (maxlayers)
      integer          layer_vary_number (maxlayers)

C  Column linearization control inputs

      logical          do_column_linearization
      integer          n_totalcolumn_wfs

C  Whole layers

      integer          ntraverse(0:maxlayers)
      double precision sunpaths(0:maxlayers,maxlayers)
      double precision radii   (0:maxlayers)
      double precision alpha_all  (0:maxlayers)

C  Fine level

      integer          ntraverse_fine(maxlayers,maxfinelayers)
      double precision sunpaths_fine (maxlayers,maxlayers,maxfinelayers)
      double precision alpha_fine    (maxlayers,maxfinelayers)

C  Partial layers

      integer          ntraverse_p(max_offgrid_usertaus)
      double precision sunpaths_p (max_offgrid_usertaus,maxlayers)
      double precision alpha_p    (max_offgrid_usertaus)

C  Extinction inputs

      double precision extinction   (maxlayers)
      double precision L_extinction (maxlayers,max_atmoswfs)

C  outputs
C  -------

      double precision multipliers  (maxlayers)
      double precision lostrans     (maxlayers)
      double precision L_multipliers(maxlayers,0:maxlayers,max_atmoswfs)
      double precision L_lostrans   (maxlayers,max_atmoswfs)

      double precision multipliers_p   (max_offgrid_usertaus)
      double precision lostrans_p      (max_offgrid_usertaus)
      double precision L_multipliers_p 
     &         (max_offgrid_usertaus,0:maxlayers,max_atmoswfs)
      double precision L_lostrans_p
     &         (max_offgrid_usertaus,max_atmoswfs)

C  local arrays
C  ------------

C  Local geoemetry arrays

      double precision csq_fine ( maxfinelayers)
      double precision cot_fine ( maxfinelayers)

C  Local attenuation factors

      double precision attn ( 0:maxlayers )
      double precision attn_fine ( maxlayers, maxfinelayers )
      double precision attn_p    ( max_offgrid_usertaus )

C  Local attenuation factors, linearizations

      double precision l_attn ( 0:maxlayers, 0:maxlayers, max_atmoswfs )
      double precision l_attn_fine ( maxlayers, maxfinelayers,
     &                               0:maxlayers, max_atmoswfs )
      double precision l_attn_p
     &       ( max_offgrid_usertaus, 0:maxlayers, max_atmoswfs  )

C  help variables
C  --------------

      integer          n, j, k, k0, q, ut, np, nfine_p

      double precision tau, sum, salpha, calpha, dfine1, raycon, kn
      double precision csq_1, csq_2, cot_1, cot_2, argm_1, argj
      double precision func_1, func_2, tran_1, tran(maxfinelayers)
      double precision term1, term2, dfine_p, step_f, step_p
      double precision l_term1, l_term2

      double precision l_func_1, l_func_2, l_tran_1, l_tran
      double precision l_sum, l_func, l_kn, l_iop, step, func

C  local optical thickness cutoff
C      (should be same as MAX_TAU_SPATH in LIDORT)

      DOUBLE PRECISION LOCAL_CUTOFF
      PARAMETER       ( LOCAL_CUTOFF = 32.0D0 )

C  initialise output
C  -----------------

C  Whole layers

      do n = 1, nlayers
        multipliers(n) = 0.0d0
        lostrans(n)    = 0.0d0
      enddo

C  Partial layers
 
      if ( do_partials ) then
        do ut = 1, n_partials
          multipliers_p(ut) = 0.0d0
          lostrans_p(ut)    = 0.0d0
        enddo
      endif

C  profile layer linearizations

      if ( do_profile_linearization ) then
        do n = 1, nlayers
          do k = 1, nlayers
            if ( layer_vary_flag(k) ) then
              do q = 1, layer_vary_number(k)
                L_multipliers(n,k,q) = zero
              enddo
            endif
          enddo
          if ( layer_vary_flag(n) ) then
            do q = 1, layer_vary_number(n)
              L_lostrans(n,q) = zero
            enddo
          endif
        enddo
      endif

      if ( do_column_linearization ) then
        k0 = 0
        do n = 1, nlayers
          do q = 1, n_totalcolumn_wfs
            L_multipliers(n,k0,q) = zero
            L_lostrans(n,q)       = zero
          enddo
        enddo
      endif

C  Partial layer linearizations

      if ( do_partials ) then
       if ( do_profile_linearization ) then
        do ut = 1, n_partials
          n = partials_idx(ut)
          do k = 1, nlayers
            if ( layer_vary_flag(k) ) then
              do q = 1, layer_vary_number(k)
                L_multipliers_p(ut,k,q) = zero
              enddo
            endif
          enddo
          if ( layer_vary_flag(n) ) then
            do q = 1, layer_vary_number(n)
              L_lostrans_p(ut,q) = zero
            enddo
          endif
        enddo
       endif
       if ( do_column_linearization ) then
        k0 = 0
        do ut = 1, n_partials
          do q = 1, n_totalcolumn_wfs
            L_multipliers_p(ut,k0,q) = zero
            L_lostrans_p(ut,q)       = zero
          enddo
        enddo
       endif
      endif

C  Create slant path optical attenuations, solar beams
C  ---------------------------------------------------

C  Ray constant at TOA

      salpha = dsin(alpha_all(0))
      raycon  = radii(0) * salpha
      dfine1 = dble(nfinelayers+1)

C  attenuation functions, whole layers

      do n = 0, nlayers
       tau     = 0.0d0
       attn(n) = 0.0d0
       do k = 1, ntraverse(n)
         tau = tau + sunpaths(n,k) * extinction(k)
       enddo
       if ( tau .le. local_cutoff ) attn(n) = dexp(-tau)
      enddo

C  Attenuations for partials

      if ( do_partials ) then
        do ut = 1, n_partials
          tau = 0.0d0
          attn_p(ut) = 0.0d0
          do k = 1, ntraverse_p(ut)
            tau = tau + sunpaths_p(ut,k) * extinction(k)
          enddo
          if ( tau.le.local_cutoff) attn_p(ut) = dexp(-tau)
        enddo
      endif

C  Attenuations for fine layer stuff

      do n = 1, nlayers
       do j = 1, nfinelayers
        tau            = 0.0d0
        attn_fine(n,j) = 0.0d0
        do k = 1, ntraverse_fine(n,j)
          tau = tau + sunpaths_fine(n,k,j) * extinction(k)
        enddo
        if ( tau .le. local_cutoff ) attn_fine(n,j) = dexp(-tau)
       enddo
      enddo

C  Linearized attenuation factors. Profile linearization

      if ( do_profile_linearization ) then
        do n = 0, nlayers
          do k = 1, ntraverse(n)
            if ( layer_vary_flag(k) ) then
              do q = 1, layer_vary_number(k)
                l_iop = - l_extinction(k,q)
                l_attn(n,k,q) = sunpaths(n,k) * attn(n) * l_iop
              enddo
            endif
          enddo
        enddo
        do n = 1, nlayers
         do j = 1, nfinelayers
          do k = 1, ntraverse_fine(n,j)
            if ( layer_vary_flag(k) ) then
              do q = 1, layer_vary_number(k)
                l_iop = - l_extinction(k,q) * attn_fine(n,j)
                l_attn_fine(n,j,k,q) = sunpaths_fine(n,k,j) * l_iop
              enddo
            endif
          enddo
         enddo
        enddo
        if ( do_partials ) then
         do ut = 1, n_partials
          do k = 1, ntraverse_p(ut)
            if ( layer_vary_flag(k) ) then
              do q = 1, layer_vary_number(k)
                l_iop = - l_extinction(k,q) * attn_p(ut)
                l_attn_p(ut,k,q) = sunpaths_p(ut,k) * l_iop
              enddo
            endif
          enddo
         enddo
        endif
      endif

C  Linearized attenuation factors. Column linearization

      if ( do_column_linearization ) then
        do n = 0, nlayers
          do q = 1, n_totalcolumn_wfs
            sum = zero
            do k = 1, ntraverse(n)
              l_iop = - l_extinction(k,q)
              sum = sum + sunpaths(n,k) * l_iop
            enddo
            l_attn(n,k0,q) = sum * attn(n) 
          enddo
        enddo
        do n = 1, nlayers
         do j = 1, nfinelayers
          do q = 1, n_totalcolumn_wfs
            sum = zero
            do k = 1, ntraverse_fine(n,j)
              l_iop = - l_extinction(k,q)
              sum = sum + sunpaths_fine(n,k,j) * l_iop
            enddo
            l_attn_fine(n,j,k0,q) = sum * attn_fine(n,j)
          enddo
         enddo
        enddo
        if ( do_partials ) then
         do ut = 1, n_partials
          do q = 1, n_totalcolumn_wfs
            sum = zero
            do k = 1, ntraverse_p(ut)
              l_iop = - l_extinction(k,q)
              sum = sum + sunpaths_p(ut,k) * l_iop
            enddo
            l_attn_p(ut,k0,q) = sum * attn_p(ut)
          enddo
         enddo
        endif
      endif

C  Work down from the top of the atmosphere
C  ========================================

C  initialise

      n = 0
      salpha = dsin(alpha_all(n))
      calpha = dcos(alpha_all(n))
      csq_1 = 1.0d0 / salpha / salpha
      cot_1 = calpha / salpha

C  Start layer loop

      do n = 1, nlayers

C  Save some quantities

        kn = raycon * extinction(n)
        salpha = dsin(alpha_all(n))
        calpha = dcos(alpha_all(n))
        csq_2 = 1.0d0 / salpha / salpha
        cot_2 = calpha / salpha
        step    = (alpha_all(n) - alpha_all(n-1))/dfine1
        do j = 1, nfinelayers
          calpha = dcos(alpha_fine(n,j))
          salpha = dsin(alpha_fine(n,j))
          cot_fine(j) = calpha / salpha
          csq_fine(j) = 1.0d0 / salpha / salpha
        enddo

C  integrated source term multiplier + transmittance
C   ---- Trapezium Rule Integration

        argm_1 = cot_1 - cot_2 
        tran_1 = dexp ( - kn * argm_1 )
        func_1 = attn(n-1) * csq_1 * tran_1
        func_2 = attn(n) * csq_2
        sum = 0.5d0 * ( func_1 + func_2 )
        do j = nfinelayers, 1, -1
          tran(j) = dexp ( -kn * ( cot_fine(j) - cot_2 ) )
          func = attn_fine(n,j) * tran(j) * csq_fine(j)
          sum = sum + func
        enddo
        lostrans(n)    = tran_1
        multipliers(n) = sum * step * kn

C  Profile Linearization.
C  Checked 30 January 2007

        if ( do_profile_linearization ) then
         do k = 1, nlayers
          if ( layer_vary_flag(k) ) then
           do q = 1, layer_vary_number(k)
            if ( k.eq.n) then
             l_kn = raycon * l_extinction(n,q)
             l_func_2 = l_attn(n,n,q) * csq_2
             l_tran_1 = -l_kn * argm_1 * tran_1
             l_func_1 = csq_1 * ( l_attn(n-1,n,q) *   tran_1 
     &                            + attn(n-1)     * l_tran_1 )
             l_sum = 0.5d0 * ( l_func_1 + l_func_2 )
             do j = nfinelayers, 1, -1
              argj = cot_fine(j) - cot_2
              l_tran = - l_kn * argj * tran(j)
              l_func = csq_fine(j) * ( l_attn_fine(n,j,n,q) * tran(j)
     &                               +   attn_fine(n,j)     * l_tran )
              l_sum = l_sum + l_func
             enddo
             l_lostrans(n,q)      = l_tran_1
             l_multipliers(n,n,q) = step * ( l_kn * sum + kn * l_sum )
            else
             l_func_1 = l_attn(n-1,k,q) * csq_1 * tran_1
             l_func_2 = csq_2 * l_attn(n,k,q)
             l_sum = 0.5d0 * ( l_func_1 + l_func_2 )
             do j = 1, nfinelayers
              l_func = csq_fine(j) * l_attn_fine(n,j,k,q) * tran(j)
              l_sum = l_sum + l_func
             enddo
             l_multipliers(n,k,q) = step * kn * l_sum
            endif
           enddo
          endif
         enddo
        endif

C  Column linearization

        if ( do_column_linearization ) then
          do q = 1, n_totalcolumn_wfs
            l_kn = raycon * l_extinction(n,q)
            l_func_2 = l_attn(n,k0,q) * csq_2
            l_tran_1 = -l_kn * argm_1 * tran_1
            l_func_1 = csq_1 * ( l_attn(n-1,k0,q) *   tran_1 
     &                           + attn(n-1)     * l_tran_1 )
            l_sum = 0.5d0 * ( l_func_1 + l_func_2 )
            do j = nfinelayers, 1, -1
              argj = cot_fine(j) - cot_2
              l_tran = - l_kn * argj * tran(j)
              l_func = csq_fine(j) * ( l_attn_fine(n,j,k0,q) * tran(j)
     &                               +   attn_fine(n,j)      * l_tran )
              l_sum = l_sum + l_func
            enddo
            l_lostrans(n,q)       = l_tran_1
            l_multipliers(n,k0,q) = step * ( l_kn * sum + kn * l_sum )
          enddo
        endif

C  update the geometry

        cot_1 = cot_2
        csq_1 = csq_2

C  Finish layer

      ENDDO

C  Finish if no partials

      if ( .not. do_partials ) return
      
C  Partial layer functions
C  =======================

C  start partials loop

      do ut = 1, n_partials
      
C  layer of occurrence and fine-layer cutoff

        np     = partials_idx(ut)
        nfine_p = partials_fineidx(ut)
        if ( nfine_p .eq. nfinelayers ) then
          step_p = (alpha_p(ut)-alpha_all(np-1))
        else if ( nfine_p.eq.0 ) then
          step_f = (alpha_all(np) - alpha_all(np-1))/dfine1
          step_p = (alpha_p(ut)-alpha_fine(np,nfine_p+1))
        else
          dfine_p = dble(nfine_p)
          step_f = (alpha_all(np) - alpha_fine(np,nfine_p))/dfine_p
          step_p = (alpha_p(ut)-alpha_fine(np,nfine_p+1))
        endif
        
C  Top of layer of occurrence

        salpha = dsin(alpha_all(np-1))
        calpha = dcos(alpha_all(np-1))
        csq_1 = 1.0d0 / salpha / salpha
        cot_1 = calpha / salpha
        kn = raycon * extinction(np)

C  End of integration = off-grid value

        salpha = dsin(alpha_p(ut))
        calpha = dcos(alpha_p(ut))
        csq_2 = 1.0d0 / salpha / salpha
        cot_2 = calpha / salpha

C  saved fine-grid quantities as far as cutoff

        do j = nfinelayers, nfine_p+1, -1
          calpha = dcos(alpha_fine(np,j))
          salpha = dsin(alpha_fine(np,j))
          cot_fine(j) = calpha / salpha
          csq_fine(j) = 1.0d0 / salpha / salpha
        enddo

C  integrated source term multiplier + transmittance
C   ---- Trapezium Rule Integration
C  Careful with the last step of the integration 

        func_2 = attn_p(ut)   * csq_2
        argm_1 = cot_1 - cot_2
        tran_1 = dexp ( - kn * argm_1 )
        func_1 = attn(np-1) * csq_1 * tran_1
        if ( nfine_p.lt.nfinelayers) then
          sum = 0.5d0 * func_1
          do j = nfinelayers, nfine_p+1, -1
            tran(j) = dexp ( -kn * ( cot_fine(j) - cot_2 ) )
            func  = attn_fine(np,j) * tran(j) * csq_fine(j)
            if ( j.gt.nfine_p+1 ) sum = sum + func
            if ( j.eq.nfine_p+1 ) sum = sum + 0.5d0*func
          enddo
          term1 = sum * step_f 
          term2 = ( func + func_2) * 0.5d0 * step_p
        else
          term1 = 0.5d0 * func_1 * step_p
          term2 = 0.5d0 * func_2 * step_p
        endif
        lostrans_p(ut)    = tran_1
        multipliers_p(ut) = ( term1 + term2) * kn

C  Profile Linearization.
C    Coded up 26 September 2007.

        if ( do_profile_linearization ) then
         do k = 1, nlayers
          if ( layer_vary_flag(k) ) then
           do q = 1, layer_vary_number(k)
            if ( k.eq.np) then
             l_kn = raycon * l_extinction(np,q)
             l_func_2 = l_attn_p(ut,np,q) * csq_2
             l_tran_1 = -l_kn * argm_1 * tran_1
             l_func_1 = csq_1 * ( l_attn(np-1,np,q) *   tran_1 
     &                            + attn(np-1)      * l_tran_1 )
             if ( nfine_p.lt.nfinelayers) then
              l_sum = 0.5d0 * l_func_1
              do j = nfinelayers, nfine_p+1, -1
               argj = cot_fine(j) - cot_2 
               l_tran = - l_kn * argj * tran(j)
               l_func = csq_fine(j) * ( l_attn_fine(np,j,np,q) * tran(j)
     &                               +   attn_fine(np,j)     * l_tran )
               if ( j.gt.nfine_p+1 ) l_sum = l_sum + l_func
               if ( j.eq.nfine_p+1 ) l_sum = l_sum + 0.5d0*l_func
              enddo
              l_term1 = l_sum * step_f 
              l_term2 = ( l_func + l_func_2 ) * 0.5d0 * step_p
             else
              l_term1 = 0.5d0 * l_func_1 * step_p
              l_term2 = 0.5d0 * l_func_2 * step_p
             endif
             l_lostrans_p(ut,q)       = l_tran_1
             l_multipliers_p(ut,np,q) =  ( term1 +   term2 ) * l_kn +
     &                                 ( l_term1 + l_term2 ) * kn
            else
             l_func_2 = l_attn_p(ut,k,q) * csq_2
             l_func_1 = csq_1 * l_attn(np-1,k,q) *   tran_1 
             if ( nfine_p.lt.nfinelayers) then
              l_sum = 0.5d0 * l_func_1
              do j = nfinelayers, nfine_p+1, -1
               l_func = csq_fine(j) * l_attn_fine(np,j,k,q) * tran(j)
               if ( j.gt.nfine_p+1 ) l_sum = l_sum + l_func
               if ( j.eq.nfine_p+1 ) l_sum = l_sum + 0.5d0*l_func
              enddo
              l_term1 = l_sum * step_f 
              l_term2 = ( l_func + l_func_2 ) * 0.5d0 * step_p
             else
              l_term1 = 0.5d0 * l_func_1 * step_p
              l_term2 = 0.5d0 * l_func_2 * step_p
             endif
             l_multipliers_p(ut,k,q) = ( l_term1 + l_term2 ) * kn
            endif
           enddo
          endif
         enddo
        endif

C  Column Linearization.

        if ( do_column_linearization ) then
          do q = 1, n_totalcolumn_wfs
            l_kn = raycon * l_extinction(np,q)
            l_func_2 = l_attn_p(ut,k0,q) * csq_2
            l_tran_1 = -l_kn * argm_1 * tran_1
            l_func_1 = csq_1 * ( l_attn(np-1,k0,q) *   tran_1 
     &                           + attn(np-1)      * l_tran_1 )
            if ( nfine_p.lt.nfinelayers) then
              l_sum = 0.5d0 * l_func_1
              do j = nfinelayers, nfine_p+1, -1
               argj = cot_fine(j) - cot_2 
               l_tran = - l_kn * argj * tran(j)
               l_func = csq_fine(j) * ( l_attn_fine(np,j,k0,q) * tran(j)
     &                               +   attn_fine(np,j)     * l_tran )
               if ( j.gt.nfine_p+1 ) l_sum = l_sum + l_func
               if ( j.eq.nfine_p+1 ) l_sum = l_sum + 0.5d0*l_func
              enddo
              l_term1 = l_sum * step_f 
              l_term2 = ( l_func + l_func_2 ) * 0.5d0 * step_p
            else
              l_term1 = 0.5d0 * l_func_1 * step_p
              l_term2 = 0.5d0 * l_func_2 * step_p
            endif
            l_lostrans_p(ut,q)       = l_tran_1
            l_multipliers_p(ut,k0,q) =  ( term1 +   term2 ) * l_kn +
     &                                ( l_term1 + l_term2 ) * kn
          enddo
        endif

C  Finish partials loop

      ENDDO

C  finish

      RETURN
      END
      
c

      SUBROUTINE LIDORT_LAP_LAMBERTIAN_DBCORRECTION(FLUXMULT)

C  Prepares Linearization of Exact Direct Beam reflection (Lambertian case)
C   Linearization with respect to atmospheric variables

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  Include file of setup variables (Input to the present module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  Include files of Linearized input and setup variables

      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'

C  Include files of Single scatter & surface reflectance variables (Input)

      INCLUDE '../includes/LIDORT_SINGSCAT.VARS'
      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  Direct beam output stored in linearized correction file

      INCLUDE '../includes/LIDORT_L_SINGSCAT.VARS'

C  input argument

      DOUBLE PRECISION FLUXMULT

C  Local variables
C  ---------------

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER          UT, UTA, UM, UA, NC, IB, V, K, Q
      DOUBLE PRECISION X0_FLUX, X0_BOA, FACTOR, L_ATTN, TR, LTR

C  first stage
C  -----------

C  return if no upwelling

      IF ( .NOT.DO_UPWELLING ) RETURN

C  Lambertian factor

      FACTOR = LAMBERTIAN_ALBEDO * FLUXMULT

C  Start weighting function loop

      DO K = 1, NLAYERS
        IF ( LAYER_VARY_FLAG(K) ) THEN
          DO Q = 1, LAYER_VARY_NUMBER(K)

C  Initialize the output results

            DO V = 1, N_GEOMETRIES
              DO UTA = 1, N_OUT_USERTAUS
                PROFILEWF_DB(Q,K,UTA,V) = ZERO
              ENDDO
            ENDDO

C  Skip calculation if albedo zero

            IF ( FACTOR .EQ. ZERO ) GO TO 6543

C  New Code  R. Spurr, 6 August 2007. RT Solutions Inc.
C  ====================================================

C  initialize cumulative source term

C  Start geometry loops

            DO UM = 1, N_USER_STREAMS
             DO IB = 1, NBEAMS
              IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
               DO UA = 1, N_USER_RELAZMS
                V = UMOFF(IB,UM) + UA

c  Beam attenuation and reflection

                IF ( DO_SSCORR_OUTGOING ) THEN
                  X0_BOA = DCOS(BEAM_SZAS_ADJUST(UM,IB,UA)*DEG_TO_RAD)
                  L_ATTN = L_BOA_ATTN(K,Q,V)
                ELSE
                 IF ( DO_REFRACTIVE_GEOMETRY ) THEN
                  X0_BOA = DCOS(SZA_LOCAL_INPUT(NLAYERS,IB)*DEG_TO_RAD)
                 ELSE
                  X0_BOA = X0(IB)
                  ENDIF 
                  L_ATTN = L_DELTAU_VERT(Q,K)*DELTAU_SLANT(NLAYERS,K,IB)
                  L_ATTN = - L_ATTN * SOLAR_BEAM_OPDEP(IB)
                ENDIF
                X0_FLUX = FOUR * X0_BOA
                L_ATTN  = L_ATTN * X0_FLUX             
                L_DB_CUMSOURCE(V) = L_ATTN * FACTOR

C  Finish loops over geometries

               ENDDO
              ENDIF
             ENDDO
            ENDDO

C  Upwelling recurrence: transmittance of exact source term
C  --------------------------------------------------------

            NC =  0

C  initialize optical depth loop

            NSTART = NLAYERS
            NUT_PREV = NSTART + 1

C  Main loop over all output optical depths

            DO UTA = N_OUT_USERTAUS, 1, -1

C  Layer index for given optical depth

              NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
              NUT    = NLEVEL + 1

C  Cumulative layer transmittance :
C    loop over layers working upwards to level NUT

              DO N = NSTART, NUT, -1
                NC = NLAYERS + 1 - N

C  If N  = K (varying layer), addition linearization of transmittance
C  If N ne K just transmittance of the linearized source

                IF ( N.EQ.K) THEN
                 DO IB = 1, NBEAMS
                  DO UM = 1, N_USER_STREAMS
                   DO UA = 1, N_USER_RELAZMS
                    V = UMOFF(IB,UM) + UA
                    IF ( DO_SSCORR_OUTGOING ) THEN
                      TR  = UP_LOSTRANS(N,V)
                      LTR = L_UP_LOSTRANS(N,Q,V)
                    ELSE
                      TR  = T_DELT_USERM(N,UM)
                      LTR = L_T_DELT_USERM(N,UM,Q)
                    ENDIF
                    L_DB_CUMSOURCE(V) = TR * L_DB_CUMSOURCE(V)
     &                                  + LTR * DB_CUMSOURCE(V,NC-1) 
                   ENDDO
                  ENDDO
                 ENDDO
                ELSE
                 DO IB = 1, NBEAMS
                  DO UM = 1, N_USER_STREAMS
                   DO UA = 1, N_USER_RELAZMS
                    V = UMOFF(IB,UM) + UA
                    IF ( DO_SSCORR_OUTGOING ) THEN
                      TR  = UP_LOSTRANS(N,V)
                    ELSE
                      TR  = T_DELT_USERM(N,UM)
                    ENDIF
                    L_DB_CUMSOURCE(V) = TR * L_DB_CUMSOURCE(V)
                   ENDDO
                  ENDDO
                 ENDDO
                ENDIF

C  end layer loop

              ENDDO

C  Offgrid output
C  ------------- 

C  Require partial layer transmittance
C  If N = K (varying layer), require additional linearization

              IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN

               UT = OFFGRID_UTAU_OUTINDEX(UTA)
               N  = OFFGRID_UTAU_LAYERIDX(UT)

               IF ( N.EQ.K) THEN
                DO IB = 1, NBEAMS
                 DO UM = 1, N_USER_STREAMS
                  DO UA = 1, N_USER_RELAZMS
                   V = UMOFF(IB,UM) + UA
                   IF ( DO_SSCORR_OUTGOING ) THEN
                    TR  = UP_LOSTRANS_UT(UT,V)
                    LTR = L_UP_LOSTRANS_UT(UT,Q,V)
                   ELSE
                    TR  = T_UTUP_USERM(UT,UM)
                    LTR = L_T_UTUP_USERM(UT,UM,Q)
                   ENDIF
                   PROFILEWF_DB(Q,K,UTA,V) = TR * L_DB_CUMSOURCE(V)
     &                                     + LTR * DB_CUMSOURCE(V,NC)
                  ENDDO
                 ENDDO
                ENDDO
               ELSE IF ( N.NE.K) THEN
                DO IB = 1, NBEAMS
                 DO UM = 1, N_USER_STREAMS
                  DO UA = 1, N_USER_RELAZMS
                   V = UMOFF(IB,UM) + UA
                   IF ( DO_SSCORR_OUTGOING ) THEN
                    TR  = UP_LOSTRANS_UT(UT,V)
                   ELSE
                    TR  = T_UTUP_USERM(UT,UM)
                   ENDIF
                   PROFILEWF_DB(Q,K,UTA,V) = TR * L_DB_CUMSOURCE(V)
                  ENDDO
                 ENDDO
                ENDDO
               ENDIF

C  Ongrid output : Set final cumulative source directly

              ELSE

               DO IB = 1, NBEAMS
                DO UM = 1, N_USER_STREAMS
                 DO UA = 1, N_USER_RELAZMS
                  V = UMOFF(IB,UM) + UA
                  PROFILEWF_DB(Q,K,UTA,V) = L_DB_CUMSOURCE(V)
                 ENDDO
                ENDDO
               ENDDO

              ENDIF

C  Check for updating the recursion 

             IF ( NUT. NE. NUT_PREV ) NSTART = NUT - 1
             NUT_PREV = NUT

C  end optical depth loop

            ENDDO

C  Continuation point for avoiding calculation

 6543       continue

C  End weighting function loop

          ENDDO
        ENDIF
      ENDDO

C  Finish

      RETURN
      END

c

      SUBROUTINE LIDORT_LAC_LAMBERTIAN_DBCORRECTION(FLUXMULT)

C  Prepares Linearization of Exact Direct Beam reflection (Lambertian case)
C   Linearization with respect to Column atmospheric variables

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  Include file of setup variables (Input to the present module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  Include files of Linearized input and setup variables

      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'

C  Include files of Single scatter & surface reflectance variables (Input)

      INCLUDE '../includes/LIDORT_SINGSCAT.VARS'
      INCLUDE '../includes/LIDORT_REFLECTANCE.VARS'

C  Direct beam output stored in linearized correction file

      INCLUDE '../includes/LIDORT_L_SINGSCAT.VARS'

C  input argument

      DOUBLE PRECISION FLUXMULT

C  Local variables
C  ---------------

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER          UT, UTA, UM, UA, NC, IB, V, K, Q, K0
      DOUBLE PRECISION X0_FLUX, X0_BOA, FACTOR, L_ATTN, TR, LTR

C  first stage
C  -----------

C  return if no upwelling

      IF ( .NOT.DO_UPWELLING ) RETURN

C  Lambertian factor

      FACTOR = LAMBERTIAN_ALBEDO * FLUXMULT

C  Weighting function index

      K0 = 0

C  Start weighting function loop

      DO Q = 1, N_TOTALCOLUMN_WFS

C  Initialize the output results

        DO V = 1, N_GEOMETRIES
          DO UTA = 1, N_OUT_USERTAUS
            COLUMNWF_DB(Q,UTA,V) = ZERO
          ENDDO
        ENDDO

C  Skip if Lambertian albedo is zero

        IF ( FACTOR .EQ. ZERO ) GO TO 6543

C  New Code  R. Spurr, 6 August 2007. RT Solutions Inc.
C  ====================================================

C  initialize cumulative source term

C  Start geometry loops

        DO UM = 1, N_USER_STREAMS
          DO IB = 1, NBEAMS
            IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
              DO UA = 1, N_USER_RELAZMS
                V = UMOFF(IB,UM) + UA

c  Beam attenuation and reflection

                IF ( DO_SSCORR_OUTGOING ) THEN
                  X0_BOA = DCOS(BEAM_SZAS_ADJUST(UM,IB,UA)*DEG_TO_RAD)
                  L_ATTN = L_BOA_ATTN(K0,Q,V)
                ELSE
                  IF ( DO_REFRACTIVE_GEOMETRY ) THEN
                   X0_BOA = DCOS(SZA_LOCAL_INPUT(NLAYERS,IB)*DEG_TO_RAD)
                  ELSE
                   X0_BOA = X0(IB)
                  ENDIF 
                  L_ATTN = ZERO
                  DO K = 1, NLAYERS
                    L_ATTN = L_ATTN - L_DELTAU_VERT(Q,K) * 
     &                                  DELTAU_SLANT(NLAYERS,K,IB)
                  ENDDO   
                  L_ATTN = L_ATTN * SOLAR_BEAM_OPDEP(IB)
                ENDIF
                X0_FLUX = FOUR * X0_BOA
                L_ATTN  = L_ATTN * X0_FLUX             
                L_DB_CUMSOURCE(V) = L_ATTN * FACTOR

C  Finish loops over geometries

              ENDDO
            ENDIF
          ENDDO
        ENDDO

C  Upwelling recurrence: transmittance of exact source term
C  --------------------------------------------------------

C  initialise

        NC =  0

C  initialize optical depth loop

        NSTART = NLAYERS
        NUT_PREV = NSTART + 1

C  Main loop over all output optical depths

        DO UTA = N_OUT_USERTAUS, 1, -1

C  Layer index for given optical depth

          NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
          NUT    = NLEVEL + 1

C  Cumulative layer transmittance :
C    loop over layers working upwards to level NUT

          DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                DO UA = 1, N_USER_RELAZMS
                  V = UMOFF(IB,UM) + UA
                  IF ( DO_SSCORR_OUTGOING ) THEN
                    TR  = UP_LOSTRANS(N,V)
                    LTR = L_UP_LOSTRANS(N,Q,V)
                  ELSE
                    TR  = T_DELT_USERM(N,UM)
                    LTR = L_T_DELT_USERM(N,UM,Q)
                  ENDIF
                  L_DB_CUMSOURCE(V) = TR * L_DB_CUMSOURCE(V)
     &                                + LTR * DB_CUMSOURCE(V,NC-1) 
                ENDDO
              ENDDO
            ENDDO
          ENDDO

C  Offgrid output
C  ------------- 

C  Require partial layer transmittance

          IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN
            UT = OFFGRID_UTAU_OUTINDEX(UTA)
            N  = OFFGRID_UTAU_LAYERIDX(UT)
            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                DO UA = 1, N_USER_RELAZMS
                  V = UMOFF(IB,UM) + UA
                  IF ( DO_SSCORR_OUTGOING ) THEN
                   TR  = UP_LOSTRANS_UT(UT,V)
                   LTR = L_UP_LOSTRANS_UT(UT,Q,V)
                  ELSE
                   TR  = T_UTUP_USERM(UT,UM)
                   LTR = L_T_UTUP_USERM(UT,UM,Q)
                  ENDIF
                  COLUMNWF_DB(Q,UTA,V) = TR * L_DB_CUMSOURCE(V)
     &                                + LTR *   DB_CUMSOURCE(V,NC)
                ENDDO
              ENDDO
            ENDDO

C  Ongrid output : Set final cumulative source directly

          ELSE

            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                DO UA = 1, N_USER_RELAZMS
                  V = UMOFF(IB,UM) + UA
                  COLUMNWF_DB(Q,UTA,V) = L_DB_CUMSOURCE(V)
                ENDDO
              ENDDO
            ENDDO

          ENDIF

C  Check for updating the recursion 

          IF ( NUT. NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT

C  end optical depth loop

        ENDDO

C  Continuation point for avoiding calculation

 6543   continue

C  End weighting function loop

      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_LS_LAMBERTIAN_DBCORRECTION

C  Prepares Linearization of Exact Direct Beam reflection (Lambertian case)
C   Linearization with respect to LAMBERTIAN amplitude variable only !

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  Include file of setup variables (Input to the present module)

      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  Include files of Linearized input and setup variables

      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_L_SETUPS.VARS'

C  Include files of Single scatter stuff

      INCLUDE '../includes/LIDORT_SINGSCAT.VARS'

C  Direct beam output stored in linearized correction file

      INCLUDE '../includes/LIDORT_L_SINGSCAT.VARS'

C  Local variables
C  ---------------

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER          UT, UTA, UM, UA, NC, IB, V, WF
      DOUBLE PRECISION TR

C  first stage
C  -----------

C  return if no upwelling

      IF ( .NOT.DO_UPWELLING ) RETURN

C  only one Surface weighting function

      WF = 1

C  initialise output

      DO V = 1, N_GEOMETRIES
        DO UTA = 1, N_OUT_USERTAUS
          SURFACEWF_DB(WF,UTA,V)  = ZERO
        ENDDO
      ENDDO

C  Transmittance of source term: upwelling recursion
C  -------------------------------------------------

C  initialize cumulative linearized source term
C    In this case, just equal to the original source term
C      ( as dependence on Lambertin albedo is linear )

      NC =  0
      DO V = 1, N_GEOMETRIES
        LS_DB_CUMSOURCE(V) = DB_CUMSOURCE(V,NC)
      ENDDO

C  initialize optical depth loop

      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

C  Main loop over all output optical depths

      DO UTA = N_OUT_USERTAUS, 1, -1

C  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
        NUT    = NLEVEL + 1

C  Cumulative layer transmittance :
C    loop over layers working upwards to level NUT

        DO N = NSTART, NUT, -1
          NC = NLAYERS + 1 - N
          DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
              DO UA = 1, N_USER_RELAZMS
                V = UMOFF(IB,UM) + UA
                IF ( DO_SSCORR_OUTGOING ) THEN
                  TR  = UP_LOSTRANS(N,V)
                ELSE IF ( DO_SSCORR_NADIR ) THEN
                  TR  = T_DELT_USERM(N,UM)
                ENDIF
                LS_DB_CUMSOURCE(V) = TR * LS_DB_CUMSOURCE(V)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

C  Offgrid output
C  -------------- 

C  Require partial layer transmittance

        IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN
          UT = OFFGRID_UTAU_OUTINDEX(UTA)
          N  = OFFGRID_UTAU_LAYERIDX(UT)
          DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
              DO UA = 1, N_USER_RELAZMS
                V = UMOFF(IB,UM) + UA
                IF ( DO_SSCORR_OUTGOING ) THEN
                  TR  = UP_LOSTRANS_UT(UT,V)
                ELSE IF ( DO_SSCORR_NADIR ) THEN
                  TR  = T_UTUP_USERM(UT,UM)
                ENDIF
                SURFACEWF_DB(WF,UTA,V) = TR * LS_DB_CUMSOURCE(V)
              ENDDO
            ENDDO
          ENDDO

C  Ongrid output : Set final cumulative source directly

        ELSE
          DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
              DO UA = 1, N_USER_RELAZMS
                V = UMOFF(IB,UM) + UA
                SURFACEWF_DB(WF,UTA,V) = LS_DB_CUMSOURCE(V)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

C  Check for updating the recursion 

        IF ( NUT. NE. NUT_PREV ) NSTART = NUT - 1
        NUT_PREV = NUT

C  end optical depth loop

      ENDDO

C  Finish

      RETURN
      END

