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
C # --------------------------                                  #
C #                                                             #
C #   1. Nadir SS-correction                                    #
C #       Version 2.4.  LC module  (column Jacobians)           #
C #                                                             #
C #            VLIDORT_LP_SSCORR_NADIR (master, renamed)        #
C #            VLIDORT_LC_SSCORR_NADIR (master, new 2.4)        #
C #              VLIDORTSS_L_ZMATRICES                          #
C #              VLIDORTSS_L_FMATRICES_REFR (placeholder)       #
C #                                                             #
C #   2. Outgoing SS-correction                                 #
C #      Version 2.2   Whole layer integration                  #
C #              2.3.  partial-layer integration                #
C #      Version 2.4.  Column Jacobians introduced              #
C #                                                             #
C #            VLIDORT_L_SSCORR_OUTGOING (master)               #
C #                 L_outgoing_integration_up                   #
C #                 L_outgoing_integration_dn                   #
C #                                                             #
C #   3. DB correction, Lambertian surfaces                     #
C #      Version 2.4.  LAC module (column Jacobians)            #
C #                                                             #
C #            VLIDORT_LAP_LAMBERTIAN_DBCORR                    #
C #            VLIDORT_LAC_LAMBERTIAN_DBCORR                    #
C #            VLIDORT_LS_LAMBERTIAN_DBCORR                     #
C #                                                             #
C #      Version 2.4R ------- Notes -------                     #
C #                                                             #
C #            Additional correction to ZMATRIX and FMATRIX     #
C #            routines for Azimuth > 180. HvdM, Eqs(94/95)     #
C #            Implemented by V. Natraj and R. Spurr, 5/1/09    #
C #                                                             #
C ###############################################################

      SUBROUTINE VLIDORT_LP_SSCORR_NADIR
     &      ( SSFLUX, LAYER_TO_VARY, NV_PARAMETERS)

C  Single scatter exact calculation. Nadir view
C   Programmed by R. Spurr, RT Solutions Inc.
C    Second Draft, April 14th 2005.
C    Third Draft,  May    6th 2005.
C       - additional code to deal with refraction.

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  Include files lineaeized input variables

      INCLUDE '../includes/VLIDORT_L_INPUTS.VARS'

C  include files of setup and multiplier variables (input)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  Include file of single scatter result variables

      INCLUDE '../includes/VLIDORT_SINGSCAT.VARS'

C  include files of linearized setup, solution and multiplier variables (input)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'

C  Include file of results (linearized single scatter weighting funcs)

      INCLUDE '../includes/VLIDORT_L_SINGSCAT.VARS'

C  Input arguments
C  ---------------

C  Flux factor

      DOUBLE PRECISION SSFLUX

C  Linearization control

      INTEGER          LAYER_TO_VARY, NV_PARAMETERS

C  local variables
C  ---------------

C  indices

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL, V
      INTEGER          UT, UTA, UM, IB, IA, NC, Q, NM1, NV, O1, O2

C  other help variables

      DOUBLE PRECISION HELP, TTT, UVAR, AVAR, FT1, FT2, DNM1
      DOUBLE PRECISION L_FINAL_SOURCE, L_SS_LAYERSOURCE
      DOUBLE PRECISION L_SSCORRECTION, L_TR_CUMSOURCE
      DOUBLE PRECISION VAR_TMS(MAX_ATMOSWFS), MUX, SUM, VIEWSIGN

C  May not require these
C  ---------------------

C  zenith angle cosines/sines, azimuth angle cosines

      DOUBLE PRECISION CTHETA (MAXLAYERS,MAXBEAMS)
      DOUBLE PRECISION STHETA (MAXLAYERS,MAXBEAMS)
      DOUBLE PRECISION CALPHA (MAX_USER_STREAMS)
      DOUBLE PRECISION SALPHA (MAX_USER_STREAMS)
      DOUBLE PRECISION CPHI   (MAX_USER_RELAZMS)

C  Set up operations
C  -----------------

      NV = LAYER_TO_VARY

C  Create Linearized TMS factor for Layer NV
C   ( Use UNSCALED linearized inputs - the Nakajima-Tanaka way)

      IF ( DO_SSCORR_TRUNCATION ) THEN
        TTT = TMS(NV) /  ( ONE - SSFDEL(NV) )
      ELSE
        TTT = TMS(NV)
      ENDIF
 
      IF ( DO_DELTAM_SCALING ) THEN
        NM1 = NMOMENTS+1
        FT1 = TRUNC_FACTOR(NV) * TTT
        FT2 = ONE + FT1
        DO Q = 1, NV_PARAMETERS
          IF ( DO_SCATMAT_VARIATION(NV,Q) ) THEN
            UVAR = L_OMEGA_TOTAL_INPUT(Q,NV)
            AVAR = UVAR + L_GREEKMAT_TOTAL_INPUT(Q,NM1,NV,1)
            VAR_TMS(Q) = UVAR + FT1 * AVAR
          ELSE
            UVAR = L_OMEGA_TOTAL_INPUT(Q,NV)
            VAR_TMS(Q) = UVAR * FT2
          ENDIF
        ENDDO
      ELSE
        DO Q = 1, NV_PARAMETERS
          UVAR = L_OMEGA_TOTAL_INPUT(Q,NV)
          VAR_TMS(Q) = UVAR
        ENDDO
      ENDIF

C  Additional Delta-M scaling
C  --------------------------

C  New section. R. Spurr, 07 September 2007.

      IF ( DO_SSCORR_TRUNCATION ) THEN
        NM1  = NGREEK_MOMENTS_INPUT
        DNM1 = DFLOAT(2*NM1+1)
        DO Q = 1, NV_PARAMETERS
          IF ( DO_SCATMAT_VARIATION(NV,Q) ) THEN
            FT1 = ONE / ( ONE- SSFDEL(NV) )
            L_SSFDEL(NV,Q) = GREEKMAT_TOTAL_INPUT(NM1,NV,1) *
     *                  L_GREEKMAT_TOTAL_INPUT(Q,NM1,NV,1) / DNM1
            VAR_TMS(Q) = VAR_TMS(Q) - FT1 * L_SSFDEL(NV,Q)
          ENDIF
        ENDDO
      ENDIF

C  save some geometrical quantities (cosines and sines)
C    Multiple layer quantities for the Refractive case

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN
        DO IB = 1, NBEAMS
          DO N = 1, NLAYERS
            MUX =  SUN_SZA_COSINES(N,IB)  
            CTHETA(N,IB) = MUX
            STHETA(N,IB) = DSQRT(ONE-MUX*MUX)
          ENDDO
        ENDDO
      ELSE
        DO IB = 1, NBEAMS
          DO N = 1, NLAYERS
            CTHETA(N,IB) = COS_SZANGLES(IB)
            STHETA(N,IB) = SIN_SZANGLES(IB)
          END DO
        END DO
      ENDIF

      DO UM = 1, N_USER_STREAMS
        CALPHA(UM) = USER_STREAMS(UM)
        SALPHA(UM) = DSQRT ( ONE - CALPHA(UM) * CALPHA(UM) )
      ENDDO

      DO IA = 1, N_USER_RELAZMS
        CPHI(IA) = DCOS ( USER_RELAZMS(IA) * DEG_TO_RAD )
      ENDDO

C  ####################
C  #    UPWELLING     #
C  ####################

      IF ( DO_UPWELLING ) THEN

C  =================================================
C  Total scattering matrix linearization (upwelling)
C  =================================================

        IF ( STERM_LAYERMASK_UP(NV)) THEN

C  Get L_Z matrices for each layer, if the scatter law is linearized

         VIEWSIGN = -1.0D0
         CALL VLIDORTSS_L_ZMATRICES
     I  ( NSTOKES, NV, NV_PARAMETERS, DO_SCATMAT_VARIATION,
     I    N_GEOMETRIES, NGREEK_MOMENTS_INPUT,
     I    NBEAMS, N_USER_STREAMS,  N_USER_RELAZMS,
     I    LAYER_MAXMOMENTS,
     I    GREEKMAT_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT,
     I    DO_SSCORR_TRUNCATION, SSFDEL, L_SSFDEL,
     I    VIEWSIGN, VZA_OFFSETS, 
     I    CTHETA, STHETA, CALPHA, SALPHA, CPHI, USER_RELAZMS,
     O    L_ZMAT_UP )

C  Loop over varying parameters Q for layer NV
C  Loop over all geometries V 
  
         DO Q = 1, NV_PARAMETERS

          DO V = 1, N_GEOMETRIES

C  Phase function moment variations
C    add TMS correction factor linearization

            IF ( DO_SCATMAT_VARIATION(NV,Q) ) THEN
              DO O1 = 1, NSTOKES
                DO O2 = 1, NSTOKES
                  L_ZMAT_UP(Q,V,NV,O1,O2) =
     &                L_ZMAT_UP(Q,V,NV,O1,O2)  * TMS(NV) +
     &                  ZMAT_UP(V,NV,O1,O2)    * VAR_TMS(Q)
                ENDDO
              ENDDO
            ELSE
              DO O1 = 1, NSTOKES
                DO O2 = 1, NSTOKES
                  L_ZMAT_UP(Q,V,NV,O1,O2) =
     &                  ZMAT_UP(V,NV,O1,O2) * VAR_TMS(Q)
                ENDDO
              ENDDO
            ENDIF

C  end parameter loop and viewing directions loop

          ENDDO
         ENDDO

C  Only if layer NV exists

        ENDIF

C  ===================================
C  Upwelling single scatter recurrence
C  ===================================

C  initialize cumulative source term

        DO V = 1, N_GEOMETRIES
          DO Q = 1, NV_PARAMETERS
            DO O1 = 1, NSTOKES
              L_SS_CUMSOURCE(Q,V,O1) = ZERO
            ENDDO
          ENDDO
        ENDDO

C  initialise optical depth loop

        NC = 0
        NSTART = NLAYERS
        NUT_PREV = NSTART + 1

C  Main loop over all output optical depths
C  ========================================

        DO UTA = N_USER_LEVELS, 1, -1

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

C  If N = NV (the layer that is varying)

            IF ( N .EQ. NV ) THEN

              DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO Q = 1, NV_PARAMETERS
                  DO O1 = 1, NSTOKES
                   SUM = ZERO
                   DO O2 = 1, NSTOKES
                    HELP = ZMAT_UP(V,N,O1,O2) * L_EMULT_UP(UM,N,NV,IB,Q)
     &                 + L_ZMAT_UP(Q,V,N,O1,O2) * EMULT_UP(UM,N,IB)
                    SUM = SUM + HELP * FLUXVEC(O2)
                   ENDDO
                   L_SS_LAYERSOURCE = SUM
                   L_SS_CUMSOURCE(Q,V,O1) = L_SS_LAYERSOURCE
     &         +   T_DELT_USERM(N,UM)   * L_SS_CUMSOURCE(Q,V,O1)
     &         + L_T_DELT_USERM(N,UM,Q) *   SS_CUMSOURCE_UP(V,O1,NC-1)
                  ENDDO
                 ENDDO
                ENDDO
               ENDDO
              ENDDO

C  Variations when N > NV

            ELSE IF ( N .GT. NV ) THEN

              DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO Q = 1, NV_PARAMETERS
                  DO O1 = 1, NSTOKES
                   SUM = ZERO
                   DO O2 = 1, NSTOKES
                    HELP = ZMAT_UP(V,N,O1,O2) * L_EMULT_UP(UM,N,NV,IB,Q)
                    SUM = SUM + HELP * FLUXVEC(O2)
                   ENDDO
                   L_SS_LAYERSOURCE = SUM
                   L_SS_CUMSOURCE(Q,V,O1) = L_SS_LAYERSOURCE
     &         +   T_DELT_USERM(N,UM)   * L_SS_CUMSOURCE(Q,V,O1)
                  ENDDO
                 ENDDO
                ENDDO
               ENDDO
              ENDDO

C  Transmittance for N < NV

            ELSE

              DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO Q = 1, NV_PARAMETERS
                  DO O1 = 1, NSTOKES
                   L_SS_CUMSOURCE(Q,V,O1) = 
     &               T_DELT_USERM(N,UM) * L_SS_CUMSOURCE(Q,V,O1)
                  ENDDO
                 ENDDO
                ENDDO
               ENDDO
              ENDDO

            ENDIF

          ENDDO

C  Offgrid output
C  --------------

C  Set final cumulative source and Single scatter Weighting function

          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)

C  If N = NV (the layer that is varying)
C    add Linearization of additional partial layer source term =
C        L_Exact_Scat(n) * Multiplier  +  Exact_Scat * L_Multiplier(n)

            IF ( N .EQ. NV ) THEN

              DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO Q = 1, NV_PARAMETERS
                  DO O1 = 1, NSTOKES
                   SUM = ZERO
                   DO O2 = 1, NSTOKES
                    HELP =
     &            ZMAT_UP(V,N,O1,O2) * L_UT_EMULT_UP(UM,UT,NV,IB,Q)
     &        + L_ZMAT_UP(Q,V,N,O1,O2) * UT_EMULT_UP(UM,UT,IB)
                    SUM = SUM + HELP * FLUXVEC(O2)
                   ENDDO
                   L_SS_LAYERSOURCE = SUM
                   L_TR_CUMSOURCE = 
     &           + L_T_UTUP_USERM(UT,UM,Q) *   SS_CUMSOURCE_UP(V,O1,NC)
     &           +   T_UTUP_USERM(UT,UM)   * L_SS_CUMSOURCE(Q,V,O1)
                   L_FINAL_SOURCE = L_TR_CUMSOURCE + L_SS_LAYERSOURCE
                   L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                   PROFILEWF_SS(Q,NV,UTA,V,O1,UPIDX) = L_SSCORRECTION
                  ENDDO
                 ENDDO
                ENDDO
               ENDDO
              ENDDO

C  Variations when N > NV
C    add Linearization of additional partial layer source term =
C         Exact_Scat * L_Multiplier(NV)

            ELSE IF ( N .GT. NV ) THEN

              DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO Q = 1, NV_PARAMETERS
                  DO O1 = 1, NSTOKES
                   SUM = ZERO
                   DO O2 = 1, NSTOKES
                    HELP =
     &            ZMAT_UP(V,N,O1,O2) * L_UT_EMULT_UP(UM,UT,NV,IB,Q)
                    SUM = SUM + HELP * FLUXVEC(O2)
                   ENDDO
                   L_SS_LAYERSOURCE = SUM
                   L_TR_CUMSOURCE = 
     &           +   T_UTUP_USERM(UT,UM)   * L_SS_CUMSOURCE(Q,V,O1)
                   L_FINAL_SOURCE = L_TR_CUMSOURCE + L_SS_LAYERSOURCE
                   L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                   PROFILEWF_SS(Q,NV,UTA,V,O1,UPIDX) = L_SSCORRECTION
                  ENDDO
                 ENDDO
                ENDDO
               ENDDO
              ENDDO

C  Transmittance for the other layers

            ELSE

              DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO Q = 1, NV_PARAMETERS
                  DO O1 = 1, NSTOKES
                   L_FINAL_SOURCE = 
     &           +   T_UTUP_USERM(UT,UM)   * L_SS_CUMSOURCE(Q,V,O1)
                   L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                   PROFILEWF_SS(Q,NV,UTA,V,O1,UPIDX) = L_SSCORRECTION
                  ENDDO
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
              DO Q = 1, NV_PARAMETERS
                DO O1 = 1, NSTOKES
                  L_FINAL_SOURCE = L_SS_CUMSOURCE(Q,V,O1)
                  L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                  PROFILEWF_SS(Q,NV,UTA,V,O1,UPIDX) = L_SSCORRECTION
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

C  =================================================
C  Total scattering matrix linearization (upwelling)
C  =================================================

        IF ( STERM_LAYERMASK_DN(NV)) THEN

C  Get L_Z matrices for each layer, if the scatter law is linearized

         VIEWSIGN = +1.0D0
         CALL VLIDORTSS_L_ZMATRICES
     I  ( NSTOKES, NV, NV_PARAMETERS, DO_SCATMAT_VARIATION,
     I    N_GEOMETRIES, NGREEK_MOMENTS_INPUT,
     I    NBEAMS, N_USER_STREAMS,  N_USER_RELAZMS,
     I    LAYER_MAXMOMENTS,
     I    GREEKMAT_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT,
     I    DO_SSCORR_TRUNCATION, SSFDEL, L_SSFDEL,
     I    VIEWSIGN, VZA_OFFSETS, 
     I    CTHETA, STHETA, CALPHA, SALPHA, CPHI, USER_RELAZMS,
     O    L_ZMAT_DN )

C  Loop over varying parameters Q for layer NV
C  Loop over all geometries V 

         DO Q = 1, NV_PARAMETERS
          DO V = 1, N_GEOMETRIES

C  Phase function moment variations
C    add TMS correction factor linearization

            IF ( DO_SCATMAT_VARIATION(NV,Q) ) THEN
              DO O1 = 1, NSTOKES
                DO O2 = 1, NSTOKES
                  L_ZMAT_DN(Q,V,NV,O1,O2) =
     &                L_ZMAT_DN(Q,V,NV,O1,O2)  * TMS(NV) +
     &                  ZMAT_DN(V,NV,O1,O2) * VAR_TMS(Q)
                ENDDO
              ENDDO
            ELSE
              DO O1 = 1, NSTOKES
                DO O2 = 1, NSTOKES
                  L_ZMAT_DN(Q,V,NV,O1,O2) =
     &                  ZMAT_DN(V,NV,O1,O2) * VAR_TMS(Q)
                ENDDO
              ENDDO
            ENDIF

C  end parameter and geometry loops

          ENDDO
         ENDDO

C  Only if layer NV exists

        ENDIF

C  =====================================
C  Downwelling single scatter recurrence
C  =====================================

C  initialize cumulative source term

        DO V = 1, N_GEOMETRIES
          DO Q = 1, NV_PARAMETERS
            DO O1 = 1, NSTOKES
              L_SS_CUMSOURCE(Q,V,O1) = ZERO
            ENDDO
          ENDDO
        ENDDO

C  initialise optical depth loop

        NC = 0
        NSTART = 1
        NUT_PREV = NSTART - 1

C  Main loop over all output optical depths
C  ========================================

        DO UTA = 1, N_USER_LEVELS

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

C  If N = NV (the layer that is varying)

            IF ( N .EQ. NV ) THEN

              DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO Q = 1, NV_PARAMETERS
                  DO O1 = 1, NSTOKES
                   SUM = ZERO
                   DO O2 = 1, NSTOKES
                    HELP = ZMAT_DN(V,N,O1,O2) * L_EMULT_DN(UM,N,NV,IB,Q)
     &                 + L_ZMAT_DN(Q,V,N,O1,O2) * EMULT_DN(UM,N,IB)
                    SUM = SUM + HELP * FLUXVEC(O2)
                   ENDDO
                   L_SS_LAYERSOURCE = SUM
                   L_SS_CUMSOURCE(Q,V,O1) = L_SS_LAYERSOURCE
     &         +   T_DELT_USERM(N,UM)   * L_SS_CUMSOURCE(Q,V,O1)
     &         + L_T_DELT_USERM(N,UM,Q) *   SS_CUMSOURCE_DN(V,O1,NC-1)
                  ENDDO
                 ENDDO
                ENDDO
               ENDDO
              ENDDO

C  Variations when N > NV

            ELSE IF ( N .GT. NV ) THEN

              DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO Q = 1, NV_PARAMETERS
                  DO O1 = 1, NSTOKES
                   SUM = ZERO
                   DO O2 = 1, NSTOKES
                    HELP = ZMAT_DN(V,N,O1,O2) * L_EMULT_DN(UM,N,NV,IB,Q)
                    SUM = SUM + HELP * FLUXVEC(O2)
                   ENDDO
                   L_SS_LAYERSOURCE = SUM
                   L_SS_CUMSOURCE(Q,V,O1) = L_SS_LAYERSOURCE
     &         +   T_DELT_USERM(N,UM)   * L_SS_CUMSOURCE(Q,V,O1)
                  ENDDO
                 ENDDO
                ENDDO
               ENDDO
              ENDDO

C  Transmittance for N < NV

            ELSE

              DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO Q = 1, NV_PARAMETERS
                  DO O1 = 1, NSTOKES
                   L_SS_CUMSOURCE(Q,V,O1) = 
     &               T_DELT_USERM(N,UM) * L_SS_CUMSOURCE(Q,V,O1)
                  ENDDO
                 ENDDO
                ENDDO
               ENDDO
              ENDDO

            ENDIF

C  end recursive layer loop

          ENDDO

C  Offgrid output
C  --------------

C  add additional partial layer source term = Exact Scat * Multiplier
C  Set final cumulative source and Correct the intensity

          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)

C  If N = NV (the layer that is varying)
C    add Linearization of additional partial layer source term =
C        L_Exact_Scat(n) * Multiplier  +  Exact_Scat * L_Multiplier(n)

            IF ( N .EQ. NV ) THEN

              DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO Q = 1, NV_PARAMETERS
                  DO O1 = 1, NSTOKES
                   SUM = ZERO
                   DO O2 = 1, NSTOKES
                    HELP =
     &            ZMAT_DN(V,N,O1,O2) * L_UT_EMULT_DN(UM,UT,NV,IB,Q)
     &        + L_ZMAT_DN(Q,V,N,O1,O2) * UT_EMULT_DN(UM,UT,IB)
                    SUM = SUM + HELP * FLUXVEC(O2)
                   ENDDO
                   L_SS_LAYERSOURCE = SUM
                   L_TR_CUMSOURCE = 
     &           + L_T_UTDN_USERM(UT,UM,Q) *   SS_CUMSOURCE_DN(V,O1,NC)
     &           +   T_UTDN_USERM(UT,UM)   * L_SS_CUMSOURCE(Q,V,O1)
                   L_FINAL_SOURCE = L_TR_CUMSOURCE + L_SS_LAYERSOURCE
                   L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                   PROFILEWF_SS(Q,NV,UTA,V,O1,DNIDX) = L_SSCORRECTION
                  ENDDO
                 ENDDO
                ENDDO
               ENDDO
              ENDDO

C  Variations when N > NV
C    add Linearization of additional partial layer source term =
C         Exact_Scat * L_Multiplier(NV)

            ELSE IF ( N .GT. NV ) THEN

              DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO Q = 1, NV_PARAMETERS
                  DO O1 = 1, NSTOKES
                   SUM = ZERO
                   DO O2 = 1, NSTOKES
                    HELP =
     &            ZMAT_DN(V,N,O1,O2) * L_UT_EMULT_DN(UM,UT,NV,IB,Q)
                    SUM = SUM + HELP * FLUXVEC(O2)
                   ENDDO
                   L_SS_LAYERSOURCE = SUM
                   L_TR_CUMSOURCE = 
     &           +   T_UTDN_USERM(UT,UM)   * L_SS_CUMSOURCE(Q,V,O1)
                   L_FINAL_SOURCE = L_TR_CUMSOURCE + L_SS_LAYERSOURCE
                   L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                   PROFILEWF_SS(Q,NV,UTA,V,O1,DNIDX) = L_SSCORRECTION
                  ENDDO
                 ENDDO
                ENDDO
               ENDDO
              ENDDO

C  Transmittance for the other layers

            ELSE

              DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO Q = 1, NV_PARAMETERS
                  DO O1 = 1, NSTOKES
                   L_FINAL_SOURCE = 
     &            +   T_UTDN_USERM(UT,UM)   * L_SS_CUMSOURCE(Q,V,O1)
                   L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                   PROFILEWF_SS(Q,NV,UTA,V,O1,DNIDX) = L_SSCORRECTION
                  ENDDO
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
              DO Q = 1, NV_PARAMETERS
                DO O1 = 1, NSTOKES
                  L_FINAL_SOURCE = L_SS_CUMSOURCE(Q,V,O1)
                  L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                  PROFILEWF_SS(Q,NV,UTA,V,O1,DNIDX) = L_SSCORRECTION
                ENDDO
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

      SUBROUTINE VLIDORT_LC_SSCORR_NADIR
     &      ( SSFLUX, NV_PARAMETERS)

C  Single scatter exact calculation. Nadir view
C   Programmed by R. Spurr, RT Solutions Inc.

C  Version 2.4. Total Column Jacobians, December 2008
C               Modeled after LIDORT code Version 3.3

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  Include files lineaeized input variables

      INCLUDE '../includes/VLIDORT_L_INPUTS.VARS'

C  include files of setup and multiplier variables (input)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  Include file of single scatter result variables

      INCLUDE '../includes/VLIDORT_SINGSCAT.VARS'

C  include files of linearized setup, solution and multiplier variables (input)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'

C  Include file of results (linearized single scatter weighting funcs)

      INCLUDE '../includes/VLIDORT_L_SINGSCAT.VARS'

C  Input arguments
C  ---------------

C  Flux factor

      DOUBLE PRECISION SSFLUX

C  Linearization control

      INTEGER          NV_PARAMETERS

C  local variables
C  ---------------

C  indices

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL, V
      INTEGER          UT, UTA, UM, IB, IA, NC, Q, NM1, NV, O1, O2

C  other help variables

      DOUBLE PRECISION HELP, TTT, UVAR, AVAR, FT1, FT2, DNM1
      DOUBLE PRECISION L_FINAL_SOURCE, L_SS_LAYERSOURCE
      DOUBLE PRECISION L_SSCORRECTION, L_TR_CUMSOURCE
      DOUBLE PRECISION VAR_TMS(MAX_ATMOSWFS,MAXLAYERS)
      DOUBLE PRECISION  MUX, SUM, VIEWSIGN

C  May not require these
C  ---------------------

C  zenith angle cosines/sines, azimuth angle cosines

      DOUBLE PRECISION CTHETA (MAXLAYERS,MAXBEAMS)
      DOUBLE PRECISION STHETA (MAXLAYERS,MAXBEAMS)
      DOUBLE PRECISION CALPHA (MAX_USER_STREAMS)
      DOUBLE PRECISION SALPHA (MAX_USER_STREAMS)
      DOUBLE PRECISION CPHI   (MAX_USER_RELAZMS)

C  Set up operations
C  -----------------

C  All layers

      DO NV = 1, NLAYERS

C  Create Linearized TMS factor for Layer NV
C   ( Use UNSCALED linearized inputs - the Nakajima-Tanaka way)

        IF ( DO_SSCORR_TRUNCATION ) THEN
          TTT = TMS(NV) /  ( ONE - SSFDEL(NV) )
        ELSE
          TTT = TMS(NV)
        ENDIF
 
        IF ( DO_DELTAM_SCALING ) THEN
          NM1 = NMOMENTS+1
          FT1 = TRUNC_FACTOR(NV) * TTT
          FT2 = ONE + FT1
          DO Q = 1, NV_PARAMETERS
            IF ( DO_SCATMAT_VARIATION(NV,Q) ) THEN
              UVAR = L_OMEGA_TOTAL_INPUT(Q,NV)
              AVAR = UVAR + L_GREEKMAT_TOTAL_INPUT(Q,NM1,NV,1)
              VAR_TMS(Q,NV) = UVAR + FT1 * AVAR
            ELSE
              UVAR = L_OMEGA_TOTAL_INPUT(Q,NV)
              VAR_TMS(Q,NV) = UVAR * FT2
            ENDIF
          ENDDO
        ELSE
          DO Q = 1, NV_PARAMETERS
            UVAR = L_OMEGA_TOTAL_INPUT(Q,NV)
            VAR_TMS(Q,NV) = UVAR
          ENDDO
        ENDIF

C  Additional Delta-M scaling
C  --------------------------

C  New section. R. Spurr, 07 September 2007.

        IF ( DO_SSCORR_TRUNCATION ) THEN
          NM1  = NGREEK_MOMENTS_INPUT
          DNM1 = DFLOAT(2*NM1+1)
          DO Q = 1, NV_PARAMETERS
            IF ( DO_SCATMAT_VARIATION(NV,Q) ) THEN
              FT1 = ONE / ( ONE- SSFDEL(NV) )
              L_SSFDEL(NV,Q) = GREEKMAT_TOTAL_INPUT(NM1,NV,1) *
     *                  L_GREEKMAT_TOTAL_INPUT(Q,NM1,NV,1) / DNM1
              VAR_TMS(Q,NV) = VAR_TMS(Q,NV) - FT1 * L_SSFDEL(NV,Q)
            ENDIF
          ENDDO
        ENDIF

C  End layer loop

      ENDDO

C  save some geometrical quantities (cosines and sines)
C    Multiple layer quantities for the Refractive case

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN
        DO IB = 1, NBEAMS
          DO N = 1, NLAYERS
            MUX =  SUN_SZA_COSINES(N,IB)  
            CTHETA(N,IB) = MUX
            STHETA(N,IB) = DSQRT(ONE-MUX*MUX)
          ENDDO
        ENDDO
      ELSE
        DO IB = 1, NBEAMS
          DO N = 1, NLAYERS
            CTHETA(N,IB) = COS_SZANGLES(IB)
            STHETA(N,IB) = SIN_SZANGLES(IB)
          END DO
        END DO
      ENDIF

      DO UM = 1, N_USER_STREAMS
        CALPHA(UM) = USER_STREAMS(UM)
        SALPHA(UM) = DSQRT ( ONE - CALPHA(UM) * CALPHA(UM) )
      ENDDO

      DO IA = 1, N_USER_RELAZMS
        CPHI(IA) = DCOS ( USER_RELAZMS(IA) * DEG_TO_RAD )
      ENDDO

C  ####################
C  #    UPWELLING     #
C  ####################

      IF ( DO_UPWELLING ) THEN

C  =================================================
C  Total scattering matrix linearization (upwelling)
C  =================================================

C  All layers

        DO NV = 1, NLAYERS
         IF ( STERM_LAYERMASK_UP(NV)) THEN

C  Get L_Z matrices for each layer, if the scatter law is linearized

          VIEWSIGN = -1.0D0
          CALL VLIDORTSS_L_ZMATRICES
     I  ( NSTOKES, NV, NV_PARAMETERS, DO_SCATMAT_VARIATION,
     I    N_GEOMETRIES, NGREEK_MOMENTS_INPUT,
     I    NBEAMS, N_USER_STREAMS,  N_USER_RELAZMS,
     I    LAYER_MAXMOMENTS,
     I    GREEKMAT_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT,
     I    DO_SSCORR_TRUNCATION, SSFDEL, L_SSFDEL,
     I    VIEWSIGN, VZA_OFFSETS, 
     I    CTHETA, STHETA, CALPHA, SALPHA, CPHI, USER_RELAZMS,
     O    L_ZMAT_UP )

C  Loop over varying parameters Q for layer NV
C  Loop over all geometries V 
  
          DO Q = 1, NV_PARAMETERS
            DO V = 1, N_GEOMETRIES

C  Phase function moment variations
C    add TMS correction factor linearization

              IF ( DO_SCATMAT_VARIATION(NV,Q) ) THEN
                DO O1 = 1, NSTOKES
                  DO O2 = 1, NSTOKES
                    L_ZMAT_UP(Q,V,NV,O1,O2) =
     &                L_ZMAT_UP(Q,V,NV,O1,O2)  * TMS(NV) +
     &                  ZMAT_UP(V,NV,O1,O2)    * VAR_TMS(Q,NV)
                  ENDDO
                ENDDO
              ELSE
                DO O1 = 1, NSTOKES
                  DO O2 = 1, NSTOKES
                    L_ZMAT_UP(Q,V,NV,O1,O2) =
     &                  ZMAT_UP(V,NV,O1,O2) * VAR_TMS(Q,NV)
                  ENDDO
                ENDDO
              ENDIF

C  end parameter loop and viewing directions loop

            ENDDO
          ENDDO

C  Only if layer NV exists

         ENDIF
        ENDDO

C  ===================================
C  Upwelling single scatter recurrence
C  ===================================

C  initialize cumulative source term

        DO V = 1, N_GEOMETRIES
          DO Q = 1, NV_PARAMETERS
            DO O1 = 1, NSTOKES
              L_SS_CUMSOURCE(Q,V,O1) = ZERO
            ENDDO
          ENDDO
        ENDDO

C  initialise optical depth loop

        NC = 0
        NSTART = NLAYERS
        NUT_PREV = NSTART + 1

C  Main loop over all output optical depths
C  ========================================

        DO UTA = N_USER_LEVELS, 1, -1

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

C  ..All layers will have some variation

            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO Q = 1, NV_PARAMETERS
                  DO O1 = 1, NSTOKES
                   SUM = ZERO
                   DO O2 = 1, NSTOKES
                    HELP = ZMAT_UP(V,N,O1,O2) * L_EMULT_UP(UM,N,0,IB,Q)
     &                 + L_ZMAT_UP(Q,V,N,O1,O2) * EMULT_UP(UM,N,IB)
                    SUM = SUM + HELP * FLUXVEC(O2)
                   ENDDO
                   L_SS_LAYERSOURCE = SUM
                   L_SS_CUMSOURCE(Q,V,O1) = L_SS_LAYERSOURCE
     &         +   T_DELT_USERM(N,UM)   * L_SS_CUMSOURCE(Q,V,O1)
     &         + L_T_DELT_USERM(N,UM,Q) *   SS_CUMSOURCE_UP(V,O1,NC-1)
                  ENDDO
                 ENDDO
                ENDDO
              ENDDO
            ENDDO

C  End layer loop

          ENDDO

C  Offgrid output
C  --------------

C  Set final cumulative source and Single scatter Weighting function

          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)

C    add Linearization of additional partial layer source term =
C        L_Exact_Scat(n) * Multiplier  +  Exact_Scat * L_Multiplier(n)

            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO Q = 1, NV_PARAMETERS
                  DO O1 = 1, NSTOKES
                   SUM = ZERO
                   DO O2 = 1, NSTOKES
                    HELP =
     &            ZMAT_UP(V,N,O1,O2) * L_UT_EMULT_UP(UM,UT,0,IB,Q)
     &        + L_ZMAT_UP(Q,V,N,O1,O2) * UT_EMULT_UP(UM,UT,IB)
                    SUM = SUM + HELP * FLUXVEC(O2)
                   ENDDO
                   L_SS_LAYERSOURCE = SUM
                   L_TR_CUMSOURCE = 
     &           + L_T_UTUP_USERM(UT,UM,Q) *   SS_CUMSOURCE_UP(V,O1,NC)
     &           +   T_UTUP_USERM(UT,UM)   * L_SS_CUMSOURCE(Q,V,O1)
                   L_FINAL_SOURCE = L_TR_CUMSOURCE + L_SS_LAYERSOURCE
                   L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                   COLUMNWF_SS(Q,UTA,V,O1,UPIDX) = L_SSCORRECTION
                  ENDDO
                 ENDDO
                ENDDO
              ENDDO
            ENDDO

C  Ongrid output
C  -------------

C  just set to the cumulative source term 

          ELSE

            DO V = 1, N_GEOMETRIES
              DO Q = 1, NV_PARAMETERS
                DO O1 = 1, NSTOKES
                  L_FINAL_SOURCE = L_SS_CUMSOURCE(Q,V,O1)
                  L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                  COLUMNWF_SS(Q,UTA,V,O1,UPIDX) = L_SSCORRECTION
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

C  =================================================
C  Total scattering matrix linearization (upwelling)
C  =================================================

C  All layers

        DO NV = 1, NLAYERS
         IF ( STERM_LAYERMASK_DN(NV)) THEN

C  Get L_Z matrices for each layer, if the scatter law is linearized

          VIEWSIGN = +1.0D0
          CALL VLIDORTSS_L_ZMATRICES
     I  ( NSTOKES, NV, NV_PARAMETERS, DO_SCATMAT_VARIATION,
     I    N_GEOMETRIES, NGREEK_MOMENTS_INPUT,
     I    NBEAMS, N_USER_STREAMS,  N_USER_RELAZMS,
     I    LAYER_MAXMOMENTS,
     I    GREEKMAT_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT,
     I    DO_SSCORR_TRUNCATION, SSFDEL, L_SSFDEL,
     I    VIEWSIGN, VZA_OFFSETS, 
     I    CTHETA, STHETA, CALPHA, SALPHA, CPHI, USER_RELAZMS,
     O    L_ZMAT_DN )

C  Loop over varying parameters Q for layer NV
C  Loop over all geometries V 

          DO Q = 1, NV_PARAMETERS
           DO V = 1, N_GEOMETRIES

C  Phase function moment variations
C    add TMS correction factor linearization

            IF ( DO_SCATMAT_VARIATION(NV,Q) ) THEN
              DO O1 = 1, NSTOKES
                DO O2 = 1, NSTOKES
                  L_ZMAT_DN(Q,V,NV,O1,O2) =
     &                L_ZMAT_DN(Q,V,NV,O1,O2)  * TMS(NV) +
     &                  ZMAT_DN(V,NV,O1,O2) * VAR_TMS(Q,NV)
                ENDDO
              ENDDO
            ELSE
              DO O1 = 1, NSTOKES
                DO O2 = 1, NSTOKES
                  L_ZMAT_DN(Q,V,NV,O1,O2) =
     &                  ZMAT_DN(V,NV,O1,O2) * VAR_TMS(Q,NV)
                ENDDO
              ENDDO
            ENDIF

C  end parameter and geometry loops

           ENDDO
          ENDDO

C  Only if layer NV exists

         ENDIF
        ENDDO

C  =====================================
C  Downwelling single scatter recurrence
C  =====================================

C  initialize cumulative source term

        DO V = 1, N_GEOMETRIES
          DO Q = 1, NV_PARAMETERS
            DO O1 = 1, NSTOKES
              L_SS_CUMSOURCE(Q,V,O1) = ZERO
            ENDDO
          ENDDO
        ENDDO

C  initialise optical depth loop

        NC = 0
        NSTART = 1
        NUT_PREV = NSTART - 1

C  Main loop over all output optical depths
C  ========================================

        DO UTA = 1, N_USER_LEVELS

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
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO Q = 1, NV_PARAMETERS
                  DO O1 = 1, NSTOKES
                   SUM = ZERO
                   DO O2 = 1, NSTOKES
                    HELP = ZMAT_DN(V,N,O1,O2) * L_EMULT_DN(UM,N,0,IB,Q)
     &                 + L_ZMAT_DN(Q,V,N,O1,O2) * EMULT_DN(UM,N,IB)
                    SUM = SUM + HELP * FLUXVEC(O2)
                   ENDDO
                   L_SS_LAYERSOURCE = SUM
                   L_SS_CUMSOURCE(Q,V,O1) = L_SS_LAYERSOURCE
     &         +   T_DELT_USERM(N,UM)   * L_SS_CUMSOURCE(Q,V,O1)
     &         + L_T_DELT_USERM(N,UM,Q) *   SS_CUMSOURCE_DN(V,O1,NC-1)
                  ENDDO
                 ENDDO
                ENDDO
              ENDDO
            ENDDO

C  end recursive layer loop

          ENDDO

C  Offgrid output
C  --------------

C  add additional partial layer source term = Exact Scat * Multiplier
C  Set final cumulative source and Correct the intensity

          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)

C    add Linearization of additional partial layer source term =
C        L_Exact_Scat(n) * Multiplier  +  Exact_Scat * L_Multiplier(n)

            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO Q = 1, NV_PARAMETERS
                  DO O1 = 1, NSTOKES
                   SUM = ZERO
                   DO O2 = 1, NSTOKES
                    HELP =
     &            ZMAT_DN(V,N,O1,O2) * L_UT_EMULT_DN(UM,UT,0,IB,Q)
     &        + L_ZMAT_DN(Q,V,N,O1,O2) * UT_EMULT_DN(UM,UT,IB)
                    SUM = SUM + HELP * FLUXVEC(O2)
                   ENDDO
                   L_SS_LAYERSOURCE = SUM
                   L_TR_CUMSOURCE = 
     &           + L_T_UTDN_USERM(UT,UM,Q) *   SS_CUMSOURCE_DN(V,O1,NC)
     &           +   T_UTDN_USERM(UT,UM)   * L_SS_CUMSOURCE(Q,V,O1)
                   L_FINAL_SOURCE = L_TR_CUMSOURCE + L_SS_LAYERSOURCE
                   L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                   COLUMNWF_SS(Q,UTA,V,O1,DNIDX) = L_SSCORRECTION
                  ENDDO
                 ENDDO
                ENDDO
              ENDDO
            ENDDO

C  Ongrid output
C  -------------

C  just set to the cumulative source term 

          ELSE

            DO V = 1, N_GEOMETRIES
              DO Q = 1, NV_PARAMETERS
                DO O1 = 1, NSTOKES
                  L_FINAL_SOURCE = L_SS_CUMSOURCE(Q,V,O1)
                  L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                  COLUMNWF_SS(Q,UTA,V,O1,DNIDX) = L_SSCORRECTION
                ENDDO
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

      SUBROUTINE VLIDORTSS_L_ZMATRICES
     I  ( NSTOKES, NV, NV_PARAMETERS, DO_SCATMAT_VARIATION,
     I    N_GEOMETRIES, NGREEKMOMS, 
     I    NBEAMS, N_USER_STREAMS,  N_USER_RELAZMS,
     I    LAYER_MAXMOMENTS, 
     I    GREEKMAT_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT,
     I    DO_SSCORR_TRUNCATION, SSFDEL, L_SSFDEL,
     I    VSIGN, VZA_OFFSETS,
     I    CTHETA, STHETA, CALPHA, SALPHA, CPHI, PHI,   
     O    L_ZMAT )

C  include file of dimensions and numbers
C     Cannot use bookkeeping vars file

      INCLUDE '../includes/VLIDORT.PARS'

C  input
C  -----

C  control

      INTEGER          NSTOKES, NV, NV_PARAMETERS
      LOGICAL          DO_SCATMAT_VARIATION (MAXLAYERS,MAX_ATMOSWFS)
      INTEGER          N_GEOMETRIES, NGREEKMOMS
      INTEGER          NBEAMS, N_USER_STREAMS,  N_USER_RELAZMS

C  +/- sign, offsets

      DOUBLE PRECISION VSIGN 
      INTEGER          VZA_OFFSETS(MAXBEAMS,MAX_USER_STREAMS)

C  scattering input information

      INTEGER          LAYER_MAXMOMENTS(MAXLAYERS)
      DOUBLE PRECISION GREEKMAT_TOTAL_INPUT
     &      ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION L_GREEKMAT_TOTAL_INPUT
     &   ( MAX_ATMOSWFS, 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )

C  Truncation control

      LOGICAL          DO_SSCORR_TRUNCATION
      DOUBLE PRECISION SSFDEL   ( MAXLAYERS )
      DOUBLE PRECISION L_SSFDEL ( MAXLAYERS,MAX_ATMOSWFS )

C  zenith angle cosines/sines, azimuth angle cosines

      DOUBLE PRECISION CTHETA (MAXLAYERS,MAXBEAMS)
      DOUBLE PRECISION STHETA (MAXLAYERS,MAXBEAMS)
      DOUBLE PRECISION CALPHA (MAX_USER_STREAMS)
      DOUBLE PRECISION SALPHA (MAX_USER_STREAMS)
      DOUBLE PRECISION CPHI   (MAX_USER_RELAZMS)

C  azimuth angles. Added V2.4R. for PHI > 180 case

      DOUBLE PRECISION PHI    (MAX_USER_RELAZMS)

C  output
C  ------

C  Z-Matrix linearized

      DOUBLE PRECISION L_ZMAT 
     &            ( MAX_ATMOSWFS, MAX_GEOMETRIES, MAXLAYERS, 4, 4 )

C  Local
C  -----

C  Rotation Angle cosines/sines

      DOUBLE PRECISION C1 (MAX_GEOMETRIES)
      DOUBLE PRECISION S1 (MAX_GEOMETRIES)
      DOUBLE PRECISION C2 (MAX_GEOMETRIES)
      DOUBLE PRECISION S2 (MAX_GEOMETRIES)

C  local F-matrix

      DOUBLE PRECISION L_FMAT ( MAX_ATMOSWFS, MAX_GEOMETRIES, 6 )

C  help variables

      INTEGER          IB, UM, IA, V, Q, GREEKMAT_INDEX(6)
      DOUBLE PRECISION COSSCAT, SINSCAT, HELP_SINSCAT
      DOUBLE PRECISION CSIG1, CSIG2, SSIG1, SSIG2, CSIG1_2, CSIG2_2

      DOUBLE PRECISION UUU(MAX_GEOMETRIES)
      DOUBLE PRECISION P00(MAX_GEOMETRIES,2)
      DOUBLE PRECISION P02(MAX_GEOMETRIES,2)

      INTEGER          K, L, LNEW, LOLD, ITMP
      INTEGER          INDEX_11, INDEX_12, INDEX_34
      INTEGER          INDEX_22, INDEX_33, INDEX_44
      DOUBLE PRECISION DL, QROOT6, FAC1, FAC2, SQL4, SQL41, TMP1, TMP2
      DOUBLE PRECISION  HELP1, HELP2, FDNL1, DNL1, FACT
      DOUBLE PRECISION  GK11, L_GK11(MAX_ATMOSWFS)
      DOUBLE PRECISION  GK12, L_GK12(MAX_ATMOSWFS)

C  Sunlight

      LOGICAL          SUNLIGHT
      PARAMETER        ( SUNLIGHT = .TRUE. )

C  indexing key

C      GREEKMAT_INDEX(1) = 1  ---> INDEX_11
C      GREEKMAT_INDEX(2) = 6  ---> INDEX_22
C      GREEKMAT_INDEX(3) = 2  ---> INDEX_12
C      GREEKMAT_INDEX(4) = 11 ---> INDEX_33
C      GREEKMAT_INDEX(5) = 12 ---> INDEX_34
C      GREEKMAT_INDEX(6) = 16 ---> INDEX_44

      GREEKMAT_INDEX(1) = 1
      GREEKMAT_INDEX(2) = 6
      GREEKMAT_INDEX(3) = 2
      GREEKMAT_INDEX(4) = 11
      GREEKMAT_INDEX(5) = 12
      GREEKMAT_INDEX(6) = 16

      INDEX_11 = 1
      INDEX_12 = 2
      INDEX_22 = 3
      INDEX_33 = 4
      INDEX_34 = 5
      INDEX_44 = 6

      SQL41 = ZERO

C  Geometrical quantities
C  ----------------------

      DO IB = 1, NBEAMS
        DO UM = 1, N_USER_STREAMS
          DO IA = 1, N_USER_RELAZMS
            V = VZA_OFFSETS(IB,UM) + IA

C  cosine scatter angle (this is valid only for non-refracting atmosphere)
C  VSIGN = -1 for upwelling, +1 for downwelling

            COSSCAT = VSIGN * CTHETA(NV,IB) * CALPHA(UM) +
     &                        STHETA(NV,IB) * SALPHA(UM) * CPHI(IA)
            UUU(V)  = COSSCAT

C  Cosine Sigma 1 and 2. H/VdM, Eqs. (99)-(101)

C  a. safety
C    Watch for sin^2(scatter angle) less than zero (machine precision)
C    R. Spurr, 16 January 2006, RT SOLUTIONS Inc.

            HELP_SINSCAT = ( ONE - COSSCAT * COSSCAT )
            IF ( HELP_SINSCAT.LE.ZERO ) THEN
              SINSCAT = 1.0D-12
            ELSE
              SINSCAT = DSQRT ( HELP_SINSCAT )
            ENDIF

C  b. necessary limit analyses - Hovenier limits.
C     R. Spurr and V. Natraj, 17 January 2006

            IF ( DABS(SINSCAT) .LE. 1.0D-12 ) THEN
              CSIG1 = ZERO
              CSIG2 = ZERO
            ELSE
             IF ( STHETA(1,IB) .EQ. ZERO ) THEN
              CSIG1 = -  CPHI(IA)
             ELSE
              CSIG1   = ( - VSIGN * CALPHA(UM) + CTHETA(NV,IB)*COSSCAT )
     &                  / SINSCAT / STHETA(NV,IB)
             ENDIF
             IF ( SALPHA(UM) .EQ. ZERO ) THEN
               CSIG2 = -  CPHI(IA)
             ELSE
               CSIG2   = ( - CTHETA(NV,IB) + VSIGN*CALPHA(UM)*COSSCAT )
     &                    / SINSCAT / SALPHA(UM)
             ENDIF
            ENDIF
c    write(7,'(4i5,1p4e15.7)')V,IB,UM,IA,CSIG1,CSIG2,sinscat

            IF ( CSIG2 .GT. ONE  ) CSIG2 = ONE
            IF ( CSIG2 .LT. -ONE ) CSIG2 = -ONE

C  output, H/VdM, Eqs. (89)-(94)
C    Rotation sines and cosines

            CSIG1_2 = TWO * CSIG1
            CSIG2_2 = TWO * CSIG2
            IF ( DABS(CSIG1-ONE).LT.1.0D-12)THEN
             SSIG1 = ZERO
            ELSE
             SSIG1 = DSQRT ( 1.0D0 - CSIG1 * CSIG1 )
            ENDIF
            IF ( DABS(CSIG2-ONE).LT.1.0D-12)THEN
             SSIG2 = ZERO
            ELSE
             SSIG2 = DSQRT ( 1.0D0 - CSIG2 * CSIG2 )
            ENDIF

C  For relazm in [180,360), need sign reversal for S1 and S2
C  See H/VdM, Eqs. 94-95. V. Natraj and R. Spurr, 01 May 2009.

            C1(V) = CSIG1_2 * CSIG1 - ONE
            C2(V) = CSIG2_2 * CSIG2 - ONE

            IF (PHI(IA) .LE. 180.D0) THEN
              S1(V) = CSIG1_2 * SSIG1
              S2(V) = CSIG2_2 * SSIG2
            ELSE
              S1(V) = -CSIG1_2 * SSIG1
              S2(V) = -CSIG2_2 * SSIG2
            ENDIF

C  End geometry loops

          ENDDO
        ENDDO
      ENDDO

C  F-matrices
C  ----------

      QROOT6 = -0.25D0 * DSQRT(6.0D0)

C initialise F-matrix

      DO V = 1, N_GEOMETRIES
        DO K = 1, 6
          DO Q = 1, NV_PARAMETERS
            L_FMAT(Q,V,K) = ZERO
          END DO
        END DO
      END DO

C  Start loop over the coefficient index l
C  first update generalized spherical functions, then calculate coefs.
C  lold and lnew are pointer-like indices used in recurrence 

      LNEW = 1
      LOLD = 2

      DO L = 0, NGREEKMOMS

C  Local Greek matrix entries
C  Set the local Greek matrix elements that you need = 11 and 12.
c   44 and 34 are not required with natural sunlight (default here)
C   Linearized 44 and 34 Have not been CODED..............!!!!!!!!!!

        DO Q = 1, NV_PARAMETERS
         IF ( DO_SCATMAT_VARIATION(NV,Q) ) THEN
          IF ( DO_SSCORR_TRUNCATION ) THEN
           DNL1  = DBLE(2*L + 1 )
           FDNL1 = SSFDEL(NV) * DNL1
           FACT  = ONE - SSFDEL(NV)
           GK11 = ( GREEKMAT_TOTAL_INPUT(L,NV,1) - FDNL1 ) / FACT
           GK12 =   GREEKMAT_TOTAL_INPUT(L,NV,2) / FACT
c           GK44 = ( GREEKMAT_TOTAL_INPUT(L,NV,16) - FDNL1 ) / FACT
c           GK34 =   GREEKMAT_TOTAL_INPUT(L,NV,12) / FACT
           HELP1 = L_GREEKMAT_TOTAL_INPUT(Q,L,NV,1) * 
     &                 GREEKMAT_TOTAL_INPUT(L,NV,1)
           HELP2 = ( GK11 - DNL1 ) * L_SSFDEL(NV,Q)
           L_GK11(Q) = ( HELP1 + HELP2 ) / FACT
           HELP1 = L_GREEKMAT_TOTAL_INPUT(Q,L,NV,2) * 
     &                 GREEKMAT_TOTAL_INPUT(L,NV,2)
           HELP2 = GK12 * L_SSFDEL(NV,Q)
           L_GK12(Q) = ( HELP1 + HELP2 ) / FACT
          ELSE
           L_GK11(Q) = L_GREEKMAT_TOTAL_INPUT(Q,L,NV,1) *
     &                   GREEKMAT_TOTAL_INPUT(L,NV,1)
           L_GK12(Q) = L_GREEKMAT_TOTAL_INPUT(Q,L,NV,2) *
     &                     GREEKMAT_TOTAL_INPUT(L,NV,2)
c           L_GK11(Q) = L_GREEKMAT_TOTAL_INPUT(Q,L,NV,1)
c     &                     GREEKMAT_TOTAL_INPUT(L,NV,1)
c           L_GK12(Q) = L_GREEKMAT_TOTAL_INPUT(Q,L,NV,2)
c     &                     GREEKMAT_TOTAL_INPUT(L,NV,2)
          ENDIF
         ENDIF 
        ENDDO

        DL   = DBLE(L)

        IF ( L .EQ. 0 ) THEN

C  Adding paper Eqs. (76) and (77) with m=0

          DO V = 1, N_GEOMETRIES
            P00(V,LOLD) = ONE
            P00(V,LNEW) = ZERO
            P02(V,LOLD) = ZERO
            P02(V,LNEW) = ZERO
          END DO

        ELSE

          FAC1 = (TWO*DL-ONE)/DL
          FAC2 = (DL-ONE)/DL

C Adding paper Eq. (81) with m=0

          DO V = 1, N_GEOMETRIES
            P00(V,LOLD) = FAC1*UUU(V)*P00(V,LNEW) - FAC2*P00(V,LOLD)
          END DO

        END IF

        IF ( L .EQ. 2 ) THEN

! Adding paper Eq. (78)  
! sql4 contains the factor dsqrt((l+1)*(l+1)-4) needed in
! the recurrence Eqs. (81) and (82)

          DO V = 1, N_GEOMETRIES
            P02(V,LOLD) = QROOT6*(ONE-UUU(V)*UUU(V))
            P02(V,LNEW) = ZERO
          END DO
          SQL41 = ZERO

        ELSE IF ( L .GT. 2) THEN

! Adding paper Eq. (82) with m=0

          SQL4  = SQL41
          SQL41 = DSQRT(DL*DL-FOUR)
          TMP1  = (TWO*DL-ONE)/SQL41
          TMP2  = SQL4/SQL41
          DO V = 1, N_GEOMETRIES
            P02(V,LOLD) = TMP1*UUU(V)*P02(V,LNEW) - TMP2*P02(V,LOLD)
          END DO

        END IF

! Switch indices so that lnew indicates the function with
! the present index value l, this mechanism prevents swapping
! of entire arrays.

        ITMP = LNEW
        LNEW = LOLD
        LOLD = ITMP

! Now add the l-th term to the scattering matrix.
! See de Haan et al. (1987) Eqs. (68)-(73).
! Remember for Mie scattering : F11 = F22 and F33 = F44

        DO Q = 1, NV_PARAMETERS
         IF ( DO_SCATMAT_VARIATION(NV,Q) ) THEN
          IF ( L.LE.LAYER_MAXMOMENTS(NV) ) THEN
           DO V = 1, N_GEOMETRIES
            L_FMAT(Q,V,INDEX_11) =
     &            L_FMAT(Q,V,INDEX_11) + L_GK11(Q) * P00(V,LNEW)
            L_FMAT(Q,V,INDEX_12) =
     &            L_FMAT(Q,V,INDEX_12) + L_GK12(Q) * P02(V,LNEW)
C  Previous code
c            L_FMAT(Q,V,INDEX_11) = L_FMAT(Q,V,INDEX_11) +
c     &       L_GREEKMAT_TOTAL_INPUT(Q,L,NV,GREEKMAT_INDEX(1))*P00(V,LNEW)
c            L_FMAT(Q,V,INDEX_12) = L_FMAT(Q,V,INDEX_12) +
c     &       L_GREEKMAT_TOTAL_INPUT(Q,L,NV,GREEKMAT_INDEX(3))*P02(V,LNEW)
c            L_FMAT(Q,V,INDEX_44) = L_FMAT(Q,V,INDEX_44) +
c     &       L_GREEKMAT_TOTAL_INPUT(Q,L,NV,GREEKMAT_INDEX(6))*P00(V,LNEW)
c            L_FMAT(Q,V,INDEX_34) = L_FMAT(Q,V,INDEX_34) +
c     &       L_GREEKMAT_TOTAL_INPUT(Q,L,NV,GREEKMAT_INDEX(5))*P02(V,LNEW)
           ENDDO
          ENDIF
         ENDIF
        END DO

C  end moment loop

      END DO

C  remaining symmetries for Mie particles

c      DO V = 1, N_GEOMETRIES
c        DO Q = 1, NV_PARAMETERS
c          L_FMAT(Q,V,INDEX_22) = L_FMAT(Q,V,INDEX_11)
c          L_FMAT(Q,V,INDEX_33) = L_FMAT(Q,V,INDEX_44)
c        END DO
c      END DO

C  Create Z matrices
C  -----------------

C  Copy for Stokes = 1
    
      IF ( NSTOKES .EQ. 1 ) THEN

       DO Q = 1, NV_PARAMETERS
         DO V = 1, N_GEOMETRIES
           L_ZMAT(Q,V,NV,1,1) = L_FMAT(Q,V,INDEX_11)
         ENDDO
       ENDDO

C  For polarized case, sunlight only !!!!!!!
C   For sunlight, only need the first Column of the Z-matrix

      ELSE IF ( NSTOKES .GT. 1 ) THEN

       IF ( SUNLIGHT ) THEN
         DO V = 1, N_GEOMETRIES
           DO Q = 1, NV_PARAMETERS
             L_ZMAT(Q,V,NV,1,1) =   L_FMAT(Q,V,INDEX_11)
             L_ZMAT(Q,V,NV,2,1) =  -L_FMAT(Q,V,INDEX_12) * C2(V)
             L_ZMAT(Q,V,NV,3,1) =   L_FMAT(Q,V,INDEX_12) * S2(V)
             L_ZMAT(Q,V,NV,4,1) =   ZERO
           ENDDO
          ENDDO
       ENDIF

      ENDIF

C  finish

      RETURN
      END

c

      SUBROUTINE VLIDORT_L_SSCORR_OUTGOING
     &        ( SSFLUX, FAIL, MESSAGE)

C  Single scatter exact calculation for the outgoing LOS
C   This one with optional linearizations
C         - NEW for Version 2.2

C   Programmed by R. Spurr, RT Solutions Inc.
C    First Draft, January 31st 2007
C   Validated against TOMRAD, 29 March 2007.
C   Partial layer output added September 2007

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include file of input variables
C  Include file of bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup and multiplier variables (input)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  Include file of single scatter result variables

      INCLUDE '../includes/VLIDORT_SINGSCAT.VARS'

C  include files of linearized inputs/setups/solution/multiplier variables

      INCLUDE '../includes/VLIDORT_L_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'

C  Include file of results (linearized single scatter weighting funcs)

      INCLUDE '../includes/VLIDORT_L_SINGSCAT.VARS'

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

      integer          ntraverse_ut(max_partlayers)
      double precision sunpaths_ut (max_partlayers,maxlayers)
      double precision radii_ut    (max_partlayers)
      double precision alpha_ut    (max_partlayers)

C  Other (incidental) geometrical output

      double precision lospaths(maxlayers)
      double precision lospaths_ut_up(max_partlayers)
      double precision lospaths_ut_dn(max_partlayers)
      double precision theta_all  (0:maxlayers)
      double precision phi_all    (0:maxlayers)
      double precision cosscat_up (0:maxlayers)
      double precision cosscat_dn (0:maxlayers)

C  Extinction and var tms

      double precision extinction   (maxlayers)
      double precision L_extinction (maxlayers,max_atmoswfs)
      double precision var_tms      (maxlayers,max_atmoswfs)

C  Partial layer heights

      double precision height_grid_ut(max_partlayers)

C  Bookkeeping

      integer          NTYPE_VARY
      logical          DO_RTSOL_VARY(maxlayers)
      integer          NPARAMS_VARY(maxlayers) 
      integer          KINDEX(maxlayers)

C  local variables
C  ---------------

C  Indices

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER          UT, UTA, UM, IA, NC, IB, V, O1, NM1
      integer          k, k_parameters, ks, q

C  help variables (double precision)

      DOUBLE PRECISION FINAL_SOURCE, HELP, SS_CUMSOURCE
      DOUBLE PRECISION SS_LAYERSOURCE, SSCORRECTION, XT
      DOUBLE PRECISION CTHETA, STHETA, CALPHA, SALPHA, CPHI, CSA
      DOUBLE PRECISION VSIGN, ZMAT_LOCAL(4), FMAT_LOCAL(6), DNM1
      DOUBLE PRECISION L_ZMAT_LOCAL ( 4, MAX_ATMOSWFS )
      DOUBLE PRECISION UVAR, AVAR, FT1, FT2, LSS, TRANS
      DOUBLE PRECISION L_FINAL_SOURCE, L_SSCORRECTION

C  Sunlight parameters
C    - assumes the Flux vector is (F,0,0,0)
C   - drop this variables when further testing has been done

      LOGICAL           SUNLIGHT
      PARAMETER         ( SUNLIGHT = .TRUE. )
      INTEGER           NPOLAR

C  Set up operations
C  -----------------

C  Set up partials flag

      do_fine = .true.
      do_partials = ( n_partlayers .gt. 0 )

C  Set npolar. For Natural light, there are only 3 Stokes parameters
C    ( no circular polarization for single scattering)

      IF (SUNLIGHT) NPOLAR = MIN(NSTOKES,3)
      IF (.NOT.SUNLIGHT) NPOLAR = NSTOKES

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
           IF ( DO_SCATMAT_VARIATION(K,Q) ) THEN
            UVAR = L_OMEGA_TOTAL_INPUT(Q,K)
            AVAR = UVAR + L_GREEKMAT_TOTAL_INPUT(Q,NM1,K,1)
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

        NM1  = NGREEK_MOMENTS_INPUT
        DNM1 = DBLE(2*NM1+1)
        DO N = 1, NLAYERS
          SSFDEL(N) = GREEKMAT_TOTAL_INPUT(NM1,N,1) / DNM1
          TMS(N) = TMS(N) * ( ONE - SSFDEL(N) )
        ENDDO

C  Linearization, Only change if SCATMAT variation is set

        IF ( DO_ATMOS_LINEARIZATION ) THEN
          DO K = 1, NLAYERS
            IF ( LAYER_VARY_FLAG(K) ) THEN
              K_PARAMETERS = LAYER_VARY_NUMBER(K)
              DO Q = 1, K_PARAMETERS
                IF ( DO_SCATMAT_VARIATION(K,Q) ) THEN
                  FT1 = ONE / ( ONE- SSFDEL(K) )
                  L_SSFDEL(K,Q) = GREEKMAT_TOTAL_INPUT(NM1,K,1) *
     *                    L_GREEKMAT_TOTAL_INPUT(Q,NM1,K,1) / DNM1
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
        do ut = 1, n_partlayers
          n = PARTLAYERS_layeridx(ut)
          xt = deltau_vert(n) - PARTAU_VERT(UT)
          height_grid_ut(ut) = height_grid(n) + xt / extinction(n)
        enddo
      endif

C  Compute and store all geometry-dependent stuff
C  ==============================================

C  Start the main loop over all solar and viewing geometries
C   Old code, not adjusted
c      DO IB = 1, NBEAMS
c        THETA_BOA = SZANGLES(IB)
c        DO UM = 1, N_USER_STREAMS
c          ALPHA_BOA = USER_VZANGLES_INPUT(UM)
c          DO IA = 1, N_USER_RELAZMS
c            PHI_BOA = USER_RELAZMS(IA)
c            V = VZA_OFFSETS(IB,UM) + IA

C  Start the main loop over all solar and viewing geometries

      DO UM = 1, N_USER_STREAMS
        ALPHA_BOA = USER_VZANGLES_ADJUST(UM)
        DO IB = 1, NBEAMS
          DO IA = 1, N_USER_RELAZMS
            THETA_BOA = SZANGLES_ADJUST(UM,IB,IA)
            PHI_BOA   = USER_RELAZMS_ADJUST(UM,IB,IA)
            V = VZA_OFFSETS(IB,UM) + IA
    
C  Call to geometry

            call outgoing_sphergeom_fine
     i  ( maxlayers, maxfinelayers, max_partlayers,
     i    do_fine, do_partials, nlayers, nfinelayers,
     i    n_partlayers, PARTLAYERS_layeridx,
     i    height_grid, height_grid_ut, earth_radius,
     i    alpha_boa, theta_boa, phi_boa,
     o    sunpaths,      radii,      ntraverse,      alpha_all, 
     o    sunpaths_fine, radii_fine, ntraverse_fine, alpha_fine,
     o    sunpaths_ut,   radii_ut,   ntraverse_ut,   alpha_ut,
     o    PARTLAYERS_layerfineidx, lospaths,
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
     i             do_partials,
     i     do_profile_linearization, layer_vary_flag, layer_vary_number,
     i     do_column_linearization, n_totalcolumn_wfs,
     i             extinction, l_extinction,
     i             n_partlayers, PARTLAYERS_layeridx,
     i             PARTLAYERS_layerfineidx,
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
c              if (v.eq.1)write(35,*)'lin',v
c              do ut = 1, n_partlayers
c                if (v.eq.1)write(35,'(1p2e18.10)')
c     &           up_lostrans_ut(ut,v),up_multipliers_ut(ut,v)
c              enddo

c              write(35,*)v
c              do n = 1, nlayers
c                write(35,'(1p2e18.10)')
c     &       up_lostrans(n,v),up_multipliers(n,v)
c              enddo
c              if ( v.eq.8)pause

c              if ( do_atmos_linearization ) then
c              do n = 1, nlayers
c                write(36,'(2i4,1p20e18.10)')v,n,
c     &       up_lostrans(n,v),up_multipliers(n,v),
c     &l_up_lostrans(n,1,v),(l_up_multipliers(n,k,1,v),k=1,nlayers)
c              enddo
c               else
c              do n = 1, nlayers
c                write(35,'(2i4,1p2e18.10)')v,n,
c     &       up_lostrans(n,v),up_multipliers(n,v)
c                enddo
c               endif

C  Start layer loop

              DO N = NLAYERS, 1, -1

C  Trigonometry

                CTHETA = DCOS(THETA_ALL(N))
                STHETA = DSIN(THETA_ALL(N))
                CALPHA = DCOS(ALPHA_ALL(N))
                SALPHA = DSIN(ALPHA_ALL(N))
                CPHI   = DCOS(PHI_ALL(N))
                CSA    = COSSCAT_UP(N)
                VSIGN  = -1.0d0

C  Get the scattering matrix and its linearization

                CALL SSCORR_OUTGOING_L_ZMATRIX
     I          ( DO_ATMOS_LINEARIZATION, DO_SCATMAT_VARIATION,
     I            N, NSTOKES, 
     I            NGREEK_MOMENTS_INPUT, LAYER_MAXMOMENTS(N),
     I            GREEKMAT_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT,
     I            DO_SSCORR_TRUNCATION, SSFDEL, L_SSFDEL,
     I            LAYER_VARY_FLAG(N), LAYER_VARY_NUMBER(N),
     I            CTHETA, STHETA, CALPHA, SALPHA, CPHI,
     I            PHI_ALL(N), CSA, VSIGN, 
     O            ZMAT_LOCAL, L_ZMAT_LOCAL, FMAT_LOCAL )

C  Phase matrix + lienarized (multiplied by TMS factor). Save them.
C     Sunlight only, the first column of the matrix

                IF ( STERM_LAYERMASK_UP(N) ) THEN
                 DO O1 = 1, NSTOKES
                  ZMAT_UP(V,N,O1,1) = ZMAT_LOCAL(O1) * TMS(N)
                 ENDDO
                 IF ( DO_ATMOS_LINEARIZATION ) THEN
                  IF ( LAYER_VARY_FLAG(N) ) THEN
                   DO Q = 1, LAYER_VARY_NUMBER(N)
                    DO O1 = 1, NSTOKES
                     L_ZMAT_UP(Q,V,N,O1,1) = 
     &                        ZMAT_UP(V,N,O1,1) * VAR_TMS(N,Q) +
     &                      L_ZMAT_LOCAL(O1,Q)  *     TMS(N)
                    ENDDO
                   ENDDO
                  ENDIF
                 ENDIF
                ENDIF

C  Finish the layer loop

              ENDDO

C  End upwelling clause

            ENDIF

C  Downwelling calculation
C  -----------------------

            IF ( DO_DNWELLING ) THEN

C  Multipliers, transmittances + linearizations

              call l_outgoing_integration_dn
     i           ( nlayers, nfinelayers,
     i             do_partials, 
     i     do_profile_linearization, layer_vary_flag, layer_vary_number,
     i     do_column_linearization, n_totalcolumn_wfs,
     i             extinction, l_extinction,
     i             n_partlayers, PARTLAYERS_layeridx,
     i             PARTLAYERS_layerfineidx,
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

C  Layer loop

              DO N = 1, NLAYERS

C  Trigonometry

                CTHETA = DCOS(THETA_ALL(N-1))
                STHETA = DSIN(THETA_ALL(N-1))
                CALPHA = DCOS(ALPHA_ALL(N-1))
                SALPHA = DSIN(ALPHA_ALL(N-1))
                CPHI   = DCOS(PHI_ALL(N-1))
                CSA    = COSSCAT_DN(N-1)
                VSIGN  = +1.0d0

C  Get the scattering matrix and its linearization

                CALL SSCORR_OUTGOING_L_ZMATRIX
     I          ( DO_ATMOS_LINEARIZATION, DO_SCATMAT_VARIATION,
     I            N, NSTOKES, 
     I            NGREEK_MOMENTS_INPUT, LAYER_MAXMOMENTS(N),
     I            GREEKMAT_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT,
     I            DO_SSCORR_TRUNCATION, SSFDEL, L_SSFDEL,
     I            LAYER_VARY_FLAG(N), LAYER_VARY_NUMBER(N),
     I            CTHETA, STHETA, CALPHA, SALPHA, CPHI,
     I            PHI_ALL(N-1), CSA, VSIGN, 
     O            ZMAT_LOCAL, L_ZMAT_LOCAL, FMAT_LOCAL )

C  Phase matrix + linearization (multiplied by TMS factor). Save them.
C     Sunlight only, the first column of the matrix

                IF ( STERM_LAYERMASK_DN(N) ) THEN
                 DO O1 = 1, NSTOKES
                  ZMAT_DN(V,N,O1,1) = ZMAT_LOCAL(O1) * TMS(N)
                 ENDDO
                 IF ( DO_ATMOS_LINEARIZATION ) THEN
                  IF ( LAYER_VARY_FLAG(N) ) THEN
                   DO Q = 1, LAYER_VARY_NUMBER(N)
                    DO O1 = 1, NSTOKES
                     L_ZMAT_DN(Q,V,N,O1,1) = 
     &                        ZMAT_DN(V,N,O1,1)  * VAR_TMS(N,Q) +
     &                      L_ZMAT_LOCAL(O1,Q)   *     TMS(N)
                    ENDDO
                   ENDDO
                  ENDIF
                 ENDIF
                ENDIF

C  Finish the layer loop

              ENDDO

C  End Downwelling clause

            ENDIF

C   Finish geometry loops

          ENDDO
        ENDDO
      ENDDO

C  Recurrence relation for the UPWELLING Stokes vector
C  ===================================================

      IF ( DO_UPWELLING ) THEN

C  initialize cumulative source term

        NC =  0
        DO V = 1, N_GEOMETRIES
          DO O1 = 1, NSTOKES
            SS_CUMSOURCE_UP(V,O1,NC) = ZERO
          ENDDO
        ENDDO

C  initialise optical depth loop

        NSTART = NLAYERS
        NUT_PREV = NSTART + 1

C  Main loop over all output optical depths

        DO UTA = N_USER_LEVELS, 1, -1

C  Layer index for given optical depth

          NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
          NUT    = NLEVEL + 1

C  Cumulative single scatter source terms :
C      For loop over layers working upwards to level NUT,
C      Get layer source terms = Exact Z-matrix * Multiplier
C      sunlight case, no circular polarization

          DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            DO V = 1, N_GEOMETRIES
              DO O1 = 1, NPOLAR
                HELP = ZMAT_UP(V,N,O1,1) * FLUXVEC(1)
                SS_LAYERSOURCE = HELP* UP_MULTIPLIERS(N,V)
                SS_CUMSOURCE_UP(V,O1,NC) = SS_LAYERSOURCE +
     &             UP_LOSTRANS(N,V)*SS_CUMSOURCE_UP(V,O1,NC-1)
              ENDDO
            ENDDO
          ENDDO

C  sunlight case, no circular polarization
C  Offgrid output-------
C    Add additional partial layer source term = Exact Phase Func * Multiplier
C    Set final cumulative source and single scatter Stokes vector
C  Ongrid output--------
C    Set final cumulative source and single scatter Stokes vector

          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)
            DO V = 1, N_GEOMETRIES
              DO O1 = 1, NPOLAR
                HELP = ZMAT_UP(V,N,O1,1) * FLUXVEC(1)
                SS_LAYERSOURCE = HELP * UP_MULTIPLIERS_UT(UT,V)
                SS_CUMSOURCE   = SS_CUMSOURCE_UP(V,O1,NC)
                TRANS = UP_LOSTRANS_UT(UT,V)
                FINAL_SOURCE   = TRANS*SS_CUMSOURCE + SS_LAYERSOURCE
                SSCORRECTION   = SSFLUX * FINAL_SOURCE
                STOKES_SS(UTA,V,O1,UPIDX) = SSCORRECTION
              ENDDO
            ENDDO
          ELSE
            DO V = 1, N_GEOMETRIES
              DO O1 = 1, NPOLAR
                FINAL_SOURCE = SS_CUMSOURCE_UP(V,O1,NC)
                SSCORRECTION = SSFLUX * FINAL_SOURCE
                STOKES_SS(UTA,V,O1,UPIDX) = SSCORRECTION
              ENDDO
            ENDDO
          ENDIF

C  Check for updating the recursion 

          IF ( NUT. NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT

C  end optical depth loop and Upwelling clause

        ENDDO
      ENDIF

C  Recurrence relation for the DOWNWELLING Stokes vector
C  =====================================================

      IF ( DO_DNWELLING ) THEN

C  initialize cumulative source term

        NC =  0
        DO V = 1, N_GEOMETRIES
          DO O1 = 1, NSTOKES
            SS_CUMSOURCE_DN(V,O1,NC) = ZERO
          ENDDO
        ENDDO

C  initialise optical depth loop

        NSTART = 1
        NUT_PREV = NSTART - 1

C  Main loop over all output optical depths

        DO UTA = 1, N_USER_LEVELS

C  Layer index for given optical depth

          NLEVEL = UTAU_LEVEL_MASK_DN(UTA)
          NUT = NLEVEL

C  Cumulative single scatter source terms :
C      For loop over layers working downwards to NUT,
C      Get layer source terms = Exact Z-matrix * Multiplier
C      sunlight case only

          DO N = NSTART, NUT
            NC = N
            DO V = 1, N_GEOMETRIES
              DO O1 = 1, NPOLAR
                HELP =  ZMAT_DN(V,N,O1,1) * FLUXVEC(1)
                SS_LAYERSOURCE = HELP* DN_MULTIPLIERS(N,V)
                SS_CUMSOURCE_DN(V,O1,NC) = SS_LAYERSOURCE +
     &             DN_LOSTRANS(N,V)*SS_CUMSOURCE_DN(V,O1,NC-1)
              ENDDO
            ENDDO
          ENDDO

C    sunlight case
C  Offgrid output :
C    add additional partial layer source term = Exact Z-matrix * Multiplier
C    Set final cumulative source and Correct Stokes vector
C  Ongrid output :
C    Set final cumulative source and Correct Stokes vector

          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)
            DO V = 1, N_GEOMETRIES
              DO O1 = 1, NPOLAR
                HELP = ZMAT_DN(V,N,O1,1) * FLUXVEC(1)
                SS_LAYERSOURCE = HELP * DN_MULTIPLIERS_UT(UT,V)
                SS_CUMSOURCE   = SS_CUMSOURCE_DN(V,O1,NC)
                TRANS = DN_LOSTRANS_UT(UT,V)
                FINAL_SOURCE   = TRANS*SS_CUMSOURCE + SS_LAYERSOURCE
                SSCORRECTION   = SSFLUX * FINAL_SOURCE
                STOKES_SS(UTA,V,O1,DNIDX) = SSCORRECTION
              ENDDO
            ENDDO
          ELSE
            DO V = 1, N_GEOMETRIES
              DO O1 = 1, NPOLAR
                FINAL_SOURCE = SS_CUMSOURCE_DN(V,O1,NC)
                SSCORRECTION = SSFLUX * FINAL_SOURCE
                STOKES_SS(UTA,V,O1,DNIDX) = SSCORRECTION
              ENDDO
            ENDDO
          ENDIF

C  Check for updating the recursion 

          IF ( NUT. NE. NUT_PREV ) NSTART = NUT + 1
          NUT_PREV = NUT

C  end optical depth loop and Downwelling clause

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
           DO O1 = 1, NSTOKES
            DO Q = 1, K_PARAMETERS
             L_SS_CUMSOURCE(Q,V,O1) = ZERO
            ENDDO
           ENDDO
          ENDDO

C  initialise optical depth loop

          NSTART = NLAYERS
          NUT_PREV = NSTART + 1

C  Main loop over all output optical depths
C  ----------------------------------------

          DO UTA = N_USER_LEVELS, 1, -1

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

             IF ( N .EQ. K .OR. KS.EQ.0 ) THEN
              DO V = 1, N_GEOMETRIES
               DO O1 = 1, NPOLAR
                DO Q = 1, K_PARAMETERS
                 LSS = ZMAT_UP(V,N,O1,1)   * L_UP_MULTIPLIERS(N,KS,Q,V)
     &             + L_ZMAT_UP(Q,V,N,O1,1) *   UP_MULTIPLIERS(N,V)
                 L_SS_CUMSOURCE(Q,V,O1) = FLUXVEC(1) * LSS
     &            +  UP_LOSTRANS(N,V)    * L_SS_CUMSOURCE(Q,V,O1)
     &            + L_UP_LOSTRANS(N,Q,V) *   SS_CUMSOURCE_UP(V,O1,NC-1)
                ENDDO
               ENDDO
              ENDDO
             ELSE
              DO V = 1, N_GEOMETRIES
               DO O1 = 1, NPOLAR
                DO Q = 1, K_PARAMETERS
                 LSS = ZMAT_UP(V,N,O1,1)   * L_UP_MULTIPLIERS(N,KS,Q,V)
                 L_SS_CUMSOURCE(Q,V,O1) = FLUXVEC(1) * LSS
     &               +  UP_LOSTRANS(N,V) * L_SS_CUMSOURCE(Q,V,O1)
                ENDDO
               ENDDO
              ENDDO
             ENDIF

C  End layer loop

            ENDDO

C  Offgrid output----------
C  Set final cumulative source and Single scatter Weighting function

           IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)

C  If N = K (the layer that is varying)
C    add Linearization of additional partial layer source term =
C        L_Exact_Scat(n) * Multiplier  +  Exact_Scat * L_Multiplier(n)
C  Variations when N > K (very occasionally N < K for tangent viewing)
C    add Linearization of additional partial layer source term =
C         Exact_Scat * L_Multiplier(k)

C  If N = K (profile) or KS = 0 (column), there are extra linearizations.
C  If N > K, N < K, transmittance + some multiplication linearizations

            IF ( N.EQ.K .OR. KS.EQ.0 ) THEN
              DO V = 1, N_GEOMETRIES
               DO O1 = 1, NPOLAR
                DO Q = 1, K_PARAMETERS
                 LSS = ZMAT_UP(V,N,O1,1) *L_UP_MULTIPLIERS_UT(UT,KS,Q,V)
     &           + L_ZMAT_UP(Q,V,N,O1,1) *  UP_MULTIPLIERS_UT(UT,V)
                 L_FINAL_SOURCE = FLUXVEC(1) * LSS
     &          +  UP_LOSTRANS_UT(UT,V)    * L_SS_CUMSOURCE(Q,V,O1)
     &         + L_UP_LOSTRANS_UT(UT,Q,V) *   SS_CUMSOURCE_UP(V,O1,NC)
                 L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                 IF ( DO_PROFILE_LINEARIZATION ) THEN
                   PROFILEWF_SS(Q,K,UTA,V,O1,UPIDX) = L_SSCORRECTION
                 ELSE IF ( DO_COLUMN_LINEARIZATION ) THEN
                   COLUMNWF_SS(Q,UTA,V,O1,UPIDX) = L_SSCORRECTION
                 ENDIF
                ENDDO
               ENDDO
              ENDDO
             ELSE
              DO V = 1, N_GEOMETRIES
               DO O1 = 1, NPOLAR
                DO Q = 1, K_PARAMETERS
                 LSS = ZMAT_UP(V,N,O1,1)*L_UP_MULTIPLIERS_UT(UT,KS,Q,V)
                 L_FINAL_SOURCE = FLUXVEC(1) * LSS
     &               +  UP_LOSTRANS_UT(UT,V) * L_SS_CUMSOURCE(Q,V,O1)
                 L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                 PROFILEWF_SS(Q,K,UTA,V,O1,UPIDX) = L_SSCORRECTION
                ENDDO
               ENDDO
              ENDDO
             ENDIF

C  Ongrid output---------
C  just set to the cumulative source term 

           ELSE

            IF ( DO_PROFILE_LINEARIZATION ) THEN
             DO V = 1, N_GEOMETRIES
              DO O1 = 1, NPOLAR
               DO Q = 1, K_PARAMETERS
                L_FINAL_SOURCE = L_SS_CUMSOURCE(Q,V,O1)
                L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                PROFILEWF_SS(Q,K,UTA,V,O1,UPIDX) = L_SSCORRECTION
               ENDDO
              ENDDO
             ENDDO
            ELSE IF ( DO_COLUMN_LINEARIZATION ) THEN
             DO V = 1, N_GEOMETRIES
              DO O1 = 1, NPOLAR
               DO Q = 1, K_PARAMETERS
                L_FINAL_SOURCE = L_SS_CUMSOURCE(Q,V,O1)
                L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                COLUMNWF_SS(Q,UTA,V,O1,UPIDX) = L_SSCORRECTION
               ENDDO
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

C  debug

c     do k = 1, nlayers
c       write(85,*)k,ATMOSWF_SS(1,K,1,1,1,UPIDX),
c    &               ATMOSWF_SS(1,K,3,1,1,UPIDX)
c     enddo

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
           DO O1 = 1, NPOLAR
            DO Q = 1, K_PARAMETERS
             L_SS_CUMSOURCE(Q,V,O1) = ZERO
            ENDDO
           ENDDO
          ENDDO

C  initialise optical depth loop

          NSTART = 1
          NUT_PREV = NSTART - 1

C  Main loop over all output optical depths

          DO UTA = 1, N_USER_LEVELS

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

             IF ( N .EQ. K .OR. KS.EQ.0 ) THEN
              DO V = 1, N_GEOMETRIES
               DO O1 = 1, NPOLAR
                DO Q = 1, K_PARAMETERS
                 LSS = ZMAT_DN(V,N,O1,1)   * L_DN_MULTIPLIERS(N,KS,Q,V)
     &             + L_ZMAT_DN(Q,V,N,O1,1) *   DN_MULTIPLIERS(N,V)
                 L_SS_CUMSOURCE(Q,V,O1) = FLUXVEC(1) * LSS
     &            +  DN_LOSTRANS(N,V)    * L_SS_CUMSOURCE(Q,V,O1)
     &            + L_DN_LOSTRANS(N,Q,V) *   SS_CUMSOURCE_DN(V,O1,NC-1)
                ENDDO
               ENDDO
              ENDDO
             ELSE
              DO V = 1, N_GEOMETRIES
               DO O1 = 1, NPOLAR
                DO Q = 1, K_PARAMETERS
                 LSS = ZMAT_DN(V,N,O1,1) * L_DN_MULTIPLIERS(N,KS,Q,V)
                 L_SS_CUMSOURCE(Q,V,O1) = FLUXVEC(1) * LSS
     &            +  DN_LOSTRANS(N,V) * L_SS_CUMSOURCE(Q,V,O1)
                ENDDO
               ENDDO
              ENDDO
             ENDIF

C  End layer loop

            ENDDO

C  Offgrid output----------
C  Set final cumulative source and Single scatter Weighting function

           IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)

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
               DO O1 = 1, NPOLAR
                DO Q = 1, K_PARAMETERS
                 LSS = ZMAT_DN(V,N,O1,1) *L_DN_MULTIPLIERS_UT(UT,KS,Q,V)
     &           + L_ZMAT_DN(Q,V,N,O1,1) *  DN_MULTIPLIERS_UT(UT,V)
                 L_FINAL_SOURCE = FLUXVEC(1) * LSS
     &          +  DN_LOSTRANS_UT(UT,V)    * L_SS_CUMSOURCE(Q,V,O1)
     &         + L_DN_LOSTRANS_UT(UT,Q,V) *   SS_CUMSOURCE_DN(V,O1,NC)
                 L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                 IF ( DO_PROFILE_LINEARIZATION ) THEN
                   PROFILEWF_SS(Q,K,UTA,V,O1,DNIDX) = L_SSCORRECTION
                 ELSE IF ( DO_COLUMN_LINEARIZATION ) THEN
                   COLUMNWF_SS(Q,UTA,V,O1,DNIDX) = L_SSCORRECTION
                 ENDIF
                ENDDO
               ENDDO
              ENDDO
            ELSE
              DO V = 1, N_GEOMETRIES
               DO O1 = 1, NPOLAR
                DO Q = 1, K_PARAMETERS
                 LSS = ZMAT_DN(V,N,O1,1)*L_DN_MULTIPLIERS_UT(UT,KS,Q,V)
                 L_FINAL_SOURCE = FLUXVEC(1) * LSS
     &               +  DN_LOSTRANS_UT(UT,V) * L_SS_CUMSOURCE(Q,V,O1)
                 L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                 PROFILEWF_SS(Q,K,UTA,V,O1,DNIDX) = L_SSCORRECTION
                ENDDO
               ENDDO
              ENDDO
            ENDIF

C  Ongrid output---------
C  just set to the cumulative source term 

           ELSE

            IF ( DO_PROFILE_LINEARIZATION ) THEN
             DO V = 1, N_GEOMETRIES
              DO O1 = 1, NPOLAR
               DO Q = 1, K_PARAMETERS
                 L_FINAL_SOURCE = L_SS_CUMSOURCE(Q,V,O1)
                 L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                 PROFILEWF_SS(Q,K,UTA,V,O1,DNIDX) = L_SSCORRECTION
               ENDDO
              ENDDO
             ENDDO
            ELSE IF ( DO_COLUMN_LINEARIZATION ) THEN
             DO V = 1, N_GEOMETRIES
              DO O1 = 1, NPOLAR
               DO Q = 1, K_PARAMETERS
                 L_FINAL_SOURCE = L_SS_CUMSOURCE(Q,V,O1)
                 L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                 COLUMNWF_SS(Q,UTA,V,O1,DNIDX) = L_SSCORRECTION
               ENDDO
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

      SUBROUTINE SSCORR_OUTGOING_L_ZMATRIX
     I  ( DO_LINEAR, DO_SCATMAT, LAYER, NSTOKES, NGREEKMOMS, 
     I    LAYER_MAXMOMENTS, GREEKMATRIX, L_GREEKMATRIX,
     I    DO_SSCORR_TRUNCATION, SSFDEL, L_SSFDEL,
     I    LAYER_VARY_FLAG, LAYER_VARY_NUMBER,
     I    CTHETA, STHETA, CALPHA, SALPHA, CPHI,
     I    PHI, COSSCAT, VSIGN,  
     O    ZMAT, L_ZMAT, FMAT )

C  Include file
C  ------------

      INCLUDE '../includes/VLIDORT.PARS'

C  input
C  -----

C  Linearization control

      LOGICAL          DO_LINEAR
      LOGICAL          DO_SCATMAT ( MAXLAYERS, MAX_ATMOSWFS )
      LOGICAL          LAYER_VARY_FLAG
      INTEGER          LAYER_VARY_NUMBER

C  control integers

      INTEGER          LAYER, NSTOKES, NGREEKMOMS

C  scattering input information

      INTEGER          LAYER_MAXMOMENTS
      DOUBLE PRECISION GREEKMATRIX ( 0:MAXMOMENTS_INPUT, MAXLAYERS, 16 )
      DOUBLE PRECISION L_GREEKMATRIX
     &      ( MAX_ATMOSWFS, 0:MAXMOMENTS_INPUT, MAXLAYERS, 16 )

C  Truncation control

      LOGICAL          DO_SSCORR_TRUNCATION
      DOUBLE PRECISION SSFDEL   ( MAXLAYERS )
      DOUBLE PRECISION L_SSFDEL ( MAXLAYERS,MAX_ATMOSWFS )

C  zenith angle, azimuth angle, scatter angle (sines/cosines)

      DOUBLE PRECISION CTHETA
      DOUBLE PRECISION STHETA
      DOUBLE PRECISION CALPHA
      DOUBLE PRECISION SALPHA
      DOUBLE PRECISION CPHI
      DOUBLE PRECISION COSSCAT
      DOUBLE PRECISION VSIGN

C  Azimuth angle. Added for Version 2.4R, required for PHI > 180.
C    R. Spurr and V. Natraj, 01 May 2009

      DOUBLE PRECISION PHI

C  output
C  ------

C  Z-Matrix, first column only. F-matrix
C    (this module only works with sunlight)

      DOUBLE PRECISION    ZMAT ( 4 )
      DOUBLE PRECISION  L_ZMAT ( 4, MAX_ATMOSWFS )
      DOUBLE PRECISION  FMAT ( 6 )

C  Local
C  -----
       
      DOUBLE PRECISION  C1, S1, C2, S2, SINSCAT, HELP_SINSCAT
      DOUBLE PRECISION  CSIG1, CSIG2, SSIG1, SSIG2, CSIG1_2, CSIG2_2
      DOUBLE PRECISION  DL, QROOT6, UUU, P00(2), P02(2)
      DOUBLE PRECISION  FAC1, FAC2, SQL4, SQL41, TMP1, TMP2
      DOUBLE PRECISION  L_FMAT ( MAX_ATMOSWFS, 6 )
      DOUBLE PRECISION  HELP1, HELP2, FDNL1, DNL1, FACT
      DOUBLE PRECISION  GK11, L_GK11(MAX_ATMOSWFS)
      DOUBLE PRECISION  GK12, L_GK12(MAX_ATMOSWFS)

      INTEGER           GREEKMAT_INDEX(6), K, L, LNEW, LOLD, ITMP, N, Q
      INTEGER           INDEX_11, INDEX_12, INDEX_34
      INTEGER           INDEX_22, INDEX_33, INDEX_44

C  indexing key

C      GREEKMAT_INDEX(1) = 1  ---> INDEX_11
C      GREEKMAT_INDEX(2) = 6  ---> INDEX_22
C      GREEKMAT_INDEX(3) = 2  ---> INDEX_12
C      GREEKMAT_INDEX(4) = 11 ---> INDEX_33
C      GREEKMAT_INDEX(5) = 12 ---> INDEX_34
C      GREEKMAT_INDEX(6) = 16 ---> INDEX_44

      GREEKMAT_INDEX(1) = 1
      GREEKMAT_INDEX(2) = 6
      GREEKMAT_INDEX(3) = 2
      GREEKMAT_INDEX(4) = 11
      GREEKMAT_INDEX(5) = 12
      GREEKMAT_INDEX(6) = 16

      INDEX_11 = 1
      INDEX_12 = 2
      INDEX_22 = 3
      INDEX_33 = 4
      INDEX_34 = 5
      INDEX_44 = 6

      N = LAYER

C  Geometrical quantities
C  ----------------------

C  cosine scatter angle (this is valid only for non-refracting atmosphere)
C  VSIGN = -1 for upwelling, +1 for downwelling
C  Should be the same as the input value CSA

      UUU  = COSSCAT

C  Cosine Sigma 1 and 2. H/VdM, Eqs. (99)-(101)

C  a. safety
C    Watch for sin^2(scatter angle) less than zero (machine precision)
C    R. Spurr, 16 January 2006, RT SOLUTIONS Inc.

      HELP_SINSCAT = ( ONE - COSSCAT * COSSCAT )
      IF ( HELP_SINSCAT.LE.ZERO ) THEN
        SINSCAT = 1.0D-12
      ELSE
        SINSCAT = DSQRT ( HELP_SINSCAT )
      ENDIF

C  b. necessary limit analyses - Hovenier limits.
C     R. Spurr and V. Natraj, 17 January 2006

      IF ( DABS(SINSCAT) .LE. 1.0D-12 ) THEN
        CSIG1 = ZERO
        CSIG2 = ZERO
      ELSE
        IF ( STHETA .EQ. ZERO ) THEN
          CSIG1 = -  CPHI
        ELSE
          CSIG1 = (-VSIGN*CALPHA+CTHETA*COSSCAT)/SINSCAT/STHETA
        ENDIF
        IF ( SALPHA .EQ. ZERO ) THEN
          CSIG2 = -  CPHI
        ELSE
          CSIG2 = (-CTHETA+VSIGN*CALPHA*COSSCAT)/SINSCAT/SALPHA
        ENDIF
      ENDIF
      IF ( CSIG2 .GT. ONE  ) CSIG2 = ONE
      IF ( CSIG2 .LT. -ONE ) CSIG2 = -ONE

C  output, H/VdM, Eqs. (89)-(94)

      CSIG1_2 = TWO * CSIG1
      CSIG2_2 = TWO * CSIG2
      IF ( DABS(CSIG1-ONE).LT.1.0D-12)THEN
        SSIG1 = ZERO
      ELSE
        SSIG1 = DSQRT ( 1.0D0 - CSIG1 * CSIG1 )
      ENDIF
      IF ( DABS(CSIG2-ONE).LT.1.0D-12)THEN
        SSIG2 = ZERO
      ELSE
        SSIG2 = DSQRT ( 1.0D0 - CSIG2 * CSIG2 )
      ENDIF

C  For relazm in [180,360), need sign reversal for S1 and S2
C  See H/VdM, Eqs. 94-95

      C1 = CSIG1_2 * CSIG1 - ONE
      C2 = CSIG2_2 * CSIG2 - ONE
      IF ( PHI .LE. 180.0d0 ) THEN
        S1 = CSIG1_2 * SSIG1
        S2 = CSIG2_2 * SSIG2
      ELSE
        S1 = -CSIG1_2 * SSIG1
        S2 = -CSIG2_2 * SSIG2
      ENDIF

C  F-matrices
C  ----------

      QROOT6 = -0.25D0 * DSQRT(6.0D0)

C initialise F-matrix and its linearization

      DO K = 1, 6
        FMAT(K) = 0.0D0
      END DO
      IF ( DO_LINEAR ) THEN
       IF ( LAYER_VARY_FLAG ) THEN
        DO Q = 1, LAYER_VARY_NUMBER
         DO K = 1, 6
          L_FMAT(Q,K) = 0.0D0
         ENDDO
        ENDDO
       ENDIF
      ENDIF

C  Start loop over the coefficient index l
C  first update generalized spherical functions, then calculate coefs.
C  lold and lnew are pointer-like indices used in recurrence 

      LNEW = 1
      LOLD = 2

      DO L = 0, NGREEKMOMS

C  Local Greek matrix entries
C  Set the local Greek matrix elements that you need = 11 and 12.
c   44 and 34 are not required with natural sunlight (default here)
C   Linearized 44 and 34 Have not been CODED

        IF ( DO_SSCORR_TRUNCATION ) THEN
          DNL1  = DBLE(2*L + 1 )
          FDNL1 = SSFDEL(N) * DNL1
          FACT  = ONE - SSFDEL(N)
          GK11 = ( GREEKMATRIX(L,N,GREEKMAT_INDEX(1)) - FDNL1 ) / FACT
          GK12 =   GREEKMATRIX(L,N,GREEKMAT_INDEX(3)) / FACT
c          GK44 = ( GREEKMATRIX(L,N,GREEKMAT_INDEX(6)) - FDNL1 ) / FACT
c          GK34 =   GREEKMATRIX(L,N,GREEKMAT_INDEX(5)) / FACT
          IF ( DO_LINEAR ) THEN
           IF ( LAYER_VARY_FLAG ) THEN
            DO Q = 1, LAYER_VARY_NUMBER
             IF ( DO_SCATMAT(N,Q) ) THEN
               HELP1 = L_GREEKMATRIX(Q,L,N,1)*GREEKMATRIX(L,N,1)
               HELP2 = ( GK11 - DNL1 ) * L_SSFDEL(N,Q)
               L_GK11(Q) = ( HELP1 + HELP2 ) / FACT
               HELP1 = L_GREEKMATRIX(Q,L,N,2)*GREEKMATRIX(L,N,2)
               HELP2 = GK12 * L_SSFDEL(N,Q)
               L_GK12(Q) = ( HELP1 + HELP2 ) / FACT
             ELSE
               L_GK11(Q) = ZERO
               L_GK12(Q) = ZERO
             ENDIF
            ENDDO
           ENDIF
          ENDIF
        ELSE
          GK11 = GREEKMATRIX(L,N,GREEKMAT_INDEX(1))
          GK12 = GREEKMATRIX(L,N,GREEKMAT_INDEX(3))
c          GK34 = GREEKMATRIX(L,N,GREEKMAT_INDEX(5))
c          GK44 = GREEKMATRIX(L,N,GREEKMAT_INDEX(6))
          IF ( DO_LINEAR ) THEN
           IF ( LAYER_VARY_FLAG ) THEN
            DO Q = 1, LAYER_VARY_NUMBER
             IF ( DO_SCATMAT(N,Q) ) THEN
               L_GK11(Q) = L_GREEKMATRIX(Q,L,N,1)*GREEKMATRIX(L,N,1)
               L_GK12(Q) = L_GREEKMATRIX(Q,L,N,2)*GREEKMATRIX(L,N,2)
             ELSE
               L_GK11(Q) = ZERO
               L_GK12(Q) = ZERO
             ENDIF
            ENDDO
           ENDIF
          ENDIF
        ENDIF

C  First moment

        IF ( L .EQ. 0 ) THEN

C  Adding paper Eqs. (76) and (77) with m=0

          P00(LOLD) = ONE
          P00(LNEW) = ZERO
          P02(LOLD) = ZERO
          P02(LNEW) = ZERO

        ELSE

          DL   = DBLE(L)
          FAC1 = (TWO*DL-ONE)/DL
          FAC2 = (DL-ONE)/DL

C Adding paper Eq. (81) with m=0

          P00(LOLD) = FAC1*UUU*P00(LNEW) - FAC2*P00(LOLD)

        END IF

        IF ( L .EQ. 2 ) THEN

! Adding paper Eq. (78)  
! sql4 contains the factor dsqrt((l+1)*(l+1)-4) needed in
! the recurrence Eqs. (81) and (82)

          P02(LOLD) = QROOT6*(ONE-UUU*UUU)
          P02(LNEW) = ZERO
          SQL41 = ZERO

        ELSE IF ( L .GT. 2) THEN

! Adding paper Eq. (82) with m=0

          SQL4  = SQL41
          SQL41 = DSQRT(DL*DL-4.0D0)
          TMP1  = (2.0D0*DL-1.0D0)/SQL41
          TMP2  = SQL4/SQL41
          P02(LOLD) = TMP1*UUU*P02(LNEW) - TMP2*P02(LOLD)

        END IF

! Switch indices so that lnew indicates the function with
! the present index value l, this mechanism prevents swapping
! of entire arrays.

        ITMP = LNEW
        LNEW = LOLD
        LOLD = ITMP

! Now add the l-th term to the scattering matrix.
! See de Haan et al. (1987) Eqs. (68)-(73).
! Remember for Mie scattering : F11 = F22 and F33 = F44
C  remaining symmetries for Mie particles. Not needed
C  do not need these with natural sunlight:
c            FMAT(INDEX_34)and FMAT(INDEX_44)

        IF ( L.LE.LAYER_MAXMOMENTS ) THEN
          FMAT(INDEX_11) = FMAT(INDEX_11) + GK11 * P00(LNEW)
          FMAT(INDEX_12) = FMAT(INDEX_12) + GK12 * P02(LNEW)
c          FMAT(INDEX_44) = FMAT(INDEX_44) + GK44 * P00(LNEW)
c          FMAT(INDEX_34) = FMAT(INDEX_34) + GK34 * P02(LNEW)

C  Previous code............................
c          FMAT(INDEX_11) = FMAT(INDEX_11) +
c     &         GREEKMATRIX(L,N,1)*P00(LNEW)
c          FMAT(INDEX_12) = FMAT(INDEX_12) +
c     &         GREEKMATRIX(L,N,2)*P02(LNEW)
c          FMAT(INDEX_44) = FMAT(INDEX_44) +
c     &         GREEKMATRIX(L,N,16)*P00(LNEW)
c          FMAT(INDEX_34) = FMAT(INDEX_34) +
c     &         GREEKMATRIX(L,N,12)*P02(LNEW)

c          FMAT(INDEX_22) = FMAT(INDEX_11)
c          FMAT(INDEX_33) = FMAT(INDEX_44)
        ENDIF

C  Linearization

        IF ( DO_LINEAR ) THEN
         IF ( LAYER_VARY_FLAG ) THEN
          DO Q = 1, LAYER_VARY_NUMBER
           IF ( DO_SCATMAT(N,Q) ) THEN
            IF ( L.LE.LAYER_MAXMOMENTS ) THEN

             L_FMAT(Q,INDEX_11) = L_FMAT(Q,INDEX_11)+L_GK11(Q)*P00(LNEW)
             L_FMAT(Q,INDEX_12) = L_FMAT(Q,INDEX_12)+L_GK12(Q)*P02(LNEW)
c             L_FMAT(Q,INDEX_34) = L_FMAT(Q,INDEX_34)+L_GK34(Q)*P02(LNEW)
c             L_FMAT(Q,INDEX_44) = L_FMAT(Q,INDEX_44)+L_GK44(Q)*P00(LNEW)

C  Previous code
c             HELP = L_GREEKMATRIX(Q,L,N,1)*GREEKMATRIX(L,N,1)
c             L_FMAT(Q,INDEX_11) = L_FMAT(Q,INDEX_11) + HELP*P00(LNEW)
c             HELP = L_GREEKMATRIX(Q,L,N,2)*GREEKMATRIX(L,N,2)
c             L_FMAT(Q,INDEX_12) = L_FMAT(Q,INDEX_12) + HELP*P02(LNEW)
c             HELP = L_GREEKMATRIX(Q,L,N,16)*GREEKMATRIX(L,N,16)
c             L_FMAT(Q,INDEX_44) = L_FMAT(Q,INDEX_44) + HELP*P00(LNEW)
c             HELP = L_GREEKMATRIX(Q,L,N,12)*GREEKMATRIX(L,N,12)
c             L_FMAT(Q,INDEX_34) = L_FMAT(Q,INDEX_34) + HELP*P02(LNEW)

c             L_FMAT(Q,INDEX_22) = L_FMAT(Q,INDEX_11)
c             L_FMAT(Q,INDEX_33) = L_FMAT(Q,INDEX_44)

            ENDIF
           ENDIF
          ENDDO
         ENDIF
        ENDIF

C  End moments

      END DO

C  Z-matrix, first column only
C  ---------------------------

      IF ( NSTOKES .EQ. 1 ) THEN
        ZMAT(1) =   FMAT(INDEX_11)
      ELSE
        ZMAT(1) =   FMAT(INDEX_11)
        ZMAT(2) =  -FMAT(INDEX_12) * C2
        ZMAT(3) =   FMAT(INDEX_12) * S2
        ZMAT(4) =   0.0D0
      ENDIF

C  Linearized Z-matrix
C  -------------------

C  Copy for Stokes = 1
C  For polarized case, sunlight only !!!!!!!
C   For sunlight, only need the first Column of the Z-matrix

      IF ( DO_LINEAR ) THEN
       DO Q = 1, LAYER_VARY_NUMBER
         IF ( NSTOKES .EQ. 1 ) THEN
          L_ZMAT(1,Q) = L_FMAT(Q,INDEX_11)
         ELSE IF ( NSTOKES .GT. 1 ) THEN
          L_ZMAT(1,Q) =   L_FMAT(Q,INDEX_11)
          L_ZMAT(2,Q) =  -L_FMAT(Q,INDEX_12) * C2
          L_ZMAT(3,Q) =   L_FMAT(Q,INDEX_12) * S2
          L_ZMAT(4,Q) =   ZERO
         ENDIF
       ENDDO
      ENDIF

C  finish

      RETURN
      END
     
C

      subroutine l_outgoing_integration_up
     i   ( nlayers, nfinelayers,
     i     do_partials,
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

      include '../includes/VLIDORT.PARS'

C  control

      logical          do_partials
      integer          nfinelayers, nlayers
      integer          n_partials
      integer          partials_idx    (max_partlayers)
      integer          partials_fineidx(max_partlayers)

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

      integer          ntraverse_p(max_partlayers)
      double precision sunpaths_p (max_partlayers,maxlayers)
      double precision alpha_p    (max_partlayers)

C  Extinction inputs

      double precision extinction   (maxlayers)
      double precision L_extinction (maxlayers,max_atmoswfs)

C  outputs
C  -------

      double precision multipliers  (maxlayers)
      double precision lostrans     (maxlayers)
      double precision L_multipliers(maxlayers,0:maxlayers,max_atmoswfs)
      double precision L_lostrans   (maxlayers,max_atmoswfs)

      double precision multipliers_p   (max_partlayers)
      double precision lostrans_p      (max_partlayers)
      double precision L_multipliers_p 
     &         (max_partlayers,0:maxlayers,max_atmoswfs)
      double precision L_lostrans_p
     &         (max_partlayers,max_atmoswfs)

      double precision boa_attn
      double precision L_boa_attn    (0:maxlayers,max_atmoswfs)

C  local arrays
C  ------------

C  Local geoemetry arrays

      double precision csq_fine ( maxfinelayers)
      double precision cot_fine ( maxfinelayers)

C  Local attenuation factors

      double precision attn ( 0:maxlayers )
      double precision attn_fine ( maxlayers, maxfinelayers )
      double precision attn_p    ( max_partlayers )

C  Local attenuation factors, linearizations

      double precision l_attn ( 0:maxlayers, 0:maxlayers, max_atmoswfs )
      double precision l_attn_fine ( maxlayers, maxfinelayers,
     &                               0:maxlayers, max_atmoswfs )
      double precision l_attn_p
     &       ( max_partlayers, 0:maxlayers, max_atmoswfs  )

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

C  Linearized attenuation factors
C     Note the extra zeroing...............bug, 01 November 2007

      if ( do_profile_linearization ) then
        do n = 0, nlayers
          do k = n+1, nlayers
            do q = 1, layer_vary_number(k)
              l_attn(n,k,q) = zero
            enddo
          enddo
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
     i     do_partials,
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

      include '../includes/VLIDORT.PARS'

C  control

      logical          do_partials
      integer          nfinelayers, nlayers
      integer          n_partials
      integer          partials_idx    (max_partlayers)
      integer          partials_fineidx(max_partlayers)

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

      integer          ntraverse_p(max_partlayers)
      double precision sunpaths_p (max_partlayers,maxlayers)
      double precision alpha_p    (max_partlayers)

C  Extinction inputs

      double precision extinction   (maxlayers)
      double precision L_extinction (maxlayers,max_atmoswfs)

C  outputs
C  -------

      double precision multipliers  (maxlayers)
      double precision lostrans     (maxlayers)
      double precision L_multipliers(maxlayers,0:maxlayers,max_atmoswfs)
      double precision L_lostrans   (maxlayers,max_atmoswfs)

      double precision multipliers_p   (max_partlayers)
      double precision lostrans_p      (max_partlayers)
      double precision L_multipliers_p 
     &         (max_partlayers,0:maxlayers,max_atmoswfs)
      double precision L_lostrans_p
     &         (max_partlayers,max_atmoswfs)

C  local arrays
C  ------------

C  Local geoemetry arrays

      double precision csq_fine ( maxfinelayers)
      double precision cot_fine ( maxfinelayers)

C  Local attenuation factors

      double precision attn ( 0:maxlayers )
      double precision attn_fine ( maxlayers, maxfinelayers )
      double precision attn_p    ( max_partlayers )

C  Local attenuation factors, linearizations

      double precision l_attn ( 0:maxlayers, 0:maxlayers, max_atmoswfs )
      double precision l_attn_fine ( maxlayers, maxfinelayers,
     &                               0:maxlayers, max_atmoswfs )
      double precision l_attn_p
     &       ( max_partlayers, 0:maxlayers, max_atmoswfs  )

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
c         if (n.eq.24)write(*,'(2i5,1p5e17.7)')
c     i   n,v,DN_LOSTRANS(N,V),DN_MULTIPLIERS(N,V)

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
C    Coded up 26 September 2007. Requires checking.....................

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

      SUBROUTINE VLIDORT_LAP_LAMBERTIAN_DBCORR(FLUXMULT)

C  Prepares Linearization of Exact Direct Beam reflection (Lambertian case)
C   Linearization with respect to PROFILE atmospheric variables

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  Include file of setup variables (Input to the present module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  Include files of Linearized input and setup variables

      INCLUDE '../includes/VLIDORT_L_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'

C  Include files of Single scatter & surface reflectance variables (Input)

      INCLUDE '../includes/VLIDORT_SINGSCAT.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'

C  Direct beam output stored in linearized correction file

      INCLUDE '../includes/VLIDORT_L_SINGSCAT.VARS'

C  input argument

      DOUBLE PRECISION FLUXMULT

C  Local variables
C  ---------------

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER          UT, UTA, UM, UA, NC, IB, V, K, Q, O1
      DOUBLE PRECISION X0_FLUX, X0_BOA, FACTOR, L_ATTN, TR, LTR

C  first stage
C  -----------

C  return if no upwelling

      IF ( .NOT.DO_UPWELLING ) RETURN

C  Only 1 Stokes component

      O1 = 1

C  Lambertian factor

      FACTOR = LAMBERTIAN_ALBEDO * FLUXMULT

C  Start weighting function loop

      DO K = 1, NLAYERS
        IF ( LAYER_VARY_FLAG(K) ) THEN
          DO Q = 1, LAYER_VARY_NUMBER(K)

C  Initialize the output results

            DO V = 1, N_GEOMETRIES
              DO UTA = 1, N_USER_LEVELS
                PROFILEWF_DB(Q,K,UTA,V,O1)      = ZERO
              ENDDO
            ENDDO

C  New Code  R. Spurr, 6 August 2007. RT Solutions Inc.
C  ====================================================

C  initialize cumulative source term

C  Start geometry loops

            DO UM = 1, N_USER_STREAMS
             DO IB = 1, NBEAMS
              IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
               DO UA = 1, N_USER_RELAZMS
                V = VZA_OFFSETS(IB,UM) + UA

c  Beam attenuation and reflection

                IF ( DO_SSCORR_OUTGOING ) THEN
                 L_ATTN = ZERO
                 IF (BOA_ATTN(V).GT.ZERO ) THEN
                  X0_BOA = DCOS(SZANGLES_ADJUST(UM,IB,UA)*DEG_TO_RAD)
                  L_ATTN = L_BOA_ATTN(K,Q,V)
                 ENDIF
                ELSE
                 IF ( DO_REFRACTIVE_GEOMETRY ) THEN
                  X0_BOA = DCOS(SZA_LOCAL_INPUT(NLAYERS,IB)*DEG_TO_RAD)
                 ELSE
                  X0_BOA = COS_SZANGLES(IB)
                 ENDIF 
                 L_ATTN = - L_DELTAU_SLANT(Q,NLAYERS,K,IB) *
     &                      SOLAR_BEAM_OPDEP(IB) 
                ENDIF
                X0_FLUX = FOUR * X0_BOA
                L_ATTN  = L_ATTN * X0_FLUX             
                L_DB_CUMSOURCE(V,O1) = L_ATTN * FACTOR

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

            DO UTA = N_USER_LEVELS, 1, -1

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
                    V = VZA_OFFSETS(IB,UM) + UA
                    IF ( DO_SSCORR_OUTGOING ) THEN
                      TR  = UP_LOSTRANS(N,V)
                      LTR = L_UP_LOSTRANS(N,Q,V)
                    ELSE
                      TR  = T_DELT_USERM(N,UM)
                      LTR = L_T_DELT_USERM(N,UM,Q)
                    ENDIF
                    L_DB_CUMSOURCE(V,O1) = TR * L_DB_CUMSOURCE(V,O1)
     &                                  + LTR * DB_CUMSOURCE(V,O1,NC-1) 
                   ENDDO
                  ENDDO
                 ENDDO
                ELSE
                 DO IB = 1, NBEAMS
                  DO UM = 1, N_USER_STREAMS
                   DO UA = 1, N_USER_RELAZMS
                    V = VZA_OFFSETS(IB,UM) + UA
                    IF ( DO_SSCORR_OUTGOING ) THEN
                      TR  = UP_LOSTRANS(N,V)
                    ELSE
                      TR  = T_DELT_USERM(N,UM)
                    ENDIF
                    L_DB_CUMSOURCE(V,O1) = TR * L_DB_CUMSOURCE(V,O1)
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

              IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

               UT = PARTLAYERS_OUTINDEX(UTA)
               N  = PARTLAYERS_LAYERIDX(UT)

               IF ( N.EQ.K) THEN
                DO IB = 1, NBEAMS
                 DO UM = 1, N_USER_STREAMS
                  DO UA = 1, N_USER_RELAZMS
                   V = VZA_OFFSETS(IB,UM) + UA
                   IF ( DO_SSCORR_OUTGOING ) THEN
                    TR  = UP_LOSTRANS_UT(UT,V)
                    LTR = L_UP_LOSTRANS_UT(UT,Q,V)
                   ELSE
                    TR  = T_UTUP_USERM(UT,UM)
                    LTR = L_T_UTUP_USERM(UT,UM,Q)
                   ENDIF
                   PROFILEWF_DB(Q,K,UTA,V,O1) = TR*L_DB_CUMSOURCE(V,O1)
     &                                + LTR*  DB_CUMSOURCE(V,O1,NC)
                  ENDDO
                 ENDDO
                ENDDO
               ELSE IF ( N.NE.K) THEN
                DO IB = 1, NBEAMS
                 DO UM = 1, N_USER_STREAMS
                  DO UA = 1, N_USER_RELAZMS
                   V = VZA_OFFSETS(IB,UM) + UA
                   IF ( DO_SSCORR_OUTGOING ) THEN
                    TR  = UP_LOSTRANS_UT(UT,V)
                   ELSE
                    TR  = T_UTUP_USERM(UT,UM)
                   ENDIF
                   PROFILEWF_DB(Q,K,UTA,V,O1) = TR*L_DB_CUMSOURCE(V,O1)
                  ENDDO
                 ENDDO
                ENDDO
               ENDIF

C  Ongrid output : Set final cumulative source directly

              ELSE

               DO IB = 1, NBEAMS
                DO UM = 1, N_USER_STREAMS
                 DO UA = 1, N_USER_RELAZMS
                  V = VZA_OFFSETS(IB,UM) + UA
                   PROFILEWF_DB(Q,K,UTA,V,O1) = L_DB_CUMSOURCE(V,O1)
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

C  Finish

      RETURN
      END

C

      SUBROUTINE VLIDORT_LAC_LAMBERTIAN_DBCORR(FLUXMULT)

C  Prepares Linearization of Exact Direct Beam reflection (Lambertian case)
C   Linearization with respect to COLUMN atmospheric variables

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  Include file of setup variables (Input to the present module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  Include files of Linearized input and setup variables

      INCLUDE '../includes/VLIDORT_L_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'

C  Include files of Single scatter & surface reflectance variables (Input)

      INCLUDE '../includes/VLIDORT_SINGSCAT.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'

C  Direct beam output stored in linearized correction file

      INCLUDE '../includes/VLIDORT_L_SINGSCAT.VARS'

C  input argument

      DOUBLE PRECISION FLUXMULT

C  Local variables
C  ---------------

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER          UT, UTA, UM, UA, NC, IB, V, K, K0, Q, O1
      DOUBLE PRECISION X0_FLUX, X0_BOA, FACTOR, L_ATTN, TR, LTR

C  first stage
C  -----------

C  return if no upwelling

      IF ( .NOT.DO_UPWELLING ) RETURN

C  Only 1 Stokes component

      O1 = 1

C  Lambertian factor

      FACTOR = LAMBERTIAN_ALBEDO * FLUXMULT

C  Weighting function index

      K0 = 0

C  Start weighting function loop

      DO Q = 1, N_TOTALCOLUMN_WFS

C  Initialize the output results

        DO V = 1, N_GEOMETRIES
          DO UTA = 1, N_USER_LEVELS
            COLUMNWF_DB(Q,UTA,V,O1)      = ZERO
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
                V = VZA_OFFSETS(IB,UM) + UA

c  Beam attenuation and reflection

                IF ( DO_SSCORR_OUTGOING ) THEN
                 L_ATTN = ZERO
                 IF (BOA_ATTN(V).GT.ZERO ) THEN
                  X0_BOA = DCOS(SZANGLES_ADJUST(UM,IB,UA)*DEG_TO_RAD)
                  L_ATTN = L_BOA_ATTN(K0,Q,V)
                 ENDIF
                ELSE
                 IF ( DO_REFRACTIVE_GEOMETRY ) THEN
                  X0_BOA = DCOS(SZA_LOCAL_INPUT(NLAYERS,IB)*DEG_TO_RAD)
                 ELSE
                  X0_BOA = COS_SZANGLES(IB)
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
                L_DB_CUMSOURCE(V,O1) = L_ATTN * FACTOR

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

        DO UTA = N_USER_LEVELS, 1, -1

C  Layer index for given optical depth

          NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
          NUT    = NLEVEL + 1

C  Cumulative layer transmittance :
C    loop over layers working upwards to level NUT
C   addition linearization of transmittance

          DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                DO UA = 1, N_USER_RELAZMS
                  V = VZA_OFFSETS(IB,UM) + UA
                  IF ( DO_SSCORR_OUTGOING ) THEN
                    TR  = UP_LOSTRANS(N,V)
                    LTR = L_UP_LOSTRANS(N,Q,V)
                  ELSE
                    TR  = T_DELT_USERM(N,UM)
                    LTR = L_T_DELT_USERM(N,UM,Q)
                  ENDIF
                  L_DB_CUMSOURCE(V,O1) = TR * L_DB_CUMSOURCE(V,O1)
     &                                + LTR * DB_CUMSOURCE(V,O1,NC-1) 
                ENDDO
              ENDDO
            ENDDO

C  end layer loop

          ENDDO

C  Offgrid output
C  ------------- 

C  Require partial layer transmittance
C  require additional linearization

          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)

            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                DO UA = 1, N_USER_RELAZMS
                  V = VZA_OFFSETS(IB,UM) + UA
                  IF ( DO_SSCORR_OUTGOING ) THEN
                   TR  = UP_LOSTRANS_UT(UT,V)
                   LTR = L_UP_LOSTRANS_UT(UT,Q,V)
                  ELSE
                   TR  = T_UTUP_USERM(UT,UM)
                   LTR = L_T_UTUP_USERM(UT,UM,Q)
                  ENDIF
                  COLUMNWF_DB(Q,UTA,V,O1) = TR * L_DB_CUMSOURCE(V,O1)
     &                                   + LTR * DB_CUMSOURCE(V,O1,NC)
                ENDDO
              ENDDO
            ENDDO

C  Ongrid output : Set final cumulative source directly

          ELSE

            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                DO UA = 1, N_USER_RELAZMS
                  V = VZA_OFFSETS(IB,UM) + UA
                  COLUMNWF_DB(Q,UTA,V,O1) = L_DB_CUMSOURCE(V,O1)
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

      SUBROUTINE VLIDORT_LS_LAMBERTIAN_DBCORR

C  Prepares Linearization of Exact Direct Beam reflection (Lambertian case)
C   Linearization with respect to LAMBERTIAN amplitude variable only !

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  Include file of setup variables (Input to the present module)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  Include files of Linearized input and setup variables

      INCLUDE '../includes/VLIDORT_L_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'

C  Include files of Single scatter stuff (Input)

      INCLUDE '../includes/VLIDORT_SINGSCAT.VARS'

C  Direct beam output stored in linearized correction file

      INCLUDE '../includes/VLIDORT_L_SINGSCAT.VARS'

C  Local variables
C  ---------------

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER          UT, UTA, UM, UA, NC, IB, V, WF, O1
      DOUBLE PRECISION TR

C  first stage
C  -----------

C  return if no upwelling

      IF ( .NOT.DO_UPWELLING ) RETURN

C  Only 1 Stokes component, only one Surface weighting function

      O1 = 1
      WF = 1

C  initialise output

      DO V = 1, N_GEOMETRIES
        DO UTA = 1, N_USER_LEVELS
          SURFACEWF_DB(WF,UTA,V,O1)  = ZERO
        ENDDO
      ENDDO

C  Transmittance of source term: upwelling recursion
C  -------------------------------------------------

C  initialize cumulative linearized source term
C    In this case, just equal to the original source term
C      ( as dependence on Lambertin albedo is linear )

      NC =  0
      DO V = 1, N_GEOMETRIES
        LS_DB_CUMSOURCE(V,O1) = DB_CUMSOURCE(V,O1,NC)
      ENDDO

C  initialize optical depth loop

      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

C  Main loop over all output optical depths

      DO UTA = N_USER_LEVELS, 1, -1

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
                V = VZA_OFFSETS(IB,UM) + UA
                IF ( DO_SSCORR_OUTGOING ) THEN
                  TR  = UP_LOSTRANS(N,V)
                ELSE IF ( DO_SSCORR_NADIR ) THEN
                  TR  = T_DELT_USERM(N,UM)
                ENDIF
                LS_DB_CUMSOURCE(V,O1) = TR * LS_DB_CUMSOURCE(V,O1)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

C  Offgrid output
C  -------------- 

C  Require partial layer transmittance

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)
          DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
              DO UA = 1, N_USER_RELAZMS
                V = VZA_OFFSETS(IB,UM) + UA
                IF ( DO_SSCORR_OUTGOING ) THEN
                  TR  = UP_LOSTRANS_UT(UT,V)
                ELSE IF ( DO_SSCORR_NADIR ) THEN
                  TR  = T_UTUP_USERM(UT,UM)
                ENDIF
                SURFACEWF_DB(WF,UTA,V,O1) = TR * LS_DB_CUMSOURCE(V,O1)
              ENDDO
            ENDDO
          ENDDO

C  Ongrid output : Set final cumulative source directly

        ELSE
          DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
              DO UA = 1, N_USER_RELAZMS
                V = VZA_OFFSETS(IB,UM) + UA
                SURFACEWF_DB(WF,UTA,V,O1) = LS_DB_CUMSOURCE(V,O1)
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


