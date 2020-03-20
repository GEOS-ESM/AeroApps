
! ###############################################################
! #                                                             #
! #                    THE VLIDORT  MODEL                       #
! #                                                             #
! #  Vectorized LInearized Discrete Ordinate Radiative Transfer #
! #  -          --         -        -        -         -        #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Author :      Robert. J. D. Spurr                          #
! #                                                             #
! #  Address :     RT Solutions, inc.                           #
! #                9 Channing Street                            #
! #                Cambridge, MA 02138, USA                     #
! #                Tel: (617) 492 1183                          #
! #                                                             #
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  Versions     :   2.0, 2.2, 2.3, 2.4, 2.4R, 2.4RT, 2.4RTC,  #
! #                   2.5, 2.6, 2.7                             #
! #  Release Date :   December 2005  (2.0)                      #
! #  Release Date :   March 2007     (2.2)                      #
! #  Release Date :   October 2007   (2.3)                      #
! #  Release Date :   December 2008  (2.4)                      #
! #  Release Date :   April 2009     (2.4R)                     #
! #  Release Date :   July 2009      (2.4RT)                    #
! #  Release Date :   October 2010   (2.4RTC)                   #
! #  Release Date :   March 2011     (2.5)                      #
! #  Release Date :   May 2012       (2.6)                      #
! #  Release Date :   August 2014    (2.7)                      #
! #                                                             #
! #       NEW: TOTAL COLUMN JACOBIANS         (2.4)             #
! #       NEW: BPDF Land-surface KERNELS      (2.4R)            #
! #       NEW: Thermal Emission Treatment     (2.4RT)           #
! #       Consolidated BRDF treatment         (2.4RTC)          #
! #       f77/f90 Release                     (2.5)             #
! #       External SS / New I/O Structures    (2.6)             #
! #                                                             #
! #       SURFACE-LEAVING / BRDF-SCALING      (2.7)             #
! #       TAYLOR Series / OMP THREADSAFE      (2.7)             #
! #                                                             #
! ###############################################################

!    #####################################################
!    #                                                   #
!    #   This Version of VLIDORT comes with a GNU-style  #
!    #   license. Please read the license carefully.     #
!    #                                                   #
!    #####################################################

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! # --------------------------                                  #
! #                                                             #
! #   1. Nadir SS-correction                                    #
! #       Version 2.4.  LC module  (column Jacobians)           #
! #                                                             #
! #            VLIDORT_LC_SSCORR_NADIR (master, new 2.4)        #
! #                                                             #
! #   2. Outgoing SS-correction                                 #
! #      Version 2.2   Whole layer integration                  #
! #              2.3.  partial-layer integration                #
! #      Version 2.4.  Column Jacobians introduced              #
! #                                                             #
! #            VLIDORT_LC_SSCORR_OUTGOING (master)              #
! #              LC_OUTGOING_INTEGRATION_UP                     #
! #              LC_OUTGOING_INTEGRATION_DN                     #
! #                                                             #
! #   3. DB correction, Lambertian and BRDF surfaces            #
! #      Version 2.4.  LAC module (column Jacobians)            #
! #                                                             #
! #            VLIDORT_LAC_DBCORRECTION                         #
! #                                                             #
! #            VFO_LCS_MASTER_INTERFACE                         #
! #                                                             #
! ###############################################################


      MODULE vlidort_lc_corrections

      USE FO_VectorSS_LinMasters

      PRIVATE
      PUBLIC :: VLIDORT_LC_SSCORR_NADIR, &
                VLIDORT_LC_SSCORR_OUTGOING, &
                VLIDORT_LAC_DBCORRECTION,&
                VFO_LCS_MASTER_INTERFACE

      CONTAINS

      SUBROUTINE VLIDORT_LC_SSCORR_NADIR ( &
        SSFLUX, NV_PARAMETERS, &
        DO_SSCORR_TRUNCATION, DO_REFRACTIVE_GEOMETRY, &
        DO_DELTAM_SCALING, DO_UPWELLING, &
        DO_DNWELLING, DO_OBSERVATION_GEOMETRY, &
        NSTOKES, NLAYERS, NGREEK_MOMENTS_INPUT, &
        N_USER_RELAZMS, USER_RELAZMS, &
        N_USER_LEVELS, GREEKMAT_TOTAL_INPUT, &
        FLUXVEC, COS_SZANGLES, &
        SIN_SZANGLES, SUN_SZA_COSINES, &
        NMOMENTS, NBEAMS, &
        N_USER_STREAMS, LAYER_MAXMOMENTS, &
        USER_STREAMS, PARTLAYERS_OUTFLAG, &
        PARTLAYERS_OUTINDEX, UTAU_LEVEL_MASK_UP, &
        UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX, &
        STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, &
        N_GEOMETRIES, VZA_OFFSETS, &
        L_OMEGA_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT, &
        TRUNC_FACTOR, T_DELT_USERM, &
        T_UTDN_USERM, T_UTUP_USERM, &
        EMULT_UP, EMULT_DN, &
        UT_EMULT_UP, UT_EMULT_DN, &
        ZMAT_UP, ZMAT_DN, &
        TMS, SSFDEL, &
        SS_CUMSOURCE_UP, SS_CUMSOURCE_DN, &
        DO_SCATMAT_VARIATION, L_T_DELT_USERM, &
        L_T_UTDN_USERM, L_T_UTUP_USERM, &
        LC_EMULT_UP, LC_EMULT_DN, &
        LC_UT_EMULT_UP, LC_UT_EMULT_DN, &
        COLUMNWF_SS )

!  SINGLE SCATTER EXACT CALCULATION. NADIR VIEW
!   PROGRAMMED BY R. SPURR, RT SOLUTIONS INC.

!  VERSION 2.4. TOTAL COLUMN JACOBIANS, DECEMBER 2008
!               MODELED AFTER LIDORT CODE VERSION 3.3

      USE VLIDORT_PARS
      USE VLIDORT_LA_CORRECTIONS

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT (IN) :: SSFLUX
      INTEGER, INTENT (IN) ::          NV_PARAMETERS
      LOGICAL, INTENT (IN) ::          DO_SSCORR_TRUNCATION
      LOGICAL, INTENT (IN) ::          DO_REFRACTIVE_GEOMETRY
      LOGICAL, INTENT (IN) ::          DO_DELTAM_SCALING
      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_DNWELLING
      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NGREEK_MOMENTS_INPUT
      INTEGER, INTENT (IN) ::          N_USER_RELAZMS
      DOUBLE PRECISION, INTENT (IN) :: USER_RELAZMS  ( MAX_USER_RELAZMS )
      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      DOUBLE PRECISION, INTENT (IN) :: GREEKMAT_TOTAL_INPUT &
          ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION, INTENT (IN) :: FLUXVEC ( MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: COS_SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) :: SIN_SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) :: SUN_SZA_COSINES &
          ( MAXLAYERS, MAX_SZANGLES )
      INTEGER, INTENT (IN) ::          NMOMENTS
      INTEGER, INTENT (IN) ::          NBEAMS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          LAYER_MAXMOMENTS ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: USER_STREAMS  ( MAX_USER_STREAMS )
      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_DN ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          N_GEOMETRIES
      INTEGER, INTENT (IN) ::          VZA_OFFSETS &
          ( MAX_SZANGLES, MAX_USER_VZANGLES )
      DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_TOTAL_INPUT &
          ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_GREEKMAT_TOTAL_INPUT &
           ( MAX_ATMOSWFS, 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION, INTENT (IN) :: TRUNC_FACTOR ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: EMULT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: EMULT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: UT_EMULT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: UT_EMULT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: ZMAT_UP &
          ( MAX_GEOMETRIES, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: ZMAT_DN &
          ( MAX_GEOMETRIES, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: TMS ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SSFDEL ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SS_CUMSOURCE_UP &
          ( MAX_GEOMETRIES, MAXSTOKES, 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SS_CUMSOURCE_DN &
          ( MAX_GEOMETRIES, MAXSTOKES, 0:MAXLAYERS )
      LOGICAL, INTENT (IN) ::          DO_SCATMAT_VARIATION &
          ( MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_USERM &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTUP_USERM &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_EMULT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_EMULT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_UT_EMULT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_UT_EMULT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (INOUT) :: COLUMNWF_SS &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, &
            MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )

!  LOCAL VARIABLES
!  ---------------

!  INDICES

      INTEGER ::          N, NUT, NSTART, NUT_PREV, NLEVEL, V, LUM
      INTEGER ::          UT, UTA, UM, IB, IA, NC, Q, NM1, NV, O1, O2

!  OTHER HELP VARIABLES

      DOUBLE PRECISION :: HELP, TTT, UVAR, AVAR, FT1, FT2, DNM1
      DOUBLE PRECISION :: L_FINAL_SOURCE, L_SS_LAYERSOURCE
      DOUBLE PRECISION :: L_SSCORRECTION, L_TR_CUMSOURCE
      DOUBLE PRECISION :: VAR_TMS(MAX_ATMOSWFS,MAXLAYERS)
      DOUBLE PRECISION :: MUX, SUM, VIEWSIGN

      DOUBLE PRECISION :: L_ZMAT_UP &
          ( MAX_ATMOSWFS, MAX_GEOMETRIES, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: L_ZMAT_DN &
          ( MAX_ATMOSWFS, MAX_GEOMETRIES, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: L_SS_CUMSOURCE &
          ( MAX_ATMOSWFS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION :: L_SSFDEL ( MAXLAYERS, MAX_ATMOSWFS )

!  MAY NOT REQUIRE THESE
!  ---------------------

!  ZENITH ANGLE COSINES/SINES, AZIMUTH ANGLE COSINES

      DOUBLE PRECISION :: CTHETA (MAXLAYERS,MAXBEAMS)
      DOUBLE PRECISION :: STHETA (MAXLAYERS,MAXBEAMS)
      DOUBLE PRECISION :: CALPHA (MAX_USER_STREAMS)
      DOUBLE PRECISION :: SALPHA (MAX_USER_STREAMS)
      DOUBLE PRECISION :: CPHI   (MAX_USER_RELAZMS)

!  SUNLIGHT PARAMETERS
!    - ASSUMES THE FLUX VECTOR IS (F,0,0,0)
!   - DROP THIS VARIABLES WHEN FURTHER TESTING HAS BEEN DONE

      LOGICAL, PARAMETER :: DO_SUNLIGHT = .TRUE.

!  SET UP OPERATIONS
!  -----------------

!  Local user index

      LUM = 1

!  ALL LAYERS

      DO NV = 1, NLAYERS

!  CREATE LINEARIZED TMS FACTOR FOR LAYER NV
!   ( USE UNSCALED LINEARIZED INPUTS - THE NAKAJIMA-TANAKA WAY)

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

!  ADDITIONAL DELTA-M SCALING
!  --------------------------

!  NEW SECTION. R. SPURR, 07 SEPTEMBER 2007.

        IF ( DO_SSCORR_TRUNCATION ) THEN
          NM1  = NGREEK_MOMENTS_INPUT
          DNM1 = DBLE(2*NM1+1)
          DO Q = 1, NV_PARAMETERS
            IF ( DO_SCATMAT_VARIATION(NV,Q) ) THEN
              FT1 = ONE / ( ONE- SSFDEL(NV) )
              L_SSFDEL(NV,Q) = GREEKMAT_TOTAL_INPUT(NM1,NV,1) * &
                        L_GREEKMAT_TOTAL_INPUT(Q,NM1,NV,1) / DNM1
              VAR_TMS(Q,NV) = VAR_TMS(Q,NV) - FT1 * L_SSFDEL(NV,Q)
            ENDIF
          ENDDO
        ENDIF

!  END LAYER LOOP

      ENDDO

!  SAVE SOME GEOMETRICAL QUANTITIES (COSINES AND SINES)
!    MULTIPLE LAYER QUANTITIES FOR THE REFRACTIVE CASE

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

      IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
        DO UM = 1, N_USER_STREAMS
          CALPHA(UM) = USER_STREAMS(UM)
          SALPHA(UM) = DSQRT ( ONE - CALPHA(UM) * CALPHA(UM) )
        ENDDO
        DO IA = 1, N_USER_RELAZMS
          CPHI(IA) = DCOS ( USER_RELAZMS(IA) * DEG_TO_RAD )
        ENDDO
      ELSE
        DO V = 1, N_GEOMETRIES
          CALPHA(V) = USER_STREAMS(V)
          SALPHA(V) = DSQRT ( ONE - CALPHA(V) * CALPHA(V) )
          CPHI(V) = DCOS ( USER_RELAZMS(V) * DEG_TO_RAD )
        ENDDO
      ENDIF

!  ####################
!  #    UPWELLING     #
!  ####################

      IF ( DO_UPWELLING ) THEN

!  =================================================
!  TOTAL SCATTERING MATRIX LINEARIZATION (UPWELLING)
!  =================================================

!  ALL LAYERS

        DO NV = 1, NLAYERS
         IF ( STERM_LAYERMASK_UP(NV)) THEN

!  GET L_Z MATRICES FOR EACH LAYER, IF THE SCATTER LAW IS LINEARIZED

          VIEWSIGN = -1.0D0

          CALL VLIDORTSS_L_ZMATRICES &
        ( DO_SUNLIGHT, DO_OBSERVATION_GEOMETRY, &
          NSTOKES, NGREEK_MOMENTS_INPUT, NV, NV_PARAMETERS, &
          N_GEOMETRIES, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, &
          DO_SCATMAT_VARIATION, LAYER_MAXMOMENTS, &
          GREEKMAT_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT, &
          DO_SSCORR_TRUNCATION, SSFDEL, L_SSFDEL, VIEWSIGN, VZA_OFFSETS, &
          CTHETA, STHETA, CALPHA, SALPHA, CPHI, USER_RELAZMS, &
          L_ZMAT_UP )

!  LOOP OVER VARYING PARAMETERS Q FOR LAYER NV
!  LOOP OVER ALL GEOMETRIES V

          DO Q = 1, NV_PARAMETERS
            DO V = 1, N_GEOMETRIES

!  PHASE FUNCTION MOMENT VARIATIONS
!    ADD TMS CORRECTION FACTOR LINEARIZATION

              IF ( DO_SCATMAT_VARIATION(NV,Q) ) THEN
                DO O1 = 1, NSTOKES
                  DO O2 = 1, NSTOKES
                    L_ZMAT_UP(Q,V,NV,O1,O2) = &
                      L_ZMAT_UP(Q,V,NV,O1,O2)  * TMS(NV) + &
                        ZMAT_UP(V,NV,O1,O2)    * VAR_TMS(Q,NV)
                  ENDDO
                ENDDO
              ELSE
                DO O1 = 1, NSTOKES
                  DO O2 = 1, NSTOKES
                    L_ZMAT_UP(Q,V,NV,O1,O2) = &
                        ZMAT_UP(V,NV,O1,O2) * VAR_TMS(Q,NV)
                  ENDDO
                ENDDO
              ENDIF

!  END PARAMETER LOOP AND VIEWING DIRECTIONS LOOP

            ENDDO
          ENDDO

!  ONLY IF LAYER NV EXISTS

         ENDIF
        ENDDO

!  ===================================
!  UPWELLING SINGLE SCATTER RECURRENCE
!  ===================================

!  INITIALIZE CUMULATIVE SOURCE TERM

        DO V = 1, N_GEOMETRIES
          DO Q = 1, NV_PARAMETERS
            DO O1 = 1, NSTOKES
              L_SS_CUMSOURCE(Q,V,O1) = ZERO
            ENDDO
          ENDDO
        ENDDO

!  INITIALISE OPTICAL DEPTH LOOP

        NC = 0
        NSTART = NLAYERS
        NUT_PREV = NSTART + 1

!  MAIN LOOP OVER ALL OUTPUT OPTICAL DEPTHS
!  ========================================

        DO UTA = N_USER_LEVELS, 1, -1

!  LAYER INDEX FOR GIVEN OPTICAL DEPTH

          NLEVEL = UTAU_LEVEL_MASK_UP(UTA)

!  FINISHING LAYER

          NUT = NLEVEL + 1

!  CUMULATIVE SINGLE SCATTER SOURCE TERMS TO LAYER NUT
!  ---------------------------------------------------

!    1. GET LAYER SOURCE TERMS = EXACT SCATTERING * MULTIPLIER
!    2. LOOP OVER LAYERS WORKING UPWARDS TO NUT

          DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N

!  ..ALL LAYERS WILL HAVE SOME VARIATION

            IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
              DO IB = 1, NBEAMS
                DO UM = 1, N_USER_STREAMS
                  DO IA = 1, N_USER_RELAZMS
                   V = VZA_OFFSETS(IB,UM) + IA
                   DO Q = 1, NV_PARAMETERS
                     DO O1 = 1, NSTOKES
                       SUM = ZERO
                       DO O2 = 1, NSTOKES
                        HELP = ZMAT_UP(V,N,O1,O2) * LC_EMULT_UP(UM,N,IB,Q) &
                           + L_ZMAT_UP(Q,V,N,O1,O2) * EMULT_UP(UM,N,IB)
                        SUM = SUM + HELP * FLUXVEC(O2)
                       ENDDO
                       L_SS_LAYERSOURCE = SUM
                       L_SS_CUMSOURCE(Q,V,O1) = L_SS_LAYERSOURCE &
                       +   T_DELT_USERM(N,UM)   * L_SS_CUMSOURCE(Q,V,O1) &
                       + L_T_DELT_USERM(N,UM,Q) *   SS_CUMSOURCE_UP(V,O1,NC-1)
                     ENDDO
                   ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ELSE
              DO V = 1, N_GEOMETRIES
                DO Q = 1, NV_PARAMETERS
                  DO O1 = 1, NSTOKES
                    SUM = ZERO
                    DO O2 = 1, NSTOKES
                      HELP = ZMAT_UP(V,N,O1,O2) * LC_EMULT_UP(LUM,N,V,Q) &
                         + L_ZMAT_UP(Q,V,N,O1,O2) * EMULT_UP(LUM,N,V)
                      SUM = SUM + HELP * FLUXVEC(O2)
                    ENDDO
                    L_SS_LAYERSOURCE = SUM
                    L_SS_CUMSOURCE(Q,V,O1) = L_SS_LAYERSOURCE &
                    +   T_DELT_USERM(N,V)   * L_SS_CUMSOURCE(Q,V,O1) &
                    + L_T_DELT_USERM(N,V,Q) *   SS_CUMSOURCE_UP(V,O1,NC-1)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF

!  END LAYER LOOP

          ENDDO

!  OFFGRID OUTPUT
!  --------------

!  SET FINAL CUMULATIVE SOURCE AND SINGLE SCATTER WEIGHTING FUNCTION

          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)

!    ADD LINEARIZATION OF ADDITIONAL PARTIAL LAYER SOURCE TERM =
!        L_EXACT_SCAT(N) * MULTIPLIER  +  EXACT_SCAT * L_MULTIPLIER(N)

            IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
              DO IB = 1, NBEAMS
                DO UM = 1, N_USER_STREAMS
                  DO IA = 1, N_USER_RELAZMS
                   V = VZA_OFFSETS(IB,UM) + IA
                   DO Q = 1, NV_PARAMETERS
                    DO O1 = 1, NSTOKES
                     SUM = ZERO
                     DO O2 = 1, NSTOKES
                      HELP = &
                          ZMAT_UP(V,N,O1,O2) * LC_UT_EMULT_UP(UM,UT,IB,Q) &
                      + L_ZMAT_UP(Q,V,N,O1,O2) * UT_EMULT_UP(UM,UT,IB)
                      SUM = SUM + HELP * FLUXVEC(O2)
                     ENDDO
                     L_SS_LAYERSOURCE = SUM
                     L_TR_CUMSOURCE = &
                     + L_T_UTUP_USERM(UT,UM,Q) *   SS_CUMSOURCE_UP(V,O1,NC) &
                     +   T_UTUP_USERM(UT,UM)   * L_SS_CUMSOURCE(Q,V,O1)
                     L_FINAL_SOURCE = L_TR_CUMSOURCE + L_SS_LAYERSOURCE
                     L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                     COLUMNWF_SS(Q,UTA,V,O1,UPIDX) = L_SSCORRECTION
                    ENDDO
                   ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ELSE
              DO V = 1, N_GEOMETRIES
                DO Q = 1, NV_PARAMETERS
                  DO O1 = 1, NSTOKES
                   SUM = ZERO
                   DO O2 = 1, NSTOKES
                    HELP = &
                        ZMAT_UP(V,N,O1,O2) * LC_UT_EMULT_UP(LUM,UT,V,Q) &
                    + L_ZMAT_UP(Q,V,N,O1,O2) * UT_EMULT_UP(LUM,UT,V)
                    SUM = SUM + HELP * FLUXVEC(O2)
                   ENDDO
                   L_SS_LAYERSOURCE = SUM
                   L_TR_CUMSOURCE = &
                   + L_T_UTUP_USERM(UT,V,Q) *   SS_CUMSOURCE_UP(V,O1,NC) &
                   +   T_UTUP_USERM(UT,V)   * L_SS_CUMSOURCE(Q,V,O1)
                   L_FINAL_SOURCE = L_TR_CUMSOURCE + L_SS_LAYERSOURCE
                   L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                   COLUMNWF_SS(Q,UTA,V,O1,UPIDX) = L_SSCORRECTION
                  ENDDO
                ENDDO
              ENDDO
            ENDIF

!  ONGRID OUTPUT
!  -------------

!  JUST SET TO THE CUMULATIVE SOURCE TERM

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

!  CHECK FOR UPDATING THE RECURSION

          IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT

!  END LOOP OVER OPTICAL DEPTH

        ENDDO

!  END UPWELLING CLAUSE

      ENDIF

!  ####################
!  #   DOWNWELLING    #
!  ####################

      IF ( DO_DNWELLING ) THEN

!  =================================================
!  TOTAL SCATTERING MATRIX LINEARIZATION (UPWELLING)
!  =================================================

!  ALL LAYERS

        DO NV = 1, NLAYERS
         IF ( STERM_LAYERMASK_DN(NV)) THEN

!  GET L_Z MATRICES FOR EACH LAYER, IF THE SCATTER LAW IS LINEARIZED

          VIEWSIGN = +1.0D0

          CALL VLIDORTSS_L_ZMATRICES &
        ( DO_SUNLIGHT, DO_OBSERVATION_GEOMETRY, &
          NSTOKES, NGREEK_MOMENTS_INPUT, NV, NV_PARAMETERS, &
          N_GEOMETRIES, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, &
          DO_SCATMAT_VARIATION, LAYER_MAXMOMENTS, &
          GREEKMAT_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT, &
          DO_SSCORR_TRUNCATION, SSFDEL, L_SSFDEL, VIEWSIGN, VZA_OFFSETS, &
          CTHETA, STHETA, CALPHA, SALPHA, CPHI, USER_RELAZMS, &
          L_ZMAT_DN )

!  LOOP OVER VARYING PARAMETERS Q FOR LAYER NV
!  LOOP OVER ALL GEOMETRIES V

          DO Q = 1, NV_PARAMETERS
           DO V = 1, N_GEOMETRIES

!  PHASE FUNCTION MOMENT VARIATIONS
!    ADD TMS CORRECTION FACTOR LINEARIZATION

            IF ( DO_SCATMAT_VARIATION(NV,Q) ) THEN
              DO O1 = 1, NSTOKES
                DO O2 = 1, NSTOKES
                  L_ZMAT_DN(Q,V,NV,O1,O2) = &
                      L_ZMAT_DN(Q,V,NV,O1,O2)  * TMS(NV) + &
                        ZMAT_DN(V,NV,O1,O2) * VAR_TMS(Q,NV)
                ENDDO
              ENDDO
            ELSE
              DO O1 = 1, NSTOKES
                DO O2 = 1, NSTOKES
                  L_ZMAT_DN(Q,V,NV,O1,O2) = &
                        ZMAT_DN(V,NV,O1,O2) * VAR_TMS(Q,NV)
                ENDDO
              ENDDO
            ENDIF

!  END PARAMETER AND GEOMETRY LOOPS

           ENDDO
          ENDDO

!  ONLY IF LAYER NV EXISTS

         ENDIF
        ENDDO

!  =====================================
!  DOWNWELLING SINGLE SCATTER RECURRENCE
!  =====================================

!  INITIALIZE CUMULATIVE SOURCE TERM

        DO V = 1, N_GEOMETRIES
          DO Q = 1, NV_PARAMETERS
            DO O1 = 1, NSTOKES
              L_SS_CUMSOURCE(Q,V,O1) = ZERO
            ENDDO
          ENDDO
        ENDDO

!  INITIALISE OPTICAL DEPTH LOOP

        NC = 0
        NSTART = 1
        NUT_PREV = NSTART - 1

!  MAIN LOOP OVER ALL OUTPUT OPTICAL DEPTHS
!  ========================================

        DO UTA = 1, N_USER_LEVELS

!  LAYER INDEX FOR GIVEN OPTICAL DEPTH

          NLEVEL = UTAU_LEVEL_MASK_DN(UTA)

!  FINISHING LAYER

          NUT = NLEVEL

!  CUMULATIVE SINGLE SCATTER SOURCE TERMS TO LAYER NUT
!  ---------------------------------------------------

!    1. GET LAYER SOURCE TERMS = EXACT SCATTERING * MULTIPLIER
!    2. LOOP OVER LAYERS WORKING UPWARDS TO NUT

          DO N = NSTART, NUT
            NC = N

            IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
              DO IB = 1, NBEAMS
                DO UM = 1, N_USER_STREAMS
                  DO IA = 1, N_USER_RELAZMS
                   V = VZA_OFFSETS(IB,UM) + IA
                   DO Q = 1, NV_PARAMETERS
                    DO O1 = 1, NSTOKES
                     SUM = ZERO
                     DO O2 = 1, NSTOKES
                      HELP = ZMAT_DN(V,N,O1,O2) * LC_EMULT_DN(UM,N,IB,Q) &
                         + L_ZMAT_DN(Q,V,N,O1,O2) * EMULT_DN(UM,N,IB)
                      SUM = SUM + HELP * FLUXVEC(O2)
                     ENDDO
                     L_SS_LAYERSOURCE = SUM
                     L_SS_CUMSOURCE(Q,V,O1) = L_SS_LAYERSOURCE &
                     +   T_DELT_USERM(N,UM)   * L_SS_CUMSOURCE(Q,V,O1) &
                     + L_T_DELT_USERM(N,UM,Q) *   SS_CUMSOURCE_DN(V,O1,NC-1)
                    ENDDO
                   ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ELSE
              DO V = 1, N_GEOMETRIES
                DO Q = 1, NV_PARAMETERS
                  DO O1 = 1, NSTOKES
                   SUM = ZERO
                   DO O2 = 1, NSTOKES
                    HELP = ZMAT_DN(V,N,O1,O2) * LC_EMULT_DN(LUM,N,V,Q) &
                       + L_ZMAT_DN(Q,V,N,O1,O2) * EMULT_DN(LUM,N,V)
                    SUM = SUM + HELP * FLUXVEC(O2)
                   ENDDO
                   L_SS_LAYERSOURCE = SUM
                   L_SS_CUMSOURCE(Q,V,O1) = L_SS_LAYERSOURCE &
                   +   T_DELT_USERM(N,V)   * L_SS_CUMSOURCE(Q,V,O1) &
                   + L_T_DELT_USERM(N,V,Q) *   SS_CUMSOURCE_DN(V,O1,NC-1)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF

!  END RECURSIVE LAYER LOOP

          ENDDO

!  OFFGRID OUTPUT
!  --------------

!  ADD ADDITIONAL PARTIAL LAYER SOURCE TERM = EXACT SCAT * MULTIPLIER
!  SET FINAL CUMULATIVE SOURCE AND CORRECT THE INTENSITY

          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)

!    ADD LINEARIZATION OF ADDITIONAL PARTIAL LAYER SOURCE TERM =
!        L_EXACT_SCAT(N) * MULTIPLIER  +  EXACT_SCAT * L_MULTIPLIER(N)

            IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
              DO IB = 1, NBEAMS
                DO UM = 1, N_USER_STREAMS
                  DO IA = 1, N_USER_RELAZMS
                   V = VZA_OFFSETS(IB,UM) + IA
                   DO Q = 1, NV_PARAMETERS
                    DO O1 = 1, NSTOKES
                     SUM = ZERO
                     DO O2 = 1, NSTOKES
                      HELP = &
                          ZMAT_DN(V,N,O1,O2) * LC_UT_EMULT_DN(UM,UT,IB,Q) &
                      + L_ZMAT_DN(Q,V,N,O1,O2) * UT_EMULT_DN(UM,UT,IB)
                      SUM = SUM + HELP * FLUXVEC(O2)
                     ENDDO
                     L_SS_LAYERSOURCE = SUM
                     L_TR_CUMSOURCE = &
                     + L_T_UTDN_USERM(UT,UM,Q) *   SS_CUMSOURCE_DN(V,O1,NC) &
                     +   T_UTDN_USERM(UT,UM)   * L_SS_CUMSOURCE(Q,V,O1)
                     L_FINAL_SOURCE = L_TR_CUMSOURCE + L_SS_LAYERSOURCE
                     L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                     COLUMNWF_SS(Q,UTA,V,O1,DNIDX) = L_SSCORRECTION
                    ENDDO
                   ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ELSE
              DO V = 1, N_GEOMETRIES
                DO Q = 1, NV_PARAMETERS
                  DO O1 = 1, NSTOKES
                   SUM = ZERO
                   DO O2 = 1, NSTOKES
                    HELP = &
                        ZMAT_DN(V,N,O1,O2) * LC_UT_EMULT_DN(LUM,UT,V,Q) &
                    + L_ZMAT_DN(Q,V,N,O1,O2) * UT_EMULT_DN(LUM,UT,V)
                    SUM = SUM + HELP * FLUXVEC(O2)
                   ENDDO
                   L_SS_LAYERSOURCE = SUM
                   L_TR_CUMSOURCE = &
                   + L_T_UTDN_USERM(UT,V,Q) *   SS_CUMSOURCE_DN(V,O1,NC) &
                   +   T_UTDN_USERM(UT,V)   * L_SS_CUMSOURCE(Q,V,O1)
                   L_FINAL_SOURCE = L_TR_CUMSOURCE + L_SS_LAYERSOURCE
                   L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                   COLUMNWF_SS(Q,UTA,V,O1,DNIDX) = L_SSCORRECTION
                  ENDDO
                ENDDO
              ENDDO
            ENDIF

!  ONGRID OUTPUT
!  -------------

!  JUST SET TO THE CUMULATIVE SOURCE TERM

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

!  CHECK FOR UPDATING THE RECURSION

          IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
          NUT_PREV = NUT

!  END LOOP OVER OPTICAL DEPTH

        ENDDO

!  END DOWNWELLING CLAUSE

      ENDIF

!  FINISH

      RETURN
      END SUBROUTINE VLIDORT_LC_SSCORR_NADIR

!

      SUBROUTINE VLIDORT_LC_SSCORR_OUTGOING ( &
        SSFLUX, DO_SSCORR_TRUNCATION, DO_DELTAM_SCALING, &
        DO_UPWELLING, DO_DNWELLING, DO_OBSERVATION_GEOMETRY, &
        NSTOKES, NLAYERS, &
        NGREEK_MOMENTS_INPUT, N_USER_RELAZMS, N_USER_LEVELS, &
        EARTH_RADIUS, HEIGHT_GRID, &
        OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT, &
        DO_TOA_CONTRIBS, FLUXVEC, NMOMENTS, &
        NBEAMS, N_USER_STREAMS, &
        LAYER_MAXMOMENTS, USER_VZANGLES_ADJUST, &
        SZANGLES_ADJUST, USER_RELAZMS_ADJUST, &
        N_PARTLAYERS, NFINELAYERS, &
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
        UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, &
        PARTLAYERS_LAYERIDX, &
        STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, &
        N_GEOMETRIES, VZA_OFFSETS, &
        DELTAU_VERT, PARTAU_VERT, TRUNC_FACTOR, &
        DO_COLUMN_LINEARIZATION, N_TOTALCOLUMN_WFS, &
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
        L_OMEGA_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT, &
        L_DELTAU_VERT, DO_SCATMAT_VARIATION, &
        UP_MULTIPLIERS, DN_MULTIPLIERS, UP_LOSTRANS, DN_LOSTRANS, &
        UP_MULTIPLIERS_UT, DN_MULTIPLIERS_UT, &
        UP_LOSTRANS_UT, DN_LOSTRANS_UT, &
        ZMAT_UP, ZMAT_DN, TMS, SSFDEL, &
        SS_CUMSOURCE_UP, SS_CUMSOURCE_DN, &
        BOA_ATTN, STOKES_SS, SS_CONTRIBS, &
        LC_UP_MULTIPLIERS, LC_DN_MULTIPLIERS, &
        L_UP_LOSTRANS, L_DN_LOSTRANS, &
        LC_UP_MULTIPLIERS_UT, LC_DN_MULTIPLIERS_UT, &
        L_UP_LOSTRANS_UT, L_DN_LOSTRANS_UT, &
        LC_BOA_ATTN, COLUMNWF_SS, &
        FAIL, MESSAGE, TRACE )

!  SINGLE SCATTER EXACT CALCULATION FOR THE OUTGOING LOS
!   THIS ONE WITH OPTIONAL LINEARIZATIONS
!         - NEW FOR VERSION 2.2

!   PROGRAMMED BY R. SPURR, RT SOLUTIONS INC.
!    FIRST DRAFT, JANUARY 31ST 2007
!   VALIDATED AGAINST TOMRAD, 29 MARCH 2007.
!   PARTIAL LAYER OUTPUT ADDED SEPTEMBER 2007

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  EXTENSION AND CHANGES 01-05 OCTOBER 2010, R. SPURR
!   - CORRECTIONS FOR THE PHI > 180 CASE
!   - ADDITION OF NON-MIE (SPHEROIDAL PARTICLE) ENTRIES TO F-MATRIX
!   - GENERALIZATION OF Z-MATRIX TO NON-SUNLIGHT SITUATIONS.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      USE VLIDORT_PARS
      USE VLIDORT_GEOMETRY
      USE VLIDORT_CORRECTIONS
      USE VLIDORT_LA_CORRECTIONS

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT (IN) :: SSFLUX
      LOGICAL, INTENT (IN) ::          DO_SSCORR_TRUNCATION
      LOGICAL, INTENT (IN) ::          DO_DELTAM_SCALING
      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_DNWELLING
      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NGREEK_MOMENTS_INPUT
      INTEGER, INTENT (IN) ::          N_USER_RELAZMS
      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      DOUBLE PRECISION, INTENT (IN) :: EARTH_RADIUS
      DOUBLE PRECISION, INTENT (IN) :: HEIGHT_GRID ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: OMEGA_TOTAL_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: GREEKMAT_TOTAL_INPUT &
            ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
      LOGICAL, INTENT (IN) ::          DO_TOA_CONTRIBS
      DOUBLE PRECISION, INTENT (IN) :: FLUXVEC ( MAXSTOKES )
      INTEGER, INTENT (IN) ::          NMOMENTS
      INTEGER, INTENT (IN) ::          NBEAMS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          LAYER_MAXMOMENTS ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: USER_VZANGLES_ADJUST &
           ( MAX_USER_VZANGLES )
      DOUBLE PRECISION, INTENT (IN) :: SZANGLES_ADJUST &
           ( MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS )
      DOUBLE PRECISION, INTENT (IN) :: USER_RELAZMS_ADJUST &
           ( MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS )
      INTEGER, INTENT(IN) ::           N_PARTLAYERS
      INTEGER, INTENT(IN) ::           NFINELAYERS
      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_DN ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          N_GEOMETRIES
      INTEGER, INTENT (IN) ::          VZA_OFFSETS &
          ( MAX_SZANGLES, MAX_USER_VZANGLES )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PARTAU_VERT ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: TRUNC_FACTOR ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::          DO_COLUMN_LINEARIZATION
      INTEGER, INTENT (IN) ::          N_TOTALCOLUMN_WFS
      LOGICAL, INTENT (IN) ::          LAYER_VARY_FLAG  ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          LAYER_VARY_NUMBER ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_TOTAL_INPUT &
          ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_GREEKMAT_TOTAL_INPUT &
          ( MAX_ATMOSWFS, 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      LOGICAL, INTENT (IN) ::          DO_SCATMAT_VARIATION &
          ( MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (OUT) ::  UP_MULTIPLIERS &
          ( MAXLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (OUT) ::  DN_MULTIPLIERS &
          ( MAXLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (OUT) ::  UP_LOSTRANS &
          ( MAXLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (OUT) ::  DN_LOSTRANS &
          ( MAXLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (OUT) ::  UP_MULTIPLIERS_UT &
          ( MAX_PARTLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (OUT) ::  DN_MULTIPLIERS_UT &
          ( MAX_PARTLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (OUT) ::  UP_LOSTRANS_UT &
          ( MAX_PARTLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (OUT) ::  DN_LOSTRANS_UT &
          ( MAX_PARTLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (OUT) ::  ZMAT_UP &
          ( MAX_GEOMETRIES, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) ::  ZMAT_DN &
          ( MAX_GEOMETRIES, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) ::  TMS ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) ::  SSFDEL ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) ::  SS_CUMSOURCE_UP &
          ( MAX_GEOMETRIES, MAXSTOKES, 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) ::  SS_CUMSOURCE_DN &
          ( MAX_GEOMETRIES, MAXSTOKES, 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) ::  BOA_ATTN ( MAX_GEOMETRIES )

      DOUBLE PRECISION, INTENT (INOUT) ::  STOKES_SS &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )

      DOUBLE PRECISION, INTENT (OUT) ::  SS_CONTRIBS &
          ( MAX_GEOMETRIES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) ::  LC_UP_MULTIPLIERS &
         ( MAXLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (OUT) ::  LC_DN_MULTIPLIERS &
         ( MAXLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (OUT) ::  L_UP_LOSTRANS &
         ( MAXLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (OUT) ::  L_DN_LOSTRANS &
         ( MAXLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (OUT) ::  LC_UP_MULTIPLIERS_UT &
          ( MAX_PARTLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (OUT) ::  LC_DN_MULTIPLIERS_UT &
          ( MAX_PARTLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (OUT) ::  L_UP_LOSTRANS_UT &
          ( MAX_PARTLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (OUT) ::  L_DN_LOSTRANS_UT &
          ( MAX_PARTLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )

      DOUBLE PRECISION, INTENT (OUT) ::  LC_BOA_ATTN &
          ( MAX_ATMOSWFS, MAX_GEOMETRIES )

      DOUBLE PRECISION, INTENT (INOUT) ::  COLUMNWF_SS &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, &
            MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )

      LOGICAL, INTENT (OUT) ::           FAIL
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  GEOMETRY ROUTINE INPUTS AND OUTPUTS
!  -----------------------------------

      LOGICAL ::          DO_FINE
      LOGICAL ::          DO_PARTIALS
      DOUBLE PRECISION :: ALPHA_BOA, THETA_BOA, PHI_BOA

!  MAIN OUTPUTS (GEOMETRY)

      INTEGER          :: NTRAVERSE(0:MAXLAYERS)
      DOUBLE PRECISION :: SUNPATHS(0:MAXLAYERS,MAXLAYERS)
      DOUBLE PRECISION :: RADII   (0:MAXLAYERS)
      DOUBLE PRECISION :: ALPHA_ALL  (0:MAXLAYERS)

!  FINE LEVEL OUTPUT (GEOMETRY)

      INTEGER          :: NTRAVERSE_FINE(MAXLAYERS,MAXFINELAYERS)
      DOUBLE PRECISION :: SUNPATHS_FINE (MAXLAYERS,MAXLAYERS,MAXFINELAYERS)
      DOUBLE PRECISION :: RADII_FINE    (MAXLAYERS,MAXFINELAYERS)
      DOUBLE PRECISION :: ALPHA_FINE    (MAXLAYERS,MAXFINELAYERS)

!  PARTIAL LAYER OUTPUT (GEOMETRY)

      INTEGER          :: NTRAVERSE_UT(MAX_PARTLAYERS)
      DOUBLE PRECISION :: SUNPATHS_UT (MAX_PARTLAYERS,MAXLAYERS)
      DOUBLE PRECISION :: RADII_UT    (MAX_PARTLAYERS)
      DOUBLE PRECISION :: ALPHA_UT    (MAX_PARTLAYERS)

!  OTHER (INCIDENTAL) GEOMETRICAL OUTPUT

      DOUBLE PRECISION :: LOSPATHS(MAXLAYERS)
      DOUBLE PRECISION :: LOSPATHS_UT_UP(MAX_PARTLAYERS)
      DOUBLE PRECISION :: LOSPATHS_UT_DN(MAX_PARTLAYERS)
      DOUBLE PRECISION :: THETA_ALL  (0:MAXLAYERS)
      DOUBLE PRECISION :: PHI_ALL    (0:MAXLAYERS)
      DOUBLE PRECISION :: COSSCAT_UP (0:MAXLAYERS)
      DOUBLE PRECISION :: COSSCAT_DN (0:MAXLAYERS)

!  EXTINCTION AND VAR TMS

      DOUBLE PRECISION :: EXTINCTION   (MAXLAYERS)
      DOUBLE PRECISION :: L_EXTINCTION (MAXLAYERS,MAX_ATMOSWFS)
      DOUBLE PRECISION :: VAR_TMS      (MAXLAYERS,MAX_ATMOSWFS)

!  PARTIAL LAYER HEIGHTS

      DOUBLE PRECISION :: HEIGHT_GRID_UT(MAX_PARTLAYERS)

!  BOOKKEEPING

      INTEGER ::          NTYPE_VARY
      LOGICAL ::          DO_RTSOL_VARY(MAXLAYERS)
      INTEGER ::          NPARAMS_VARY(MAXLAYERS)

!  LOCAL VARIABLES
!  ---------------

!  CUMULATIVE TRANSMITTANCES
!    NEW, 27 JANUARY 2010. TOA CONTRIBUTION FUNCTIONS CODE.

      DOUBLE PRECISION :: OGCTRANS(MAXLAYERS,MAX_GEOMETRIES)

!  INDICES

      INTEGER ::          N, NUT, NSTART, NUT_PREV, NLEVEL, LUM, LIB, LUA
      INTEGER ::          UT, UTA, UM, IA, NC, IB, V, O1, O2, NM1
      INTEGER ::          K, K_PARAMETERS, Q

!  HELP VARIABLES

      DOUBLE PRECISION :: FINAL_SOURCE, HELP, SS_CUMSOURCE
      DOUBLE PRECISION :: SS_LAYERSOURCE, SSCORRECTION, XT
      DOUBLE PRECISION :: CTHETA, STHETA, CALPHA, SALPHA, CPHI, CSA
      DOUBLE PRECISION :: VSIGN, ZMAT_LOCAL(4,4), FMAT_LOCAL(6), DNM1
      DOUBLE PRECISION :: L_ZMAT_LOCAL(4,4,MAX_ATMOSWFS), LMULT
      DOUBLE PRECISION :: UVAR, AVAR, FT1, FT2, LSS, TRANS
      DOUBLE PRECISION :: L_FINAL_SOURCE, L_SSCORRECTION

      CHARACTER (LEN=3) :: CV

      INTEGER ::          PARTLAYERS_LAYERFINEIDX ( MAX_PARTLAYERS )

      DOUBLE PRECISION ::  L_ZMAT_UP &
          ( MAX_ATMOSWFS, MAX_GEOMETRIES, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION ::  L_ZMAT_DN &
          ( MAX_ATMOSWFS, MAX_GEOMETRIES, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION ::  L_SS_CUMSOURCE &
          ( MAX_ATMOSWFS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION ::  L_SSFDEL ( MAXLAYERS, MAX_ATMOSWFS )

!  SUNLIGHT PARAMETERS
!    - ASSUMES THE FLUX VECTOR IS (F,0,0,0)
!   - DROP THIS VARIABLES WHEN FURTHER TESTING HAS BEEN DONE

      LOGICAL, PARAMETER :: SUNLIGHT = .TRUE.
      INTEGER ::            NPOLAR

!  SET UP OPERATIONS
!  -----------------

!  INITIALIZE EXCEPTION HANDLING

      FAIL    = .FALSE.
      MESSAGE = ' '
      TRACE   = ' '

!  LOCAL USER INDICES

      LUM = 1
      LIB = 1
      LUA = 1

!  SET UP PARTIALS FLAG

      DO_FINE = .TRUE.
      DO_PARTIALS = ( N_PARTLAYERS .GT. 0 )

!  SET NPOLAR. FOR NATURAL LIGHT, THERE ARE ONLY 3 STOKES PARAMETERS
!    ( NO CIRCULAR POLARIZATION FOR SINGLE SCATTERING)

      IF (SUNLIGHT) NPOLAR = MIN(NSTOKES,3)
      IF (.NOT.SUNLIGHT) NPOLAR = NSTOKES

!  CREATE TMS FACTORS, THESE GET STORED
!    DELTA-M SCALING INTRODUCED APRIL 2005.

      DO N = 1, NLAYERS
        IF ( DO_DELTAM_SCALING ) THEN
          HELP   = ONE - TRUNC_FACTOR(N) * OMEGA_TOTAL_INPUT(N)
          TMS(N) = OMEGA_TOTAL_INPUT(N) / HELP
        ELSE
          TMS(N) = OMEGA_TOTAL_INPUT(N)
        ENDIF
      ENDDO

!  LINEARIZATIONS OF THE TMS FACTORS
!  ---------------------------------

!  CREATE LINEARIZED TMS FACTOR FOR EACH LAYER
!   ( USE UNSCALED LINEARIZED INPUTS - THE NAKAJIMA-TANAKA WAY)
!    DISTINGUISH BETWEEN DELTAM CASE OR NOT.
!  ONLY IF REQUIRED

      IF ( DO_COLUMN_LINEARIZATION ) THEN

!  LAYER LOOP

       DO K = 1, NLAYERS
        IF ( LAYER_VARY_FLAG(K) ) THEN
         K_PARAMETERS = LAYER_VARY_NUMBER(K)

!  DELTAM_SCALING. EXTRA CONTRIBUTIONS IF PHASE FUNCTION MOMENTS ARE VAR

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

!  NO DELTA-M SCALING, JUST COPY

         IF ( .NOT. DO_DELTAM_SCALING ) THEN
          DO Q = 1, K_PARAMETERS
           UVAR = L_OMEGA_TOTAL_INPUT(Q,K)
           VAR_TMS(K,Q) = UVAR
          ENDDO
         ENDIF

!  END LAYER LOOP

        ENDIF
       ENDDO

!  END LINEARIZATION CLAUSE

      ENDIF

!  ADDITIONAL DELTA-M SCALING
!  --------------------------

!  NEW SECTION. R. SPURR, 07 SEPTEMBER 2007.
!   TMS GETS MODIFIED BY (1-F). SAVE THE TRUNCATION FACTOR.
!   PHASE FUNCTION MOMENTS ARE MODIFIED LATER ON.

      IF ( DO_SSCORR_TRUNCATION ) THEN

!  BASIC SCALING

        NM1  = NGREEK_MOMENTS_INPUT
        DNM1 = DBLE(2*NM1+1)
        DO N = 1, NLAYERS
          SSFDEL(N) = GREEKMAT_TOTAL_INPUT(NM1,N,1) / DNM1
          TMS(N) = TMS(N) * ( ONE - SSFDEL(N) )
        ENDDO

!  LINEARIZATION, ONLY CHANGE IF SCATMAT VARIATION IS SET

        IF ( DO_COLUMN_LINEARIZATION ) THEN
          DO K = 1, NLAYERS
            IF ( LAYER_VARY_FLAG(K) ) THEN
              K_PARAMETERS = LAYER_VARY_NUMBER(K)
              DO Q = 1, K_PARAMETERS
                IF ( DO_SCATMAT_VARIATION(K,Q) ) THEN
                  FT1 = ONE / ( ONE- SSFDEL(K) )
                  L_SSFDEL(K,Q) = GREEKMAT_TOTAL_INPUT(NM1,K,1) * &
                          L_GREEKMAT_TOTAL_INPUT(Q,NM1,K,1) / DNM1
                  VAR_TMS(K,Q) = VAR_TMS(K,Q) - FT1 * L_SSFDEL(K,Q)
                ENDIF
              ENDDO
            ENDIF
          ENDDO
        ENDIF

      ENDIF

!  PARAMETER CONTROL
!  -----------------

      IF ( DO_COLUMN_LINEARIZATION ) THEN
        NTYPE_VARY = 1
        DO_RTSOL_VARY(1) = .TRUE.
        NPARAMS_VARY(1)  = N_TOTALCOLUMN_WFS
      ENDIF

!  CREATE EXTINCTIONS
!  ------------------

!  USE BASIC DEFINITIONS

      DO N = 1, NLAYERS
        HELP = HEIGHT_GRID(N-1) - HEIGHT_GRID(N)
        EXTINCTION(N) = DELTAU_VERT(N) / HELP
      ENDDO

!  LINEARIZED EXTINCTIONS

      IF ( DO_COLUMN_LINEARIZATION ) THEN
       DO K = 1, NLAYERS
        IF ( LAYER_VARY_FLAG(K) ) THEN
         K_PARAMETERS = LAYER_VARY_NUMBER(K)
         DO Q = 1, K_PARAMETERS
           L_EXTINCTION(K,Q) = L_DELTAU_VERT(Q,K) * EXTINCTION(K)
         ENDDO
        ENDIF
       ENDDO
      ENDIF

!  CREATE THE HEIGHTS OF PARTIAL LAYERS

      IF ( DO_PARTIALS ) THEN
        DO UT = 1, N_PARTLAYERS
          N = PARTLAYERS_LAYERIDX(UT)
          XT = DELTAU_VERT(N) - PARTAU_VERT(UT)
          HEIGHT_GRID_UT(UT) = HEIGHT_GRID(N) + XT / EXTINCTION(N)
        ENDDO
      ENDIF

!  COMPUTE AND STORE ALL GEOMETRY-DEPENDENT STUFF
!  ==============================================

!  Lattice calculation

      IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN

!  START THE MAIN LOOP OVER ALL SOLAR AND VIEWING GEOMETRIES
!   OLD CODE, NOT ADJUSTED
!      DO IB = 1, NBEAMS
!        THETA_BOA = SZANGLES(IB)
!        DO UM = 1, N_USER_STREAMS
!          ALPHA_BOA = USER_VZANGLES_INPUT(UM)
!          DO IA = 1, N_USER_RELAZMS
!            PHI_BOA = USER_RELAZMS(IA)
!            V = VZA_OFFSETS(IB,UM) + IA

!  START THE MAIN LOOP OVER ALL SOLAR AND VIEWING GEOMETRIES

        DO UM = 1, N_USER_STREAMS
          ALPHA_BOA = USER_VZANGLES_ADJUST(UM)
          DO IB = 1, NBEAMS
            DO IA = 1, N_USER_RELAZMS
              THETA_BOA = SZANGLES_ADJUST(UM,IB,IA)
              PHI_BOA   = USER_RELAZMS_ADJUST(UM,IB,IA)
              V = VZA_OFFSETS(IB,UM) + IA

!  UPWELLING CALCULATION
!  ---------------------

              IF ( DO_UPWELLING ) THEN

!  CALL TO GEOMETRY

                CALL OUTGOING_SPHERGEOM_FINE_UP ( &
                  MAXLAYERS, MAXFINELAYERS, MAX_PARTLAYERS, &
                  DO_FINE, DO_PARTIALS, NLAYERS, NFINELAYERS, &
                  N_PARTLAYERS, PARTLAYERS_LAYERIDX, &
                  HEIGHT_GRID, HEIGHT_GRID_UT, EARTH_RADIUS, &
                  ALPHA_BOA, THETA_BOA, PHI_BOA, &
                  SUNPATHS,      RADII,      NTRAVERSE,      ALPHA_ALL, &
                  SUNPATHS_FINE, RADII_FINE, NTRAVERSE_FINE, ALPHA_FINE, &
                  SUNPATHS_UT,   RADII_UT,   NTRAVERSE_UT,   ALPHA_UT, &
                  PARTLAYERS_LAYERFINEIDX, LOSPATHS, LOSPATHS_UT_UP, &
                  THETA_ALL, PHI_ALL, COSSCAT_UP, &
                  FAIL, MESSAGE )

!  EXCEPTION HANDLING

                IF ( FAIL ) THEN
                  WRITE(CV,'(I3)')V
                  TRACE = 'ERROR FROM OUTGOING_SPHERGEOM_FINE_UP, '// &
                          ' GEOMETRY '//CV
                  RETURN
                ENDIF

!  DEBUG GEOMETRY
!            DO N = 0, NLAYERS
!              WRITE(45,'(I4,101F10.5)')N,(SUNPATHS(N,V),V=1,NLAYERS)
!            ENDDO
!            PAUSE
!  MULTIPLIERS AND TRANSMITTANCES + LINEARIZATIONS

                CALL LC_OUTGOING_INTEGRATION_UP ( &
                  NLAYERS, NFINELAYERS, &
                  DO_PARTIALS, &
                  LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
                  DO_COLUMN_LINEARIZATION, N_TOTALCOLUMN_WFS, &
                  EXTINCTION, L_EXTINCTION, &
                  N_PARTLAYERS, PARTLAYERS_LAYERIDX, &
                  PARTLAYERS_LAYERFINEIDX, &
                  SUNPATHS, RADII, NTRAVERSE, ALPHA_ALL, &
                  SUNPATHS_UT, RADII_UT,  NTRAVERSE_UT, ALPHA_UT, &
                  SUNPATHS_FINE, NTRAVERSE_FINE, ALPHA_FINE, &
                  UP_MULTIPLIERS(1,V), UP_LOSTRANS(1,V), BOA_ATTN(V), &
                  UP_MULTIPLIERS_UT(1,V), UP_LOSTRANS_UT(1,V), &
                  LC_UP_MULTIPLIERS(1,1,V), &
                  L_UP_LOSTRANS(1,1,V), &
                  LC_BOA_ATTN(1,V), &
                  LC_UP_MULTIPLIERS_UT(1,1,V), &
                  L_UP_LOSTRANS_UT(1,1,V) )

!  DEBUG
!                IF (V.EQ.1)WRITE(35,*)'LIN',V
!                DO UT = 1, N_PARTLAYERS
!                  IF (V.EQ.1)WRITE(35,'(1P2E18.10)')
!     &             UP_LOSTRANS_UT(UT,V),UP_MULTIPLIERS_UT(UT,V)
!                ENDDO

!                WRITE(35,*)V
!                DO N = 1, NLAYERS
!                  WRITE(35,'(1P2E18.10)')
!     &         UP_LOSTRANS(N,V),UP_MULTIPLIERS(N,V)
!                ENDDO
!                IF ( V.EQ.8)PAUSE

!                IF ( DO_COLUMN_LINEARIZATION ) THEN
!                 DO N = 1, NLAYERS
!                   WRITE(36,'(2I4,1P20E18.10)')V,N,
!     &             UP_LOSTRANS(N,V),UP_MULTIPLIERS(N,V),
!     &             L_UP_LOSTRANS(N,1,V),(LC_UP_MULTIPLIERS(N,1,V),K=1,NLAYERS)
!                 ENDDO
!                ELSE
!                 DO N = 1, NLAYERS
!                   WRITE(35,'(2I4,1P2E18.10)')V,N,
!     &             UP_LOSTRANS(N,V),UP_MULTIPLIERS(N,V)
!                 ENDDO
!                ENDIF

!  START LAYER LOOP

                DO N = NLAYERS, 1, -1

!  TRIGONOMETRY

                  CTHETA = DCOS(THETA_ALL(N))
                  STHETA = DSIN(THETA_ALL(N))
                  CALPHA = DCOS(ALPHA_ALL(N))
                  SALPHA = DSIN(ALPHA_ALL(N))
                  CPHI   = DCOS(PHI_ALL(N))
                  CSA    = COSSCAT_UP(N)
                  VSIGN  = -1.0D0

!  GET THE SCATTERING MATRIX AND ITS LINEARIZATION

                  CALL SSCORR_OUTGOING_L_ZMATRIX ( &
                    DO_SSCORR_TRUNCATION, SUNLIGHT, DO_COLUMN_LINEARIZATION, &
                    N, NSTOKES, NGREEK_MOMENTS_INPUT, &
                    LAYER_VARY_FLAG(N), LAYER_VARY_NUMBER(N), &
                    LAYER_MAXMOMENTS(N),  GREEKMAT_TOTAL_INPUT,   SSFDEL, &
                    DO_SCATMAT_VARIATION, L_GREEKMAT_TOTAL_INPUT, L_SSFDEL, &
                    CTHETA, STHETA, CALPHA, SALPHA, CPHI, &
                    PHI_ALL(N), CSA, VSIGN, &
                    ZMAT_LOCAL, L_ZMAT_LOCAL, FMAT_LOCAL )

!  PHASE MATRIX + LINEARIZED (MULTIPLIED BY TMS FACTOR). SAVE THEM.
!     SUNLIGHT ONLY, THE FIRST COLUMN OF THE MATRIX
!    GENERAL CASE, CODE PROGRAMMED 05 OCTOBER 2010

                  IF ( STERM_LAYERMASK_UP(N) ) THEN
                   DO O1 = 1, NSTOKES
                    DO O2 = 1, NSTOKES
                     ZMAT_UP(V,N,O1,O2) = ZMAT_LOCAL(O1,O2) * TMS(N)
                    ENDDO
                   ENDDO
                   IF ( DO_COLUMN_LINEARIZATION ) THEN
                    IF ( LAYER_VARY_FLAG(N) ) THEN
                     DO Q = 1, LAYER_VARY_NUMBER(N)
                      DO O1 = 1, NSTOKES
                       DO O2 = 1, NSTOKES
                        L_ZMAT_UP(Q,V,N,O1,O2) = &
                                ZMAT_UP(V,N,O1,O2)   * VAR_TMS(N,Q) + &
                              L_ZMAT_LOCAL(O1,O2,Q)  *     TMS(N)
                       ENDDO
                      ENDDO
                     ENDDO
                    ENDIF
                   ENDIF
                  ENDIF

!  FINISH THE LAYER LOOP

                ENDDO

!  END UPWELLING CLAUSE

              ENDIF

!  DOWNWELLING CALCULATION
!  -----------------------

              IF ( DO_DNWELLING ) THEN

!  CALL TO GEOMETRY

                CALL OUTGOING_SPHERGEOM_FINE_DN ( &
                  MAXLAYERS, MAXFINELAYERS, MAX_PARTLAYERS, &
                  DO_FINE, DO_PARTIALS, NLAYERS, NFINELAYERS, &
                  N_PARTLAYERS, PARTLAYERS_LAYERIDX, &
                  HEIGHT_GRID, HEIGHT_GRID_UT, EARTH_RADIUS, &
                  ALPHA_BOA, THETA_BOA, PHI_BOA, &
                  SUNPATHS,      RADII,      NTRAVERSE,      ALPHA_ALL, &
                  SUNPATHS_FINE, RADII_FINE, NTRAVERSE_FINE, ALPHA_FINE, &
                  SUNPATHS_UT,   RADII_UT,   NTRAVERSE_UT,   ALPHA_UT, &
                  PARTLAYERS_LAYERFINEIDX, LOSPATHS, LOSPATHS_UT_DN, &
                  THETA_ALL, PHI_ALL, COSSCAT_DN, &
                  FAIL, MESSAGE )

!  EXCEPTION HANDLING

                IF ( FAIL ) THEN
                  WRITE(CV,'(I3)')V
                  TRACE = 'ERROR FROM OUTGOING_SPHERGEOM_FINE_DN, '// &
                          ' GEOMETRY '//CV
                  RETURN
                ENDIF

!  MULTIPLIERS, TRANSMITTANCES + LINEARIZATIONS

                CALL LC_OUTGOING_INTEGRATION_DN ( &
                  NLAYERS, NFINELAYERS, DO_PARTIALS, &
                  LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
                  DO_COLUMN_LINEARIZATION, N_TOTALCOLUMN_WFS, &
                  EXTINCTION, L_EXTINCTION, &
                  N_PARTLAYERS, PARTLAYERS_LAYERIDX, &
                  PARTLAYERS_LAYERFINEIDX, &
                  SUNPATHS, RADII, NTRAVERSE, ALPHA_ALL, &
                  SUNPATHS_UT,   NTRAVERSE_UT,   ALPHA_UT, &
                  SUNPATHS_FINE, NTRAVERSE_FINE, ALPHA_FINE, &
                  DN_MULTIPLIERS(1,V), DN_LOSTRANS(1,V), &
                  DN_MULTIPLIERS_UT(1,V), DN_LOSTRANS_UT(1,V), &
                  LC_DN_MULTIPLIERS(1,1,V), &
                  L_DN_LOSTRANS(1,1,V), &
                  LC_DN_MULTIPLIERS_UT(1,1,V), &
                  L_DN_LOSTRANS_UT(1,1,V) )

!  DEBUG
!                DO N = 1, NLAYERS
!                  WRITE(88,*)N,DN_LOSTRANS(N,1),DN_MULTIPLIERS(N,1)
!                ENDDO

!  LAYER LOOP

                DO N = 1, NLAYERS

!  TRIGONOMETRY

                  CTHETA = DCOS(THETA_ALL(N-1))
                  STHETA = DSIN(THETA_ALL(N-1))
                  CALPHA = DCOS(ALPHA_ALL(N-1))
                  SALPHA = DSIN(ALPHA_ALL(N-1))
                  CPHI   = DCOS(PHI_ALL(N-1))
                  CSA    = COSSCAT_DN(N-1)
                  VSIGN  = +1.0D0

!  GET THE SCATTERING MATRIX AND ITS LINEARIZATION

                  CALL SSCORR_OUTGOING_L_ZMATRIX ( &
                    DO_SSCORR_TRUNCATION, SUNLIGHT, DO_COLUMN_LINEARIZATION, &
                    N, NSTOKES, NGREEK_MOMENTS_INPUT, &
                    LAYER_VARY_FLAG(N), LAYER_VARY_NUMBER(N), &
                    LAYER_MAXMOMENTS(N),  GREEKMAT_TOTAL_INPUT,   SSFDEL, &
                    DO_SCATMAT_VARIATION, L_GREEKMAT_TOTAL_INPUT, L_SSFDEL, &
                    CTHETA, STHETA, CALPHA, SALPHA, CPHI, &
                    PHI_ALL(N-1), CSA, VSIGN, &
                    ZMAT_LOCAL, L_ZMAT_LOCAL, FMAT_LOCAL )

!  PHASE MATRIX + LINEARIZATION (MULTIPLIED BY TMS FACTOR). SAVE THEM.
!     SUNLIGHT ONLY, THE FIRST COLUMN OF THE MATRIX
!    GENERAL CASE, CODE PROGRAMMED 05 OCTOBER 2010

                  IF ( STERM_LAYERMASK_DN(N) ) THEN
                   DO O1 = 1, NSTOKES
                    DO O2 = 1, NSTOKES
                     ZMAT_DN(V,N,O1,O2) = ZMAT_LOCAL(O1,O2) * TMS(N)
                    ENDDO
                   ENDDO
                   IF ( DO_COLUMN_LINEARIZATION ) THEN
                    IF ( LAYER_VARY_FLAG(N) ) THEN
                     DO Q = 1, LAYER_VARY_NUMBER(N)
                      DO O1 = 1, NSTOKES
                       DO O2 = 1, NSTOKES
                        L_ZMAT_DN(Q,V,N,O1,O2) = &
                                ZMAT_DN(V,N,O1,O2)   * VAR_TMS(N,Q) + &
                              L_ZMAT_LOCAL(O1,O2,Q)  *     TMS(N)
                       ENDDO
                      ENDDO
                     ENDDO
                    ENDIF
                   ENDIF
                  ENDIF

!                  IF ( N.GT.23) WRITE(*,*)V,N,L_ZMAT_DN(1,V,N,1,1)

!  FINISH THE LAYER LOOP

                ENDDO

!  END DOWNWELLING CLAUSE

              ENDIF

!   FINISH GEOMETRY LOOPS

            ENDDO
          ENDDO
        ENDDO

!  End lattice calculation

      ENDIF

!  Observational geometry calculation

      IF ( DO_OBSERVATION_GEOMETRY ) THEN

!  START THE MAIN LOOP OVER ALL SOLAR AND VIEWING GEOMETRIES
        DO V = 1, N_GEOMETRIES
          ALPHA_BOA = USER_VZANGLES_ADJUST(V)
          THETA_BOA = SZANGLES_ADJUST(LUM,V,LUA)
          PHI_BOA   = USER_RELAZMS_ADJUST(LUM,LIB,V)

!  UPWELLING CALCULATION
!  ---------------------

          IF ( DO_UPWELLING ) THEN

!  CALL TO GEOMETRY

            CALL OUTGOING_SPHERGEOM_FINE_UP ( &
              MAXLAYERS, MAXFINELAYERS, MAX_PARTLAYERS, &
              DO_FINE, DO_PARTIALS, NLAYERS, NFINELAYERS, &
              N_PARTLAYERS, PARTLAYERS_LAYERIDX, &
              HEIGHT_GRID, HEIGHT_GRID_UT, EARTH_RADIUS, &
              ALPHA_BOA, THETA_BOA, PHI_BOA, &
              SUNPATHS,      RADII,      NTRAVERSE,      ALPHA_ALL, &
              SUNPATHS_FINE, RADII_FINE, NTRAVERSE_FINE, ALPHA_FINE, &
              SUNPATHS_UT,   RADII_UT,   NTRAVERSE_UT,   ALPHA_UT, &
              PARTLAYERS_LAYERFINEIDX, LOSPATHS, LOSPATHS_UT_UP, &
              THETA_ALL, PHI_ALL, COSSCAT_UP, &
              FAIL, MESSAGE )

!  EXCEPTION HANDLING

            IF ( FAIL ) THEN
              WRITE(CV,'(I3)')V
              TRACE = 'ERROR FROM OUTGOING_SPHERGEOM_FINE_UP, '// &
                      ' GEOMETRY '//CV
              RETURN
            ENDIF

!  DEBUG GEOMETRY
!            DO N = 0, NLAYERS
!              WRITE(45,'(I4,101F10.5)')N,(SUNPATHS(N,V),V=1,NLAYERS)
!            ENDDO
!            PAUSE
!  MULTIPLIERS AND TRANSMITTANCES + LINEARIZATIONS

            CALL LC_OUTGOING_INTEGRATION_UP ( &
              NLAYERS, NFINELAYERS, &
              DO_PARTIALS, &
              LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
              DO_COLUMN_LINEARIZATION, N_TOTALCOLUMN_WFS, &
              EXTINCTION, L_EXTINCTION, &
              N_PARTLAYERS, PARTLAYERS_LAYERIDX, &
              PARTLAYERS_LAYERFINEIDX, &
              SUNPATHS, RADII, NTRAVERSE, ALPHA_ALL, &
              SUNPATHS_UT, RADII_UT,  NTRAVERSE_UT, ALPHA_UT, &
              SUNPATHS_FINE, NTRAVERSE_FINE, ALPHA_FINE, &
              UP_MULTIPLIERS(1,V), UP_LOSTRANS(1,V), BOA_ATTN(V), &
              UP_MULTIPLIERS_UT(1,V), UP_LOSTRANS_UT(1,V), &
              LC_UP_MULTIPLIERS(1,1,V), &
              L_UP_LOSTRANS(1,1,V), &
              LC_BOA_ATTN(1,V), &
              LC_UP_MULTIPLIERS_UT(1,1,V), &
              L_UP_LOSTRANS_UT(1,1,V) )

!  DEBUG
!              IF (V.EQ.1)WRITE(35,*)'LIN',V
!              DO UT = 1, N_PARTLAYERS
!                IF (V.EQ.1)WRITE(35,'(1P2E18.10)')
!     &           UP_LOSTRANS_UT(UT,V),UP_MULTIPLIERS_UT(UT,V)
!              ENDDO

!              WRITE(35,*)V
!              DO N = 1, NLAYERS
!                WRITE(35,'(1P2E18.10)')
!     &       UP_LOSTRANS(N,V),UP_MULTIPLIERS(N,V)
!              ENDDO
!              IF ( V.EQ.8)PAUSE

!              IF ( DO_COLUMN_LINEARIZATION ) THEN
!              DO N = 1, NLAYERS
!                WRITE(36,'(2I4,1P20E18.10)')V,N,
!     &       UP_LOSTRANS(N,V),UP_MULTIPLIERS(N,V),
!     &L_UP_LOSTRANS(N,1,V),(LC_UP_MULTIPLIERS(N,1,V),K=1,NLAYERS)
!              ENDDO
!               ELSE
!              DO N = 1, NLAYERS
!                WRITE(35,'(2I4,1P2E18.10)')V,N,
!     &       UP_LOSTRANS(N,V),UP_MULTIPLIERS(N,V)
!                ENDDO
!               ENDIF

!  START LAYER LOOP

            DO N = NLAYERS, 1, -1

!  TRIGONOMETRY

              CTHETA = DCOS(THETA_ALL(N))
              STHETA = DSIN(THETA_ALL(N))
              CALPHA = DCOS(ALPHA_ALL(N))
              SALPHA = DSIN(ALPHA_ALL(N))
              CPHI   = DCOS(PHI_ALL(N))
              CSA    = COSSCAT_UP(N)
              VSIGN  = -1.0D0

!  GET THE SCATTERING MATRIX AND ITS LINEARIZATION

              CALL SSCORR_OUTGOING_L_ZMATRIX ( &
                DO_SSCORR_TRUNCATION, SUNLIGHT, DO_COLUMN_LINEARIZATION, &
                N, NSTOKES, NGREEK_MOMENTS_INPUT, &
                LAYER_VARY_FLAG(N), LAYER_VARY_NUMBER(N), &
                LAYER_MAXMOMENTS(N),  GREEKMAT_TOTAL_INPUT,   SSFDEL, &
                DO_SCATMAT_VARIATION, L_GREEKMAT_TOTAL_INPUT, L_SSFDEL, &
                CTHETA, STHETA, CALPHA, SALPHA, CPHI, &
                PHI_ALL(N), CSA, VSIGN, &
                ZMAT_LOCAL, L_ZMAT_LOCAL, FMAT_LOCAL )

!  PHASE MATRIX + LINEARIZED (MULTIPLIED BY TMS FACTOR). SAVE THEM.
!     SUNLIGHT ONLY, THE FIRST COLUMN OF THE MATRIX
!    GENERAL CASE, CODE PROGRAMMED 05 OCTOBER 2010

              IF ( STERM_LAYERMASK_UP(N) ) THEN
               DO O1 = 1, NSTOKES
                DO O2 = 1, NSTOKES
                 ZMAT_UP(V,N,O1,O2) = ZMAT_LOCAL(O1,O2) * TMS(N)
                ENDDO
               ENDDO
               IF ( DO_COLUMN_LINEARIZATION ) THEN
                IF ( LAYER_VARY_FLAG(N) ) THEN
                 DO Q = 1, LAYER_VARY_NUMBER(N)
                  DO O1 = 1, NSTOKES
                   DO O2 = 1, NSTOKES
                    L_ZMAT_UP(Q,V,N,O1,O2) = &
                            ZMAT_UP(V,N,O1,O2)   * VAR_TMS(N,Q) + &
                          L_ZMAT_LOCAL(O1,O2,Q)  *     TMS(N)
                   ENDDO
                  ENDDO
                 ENDDO
                ENDIF
               ENDIF
              ENDIF

!  FINISH THE LAYER LOOP

            ENDDO

!  END UPWELLING CLAUSE

          ENDIF

!  DOWNWELLING CALCULATION
!  -----------------------

          IF ( DO_DNWELLING ) THEN

!  CALL TO GEOMETRY

            CALL OUTGOING_SPHERGEOM_FINE_DN ( &
              MAXLAYERS, MAXFINELAYERS, MAX_PARTLAYERS, &
              DO_FINE, DO_PARTIALS, NLAYERS, NFINELAYERS, &
              N_PARTLAYERS, PARTLAYERS_LAYERIDX, &
              HEIGHT_GRID, HEIGHT_GRID_UT, EARTH_RADIUS, &
              ALPHA_BOA, THETA_BOA, PHI_BOA, &
              SUNPATHS,      RADII,      NTRAVERSE,      ALPHA_ALL, &
              SUNPATHS_FINE, RADII_FINE, NTRAVERSE_FINE, ALPHA_FINE, &
              SUNPATHS_UT,   RADII_UT,   NTRAVERSE_UT,   ALPHA_UT, &
              PARTLAYERS_LAYERFINEIDX, LOSPATHS, LOSPATHS_UT_DN, &
              THETA_ALL, PHI_ALL, COSSCAT_DN, &
              FAIL, MESSAGE )

!  EXCEPTION HANDLING

            IF ( FAIL ) THEN
              WRITE(CV,'(I3)')V
              TRACE = 'ERROR FROM OUTGOING_SPHERGEOM_FINE_DN, '// &
                      ' GEOMETRY '//CV
              RETURN
            ENDIF

!  MULTIPLIERS, TRANSMITTANCES + LINEARIZATIONS

            CALL LC_OUTGOING_INTEGRATION_DN ( &
              NLAYERS, NFINELAYERS, DO_PARTIALS, &
              LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
              DO_COLUMN_LINEARIZATION, N_TOTALCOLUMN_WFS, &
              EXTINCTION, L_EXTINCTION, &
              N_PARTLAYERS, PARTLAYERS_LAYERIDX, &
              PARTLAYERS_LAYERFINEIDX, &
              SUNPATHS, RADII, NTRAVERSE, ALPHA_ALL, &
              SUNPATHS_UT,   NTRAVERSE_UT,   ALPHA_UT, &
              SUNPATHS_FINE, NTRAVERSE_FINE, ALPHA_FINE, &
              DN_MULTIPLIERS(1,V), DN_LOSTRANS(1,V), &
              DN_MULTIPLIERS_UT(1,V), DN_LOSTRANS_UT(1,V), &
              LC_DN_MULTIPLIERS(1,1,V), &
              L_DN_LOSTRANS(1,1,V), &
              LC_DN_MULTIPLIERS_UT(1,1,V), &
              L_DN_LOSTRANS_UT(1,1,V) )

!  DEBUG
!              DO N = 1, NLAYERS
!                WRITE(88,*)N,DN_LOSTRANS(N,1),DN_MULTIPLIERS(N,1)
!              ENDDO

!  LAYER LOOP

            DO N = 1, NLAYERS

!  TRIGONOMETRY

              CTHETA = DCOS(THETA_ALL(N-1))
              STHETA = DSIN(THETA_ALL(N-1))
              CALPHA = DCOS(ALPHA_ALL(N-1))
              SALPHA = DSIN(ALPHA_ALL(N-1))
              CPHI   = DCOS(PHI_ALL(N-1))
              CSA    = COSSCAT_DN(N-1)
              VSIGN  = +1.0D0

!  GET THE SCATTERING MATRIX AND ITS LINEARIZATION

              CALL SSCORR_OUTGOING_L_ZMATRIX ( &
                DO_SSCORR_TRUNCATION, SUNLIGHT, DO_COLUMN_LINEARIZATION, &
                N, NSTOKES, NGREEK_MOMENTS_INPUT, &
                LAYER_VARY_FLAG(N), LAYER_VARY_NUMBER(N), &
                LAYER_MAXMOMENTS(N),  GREEKMAT_TOTAL_INPUT,   SSFDEL, &
                DO_SCATMAT_VARIATION, L_GREEKMAT_TOTAL_INPUT, L_SSFDEL, &
                CTHETA, STHETA, CALPHA, SALPHA, CPHI, &
                PHI_ALL(N-1), CSA, VSIGN, &
                ZMAT_LOCAL, L_ZMAT_LOCAL, FMAT_LOCAL )

!  PHASE MATRIX + LINEARIZATION (MULTIPLIED BY TMS FACTOR). SAVE THEM.
!     SUNLIGHT ONLY, THE FIRST COLUMN OF THE MATRIX
!    GENERAL CASE, CODE PROGRAMMED 05 OCTOBER 2010

              IF ( STERM_LAYERMASK_DN(N) ) THEN
               DO O1 = 1, NSTOKES
                DO O2 = 1, NSTOKES
                 ZMAT_DN(V,N,O1,O2) = ZMAT_LOCAL(O1,O2) * TMS(N)
                ENDDO
               ENDDO
               IF ( DO_COLUMN_LINEARIZATION ) THEN
                IF ( LAYER_VARY_FLAG(N) ) THEN
                 DO Q = 1, LAYER_VARY_NUMBER(N)
                  DO O1 = 1, NSTOKES
                   DO O2 = 1, NSTOKES
                    L_ZMAT_DN(Q,V,N,O1,O2) = &
                            ZMAT_DN(V,N,O1,O2)   * VAR_TMS(N,Q) + &
                          L_ZMAT_LOCAL(O1,O2,Q)  *     TMS(N)
                   ENDDO
                  ENDDO
                 ENDDO
                ENDIF
               ENDIF
              ENDIF

!              IF ( N.GT.23) WRITE(*,*)V,N,L_ZMAT_DN(1,V,N,1,1)

!  FINISH THE LAYER LOOP

            ENDDO

!  END DOWNWELLING CLAUSE

          ENDIF

!   FINISH GEOMETRY LOOP

        ENDDO

!  End observationl geometry calculation

      ENDIF

!  RECURRENCE RELATION FOR THE UPWELLING STOKES VECTOR
!  ===================================================

      IF ( DO_UPWELLING ) THEN

!  CUMULATIVE TRANMSITTANCES (TOA CONTRIBUTION FUNCTIONS)
!     NEW SECTION, 27 JANAURY 2010

        IF ( DO_TOA_CONTRIBS ) THEN
          DO V = 1, N_GEOMETRIES
            OGCTRANS(1,V) = ONE
            DO N = 1, NLAYERS - 1
              OGCTRANS(N+1,V) = OGCTRANS(N,V) * UP_LOSTRANS(N,V)
            ENDDO
          ENDDO
        ENDIF

!  INITIALIZE CUMULATIVE SOURCE TERM

        NC =  0
        DO V = 1, N_GEOMETRIES
          DO O1 = 1, NSTOKES
            SS_CUMSOURCE_UP(V,O1,NC) = ZERO
          ENDDO
        ENDDO

!  INITIALISE OPTICAL DEPTH LOOP

        NSTART = NLAYERS
        NUT_PREV = NSTART + 1

!  MAIN LOOP OVER ALL OUTPUT OPTICAL DEPTHS

        DO UTA = N_USER_LEVELS, 1, -1

!  LAYER INDEX FOR GIVEN OPTICAL DEPTH

          NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
          NUT    = NLEVEL + 1

!  CUMULATIVE SINGLE SCATTER SOURCE TERMS :
!      FOR LOOP OVER LAYERS WORKING UPWARDS TO LEVEL NUT,
!      GET LAYER SOURCE TERMS = EXACT Z-MATRIX * MULTIPLIER
!      PROGRAMMED FOR GENERAL CASE (INCLUDING SUNLIGHT)

          DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            DO V = 1, N_GEOMETRIES
              DO O1 = 1, NPOLAR
                HELP = ZERO
                DO O2 = 1, NSTOKES
                 HELP = HELP + ZMAT_UP(V,N,O1,O2) * FLUXVEC(O2)
                ENDDO
                SS_LAYERSOURCE = HELP* UP_MULTIPLIERS(N,V)
                IF (DO_TOA_CONTRIBS) THEN
                  SS_CONTRIBS(V,O1,N) = OGCTRANS(N,V) * SS_LAYERSOURCE
                  SS_CONTRIBS(V,O1,N) = SSFLUX * SS_CONTRIBS(V,O1,N)
                ENDIF
                SS_CUMSOURCE_UP(V,O1,NC) = SS_LAYERSOURCE + &
                   UP_LOSTRANS(N,V)*SS_CUMSOURCE_UP(V,O1,NC-1)
              ENDDO
            ENDDO
          ENDDO

!  OFFGRID OUTPUT-------
!    ADD ADDITIONAL PARTIAL LAYER SOURCE TERM = EXACT PHASE FUNC * MULTI
!    SET FINAL CUMULATIVE SOURCE AND SINGLE SCATTER STOKES VECTOR
!  ONGRID OUTPUT--------
!    SET FINAL CUMULATIVE SOURCE AND SINGLE SCATTER STOKES VECTOR
!      PROGRAMMED FOR GENERAL CASE (INCLUDING SUNLIGHT)

          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)
            DO V = 1, N_GEOMETRIES
              DO O1 = 1, NPOLAR
                HELP = ZERO
                DO O2 = 1, NSTOKES
                 HELP = HELP + ZMAT_UP(V,N,O1,O2) * FLUXVEC(O2)
                ENDDO
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

!  CHECK FOR UPDATING THE RECURSION

          IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT

!  END OPTICAL DEPTH LOOP AND UPWELLING CLAUSE

        ENDDO
      ENDIF

!  RECURRENCE RELATION FOR THE DOWNWELLING STOKES VECTOR
!  =====================================================

      IF ( DO_DNWELLING ) THEN

!  INITIALIZE CUMULATIVE SOURCE TERM

        NC =  0
        DO V = 1, N_GEOMETRIES
          DO O1 = 1, NSTOKES
            SS_CUMSOURCE_DN(V,O1,NC) = ZERO
          ENDDO
        ENDDO

!  INITIALISE OPTICAL DEPTH LOOP

        NSTART = 1
        NUT_PREV = NSTART - 1

!  MAIN LOOP OVER ALL OUTPUT OPTICAL DEPTHS

        DO UTA = 1, N_USER_LEVELS

!  LAYER INDEX FOR GIVEN OPTICAL DEPTH

          NLEVEL = UTAU_LEVEL_MASK_DN(UTA)
          NUT = NLEVEL

!  CUMULATIVE SINGLE SCATTER SOURCE TERMS :
!      FOR LOOP OVER LAYERS WORKING DOWNWARDS TO NUT,
!      GET LAYER SOURCE TERMS = EXACT Z-MATRIX * MULTIPLIER
!      SUNLIGHT CASE ONLY

          DO N = NSTART, NUT
            NC = N
            DO V = 1, N_GEOMETRIES
              DO O1 = 1, NPOLAR
                HELP = ZERO
                DO O2 = 1, NSTOKES
                 HELP = HELP + ZMAT_DN(V,N,O1,O2) * FLUXVEC(O2)
                ENDDO
                SS_LAYERSOURCE = HELP* DN_MULTIPLIERS(N,V)
                SS_CUMSOURCE_DN(V,O1,NC) = SS_LAYERSOURCE + &
                   DN_LOSTRANS(N,V)*SS_CUMSOURCE_DN(V,O1,NC-1)
              ENDDO
            ENDDO
          ENDDO

!    SUNLIGHT CASE
!  OFFGRID OUTPUT :
!    ADD ADDITIONAL PARTIAL LAYER SOURCE TERM = EXACT Z-MATRIX * MULTIPL
!    SET FINAL CUMULATIVE SOURCE AND CORRECT STOKES VECTOR
!  ONGRID OUTPUT :
!    SET FINAL CUMULATIVE SOURCE AND CORRECT STOKES VECTOR

          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)
            DO V = 1, N_GEOMETRIES
              DO O1 = 1, NPOLAR
                HELP = ZERO
                DO O2 = 1, NSTOKES
                 HELP = HELP + ZMAT_DN(V,N,O1,O2) * FLUXVEC(O2)
                ENDDO
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

!  CHECK FOR UPDATING THE RECURSION

          IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
          NUT_PREV = NUT

!  END OPTICAL DEPTH LOOP AND DOWNWELLING CLAUSE

        ENDDO
      ENDIF

!  RECURRENCE RELATION FOR THE UPWELLING JACOBIANS
!  ===============================================

      IF ( DO_UPWELLING .AND. DO_COLUMN_LINEARIZATION ) THEN

!  START THE MAIN LAYER VARIATION LOOP
!  -----------------------------------

        DO K = 1, NTYPE_VARY
         IF ( DO_RTSOL_VARY(K) ) THEN
          K_PARAMETERS = NPARAMS_VARY(K)

!  INITIALIZE CUMULATIVE SOURCE TERM

          NC = 0
          DO V = 1, N_GEOMETRIES
           DO O1 = 1, NSTOKES
            DO Q = 1, K_PARAMETERS
             L_SS_CUMSOURCE(Q,V,O1) = ZERO
            ENDDO
           ENDDO
          ENDDO

!  INITIALISE OPTICAL DEPTH LOOP

          NSTART = NLAYERS
          NUT_PREV = NSTART + 1

!  MAIN LOOP OVER ALL OUTPUT OPTICAL DEPTHS
!  ----------------------------------------

          DO UTA = N_USER_LEVELS, 1, -1

!  LAYER INDEX FOR GIVEN OPTICAL DEPTH

           NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
           NUT = NLEVEL + 1

!  CUMULATIVE SINGLE SCATTER SOURCE TERMS TO LAYER NUT
!  ---------------------------------------------------

!    1. GET LAYER SOURCE TERMS = EXACT SCATTERING * MULTIPLIER
!    2. LOOP OVER LAYERS WORKING UPWARDS TO NUT

           DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            DO V = 1, N_GEOMETRIES
             DO Q = 1, K_PARAMETERS
              DO O1 = 1, NPOLAR
               HELP = ZERO
               DO O2 = 1, NSTOKES
                 LSS = ZMAT_UP(V,N,O1,O2) * LC_UP_MULTIPLIERS(N,Q,V) &
                   + L_ZMAT_UP(Q,V,N,O1,O2) * UP_MULTIPLIERS(N,V)
                 HELP = HELP + FLUXVEC(O2) * LSS
               ENDDO
               L_SS_CUMSOURCE(Q,V,O1) = HELP &
                +  UP_LOSTRANS(N,V)    * L_SS_CUMSOURCE(Q,V,O1) &
                + L_UP_LOSTRANS(N,Q,V) *   SS_CUMSOURCE_UP(V,O1,NC-1)
              ENDDO
             ENDDO
            ENDDO

!  END LAYER LOOP

           ENDDO

!  OFFGRID OUTPUT----------
!  SET FINAL CUMULATIVE SOURCE AND SINGLE SCATTER WEIGHTING FUNCTION

           IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

             UT = PARTLAYERS_OUTINDEX(UTA)
             N  = PARTLAYERS_LAYERIDX(UT)

             DO V = 1, N_GEOMETRIES
              DO Q = 1, K_PARAMETERS
               LMULT = LC_UP_MULTIPLIERS_UT(UT,Q,V)
               DO O1 = 1, NPOLAR
                HELP = ZERO
                DO O2 = 1, NSTOKES
                 LSS = ZMAT_UP(V,N,O1,O2)   * LMULT &
                   + L_ZMAT_UP(Q,V,N,O1,O2) * UP_MULTIPLIERS_UT(UT,V)
                 HELP = HELP + FLUXVEC(O2) * LSS
                ENDDO
                L_FINAL_SOURCE = HELP &
                 +   UP_LOSTRANS_UT(UT,V)   * L_SS_CUMSOURCE(Q,V,O1) &
                 + L_UP_LOSTRANS_UT(UT,Q,V) *   SS_CUMSOURCE_UP(V,O1,NC)
                L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                IF ( DO_COLUMN_LINEARIZATION ) THEN
                  COLUMNWF_SS(Q,UTA,V,O1,UPIDX) = L_SSCORRECTION
                ENDIF
               ENDDO
              ENDDO
             ENDDO

!  ONGRID OUTPUT---------
!  JUST SET TO THE CUMULATIVE SOURCE TERM

           ELSE

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

!  CHECK FOR UPDATING THE RECURSION

           IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
           NUT_PREV = NUT

!  END LOOP OVER OUTPUT OPTICAL DEPTHS

          ENDDO

!  END LOOP OVER VARYING LAYERS

         ENDIF
        ENDDO

!  END UPWELLING JACOBIAN CLAUSE

      ENDIF

!  DEBUG

!     DO K = 1, NLAYERS
!       WRITE(85,*)K,ATMOSWF_SS(1,K,1,1,1,UPIDX),
!    &               ATMOSWF_SS(1,K,3,1,1,UPIDX)
!     ENDDO

!  RECURRENCE RELATION FOR THE DOWNWELLING JACOBIANS
!  =================================================

      IF ( DO_DNWELLING .AND. DO_COLUMN_LINEARIZATION ) THEN

!  START THE MAIN LAYER VARIATION LOOP
!  -----------------------------------

        DO K = 1, NTYPE_VARY
         IF ( DO_RTSOL_VARY(K) ) THEN
          K_PARAMETERS = NPARAMS_VARY(K)

!  INITIALIZE CUMULATIVE SOURCE TERM

          NC = 0
          DO V = 1, N_GEOMETRIES
           DO O1 = 1, NPOLAR
            DO Q = 1, K_PARAMETERS
             L_SS_CUMSOURCE(Q,V,O1) = ZERO
            ENDDO
           ENDDO
          ENDDO

!  INITIALISE OPTICAL DEPTH LOOP

          NSTART = 1
          NUT_PREV = NSTART - 1

!  MAIN LOOP OVER ALL OUTPUT OPTICAL DEPTHS

          DO UTA = 1, N_USER_LEVELS

!  LAYER INDEX FOR GIVEN OPTICAL DEPTH

            NLEVEL = UTAU_LEVEL_MASK_DN(UTA)
            NUT = NLEVEL

!  CUMULATIVE SINGLE SCATTER SOURCE TERMS TO LAYER NUT
!  ---------------------------------------------------

!    1. GET LAYER SOURCE TERMS = EXACT SCATTERING * MULTIPLIER
!    2. LOOP OVER LAYERS WORKING DOWNWARDS TO NUT

            DO N = NSTART, NUT
              NC = N

              DO V = 1, N_GEOMETRIES
               DO Q = 1, K_PARAMETERS
                DO O1 = 1, NPOLAR
                 HELP = ZERO
                 DO O2 = 1, NSTOKES
                   LSS = ZMAT_DN(V,N,O1,O2) * LC_DN_MULTIPLIERS(N,Q,V) &
                     + L_ZMAT_DN(Q,V,N,O1,O2) * DN_MULTIPLIERS(N,V)
                   HELP = HELP + FLUXVEC(O2) * LSS
                 ENDDO
                 L_SS_CUMSOURCE(Q,V,O1) = HELP &
                  +  DN_LOSTRANS(N,V)    * L_SS_CUMSOURCE(Q,V,O1) &
                  + L_DN_LOSTRANS(N,Q,V) *   SS_CUMSOURCE_DN(V,O1,NC-1)
                ENDDO
               ENDDO
              ENDDO

!  END LAYER LOOP

            ENDDO

!  OFFGRID OUTPUT----------
!  SET FINAL CUMULATIVE SOURCE AND SINGLE SCATTER WEIGHTING FUNCTION

           IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

             UT = PARTLAYERS_OUTINDEX(UTA)
             N  = PARTLAYERS_LAYERIDX(UT)

             DO V = 1, N_GEOMETRIES
              DO Q = 1, K_PARAMETERS
               LMULT = LC_DN_MULTIPLIERS_UT(UT,Q,V)
               DO O1 = 1, NPOLAR
                HELP = ZERO
                DO O2 = 1, NSTOKES
                 LSS = ZMAT_DN(V,N,O1,O2)   * LMULT &
                   + L_ZMAT_DN(Q,V,N,O1,O2) * DN_MULTIPLIERS_UT(UT,V)
                 HELP = HELP + FLUXVEC(O2) * LSS
                ENDDO
                L_FINAL_SOURCE = HELP &
                 +   DN_LOSTRANS_UT(UT,V)   * L_SS_CUMSOURCE(Q,V,O1) &
                 + L_DN_LOSTRANS_UT(UT,Q,V) *   SS_CUMSOURCE_DN(V,O1,NC)
                L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                IF ( DO_COLUMN_LINEARIZATION ) THEN
                  COLUMNWF_SS(Q,UTA,V,O1,DNIDX) = L_SSCORRECTION
                ENDIF
               ENDDO
              ENDDO
             ENDDO

!  ONGRID OUTPUT---------
!  JUST SET TO THE CUMULATIVE SOURCE TERM

           ELSE

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

!  CHECK FOR UPDATING THE RECURSION

           IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
           NUT_PREV = NUT

!  END LOOP OVER OUTPUT OPTICAL DEPTHS

          ENDDO

!  END LOOP OVER VARYING LAYERS

         ENDIF
        ENDDO

!  END DOWNWELLING JACOBIAN CLAUSE

      ENDIF

!  FINISH

      RETURN
      END SUBROUTINE VLIDORT_LC_SSCORR_OUTGOING

!

      SUBROUTINE LC_OUTGOING_INTEGRATION_UP ( &
        NLAYERS, NFINELAYERS, &
        DO_PARTIALS, &
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
        DO_COLUMN_LINEARIZATION, N_TOTALCOLUMN_WFS, &
        EXTINCTION, L_EXTINCTION, &
        N_PARTIALS, PARTIALS_IDX, PARTIALS_FINEIDX, &
        SUNPATHS, RADII, NTRAVERSE, ALPHA_ALL, &
        SUNPATHS_P, RADII_P, NTRAVERSE_P,    ALPHA_P, &
        SUNPATHS_FINE, NTRAVERSE_FINE, ALPHA_FINE, &
        MULTIPLIERS,     LOSTRANS,     BOA_ATTN, &
        MULTIPLIERS_P,   LOSTRANS_P, &
        LC_MULTIPLIERS,   L_LOSTRANS,   LC_BOA_ATTN, &
        LC_MULTIPLIERS_P, L_LOSTRANS_P )

!  DOES THE OPTICAL DEPTH INTEGRATION OVER LAYERS.
!  PARTIAL LAYER INTEGRATION ADDED SEPTEMBER 2007.

      USE VLIDORT_PARS

      IMPLICIT NONE

!  INPUTS
!  ------

!  CONTROL

      LOGICAL, INTENT (IN) ::          DO_PARTIALS
      INTEGER, INTENT (IN) ::          NFINELAYERS, NLAYERS
      INTEGER, INTENT (IN) ::          N_PARTIALS
      INTEGER, INTENT (IN) ::          PARTIALS_IDX    (MAX_PARTLAYERS)
      INTEGER, INTENT (IN) ::          PARTIALS_FINEIDX(MAX_PARTLAYERS)

!  PROFILE LINEARIZATION CONTROL INPUTS

      LOGICAL, INTENT (IN) ::          LAYER_VARY_FLAG   (MAXLAYERS)
      INTEGER, INTENT (IN) ::          LAYER_VARY_NUMBER (MAXLAYERS)

!  COLUMN LINEARIZATION CONTROL INPUTS

      LOGICAL, INTENT (IN) ::          DO_COLUMN_LINEARIZATION
      INTEGER, INTENT (IN) ::          N_TOTALCOLUMN_WFS

!  WHOLE LAYERS

      INTEGER, INTENT (IN)          :: NTRAVERSE(0:MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: SUNPATHS(0:MAXLAYERS,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: RADII   (0:MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: ALPHA_ALL  (0:MAXLAYERS)

!  FINE LEVEL

      INTEGER, INTENT (IN)          :: NTRAVERSE_FINE(MAXLAYERS,MAXFINELAYERS)
      DOUBLE PRECISION, INTENT (IN) :: SUNPATHS_FINE &
          (MAXLAYERS,MAXLAYERS,MAXFINELAYERS)
      DOUBLE PRECISION, INTENT (IN) :: ALPHA_FINE (MAXLAYERS,MAXFINELAYERS)

!  PARTIAL LAYERS

      INTEGER, INTENT (IN)          :: NTRAVERSE_P(MAX_PARTLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: SUNPATHS_P (MAX_PARTLAYERS,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: ALPHA_P    (MAX_PARTLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: RADII_P    (MAX_PARTLAYERS)

!  EXTINCTION INPUTS

      DOUBLE PRECISION, INTENT (IN) :: EXTINCTION   (MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: L_EXTINCTION (MAXLAYERS,MAX_ATMOSWFS)

!  OUTPUTS
!  -------

      DOUBLE PRECISION, INTENT (OUT) :: MULTIPLIERS  (MAXLAYERS)
      DOUBLE PRECISION, INTENT (OUT) :: LOSTRANS     (MAXLAYERS)
      DOUBLE PRECISION, INTENT (OUT) :: LC_MULTIPLIERS &
          (MAXLAYERS,MAX_ATMOSWFS)
      DOUBLE PRECISION, INTENT (OUT) :: L_LOSTRANS   (MAXLAYERS,MAX_ATMOSWFS)

      DOUBLE PRECISION, INTENT (OUT) :: MULTIPLIERS_P(MAX_PARTLAYERS)
      DOUBLE PRECISION, INTENT (OUT) :: LOSTRANS_P   (MAX_PARTLAYERS)
      DOUBLE PRECISION, INTENT (OUT) :: LC_MULTIPLIERS_P &
          (MAX_PARTLAYERS,MAX_ATMOSWFS)
      DOUBLE PRECISION, INTENT (OUT) :: L_LOSTRANS_P &
          (MAX_PARTLAYERS,MAX_ATMOSWFS)

      DOUBLE PRECISION, INTENT (OUT) :: BOA_ATTN
      DOUBLE PRECISION, INTENT (OUT) :: LC_BOA_ATTN (MAX_ATMOSWFS)

!  LOCAL ARRAYS
!  ------------

!  LOCAL GEOEMETRY ARRAYS

      DOUBLE PRECISION :: CSQ_FINE ( MAXFINELAYERS)
      DOUBLE PRECISION :: COT_FINE ( MAXFINELAYERS)

!  LOCAL ATTENUATION FACTORS

      DOUBLE PRECISION :: ATTN ( 0:MAXLAYERS )
      DOUBLE PRECISION :: ATTN_FINE ( MAXLAYERS, MAXFINELAYERS )
      DOUBLE PRECISION :: ATTN_P    ( MAX_PARTLAYERS )

!  LOCAL ATTENUATION FACTORS, LINEARIZATIONS

      DOUBLE PRECISION :: L_ATTN ( 0:MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_ATTN_FINE &
          ( MAXLAYERS, MAXFINELAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_ATTN_P ( MAX_PARTLAYERS, MAX_ATMOSWFS  )

!  HELP VARIABLES
!  --------------

      INTEGER ::          N, J, K, Q, UT, NP, NFINE_P

      DOUBLE PRECISION :: TAU, SUM, SALPHA, CALPHA, DFINE1, RAYCON, KN
      DOUBLE PRECISION :: CSQ_1, CSQ_2, COT_1, COT_2, ARGM_1, ARGJ
      DOUBLE PRECISION :: FUNC_1, FUNC_2, TRAN_1, TRAN(MAXFINELAYERS)
      DOUBLE PRECISION :: TERM1, TERM2, DFINE_P, STEP_F, STEP_P
      DOUBLE PRECISION :: L_TERM1, L_TERM2

      DOUBLE PRECISION :: L_FUNC_1, L_FUNC_2, L_TRAN_1, L_TRAN
      DOUBLE PRECISION :: L_SUM, L_FUNC, L_KN, L_IOP, STEP, FUNC

      LOGICAL          :: DO_SPECIAL
      DOUBLE PRECISION :: HN, HJ(MAXFINELAYERS), HP

!  LOCAL OPTICAL THICKNESS CUTOFF
!      (SHOULD BE SAME AS MAX_TAU_SPATH IN LIDORT)

      DOUBLE PRECISION, PARAMETER :: LOCAL_CUTOFF = 32.0D0

!  INITIALISE OUTPUT
!  -----------------

!  WHOLE LAYERS

      DO N = 1, NLAYERS
        MULTIPLIERS(N) = 0.0D0
        LOSTRANS(N)    = 0.0D0
      ENDDO

!  PARTIAL LAYERS

      IF ( DO_PARTIALS ) THEN
        DO UT = 1, N_PARTIALS
          MULTIPLIERS_P(UT) = 0.0D0
          LOSTRANS_P(UT)    = 0.0D0
        ENDDO
      ENDIF

!  WHOLE LAYER LINEARIZATIONS

      IF ( DO_COLUMN_LINEARIZATION ) THEN
        DO N = 1, NLAYERS
          DO Q = 1, N_TOTALCOLUMN_WFS
            LC_MULTIPLIERS(N,Q) = ZERO
            L_LOSTRANS(N,Q)     = ZERO
          ENDDO
        ENDDO
      ENDIF

!  PARTIAL LAYER LINEARIZATIONS

      IF ( DO_PARTIALS ) THEN
       IF ( DO_COLUMN_LINEARIZATION ) THEN
        DO UT = 1, N_PARTIALS
          DO Q = 1, N_TOTALCOLUMN_WFS
            LC_MULTIPLIERS_P(UT,Q) = ZERO
            L_LOSTRANS_P(UT,Q)     = ZERO
          ENDDO
        ENDDO
       ENDIF
      ENDIF

!  CREATE SLANT PATH OPTICAL ATTENUATIONS, SOLAR BEAMS
!  ---------------------------------------------------

!  RAY CONSTANT AT TOA

      IF ( ALPHA_ALL(0).EQ.0.0D0 ) THEN
        DO_SPECIAL = .TRUE.
      ELSE
        DO_SPECIAL = .FALSE.
        SALPHA = DSIN(ALPHA_ALL(0))
        RAYCON  = RADII(0) * SALPHA
      ENDIF
      DFINE1 = DBLE(NFINELAYERS+1)

!  ZERO THE ATTENUATION FACTORS.  EXTREMELY IMPORTANT

      ATTN = 0.0D0
      IF ( NFINELAYERS.GT.0) ATTN_FINE = 0.0D0
      IF ( DO_PARTIALS )     ATTN_P    = 0.0D0
      IF ( DO_COLUMN_LINEARIZATION ) THEN
        L_ATTN = 0.0D0
        IF ( NFINELAYERS.GT.0)L_ATTN_FINE = 0.0D0
        IF ( DO_PARTIALS )    L_ATTN_P    = 0.0D0
      ENDIF

!  ATTENUATION FUNCTIONS, WHOLE LAYERS

      DO N = 0, NLAYERS
       TAU     = 0.0D0
       DO K = 1, NTRAVERSE(N)
         TAU = TAU + SUNPATHS(N,K) * EXTINCTION(K)
       ENDDO
       IF ( TAU .LE. LOCAL_CUTOFF ) ATTN(N) = DEXP(-TAU)
      ENDDO
      BOA_ATTN = ATTN(NLAYERS)

!  ATTENUATIONS FOR PARTIALS

      IF ( DO_PARTIALS ) THEN
        DO UT = 1, N_PARTIALS
          TAU = 0.0D0
          DO K = 1, NTRAVERSE_P(UT)
            TAU = TAU + SUNPATHS_P(UT,K) * EXTINCTION(K)
          ENDDO
          IF ( TAU.LE.LOCAL_CUTOFF) ATTN_P(UT) = DEXP(-TAU)
        ENDDO
      ENDIF

!  ATTENUATIONS FOR FINE LAYER STUFF

      DO N = 1, NLAYERS
       DO J = 1, NFINELAYERS
        TAU            = 0.0D0
        DO K = 1, NTRAVERSE_FINE(N,J)
          TAU = TAU + SUNPATHS_FINE(N,K,J) * EXTINCTION(K)
        ENDDO
        IF ( TAU .LE. LOCAL_CUTOFF ) ATTN_FINE(N,J) = DEXP(-TAU)
       ENDDO
      ENDDO

!  LINEARIZED ATTENUATION FACTORS. COLUMN LINEARIZATION

      IF ( DO_COLUMN_LINEARIZATION ) THEN
        DO N = 0, NLAYERS
          DO Q = 1, N_TOTALCOLUMN_WFS
            SUM = ZERO
            DO K = 1, NTRAVERSE(N)
              L_IOP = - L_EXTINCTION(K,Q)
              SUM = SUM + SUNPATHS(N,K) * L_IOP
            ENDDO
            L_ATTN(N,Q) = SUM * ATTN(N)
            IF (N.EQ.NLAYERS)LC_BOA_ATTN(Q)=L_ATTN(N,Q)
          ENDDO
        ENDDO
        DO N = 1, NLAYERS
         DO J = 1, NFINELAYERS
          DO Q = 1, N_TOTALCOLUMN_WFS
            SUM = ZERO
            DO K = 1, NTRAVERSE_FINE(N,J)
              L_IOP = - L_EXTINCTION(K,Q)
              SUM = SUM + SUNPATHS_FINE(N,K,J) * L_IOP
            ENDDO
            L_ATTN_FINE(N,J,Q) = SUM * ATTN_FINE(N,J)
          ENDDO
         ENDDO
        ENDDO
        IF ( DO_PARTIALS ) THEN
         DO UT = 1, N_PARTIALS
          DO Q = 1, N_TOTALCOLUMN_WFS
            SUM = ZERO
            DO K = 1, NTRAVERSE_P(UT)
              L_IOP = - L_EXTINCTION(K,Q)
              SUM = SUM + SUNPATHS_P(UT,K) * L_IOP
            ENDDO
            L_ATTN_P(UT,Q) = SUM * ATTN_P(UT)
          ENDDO
         ENDDO
        ENDIF
      ENDIF

!   SPECIAL CASE (NADIR VIEWING)
!   ============================

      IF ( DO_SPECIAL ) THEN

!  START LAYER LOOP

        DO N = NLAYERS, 1, -1

!  MULTIPLIERS AND TRANSMITTANCES

          HN = RADII(N-1) - RADII(N)
          KN = EXTINCTION(N)
          TRAN_1 = DEXP ( - HN * KN)
          STEP   = HN/DFINE1
          FUNC_2 = ATTN(N-1)
          FUNC_1 = ATTN(N) * TRAN_1
          SUM = 0.5D0 * ( FUNC_1 + FUNC_2 )
          DO J = 1, NFINELAYERS
            HJ(J) = HN - STEP * DBLE(J)
            TRAN(J) = DEXP ( - KN * HJ(J) )
            FUNC = ATTN_FINE(N,J) * TRAN(J)
            SUM = SUM + FUNC
          ENDDO
          LOSTRANS(N)    = TRAN_1
          MULTIPLIERS(N) = SUM * STEP * KN

!  COLUMN LINEARIZATION.

          IF ( DO_COLUMN_LINEARIZATION ) THEN
           DO Q = 1, N_TOTALCOLUMN_WFS
            L_KN = L_EXTINCTION(N,Q)
            L_FUNC_2 = L_ATTN(N-1,Q)
            L_TRAN_1 = -L_KN * HN * TRAN_1
            L_FUNC_1 = L_ATTN(N,Q) * TRAN_1 + ATTN(N) * L_TRAN_1
            L_SUM = 0.5D0 * ( L_FUNC_1 + L_FUNC_2 )
            DO J = 1, NFINELAYERS
              L_TRAN = - L_KN * HJ(J) * TRAN(J)
              L_FUNC = L_ATTN_FINE(N,J,Q) * TRAN(J) &
                     +   ATTN_FINE(N,J)   * L_TRAN
              L_SUM = L_SUM + L_FUNC
            ENDDO
            L_LOSTRANS(N,Q)      = L_TRAN_1
            LC_MULTIPLIERS(N,Q) = STEP * ( L_KN * SUM + KN * L_SUM )
           ENDDO
          ENDIF

!  END LAYER LOOP

        ENDDO

!  FINISH IF NO PARTIALS

        IF ( .NOT. DO_PARTIALS ) RETURN

!  START PARTIALS LOOP

        DO UT = 1, N_PARTIALS

!  LAYER OF OCCURRENCE AND FINE-LAYER CUTOFF

          NP      = PARTIALS_IDX(UT)
          NFINE_P = PARTIALS_FINEIDX(UT)
          HN = RADII(NP-1) - RADII(NP)
          HP = RADII_P(UT) - RADII(NP)

          IF ( NFINE_P .GT. 0 ) THEN
            DFINE_P = DBLE(NFINE_P)
            STEP_F = HN/DFINE1
            STEP_P = HP - DFINE_P * STEP_F
          ELSE
            STEP_F = 0.0D0
            STEP_P = HP
          ENDIF

!  INTEGRATED SOURCE TERM MULTIPLIER + TRANSMITTANCE
!   ---- TRAPEZIUM RULE INTEGRATION
!  CAREFUL WITH THE LAST STEP OF THE INTEGRATION

          KN = EXTINCTION(NP)
          FUNC_2 = ATTN_P(UT)
          TRAN_1 = DEXP ( - KN * HP )
          FUNC_1 = ATTN(NP) *  TRAN_1
          IF ( NFINE_P.GT.0) THEN
            SUM = 0.5D0 * FUNC_1
            DO J = 1, NFINE_P
              HJ(J) = HP - STEP_F * DBLE(J)
              TRAN(J) = DEXP ( - KN * HJ(J) )
              FUNC = ATTN_FINE(NP,J) * TRAN(J)
              IF ( J.LT.NFINE_P ) SUM = SUM + FUNC
              IF ( J.EQ.NFINE_P ) SUM = SUM + 0.5D0*FUNC
            ENDDO
            TERM1 = SUM * STEP_F
            TERM2 = ( FUNC + FUNC_2) * 0.5D0 * STEP_P
          ELSE
            TERM1 = 0.5D0 * FUNC_1 * STEP_P
            TERM2 = 0.5D0 * FUNC_2 * STEP_P
          ENDIF
          LOSTRANS_P(UT)    = TRAN_1
          MULTIPLIERS_P(UT) = ( TERM1 + TERM2 ) * KN

!        WRITE(*,*)'SPECIAL PARTIAL UP',UT,NP,NFINE_P,
!     &              MULTIPLIERS_P(UT),LOSTRANS_P(UT)

!  COLUMN LINEARIZATION

          IF ( DO_COLUMN_LINEARIZATION ) THEN
            DO Q = 1, N_TOTALCOLUMN_WFS
             L_KN = L_EXTINCTION(NP,Q)
             L_FUNC_2 = L_ATTN_P(UT,Q)
             L_TRAN_1 = -L_KN * HP * TRAN_1
             L_FUNC_1 = L_ATTN(NP,Q) * TRAN_1 + ATTN(NP) * L_TRAN_1
             IF ( NFINE_P.GT.0) THEN
              L_SUM = 0.5D0 * L_FUNC_1
              DO J = 1, NFINE_P
                L_TRAN = - L_KN * HJ(J) * TRAN(J)
                L_FUNC = L_ATTN_FINE(NP,J,Q) * TRAN(J) &
                       +   ATTN_FINE(NP,J)   * L_TRAN
                IF ( J.LT.NFINE_P ) L_SUM = L_SUM + L_FUNC
                IF ( J.EQ.NFINE_P ) L_SUM = L_SUM + 0.5D0*L_FUNC
              ENDDO
              L_TERM1 = L_SUM * STEP_F
              L_TERM2 = ( L_FUNC + L_FUNC_2 ) * 0.5D0 * STEP_P
             ELSE
              L_TERM1 = 0.5D0 * L_FUNC_1 * STEP_P
              L_TERM2 = 0.5D0 * L_FUNC_2 * STEP_P
             ENDIF
             L_LOSTRANS_P(UT,Q)       = L_TRAN_1
             LC_MULTIPLIERS_P(UT,Q) =  ( TERM1 +   TERM2 ) * L_KN + &
                                       ( L_TERM1 + L_TERM2 ) * KN
            ENDDO
          ENDIF

!  END PARTIALS

        ENDDO

!  RETURN AND END SPECIALS CLAUSE

        RETURN

      ENDIF

!  GENERAL CASE: WORK UP FROM THE BOTTOM OF THE ATMOSPHERE
!  =======================================================

!  INITIALISE

      N = NLAYERS
      SALPHA = DSIN(ALPHA_ALL(N))
      CALPHA = DCOS(ALPHA_ALL(N))
      CSQ_1 = 1.0D0 / SALPHA / SALPHA
      COT_1 = CALPHA / SALPHA

!  START LAYER LOOP

      DO N = NLAYERS, 1, -1

!  SAVE SOME QUANTITIES

        KN = RAYCON * EXTINCTION(N)
        SALPHA = DSIN(ALPHA_ALL(N-1))
        CALPHA = DCOS(ALPHA_ALL(N-1))
        CSQ_2 = 1.0D0 / SALPHA / SALPHA
        COT_2 = CALPHA / SALPHA
        STEP    = (ALPHA_ALL(N) - ALPHA_ALL(N-1))/DFINE1

        DO J = 1, NFINELAYERS
          CALPHA = DCOS(ALPHA_FINE(N,J))
          SALPHA = DSIN(ALPHA_FINE(N,J))
          COT_FINE(J) = CALPHA / SALPHA
          CSQ_FINE(J) = 1.0D0 / SALPHA / SALPHA
        ENDDO

!  INTEGRATED SOURCE TERM MULTIPLIER + TRANSMITTANCE
!   ---- TRAPEZIUM RULE INTEGRATION

!   SURELY THIS IS WRONG.................???? 01 OCTOBER 2007
!        FUNC_1 = ATTN(N-1) * CSQ_1
!        ARGM_2 = COT_2 - COT_1
!        TRAN_2 = DEXP ( - KN * ARGM_2 )
!        FUNC_2 = ATTN(N) * CSQ_2 * TRAN_2

        FUNC_2 = ATTN(N-1)   * CSQ_2
        ARGM_1 = COT_2 - COT_1
        TRAN_1 = DEXP ( - KN * ARGM_1 )
        FUNC_1 = ATTN(N) * CSQ_1 * TRAN_1
        SUM = 0.5D0 * ( FUNC_1 + FUNC_2 )
        DO J = 1, NFINELAYERS
          TRAN(J) = DEXP ( -KN * ( COT_2 - COT_FINE(J) ) )
          FUNC = ATTN_FINE(N,J) * TRAN(J) * CSQ_FINE(J)
          SUM = SUM + FUNC
        ENDDO
        LOSTRANS(N)    = TRAN_1
        MULTIPLIERS(N) = SUM * STEP * KN

!  COLUMN LINEARIZATION.

        IF ( DO_COLUMN_LINEARIZATION ) THEN
          DO Q = 1, N_TOTALCOLUMN_WFS
            L_KN = RAYCON * L_EXTINCTION(N,Q)
            L_FUNC_2 = L_ATTN(N-1,Q) * CSQ_2
            L_TRAN_1 = -L_KN * ARGM_1 * TRAN_1
            L_FUNC_1 = CSQ_1 * ( L_ATTN(N,Q) *   TRAN_1 &
                                 + ATTN(N)   * L_TRAN_1 )
            L_SUM = 0.5D0 * ( L_FUNC_1 + L_FUNC_2 )
            DO J = 1, NFINELAYERS
              ARGJ = COT_2 - COT_FINE(J)
              L_TRAN = - L_KN * ARGJ * TRAN(J)
              L_FUNC = CSQ_FINE(J) * ( L_ATTN_FINE(N,J,Q) * TRAN(J) &
                                     +   ATTN_FINE(N,J)   * L_TRAN )
              L_SUM = L_SUM + L_FUNC
            ENDDO
            L_LOSTRANS(N,Q)      = L_TRAN_1
            LC_MULTIPLIERS(N,Q) = STEP * ( L_KN * SUM + KN * L_SUM )
          ENDDO
        ENDIF

!  UPDATE THE GEOMETRY

        COT_1 = COT_2
        CSQ_1 = CSQ_2

!  FINISH LAYER

      ENDDO

!  FINISH IF NO PARTIALS

      IF ( .NOT. DO_PARTIALS ) RETURN

!  PARTIAL LAYER FUNCTIONS
!  =======================

!  START PARTIALS LOOP

      DO UT = 1, N_PARTIALS

!  LAYER OF OCCURRENCE AND FINE-LAYER CUTOFF

        NP     = PARTIALS_IDX(UT)
        NFINE_P = PARTIALS_FINEIDX(UT)
        IF ( NFINE_P .GT. 0 ) THEN
          DFINE_P = DBLE(NFINE_P)
          STEP_F = (ALPHA_ALL(NP) - ALPHA_FINE(NP,NFINE_P))/DFINE_P
          STEP_P = (ALPHA_FINE(NP,NFINE_P)-ALPHA_P(UT))
        ELSE
          STEP_P = (ALPHA_ALL(NP)-ALPHA_P(UT))
        ENDIF

!  BOTTOM OF LAYER OF OCCURRENCE

        SALPHA = DSIN(ALPHA_ALL(NP))
        CALPHA = DCOS(ALPHA_ALL(NP))
        CSQ_1 = 1.0D0 / SALPHA / SALPHA
        COT_1 = CALPHA / SALPHA
        KN = RAYCON * EXTINCTION(NP)

!  TOP OF INTEGRATION = OFF-GRID VALUE

        SALPHA = DSIN(ALPHA_P(UT))
        CALPHA = DCOS(ALPHA_P(UT))
        CSQ_2 = 1.0D0 / SALPHA / SALPHA
        COT_2 = CALPHA / SALPHA

!  SAVED FINE-GRID QUANTITIES AS FAR AS CUTOFF

        DO J = 1, NFINE_P
          CALPHA = DCOS(ALPHA_FINE(NP,J))
          SALPHA = DSIN(ALPHA_FINE(NP,J))
          COT_FINE(J) = CALPHA / SALPHA
          CSQ_FINE(J) = 1.0D0 / SALPHA / SALPHA
        ENDDO

!  INTEGRATED SOURCE TERM MULTIPLIER + TRANSMITTANCE
!   ---- TRAPEZIUM RULE INTEGRATION
!  CAREFUL WITH THE LAST STEP OF THE INTEGRATION

        FUNC_2 = ATTN_P(UT)   * CSQ_2
        ARGM_1 = COT_2 - COT_1
        TRAN_1 = DEXP ( - KN * ARGM_1  )
        FUNC_1 = ATTN(NP) * CSQ_1 * TRAN_1

        IF ( NFINE_P.GT.0) THEN
          SUM = 0.5D0 * FUNC_1
          DO J = 1, NFINE_P
            TRAN(J) = DEXP ( -KN * ( COT_2 - COT_FINE(J) ) )
            FUNC = ATTN_FINE(NP,J) * TRAN(J) * CSQ_FINE(J)
            IF ( J.LT.NFINE_P ) SUM = SUM + FUNC
            IF ( J.EQ.NFINE_P ) SUM = SUM + 0.5D0*FUNC
          ENDDO
          TERM1 = SUM * STEP_F
          TERM2 = (FUNC + FUNC_2) * 0.5D0 * STEP_P
        ELSE
          TERM1 = 0.5D0 * FUNC_1 * STEP_P
          TERM2 = 0.5D0 * FUNC_2 * STEP_P
        ENDIF
        LOSTRANS_P(UT)    = TRAN_1
        MULTIPLIERS_P(UT) = ( TERM1 + TERM2 ) * KN

!  COLUMN LINEARIZATION

        IF ( DO_COLUMN_LINEARIZATION ) THEN
          DO Q = 1, N_TOTALCOLUMN_WFS
            L_KN = RAYCON * L_EXTINCTION(NP,Q)
            L_FUNC_2 = L_ATTN_P(UT,Q) * CSQ_2
            L_TRAN_1 = -L_KN * ARGM_1 * TRAN_1
            L_FUNC_1 = CSQ_1 * ( L_ATTN(NP,Q) *   TRAN_1 &
                                 + ATTN(NP)   * L_TRAN_1 )
            IF ( NFINE_P.GT.0) THEN
              L_SUM = 0.5D0 * L_FUNC_1
              DO J = 1, NFINE_P
                ARGJ = COT_2 - COT_FINE(J)
                L_TRAN = - L_KN * ARGJ * TRAN(J)
                L_FUNC = CSQ_FINE(J) * ( L_ATTN_FINE(NP,J,Q) * TRAN(J) &
                                     +   ATTN_FINE(NP,J)     * L_TRAN )
                IF ( J.LT.NFINE_P ) L_SUM = L_SUM + L_FUNC
                IF ( J.EQ.NFINE_P ) L_SUM = L_SUM + 0.5D0*L_FUNC
              ENDDO
              L_TERM1 = L_SUM * STEP_F
              L_TERM2 = ( L_FUNC + L_FUNC_2 ) * 0.5D0 * STEP_P
            ELSE
              L_TERM1 = 0.5D0 * L_FUNC_1 * STEP_P
              L_TERM2 = 0.5D0 * L_FUNC_2 * STEP_P
            ENDIF
            L_LOSTRANS_P(UT,Q)     = L_TRAN_1
            LC_MULTIPLIERS_P(UT,Q) =  ( TERM1 +   TERM2 ) * L_KN + &
                                      ( L_TERM1 + L_TERM2 ) * KN
          ENDDO
        ENDIF

!  FINISH PARTIALS LOOP

      ENDDO

!  FINISH

      RETURN
      END SUBROUTINE LC_OUTGOING_INTEGRATION_UP

!

      SUBROUTINE LC_OUTGOING_INTEGRATION_DN ( &
        NLAYERS, NFINELAYERS, DO_PARTIALS, &
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
        DO_COLUMN_LINEARIZATION, N_TOTALCOLUMN_WFS, &
        EXTINCTION, L_EXTINCTION, &
        N_PARTIALS, PARTIALS_IDX, PARTIALS_FINEIDX, &
        SUNPATHS, RADII, NTRAVERSE, ALPHA_ALL, &
        SUNPATHS_P,    NTRAVERSE_P,    ALPHA_P, &
        SUNPATHS_FINE, NTRAVERSE_FINE, ALPHA_FINE, &
        MULTIPLIERS,     LOSTRANS, &
        MULTIPLIERS_P,   LOSTRANS_P, &
        LC_MULTIPLIERS,   L_LOSTRANS, &
        LC_MULTIPLIERS_P, L_LOSTRANS_P )

!  DOES THE OPTICAL DEPTH INTEGRATION OVER LAYERS.
!  PARTIAL LAYER INTEGRATION ADDED SEPTEMBER 2007.

      USE VLIDORT_PARS

      IMPLICIT NONE

!  INPUTS
!  ------

!  CONTROL

      LOGICAL, INTENT (IN) ::          DO_PARTIALS
      INTEGER, INTENT (IN) ::          NFINELAYERS, NLAYERS
      INTEGER, INTENT (IN) ::          N_PARTIALS
      INTEGER, INTENT (IN) ::          PARTIALS_IDX    (MAX_PARTLAYERS)
      INTEGER, INTENT (IN) ::          PARTIALS_FINEIDX(MAX_PARTLAYERS)

!  PROFILE LINEARIZATION CONTROL INPUTS

      LOGICAL, INTENT (IN) ::          LAYER_VARY_FLAG   (MAXLAYERS)
      INTEGER, INTENT (IN) ::          LAYER_VARY_NUMBER (MAXLAYERS)

!  COLUMN LINEARIZATION CONTROL INPUTS

      LOGICAL, INTENT (IN) ::          DO_COLUMN_LINEARIZATION
      INTEGER, INTENT (IN) ::          N_TOTALCOLUMN_WFS

!  WHOLE LAYERS

      INTEGER, INTENT (IN) ::          NTRAVERSE(0:MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: SUNPATHS(0:MAXLAYERS,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: RADII   (0:MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: ALPHA_ALL  (0:MAXLAYERS)

!  FINE LEVEL

      INTEGER, INTENT (IN) ::          NTRAVERSE_FINE(MAXLAYERS,MAXFINELAYERS)
      DOUBLE PRECISION, INTENT (IN) :: SUNPATHS_FINE &
          (MAXLAYERS,MAXLAYERS,MAXFINELAYERS)
      DOUBLE PRECISION, INTENT (IN) :: ALPHA_FINE    (MAXLAYERS,MAXFINELAYERS)

!  PARTIAL LAYERS

      INTEGER, INTENT (IN) ::          NTRAVERSE_P(MAX_PARTLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: SUNPATHS_P (MAX_PARTLAYERS,MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: ALPHA_P    (MAX_PARTLAYERS)

!  EXTINCTION INPUTS

      DOUBLE PRECISION, INTENT (IN) :: EXTINCTION   (MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: L_EXTINCTION (MAXLAYERS,MAX_ATMOSWFS)

!  OUTPUTS
!  -------

      DOUBLE PRECISION, INTENT (OUT) :: MULTIPLIERS  (MAXLAYERS)
      DOUBLE PRECISION, INTENT (OUT) :: LOSTRANS     (MAXLAYERS)
      DOUBLE PRECISION, INTENT (OUT) :: LC_MULTIPLIERS &
          (MAXLAYERS,MAX_ATMOSWFS)
      DOUBLE PRECISION, INTENT (OUT) :: L_LOSTRANS   (MAXLAYERS,MAX_ATMOSWFS)

      DOUBLE PRECISION, INTENT (OUT) :: MULTIPLIERS_P   (MAX_PARTLAYERS)
      DOUBLE PRECISION, INTENT (OUT) :: LOSTRANS_P      (MAX_PARTLAYERS)
      DOUBLE PRECISION, INTENT (OUT) :: LC_MULTIPLIERS_P &
          (MAX_PARTLAYERS,MAX_ATMOSWFS)
      DOUBLE PRECISION, INTENT (OUT) :: L_LOSTRANS_P &
          (MAX_PARTLAYERS,MAX_ATMOSWFS)

!  LOCAL ARRAYS
!  ------------

!  LOCAL GEOEMETRY ARRAYS

      DOUBLE PRECISION :: CSQ_FINE ( MAXFINELAYERS)
      DOUBLE PRECISION :: COT_FINE ( MAXFINELAYERS)

!  LOCAL ATTENUATION FACTORS

      DOUBLE PRECISION :: ATTN ( 0:MAXLAYERS )
      DOUBLE PRECISION :: ATTN_FINE ( MAXLAYERS, MAXFINELAYERS )
      DOUBLE PRECISION :: ATTN_P    ( MAX_PARTLAYERS )

!  LOCAL ATTENUATION FACTORS, LINEARIZATIONS

      DOUBLE PRECISION :: L_ATTN ( 0:MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_ATTN_FINE &
          ( MAXLAYERS, MAXFINELAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_ATTN_P ( MAX_PARTLAYERS, MAX_ATMOSWFS  )

!  HELP VARIABLES
!  --------------

      INTEGER ::          N, J, K, Q, UT, NP, NFINE_P

      DOUBLE PRECISION :: TAU, SUM, SALPHA, CALPHA, DFINE1, RAYCON, KN
      DOUBLE PRECISION :: CSQ_1, CSQ_2, COT_1, COT_2, ARGM_1, ARGJ
      DOUBLE PRECISION :: FUNC_1, FUNC_2, TRAN_1, TRAN(MAXFINELAYERS)
      DOUBLE PRECISION :: TERM1, TERM2, DFINE_P, STEP_F, STEP_P
      DOUBLE PRECISION :: L_TERM1, L_TERM2

      DOUBLE PRECISION :: L_FUNC_1, L_FUNC_2, L_TRAN_1, L_TRAN
      DOUBLE PRECISION :: L_SUM, L_FUNC, L_KN, L_IOP, STEP, FUNC

!  LOCAL OPTICAL THICKNESS CUTOFF
!      (SHOULD BE SAME AS MAX_TAU_SPATH IN LIDORT)

      DOUBLE PRECISION, PARAMETER :: LOCAL_CUTOFF = 32.0D0

!  INITIALISE OUTPUT
!  -----------------

!  WHOLE LAYERS

      DO N = 1, NLAYERS
        MULTIPLIERS(N) = 0.0D0
        LOSTRANS(N)    = 0.0D0
      ENDDO

!  PARTIAL LAYERS

      IF ( DO_PARTIALS ) THEN
        DO UT = 1, N_PARTIALS
          MULTIPLIERS_P(UT) = 0.0D0
          LOSTRANS_P(UT)    = 0.0D0
        ENDDO
      ENDIF

!  LINEARIZATIONS

      IF ( DO_COLUMN_LINEARIZATION ) THEN
        DO N = 1, NLAYERS
          DO Q = 1, N_TOTALCOLUMN_WFS
            LC_MULTIPLIERS(N,Q) = ZERO
            L_LOSTRANS(N,Q)     = ZERO
          ENDDO
        ENDDO
      ENDIF

!  PARTIAL LAYER LINEARIZATIONS

      IF ( DO_PARTIALS ) THEN
       IF ( DO_COLUMN_LINEARIZATION ) THEN
        DO UT = 1, N_PARTIALS
          DO Q = 1, N_TOTALCOLUMN_WFS
            LC_MULTIPLIERS_P(UT,Q) = ZERO
            L_LOSTRANS_P(UT,Q)     = ZERO
          ENDDO
        ENDDO
       ENDIF
      ENDIF

!  CREATE SLANT PATH OPTICAL ATTENUATIONS, SOLAR BEAMS
!  ---------------------------------------------------

!  RAY CONSTANT AT TOA

      SALPHA = DSIN(ALPHA_ALL(0))
      RAYCON  = RADII(0) * SALPHA
      DFINE1 = DBLE(NFINELAYERS+1)

!  ZERO THE ATTENUATION FACTORS.  EXTREMELY IMPORTANT

      ATTN = 0.0D0
      IF ( NFINELAYERS.GT.0) ATTN_FINE = 0.0D0
      IF ( DO_PARTIALS )     ATTN_P    = 0.0D0
      IF ( DO_COLUMN_LINEARIZATION ) THEN
        L_ATTN = 0.0D0
        IF ( NFINELAYERS.GT.0)L_ATTN_FINE = 0.0D0
        IF ( DO_PARTIALS )    L_ATTN_P    = 0.0D0
      ENDIF

!  ATTENUATION FUNCTIONS, WHOLE LAYERS

      DO N = 0, NLAYERS
       TAU     = 0.0D0
       DO K = 1, NTRAVERSE(N)
         TAU = TAU + SUNPATHS(N,K) * EXTINCTION(K)
       ENDDO
       IF ( TAU .LE. LOCAL_CUTOFF ) ATTN(N) = DEXP(-TAU)
      ENDDO

!  ATTENUATIONS FOR PARTIALS

      IF ( DO_PARTIALS ) THEN
        DO UT = 1, N_PARTIALS
          TAU = 0.0D0
          DO K = 1, NTRAVERSE_P(UT)
            TAU = TAU + SUNPATHS_P(UT,K) * EXTINCTION(K)
          ENDDO
          IF ( TAU.LE.LOCAL_CUTOFF) ATTN_P(UT) = DEXP(-TAU)
        ENDDO
      ENDIF

!  ATTENUATIONS FOR FINE LAYER STUFF

      DO N = 1, NLAYERS
       DO J = 1, NFINELAYERS
        TAU            = 0.0D0
        DO K = 1, NTRAVERSE_FINE(N,J)
          TAU = TAU + SUNPATHS_FINE(N,K,J) * EXTINCTION(K)
        ENDDO
        IF ( TAU .LE. LOCAL_CUTOFF ) ATTN_FINE(N,J) = DEXP(-TAU)
       ENDDO
      ENDDO

!  LINEARIZED ATTENUATION FACTORS. COLUMN LINEARIZATION

      IF ( DO_COLUMN_LINEARIZATION ) THEN
        DO N = 0, NLAYERS
          DO Q = 1, N_TOTALCOLUMN_WFS
            SUM = ZERO
            DO K = 1, NTRAVERSE(N)
              L_IOP = - L_EXTINCTION(K,Q)
              SUM = SUM + SUNPATHS(N,K) * L_IOP
            ENDDO
            L_ATTN(N,Q) = SUM * ATTN(N)
          ENDDO
        ENDDO
        DO N = 1, NLAYERS
         DO J = 1, NFINELAYERS
          DO Q = 1, N_TOTALCOLUMN_WFS
            SUM = ZERO
            DO K = 1, NTRAVERSE_FINE(N,J)
              L_IOP = - L_EXTINCTION(K,Q)
              SUM = SUM + SUNPATHS_FINE(N,K,J) * L_IOP
            ENDDO
            L_ATTN_FINE(N,J,Q) = SUM * ATTN_FINE(N,J)
          ENDDO
         ENDDO
        ENDDO
        IF ( DO_PARTIALS ) THEN
         DO UT = 1, N_PARTIALS
          DO Q = 1, N_TOTALCOLUMN_WFS
            SUM = ZERO
            DO K = 1, NTRAVERSE_P(UT)
              L_IOP = - L_EXTINCTION(K,Q)
              SUM = SUM + SUNPATHS_P(UT,K) * L_IOP
            ENDDO
            L_ATTN_P(UT,Q) = SUM * ATTN_P(UT)
          ENDDO
         ENDDO
        ENDIF
      ENDIF

!  WORK DOWN FROM THE TOP OF THE ATMOSPHERE
!  ========================================

!  INITIALISE

      N = 0
      SALPHA = DSIN(ALPHA_ALL(N))
      CALPHA = DCOS(ALPHA_ALL(N))
      CSQ_1 = 1.0D0 / SALPHA / SALPHA
      COT_1 = CALPHA / SALPHA

!  START LAYER LOOP

      DO N = 1, NLAYERS

!  SAVE SOME QUANTITIES

        KN = RAYCON * EXTINCTION(N)
        SALPHA = DSIN(ALPHA_ALL(N))
        CALPHA = DCOS(ALPHA_ALL(N))
        CSQ_2 = 1.0D0 / SALPHA / SALPHA
        COT_2 = CALPHA / SALPHA
        STEP    = (ALPHA_ALL(N) - ALPHA_ALL(N-1))/DFINE1
        DO J = 1, NFINELAYERS
          CALPHA = DCOS(ALPHA_FINE(N,J))
          SALPHA = DSIN(ALPHA_FINE(N,J))
          COT_FINE(J) = CALPHA / SALPHA
          CSQ_FINE(J) = 1.0D0 / SALPHA / SALPHA
        ENDDO

!  INTEGRATED SOURCE TERM MULTIPLIER + TRANSMITTANCE
!   ---- TRAPEZIUM RULE INTEGRATION

        ARGM_1 = COT_1 - COT_2
        TRAN_1 = DEXP ( - KN * ARGM_1 )
        FUNC_1 = ATTN(N-1) * CSQ_1 * TRAN_1
        FUNC_2 = ATTN(N) * CSQ_2
        SUM = 0.5D0 * ( FUNC_1 + FUNC_2 )
        DO J = NFINELAYERS, 1, -1
          TRAN(J) = DEXP ( -KN * ( COT_FINE(J) - COT_2 ) )
          FUNC = ATTN_FINE(N,J) * TRAN(J) * CSQ_FINE(J)
          SUM = SUM + FUNC
        ENDDO
        LOSTRANS(N)    = TRAN_1
        MULTIPLIERS(N) = SUM * STEP * KN

!  COLUMN LINEARIZATION

        IF ( DO_COLUMN_LINEARIZATION ) THEN
          DO Q = 1, N_TOTALCOLUMN_WFS
            L_KN = RAYCON * L_EXTINCTION(N,Q)
            L_FUNC_2 = L_ATTN(N,Q) * CSQ_2
            L_TRAN_1 = -L_KN * ARGM_1 * TRAN_1
            L_FUNC_1 = CSQ_1 * ( L_ATTN(N-1,Q) *   TRAN_1 &
                                 + ATTN(N-1)   * L_TRAN_1 )
            L_SUM = 0.5D0 * ( L_FUNC_1 + L_FUNC_2 )
            DO J = NFINELAYERS, 1, -1
              ARGJ = COT_FINE(J) - COT_2
              L_TRAN = - L_KN * ARGJ * TRAN(J)
              L_FUNC = CSQ_FINE(J) * ( L_ATTN_FINE(N,J,Q) * TRAN(J) &
                                     +   ATTN_FINE(N,J)   * L_TRAN )
              L_SUM = L_SUM + L_FUNC
            ENDDO
            L_LOSTRANS(N,Q)       = L_TRAN_1
            LC_MULTIPLIERS(N,Q) = STEP * ( L_KN * SUM + KN * L_SUM )
          ENDDO
        ENDIF

!  UPDATE THE GEOMETRY

        COT_1 = COT_2
        CSQ_1 = CSQ_2

!  FINISH LAYER

      ENDDO

!  FINISH IF NO PARTIALS

      IF ( .NOT. DO_PARTIALS ) RETURN

!  PARTIAL LAYER FUNCTIONS
!  =======================

!  START PARTIALS LOOP

      DO UT = 1, N_PARTIALS

!  LAYER OF OCCURRENCE AND FINE-LAYER CUTOFF

        NP     = PARTIALS_IDX(UT)
        NFINE_P = PARTIALS_FINEIDX(UT)
        IF ( NFINE_P .EQ. NFINELAYERS ) THEN
          STEP_P = (ALPHA_P(UT)-ALPHA_ALL(NP-1))
        ELSE IF ( NFINE_P.EQ.0 ) THEN
          STEP_F = (ALPHA_ALL(NP) - ALPHA_ALL(NP-1))/DFINE1
          STEP_P = (ALPHA_P(UT)-ALPHA_FINE(NP,NFINE_P+1))
        ELSE
          DFINE_P = DBLE(NFINE_P)
          STEP_F = (ALPHA_ALL(NP) - ALPHA_FINE(NP,NFINE_P))/DFINE_P
          STEP_P = (ALPHA_P(UT)-ALPHA_FINE(NP,NFINE_P+1))
        ENDIF

!  TOP OF LAYER OF OCCURRENCE

        SALPHA = DSIN(ALPHA_ALL(NP-1))
        CALPHA = DCOS(ALPHA_ALL(NP-1))
        CSQ_1 = 1.0D0 / SALPHA / SALPHA
        COT_1 = CALPHA / SALPHA
        KN = RAYCON * EXTINCTION(NP)

!  END OF INTEGRATION = OFF-GRID VALUE

        SALPHA = DSIN(ALPHA_P(UT))
        CALPHA = DCOS(ALPHA_P(UT))
        CSQ_2 = 1.0D0 / SALPHA / SALPHA
        COT_2 = CALPHA / SALPHA

!  SAVED FINE-GRID QUANTITIES AS FAR AS CUTOFF

        DO J = NFINELAYERS, NFINE_P+1, -1
          CALPHA = DCOS(ALPHA_FINE(NP,J))
          SALPHA = DSIN(ALPHA_FINE(NP,J))
          COT_FINE(J) = CALPHA / SALPHA
          CSQ_FINE(J) = 1.0D0 / SALPHA / SALPHA
        ENDDO

!  INTEGRATED SOURCE TERM MULTIPLIER + TRANSMITTANCE
!   ---- TRAPEZIUM RULE INTEGRATION
!  CAREFUL WITH THE LAST STEP OF THE INTEGRATION

        FUNC_2 = ATTN_P(UT)   * CSQ_2
        ARGM_1 = COT_1 - COT_2
        TRAN_1 = DEXP ( - KN * ARGM_1 )
        FUNC_1 = ATTN(NP-1) * CSQ_1 * TRAN_1
        IF ( NFINE_P.LT.NFINELAYERS) THEN
          SUM = 0.5D0 * FUNC_1
          DO J = NFINELAYERS, NFINE_P+1, -1
            TRAN(J) = DEXP ( -KN * ( COT_FINE(J) - COT_2 ) )
            FUNC  = ATTN_FINE(NP,J) * TRAN(J) * CSQ_FINE(J)
            IF ( J.GT.NFINE_P+1 ) SUM = SUM + FUNC
            IF ( J.EQ.NFINE_P+1 ) SUM = SUM + 0.5D0*FUNC
          ENDDO
          TERM1 = SUM * STEP_F
          TERM2 = ( FUNC + FUNC_2) * 0.5D0 * STEP_P
        ELSE
          TERM1 = 0.5D0 * FUNC_1 * STEP_P
          TERM2 = 0.5D0 * FUNC_2 * STEP_P
        ENDIF
        LOSTRANS_P(UT)    = TRAN_1
        MULTIPLIERS_P(UT) = ( TERM1 + TERM2) * KN

!  COLUMN LINEARIZATION.

        IF ( DO_COLUMN_LINEARIZATION ) THEN
          DO Q = 1, N_TOTALCOLUMN_WFS
            L_KN = RAYCON * L_EXTINCTION(NP,Q)
            L_FUNC_2 = L_ATTN_P(UT,Q) * CSQ_2
            L_TRAN_1 = -L_KN * ARGM_1 * TRAN_1
            L_FUNC_1 = CSQ_1 * ( L_ATTN(NP-1,Q) *   TRAN_1 &
                                 + ATTN(NP-1)   * L_TRAN_1 )
            IF ( NFINE_P.LT.NFINELAYERS) THEN
              L_SUM = 0.5D0 * L_FUNC_1
              DO J = NFINELAYERS, NFINE_P+1, -1
               ARGJ = COT_FINE(J) - COT_2
               L_TRAN = - L_KN * ARGJ * TRAN(J)
               L_FUNC = CSQ_FINE(J) * ( L_ATTN_FINE(NP,J,Q) * TRAN(J) &
                                     +    ATTN_FINE(NP,J)   * L_TRAN )
               IF ( J.GT.NFINE_P+1 ) L_SUM = L_SUM + L_FUNC
               IF ( J.EQ.NFINE_P+1 ) L_SUM = L_SUM + 0.5D0*L_FUNC
              ENDDO
              L_TERM1 = L_SUM * STEP_F
              L_TERM2 = ( L_FUNC + L_FUNC_2 ) * 0.5D0 * STEP_P
            ELSE
              L_TERM1 = 0.5D0 * L_FUNC_1 * STEP_P
              L_TERM2 = 0.5D0 * L_FUNC_2 * STEP_P
            ENDIF
            L_LOSTRANS_P(UT,Q)       = L_TRAN_1
            LC_MULTIPLIERS_P(UT,Q) =  ( TERM1 +   TERM2 ) * L_KN + &
                                      ( L_TERM1 + L_TERM2 ) * KN
          ENDDO
        ENDIF

!  FINISH PARTIALS LOOP

      ENDDO

!  FINISH

      RETURN
      END SUBROUTINE LC_OUTGOING_INTEGRATION_DN

!

      SUBROUTINE VLIDORT_LAC_DBCORRECTION ( &
        FLUXMULT, &
        DO_SSCORR_OUTGOING, DO_REFRACTIVE_GEOMETRY, &
        DO_UPWELLING, DO_OBSERVATION_GEOMETRY, &
        NSTOKES, NLAYERS, &
        SZA_LOCAL_INPUT, N_USER_RELAZMS, &
        N_USER_LEVELS, &
        DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
        FLUXVEC, COS_SZANGLES, NBEAMS, &
        N_USER_STREAMS, SZANGLES_ADJUST, &
        PARTLAYERS_OUTFLAG, &
        PARTLAYERS_OUTINDEX, UTAU_LEVEL_MASK_UP, &
        PARTLAYERS_LAYERIDX, N_GEOMETRIES, &
        VZA_OFFSETS, DELTAU_SLANT, &
        TRANS_SOLAR_BEAM, DO_REFLECTED_DIRECTBEAM, &
        T_DELT_USERM, T_UTUP_USERM, &
        N_TOTALCOLUMN_WFS, &
        L_DELTAU_VERT, L_T_DELT_USERM, &
        L_T_UTUP_USERM, &
        DB_CUMSOURCE, EXACTDB_SOURCE, &
        UP_LOSTRANS, UP_LOSTRANS_UT, &
        BOA_ATTN, ATTN_DB_SAVE, &
        L_UP_LOSTRANS, L_UP_LOSTRANS_UT, &
        LC_BOA_ATTN, L_EXACTDBC_SOURCE, &
        L_DB_CUMSOURCE, COLUMNWF_DB )

!  PREPARES LINEARIZATION OF EXACT DIRECT BEAM REFLECTION
!   LINEARIZATION WITH RESPECT TO COLUMN ATMOSPHERIC VARIABLES

      USE VLIDORT_PARS

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT (IN) :: FLUXMULT
      LOGICAL, INTENT (IN) ::          DO_SSCORR_OUTGOING
      LOGICAL, INTENT (IN) ::          DO_REFRACTIVE_GEOMETRY
      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NLAYERS
      DOUBLE PRECISION, INTENT (IN) :: SZA_LOCAL_INPUT &
          ( 0:MAXLAYERS, MAX_SZANGLES )
      INTEGER, INTENT (IN) ::          N_USER_RELAZMS
      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      DOUBLE PRECISION, INTENT (IN) :: LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: FLUXVEC ( MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: COS_SZANGLES ( MAX_SZANGLES )
      INTEGER, INTENT (IN) ::          NBEAMS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      DOUBLE PRECISION, INTENT (IN) :: SZANGLES_ADJUST &
           ( MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS )
      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      INTEGER, INTENT (IN) ::          N_GEOMETRIES
      INTEGER, INTENT (IN) ::          VZA_OFFSETS &
          ( MAX_SZANGLES, MAX_USER_VZANGLES )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_SLANT &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: TRANS_SOLAR_BEAM ( MAXBEAMS )
      LOGICAL, INTENT (IN) ::          DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      INTEGER, INTENT (IN) ::          N_TOTALCOLUMN_WFS
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTUP_USERM &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: DB_CUMSOURCE &
          ( MAX_GEOMETRIES, MAXSTOKES, 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: EXACTDB_SOURCE &
          ( MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: UP_LOSTRANS ( MAXLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (IN) :: UP_LOSTRANS_UT &
          ( MAX_PARTLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (IN) :: BOA_ATTN ( MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (IN) :: ATTN_DB_SAVE ( MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (IN) :: L_UP_LOSTRANS &
         ( MAXLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (IN) :: L_UP_LOSTRANS_UT &
          ( MAX_PARTLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (IN) :: LC_BOA_ATTN &
         ( MAX_ATMOSWFS, MAX_GEOMETRIES )

      DOUBLE PRECISION, INTENT (OUT) :: L_EXACTDBC_SOURCE &
          ( MAX_ATMOSWFS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: L_DB_CUMSOURCE &
          ( MAX_GEOMETRIES, MAXSTOKES )

      DOUBLE PRECISION, INTENT (INOUT) :: COLUMNWF_DB &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  LOCAL VARIABLES
!  ---------------

      INTEGER ::          N, NUT, NSTART, NUT_PREV, NLEVEL, NELEMENTS
      INTEGER ::          UT, UTA, UM, UA, NC, IB, V, K, Q, O1, LUM, LUA
      DOUBLE PRECISION :: X0_FLUX, X0_BOA, FACTOR, L_ATTN, TR, LTR, RT

!  FIRST STAGE
!  -----------

!  RETURN IF NO UPWELLING

      IF ( .NOT.DO_UPWELLING ) RETURN

!  DEFINE LOCAL INDICES

      LUM = 1
      LUA = 1

!  START WEIGHTING FUNCTION LOOP

      DO Q = 1, N_TOTALCOLUMN_WFS

!  INITIALIZE THE OUTPUT RESULTS

        DO V = 1, N_GEOMETRIES
          DO O1 = 1, NSTOKES
            L_EXACTDBC_SOURCE(Q,V,O1) = ZERO
            DO UTA = 1, N_USER_LEVELS
              COLUMNWF_DB(Q,UTA,V,O1) = ZERO
            ENDDO
          ENDDO
        ENDDO

!  NEW CODE  R. SPURR, 6 AUGUST 2007. RT SOLUTIONS INC.
!  ====================================================

!  INITIALIZE CUMULATIVE SOURCE TERM

!  Lattice calculation

        IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN

!  START GEOMETRY LOOPS

          DO UM = 1, N_USER_STREAMS
            DO IB = 1, NBEAMS
              IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
                DO UA = 1, N_USER_RELAZMS
                  V = VZA_OFFSETS(IB,UM) + UA

!  BEAM ATTENUATION AND REFLECTION
!  INITIALIZE CUMULATIVE SOURCE TERM
!    (ALREADY FLUX-MULTIPLIED, BECAUSE EXACTBD IS A FINAL RESULT)

                  IF ( DO_SSCORR_OUTGOING ) THEN
                   L_ATTN = ZERO
                   IF (BOA_ATTN(V).GT.ZERO ) THEN
                    X0_BOA = DCOS(SZANGLES_ADJUST(UM,IB,UA)*DEG_TO_RAD)
                    L_ATTN = LC_BOA_ATTN(Q,V)
                   ENDIF
                  ELSE
                   IF ( DO_REFRACTIVE_GEOMETRY ) THEN
                    X0_BOA = DCOS(SZA_LOCAL_INPUT(NLAYERS,IB)*DEG_TO_RAD)
                   ELSE
                    X0_BOA = COS_SZANGLES(IB)
                   ENDIF
                   L_ATTN = ZERO
                   DO K = 1, NLAYERS
                     L_ATTN = L_ATTN - L_DELTAU_VERT(Q,K) * &
                                          DELTAU_SLANT(NLAYERS,K,IB)
                   ENDDO
                   L_ATTN = L_ATTN * TRANS_SOLAR_BEAM(IB)
                  ENDIF
                  X0_FLUX = FOUR * X0_BOA
                  L_ATTN  = L_ATTN * X0_FLUX

!  LOOP OVER ALBEDO KERNELS, ASSIGN REFLECTION

                  IF ( DO_LAMBERTIAN_SURFACE ) THEN
                    O1 = 1
                    FACTOR = LAMBERTIAN_ALBEDO * L_ATTN
!mick fix 8/16/2013
                    !L_EXACTDBC_SOURCE(Q,V,O1) = FACTOR * FLUXMULT
                    L_EXACTDBC_SOURCE(Q,V,O1) = FACTOR * FLUXMULT * FLUXVEC(O1)
                  ELSE
                    IF ( ATTN_DB_SAVE(V).NE.ZERO ) THEN
                      RT = L_ATTN / ATTN_DB_SAVE(V)
                    ELSE
                      RT = ZERO
                    ENDIF
                    DO O1 = 1, NSTOKES
                      FACTOR = EXACTDB_SOURCE(V,O1) * RT
!mick fix 8/16/2013
                      !L_EXACTDBC_SOURCE(Q,V,O1) = FACTOR * FLUXMULT
                      L_EXACTDBC_SOURCE(Q,V,O1) = FACTOR * FLUXMULT * FLUXVEC(O1)
                    ENDDO
                  ENDIF

!  FINISH LOOPS OVER GEOMETRIES

                ENDDO
              ENDIF
            ENDDO
          ENDDO

!  End lattice calculation

        ENDIF

!  Observational geometry calculation

        IF ( DO_OBSERVATION_GEOMETRY ) THEN

!  START GEOMETRY LOOP

          DO V = 1, N_GEOMETRIES
            IF ( DO_REFLECTED_DIRECTBEAM(V) ) THEN

!  BEAM ATTENUATION AND REFLECTION
!  INITIALIZE CUMULATIVE SOURCE TERM
!    (ALREADY FLUX-MULTIPLIED, BECAUSE EXACTBD IS A FINAL RESULT)

              IF ( DO_SSCORR_OUTGOING ) THEN
                L_ATTN = ZERO
                IF (BOA_ATTN(V).GT.ZERO ) THEN
                  X0_BOA = DCOS(SZANGLES_ADJUST(LUM,V,LUA)*DEG_TO_RAD)
                  L_ATTN = LC_BOA_ATTN(Q,V)
                ENDIF
              ELSE
                IF ( DO_REFRACTIVE_GEOMETRY ) THEN
                  X0_BOA = DCOS(SZA_LOCAL_INPUT(NLAYERS,V)*DEG_TO_RAD)
                ELSE
                  X0_BOA = COS_SZANGLES(V)
                ENDIF
                L_ATTN = ZERO
                DO K = 1, NLAYERS
                  L_ATTN = L_ATTN - L_DELTAU_VERT(Q,K) * &
                                      DELTAU_SLANT(NLAYERS,K,V)
                ENDDO
                L_ATTN = L_ATTN * TRANS_SOLAR_BEAM(V)
              ENDIF
              X0_FLUX = FOUR * X0_BOA
              L_ATTN  = L_ATTN * X0_FLUX

!  LOOP OVER ALBEDO KERNELS, ASSIGN REFLECTION

              IF ( DO_LAMBERTIAN_SURFACE ) THEN
                O1 = 1
                FACTOR = LAMBERTIAN_ALBEDO * L_ATTN
!mick fix 8/16/2013
                !L_EXACTDBC_SOURCE(Q,V,O1) = FACTOR * FLUXMULT
                L_EXACTDBC_SOURCE(Q,V,O1) = FACTOR * FLUXMULT * FLUXVEC(O1)
              ELSE
                IF ( ATTN_DB_SAVE(V).NE.ZERO ) THEN
                  RT = L_ATTN / ATTN_DB_SAVE(V)
                ELSE
                  RT = ZERO
                ENDIF
                DO O1 = 1, NSTOKES
                  FACTOR = EXACTDB_SOURCE(V,O1) * RT
!mick fix 8/16/2013
                  !L_EXACTDBC_SOURCE(Q,V,O1) = FACTOR * FLUXMULT
                  L_EXACTDBC_SOURCE(Q,V,O1) = FACTOR * FLUXMULT * FLUXVEC(O1)
                ENDDO
              ENDIF

!  FINISH LOOP OVER GEOMETRIES

            ENDIF
          ENDDO

!  End observationl geometry calculation

        ENDIF

!  UPWELLING RECURRENCE: TRANSMITTANCE OF EXACT SOURCE TERM
!  --------------------------------------------------------

!  NUMBER OF ELEMENTS

        NELEMENTS = NSTOKES
        IF ( DO_LAMBERTIAN_SURFACE ) NELEMENTS = 1

!  INITIALIZE CUMULATIVE SOURCE TERM
!    (ALREADY FLUX-MULTIPLIED, BECAUSE EXACTBD IS A FINAL RESULT)

        NC =  0
        DO V = 1, N_GEOMETRIES
          DO O1 = 1, NELEMENTS
            L_DB_CUMSOURCE(V,O1) = L_EXACTDBC_SOURCE(Q,V,O1)
          ENDDO
        ENDDO

!  INITIALIZE OPTICAL DEPTH LOOP

        NSTART = NLAYERS
        NUT_PREV = NSTART + 1

!  MAIN LOOP OVER ALL OUTPUT OPTICAL DEPTHS

        DO UTA = N_USER_LEVELS, 1, -1

!  LAYER INDEX FOR GIVEN OPTICAL DEPTH

          NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
          NUT    = NLEVEL + 1

!  CUMULATIVE LAYER TRANSMITTANCE :
!    LOOP OVER LAYERS WORKING UPWARDS TO LEVEL NUT
!   ADDITION LINEARIZATION OF TRANSMITTANCE

          DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N

            IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
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
                    DO O1 = 1, NELEMENTS
                     L_DB_CUMSOURCE(V,O1) = TR * L_DB_CUMSOURCE(V,O1) &
                            + LTR * DB_CUMSOURCE(V,O1,NC-1)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ELSE
              DO V = 1, N_GEOMETRIES
                IF ( DO_SSCORR_OUTGOING ) THEN
                  TR  = UP_LOSTRANS(N,V)
                  LTR = L_UP_LOSTRANS(N,Q,V)
                ELSE
                  TR  = T_DELT_USERM(N,V)
                  LTR = L_T_DELT_USERM(N,V,Q)
                ENDIF
                DO O1 = 1, NELEMENTS
                 L_DB_CUMSOURCE(V,O1) = TR * L_DB_CUMSOURCE(V,O1) &
                        + LTR * DB_CUMSOURCE(V,O1,NC-1)
                ENDDO
              ENDDO
            ENDIF

!  END LAYER LOOP

          ENDDO

!  OFFGRID OUTPUT
!  -------------

!  REQUIRE PARTIAL LAYER TRANSMITTANCE
!  REQUIRE ADDITIONAL LINEARIZATION

          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)

            IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
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
                    DO O1 = 1, NELEMENTS
                      COLUMNWF_DB(Q,UTA,V,O1) = TR*L_DB_CUMSOURCE(V,O1) &
                                             + LTR*DB_CUMSOURCE(V,O1,NC)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ELSE
              DO V = 1, N_GEOMETRIES
                IF ( DO_SSCORR_OUTGOING ) THEN
                 TR  = UP_LOSTRANS_UT(UT,V)
                 LTR = L_UP_LOSTRANS_UT(UT,Q,V)
                ELSE
                 TR  = T_UTUP_USERM(UT,V)
                 LTR = L_T_UTUP_USERM(UT,V,Q)
                ENDIF
                DO O1 = 1, NELEMENTS
                  COLUMNWF_DB(Q,UTA,V,O1) = TR*L_DB_CUMSOURCE(V,O1) &
                                         + LTR*DB_CUMSOURCE(V,O1,NC)
                ENDDO
              ENDDO
            ENDIF

!  ONGRID OUTPUT : SET FINAL CUMULATIVE SOURCE DIRECTLY

          ELSE

            DO V = 1, N_GEOMETRIES
              DO O1 = 1, NELEMENTS
                COLUMNWF_DB(Q,UTA,V,O1) = L_DB_CUMSOURCE(V,O1)
              ENDDO
            ENDDO

          ENDIF

!  CHECK FOR UPDATING THE RECURSION

          IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT

!  END OPTICAL DEPTH LOOP

        ENDDO

!  END WEIGHTING FUNCTION LOOP

      ENDDO

!  FINISH

      RETURN
      END SUBROUTINE VLIDORT_LAC_DBCORRECTION

!

      SUBROUTINE VFO_LCS_MASTER_INTERFACE ( &
        DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,    & ! Input flags
        DO_PLANE_PARALLEL, DO_SSCORR_NADIR, DO_SSCORR_OUTGOING,        & ! Input flags
        DO_DELTAM_SCALING, DO_UPWELLING, DO_DNWELLING,                 & ! Input flags
        DO_LAMBERTIAN_SURFACE, DO_OBSERVATION_GEOMETRY,                & ! Input flags
        DO_COLUMN_LINEARIZATION, DO_SURFACE_LINEARIZATION,             & ! Input flags
        NSTOKES, NLAYERS, NFINELAYERS, NMOMENTS, NSTREAMS,             & ! Input numbers
        NGREEK_MOMENTS_INPUT, N_TOTALCOLUMN_WFS, N_SURFACE_WFS,        & ! Input numbers
        N_SZANGLES, SZANGLES, N_USER_VZANGLES, USER_VZANGLES,          & ! Input geometry
        N_USER_RELAZMS, USER_RELAZMS,                                  & ! Input geometry
        N_USER_LEVELS, USER_LEVELS, EARTH_RADIUS, HEIGHT_GRID,         & ! Input other
        SS_FLUX_MULTIPLIER, FLUXVEC,                                   & ! Inputs (Optical - Atmos)
        DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,    & ! Inputs (Optical - Atmos)
        DELTAU_VERT, OMEGA_TOTAL, GREEKMAT_TOTAL, TRUNC_FACTOR,        & ! Inputs (Optical - Atmos)
        THERMAL_BB_INPUT,                                              & ! Inputs (Optical - Atmos)
        LAMBERTIAN_ALBEDO, EXACTDB_BRDFUNC, SURFBB, USER_EMISSIVITY,   & ! Inputs (Optical - Surf)
        L_DELTAU_VERT_INPUT, L_OMEGA_TOTAL_INPUT,                      & ! Inputs (Optical - Lin Atmos)
        L_GREEKMAT_TOTAL_INPUT, L_DELTAU_VERT, L_TRUNC_FACTOR,         & ! Inputs (Optical - Lin Atmos)
        LS_EXACTDB_BRDFUNC, LS_USER_EMISSIVITY,                        & ! Inputs (Optical - Lin Surf)
        FO_STOKES_SS, FO_STOKES_DB, FO_STOKES_DTA, FO_STOKES_DTS,      & ! Output - Std
        FO_STOKES,                                                     & ! Output - Std
        FO_COLUMNWF_SS,  FO_COLUMNWF_DB,  FO_SURFACEWF_DB,             & ! Output - Lin
        FO_COLUMNWF_DTA, FO_COLUMNWF_DTS, FO_SURFACEWF_DTS,            & ! Output - Lin
        FO_COLUMNWF, FO_SURFACEWF,                                     & ! Output - Lin
        FAIL, MESSAGE, TRACE )                                           ! Output

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Parameter argument
!  ------------------

      INTEGER, PARAMETER  :: FFP = SELECTED_REAL_KIND(15)

!  Inputs
!  ======

      LOGICAL, INTENT(IN) ::            DO_SOLAR_SOURCES
      LOGICAL, INTENT(IN) ::            DO_THERMAL_EMISSION
      LOGICAL, INTENT(IN) ::            DO_SURFACE_EMISSION
      LOGICAL, INTENT(IN) ::            DO_PLANE_PARALLEL
      LOGICAL, INTENT(IN) ::            DO_SSCORR_NADIR
      LOGICAL, INTENT(IN) ::            DO_SSCORR_OUTGOING
      LOGICAL, INTENT(IN) ::            DO_DELTAM_SCALING
      LOGICAL, INTENT(IN) ::            DO_UPWELLING
      LOGICAL, INTENT(IN) ::            DO_DNWELLING
      LOGICAL, INTENT(IN) ::            DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT(IN) ::            DO_OBSERVATION_GEOMETRY
      LOGICAL, INTENT(IN) ::            DO_COLUMN_LINEARIZATION
      LOGICAL, INTENT(IN) ::            DO_SURFACE_LINEARIZATION

      INTEGER, INTENT(IN) ::            NSTOKES
      INTEGER, INTENT(IN) ::            NLAYERS
      INTEGER, INTENT(IN) ::            NFINELAYERS

      INTEGER, INTENT(IN) ::            NMOMENTS
      INTEGER, INTENT(IN) ::            NSTREAMS
      INTEGER, INTENT(IN) ::            NGREEK_MOMENTS_INPUT
      INTEGER, INTENT(IN) ::            N_TOTALCOLUMN_WFS
      INTEGER, INTENT(IN) ::            N_SURFACE_WFS

      INTEGER, INTENT(IN) ::            N_SZANGLES
      DOUBLE PRECISION, INTENT(IN) ::   SZANGLES ( MAX_SZANGLES )
      INTEGER, INTENT(IN) ::            N_USER_VZANGLES
      DOUBLE PRECISION, INTENT(IN) ::   USER_VZANGLES ( MAX_USER_VZANGLES )
      INTEGER, INTENT(IN) ::            N_USER_RELAZMS
      DOUBLE PRECISION, INTENT(IN) ::   USER_RELAZMS  ( MAX_USER_RELAZMS )
      INTEGER, INTENT(IN) ::            N_USER_LEVELS
      DOUBLE PRECISION, INTENT(IN) ::   USER_LEVELS ( MAX_USER_LEVELS )

      DOUBLE PRECISION, INTENT(IN) ::   EARTH_RADIUS
      DOUBLE PRECISION, INTENT(IN) ::   HEIGHT_GRID ( 0:MAXLAYERS )

      DOUBLE PRECISION, INTENT(IN) ::   SS_FLUX_MULTIPLIER
      DOUBLE PRECISION, INTENT(IN) ::   FLUXVEC ( MAXSTOKES )

      DOUBLE PRECISION, INTENT(IN) ::   DELTAU_VERT_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) ::   OMEGA_TOTAL_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) ::   GREEKMAT_TOTAL_INPUT &
          ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )

      DOUBLE PRECISION, INTENT(IN) ::   DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) ::   OMEGA_TOTAL ( MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) ::   GREEKMAT_TOTAL &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION, INTENT(IN) ::   TRUNC_FACTOR ( MAXLAYERS )

      DOUBLE PRECISION, INTENT(IN) ::   THERMAL_BB_INPUT ( 0:MAXLAYERS )

      DOUBLE PRECISION, INTENT(IN) ::   LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT(IN) ::   EXACTDB_BRDFUNC &
          ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

      DOUBLE PRECISION, INTENT(IN) ::   SURFBB
      DOUBLE PRECISION, INTENT(IN) ::   USER_EMISSIVITY &
          ( MAXSTOKES, MAX_USER_STREAMS )

      DOUBLE PRECISION, INTENT(IN) ::   L_DELTAU_VERT_INPUT &
          ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) ::   L_OMEGA_TOTAL_INPUT &
          ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) ::   L_GREEKMAT_TOTAL_INPUT &
          ( MAX_ATMOSWFS, 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION, INTENT(IN) ::   L_DELTAU_VERT &
          ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) ::   L_TRUNC_FACTOR &
          ( MAX_ATMOSWFS, MAXLAYERS )

      DOUBLE PRECISION, INTENT(IN) ::   LS_EXACTDB_BRDFUNC &
          ( MAX_SURFACEWFS, MAXSTOKES_SQ, &
            MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT(IN) ::   LS_USER_EMISSIVITY &
          ( MAX_SURFACEWFS, MAXSTOKES, MAX_USER_STREAMS )

!  Outputs
!  =======

!  Standard
!  --------

!  Solar

      DOUBLE PRECISION, INTENT (INOUT)  :: FO_STOKES_SS &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT)  :: FO_STOKES_DB &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  Thermal

      DOUBLE PRECISION, INTENT (INOUT)  :: FO_STOKES_DTA &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT)  :: FO_STOKES_DTS &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  Composite

      DOUBLE PRECISION, INTENT (INOUT)  :: FO_STOKES &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )

!  Linearized
!  ----------

!  Solar

      DOUBLE PRECISION, INTENT (INOUT) :: FO_COLUMNWF_SS &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS,&
            MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) :: FO_COLUMNWF_DB &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) :: FO_SURFACEWF_DB &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  Thermal

      DOUBLE PRECISION, INTENT (INOUT) :: FO_COLUMNWF_DTA &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS,&
            MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) :: FO_COLUMNWF_DTS &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) :: FO_SURFACEWF_DTS &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  Composite

      DOUBLE PRECISION, INTENT (INOUT) :: FO_COLUMNWF &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS,&
            MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) :: FO_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  Error Handling

      LOGICAL, INTENT (OUT)             :: FAIL
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  Local variables
!  ===============

!  Max dimensions
!  --------------

      INTEGER   :: MAXGEOMS, MAXFINE!, MAXLAYERS, MAXMOMENTS_INPUT, MAX_USER_LEVELS

!  Dimensions
!  ----------

!  Layer and geometry control. Finelayer divisions may be changed

      INTEGER   :: NGEOMS!, NLAYERS
      INTEGER   :: NFINEDIVS (MAXLAYERS, MAX_GEOMETRIES )
      INTEGER   :: NMOMENTS_INPUT
      !INTEGER   :: N_USER_LEVELS

!  Number of Stokes components

      !INTEGER   :: NSTOKES

!  Number of column & surface weighting functions

      INTEGER   :: N_COLUMNWFS
      INTEGER   :: N_SURFACEWFS

!  Configuration inputs
!  --------------------

!  Sources control, including thermal

      !LOGICAL   :: DO_SOLAR_SOURCES
      !LOGICAL   :: DO_THERMAL_EMISSION
      !LOGICAL   :: DO_SURFACE_EMISSION

!  Flags (sphericity flags should be mutually exclusive)

      LOGICAL   :: DO_PLANPAR
      LOGICAL   :: DO_REGULAR_PS
      LOGICAL   :: DO_ENHANCED_PS

!  Delta-m scaling flag

      !LOGICAL   :: DO_DELTAM_SCALING

!  Directional Flags

      !LOGICAL   :: DO_UPWELLING, DO_DNWELLING

!  Lambertian surface flag

      LOGICAL   :: DO_LAMBERTIAN

!  Vector sunlight flag

      LOGICAL   :: DO_SUNLIGHT

!  Linearization flags

      LOGICAL   :: DO_COLUMNWFS, DO_SURFACEWFS

!  General inputs
!  --------------

!  DTR = degrees-to-Radians

      REAL(FFP) :: DTR!, VSIGN, PIE

!  Critical adjustment for cloud layers

      LOGICAL   :: DoCrit
      REAL(FFP) :: Acrit

!  Earth radius + heights

      REAL(FFP) :: ERADIUS
      REAL(FFP) :: HEIGHTS ( 0:MAXLAYERS )

!  Output levels

      INTEGER   :: FO_USER_LEVELS ( MAX_USER_LEVELS )

!  Geometry inputs
!  ---------------

!  Input angles (Degrees)

      REAL(FFP) :: theta_boa ( MAX_GEOMETRIES ) !SZA
      REAL(FFP) :: alpha_boa ( MAX_GEOMETRIES ) !UZA
      REAL(FFP) :: phi_boa   ( MAX_GEOMETRIES ) !RAA

!     Mu0 = cos(theta_boa), required for surface term (both regular & enhanced)
!     Mu1 = cos(alpha_boa), required for the Regular PS only

      REAL(FFP) :: Mu0 ( MAX_GEOMETRIES )
      REAL(FFP) :: Mu1 ( MAX_GEOMETRIES )

!  Flag for the Nadir case. Intent(inout), input if DO_LOSpaths set

      LOGICAL   :: DoNadir ( MAX_GEOMETRIES )

!  Optical inputs
!  --------------

!  Solar flux

      REAL(FFP) :: FLUX
      !REAL(FFP) :: FLUXVEC ( MAXSTOKES )

!  Atmosphere

      REAL(FFP) :: EXTINCTION  ( MAXLAYERS )
      REAL(FFP) :: DELTAUS     ( MAXLAYERS )
      REAL(FFP) :: OMEGA       ( MAXLAYERS )
      REAL(FFP) :: GREEKMAT    ( MAXLAYERS, 0:MAXMOMENTS_INPUT, MAXSTOKES, MAXSTOKES )

!  For TMS correction

      REAL(FFP) :: TRUNCFAC    ( MAXLAYERS )

!  Thermal inputs

      REAL(FFP) :: BB_INPUT ( 0:MAXLAYERS )

!  Surface properties - reflective (could be the albedo)

      REAL(FFP) :: REFLEC ( MAXSTOKES, MAXSTOKES, MAX_GEOMETRIES )

!  Surface properties - emissive

      !REAL(FFP) :: SURFBB
      REAL(FFP) :: EMISS ( MAX_GEOMETRIES )

!  Linearized inputs
!  -----------------

!  Linearization control

      LOGICAL   :: LVARYFLAGS ( MAXLAYERS )
      INTEGER   :: LVARYNUMS  ( MAXLAYERS )
      LOGICAL   :: LVARYMOMS  ( MAXLAYERS, MAX_ATMOSWFS )

!  Linearized optical inputs

      REAL(FFP) :: L_EXTINCTION ( MAXLAYERS, MAX_ATMOSWFS )
      REAL(FFP) :: L_DELTAUS    ( MAXLAYERS, MAX_ATMOSWFS )
      REAL(FFP) :: L_OMEGA      ( MAXLAYERS, MAX_ATMOSWFS )
      REAL(FFP) :: L_GREEKMAT   ( MAXLAYERS, 0:MAXMOMENTS_INPUT, &
                                  MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )

!  Linearized TMS correction

      REAL(FFP) :: L_TRUNCFAC ( MAXLAYERS, MAX_ATMOSWFS )

!  Surface properties - reflective

      REAL(FFP) :: LS_REFLEC ( MAXSTOKES, MAXSTOKES, &
                               MAX_GEOMETRIES, MAX_SURFACEWFS )

!  Surface properties - emissive

      REAL(FFP) :: LS_EMISS  ( MAX_GEOMETRIES, MAX_SURFACEWFS )

!  Local variables
!  ---------------

      INTEGER   :: G, GK, L, N, NS2, O1, O2, PAR, SPAR
      INTEGER   :: LUM=1, LUA=1
      REAL(FFP) :: DNM1

!  Define VFO_MASTER inputs
!  ========================

!  Note: argument passing for variables defined elsewhere
!        commented out here

!  Max dimensions

      maxgeoms          = MAX_GEOMETRIES
      !maxlayers        =
      maxfine           = MAXFINELAYERS
      !maxmoments_input =
      !max_user_levels  =

      !max_atmoswfs      =
      !max_surfacewfs    =

!  Dimensions

      !Note: set for ObsGeo mode only at present
      ngeoms          = N_SZANGLES

      !nlayers        =
      nfinedivs       = NFINELAYERS

      !recall NMOMENTS = MIN(2*NSTREAMS-1, NGREEK_MOMENTS_INPUT)
      !nmoments_input  = NMOMENTS
      nmoments_input  = NGREEK_MOMENTS_INPUT

      !n_user_levels  =
      !nstokes        =

      n_columnwfs     = N_TOTALCOLUMN_WFS
      n_surfacewfs    = N_SURFACE_WFS

!  Configuration inputs

      !do_solar_sources    =
      !do_thermal_emission =
      !do_surface_emission =

      if (DO_PLANE_PARALLEL) then
        do_planpar     = .TRUE.
        do_regular_ps  = .FALSE.
        do_enhanced_ps = .FALSE.
      else
        do_planpar = .FALSE.
        if (DO_SSCORR_NADIR) then
          do_regular_ps  = .TRUE.
          do_enhanced_ps = .FALSE.
        else if (DO_SSCORR_OUTGOING) then
          do_regular_ps  = .FALSE.
          do_enhanced_ps = .TRUE.
        end if
      end if

      !do_deltam_scaling  =
      !do_upwelling       =
      !do_dnwelling       =

      do_lambertian      = DO_LAMBERTIAN_SURFACE
      do_sunlight        = .TRUE. !set for now

      do_columnwfs       = DO_COLUMN_LINEARIZATION
      do_surfacewfs      = DO_SURFACE_LINEARIZATION

!  General inputs

      dtr    = DEG_TO_RAD
      !Pie   =
      DoCrit = .FALSE. !set for now
      Acrit  = ZERO    !set for now

      !defined within VFO_MASTER now
      !if ( do_upwelling ) then
      !  vsign =  one
      !else if ( do_dnwelling ) then
      !  vsign = -one
      !end if

      eradius             = EARTH_RADIUS
      heights(0:nlayers)  = HEIGHT_GRID(0:NLAYERS)

      fo_user_levels(1:n_user_levels) = &
        INT(user_levels(1:n_user_levels))

!  Geometry inputs

      !Note: set for ObsGeo mode only at present

      !angles in deg
      theta_boa(1:ngeoms) = SZANGLES(1:N_SZANGLES)
      alpha_boa(1:ngeoms) = USER_VZANGLES(1:N_USER_VZANGLES)
      phi_boa  (1:ngeoms) = USER_RELAZMS (1:N_USER_RELAZMS)

      !For regular & enhanced PS
      Mu0(1:ngeoms) = cos( theta_boa(1:ngeoms) * dtr ) !cos(SZA)
      !For regular PS only
      Mu1(1:ngeoms) = cos( alpha_boa(1:ngeoms) * dtr ) !cos(UZA)

      !For enhanced PS only
      do g=1,ngeoms
        if (alpha_boa(g) == ZERO) then
          DoNadir(g) = .TRUE.
        else
          DoNadir(g) = .FALSE.
        end if
      end do

!  Optical inputs

      flux     = SS_FLUX_MULTIPLIER !=FLUX_FACTOR/(4.0d0*PIE)
      !fluxvec

      if (do_deltam_scaling) then
        extinction(1:nlayers) = DELTAU_VERT(1:NLAYERS)/&
          (HEIGHT_GRID(0:NLAYERS-1)-HEIGHT_GRID(1:NLAYERS))
        deltaus(1:nlayers)    = DELTAU_VERT(1:NLAYERS)
      else
        extinction(1:nlayers) = DELTAU_VERT_INPUT(1:NLAYERS)/&
          (HEIGHT_GRID(0:NLAYERS-1)-HEIGHT_GRID(1:NLAYERS))
        deltaus(1:nlayers)    = DELTAU_VERT_INPUT(1:NLAYERS)
      endif

      omega(1:nlayers)      = OMEGA_TOTAL_INPUT(1:NLAYERS)
      do o1=1,nstokes
        do o2=1,nstokes
          gk = 4*(o1-1) + o2
          do l=0,nmoments_input
            do n=1,nlayers
              greekmat(n,l,o1,o2) = GREEKMAT_TOTAL_INPUT(L,N,GK)
            end do
          end do
        end do
      end do

      !For TMS correction only
      if (do_deltam_scaling) then
        do n=1,nlayers
          truncfac(n) = TRUNC_FACTOR(N)
        end do
      end if

      bb_input(0:nlayers) = THERMAL_BB_INPUT(0:NLAYERS)

      reflec = ZERO
      if (do_lambertian) then
        reflec(1,1,1:ngeoms) = LAMBERTIAN_ALBEDO
      else
        !Note: set for ObsGeo mode only at present
        do g=1,ngeoms
          do o1=1,nstokes
            do o2=1,nstokes
              gk = 4*(o1-1) + o2
              reflec(o1,o2,g) = EXACTDB_BRDFUNC(GK,LUM,LUA,G)
            end do
          end do
        end do
      end if

      !surfbb =
      do g=1,ngeoms
        emiss(g) = USER_EMISSIVITY(1,G)
      end do

!  Linearized control

      !Default setting: see below for modification when
      !Greek moments vary
      Lvarymoms(1:nlayers,1:max_atmoswfs) = .FALSE.

!  Linearized optical

      !Atmospheric quantities
      L_extinction = ZERO
      L_deltaus    = ZERO
      L_omega      = ZERO
      L_truncfac   = ZERO
      L_greekmat   = ZERO

      !Recall that VLIDORT takes FULLY NORMALIZED linearized
      !atmospheric inputs
      !                      (x/y)*(dy/dx)
      !whereas the FO code only takes PARTIALLY NORMALIZED linearized
      !atmospheric inputs
      !                        x*(dy/dx)
      !thus, we must compensate for this difference when formulating the
      !inputs here!
      do n=1,nlayers
        if (do_columnwfs) then
          do par=1,n_columnwfs

            if (do_deltam_scaling) then
              L_extinction(n,par) = (DELTAU_VERT(N)/&
                (HEIGHT_GRID(N-1)-HEIGHT_GRID(N)))*L_DELTAU_VERT(PAR,N)
              L_deltaus(n,par)    = DELTAU_VERT(N)*L_DELTAU_VERT(PAR,N)
            else
              L_extinction(n,par) = (DELTAU_VERT_INPUT(N)/&
                (HEIGHT_GRID(N-1)-HEIGHT_GRID(N)))*L_DELTAU_VERT_INPUT(PAR,N)
              L_deltaus(n,par)    = DELTAU_VERT_INPUT(N)*L_DELTAU_VERT_INPUT(PAR,N)
            endif

            L_omega(n,par)      = OMEGA_TOTAL_INPUT(N)*L_OMEGA_TOTAL_INPUT(PAR,N)

            if (do_deltam_scaling) &
              L_truncfac(n,par) = L_TRUNC_FACTOR(PAR,N)

            do o1=1,nstokes
              do o2=1,nstokes
                gk = 4*(o1-1) + o2
                do l=0,nmoments_input
                  L_greekmat(n,l,o1,o2,par) = &
                    GREEKMAT_TOTAL_INPUT(L,N,GK) &
                    *L_GREEKMAT_TOTAL_INPUT(PAR,L,N,GK)

                  !Check for variation of Greek moments associated with
                  !Jacobian wrt current atmospheric parameter
                  if ( (l==1) .and. (gk==1) .and. &
                       (L_GREEKMAT_TOTAL_INPUT(PAR,L,N,GK) > 1.0d-8) ) &
                    Lvarymoms(n,par) = .TRUE.
                end do
              end do
            end do

          end do
        end if
      end do

      !Surface quantities
      LS_reflec = ZERO
      LS_emiss  = ONE
      if (do_surfacewfs) then
        if (do_lambertian) then
          do spar=1,n_surfacewfs
            do g=1,ngeoms
              LS_reflec(1,1,1:ngeoms,1:n_surfacewfs) = ONE
            end do
          end do
        else
          !Note: set for ObsGeo mode only at present
          do spar=1,n_surfacewfs
            do g=1,ngeoms
              do o1=1,nstokes
                do o2=1,nstokes
                  gk = 4*(o1-1) + o2
                  LS_reflec(o1,o2,g,spar) = LS_EXACTDB_BRDFUNC(SPAR,GK,LUM,LUA,G)
                end do
              end do
            end do
          end do
        end if

        do spar=1,n_surfacewfs
          do g=1,ngeoms
            LS_emiss(g,spar) = LS_USER_EMISSIVITY(SPAR,1,G)
          end do
        end do
      end if

!  Call VFO_MASTER
!  ===============

      CALL VFO_LCS_MASTER &
       ( maxgeoms, maxlayers, maxfine, maxmoments_input, max_user_levels,    & ! Input max dims
         max_atmoswfs, max_surfacewfs,                                       & ! Input max dims
         ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, nstokes, & ! Input dims
         n_columnwfs, n_surfacewfs,                                          & ! Input dims
         do_solar_sources, do_thermal_emission, do_surface_emission,         & ! Input flags
         do_planpar, do_regular_ps, do_enhanced_ps, do_deltam_scaling,       & ! Input flags
         do_upwelling, do_dnwelling, do_lambertian, do_sunlight,             & ! Input flags
         do_columnwfs, do_surfacewfs,                                        & ! Input flags
         dtr, Pie, doCrit, Acrit, eradius, heights, fo_user_levels,          & ! Input general
         theta_boa, alpha_boa, phi_boa, Mu0, Mu1, doNadir,                   & ! Input geometry
         flux, fluxvec, extinction, deltaus, omega, greekmat, truncfac,      & ! Input - std atmos optical
         bb_input, reflec, surfbb, emiss,                                    & ! Input - std surf optical
         Lvarymoms,                                                          & ! Input - lin control
         L_extinction, L_deltaus, L_omega, L_greekmat, L_truncfac,           & ! Input - lin atmos optical
         LS_reflec, LS_emiss,                                                & ! Input - lin surf optical
         fo_stokes_ss, fo_stokes_db, fo_stokes_dta, fo_stokes_dts,           & ! Output - std
         fo_stokes,                                                          & ! Output - std
         fo_columnwf_ss,  fo_columnwf_db,  fo_surfacewf_db,                  & ! Output - lin
         fo_columnwf_dta, fo_columnwf_dts, fo_surfacewf_dts,                 & ! Output - lin
         fo_columnwf, fo_surfacewf,                                          & ! Output - lin
         fail, message, trace )                                                ! Output - error

      END SUBROUTINE VFO_LCS_MASTER_INTERFACE

      END MODULE vlidort_lc_corrections

