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
! #                   2.5, 2.6                                  #
! #  Release Date :   December 2005  (2.0)                      #
! #  Release Date :   March 2007     (2.2)                      #
! #  Release Date :   October 2007   (2.3)                      #
! #  Release Date :   December 2008  (2.4)                      #
! #  Release Date :   April 2009     (2.4R)                     #
! #  Release Date :   July 2009      (2.4RT)                    #
! #  Release Date :   October 2010   (2.4RTC)                   #
! #  Release Date :   March 2011     (2.5)                      #
! #  Release Date :   May 2012       (2.6)                      #
! #                                                             #
! #       NEW: TOTAL COLUMN JACOBIANS         (2.4)             #
! #       NEW: BPDF Land-surface KERNELS      (2.4R)            #
! #       NEW: Thermal Emission Treatment     (2.4RT)           #
! #       Consolidated BRDF treatment         (2.4RTC)          #
! #       f77/f90 Release                     (2.5)             #
! #       External SS / New I/O Structures    (2.6)             #
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
! #            VLIDORT_LP_SSCORR_NADIR (master, renamed)        #
! #            VLIDORT_LC_SSCORR_NADIR (master, new 2.4)        #
! #              VLIDORTSS_L_ZMATRICES                          #
! #              VLIDORTSS_L_FMATRICES_REFR (placeholder)       #
! #                                                             #
! #   2. Outgoing SS-correction                                 #
! #      Version 2.2   Whole layer integration                  #
! #              2.3.  partial-layer integration                #
! #      Version 2.4.  Column Jacobians introduced              #
! #                                                             #
! #            VLIDORT_L_SSCORR_OUTGOING (master)               #
! #              L_OUTGOING_INTEGRATION_UP                      #
! #              L_OUTGOING_INTEGRATION_DN                      #
! #                                                             #
! #   3. DB correction, Lambertian and BRDF surfaces            #
! #      Version 2.4.  LAC module (column Jacobians)            #
! #                                                             #
! #            VLIDORT_LAP_DBCORRECTION                         #
! #            VLIDORT_LAC_DBCORRECTION                         #
! #            VLIDORT_LS_DBCORRECTION                          #
! #                                                             #
! #      Version 2.4R ------- Notes -------                     #
! #                                                             #
! #            Additional correction to ZMATRIX and FMATRIX     #
! #            routines for Azimuth > 180. HvdM, Eqs(94/95)     #
! #            Implemented by V. Natraj and R. Spurr, 5/1/09    #
! #                                                             #
! ###############################################################


      MODULE vlidort_l_corrections

      PRIVATE
      PUBLIC :: VLIDORT_LP_SSCORR_NADIR, &
                VLIDORT_LC_SSCORR_NADIR, &
                VLIDORT_L_SSCORR_OUTGOING, &
                VLIDORT_LAP_DBCORRECTION, &
                VLIDORT_LAC_DBCORRECTION, &
                VLIDORT_LS_DBCORRECTION

      CONTAINS

      SUBROUTINE VLIDORT_LP_SSCORR_NADIR ( &
        SSFLUX, LAYER_TO_VARY, NV_PARAMETERS, &
        DO_SSCORR_TRUNCATION, DO_REFRACTIVE_GEOMETRY, &
        DO_DELTAM_SCALING, DO_UPWELLING, &
        DO_DNWELLING, NSTOKES, &
        NLAYERS, NGREEK_MOMENTS_INPUT, &
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
        L_EMULT_UP, L_EMULT_DN, &
        L_UT_EMULT_UP, L_UT_EMULT_DN, &
        PROFILEWF_SS )

!  SINGLE SCATTER EXACT CALCULATION. NADIR VIEW
!   PROGRAMMED BY R. SPURR, RT SOLUTIONS INC.
!    SECOND DRAFT, APRIL 14TH 2005.
!    THIRD DRAFT,  MAY    6TH 2005.
!       - ADDITIONAL CODE TO DEAL WITH REFRACTION.

      USE VLIDORT_PARS

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT (IN) :: SSFLUX
      INTEGER, INTENT (IN) ::          LAYER_TO_VARY
      INTEGER, INTENT (IN) ::          NV_PARAMETERS
      LOGICAL, INTENT (IN) ::          DO_SSCORR_TRUNCATION
      LOGICAL, INTENT (IN) ::          DO_REFRACTIVE_GEOMETRY
      LOGICAL, INTENT (IN) ::          DO_DELTAM_SCALING
      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_DNWELLING
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
      DOUBLE PRECISION, INTENT (IN) :: L_EMULT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_EMULT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UT_EMULT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, 0:MAXLAYERS, &
            MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UT_EMULT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, 0:MAXLAYERS, &
            MAXBEAMS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (INOUT) :: PROFILEWF_SS &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
            MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )

!  LOCAL VARIABLES
!  ---------------

!  INDICES

      INTEGER ::          N, NUT, NSTART, NUT_PREV, NLEVEL, V
      INTEGER ::          UT, UTA, UM, IB, IA, NC, Q, NM1, NV, O1, O2

!  OTHER HELP VARIABLES

      DOUBLE PRECISION :: HELP, TTT, UVAR, AVAR, FT1, FT2, DNM1
      DOUBLE PRECISION :: L_FINAL_SOURCE, L_SS_LAYERSOURCE
      DOUBLE PRECISION :: L_SSCORRECTION, L_TR_CUMSOURCE
      DOUBLE PRECISION :: VAR_TMS(MAX_ATMOSWFS), MUX, SUM, VIEWSIGN

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
!    - DROP THESE VARIABLES WHEN FURTHER TESTING HAS BEEN DONE

      LOGICAL, PARAMETER :: DO_SUNLIGHT = .TRUE.

!  SET UP OPERATIONS
!  -----------------

      NV = LAYER_TO_VARY

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
            VAR_TMS(Q) = VAR_TMS(Q) - FT1 * L_SSFDEL(NV,Q)
          ENDIF
        ENDDO
      ENDIF

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

      DO UM = 1, N_USER_STREAMS
        CALPHA(UM) = USER_STREAMS(UM)
        SALPHA(UM) = DSQRT ( ONE - CALPHA(UM) * CALPHA(UM) )
      ENDDO

      DO IA = 1, N_USER_RELAZMS
        CPHI(IA) = DCOS ( USER_RELAZMS(IA) * DEG_TO_RAD )
      ENDDO

!  ####################
!  #    UPWELLING     #
!  ####################

      IF ( DO_UPWELLING ) THEN

!  =================================================
!  TOTAL SCATTERING MATRIX LINEARIZATION (UPWELLING)
!  =================================================

        IF ( STERM_LAYERMASK_UP(NV)) THEN

!  GET L_Z MATRICES FOR EACH LAYER, IF THE SCATTER LAW IS LINEARIZED

         VIEWSIGN = -1.0D0

         CALL VLIDORTSS_L_ZMATRICES &
        ( DO_SUNLIGHT, NSTOKES, NGREEK_MOMENTS_INPUT, NV, NV_PARAMETERS, &
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
                        ZMAT_UP(V,NV,O1,O2)    * VAR_TMS(Q)
                ENDDO
              ENDDO
            ELSE
              DO O1 = 1, NSTOKES
                DO O2 = 1, NSTOKES
                  L_ZMAT_UP(Q,V,NV,O1,O2) = &
                        ZMAT_UP(V,NV,O1,O2) * VAR_TMS(Q)
                ENDDO
              ENDDO
            ENDIF

!  END PARAMETER LOOP AND VIEWING DIRECTIONS LOOP

          ENDDO
         ENDDO

!  ONLY IF LAYER NV EXISTS

        ENDIF

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

!  IF N = NV (THE LAYER THAT IS VARYING)

            IF ( N .EQ. NV ) THEN

              DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO Q = 1, NV_PARAMETERS
                  DO O1 = 1, NSTOKES
                   SUM = ZERO
                   DO O2 = 1, NSTOKES
                    HELP = ZMAT_UP(V,N,O1,O2) * L_EMULT_UP(UM,N,NV,IB,Q) &
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

!  VARIATIONS WHEN N > NV

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
                   L_SS_CUMSOURCE(Q,V,O1) = L_SS_LAYERSOURCE &
               +   T_DELT_USERM(N,UM)   * L_SS_CUMSOURCE(Q,V,O1)
                  ENDDO
                 ENDDO
                ENDDO
               ENDDO
              ENDDO

!  TRANSMITTANCE FOR N < NV

            ELSE

              DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO Q = 1, NV_PARAMETERS
                  DO O1 = 1, NSTOKES
                   L_SS_CUMSOURCE(Q,V,O1) = &
                     T_DELT_USERM(N,UM) * L_SS_CUMSOURCE(Q,V,O1)
                  ENDDO
                 ENDDO
                ENDDO
               ENDDO
              ENDDO

            ENDIF

          ENDDO

!  OFFGRID OUTPUT
!  --------------

!  SET FINAL CUMULATIVE SOURCE AND SINGLE SCATTER WEIGHTING FUNCTION

          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)

!  IF N = NV (THE LAYER THAT IS VARYING)
!    ADD LINEARIZATION OF ADDITIONAL PARTIAL LAYER SOURCE TERM =
!        L_EXACT_SCAT(N) * MULTIPLIER  +  EXACT_SCAT * L_MULTIPLIER(N)

            IF ( N .EQ. NV ) THEN

              DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO Q = 1, NV_PARAMETERS
                  DO O1 = 1, NSTOKES
                   SUM = ZERO
                   DO O2 = 1, NSTOKES
                    HELP = &
                  ZMAT_UP(V,N,O1,O2) * L_UT_EMULT_UP(UM,UT,NV,IB,Q) &
              + L_ZMAT_UP(Q,V,N,O1,O2) * UT_EMULT_UP(UM,UT,IB)
                    SUM = SUM + HELP * FLUXVEC(O2)
                   ENDDO
                   L_SS_LAYERSOURCE = SUM
                   L_TR_CUMSOURCE = &
                 + L_T_UTUP_USERM(UT,UM,Q) *   SS_CUMSOURCE_UP(V,O1,NC) &
                 +   T_UTUP_USERM(UT,UM)   * L_SS_CUMSOURCE(Q,V,O1)
                   L_FINAL_SOURCE = L_TR_CUMSOURCE + L_SS_LAYERSOURCE
                   L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                   PROFILEWF_SS(Q,NV,UTA,V,O1,UPIDX) = L_SSCORRECTION
                  ENDDO
                 ENDDO
                ENDDO
               ENDDO
              ENDDO

!  VARIATIONS WHEN N > NV
!    ADD LINEARIZATION OF ADDITIONAL PARTIAL LAYER SOURCE TERM =
!         EXACT_SCAT * L_MULTIPLIER(NV)

            ELSE IF ( N .GT. NV ) THEN

              DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO Q = 1, NV_PARAMETERS
                  DO O1 = 1, NSTOKES
                   SUM = ZERO
                   DO O2 = 1, NSTOKES
                    HELP = &
                      ZMAT_UP(V,N,O1,O2) * L_UT_EMULT_UP(UM,UT,NV,IB,Q)
                    SUM = SUM + HELP * FLUXVEC(O2)
                   ENDDO
                   L_SS_LAYERSOURCE = SUM
                   L_TR_CUMSOURCE = &
                 +   T_UTUP_USERM(UT,UM) * L_SS_CUMSOURCE(Q,V,O1)
                   L_FINAL_SOURCE = L_TR_CUMSOURCE + L_SS_LAYERSOURCE
                   L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                   PROFILEWF_SS(Q,NV,UTA,V,O1,UPIDX) = L_SSCORRECTION
                  ENDDO
                 ENDDO
                ENDDO
               ENDDO
              ENDDO

!  TRANSMITTANCE FOR THE OTHER LAYERS

            ELSE

              DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO Q = 1, NV_PARAMETERS
                  DO O1 = 1, NSTOKES
                   L_FINAL_SOURCE = &
                 +   T_UTUP_USERM(UT,UM)   * L_SS_CUMSOURCE(Q,V,O1)
                   L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                   PROFILEWF_SS(Q,NV,UTA,V,O1,UPIDX) = L_SSCORRECTION
                  ENDDO
                 ENDDO
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
                  PROFILEWF_SS(Q,NV,UTA,V,O1,UPIDX) = L_SSCORRECTION
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

        IF ( STERM_LAYERMASK_DN(NV)) THEN

!  GET L_Z MATRICES FOR EACH LAYER, IF THE SCATTER LAW IS LINEARIZED

         VIEWSIGN = +1.0D0

         CALL VLIDORTSS_L_ZMATRICES &
        ( DO_SUNLIGHT, NSTOKES, NGREEK_MOMENTS_INPUT, NV, NV_PARAMETERS, &
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
                        ZMAT_DN(V,NV,O1,O2) * VAR_TMS(Q)
                ENDDO
              ENDDO
            ELSE
              DO O1 = 1, NSTOKES
                DO O2 = 1, NSTOKES
                  L_ZMAT_DN(Q,V,NV,O1,O2) = &
                        ZMAT_DN(V,NV,O1,O2) * VAR_TMS(Q)
                ENDDO
              ENDDO
            ENDIF

!  END PARAMETER AND GEOMETRY LOOPS

          ENDDO
         ENDDO

!  ONLY IF LAYER NV EXISTS

        ENDIF

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

!  IF N = NV (THE LAYER THAT IS VARYING)

            IF ( N .EQ. NV ) THEN

              DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO Q = 1, NV_PARAMETERS
                  DO O1 = 1, NSTOKES
                   SUM = ZERO
                   DO O2 = 1, NSTOKES
                    HELP = ZMAT_DN(V,N,O1,O2) * L_EMULT_DN(UM,N,NV,IB,Q) &
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

!  VARIATIONS WHEN N > NV

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
                   L_SS_CUMSOURCE(Q,V,O1) = L_SS_LAYERSOURCE &
               +   T_DELT_USERM(N,UM)   * L_SS_CUMSOURCE(Q,V,O1)
                  ENDDO
                 ENDDO
                ENDDO
               ENDDO
              ENDDO

!  TRANSMITTANCE FOR N < NV

            ELSE

              DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO Q = 1, NV_PARAMETERS
                  DO O1 = 1, NSTOKES
                   L_SS_CUMSOURCE(Q,V,O1) = &
                     T_DELT_USERM(N,UM) * L_SS_CUMSOURCE(Q,V,O1)
                  ENDDO
                 ENDDO
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

!  IF N = NV (THE LAYER THAT IS VARYING)
!    ADD LINEARIZATION OF ADDITIONAL PARTIAL LAYER SOURCE TERM =
!        L_EXACT_SCAT(N) * MULTIPLIER  +  EXACT_SCAT * L_MULTIPLIER(N)

            IF ( N .EQ. NV ) THEN

              DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO Q = 1, NV_PARAMETERS
                  DO O1 = 1, NSTOKES
                   SUM = ZERO
                   DO O2 = 1, NSTOKES
                    HELP = &
                  ZMAT_DN(V,N,O1,O2) * L_UT_EMULT_DN(UM,UT,NV,IB,Q) &
              + L_ZMAT_DN(Q,V,N,O1,O2) * UT_EMULT_DN(UM,UT,IB)
                    SUM = SUM + HELP * FLUXVEC(O2)
                   ENDDO
                   L_SS_LAYERSOURCE = SUM
                   L_TR_CUMSOURCE = &
                 + L_T_UTDN_USERM(UT,UM,Q) *   SS_CUMSOURCE_DN(V,O1,NC) &
                 +   T_UTDN_USERM(UT,UM)   * L_SS_CUMSOURCE(Q,V,O1)
                   L_FINAL_SOURCE = L_TR_CUMSOURCE + L_SS_LAYERSOURCE
                   L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                   PROFILEWF_SS(Q,NV,UTA,V,O1,DNIDX) = L_SSCORRECTION
                  ENDDO
                 ENDDO
                ENDDO
               ENDDO
              ENDDO

!  VARIATIONS WHEN N > NV
!    ADD LINEARIZATION OF ADDITIONAL PARTIAL LAYER SOURCE TERM =
!         EXACT_SCAT * L_MULTIPLIER(NV)

            ELSE IF ( N .GT. NV ) THEN

              DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO Q = 1, NV_PARAMETERS
                  DO O1 = 1, NSTOKES
                   SUM = ZERO
                   DO O2 = 1, NSTOKES
                    HELP = &
                  ZMAT_DN(V,N,O1,O2) * L_UT_EMULT_DN(UM,UT,NV,IB,Q)
                    SUM = SUM + HELP * FLUXVEC(O2)
                   ENDDO
                   L_SS_LAYERSOURCE = SUM
                   L_TR_CUMSOURCE = &
                 +   T_UTDN_USERM(UT,UM)   * L_SS_CUMSOURCE(Q,V,O1)
                   L_FINAL_SOURCE = L_TR_CUMSOURCE + L_SS_LAYERSOURCE
                   L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                   PROFILEWF_SS(Q,NV,UTA,V,O1,DNIDX) = L_SSCORRECTION
                  ENDDO
                 ENDDO
                ENDDO
               ENDDO
              ENDDO

!  TRANSMITTANCE FOR THE OTHER LAYERS

            ELSE

              DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO Q = 1, NV_PARAMETERS
                  DO O1 = 1, NSTOKES
                   L_FINAL_SOURCE = &
                  +   T_UTDN_USERM(UT,UM)   * L_SS_CUMSOURCE(Q,V,O1)
                   L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                   PROFILEWF_SS(Q,NV,UTA,V,O1,DNIDX) = L_SSCORRECTION
                  ENDDO
                 ENDDO
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
                  PROFILEWF_SS(Q,NV,UTA,V,O1,DNIDX) = L_SSCORRECTION
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
      END SUBROUTINE VLIDORT_LP_SSCORR_NADIR

!

      SUBROUTINE VLIDORT_LC_SSCORR_NADIR ( &
        SSFLUX, NV_PARAMETERS, &
        DO_SSCORR_TRUNCATION, DO_REFRACTIVE_GEOMETRY, &
        DO_DELTAM_SCALING, DO_UPWELLING, &
        DO_DNWELLING, NSTOKES, &
        NLAYERS, NGREEK_MOMENTS_INPUT, &
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
        L_EMULT_UP, L_EMULT_DN, &
        L_UT_EMULT_UP, L_UT_EMULT_DN, &
        COLUMNWF_SS )

!  SINGLE SCATTER EXACT CALCULATION. NADIR VIEW
!   PROGRAMMED BY R. SPURR, RT SOLUTIONS INC.

!  VERSION 2.4. TOTAL COLUMN JACOBIANS, DECEMBER 2008
!               MODELED AFTER LIDORT CODE VERSION 3.3

      USE VLIDORT_PARS

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT (IN) :: SSFLUX
      INTEGER, INTENT (IN) ::          NV_PARAMETERS
      LOGICAL, INTENT (IN) ::          DO_SSCORR_TRUNCATION
      LOGICAL, INTENT (IN) ::          DO_REFRACTIVE_GEOMETRY
      LOGICAL, INTENT (IN) ::          DO_DELTAM_SCALING
      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_DNWELLING
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
      DOUBLE PRECISION, INTENT (IN) :: L_EMULT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_EMULT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UT_EMULT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, &
            0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UT_EMULT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, &
            0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (INOUT) :: COLUMNWF_SS &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, &
            MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )

!  LOCAL VARIABLES
!  ---------------

!  INDICES

      INTEGER ::          N, NUT, NSTART, NUT_PREV, NLEVEL, V
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

      DO UM = 1, N_USER_STREAMS
        CALPHA(UM) = USER_STREAMS(UM)
        SALPHA(UM) = DSQRT ( ONE - CALPHA(UM) * CALPHA(UM) )
      ENDDO

      DO IA = 1, N_USER_RELAZMS
        CPHI(IA) = DCOS ( USER_RELAZMS(IA) * DEG_TO_RAD )
      ENDDO

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
        ( DO_SUNLIGHT, NSTOKES, NGREEK_MOMENTS_INPUT, NV, NV_PARAMETERS, &
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

            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO Q = 1, NV_PARAMETERS
                  DO O1 = 1, NSTOKES
                   SUM = ZERO
                   DO O2 = 1, NSTOKES
                    HELP = ZMAT_UP(V,N,O1,O2) * L_EMULT_UP(UM,N,0,IB,Q) &
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

            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO Q = 1, NV_PARAMETERS
                  DO O1 = 1, NSTOKES
                   SUM = ZERO
                   DO O2 = 1, NSTOKES
                    HELP = &
                  ZMAT_UP(V,N,O1,O2) * L_UT_EMULT_UP(UM,UT,0,IB,Q) &
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
        ( DO_SUNLIGHT, NSTOKES, NGREEK_MOMENTS_INPUT, NV, NV_PARAMETERS, &
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
            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO Q = 1, NV_PARAMETERS
                  DO O1 = 1, NSTOKES
                   SUM = ZERO
                   DO O2 = 1, NSTOKES
                    HELP = ZMAT_DN(V,N,O1,O2) * L_EMULT_DN(UM,N,0,IB,Q) &
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

            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO Q = 1, NV_PARAMETERS
                  DO O1 = 1, NSTOKES
                   SUM = ZERO
                   DO O2 = 1, NSTOKES
                    HELP = &
                  ZMAT_DN(V,N,O1,O2) * L_UT_EMULT_DN(UM,UT,0,IB,Q) &
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

      SUBROUTINE VLIDORTSS_L_ZMATRICES ( &
        DO_SUNLIGHT, NSTOKES, NGREEKMOMS, NV, NV_PARAMETERS, &
        N_GEOMETRIES, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, &
        DO_SCATMAT_VARIATION, LAYER_MAXMOMENTS, GREEKMAT, L_GREEKMAT, &
        DO_SSCORR_TRUNCATION, SSFDEL, L_SSFDEL, VSIGN, VZA_OFFSETS, &
        CTHETA, STHETA, CALPHA, SALPHA, CPHI, PHI, &
        L_ZMAT )

      USE VLIDORT_PARS

      IMPLICIT NONE

!  INPUT
!  -----

!  CONTROL

      LOGICAL, INTENT (IN) ::          DO_SUNLIGHT
      INTEGER, INTENT (IN) ::          NSTOKES, NGREEKMOMS, NV, NV_PARAMETERS
      INTEGER, INTENT (IN) ::          N_GEOMETRIES
      INTEGER, INTENT (IN) ::          NBEAMS, N_USER_STREAMS, N_USER_RELAZMS

!  +/- SIGN, OFFSETS

      DOUBLE PRECISION, INTENT (IN) :: VSIGN
      INTEGER, INTENT (IN) ::          VZA_OFFSETS &
          ( MAXBEAMS, MAX_USER_STREAMS )

!  SCATTERING INPUT INFORMATION

      LOGICAL, INTENT (IN) ::          DO_SCATMAT_VARIATION &
          ( MAXLAYERS, MAX_ATMOSWFS )
      INTEGER, INTENT (IN) ::          LAYER_MAXMOMENTS ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: GREEKMAT &
          ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION, INTENT (IN) :: L_GREEKMAT &
          ( MAX_ATMOSWFS, 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )

!  TRUNCATION CONTROL

      LOGICAL, INTENT (IN) ::          DO_SSCORR_TRUNCATION
      DOUBLE PRECISION, INTENT (IN) :: SSFDEL   ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_SSFDEL ( MAXLAYERS, MAX_ATMOSWFS )

!  ZENITH ANGLE COSINES/SINES, AZIMUTH ANGLE COSINES

      DOUBLE PRECISION, INTENT (IN) :: CTHETA (MAXLAYERS,MAXBEAMS)
      DOUBLE PRECISION, INTENT (IN) :: STHETA (MAXLAYERS,MAXBEAMS)
      DOUBLE PRECISION, INTENT (IN) :: CALPHA (MAX_USER_STREAMS)
      DOUBLE PRECISION, INTENT (IN) :: SALPHA (MAX_USER_STREAMS)
      DOUBLE PRECISION, INTENT (IN) :: CPHI   (MAX_USER_RELAZMS)

!  AZIMUTH ANGLES. ADDED V2.4R. FOR PHI > 180 CASE

      DOUBLE PRECISION, INTENT (IN) :: PHI    (MAX_USER_RELAZMS)

!  OUTPUT
!  ------

!  Z-MATRIX LINEARIZED

      DOUBLE PRECISION, INTENT (INOUT) :: L_ZMAT &
          ( MAX_ATMOSWFS, MAX_GEOMETRIES, MAXLAYERS, 4, 4 )

!  LOCAL
!  -----

!  ROTATION ANGLE COSINES/SINES

      DOUBLE PRECISION :: C1 (MAX_GEOMETRIES)
      DOUBLE PRECISION :: S1 (MAX_GEOMETRIES)
      DOUBLE PRECISION :: C2 (MAX_GEOMETRIES)
      DOUBLE PRECISION :: S2 (MAX_GEOMETRIES)

!  LOCAL F-MATRIX

      DOUBLE PRECISION :: L_FMAT ( MAX_ATMOSWFS, MAX_GEOMETRIES, 6 )

!  HELP VARIABLES

      INTEGER          :: IB, UM, IA, V, Q, GREEKMAT_INDEX(6)
      DOUBLE PRECISION :: COSSCAT, SINSCAT, HELP_SINSCAT
      DOUBLE PRECISION :: CSIG1, CSIG2, SSIG1, SSIG2, CSIG1_2, CSIG2_2

      DOUBLE PRECISION :: UUU(MAX_GEOMETRIES)

!  GENERALIZED SPHERICAL FUNCTIONS
!    P2P2 AND P2M2 ADDED 20 MARCH 2006

      DOUBLE PRECISION :: P00(MAX_GEOMETRIES,2)
      DOUBLE PRECISION :: P02(MAX_GEOMETRIES,2)
      DOUBLE PRECISION :: P2P2(MAX_GEOMETRIES,2)
      DOUBLE PRECISION :: P2M2(MAX_GEOMETRIES,2)

      INTEGER          :: K, L, LNEW, LOLD, ITMP
      INTEGER          :: INDEX_11, INDEX_12, INDEX_34
      INTEGER          :: INDEX_22, INDEX_33, INDEX_44
      DOUBLE PRECISION :: DL, QROOT6, FAC1, FAC2, SQL4, SQL41
      DOUBLE PRECISION :: TMP1, TMP2, SUM23, DIF23
      DOUBLE PRECISION :: FL2, FLL1, PERLL4, WFACT, QFAC, DL1
      DOUBLE PRECISION :: HELP1, HELP2, FDNL1, DNL1, FACT
      DOUBLE PRECISION :: HELP2C1, HELP2S1, HELP3C1, HELP3S1

      DOUBLE PRECISION :: GK11, L_GK11(MAX_ATMOSWFS)
      DOUBLE PRECISION :: GK12, L_GK12(MAX_ATMOSWFS)
      DOUBLE PRECISION :: GK44, L_GK44(MAX_ATMOSWFS)
      DOUBLE PRECISION :: GK34, L_GK34(MAX_ATMOSWFS)
      DOUBLE PRECISION :: GK22, L_GK22(MAX_ATMOSWFS)
      DOUBLE PRECISION :: GK33, L_GK33(MAX_ATMOSWFS)

!  INDEXING KEY

!      GREEKMAT_INDEX(1) = 1  ---> INDEX_11
!      GREEKMAT_INDEX(2) = 6  ---> INDEX_22
!      GREEKMAT_INDEX(3) = 2  ---> INDEX_12
!      GREEKMAT_INDEX(4) = 11 ---> INDEX_33
!      GREEKMAT_INDEX(5) = 12 ---> INDEX_34
!      GREEKMAT_INDEX(6) = 16 ---> INDEX_44

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

!  GEOMETRICAL QUANTITIES
!  ----------------------

      DO IB = 1, NBEAMS
        DO UM = 1, N_USER_STREAMS
          DO IA = 1, N_USER_RELAZMS
            V = VZA_OFFSETS(IB,UM) + IA

!  COSINE SCATTER ANGLE (THIS IS VALID ONLY FOR NON-REFRACTING ATMOSPHER
!  VSIGN = -1 FOR UPWELLING, +1 FOR DOWNWELLING

            COSSCAT = VSIGN * CTHETA(NV,IB) * CALPHA(UM) + &
                              STHETA(NV,IB) * SALPHA(UM) * CPHI(IA)
            UUU(V)  = COSSCAT

!  COSINE SIGMA 1 AND 2. H/VDM, EQS. (99)-(101)

!  A. SAFETY
!    WATCH FOR SIN^2(SCATTER ANGLE) LESS THAN ZERO (MACHINE PRECISION)
!    R. SPURR, 16 JANUARY 2006, RT SOLUTIONS INC.

            HELP_SINSCAT = ( ONE - COSSCAT * COSSCAT )
            IF ( HELP_SINSCAT.LE.ZERO ) THEN
              SINSCAT = 1.0D-12
            ELSE
              SINSCAT = DSQRT ( HELP_SINSCAT )
            ENDIF

!  B. NECESSARY LIMIT ANALYSES - HOVENIER LIMITS.
!     R. SPURR AND V. NATRAJ, 17 JANUARY 2006

            IF ( DABS(SINSCAT) .LE. 1.0D-12 ) THEN
              CSIG1 = ZERO
              CSIG2 = ZERO
            ELSE
             IF ( STHETA(1,IB) .EQ. ZERO ) THEN
              CSIG1 = -  CPHI(IA)
             ELSE
              CSIG1   = ( - VSIGN * CALPHA(UM) + CTHETA(NV,IB)*COSSCAT ) &
                        / SINSCAT / STHETA(NV,IB)
             ENDIF
             IF ( SALPHA(UM) .EQ. ZERO ) THEN
               CSIG2 = -  CPHI(IA)
             ELSE
               CSIG2   = ( - CTHETA(NV,IB) + VSIGN*CALPHA(UM)*COSSCAT ) &
                          / SINSCAT / SALPHA(UM)
             ENDIF
            ENDIF
!    WRITE(7,'(4I5,1P4E15.7)')V,IB,UM,IA,CSIG1,CSIG2,SINSCAT

!  THESE LINES ARE NECESSARY TO AVOID BAD VALUES

            IF ( CSIG1 .GT. ONE  ) CSIG1 = ONE
            IF ( CSIG1 .LT. -ONE ) CSIG1 = -ONE

            IF ( CSIG2 .GT. ONE  ) CSIG2 = ONE
            IF ( CSIG2 .LT. -ONE ) CSIG2 = -ONE

!  OUTPUT, H/VDM, EQS. (89)-(94)
!    ROTATION SINES AND COSINES

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

!  FOR RELAZM IN [180,360), NEED SIGN REVERSAL FOR S1 AND S2
!  SEE H/VDM, EQS. 94-95. V. NATRAJ AND R. SPURR, 01 MAY 2009.

            C1(V) = CSIG1_2 * CSIG1 - ONE
            C2(V) = CSIG2_2 * CSIG2 - ONE

            IF (PHI(IA) .LE. 180.D0) THEN
              S1(V) = CSIG1_2 * SSIG1
              S2(V) = CSIG2_2 * SSIG2
            ELSE
              S1(V) = -CSIG1_2 * SSIG1
              S2(V) = -CSIG2_2 * SSIG2
            ENDIF

!  END GEOMETRY LOOPS

          ENDDO
        ENDDO
      ENDDO

!  LINEARIZED F-MATRICES
!  ---------------------

      QROOT6 = -0.25D0 * DSQRT(6.0D0)

! INITIALISE F-MATRIX

      DO V = 1, N_GEOMETRIES
        DO K = 1, 6
          DO Q = 1, NV_PARAMETERS
            L_FMAT(Q,V,K) = ZERO
          END DO
        END DO
      END DO

!  START LOOP OVER THE COEFFICIENT INDEX L
!  FIRST UPDATE GENERALIZED SPHERICAL FUNCTIONS, THEN CALCULATE COEFS.
!  LOLD AND LNEW ARE POINTER-LIKE INDICES USED IN RECURRENCE

      LNEW = 1
      LOLD = 2

      DO L = 0, NGREEKMOMS

        DL   = DBLE(L)
        DL1  = DL - ONE

!  SET THE LOCAL GREEK MATRIX ELEMENTS THAT YOU NEED
!   44 AND 34 ARE NOT REQUIRED WITH NATURAL SUNLIGHT (DEFAULT HERE)
!   22 AND 33 REQUIRED FOR NON-MIE SPHEROIDAL PARTICLES

        DO Q = 1, NV_PARAMETERS

         IF ( DO_SCATMAT_VARIATION(NV,Q) ) THEN
          IF ( DO_SSCORR_TRUNCATION ) THEN

           DNL1  = DBLE(2*L + 1 )
           FDNL1 = SSFDEL(NV) * DNL1
           FACT  = ONE - SSFDEL(NV)

!  (1,1) AND (1,2) COMPONENTS

           GK11 = ( GREEKMAT(L,NV,1) - FDNL1 ) / FACT
           HELP1 = L_GREEKMAT(Q,L,NV,1) * GREEKMAT(L,NV,1)
           HELP2 = ( GK11 - DNL1 ) * L_SSFDEL(NV,Q)
           L_GK11(Q) = ( HELP1 + HELP2 ) / FACT
           IF ( NSTOKES.GT.1 ) THEN
             GK12 =   GREEKMAT(L,NV,2) / FACT
             HELP1 = L_GREEKMAT(Q,L,NV,2) * GREEKMAT(L,NV,2)
             HELP2 = GK12 * L_SSFDEL(NV,Q)
             L_GK12(Q) = ( HELP1 + HELP2 ) / FACT
           ENDIF

!  NON-SUNLIGHT CASES. (4,4) AND (3,4)

!mick fix
!           IF ( .NOT. DO_SUNLIGHT .AND. NSTOKES.EQ.4 ) THEN
!             GK44 = ( GREEKMAT(L,NV,16) - FDNL1 ) / FACT
!             GK34 =   GREEKMAT(L,NV,12) / FACT
!             HELP1 = L_GREEKMAT(Q,L,NV,16) * GREEKMAT(L,NV,16)
!             HELP2 = ( GK44 - DNL1 ) * L_SSFDEL(NV,Q)
!             L_GK44(Q) = ( HELP1 + HELP2 ) / FACT
!             HELP1 = L_GREEKMAT(Q,L,NV,12) * GREEKMAT(L,NV,12)
!             HELP2 = GK34 * L_SSFDEL(NV,Q)
!             L_GK34(Q) = ( HELP1 + HELP2 ) / FACT
!           ENDIF
           IF ( .NOT. DO_SUNLIGHT .AND. NSTOKES.EQ.4 ) THEN
             GK44 = ( GREEKMAT(L,NV,16) - FDNL1 ) / FACT
             GK34 =   GREEKMAT(L,NV,12) / FACT
             HELP1 = L_GREEKMAT(Q,L,NV,16) * GREEKMAT(L,NV,16)
             HELP2 = ( GK44 - DNL1 ) * L_SSFDEL(NV,Q)
             L_GK44(Q) = ( HELP1 + HELP2 ) / FACT
             HELP1 = L_GREEKMAT(Q,L,NV,12) * GREEKMAT(L,NV,12)
             HELP2 = GK34 * L_SSFDEL(NV,Q)
             L_GK34(Q) = ( HELP1 + HELP2 ) / FACT
           ELSE
             L_GK44(Q) = ZERO
             L_GK34(Q) = ZERO
           ENDIF

!  NON-SUNLIGHT CASES. (2,2) AND (3,3)

!mick fix
!           IF ( .NOT. DO_SUNLIGHT  ) THEN
!             GK22 = ( GREEKMAT(L,NV,6)  - FDNL1 ) / FACT
!             HELP1 = L_GREEKMAT(Q,L,NV,6) * GREEKMAT(L,NV,6)
!             HELP2 = ( GK22 - DNL1 ) * L_SSFDEL(NV,Q)
!             L_GK22(Q) = ( HELP1 + HELP2 ) / FACT
!             GK33 = ( GREEKMAT(L,NV,11) - FDNL1 ) / FACT
!             HELP1 = L_GREEKMAT(Q,L,NV,11) * GREEKMAT(L,NV,11)
!             HELP2 = ( GK33 - DNL1 ) * L_SSFDEL(NV,Q)
!             L_GK33(Q) = ( HELP1 + HELP2 ) / FACT
!           ENDIF
           IF ( .NOT. DO_SUNLIGHT  ) THEN
             GK22 = ( GREEKMAT(L,NV,6)  - FDNL1 ) / FACT
             HELP1 = L_GREEKMAT(Q,L,NV,6) * GREEKMAT(L,NV,6)
             HELP2 = ( GK22 - DNL1 ) * L_SSFDEL(NV,Q)
             L_GK22(Q) = ( HELP1 + HELP2 ) / FACT
             GK33 = ( GREEKMAT(L,NV,11) - FDNL1 ) / FACT
             HELP1 = L_GREEKMAT(Q,L,NV,11) * GREEKMAT(L,NV,11)
             HELP2 = ( GK33 - DNL1 ) * L_SSFDEL(NV,Q)
             L_GK33(Q) = ( HELP1 + HELP2 ) / FACT
           ELSE
             L_GK22(Q) = ZERO
             L_GK33(Q) = ZERO
           ENDIF

!  USUAL CASE WITH NO TRUNCATION

          ELSE

!  (1,1) AND (1,2) COMPONENTS

           L_GK11(Q) = L_GREEKMAT(Q,L,NV,1) * GREEKMAT(L,NV,1)

!mick fix
!           IF ( NSTOKES.GT.1 ) THEN
!             L_GK12(Q) = L_GREEKMAT(Q,L,NV,2) * GREEKMAT(L,NV,2)
!           ENDIF
           IF ( NSTOKES.GT.1 ) THEN
             L_GK12(Q) = L_GREEKMAT(Q,L,NV,2) * GREEKMAT(L,NV,2)
           ELSE
             L_GK12(Q) = ZERO
           ENDIF

!  NON-SUNLIGHT CASES. (4,4) AND (3,4)

!mick fix
!           IF ( .NOT. DO_SUNLIGHT .AND. NSTOKES.EQ.4 ) THEN
!             L_GK44(Q) = L_GREEKMAT(Q,L,NV,16) * GREEKMAT(L,NV,16)
!             L_GK34(Q) = L_GREEKMAT(Q,L,NV,12) * GREEKMAT(L,NV,12)
!           ENDIF
           IF ( .NOT. DO_SUNLIGHT .AND. NSTOKES.EQ.4 ) THEN
             L_GK44(Q) = L_GREEKMAT(Q,L,NV,16) * GREEKMAT(L,NV,16)
             L_GK34(Q) = L_GREEKMAT(Q,L,NV,12) * GREEKMAT(L,NV,12)
           ELSE
             L_GK44(Q) = ZERO
             L_GK34(Q) = ZERO
           ENDIF

!  NON-SUNLIGHT CASES. (2,2) AND (3,3)

!mick fix
!           IF ( .NOT. DO_SUNLIGHT  ) THEN
!             L_GK22(Q) = L_GREEKMAT(Q,L,NV,6)  * GREEKMAT(L,NV,6)
!             L_GK33(Q) = L_GREEKMAT(Q,L,NV,11) * GREEKMAT(L,NV,11)
!           ENDIF
           IF ( .NOT. DO_SUNLIGHT  ) THEN
             L_GK22(Q) = L_GREEKMAT(Q,L,NV,6)  * GREEKMAT(L,NV,6)
             L_GK33(Q) = L_GREEKMAT(Q,L,NV,11) * GREEKMAT(L,NV,11)
           ELSE
             L_GK22(Q) = ZERO
             L_GK33(Q) = ZERO
           ENDIF

          ENDIF
         ENDIF
        ENDDO

        IF ( L .EQ. 0 ) THEN

!  ADDING PAPER EQS. (76) AND (77) WITH M=0
!   ADDITIONAL FUNCTIONS P2M2 AND P2P2 ZERO FOR M = 0

          DO V = 1, N_GEOMETRIES
            P00(V,LOLD) = ONE
            P00(V,LNEW) = ZERO
            P02(V,LOLD) = ZERO
            P02(V,LNEW) = ZERO
            P2P2(V,LOLD) = ZERO
            P2P2(V,LNEW) = ZERO
            P2M2(V,LOLD) = ZERO
            P2M2(V,LNEW) = ZERO
          END DO

        ELSE

          FAC1 = (TWO*DL-ONE)/DL
          FAC2 = DL1/DL

! ADDING PAPER EQ. (81) WITH M=0

          DO V = 1, N_GEOMETRIES
            P00(V,LOLD) = FAC1*UUU(V)*P00(V,LNEW) - FAC2*P00(V,LOLD)
          END DO

        END IF

        IF ( L .EQ. 2 ) THEN

! ADDING PAPER EQ. (78)
! SQL4 CONTAINS THE FACTOR DSQRT((L+1)*(L+1)-4) NEEDED IN
! THE RECURRENCE EQS. (81) AND (82)

          DO V = 1, N_GEOMETRIES
            P02(V,LOLD) = QROOT6*(ONE-UUU(V)*UUU(V))
            P02(V,LNEW) = ZERO
          END DO
          SQL41 = ZERO

!  INTRODUCE THE P2P2 AND P2M2 FUNCTIONS FOR L = 2

          DO V = 1, N_GEOMETRIES
            P2P2(V,LOLD)= 0.25D0*(ONE+UUU(V))*(ONE+UUU(V))
            P2M2(V,LOLD)= 0.25D0*(ONE-UUU(V))*(ONE-UUU(V))
          ENDDO

        ELSE IF ( L .GT. 2) THEN

! ADDING PAPER EQ. (82) WITH M=0

          SQL4  = SQL41
          SQL41 = DSQRT(DL*DL-FOUR)
          TMP1  = (TWO*DL-ONE)/SQL41
          TMP2  = SQL4/SQL41
          DO V = 1, N_GEOMETRIES
            P02(V,LOLD) = TMP1*UUU(V)*P02(V,LNEW) - TMP2*P02(V,LOLD)
          END DO

!  INTRODUCE THE P2P2 AND P2M2 FUNCTIONS FOR L > 2

          FL2    = TWO * DL - ONE
          FLL1   = DL * DL1
          PERLL4 = ONE/(DL1*SQL41**2)
          QFAC   = DL  * ( DL1*DL1 - FOUR)
          DO V = 1, N_GEOMETRIES
           WFACT = FL2 * ( FLL1 * UUU(V) - FOUR )
           P2P2(V,LOLD) = (WFACT*P2P2(V,LNEW) - &
                QFAC*P2P2(V,LOLD)) * PERLL4
           WFACT = FL2 * ( FLL1 * UUU(V) + FOUR )
           P2M2(V,LOLD) = (WFACT*P2M2(V,LNEW) - &
                QFAC*P2M2(V,LOLD)) * PERLL4
          ENDDO

        END IF

! SWITCH INDICES SO THAT LNEW INDICATES THE FUNCTION WITH
! THE PRESENT INDEX VALUE L, THIS MECHANISM PREVENTS SWAPPING
! OF ENTIRE ARRAYS.

        ITMP = LNEW
        LNEW = LOLD
        LOLD = ITMP

! NOW ADD THE L-TH TERM TO THE SCATTERING MATRIX.
! SEE DE HAAN ET AL. (1987) EQS. (68)-(73).

! SECTION FOR RANDOMLY-ORIENTED SPHEROIDS, ADDED 20 MARCH 2006
!  R. SPURR AND V. NATRAJ

        DO Q = 1, NV_PARAMETERS

         IF ( DO_SCATMAT_VARIATION(NV,Q) ) THEN
          IF ( L.LE.LAYER_MAXMOMENTS(NV) ) THEN
           DO V = 1, N_GEOMETRIES
            L_FMAT(Q,V,INDEX_11) = &
                  L_FMAT(Q,V,INDEX_11) + L_GK11(Q) * P00(V,LNEW)
            L_FMAT(Q,V,INDEX_12) = &
                  L_FMAT(Q,V,INDEX_12) + L_GK12(Q) * P02(V,LNEW)

!  THIS SECTION ADDS F-MATRIX TERMS NEEDED FOR RANDOMLY-ORIENTED SPHEROIDS
            SUM23 = L_GK22(Q) + L_GK33(Q)
            DIF23 = L_GK22(Q) - L_GK33(Q)
            L_FMAT(Q,V,INDEX_22) = &
                  L_FMAT(Q,V,INDEX_22) + SUM23 * P2P2(V,LNEW)
            L_FMAT(Q,V,INDEX_33) = &
                  L_FMAT(Q,V,INDEX_33) + DIF23 * P2M2(V,LNEW)
!  END OF NEW SECTION---------------------------------------------------

            L_FMAT(Q,V,INDEX_44) = &
                  L_FMAT(Q,V,INDEX_44) + L_GK44(Q) * P00(V,LNEW)
            L_FMAT(Q,V,INDEX_34) = &
                  L_FMAT(Q,V,INDEX_34) + L_GK34(Q) * P02(V,LNEW)

           ENDDO
          ENDIF
         ENDIF
        END DO

!  END MOMENT LOOP

      END DO

!   THIS MUST BE DONE AFTER THE MOMENT LOOP.

      DO V = 1, N_GEOMETRIES
        DO Q = 1, NV_PARAMETERS
          L_FMAT(Q,V,INDEX_22) = HALF * &
             ( L_FMAT(Q,V,INDEX_22) + L_FMAT(Q,V,INDEX_33) )
          L_FMAT(Q,V,INDEX_33) = &
             ( L_FMAT(Q,V,INDEX_22) - L_FMAT(Q,V,INDEX_33) )
        END DO
      END DO

! REMEMBER FOR MIE SCATTERING : F11 = F22 AND F33 = F44
!  THIS CODE IS NO LONGER REQUIRED, AS WE HAVE INTRODUCED CODE NOW
!   FOR RANDOMLY ORIENTED SPHEROIDS. THE SYMMETRY SHOULD STILL OF
!   COURSE BE PRESENT FOR THE MIE PARTICLES, SO THIS WILL BE A
!   CHECK ON THE NEW CODE.
!   R. SPURR AND V. NATRAJ, 20 MARCH 2006

!      DO V = 1, N_GEOMETRIES
!        DO Q = 1, NV_PARAMETERS
!          L_FMAT(Q,V,INDEX_22) = L_FMAT(Q,V,INDEX_11)
!          L_FMAT(Q,V,INDEX_33) = L_FMAT(Q,V,INDEX_44)
!        END DO
!      END DO

!  CREATE LINEARIZED Z MATRICES
!  ----------------------------

!  COPY FOR STOKES = 1

      IF ( NSTOKES .EQ. 1 ) THEN
        DO Q = 1, NV_PARAMETERS
          DO V = 1, N_GEOMETRIES
            L_ZMAT(Q,V,NV,1,1) = L_FMAT(Q,V,INDEX_11)
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  FOR POLARIZED CASE, SUNLIGHT ONLY
!    ** ONLY NEED THE FIRST COLUMN OF THE Z-MATRIX

      IF ( DO_SUNLIGHT .AND. NSTOKES .GT. 1 ) THEN
        DO V = 1, N_GEOMETRIES
          DO Q = 1, NV_PARAMETERS
            L_ZMAT(Q,V,NV,1,1) =   L_FMAT(Q,V,INDEX_11)
            L_ZMAT(Q,V,NV,2,1) =  -L_FMAT(Q,V,INDEX_12) * C2(V)
            L_ZMAT(Q,V,NV,3,1) =   L_FMAT(Q,V,INDEX_12) * S2(V)
            L_ZMAT(Q,V,NV,4,1) =   ZERO
          ENDDO
        ENDDO
      ENDIF

!  FOR POLARIZED CASE, NON-SUNLIGHT GENERAL CASE

      IF ( .NOT. DO_SUNLIGHT .AND. NSTOKES .GT. 1 ) THEN
        DO V = 1, N_GEOMETRIES
          DO Q = 1, NV_PARAMETERS
            L_ZMAT(Q,V,NV,1,1) =   L_FMAT(Q,V,INDEX_11)
            L_ZMAT(Q,V,NV,2,1) =  -L_FMAT(Q,V,INDEX_12) * C2(V)
            L_ZMAT(Q,V,NV,3,1) =   L_FMAT(Q,V,INDEX_12) * S2(V)
            L_ZMAT(Q,V,NV,4,1) =   ZERO
!  THIS CODE UNTESTED
            HELP2C1 = L_FMAT(Q,V,INDEX_22) * C1(V)
            HELP2S1 = L_FMAT(Q,V,INDEX_22) * S1(V)
            HELP3C1 = L_FMAT(Q,V,INDEX_33) * C1(V)
            HELP3S1 = L_FMAT(Q,V,INDEX_33) * S1(V)
            L_ZMAT(Q,V,NV,1,2) =   L_FMAT(Q,V,INDEX_12) * C1(V)
            L_ZMAT(Q,V,NV,1,3) = - L_FMAT(Q,V,INDEX_12) * S1(V)
            L_ZMAT(Q,V,NV,1,4) =   ZERO
            L_ZMAT(Q,V,NV,2,2) =   C2(V) * HELP2C1 - S2(V) * HELP3S1
            L_ZMAT(Q,V,NV,2,3) = - C2(V) * HELP2S1 - S2(V) * HELP3C1
            L_ZMAT(Q,V,NV,2,4) = - L_FMAT(Q,V,INDEX_34) * S2(V)
            L_ZMAT(Q,V,NV,3,2) =   S2(V) * HELP2C1 + C2(V) * HELP3S1
            L_ZMAT(Q,V,NV,3,3) = - S2(V) * HELP2S1 + C2(V) * HELP3C1
            L_ZMAT(Q,V,NV,3,4) =   L_FMAT(Q,V,INDEX_34) * C2(V)
            L_ZMAT(Q,V,NV,4,2) = - L_FMAT(Q,V,INDEX_34) * S1(V)
            L_ZMAT(Q,V,NV,4,3) = - L_FMAT(Q,V,INDEX_34) * C1(V)
            L_ZMAT(Q,V,NV,4,4) =   L_FMAT(Q,V,INDEX_44)
          ENDDO
        ENDDO
      ENDIF

!  FINISH

      RETURN
      END SUBROUTINE VLIDORTSS_L_ZMATRICES

!

      SUBROUTINE VLIDORT_L_SSCORR_OUTGOING ( &
        SSFLUX, DO_SSCORR_TRUNCATION, DO_DELTAM_SCALING, &
        DO_UPWELLING, DO_DNWELLING, NSTOKES, NLAYERS, &
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
        DO_PROFILE_LINEARIZATION, DO_COLUMN_LINEARIZATION, &
        DO_ATMOS_LINEARIZATION, N_TOTALCOLUMN_WFS, &
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
        L_OMEGA_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT, &
        L_DELTAU_VERT, DO_SCATMAT_VARIATION, &
        UP_MULTIPLIERS, DN_MULTIPLIERS, UP_LOSTRANS, DN_LOSTRANS, &
        UP_MULTIPLIERS_UT, DN_MULTIPLIERS_UT, &
        UP_LOSTRANS_UT, DN_LOSTRANS_UT, &
        ZMAT_UP, ZMAT_DN, TMS, SSFDEL, &
        SS_CUMSOURCE_UP, SS_CUMSOURCE_DN, &
        BOA_ATTN, STOKES_SS, SS_CONTRIBS, &
        L_UP_MULTIPLIERS, L_DN_MULTIPLIERS, &
        L_UP_LOSTRANS, L_DN_LOSTRANS, &
        L_UP_MULTIPLIERS_UT, L_DN_MULTIPLIERS_UT, &
        L_UP_LOSTRANS_UT, L_DN_LOSTRANS_UT, &
        L_BOA_ATTN, &
        PROFILEWF_SS, COLUMNWF_SS, &
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
      USE VLIDORT_CORRECTIONS

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT (IN) :: SSFLUX
      LOGICAL, INTENT (IN) ::          DO_SSCORR_TRUNCATION
      LOGICAL, INTENT (IN) ::          DO_DELTAM_SCALING
      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_DNWELLING
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
      LOGICAL, INTENT (IN) ::          DO_PROFILE_LINEARIZATION
      LOGICAL, INTENT (IN) ::          DO_COLUMN_LINEARIZATION
      LOGICAL, INTENT (IN) ::          DO_ATMOS_LINEARIZATION
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
      DOUBLE PRECISION, INTENT (OUT) ::  L_UP_MULTIPLIERS &
         ( MAXLAYERS, 0:MAXLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (OUT) ::  L_DN_MULTIPLIERS &
         ( MAXLAYERS, 0:MAXLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (OUT) ::  L_UP_LOSTRANS &
         ( MAXLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (OUT) ::  L_DN_LOSTRANS &
         ( MAXLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (OUT) ::  L_UP_MULTIPLIERS_UT &
          ( MAX_PARTLAYERS, 0:MAXLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (OUT) ::  L_DN_MULTIPLIERS_UT &
          ( MAX_PARTLAYERS, 0:MAXLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (OUT) ::  L_UP_LOSTRANS_UT &
          ( MAX_PARTLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (OUT) ::  L_DN_LOSTRANS_UT &
          ( MAX_PARTLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )

      DOUBLE PRECISION, INTENT (OUT) ::  L_BOA_ATTN &
          ( 0:MAXLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )

      DOUBLE PRECISION, INTENT (INOUT) ::  PROFILEWF_SS &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
            MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
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
      INTEGER ::          KINDEX(MAXLAYERS)

!  LOCAL VARIABLES
!  ---------------

!  CUMULATIVE TRANSMITTANCES
!    NEW, 27 JANUARY 2010. TOA CONTRIBUTION FUNCTIONS CODE.

      DOUBLE PRECISION :: OGCTRANS(MAXLAYERS,MAX_GEOMETRIES)

!  INDICES

      INTEGER ::          N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER ::          UT, UTA, UM, IA, NC, IB, V, O1, O2, NM1
      INTEGER ::          K, K_PARAMETERS, KS, Q

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

      IF ( DO_ATMOS_LINEARIZATION ) THEN

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

        IF ( DO_ATMOS_LINEARIZATION ) THEN
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

      IF ( DO_PROFILE_LINEARIZATION ) THEN
        NTYPE_VARY = NLAYERS
        DO N = 1, NLAYERS
         DO_RTSOL_VARY(N) = LAYER_VARY_FLAG(N)
         NPARAMS_VARY(N)  = LAYER_VARY_NUMBER(N)
         KINDEX(N) = N
        ENDDO
      ELSE IF ( DO_COLUMN_LINEARIZATION ) THEN
        NTYPE_VARY = 1
        DO_RTSOL_VARY(1) = .TRUE.
        NPARAMS_VARY(1)  = N_TOTALCOLUMN_WFS
        KINDEX(1) = 0
      ENDIF

!  CREATE EXTINCTIONS
!  ------------------

!  USE BASIC DEFINITIONS

      DO N = 1, NLAYERS
        HELP = HEIGHT_GRID(N-1) - HEIGHT_GRID(N)
        EXTINCTION(N) = DELTAU_VERT(N) / HELP
      ENDDO

!  LINEARIZED EXTINCTIONS

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

              CALL L_OUTGOING_INTEGRATION_UP ( &
                NLAYERS, NFINELAYERS, &
                DO_PARTIALS, &
                DO_PROFILE_LINEARIZATION, LAYER_VARY_FLAG, &
                LAYER_VARY_NUMBER, &
                DO_COLUMN_LINEARIZATION, N_TOTALCOLUMN_WFS, &
                EXTINCTION, L_EXTINCTION, &
                N_PARTLAYERS, PARTLAYERS_LAYERIDX, &
                PARTLAYERS_LAYERFINEIDX, &
                SUNPATHS, RADII, NTRAVERSE, ALPHA_ALL, &
                SUNPATHS_UT, RADII_UT,  NTRAVERSE_UT, ALPHA_UT, &
                SUNPATHS_FINE, NTRAVERSE_FINE, ALPHA_FINE, &
                UP_MULTIPLIERS(1,V), UP_LOSTRANS(1,V), BOA_ATTN(V), &
                UP_MULTIPLIERS_UT(1,V), UP_LOSTRANS_UT(1,V), &
                L_UP_MULTIPLIERS(1,0,1,V), &
                L_UP_LOSTRANS(1,1,V), &
                L_BOA_ATTN(0,1,V), &
                L_UP_MULTIPLIERS_UT(1,0,1,V), &
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

!              IF ( DO_ATMOS_LINEARIZATION ) THEN
!              DO N = 1, NLAYERS
!                WRITE(36,'(2I4,1P20E18.10)')V,N,
!     &       UP_LOSTRANS(N,V),UP_MULTIPLIERS(N,V),
!     &L_UP_LOSTRANS(N,1,V),(L_UP_MULTIPLIERS(N,K,1,V),K=1,NLAYERS)
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
                  DO_SSCORR_TRUNCATION, SUNLIGHT, DO_ATMOS_LINEARIZATION, &
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
                 IF ( DO_ATMOS_LINEARIZATION ) THEN
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

              CALL L_OUTGOING_INTEGRATION_DN ( &
                NLAYERS, NFINELAYERS, DO_PARTIALS, &
                DO_PROFILE_LINEARIZATION, LAYER_VARY_FLAG, &
                LAYER_VARY_NUMBER, &
                DO_COLUMN_LINEARIZATION, N_TOTALCOLUMN_WFS, &
                EXTINCTION, L_EXTINCTION, &
                N_PARTLAYERS, PARTLAYERS_LAYERIDX, &
                PARTLAYERS_LAYERFINEIDX, &
                SUNPATHS, RADII, NTRAVERSE, ALPHA_ALL, &
                SUNPATHS_UT,   NTRAVERSE_UT,   ALPHA_UT, &
                SUNPATHS_FINE, NTRAVERSE_FINE, ALPHA_FINE, &
                DN_MULTIPLIERS(1,V), DN_LOSTRANS(1,V), &
                DN_MULTIPLIERS_UT(1,V), DN_LOSTRANS_UT(1,V), &
                L_DN_MULTIPLIERS(1,0,1,V), &
                L_DN_LOSTRANS(1,1,V), &
                L_DN_MULTIPLIERS_UT(1,0,1,V), &
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
                  DO_SSCORR_TRUNCATION, SUNLIGHT, DO_ATMOS_LINEARIZATION, &
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
                 IF ( DO_ATMOS_LINEARIZATION ) THEN
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

!                IF ( N.GT.23) WRITE(*,*)V,N,L_ZMAT_DN(1,V,N,1,1)

!  FINISH THE LAYER LOOP

              ENDDO

!  END DOWNWELLING CLAUSE

            ENDIF

!   FINISH GEOMETRY LOOPS

          ENDDO
        ENDDO
      ENDDO

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

      IF ( DO_UPWELLING .AND. DO_ATMOS_LINEARIZATION ) THEN

!  START THE MAIN LAYER VARIATION LOOP
!  -----------------------------------

        DO K = 1, NTYPE_VARY
         IF ( DO_RTSOL_VARY(K) ) THEN
          K_PARAMETERS = NPARAMS_VARY(K)
          KS = KINDEX(K)

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

!  IF N = K (PROFILE) OR KS = 0 (COLUMN), THERE ARE EXTRA LINEARIZATIONS
!  IF N > K, N < K, TRANSMITTANCE + SOME MULTIPLICATION LINEARIZATIONS

             IF ( N .EQ. K .OR. KS.EQ.0 ) THEN
              DO V = 1, N_GEOMETRIES
               DO Q = 1, K_PARAMETERS
                DO O1 = 1, NPOLAR
                 HELP = ZERO
                 DO O2 = 1, NSTOKES
                   LSS = ZMAT_UP(V,N,O1,O2) * L_UP_MULTIPLIERS(N,KS,Q,V) &
                     + L_ZMAT_UP(Q,V,N,O1,O2) * UP_MULTIPLIERS(N,V)
                   HELP = HELP + FLUXVEC(O2) * LSS
                 ENDDO
                 L_SS_CUMSOURCE(Q,V,O1) = HELP &
                  +  UP_LOSTRANS(N,V)    * L_SS_CUMSOURCE(Q,V,O1) &
                  + L_UP_LOSTRANS(N,Q,V) *   SS_CUMSOURCE_UP(V,O1,NC-1)
                ENDDO
               ENDDO
              ENDDO
             ELSE
              DO V = 1, N_GEOMETRIES
               DO Q = 1, K_PARAMETERS
                DO O1 = 1, NPOLAR
                 HELP = ZERO
                 DO O2 = 1, NSTOKES
                   LSS = ZMAT_UP(V,N,O1,O2) * L_UP_MULTIPLIERS(N,KS,Q,V)
                   HELP = HELP + FLUXVEC(O2) * LSS
                 ENDDO
                 L_SS_CUMSOURCE(Q,V,O1) = HELP &
                     +  UP_LOSTRANS(N,V) * L_SS_CUMSOURCE(Q,V,O1)
                ENDDO
               ENDDO
              ENDDO
             ENDIF

!  END LAYER LOOP

            ENDDO

!  OFFGRID OUTPUT----------
!  SET FINAL CUMULATIVE SOURCE AND SINGLE SCATTER WEIGHTING FUNCTION

           IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)

!  IF N = K (THE LAYER THAT IS VARYING)
!    ADD LINEARIZATION OF ADDITIONAL PARTIAL LAYER SOURCE TERM =
!        L_EXACT_SCAT(N) * MULTIPLIER  +  EXACT_SCAT * L_MULTIPLIER(N)
!  VARIATIONS WHEN N > K (VERY OCCASIONALLY N < K FOR TANGENT VIEWING)
!    ADD LINEARIZATION OF ADDITIONAL PARTIAL LAYER SOURCE TERM =
!         EXACT_SCAT * L_MULTIPLIER(K)

!  IF N = K (PROFILE) OR KS = 0 (COLUMN), THERE ARE EXTRA LINEARIZATIONS
!  IF N > K, N < K, TRANSMITTANCE + SOME MULTIPLICATION LINEARIZATIONS

            IF ( N.EQ.K .OR. KS.EQ.0 ) THEN
              DO V = 1, N_GEOMETRIES
               DO Q = 1, K_PARAMETERS
                LMULT = L_UP_MULTIPLIERS_UT(UT,KS,Q,V)
                DO O1 = 1, NPOLAR
                 HELP = ZERO
                 DO O2 = 1, NSTOKES
                  LSS = ZMAT_UP(V,N,O1,O2)   * LMULT &
                    + L_ZMAT_UP(Q,V,N,O1,O2) * UP_MULTIPLIERS_UT(UT,V)
                  HELP = HELP + FLUXVEC(O2) * LSS
                 ENDDO
                 L_FINAL_SOURCE = HELP &
                +  UP_LOSTRANS_UT(UT,V)    * L_SS_CUMSOURCE(Q,V,O1) &
               + L_UP_LOSTRANS_UT(UT,Q,V) *   SS_CUMSOURCE_UP(V,O1,NC)
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
               DO Q = 1, K_PARAMETERS
                LMULT = L_UP_MULTIPLIERS_UT(UT,KS,Q,V)
                DO O1 = 1, NPOLAR
                 HELP = ZERO
                 DO O2 = 1, NSTOKES
                  HELP = HELP + FLUXVEC(O2) * ZMAT_UP(V,N,O1,O2)
                 ENDDO
                 L_FINAL_SOURCE = HELP * LMULT &
                     +  UP_LOSTRANS_UT(UT,V) * L_SS_CUMSOURCE(Q,V,O1)
                 L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                 PROFILEWF_SS(Q,K,UTA,V,O1,UPIDX) = L_SSCORRECTION
                ENDDO
               ENDDO
              ENDDO
             ENDIF

!  ONGRID OUTPUT---------
!  JUST SET TO THE CUMULATIVE SOURCE TERM

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

      IF ( DO_DNWELLING .AND. DO_ATMOS_LINEARIZATION ) THEN

!  START THE MAIN LAYER VARIATION LOOP
!  -----------------------------------

        DO K = 1, NTYPE_VARY
         IF ( DO_RTSOL_VARY(K) ) THEN
          K_PARAMETERS = NPARAMS_VARY(K)
          KS = KINDEX(K)

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

!  IF N = K (PROFILE) OR KS = 0 (COLUMN), THERE ARE EXTRA LINEARIZATIONS
!  IF N > K, N < K, TRANSMITTANCE + SOME MULTIPLICATION LINEARIZATIONS

             IF ( N .EQ. K .OR. KS.EQ.0 ) THEN
              DO V = 1, N_GEOMETRIES
               DO Q = 1, K_PARAMETERS
                DO O1 = 1, NPOLAR
                 HELP = ZERO
                 DO O2 = 1, NSTOKES
                   LSS = ZMAT_DN(V,N,O1,O2) * L_DN_MULTIPLIERS(N,KS,Q,V) &
                     + L_ZMAT_DN(Q,V,N,O1,O2) * DN_MULTIPLIERS(N,V)
                   HELP = HELP + FLUXVEC(O2) * LSS
                 ENDDO
                 L_SS_CUMSOURCE(Q,V,O1) = HELP &
                  +  DN_LOSTRANS(N,V)    * L_SS_CUMSOURCE(Q,V,O1) &
                  + L_DN_LOSTRANS(N,Q,V) *   SS_CUMSOURCE_DN(V,O1,NC-1)
                ENDDO
               ENDDO
              ENDDO
             ELSE
              DO V = 1, N_GEOMETRIES
               DO Q = 1, K_PARAMETERS
                DO O1 = 1, NPOLAR
                 HELP = ZERO
                 DO O2 = 1, NSTOKES
                   LSS = ZMAT_DN(V,N,O1,O2) * L_DN_MULTIPLIERS(N,KS,Q,V)
                   HELP = HELP + FLUXVEC(O2) * LSS
                 ENDDO
                 L_SS_CUMSOURCE(Q,V,O1) = HELP &
                  +  DN_LOSTRANS(N,V) * L_SS_CUMSOURCE(Q,V,O1)
                ENDDO
               ENDDO
              ENDDO
             ENDIF

!  END LAYER LOOP

            ENDDO

!  OFFGRID OUTPUT----------
!  SET FINAL CUMULATIVE SOURCE AND SINGLE SCATTER WEIGHTING FUNCTION

           IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)

!  IF N = K (THE LAYER THAT IS VARYING)
!    ADD LINEARIZATION OF ADDITIONAL PARTIAL LAYER SOURCE TERM =
!        L_EXACT_SCAT(N) * MULTIPLIER  +  EXACT_SCAT * L_MULTIPLIER(N)
!  VARIATIONS WHEN N > K
!    ADD LINEARIZATION OF ADDITIONAL PARTIAL LAYER SOURCE TERM =
!         EXACT_SCAT * L_MULTIPLIER(K)

!  IF N = K (PROFILE) OR KS = 0 (COLUMN), THERE ARE EXTRA LINEARIZATIONS
!  IF N > K, N < K, TRANSMITTANCE + SOME MULTIPLICATION LINEARIZATIONS

            IF ( N.EQ.K .OR. KS.EQ.0 ) THEN
              DO V = 1, N_GEOMETRIES
               DO Q = 1, K_PARAMETERS
                LMULT = L_DN_MULTIPLIERS_UT(UT,KS,Q,V)
                DO O1 = 1, NPOLAR
                 HELP = ZERO
                 DO O2 = 1, NSTOKES
                  LSS = ZMAT_DN(V,N,O1,O2)   * LMULT &
                    + L_ZMAT_DN(Q,V,N,O1,O2) * DN_MULTIPLIERS_UT(UT,V)
                  HELP = HELP + FLUXVEC(O2) * LSS
                 ENDDO
                 L_FINAL_SOURCE = HELP &
                +  DN_LOSTRANS_UT(UT,V)    * L_SS_CUMSOURCE(Q,V,O1) &
               + L_DN_LOSTRANS_UT(UT,Q,V) *   SS_CUMSOURCE_DN(V,O1,NC)
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
               DO Q = 1, K_PARAMETERS
                LMULT = L_DN_MULTIPLIERS_UT(UT,KS,Q,V)
                DO O1 = 1, NPOLAR
                 HELP = ZERO
                 DO O2 = 1, NSTOKES
                  HELP = HELP + FLUXVEC(O2) * ZMAT_DN(V,N,O1,O2)
                 ENDDO
                 L_FINAL_SOURCE = HELP * LMULT &
                     +  DN_LOSTRANS_UT(UT,V) * L_SS_CUMSOURCE(Q,V,O1)
                 L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                 PROFILEWF_SS(Q,K,UTA,V,O1,DNIDX) = L_SSCORRECTION
                ENDDO
               ENDDO
              ENDDO
            ENDIF

!  ONGRID OUTPUT---------
!  JUST SET TO THE CUMULATIVE SOURCE TERM

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
      END SUBROUTINE VLIDORT_L_SSCORR_OUTGOING

!

      SUBROUTINE SSCORR_OUTGOING_L_ZMATRIX ( &
        DO_SSCORR_TRUNCATION, DO_SUNLIGHT, DO_LINEAR, &
        LAYER, NSTOKES, NGREEKMOMS, &
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
        LAYER_MAXMOMENTS, GREEKMATRIX, SSFDEL, &
        DO_SCATMAT, L_GREEKMATRIX, L_SSFDEL, &
        CTHETA, STHETA, CALPHA, SALPHA, CPHI, &
        PHI, COSSCAT, VSIGN, &
        ZMAT, L_ZMAT, FMAT )

      USE VLIDORT_PARS

      IMPLICIT NONE

!  INPUT
!  -----

!  LINEARIZATION CONTROL

      LOGICAL, INTENT (IN) ::          DO_LINEAR
      LOGICAL, INTENT (IN) ::          DO_SCATMAT ( MAXLAYERS, MAX_ATMOSWFS )
      LOGICAL, INTENT (IN) ::          LAYER_VARY_FLAG
      INTEGER, INTENT (IN) ::          LAYER_VARY_NUMBER

!  CONTROL INTEGERS AND SUNLIGHT FLAG

      LOGICAL, INTENT (IN) ::          DO_SUNLIGHT
      INTEGER, INTENT (IN) ::          LAYER, NSTOKES, NGREEKMOMS

!  SCATTERING INPUT INFORMATION

      INTEGER, INTENT (IN) ::          LAYER_MAXMOMENTS
      DOUBLE PRECISION, INTENT (IN) :: GREEKMATRIX &
          ( 0:MAXMOMENTS_INPUT, MAXLAYERS, 16 )
      DOUBLE PRECISION, INTENT (IN) :: L_GREEKMATRIX &
          ( MAX_ATMOSWFS, 0:MAXMOMENTS_INPUT, MAXLAYERS, 16 )

!  TRUNCATION CONTROL

      LOGICAL, INTENT (IN) ::          DO_SSCORR_TRUNCATION
      DOUBLE PRECISION, INTENT (IN) :: SSFDEL   ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_SSFDEL ( MAXLAYERS,MAX_ATMOSWFS )

!  ZENITH ANGLE, AZIMUTH ANGLE, SCATTER ANGLE (SINES/COSINES)

      DOUBLE PRECISION, INTENT (IN) :: CTHETA
      DOUBLE PRECISION, INTENT (IN) :: STHETA
      DOUBLE PRECISION, INTENT (IN) :: CALPHA
      DOUBLE PRECISION, INTENT (IN) :: SALPHA
      DOUBLE PRECISION, INTENT (IN) :: CPHI
      DOUBLE PRECISION, INTENT (IN) :: COSSCAT
      DOUBLE PRECISION, INTENT (IN) :: VSIGN

!  AZIMUTH ANGLE. ADDED FOR VERSION 2.4R, REQUIRED FOR PHI > 180.
!    R. SPURR AND V. NATRAJ, 01 MAY 2009

      DOUBLE PRECISION, INTENT (IN) :: PHI

!  OUTPUT
!  ------

!  Z-MATRIX, FIRST COLUMN ONLY. F-MATRIX ( SUNLIGHT CASE )
!  GENERALIZED, 05 OCTOBER 2010, NON-TESTED

      DOUBLE PRECISION, INTENT (OUT) :: ZMAT ( 4, 4 )
      DOUBLE PRECISION, INTENT (OUT) :: L_ZMAT ( 4, 4, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: FMAT ( 6 )

!  LOCAL
!  -----

      DOUBLE PRECISION :: C1, S1, C2, S2, SINSCAT, HELP_SINSCAT
      DOUBLE PRECISION :: CSIG1, CSIG2, SSIG1, SSIG2, CSIG1_2, CSIG2_2

      DOUBLE PRECISION :: DL, QROOT6, UUU, PHIDEG, DNL1, FDNL1, FACT
      DOUBLE PRECISION :: FAC1, FAC2, SQL4, SQL41
      DOUBLE PRECISION :: TMP1, TMP2, SUM23, DIF23, HELP1, HELP2
      DOUBLE PRECISION :: GK11, GK12, GK34, GK44, GK22, GK33
      DOUBLE PRECISION :: L_GK11(MAX_ATMOSWFS),L_GK12(MAX_ATMOSWFS)
      DOUBLE PRECISION :: L_GK44(MAX_ATMOSWFS),L_GK34(MAX_ATMOSWFS)
      DOUBLE PRECISION :: L_GK22(MAX_ATMOSWFS),L_GK33(MAX_ATMOSWFS)
      DOUBLE PRECISION :: FL2, FLL1, PERLL4, QFAC, WFACT, DL1
      DOUBLE PRECISION :: HELP2C1, HELP2S1, HELP3C1, HELP3S1

      DOUBLE PRECISION :: L_FMAT ( MAX_ATMOSWFS, 6 )

      INTEGER ::          GREEKMAT_INDEX(6), K, L, LNEW, LOLD, ITMP, N, Q
      INTEGER ::          INDEX_11, INDEX_12, INDEX_34
      INTEGER ::          INDEX_22, INDEX_33, INDEX_44

!  GENERALIZED SPHERICAL FUNCTIONS
!    P2P2 AND P2M2 ADDED, 05 OCTOBER 2010

      DOUBLE PRECISION :: P00(2), P02(2), P2P2(2), P2M2(2)

!  INDEXING KEY

!      GREEKMAT_INDEX(1) = 1  ---> INDEX_11
!      GREEKMAT_INDEX(2) = 6  ---> INDEX_22
!      GREEKMAT_INDEX(3) = 2  ---> INDEX_12
!      GREEKMAT_INDEX(4) = 11 ---> INDEX_33
!      GREEKMAT_INDEX(5) = 12 ---> INDEX_34
!      GREEKMAT_INDEX(6) = 16 ---> INDEX_44

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

!      IF ( V.EQ.42.AND.N.EQ.1.AND.VSIGN.GT.0.0D0 ) THEN
!         WRITE(*,*)'JJJ'
!         WRITE(*,*)CTHETA, STHETA, CALPHA, SALPHA, CPHI,
!     I    PHI, COSSCAT, VSIGN
!      ENDIF

!  GEOMETRICAL QUANTITIES
!  ----------------------

!  COSINE SCATTER ANGLE (THIS IS VALID ONLY FOR NON-REFRACTING ATMOSPHER
!  VSIGN = -1 FOR UPWELLING, +1 FOR DOWNWELLING
!  SHOULD BE THE SAME AS THE INPUT VALUE CSA

      UUU  = COSSCAT

!  COSINE SIGMA 1 AND 2. H/VDM, EQS. (99)-(101)

!  A. SAFETY
!    WATCH FOR SIN^2(SCATTER ANGLE) LESS THAN ZERO (MACHINE PRECISION)
!    R. SPURR, 16 JANUARY 2006, RT SOLUTIONS INC.

      HELP_SINSCAT = ( ONE - COSSCAT * COSSCAT )
      IF ( HELP_SINSCAT.LE.ZERO ) THEN
        SINSCAT = 1.0D-12
      ELSE
        SINSCAT = DSQRT ( HELP_SINSCAT )
      ENDIF

!  B. NECESSARY LIMIT ANALYSES - HOVENIER LIMITS.
!     R. SPURR AND V. NATRAJ, 17 JANUARY 2006

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

!  THESE LINES ARE NECESSARY TO AVOID BAD VALUES

      IF ( CSIG1 .GT. ONE  ) CSIG1 = ONE
      IF ( CSIG1 .LT. -ONE ) CSIG1 = -ONE

      IF ( CSIG2 .GT. ONE  ) CSIG2 = ONE
      IF ( CSIG2 .LT. -ONE ) CSIG2 = -ONE

!  OUTPUT, H/VDM, EQS. (89)-(94)

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

!  FOR RELAZM IN [180,360), NEED SIGN REVERSAL FOR S1 AND S2
!  SEE H/VDM, EQS. 94-95

      C1 = CSIG1_2 * CSIG1 - ONE
      C2 = CSIG2_2 * CSIG2 - ONE

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 !  BUG CORRECTED 01 OCTOBER 2010------------------------- RTS, R. SPURR

!  PHI IS INPUT IN RADIANS HERE. MUST CONVERT TO DEGREES
!     WE WERE USING DEGREES FOR THE S1/S2 SIGN REVERSLA

      PHIDEG = PHI / DEG_TO_RAD
!      IF ( PHI    .LE. 180.0D0 ) THEN            ! OLD
      IF ( PHIDEG .LE. 180.0D0 ) THEN             ! NEW
        S1 = CSIG1_2 * SSIG1
        S2 = CSIG2_2 * SSIG2
      ELSE
        S1 = -CSIG1_2 * SSIG1
        S2 = -CSIG2_2 * SSIG2
      ENDIF

 !  B G CORRECTED 01 OCTOBER 2010------------------------- RTS, R. SPURR
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  F-MATRICES AND LINEARIZED F-MATRICES
!  ------------------------------------

      QROOT6 = -0.25D0 * DSQRT(6.0D0)

! INITIALISE F-MATRIX AND ITS LINEARIZATION

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

!  START LOOP OVER THE COEFFICIENT INDEX L
!  FIRST UPDATE GENERALIZED SPHERICAL FUNCTIONS, THEN CALCULATE COEFS.
!  LOLD AND LNEW ARE POINTER-LIKE INDICES USED IN RECURRENCE

      LNEW = 1
      LOLD = 2

      DO L = 0, NGREEKMOMS

        DL   = DBLE(L)
        DL1  = DL - ONE

!  SET THE LOCAL GREEK MATRIX ELEMENTS THAT YOU NEED
!   44 AND 34 ARE NOT REQUIRED WITH NATURAL SUNLIGHT (DEFAULT HERE)
!   22 AND 33 REQUIRED FOR NON-MIE SPHEROIDAL PARTICLES

        IF ( DO_SSCORR_TRUNCATION ) THEN
          DNL1  = DBLE(2*L + 1 )
          FDNL1 = SSFDEL(N) * DNL1
          FACT  = ONE - SSFDEL(N)

!  REGULAR COEFFICIENTS

          GK11 = ( GREEKMATRIX(L,N,GREEKMAT_INDEX(1)) - FDNL1 ) / FACT
          GK12 =   GREEKMATRIX(L,N,GREEKMAT_INDEX(3)) / FACT
          IF ( .NOT. DO_SUNLIGHT ) THEN
           GK44 = ( GREEKMATRIX(L,N,GREEKMAT_INDEX(6)) - FDNL1 ) / FACT
           GK34 =   GREEKMATRIX(L,N,GREEKMAT_INDEX(5)) / FACT
           GK22 = ( GREEKMATRIX(L,N,GREEKMAT_INDEX(2)) - FDNL1 ) / FACT
           GK33 = ( GREEKMATRIX(L,N,GREEKMAT_INDEX(4)) - FDNL1 ) / FACT
          ENDIF

!  LINEARIZED COEFFICIENTS

          IF ( DO_LINEAR ) THEN
           IF ( LAYER_VARY_FLAG ) THEN
            DO Q = 1, LAYER_VARY_NUMBER
             IF ( DO_SCATMAT(N,Q) ) THEN

!  (1,1) AND (1,2) CASES

               HELP1 = L_GREEKMATRIX(Q,L,N,1)*GREEKMATRIX(L,N,1)
               HELP2 = ( GK11 - DNL1 ) * L_SSFDEL(N,Q)
               L_GK11(Q) = ( HELP1 + HELP2 ) / FACT
               HELP1 = L_GREEKMATRIX(Q,L,N,2)*GREEKMATRIX(L,N,2)
               HELP2 = GK12 * L_SSFDEL(N,Q)
               L_GK12(Q) = ( HELP1 + HELP2 ) / FACT

!  OTHER NON-SUNLIGHT CASES

               IF ( .NOT. DO_SUNLIGHT ) THEN

                 HELP1 = L_GREEKMATRIX(Q,L,N,16) * GREEKMATRIX(L,N,16)
                 HELP2 = ( GK44 - DNL1 ) * L_SSFDEL(N,Q)
                 L_GK44(Q) = ( HELP1 + HELP2 ) / FACT
                 HELP1 = L_GREEKMATRIX(Q,L,N,12) * GREEKMATRIX(L,N,12)
                 HELP2 = GK34 * L_SSFDEL(N,Q)
                 L_GK34(Q) = ( HELP1 + HELP2 ) / FACT

                 HELP1 = L_GREEKMATRIX(Q,L,N,6)  * GREEKMATRIX(L,N,6)
                 HELP2 = ( GK22 - DNL1 ) * L_SSFDEL(N,Q)
                 L_GK22(Q) = ( HELP1 + HELP2 ) / FACT
                 HELP1 = L_GREEKMATRIX(Q,L,N,11) * GREEKMATRIX(L,N,11)
                 HELP2 = ( GK33 - DNL1 ) * L_SSFDEL(N,Q)
                 L_GK33(Q) = ( HELP1 + HELP2 ) / FACT

               ENDIF

!  NO SCAT-MAT VARIATION

             ELSE
               L_GK11(Q) = ZERO
               L_GK12(Q) = ZERO
               IF ( .NOT. DO_SUNLIGHT ) THEN
                L_GK44(Q) = ZERO
                L_GK34(Q) = ZERO
                L_GK22(Q) = ZERO
                L_GK33(Q) = ZERO
               ENDIF
             ENDIF

!  END LINEARIZATION CONTROL

            ENDDO
           ENDIF
          ENDIF

!  NO SSCORR TRUNCATION, JUST COPY

        ELSE

!  REGULAR COEFFICIENTS

          GK11 = GREEKMATRIX(L,N,GREEKMAT_INDEX(1))
          GK12 = GREEKMATRIX(L,N,GREEKMAT_INDEX(3))
          IF ( .NOT. DO_SUNLIGHT ) THEN
           GK44 = GREEKMATRIX(L,N,GREEKMAT_INDEX(6))
           GK34 = GREEKMATRIX(L,N,GREEKMAT_INDEX(5))
           GK22 = GREEKMATRIX(L,N,GREEKMAT_INDEX(2))
           GK33 = GREEKMATRIX(L,N,GREEKMAT_INDEX(4))
          ENDIF

! LINAERIZED COEFFICIENTS

          IF ( DO_LINEAR ) THEN
           IF ( LAYER_VARY_FLAG ) THEN
            DO Q = 1, LAYER_VARY_NUMBER
             IF ( DO_SCATMAT(N,Q) ) THEN
               L_GK11(Q) = L_GREEKMATRIX(Q,L,N,1)*GREEKMATRIX(L,N,1)
               L_GK12(Q) = L_GREEKMATRIX(Q,L,N,2)*GREEKMATRIX(L,N,2)
               IF ( .NOT. DO_SUNLIGHT ) THEN
                L_GK44(Q) = L_GREEKMATRIX(Q,L,N,16)*GREEKMATRIX(L,N,16)
                L_GK34(Q) = L_GREEKMATRIX(Q,L,N,12)*GREEKMATRIX(L,N,12)
                L_GK22(Q) = L_GREEKMATRIX(Q,L,N,6) *GREEKMATRIX(L,N,6)
                L_GK33(Q) = L_GREEKMATRIX(Q,L,N,11)*GREEKMATRIX(L,N,11)
               ENDIF
             ELSE
               L_GK11(Q) = ZERO
               L_GK12(Q) = ZERO
               IF ( .NOT. DO_SUNLIGHT ) THEN
                L_GK44(Q) = ZERO
                L_GK34(Q) = ZERO
                L_GK22(Q) = ZERO
                L_GK33(Q) = ZERO
               ENDIF
             ENDIF
            ENDDO
           ENDIF
          ENDIF
        ENDIF

!  FIRST MOMENT

        IF ( L .EQ. 0 ) THEN

!  ADDING PAPER EQS. (76) AND (77) WITH M=0
!   ADDITIONAL FUNCTIONS P2M2 AND P2P2 ZERO FOR M = 0

          P00(LOLD) = ONE
          P00(LNEW) = ZERO
          P02(LOLD) = ZERO
          P02(LNEW) = ZERO
          P2P2(LOLD) = ZERO
          P2P2(LNEW) = ZERO
          P2M2(LOLD) = ZERO
          P2M2(LNEW) = ZERO

        ELSE

          FAC1 = (TWO*DL-ONE)/DL
          FAC2 = DL1/DL

! ADDING PAPER EQ. (81) WITH M=0

          P00(LOLD) = FAC1*UUU*P00(LNEW) - FAC2*P00(LOLD)

        END IF

        IF ( L .EQ. 2 ) THEN

! ADDING PAPER EQ. (78)
! SQL4 CONTAINS THE FACTOR DSQRT((L+1)*(L+1)-4) NEEDED IN
! THE RECURRENCE EQS. (81) AND (82)

          P02(LOLD) = QROOT6*(ONE-UUU*UUU)
          P02(LNEW) = ZERO
          SQL41 = ZERO

!  INTRODUCE THE P2P2 AND P2M2 FUNCTIONS FOR L = 2

          P2P2(LOLD)= 0.25D0*(ONE+UUU)*(ONE+UUU)
          P2M2(LOLD)= 0.25D0*(ONE-UUU)*(ONE-UUU)

        ELSE IF ( L .GT. 2) THEN

! ADDING PAPER EQ. (82) WITH M=0

          SQL4  = SQL41
          SQL41 = DSQRT(DL*DL-4.0D0)
          TMP1  = (2.0D0*DL-1.0D0)/SQL41
          TMP2  = SQL4/SQL41
          P02(LOLD) = TMP1*UUU*P02(LNEW) - TMP2*P02(LOLD)

!  INTRODUCE THE P2P2 AND P2M2 FUNCTIONS FOR L > 2

          FL2 = TWO * DL - ONE
          FLL1 = DL * DL1
          PERLL4=ONE/(DL1*SQL41**2)
          QFAC  = DL  * ( DL1*DL1 - FOUR)
          WFACT = FL2 * ( FLL1 * UUU - FOUR )
          P2P2(LOLD) = (WFACT*P2P2(LNEW) - QFAC*P2P2(LOLD)) * PERLL4
          WFACT = FL2 * ( FLL1 * UUU + FOUR )
          P2M2(LOLD) = (WFACT*P2M2(LNEW) - QFAC*P2M2(LOLD)) * PERLL4

        END IF

! SWITCH INDICES SO THAT LNEW INDICATES THE FUNCTION WITH
! THE PRESENT INDEX VALUE L, THIS MECHANISM PREVENTS SWAPPING
! OF ENTIRE ARRAYS.

        ITMP = LNEW
        LNEW = LOLD
        LOLD = ITMP

! NOW ADD THE L-TH TERM TO THE SCATTERING MATRIX.
! SEE DE HAAN ET AL. (1987) EQS. (68)-(73).

! REMEMBER FOR MIE SCATTERING : F11 = F22 AND F33 = F44
!  DO NOT NEED THESE WITH NATURAL SUNLIGHT:
!            FMAT(INDEX_34)AND FMAT(INDEX_44)

! SECTION FOR RANDOMLY-ORIENTED SPHEROIDS, ADDED 05 OCTOBER 2010
!  R. SPURR AND V. NATRAJ

        IF ( L.LE.LAYER_MAXMOMENTS ) THEN

          FMAT(INDEX_11) = FMAT(INDEX_11) + GK11 * P00(LNEW)
          FMAT(INDEX_12) = FMAT(INDEX_12) + GK12 * P02(LNEW)

          IF ( .NOT. DO_SUNLIGHT ) THEN
            SUM23 = GK22 + GK33
            DIF23 = GK22 - GK33
            FMAT(INDEX_22) = FMAT(INDEX_22) + SUM23 * P2P2(LNEW)
            FMAT(INDEX_33) = FMAT(INDEX_33) + DIF23 * P2M2(LNEW)
            FMAT(INDEX_44) = FMAT(INDEX_44) + GK44 * P00(LNEW)
            FMAT(INDEX_34) = FMAT(INDEX_34) + GK34 * P02(LNEW)
          ENDIF

        ENDIF

!  LINEARIZATION

        IF ( DO_LINEAR ) THEN
         IF ( LAYER_VARY_FLAG ) THEN
          DO Q = 1, LAYER_VARY_NUMBER
           IF ( DO_SCATMAT(N,Q) ) THEN
            IF ( L.LE.LAYER_MAXMOMENTS ) THEN

             L_FMAT(Q,INDEX_11) = L_FMAT(Q,INDEX_11)+L_GK11(Q)*P00(LNEW)
             L_FMAT(Q,INDEX_12) = L_FMAT(Q,INDEX_12)+L_GK12(Q)*P02(LNEW)
             IF ( .NOT. DO_SUNLIGHT ) THEN
              SUM23 = L_GK22(Q) + L_GK33(Q)
              DIF23 = L_GK22(Q) - L_GK33(Q)
              L_FMAT(Q,INDEX_22) = L_FMAT(Q,INDEX_22)+SUM23*P2P2(LNEW)
              L_FMAT(Q,INDEX_33) = L_FMAT(Q,INDEX_33)+DIF23*P2M2(LNEW)
              L_FMAT(Q,INDEX_44) =L_FMAT(Q,INDEX_44)+L_GK44(Q)*P00(LNEW)
              L_FMAT(Q,INDEX_34) =L_FMAT(Q,INDEX_34)+L_GK34(Q)*P02(LNEW)
             ENDIF

            ENDIF
           ENDIF
          ENDDO
         ENDIF
        ENDIF

!  END MOMENTS

      END DO

!   THIS MUST BE DONE AFTER THE MOMENT LOOP. (NON-SUNLIGHT ONLY)

      IF ( .NOT. DO_SUNLIGHT ) THEN

        FMAT(INDEX_22) = HALF * ( FMAT(INDEX_22) + FMAT(INDEX_33) )
        FMAT(INDEX_33) = FMAT(INDEX_22) - FMAT(INDEX_33)

        IF ( DO_LINEAR ) THEN
          IF ( LAYER_VARY_FLAG ) THEN
            DO Q = 1, LAYER_VARY_NUMBER
              IF ( DO_SCATMAT(N,Q) ) THEN
                L_FMAT(Q,INDEX_22) = HALF * &
                   ( L_FMAT(Q,INDEX_22) + L_FMAT(Q,INDEX_33) )
                L_FMAT(Q,INDEX_33) = &
                   ( L_FMAT(Q,INDEX_22) - L_FMAT(Q,INDEX_33) )
              ENDIF
            ENDDO
          ENDIF
        ENDIF

      ENDIF

! REMEMBER FOR MIE SCATTERING : F11 = F22 AND F33 = F44
!  THIS CODE IS NO LONGER REQUIRED, AS WE HAVE INTRODUCED CODE NOW
!   FOR RANDOMLY ORIENTED SPHEROIDS. THE SYMMETRY SHOULD STILL OF
!   COURSE BE PRESENT FOR THE MIE PARTICLES, SO THIS WILL BE A
!   CHECK ON THE NEW CODE. R. SPURR AND V. NATRAJ,, 20 MARCH 2006

!      IF ( .NOT. DO_SUNLIGHT ) THEN
!        FMAT(INDEX_22) = FMAT(INDEX_11)
!        FMAT(INDEX_33) = FMAT(INDEX_44)
!        IF ( DO_LINEAR ) THEN
!          IF ( LAYER_VARY_FLAG ) THEN
!            DO Q = 1, LAYER_VARY_NUMBER
!              IF ( DO_SCATMAT(N,Q) ) THEN
!                L_FMAT(Q,INDEX_22) = L_FMAT(Q,INDEX_11)
!                L_FMAT(Q,INDEX_33) = L_FMAT(Q,INDEX_44)
!              ENDIF
!            ENDDO
!          ENDIF
!        ENDIF
!      ENDIF

!  Z-MATRIX CALCULATION
!  ====================

!  INITIALIZE Z-MATRIX + LINEARIZATION

      ZMAT   = ZERO
      L_ZMAT = ZERO

!  SCALAR CASE

      IF ( NSTOKES .EQ. 1 ) THEN
        ZMAT(1,1) =   FMAT(INDEX_11)
      ENDIF

!  NATURAL LIGHT, FIRST COLUMN ONLY

      IF ( DO_SUNLIGHT .AND. NSTOKES.GT.1 ) THEN
        ZMAT(1,1) =   FMAT(INDEX_11)
        ZMAT(2,1) =  -FMAT(INDEX_12) * C2
        ZMAT(3,1) =   FMAT(INDEX_12) * S2
        ZMAT(4,1) =   ZERO
      ENDIF

!  OTHER SOURCES, USE FULL 4X4 MATRIX
!    CODE INTRODUCED BUT NOT TESTED, 05 OCTOBER 2010

      IF ( .NOT. DO_SUNLIGHT .AND. NSTOKES.GT.1 ) THEN
         HELP2C1 = FMAT(INDEX_22) * C1
         HELP2S1 = FMAT(INDEX_22) * S1
         HELP3C1 = FMAT(INDEX_33) * C1
         HELP3S1 = FMAT(INDEX_33) * S1
         ZMAT(1,2) =   FMAT(INDEX_12) * C1
         ZMAT(1,3) = - FMAT(INDEX_12) * S1
         ZMAT(1,4) =   ZERO
         ZMAT(2,2) =   C2 * HELP2C1 - S2 * HELP3S1
         ZMAT(2,3) = - C2 * HELP2S1 - S2 * HELP3C1
         ZMAT(2,4) = - FMAT(INDEX_34) * S2
         ZMAT(3,2) =   S2 * HELP2C1 + C2 * HELP3S1
         ZMAT(3,3) = - S2 * HELP2S1 + C2 * HELP3C1
         ZMAT(3,4) =   FMAT(INDEX_34) * C2
         ZMAT(4,2) = - FMAT(INDEX_34) * S1
         ZMAT(4,3) = - FMAT(INDEX_34) * C1
         ZMAT(4,4) =   FMAT(INDEX_44)
      ENDIF

!  RETURN IF NO LINEARIZATION

      IF ( .NOT. DO_LINEAR ) RETURN

!  LINEARIZED Z-MATRIX
!  -------------------

!  COPY FOR STOKES = 1
!  FOR POLARIZED CASE, SUNLIGHT ONLY !!!!!!!
!   FOR SUNLIGHT, ONLY NEED THE FIRST COLUMN OF THE Z-MATRIX

      DO Q = 1, LAYER_VARY_NUMBER

        IF ( NSTOKES .EQ. 1 ) THEN
          L_ZMAT(1,1,Q) = L_FMAT(Q,INDEX_11)
        ENDIF

        IF ( DO_SUNLIGHT .AND. NSTOKES .GT. 1 ) THEN
          L_ZMAT(1,1,Q) =   L_FMAT(Q,INDEX_11)
          L_ZMAT(2,1,Q) =  -L_FMAT(Q,INDEX_12) * C2
          L_ZMAT(3,1,Q) =   L_FMAT(Q,INDEX_12) * S2
          L_ZMAT(4,1,Q) =   ZERO
        ENDIF

!  OTHER SOURCES, USE FULL 4X4 MATRIX
!    CODE INTRODUCED BUT NOT TESTED, 05 OCTOBER 2010

        IF ( .NOT. DO_SUNLIGHT .AND. NSTOKES .GT. 1 ) THEN
           HELP2C1 = L_FMAT(Q,INDEX_22) * C1
           HELP2S1 = L_FMAT(Q,INDEX_22) * S1
           HELP3C1 = L_FMAT(Q,INDEX_33) * C1
           HELP3S1 = L_FMAT(Q,INDEX_33) * S1
           L_ZMAT(1,2,Q) =   L_FMAT(Q,INDEX_12) * C1
           L_ZMAT(1,3,Q) = - L_FMAT(Q,INDEX_12) * S1
           L_ZMAT(1,4,Q) =   ZERO
           L_ZMAT(2,2,Q) =   C2 * HELP2C1 - S2 * HELP3S1
           L_ZMAT(2,3,Q) = - C2 * HELP2S1 - S2 * HELP3C1
           L_ZMAT(2,4,Q) = - L_FMAT(Q,INDEX_34) * S2
           L_ZMAT(3,2,Q) =   S2 * HELP2C1 + C2 * HELP3S1
           L_ZMAT(3,3,Q) = - S2 * HELP2S1 + C2 * HELP3C1
           L_ZMAT(3,4,Q) =   L_FMAT(Q,INDEX_34) * C2
           L_ZMAT(4,2,Q) = - L_FMAT(Q,INDEX_34) * S1
           L_ZMAT(4,3,Q) = - L_FMAT(Q,INDEX_34) * C1
           L_ZMAT(4,4,Q) =   L_FMAT(Q,INDEX_44)
        ENDIF

      ENDDO

!  FINISH

      RETURN
      END SUBROUTINE SSCORR_OUTGOING_L_ZMATRIX

!

      SUBROUTINE L_OUTGOING_INTEGRATION_UP ( &
        NLAYERS, NFINELAYERS, &
        DO_PARTIALS, &
        DO_PROFILE_LINEARIZATION, LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
        DO_COLUMN_LINEARIZATION, N_TOTALCOLUMN_WFS, &
        EXTINCTION, L_EXTINCTION, &
        N_PARTIALS, PARTIALS_IDX, PARTIALS_FINEIDX, &
        SUNPATHS, RADII, NTRAVERSE, ALPHA_ALL, &
        SUNPATHS_P, RADII_P, NTRAVERSE_P,    ALPHA_P, &
        SUNPATHS_FINE, NTRAVERSE_FINE, ALPHA_FINE, &
        MULTIPLIERS,     LOSTRANS,     BOA_ATTN, &
        MULTIPLIERS_P,   LOSTRANS_P, &
        L_MULTIPLIERS,   L_LOSTRANS,   L_BOA_ATTN, &
        L_MULTIPLIERS_P, L_LOSTRANS_P )

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

      LOGICAL, INTENT (IN) ::          DO_PROFILE_LINEARIZATION
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
      DOUBLE PRECISION, INTENT (OUT) :: L_MULTIPLIERS &
          (MAXLAYERS,0:MAXLAYERS,MAX_ATMOSWFS)
      DOUBLE PRECISION, INTENT (OUT) :: L_LOSTRANS   (MAXLAYERS,MAX_ATMOSWFS)

      DOUBLE PRECISION, INTENT (OUT) :: MULTIPLIERS_P(MAX_PARTLAYERS)
      DOUBLE PRECISION, INTENT (OUT) :: LOSTRANS_P   (MAX_PARTLAYERS)
      DOUBLE PRECISION, INTENT (OUT) :: L_MULTIPLIERS_P &
          (MAX_PARTLAYERS,0:MAXLAYERS,MAX_ATMOSWFS)
      DOUBLE PRECISION, INTENT (OUT) :: L_LOSTRANS_P &
          (MAX_PARTLAYERS,MAX_ATMOSWFS)

      DOUBLE PRECISION, INTENT (OUT) :: BOA_ATTN
      DOUBLE PRECISION, INTENT (OUT) :: L_BOA_ATTN (0:MAXLAYERS,MAX_ATMOSWFS)

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

      DOUBLE PRECISION :: L_ATTN ( 0:MAXLAYERS, 0:MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_ATTN_FINE ( MAXLAYERS, MAXFINELAYERS, &
                                     0:MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_ATTN_P &
          ( MAX_PARTLAYERS, 0:MAXLAYERS, MAX_ATMOSWFS  )

!  HELP VARIABLES
!  --------------

      INTEGER ::          N, J, K, K0, Q, UT, NP, NFINE_P

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

      IF ( DO_PROFILE_LINEARIZATION ) THEN
        DO N = 1, NLAYERS
          DO K = 1, NLAYERS
            IF ( LAYER_VARY_FLAG(K) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(K)
                L_MULTIPLIERS(N,K,Q) = ZERO
              ENDDO
            ENDIF
          ENDDO
          IF ( LAYER_VARY_FLAG(N) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_LOSTRANS(N,Q) = ZERO
            ENDDO
          ENDIF
        ENDDO
      ENDIF

      IF ( DO_COLUMN_LINEARIZATION ) THEN
        K0 = 0
        DO N = 1, NLAYERS
          DO Q = 1, N_TOTALCOLUMN_WFS
            L_MULTIPLIERS(N,K0,Q) = ZERO
            L_LOSTRANS(N,Q)       = ZERO
          ENDDO
        ENDDO
      ENDIF

!  PARTIAL LAYER LINEARIZATIONS

      IF ( DO_PARTIALS ) THEN
       IF ( DO_PROFILE_LINEARIZATION ) THEN
        DO UT = 1, N_PARTIALS
          N = PARTIALS_IDX(UT)
          DO K = 1, NLAYERS
            IF ( LAYER_VARY_FLAG(K) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(K)
                L_MULTIPLIERS_P(UT,K,Q) = ZERO
              ENDDO
            ENDIF
          ENDDO
          IF ( LAYER_VARY_FLAG(N) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_LOSTRANS_P(UT,Q) = ZERO
            ENDDO
          ENDIF
        ENDDO
       ENDIF
       IF ( DO_COLUMN_LINEARIZATION ) THEN
        K0 = 0
        DO UT = 1, N_PARTIALS
          DO Q = 1, N_TOTALCOLUMN_WFS
            L_MULTIPLIERS_P(UT,K0,Q) = ZERO
            L_LOSTRANS_P(UT,Q)       = ZERO
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
      IF ( DO_PROFILE_LINEARIZATION.OR.DO_COLUMN_LINEARIZATION ) THEN
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

!  LINEARIZED ATTENUATION FACTORS
!     NOTE THE EXTRA ZEROING...............BUG, 01 NOVEMBER 2007

      IF ( DO_PROFILE_LINEARIZATION ) THEN
        DO N = 0, NLAYERS
          DO K = N+1, NLAYERS
            DO Q = 1, LAYER_VARY_NUMBER(K)
              L_ATTN(N,K,Q) = ZERO
            ENDDO
          ENDDO
          DO K = 1, NTRAVERSE(N)
            IF ( LAYER_VARY_FLAG(K) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(K)
                L_IOP = - L_EXTINCTION(K,Q)
                L_ATTN(N,K,Q) = SUNPATHS(N,K) * ATTN(N) * L_IOP
                IF (N.EQ.NLAYERS)L_BOA_ATTN(K,Q)=L_ATTN(N,K,Q)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
        DO N = 1, NLAYERS
         DO J = 1, NFINELAYERS
          DO K = 1, NTRAVERSE_FINE(N,J)
            IF ( LAYER_VARY_FLAG(K) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(K)
                L_IOP = - L_EXTINCTION(K,Q) * ATTN_FINE(N,J)
                L_ATTN_FINE(N,J,K,Q) = SUNPATHS_FINE(N,K,J) * L_IOP
              ENDDO
            ENDIF
          ENDDO
         ENDDO
        ENDDO
        IF ( DO_PARTIALS ) THEN
         DO UT = 1, N_PARTIALS
          DO K = 1, NTRAVERSE_P(UT)
            IF ( LAYER_VARY_FLAG(K) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(K)
                L_IOP = - L_EXTINCTION(K,Q) * ATTN_P(UT)
                L_ATTN_P(UT,K,Q) = SUNPATHS_P(UT,K) * L_IOP
              ENDDO
            ENDIF
          ENDDO
         ENDDO
        ENDIF
      ENDIF

!  LINEARIZED ATTENUATION FACTORS. COLUMN LINEARIZATION

      IF ( DO_COLUMN_LINEARIZATION ) THEN
        DO N = 0, NLAYERS
          DO Q = 1, N_TOTALCOLUMN_WFS
            SUM = ZERO
            DO K = 1, NTRAVERSE(N)
              L_IOP = - L_EXTINCTION(K,Q)
              SUM = SUM + SUNPATHS(N,K) * L_IOP
            ENDDO
            L_ATTN(N,K0,Q) = SUM * ATTN(N)
            IF (N.EQ.NLAYERS)L_BOA_ATTN(K0,Q)=L_ATTN(N,K0,Q)
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
            L_ATTN_FINE(N,J,K0,Q) = SUM * ATTN_FINE(N,J)
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
            L_ATTN_P(UT,K0,Q) = SUM * ATTN_P(UT)
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

!  PROFILE LINEARIZATIONS

          IF ( DO_PROFILE_LINEARIZATION ) THEN
           DO K = 1, NLAYERS
            IF ( LAYER_VARY_FLAG(K) ) THEN
             DO Q = 1, LAYER_VARY_NUMBER(K)
              IF ( K.EQ.N) THEN
               L_KN = L_EXTINCTION(N,Q)
               L_FUNC_2 = L_ATTN(N-1,N,Q)
               L_TRAN_1 = - L_KN * HN * TRAN_1
               L_FUNC_1 = L_ATTN(N,N,Q) * TRAN_1 + ATTN(N) * L_TRAN_1
               L_SUM = 0.5D0 * ( L_FUNC_1 + L_FUNC_2 )
               DO J = 1, NFINELAYERS
                L_TRAN = - L_KN * HJ(J) * TRAN(J)
                L_FUNC = L_ATTN_FINE(N,J,N,Q) * TRAN(J) + &
                           ATTN_FINE(N,J) * L_TRAN
                L_SUM = L_SUM + L_FUNC
               ENDDO
               L_LOSTRANS(N,Q)      = L_TRAN_1
               L_MULTIPLIERS(N,N,Q) = STEP * ( L_KN * SUM + KN * L_SUM )
              ELSE
               L_FUNC_2 = L_ATTN(N-1,K,Q)
               L_FUNC_1 = L_ATTN(N,K,Q) * TRAN_1
               L_SUM = 0.5D0 * ( L_FUNC_1 + L_FUNC_2 )
               DO J = 1, NFINELAYERS
                L_FUNC = L_ATTN_FINE(N,J,K,Q) * TRAN(J)
                L_SUM = L_SUM + L_FUNC
               ENDDO
               L_MULTIPLIERS(N,K,Q) = STEP * KN * L_SUM
              ENDIF
             ENDDO
            ENDIF
           ENDDO

!  DEBUG SPECIAL CASE
!           WRITE(65,'(I4,1P3E19.9)')N,LOSTRANS(N),MULTIPLIERS(N),
!     &                              L_LOSTRANS(N,1)
!           DO K = 1, NLAYERS
!             WRITE(65,'(4X,I4,1PE19.9)')K,L_MULTIPLIERS(N,K,1)
!           ENDDO

          ENDIF

!  COLUMN LINEARIZATION.

          IF ( DO_COLUMN_LINEARIZATION ) THEN
           DO Q = 1, N_TOTALCOLUMN_WFS
            L_KN = L_EXTINCTION(N,Q)
            L_FUNC_2 = L_ATTN(N-1,K0,Q)
            L_TRAN_1 = -L_KN * HN * TRAN_1
            L_FUNC_1 = L_ATTN(N,K0,Q) * TRAN_1 + ATTN(N) * L_TRAN_1
            L_SUM = 0.5D0 * ( L_FUNC_1 + L_FUNC_2 )
            DO J = 1, NFINELAYERS
              L_TRAN = - L_KN * HJ(J) * TRAN(J)
              L_FUNC = L_ATTN_FINE(N,J,K0,Q) * TRAN(J) &
                     +   ATTN_FINE(N,J)      * L_TRAN
              L_SUM = L_SUM + L_FUNC
            ENDDO
            L_LOSTRANS(N,Q)      = L_TRAN_1
            L_MULTIPLIERS(N,K0,Q) = STEP * ( L_KN * SUM + KN * L_SUM )
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

!  PROFILE LINEARIZATION.

          IF ( DO_PROFILE_LINEARIZATION ) THEN

           DO K = 1, NLAYERS
            IF ( LAYER_VARY_FLAG(K) ) THEN
             DO Q = 1, LAYER_VARY_NUMBER(K)
              IF ( K.EQ.NP) THEN
               L_KN = L_EXTINCTION(NP,Q)
               L_FUNC_2 = L_ATTN_P(UT,NP,Q)
               L_TRAN_1 = -L_KN * HP * TRAN_1
               L_FUNC_1 = L_ATTN(NP,NP,Q) * TRAN_1 + ATTN(NP) * L_TRAN_1
               IF ( NFINE_P.GT.0) THEN
                L_SUM = 0.5D0 * L_FUNC_1
                DO J = 1, NFINE_P
                 L_TRAN = - L_KN * HJ(J) * TRAN(J)
                 L_FUNC = L_ATTN_FINE(NP,J,NP,Q) * TRAN(J) &
                                     +   ATTN_FINE(NP,J) * L_TRAN
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
               L_MULTIPLIERS_P(UT,NP,Q) =  ( TERM1 +   TERM2 ) * L_KN + &
                                         ( L_TERM1 + L_TERM2 ) * KN
              ELSE
               L_FUNC_2 = L_ATTN_P(UT,K,Q)
               L_FUNC_1 = L_ATTN(NP,K,Q) * TRAN_1
               IF ( NFINE_P.GT.0) THEN
                L_SUM = 0.5D0 * L_FUNC_1
                DO J = 1, NFINE_P
                 L_FUNC = L_ATTN_FINE(NP,J,K,Q) * TRAN(J)
                 IF ( J.LT.NFINE_P ) L_SUM = L_SUM + L_FUNC
                 IF ( J.EQ.NFINE_P ) L_SUM = L_SUM + 0.5D0*L_FUNC
                ENDDO
                L_TERM1 = L_SUM * STEP_F
                L_TERM2 = ( L_FUNC + L_FUNC_2 ) * 0.5D0 * STEP_P
               ELSE
                L_TERM1 = 0.5D0 * L_FUNC_1 * STEP_P
                L_TERM2 = 0.5D0 * L_FUNC_2 * STEP_P
               ENDIF
               L_MULTIPLIERS_P(UT,K,Q) =  ( L_TERM1 + L_TERM2 ) * KN
              ENDIF
             ENDDO
            ENDIF
           ENDDO

!  DEBUG SPECIAL CASE
!        WRITE(78,*)'SPECIAL PARTIAL UP',UT,NP,NFINE_P,
!     &          MULTIPLIERS_P(UT),LOSTRANS_P(UT),L_LOSTRANS_P(UT,1)
!        DO K = 1, NLAYERS
!           WRITE(78,*)K, L_MULTIPLIERS_P(UT,K,1)
!        ENDDO

          ENDIF

!  COLUMN LINEARIZATION

          IF ( DO_COLUMN_LINEARIZATION ) THEN
            DO Q = 1, N_TOTALCOLUMN_WFS
             L_KN = L_EXTINCTION(NP,Q)
             L_FUNC_2 = L_ATTN_P(UT,K0,Q)
             L_TRAN_1 = -L_KN * HP * TRAN_1
             L_FUNC_1 = L_ATTN(NP,K0,Q) * TRAN_1 + ATTN(NP) * L_TRAN_1
             IF ( NFINE_P.GT.0) THEN
              L_SUM = 0.5D0 * L_FUNC_1
              DO J = 1, NFINE_P
                L_TRAN = - L_KN * HJ(J) * TRAN(J)
                L_FUNC = L_ATTN_FINE(NP,J,K0,Q)* TRAN(J) &
                       +   ATTN_FINE(NP,J)     * L_TRAN
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
             L_MULTIPLIERS_P(UT,K0,Q) =  ( TERM1 +   TERM2 ) * L_KN + &
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

!  LINEARIZATION.
!  CHECKED 30 JANUARY 2007

        IF ( DO_PROFILE_LINEARIZATION ) THEN
         DO K = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(K) ) THEN
           DO Q = 1, LAYER_VARY_NUMBER(K)
            IF ( K.EQ.N) THEN
             L_KN = RAYCON * L_EXTINCTION(N,Q)
             L_FUNC_2 = L_ATTN(N-1,N,Q) * CSQ_2
             L_TRAN_1 = -L_KN * ARGM_1 * TRAN_1
             L_FUNC_1 = CSQ_1 * ( L_ATTN(N,N,Q) *   TRAN_1 &
                                  + ATTN(N)       * L_TRAN_1 )
             L_SUM = 0.5D0 * ( L_FUNC_1 + L_FUNC_2 )
             DO J = 1, NFINELAYERS
              ARGJ = COT_2 - COT_FINE(J)
              L_TRAN = - L_KN * ARGJ * TRAN(J)
              L_FUNC = CSQ_FINE(J) * ( L_ATTN_FINE(N,J,N,Q) * TRAN(J) &
                                     +   ATTN_FINE(N,J)     * L_TRAN )
              L_SUM = L_SUM + L_FUNC
             ENDDO
             L_LOSTRANS(N,Q)      = L_TRAN_1
             L_MULTIPLIERS(N,N,Q) = STEP * ( L_KN * SUM + KN * L_SUM )
            ELSE
             L_FUNC_2 = L_ATTN(N-1,K,Q) * CSQ_2
             L_FUNC_1 = CSQ_1 * L_ATTN(N,K,Q) * TRAN_1
             L_SUM = 0.5D0 * ( L_FUNC_1 + L_FUNC_2 )
             DO J = 1, NFINELAYERS
              L_FUNC = CSQ_FINE(J) * L_ATTN_FINE(N,J,K,Q) * TRAN(J)
              L_SUM = L_SUM + L_FUNC
             ENDDO
             L_MULTIPLIERS(N,K,Q) = STEP * KN * L_SUM
            ENDIF
           ENDDO
          ENDIF
         ENDDO

!  DEBUG GENERAL CASE
!           WRITE(66,'(I4,1P3E19.9)')N,LOSTRANS(N),MULTIPLIERS(N), &
!                                    L_LOSTRANS(N,1)
!           DO K = 1, NLAYERS
!             WRITE(66,'(4X,I4,1PE19.9)')K,L_MULTIPLIERS(N,K,1)
!           ENDDO

        ENDIF

!  COLUMN LINEARIZATION.

        IF ( DO_COLUMN_LINEARIZATION ) THEN
          DO Q = 1, N_TOTALCOLUMN_WFS
            L_KN = RAYCON * L_EXTINCTION(N,Q)
            L_FUNC_2 = L_ATTN(N-1,K0,Q) * CSQ_2
            L_TRAN_1 = -L_KN * ARGM_1 * TRAN_1
            L_FUNC_1 = CSQ_1 * ( L_ATTN(N,K0,Q) *   TRAN_1 &
                                 + ATTN(N)       * L_TRAN_1 )
            L_SUM = 0.5D0 * ( L_FUNC_1 + L_FUNC_2 )
            DO J = 1, NFINELAYERS
              ARGJ = COT_2 - COT_FINE(J)
              L_TRAN = - L_KN * ARGJ * TRAN(J)
              L_FUNC = CSQ_FINE(J) * ( L_ATTN_FINE(N,J,K0,Q) * TRAN(J) &
                                     +   ATTN_FINE(N,J)      * L_TRAN )
              L_SUM = L_SUM + L_FUNC
            ENDDO
            L_LOSTRANS(N,Q)      = L_TRAN_1
            L_MULTIPLIERS(N,K0,Q) = STEP * ( L_KN * SUM + KN * L_SUM )
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

!  LINEARIZATION.
!    CODED UP 26 SEPTEMBER 2007.

        IF ( DO_PROFILE_LINEARIZATION ) THEN

         DO K = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(K) ) THEN
           DO Q = 1, LAYER_VARY_NUMBER(K)
            IF ( K.EQ.NP) THEN
             L_KN = RAYCON * L_EXTINCTION(NP,Q)
             L_FUNC_2 = L_ATTN_P(UT,NP,Q) * CSQ_2
             L_TRAN_1 = -L_KN * ARGM_1 * TRAN_1
             L_FUNC_1 = CSQ_1 * ( L_ATTN(NP,NP,Q) *   TRAN_1 &
                                  + ATTN(NP)      * L_TRAN_1 )
             IF ( NFINE_P.GT.0) THEN
              L_SUM = 0.5D0 * L_FUNC_1
              DO J = 1, NFINE_P
                ARGJ = COT_2 - COT_FINE(J)
                L_TRAN = - L_KN * ARGJ * TRAN(J)
                L_FUNC = CSQ_FINE(J) * ( L_ATTN_FINE(NP,J,NP,Q)* TRAN(J) &
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
             L_LOSTRANS_P(UT,Q)       = L_TRAN_1
             L_MULTIPLIERS_P(UT,NP,Q) =  ( TERM1 +   TERM2 ) * L_KN + &
                                       ( L_TERM1 + L_TERM2 ) * KN
            ELSE
             L_FUNC_2 = L_ATTN_P(UT,K,Q) * CSQ_2
             L_FUNC_1 = CSQ_1 * L_ATTN(NP,K,Q) * TRAN_1
             IF ( NFINE_P.GT.0) THEN
              L_SUM = 0.5D0 * L_FUNC_1
              DO J = 1, NFINE_P
               L_FUNC = CSQ_FINE(J) * L_ATTN_FINE(NP,J,K,Q) * TRAN(J)
               IF ( J.LT.NFINE_P ) L_SUM = L_SUM + L_FUNC
               IF ( J.EQ.NFINE_P ) L_SUM = L_SUM + 0.5D0*L_FUNC
              ENDDO
              L_TERM1 = L_SUM * STEP_F
              L_TERM2 = ( L_FUNC + L_FUNC_2 ) * 0.5D0 * STEP_P
             ELSE
              L_TERM1 = 0.5D0 * L_FUNC_1 * STEP_P
              L_TERM2 = 0.5D0 * L_FUNC_2 * STEP_P
             ENDIF
             L_MULTIPLIERS_P(UT,K,Q) =  ( L_TERM1 + L_TERM2 ) * KN
            ENDIF
           ENDDO
          ENDIF
         ENDDO

!  DEBUG GENERAL CASE
!        WRITE(79,*)'GENERAL PARTIAL UP',UT,NP,NFINE_P, &
!                MULTIPLIERS_P(UT),LOSTRANS_P(UT),L_LOSTRANS_P(UT,1)
!        DO K = 1, NLAYERS
!           WRITE(79,*)K, L_MULTIPLIERS_P(UT,K,1)
!        ENDDO

        ENDIF

!  COLUMN LINEARIZATION

        IF ( DO_COLUMN_LINEARIZATION ) THEN
          DO Q = 1, N_TOTALCOLUMN_WFS
            L_KN = RAYCON * L_EXTINCTION(NP,Q)
            L_FUNC_2 = L_ATTN_P(UT,K0,Q) * CSQ_2
            L_TRAN_1 = -L_KN * ARGM_1 * TRAN_1
            L_FUNC_1 = CSQ_1 * ( L_ATTN(NP,K0,Q) *   TRAN_1 &
                                 + ATTN(NP)      * L_TRAN_1 )
            IF ( NFINE_P.GT.0) THEN
              L_SUM = 0.5D0 * L_FUNC_1
              DO J = 1, NFINE_P
                ARGJ = COT_2 - COT_FINE(J)
                L_TRAN = - L_KN * ARGJ * TRAN(J)
                L_FUNC = CSQ_FINE(J) * ( L_ATTN_FINE(NP,J,K0,Q)* TRAN(J) &
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
            L_LOSTRANS_P(UT,Q)       = L_TRAN_1
            L_MULTIPLIERS_P(UT,K0,Q) =  ( TERM1 +   TERM2 ) * L_KN + &
                                       ( L_TERM1 + L_TERM2 ) * KN
          ENDDO
        ENDIF

!  FINISH PARTIALS LOOP

      ENDDO

!  FINISH

      RETURN
      END SUBROUTINE L_OUTGOING_INTEGRATION_UP

!

      SUBROUTINE L_OUTGOING_INTEGRATION_DN ( &
        NLAYERS, NFINELAYERS, DO_PARTIALS, &
        DO_PROFILE_LINEARIZATION, LAYER_VARY_FLAG, &
        LAYER_VARY_NUMBER, &
        DO_COLUMN_LINEARIZATION, N_TOTALCOLUMN_WFS, &
        EXTINCTION, L_EXTINCTION, &
        N_PARTIALS, PARTIALS_IDX, PARTIALS_FINEIDX, &
        SUNPATHS, RADII, NTRAVERSE, ALPHA_ALL, &
        SUNPATHS_P,    NTRAVERSE_P,    ALPHA_P, &
        SUNPATHS_FINE, NTRAVERSE_FINE, ALPHA_FINE, &
        MULTIPLIERS,     LOSTRANS, &
        MULTIPLIERS_P,   LOSTRANS_P, &
        L_MULTIPLIERS,   L_LOSTRANS, &
        L_MULTIPLIERS_P, L_LOSTRANS_P )

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

      LOGICAL, INTENT (IN) ::          DO_PROFILE_LINEARIZATION
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
      DOUBLE PRECISION, INTENT (OUT) :: L_MULTIPLIERS &
          (MAXLAYERS,0:MAXLAYERS,MAX_ATMOSWFS)
      DOUBLE PRECISION, INTENT (OUT) :: L_LOSTRANS   (MAXLAYERS,MAX_ATMOSWFS)

      DOUBLE PRECISION, INTENT (OUT) :: MULTIPLIERS_P   (MAX_PARTLAYERS)
      DOUBLE PRECISION, INTENT (OUT) :: LOSTRANS_P      (MAX_PARTLAYERS)
      DOUBLE PRECISION, INTENT (OUT) :: L_MULTIPLIERS_P &
          (MAX_PARTLAYERS,0:MAXLAYERS,MAX_ATMOSWFS)
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

      DOUBLE PRECISION :: L_ATTN ( 0:MAXLAYERS, 0:MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_ATTN_FINE &
          ( MAXLAYERS, MAXFINELAYERS, 0:MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_ATTN_P &
          ( MAX_PARTLAYERS, 0:MAXLAYERS, MAX_ATMOSWFS  )

!  HELP VARIABLES
!  --------------

      INTEGER ::          N, J, K, K0, Q, UT, NP, NFINE_P

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

!  PROFILE LAYER LINEARIZATIONS

      IF ( DO_PROFILE_LINEARIZATION ) THEN
        DO N = 1, NLAYERS
          DO K = 1, NLAYERS
            IF ( LAYER_VARY_FLAG(K) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(K)
                L_MULTIPLIERS(N,K,Q) = ZERO
              ENDDO
            ENDIF
          ENDDO
          IF ( LAYER_VARY_FLAG(N) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_LOSTRANS(N,Q) = ZERO
            ENDDO
          ENDIF
        ENDDO
      ENDIF

      IF ( DO_COLUMN_LINEARIZATION ) THEN
        K0 = 0
        DO N = 1, NLAYERS
          DO Q = 1, N_TOTALCOLUMN_WFS
            L_MULTIPLIERS(N,K0,Q) = ZERO
            L_LOSTRANS(N,Q)       = ZERO
          ENDDO
        ENDDO
      ENDIF

!  PARTIAL LAYER LINEARIZATIONS

      IF ( DO_PARTIALS ) THEN
       IF ( DO_PROFILE_LINEARIZATION ) THEN
        DO UT = 1, N_PARTIALS
          N = PARTIALS_IDX(UT)
          DO K = 1, NLAYERS
            IF ( LAYER_VARY_FLAG(K) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(K)
                L_MULTIPLIERS_P(UT,K,Q) = ZERO
              ENDDO
            ENDIF
          ENDDO
          IF ( LAYER_VARY_FLAG(N) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_LOSTRANS_P(UT,Q) = ZERO
            ENDDO
          ENDIF
        ENDDO
       ENDIF
       IF ( DO_COLUMN_LINEARIZATION ) THEN
        K0 = 0
        DO UT = 1, N_PARTIALS
          DO Q = 1, N_TOTALCOLUMN_WFS
            L_MULTIPLIERS_P(UT,K0,Q) = ZERO
            L_LOSTRANS_P(UT,Q)       = ZERO
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
      IF ( DO_PROFILE_LINEARIZATION.OR.DO_COLUMN_LINEARIZATION ) THEN
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

!  LINEARIZED ATTENUATION FACTORS. PROFILE LINEARIZATION

      IF ( DO_PROFILE_LINEARIZATION ) THEN
        DO N = 0, NLAYERS
          DO K = 1, NTRAVERSE(N)
            IF ( LAYER_VARY_FLAG(K) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(K)
                L_IOP = - L_EXTINCTION(K,Q)
                L_ATTN(N,K,Q) = SUNPATHS(N,K) * ATTN(N) * L_IOP
              ENDDO
            ENDIF
          ENDDO
        ENDDO
        DO N = 1, NLAYERS
         DO J = 1, NFINELAYERS
          DO K = 1, NTRAVERSE_FINE(N,J)
            IF ( LAYER_VARY_FLAG(K) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(K)
                L_IOP = - L_EXTINCTION(K,Q) * ATTN_FINE(N,J)
                L_ATTN_FINE(N,J,K,Q) = SUNPATHS_FINE(N,K,J) * L_IOP
              ENDDO
            ENDIF
          ENDDO
         ENDDO
        ENDDO
        IF ( DO_PARTIALS ) THEN
         DO UT = 1, N_PARTIALS
          DO K = 1, NTRAVERSE_P(UT)
            IF ( LAYER_VARY_FLAG(K) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(K)
                L_IOP = - L_EXTINCTION(K,Q) * ATTN_P(UT)
                L_ATTN_P(UT,K,Q) = SUNPATHS_P(UT,K) * L_IOP
              ENDDO
            ENDIF
          ENDDO
         ENDDO
        ENDIF
      ENDIF

!  LINEARIZED ATTENUATION FACTORS. COLUMN LINEARIZATION

      IF ( DO_COLUMN_LINEARIZATION ) THEN
        DO N = 0, NLAYERS
          DO Q = 1, N_TOTALCOLUMN_WFS
            SUM = ZERO
            DO K = 1, NTRAVERSE(N)
              L_IOP = - L_EXTINCTION(K,Q)
              SUM = SUM + SUNPATHS(N,K) * L_IOP
            ENDDO
            L_ATTN(N,K0,Q) = SUM * ATTN(N)
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
            L_ATTN_FINE(N,J,K0,Q) = SUM * ATTN_FINE(N,J)
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
            L_ATTN_P(UT,K0,Q) = SUM * ATTN_P(UT)
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

!  PROFILE LINEARIZATION.
!  CHECKED 30 JANUARY 2007

        IF ( DO_PROFILE_LINEARIZATION ) THEN
         DO K = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(K) ) THEN
           DO Q = 1, LAYER_VARY_NUMBER(K)
            IF ( K.EQ.N) THEN
             L_KN = RAYCON * L_EXTINCTION(N,Q)
             L_FUNC_2 = L_ATTN(N,N,Q) * CSQ_2
             L_TRAN_1 = -L_KN * ARGM_1 * TRAN_1
             L_FUNC_1 = CSQ_1 * ( L_ATTN(N-1,N,Q) *   TRAN_1 &
                                  + ATTN(N-1)     * L_TRAN_1 )
             L_SUM = 0.5D0 * ( L_FUNC_1 + L_FUNC_2 )
             DO J = NFINELAYERS, 1, -1
              ARGJ = COT_FINE(J) - COT_2
              L_TRAN = - L_KN * ARGJ * TRAN(J)
              L_FUNC = CSQ_FINE(J) * ( L_ATTN_FINE(N,J,N,Q) * TRAN(J) &
                                     +   ATTN_FINE(N,J)     * L_TRAN )
              L_SUM = L_SUM + L_FUNC
             ENDDO
             L_LOSTRANS(N,Q)      = L_TRAN_1
             L_MULTIPLIERS(N,N,Q) = STEP * ( L_KN * SUM + KN * L_SUM )
            ELSE
             L_FUNC_1 = L_ATTN(N-1,K,Q) * CSQ_1 * TRAN_1
             L_FUNC_2 = CSQ_2 * L_ATTN(N,K,Q)
             L_SUM = 0.5D0 * ( L_FUNC_1 + L_FUNC_2 )
             DO J = 1, NFINELAYERS
              L_FUNC = CSQ_FINE(J) * L_ATTN_FINE(N,J,K,Q) * TRAN(J)
              L_SUM = L_SUM + L_FUNC
             ENDDO
             L_MULTIPLIERS(N,K,Q) = STEP * KN * L_SUM
            ENDIF
           ENDDO
          ENDIF
         ENDDO
        ENDIF

!  COLUMN LINEARIZATION

        IF ( DO_COLUMN_LINEARIZATION ) THEN
          DO Q = 1, N_TOTALCOLUMN_WFS
            L_KN = RAYCON * L_EXTINCTION(N,Q)
            L_FUNC_2 = L_ATTN(N,K0,Q) * CSQ_2
            L_TRAN_1 = -L_KN * ARGM_1 * TRAN_1
            L_FUNC_1 = CSQ_1 * ( L_ATTN(N-1,K0,Q) *   TRAN_1 &
                                 + ATTN(N-1)     * L_TRAN_1 )
            L_SUM = 0.5D0 * ( L_FUNC_1 + L_FUNC_2 )
            DO J = NFINELAYERS, 1, -1
              ARGJ = COT_FINE(J) - COT_2
              L_TRAN = - L_KN * ARGJ * TRAN(J)
              L_FUNC = CSQ_FINE(J) * ( L_ATTN_FINE(N,J,K0,Q) * TRAN(J) &
                                     +   ATTN_FINE(N,J)      * L_TRAN )
              L_SUM = L_SUM + L_FUNC
            ENDDO
            L_LOSTRANS(N,Q)       = L_TRAN_1
            L_MULTIPLIERS(N,K0,Q) = STEP * ( L_KN * SUM + KN * L_SUM )
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

!  PROFILE LINEARIZATION.
!    CODED UP 26 SEPTEMBER 2007. REQUIRES CHECKING.....................

        IF ( DO_PROFILE_LINEARIZATION ) THEN
         DO K = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(K) ) THEN
           DO Q = 1, LAYER_VARY_NUMBER(K)
            IF ( K.EQ.NP) THEN
             L_KN = RAYCON * L_EXTINCTION(NP,Q)
             L_FUNC_2 = L_ATTN_P(UT,NP,Q) * CSQ_2
             L_TRAN_1 = -L_KN * ARGM_1 * TRAN_1
             L_FUNC_1 = CSQ_1 * ( L_ATTN(NP-1,NP,Q) *   TRAN_1 &
                                  + ATTN(NP-1)      * L_TRAN_1 )
             IF ( NFINE_P.LT.NFINELAYERS) THEN
              L_SUM = 0.5D0 * L_FUNC_1
              DO J = NFINELAYERS, NFINE_P+1, -1
               ARGJ = COT_FINE(J) - COT_2
               L_TRAN = - L_KN * ARGJ * TRAN(J)
               L_FUNC = CSQ_FINE(J) * ( L_ATTN_FINE(NP,J,NP,Q) * TRAN(J) &
                                     +   ATTN_FINE(NP,J)     * L_TRAN )
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
             L_MULTIPLIERS_P(UT,NP,Q) =  ( TERM1 +   TERM2 ) * L_KN + &
                                       ( L_TERM1 + L_TERM2 ) * KN
            ELSE
             L_FUNC_2 = L_ATTN_P(UT,K,Q) * CSQ_2
             L_FUNC_1 = CSQ_1 * L_ATTN(NP-1,K,Q) *   TRAN_1
             IF ( NFINE_P.LT.NFINELAYERS) THEN
              L_SUM = 0.5D0 * L_FUNC_1
              DO J = NFINELAYERS, NFINE_P+1, -1
               L_FUNC = CSQ_FINE(J) * L_ATTN_FINE(NP,J,K,Q) * TRAN(J)
               IF ( J.GT.NFINE_P+1 ) L_SUM = L_SUM + L_FUNC
               IF ( J.EQ.NFINE_P+1 ) L_SUM = L_SUM + 0.5D0*L_FUNC
              ENDDO
              L_TERM1 = L_SUM * STEP_F
              L_TERM2 = ( L_FUNC + L_FUNC_2 ) * 0.5D0 * STEP_P
             ELSE
              L_TERM1 = 0.5D0 * L_FUNC_1 * STEP_P
              L_TERM2 = 0.5D0 * L_FUNC_2 * STEP_P
             ENDIF
             L_MULTIPLIERS_P(UT,K,Q) = ( L_TERM1 + L_TERM2 ) * KN
            ENDIF
           ENDDO
          ENDIF
         ENDDO
        ENDIF

!  COLUMN LINEARIZATION.

        IF ( DO_COLUMN_LINEARIZATION ) THEN
          DO Q = 1, N_TOTALCOLUMN_WFS
            L_KN = RAYCON * L_EXTINCTION(NP,Q)
            L_FUNC_2 = L_ATTN_P(UT,K0,Q) * CSQ_2
            L_TRAN_1 = -L_KN * ARGM_1 * TRAN_1
            L_FUNC_1 = CSQ_1 * ( L_ATTN(NP-1,K0,Q) *   TRAN_1 &
                                 + ATTN(NP-1)      * L_TRAN_1 )
            IF ( NFINE_P.LT.NFINELAYERS) THEN
              L_SUM = 0.5D0 * L_FUNC_1
              DO J = NFINELAYERS, NFINE_P+1, -1
               ARGJ = COT_FINE(J) - COT_2
               L_TRAN = - L_KN * ARGJ * TRAN(J)
               L_FUNC = CSQ_FINE(J) * ( L_ATTN_FINE(NP,J,K0,Q) * TRAN(J) &
                                     +   ATTN_FINE(NP,J)     * L_TRAN )
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
            L_MULTIPLIERS_P(UT,K0,Q) =  ( TERM1 +   TERM2 ) * L_KN + &
                                      ( L_TERM1 + L_TERM2 ) * KN
          ENDDO
        ENDIF

!  FINISH PARTIALS LOOP

      ENDDO

!  FINISH

      RETURN
      END SUBROUTINE L_OUTGOING_INTEGRATION_DN

!

      SUBROUTINE VLIDORT_LAP_DBCORRECTION ( &
        FLUXMULT, &
        DO_SSCORR_OUTGOING, DO_REFRACTIVE_GEOMETRY, &
        DO_UPWELLING, NSTOKES, NLAYERS, &
        SZA_LOCAL_INPUT, N_USER_RELAZMS, &
        N_USER_LEVELS, &
        DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
        COS_SZANGLES, NBEAMS, &
        N_USER_STREAMS, SZANGLES_ADJUST, &
        PARTLAYERS_OUTFLAG, &
        PARTLAYERS_OUTINDEX, UTAU_LEVEL_MASK_UP, &
        PARTLAYERS_LAYERIDX, N_GEOMETRIES, &
        VZA_OFFSETS, SOLAR_BEAM_OPDEP, &
        DO_REFLECTED_DIRECTBEAM, T_DELT_USERM, &
        T_UTUP_USERM, &
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
        L_DELTAU_SLANT, L_T_DELT_USERM, &
        L_T_UTUP_USERM, &
        DB_CUMSOURCE, EXACTDB_SOURCE, &
        UP_LOSTRANS, UP_LOSTRANS_UT, &
        BOA_ATTN, ATTN_DB_SAVE, &
        L_UP_LOSTRANS, L_UP_LOSTRANS_UT, &
        L_BOA_ATTN, L_EXACTDB_SOURCE, &
        L_DB_CUMSOURCE, PROFILEWF_DB )

!  PREPARES LINEARIZATION OF EXACT DIRECT BEAM REFLECTION
!   LINEARIZATION WITH RESPECT TO PROFILE ATMOSPHERIC VARIABLES

      USE VLIDORT_PARS

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT (IN) :: FLUXMULT
      LOGICAL, INTENT (IN) ::          DO_SSCORR_OUTGOING
      LOGICAL, INTENT (IN) ::          DO_REFRACTIVE_GEOMETRY
      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NLAYERS
      DOUBLE PRECISION, INTENT (IN) :: SZA_LOCAL_INPUT &
          ( 0:MAXLAYERS, MAX_SZANGLES )
      INTEGER, INTENT (IN) ::          N_USER_RELAZMS
      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      DOUBLE PRECISION, INTENT (IN) :: LAMBERTIAN_ALBEDO
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
      DOUBLE PRECISION, INTENT (IN) :: SOLAR_BEAM_OPDEP ( MAXBEAMS )
      LOGICAL, INTENT (IN) ::          DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      LOGICAL, INTENT (IN) ::          LAYER_VARY_FLAG  ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          LAYER_VARY_NUMBER ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_SLANT &
          ( MAX_ATMOSWFS, MAXLAYERS, MAXLAYERS, MAXBEAMS )
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
      DOUBLE PRECISION, INTENT (IN) :: L_BOA_ATTN &
         ( 0:MAXLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )

      DOUBLE PRECISION, INTENT (OUT) :: L_EXACTDB_SOURCE &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: L_DB_CUMSOURCE &
          ( MAX_GEOMETRIES, MAXSTOKES )

      DOUBLE PRECISION, INTENT (INOUT) :: PROFILEWF_DB &
          ( MAX_ATMOSWFS, MAXLAYERS, &
            MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  LOCAL VARIABLES
!  ---------------

      INTEGER ::          N, NUT, NSTART, NUT_PREV, NLEVEL, NELEMENTS
      INTEGER ::          UT, UTA, UM, UA, NC, IB, V, K, Q, O1
      DOUBLE PRECISION :: X0_FLUX, X0_BOA, FACTOR, L_ATTN, TR, LTR, RT

!  FIRST STAGE
!  -----------

!  RETURN IF NO UPWELLING

      IF ( .NOT.DO_UPWELLING ) RETURN

!  START WEIGHTING FUNCTION LOOP

      DO K = 1, NLAYERS
        IF ( LAYER_VARY_FLAG(K) ) THEN
          DO Q = 1, LAYER_VARY_NUMBER(K)

!  INITIALIZE THE OUTPUT RESULTS

            DO V = 1, N_GEOMETRIES
              DO O1 = 1, NSTOKES
                L_EXACTDB_SOURCE(Q,K,V,O1) = ZERO
                DO UTA = 1, N_USER_LEVELS
                  PROFILEWF_DB(Q,K,UTA,V,O1)      = ZERO
                ENDDO
              ENDDO
            ENDDO

!  NEW CODE  R. SPURR, 6 AUGUST 2007. RT SOLUTIONS INC.
!  ====================================================

!  INITIALIZE CUMULATIVE SOURCE TERM

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
                  L_ATTN = L_BOA_ATTN(K,Q,V)
                 ENDIF
                ELSE
                 IF ( DO_REFRACTIVE_GEOMETRY ) THEN
                  X0_BOA = DCOS(SZA_LOCAL_INPUT(NLAYERS,IB)*DEG_TO_RAD)
                 ELSE
                  X0_BOA = COS_SZANGLES(IB)
                 ENDIF
                 L_ATTN = - L_DELTAU_SLANT(Q,NLAYERS,K,IB) * &
                            SOLAR_BEAM_OPDEP(IB)
                ENDIF
                X0_FLUX = FOUR * X0_BOA
                L_ATTN  = L_ATTN * X0_FLUX

!  LOOP OVER ALBEDO KERNELS, ASSIGN REFLECTION

                IF ( DO_LAMBERTIAN_SURFACE ) THEN
                  O1 = 1
                  FACTOR = LAMBERTIAN_ALBEDO * L_ATTN
                  L_EXACTDB_SOURCE(Q,K,V,O1) = FACTOR * FLUXMULT
                ELSE
                  IF ( ATTN_DB_SAVE(V).NE.ZERO ) THEN
                    RT = L_ATTN / ATTN_DB_SAVE(V)
                  ELSE
                    RT = ZERO
                  ENDIF
                  DO O1 = 1, NSTOKES
                    FACTOR = EXACTDB_SOURCE(V,O1) * RT
                    L_EXACTDB_SOURCE(Q,K,V,O1) = FACTOR * FLUXMULT
                  ENDDO
                ENDIF

!  FINISH LOOPS OVER GEOMETRIES

               ENDDO
              ENDIF
             ENDDO
            ENDDO

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
                L_DB_CUMSOURCE(V,O1) = L_EXACTDB_SOURCE(Q,K,V,O1)
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

              DO N = NSTART, NUT, -1
                NC = NLAYERS + 1 - N

!  IF N  = K (VARYING LAYER), ADDITION LINEARIZATION OF TRANSMITTANCE
!  IF N NE K JUST TRANSMITTANCE OF THE LINEARIZED SOURCE

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
                    DO O1 = 1, NELEMENTS
                     L_DB_CUMSOURCE(V,O1) = TR * L_DB_CUMSOURCE(V,O1) &
                          + LTR * DB_CUMSOURCE(V,O1,NC-1)
                    ENDDO
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
                    DO O1 = 1, NELEMENTS
                      L_DB_CUMSOURCE(V,O1) = TR * L_DB_CUMSOURCE(V,O1)
                    ENDDO
                   ENDDO
                  ENDDO
                 ENDDO
                ENDIF

!  END LAYER LOOP

              ENDDO

!  OFFGRID OUTPUT
!  -------------

!  REQUIRE PARTIAL LAYER TRANSMITTANCE
!  IF N = K (VARYING LAYER), REQUIRE ADDITIONAL LINEARIZATION

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
                   DO O1 = 1, NELEMENTS
                    PROFILEWF_DB(Q,K,UTA,V,O1) = TR*L_DB_CUMSOURCE(V,O1) &
                                           + LTR*DB_CUMSOURCE(V,O1,NC)
                   ENDDO
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
                   DO O1 = 1, NELEMENTS
                    PROFILEWF_DB(Q,K,UTA,V,O1) = TR*L_DB_CUMSOURCE(V,O1)
                   ENDDO
                  ENDDO
                 ENDDO
                ENDDO
               ENDIF

!  ONGRID OUTPUT : SET FINAL CUMULATIVE SOURCE DIRECTLY

              ELSE

               DO IB = 1, NBEAMS
                DO UM = 1, N_USER_STREAMS
                 DO UA = 1, N_USER_RELAZMS
                  V = VZA_OFFSETS(IB,UM) + UA
                  DO O1 = 1, NELEMENTS
                   PROFILEWF_DB(Q,K,UTA,V,O1) = L_DB_CUMSOURCE(V,O1)
                  ENDDO
                 ENDDO
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
        ENDIF
      ENDDO

!  FINISH

      RETURN
      END SUBROUTINE VLIDORT_LAP_DBCORRECTION

!

      SUBROUTINE VLIDORT_LAC_DBCORRECTION ( &
        FLUXMULT, &
        DO_SSCORR_OUTGOING, DO_REFRACTIVE_GEOMETRY, &
        DO_UPWELLING, NSTOKES, NLAYERS, &
        SZA_LOCAL_INPUT, N_USER_RELAZMS, &
        N_USER_LEVELS, &
        DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
        COS_SZANGLES, NBEAMS, &
        N_USER_STREAMS, SZANGLES_ADJUST, &
        PARTLAYERS_OUTFLAG, &
        PARTLAYERS_OUTINDEX, UTAU_LEVEL_MASK_UP, &
        PARTLAYERS_LAYERIDX, N_GEOMETRIES, &
        VZA_OFFSETS, DELTAU_SLANT, &
        SOLAR_BEAM_OPDEP, DO_REFLECTED_DIRECTBEAM, &
        T_DELT_USERM, T_UTUP_USERM, &
        N_TOTALCOLUMN_WFS, &
        L_DELTAU_VERT, L_T_DELT_USERM, &
        L_T_UTUP_USERM, &
        DB_CUMSOURCE, EXACTDB_SOURCE, &
        UP_LOSTRANS, UP_LOSTRANS_UT, &
        BOA_ATTN, ATTN_DB_SAVE, &
        L_UP_LOSTRANS, L_UP_LOSTRANS_UT, &
        L_BOA_ATTN, L_EXACTDBC_SOURCE, &
        L_DB_CUMSOURCE, COLUMNWF_DB )

!  PREPARES LINEARIZATION OF EXACT DIRECT BEAM REFLECTION
!   LINEARIZATION WITH RESPECT TO COLUMN ATMOSPHERIC VARIABLES

      USE VLIDORT_PARS

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT (IN) :: FLUXMULT
      LOGICAL, INTENT (IN) ::          DO_SSCORR_OUTGOING
      LOGICAL, INTENT (IN) ::          DO_REFRACTIVE_GEOMETRY
      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NLAYERS
      DOUBLE PRECISION, INTENT (IN) :: SZA_LOCAL_INPUT &
          ( 0:MAXLAYERS, MAX_SZANGLES )
      INTEGER, INTENT (IN) ::          N_USER_RELAZMS
      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      DOUBLE PRECISION, INTENT (IN) :: LAMBERTIAN_ALBEDO
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
      DOUBLE PRECISION, INTENT (IN) :: SOLAR_BEAM_OPDEP ( MAXBEAMS )
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
      DOUBLE PRECISION, INTENT (IN) :: L_BOA_ATTN &
         ( 0:MAXLAYERS, MAX_ATMOSWFS, MAX_GEOMETRIES )

      DOUBLE PRECISION, INTENT (OUT) :: L_EXACTDBC_SOURCE &
          ( MAX_ATMOSWFS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: L_DB_CUMSOURCE &
          ( MAX_GEOMETRIES, MAXSTOKES )

      DOUBLE PRECISION, INTENT (INOUT) :: COLUMNWF_DB &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  LOCAL VARIABLES
!  ---------------

      INTEGER ::          N, NUT, NSTART, NUT_PREV, NLEVEL, NELEMENTS
      INTEGER ::          UT, UTA, UM, UA, NC, IB, V, K, K0, Q, O1
      DOUBLE PRECISION :: X0_FLUX, X0_BOA, FACTOR, L_ATTN, TR, LTR, RT

!  FIRST STAGE
!  -----------

!  RETURN IF NO UPWELLING

      IF ( .NOT.DO_UPWELLING ) RETURN

!  WEIGHTING FUNCTION INDEX

      K0 = 0

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

!  START GEOMETRY LOOPS

        DO UM = 1, N_USER_STREAMS
          DO IB = 1, NBEAMS
            IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
              DO UA = 1, N_USER_RELAZMS
                V = VZA_OFFSETS(IB,UM) + UA

!  BEAM ATTENUATION AND REFLECTION

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
                   L_ATTN = L_ATTN - L_DELTAU_VERT(Q,K) * &
                                        DELTAU_SLANT(NLAYERS,K,IB)
                 ENDDO
                 L_ATTN = L_ATTN * SOLAR_BEAM_OPDEP(IB)
                ENDIF
                X0_FLUX = FOUR * X0_BOA
                L_ATTN  = L_ATTN * X0_FLUX


!  ASSIGN SURFACE

                IF ( DO_LAMBERTIAN_SURFACE ) THEN
                  O1 = 1
                  FACTOR = LAMBERTIAN_ALBEDO * L_ATTN
                  L_EXACTDBC_SOURCE(Q,V,O1) = FACTOR * FLUXMULT
                ELSE
                  IF ( ATTN_DB_SAVE(V).NE.ZERO ) THEN
                    RT = L_ATTN / ATTN_DB_SAVE(V)
                  ELSE
                    RT = ZERO
                  ENDIF
                  DO O1 = 1, NSTOKES
                    FACTOR = EXACTDB_SOURCE(V,O1) * RT
                    L_EXACTDBC_SOURCE(Q,V,O1) = FACTOR * FLUXMULT
                  ENDDO
                ENDIF

!  FINISH LOOPS OVER GEOMETRIES

              ENDDO
            ENDIF
          ENDDO
        ENDDO

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

!  END LAYER LOOP

          ENDDO

!  OFFGRID OUTPUT
!  -------------

!  REQUIRE PARTIAL LAYER TRANSMITTANCE
!  REQUIRE ADDITIONAL LINEARIZATION

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
                  DO O1 = 1, NELEMENTS
                    COLUMNWF_DB(Q,UTA,V,O1) = TR*L_DB_CUMSOURCE(V,O1) &
                                           + LTR*DB_CUMSOURCE(V,O1,NC)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO

!  ONGRID OUTPUT : SET FINAL CUMULATIVE SOURCE DIRECTLY

          ELSE

            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                DO UA = 1, N_USER_RELAZMS
                  V = VZA_OFFSETS(IB,UM) + UA
                  DO O1 = 1, NELEMENTS
                    COLUMNWF_DB(Q,UTA,V,O1) = L_DB_CUMSOURCE(V,O1)
                  ENDDO
                ENDDO
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

      SUBROUTINE VLIDORT_LS_DBCORRECTION ( &
        FLUXMULT, &
        DO_SSCORR_OUTGOING, DO_UPWELLING, &
        NSTOKES, NLAYERS, &
        N_USER_RELAZMS, N_USER_LEVELS, &
        DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
        FLUXVEC, NBEAMS, &
        N_USER_STREAMS, MUELLER_INDEX, &
        PARTLAYERS_OUTFLAG, &
        PARTLAYERS_OUTINDEX, UTAU_LEVEL_MASK_UP, &
        PARTLAYERS_LAYERIDX, N_GEOMETRIES, &
        VZA_OFFSETS, &
        DO_REFLECTED_DIRECTBEAM, T_DELT_USERM, &
        T_UTUP_USERM, ATTN_DB_SAVE, &
        N_SURFACE_WFS, LS_EXACTDB_BRDFUNC, &
        UP_LOSTRANS, UP_LOSTRANS_UT, &
        LS_EXACTDB_SOURCE, LS_DB_CUMSOURCE, &
        SURFACEWF_DB )

!  PREPARES LINEARIZATION OF EXACT DIRECT BEAM REFLECTION
!   LINEARIZATION WITH RESPECT TO ALL SURFACE PARAMETERS

      USE VLIDORT_PARS

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT (IN) :: FLUXMULT
      LOGICAL, INTENT (IN) ::          DO_SSCORR_OUTGOING
      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_USER_RELAZMS
      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      DOUBLE PRECISION, INTENT (IN) :: LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: FLUXVEC ( MAXSTOKES )
      INTEGER, INTENT (IN) ::          NBEAMS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      INTEGER, INTENT (IN) ::          N_GEOMETRIES
      INTEGER, INTENT (IN) ::          VZA_OFFSETS &
          ( MAX_SZANGLES, MAX_USER_VZANGLES )
      LOGICAL, INTENT (IN) ::          DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: ATTN_DB_SAVE ( MAX_GEOMETRIES )
      INTEGER, INTENT (IN) ::          N_SURFACE_WFS
      DOUBLE PRECISION, INTENT (IN) :: LS_EXACTDB_BRDFUNC &
          ( MAX_SURFACEWFS, MAXSTOKES_SQ, &
            MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: UP_LOSTRANS ( MAXLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (IN) :: UP_LOSTRANS_UT &
          ( MAX_PARTLAYERS, MAX_GEOMETRIES )

      DOUBLE PRECISION, INTENT (OUT) :: LS_EXACTDB_SOURCE &
          ( MAX_SURFACEWFS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: LS_DB_CUMSOURCE &
          ( MAX_GEOMETRIES, MAXSTOKES )

      DOUBLE PRECISION, INTENT (INOUT) :: SURFACEWF_DB &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  LOCAL VARIABLES
!  ---------------

      INTEGER ::          N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER ::          NSURFWFS, NELEMENTS
      INTEGER ::          UT, UTA, UM, UA, NC, IB, V, O1, O2, Q, OM
      DOUBLE PRECISION :: TR, SUM, L_KER, REFL_ATTN

!  FIRST STAGE
!  -----------

!  RETURN IF NO UPWELLING

      IF ( .NOT.DO_UPWELLING ) RETURN

!  LAMBERTIAN: ONLY 1 STOKES COMPONENT, ONLY ONE SURFACE WEIGHTING FUNCTION

      IF ( DO_LAMBERTIAN_SURFACE ) THEN
        NELEMENTS = 1
        NSURFWFS  = 1
      ELSE
        NELEMENTS = NSTOKES
        NSURFWFS  = N_SURFACE_WFS
      ENDIF

!  INITIALISE OUTPUT

      DO V = 1, N_GEOMETRIES
        DO UTA = 1, N_USER_LEVELS
          DO Q = 1, NSURFWFS
            DO O1 = 1, NELEMENTS
              SURFACEWF_DB(Q,UTA,V,O1)  = ZERO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!  GET THE SOURCE FUNCTION LINEARIZATION

      DO UM = 1, N_USER_STREAMS
        DO IB = 1, NBEAMS
          DO UA = 1, N_USER_RELAZMS
            V = VZA_OFFSETS(IB,UM) + UA
            IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
              REFL_ATTN = ATTN_DB_SAVE(V)
              DO Q = 1, NSURFWFS
                DO O1 = 1, NELEMENTS
                  IF ( DO_LAMBERTIAN_SURFACE ) THEN
!  @@@ Rob fix       LS_EXACTDB_SOURCE(Q,V,O1) = REFL_ATTN * LAMBERTIAN_ALBEDO
!      (Albedo WF must be Unnormalized)
                    LS_EXACTDB_SOURCE(Q,V,O1) = REFL_ATTN
                  ELSE
                    SUM = ZERO
                    DO O2 = 1, NSTOKES
                      OM = MUELLER_INDEX(O1,O2)
                      L_KER = LS_EXACTDB_BRDFUNC(Q,O1,UM,UA,IB)
                      SUM = SUM + L_KER * FLUXVEC(O2)
                    ENDDO
                    LS_EXACTDB_SOURCE(Q,V,O1) = REFL_ATTN * SUM
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO

!  ====================================================
!  2. TRANSMITTANCE OF SOURCE TERM: UPWELLING RECURSION
!  ====================================================

      DO Q = 1, NSURFWFS

!  INITIALIZE CUMULATIVE SOURCE TERM

        NC =  0
        DO V = 1, N_GEOMETRIES
         DO O1 = 1, NELEMENTS
          LS_DB_CUMSOURCE(V,O1) = LS_EXACTDB_SOURCE(Q,V,O1)*FLUXMULT
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

          DO N = NSTART, NUT, -1
           NC = NLAYERS + 1 - N
           DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
             DO UA = 1, N_USER_RELAZMS
              V = VZA_OFFSETS(IB,UM) + UA
              IF ( DO_SSCORR_OUTGOING ) THEN
                TR  = UP_LOSTRANS(N,V)
              ELSE
                TR  = T_DELT_USERM(N,UM)
              ENDIF
              DO O1 = 1, NELEMENTS
                LS_DB_CUMSOURCE(V,O1) = TR * LS_DB_CUMSOURCE(V,O1)
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO

!  OFFGRID OUTPUT
!  --------------

!  REQUIRE PARTIAL LAYER TRANSMITTANCE

          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
           UT = PARTLAYERS_OUTINDEX(UTA)
           N  = PARTLAYERS_LAYERIDX(UT)
           DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
             DO UA = 1, N_USER_RELAZMS
              V = VZA_OFFSETS(IB,UM) + UA
              IF ( DO_SSCORR_OUTGOING ) THEN
                TR  = UP_LOSTRANS_UT(UT,V)
              ELSE
                TR  = T_UTUP_USERM(UT,UM)
              ENDIF
              DO O1 = 1, NELEMENTS
               SURFACEWF_DB(Q,UTA,V,O1) = TR * LS_DB_CUMSOURCE(V,O1)
              ENDDO
             ENDDO
            ENDDO
           ENDDO

!  ONGRID OUTPUT : SET FINAL CUMULATIVE SOURCE DIRECTLY

          ELSE
           DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
             DO UA = 1, N_USER_RELAZMS
              V = VZA_OFFSETS(IB,UM) + UA
              DO O1 = 1, NELEMENTS
               SURFACEWF_DB(Q,UTA,V,O1) = LS_DB_CUMSOURCE(V,O1)
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDIF

!  CHECK FOR UPDATING THE RECURSION

          IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT

!  END OPTICAL DEPTH LOOP

        ENDDO

!  END LOOP OVER ALL SURFACE  WEIGHTING FUNCTIONS

      ENDDO

!  FINISH

      RETURN
      END SUBROUTINE VLIDORT_LS_DBCORRECTION

      END MODULE vlidort_l_corrections

