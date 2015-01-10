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

! ##########################################################
! #                                                        #
! # Subroutines in this Module                             #
! #                                                        #
! #     Top level routines--------------                   #
! #                                                        #
! #            VLIDORT_UPUSER_ATMOSWF                      #
! #            VLIDORT_DNUSER_ATMOSWF                      #
! #                                                        #
! #     ------Calling                                      #
! #                                                        #
! #            GET_L_TOA_SOURCE                            #
! #            GET_L_BOA_SOURCE                            #
! #                                                        #
! #            L_WHOLELAYER_STERM_UP                       #
! #            L_WHOLELAYER_STERM_DN                       #
! #            L_PARTLAYER_STERM_UP                        #
! #            L_PARTLAYER_STERM_DN                        #
! #                                                        #
! #       VLIDORT_L_INTEGRATED_OUTPUT (master)             #
! #                                                        #
! #     ------Calling                                      #
! #                                                        #
! #            QUADATMOSWF_LEVEL_UP                        #
! #            QUADATMOSWF_LEVEL_DN                        #
! #            QUADATMOSWF_OFFGRID_UP                      #
! #            QUADATMOSWF_OFFGRID_DN                      #
! #                                                        #
! ##########################################################


      MODULE vlidort_l_wfatmos

      PRIVATE
      PUBLIC :: VLIDORT_UPUSER_ATMOSWF, &
                VLIDORT_DNUSER_ATMOSWF, &
                VLIDORT_L_INTEGRATED_OUTPUT

      CONTAINS

      SUBROUTINE VLIDORT_UPUSER_ATMOSWF ( &
        DO_INCLUDE_SURFACE, DO_INCLUDE_THERMEMISS, &
        DO_INCLUDE_MVOUTPUT, DO_INCLUDE_DIRECTBEAM, &
        SURFACE_FACTOR, FLUX_MULTIPLIER, FOURIER_COMPONENT, &
        IBEAM, LAYER_TO_VARY, NV_PARAMETERS, &
        NSTOKES, NLAYERS, N_USER_LEVELS, &
        N_USER_STREAMS, DO_LAYER_SCATTERING, &
        DO_USER_STREAMS, LOCAL_UM_START, &
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
        UTAU_LEVEL_MASK_UP, PARTLAYERS_LAYERIDX, &
        T_DELT_USERM, T_UTUP_USERM, CUMSOURCE_UP, &
        L_T_DELT_USERM, L_T_UTUP_USERM, &
        DO_SOLAR_SOURCES, DO_QUAD_OUTPUT, &
        NSTREAMS, DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
        BRDF_F, USER_BRDF_F, DO_THERMAL_TRANSONLY, DO_DBCORRECTION, &
        QUAD_WEIGHTS, QUAD_STRMWTS, MUELLER_INDEX, &
        DELTAU_SLANT, T_DELT_DISORDS, T_DELT_EIGEN, &
        K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, &
        LCON, MCON, T_WLOWER, USER_DIRECT_BEAM, &
        L_DELTAU_VERT, L_T_DELT_EIGEN, L_T_DELT_DISORDS, &
        L_SOLA_XPOS, L_SOLB_XNEG, &
        L_WLOWER, NCON, PCON, L_T_WLOWER, &
        DO_MSMODE_VLIDORT, &
        UHOM_UPDN, UHOM_UPUP, UPAR_UP_1, UPAR_UP_2, &
        HMULT_1, HMULT_2, EMULT_UP, &
        L_UHOM_UPDN, L_UHOM_UPUP, L_UPAR_UP_1, L_UPAR_UP_2, &
        L_HMULT_1, L_HMULT_2, L_EMULT_UP, &
        L_LAYER_TSUP_UP, &
        UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP, &
        L_UT_HMULT_UU, L_UT_HMULT_UD, L_UT_EMULT_UP, &
        L_LAYER_TSUP_UTUP, &
        L_BOA_THTONLY_SOURCE, ATMOSWF_F )

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_MVOUTPUT
      INTEGER, INTENT (IN) ::          LAYER_TO_VARY, NV_PARAMETERS
      INTEGER, INTENT (IN) ::          FOURIER_COMPONENT, IBEAM
      DOUBLE PRECISION, INTENT (IN) :: SURFACE_FACTOR
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      LOGICAL, INTENT (IN) ::          DO_LAYER_SCATTERING &
          ( 0:MAXMOMENTS, MAXLAYERS )
      LOGICAL, INTENT (IN) ::          DO_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: CUMSOURCE_UP &
          ( MAX_USER_STREAMS, MAXSTOKES, 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTUP_USERM &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::          DO_QUAD_OUTPUT
      INTEGER, INTENT (IN) ::          NSTREAMS
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      DOUBLE PRECISION, INTENT (IN) :: LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: BRDF_F &
          ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: USER_BRDF_F &
          ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_DBCORRECTION
      DOUBLE PRECISION, INTENT (IN) :: QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_SLANT &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_WLOWER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: USER_DIRECT_BEAM &
          ( MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_EIGEN &
          ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_DISORDS &
          ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_WLOWER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: NCON &
          ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PCON &
          ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_WLOWER &
          ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )
      LOGICAL, INTENT (IN) ::          DO_MSMODE_VLIDORT
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPDN &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPUP &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UPAR_UP_1 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UPAR_UP_2 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_1 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_2 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: EMULT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_UHOM_UPDN ( MAX_USER_STREAMS, &
          MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UHOM_UPUP ( MAX_USER_STREAMS, &
          MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UPAR_UP_1 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UPAR_UP_2 ( MAX_USER_STREAMS, &
          MAXSTOKES, MAXLAYERS, 0:MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_HMULT_1 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_HMULT_2 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_EMULT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_LAYER_TSUP_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_UU &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_UD &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_EMULT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_UT_HMULT_UU &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UT_HMULT_UD &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UT_EMULT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, 0:MAXLAYERS, &
            MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_LAYER_TSUP_UTUP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (OUT) :: L_BOA_THTONLY_SOURCE &
          ( MAXSTREAMS, MAX_ATMOSWFS )

!  Linearized output (Ostensibly output)

      DOUBLE PRECISION, INTENT (INOUT) :: ATMOSWF_F &
          ( MAX_ATMOSWFS, 0:MAXLAYERS, MAX_USER_LEVELS, &
            MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  local variables
!  ---------------

      LOGICAL ::          SFLAG
      INTEGER ::          N, NUT, NSTART, NUT_PREV, NLEVEL, O1
      INTEGER ::          UTA, UM, NV, Q, NC, UT, IB

      DOUBLE PRECISION :: L_CUMUL_SOURCE &
          ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

      DOUBLE PRECISION :: L_BOA_MSSOURCE &
          ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_BOA_DBSOURCE &
          ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

      DOUBLE PRECISION :: L_LAYER_SOURCE &
          ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

      DOUBLE PRECISION :: L_FINAL_SOURCE

!  index

      NV = LAYER_TO_VARY
      IB = IBEAM

!  Zero all Fourier components - New rule, better for safety
!    Only did this for components close to zenith (formerly)

      IF ( DO_USER_STREAMS ) THEN
        DO UTA = 1, N_USER_LEVELS
          DO Q = 1, NV_PARAMETERS
            DO UM = 1, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                ATMOSWF_F(Q,NV,UTA,UM,IB,O1,UPIDX) = ZERO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Initialize post-processing recursion
!  ====================================

      IF ( DO_USER_STREAMS ) THEN

!  Get the linearized BOA source terms (diffuse and direct)

        CALL GET_L_BOA_SOURCE ( &
          DO_INCLUDE_SURFACE, DO_INCLUDE_THERMEMISS, &
          DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_MVOUTPUT, &
          SURFACE_FACTOR, FOURIER_COMPONENT, &
          IBEAM, NV, NV_PARAMETERS, &
          DO_SOLAR_SOURCES, DO_QUAD_OUTPUT, &
          NSTOKES, NSTREAMS, &
          NLAYERS, DO_LAMBERTIAN_SURFACE, &
          LAMBERTIAN_ALBEDO, &
          USER_BRDF_F, DO_THERMAL_TRANSONLY, &
          DO_DBCORRECTION, QUAD_WEIGHTS, &
          QUAD_STRMWTS, N_USER_STREAMS, &
          MUELLER_INDEX, DO_USER_STREAMS, &
          LOCAL_UM_START, DELTAU_SLANT, &
          T_DELT_DISORDS, T_DELT_EIGEN, &
          K_REAL, K_COMPLEX, &
          SOLA_XPOS, SOLB_XNEG, &
          LCON, MCON, &
          T_WLOWER, USER_DIRECT_BEAM, BRDF_F, &
          L_DELTAU_VERT, L_T_DELT_EIGEN, &
          L_T_DELT_DISORDS, &
          L_SOLA_XPOS, L_SOLB_XNEG, &
          L_WLOWER, NCON, PCON, &
          L_T_WLOWER, L_BOA_THTONLY_SOURCE, &
          L_BOA_MSSOURCE, L_BOA_DBSOURCE )

!  Set the cumulative source term equal to the BOA sum

        DO UM = LOCAL_UM_START, N_USER_STREAMS
         DO O1 = 1, NSTOKES
          DO Q = 1, NV_PARAMETERS
            L_CUMUL_SOURCE(UM,O1,Q) = L_BOA_MSSOURCE(UM,O1,Q) + &
                                      L_BOA_DBSOURCE(UM,O1,Q)
          ENDDO
         ENDDO
        ENDDO

      ENDIF

!  Recursion Loop for linearized Post-processing
!  =============================================

!  initialise cumulative source term loop

      NC  = 0
      NUT = 0
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = N_USER_LEVELS, 1, -1

!  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_UP(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only
!    1. Get layer source terms
!    2. Find cumulative source term
!    3. Set multiple scatter source term (MSST) output if flagged

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL + 1
          DO N = NSTART, NUT, -1
            SFLAG = DO_LAYER_SCATTERING(FOURIER_COMPONENT,N)
            NC = NLAYERS + 1 - N

            CALL L_WHOLELAYER_STERM_UP ( &
              IB, N, &
              SFLAG, NV, &
              NV_PARAMETERS, DO_INCLUDE_THERMEMISS, &
              DO_SOLAR_SOURCES, NSTOKES, &
              DO_THERMAL_TRANSONLY, &
              DO_MSMODE_VLIDORT, N_USER_STREAMS, &
              LOCAL_UM_START, K_REAL, K_COMPLEX, &
              LCON, MCON, &
              UHOM_UPDN, UHOM_UPUP, &
              UPAR_UP_1, UPAR_UP_2, &
              HMULT_1, HMULT_2, &
              EMULT_UP, NCON, PCON, &
              L_UHOM_UPDN, L_UHOM_UPUP, &
              L_UPAR_UP_1, L_UPAR_UP_2, &
              L_HMULT_1, L_HMULT_2, &
              L_EMULT_UP, &
              L_LAYER_TSUP_UP, &
              L_LAYER_SOURCE )

            IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO O1 = 1, NSTOKES
                  DO Q = 1, NV_PARAMETERS
                    L_CUMUL_SOURCE(UM,O1,Q) = L_LAYER_SOURCE(UM,O1,Q) &
                     +   T_DELT_USERM(N,UM)   * L_CUMUL_SOURCE(UM,O1,Q) &
                     + L_T_DELT_USERM(N,UM,Q) * CUMSOURCE_UP(UM,O1,NC-1)
                  ENDDO
                ENDDO
              ENDDO
            ELSE IF ( N.NE.NV.AND.NV.NE.0 ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO O1 = 1, NSTOKES
                  DO Q = 1, NV_PARAMETERS
                    L_CUMUL_SOURCE(UM,O1,Q) = L_LAYER_SOURCE(UM,O1,Q) &
                     +   T_DELT_USERM(N,UM)   * L_CUMUL_SOURCE(UM,O1,Q)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDIF

!  Offgrid output
!  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT    = PARTLAYERS_OUTINDEX(UTA)
          N     = PARTLAYERS_LAYERIDX(UT)
          SFLAG = DO_LAYER_SCATTERING(FOURIER_COMPONENT,N)

!  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN
            CALL L_PARTLAYER_STERM_UP ( &
              IB, UT, N, SFLAG, NV, &
              NV_PARAMETERS, DO_INCLUDE_THERMEMISS, &
              DO_SOLAR_SOURCES, NSTOKES, &
              DO_THERMAL_TRANSONLY, &
              DO_MSMODE_VLIDORT, N_USER_STREAMS, &
              LOCAL_UM_START, &
              K_REAL, K_COMPLEX, &
              LCON, MCON, &
              UHOM_UPDN, UHOM_UPUP, &
              UPAR_UP_1, UPAR_UP_2, &
              UT_HMULT_UU, UT_HMULT_UD, &
              UT_EMULT_UP, &
              NCON, PCON, &
              L_UHOM_UPDN, L_UHOM_UPUP, &
              L_UPAR_UP_1, L_UPAR_UP_2, &
              L_UT_HMULT_UU, L_UT_HMULT_UD, &
              L_UT_EMULT_UP, &
              L_LAYER_TSUP_UTUP, &
              L_LAYER_SOURCE )

            IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO O1 = 1, NSTOKES
                  DO Q = 1, NV_PARAMETERS
                    L_FINAL_SOURCE = L_LAYER_SOURCE(UM,O1,Q) &
                    +   T_UTUP_USERM(UT,UM)   * L_CUMUL_SOURCE(UM,O1,Q) &
                    + L_T_UTUP_USERM(UT,UM,Q) *   CUMSOURCE_UP(UM,O1,NC)
                    ATMOSWF_F(Q,NV,UTA,UM,IB,O1,UPIDX) = &
                         FLUX_MULTIPLIER * L_FINAL_SOURCE
                  ENDDO
                ENDDO
              ENDDO
            ELSE IF ( N.NE.NV .AND. NV.NE.0 ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO O1 = 1, NSTOKES
                  DO Q = 1, NV_PARAMETERS
                    L_FINAL_SOURCE = L_LAYER_SOURCE(UM,O1,Q) &
                    +   T_UTUP_USERM(UT,UM)   * L_CUMUL_SOURCE(UM,O1,Q)
                    ATMOSWF_F(Q,NV,UTA,UM,IB,O1,UPIDX) = &
                        FLUX_MULTIPLIER * L_FINAL_SOURCE
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDIF

!  Ongrid output
!  -------------

        ELSE

!  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                DO Q = 1, NV_PARAMETERS
                  L_FINAL_SOURCE = &
                         FLUX_MULTIPLIER * L_CUMUL_SOURCE(UM,O1,Q)
!                  if (DABS(L_FINAL_SOURCE).GT.1.0d-12 ) then
                    ATMOSWF_F(Q,NV,UTA,UM,IB,O1,UPIDX) = L_FINAL_SOURCE
!                  endif
                ENDDO
              ENDDO
            ENDDO
          ENDIF

        ENDIF

!  Check for updating the recursion

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT
        ENDIF

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_UPUSER_ATMOSWF

!

      SUBROUTINE VLIDORT_DNUSER_ATMOSWF ( &
        DO_INCLUDE_THERMEMISS, FLUX_MULTIPLIER, &
        FOURIER_COMPONENT, IBEAM, &
        LAYER_TO_VARY, NV_PARAMETERS, &
        NSTOKES, N_USER_LEVELS, &
        N_USER_STREAMS, DO_LAYER_SCATTERING, &
        DO_USER_STREAMS, LOCAL_UM_START, &
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
        UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX, &
        T_DELT_USERM, T_UTDN_USERM, CUMSOURCE_DN, &
        L_T_DELT_USERM, L_T_UTDN_USERM, &
        DO_SOLAR_SOURCES, DO_THERMAL_TRANSONLY, &
        DO_MSMODE_VLIDORT, &
        K_REAL, K_COMPLEX, &
        LCON, MCON, &
        UHOM_DNDN, UHOM_DNUP, &
        UPAR_DN_1, UPAR_DN_2, &
        HMULT_1, HMULT_2, EMULT_DN, &
        NCON, PCON, &
        L_UHOM_DNDN, L_UHOM_DNUP, &
        L_UPAR_DN_1, L_UPAR_DN_2, &
        L_HMULT_1, L_HMULT_2, L_EMULT_DN, &
        L_LAYER_TSUP_DN, &
        UT_HMULT_DU, UT_HMULT_DD, UT_EMULT_DN, &
        L_UT_HMULT_DU, L_UT_HMULT_DD, &
        L_UT_EMULT_DN, &
        L_LAYER_TSUP_UTDN, &
        ATMOSWF_F )

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS
      INTEGER, INTENT (IN) ::          FOURIER_COMPONENT, IBEAM
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER
      INTEGER, INTENT (IN) ::          LAYER_TO_VARY, NV_PARAMETERS
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      LOGICAL, INTENT (IN) ::          DO_LAYER_SCATTERING &
          ( 0:MAXMOMENTS, MAXLAYERS )
      LOGICAL, INTENT (IN) ::          DO_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: CUMSOURCE_DN &
          ( MAX_USER_STREAMS, MAXSTOKES, 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_USERM &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_MSMODE_VLIDORT
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNDN &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNUP &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UPAR_DN_1 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UPAR_DN_2 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_1 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_2 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: EMULT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: NCON &
          ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PCON &
          ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UHOM_DNDN ( MAX_USER_STREAMS, &
          MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UHOM_DNUP ( MAX_USER_STREAMS, &
          MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UPAR_DN_1 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UPAR_DN_2 ( MAX_USER_STREAMS, &
          MAXSTOKES, MAXLAYERS, 0:MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_HMULT_1 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_HMULT_2 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_EMULT_DN ( MAX_USER_STREAMS, &
                 MAXLAYERS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_LAYER_TSUP_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_DU &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_DD &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_EMULT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_UT_HMULT_DU &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UT_HMULT_DD &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UT_EMULT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, 0:MAXLAYERS, &
            MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_LAYER_TSUP_UTDN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Linearized output (Ostensibly output)

      DOUBLE PRECISION, INTENT (INOUT) :: ATMOSWF_F &
         ( MAX_ATMOSWFS, 0:MAXLAYERS, MAX_USER_LEVELS, &
           MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  local variables
!  ---------------

      LOGICAL ::          SFLAG
      INTEGER ::          N, NUT, NSTART, NUT_PREV, NLEVEL, O1
      INTEGER ::          UTA, UM, NV, Q, NC, UT, IB, M

      DOUBLE PRECISION :: L_CUMUL_SOURCE &
          ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

      DOUBLE PRECISION :: L_TOA_SOURCE &
          ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

      DOUBLE PRECISION :: L_LAYER_SOURCE &
          ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

      DOUBLE PRECISION :: L_FINAL_SOURCE

!  Initialise

      NV = LAYER_TO_VARY
      IB = IBEAM
      M  = FOURIER_COMPONENT

!  Zero all Fourier component output

      IF ( DO_USER_STREAMS ) THEN
        DO UTA = 1, N_USER_LEVELS
          DO Q = 1, NV_PARAMETERS
            DO UM = 1, LOCAL_UM_START
              DO O1 = 1, NSTOKES
                ATMOSWF_F(Q,NV,UTA,UM,IB,O1,DNIDX) = ZERO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Initialize post-processing recursion
!  ====================================

!  Get the linearized TOA source terms

      IF ( DO_USER_STREAMS ) THEN
        CALL GET_L_TOA_SOURCE ( &
          NV_PARAMETERS, NSTOKES, &
          N_USER_STREAMS, LOCAL_UM_START, &
          L_TOA_SOURCE )

        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            DO Q = 1, NV_PARAMETERS
              L_CUMUL_SOURCE(UM,O1,Q) = L_TOA_SOURCE(UM,O1,Q)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Recursion Loop for linearized Post-processing
!  =============================================

!  initialise cumulative source term loop

      NC  = 0
      NUT = 0
      NSTART = 1
      NUT_PREV = NSTART - 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = 1, N_USER_LEVELS

!  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_DN(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only
!    1. Get layer source terms
!    2. Find cumulative source term
!    3. Set multiple scatter source term output if flagged

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL
          DO N = NSTART, NUT
            SFLAG = DO_LAYER_SCATTERING(FOURIER_COMPONENT,N)
            NC = N

            CALL L_WHOLELAYER_STERM_DN ( &
              IB, N, SFLAG, NV, &
              NV_PARAMETERS, DO_INCLUDE_THERMEMISS, &
              DO_SOLAR_SOURCES, NSTOKES, &
              DO_THERMAL_TRANSONLY, &
              DO_MSMODE_VLIDORT, N_USER_STREAMS, &
              LOCAL_UM_START, &
              K_REAL, K_COMPLEX, &
              LCON, MCON, &
              UHOM_DNDN, UHOM_DNUP, &
              UPAR_DN_1, UPAR_DN_2, &
              HMULT_1, HMULT_2, EMULT_DN, &
              NCON, PCON, &
              L_UHOM_DNDN, L_UHOM_DNUP, &
              L_UPAR_DN_1, L_UPAR_DN_2, &
              L_HMULT_1, L_HMULT_2, L_EMULT_DN, &
              L_LAYER_TSUP_DN, &
              L_LAYER_SOURCE )

            IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO O1 = 1, NSTOKES
                  DO Q = 1, NV_PARAMETERS
                    L_CUMUL_SOURCE(UM,O1,Q) = L_LAYER_SOURCE(UM,O1,Q) &
                     +   T_DELT_USERM(N,UM)   * L_CUMUL_SOURCE(UM,O1,Q) &
                     + L_T_DELT_USERM(N,UM,Q) * CUMSOURCE_DN(UM,O1,NC-1)
                  ENDDO
                ENDDO
              ENDDO
            ELSE IF ( N.NE.NV.AND.NV.NE.0 ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO O1 = 1, NSTOKES
                  DO Q = 1, NV_PARAMETERS
                    L_CUMUL_SOURCE(UM,O1,Q) = L_LAYER_SOURCE(UM,O1,Q) &
                     +   T_DELT_USERM(N,UM) * L_CUMUL_SOURCE(UM,O1,Q)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDIF

!  Offgrid output
!  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT    = PARTLAYERS_OUTINDEX(UTA)
          N     = PARTLAYERS_LAYERIDX(UT)
          SFLAG = DO_LAYER_SCATTERING(FOURIER_COMPONENT,N)

!  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN
            CALL L_PARTLAYER_STERM_DN ( &
              IB, UT, N, SFLAG, NV, &
              NV_PARAMETERS, DO_INCLUDE_THERMEMISS, &
              DO_SOLAR_SOURCES, NSTOKES, &
              DO_THERMAL_TRANSONLY, &
              DO_MSMODE_VLIDORT, N_USER_STREAMS, &
              LOCAL_UM_START, K_REAL, K_COMPLEX, &
              LCON, MCON, &
              UHOM_DNDN, UHOM_DNUP, &
              UPAR_DN_1, UPAR_DN_2, &
              UT_HMULT_DU, UT_HMULT_DD, UT_EMULT_DN, &
              NCON, PCON, &
              L_UHOM_DNDN, L_UHOM_DNUP, &
              L_UPAR_DN_1, L_UPAR_DN_2, &
              L_UT_HMULT_DU, L_UT_HMULT_DD, &
              L_UT_EMULT_DN, &
              L_LAYER_TSUP_UTDN, &
              L_LAYER_SOURCE )

            IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO O1 = 1, NSTOKES
                  DO Q = 1, NV_PARAMETERS
                    L_FINAL_SOURCE = L_LAYER_SOURCE(UM,O1,Q) &
                    +   T_UTDN_USERM(UT,UM)   * L_CUMUL_SOURCE(UM,O1,Q) &
                    + L_T_UTDN_USERM(UT,UM,Q) *   CUMSOURCE_DN(UM,O1,NC)
                    ATMOSWF_F(Q,NV,UTA,UM,IB,O1,DNIDX) = &
                      FLUX_MULTIPLIER * L_FINAL_SOURCE
                  ENDDO
                ENDDO
              ENDDO
            ELSE IF ( N.NE.NV.AND.NV.NE.0 ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO O1 = 1, NSTOKES
                  DO Q = 1, NV_PARAMETERS
                    L_FINAL_SOURCE = L_LAYER_SOURCE(UM,O1,Q) &
                    +   T_UTDN_USERM(UT,UM)   * L_CUMUL_SOURCE(UM,O1,Q)
                    ATMOSWF_F(Q,NV,UTA,UM,IB,O1,DNIDX) = &
                      FLUX_MULTIPLIER * L_FINAL_SOURCE
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDIF

!  Ongrid output
!  -------------

        ELSE

!  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                DO Q = 1, NV_PARAMETERS
                  ATMOSWF_F(Q,NV,UTA,UM,IB,O1,DNIDX) = &
                          FLUX_MULTIPLIER * L_CUMUL_SOURCE(UM,O1,Q)
                ENDDO
              ENDDO
            ENDDO
          ENDIF

        ENDIF

!  Check for updating the recursion

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT.NE. NUT_PREV ) NSTART = NUT + 1
          NUT_PREV = NUT
        ENDIF

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_DNUSER_ATMOSWF

!

      SUBROUTINE GET_L_TOA_SOURCE ( &
        NV_PARAMETERS, NSTOKES, &
        N_USER_STREAMS, LOCAL_UM_START, &
        L_TOA_SOURCE )

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          NV_PARAMETERS
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START

      DOUBLE PRECISION, INTENT (OUT) :: L_TOA_SOURCE &
          ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

!  local variables
!  ---------------

      INTEGER :: UM, Q, O1

!  initialise TOA source function
!  ------------------------------

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO Q = 1, NV_PARAMETERS
          DO O1 = 1, NSTOKES
            L_TOA_SOURCE(UM,O1,Q) = ZERO
          ENDDO
        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE GET_L_TOA_SOURCE

!

      SUBROUTINE GET_L_BOA_SOURCE ( &
        DO_INCLUDE_SURFACE, DO_INCLUDE_THERMEMISS, &
        DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_MVOUTPUT, &
        SURFACE_FACTOR, FOURIER_COMPONENT, &
        IBEAM, NV, NV_PARAMETERS, &
        DO_SOLAR_SOURCES, DO_QUAD_OUTPUT, &
        NSTOKES, NSTREAMS, &
        NLAYERS, DO_LAMBERTIAN_SURFACE, &
        LAMBERTIAN_ALBEDO, &
        USER_BRDF_F, DO_THERMAL_TRANSONLY, &
        DO_DBCORRECTION, QUAD_WEIGHTS, &
        QUAD_STRMWTS, N_USER_STREAMS, &
        MUELLER_INDEX, DO_USER_STREAMS, &
        LOCAL_UM_START, DELTAU_SLANT, &
        T_DELT_DISORDS, T_DELT_EIGEN, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        LCON, MCON, &
        T_WLOWER, USER_DIRECT_BEAM, BRDF_F, &
        L_DELTAU_VERT, L_T_DELT_EIGEN, &
        L_T_DELT_DISORDS, &
        L_SOLA_XPOS, L_SOLB_XNEG, &
        L_WLOWER, NCON, PCON, &
        L_T_WLOWER, L_BOA_THTONLY_SOURCE, &
        L_BOA_MSSOURCE, L_BOA_DBSOURCE )

!  Linearized Bottom of the atmosphere source term

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_MVOUTPUT
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS
      DOUBLE PRECISION, INTENT (IN) :: SURFACE_FACTOR
      INTEGER, INTENT (IN) ::          FOURIER_COMPONENT, IBEAM
      INTEGER, INTENT (IN) ::          NV, NV_PARAMETERS
      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::          DO_QUAD_OUTPUT
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      DOUBLE PRECISION, INTENT (IN) :: LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: USER_BRDF_F &
          ( 0:MAXMOMENTS,MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_DBCORRECTION
      DOUBLE PRECISION, INTENT (IN) :: QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      LOGICAL, INTENT (IN) ::          DO_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_SLANT &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_WLOWER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: USER_DIRECT_BEAM &
          ( MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: BRDF_F &
          ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_EIGEN &
          ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_DISORDS &
          ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_WLOWER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: NCON &
          ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PCON &
          ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_WLOWER &
          ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (OUT) :: L_BOA_THTONLY_SOURCE &
          ( MAXSTREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_BOA_MSSOURCE &
          ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_BOA_DBSOURCE &
          ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

!  local variables
!  ---------------

      INTEGER ::          M, N, NN, J, I, UM, NELEMENTS
      INTEGER ::          Q, IB, O1, O2, OM, O11
      INTEGER ::          K, KO1, K0, K1, K2

      LOGICAL ::          DO_QTHTONLY
      DOUBLE PRECISION :: DOWN (MAXSTREAMS,MAXSTOKES)
      DOUBLE PRECISION :: L_DOWN (MAXSTREAMS,MAXSTOKES,MAX_ATMOSWFS)

      DOUBLE PRECISION :: REFLEC, S_REFLEC, L_BEAM, FAC, KMULT
      DOUBLE PRECISION :: SPAR, SHOM_R, SHOM_CR, SHOM
      DOUBLE PRECISION :: HOM1, HOM2, HOM1CR, HOM2CR, HOM3CR
      DOUBLE PRECISION :: LXR, NXR, PXR, LLXR, MLXR
      DOUBLE PRECISION :: LXR1, NXR1, PXR1, LLXR1, MLXR1
      DOUBLE PRECISION :: LXR2, NXR2, LLXR2

!  Starting section
!  ----------------

!  Fourier number, layer number

      M   = FOURIER_COMPONENT
      N   = NLAYERS
      KO1 = K_REAL(N) + 1
      IB  = IBEAM

!  Special flag

      DO_QTHTONLY = ( DO_THERMAL_TRANSONLY ) .AND. &
            ( DO_QUAD_OUTPUT .OR. DO_INCLUDE_MVOUTPUT )

!  initialise linearized BOA source functions

      IF ( DO_USER_STREAMS ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO Q = 1, NV_PARAMETERS
            DO O1 = 1, NSTOKES
              L_BOA_MSSOURCE(UM,O1,Q) = ZERO
              L_BOA_DBSOURCE(UM,O1,Q) = ZERO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Thermal tranmsittance only, special term

      IF ( DO_QTHTONLY ) THEN
        DO I = 1, NSTREAMS
          DO Q = 1, NV_PARAMETERS
            L_BOA_THTONLY_SOURCE(I,Q) = ZERO
          ENDDO
        ENDDO
      ENDIF

!  Number of Elements
!    Only want the (1,1) component for Lambertian

      NELEMENTS = 1
      IF ( .not. DO_LAMBERTIAN_SURFACE ) NELEMENTS = NSTOKES

!  reflectance integrand  a(j).x(j).L(downwelling)(-j)

      IF ( DO_INCLUDE_SURFACE ) THEN

!  1. Thermal Transmittance only
!     %%%%%%%%%%%%%%%%%%%%%%%%%%

!  Thermal transmittance solution, build from TOA downwards

        IF ( DO_THERMAL_TRANSONLY ) THEN

!  Initialise

          O1 = 1
          DO I = 1, NSTREAMS
            DOWN(I,O1) = ZERO
            DO Q = 1, NV_PARAMETERS
              L_DOWN(I,O1,Q) = ZERO
            ENDDO
          ENDDO

!  Build

          DO NN = 1, NLAYERS
            IF ( NV.EQ.0 .OR. NV.EQ.NN ) THEN
              DO I = 1, NSTREAMS
                DO Q = 1, NV_PARAMETERS
                  L_DOWN(I,O1,Q) = &
                         L_DOWN(I,O1,Q) * T_DELT_DISORDS(I,NN) &
                       +   DOWN(I,O1)   * L_T_DELT_DISORDS(I,NN,Q) &
                           + L_T_WLOWER(I,NN,Q)
                ENDDO
              ENDDO
            ELSE
              DO I = 1, NSTREAMS
                DO Q = 1, NV_PARAMETERS
                  L_DOWN(I,O1,Q) = &
                    L_DOWN(I,O1,Q) * T_DELT_DISORDS(I,NN)
                ENDDO
              ENDDO
            ENDIF
            DO I = 1, NSTREAMS
              DOWN(I,O1) = &
                 DOWN(I,O1)*T_DELT_DISORDS(I,NN) + T_WLOWER(I,NN)
            ENDDO
          ENDDO

!  Finalize

          DO I = 1, NSTREAMS
            DO Q = 1, NV_PARAMETERS
              L_DOWN(I,O1,Q) = QUAD_WEIGHTS(I) * L_DOWN(I,O1,Q)
            ENDDO
          ENDDO

!  Continuation point for avoiding the next section

          GO TO 5432

        ENDIF

!  2. Scattering solutions
!     %%%%%%%%%%%%%%%%%%%%

!  Two cases:
!  (a) If  NV = N, this is also the layer that is varying --> Extras!
!      If  NV = 0, same thing
!  (b) If  N > NV with variations in layer NV above N

        IF ( NV .EQ. N .or. NV .eq. 0 ) THEN

!  stream and stokes loops

          DO I = 1, NSTREAMS
            DO O1 = 1, NELEMENTS

!  Parameter loop

              DO Q = 1, NV_PARAMETERS

!  Real homogeneous solutions

               SHOM_R = ZERO
               DO K = 1, K_REAL(N)
                LXR  = LCON(K,N)   *   SOLA_XPOS(I,O1,K,N)
                LLXR = LCON(K,N)   * L_SOLA_XPOS(I,O1,K,N,Q)
                MLXR = MCON(K,N)   * L_SOLB_XNEG(I,O1,K,N,Q)
                NXR  = NCON(K,N,Q) *   SOLA_XPOS(I,O1,K,N)
                PXR  = PCON(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
                HOM1 = ( NXR + LLXR ) *   T_DELT_EIGEN(K,N) &
                             +  LXR   * L_T_DELT_EIGEN(K,N,Q)
                HOM2 = PXR + MLXR
                SHOM_R = SHOM_R + HOM1 + HOM2
               ENDDO

!  Complex homogeneous solutions

               SHOM_CR = ZERO
               KO1 = K_REAL(N) + 1
               DO K = 1, K_COMPLEX(N)
                K0 = 2 * K - 2
                K1 = KO1 + K0
                K2 = K1  + 1

                NXR1  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K1,N) &
                        - NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K2,N)
                NXR2  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K2,N) &
                        + NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K1,N)
                PXR1  =   PCON(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N) &
                        - PCON(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)

                LXR1  =   LCON(K1,N) *   SOLA_XPOS(I,O1,K1,N) &
                        - LCON(K2,N) *   SOLA_XPOS(I,O1,K2,N)
                LXR2  =   LCON(K1,N) *   SOLA_XPOS(I,O1,K2,N) &
                        + LCON(K2,N) *   SOLA_XPOS(I,O1,K1,N)

                LLXR1  =   LCON(K1,N) * L_SOLA_XPOS(I,O1,K1,N,Q) &
                         - LCON(K2,N) * L_SOLA_XPOS(I,O1,K2,N,Q)
                LLXR2  =   LCON(K1,N) * L_SOLA_XPOS(I,O1,K2,N,Q) &
                         + LCON(K2,N) * L_SOLA_XPOS(I,O1,K1,N,Q)
                MLXR1  =   MCON(K1,N) * L_SOLB_XNEG(I,O1,K1,N,Q) &
                         - MCON(K2,N) * L_SOLB_XNEG(I,O1,K2,N,Q)

                HOM1CR =   ( NXR1 + LLXR1 ) *   T_DELT_EIGEN(K1,N) &
                         - ( NXR2 + LLXR2 ) *   T_DELT_EIGEN(K2,N)
                HOM2CR =             LXR1   * L_T_DELT_EIGEN(K1,N,Q) &
                                   - LXR2   * L_T_DELT_EIGEN(K2,N,Q)
                HOM3CR = PXR1 + MLXR1
                SHOM_CR = SHOM_CR + HOM1CR + HOM2CR + HOM3CR

               ENDDO

!  real part and add particular solution

               SHOM = SHOM_R + SHOM_CR


!  Add Particular integral linearization
!    Does not exist if thermal only and N > K, K = 0

               SPAR = ZERO
               IF ( DO_SOLAR_SOURCES .OR. &
            (DO_INCLUDE_THERMEMISS.AND.N.EQ.NV) .OR. &
            (DO_INCLUDE_THERMEMISS.AND.NV.EQ.0) ) THEN
                 SPAR = L_WLOWER(I,O1,N,Q)
               ENDIF

!  Final result

               L_DOWN(I,O1,Q) = SPAR + SHOM

!  End loops over Q, O1 and I

              ENDDO
            ENDDO
          ENDDO

!  otherwise the varying layer is above the boundary layer

        ELSE IF ( NV.LT.N .and. NV.NE.0 ) THEN

!  stream and stokes loops

          DO I = 1, NSTREAMS
            DO O1 = 1, NELEMENTS

!  Parameter loop

              DO Q = 1, NV_PARAMETERS

!  Real homogeneous solutions

               SHOM_R = ZERO
               DO K = 1, K_REAL(N)
                NXR = NCON(K,N,Q) *   SOLA_XPOS(I,O1,K,N)
                PXR = PCON(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
                HOM1 = NXR * T_DELT_EIGEN(K,N)
                HOM2 = PXR
                SHOM_R = SHOM_R + HOM1 + HOM2
               ENDDO

!  Complex homogeneous solutions

               SHOM_CR = ZERO
               KO1 = K_REAL(N) + 1
               DO K = 1, K_COMPLEX(N)
                K0 = 2 * K - 2
                K1 = KO1 + K0
                K2 = K1  + 1
                NXR1  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K1,N) &
                        - NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K2,N)
                NXR2  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K2,N) &
                        + NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K1,N)
                PXR1  =   PCON(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N) &
                        - PCON(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)
                HOM1CR =   NXR1 * T_DELT_EIGEN(K1,N) &
                         - NXR2 * T_DELT_EIGEN(K2,N)
                HOM2CR = PXR1
                SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
               ENDDO

!  real part homogeneous solution

               SHOM = SHOM_R + SHOM_CR

!  Add Particular integral linearization
!    Does not exist if thermal only and N > K, K = 0

               SPAR = ZERO
               IF ( DO_SOLAR_SOURCES .OR. &
            (DO_INCLUDE_THERMEMISS.AND.N.EQ.NV) .OR. &
            (DO_INCLUDE_THERMEMISS.AND.NV.EQ.0) ) THEN
                 SPAR = L_WLOWER(I,O1,N,Q)
               ENDIF

!  Final result

               L_DOWN(I,O1,Q) = SPAR + SHOM

!  End loops over Q, O1 and I

              ENDDO
            ENDDO
          ENDDO

!  End cases

        ENDIF

!  Continuation point

 5432   CONTINUE

!  reflectance integrand  a(j).x(j).L_DOWN(-j)

        DO O1 = 1, NELEMENTS
          DO Q = 1, NV_PARAMETERS
            DO I = 1, NSTREAMS
              L_DOWN(I,O1,Q) = L_DOWN(I,O1,Q) * QUAD_STRMWTS(I)
            ENDDO
          ENDDO
        ENDDO

!  reflected multiple scatter intensity at user defined-angles
!  -----------------------------------------------------------

!  ###### Lambertian reflectance (same for all user-streams)

        IF ( DO_LAMBERTIAN_SURFACE ) THEN
          KMULT = SURFACE_FACTOR * LAMBERTIAN_ALBEDO
          IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
            O1 = 1
            DO Q = 1, NV_PARAMETERS
              REFLEC = ZERO
              DO J = 1, NSTREAMS
                REFLEC = REFLEC + L_DOWN(J,O1,Q)
              ENDDO
              REFLEC = KMULT * REFLEC
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                L_BOA_MSSOURCE(UM,O1,Q) = L_BOA_MSSOURCE(UM,O1,Q) &
                                 + REFLEC
              ENDDO
              IF ( DO_QTHTONLY ) THEN
                DO I = 1, NSTREAMS
                 L_BOA_THTONLY_SOURCE(I,Q) = &
                     L_BOA_THTONLY_SOURCE(I,Q) + REFLEC
                ENDDO
              ENDIF
            ENDDO
          ENDIF
        ENDIF

!  .. integrate with BRDF reflectance function at user angles
!  @@@@@@@@ Rob Fix, 2/9/11, DO_QTHTONLY clause was absent

        IF ( .not. DO_LAMBERTIAN_SURFACE ) THEN
          KMULT = SURFACE_FACTOR
          DO Q = 1, NV_PARAMETERS
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                REFLEC = ZERO
                DO J = 1, NSTREAMS
                  S_REFLEC = ZERO
                  DO O2 = 1, NSTOKES
                    OM = MUELLER_INDEX(O1,O2)
                    S_REFLEC = S_REFLEC + L_DOWN(J,O2,Q) * &
                           USER_BRDF_F(M,OM,UM,J)
                  ENDDO
                  REFLEC = REFLEC + S_REFLEC
                ENDDO
                REFLEC = KMULT * REFLEC
                L_BOA_MSSOURCE(UM,O1,Q) = &
                       L_BOA_MSSOURCE(UM,O1,Q) + REFLEC
              ENDDO
            ENDDO
            IF ( DO_QTHTONLY ) THEN
              O11 = 1
              DO I = 1, NSTREAMS
                REFLEC = ZERO
                DO J = 1, NSTREAMS
                  REFLEC = REFLEC + L_DOWN(J,O11,Q) * BRDF_F(M,O11,I,J)
                ENDDO
                REFLEC = KMULT * REFLEC
                L_BOA_THTONLY_SOURCE(I,Q) = &
                    L_BOA_THTONLY_SOURCE(I,Q) + REFLEC
              ENDDO
            ENDIF
          ENDDO
        ENDIF

!  Add direct beam if flagged
!    For NV > 0, profile weighting functions
!    For NV = 0, Bulk (column) weighting functions

        IF ( DO_INCLUDE_DIRECTBEAM.AND..NOT.DO_DBCORRECTION ) THEN
          IF ( NV .NE. 0 ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NELEMENTS
                FAC = - USER_DIRECT_BEAM(UM,IB,O1)*DELTAU_SLANT(N,NV,IB)
                DO Q = 1, NV_PARAMETERS
                  L_BEAM = L_DELTAU_VERT(Q,NV) * FAC
                  L_BOA_DBSOURCE(UM,O1,Q) = L_BEAM
                ENDDO
              ENDDO
            ENDDO
          ELSE
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NELEMENTS
                FAC = - USER_DIRECT_BEAM(UM,IB,O1)
                DO Q = 1, NV_PARAMETERS
                  L_BEAM = ZERO
                  DO K1 = 1, N
                    L_BEAM = L_BEAM + L_DELTAU_VERT(Q,K1) &
                                     * DELTAU_SLANT(N,K1,IB)
                   ENDDO
                  L_BOA_DBSOURCE(UM,O1,Q) = L_BEAM * FAC
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDIF

!  End inclusion of surface terms

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE GET_L_BOA_SOURCE

!

      SUBROUTINE L_WHOLELAYER_STERM_UP ( &
        IBEAM, GIVEN_LAYER, &
        SOURCETERM_FLAG, LAYER_TO_VARY, &
        NV_PARAMETERS, DO_INCLUDE_THERMEMISS, &
        DO_SOLAR_SOURCES, NSTOKES, &
        DO_THERMAL_TRANSONLY, &
        DO_MSMODE_VLIDORT, N_USER_STREAMS, &
        LOCAL_UM_START, K_REAL, K_COMPLEX, &
        LCON, MCON, &
        UHOM_UPDN, UHOM_UPUP, &
        UPAR_UP_1, UPAR_UP_2, &
        HMULT_1, HMULT_2, &
        EMULT_UP, NCON, PCON, &
        L_UHOM_UPDN, L_UHOM_UPUP, &
        L_UPAR_UP_1, L_UPAR_UP_2, &
        L_HMULT_1, L_HMULT_2, &
        L_EMULT_UP, &
        L_LAYER_TSUP_UP, &
        L_LAYERSOURCE )

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          IBEAM, GIVEN_LAYER
      LOGICAL, INTENT (IN) ::          SOURCETERM_FLAG
      INTEGER, INTENT (IN) ::          LAYER_TO_VARY, NV_PARAMETERS
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      INTEGER, INTENT (IN) ::          NSTOKES
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_MSMODE_VLIDORT
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPDN &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPUP &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UPAR_UP_1 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UPAR_UP_2 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_1 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_2 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: EMULT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: NCON &
          ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PCON &
          ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UHOM_UPDN ( MAX_USER_STREAMS, &
          MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UHOM_UPUP ( MAX_USER_STREAMS, &
          MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UPAR_UP_1 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UPAR_UP_2 ( MAX_USER_STREAMS, &
          MAXSTOKES, MAXLAYERS, 0:MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_HMULT_1 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_HMULT_2 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_EMULT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_LAYER_TSUP_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (OUT) :: L_LAYERSOURCE &
          ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

!  local variables
!  ---------------

      INTEGER ::          N, NV, UM, O1, Q, IB
      INTEGER ::          K, KO1, K0, K1, K2
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, SPAR, H1, H2, TM
      DOUBLE PRECISION :: LUXR, MUXR, NUXR, PUXR, LLUXR, MLUXR
      DOUBLE PRECISION :: LUXR1, MUXR1, NUXR1, PUXR1, LLUXR1, MLUXR1
      DOUBLE PRECISION :: LUXR2, MUXR2, NUXR2, PUXR2, LLUXR2, MLUXR2

!   Very important to zero both output terms (bug solved 12/29/05)

      DO Q = 1, NV_PARAMETERS
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            L_LAYERSOURCE(UM,O1,Q)  = ZERO
          ENDDO
       ENDDO
      ENDDO

!  return if no source term
!      Need to go on if thermal transmittances only

      IF ( .NOT. SOURCETERM_FLAG .and. &
           .NOT. DO_THERMAL_TRANSONLY ) RETURN

!  local indices

      N   = GIVEN_LAYER
      KO1 = K_REAL(N) + 1
      NV  = LAYER_TO_VARY
      IB  = IBEAM

!  Avoid this section if thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789

!  Homogeneous solutions
!  =====================

!  Special case when N = NV, or NV = 0 (bulk)
!  ------------------------------------------

      IF ( N .EQ. NV .OR. NV.EQ.0 ) THEN

!  Loop over user angles and STokes

        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES

!  parameter loop

           DO Q = 1, NV_PARAMETERS

!  Real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(N)
              LUXR  = LCON(K,N)   *   UHOM_UPDN(UM,O1,K,N)
              MUXR  = MCON(K,N)   *   UHOM_UPUP(UM,O1,K,N)
              LLUXR = LCON(K,N)   * L_UHOM_UPDN(UM,O1,K,N,Q)
              MLUXR = MCON(K,N)   * L_UHOM_UPUP(UM,O1,K,N,Q)
              NUXR  = NCON(K,N,Q) *   UHOM_UPDN(UM,O1,K,N)
              PUXR  = PCON(K,N,Q) *   UHOM_UPUP(UM,O1,K,N)
              H1 = ( NUXR + LLUXR ) *   HMULT_2(K,UM,N) &
                         +  LUXR    * L_HMULT_2(K,UM,N,Q)
              H2 = ( PUXR + MLUXR ) *   HMULT_1(K,UM,N) &
                         +  MUXR    * L_HMULT_1(K,UM,N,Q)
              SHOM_R = SHOM_R + H1 + H2
            ENDDO

!  Complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K-2
              K1 = KO1 + K0
              K2 = K1  + 1

              LUXR1 =  LCON(K1,N)   *   UHOM_UPDN(UM,O1,K1,N) - &
                       LCON(K2,N)   *   UHOM_UPDN(UM,O1,K2,N)
              LUXR2 =  LCON(K1,N)   *   UHOM_UPDN(UM,O1,K2,N) + &
                       LCON(K2,N)   *   UHOM_UPDN(UM,O1,K1,N)
              MUXR1 =  MCON(K1,N)   *   UHOM_UPUP(UM,O1,K1,N) - &
                       MCON(K2,N)   *   UHOM_UPUP(UM,O1,K2,N)
              MUXR2 =  MCON(K1,N)   *   UHOM_UPUP(UM,O1,K2,N) + &
                       MCON(K2,N)   *   UHOM_UPUP(UM,O1,K1,N)

              LLUXR1 = LCON(K1,N)   * L_UHOM_UPDN(UM,O1,K1,N,Q) - &
                       LCON(K2,N)   * L_UHOM_UPDN(UM,O1,K2,N,Q)
              LLUXR2 = LCON(K1,N)   * L_UHOM_UPDN(UM,O1,K2,N,Q) + &
                       LCON(K2,N)   * L_UHOM_UPDN(UM,O1,K1,N,Q)
              MLUXR1 = MCON(K1,N)   * L_UHOM_UPUP(UM,O1,K1,N,Q) - &
                       MCON(K2,N)   * L_UHOM_UPUP(UM,O1,K2,N,Q)
              MLUXR2 = MCON(K1,N)   * L_UHOM_UPUP(UM,O1,K2,N,Q) + &
                       MCON(K2,N)   * L_UHOM_UPUP(UM,O1,K1,N,Q)

              NUXR1 =  NCON(K1,N,Q) *   UHOM_UPDN(UM,O1,K1,N) - &
                       NCON(K2,N,Q) *   UHOM_UPDN(UM,O1,K2,N)
              NUXR2 =  NCON(K1,N,Q) *   UHOM_UPDN(UM,O1,K2,N) + &
                       NCON(K2,N,Q) *   UHOM_UPDN(UM,O1,K1,N)
              PUXR1 =  PCON(K1,N,Q) *   UHOM_UPUP(UM,O1,K1,N) - &
                       PCON(K2,N,Q) *   UHOM_UPUP(UM,O1,K2,N)
              PUXR2 =  PCON(K1,N,Q) *   UHOM_UPUP(UM,O1,K2,N) + &
                       PCON(K2,N,Q) *   UHOM_UPUP(UM,O1,K1,N)

              H1 =    ( NUXR1 + LLUXR1 ) * HMULT_2(K1,UM,N) &
                    - ( NUXR2 + LLUXR2 ) * HMULT_2(K2,UM,N) &
                          +  LUXR1       * L_HMULT_2(K1,UM,N,Q) &
                          -  LUXR2       * L_HMULT_2(K2,UM,N,Q)
              H2 =    ( PUXR1 + MLUXR1 ) * HMULT_1(K1,UM,N) &
                    - ( PUXR2 + MLUXR2 ) * HMULT_1(K2,UM,N) &
                          +  MUXR1       * L_HMULT_1(K1,UM,N,Q) &
                          -  MUXR2       * L_HMULT_1(K2,UM,N,Q)

              SHOM_CR = SHOM_CR + H1 + H2

            ENDDO

!  homogeneous contribution

            L_LAYERSOURCE(UM,O1,Q) = SHOM_R + SHOM_CR
           ENDDO
          ENDDO
        ENDDO

!  Other cases when N not equal to NV (only variation of Integ-Cons)
!  ----------------------------------

      ELSE IF ( N.NE.NV .AND. NV.NE.0 ) THEN

!  Loop over user angles and STokes

        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES

!  parameter loop

           DO Q = 1, NV_PARAMETERS

!  Real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(N)
              NUXR  = NCON(K,N,Q) *   UHOM_UPDN(UM,O1,K,N)
              PUXR  = PCON(K,N,Q) *   UHOM_UPUP(UM,O1,K,N)
              H1 = NUXR *   HMULT_2(K,UM,N)
              H2 = PUXR *   HMULT_1(K,UM,N)
              SHOM_R = SHOM_R + H1 + H2
            ENDDO

!  Complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K-2
              K1 = KO1 + K0
              K2 = K1  + 1
              NUXR1 =  NCON(K1,N,Q) *   UHOM_UPDN(UM,O1,K1,N) - &
                       NCON(K2,N,Q) *   UHOM_UPDN(UM,O1,K2,N)
              NUXR2 =  NCON(K1,N,Q) *   UHOM_UPDN(UM,O1,K2,N) + &
                       NCON(K2,N,Q) *   UHOM_UPDN(UM,O1,K1,N)
              PUXR1 =  PCON(K1,N,Q) *   UHOM_UPUP(UM,O1,K1,N) - &
                       PCON(K2,N,Q) *   UHOM_UPUP(UM,O1,K2,N)
              PUXR2 =  PCON(K1,N,Q) *   UHOM_UPUP(UM,O1,K2,N) + &
                       PCON(K2,N,Q) *   UHOM_UPUP(UM,O1,K1,N)

              H1 =    NUXR1 * HMULT_2(K1,UM,N) &
                    - NUXR2 * HMULT_2(K2,UM,N)
              H2 =    PUXR1 * HMULT_1(K1,UM,N) &
                    - PUXR2 * HMULT_1(K2,UM,N)

              SHOM_CR = SHOM_CR + H1 + H2
            ENDDO

!  homogeneous contribution

            L_LAYERSOURCE(UM,O1,Q) = SHOM_R + SHOM_CR

!  End loops over Q, O1 and UM

           ENDDO
          ENDDO
        ENDDO

      ENDIF

!  Continuation point

 6789 continue

!  Add thermal emission term (direct and diffuse)
!     -----Modulus 1 if solar sources are included (taken care of earlier)
!     ----- Linearization only exists if N = K

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        O1 = 1
        TM = ONE
        IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
        IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO Q = 1, NV_PARAMETERS
              L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) &
                     + L_LAYER_TSUP_UP(UM,N,Q)*TM
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Particular and single scatter contributions
!  ===========================================

!  Special case when N = NV, or NV = 0 (bulk)
!  ------------------------------------------

      IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN

!  add particular solution

        DO UM = LOCAL_UM_START, N_USER_STREAMS
         DO O1 = 1, NSTOKES
          DO Q = 1, NV_PARAMETERS
           SPAR = L_UPAR_UP_2(UM,O1,N,NV,Q) *   EMULT_UP(UM,N,IB) + &
                    UPAR_UP_2(UM,O1,N)      * L_EMULT_UP(UM,N,NV,IB,Q)
!           if (nv.eq.22.and.ibeam.eq.1)write(447,*)um,q,&
!              L_LAYERSOURCE(UM,O1,Q),SPAR
           L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SPAR
          ENDDO
         ENDDO
        ENDDO

!  Add single scatter term if flagged

        IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
           DO O1 = 1, NSTOKES
            DO Q = 1, NV_PARAMETERS
             SPAR = L_UPAR_UP_1(UM,O1,N,Q) *   EMULT_UP(UM,N,IB) + &
                      UPAR_UP_1(UM,O1,N)   * L_EMULT_UP(UM,N,NV,IB,Q)
             L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SPAR
            ENDDO
           ENDDO
          ENDDO
        ENDIF

!  Other cases when N > NV
!  -----------------------

      ELSE IF ( N .GT. NV .AND. NV.NE.0 ) THEN

!  add particular solution

        DO UM = LOCAL_UM_START, N_USER_STREAMS
         DO O1 = 1, NSTOKES
          DO Q = 1, NV_PARAMETERS
           SPAR = L_UPAR_UP_2(UM,O1,N,NV,Q) *   EMULT_UP(UM,N,IB) + &
                    UPAR_UP_2(UM,O1,N)      * L_EMULT_UP(UM,N,NV,IB,Q)
           L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SPAR
          ENDDO
         ENDDO
        ENDDO

!  Add single scatter term if flagged

        IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
           DO O1 = 1, NSTOKES
            DO Q = 1, NV_PARAMETERS
             SPAR =  UPAR_UP_1(UM,O1,N) * L_EMULT_UP(UM,N,NV,IB,Q)
             L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SPAR
            ENDDO
           ENDDO
          ENDDO
        ENDIF

!  End layer clause

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE L_WHOLELAYER_STERM_UP

!

      SUBROUTINE L_WHOLELAYER_STERM_DN ( &
        IBEAM, GIVEN_LAYER, &
        SOURCETERM_FLAG, LAYER_TO_VARY, &
        NV_PARAMETERS, DO_INCLUDE_THERMEMISS, &
        DO_SOLAR_SOURCES, NSTOKES, &
        DO_THERMAL_TRANSONLY, &
        DO_MSMODE_VLIDORT, N_USER_STREAMS, &
        LOCAL_UM_START, &
        K_REAL, K_COMPLEX, &
        LCON, MCON, &
        UHOM_DNDN, UHOM_DNUP, &
        UPAR_DN_1, UPAR_DN_2, &
        HMULT_1, HMULT_2, EMULT_DN, &
        NCON, PCON, &
        L_UHOM_DNDN, L_UHOM_DNUP, &
        L_UPAR_DN_1, L_UPAR_DN_2, &
        L_HMULT_1, L_HMULT_2, L_EMULT_DN, &
        L_LAYER_TSUP_DN, &
        L_LAYERSOURCE )

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          IBEAM, GIVEN_LAYER
      LOGICAL, INTENT (IN) ::          SOURCETERM_FLAG
      INTEGER, INTENT (IN) ::          LAYER_TO_VARY, NV_PARAMETERS
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      INTEGER, INTENT (IN) ::          NSTOKES
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_MSMODE_VLIDORT
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNDN &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNUP &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UPAR_DN_1 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UPAR_DN_2 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_1 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_2 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: EMULT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: NCON &
          ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PCON &
          ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UHOM_DNDN ( MAX_USER_STREAMS, &
          MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UHOM_DNUP ( MAX_USER_STREAMS, &
          MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UPAR_DN_1 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UPAR_DN_2 ( MAX_USER_STREAMS, &
          MAXSTOKES, MAXLAYERS, 0:MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_HMULT_1 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_HMULT_2 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_EMULT_DN ( MAX_USER_STREAMS, &
                 MAXLAYERS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_LAYER_TSUP_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (OUT) :: L_LAYERSOURCE &
          ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

!  local variables
!  ---------------

      INTEGER ::          N, NV, UM, O1, Q, IB
      INTEGER ::          K, KO1, K0, K1, K2
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, SPAR, H1, H2, TM
      DOUBLE PRECISION :: LUXR, MUXR, NUXR, PUXR, LLUXR, MLUXR
      DOUBLE PRECISION :: LUXR1, MUXR1, NUXR1, PUXR1, LLUXR1, MLUXR1
      DOUBLE PRECISION :: LUXR2, MUXR2, NUXR2, PUXR2, LLUXR2, MLUXR2

!   Very important to zero both output terms (bug solved 12/29/05)

      DO Q = 1, NV_PARAMETERS
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            L_LAYERSOURCE(UM,O1,Q)       = ZERO
          ENDDO
        ENDDO
      ENDDO

!  return if no source term
!      Need to go on if thermal transmittances only

      IF ( .NOT. SOURCETERM_FLAG .and. &
           .NOT. DO_THERMAL_TRANSONLY ) RETURN

!  local indices

      N   = GIVEN_LAYER
      KO1 = K_REAL(N) + 1
      NV  = LAYER_TO_VARY
      IB  = IBEAM

!  Avoid this section if thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789

!  Homogeneous solutions
!  =====================

!  Special case when N = NV, or NV = 0 (bulk)
!  ------------------------------------------

      IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN

!  Loop over user angles and STokes

        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES

!  parameter loop

           DO Q = 1, NV_PARAMETERS

!  Real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(N)
              LUXR  = LCON(K,N)   *   UHOM_DNDN(UM,O1,K,N)
              MUXR  = MCON(K,N)   *   UHOM_DNUP(UM,O1,K,N)
              LLUXR = LCON(K,N)   * L_UHOM_DNDN(UM,O1,K,N,Q)
              MLUXR = MCON(K,N)   * L_UHOM_DNUP(UM,O1,K,N,Q)
              NUXR  = NCON(K,N,Q) *   UHOM_DNDN(UM,O1,K,N)
              PUXR  = PCON(K,N,Q) *   UHOM_DNUP(UM,O1,K,N)
              H1 = ( NUXR + LLUXR ) *   HMULT_1(K,UM,N) &
                         +  LUXR    * L_HMULT_1(K,UM,N,Q)
              H2 = ( PUXR + MLUXR ) *   HMULT_2(K,UM,N) &
                         +  MUXR    * L_HMULT_2(K,UM,N,Q)
              SHOM_R = SHOM_R + H1 + H2
            ENDDO

!  Complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K-2
              K1 = KO1 + K0
              K2 = K1  + 1

              LUXR1 =  LCON(K1,N)   *   UHOM_DNDN(UM,O1,K1,N) - &
                       LCON(K2,N)   *   UHOM_DNDN(UM,O1,K2,N)
              LUXR2 =  LCON(K1,N)   *   UHOM_DNDN(UM,O1,K2,N) + &
                       LCON(K2,N)   *   UHOM_DNDN(UM,O1,K1,N)
              MUXR1 =  MCON(K1,N)   *   UHOM_DNUP(UM,O1,K1,N) - &
                       MCON(K2,N)   *   UHOM_DNUP(UM,O1,K2,N)
              MUXR2 =  MCON(K1,N)   *   UHOM_DNUP(UM,O1,K2,N) + &
                       MCON(K2,N)   *   UHOM_DNUP(UM,O1,K1,N)

              LLUXR1 = LCON(K1,N)   * L_UHOM_DNDN(UM,O1,K1,N,Q) - &
                       LCON(K2,N)   * L_UHOM_DNDN(UM,O1,K2,N,Q)
              LLUXR2 = LCON(K1,N)   * L_UHOM_DNDN(UM,O1,K2,N,Q) + &
                       LCON(K2,N)   * L_UHOM_DNDN(UM,O1,K1,N,Q)
              MLUXR1 = MCON(K1,N)   * L_UHOM_DNUP(UM,O1,K1,N,Q) - &
                       MCON(K2,N)   * L_UHOM_DNUP(UM,O1,K2,N,Q)
              MLUXR2 = MCON(K1,N)   * L_UHOM_DNUP(UM,O1,K2,N,Q) + &
                       MCON(K2,N)   * L_UHOM_DNUP(UM,O1,K1,N,Q)

              NUXR1 =  NCON(K1,N,Q) *   UHOM_DNDN(UM,O1,K1,N) - &
                       NCON(K2,N,Q) *   UHOM_DNDN(UM,O1,K2,N)
              NUXR2 =  NCON(K1,N,Q) *   UHOM_DNDN(UM,O1,K2,N) + &
                       NCON(K2,N,Q) *   UHOM_DNDN(UM,O1,K1,N)
              PUXR1 =  PCON(K1,N,Q) *   UHOM_DNUP(UM,O1,K1,N) - &
                       PCON(K2,N,Q) *   UHOM_DNUP(UM,O1,K2,N)
              PUXR2 =  PCON(K1,N,Q) *   UHOM_DNUP(UM,O1,K2,N) + &
                       PCON(K2,N,Q) *   UHOM_DNUP(UM,O1,K1,N)

              H1 =    ( NUXR1 + LLUXR1 ) * HMULT_1(K1,UM,N) &
                    - ( NUXR2 + LLUXR2 ) * HMULT_1(K2,UM,N) &
                          +  LUXR1       * L_HMULT_1(K1,UM,N,Q) &
                          -  LUXR2       * L_HMULT_1(K2,UM,N,Q)
              H2 =    ( PUXR1 + MLUXR1 ) * HMULT_2(K1,UM,N) &
                    - ( PUXR2 + MLUXR2 ) * HMULT_2(K2,UM,N) &
                          +  MUXR1       * L_HMULT_2(K1,UM,N,Q) &
                          -  MUXR2       * L_HMULT_2(K2,UM,N,Q)

              SHOM_CR = SHOM_CR + H1 + H2

            ENDDO

!  homogeneous contribution

            L_LAYERSOURCE(UM,O1,Q) = SHOM_R + SHOM_CR

!  End loops over Q, O1 and UM

           ENDDO
          ENDDO
        ENDDO

!  Other cases when N not equal to NV (only variation of Integ-Cons)
!  ----------------------------------

      ELSE IF ( NV.NE.0 .AND. N.NE.NV ) THEN

!  Loop over user angles and STokes

        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES

!  parameter loop

           DO Q = 1, NV_PARAMETERS

!  Real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(N)
              NUXR  = NCON(K,N,Q) *   UHOM_DNDN(UM,O1,K,N)
              PUXR  = PCON(K,N,Q) *   UHOM_DNUP(UM,O1,K,N)
              H1 = NUXR * HMULT_1(K,UM,N)
              H2 = PUXR * HMULT_2(K,UM,N)
              SHOM_R = SHOM_R + H1 + H2
            ENDDO

!  Complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K-2
              K1 = KO1 + K0
              K2 = K1  + 1

              NUXR1 =  NCON(K1,N,Q) *   UHOM_DNDN(UM,O1,K1,N) - &
                       NCON(K2,N,Q) *   UHOM_DNDN(UM,O1,K2,N)
              NUXR2 =  NCON(K1,N,Q) *   UHOM_DNDN(UM,O1,K2,N) + &
                       NCON(K2,N,Q) *   UHOM_DNDN(UM,O1,K1,N)
              PUXR1 =  PCON(K1,N,Q) *   UHOM_DNUP(UM,O1,K1,N) - &
                       PCON(K2,N,Q) *   UHOM_DNUP(UM,O1,K2,N)
              PUXR2 =  PCON(K1,N,Q) *   UHOM_DNUP(UM,O1,K2,N) + &
                       PCON(K2,N,Q) *   UHOM_DNUP(UM,O1,K1,N)

              H1 =     NUXR1 * HMULT_1(K1,UM,N) &
                    -  NUXR2 * HMULT_1(K2,UM,N)
              H2 =     PUXR1 * HMULT_2(K1,UM,N) &
                    -  PUXR2 * HMULT_2(K2,UM,N)

              SHOM_CR = SHOM_CR + H1 + H2

            ENDDO

!  homogeneous contribution

            L_LAYERSOURCE(UM,O1,Q) = SHOM_R + SHOM_CR

!  End loops over Q, O1 and UM

           ENDDO
          ENDDO
        ENDDO

!  End layer clause

      ENDIF

!  Continuation point

 6789 continue

!  Add thermal emission term (direct and diffuse)
!     ----- Modulus 1 if solar sources are included (taken care of earlier)
!     ----- Linearization only exists if N = K, or K = 0 (bulk)

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        O1 = 1
        TM = ONE
        IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
        IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO Q = 1, NV_PARAMETERS
              L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) &
                     + L_LAYER_TSUP_DN(UM,N,Q)*TM
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Particular and single scatter contributions
!  ===========================================

!  Special case when N = NV, or NV = 0 (bulk)
!  -----------------------------------------

      IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN

!  add particular solution

        DO UM = LOCAL_UM_START, N_USER_STREAMS
         DO O1 = 1, NSTOKES
          DO Q = 1, NV_PARAMETERS
           SPAR = L_UPAR_DN_2(UM,O1,N,NV,Q) *   EMULT_DN(UM,N,IB) + &
                    UPAR_DN_2(UM,O1,N)      * L_EMULT_DN(UM,N,NV,IB,Q)
           L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SPAR
          ENDDO
         ENDDO
        ENDDO

!  Add single scatter term if flagged

        IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
           DO O1 = 1, NSTOKES
            DO Q = 1, NV_PARAMETERS
             SPAR = L_UPAR_DN_1(UM,O1,N,Q) *   EMULT_DN(UM,N,IB) + &
                      UPAR_DN_1(UM,O1,N)   * L_EMULT_DN(UM,N,NV,IB,Q)
             L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SPAR
            ENDDO
           ENDDO
          ENDDO
        ENDIF

!  Other cases when N > NV
!  -----------------------

      ELSE IF ( N.GT.NV .AND. NV.NE.0) THEN

!  add particular solution

        DO UM = LOCAL_UM_START, N_USER_STREAMS
         DO O1 = 1, NSTOKES
          DO Q = 1, NV_PARAMETERS
          SPAR = L_UPAR_DN_2(UM,O1,N,NV,Q) *   EMULT_DN(UM,N,IB) + &
                   UPAR_DN_2(UM,O1,N)      * L_EMULT_DN(UM,N,NV,IB,Q)
          L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SPAR
          ENDDO
         ENDDO
        ENDDO

!  Add single scatter term if flagged

        IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
           DO O1 = 1, NSTOKES
            DO Q = 1, NV_PARAMETERS
             SPAR =  UPAR_DN_1(UM,O1,N) * L_EMULT_DN(UM,N,NV,IB,Q)
             L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SPAR
            ENDDO
           ENDDO
          ENDDO
        ENDIF

!  End layer clause

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE L_WHOLELAYER_STERM_DN

!

      SUBROUTINE L_PARTLAYER_STERM_UP ( &
        IBEAM, OFFGRID_INDEX, GIVEN_LAYER, &
        SOURCETERM_FLAG, LAYER_TO_VARY, &
        NV_PARAMETERS, DO_INCLUDE_THERMEMISS, &
        DO_SOLAR_SOURCES, NSTOKES, &
        DO_THERMAL_TRANSONLY, &
        DO_MSMODE_VLIDORT, N_USER_STREAMS, &
        LOCAL_UM_START, &
        K_REAL, K_COMPLEX, &
        LCON, MCON, &
        UHOM_UPDN, UHOM_UPUP, &
        UPAR_UP_1, UPAR_UP_2, &
        UT_HMULT_UU, UT_HMULT_UD, &
        UT_EMULT_UP, &
        NCON, PCON, &
        L_UHOM_UPDN, L_UHOM_UPUP, &
        L_UPAR_UP_1, L_UPAR_UP_2, &
        L_UT_HMULT_UU, L_UT_HMULT_UD, &
        L_UT_EMULT_UP, &
        L_LAYER_TSUP_UTUP, &
        L_LAYERSOURCE )

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          GIVEN_LAYER, IBEAM
      LOGICAL, INTENT (IN) ::          SOURCETERM_FLAG
      INTEGER, INTENT (IN) ::          OFFGRID_INDEX
      INTEGER, INTENT (IN) ::          LAYER_TO_VARY, NV_PARAMETERS
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      INTEGER, INTENT (IN) ::          NSTOKES
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_MSMODE_VLIDORT
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPDN &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPUP &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UPAR_UP_1 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UPAR_UP_2 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_UU &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_UD &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_EMULT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: NCON &
          ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PCON &
          ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UHOM_UPDN ( MAX_USER_STREAMS, &
          MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UHOM_UPUP ( MAX_USER_STREAMS, &
          MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UPAR_UP_1 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UPAR_UP_2 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, 0:MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UT_HMULT_UU &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UT_HMULT_UD &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UT_EMULT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, 0:MAXLAYERS, &
            MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_LAYER_TSUP_UTUP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (OUT) :: L_LAYERSOURCE &
          ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

!  local variables
!  ---------------

      INTEGER ::          N, NV, UM, O1, Q, IB, UT
      INTEGER ::          K, KO1, K0, K1, K2
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, SP, H1, H2, TM
      DOUBLE PRECISION :: LUXR, MUXR, NUXR, PUXR, LLUXR, MLUXR
      DOUBLE PRECISION :: LUXR1, MUXR1, NUXR1, PUXR1, LLUXR1, MLUXR1
      DOUBLE PRECISION :: LUXR2, MUXR2, NUXR2, PUXR2, LLUXR2, MLUXR2

!   Very important to zero both output terms (bug solved 12/29/05)

      DO Q = 1, NV_PARAMETERS
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            L_LAYERSOURCE(UM,O1,Q)       = ZERO
          ENDDO
        ENDDO
      ENDDO

!  return if no source term
!      Need to go on if thermal transmittances only

      IF ( .NOT. SOURCETERM_FLAG .and. &
           .NOT. DO_THERMAL_TRANSONLY ) RETURN

!  local indices

      N   = GIVEN_LAYER
      KO1 = K_REAL(N) + 1
      NV  = LAYER_TO_VARY
      UT  = OFFGRID_INDEX
      IB = IBEAM

!  Avoid this section if thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789

!  Partial layer source function ( Homogeneous/constants variation )
!  =================================================================

!  Special case when N = NV, or NV = 0 (bulk)
!  ------------------------------------------

      IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN

!  Loop over user angles and STokes

        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES

!  parameter loop

           DO Q = 1, NV_PARAMETERS

!  Real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(N)

              LUXR  = LCON(K,N)   *   UHOM_UPDN(UM,O1,K,N)
              MUXR  = MCON(K,N)   *   UHOM_UPUP(UM,O1,K,N)
              LLUXR = LCON(K,N)   * L_UHOM_UPDN(UM,O1,K,N,Q)
              MLUXR = MCON(K,N)   * L_UHOM_UPUP(UM,O1,K,N,Q)
              NUXR  = NCON(K,N,Q) *   UHOM_UPDN(UM,O1,K,N)
              PUXR  = PCON(K,N,Q) *   UHOM_UPUP(UM,O1,K,N)

              H1 = ( NUXR + LLUXR ) *   UT_HMULT_UD(K,UM,UT) &
                         +  LUXR    * L_UT_HMULT_UD(K,UM,UT,Q)
              H2 = ( PUXR + MLUXR ) *   UT_HMULT_UU(K,UM,UT) &
                         +  MUXR    * L_UT_HMULT_UU(K,UM,UT,Q)
              SHOM_R = SHOM_R + H1 + H2

            ENDDO

!  Complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K-2
              K1 = KO1 + K0
              K2 = K1  + 1

              LUXR1 =  LCON(K1,N)   *   UHOM_UPDN(UM,O1,K1,N) - &
                       LCON(K2,N)   *   UHOM_UPDN(UM,O1,K2,N)
              LUXR2 =  LCON(K1,N)   *   UHOM_UPDN(UM,O1,K2,N) + &
                       LCON(K2,N)   *   UHOM_UPDN(UM,O1,K1,N)
              MUXR1 =  MCON(K1,N)   *   UHOM_UPUP(UM,O1,K1,N) - &
                       MCON(K2,N)   *   UHOM_UPUP(UM,O1,K2,N)
              MUXR2 =  MCON(K1,N)   *   UHOM_UPUP(UM,O1,K2,N) + &
                       MCON(K2,N)   *   UHOM_UPUP(UM,O1,K1,N)

              LLUXR1 = LCON(K1,N)   * L_UHOM_UPDN(UM,O1,K1,N,Q) - &
                       LCON(K2,N)   * L_UHOM_UPDN(UM,O1,K2,N,Q)
              LLUXR2 = LCON(K1,N)   * L_UHOM_UPDN(UM,O1,K2,N,Q) + &
                       LCON(K2,N)   * L_UHOM_UPDN(UM,O1,K1,N,Q)
              MLUXR1 = MCON(K1,N)   * L_UHOM_UPUP(UM,O1,K1,N,Q) - &
                       MCON(K2,N)   * L_UHOM_UPUP(UM,O1,K2,N,Q)
              MLUXR2 = MCON(K1,N)   * L_UHOM_UPUP(UM,O1,K2,N,Q) + &
                       MCON(K2,N)   * L_UHOM_UPUP(UM,O1,K1,N,Q)

              NUXR1 =  NCON(K1,N,Q) *   UHOM_UPDN(UM,O1,K1,N) - &
                       NCON(K2,N,Q) *   UHOM_UPDN(UM,O1,K2,N)
              NUXR2 =  NCON(K1,N,Q) *   UHOM_UPDN(UM,O1,K2,N) + &
                       NCON(K2,N,Q) *   UHOM_UPDN(UM,O1,K1,N)
              PUXR1 =  PCON(K1,N,Q) *   UHOM_UPUP(UM,O1,K1,N) - &
                       PCON(K2,N,Q) *   UHOM_UPUP(UM,O1,K2,N)
              PUXR2 =  PCON(K1,N,Q) *   UHOM_UPUP(UM,O1,K2,N) + &
                       PCON(K2,N,Q) *   UHOM_UPUP(UM,O1,K1,N)

              H1 =    ( NUXR1 + LLUXR1 ) * UT_HMULT_UD(K1,UM,UT) &
                    - ( NUXR2 + LLUXR2 ) * UT_HMULT_UD(K2,UM,UT) &
                          +  LUXR1       * L_UT_HMULT_UD(K1,UM,UT,Q) &
                          -  LUXR2       * L_UT_HMULT_UD(K2,UM,UT,Q)
              H2 =    ( PUXR1 + MLUXR1 ) * UT_HMULT_UU(K1,UM,UT) &
                    - ( PUXR2 + MLUXR2 ) * UT_HMULT_UU(K2,UM,UT) &
                          +  MUXR1       * L_UT_HMULT_UU(K1,UM,UT,Q) &
                          -  MUXR2       * L_UT_HMULT_UU(K2,UM,UT,Q)

              SHOM_CR = SHOM_CR + H1 + H2

            ENDDO

!  homogeneous contribution

            L_LAYERSOURCE(UM,O1,Q) = SHOM_R + SHOM_CR

!  End loops over Q, O1 and UM

           ENDDO
          ENDDO
        ENDDO

!  Other cases when N not equal to NV (only variation of Integ-Cons)
!  ----------------------------------

      ELSE IF ( N.NE.NV .AND. NV.NE.0 ) THEN

!  Loop over user angles and STokes

        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES

!  parameter loop

           DO Q = 1, NV_PARAMETERS

!  Real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(N)
              NUXR  = NCON(K,N,Q) *   UHOM_UPDN(UM,O1,K,N)
              PUXR  = PCON(K,N,Q) *   UHOM_UPUP(UM,O1,K,N)
              H1 = NUXR *   UT_HMULT_UD(K,UM,UT)
              H2 = PUXR *   UT_HMULT_UU(K,UM,UT)
              SHOM_R = SHOM_R + H1 + H2
            ENDDO

!  Complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K-2
              K1 = KO1 + K0
              K2 = K1  + 1

              NUXR1 =  NCON(K1,N,Q) *   UHOM_UPDN(UM,O1,K1,N) - &
                       NCON(K2,N,Q) *   UHOM_UPDN(UM,O1,K2,N)
              NUXR2 =  NCON(K1,N,Q) *   UHOM_UPDN(UM,O1,K2,N) + &
                       NCON(K2,N,Q) *   UHOM_UPDN(UM,O1,K1,N)
              PUXR1 =  PCON(K1,N,Q) *   UHOM_UPUP(UM,O1,K1,N) - &
                       PCON(K2,N,Q) *   UHOM_UPUP(UM,O1,K2,N)
              PUXR2 =  PCON(K1,N,Q) *   UHOM_UPUP(UM,O1,K2,N) + &
                       PCON(K2,N,Q) *   UHOM_UPUP(UM,O1,K1,N)

              H1 =    NUXR1 * UT_HMULT_UD(K1,UM,UT) &
                    - NUXR2 * UT_HMULT_UD(K2,UM,UT)
              H2 =    PUXR1 * UT_HMULT_UU(K1,UM,UT) &
                    - PUXR2 * UT_HMULT_UU(K2,UM,UT)

              SHOM_CR = SHOM_CR + H1 + H2

            ENDDO

!  homogeneous contribution

            L_LAYERSOURCE(UM,O1,Q) = SHOM_R + SHOM_CR

!  End loops over Q, O1 and UM

           ENDDO
          ENDDO
        ENDDO

!  End layer clause

      ENDIF

!  Continuation point

 6789 continue

!  Add thermal emission term (direct and diffuse)
!     ----- Modulus 1.0 if solar sources are included (taken care of earlier)
!     ----- Linearization only exists if N = K, or K = 0 (bulk)

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        O1 = 1
        TM = ONE
        IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
        IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO Q = 1, NV_PARAMETERS
              L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) &
                     + L_LAYER_TSUP_UTUP(UM,UT,Q)*TM
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Particular and single scatter contributions
!  ===========================================

!  Special case when N = NV, or NV = 0 (bulk)
!  ------------------------------------------

      IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN

!  add particular solution

        DO UM = LOCAL_UM_START, N_USER_STREAMS
         DO O1 = 1, NSTOKES
          DO Q = 1, NV_PARAMETERS
           SP = L_UPAR_UP_2(UM,O1,N,NV,Q) *   UT_EMULT_UP(UM,UT,IB) + &
                  UPAR_UP_2(UM,O1,N)      * L_UT_EMULT_UP(UM,UT,NV,IB,Q)
           L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SP
          ENDDO
         ENDDO
        ENDDO

!  Add single scatter term if flagged

        IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
           DO O1 = 1, NSTOKES
            DO Q = 1, NV_PARAMETERS
             SP = L_UPAR_UP_1(UM,O1,N,Q) *   UT_EMULT_UP(UM,UT,IB) + &
                    UPAR_UP_1(UM,O1,N)   * L_UT_EMULT_UP(UM,UT,NV,IB,Q)
             L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SP
            ENDDO
           ENDDO
          ENDDO
        ENDIF

!  Other cases when N > NV
!  -----------------------

      ELSE IF ( N.GT.NV .AND. NV.NE.0 ) THEN

!  add particular solution

        DO UM = LOCAL_UM_START, N_USER_STREAMS
         DO O1 = 1, NSTOKES
          DO Q = 1, NV_PARAMETERS
           SP = L_UPAR_UP_2(UM,O1,N,NV,Q) *   UT_EMULT_UP(UM,UT,IB) + &
                  UPAR_UP_2(UM,O1,N)      * L_UT_EMULT_UP(UM,UT,NV,IB,Q)
           L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SP
          ENDDO
         ENDDO
        ENDDO

!  Add single scatter term if flagged

        IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
           DO O1 = 1, NSTOKES
            DO Q = 1, NV_PARAMETERS
             SP =  UPAR_UP_1(UM,O1,N) * L_UT_EMULT_UP(UM,UT,NV,IB,Q)
             L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SP
            ENDDO
           ENDDO
          ENDDO
        ENDIF

!  End layer clause

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE L_PARTLAYER_STERM_UP

!

      SUBROUTINE L_PARTLAYER_STERM_DN ( &
        IBEAM, OFFGRID_INDEX, GIVEN_LAYER, &
        SOURCETERM_FLAG, LAYER_TO_VARY, &
        NV_PARAMETERS, DO_INCLUDE_THERMEMISS, &
        DO_SOLAR_SOURCES, NSTOKES, &
        DO_THERMAL_TRANSONLY, &
        DO_MSMODE_VLIDORT, N_USER_STREAMS, &
        LOCAL_UM_START, K_REAL, K_COMPLEX, &
        LCON, MCON, &
        UHOM_DNDN, UHOM_DNUP, &
        UPAR_DN_1, UPAR_DN_2, &
        UT_HMULT_DU, UT_HMULT_DD, UT_EMULT_DN, &
        NCON, PCON, &
        L_UHOM_DNDN, L_UHOM_DNUP, &
        L_UPAR_DN_1, L_UPAR_DN_2, &
        L_UT_HMULT_DU, L_UT_HMULT_DD, &
        L_UT_EMULT_DN, &
        L_LAYER_TSUP_UTDN, &
        L_LAYERSOURCE )

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          GIVEN_LAYER, IBEAM
      LOGICAL, INTENT (IN) ::          SOURCETERM_FLAG
      INTEGER, INTENT (IN) ::          OFFGRID_INDEX
      INTEGER, INTENT (IN) ::          LAYER_TO_VARY, NV_PARAMETERS
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      INTEGER, INTENT (IN) ::          NSTOKES
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_MSMODE_VLIDORT
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNDN &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNUP &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UPAR_DN_1 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UPAR_DN_2 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_DU &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_DD &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_EMULT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: NCON &
          ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PCON &
          ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UHOM_DNDN &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UHOM_DNUP &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UPAR_DN_1 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UPAR_DN_2 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, 0:MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UT_HMULT_DU &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UT_HMULT_DD &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UT_EMULT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, 0:MAXLAYERS, &
            MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_LAYER_TSUP_UTDN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (OUT) :: L_LAYERSOURCE &
          ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

!  local variables
!  ---------------

      INTEGER ::          N, NV, UM, O1, Q, IB, UT
      INTEGER ::          K, KO1, K0, K1, K2
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, SP, H1, H2, TM
      DOUBLE PRECISION :: LUXR, MUXR, NUXR, PUXR, LLUXR, MLUXR
      DOUBLE PRECISION :: LUXR1, MUXR1, NUXR1, PUXR1, LLUXR1, MLUXR1
      DOUBLE PRECISION :: LUXR2, MUXR2, NUXR2, PUXR2, LLUXR2, MLUXR2

!   Very important to zero both output terms (bug solved 12/29/05)

      DO Q = 1, NV_PARAMETERS
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            L_LAYERSOURCE(UM,O1,Q)       = ZERO
          ENDDO
        ENDDO
      ENDDO

!  return if no source term
!      Need to go on if thermal transmittances only

      IF ( .NOT. SOURCETERM_FLAG .and. &
           .NOT. DO_THERMAL_TRANSONLY ) RETURN

!  local indices

      N   = GIVEN_LAYER
      KO1 = K_REAL(N) + 1
      NV  = LAYER_TO_VARY
      UT  = OFFGRID_INDEX
      IB  = IBEAM

!  Avoid this section if thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789

!  Partial layer source function ( Homogeneous/constants variation )
!  =================================================================

!  Special case when N = NV, or NV = 0 (bulk)
!  ------------------------------------------

      IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN

!  Loop over user angles and STokes

        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES

!  parameter loop

           DO Q = 1, NV_PARAMETERS

!  Real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(N)

              LUXR  = LCON(K,N)   *   UHOM_DNDN(UM,O1,K,N)
              MUXR  = MCON(K,N)   *   UHOM_DNUP(UM,O1,K,N)
              LLUXR = LCON(K,N)   * L_UHOM_DNDN(UM,O1,K,N,Q)
              MLUXR = MCON(K,N)   * L_UHOM_DNUP(UM,O1,K,N,Q)
              NUXR  = NCON(K,N,Q) *   UHOM_DNDN(UM,O1,K,N)
              PUXR  = PCON(K,N,Q) *   UHOM_DNUP(UM,O1,K,N)

              H1 = ( NUXR + LLUXR ) *   UT_HMULT_DD(K,UM,UT) &
                         +  LUXR    * L_UT_HMULT_DD(K,UM,UT,Q)
              H2 = ( PUXR + MLUXR ) *   UT_HMULT_DU(K,UM,UT) &
                         +  MUXR    * L_UT_HMULT_DU(K,UM,UT,Q)
              SHOM_R = SHOM_R + H1 + H2

            ENDDO

!  Complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K-2
              K1 = KO1 + K0
              K2 = K1  + 1

              LUXR1 =  LCON(K1,N)   *   UHOM_DNDN(UM,O1,K1,N) - &
                       LCON(K2,N)   *   UHOM_DNDN(UM,O1,K2,N)
              LUXR2 =  LCON(K1,N)   *   UHOM_DNDN(UM,O1,K2,N) + &
                       LCON(K2,N)   *   UHOM_DNDN(UM,O1,K1,N)
              MUXR1 =  MCON(K1,N)   *   UHOM_DNUP(UM,O1,K1,N) - &
                       MCON(K2,N)   *   UHOM_DNUP(UM,O1,K2,N)
              MUXR2 =  MCON(K1,N)   *   UHOM_DNUP(UM,O1,K2,N) + &
                       MCON(K2,N)   *   UHOM_DNUP(UM,O1,K1,N)

              LLUXR1 = LCON(K1,N)   * L_UHOM_DNDN(UM,O1,K1,N,Q) - &
                       LCON(K2,N)   * L_UHOM_DNDN(UM,O1,K2,N,Q)
              LLUXR2 = LCON(K1,N)   * L_UHOM_DNDN(UM,O1,K2,N,Q) + &
                       LCON(K2,N)   * L_UHOM_DNDN(UM,O1,K1,N,Q)
              MLUXR1 = MCON(K1,N)   * L_UHOM_DNUP(UM,O1,K1,N,Q) - &
                       MCON(K2,N)   * L_UHOM_DNUP(UM,O1,K2,N,Q)
              MLUXR2 = MCON(K1,N)   * L_UHOM_DNUP(UM,O1,K2,N,Q) + &
                       MCON(K2,N)   * L_UHOM_DNUP(UM,O1,K1,N,Q)

              NUXR1 =  NCON(K1,N,Q) *   UHOM_DNDN(UM,O1,K1,N) - &
                       NCON(K2,N,Q) *   UHOM_DNDN(UM,O1,K2,N)
              NUXR2 =  NCON(K1,N,Q) *   UHOM_DNDN(UM,O1,K2,N) + &
                       NCON(K2,N,Q) *   UHOM_DNDN(UM,O1,K1,N)
              PUXR1 =  PCON(K1,N,Q) *   UHOM_DNUP(UM,O1,K1,N) - &
                       PCON(K2,N,Q) *   UHOM_DNUP(UM,O1,K2,N)
              PUXR2 =  PCON(K1,N,Q) *   UHOM_DNUP(UM,O1,K2,N) + &
                       PCON(K2,N,Q) *   UHOM_DNUP(UM,O1,K1,N)

              H1 =    ( NUXR1 + LLUXR1 ) * UT_HMULT_DD(K1,UM,UT) &
                    - ( NUXR2 + LLUXR2 ) * UT_HMULT_DD(K2,UM,UT) &
                          +  LUXR1       * L_UT_HMULT_DD(K1,UM,UT,Q) &
                          -  LUXR2       * L_UT_HMULT_DD(K2,UM,UT,Q)
              H2 =    ( PUXR1 + MLUXR1 ) * UT_HMULT_DU(K1,UM,UT) &
                    - ( PUXR2 + MLUXR2 ) * UT_HMULT_DU(K2,UM,UT) &
                          +  MUXR1       * L_UT_HMULT_DU(K1,UM,UT,Q) &
                          -  MUXR2       * L_UT_HMULT_DU(K2,UM,UT,Q)

              SHOM_CR = SHOM_CR + H1 + H2

            ENDDO

!  homogeneous contribution

            L_LAYERSOURCE(UM,O1,Q) = SHOM_R + SHOM_CR

!  End loops over Q, O1 and UM

           ENDDO
          ENDDO
        ENDDO

!  Other cases when N not equal to NV (only variation of Integ-Cons)
!  ----------------------------------

      ELSE IF ( N.NE.NV .AND. NV.NE.0 ) THEN

!  Loop over user angles and STokes

        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES

!  parameter loop

           DO Q = 1, NV_PARAMETERS

!  Real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(N)
              NUXR  = NCON(K,N,Q) *   UHOM_DNDN(UM,O1,K,N)
              PUXR  = PCON(K,N,Q) *   UHOM_DNUP(UM,O1,K,N)
              H1 = NUXR *   UT_HMULT_DD(K,UM,UT)
              H2 = PUXR *   UT_HMULT_DU(K,UM,UT)
              SHOM_R = SHOM_R + H1 + H2
            ENDDO

!  Complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K-2
              K1 = KO1 + K0
              K2 = K1  + 1

              NUXR1 =  NCON(K1,N,Q) *   UHOM_DNDN(UM,O1,K1,N) - &
                       NCON(K2,N,Q) *   UHOM_DNDN(UM,O1,K2,N)
              NUXR2 =  NCON(K1,N,Q) *   UHOM_DNDN(UM,O1,K2,N) + &
                       NCON(K2,N,Q) *   UHOM_DNDN(UM,O1,K1,N)
              PUXR1 =  PCON(K1,N,Q) *   UHOM_DNUP(UM,O1,K1,N) - &
                       PCON(K2,N,Q) *   UHOM_DNUP(UM,O1,K2,N)
              PUXR2 =  PCON(K1,N,Q) *   UHOM_DNUP(UM,O1,K2,N) + &
                       PCON(K2,N,Q) *   UHOM_DNUP(UM,O1,K1,N)

              H1 =    NUXR1 * UT_HMULT_DD(K1,UM,UT) &
                    - NUXR2 * UT_HMULT_DD(K2,UM,UT)
              H2 =    PUXR1 * UT_HMULT_DU(K1,UM,UT) &
                    - PUXR2 * UT_HMULT_DU(K2,UM,UT)

              SHOM_CR = SHOM_CR + H1 + H2

            ENDDO

!  homogeneous contribution

            L_LAYERSOURCE(UM,O1,Q) = SHOM_R + SHOM_CR

!  End loops over Q, O1 and UM

           ENDDO
          ENDDO
        ENDDO

!  End layer clause

      ENDIF

!  Continuation point

 6789 continue

!  Add thermal emission term (direct and diffuse)
!     ----- Modulus 1.0 if solar sources are included (taken care of earlier)
!     ----- Linearization only exists if N = K, or K = 0 (bulk)

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        O1 = 1
        TM = ONE
        IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
        IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO Q = 1, NV_PARAMETERS
              L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) &
                     + L_LAYER_TSUP_UTDN(UM,UT,Q)*TM
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Particular and single scatter contributions
!  ===========================================

!  Special case when N = NV, or NV = 0 (bulk)
!  ------------------------------------------

      IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN

!  add particular solution

        DO UM = LOCAL_UM_START, N_USER_STREAMS
         DO O1 = 1, NSTOKES
          DO Q = 1, NV_PARAMETERS
           SP = L_UPAR_DN_2(UM,O1,N,NV,Q) *   UT_EMULT_DN(UM,UT,IB) + &
                  UPAR_DN_2(UM,O1,N)      * L_UT_EMULT_DN(UM,UT,NV,IB,Q)
           L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SP
!          if (ut.eq.1)write(*,*)um,q,
!     &        L_UPAR_DN_2(UM,O1,N,NV,Q),L_UT_EMULT_DN(UM,UT,NV,IB,Q)
          ENDDO
         ENDDO
        ENDDO

!  Add single scatter term if flagged

        IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
           DO O1 = 1, NSTOKES
            DO Q = 1, NV_PARAMETERS
             SP = L_UPAR_DN_1(UM,O1,N,Q) *   UT_EMULT_DN(UM,UT,IB) + &
                    UPAR_DN_1(UM,O1,N)   * L_UT_EMULT_DN(UM,UT,NV,IB,Q)
             L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SP
!            write(*,*)ib,n,nv,ut,um,q,'par2'
            ENDDO
           ENDDO
          ENDDO
        ENDIF

!  Other cases when N > NV
!  -----------------------

      ELSE IF ( N.GT.NV .AND. NV.NE.0 ) THEN

!  add particular solution

        DO UM = LOCAL_UM_START, N_USER_STREAMS
         DO O1 = 1, NSTOKES
          DO Q = 1, NV_PARAMETERS
           SP = L_UPAR_DN_2(UM,O1,N,NV,Q) *   UT_EMULT_DN(UM,UT,IB) + &
                  UPAR_DN_2(UM,O1,N)      * L_UT_EMULT_DN(UM,UT,NV,IB,Q)
           L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SP
          ENDDO
         ENDDO
        ENDDO

!  Add single scatter term if flagged

        IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
           DO O1 = 1, NSTOKES
            DO Q = 1, NV_PARAMETERS
             SP =  UPAR_DN_1(UM,O1,N) * L_UT_EMULT_DN(UM,UT,NV,IB,Q)
             L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SP
            ENDDO
           ENDDO
          ENDDO
        ENDIF

!  End layer clause

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE L_PARTLAYER_STERM_DN

!
! @@@ Rob fix 1/31/11, - added FLUX_FACTOR argument

      SUBROUTINE VLIDORT_L_INTEGRATED_OUTPUT ( &
        DO_INCLUDE_MVOUTPUT, DO_INCLUDE_DIRECTBEAM, &
        DO_INCLUDE_THERMEMISS, FLUX_MULTIPLIER, FLUX_FACTOR, & !  Added
        IBEAM, NV, NV_PARAMETERS, &
        NSTOKES, NSTREAMS, N_USER_LEVELS, &
        FLUXVEC, LAYER_PIS_CUTOFF, &
        QUAD_WEIGHTS, QUAD_STRMWTS, &
        N_DIRECTIONS, WHICH_DIRECTIONS, &
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
        UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, &
        PARTLAYERS_LAYERIDX, &
        T_DELT_MUBAR, T_UTDN_MUBAR, &
        INITIAL_TRANS, LOCAL_CSZA, &
        L_INITIAL_TRANS, L_T_DELT_MUBAR, &
        L_T_UTDN_MUBAR, &
        DO_SOLAR_SOURCES, &
        NLAYERS, DO_THERMAL_TRANSONLY, &
        DO_CLASSICAL_SOLUTION, QUAD_STREAMS, &
        T_DELT_DISORDS, T_DISORDS_UTUP, &
        T_UTUP_EIGEN, T_UTDN_EIGEN, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        WUPPER, LCON, MCON, T_WUPPER, &
        BOA_THTONLY_SOURCE, &
        L_T_UTUP_EIGEN, L_T_UTDN_EIGEN, &
        L_T_DELT_DISORDS, L_T_DISORDS_UTUP, &
        L_SOLA_XPOS, L_SOLB_XNEG, &
        L_WUPPER, NCON, PCON, &
        L_T_WUPPER, L_UT_T_PARTIC, &
        L_BOA_THTONLY_SOURCE, &
        T_DELT_EIGEN, L_T_DELT_EIGEN, &
        L_WLOWER, &
        T_DISORDS_UTDN, L_T_DISORDS_UTDN, &
        L_T_WLOWER, T_WLOWER, &
        MINT_ATMOSWF, MINT_ATMOSWF_DIRECT, &
        FLUX_ATMOSWF, FLUX_ATMOSWF_DIRECT )

!  Quadrature output at offgrid or ongrid optical depths
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_MVOUTPUT
      INTEGER, INTENT (IN) ::          NV, NV_PARAMETERS
      INTEGER, INTENT (IN) ::          IBEAM
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER

! @@@ Rob fix 1/31/11, - added FLUX_FACTOR argument 
      DOUBLE PRECISION, INTENT (IN) :: FLUX_FACTOR

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      DOUBLE PRECISION, INTENT (IN) :: FLUXVEC ( MAXSTOKES )
      INTEGER, INTENT (IN) ::          LAYER_PIS_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::          N_DIRECTIONS
      INTEGER, INTENT (IN) ::          WHICH_DIRECTIONS ( MAX_DIRECTIONS )
      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: LOCAL_CSZA ( 0:MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_INITIAL_TRANS &
          ( MAXLAYERS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_MUBAR &
          ( MAXLAYERS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_MUBAR &
          ( MAX_USER_LEVELS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      INTEGER, INTENT (IN) ::          NLAYERS
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTUP &
          ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN &
          ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN &
          ( MAXEVALUES, MAX_PARTLAYERS )
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: WUPPER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_WUPPER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: BOA_THTONLY_SOURCE &
          ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTUP_EIGEN &
          ( MAXEVALUES, MAX_USER_LEVELS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_EIGEN &
          ( MAXEVALUES, MAX_USER_LEVELS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_DISORDS &
          ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DISORDS_UTUP &
          ( MAXSTREAMS, MAX_USER_LEVELS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_WUPPER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: NCON &
          ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PCON &
          ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_WUPPER &
          ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UT_T_PARTIC &
          ( MAXSTREAMS_2, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_BOA_THTONLY_SOURCE &
          (MAXSTREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_EIGEN &
          ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_WLOWER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTDN &
          ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DISORDS_UTDN &
          ( MAXSTREAMS, MAX_USER_LEVELS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_WLOWER &
          ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: T_WLOWER ( MAXSTREAMS_2, MAXLAYERS )

!  Linearized output (Has to be INOUT)

      DOUBLE PRECISION, INTENT (INOUT) :: MINT_ATMOSWF &
          ( MAX_ATMOSWFS, 0:MAXLAYERS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) :: MINT_ATMOSWF_DIRECT &
          ( MAX_ATMOSWFS, 0:MAXLAYERS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) :: FLUX_ATMOSWF &
          ( MAX_ATMOSWFS, 0:MAXLAYERS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) :: FLUX_ATMOSWF_DIRECT &
          ( MAX_ATMOSWFS, 0:MAXLAYERS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES )

!  local variables
!  ---------------

!  Local quadrature output (for debug)

      DOUBLE PRECISION :: QATMOSWF_F &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  Help variables

      INTEGER ::          I, IDIR, WDIR, UTA, UT, Q, N, O1, NLEVEL
      DOUBLE PRECISION :: SMI, SFX, FTRANS
      DOUBLE PRECISION :: L_TRANS, L_DIRECT_FLUX, L_DIRECT_MEANI

!  direction loop

      DO IDIR = 1, N_DIRECTIONS
        WDIR = WHICH_DIRECTIONS(IDIR)

!  Upwelling Jacobian output at Quadrature angles

        IF ( WDIR .EQ. UPIDX ) THEN
          DO UTA = 1, N_USER_LEVELS
            NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
              UT = PARTLAYERS_OUTINDEX(UTA)
              N  = PARTLAYERS_LAYERIDX(UT)

              CALL QUADATMOSWF_OFFGRID_UP ( &
                IBEAM, UTA, UT, N, NV, NV_PARAMETERS, &
                DO_INCLUDE_THERMEMISS, FLUX_MULTIPLIER, &
                DO_SOLAR_SOURCES, NSTOKES, &
                NSTREAMS, NLAYERS, &
                DO_THERMAL_TRANSONLY, &
                DO_CLASSICAL_SOLUTION, QUAD_STREAMS, &
                T_DELT_DISORDS, T_DISORDS_UTUP, &
                T_UTUP_EIGEN, T_UTDN_EIGEN, &
                T_UTDN_MUBAR, &
                K_REAL, K_COMPLEX, &
                SOLA_XPOS, SOLB_XNEG, &
                WUPPER, LCON, MCON, &
                T_WUPPER, &
                BOA_THTONLY_SOURCE, &
                L_T_UTUP_EIGEN, L_T_UTDN_EIGEN, &
                L_T_DELT_DISORDS, L_T_DISORDS_UTUP, &
                L_T_UTDN_MUBAR, &
                L_SOLA_XPOS, L_SOLB_XNEG, &
                L_WUPPER, NCON, PCON, &
                L_T_WUPPER, L_UT_T_PARTIC, &
                L_BOA_THTONLY_SOURCE, &
                QATMOSWF_F )

            ELSE
              CALL QUADATMOSWF_LEVEL_UP ( &
                UTA, NLEVEL, NV, NV_PARAMETERS, &
                FLUX_MULTIPLIER, &
                NSTOKES, NSTREAMS, NLAYERS, &
                DO_THERMAL_TRANSONLY, &
                QUAD_STREAMS, &
                T_DELT_DISORDS, T_DELT_EIGEN, &
                K_REAL, K_COMPLEX, &
                SOLA_XPOS, SOLB_XNEG, &
                LCON, MCON, T_WUPPER, &
                BOA_THTONLY_SOURCE, &
                L_T_DELT_EIGEN, L_T_DELT_DISORDS, &
                L_SOLA_XPOS, L_SOLB_XNEG, &
                L_WUPPER, L_WLOWER, &
                NCON, PCON, &
                L_T_WUPPER, L_BOA_THTONLY_SOURCE, &
                QATMOSWF_F )

            ENDIF
          ENDDO
        ENDIF

!  Downwelling Jacobian output at Quadrature angles

        IF ( WDIR .EQ. DNIDX ) THEN
          DO UTA = 1, N_USER_LEVELS
            NLEVEL = UTAU_LEVEL_MASK_DN(UTA)
            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
              UT = PARTLAYERS_OUTINDEX(UTA)
              N  = PARTLAYERS_LAYERIDX(UT)

              CALL QUADATMOSWF_OFFGRID_DN ( &
                IBEAM, UTA, UT, N, NV, NV_PARAMETERS, &
                DO_INCLUDE_THERMEMISS, FLUX_MULTIPLIER, &
                DO_SOLAR_SOURCES, NSTOKES, &
                NSTREAMS, DO_THERMAL_TRANSONLY, &
                DO_CLASSICAL_SOLUTION, QUAD_STREAMS, &
                T_DELT_DISORDS, T_DISORDS_UTDN, &
                T_UTUP_EIGEN, T_UTDN_EIGEN, &
                T_UTDN_MUBAR, &
                K_REAL, K_COMPLEX, &
                SOLA_XPOS, SOLB_XNEG, &
                WUPPER, LCON, MCON, &
                T_WUPPER, T_WLOWER, &
                L_T_UTUP_EIGEN, L_T_UTDN_EIGEN, &
                L_T_DELT_DISORDS, L_T_DISORDS_UTDN, &
                L_T_UTDN_MUBAR, &
                L_SOLA_XPOS, L_SOLB_XNEG, &
                L_WUPPER, NCON, PCON, &
                L_T_WUPPER, L_T_WLOWER, L_UT_T_PARTIC, &
                QATMOSWF_F )

            ELSE
              CALL QUADATMOSWF_LEVEL_DN ( &
                UTA, NLEVEL, NV, NV_PARAMETERS, &
                FLUX_MULTIPLIER, &
                NSTOKES, NSTREAMS, &
                DO_THERMAL_TRANSONLY, &
                QUAD_STREAMS, &
                T_DELT_DISORDS, T_DELT_EIGEN, &
                K_REAL, K_COMPLEX, &
                SOLA_XPOS, SOLB_XNEG, &
                LCON, MCON, &
                T_WLOWER, &
                L_T_DELT_EIGEN, L_T_DELT_DISORDS, &
                L_SOLA_XPOS, L_SOLB_XNEG, &
                L_WLOWER, NCON, PCON, &
                L_T_WLOWER, &
                QATMOSWF_F )

            ENDIF
          ENDDO
        ENDIF

!  Mean Intensity and Flux  output
!  -------------------------------

        IF ( DO_INCLUDE_MVOUTPUT ) THEN

!  Diffuse term integrated output

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

!  nothing further to do if no solar sources

          IF ( .NOT. DO_INCLUDE_DIRECTBEAM ) GO TO 455

!  For the downward direction, add the direct beam contributions

          IF ( WDIR .EQ. DNIDX ) THEN

!  loop over all the output optical depths

           DO UTA = 1, N_USER_LEVELS

!  For the offgrid values

            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
              UT = PARTLAYERS_OUTINDEX(UTA)
              N  = PARTLAYERS_LAYERIDX(UT)

!  For the offgrid values.......
!     .....Only contributions for layers above the PI cutoff
!    L_INITIAL_TRANS is a logarithmic derivative

              IF ( N .LE. LAYER_PIS_CUTOFF(IBEAM) ) THEN
                IF ( NV.LE.N .OR.NV.EQ.0 ) THEN
                 DO Q = 1, NV_PARAMETERS
                  L_TRANS = L_T_UTDN_MUBAR(UT,NV,IBEAM,Q) + &
                 L_INITIAL_TRANS(N,NV,IBEAM,Q) * T_UTDN_MUBAR(UT,IBEAM)
                  L_TRANS = L_TRANS * INITIAL_TRANS(N,IBEAM)
                  DO O1 = 1, NSTOKES

!  @@@ Rob fix (Forgot the Flux factor)
!                  FTRANS = FLUXVEC(O1) * L_TRANS
!                  FTRANS = FLUX_FACTOR * FLUXVEC(O1) * L_TRANS
!mick fix - changed FLUX_FACTOR to FLUX_MULTIPLIER
!                   FTRANS = FLUX_MULTIPLIER * FLUXVEC(O1) * L_TRANS
! @@@ Rob fix 1/31/11, - changed FLUX_MULTIPLIER to FLUX_FACTOR
                   FTRANS = FLUX_FACTOR * FLUXVEC(O1) * L_TRANS

                   L_DIRECT_MEANI = FTRANS / PI4
                   L_DIRECT_FLUX  = FTRANS * LOCAL_CSZA(N,IBEAM)
                   MINT_ATMOSWF_DIRECT(Q,NV,UTA,IBEAM,O1)=L_DIRECT_MEANI
                   FLUX_ATMOSWF_DIRECT(Q,NV,UTA,IBEAM,O1)=L_DIRECT_FLUX
                   MINT_ATMOSWF(Q,NV,UTA,IBEAM,O1,WDIR) = &
                   MINT_ATMOSWF(Q,NV,UTA,IBEAM,O1,WDIR) + L_DIRECT_MEANI
                   FLUX_ATMOSWF(Q,NV,UTA,IBEAM,O1,WDIR) = &
                   FLUX_ATMOSWF(Q,NV,UTA,IBEAM,O1,WDIR) + L_DIRECT_FLUX
                  ENDDO
                 ENDDO
                ENDIF
              ENDIF


!  For the on-grid values
!    L_INITIAL_TRANS is a logarithmic derivative

             ELSE
              N = UTAU_LEVEL_MASK_DN(UTA)
              IF ( N .LE. LAYER_PIS_CUTOFF(IBEAM) ) THEN
               IF ( N.GT.0 ) THEN
                IF ( NV.LE.N .OR.NV.EQ.0 ) THEN
                 DO Q = 1, NV_PARAMETERS
                  L_TRANS = L_T_DELT_MUBAR(N,NV,IBEAM,Q) + &
                   L_INITIAL_TRANS(N,NV,IBEAM,Q) * T_DELT_MUBAR(N,IBEAM)
                  L_TRANS = L_TRANS * INITIAL_TRANS(N,IBEAM)
                  DO O1 = 1, NSTOKES

!  @@@ Rob fix (Forgot the Flux factor)
!                  FTRANS = FLUXVEC(O1) * L_TRANS
!                  FTRANS = FLUX_FACTOR * FLUXVEC(O1) * L_TRANS
!mick fix - changed FLUX_FACTOR to FLUX_MULTIPLIER
!                   FTRANS = FLUX_MULTIPLIER * FLUXVEC(O1) * L_TRANS
! @@@ Rob fix 1/31/11, - changed FLUX_MULTIPLIER to FLUX_FACTOR
                   FTRANS = FLUX_FACTOR * FLUXVEC(O1) * L_TRANS

                   L_DIRECT_MEANI = FTRANS / PI4
                   L_DIRECT_FLUX  = FTRANS * LOCAL_CSZA(N,IBEAM)
                   MINT_ATMOSWF_DIRECT(Q,NV,UTA,IBEAM,O1)=L_DIRECT_MEANI
                   FLUX_ATMOSWF_DIRECT(Q,NV,UTA,IBEAM,O1)=L_DIRECT_FLUX
                   MINT_ATMOSWF(Q,NV,UTA,IBEAM,O1,WDIR) = &
                   MINT_ATMOSWF(Q,NV,UTA,IBEAM,O1,WDIR) + L_DIRECT_MEANI
                   FLUX_ATMOSWF(Q,NV,UTA,IBEAM,O1,WDIR) = &
                   FLUX_ATMOSWF(Q,NV,UTA,IBEAM,O1,WDIR) + L_DIRECT_FLUX
                  ENDDO
                 ENDDO
                ENDIF
               ENDIF
              ENDIF
             ENDIF

!  End UTA loop

           ENDDO

!  Finish downwelling direct contribution

          ENDIF

!  Continuation point for avoiding direct beam calculation

 455      CONTINUE

!  Finish MV output

        ENDIF

!  end direction loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_L_INTEGRATED_OUTPUT

!

      SUBROUTINE QUADATMOSWF_LEVEL_UP ( &
        UTA, NL, NV, NV_PARAMETERS, &
        FLUX_MULTIPLIER, &
        NSTOKES, NSTREAMS, NLAYERS, &
        DO_THERMAL_TRANSONLY, &
        QUAD_STREAMS, &
        T_DELT_DISORDS, T_DELT_EIGEN, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        LCON, MCON, T_WUPPER, &
        BOA_THTONLY_SOURCE, &
        L_T_DELT_EIGEN, L_T_DELT_DISORDS, &
        L_SOLA_XPOS, L_SOLB_XNEG, &
        L_WUPPER, L_WLOWER, &
        NCON, PCON, &
        L_T_WUPPER, L_BOA_THTONLY_SOURCE, &
        QATMOSWF_F )

!  Upwelling weighting function Fourier components at level boundary NL
!  Quadrature angles only

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          NV, NV_PARAMETERS
      INTEGER, INTENT (IN) ::          UTA, NL
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_WUPPER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: BOA_THTONLY_SOURCE &
          ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_EIGEN &
          ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_DISORDS &
          ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_WUPPER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_WLOWER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: NCON &
          ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PCON &
          ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_WUPPER &
          ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_BOA_THTONLY_SOURCE &
          (MAXSTREAMS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (INOUT) :: QATMOSWF_F &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          N, I, I1, Q, O1, K, KO1, K0, K1, K2, LAY
      DOUBLE PRECISION :: SPAR, SHOM_R, SHOM_CR, SHOM
      DOUBLE PRECISION :: HOM1, HOM2, HOM1CR, HOM2CR, HOM3CR
      DOUBLE PRECISION :: NXR, PXR, NXR1, NXR2, PXR1, PXR2
      DOUBLE PRECISION :: LXR, MXR, LXR1, MXR1, LXR2, MXR2
      DOUBLE PRECISION :: LLXR, MLXR, LLXR1, MLXR1, LLXR2, MLXR2
      DOUBLE PRECISION :: TPROP, THELP, L_TPROP, L_THELP, FM

!  homogeneous and particular solution contributions SHOM and SPAR

!  This depends on the level mask - if this is 0 to NLAYERS - 1, then we
!  looking at the perturbation field at the top of these layers. The
!  case where the level mask = NLAYERS is the upwelling perturbed fields
!  at the bottom of the atmosphere (treated separately).

      N  = NL + 1
      FM = FLUX_MULTIPLIER

!  For the lowest level
!  ====================

!  Thermal transmittance only

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

!  For the lowest level, scattering solution
!  -----------------------------------------

      IF ( NL .EQ. NLAYERS ) THEN

!  If this is also the layer that is varying, extra contributions

        IF ( NV .EQ. NL .OR. NV.EQ. 0 ) THEN

!  Stokes and streams loops

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES

!  Parameter loop

              DO Q = 1, NV_PARAMETERS

!  Real homogeneous solutions

               SHOM_R = ZERO
               DO K = 1, K_REAL(NL)

                LXR  = LCON(K,NL)   *   SOLA_XPOS(I1,O1,K,NL)
                LLXR = LCON(K,NL)   * L_SOLA_XPOS(I1,O1,K,NL,Q)
                MLXR = MCON(K,NL)   * L_SOLB_XNEG(I1,O1,K,NL,Q)
                NXR  = NCON(K,NL,Q) *   SOLA_XPOS(I1,O1,K,NL)
                PXR  = PCON(K,NL,Q) *   SOLB_XNEG(I1,O1,K,NL)

                HOM1 = ( LLXR + NXR )  *   T_DELT_EIGEN(K,NL) &
                             +  LXR    * L_T_DELT_EIGEN(K,NL,Q)
                HOM2 = PXR + MLXR
                SHOM_R = SHOM_R + HOM1 + HOM2
               ENDDO

!  Complex homogeneous solutions

               SHOM_CR = ZERO
               KO1 = K_REAL(NL) + 1
               DO K = 1, K_COMPLEX(NL)
                K0 = 2 * K - 2
                K1 = KO1 + K0
                K2 = K1  + 1

                NXR1  =   NCON(K1,NL,Q) *   SOLA_XPOS(I1,O1,K1,NL) &
                        - NCON(K2,NL,Q) *   SOLA_XPOS(I1,O1,K2,NL)
                NXR2  =   NCON(K1,NL,Q) *   SOLA_XPOS(I1,O1,K2,NL) &
                        + NCON(K2,NL,Q) *   SOLA_XPOS(I1,O1,K1,NL)
                PXR1  =   PCON(K1,NL,Q) *   SOLB_XNEG(I1,O1,K1,NL) &
                        - PCON(K2,NL,Q) *   SOLB_XNEG(I1,O1,K2,NL)

                LXR1  =   LCON(K1,NL) *   SOLA_XPOS(I1,O1,K1,NL) &
                        - LCON(K2,NL) *   SOLA_XPOS(I1,O1,K2,NL)
                LXR2  =   LCON(K1,NL) *   SOLA_XPOS(I1,O1,K2,NL) &
                        + LCON(K2,NL) *   SOLA_XPOS(I1,O1,K1,NL)

                LLXR1  =   LCON(K1,NL) * L_SOLA_XPOS(I1,O1,K1,NL,Q) &
                         - LCON(K2,NL) * L_SOLA_XPOS(I1,O1,K2,NL,Q)
                LLXR2  =   LCON(K1,NL) * L_SOLA_XPOS(I1,O1,K2,NL,Q) &
                         + LCON(K2,NL) * L_SOLA_XPOS(I1,O1,K1,NL,Q)
                MLXR1  =   MCON(K1,NL) * L_SOLB_XNEG(I1,O1,K1,NL,Q) &
                         - MCON(K2,NL) * L_SOLB_XNEG(I1,O1,K2,NL,Q)

                HOM1CR =   ( NXR1 + LLXR1 ) *   T_DELT_EIGEN(K1,NL) &
                         - ( NXR2 + LLXR2 ) *   T_DELT_EIGEN(K2,NL)
                HOM2CR =             LXR1   * L_T_DELT_EIGEN(K1,NL,Q) &
                                   - LXR2   * L_T_DELT_EIGEN(K2,NL,Q)
                HOM3CR = PXR1 + MLXR1
                SHOM_CR = SHOM_CR + HOM1CR + HOM2CR + HOM3CR
               ENDDO

!  real part, add particular solution, complete result

               SHOM = SHOM_R + SHOM_CR
               SPAR = L_WLOWER(I1,O1,NL,Q)
               QATMOSWF_F(Q,UTA,I,O1) = FM * ( SPAR + SHOM )

!  Finish Q, I and O1 loops

              ENDDO
            ENDDO
          ENDDO

!  non-varying lowest layer

        ELSE IF ( NV.LT.NL .AND. NV.NE.0 ) THEN

!  Stokes and streams loops

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES

!  Parameter loop

              DO Q = 1, NV_PARAMETERS

!  Real homogeneous solutions

               SHOM_R = ZERO
               DO K = 1, K_REAL(NL)
                NXR = NCON(K,NL,Q) *   SOLA_XPOS(I1,O1,K,NL)
                PXR = PCON(K,NL,Q) *   SOLB_XNEG(I1,O1,K,NL)
                HOM1 = NXR * T_DELT_EIGEN(K,NL)
                HOM2 = PXR
                SHOM_R = SHOM_R + HOM1 + HOM2
               ENDDO

!  Complex homogeneous solutions

               SHOM_CR = ZERO
               KO1 = K_REAL(NL) + 1
               DO K = 1, K_COMPLEX(NL)
                K0 = 2 * K - 2
                K1 = KO1 + K0
                K2 = K1  + 1
                NXR1  =   NCON(K1,NL,Q) *   SOLA_XPOS(I1,O1,K1,NL) &
                        - NCON(K2,NL,Q) *   SOLA_XPOS(I1,O1,K2,NL)
                NXR2  =   NCON(K1,NL,Q) *   SOLA_XPOS(I1,O1,K2,NL) &
                        + NCON(K2,NL,Q) *   SOLA_XPOS(I1,O1,K1,NL)
                PXR1  =   PCON(K1,NL,Q) *   SOLB_XNEG(I1,O1,K1,NL) &
                        - PCON(K2,NL,Q) *   SOLB_XNEG(I1,O1,K2,NL)
                HOM1CR =   NXR1 * T_DELT_EIGEN(K1,NL) &
                         - NXR2 * T_DELT_EIGEN(K2,NL)
                HOM2CR = PXR1
                SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
               ENDDO

!  real part, add particular solution, complete result

               SHOM = SHOM_R + SHOM_CR
               SPAR = L_WLOWER(I1,O1,NL,Q)
               QATMOSWF_F(Q,UTA,I,O1) = FM * ( SPAR + SHOM )

!  Finish Q, I and O1 loops

              ENDDO
            ENDDO
          ENDDO

        ENDIF

!  End lowest level clause

      ENDIF

!  For other levels in the atmosphere
!  ==================================

!  For other levels, thermal transmittance only

      IF ( NL .NE. NLAYERS .and. DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO Q = 1, NV_PARAMETERS
            L_THELP = L_BOA_THTONLY_SOURCE(I,Q)
            THELP   =   BOA_THTONLY_SOURCE(I,O1)
            DO LAY = NLAYERS, N, -1
              L_THELP = L_THELP * T_DELT_DISORDS(I,LAY)
              IF ( LAY.EQ.NV .OR. NV.EQ.0 ) THEN
                L_TPROP = L_T_WUPPER(I1,LAY,Q) / QUAD_STREAMS(I)
                L_THELP = L_THELP + L_TPROP &
                          + THELP * L_T_DELT_DISORDS(I,LAY,Q)
              ENDIF
              TPROP = T_WUPPER(I1,LAY) / QUAD_STREAMS(I)
              THELP = THELP * T_DELT_DISORDS(I,LAY) + TPROP
            ENDDO
            QATMOSWF_F(Q,UTA,I,O1) = FM * L_THELP
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  For other levels, scattering solution
!  -------------------------------------

      IF ( NL .NE. NLAYERS ) THEN

!  If this is also the layer that is varying, extra contributions

        IF ( NV .EQ. N .OR. NV.EQ. 0 ) THEN

!  Stokes and streams loops

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES

!  Parameter loop

              DO Q = 1, NV_PARAMETERS

!  Real homogeneous solutions

               SHOM_R = ZERO
               DO K = 1, K_REAL(N)
                LXR  = LCON(K,N)   *   SOLA_XPOS(I1,O1,K,N)
                MXR  = MCON(K,N)   *   SOLB_XNEG(I1,O1,K,N)
                LLXR = LCON(K,N)   * L_SOLA_XPOS(I1,O1,K,N,Q)
                MLXR = MCON(K,N)   * L_SOLB_XNEG(I1,O1,K,N,Q)
                NXR  = NCON(K,N,Q) *   SOLA_XPOS(I1,O1,K,N)
                PXR  = PCON(K,N,Q) *   SOLB_XNEG(I1,O1,K,N)
                HOM1 =   NXR + LLXR
                HOM2 = ( PXR + MLXR ) *   T_DELT_EIGEN(K,N) &
                             + MXR    * L_T_DELT_EIGEN(K,N,Q)
                SHOM_R = SHOM_R + HOM1 + HOM2
               ENDDO

!  Complex homogeneous solutions

               SHOM_CR = ZERO
               KO1 = K_REAL(N) + 1
               DO K = 1, K_COMPLEX(N)
                K0 = 2 * K - 2
                K1 = KO1 + K0
                K2 = K1  + 1

                NXR1  =   NCON(K1,N,Q) *   SOLA_XPOS(I1,O1,K1,N) &
                        - NCON(K2,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
                PXR1  =   PCON(K1,N,Q) *   SOLB_XNEG(I1,O1,K1,N) &
                        - PCON(K2,N,Q) *   SOLB_XNEG(I1,O1,K2,N)
                PXR2  =   PCON(K1,N,Q) *   SOLB_XNEG(I1,O1,K2,N) &
                        + PCON(K2,N,Q) *   SOLB_XNEG(I1,O1,K1,N)

                MXR1  =   MCON(K1,N) *   SOLB_XNEG(I1,O1,K1,N) &
                        - MCON(K2,N) *   SOLB_XNEG(I1,O1,K2,N)
                MXR2  =   MCON(K1,N) *   SOLB_XNEG(I1,O1,K2,N) &
                        + MCON(K2,N) *   SOLB_XNEG(I1,O1,K1,N)

                LLXR1  =   LCON(K1,N) * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                         - LCON(K2,N) * L_SOLA_XPOS(I1,O1,K2,N,Q)
                MLXR1  =   MCON(K1,N) * L_SOLB_XNEG(I1,O1,K1,N,Q) &
                         - MCON(K2,N) * L_SOLB_XNEG(I1,O1,K2,N,Q)
                MLXR2  =   MCON(K1,N) * L_SOLB_XNEG(I1,O1,K2,N,Q) &
                         + MCON(K2,N) * L_SOLB_XNEG(I1,O1,K1,N,Q)

                HOM1CR =   ( PXR1 + MLXR1 ) *   T_DELT_EIGEN(K1,N) &
                         - ( PXR2 + MLXR2 ) *   T_DELT_EIGEN(K2,N)
                HOM2CR =             MXR1   * L_T_DELT_EIGEN(K1,N,Q) &
                                   - MXR2   * L_T_DELT_EIGEN(K2,N,Q)
                HOM3CR = NXR1 + LLXR1
                SHOM_CR = SHOM_CR + HOM1CR + HOM2CR + HOM3CR

               ENDDO

!  real part, add particular solution, complete result

               SHOM = SHOM_R + SHOM_CR
               SPAR = L_WUPPER(I1,O1,N,Q)
               QATMOSWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER * (SPAR+SHOM)

!  Finish Q, I and O1 loops

              ENDDO
            ENDDO
          ENDDO

!  For layers N beneath or above the varying layer NV

        ELSE IF ( NV.NE.N .AND. NV.NE.0 ) THEN

!  Stokes and streams loops

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES

!  Parameter loop

              DO Q = 1, NV_PARAMETERS

!  Real homogeneous solutions

               SHOM_R = ZERO
               DO K = 1, K_REAL(N)
                NXR  = NCON(K,N,Q) *   SOLA_XPOS(I1,O1,K,N)
                PXR  = PCON(K,N,Q) *   SOLB_XNEG(I1,O1,K,N)
                HOM1 = NXR
                HOM2 = PXR * T_DELT_EIGEN(K,N)
                SHOM_R = SHOM_R + HOM1 + HOM2
               ENDDO

!  Complex homogeneous solutions

               SHOM_CR = ZERO
               KO1 = K_REAL(N) + 1
               DO K = 1, K_COMPLEX(N)
                K0 = 2 * K - 2
                K1 = KO1 + K0
                K2 = K1  + 1
                NXR1  =   NCON(K1,N,Q) *   SOLA_XPOS(I1,O1,K1,N) &
                        - NCON(K2,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
                PXR1  =   PCON(K1,N,Q) *   SOLB_XNEG(I1,O1,K1,N) &
                        - PCON(K2,N,Q) *   SOLB_XNEG(I1,O1,K2,N)
                PXR2  =   PCON(K1,N,Q) *   SOLB_XNEG(I1,O1,K2,N) &
                        + PCON(K2,N,Q) *   SOLB_XNEG(I1,O1,K1,N)
                HOM1CR =   PXR1 * T_DELT_EIGEN(K1,N) &
                         - PXR2 * T_DELT_EIGEN(K2,N)
                HOM2CR = NXR1
                SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
               ENDDO

!  real part, add particular solution (only if N>NV), complete result

               SHOM = SHOM_R + SHOM_CR
               SPAR = ZERO
               IF (NV.LT.N)SPAR = L_WUPPER(I1,O1,N,Q)
               QATMOSWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER * (SPAR+SHOM)

!  Finish Q, I and O1 loops

              ENDDO
            ENDDO
          ENDDO

!  End variability clauses

        ENDIF
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE QUADATMOSWF_LEVEL_UP

!

      SUBROUTINE QUADATMOSWF_LEVEL_DN ( &
        UTA, NL, NV, NV_PARAMETERS, &
        FLUX_MULTIPLIER, &
        NSTOKES, NSTREAMS, &
        DO_THERMAL_TRANSONLY, &
        QUAD_STREAMS, &
        T_DELT_DISORDS, T_DELT_EIGEN, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        LCON, MCON, &
        T_WLOWER, &
        L_T_DELT_EIGEN, L_T_DELT_DISORDS, &
        L_SOLA_XPOS, L_SOLB_XNEG, &
        L_WLOWER, NCON, PCON, &
        L_T_WLOWER, &
        QATMOSWF_F )

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          NV, NV_PARAMETERS
      INTEGER, INTENT (IN) ::          UTA, NL
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_WLOWER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_EIGEN &
          ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_DISORDS &
          ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_WLOWER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: NCON &
          ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PCON &
          ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_WLOWER &
          ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS)

      DOUBLE PRECISION, INTENT (INOUT) :: QATMOSWF_F &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          N, I, Q, O1, K, KO1, K0, K1, K2, LAY
      DOUBLE PRECISION :: SPAR, SHOM_R, SHOM_CR, SHOM
      DOUBLE PRECISION :: HOM1, HOM2, HOM1CR, HOM2CR, HOM3CR
      DOUBLE PRECISION :: NXR, PXR, NXR1, NXR2, PXR1
      DOUBLE PRECISION :: LXR, LXR1, LXR2
      DOUBLE PRECISION :: LLXR, MLXR, LLXR1, MLXR1, LLXR2
      DOUBLE PRECISION :: TPROP, THELP, L_TPROP, L_THELP

!  Downwelling weighting function at TOA ( or N = 0 ) is zero
!    Zero and return

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

!  For other levels in the atmosphere
!  ----------------------------------

!  Other levels, Thermal transmittance-only solution
!  Thermal transmittance solution, build from TOA downwards
!  Scattering solution, use the Discrete Ordinate solution

      IF ( NL.NE.0 .and. DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          DO Q = 1, NV_PARAMETERS
            L_THELP = ZERO
            THELP   = ZERO
            DO LAY = 1, NL
              L_THELP = L_THELP * T_DELT_DISORDS(I,LAY)
              IF ( LAY.EQ.NV .OR. NV.EQ.0 ) THEN
                L_TPROP = L_T_WLOWER(I,LAY,Q) / QUAD_STREAMS(I)
                L_THELP = L_THELP +  L_TPROP &
                          + THELP * L_T_DELT_DISORDS(I,LAY,Q)
              ENDIF
              TPROP = T_WLOWER(I,LAY) / QUAD_STREAMS(I)
              THELP = THELP * T_DELT_DISORDS(I,LAY) + TPROP
            ENDDO
            QATMOSWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER * L_THELP
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Scattering solutions
!  --------------------

      IF ( NL.NE.0 ) THEN

!  Shorthand

        N = NL

!  If this is also the layer that is varying, extra contributions

        IF ( NV .EQ. N .OR. NV.EQ.0 ) THEN

!  Stokes and streams loops

          DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES

!  Parameter loop

              DO Q = 1, NV_PARAMETERS

!  Real homogeneous solutions

               SHOM_R = ZERO
               DO K = 1, K_REAL(N)
                LXR  = LCON(K,N)   *   SOLA_XPOS(I,O1,K,N)
                LLXR = LCON(K,N)   * L_SOLA_XPOS(I,O1,K,N,Q)
                MLXR = MCON(K,N)   * L_SOLB_XNEG(I,O1,K,N,Q)
                NXR  = NCON(K,N,Q) *   SOLA_XPOS(I,O1,K,N)
                PXR  = PCON(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
                HOM1 = ( NXR + LLXR ) *   T_DELT_EIGEN(K,N) &
                             +  LXR   * L_T_DELT_EIGEN(K,N,Q)
                HOM2 = PXR + MLXR
                SHOM_R = SHOM_R + HOM1 + HOM2
               ENDDO

!  Complex homogeneous solutions

               SHOM_CR = ZERO
               KO1 = K_REAL(N) + 1
               DO K = 1, K_COMPLEX(N)
                K0 = 2 * K - 2
                K1 = KO1 + K0
                K2 = K1  + 1

                NXR1  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K1,N) &
                        - NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K2,N)
                NXR2  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K2,N) &
                        + NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K1,N)
                PXR1  =   PCON(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N) &
                        - PCON(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)

                LXR1  =   LCON(K1,N) *   SOLA_XPOS(I,O1,K1,N) &
                        - LCON(K2,N) *   SOLA_XPOS(I,O1,K2,N)
                LXR2  =   LCON(K1,N) *   SOLA_XPOS(I,O1,K2,N) &
                        + LCON(K2,N) *   SOLA_XPOS(I,O1,K1,N)

                LLXR1  =   LCON(K1,N) * L_SOLA_XPOS(I,O1,K1,N,Q) &
                         - LCON(K2,N) * L_SOLA_XPOS(I,O1,K2,N,Q)
                LLXR2  =   LCON(K1,N) * L_SOLA_XPOS(I,O1,K2,N,Q) &
                         + LCON(K2,N) * L_SOLA_XPOS(I,O1,K1,N,Q)
                MLXR1  =   MCON(K1,N) * L_SOLB_XNEG(I,O1,K1,N,Q) &
                         - MCON(K2,N) * L_SOLB_XNEG(I,O1,K2,N,Q)

                HOM1CR =   ( NXR1 + LLXR1 ) *   T_DELT_EIGEN(K1,N) &
                         - ( NXR2 + LLXR2 ) *   T_DELT_EIGEN(K2,N)
                HOM2CR =             LXR1   * L_T_DELT_EIGEN(K1,N,Q) &
                                   - LXR2   * L_T_DELT_EIGEN(K2,N,Q)
                HOM3CR = PXR1 + MLXR1
                SHOM_CR = SHOM_CR + HOM1CR + HOM2CR + HOM3CR

               ENDDO

!  real part, add particular solution, complete result

               SHOM = SHOM_R + SHOM_CR
               SPAR = L_WLOWER(I,O1,N,Q)
               QATMOSWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER * (SPAR+SHOM)

!  Finish Q, I and O1 loops

              ENDDO
            ENDDO
          ENDDO

!  varying layer above or below active layer

        ELSE IF ( NV.NE.N .AND. NV.NE.0 ) THEN

!  Stokes and streams loops

          DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES

!  Parameter loop

              DO Q = 1, NV_PARAMETERS

!  Real homogeneous solutions

               SHOM_R = ZERO
               DO K = 1, K_REAL(N)
                NXR = NCON(K,N,Q) *   SOLA_XPOS(I,O1,K,N)
                PXR = PCON(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
                HOM1 = NXR * T_DELT_EIGEN(K,N)
                HOM2 = PXR
                SHOM_R = SHOM_R + HOM1 + HOM2
               ENDDO

!  Complex homogeneous solutions

               SHOM_CR = ZERO
               KO1 = K_REAL(N) + 1
               DO K = 1, K_COMPLEX(N)
                K0 = 2 * K - 2
                K1 = KO1 + K0
                K2 = K1  + 1
                NXR1  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K1,N) &
                        - NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K2,N)
                NXR2  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K2,N) &
                        + NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K1,N)
                PXR1  =   PCON(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N) &
                        - PCON(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)
                HOM1CR =   NXR1 * T_DELT_EIGEN(K1,N) &
                         - NXR2 * T_DELT_EIGEN(K2,N)
                HOM2CR = PXR1
                SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
               ENDDO

!  real part, add particular solution (only if N > NV), complete result

               SHOM = SHOM_R + SHOM_CR
               SPAR = ZERO
               IF ( N .GT. NV ) SPAR = L_WLOWER(I,O1,N,Q)
               QATMOSWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER * (SPAR+SHOM)

!  Finish Q, I and O1 loops

              ENDDO
            ENDDO
          ENDDO

!  Finish variability clauses

        ENDIF
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE QUADATMOSWF_LEVEL_DN

!

      SUBROUTINE QUADATMOSWF_OFFGRID_UP ( &
        IB, UTA, UT, N, NV, NV_PARAMETERS, &
        DO_INCLUDE_THERMEMISS, FLUX_MULTIPLIER, &
        DO_SOLAR_SOURCES, NSTOKES, &
        NSTREAMS, NLAYERS, &
        DO_THERMAL_TRANSONLY, &
        DO_CLASSICAL_SOLUTION, QUAD_STREAMS, &
        T_DELT_DISORDS, T_DISORDS_UTUP, &
        T_UTUP_EIGEN, T_UTDN_EIGEN, &
        T_UTDN_MUBAR, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        WUPPER, LCON, MCON, &
        T_WUPPER, &
        BOA_THTONLY_SOURCE, &
        L_T_UTUP_EIGEN, L_T_UTDN_EIGEN, &
        L_T_DELT_DISORDS, L_T_DISORDS_UTUP, &
        L_T_UTDN_MUBAR, &
        L_SOLA_XPOS, L_SOLB_XNEG, &
        L_WUPPER, NCON, PCON, &
        L_T_WUPPER, L_UT_T_PARTIC, &
        L_BOA_THTONLY_SOURCE, &
        QATMOSWF_F )

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          NV, NV_PARAMETERS
      INTEGER, INTENT (IN) ::          IB, UTA, UT, N
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTUP &
          ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN &
          ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN &
          ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: WUPPER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_WUPPER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: BOA_THTONLY_SOURCE &
          ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTUP_EIGEN &
          ( MAXEVALUES, MAX_USER_LEVELS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_EIGEN &
          ( MAXEVALUES, MAX_USER_LEVELS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_DISORDS &
          ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DISORDS_UTUP &
          ( MAXSTREAMS, MAX_USER_LEVELS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_MUBAR &
          ( MAX_USER_LEVELS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)
      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_WUPPER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: NCON &
          ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PCON &
          ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_WUPPER &
          ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UT_T_PARTIC &
          ( MAXSTREAMS_2, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_BOA_THTONLY_SOURCE &
          (MAXSTREAMS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (INOUT) :: QATMOSWF_F &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  local variables
!  ---------------

!  @@@@@@@ Rob fix, 2/9/11, added variables VAREX, SUNP, L_SUNP

      LOGICAL ::          VAREX
      INTEGER ::          I, I1, Q, O1, K, KO1, K0, K1, K2, LAY
      DOUBLE PRECISION :: SPAR, SHOM_R, SHOM_CR, SHOM, SUNP, L_SUNP
      DOUBLE PRECISION :: HOM1, HOM2, HOM1CR, HOM2CR, HOM3CR, HOM4CR
      DOUBLE PRECISION :: NXR, PXR, NXR1, NXR2, PXR1, PXR2
      DOUBLE PRECISION :: LXR, MXR, LXR1, MXR1, LXR2, MXR2
      DOUBLE PRECISION :: LLXR, MLXR, LLXR1, MLXR1, LLXR2, MLXR2
      DOUBLE PRECISION :: TPROP, THELP, L_TPROP, L_THELP, FMULT

!  short hand

      FMULT = FLUX_MULTIPLIER

!  Thermal Transmittance only
!  --------------------------

      IF ( DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO Q = 1, NV_PARAMETERS
            THELP     =   BOA_THTONLY_SOURCE(I,O1)
            L_THELP   = L_BOA_THTONLY_SOURCE(I,Q)
            DO LAY = NLAYERS, N+1, -1
              L_THELP = L_THELP *   T_DELT_DISORDS(I,LAY)
              IF ( LAY.EQ.NV .OR. NV.EQ.0 ) THEN
                L_TPROP = L_T_WUPPER(I1,LAY,Q) / QUAD_STREAMS(I)
                L_THELP = L_THELP + L_TPROP &
                          + THELP * L_T_DELT_DISORDS(I,LAY,Q)
              ENDIF
              TPROP = T_WUPPER(I1,Q) / QUAD_STREAMS(I)
              THELP = THELP * T_DELT_DISORDS(I,LAY) + TPROP
            ENDDO
            L_THELP = L_THELP * T_DISORDS_UTUP(I,UT)
            IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN
              L_TPROP = L_UT_T_PARTIC(I1,UT,Q) / QUAD_STREAMS(I)
              L_THELP = L_THELP + L_TPROP &
                        + THELP * L_T_DISORDS_UTUP(I,UT,Q)
            ENDIF
            QATMOSWF_F(Q,UTA,I,O1) = FMULT * L_THELP
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  For those optical depths at off-grid levels
!  ###########################################

!  Homogeneous
!  -----------

!  solution for N being the varying layer NV

      IF ( N.EQ. NV .OR. NV.EQ.0 ) THEN

!  stream directions

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES

!  Parameter loop

            DO Q = 1, NV_PARAMETERS

!  real homogeneous solutions

              SHOM_R = ZERO
              DO K = 1, K_REAL(N)
                LXR  = LCON(K,N)   *   SOLA_XPOS(I1,O1,K,N)
                MXR  = MCON(K,N)   *   SOLB_XNEG(I1,O1,K,N)
                LLXR = LCON(K,N)   * L_SOLA_XPOS(I1,O1,K,N,Q)
                MLXR = MCON(K,N)   * L_SOLB_XNEG(I1,O1,K,N,Q)
                NXR  = NCON(K,N,Q) *   SOLA_XPOS(I1,O1,K,N)
                PXR  = PCON(K,N,Q) *   SOLB_XNEG(I1,O1,K,N)
                HOM1 = ( NXR + LLXR ) *   T_UTDN_EIGEN(K,UT) &
                             +  LXR   * L_T_UTDN_EIGEN(K,UT,Q)
                HOM2 = ( PXR + MLXR ) *   T_UTUP_EIGEN(K,UT) &
                             +  MXR   * L_T_UTUP_EIGEN(K,UT,Q)
                SHOM_R = SHOM_R + HOM1 + HOM2
              ENDDO

!  complex homogeneous solutions

              SHOM_CR = ZERO
              KO1 = K_REAL(N) + 1
              DO K = 1, K_COMPLEX(N)
                K0 = 2 * K - 2
                K1 = KO1 + K0
                K2 = K1  + 1

                NXR1  =   NCON(K1,N,Q) *   SOLA_XPOS(I1,O1,K1,N) &
                        - NCON(K2,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
                NXR2  =   NCON(K1,N,Q) *   SOLA_XPOS(I1,O1,K2,N) &
                        + NCON(K2,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
                PXR1  =   PCON(K1,N,Q) *   SOLB_XNEG(I1,O1,K1,N) &
                        - PCON(K2,N,Q) *   SOLB_XNEG(I1,O1,K2,N)
                PXR2  =   PCON(K1,N,Q) *   SOLB_XNEG(I1,O1,K2,N) &
                        + PCON(K2,N,Q) *   SOLB_XNEG(I1,O1,K1,N)

                LXR1  =   LCON(K1,N) *   SOLA_XPOS(I1,O1,K1,N) &
                        - LCON(K2,N) *   SOLA_XPOS(I1,O1,K2,N)
                LXR2  =   LCON(K1,N) *   SOLA_XPOS(I1,O1,K2,N) &
                        + LCON(K2,N) *   SOLA_XPOS(I1,O1,K1,N)
                MXR1  =   MCON(K1,N) *   SOLB_XNEG(I1,O1,K1,N) &
                        - MCON(K2,N) *   SOLB_XNEG(I1,O1,K2,N)
                MXR2  =   MCON(K1,N) *   SOLB_XNEG(I1,O1,K2,N) &
                        + MCON(K2,N) *   SOLB_XNEG(I1,O1,K1,N)

                LLXR1  =   LCON(K1,N) * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                         - LCON(K2,N) * L_SOLA_XPOS(I1,O1,K2,N,Q)
                LLXR2  =   LCON(K1,N) * L_SOLA_XPOS(I1,O1,K2,N,Q) &
                         + LCON(K2,N) * L_SOLA_XPOS(I1,O1,K1,N,Q)
                MLXR1  =   MCON(K1,N) * L_SOLB_XNEG(I1,O1,K1,N,Q) &
                         - MCON(K2,N) * L_SOLB_XNEG(I1,O1,K2,N,Q)
                MLXR2  =   MCON(K1,N) * L_SOLB_XNEG(I1,O1,K2,N,Q) &
                         + MCON(K2,N) * L_SOLB_XNEG(I1,O1,K1,N,Q)

                HOM1CR =   ( NXR1 + LLXR1 ) *   T_UTDN_EIGEN(K1,UT) &
                         - ( NXR2 + LLXR2 ) *   T_UTDN_EIGEN(K2,UT)
                HOM2CR =             LXR1   * L_T_UTDN_EIGEN(K1,UT,Q) &
                                   - LXR2   * L_T_UTDN_EIGEN(K2,UT,Q)
                HOM3CR =   ( PXR1 + MLXR1 ) *   T_UTUP_EIGEN(K1,UT) &
                         - ( PXR2 + MLXR2 ) *   T_UTUP_EIGEN(K2,UT)
                HOM4CR =             MXR1   * L_T_UTUP_EIGEN(K1,UT,Q) &
                                   - MXR2   * L_T_UTUP_EIGEN(K2,UT,Q)

                SHOM_CR = SHOM_CR + HOM1CR + HOM2CR + HOM3CR + HOM4CR

              ENDDO

!  real part

              SHOM = SHOM_R + SHOM_CR
              QATMOSWF_F(Q,UTA,I,O1) = FMULT * SHOM

!  Finish Q, I and O1 loops

            ENDDO
          ENDDO
        ENDDO

!  Solution for N  below/above the layer NV that varies

      ELSE IF ( NV.NE.N .AND. NV.NE.0 ) THEN

!  stream/stokes directions

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES

!  Parameter loop

            DO Q = 1, NV_PARAMETERS

!  real homogeneous solutions

              SHOM_R = ZERO
              DO K = 1, K_REAL(N)
                NXR  = NCON(K,N,Q) *   SOLA_XPOS(I1,O1,K,N)
                PXR  = PCON(K,N,Q) *   SOLB_XNEG(I1,O1,K,N)
                HOM1 = NXR * T_UTDN_EIGEN(K,UT)
                HOM2 = PXR * T_UTUP_EIGEN(K,UT)
                SHOM_R = SHOM_R + HOM1 + HOM2
              ENDDO

!  complex homogeneous solutions

              SHOM_CR = ZERO
              KO1 = K_REAL(N) + 1
              DO K = 1, K_COMPLEX(N)
                K0 = 2 * K - 2
                K1 = KO1 + K0
                K2 = K1  + 1
                NXR1  =   NCON(K1,N,Q) *   SOLA_XPOS(I1,O1,K1,N) &
                        - NCON(K2,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
                NXR2  =   NCON(K1,N,Q) *   SOLA_XPOS(I1,O1,K2,N) &
                        + NCON(K2,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
                PXR1  =   PCON(K1,N,Q) *   SOLB_XNEG(I1,O1,K1,N) &
                        - PCON(K2,N,Q) *   SOLB_XNEG(I1,O1,K2,N)
                PXR2  =   PCON(K1,N,Q) *   SOLB_XNEG(I1,O1,K2,N) &
                        + PCON(K2,N,Q) *   SOLB_XNEG(I1,O1,K1,N)
                HOM1CR =   NXR1  * T_UTDN_EIGEN(K1,UT) &
                         - NXR2  * T_UTDN_EIGEN(K2,UT)
                HOM2CR =   PXR1  * T_UTUP_EIGEN(K1,UT) &
                         - PXR2  * T_UTUP_EIGEN(K2,UT)
                SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
              ENDDO

!  real part

              SHOM = SHOM_R + SHOM_CR
              QATMOSWF_F(Q,UTA,I,O1) = FMULT * SHOM

!  Finish Q, I and O1 loops

            ENDDO
          ENDDO
        ENDDO

!  end variability clause

      ENDIF

!  Add the linearized thermal solution  (if flagged)
!    ---Only present if N = NV, or NV = 0 (column linearization)
!   THIS is the solution with scattering

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        O1 = 1
        IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO Q = 1, NV_PARAMETERS
              SPAR = L_UT_T_PARTIC(I1,UT,Q)
              QATMOSWF_F(Q,UTA,I,O1) = &
                  QATMOSWF_F(Q,UTA,I,O1) + FMULT * SPAR
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  Finished if no solar terms

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Add Linearized Classical beam solution
!  Green function solution is still a Placeholder....!
!    solutions exist only for N = NV or N > NV or NV = 0

! @@@@@@@@@@ Rob fix, 2/9/11, Replacement of commented section
!      IF ( DO_CLASSICAL_SOLUTION ) THEN
!        IF ( N.GE.NV .OR. NV.EQ.0 ) THEN
!          DO I = 1, NSTREAMS
!            I1 = I + NSTREAMS
!            DO O1 = 1, NSTOKES
!              DO Q = 1, NV_PARAMETERS
!                SPAR = L_WUPPER(I1,O1,N,Q) *   T_UTDN_MUBAR(UT,IB) + &
!                         WUPPER(I1,O1,N)   * L_T_UTDN_MUBAR(UT,NV,IB,Q)
!                QATMOSWF_F(Q,UTA,I,O1) = &
!                  QATMOSWF_F(Q,UTA,I,O1) + FMULT * SPAR
!              ENDDO
!            ENDDO
!          ENDDO
!        ENDIF
!      ELSE
!                P L A C E H O L D E R
!      ENDIF

      IF ( DO_CLASSICAL_SOLUTION ) THEN
        IF ( DO_INCLUDE_THERMEMISS ) THEN
          VAREX = ( N.EQ.NV .OR. NV.EQ.0 )
          IF ( N.GT.NV .or. VAREX ) THEN
            DO I = 1, NSTREAMS
              I1 = I + NSTREAMS
              DO O1 = 1, NSTOKES
                DO Q = 1, NV_PARAMETERS
                  SUNP   = WUPPER(I1,O1,N)
                  L_SUNP = L_WUPPER(I1,O1,N,Q)
                  IF ( O1 .eq. 1 ) THEN
                    SUNP   = SUNP - T_WUPPER(I1,N)
                    IF (VAREX) L_SUNP = L_SUNP - L_T_WUPPER(I1,N,Q)
                  ENDIF
                  SPAR = L_SUNP *   T_UTDN_MUBAR(UT,IB) + &
                           SUNP * L_T_UTDN_MUBAR(UT,NV,IB,Q)
                  QATMOSWF_F(Q,UTA,I,O1) = &
                  QATMOSWF_F(Q,UTA,I,O1) + FMULT * SPAR
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ELSE
          IF ( N.GE.NV .OR. NV.EQ.0 ) THEN
            DO I = 1, NSTREAMS
              I1 = I + NSTREAMS
              DO O1 = 1, NSTOKES
                DO Q = 1, NV_PARAMETERS
                  SPAR = L_WUPPER(I1,O1,N,Q) *   T_UTDN_MUBAR(UT,IB) + &
                           WUPPER(I1,O1,N)   * L_T_UTDN_MUBAR(UT,NV,IB,Q)
                  QATMOSWF_F(Q,UTA,I,O1) = &
                    QATMOSWF_F(Q,UTA,I,O1) + FMULT * SPAR
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDIF
      ENDIF
! @@@@@@@@@@ Rob fix, 2/9/11, End Replacement section

!  Finish

      RETURN
      END SUBROUTINE QUADATMOSWF_OFFGRID_UP

!

      SUBROUTINE QUADATMOSWF_OFFGRID_DN ( &
        IB, UTA, UT, N, NV, NV_PARAMETERS, &
        DO_INCLUDE_THERMEMISS, FLUX_MULTIPLIER, &
        DO_SOLAR_SOURCES, NSTOKES, &
        NSTREAMS, DO_THERMAL_TRANSONLY, &
        DO_CLASSICAL_SOLUTION, QUAD_STREAMS, &
        T_DELT_DISORDS, T_DISORDS_UTDN, &
        T_UTUP_EIGEN, T_UTDN_EIGEN, &
        T_UTDN_MUBAR, &
        K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, &
        WUPPER, LCON, MCON, &
        T_WUPPER, T_WLOWER, &
        L_T_UTUP_EIGEN, L_T_UTDN_EIGEN, &
        L_T_DELT_DISORDS, L_T_DISORDS_UTDN, &
        L_T_UTDN_MUBAR, &
        L_SOLA_XPOS, L_SOLB_XNEG, &
        L_WUPPER, NCON, PCON, &
        L_T_WUPPER, L_T_WLOWER, L_UT_T_PARTIC, &
        QATMOSWF_F )

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          NV, NV_PARAMETERS
      INTEGER, INTENT (IN) ::          IB, UTA, UT, N
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTDN &
          ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN &
          ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN &
          ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: WUPPER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_WUPPER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_WLOWER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTUP_EIGEN &
          ( MAXEVALUES, MAX_USER_LEVELS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_EIGEN &
          ( MAXEVALUES, MAX_USER_LEVELS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_DISORDS &
          ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DISORDS_UTDN &
          ( MAXSTREAMS, MAX_USER_LEVELS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_MUBAR &
          ( MAX_USER_LEVELS, 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)
      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_WUPPER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: NCON &
          ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PCON &
          ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_WUPPER &
          ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_WLOWER &
          ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UT_T_PARTIC &
          ( MAXSTREAMS_2, MAX_PARTLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (INOUT) :: QATMOSWF_F &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  local variables
!  ---------------

!  @@@@@@@ Rob fix, 2/9/11, added variables VAREX, SUNP, L_SUNP

      LOGICAL ::          VAREX
      INTEGER ::          I, Q, O1, K, KO1, K0, K1, K2, LAY
      DOUBLE PRECISION :: SPAR, SHOM_R, SHOM_CR, SHOM, SUNP, L_SUNP
      DOUBLE PRECISION :: HOM1, HOM2, HOM1CR, HOM2CR, HOM3CR, HOM4CR
      DOUBLE PRECISION :: NXR, PXR, NXR1, NXR2, PXR1, PXR2
      DOUBLE PRECISION :: LXR, MXR, LXR1, MXR1, LXR2, MXR2
      DOUBLE PRECISION :: LLXR, MLXR, LLXR1, MLXR1, LLXR2, MLXR2
      DOUBLE PRECISION :: TPROP, THELP, L_TPROP, L_THELP, FMULT

!  Short hand

      FMULT = FLUX_MULTIPLIER

!  Thermal Transmittance only
!  --------------------------

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
                L_THELP = L_THELP + L_TPROP &
                          + THELP * L_T_DELT_DISORDS(I,LAY,Q)
              ENDIF
              TPROP = T_WLOWER(I,LAY) / QUAD_STREAMS(I)
              THELP = THELP * T_DELT_DISORDS(I,LAY) + TPROP
            ENDDO
            L_THELP = L_THELP * T_DISORDS_UTDN(I,UT)
            IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN
              L_TPROP = L_UT_T_PARTIC(I,UT,Q) / QUAD_STREAMS(I)
              L_THELP = L_THELP + L_TPROP &
                        + THELP * L_T_DISORDS_UTDN(I,UT,Q)
            ENDIF
            QATMOSWF_F(Q,UTA,I,O1) = FMULT * L_THELP
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  For those optical depths at off-grid levels
!  ###########################################

!  Homogeneous
!  -----------

!  solution for N being the varying layer NV

      IF ( N .EQ. NV .OR. NV.EQ.0 ) THEN

!  stream and Stokes directions

        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES

!  Parameter loop

            DO Q = 1, NV_PARAMETERS

!  real homogeneous solutions

              SHOM_R = ZERO
              DO K = 1, K_REAL(N)
                LXR  = LCON(K,N)   *   SOLA_XPOS(I,O1,K,N)
                MXR  = MCON(K,N)   *   SOLB_XNEG(I,O1,K,N)
                LLXR = LCON(K,N)   * L_SOLA_XPOS(I,O1,K,N,Q)
                MLXR = MCON(K,N)   * L_SOLB_XNEG(I,O1,K,N,Q)
                NXR  = NCON(K,N,Q) *   SOLA_XPOS(I,O1,K,N)
                PXR  = PCON(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
                HOM1 = ( NXR + LLXR ) *   T_UTDN_EIGEN(K,UT) &
                             +  LXR   * L_T_UTDN_EIGEN(K,UT,Q)
                HOM2 = ( PXR + MLXR ) *   T_UTUP_EIGEN(K,UT) &
                             +  MXR   * L_T_UTUP_EIGEN(K,UT,Q)
                SHOM_R = SHOM_R + HOM1 + HOM2
              ENDDO

!  complex homogeneous solutions

              SHOM_CR = ZERO
              KO1 = K_REAL(N) + 1
              DO K = 1, K_COMPLEX(N)
                K0 = 2 * K - 2
                K1 = KO1 + K0
                K2 = K1  + 1

                NXR1  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K1,N) &
                        - NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K2,N)
                NXR2  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K2,N) &
                        + NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K1,N)
                PXR1  =   PCON(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N) &
                        - PCON(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)
                PXR2  =   PCON(K1,N,Q) *   SOLB_XNEG(I,O1,K2,N) &
                        + PCON(K2,N,Q) *   SOLB_XNEG(I,O1,K1,N)

                LXR1  =   LCON(K1,N) *   SOLA_XPOS(I,O1,K1,N) &
                        - LCON(K2,N) *   SOLA_XPOS(I,O1,K2,N)
                LXR2  =   LCON(K1,N) *   SOLA_XPOS(I,O1,K2,N) &
                        + LCON(K2,N) *   SOLA_XPOS(I,O1,K1,N)
                MXR1  =   MCON(K1,N) *   SOLB_XNEG(I,O1,K1,N) &
                        - MCON(K2,N) *   SOLB_XNEG(I,O1,K2,N)
                MXR2  =   MCON(K1,N) *   SOLB_XNEG(I,O1,K2,N) &
                        + MCON(K2,N) *   SOLB_XNEG(I,O1,K1,N)

                LLXR1  =   LCON(K1,N) * L_SOLA_XPOS(I,O1,K1,N,Q) &
                         - LCON(K2,N) * L_SOLA_XPOS(I,O1,K2,N,Q)
                LLXR2  =   LCON(K1,N) * L_SOLA_XPOS(I,O1,K2,N,Q) &
                         + LCON(K2,N) * L_SOLA_XPOS(I,O1,K1,N,Q)
                MLXR1  =   MCON(K1,N) * L_SOLB_XNEG(I,O1,K1,N,Q) &
                         - MCON(K2,N) * L_SOLB_XNEG(I,O1,K2,N,Q)
                MLXR2  =   MCON(K1,N) * L_SOLB_XNEG(I,O1,K2,N,Q) &
                         + MCON(K2,N) * L_SOLB_XNEG(I,O1,K1,N,Q)

                HOM1CR =   ( NXR1 + LLXR1 ) *   T_UTDN_EIGEN(K1,UT) &
                         - ( NXR2 + LLXR2 ) *   T_UTDN_EIGEN(K2,UT)
                HOM2CR =             LXR1   * L_T_UTDN_EIGEN(K1,UT,Q) &
                                   - LXR2   * L_T_UTDN_EIGEN(K2,UT,Q)
                HOM3CR =   ( PXR1 + MLXR1 ) *   T_UTUP_EIGEN(K1,UT) &
                         - ( PXR2 + MLXR2 ) *   T_UTUP_EIGEN(K2,UT)
                HOM4CR =             MXR1   * L_T_UTUP_EIGEN(K1,UT,Q) &
                                   - MXR2   * L_T_UTUP_EIGEN(K2,UT,Q)

                SHOM_CR = SHOM_CR + HOM1CR + HOM2CR + HOM3CR + HOM4CR
              ENDDO

!  real part

              SHOM = SHOM_R + SHOM_CR
              QATMOSWF_F(Q,UTA,I,O1) = FMULT * SHOM

!  Finish Q, I and O1 loops

            ENDDO
          ENDDO
        ENDDO

!  Solution for N  below/above the layer NV that varies

      ELSE IF ( NV.NE.N .AND. NV.NE.0 ) THEN

!  stream/stokes directions

        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES

!  Parameter loop

            DO Q = 1, NV_PARAMETERS

!  real homogeneous solutions

              SHOM_R = ZERO
              DO K = 1, K_REAL(N)
                NXR  = NCON(K,N,Q) *   SOLA_XPOS(I,O1,K,N)
                PXR  = PCON(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
                HOM1 = NXR * T_UTDN_EIGEN(K,UT)
                HOM2 = PXR * T_UTUP_EIGEN(K,UT)
                SHOM_R = SHOM_R + HOM1 + HOM2
              ENDDO

!  complex homogeneous solutions

              SHOM_CR = ZERO
              KO1 = K_REAL(N) + 1
              DO K = 1, K_COMPLEX(N)
                K0 = 2 * K - 2
                K1 = KO1 + K0
                K2 = K1  + 1
                NXR1  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K1,N) &
                        - NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K2,N)
                NXR2  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K2,N) &
                        + NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K1,N)
                PXR1  =   PCON(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N) &
                        - PCON(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)
                PXR2  =   PCON(K1,N,Q) *   SOLB_XNEG(I,O1,K2,N) &
                        + PCON(K2,N,Q) *   SOLB_XNEG(I,O1,K1,N)
                HOM1CR =   NXR1  * T_UTDN_EIGEN(K1,UT) &
                         - NXR2  * T_UTDN_EIGEN(K2,UT)
                HOM2CR =   PXR1  * T_UTUP_EIGEN(K1,UT) &
                         - PXR2  * T_UTUP_EIGEN(K2,UT)
                SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
              ENDDO

!  real part

              SHOM = SHOM_R + SHOM_CR
              QATMOSWF_F(Q,UTA,I,O1) = FMULT * SHOM

!  Finish Q, I and O1 loops

            ENDDO
          ENDDO
        ENDDO

!  end variability clause

      ENDIF

!  Add the linearized thermal solution  (if flagged)
!    ---Only present if N = NV, or NV = 0
!   THIS IS THE SOLUTION with scattering

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        O1 = 1
        IF ( N.EQ.NV .OR. NV.EQ.0 ) THEN
          DO I = 1, NSTREAMS
            DO Q = 1, NV_PARAMETERS
              SPAR = L_UT_T_PARTIC(I,UT,Q)
              QATMOSWF_F(Q,UTA,I,O1) = &
                 QATMOSWF_F(Q,UTA,I,O1) + FMULT * SPAR
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  Finished if no solar terms

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Add Linearized Classical beam solution
!  Green function solution is still a Placeholder....!
!    solutions exist only for N = NV or N > NV or NV = 0

! @@@@@@@@@@ Rob fix, 2/9/11, Replacement of commented section
!      IF ( DO_CLASSICAL_SOLUTION ) THEN
!        IF ( N.GE.NV .OR. NV.EQ.0 ) THEN
!          DO I = 1, NSTREAMS
!            DO O1 = 1, NSTOKES
!              DO Q = 1, NV_PARAMETERS
!                SPAR = L_WUPPER(I,O1,N,Q) *   T_UTDN_MUBAR(UT,IB) + &
!                         WUPPER(I,O1,N)   * L_T_UTDN_MUBAR(UT,NV,IB,Q)
!                QATMOSWF_F(Q,UTA,I,O1) = &
!                  QATMOSWF_F(Q,UTA,I,O1) + FMULT * SPAR
!              ENDDO
!            ENDDO
!          ENDDO
!        ENDIF
!      ELSE
!                P L A C E H O L D E R
!      ENDIF

      IF ( DO_CLASSICAL_SOLUTION ) THEN
        IF ( DO_INCLUDE_THERMEMISS ) THEN
          VAREX = ( N.EQ.NV .OR. NV.EQ.0 )
          IF ( N.GT.NV .or. VAREX ) THEN
            DO I = 1, NSTREAMS
              DO O1 = 1, NSTOKES
                DO Q = 1, NV_PARAMETERS
                  SUNP   = WUPPER(I,O1,N)
                  L_SUNP = L_WUPPER(I,O1,N,Q)
                  IF ( O1 .eq. 1 ) THEN
                    SUNP   = SUNP - T_WUPPER(I,N)
                    IF (VAREX) L_SUNP = L_SUNP - L_T_WUPPER(I,N,Q)
                  ENDIF
                  SPAR = L_SUNP *   T_UTDN_MUBAR(UT,IB) + &
                           SUNP * L_T_UTDN_MUBAR(UT,NV,IB,Q)
                  QATMOSWF_F(Q,UTA,I,O1) = &
                  QATMOSWF_F(Q,UTA,I,O1) + FMULT * SPAR
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ELSE
          IF ( N.GE.NV .OR. NV.EQ.0 ) THEN
            DO I = 1, NSTREAMS
              DO O1 = 1, NSTOKES
                DO Q = 1, NV_PARAMETERS
                  SPAR = L_WUPPER(I,O1,N,Q) *   T_UTDN_MUBAR(UT,IB) + &
                           WUPPER(I,O1,N)   * L_T_UTDN_MUBAR(UT,NV,IB,Q)
                  QATMOSWF_F(Q,UTA,I,O1) = &
                    QATMOSWF_F(Q,UTA,I,O1) + FMULT * SPAR
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDIF
      ENDIF
! @@@@@@@@@@ Rob fix, 2/9/11, End Replacement section

!  Finish

      RETURN
      END SUBROUTINE QUADATMOSWF_OFFGRID_DN


      END MODULE vlidort_l_wfatmos

