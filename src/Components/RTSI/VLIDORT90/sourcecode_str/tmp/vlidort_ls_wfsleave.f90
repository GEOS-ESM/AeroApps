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

! ##########################################################
! #                                                        #
! # Subroutines in this Module                             #
! #                                                        #
! #     Top level PUBLIC routines--------------            #
! #            VLIDORT_LSSL_DBCORRECTION (direct)          #
! #            VLIDORT_LSSL_DBSETUPS (diffuse)             #
! #            VLIDORT_LSSL_WFS   (diffuse, master)        #
! #                                                        #
! #     High level Jacobian routines  ---------            #
! #            LSSL_UPUSER_SURFACEWF                       #
! #            LSSL_DNUSER_SURFACEWF                       #
! #                                                        #
! #     Post-processing at user angles --------            #
! #            LSSL_WHOLELAYER_STERM_UP                    #
! #            LSSL_WHOLELAYER_STERM_DN                    #
! #            LSSL_PARTLAYER_STERM_UP                     #
! #            LSSL_PARTLAYER_STERM_DN                     #
! #                                                        #
! ##########################################################
! #                                                        #
! #     High level Jacobian routine   ---------            #
! #            LSSL_INTEGRATED_OUTPUT (master)             #
! #                                                        #
! #     Post-processing at quad angles --------            #
! #            LSSL_QUADWF_LEVEL_UP                        #
! #            LSSL_QUADWF_LEVEL_DN                        #
! #            LSSL_QUADWF_OFFGRID_UP                      #
! #            LSSL_QUADWF_OFFGRID_DN                      #
! #                                                        #
! ##########################################################

      MODULE vlidort_ls_wfsleave

!  Only 3 routines available publicly to the rest of VLIDORT...

      PRIVATE
      PUBLIC :: VLIDORT_LSSL_DBCORRECTION, &
                VLIDORT_LSSL_DBSETUPS, &
                VLIDORT_LSSL_WFS

      CONTAINS

!

      SUBROUTINE VLIDORT_LSSL_DBCORRECTION ( &
        DO_SSCORR_OUTGOING, DO_UPWELLING, DO_SL_ISOTROPIC,       &
        DO_OBSERVATION_GEOMETRY,                                 &
        NSTOKES, NLAYERS, NBEAMS, N_SLEAVE_WFS, N_REFLEC_WFS,    &
        N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS,           &
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                 &
        PARTLAYERS_LAYERIDX, UTAU_LEVEL_MASK_UP,                 &
        N_GEOMETRIES, VZA_OFFSETS, DO_REFLECTED_DIRECTBEAM,      &
        FLUXMULT, LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_USERANGLES, &
        T_DELT_USERM, T_UTUP_USERM, UP_LOSTRANS, UP_LOSTRANS_UT, &
        SURFACEWF_DB )

!  PREPARES LINEARIZATION OF EXACT DIRECT BEAM REFLECTION
!   LINEARIZATION WITH RESPECT TO ALL SURFACE PARAMETERS

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_SL_ISOTROPIC
      LOGICAL, INTENT (IN) ::          DO_SSCORR_OUTGOING
      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NBEAMS
      INTEGER, INTENT (IN) ::          N_SLEAVE_WFS
      INTEGER, INTENT (IN) ::          N_REFLEC_WFS

      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          N_USER_RELAZMS
      INTEGER, INTENT (IN) ::          N_USER_LEVELS

      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

      INTEGER, INTENT (IN) ::          N_GEOMETRIES
      INTEGER, INTENT (IN) ::          VZA_OFFSETS ( MAX_SZANGLES, MAX_USER_VZANGLES )
      LOGICAL, INTENT (IN) ::          DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )

      DOUBLE PRECISION, INTENT (IN) :: FLUXMULT
      DOUBLE PRECISION, INTENT (IN) ::   LSSL_SLTERM_ISOTROPIC &
          ( MAX_SLEAVEWFS, MAXSTOKES, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::   LSSL_SLTERM_USERANGLES &
          ( MAX_SLEAVEWFS, MAXSTOKES, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS  )

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: UP_LOSTRANS    ( MAXLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (IN) :: UP_LOSTRANS_UT ( MAX_PARTLAYERS, MAX_GEOMETRIES )

!  Output
!  ------

      DOUBLE PRECISION, INTENT (INOUT) :: SURFACEWF_DB &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  LOCAL VARIABLES
!  ---------------

!  Local lsource arrays (not saved this time)

      DOUBLE PRECISION  :: LSSL_EXACTDB_SOURCE ( MAX_SLEAVEWFS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION  :: LSSL_DB_CUMSOURCE   ( MAX_GEOMETRIES, MAXSTOKES )

!  Help

      INTEGER ::          N, NUT, NSTART, NUT_PREV, NLEVEL, NELEMENTS
      INTEGER ::          UT, UTA, UM, UA, NC, IB, V, O1, Q, Q1, LUM, LUA
      DOUBLE PRECISION :: TR

!  Local user indices

      LUM = 1
      LUA = 1

!  FIRST STAGE
!  -----------

!  RETURN IF NO UPWELLING

      IF ( .NOT.DO_UPWELLING ) RETURN

!  LAMBERTIAN: ONLY 1 STOKES COMPONENT, ONLY ONE SURFACE WEIGHTING FUNCTION

      IF ( DO_SL_ISOTROPIC ) THEN
        NELEMENTS = 1
      ELSE
        NELEMENTS = NSTOKES
      ENDIF

!  INITIALISE OUTPUT

      DO V = 1, N_GEOMETRIES
        DO UTA = 1, N_USER_LEVELS
          DO Q = 1, N_SLEAVE_WFS
            Q1 = Q + N_REFLEC_WFS
            DO O1 = 1, NELEMENTS
              SURFACEWF_DB(Q1,UTA,V,O1)  = ZERO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!  GET THE SLEAVING TERM LINEARIZATION

      IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
        DO UM = 1, N_USER_STREAMS
          DO IB = 1, NBEAMS
            DO UA = 1, N_USER_RELAZMS
              V = VZA_OFFSETS(IB,UM) + UA
              IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
                DO Q = 1, N_SLEAVE_WFS
                  DO O1 = 1, NELEMENTS
                    IF ( DO_SL_ISOTROPIC ) THEN
                      LSSL_EXACTDB_SOURCE(Q,V,O1) = LSSL_SLTERM_ISOTROPIC(Q,O1,IB)*PI4
                    ELSE
                      LSSL_EXACTDB_SOURCE(Q,V,O1) = LSSL_SLTERM_USERANGLES(Q,O1,UM,UA,IB)*PI4
                    ENDIF
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO V = 1, N_GEOMETRIES
          IF ( DO_REFLECTED_DIRECTBEAM(V) ) THEN
            DO Q = 1, N_SLEAVE_WFS
              DO O1 = 1, NELEMENTS
                IF ( DO_SL_ISOTROPIC ) THEN
                  LSSL_EXACTDB_SOURCE(Q,V,O1) = LSSL_SLTERM_ISOTROPIC(Q,O1,V)*PI4
                ELSE
                  LSSL_EXACTDB_SOURCE(Q,V,O1) = LSSL_SLTERM_USERANGLES(Q,O1,LUM,LUA,V)*PI4
                ENDIF
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDIF

!  ====================================================
!  2. TRANSMITTANCE OF SLEAVE TERM: UPWELLING RECURSION
!  ====================================================

      DO Q = 1, N_SLEAVE_WFS

!  Offset

        Q1 = Q + N_REFLEC_WFS

!  INITIALIZE CUMULATIVE SOURCE TERM

        NC =  0
        DO V = 1, N_GEOMETRIES
         DO O1 = 1, NELEMENTS
          LSSL_DB_CUMSOURCE(V,O1) = LSSL_EXACTDB_SOURCE(Q,V,O1)*FLUXMULT
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
           IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
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
                  LSSL_DB_CUMSOURCE(V,O1) = TR * LSSL_DB_CUMSOURCE(V,O1)
                ENDDO
               ENDDO
              ENDDO
             ENDDO
           ELSE
            DO V = 1, N_GEOMETRIES
              IF ( DO_SSCORR_OUTGOING ) THEN
                TR  = UP_LOSTRANS(N,V)
              ELSE
                TR  = T_DELT_USERM(N,V)
              ENDIF
              DO O1 = 1, NELEMENTS
                LSSL_DB_CUMSOURCE(V,O1) = TR * LSSL_DB_CUMSOURCE(V,O1)
              ENDDO
            ENDDO
           ENDIF
          ENDDO

!  OFFGRID OUTPUT
!  --------------

!  REQUIRE PARTIAL LAYER TRANSMITTANCE

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
                ELSE
                  TR  = T_UTUP_USERM(UT,UM)
                ENDIF
                DO O1 = 1, NELEMENTS
                 SURFACEWF_DB(Q1,UTA,V,O1) = TR * LSSL_DB_CUMSOURCE(V,O1)
                ENDDO
               ENDDO
              ENDDO
             ENDDO
           ELSE
             DO V = 1, N_GEOMETRIES
              IF ( DO_SSCORR_OUTGOING ) THEN
                TR  = UP_LOSTRANS_UT(UT,V)
              ELSE
                TR  = T_UTUP_USERM(UT,V)
              ENDIF
              DO O1 = 1, NELEMENTS
               SURFACEWF_DB(Q1,UTA,V,O1) = TR * LSSL_DB_CUMSOURCE(V,O1)
              ENDDO
             ENDDO
           ENDIF


!  ONGRID OUTPUT : SET FINAL CUMULATIVE SOURCE DIRECTLY

          ELSE
           IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
             DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
               DO UA = 1, N_USER_RELAZMS
                V = VZA_OFFSETS(IB,UM) + UA
                DO O1 = 1, NELEMENTS
                 SURFACEWF_DB(Q1,UTA,V,O1) = LSSL_DB_CUMSOURCE(V,O1)
                ENDDO
               ENDDO
              ENDDO
             ENDDO
           ELSE
             DO V = 1, N_GEOMETRIES
              DO O1 = 1, NELEMENTS
               SURFACEWF_DB(Q1,UTA,V,O1) = LSSL_DB_CUMSOURCE(V,O1)
              ENDDO
             ENDDO
           ENDIF
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
      END SUBROUTINE VLIDORT_LSSL_DBCORRECTION

!

      SUBROUTINE VLIDORT_LSSL_DBSETUPS ( &
        DO_OBSERVATION_GEOMETRY, DO_SL_ISOTROPIC, DO_USER_STREAMS,    &
        DO_REFLECTED_DIRECTBEAM, FOURIER_COMPONENT,                   &
        NSTOKES, NSTREAMS, NBEAMS, N_USER_STREAMS, N_SLEAVE_WFS,      &
        FLUX_FACTOR, DELTA_FACTOR,                                    &
        LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0, LSSL_USER_SLTERM_F_0, &
        LSSL_DIRECT_BEAM, LSSL_USER_DIRECT_BEAM )

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::           DO_OBSERVATION_GEOMETRY
      LOGICAL, INTENT (IN) ::           DO_SL_ISOTROPIC
      LOGICAL, INTENT (IN) ::           DO_USER_STREAMS
      LOGICAL, INTENT (IN) ::           DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )

      INTEGER, INTENT (IN) ::           FOURIER_COMPONENT
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NBEAMS
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      INTEGER, INTENT (IN) ::           N_SLEAVE_WFS

      DOUBLE PRECISION, INTENT (IN) ::  FLUX_FACTOR, DELTA_FACTOR
      DOUBLE PRECISION, INTENT (IN) ::   LSSL_SLTERM_ISOTROPIC &
          ( MAX_SLEAVEWFS, MAXSTOKES, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::   LSSL_SLTERM_F_0 &
          ( MAX_SLEAVEWFS, 0:MAXMOMENTS, MAXSTOKES, MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::   LSSL_USER_SLTERM_F_0 &
          ( MAX_SLEAVEWFS, 0:MAXMOMENTS, MAXSTOKES, MAX_USER_STREAMS, MAXBEAMS )

!  Outputs

      DOUBLE PRECISION, INTENT (OUT) :: LSSL_DIRECT_BEAM &
          ( MAX_SLEAVEWFS, MAXSTREAMS, MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: LSSL_USER_DIRECT_BEAM &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )

!  Local variables
!  ---------------

      DOUBLE PRECISION :: SL, HELP
      INTEGER          :: I, UI, O1, IB, M, Q, LUM

!  Initialize
!  ----------

!  Local user index

      LUM = 1

!  Safety first!

      DO IB = 1, NBEAMS
        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES
           LSSL_DIRECT_BEAM(1:N_SLEAVE_WFS,I,IB,O1) = ZERO
          ENDDO
        ENDDO
        IF ( DO_USER_STREAMS ) THEN
          IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
            DO UI = 1, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                LSSL_USER_DIRECT_BEAM(1:N_SLEAVE_WFS,UI,IB,O1) = ZERO
              ENDDO
            ENDDO
          ELSE
            DO O1 = 1, NSTOKES
              LSSL_USER_DIRECT_BEAM(1:N_SLEAVE_WFS,LUM,IB,O1) = ZERO
            ENDDO
          ENDIF
        ENDIF
      ENDDO

!  Fourier component

      M = FOURIER_COMPONENT

!  Attenuation of solar beam
!  -------------------------

!  New code to deal with refractive geometry case
!    R. Spurr, 7 May 2005. RT Solutions Inc.

      DO IB = 1, NBEAMS
        IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN

!  Corrected implementation, 30 July 2012
!    Normalized to Flux-factor / DELTA_Factor
!    Delta_Factor = 1.0 for the Isotropic or non-iso Fourier = 0 cases

          HELP = FLUX_FACTOR / DELTA_FACTOR
          IF ( DO_SL_ISOTROPIC .and. M.EQ.0 ) THEN
            DO Q = 1, N_SLEAVE_WFS
              DO O1 = 1, NSTOKES
                SL = LSSL_SLTERM_ISOTROPIC(Q,O1,IB) * HELP
                LSSL_DIRECT_BEAM(Q,1:NSTREAMS,IB,O1) = SL
                IF ( DO_USER_STREAMS ) THEN
                  LSSL_USER_DIRECT_BEAM(Q,1:N_USER_STREAMS,IB,O1) = SL
                ELSE
                  LSSL_USER_DIRECT_BEAM(Q,LUM,IB,O1) = SL
                ENDIF
              ENDDO
            ENDDO
          ELSE
            DO Q = 1, N_SLEAVE_WFS
              DO O1 = 1, NSTOKES
                DO I = 1, NSTREAMS
                  SL = LSSL_SLTERM_F_0(Q,M,O1,I,IB) * HELP
                  LSSL_DIRECT_BEAM(Q,I,IB,O1) = SL
                ENDDO
                IF ( DO_USER_STREAMS ) THEN
                  IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
                    DO UI = 1, N_USER_STREAMS
                      SL = LSSL_USER_SLTERM_F_0(Q,M,O1,UI,IB) * HELP
                      LSSL_USER_DIRECT_BEAM(Q,UI,IB,O1) = SL
                    ENDDO
                  ELSE
                    SL = LSSL_USER_SLTERM_F_0(Q,M,O1,LUM,IB) * HELP
                    LSSL_USER_DIRECT_BEAM(Q,LUM,IB,O1) = SL
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDIF

!  end direct beam calculation

       ENDIF
      ENDDO

!  finish

      RETURN
      END SUBROUTINE VLIDORT_LSSL_DBSETUPS

!

      SUBROUTINE VLIDORT_LSSL_WFS ( &
        DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_MVOUTPUT,                  &
        DO_QUAD_OUTPUT, DO_USER_STREAMS, DO_UPWELLING, DO_DNWELLING, &
        DO_SL_ISOTROPIC, DO_OBSERVATION_GEOMETRY,                    &
        N_SLEAVE_WFS, N_REFLEC_WFS, NSTOKES, NSTREAMS, NLAYERS,      &
        N_USER_STREAMS, LOCAL_UM_START, FOURIER_COMPONENT, IBEAM,    &
        N_USER_LEVELS, N_DIRECTIONS, WHICH_DIRECTIONS,               &
        NTOTAL, N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2,  &
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                     &
        UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX, &
        DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, USER_BRDF_F,       & ! @@@ New Line
        SURFACE_FACTOR, FLUX_MULTIPLIER, QUAD_WEIGHTS, QUAD_STRMWTS, &
        LSSL_DIRECT_BEAM, LSSL_USER_DIRECT_BEAM,            &
        K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,            &
        T_DELT_EIGEN, BANDMAT2, IPIVOT, SMAT2, SIPIVOT,     &
        T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,           &
        T_UTUP_EIGEN, T_UTDN_EIGEN,  HMULT_1, HMULT_2,      &
        UHOM_UPDN, UHOM_UPUP, UHOM_DNDN, UHOM_DNUP,         &
        UT_HMULT_UU, UT_HMULT_UD, UT_HMULT_DU, UT_HMULT_DD, &
        SURFACEWF_F, MINT_SURFACEWF, FLUX_SURFACEWF,        &
        STATUS, MESSAGE, TRACE )

      USE VLIDORT_PARS
      USE LAPACK_TOOLS

      IMPLICIT NONE

!  Control

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_MVOUTPUT
      LOGICAL, INTENT (IN) ::          DO_QUAD_OUTPUT

      LOGICAL, INTENT (IN) ::          DO_USER_STREAMS
      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_DNWELLING
      LOGICAL, INTENT (IN) ::          DO_SL_ISOTROPIC
      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY

!      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
!      DOUBLE PRECISION, INTENT (IN) :: FLUXVEC ( MAXSTOKES )

      INTEGER, INTENT (IN) ::          N_SLEAVE_WFS, N_REFLEC_WFS
      INTEGER, INTENT (IN) ::          FOURIER_COMPONENT, IBEAM
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS

      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START

      INTEGER, INTENT (IN) ::          N_DIRECTIONS
      INTEGER, INTENT (IN) ::          WHICH_DIRECTIONS ( MAX_DIRECTIONS )

      INTEGER, INTENT (IN) ::          NTOTAL
      INTEGER, INTENT (IN) ::          N_SUBDIAG
      INTEGER, INTENT (IN) ::          N_SUPDIAG
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS_2

      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

!  @@@@@@@@@@@@@@@@ Rob Fix @@@ 11 Sep 12 ADDITIONAL CODE  @@@@@@@@@@@
!  @@@@@@@@@@@@@@@@ START OF BLOCK  @@@@@@@@@@@
!    Need extra line of input variables
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      DOUBLE PRECISION, INTENT (IN) :: LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: USER_BRDF_F &
          ( 0:MAXMOMENTS,MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )
!  @@@@@@@@@@@@@@@@ END   OF BLOCK  @@@@@@@@@@@
!  @@@@@@@@@@@@@@@@ Rob Fix @@@ 11 Sep 12 ADDITIONAL CODE  @@@@@@@@@@@

!  Solutions

      DOUBLE PRECISION, INTENT (IN) :: SURFACE_FACTOR
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER

      DOUBLE PRECISION, INTENT (IN) :: QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )

      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER, INTENT (IN) ::          IPIVOT ( MAXTOTAL )
      DOUBLE PRECISION, INTENT (IN) :: SMAT2 ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )
      INTEGER, INTENT (IN) ::          SIPIVOT ( MAXSTRMSTKS_2 )

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN &
          ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN &
          ( MAXEVALUES, MAX_PARTLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!mick fix 9/7/2012 - moved MAX_SLEAVEWFS to 1st dimension
      DOUBLE PRECISION, INTENT (IN) :: LSSL_DIRECT_BEAM &
          ( MAX_SLEAVEWFS, MAXSTREAMS, MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: LSSL_USER_DIRECT_BEAM &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )

      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPDN &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPUP &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNDN &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNUP &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: HMULT_1 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_2 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_UU &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_UD &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_DU &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_DD &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )


!  Linearized surface output

      DOUBLE PRECISION, INTENT (INOUT) ::  SURFACEWF_F &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_VZANGLES, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) ::  MINT_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT)  :: FLUX_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  Exception handling

      INTEGER, INTENT (OUT) ::             STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  Local variables
!  ---------------

!  Linearized BOA terms

      DOUBLE PRECISION :: LSSL_BOA_SOURCE &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAXSTOKES )

!  Linearized BVP solution

      DOUBLE PRECISION :: COL2_WFSLEAVE  ( MAXTOTAL, MAX_SLEAVEWFS )
      DOUBLE PRECISION :: SCOL2_WFSLEAVE ( MAXSTRMSTKS_2, MAX_SLEAVEWFS )

      DOUBLE PRECISION :: NCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION :: PCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )

!  Other local variables

      INTEGER ::           INFO, NSTOKES_ACTUAL, IB, M, N, Q, O1
      INTEGER ::           K, K0, K1, K2, KO1, C0, CM, I, UM, LUM
      INTEGER ::           IR, IROW, IROW1, IROW_S, IROW1_S
      DOUBLE PRECISION  :: KS
      CHARACTER (LEN=3) :: CI

!  @@@@@@@@@@@@@@@@ Rob Fix @@@ 11 Sep 12 ADDITIONAL CODE  @@@@@@@@@@@
!  @@@@@@@@@@@@@@@@ START OF BLOCK  @@@@@@@@@@@
     INTEGER          :: O2, J, OM
      DOUBLE PRECISION :: INTEGRAND ( MAX_SLEAVEWFS, MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION :: SUM_R, SUM_CR, REFLEC, S_REFLEC, REFL_ATTN
      DOUBLE PRECISION :: H1, H2, NXR, PXR, NXR1, NXR2, PXR1
!  @@@@@@@@@@@@@@@@ END   OF BLOCK  @@@@@@@@@@@
!  @@@@@@@@@@@@@@@@ Rob Fix @@@ 11 Sep 12 ADDITIONAL CODE  @@@@@@@@@@@

!  Initialise status

      STATUS  = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Short-hand

      IB  = IBEAM
      M   = FOURIER_COMPONENT
      LUM = 1

!  Actual number of Stokes calculations required

      NSTOKES_ACTUAL = NSTOKES
      IF ( DO_SL_ISOTROPIC ) NSTOKES_ACTUAL = 1

!  Nothing to do if M > 0 and Isotropic

!mick fix - turned off
      !IF ( .not. DO_INCLUDE_DIRECTBEAM  ) RETURN
      IF ( M.gt.0 .and. DO_SL_ISOTROPIC ) RETURN

!  BVP solution for perturbed integration constants
!  ------------------------------------------------

!  Compute the main column B' where AX = B'
!     Regular BVP Solution --->  NO TELESCOPING HERE

!  initialise. Vitally necessary

      DO Q = 1, N_SLEAVE_WFS
        DO I = 1, NTOTAL
          COL2_WFSLEAVE(I,Q) = ZERO
        ENDDO
      ENDDO

!  Last layer : Add direct beam variation due to surface leaving 

      N  = NLAYERS
      C0 = N*NSTKS_NSTRMS_2 - NSTKS_NSTRMS
      DO Q = 1, N_SLEAVE_WFS
        DO O1 = 1, NSTOKES_ACTUAL
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            CM   = C0 + IR + O1
!mick fix 9/7/2012 - moved Q to 1st dimension
! Rob Fix @@@ 11 Sep 12, remove SURFACE_FACTOR
            !COL2_WFSLEAVE(CM,Q) = SURFACE_FACTOR * LSSL_DIRECT_BEAM(I,IB,O1,Q)
            !COL2_WFSLEAVE(CM,Q) = SURFACE_FACTOR * LSSL_DIRECT_BEAM(Q,I,IB,O1)
            COL2_WFSLEAVE(CM,Q) =  LSSL_DIRECT_BEAM(Q,I,IB,O1)
          ENDDO
        ENDDO
      ENDDO

!  Copy for the single layer case

      IF ( NLAYERS .EQ. 1 ) THEN
!mick fix 9/17/2012 - replaced head of do loop
        !DO O1 = 1, NSTOKES_ACTUAL
        DO Q = 1, N_SLEAVE_WFS
          DO N = 1, NTOTAL
            SCOL2_WFSLEAVE(N,Q) = COL2_WFSLEAVE(N,Q)
          ENDDO
        ENDDO
      ENDIF

!  BVP back-substitution: With compression (multilayers)
!  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  LAPACK substitution (DGBTRS) using RHS column vector COL2_WF
!  BV solution for perturbed integration constants
!    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS &
           ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, N_SLEAVE_WFS, &
              BANDMAT2, MAXBANDTOTAL, IPIVOT, &
              COL2_WFSLEAVE, MAXTOTAL, INFO )

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRS call (multilayer) in VLIDORT_SLEAVE WFS'
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  Set Linearized integration constants NCON_SLEAVE and PCON_SLEAVE, all layer

        DO Q = 1, N_SLEAVE_WFS
          DO N = 1, NLAYERS
            C0 = (N-1)*NSTKS_NSTRMS_2
            DO K = 1, K_REAL(N)
              IROW = K
              IROW1 = IROW + NSTKS_NSTRMS
              NCON_SLEAVE(Q,K,N) = COL2_WFSLEAVE(C0+IROW,Q)
              PCON_SLEAVE(Q,K,N) = COL2_WFSLEAVE(C0+IROW1,Q)
            ENDDO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2 * K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              IROW    = K + K_REAL(N)
              IROW1   = IROW + NSTKS_NSTRMS
              IROW_S  = IROW + K_COMPLEX(N)
              IROW1_S = IROW_S + NSTKS_NSTRMS
              NCON_SLEAVE(Q,K1,N) = COL2_WFSLEAVE(C0+IROW,   Q)
              NCON_SLEAVE(Q,K2,N) = COL2_WFSLEAVE(C0+IROW_S, Q)
              PCON_SLEAVE(Q,K1,N) = COL2_WFSLEAVE(C0+IROW1,  Q)
              PCON_SLEAVE(Q,K2,N) = COL2_WFSLEAVE(C0+IROW1_S,Q)
            ENDDO
          ENDDO
        ENDDO

!  Solve the boundary problem: No compression, Single Layer only
!  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2_WFSLEAVE

        CALL DGETRS &
           ( 'N', NTOTAL, N_SLEAVE_WFS, SMAT2, MAXSTRMSTKS_2, &
              SIPIVOT, SCOL2_WFSLEAVE, MAXSTRMSTKS_2, INFO )

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRS call (Reg. 1 layer) in VLIDORT_SLEAVE_WFS'
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  Set Linearized integration constants NCON_SLEAVE and PCON_SLEAVE, 1 layer

        DO Q = 1, N_SLEAVE_WFS
          N = 1
          DO K = 1, K_REAL(N)
            IROW = K
            IROW1 = IROW + NSTKS_NSTRMS
            NCON_SLEAVE(Q,K,N) = SCOL2_WFSLEAVE(IROW,Q)
            PCON_SLEAVE(Q,K,N) = SCOL2_WFSLEAVE(IROW1,Q)
          ENDDO
          KO1 = K_REAL(N) + 1
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            IROW    = K + K_REAL(N)
            IROW1   = IROW + NSTKS_NSTRMS
            IROW_S  = IROW + K_COMPLEX(N)
            IROW1_S = IROW_S + NSTKS_NSTRMS
            NCON_SLEAVE(Q,K1,N) = SCOL2_WFSLEAVE(IROW,   Q)
            NCON_SLEAVE(Q,K2,N) = SCOL2_WFSLEAVE(IROW_S, Q)
            PCON_SLEAVE(Q,K1,N) = SCOL2_WFSLEAVE(IROW1,  Q)
            PCON_SLEAVE(Q,K2,N) = SCOL2_WFSLEAVE(IROW1_S,Q)
          ENDDO
        ENDDO

!  end clause

      ENDIF

!  debug------------------------------------------
!        if ( do_debug_write.and.fourier_component.eq.0 ) then
!         DO N = 1, NLAYERS
!          DO K = 1, K_REAL(N)
!           write(86,'(3i2,1p6e13.5)')FOURIER_COMPONENT,N,K,
!     &                LCON(K,N), MCON(K,N),
!     &                NCON_SLEAVE(1,K,N),PCON_SLEAVE(1,K,N)
!          ENDDO
!         ENDDO
!        ENDIF

!  Get the Post-processed weighting functions
!  ==========================================

!  Upwelling weighting functions
!  -----------------------------

      IF ( DO_UPWELLING .and. DO_USER_STREAMS ) THEN

!  Derivative of BOA source function

!        KS = SURFACE_FACTOR. Bug rob fix 9/9/14. Should be 1.0 always
        KS = one

!mick fix 9/11/2012 - modified IF condition and added ELSE
!        IF ( DO_INCLUDE_DIRECTBEAM ) THEN
!          IF ( (DO_SL_ISOTROPIC .AND. M.EQ.0) .OR. .NOT.DO_SL_ISOTROPIC ) THEN
        IF ( DO_INCLUDE_DIRECTBEAM .AND. &
             ((DO_SL_ISOTROPIC .AND. M.EQ.0) .OR. .NOT.DO_SL_ISOTROPIC) ) THEN

          DO O1 = 1, NSTOKES_ACTUAL
            IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
              DO Q = 1, N_SLEAVE_WFS
! Rob Fix @@@ 11 Sep 12, Wrong range for UM do loop
!               DO UM = 1, LOCAL_UM_START
                DO UM = LOCAL_UM_START, N_USER_STREAMS
!mick fix 9/7/2012 - moved Q to 1st dimension
                  !LSSL_BOA_SOURCE(Q,UM,O1) = KS * LSSL_USER_DIRECT_BEAM(UM,IB,O1,Q)
                  LSSL_BOA_SOURCE(Q,UM,O1) = KS * LSSL_USER_DIRECT_BEAM(Q,UM,IB,O1)
                ENDDO
              ENDDO
            ELSE
              DO Q = 1, N_SLEAVE_WFS
                LSSL_BOA_SOURCE(Q,IB,O1) = KS * LSSL_USER_DIRECT_BEAM(Q,LUM,IB,O1)
              ENDDO
            ENDIF
          ENDDO

        ELSE
          LSSL_BOA_SOURCE = ZERO
        ENDIF

!  @@@@@@@@@@@@@@@@ Rob Fix @@@ 11 Sep 12 ADDITIONAL CODE  @@@@@@@@@@@
!  @@@@@@@@@@@@@@@@ START OF BLOCK  @@@@@@@@@@@

!  Diffuse Term: Contribution due to derivatives of BVP constants
!  -------------------------------------------------------------

!  First compute derivative of downward intensity Integrand at stream an
!        .. reflectance integrand  = a(j).x(j).dI_DOWN(-j)/dS

!  start loops

        N = NLAYERS
        DO Q = 1, N_SLEAVE_WFS
         DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES_ACTUAL

!  Real homogeneous solutions

           SUM_R = ZERO
           DO K = 1, K_REAL(N)
            NXR = NCON_SLEAVE(Q,K,N) * SOLA_XPOS(I,O1,K,N)
            PXR = PCON_SLEAVE(Q,K,N) * SOLB_XNEG(I,O1,K,N)
            SUM_R = SUM_R + NXR*T_DELT_EIGEN(K,N) + PXR
           ENDDO

!  Complex solutions

           SUM_CR = ZERO
           DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            NXR1 =   NCON_SLEAVE(Q,K1,N) * SOLA_XPOS(I,O1,K1,N) &
                   - NCON_SLEAVE(Q,K2,N) * SOLA_XPOS(I,O1,K2,N)
            NXR2 =   NCON_SLEAVE(Q,K1,N) * SOLA_XPOS(I,O1,K2,N) &
                   + NCON_SLEAVE(Q,K2,N) * SOLA_XPOS(I,O1,K1,N)
            PXR1 =   PCON_SLEAVE(Q,K1,N) * SOLB_XNEG(I,O1,K1,N) &
                   - PCON_SLEAVE(Q,K2,N) * SOLB_XNEG(I,O1,K2,N)
            H1 =  NXR1 * T_DELT_EIGEN(K1,N) &
                 -NXR2 * T_DELT_EIGEN(K2,N)
            H2 =  PXR1
            SUM_CR = SUM_CR + H1 + H2
           ENDDO

!  Final result

           INTEGRAND(Q,I,O1) = QUAD_STRMWTS(I) * ( SUM_R + SUM_CR )

!  end loops

          ENDDO
         ENDDO
        ENDDO

!  integrated reflectance term
!  ---------------------------

!  Lambertian case, same for all user-streams

        IF ( DO_LAMBERTIAN_SURFACE ) THEN
         IF ( FOURIER_COMPONENT.EQ.0 ) THEN
          O1 = 1
          DO Q = 1, N_SLEAVE_WFS
           REFLEC = ZERO
           DO J = 1, NSTREAMS
            REFLEC = REFLEC + INTEGRAND(Q,J,O1)
           ENDDO
           REFLEC = SURFACE_FACTOR * REFLEC * LAMBERTIAN_ALBEDO
           IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
             DO UM = LOCAL_UM_START, N_USER_STREAMS
              LSSL_BOA_SOURCE(Q,UM,O1) = LSSL_BOA_SOURCE(Q,UM,O1) + REFLEC
             ENDDO
           ELSE
             LSSL_BOA_SOURCE(Q,IB,O1) = LSSL_BOA_SOURCE(Q,IB,O1) + REFLEC
           ENDIF
          ENDDO
         ENDIF
        ENDIF

!  BRDF case

        IF ( .not. DO_LAMBERTIAN_SURFACE ) THEN
         DO Q = 1, N_SLEAVE_WFS
          IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
             DO O1 = 1, NSTOKES_ACTUAL
              REFLEC = ZERO
              DO J = 1, NSTREAMS
               S_REFLEC = ZERO
               DO O2 = 1, NSTOKES
                OM = MAXSTOKES*(O1-1) + O2
                S_REFLEC = S_REFLEC + INTEGRAND(Q,J,O2) * &
                                      USER_BRDF_F(M,OM,UM,J)
               ENDDO
               REFLEC = REFLEC + S_REFLEC
              ENDDO
              LSSL_BOA_SOURCE(Q,UM,O1) = LSSL_BOA_SOURCE(Q,UM,O1) + &
                                         REFLEC * SURFACE_FACTOR
             ENDDO
            ENDDO
          ELSE
           DO O1 = 1, NSTOKES_ACTUAL
            REFLEC = ZERO
            DO J = 1, NSTREAMS
             S_REFLEC = ZERO
             DO O2 = 1, NSTOKES
              OM = MAXSTOKES*(O1-1) + O2
              S_REFLEC = S_REFLEC + INTEGRAND(Q,J,O2) * &
                                    USER_BRDF_F(M,OM,IB,J)
             ENDDO
             REFLEC = REFLEC + S_REFLEC
            ENDDO
            LSSL_BOA_SOURCE(Q,IB,O1) = LSSL_BOA_SOURCE(Q,IB,O1) + &
                                          REFLEC * SURFACE_FACTOR
           ENDDO
          ENDIF
         ENDDO
        ENDIF

!  @@@@@@@@@@@@@@@@ END OF BLOCK    @@@@@@@@@@@
!  @@@@@@@@@@@@@@@@ Rob Fix @@@ 11 Sep 12 ADDITIONAL CODE  @@@@@@@@@@@


!  Upwelling Surface WF contribution for SLEAVE

        CALL LSSL_UPUSER_SURFACEWF ( &
          DO_OBSERVATION_GEOMETRY,                       &
          N_SLEAVE_WFS, N_REFLEC_WFS, NSTOKES, NLAYERS,  &
          N_USER_LEVELS, N_USER_STREAMS, LOCAL_UM_START, &
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,       &
          UTAU_LEVEL_MASK_UP, PARTLAYERS_LAYERIDX,       &
          IBEAM, FLUX_MULTIPLIER, LSSL_BOA_SOURCE,       &
          T_DELT_USERM, T_UTUP_USERM,                    &
          K_REAL, K_COMPLEX, UHOM_UPDN, UHOM_UPUP,       &
          HMULT_1, HMULT_2, UT_HMULT_UU, UT_HMULT_UD,    &
          NCON_SLEAVE, PCON_SLEAVE,                      &
          SURFACEWF_F )

!mick temp fix 9/7/2012 - ELSE added
      ELSE
        IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
          SURFACEWF_F(N_REFLEC_WFS+1:N_REFLEC_WFS+N_SLEAVE_WFS,&
                      1:N_USER_LEVELS,1:N_USER_STREAMS,&
                      IBEAM,1:NSTOKES,UPIDX) = ZERO
        ELSE
          SURFACEWF_F(N_REFLEC_WFS+1:N_REFLEC_WFS+N_SLEAVE_WFS,&
                      1:N_USER_LEVELS,LUM,&
                      IBEAM,1:NSTOKES,UPIDX) = ZERO
        ENDIF
      ENDIF

!  Downwelling Albedo weighting functions
!  --------------------------------------

      IF ( DO_DNWELLING ) THEN
        CALL LSSL_DNUSER_SURFACEWF ( &
          DO_OBSERVATION_GEOMETRY,                       &
          N_SLEAVE_WFS, N_REFLEC_WFS, NSTOKES,           &
          N_USER_LEVELS, N_USER_STREAMS, LOCAL_UM_START, &
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,       &
          UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX,       &
          IBEAM, FLUX_MULTIPLIER,                        &
          T_DELT_USERM, T_UTDN_USERM,                    &
          K_REAL, K_COMPLEX, UHOM_DNDN, UHOM_DNUP,       &
          HMULT_1, HMULT_2, UT_HMULT_DU, UT_HMULT_DD,    &
          NCON_SLEAVE, PCON_SLEAVE,                      &
          SURFACEWF_F )
!mick temp fix 9/7/2012 - ELSE added
      ELSE
        IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
          SURFACEWF_F(N_REFLEC_WFS+1:N_REFLEC_WFS+N_SLEAVE_WFS,&
                      1:N_USER_LEVELS,1:N_USER_STREAMS,&
                      IBEAM,1:NSTOKES,DNIDX) = ZERO
        ELSE
          SURFACEWF_F(N_REFLEC_WFS+1:N_REFLEC_WFS+N_SLEAVE_WFS,&
                      1:N_USER_LEVELS,LUM,&
                      IBEAM,1:NSTOKES,DNIDX) = ZERO
        ENDIF
      ENDIF

!  mean value output
!  -----------------

      IF ( DO_INCLUDE_MVOUTPUT.OR.DO_QUAD_OUTPUT ) THEN
        CALL LSSL_INTEGRATED_OUTPUT ( &
          DO_INCLUDE_MVOUTPUT, N_SLEAVE_WFS, N_REFLEC_WFS, &
          NSTOKES, NSTREAMS, NLAYERS, N_USER_LEVELS,       &
          N_DIRECTIONS, WHICH_DIRECTIONS, IBEAM,           &
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,         &
          UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,          &
          PARTLAYERS_LAYERIDX,                             &
          FLUX_MULTIPLIER, QUAD_WEIGHTS, QUAD_STRMWTS,     &
          T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,        &
          K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,         &
          NCON_SLEAVE, PCON_SLEAVE,                        &
          MINT_SURFACEWF, FLUX_SURFACEWF )
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_LSSL_WFS

!

      SUBROUTINE LSSL_UPUSER_SURFACEWF ( &
        DO_OBSERVATION_GEOMETRY,                       &
        N_SLEAVE_WFS, N_REFLEC_WFS, NSTOKES, NLAYERS,  &
        N_USER_LEVELS, N_USER_STREAMS, LOCAL_UM_START, &
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,       &
        UTAU_LEVEL_MASK_UP, PARTLAYERS_LAYERIDX,       &
        IBEAM, FLUX_MULTIPLIER, LSSL_BOA_SOURCE,       &
        T_DELT_USERM, T_UTUP_USERM,                    &
        K_REAL, K_COMPLEX, UHOM_UPDN, UHOM_UPUP,       &
        HMULT_1, HMULT_2, UT_HMULT_UU, UT_HMULT_UD,    &
        NCON_SLEAVE, PCON_SLEAVE,                      &
        SURFACEWF_F )

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Input control

      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY

      INTEGER, INTENT (IN) ::          N_REFLEC_WFS
      INTEGER, INTENT (IN) ::          N_SLEAVE_WFS
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NLAYERS

      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START

      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

!  Solution inputs

      INTEGER, INTENT (IN) ::          IBEAM
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER
      DOUBLE PRECISION, INTENT (IN) :: LSSL_BOA_SOURCE &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAXSTOKES )

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )

      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPDN &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPUP &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: HMULT_1 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_2 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_UU &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_UD &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: NCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )

!  output

      DOUBLE PRECISION, INTENT (INOUT) :: SURFACEWF_F &
         ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_VZANGLES, &
           MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  local variables
!  ---------------

      INTEGER ::          N, NUT, NSTART, NUT_PREV, NLEVEL, O1
      INTEGER ::          UTA, UM, Q, Q1, UT, IB, LUM

      DOUBLE PRECISION :: LSSL_CUMUL_SOURCE &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION :: LSSL_LAYER_SOURCE &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION :: LSSL_FINAL_SOURCE

!  index

      IB = IBEAM
      LUM = 1

!  Zero all Fourier components - New rule, better for safety
!    Only did this for components close to zenith (formerly)

      DO UTA = 1, N_USER_LEVELS
        IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, N_SLEAVE_WFS
              Q1 = Q + N_REFLEC_WFS
              DO O1 = 1, NSTOKES
                SURFACEWF_F(Q1,UTA,UM,IB,O1,UPIDX) = ZERO
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO Q = 1, N_SLEAVE_WFS
            Q1 = Q + N_REFLEC_WFS
            DO O1 = 1, NSTOKES
              SURFACEWF_F(Q1,UTA,LUM,IB,O1,UPIDX) = ZERO
            ENDDO
          ENDDO
        ENDIF
      ENDDO

!  Initialize post-processing recursion
!  ====================================

!  Set the cumulative source term equal to the BOA sum

      IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO Q = 1, N_SLEAVE_WFS
            DO O1 = 1, NSTOKES
              LSSL_CUMUL_SOURCE(Q,UM,O1) = LSSL_BOA_SOURCE(Q,UM,O1)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO Q = 1, N_SLEAVE_WFS
          DO O1 = 1, NSTOKES
            LSSL_CUMUL_SOURCE(Q,IB,O1) = LSSL_BOA_SOURCE(Q,IB,O1)
          ENDDO
        ENDDO
      ENDIF

!  Recursion Loop for linearized Post-processing
!  =============================================

!  initialise cumulative source term loop

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

        NUT = NLEVEL + 1
        DO N = NSTART, NUT, -1
          CALL LSSL_WHOLELAYER_STERM_UP ( &
              DO_OBSERVATION_GEOMETRY,                     &
              N_SLEAVE_WFS, NSTOKES, IB, N,                &
              N_USER_STREAMS, LOCAL_UM_START,              &
              K_REAL, K_COMPLEX, NCON_SLEAVE, PCON_SLEAVE, &
              UHOM_UPDN, UHOM_UPUP, HMULT_1, HMULT_2,      &
              LSSL_LAYER_SOURCE )

          IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO Q = 1, N_SLEAVE_WFS
                DO O1 = 1, NSTOKES
                  LSSL_CUMUL_SOURCE(Q,UM,O1) = LSSL_LAYER_SOURCE(Q,UM,O1) &
                       + T_DELT_USERM(N,UM) * LSSL_CUMUL_SOURCE(Q,UM,O1)
                ENDDO
              ENDDO
            ENDDO
          ELSE
            DO Q = 1, N_SLEAVE_WFS
              DO O1 = 1, NSTOKES
                LSSL_CUMUL_SOURCE(Q,IB,O1) = LSSL_LAYER_SOURCE(Q,IB,O1) &
                     + T_DELT_USERM(N,IB) * LSSL_CUMUL_SOURCE(Q,IB,O1)
              ENDDO
            ENDDO
          ENDIF
        ENDDO

!  Offgrid output
!  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)

          CALL LSSL_PARTLAYER_STERM_UP ( &
              DO_OBSERVATION_GEOMETRY,                        &
              N_SLEAVE_WFS, NSTOKES, IB, N,                   &
              UT, N_USER_STREAMS, LOCAL_UM_START,             &
              K_REAL, K_COMPLEX, NCON_SLEAVE, PCON_SLEAVE,    &
              UHOM_UPDN, UHOM_UPUP, UT_HMULT_UU, UT_HMULT_UD, &
              LSSL_LAYER_SOURCE )

          IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO Q = 1, N_SLEAVE_WFS
                Q1 = Q + N_REFLEC_WFS
                DO O1 = 1, NSTOKES
                  LSSL_FINAL_SOURCE = LSSL_LAYER_SOURCE(Q,UM,O1) &
                  + T_UTUP_USERM(UT,UM) * LSSL_CUMUL_SOURCE(Q,UM,O1)
                  SURFACEWF_F(Q1,UTA,UM,IB,O1,UPIDX) = &
                           FLUX_MULTIPLIER * LSSL_FINAL_SOURCE
                ENDDO
              ENDDO
            ENDDO
          ELSE
            DO Q = 1, N_SLEAVE_WFS
              Q1 = Q + N_REFLEC_WFS
              DO O1 = 1, NSTOKES
                LSSL_FINAL_SOURCE = LSSL_LAYER_SOURCE(Q,IB,O1) &
                + T_UTUP_USERM(UT,IB) * LSSL_CUMUL_SOURCE(Q,IB,O1)
                SURFACEWF_F(Q1,UTA,LUM,IB,O1,UPIDX) = &
                         FLUX_MULTIPLIER * LSSL_FINAL_SOURCE
              ENDDO
            ENDDO
          ENDIF

!  Ongrid output
!  -------------

        ELSE
          IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO Q = 1, N_SLEAVE_WFS
                Q1 = Q + N_REFLEC_WFS
                DO O1 = 1, NSTOKES
                  SURFACEWF_F(Q1,UTA,UM,IB,O1,UPIDX) = &
                             FLUX_MULTIPLIER * LSSL_CUMUL_SOURCE(Q,UM,O1)
                ENDDO
              ENDDO
            ENDDO
          ELSE
            DO Q = 1, N_SLEAVE_WFS
              Q1 = Q + N_REFLEC_WFS
              DO O1 = 1, NSTOKES
                SURFACEWF_F(Q1,UTA,LUM,IB,O1,UPIDX) = &
                           FLUX_MULTIPLIER * LSSL_CUMUL_SOURCE(Q,IB,O1)
              ENDDO
            ENDDO
          ENDIF
        ENDIF

!  Check for updating the recursion

        IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
        NUT_PREV = NUT

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSSL_UPUSER_SURFACEWF

!

      SUBROUTINE LSSL_DNUSER_SURFACEWF ( &
        DO_OBSERVATION_GEOMETRY,                       &
        N_SLEAVE_WFS, N_REFLEC_WFS, NSTOKES,           &
        N_USER_LEVELS, N_USER_STREAMS, LOCAL_UM_START, &
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,       &
        UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX,       &
        IBEAM, FLUX_MULTIPLIER,                        & 
        T_DELT_USERM, T_UTDN_USERM,                    &
        K_REAL, K_COMPLEX, UHOM_DNDN, UHOM_DNUP,       &
        HMULT_1, HMULT_2, UT_HMULT_DU, UT_HMULT_DD,    &
        NCON_SLEAVE, PCON_SLEAVE,                      &
        SURFACEWF_F )

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Input control

      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY

      INTEGER, INTENT (IN) ::          N_REFLEC_WFS
      INTEGER, INTENT (IN) ::          N_SLEAVE_WFS
      INTEGER, INTENT (IN) ::          NSTOKES

      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START

      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

!  Solution input variables

      INTEGER, INTENT (IN) ::          IBEAM
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )

      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNDN &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNUP &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: HMULT_1 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_2 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: NCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_DU &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_DD &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )

!  output

      DOUBLE PRECISION, INTENT (INOUT) :: SURFACEWF_F &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_VZANGLES, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  local variables
!  ---------------

      INTEGER ::          N, NUT, NSTART, NUT_PREV, NLEVEL, O1
      INTEGER ::          UTA, UM, Q, Q1,UT, IB, LUM

      DOUBLE PRECISION :: LSSL_CUMUL_SOURCE &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION :: LSSL_LAYER_SOURCE &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION :: LSSL_TOA_SOURCE &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION :: LSSL_FINAL_SOURCE

!  Initialise

      IB  = IBEAM
      LUM = 1

!  Zero all Fourier component output

      DO UTA = 1, N_USER_LEVELS
        IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
!mick fix 9/17/2012 - changed limit of do loop
          !DO UM = 1, LOCAL_UM_START
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, N_SLEAVE_WFS
              Q1 = Q + N_REFLEC_WFS
               DO O1 = 1, NSTOKES
                  SURFACEWF_F(Q1,UTA,UM,IB,O1,DNIDX) = ZERO
               ENDDO
            ENDDO
          ENDDO
        ELSE
          DO Q = 1, N_SLEAVE_WFS
            Q1 = Q + N_REFLEC_WFS
             DO O1 = 1, NSTOKES
                SURFACEWF_F(Q1,UTA,LUM,IB,O1,DNIDX) = ZERO
             ENDDO
          ENDDO
        ENDIF
      ENDDO

!  Initialize post-processing recursion
!  ====================================

!  Get the linearized TOA source terms

      IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO Q = 1, N_SLEAVE_WFS
            DO O1 = 1, NSTOKES
              LSSL_TOA_SOURCE(Q,UM,O1)   = ZERO
              LSSL_CUMUL_SOURCE(Q,UM,O1) = LSSL_TOA_SOURCE(Q,UM,O1)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO Q = 1, N_SLEAVE_WFS
          DO O1 = 1, NSTOKES
            LSSL_TOA_SOURCE(Q,IB,O1)   = ZERO
            LSSL_CUMUL_SOURCE(Q,IB,O1) = LSSL_TOA_SOURCE(Q,IB,O1)
          ENDDO
        ENDDO
      ENDIF

!  Recursion Loop for linearized Post-processing
!  =============================================

!  initialise cumulative source term loop

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

        NUT = NLEVEL
        DO N = NSTART, NUT

          CALL LSSL_WHOLELAYER_STERM_DN ( &
              DO_OBSERVATION_GEOMETRY,                     &
              N_SLEAVE_WFS, NSTOKES, IB, N,                &
              N_USER_STREAMS, LOCAL_UM_START,              &
              K_REAL, K_COMPLEX, NCON_SLEAVE, PCON_SLEAVE, &
              UHOM_DNDN, UHOM_DNUP, HMULT_1, HMULT_2,      &
              LSSL_LAYER_SOURCE )

          IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO Q = 1, N_SLEAVE_WFS
                DO O1 = 1, NSTOKES
                  LSSL_CUMUL_SOURCE(Q,UM,O1) = LSSL_LAYER_SOURCE(Q,UM,O1) &
                       + T_DELT_USERM(N,UM) * LSSL_CUMUL_SOURCE(Q,UM,O1)
                ENDDO
              ENDDO
            ENDDO
          ELSE
            DO Q = 1, N_SLEAVE_WFS
              DO O1 = 1, NSTOKES
                LSSL_CUMUL_SOURCE(Q,IB,O1) = LSSL_LAYER_SOURCE(Q,IB,O1) &
                     + T_DELT_USERM(N,IB) * LSSL_CUMUL_SOURCE(Q,IB,O1)
              ENDDO
            ENDDO
          ENDIF

        ENDDO

!  Offgrid output
!  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)

!  User-defined stream output, add additional partial layer source term

          CALL LSSL_PARTLAYER_STERM_DN ( &
              DO_OBSERVATION_GEOMETRY,                        &
              N_SLEAVE_WFS, NSTOKES, IB, N,                   &
              UT, N_USER_STREAMS, LOCAL_UM_START,             &
              K_REAL, K_COMPLEX, NCON_SLEAVE, PCON_SLEAVE,    &
              UHOM_DNDN, UHOM_DNUP, UT_HMULT_DU, UT_HMULT_DD, &
              LSSL_LAYER_SOURCE )

          IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO Q = 1, N_SLEAVE_WFS
                Q1 = Q + N_REFLEC_WFS
                DO O1 = 1, NSTOKES
                  LSSL_FINAL_SOURCE = LSSL_LAYER_SOURCE(Q,UM,O1) &
                    + T_UTDN_USERM(UT,UM) * LSSL_CUMUL_SOURCE(Q,UM,O1)
                  SURFACEWF_F(Q1,UTA,UM,IB,O1,DNIDX) = &
                    FLUX_MULTIPLIER * LSSL_FINAL_SOURCE
                ENDDO
              ENDDO
            ENDDO
          ELSE
            DO Q = 1, N_SLEAVE_WFS
              Q1 = Q + N_REFLEC_WFS
              DO O1 = 1, NSTOKES
                LSSL_FINAL_SOURCE = LSSL_LAYER_SOURCE(Q,IB,O1) &
                  + T_UTDN_USERM(UT,IB) * LSSL_CUMUL_SOURCE(Q,IB,O1)
                SURFACEWF_F(Q1,UTA,LUM,IB,O1,DNIDX) = &
                  FLUX_MULTIPLIER * LSSL_FINAL_SOURCE
              ENDDO
            ENDDO
          ENDIF

!  Ongrid output
!  -------------

        ELSE

!  User-defined stream output, just set to the cumulative source term

          IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO Q = 1, N_SLEAVE_WFS
                Q1 = Q + N_REFLEC_WFS
                DO O1 = 1, NSTOKES
                  SURFACEWF_F(Q1,UTA,UM,IB,O1,DNIDX) = &
                            FLUX_MULTIPLIER * LSSL_CUMUL_SOURCE(Q,UM,O1)
                ENDDO
              ENDDO
            ENDDO
          ELSE
            DO Q = 1, N_SLEAVE_WFS
              Q1 = Q + N_REFLEC_WFS
              DO O1 = 1, NSTOKES
                SURFACEWF_F(Q1,UTA,LUM,IB,O1,DNIDX) = &
                          FLUX_MULTIPLIER * LSSL_CUMUL_SOURCE(Q,IB,O1)
              ENDDO
            ENDDO
          ENDIF

        ENDIF

!  Check for updating the recursion

        IF ( NUT.NE. NUT_PREV ) NSTART = NUT + 1
        NUT_PREV = NUT

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSSL_DNUSER_SURFACEWF

!

      SUBROUTINE LSSL_WHOLELAYER_STERM_UP ( &
        DO_OBSERVATION_GEOMETRY,                     &
        N_SLEAVE_WFS, NSTOKES, IBEAM, GIVEN_LAYER,   &
        N_USER_STREAMS, LOCAL_UM_START,              &
        K_REAL, K_COMPLEX, NCON_SLEAVE, PCON_SLEAVE, &
        UHOM_UPDN, UHOM_UPUP, HMULT_1, HMULT_2,      &
        LSSL_LAYERSOURCE )

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Control inputs

      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY

      INTEGER, INTENT (IN) ::          N_SLEAVE_WFS
      INTEGER, INTENT (IN) ::          IBEAM, GIVEN_LAYER
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START

!  Solution inputs

      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPDN &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPUP &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_1 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_2 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: NCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )

!  outputs

      DOUBLE PRECISION, INTENT (OUT) :: LSSL_LAYERSOURCE &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          N, UM, O1, IB, K, KO1, K0, K1, K2, Q
      INTEGER ::          UM_START, UM_END
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, H1, H2
      DOUBLE PRECISION :: NUXR, PUXR, NUXR1, NUXR2, PUXR1, PUXR2


!  local indices

      N   = GIVEN_LAYER
      KO1 = K_REAL(N) + 1
      IB  = IBEAM

      IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
        UM_START = LOCAL_UM_START
        UM_END   = N_USER_STREAMS
      ELSE
        UM_START = IB
        UM_END   = IB
      ENDIF

!  Homogeneous solutions
!  =====================

!  Loops over user angles, weighting functions and Stokes

      DO UM = UM_START,UM_END
       DO Q = 1, N_SLEAVE_WFS
        DO O1 = 1, NSTOKES

!  Real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, K_REAL(N)
            NUXR = NCON_SLEAVE(Q,K,N)*UHOM_UPDN(UM,O1,K,N)
            PUXR = PCON_SLEAVE(Q,K,N)*UHOM_UPUP(UM,O1,K,N)
            H1 =  NUXR * HMULT_2(K,UM,N)
            H2 =  PUXR * HMULT_1(K,UM,N)
            SHOM_R = SHOM_R + H1 + H2
          ENDDO

!  Complex homogeneous solutions

          SHOM_CR = ZERO
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            NUXR1 =   NCON_SLEAVE(Q,K1,N) * UHOM_UPDN(UM,O1,K1,N) &
                    - NCON_SLEAVE(Q,K2,N) * UHOM_UPDN(UM,O1,K2,N)
            NUXR2 =   NCON_SLEAVE(Q,K1,N) * UHOM_UPDN(UM,O1,K2,N) &
                    + NCON_SLEAVE(Q,K2,N) * UHOM_UPDN(UM,O1,K1,N)
            PUXR1 =   PCON_SLEAVE(Q,K1,N) * UHOM_UPUP(UM,O1,K1,N) &
                    - PCON_SLEAVE(Q,K2,N) * UHOM_UPUP(UM,O1,K2,N)
            PUXR2 =   PCON_SLEAVE(Q,K1,N) * UHOM_UPUP(UM,O1,K2,N) &
                    + PCON_SLEAVE(Q,K2,N) * UHOM_UPUP(UM,O1,K1,N)
            H1 =   NUXR1 * HMULT_2(K1,UM,N) &
                 - NUXR2 * HMULT_2(K2,UM,N)
            H2 =   PUXR1 * HMULT_1(K1,UM,N) &
                 - PUXR2 * HMULT_1(K2,UM,N)
            SHOM_CR = SHOM_CR + H1 + H2
          ENDDO

!  homogeneous contribution

          LSSL_LAYERSOURCE(Q,UM,O1) = SHOM_R + SHOM_CR

!  End loops over Q, O1 and UM

        ENDDO
       ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSSL_WHOLELAYER_STERM_UP

!

      SUBROUTINE LSSL_WHOLELAYER_STERM_DN ( &
        DO_OBSERVATION_GEOMETRY,                     &
        N_SLEAVE_WFS, NSTOKES, IBEAM, GIVEN_LAYER,   &
        N_USER_STREAMS, LOCAL_UM_START,              &
        K_REAL, K_COMPLEX, NCON_SLEAVE, PCON_SLEAVE, &
        UHOM_DNDN, UHOM_DNUP, HMULT_1, HMULT_2,      &
        LSSL_LAYERSOURCE )

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Control inputs

      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY

      INTEGER, INTENT (IN) ::          N_SLEAVE_WFS
      INTEGER, INTENT (IN) ::          IBEAM, GIVEN_LAYER
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START

!  Solution inputs

      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNDN &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNUP &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_1 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_2 &
          ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: NCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )

!  output

      DOUBLE PRECISION, INTENT (OUT) :: LSSL_LAYERSOURCE &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          N, UM, O1, IB, K, KO1, K0, K1, K2, Q
      INTEGER ::          UM_START, UM_END
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, H1, H2
      DOUBLE PRECISION :: NUXR, PUXR, NUXR1, NUXR2, PUXR1, PUXR2

!  local indices

      N   = GIVEN_LAYER
      KO1 = K_REAL(N) + 1
      IB  = IBEAM

      IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
        UM_START = LOCAL_UM_START
        UM_END   = N_USER_STREAMS
      ELSE
        UM_START = IB
        UM_END   = IB
      ENDIF

!  Homogeneous solutions
!  =====================

!  Loop over user angles, weighting functions and STokes

      DO UM = UM_START,UM_END
       DO Q = 1, N_SLEAVE_WFS
        DO O1 = 1, NSTOKES

!  Real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, K_REAL(N)
            NUXR = NCON_SLEAVE(Q,K,N)*UHOM_DNDN(UM,O1,K,N)
            PUXR = PCON_SLEAVE(Q,K,N)*UHOM_DNUP(UM,O1,K,N)
            H1 =  NUXR * HMULT_1(K,UM,N)
            H2 =  PUXR * HMULT_2(K,UM,N)
            SHOM_R = SHOM_R + H1 + H2
          ENDDO

!  Complex homogeneous solutions

          SHOM_CR = ZERO
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            NUXR1 =   NCON_SLEAVE(Q,K1,N) * UHOM_DNDN(UM,O1,K1,N) &
                    - NCON_SLEAVE(Q,K2,N) * UHOM_DNDN(UM,O1,K2,N)
            NUXR2 =   NCON_SLEAVE(Q,K1,N) * UHOM_DNDN(UM,O1,K2,N) &
                    + NCON_SLEAVE(Q,K2,N) * UHOM_DNDN(UM,O1,K1,N)
            PUXR1 =   PCON_SLEAVE(Q,K1,N) * UHOM_DNUP(UM,O1,K1,N) &
                    - PCON_SLEAVE(Q,K2,N) * UHOM_DNUP(UM,O1,K2,N)
            PUXR2 =   PCON_SLEAVE(Q,K1,N) * UHOM_DNUP(UM,O1,K2,N) &
                    + PCON_SLEAVE(Q,K2,N) * UHOM_DNUP(UM,O1,K1,N)
            H1 =   NUXR1 * HMULT_1(K1,UM,N) &
                 - NUXR2 * HMULT_1(K2,UM,N)
            H2 =   PUXR1 * HMULT_2(K1,UM,N) &
                 - PUXR2 * HMULT_2(K2,UM,N)
            SHOM_CR = SHOM_CR + H1 + H2
          ENDDO

!  homogeneous contribution

          LSSL_LAYERSOURCE(Q,UM,O1) = SHOM_R + SHOM_CR

!  End loops over Q, O1 and UM

        ENDDO
       ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSSL_WHOLELAYER_STERM_DN

!

      SUBROUTINE LSSL_PARTLAYER_STERM_UP ( &
        DO_OBSERVATION_GEOMETRY, &
        N_SLEAVE_WFS, NSTOKES, IBEAM, GIVEN_LAYER,      &
        OFFGRID_INDEX, N_USER_STREAMS, LOCAL_UM_START,  &
        K_REAL, K_COMPLEX, NCON_SLEAVE, PCON_SLEAVE,    &
        UHOM_UPDN, UHOM_UPUP, UT_HMULT_UU, UT_HMULT_UD, &
        LSSL_LAYERSOURCE )

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Control inputs

      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY

      INTEGER, INTENT (IN) ::          N_SLEAVE_WFS
      INTEGER, INTENT (IN) ::          GIVEN_LAYER, IBEAM
      INTEGER, INTENT (IN) ::          OFFGRID_INDEX
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START

!  Solution input variables

      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPDN &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPUP &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_UU &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_UD &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: NCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )

!  output

      DOUBLE PRECISION, INTENT (OUT) :: LSSL_LAYERSOURCE &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          N, UM, O1, IB, UT, K, KO1, K0, K1, K2, Q
      INTEGER ::          UM_START, UM_END
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, H1, H2
      DOUBLE PRECISION :: NUXR, PUXR, NUXR1, NUXR2, PUXR1, PUXR2

!  local indices

      N   = GIVEN_LAYER
      KO1 = K_REAL(N) + 1
      UT  = OFFGRID_INDEX
      IB  = IBEAM

      IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
        UM_START = LOCAL_UM_START
        UM_END   = N_USER_STREAMS
      ELSE
        UM_START = IB
        UM_END   = IB
      ENDIF

!  Partial layer source function ( Homogeneous/constants variation )
!  =================================================================

!  Loop over user angles, weighting functions and STokes

      DO UM = UM_START,UM_END
       DO Q = 1, N_SLEAVE_WFS
        DO O1 = 1, NSTOKES

!  Real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, K_REAL(N)
            NUXR = NCON_SLEAVE(Q,K,N) * UHOM_UPDN(UM,O1,K,N)
            PUXR = PCON_SLEAVE(Q,K,N) * UHOM_UPUP(UM,O1,K,N)
            H1 =  NUXR * UT_HMULT_UD(K,UM,UT)
            H2 =  PUXR * UT_HMULT_UU(K,UM,UT)
            SHOM_R = SHOM_R + H1 + H2
          ENDDO

!  Complex homogeneous solutions
!     Rob, 12/21/10, Changed N -> UT in UT_HMULT_DD, UT_HMULT_DU

          SHOM_CR = ZERO
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            NUXR1 =   NCON_SLEAVE(Q,K1,N) * UHOM_UPDN(UM,O1,K1,N) &
                    - NCON_SLEAVE(Q,K2,N) * UHOM_UPDN(UM,O1,K2,N)
            NUXR2 =   NCON_SLEAVE(Q,K1,N) * UHOM_UPDN(UM,O1,K2,N) &
                    + NCON_SLEAVE(Q,K2,N) * UHOM_UPDN(UM,O1,K1,N)
            PUXR1 =   PCON_SLEAVE(Q,K1,N) * UHOM_UPUP(UM,O1,K1,N) &
                    - PCON_SLEAVE(Q,K2,N) * UHOM_UPUP(UM,O1,K2,N)
            PUXR2 =   PCON_SLEAVE(Q,K1,N) * UHOM_UPUP(UM,O1,K2,N) &
                    + PCON_SLEAVE(Q,K2,N) * UHOM_UPUP(UM,O1,K1,N)
            H1 =   NUXR1 * UT_HMULT_UD(K1,UM,UT) &
                 - NUXR2 * UT_HMULT_UD(K2,UM,UT)
            H2 =   PUXR1 * UT_HMULT_UU(K1,UM,UT) &
                 - PUXR2 * UT_HMULT_UU(K2,UM,UT)
            SHOM_CR = SHOM_CR + H1 + H2
          ENDDO

!  homogeneous contribution

          LSSL_LAYERSOURCE(Q,UM,O1) = SHOM_R + SHOM_CR

!  End loops over Q, O1 and UM

        ENDDO
       ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSSL_PARTLAYER_STERM_UP

!

      SUBROUTINE LSSL_PARTLAYER_STERM_DN ( &
        DO_OBSERVATION_GEOMETRY, &
        N_SLEAVE_WFS, NSTOKES, IBEAM, GIVEN_LAYER,      &
        OFFGRID_INDEX, N_USER_STREAMS, LOCAL_UM_START,  &
        K_REAL, K_COMPLEX, NCON_SLEAVE, PCON_SLEAVE,    &
        UHOM_DNDN, UHOM_DNUP, UT_HMULT_DU, UT_HMULT_DD, &
        LSSL_LAYERSOURCE )

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Control inputs

      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY

      INTEGER, INTENT (IN) ::          N_SLEAVE_WFS
      INTEGER, INTENT (IN) ::          GIVEN_LAYER, IBEAM
      INTEGER, INTENT (IN) ::          OFFGRID_INDEX
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START

!  Solution input variables

      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNDN &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNUP &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_DU &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_DD &
          ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: NCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )

!  output

      DOUBLE PRECISION, INTENT (OUT) :: LSSL_LAYERSOURCE &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          N, UM, O1, IB, UT, K, KO1, K0, K1, K2, Q
      INTEGER ::          UM_START, UM_END
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, H1, H2
      DOUBLE PRECISION :: NUXR, PUXR, NUXR1, NUXR2, PUXR1, PUXR2

!  local indices

      N   = GIVEN_LAYER
      KO1 = K_REAL(N) + 1
      UT  = OFFGRID_INDEX
      IB  = IBEAM

      IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
        UM_START = LOCAL_UM_START
        UM_END   = N_USER_STREAMS
      ELSE
        UM_START = IB
        UM_END   = IB
      ENDIF

!  Partial layer source function ( Homogeneous/constants variation )
!  =================================================================

!  Loop over user angles, weighting functions and STokes

      DO UM = UM_START,UM_END
       DO Q = 1, N_SLEAVE_WFS
        DO O1 = 1, NSTOKES

!  Real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, K_REAL(N)
            NUXR = NCON_SLEAVE(Q,K,N) * UHOM_DNDN(UM,O1,K,N)
            PUXR = PCON_SLEAVE(Q,K,N) * UHOM_DNUP(UM,O1,K,N)
            H1 =  NUXR * UT_HMULT_DD(K,UM,UT)
            H2 =  PUXR * UT_HMULT_DU(K,UM,UT)
            SHOM_R = SHOM_R + H1 + H2
          ENDDO

!  Complex homogeneous solutions
!     Rob, 12/21/10, Changed N -> UT in UT_HMULT_DD, UT_HMULT_DU

          SHOM_CR = ZERO
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            NUXR1 =   NCON_SLEAVE(Q,K1,N) * UHOM_DNDN(UM,O1,K1,N) &
                    - NCON_SLEAVE(Q,K2,N) * UHOM_DNDN(UM,O1,K2,N)
            NUXR2 =   NCON_SLEAVE(Q,K1,N) * UHOM_DNDN(UM,O1,K2,N) &
                    + NCON_SLEAVE(Q,K2,N) * UHOM_DNDN(UM,O1,K1,N)
            PUXR1 =   PCON_SLEAVE(Q,K1,N) * UHOM_DNUP(UM,O1,K1,N) &
                    - PCON_SLEAVE(Q,K2,N) * UHOM_DNUP(UM,O1,K2,N)
            PUXR2 =   PCON_SLEAVE(Q,K1,N) * UHOM_DNUP(UM,O1,K2,N) &
                    + PCON_SLEAVE(Q,K2,N) * UHOM_DNUP(UM,O1,K1,N)
            H1 =   NUXR1 * UT_HMULT_DD(K1,UM,UT) &
                 - NUXR2 * UT_HMULT_DD(K2,UM,UT)
            H2 =   PUXR1 * UT_HMULT_DU(K1,UM,UT) &
                 - PUXR2 * UT_HMULT_DU(K2,UM,UT)
            SHOM_CR = SHOM_CR + H1 + H2
          ENDDO

!  homogeneous contribution

          LSSL_LAYERSOURCE(Q,UM,O1) = SHOM_R + SHOM_CR

!  End loops over Q, O1 and UM

        ENDDO
       ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSSL_PARTLAYER_STERM_DN

!

      SUBROUTINE LSSL_INTEGRATED_OUTPUT ( &
        DO_INCLUDE_MVOUT, N_SLEAVE_WFS, N_REFLEC_WFS, &
        NSTOKES, NSTREAMS, NLAYERS, N_USER_LEVELS,    &
        N_DIRECTIONS, WHICH_DIRECTIONS, IBEAM,        &
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,      &
        UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,       &
        PARTLAYERS_LAYERIDX,                          &
        FLUX_MULTIPLIER, QUAD_WEIGHTS, QUAD_STRMWTS,  &
        T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,     &
        K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,      &
        NCON_SLEAVE, PCON_SLEAVE,                     &
        MINT_SURFACEWF, FLUX_SURFACEWF )

!  Quadrature output at offgrid or ongrid optical depths
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

      USE VLIDORT_PARS

      IMPLICIT NONE

!  INPUT
!  =====

!  Top level flag

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_MVOUT

!  Basic munbers

      INTEGER, INTENT (IN) ::          IBEAM
      INTEGER, INTENT (IN) ::          N_SLEAVE_WFS, N_REFLEC_WFS
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS

!  Directional and Level output controls

      INTEGER, INTENT (IN) ::          N_DIRECTIONS
      INTEGER, INTENT (IN) ::          WHICH_DIRECTIONS ( MAX_DIRECTIONS )

      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

!  stream and multiplier input

      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER
      DOUBLE PRECISION, INTENT (IN) :: QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )

!  Solution variables

      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )

      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: NCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )

!  output
!  ======

      DOUBLE PRECISION, INTENT (INOUT) :: MINT_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) :: FLUX_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  Local variables
!  ===============

!  Local quadrature output

      DOUBLE PRECISION :: QSLEAVEWF_F &
          ( MAX_SLEAVEWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  Help variables

      INTEGER ::          I, IDIR, WDIR, UTA, UT, N, NLEVEL, O1, Q, Q1
      DOUBLE PRECISION :: SMI, SFX

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
              CALL LSSL_QUADWF_OFFGRID_UP ( &
                   UTA, UT, N, N_SLEAVE_WFS, NSTOKES, NSTREAMS, &
                   FLUX_MULTIPLIER, T_UTUP_EIGEN, T_UTDN_EIGEN, & 
                   K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,     &
                   NCON_SLEAVE, PCON_SLEAVE, QSLEAVEWF_F )
            ELSE
              CALL LSSL_QUADWF_LEVEL_UP  (  &
                 UTA, NLEVEL, N_SLEAVE_WFS,                &
                 NSTOKES, NSTREAMS, NLAYERS,               &
                 FLUX_MULTIPLIER, T_DELT_EIGEN,            &
                 K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,  &
                 NCON_SLEAVE, PCON_SLEAVE, QSLEAVEWF_F )
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
              CALL LSSL_QUADWF_OFFGRID_DN ( &
                   UTA, UT, N, N_SLEAVE_WFS, NSTOKES, NSTREAMS, &
                   FLUX_MULTIPLIER, T_UTUP_EIGEN, T_UTDN_EIGEN, & 
                   K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,     &
                   NCON_SLEAVE, PCON_SLEAVE, QSLEAVEWF_F )
            ELSE
              CALL LSSL_QUADWF_LEVEL_DN  (  &
                 UTA, NLEVEL, N_SLEAVE_WFS,                &
                 NSTOKES, NSTREAMS,                        &
                 FLUX_MULTIPLIER, T_DELT_EIGEN,            &
                 K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,  &
                 NCON_SLEAVE, PCON_SLEAVE, QSLEAVEWF_F )
            ENDIF
          ENDDO
        ENDIF

!  Mean Intensity and Flux output (diffuse term only). If flagged

        IF ( DO_INCLUDE_MVOUT ) THEN
          DO UTA = 1, N_USER_LEVELS
           DO Q = 1, N_SLEAVE_WFS
            Q1 = Q + N_REFLEC_WFS
            DO O1 = 1, NSTOKES
              SMI = ZERO
              SFX = ZERO
              DO I = 1, NSTREAMS
                SMI = SMI + QUAD_WEIGHTS(I) * QSLEAVEWF_F(Q,UTA,I,O1)
                SFX = SFX + QUAD_STRMWTS(I) * QSLEAVEWF_F(Q,UTA,I,O1)
              ENDDO
              MINT_SURFACEWF(Q1,UTA,IBEAM,O1,WDIR) = SMI * HALF
              FLUX_SURFACEWF(Q1,UTA,IBEAM,O1,WDIR) = SFX * PI2
            ENDDO
           ENDDO
          ENDDO
        ENDIF

!  End directions loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSSL_INTEGRATED_OUTPUT

!

      SUBROUTINE LSSL_QUADWF_LEVEL_UP (                                     &
        UTA, NL, N_SLEAVE_WFS, NSTOKES, NSTREAMS, NLAYERS, FLUX_MULTIPLIER, &
        T_DELT_EIGEN, K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,              &
        NCON_SLEAVE, PCON_SLEAVE, QSURFACEWF_F )

!  Upwelling weighting function Fourier components at level boundary NL
!  Quadrature angles only

      USE VLIDORT_PARS

      IMPLICIT NONE

!  input

      INTEGER, INTENT (IN) ::          UTA, NL, N_SLEAVE_WFS
      INTEGER, INTENT (IN) ::          NSTOKES, NSTREAMS, NLAYERS
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: NCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )

!  output

      DOUBLE PRECISION, INTENT (INOUT) :: QSURFACEWF_F &
          ( MAX_SLEAVEWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          N, I, I1, O1, K, KO1, K0, K1, K2, Q
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, H1, H2
      DOUBLE PRECISION :: NXR, NXR1, NXR2, PXR, PXR1, PXR2

!  This depends on the level mask - if this is 0 to NLAYERS - 1, then we
!  looking at the perturbation field at the top of these layers. The
!  case where the level mask = NLAYERS is the upwelling perturbed fields
!  at the bottom of the atmosphere (treated separately).

      N = NL + 1

!  For the lowest level
!  ====================

      IF ( NL .EQ. NLAYERS ) THEN

!  Offset

       KO1 = K_REAL(NL) + 1

!  parameter loop, Stokes and streams loops

       DO Q = 1, N_SLEAVE_WFS
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES

!  Real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(NL)
              NXR = NCON_SLEAVE(Q,K,NL) * SOLA_XPOS(I1,O1,K,NL)
              PXR = PCON_SLEAVE(Q,K,NL) * SOLB_XNEG(I1,O1,K,NL)
              SHOM_R = SHOM_R + NXR*T_DELT_EIGEN(K,NL) + PXR
            ENDDO

!  Complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(NL)
             K0 = 2 * K - 2
             K1 = KO1 + K0
             K2 = K1  + 1
             NXR1 =   NCON_SLEAVE(Q,K1,NL) * SOLA_XPOS(I1,O1,K1,NL) &
                    - NCON_SLEAVE(Q,K2,NL) * SOLA_XPOS(I1,O1,K2,NL)
             NXR2 =   NCON_SLEAVE(Q,K1,NL) * SOLA_XPOS(I1,O1,K2,NL) &
                    + NCON_SLEAVE(Q,K2,NL) * SOLA_XPOS(I1,O1,K1,NL)
             PXR1 =   PCON_SLEAVE(Q,K1,NL) * SOLB_XNEG(I1,O1,K1,NL) &
                    - PCON_SLEAVE(Q,K2,NL) * SOLB_XNEG(I1,O1,K2,NL)
             H1 =   NXR1 * T_DELT_EIGEN(K1,NL) &
                  - NXR2 * T_DELT_EIGEN(K2,NL)
             H2 =  PXR1
             SHOM_CR = SHOM_CR + H1 + H2
            ENDDO

!  collect solution

            QSURFACEWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER*(SHOM_R+SHOM_CR)

!  Finish loops

          ENDDO
        ENDDO
       ENDDO

      ENDIF

!  For other levels in the atmosphere
!  ==================================

      IF ( NL .NE. NLAYERS ) THEN

!  Offset

       KO1 = K_REAL(N) + 1

!  parameter loop, Stokes and streams loops

       DO Q = 1, N_SLEAVE_WFS
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES

!  Real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(N)
              NXR = NCON_SLEAVE(Q,K,N) * SOLA_XPOS(I1,O1,K,N)
              PXR = PCON_SLEAVE(Q,K,N) * SOLB_XNEG(I1,O1,K,N)
              SHOM_R = SHOM_R + PXR*T_DELT_EIGEN(K,N) + NXR
            ENDDO

!  Complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(N)
             K0 = 2 * K - 2
             K1 = KO1 + K0
             K2 = K1  + 1
             NXR1 =   NCON_SLEAVE(Q,K1,N) * SOLA_XPOS(I1,O1,K1,N) &
                    - NCON_SLEAVE(Q,K2,N) * SOLA_XPOS(I1,O1,K2,N)
             PXR1 =   PCON_SLEAVE(Q,K1,N) * SOLB_XNEG(I1,O1,K1,N) &
                    - PCON_SLEAVE(Q,K2,N) * SOLB_XNEG(I1,O1,K2,N)
             PXR2 =   PCON_SLEAVE(Q,K1,N) * SOLB_XNEG(I1,O1,K2,N) &
                    + PCON_SLEAVE(Q,K2,N) * SOLB_XNEG(I1,O1,K1,N)
             H1 =   PXR1 * T_DELT_EIGEN(K1,N) &
                  - PXR2 * T_DELT_EIGEN(K2,N)
             H2 =  NXR1
             SHOM_CR = SHOM_CR + H1 + H2
            ENDDO

!  collect solution

           QSURFACEWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER*(SHOM_R+SHOM_CR)

!  Finish loops

          ENDDO
        ENDDO
       ENDDO

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LSSL_QUADWF_LEVEL_UP

!

      SUBROUTINE LSSL_QUADWF_LEVEL_DN (                           &
        UTA, NL, N_SLEAVE_WFS, NSTOKES, NSTREAMS, FLUX_MULTIPLIER,  &
        T_DELT_EIGEN,  K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,     &
        NCON_SLEAVE, PCON_SLEAVE, QSURFACEWF_F )

      USE VLIDORT_PARS

      IMPLICIT NONE

!  input

      INTEGER, INTENT (IN) ::          UTA, NL, N_SLEAVE_WFS
      INTEGER, INTENT (IN) ::          NSTOKES, NSTREAMS
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: NCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )

!  output

      DOUBLE PRECISION, INTENT (INOUT) :: QSURFACEWF_F &
          ( MAX_SLEAVEWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          N, I, O1, K, KO1, K0, K1, K2, Q
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, H1, H2
      DOUBLE PRECISION :: NXR, NXR1, NXR2, PXR, PXR1

!  Downwelling weighting functions at TOA ( or N = 0 ) are zero
!  ============================================================

      IF ( NL .EQ. 0 ) THEN
        DO Q = 1, N_SLEAVE_WFS
          DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES
              QSURFACEWF_F(Q,UTA,I,O1) = ZERO
            ENDDO
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Downwelling weighting functions at other levels
!  ===============================================

!  Offset

      N = NL
      KO1 = K_REAL(N) + 1

!  parameter loop,  Stokes and streams loops

      DO Q = 1, N_SLEAVE_WFS
        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES

!  Real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(N)
              NXR = NCON_SLEAVE(Q,K,N) * SOLA_XPOS(I,O1,K,N)
              PXR = PCON_SLEAVE(Q,K,N) * SOLB_XNEG(I,O1,K,N)
              SHOM_R = SHOM_R + NXR*T_DELT_EIGEN(K,N) + PXR
            ENDDO

!  Complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(N)
             K0 = 2 * K - 2
             K1 = KO1 + K0
             K2 = K1  + 1
             NXR1 =   NCON_SLEAVE(Q,K1,N) * SOLA_XPOS(I,O1,K1,N) &
                    - NCON_SLEAVE(Q,K2,N) * SOLA_XPOS(I,O1,K2,N)
             NXR2 =   NCON_SLEAVE(Q,K1,N) * SOLA_XPOS(I,O1,K2,N) &
                    + NCON_SLEAVE(Q,K2,N) * SOLA_XPOS(I,O1,K1,N)
             PXR1 =   PCON_SLEAVE(Q,K1,N) * SOLB_XNEG(I,O1,K1,N) &
                    - PCON_SLEAVE(Q,K2,N) * SOLB_XNEG(I,O1,K2,N)
             H1 =   NXR1 * T_DELT_EIGEN(K1,N) &
                  - NXR2 * T_DELT_EIGEN(K2,N)
             H2 =  PXR1
             SHOM_CR = SHOM_CR + H1 + H2
            ENDDO

!  collect solution

            QSURFACEWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER*(SHOM_R+SHOM_CR)

!  Finish loops

          ENDDO
        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSSL_QUADWF_LEVEL_DN

!

      SUBROUTINE LSSL_QUADWF_OFFGRID_UP ( &
        UTA, UT, N, N_SLEAVE_WFS, NSTOKES, NSTREAMS,  &
        FLUX_MULTIPLIER, T_UTUP_EIGEN, T_UTDN_EIGEN,  &
        K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,      &
        NCON_SLEAVE, PCON_SLEAVE, QSURFACEWF_F )

      USE VLIDORT_PARS

      IMPLICIT NONE

!  input

      INTEGER, INTENT (IN) ::          UTA, UT, N, N_SLEAVE_WFS
      INTEGER, INTENT (IN) ::          NSTOKES, NSTREAMS
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER

      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )

      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: NCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )

!  output

      DOUBLE PRECISION, INTENT (INOUT) :: QSURFACEWF_F &
          ( MAX_SLEAVEWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          I, I1, O1, K, KO1, K0, K1, K2, Q
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, H1, H2
      DOUBLE PRECISION :: NXR, NXR1, NXR2, PXR, PXR1, PXR2

!  For those optical depths at off-grid levels
!  -------------------------------------------

!  Offset

      KO1 = K_REAL(N) + 1

!  parameter loop, stream and Stokes loops

      DO Q = 1, N_SLEAVE_WFS
       DO I = 1, NSTREAMS
        I1 = I + NSTREAMS
        DO O1 = 1, NSTOKES

!  real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, K_REAL(N)
            NXR = NCON_SLEAVE(Q,K,N) * SOLA_XPOS(I1,O1,K,N)
            PXR = PCON_SLEAVE(Q,K,N) * SOLB_XNEG(I1,O1,K,N)
            SHOM_R = SHOM_R + NXR * T_UTDN_EIGEN(K,UT) &
                            + PXR * T_UTUP_EIGEN(K,UT)
          ENDDO

!  complex homogeneous solutions

          SHOM_CR = ZERO
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            NXR1 =   NCON_SLEAVE(Q,K1,N) * SOLA_XPOS(I1,O1,K1,N) &
                   - NCON_SLEAVE(Q,K2,N) * SOLA_XPOS(I1,O1,K2,N)
            NXR2 =   NCON_SLEAVE(Q,K1,N) * SOLA_XPOS(I1,O1,K2,N) &
                   + NCON_SLEAVE(Q,K2,N) * SOLA_XPOS(I1,O1,K1,N)
            PXR1 =   PCON_SLEAVE(Q,K1,N) * SOLB_XNEG(I1,O1,K1,N) &
                   - PCON_SLEAVE(Q,K2,N) * SOLB_XNEG(I1,O1,K2,N)
            PXR2 =   PCON_SLEAVE(Q,K1,N) * SOLB_XNEG(I1,O1,K2,N) &
                   + PCON_SLEAVE(Q,K2,N) * SOLB_XNEG(I1,O1,K1,N)
            H1 =   NXR1 * T_UTDN_EIGEN(K1,UT) &
                 - NXR2 * T_UTDN_EIGEN(K2,UT)
            H2 =   PXR1 * T_UTUP_EIGEN(K1,UT) &
                 - PXR2 * T_UTUP_EIGEN(K2,UT)
            SHOM_CR = SHOM_CR + H1 + H2
          ENDDO

!  set result

          QSURFACEWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER*(SHOM_R+SHOM_CR)

!  end O1 and I loops, end Q loop

        ENDDO
       ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSSL_QUADWF_OFFGRID_UP

!

      SUBROUTINE LSSL_QUADWF_OFFGRID_DN ( &
        UTA, UT, N, N_SLEAVE_WFS, NSTOKES, NSTREAMS,  &
        FLUX_MULTIPLIER, T_UTUP_EIGEN, T_UTDN_EIGEN,  &
        K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,      &
        NCON_SLEAVE, PCON_SLEAVE, QSURFACEWF_F )

      USE VLIDORT_PARS

      IMPLICIT NONE

!  input

      INTEGER, INTENT (IN) ::          UTA, UT, N, N_SLEAVE_WFS
      INTEGER, INTENT (IN) ::          NSTOKES, NSTREAMS
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER

      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )

      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: NCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )

!  output

      DOUBLE PRECISION, INTENT (INOUT) :: QSURFACEWF_F &
          ( MAX_SLEAVEWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          I, O1, K, KO1, K0, K1, K2, Q
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, H1, H2
      DOUBLE PRECISION :: NXR, NXR1, NXR2, PXR, PXR1, PXR2

!  For those optical depths at off-grid levels
!  -------------------------------------------

!  Offset

      KO1 = K_REAL(N) + 1

!  parameter loop, stream and Stokes loops

      DO Q = 1, N_SLEAVE_WFS
       DO I = 1, NSTREAMS
        DO O1 = 1, NSTOKES

!  real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, K_REAL(N)
            NXR = NCON_SLEAVE(Q,K,N) * SOLA_XPOS(I,O1,K,N)
            PXR = PCON_SLEAVE(Q,K,N) * SOLB_XNEG(I,O1,K,N)
            SHOM_R = SHOM_R + NXR * T_UTDN_EIGEN(K,UT) &
                            + PXR * T_UTUP_EIGEN(K,UT)
          ENDDO

!  complex homogeneous solutions

          SHOM_CR = ZERO
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            NXR1 =   NCON_SLEAVE(Q,K1,N) * SOLA_XPOS(I,O1,K1,N) &
                   - NCON_SLEAVE(Q,K2,N) * SOLA_XPOS(I,O1,K2,N)
            NXR2 =   NCON_SLEAVE(Q,K1,N) * SOLA_XPOS(I,O1,K2,N) &
                   + NCON_SLEAVE(Q,K2,N) * SOLA_XPOS(I,O1,K1,N)
            PXR1 =   PCON_SLEAVE(Q,K1,N) * SOLB_XNEG(I,O1,K1,N) &
                   - PCON_SLEAVE(Q,K2,N) * SOLB_XNEG(I,O1,K2,N)
            PXR2 =   PCON_SLEAVE(Q,K1,N) * SOLB_XNEG(I,O1,K2,N) &
                   + PCON_SLEAVE(Q,K2,N) * SOLB_XNEG(I,O1,K1,N)
            H1 =   NXR1 * T_UTDN_EIGEN(K1,UT) &
                 - NXR2 * T_UTDN_EIGEN(K2,UT)
            H2 =   PXR1 * T_UTUP_EIGEN(K1,UT) &
                 - PXR2 * T_UTUP_EIGEN(K2,UT)
            SHOM_CR = SHOM_CR + H1 + H2
          ENDDO

!  set result

          QSURFACEWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER*(SHOM_R+SHOM_CR)

!  end O1 and I loops, end Q loop

        ENDDO
       ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSSL_QUADWF_OFFGRID_DN

!  Finish module

      END MODULE vlidort_ls_wfsleave

