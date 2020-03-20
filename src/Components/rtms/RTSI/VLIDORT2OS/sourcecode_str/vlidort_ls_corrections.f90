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
! #            VLIDORT_LS_DBCORRECTION                          #
! #                                                             #
! ###############################################################


      MODULE vlidort_ls_corrections

      PRIVATE
      PUBLIC :: VLIDORT_LS_DBCORRECTION

      CONTAINS

      SUBROUTINE VLIDORT_LS_DBCORRECTION ( &
        FLUXMULT, &
        DO_SSCORR_OUTGOING, DO_UPWELLING, &
        DO_OBSERVATION_GEOMETRY, &
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
      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY
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
      INTEGER ::          UT, UTA, UM, UA, NC, IB, V, O1, O2, Q, OM, LUM, LUA
      DOUBLE PRECISION :: TR, SUM, L_KER, REFL_ATTN

!  Local user indices

      LUM = 1
      LUA = 1

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

      IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
        DO UM = 1, N_USER_STREAMS
          DO IB = 1, NBEAMS
            DO UA = 1, N_USER_RELAZMS
              V = VZA_OFFSETS(IB,UM) + UA
              IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
                REFL_ATTN = ATTN_DB_SAVE(V)
                DO Q = 1, NSURFWFS
                  DO O1 = 1, NELEMENTS
                    IF ( DO_LAMBERTIAN_SURFACE ) THEN
!  @@@ Rob fix        LS_EXACTDB_SOURCE(Q,V,O1) = REFL_ATTN * LAMBERTIAN_ALBEDO
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
      ELSE
        DO V = 1, N_GEOMETRIES
          IF ( DO_REFLECTED_DIRECTBEAM(V) ) THEN
            REFL_ATTN = ATTN_DB_SAVE(V)
            DO Q = 1, NSURFWFS
              DO O1 = 1, NELEMENTS
                IF ( DO_LAMBERTIAN_SURFACE ) THEN
!  @@@ Rob fix     LS_EXACTDB_SOURCE(Q,V,O1) = REFL_ATTN * LAMBERTIAN_ALBEDO
!      (Albedo WF must be Unnormalized)
                  LS_EXACTDB_SOURCE(Q,V,O1) = REFL_ATTN
                ELSE
                  SUM = ZERO
                  DO O2 = 1, NSTOKES
                    OM = MUELLER_INDEX(O1,O2)
                    L_KER = LS_EXACTDB_BRDFUNC(Q,O1,LUM,LUA,V)
                    SUM = SUM + L_KER * FLUXVEC(O2)
                  ENDDO
                  LS_EXACTDB_SOURCE(Q,V,O1) = REFL_ATTN * SUM
                ENDIF
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDIF

!  ====================================================
!  2. TRANSMITTANCE OF SOURCE TERM: UPWELLING RECURSION
!  ====================================================

      DO Q = 1, NSURFWFS

!  INITIALIZE CUMULATIVE SOURCE TERM

!mick fix  8/8/2013 - added if structure and else clause
!mick fix 11/8/2013 - modifed if structure
        NC =  0
        IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
          DO UM = 1, N_USER_STREAMS
            DO IB = 1, NBEAMS
              DO UA = 1, N_USER_RELAZMS
                V = VZA_OFFSETS(IB,UM) + UA
                IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
                  DO O1 = 1, NELEMENTS
                    LS_DB_CUMSOURCE(V,O1) = LS_EXACTDB_SOURCE(Q,V,O1)*FLUXMULT
                  ENDDO
                ELSE
                  DO O1 = 1, NELEMENTS
                    LS_DB_CUMSOURCE(V,O1) = ZERO
                  ENDDO
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO V = 1, N_GEOMETRIES
            IF ( DO_REFLECTED_DIRECTBEAM(V) ) THEN
              DO O1 = 1, NELEMENTS
                LS_DB_CUMSOURCE(V,O1) = LS_EXACTDB_SOURCE(Q,V,O1)*FLUXMULT
              ENDDO
            ELSE
              DO O1 = 1, NELEMENTS
                LS_DB_CUMSOURCE(V,O1) = ZERO
              ENDDO
            ENDIF
          ENDDO
        ENDIF

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
                  LS_DB_CUMSOURCE(V,O1) = TR * LS_DB_CUMSOURCE(V,O1)
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
                LS_DB_CUMSOURCE(V,O1) = TR * LS_DB_CUMSOURCE(V,O1)
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

           IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
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
           ELSE
             DO V = 1, N_GEOMETRIES
              IF ( DO_SSCORR_OUTGOING ) THEN
                TR  = UP_LOSTRANS_UT(UT,V)
              ELSE
                TR  = T_UTUP_USERM(UT,V)
              ENDIF
              DO O1 = 1, NELEMENTS
               SURFACEWF_DB(Q,UTA,V,O1) = TR * LS_DB_CUMSOURCE(V,O1)
              ENDDO
             ENDDO
           ENDIF

!  ONGRID OUTPUT : SET FINAL CUMULATIVE SOURCE DIRECTLY

          ELSE
           IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
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
           ELSE
             DO V = 1, N_GEOMETRIES
              DO O1 = 1, NELEMENTS
               SURFACEWF_DB(Q,UTA,V,O1) = LS_DB_CUMSOURCE(V,O1)
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
      END SUBROUTINE VLIDORT_LS_DBCORRECTION

      END MODULE vlidort_ls_corrections

