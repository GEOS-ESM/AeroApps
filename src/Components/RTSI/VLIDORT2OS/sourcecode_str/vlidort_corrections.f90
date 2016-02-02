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
! #                                                             #
! #      Version 2.0. Nadir single scatter correction           #
! #                                                             #
! #            VLIDORT_SSCORR_NADIR (master)                    #
! #              VLIDORTSS_FMATRICES                            #
! #              VLIDORTSS_FMATRICES_MULTI                      #
! #                                                             #
! #      Version 2.2. outgoing sphericity correction            #
! #              2.3. partial-layer integration                 #
! #                                                             #
! #            VLIDORT_SSCORR_OUTGOING (master)                 #
! #                 SSCORR_OUTGOING_ZMATRIX                     #
! #                 OUTGOING_INTEGRATION_UP                     #
! #                 OUTGOING_INTEGRATION_DN                     #
! #                                                             #
! #            VLIDORT_DBCORRECTION                             #
! #                                                             #
! #            VFO_MASTER_INTERFACE                             #
! #                                                             #
! #      Version 2.4RTC ------- Notes -------                   #
! #                                                             #
! #            Additional correction to ZMATRIX and FMATRIX     #
! #            routines for Azimuth > 180. HvdM, Eqs(94/95)     #
! #            Implemented by V. Natraj and R. Spurr, 5/1/09    #
! #            Correctly done, October 2010.                    #
! #                                                             #
! #            Correction for Downwelling SSCORR_OUTGOING       #
! #            Separate Geometry routine implemented.           #
! #            Implemented by R. Spurr, 10/06/10                #
! #                                                             #
! #            Consolidated DB correction (Lambertian/BRDF)     #
! #                                                             #
! ###############################################################


      MODULE vlidort_corrections

      USE FO_VectorSS_masters

      PRIVATE
      PUBLIC :: VLIDORT_SSCORR_NADIR, &
                VLIDORT_SSCORR_OUTGOING, &
                VLIDORT_DBCORRECTION, &
                VFO_MASTER_INTERFACE

      CONTAINS

      SUBROUTINE VLIDORT_SSCORR_NADIR ( &
        SSFLUX, &
        DO_SSCORR_TRUNCATION, DO_REFRACTIVE_GEOMETRY, &
        DO_DELTAM_SCALING, DO_UPWELLING, &
        DO_DNWELLING, DO_OBSERVATION_GEOMETRY, NSTOKES, &
        NLAYERS, NGREEK_MOMENTS_INPUT, &
        N_USER_RELAZMS, USER_RELAZMS, &
        N_USER_LEVELS, &
        OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT, &
        DO_TOA_CONTRIBS, &
        FLUXVEC, COS_SZANGLES, &
        SIN_SZANGLES, SUN_SZA_COSINES, &
        NBEAMS, N_USER_STREAMS, &
        LAYER_MAXMOMENTS, USER_STREAMS, &
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
        UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, &
        PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, &
        STERM_LAYERMASK_DN, N_GEOMETRIES, &
        VZA_OFFSETS, GREEKMAT_TOTAL, &
        TRUNC_FACTOR, T_DELT_USERM, &
        T_UTDN_USERM, T_UTUP_USERM, &
        CUMTRANS, &
        EMULT_UP, EMULT_DN, &
        UT_EMULT_UP, UT_EMULT_DN, &
        ZMAT_UP, ZMAT_DN, TMS, &
        SSFDEL, SS_CUMSOURCE_UP, SS_CUMSOURCE_DN, &
        STOKES_SS, SS_CONTRIBS )

!  Single scatter exact calculation. Nadir viewing output
!   Programmed by R. Spurr, RT Solutions Inc.
!    Second Draft, April 14th 2005.
!    Third Draft,  May    6th 2005.
!       - only generates SS Stokes vector.
!       - additional code to deal with refraction.

      USE VLIDORT_PARS

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT (IN) ::  SSFLUX

      LOGICAL, INTENT (IN) ::           DO_SSCORR_TRUNCATION
      LOGICAL, INTENT (IN) ::           DO_REFRACTIVE_GEOMETRY
      LOGICAL, INTENT (IN) ::           DO_DELTAM_SCALING
      LOGICAL, INTENT (IN) ::           DO_UPWELLING
      LOGICAL, INTENT (IN) ::           DO_DNWELLING
      LOGICAL, INTENT (IN) ::           DO_OBSERVATION_GEOMETRY
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           NGREEK_MOMENTS_INPUT
      INTEGER, INTENT (IN) ::           N_USER_RELAZMS
      DOUBLE PRECISION, INTENT (IN) ::  USER_RELAZMS  ( MAX_USER_RELAZMS )
      INTEGER, INTENT (IN) ::           N_USER_LEVELS
      DOUBLE PRECISION, INTENT (IN) ::  OMEGA_TOTAL_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  GREEKMAT_TOTAL_INPUT &
          ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
      LOGICAL, INTENT (IN) ::           DO_TOA_CONTRIBS
      DOUBLE PRECISION, INTENT (IN) ::  FLUXVEC ( MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  COS_SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) ::  SIN_SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) ::  SUN_SZA_COSINES &
          ( MAXLAYERS, MAX_SZANGLES )
      INTEGER, INTENT (IN) ::           NBEAMS
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      INTEGER, INTENT (IN) ::           LAYER_MAXMOMENTS ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  USER_STREAMS  ( MAX_USER_STREAMS )
      LOGICAL, INTENT (IN) ::           PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      LOGICAL, INTENT (IN) ::           STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::           STERM_LAYERMASK_DN ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           N_GEOMETRIES
      INTEGER, INTENT (IN) ::           VZA_OFFSETS &
          ( MAX_SZANGLES, MAX_USER_VZANGLES )
      DOUBLE PRECISION, INTENT (IN) ::  GREEKMAT_TOTAL &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION, INTENT (IN) ::  TRUNC_FACTOR ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_UTDN_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_UTUP_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  CUMTRANS ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  EMULT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  EMULT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_EMULT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_EMULT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )

      DOUBLE PRECISION, INTENT (OUT) :: ZMAT_UP &
          ( MAX_GEOMETRIES, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: ZMAT_DN &
          ( MAX_GEOMETRIES, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: TMS ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: SSFDEL ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: SS_CUMSOURCE_UP &
          ( MAX_GEOMETRIES, MAXSTOKES, 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: SS_CUMSOURCE_DN &
          ( MAX_GEOMETRIES, MAXSTOKES, 0:MAXLAYERS )

      DOUBLE PRECISION, INTENT (INOUT) :: STOKES_SS &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )

      DOUBLE PRECISION, INTENT (OUT) :: SS_CONTRIBS &
          ( MAX_GEOMETRIES, MAXSTOKES, MAXLAYERS )

!  local variables
!  ---------------

!  Indices

      INTEGER ::           N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER ::           UT, UTA, UM, IA, NC, IB, V, O1, O2, NM1, LUM
      INTEGER ::           INDEX_11, INDEX_12, INDEX_34
      INTEGER ::           INDEX_22, INDEX_33, INDEX_44

!  help variables (double precision)

      DOUBLE PRECISION ::  VIEWSIGN, MUX, DNM1
      DOUBLE PRECISION ::  FINAL_SOURCE, HELP, SS_LAYERSOURCE
      DOUBLE PRECISION ::  MULTIPLIER, SSCORRECTION
      DOUBLE PRECISION ::  TR_CUMSOURCE, SS_CUMSOURCE
      DOUBLE PRECISION ::  HELP2C1, HELP2S1, HELP3C1, HELP3S1

      DOUBLE PRECISION :: FMAT_UP ( MAX_GEOMETRIES, MAXLAYERS, 6 )
      DOUBLE PRECISION :: FMAT_DN ( MAX_GEOMETRIES, MAXLAYERS, 6 )

!  zenith angle cosines/sines, azimuth angle cosines

      DOUBLE PRECISION ::  CTHETA (MAXLAYERS,MAXBEAMS)
      DOUBLE PRECISION ::  STHETA (MAXLAYERS,MAXBEAMS)

      DOUBLE PRECISION ::  CALPHA (MAX_USER_STREAMS)
      DOUBLE PRECISION ::  SALPHA (MAX_USER_STREAMS)
      DOUBLE PRECISION ::  CPHI   (MAX_USER_RELAZMS)

!  output Rotation cosines/sines from F-Matrix subroutine

      DOUBLE PRECISION ::  C1 (MAX_GEOMETRIES,MAXLAYERS)
      DOUBLE PRECISION ::  S1 (MAX_GEOMETRIES,MAXLAYERS)
      DOUBLE PRECISION ::  C2 (MAX_GEOMETRIES,MAXLAYERS)
      DOUBLE PRECISION ::  S2 (MAX_GEOMETRIES,MAXLAYERS)

!  Sunlight parameters
!   - assumes the Flux vector is (F,0,0,0)
!   - drop this variable when further testing has been done

      LOGICAL, PARAMETER :: SUNLIGHT = .TRUE.
      INTEGER ::            NPOLAR

!  Set up operations
!  -----------------

!  Set npolar. For Natural light, there are only 3 Stokes parameters
!    ( no circular polarization for single scattering)

      IF (SUNLIGHT) NPOLAR = MIN(NSTOKES,3)
      IF (.NOT.SUNLIGHT) NPOLAR = NSTOKES

!  Local user index

      LUM = 1

!  indexing keys

      INDEX_11 = 1
      INDEX_12 = 2
      INDEX_22 = 3
      INDEX_33 = 4
      INDEX_34 = 5
      INDEX_44 = 6

!  Create TMS factors, these get stored
!    Delta-M Scaling introduced April 2005.

      DO N = 1, NLAYERS
        IF ( DO_DELTAM_SCALING ) THEN
          HELP   = ONE - TRUNC_FACTOR(N) * OMEGA_TOTAL_INPUT(N)
          TMS(N) = OMEGA_TOTAL_INPUT(N) / HELP
        ELSE
          TMS(N) = OMEGA_TOTAL_INPUT(N)
        ENDIF
      ENDDO

!  Additional Delta-M scaling
!  New section. R. Spurr, 07 September 2007.
!   TMS gets modified by (1-F). Save the truncation factor.
!   Phase function moments are modified later on.

      IF ( DO_SSCORR_TRUNCATION ) THEN
        NM1  = NGREEK_MOMENTS_INPUT
        DNM1 = DBLE(2*NM1+1)
        DO N = 1, NLAYERS
          SSFDEL(N) = GREEKMAT_TOTAL_INPUT(NM1,N,1) / DNM1
          TMS(N) = TMS(N) * ( ONE - SSFDEL(N) )
        ENDDO
      ENDIF

!  save some geometrical quantities (cosines and sines)
!    Multiple layer quantities for the Refractive case

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

!  Upwelling single scatter Phase matrices
!  ---------------------------------------

      IF ( DO_UPWELLING ) THEN

!  Get Rotation cosine/sines and F matrices for all layers
!    Distinguish between refractive case and straightline case

        VIEWSIGN = -1.0D0
        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
          CALL VLIDORTSS_FMATRICES_MULTI ( &
            DO_SSCORR_TRUNCATION, DO_OBSERVATION_GEOMETRY, NLAYERS, &
            NSTOKES, N_GEOMETRIES, NGREEK_MOMENTS_INPUT, &
            NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, &
            LAYER_MAXMOMENTS, GREEKMAT_TOTAL_INPUT, SSFDEL, &
            VIEWSIGN, VZA_OFFSETS, &
            CTHETA, STHETA, CALPHA, SALPHA, CPHI, USER_RELAZMS, &
            C1, S1, C2, S2, FMAT_UP )
        ELSE
          CALL VLIDORTSS_FMATRICES ( &
            DO_SSCORR_TRUNCATION, DO_OBSERVATION_GEOMETRY, NLAYERS, &
            NSTOKES, N_GEOMETRIES, NGREEK_MOMENTS_INPUT, &
            NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, &
            LAYER_MAXMOMENTS, GREEKMAT_TOTAL_INPUT, SSFDEL, &
            VIEWSIGN, VZA_OFFSETS, &
            CTHETA, STHETA, CALPHA, SALPHA, CPHI, USER_RELAZMS, &
            C1, S1, C2, S2, FMAT_UP )
        ENDIF

!   Z matrices, Hovenier and van der Mee (1983), Eq 88.
!   ---------------------------------------------------

!mick fix 1/16/2013 - initialize all entries
        ZMAT_UP = ZERO

        IF ( NSTOKES .GT. 1 ) THEN

!  For sunlight, only need the first Column of the Z-matrix

         IF ( SUNLIGHT ) THEN
          DO V = 1, N_GEOMETRIES
           DO N = 1, NLAYERS
            IF ( STERM_LAYERMASK_UP(N)) THEN
             ZMAT_UP(V,N,1,1) = FMAT_UP(V,N,INDEX_11)
             ZMAT_UP(V,N,2,1) =  -FMAT_UP(V,N,INDEX_12) * C2(V,N)
             ZMAT_UP(V,N,3,1) =   FMAT_UP(V,N,INDEX_12) * S2(V,N)
             ZMAT_UP(V,N,4,1) =   ZERO
            ENDIF
           ENDDO
          ENDDO
         ENDIF

!  For general case. Not tested as of 15 April 2005.

         IF ( .NOT. SUNLIGHT ) THEN
          DO V = 1, N_GEOMETRIES
           DO N = 1, NLAYERS
            IF ( STERM_LAYERMASK_UP(N)) THEN
             ZMAT_UP(V,N,1,1) = FMAT_UP(V,N,INDEX_11)
             ZMAT_UP(V,N,2,1) =  -FMAT_UP(V,N,INDEX_12) * C2(V,N)
             ZMAT_UP(V,N,3,1) =   FMAT_UP(V,N,INDEX_12) * S2(V,N)
             ZMAT_UP(V,N,4,1) =   ZERO
!  this next part remains untested. Need to check signs.
             HELP2C1 = FMAT_UP(V,N,INDEX_22) * C1(V,N)
             HELP2S1 = FMAT_UP(V,N,INDEX_22) * S1(V,N)
             HELP3C1 = FMAT_UP(V,N,INDEX_33) * C1(V,N)
             HELP3S1 = FMAT_UP(V,N,INDEX_33) * S1(V,N)
             ZMAT_UP(V,N,1,2) =   FMAT_UP(V,N,INDEX_12) * C1(V,N)
             ZMAT_UP(V,N,1,3) = - FMAT_UP(V,N,INDEX_12) * S1(V,N)
             ZMAT_UP(V,N,1,4) =   ZERO
             ZMAT_UP(V,N,2,2) =   C2(V,N) * HELP2C1 - S2(V,N) * HELP3S1
             ZMAT_UP(V,N,2,3) = - C2(V,N) * HELP2S1 - S2(V,N) * HELP3C1
             ZMAT_UP(V,N,2,4) = - FMAT_UP(V,N,INDEX_34) * S2(V,N)
             ZMAT_UP(V,N,3,2) =   S2(V,N) * HELP2C1 + C2(V,N) * HELP3S1
             ZMAT_UP(V,N,3,3) = - S2(V,N) * HELP2S1 + C2(V,N) * HELP3C1
             ZMAT_UP(V,N,3,4) =   FMAT_UP(V,N,INDEX_34) * C2(V,N)
             ZMAT_UP(V,N,4,2) = - FMAT_UP(V,N,INDEX_34) * S1(V,N)
             ZMAT_UP(V,N,4,3) = - FMAT_UP(V,N,INDEX_34) * C1(V,N)
             ZMAT_UP(V,N,4,4) =   FMAT_UP(V,N,INDEX_44)
            ENDIF
           ENDDO
          ENDDO
         ENDIF

!  end NSTOKES > 1

        ENDIF

!  debug
!         write(97,'(2i4,1pe15.7)')V,21,ZMAT_UP(V,21,1,1)
!         write(97,'(2i4,1pe15.7)')V,22,ZMAT_UP(V,22,1,1)
!        pause

!  add TMS correction factor, scalar case

        IF ( NSTOKES .EQ. 1 ) THEN
         DO N = 1, NLAYERS
          IF ( STERM_LAYERMASK_UP(N)) THEN
           DO V = 1, N_GEOMETRIES
            ZMAT_UP(V,N,1,1) = FMAT_UP(V,N,INDEX_11) * TMS(N)
!            write(34,'(2i5,1pe22.12)')V,N,ZMAT_UP(V,N,1,1)
           ENDDO
          ENDIF
         ENDDO
        ENDIF

!  add TMS correction factor, vector case

        IF ( NSTOKES .GT. 1 ) THEN

!  sunlight only

         IF ( SUNLIGHT ) THEN
          DO N = 1, NLAYERS
           IF ( STERM_LAYERMASK_UP(N)) THEN
            DO V = 1, N_GEOMETRIES
             DO O1 = 1, NSTOKES
!        if ( n.gt.nlayers-2)write(*,*)'Reg',N,V,O1,ZMAT_UP(V,N,O1,1)
              ZMAT_UP(V,N,O1,1) = ZMAT_UP(V,N,O1,1) * TMS(N)
             ENDDO
            ENDDO
           ENDIF
          ENDDO

!  general case

         ELSE
          DO N = 1, NLAYERS
           IF ( STERM_LAYERMASK_UP(N)) THEN
            DO V = 1, N_GEOMETRIES
             DO O1 = 1, NSTOKES
              DO O2 = 1, NSTOKES
               ZMAT_UP(V,N,O1,O2) = ZMAT_UP(V,N,O1,O2) * TMS(N)
              ENDDO
             ENDDO
            ENDDO
           ENDIF
          ENDDO
         ENDIF

!  finish this section

        ENDIF
      ENDIF

!  Upwelling single scatter recurrence
!  -----------------------------------

      IF ( DO_UPWELLING ) THEN

!  initialize cumulative source term

        NC =  0
        DO V = 1, N_GEOMETRIES
          DO O1 = 1, NSTOKES
            SS_CUMSOURCE_UP(V,O1,NC) = ZERO
          ENDDO
        ENDDO

!  Lattice calculation

        IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN

!  initialise optical depth loop

          NSTART = NLAYERS
          NUT_PREV = NSTART + 1

!  Main loop over all output optical depths

          DO UTA = N_USER_LEVELS, 1, -1

!  Layer index for given optical depth

            NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
            NUT    = NLEVEL + 1

!  Cumulative single scatter source terms :
!      For loop over layers working upwards to level NUT,
!      Get layer source terms = Exact Z-matrix * Multiplier

            DO N = NSTART, NUT, -1
             NC = NLAYERS + 1 - N

!  sunlight case, no circular polarization
!    Addition of Single scatter Contribution functions

             IF ( SUNLIGHT ) THEN
              DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                MULTIPLIER = EMULT_UP(UM,N,IB)
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO O1 = 1, NPOLAR
                  HELP = ZMAT_UP(V,N,O1,1) * FLUXVEC(1)
                  SS_LAYERSOURCE = HELP* MULTIPLIER
                  IF (DO_TOA_CONTRIBS) THEN
                    SS_CONTRIBS(V,O1,N) = CUMTRANS(N,UM) * SS_LAYERSOURCE
                    SS_CONTRIBS(V,O1,N) = SSFLUX * SS_CONTRIBS(V,O1,N)
                  ENDIF
                  SS_CUMSOURCE_UP(V,O1,NC) = SS_LAYERSOURCE + &
                     T_DELT_USERM(N,UM)*SS_CUMSOURCE_UP(V,O1,NC-1)
                 ENDDO
                ENDDO
               ENDDO
              ENDDO

!  general case

             ELSE
              DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                MULTIPLIER = EMULT_UP(UM,N,IB)
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO O1 = 1, NSTOKES
                  HELP = ZERO
                  DO O2 = 1, NSTOKES
                   HELP = HELP + ZMAT_UP(V,N,O1,O2) * FLUXVEC(O2)
                  ENDDO
                  SS_LAYERSOURCE = HELP* MULTIPLIER
                  IF (DO_TOA_CONTRIBS) THEN
                    SS_CONTRIBS(V,O1,N) = CUMTRANS(N,UM) * SS_LAYERSOURCE
                    SS_CONTRIBS(V,O1,N) = SSFLUX * SS_CONTRIBS(V,O1,N)
                  ENDIF
                  SS_CUMSOURCE_UP(V,O1,NC) = SS_LAYERSOURCE + &
                     T_DELT_USERM(N,UM)*SS_CUMSOURCE_UP(V,O1,NC-1)
                 ENDDO
                ENDDO
               ENDDO
              ENDDO
             ENDIF

!  end layer loop

            ENDDO

!  Offgrid output :
!    add additional partial layer source term = Exact Z-matrix * Multipl
!    Set final cumulative source and Correct the intensity

            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
              UT = PARTLAYERS_OUTINDEX(UTA)
              N  = PARTLAYERS_LAYERIDX(UT)

!  sunlight case, no circular polarization

              IF ( SUNLIGHT ) THEN

               DO IB = 1, NBEAMS
                DO UM = 1, N_USER_STREAMS
                 MULTIPLIER = UT_EMULT_UP(UM,UT,IB)
                 DO IA = 1, N_USER_RELAZMS
                  V = VZA_OFFSETS(IB,UM) + IA
                  DO O1 = 1, NPOLAR
                   HELP = ZMAT_UP(V,N,O1,1) * FLUXVEC(1)
                   SS_LAYERSOURCE = HELP * MULTIPLIER
                   SS_CUMSOURCE   = SS_CUMSOURCE_UP(V,O1,NC)
                   TR_CUMSOURCE   = T_UTUP_USERM(UT,UM) * SS_CUMSOURCE
                   FINAL_SOURCE   = TR_CUMSOURCE + SS_LAYERSOURCE
                   SSCORRECTION   = SSFLUX * FINAL_SOURCE
                   STOKES_SS(UTA,V,O1,UPIDX) = SSCORRECTION
                  ENDDO
                 ENDDO
                ENDDO
               ENDDO

!  general case

              ELSE

               DO IB = 1, NBEAMS
                DO UM = 1, N_USER_STREAMS
                 MULTIPLIER = UT_EMULT_UP(UM,UT,IB)
                 DO IA = 1, N_USER_RELAZMS
                  V = VZA_OFFSETS(IB,UM) + IA
                  DO O1 = 1, NSTOKES
                   HELP = ZERO
                   DO O2 = 1, NSTOKES
                    HELP = HELP + ZMAT_UP(V,N,O1,O2) * FLUXVEC(O2)
                   ENDDO
                   SS_LAYERSOURCE = HELP * MULTIPLIER
                   SS_CUMSOURCE   = SS_CUMSOURCE_UP(V,O1,NC)
                   TR_CUMSOURCE   = T_UTUP_USERM(UT,UM) * SS_CUMSOURCE
                   FINAL_SOURCE   = TR_CUMSOURCE + SS_LAYERSOURCE
                   SSCORRECTION   = SSFLUX * FINAL_SOURCE
                   STOKES_SS(UTA,V,O1,UPIDX) = SSCORRECTION
                  ENDDO
                 ENDDO
                ENDDO
               ENDDO

              ENDIF

!  debug code to be inserted above
!                     write(98,'(a,3i4,1p2e15.7)')
!     &                 'up  ',UTA,IB,V,UNCORRECTED,SSCORRECTION

!  Ongrid output :
!     Set final cumulative source and correct Stokes vector

            ELSE
              DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO O1 = 1, NPOLAR
                  FINAL_SOURCE = SS_CUMSOURCE_UP(V,O1,NC)
                  SSCORRECTION = SSFLUX * FINAL_SOURCE
                  STOKES_SS(UTA,V,O1,UPIDX) = SSCORRECTION
!              write(34,'(2i5,1pe22.12)')V,UTA,STOKES_SS(UTA,V,O1,UPIDX)
                 ENDDO
                ENDDO
               ENDDO
              ENDDO
            ENDIF

!  insert this debug code in the above loops
!                     write(97,'(a,3i4,1p2e15.7)')
!     &                 'up  ',UTA,IB,V,UNCORRECTED,SSCORRECTION

!  Check for updating the recursion

            IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
            NUT_PREV = NUT

!  end optical depth loop

          ENDDO

!  End lattice calculation

        ENDIF

!  Observational geometry calculation

        IF ( DO_OBSERVATION_GEOMETRY ) THEN

!  initialise optical depth loop

          NSTART = NLAYERS
          NUT_PREV = NSTART + 1

!  Main loop over all output optical depths

          DO UTA = N_USER_LEVELS, 1, -1

!  Layer index for given optical depth

            NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
            NUT    = NLEVEL + 1

!  Cumulative single scatter source terms :
!      For loop over layers working upwards to level NUT,
!      Get layer source terms = Exact Z-matrix * Multiplier

            DO N = NSTART, NUT, -1
             NC = NLAYERS + 1 - N

!  sunlight case, no circular polarization
!    Addition of Single scatter Contribution functions

             IF ( SUNLIGHT ) THEN
               DO V = 1, N_GEOMETRIES
                 MULTIPLIER = EMULT_UP(LUM,N,V)
                 DO O1 = 1, NPOLAR
                   HELP = ZMAT_UP(V,N,O1,1) * FLUXVEC(1)
                   SS_LAYERSOURCE = HELP* MULTIPLIER
                   IF (DO_TOA_CONTRIBS) THEN
                     SS_CONTRIBS(V,O1,N) = CUMTRANS(N,V) * SS_LAYERSOURCE
                     SS_CONTRIBS(V,O1,N) = SSFLUX * SS_CONTRIBS(V,O1,N)
                   ENDIF
                   SS_CUMSOURCE_UP(V,O1,NC) = SS_LAYERSOURCE + &
                      T_DELT_USERM(N,V)*SS_CUMSOURCE_UP(V,O1,NC-1)
                 ENDDO
               ENDDO

!  general case

             ELSE
               DO V = 1, N_GEOMETRIES
                 MULTIPLIER = EMULT_UP(LUM,N,V)
                 DO O1 = 1, NSTOKES
                   HELP = ZERO
                   DO O2 = 1, NSTOKES
                     HELP = HELP + ZMAT_UP(V,N,O1,O2) * FLUXVEC(O2)
                   ENDDO
                   SS_LAYERSOURCE = HELP* MULTIPLIER
                   IF (DO_TOA_CONTRIBS) THEN
                     SS_CONTRIBS(V,O1,N) = CUMTRANS(N,V) * SS_LAYERSOURCE
                     SS_CONTRIBS(V,O1,N) = SSFLUX * SS_CONTRIBS(V,O1,N)
                   ENDIF
                   SS_CUMSOURCE_UP(V,O1,NC) = SS_LAYERSOURCE + &
                      T_DELT_USERM(N,V)*SS_CUMSOURCE_UP(V,O1,NC-1)
                 ENDDO
               ENDDO
             ENDIF

!  end layer loop

            ENDDO

!  Offgrid output :
!    add additional partial layer source term = Exact Z-matrix * Multipl
!    Set final cumulative source and Correct the intensity

            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
              UT = PARTLAYERS_OUTINDEX(UTA)
              N  = PARTLAYERS_LAYERIDX(UT)

!  sunlight case, no circular polarization

              IF ( SUNLIGHT ) THEN

               DO V = 1, N_GEOMETRIES
                 MULTIPLIER = UT_EMULT_UP(LUM,UT,V)
                 DO O1 = 1, NPOLAR
                   HELP = ZMAT_UP(V,N,O1,1) * FLUXVEC(1)
                   SS_LAYERSOURCE = HELP * MULTIPLIER
                   SS_CUMSOURCE   = SS_CUMSOURCE_UP(V,O1,NC)
                   TR_CUMSOURCE   = T_UTUP_USERM(UT,V) * SS_CUMSOURCE
                   FINAL_SOURCE   = TR_CUMSOURCE + SS_LAYERSOURCE
                   SSCORRECTION   = SSFLUX * FINAL_SOURCE
                   STOKES_SS(UTA,V,O1,UPIDX) = SSCORRECTION
                 ENDDO
               ENDDO

!  general case

              ELSE

               DO V = 1, N_GEOMETRIES
                 MULTIPLIER = UT_EMULT_UP(LUM,UT,V)
                 DO O1 = 1, NSTOKES
                   HELP = ZERO
                   DO O2 = 1, NSTOKES
                    HELP = HELP + ZMAT_UP(V,N,O1,O2) * FLUXVEC(O2)
                   ENDDO
                   SS_LAYERSOURCE = HELP * MULTIPLIER
                   SS_CUMSOURCE   = SS_CUMSOURCE_UP(V,O1,NC)
                   TR_CUMSOURCE   = T_UTUP_USERM(UT,V) * SS_CUMSOURCE
                   FINAL_SOURCE   = TR_CUMSOURCE + SS_LAYERSOURCE
                   SSCORRECTION   = SSFLUX * FINAL_SOURCE
                   STOKES_SS(UTA,V,O1,UPIDX) = SSCORRECTION
                 ENDDO
               ENDDO

              ENDIF

!  debug code to be inserted above
!                     write(98,'(a,3i4,1p2e15.7)')
!     &                 'up  ',UTA,IB,V,UNCORRECTED,SSCORRECTION

!  Ongrid output :
!     Set final cumulative source and correct Stokes vector

            ELSE
              DO V = 1, N_GEOMETRIES
               DO O1 = 1, NPOLAR
                FINAL_SOURCE = SS_CUMSOURCE_UP(V,O1,NC)
                SSCORRECTION = SSFLUX * FINAL_SOURCE
                STOKES_SS(UTA,V,O1,UPIDX) = SSCORRECTION
               ENDDO
              ENDDO

            ENDIF

!  insert this debug code in the above loops
!                     write(97,'(a,3i4,1p2e15.7)')
!     &                 'up  ',UTA,IB,V,UNCORRECTED,SSCORRECTION

!  Check for updating the recursion

          IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT

!  end optical depth loop

        ENDDO

!  End observationl geometry calculation

        ENDIF

!  end Upwelling clause

      ENDIF

!  Downwelling single scatter Phase matrices
!  -----------------------------------------

      IF ( DO_DNWELLING ) THEN

!  Get Rotation cosine/sines and F matrices for all layers
!    Distinguish between refractive case and straightline case

        VIEWSIGN = +1.0D0
        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
          CALL VLIDORTSS_FMATRICES_MULTI &
        ( DO_SSCORR_TRUNCATION, DO_OBSERVATION_GEOMETRY, NLAYERS, &
          NSTOKES, N_GEOMETRIES, NGREEK_MOMENTS_INPUT, &
          NBEAMS, N_USER_STREAMS,  N_USER_RELAZMS, &
          LAYER_MAXMOMENTS,  GREEKMAT_TOTAL_INPUT, SSFDEL, &
          VIEWSIGN, VZA_OFFSETS, &
          CTHETA, STHETA, CALPHA, SALPHA, CPHI, USER_RELAZMS, &
          C1, S1, C2, S2, FMAT_DN )
        ELSE
          CALL VLIDORTSS_FMATRICES &
        ( DO_SSCORR_TRUNCATION, DO_OBSERVATION_GEOMETRY, NLAYERS, &
          NSTOKES, N_GEOMETRIES, NGREEK_MOMENTS_INPUT, &
          NBEAMS, N_USER_STREAMS,  N_USER_RELAZMS, &
          LAYER_MAXMOMENTS,  GREEKMAT_TOTAL_INPUT, SSFDEL, &
          VIEWSIGN, VZA_OFFSETS, &
          CTHETA, STHETA, CALPHA, SALPHA, CPHI, USER_RELAZMS, &
          C1, S1, C2, S2, FMAT_DN )
        ENDIF

!   Z matrices, Hovenier and van der Mee (1983), Eq 88.

!mick fix 1/16/2013 - initialize all entries
        ZMAT_DN = ZERO

        IF ( NSTOKES .GT. 1 ) THEN

!  For sunlight, only need the first column of Z matrix

         IF ( SUNLIGHT ) THEN
          DO V = 1, N_GEOMETRIES
           DO N = 1, NLAYERS
            IF ( STERM_LAYERMASK_DN(N)) THEN
             ZMAT_DN(V,N,1,1) = FMAT_DN(V,N,INDEX_11)
             ZMAT_DN(V,N,2,1) =   -FMAT_DN(V,N,INDEX_12) * C2(V,N)
             ZMAT_DN(V,N,3,1) =   FMAT_DN(V,N,INDEX_12) * S2(V,N)
             ZMAT_DN(V,N,4,1) =   ZERO
            ENDIF
           ENDDO
          ENDDO
         ENDIF

!  For general case. Not tested as of 15 April 2005.

         IF ( .NOT. SUNLIGHT ) THEN

          DO V = 1, N_GEOMETRIES
           DO N = 1, NLAYERS
            IF ( STERM_LAYERMASK_DN(N)) THEN
             ZMAT_DN(V,N,1,1) = FMAT_DN(V,N,INDEX_11)
             ZMAT_DN(V,N,2,1) =   -FMAT_DN(V,N,INDEX_12) * C2(V,N)
             ZMAT_DN(V,N,3,1) =   FMAT_DN(V,N,INDEX_12) * S2(V,N)
             ZMAT_DN(V,N,4,1) =   ZERO
!  This code untested
             HELP2C1 = FMAT_DN(V,N,INDEX_22) * C1(V,N)
             HELP2S1 = FMAT_DN(V,N,INDEX_22) * S1(V,N)
             HELP3C1 = FMAT_DN(V,N,INDEX_33) * C1(V,N)
             HELP3S1 = FMAT_DN(V,N,INDEX_33) * S1(V,N)
             ZMAT_DN(V,N,1,2) =   FMAT_DN(V,N,INDEX_12) * C1(V,N)
             ZMAT_DN(V,N,1,3) = - FMAT_DN(V,N,INDEX_12) * S1(V,N)
             ZMAT_DN(V,N,1,4) =   ZERO
             ZMAT_DN(V,N,2,2) =   C2(V,N) * HELP2C1 - S2(V,N) * HELP3S1
             ZMAT_DN(V,N,2,3) = - C2(V,N) * HELP2S1 - S2(V,N) * HELP3C1
             ZMAT_DN(V,N,2,4) = - FMAT_DN(V,N,INDEX_34) * S2(V,N)
             ZMAT_DN(V,N,3,2) =   S2(V,N) * HELP2C1 + C2(V,N) * HELP3S1
             ZMAT_DN(V,N,3,3) = - S2(V,N) * HELP2S1 + C2(V,N) * HELP3C1
             ZMAT_DN(V,N,3,4) =   FMAT_DN(V,N,INDEX_34) * C2(V,N)
             ZMAT_DN(V,N,4,2) = - FMAT_DN(V,N,INDEX_34) * S1(V,N)
             ZMAT_DN(V,N,4,3) = - FMAT_DN(V,N,INDEX_34) * C1(V,N)
             ZMAT_DN(V,N,4,4) =   FMAT_DN(V,N,INDEX_44)
            ENDIF
           ENDDO
          ENDDO
         ENDIF
        ENDIF

!  add TMS correction factor, scalar case

        IF ( NSTOKES .EQ. 1 ) THEN
         DO N = 1, NLAYERS
          IF ( STERM_LAYERMASK_DN(N)) THEN
           DO V = 1, N_GEOMETRIES
            ZMAT_DN(V,N,1,1) = FMAT_DN(V,N,INDEX_11) * TMS(N)
           ENDDO
          ENDIF
         ENDDO
        ENDIF

!  add TMS correction factor, vector case

        IF ( NSTOKES .GT. 1 ) THEN
         IF ( SUNLIGHT ) THEN

!  sunlight only

          DO N = 1, NLAYERS
           IF ( STERM_LAYERMASK_DN(N)) THEN
            DO V = 1, N_GEOMETRIES
             DO O1 = 1, NSTOKES
              ZMAT_DN(V,N,O1,1) = ZMAT_DN(V,N,O1,1) * TMS(N)
             ENDDO
            ENDDO
           ENDIF
          ENDDO

!  general case

         ELSE
          DO N = 1, NLAYERS
           IF ( STERM_LAYERMASK_DN(N)) THEN
            DO V = 1, N_GEOMETRIES
             DO O1 = 1, NSTOKES
              DO O2 = 1, NSTOKES
               ZMAT_DN(V,N,O1,O2) = ZMAT_DN(V,N,O1,O2) * TMS(N)
              ENDDO
             ENDDO
            ENDDO
           ENDIF
          ENDDO
         ENDIF

!  finish this section

        ENDIF
      ENDIF

!  Downwelling single scatter recurrence
!  -------------------------------------

      IF ( DO_DNWELLING ) THEN

!  initialize cumulative source term

        NC =  0
        DO V = 1, N_GEOMETRIES
          DO O1 = 1, NSTOKES
            SS_CUMSOURCE_DN(V,O1,NC) = ZERO
          ENDDO
        ENDDO

!  Lattice calculation

        IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN

!  initialise optical depth loop

          NSTART = 1
          NUT_PREV = NSTART - 1

!  Main loop over all output optical depths

          DO UTA = 1, N_USER_LEVELS

!  Layer index for given optical depth

            NLEVEL = UTAU_LEVEL_MASK_DN(UTA)
            NUT = NLEVEL

!  Cumulative single scatter source terms :
!      For loop over layers working downwards to NUT,
!      Get layer source terms = Exact Z-matrix * Multiplier

            DO N = NSTART, NUT
             NC = N

!  sunlight case

             IF ( SUNLIGHT ) THEN

              DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                MULTIPLIER = EMULT_DN(UM,N,IB)
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO O1 = 1, NPOLAR
                  HELP =  ZMAT_DN(V,N,O1,1) * FLUXVEC(1)
                  SS_LAYERSOURCE = HELP* MULTIPLIER
                  SS_CUMSOURCE_DN(V,O1,NC) = SS_LAYERSOURCE + &
                     T_DELT_USERM(N,UM)*SS_CUMSOURCE_DN(V,O1,NC-1)
                 ENDDO
                ENDDO
               ENDDO
              ENDDO

!  general case

             ELSE

              DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                MULTIPLIER = EMULT_DN(UM,N,IB)
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO O1 = 1, NSTOKES
                  HELP = ZERO
                  DO O2 = 1, NSTOKES
                   HELP = HELP + ZMAT_DN(V,N,O1,O2) * FLUXVEC(O2)
                  ENDDO
                  SS_LAYERSOURCE = HELP* MULTIPLIER
                  SS_CUMSOURCE_DN(V,O1,NC) = SS_LAYERSOURCE + &
                     T_DELT_USERM(N,UM)*SS_CUMSOURCE_DN(V,O1,NC-1)
                 ENDDO
                ENDDO
               ENDDO
              ENDDO

             ENDIF

!  end layer loop

            ENDDO

!  Offgrid output :
!    add additional partial layer source term = Exact Z-matrix * Multipl
!    Set final cumulative source and Correct the intensity

            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
              UT = PARTLAYERS_OUTINDEX(UTA)
              N  = PARTLAYERS_LAYERIDX(UT)

!  sunlight case

             IF ( SUNLIGHT ) THEN

              DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                MULTIPLIER = UT_EMULT_DN(UM,UT,IB)
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO O1 = 1, NSTOKES
                  HELP = ZMAT_DN(V,N,O1,1) * FLUXVEC(1)
                  SS_LAYERSOURCE = HELP * MULTIPLIER
                  SS_CUMSOURCE   = SS_CUMSOURCE_DN(V,O1,NC)
                  TR_CUMSOURCE   = T_UTDN_USERM(UT,UM) * SS_CUMSOURCE
                  FINAL_SOURCE   = TR_CUMSOURCE + SS_LAYERSOURCE
                  SSCORRECTION   = SSFLUX * FINAL_SOURCE
                  STOKES_SS(UTA,V,O1,DNIDX) = SSCORRECTION
                 ENDDO
                ENDDO
                 ENDDO
              ENDDO

!  general case

             ELSE

              DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                MULTIPLIER = UT_EMULT_DN(UM,UT,IB)
                DO IA = 1, N_USER_RELAZMS
                 V = VZA_OFFSETS(IB,UM) + IA
                 DO O1 = 1, NSTOKES
                  HELP = ZERO
                  DO O2 = 1, NSTOKES
                   HELP = HELP + ZMAT_DN(V,N,O1,O2) * FLUXVEC(O2)
                  ENDDO
                  SS_LAYERSOURCE = HELP * MULTIPLIER
                  SS_CUMSOURCE   = SS_CUMSOURCE_DN(V,O1,NC)
                  TR_CUMSOURCE   = T_UTDN_USERM(UT,UM) * SS_CUMSOURCE
                  FINAL_SOURCE   = TR_CUMSOURCE + SS_LAYERSOURCE
                  SSCORRECTION   = SSFLUX * FINAL_SOURCE
                  STOKES_SS(UTA,V,O1,DNIDX) = SSCORRECTION
                 ENDDO
                ENDDO
               ENDDO
              ENDDO

             ENDIF

!  debug code to be inserted
!                     write(98,'(a,3i4,1p2e15.7)')
!     &                 'down',UTA,IB,V,UNCORRECTED,SSCORRECTION

!  Ongrid output :
!     Set final cumulative source and correct Stokes vector

            ELSE

             DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
               DO IA = 1, N_USER_RELAZMS
                V = VZA_OFFSETS(IB,UM) + IA
                DO O1 = 1, NSTOKES
                 FINAL_SOURCE = SS_CUMSOURCE_DN(V,O1,NC)
                 SSCORRECTION = SSFLUX * FINAL_SOURCE
                 STOKES_SS(UTA,V,O1,DNIDX) = SSCORRECTION
                ENDDO
               ENDDO
              ENDDO
             ENDDO

!  debug code to be inserted above
!         write(97,'(a,3i4,1p2e15.7)')
!     &        'down',UTA,IB,V,UNCORRECTED,SSCORRECTION

            ENDIF

!  Check for updating the recursion

            IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
            NUT_PREV = NUT

!  end optical depth loop

          ENDDO

!  End lattice calculation

        ENDIF

!  Observational geometry calculation

        IF ( DO_OBSERVATION_GEOMETRY ) THEN

!  initialise optical depth loop

          NSTART = 1
          NUT_PREV = NSTART - 1

!  Main loop over all output optical depths

          DO UTA = 1, N_USER_LEVELS

!  Layer index for given optical depth

            NLEVEL = UTAU_LEVEL_MASK_DN(UTA)
            NUT = NLEVEL

!  Cumulative single scatter source terms :
!      For loop over layers working downwards to NUT,
!      Get layer source terms = Exact Z-matrix * Multiplier

            DO N = NSTART, NUT
             NC = N

!  sunlight case

             IF ( SUNLIGHT ) THEN
              DO V = 1, N_GEOMETRIES
                MULTIPLIER = EMULT_DN(LUM,N,V)
                DO O1 = 1, NPOLAR
                  HELP =  ZMAT_DN(V,N,O1,1) * FLUXVEC(1)
                  SS_LAYERSOURCE = HELP* MULTIPLIER
                  SS_CUMSOURCE_DN(V,O1,NC) = SS_LAYERSOURCE + &
                     T_DELT_USERM(N,V)*SS_CUMSOURCE_DN(V,O1,NC-1)
                ENDDO
              ENDDO

!  general case

             ELSE

              DO V = 1, N_GEOMETRIES
                MULTIPLIER = EMULT_DN(LUM,N,V)
                DO O1 = 1, NSTOKES
                  HELP = ZERO
                  DO O2 = 1, NSTOKES
                   HELP = HELP + ZMAT_DN(V,N,O1,O2) * FLUXVEC(O2)
                  ENDDO
                  SS_LAYERSOURCE = HELP* MULTIPLIER
                  SS_CUMSOURCE_DN(V,O1,NC) = SS_LAYERSOURCE + &
                     T_DELT_USERM(N,V)*SS_CUMSOURCE_DN(V,O1,NC-1)
                ENDDO
              ENDDO

             ENDIF

!  end layer loop

            ENDDO

!  Offgrid output :
!    add additional partial layer source term = Exact Z-matrix * Multiple
!    Set final cumulative source and correct the intensity

            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
              UT = PARTLAYERS_OUTINDEX(UTA)
              N  = PARTLAYERS_LAYERIDX(UT)

!  sunlight case

             IF ( SUNLIGHT ) THEN

              DO V = 1, N_GEOMETRIES
                MULTIPLIER = UT_EMULT_DN(LUM,UT,V)
                DO O1 = 1, NSTOKES
                  HELP = ZMAT_DN(V,N,O1,1) * FLUXVEC(1)
                  SS_LAYERSOURCE = HELP * MULTIPLIER
                  SS_CUMSOURCE   = SS_CUMSOURCE_DN(V,O1,NC)
                  TR_CUMSOURCE   = T_UTDN_USERM(UT,V) * SS_CUMSOURCE
                  FINAL_SOURCE   = TR_CUMSOURCE + SS_LAYERSOURCE
                  SSCORRECTION   = SSFLUX * FINAL_SOURCE
                  STOKES_SS(UTA,V,O1,DNIDX) = SSCORRECTION
                ENDDO
              ENDDO

!  general case

             ELSE

              DO V = 1, N_GEOMETRIES
                MULTIPLIER = UT_EMULT_DN(LUM,UT,V)
                DO O1 = 1, NSTOKES
                  HELP = ZERO
                  DO O2 = 1, NSTOKES
                   HELP = HELP + ZMAT_DN(V,N,O1,O2) * FLUXVEC(O2)
                  ENDDO
                  SS_LAYERSOURCE = HELP * MULTIPLIER
                  SS_CUMSOURCE   = SS_CUMSOURCE_DN(V,O1,NC)
                  TR_CUMSOURCE   = T_UTDN_USERM(UT,V) * SS_CUMSOURCE
                  FINAL_SOURCE   = TR_CUMSOURCE + SS_LAYERSOURCE
                  SSCORRECTION   = SSFLUX * FINAL_SOURCE
                  STOKES_SS(UTA,V,O1,DNIDX) = SSCORRECTION
                ENDDO
              ENDDO

             ENDIF

!  debug code to be inserted
!                     write(98,'(a,3i4,1p2e15.7)')
!     &                 'down',UTA,IB,V,UNCORRECTED,SSCORRECTION

!  Ongrid output :
!     Set final cumulative source and correct Stokes vector

            ELSE

             DO V = 1, N_GEOMETRIES
               DO O1 = 1, NSTOKES
                 FINAL_SOURCE = SS_CUMSOURCE_DN(V,O1,NC)
                 SSCORRECTION = SSFLUX * FINAL_SOURCE
                 STOKES_SS(UTA,V,O1,DNIDX) = SSCORRECTION
               ENDDO
             ENDDO

!  debug code to be inserted above
!         write(97,'(a,3i4,1p2e15.7)')
!     &        'down',UTA,IB,V,UNCORRECTED,SSCORRECTION

            ENDIF

!  Check for updating the recursion

            IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
            NUT_PREV = NUT

!  end optical depth loop

          ENDDO

!  End observationl geometry calculation

        ENDIF

!  end Downwelling clause

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_SSCORR_NADIR

!

      SUBROUTINE VLIDORTSS_FMATRICES &
        ( DO_SSCORR_TRUNCATION, DO_OBSERVATION_GEOMETRY, NLAYERS, &
          NSTOKES, N_GEOMETRIES, NGREEKMOMS, &
          NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, &
          LAYER_MAXMOMENTS, GREEKMAT_TOTAL_INPUT, SSFDEL, &
          VSIGN, VZA_OFFSETS, &
          CTHETA, STHETA, CALPHA, SALPHA, CPHI, PHI, &
          C1, S1, C2, S2, FMAT )

!  include file of dimensions and numbers
!     Cannot use bookkeeping vars file

      USE VLIDORT_PARS

      IMPLICIT NONE

!  input
!  -----

!  Additional control for the truncation

      LOGICAL, INTENT(IN) ::          DO_SSCORR_TRUNCATION

!  Observational geometry control

      LOGICAL, INTENT(IN) ::          DO_OBSERVATION_GEOMETRY

!  control integers

      INTEGER, INTENT(IN) ::          NLAYERS, NSTOKES, N_GEOMETRIES, NGREEKMOMS
      INTEGER, INTENT(IN) ::          NBEAMS, N_USER_STREAMS, N_USER_RELAZMS

!  +/- sign, offsets

      DOUBLE PRECISION, INTENT(IN) :: VSIGN
      INTEGER, INTENT(IN) ::          VZA_OFFSETS(MAXBEAMS,MAX_USER_STREAMS)

!  scattering input information

      INTEGER, INTENT(IN) ::          LAYER_MAXMOMENTS(MAXLAYERS)
      DOUBLE PRECISION, INTENT(IN) :: GREEKMAT_TOTAL_INPUT &
          ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION, INTENT(IN) :: SSFDEL( MAXLAYERS )

!  zenith angle cosines/sines, azimuth angle cosines

      DOUBLE PRECISION, INTENT(IN) :: CTHETA (MAXLAYERS,MAXBEAMS)
      DOUBLE PRECISION, INTENT(IN) :: STHETA (MAXLAYERS,MAXBEAMS)
      DOUBLE PRECISION, INTENT(IN) :: CALPHA (MAX_USER_STREAMS)
      DOUBLE PRECISION, INTENT(IN) :: SALPHA (MAX_USER_STREAMS)
      DOUBLE PRECISION, INTENT(IN) :: CPHI   (MAX_USER_RELAZMS)

!  azimuth angles. Added V2.4R. for PHI > 180 case

      DOUBLE PRECISION, INTENT(IN) :: PHI (MAX_USER_RELAZMS)

!  output
!  ------

!  Rotation Angle cosines/sines

      DOUBLE PRECISION, INTENT(OUT) ::  C1 (MAX_GEOMETRIES,MAXLAYERS)
      DOUBLE PRECISION, INTENT(OUT) ::  S1 (MAX_GEOMETRIES,MAXLAYERS)
      DOUBLE PRECISION, INTENT(OUT) ::  C2 (MAX_GEOMETRIES,MAXLAYERS)
      DOUBLE PRECISION, INTENT(OUT) ::  S2 (MAX_GEOMETRIES,MAXLAYERS)

!  F-Matrix

      DOUBLE PRECISION, INTENT(OUT) ::  FMAT ( MAX_GEOMETRIES, MAXLAYERS, 6 )

!  Local
!  -----

      INTEGER ::           IB, UM, IA, V, N, GREEKMAT_INDEX(6)
      DOUBLE PRECISION ::  COSSCAT, SINSCAT, HELP_SINSCAT, F, FT, DNL1
      DOUBLE PRECISION ::  CSIG1, CSIG2, SSIG1, SSIG2, CSIG1_2, CSIG2_2

      DOUBLE PRECISION ::  UUU(MAX_GEOMETRIES)

!  Generalized spherical functions
!    P2P2 and P2M2 added 20 March 2006

      DOUBLE PRECISION ::  P00(MAX_GEOMETRIES,2)
      DOUBLE PRECISION ::  P02(MAX_GEOMETRIES,2)
      DOUBLE PRECISION ::  P2P2(MAX_GEOMETRIES,2)
      DOUBLE PRECISION ::  P2M2(MAX_GEOMETRIES,2)

      INTEGER ::           K, L, LNEW, LOLD, ITMP
      INTEGER ::           INDEX_11, INDEX_12, INDEX_34
      INTEGER ::           INDEX_22, INDEX_33, INDEX_44
      DOUBLE PRECISION ::  DL, QROOT6, FAC1, FAC2, SQL4, SQL41
      DOUBLE PRECISION ::  TMP1, TMP2, SUM23, DIF23
      DOUBLE PRECISION ::  FL2, FLL1, PERLL4, Q, WFACT, DL1

!  Greek matrix components (local usage)

      DOUBLE PRECISION ::  GK11(MAXLAYERS), GK12(MAXLAYERS)
      DOUBLE PRECISION ::  GK44(MAXLAYERS), GK34(MAXLAYERS)
      DOUBLE PRECISION ::  GK22(MAXLAYERS), GK33(MAXLAYERS)  !! Add 10/05/10

!  indexing key

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

!  Geometrical quantities
!  ----------------------

      IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
        DO IB = 1, NBEAMS
          DO UM = 1, N_USER_STREAMS
            DO IA = 1, N_USER_RELAZMS
              V = VZA_OFFSETS(IB,UM) + IA

!  cosine scatter angle (this is valid only for non-refracting atmospher
!  VSIGN = -1 for upwelling, +1 for downwelling

              COSSCAT = VSIGN * CTHETA(1,IB) * CALPHA(UM) + &
                                STHETA(1,IB) * SALPHA(UM) * CPHI(IA)

              UUU(V)  = COSSCAT

!  Cosine Sigma 1 and 2. H/VdM, Eqs. (99)-(101)

!  a. safety
!    Watch for sin^2(scatter angle) less than zero (machine precision)
!    R. Spurr, 16 January 2006, RT SOLUTIONS Inc.

              HELP_SINSCAT = ( ONE - COSSCAT * COSSCAT )
              IF ( HELP_SINSCAT.LE.ZERO ) THEN
                SINSCAT = 1.0D-12
              ELSE
                SINSCAT = DSQRT ( HELP_SINSCAT )
              ENDIF

!  b. necessary limit analyses - Hovenier limits.
!     R. Spurr and V. Natraj, 17 January 2006

              IF ( DABS(SINSCAT) .LE. 1.0D-12 ) THEN
                CSIG1 = ZERO
                CSIG2 = ZERO
              ELSE
               IF ( STHETA(1,IB) .EQ. ZERO ) THEN
                CSIG1 = -  CPHI(IA)
               ELSE
                CSIG1   = ( - VSIGN * CALPHA(UM) + CTHETA(1,IB)*COSSCAT ) &
                          / SINSCAT / STHETA(1,IB)
               ENDIF
               IF ( SALPHA(UM) .EQ. ZERO ) THEN
                 CSIG2 = -  CPHI(IA)
               ELSE
                 CSIG2   = ( - CTHETA(1,IB) + VSIGN*CALPHA(UM)*COSSCAT ) &
                            / SINSCAT / SALPHA(UM)
               ENDIF
              ENDIF

!  Debugged again, 06 October 2010. Everything OK.
!       write(*,'(2i4,f6.2,1p3e15.7)')V,1,VSIGN,CSIG1,CSIG2,PHI(IA)
!       write(*,'(i4,2f6.1,1p6e15.7)')V,VSIGN,PHI(IA),
!     &   SALPHA(UM),CALPHA(UM),STHETA(1,IB),CTHETA(1,IB),CSIG1,CSIG2

!  These lines are necessary to avoid bad values

              IF ( CSIG1 .GT. ONE  ) CSIG1 = ONE
              IF ( CSIG1 .LT. -ONE ) CSIG1 = -ONE

              IF ( CSIG2 .GT. ONE  ) CSIG2 = ONE
              IF ( CSIG2 .LT. -ONE ) CSIG2 = -ONE

!  output, H/VdM, Eqs. (89)-(94)
!    Rotation sines and cosines. Same for all layers

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

!  For relazm in [180,360), need sign reversal for S1 and S2
!  See H/VdM, Eqs. 94-95. V. Natraj and R. Spurr, 01 May 2009.

              C1(V,1) = CSIG1_2 * CSIG1 - ONE
              C2(V,1) = CSIG2_2 * CSIG2 - ONE

              IF (PHI(IA) .LE. 180.D0) THEN
                S1(V,1) = CSIG1_2 * SSIG1
                S2(V,1) = CSIG2_2 * SSIG2
              ELSE
                S1(V,1) = -CSIG1_2 * SSIG1
                S2(V,1) = -CSIG2_2 * SSIG2
              ENDIF

!  Copy to all other layers

              DO N = 2, NLAYERS
                C1(V,N) = C1(V,1)
                S1(V,N) = S1(V,1)
                C2(V,N) = C2(V,1)
                S2(V,N) = S2(V,1)
              ENDDO
!      write(8,'(i5,1p4e15.7)')V,C1(V,1),S1(V,1),C2(V,1),S2(V,1)

!  End geometry loops

            ENDDO
          ENDDO
        ENDDO

      ELSE

        DO V = 1, N_GEOMETRIES

!  cosine scatter angle (this is valid only for non-refracting atmospher
!  VSIGN = -1 for upwelling, +1 for downwelling

          COSSCAT = VSIGN * CTHETA(1,V) * CALPHA(V) + &
                            STHETA(1,V) * SALPHA(V) * CPHI(V)

          UUU(V)  = COSSCAT

!  Cosine Sigma 1 and 2. H/VdM, Eqs. (99)-(101)

!  a. safety
!    Watch for sin^2(scatter angle) less than zero (machine precision)
!    R. Spurr, 16 January 2006, RT SOLUTIONS Inc.

          HELP_SINSCAT = ( ONE - COSSCAT * COSSCAT )
          IF ( HELP_SINSCAT.LE.ZERO ) THEN
            SINSCAT = 1.0D-12
          ELSE
            SINSCAT = DSQRT ( HELP_SINSCAT )
          ENDIF

!  b. necessary limit analyses - Hovenier limits.
!     R. Spurr and V. Natraj, 17 January 2006

          IF ( DABS(SINSCAT) .LE. 1.0D-12 ) THEN
            CSIG1 = ZERO
            CSIG2 = ZERO
          ELSE
           IF ( STHETA(1,V) .EQ. ZERO ) THEN
            CSIG1 = -  CPHI(V)
           ELSE
            CSIG1   = ( - VSIGN * CALPHA(V) + CTHETA(1,V)*COSSCAT ) &
                      / SINSCAT / STHETA(1,V)
           ENDIF
           IF ( SALPHA(V) .EQ. ZERO ) THEN
             CSIG2 = -  CPHI(V)
           ELSE
             CSIG2   = ( - CTHETA(1,V) + VSIGN*CALPHA(V)*COSSCAT ) &
                        / SINSCAT / SALPHA(V)
           ENDIF
          ENDIF

!  Debugged again, 06 October 2010. Everything OK.
!       write(*,'(2i4,f6.2,1p3e15.7)')V,1,VSIGN,CSIG1,CSIG2,PHI(IA)
!       write(*,'(i4,2f6.1,1p6e15.7)')V,VSIGN,PHI(IA),
!     &   SALPHA(UM),CALPHA(UM),STHETA(1,IB),CTHETA(1,IB),CSIG1,CSIG2

!  These lines are necessary to avoid bad values

          IF ( CSIG1 .GT. ONE  ) CSIG1 = ONE
          IF ( CSIG1 .LT. -ONE ) CSIG1 = -ONE

          IF ( CSIG2 .GT. ONE  ) CSIG2 = ONE
          IF ( CSIG2 .LT. -ONE ) CSIG2 = -ONE

!  output, H/VdM, Eqs. (89)-(94)
!    Rotation sines and cosines. Same for all layers

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

!  For relazm in [180,360), need sign reversal for S1 and S2
!  See H/VdM, Eqs. 94-95. V. Natraj and R. Spurr, 01 May 2009.

          C1(V,1) = CSIG1_2 * CSIG1 - ONE
          C2(V,1) = CSIG2_2 * CSIG2 - ONE

          IF (PHI(V) .LE. 180.D0) THEN
            S1(V,1) = CSIG1_2 * SSIG1
            S2(V,1) = CSIG2_2 * SSIG2
          ELSE
            S1(V,1) = -CSIG1_2 * SSIG1
            S2(V,1) = -CSIG2_2 * SSIG2
          ENDIF

!  Copy to all other layers

          DO N = 2, NLAYERS
            C1(V,N) = C1(V,1)
            S1(V,N) = S1(V,1)
            C2(V,N) = C2(V,1)
            S2(V,N) = S2(V,1)
          ENDDO

!         write(8,'(i5,1p4e15.7)')V,C1(V,1),S1(V,1),C2(V,1),S2(V,1)

!  End geometry loop

        ENDDO

      ENDIF

!  F-matrices
!  ----------

      QROOT6 = -0.25D0 * DSQRT(6.0D0)

! initialise F-matrix

      DO V = 1, N_GEOMETRIES
        DO K = 1, 6
          DO N = 1, NLAYERS
            FMAT(V,N,K) = ZERO
          END DO
        END DO
      END DO

!  Start loop over the coefficient index l
!  first update generalized spherical functions, then calculate coefs.
!  lold and lnew are pointer-like indices used in recurrence

      LNEW = 1
      LOLD = 2

      DO L = 0, NGREEKMOMS

        DL   = DBLE(L)
        DL1  = DL - ONE

!  Set the local Greek matrix elements that you need
!   44 and 34 are not required with natural sunlight (default here)
!   22 and 33 required for non-Mie spheroidal particles

        IF ( DO_SSCORR_TRUNCATION ) THEN
          DO N = 1, NLAYERS
           DNL1 = DBLE(2*L + 1 )
           F    = SSFDEL(N) * DNL1
           FT   = ONE - SSFDEL(N)
           GK11(N) = (GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(1))-F)/FT
!mick fix 1/21/2013 - added IF structure and ELSE section
           IF (NSTOKES .GT. 1) THEN
             GK12(N) =  GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(3))   /FT
             GK44(N) = (GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(6))-F)/FT
             GK34(N) =  GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(5))   /FT
             GK22(N) = (GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(2))-F)/FT
             GK33(N) = (GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(4))-F)/FT
           ELSE
             GK12(N) = ZERO
             GK44(N) = ZERO
             GK34(N) = ZERO
             GK22(N) = ZERO
             GK33(N) = ZERO
           END IF
          ENDDO
        ELSE
          DO N = 1, NLAYERS
           GK11(N) = GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(1))
!mick fix 1/21/2013 - added IF structure and ELSE section
           IF (NSTOKES .GT. 1) THEN
             GK12(N) = GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(3))
             GK44(N) = GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(6))
             GK34(N) = GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(5))
             GK22(N) = GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(2))
             GK33(N) = GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(4))
           ELSE
             GK12(N) = ZERO
             GK44(N) = ZERO
             GK34(N) = ZERO
             GK22(N) = ZERO
             GK33(N) = ZERO
           ENDIF
          ENDDO
        ENDIF

!  First moment

        IF ( L .EQ. 0 ) THEN

!  Adding paper Eqs. (76) and (77) with m=0
!   Additional functions P2M2 and P2P2 zero for m = 0

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

! Adding paper Eq. (81) with m=0

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

!  introduce the P2P2 and P2M2 functions for L = 2

          DO V = 1, N_GEOMETRIES
            P2P2(V,LOLD)= 0.25D0*(ONE+UUU(V))*(ONE+UUU(V))
            P2M2(V,LOLD)= 0.25D0*(ONE-UUU(V))*(ONE-UUU(V))
          ENDDO

        ELSE IF ( L .GT. 2) THEN

! Adding paper Eq. (82) with m=0

          SQL4  = SQL41
          SQL41 = DSQRT(DL*DL-FOUR)
          TMP1  = (TWO*DL-ONE)/SQL41
          TMP2  = SQL4/SQL41
          DO V = 1, N_GEOMETRIES
            P02(V,LOLD) = TMP1*UUU(V)*P02(V,LNEW) - TMP2*P02(V,LOLD)
          END DO

!  introduce the P2P2 and P2M2 functions for L > 2

          FL2 = TWO * DL - ONE
          FLL1 = DL * DL1
          PERLL4=ONE/(DL1*SQL41**2)
          Q     = DL  * ( DL1*DL1 - FOUR)
          DO V = 1, N_GEOMETRIES
           WFACT = FL2 * ( FLL1 * UUU(V) - FOUR )
           P2P2(V,LOLD) = (WFACT*P2P2(V,LNEW) - &
                Q*P2P2(V,LOLD)) * PERLL4
           WFACT = FL2 * ( FLL1 * UUU(V) + FOUR )
           P2M2(V,LOLD) = (WFACT*P2M2(V,LNEW) - &
                Q*P2M2(V,LOLD)) * PERLL4
          ENDDO

        END IF

! Switch indices so that lnew indicates the function with
! the present index value l, this mechanism prevents swapping
! of entire arrays.

        ITMP = LNEW
        LNEW = LOLD
        LOLD = ITMP

! Now add the l-th term to the scattering matrix.
! See de Haan et al. (1987) Eqs. (68)-(73).

! Section for randomly-oriented spheroids, added 20 march 2006
!  R. Spurr and V. Natraj

        DO N = 1, NLAYERS
         IF ( L.LE.LAYER_MAXMOMENTS(N) ) THEN
          DO V = 1, N_GEOMETRIES

           FMAT(V,N,INDEX_11) = FMAT(V,N,INDEX_11) + GK11(N)*P00(V,LNEW)
           FMAT(V,N,INDEX_12) = FMAT(V,N,INDEX_12) + GK12(N)*P02(V,LNEW)

!  This section adds F-matrix terms needed for randomly-oriented spheroi
           SUM23 = GK22(N) + GK33(N)
           DIF23 = GK22(N) - GK33(N)
           FMAT(V,N,INDEX_22) = FMAT(V,N,INDEX_22) + SUM23*P2P2(V,LNEW)
           FMAT(V,N,INDEX_33) = FMAT(V,N,INDEX_33) + DIF23*P2M2(V,LNEW)
!  end of new section---------------------------------------------------

           FMAT(V,N,INDEX_44) = FMAT(V,N,INDEX_44) + GK44(N)*P00(V,LNEW)
           FMAT(V,N,INDEX_34) = FMAT(V,N,INDEX_34) + GK34(N)*P02(V,LNEW)

          ENDDO
         ENDIF
        END DO

!  end moment loop

      END DO

!   This must be done AFTER the moment loop.

      DO V = 1, N_GEOMETRIES
        DO N = 1, NLAYERS
           FMAT(V,N,INDEX_22) = HALF*(FMAT(V,N,INDEX_22)+ &
               FMAT(V,N,INDEX_33))
           FMAT(V,N,INDEX_33) = FMAT(V,N,INDEX_22)- &
              FMAT(V,N,INDEX_33)
        END DO
      END DO

! Remember for Mie scattering : F11 = F22 and F33 = F44
!  This code is no longer required, as we have introduced code now
!   for randomly oriented spheroids. The symmetry should still of
!   course be present for the Mie particles, so this will be a
!   check on the new code.
!   R. Spurr and V. Natraj, 20 March 2006

!      DO V = 1, N_GEOMETRIES
!        DO N = 1, NLAYERS
!          FMAT(V,N,INDEX_22) = FMAT(V,N,INDEX_11)
!          FMAT(V,N,INDEX_33) = FMAT(V,N,INDEX_44)
!        END DO
!      END DO

!  finish

      RETURN
      END SUBROUTINE VLIDORTSS_FMATRICES

!

      SUBROUTINE VLIDORTSS_FMATRICES_MULTI &
        ( DO_SSCORR_TRUNCATION, DO_OBSERVATION_GEOMETRY, NLAYERS, &
          NSTOKES, N_GEOMETRIES, NGREEKMOMS, &
          NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, &
          LAYER_MAXMOMENTS, GREEKMAT_TOTAL_INPUT, SSFDEL, &
          VSIGN, VZA_OFFSETS, &
          CTHETA, STHETA, CALPHA, SALPHA, CPHI, PHI, &
          C1, S1, C2, S2, FMAT )

!  F-matrix routine with multiple-solar zenith angles (one for each laye
!  Applicable with the Refractive Geometry case.

!  include file of dimensions and numbers

      USE VLIDORT_PARS

      IMPLICIT NONE

!  input
!  -----

!  Additional control for the truncation

      LOGICAL, INTENT(IN) ::           DO_SSCORR_TRUNCATION

!  Observational geometry control

      LOGICAL, INTENT(IN) ::           DO_OBSERVATION_GEOMETRY

!  control integers

      INTEGER, INTENT(IN) ::           NLAYERS, NSTOKES, N_GEOMETRIES, NGREEKMOMS
      INTEGER, INTENT(IN) ::           NBEAMS, N_USER_STREAMS,  N_USER_RELAZMS

!  +/- sign, offsets

      DOUBLE PRECISION, INTENT(IN) ::  VSIGN
      INTEGER, INTENT(IN) ::           VZA_OFFSETS(MAXBEAMS,MAX_USER_STREAMS)

!  scattering input information

      INTEGER, INTENT(IN) ::           LAYER_MAXMOMENTS(MAXLAYERS)
      DOUBLE PRECISION, INTENT(IN) ::  GREEKMAT_TOTAL_INPUT &
          ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION, INTENT(IN) ::  SSFDEL( MAXLAYERS )

!  zenith angle cosines/sines, azimuth angle cosines

      DOUBLE PRECISION, INTENT(IN) ::  CTHETA (MAXLAYERS,MAXBEAMS)
      DOUBLE PRECISION, INTENT(IN) ::  STHETA (MAXLAYERS,MAXBEAMS)
      DOUBLE PRECISION, INTENT(IN) ::  CALPHA (MAX_USER_STREAMS)
      DOUBLE PRECISION, INTENT(IN) ::  SALPHA (MAX_USER_STREAMS)
      DOUBLE PRECISION, INTENT(IN) ::  CPHI   (MAX_USER_RELAZMS)

!  azimuth angles. Added V2.4R. for PHI > 180 case

      DOUBLE PRECISION, INTENT(IN) ::  PHI (MAX_USER_RELAZMS)

!  output
!  ------

!  Rotation Angle cosines/sines

      DOUBLE PRECISION, INTENT(OUT) ::  C1 (MAX_GEOMETRIES,MAXLAYERS)
      DOUBLE PRECISION, INTENT(OUT) ::  S1 (MAX_GEOMETRIES,MAXLAYERS)
      DOUBLE PRECISION, INTENT(OUT) ::  C2 (MAX_GEOMETRIES,MAXLAYERS)
      DOUBLE PRECISION, INTENT(OUT) ::  S2 (MAX_GEOMETRIES,MAXLAYERS)

!  F-Matrix

      DOUBLE PRECISION, INTENT(OUT) ::   FMAT ( MAX_GEOMETRIES, MAXLAYERS, 6 )

!  Local
!  -----

      INTEGER ::           IB, UM, IA, V, N, GREEKMAT_INDEX(6)
      DOUBLE PRECISION ::  COSSCAT, SINSCAT, HELP_SINSCAT, F, FT, DNL1
      DOUBLE PRECISION ::  CSIG1, CSIG2, SSIG1, SSIG2, CSIG1_2, CSIG2_2

!  now have layer dimensioning in the help arrays

      DOUBLE PRECISION ::  UUU(MAX_GEOMETRIES,MAXLAYERS)

!  Generalized spherical functions
!    P2P2 and P2M2 added 20 March 2006

      DOUBLE PRECISION ::  P00(MAX_GEOMETRIES,MAXLAYERS,2)
      DOUBLE PRECISION ::  P02(MAX_GEOMETRIES,MAXLAYERS,2)
      DOUBLE PRECISION ::  P2P2(MAX_GEOMETRIES,MAXLAYERS,2)
      DOUBLE PRECISION ::  P2M2(MAX_GEOMETRIES,MAXLAYERS,2)

      INTEGER ::           K, L, LNEW, LOLD, ITMP
      INTEGER ::           INDEX_11, INDEX_12, INDEX_34
      INTEGER ::           INDEX_22, INDEX_33, INDEX_44
      DOUBLE PRECISION ::  DL, QROOT6, FAC1, FAC2, SQL4, SQL41
      DOUBLE PRECISION ::  TMP1, TMP2, SUM23, DIF23
      DOUBLE PRECISION ::  FL2, FLL1, PERLL4, Q, WFACT, DL1
      DOUBLE PRECISION ::  GK11(MAXLAYERS), GK12(MAXLAYERS)
      DOUBLE PRECISION ::  GK44(MAXLAYERS), GK34(MAXLAYERS)
      DOUBLE PRECISION ::  GK22(MAXLAYERS), GK33(MAXLAYERS)

!  indexing key

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

!  Geometrical quantities
!  ----------------------

      IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
        DO IB = 1, NBEAMS
          DO UM = 1, N_USER_STREAMS
            DO IA = 1, N_USER_RELAZMS
              V = VZA_OFFSETS(IB,UM) + IA
              DO N = 1, NLAYERS

!  cosine scatter angle (this is valid only for non-refracting atmospher
!  VSIGN = -1 for upwelling, +1 for downwelling

                COSSCAT = VSIGN * CTHETA(N,IB) * CALPHA(UM) + &
                                  STHETA(N,IB) * SALPHA(UM) * CPHI(IA)

                UUU(V,N)  = COSSCAT

!  Cosine Sigma 1 and 2. H/VdM, Eqs. (99)-(101)

!  a. safety
!    Watch for sin^2(scatter angle) less than zero (machine precision)
!    R. Spurr, 16 January 2006, RT SOLUTIONS Inc.

                HELP_SINSCAT = ( ONE - COSSCAT * COSSCAT )
                IF ( HELP_SINSCAT.LE.ZERO ) THEN
                  SINSCAT = 1.0D-12
                ELSE
                  SINSCAT = DSQRT ( HELP_SINSCAT )
                ENDIF

!  b. necessary limit analyses - Hovenier limits.
!     R. Spurr and V. Natraj, 17 January 2006

                IF ( DABS(SINSCAT) .LE. 1.0D-12 ) THEN
                 CSIG1 = ZERO
                 CSIG2 = ZERO
                ELSE
                 IF ( STHETA(1,IB) .EQ. ZERO ) THEN
                  CSIG1 = -  CPHI(IA)
                 ELSE
                  CSIG1 = ( - VSIGN*CALPHA(UM) + CTHETA(N,IB)*COSSCAT ) &
                               / SINSCAT / STHETA(N,IB)
                 ENDIF
                 IF ( SALPHA(UM) .EQ. ZERO ) THEN
                  CSIG2 = -  CPHI(IA)
                 ELSE
                  CSIG2 = ( - CTHETA(N,IB) + VSIGN*CALPHA(UM)*COSSCAT ) &
                               / SINSCAT / SALPHA(UM)
                 ENDIF
                ENDIF

!  These lines are necessary to avoid bad values

                IF ( CSIG1 .GT. ONE  ) CSIG1 = ONE
                IF ( CSIG1 .LT. -ONE ) CSIG1 = -ONE

                IF ( CSIG2 .GT. ONE  ) CSIG2 = ONE
                IF ( CSIG2 .LT. -ONE ) CSIG2 = -ONE

!  output, H/VdM, Eqs. (89)-(94)

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

!  For relazm in [180,360), need sign reversal for S1 and S2
!  See H/VdM, Eqs. 94-95. V. Natraj and R. Spurr, 01 May 2009.

                C1(V,N) = CSIG1_2 * CSIG1 - ONE
                C2(V,N) = CSIG2_2 * CSIG2 - ONE

                IF (PHI(IA) .LE. 180.D0) THEN
                  S1(V,N) = CSIG1_2 * SSIG1
                  S2(V,N) = CSIG2_2 * SSIG2
                ELSE
                  S1(V,N) = -CSIG1_2 * SSIG1
                  S2(V,N) = -CSIG2_2 * SSIG2
                ENDIF

!  End layer loop

              ENDDO

!  End geometry loops

            ENDDO
          ENDDO
        ENDDO

      ELSE

        DO V = 1, N_GEOMETRIES
          DO N = 1, NLAYERS

!  cosine scatter angle (this is valid only for non-refracting atmospher
!  VSIGN = -1 for upwelling, +1 for downwelling

            COSSCAT = VSIGN * CTHETA(N,V) * CALPHA(V) + &
                              STHETA(N,V) * SALPHA(V) * CPHI(V)

            UUU(V,N)  = COSSCAT

!  Cosine Sigma 1 and 2. H/VdM, Eqs. (99)-(101)

!  a. safety
!    Watch for sin^2(scatter angle) less than zero (machine precision)
!    R. Spurr, 16 January 2006, RT SOLUTIONS Inc.

            HELP_SINSCAT = ( ONE - COSSCAT * COSSCAT )
            IF ( HELP_SINSCAT.LE.ZERO ) THEN
              SINSCAT = 1.0D-12
            ELSE
              SINSCAT = DSQRT ( HELP_SINSCAT )
            ENDIF

!  b. necessary limit analyses - Hovenier limits.
!     R. Spurr and V. Natraj, 17 January 2006

            IF ( DABS(SINSCAT) .LE. 1.0D-12 ) THEN
             CSIG1 = ZERO
             CSIG2 = ZERO
            ELSE
             IF ( STHETA(1,V) .EQ. ZERO ) THEN
              CSIG1 = -  CPHI(V)
             ELSE
              CSIG1 = ( - VSIGN*CALPHA(V) + CTHETA(N,V)*COSSCAT ) &
                           / SINSCAT / STHETA(N,V)
             ENDIF
             IF ( SALPHA(V) .EQ. ZERO ) THEN
              CSIG2 = -  CPHI(V)
             ELSE
              CSIG2 = ( - CTHETA(N,V) + VSIGN*CALPHA(V)*COSSCAT ) &
                           / SINSCAT / SALPHA(V)
             ENDIF
            ENDIF

!  These lines are necessary to avoid bad values

            IF ( CSIG1 .GT. ONE  ) CSIG1 = ONE
            IF ( CSIG1 .LT. -ONE ) CSIG1 = -ONE

            IF ( CSIG2 .GT. ONE  ) CSIG2 = ONE
            IF ( CSIG2 .LT. -ONE ) CSIG2 = -ONE

!  output, H/VdM, Eqs. (89)-(94)

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

!  For relazm in [180,360), need sign reversal for S1 and S2
!  See H/VdM, Eqs. 94-95. V. Natraj and R. Spurr, 01 May 2009.

            C1(V,N) = CSIG1_2 * CSIG1 - ONE
            C2(V,N) = CSIG2_2 * CSIG2 - ONE

            IF (PHI(V) .LE. 180.D0) THEN
              S1(V,N) = CSIG1_2 * SSIG1
              S2(V,N) = CSIG2_2 * SSIG2
            ELSE
              S1(V,N) = -CSIG1_2 * SSIG1
              S2(V,N) = -CSIG2_2 * SSIG2
            ENDIF

!  End layer loop

          ENDDO

!  End geometry loop
        ENDDO

      ENDIF

!  F-matrices
!  ----------

      QROOT6 = -0.25D0 * DSQRT(6.0D0)

! initialise F-matrix

      DO V = 1, N_GEOMETRIES
        DO K = 1, 6
          DO N = 1, NLAYERS
            FMAT(V,N,K) = ZERO
          END DO
        END DO
      END DO

!  Start loop over the coefficient index l
!  first update generalized spherical functions, then calculate coefs.
!  lold and lnew are pointer-like indices used in recurrence

      LNEW = 1
      LOLD = 2

      DO L = 0, NGREEKMOMS

        DL   = DBLE(L)
        DL1  = DL - ONE

!  Set the local Greek matrix elements that you need
!   44 and 34 are not required with natural sunlight (default here)
!   22 and 33 required for non-Mie spheroidal particles

        IF ( DO_SSCORR_TRUNCATION ) THEN
          DO N = 1, NLAYERS
           DNL1 = DBLE(2*L + 1 )
           F    = SSFDEL(N) * DNL1
           FT   = ONE - SSFDEL(N)
           GK11(N) = (GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(1))-F)/FT
!mick fix 1/21/2013 - added IF structure and ELSE section
           IF (NSTOKES .GT. 1) THEN
             GK12(N) =  GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(3))   /FT
             GK44(N) = (GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(6))-F)/FT
             GK34(N) =  GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(5))   /FT
             GK22(N) = (GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(2))-F)/FT
             GK33(N) = (GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(4))-F)/FT
           ELSE
             GK12(N) = ZERO
             GK44(N) = ZERO
             GK34(N) = ZERO
             GK22(N) = ZERO
             GK33(N) = ZERO
           END IF
          ENDDO
        ELSE
          DO N = 1, NLAYERS
           GK11(N) = GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(1))
!mick fix 1/21/2013 - added IF structure and ELSE section
           IF (NSTOKES .GT. 1) THEN
             GK12(N) = GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(3))
             GK44(N) = GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(6))
             GK34(N) = GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(5))
             GK22(N) = GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(2))
             GK33(N) = GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(4))
           ELSE
             GK12(N) = ZERO
             GK44(N) = ZERO
             GK34(N) = ZERO
             GK22(N) = ZERO
             GK33(N) = ZERO
           ENDIF
          ENDDO
        ENDIF

        IF ( L .EQ. 0 ) THEN

!  Adding paper Eqs. (76) and (77) with m=0

          DO V = 1, N_GEOMETRIES
            DO N = 1, NLAYERS
              P00(V,N,LOLD) = ONE
              P00(V,N,LNEW) = ZERO
              P02(V,N,LOLD) = ZERO
              P02(V,N,LNEW) = ZERO
              P2P2(V,N,LOLD) = ZERO
              P2P2(V,N,LNEW) = ZERO
              P2M2(V,N,LOLD) = ZERO
              P2M2(V,N,LNEW) = ZERO
            END DO
          END DO

        ELSE

          FAC1 = (TWO*DL-ONE)/DL
          FAC2 = DL1/DL

! Adding paper Eq. (81) with m=0

          DO V = 1, N_GEOMETRIES
            DO N = 1, NLAYERS
              P00(V,N,LOLD) = FAC1*UUU(V,N)*P00(V,N,LNEW) &
                              - FAC2*P00(V,N,LOLD)
            END DO
          END DO

        END IF

        IF ( L .EQ. 2 ) THEN

! Adding paper Eq. (78)
! sql4 contains the factor dsqrt((l+1)*(l+1)-4) needed in
! the recurrence Eqs. (81) and (82)

          DO V = 1, N_GEOMETRIES
            DO N = 1, NLAYERS
              P02(V,N,LOLD) = QROOT6*(ONE-UUU(V,N)*UUU(V,N))
              P02(V,N,LNEW) = ZERO
            END DO
          END DO
          SQL41 = ZERO

!  introduce the P2P2 and P2M2 functions for L = 2

          DO V = 1, N_GEOMETRIES
            DO N = 1, NLAYERS
              P2P2(V,N,LOLD)= 0.25D0*(ONE+UUU(V,N))*(ONE+UUU(V,N))
              P2M2(V,N,LOLD)= 0.25D0*(ONE-UUU(V,N))*(ONE-UUU(V,N))
            ENDDO
          ENDDO

        ELSE IF ( L .GT. 2) THEN

! Adding paper Eq. (82) with m=0

          SQL4  = SQL41
          SQL41 = DSQRT(DL*DL-FOUR)
          TMP1  = (TWO*DL-ONE)/SQL41
          TMP2  = SQL4/SQL41
          DO V = 1, N_GEOMETRIES
            DO N = 1, NLAYERS
              P02(V,N,LOLD) = TMP1*UUU(V,N)*P02(V,N,LNEW) &
                              - TMP2*P02(V,N,LOLD)
            END DO
          END DO

!  introduce the P2P2 and P2M2 functions for L > 2

          FL2 = TWO * DL - ONE
          FLL1 = DL * DL1
          PERLL4=ONE/(DL1*SQL41**2)
          Q     = DL  * ( DL1*DL1 - FOUR)
          DO V = 1, N_GEOMETRIES
           DO N = 1, NLAYERS
            WFACT = FL2 * ( FLL1*UUU(V,N) - FOUR )
            P2P2(V,N,LOLD) = (WFACT*P2P2(V,N,LNEW) - &
                Q*P2P2(V,N,LOLD)) * PERLL4
            WFACT = FL2 * ( FLL1 * UUU(V,N) + FOUR )
            P2M2(V,N,LOLD) = (WFACT*P2M2(V,N,LNEW) - &
                Q*P2M2(V,N,LOLD)) * PERLL4
           ENDDO
          ENDDO

        END IF

! Switch indices so that lnew indicates the function with
! the present index value l, this mechanism prevents swapping
! of entire arrays.

        ITMP = LNEW
        LNEW = LOLD
        LOLD = ITMP

! Now add the l-th term to the scattering matrix.
! See de Haan et al. (1987) Eqs. (68)-(73).

! Section for randomly-oriented spheroids, added 20 march 2006
!  R. Spurr and V. Natraj

        DO N = 1, NLAYERS
         IF ( L.LE.LAYER_MAXMOMENTS(N) ) THEN
          DO V = 1, N_GEOMETRIES

           FMAT(V,N,INDEX_11) = FMAT(V,N,INDEX_11)+GK11(N)*P00(V,N,LNEW)
           FMAT(V,N,INDEX_12) = FMAT(V,N,INDEX_12)+GK12(N)*P02(V,N,LNEW)

!  This section adds F-matrix terms needed for randomly-oriented spheroi
           SUM23 = GK22(N) + GK33(N)
           DIF23 = GK22(N) - GK33(N)
           FMAT(V,N,INDEX_22) = FMAT(V,N,INDEX_22)+SUM23*P2P2(V,N,LNEW)
           FMAT(V,N,INDEX_33) = FMAT(V,N,INDEX_33)+DIF23*P2M2(V,N,LNEW)
!  end of new section---------------------------------------------------

           FMAT(V,N,INDEX_44) = FMAT(V,N,INDEX_44)+GK44(N)*P00(V,N,LNEW)
           FMAT(V,N,INDEX_34) = FMAT(V,N,INDEX_34)+GK34(N)*P02(V,N,LNEW)

          ENDDO
         ENDIF
        END DO

!  end moment loop

      END DO

!   This must be done AFTER the moment loop.

      DO V = 1, N_GEOMETRIES
        DO N = 1, NLAYERS
           FMAT(V,N,INDEX_22) = HALF*(FMAT(V,N,INDEX_22)+ &
               FMAT(V,N,INDEX_33))
           FMAT(V,N,INDEX_33) = FMAT(V,N,INDEX_22)- &
              FMAT(V,N,INDEX_33)
        END DO
      END DO

! Remember for Mie scattering : F11 = F22 and F33 = F44
!  This code is no longer required, as we have introduced code now
!   for randomly oriented spheroids. The symmetry should still of
!   course be present for the Mie particles, so this will be a
!   check on the new code.
!   R. Spurr and V. Natraj,, 20 March 2006

!      DO V = 1, N_GEOMETRIES
!        DO N = 1, NLAYERS
!          FMAT(V,N,INDEX_22) = FMAT(V,N,INDEX_11)
!          FMAT(V,N,INDEX_33) = FMAT(V,N,INDEX_44)
!        END DO
!      END DO

!  finish

      RETURN
      END SUBROUTINE VLIDORTSS_FMATRICES_MULTI

!

      SUBROUTINE VLIDORT_SSCORR_OUTGOING ( &
        SSFLUX, &
        DO_SSCORR_TRUNCATION, DO_DELTAM_SCALING, &
        DO_UPWELLING, DO_DNWELLING, DO_OBSERVATION_GEOMETRY, &
        NSTOKES, NLAYERS, &
        NGREEK_MOMENTS_INPUT, &
        N_USER_RELAZMS, N_USER_LEVELS, &
        EARTH_RADIUS, HEIGHT_GRID, &
        OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT, &
        DO_TOA_CONTRIBS, &
        FLUXVEC, NBEAMS, &
        N_USER_STREAMS, LAYER_MAXMOMENTS, &
        USER_VZANGLES_ADJUST, SZANGLES_ADJUST, &
        USER_RELAZMS_ADJUST, &
        UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, &
        STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, &
        N_GEOMETRIES, VZA_OFFSETS, &
        DELTAU_VERT, TRUNC_FACTOR, &
        N_PARTLAYERS, NFINELAYERS, &
        PARTLAYERS_LAYERIDX, PARTAU_VERT, &
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, &
        UP_MULTIPLIERS, DN_MULTIPLIERS, &
        UP_MULTIPLIERS_UT, DN_MULTIPLIERS_UT, &
        UP_LOSTRANS, DN_LOSTRANS, &
        UP_LOSTRANS_UT, DN_LOSTRANS_UT, &
        ZMAT_UP, ZMAT_DN, TMS, &
        SSFDEL, SS_CUMSOURCE_UP, SS_CUMSOURCE_DN, &
        BOA_ATTN, STOKES_SS, SS_CONTRIBS, &
        FAIL, MESSAGE, TRACE)

!  SINGLE SCATTER EXACT CALCULATION FOR THE OUTGOING LOS
!         - NEW FOR VERSION 2.2

!   PROGRAMMED BY R. SPURR, RT SOLUTIONS INC.
!    FIRST DRAFT, JANUARY 30TH 2007.
!    VALIDATED AGAINST TOMRAD. 29 MARCH 2007.
!   PARTIAL LAYER OUTPUT ADDED SEPTEMBER 2007

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  EXTENSION AND CHANGES 01-05 OCTOBER 2010, R. SPURR
!   - CORRECTIONS FOR THE PHI > 180 CASE
!   - ADDITION OF NON-MIE (SPHEROIDAL PARTICLE) ENTRIES TO F-MATRIX
!   - GENERALIZATION OF Z-MATRIX TO NON-SUNLIGHT SITUATIONS.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      USE VLIDORT_PARS
      USE VLIDORT_GEOMETRY

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT (IN) ::  SSFLUX

      LOGICAL, INTENT (IN) ::           DO_SSCORR_TRUNCATION
      LOGICAL, INTENT (IN) ::           DO_DELTAM_SCALING
      LOGICAL, INTENT (IN) ::           DO_UPWELLING
      LOGICAL, INTENT (IN) ::           DO_DNWELLING
      LOGICAL, INTENT (IN) ::           DO_OBSERVATION_GEOMETRY
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           NGREEK_MOMENTS_INPUT
      INTEGER, INTENT (IN) ::           N_USER_RELAZMS
      INTEGER, INTENT (IN) ::           N_USER_LEVELS
      DOUBLE PRECISION, INTENT (IN) ::  EARTH_RADIUS
      DOUBLE PRECISION, INTENT (IN) ::  HEIGHT_GRID ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  OMEGA_TOTAL_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  GREEKMAT_TOTAL_INPUT &
          ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
      LOGICAL, INTENT (IN) ::           DO_TOA_CONTRIBS
      DOUBLE PRECISION, INTENT (IN) ::  FLUXVEC ( MAXSTOKES )
      INTEGER, INTENT (IN) ::           NBEAMS
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      INTEGER, INTENT (IN) ::           LAYER_MAXMOMENTS ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  USER_VZANGLES_ADJUST &
          ( MAX_USER_VZANGLES )
      DOUBLE PRECISION, INTENT (IN) ::  SZANGLES_ADJUST &
          ( MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS )
      DOUBLE PRECISION, INTENT (IN) ::  USER_RELAZMS_ADJUST &
          ( MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS )
      INTEGER, INTENT (IN) ::           UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )
      LOGICAL, INTENT (IN) ::           STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::           STERM_LAYERMASK_DN ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           N_GEOMETRIES
      INTEGER, INTENT (IN) ::           VZA_OFFSETS &
          ( MAX_SZANGLES, MAX_USER_VZANGLES )
      DOUBLE PRECISION, INTENT (IN) ::  DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  TRUNC_FACTOR ( MAXLAYERS )
      INTEGER, INTENT(IN) ::            N_PARTLAYERS
      INTEGER, INTENT(IN) ::            NFINELAYERS
      INTEGER, INTENT(IN) ::            PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT(IN) ::   PARTAU_VERT ( MAX_PARTLAYERS )
      LOGICAL, INTENT(IN) ::            PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT(IN) ::            PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )

      DOUBLE PRECISION, INTENT (OUT) ::  UP_MULTIPLIERS &
          ( MAXLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (OUT) ::  DN_MULTIPLIERS &
          ( MAXLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (OUT) ::  UP_MULTIPLIERS_UT &
          ( MAX_PARTLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (OUT) ::  DN_MULTIPLIERS_UT &
          ( MAX_PARTLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (OUT) ::  UP_LOSTRANS &
          ( MAXLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (OUT) ::  DN_LOSTRANS &
          ( MAXLAYERS, MAX_GEOMETRIES )
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
      LOGICAL, INTENT (OUT) ::           FAIL
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  GEOMETRY ROUTINE INPUTS AND OUTPUTS
!  -----------------------------------

!  CONTROL

      LOGICAL ::           DO_FINE
      LOGICAL ::           DO_PARTIALS
      DOUBLE PRECISION ::  ALPHA_BOA, THETA_BOA, PHI_BOA

!  MAIN OUTPUTS (GEOMETRY)

      INTEGER ::           NTRAVERSE(0:MAXLAYERS)
      DOUBLE PRECISION ::  SUNPATHS(0:MAXLAYERS,MAXLAYERS)
      DOUBLE PRECISION ::  RADII   (0:MAXLAYERS)
      DOUBLE PRECISION ::  ALPHA_ALL  (0:MAXLAYERS)

!  FINE LEVEL OUTPUT (GEOMETRY)

      INTEGER ::           NTRAVERSE_FINE(MAXLAYERS,MAXFINELAYERS)
      DOUBLE PRECISION ::  SUNPATHS_FINE (MAXLAYERS,MAXLAYERS,MAXFINELAYERS)
      DOUBLE PRECISION ::  RADII_FINE    (MAXLAYERS,MAXFINELAYERS)
      DOUBLE PRECISION ::  ALPHA_FINE    (MAXLAYERS,MAXFINELAYERS)

!  PARTIAL LAYER OUTPUT (GEOMETRY)

      INTEGER ::           NTRAVERSE_UT(MAX_PARTLAYERS)
      DOUBLE PRECISION ::  SUNPATHS_UT (MAX_PARTLAYERS,MAXLAYERS)
      DOUBLE PRECISION ::  RADII_UT    (MAX_PARTLAYERS)
      DOUBLE PRECISION ::  ALPHA_UT    (MAX_PARTLAYERS)

!  OTHER (INCIDENTAL) GEOMETRICAL OUTPUT

      DOUBLE PRECISION ::  LOSPATHS(MAXLAYERS)
      DOUBLE PRECISION ::  LOSPATHS_UT_UP(MAX_PARTLAYERS)
      DOUBLE PRECISION ::  LOSPATHS_UT_DN(MAX_PARTLAYERS)
      DOUBLE PRECISION ::  THETA_ALL  (0:MAXLAYERS)
      DOUBLE PRECISION ::  PHI_ALL    (0:MAXLAYERS)
      DOUBLE PRECISION ::  COSSCAT_UP (0:MAXLAYERS)
      DOUBLE PRECISION ::  COSSCAT_DN (0:MAXLAYERS)

!  EXTINCTION

      DOUBLE PRECISION ::  EXTINCTION ( MAXLAYERS)

!  PARTIAL LAYER HEIGHTS

      DOUBLE PRECISION ::  HEIGHT_GRID_UT(MAX_PARTLAYERS)

!  LOCAL VARIABLES
!  ---------------

!  CUMULATIVE TRANSMITTANCES
!    NEW, 27 JANUARY 2010. TOA CONTRIBUTION FUNCTIONS CODE.

      DOUBLE PRECISION ::  OGCTRANS(MAXLAYERS,MAX_GEOMETRIES)

!  INDICES

      INTEGER ::           N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER ::           UT, UTA, UM, IA, NC, IB, V, O1, O2, NM1, &
                           LUM, LIB, LUA

!  HELP VARIABLES (DOUBLE PRECISION)

      DOUBLE PRECISION ::  FINAL_SOURCE, HELP, SS_CUMSOURCE
      DOUBLE PRECISION ::  SS_LAYERSOURCE, SSCORRECTION, XT
      DOUBLE PRECISION ::  CTHETA, STHETA, CALPHA, SALPHA, CPHI, CSA, TRANS
      DOUBLE PRECISION ::  VSIGN, ZMAT_LOCAL(4,4), FMAT_LOCAL(6), DNM1
      CHARACTER (LEN=3) :: CV

      INTEGER ::           PARTLAYERS_LAYERFINEIDX ( MAX_PARTLAYERS )

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

!  Local user indices

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

!  CREATE EXTINCTIONS

      DO N = 1, NLAYERS
        HELP = HEIGHT_GRID(N-1) - HEIGHT_GRID(N)
        EXTINCTION(N) = DELTAU_VERT(N) / HELP
      ENDDO

!  CREATE THE HEIGHTS OF PARTIAL LAYERS

      IF ( DO_PARTIALS ) THEN
        DO UT = 1, N_PARTLAYERS
          N = PARTLAYERS_LAYERIDX(UT)
          XT = DELTAU_VERT(N) - PARTAU_VERT(UT)
          HEIGHT_GRID_UT(UT) = HEIGHT_GRID(N) + XT / EXTINCTION(N)
        ENDDO
      ENDIF

!  ADDITIONAL DELTA-M SCALING
!  --------------------------

!  NEW SECTION. R. SPURR, 07 SEPTEMBER 2007.
!   TMS GETS MODIFIED BY (1-F). SAVE THE TRUNCATION FACTOR.
!   PHASE FUNCTION MOMENTS ARE MODIFIED LATER ON.

      IF ( DO_SSCORR_TRUNCATION ) THEN
        NM1  = NGREEK_MOMENTS_INPUT
        DNM1 = DBLE(2*NM1+1)
        DO N = 1, NLAYERS
          SSFDEL(N) = GREEKMAT_TOTAL_INPUT(NM1,N,1) / DNM1
          TMS(N) = TMS(N) * ( ONE - SSFDEL(N) )
        ENDDO
      ENDIF

!  SOURCE TERM CALCULATIONS
!  ========================

!  Lattice calculation

      IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN

!  START THE MAIN LOOP OVER ALL SOLAR AND VIEWING GEOMETRIES
!   USE THE ADJUSTED VALUES OF THE ANGLES

        DO UM = 1, N_USER_STREAMS
          ALPHA_BOA = USER_VZANGLES_ADJUST(UM)
          DO IB = 1, NBEAMS
            DO IA = 1, N_USER_RELAZMS
              THETA_BOA = SZANGLES_ADJUST(UM,IB,IA)
              PHI_BOA   = USER_RELAZMS_ADJUST(UM,IB,IA)
              V = VZA_OFFSETS(IB,UM) + IA

!  UPWELLING SOURCE TERMS
!  ----------------------

              IF ( DO_UPWELLING ) THEN

!  CALL TO GEOMETRY ROUTINE, PATH DISTANCES ETC....

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
!              DO N = 0, NLAYERS
!                WRITE(45,'(I4,101F10.5)')N,(SUNPATHS(N,V),V=1,NLAYERS)
!              ENDDO
!              PAUSE

!  MULTIPLIERS, TRANSMITTANCES

                CALL OUTGOING_INTEGRATION_UP &
                 ( NLAYERS, NFINELAYERS, DO_PARTIALS, EXTINCTION, &
                   N_PARTLAYERS, PARTLAYERS_LAYERIDX, &
                   PARTLAYERS_LAYERFINEIDX, &
                   SUNPATHS, RADII, NTRAVERSE, ALPHA_ALL, &
                   SUNPATHS_UT, RADII_UT, NTRAVERSE_UT, ALPHA_UT, &
                   SUNPATHS_FINE, NTRAVERSE_FINE, ALPHA_FINE, &
                   UP_MULTIPLIERS(1,V), UP_LOSTRANS(1,V), BOA_ATTN(V), &
                   UP_MULTIPLIERS_UT(1,V), UP_LOSTRANS_UT(1,V) )

!        DO N = 1, NLAYERS
!          WRITE(85,'(I4,1P2E18.10)')N,UP_LOSTRANS(N,1),UP_MULTIPLIERS(N,1)
!        ENDDO

!  DEBUG
!                WRITE(45,*)V
!                DO N = 1, NLAYERS
!                  WRITE(45,'(1P2E18.10)')
!     &             UP_LOSTRANS(N,V),UP_MULTIPLIERS(N,V)
!                ENDDO
!                PAUSE

!                IF (V.EQ.1)WRITE(35,*)V
!                DO UT = 1, N_PARTLAYERS
!                  IF (V.EQ.1)WRITE(35,'(1P2E18.10)')
!     &             UP_LOSTRANS_UT(UT,V),UP_MULTIPLIERS_UT(UT,V)
!                ENDDO

!  DEBUG: MULTIPLIERS (ALL LAYERS), UPWELLING OLD WAY
!                CALL OUTGOING_INTEGRATION_OLD_UP
!     I             ( NLAYERS, EXTINCTION, DELTAU_VERT,
!     I               LOSPATHS, SUNPATHS, NTRAVERSE,
!     O               UP_MULTIPLIERS(1,V), UP_LOSTRANS(1,V) )
!                WRITE(46,*)V
!                DO N = 1, NLAYERS
!                  WRITE(46,*)UP_LOSTRANS(N,V),UP_MULTIPLIERS(N,V)
!                ENDDO

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

!  CALL TO SCATTERING LAW FOR PHASE MATRIX ZMAT

                  CALL SSCORR_OUTGOING_ZMATRIX &
                  ( DO_SSCORR_TRUNCATION, SUNLIGHT, N, NSTOKES, &
                    NGREEK_MOMENTS_INPUT, LAYER_MAXMOMENTS(N), &
                    GREEKMAT_TOTAL_INPUT, SSFDEL(N), &
                    CTHETA, STHETA, CALPHA, SALPHA, CPHI, &
                    PHI_ALL(N), CSA, VSIGN, &
                    ZMAT_LOCAL, FMAT_LOCAL )

!  PHASE MATRIX (MULTIPLIED BY TMS FACTOR). SAVE THEM.
!     SUNLIGHT ONLY, THE FIRST COLUMN OF THE MATRIX
!    GENERAL CASE, CODE PROGRAMMED 05 OCTOBER 2010

                  IF ( STERM_LAYERMASK_UP(N) ) THEN
                   DO O1 = 1, NSTOKES
                    DO O2 = 1, NSTOKES
                     ZMAT_UP(V,N,O1,O2) = ZMAT_LOCAL(O1,O2) * TMS(N)
                    ENDDO
                   ENDDO
                  ENDIF

!  DEBUG FOR 2OS GEOMTEST
!      WRITE(46,'(I4,1P3E20.10)') &
!              N,ZMAT_LOCAL(1,1),ZMAT_LOCAL(2,1),ZMAT_LOCAL(3,1)

!  FINISH THE LAYER LOOP AND END UPWELLING

                ENDDO
              ENDIF

!  DOWNWELLING SOURCE TERMS: PHASE MATRIX, MULTIPLIERS, TRANSMITTANCES
!  -------------------------------------------------------------------

              IF ( DO_DNWELLING ) THEN

!  CALL TO GEOMETRY ROUTINE, PATH DISTANCES ETC....

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

!  MULTIPLIERS AND TRANSMITTANCES

                CALL OUTGOING_INTEGRATION_DN &
                 ( NLAYERS, NFINELAYERS, &
                   DO_PARTIALS, EXTINCTION, &
                   N_PARTLAYERS, PARTLAYERS_LAYERIDX, &
                   PARTLAYERS_LAYERFINEIDX, &
                   SUNPATHS, RADII, NTRAVERSE, ALPHA_ALL, &
                   SUNPATHS_UT,   NTRAVERSE_UT,   ALPHA_UT, &
                   SUNPATHS_FINE, NTRAVERSE_FINE, ALPHA_FINE, &
                   DN_MULTIPLIERS(1,V), DN_LOSTRANS(1,V), &
                   DN_MULTIPLIERS_UT(1,V), DN_LOSTRANS_UT(1,V) )

!  DEBUG
!                WRITE(65,*)V
!                DO N = 1, NLAYERS
!                  WRITE(65,'(1P2E18.10)')
!     &             DN_LOSTRANS(N,V),DN_MULTIPLIERS(N,V)
!                ENDDO

!  DEBUG: MULTIPLIERS (ALL LAYERS), DOWN WELLING OLD WAY
!                CALL OUTGOING_INTEGRATION_OLD_DN
!     I             ( NLAYERS, EXTINCTION, DELTAU_VERT,
!     I               LOSPATHS, SUNPATHS, NTRAVERSE,
!     O               DN_MULTIPLIERS(1,V), DN_LOSTRANS(1,V) )
!                WRITE(66,*)V
!                DO N = 1, NLAYERS
!                  WRITE(66,*)DN_LOSTRANS(N,V),DN_MULTIPLIERS(N,V)
!                ENDDO

!  START LAYER LOOP

                DO N = 1, NLAYERS

!  TRIGONOMETRY

                  CTHETA = DCOS(THETA_ALL(N-1))
                  STHETA = DSIN(THETA_ALL(N-1))
                  CALPHA = DCOS(ALPHA_ALL(N-1))
                  SALPHA = DSIN(ALPHA_ALL(N-1))
                  CPHI   = DCOS(PHI_ALL(N-1))
                  CSA    = COSSCAT_DN(N-1)
                  VSIGN  = +1.0D0

!  CALL TO SCATTERING LAW FOR PHASE MATRIX ZMAT

                  CALL SSCORR_OUTGOING_ZMATRIX &
                  ( DO_SSCORR_TRUNCATION, SUNLIGHT, N, NSTOKES, &
                    NGREEK_MOMENTS_INPUT, LAYER_MAXMOMENTS(N), &
                    GREEKMAT_TOTAL_INPUT, SSFDEL(N), &
                    CTHETA, STHETA, CALPHA, SALPHA, CPHI, &
                    PHI_ALL(N-1), CSA, VSIGN, &
                    ZMAT_LOCAL, FMAT_LOCAL )

!  PHASE MATRIX (MULTIPLIED BY TMS FACTOR). SAVE THEM.
!     SUNLIGHT ONLY, THE FIRST COLUMN OF THE MATRIX
!    GENERAL CASE, CODE PROGRAMMED 05 OCTOBER 2010

                  IF ( STERM_LAYERMASK_DN(N) ) THEN
                   DO O1 = 1, NSTOKES
                    DO O2 = 1, NSTOKES
                     ZMAT_DN(V,N,O1,O2) = ZMAT_LOCAL(O1,O2) * TMS(N)
                    ENDDO
                   ENDDO
                  ENDIF

!  END LAYER LOOP AND DOWNWELLING

                ENDDO
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
!   USE THE ADJUSTED VALUES OF THE ANGLES

        DO V = 1, N_GEOMETRIES
          ALPHA_BOA = USER_VZANGLES_ADJUST(V)
          THETA_BOA = SZANGLES_ADJUST(LUM,V,LUA)
          PHI_BOA   = USER_RELAZMS_ADJUST(LUM,LIB,V)

!  UPWELLING SOURCE TERMS
!  ----------------------

          IF ( DO_UPWELLING ) THEN

!  CALL TO GEOMETRY ROUTINE, PATH DISTANCES ETC....

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
!           DO N = 0, NLAYERS
!             WRITE(45,'(I4,101F10.5)')N,(SUNPATHS(N,V),V=1,NLAYERS)
!           ENDDO
!           PAUSE

!  MULTIPLIERS, TRANSMITTANCES

            CALL OUTGOING_INTEGRATION_UP &
               ( NLAYERS, NFINELAYERS, DO_PARTIALS, EXTINCTION, &
                 N_PARTLAYERS, PARTLAYERS_LAYERIDX, &
                 PARTLAYERS_LAYERFINEIDX, &
                 SUNPATHS, RADII, NTRAVERSE, ALPHA_ALL, &
                 SUNPATHS_UT, RADII_UT, NTRAVERSE_UT, ALPHA_UT, &
                 SUNPATHS_FINE, NTRAVERSE_FINE, ALPHA_FINE, &
                 UP_MULTIPLIERS(1,V), UP_LOSTRANS(1,V), BOA_ATTN(V), &
                 UP_MULTIPLIERS_UT(1,V), UP_LOSTRANS_UT(1,V) )

!           DO N = 1, NLAYERS
!             WRITE(85,'(I4,1P2E18.10)')N,UP_LOSTRANS(N,1),UP_MULTIPLIERS(N,1)
!           ENDDO

!  DEBUG
!           WRITE(45,*)V
!           DO N = 1, NLAYERS
!             WRITE(45,'(1P2E18.10)')
!     &        UP_LOSTRANS(N,V),UP_MULTIPLIERS(N,V)
!           ENDDO
!           PAUSE

!           IF (V.EQ.1)WRITE(35,*)V
!           DO UT = 1, N_PARTLAYERS
!             IF (V.EQ.1)WRITE(35,'(1P2E18.10)')
!     &         UP_LOSTRANS_UT(UT,V),UP_MULTIPLIERS_UT(UT,V)
!           ENDDO

!  DEBUG: MULTIPLIERS (ALL LAYERS), UPWELLING OLD WAY
!           CALL OUTGOING_INTEGRATION_OLD_UP
!     I        ( NLAYERS, EXTINCTION, DELTAU_VERT,
!     I          LOSPATHS, SUNPATHS, NTRAVERSE,
!     O          UP_MULTIPLIERS(1,V), UP_LOSTRANS(1,V) )
!           WRITE(46,*)V
!           DO N = 1, NLAYERS
!             WRITE(46,*)UP_LOSTRANS(N,V),UP_MULTIPLIERS(N,V)
!           ENDDO

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

!  CALL TO SCATTERING LAW FOR PHASE MATRIX ZMAT

              CALL SSCORR_OUTGOING_ZMATRIX &
                ( DO_SSCORR_TRUNCATION, SUNLIGHT, N, NSTOKES, &
                  NGREEK_MOMENTS_INPUT, LAYER_MAXMOMENTS(N), &
                  GREEKMAT_TOTAL_INPUT, SSFDEL(N), &
                  CTHETA, STHETA, CALPHA, SALPHA, CPHI, &
                  PHI_ALL(N), CSA, VSIGN, &
                  ZMAT_LOCAL, FMAT_LOCAL )

!  PHASE MATRIX (MULTIPLIED BY TMS FACTOR). SAVE THEM.
!     SUNLIGHT ONLY, THE FIRST COLUMN OF THE MATRIX
!    GENERAL CASE, CODE PROGRAMMED 05 OCTOBER 2010

              IF ( STERM_LAYERMASK_UP(N) ) THEN
                DO O1 = 1, NSTOKES
                  DO O2 = 1, NSTOKES
                    ZMAT_UP(V,N,O1,O2) = ZMAT_LOCAL(O1,O2) * TMS(N)
                  ENDDO
                ENDDO
              ENDIF

!  DEBUG FOR 2OS GEOMTEST
!      WRITE(46,'(I4,1P3E20.10)') &
!              N,ZMAT_LOCAL(1,1),ZMAT_LOCAL(2,1),ZMAT_LOCAL(3,1)

!  FINISH THE LAYER LOOP AND END UPWELLING

            ENDDO
          ENDIF

!  DOWNWELLING SOURCE TERMS: PHASE MATRIX, MULTIPLIERS, TRANSMITTANCES
!  -------------------------------------------------------------------

          IF ( DO_DNWELLING ) THEN

!  CALL TO GEOMETRY ROUTINE, PATH DISTANCES ETC....

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

!  MULTIPLIERS AND TRANSMITTANCES

            CALL OUTGOING_INTEGRATION_DN &
               ( NLAYERS, NFINELAYERS, &
                 DO_PARTIALS, EXTINCTION, &
                 N_PARTLAYERS, PARTLAYERS_LAYERIDX, &
                 PARTLAYERS_LAYERFINEIDX, &
                 SUNPATHS, RADII, NTRAVERSE, ALPHA_ALL, &
                 SUNPATHS_UT,   NTRAVERSE_UT,   ALPHA_UT, &
                 SUNPATHS_FINE, NTRAVERSE_FINE, ALPHA_FINE, &
                 DN_MULTIPLIERS(1,V), DN_LOSTRANS(1,V), &
                 DN_MULTIPLIERS_UT(1,V), DN_LOSTRANS_UT(1,V) )

!  DEBUG
!           WRITE(65,*)V
!           DO N = 1, NLAYERS
!             WRITE(65,'(1P2E18.10)')
!     &        DN_LOSTRANS(N,V),DN_MULTIPLIERS(N,V)
!           ENDDO

!  DEBUG: MULTIPLIERS (ALL LAYERS), DOWN WELLING OLD WAY
!           CALL OUTGOING_INTEGRATION_OLD_DN
!     I        ( NLAYERS, EXTINCTION, DELTAU_VERT,
!     I          LOSPATHS, SUNPATHS, NTRAVERSE,
!     O          DN_MULTIPLIERS(1,V), DN_LOSTRANS(1,V) )
!           WRITE(66,*)V
!           DO N = 1, NLAYERS
!             WRITE(66,*)DN_LOSTRANS(N,V),DN_MULTIPLIERS(N,V)
!           ENDDO

!  START LAYER LOOP

            DO N = 1, NLAYERS

!  TRIGONOMETRY

              CTHETA = DCOS(THETA_ALL(N-1))
              STHETA = DSIN(THETA_ALL(N-1))
              CALPHA = DCOS(ALPHA_ALL(N-1))
              SALPHA = DSIN(ALPHA_ALL(N-1))
              CPHI   = DCOS(PHI_ALL(N-1))
              CSA    = COSSCAT_DN(N-1)
              VSIGN  = +1.0D0

!  CALL TO SCATTERING LAW FOR PHASE MATRIX ZMAT

              CALL SSCORR_OUTGOING_ZMATRIX &
                ( DO_SSCORR_TRUNCATION, SUNLIGHT, N, NSTOKES, &
                  NGREEK_MOMENTS_INPUT, LAYER_MAXMOMENTS(N), &
                  GREEKMAT_TOTAL_INPUT, SSFDEL(N), &
                  CTHETA, STHETA, CALPHA, SALPHA, CPHI, &
                  PHI_ALL(N-1), CSA, VSIGN, &
                  ZMAT_LOCAL, FMAT_LOCAL )

!  PHASE MATRIX (MULTIPLIED BY TMS FACTOR). SAVE THEM.
!     SUNLIGHT ONLY, THE FIRST COLUMN OF THE MATRIX
!    GENERAL CASE, CODE PROGRAMMED 05 OCTOBER 2010

              IF ( STERM_LAYERMASK_DN(N) ) THEN
                DO O1 = 1, NSTOKES
                  DO O2 = 1, NSTOKES
                    ZMAT_DN(V,N,O1,O2) = ZMAT_LOCAL(O1,O2) * TMS(N)
                  ENDDO
                ENDDO
              ENDIF

!  END LAYER LOOP AND DOWNWELLING

            ENDDO
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
!     SET FINAL CUMULATIVE SOURCE AND SINGLE SCATTER STOKES VECTOR
!  ONGRID OUTPUT--------
!     SET FINAL CUMULATIVE SOURCE AND SINGLE SCATTER STOKES VECTOR
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
!      PROGRAMMED FOR GENERAL CASE (INCLUDING SUNLIGHT)

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

!  OFFGRID OUTPUT :
!    ADD ADDITIONAL PARTIAL LAYER SOURCE TERM = EXACT Z-MATRIX * MULTIPL
!     SET FINAL CUMULATIVE SOURCE AND CORRECT STOKES VECTOR
!  ONGRID OUTPUT :
!     SET FINAL CUMULATIVE SOURCE AND CORRECT STOKES VECTOR
!      PROGRAMMED FOR GENERAL CASE (INCLUDING SUNLIGHT)

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

!  FINISH

      RETURN
      END SUBROUTINE VLIDORT_SSCORR_OUTGOING

!

      SUBROUTINE SSCORR_OUTGOING_ZMATRIX &
        ( DO_SSCORR_TRUNCATION, DO_SUNLIGHT, LAYER, NSTOKES, &
          NGREEKMOMS,  LAYER_MAXMOMENTS, GREEKMATRIX, TRUNCFACTOR, &
          CTHETA, STHETA, CALPHA, SALPHA, CPHI, PHI, COSSCAT, VSIGN, &
          ZMAT, FMAT )

!  INCLUDE FILE
!  ------------

      USE VLIDORT_PARS

      IMPLICIT NONE


!  INPUT
!  -----

!  CONTROL FOR USING NATURAL LIGHT CASE

      LOGICAL, INTENT(IN) ::           DO_SUNLIGHT

!  CONTROL FOR USING THE LOCAL TRUNCATION

      LOGICAL, INTENT(IN) ::           DO_SSCORR_TRUNCATION

!  CONTROL INTEGERS

      INTEGER, INTENT(IN) ::           LAYER, NSTOKES, NGREEKMOMS

!  SCATTERING INPUT INFORMATION

      INTEGER, INTENT(IN) ::           LAYER_MAXMOMENTS
      DOUBLE PRECISION, INTENT(IN) ::  GREEKMATRIX &
          (0:MAXMOMENTS_INPUT,MAXLAYERS,16)
      DOUBLE PRECISION, INTENT(IN) ::  TRUNCFACTOR

!  ZENITH ANGLE, AZIMUTH ANGLE, SCATTER ANGLE (SINES/COSINES)

      DOUBLE PRECISION, INTENT(IN) ::  CTHETA
      DOUBLE PRECISION, INTENT(IN) ::  STHETA
      DOUBLE PRECISION, INTENT(IN) ::  CALPHA
      DOUBLE PRECISION, INTENT(IN) ::  SALPHA
      DOUBLE PRECISION, INTENT(IN) ::  CPHI
      DOUBLE PRECISION, INTENT(IN) ::  COSSCAT
      DOUBLE PRECISION, INTENT(IN) ::  VSIGN

!  AZIMUTH ANGLE. ADDED FOR VERSION 2.4R, REQUIRED FOR PHI > 180.
!    R. SPURR AND V. NATRAJ, 01 MAY 2009
!    PROPER USAGE, 01 OCTOBER 2010. R. SPURR

      DOUBLE PRECISION, INTENT(IN) ::  PHI

!  OUTPUT
!  ------

!  Z-MATRIX, FIRST COLUMN ONLY. F-MATRIX ( SUNLIGHT CASE )
!  GENERALIZED, 05 OCTOBER 2010, NON-TESTED

      DOUBLE PRECISION, INTENT(OUT) ::   ZMAT ( 4, 4 )
      DOUBLE PRECISION, INTENT(OUT) ::   FMAT ( 6 )

!  LOCAL
!  -----

      DOUBLE PRECISION ::   C1, S1, C2, S2, SINSCAT, HELP_SINSCAT
      DOUBLE PRECISION ::   CSIG1, CSIG2, SSIG1, SSIG2, CSIG1_2, CSIG2_2

      DOUBLE PRECISION ::   DL, QROOT6, UUU, PHIDEG, DNL1, FDNL1, FACT
      DOUBLE PRECISION ::   FAC1, FAC2, SQL4, SQL41
      DOUBLE PRECISION ::   TMP1, TMP2, SUM23, DIF23
      DOUBLE PRECISION ::   GK11, GK12, GK34, GK44, GK22, GK33
      DOUBLE PRECISION ::   FL2, FLL1, PERLL4, Q, WFACT, DL1
      DOUBLE PRECISION ::   HELP2C1, HELP2S1, HELP3C1, HELP3S1

      INTEGER ::            GREEKMAT_INDEX(6), K, L, LNEW, LOLD, ITMP, N
      INTEGER ::            INDEX_11, INDEX_12, INDEX_34
      INTEGER ::            INDEX_22, INDEX_33, INDEX_44

!  GENERALIZED SPHERICAL FUNCTIONS
!    P2P2 AND P2M2 ADDED, 05 OCTOBER 2010

      DOUBLE PRECISION ::   P00(2), P02(2), P2P2(2), P2M2(2)

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

!  DEBUG, 06 OCTOBER 2010. FOUND DOWNWELLING BUG, NEED SEPARATE GEOMETRY
!      WRITE(7,'(I4,2F6.1,1P6E15.7)')N,VSIGN,PHI/DEG_TO_RAD,
!     &      SALPHA, CALPHA, STHETA, CTHETA, CSIG1,CSIG2

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
!  SEE H/VDM, EQS. 94-95.

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

!  F-MATRICES
!  ----------

      QROOT6 = -0.25D0 * DSQRT(6.0D0)

! INITIALISE F-MATRIX

      DO K = 1, 6
        FMAT(K) = ZERO
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

        IF ( DO_SSCORR_TRUNCATION ) THEN
          DNL1  = DBLE(2*L + 1 )
          FDNL1 = TRUNCFACTOR * DNL1
          FACT  = ONE - TRUNCFACTOR
          GK11 = ( GREEKMATRIX(L,N,GREEKMAT_INDEX(1)) - FDNL1 ) / FACT
!mick fix 1/21/2013 - added IF structure and ELSE section
          IF ( NSTOKES .GT. 1 ) THEN
            GK12 =   GREEKMATRIX(L,N,GREEKMAT_INDEX(3)) / FACT
            IF ( .NOT. DO_SUNLIGHT ) THEN
             GK44 = ( GREEKMATRIX(L,N,GREEKMAT_INDEX(6)) - FDNL1 ) / FACT
             GK34 =   GREEKMATRIX(L,N,GREEKMAT_INDEX(5)) / FACT
             GK22 = ( GREEKMATRIX(L,N,GREEKMAT_INDEX(2)) - FDNL1 ) / FACT
             GK33 = ( GREEKMATRIX(L,N,GREEKMAT_INDEX(4)) - FDNL1 ) / FACT
            ENDIF
          ELSE
             GK12 = ZERO
             IF ( .NOT. DO_SUNLIGHT ) THEN
              GK44 = ZERO
              GK34 = ZERO
              GK22 = ZERO
              GK33 = ZERO
             ENDIF
          ENDIF
        ELSE
          GK11 = GREEKMATRIX(L,N,GREEKMAT_INDEX(1))
!mick fix 1/21/2013 - added IF structure and ELSE section
          IF ( NSTOKES .GT. 1 ) THEN
            GK12 = GREEKMATRIX(L,N,GREEKMAT_INDEX(3))
            IF ( .NOT. DO_SUNLIGHT ) THEN
             GK44 = GREEKMATRIX(L,N,GREEKMAT_INDEX(6))
             GK34 = GREEKMATRIX(L,N,GREEKMAT_INDEX(5))
             GK22 = GREEKMATRIX(L,N,GREEKMAT_INDEX(2))
             GK33 = GREEKMATRIX(L,N,GREEKMAT_INDEX(4))
            ENDIF
          ELSE
             GK12 = ZERO
             IF ( .NOT. DO_SUNLIGHT ) THEN
              GK44 = ZERO
              GK34 = ZERO
              GK22 = ZERO
              GK33 = ZERO
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
          Q     = DL  * ( DL1*DL1 - FOUR)
          WFACT = FL2 * ( FLL1 * UUU - FOUR )
          P2P2(LOLD) = (WFACT*P2P2(LNEW) - Q*P2P2(LOLD)) * PERLL4
          WFACT = FL2 * ( FLL1 * UUU + FOUR )
          P2M2(LOLD) = (WFACT*P2M2(LNEW) - Q*P2M2(LOLD)) * PERLL4

        END IF

! SWITCH INDICES SO THAT LNEW INDICATES THE FUNCTION WITH
! THE PRESENT INDEX VALUE L, THIS MECHANISM PREVENTS SWAPPING
! OF ENTIRE ARRAYS.

        ITMP = LNEW
        LNEW = LOLD
        LOLD = ITMP

! NOW ADD THE L-TH TERM TO THE SCATTERING MATRIX.
! SEE DE HAAN ET AL. (1987) EQS. (68)-(73).

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

!  END MOMENTS

      END DO

!   THIS MUST BE DONE AFTER THE MOMENT LOOP.

      IF ( .NOT. DO_SUNLIGHT ) THEN
        FMAT(INDEX_22) = HALF * ( FMAT(INDEX_22) + FMAT(INDEX_33) )
        FMAT(INDEX_33) = FMAT(INDEX_22) - FMAT(INDEX_33)
      ENDIF

! REMEMBER FOR MIE SCATTERING : F11 = F22 AND F33 = F44
!  THIS CODE IS NO LONGER REQUIRED, AS WE HAVE INTRODUCED CODE NOW
!   FOR RANDOMLY ORIENTED SPHEROIDS. THE SYMMETRY SHOULD STILL OF
!   COURSE BE PRESENT FOR THE MIE PARTICLES, SO THIS WILL BE A
!   CHECK ON THE NEW CODE. R. SPURR AND V. NATRAJ,, 20 MARCH 2006

!      IF ( .NOT. DO_SUNLIGHT ) THEN
!        FMAT(INDEX_22) = FMAT(INDEX_11)
!        FMAT(INDEX_33) = FMAT(INDEX_44)
!      ENDIF

!  Z-MATRIX CALCULATION
!  ====================

!  INITIALIZE Z-MATRIX

      ZMAT = ZERO

!  SCALAR CASE

      IF ( NSTOKES .EQ. 1 ) THEN
        ZMAT(1,1) =   FMAT(INDEX_11)
        RETURN
      ENDIF

!  NATURAL LIGHT, FIRST COLUMN ONLY

      IF ( DO_SUNLIGHT ) THEN
        ZMAT(1,1) =   FMAT(INDEX_11)
        ZMAT(2,1) =  -FMAT(INDEX_12) * C2
        ZMAT(3,1) =   FMAT(INDEX_12) * S2
        ZMAT(4,1) =   ZERO
      ENDIF

!  OTHER SOURCES, USE FULL 4X4 MATRIX
!    CODE INTRODUCED BUT NOT TESTED, 05 OCTOBER 2010

      IF ( .NOT. DO_SUNLIGHT ) THEN
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

!  FINISH

      RETURN
      END SUBROUTINE SSCORR_OUTGOING_ZMATRIX

!

      SUBROUTINE OUTGOING_INTEGRATION_UP &
        ( NLAYERS, NFINELAYERS, DO_PARTIALS, EXTINCTION, &
          N_PARTIALS, PARTIALS_IDX, PARTIALS_FINEIDX, &
          SUNPATHS, RADII,  NTRAVERSE, ALPHA_ALL, &
          SUNPATHS_P, RADII_P, NTRAVERSE_P,    ALPHA_P, &
          SUNPATHS_FINE, NTRAVERSE_FINE, ALPHA_FINE, &
          MULTIPLIERS, LOSTRANS, BOA_ATTN, &
          MULTIPLIERS_P, LOSTRANS_P )

!  DOES THE OPTICAL DEPTH INTEGRATION OVER LAYERS.
!  PARTIAL LAYER INTEGRATION ADDED SEPTEMBER 2007.

!  INCLUDE FILES

      USE VLIDORT_PARS

      IMPLICIT NONE


!  GEOMETRY ROUTINE INPUTS
!  -----------------------

!  CONTROL

      LOGICAL, INTENT(IN) ::           DO_PARTIALS
      INTEGER, INTENT(IN) ::           NFINELAYERS, NLAYERS
      INTEGER, INTENT(IN) ::           N_PARTIALS
      INTEGER, INTENT(IN) ::           PARTIALS_IDX    (MAX_PARTLAYERS)
      INTEGER, INTENT(IN) ::           PARTIALS_FINEIDX(MAX_PARTLAYERS)

!  WHOLE LAYERS

      INTEGER, INTENT(IN) ::           NTRAVERSE(0:MAXLAYERS)
      DOUBLE PRECISION, INTENT(IN) ::  SUNPATHS(0:MAXLAYERS,MAXLAYERS)
      DOUBLE PRECISION, INTENT(IN) ::  RADII   (0:MAXLAYERS)
      DOUBLE PRECISION, INTENT(IN) ::  ALPHA_ALL  (0:MAXLAYERS)

!  FINE LEVEL

      INTEGER, INTENT(IN) ::           NTRAVERSE_FINE(MAXLAYERS,MAXFINELAYERS)
      DOUBLE PRECISION, INTENT(IN) ::  SUNPATHS_FINE &
          (MAXLAYERS,MAXLAYERS,MAXFINELAYERS)
      DOUBLE PRECISION, INTENT(IN) ::  ALPHA_FINE    (MAXLAYERS,MAXFINELAYERS)

!  PARTIAL LAYERS

      INTEGER, INTENT(IN) ::           NTRAVERSE_P(MAX_PARTLAYERS)
      DOUBLE PRECISION, INTENT(IN) ::  SUNPATHS_P (MAX_PARTLAYERS,MAXLAYERS)
      DOUBLE PRECISION, INTENT(IN) ::  ALPHA_P    (MAX_PARTLAYERS)
      DOUBLE PRECISION, INTENT(IN) ::  RADII_P    (MAX_PARTLAYERS)

!  EXTINCTION

      DOUBLE PRECISION, INTENT(IN) ::  EXTINCTION (MAXLAYERS)

!  OUTPUTS
!  -------

      DOUBLE PRECISION, INTENT(OUT) ::  MULTIPLIERS (MAXLAYERS)
      DOUBLE PRECISION, INTENT(OUT) ::  LOSTRANS    (MAXLAYERS)
      DOUBLE PRECISION, INTENT(OUT) ::  MULTIPLIERS_P (MAX_PARTLAYERS)
      DOUBLE PRECISION, INTENT(OUT) ::  LOSTRANS_P    (MAX_PARTLAYERS)
      DOUBLE PRECISION, INTENT(OUT) ::  BOA_ATTN

!  LOCAL ARRAYS
!  ------------

!  LOCAL GEOEMETRY ARRAYS

      DOUBLE PRECISION ::  CSQ_FINE ( MAXFINELAYERS)
      DOUBLE PRECISION ::  COT_FINE ( MAXFINELAYERS)

!  LOCAL ATTENUATION FACTORS

      DOUBLE PRECISION ::  ATTN      ( 0:MAXLAYERS )
      DOUBLE PRECISION ::  ATTN_FINE ( MAXLAYERS, MAXFINELAYERS )
      DOUBLE PRECISION ::  ATTN_P    ( MAX_PARTLAYERS )

!  HELP VARIABLES
!  --------------

      INTEGER ::           N, J, K, UT, NP, NFINE_P
      DOUBLE PRECISION ::  TAU, SUM, SALPHA, CALPHA, DFINE1, RAYCON, KN
      DOUBLE PRECISION ::  CSQ_1, CSQ_2, COT_1, COT_2
      DOUBLE PRECISION ::  FUNC_1, FUNC_2, TRAN_1, TRAN, STEP, FUNC
      DOUBLE PRECISION ::  TERM1, TERM2, DFINE_P, STEP_F, STEP_P

      LOGICAL ::           DO_SPECIAL
      DOUBLE PRECISION ::  HN, HJ, HP

!  LOCAL OPTICAL THICKNESS CUTOFF
!      (SHOULD BE SAME AS MAX_TAU_SPATH IN VLIDORT)

      DOUBLE PRECISION, PARAMETER ::  LOCAL_CUTOFF = 32.0D0

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

!  ATTENUATION FUNCTIONS, WHOLE LAYERS

      DO N = 0, NLAYERS
       TAU     = 0.0D0
       ATTN(N) = 0.0D0
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
          ATTN_P(UT) = 0.0D0
          DO K = 1, NTRAVERSE_P(UT)
            TAU = TAU + SUNPATHS_P(UT,K) * EXTINCTION(K)
          ENDDO
          IF ( TAU.LE.LOCAL_CUTOFF) ATTN_P(UT) = DEXP(-TAU)
        ENDDO
      ENDIF

!  FINE GRID ATTENUATIONS

      DO N = 1, NLAYERS
       DO J = 1, NFINELAYERS
        TAU            = 0.0D0
        ATTN_FINE(N,J) = 0.0D0
        DO K = 1, NTRAVERSE_FINE(N,J)
          TAU = TAU + SUNPATHS_FINE(N,K,J) * EXTINCTION(K)
        ENDDO
        IF ( TAU .LE. LOCAL_CUTOFF ) ATTN_FINE(N,J) = DEXP(-TAU)
       ENDDO
      ENDDO

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
            HJ = HN - STEP * DBLE(J)
            TRAN = DEXP ( - KN * HJ )
            FUNC = ATTN_FINE(N,J) * TRAN
            SUM = SUM + FUNC
          ENDDO
          LOSTRANS(N)    = TRAN_1
          MULTIPLIERS(N) = SUM * STEP * KN

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
              HJ = HP - STEP_F * DBLE(J)
              TRAN = DEXP ( - KN * HJ )
              FUNC = ATTN_FINE(NP,J) * TRAN
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

!        WRITE(*,*)'SPECIAL PARTIAL UP',UT,NP,NFINE_P, &
!                    MULTIPLIERS_P(UT),LOSTRANS_P(UT)

!  END PARTIALS

        ENDDO

!  RETURN AND END SPECIALS CLAUSE

        RETURN

      ENDIF

!  WORK UP FROM THE BOTTOM OF THE ATMOSPHERE
!  =========================================

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

!   SURELY THIS IS WRONG.................???? 27 SEPTEMBER 2007
!        FUNC_1 = ATTN(N-1)   * CSQ_1
!        TRAN_2 = DEXP ( - KN * ( COT_2 - COT_1 ) )
!        FUNC_2 = ATTN(N) * CSQ_2 * TRAN_2

        FUNC_2 = ATTN(N-1)   * CSQ_2
        TRAN_1 = DEXP ( - KN * ( COT_2 - COT_1 ) )
        FUNC_1 = ATTN(N) * CSQ_1 * TRAN_1
        SUM = 0.5D0 * ( FUNC_1 + FUNC_2 )
        DO J = 1, NFINELAYERS
          TRAN = DEXP ( -KN * ( COT_2 - COT_FINE(J) ) )
          FUNC = ATTN_FINE(N,J) * TRAN * CSQ_FINE(J)
          SUM = SUM + FUNC
        ENDDO
        LOSTRANS(N)    = TRAN_1
        MULTIPLIERS(N) = SUM * STEP * KN

!        IF (N.EQ.6)WRITE(19,*)'WHOLE',N,MULTIPLIERS(N),LOSTRANS(N)

!  UPDATE GEOMETRY

        CSQ_1 = CSQ_2
        COT_1 = COT_2

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

!  INITIALIZE FOR LAYER OF OCCURRENCE

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
        TRAN_1 = DEXP ( - KN * ( COT_2 - COT_1 ) )
        FUNC_1 = ATTN(NP) * CSQ_1 * TRAN_1

        IF ( NFINE_P.GT.0) THEN
          SUM = 0.5D0 * FUNC_1
          DO J = 1, NFINE_P
            TRAN = DEXP ( -KN * ( COT_2 - COT_FINE(J) ) )
            FUNC = ATTN_FINE(NP,J) * TRAN * CSQ_FINE(J)
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

!        WRITE(*,*)'GENERAL PARTIAL UP',UT,NP,NFINE_P, &
!                    MULTIPLIERS_P(UT),LOSTRANS_P(UT)

!  FINISH PARTIALS LOOP

      ENDDO

!  FINISH

      RETURN
      END SUBROUTINE OUTGOING_INTEGRATION_UP

!

      SUBROUTINE OUTGOING_INTEGRATION_DN &
        ( NLAYERS, NFINELAYERS, DO_PARTIALS, EXTINCTION, &
          N_PARTIALS, PARTIALS_IDX, PARTIALS_FINEIDX, &
          SUNPATHS, RADII, NTRAVERSE, ALPHA_ALL, &
          SUNPATHS_P,   NTRAVERSE_P,   ALPHA_P, &
          SUNPATHS_FINE, NTRAVERSE_FINE, ALPHA_FINE, &
          MULTIPLIERS, LOSTRANS, &
          MULTIPLIERS_P, LOSTRANS_P )

!  DOES THE OPTICAL DEPTH INTEGRATION OVER LAYERS.
!  PARTIAL LAYER INTEGRATION ADDED SEPTEMBER 2007.

!  INCLUDE FILES

      USE VLIDORT_PARS

      IMPLICIT NONE


!  GEOMETRY ROUTINE INPUTS
!  -----------------------

!  CONTROL

      LOGICAL, INTENT(IN) ::           DO_PARTIALS
      INTEGER, INTENT(IN) ::           NFINELAYERS, NLAYERS
      INTEGER, INTENT(IN) ::           N_PARTIALS
      INTEGER, INTENT(IN) ::           PARTIALS_IDX    (MAX_PARTLAYERS)
      INTEGER, INTENT(IN) ::           PARTIALS_FINEIDX(MAX_PARTLAYERS)

!  WHOLE LAYERS

      INTEGER, INTENT(IN) ::           NTRAVERSE(0:MAXLAYERS)
      DOUBLE PRECISION, INTENT(IN) ::  SUNPATHS(0:MAXLAYERS,MAXLAYERS)
      DOUBLE PRECISION, INTENT(IN) ::  RADII   (0:MAXLAYERS)
      DOUBLE PRECISION, INTENT(IN) ::  ALPHA_ALL  (0:MAXLAYERS)

!  FINE LEVEL

      INTEGER, INTENT(IN) ::           NTRAVERSE_FINE(MAXLAYERS,MAXFINELAYERS)
      DOUBLE PRECISION, INTENT(IN) ::  SUNPATHS_FINE &
          (MAXLAYERS,MAXLAYERS,MAXFINELAYERS)
      DOUBLE PRECISION, INTENT(IN) ::  ALPHA_FINE    (MAXLAYERS,MAXFINELAYERS)

!  PARTIAL LAYERS

      INTEGER, INTENT(IN) ::           NTRAVERSE_P(MAX_PARTLAYERS)
      DOUBLE PRECISION, INTENT(IN) ::  SUNPATHS_P (MAX_PARTLAYERS,MAXLAYERS)
      DOUBLE PRECISION, INTENT(IN) ::  ALPHA_P    (MAX_PARTLAYERS)

!  EXTINCTION

      DOUBLE PRECISION, INTENT(IN) ::  EXTINCTION (MAXLAYERS)

!  OUTPUTS
!  -------

      DOUBLE PRECISION, INTENT(OUT) ::  MULTIPLIERS (MAXLAYERS)
      DOUBLE PRECISION, INTENT(OUT) ::  LOSTRANS    (MAXLAYERS)
      DOUBLE PRECISION, INTENT(OUT) ::  MULTIPLIERS_P (MAX_PARTLAYERS)
      DOUBLE PRECISION, INTENT(OUT) ::  LOSTRANS_P    (MAX_PARTLAYERS)

!  LOCAL ARRAYS
!  ------------

!  LOCAL GEOEMETRY ARRAYS

      DOUBLE PRECISION ::  CSQ_FINE ( MAXFINELAYERS)
      DOUBLE PRECISION ::  COT_FINE ( MAXFINELAYERS)

!  LOCAL ATTENUATION FACTORS

      DOUBLE PRECISION ::  ATTN      ( 0:MAXLAYERS )
      DOUBLE PRECISION ::  ATTN_FINE ( MAXLAYERS, MAXFINELAYERS )
      DOUBLE PRECISION ::  ATTN_P    ( MAX_PARTLAYERS )

!  HELP VARIABLES
!  --------------

      INTEGER ::           N, J, K, UT, NP, NFINE_P
      DOUBLE PRECISION ::  TAU, SUM, SALPHA, CALPHA, DFINE1, RAYCON, KN
      DOUBLE PRECISION ::  CSQ_1, CSQ_2, COT_1, COT_2
      DOUBLE PRECISION ::  FUNC_1, FUNC_2, TRAN_1, TRAN, STEP, FUNC
      DOUBLE PRECISION ::  TERM1, TERM2, DFINE_P, STEP_F, STEP_P

!  LOCAL OPTICAL THICKNESS CUTOFF
!      (SHOULD BE SAME AS MAX_TAU_SPATH IN VLIDORT)

      DOUBLE PRECISION, PARAMETER ::  LOCAL_CUTOFF = 32.0D0

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

!  CREATE SLANT PATH OPTICAL ATTENUATIONS, SOLAR BEAMS
!  ---------------------------------------------------

!  RAY CONSTANT AT TOA

      SALPHA = DSIN(ALPHA_ALL(0))
      RAYCON  = RADII(0) * SALPHA
      DFINE1 = DBLE(NFINELAYERS+1)

!  ATTENUATION FUNCTIONS, WHOLE LAYERS

      DO N = 0, NLAYERS
       TAU     = 0.0D0
       ATTN(N) = 0.0D0
       DO K = 1, NTRAVERSE(N)
         TAU = TAU + SUNPATHS(N,K) * EXTINCTION(K)
       ENDDO
       IF ( TAU .LE. LOCAL_CUTOFF ) ATTN(N) = DEXP(-TAU)
      ENDDO

!  ATTENUATIONS FOR PARTIALS

      IF ( DO_PARTIALS ) THEN
        DO UT = 1, N_PARTIALS
          TAU = 0.0D0
          ATTN_P(UT) = 0.0D0
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
        ATTN_FINE(N,J) = 0.0D0
        DO K = 1, NTRAVERSE_FINE(N,J)
          TAU = TAU + SUNPATHS_FINE(N,K,J) * EXTINCTION(K)
        ENDDO
        IF ( TAU .LE. LOCAL_CUTOFF ) ATTN_FINE(N,J) = DEXP(-TAU)
       ENDDO
      ENDDO

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
        DO J = NFINELAYERS, 1, -1
          CALPHA = DCOS(ALPHA_FINE(N,J))
          SALPHA = DSIN(ALPHA_FINE(N,J))
          COT_FINE(J) = CALPHA / SALPHA
          CSQ_FINE(J) = 1.0D0 / SALPHA / SALPHA
        ENDDO

!  INTEGRATED SOURCE TERM MULTIPLIER + TRANSMITTANCE
!   ---- TRAPEZIUM RULE INTEGRATION

        TRAN_1 = DEXP ( - KN * ( COT_1 - COT_2 ) )
        FUNC_1 = ATTN(N-1) * CSQ_1  * TRAN_1
        FUNC_2 = ATTN(N)   * CSQ_2
        SUM = 0.5D0 * ( FUNC_1 + FUNC_2 )
        DO J = NFINELAYERS, 1, -1
          TRAN = DEXP ( -KN * ( COT_FINE(J) - COT_2 ) )
          FUNC = ATTN_FINE(N,J) * TRAN * CSQ_FINE(J)
          SUM = SUM + FUNC
        ENDDO
        LOSTRANS(N)    = TRAN_1
        MULTIPLIERS(N) = SUM * STEP * KN

!        IF (N.EQ.6)WRITE(20,*)'WHOLE',N,MULTIPLIERS(N),LOSTRANS(N)

!  UPDATE GEOMETRY

        CSQ_1 = CSQ_2
        COT_1 = COT_2

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
        TRAN_1 = DEXP ( - KN * ( COT_1 - COT_2 ) )
        FUNC_1 = ATTN(NP-1) * CSQ_1 * TRAN_1
        IF ( NFINE_P.LT.NFINELAYERS) THEN
          SUM = 0.5D0 * FUNC_1
          DO J = NFINELAYERS, NFINE_P+1, -1
            TRAN = DEXP ( -KN * ( COT_FINE(J) - COT_2 ) )
            FUNC = ATTN_FINE(NP,J) * TRAN * CSQ_FINE(J)
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
        MULTIPLIERS_P(UT) = ( TERM1 + TERM2 ) * KN

!        WRITE(20,*)'PARTIAL DN',UT,NP,NFINE_P,
!     &           MULTIPLIERS_P(UT),LOSTRANS_P(UT)

!  FINISH PARTIALS LOOP

      ENDDO

!  FINISH

      RETURN
      END SUBROUTINE OUTGOING_INTEGRATION_DN

!

      SUBROUTINE VLIDORT_DBCORRECTION ( &
        FLUXMULT, &
        DO_SSCORR_OUTGOING, DO_REFRACTIVE_GEOMETRY, &
        DO_UPWELLING, DO_OBSERVATION_GEOMETRY, &
        NSTOKES, NLAYERS, &
        SZA_LOCAL_INPUT, N_USER_RELAZMS, &
        N_USER_LEVELS, &
        DO_LAMBERTIAN_SURFACE, &
        LAMBERTIAN_ALBEDO, EXACTDB_BRDFUNC, &
        DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, &
        SLTERM_ISOTROPIC, SLTERM_USERANGLES, &
        FLUXVEC, COS_SZANGLES, &
        NBEAMS, N_USER_STREAMS, &
        MUELLER_INDEX, SZANGLES_ADJUST, &
        PARTLAYERS_OUTFLAG, &
        PARTLAYERS_OUTINDEX, UTAU_LEVEL_MASK_UP, &
        PARTLAYERS_LAYERIDX, N_GEOMETRIES, &
        VZA_OFFSETS, &
        TRANS_SOLAR_BEAM, DO_REFLECTED_DIRECTBEAM, &
        T_DELT_USERM, T_UTUP_USERM, &
        UP_LOSTRANS, UP_LOSTRANS_UT, BOA_ATTN, &
        DB_CUMSOURCE, EXACTDB_SOURCE, &
        ATTN_DB_SAVE, STOKES_DB )

!  PREPARES EXACT DIRECT BEAM REFLECTION

!  Surface leaving option and variables added 17 May 2012

      USE VLIDORT_PARS

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT (IN) ::  FLUXMULT

      LOGICAL, INTENT (IN) ::           DO_SSCORR_OUTGOING
      LOGICAL, INTENT (IN) ::           DO_REFRACTIVE_GEOMETRY
      LOGICAL, INTENT (IN) ::           DO_UPWELLING
      LOGICAL, INTENT (IN) ::           DO_OBSERVATION_GEOMETRY
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NLAYERS
      DOUBLE PRECISION, INTENT (IN) ::  SZA_LOCAL_INPUT &
          ( 0:MAXLAYERS, MAX_SZANGLES )
      INTEGER, INTENT (IN) ::           N_USER_RELAZMS
      INTEGER, INTENT (IN) ::           N_USER_LEVELS
      LOGICAL, INTENT (IN) ::           DO_LAMBERTIAN_SURFACE
      DOUBLE PRECISION, INTENT (IN) ::  LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) ::  EXACTDB_BRDFUNC &
          ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  New Surface-Leaving stuff 17 May 2012
      LOGICAL, INTENT (IN) ::            DO_SURFACE_LEAVING
      LOGICAL, INTENT (IN) ::            DO_SL_ISOTROPIC
      DOUBLE PRECISION, INTENT (IN) ::   SLTERM_ISOTROPIC &
          ( MAXSTOKES, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::   SLTERM_USERANGLES &
          ( MAXSTOKES, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS  )

      DOUBLE PRECISION, INTENT (IN) ::  FLUXVEC ( MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  COS_SZANGLES ( MAX_SZANGLES )
      INTEGER, INTENT (IN) ::           NBEAMS
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      INTEGER, INTENT (IN) ::           MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  SZANGLES_ADJUST &
          ( MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS )
      LOGICAL, INTENT (IN) ::           PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      INTEGER, INTENT (IN) ::           N_GEOMETRIES
      INTEGER, INTENT (IN) ::           VZA_OFFSETS &
          ( MAX_SZANGLES, MAX_USER_VZANGLES )
      DOUBLE PRECISION, INTENT (IN) ::  TRANS_SOLAR_BEAM ( MAXBEAMS )
      LOGICAL, INTENT (IN) ::           DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_UTUP_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  UP_LOSTRANS &
          ( MAXLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (IN) ::  UP_LOSTRANS_UT &
          ( MAX_PARTLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION, INTENT (IN) ::  BOA_ATTN ( MAX_GEOMETRIES )

      DOUBLE PRECISION, INTENT (OUT) :: DB_CUMSOURCE &
          ( MAX_GEOMETRIES, MAXSTOKES, 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) :: EXACTDB_SOURCE &
          ( MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: ATTN_DB_SAVE ( MAX_GEOMETRIES )

      DOUBLE PRECISION, INTENT (INOUT) :: STOKES_DB &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  LOCAL VARIABLES
!  ---------------

      INTEGER ::           N, NUT, NSTART, NUT_PREV, NLEVEL, NELEMENTS
      INTEGER ::           UT, UTA, UM, UA, NC, IB, V, O1, O2, OM, LUM, LUA
      DOUBLE PRECISION ::  FINAL_SOURCE, TR, X0_FLUX, X0_BOA, ATTN, SUM

!  FIRST STAGE
!  -----------

!  Local user indices

      LUM = 1
      LUA = 1

!  INITIALIZE

      DO V = 1, N_GEOMETRIES
        DO O1 = 1, NSTOKES
          EXACTDB_SOURCE(V,O1) = ZERO
          DO UTA = 1, N_USER_LEVELS
            STOKES_DB(UTA,V,O1)  = ZERO
          ENDDO
        ENDDO
      ENDDO

!  RETURN IF NO UPWELLING

      IF ( .NOT.DO_UPWELLING ) RETURN

!  REFLECTION OF SOLAR BEAM
!  ------------------------

!  MUST USE ADJUSTED VALUES HERE.

!  NEW CODE  R. SPURR, 6 AUGUST 2007. RT SOLUTIONS INC.
!  ====================================================

!  Lattice calculation

      IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN

!  START GEOMETRY LOOPS

        DO UM = 1, N_USER_STREAMS
          DO IB = 1, NBEAMS
            IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
              DO UA = 1, N_USER_RELAZMS
                V = VZA_OFFSETS(IB,UM) + UA

!  BEAM ATTENUATION

                IF ( DO_SSCORR_OUTGOING ) THEN
                  X0_BOA = DCOS(SZANGLES_ADJUST(UM,IB,UA)*DEG_TO_RAD)
                  ATTN = BOA_ATTN(V)
                ELSE
                  IF ( DO_REFRACTIVE_GEOMETRY ) THEN
                    X0_BOA = DCOS(SZA_LOCAL_INPUT(NLAYERS,IB)*DEG_TO_RAD)
                  ELSE
                    X0_BOA = COS_SZANGLES(IB)
                  ENDIF
                  ATTN = TRANS_SOLAR_BEAM(IB)
                ENDIF
                X0_FLUX = FOUR * X0_BOA
                ATTN    = ATTN * X0_FLUX
                ATTN_DB_SAVE(V) = ATTN

!  DEFINE EXACT DB SOURCE TERM

                IF ( .NOT. DO_LAMBERTIAN_SURFACE ) THEN
                  DO O1 = 1, NSTOKES
                    SUM = ZERO
                    DO O2 = 1, NSTOKES
                      OM = MUELLER_INDEX(O1,O2)
                      SUM = SUM + FLUXVEC(O2)*EXACTDB_BRDFUNC(OM,UM,UA,IB)
                    ENDDO
                    EXACTDB_SOURCE(V,O1) = ATTN * SUM
                  ENDDO
                ELSE
                  O1 = 1
                  EXACTDB_SOURCE(V,O1) = ATTN * LAMBERTIAN_ALBEDO * FLUXVEC(O1)
                ENDIF

!  New Surface-Leaving stuff 17 May 2012
!   Multiply terms by 4pi to get correct normalization. 7/30/12

                IF ( DO_SURFACE_LEAVING ) THEN
                  IF ( DO_SL_ISOTROPIC ) THEN
                    O1 = 1
!                    EXACTDB_SOURCE(V,O1) = EXACTDB_SOURCE(V,O1) &
!                      + SLTERM_ISOTROPIC(O1,IB)
                    EXACTDB_SOURCE(V,O1) = EXACTDB_SOURCE(V,O1) &
                      + SLTERM_ISOTROPIC(O1,IB) * PI4
                  ELSE
                    DO O1 = 1, NSTOKES
                      EXACTDB_SOURCE(V,O1) = EXACTDB_SOURCE(V,O1) &
                        + SLTERM_USERANGLES(O1,UM,UA,IB) * PI4
                    ENDDO
                  ENDIF
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

!  START GEOMETRY LOOPS

        DO V = 1, N_GEOMETRIES
          IF ( DO_REFLECTED_DIRECTBEAM(V) ) THEN

!  BEAM ATTENUATION

            IF ( DO_SSCORR_OUTGOING ) THEN
              X0_BOA = DCOS(SZANGLES_ADJUST(LUM,V,LUA)*DEG_TO_RAD)
              ATTN = BOA_ATTN(V)
            ELSE
              IF ( DO_REFRACTIVE_GEOMETRY ) THEN
                X0_BOA = DCOS(SZA_LOCAL_INPUT(NLAYERS,V)*DEG_TO_RAD)
              ELSE
                X0_BOA = COS_SZANGLES(V)
              ENDIF
              ATTN = TRANS_SOLAR_BEAM(V)
            ENDIF
            X0_FLUX = FOUR * X0_BOA
            ATTN    = ATTN * X0_FLUX
            ATTN_DB_SAVE(V) = ATTN

!  DEFINE EXACT DB SOURCE TERM

            IF ( .NOT. DO_LAMBERTIAN_SURFACE ) THEN
              DO O1 = 1, NSTOKES
                SUM = ZERO
                DO O2 = 1, NSTOKES
                  OM = MUELLER_INDEX(O1,O2)
                  SUM = SUM + FLUXVEC(O2)*EXACTDB_BRDFUNC(OM,LUM,LUA,V)
                ENDDO
                EXACTDB_SOURCE(V,O1) = ATTN * SUM
              ENDDO
            ELSE
              O1 = 1
              EXACTDB_SOURCE(V,O1) = ATTN * LAMBERTIAN_ALBEDO * FLUXVEC(O1)
            ENDIF

!  New Surface-Leaving stuff 17 May 2012
!   Multiply terms by 4pi to get correct normalization. 7/30/12

            IF ( DO_SURFACE_LEAVING ) THEN
              IF ( DO_SL_ISOTROPIC ) THEN
                O1 = 1
!                EXACTDB_SOURCE(V,O1) = EXACTDB_SOURCE(V,O1) &
!                  + SLTERM_ISOTROPIC(O1,V)
                EXACTDB_SOURCE(V,O1) = EXACTDB_SOURCE(V,O1) &
                  + SLTERM_ISOTROPIC(O1,V)*PI4
              ELSE
                DO O1 = 1, NSTOKES
                  EXACTDB_SOURCE(V,O1) = EXACTDB_SOURCE(V,O1) &
                    + SLTERM_USERANGLES(O1,LUM,LUA,V) * PI4
                ENDDO
              ENDIF
            ENDIF

!  FINISH LOOP OVER GEOMETRIES

          ENDIF
        ENDDO

!  End observationl geometry calculation

      ENDIF

!  UPWELLING RECURRENCE: TRANSMITTANCE OF EXACT SOURCE TERM
!  --------------------------------------------------------

!  INITIALIZE CUMULATIVE SOURCE TERM = REFLEC * F.MU_0.T/PI
!    T = ATTENUATION OF DIRECT BEAM TO BOA, F = FLUX, REFLEC = ALBEDO/BRDF

      NC = 0
      DO O1 = 1, NSTOKES
        DO V = 1, N_GEOMETRIES
          DB_CUMSOURCE(V,O1,NC) = FLUXMULT * EXACTDB_SOURCE(V,O1)
        ENDDO
      ENDDO

!  NUMBER OF ELEMENTS

      IF ( DO_LAMBERTIAN_SURFACE ) THEN
        NELEMENTS = 1
      ELSE
        NELEMENTS = NSTOKES
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
                    TR = UP_LOSTRANS(N,V)
                  ELSE
                    TR = T_DELT_USERM(N,UM)
                  ENDIF
                  DO O1 = 1, NELEMENTS
                    DB_CUMSOURCE(V,O1,NC) = TR * DB_CUMSOURCE(V,O1,NC-1)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ELSE
            DO V = 1, N_GEOMETRIES
              IF ( DO_SSCORR_OUTGOING ) THEN
                TR = UP_LOSTRANS(N,V)
              ELSE
                TR = T_DELT_USERM(N,V)
              ENDIF
              DO O1 = 1, NELEMENTS
                DB_CUMSOURCE(V,O1,NC) = TR * DB_CUMSOURCE(V,O1,NC-1)
              ENDDO
            ENDDO
          ENDIF

!  END LAYER LOOP

        ENDDO

!  OFFGRID OUTPUT : PARTIAL LAYER TRANSMITTANCE, THEN SET RESULT

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)

          IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                DO UA = 1, N_USER_RELAZMS
                  V = VZA_OFFSETS(IB,UM) + UA
                  IF ( DO_SSCORR_OUTGOING ) THEN
                    TR = UP_LOSTRANS_UT(UT,V)
                  ELSE
                    TR = T_UTUP_USERM(UT,UM)
                  ENDIF
                  DO O1 = 1, NELEMENTS
                    STOKES_DB(UTA,V,O1) = TR * DB_CUMSOURCE(V,O1,NC)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ELSE
            DO V = 1, N_GEOMETRIES
              IF ( DO_SSCORR_OUTGOING ) THEN
                TR = UP_LOSTRANS_UT(UT,V)
              ELSE
                TR = T_UTUP_USERM(UT,V)
              ENDIF
              DO O1 = 1, NELEMENTS
                STOKES_DB(UTA,V,O1) = TR * DB_CUMSOURCE(V,O1,NC)
              ENDDO
            ENDDO
          ENDIF

!  ONGRID OUTPUT : SET FINAL CUMULATIVE SOURCE DIRECTLY

        ELSE

          DO V = 1, N_GEOMETRIES
            DO O1 = 1, NELEMENTS
              FINAL_SOURCE = DB_CUMSOURCE(V,O1,NC)
              STOKES_DB(UTA,V,O1) = FINAL_SOURCE
            ENDDO
          ENDDO

!          DO IB = 1, NBEAMS
!            DO UM = 1, N_USER_STREAMS
!              DO UA = 1, N_USER_RELAZMS
!                V = VZA_OFFSETS(IB,UM) + UA
!                DO O1 = 1, NELEMENTS
!                  FINAL_SOURCE = DB_CUMSOURCE(V,O1,NC)
!                  STOKES_DB(UTA,V,O1) = FINAL_SOURCE
!               ENDDO
!              ENDDO
!            ENDDO
!          ENDDO

        ENDIF

!  CHECK FOR UPDATING THE RECURSION

        IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
        NUT_PREV = NUT

!  END OPTICAL DEPTH LOOP

      ENDDO

!  FINISH

      RETURN
      END SUBROUTINE VLIDORT_DBCORRECTION

!

      SUBROUTINE VFO_MASTER_INTERFACE ( &
        DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,    & ! Input flags
        DO_PLANE_PARALLEL, DO_SSCORR_NADIR, DO_SSCORR_OUTGOING,        & ! Input flags
        DO_DELTAM_SCALING, DO_UPWELLING, DO_DNWELLING,                 & ! Input flags
        DO_LAMBERTIAN_SURFACE, DO_OBSERVATION_GEOMETRY,                & ! Input flags
        NSTOKES, NLAYERS, NFINELAYERS, NMOMENTS, NSTREAMS,             & ! Input numbers
        NGREEK_MOMENTS_INPUT,                                          & ! Input numbers
        N_SZANGLES, SZANGLES, N_USER_VZANGLES, USER_VZANGLES,          & ! Input geometry
        N_USER_RELAZMS, USER_RELAZMS,                                  & ! Input geometry
        N_USER_LEVELS, USER_LEVELS, EARTH_RADIUS, HEIGHT_GRID,         & ! Input other
        SS_FLUX_MULTIPLIER, FLUXVEC,                                   & ! Inputs (Optical - Regular)
        DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,    & ! Inputs (Optical - Regular)
        DELTAU_VERT, OMEGA_TOTAL, GREEKMAT_TOTAL, TRUNC_FACTOR,        & ! Inputs (Optical - Regular)
        THERMAL_BB_INPUT,                                              & ! Inputs (Optical - Regular)
        LAMBERTIAN_ALBEDO, EXACTDB_BRDFUNC, SURFBB, USER_EMISSIVITY,   & ! Inputs (Optical - Surface)
        FO_STOKES_SS, FO_STOKES_DB, FO_STOKES_DTA, FO_STOKES_DTS,      & ! Output
        FO_STOKES, FAIL, MESSAGE, TRACE )                                ! Output

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

      INTEGER, INTENT(IN) ::            NSTOKES
      INTEGER, INTENT(IN) ::            NLAYERS
      INTEGER, INTENT(IN) ::            NFINELAYERS

      INTEGER, INTENT(IN) ::            NMOMENTS
      INTEGER, INTENT(IN) ::            NSTREAMS
      INTEGER, INTENT(IN) ::            NGREEK_MOMENTS_INPUT

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

!  Outputs
!  =======

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

!  Subroutine outputs
!  ------------------

      !REAL(FFP) :: FO_STOKES ( MAXSTOKES, MAX_GEOMETRIES )

!  Other variables

      INTEGER   :: G, GK, L, N, NS2, O1, O2
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

! @@ Rob 7/30/13
!   If deltam-scaling, must use scaled deltau !!!

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
        do g=1,ngeoms
          reflec(1,1,g) = LAMBERTIAN_ALBEDO
        end do
      else
        !Note: set for ObsGeo mode only at present
        do g=1,ngeoms
          do o1=1,nstokes
            do o2=1,nstokes
              gk = 4*(o1-1) + o2
              reflec(o1,o2,g) = EXACTDB_BRDFUNC(gk,lum,lua,g)
            end do
          end do
        end do
      end if

      !surfbb =
      do g=1,ngeoms
        emiss(g) = USER_EMISSIVITY(1,g)
      end do

!  Call VFO_MASTER
!  ===============

      CALL VFO_MASTER &
       ( maxgeoms, maxlayers, maxfine, maxmoments_input, max_user_levels,    & ! Input max dims
         ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, nstokes, & ! Input dims
         do_solar_sources, do_thermal_emission, do_surface_emission,         & ! Input flags
         do_planpar, do_regular_ps, do_enhanced_ps, do_deltam_scaling,       & ! Input flags
         do_upwelling, do_dnwelling, do_lambertian, do_sunlight,             & ! Input flags
         dtr, Pie, doCrit, Acrit, eradius, heights, fo_user_levels,          & ! Input general
         theta_boa, alpha_boa, phi_boa, Mu0, Mu1, doNadir,                   & ! Input geometry
         flux, fluxvec, extinction, deltaus, omega, greekmat, truncfac,      & ! Input (Optical - Regular)
         bb_input, reflec, surfbb, emiss,                                    & ! Input (Optical - Surface)
         fo_stokes_ss, fo_stokes_db, fo_stokes_dta, fo_stokes_dts,           & ! Output
         fo_stokes, fail, message, trace )                                     ! Output

      END SUBROUTINE VFO_MASTER_INTERFACE

      END MODULE vlidort_corrections

