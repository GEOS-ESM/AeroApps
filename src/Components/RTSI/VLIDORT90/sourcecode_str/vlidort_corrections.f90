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
! #                                                             #
! #      Version 2.0. Nadir single scatter correction           #
! #                                                             #
! #            VLIDORT_SSCORR_NADIR (master)                    #
! #              VLIDORTSS_FMATRICES                            #
! #              VLIDORTSS_FMATRICES_MULTI                      #
! #                                                             #
! #      Version 2.2. outgoing sphericity correction            #
! #              2.3.  partial-layer integration                #
! #                                                             #
! #            VLIDORT_SSCORR_OUTGOING (master)                 #
! #                 SSCORR_OUTGOING_ZMATRIX                     #
! #                 OUTGOING_SPHERGEOM_FINE_UP                  #
! #                 OUTGOING_SPHERGEOM_FINE_DN                  #
! #                 OUTGOING_INTEGRATION_UP                     #
! #                 OUTGOING_INTEGRATION_DN                     #
! #                 MULTI_OUTGOING_ADJUSTGEOM                   #
! #                                                             #
! #            VLIDORT_DBCORRECTION                             #
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

      PRIVATE
      PUBLIC :: VLIDORT_SSCORR_NADIR, &
                VLIDORT_SSCORR_OUTGOING, &
                VLIDORT_DBCORRECTION, &
                OUTGOING_SPHERGEOM_FINE_UP, &
                OUTGOING_SPHERGEOM_FINE_DN, &
                MULTI_OUTGOING_ADJUSTGEOM

      CONTAINS

      SUBROUTINE VLIDORT_SSCORR_NADIR ( &
        SSFLUX, &
        DO_SSCORR_TRUNCATION, DO_REFRACTIVE_GEOMETRY, &
        DO_DELTAM_SCALING, DO_UPWELLING, &
        DO_DNWELLING, NSTOKES, &
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
      INTEGER ::           UT, UTA, UM, IA, NC, IB, V, O1, O2, NM1
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
!    - assumes the Flux vector is (F,0,0,0)
!   - drop this variables when further testing has been done

      LOGICAL, PARAMETER :: SUNLIGHT = .TRUE.
      INTEGER ::            NPOLAR

!  Set up operations
!  -----------------

!  Set npolar. For Natural light, there are only 3 Stokes parameters
!    ( no circular polarization for single scattering)

      IF (SUNLIGHT) NPOLAR = MIN(NSTOKES,3)
      IF (.NOT.SUNLIGHT) NPOLAR = NSTOKES

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

      DO UM = 1, N_USER_STREAMS
        CALPHA(UM) = USER_STREAMS(UM)
        SALPHA(UM) = DSQRT ( ONE - CALPHA(UM) * CALPHA(UM) )
      ENDDO

      DO IA = 1, N_USER_RELAZMS
        CPHI(IA) = DCOS ( USER_RELAZMS(IA) * DEG_TO_RAD )
      ENDDO

!  Upwelling single scatter Phase matrices
!  ---------------------------------------

      IF ( DO_UPWELLING ) THEN

!  Get Rotation cosine/sines and F matrices for all layers
!    Distinguish between refractive case and straightline case

        VIEWSIGN = -1.0D0
        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
          CALL VLIDORTSS_FMATRICES_MULTI ( &
            DO_SSCORR_TRUNCATION, NLAYERS, &
            N_GEOMETRIES, NGREEK_MOMENTS_INPUT, &
            NBEAMS, N_USER_STREAMS,  N_USER_RELAZMS, &
            LAYER_MAXMOMENTS,  GREEKMAT_TOTAL_INPUT, SSFDEL, &
            VIEWSIGN, VZA_OFFSETS, &
            CTHETA, STHETA, CALPHA, SALPHA, CPHI, USER_RELAZMS, &
            C1, S1, C2, S2, FMAT_UP )
        ELSE
          CALL VLIDORTSS_FMATRICES ( &
            DO_SSCORR_TRUNCATION, NLAYERS, &
            N_GEOMETRIES, NGREEK_MOMENTS_INPUT, &
            NBEAMS, N_USER_STREAMS,  N_USER_RELAZMS, &
            LAYER_MAXMOMENTS,  GREEKMAT_TOTAL_INPUT, SSFDEL, &
            VIEWSIGN, VZA_OFFSETS, &
            CTHETA, STHETA, CALPHA, SALPHA, CPHI, USER_RELAZMS, &
            C1, S1, C2, S2, FMAT_UP )
        ENDIF

!   Z matrices, Hovenier and van der Mee (1983), Eq 88.
!   ---------------------------------------------------

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
!  this next part remain untested. Need to check signs.
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
!         write(97,'(2i4,1pe15.7)')V,22,ZMAT_UP(V,22,1,1
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
!            write(34,'(2i5,1pe22.12)')V,UTA,STOKES_SS(UTA,V,O1,UPIDX)
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

!  end optical depth loop and Upwelling clause

        ENDDO
      ENDIF

!  Downwelling single scatter Phase matrices
!  -----------------------------------------

      IF ( DO_DNWELLING ) THEN

!  Get Rotation cosine/sines and F matrices for all layers
!    Distinguish between refractive case and straightline case

        VIEWSIGN = +1.0D0
        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
          CALL VLIDORTSS_FMATRICES_MULTI &
        ( DO_SSCORR_TRUNCATION, NLAYERS, &
          N_GEOMETRIES, NGREEK_MOMENTS_INPUT, &
          NBEAMS, N_USER_STREAMS,  N_USER_RELAZMS, &
          LAYER_MAXMOMENTS,  GREEKMAT_TOTAL_INPUT, SSFDEL, &
          VIEWSIGN, VZA_OFFSETS, &
          CTHETA, STHETA, CALPHA, SALPHA, CPHI, USER_RELAZMS, &
          C1, S1, C2, S2, FMAT_DN )
        ELSE
          CALL VLIDORTSS_FMATRICES &
        ( DO_SSCORR_TRUNCATION, NLAYERS, &
          N_GEOMETRIES, NGREEK_MOMENTS_INPUT, &
          NBEAMS, N_USER_STREAMS,  N_USER_RELAZMS, &
          LAYER_MAXMOMENTS,  GREEKMAT_TOTAL_INPUT, SSFDEL, &
          VIEWSIGN, VZA_OFFSETS, &
          CTHETA, STHETA, CALPHA, SALPHA, CPHI, USER_RELAZMS, &
          C1, S1, C2, S2, FMAT_DN )
        ENDIF

!   Z matrices, Hovenier and van der Mee (1983), Eq 88.

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

!  end optical depth loop and Downwelling clause

        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_SSCORR_NADIR

!

      SUBROUTINE VLIDORTSS_FMATRICES &
        ( DO_SSCORR_TRUNCATION, NLAYERS, &
          N_GEOMETRIES, NGREEKMOMS, &
          NBEAMS, N_USER_STREAMS,  N_USER_RELAZMS, &
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

!  control integers

      INTEGER, INTENT(IN) ::          NLAYERS, N_GEOMETRIES, NGREEKMOMS
      INTEGER, INTENT(IN) ::          NBEAMS, N_USER_STREAMS,  N_USER_RELAZMS

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
!    write(8,'(i5,1p4e15.7)')V,C1(V,1),S1(V,1),C2(V,1),S2(V,1)

!  End geometry loops

          ENDDO
        ENDDO
      ENDDO

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
           GK12(N) =  GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(3))   /FT
           GK44(N) = (GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(6))-F)/FT
           GK34(N) =  GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(5))   /FT
           GK22(N) = (GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(2))-F)/FT
           GK33(N) = (GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(4))-F)/FT
          ENDDO
        ELSE
          DO N = 1, NLAYERS
           GK11(N) = GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(1))
           GK12(N) = GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(3))
           GK44(N) = GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(6))
           GK34(N) = GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(5))
           GK22(N) = GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(2))
           GK33(N) = GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(4))
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
        ( DO_SSCORR_TRUNCATION, NLAYERS, &
          N_GEOMETRIES, NGREEKMOMS, &
          NBEAMS, N_USER_STREAMS,  N_USER_RELAZMS, &
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

!  control integers

      INTEGER, INTENT(IN) ::           NLAYERS, N_GEOMETRIES, NGREEKMOMS
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
           GK12(N) =  GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(3))   /FT
           GK44(N) = (GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(6))-F)/FT
           GK34(N) =  GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(5))   /FT
           GK22(N) = (GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(2))-F)/FT
           GK33(N) = (GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(4))-F)/FT
          ENDDO
        ELSE
          DO N = 1, NLAYERS
           GK11(N) = GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(1))
           GK12(N) = GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(3))
           GK44(N) = GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(6))
           GK34(N) = GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(5))
           GK22(N) = GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(2))
           GK33(N) = GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(4))
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
        DO_UPWELLING, DO_DNWELLING, &
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

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT (IN) ::  SSFLUX

      LOGICAL, INTENT (IN) ::           DO_SSCORR_TRUNCATION
      LOGICAL, INTENT (IN) ::           DO_DELTAM_SCALING
      LOGICAL, INTENT (IN) ::           DO_UPWELLING
      LOGICAL, INTENT (IN) ::           DO_DNWELLING
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
      INTEGER ::           UT, UTA, UM, IA, NC, IB, V, O1, O2, NM1

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
!            DO N = 0, NLAYERS
!              WRITE(45,'(I4,101F10.5)')N,(SUNPATHS(N,V),V=1,NLAYERS)
!            ENDDO
!            PAUSE

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

!      DO N = 1, NLAYERS
!        WRITE(85,'(I4,1P2E18.10)')N,UP_LOSTRANS(N,1),UP_MULTIPLIERS(N,1)
!      ENDDO

!  DEBUG
!              WRITE(45,*)V
!              DO N = 1, NLAYERS
!                WRITE(45,'(1P2E18.10)')
!     &           UP_LOSTRANS(N,V),UP_MULTIPLIERS(N,V)
!              ENDDO
!               PAUSE

!              IF (V.EQ.1)WRITE(35,*)V
!              DO UT = 1, N_PARTLAYERS
!                IF (V.EQ.1)WRITE(35,'(1P2E18.10)')
!     &           UP_LOSTRANS_UT(UT,V),UP_MULTIPLIERS_UT(UT,V)
!              ENDDO

!  DEBUG: MULTIPLIERS (ALL LAYERS), UPWELLING OLD WAY
!              CALL OUTGOING_INTEGRATION_OLD_UP
!     I           ( NLAYERS, EXTINCTION, DELTAU_VERT,
!     I             LOSPATHS, SUNPATHS, NTRAVERSE,
!     O             UP_MULTIPLIERS(1,V), UP_LOSTRANS(1,V) )
!              WRITE(46,*)V
!              DO N = 1, NLAYERS
!                WRITE(46,*)UP_LOSTRANS(N,V),UP_MULTIPLIERS(N,V)
!              ENDDO

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
!              WRITE(65,*)V
!              DO N = 1, NLAYERS
!                WRITE(65,'(1P2E18.10)')
!     &           DN_LOSTRANS(N,V),DN_MULTIPLIERS(N,V)
!              ENDDO

!  DEBUG: MULTIPLIERS (ALL LAYERS), DOWN WELLING OLD WAY
!              CALL OUTGOING_INTEGRATION_OLD_DN
!     I           ( NLAYERS, EXTINCTION, DELTAU_VERT,
!     I             LOSPATHS, SUNPATHS, NTRAVERSE,
!     O             DN_MULTIPLIERS(1,V), DN_LOSTRANS(1,V) )
!              WRITE(66,*)V
!              DO N = 1, NLAYERS
!                WRITE(66,*)DN_LOSTRANS(N,V),DN_MULTIPLIERS(N,V)
!              ENDDO

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
          GK12 =   GREEKMATRIX(L,N,GREEKMAT_INDEX(3)) / FACT
          IF ( .NOT. DO_SUNLIGHT ) THEN
           GK44 = ( GREEKMATRIX(L,N,GREEKMAT_INDEX(6)) - FDNL1 ) / FACT
           GK34 =   GREEKMATRIX(L,N,GREEKMAT_INDEX(5)) / FACT
           GK22 = ( GREEKMATRIX(L,N,GREEKMAT_INDEX(2)) - FDNL1 ) / FACT
           GK33 = ( GREEKMATRIX(L,N,GREEKMAT_INDEX(4)) - FDNL1 ) / FACT
          ENDIF
        ELSE
          GK11 = GREEKMATRIX(L,N,GREEKMAT_INDEX(1))
          GK12 = GREEKMATRIX(L,N,GREEKMAT_INDEX(3))
          IF ( .NOT. DO_SUNLIGHT ) THEN
           GK44 = GREEKMATRIX(L,N,GREEKMAT_INDEX(6))
           GK34 = GREEKMATRIX(L,N,GREEKMAT_INDEX(5))
           GK22 = GREEKMATRIX(L,N,GREEKMAT_INDEX(2))
           GK33 = GREEKMATRIX(L,N,GREEKMAT_INDEX(4))
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

      SUBROUTINE OUTGOING_SPHERGEOM_FINE_UP &
        ( MAXLAYERS, MAXFINE, MAXPARTIALS, &
          DO_FINE, DO_PARTIALS, NLAYERS, NFINE, &
          N_PARTIALS, PARTIALS_IDX, &
          HEIGHTS, HEIGHTS_P, ERADIUS, &
          ALPHA_BOA, THETA_BOA, PHI_BOA, &
          SUNPATHS,      RADII,      NTRAVERSE,      ALPHA_ALL, &
          SUNPATHS_FINE, RADII_FINE, NTRAVERSE_FINE, ALPHA_FINE, &
          SUNPATHS_P,    RADII_P,    NTRAVERSE_P,    ALPHA_P, &
          PARTIALS_FINEIDX, LOSPATHS, LOSPATHS_P_UP, &
          THETA_ALL, PHI_ALL, COSSCAT_UP, &
          FAIL, MESSAGE )

      IMPLICIT NONE

!  COMPLETELY STAND-ALONE GEOMETRY ROUTINE FOR THE OUTGOING CORRECTION
!   THIS IS APPLICABLE TO THE UPWELLING PATH GEOEMTRY

!    STARTING INPUTS ARE THE BOA VALUES OF SZA, VZA AND PHI
!    NEED ALSO THE HEIGHT GRIDS, EARTH RADIUS AND CONTROL

!  THIS ROUTINE HAS THE FINE GRIDDING TREATMENT
!  VERSION 2.3. SEPTEMBER 2007, PARTIAL LAYER GEOMETRIES ADDED

!  INPUTS

      INTEGER, INTENT(IN) ::           MAXLAYERS, MAXFINE, MAXPARTIALS
      INTEGER, INTENT(IN) ::           NLAYERS, NFINE
      LOGICAL, INTENT(IN) ::           DO_FINE, DO_PARTIALS
      INTEGER, INTENT(IN) ::           N_PARTIALS
      INTEGER, INTENT(IN) ::           PARTIALS_IDX (MAXPARTIALS)
      DOUBLE PRECISION, INTENT(IN) ::  ERADIUS, HEIGHTS (0:MAXLAYERS)
      DOUBLE PRECISION, INTENT(IN) ::  HEIGHTS_P(MAXPARTIALS)
      DOUBLE PRECISION, INTENT(IN) ::  ALPHA_BOA, THETA_BOA

!  INPUT/OUTPUT

      DOUBLE PRECISION, INTENT(INOUT) ::  PHI_BOA

!  MAIN OUTPUTS (GEOMETRY)

      INTEGER, INTENT(OUT) ::           NTRAVERSE   (0:MAXLAYERS)
      DOUBLE PRECISION, INTENT(OUT) ::  SUNPATHS   (0:MAXLAYERS,MAXLAYERS)
      DOUBLE PRECISION, INTENT(OUT) ::  RADII      (0:MAXLAYERS)
      DOUBLE PRECISION, INTENT(OUT) ::  ALPHA_ALL  (0:MAXLAYERS)

!  FINE LEVEL OUTPUT (GEOMETRY)

      INTEGER, INTENT(OUT) ::           NTRAVERSE_FINE(MAXLAYERS,MAXFINE)
      DOUBLE PRECISION, INTENT(OUT) ::  SUNPATHS_FINE & 
          (MAXLAYERS,MAXLAYERS,MAXFINE)
      DOUBLE PRECISION, INTENT(OUT) ::  RADII_FINE    (MAXLAYERS,MAXFINE)
      DOUBLE PRECISION, INTENT(OUT) ::  ALPHA_FINE    (MAXLAYERS,MAXFINE)

!  PARTIAL LAYER OUTPUT (GEOMETRY)

      INTEGER, INTENT(OUT) ::           NTRAVERSE_P(MAXPARTIALS)
      DOUBLE PRECISION, INTENT(OUT) ::  SUNPATHS_P (MAXPARTIALS,MAXLAYERS)
      DOUBLE PRECISION, INTENT(OUT) ::  RADII_P    (MAXPARTIALS)
      DOUBLE PRECISION, INTENT(OUT) ::  ALPHA_P    (MAXPARTIALS)

!  OTHER (INCIDENTAL) GEOMETRICAL OUTPUT

      DOUBLE PRECISION, INTENT(OUT) ::  LOSPATHS(MAXLAYERS)
      DOUBLE PRECISION, INTENT(OUT) ::  LOSPATHS_P_UP(MAXPARTIALS)
      INTEGER, INTENT(OUT) ::           PARTIALS_FINEIDX(MAXPARTIALS)

      DOUBLE PRECISION, INTENT(OUT) ::  THETA_ALL  (0:MAXLAYERS)
      DOUBLE PRECISION, INTENT(OUT) ::  PHI_ALL    (0:MAXLAYERS)
      DOUBLE PRECISION, INTENT(OUT) ::  COSSCAT_UP (0:MAXLAYERS)

!  STATUS OUTPUT

      LOGICAL, INTENT(OUT) ::           FAIL
      CHARACTER (LEN=*), INTENT(INOUT) :: MESSAGE

!  LOCAL

      LOGICAL ::           DIRECT_SUN
      INTEGER ::           N, K, KRAD, N1, UT, NP
      DOUBLE PRECISION ::  DEG_TO_RAD, EX, EY, EZ, PX, PY, PZ
      DOUBLE PRECISION ::  SALPHA_BOA, CALPHA_BOA, SPHI_BOA, PIE, PI2
      DOUBLE PRECISION ::  STHETA_BOA, CTHETA_BOA, CPHI_BOA
      DOUBLE PRECISION ::  KSI, CKSI, SKSI, XICUM, TANGR, FAC
      DOUBLE PRECISION ::  CTHETA, STHETA, CALPHA, SALPHA, CPHI
      DOUBLE PRECISION ::  B, STH0, TH0, KS1, STH1, TH1, THETA_P

!  LOCAL ARRAYS ASSOCIATED WITH FINE GRID OUTPUT

      INTEGER ::            J
      INTEGER, PARAMETER :: MAXLOCALFINE = 20
      LOGICAL ::            DIRECT_SUNF(MAXLOCALFINE)
      DOUBLE PRECISION ::   DIFZ, DFINE1, SAF, XICUM0, PATH
      DOUBLE PRECISION ::   THETAF(MAXLOCALFINE), XICUMF, DIFA
      DOUBLE PRECISION ::   CTHETAF(MAXLOCALFINE)
      DOUBLE PRECISION ::   STHETAF(MAXLOCALFINE)
      DOUBLE PRECISION ::   KSIF(MAXLOCALFINE)

!  INITIALISE OUTPUT

      FAIL = .FALSE.
      MESSAGE = ' '

!  CHECK RANGE OF INPUTS

      IF ( ALPHA_BOA.GE.90.0D0.OR.ALPHA_BOA.LT.0.0D0 ) THEN
        MESSAGE = 'BOA LOS ANGLE OUTSIDE RANGE [0,90])'
        FAIL    = .TRUE.
        RETURN
      ENDIF
      IF ( PHI_BOA.LT.0.0D0 )   PHI_BOA = - PHI_BOA

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  BUG CORRECTED 01 OCTOBER 2010------------------------- RTS, R. SPURR

!  IMPORTANT, DO NOT LIMIT THE AZIMUTH TO < 180 DEGREES

!      IF ( PHI_BOA.GT.180.0D0 ) PHI_BOA = 360.0D0 - PHI_BOA    ! OLD
       IF ( PHI_BOA.GT.360.0D0 ) PHI_BOA = PHI_BOA - 360.0D0    ! NEW

!  BUG CORRECTED 01 OCTOBER 2010------------------------- RTS, R. SPURR
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      IF ( THETA_BOA.GE.90.0D0.OR.THETA_BOA.LT.0.0D0 ) THEN
        MESSAGE = 'BOA SZA ANGLE OUTSIDE RANGE [0,90])'
        FAIL    = .TRUE.
        RETURN
      ENDIF
      IF ( DO_FINE ) THEN
       IF ( NFINE.GT.MAXLOCALFINE ) THEN
         MESSAGE = 'LOCAL FINELAYER DIMENSIONING INSUFFICIENT'
         FAIL    = .TRUE.
         RETURN
       ENDIF
      ENDIF

!  ZERO THE SUN PATHS
!  INITIALIZE NUMBER OF LAYERS TRAVERSED  (NOMINAL CONDITIONS)

      DO N = 0, NLAYERS
        NTRAVERSE(N) = N
        DO K = 1, NLAYERS
         SUNPATHS(N,K) = 0.0D0
        ENDDO
      ENDDO

!  ZERO THE PARTIAL PATHS, INITIALIZE THE TRAVERSE NUMBER (NOMINAL)

      IF ( DO_PARTIALS ) THEN
        DO UT = 1, N_PARTIALS
          NP = PARTIALS_IDX(UT)
          NTRAVERSE_P(UT) = NP
          DO K = 1, NP
            SUNPATHS_P(UT,K) = 0.0D0
          ENDDO
        ENDDO
      ENDIF

!  ZERO THE FINE DATA PATHS

      IF ( DO_FINE ) THEN
       DFINE1 = DBLE(NFINE) + 1
       DO N = 1, NLAYERS
        DO J = 1, NFINE
         NTRAVERSE_FINE(N,J) = N
         DO K = 1, NLAYERS
          SUNPATHS_FINE(N,K,J) = 0.0D0
         ENDDO
        ENDDO
       ENDDO
      ENDIF

!  START AT BOA

      PIE        = DACOS(-1.0D0)
      DEG_TO_RAD = PIE / 180.0D0
      PI2        = 2.0D0 * PIE

      ALPHA_ALL(NLAYERS) = ALPHA_BOA * DEG_TO_RAD
      THETA_ALL(NLAYERS) = THETA_BOA * DEG_TO_RAD
      PHI_ALL(NLAYERS)   = PHI_BOA   * DEG_TO_RAD

!  COSINE OF SCATTERING ANGLE AT BOA

      SALPHA_BOA = DSIN(ALPHA_ALL(NLAYERS))
      CALPHA_BOA = DCOS(ALPHA_ALL(NLAYERS))
      STHETA_BOA = DSIN(THETA_ALL(NLAYERS))
      CTHETA_BOA = DCOS(THETA_ALL(NLAYERS))
      CPHI_BOA   = DCOS(PHI_ALL(NLAYERS))
      SPHI_BOA   = DSIN(PHI_ALL(NLAYERS))

      COSSCAT_UP (NLAYERS) = - CALPHA_BOA * CTHETA_BOA + &
                               SALPHA_BOA * STHETA_BOA * CPHI_BOA

!  RADII
!  -----

!  LAYER LEVELS

      DO N = 0, NLAYERS
        RADII(N) = ERADIUS + HEIGHTS(N)
      ENDDO

!  FINE LEVELS

      IF ( DO_FINE ) THEN
        DO N = 1, NLAYERS
          DIFZ = (RADII(N-1)-RADII(N))/DFINE1
          DO J = 1, NFINE
            RADII_FINE(N,J) = RADII(N) + DIFZ * DBLE(J)
          ENDDO
        ENDDO
      ENDIF

!  PARTIAL LEVELS
!   FIND THE FINE LEVEL JUST BELOW PARTIAL POINT = (0,1,2,...NFINE)

      IF ( DO_PARTIALS ) THEN
        DO UT = 1, N_PARTIALS
          NP = PARTIALS_IDX(UT)
          RADII_P(UT) = ERADIUS + HEIGHTS_P(UT)
          J = 1
          DO WHILE (RADII_P(UT).GT.RADII_FINE(NP,J))
             J = J + 1
             IF (J.GT.NFINE) GO TO 5778
          ENDDO
 5778     CONTINUE
          PARTIALS_FINEIDX(UT) = J - 1
        ENDDO
      ENDIF
!          J = 1
!          DO WHILE (RADII_P(UT).GT.RADII_FINE(NP,J).AND.J.LE.NFINE)
!           J = J + 1
!          ENDDO

!  SPECIAL CASE. DIRECT NADIR VIEWING
!  ==================================

!  COMPUTE EVERYTHING AND EXIT.
!    (THIS IS THE SAME AS THE REGULAR PSEUDO-SPHERICAL )

      IF ( SALPHA_BOA.EQ.0.0D0 ) THEN

!  WHOLE LAYER AND FINE DIVISIONS
!  ------------------------------

!  START LAYER LOOP, WORKING UPWARDS

        DO N = NLAYERS,1,-1

!  SET MAIN OUTPUT.

          ALPHA_ALL(N-1)   = ALPHA_ALL(N)
          THETA_ALL(N-1)   = THETA_ALL(N)
          PHI_ALL(N-1)     = PHI_ALL(N)
          COSSCAT_UP(N-1) = COSSCAT_UP(N)
          LOSPATHS(N) = RADII(N-1)-RADII(N)
          IF ( DO_FINE ) THEN
            DO J = 1, NFINE
              ALPHA_FINE(N,J) = 0.0D0
            ENDDO
          ENDIF

!  OVERHEAD SUN

          IF (STHETA_BOA.EQ.0.0D0 ) THEN
            DO K = N, 1, -1
              SUNPATHS(N,K) = RADII(K-1)-RADII(K)
            ENDDO
            IF ( DO_FINE ) THEN
              DO J = 1, NFINE
                DO K = N - 1, 1, -1
                  SUNPATHS_FINE(N,K,J) = RADII(K-1)-RADII(K)
                ENDDO
                SUNPATHS_FINE(N,N,J) = RADII(N-1)-RADII_FINE(N,J)
              ENDDO
            ENDIF
          ENDIF

!  NON-OVERHEAD SUN
!  MAIN OUTPUT OF SOLAR PATHS
!  SOLAR PATH DISTANCES FOR FINE OUTPUT

          IF (STHETA_BOA.GT.0.0D0 ) THEN
            STH0 = STHETA_BOA
            TH0  = THETA_ALL(N)
            DO K = N, 1, -1
              STH1 = STH0*RADII(K)/RADII(K-1)
              TH1  = DASIN(STH1)
              KS1  = TH0-TH1
              SUNPATHS(N,K) = DSIN(KS1)*RADII(K)/STH1
              STH0 = STH1
              TH0  = TH1
            ENDDO
            IF ( DO_FINE ) THEN
              DO J = 1, NFINE
                STH0 = STHETA_BOA
                TH0  = THETA_ALL(N)
                STH1 = STH0*RADII_FINE(N,J)/RADII(N-1)
                TH1  = DASIN(STH1)
                KS1  = TH0-TH1
                SUNPATHS_FINE(N,N,J) = DSIN(KS1)*RADII_FINE(N,J)/STH1
                STH0 = STH1
                TH0  = TH1
                DO K = N-1, 1, -1
                  STH1 = STH0*RADII(K)/RADII(K-1)
                  TH1  = DASIN(STH1)
                  KS1  = TH0-TH1
                  SUNPATHS_FINE(N,K,J) = DSIN(KS1)*RADII(K)/STH1
                  STH0 = STH1
                  TH0  = TH1
                ENDDO
              ENDDO
            ENDIF
          ENDIF

!  END MAIN LAYER LOOP

        ENDDO

!  PARTIAL LAYERS
!  --------------

        IF ( DO_PARTIALS ) THEN
          DO UT = 1, N_PARTIALS
            NP = PARTIALS_IDX(UT)

!  LOS ANGLE AND PATHS

            ALPHA_P(UT) = ALPHA_ALL(NLAYERS)
            LOSPATHS_P_UP(UT) = RADII_P(UT)-RADII(NP)

!  OVERHEAD SUN

            IF (STHETA_BOA.EQ.0.0D0 ) THEN
              DO K = 1, NP-1
                SUNPATHS_P(UT,K) = RADII(K-1)-RADII(K)
              ENDDO
              SUNPATHS_P(UT,NP) = RADII(NP-1)-RADII_P(UT)
            ENDIF

!  NON-OVERHEAD SUN

            IF (STHETA_BOA.GT.0.0D0 ) THEN
              STH0 = STHETA_BOA
              TH0  = THETA_ALL(NP)
              STH1 = STH0*RADII_P(UT)/RADII(NP-1)
              TH1  = DASIN(STH1)
              KS1  = TH0-TH1
              SUNPATHS_P(UT,NP) = DSIN(KS1)*RADII_P(UT)/STH1
              STH0 = STH1
              TH0  = TH1
              DO K = NP-1, 1, -1
                STH1 = STH0*RADII(K)/RADII(K-1)
                TH1  = DASIN(STH1)
                KS1  = TH0-TH1
                SUNPATHS_P(UT,K) = DSIN(KS1)*RADII(K)/STH1
                STH0 = STH1
                TH0  = TH1
              ENDDO
            ENDIF

!  FINISH PARTIAL LAYER STUFF

          ENDDO
        ENDIF

!  RETURN, AS EVERYTHING NOW DONE

        RETURN

!  END REGULAR PSEUDO-SPHERICAL CLAUSE, LOS IS ZERO

      ENDIF

!  OUTGOING SPHERICITY GEOMETRY
!  ============================

!  DEFINE UNIT SOLAR VECTOR AT BOA

      EX = - STHETA_BOA * CPHI_BOA
      EY = - STHETA_BOA * SPHI_BOA
      EZ = - CTHETA_BOA

!  SUN PATHS, BOA GEOMETRY, ALWAYS DIRECTLY ILLUMINATED

      IF ( STHETA_BOA.EQ.0.0D0 ) THEN
        DO K = NLAYERS, 1, -1
          SUNPATHS(NLAYERS,K) = RADII(K-1)-RADII(K)
        ENDDO
      ELSE
        STH0 = STHETA_BOA
        TH0  = THETA_ALL(NLAYERS)
        DO K = NLAYERS, 1, -1
          STH1 = STH0*RADII(K)/RADII(K-1)
          TH1  = DASIN(STH1)
          KS1  = TH0-TH1
          SUNPATHS(NLAYERS,K) = DSIN(KS1)*RADII(K)/STH1
          STH0 = STH1
          TH0  = TH1
        ENDDO
      ENDIF

!  CHECK SINGLE ILLUMINATION
!      --NOT REQUIRED, NOW WE HAVE THE TANGENT POINT TREATMENT
!      IF (STHETA_BOA.GT.0.0D0 ) THEN
!        XICUM = DASIN(RADII(NLAYERS)*SALPHA_BOA/RADII(0))
!        XICUM = ALPHA_ALL(NLAYERS)-XICUM
!        PX = - RADII(0) * DSIN(XICUM)
!        PY = 0.0D0
!        PZ =   RADII(0) * DCOS(XICUM)
!        B = EX*PX + EY*PY + EZ*PZ
!        CTHETA = -B/RADII(0)
!        IF ( CTHETA.LE.0.0D0 ) THEN
!          WRITE(*,*)'LIMIT VALUE = ',90.0D0-XICUM/DEG_TO_RAD
!        ENDIF
!      ENDIF

!  INITIALISE LOS CUMULATIVE ANGLE

      XICUM  = 0.0D0

!  SET TOA DIRECT ILLUMINATION FLAG

      DIRECT_SUN = .TRUE.
      IF ( DO_FINE ) THEN
        DO J = 1, NFINE
          DIRECT_SUNF(J) = .TRUE.
        ENDDO
      ENDIF

!  START LOOP OVER POSITIONS (LAYER UPPER BOUNDARIES)

      DO N = NLAYERS - 1, 0, -1

!  NEXT LEVEL UP

        N1 = N + 1

!  LOS ANGLES AT LEVEL BOUNDARIES

        SALPHA = RADII(NLAYERS) * SALPHA_BOA / RADII(N)
        ALPHA_ALL(N)  = DASIN(SALPHA)
        CALPHA = DCOS(ALPHA_ALL(N))

!  LOSPATHS

        KSI = ALPHA_ALL(N1) - ALPHA_ALL(N)
        SKSI = DSIN(KSI)
        CKSI = DCOS(KSI)
        LOSPATHS(N1) = SKSI * RADII(N1) / SALPHA
        XICUM0 = XICUM
        XICUM  = XICUM + KSI

!  FINE GRID LOSPATH OUTPUT (ANGLE AND RADIUS)
!    LOCALLY SAVE THE EARTH-CENTER ANGLE KSIF

        IF ( DO_FINE ) THEN
          DIFA = (ALPHA_ALL(N1)-ALPHA_ALL(N))/DFINE1
          DO J = 1, NFINE
            ALPHA_FINE(N1,J) = ALPHA_ALL(N1) - DIFA * DBLE(J)
            SAF = DSIN(ALPHA_FINE(N1,J))
            RADII_FINE(N1,J) = SALPHA_BOA * RADII(NLAYERS) / SAF
            KSIF(J) = ALPHA_ALL(N1) - ALPHA_FINE(N1,J)
          ENDDO
        ENDIF

!  SUN ANGLES FOR THE DIRECT NADIR CASE

        IF (STHETA_BOA.EQ.0.0D0 ) THEN
         THETA_ALL(N) = XICUM
         CTHETA = DCOS(THETA_ALL(N))
         STHETA = DSQRT(1.0D0-CTHETA*CTHETA)
         IF ( DO_FINE ) THEN
           DO J = 1, NFINE
             THETAF(J)  = XICUM0 + KSIF(J)
             CTHETAF(J) = DCOS(THETAF(J))
             STHETAF(J) = DSQRT(1.0D0-CTHETA*CTHETA)
           ENDDO
         ENDIF
        ENDIF

!  SUN ANGLES FOR THE GENERAL CASE
!    LOCAL SAVE OF ANGLES, COSINES, SINES AND  ILLUMINATION FLAGS

        IF (STHETA_BOA.GT.0.0D0 ) THEN
         PX = - RADII(N) * DSIN(XICUM)
         PY = 0.0D0
         PZ =   RADII(N) * DCOS(XICUM)
         B = EX*PX + EY*PY + EZ*PZ
         CTHETA = -B/RADII(N)
         DIRECT_SUN = (DIRECT_SUN.AND.CTHETA.GE.0.D0)
         STHETA = DSQRT(1.0D0-CTHETA*CTHETA)
         THETA_ALL(N) = DACOS(CTHETA)
         IF ( DO_FINE ) THEN
           DO J = 1, NFINE
             XICUMF  = XICUM0 + KSIF(J)
             PX = - RADII_FINE(N1,J) * DSIN(XICUMF)
             PY = 0.0D0
             PZ =   RADII_FINE(N1,J) * DCOS(XICUMF)
             B  = EX*PX + EY*PY + EZ*PZ
             CTHETAF(J) = -B/RADII_FINE(N1,J)
             DIRECT_SUNF(J) = (DIRECT_SUNF(J).AND.CTHETAF(J).GE.0.D0)
             STHETAF(J) = DSQRT(1.0D0-CTHETAF(J)*CTHETAF(J))
             THETAF(J)  = DACOS(CTHETAF(J))
           ENDDO
         ENDIF
        ENDIF

!  UNIT VECTOR F2(I) PERPENDICULAR TO OP BUT IN PLANE OF PATH
!  PROJECTION OF F2(I) ON SOLAR PATH GIVES THE RELATIVE AZIMUTH AT P
!        F2X = DSIN(XICUM)
!        F2Y = 0.0D0
!        F2Z = DCOS(XICUM)
!        CPHI = - (EX*F2X + EY*F2Y + EZ*F2Z ) / STHETA
!        CPHI = - (EX*F2X + EY*F2Y + EZ*F2Z ) / STHETA
!        IF ( CPHI.GT.1.0D0)  CPHI = 1.0D0
!        IF ( CPHI.LT.-1.0D0) CPHI = -1.0D0
! ********************************************* APPARENTLY NOT CORRECT

!  FIX PHI BY USING CONSTANCY OF SCATTER ANGLE

        COSSCAT_UP(N) = COSSCAT_UP(N+1)
        IF (STHETA_BOA.EQ.0.0D0 ) THEN
          PHI_ALL(N)     = PHI_ALL(N+1)
        ELSE
         CPHI = (COSSCAT_UP(N)+CALPHA*CTHETA)/STHETA/SALPHA
         IF ( CPHI.GT.1.0D0) CPHI = 1.0D0
         IF ( CPHI.LT.-1.0D0) CPHI = -1.0D0
         PHI_ALL(N)     = DACOS(CPHI)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 !  BUG CORRECTED 01 OCTOBER 2010------------------------- RTS, R. SPURR

!  IF AZIMUTH > 180 DEGREES. MUST ENSURE CONSISTENCY, ADD THIS LINE.

         IF ( PHI_BOA.GT.180.0D0 ) PHI_ALL(N) = PI2 - PHI_ALL(N)

 !  B G CORRECTED 01 OCTOBER 2010------------------------- RTS, R. SPURR
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        ENDIF

!  SUN PATHS, DIRECT SUN AT LAYER TOP
!  ==================================

!   MEANS THAT THE SZA AT LAYER TOP IS < 90.
!    ===> SZA < 90 FOR ALL FINE POINTS IN LAYER BENEATH

        IF ( DIRECT_SUN ) THEN

!  WORK UP FROM LEVEL N TO TOA
!    LAYER TOP CALCULATION GETS LEFT OUT AT TOA

         IF ( N .GT. 0 ) THEN
          STH0 = STHETA
          TH0  = THETA_ALL(N)
          DO K = N, 1, -1
           STH1 = STH0*RADII(K)/RADII(K-1)
           TH1  = DASIN(STH1)
           KS1  = TH0-TH1
           SUNPATHS(N,K) = DSIN(KS1)*RADII(K)/STH1
           STH0 = STH1
           TH0  = TH1
          ENDDO
         ENDIF

! DBG         WRITE(*,*)'REGULAR',N1,(SUNPATHS(N1,K),K=N1,1,-1)

!  FINE GRID CALCULATION IS ALWAYS REQUIRED
!  ----------------------------------------

!   START AT THE GRID POINT ON LOS PATH AND WORK UPWARDS,
!     FIRST ACROSS PARTIAL LAYER TO UPPER BOUNDARY, THEN CONTINUE
!     UPWARDS TO TOA ON A WHOLE LAYER BASIS. SINE RULE.

         IF ( DO_FINE ) THEN
          DO J = 1, NFINE
           STH0 = STHETAF(J)
           TH0  = THETAF(J)
           STH1 = STH0*RADII_FINE(N1,J)/RADII(N)
           TH1  = DASIN(STH1)
           KS1  = TH0-TH1
           SUNPATHS_FINE(N1,N1,J) = DSIN(KS1)*RADII_FINE(N1,J)/STH1
           STH0 = STH1
           TH0  = TH1
           DO K = N, 1, -1
            STH1 = STH0*RADII(K)/RADII(K-1)
            TH1  = DASIN(STH1)
            KS1  = TH0-TH1
            SUNPATHS_FINE(N1,K,J) = DSIN(KS1)*RADII(K)/STH1
            STH0 = STH1
            TH0  = TH1
           ENDDO
!  DBG           WRITE(*,*)'REGULAR',N1,J,(SUNPATHS_FINE(N1,K,J),K=N1,1,
          ENDDO

!  DBG         WRITE(*,*)'REGULAR',N,(SUNPATHS(N,K),K=N,1,-1)

         ENDIF

!  COMPLETE DIRECT SUN COMPUTATIONS

        ENDIF

!  SUN PATHS, NOT DIRECT SUN , WITH TANGENT POINT
!  ==============================================

!  ALTHOUGH LAYER TOP HAS A TANGENT POINT, NOT ALL OF THE FINE-GRID
!   POINTS WILL HAVE A TANGENT POINT.

        IF (.NOT.DIRECT_SUN ) THEN

!  FIRST DO THE LAYER-TOP CALCULATION.
!  -----------------------------------

!  TANGR = TANGENT POINT RADIUS.

         TANGR = STHETA*RADII(N)

!  NTRAVERSE(N) IS THE NUMBER OF LAYERS TRAVERSED BY RAY.

         KRAD = NLAYERS
         DO WHILE (TANGR.GT.RADII(KRAD))
          KRAD = KRAD - 1
         ENDDO
         NTRAVERSE(N) = KRAD + 1

!  START AT THE TOA ANGLES

         STH0 = TANGR/RADII(0)
         TH0 = DASIN(STH0)

!  WORK DOWNWARDS FROM TOA (SINE RULE) TO LEVEL IMMEDIATELY ABOVE
!  THE TANGENT LAYER. DON'T FORGET TO DOUBLE THE PATH LENGTH FOR
!  ANY LAYERS WHICH ARE TRAVERSED TWICE.

         DO K = 1, KRAD
          STH1 = RADII(K-1)*STH0/RADII(K)
          TH1 = DASIN(STH1)
          KS1 = TH1-TH0
          FAC = 1.0D0
          IF ( K.GT.N) FAC = 2.0D0
          SUNPATHS(N,K) = FAC*DSIN(KS1)*RADII(K-1)/STH1
          STH0 = STH1
          TH0  = TH1
         ENDDO

!  TANGENT LAYER PATH LENGTH. TWICE AGAIN. THE FOLLOWING CHECK IS GOOD.
!  CHECK       WRITE(*,*)TANGR/DTAN(TH1),RADII(KRAD)*DCOS(TH1)

         SUNPATHS(N,KRAD+1)=2.0D0*RADII(KRAD)*DCOS(TH1)

! DBG         WRITE(*,*)'--TANGENT',N1,(SUNPATHS(N1,K),K=NTRAVERSE(N1),1

!  FINE LAYER DEVELOPMENT (SLIGHTLY DIFFERENT)
!  ---------------------

         IF ( DO_FINE ) THEN
          DO J = 1, NFINE

!  IF THERE IS NO TANGENT POINT, REPEAT CALCULATION ABOVE.

           IF ( DIRECT_SUNF(J) ) THEN
            STH0 = STHETAF(J)
            TH0  = THETAF(J)
            STH1 = STH0*RADII_FINE(N1,J)/RADII(N)
            TH1  = DASIN(STH1)
            KS1  = TH0-TH1
            SUNPATHS_FINE(N1,N1,J) = DSIN(KS1)*RADII_FINE(N1,J)/STH1
            STH0 = STH1
            TH0  = TH1
            DO K = N, 1, -1
             STH1 = STH0*RADII(K)/RADII(K-1)
             TH1  = DASIN(STH1)
             KS1  = TH0-TH1
             SUNPATHS_FINE(N1,K,J) = DSIN(KS1)*RADII(K)/STH1
             STH0 = STH1
             TH0  = TH1
            ENDDO

! DBG           WRITE(*,*)'REGULAR',N1,J,
!     &       (SUNPATHS_FINE(N1,K,J),K=NTRAVERSE_FINE(N1,J),1,-1)

!  FINE GRID CALCULATION WITH TANGENT POINT
!  ----------------------------------------

           ELSE

!  LOCAL TANGENT RADIUS AND NUMBER OF LAYERS TRAVERSED

            TANGR = STHETAF(J)*RADII_FINE(N1,J)
            KRAD = NLAYERS
            DO WHILE (TANGR.GT.RADII(KRAD))
             KRAD = KRAD - 1
            ENDDO
            NTRAVERSE_FINE(N1,J) = KRAD + 1

!  START AGAIN AT TOA

            STH0 = TANGR/RADII(0)
            TH0  = DASIN(STH0)

!  WORK DOWN TO THE LEVEL N

            DO K = 1, N
             STH1 = RADII(K-1)*STH0/RADII(K)
             TH1 = DASIN(STH1)
             KS1 = TH1-TH0
             SUNPATHS_FINE(N1,K,J) = DSIN(KS1)*RADII(K-1)/STH1
             STH0 = STH1
             TH0  = TH1
            ENDDO

!  IN LAYER BELOW LEVEL N, WORK DOWN TO LEVEL OF FINE GRID RADIUS
!     (SINGLE CONTRIBUTION)

            STH1 = RADII(N)*STH0/RADII_FINE(N1,J)
            TH1 = DASIN(STH1)
            KS1 = TH1-TH0
            PATH = DSIN(KS1)*RADII(N)/STH1
            STH0 = STH1
            TH0  = TH1

!  IN LAYER BELOW LEVEL N, COMPLETE THE PATH DOWN TO THE TANGENT POINT
!    (DOUBLE CONTRIBUTION). FINISH.

            IF ( KRAD.EQ.N ) THEN

              PATH = PATH + 2.0D0*RADII_FINE(N1,J)*DCOS(TH1)
              SUNPATHS_FINE(N1,KRAD+1,J) = PATH

!  CHECK    WRITE(*,*)'1',TANGR/DTAN(TH1),RADII_FINE(N1,J)*DCOS(TH1)

!  IN LAYER BELOW LEVEL N, NEED TO GO DOWN FURTHER.
!    (FROM NOW ON, WE NEED DOUBLE CONTRIBUTIONS)
!     --- CONTINUE THE PATH TO THE BOTTOM OF THIS LAYER
!     --- CONTINUE DOWN BY WHOLE-LAYER STEPS UNTIL REACH LEVEL ABOVE TAN
!     --- COMPLETE THE PATH DOWN TO THE TANGENT POINT

            ELSE
              STH1 = RADII_FINE(N1,J)*STH0/RADII(N1)
              TH1 = DASIN(STH1)
              KS1 = TH1-TH0
              PATH = PATH + 2.0D0 * DSIN(KS1)*RADII_FINE(N1,J)/STH1
              SUNPATHS_FINE(N1,N1,J) = PATH
              STH0 = STH1
              TH0  = TH1
              DO K = N1 + 1, KRAD
               STH1 = RADII(K-1)*STH0/RADII(K)
               TH1 = DASIN(STH1)
               KS1 = TH1-TH0
               SUNPATHS_FINE(N1,K,J) = DSIN(KS1)*RADII(K-1)/STH1
               STH0 = STH1
               TH0  = TH1
              ENDDO
              SUNPATHS_FINE(N1,KRAD+1,J)=2.0D0*RADII(KRAD)*DCOS(TH1)
!  CHECK 2    WRITE(*,*)'2',TANGR/DTAN(TH1),RADII(KRAD)*DCOS(TH1)
            ENDIF

! DBG           WRITE(*,*)'TANGENT',N1,J,
!     &       (SUNPATHS_FINE(N1,K,J),K=NTRAVERSE_FINE(N1,J),1,-1)

!  COMPLETE TANGENT POINT CLAUSE

           ENDIF

!  COMPLETE FINE-GRID LOOP

          ENDDO

!  DBG        WRITE(*,*)'TANGENT',N,(SUNPATHS(N,K),K=NTRAVERSE(N),1,-1)

!  COMPLETE TANGENT POINT CALCULATION

         ENDIF

        ENDIF

!  END LAYER LOOP

      ENDDO

!  PARTIALS SECTION
!  ----------------

!  PARTIAL ANGLES AND LOSPATHS

      IF ( DO_PARTIALS ) THEN
        DO UT = 1, N_PARTIALS
          NP = PARTIALS_IDX(UT)

!  ESTABLISH THE LOS ANGLE, AND LOSPATHS

          TH0 = ALPHA_ALL(NP)
          STH0 = DSIN(TH0)
          SALPHA = RADII(NP)*STH0 / RADII_P(UT)
          ALPHA_P(UT) = DASIN(SALPHA)
          CALPHA = DSQRT(1.0D0-SALPHA*SALPHA)

          KSI = ALPHA_ALL(NP) - ALPHA_P(UT)
          LOSPATHS_P_UP(UT) = DSIN(KSI)*RADII(NP)/SALPHA

!          KSI = ALPHA_P(UT) - ALPHA_ALL(NP-1)
!          LOSPATHS_P_DN(UT) = DSIN(KSI)*RADII(NP-1)/SALPHA

!  CHECK DEBUG
!          LOSPATHS_P_DN(UT) = LOSPATHS(NP)-LOSPATHS_P_UP(UT)
!          WRITE(*,*)LOSPATHS_P_UP(UT),LOSPATHS_P_DN(UT),
!     &              RADII_P(UT)-RADII(NP)
!          PAUSE

!  ESTABLISH THE SOLAR ANGLE, BOTH NAIDR =SUN-AT-BOA AND GENERALLY

          XICUM = ALPHA_ALL(NLAYERS) - ALPHA_P(UT)
          DIRECT_SUN = .TRUE.
          IF (STHETA_BOA.EQ.0.0D0 ) THEN
            THETA_P = XICUM
            CTHETA = DCOS(THETA_P)
            STHETA = DSQRT(1.0D0-CTHETA*CTHETA)
          ELSE
            PX = - RADII_P(UT) * DSIN(XICUM)
            PY = 0.0D0
            PZ =   RADII_P(UT) * DCOS(XICUM)
            B = EX*PX + EY*PY + EZ*PZ
            CTHETA = -B/RADII_P(UT)
            DIRECT_SUN = (DIRECT_SUN.AND.CTHETA.GE.0.D0)
            STHETA = DSQRT(1.0D0-CTHETA*CTHETA)
            THETA_P = DACOS(CTHETA)
          ENDIF

!         WRITE(*,*)DIRECT_SUN,THETA_ALL(NP-1),THETA_P,THETA_ALL(NP)

!  SUN PATHS, DIRECT SUN AT LAYER TOP
!   FIRST CALCULATE SUNPATH FOR LAYER CONTAINING PARTIAL POINT
!   THEN WORK UPWARDS TO TOA USING THE SINE RULE.

          IF ( DIRECT_SUN ) THEN
            STH0 = STHETA
            TH0  = THETA_P
            STH1 = STH0*RADII_P(UT)/RADII(NP-1)
            TH1  = DASIN(STH1)
            KS1  = TH0-TH1
            SUNPATHS_P(UT,NP) = DSIN(KS1)*RADII_P(UT)/STH1
            STH0 = STH1
            TH0  = TH1
            DO K = NP-1, 1, -1
              STH1 = STH0*RADII(K)/RADII(K-1)
              TH1  = DASIN(STH1)
              KS1  = TH0-TH1
              SUNPATHS_P(UT,K) = DSIN(KS1)*RADII(K)/STH1
              STH0 = STH1
              TH0  = TH1
            ENDDO
          ENDIF

!  DEBUG
! 14       FORMAT(A11,1X,10F10.6)
!          J = PARTIALS_FINEIDX(UT)
!          WRITE(*,*)UT,NTRAVERSE_P(UT)
!          WRITE(*,*)UT,J,NP,LOSPATHS_P_UP(UT)/LOSPATHS(NP)
!           WRITE(*,14)'REGULAR FEB',(SUNPATHS(NP,K),K=NP,1,-1)
!           DO J = 1, NFINE
!             WRITE(*,14)'REGULAR FEB',(SUNPATHS_FINE(NP,K,J),K=NP,1,-1)
!             IF (J.EQ.PARTIALS_FINEIDX(UT))
!     &        WRITE(*,14)'REGULAR UT ',(SUNPATHS_P(UT,K),K=NP,1,-1)
!           ENDDO
!           WRITE(*,14)'REGULAR FEB',0.0D0,(SUNPATHS(NP-1,K),K=NP-1,1,-1
!
!  SUN PATHS, NOT DIRECT SUN , WITH TANGENT POINT. CODE SIMILAR TO ABOVE
!  TANGR = TANGENT POINT RADIUS FOR THIS RAY
!   NTRAVERSE_P(UT) = NUMBER OF LAYERS TRAVERSED BY RAY.

!  WORK DOWNWARDS FROM TOA (SINE RULE) TO LEVEL IMMEDIATELY ABOVE
!  THE TANGENT LAYER. DON'T FORGET TO DOUBLE THE PATH LENGTH FOR
!  ANY LAYERS WHICH ARE TRAVERSED TWICE.

!  TANGENT LAYER PATH LENGTH. TWICE AGAIN. THE FOLLOWING CHECK IS GOOD.
!  CHECK       WRITE(*,*)TANGR/DTAN(TH1),RADII(KRAD)*DCOS(TH1)

          IF (.NOT.DIRECT_SUN ) THEN
            TANGR = STHETA*RADII(N)
            KRAD = NLAYERS
            DO WHILE (TANGR.GT.RADII(KRAD))
              KRAD = KRAD - 1
            ENDDO
            NTRAVERSE_P(UT) = KRAD + 1
            STH0 = TANGR/RADII(0)
            TH0 = DASIN(STH0)
            DO K = 1, KRAD
              STH1 = RADII(K-1)*STH0/RADII(K)
              TH1 = DASIN(STH1)
              KS1 = TH1-TH0
              FAC = 1.0D0
              IF ( K.GT.NP) FAC = 2.0D0
              SUNPATHS_P(UT,K) = FAC*DSIN(KS1)*RADII(K-1)/STH1
              STH0 = STH1
              TH0  = TH1
            ENDDO
            SUNPATHS_P(UT,KRAD+1)=2.0D0*RADII_P(UT)*DCOS(TH1)
          ENDIF

!  FINISH PARTIALS LOOP

        ENDDO
      ENDIF

!  FINISH

      RETURN
      END SUBROUTINE OUTGOING_SPHERGEOM_FINE_UP

!

      SUBROUTINE OUTGOING_SPHERGEOM_FINE_DN &
        ( MAXLAYERS, MAXFINE, MAXPARTIALS, &
          DO_FINE, DO_PARTIALS, NLAYERS, NFINE, &
          N_PARTIALS, PARTIALS_IDX, &
          HEIGHTS, HEIGHTS_P, ERADIUS, &
          ALPHA_BOA, THETA_BOA, PHI_BOA, &
          SUNPATHS,      RADII,      NTRAVERSE,      ALPHA_ALL, &
          SUNPATHS_FINE, RADII_FINE, NTRAVERSE_FINE, ALPHA_FINE, &
          SUNPATHS_P,    RADII_P,    NTRAVERSE_P,    ALPHA_P, &
          PARTIALS_FINEIDX, LOSPATHS, LOSPATHS_P_DN, &
          THETA_ALL, PHI_ALL, COSSCAT_DN, &
          FAIL, MESSAGE )

      IMPLICIT NONE


!  COMPLETELY STAND-ALONE GEOMETRY ROUTINE FOR THE OUTGOING CORRECTION
!   THIS IS APPLICABLE TO THE DOWNWELLING PATH GEOEMTRY

!    STARTING INPUTS ARE THE BOA VALUES OF SZA, VZA AND PHI
!    NEED ALSO THE HEIGHT GRIDS, EARTH RADIUS AND CONTROL

!  THIS ROUTINE HAS THE FINE GRIDDING TREATMENT
!  VERSION 2.3. SEPTEMBER 2007, PARTIAL LAYER GEOMETRIES ADDED

!  INPUTS

      INTEGER, INTENT(IN) ::           MAXLAYERS, MAXFINE, MAXPARTIALS
      INTEGER, INTENT(IN) ::           NLAYERS, NFINE
      LOGICAL, INTENT(IN) ::           DO_FINE, DO_PARTIALS
      INTEGER, INTENT(IN) ::           N_PARTIALS
      INTEGER, INTENT(IN) ::           PARTIALS_IDX    (MAXPARTIALS)
      DOUBLE PRECISION, INTENT(IN) ::  ERADIUS, HEIGHTS (0:MAXLAYERS)
      DOUBLE PRECISION, INTENT(IN) ::  HEIGHTS_P(MAXPARTIALS)
      DOUBLE PRECISION, INTENT(IN) ::  ALPHA_BOA, THETA_BOA

!  INPUT/OUTPUT

      DOUBLE PRECISION, INTENT(INOUT) :: PHI_BOA

!  MAIN OUTPUTS (GEOMETRY)

      INTEGER, INTENT(OUT) ::           NTRAVERSE   (0:MAXLAYERS)
      DOUBLE PRECISION, INTENT(OUT) ::  SUNPATHS   (0:MAXLAYERS,MAXLAYERS)
      DOUBLE PRECISION, INTENT(OUT) ::  RADII      (0:MAXLAYERS)
      DOUBLE PRECISION, INTENT(OUT) ::  ALPHA_ALL  (0:MAXLAYERS)

!  FINE LEVEL OUTPUT (GEOMETRY)

      INTEGER, INTENT(OUT) ::           NTRAVERSE_FINE(MAXLAYERS,MAXFINE)
      DOUBLE PRECISION, INTENT(OUT) ::  SUNPATHS_FINE & 
          (MAXLAYERS,MAXLAYERS,MAXFINE)
      DOUBLE PRECISION, INTENT(OUT) ::  RADII_FINE    (MAXLAYERS,MAXFINE)
      DOUBLE PRECISION, INTENT(OUT) ::  ALPHA_FINE    (MAXLAYERS,MAXFINE)

!  PARTIAL LAYER OUTPUT (GEOMETRY)

      INTEGER, INTENT(OUT) ::           NTRAVERSE_P(MAXPARTIALS)
      DOUBLE PRECISION, INTENT(OUT) ::  SUNPATHS_P (MAXPARTIALS,MAXLAYERS)
      DOUBLE PRECISION, INTENT(OUT) ::  RADII_P    (MAXPARTIALS)
      DOUBLE PRECISION, INTENT(OUT) ::  ALPHA_P    (MAXPARTIALS)

!  OTHER (INCIDENTAL) GEOMETRICAL OUTPUT

      DOUBLE PRECISION, INTENT(OUT) ::  LOSPATHS(MAXLAYERS)
      DOUBLE PRECISION, INTENT(OUT) ::  LOSPATHS_P_DN(MAXPARTIALS)
      INTEGER, INTENT(OUT) ::           PARTIALS_FINEIDX(MAXPARTIALS)

      DOUBLE PRECISION, INTENT(OUT) ::  THETA_ALL  (0:MAXLAYERS)
      DOUBLE PRECISION, INTENT(OUT) ::  PHI_ALL    (0:MAXLAYERS)
      DOUBLE PRECISION, INTENT(OUT) ::  COSSCAT_DN (0:MAXLAYERS)

!  STATUS OUTPUT

      LOGICAL, INTENT(OUT) ::           FAIL
      CHARACTER (LEN=*), INTENT(INOUT) :: MESSAGE

!  LOCAL

      LOGICAL ::           DIRECT_SUN
      INTEGER ::           N, K, KRAD, N1, UT, NP
      DOUBLE PRECISION ::  DEG_TO_RAD, EX, EY, EZ, PX, PY, PZ
      DOUBLE PRECISION ::  SALPHA_BOA, CALPHA_BOA, SPHI_BOA, PIE, PI2
      DOUBLE PRECISION ::  STHETA_BOA, CTHETA_BOA, CPHI_BOA
      DOUBLE PRECISION ::  KSI, CKSI, SKSI, XICUM, TANGR, FAC
      DOUBLE PRECISION ::  CTHETA, STHETA, CALPHA, SALPHA, CPHI
      DOUBLE PRECISION ::  B, STH0, TH0, KS1, STH1, TH1, THETA_P

!  LOCAL ARRAYS ASSOCIATED WITH FINE GRID OUTPUT

      INTEGER ::            J
      INTEGER, PARAMETER :: MAXLOCALFINE = 20
      LOGICAL ::            DIRECT_SUNF(MAXLOCALFINE)
      DOUBLE PRECISION ::   DIFZ, DFINE1, SAF, XICUM0, PATH
      DOUBLE PRECISION ::   THETAF(MAXLOCALFINE), XICUMF, DIFA
      DOUBLE PRECISION ::   CTHETAF(MAXLOCALFINE)
      DOUBLE PRECISION ::   STHETAF(MAXLOCALFINE)
      DOUBLE PRECISION ::   KSIF(MAXLOCALFINE)

!  INITIALISE OUTPUT

      FAIL = .FALSE.
      MESSAGE = ' '

!  CHECK RANGE OF INPUTS

      IF ( ALPHA_BOA.GE.90.0D0.OR.ALPHA_BOA.LT.0.0D0 ) THEN
        MESSAGE = 'BOA LOS ANGLE OUTSIDE RANGE [0,90])'
        FAIL    = .TRUE.
        RETURN
      ENDIF
      IF ( PHI_BOA.LT.0.0D0 )   PHI_BOA = - PHI_BOA

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  BUG CORRECTED 01 OCTOBER 2010------------------------- RTS, R. SPURR

!  IMPORTANT, DO NOT LIMIT THE AZIMUTH TO < 180 DEGREES

!      IF ( PHI_BOA.GT.180.0D0 ) PHI_BOA = 360.0D0 - PHI_BOA    ! OLD
       IF ( PHI_BOA.GT.360.0D0 ) PHI_BOA = PHI_BOA - 360.0D0    ! NEW

!  BUG CORRECTED 01 OCTOBER 2010------------------------- RTS, R. SPURR
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      IF ( THETA_BOA.GE.90.0D0.OR.THETA_BOA.LT.0.0D0 ) THEN
        MESSAGE = 'BOA SZA ANGLE OUTSIDE RANGE [0,90])'
        FAIL    = .TRUE.
        RETURN
      ENDIF
      IF ( DO_FINE ) THEN
       IF ( NFINE.GT.MAXLOCALFINE ) THEN
         MESSAGE = 'LOCAL FINELAYER DIMENSIONING INSUFFICIENT'
         FAIL    = .TRUE.
         RETURN
       ENDIF
      ENDIF

!  ZERO THE SUN PATHS
!  INITIALIZE NUMBER OF LAYERS TRAVERSED  (NOMINAL CONDITIONS)

      DO N = 0, NLAYERS
        NTRAVERSE(N) = N
        DO K = 1, NLAYERS
         SUNPATHS(N,K) = 0.0D0
        ENDDO
      ENDDO

!  ZERO THE PARTIAL PATHS, INITIALIZE THE TRAVERSE NUMBER (NOMINAL)

      IF ( DO_PARTIALS ) THEN
        DO UT = 1, N_PARTIALS
          NP = PARTIALS_IDX(UT)
          NTRAVERSE_P(UT) = NP
          DO K = 1, NP
            SUNPATHS_P(UT,K) = 0.0D0
          ENDDO
        ENDDO
      ENDIF

!  ZERO THE FINE DATA PATHS

      IF ( DO_FINE ) THEN
       DFINE1 = DBLE(NFINE) + 1
       DO N = 1, NLAYERS
        DO J = 1, NFINE
         NTRAVERSE_FINE(N,J) = N
         DO K = 1, NLAYERS
          SUNPATHS_FINE(N,K,J) = 0.0D0
         ENDDO
        ENDDO
       ENDDO
      ENDIF

!  START AT BOA

      PIE        = DACOS(-1.0D0)
      DEG_TO_RAD = PIE / 180.0D0
      PI2        = 2.0D0 * PIE

      ALPHA_ALL(NLAYERS) = ALPHA_BOA * DEG_TO_RAD
      THETA_ALL(NLAYERS) = THETA_BOA * DEG_TO_RAD
      PHI_ALL(NLAYERS)   = PHI_BOA   * DEG_TO_RAD

!  COSINE OF SCATTERING ANGLE AT BOA

      SALPHA_BOA = DSIN(ALPHA_ALL(NLAYERS))
      CALPHA_BOA = DCOS(ALPHA_ALL(NLAYERS))
      STHETA_BOA = DSIN(THETA_ALL(NLAYERS))
      CTHETA_BOA = DCOS(THETA_ALL(NLAYERS))
      CPHI_BOA   = DCOS(PHI_ALL(NLAYERS))
      SPHI_BOA   = DSIN(PHI_ALL(NLAYERS))

      COSSCAT_DN (NLAYERS) = + CALPHA_BOA * CTHETA_BOA + &
                               SALPHA_BOA * STHETA_BOA * CPHI_BOA

!  RADII
!  -----

!  LAYER LEVELS

      DO N = 0, NLAYERS
        RADII(N) = ERADIUS + HEIGHTS(N)
      ENDDO

!  FINE LEVELS

      IF ( DO_FINE ) THEN
        DO N = 1, NLAYERS
          DIFZ = (RADII(N-1)-RADII(N))/DFINE1
          DO J = 1, NFINE
            RADII_FINE(N,J) = RADII(N) + DIFZ * DBLE(J)
          ENDDO
        ENDDO
      ENDIF

!  PARTIAL LEVELS
!   FIND THE FINE LEVEL JUST BELOW PARTIAL POINT = (0,1,2,...NFINE)

      IF ( DO_PARTIALS ) THEN
        DO UT = 1, N_PARTIALS
          NP = PARTIALS_IDX(UT)
          RADII_P(UT) = ERADIUS + HEIGHTS_P(UT)
          J = 1
          DO WHILE (RADII_P(UT).GT.RADII_FINE(NP,J))
             J = J + 1
             IF (J.GT.NFINE) GO TO 5778
          ENDDO
 5778     CONTINUE
          PARTIALS_FINEIDX(UT) = J - 1
        ENDDO
      ENDIF

!  SPECIAL CASE. DIRECT NADIR VIEWING
!  ==================================

!  COMPUTE EVERYTHING AND EXIT.
!    (THIS IS THE SAME AS THE REGULAR PSEUDO-SPHERICAL )

      IF ( SALPHA_BOA.EQ.0.0D0 ) THEN

!  WHOLE LAYER AND FINE DIVISIONS
!  ------------------------------

!  START LAYER LOOP, WORKING UPWARDS

        DO N = NLAYERS,1,-1

!  SET MAIN OUTPUT.

          ALPHA_ALL(N-1)   = ALPHA_ALL(N)
          THETA_ALL(N-1)   = THETA_ALL(N)
          PHI_ALL(N-1)     = PHI_ALL(N)
          COSSCAT_DN(N-1) = COSSCAT_DN(N)
          LOSPATHS(N) = RADII(N-1)-RADII(N)
          IF ( DO_FINE ) THEN
            DO J = 1, NFINE
              ALPHA_FINE(N,J) = 0.0D0
            ENDDO
          ENDIF

!  OVERHEAD SUN

          IF (STHETA_BOA.EQ.0.0D0 ) THEN
            DO K = N, 1, -1
              SUNPATHS(N,K) = RADII(K-1)-RADII(K)
            ENDDO
            IF ( DO_FINE ) THEN
              DO J = 1, NFINE
                DO K = N - 1, 1, -1
                  SUNPATHS_FINE(N,K,J) = RADII(K-1)-RADII(K)
                ENDDO
                SUNPATHS_FINE(N,N,J) = RADII(N-1)-RADII_FINE(N,J)
              ENDDO
            ENDIF
          ENDIF

!  NON-OVERHEAD SUN
!  MAIN OUTPUT OF SOLAR PATHS
!  SOLAR PATH DISTANCES FOR FINE OUTPUT

          IF (STHETA_BOA.GT.0.0D0 ) THEN
            STH0 = STHETA_BOA
            TH0  = THETA_ALL(N)
            DO K = N, 1, -1
              STH1 = STH0*RADII(K)/RADII(K-1)
              TH1  = DASIN(STH1)
              KS1  = TH0-TH1
              SUNPATHS(N,K) = DSIN(KS1)*RADII(K)/STH1
              STH0 = STH1
              TH0  = TH1
            ENDDO
            IF ( DO_FINE ) THEN
              DO J = 1, NFINE
                STH0 = STHETA_BOA
                TH0  = THETA_ALL(N)
                STH1 = STH0*RADII_FINE(N,J)/RADII(N-1)
                TH1  = DASIN(STH1)
                KS1  = TH0-TH1
                SUNPATHS_FINE(N,N,J) = DSIN(KS1)*RADII_FINE(N,J)/STH1
                STH0 = STH1
                TH0  = TH1
                DO K = N-1, 1, -1
                  STH1 = STH0*RADII(K)/RADII(K-1)
                  TH1  = DASIN(STH1)
                  KS1  = TH0-TH1
                  SUNPATHS_FINE(N,K,J) = DSIN(KS1)*RADII(K)/STH1
                  STH0 = STH1
                  TH0  = TH1
                ENDDO
              ENDDO
            ENDIF
          ENDIF

!  END MAIN LAYER LOOP

        ENDDO

!  PARTIAL LAYERS
!  --------------

        IF ( DO_PARTIALS ) THEN
          DO UT = 1, N_PARTIALS
            NP = PARTIALS_IDX(UT)

!  LOS ANGLE AND PATHS

            ALPHA_P(UT) = ALPHA_ALL(NLAYERS)
            LOSPATHS_P_DN(UT) = RADII(NP-1)-RADII_P(UT)

!  OVERHEAD SUN

            IF (STHETA_BOA.EQ.0.0D0 ) THEN
              DO K = 1, NP-1
                SUNPATHS_P(UT,K) = RADII(K-1)-RADII(K)
              ENDDO
              SUNPATHS_P(UT,NP) = RADII(NP-1)-RADII_P(UT)
            ENDIF

!  NON-OVERHEAD SUN

            IF (STHETA_BOA.GT.0.0D0 ) THEN
              STH0 = STHETA_BOA
              TH0  = THETA_ALL(NP)
              STH1 = STH0*RADII_P(UT)/RADII(NP-1)
              TH1  = DASIN(STH1)
              KS1  = TH0-TH1
              SUNPATHS_P(UT,NP) = DSIN(KS1)*RADII_P(UT)/STH1
              STH0 = STH1
              TH0  = TH1
              DO K = NP-1, 1, -1
                STH1 = STH0*RADII(K)/RADII(K-1)
                TH1  = DASIN(STH1)
                KS1  = TH0-TH1
                SUNPATHS_P(UT,K) = DSIN(KS1)*RADII(K)/STH1
                STH0 = STH1
                TH0  = TH1
              ENDDO
            ENDIF

!  FINISH PARTIAL LAYER STUFF

          ENDDO
        ENDIF

!  RETURN, AS EVERYTHING NOW DONE

        RETURN

!  END REGULAR PSEUDO-SPHERICAL CLAUSE, LOS IS ZERO

      ENDIF

!  OUTGOING SPHERICITY GEOMETRY
!  ============================

!  DEFINE UNIT SOLAR VECTOR AT BOA
!   NOTE CHANGE OF SIGN FOR EX HERE.

      EX = + STHETA_BOA * CPHI_BOA
      EY = - STHETA_BOA * SPHI_BOA
      EZ = - CTHETA_BOA

!  SUN PATHS, BOA GEOMETRY, ALWAYS DIRECTLY ILLUMINATED

      IF ( STHETA_BOA.EQ.0.0D0 ) THEN
        DO K = NLAYERS, 1, -1
          SUNPATHS(NLAYERS,K) = RADII(K-1)-RADII(K)
        ENDDO
      ELSE
        STH0 = STHETA_BOA
        TH0  = THETA_ALL(NLAYERS)
        DO K = NLAYERS, 1, -1
          STH1 = STH0*RADII(K)/RADII(K-1)
          TH1  = DASIN(STH1)
          KS1  = TH0-TH1
          SUNPATHS(NLAYERS,K) = DSIN(KS1)*RADII(K)/STH1
          STH0 = STH1
          TH0  = TH1
        ENDDO
      ENDIF

!  CHECK SINGLE ILLUMINATION
!      --NOT REQUIRED, NOW WE HAVE THE TANGENT POINT TREATMENT
!      IF (STHETA_BOA.GT.0.0D0 ) THEN
!        XICUM = DASIN(RADII(NLAYERS)*SALPHA_BOA/RADII(0))
!        XICUM = ALPHA_ALL(NLAYERS)-XICUM
!        PX = - RADII(0) * DSIN(XICUM)
!        PY = 0.0D0
!        PZ =   RADII(0) * DCOS(XICUM)
!        B = EX*PX + EY*PY + EZ*PZ
!        CTHETA = -B/RADII(0)
!        IF ( CTHETA.LE.0.0D0 ) THEN
!          WRITE(*,*)'LIMIT VALUE = ',90.0D0-XICUM/DEG_TO_RAD
!        ENDIF
!      ENDIF

!  INITIALISE LOS CUMULATIVE ANGLE

      XICUM  = 0.0D0

!  SET TOA DIRECT ILLUMINATION FLAG

      DIRECT_SUN = .TRUE.
      IF ( DO_FINE ) THEN
        DO J = 1, NFINE
          DIRECT_SUNF(J) = .TRUE.
        ENDDO
      ENDIF

!  START LOOP OVER POSITIONS (LAYER UPPER BOUNDARIES)

      DO N = NLAYERS - 1, 0, -1

!  NEXT LEVEL UP

        N1 = N + 1

!  LOS ANGLES AT LEVEL BOUNDARIES

        SALPHA = RADII(NLAYERS) * SALPHA_BOA / RADII(N)
        ALPHA_ALL(N)  = DASIN(SALPHA)
        CALPHA = DCOS(ALPHA_ALL(N))

!  LOSPATHS

        KSI = ALPHA_ALL(N1) - ALPHA_ALL(N)
        SKSI = DSIN(KSI)
        CKSI = DCOS(KSI)
        LOSPATHS(N1) = SKSI * RADII(N1) / SALPHA
        XICUM0 = XICUM
        XICUM  = XICUM + KSI

!  FINE GRID LOSPATH OUTPUT (ANGLE AND RADIUS)
!    LOCALLY SAVE THE EARTH-CENTER ANGLE KSIF

        IF ( DO_FINE ) THEN
          DIFA = (ALPHA_ALL(N1)-ALPHA_ALL(N))/DFINE1
          DO J = 1, NFINE
            ALPHA_FINE(N1,J) = ALPHA_ALL(N1) - DIFA * DBLE(J)
            SAF = DSIN(ALPHA_FINE(N1,J))
            RADII_FINE(N1,J) = SALPHA_BOA * RADII(NLAYERS) / SAF
            KSIF(J) = ALPHA_ALL(N1) - ALPHA_FINE(N1,J)
          ENDDO
        ENDIF

!  SUN ANGLES FOR THE DIRECT NADIR CASE

        IF (STHETA_BOA.EQ.0.0D0 ) THEN
         THETA_ALL(N) = XICUM
         CTHETA = DCOS(THETA_ALL(N))
         STHETA = DSQRT(1.0D0-CTHETA*CTHETA)
         IF ( DO_FINE ) THEN
           DO J = 1, NFINE
             THETAF(J)  = XICUM0 + KSIF(J)
             CTHETAF(J) = DCOS(THETAF(J))
             STHETAF(J) = DSQRT(1.0D0-CTHETA*CTHETA)
           ENDDO
         ENDIF
        ENDIF

!  SUN ANGLES FOR THE GENERAL CASE
!    LOCAL SAVE OF ANGLES, COSINES, SINES AND  ILLUMINATION FLAGS

        IF (STHETA_BOA.GT.0.0D0 ) THEN
         PX = - RADII(N) * DSIN(XICUM)
         PY = 0.0D0
         PZ =   RADII(N) * DCOS(XICUM)
         B = EX*PX + EY*PY + EZ*PZ
         CTHETA = -B/RADII(N)
         DIRECT_SUN = (DIRECT_SUN.AND.CTHETA.GE.0.D0)
         STHETA = DSQRT(1.0D0-CTHETA*CTHETA)
         THETA_ALL(N) = DACOS(CTHETA)
         IF ( DO_FINE ) THEN
           DO J = 1, NFINE
             XICUMF  = XICUM0 + KSIF(J)
             PX = - RADII_FINE(N1,J) * DSIN(XICUMF)
             PY = 0.0D0
             PZ =   RADII_FINE(N1,J) * DCOS(XICUMF)
             B  = EX*PX + EY*PY + EZ*PZ
             CTHETAF(J) = -B/RADII_FINE(N1,J)
             DIRECT_SUNF(J) = (DIRECT_SUNF(J).AND.CTHETAF(J).GE.0.D0)
             STHETAF(J) = DSQRT(1.0D0-CTHETAF(J)*CTHETAF(J))
             THETAF(J)  = DACOS(CTHETAF(J))
           ENDDO
         ENDIF
        ENDIF

!  UNIT VECTOR F2(I) PERPENDICULAR TO OP BUT IN PLANE OF PATH
!  PROJECTION OF F2(I) ON SOLAR PATH GIVES THE RELATIVE AZIMUTH AT P
!        F2X = DSIN(XICUM)
!        F2Y = 0.0D0
!        F2Z = DCOS(XICUM)
!        CPHI = - (EX*F2X + EY*F2Y + EZ*F2Z ) / STHETA
!        CPHI = - (EX*F2X + EY*F2Y + EZ*F2Z ) / STHETA
!        IF ( CPHI.GT.1.0D0)  CPHI = 1.0D0
!        IF ( CPHI.LT.-1.0D0) CPHI = -1.0D0
! ********************************************* APPARENTLY NOT CORRECT

!  FIX PHI BY USING CONSTANCY OF SCATTER ANGLE

        COSSCAT_DN(N) = COSSCAT_DN(N+1)
        IF (STHETA_BOA.EQ.0.0D0 ) THEN
          PHI_ALL(N)     = PHI_ALL(N+1)
        ELSE
         CPHI = (COSSCAT_DN(N)-CALPHA*CTHETA)/STHETA/SALPHA
         IF ( CPHI.GT.1.0D0) CPHI = 1.0D0
         IF ( CPHI.LT.-1.0D0) CPHI = -1.0D0
         PHI_ALL(N)     = DACOS(CPHI)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 !  BUG CORRECTED 01 OCTOBER 2010------------------------- RTS, R. SPURR

!  IF AZIMUTH > 180 DEGREES. MUST ENSURE CONSISTENCY, ADD THIS LINE.

         IF ( PHI_BOA.GT.180.0D0 ) PHI_ALL(N) = PI2 - PHI_ALL(N)

 !  B G CORRECTED 01 OCTOBER 2010------------------------- RTS, R. SPURR
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        ENDIF

!  SUN PATHS, DIRECT SUN AT LAYER TOP
!  ==================================

!   MEANS THAT THE SZA AT LAYER TOP IS < 90.
!    ===> SZA < 90 FOR ALL FINE POINTS IN LAYER BENEATH

        IF ( DIRECT_SUN ) THEN

!  WORK UP FROM LEVEL N TO TOA
!    LAYER TOP CALCULATION GETS LEFT OUT AT TOA

         IF ( N .GT. 0 ) THEN
          STH0 = STHETA
          TH0  = THETA_ALL(N)
          DO K = N, 1, -1
           STH1 = STH0*RADII(K)/RADII(K-1)
           TH1  = DASIN(STH1)
           KS1  = TH0-TH1
           SUNPATHS(N,K) = DSIN(KS1)*RADII(K)/STH1
           STH0 = STH1
           TH0  = TH1
          ENDDO
         ENDIF

! DBG         WRITE(*,*)'REGULAR',N1,(SUNPATHS(N1,K),K=N1,1,-1)

!  FINE GRID CALCULATION IS ALWAYS REQUIRED
!  ----------------------------------------

!   START AT THE GRID POINT ON LOS PATH AND WORK UPWARDS,
!     FIRST ACROSS PARTIAL LAYER TO UPPER BOUNDARY, THEN CONTINUE
!     UPWARDS TO TOA ON A WHOLE LAYER BASIS. SINE RULE.

         IF ( DO_FINE ) THEN
          DO J = 1, NFINE
           STH0 = STHETAF(J)
           TH0  = THETAF(J)
           STH1 = STH0*RADII_FINE(N1,J)/RADII(N)
           TH1  = DASIN(STH1)
           KS1  = TH0-TH1
           SUNPATHS_FINE(N1,N1,J) = DSIN(KS1)*RADII_FINE(N1,J)/STH1
           STH0 = STH1
           TH0  = TH1
           DO K = N, 1, -1
            STH1 = STH0*RADII(K)/RADII(K-1)
            TH1  = DASIN(STH1)
            KS1  = TH0-TH1
            SUNPATHS_FINE(N1,K,J) = DSIN(KS1)*RADII(K)/STH1
            STH0 = STH1
            TH0  = TH1
           ENDDO
!  DBG           WRITE(*,*)'REGULAR',N1,J,(SUNPATHS_FINE(N1,K,J),K=N1,1,
          ENDDO

!  DBG         WRITE(*,*)'REGULAR',N,(SUNPATHS(N,K),K=N,1,-1)

         ENDIF

!  COMPLETE DIRECT SUN COMPUTATIONS

        ENDIF

!  SUN PATHS, NOT DIRECT SUN , WITH TANGENT POINT
!  ==============================================

!  ALTHOUGH LAYER TOP HAS A TANGENT POINT, NOT ALL OF THE FINE-GRID
!   POINTS WILL HAVE A TANGENT POINT.

        IF (.NOT.DIRECT_SUN ) THEN

!  FIRST DO THE LAYER-TOP CALCULATION.
!  -----------------------------------

!  TANGR = TANGENT POINT RADIUS.

         TANGR = STHETA*RADII(N)

!  NTRAVERSE(N) IS THE NUMBER OF LAYERS TRAVERSED BY RAY.

         KRAD = NLAYERS
         DO WHILE (TANGR.GT.RADII(KRAD))
          KRAD = KRAD - 1
         ENDDO
         NTRAVERSE(N) = KRAD + 1

!  START AT THE TOA ANGLES

         STH0 = TANGR/RADII(0)
         TH0 = DASIN(STH0)

!  WORK DOWNWARDS FROM TOA (SINE RULE) TO LEVEL IMMEDIATELY ABOVE
!  THE TANGENT LAYER. DON'T FORGET TO DOUBLE THE PATH LENGTH FOR
!  ANY LAYERS WHICH ARE TRAVERSED TWICE.

         DO K = 1, KRAD
          STH1 = RADII(K-1)*STH0/RADII(K)
          TH1 = DASIN(STH1)
          KS1 = TH1-TH0
          FAC = 1.0D0
          IF ( K.GT.N) FAC = 2.0D0
          SUNPATHS(N,K) = FAC*DSIN(KS1)*RADII(K-1)/STH1
          STH0 = STH1
          TH0  = TH1
         ENDDO

!  TANGENT LAYER PATH LENGTH. TWICE AGAIN. THE FOLLOWING CHECK IS GOOD.
!  CHECK       WRITE(*,*)TANGR/DTAN(TH1),RADII(KRAD)*DCOS(TH1)

         SUNPATHS(N,KRAD+1)=2.0D0*RADII(KRAD)*DCOS(TH1)

! DBG         WRITE(*,*)'--TANGENT',N1,(SUNPATHS(N1,K),K=NTRAVERSE(N1),1

!  FINE LAYER DEVELOPMENT (SLIGHTLY DIFFERENT)
!  ---------------------

         IF ( DO_FINE ) THEN
          DO J = 1, NFINE

!  IF THERE IS NO TANGENT POINT, REPEAT CALCULATION ABOVE.

           IF ( DIRECT_SUNF(J) ) THEN
            STH0 = STHETAF(J)
            TH0  = THETAF(J)
            STH1 = STH0*RADII_FINE(N1,J)/RADII(N)
            TH1  = DASIN(STH1)
            KS1  = TH0-TH1
            SUNPATHS_FINE(N1,N1,J) = DSIN(KS1)*RADII_FINE(N1,J)/STH1
            STH0 = STH1
            TH0  = TH1
            DO K = N, 1, -1
             STH1 = STH0*RADII(K)/RADII(K-1)
             TH1  = DASIN(STH1)
             KS1  = TH0-TH1
             SUNPATHS_FINE(N1,K,J) = DSIN(KS1)*RADII(K)/STH1
             STH0 = STH1
             TH0  = TH1
            ENDDO

! DBG           WRITE(*,*)'REGULAR',N1,J,
!     &       (SUNPATHS_FINE(N1,K,J),K=NTRAVERSE_FINE(N1,J),1,-1)

!  FINE GRID CALCULATION WITH TANGENT POINT
!  ----------------------------------------

           ELSE

!  LOCAL TANGENT RADIUS AND NUMBER OF LAYERS TRAVERSED

            TANGR = STHETAF(J)*RADII_FINE(N1,J)
            KRAD = NLAYERS
            DO WHILE (TANGR.GT.RADII(KRAD))
             KRAD = KRAD - 1
            ENDDO
            NTRAVERSE_FINE(N1,J) = KRAD + 1

!  START AGAIN AT TOA

            STH0 = TANGR/RADII(0)
            TH0  = DASIN(STH0)

!  WORK DOWN TO THE LEVEL N

            DO K = 1, N
             STH1 = RADII(K-1)*STH0/RADII(K)
             TH1 = DASIN(STH1)
             KS1 = TH1-TH0
             SUNPATHS_FINE(N1,K,J) = DSIN(KS1)*RADII(K-1)/STH1
             STH0 = STH1
             TH0  = TH1
            ENDDO

!  IN LAYER BELOW LEVEL N, WORK DOWN TO LEVEL OF FINE GRID RADIUS
!     (SINGLE CONTRIBUTION)

            STH1 = RADII(N)*STH0/RADII_FINE(N1,J)
            TH1 = DASIN(STH1)
            KS1 = TH1-TH0
            PATH = DSIN(KS1)*RADII(N)/STH1
            STH0 = STH1
            TH0  = TH1

!  IN LAYER BELOW LEVEL N, COMPLETE THE PATH DOWN TO THE TANGENT POINT
!    (DOUBLE CONTRIBUTION). FINISH.

            IF ( KRAD.EQ.N ) THEN

              PATH = PATH + 2.0D0*RADII_FINE(N1,J)*DCOS(TH1)
              SUNPATHS_FINE(N1,KRAD+1,J) = PATH

!  CHECK    WRITE(*,*)'1',TANGR/DTAN(TH1),RADII_FINE(N1,J)*DCOS(TH1)

!  IN LAYER BELOW LEVEL N, NEED TO GO DOWN FURTHER.
!    (FROM NOW ON, WE NEED DOUBLE CONTRIBUTIONS)
!     --- CONTINUE THE PATH TO THE BOTTOM OF THIS LAYER
!     --- CONTINUE DOWN BY WHOLE-LAYER STEPS UNTIL REACH LEVEL ABOVE TAN
!     --- COMPLETE THE PATH DOWN TO THE TANGENT POINT

            ELSE
              STH1 = RADII_FINE(N1,J)*STH0/RADII(N1)
              TH1 = DASIN(STH1)
              KS1 = TH1-TH0
              PATH = PATH + 2.0D0 * DSIN(KS1)*RADII_FINE(N1,J)/STH1
              SUNPATHS_FINE(N1,N1,J) = PATH
              STH0 = STH1
              TH0  = TH1
              DO K = N1 + 1, KRAD
               STH1 = RADII(K-1)*STH0/RADII(K)
               TH1 = DASIN(STH1)
               KS1 = TH1-TH0
               SUNPATHS_FINE(N1,K,J) = DSIN(KS1)*RADII(K-1)/STH1
               STH0 = STH1
               TH0  = TH1
              ENDDO
              SUNPATHS_FINE(N1,KRAD+1,J)=2.0D0*RADII(KRAD)*DCOS(TH1)
!  CHECK 2    WRITE(*,*)'2',TANGR/DTAN(TH1),RADII(KRAD)*DCOS(TH1)
            ENDIF

! DBG           WRITE(*,*)'TANGENT',N1,J,
!     &       (SUNPATHS_FINE(N1,K,J),K=NTRAVERSE_FINE(N1,J),1,-1)

!  COMPLETE TANGENT POINT CLAUSE

           ENDIF

!  COMPLETE FINE-GRID LOOP

          ENDDO

!  DBG        WRITE(*,*)'TANGENT',N,(SUNPATHS(N,K),K=NTRAVERSE(N),1,-1)

!  COMPLETE TANGENT POINT CALCULATION

         ENDIF

        ENDIF

!  END LAYER LOOP

      ENDDO

!  PARTIALS SECTION
!  ----------------

!  PARTIAL ANGLES AND LOSPATHS

      IF ( DO_PARTIALS ) THEN
        DO UT = 1, N_PARTIALS
          NP = PARTIALS_IDX(UT)

!  ESTABLISH THE LOS ANGLE, AND LOSPATHS

          TH0 = ALPHA_ALL(NP)
          STH0 = DSIN(TH0)
          SALPHA = RADII(NP)*STH0 / RADII_P(UT)
          ALPHA_P(UT) = DASIN(SALPHA)
          CALPHA = DSQRT(1.0D0-SALPHA*SALPHA)

!          KSI = ALPHA_ALL(NP) - ALPHA_P(UT)
!          LOSPATHS_P_UP(UT) = DSIN(KSI)*RADII(NP)/SALPHA

          KSI = ALPHA_P(UT) - ALPHA_ALL(NP-1)
          LOSPATHS_P_DN(UT) = DSIN(KSI)*RADII(NP-1)/SALPHA

!  ESTABLISH THE SOLAR ANGLE, BOTH NAIDR =SUN-AT-BOA AND GENERALLY

          XICUM = ALPHA_ALL(NLAYERS) - ALPHA_P(UT)
          DIRECT_SUN = .TRUE.
          IF (STHETA_BOA.EQ.0.0D0 ) THEN
            THETA_P = XICUM
            CTHETA = DCOS(THETA_P)
            STHETA = DSQRT(1.0D0-CTHETA*CTHETA)
          ELSE
            PX = - RADII_P(UT) * DSIN(XICUM)
            PY = 0.0D0
            PZ =   RADII_P(UT) * DCOS(XICUM)
            B = EX*PX + EY*PY + EZ*PZ
            CTHETA = -B/RADII_P(UT)
            DIRECT_SUN = (DIRECT_SUN.AND.CTHETA.GE.0.D0)
            STHETA = DSQRT(1.0D0-CTHETA*CTHETA)
            THETA_P = DACOS(CTHETA)
          ENDIF

!         WRITE(*,*)DIRECT_SUN,THETA_ALL(NP-1),THETA_P,THETA_ALL(NP)

!  SUN PATHS, DIRECT SUN AT LAYER TOP
!   FIRST CALCULATE SUNPATH FOR LAYER CONTAINING PARTIAL POINT
!   THEN WORK UPWARDS TO TOA USING THE SINE RULE.

          IF ( DIRECT_SUN ) THEN
            STH0 = STHETA
            TH0  = THETA_P
            STH1 = STH0*RADII_P(UT)/RADII(NP-1)
            TH1  = DASIN(STH1)
            KS1  = TH0-TH1
            SUNPATHS_P(UT,NP) = DSIN(KS1)*RADII_P(UT)/STH1
            STH0 = STH1
            TH0  = TH1
            DO K = NP-1, 1, -1
              STH1 = STH0*RADII(K)/RADII(K-1)
              TH1  = DASIN(STH1)
              KS1  = TH0-TH1
              SUNPATHS_P(UT,K) = DSIN(KS1)*RADII(K)/STH1
              STH0 = STH1
              TH0  = TH1
            ENDDO
          ENDIF

!  DEBUG
! 14       FORMAT(A11,1X,10F10.6)
!          J = PARTIALS_FINEIDX(UT)
!          WRITE(*,*)UT,NTRAVERSE_P(UT)
!          WRITE(*,*)UT,J,NP,LOSPATHS_P_UP(UT)/LOSPATHS(NP)
!           WRITE(*,14)'REGULAR FEB',(SUNPATHS(NP,K),K=NP,1,-1)
!           DO J = 1, NFINE
!             WRITE(*,14)'REGULAR FEB',(SUNPATHS_FINE(NP,K,J),K=NP,1,-1)
!             IF (J.EQ.PARTIALS_FINEIDX(UT))
!     &        WRITE(*,14)'REGULAR UT ',(SUNPATHS_P(UT,K),K=NP,1,-1)
!           ENDDO
!           WRITE(*,14)'REGULAR FEB',0.0D0,(SUNPATHS(NP-1,K),K=NP-1,1,-1
!
!  SUN PATHS, NOT DIRECT SUN , WITH TANGENT POINT. CODE SIMILAR TO ABOVE
!  TANGR = TANGENT POINT RADIUS FOR THIS RAY
!   NTRAVERSE_P(UT) = NUMBER OF LAYERS TRAVERSED BY RAY.

!  WORK DOWNWARDS FROM TOA (SINE RULE) TO LEVEL IMMEDIATELY ABOVE
!  THE TANGENT LAYER. DON'T FORGET TO DOUBLE THE PATH LENGTH FOR
!  ANY LAYERS WHICH ARE TRAVERSED TWICE.

!  TANGENT LAYER PATH LENGTH. TWICE AGAIN. THE FOLLOWING CHECK IS GOOD.
!  CHECK       WRITE(*,*)TANGR/DTAN(TH1),RADII(KRAD)*DCOS(TH1)

          IF (.NOT.DIRECT_SUN ) THEN
            TANGR = STHETA*RADII(N)
            KRAD = NLAYERS
            DO WHILE (TANGR.GT.RADII(KRAD))
              KRAD = KRAD - 1
            ENDDO
            NTRAVERSE_P(UT) = KRAD + 1
            STH0 = TANGR/RADII(0)
            TH0 = DASIN(STH0)
            DO K = 1, KRAD
              STH1 = RADII(K-1)*STH0/RADII(K)
              TH1 = DASIN(STH1)
              KS1 = TH1-TH0
              FAC = 1.0D0
              IF ( K.GT.NP) FAC = 2.0D0
              SUNPATHS_P(UT,K) = FAC*DSIN(KS1)*RADII(K-1)/STH1
              STH0 = STH1
              TH0  = TH1
            ENDDO
            SUNPATHS_P(UT,KRAD+1)=2.0D0*RADII_P(UT)*DCOS(TH1)
          ENDIF

!  FINISH PARTIALS LOOP

        ENDDO
      ENDIF

!  FINISH

      RETURN
      END SUBROUTINE OUTGOING_SPHERGEOM_FINE_DN

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

      SUBROUTINE MULTI_OUTGOING_ADJUSTGEOM &
         ( MAX_VZA, MAX_SZA, MAX_AZM, N_VZA, N_SZA, N_AZM, &
           HSURFACE, ERADIUS, DO_ADJUST_SURFACE, &
           ALPHA_BOA, THETA_BOA, PHI_BOA, &
           ALPHA_SSA, THETA_SSA, PHI_SSA, FAIL, MESSAGE, TRACE )

      IMPLICIT NONE

! STAND-ALONE GEOMETRY ROUTINE FOR ADJUSTING THE OUTGOING CORRECTION
!    STARTING INPUTS ARE THE BOA VALUES OF SZA, VZA AND PHI
!    NEED ALSO THE HEIGHT OF NEW SURFACE, EARTH RADIUS.

!  HEIGHT GRID HERE IS ARTIFICIAL

!  INPUTS

      INTEGER, INTENT(IN) ::           MAX_VZA, MAX_SZA, MAX_AZM
      INTEGER, INTENT(IN) ::           N_VZA, N_SZA, N_AZM
      DOUBLE PRECISION, INTENT(IN) ::  ERADIUS, HSURFACE
      LOGICAL, INTENT(IN) ::           DO_ADJUST_SURFACE
      DOUBLE PRECISION, INTENT(IN) ::  ALPHA_BOA (MAX_VZA)
      DOUBLE PRECISION, INTENT(IN) ::  THETA_BOA (MAX_SZA)

!  INPUT/OUTPUT

      DOUBLE PRECISION, INTENT(INOUT) ::  PHI_BOA   (MAX_AZM)

!  OUTPUTS

      DOUBLE PRECISION, INTENT(OUT) ::  ALPHA_SSA (MAX_VZA)
      DOUBLE PRECISION, INTENT(OUT) ::  THETA_SSA (MAX_VZA,MAX_SZA,MAX_AZM)
      DOUBLE PRECISION, INTENT(OUT) ::  PHI_SSA   (MAX_VZA,MAX_SZA,MAX_AZM)
      LOGICAL, INTENT(OUT) ::           FAIL
      CHARACTER (LEN=*), INTENT(INOUT) :: MESSAGE, TRACE

!  LOCAL

      INTEGER ::           J, I, K
      DOUBLE PRECISION ::  DEG_TO_RAD, EX, EY, EZ, PX, PY, PZ, PIE, PI2
      DOUBLE PRECISION ::  SALPHA_BOA, CALPHA_BOA, SPHI_BOA
      DOUBLE PRECISION ::  STHETA_BOA, CTHETA_BOA, CPHI_BOA
      DOUBLE PRECISION ::  KSI, CKSI, SKSI, XICUM, COSSCAT_UP
      DOUBLE PRECISION ::  PHI_ALL, ALPHA_ALL, THETA_ALL
      DOUBLE PRECISION ::  CTHETA, STHETA, CALPHA, SALPHA, CPHI
      DOUBLE PRECISION ::  B,RSSA
      CHARACTER (LEN=2) :: C2

!  INITIALISE OUTPUT

      FAIL    = .FALSE.
      MESSAGE = ' '
      TRACE   = ' '

!  CHECK RANGE OF INPUTS

      DO J = 1, N_VZA
       IF ( ALPHA_BOA(J).GE.90.0D0.OR.ALPHA_BOA(J).LT.0.0D0 ) THEN
        WRITE(C2,'(I2)')J
        MESSAGE = 'BOA LOS ANGLE OUTSIDE RANGE [0,90])'
        TRACE   = 'CHANGE BOA LOS ANGLE, NUMBER '//C2
        FAIL    = .TRUE.
        RETURN
       ENDIF
      ENDDO

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  BUG CORRECTED 01 OCTOBER 2010------------------------- RTS, R. SPURR

!  IMPORTANT, DO NOT LIMIT THE AZIMUTH TO < 180 DEGREES

      DO K = 1, N_AZM
        IF ( PHI_BOA(K).LT.0.0D0   ) PHI_BOA(K) = - PHI_BOA(K)
!        IF ( PHI_BOA(K).GT.180.0D0 ) PHI_BOA(K) = 360.0D0 - PHI_BOA(K)
        IF ( PHI_BOA(K).GT.360.0D0 ) PHI_BOA(K) = PHI_BOA(K) - 360.0D0
      ENDDO

!  BUG CORRECTED 01 OCTOBER 2010------------------------- RTS, R. SPURR
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      DO I = 1, N_SZA
       IF ( THETA_BOA(I).GE.90.0D0.OR.THETA_BOA(I).LT.0.0D0 ) THEN
        WRITE(C2,'(I2)')I
        MESSAGE = 'BOA SZA ANGLE OUTSIDE RANGE [0,90])'
        TRACE   = 'CHANGE BOA SZA ANGLE, NUMBER '//C2
        FAIL    = .TRUE.
        RETURN
       ENDIF
      ENDDO

!  NO ADJUSTMENT, JUST COPY AND EXIT

      IF ( .NOT. DO_ADJUST_SURFACE ) THEN
        DO J = 1, N_VZA
          ALPHA_SSA(J)   = ALPHA_BOA(J)
          DO I = 1, N_SZA
            DO K = 1, N_AZM
              THETA_SSA(J,I,K) = THETA_BOA(I)
              PHI_SSA(J,I,K)   = PHI_BOA(K)
            ENDDO
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  CONVERSION

      PIE = DACOS(-1.0D0)
      PI2 = 2.0D0 * PIE
      DEG_TO_RAD = PIE / 180.0D0

!  RADIUS OF SURFACE

      RSSA = HSURFACE + ERADIUS

!  START VZA LOOP

      DO J = 1, N_VZA

        ALPHA_ALL = ALPHA_BOA(J) * DEG_TO_RAD
        SALPHA_BOA = DSIN(ALPHA_ALL)
        CALPHA_BOA = DCOS(ALPHA_ALL)

!  SPECIAL CASE. DIRECT NADIR VIEWING. COMPUTE EVERYTHING AND EXIT.
!    (THIS IS THE SAME AS THE REGULAR PSEUDO-SPHERICAL )

        IF ( SALPHA_BOA.EQ.0.0D0 ) THEN
          ALPHA_SSA(J)   = ALPHA_BOA(J)
          DO I = 1, N_SZA
            DO K = 1, N_AZM
              THETA_SSA(J,I,K) = THETA_BOA(I)
              PHI_SSA(J,I,K)   = PHI_BOA(K)
            ENDDO
          ENDDO
          GO TO 567
        ENDIF

!  LOS ANGLE

        SALPHA       = ERADIUS * SALPHA_BOA / RSSA
        ALPHA_SSA(J) = DASIN(SALPHA)
        CALPHA       = DCOS(ALPHA_SSA(J))

!  LOSPATHS

        KSI  = ALPHA_ALL - ALPHA_SSA(J)
        SKSI = DSIN(KSI)
        CKSI = DCOS(KSI)
        XICUM = KSI

!  OUTPUT ANGLE IN DEGREES

        ALPHA_SSA(J) = ALPHA_SSA(J) / DEG_TO_RAD

!  LOS VECTOR

        PX = - RSSA * DSIN(XICUM)
        PY = 0.0D0
        PZ =   RSSA * DCOS(XICUM)

!  START SZA LOOP

        DO I = 1, N_SZA

          THETA_ALL = THETA_BOA(I) * DEG_TO_RAD
          STHETA_BOA = DSIN(THETA_ALL)
          CTHETA_BOA = DCOS(THETA_ALL)

!  START AZIMUTH LOOP

          DO K = 1, N_AZM

            PHI_ALL   = PHI_BOA(K)   * DEG_TO_RAD
            CPHI_BOA  = DCOS(PHI_ALL)
            SPHI_BOA  = DSIN(PHI_ALL)

!  DEFINE UNIT SOLAR VECTOR

            EX = - STHETA_BOA * CPHI_BOA
            EY = - STHETA_BOA * SPHI_BOA
            EZ = - CTHETA_BOA

!  SUN ANGLE

            B = EX*PX + EY*PY + EZ*PZ
            CTHETA = -B/RSSA
            STHETA = DSQRT(1.0D0-CTHETA*CTHETA)
            THETA_SSA(J,I,K) = DACOS(CTHETA)/DEG_TO_RAD
            IF ( CTHETA.LT.0.0D0 ) THEN
              WRITE(C2,'(I2)')J
              MESSAGE = 'LOS-PATH SZA ANGLE OUTSIDE RANGE [0,90])'
              TRACE   = 'CHECK INPUTS FOR LOS ANGLE '//C2
              FAIL    = .TRUE.
              RETURN
            ENDIF

!  SCATTERING ANGLE

            COSSCAT_UP  = - CALPHA_BOA * CTHETA_BOA + &
                            SALPHA_BOA * STHETA_BOA * CPHI_BOA

!  FIX PHI BY USING CONSTANCY OF SCATTER ANGLE

            IF ( PHI_BOA(K).EQ.180.0D0 ) THEN
              PHI_SSA(J,I,K) = PHI_ALL
            ELSE IF ( PHI_BOA(K) .EQ. 0.0D0 ) THEN
              PHI_SSA(J,I,K) = 0.0D0
            ELSE
              CPHI = (COSSCAT_UP+CALPHA*CTHETA)/STHETA/SALPHA
              IF ( CPHI.GT.1.0D0) CPHI = 1.0D0
              IF ( CPHI.LT.-1.0D0) CPHI = -1.0D0
              PHI_SSA(J,I,K) = DACOS(CPHI)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  BUG CORRECTED 05 OCTOBER 2010----------------- V. NATRAJ, R. SPURR
!  IF AZIMUTH > 180 DEGREES. MUST ENSURE CONSISTENCY, ADD THIS LINE.
              IF ( PHI_BOA(K).GT.180.0D0 ) THEN
                 PHI_SSA(J,I,K)= PI2 - PHI_SSA(J,I,K)
              ENDIF
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

            ENDIF
            PHI_SSA(J,I,K) = PHI_SSA(J,I,K) / DEG_TO_RAD

!  END AZIMUTH AND SOLAR LOOPS

          ENDDO
        ENDDO

!  CONTINUATION POINT

 567    CONTINUE

!  END LOS LOOP

      ENDDO

!  FINISH

      RETURN
      END SUBROUTINE MULTI_OUTGOING_ADJUSTGEOM

!

      SUBROUTINE VLIDORT_DBCORRECTION ( &
        FLUXMULT, &
        DO_SSCORR_OUTGOING, DO_REFRACTIVE_GEOMETRY, &
        DO_UPWELLING, NSTOKES, &
        NLAYERS, &
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
        SOLAR_BEAM_OPDEP, DO_REFLECTED_DIRECTBEAM, &
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
      DOUBLE PRECISION, INTENT (IN) ::   SLTERM_ISOTROPIC ( MAXSTOKES )
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
      DOUBLE PRECISION, INTENT (IN) ::  SOLAR_BEAM_OPDEP ( MAXBEAMS )
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
      INTEGER ::           UT, UTA, UM, UA, NC, IB, V, O1, O2, OM
      DOUBLE PRECISION ::  FINAL_SOURCE, TR, X0_FLUX, X0_BOA, ATTN, SUM

!  FIRST STAGE
!  -----------

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
                ATTN = SOLAR_BEAM_OPDEP(IB)
              ENDIF
              X0_FLUX = FOUR * X0_BOA
              ATTN    = ATTN * X0_FLUX
              ATTN_DB_SAVE(V) = ATTN

!  ASSIGN REFLECTION TERM
!  INITIALIZE CUMULATIVE SOURCE TERM = REFLEC * F.MU_0.T/PI
!    T = ATTENUATION OF DIRECT BEAM TO BOA, F = FLUX

              IF ( .NOT. DO_LAMBERTIAN_SURFACE ) THEN
                DO O1 = 1, NSTOKES
                  SUM = ZERO
                  DO O2 = 1, NSTOKES
                    OM = MUELLER_INDEX(O1,O2)
                    SUM = SUM + &
                    FLUXVEC(O2)*EXACTDB_BRDFUNC(OM,UM,UA,IB)
                  ENDDO
                  EXACTDB_SOURCE(V,O1) = ATTN * SUM
                ENDDO
              ELSE
                O1 = 1
                EXACTDB_SOURCE(V,O1) = ATTN * LAMBERTIAN_ALBEDO
              ENDIF

!  New Surface-Leaving stuff 17 May 2012

              IF ( DO_SURFACE_LEAVING ) THEN
                IF ( DO_SL_ISOTROPIC ) THEN
                   O1 = 1
                   EXACTDB_SOURCE(V,O1) = &
                      EXACTDB_SOURCE(V,O1) + SLTERM_ISOTROPIC(O1)
                ELSE
                  DO O1 = 1, NSTOKES
                    EXACTDB_SOURCE(V,O1) = &
                      EXACTDB_SOURCE(V,O1) + SLTERM_USERANGLES(O1,UM,UA,IB)
                  ENDDO
                ENDIF
              ENDIF

!  FINISH LOOPS OVER GEOMETRIES

            ENDDO
          ENDIF
        ENDDO
      ENDDO

!  UPWELLING RECURRENCE: TRANSMITTANCE OF EXACT SOURCE TERM
!  --------------------------------------------------------

!  INITIALIZE CUMULATIVE COUNT
!    MULTIPLY BY FLUX MULTIPLIER

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

!  END LAYER LOOP

        ENDDO

!  OFFGRID OUTPUT : PARTIAL LAYER TRANSMITTANCE, THEN SET RESULT

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)
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

!  ONGRID OUTPUT : SET FINAL CUMULATIVE SOURCE DIRECTLY

        ELSE

          DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
              DO UA = 1, N_USER_RELAZMS
                V = VZA_OFFSETS(IB,UM) + UA
                DO O1 = 1, NELEMENTS
                  FINAL_SOURCE = DB_CUMSOURCE(V,O1,NC)
                  STOKES_DB(UTA,V,O1) = FINAL_SOURCE
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

!  FINISH

      RETURN
      END SUBROUTINE VLIDORT_DBCORRECTION

      END MODULE vlidort_corrections

