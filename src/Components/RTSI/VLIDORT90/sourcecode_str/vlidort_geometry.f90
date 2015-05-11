! ###############################################################
! #                                                             #
! #                    THE VECTOR LIDORT MODEL                  #
! #                                                             #
! #  (Vector LInearized Discrete Ordinate Radiative Transfer)   #
! #   -      --         -        -        -         -           #
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
! #            VLIDORT_CHAPMAN  (called by VLIDORT_MASTER)      #
! #              BEAM_GEOMETRY_PREPARE                          #
! #                                                             #
! #            OUTGOING_SPHERGEOM_FINE_UP                       #
! #            OUTGOING_SPHERGEOM_FINE_DN                       #
! #            MULTI_OUTGOING_ADJUSTGEOM                        #
! #            OBSGEOM_OUTGOING_ADJUSTGEOM                      #
! #                                                             #
! #      Version 2.6. Modified geometry routine for when        #
! #                   using thermal emission at night           #
! #                                                             #
! #            LOSONLY_OUTGOING_ADJUSTGEOM                      #
! #                                                             #
! ###############################################################

      MODULE vlidort_geometry

      PRIVATE
      PUBLIC :: VLIDORT_CHAPMAN,&
                OUTGOING_SPHERGEOM_FINE_UP,&
                OUTGOING_SPHERGEOM_FINE_DN,&
                MULTI_OUTGOING_ADJUSTGEOM,&
                OBSGEOM_OUTGOING_ADJUSTGEOM,&
                LOSONLY_OUTGOING_ADJUSTGEOM

      CONTAINS

      SUBROUTINE VLIDORT_CHAPMAN ( &
        DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, &
        NLAYERS, N_SZANGLES, &
        SZANGLES, &
        EARTH_RADIUS, RFINDEX_PARAMETER, &
        HEIGHT_GRID, PRESSURE_GRID, &
        TEMPERATURE_GRID, FINEGRID, &
        SZA_LOCAL_INPUT, CHAPMAN_FACTORS, &
        FAIL, MESSAGE, TRACE )

!  This is the internal Chapman function calculation of the slant path
!  Chapman Factors required for slant-path optical thickness values.

!  This module calculates CHAPMAN_FACTORS internally inside VLIDORT,
!  saving you the job of doing it yourself, though you can input these
!  quantities, if the DO_CHAPMAN_FACTORS flag is not set.

!  The following options apply:

!   1. If the plane-parallel flag is on, no further inputs are required

!   2. If the plane-parallel flag is off, then Pseudo-spherical:
!       (a) Straight line geometry, must specify
!               Earth_radius, height grid
!       (b) Refractive geometry, must specify
!               Earth_radius, height grid
!               pressure grid, temperature grid

!  The logic will be checked before the module is called.

!  Newly programmed by R. Spurr, RT SOLUTIONS Inc. 5/5/05.

!    Based round a call to a pure geometry module which returns slant
!    path distances which was adapted for use in the Radiant model by
!    R. Spurr during an OCO L2 intensive April 24-29, 2005.

      USE VLIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::           DO_PLANE_PARALLEL
      LOGICAL, INTENT (IN) ::           DO_REFRACTIVE_GEOMETRY
      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           N_SZANGLES
      DOUBLE PRECISION, INTENT (IN) ::  SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) ::  EARTH_RADIUS
      DOUBLE PRECISION, INTENT (IN) ::  RFINDEX_PARAMETER
      DOUBLE PRECISION, INTENT (IN) ::  HEIGHT_GRID ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  PRESSURE_GRID ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  TEMPERATURE_GRID ( 0:MAXLAYERS )
      INTEGER, INTENT (IN) ::           FINEGRID ( MAXLAYERS )

      DOUBLE PRECISION, INTENT (OUT) :: SZA_LOCAL_INPUT &
          ( 0:MAXLAYERS, MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (OUT) :: CHAPMAN_FACTORS &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      LOGICAL, INTENT (OUT)           :: FAIL
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  Local variables
!  ---------------

!  number of iterations (refractive case only)
!      This is debug output

      INTEGER          ::  ITERSAVE(MAXLAYERS)

!  other local variables

      INTEGER          :: IB
      DOUBLE PRECISION :: SUN0

!  get spherical optical depths
!  ----------------------------

!  start beam loop

!mick - added initialization
      CHAPMAN_FACTORS = ZERO
      SZA_LOCAL_INPUT = ZERO

      DO IB = 1, N_SZANGLES

        SUN0 = SZANGLES(IB)

        CALL BEAM_GEOMETRY_PREPARE &
           ( MAX_SZANGLES, MAXLAYERS, IB, NLAYERS, FINEGRID, &
             SUN0, EARTH_RADIUS, RFINDEX_PARAMETER, &
             DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, &
             HEIGHT_GRID, PRESSURE_GRID, TEMPERATURE_GRID, &
             CHAPMAN_FACTORS, SZA_LOCAL_INPUT, &
             ITERSAVE, FAIL, MESSAGE )

!  return if failed

        IF ( FAIL ) THEN
          TRACE = 'Geometry failure in VLIDORT_CHAPMAN'
          RETURN
        ENDIF

!  end beam loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_CHAPMAN

!

      SUBROUTINE BEAM_GEOMETRY_PREPARE ( &
        MAXBEAMS, MAXLAYERS, IBEAM, NLAYERS, FINEGRID, &
        SZA_GEOM_TRUE, REARTH, RFINDEX_PARAMETER, &
        DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, &
        HEIGHTS, PRESSURES, TEMPERATURES, &
        CHAPMAN_FACTORS, SZA_LEVEL_OUTPUT, &
        ITERSAVE, FAIL, MESSAGE )

      USE VLIDORT_PARS, ONLY: ZERO, ONE

      IMPLICIT NONE

!  Generate path CHAPMAN_FACTORS and SZA angles SZA_LEVEL_OUTPUT
!  for a curved ray-traced beam through a multilayer atmosphere.

!  Coarse layering is input to the module. Values of Z, P, T  are
!  given at the layer boundaries, with the first value (index 0) at TOA.

!  The refractive geometry is assumed to start at the TOA level.

!  We also require the earth radius and the refractive index parameter
!   (For the Born-Wolf approximation)

!  There is no refraction if the flag DO_REFRACTIVE_GEOMETRY is not set.
!  In this case we do not require pressure and temperature information.
!  The calculation will then be for geometric rays.

!  The plane parallel Flag and the refractive geometry flag should not
!  both be true - this should be checked outside.

!  In the refracting case, fine-gridding of pressure and temperature is
!  done internally, temperature is interpolated linearly with height,
!  and pressure log-linearly. The refraction uses Snell's law rule.
!  Finelayer gridding assumes equidistant heights within coarse layers
!  but the number of fine layers can be varied

!  Output is specified at coarse layer boundaries

!  Module is stand-alone.

!  Reprogrammed for the OCO L2 algorithm
!   R. Spurr, RT Solutions, Inc.   April 27, 2005

!  Intended use in LIDORT and Radiant RT models.

!  Input arguments
!  ===============

!  input dimensioning

      INTEGER, INTENT (IN) ::          MAXLAYERS, MAXBEAMS

!  Beam index

      INTEGER, INTENT (IN) ::          IBEAM

!  number of coarse layers

      INTEGER, INTENT (IN) ::          NLAYERS

!  number of fine layers within coarse layers

      INTEGER, INTENT (IN) ::          FINEGRID(MAXLAYERS)

!  True solar zenith angle (degrees)

      DOUBLE PRECISION, INTENT (IN) :: SZA_GEOM_TRUE

!  Earth radius (km)

      DOUBLE PRECISION, INTENT (IN) :: REARTH

!  Refractive index parametaer (Born-Wolf approximation)

      DOUBLE PRECISION, INTENT (IN) :: RFINDEX_PARAMETER

!  flag for plane parallel case

      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL

!  flag for refractive geometry

      LOGICAL, INTENT (IN) ::          DO_REFRACTIVE_GEOMETRY

!  Coarse grids of heights, pressures and temperatures

      DOUBLE PRECISION, INTENT (IN) :: HEIGHTS (0:MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: PRESSURES (0:MAXLAYERS)
      DOUBLE PRECISION, INTENT (IN) :: TEMPERATURES (0:MAXLAYERS)

!  Output arguments
!  ================

!  Path segments distances (km)

!mick
!      DOUBLE PRECISION, INTENT (OUT) :: &
!        CHAPMAN_FACTORS (MAXLAYERS, MAXLAYERS, MAXBEAMS)
      DOUBLE PRECISION, INTENT (INOUT) :: &
        CHAPMAN_FACTORS (MAXLAYERS, MAXLAYERS, MAXBEAMS)

!  solar zenith angles at nadir

!mick
!      DOUBLE PRECISION, INTENT (OUT) :: &
!        SZA_LEVEL_OUTPUT (0:MAXLAYERS,MAXBEAMS)
      DOUBLE PRECISION, INTENT (INOUT) :: &
        SZA_LEVEL_OUTPUT (0:MAXLAYERS,MAXBEAMS)

!  number of iterations (refractive case only)
!   This is debug output

      INTEGER, INTENT (OUT)         :: ITERSAVE(MAXLAYERS)

!  output status

      LOGICAL, INTENT (OUT)           :: FAIL
      CHARACTER (LEN=*), INTENT (OUT) :: MESSAGE

!  Local variables
!  ===============

!  local dimensioning

      INTEGER, PARAMETER :: &
        LOCAL_MAXLAYERS = 101, &
        LOCAL_MAXFINELAYERS = 20

!  fine layer gridding for refraction

      DOUBLE PRECISION :: ZRFINE(LOCAL_MAXLAYERS,0:LOCAL_MAXFINELAYERS)
      DOUBLE PRECISION :: PRFINE(LOCAL_MAXLAYERS,0:LOCAL_MAXFINELAYERS)
      DOUBLE PRECISION :: TRFINE(LOCAL_MAXLAYERS,0:LOCAL_MAXFINELAYERS)

!  local height arrays

      DOUBLE PRECISION :: H(0:LOCAL_MAXLAYERS)
      DOUBLE PRECISION :: DELZ(LOCAL_MAXLAYERS)

!  help variables

      INTEGER          :: N, J, NRFINE, K, ITER, MAXF,IB
      LOGICAL          :: LOOP
      DOUBLE PRECISION :: GM_TOA, TH_TOA, MU_TOA, MU_NEXT
      DOUBLE PRECISION :: Z1, Z0, Z, T1, T0, T, P1, P0, Q1, Q0, Q
      DOUBLE PRECISION :: FU, FL, DEG_TO_RAD

      DOUBLE PRECISION :: LAYER_DIST, &
            MU_PREV, STH1, SINTH1, STH2, SINTH2, LOCAL_SUBTHICK, &
            PHI, PHI_0, PHI_CUM, SINPHI, DELPHI, REFRAC, RATIO, &
            RE_LOWER, RE_UPPER, DIST, STH2D, SINTH2D, SNELL

!  Standard temperature (K) and pressure (mbar).

      DOUBLE PRECISION, PARAMETER :: &
        T_STANDARD = 273.16D0, &
        P_STANDARD = 1013.25D0, &
        STP_RATIO  = T_STANDARD / P_STANDARD

!  Loschmidt's number (particles/cm2/km).

      DOUBLE PRECISION, PARAMETER :: &
        RHO_STANDARD = 2.68675D+24

!  Some setup operations
!  =====================

!  initialise output

      IB = IBEAM
      SZA_LEVEL_OUTPUT(0,IB) = 0.0D0
      DO N = 1, NLAYERS
        SZA_LEVEL_OUTPUT(N,IB) = 0.0D0
        DO K = 1, NLAYERS
          CHAPMAN_FACTORS(N,K,IB) = 0.0D0
        ENDDO
      ENDDO
      FAIL    = .FALSE.
      MESSAGE = ' '

!  check local dimensioning

      IF ( LOCAL_MAXLAYERS .LT. NLAYERS ) THEN
        MESSAGE = 'local coarse layer dimensioning insufficient'
        FAIL = .TRUE.
        RETURN
      ENDIF

!  earth radii and heights differences

      DO N = 0, NLAYERS
        H(N) = HEIGHTS(N) + REARTH
      ENDDO

      DO N = 1, NLAYERS
        DELZ(N) = HEIGHTS(N-1)-HEIGHTS(N)
      ENDDO

!  TOA values

      SZA_LEVEL_OUTPUT(0,IB) = SZA_GEOM_TRUE
      DEG_TO_RAD = DATAN(1.0D0) / 45.0D0
      TH_TOA = SZA_GEOM_TRUE * DEG_TO_RAD
      MU_TOA = DCOS(TH_TOA)
      GM_TOA = DSQRT ( 1.0D0 - MU_TOA * MU_TOA )
      STH2D  = 0.0D0

!  derive the fine values

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN

!mick fix - moved inside if block from above
        !  Check fine layers do not exceed local dimensions assigned
        MAXF = 0
        DO N = 1, NLAYERS
          MAXF = MAX(MAXF,FINEGRID(N))
        ENDDO
        IF ( LOCAL_MAXFINELAYERS .LT. MAXF ) THEN
          MESSAGE = 'local fine layer dimensioning insufficient'
          FAIL = .TRUE.
          RETURN
        ENDIF

        Z0 = HEIGHTS(0)
        P0 = PRESSURES(0)
        T0 = TEMPERATURES(0)
        Q0 = DLOG(P0)
        DO N = 1, NLAYERS
          NRFINE = FINEGRID(N)
          LOCAL_SUBTHICK = DELZ(N) / DBLE(NRFINE)
          P1 = PRESSURES(N)
          Z1 = HEIGHTS(N)
          T1 = TEMPERATURES(N)
          Q1 = DLOG(P1)
          ZRFINE(N,0) = Z0
          PRFINE(N,0) = P0
          TRFINE(N,0) = T0
          DO J = 1, NRFINE - 1
            Z  = Z0 - DBLE(J)*LOCAL_SUBTHICK
            FL = ( Z0 - Z ) / DELZ(N)
            FU = 1.0d0 - FL
            Q  = FL * Q1 + FU * Q0
            T  = FL * T0 + FU * T1
            PRFINE(N,J) = DEXP (  Q )
            TRFINE(N,J) = T
            ZRFINE(N,J) = Z
          ENDDO
          PRFINE(N,NRFINE) = P1
          TRFINE(N,NRFINE) = T1
          ZRFINE(N,NRFINE) = Z1
!              write(*,'(i3,11F10.4)')N,(PRFINE(N,J),J=0,NRFINE)
          Z0 = Z1
          P0 = P1
          T0 = T1
          Q0 = Q1
        ENDDO
      ENDIF

!  plane-parallel case
!  ===================

      IF ( DO_PLANE_PARALLEL ) THEN
        DO N = 1, NLAYERS
          SZA_LEVEL_OUTPUT(N,IB) = SZA_GEOM_TRUE
          DO K = 1, N
            CHAPMAN_FACTORS(N,K,IB) = 1.0D0 / MU_TOA
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Refractive Geometry case
!  ========================

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN

!  Rob fix 9/9/14. Sun zero case

        IF ( GM_TOA .eq. zero ) then
          DO N = 1, NLAYERS
            CHAPMAN_FACTORS(N,1:N,IB) = one
          ENDDO
          RETURN
        ENDIF

!  Non-zero case....(resume code)

!  Start value of SZA cosine

        MU_PREV = MU_TOA

!  start layer loop

        DO N = 1, NLAYERS

!  start values

          SINTH1 = GM_TOA * H(N) / H(0)
          STH1 = DASIN(SINTH1)
          PHI_0 = TH_TOA - STH1
          NRFINE = FINEGRID(N)

!  iteration loop

          ITER = 0
          LOOP = .TRUE.
          DO WHILE (LOOP.AND.ITER.LT.100)
            ITER = ITER + 1
            PHI_CUM = 0.0D0
            RE_UPPER = ZRFINE(1,0) + REARTH
            RATIO  = PRFINE(1,0) * STP_RATIO / TRFINE(1,0)
            REFRAC = 1.0D0 + RFINDEX_PARAMETER * RATIO
            SNELL = REFRAC * RE_UPPER * SINTH1
            DO K = 1, N
              LAYER_DIST = 0.0D0
              LOCAL_SUBTHICK = DELZ(K) / DBLE(NRFINE)
              DO J = 0, NRFINE - 1
                RATIO  = PRFINE(K,J) * STP_RATIO / TRFINE(K,J)
                REFRAC = 1.0D0 + RFINDEX_PARAMETER * RATIO
                RE_LOWER = RE_UPPER - LOCAL_SUBTHICK
                SINTH2 = SNELL/ (REFRAC * RE_UPPER )
                IF ( SINTH2.GT.1.0D0 ) SINTH2 = 1.0D0
                STH2 = DASIN(SINTH2)
                SINTH2D = RE_UPPER * SINTH2 / RE_LOWER
                IF ( SINTH2D .GT. 1.0D0 ) THEN
                  MESSAGE = 'refraction yields angles > 90 some levels'
                  FAIL = .TRUE.
                  RETURN
                ENDIF
                STH2D = DASIN(SINTH2D)
                PHI = STH2D - STH2
                SINPHI = DSIN(PHI)
                PHI_CUM = PHI_CUM + PHI
                DIST = RE_UPPER * SINPHI / SINTH2D
                LAYER_DIST = LAYER_DIST +  DIST
                RE_UPPER = RE_LOWER
              ENDDO
              CHAPMAN_FACTORS(N,K,IB) = LAYER_DIST / DELZ(K)
            ENDDO

!  examine convergence

            DELPHI = PHI_0 - PHI_CUM
            LOOP = (DABS(DELPHI/PHI_CUM).GT.1.0D-4)

!  Fudge factors to speed up the iteration

            IF ( SZA_GEOM_TRUE .GT. 88.7D0 ) THEN
              STH1 = STH1 + 0.1 * DELPHI
              PHI_0 = TH_TOA - STH1
            ELSE IF ( SZA_GEOM_TRUE .LT. 80.0D0 ) THEN
              PHI_0 = PHI_CUM
              STH1 = TH_TOA - PHI_0
            ELSE
              STH1 = STH1 + 0.3 * DELPHI
              PHI_0 = TH_TOA - STH1
            ENDIF
            SINTH1 = DSIN(STH1)

          ENDDO

!  failure

          IF ( LOOP ) THEN
            MESSAGE = 'refractive iteration not converged'
            FAIL = .TRUE.
            RETURN
          ENDIF

!  Update and save angle output

          MU_NEXT = DCOS(STH2D)
          MU_PREV = MU_NEXT
          SZA_LEVEL_OUTPUT(N,IB) = DACOS(MU_NEXT) / DEG_TO_RAD
          ITERSAVE(N) = ITER

        ENDDO

!  Straight line geometry
!  ======================

      ELSE

!  Rob fix 9/9/14. Sun zero case

        IF ( GM_TOA .eq. zero ) then
          DO N = 1, NLAYERS
            CHAPMAN_FACTORS(N,1:N,IB) = one
          ENDDO
          RETURN
        ENDIF

!  Non-zero case....(resume code)

        DO N = 1, NLAYERS

!  start values

          SINTH1 = GM_TOA * H(N) / H(0)
          STH1   = DASIN(SINTH1)
          RE_UPPER = H(0)

!  solar zenith angles are all the same = input value

          SZA_LEVEL_OUTPUT(N,IB) = SZA_GEOM_TRUE

! loop over layers K from 1 to layer N

          DO K = 1, N

!  sine-rule; PHI = earth-centered angle

            RE_LOWER = RE_UPPER - DELZ(K)
            SINTH2 = RE_UPPER * SINTH1 / RE_LOWER
            STH2   = DASIN(SINTH2)
            PHI    = STH2 - STH1
            SINPHI = DSIN(PHI)
            DIST = RE_UPPER * SINPHI / SINTH2
            CHAPMAN_FACTORS(N,K,IB) = DIST / DELZ(K)

!  re-set

            RE_UPPER = RE_LOWER
            SINTH1 = SINTH2
            STH1   = STH2

          ENDDO

!  finish main layer loop

        ENDDO

!  Finish

      ENDIF

!  end of routine

      RETURN
      END SUBROUTINE BEAM_GEOMETRY_PREPARE

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

      SUBROUTINE OBSGEOM_OUTGOING_ADJUSTGEOM                      &
         ( MAX_VZA, MAX_SZA, MAX_AZM, N_VZA, N_SZA, N_AZM,        & ! INPUT
          HSURFACE, ERADIUS, DO_ADJUST_SURFACE,                   & ! INPUT
          ALPHA_BOA, THETA_BOA, PHI_BOA,                          & ! INPUT
          ALPHA_SSA, THETA_SSA, PHI_SSA, FAIL, MESSAGE, TRACE )     ! OUTPUT

! STAND-ALONE GEOMETRY ROUTINE FOR ADJUSTING THE OUTGOING CORRECTION
!    STARTING INPUTS ARE THE BOA VALUES OF SZA, VZA AND PHI
!    NEED ALSO THE HEIGHT OF NEW SURFACE, EARTH RADIUS.

!  HEIGHT GRID HERE IS ARTIFICIAL

      IMPLICIT NONE

!  INPUTS

      INTEGER          :: MAX_VZA, MAX_SZA, MAX_AZM
      INTEGER          :: N_VZA, N_SZA, N_AZM
      DOUBLE PRECISION :: ERADIUS, HSURFACE
      LOGICAL          :: DO_ADJUST_SURFACE
      DOUBLE PRECISION :: ALPHA_BOA (MAX_VZA)
      DOUBLE PRECISION :: THETA_BOA (MAX_SZA)
      DOUBLE PRECISION :: PHI_BOA   (MAX_AZM)

!  OUTPUTS

      DOUBLE PRECISION :: ALPHA_SSA (MAX_VZA)
      DOUBLE PRECISION :: THETA_SSA (MAX_VZA,MAX_SZA,MAX_AZM)
      DOUBLE PRECISION :: PHI_SSA   (MAX_VZA,MAX_SZA,MAX_AZM)
      LOGICAL          :: FAIL
      CHARACTER*(*)    :: MESSAGE, TRACE

!  LOCAL

      INTEGER          :: J, I, K
      DOUBLE PRECISION :: PIE, PI2, DEG_TO_RAD, EX, EY, EZ, PX, PY, PZ
      DOUBLE PRECISION :: SALPHA_BOA, CALPHA_BOA, SPHI_BOA
      DOUBLE PRECISION :: STHETA_BOA, CTHETA_BOA, CPHI_BOA
      DOUBLE PRECISION :: KSI, CKSI, SKSI, XICUM, COSSCAT_UP
      DOUBLE PRECISION :: PHI_ALL, ALPHA_ALL, THETA_ALL
      DOUBLE PRECISION :: CTHETA, STHETA, CALPHA, SALPHA, CPHI
      DOUBLE PRECISION :: B,RSSA
      CHARACTER*2      :: C2

!  INITIALISE OUTPUT

      FAIL = .FALSE.
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

      DO J = 1, N_AZM
        IF ( PHI_BOA(J).LT.0.0D0   ) PHI_BOA(J) = - PHI_BOA(J)
!        IF ( PHI_BOA(J).GT.180.0D0 ) PHI_BOA(J) = 360.0D0 - PHI_BOA(J)
        IF ( PHI_BOA(J).GT.360.0D0 ) PHI_BOA(J) = PHI_BOA(J) - 360.0D0
      ENDDO

!  BUG CORRECTED 01 OCTOBER 2010------------------------- RTS, R. SPURR
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!      DO J = 1, N_AZM
!        IF ( PHI_BOA(J).LT.0.0D0   ) PHI_BOA(J) = - PHI_BOA(J)
!        IF ( PHI_BOA(J).GT.180.0D0 ) PHI_BOA(J) = 360.0D0 - PHI_BOA(J)
!      ENDDO

      DO J = 1, N_SZA
       IF ( THETA_BOA(J).GE.90.0D0.OR.THETA_BOA(J).LT.0.0D0 ) THEN
        WRITE(C2,'(I2)')J
        MESSAGE = 'BOA SZA ANGLE OUTSIDE RANGE [0,90])'
        TRACE   = 'CHANGE BOA SZA ANGLE, NUMBER '//C2
        FAIL    = .TRUE. ; RETURN
       ENDIF
      ENDDO

!  NO ADJUSTMENT, JUST COPY AND EXIT

      IF ( .NOT. DO_ADJUST_SURFACE ) THEN
        DO J = 1, N_VZA
          ALPHA_SSA(J)     = ALPHA_BOA(J)
          THETA_SSA(1,J,1) = THETA_BOA(J)
          PHI_SSA(1,1,J)   = PHI_BOA(J)
        ENDDO
        RETURN
      ENDIF

!  CONVERSION

      PIE = DACOS(-1.0D0)
      PI2 = 2.0D0 * PIE
      DEG_TO_RAD = PIE / 180.0D0

!  RADIUS OF SURFACE

      RSSA = HSURFACE + ERADIUS

!  START GEOMETRY LOOP (VZA INDEX WILL BE THE DRIVER.....)

      DO J = 1, N_VZA

        ALPHA_ALL = ALPHA_BOA(J) * DEG_TO_RAD
        SALPHA_BOA = DSIN(ALPHA_ALL)
        CALPHA_BOA = DCOS(ALPHA_ALL)

!  SPECIAL CASE. DIRECT NADIR VIEWING. COMPUTE EVERYTHING AND EXIT.
!    (THIS IS THE SAME AS THE REGULAR PSEUDO-SPHERICAL )

        IF ( SALPHA_BOA.EQ.0.0D0 ) THEN
          ALPHA_SSA(J)     = ALPHA_BOA(J)
          THETA_SSA(1,J,1) = THETA_BOA(J)
          PHI_SSA(1,1,J)   = PHI_BOA(J)
        ENDIF

!  LOS ANGLE, GENERAL CASE

        IF ( SALPHA_BOA.NE.0.0D0 ) THEN

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

!  SZA CALCULATION

         THETA_ALL = THETA_BOA(J) * DEG_TO_RAD
         STHETA_BOA = DSIN(THETA_ALL)
         CTHETA_BOA = DCOS(THETA_ALL)

!  AZIMUTH CALCULATION

         PHI_ALL   = PHI_BOA(J)   * DEG_TO_RAD
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
         THETA_SSA(1,J,1) = DACOS(CTHETA)/DEG_TO_RAD
         IF ( CTHETA.LT.0.0D0 ) THEN
            WRITE(C2,'(I2)')J
            MESSAGE = 'LOS-PATH SZA ANGLE OUTSIDE RANGE [0,90])'
            TRACE   = 'CHECK INPUTS FOR LOS ANGLE '//C2
            FAIL    = .TRUE. ; RETURN
         ENDIF

!  SCATTERING ANGLE

         COSSCAT_UP  = - CALPHA_BOA * CTHETA_BOA + &
                         SALPHA_BOA * STHETA_BOA * CPHI_BOA 

!  FIX PHI BY USING CONSTANCY OF SCATTER ANGLE

         IF ( PHI_BOA(J).EQ.180.0D0 ) THEN
            PHI_SSA(1,1,J) = PHI_ALL
         ELSE IF ( PHI_BOA(J) .EQ. 0.0D0 ) THEN
            PHI_SSA(1,1,J) = 0.0D0
         ELSE
            CPHI = (COSSCAT_UP+CALPHA*CTHETA)/STHETA/SALPHA
            IF ( CPHI.GT.1.0D0) CPHI = 1.0D0
            IF ( CPHI.LT.-1.0D0) CPHI = -1.0D0
            PHI_SSA(1,1,J) = DACOS(CPHI)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  BUG CORRECTED 05 OCTOBER 2010----------------- V. NATRAJ, R. SPURR
!  IF AZIMUTH > 180 DEGREES. MUST ENSURE CONSISTENCY, ADD THIS LINE.
            IF ( PHI_BOA(J).GT.180.0D0 ) THEN
                 PHI_SSA(1,1,J)= PI2 - PHI_SSA(1,1,J)
            ENDIF
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         ENDIF
         PHI_SSA(1,1,J) = PHI_SSA(1,1,J) / DEG_TO_RAD

!  END OF REGULAR LOS CASE

        ENDIF

!  END LOS LOOP

      ENDDO

!  FINISH

      RETURN
      END SUBROUTINE OBSGEOM_OUTGOING_ADJUSTGEOM

!

      SUBROUTINE LOSONLY_OUTGOING_ADJUSTGEOM &
       ( MAX_VZA, N_VZA, HSURFACE, ERADIUS, DO_ADJUST_SURFACE, ALPHA_BOA, &
         ALPHA_SSA, FAIL, MESSAGE, TRACE )

      IMPLICIT NONE

!  THERMAL ONLY (No solar)

! STAND-ALONE GEOMETRY ROUTINE FOR ADJUSTING THE OUTGOING CORRECTION
!    STARTING INPUTS ARE THE BOA VALUES OF VZA
!    NEED ALSO THE HEIGHT OF NEW SURFACE, EARTH RADIUS.

!  HEIGHT GRID HERE IS ARTIFICIAL

!  INPUTS

      INTEGER, INTENT(IN) ::           MAX_VZA
      INTEGER, INTENT(IN) ::           N_VZA
      DOUBLE PRECISION, INTENT(IN) ::  ERADIUS, HSURFACE
      LOGICAL, INTENT(IN) ::           DO_ADJUST_SURFACE
      DOUBLE PRECISION, INTENT(IN) ::  ALPHA_BOA (MAX_VZA)

!  OUTPUTS

      DOUBLE PRECISION, INTENT(OUT) ::  ALPHA_SSA (MAX_VZA)
      LOGICAL, INTENT(OUT) ::           FAIL
      CHARACTER (LEN=*), INTENT(INOUT) :: MESSAGE, TRACE

!  LOCAL

      INTEGER ::           J
      DOUBLE PRECISION ::  DEG_TO_RAD, PIE, PI2
      DOUBLE PRECISION ::  SALPHA_BOA, CALPHA_BOA
      DOUBLE PRECISION ::  ALPHA_ALL,  CALPHA, SALPHA, RSSA
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

!  NO ADJUSTMENT, JUST COPY AND EXIT

      IF ( .NOT. DO_ADJUST_SURFACE ) THEN
        DO J = 1, N_VZA
          ALPHA_SSA(J)   = ALPHA_BOA(J)
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
        ELSE
           SALPHA       = ERADIUS * SALPHA_BOA / RSSA
           ALPHA_SSA(J) = DASIN(SALPHA)
           CALPHA       = DCOS(ALPHA_SSA(J))
           ALPHA_SSA(J) = ALPHA_SSA(J) / DEG_TO_RAD
        ENDIF

!  END LOS LOOP

      ENDDO

!  FINISH

      RETURN
      END SUBROUTINE LOSONLY_OUTGOING_ADJUSTGEOM

      END MODULE vlidort_geometry
