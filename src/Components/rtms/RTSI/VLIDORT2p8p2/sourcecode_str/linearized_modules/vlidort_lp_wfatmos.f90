
! ###############################################################
! #                                                             #
! #                       VLIDORT_2p8p2                         #
! #                                                             #
! #  Vectorized LInearized Discrete Ordinate Radiative Transfer #
! #  -          --         -        -        -         -        #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Authors :     Robert. J. D. Spurr (1)                      #
! #                Matt Christi                                 #
! #                                                             #
! #  Address (1) : RT Solutions, inc.                           #
! #                9 Channing Street                            #
! #                Cambridge, MA 02138, USA                     #
! #                                                             #
! #  Tel:          (617) 492 1183                               #
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  This Version :   VLIDORT_2p8p2                             #
! #  Release Date :   15 April 2020                             #
! #                                                             #
! #  Previous VLIDORT Versions under Standard GPL 3.0:          #
! #  ------------------------------------------------           #
! #                                                             #
! #      2.7   F90, released August 2014                        #
! #      2.8   F90, released May    2017                        #
! #      2.8.1 F90, released August 2019                        # 
! #                                                             #
! #  Features Summary of Recent VLIDORT Versions:               #
! #  -------------------------------------------                #
! #                                                             #
! #      NEW: TOTAL COLUMN JACOBIANS         (2.4)              #
! #      NEW: BPDF Land-surface KERNELS      (2.4R)             #
! #      NEW: Thermal Emission Treatment     (2.4RT)            #
! #      Consolidated BRDF treatment         (2.4RTC)           #
! #      f77/f90 Release                     (2.5)              #
! #      External SS / New I/O Structures    (2.6)              #
! #                                                             #
! #      SURFACE-LEAVING / BRDF-SCALING      (2.7)              #
! #      TAYLOR Series / OMP THREADSAFE      (2.7)              #
! #      New Water-Leaving Treatment         (2.8)              #
! #      LBBF & BRDF-Telescoping, enabled    (2.8)              #
! #      Several Performance Enhancements    (2.8)              #
! #      Water-leaving coupled code          (2.8.1)            #
! #      Planetary problem, media properties (2.8.1)            #
! #                                                             #
! #  Features Summary of This VLIDORT Version                   #
! #  ----------------------------------------                   #
! #                                                             #
! #   2.8.2, released 15 April 2020.                            #
! #     ==> Geometry (FO/MS), check/derive separation           #
! #     ==> New setup_master for Geometry/Check/Derive          #
! #     ==> Reduction of zeroing, some dynamic memory           #
! #     ==> Use of F-matrixes only in FO code                   #
! #     ==> Use I/O type structures directly                    #
! #     ==> Doublet geometry post-processing option             #
! #                                                             #
! ###############################################################

! ###################################################################
! #                                                                 #
! # This is Version 2.8.2 of the VLIDORT_2p8 software library.      #
! # This library comes with the Standard GNU General Public License,#
! # Version 3.0, 29 June 2007. Please read this license carefully.  #
! #                                                                 #
! #      VLIDORT Copyright (c) 2003-2020.                           #
! #          Robert Spurr, RT Solutions, Inc.                       #
! #          9 Channing Street, Cambridge, MA 02138, USA.           #
! #                                                                 #
! # This file is part of VLIDORT_2p8p2 ( Version 2.8.2 )            #
! #                                                                 #
! # VLIDORT_2p8p2 is free software: you can redistribute it         #
! # and/or modify it under the terms of the Standard GNU GPL        #
! # (General Public License) as published by the Free Software      #
! # Foundation, either version 3.0 of the License, or any           #
! # later version.                                                  #
! #                                                                 #
! # VLIDORT_2p8p2 is distributed in the hope that it will be        #
! # useful, but WITHOUT ANY WARRANTY; without even the implied      #
! # warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR         #
! # PURPOSE. See the Standard GNU General Public License (GPL)      #
! # for more details.                                               #
! #                                                                 #
! # You should have received a copy of the Standard GNU General     #
! # Public License (GPL) Version 3.0, along with the VLIDORT_2p8p2  #
! # code package. If not, see <http://www.gnu.org/licenses/>.       #
! #                                                                 #
! ###################################################################

! ##########################################################
! #                                                        #
! # Subroutines in this Module                             #
! #                                                        #
! #     Top level routines--------------                   #
! #                                                        #
! #       UPUSER_PROFILEWF                                 #
! #       DNUSER_PROFILEWF                                 #
! #                                                        #
! #       ------Calling                                    #
! #                                                        #
! #            GET_LP_TOASOURCE                            #
! #            GET_LP_BOASOURCE                            #
! #                                                        #
! #            LP_WHOLELAYER_STERM_UP                      #
! #            LP_WHOLELAYER_STERM_DN                      #
! #            LP_PARTLAYER_STERM_UP                       #
! #            LP_PARTLAYER_STERM_DN                       #
! #                                                        #
! #       MIFLUX_PROFILEWF                                 #
! #                                                        #
! #       ------Calling                                    #
! #                                                        #
! #            QUADPROFILEWF_LEVEL_UP                      #
! #            QUADPROFILEWF_LEVEL_DN                      #
! #            QUADPROFILEWF_OFFGRID_UP                    #
! #            QUADPROFILEWF_OFFGRID_DN                    #
! #                                                        #
! ##########################################################

!  4/15/20. Version 2.8.2. Separate module created for the LC_CONVERGE routines (removed from here)

      MODULE vlidort_lp_wfatmos_m

      PRIVATE :: GET_LP_TOASOURCE,        GET_LP_BOASOURCE,         & 
                 LP_WHOLELAYER_STERM_UP,  LP_WHOLELAYER_STERM_DN,   &
                 LP_PARTLAYER_STERM_UP,   LP_PARTLAYER_STERM_DN,    &
                 QUADPROFILEWF_LEVEL_UP  , QUADPROFILEWF_LEVEL_DN,  &
                 QUADPROFILEWF_OFFGRID_UP, QUADPROFILEWF_OFFGRID_DN

      PUBLIC  :: UPUSER_PROFILEWF,     &
                 DNUSER_PROFILEWF,     &
                 MIFLUX_PROFILEWF

      CONTAINS

      SUBROUTINE UPUSER_PROFILEWF ( &
        DO_USER_STREAMS,    DO_OBSERVATION_GEOMETRY, DO_INCLUDE_MVOUTPUT,       & ! Input flags (RT mode)
        DO_SOLAR_SOURCES,   DO_INCLUDE_THERMEMISS,   DO_THERMAL_TRANSONLY,      & ! Input flags (sources)
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE,   DO_MSMODE_VLIDORT,         & ! Input flags (Surface)
        DO_DBCORRECTION,    DO_INCLUDE_DIRECTRF,     DO_INCLUDE_DIRECTSL,       & ! Input flags (Beam/scattering)
        DO_LAYER_SCATTERING, FOURIER, IBEAM, NSTOKES, NSTREAMS, NLAYERS,        & ! Input numbers (basic)
        N_USER_STREAMS, N_USER_LEVELS, LAYER_TO_VARY, NV_PARAMETERS, LOCAL_UM_START, & ! Input numbers (basic)
        UTAU_LEVEL_MASK_UP, MUELLER_INDEX,                                   & ! Input bookkeeping + levels
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,        & ! Input partial-layer control
        FLUX_MULTIPLIER, SURFACE_FACTOR, QUAD_WEIGHTS, QUAD_STRMWTS,         & ! Input Flux and quadrature
        T_DELT_DISORDS,  T_DELT_USERM, T_UTUP_USERM, CUMSOURCE_UP,           & ! Input Transmittances
        ALBEDO, BRDF_F, USER_BRDF_F, DELTAU_SLANT,                           & ! Input Surface BRDF.
        RF_USER_DIRECT_BEAM, SL_USERTERM, LP_TRANS_ATMOS_FINAL,              & ! Input surface radiances (beam/sleave)
        K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, LCON, MCON,   & ! Input Homog. RTE Soln.
        T_WLOWER, UHOM_UPDN, UHOM_UPUP, UPAR_UP_1, UPAR_UP_2,                & ! Input thermal and User solutions
        HMULT_1, HMULT_2, EMULT_UP, UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP,   & ! Input multiplier  
        L_DELTAU_VERT, L_T_DELT_USERM, L_T_UTUP_USERM,  L_T_DELT_DISORDS,    & ! Linearized delta and transmittances
        L_T_DELT_EIGEN, L_SOLA_XPOS, L_SOLB_XNEG, L_WLOWER, NCON, PCON,      & ! Linearized Discr.Ord. solutions.
        L_T_WLOWER, L_LAYER_TSUP_UP, L_LAYER_TSUP_UTUP,                      & ! Linearized thermal solutions
        L_UHOM_UPDN, L_UHOM_UPUP, L_UPAR_UP_1, LP_UPAR_UP_2,                 & ! Linearized homog. and beam solutions
        L_HMULT_1, L_HMULT_2, LP_EMULT_UP,                                   & ! Linearized Homog. and Beam multipliers
        L_UT_HMULT_UU, L_UT_HMULT_UD, LP_UT_EMULT_UP,                        & ! Linearized Homog. and Beam multipliers  
        L_BOA_THTONLY_SOURCE, PROFILEWF_F )                                    ! OUTPUT

!   Streamlined for Version 2.8. 7/6/16
!  4/15/20. Version 2.8.2. BRDF arrays defined locally, each Fourier

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXSTREAMS, MAXBEAMS, MAXSTOKES, MAXLAYERS, MAX_PARTLAYERS, MAX_ATMOSWFS, &
                                 MAX_USER_STREAMS, MAX_USER_LEVELS, MAX_USER_VZANGLES, MAX_SZANGLES, MAX_DIRECTIONS,   &
                                 MAXSTOKES_SQ, MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, ZERO, UPIDX

      IMPLICIT NONE

!  Inputs
!  ======

!  Input flags

      LOGICAL, INTENT (IN) ::           DO_USER_STREAMS
      LOGICAL, INTENT (IN) ::           DO_OBSERVATION_GEOMETRY
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_MVOUTPUT

      LOGICAL, INTENT (IN) ::           DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::           DO_THERMAL_TRANSONLY

      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::           DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT (IN) ::           DO_DBCORRECTION

!  4/9/19. Replacement of DIRECTBEAM by two separate flags
      
      LOGICAL, intent(in)  ::           DO_INCLUDE_DIRECTRF
      LOGICAL, intent(in)  ::           DO_INCLUDE_DIRECTSL
!      LOGICAL, INTENT (IN) ::           DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT (IN) ::           DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )
      LOGICAL, INTENT (IN) ::           DO_MSMODE_VLIDORT

!  removed, Version 2.8, 7.17.16
!      LOGICAL, INTENT (IN) ::           DO_QUAD_OUTPUT

!  numbers

      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           IBEAM
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      INTEGER, INTENT (IN) ::           LAYER_TO_VARY, NV_PARAMETERS

!  Level output control + bookkeeping

      INTEGER, INTENT (IN) ::           LOCAL_UM_START
      INTEGER, INTENT (IN) ::           N_USER_LEVELS
      INTEGER, INTENT (IN) ::           UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )

!  Partial-lauyer control

      LOGICAL, INTENT (IN) ::           PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

!  Flux multipliers and quadrature

      DOUBLE PRECISION, INTENT (IN) ::  FLUX_MULTIPLIER
      DOUBLE PRECISION, INTENT (IN) ::  SURFACE_FACTOR
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_STRMWTS ( MAXSTREAMS )

!  User-stream and discrete ordinate transmittances

      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_USERM   ( MAXLAYERS,      MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_UTUP_USERM   ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Cumulative source term

      DOUBLE PRECISION, INTENT (IN) ::  CUMSOURCE_UP   ( MAX_USER_STREAMS, MAXSTOKES, 0:MAXLAYERS )

!  Surface reflectance variables
!  4/15/20. Version 2.8.2. BRDF arrays defined locally, each Fourier, remove MAXMOMENTS dimension

      DOUBLE PRECISION, INTENT (IN) ::  ALBEDO
      DOUBLE PRECISION, INTENT (IN) ::  BRDF_F      ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  USER_BRDF_F ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS ) 
      
!  Reflected Direct beam and surface-leaving solutions

      DOUBLE PRECISION, INTENT (IN) ::  RF_USER_DIRECT_BEAM ( MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  SL_USERTERM         ( MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )
!      DOUBLE PRECISION, INTENT (IN) ::  USER_DIRECT_BEAM ( MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )

!  Linearized transmittance flux for water-leaving

      DOUBLE PRECISION, INTENT (IN) ::  LP_TRANS_ATMOS_FINAL (MAXBEAMS,MAXLAYERS,MAX_ATMOSWFS)
      
!  RTE solutions (Homog.)

      INTEGER, INTENT (IN) ::           K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  LCON    ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  MCON    ( MAXSTRMSTKS, MAXLAYERS )

!  Thermal solutions

      DOUBLE PRECISION, INTENT (IN) ::  T_WLOWER ( MAXSTREAMS_2, MAXLAYERS )

!  User stream RTE solutions

      DOUBLE PRECISION, INTENT (IN) ::  UHOM_UPDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UHOM_UPUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) ::  UPAR_UP_1 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UPAR_UP_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )

!  Multipliers

      DOUBLE PRECISION, INTENT (IN) ::  HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_HMULT_UU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_HMULT_UD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )

      DOUBLE PRECISION, INTENT (IN) ::  EMULT_UP    ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )

!  Linearized transmittances

      DOUBLE PRECISION, INTENT (IN) ::  L_T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  L_T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  L_T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  L_DELTAU_VERT    ( MAX_ATMOSWFS, MAXLAYERS )

!  Linearized homogeneous and sdiscrete ordinate solution

      DOUBLE PRECISION, INTENT (IN) ::  L_T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  L_SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  L_SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  L_WLOWER    ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  NCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  PCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized thermal solutions

      DOUBLE PRECISION, INTENT (IN) ::  L_T_WLOWER        ( MAXSTREAMS_2,     MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  L_LAYER_TSUP_UP   ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) ::  L_LAYER_TSUP_UTUP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  linearized User homogeneous and beam solutions

      DOUBLE PRECISION, INTENT (IN) :: L_UHOM_UPDN  ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UHOM_UPUP  ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UPAR_UP_1  ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_UPAR_UP_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized Homog. solution multipliers

      DOUBLE PRECISION, INTENT (IN) :: L_HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UT_HMULT_UU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UT_HMULT_UD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Linearized Beam solution multipliers

      DOUBLE PRECISION, INTENT (IN) :: LP_EMULT_UP    ( MAX_USER_STREAMS, MAXLAYERS,      MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  OUTPUT
!  ======

!  Linearized BOA source term for Thermal-transonly

      DOUBLE PRECISION, INTENT (OUT) :: L_BOA_THTONLY_SOURCE ( MAXSTREAMS, MAX_ATMOSWFS )

!  Linearized output (Ostensibly output)

      DOUBLE PRECISION, INTENT (INOUT) :: PROFILEWF_F ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                                        MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  local variables
!  ---------------

!  Help

      LOGICAL ::          SFLAG
      INTEGER ::          M, N, NUT, NSTART, NUT_PREV, NLEVEL, O1
      INTEGER ::          UTA, UM, NV, NVPARS, Q, NC, UT, IB, LUM
      DOUBLE PRECISION :: L_FINAL_SOURCE

!  Local arrays

      DOUBLE PRECISION :: L_CUMUL_SOURCE ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_BOA_MSSOURCE ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_BOA_DBSOURCE ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_LAYER_SOURCE ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

!  Proxies

      NV     = LAYER_TO_VARY
      NVPARS = NV_PARAMETERS
      M   = FOURIER
      IB  = IBEAM
      LUM = 1

!  Zero all Fourier components - New rule, better for safety
!    Only did this for components close to zenith (formerly)

      IF ( DO_USER_STREAMS ) THEN
        DO UTA = 1, N_USER_LEVELS
          DO Q = 1, NV_PARAMETERS
            IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
              DO UM = 1, N_USER_STREAMS
                DO O1 = 1, NSTOKES
                  PROFILEWF_F(Q,NV,UTA,UM,IB,O1,UPIDX) = ZERO
                ENDDO
              ENDDO
            ELSE
              DO O1 = 1, NSTOKES
                PROFILEWF_F(Q,NV,UTA,LUM,IB,O1,UPIDX) = ZERO
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDIF

!  Initialize post-processing recursion
!  ====================================

      IF ( DO_USER_STREAMS ) THEN

!  Get the linearized BOA source terms (diffuse and direct)

        CALL GET_LP_BOASOURCE ( &
          DO_USER_STREAMS, DO_OBSERVATION_GEOMETRY, DO_INCLUDE_MVOUTPUT,   & ! Input flags
          DO_SOLAR_SOURCES,   DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY, & ! Input flags
          DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, DO_DBCORRECTION,      & ! Input flags
          DO_INCLUDE_DIRECTRF, DO_INCLUDE_DIRECTSL, NV, NV_PARAMETERS, FOURIER, & ! Input flags, linearization control
          IBEAM, NSTOKES, NSTREAMS, NLAYERS, N_USER_STREAMS, LOCAL_UM_START,    & ! Input Numbers
          DELTAU_SLANT, QUAD_WEIGHTS, QUAD_STRMWTS, MUELLER_INDEX,         & ! Input bookkeeping
          SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,                       & ! Input Direct beam, emiss.
          RF_USER_DIRECT_BEAM, SL_USERTERM, LP_TRANS_ATMOS_FINAL,            & ! Input surface radiance
          K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN,           & ! Input Homog solutions
          LCON, MCON, T_DELT_DISORDS, T_WLOWER,                            & ! Input Thermal and PI
          L_DELTAU_VERT, L_T_DELT_DISORDS, L_T_WLOWER,                     & ! Linearized thermal and misc.
          L_SOLA_XPOS, L_SOLB_XNEG,  L_T_DELT_EIGEN, L_WLOWER, NCON, PCON, & ! Linearized DsOr solutions
          L_BOA_THTONLY_SOURCE, L_BOA_MSSOURCE, L_BOA_DBSOURCE )             ! OUTPUT BOA terms (linearized)

!  Set the cumulative source term equal to the BOA sum

        IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
           DO O1 = 1, NSTOKES
            DO Q = 1, NV_PARAMETERS
              L_CUMUL_SOURCE(UM,O1,Q) = L_BOA_MSSOURCE(UM,O1,Q) + L_BOA_DBSOURCE(UM,O1,Q)
            ENDDO
           ENDDO
          ENDDO
        ELSE
          DO O1 = 1, NSTOKES
           DO Q = 1, NV_PARAMETERS
             L_CUMUL_SOURCE(IB,O1,Q) = L_BOA_MSSOURCE(IB,O1,Q) + L_BOA_DBSOURCE(IB,O1,Q)
           ENDDO
          ENDDO
        ENDIF

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
            SFLAG = DO_LAYER_SCATTERING(FOURIER,N)
            NC = NLAYERS + 1 - N

            CALL LP_WHOLELAYER_STERM_UP ( &
              DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,   & ! Input flags
              DO_OBSERVATION_GEOMETRY, DO_MSMODE_VLIDORT, SFLAG,               & ! Input flags
              N, M, IB, NV, NVPARS, NSTOKES, N_USER_STREAMS, LOCAL_UM_START,   & ! Input numbers
              K_REAL, K_COMPLEX, LCON, MCON, HMULT_1, HMULT_2, EMULT_UP,       & ! Input solutions and multipliers
              UHOM_UPDN, UHOM_UPUP, UPAR_UP_1, UPAR_UP_2,                      & ! Input user-solutions
              NCON, PCON, L_UHOM_UPDN, L_UHOM_UPUP, L_UPAR_UP_1, LP_UPAR_UP_2, & ! Input linearized solutions
              L_HMULT_1, L_HMULT_2, LP_EMULT_UP, L_LAYER_TSUP_UP,              & ! Input linearized multiplieres
              L_LAYER_SOURCE )                                                   ! OUTPUT linearized source term

            IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
              IF ( N.EQ.NV ) THEN
                DO UM = LOCAL_UM_START, N_USER_STREAMS
                  DO O1 = 1, NSTOKES
                    DO Q = 1, NV_PARAMETERS
                      L_CUMUL_SOURCE(UM,O1,Q) = L_LAYER_SOURCE(UM,O1,Q) +   T_DELT_USERM(N,UM)   * L_CUMUL_SOURCE(UM,O1,Q) &
                                                                        + L_T_DELT_USERM(N,UM,Q) * CUMSOURCE_UP(UM,O1,NC-1)
                    ENDDO
                  ENDDO
                ENDDO
              ELSE IF ( N.NE.NV ) THEN
                DO UM = LOCAL_UM_START, N_USER_STREAMS
                  DO O1 = 1, NSTOKES
                    DO Q = 1, NV_PARAMETERS
                      L_CUMUL_SOURCE(UM,O1,Q) = L_LAYER_SOURCE(UM,O1,Q) + T_DELT_USERM(N,UM) * L_CUMUL_SOURCE(UM,O1,Q)
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF
            ELSE
              IF ( N.EQ.NV ) THEN
                DO O1 = 1, NSTOKES
                  DO Q = 1, NV_PARAMETERS
                    L_CUMUL_SOURCE(IB,O1,Q) = L_LAYER_SOURCE(IB,O1,Q) +   T_DELT_USERM(N,IB)   * L_CUMUL_SOURCE(IB,O1,Q) &
                                                                      + L_T_DELT_USERM(N,IB,Q) * CUMSOURCE_UP(IB,O1,NC-1)
                  ENDDO
                ENDDO
              ELSE IF ( N.NE.NV ) THEN
                DO O1 = 1, NSTOKES
                  DO Q = 1, NV_PARAMETERS
                    L_CUMUL_SOURCE(IB,O1,Q) = L_LAYER_SOURCE(IB,O1,Q) +  T_DELT_USERM(N,IB) * L_CUMUL_SOURCE(IB,O1,Q)
                  ENDDO
                ENDDO
              ENDIF
            ENDIF

          ENDDO
        ENDIF

!  Offgrid output
!  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT    = PARTLAYERS_OUTINDEX(UTA)
          N     = PARTLAYERS_LAYERIDX(UT)
          SFLAG = DO_LAYER_SCATTERING(FOURIER,N)

!  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN

            CALL LP_PARTLAYER_STERM_UP ( &
              DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,   & ! Input flags
              DO_OBSERVATION_GEOMETRY, DO_MSMODE_VLIDORT, SFLAG,               & ! Input flags
              N, UT, IB, NV, NVPARS, NSTOKES, N_USER_STREAMS, LOCAL_UM_START,  & ! Input numbers
              K_REAL, K_COMPLEX, LCON, MCON, UHOM_UPDN, UHOM_UPUP, UPAR_UP_1,  & ! Input user-solutions
              UPAR_UP_2, UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP,                & ! Input multipliers
              NCON, PCON, L_UHOM_UPDN, L_UHOM_UPUP, L_UPAR_UP_1, LP_UPAR_UP_2, & ! Linearized solutions
              L_UT_HMULT_UU, L_UT_HMULT_UD, LP_UT_EMULT_UP, L_LAYER_TSUP_UTUP, & ! Linearized multipliers and thermal
              L_LAYER_SOURCE )                                                   ! OUTPUT

            IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
              IF ( N.EQ.NV ) THEN
                DO UM = LOCAL_UM_START, N_USER_STREAMS
                  DO O1 = 1, NSTOKES
                    DO Q = 1, NV_PARAMETERS
                      L_FINAL_SOURCE = L_LAYER_SOURCE(UM,O1,Q) &
                                            +   T_UTUP_USERM(UT,UM)   * L_CUMUL_SOURCE(UM,O1,Q) &
                                            + L_T_UTUP_USERM(UT,UM,Q) *   CUMSOURCE_UP(UM,O1,NC)
                      PROFILEWF_F(Q,NV,UTA,UM,IB,O1,UPIDX) = FLUX_MULTIPLIER * L_FINAL_SOURCE
                    ENDDO
                  ENDDO
                ENDDO
              ELSE IF ( N.NE.NV ) THEN
                DO UM = LOCAL_UM_START, N_USER_STREAMS
                  DO O1 = 1, NSTOKES
                    DO Q = 1, NV_PARAMETERS
                      L_FINAL_SOURCE = L_LAYER_SOURCE(UM,O1,Q) + T_UTUP_USERM(UT,UM)  * L_CUMUL_SOURCE(UM,O1,Q)
                      PROFILEWF_F(Q,NV,UTA,UM,IB,O1,UPIDX) =  FLUX_MULTIPLIER * L_FINAL_SOURCE
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF
            ELSE
              IF ( N.EQ.NV ) THEN
                DO O1 = 1, NSTOKES
                  DO Q = 1, NV_PARAMETERS
                    L_FINAL_SOURCE = L_LAYER_SOURCE(IB,O1,Q) &
                                          +   T_UTUP_USERM(UT,IB)   * L_CUMUL_SOURCE(IB,O1,Q) &
                                          + L_T_UTUP_USERM(UT,IB,Q) *   CUMSOURCE_UP(IB,O1,NC)
                    PROFILEWF_F(Q,NV,UTA,LUM,IB,O1,UPIDX) = FLUX_MULTIPLIER * L_FINAL_SOURCE
                  ENDDO
                ENDDO
              ELSE IF ( N.NE.NV ) THEN
                DO O1 = 1, NSTOKES
                  DO Q = 1, NV_PARAMETERS
                    L_FINAL_SOURCE = L_LAYER_SOURCE(IB,O1,Q) + T_UTUP_USERM(UT,IB) * L_CUMUL_SOURCE(IB,O1,Q)
                    PROFILEWF_F(Q,NV,UTA,LUM,IB,O1,UPIDX) = FLUX_MULTIPLIER * L_FINAL_SOURCE
                  ENDDO
                ENDDO
              ENDIF
            ENDIF

          ENDIF

!  Ongrid output
!  -------------

        ELSE

!  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
            IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO O1 = 1, NSTOKES
                  DO Q = 1, NV_PARAMETERS
                    L_FINAL_SOURCE =  FLUX_MULTIPLIER * L_CUMUL_SOURCE(UM,O1,Q)
!                    if (DABS(L_FINAL_SOURCE).GT.1.0d-12 ) then
                      PROFILEWF_F(Q,NV,UTA,UM,IB,O1,UPIDX) = L_FINAL_SOURCE
!                    endif
                  ENDDO
                ENDDO
              ENDDO
            ELSE
              DO O1 = 1, NSTOKES
                DO Q = 1, NV_PARAMETERS
                  L_FINAL_SOURCE =  FLUX_MULTIPLIER * L_CUMUL_SOURCE(IB,O1,Q)
!                  if (DABS(L_FINAL_SOURCE).GT.1.0d-12 ) then
                    PROFILEWF_F(Q,NV,UTA,LUM,IB,O1,UPIDX) = L_FINAL_SOURCE
!                  endif
                ENDDO
              ENDDO
            ENDIF
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
      END SUBROUTINE UPUSER_PROFILEWF

!

      SUBROUTINE DNUSER_PROFILEWF ( &
        DO_USER_STREAMS,  DO_OBSERVATION_GEOMETRY, DO_MSMODE_VLIDORT,              & ! Input flags (RT mode)
        DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS,   DO_THERMAL_TRANSONLY,           & ! Input flags (sources)
        DO_LAYER_SCATTERING, FOURIER, IBEAM, LVARY, NV_PARAMETERS, NSTOKES,        & ! Input numbers (basic)
        N_USER_STREAMS, LOCAL_UM_START, N_USER_LEVELS, UTAU_LEVEL_MASK_DN,         & ! Input bookkeeping + levels
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,              & ! Input partial-layer control
        FLUX_MULTIPLIER, T_DELT_USERM, T_UTDN_USERM, CUMSOURCE_DN,                 & ! Input Transmittances, Flux
        K_REAL, K_COMPLEX, LCON, MCON, UHOM_DNDN, UHOM_DNUP, UPAR_DN_1, UPAR_DN_2, & ! Input RT and User solutions
        HMULT_1, HMULT_2, EMULT_DN, UT_HMULT_DD, UT_HMULT_DU, UT_EMULT_DN,         & ! Input multipliers
        L_T_DELT_USERM, L_T_UTDN_USERM, L_LAYER_TSUP_DN,  L_LAYER_TSUP_UTDN,       & ! Input linearized
        NCON, PCON, L_UHOM_DNDN, L_UHOM_DNUP, L_UPAR_DN_1, LP_UPAR_DN_2,           & ! Linearized homog. and beam solutions
        L_HMULT_1, L_HMULT_2, LP_EMULT_DN,                                         & ! Linearized Homog. and Beam multipliers
        L_UT_HMULT_DD, L_UT_HMULT_DU, LP_UT_EMULT_DN,                              & ! Linearized Homog. and Beam multipliers
        PROFILEWF_F )                                                                ! OUTPUT
       
!   Streamlined for Version 2.8. 7/6/16

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXSTREAMS, MAXBEAMS, MAXSTOKES, MAXLAYERS, MAX_PARTLAYERS, MAX_ATMOSWFS, &
                                 MAX_USER_STREAMS, MAX_USER_LEVELS, MAX_USER_VZANGLES, MAX_SZANGLES, MAX_DIRECTIONS,   &
                                 MAXSTOKES_SQ, MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, ZERO, DNIDX

      IMPLICIT NONE

!  Input flags

      LOGICAL, INTENT (IN) ::           DO_USER_STREAMS
      LOGICAL, INTENT (IN) ::           DO_OBSERVATION_GEOMETRY

      LOGICAL, INTENT (IN) ::           DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::           DO_THERMAL_TRANSONLY

      LOGICAL, INTENT (IN) ::           DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )
      LOGICAL, INTENT (IN) ::           DO_MSMODE_VLIDORT

!  removed, Version 2.8, 7.17.16
!      LOGICAL, INTENT (IN) ::           DO_QUAD_OUTPUT


!  numbers

      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           IBEAM
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      INTEGER, INTENT (IN) ::           LVARY, NV_PARAMETERS

!  Level output control + bookkeeping

      INTEGER, INTENT (IN) ::           LOCAL_UM_START
      INTEGER, INTENT (IN) ::           N_USER_LEVELS
      INTEGER, INTENT (IN) ::           UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )

!  Partial-lauyer control

      LOGICAL, INTENT (IN) ::           PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::           PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

!  Flux multipliers and User-stream transmittances

      DOUBLE PRECISION, INTENT (IN) ::  FLUX_MULTIPLIER
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_USERM   ( MAXLAYERS,      MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  T_UTDN_USERM   ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Fourier contributions to radiance

      DOUBLE PRECISION, INTENT (IN) ::  CUMSOURCE_DN ( MAX_USER_STREAMS, MAXSTOKES, 0:MAXLAYERS )

!  RTE solutions

      INTEGER, INTENT (IN) ::           K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  LCON    ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  MCON    ( MAXSTRMSTKS, MAXLAYERS )

!  User stream RTE solutions

      DOUBLE PRECISION, INTENT (IN) ::  UHOM_DNDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UHOM_DNUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) ::  UPAR_DN_1 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UPAR_DN_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )

!  Multipliers

      DOUBLE PRECISION, INTENT (IN) ::  HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_HMULT_DD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_HMULT_DU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )

      DOUBLE PRECISION, INTENT (IN) ::  EMULT_DN    ( MAX_USER_STREAMS, MAXLAYERS,      MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )

!  Linearized inputs

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_LAYER_TSUP_DN     ( MAX_USER_STREAMS, MAXLAYERS,      MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_LAYER_TSUP_UTDN   ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
!  Linearized solutions

      DOUBLE PRECISION, INTENT (IN) :: NCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UHOM_DNDN  ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UHOM_DNUP  ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UPAR_DN_1  ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_UPAR_DN_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized Multipliers

      DOUBLE PRECISION, INTENT (IN) :: L_HMULT_1     ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_HMULT_2     ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UT_HMULT_DU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UT_HMULT_DD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (IN) :: LP_EMULT_DN    ( MAX_USER_STREAMS, MAXLAYERS,      MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized output (Ostensibly output)

      DOUBLE PRECISION, INTENT (INOUT) :: PROFILEWF_F ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                                        MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  local variables
!  ---------------

!  help

      LOGICAL ::          SFLAG
      INTEGER ::          N, NUT, NSTART, NUT_PREV, NLEVEL, O1
      INTEGER ::          UTA, UM, NV, NVPARS, Q, NC, UT, IB, M, LUM
      DOUBLE PRECISION :: L_FINAL_SOURCE

!  Local arrays

      DOUBLE PRECISION :: L_CUMUL_SOURCE ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_TOA_SOURCE   ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_LAYER_SOURCE ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

!  Proxies

      NV     = LVARY
      NVPARS = NV_PARAMETERS
      IB  = IBEAM
      M   = FOURIER
      LUM = 1

!  Zero all Fourier component output

      IF ( DO_USER_STREAMS ) THEN
        DO UTA = 1, N_USER_LEVELS
          DO Q = 1, NV_PARAMETERS
            IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
!mick fix 1/9/2013 - changed limit of do loop
              !DO UM = 1, LOCAL_UM_START
              DO UM = 1, N_USER_STREAMS
                DO O1 = 1, NSTOKES
                  PROFILEWF_F(Q,NV,UTA,UM,IB,O1,DNIDX) = ZERO
                ENDDO
              ENDDO
            ELSE
              DO O1 = 1, NSTOKES
                PROFILEWF_F(Q,NV,UTA,LUM,IB,O1,DNIDX) = ZERO
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDIF

!  Initialize post-processing recursion
!  ====================================

!  Get the linearized TOA source terms

      IF ( DO_USER_STREAMS ) THEN
        CALL GET_LP_TOASOURCE ( &
          DO_OBSERVATION_GEOMETRY, IB, NV_PARAMETERS, &
          NSTOKES, N_USER_STREAMS, LOCAL_UM_START,    &
          L_TOA_SOURCE )

        IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              DO Q = 1, NV_PARAMETERS
                L_CUMUL_SOURCE(UM,O1,Q) = L_TOA_SOURCE(UM,O1,Q)
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO O1 = 1, NSTOKES
            DO Q = 1, NV_PARAMETERS
              L_CUMUL_SOURCE(IB,O1,Q) = L_TOA_SOURCE(IB,O1,Q)
            ENDDO
          ENDDO
        ENDIF
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
            SFLAG = DO_LAYER_SCATTERING(FOURIER,N)
            NC = N

            CALL LP_WHOLELAYER_STERM_DN ( &
              DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,   & ! Input flags
              DO_OBSERVATION_GEOMETRY, DO_MSMODE_VLIDORT, SFLAG,               & ! Input flags
              N, M, IB, NV, NVPARS, NSTOKES, N_USER_STREAMS, LOCAL_UM_START,   & ! Input numbers
              K_REAL, K_COMPLEX, LCON, MCON, HMULT_1, HMULT_2, EMULT_DN,       & ! Input solutions and multipliers
              UHOM_DNDN, UHOM_DNUP, UPAR_DN_1, UPAR_DN_2,                      & ! Input user-solutions
              NCON, PCON, L_UHOM_DNDN, L_UHOM_DNUP, L_UPAR_DN_1, LP_UPAR_DN_2, & ! Input linearized solutions
              L_HMULT_1, L_HMULT_2, LP_EMULT_DN, L_LAYER_TSUP_DN,              & ! Input linearized multiplieres
              L_LAYER_SOURCE )                                                   ! OUTPUT

            IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
              IF ( N.EQ.NV ) THEN
                DO UM = LOCAL_UM_START, N_USER_STREAMS
                  DO O1 = 1, NSTOKES
                    DO Q = 1, NV_PARAMETERS
                      L_CUMUL_SOURCE(UM,O1,Q) = L_LAYER_SOURCE(UM,O1,Q) + T_DELT_USERM(N,UM)   * L_CUMUL_SOURCE(UM,O1,Q) &
                                                                      + L_T_DELT_USERM(N,UM,Q) * CUMSOURCE_DN(UM,O1,NC-1)
                    ENDDO
                  ENDDO
                ENDDO
              ELSE IF ( N.NE.NV ) THEN
                DO UM = LOCAL_UM_START, N_USER_STREAMS
                  DO O1 = 1, NSTOKES
                    DO Q = 1, NV_PARAMETERS
                      L_CUMUL_SOURCE(UM,O1,Q) = L_LAYER_SOURCE(UM,O1,Q) + T_DELT_USERM(N,UM) * L_CUMUL_SOURCE(UM,O1,Q)
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF
            ELSE
              IF ( N.EQ.NV ) THEN
                DO O1 = 1, NSTOKES
                  DO Q = 1, NV_PARAMETERS
                    L_CUMUL_SOURCE(IB,O1,Q) = L_LAYER_SOURCE(IB,O1,Q) + T_DELT_USERM(N,IB)   * L_CUMUL_SOURCE(IB,O1,Q) &
                                                                    + L_T_DELT_USERM(N,IB,Q) * CUMSOURCE_DN(IB,O1,NC-1)
                  ENDDO
                ENDDO
              ELSE IF ( N.NE.NV ) THEN
                DO O1 = 1, NSTOKES
                  DO Q = 1, NV_PARAMETERS
                    L_CUMUL_SOURCE(IB,O1,Q) = L_LAYER_SOURCE(IB,O1,Q) + T_DELT_USERM(N,IB) * L_CUMUL_SOURCE(IB,O1,Q)
                  ENDDO
                ENDDO
              ENDIF
            ENDIF

          ENDDO
        ENDIF

!  Offgrid output
!  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT    = PARTLAYERS_OUTINDEX(UTA)
          N     = PARTLAYERS_LAYERIDX(UT)
          SFLAG = DO_LAYER_SCATTERING(FOURIER,N)

!  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN
            CALL LP_PARTLAYER_STERM_DN ( &
              DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,   & ! Input flags
              DO_OBSERVATION_GEOMETRY, DO_MSMODE_VLIDORT, SFLAG,               & ! Input flags
              N, UT, IB, NV, NVPARS, NSTOKES, N_USER_STREAMS, LOCAL_UM_START,  & ! Input numbers
              K_REAL, K_COMPLEX, LCON, MCON, UHOM_DNDN, UHOM_DNUP, UPAR_DN_1,  & ! Input user-solutions
              UPAR_DN_2, UT_HMULT_DD, UT_HMULT_DU, UT_EMULT_DN,                & ! Input multipliers
              NCON, PCON, L_UHOM_DNDN, L_UHOM_DNUP, L_UPAR_DN_1, LP_UPAR_DN_2, & ! Input linearized solutions
              L_UT_HMULT_DU, L_UT_HMULT_DD, LP_UT_EMULT_DN, L_LAYER_TSUP_UTDN, & ! Input linearized Multipliers
              L_LAYER_SOURCE )                                                   ! OUTPUT

            IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
              IF ( N.EQ.NV ) THEN
                DO UM = LOCAL_UM_START, N_USER_STREAMS
                  DO O1 = 1, NSTOKES
                    DO Q = 1, NV_PARAMETERS
                      L_FINAL_SOURCE = L_LAYER_SOURCE(UM,O1,Q) + T_UTDN_USERM(UT,UM)   * L_CUMUL_SOURCE(UM,O1,Q) &
                                                             + L_T_UTDN_USERM(UT,UM,Q) *   CUMSOURCE_DN(UM,O1,NC)
                      PROFILEWF_F(Q,NV,UTA,UM,IB,O1,DNIDX) = FLUX_MULTIPLIER * L_FINAL_SOURCE
                    ENDDO
                  ENDDO
                ENDDO
              ELSE IF ( N.NE.NV ) THEN
                DO UM = LOCAL_UM_START, N_USER_STREAMS
                  DO O1 = 1, NSTOKES
                    DO Q = 1, NV_PARAMETERS
                      L_FINAL_SOURCE = L_LAYER_SOURCE(UM,O1,Q) + T_UTDN_USERM(UT,UM) * L_CUMUL_SOURCE(UM,O1,Q)
                      PROFILEWF_F(Q,NV,UTA,UM,IB,O1,DNIDX) = FLUX_MULTIPLIER * L_FINAL_SOURCE
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF
            ELSE
              IF ( N.EQ.NV ) THEN
                DO O1 = 1, NSTOKES
                  DO Q = 1, NV_PARAMETERS
                    L_FINAL_SOURCE = L_LAYER_SOURCE(IB,O1,Q) + T_UTDN_USERM(UT,IB)   * L_CUMUL_SOURCE(IB,O1,Q) &
                                                           + L_T_UTDN_USERM(UT,IB,Q) *   CUMSOURCE_DN(IB,O1,NC)
                    PROFILEWF_F(Q,NV,UTA,LUM,IB,O1,DNIDX) = FLUX_MULTIPLIER * L_FINAL_SOURCE
                  ENDDO
                ENDDO
              ELSE IF ( N.NE.NV ) THEN
                DO O1 = 1, NSTOKES
                  DO Q = 1, NV_PARAMETERS
                    L_FINAL_SOURCE = L_LAYER_SOURCE(IB,O1,Q) + T_UTDN_USERM(UT,IB) * L_CUMUL_SOURCE(IB,O1,Q)
                    PROFILEWF_F(Q,NV,UTA,LUM,IB,O1,DNIDX) = FLUX_MULTIPLIER * L_FINAL_SOURCE
                  ENDDO
                ENDDO
              ENDIF
            ENDIF

          ENDIF

!  Ongrid output
!  -------------

        ELSE

!  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
            IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO O1 = 1, NSTOKES
                  DO Q = 1, NV_PARAMETERS
                    PROFILEWF_F(Q,NV,UTA,UM,IB,O1,DNIDX) = FLUX_MULTIPLIER * L_CUMUL_SOURCE(UM,O1,Q)
                  ENDDO
                ENDDO
              ENDDO
            ELSE
              DO O1 = 1, NSTOKES
                DO Q = 1, NV_PARAMETERS
                  PROFILEWF_F(Q,NV,UTA,LUM,IB,O1,DNIDX) = FLUX_MULTIPLIER * L_CUMUL_SOURCE(IB,O1,Q)
                ENDDO
              ENDDO
            ENDIF
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
      END SUBROUTINE DNUSER_PROFILEWF

!

      SUBROUTINE GET_LP_TOASOURCE ( &
        DO_OBSERVATION_GEOMETRY, IBEAM, NV_PARAMETERS, &
        NSTOKES, N_USER_STREAMS, LOCAL_UM_START,       &
        L_TOA_SOURCE )

      USE VLIDORT_PARS_m, only : MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS, ZERO

      IMPLICIT NONE

!  Input

      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY
      INTEGER, INTENT (IN) ::          IBEAM
      INTEGER, INTENT (IN) ::          NV_PARAMETERS
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START

!  Output

      DOUBLE PRECISION, INTENT (OUT) :: L_TOA_SOURCE ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

!  local variables

      INTEGER :: UM, Q, O1

!  initialise TOA source function

      IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO Q = 1, NV_PARAMETERS
            DO O1 = 1, NSTOKES
              L_TOA_SOURCE(UM,O1,Q) = ZERO
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO Q = 1, NV_PARAMETERS
          DO O1 = 1, NSTOKES
            L_TOA_SOURCE(IBEAM,O1,Q) = ZERO
          ENDDO
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE GET_LP_TOASOURCE

!

      SUBROUTINE GET_LP_BOASOURCE ( &
          DO_USER_STREAMS, DO_OBSERVATION_GEOMETRY, DO_INCLUDE_MVOUTPUT,   & ! Input flags
          DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,   & ! Input flags
          DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, DO_DBCORRECTION,      & ! Input flags
          DO_INCLUDE_DIRECTRF, DO_INCLUDE_DIRECTSL, NV, NV_PARAMETERS, FOURIER, & ! Input flags, linearization control
          IBEAM, NSTOKES, NSTREAMS, NLAYERS, N_USER_STREAMS, LOCAL_UM_START,    & ! Input Numbers
          DELTAU_SLANT, QUAD_WEIGHTS, QUAD_STRMWTS, MUELLER_INDEX,         & ! Input bookkeeping
          SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,                       & ! Input Direct beam, emiss.
          RF_USER_DIRECT_BEAM, SL_USERTERM, LP_TRANS_ATMOS_FINAL,            & ! Input surface radiance
          K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN,           & ! Input Homog solutions
          LCON, MCON, T_DELT_DISORDS, T_WLOWER,                            & ! Input Thermal and PI
          L_DELTAU_VERT, L_T_DELT_DISORDS, L_T_WLOWER,                     & ! Linearized thermal and misc.
          L_SOLA_XPOS, L_SOLB_XNEG,  L_T_DELT_EIGEN, L_WLOWER, NCON, PCON, & ! Linearized DsOr solutions
          L_BOA_THTONLY_SOURCE, L_BOA_MSSOURCE, L_BOA_DBSOURCE )             ! OUTPUT BOA terms (linearized)
        
!  Linearized Bottom of the atmosphere source term

!  4/15/20. Version 2.8.2. BRDF arrays are defined locally, each Fourier.

!  module, dimensions and numbers
!  4/15/20. Version 2.8.2. Remove MAXMOMENTS dimension.

      USE VLIDORT_PARS_m, Only : MAXSTREAMS, MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES, MAXLAYERS,  &
                                 MAX_ATMOSWFS, MAXSTOKES_SQ, MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, ZERO, ONE

      IMPLICIT NONE
      
!  INPUTS
!  ======

!  Flags

      LOGICAL, INTENT (IN) ::           DO_USER_STREAMS
      LOGICAL, INTENT (IN) ::           DO_OBSERVATION_GEOMETRY
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_MVOUTPUT

      LOGICAL, INTENT (IN) ::           DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::           DO_THERMAL_TRANSONLY

      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::           DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT (IN) ::           DO_DBCORRECTION

!  4/9/19. Replacement of DIRECTBEAM by two separate flags
      
      LOGICAL, intent(in)  ::           DO_INCLUDE_DIRECTRF
      LOGICAL, intent(in)  ::           DO_INCLUDE_DIRECTSL
!      LOGICAL, INTENT (IN) ::           DO_INCLUDE_DIRECTBEAM
      
!  Numbers

      INTEGER, INTENT (IN) ::           NV, NV_PARAMETERS
      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           IBEAM
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      INTEGER, INTENT (IN) ::           LOCAL_UM_START

!  Bookkeeping (quadrature etc)

      DOUBLE PRECISION, INTENT (IN) ::  DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS ) 
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_STRMWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::           MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )

!  Surface reflectance
!  4/15/20. Version 2.8.2. BRDF arrays defined locally, each Fourier, remove MAXMOMENTS dimension

      DOUBLE PRECISION, INTENT (IN) ::  SURFACE_FACTOR
      DOUBLE PRECISION, INTENT (IN) ::  ALBEDO
      DOUBLE PRECISION, INTENT (IN) ::  BRDF_F      ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  USER_BRDF_F ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )

!  Reflected Direct beam and surface-leaving solutions

      DOUBLE PRECISION, INTENT (IN) ::  RF_USER_DIRECT_BEAM ( MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  SL_USERTERM         ( MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )
!      DOUBLE PRECISION, INTENT (IN) ::  USER_DIRECT_BEAM ( MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )

!  Linearized transmittance flux for water-leaving

      DOUBLE PRECISION, INTENT (IN) ::  LP_TRANS_ATMOS_FINAL (MAXBEAMS,MAXLAYERS,MAX_ATMOSWFS)
      
!  Homogeneous solutions

      INTEGER, INTENT (IN) ::           K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )

!  Particular integrals and thermal

      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_WLOWER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  MCON ( MAXSTRMSTKS, MAXLAYERS )

!  Linearized inputs

      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT  ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_EIGEN  ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_DISORDS( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_WLOWER    ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: NCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_WLOWER  ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )

!  local variables
!  ---------------

!  Output
!  ======

!  Quad-angle output

      DOUBLE PRECISION, INTENT (OUT) :: L_BOA_THTONLY_SOURCE ( MAXSTREAMS, MAX_ATMOSWFS )

!  User-angle output

      DOUBLE PRECISION, INTENT (OUT) :: L_BOA_MSSOURCE ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_BOA_DBSOURCE ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

!  local variables
!  ---------------


      INTEGER ::          M, N, NN, J, I, UM, NELEMENTS
      INTEGER ::          Q, IB, O1, O2, OM, O11
      INTEGER ::          K, KO1, K0, K1, K2, LUM

      LOGICAL ::          DO_QTHTONLY
      DOUBLE PRECISION :: DOWN   (MAXSTREAMS,MAXSTOKES)
      DOUBLE PRECISION :: L_DOWN (MAXSTREAMS,MAXSTOKES,MAX_ATMOSWFS)

      DOUBLE PRECISION :: REFLEC, S_REFLEC, FAC, KMULT
      DOUBLE PRECISION :: SPAR, SHOM_R, SHOM_CR, SHOM
      DOUBLE PRECISION :: HOM1, HOM2, HOM1CR, HOM2CR, HOM3CR
      DOUBLE PRECISION :: LXR, NXR, PXR, LLXR, MLXR
      DOUBLE PRECISION :: LXR1, NXR1, PXR1, LLXR1, MLXR1
      DOUBLE PRECISION :: LXR2, NXR2, LLXR2

!  Starting section
!  ----------------

!  Fourier number, layer number

      M   = FOURIER
      N   = NLAYERS
      KO1 = K_REAL(N) + 1
      IB  = IBEAM
      LUM = 1

!  Special flag
!    Version 2.8. QUAD_OUTPUT flag removed
!      DO_QTHTONLY = ( DO_THERMAL_TRANSONLY ) .AND. ( DO_QUAD_OUTPUT .OR. DO_INCLUDE_MVOUTPUT )

      DO_QTHTONLY = DO_THERMAL_TRANSONLY .AND. DO_INCLUDE_MVOUTPUT

!  initialise linearized BOA source functions

      IF ( DO_USER_STREAMS ) THEN
        IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO Q = 1, NV_PARAMETERS
              DO O1 = 1, NSTOKES
                L_BOA_MSSOURCE(UM,O1,Q) = ZERO
                L_BOA_DBSOURCE(UM,O1,Q) = ZERO
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO Q = 1, NV_PARAMETERS
            DO O1 = 1, NSTOKES
              L_BOA_MSSOURCE(IB,O1,Q) = ZERO
              L_BOA_DBSOURCE(IB,O1,Q) = ZERO
            ENDDO
          ENDDO
        ENDIF
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

!  Add direct beam if flagged
!    For NV > 0, profile weighting functions

      IF ( DO_INCLUDE_DIRECTRF.AND..NOT.DO_DBCORRECTION .AND.DO_USER_STREAMS) THEN
         IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
               DO O1 = 1, NELEMENTS
                  FAC = - RF_USER_DIRECT_BEAM(UM,IB,O1) * DELTAU_SLANT(N,NV,IB)
                  DO Q = 1, NV_PARAMETERS
                     L_BOA_DBSOURCE(UM,O1,Q) = L_BOA_DBSOURCE(UM,O1,Q) + L_DELTAU_VERT(Q,NV) * FAC
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            DO O1 = 1, NELEMENTS
               FAC = - RF_USER_DIRECT_BEAM(LUM,IB,O1) * DELTAU_SLANT(N,NV,IB)
               DO Q = 1, NV_PARAMETERS
                  L_BOA_DBSOURCE(IB,O1,Q) = L_BOA_DBSOURCE(IB,O1,Q) + L_DELTAU_VERT(Q,NV) * FAC
               ENDDO
            ENDDO
         ENDIF
      ENDIF

!  Add surface-term linearization if flagged. O1 = 1 (I-component only)
      
      IF ( DO_INCLUDE_DIRECTSL .AND. DO_USER_STREAMS ) THEN
         O1 = 1
         IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
               FAC = SL_USERTERM(UM,IB,O1)
               DO Q = 1, NV_PARAMETERS
                  L_BOA_DBSOURCE(UM,O1,Q) = L_BOA_DBSOURCE(UM,O1,Q) +  LP_TRANS_ATMOS_FINAL(IB,NV,Q) * FAC
               ENDDO
            ENDDO
         ELSE
            FAC = SL_USERTERM(LUM,IB,O1)
            DO Q = 1, NV_PARAMETERS
               L_BOA_DBSOURCE(IB,O1,Q) = L_BOA_DBSOURCE(IB,O1,Q) +  LP_TRANS_ATMOS_FINAL(IB,NV,Q) * FAC
            ENDDO
         ENDIF
      ENDIF

!  Exit if no surface

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  Include surface only, from now on
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
            IF ( NV.EQ.NN ) THEN
              DO I = 1, NSTREAMS
                DO Q = 1, NV_PARAMETERS
                  L_DOWN(I,O1,Q) = L_DOWN(I,O1,Q) *   T_DELT_DISORDS(I,NN)     &
                                 +   DOWN(I,O1)   * L_T_DELT_DISORDS(I,NN,Q) + L_T_WLOWER(I,NN,Q)
                ENDDO
              ENDDO
            ELSE
              DO I = 1, NSTREAMS
                DO Q = 1, NV_PARAMETERS
                  L_DOWN(I,O1,Q) = L_DOWN(I,O1,Q) * T_DELT_DISORDS(I,NN)
                ENDDO
              ENDDO
            ENDIF
            DO I = 1, NSTREAMS
              DOWN(I,O1) = DOWN(I,O1)*T_DELT_DISORDS(I,NN) + T_WLOWER(I,NN)
            ENDDO
          ENDDO

!  Finalize

          DO I = 1, NSTREAMS
            DO Q = 1, NV_PARAMETERS
              L_DOWN(I,O1,Q) = QUAD_WEIGHTS(I) * L_DOWN(I,O1,Q)
            ENDDO
          ENDDO

!  Continuation point for avoiding the next section
!  Version 2.8 remove GOTO -->        GO TO 5432

!  2. Scattering solutions
!     %%%%%%%%%%%%%%%%%%%%

        ELSE

!  Two cases:
!  (a) If  NV = N, this is also the layer that is varying --> Extras!
!  (b) If  N > NV with variations in layer NV above N

          IF ( NV .EQ. N ) THEN

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
                    K0 = 2 * K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
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
                  IF ( DO_SOLAR_SOURCES .OR. (DO_INCLUDE_THERMEMISS .AND. N.EQ.NV) ) THEN
                    SPAR = L_WLOWER(I,O1,N,Q)
                  ENDIF

!  Final result

                  L_DOWN(I,O1,Q) = SPAR + SHOM

!  End loops over Q, O1 and I

                ENDDO
              ENDDO
            ENDDO

!  otherwise the varying layer is above the boundary layer

          ELSE IF ( NV.LT.N ) THEN

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
                    K0 = 2 * K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1

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
                  IF ( DO_SOLAR_SOURCES .OR. (DO_INCLUDE_THERMEMISS .AND. N.EQ.NV) ) THEN
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

!  Continuation point. Removed Version 2.8
! 5432   CONTINUE

!  End scattering vs. thermal-only

        ENDIF

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
          KMULT = SURFACE_FACTOR * ALBEDO
          IF ( FOURIER .EQ. 0 ) THEN
            O1 = 1
            DO Q = 1, NV_PARAMETERS
              REFLEC = KMULT * SUM(L_DOWN(1:NSTREAMS,O1,Q))
              IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
                DO UM = LOCAL_UM_START, N_USER_STREAMS
                  L_BOA_MSSOURCE(UM,O1,Q) = L_BOA_MSSOURCE(UM,O1,Q) + REFLEC
                ENDDO
              ELSE
                L_BOA_MSSOURCE(IB,O1,Q) = L_BOA_MSSOURCE(IB,O1,Q) + REFLEC
              ENDIF
              IF ( DO_QTHTONLY ) THEN
                DO I = 1, NSTREAMS
                  L_BOA_THTONLY_SOURCE(I,Q) = L_BOA_THTONLY_SOURCE(I,Q) + REFLEC
                ENDDO
              ENDIF
            ENDDO
          ENDIF
        ENDIF

!  .. integrate with BRDF reflectance function at user angles
!  @@@@@@@@ Rob Fix, 2/9/11, DO_QTHTONLY clause was absent
!  4/15/20. Version 2.8.2. BRDF arrays defined locally, remove M=Fourier index

        IF ( .not. DO_LAMBERTIAN_SURFACE ) THEN
          KMULT = SURFACE_FACTOR
          DO Q = 1, NV_PARAMETERS
            IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                DO O1 = 1, NSTOKES
                  REFLEC = ZERO
                  DO J = 1, NSTREAMS
                    S_REFLEC = ZERO
                    DO O2 = 1, NSTOKES
                      OM = MUELLER_INDEX(O1,O2)
                      S_REFLEC = S_REFLEC + L_DOWN(J,O2,Q) * USER_BRDF_F(OM,UM,J)
                    ENDDO
                    REFLEC = REFLEC + S_REFLEC
                  ENDDO
                  REFLEC = KMULT * REFLEC
                  L_BOA_MSSOURCE(UM,O1,Q) = L_BOA_MSSOURCE(UM,O1,Q) + REFLEC
                ENDDO
              ENDDO
            ELSE
              DO O1 = 1, NSTOKES
                REFLEC = ZERO
                DO J = 1, NSTREAMS
                  S_REFLEC = ZERO
                  DO O2 = 1, NSTOKES
                    OM = MUELLER_INDEX(O1,O2)
                    S_REFLEC = S_REFLEC + L_DOWN(J,O2,Q) * USER_BRDF_F(OM,IB,J)
                  ENDDO
                  REFLEC = REFLEC + S_REFLEC
                ENDDO
                REFLEC = KMULT * REFLEC
                L_BOA_MSSOURCE(IB,O1,Q) = L_BOA_MSSOURCE(IB,O1,Q) + REFLEC
              ENDDO
            ENDIF

            IF ( DO_QTHTONLY ) THEN
              O11 = 1
              DO I = 1, NSTREAMS
                REFLEC = ZERO
                DO J = 1, NSTREAMS
                  REFLEC = REFLEC + L_DOWN(J,O11,Q) * BRDF_F(O11,I,J)
                ENDDO
                REFLEC = KMULT * REFLEC
                L_BOA_THTONLY_SOURCE(I,Q) = L_BOA_THTONLY_SOURCE(I,Q) + REFLEC
              ENDDO
            ENDIF
          ENDDO
        ENDIF

!  End inclusion of surface terms

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE GET_LP_BOASOURCE

!

      SUBROUTINE LP_WHOLELAYER_STERM_UP ( &
          DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,   & ! Input flags
          DO_OBSERVATION_GEOMETRY, DO_MSMODE_VLIDORT, SOURCETERM_FLAG,     & ! Input flags
          GIVEN_LAYER, FOURIER, IBEAM, LAYER_TO_VARY, NV_PARAMETERS,       & ! Input numbers
          NSTOKES, N_USER_STREAMS, LOCAL_UM_START,                         & ! Input numbers
          K_REAL, K_COMPLEX, LCON, MCON, HMULT_1, HMULT_2, EMULT_UP,       & ! Input solutions and multipliers
          UHOM_UPDN, UHOM_UPUP, UPAR_UP_1, UPAR_UP_2,                      & ! Input user-solutions
          NCON, PCON, L_UHOM_UPDN, L_UHOM_UPUP, L_UPAR_UP_1, LP_UPAR_UP_2, & ! Input linearized solutions
          L_HMULT_1, L_HMULT_2, LP_EMULT_UP, L_LAYER_TSUP_UP,              & ! Input linearized multiplieres
          L_LAYERSOURCE )                                                    ! OUTPUT linearized source term

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXSTREAMS, MAXBEAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS, &
                                 MAX_USER_STREAMS, MAXEVALUES, MAXSTRMSTKS, ZERO, ONE, PI4

      IMPLICIT NONE

!  Inputs
!  ======

!  flags

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY

      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY
      LOGICAL, INTENT (IN) ::          DO_MSMODE_VLIDORT
      LOGICAL, INTENT (IN) ::          SOURCETERM_FLAG

!  Numbers

      INTEGER, INTENT (IN) ::          IBEAM, GIVEN_LAYER, FOURIER
      INTEGER, INTENT (IN) ::          LAYER_TO_VARY, NV_PARAMETERS
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START

!  Homog solutions

      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )

!  User solutions

      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UPAR_UP_1 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UPAR_UP_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )

!  Multipliers

      DOUBLE PRECISION, INTENT (IN) :: HMULT_1  ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_2  ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: EMULT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )

!  Linearized solutions

      DOUBLE PRECISION, INTENT (IN) :: NCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (IN) :: L_UHOM_UPDN  ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UHOM_UPUP  ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (IN) :: L_UPAR_UP_1  ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_UPAR_UP_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (IN) :: L_LAYER_TSUP_UP  ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized Multipliers

      DOUBLE PRECISION, INTENT (IN) :: L_HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_EMULT_UP( MAX_USER_STREAMS, MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Output
!  ------

      DOUBLE PRECISION, INTENT (OUT) :: L_LAYERSOURCE ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

!  local variables
!  ---------------

      INTEGER ::          N, NV, UM, O1, Q, IB, M
      INTEGER ::          K, KO1, K0, K1, K2, LUM
      INTEGER ::          UM_START, UM_END
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, SPAR, H1, H2, TM
      DOUBLE PRECISION :: LUXR, MUXR, NUXR, PUXR, LLUXR, MLUXR
      DOUBLE PRECISION :: LUXR1, MUXR1, NUXR1, PUXR1, LLUXR1, MLUXR1
      DOUBLE PRECISION :: LUXR2, MUXR2, NUXR2, PUXR2, LLUXR2, MLUXR2

!  Proxies

      N   = GIVEN_LAYER
      KO1 = K_REAL(N) + 1
      NV  = LAYER_TO_VARY
      IB  = IBEAM
      M   = FOURIER        ! Only for debugging
      LUM = 1

!  Local ranges

      IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
        UM_START = LOCAL_UM_START
        UM_END   = N_USER_STREAMS
      ELSE
        UM_START = IB
        UM_END   = IB
      ENDIF

!   Very important to zero both output terms (bug solved 12/29/05)

      DO Q = 1, NV_PARAMETERS
        DO UM = UM_START,UM_END
          DO O1 = 1, NSTOKES
            L_LAYERSOURCE(UM,O1,Q)  = ZERO
          ENDDO
       ENDDO
      ENDDO

!  return if no source term
!      Need to go on if thermal transmittances only

      IF ( .NOT. SOURCETERM_FLAG .and. .NOT. DO_THERMAL_TRANSONLY ) RETURN

!  Avoid this section if thermal transmittance only
!      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789
!      Version 2.8, GOTO statement removed, 19 July 2016

      IF ( .not. DO_THERMAL_TRANSONLY ) THEN

!  Homogeneous solutions
!  =====================

!  Special case when N = NV
!  ------------------------

        IF ( N .EQ. NV ) THEN

!  Loop over user angles and STokes

          DO UM = UM_START,UM_END
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
                  K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1

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

        ELSE IF ( N.NE.NV ) THEN

!  Loop over user angles and STokes

          DO UM = UM_START,UM_END
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
                  K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1

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

!  End clause N.eq.NV

        ENDIF

!  Continuation point
! 6789 continue

!  End scattering-only solutions

      ENDIF

!  Add thermal emission term (direct and diffuse)
!     -----Modulus 1 if solar sources are included (taken care of earlier)
!     ----- Linearization only exists if N = K

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        O1 = 1 ; TM = ONE
        IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
        IF ( N.EQ.NV ) THEN
          DO UM = UM_START,UM_END
            DO Q = 1, NV_PARAMETERS
              L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + L_LAYER_TSUP_UP(UM,N,Q)*TM
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Particular and single scatter contributions
!  ===========================================

!  Lattice calculation

      IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN

!  Special case when N = NV
!  ------------------------

        IF ( N.EQ.NV ) THEN

!  add particular solution

          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              DO Q = 1, NV_PARAMETERS
                SPAR = LP_UPAR_UP_2(UM,O1,N,NV,Q) *   EMULT_UP(UM,N,IB) + &
                          UPAR_UP_2(UM,O1,N)      * LP_EMULT_UP(UM,N,NV,IB,Q)
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
                           UPAR_UP_1(UM,O1,N)   * LP_EMULT_UP(UM,N,NV,IB,Q)
                  L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SPAR
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!  Other cases when N > NV
!  -----------------------

        ELSE IF ( N .GT. NV ) THEN

!  add particular solution

          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              DO Q = 1, NV_PARAMETERS
                SPAR = LP_UPAR_UP_2(UM,O1,N,NV,Q) *   EMULT_UP(UM,N,IB) + &
                          UPAR_UP_2(UM,O1,N)      * LP_EMULT_UP(UM,N,NV,IB,Q)
                L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SPAR
              ENDDO
            ENDDO
          ENDDO

!  Add single scatter term if flagged

          IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                DO Q = 1, NV_PARAMETERS
                  SPAR =  UPAR_UP_1(UM,O1,N) * LP_EMULT_UP(UM,N,NV,IB,Q)
                  L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SPAR
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!  End layer clause

        ENDIF

!  Observational geometry calculation

      ELSE

!  Special case when N = NV
!  ------------------------

        IF ( N.EQ.NV ) THEN

!  add particular solution

          DO O1 = 1, NSTOKES
            DO Q = 1, NV_PARAMETERS
              SPAR = LP_UPAR_UP_2(IB,O1,N,NV,Q) *   EMULT_UP(LUM,N,IB) + &
                        UPAR_UP_2(IB,O1,N)      * LP_EMULT_UP(LUM,N,NV,IB,Q)
              L_LAYERSOURCE(IB,O1,Q) = L_LAYERSOURCE(IB,O1,Q) + SPAR
            ENDDO
          ENDDO

!  Add single scatter term if flagged

          IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
            DO O1 = 1, NSTOKES
              DO Q = 1, NV_PARAMETERS
                SPAR = L_UPAR_UP_1(IB,O1,N,Q) *   EMULT_UP(LUM,N,IB) + &
                         UPAR_UP_1(IB,O1,N)   * LP_EMULT_UP(LUM,N,NV,IB,Q)
                L_LAYERSOURCE(IB,O1,Q) = L_LAYERSOURCE(IB,O1,Q) + SPAR
              ENDDO
            ENDDO
          ENDIF

!  Other cases when N > NV
!  -----------------------

        ELSE IF ( N .GT. NV ) THEN

!  add particular solution

          DO O1 = 1, NSTOKES
            DO Q = 1, NV_PARAMETERS
              SPAR = LP_UPAR_UP_2(IB,O1,N,NV,Q) *   EMULT_UP(LUM,N,IB) + &
                        UPAR_UP_2(IB,O1,N)      * LP_EMULT_UP(LUM,N,NV,IB,Q)
              L_LAYERSOURCE(IB,O1,Q) = L_LAYERSOURCE(IB,O1,Q) + SPAR
            ENDDO
          ENDDO

!  Add single scatter term if flagged

          IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
            DO O1 = 1, NSTOKES
              DO Q = 1, NV_PARAMETERS
                SPAR =  UPAR_UP_1(IB,O1,N) * LP_EMULT_UP(LUM,N,NV,IB,Q)
                L_LAYERSOURCE(IB,O1,Q) = L_LAYERSOURCE(IB,O1,Q) + SPAR
              ENDDO
            ENDDO
          ENDIF

!  End layer clause

        ENDIF

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LP_WHOLELAYER_STERM_UP

!

      SUBROUTINE LP_WHOLELAYER_STERM_DN ( &
          DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,   & ! Input flags
          DO_OBSERVATION_GEOMETRY, DO_MSMODE_VLIDORT, SOURCETERM_FLAG,     & ! Input flags
          GIVEN_LAYER, FOURIER, IBEAM, LAYER_TO_VARY, NV_PARAMETERS,       & ! Input numbers
          NSTOKES, N_USER_STREAMS, LOCAL_UM_START,                         & ! Input numbers
          K_REAL, K_COMPLEX, LCON, MCON, HMULT_1, HMULT_2, EMULT_DN,       & ! Input solutions and multipliers
          UHOM_DNDN, UHOM_DNUP, UPAR_DN_1, UPAR_DN_2,                      & ! Input user-solutions
          NCON, PCON, L_UHOM_DNDN, L_UHOM_DNUP, L_UPAR_DN_1, LP_UPAR_DN_2, & ! Input linearized solutions
          L_HMULT_1, L_HMULT_2, LP_EMULT_DN, L_LAYER_TSUP_DN,              & ! Input linearized multiplieres
          L_LAYERSOURCE )                                                    ! OUTPUT linearized source term

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXSTREAMS, MAXBEAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS, &
                                 MAX_USER_STREAMS, MAXEVALUES, MAXSTRMSTKS, ZERO, ONE, PI4

      IMPLICIT NONE

!  Inputs
!  ======

!  flags

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY

      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY
      LOGICAL, INTENT (IN) ::          DO_MSMODE_VLIDORT
      LOGICAL, INTENT (IN) ::          SOURCETERM_FLAG

!  Numbers

      INTEGER, INTENT (IN) ::          IBEAM, GIVEN_LAYER, FOURIER
      INTEGER, INTENT (IN) ::          LAYER_TO_VARY, NV_PARAMETERS
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          LOCAL_UM_START

!  Homog solutions

      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )

!  User solutions

      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UPAR_DN_1 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UPAR_DN_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )

!  Multipliers

      DOUBLE PRECISION, INTENT (IN) :: HMULT_1  ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_2  ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: EMULT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )

!  Linearized solutions

      DOUBLE PRECISION, INTENT (IN) :: NCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (IN) :: L_UHOM_DNDN  ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UHOM_DNUP  ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (IN) :: L_UPAR_DN_1  ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_UPAR_DN_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (IN) :: L_LAYER_TSUP_DN  ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized Multipliers

      DOUBLE PRECISION, INTENT (IN) :: L_HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_EMULT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Output
!  ------

      DOUBLE PRECISION, INTENT (OUT) :: L_LAYERSOURCE ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

!  local variables
!  ---------------

      INTEGER ::          N, NV, UM, O1, Q, IB, M
      INTEGER ::          K, KO1, K0, K1, K2, LUM
      INTEGER ::          UM_START, UM_END
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, SPAR, H1, H2, TM
      DOUBLE PRECISION :: LUXR, MUXR, NUXR, PUXR, LLUXR, MLUXR
      DOUBLE PRECISION :: LUXR1, MUXR1, NUXR1, PUXR1, LLUXR1, MLUXR1
      DOUBLE PRECISION :: LUXR2, MUXR2, NUXR2, PUXR2, LLUXR2, MLUXR2

!  Proxies

      N   = GIVEN_LAYER
      KO1 = K_REAL(N) + 1
      NV  = LAYER_TO_VARY
      IB  = IBEAM
      M   = FOURIER   ! Only for debugging
      LUM = 1

!  Local Range

      IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
        UM_START = LOCAL_UM_START
        UM_END   = N_USER_STREAMS
      ELSE
        UM_START = IB
        UM_END   = IB
      ENDIF

!   Very important to zero both output terms (bug solved 12/29/05)

      DO Q = 1, NV_PARAMETERS
        DO UM = UM_START,UM_END
          DO O1 = 1, NSTOKES
            L_LAYERSOURCE(UM,O1,Q) = ZERO
          ENDDO
        ENDDO
      ENDDO

!  return if no source term
!      Need to go on if thermal transmittances only

      IF ( .NOT. SOURCETERM_FLAG .and. .NOT. DO_THERMAL_TRANSONLY ) RETURN

!  Avoid this section if thermal transmittance only
!      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789
!      Version 2.8, GOTO statement removed, 19 July 2016

      IF ( .not. DO_THERMAL_TRANSONLY ) THEN

!  Homogeneous solutions
!  =====================

!  Special case when N = NV
!  ------------------------

        IF ( N.EQ.NV ) THEN

!  Loop over user angles and STokes

          DO UM = UM_START,UM_END
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
                  K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1

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

        ELSE IF ( N.NE.NV ) THEN

!  Loop over user angles and STokes

          DO UM = UM_START,UM_END
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
                  K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1

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
! 6789 continue

!  End scattering clause

      ENDIF

!  Add thermal emission term (direct and diffuse)
!     ----- Modulus 1 if solar sources are included (taken care of earlier)
!     ----- Linearization only exists if N = K

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        O1 = 1 ; TM = ONE
        IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
        IF ( N.EQ.NV ) THEN
          DO UM = UM_START,UM_END
            DO Q = 1, NV_PARAMETERS
              L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + L_LAYER_TSUP_DN(UM,N,Q)*TM
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Particular and single scatter contributions
!  ===========================================

!  Lattice calculation

      IF (.not. DO_OBSERVATION_GEOMETRY ) THEN

!  Special case when N = NV
!  ------------------------

        IF ( N.EQ.NV ) THEN

!  add particular solution

          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              DO Q = 1, NV_PARAMETERS
                SPAR = LP_UPAR_DN_2(UM,O1,N,NV,Q) *   EMULT_DN(UM,N,IB) + &
                          UPAR_DN_2(UM,O1,N)      * LP_EMULT_DN(UM,N,NV,IB,Q)
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
                           UPAR_DN_1(UM,O1,N)   * LP_EMULT_DN(UM,N,NV,IB,Q)
                  L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SPAR
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!  Other cases when N > NV
!  -----------------------

        ELSE IF ( N.GT.NV ) THEN

!  add particular solution

          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              DO Q = 1, NV_PARAMETERS
              SPAR = LP_UPAR_DN_2(UM,O1,N,NV,Q) *   EMULT_DN(UM,N,IB) + &
                        UPAR_DN_2(UM,O1,N)      * LP_EMULT_DN(UM,N,NV,IB,Q)
               L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SPAR
              ENDDO
            ENDDO
          ENDDO

!  Add single scatter term if flagged

          IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                DO Q = 1, NV_PARAMETERS
                  SPAR =  UPAR_DN_1(UM,O1,N) * LP_EMULT_DN(UM,N,NV,IB,Q)
                  L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SPAR
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!  End layer clause

        ENDIF

!  Observational geometry calculation

      ELSE

!  Special case when N = NV
!  ------------------------

        IF ( N.EQ.NV ) THEN

!  add particular solution

          DO O1 = 1, NSTOKES
            DO Q = 1, NV_PARAMETERS
              SPAR = LP_UPAR_DN_2(IB,O1,N,NV,Q) *    EMULT_DN(LUM,N,IB) + &
                        UPAR_DN_2(IB,O1,N)      * LP_EMULT_DN(LUM,N,NV,IB,Q)
              L_LAYERSOURCE(IB,O1,Q) = L_LAYERSOURCE(IB,O1,Q) + SPAR
            ENDDO
          ENDDO

!  Add single scatter term if flagged

          IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
           DO O1 = 1, NSTOKES
             DO Q = 1, NV_PARAMETERS
               SPAR = L_UPAR_DN_1(IB,O1,N,Q) *   EMULT_DN(LUM,N,IB) + &
                        UPAR_DN_1(IB,O1,N)   * LP_EMULT_DN(LUM,N,NV,IB,Q)
               L_LAYERSOURCE(IB,O1,Q) = L_LAYERSOURCE(IB,O1,Q) + SPAR
             ENDDO
           ENDDO
          ENDIF

!  Other cases when N > NV
!  -----------------------

        ELSE IF ( N.GT.NV ) THEN

!  add particular solution

          DO O1 = 1, NSTOKES
            DO Q = 1, NV_PARAMETERS
              SPAR = LP_UPAR_DN_2(IB,O1,N,NV,Q) *   EMULT_DN(LUM,N,IB) + &
                        UPAR_DN_2(IB,O1,N)      * LP_EMULT_DN(LUM,N,NV,IB,Q)
              L_LAYERSOURCE(IB,O1,Q) = L_LAYERSOURCE(IB,O1,Q) + SPAR
            ENDDO
          ENDDO

!  Add single scatter term if flagged

          IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
            DO O1 = 1, NSTOKES
              DO Q = 1, NV_PARAMETERS
                SPAR =  UPAR_DN_1(IB,O1,N) * LP_EMULT_DN(LUM,N,NV,IB,Q)
                L_LAYERSOURCE(IB,O1,Q) = L_LAYERSOURCE(IB,O1,Q) + SPAR
              ENDDO
            ENDDO
          ENDIF

!  End layer clause

        ENDIF

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LP_WHOLELAYER_STERM_DN

!

      SUBROUTINE LP_PARTLAYER_STERM_UP ( &
           DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,   & ! Input flags
           DO_OBSERVATION_GEOMETRY, DO_MSMODE_VLIDORT, SOURCETERM_FLAG,     & ! Input flags
           GIVEN_LAYER, OFFGRID_INDEX, IBEAM, LAYER_TO_VARY, NV_PARAMETERS, & ! Input numbers
           NSTOKES, N_USER_STREAMS, LOCAL_UM_START,                         & ! Input numbers
           K_REAL, K_COMPLEX, LCON, MCON, UHOM_UPDN, UHOM_UPUP, UPAR_UP_1,  & ! Input user-solutions
           UPAR_UP_2, UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP,                & ! Input multipliers
           NCON, PCON, L_UHOM_UPDN, L_UHOM_UPUP, L_UPAR_UP_1, LP_UPAR_UP_2, & ! Input linearized solutions
           L_UT_HMULT_UU, L_UT_HMULT_UD, LP_UT_EMULT_UP, L_LAYER_TSUP_UTUP, & ! Linearized multipliers and thermal
           L_LAYERSOURCE )                                                    ! OUTPUT

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXSTREAMS, MAXBEAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS, &
                                 MAX_USER_STREAMS, MAX_PARTLAYERS, MAXEVALUES, MAXSTRMSTKS, ZERO, ONE, PI4

      IMPLICIT NONE

!  Input flags

      LOGICAL, INTENT (IN) ::           DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::           DO_THERMAL_TRANSONLY

      LOGICAL, INTENT (IN) ::           DO_OBSERVATION_GEOMETRY
      LOGICAL, INTENT (IN) ::           DO_MSMODE_VLIDORT
      LOGICAL, INTENT (IN) ::           SOURCETERM_FLAG

!  numbers

      INTEGER, INTENT (IN) ::           GIVEN_LAYER, OFFGRID_INDEX, IBEAM
      INTEGER, INTENT (IN) ::           NSTOKES, N_USER_STREAMS, LOCAL_UM_START
      INTEGER, INTENT (IN) ::           LAYER_TO_VARY, NV_PARAMETERS

!  RTE solutions

      INTEGER, INTENT (IN) ::           K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  LCON    ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  MCON    ( MAXSTRMSTKS, MAXLAYERS )

!  User stream RTE solutions

      DOUBLE PRECISION, INTENT (IN) ::  UHOM_UPDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UHOM_UPUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) ::  UPAR_UP_1 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UPAR_UP_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )

!  Multipliers

      DOUBLE PRECISION, INTENT (IN) ::  UT_HMULT_UD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_HMULT_UU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )

!  Linearized inputs

      DOUBLE PRECISION, INTENT (IN) :: L_LAYER_TSUP_UTUP   ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Linearized solutions

      DOUBLE PRECISION, INTENT (IN) :: NCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UHOM_UPDN  ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UHOM_UPUP  ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UPAR_UP_1  ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_UPAR_UP_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized Multipliers

      DOUBLE PRECISION, INTENT (IN) :: L_UT_HMULT_UU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UT_HMULT_UD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized output

      DOUBLE PRECISION, INTENT (OUT) :: L_LAYERSOURCE ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

!  local variables
!  ---------------

      INTEGER ::          N, NV, UM, O1, Q, IB, UT
      INTEGER ::          K, KO1, K0, K1, K2, LUM
      INTEGER ::          UM_START, UM_END
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, SP, H1, H2, TM
      DOUBLE PRECISION :: LUXR, MUXR, NUXR, PUXR, LLUXR, MLUXR
      DOUBLE PRECISION :: LUXR1, MUXR1, NUXR1, PUXR1, LLUXR1, MLUXR1
      DOUBLE PRECISION :: LUXR2, MUXR2, NUXR2, PUXR2, LLUXR2, MLUXR2

!  Proxies

      N   = GIVEN_LAYER
      KO1 = K_REAL(N) + 1
      NV  = LAYER_TO_VARY
      UT  = OFFGRID_INDEX
      IB  = IBEAM
      LUM = 1

!  Local Range

      IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
        UM_START = LOCAL_UM_START
        UM_END   = N_USER_STREAMS
      ELSE
        UM_START = IB
        UM_END   = IB
      ENDIF

!   Very important to zero both output terms (bug solved 12/29/05)

      DO Q = 1, NV_PARAMETERS
        DO UM = UM_START,UM_END
          DO O1 = 1, NSTOKES
            L_LAYERSOURCE(UM,O1,Q) = ZERO
          ENDDO
        ENDDO
      ENDDO

!  return if no source term
!      Need to go on if thermal transmittances only

      IF ( .NOT. SOURCETERM_FLAG .and. .NOT. DO_THERMAL_TRANSONLY ) RETURN

!  Avoid this section if thermal transmittance only
!  Version 2.8. remove GOTO Statement
!      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789

      IF ( .not. DO_THERMAL_TRANSONLY ) THEN

!  Partial layer source function ( Homogeneous/constants variation )
!  =================================================================

!  Special case when N = NV
!  ------------------------

        IF ( N.EQ.NV ) THEN

!  Loop over user angles and STokes

          DO UM = UM_START,UM_END
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
                  K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1

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

        ELSE IF ( N.NE.NV ) THEN

!  Loop over user angles and STokes

          DO UM = UM_START,UM_END
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
                  K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1

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
! 6789 continue

!  End scattering clause

      ENDIF

!  Add thermal emission term (direct and diffuse)
!     ----- Modulus 1.0 if solar sources are included (taken care of earlier)
!     ----- Linearization only exists if N = K

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        O1 = 1 ; TM = ONE
        IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
        IF ( N.EQ.NV ) THEN
          DO UM = UM_START,UM_END
            DO Q = 1, NV_PARAMETERS
              L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + L_LAYER_TSUP_UTUP(UM,UT,Q)*TM
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Particular and single scatter contributions
!  ===========================================

!  Lattice calculation

      IF (.not. DO_OBSERVATION_GEOMETRY ) THEN

!  Special case when N = NV
!  ------------------------

        IF ( N.EQ.NV ) THEN

!  add particular solution

          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              DO Q = 1, NV_PARAMETERS
                SP = LP_UPAR_UP_2(UM,O1,N,NV,Q) *   UT_EMULT_UP(UM,UT,IB) + &
                        UPAR_UP_2(UM,O1,N)      * LP_UT_EMULT_UP(UM,UT,NV,IB,Q)
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
                         UPAR_UP_1(UM,O1,N)   * LP_UT_EMULT_UP(UM,UT,NV,IB,Q)
                  L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SP
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!  Other cases when N > NV
!  -----------------------

        ELSE IF ( N.GT.NV ) THEN

!  add particular solution

          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              DO Q = 1, NV_PARAMETERS
               SP = LP_UPAR_UP_2(UM,O1,N,NV,Q) *   UT_EMULT_UP(UM,UT,IB) + &
                       UPAR_UP_2(UM,O1,N)      * LP_UT_EMULT_UP(UM,UT,NV,IB,Q)
                L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SP
              ENDDO
            ENDDO
          ENDDO

!  Add single scatter term if flagged

          IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                DO Q = 1, NV_PARAMETERS
                  SP =  UPAR_UP_1(UM,O1,N) * LP_UT_EMULT_UP(UM,UT,NV,IB,Q)
                  L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SP
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!  End layer clause

        ENDIF

!  Observational geometry calculation

      ELSE

!  Special case when N = NV
!  ------------------------

        IF ( N.EQ.NV ) THEN

!  add particular solution

          DO O1 = 1, NSTOKES
            DO Q = 1, NV_PARAMETERS
              SP = LP_UPAR_UP_2(IB,O1,N,NV,Q) *   UT_EMULT_UP(LUM,UT,IB) + &
                      UPAR_UP_2(IB,O1,N)      * LP_UT_EMULT_UP(LUM,UT,NV,IB,Q)
              L_LAYERSOURCE(IB,O1,Q) = L_LAYERSOURCE(IB,O1,Q) + SP
            ENDDO
          ENDDO

!  Add single scatter term if flagged

          IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
            DO O1 = 1, NSTOKES
              DO Q = 1, NV_PARAMETERS
                SP = L_UPAR_UP_1(IB,O1,N,Q) *   UT_EMULT_UP(LUM,UT,IB) + &
                       UPAR_UP_1(IB,O1,N)   * LP_UT_EMULT_UP(LUM,UT,NV,IB,Q)
                L_LAYERSOURCE(IB,O1,Q) = L_LAYERSOURCE(IB,O1,Q) + SP
              ENDDO
            ENDDO
          ENDIF

!  Other cases when N > NV
!  -----------------------

        ELSE IF ( N.GT.NV ) THEN

!  add particular solution

          DO O1 = 1, NSTOKES
            DO Q = 1, NV_PARAMETERS
              SP = LP_UPAR_UP_2(IB,O1,N,NV,Q) *   UT_EMULT_UP(LUM,UT,IB) + &
                      UPAR_UP_2(IB,O1,N)      * LP_UT_EMULT_UP(LUM,UT,NV,IB,Q)
              L_LAYERSOURCE(IB,O1,Q) = L_LAYERSOURCE(IB,O1,Q) + SP
            ENDDO
          ENDDO

!  Add single scatter term if flagged

          IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
            DO O1 = 1, NSTOKES
              DO Q = 1, NV_PARAMETERS
                SP =  UPAR_UP_1(IB,O1,N) * LP_UT_EMULT_UP(LUM,UT,NV,IB,Q)
                L_LAYERSOURCE(IB,O1,Q) = L_LAYERSOURCE(IB,O1,Q) + SP
              ENDDO
            ENDDO
          ENDIF

!  End layer clause

        ENDIF

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LP_PARTLAYER_STERM_UP

!

      SUBROUTINE LP_PARTLAYER_STERM_DN ( &
           DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,   & ! Input flags
           DO_OBSERVATION_GEOMETRY, DO_MSMODE_VLIDORT, SOURCETERM_FLAG,     & ! Input flags
           GIVEN_LAYER, OFFGRID_INDEX, IBEAM, LAYER_TO_VARY, NV_PARAMETERS, & ! Input numbers
           NSTOKES, N_USER_STREAMS, LOCAL_UM_START,                         & ! Input numbers
           K_REAL, K_COMPLEX, LCON, MCON, UHOM_DNDN, UHOM_DNUP, UPAR_DN_1,  & ! Input user-solutions
           UPAR_DN_2, UT_HMULT_DD, UT_HMULT_DU, UT_EMULT_DN,                & ! Input multipliers
           NCON, PCON, L_UHOM_DNDN, L_UHOM_DNUP, L_UPAR_DN_1, LP_UPAR_DN_2, & ! Input linearized solutions
           L_UT_HMULT_DU, L_UT_HMULT_DD, LP_UT_EMULT_DN, L_LAYER_TSUP_UTDN, & ! Input linearized Multipliers
           L_LAYERSOURCE )                                                    ! OUTPUT

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXSTREAMS, MAXBEAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS, &
                                 MAX_USER_STREAMS, MAX_PARTLAYERS, MAXEVALUES, MAXSTRMSTKS, ZERO, ONE, PI4

      IMPLICIT NONE

!  Input flags

      LOGICAL, INTENT (IN) ::           DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::           DO_THERMAL_TRANSONLY

      LOGICAL, INTENT (IN) ::           DO_OBSERVATION_GEOMETRY
      LOGICAL, INTENT (IN) ::           DO_MSMODE_VLIDORT
      LOGICAL, INTENT (IN) ::           SOURCETERM_FLAG

!  numbers

      INTEGER, INTENT (IN) ::           GIVEN_LAYER, OFFGRID_INDEX, IBEAM
      INTEGER, INTENT (IN) ::           NSTOKES, N_USER_STREAMS, LOCAL_UM_START
      INTEGER, INTENT (IN) ::           LAYER_TO_VARY, NV_PARAMETERS

!  RTE solutions

      INTEGER, INTENT (IN) ::           K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  LCON    ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  MCON    ( MAXSTRMSTKS, MAXLAYERS )

!  User stream RTE solutions

      DOUBLE PRECISION, INTENT (IN) ::  UHOM_DNDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UHOM_DNUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) ::  UPAR_DN_1 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UPAR_DN_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )

!  Multipliers

      DOUBLE PRECISION, INTENT (IN) ::  UT_HMULT_DD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_HMULT_DU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )

!  Linearized inputs

      DOUBLE PRECISION, INTENT (IN) :: L_LAYER_TSUP_UTDN   ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Linearized solutions

      DOUBLE PRECISION, INTENT (IN) :: NCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UHOM_DNDN  ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UHOM_DNUP  ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UPAR_DN_1  ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_UPAR_DN_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized Multipliers

      DOUBLE PRECISION, INTENT (IN) :: L_UT_HMULT_DU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_UT_HMULT_DD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized output

      DOUBLE PRECISION, INTENT (OUT) :: L_LAYERSOURCE ( MAX_USER_STREAMS, MAXSTOKES, MAX_ATMOSWFS )

!  local variables
!  ---------------

      INTEGER ::          N, NV, UM, O1, Q, IB, UT
      INTEGER ::          K, KO1, K0, K1, K2, LUM
      INTEGER ::          UM_START, UM_END
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, SP, H1, H2, TM
      DOUBLE PRECISION :: LUXR, MUXR, NUXR, PUXR, LLUXR, MLUXR
      DOUBLE PRECISION :: LUXR1, MUXR1, NUXR1, PUXR1, LLUXR1, MLUXR1
      DOUBLE PRECISION :: LUXR2, MUXR2, NUXR2, PUXR2, LLUXR2, MLUXR2

!  Proxies

      N   = GIVEN_LAYER
      KO1 = K_REAL(N) + 1
      NV  = LAYER_TO_VARY
      UT  = OFFGRID_INDEX
      IB  = IBEAM
      LUM = 1

!  Local Range

      IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
        UM_START = LOCAL_UM_START
        UM_END   = N_USER_STREAMS
      ELSE
        UM_START = IB
        UM_END   = IB
      ENDIF

!   Very important to zero both output terms (bug solved 12/29/05)

      DO Q = 1, NV_PARAMETERS
        DO UM = UM_START,UM_END
          DO O1 = 1, NSTOKES
            L_LAYERSOURCE(UM,O1,Q) = ZERO
          ENDDO
        ENDDO
      ENDDO

!  return if no source term
!      Need to go on if thermal transmittances only

      IF ( .NOT. SOURCETERM_FLAG .and. .NOT. DO_THERMAL_TRANSONLY ) RETURN

!  Avoid this section if thermal transmittance only
!  Version 2.8, remove GOTO statement
!      IF ( DO_THERMAL_TRANSONLY ) GO TO 6789

      IF ( .not. DO_THERMAL_TRANSONLY ) THEN

!  Partial layer source function ( Homogeneous/constants variation )
!  =================================================================

!  Special case when N = NV
!  ------------------------

        IF ( N.EQ.NV ) THEN

!  Loop over user angles and STokes

          DO UM = UM_START,UM_END
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
                  K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1

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

        ELSE IF ( N.NE.NV ) THEN

!  Loop over user angles and STokes

          DO UM = UM_START,UM_END
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
                  K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1

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
! 6789 continue

!  End scattering section

      ENDIF

!  Add thermal emission term (direct and diffuse)
!     ----- Modulus 1.0 if solar sources are included (taken care of earlier)
!     ----- Linearization only exists if N = K

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        O1 = 1 ; TM = ONE
        IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
        IF ( N.EQ.NV ) THEN
          DO UM = UM_START,UM_END
            DO Q = 1, NV_PARAMETERS
              L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + L_LAYER_TSUP_UTDN(UM,UT,Q)*TM
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Particular and single scatter contributions
!  ===========================================

!  Lattice calculation

      IF (.not. DO_OBSERVATION_GEOMETRY ) THEN

!  Special case when N = NV
!  ------------------------

        IF ( N.EQ.NV ) THEN

!  add particular solution

          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              DO Q = 1, NV_PARAMETERS
                SP = LP_UPAR_DN_2(UM,O1,N,NV,Q) *   UT_EMULT_DN(UM,UT,IB) + &
                        UPAR_DN_2(UM,O1,N)      * LP_UT_EMULT_DN(UM,UT,NV,IB,Q)
                L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SP
              ENDDO
            ENDDO
          ENDDO

!  Add single scatter term if flagged

          IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                DO Q = 1, NV_PARAMETERS
                  SP = L_UPAR_DN_1(UM,O1,N,Q) *   UT_EMULT_DN(UM,UT,IB) + &
                         UPAR_DN_1(UM,O1,N)   * LP_UT_EMULT_DN(UM,UT,NV,IB,Q)
                  L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SP
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!  Other cases when N > NV
!  -----------------------

        ELSE IF ( N.GT.NV ) THEN

!  add particular solution

          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              DO Q = 1, NV_PARAMETERS
                SP = LP_UPAR_DN_2(UM,O1,N,NV,Q) *   UT_EMULT_DN(UM,UT,IB) + &
                        UPAR_DN_2(UM,O1,N)      * LP_UT_EMULT_DN(UM,UT,NV,IB,Q)
                L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SP
              ENDDO
            ENDDO
          ENDDO

!  Add single scatter term if flagged

          IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                DO Q = 1, NV_PARAMETERS
                  SP =  UPAR_DN_1(UM,O1,N) * LP_UT_EMULT_DN(UM,UT,NV,IB,Q)
                  L_LAYERSOURCE(UM,O1,Q) = L_LAYERSOURCE(UM,O1,Q) + SP
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!  End layer clause

        ENDIF

!  Observational geometry calculation

      ELSE

!  Special case when N = NV
!  ------------------------

        IF ( N.EQ.NV ) THEN

!  add particular solution

          DO O1 = 1, NSTOKES
            DO Q = 1, NV_PARAMETERS
              SP = LP_UPAR_DN_2(IB,O1,N,NV,Q) *   UT_EMULT_DN(LUM,UT,IB) + &
                      UPAR_DN_2(IB,O1,N)      * LP_UT_EMULT_DN(LUM,UT,NV,IB,Q)
              L_LAYERSOURCE(IB,O1,Q) = L_LAYERSOURCE(IB,O1,Q) + SP
            ENDDO
          ENDDO

!  Add single scatter term if flagged

          IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
           DO O1 = 1, NSTOKES
             DO Q = 1, NV_PARAMETERS
               SP = L_UPAR_DN_1(IB,O1,N,Q) *   UT_EMULT_DN(LUM,UT,IB) + &
                      UPAR_DN_1(IB,O1,N)   * LP_UT_EMULT_DN(LUM,UT,NV,IB,Q)
               L_LAYERSOURCE(IB,O1,Q) = L_LAYERSOURCE(IB,O1,Q) + SP
             ENDDO
           ENDDO
          ENDIF

!  Other cases when N > NV
!  -----------------------

        ELSE IF ( N.GT.NV ) THEN

!  add particular solution

          DO O1 = 1, NSTOKES
            DO Q = 1, NV_PARAMETERS
              SP = LP_UPAR_DN_2(IB,O1,N,NV,Q) *   UT_EMULT_DN(LUM,UT,IB) + &
                      UPAR_DN_2(IB,O1,N)      * LP_UT_EMULT_DN(LUM,UT,NV,IB,Q)
              L_LAYERSOURCE(IB,O1,Q) = L_LAYERSOURCE(IB,O1,Q) + SP
            ENDDO
          ENDDO

!  Add single scatter term if flagged

          IF ( .NOT. DO_MSMODE_VLIDORT ) THEN
            DO O1 = 1, NSTOKES
              DO Q = 1, NV_PARAMETERS
                SP =  UPAR_DN_1(IB,O1,N) * LP_UT_EMULT_DN(LUM,UT,NV,IB,Q)
                L_LAYERSOURCE(IB,O1,Q) = L_LAYERSOURCE(IB,O1,Q) + SP
              ENDDO
            ENDDO
          ENDIF

!  End layer clause

        ENDIF

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LP_PARTLAYER_STERM_DN

!

      SUBROUTINE MIFLUX_PROFILEWF ( &
        DO_INCLUDE_MVOUTPUT, DO_INCLUDE_DIRECTBEAM,                    & ! Input flags
        DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY, & ! Input flags
        IBEAM, NV, NV_PARAMETERS, NSTOKES, NSTREAMS, NLAYERS,          & ! Input numbers
        N_DIRECTIONS, WHICH_DIRECTIONS,                                & ! Input directions
        N_USER_LEVELS, UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,         & ! Level output control
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,  & ! Partial layer output control
        FLUX_MULTIPLIER, FLUX_FACTOR, FLUXVEC, LOCAL_CSZA,             & ! Input Flux
        QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS,                      & ! Input quadrature
        INITIAL_TRANS, BEAM_CUTOFF, T_DELT_MUBAR, T_UTDN_MUBAR,        & ! Beam parameterization
        T_DELT_DISORDS, T_DISORDS_UTDN, T_DISORDS_UTUP,                & ! Discrete Ordinate Transmittances
        T_DELT_EIGEN,  T_UTDN_EIGEN, T_UTUP_EIGEN,                     & ! Homog solution    Transmittances
        K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, WUPPER,               & ! Discrete Ordinate Solutons
        LCON, MCON, T_WLOWER, T_WUPPER, BOA_THTONLY_SOURCE,            & ! Thermal solutions
        LP_INITIAL_TRANS, LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR,            & ! Lin. Beam parameterization
        LP_LEVELS_SOLARTRANS, LP_PARTIALS_SOLARTRANS,                  & ! Lin. Beam transmittances
        L_T_DELT_DISORDS, L_T_DISORDS_UTDN,  L_T_DISORDS_UTUP,         & ! Lin. Discrete Ordinate Transmittances
        L_T_DELT_EIGEN, L_T_UTDN_EIGEN, L_T_UTUP_EIGEN,                & ! Lin. Homog solution    Transmittances
        L_SOLA_XPOS, L_SOLB_XNEG, L_WLOWER, L_WUPPER, NCON, PCON,      & ! Lin. Discrete Ordinate Solutons
        L_UT_T_PARTIC, L_T_WLOWER, L_T_WUPPER, L_BOA_THTONLY_SOURCE,   & ! Lin. Thermal solutions
        MEANST_DIFFUSE_PROFWF, DNMEANST_DIRECT_PROFWF,                 & ! OUTPUT Actinic Fluxes linearized
        FLUX_DIFFUSE_PROFWF, DNFLUX_DIRECT_PROFWF )                      ! OUTPUT Regular Fluxes linearized

!  Quadrature output at offgrid or ongrid optical depths
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

! @@@ Rob fix 1/31/11, - added FLUX_FACTOR argument
!mick mod 9/19/2017 - output flux variables renamed: distinguish between diffuse and direct
!mick fix 9/19/2017 - added LP_LEVELS_SOLARTRANS & LP_PARTIALS_SOLARTRANS to facilitate correction
!                     of linearized direct flux   

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXSTREAMS, MAXBEAMS, MAXSTOKES, MAXLAYERS, MAX_PARTLAYERS, &
                                 MAX_USER_LEVELS, MAX_SZANGLES, MAX_DIRECTIONS, MAX_ATMOSWFS,            &
                                 MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, ZERO, HALF, PI2, PI4, UPIDX, DNIDX

      IMPLICIT NONE

! @@@ Rob fix 1/31/11, - added FLUX_FACTOR argument 

!  Input flags

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_MVOUTPUT
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_DIRECTBEAM

      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY

!  Input basic numbers

      INTEGER, INTENT (IN) ::          IBEAM, NV, NV_PARAMETERS
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS

!  input directions

      INTEGER, INTENT (IN) ::          N_DIRECTIONS
      INTEGER, INTENT (IN) ::          WHICH_DIRECTIONS ( MAX_DIRECTIONS )

!  Level output control

      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )

!  Partial output control

      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

!  Flux

      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER
      DOUBLE PRECISION, INTENT (IN) :: FLUX_FACTOR
      DOUBLE PRECISION, INTENT (IN) :: FLUXVEC ( MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: LOCAL_CSZA ( 0:MAXLAYERS, MAXBEAMS )

!  Quadrature

      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )

!  Beam parameterization

      INTEGER, INTENT (IN) ::          BEAM_CUTOFF    ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )

!  Discrete-ordinate transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTDN ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS )

!  Eigenstream transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN  ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN  ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN  ( MAXEVALUES, MAX_PARTLAYERS )

!  Discrete Ordinate solutions

      INTEGER, INTENT (IN) ::          K_REAL     ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX  ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS  ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG  ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: WUPPER     ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON       ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON       ( MAXSTRMSTKS, MAXLAYERS )

!  Thermal solutions

      DOUBLE PRECISION, INTENT (IN) :: T_WLOWER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_WUPPER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: BOA_THTONLY_SOURCE ( MAXSTREAMS, MAXSTOKES )

!  Linearized Input
!  ================

!  Beam parameterization

      DOUBLE PRECISION, INTENT (IN) :: LP_INITIAL_TRANS ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_T_DELT_MUBAR  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_T_UTDN_MUBAR  ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Derived Solar-beam Linearized Transmittance at all levels

      DOUBLE PRECISION, INTENT (IN) :: LP_LEVELS_SOLARTRANS ( 0:MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_PARTIALS_SOLARTRANS ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Discrete-ordinate transmittances

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_DISORDS  ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DISORDS_UTDN  ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DISORDS_UTUP  ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Eigenstream transmittances

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTUP_EIGEN ( MAXEVALUES, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_EIGEN ( MAXEVALUES, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Discrete Ordinate solutions

      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_WLOWER    ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_WUPPER    ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: NCON  ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PCON  ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )

!  Thermal solutions

      DOUBLE PRECISION, INTENT (IN) :: L_UT_T_PARTIC ( MAXSTREAMS_2, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_WLOWER    ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_WUPPER    ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_BOA_THTONLY_SOURCE  ( MAXSTREAMS, MAX_ATMOSWFS )

!  Linearized output (Has to be INOUT)

      DOUBLE PRECISION, INTENT (INOUT) :: MEANST_DIFFUSE_PROFWF &
             ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) :: DNMEANST_DIRECT_PROFWF &
             ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES )

      DOUBLE PRECISION, INTENT (INOUT) :: FLUX_DIFFUSE_PROFWF &
             ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) :: DNFLUX_DIRECT_PROFWF &
             ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES )

!  local variables
!  ---------------

!  Local quadrature output (for debug)

      DOUBLE PRECISION :: QATMOSWF_F ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  Help variables

      INTEGER ::          I, IDIR, WDIR, UTA, UT, Q, N, O1, NLEVEL
      DOUBLE PRECISION :: SMI, SFX, L_TRANS, L_FTRANS, L_DNDIRECT_MEANST, L_DNDIRECT_FLUX
!mick fix 9/19/2017 - added to facilitate correction of direct flux
      DOUBLE PRECISION :: HELP, L_TRANS_SCALED, L_FTRANS_SCALED, L_DNDIRECT_MEANST_SCALED, L_DNDIRECT_FLUX_SCALED

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

              CALL QUADPROFILEWF_OFFGRID_UP ( &
                DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,       & ! Input flags
                N, UTA, UT, IBEAM, NV, NV_PARAMETERS, NSTOKES, NSTREAMS, NLAYERS,    & ! Input numbers
                FLUX_MULTIPLIER, QUAD_STREAMS, T_DELT_DISORDS, T_DISORDS_UTUP,       & ! flux/quad/D.O.trans
                K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_UTUP_EIGEN, T_UTDN_EIGEN, & ! Input Homog. solutions
                T_UTDN_MUBAR, WUPPER, LCON, MCON, T_WUPPER, BOA_THTONLY_SOURCE,      & ! Input thermal and PI solutions
                L_T_DELT_DISORDS, L_T_DISORDS_UTUP, LP_T_UTDN_MUBAR,                 & ! Linearized
                L_SOLA_XPOS, L_SOLB_XNEG, L_T_UTUP_EIGEN, L_T_UTDN_EIGEN,  L_WUPPER, & ! Linearized Homog. solutions
                NCON, PCON, L_T_WUPPER, L_UT_T_PARTIC, L_BOA_THTONLY_SOURCE,         & ! Linearized thermal and PI solutions
                QATMOSWF_F )                                                           ! OUTPUT

            ELSE

              CALL QUADPROFILEWF_LEVEL_UP ( &
                DO_THERMAL_TRANSONLY, NLEVEL, UTA, NV, NV_PARAMETERS, NSTOKES,     & ! Input flag/numbers
                NSTREAMS, NLAYERS, FLUX_MULTIPLIER, QUAD_STREAMS, T_DELT_DISORDS,  & ! Input numbers/flux/quad/D.O.trans
                K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN,             & ! Input Homog. solutions
                LCON, MCON, T_WUPPER, BOA_THTONLY_SOURCE,                          & ! Input thermal and PI solutions
                L_T_DELT_DISORDS, L_SOLA_XPOS, L_SOLB_XNEG,  L_T_DELT_EIGEN,       & ! Linearized Homog. solutions
                L_WLOWER, L_WUPPER, NCON, PCON, L_T_WUPPER, L_BOA_THTONLY_SOURCE,  & ! Linearized thermal and PI solutions
                QATMOSWF_F )                                                         ! OUTPUT

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

              CALL QUADPROFILEWF_OFFGRID_DN ( &
                DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,       & ! Input flags
                N, UTA, UT, IBEAM, NV, NV_PARAMETERS, NSTOKES, NSTREAMS,             & ! Input numbers
                FLUX_MULTIPLIER, QUAD_STREAMS, T_DELT_DISORDS, T_DISORDS_UTDN,       & ! flux/quad/D.O.trans
                K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_UTUP_EIGEN, T_UTDN_EIGEN, & ! Input Homog. solutions
                T_UTDN_MUBAR, WUPPER, LCON, MCON, T_WUPPER, T_WLOWER,                & ! Input thermal and PI solutions
                L_T_DELT_DISORDS, L_T_DISORDS_UTDN, LP_T_UTDN_MUBAR,                 & ! Linearized
                L_SOLA_XPOS, L_SOLB_XNEG, L_T_UTUP_EIGEN, L_T_UTDN_EIGEN, L_WUPPER,  & ! Linearized Homog. solutions
                NCON, PCON, L_T_WUPPER, L_T_WLOWER, L_UT_T_PARTIC,                   & ! Linearized thermal and PI solutions
                QATMOSWF_F )                                                           ! OUTPUT

            ELSE

              CALL QUADPROFILEWF_LEVEL_DN ( &
                DO_THERMAL_TRANSONLY, NLEVEL, UTA, NV, NV_PARAMETERS, NSTOKES,     & ! Input flag/numbers
                NSTREAMS, FLUX_MULTIPLIER, QUAD_STREAMS, T_DELT_DISORDS,           & ! Input numbers/flux/quad/D.O.trans
                K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, LCON, MCON, & ! Input Homog. solutions
                T_WLOWER, L_T_DELT_DISORDS, L_SOLA_XPOS, L_SOLB_XNEG,              & ! Linearized Homog. solutions
                L_T_DELT_EIGEN, L_WLOWER, NCON, PCON, L_T_WLOWER,                  & ! Linearized thermal and PI solutions
                QATMOSWF_F )                                                         ! OUTPUT

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
                MEANST_DIFFUSE_PROFWF(Q,NV,UTA,IBEAM,O1,WDIR) = SMI * HALF
                FLUX_DIFFUSE_PROFWF(Q,NV,UTA,IBEAM,O1,WDIR)   = SFX * PI2
              ENDDO
            ENDDO
          ENDDO

!  nothing further to do if no solar sources
!   Version 2.8, remove GOTO statment, 7/19/16
!          IF ( .NOT. DO_INCLUDE_DIRECTBEAM ) GO TO 455

          IF ( DO_INCLUDE_DIRECTBEAM ) THEN

!  For the downward direction, add the direct beam contributions
!mick fix 9/19/2017 - use unscaled optical thicknesses for calculation of linearized
!                     direct flux (following LIDORT)

            IF ( WDIR .EQ. DNIDX ) THEN

!  loop over all the output optical depths

              DO UTA = 1, N_USER_LEVELS

!  For the offgrid values

                IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
                  UT = PARTLAYERS_OUTINDEX(UTA)
                  N  = PARTLAYERS_LAYERIDX(UT)

!  For the offgrid values.......
!     .....Only contributions for layers above the PI cutoff
!    LP_INITIAL_TRANS is a logarithmic derivative

!                  IF ( N .LE. BEAM_CUTOFF(IBEAM) ) THEN
!                    IF ( NV.LE.N ) THEN
!                      DO Q = 1, NV_PARAMETERS
!                        L_TRANS = LP_T_UTDN_MUBAR(UT,NV,IBEAM,Q) + &
!                                  LP_INITIAL_TRANS(N,NV,IBEAM,Q) * T_UTDN_MUBAR(UT,IBEAM)
!                        L_TRANS = L_TRANS * INITIAL_TRANS(N,IBEAM)
!                        DO O1 = 1, NSTOKES
!
!!  @@@ Rob fix (Forgot the Flux factor)
!!                  FTRANS = FLUXVEC(O1) * L_TRANS
!!                  FTRANS = FLUX_FACTOR * FLUXVEC(O1) * L_TRANS
!!mick fix - changed FLUX_FACTOR to FLUX_MULTIPLIER
!!                   FTRANS = FLUX_MULTIPLIER * FLUXVEC(O1) * L_TRANS
!! @@@ Rob fix 1/31/11, - changed FLUX_MULTIPLIER to FLUX_FACTOR
!
!                          FTRANS = FLUX_FACTOR * FLUXVEC(O1) * L_TRANS
!                          L_DIRECT_MEANI = FTRANS / PI4
!                          L_DIRECT_FLUX  = FTRANS * LOCAL_CSZA(N,IBEAM)
!                          MINT_PROFILEWF_DIRECT(Q,NV,UTA,IBEAM,O1)=L_DIRECT_MEANI
!                          FLUX_PROFILEWF_DIRECT(Q,NV,UTA,IBEAM,O1)=L_DIRECT_FLUX
!                          MINT_PROFILEWF(Q,NV,UTA,IBEAM,O1,WDIR) = &
!                                  MINT_PROFILEWF(Q,NV,UTA,IBEAM,O1,WDIR) + L_DIRECT_MEANI
!                          FLUX_PROFILEWF(Q,NV,UTA,IBEAM,O1,WDIR) = &
!                                  FLUX_PROFILEWF(Q,NV,UTA,IBEAM,O1,WDIR) + L_DIRECT_FLUX
!                        ENDDO
!                      ENDDO
!                    ENDIF
!                  ENDIF

                  IF ( N .LE. BEAM_CUTOFF(IBEAM) ) THEN
                    IF ( NV.LE.N ) THEN
                      DO Q = 1, NV_PARAMETERS

!  Linearized direct transmittances, scaled and unscaled

                        L_TRANS = LP_PARTIALS_SOLARTRANS(UT,NV,IBEAM,Q)
                        HELP    = LP_T_UTDN_MUBAR(UT,NV,IBEAM,Q) + &
                                  LP_INITIAL_TRANS(N,NV,IBEAM,Q) * T_UTDN_MUBAR(UT,IBEAM)
                        L_TRANS_SCALED = HELP * INITIAL_TRANS(N,IBEAM)

                        DO O1 = 1, NSTOKES

!  Direct calculation with non-scaled linearized transmittances

                          L_FTRANS = FLUX_FACTOR * FLUXVEC(O1) * L_TRANS
                          L_DNDIRECT_MEANST = L_FTRANS / PI4
                          L_DNDIRECT_FLUX   = L_FTRANS * LOCAL_CSZA(N,IBEAM)
                          DNMEANST_DIRECT_PROFWF(Q,NV,UTA,IBEAM,O1) = L_DNDIRECT_MEANST
                          DNFLUX_DIRECT_PROFWF(Q,NV,UTA,IBEAM,O1)   = L_DNDIRECT_FLUX

!  Diffuse calculation

                          L_FTRANS_SCALED = FLUX_FACTOR * FLUXVEC(O1) * L_TRANS_SCALED
                          L_DNDIRECT_MEANST_SCALED = L_FTRANS_SCALED / PI4
                          L_DNDIRECT_FLUX_SCALED   = L_FTRANS_SCALED * LOCAL_CSZA(N,IBEAM)
                          MEANST_DIFFUSE_PROFWF(Q,NV,UTA,IBEAM,O1,WDIR) = &
                                MEANST_DIFFUSE_PROFWF(Q,NV,UTA,IBEAM,O1,WDIR) + ( L_DNDIRECT_MEANST_SCALED - L_DNDIRECT_MEANST )
                          FLUX_DIFFUSE_PROFWF(Q,NV,UTA,IBEAM,O1,WDIR) = &
                                FLUX_DIFFUSE_PROFWF(Q,NV,UTA,IBEAM,O1,WDIR)   + ( L_DNDIRECT_FLUX_SCALED - L_DNDIRECT_FLUX )
                        ENDDO
                      ENDDO
                    ENDIF
                  ENDIF

!  For the on-grid values
!    LP_INITIAL_TRANS is a logarithmic derivative

                ELSE
                  N = UTAU_LEVEL_MASK_DN(UTA)

!                  IF ( N .LE. BEAM_CUTOFF(IBEAM) ) THEN
!                    IF ( N.GT.0 ) THEN
!                      IF ( NV.LE.N ) THEN
!                        DO Q = 1, NV_PARAMETERS
!                          L_TRANS = LP_T_DELT_MUBAR(N,NV,IBEAM,Q) + &
!                                    LP_INITIAL_TRANS(N,NV,IBEAM,Q) * T_DELT_MUBAR(N,IBEAM)
!                          L_TRANS = L_TRANS * INITIAL_TRANS(N,IBEAM)
!                          DO O1 = 1, NSTOKES
!
!!  @@@ Rob fix (Forgot the Flux factor)
!!                  FTRANS = FLUXVEC(O1) * L_TRANS
!!                  FTRANS = FLUX_FACTOR * FLUXVEC(O1) * L_TRANS
!!mick fix - changed FLUX_FACTOR to FLUX_MULTIPLIER
!!                   FTRANS = FLUX_MULTIPLIER * FLUXVEC(O1) * L_TRANS
!! @@@ Rob fix 1/31/11, - changed FLUX_MULTIPLIER to FLUX_FACTOR
!
!                            FTRANS = FLUX_FACTOR * FLUXVEC(O1) * L_TRANS
!                            L_DIRECT_MEANI = FTRANS / PI4
!                            L_DIRECT_FLUX  = FTRANS * LOCAL_CSZA(N,IBEAM)
!                            MINT_PROFILEWF_DIRECT(Q,NV,UTA,IBEAM,O1)=L_DIRECT_MEANI
!                            FLUX_PROFILEWF_DIRECT(Q,NV,UTA,IBEAM,O1)=L_DIRECT_FLUX
!                            MINT_PROFILEWF(Q,NV,UTA,IBEAM,O1,WDIR) = &
!                                     MINT_PROFILEWF(Q,NV,UTA,IBEAM,O1,WDIR) + L_DIRECT_MEANI
!                            FLUX_PROFILEWF(Q,NV,UTA,IBEAM,O1,WDIR) = &
!                                     FLUX_PROFILEWF(Q,NV,UTA,IBEAM,O1,WDIR) + L_DIRECT_FLUX
!                          ENDDO
!                        ENDDO
!                      ENDIF
!                    ENDIF
!                  ENDIF

                  IF ( N .LE. BEAM_CUTOFF(IBEAM) ) THEN
                    IF ( N.GT.0 ) THEN
                      IF ( NV.LE.N ) THEN
                        DO Q = 1, NV_PARAMETERS

!  Linearized direct transmittances, scaled and unscaled

                          L_TRANS = LP_LEVELS_SOLARTRANS(N,NV,IBEAM,Q)
                          HELP    = LP_T_DELT_MUBAR(N,NV,IBEAM,Q) + &
                                    LP_INITIAL_TRANS(N,NV,IBEAM,Q) * T_DELT_MUBAR(N,IBEAM)
                          L_TRANS_SCALED = HELP * INITIAL_TRANS(N,IBEAM)

                          DO O1 = 1, NSTOKES

!  Direct calculation with non-scaled linearized transmittances

                            L_FTRANS = FLUX_FACTOR * FLUXVEC(O1) * L_TRANS
                            L_DNDIRECT_MEANST = L_FTRANS / PI4
                            L_DNDIRECT_FLUX   = L_FTRANS * LOCAL_CSZA(N,IBEAM)
                            DNMEANST_DIRECT_PROFWF(Q,NV,UTA,IBEAM,O1) = L_DNDIRECT_MEANST
                            DNFLUX_DIRECT_PROFWF(Q,NV,UTA,IBEAM,O1)   = L_DNDIRECT_FLUX

!  Diffuse calculation

                            L_FTRANS_SCALED = FLUX_FACTOR * FLUXVEC(O1) * L_TRANS_SCALED
                            L_DNDIRECT_MEANST_SCALED = L_FTRANS_SCALED / PI4
                            L_DNDIRECT_FLUX_SCALED   = L_FTRANS_SCALED * LOCAL_CSZA(N,IBEAM)
                            MEANST_DIFFUSE_PROFWF(Q,NV,UTA,IBEAM,O1,WDIR) = &
                                   MEANST_DIFFUSE_PROFWF(Q,NV,UTA,IBEAM,O1,WDIR) + ( L_DNDIRECT_MEANST_SCALED - L_DNDIRECT_MEANST )
                            FLUX_DIFFUSE_PROFWF(Q,NV,UTA,IBEAM,O1,WDIR) = &
                                   FLUX_DIFFUSE_PROFWF(Q,NV,UTA,IBEAM,O1,WDIR)   + ( L_DNDIRECT_FLUX_SCALED - L_DNDIRECT_FLUX )
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
! 455      CONTINUE. Version 2.8, GOTO REMOVED 7/19/16. Replaced by IF clause

          ENDIF

!  Finish MV output

        ENDIF

!  end direction loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE MIFLUX_PROFILEWF

!

      SUBROUTINE QUADPROFILEWF_LEVEL_UP ( &
           DO_THERMAL_TRANSONLY, NL, UTA, NV, NV_PARAMETERS, NSTOKES,         & ! Input flag/numbers
           NSTREAMS, NLAYERS, FLUX_MULTIPLIER, QUAD_STREAMS, T_DELT_DISORDS,  & ! Input numbers/flux/quad/D.O.trans
           K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN,             & ! Input Homog. solutions
           LCON, MCON, T_WUPPER, BOA_THTONLY_SOURCE,                          & ! Input thermal and PI solutions
           L_T_DELT_DISORDS, L_SOLA_XPOS, L_SOLB_XNEG,  L_T_DELT_EIGEN,       & ! Linearized Homog. solutions
           L_WLOWER, L_WUPPER, NCON, PCON, L_T_WUPPER, L_BOA_THTONLY_SOURCE,  & ! Linearized thermal and PI solutions
           QATMOSWF_F )                                                         ! OUTPUT

!  Upwelling weighting function Fourier components at level boundary NL
!  Quadrature angles only

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXLAYERS, &
                                 MAX_USER_LEVELS, MAX_SZANGLES, MAX_ATMOSWFS,  &
                                 MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, ZERO 

      IMPLICIT NONE

!  Input flags

      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY

!  Input basic numbers

      INTEGER, INTENT (IN) ::          NL, UTA, NV, NV_PARAMETERS
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS

!  Flux

      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER

!  Quadrature

      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )

!  Discrete-ordinate transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )

!  Discrete Ordinate solutions

      INTEGER, INTENT (IN) ::          K_REAL     ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX  ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS  ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG  ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON       ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON       ( MAXSTRMSTKS, MAXLAYERS )

!  Thermal solutions

      DOUBLE PRECISION, INTENT (IN) :: T_WUPPER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: BOA_THTONLY_SOURCE ( MAXSTREAMS, MAXSTOKES )

!  Linearized Input
!  ================

!  Discrete-ordinate transmittances

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_DISORDS  ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )

!  Eigenstream transmittances

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

!  Discrete Ordinate solutions

      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_WLOWER    ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_WUPPER    ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: NCON  ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PCON  ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )

!  Thermal solutions

      DOUBLE PRECISION, INTENT (IN) :: L_T_WUPPER    ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_BOA_THTONLY_SOURCE  ( MAXSTREAMS, MAX_ATMOSWFS )

!  Output
!  ------

      DOUBLE PRECISION, INTENT (INOUT) :: QATMOSWF_F ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

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

        IF ( NV .EQ. NL ) THEN

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
                  K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1

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

        ELSE IF ( NV.LT.NL ) THEN

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
                  K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1

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
              IF ( LAY.EQ.NV ) THEN
                L_TPROP = L_T_WUPPER(I1,LAY,Q) / QUAD_STREAMS(I)
                L_THELP = L_THELP + L_TPROP + THELP * L_T_DELT_DISORDS(I,LAY,Q)
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

        IF ( NV .EQ. N ) THEN

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
                  K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1

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

        ELSE IF ( NV.NE.N ) THEN

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
                  K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1

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
                IF ( NV .LT. N ) SPAR = L_WUPPER(I1,O1,N,Q)
                QATMOSWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER * ( SPAR + SHOM )

!  Finish Q, I and O1 loops

              ENDDO
            ENDDO
          ENDDO

!  End variability clauses

        ENDIF
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE QUADPROFILEWF_LEVEL_UP

!

      SUBROUTINE QUADPROFILEWF_LEVEL_DN ( &
           DO_THERMAL_TRANSONLY, NL, UTA, NV, NV_PARAMETERS, NSTOKES,         & ! Input flag/numbers
           NSTREAMS, FLUX_MULTIPLIER, QUAD_STREAMS, T_DELT_DISORDS,           & ! Input numbers/flux/quad/D.O.trans
           K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, LCON, MCON, & ! Input RT solutions
           T_WLOWER, L_T_DELT_DISORDS, L_SOLA_XPOS, L_SOLB_XNEG,              & ! Linearized Homog. solutions
           L_T_DELT_EIGEN, L_WLOWER,  NCON, PCON, L_T_WLOWER,                 & ! Linearized thermal and PI solutions
           QATMOSWF_F )                                                         ! OUTPUT

!  Upwelling weighting function Fourier components at level boundary NL
!  Quadrature angles only

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXLAYERS, &
                                 MAX_USER_LEVELS, MAX_SZANGLES, MAX_ATMOSWFS,  &
                                 MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, ZERO 

      IMPLICIT NONE

!  Input flags

      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY

!  Input basic numbers

      INTEGER, INTENT (IN) ::          NL, UTA, NV, NV_PARAMETERS
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS

!  Flux

      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER

!  Quadrature

      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )

!  Discrete-ordinate transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )

!  Discrete Ordinate solutions

      INTEGER, INTENT (IN) ::          K_REAL     ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX  ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS  ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG  ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN  ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON       ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON       ( MAXSTRMSTKS, MAXLAYERS )

!  Thermal solutions

      DOUBLE PRECISION, INTENT (IN) :: T_WLOWER ( MAXSTREAMS_2, MAXLAYERS )

!  Linearized Input
!  ================

!  Discrete-ordinate transmittances

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_DISORDS  ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )

!  Eigenstream transmittances

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

!  Discrete Ordinate solutions

      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_WLOWER    ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: NCON  ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PCON  ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )

!  Thermal solutions

      DOUBLE PRECISION, INTENT (IN) :: L_T_WLOWER    ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )

!  Output
!  ------

      DOUBLE PRECISION, INTENT (INOUT) :: QATMOSWF_F ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

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
              IF ( LAY.EQ.NV ) THEN
                L_TPROP = L_T_WLOWER(I,LAY,Q) / QUAD_STREAMS(I)
                L_THELP = L_THELP +  L_TPROP + THELP * L_T_DELT_DISORDS(I,LAY,Q)
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

        IF ( NV .EQ. N ) THEN

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
                  K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1

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
                QATMOSWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER * ( SPAR + SHOM )

!  Finish Q, I and O1 loops

              ENDDO
            ENDDO
          ENDDO

!  varying layer above or below active layer

        ELSE IF ( NV.NE.N ) THEN

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
                  K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1
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
                QATMOSWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER * ( SPAR + SHOM )

!  Finish Q, I and O1 loops

              ENDDO
            ENDDO
          ENDDO

!  Finish variability clauses

        ENDIF
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE QUADPROFILEWF_LEVEL_DN

!

      SUBROUTINE QUADPROFILEWF_OFFGRID_UP ( &
           DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,       & ! Input flags
           N, UTA, UT, IB, NV, NV_PARAMETERS, NSTOKES, NSTREAMS, NLAYERS,       & ! Input numbers
           FLUX_MULTIPLIER, QUAD_STREAMS, T_DELT_DISORDS, T_DISORDS_UTUP,       & ! flux/quad/D.O.trans
           K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_UTUP_EIGEN, T_UTDN_EIGEN, & ! Input Homog. solutions
           T_UTDN_MUBAR, WUPPER, LCON, MCON, T_WUPPER, BOA_THTONLY_SOURCE,      & ! Input thermal and PI solutions
           L_T_DELT_DISORDS, L_T_DISORDS_UTUP, LP_T_UTDN_MUBAR,                 & ! Linearized
           L_SOLA_XPOS, L_SOLB_XNEG, L_T_UTUP_EIGEN, L_T_UTDN_EIGEN, L_WUPPER,  & ! Linearized Homog. solutions
           NCON, PCON, L_T_WUPPER, L_UT_T_PARTIC, L_BOA_THTONLY_SOURCE,         & ! Linearized thermal and PI solutions
           QATMOSWF_F )                                                           ! OUTPUT

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXBEAMS, MAXLAYERS, &
                                 MAX_PARTLAYERS, MAX_USER_LEVELS, MAX_ATMOSWFS,          &
                                 MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, ZERO

      IMPLICIT NONE

!  Input flags

      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY

!  Input basic numbers

      INTEGER, INTENT (IN) ::          N, UTA, UT, IB, NV, NV_PARAMETERS
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS

!  Flux

      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER

!  Quadrature

      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )

!  Discrete-ordinate and beam transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXBEAMS )

!  Eigenstream transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN  ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN  ( MAXEVALUES, MAX_PARTLAYERS )

!  Discrete Ordinate solutions

      INTEGER, INTENT (IN) ::          K_REAL     ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX  ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS  ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG  ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: WUPPER     ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON       ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON       ( MAXSTRMSTKS, MAXLAYERS )

!  Thermal solutions

      DOUBLE PRECISION, INTENT (IN) :: T_WUPPER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: BOA_THTONLY_SOURCE ( MAXSTREAMS, MAXSTOKES )

!  Linearized Input
!  ================

!  Discrete-ordinate and beam  transmittances

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_DISORDS  ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DISORDS_UTUP  ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_T_UTDN_MUBAR  ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Eigenstream transmittances

      DOUBLE PRECISION, INTENT (IN) :: L_T_UTUP_EIGEN ( MAXEVALUES, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_EIGEN ( MAXEVALUES, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Discrete Ordinate solutions

      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_WUPPER    ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: NCON  ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PCON  ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )

!  Thermal solutions

      DOUBLE PRECISION, INTENT (IN) :: L_UT_T_PARTIC ( MAXSTREAMS_2, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_WUPPER    ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_BOA_THTONLY_SOURCE  ( MAXSTREAMS, MAX_ATMOSWFS )

!  output
!  ------

      DOUBLE PRECISION, INTENT (INOUT) :: QATMOSWF_F ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

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
              IF ( LAY.EQ.NV ) THEN
                L_TPROP = L_T_WUPPER(I1,LAY,Q) / QUAD_STREAMS(I)
                L_THELP = L_THELP + L_TPROP + THELP * L_T_DELT_DISORDS(I,LAY,Q)
              ENDIF
              TPROP = T_WUPPER(I1,Q) / QUAD_STREAMS(I)
              THELP = THELP * T_DELT_DISORDS(I,LAY) + TPROP
            ENDDO
            L_THELP = L_THELP * T_DISORDS_UTUP(I,UT)
            IF ( N.EQ.NV ) THEN
              L_TPROP = L_UT_T_PARTIC(I1,UT,Q) / QUAD_STREAMS(I)
              L_THELP = L_THELP + L_TPROP + THELP * L_T_DISORDS_UTUP(I,UT,Q)
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

      IF ( N.EQ. NV ) THEN

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
                K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1

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

      ELSE IF ( NV.NE.N ) THEN

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
                K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1

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
!    ---Only present if N = NV
!   THIS is the solution with scattering

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        O1 = 1
        IF ( N.EQ.NV ) THEN
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO Q = 1, NV_PARAMETERS
              SPAR = L_UT_T_PARTIC(I1,UT,Q)
              QATMOSWF_F(Q,UTA,I,O1) = QATMOSWF_F(Q,UTA,I,O1) + FMULT * SPAR
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  Finished if no solar terms

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Add Linearized Classical beam solution
!  Green function solution is still a Placeholder....!
!    solutions exist only for N = NV or N > NV

! @@@@@@@@@@ Rob fix, 2/9/11, Replacement of commented section
!      IF ( DO_CLASSICAL_SOLUTION ) THEN
!        IF ( N.GE.NV ) THEN
!          DO I = 1, NSTREAMS
!            I1 = I + NSTREAMS
!            DO O1 = 1, NSTOKES
!              DO Q = 1, NV_PARAMETERS
!                SPAR = L_WUPPER(I1,O1,N,Q) *   T_UTDN_MUBAR(UT,IB) + &
!                         WUPPER(I1,O1,N)   * LP_T_UTDN_MUBAR(UT,NV,IB,Q)
!                QATMOSWF_F(Q,UTA,I,O1) = &
!                  QATMOSWF_F(Q,UTA,I,O1) + FMULT * SPAR
!              ENDDO
!            ENDDO
!          ENDDO
!        ENDIF
!      ELSE
!                P L A C E H O L D E R
!      ENDIF
! @@@@@@@@@@ Rob fix, 2/9/11, End Replaced section

!  Version 2.8. Remove Classical solution flag

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        VAREX = ( N.EQ.NV )
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
                         SUNP * LP_T_UTDN_MUBAR(UT,NV,IB,Q)
                QATMOSWF_F(Q,UTA,I,O1) = QATMOSWF_F(Q,UTA,I,O1) + FMULT * SPAR
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ELSE
        IF ( N.GE.NV ) THEN
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES
              DO Q = 1, NV_PARAMETERS
                SPAR = L_WUPPER(I1,O1,N,Q) *   T_UTDN_MUBAR(UT,IB) + &
                         WUPPER(I1,O1,N)   * LP_T_UTDN_MUBAR(UT,NV,IB,Q)
                QATMOSWF_F(Q,UTA,I,O1) = QATMOSWF_F(Q,UTA,I,O1) + FMULT * SPAR
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE QUADPROFILEWF_OFFGRID_UP

!

      SUBROUTINE QUADPROFILEWF_OFFGRID_DN ( &
           DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,       & ! Input flags
           N, UTA, UT, IB, NV, NV_PARAMETERS, NSTOKES, NSTREAMS,                & ! Input numbers
           FLUX_MULTIPLIER, QUAD_STREAMS, T_DELT_DISORDS, T_DISORDS_UTDN,       & ! flux/quad/D.O.trans
           K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_UTUP_EIGEN, T_UTDN_EIGEN, & ! Input Homog. solutions
           T_UTDN_MUBAR, WUPPER, LCON, MCON, T_WUPPER, T_WLOWER,                & ! Input thermal and PI solutions
           L_T_DELT_DISORDS, L_T_DISORDS_UTDN, LP_T_UTDN_MUBAR,                 & ! Linearized
           L_SOLA_XPOS, L_SOLB_XNEG, L_T_UTUP_EIGEN, L_T_UTDN_EIGEN, L_WUPPER,  & ! Linearized Homog. solutions
           NCON, PCON, L_T_WUPPER, L_T_WLOWER, L_UT_T_PARTIC,                   & ! Linearized thermal and PI solutions
           QATMOSWF_F )                                                           ! OUTPUT

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXBEAMS, MAXLAYERS, &
                                 MAX_PARTLAYERS, MAX_USER_LEVELS, MAX_ATMOSWFS,          &
                                 MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, ZERO

      IMPLICIT NONE

!  Input flags

      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY

!  Input basic numbers

      INTEGER, INTENT (IN) ::          N, UTA, UT, IB, NV, NV_PARAMETERS
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
!  Flux

      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER

!  Quadrature

      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )

!  Discrete-ordinate and beam transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTDN ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXBEAMS )

!  Eigenstream transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN  ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN  ( MAXEVALUES, MAX_PARTLAYERS )

!  Discrete Ordinate solutions

      INTEGER, INTENT (IN) ::          K_REAL     ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX  ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS  ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG  ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: WUPPER     ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON       ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON       ( MAXSTRMSTKS, MAXLAYERS )

!  Thermal solutions

      DOUBLE PRECISION, INTENT (IN) :: T_WLOWER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_WUPPER ( MAXSTREAMS_2, MAXLAYERS )

!  Linearized Input
!  ================

!  Discrete-ordinate and beam  transmittances

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_DISORDS  ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DISORDS_UTDN  ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LP_T_UTDN_MUBAR  ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Eigenstream transmittances

      DOUBLE PRECISION, INTENT (IN) :: L_T_UTUP_EIGEN ( MAXEVALUES, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_UTDN_EIGEN ( MAXEVALUES, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Discrete Ordinate solutions

      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_WUPPER    ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: NCON  ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: PCON  ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )

!  Thermal solutions

      DOUBLE PRECISION, INTENT (IN) :: L_UT_T_PARTIC ( MAXSTREAMS_2, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_WLOWER    ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_WUPPER    ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )

!  output
!  ------

      DOUBLE PRECISION, INTENT (INOUT) :: QATMOSWF_F ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

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
              IF ( LAY.EQ.NV ) THEN
                L_TPROP = L_T_WLOWER(I,LAY,Q) / QUAD_STREAMS(I)
                L_THELP = L_THELP + L_TPROP + THELP * L_T_DELT_DISORDS(I,LAY,Q)
              ENDIF
              TPROP = T_WLOWER(I,LAY) / QUAD_STREAMS(I)
              THELP = THELP * T_DELT_DISORDS(I,LAY) + TPROP
            ENDDO
            L_THELP = L_THELP * T_DISORDS_UTDN(I,UT)
            IF ( N.EQ.NV ) THEN
              L_TPROP = L_UT_T_PARTIC(I,UT,Q) / QUAD_STREAMS(I)
              L_THELP = L_THELP + L_TPROP  + THELP * L_T_DISORDS_UTDN(I,UT,Q)
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

      IF ( N .EQ. NV ) THEN

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
                K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1

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

      ELSE IF ( NV.NE.N ) THEN

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
                K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1

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
!    ---Only present if N = NV
!   THIS IS THE SOLUTION with scattering

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        O1 = 1
        IF ( N.EQ.NV ) THEN
          DO I = 1, NSTREAMS
            DO Q = 1, NV_PARAMETERS
              SPAR = L_UT_T_PARTIC(I,UT,Q)
              QATMOSWF_F(Q,UTA,I,O1) = QATMOSWF_F(Q,UTA,I,O1) + FMULT * SPAR
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  Finished if no solar terms

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Add Linearized Classical beam solution
!  Green function solution is still a Placeholder....!
!    solutions exist only for N = NV or N > NV

! @@@@@@@@@@ Rob fix, 2/9/11, Replacement of commented section
!      IF ( DO_CLASSICAL_SOLUTION ) THEN
!        IF ( N.GE.NV ) THEN
!          DO I = 1, NSTREAMS
!            DO O1 = 1, NSTOKES
!              DO Q = 1, NV_PARAMETERS
!                SPAR = L_WUPPER(I,O1,N,Q) *   T_UTDN_MUBAR(UT,IB) + &
!                         WUPPER(I,O1,N)   * LP_T_UTDN_MUBAR(UT,NV,IB,Q)
!                QATMOSWF_F(Q,UTA,I,O1) = &
!                  QATMOSWF_F(Q,UTA,I,O1) + FMULT * SPAR
!              ENDDO
!            ENDDO
!          ENDDO
!        ENDIF
!      ELSE
!                P L A C E H O L D E R
!      ENDIF
! @@@@@@@@@@ Rob fix, 2/9/11, End Replaced section

!  Version 2.8. Remove Classical solution flag

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        VAREX = ( N.EQ.NV )
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
                         SUNP * LP_T_UTDN_MUBAR(UT,NV,IB,Q)
                QATMOSWF_F(Q,UTA,I,O1) = QATMOSWF_F(Q,UTA,I,O1) + FMULT * SPAR
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ELSE
        IF ( N.GE.NV ) THEN
          DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES
              DO Q = 1, NV_PARAMETERS
                SPAR = L_WUPPER(I,O1,N,Q) *   T_UTDN_MUBAR(UT,IB) + &
                         WUPPER(I,O1,N)   * LP_T_UTDN_MUBAR(UT,NV,IB,Q)
                QATMOSWF_F(Q,UTA,I,O1) = QATMOSWF_F(Q,UTA,I,O1) + FMULT * SPAR
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE QUADPROFILEWF_OFFGRID_DN

!  End Module

      END MODULE vlidort_lp_wfatmos_m
