
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
! #  Author :      Robert. J. D. Spurr                          #
! #                                                             #
! #  Address :     RT Solutions, Inc.                           #
! #                9 Channing Street                            #
! #                Cambridge, MA 02138, USA                     #
! #                                                             #
! #  Tel:          (617) 492 1183                               #
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #              FIRST-ORDER SCALAR/VECTOR MODEL                #
! #     (EXACT SINGLE-SCATTERING and DIRECT-THERMAL)            #
! #                                                             #
! #  This Version :   FO_1p5 (Version 1.5.2)                    #
! #  Release Date :   15 April 2020                             #
! #                                                             #
! #  Previous FO CODE Versions under Standard GPL 3.0:          #
! #  -------------------------------------------------          #
! #                                                             #
! #   Version 1.1,  13 February 2012, First Code                #
! #   Version 1.2,  01 June     2012, Modularization            #
! #   Version 1.3a, 29 October  2012, Obsgeom Multi-geom.       #
! #   Version 1.3b, 24 January  2013, BRDF/SL Supplements       #
! #   Version 1.4,  31 July     2013, Lattice Multi-geom.       #
! #   Version 1.5,  7  July     2016. Use Fmatrix/Phasfunc      #
! #   Version 1.5,  22 August   2016. Partial-layer output.     #
! #   Version 1.5,  30 April    2017. Shakedown completed       #
! #                                                             #
! #   FO Version 1.5   coincides (V)LIDORT Version (2.8)3.8     #
! #   FO Version 1.5.1 coincides (V)LIDORT Version (2.8.1)3.8.1 #
! #                                                             #
! #  This Version                                               #
! #  ------------                                               #
! #                                                             #
! #   FO_1.5.2, released April 15 2020.                         #
! #     ==> Geometry (FO), separated from RT Calculations       #
! #     ==> Use of phase functions only in FO code              #
! #     ==> Use I/O Geometry type structures directly           #
! #     ==> Optional for the Doublet Geometry                   #
! #                                                             #
! #   FO Version 1.5.2 coincides (V)LIDORT Version (2.8.2)3.8.2 #
! #                                                             #
! ###############################################################

! ##################################################################
! #                                                                #
! # This is Version 1.5.2 of the FO_1p5 software library. This     #
! # library comes with the Standard GNU General Public License,    #
! # Version 3.0, 29 June 2007. Please read this license carefully. #
! #                                                                #
! #      FO CODE Copyright (c) 2010-2020.                          #
! #          Robert Spurr, RT Solutions Inc.                       #
! #          9 Channing Street, Cambridge, MA 02138, USA.          #
! #                                                                #
! # This file is part of FO_1p5 Version 1.5.2.                     #
! #                                                                #
! # FO_1p5 is free software: you can redistribute it               #
! # and/or modify it under the terms of the GNU GPL (Standard)     #
! # General Public License) as published by the Free Software      #
! # Foundation, either version 3.0 of this License, or any         #
! # later version.                                                 #
! #                                                                #
! # FO_1p5 is distributed in the hope that it will be              #
! # useful, but WITHOUT ANY WARRANTY; without even the implied     #
! # warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR        #
! # PURPOSE. See the GNU Standard General Public License (GPL)     #
! # for more details.                                              #
! #                                                                #
! # You should have received a copy of the GNU Standard General    #
! # Public License (GPL) Version 3.0, along with FO_1p5_           #
! # Version 1.5.2. If not, see <http://www.gnu.org/licenses/>.     #
! #                                                                #
! ##################################################################

!  4/15/20. Version 2.8.2
!  ======================

!  Geometry separation.

!  FO Version history
!  ------------------

!  Versions to 1.4, without Partials. Code is stand alone with no dependencies.
!  Version     1.5, with optional phase function, Surface Leaving  and partials.

!    Version 1a, 01 December 2011, R. Spurr, RT Solutions Inc.
!    Version 1b, 13 February 2012, R. Spurr, RT Solutions Inc.
!    Version 2,  01 June     2012, R. Spurr, RT Solutions Inc.
!    Version 3,  29 October  2012, Extension to Observational multiple geometries
!    Version 4,  31 July     2013, Lattice Multi-geometry
!    Version 5,  07 July     2016, Optional F-matrix usage
!    Version 5,  02 August   2016. Surface leaving  + Jacobians
!    Version 5,  19 August   2016, Partial-layer output

!  VLIDORT Interface history
!  ------------------------

!    FO Version 1.4: This module is interface with VLIDORT V2.7. R.Spurr 3/19/15
!    FO Version 1.5, 7/7/16 and 8/2/16 and 8/20/16 ( partials)
!       - Optional calculation using F Matrices directly.
!       - inclusion of Surface-leaving terms + LSSL weighting functions.
!       - Partial-layer output introduced
!       - Introduce water-leaving flag, output cumtrans (for VLIDORT 2.8.1 )
!    FO Version 1.5.2: Interface module upgraded to Version 2.8.2. R.Spurr, 4/15/20

module VFO_Geometry_Master_m

!  Use statements

   Use VLIDORT_pars_m
   USE VLIDORT_Inputs_def_m
   USE VLIDORT_Setups_def_m
   
!  Use modules

   USE FO_SSWPGeometry_Master_m
   USE FO_DTWPGeometry_Master_m
   USE FO_VectorSS_spherfuncs_m

!  Main subroutine public

public

contains

subroutine VFO_GEOMETRY_MASTER &
     ( VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Bookkeep,     & ! Inputs and Bookkeeping
       VLIDORT_FOGeom, Gfail, message, trace_1, trace_2 )    ! Output FOGeometry and Exception-Handling
  
!         extinction_crit             & ! Input thermal/surf optical
  
   implicit none

!  Subroutine Arguments
!  ====================
   
!  Type structures, inputs and bookkeeping

   TYPE(VLIDORT_Fixed_Inputs)   , INTENT(IN)    :: VLIDORT_FixIn
   TYPE(VLIDORT_Modified_Inputs), INTENT(IN)    :: VLIDORT_ModIn
   TYPE(VLIDORT_Bookkeeping)    , INTENT(INOUT) :: VLIDORT_Bookkeep

!  Type structures, FO Geometry output

   TYPE(VLIDORT_Geometry_FO), INTENT(INOUT)    :: VLIDORT_FOGeom

!  Exception handling

   logical, intent(out)           :: Gfail
   character (len=*), intent(out) :: message
   character (len=*), intent(out) :: trace_1, trace_2

!  LOCAL VARIABLES (PROXY INPUTS)
!  ==============================
   
!  Max dimensions

   integer  :: maxgeoms, maxszas, maxvzas, maxazms, maxpartials, maxfine

!  Directional Flags

   logical  :: do_upwelling, do_dnwelling

!  Obsgeom/Doublet flags

   logical  :: do_Obsgeom, do_doublet

!  Flags (sphericity flags are mutually exclusive). Regular PS now removed, version 1.5

   logical  :: do_planpar, do_enhanced_ps

!  Numbers
!  -------

!  Layer and geometry control. Finelayer divisions may be changed

   integer :: ngeoms, nszas, nvzas, nazms, nlayers, nfine
   integer :: nmoments_input, nstokes

!  Control for partial-layer output, added 8/25/16

   logical :: do_Partials
   integer :: Npartials
   integer :: partial_layeridx(MAX_PARTLAYERS)

!  General inputs
!  --------------

!  Critical adjustment for cloud layers. Not enabled. 9/17/16

!   logical   :: doCrit
!   real(fpk) :: Acrit

!  Earth radius + heights. Partials added 9/17/16.

   real(fpk)   :: eradius, heights (0:MAXLAYERS)
   real(fpk)   :: partial_heights (MAX_PARTLAYERS)

!  Geometry inputs
!  ---------------

!  input angles (Degrees). Enough information for Lattice or Obsgeom.
!   Convention for ObsGeom = same as VLIDORT/LIDORT (1=sza,2=vza,3=azm)
!    In both cases, the Phi angle may be changed.....

   real(fpk)  :: Obsgeom_boa(MAX_GEOMETRIES,3)
   real(fpk)  :: alpha_boa(MAX_USER_STREAMS), theta_boa(MAXBEAMS), phi_boa(MAX_USER_RELAZMS)

!  Geometry routine: proxy outputs
!  ===============================

!  VSIGN = +1 (Up); -1(Down)

   real(fpk), parameter  :: vsign_up = ONE
   real(fpk), parameter  :: vsign_dn = MINUS_ONE
  
!  Alphas,  Radii, Ray constant. 

   real(fpk)  :: radii    (0:MAXLAYERS)
   real(fpk)  :: alpha    (0:MAXLAYERS,MAX_GEOMETRIES)
   real(fpk)  :: cosa     (0:MAXLAYERS,MAX_GEOMETRIES)
   real(fpk)  :: sina     (0:MAXLAYERS,MAX_GEOMETRIES)

   real(fpk)  :: radii_p    (MAX_PARTLAYERS)
   real(fpk)  :: alpha_p    (MAX_PARTLAYERS,MAX_GEOMETRIES)
   real(fpk)  :: cosa_p     (MAX_PARTLAYERS,MAX_GEOMETRIES)
   real(fpk)  :: sina_p     (MAX_PARTLAYERS,MAX_GEOMETRIES)

!  Critical layer. Not yet active 9/17/16.
!   integer    :: Ncrit(maxgeoms)
!   real(fpk)  :: RadCrit(maxgeoms), CotCrit(maxgeoms)

!  LOS Geometries, whole and partial

   real(fpk)  :: radiifine (MAXLAYERS,MAXFINELAYERS,MAX_GEOMETRIES)
   real(fpk)  :: alphafine (MAXLAYERS,MAXFINELAYERS,MAX_GEOMETRIES)
   real(fpk)  :: sinfine   (MAXLAYERS,MAXFINELAYERS,MAX_GEOMETRIES)
   real(fpk)  :: cosfine   (MAXLAYERS,MAXFINELAYERS,MAX_GEOMETRIES)

   real(fpk)  :: radiifine_p (MAX_PARTLAYERS,MAXFINELAYERS,MAX_GEOMETRIES)
   real(fpk)  :: alphafine_p (MAX_PARTLAYERS,MAXFINELAYERS,MAX_GEOMETRIES)
   real(fpk)  :: sinfine_p   (MAX_PARTLAYERS,MAXFINELAYERS,MAX_GEOMETRIES)
   real(fpk)  :: cosfine_p   (MAX_PARTLAYERS,MAXFINELAYERS,MAX_GEOMETRIES)

!  Help variables for the generalized spherical functions

   REAL(fpk) :: GSHELP(7,MAXMOMENTS_INPUT)

!  LOS VARIABLES (THERMAL SOLUTION)
!  --------------------------------

!  Geometry.

   real(fpk)  :: alpha_LOS    (0:MAXLAYERS,MAX_USER_STREAMS)
   real(fpk)  :: cosa_LOS     (0:MAXLAYERS,MAX_USER_STREAMS)
   real(fpk)  :: sina_LOS     (0:MAXLAYERS,MAX_USER_STREAMS)

   real(fpk)  :: alpha_p_LOS    (MAX_PARTLAYERS,MAX_USER_STREAMS)
   real(fpk)  :: cosa_p_LOS     (MAX_PARTLAYERS,MAX_USER_STREAMS)
   real(fpk)  :: sina_p_LOS     (MAX_PARTLAYERS,MAX_USER_STREAMS)

!  LOS Quadratures for Enhanced PS. Partials added 8/25/16.
!  5/22/20. Version 2.8.2 Upgrades. Remove radiifine_LOS and radiifine_p_LOS
!   ==> These are replaced by hfine/hfine_p variables in FOGeometry.

   real(fpk)  :: cosfine_LOS   (MAXLAYERS,MAXFINELAYERS,MAX_USER_STREAMS)
   real(fpk)  :: sinfine_LOS   (MAXLAYERS,MAXFINELAYERS,MAX_USER_STREAMS)
   real(fpk)  :: alphafine_LOS (MAXLAYERS,MAXFINELAYERS,MAX_USER_STREAMS)
!   real(fpk)  :: radiifine_LOS (MAXLAYERS,MAXFINELAYERS,MAX_USER_STREAMS)

   real(fpk)  :: cosfine_p_LOS   (MAX_PARTLAYERS,MAXFINELAYERS,MAX_USER_STREAMS)
   real(fpk)  :: sinfine_p_LOS   (MAX_PARTLAYERS,MAXFINELAYERS,MAX_USER_STREAMS)
   real(fpk)  :: alphafine_p_LOS (MAX_PARTLAYERS,MAXFINELAYERS,MAX_USER_STREAMS)
!   real(fpk)  :: radiifine_p_LOS (MAX_PARTLAYERS,MAXFINELAYERS,MAX_USER_STREAMS)

!  No criticality yet. 9/17/16
!   integer    :: Ncrit_LOS(MAX_USER_STREAMS)
!   real(fpk)  :: RadCrit_LOS(MAX_USER_STREAMS), CotCrit_LOS(MAX_USER_STREAMS)

!  LOCAL HELP VARIABLES
!  ====================

   integer   :: na_offset(MAXBEAMS,MAX_USER_STREAMS)
   integer   :: nd_offset(MAXBEAMS)
   logical   :: STARTER, do_Thermset, fail, do_Chapman, Do_spherfunc, DO_Sunlight

!  Initialize exception handling

      Gfail   = .false.
      message = ' '
      trace_1 = ' '
      trace_2 = ' '

!  Overall Flags to be set for each calculation (safety)

      do_Chapman  = VLIDORT_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION .and. VLIDORT_FixIn%Bool%TS_DO_FULLRAD_MODE
      do_Thermset = VLIDORT_FixIn%Bool%TS_DO_THERMAL_EMISSION
      starter     = VLIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES

!mick fix 3/25/2015 - added initialization

      VLIDORT_FOGeom%doNadir     = .false.

!  Phase matrices will be constructed from coefficients!
   
      Do_spherfunc = .true.
      DO_Sunlight  = .true.

!  Local Copy
!  ==========

!  Max dimensions.

      maxpartials       = MAX_PARTLAYERS
      maxgeoms          = MAX_GEOMETRIES
      maxfine           = MAXFINELAYERS
      maxszas           = MAX_SZANGLES
      maxvzas           = MAX_USER_VZANGLES
      maxazms           = MAX_USER_RELAZMS

!  directions

      do_upwelling = VLIDORT_FixIn%Bool%TS_DO_UPWELLING
      do_dnwelling = VLIDORT_FixIn%Bool%TS_DO_DNWELLING

!  geometry Modes

      do_ObsGeom  = VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY
      do_Doublet  = VLIDORT_ModIn%MBool%TS_DO_DOUBLET_GEOMETRY

!  Numbers

      nstokes         = VLIDORT_FixIn%Cont%TS_NSTOKES
      nlayers         = VLIDORT_FixIn%Cont%TS_NLAYERS
      nfine           = VLIDORT_FixIn%Cont%TS_NFINELAYERS
      nmoments_input  = VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT

!  Partials
!mick fix 3/2/2020 - trim dimensions on "partial_layeridx"

      npartials        = VLIDORT_Bookkeep%N_PARTLAYERS
      partial_layeridx(1:npartials) = VLIDORT_Bookkeep%PARTLAYERS_LAYERIDX(1:npartials)
      do_partials = ( npartials .ne. 0 )

!  Geometry numbers

      nszas   = VLIDORT_ModIn%MSunrays%TS_N_SZANGLES
      nvzas   = VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES
      nazms   = VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS
      ngeoms  = nszas ; if ( .not. do_ObsGeom ) ngeoms = ngeoms * nvzas * nazms

!  Geometry Offsets
!mick fix 3/2/2020 - trim dimensions on "na_offset"

      na_offset = 0 ; nd_offset = 0
      if ( do_doublet ) then
        nd_offset(1:nszas)         = VLIDORT_Bookkeep%SZD_OFFSETS(1:nszas)
      else if ( .not. do_obsgeom ) then
        na_offset(1:nszas,1:nvzas) = VLIDORT_Bookkeep%VZA_OFFSETS(1:nszas,1:nvzas)
      endif

!  not yet enabled
!      DoCrit = .TRUE.     
!      Acrit  = 1.0d-12

!  Configuration inputs. Removed "regular_ps" flag 9/17/16 (Redundant)

      if (VLIDORT_FixIn%Bool%TS_DO_PLANE_PARALLEL) then
        do_planpar     = .TRUE.
        do_enhanced_ps = .FALSE.
      else
        do_planpar = .FALSE.
        if (VLIDORT_ModIn%MBool%TS_DO_FOCORR_NADIR) then
          do_enhanced_ps = .FALSE.
        else if (VLIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING) then
          do_enhanced_ps = .TRUE.
        end if
      end if

!  Earth Radius and heights. Partials added, 9/17/16
!mick fix 3/2/2020 - trim dimensions on "partial_heights"

      eradius               = VLIDORT_ModIn%MCont%TS_EARTH_RADIUS
      heights(0:maxlayers)  = VLIDORT_FixIn%Chapman%TS_HEIGHT_GRID(0:maxlayers)
      partial_heights(1:npartials) = VLIDORT_Bookkeep%PARTLAYERS_HEIGHTS(1:npartials)

!  Geometry inputs. all angles in Degrees.
!  These are good for Lattice and ObsGeom

      obsgeom_boa = zero ; theta_boa = zero ; ALPHA_boa = zero ; phi_boa = zero
      if ( do_ObsGeom) then
         obsgeom_boa(1:ngeoms,1) = VLIDORT_ModIn%MSunrays%TS_SZANGLES            (1:ngeoms)
         obsgeom_boa(1:ngeoms,2) = VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT (1:ngeoms)
         obsgeom_boa(1:ngeoms,3) = VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS        (1:ngeoms)
      endif
      theta_boa(1:nszas) = VLIDORT_ModIn%MSunrays%TS_SZANGLES           (1:nszas)
      alpha_boa(1:nvzas) = VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(1:nvzas)
      phi_boa  (1:nazms) = VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS       (1:nazms)

!  Solar sources run (NO THERMAL)
!  ------------------------------

      if (  VLIDORT_ModIn%MBool%TS_do_solar_sources ) then

!  Upwelling

         if ( do_upwelling ) then

!  Geometry call. Updated 9/17/16.

            call FO_SSWPGeometry_Master &
              ( maxgeoms, maxszas, maxvzas, maxazms, maxlayers, maxpartials, maxfine, deg_to_rad, pie, & ! Input dims/constants
                vsign_up, eradius, do_obsgeom, do_doublet, do_Chapman, do_planpar, do_enhanced_ps,     & ! Input flags
                do_Partials, ngeoms, nszas, nvzas, nazms, nlayers, npartials, nfine, partial_layeridx, & ! Input control
                heights, partial_heights, obsgeom_boa, alpha_boa, theta_boa, phi_boa,                  & ! Input heights/geometry
                VLIDORT_FOGeom%doNadir, VLIDORT_FOGeom%Raycon_up, VLIDORT_FOGeom%Mu0_up,               & ! Output geometry
                VLIDORT_FOGeom%Mu1_up, VLIDORT_FOGeom%cosscat_up,                                      & ! Output geometry
                Radii,   VLIDORT_FOGeom%LosW_paths, alpha,   sina,   cosa,                             & ! Output (Layer bounds)
                VLIDORT_FOGeom%sunpaths_up, VLIDORT_FOGeom%ntraverse_up, VLIDORT_FOGeom%chapfacs,      & ! Output (Layer bounds)
                Radii_p, VLIDORT_FOGeom%LosP_paths, alpha_p, sina_p, cosa_p,                           & ! Output (partial levs)
                VLIDORT_FOGeom%sunpaths_p_up, VLIDORT_FOGeom%ntraverse_p_up, VLIDORT_FOGeom%chapfacs_p,& ! Output (partial levs)
                VLIDORT_FOGeom%nfinedivs, VLIDORT_FOGeom%xfine, VLIDORT_FOGeom%wfine,                  & ! Output Wholelayer
                radiifine,   alphafine,   sinfine,  cosfine,                                           & ! Output Wholelayer
                VLIDORT_FOGeom%sunpathsfine_up,   VLIDORT_FOGeom%ntraversefine_up,                     & ! Output Wholelayer      
                VLIDORT_FOGeom%nfinedivs_p_up, VLIDORT_FOGeom%xfine_p_up, VLIDORT_FOGeom%wfine_p_up,   & ! Output partial up
                radiifine_p, alphafine_p, sinfine_p, cosfine_p,                                        & ! Output partial up
                VLIDORT_FOGeom%sunpathsfine_p_up, VLIDORT_FOGeom%ntraversefine_p_up,                   & ! Output partial up
                fail, message, trace_1 )                                                                 ! Output (Status)

            if ( fail ) then
               trace_2 = 'Failure from FO_SSWPGeometry_Master, Solar Sources, Upwelling calculation'
               Gfail = .true. ; return
            endif

!  Spherical functions call. Updated 9/17/16.

            call FO_VectorSS_spherfuncs &
              ( MAXMOMENTS_INPUT, MAXGEOMS, MAXSZAS, MAXVZAS, MAXAZMS, deg_to_rad,       & ! Inputs
                NMOMENTS_INPUT, NGEOMS, NSZAS, NVZAS, NAZMS, NSTOKES, STARTER,           & ! Inputs
                DO_OBSGEOM, DO_DOUBLET, DO_SPHERFUNC, DO_SUNLIGHT, NA_OFFSET, ND_OFFSET, & ! Inputs
                THETA_BOA, ALPHA_BOA, PHI_BOA, VLIDORT_FOGeom%cosscat_up, vsign_up,      & ! Inputs
                VLIDORT_FOGeom%ROTATIONS_up, GSHELP, VLIDORT_FOGeom%GENSPHER_up )          ! Outputs

!  End upwelling geometry

         endif

!  Donwelling

         if ( do_dnwelling ) then

!  Geometry call. Updated 9/17/16.

            call FO_SSWPGeometry_Master &
              ( maxgeoms, maxszas, maxvzas, maxazms, maxlayers, maxpartials, maxfine, deg_to_rad, pie, & ! Input dims/constants
                vsign_dn, eradius, do_obsgeom, do_doublet, do_Chapman, do_planpar, do_enhanced_ps,     & ! Input flags
                do_Partials, ngeoms, nszas, nvzas, nazms, nlayers, npartials, nfine, partial_layeridx, & ! Input control
                heights, partial_heights, obsgeom_boa, alpha_boa, theta_boa, phi_boa,                  & ! Input heights/geometry
                VLIDORT_FOGeom%doNadir, VLIDORT_FOGeom%Raycon_dn, VLIDORT_FOGeom%Mu0_dn,               & ! Output geometry
                VLIDORT_FOGeom%Mu1_dn, VLIDORT_FOGeom%cosscat_dn,                                      & ! Output geometry
                Radii,   VLIDORT_FOGeom%LosW_paths, alpha,   sina,   cosa,                             & ! Output (Layer bounds)
                VLIDORT_FOGeom%sunpaths_dn, VLIDORT_FOGeom%ntraverse_dn, VLIDORT_FOGeom%chapfacs,      & ! Output (Layer bounds)
                Radii_p, VLIDORT_FOGeom%LosP_paths, alpha_p, sina_p, cosa_p,                           & ! Output (partial levs)
                VLIDORT_FOGeom%sunpaths_p_dn, VLIDORT_FOGeom%ntraverse_p_dn, VLIDORT_FOGeom%chapfacs_p,& ! Output (partial levs)
                VLIDORT_FOGeom%nfinedivs, VLIDORT_FOGeom%xfine, VLIDORT_FOGeom%wfine,                  & ! Output Wholelayer
                radiifine,   alphafine,   sinfine,  cosfine,                                           & ! Output Wholelayer
                VLIDORT_FOGeom%sunpathsfine_dn,   VLIDORT_FOGeom%ntraversefine_dn,                     & ! Output Wholelayer      
                VLIDORT_FOGeom%nfinedivs_p_dn, VLIDORT_FOGeom%xfine_p_dn, VLIDORT_FOGeom%wfine_p_dn,   & ! Output partial dn
                radiifine_p, alphafine_p, sinfine_p, cosfine_p,                                        & ! Output partial dn
                VLIDORT_FOGeom%sunpathsfine_p_dn, VLIDORT_FOGeom%ntraversefine_p_dn,                   & ! Output partial dn
                fail, message, trace_1 )                                                                 ! Output (Status)

            if ( fail ) then
               trace_2 = 'Failure from FO_SSWPGeometry_Master, Solar Sources, Downwelling calculation'
               Gfail = .true. ; return
            endif

!  Spherical functions call. Updated 9/17/16.

            call FO_VectorSS_spherfuncs &
              ( MAXMOMENTS_INPUT, MAXGEOMS, MAXSZAS, MAXVZAS, MAXAZMS, deg_to_rad,       & ! Inputs
                NMOMENTS_INPUT, NGEOMS, NSZAS, NVZAS, NAZMS, NSTOKES, STARTER,           & ! Inputs
                DO_OBSGEOM, DO_DOUBLET, DO_SPHERFUNC, DO_SUNLIGHT, NA_OFFSET, ND_OFFSET, & ! Inputs
                THETA_BOA, ALPHA_BOA, PHI_BOA, VLIDORT_FOGeom%cosscat_dn, vsign_dn,      & ! Inputs
                VLIDORT_FOGeom%ROTATIONS_dn, GSHELP, VLIDORT_FOGeom%GENSPHER_dn )          ! Outputs

!  End downwelling

         endif

!  End solar only

      endif

!  Thermal sources, geometry
!  -------------------------

      if ( do_thermset ) then

!  Upwelling DT Geometry call. Updated 9/17/16.
!  5/22/20. Version 2.8.2. Upgrades.
!     ==> Must use vertical distances in Thermal source terms (not path distances, bug corrected)
!     ==> Replace radiifine/radiifine_p output with VLIDORT_FOGeom%hfine_LOS and VLIDORT_FOGeom%hfine_p_LOS
!     ==> Distinguish between the upwelling and downwelling values for this inpu

         if ( do_upwelling ) then
           call FO_DTWPGeometry_Master  &
            (maxvzas, maxlayers, maxpartials, maxfine, deg_to_rad, eradius,                                & ! Input dims/constants
             .true., do_planpar, do_enhanced_ps, do_Partials,                                              & ! Input flags
             nvzas, nlayers, npartials, nfine, partial_layeridx,                                           & ! Input control
             heights, alpha_boa, partial_heights, VLIDORT_FOGeom%Mu1_LOS,                                  & ! Input heights/geom
             Radii,   VLIDORT_FOGeom%LosW_paths_LOS, alpha_LOS,   sina_LOS,   cosa_LOS,                    & ! Output (Layer bounds)
             Radii_p, VLIDORT_FOGeom%LosP_paths_LOS, alpha_p_LOS, sina_p_LOS, cosa_p_LOS,                  & ! Output (partial levs)
             VLIDORT_FOGeom%nfinedivs_LOS, VLIDORT_FOGeom%xfine_LOS, VLIDORT_FOGeom%wfine_LOS,             & ! Output Wholelayer
             VLIDORT_FOGeom%hfine_LOS_up, alphafine_LOS, sinfine_LOS, cosfine_LOS,                         & ! Output Wholelayer
             VLIDORT_FOGeom%nfinedivs_p_LOS_up,VLIDORT_FOGeom%xfine_p_LOS_up,VLIDORT_FOGeom%wfine_p_LOS_up,& ! Output partial
             VLIDORT_FOGeom%hfine_p_LOS_up, alphafine_p_LOS, sinfine_p_LOS, cosfine_p_LOS,                 & ! Output partial
             fail, message, trace_1)                                                                         ! Output (Status)
           if ( fail ) then
              trace_2 = 'Failure from FO_DTWPGeometry_Master, upwelling call'
              Gfail = .true. ; return
           endif
         endif

!  Downwelling DT Geometry call. Updated 9/17/16.

         if ( do_dnwelling ) then
           call FO_DTWPGeometry_Master  &
            (maxvzas, maxlayers, maxpartials, maxfine, deg_to_rad, eradius,                                & ! Input dims/constants
             .false., do_planpar, do_enhanced_ps, do_Partials,                                             & ! Input flags
             nvzas, nlayers, npartials, nfine, partial_layeridx,                                           & ! Input control
             heights, alpha_boa, partial_heights, VLIDORT_FOGeom%Mu1_LOS,                                  & ! Input heights/geom
             Radii,   VLIDORT_FOGeom%LosW_paths_LOS, alpha_LOS,   sina_LOS,   cosa_LOS,                    & ! Output (Layer bounds)
             Radii_p, VLIDORT_FOGeom%LosP_paths_LOS, alpha_p_LOS, sina_p_LOS, cosa_p_LOS,                  & ! Output (partial levs)
             VLIDORT_FOGeom%nfinedivs_LOS, VLIDORT_FOGeom%xfine_LOS, VLIDORT_FOGeom%wfine_LOS,             & ! Output Wholelayer
             VLIDORT_FOGeom%Hfine_LOS_dn, alphafine_LOS, sinfine_LOS, cosfine_LOS,                         & ! Output Wholelayer
             VLIDORT_FOGeom%nfinedivs_p_LOS_dn,VLIDORT_FOGeom%xfine_p_LOS_dn,VLIDORT_FOGeom%wfine_p_LOS_dn,& ! Output partial
             VLIDORT_FOGeom%hfine_p_LOS_dn, alphafine_p_LOS, sinfine_p_LOS, cosfine_p_LOS,                 & ! Output partial
             fail, message, trace_1)                                                                         ! Output (Status)
           if ( fail ) then
              trace_2 = 'Failure from FO_DTWPGeometry_Master, downwelling call'
              Gfail = .true. ; return
           endif
         endif

!  End thermal run

      endif

!  Finish

      return
end subroutine VFO_GEOMETRY_MASTER

!  End Module

end module VFO_Geometry_Master_m


