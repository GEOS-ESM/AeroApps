
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

!  FO Version history
!  ------------------

!  Versions to 1.4, without Partials. Code is stand alone with no dependencies.
!  Version     1.5, with optional phase matrices, Surface Leaving and partials.

!    Version 1a, 01 December 2011, R. Spurr, RT Solutions Inc.
!    Version 1b, 13 February 2012, R. Spurr, RT Solutions Inc.
!    Version 2,  01 June     2012, R. Spurr, RT Solutions Inc.
!    Version 3,  29 October  2012, Extension to Observational multiple geometries
!    Version 4,  31 July     2013, Lattice Multi-geometry
!    Version 5,  07 July     2016, Optional F-matrix usage
!    Version 5,  02 August   2016. Sleave + Jacobians
!    Version 5,  19 August   2016, Partial-layer output
!    Version 5,  11 December 2017, Optional code for polarized emissivity

!   4/09/19. Version 1.5.1. Add CUMTRANS output, and Waterleaving input control
!   4/15/20. Version 1.5.2. Add doublet geometry option

!  VLIDORT Interface history
!  -------------------------

!    FO Version 1.4: This module is interface with VLIDORT V2.7. R.Spurr 3/19/15
!    FO Version 1.5: Interface module upgraded to  VLIDORT V2.8. R.Spurr 7/7/16, 9/17/16

!  5/22/20. Version 2.8.2 Upgrades. Direct-Thermal upgrade. 
!    ==> Geometrical calculation for direct thermal outgoing downwelling and upwelling was corrected.
!    ==> Water-leaving cumulative tranmsittance was properly initialized.

module VFO_RTCalc_Master_m

!  Use statements

   Use VLIDORT_pars_m
   USE VLIDORT_Setups_def_m

public

contains

subroutine VFO_RTCALC_MASTER &
       ( do_solar_sources, do_thermal_emission, do_surface_emission, do_Polarized_Emissivity, & ! Input flags (sources)
         do_upwelling, do_dnwelling, do_lambertian, do_surface_leaving, do_water_leaving,     & ! Input flags (Up/dn/surface)
         do_sunlight, do_deltam_scaling, do_obsgeom, do_doublet, do_Partials, do_planpar,     & ! Input flags (RT Control)
         do_enhanced_ps, nstokes, ngeoms, nszas, nvzas, nazms, na_offset, nd_offset, nlayers, & ! Inputs (control-numbers)
         n_user_levels, user_level_mask_up, user_level_mask_dn,                               & ! Inputs (control-levels)
         npartials, partial_outindex, partial_outflag, partial_layeridx, FOGeometry,          & ! Inputs (control-partial)
         flux, fluxvec, extinction, deltaus, omega, truncfac, fmatrix_up, fmatrix_dn,         & ! Input atmos optical
         bb_input, surfbb, emiss, reflec, slterm,                                             & ! Input thermal/surface/Geometry
         fo_stokes_ss, fo_stokes_db, fo_stokes_dta, fo_stokes_dts,                            & ! Output
         fo_stokes_atmos, fo_stokes_surf, fo_stokes, cumtrans )                                 ! Output

!  4/9/19. Add CUMTRANS output, and Waterleaving input control
!  4/15/20. Version 1.5.2. Add doublet flag and offsets

!  Use modules

   USE FO_Vector_SSRT_m
   USE FO_Thermal_DTRT_m

   implicit none

!  parameter arguments

   integer, parameter :: ffp = selected_real_kind(15)

!  Subroutine inputs
!  =================

!  Input/setup Proxies (input to this Master routine)
!  --------------------------------------------------

!  Sources control, including thermal

   logical, intent(in)  :: do_solar_sources, do_thermal_emission, do_surface_emission

!  Polarized Emissivity flag, 12/11/17 Rob added
!      4/15/20. Version 1.5.2. Introduced as an input here

   LOGICAL, Intent(in)  :: do_Polarized_Emissivity

!  Directional Flags

   logical, intent(in)  :: do_upwelling, do_dnwelling

!  Sleave flag added, 8/2/16, water-leaving 4/9/19

   logical, intent(in)  :: do_Lambertian, do_surface_leaving, do_water_leaving

!  Sunlight and deltam scaling flags

   logical, intent(in)  :: do_sunlight, do_deltam_scaling

!  Geometry Flags (sphericity flags are mutually exclusive). Regular PS now removed, version 1.5
!   ==> 4/15/20. Version 1.5.2. Add doublet flag

   logical, intent(in)  :: do_Obsgeom, do_doublet, do_planpar, do_enhanced_ps

!  Layer and geometry control. 
!    ==> 4/15/20. Version 1.5.2. Add doublet offsets

   integer, intent(in)  :: nstokes, ngeoms, nszas, nvzas, nazms, nlayers
   integer, intent(in)  :: na_offset(MAX_SZANGLES,MAX_USER_VZANGLES)
   integer, intent(in)  :: nd_offset(MAX_SZANGLES)

!  Output levels. Use masking as in main codes, 9/17/16

   integer, intent(in)  :: n_user_levels
   integer, intent(in)  :: user_level_mask_up ( max_user_levels )
   integer, intent(in)  :: user_level_mask_dn ( max_user_levels )

!  Control for partial-layer output, added 8/25/16

   logical, intent(in)  :: do_Partials
   integer, intent(in)  :: Npartials
   integer, intent(in)  :: partial_layeridx( max_partlayers  )
   logical, intent(in)  :: partial_outflag ( MAX_USER_LEVELS )
   integer, intent(in)  :: partial_outindex( MAX_USER_LEVELS )

!  Optical Inputs
!  --------------

!  Solar flux

   real(ffp), intent(in) :: flux, fluxvec(maxstokes)

!  Atmosphere
!   4/15/20. Version 1.5.2. Greekmat removed; Now we must use FMATRICES

   real(ffp), intent(in) :: extinction  ( maxlayers )
   real(ffp), intent(in) :: deltaus     ( maxlayers )
   real(ffp), intent(in) :: omega       ( maxlayers )
   real(ffp), intent(in) :: truncfac    ( maxlayers )

!  Fmatrix input. (FO 1.5, VLIDORT 2.8). Introduced 7/7/16

   real(ffp), intent(in) :: fmatrix_up  ( maxlayers, Max_geometries, 6 )
   real(ffp), intent(in) :: fmatrix_dn  ( maxlayers, Max_geometries, 6 )

!  Thermal inputs, surface emissivity.
!   Emissivity is now polarized. 12/11/17 Rob add.

   real(ffp), intent(in) :: bb_input ( 0:maxlayers )
   real(ffp), intent(in) :: surfbb
   real(ffp), intent(in) :: emiss ( maxstokes, max_user_vzangles )

!  Surface properties - reflective (could be the albedo), surface leaving added 8/2/16

   real(ffp), intent(in) :: reflec(maxstokes,maxstokes,Max_geometries)
   real(ffp), intent(in) :: slterm(maxstokes,Max_geometries)

!  Type structures, FO Geometry output
!  -----------------------------------

   TYPE(VLIDORT_Geometry_FO), intent(in)  :: FOGeometry

!  Subroutine outputs
!  ==================

!  Solar

   real(ffp), intent(out) :: fo_stokes_ss ( max_user_levels,MAX_GEOMETRIES,maxstokes,max_directions )
   real(ffp), intent(out) :: fo_stokes_db ( max_user_levels,MAX_GEOMETRIES,maxstokes )

!  Thermal

   real(ffp), intent(out) :: fo_stokes_dta ( max_user_levels,MAX_GEOMETRIES,maxstokes,max_directions )
   real(ffp), intent(out) :: fo_stokes_dts ( max_user_levels,MAX_GEOMETRIES,maxstokes )

!  Composite
!mick mod 9/19/2017 - added "fo_stokes_atmos" & "fo_stokes_surf" to vector FO output to
!                     keep FO scalar & vector subroutines on an equal footing
                     
   real(ffp), intent(out) :: fo_stokes_atmos ( max_user_levels,MAX_GEOMETRIES,maxstokes,max_directions )
   real(ffp), intent(out) :: fo_stokes_surf  ( max_user_levels,MAX_GEOMETRIES,maxstokes )
   real(ffp), intent(out) :: fo_stokes       ( max_user_levels,MAX_GEOMETRIES,maxstokes,max_directions )

!  4/9/19. Additional output for the sleave correction

   real(ffp), Intent(out) :: CUMTRANS ( max_user_levels, MAX_GEOMETRIES )

!  Other variables
!  ===============

!  Source flags
!  ------------

!  Existence flags (Solar). 8/19/16. Criticality enters here

   logical    :: do_sources_up       (MAXLAYERS,MAX_GEOMETRIES)
   logical    :: do_sources_dn       (MAXLAYERS,MAX_GEOMETRIES)
   logical    :: do_sources_up_p     (MAX_PARTLAYERS,MAX_GEOMETRIES)
   logical    :: do_sources_dn_p     (MAX_PARTLAYERS,MAX_GEOMETRIES)

!  Existence flags (Thermal). 8/19/16. Criticality enters here

   logical    :: do_Tsources_up       (MAXLAYERS,MAX_USER_STREAMS)
   logical    :: do_Tsources_dn       (MAXLAYERS,MAX_USER_STREAMS)
   logical    :: do_Tsources_up_p     (MAX_PARTLAYERS,MAX_USER_STREAMS)
   logical    :: do_Tsources_dn_p     (MAX_PARTLAYERS,MAX_USER_STREAMS)

!  RT calculation outputs
!  ----------------------

!  Solar routines

   real(ffp)  :: stokes_up    ( max_user_levels,maxstokes,MAX_GEOMETRIES )
   real(ffp)  :: stokes_dn    ( max_user_levels,maxstokes,MAX_GEOMETRIES )
   real(ffp)  :: stokes_db    ( max_user_levels,maxstokes,MAX_GEOMETRIES )

!  Thermal routines (Scalar, no polarization). SEE LOS VARIABLES (THERMAL), below
!   real(ffp)  :: intensity_dta_up ( max_user_levels,MAX_GEOMETRIES )
!   real(ffp)  :: intensity_dta_dn ( max_user_levels,MAX_GEOMETRIES )
!   real(ffp)  :: intensity_dts    ( max_user_levels,MAX_GEOMETRIES )

!  Intermediate RT products
!  ------------------------

!  Composite
!mick mod 9/19/2017 - added "fo_stokes_atmos" & "fo_stokes_surf" to vector FO output to
!                     keep FO scalar & vector subroutines on an equal footing

   !real(ffp)  :: fo_stokes_atmos ( max_user_levels,MAX_GEOMETRIES,maxstokes,max_directions )
   !real(ffp)  :: fo_stokes_surf  ( max_user_levels,MAX_GEOMETRIES,maxstokes )

!  LOS VARIABLES (THERMAL SOLUTION)
!  --------------------------------

!  Output

   real(ffp)  :: intensity_dta_up_LOS ( max_user_levels, MAX_USER_VZANGLES )
   real(ffp)  :: intensity_dta_dn_LOS ( max_user_levels, MAX_USER_VZANGLES )
   real(ffp)  :: intensity_dts_LOS    ( max_user_levels, MAX_USER_VZANGLES )
   real(ffp)  :: cumsource_up_LOS     ( 0:maxlayers, MAX_USER_VZANGLES )
   real(ffp)  :: cumsource_dn_LOS     ( 0:maxlayers, MAX_USER_VZANGLES )

!  5/22/20. Version 2.8.2 Upgrades.
!    ==> LOSTRANS output might be used again for the (upwelling) thermal-NoScattering contribution??

   real(ffp)  :: lostrans_up_LOS      ( maxlayers     ,MAX_USER_VZANGLES )
   real(ffp)  :: lostrans_up_p_LOS    ( MAX_PARTLAYERS,MAX_USER_VZANGLES )

!  StokesQUV_dts_LOS, contribution from Polarized emissivity. 12/11/17 Rob add.

   real(ffp)  :: StokesQUV_dts_LOS    ( max_user_levels, 3, MAX_USER_VZANGLES )

!  Other products
!  --------------

!  Thermal setup

   real(ffp)  :: tcom1(maxlayers,2)

!  Dummies

   real(ffp)  :: SScumsource_up ( 0:maxlayers,maxstokes,MAX_GEOMETRIES )
   real(ffp)  :: SScumsource_dn ( 0:maxlayers,maxstokes,MAX_GEOMETRIES )

!   real(ffp)  :: DTcumsource_up ( 0:maxlayers,MAX_GEOMETRIES )
!   real(ffp)  :: DTcumsource_dn ( 0:maxlayers,MAX_GEOMETRIES )

!  LOCAL HELP VARIABLES
!  --------------------

   integer   :: ns, nv, na, g, lev, o1

!  Initialize output
!mick mod 9/19/2017 - initialized "fo_stokes_atmos" & "fo_stokes_surf"

   fo_stokes_ss    = zero
   fo_stokes_db    = zero
   fo_stokes_dta   = zero
   fo_stokes_dts   = zero
   fo_stokes_atmos = zero
   fo_stokes_surf  = zero
   fo_stokes       = zero

!  Set do_Polarized_Emissivity flag. 12/11/17 Rob add.
!    Now set outside this module , using the following (commented out) code
!   do_Polarized_Emissivity = .false.
!   if ( nstokes.gt.1 .and. do_surface_emission ) then
!      do nv = 1, nvzas
!        if ( SUM(emiss(2:4,nv)).ne.zero ) do_Polarized_Emissivity = .true.
!      enddo
!   endif

!  Solar sources run (NO THERMAL)
!  ------------------------------

   if ( do_solar_sources ) then

!  Upwelling

     if ( do_upwelling ) then

!  Sources, temporary fix until criticality realized. 9/17/16

       do_sources_up   = .true. ; do_sources_up_p = .true.

!  RT Call Solar only
!  - Updated to include surface leaving, 8/2/16. Updated 9/17/16.
!  4/9/19. water-leaving control, Additional cumtrans output for the sleave correction

      Call FO_Vector_SSRT_UP &
        ( do_sunlight, do_deltam_scaling, do_Lambertian, do_surface_leaving, do_water_leaving, & ! Inputs (Flags-General/Surface)
          do_Partials, do_PlanPar, do_enhanced_ps, do_sources_up, do_sources_up_p,             & ! Inputs(Flags/criticality)
          nstokes, ngeoms, nlayers, n_user_levels, user_level_mask_up,                         & ! Inputs (control/flux)
          npartials, partial_outindex, partial_outflag, partial_layeridx, FOGeometry,          & ! Inputs (control-partial)
          flux, fluxvec, extinction, deltaus, omega, truncfac, fmatrix_up, Reflec, Slterm,     & ! Inputs (Geometry/Optical/surface)
          stokes_up, stokes_db, SScumsource_up, cumtrans )                                       ! Outputs

!  Save results
!mick fix 9/19/2017 - turned off "fo_stokes" (defined later)

       do o1 = 1, nstokes
         do g = 1, ngeoms
           do lev=1,n_user_levels
             fo_stokes_ss(lev,g,o1,upidx) = stokes_up(lev,o1,g)
             fo_stokes_db(lev,g,o1)       = stokes_db(lev,o1,g)
             !fo_stokes(lev,g,o1,upidx)    = fo_stokes_ss(lev,g,o1,upidx) + fo_stokes_db(lev,g,o1)
           enddo
         enddo
       enddo

!  End upwelling

     endif

!  Donwelling

     if ( do_dnwelling ) then

!  Sources, temporary fix until criticality realized. 9/17/16

       do_sources_dn   = .true. ; do_sources_dn_p = .true.

!  RT Call Solar only. Updated 9/17/16.

       Call FO_Vector_SSRT_DN &
        ( do_sunlight, do_deltam_scaling, do_Partials,                                & ! Inputs (Flags)
          do_PlanPar, do_enhanced_ps, do_sources_dn, do_sources_dn_p,                 & ! Inputs (Flags)
          nstokes, ngeoms, nlayers, n_user_levels, user_level_mask_dn,                & ! Inputs (control/flux)
          npartials, partial_outindex, partial_outflag, partial_layeridx, FOGeometry, & ! Inputs (control-partial)
          flux, fluxvec, extinction, deltaus, omega, truncfac, fmatrix_dn,            & ! Inputs (Optical)
          stokes_dn, SScumsource_dn )                                                   ! Outputs

!  Save results
!mick fix 9/19/2017 - turned off "fo_stokes" (defined later)

       do o1 = 1, nstokes
         do g = 1, ngeoms
           do lev=1,n_user_levels
             fo_stokes_ss(lev,g,o1,dnidx) = stokes_dn(lev,o1,g)
             !fo_stokes(lev,g,o1,dnidx)    = fo_stokes_ss(lev,g,o1,dnidx)
           enddo
         enddo
       enddo

!  End downwelling

     endif

!  End solar only

   endif

!  Thermal sources run
!  -------------------

   if ( do_thermal_emission.and.do_surface_emission ) then

!  Upwelling
!  ---------

     if ( do_upwelling ) then

!  Sources, temporary fix until criticality realized. 9/17/16

       do_Tsources_up   = .true. ; do_Tsources_up_p = .true.

!  Direct thermal, calculate. Updated 9/17/16
!    -- If Polarized Emissivity flag present, then use optional call. 12/11/17 Rob add.
!    -- Array "Emiss" now has vector dimension.
!mick fix 3/2/2020 - replaced input "ngeoms" with "nvzas"

!  5/22/20. Version 2.8.2 Upgrades.
!        ==> Add lostrans_up_LOS + lostrans_up_p_LOS to argument list

       call FO_Thermal_DTRT_UP &
         ( do_deltam_scaling, do_Partials, Do_Polarized_Emissivity,        & ! Inputs (Flags)
           do_PlanPar, do_enhanced_ps, do_Tsources_up, do_Tsources_up_p,   & ! Inputs (Flags)
           nstokes, nvzas, nlayers, n_user_levels, user_level_mask_up,     & ! Inputs (control output)
           npartials, partial_outindex, partial_outflag, partial_layeridx, & ! Inputs (control-partial)
           FOGeometry, bb_input, surfbb, emiss(1,:), emiss(2:4,:),         & ! Inputs (Thermal)
           extinction, deltaus, omega, truncfac,                           & ! Inputs (Optical)
           intensity_dta_up_LOS, intensity_dts_LOS, StokesQUV_dts_LOS,     & ! Main Outputs
           cumsource_up_LOS, tcom1, lostrans_up_LOS, lostrans_up_p_LOS )     ! Other Outputs

!  Save results
!mick fix 9/19/2017 - turned off "fo_stokes" (defined later)
!  4/15/20. Version 1.5.2.  Add the Doublet goemetry option

       o1=1
       if ( do_obsgeom ) then
         do g = 1, nvzas
           do lev=1,n_user_levels
             fo_stokes_dta(lev,g,o1,upidx) = intensity_dta_up_LOS(lev,g)
             fo_stokes_dts(lev,g,o1)       = intensity_dts_LOS(lev,g)
           enddo
         enddo
       else if ( do_Doublet ) then
          do nv = 1, nvzas ; do ns = 1, nszas
             g = nd_offset(ns) + nv
             do lev=1,n_user_levels
                fo_stokes_dta(lev,g,o1,upidx) = intensity_dta_up_LOS(lev,nv)
                fo_stokes_dts(lev,g,o1)       = intensity_dts_LOS(lev,nv)
             enddo
          enddo ; enddo
       else
          do nv = 1, nvzas ; do ns = 1, nszas ; do na = 1, nazms
             g = na_offset(ns,nv) + na
             do lev=1,n_user_levels
                fo_stokes_dta(lev,g,o1,upidx) = intensity_dta_up_LOS(lev,nv)
                fo_stokes_dts(lev,g,o1)       = intensity_dts_LOS(lev,nv)
                !fo_stokes(lev,g,o1,upidx)     = fo_stokes_dta(lev,g,o1,upidx) + fo_stokes_dts(lev,g,o1)
             enddo
          enddo ; enddo ; enddo
       endif

!  Save polarized Emissivity results. 12/11/17  Rob add.
!  4/15/20. Version 1.5.2.  Add the Doublet goemetry option

       if ( do_Polarized_Emissivity.and.nstokes.gt.1 ) then
         if ( do_obsgeom ) then
           do g = 1, nvzas
             do lev=1,n_user_levels
               fo_stokes_dts(lev,g,2:nstokes) = StokesQUV_dts_LOS(lev,1:nstokes-1,g)
             enddo
           enddo
         else if ( do_doublet ) then
            do nv = 1, nvzas ; do ns = 1, nszas 
               g = nd_offset(ns) + nv
               do lev=1,n_user_levels
                  fo_stokes_dts(lev,g,2:nstokes) = StokesQUV_dts_LOS(lev,1:nstokes-1,nv)
               enddo
            enddo ; enddo
         else
            do nv = 1, nvzas ; do ns = 1, nszas ; do na = 1, nazms
               g = na_offset(ns,nv) + na
               do lev=1,n_user_levels
                  fo_stokes_dts(lev,g,2:nstokes) = StokesQUV_dts_LOS(lev,1:nstokes-1,nv)
               enddo
            enddo ; enddo ; enddo
         endif
       endif

!  End upwelling

     endif

!  Downwelling
!  -----------

     if ( do_dnwelling ) then

!  Sources, temporary fix until criticality realized. 9/17/16

       do_Tsources_dn   = .true. ; do_Tsources_dn_p = .true.

!  Direct thermal, calculate. Updated, 9/17/16.
!mick fix 3/2/2020 - replaced input "ngeoms" with "nvzas"

       call FO_Thermal_DTRT_DN &
        ( do_deltam_scaling, do_Partials, do_PlanPar,                   & ! Inputs (Flags)
          do_enhanced_ps, do_Tsources_dn, do_Tsources_dn_p,             & ! Inputs (Flags)
          nvzas, nlayers, n_user_levels, user_level_mask_dn, npartials, & ! Inputs (control output)
          partial_outindex, partial_outflag, partial_layeridx,          & ! Inputs (control-partial)
          FOGeometry, BB_input, extinction, deltaus, omega, truncfac,   & ! Inputs (Optical)
          intensity_dta_dn_LOS, cumsource_dn_LOS, tcom1 )                 ! Outputs

!  Save results
!mick fix 9/19/2017 - turned off "fo_stokes" (defined later)
!  4/15/20. Version 1.5.2.  Add the Doublet goemetry option

       o1=1
       if ( do_obsgeom ) then
          do g = 1, nvzas
             do lev=1,n_user_levels
                fo_stokes_dta(lev,g,o1,dnidx) = intensity_dta_dn_LOS(lev,g)
                !fo_stokes(lev,g,o1,dnidx)     = fo_stokes_dta(lev,g,o1,dnidx)
             enddo
          enddo
       else if ( do_doublet ) then
          do nv = 1, nvzas ; do ns = 1, nszas
             g = nd_offset(ns) + nv
             do lev=1,n_user_levels
                fo_stokes_dta(lev,g,o1,dnidx) = intensity_dta_dn_LOS(lev,nv)
                !fo_stokes(lev,g,o1,dnidx)     = fo_stokes_dta(lev,g,o1,dnidx)
             enddo
          enddo ; enddo
       else
          do nv = 1, nvzas ; do ns = 1, nszas ; do na = 1, nazms
             g = na_offset(ns,nv) + na
             do lev=1,n_user_levels
                fo_stokes_dta(lev,g,o1,dnidx) = intensity_dta_dn_LOS(lev,nv)
                !fo_stokes(lev,g,o1,dnidx)     = fo_stokes_dta(lev,g,o1,dnidx)
             enddo
          enddo ; enddo ; enddo
       endif

!  end downwelling

     endif

!  End thermal run

   endif

!  Final computation of Composites
!  -------------------------------
!mick fix 9/19/2017 - this section overhauled to avoid foreseen bugs & simplify computations

   if ( do_upwelling ) then
     do o1 = 1, nstokes
       do g = 1, ngeoms
         do lev=1,n_user_levels
           fo_stokes_atmos(lev,g,o1,upidx) = fo_stokes_ss(lev,g,o1,upidx)    + fo_stokes_dta(lev,g,o1,upidx)
           fo_stokes_surf(lev,g,o1)        = fo_stokes_db(lev,g,o1)          + fo_stokes_dts(lev,g,o1)
           fo_stokes(lev,g,o1,upidx)       = fo_stokes_atmos(lev,g,o1,upidx) + fo_stokes_surf(lev,g,o1)
         enddo
       enddo
     enddo
   endif

   if ( do_dnwelling ) then
     do o1 = 1, nstokes
       do g = 1, ngeoms
         do lev=1,n_user_levels
           fo_stokes_atmos(lev,g,o1,dnidx) = fo_stokes_ss(lev,g,o1,dnidx) + fo_stokes_dta(lev,g,o1,dnidx)
           fo_stokes(lev,g,o1,dnidx)       = fo_stokes_atmos(lev,g,o1,dnidx)
         enddo
       enddo
     enddo
   endif

!  Finish

   return
end subroutine VFO_RTCALC_MASTER

end module VFO_RTCalc_Master_m

