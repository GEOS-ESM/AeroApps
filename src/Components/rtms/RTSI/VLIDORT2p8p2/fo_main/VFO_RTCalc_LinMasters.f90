
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
!  Version     1.5, with optional F-matrices, Surface leaving, and partials.

!    Version 1a, 01 December 2011, R. Spurr, RT Solutions Inc.
!    Version 1b, 13 February 2012, R. Spurr, RT Solutions Inc.
!    Version 2,  01 June     2012, R. Spurr, RT Solutions Inc.
!    Version 3,  29 October  2012, Extension to Observational multiple geometries
!    Version 4,  31 July     2013, Lattice Multi-geometry
!    Version 5,  07 July     2016, Optional F-matrix usage
!    Version 5,  02 August   2016. Surface leaving and Sleave Jacobians
!    Version 5,  25 August   2016, Partial-layer output
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

module VFO_RTCalc_LinMasters_m

!  Use statements

   Use VLIDORT_pars_m
   USE VLIDORT_Setups_def_m

!  All subroutines public

public  :: VFO_RTCALC_LPS_MASTER, &
           VFO_RTCALC_LCS_MASTER


contains

subroutine VFO_RTCALC_LPS_MASTER &
       ( do_solar_sources, do_thermal_emission, do_surface_emission, do_Polarized_Emissivity, & ! Input flags (sources)
         do_upwelling, do_dnwelling, do_lambertian, do_surface_leaving, do_water_leaving,     & ! Input flags (dir/surface)
         do_sunlight, do_deltam_scaling, do_obsgeom, do_doublet, do_Partials, do_planpar,     & ! Input flags (surface/geoms)
         do_enhanced_ps, do_profilewfs, do_surfacewfs, do_sleavewfs,                          & ! Input Lin flags
         nstokes, ngeoms, nszas, nvzas, nazms, na_offset, nd_offset, nlayers,                 & ! Input Numbers/offsets
         n_user_levels, user_level_mask_up, user_level_mask_dn,                           & ! Input control-levels
         npartials, partial_outindex, partial_outflag, partial_layeridx, FOGeometry,      & ! Input control-partial/Geometry
         Lvaryflags, Lvarynums, LvaryFmat, n_reflecwfs, n_sleavewfs, n_surfacewfs,        & ! Input Lin control
         flux, fluxvec, extinction, deltaus, omega, truncfac, fmatrix_up, fmatrix_dn,     & ! Input atmos optical
         bb_input, surfbb, emiss, reflec, slterm,                                         & ! Input thermal/surface
         L_extinction, L_deltaus, L_omega, L_truncfac, L_fmatrix_up, L_fmatrix_dn,        & ! Input Lin atmos optical
         LS_emiss, LS_reflec, LSSL_slterm,                                                & ! Input Lin thermal/surf
         fo_stokes_ss, fo_stokes_db, fo_stokes_dta, fo_stokes_dts,                        & ! Output - Stokes   Vectors
         fo_profilewf_ss,  fo_profilewf_db,  fo_profilewf_dta, fo_profilewf_dts,          & ! Output - Profile  Jacobians
         fo_surfacewf_db,  fo_surfacewf_dts,                                              & ! Output - Surface  Jacobians
         fo_stokes_atmos,  fo_stokes_surf, fo_stokes,                                     & ! Output - Stokes   composites
         fo_profilewf_atmos, fo_profilewf_surf, fo_profilewf, fo_surfacewf,               & ! Output - Jacobian composites
         cumtrans, LP_cumtrans )                                                            ! Output - ancillary

!  4/9/19. Add CUMTRANS and LP_cumtrans output, and Waterleaving input control
!  4/15/20. Version 1.5.2. Add doublet flag and offsets

!  Use modules

   USE FO_Vector_SSRT_ILPS_m
   USE FO_Thermal_DTRT_ILPS_m

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
!      4/15/20/20. Version 2.8.2. Introduced as an input here

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
!   ==> 4/15/20. Version 1.5.2. Add doublet offsets

   integer, intent(in)  :: nstokes, ngeoms, nszas, nvzas, nazms, nlayers
   integer, intent(in)  :: na_offset(MAX_SZANGLES,MAX_USER_VZANGLES)
   integer, intent(in)  :: nd_offset(MAX_SZANGLES)

!  Output levels. Use masking as in main codes, 9/17/16

   integer, intent(in)  :: n_user_levels
   integer, intent(in)  :: user_level_mask_up ( MAX_USER_LEVELS )
   integer, intent(in)  :: user_level_mask_dn ( MAX_USER_LEVELS )

!  Control for partial-layer output, added 8/25/16

   logical, intent(in)  :: do_Partials
   integer, intent(in)  :: Npartials
   integer, intent(in)  :: partial_layeridx( max_partlayers  )
   logical, intent(in)  :: partial_outflag ( MAX_USER_LEVELS )
   integer, intent(in)  :: partial_outindex( MAX_USER_LEVELS )

!  Jacobian Flags. do_sleavewfs added 8/2/16

   LOGICAL, Intent(in) :: do_profilewfs
   LOGICAL, Intent(in) :: do_surfacewfs
   LOGICAL, Intent(in) :: do_sleavewfs

!  Jacobian control. Reflec and sleave numbers added, 8/2/16
!    Note that n_surfacewfs = n_reflecwfs + n_sleavewfs

   INTEGER, Intent(in) :: n_reflecwfs
   INTEGER, Intent(in) :: n_sleavewfs
   INTEGER, Intent(in) :: n_surfacewfs
   LOGICAL, Intent(in) :: Lvaryflags(maxlayers)
   INTEGER, Intent(in) :: Lvarynums (maxlayers)
   LOGICAL, Intent(in) :: LvaryFmat (maxlayers,max_atmoswfs)

!  Type structures, FO Geometry output
!  -----------------------------------

   TYPE(VLIDORT_Geometry_FO), intent(in)  :: FOGeometry

!  Optical inputs
!  --------------

!  Solar flux

   real(ffp), intent(in) :: flux, fluxvec(maxstokes)

!  Atmosphere opticals
!  4/15/20. Version 1.5.2. GREEKMAT removed; Now we must use FMATRIX

   real(ffp), intent(in) :: extinction ( maxlayers )
   real(ffp), intent(in) :: deltaus    ( maxlayers )
   real(ffp), intent(in) :: omega      ( maxlayers )
   real(ffp), intent(in) :: truncfac   ( maxlayers )
   real(ffp), intent(in) :: fmatrix_up ( maxlayers, MAX_GEOMETRIES, 6 )
   real(ffp), intent(in) :: fmatrix_dn ( maxlayers, MAX_GEOMETRIES, 6 )

!  Thermal inputs, surface emissivity
!mick fix 4/3/2015 - fix dimension
   !real(ffp), intent(in) :: emiss ( MAX_GEOMETRIES )
   !real(ffp), intent(in) :: LS_emiss ( MAX_GEOMETRIES, max_surfacewfs )

!   Emissivity is now polarized. 12/11/17 Rob add.
!   real(ffp), intent(in) :: emiss ( MAX_USER_VZANGLES )
!   real(ffp), intent(in) :: LS_emiss ( MAX_USER_VZANGLES, max_surfacewfs )

   real(ffp), intent(in) :: bb_input ( 0:maxlayers )
   real(ffp), intent(in) :: surfbb
   real(ffp), intent(in) :: emiss    ( maxstokes, MAX_USER_VZANGLES )
   real(ffp), intent(in) :: LS_emiss ( maxstokes, MAX_USER_VZANGLES, max_surfacewfs )

!  Surface reflectivity (Could be the albedo) + linearizations
!    Surface leaving input added 8/2/16

   real(ffp), intent(in) :: reflec(maxstokes,maxstokes,MAX_GEOMETRIES)
   real(ffp), intent(in) :: slterm(maxstokes,MAX_GEOMETRIES)
   real(ffp), Intent(in) :: ls_reflec   ( maxstokes, maxstokes, MAX_GEOMETRIES, max_surfacewfs )
   real(ffp), Intent(in) :: lssl_slterm ( maxstokes, MAX_GEOMETRIES, max_sleavewfs  )

!  Linearized optical inputs
!  -------------------------

!  4/15/20. Version 1.5.2. L_GREEKMAT removed; Now we must use L_FMATRIX

   real(ffp), intent(in) :: L_extinction  ( maxlayers, max_atmoswfs )
   real(ffp), intent(in) :: L_deltaus     ( maxlayers, max_atmoswfs )
   real(ffp), intent(in) :: L_omega       ( maxlayers, max_atmoswfs )
   real(ffp), intent(in) :: L_truncfac    ( maxlayers, max_atmoswfs )

   real(ffp), intent(in) :: L_fmatrix_up  ( maxlayers, MAX_GEOMETRIES, 6, max_atmoswfs )
   real(ffp), intent(in) :: L_fmatrix_dn  ( maxlayers, MAX_GEOMETRIES, 6, max_atmoswfs )

!  Subroutine outputs
!  ==================

!  Stokes vectors
!  --------------

!  Solar

   real(ffp), intent(out) :: fo_stokes_ss ( max_user_levels,MAX_GEOMETRIES,maxstokes,max_directions )
   real(ffp), intent(out) :: fo_stokes_db ( max_user_levels,MAX_GEOMETRIES,maxstokes )

!  Thermal

   real(ffp), intent(out) :: fo_stokes_dta ( max_user_levels,MAX_GEOMETRIES,maxstokes,max_directions )
   real(ffp), intent(out) :: fo_stokes_dts ( max_user_levels,MAX_GEOMETRIES,maxstokes )

!  Composite

   real(ffp), intent(out) :: fo_stokes_atmos ( max_user_levels,MAX_GEOMETRIES,maxstokes,max_directions )
   real(ffp), intent(out) :: fo_stokes_surf  ( max_user_levels,MAX_GEOMETRIES,maxstokes )
   real(ffp), intent(out) :: fo_stokes       ( max_user_levels,MAX_GEOMETRIES,maxstokes,max_directions )

!  Jacobians
!  ---------

!  Solar

   real(ffp), intent(out) :: fo_profilewf_ss ( max_atmoswfs,maxlayers,max_user_levels,&
                                               MAX_GEOMETRIES,maxstokes,max_directions )
   real(ffp), intent(out) :: fo_profilewf_db ( max_atmoswfs,maxlayers,max_user_levels,&
                                               MAX_GEOMETRIES,maxstokes )
   real(ffp), intent(out) :: fo_surfacewf_db ( max_surfacewfs,max_user_levels,&
                                               MAX_GEOMETRIES,maxstokes )

!  Thermal

   real(ffp), intent(out) :: fo_profilewf_dta ( max_atmoswfs,maxlayers,max_user_levels,&
                                                MAX_GEOMETRIES,maxstokes,max_directions )
   real(ffp), intent(out) :: fo_profilewf_dts ( max_atmoswfs,maxlayers,max_user_levels,&
                                                MAX_GEOMETRIES,maxstokes )
   real(ffp), intent(out) :: fo_surfacewf_dts ( max_surfacewfs,max_user_levels,&
                                                MAX_GEOMETRIES,maxstokes )

!  Composite

   real(ffp), intent(out) :: fo_profilewf_atmos ( max_atmoswfs,maxlayers,max_user_levels,&
                                                  MAX_GEOMETRIES,maxstokes,max_directions )
   real(ffp), intent(out) :: fo_profilewf_surf  ( max_atmoswfs,maxlayers,max_user_levels,&
                                                  MAX_GEOMETRIES,maxstokes )
   real(ffp), intent(out) :: fo_profilewf       ( max_atmoswfs,maxlayers,max_user_levels,&
                                                  MAX_GEOMETRIES,maxstokes,max_directions )
   real(ffp), intent(out) :: fo_surfacewf       ( max_surfacewfs,max_user_levels,&
                                                  MAX_GEOMETRIES,maxstokes )

!  4/9/19. Additional output for the sleave correction

   real(ffp), Intent(out) :: CUMTRANS    ( max_user_levels, MAX_GEOMETRIES )
   real(ffp), Intent(out) :: LP_CUMTRANS ( max_user_levels, MAX_GEOMETRIES, maxlayers, max_atmoswfs )

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

   logical    :: do_Tsources_up       (MAXLAYERS,MAX_USER_VZANGLES)
   logical    :: do_Tsources_dn       (MAXLAYERS,MAX_USER_VZANGLES)
   logical    :: do_Tsources_up_p     (MAX_PARTLAYERS,MAX_USER_VZANGLES)
   logical    :: do_Tsources_dn_p     (MAX_PARTLAYERS,MAX_USER_VZANGLES)

!  RT Calculation outputs
!  ----------------------

!  SS routines output

   real(ffp)  :: stokes_up    ( max_user_levels,maxstokes,MAX_GEOMETRIES )
   real(ffp)  :: stokes_dn    ( max_user_levels,maxstokes,MAX_GEOMETRIES )
   real(ffp)  :: stokes_db    ( max_user_levels,maxstokes,MAX_GEOMETRIES )

   real(ffp)  :: LP_Jacobians_up  ( max_user_levels, maxstokes, MAX_GEOMETRIES, maxlayers, max_atmoswfs )
   real(ffp)  :: LP_Jacobians_dn  ( max_user_levels, maxstokes, MAX_GEOMETRIES, maxlayers, max_atmoswfs )
   real(ffp)  :: LP_Jacobians_db  ( max_user_levels, maxstokes, MAX_GEOMETRIES, maxlayers, max_atmoswfs )

   real(ffp)  :: LS_Jacobians_db  ( max_user_levels, maxstokes, MAX_GEOMETRIES, max_surfacewfs )

!  Thermal routines output

   real(ffp)  :: intensity_dta_up_LOS ( max_user_levels,MAX_USER_VZANGLES )
   real(ffp)  :: intensity_dta_dn_LOS ( max_user_levels,MAX_USER_VZANGLES )
   real(ffp)  :: intensity_dts_LOS    ( max_user_levels,MAX_USER_VZANGLES )

   real(ffp)  :: LP_Jacobians_dta_up_LOS  ( max_user_levels, MAX_USER_VZANGLES, maxlayers, max_atmoswfs )
   real(ffp)  :: LP_Jacobians_dta_dn_LOS  ( max_user_levels, MAX_USER_VZANGLES, maxlayers, max_atmoswfs )
   real(ffp)  :: LP_Jacobians_dts_LOS     ( max_user_levels, MAX_USER_VZANGLES, maxlayers, max_atmoswfs )
   real(ffp)  :: LS_Jacobians_dts_LOS     ( max_user_levels, MAX_USER_VZANGLES, max_surfacewfs )

!  5/22/20. Version 2.8.2 Upgrades.
!    ==> LOSTRANS output might be used again for the (upwelling) thermal-NoScattering contribution??

   real(ffp)  :: lostrans_up_LOS      ( maxlayers     , MAX_USER_VZANGLES )
   real(ffp)  :: lostrans_up_p_LOS    ( MAX_PARTLAYERS, MAX_USER_VZANGLES )
   real(ffp)  :: L_lostrans_up_LOS    ( maxlayers     , MAX_USER_VZANGLES, max_atmoswfs )
   real(ffp)  :: L_lostrans_up_p_LOS  ( MAX_PARTLAYERS, MAX_USER_VZANGLES, max_atmoswfs )

!  StokesQUV_dts_LOS, contribution from Polarized emissivity. 12/11/17 Rob add.
!   Also the Jacobian arrays LP_JacobiansQUV_dts_LOS and LS_JacobiansQUV_dts_LOS

   real(ffp)  :: StokesQUV_dts_LOS       ( max_user_levels, 3, MAX_USER_VZANGLES )
   real(ffp)  :: LP_JacobiansQUV_dts_LOS ( max_user_levels, 3, MAX_USER_VZANGLES, maxlayers, max_atmoswfs )
   real(ffp)  :: LS_JacobiansQUV_dts_LOS ( max_user_levels, 3, MAX_USER_VZANGLES, max_surfacewfs )

!  Other products
!  --------------

!  Thermal setup and linearization

   real(ffp)  :: tcom1(maxlayers,2)
   real(ffp)  :: L_tcom1(maxlayers,2,max_atmoswfs)

!  Dummies

!   real(ffp)  :: SScumsource_up     ( 0:maxlayers,maxstokes,MAX_GEOMETRIES )
!   real(ffp)  :: SScumsource_dn     ( 0:maxlayers,maxstokes,MAX_GEOMETRIES )
!   real(ffp)  :: DTcumsource_up     ( 0:maxlayers,MAX_GEOMETRIES )
!   real(ffp)  :: DTcumsource_dn     ( 0:maxlayers,MAX_GEOMETRIES )

!  LOCAL HELP VARIABLES
!  --------------------

!  help variables. 

   integer   :: ns, nv, na, g, n, spar, par, lev, o1

!mick fix 9/19/2017 - initialized "atmos" and "surf" components of both
!                     intensity & profilewf quantities

!  Initialize Intensity output. Including composites (3/9/17)

   fo_stokes_ss       = zero
   fo_stokes_db       = zero
   fo_stokes_dta      = zero
   fo_stokes_dts      = zero

   fo_stokes_atmos    = zero
   fo_stokes_surf     = zero
   fo_stokes          = zero

!  Initialize Jacobian output. Including composites (3/9/17)

   fo_profilewf_ss    = zero
   fo_profilewf_db    = zero
   fo_surfacewf_db    = zero

   fo_profilewf_dta   = zero
   fo_profilewf_dts   = zero
   fo_surfacewf_dts   = zero

   fo_profilewf_atmos = zero
   fo_profilewf_surf  = zero
   fo_profilewf       = zero
   fo_surfacewf       = zero

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
!  - 4/9/19. Add the CUMTRANS/LP_CUMTRANS output, add water-leaving control

       call FO_Vector_SSRT_ILPS_UP &
        ( do_sunlight, do_deltam_scaling, do_Lambertian, do_surface_leaving, do_water_leaving,     & ! Inputs (Flags-General/Surface)
          do_Partials, do_PlanPar, do_enhanced_ps, do_sources_up, do_sources_up_p, do_profilewfs,  & ! Inputs(Flags/sources)
          do_surfacewfs, do_sleavewfs, Lvaryflags, Lvarynums, LvaryFmat, n_reflecwfs, n_sleavewfs, & ! Inputs (Control, Jacobian)
          n_surfacewfs, nstokes, ngeoms, nlayers, n_user_levels, user_level_mask_up,               & ! Inputs (control, output)
          npartials, partial_outindex, partial_outflag, partial_layeridx, FOGeometry,              & ! Inputs (patial/Geometry)
          flux, fluxvec, extinction, deltaus, omega, truncfac, fmatrix_up, reflec, slterm,         & ! Inputs (Optical)
          L_extinction, L_deltaus, L_omega, L_truncfac, L_fmatrix_up, LS_reflec, LSSL_slterm,      & ! Inputs (Linearized)
          Stokes_up, Stokes_db, LP_Jacobians_up, LP_Jacobians_db, LS_Jacobians_db, cumtrans, LP_cumtrans )   ! Output

!  Save results
!mick mod 9/19/2017 - turned off "fo_stokes", "fo_profilewf", & "fo_surfacewf" (defined later)

       do o1 = 1, nstokes
         do g = 1, ngeoms
           do lev=1,n_user_levels
             fo_stokes_ss(lev,g,o1,upidx) = stokes_up(lev,o1,g)
             fo_stokes_db(lev,g,o1)       = stokes_db(lev,o1,g)
           enddo
         enddo
       enddo

       if ( do_profilewfs ) then
         do o1 = 1, nstokes
           do g = 1, ngeoms
             do lev=1,n_user_levels
               do n = 1, nlayers
                 if ( Lvaryflags(n) ) then
                   do par=1,Lvarynums(n)
                     fo_profilewf_ss(par,n,lev,g,o1,upidx) = LP_Jacobians_up(lev,o1,g,n,par)
                     fo_profilewf_db(par,n,lev,g,o1)       = LP_Jacobians_db(lev,o1,g,n,par)
                   enddo
                 endif
               enddo
             enddo
           enddo
         enddo
       endif

       if ( do_surfacewfs ) then
         do o1 = 1, nstokes
           do g = 1, ngeoms
             do lev=1,n_user_levels
               do spar=1,n_surfacewfs
                 fo_surfacewf_db(spar,lev,g,o1) = LS_Jacobians_db(lev,o1,g,spar)

               enddo
             enddo
           enddo
         enddo
       endif

! End upwelling

     end if

!  Donwelling

     if ( do_dnwelling ) then

!  Sources, temporary fix until criticality realized. 9/17/16

       do_sources_dn   = .true. ; do_sources_dn_p = .true.

! RT call - solar only. Updated 9/17/16.

       call FO_Vector_SSRT_ILPS_DN &
        ( do_sunlight, do_deltam_scaling, do_Partials, do_PlanPar, do_enhanced_ps,    & ! Inputs (Flags/flux)
          do_sources_dn, do_sources_dn_p, do_profilewfs, Lvaryflags, Lvarynums,       & ! Inputs (Control, Lin )
          LvaryFmat, nstokes, ngeoms, nlayers, n_user_levels, user_level_mask_dn,     & ! Inputs (Control, Output)
          npartials, partial_outindex, partial_outflag, partial_layeridx, FOGeometry, & ! Inputs (Control, Partial)
          flux, fluxvec, extinction, deltaus, omega, truncfac, fmatrix_dn,            & ! Inputs (Optical)
          L_extinction, L_deltaus, L_omega, L_truncfac, L_fmatrix_dn,                 & ! Inputs (Optical - Lin)
          Stokes_dn, LP_Jacobians_dn )                                                 ! Output

!  Save results
!mick mod 9/19/2017 - turned off "fo_stokes" & "fo_profilewf" (defined later)

       do o1 = 1, nstokes
         do g = 1, ngeoms
           do lev=1,n_user_levels
             fo_stokes_ss(lev,g,o1,dnidx) = stokes_dn(lev,o1,g)
           enddo
         enddo
       enddo

       if ( do_profilewfs ) then
         do o1 = 1, nstokes
           do g = 1, ngeoms
             do lev=1,n_user_levels
               do n = 1, nlayers
                 if ( Lvaryflags(n) ) then
                   do par=1,Lvarynums(n)
                     fo_profilewf_ss(par,n,lev,g,o1,dnidx) = LP_Jacobians_dn(lev,o1,g,n,par)
                   enddo
                 endif
               enddo
             enddo
           enddo
         enddo
       endif

!  End downwelling

     endif

!  End solar run

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
!    -- Arrays "Emiss" and "LS_emiss" now have vector dimension.

!  5/22/20. Version 2.8.2 Upgrades.
!        ==> Add lostrans_up_LOS   + lostrans_up_p_LOS to argument list
!        ==> Add L_lostrans_up_LOS + L_lostrans_up_p_LOS to argument list, similarly

       call FO_Thermal_DTRT_ILPS_UP &
         ( do_deltam_scaling, do_Partials, do_PlanPar, do_enhanced_ps,            & ! Inputs (Flags)
           Do_Polarized_Emissivity, do_Tsources_up, do_Tsources_up_p,             & ! Inputs (Flags)
           do_profilewfs, do_surfacewfs, Lvaryflags, Lvarynums, n_surfacewfs,     & ! Inputs (Control, Jacobians)
           nstokes, nvzas, nlayers, n_user_levels, user_level_mask_up, npartials, & ! Inputs (Control output)
           partial_outindex, partial_outflag, partial_layeridx, FOGeometry,       & ! Inputs (Partial/Geometry)
           extinction, deltaus, omega, truncfac, bb_input,                        & ! Inputs (Optical/thermal)
           surfbb, emiss(1,:), emiss(2:4,:),                                      & ! Inputs (Surface)
           L_extinction, L_deltaus, L_omega, L_truncfac,                          & ! Inputs (Optical - Linearized)
           LS_emiss(1,:,:), LS_emiss(2:4,:,:),                                    & ! Inputs (Surface - Linearized)
           intensity_dta_up_LOS, intensity_dts_LOS, LP_Jacobians_dta_up_LOS,      & ! Main Outputs
           LP_Jacobians_dts_LOS, LS_Jacobians_dts_LOS, tcom1, L_tcom1,                 & ! Main Outputs
           lostrans_up_LOS, lostrans_up_p_LOS, L_lostrans_up_LOS, L_lostrans_up_p_LOS, & ! Other Outputs
           StokesQUV_dts_LOS, LP_JacobiansQUV_dts_LOS, LS_JacobiansQUV_dts_LOS )         ! Optional Output. 12/11/17 Rob Add.

!  Save results
!mick mod 9/19/2017 - turned off "fo_stokes", "fo_profilewf", & "fo_surfacewf" (defined later)
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

       if ( do_profilewfs ) then
         if ( do_ObsGeom ) then
           do g = 1, ngeoms
             do lev=1,n_user_levels
               do n = 1, nlayers
                 if ( Lvaryflags(n) ) then
                   do par=1,Lvarynums(n)
                     fo_profilewf_dta(par,n,lev,g,o1,upidx) = LP_Jacobians_dta_up_LOS(lev,g,n,par)
                     fo_profilewf_dts(par,n,lev,g,o1) = LP_Jacobians_dts_LOS(lev,g,n,par)
                   enddo
                 endif
               enddo
             enddo
           enddo
         else if ( do_doublet ) then
           do nv = 1, nvzas ; do ns = 1, nszas
             g = nd_offset(ns) + nv
             do lev=1,n_user_levels
               do n = 1, nlayers
                 if ( Lvaryflags(n) ) then
                   do par=1,Lvarynums(n)
                     fo_profilewf_dta(par,n,lev,g,o1,upidx) = LP_Jacobians_dta_up_LOS(lev,nv,n,par)
                     fo_profilewf_dts(par,n,lev,g,o1) = LP_Jacobians_dts_LOS(lev,nv,n,par)
                   enddo
                 endif
               enddo
             enddo
           enddo ; enddo
         else
           do nv = 1, nvzas ; do ns = 1, nszas ; do na = 1, nazms
             g = na_offset(ns,nv) + na
             do lev=1,n_user_levels
               do n = 1, nlayers
                 if ( Lvaryflags(n) ) then
                   do par=1,Lvarynums(n)
                     fo_profilewf_dta(par,n,lev,g,o1,upidx) = LP_Jacobians_dta_up_LOS(lev,nv,n,par)
                     fo_profilewf_dts(par,n,lev,g,o1) = LP_Jacobians_dts_LOS(lev,nv,n,par)
                   enddo
                 endif
               enddo
             enddo
           enddo ; enddo ; enddo
         endif
       endif

       if ( do_surfacewfs ) then
         if ( do_ObsGeom ) then
           do g = 1, ngeoms
             do lev=1,n_user_levels
               do spar=1,n_surfacewfs
                 fo_surfacewf_dts(spar,lev,g,o1) = LS_Jacobians_dts_LOS(lev,g,spar)
               enddo
             enddo
           enddo
         else if ( do_doublet ) then
           do nv = 1, nvzas ; do ns = 1, nszas
             g = nd_offset(ns) + nv
             do spar=1,n_surfacewfs
               fo_surfacewf_dts(spar,lev,g,o1) = LS_Jacobians_dts_LOS(lev,nv,spar)
             enddo
           enddo ; enddo
         else
           do nv = 1, nvzas ; do ns = 1, nszas ; do na = 1, nazms
             g = na_offset(ns,nv) + na
!mick fix 3/22/2017 - added lev loop
             do lev=1,n_user_levels
               do spar=1,n_surfacewfs
                 fo_surfacewf_dts(spar,lev,g,o1) = LS_Jacobians_dts_LOS(lev,nv,spar)
               enddo
             enddo
           enddo ; enddo ; enddo
         endif
       endif

!  Save polarized Emissivity results. 12/11/17  Rob add.
!  ----------------------------------------------------

!  4/15/20. Version 1.5.2.  Add the Doublet goemetry option

       if ( do_Polarized_Emissivity.and.nstokes.gt.1 ) then

!  Stokes vector

         if ( do_obsgeom ) then
           do g = 1, nvzas
             do lev=1,n_user_levels
               fo_stokes_dts(lev,g,2:nstokes) = StokesQUV_dts_LOS(lev,1:nstokes-1,g)
             enddo
           enddo
         else
           do nv = 1, nvzas
             do ns = 1, nszas
               do na = 1, nazms
                 g = na_offset(ns,nv) + na
                 do lev=1,n_user_levels
                   fo_stokes_dts(lev,g,2:nstokes) = StokesQUV_dts_LOS(lev,1:nstokes-1,nv)
                 enddo
               enddo
             enddo
           enddo
         endif

!  Profile Jacobians

         if ( do_profilewfs ) then
           if ( do_ObsGeom ) then
             do g = 1, ngeoms
               do lev=1,n_user_levels
                 do n = 1, nlayers
                   if ( Lvaryflags(n) ) then
                     do par=1,Lvarynums(n)
                       fo_profilewf_dts(par,n,lev,g,2:nstokes) = LP_JacobiansQUV_dts_LOS(lev,1:nstokes-1,g,n,par)
                     enddo
                   endif
                 enddo
               enddo
             enddo
           else
             do nv = 1, nvzas
               do ns = 1, nszas
                 do na = 1, nazms
                   g = na_offset(ns,nv) + na
                   do lev=1,n_user_levels
                     do n = 1, nlayers
                       if ( Lvaryflags(n) ) then
                         do par=1,Lvarynums(n)
                           fo_profilewf_dts(par,n,lev,g,2:nstokes) = LP_JacobiansQUV_dts_LOS(lev,1:nstokes-1,nv,n,par)
                         enddo
                       endif
                     enddo
                   enddo
                 enddo
               enddo
             enddo
           endif
         endif

!  Surface Jacobians

         if ( do_surfacewfs ) then
           if ( do_ObsGeom ) then
             do g = 1, ngeoms
               do lev=1,n_user_levels
                 do spar=1,n_surfacewfs
                   fo_surfacewf_dts(spar,lev,g,2:nstokes) = LS_JacobiansQUV_dts_LOS(lev,1:nstokes-1,g,spar)
                 enddo
               enddo
             enddo
           else
             do nv = 1, nvzas
               do ns = 1, nszas
                 do na = 1, nazms
                   g = na_offset(ns,nv) + na
                   do lev=1,n_user_levels
                     do spar=1,n_surfacewfs
                       fo_surfacewf_dts(spar,lev,g,2:nstokes) = LS_JacobiansQUV_dts_LOS(lev,1:nstokes-1,nv,spar)
                     enddo
                   enddo
                 enddo
               enddo
             enddo
           endif
         endif

!  End polarized emissivity clause

       endif

!  End upwelling

     endif

!  Downwelling
!  -----------

     if ( do_dnwelling ) then

!  Sources, temporary fix until criticality realized. 9/17/16

       do_Tsources_dn   = .true. ; do_Tsources_dn_p = .true.

!  Direct thermal, calculate. Updated 9/17/16

        call FO_Thermal_DTRT_ILPS_DN &
         ( do_deltam_scaling, do_Partials, do_PlanPar, do_enhanced_ps,                 & ! Inputs (Flags)
           do_Tsources_dn, do_Tsources_dn_p, do_profilewfs, Lvaryflags, Lvarynums,     & ! Inputs (Flags/Jac-control)
           nvzas, nlayers, n_user_levels, user_level_mask_dn,                          & ! Inputs (control output)
           npartials, partial_outindex, partial_outflag, partial_layeridx, FOGeometry, & ! Inputs (control-partial)
           bb_input, extinction, deltaus, omega, truncfac,                             & ! Inputs (Optical - Regular)
           L_extinction, L_deltaus, L_omega, L_truncfac,                               & ! Inputs (Optical - Linearized)
           intensity_dta_dn_LOS, LP_Jacobians_dta_dn_LOS, tcom1, L_tcom1 )               ! Output

!  Save results
!mick mod 9/19/2017 - turned off "fo_stokes" & "fo_profilewf" (defined later)
!  4/15/20. Version 1.5.2.  Add the Doublet goemetry option

       o1=1
       if ( do_obsgeom ) then
         do g = 1, nvzas
           do lev=1,n_user_levels
             fo_stokes_dta(lev,g,o1,dnidx) = intensity_dta_dn_LOS(lev,g)
           enddo
         enddo
       else
         do nv = 1, nvzas
           do ns = 1, nszas
             do na = 1, nazms
               g = na_offset(ns,nv) + na
               do lev=1,n_user_levels
                 fo_stokes_dta(lev,g,o1,dnidx) = intensity_dta_dn_LOS(lev,nv)
               enddo
             enddo
           enddo
         enddo
       endif

       if ( do_profilewfs ) then
          if (do_ObsGeom ) then
           do g = 1, ngeoms
             do lev=1,n_user_levels
               do n = 1, nlayers
                 if ( Lvaryflags(n) ) then
                   do par=1,Lvarynums(n)
                     fo_profilewf_dta(par,n,lev,g,o1,dnidx) = LP_Jacobians_dta_dn_LOS(lev,g,n,par)
                   enddo
                 endif
               enddo
             enddo
           enddo
         else if ( do_doublet ) then
           do nv = 1, nvzas ; do ns = 1, nszas
             g = nd_offset(ns) + nv
             do lev=1,n_user_levels
               do n = 1, nlayers
                 if ( Lvaryflags(n) ) then
                   do par=1,Lvarynums(n)
                     fo_profilewf_dta(par,n,lev,g,o1,dnidx) = LP_Jacobians_dta_dn_LOS(lev,nv,n,par)
                   enddo
                 endif
               enddo
             enddo
           enddo ; enddo
         else
           do nv = 1, nvzas ; do ns = 1, nszas ; do na = 1, nazms
             g = na_offset(ns,nv) + na
             do lev=1,n_user_levels
               do n = 1, nlayers
                 if ( Lvaryflags(n) ) then
                   do par=1,Lvarynums(n)
                     fo_profilewf_dta(par,n,lev,g,o1,dnidx) = LP_Jacobians_dta_dn_LOS(lev,nv,n,par)
                   enddo
                 endif
               enddo
             enddo
           enddo ; enddo ; enddo
         endif
       endif

!  end downwelling

      endif

!  End Thermal run

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
     if ( do_profilewfs ) then
       do o1 = 1, nstokes
         do g = 1, ngeoms
           do lev=1,n_user_levels
             do n = 1, nlayers
               if ( Lvaryflags(n) ) then
                 do par=1,Lvarynums(n)
                   fo_profilewf_atmos(par,n,lev,g,o1,upidx) = fo_profilewf_ss(par,n,lev,g,o1,upidx) &
                                                            + fo_profilewf_dta(par,n,lev,g,o1,upidx)
                   fo_profilewf_surf(par,n,lev,g,o1)        = fo_profilewf_db(par,n,lev,g,o1) &
                                                            + fo_profilewf_dts(par,n,lev,g,o1)
                   fo_profilewf(par,n,lev,g,o1,upidx)       = fo_profilewf_atmos(par,n,lev,g,o1,upidx) &
                                                            + fo_profilewf_surf(par,n,lev,g,o1)
                 enddo
               endif
             enddo
           enddo
         enddo
       enddo
     endif
     if ( do_surfacewfs ) then
       do o1 = 1, nstokes
         do g = 1, ngeoms
           do lev=1,n_user_levels
             do spar=1,n_surfacewfs
               fo_surfacewf(spar,lev,g,o1) = fo_surfacewf_db(spar,lev,g,o1) + fo_surfacewf_dts(spar,lev,g,o1)
             enddo
           enddo
         enddo
       enddo
     endif
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
     if ( do_profilewfs ) then
       do o1 = 1, nstokes
         do g = 1, ngeoms
           do lev=1,n_user_levels
             do n = 1, nlayers
               if ( Lvaryflags(n) ) then
                 do par=1,Lvarynums(n)
                     fo_profilewf_atmos(par,n,lev,g,o1,dnidx) = fo_profilewf_ss(par,n,lev,g,o1,dnidx) &
                                                              + fo_profilewf_dta(par,n,lev,g,o1,dnidx)
                     fo_profilewf(par,n,lev,g,o1,dnidx)       = fo_profilewf_atmos(par,n,lev,g,o1,dnidx)
                 enddo
               endif
             enddo
           enddo
         enddo
       enddo
     endif
   endif

!  Finish

   return
end subroutine VFO_RTCALC_LPS_MASTER

!

subroutine VFO_RTCALC_LCS_MASTER &
       ( do_solar_sources, do_thermal_emission, do_surface_emission, do_Polarized_Emissivity, & ! Input flags (sources)
         do_upwelling, do_dnwelling, do_lambertian, do_surface_leaving, do_water_leaving,     & ! Input flags (dir/surface)
         do_sunlight, do_deltam_scaling, do_obsgeom, do_doublet, do_Partials, do_planpar,     & ! Input flags (surface/geoms)
         do_enhanced_ps, do_columnwfs, do_surfacewfs, do_sleavewfs,                           & ! Input Lin flags
         nstokes, ngeoms, nszas, nvzas, nazms, na_offset, nd_offset, nlayers,                 & ! Input Numbers/Offsets
         n_user_levels, user_level_mask_up, user_level_mask_dn,                               & ! Input control-levels
         npartials, partial_outindex, partial_outflag, partial_layeridx, FOGeometry,          & ! Input control-partial/Geometry
         n_columnwfs, LvaryFmat, n_reflecwfs, n_sleavewfs, n_surfacewfs,                      & ! Input Lin control
         flux, fluxvec, extinction, deltaus, omega, truncfac, fmatrix_up, fmatrix_dn,         & ! Input atmos optical
         bb_input, surfbb, emiss, reflec, slterm,                                             & ! Input thermal/surface
         L_extinction, L_deltaus, L_omega, L_truncfac, L_fmatrix_up, L_fmatrix_dn,            & ! Input Lin atmos optical
         LS_emiss, LS_reflec, LSSL_slterm,                                                    & ! Input Lin thermal/surf
         fo_stokes_ss, fo_stokes_db, fo_stokes_dta, fo_stokes_dts,                            & ! Output - Stokes   Vectors
         fo_columnwf_ss,  fo_columnwf_db,  fo_columnwf_dta, fo_columnwf_dts,                  & ! Output - Column Jacobians
         fo_surfacewf_db,  fo_surfacewf_dts,                                                  & ! Output - Surface  Jacobians
         fo_stokes_atmos,  fo_stokes_surf, fo_stokes,                                         & ! Output - Stokes   composites
         fo_columnwf_atmos, fo_columnwf_surf, fo_columnwf, fo_surfacewf,                      & ! Output - Jacobian composites
         cumtrans, LC_cumtrans )                                                                ! Output - ancillary

!  4/9/19. Add CUMTRANS, LC_Cumtrans output, and Waterleaving input control
!  4/15/20. Version 1.5.2. Add doublet flag and offsets

!  Use modules

   USE FO_Vector_SSRT_ILCS_m
   USE FO_Thermal_DTRT_ILCS_m

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
   integer, intent(in)  :: user_level_mask_up ( MAX_USER_LEVELS )
   integer, intent(in)  :: user_level_mask_dn ( MAX_USER_LEVELS )

!  Control for partial-layer output, added 8/25/16

   logical, intent(in)  :: do_Partials
   integer, intent(in)  :: Npartials
   integer, intent(in)  :: partial_layeridx( max_partlayers  )
   logical, intent(in)  :: partial_outflag ( MAX_USER_LEVELS )
   integer, intent(in)  :: partial_outindex( MAX_USER_LEVELS )

!  Jacobian Flags. do_sleavewfs added 8/2/16

   LOGICAL, Intent(in) :: do_columnwfs
   LOGICAL, Intent(in) :: do_surfacewfs
   LOGICAL, Intent(in) :: do_sleavewfs

!  Jacobian control. Reflec and sleave numbers added, 8/2/16
!    Note that n_surfacewfs = n_reflecwfs + n_sleavewfs

   INTEGER, Intent(in) :: n_reflecwfs
   INTEGER, Intent(in) :: n_sleavewfs
   INTEGER, Intent(in) :: n_surfacewfs
   INTEGER, Intent(in) :: n_columnwfs
   LOGICAL, Intent(in) :: LvaryFmat (maxlayers,max_atmoswfs)

!  Type structures, FO Geometry output
!  -----------------------------------

   TYPE(VLIDORT_Geometry_FO), intent(in)  :: FOGeometry

!  Optical inputs
!  --------------

!  Solar flux

   real(ffp), intent(in) :: flux, fluxvec(maxstokes)

!  Atmosphere opticals
!  4/15/20. Version 1.5.2. GREEKMAT removed; Now we must use FMATRIX

   real(ffp), intent(in) :: extinction  ( maxlayers )
   real(ffp), intent(in) :: deltaus     ( maxlayers )
   real(ffp), intent(in) :: omega       ( maxlayers )
   real(ffp), intent(in) :: truncfac ( maxlayers )
   real(ffp), intent(in) :: fmatrix_up  ( maxlayers, MAX_GEOMETRIES, 6 )
   real(ffp), intent(in) :: fmatrix_dn  ( maxlayers, MAX_GEOMETRIES, 6 )

!  Thermal inputs, surface emissivity
!    Emissivity is now polarized. 12/11/17 Rob add.

   real(ffp), intent(in) :: bb_input ( 0:maxlayers )
   real(ffp), intent(in) :: surfbb
   real(ffp), intent(in) :: emiss    ( maxstokes, MAX_USER_VZANGLES )
   real(ffp), intent(in) :: LS_emiss ( maxstokes, MAX_USER_VZANGLES, max_surfacewfs )

!  Surface reflectivity (Could be the albedo) + linearizations
!    Surface leaving input added 8/2/16

   real(ffp), intent(in) :: reflec(maxstokes,maxstokes,MAX_GEOMETRIES)
   real(ffp), intent(in) :: slterm(maxstokes,MAX_GEOMETRIES)
   real(ffp), Intent(in) :: ls_reflec   ( maxstokes, maxstokes, MAX_GEOMETRIES, max_surfacewfs )
   real(ffp), Intent(in) :: lssl_slterm ( maxstokes, MAX_GEOMETRIES, max_sleavewfs  )

!  Linearized optical inputs
!  -------------------------

!  4/15/20. Version 1.5.2. L_GREEKMAT removed; Now we must use L_FMATRIX

   real(ffp), intent(in) :: L_extinction  ( maxlayers, max_atmoswfs )
   real(ffp), intent(in) :: L_deltaus     ( maxlayers, max_atmoswfs )
   real(ffp), intent(in) :: L_omega       ( maxlayers, max_atmoswfs )
   real(ffp), intent(in) :: L_truncfac    ( maxlayers, max_atmoswfs )

   real(ffp), intent(in) :: L_fmatrix_up  ( maxlayers, MAX_GEOMETRIES, 6, max_atmoswfs )
   real(ffp), intent(in) :: L_fmatrix_dn  ( maxlayers, MAX_GEOMETRIES, 6, max_atmoswfs )

!  Subroutine outputs
!  ==================

!  Stokes Vectors
!  --------------

!  Solar

   real(ffp), intent(out) :: fo_stokes_ss ( max_user_levels,MAX_GEOMETRIES,maxstokes,max_directions )
   real(ffp), intent(out) :: fo_stokes_db ( max_user_levels,MAX_GEOMETRIES,maxstokes )

!  Thermal

   real(ffp), intent(out) :: fo_stokes_dta ( max_user_levels,MAX_GEOMETRIES,maxstokes,max_directions )
   real(ffp), intent(out) :: fo_stokes_dts ( max_user_levels,MAX_GEOMETRIES,maxstokes )

!  Composite

   real(ffp), intent(out) :: fo_stokes_atmos ( max_user_levels,MAX_GEOMETRIES,maxstokes,max_directions )
   real(ffp), intent(out) :: fo_stokes_surf  ( max_user_levels,MAX_GEOMETRIES,maxstokes )
   real(ffp), intent(out) :: fo_stokes       ( max_user_levels,MAX_GEOMETRIES,maxstokes,max_directions )

!  Jacobians
!  ---------

!  Solar

   real(ffp), intent(out) :: fo_columnwf_ss  ( max_atmoswfs,max_user_levels,MAX_GEOMETRIES,maxstokes,max_directions )
   real(ffp), intent(out) :: fo_columnwf_db  ( max_atmoswfs,max_user_levels,MAX_GEOMETRIES,maxstokes )
   real(ffp), intent(out) :: fo_surfacewf_db ( max_surfacewfs,max_user_levels,MAX_GEOMETRIES,maxstokes )

!  Thermal

   real(ffp), intent(out) :: fo_columnwf_dta  ( max_atmoswfs,max_user_levels,MAX_GEOMETRIES,maxstokes,max_directions )
   real(ffp), intent(out) :: fo_columnwf_dts  ( max_atmoswfs,max_user_levels,MAX_GEOMETRIES,maxstokes )
   real(ffp), intent(out) :: fo_surfacewf_dts ( max_surfacewfs,max_user_levels,MAX_GEOMETRIES,maxstokes )

!  Composite

   real(ffp), intent(out) :: fo_columnwf_atmos  ( max_atmoswfs,max_user_levels,MAX_GEOMETRIES,maxstokes,max_directions )
   real(ffp), intent(out) :: fo_columnwf_surf   ( max_atmoswfs,max_user_levels,MAX_GEOMETRIES,maxstokes )
   real(ffp), intent(out) :: fo_columnwf        ( max_atmoswfs,max_user_levels,MAX_GEOMETRIES,maxstokes,max_directions )
   real(ffp), intent(out) :: fo_surfacewf       ( max_surfacewfs,max_user_levels,MAX_GEOMETRIES,maxstokes )

!  4/9/19. Additional output for the sleave correction

   real(ffp), Intent(out) :: CUMTRANS    ( max_user_levels, MAX_GEOMETRIES )
   real(ffp), Intent(out) :: LC_CUMTRANS ( max_user_levels, MAX_GEOMETRIES,max_atmoswfs )

!  Other variables
!  ===============

!  Existence flags. 8/19/16. Criticality enters here

   logical    :: do_sources_up       (maxlayers,MAX_GEOMETRIES)
   logical    :: do_sources_dn       (maxlayers,MAX_GEOMETRIES)
   logical    :: do_sources_up_p     (MAX_PARTLAYERS,MAX_GEOMETRIES)
   logical    :: do_sources_dn_p     (MAX_PARTLAYERS,MAX_GEOMETRIES)

!  Existence flags (Thermal). 8/19/16. Criticality enters here

   logical    :: do_Tsources_up       (MAXLAYERS,MAX_USER_VZANGLES)
   logical    :: do_Tsources_dn       (MAXLAYERS,MAX_USER_VZANGLES)
   logical    :: do_Tsources_up_p     (MAX_PARTLAYERS,MAX_USER_VZANGLES)
   logical    :: do_Tsources_dn_p     (MAX_PARTLAYERS,MAX_USER_VZANGLES)


!  RT Calculation outputs
!  ----------------------

!  SS routines output

   real(ffp)  :: stokes_up    ( max_user_levels,maxstokes,MAX_GEOMETRIES )
   real(ffp)  :: stokes_dn    ( max_user_levels,maxstokes,MAX_GEOMETRIES )
   real(ffp)  :: stokes_db    ( max_user_levels,maxstokes,MAX_GEOMETRIES )

   real(ffp)  :: LC_Jacobians_up  ( max_user_levels, maxstokes, MAX_GEOMETRIES, max_atmoswfs )
   real(ffp)  :: LC_Jacobians_dn  ( max_user_levels, maxstokes, MAX_GEOMETRIES, max_atmoswfs )
   real(ffp)  :: LC_Jacobians_db  ( max_user_levels, maxstokes, MAX_GEOMETRIES, max_atmoswfs )

   real(ffp)  :: LS_Jacobians_db  ( max_user_levels, maxstokes, MAX_GEOMETRIES, max_surfacewfs )

!  LOS VARIABLES (THERMAL SOLUTION)
!  --------------------------------

   real(ffp)  :: intensity_dta_up_LOS ( max_user_levels,MAX_USER_VZANGLES )
   real(ffp)  :: intensity_dta_dn_LOS ( max_user_levels,MAX_USER_VZANGLES )
   real(ffp)  :: intensity_dts_LOS    ( max_user_levels,MAX_USER_VZANGLES )

   real(ffp)  :: LC_Jacobians_dta_up_LOS  ( max_user_levels, MAX_USER_VZANGLES, max_atmoswfs )
   real(ffp)  :: LC_Jacobians_dta_dn_LOS  ( max_user_levels, MAX_USER_VZANGLES, max_atmoswfs )
   real(ffp)  :: LC_Jacobians_dts_LOS     ( max_user_levels, MAX_USER_VZANGLES, max_atmoswfs )
   real(ffp)  :: LS_Jacobians_dts_LOS     ( max_user_levels, MAX_USER_VZANGLES, max_surfacewfs )

!  5/22/20. Version 2.8.2 Upgrades.
!    ==> LOSTRANS output might be used again for the (upwelling) thermal-NoScattering contribution??

   real(ffp)  :: lostrans_up_LOS      ( maxlayers     , MAX_USER_VZANGLES )
   real(ffp)  :: lostrans_up_p_LOS    ( MAX_PARTLAYERS, MAX_USER_VZANGLES )
   real(ffp)  :: L_lostrans_up_LOS    ( maxlayers     , MAX_USER_VZANGLES, max_atmoswfs )
   real(ffp)  :: L_lostrans_up_p_LOS  ( MAX_PARTLAYERS, MAX_USER_VZANGLES, max_atmoswfs )

!  StokesQUV_dts_LOS, contribution from Polarized emissivity. 12/11/17 Rob add.
!   Also the Jacobian arrays LC_JacobiansQUV_dts_LOS and LS_JacobiansQUV_dts_LOS

   real(ffp)  :: StokesQUV_dts_LOS       ( max_user_levels, 3, MAX_USER_VZANGLES )
   real(ffp)  :: LC_JacobiansQUV_dts_LOS ( max_user_levels, 3, MAX_USER_VZANGLES, max_atmoswfs )
   real(ffp)  :: LS_JacobiansQUV_dts_LOS ( max_user_levels, 3, MAX_USER_VZANGLES, max_surfacewfs )

!  Other products
!  --------------

!  Thermal setup and linearization

   real(ffp)  :: tcom1(maxlayers,2)
   real(ffp)  :: L_tcom1(maxlayers,2,max_atmoswfs)

!  LOCAL HELP VARIABLES
!  --------------------

!  help variables.

   integer   :: ns, nv, na, g, par, spar, lev, o1

!mick fix 9/19/2017 - initialized "atmos" and "surf" components of both
!                     intensity & columnwf quantities

!  Initialize Intensity output. Including composites (3/9/17)

   fo_stokes_ss      = zero
   fo_stokes_db      = zero
   fo_stokes_dta     = zero
   fo_stokes_dts     = zero

   fo_stokes_atmos   = zero
   fo_stokes_surf    = zero
   fo_stokes         = zero

!  Initialize Jacobian output. Including composites (3/9/17)

   fo_columnwf_ss    = zero
   fo_columnwf_db    = zero
   fo_surfacewf_db   = zero

   fo_columnwf_dta   = zero
   fo_columnwf_dts   = zero
   fo_surfacewf_dts  = zero

   fo_columnwf_atmos = zero
   fo_columnwf_surf  = zero
   fo_columnwf       = zero
   fo_surfacewf      = zero

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

!  RT Call Solar only. Updated 9/17/16.
!  - 4/9/19. Add the CUMTRANS output, add water-leaving control

       call FO_Vector_SSRT_ILCS_UP &
         ( do_sunlight, do_deltam_scaling, do_Lambertian, do_surface_leaving, do_water_leaving,    & ! Inputs (Flags-General/Surface)
           do_Partials, do_PlanPar, do_enhanced_ps, do_sources_up, do_sources_up_p,                & ! Inputs(Flags/criticality)
           do_columnwfs, do_surfacewfs, do_sleavewfs, n_columnwfs, LvaryFmat, n_reflecwfs,         & ! Inputs (Control, Jacobian)
           n_sleavewfs, n_surfacewfs, nstokes, ngeoms, nlayers, n_user_levels, user_level_mask_up, & ! Inputs (control, output)
           npartials, partial_outindex, partial_outflag, partial_layeridx, FOGeometry,             & ! Inputs (patial/Geometry)
           flux, fluxvec, extinction, deltaus, omega, truncfac, fmatrix_up, reflec, slterm,        & ! Inputs (Optical)
           L_extinction, L_deltaus, L_omega, L_truncfac, L_fmatrix_up, LS_reflec, LSSL_slterm,     & ! Inputs (Linearized)
           Stokes_up, Stokes_db, LC_Jacobians_up, LC_Jacobians_db, LS_Jacobians_db, cumtrans, LC_cumtrans )  ! Output &

!  Save results
!mick mod 9/19/2017 - turned off "fo_stokes", "fo_columnwf", & "fo_surfacewf" (defined later)

       do o1 = 1, nstokes
         do g = 1, ngeoms
           do lev=1,n_user_levels
             fo_stokes_ss(lev,g,o1,upidx) = stokes_up(lev,o1,g)
             fo_stokes_db(lev,g,o1)       = stokes_db(lev,o1,g)
           enddo
         enddo
       enddo

       if ( do_columnwfs ) then
         do o1 = 1, nstokes
           do g = 1, ngeoms
             do lev=1,n_user_levels
               do par=1,n_columnwfs
                 fo_columnwf_ss(par,lev,g,o1,upidx) = LC_Jacobians_up(lev,o1,g,par)
                 fo_columnwf_db(par,lev,g,o1)       = LC_Jacobians_db(lev,o1,g,par)
               enddo
             enddo
           enddo
         enddo
       endif

       if ( do_surfacewfs ) then
         do o1 = 1, nstokes
           do g = 1, ngeoms
             do lev=1,n_user_levels
               do spar=1,n_surfacewfs
                 fo_surfacewf_db(spar,lev,g,o1) = LS_Jacobians_db(lev,o1,g,spar)
               enddo
             enddo
           enddo
         enddo
       endif

! End upwelling

     end if

!  Donwelling

     if ( do_dnwelling ) then

!  Sources, temporary fix until criticality realized. 9/17/16

       do_sources_dn   = .true. ; do_sources_dn_p = .true.

! RT call - solar only. Updated 9/17/16.

       call FO_Vector_SSRT_ILCS_DN &
         ( do_sunlight, do_deltam_scaling, do_Partials, do_PlanPar, do_enhanced_ps,    & ! Inputs (Flags/flux)
           do_sources_dn, do_sources_dn_p, do_columnwfs, n_columnwfs, LvaryFmat,       & ! Inputs (control, Jacobian )
           nstokes, ngeoms, nlayers, n_user_levels, user_level_mask_dn,                & ! Inputs (control)
           npartials, partial_outindex, partial_outflag, partial_layeridx, FOGeometry, & ! Inputs (control-partial)
           flux, fluxvec, extinction, deltaus, omega, truncfac, fmatrix_dn,            & ! Inputs (Optical)
           L_extinction, L_deltaus, L_omega, L_truncfac, L_fmatrix_dn,                 & ! Inputs (Optical - Linearized)
           Stokes_dn, LC_Jacobians_dn )                                                  ! Output

!  Save results
!mick mod 9/19/2017 - turned off "fo_stokes" & "fo_columnwf" (defined later)

       do o1 = 1, nstokes
         do g = 1, ngeoms
           do lev=1,n_user_levels
             fo_stokes_ss(lev,g,o1,dnidx) = stokes_dn(lev,o1,g)
           enddo
         enddo
       enddo

       if ( do_columnwfs ) then
         do o1 = 1, nstokes
           do g = 1, ngeoms
             do lev=1,n_user_levels
               do par=1,n_columnwfs
                 fo_columnwf_ss(par,lev,g,o1,dnidx) = LC_Jacobians_dn(lev,o1,g,par)
               enddo
             enddo
           enddo
         enddo
       endif

!  End downwelling

     endif

!  End solar run

   endif

!  Thermal sources run
!  -------------------

   if ( do_thermal_emission.and.do_surface_emission ) then


!  Upwelling
!  ---------

     if ( do_upwelling ) then

!  Sources, temporary fix until criticality realized. 9/17/16

       do_Tsources_up   = .true. ; do_Tsources_up_p = .true.

!  Direct thermal, calculate. Updated 9/17/16.
!    -- If Polarized Emissivity flag present, then use optional call. 12/11/17 Rob add.
!    -- Array "Emiss" now has vector dimension.

       call FO_Thermal_DTRT_ILCS_UP &
         ( do_deltam_scaling, do_Partials, do_PlanPar, do_enhanced_ps,                         & ! Inputs (Flags)
           Do_Polarized_Emissivity, do_Tsources_up, do_Tsources_up_p,                          & ! Inputs (Flags)
           do_columnwfs, do_surfacewfs, n_columnwfs, n_surfacewfs,                             & ! Inputs (Control, Jacobians)
           nstokes, nvzas, nlayers, n_user_levels, user_level_mask_up, npartials,              & ! Inputs (Control output)
           partial_outindex, partial_outflag, partial_layeridx, FOGeometry,                    & ! Inputs (Partial/Geometry)
           extinction, deltaus, omega, truncfac, bb_input,                                     & ! Inputs (Optical/Thermal)
           surfbb, emiss(1,:), emiss(2:4,:),                                                   & ! Inputs (Surface)
           L_extinction, L_deltaus, L_omega, L_truncfac,                                       & ! Inputs (Optical - Linearized)
           LS_emiss(1,:,:), LS_emiss(2:4,:,:),                                                 & ! Inputs (Surface - Linearized)
           intensity_dta_up_LOS, intensity_dts_LOS, LC_Jacobians_dta_up_LOS,                   & ! Main Outputs
           LC_Jacobians_dts_LOS, LS_Jacobians_dts_LOS, tcom1, L_tcom1,                         & ! Main Outputs
           lostrans_up_LOS, lostrans_up_p_LOS, L_lostrans_up_LOS, L_lostrans_up_p_LOS,         & ! Other Outputs
           StokesQUV_dts_LOS, LC_JacobiansQUV_dts_LOS, LS_JacobiansQUV_dts_LOS )                 ! Optional Output. 12/11/17 Rob Add.

!  Save results
!mick mod 9/19/2017 - turned off "fo_stokes", "fo_columnwf", & "fo_surfacewf" (defined later)
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

       if ( do_columnwfs ) then
         if ( do_ObsGeom ) then
           do g = 1, ngeoms
             do lev=1,n_user_levels
               do par=1,n_columnwfs
                 fo_columnwf_dta(par,lev,g,o1,upidx) = LC_Jacobians_dta_up_LOS(lev,g,par)
                 fo_columnwf_dts(par,lev,g,o1) = LC_Jacobians_dts_LOS(lev,g,par)
               enddo
             enddo
           enddo
         else if ( do_doublet ) then
           do nv = 1, nvzas ; do ns = 1, nszas
             g = nd_offset(ns) + nv
             do lev=1,n_user_levels
               do par=1,n_columnwfs
                 fo_columnwf_dta(par,lev,g,o1,upidx) = LC_Jacobians_dta_up_LOS(lev,nv,par)
                 fo_columnwf_dts(par,lev,g,o1)       = LC_Jacobians_dts_LOS(lev,nv,par)
               enddo
             enddo
           enddo ; enddo
         else
           do nv = 1, nvzas ; do ns = 1, nszas ; do na = 1, nazms
             g = na_offset(ns,nv) + na
             do lev=1,n_user_levels
               do par=1,n_columnwfs
                 fo_columnwf_dta(par,lev,g,o1,upidx) = LC_Jacobians_dta_up_LOS(lev,nv,par)
                 fo_columnwf_dts(par,lev,g,o1)       = LC_Jacobians_dts_LOS(lev,nv,par)
               enddo
             enddo
           enddo ; enddo ; enddo
         endif
       endif

       if ( do_surfacewfs ) then
         if ( do_ObsGeom ) then
           do g = 1, ngeoms
             do lev=1,n_user_levels
               do spar=1,n_surfacewfs
                 fo_surfacewf_dts(spar,lev,g,o1) = LS_Jacobians_dts_LOS(lev,g,spar)
               enddo
             enddo
           enddo
         else if ( do_doublet ) then
           do nv = 1, nvzas ; do ns = 1, nszas
             g = nd_offset(ns) + nv
!mick fix 9/19/2017 - added lev loop
             do lev=1,n_user_levels
               do spar=1,n_surfacewfs
                 fo_surfacewf_dts(spar,lev,g,o1) = LS_Jacobians_dts_LOS(lev,nv,spar)
               enddo
             enddo
           enddo ; enddo
         else
           do nv = 1, nvzas ; do ns = 1, nszas ; do na = 1, nazms
             g = na_offset(ns,nv) + na
!mick fix 9/19/2017 - added lev loop
             do lev=1,n_user_levels
               do spar=1,n_surfacewfs
                 fo_surfacewf_dts(spar,lev,g,o1) = LS_Jacobians_dts_LOS(lev,nv,spar)
               enddo
             enddo
           enddo ; enddo ; enddo
         endif
       endif

!  Save polarized Emissivity results. 12/11/17  Rob add.
!  -----------------------------------------------------

!  4/15/20. Version 1.5.2.  Add the Doublet goemetry option

       if ( do_Polarized_Emissivity.and.nstokes.gt.1 ) then

!  stokes vector

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

!  Column jacobians

         if ( do_columnwfs ) then
           if ( do_ObsGeom ) then
             do g = 1, ngeoms
               do lev=1,n_user_levels
                 do par=1,n_columnwfs
                   fo_columnwf_dts(par,lev,g,2:nstokes) = LC_JacobiansQUV_dts_LOS(lev,1:nstokes-1,g,par)
                 enddo
               enddo
             enddo
           else if ( do_doublet ) then
             do nv = 1, nvzas ; do ns = 1, nszas
               g = nd_offset(ns) + nv
               do lev=1,n_user_levels
                 do par=1,n_columnwfs
                   fo_columnwf_dts(par,lev,g,2:nstokes)       = LC_JacobiansQUV_dts_LOS(lev,1:nstokes-1,nv,par)
                 enddo
               enddo
             enddo ; enddo
           else
             do nv = 1, nvzas ; do ns = 1, nszas ; do na = 1, nazms
               g = na_offset(ns,nv) + na
               do lev=1,n_user_levels
                 do par=1,n_columnwfs
                   fo_columnwf_dts(par,lev,g,2:nstokes)       = LC_JacobiansQUV_dts_LOS(lev,1:nstokes-1,nv,par)
                 enddo
               enddo
             enddo ; enddo ; enddo
           endif
         endif

!  surface Jacobians

         if ( do_surfacewfs ) then
           if ( do_ObsGeom ) then
             do g = 1, ngeoms
               do lev=1,n_user_levels
                 do spar=1,n_surfacewfs
                   fo_surfacewf_dts(spar,lev,g,2:nstokes) = LS_JacobiansQUV_dts_LOS(lev,1:nstokes-1,g,spar)
                 enddo
               enddo
             enddo
           else if ( do_doublet ) then
             do nv = 1, nvzas ; do ns = 1, nszas
               g = nd_offset(ns) + nv
               do lev=1,n_user_levels
                 do spar=1,n_surfacewfs
                   fo_surfacewf_dts(spar,lev,g,2:nstokes) = LS_JacobiansQUV_dts_LOS(lev,1:nstokes-1,nv,spar)
                 enddo
               enddo
             enddo ; enddo
           else
             do nv = 1, nvzas ; do ns = 1, nszas ; do na = 1, nazms
               g = na_offset(ns,nv) + na
               do lev=1,n_user_levels
                 do spar=1,n_surfacewfs
                   fo_surfacewf_dts(spar,lev,g,2:nstokes) = LS_JacobiansQUV_dts_LOS(lev,1:nstokes-1,nv,spar)
                 enddo
               enddo
             enddo ; enddo ; enddo
           endif
         endif

!  End Polarized emissivity clause

       endif

!  End upwelling

     endif

!  Downwelling
!  -----------

     if ( do_dnwelling ) then

!  Sources, temporary fix until criticality realized. 9/17/16

       do_Tsources_dn   = .true. ; do_Tsources_dn_p = .true.

!  Direct thermal, calculate. Updated 9/17/16

        call FO_Thermal_DTRT_ILCS_DN &
         ( do_deltam_scaling, do_Partials, do_PlanPar, do_enhanced_ps,                 & ! Inputs (Flags)
           do_Tsources_dn, do_Tsources_dn_p, do_columnwfs, n_columnwfs,                & ! Inputs (Flags/Jac-control)
           nvzas, nlayers, n_user_levels, user_level_mask_dn,                          & ! Inputs (control output)
           npartials, partial_outindex, partial_outflag, partial_layeridx, FOGeometry, & ! Inputs (control-partial)
           bb_input, extinction, deltaus, omega, truncfac,                             & ! Inputs (Optical - Regular)
           L_extinction, L_deltaus, L_omega, L_truncfac,                               & ! Inputs (Optical - Linearized)
           intensity_dta_dn_LOS, LC_Jacobians_dta_dn_LOS, tcom1, L_tcom1 )               ! Output

!  Save results
!mick mod 9/19/2017 - turned off "fo_stokes" & "fo_columnwf" (defined later)
!  4/15/20. Version 1.5.2.  Add the Doublet goemetry option

       o1=1
       if ( do_obsgeom ) then
         do g = 1, nvzas
           do lev=1,n_user_levels
             fo_stokes_dta(lev,g,o1,dnidx) = intensity_dta_dn_LOS(lev,g)
           enddo
         enddo
       else if ( do_doublet ) then
         do nv = 1, nvzas ; do ns = 1, nszas
           g = nd_offset(ns) + nv
           do lev=1,n_user_levels
              fo_stokes_dta(lev,g,o1,dnidx) = intensity_dta_dn_LOS(lev,nv)
           enddo
         enddo ; enddo
       else
         do nv = 1, nvzas ; do ns = 1, nszas ; do na = 1, nazms
           g = na_offset(ns,nv) + na
           do lev=1,n_user_levels
             fo_stokes_dta(lev,g,o1,dnidx) = intensity_dta_dn_LOS(lev,nv)
           enddo
         enddo ; enddo ; enddo
       endif

       if ( do_columnwfs ) then
          if (do_ObsGeom ) then
           do g = 1, ngeoms
             do lev=1,n_user_levels
               do par=1,n_columnwfs
                 fo_columnwf_dta(par,lev,g,o1,dnidx) = LC_Jacobians_dta_dn_LOS(lev,g,par)
               enddo
             enddo
           enddo
         else if ( do_doublet ) then
           do nv = 1, nvzas ; do ns = 1, nszas
             g = nd_offset(ns) + nv
             do lev=1,n_user_levels
               do par=1,n_columnwfs
                 fo_columnwf_dta(par,lev,g,o1,dnidx) = LC_Jacobians_dta_dn_LOS(lev,nv,par)
               enddo
             enddo
           enddo ; enddo
         else
           do nv = 1, nvzas ; do ns = 1, nszas ; do na = 1, nazms
             g = na_offset(ns,nv) + na
             do lev=1,n_user_levels
               do par=1,n_columnwfs
                 fo_columnwf_dta(par,lev,g,o1,dnidx) = LC_Jacobians_dta_dn_LOS(lev,nv,par)
               enddo
             enddo
           enddo ; enddo ; enddo
         endif
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
     if ( do_columnwfs ) then
       do o1 = 1, nstokes
         do g = 1, ngeoms
           do lev=1,n_user_levels
             do par=1,n_columnwfs
               fo_columnwf_atmos(par,lev,g,o1,upidx) = fo_columnwf_ss(par,lev,g,o1,upidx) &
                                                     + fo_columnwf_dta(par,lev,g,o1,upidx)
               fo_columnwf_surf(par,lev,g,o1)        = fo_columnwf_db(par,lev,g,o1) &
                                                     + fo_columnwf_dts(par,lev,g,o1)
               fo_columnwf(par,lev,g,o1,upidx)       = fo_columnwf_atmos(par,lev,g,o1,upidx) &
                                                     + fo_columnwf_surf(par,lev,g,o1)
             enddo
           enddo
         enddo
       enddo
     endif
     if ( do_surfacewfs ) then
       do o1 = 1, nstokes
         do g = 1, ngeoms
           do lev=1,n_user_levels
             do spar=1,n_surfacewfs
               fo_surfacewf(spar,lev,g,o1) = fo_surfacewf_db(spar,lev,g,o1) + fo_surfacewf_dts(spar,lev,g,o1)
             enddo
           enddo
         enddo
       enddo
     endif
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
     if ( do_columnwfs ) then
       do o1 = 1, nstokes
         do g = 1, ngeoms
           do lev=1,n_user_levels
             do par=1,n_columnwfs
               fo_columnwf_atmos(par,lev,g,o1,dnidx) = fo_columnwf_ss(par,lev,g,o1,dnidx) &
                                                     + fo_columnwf_dta(par,lev,g,o1,dnidx)
               fo_columnwf(par,lev,g,o1,dnidx)       = fo_columnwf_atmos(par,lev,g,o1,dnidx)
             enddo
           enddo
         enddo
       enddo
     endif
   endif

!  Finish

   return
end subroutine VFO_RTCALC_LCS_MASTER

!  End module

end module VFO_RTCalc_LinMasters_m

