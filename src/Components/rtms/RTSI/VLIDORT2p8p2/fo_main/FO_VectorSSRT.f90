
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

module FO_Vector_SSRT_m

!  FUNCTION
!  ========

!  For given wavelength, Module calculates First-Order upwelling+downwelling Stokes vectors

!  1. For the Atmospheric Solar Single-scatter and Surface Direct-Beam (SS) sources.
!  2. This is based on Precalculated Geometrical quantities and appropriate Optical properties.
!  3. This will perform Enhanced-PS calculations (incoming solar and outgoing LOS-path sphericity) 
!  4. This will perform Regular-PS  calculations (plane-parallel or incoming solar pseudo-spherical)

!  THIS VERSION
!  ============

!  4/15/20. Version 1.5.2, substantial rewrite for VLIDORT 2.8.2
!   - Geometry type structure introduced as input. No longer stand-alone

!  HISTORICAL NOTES, up to Version 1.5.1
!  =====================================

!  Versions to 1.4, without Partials. Code is stand alone with no dependencies.
!  Version     1.5, with optional phase function, surface leaving and partials.

!    Version 1.1a, 01 December 2011, R. Spurr, RT Solutions Inc.
!    Version 1.1b, 13 February 2012, R. Spurr, RT Solutions Inc.
!    Version 1.2,  01 June     2012, R. Spurr, RT Solutions Inc.
!    Version 1.3,  29 October  2012, Extension to Observational multiple geometries
!    Version 1.4,  31 July     2013, Extension to Lattice       multiple geometries
!    Version 1.5,  07 July     2016, Optional calculation using Phase functions directly
!    Version 1.5,  02 August   2016. Inclusion of Surface-leaving terms + LSSL weighting functions
!    Version 1.5,  17 August   2016, Partial-layer output introduced
!    Version 1.5.1 09 April    2019, Add the CUMTRANS output, add water-leaving control

!  SUBROUTINES
!  ===========

!  For Solar sources, the subroutines are
!       FO_Vector_SSRT_UP   (Upwelling only)
!       FO_Vector_SSRT_DN   (Downwelling only)
!       FO_Vector_SSRT_UPDN (Upwelling and Downwelling)

!  Dependencies
!  ============

   use VLIDORT_PARS_m          , only : zero, one, Expcutoff, MAX_USER_LEVELS, MAXLAYERS, &
                                        MAXSTOKES, MAX_PARTLAYERS, MAXFINELAYERS, MAX_GEOMETRIES
   use VLIDORT_Setups_def_m

!  All three subroutines public

public

contains

subroutine FO_Vector_SSRT_UP &
   ( do_sunlight, do_deltam_scaling, do_Lambertian, do_surface_leaving, do_water_leaving, & ! Inputs (Flags-General/Surface)
     do_Partials, do_PlanPar, do_enhanced_ps, do_sources_up, do_sources_up_p,             & ! Inputs(Flags/criticality)
     nstokes, ngeoms, nlayers, n_user_levels, user_levels,                                & ! Inputs (control/flux)
     npartials, partial_outindex, partial_outflag, partial_layeridx, FOGeometry,          & ! Inputs (control-partial)
     flux, fluxvec, extinction, deltaus, omega, truncfac, fmatrix_up, Reflec, Slterm,     & ! Inputs (Geometry/Optical/surface)
     stokes_up, stokes_db, cumsource_up, cumtrans )                                         ! Outputs

!  FO Routine for calculation of Upwelling Solar-beam Single-scatter (SS) radiation field.
!    Computation of Stokes-vector. Inputs: Control, geometry, optical properties.

   implicit none         

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ======

!  flags.
!  Version 1.5: --> F-matrix flag added 7/7/16; surface-leaving flag 8/2/16; Partials 8/20/16
!  4/15/20. Version 2.8.2. DO_FMATRIX flag removed; Now we must use FMATRIX input (no choice)

   LOGICAL, Intent(in) :: DO_SUNLIGHT
   LOGICAL, Intent(in) :: DO_DELTAM_SCALING

   LOGICAL, Intent(in) :: DO_LAMBERTIAN
   LOGICAL, Intent(in) :: DO_SURFACE_LEAVING
   logical, Intent(in) :: DO_WATER_LEAVING

   logical, Intent(in) :: DO_Partials
   LOGICAL, Intent(in) :: DO_PLANPAR
   LOGICAL, Intent(in) :: DO_ENHANCED_PS

!  Existence flags. 8/19/16. Criticality enters here

   logical, Intent(in)    :: do_sources_up       (maxlayers,MAX_GEOMETRIES)
   logical, Intent(in)    :: do_sources_up_p     (MAX_PARTLAYERS,MAX_GEOMETRIES)

!  Numbers
!  4/15/20. Version 2.8.2. NMOMSINP removed; Now we must use FMATRIX input (no choice)

   INTEGER, Intent(in) :: NSTOKES, NGEOMS, NLAYERS, N_USER_LEVELS
   INTEGER, Intent(in) :: USER_LEVELS ( MAX_USER_LEVELS )

!  Solar Flux 

   real(ffp), Intent(in) :: FLUX, fluxvec(4)

!  Numbers for Version 1.5: -->  Partial Control added, 8/20/16

   integer, Intent(in) :: Npartials
   integer, Intent(in) :: partial_layeridx(MAX_PARTLAYERS)
   logical, Intent(in) :: partial_outflag ( MAX_USER_LEVELS )
   integer, Intent(in) :: partial_outindex( MAX_USER_LEVELS )

!  Geometrical inputs
!  ==================

   Type(VLIDORT_Geometry_FO), Intent(in) :: FOGeometry

!  optical inputs
!  --------------

!  Atmosphere. Fmatrix input added 7/7/16
!  4/15/20. Version 2.8.2. GREEKMAT removed; Now we must use FMATRIX input (no choice)

   REAL(ffp), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   REAL(ffp), Intent(in) :: DELTAUS     ( MAXLAYERS )
   REAL(ffp), Intent(in) :: OMEGA       ( MAXLAYERS )
   REAL(ffp), Intent(in) :: TRUNCFAC    ( MAXLAYERS )
   REAL(ffp), Intent(in) :: FMATRIX_UP  ( MAXLAYERS, MAX_GEOMETRIES, 6 )

!  Surface reflectivity (Could be the albedo)
!    Surface leaving input added 8/2/16

   real(ffp), Intent(in) :: REFLEC ( 4, 4, MAX_GEOMETRIES )
   real(ffp), Intent(in) :: SLTERM ( 4,    MAX_GEOMETRIES )

!  outputs
!  -------

   REAL(ffp), Intent(Out)  :: stokes_up     ( max_user_levels, 4, MAX_GEOMETRIES )
   REAL(ffp), Intent(Out)  :: stokes_db     ( max_user_levels, 4, MAX_GEOMETRIES )
   REAL(ffp), Intent(Out)  :: cumsource_up  ( 0:maxlayers,     4, MAX_GEOMETRIES )

!  4/9/19. Additional output for the sleave correction

   real(ffp), Intent(out)  :: CUMTRANS ( max_user_levels, MAX_GEOMETRIES )

!  LOCAL
!  -----

!  Attenuations. Partials added, 8/20/16

   real(ffp)  :: attenuations       (0:maxlayers)
   real(ffp)  :: attenuationsfine   (maxlayers,MAXFINELAYERS)
   real(ffp)  :: attenuations_p     (MAX_PARTLAYERS)
   real(ffp)  :: attenuationsfine_p (MAX_PARTLAYERS,MAXFINELAYERS)

!  Scattering

   REAL(ffp)  :: tms (maxlayers)
   REAL(ffp)  :: exactscat_up (maxlayers,4,4)

!  Source function integration results

   real(ffp)  :: sources_up  (maxlayers, 4), sources_up_p  (MAX_PARTLAYERS, 4)
   real(ffp)  :: lostrans_up (maxlayers), lostrans_up_p (MAX_PARTLAYERS)

   real(ffp)  :: multiplier       ( maxlayers )
   real(ffp)  :: multiplier_p     ( MAX_PARTLAYERS )

!  Regular_PS or plane-parallel flag

   LOGICAL    :: do_RegPSorPP

!  Help

   LOGICAL    :: do_regular_ps, layermask_up(maxlayers)
   INTEGER    :: n, ns, uta, nstart, nc, nut, nut_prev, j, o1, v, nt, np, ut

   REAL(ffp)  :: sumd, help, sum, tran, factor1, factor2, kn, pi4, dj, path_up, Shelp(4)
   REAL(ffp)  :: cumsource_db(4), ctrans, CUMSOURCE_DB_START(4)
   REAL(ffp)  :: help3c1, help3s1, help4c1, help4s1, m4
   real(ffp)  :: suntau(0:maxlayers), suntau_p(MAX_PARTLAYERS), lostau

!  Number

!mick fix 9/19/2017 - define pi4 as in VLIDORT_PARS
   !pi4 = acos(-one)/4.0_ffp
   pi4 = acos(-one)*4.0_ffp

!  Zero the output. 4/9/19 include CUMTRANS

   CUMSOURCE_UP  = zero ; STOKES_UP  = zero ; STOKES_DB    = zero ; cumtrans = zero

!  Regular_PS or plane-parallel flag

   do_regular_ps = .false.
   if ( .not.do_Planpar ) do_regular_ps = .not. do_enhanced_ps
   do_RegPSorPP = (do_regular_ps .or. do_PlanPar)

!  Bookkeeping

   ns = nstokes
   NUT = USER_LEVELS(1) + 1
   LAYERMASK_UP = .false.
   LAYERMASK_UP(NUT:NLAYERS) = .true.

!  TMS factors

   do n = 1, nlayers
      if ( do_deltam_scaling ) then
         help = one - truncfac(n) * omega(n)
         tms(n) = omega(n) / help
      else
         tms(n) = omega(n)
      endif
   enddo

!  Start geometry loop
!  -------------------

   do v = 1, ngeoms

!  Zero the local sources

      lostrans_up   = zero  ; sources_up   = zero ; exactscat_up = zero
      lostrans_up_p = zero  ; sources_up_p = zero

!  Scattering functions
!  ====================

!  Version 1.5, F-matrix option introduced. 7/7/16
!mick note 9/19/2017 - defining F-matrix using either (1) F-matrix input ("do_fmatrix")
!                                                  or (2) Greek matrix moment input

!  4/15/20. Version 1.5.2 for VLIDORT_2.8.2, Use F-matrices only.

!  Scalar only
!  -----------

      if ( nstokes .eq. 1 ) then
        do n = 1, nlayers
          if ( layermask_up(n) ) then
            exactscat_up(n,1,1) = fmatrix_up(n,v,1) * tms(n)
          endif
        enddo
      endif

!  Vector with Sunlight
!  --------------------

!   Apply rotations for Z-matrix, multiply by TMS

      if ( nstokes .gt. 1 .and. do_sunlight ) then
        do n = 1, nlayers
          if ( layermask_up(n) ) then
            exactscat_up(n,1,1) = + fmatrix_up(n,v,1)
            exactscat_up(n,2,1) = - fmatrix_up(n,v,2) * FOGeometry%Rotations_up(3,v)
            exactscat_up(n,3,1) = + fmatrix_up(n,v,2) * FOGeometry%Rotations_up(4,v)
            exactscat_up(n,1:ns,1) = tms(n) * exactscat_up(n,1:ns,1) 
          endif
        enddo
      endif

!  Vector General case
!  -------------------

! USE FULL 4X4 MATRIX; CODE INTRODUCED BUT NOT TESTED, 05 OCTOBER 2010

!   Apply rotations for Z-matrix, multiply by TMS

      if ( nstokes .gt. 1 .and. .not. do_sunlight ) then
        do n = 1, nlayers
          if ( layermask_up(n) ) then
            help3c1 = fmatrix_up(n,v,3) * FOGeometry%Rotations_up(1,v)
            help3s1 = fmatrix_up(n,v,3) * FOGeometry%Rotations_up(2,v)
            help4c1 = fmatrix_up(n,v,4) * FOGeometry%Rotations_up(1,v)
            help4s1 = fmatrix_up(n,v,4) * FOGeometry%Rotations_up(2,v)
            exactscat_up(n,1,1) = + fmatrix_up(n,v,1)
            exactscat_up(n,2,1) = - fmatrix_up(n,v,2) * FOGeometry%Rotations_up(3,v)
            exactscat_up(n,3,1) = + fmatrix_up(n,v,2) * FOGeometry%Rotations_up(4,v)
            exactscat_up(n,1,2) = + fmatrix_up(n,v,2) * FOGeometry%Rotations_up(1,v)
            exactscat_up(n,1,3) = - fmatrix_up(n,v,2) * FOGeometry%Rotations_up(2,v)
            exactscat_up(n,2,2) = + help3c1 * FOGeometry%Rotations_up(3,v) - help4s1 * FOGeometry%Rotations_up(4,v)
            exactscat_up(n,2,3) = - help3s1 * FOGeometry%Rotations_up(3,v) - help4c1 * FOGeometry%Rotations_up(4,v)
            exactscat_up(n,3,2) = + help3c1 * FOGeometry%Rotations_up(4,v) + help4s1 * FOGeometry%Rotations_up(3,v)
            exactscat_up(n,3,3) = - help3s1 * FOGeometry%Rotations_up(4,v) + help4c1 * FOGeometry%Rotations_up(3,v)
            if ( nstokes .eq. 4 ) then
               exactscat_up(n,2,4) = - fmatrix_up(n,v,5) * FOGeometry%Rotations_up(4,v) 
               exactscat_up(n,4,2) = - fmatrix_up(n,v,5) * FOGeometry%Rotations_up(2,v) 
               exactscat_up(n,3,4) = + fmatrix_up(n,v,5) * FOGeometry%Rotations_up(3,v) 
               exactscat_up(n,4,3) = - fmatrix_up(n,v,5) * FOGeometry%Rotations_up(1,v) 
               exactscat_up(n,4,4) = + fmatrix_up(n,v,6)
            endif
            exactscat_up(n,1:ns,1:ns) = tms(n)*exactscat_up(n,1:ns,1:ns)
          endif
        enddo
      endif

!  Attenuations and Solar solutions
!  ================================

!  Initialize, only to layer Ncrit if applicable

      Attenuations   = zero ; suntau = zero   ; Attenuationsfine   = zero
      Attenuations_p = zero ; suntau_p = zero ; Attenuationsfine_p = zero
      nstart = nlayers

!  Critical removed TRC
!  Initialize, only to layer Ncrit if applicable
!     if (Ncrit(v).ne.0) nstart = nCrit(v)

!  Attenuations to End points (including TOA). Both PS representations
!    MUST go all the way to NLAYERS (surface term required)

      do n = 0, nlayers
         nt = FOGeometry%ntraverse_up(n,v) ; sumd = dot_product(extinction(1:nt),FOGeometry%sunpaths_up(n,1:nt,v))
         suntau(n) = sumd    ; If (sumd .lt. Expcutoff ) Attenuations(n) = exp( - sumd )
      enddo

!  RobFix 8/20/16. Attenuations to partial-layer points

      if ( do_Partials ) then
        do ut = 1, npartials
          nt = FOGeometry%ntraverse_p_up(ut,v) ; sumd = dot_product(extinction(1:nt),FOGeometry%sunpaths_p_up(ut,1:nt,v))
          suntau_p(ut) = sumd    ; If (sumd .lt. Expcutoff ) Attenuations_p(ut) = exp( - sumd )
        enddo
      endif

!  Enhanced-spherical, fine-layer attenuations, Whole-layer integration

      if ( do_enhanced_ps ) then
        do n = 1, nlayers
          if ( layermask_up(n) .and. do_sources_up(n,v) ) then
            do j = 1, FOGeometry%nfinedivs(n,v)
              nt = FOGeometry%ntraversefine_up(n,j,v)
              sumd = dot_product(extinction(1:nt),FOGeometry%sunpathsfine_up(n,1:nt,j,v))
              if (sumd .lt. Expcutoff ) Attenuationsfine(n,j) = exp( - sumd )
            enddo
          endif
        enddo
      endif

!  RobFix 8/20/16. Enhanced-spherical, fine-layer attenuations, Partial layer integration

      if ( do_enhanced_ps .and. do_Partials ) then
        do ut = 1, npartials
          if ( do_sources_up_p(ut,v) ) then
            np = partial_layeridx(ut)
            do j = 1, FOGeometry%nfinedivs_p_up(ut,v)
              nt = FOGeometry%ntraversefine_p_up(ut,j,v)
              sumd = dot_product(extinction(1:nt),FOGeometry%sunpathsfine_p_up(ut,1:nt,j,v))
              If (sumd .lt. Expcutoff ) Attenuationsfine_p(ut,j) = exp( - sumd )
            enddo
          endif
        enddo
      endif

!  Layer integrated Solar sources
!  ==============================

!  Plane/Parallel or Regular-PS, Whole-layer source terms
!  ------------------------------------------------------

      if ( do_RegPSorPP ) then
        do n = nlayers, 1, -1
          factor1 = zero ; factor2 = zero
          if ( layermask_up(n) .and. do_sources_up(n,v)  ) then
            if ( FOGeometry%Mu1_up(v) .gt. zero ) then
              lostau = deltaus(n)  / FOGeometry%Mu1_up(v)
              if ( lostau .lt. Expcutoff ) lostrans_up(n) = exp( - lostau )
              factor1 = Attenuations(n-1) - Attenuations(n)*lostrans_up(n)
              factor2 = one + (suntau(n) - suntau(n-1))/lostau
              multiplier(n) = factor1 / factor2
            endif
          endif
        enddo
      endif

!  RobFix 8/20/16. Plane/Parallel or Regular-PS, Partial-layer output

      if ( do_RegPSorPP .and.do_Partials ) then
        do ut = 1, npartials
          if ( do_sources_up_p(ut,v) ) then
            np = Partial_layeridx(ut) ; kn = extinction(np)
            path_up = FOGeometry%LosW_paths(np,v) - FOGeometry%LosP_paths(ut,v)
            factor1 = zero ; factor2 = zero
            if ( FOGeometry%Mu1_up(v) .gt. zero ) then
              lostau = kn * path_up
              if ( lostau .lt. Expcutoff ) lostrans_up_p(ut) = exp( - lostau )
              factor1 = Attenuations_p(ut) - Attenuations(np)*lostrans_up_p(ut)
              factor2 = one + (suntau(np) - suntau_p(ut))/lostau
              multiplier_p(ut) = factor1 / factor2
            endif
          endif
        enddo
      endif

!  Enhanced PS: General case, whole layers. 
!  ----------------------------------------

!     RobFix 8/20/16 streamlined code using distances
!      Quadratures from Bottom of the layer

      if ( do_enhanced_ps ) then
        do n = nlayers, 1, -1
          if ( layermask_up(n) .and. do_sources_up(n,v)   ) then
!mick fix 3/22/2017 - replaced index "np" with "n" in "LosW_paths"
            kn = extinction(n) ; path_up = FOGeometry%LosW_paths(n,v)
            lostau = kn * path_up ; if( lostau.lt.Expcutoff ) lostrans_up(n) = exp ( - lostau )
            sum = zero
            do j = 1, FOGeometry%nfinedivs(n,v)
              dj = FOGeometry%LosW_paths(n,v) - FOGeometry%xfine(n,j,v) ; tran = exp ( - kn * dj )
              sum  = sum + Attenuationsfine(n,j) * tran * FOGeometry%wfine(n,j,v)
            enddo
            multiplier(n) = sum * kn
          endif
        enddo
      endif

!  Enhanced PS: General case, partial layers. 
!  -----------------------------------------

!     RobFix 8/02/16 streamlined code using distances
!      Quadratures from Bottom of the layer

      if ( do_enhanced_ps .and. do_partials ) then
        do ut = 1, npartials
          if ( do_sources_up_p(ut,v) ) then
            np = partial_layeridx(ut) ; kn = extinction(np)
            path_up = FOGeometry%LosW_paths(np,v)- FOGeometry%LosP_paths(ut,v)
            lostau = kn * path_up ; if ( lostau.lt.Expcutoff ) lostrans_up_p(ut) = exp ( - lostau )
            sum = zero
            do j = 1, FOGeometry%nfinedivs_p_up(ut,v)
              dj = path_up - FOGeometry%xfine_p_up(ut,j,v) ; tran = exp ( - kn * dj )     ! Correct
              sum  = sum + Attenuationsfine_p(ut,j) * tran * FOGeometry%wfine_p_up(ut,j,v)
            enddo
            multiplier_p(ut) = sum * kn
          endif
        enddo
      endif

!  Layer integrated Solar sources
!  ==============================

!  General case, Whole layers
!mick fix 9/19/2017 - defined "sources_up" for all stokes vector elements

      do n = nlayers, 1, -1
        if ( layermask_up(n) .and. do_sources_up(n,v)  ) then
          if ( do_sunlight ) then
             shelp(1:ns) = exactscat_up(n,1:ns,1) * fluxvec(1)
          else
            do o1 = 1, nstokes
              shelp(o1) = dot_product(exactscat_up(n,o1,1:ns),fluxvec(1:ns))
            enddo
          endif
          !sources_up(n,o1) = shelp(o1) * multiplier(n)
          sources_up(n,1:ns) = shelp(1:ns) * multiplier(n)
        endif
      enddo

!  Partials case

      if ( do_partials ) then
        do ut = 1, npartials
          if ( do_sources_up_p(ut,v)  ) then
            np = partial_layeridx(ut)
            if ( do_sunlight ) then
              shelp(1:ns) = exactscat_up(np,1:ns,1) * fluxvec(1)
            else
              do o1 = 1, nstokes
                shelp(o1) = dot_product(exactscat_up(np,o1,1:ns),fluxvec(1:ns))
              enddo
            endif
            sources_up_p(ut,1:ns) = shelp(1:ns)* multiplier_p(ut)
          endif
        enddo
      endif

!  Source function integration
!  ===========================

!  initialize recursion ( For Direct Beam, use PI.mu0.R.Atten )
!  4/9/19. Add start of CUMTRANS recursion (CTRANS = 1.0). Rename CUMSOURCE_DB_START

      NC =  0
      CUMSOURCE_UP(NC,:,V) = zero
      CUMSOURCE_DB_START = zero
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  Surface term

      M4 = 4.0d0 * FOGeometry%Mu0_up(V)
      if ( DO_LAMBERTIAN ) then
         CUMSOURCE_DB_START(1) = M4 * REFLEC(1,1,V) * attenuations(nlayers) * fluxvec(1)
      else
         do o1 = 1, nstokes
            sum = dot_product(REFLEC(O1,1:nstokes,V),fluxvec(1:nstokes))
            CUMSOURCE_DB_START(o1) = M4 * sum * attenuations(nlayers)      ! Bug 3/23/15 attenuations was left out
         enddo
      endif

!  surface-leaving term. Added, 8/2/16
!   -- (modeled after the DBCORRECTION code in Version 2.7)
!   -- 4/9/19. Not done for water-leaving, as need to use adjusted values

     IF ( DO_SURFACE_LEAVING .and. .not. DO_WATER_LEAVING ) THEN
        do o1 = 1, nstokes
           CUMSOURCE_DB_START(o1) = CUMSOURCE_DB_START(o1) + PI4 * SLTERM(o1,v)
        enddo
     ENDIF

!  Main loop over all output optical depths
!     NLEVEL = Layer index for given optical depth
!     Cumulative source terms : Loop over layers working upwards from NSTART to level NUT,
!     Check for updating the recursion

      DO UTA = N_USER_LEVELS, 1, -1
         NUT    = USER_LEVELS(UTA) + 1
         DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            do o1 = 1, nstokes
               CUMSOURCE_DB_START(O1) = LOSTRANS_UP(N) * CUMSOURCE_DB_START(O1)
               CUMSOURCE_UP(NC,O1,V)  = LOSTRANS_UP(N) * CUMSOURCE_UP(NC-1,O1,V) + SOURCES_UP(N,O1)
            enddo
         ENDDO
         CUMSOURCE_DB(1:nstokes)     = CUMSOURCE_DB_START(1:nstokes)
         IF ( Partial_OUTFLAG(UTA) ) THEN
            UT = Partial_OUTINDEX(UTA)
            do o1 = 1, nstokes
              STOKES_UP(UTA,O1,V) = FLUX * ( CUMSOURCE_UP(NC,O1,V) * LOSTRANS_UP_p(UT) + SOURCES_UP_p(UT,O1) )
              STOKES_DB(UTA,O1,V) = FLUX * CUMSOURCE_DB(O1) * LOSTRANS_UP_p(UT)
            enddo
         ELSE
           do o1 = 1, nstokes
             STOKES_UP(UTA,O1,V) = FLUX * CUMSOURCE_UP(NC,O1,V)
             STOKES_DB(UTA,O1,V) = FLUX * CUMSOURCE_DB(O1)
           enddo
         ENDIF
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
         NUT_PREV = NUT
      ENDDO

!  4/9/19.  WATERLEAVING CASE. Add CUMTRANS calculation
!    Add start of CUMTRANS recursion (CTRANS = 1.0).
      
!  5/22/20. Version 2.8.2 Upgrades. CUMTRANS calculation was not properly initialized

      if ( do_water_leaving ) then
         NSTART = NLAYERS                ! 5/22/20 Needed for initializing
         NUT_PREV = NSTART + 1           ! 5/22/20 Needed for initializing
         ctrans = one
         DO UTA = N_USER_LEVELS, 1, -1
            NUT    = USER_LEVELS(UTA) + 1
            DO N = NSTART, NUT, -1
               CTRANS = CTRANS * LOSTRANS_UP(N)
            ENDDO
            IF ( Partial_OUTFLAG(UTA) ) THEN
               UT = Partial_OUTINDEX(UTA)
               CUMTRANS(UTA,V)     = CTRANS * LOSTRANS_UP_p(UT)
            ELSE
               CUMTRANS(UTA,V)     = CTRANS
            ENDIF
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
            NUT_PREV = NUT
         ENDDO
      endif
      
!  Finish geometry loop

   enddo

!  Finish

   return
end subroutine FO_Vector_SSRT_UP

!

subroutine FO_Vector_SSRT_DN &
   ( do_sunlight, do_deltam_scaling, do_Partials,                                & ! Inputs (Flags)
     do_PlanPar, do_enhanced_ps, do_sources_dn, do_sources_dn_p,                 & ! Inputs (Flags)
     nstokes, ngeoms, nlayers, n_user_levels, user_levels,                       & ! Inputs (control/flux)
     npartials, partial_outindex, partial_outflag, partial_layeridx, FOGeometry, & ! Inputs (control-partial)
     flux, fluxvec, extinction, deltaus, omega, truncfac, fmatrix_dn,            & ! Inputs (Optical)
     stokes_dn, cumsource_dn )                                                     ! Outputs

!  FO Routine for calculation of Downwelling Solar-beam Single-scatter (SS) radiation field.
!    computation of Stokes vector. Inputs: geometry, spherical functions, optical properties.

   implicit none         

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ======


!  flags. F-matrix flag added, 7/7/16
!  4/15/20. Version 2.8.2. DO_FMATRIX flag removed; Now we must use FMATRIX (no choice)

   LOGICAL, Intent(in) ::  DO_SUNLIGHT
   LOGICAL, Intent(in) ::  DO_DELTAM_SCALING

   logical, Intent(in) ::  DO_Partials
   LOGICAL, Intent(in) ::  DO_PLANPAR
   LOGICAL, Intent(in) ::  DO_ENHANCED_PS

!  Existence flags. 8/20/16. Criticality enters here

   logical, Intent(in)    :: do_sources_dn       (maxlayers,MAX_GEOMETRIES)
   logical, Intent(in)    :: do_sources_dn_p     (MAX_PARTLAYERS,MAX_GEOMETRIES)

!  Numbers

   INTEGER, Intent(in) ::  NSTOKES, NGEOMS, NLAYERS, N_USER_LEVELS
   INTEGER, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  Numbers for Version 1.5: -->  Partial Control added, 8/20/16

   integer, Intent(in) :: Npartials
   integer, Intent(in) :: partial_layeridx(MAX_PARTLAYERS)
   logical, Intent(in) :: partial_outflag ( MAX_USER_LEVELS )
   integer, Intent(in) :: partial_outindex( MAX_USER_LEVELS )

!  Geometrical inputs
!  ==================

   Type(VLIDORT_Geometry_FO), Intent(in) :: FOGeometry

!  optical inputs
!  --------------

!  Solar Flux 

   REAL(ffp), Intent(in) :: FLUX, FLUXVEC(4)

!  Atmosphere. Fmatrix input added 7/7/16
!  4/15/20. Version 2.8.2. GREEKMAT removed; Now we must use FMATRIX_DN

   REAL(ffp), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   REAL(ffp), Intent(in) :: DELTAUS     ( MAXLAYERS )
   REAL(ffp), Intent(in) :: OMEGA       ( MAXLAYERS )
   REAL(ffp), Intent(in) :: TRUNCFAC    ( MAXLAYERS )
   REAL(ffp), Intent(in) :: FMATRIX_DN  ( MAXLAYERS, MAX_GEOMETRIES, 6 )

!  outputs
!  -------

   REAL(ffp), Intent(Out)  :: stokes_dn     ( max_user_levels, 4, MAX_GEOMETRIES )
   REAL(ffp), Intent(Out)  :: cumsource_dn  ( 0:maxlayers,     4, MAX_GEOMETRIES )

!  LOCAL
!  -----

!  Attenuations. Partials added, 8/20/16

   real(ffp)  :: attenuations       (0:maxlayers)
   real(ffp)  :: attenuationsfine   (maxlayers,MAXFINELAYERS)
   real(ffp)  :: attenuations_p     (MAX_PARTLAYERS)
   real(ffp)  :: attenuationsfine_p (MAX_PARTLAYERS,MAXFINELAYERS)

!  Scattering

   REAL(ffp)  :: tms (maxlayers)
   REAL(ffp)  :: exactscat_dn (maxlayers,4,4)

!  Source function integration results

   real(ffp)  :: sources_dn  (maxlayers, 4), sources_dn_p  (MAX_PARTLAYERS, 4)
   real(ffp)  :: lostrans_dn (maxlayers), lostrans_dn_p (MAX_PARTLAYERS)

   real(ffp)  :: multiplier       ( maxlayers )
   real(ffp)  :: multiplier_p     ( MAX_PARTLAYERS )

!  Regular_PS or plane-parallel flag

   LOGICAL    :: do_RegPSorPP

!  Help

   INTEGER    :: n, ns, uta, nstart, nc, nut, nut_prev, j, O1, v, nt, ut, np
   logical    :: do_regular_ps, layermask_dn(maxlayers)

   REAL(ffp)  :: sumd, help, sum, tran, factor1, factor2, kn, dj
   real(ffp)  :: suntau(0:maxlayers), suntau_p(MAX_PARTLAYERS), lostau, path_dn, Shelp(4)
   REAL(ffp)  :: help3c1, help3s1, help4c1, help4s1

!  Zero the output

   CUMSOURCE_DN = zero ; STOKES_DN  = zero

!  Regular_PS or plane-parallel flag

   do_regular_ps = .false.
   if ( .not.do_Planpar ) do_regular_ps = .not. do_enhanced_ps
   do_RegPSorPP = (do_regular_ps .or. do_PlanPar)

!  Bookkeeping

   ns = nstokes 
   NUT = USER_LEVELS(N_USER_LEVELS) + 1
   IF ( NUT > NLAYERS ) NUT = NLAYERS
   LAYERMASK_DN = .false.
   LAYERMASK_DN(1:NUT) = .true.

!  TMS factors

   do n = 1, nlayers
      if ( do_deltam_scaling ) then
         help = one - truncfac(n) * omega(n)
         tms(n) = omega(n) / help
      else
         tms(n) = omega(n)
      endif
   enddo

!  Start geometry loop
!  ===================

   do v = 1, ngeoms

!  Zero the local sources

      lostrans_dn   = zero  ; sources_dn   = zero ; exactscat_dn = zero
      lostrans_dn_p = zero  ; sources_dn_p = zero

!  Scattering functions
!  ====================

!  Version 1.5, F-matrix option introduced. 7/7/16

!  Scalar only
!  -----------

!  4/15/20. Version 1.5.2 for LIDORT_2.8.2, Use Phase function only.

      if ( nstokes .eq. 1 ) then
        do n = 1, nlayers
          if ( layermask_dn(n) ) then
            exactscat_dn(n,1,1) = fmatrix_dn(n,v,1) * tms(n)
          endif
        enddo
      endif

!  Vector with Sunlight
!  --------------------

!   Apply rotations for Z-matrix, multiply by TMS
!  4/15/20. Version 1.5.2 for LIDORT_2.8.2, Use F-matrix only.

      if ( nstokes .gt. 1 .and. do_sunlight ) then
        do n = 1, nlayers
          if ( layermask_dn(n) ) then
            exactscat_dn(n,1,1) = + fmatrix_dn(n,v,1)
            exactscat_dn(n,2,1) = - fmatrix_dn(n,v,2) * FOGeometry%Rotations_dn(3,v)
            exactscat_dn(n,3,1) = + fmatrix_dn(n,v,2) * FOGeometry%Rotations_dn(4,v)
            exactscat_dn(n,1:ns,1) = tms(n) * exactscat_dn(n,1:ns,1) 
          endif
        enddo
      endif

!  Vector General case
!  -------------------

! USE FULL 4X4 MATRIX; CODE INTRODUCED BUT NOT TESTED, 05 OCTOBER 2010
!   Apply rotations for Z-matrix, multiply by TMS
!  4/15/20. Version 1.5.2 for LIDORT_2.8.2, Use F-matrix only.

      if ( nstokes .gt. 1 .and. .not. do_sunlight ) then
        do n = 1, nlayers
          if ( layermask_dn(n) ) then
            help3c1 = fmatrix_dn(n,v,3) * FOGeometry%Rotations_dn(1,v)
            help3s1 = fmatrix_dn(n,v,3) * FOGeometry%Rotations_dn(2,v)
            help4c1 = fmatrix_dn(n,v,4) * FOGeometry%Rotations_dn(1,v)
            help4s1 = fmatrix_dn(n,v,4) * FOGeometry%Rotations_dn(2,v)
            exactscat_dn(n,1,1) = + fmatrix_dn(n,v,1)
            exactscat_dn(n,2,1) = - fmatrix_dn(n,v,2) * FOGeometry%Rotations_dn(3,v)
            exactscat_dn(n,3,1) = + fmatrix_dn(n,v,2) * FOGeometry%Rotations_dn(4,v)
            exactscat_dn(n,1,2) = + fmatrix_dn(n,v,2) * FOGeometry%Rotations_dn(1,v)
            exactscat_dn(n,1,3) = - fmatrix_dn(n,v,2) * FOGeometry%Rotations_dn(2,v)
            exactscat_dn(n,2,2) = + help3c1 * FOGeometry%Rotations_dn(3,v) - help4s1 * FOGeometry%Rotations_dn(4,v)
            exactscat_dn(n,2,3) = - help3s1 * FOGeometry%Rotations_dn(3,v) - help4c1 * FOGeometry%Rotations_dn(4,v)
            exactscat_dn(n,3,2) = + help3c1 * FOGeometry%Rotations_dn(4,v) + help4s1 * FOGeometry%Rotations_dn(3,v)
            exactscat_dn(n,3,3) = - help3s1 * FOGeometry%Rotations_dn(4,v) + help4c1 * FOGeometry%Rotations_dn(3,v)
            if ( nstokes .eq. 4 ) then
              exactscat_dn(n,2,4) = - fmatrix_dn(n,v,5) * FOGeometry%Rotations_dn(4,v) 
              exactscat_dn(n,4,2) = - fmatrix_dn(n,v,5) * FOGeometry%Rotations_dn(2,v) 
              exactscat_dn(n,3,4) = + fmatrix_dn(n,v,5) * FOGeometry%Rotations_dn(3,v) 
              exactscat_dn(n,4,3) = - fmatrix_dn(n,v,5) * FOGeometry%Rotations_dn(1,v) 
              exactscat_dn(n,4,4) = + fmatrix_dn(n,v,6)
            endif
            exactscat_dn(n,1:ns,1:ns) = tms(n)*exactscat_dn(n,1:ns,1:ns)
          endif
        enddo
      endif

!  Attenuations and Solar solutions
!  ================================

!  Initialize

      Attenuations   = zero ; suntau = zero   ; Attenuationsfine   = zero
      Attenuations_p = zero ; suntau_p = zero ; Attenuationsfine_p = zero
      nstart = nlayers

!  Critical removed TRC
!  Initialize, only to layer Ncrit if applicable
!     if (Ncrit(v).ne.0) nstart = nCrit(v)

!  Attenuations to End points (including TOA). Both PS representations
!    MUST go all the way to NLAYERS (surface term required)

      do n = 0, nlayers
         nt = FOGeometry%ntraverse_dn(n,v) ; sumd = dot_product(extinction(1:nt),FOGeometry%sunpaths_dn(n,1:nt,v))
         suntau(n) = sumd    ; If (sumd .lt. Expcutoff ) Attenuations(n) = exp( - sumd )
      enddo

!  RobFix 8/20/16. Attenuations to partial-layer points

      if ( do_Partials ) then
        do ut = 1, npartials
          nt = FOGeometry%ntraverse_p_dn(ut,v) ; sumd = dot_product(extinction(1:nt),FOGeometry%sunpaths_p_dn(ut,1:nt,v))
          suntau_p(ut) = sumd    ; If (sumd .lt. Expcutoff ) Attenuations_p(ut) = exp( - sumd )
        enddo
      endif

!  Adjust nstart

      do n = 1, nlayers
         if ( layermask_dn(n) .and. attenuations(n-1).ne.zero )  nstart = n
      enddo

!  Enhanced-spherical, fine-layer attenuations, Whole-layer integration

      if ( do_enhanced_ps ) then
        do n = 1, nlayers
          if ( layermask_dn(n) .and. do_sources_dn(n,v) ) then
            do j = 1, FOGeometry%nfinedivs(n,v)
              nt = FOGeometry%ntraversefine_dn(n,j,v)
              sumd = dot_product(extinction(1:nt),FOGeometry%sunpathsfine_dn(n,1:nt,j,v))
              if (sumd .lt. Expcutoff ) Attenuationsfine(n,j) = exp( - sumd )
            enddo
          endif
        enddo
      endif

!  RobFix 8/20/16. Enhanced-spherical, fine-layer attenuations, Partial layer integration

      if ( do_enhanced_ps .and. do_Partials ) then
        do ut = 1, npartials
          if ( do_sources_dn_p(ut,v) ) then
            np = partial_layeridx(ut)
            do j = 1, FOGeometry%nfinedivs_p_dn(ut,v)
              nt = FOGeometry%ntraversefine_p_dn(ut,j,v)
              sumd = dot_product(extinction(1:nt),FOGeometry%sunpathsfine_p_dn(ut,1:nt,j,v))
              If (sumd .lt. Expcutoff ) Attenuationsfine_p(ut,j) = exp( - sumd )
            enddo
          endif
        enddo
      endif

!  Layer integrated Solar sources
!  ==============================

!  Plane/Parallel or Regular-PS (Average secant formulation)
!    Special treatment for the horizonal case --> Factor2 = 0, lostrans = 0

      if ( do_RegPSorPP ) then
        do n = nlayers, 1, -1
          factor1 = zero ; factor2 = zero
          if ( layermask_dn(n) .and. do_sources_dn(n,v)  ) then
            if ( FOGeometry%Mu1_dn(v) .gt. zero ) then
              lostau = deltaus(n)  / FOGeometry%Mu1_dn(v)
              if ( lostau .lt. Expcutoff ) lostrans_dn(n) = exp( - lostau )
              factor1 = Attenuations(n-1)*lostrans_dn(n) - Attenuations(n)
              factor2 = ((suntau(n) - suntau(n-1))/lostau) - one

!  Patch for PC. Use Taylor series expansion, keep only the first-term
              if ( abs(factor2).lt.1.0d-10 ) then
                 multiplier(n) = Attenuations(n) * lostau
              else
                 multiplier(n) = factor1 / factor2
              endif
!  End patch for PC

            endif
          endif
        enddo
      endif


!  RobFix 8/20/16. Plane/Parallel or Regular-PS, Partial-layer output

      if ( do_RegPSorPP .and.do_Partials ) then
        do ut = 1, npartials
          if ( do_sources_dn_p(ut,v) ) then
            np = Partial_layeridx(ut) ; kn = extinction(np)
            path_dn = FOGeometry%LosP_paths(ut,v)
            factor1 = zero ; factor2 = zero
            if ( FOGeometry%Mu1_dn(v) .gt. zero ) then
              lostau = kn * path_dn
              if ( lostau .lt. Expcutoff ) lostrans_dn_p(ut) = exp( - lostau )
              factor1 = Attenuations(np-1)*lostrans_dn_p(ut) - Attenuations_p(ut)
              factor2 = ((suntau_p(ut) - suntau(np-1))/lostau) - one

!  Patch for PC. Use Taylor series expansion, keep only the first-term
              if ( abs(factor2).lt.1.0d-10 ) then
                 multiplier_p(ut) =  Attenuations_p(ut) * lostau
              else
                 multiplier_p(ut) = factor1 / factor2
              endif
!  End patch for PC

            endif
          endif
        enddo
      endif

!  Enhanced PS: General case. RobFix 8/20/16 streamlined code using distances
!   Quadratures measured from the bottom

      if ( do_enhanced_ps ) then
        do n = nlayers, 1, -1
          if ( layermask_dn(n) .and. do_sources_dn(n,v)  ) then
!mick fix 3/22/2017 - replaced index "np" with "n" in "LosW_paths"
            kn = extinction(n) ; path_dn = FOGeometry%LosW_paths(n,v)
            lostau = kn * path_dn ; if( lostau.lt.Expcutoff ) lostrans_dn(n) = exp ( - lostau )
            sum = zero
            do j = 1, FOGeometry%nfinedivs(n,v)
              dj = FOGeometry%LosW_paths(n,v) - FOGeometry%xfine(n,j,v) ; tran = exp ( - kn * dj )
              sum  = sum + attenuationsfine(n,j) * tran * FOGeometry%wfine(n,j,v)
            enddo 
            multiplier(n) = sum * kn
          endif
        enddo
      endif

!  Enhanced PS:  RobFix 8/20/16 Partials

      if ( do_enhanced_ps .and. do_partials ) then
        do ut = 1, npartials
          if ( do_sources_dn_p(ut,v) ) then
            np = partial_layeridx(ut) ; kn = extinction(np)
            path_dn = FOGeometry%LosP_paths(ut,v)
            lostau = kn * path_dn ; if ( lostau.lt.Expcutoff ) lostrans_dn_p(ut) = exp ( - lostau )
            sum = zero
            do j = 1, FOGeometry%nfinedivs_p_dn(ut,v)
              dj = path_dn - FOGeometry%xfine_p_dn(ut,j,v) ; tran = exp ( - kn * dj )
              sum  = sum + attenuationsfine_p(ut,j) * tran * FOGeometry%wfine_p_dn(ut,j,v)
            enddo
            multiplier_p(ut) = sum * kn
          endif
        enddo
      endif

!  Layer sources
!  -------------

!  General case, Whole layers
!mick fix 9/19/2017 - defined "sources_dn" for all stokes vector elements

      do n = 1, nlayers, 1
        if ( layermask_dn(n) .and. do_sources_dn(n,v)  ) then
          if ( do_sunlight ) then
             shelp(1:ns) = exactscat_dn(n,1:ns,1) * fluxvec(1)
          else
            do o1 = 1, nstokes
              shelp(o1) = dot_product(exactscat_dn(n,o1,1:ns),fluxvec(1:ns))
            enddo
          endif
          !sources_dn(n,o1) = shelp(o1) * multiplier(n)
          sources_dn(n,1:ns) = shelp(1:ns) * multiplier(n)
        endif
      enddo

!  Partials case

      if ( do_partials ) then
        do ut = 1, npartials
          if ( do_sources_dn_p(ut,v)  ) then
            np = partial_layeridx(ut)
            if ( do_sunlight ) then
              shelp(1:ns) = exactscat_dn(np,1:ns,1) * fluxvec(1)
            else
              do o1 = 1, nstokes
                shelp(o1) = dot_product(exactscat_dn(np,o1,1:ns),fluxvec(1:ns))
              enddo
            endif
            sources_dn_p(ut,1:ns) = shelp(1:ns)* multiplier_p(ut)
          endif
        enddo
      endif

!  Source function integration
!  ===========================

!  start recursion

      NC =  0
      CUMSOURCE_DN(NC,:,v) = zero
      NSTART = 1 ; NUT_PREV = NSTART - 1

!  Main loop over all output optical depths
!     NLEVEL = Layer index for given optical depth
!     Cumulative source terms : Loop over layers working Downn from NSTART to NUT
!     Check for dndating the recursion

      DO UTA = 1, N_USER_LEVELS
         NUT    = USER_LEVELS(UTA)
         DO N = NSTART, NUT
            NC = N
            do o1 = 1, nstokes
               CUMSOURCE_DN(NC,O1,v)  = LOSTRANS_DN(N) * CUMSOURCE_DN(NC-1,O1,v) + SOURCES_DN(N,O1)
            enddo
         ENDDO
         IF ( Partial_OUTFLAG(UTA) ) THEN
           UT = Partial_OUTINDEX(UTA)
           do o1 = 1, nstokes
             STOKES_DN(UTA,O1,v) = FLUX * ( CUMSOURCE_DN(NC,O1,V) * LOSTRANS_DN_p(UT) + SOURCES_DN_p(UT,O1) )
           enddo
         ELSE
           do o1 = 1, nstokes
             STOKES_DN(UTA,O1,v) = FLUX * CUMSOURCE_DN(NC,O1,v)
           enddo
         ENDIF
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
         NUT_PREV = NUT
      ENDDO

!  Finish geometry loop

   enddo

!  Finish

   return
end subroutine FO_Vector_SSRT_DN

!

subroutine FO_Vector_SSRT_UPDN  &
   ( do_upwelling, do_dnwelling, do_sunlight, do_deltam_scaling, do_lambertian,             & ! Inputs (Flags)
     do_surface_leaving, do_water_leaving, do_Partials, do_PlanPar, do_enhanced_ps,         & ! Inputs (Flags)
     do_sources_up, do_sources_up_p, do_sources_dn, do_sources_dn_p,                        & ! Inputs (Flags/sources)
     nstokes, ngeoms, nlayers, n_user_levels, user_level_mask_up, user_level_mask_dn,       & ! Inputs (control)
     npartials, partial_outindex, partial_outflag, partial_layeridx, FOGeometry, flux,      & ! Inputs (partial/geometry)
     fluxvec, extinction, deltaus, omega, truncfac, fmatrix_up, fmatrix_dn, reflec, slterm, & ! Inputs (Optical)
     stokes_up, stokes_db, cumsource_up, stokes_dn, cumsource_dn, cumtrans )                  ! Outputs

!  Stand-alone routine for Upwelling and Downwelling Solar-beam Single-scatter (SS)
!    computation of Stokes-vector. Inputs: geometry, spherical functions, optical properties.

   implicit none         

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ======

!  flags
!  Version 1.5: --> F-matrix flag added 7/7/16; surface-leaving flag 8/2/16; Partials 8/20/16
!  4/15/20. Version 2.8.2. DO_FMATRIX flag removed; Now we must use FMATRIX (no choice)

   LOGICAL, Intent(in) :: DO_UPWELLING
   LOGICAL, Intent(in) :: DO_DNWELLING

   LOGICAL, Intent(in) :: DO_SUNLIGHT
   LOGICAL, Intent(in) :: DO_DELTAM_SCALING

   LOGICAL, Intent(in) :: DO_LAMBERTIAN
   LOGICAL, Intent(in) :: DO_SURFACE_LEAVING
   LOGICAL, Intent(in) :: DO_WATER_LEAVING

   logical, Intent(in) :: DO_Partials
   LOGICAL, Intent(in) :: DO_PLANPAR
   LOGICAL, Intent(in) :: DO_ENHANCED_PS

!  Existence flags. 8/19/16. Criticality enters here

   logical, Intent(in)    :: do_sources_up       (maxlayers,MAX_GEOMETRIES)
   logical, Intent(in)    :: do_sources_up_p     (MAX_PARTLAYERS,MAX_GEOMETRIES)
   logical, Intent(in)    :: do_sources_dn       (maxlayers,MAX_GEOMETRIES)
   logical, Intent(in)    :: do_sources_dn_p     (MAX_PARTLAYERS,MAX_GEOMETRIES)

!  Numbers

   INTEGER, Intent(in) :: NSTOKES, NGEOMS, NLAYERS, N_USER_LEVELS
   integer, Intent(in) :: USER_LEVEL_MASK_UP ( MAX_USER_LEVELS )
   integer, Intent(in) :: USER_LEVEL_MASK_DN ( MAX_USER_LEVELS )

!  Numbers for Version 1.5: -->  Partial Control added, 8/20/16

   integer, Intent(in) :: Npartials
   integer, Intent(in) :: partial_layeridx(MAX_PARTLAYERS)
   logical, Intent(in) :: partial_outflag ( MAX_USER_LEVELS )
   integer, Intent(in) :: partial_outindex( MAX_USER_LEVELS )

!  Geometrical inputs
!  ==================

   Type(VLIDORT_Geometry_FO), Intent(in) :: FOGeometry

!  optical inputs
!  --------------

!  Solar Flux 

   real(ffp), Intent(in) :: FLUX, fluxvec(4)

!  Atmosphere. Fmatrix input added 7/7/16
!  4/15/20. Version 2.8.2. GREEKMAT removed; Now we must use FMATRIX_UP, FMATRIX_DN

   REAL(ffp), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   REAL(ffp), Intent(in) :: DELTAUS     ( MAXLAYERS )
   REAL(ffp), Intent(in) :: OMEGA       ( MAXLAYERS )
   REAL(ffp), Intent(in) :: TRUNCFAC    ( MAXLAYERS )
   REAL(ffp), Intent(in) :: FMATRIX_UP  ( MAXLAYERS, MAX_GEOMETRIES, 6 )
   REAL(ffp), Intent(in) :: FMATRIX_DN  ( MAXLAYERS, MAX_GEOMETRIES, 6 )

!  Surface reflectivity (Could be the albedo)
!    Surface leaving input added 8/2/16

   real(ffp), Intent(in) :: REFLEC ( 4, 4, MAX_GEOMETRIES )
   real(ffp), Intent(in) :: SLTERM ( 4,    MAX_GEOMETRIES )

!  outputs
!  -------

!  4/9/19. Additional output for the sleave correction

   REAL(ffp), Intent(Out)  :: stokes_up     ( max_user_levels, 4, MAX_GEOMETRIES )
   REAL(ffp), Intent(Out)  :: stokes_db     ( max_user_levels, 4, MAX_GEOMETRIES )
   REAL(ffp), Intent(Out)  :: cumsource_up  ( 0:maxlayers, 4, MAX_GEOMETRIES )
   REAL(ffp), Intent(Out)  :: stokes_dn     ( max_user_levels, 4, MAX_GEOMETRIES )
   REAL(ffp), Intent(Out)  :: cumsource_dn  ( 0:maxlayers, 4, MAX_GEOMETRIES )
   
   REAL(ffp), Intent(Out)  :: cumtrans     ( max_user_levels, MAX_GEOMETRIES )

!  Upwelling
!  ---------

   if ( do_upwelling ) then

      Call FO_Vector_SSRT_UP &
        ( do_sunlight, do_deltam_scaling, do_Lambertian, do_surface_leaving, do_water_leaving, & ! Inputs (Flags-General/Surface)
          do_Partials, do_PlanPar, do_enhanced_ps, do_sources_up, do_sources_up_p,             & ! Inputs(Flags/criticality)
          nstokes, ngeoms, nlayers, n_user_levels, user_level_mask_up,                         & ! Inputs (control/flux)
          npartials, partial_outindex, partial_outflag, partial_layeridx, FOGeometry,          & ! Inputs (control-partial)
          flux, fluxvec, extinction, deltaus, omega, truncfac, fmatrix_up, Reflec, Slterm,     & ! Inputs (Geometry/Optical/surface)
          stokes_up, stokes_db, cumsource_up, cumtrans )                                         ! Outputs

   endif

   if ( do_dnwelling ) then

      Call FO_Vector_SSRT_DN &
        ( do_sunlight, do_deltam_scaling, do_Partials,                                & ! Inputs (Flags)
          do_PlanPar, do_enhanced_ps, do_sources_dn, do_sources_dn_p,                 & ! Inputs (Flags)
          nstokes, ngeoms, nlayers, n_user_levels, user_level_mask_dn,                & ! Inputs (control/flux)
          npartials, partial_outindex, partial_outflag, partial_layeridx, FOGeometry, & ! Inputs (control-partial)
          flux, fluxvec, extinction, deltaus, omega, truncfac, fmatrix_dn,            & ! Inputs (Optical)
          stokes_dn, cumsource_dn )                                                     ! Outputs

   endif

!  Finish

   return
end subroutine FO_Vector_SSRT_UPDN

!  End module

end module FO_Vector_SSRT_m

