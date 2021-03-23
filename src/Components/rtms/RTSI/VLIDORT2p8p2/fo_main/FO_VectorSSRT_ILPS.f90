
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

module FO_Vector_SSRT_ILPS_m

!  FUNCTION
!  ========

!  For a given wavelength, this Module will calculate First-Order upwelling+downwelling
!  First Order Stokes vectors, and any number of LPS Jacobians (profile/surface)

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
!  Version     1.5, with optional F-matrices, Surface leaving, and partials.

!    Version 1.1a, 01 December 2011, R. Spurr, RT Solutions Inc.
!    Version 1.1b, 13 February 2012, R. Spurr, RT Solutions Inc.
!    Version 1.2,  01 June     2012, R. Spurr, RT Solutions Inc.
!    Version 1.3,  29 October  2012, Extension to Observational multiple geometries
!    Version 1.4,  31 July     2013, Extension to Lattice       multiple geometries
!    Version 1.5,  07 July     2016, Optional calculation using F-matrices directly
!    Version 1.5,  02 August   2016. Inclusion of Surface-leaving terms + LSSL weighting functions
!    Version 1.5,  22 August   2016, Partial-layer output
!    Version 1.5.1 09 April    2019, Add the CUMTRANS output, add water-leaving control

!  SUBROUTINES
!  ===========

!  For Solar sources, the subroutines are
!       FO_Vector_SSRT_ILPS_UP   (Upwelling only)
!       FO_Vector_SSRT_ILPS_DN   (Downwelling only)
!       FO_Vector_SSRT_ILPS_UPDN (Upwelling and Downwelling)

!  Dependencies
!  ============

   use VLIDORT_PARS_m          , only : fpk, zero, one, Expcutoff, MAX_USER_LEVELS, MAXLAYERS, &
                                        MAX_PARTLAYERS, MAXFINELAYERS, MAX_GEOMETRIES,         &
                                        MAX_ATMOSWFS, MAX_SURFACEWFS, MAX_SLEAVEWFS
   use VLIDORT_Setups_def_m

!  All three subroutines public

public

contains

subroutine FO_Vector_SSRT_ILPS_UP &
        ( do_sunlight, do_deltam_scaling, do_Lambertian, do_surface_leaving, do_water_leaving,     & ! Inputs (Flags-General/Surface)
          do_Partials, do_PlanPar, do_enhanced_ps, do_sources_up, do_sources_up_p, do_profilewfs,  & ! Inputs(Flags/sources)
          do_surfacewfs, do_sleavewfs, Lvaryflags, Lvarynums, LvaryFmat, n_reflecwfs, n_sleavewfs, & ! Inputs (Control, Jacobian)
          n_surfacewfs, nstokes, ngeoms, nlayers, n_user_levels, user_levels,                      & ! Inputs (control, output)
          npartials, partial_outindex, partial_outflag, partial_layeridx, FOGeometry,              & ! Inputs (patial/Geometry)
          flux, fluxvec, extinction, deltaus, omega, truncfac, fmatrix_up, reflec, slterm,         & ! Inputs (Optical)
          L_extinction, L_deltaus, L_omega, L_truncfac, L_fmatrix_up, LS_reflec, LSSL_slterm,      & ! Inputs (Linearized)
          Stokes_up, Stokes_db, LP_Jacobians_up, LP_Jacobians_db, LS_Jacobians_db, cumtrans, LP_cumtrans )   ! Output

!  FO routine for Upwelling Solar-beam Single-scatter (SS)
!    computation of Stokes vectors and LPS Jacobians. Inputs: geometry, spherical functions, optical properties.

   implicit none         

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ======

!  flags
!  Version 1.5: --> F-matrix Flag added 7/7/16; surface-leaving flag 8/2/16; Partials 8/22/16
!  4/15/20. Version 2.8.2. DO_FMATRIX flag removed; Now we must use FMATRIX input (no choice)

   LOGICAL, Intent(in) :: DO_SUNLIGHT
   LOGICAL, Intent(in) :: DO_DELTAM_SCALING

   LOGICAL, Intent(in) :: DO_LAMBERTIAN
   LOGICAL, Intent(in) :: DO_SURFACE_LEAVING
   logical, Intent(in) :: DO_WATER_LEAVING    ! 4/9/19 added

   logical, Intent(in) :: DO_Partials
   logical, Intent(in) :: DO_PLANPAR
   logical, Intent(in) :: DO_ENHANCED_PS

!  Existence flags. 8/22/16. Criticality enters here

   logical, Intent(in)    :: do_sources_up       (maxlayers,MAX_GEOMETRIES)
   logical, Intent(in)    :: do_sources_up_p     (MAX_PARTLAYERS,MAX_GEOMETRIES)

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

!  Numbers
!  4/15/20. Version 2.8.2. NGREEK_MOMENTS_INPUT removed; Now we must use FMATRIX (no choice)

   INTEGER, Intent(in) :: NSTOKES, NGEOMS, NLAYERS, N_USER_LEVELS
   INTEGER, Intent(in) :: USER_LEVELS ( MAX_USER_LEVELS )

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
!  4/15/20. Version 2.8.2. GREEKMAT removed; Now we must use FMATRIX input (no choice)

   REAL(ffp), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   REAL(ffp), Intent(in) :: DELTAUS     ( MAXLAYERS )
   REAL(ffp), Intent(in) :: OMEGA       ( MAXLAYERS )
   REAL(ffp), Intent(in) :: TRUNCFAC    ( MAXLAYERS )

   REAL(ffp), Intent(in) :: FMATRIX_UP  ( MAXLAYERS, MAX_GEOMETRIES, 6 )

!  Surface reflectivity (Could be the albedo) + linearizations
!    Surface leaving input added 8/2/16

   real(ffp), Intent(in) :: REFLEC ( 4, 4, MAX_GEOMETRIES )
   real(ffp), Intent(in) :: SLTERM ( 4,    MAX_GEOMETRIES )
   real(ffp), Intent(in) :: LS_REFLEC   ( 4, 4, MAX_GEOMETRIES, max_surfacewfs )
   real(ffp), Intent(in) :: LSSL_SLTERM ( 4,    MAX_GEOMETRIES, max_sleavewfs  )

!  Linearized optical inputs. Fmatrix input added 7/7/16
!  4/15/20. Version 2.8.2. L_GREEKMAT removed; Now we must use L_FMATRIX input (no choice)
!  mick fix 3/2/2020 - moved dimension upper bound "max_atmoswfs" from dim 4 to dim 1 in L_FMATRIX_UP

   real(ffp), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_OMEGA       ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_TRUNCFAC    ( MAXLAYERS, max_atmoswfs )

   real(ffp), Intent(in) :: L_FMATRIX_UP  ( max_atmoswfs, MAXLAYERS, MAX_GEOMETRIES, 6 )

!  outputs
!  -------

   REAL(ffp), Intent(Out)  :: stokes_up     ( max_user_levels, 4, MAX_GEOMETRIES )
   REAL(ffp), Intent(Out)  :: stokes_db     ( max_user_levels, 4, MAX_GEOMETRIES )
   real(ffp), Intent(Out)  :: LP_Jacobians_up  ( max_user_levels, 4, MAX_GEOMETRIES, maxlayers, max_atmoswfs )
   real(ffp), Intent(Out)  :: LP_Jacobians_db  ( max_user_levels, 4, MAX_GEOMETRIES, maxlayers, max_atmoswfs )
   real(ffp), Intent(Out)  :: LS_Jacobians_db  ( max_user_levels, 4, MAX_GEOMETRIES, max_surfacewfs )

!  4/9/19. Additional output for the sleave correction

   real(ffp), Intent(out)  :: CUMTRANS    ( max_user_levels, MAX_GEOMETRIES )
   real(ffp), Intent(out)  :: LP_CUMTRANS ( max_user_levels, MAX_GEOMETRIES, maxlayers, max_atmoswfs )

!  LOCAL
!  -----

!  Attenuations. Partials added, 8/22/16

   real(ffp)  :: attenuations      (0:maxlayers)
   real(ffp)  :: LP_attenuations   (0:maxlayers,maxlayers,max_atmoswfs)

   real(ffp)  :: attenuations_p      (MAX_PARTLAYERS)
   real(ffp)  :: LP_attenuations_p   (MAX_PARTLAYERS,maxlayers,max_atmoswfs)

   real(ffp)  :: Attenuationsfine    (maxlayers,maxfinelayers)
   real(ffp)  :: LP_Attenuationsfine (maxlayers,maxfinelayers,maxlayers,max_atmoswfs)

   real(ffp)  :: Attenuationsfine_p    (MAX_PARTLAYERS,maxfinelayers)
   real(ffp)  :: LP_Attenuationsfine_p (MAX_PARTLAYERS,maxfinelayers,maxlayers,max_atmoswfs)

!  Scattering

   real(ffp)  :: tms            (maxlayers)
   real(ffp)  :: exactscat_up   (maxlayers,4,4)
   real(ffp)  :: L_tms          (maxlayers,max_atmoswfs)
   real(ffp)  :: L_exactscat_up (maxlayers,4,4, max_atmoswfs)

!  Source function integration results

   real(ffp)  :: sources_up  (maxlayers,4), sources_up_p  (MAX_PARTLAYERS,4)
   real(ffp)  :: lostrans_up (maxlayers),   lostrans_up_p (MAX_PARTLAYERS)
   real(ffp)  :: LP_sources_up    (maxlayers,4,maxlayers,max_atmoswfs)
   real(ffp)  :: LP_sources_up_p  (MAX_PARTLAYERS,4,maxlayers,max_atmoswfs)
   real(ffp)  :: LP_lostrans_up (maxlayers,max_atmoswfs), LP_lostrans_up_p (MAX_PARTLAYERS,max_atmoswfs)

   real(ffp)  :: multiplier       ( maxlayers )
   real(ffp)  :: LP_multiplier    ( maxlayers, maxlayers, max_atmoswfs )
   real(ffp)  :: multiplier_p     ( MAX_PARTLAYERS )
   real(ffp)  :: LP_multiplier_p  ( MAX_PARTLAYERS, maxlayers, max_atmoswfs )

!  Local cumulative source terms

   real(ffp)  :: cumsource_db      ( 0:maxlayers, 4 )
   real(ffp)  :: cumsource_up      ( 0:maxlayers, 4 )
   real(ffp)  :: L_cumsource       ( 4, max_atmoswfs )
   real(ffp)  :: LS_cumsource      ( 4, max_surfacewfs )

!  Regular_PS or plane-parallel flag

   logical    :: do_RegPSorPP

!  Help

   integer    :: n, ns, k, j, q, q1, o1, v, uta, nstart, nc, nut, nut_prev, Qnums(maxlayers), nt, np, ut
   logical    :: do_regular_ps, layermask_up(maxlayers), Qvary(maxlayers)

   real(ffp)  :: help, sum, kn, tran, factor1, factor2, m4, m4a, rhelp(4), shelp(4), pi4, term1(4), dj, path_up
   real(ffp)  :: L_help, L_sum(maxlayers,max_atmoswfs), L_tran, L_func, L_factor1, L_factor2, sumd, L_sumd
   real(ffp)  :: lostau, L_lostau, L_Shelp

   real(ffp)  :: help3c1, help3s1, help4c1, help4s1
   real(ffp)  :: suntau(0:maxlayers), suntau_p(MAX_PARTLAYERS)
   real(ffp)  :: LP_suntau(0:maxlayers,maxlayers,max_atmoswfs),LP_suntau_p(MAX_PARTLAYERS,maxlayers,max_atmoswfs)
   real(ffp)  :: ctrans, LP_ctrans(max_atmoswfs)

!  Number

!mick fix 9/19/2017 - define pi4 as in VLIDORT_PARS
   !pi4 = acos(-one)/4.0_ffp
   pi4 = acos(-one)*4.0_ffp

!  Zero the output. 4/9/19 include CUMTRANS and LP_CUMTRANS

   STOKES_UP       = zero ; STOKES_DB       = zero ; cumtrans = zero
   LP_JACOBIANS_UP = zero ; LP_JACOBIANS_DB = zero ; LS_JACOBIANS_DB = zero ; LP_cumtrans = zero

!  Regular_PS or plane-parallel flag

   do_regular_ps = .false.
   if ( .not.do_Planpar ) do_regular_ps = .not. do_enhanced_ps
   do_RegPSorPP = (do_regular_ps .or. do_PlanPar)

!  Bookkeeping

   ns = nstokes

!mick fix 3/22/2017 - turned on all values of LAYERMASK_UP
   !NUT = USER_LEVELS(1) + 1
   !LAYERMASK_UP = .false.
   !LAYERMASK_UP(NUT:NLAYERS) = .true.
   LAYERMASK_UP = .true.

!  Linearization bookkeeping

   Qvary = .false. ; QNums = 0
   if ( do_profilewfs ) then
      Qvary(1:nlayers) = Lvaryflags(1:nlayers)
      QNums(1:nlayers) = Lvarynums (1:nlayers)
   endif

!  TMS factors and linearizations

   if ( do_deltam_scaling ) then
      do n = 1, nlayers
         help = one - truncfac(n) * omega(n)
         tms(n) = omega(n) / help
         if ( Qvary(n) ) then
            do q = 1, Qnums(n)
               L_help = - L_truncfac(n,q)*omega(n) - truncfac(n) * L_omega(n,q)
               L_tms(n,q) = ( L_omega(n,q) - tms(n)*L_help ) / help
            enddo
         endif
      enddo
   else
      do n = 1, nlayers
         tms(n) = omega(n)
         if ( Qvary(n) ) then
            do q = 1, Qnums(n)
               L_tms(n,q) = L_omega(n,q)
            enddo
         endif
      enddo
   endif

!  Start Geometry loop

   do v = 1, ngeoms

!  Zero local sources

      lostrans_up    = zero  ; sources_up    = zero ; exactscat_up   = zero ; cumsource_up = zero
      LP_lostrans_up = zero  ; LP_sources_up = zero ; L_exactscat_up = zero

      lostrans_up_p    = zero  ; sources_up_p    = zero
      LP_lostrans_up_p = zero  ; LP_sources_up_p = zero 

!  Scattering functions and Linearization
!  ======================================

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
            if ( Qvary(n) ) then
              do q = 1, Qnums(n)
                if ( LvaryFmat(n,q) ) then
                  L_exactscat_up(n,1,1,q) = L_fmatrix_up(q,n,v,1) * tms(n) + fmatrix_up(n,v,1) * L_tms(n,q)
                else
                  L_exactscat_up(n,1,1,q) = fmatrix_up(n,v,1) * L_tms(n,q)
                endif
              enddo
            endif
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
            if ( Qvary(n) ) then
              do q = 1, Qnums(n)
                if ( LvaryFmat(n,q) ) then
                  L_exactscat_up(n,1,1,q) = + L_fmatrix_up(q,n,v,1)
                  L_exactscat_up(n,2,1,q) = - L_fmatrix_up(q,n,v,2) * FOGeometry%Rotations_up(3,V)
                  L_exactscat_up(n,3,1,q) = + L_fmatrix_up(q,n,v,2) * FOGeometry%Rotations_up(4,V)
                endif
                L_exactscat_up(n,1:ns,1,q) = L_exactscat_up(n,1:ns,1,q) *   tms(n) &
                                             + exactscat_up(n,1:ns,1)   * L_tms(n,q)
              enddo
            endif               
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
            if ( Qvary(n) ) then
              do q = 1, Qnums(n)
                if ( LvaryFmat(n,q) ) then
                  help3c1 = L_fmatrix_up(q,n,v,3) * FOGeometry%Rotations_up(1,V)
                  help3s1 = L_fmatrix_up(q,n,v,3) * FOGeometry%Rotations_up(2,V)
                  help4c1 = L_fmatrix_up(q,n,v,4) * FOGeometry%Rotations_up(1,V)
                  help4s1 = L_fmatrix_up(q,n,v,4) * FOGeometry%Rotations_up(2,V)
                  L_exactscat_up(n,1,1,q) = + L_fmatrix_up(q,n,v,1)
                  L_exactscat_up(n,2,1,q) = - L_fmatrix_up(q,n,v,2) * FOGeometry%Rotations_up(3,V)
                  L_exactscat_up(n,3,1,q) = + L_fmatrix_up(q,n,v,2) * FOGeometry%Rotations_up(4,V)
                  L_exactscat_up(n,2,2,q) = + help3c1 * FOGeometry%Rotations_up(3,V) - help4s1 * FOGeometry%Rotations_up(4,V)
                  L_exactscat_up(n,2,3,q) = - help3s1 * FOGeometry%Rotations_up(3,V) - help4c1 * FOGeometry%Rotations_up(4,V)
                  L_exactscat_up(n,3,2,q) = + help3c1 * FOGeometry%Rotations_up(4,V) + help4s1 * FOGeometry%Rotations_up(3,V)
                  L_exactscat_up(n,3,3,q) = - help3s1 * FOGeometry%Rotations_up(4,V) + help4c1 * FOGeometry%Rotations_up(3,V)
                  if ( nstokes .eq. 4 ) then
                    L_exactscat_up(n,2,4,q) = - L_fmatrix_up(q,n,v,5) * FOGeometry%Rotations_up(4,V) 
                    L_exactscat_up(n,4,2,q) = - L_fmatrix_up(q,n,v,5) * FOGeometry%Rotations_up(2,V) 
                    L_exactscat_up(n,3,4,q) = + L_fmatrix_up(q,n,v,5) * FOGeometry%Rotations_up(3,V) 
                    L_exactscat_up(n,4,3,q) = - L_fmatrix_up(q,n,v,5) * FOGeometry%Rotations_up(1,V) 
                    L_exactscat_up(n,4,4,q) = + L_fmatrix_up(q,n,v,6)
                  endif
                endif
                L_exactscat_up(n,1:ns,1:ns,q) = L_exactscat_up(n,1:ns,1:ns,q) *   tms(n) &
                                                + exactscat_up(n,1:ns,1:ns)   * L_tms(n,q)
              enddo
            endif               
            exactscat_up(n,1:ns,1:ns) = tms(n)*exactscat_up(n,1:ns,1:ns)

!  End calculation

          endif
        enddo
      endif

!  Attenuations
!  ============

!  Initialize

      Attenuations   = zero ; Attenuationsfine    = zero ; suntau = zero
      Attenuations_p = zero ; Attenuationsfine_p  = zero ; suntau_p = zero

      LP_Attenuations   = zero ; LP_Attenuationsfine    = zero ; LP_suntau = zero
      LP_Attenuations_p = zero ; LP_Attenuationsfine_p  = zero ; LP_suntau_p = zero

!mick fix 9/19/2017 - initialize these also
      LP_multiplier   = zero
      LP_multiplier_p = zero

!  Attenuations to End points (including TOA). All representations
!    MUST go all the way to NLAYERS (surface term required)

      do n = 0, nlayers
         nt = FOGeometry%ntraverse_up(n,v) ; sumd = dot_product(extinction(1:nt),FOGeometry%sunpaths_up(n,1:nt,v))
        suntau(n) = sumd    ; If (sumd .lt. Expcutoff ) Attenuations(n) = exp( - sumd )
        if ( do_profilewfs ) then
          do k = 1, nlayers
            if ( Qvary(k) .and. k.le.nt ) then
              do q = 1, Qnums(k)
                LP_suntau(n,k,q) = L_extinction(k,q) * FOGeometry%sunpaths_up(n,k,v)
                LP_Attenuations(n,k,q) = - Attenuations(n) * LP_suntau(n,k,q)
              enddo
            endif
          enddo
        endif
      enddo

!  RobFix 8/22/16. Attenuations to partial-layer points

      if ( do_Partials ) then
        do ut = 1, npartials
          nt = FOGeometry%ntraverse_p_up(ut,v) ; sumd = dot_product(extinction(1:nt),FOGeometry%sunpaths_p_up(ut,1:nt,v))
          suntau_p(ut) = sumd    ; if (sumd .lt. Expcutoff ) Attenuations_p(ut) = exp( - sumd )
          if ( do_profilewfs ) then
            do k = 1, nlayers
              if ( Qvary(k) .and. k.le.nt ) then
                do q = 1, Qnums(k)
                  LP_suntau_p(ut,k,q) = L_extinction(k,q) * FOGeometry%sunpaths_p_up(ut,k,v)
                  LP_Attenuations_p(ut,k,q) = - Attenuations_p(ut) * LP_suntau_p(ut,k,q)
                enddo
              endif
            enddo
          endif
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
              if ( do_profilewfs ) then
                do k = 1, nlayers
                  if ( Qvary(k) .and. k.le.nt ) then
                    do q = 1, Qnums(k)
                      L_sumd = L_extinction(k,q) * FOGeometry%sunpathsfine_up(n,k,j,v)
                      LP_Attenuationsfine(n,j,k,q) = - Attenuationsfine(n,j) * L_sumd 
                    enddo
                  endif
                enddo
              endif
            enddo
          endif
        enddo
      endif

!  RobFix 8/22/16. Enhanced-spherical, fine-layer attenuations, Partial layer integration

      if ( do_enhanced_ps .and. do_Partials ) then
        do ut = 1, npartials
          if ( do_sources_up_p(ut,v) ) then
            np = partial_layeridx(ut)
            do j = 1, FOGeometry%nfinedivs_p_up(ut,v)
              nt = FOGeometry%ntraversefine_p_up(ut,j,v)
              sumd = dot_product(extinction(1:nt),FOGeometry%sunpathsfine_p_up(ut,1:nt,j,v))
              If (sumd .lt. Expcutoff ) Attenuationsfine_p(ut,j) = exp( - sumd )
              if ( do_profilewfs ) then
                do k = 1, nlayers
                  if ( Qvary(k) .and. k.le.nt ) then
                    do q = 1, Qnums(k)
                      L_sumd = L_extinction(k,q) * FOGeometry%sunpathsfine_p_up(ut,k,j,v)
                      LP_Attenuationsfine_p(ut,j,k,q) = - Attenuationsfine_p(ut,j) * L_sumd 
                    enddo
                  endif
                enddo
              endif
            enddo
          endif
        enddo
      endif

!  Layer integrated solar sources
!  ==============================

!  Plane/Parallel or Regular-PS, Whole-layer source terms
!  ------------------------------------------------------

      if ( do_RegPSorPP ) then
        do n = nlayers, 1, -1
          factor1 = zero ; factor2 = zero
          if ( layermask_up(n) .and. do_sources_up(n,v)  ) then

 !  Sources, general case

            if ( FOGeometry%Mu1_up(v) .gt. zero ) then
              lostau = deltaus(n)  / FOGeometry%Mu1_up(v)
              if ( lostau .lt. Expcutoff ) lostrans_up(n) = exp( - lostau )
              factor1 = Attenuations(n-1) - Attenuations(n)*lostrans_up(n)
              factor2 = one + (suntau(n) - suntau(n-1))/lostau
              multiplier(n) = factor1 / factor2
!mick note: where is statement "sources_up(n) = exactscat_up(n) * multiplier"?
              if ( do_profilewfs ) then
                do k = 1, nlayers
                  if ( Qvary(k) .and. k.le.FOGeometry%ntraverse_up(n,v) ) then
                    do q = 1, Qnums(k)
                      L_factor1 = LP_Attenuations(n-1,k,q) - LP_Attenuations(n,k,q)*lostrans_up(n)
                      L_factor2 = ( LP_suntau(n,k,q) - LP_suntau(n-1,k,q) ) / lostau
                      if ( k.eq.n ) then
                        L_lostau            = L_deltaus(n,q) / FOGeometry%Mu1_up(v)
                        LP_lostrans_up(n,q) = - L_lostau * lostrans_up(n)
                        L_factor1 = L_factor1 - Attenuations(n) * LP_lostrans_up(n,q)
!mick fix 9/19/2017 - replaced following line
                        !L_factor2 = L_factor2 - factor2 *  L_lostau / lostau
                        L_factor2 = L_factor2 - (factor2 - one) *  L_lostau / lostau
                      endif
                      LP_multiplier(n,k,q) = ( L_factor1 - multiplier(n) * L_factor2 ) / factor2
                    enddo
                  endif
                enddo
              endif
            endif

!  End whole layers and regular-PS or plane-parallel formulation

          endif
        enddo
      endif

!  RobFix 8/22/16. Plane/Parallel or Regular-PS, Partial-layer output
!  ------------------------------------------------------------------

      if ( do_RegPSorPP .and.do_Partials ) then
        do ut = 1, npartials
          if ( do_sources_up_p(ut,v) ) then          
            np = Partial_layeridx(ut) ; kn = extinction(np)
            path_up = FOGeometry%LosW_paths(np,v) - FOGeometry%LosP_paths(ut,v)
            factor1 = zero ; factor2 = zero

!  Sources, general case

            if ( FOGeometry%Mu1_up(v) .gt. zero ) then
              lostau = kn * path_up
              if ( lostau .lt. Expcutoff ) lostrans_up_p(ut) = exp( - lostau )
              factor1 = Attenuations_p(ut) - Attenuations(np)*lostrans_up_p(ut)
              factor2 = one + (suntau(np) - suntau_p(ut))/lostau
              multiplier_p(ut) = factor1 / factor2
!mick note: where is statement "sources_up_p(ut) = exactscat_up(np) * multiplier"?
              if ( do_profilewfs ) then
                do k = 1, nlayers
                  if ( Qvary(k) .and. k.le.FOGeometry%ntraverse_p_up(ut,v) ) then
                    do q = 1, Qnums(k)
                      L_factor1 = LP_Attenuations_p(ut,k,q) - LP_Attenuations(np,k,q)*lostrans_up_p(ut)
                      L_factor2 = ( LP_suntau(np,k,q) - LP_suntau_p(ut,k,q) ) / lostau
!mick fix 9/19/2017 - replaced index "n" with "np" in IF condition below
                      if ( k.eq.np ) then
                        L_lostau = L_extinction(np,q) * path_up
                        LP_lostrans_up_p(ut,q) = - L_lostau * lostrans_up_p(ut)
                        L_factor1 = L_factor1 - Attenuations(np) * LP_lostrans_up_p(ut,q)
!mick fix 9/19/2017 - replaced following line
                        !L_factor2 = L_factor2 - factor2 *  L_lostau / lostau
                        L_factor2 = L_factor2 - (factor2 - one) *  L_lostau / lostau
                      endif
                      LP_multiplier_p(ut,k,q) = ( L_factor1 - multiplier_p(ut) * L_factor2 ) / factor2
                    enddo
                  endif
                enddo
              endif
            endif

!  End partial layers and regular-PS or plane-parallel formulation

          endif
        enddo
      endif

!  Enhanced PS: General case, whole layers. 
!  ----------------------------------------

!     RobFix 8/22/16 streamlined code using distances
!      Quadratures from Bottom of the layer

      if ( do_enhanced_ps ) then
        do n = nlayers, 1, -1
          if ( layermask_up(n) .and. do_sources_up(n,v)  ) then
!mick fix 3/22/2017 - replaced index "np" with "n" in "LosW_paths"
            kn = extinction(n) ; path_up = FOGeometry%LosW_paths(n,v)
            lostau = kn * path_up ; if( lostau.lt.Expcutoff ) lostrans_up(n) = exp ( - lostau )
            if ( do_profilewfs ) then
              if ( Qvary(n) ) then
                do q = 1, Qnums(n)
                  L_lostau = L_extinction(n,q) * path_up
                  LP_lostrans_up(n,q) = - L_lostau * lostrans_up(n)
                enddo
              endif
            endif
            sum = zero ; L_sum = zero 
            do j = 1, FOGeometry%nfinedivs(n,v)
              dj = FOGeometry%LosW_paths(n,v) - FOGeometry%xfine(n,j,v) ; tran = exp ( - kn * dj )
              sum  = sum + attenuationsfine(n,j) * tran * FOGeometry%wfine(n,j,v)
              if ( do_profilewfs ) then
                do k = 1, nlayers
                  if ( Qvary(k) .and. k.le.FOGeometry%ntraverse_up(n,v) ) then
                    do q = 1, Qnums(k)
                      if ( k.eq.n ) then
                        L_tran = - dj * L_extinction(n,q)
                        L_func = ( LP_attenuationsfine(n,j,k,q) + L_tran * attenuationsfine(n,j) ) * tran
                        L_sum(k,q) = L_sum(k,q) + L_func * FOGeometry%wfine(n,j,v)
                      else
                        L_func = LP_attenuationsfine(n,j,k,q) * tran
                        L_sum(k,q) = L_sum(k,q) + L_func * FOGeometry%wfine(n,j,v)
                      endif
                    enddo
                  endif
                enddo
              endif
            enddo
            multiplier(n) = sum * kn
            if ( do_profilewfs ) then
              do k = 1, nlayers
                if ( Qvary(k) .and. k.le.FOGeometry%ntraverse_up(n,v) ) then
                  do q = 1, Qnums(k)
                    LP_multiplier(n,k,q) = kn * L_sum(k,q)
                    if ( k.eq.n ) LP_multiplier(n,k,q) = LP_multiplier(n,k,q) + sum * L_extinction(n,q)
                  enddo
                endif
              enddo
            endif        
          endif
        enddo
      endif

!  Enhanced PS: General case, partial layers. 
!  -----------------------------------------

!     RobFix 8/22/16 streamlined code using distances
!      Quadratures from Bottom of the layer

      if ( do_enhanced_ps .and. do_partials ) then
        do ut = 1, npartials
          if ( do_sources_up_p(ut,v) ) then
            np = partial_layeridx(ut) ; kn = extinction(np)
            path_up = FOGeometry%LosW_paths(np,v)- FOGeometry%LosP_paths(ut,v)
            lostau = kn * path_up ; if ( lostau.lt.Expcutoff ) lostrans_up_p(ut) = exp ( - lostau )
            if ( do_profilewfs ) then
              if ( Qvary(np) ) then
                do q = 1, Qnums(np)
                  L_lostau = L_extinction(np,q) * path_up
                  LP_lostrans_up_p(ut,q) = - L_lostau * lostrans_up_p(ut)
                enddo
              endif
            endif
            sum = zero ; L_sum = zero 
            do j = 1, FOGeometry%nfinedivs_p_up(ut,v)
              dj = path_up - FOGeometry%xfine_p_up(ut,j,v) ; tran = exp ( - kn * dj )     ! Correct
              sum  = sum + Attenuationsfine_p(ut,j) * tran * FOGeometry%wfine_p_up(ut,j,v)
              if ( do_profilewfs ) then
                do k = 1, nlayers
                  if ( Qvary(k) .and. k.le.FOGeometry%ntraverse_p_up(ut,v) ) then
                    do q = 1, Qnums(k)
                      if ( k.eq.np ) then
                        L_tran = - dj * L_extinction(np,q)
                        L_func = ( LP_attenuationsfine_p(ut,j,k,q) + L_tran * attenuationsfine_p(ut,j) ) * tran
                        L_sum(k,q) = L_sum(k,q) + L_func * FOGeometry%wfine_p_up(ut,j,v)
                      else
                        L_func = LP_attenuationsfine_p(ut,j,k,q) * tran
                        L_sum(k,q) = L_sum(k,q) + L_func * FOGeometry%wfine_p_up(ut,j,v)
                      endif
                    enddo
                  endif
                enddo
              endif
            enddo
            multiplier_p(ut) = sum * kn
            if ( do_profilewfs ) then
              do k = 1, nlayers
                if ( Qvary(k) .and. k.le.FOGeometry%ntraverse_p_up(ut,v) ) then
                  do q = 1, Qnums(k)
                    LP_multiplier_p(ut,k,q) = kn * L_sum(k,q)
                    if ( k.eq.np ) LP_multiplier_p(ut,k,q) = LP_multiplier_p(ut,k,q) + sum * L_extinction(np,q)
                  enddo
                endif
              enddo
            endif
          endif
        enddo
      endif

!  Layer integrated solar sources
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

          if ( do_profilewfs ) then
            do k = 1, nlayers
              if ( Qvary(k) ) then
                do q = 1, Qnums(k)
                  do o1 = 1, nstokes
                    LP_sources_up(n,o1,k,q) = shelp(o1) * LP_multiplier(n,k,q)
                    if ( k.eq.n ) then
                      L_Shelp = dot_product(L_exactscat_up(n,o1,1:ns,q),fluxvec(1:ns))
                      LP_sources_up(n,o1,k,q) = LP_sources_up(n,o1,k,q) + L_Shelp * multiplier(n)
                    endif
                  enddo
                enddo
              endif
            enddo
          endif
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

            if ( do_profilewfs ) then
              do k = 1, nlayers
                if ( Qvary(k) ) then
                  do q = 1, Qnums(k)
                    do o1 = 1, nstokes
                      LP_sources_up_p(ut,o1,k,q) = shelp(o1) * LP_multiplier_p(ut,k,q)
                      if ( k.eq.np ) then
                        L_Shelp = dot_product(L_exactscat_up(np,o1,1:ns,q),fluxvec(1:ns))
                        LP_sources_up_p(ut,o1,k,q) = LP_sources_up_p(ut,o1,k,q) + L_Shelp * multiplier_p(ut)
                      endif
                    enddo
                  enddo
                endif
              enddo
            endif
          endif
        enddo
      endif

!  Source function integration
!  ===========================

!  NLEVEL = Layer index for given optical depth
!  Cumulative source terms : Loop over layers working upwards from NSTART to level NUT,
!  Check for updating the recursion

!  Stokes-vector Main loop over all output optical depths
!          Cumulative source term will be saved

      NC = 0 
      CUMSOURCE_UP(NC,:) = zero
      CUMSOURCE_DB(NC,:) = zero
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  Surface term

      RHELP = zero; M4 = 4.0_ffp * FOGeometry%Mu0_up(V) ; M4A = M4 * attenuations(nlayers)
      if ( DO_LAMBERTIAN ) then
         RHELP(1) = M4 * REFLEC(1,1,V) * fluxvec(1)
         CUMSOURCE_DB(NC,1) = RHELP(1) * attenuations(nlayers)
      else
         do o1 = 1, nstokes
            RHELP(O1) = M4 * dot_product(REFLEC(O1,1:ns,V),fluxvec(1:ns))
            CUMSOURCE_DB(NC,o1) = RHELP(O1) * attenuations(nlayers)
         enddo
      endif

!  Surface-leaving term. Added, 8/2/16
!   -- (modeled after the DBCORRECTION code in Version 2.7)
!   -- 4/9/19. Not done for water-leaving, as need to use adjusted values

     IF ( DO_SURFACE_LEAVING .and. .not. DO_WATER_LEAVING ) THEN
        do o1 = 1, nstokes
           CUMSOURCE_DB(NC,o1) = CUMSOURCE_DB(NC,o1) + PI4 * SLTERM(o1,v)
        enddo
     ENDIF

!  Main loop over all output optical depths
!     NLEVEL = Layer index for given optical depth
!     Cumulative source terms : Loop over layers working upwards from NSTART to level NUT,
!     Check for updating the recursion. RobFix 8/22/16 Partials.

      DO UTA = N_USER_LEVELS, 1, -1
         NUT    = USER_LEVELS(UTA) + 1
         DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            do o1 = 1, nstokes
               CUMSOURCE_DB(NC,O1) = LOSTRANS_UP(N) * CUMSOURCE_DB(NC-1,O1)
               CUMSOURCE_UP(NC,O1) = LOSTRANS_UP(N) * CUMSOURCE_UP(NC-1,O1) + SOURCES_UP(N,O1)
            enddo
         ENDDO
         IF ( Partial_OUTFLAG(UTA) ) THEN
           UT = Partial_OUTINDEX(UTA)
           STOKES_UP(UTA,1:NS,V) = FLUX * ( CUMSOURCE_UP(NC,1:NS) * LOSTRANS_UP_p(UT) + SOURCES_UP_p(UT,1:NS) )
           STOKES_DB(UTA,1:NS,V) = FLUX * CUMSOURCE_DB(NC,1:NS) * LOSTRANS_UP_p(UT)
         ELSE
           STOKES_UP(UTA,1:NS,V) = FLUX * CUMSOURCE_UP(NC,1:NS)
           STOKES_DB(UTA,1:NS,V) = FLUX * CUMSOURCE_DB(NC,1:NS)
         ENDIF
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1 ; NUT_PREV = NUT
      ENDDO

!  Surface WFs. Change to n_reflecwfs, 8/2/16
!mick fix 9/19/2017 - swapped indices for last two dimemsions in LS_REFLEC for
!                     non-Lambertain case

      if ( do_surfacewfs ) then
         LS_CUMSOURCE = zero
         if ( DO_LAMBERTIAN ) then
            LS_cumsource(1,1:n_reflecwfs) = M4A * LS_REFLEC(1,1,v,1:n_reflecwfs)
         else
            do q = 1, n_reflecwfs
               do o1 = 1, nstokes
                  LS_cumsource(o1,q) = M4A * dot_product(LS_REFLEC(o1,1:ns,v,q),fluxvec(1:ns))
               enddo
            enddo
         endif
       endif

!  Sleave WFs. This section added, 8/2/16
!   -- (modeled after the LSSL_DBCORRECTION code in Version 2.7)

      if ( do_surface_leaving .and. do_sleavewfs ) then
         do q = 1, n_sleavewfs
            q1  = q + n_reflecwfs
            do o1 = 1, nstokes
               LS_cumsource(o1,q1) = pi4 * lssl_slterm(o1,v,q)
            enddo
         enddo
      endif

!  Propagation of surface+sleave weighting functions

      if ( do_surfacewfs .or. ( do_surface_leaving.and.do_sleavewfs) ) then
         NSTART = NLAYERS ; NUT_PREV = NSTART + 1
         DO UTA = N_USER_LEVELS, 1, -1
            NUT    = USER_LEVELS(UTA) + 1
            DO N = NSTART, NUT, -1
               do q = 1, n_surfacewfs
                  LS_cumsource(1:ns,q) = LOSTRANS_UP(N) * LS_CUMSOURCE(1:ns,q)
               enddo
            ENDDO
            IF ( Partial_OUTFLAG(UTA) ) THEN
              UT = Partial_OUTINDEX(UTA)
              do q = 1, n_surfacewfs
                LS_JACOBIANS_DB(UTA,1:NS,V,Q) = FLUX * LS_CUMSOURCE(1:NS,Q) * LOSTRANS_UP_P(UT)
              enddo
            ELSE
              do q = 1, n_surfacewfs
                LS_JACOBIANS_DB(UTA,1:NS,V,Q) = FLUX * LS_CUMSOURCE(1:NS,Q)
              enddo
            ENDIF
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1 ; NUT_PREV = NUT
         ENDDO
      endif

!  Profile Wfs (Atmospheric term)

      if ( do_profilewfs ) then
        do k = 1, nlayers
          if ( Qvary(k) ) then
!mick fix 3/22/2017 - initialize NC
            NC = 0
            L_CUMSOURCE = zero
            NSTART = NLAYERS ; NUT_PREV = NSTART + 1
            DO UTA = N_USER_LEVELS, 1, -1
              NUT = USER_LEVELS(UTA) + 1
              DO N = NSTART, NUT, -1
                NC = NLAYERS + 1 - N
                if ( k.eq.n ) then
                  do q = 1, Qnums(k)
                    L_CUMSOURCE(1:ns,Q) = LP_SOURCES_UP(N,1:ns,K,Q) + &
                       LP_LOSTRANS_UP(N,Q) * CUMSOURCE_UP(NC-1,1:ns) + LOSTRANS_UP(N) * L_CUMSOURCE(1:ns,Q)
                  enddo
                else
                  do q = 1, Qnums(k)
                    L_CUMSOURCE(1:ns,Q) = LP_SOURCES_UP(N,1:ns,K,Q) + LOSTRANS_UP(N) * L_CUMSOURCE(1:ns,Q)
                  enddo
                endif

              ENDDO
              IF ( Partial_OUTFLAG(UTA) ) THEN
!mick fix 9/19/2017 - added PARTIAL_LAYERIDX statement ; changed "n" to "np" in IF condition 
                UT = Partial_OUTINDEX(UTA) ; np = PARTIAL_LAYERIDX(UT)
                if ( k.eq.np ) then
                  do q = 1, Qnums(k)
                    Term1(1:NS) = CUMSOURCE_UP(NC,1:ns) * LP_LOSTRANS_UP_P(UT,q) + &
                                   L_CUMSOURCE(1:ns,Q) * LOSTRANS_UP_p(UT)
                    LP_JACOBIANS_UP(UTA,1:ns,V,K,Q) = FLUX * ( TERM1(1:NS) + LP_SOURCES_UP_p(UT,1:ns,k,Q) )
                  enddo
                else
                  do q = 1, Qnums(k)
                    Term1(1:NS) = L_CUMSOURCE(1:ns,Q) * LOSTRANS_UP_p(UT)
                    LP_JACOBIANS_UP(UTA,1:ns,V,K,Q) = FLUX * ( TERM1(1:NS) + LP_SOURCES_UP_p(UT,1:ns,k,Q) )
                  enddo
                endif
              ELSE
                do q = 1, Qnums(k)
                  LP_JACOBIANS_UP(UTA,1:ns,V,K,Q) = FLUX * L_CUMSOURCE(1:ns,Q)
                enddo
              ENDIF
              IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1 ; NUT_PREV = NUT
            ENDDO
          endif
        enddo
      endif

!  Profile Wfs (Direct beam term)

      if ( do_profilewfs ) then
        do k = 1, nlayers
          if ( Qvary(k) ) then
!mick fix 3/22/2017 - initialize NC
            NC = 0
            do q = 1, Qnums(k)
              L_CUMSOURCE(1:ns,q) = RHELP(1:ns) * LP_attenuations(nlayers,k,q)
            enddo
            NSTART = NLAYERS ; NUT_PREV = NSTART + 1
            DO UTA = N_USER_LEVELS, 1, -1
              NUT = USER_LEVELS(UTA) + 1
              DO N = NSTART, NUT, -1
                NC = NLAYERS + 1 - N
                if ( k.eq.n ) then
                  do q = 1, Qnums(k)
                    L_cumsource(1:ns,q) =  LP_LOSTRANS_UP(N,Q) * CUMSOURCE_DB(NC-1,1:ns) + &
                                              LOSTRANS_UP(N)   * L_CUMSOURCE(1:ns,Q)
                  enddo
                else
                  do q = 1, Qnums(k)
                    L_cumsource(1:ns,q) = LOSTRANS_UP(N) * L_CUMSOURCE(1:ns,Q)
                  enddo
                endif
              ENDDO
              IF ( Partial_OUTFLAG(UTA) ) THEN
!mick fix 9/19/2017 - added PARTIAL_LAYERIDX statement ; changed "n" to "np" in IF condition 
                UT = Partial_OUTINDEX(UTA) ; np = PARTIAL_LAYERIDX(UT)
                if ( k.eq.np ) then
                  do q = 1, Qnums(k)
                    Term1(1:NS) = CUMSOURCE_DB(NC,1:ns) * LP_LOSTRANS_UP_P(UT,q) + &
                                    L_CUMSOURCE(1:ns,Q) * LOSTRANS_UP_P(UT)
!mick fix 9/19/2017 - removed LP_SOURCES_UP_p(UT,1:ns,k,Q) term here
                    !LP_JACOBIANS_DB(UTA,1:ns,V,K,Q) = FLUX * ( TERM1(1:NS) + LP_SOURCES_UP_p(UT,1:ns,k,Q) )
                    LP_JACOBIANS_DB(UTA,1:ns,V,K,Q) = FLUX * TERM1(1:NS)
                  enddo
                else
                  do q = 1, Qnums(k)
                    Term1(1:NS) = L_CUMSOURCE(1:ns,Q) * LOSTRANS_UP_p(UT)
                    LP_JACOBIANS_DB(UTA,1:ns,V,k,Q) = FLUX * TERM1(1:NS)
                  enddo
                endif
              ELSE
                do q = 1, Qnums(k)
                  LP_JACOBIANS_DB(UTA,1:ns,V,K,Q) = FLUX * L_CUMSOURCE(1:ns,Q)
                enddo
              ENDIF
              IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1 ; NUT_PREV = NUT
            ENDDO
          endif
        enddo
      endif

!  4/9/19. WATER_LEAVING CASE. Add CUMTRANS calculation
!    Add start of CUMTRANS recursion (CTRANS = 1.0).

!  5/22/20. Version 2.8.2 Upgrades. CUMTRANS calculation was not properly initialized (no Jacobians)
!                                   Was properly initialized for the Profile jacobian output.
   
      if ( do_water_leaving ) then
         if ( do_profilewfs ) then
            do k = 1, nlayers
               if ( Qvary(k) ) then
                  NC = 0 ; LP_Ctrans = zero ; ctrans = one
                  NSTART = NLAYERS ; NUT_PREV = NSTART + 1
                  DO UTA = N_USER_LEVELS, 1, -1
                     NUT = USER_LEVELS(UTA) + 1
                     DO N = NSTART, NUT, -1
                        if ( k.eq.n ) then
                           do q = 1, Qnums(k)
                              LP_Ctrans(q) =  LP_LOSTRANS_UP(N,Q) * CTRANS + LOSTRANS_UP(N) * LP_CTRANS(Q)
                           enddo
                        else
                           do q = 1, Qnums(k)
                              LP_Ctrans(q) = LOSTRANS_UP(N) * LP_CTRANS(Q)
                           enddo
                        endif
                        ctrans = ctrans * LOSTRANS_UP(N)
                     ENDDO
                     IF ( Partial_OUTFLAG(UTA) ) THEN
                        UT = Partial_OUTINDEX(UTA) ; np = PARTIAL_LAYERIDX(ut)
                        if ( k.eq.np ) then
                           do q = 1, Qnums(k)
                              LP_CUMTRANS(uta,v,k,q) = LP_CTRANS(q) * LOSTRANS_UP_p(UT) + CTRANS * LP_LOSTRANS_UP_p(UT,q)
                           enddo
                        else
                           do q = 1, Qnums(k)
                              LP_CUMTRANS(uta,v,k,q) = LP_CTRANS(q) * LOSTRANS_UP_p(UT)
                           enddo
                        endif
                        CUMTRANS(UTA,V) = CTRANS * LOSTRANS_UP_p(UT)
                     ELSE
                        do q = 1, Qnums(k)
                           LP_CUMTRANS(uta,v,k,q) = LP_CTRANS(q) 
                        enddo
                        CUMTRANS(UTA,V) = CTRANS
                     ENDIF
                     IF  ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
                     NUT_PREV = NUT
                  ENDDO !uta loop
               endif !Qvary(k) if
            enddo !k loop
         else
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
                  CUMTRANS(UTA,V) = CTRANS * LOSTRANS_UP_p(UT)
               ELSE
                  CUMTRANS(UTA,V) = CTRANS
               ENDIF
               IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
               NUT_PREV = NUT
            ENDDO
         endif
      endif
      
!  End geometry loop

   enddo

!  Finish

   return
end subroutine FO_Vector_SSRT_ILPS_UP

!

subroutine FO_Vector_SSRT_ILPS_DN &
   ( do_sunlight, do_deltam_scaling, do_Partials, do_PlanPar, do_enhanced_ps,    & ! Inputs (Flags/flux)
     do_sources_dn, do_sources_dn_p, do_profilewfs, Lvaryflags, Lvarynums,       & ! Inputs (Control, Lin )
     LvaryFmat, nstokes, ngeoms, nlayers, n_user_levels, user_levels,            & ! Inputs (Control, Output)
     npartials, partial_outindex, partial_outflag, partial_layeridx, FOGeometry, & ! Inputs (Control, Partial)
     flux, fluxvec, extinction, deltaus, omega, truncfac, fmatrix_dn,            & ! Inputs (Optical)
     L_extinction, L_deltaus, L_omega, L_truncfac, L_fmatrix_dn,                 & ! Inputs (Optical - Lin)
     Stokes_dn, LP_Jacobians_dn )                                                  ! Output

!  FO routine for Downwelling Solar-beam Single-scatter (SS)
!    computation of stokes-vectors and LPS Jacobians. Inputs: geometry, spherical functions, optical properties.

   implicit none         

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ======

!  flags. F-matrix flag added, 7/7/16
!  4/15/20. Version 2.8.2. DO_FMATRIX flag removed; Now we must use FMATRIX input (no choice)

   LOGICAL, Intent(in) :: DO_SUNLIGHT
   LOGICAL, Intent(in) :: DO_DELTAM_SCALING

   logical, Intent(in) :: DO_Partials
   logical, Intent(in) :: DO_PLANPAR
   logical, Intent(in) :: DO_ENHANCED_PS

!  Existence flags. 8/19/16. Criticality enters here

   logical, Intent(in)    :: do_sources_dn       (maxlayers,MAX_GEOMETRIES)
   logical, Intent(in)    :: do_sources_dn_p     (MAX_PARTLAYERS,MAX_GEOMETRIES)

!  Jacobian control

   LOGICAL, Intent(in) :: do_profilewfs
   LOGICAL, Intent(in) :: Lvaryflags(maxlayers)
   INTEGER, Intent(in) :: Lvarynums (maxlayers)
   LOGICAL, Intent(in) :: LvaryFmat (maxlayers,max_atmoswfs)

!  Numbers
!  4/15/20. Version 2.8.2. NGREEK_MOMENTS_INPUT removed; Now we must use FMATRIX (no choice)

   INTEGER, Intent(in) :: NSTOKES, NGEOMS, NLAYERS, N_USER_LEVELS
   INTEGER, Intent(in) :: USER_LEVELS ( MAX_USER_LEVELS )

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
!  4/15/20. Version 2.8.2. GREEKMAT removed; Now we must use FMATRIX input (no choice)

   REAL(ffp), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   REAL(ffp), Intent(in) :: DELTAUS     ( MAXLAYERS )
   REAL(ffp), Intent(in) :: OMEGA       ( MAXLAYERS )
   REAL(ffp), Intent(in) :: TRUNCFAC    ( MAXLAYERS )

   REAL(ffp), Intent(in) :: FMATRIX_DN  ( MAXLAYERS, MAX_GEOMETRIES, 6 )

!  Linearized optical inputs. Fmatrix input added 7/7/16
!  4/15/20. Version 2.8.2. L_GREEKMAT removed; Now we must use L_FMATRIX input (no choice)
!  mick fix 3/2/2020 - moved dimension upper bound "max_atmoswfs" from dim 4 to dim 1 in L_FMATRIX_DN

   real(ffp), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_OMEGA       ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_TRUNCFAC    ( MAXLAYERS, max_atmoswfs )

   REAL(ffp), Intent(in) :: L_FMATRIX_DN  ( max_atmoswfs, MAXLAYERS, MAX_GEOMETRIES, 6 )

!  outputs
!  -------

   real(ffp), Intent(Out)  :: stokes_dn        ( max_user_levels, 4, MAX_GEOMETRIES )
   real(ffp), Intent(Out)  :: LP_Jacobians_dn  ( max_user_levels, 4, MAX_GEOMETRIES, maxlayers, max_atmoswfs )

!  LOCAL
!  -----

!  Attenuations. Partials added, 8/22/16

   real(ffp)  :: attenuations      (0:maxlayers)
   real(ffp)  :: LP_attenuations   (0:maxlayers,maxlayers,max_atmoswfs)

   real(ffp)  :: attenuations_p      (MAX_PARTLAYERS)
   real(ffp)  :: LP_attenuations_p   (MAX_PARTLAYERS,maxlayers,max_atmoswfs)

   real(ffp)  :: Attenuationsfine    (maxlayers,maxfinelayers)
   real(ffp)  :: LP_Attenuationsfine (maxlayers,maxfinelayers,maxlayers,max_atmoswfs)

   real(ffp)  :: Attenuationsfine_p    (MAX_PARTLAYERS,maxfinelayers)
   real(ffp)  :: LP_Attenuationsfine_p (MAX_PARTLAYERS,maxfinelayers,maxlayers,max_atmoswfs)

!  Scattering

   real(ffp)  :: tms            (maxlayers)
   real(ffp)  :: exactscat_dn   (maxlayers,4,4)
   real(ffp)  :: L_tms          (maxlayers,max_atmoswfs)
   real(ffp)  :: L_exactscat_dn (maxlayers,4,4, max_atmoswfs)

!  Source function integration results

   real(ffp)  :: sources_dn  (maxlayers,4), sources_dn_p  (MAX_PARTLAYERS,4)
   real(ffp)  :: lostrans_dn (maxlayers),   lostrans_dn_p (MAX_PARTLAYERS)
   real(ffp)  :: LP_sources_dn    (maxlayers,4,maxlayers,max_atmoswfs)
   real(ffp)  :: LP_sources_dn_p  (MAX_PARTLAYERS,4,maxlayers,max_atmoswfs)
   real(ffp)  :: LP_lostrans_dn (maxlayers,max_atmoswfs), LP_lostrans_dn_p (MAX_PARTLAYERS,max_atmoswfs)

   real(ffp)  :: multiplier       ( maxlayers )
   real(ffp)  :: LP_multiplier    ( maxlayers, maxlayers, max_atmoswfs )
   real(ffp)  :: multiplier_p     ( MAX_PARTLAYERS )
   real(ffp)  :: LP_multiplier_p  ( MAX_PARTLAYERS, maxlayers, max_atmoswfs )

!  Local cumulative source terms

   real(ffp)  :: cumsource_dn      ( 0:maxlayers, 4 )
   real(ffp)  :: L_cumsource       ( 4, max_atmoswfs )

!  Regular_PS or plane-parallel flag

   logical    :: do_RegPSorPP

!  Help

   integer    :: n, ns, k, j, q, o1, v, uta, nstart, nc, nut, nut_prev, Qnums(maxlayers), nt, np, ut
   logical    :: do_regular_ps, layermask_dn(maxlayers), Qvary(maxlayers)

   real(ffp)  :: help, sum, kn, tran, factor1, factor2, shelp(4), term1(4), dj, path_dn
   real(ffp)  :: L_help, L_sum(maxlayers,max_atmoswfs), L_tran, L_func, L_factor1, L_factor2, sumd, L_sumd
   real(ffp)  :: lostau, L_lostau, L_Shelp

   real(ffp)  :: help3c1, help3s1, help4c1, help4s1
   real(ffp)  :: suntau(0:maxlayers), suntau_p(MAX_PARTLAYERS)
   real(ffp)  :: LP_suntau(0:maxlayers,maxlayers,max_atmoswfs),LP_suntau_p(MAX_PARTLAYERS,maxlayers,max_atmoswfs)

!  Zero the output

   STOKES_DN       = zero
   LP_JACOBIANS_DN = zero

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

!  Linearization bookkeeping

   Qvary = .false. ; QNums = 0
   if ( do_profilewfs ) then
      Qvary(1:nlayers) = Lvaryflags(1:nlayers)
      QNums(1:nlayers) = Lvarynums (1:nlayers)
   endif

!  TMS factors and linearizations

   if ( do_deltam_scaling ) then
      do n = 1, nlayers
         help = one - truncfac(n) * omega(n)
         tms(n) = omega(n) / help
         if ( Qvary(n) ) then
            do q = 1, Qnums(n)
               L_help = - L_truncfac(n,q)*omega(n) - truncfac(n) * L_omega(n,q)
               L_tms(n,q) = ( L_omega(n,q) - tms(n)*L_help ) / help
            enddo
         endif
      enddo
   else
      do n = 1, nlayers
         tms(n) = omega(n)
         if ( Qvary(n) ) then
            do q = 1, Qnums(n)
               L_tms(n,q) = L_omega(n,q)
            enddo
         endif
      enddo
   endif

!  Start Geometry loop

   do v = 1, ngeoms

!  Zero local sources

      lostrans_dn    = zero  ; sources_dn    = zero ; exactscat_dn   = zero ; cumsource_dn = zero
      LP_lostrans_dn = zero  ; LP_sources_dn = zero ; L_exactscat_dn = zero

      lostrans_dn_p    = zero  ; sources_dn_p    = zero
      LP_lostrans_dn_p = zero  ; LP_sources_dn_p = zero 
 
!  Scattering functions and Linearization
!  ======================================

!  Version 1.5, F-matrix option introduced. 7/7/16

!  Scalar only
!  -----------

      if ( nstokes .eq. 1 ) then
        do n = 1, nlayers
          if ( layermask_dn(n) ) then
            exactscat_dn(n,1,1) = fmatrix_dn(n,v,1) * tms(n)
            if ( Qvary(n) ) then
              do q = 1, Qnums(n)
                if ( LvaryFmat(n,q) ) then
                  L_exactscat_dn(n,1,1,q) = L_fmatrix_dn(q,n,v,1)* tms(n) + fmatrix_dn(n,v,1) * L_tms(n,q)
                else
                  L_exactscat_dn(n,1,1,q) = fmatrix_dn(n,v,1) * L_tms(n,q)
                endif
              enddo
            endif    
          endif
        enddo
      endif

!  Vector with Sunlight
!  --------------------

!   Apply rotations for Z-matrix, multiply by TMS

      if ( nstokes .gt. 1 .and. do_sunlight ) then
        do n = 1, nlayers
          if ( layermask_dn(n) ) then
            exactscat_dn(n,1,1) = + fmatrix_dn(n,v,1)
            exactscat_dn(n,2,1) = - fmatrix_dn(n,v,2) * FOGeometry%Rotations_dn(3,V)
            exactscat_dn(n,3,1) = + fmatrix_dn(n,v,2) * FOGeometry%Rotations_dn(4,V)
            if ( Qvary(n) ) then
              do q = 1, Qnums(n)
                if ( LvaryFmat(n,q) ) then
                  L_exactscat_dn(n,1,1,q) = + L_fmatrix_dn(q,n,v,1)
                  L_exactscat_dn(n,2,1,q) = - L_fmatrix_dn(q,n,v,2) * FOGeometry%Rotations_dn(3,V)
                  L_exactscat_dn(n,3,1,q) = + L_fmatrix_dn(q,n,v,2) * FOGeometry%Rotations_dn(4,V)
                endif
                L_exactscat_dn(n,1:ns,1,q) = L_exactscat_dn(n,1:ns,1,q) *   tms(n) &
                                             + exactscat_dn(n,1:ns,1)   * L_tms(n,q)
              enddo
            endif               
            exactscat_dn(n,1:ns,1) = tms(n) * exactscat_dn(n,1:ns,1) 
          endif
        enddo
      endif

!  Vector General case
!  -------------------

! USE FULL 4X4 MATRIX; CODE INTRODUCED BUT NOT TESTED, 05 OCTOBER 2010
!   Apply rotations for Z-matrix, multiply by TMS

      if ( nstokes .gt. 1 .and. .not. do_sunlight ) then
        do n = 1, nlayers
          if ( layermask_dn(n) ) then
            help3c1 = fmatrix_dn(n,v,3) * FOGeometry%Rotations_dn(1,V)
            help3s1 = fmatrix_dn(n,v,3) * FOGeometry%Rotations_dn(2,V)
            help4c1 = fmatrix_dn(n,v,4) * FOGeometry%Rotations_dn(1,V)
            help4s1 = fmatrix_dn(n,v,4) * FOGeometry%Rotations_dn(2,V)
            exactscat_dn(n,1,1) = + fmatrix_dn(n,v,1)
            exactscat_dn(n,2,1) = - fmatrix_dn(n,v,2) * FOGeometry%Rotations_dn(3,V)
            exactscat_dn(n,3,1) = + fmatrix_dn(n,v,2) * FOGeometry%Rotations_dn(4,V)
            exactscat_dn(n,2,2) = + help3c1 * FOGeometry%Rotations_dn(3,V) - help4s1 * FOGeometry%Rotations_dn(4,V)
            exactscat_dn(n,2,3) = - help3s1 * FOGeometry%Rotations_dn(3,V) - help4c1 * FOGeometry%Rotations_dn(4,V)
            exactscat_dn(n,3,2) = + help3c1 * FOGeometry%Rotations_dn(4,V) + help4s1 * FOGeometry%Rotations_dn(3,V)
            exactscat_dn(n,3,3) = - help3s1 * FOGeometry%Rotations_dn(4,V) + help4c1 * FOGeometry%Rotations_dn(3,V)
            if ( nstokes .eq. 4 ) then
              exactscat_dn(n,2,4) = - fmatrix_dn(n,v,5) * FOGeometry%Rotations_dn(4,V) 
              exactscat_dn(n,4,2) = - fmatrix_dn(n,v,5) * FOGeometry%Rotations_dn(2,V) 
              exactscat_dn(n,3,4) = + fmatrix_dn(n,v,5) * FOGeometry%Rotations_dn(3,V) 
              exactscat_dn(n,4,3) = - fmatrix_dn(n,v,5) * FOGeometry%Rotations_dn(1,V) 
              exactscat_dn(n,4,4) = + fmatrix_dn(n,v,6)
            endif
            if ( Qvary(n) ) then
              do q = 1, Qnums(n)
                if ( LvaryFmat(n,q) ) then
                  help3c1 = L_fmatrix_dn(q,n,v,3) * FOGeometry%Rotations_dn(1,V)
                  help3s1 = L_fmatrix_dn(q,n,v,3) * FOGeometry%Rotations_dn(2,V)
                  help4c1 = L_fmatrix_dn(q,n,v,4) * FOGeometry%Rotations_dn(1,V)
                  help4s1 = L_fmatrix_dn(q,n,v,4) * FOGeometry%Rotations_dn(2,V)
                  L_exactscat_dn(n,1,1,q) = + L_fmatrix_dn(q,n,v,1)
                  L_exactscat_dn(n,2,1,q) = - L_fmatrix_dn(q,n,v,2) * FOGeometry%Rotations_dn(3,V)
                  L_exactscat_dn(n,3,1,q) = + L_fmatrix_dn(q,n,v,2) * FOGeometry%Rotations_dn(4,V)
                  L_exactscat_dn(n,2,2,q) = + help3c1 * FOGeometry%Rotations_dn(3,V) - help4s1 * FOGeometry%Rotations_dn(4,V)
                  L_exactscat_dn(n,2,3,q) = - help3s1 * FOGeometry%Rotations_dn(3,V) - help4c1 * FOGeometry%Rotations_dn(4,V)
                  L_exactscat_dn(n,3,2,q) = + help3c1 * FOGeometry%Rotations_dn(4,V) + help4s1 * FOGeometry%Rotations_dn(3,V)
                  L_exactscat_dn(n,3,3,q) = - help3s1 * FOGeometry%Rotations_dn(4,V) + help4c1 * FOGeometry%Rotations_dn(3,V)
                  if ( nstokes .eq. 4 ) then
                    L_exactscat_dn(n,2,4,q) = - L_fmatrix_dn(q,n,v,5) * FOGeometry%Rotations_dn(4,V) 
                    L_exactscat_dn(n,4,2,q) = - L_fmatrix_dn(q,n,v,5) * FOGeometry%Rotations_dn(2,V) 
                    L_exactscat_dn(n,3,4,q) = + L_fmatrix_dn(q,n,v,5) * FOGeometry%Rotations_dn(3,V) 
                    L_exactscat_dn(n,4,3,q) = - L_fmatrix_dn(q,n,v,5) * FOGeometry%Rotations_dn(1,V) 
                    L_exactscat_dn(n,4,4,q) = + L_fmatrix_dn(q,n,v,6)
                  endif
                endif
                L_exactscat_dn(n,1:ns,1:ns,q) = L_exactscat_dn(n,1:ns,1:ns,q) *   tms(n) &
                                                + exactscat_dn(n,1:ns,1:ns)   * L_tms(n,q)
              enddo
            endif               
            exactscat_dn(n,1:ns,1:ns) = tms(n)*exactscat_dn(n,1:ns,1:ns)
          endif
        enddo
      endif

!  Attenuations
!  ============

!  Initialize

      Attenuations   = zero ; Attenuationsfine    = zero ; suntau = zero
      Attenuations_p = zero ; Attenuationsfine_p  = zero ; suntau_p = zero

      LP_Attenuations   = zero ; LP_Attenuationsfine    = zero ; LP_suntau = zero
      LP_Attenuations_p = zero ; LP_Attenuationsfine_p  = zero ; LP_suntau_p = zero

!mick fix 9/19/2017 - initialize these also
      LP_multiplier   = zero
      LP_multiplier_p = zero

!  Attenuations to End points (including TOA). All representations
!    MUST go all the way to NLAYERS (surface term required)

      do n = 0, nlayers
        nt   = FOGeometry%ntraverse_dn(n,v)
        sumd = dot_product(extinction(1:nt),FOGeometry%sunpaths_dn(n,1:nt,v))
        suntau(n) = sumd    ; If (sumd .lt. Expcutoff ) Attenuations(n) = exp( - sumd )
        if ( do_profilewfs ) then
          do k = 1, nlayers
            if ( Qvary(k) .and. k.le.nt ) then
              do q = 1, Qnums(k)
                LP_suntau(n,k,q) = L_extinction(k,q) * FOGeometry%sunpaths_dn(n,k,v)
                LP_Attenuations(n,k,q) = - Attenuations(n) * LP_suntau(n,k,q)
              enddo
            endif
          enddo
        endif
      enddo

!  RobFix 8/22/16. Attenuations to partial-layer points

      if ( do_Partials ) then
        do ut = 1, npartials
          nt = FOGeometry%ntraverse_p_dn(ut,v) ; sumd = dot_product(extinction(1:nt),FOGeometry%sunpaths_p_dn(ut,1:nt,v))
          suntau_p(ut) = sumd    ; If (sumd .lt. Expcutoff ) Attenuations_p(ut) = exp( - sumd )
          if ( do_profilewfs ) then
            do k = 1, nlayers
              if ( Qvary(k) .and. k.le.nt ) then
                do q = 1, Qnums(k)
                  LP_suntau_p(ut,k,q) = L_extinction(k,q) * FOGeometry%sunpaths_p_dn(ut,k,v)
                  LP_Attenuations_p(ut,k,q) = - Attenuations_p(ut) * LP_suntau_p(ut,k,q)
                enddo
              endif
            enddo
          endif
        enddo
      endif

!  Enhanced-spherical, fine-layer attenuations, Whole-layer integration

      if ( do_enhanced_ps ) then
        do n = 1, nlayers
          if ( layermask_dn(n) .and. do_sources_dn(n,v) ) then
            do j = 1, FOGeometry%nfinedivs(n,v)
              nt = FOGeometry%ntraversefine_dn(n,j,v) ; sumd = dot_product(extinction(1:nt),FOGeometry%sunpathsfine_dn(n,1:nt,j,v))
              if (sumd .lt. Expcutoff ) Attenuationsfine(n,j) = exp( - sumd )
              if ( do_profilewfs ) then
                do k = 1, nlayers
                  if ( Qvary(k) .and. k.le.nt ) then
                    do q = 1, Qnums(k)
                      L_sumd = L_extinction(k,q) * FOGeometry%sunpathsfine_dn(n,k,j,v)
                      LP_Attenuationsfine(n,j,k,q) = - Attenuationsfine(n,j) * L_sumd 
                    enddo
                  endif
                enddo
              endif
            enddo
          endif
        enddo
      endif

!  RobFix 8/22/16. Enhanced-spherical, fine-layer attenuations, Partial layer integration

      if ( do_enhanced_ps .and. do_Partials ) then
        do ut = 1, npartials
          if ( do_sources_dn_p(ut,v) ) then
            np = partial_layeridx(ut)
            do j = 1, FOGeometry%nfinedivs_p_dn(ut,v)
              nt   = FOGeometry%ntraversefine_p_dn(ut,j,v)
              sumd = dot_product(extinction(1:nt),FOGeometry%sunpathsfine_p_dn(ut,1:nt,j,v))
              If (sumd .lt. Expcutoff ) Attenuationsfine_p(ut,j) = exp( - sumd )
              if ( do_profilewfs ) then
                do k = 1, nlayers
                  if ( Qvary(k) .and. k.le.nt ) then
                    do q = 1, Qnums(k)
                      L_sumd = L_extinction(k,q) * FOGeometry%sunpathsfine_p_dn(ut,k,j,v)
                      LP_Attenuationsfine_p(ut,j,k,q) = - Attenuationsfine_p(ut,j) * L_sumd 
                    enddo
                  endif
                enddo
              endif
            enddo
          endif
        enddo
      endif

!  Layer integrated solar sources
!  ==============================

!  Plane/Parallel or Regular-PS, Whole-layer source terms
!  ------------------------------------------------------

      if ( do_RegPSorPP ) then
        do n = nlayers, 1, -1
          factor1 = zero ; factor2 = zero
          if ( layermask_dn(n) .and. do_sources_dn(n,v)  ) then

 !  Sources, general case

            if ( FOGeometry%Mu1_dn(v) .gt. zero ) then
              lostau = deltaus(n)  / FOGeometry%Mu1_dn(v)
              if ( lostau .lt. Expcutoff ) lostrans_dn(n) = exp( - lostau )
              factor1 = Attenuations(n-1)*lostrans_dn(n) - Attenuations(n)
              factor2 = ((suntau(n) - suntau(n-1))/lostau) - one
              multiplier(n) = factor1 / factor2
              if ( do_profilewfs ) then
                do k = 1, nlayers
                  if ( Qvary(k) .and. k.le.FOGeometry%ntraverse_dn(n,v) ) then
                    do q = 1, Qnums(k)
                      L_factor1 = LP_Attenuations(n-1,k,q)*lostrans_dn(n) - LP_Attenuations(n,k,q)
                      L_factor2 = (LP_suntau(n,k,q) - LP_suntau(n-1,k,q))/lostau
                      if ( k.eq.n) then
                        L_lostau            = L_deltaus(n,q) / FOGeometry%Mu1_dn(v)
                        LP_lostrans_dn(n,q) = - L_lostau * lostrans_dn(n)
                        L_factor1 = L_factor1 + Attenuations(n-1)*LP_lostrans_dn(n,q)
                        L_factor2 = L_factor2 - (factor2 + one)*L_lostau/lostau 
                        LP_multiplier(n,k,q) = ( L_factor1 - multiplier(n)*L_factor2 ) / factor2
                      else
                        LP_multiplier(n,k,q) = ( L_factor1 - multiplier(n)*L_factor2 ) / factor2
                      endif
                    enddo
                  endif
                enddo
              endif
            endif

!  End whole layers and regular-PS or plane-parallel formulation

          endif
        enddo
      endif

!  RobFix 8/24/16. Plane/Parallel or Regular-PS, Partial-layer output
!  ------------------------------------------------------------------

      if ( do_RegPSorPP .and.do_Partials ) then
        do ut = 1, npartials
          if ( do_sources_dn_p(ut,v) ) then
            np = Partial_layeridx(ut) ; kn = extinction(np)
            path_dn = FOGeometry%LosP_paths(ut,v)
            factor1 = zero ; factor2 = zero

 !  Sources, general case

            if ( FOGeometry%Mu1_dn(v) .gt. zero ) then
              lostau = kn * path_dn
              if ( lostau .lt. Expcutoff ) lostrans_dn_p(ut) = exp( - lostau )
              factor1 = Attenuations(np-1)*lostrans_dn_p(ut) - Attenuations_p(ut)
              factor2 = ((suntau_p(ut) - suntau(np-1))/lostau) - one
              multiplier_p(ut) = factor1 / factor2
              if ( do_profilewfs ) then
                do k = 1, nlayers
                  if ( Qvary(k) .and. k.le.FOGeometry%ntraverse_p_dn(ut,v) ) then
                    do q = 1, Qnums(k)
                      L_factor1 = LP_Attenuations(np-1,k,q)*lostrans_dn_p(ut) - LP_Attenuations_p(ut,k,q)
                      L_factor2 = (LP_suntau_p(ut,k,q) - LP_suntau(np-1,k,q))/lostau
                      if ( k.eq.np) then
                        L_lostau            = L_extinction(np,q) * path_dn
                        LP_lostrans_dn_p(ut,q) = - L_lostau * lostrans_dn_p(ut)
                        L_factor1 = L_factor1 + Attenuations(np-1)*LP_lostrans_dn_p(ut,q)
                        L_factor2 = L_factor2 - (factor2 + one)*L_lostau/lostau 
                        LP_multiplier_p(ut,k,q) = ( L_factor1 - multiplier_p(ut)*L_factor2 ) / factor2
                      else
                        LP_multiplier_p(ut,k,q) = ( L_factor1 - multiplier_p(ut) * L_factor2 ) / factor2
                      endif
                    enddo
                  endif
                enddo
              endif
            endif

!  End partial layers and regular-PS or plane-parallel formulation

          endif
        enddo
      endif

!  Enhanced PS: General case, whole layers. 
!  ----------------------------------------

!     RobFix 8/24/16 streamlined code using distances
!      Quadratures from Top of the layer

      if ( do_enhanced_ps ) then
        do n = nlayers, 1, -1
          if ( layermask_dn(n) .and. do_sources_dn(n,v)  ) then
!mick fix 3/22/2017 - replaced index "np" with "n" in "LosW_paths"
            kn = extinction(n) ; path_dn = FOGeometry%LosW_paths(n,v)
            lostau = kn * path_dn ; if( lostau.lt.Expcutoff ) lostrans_dn(n) = exp ( - lostau )
            if ( do_profilewfs ) then
              if ( Qvary(n) ) then
                do q = 1, Qnums(n)
                  L_lostau        = L_extinction(n,q) * path_dn
                  LP_lostrans_dn(n,q) = - L_lostau * lostrans_dn(n)
                enddo
              endif
            endif
            sum = zero ; L_sum = zero 
            do j = 1, FOGeometry%nfinedivs(n,v)
              dj = FOGeometry%LosW_paths(n,v) - FOGeometry%xfine(n,j,v) ; tran = exp ( - kn * dj )
              sum  = sum + attenuationsfine(n,j) * tran * FOGeometry%wfine(n,j,v)
              if ( do_profilewfs ) then
                do k = 1, nlayers
                  if ( Qvary(k) .and. k.le.FOGeometry%ntraverse_dn(n,v) ) then
                    do q = 1, Qnums(k)
                      if ( k.eq.n ) then
                        L_tran = - dj * L_extinction(n,q)
                        L_func = ( LP_attenuationsfine(n,j,k,q) + L_tran * attenuationsfine(n,j) ) * tran
                        L_sum(k,q) = L_sum(k,q) + L_func * FOGeometry%wfine(n,j,v)
                      else
                        L_func = LP_attenuationsfine(n,j,k,q)  * tran
                        L_sum(k,q) = L_sum(k,q) + L_func * FOGeometry%wfine(n,j,v)
                      endif
                    enddo
                  endif
                enddo
              endif
            enddo 
            multiplier(n) = sum * kn
            if ( do_profilewfs ) then
              do k = 1, nlayers
                if ( Qvary(k) .and. k.le.FOGeometry%ntraverse_dn(n,v) ) then
                  do q = 1, Qnums(k)
                    LP_multiplier(n,k,q) = kn * L_sum(k,q) 
                    if ( k.eq.n ) LP_multiplier(n,k,q) = LP_multiplier(n,k,q) + sum * L_extinction(n,q)
                  enddo
                endif
              enddo
            endif
          endif
        enddo
      endif

!  Enhanced PS: General case, partial layers. 
!  -----------------------------------------

!     RobFix 8/24/16 streamlined code using distances
!      Quadratures from Bottom of the layer

      if ( do_enhanced_ps .and. do_partials ) then
        do ut = 1, npartials
          if ( do_sources_dn_p(ut,v) ) then
            np = partial_layeridx(ut) ; kn = extinction(np)
            path_dn = FOGeometry%LosP_paths(ut,v)
            lostau = kn * path_dn ; if ( lostau.lt.Expcutoff ) lostrans_dn_p(ut) = exp ( - lostau )
            if ( do_profilewfs ) then
!mick fix 3/22/2017 - replaced index "n" with "np" in "Qvary" and "Qnums"
              if ( Qvary(np) ) then
                do q = 1, Qnums(np)
                  L_lostau        = L_extinction(np,q) * path_dn
                  LP_lostrans_dn_p(ut,q) = - L_lostau * lostrans_dn_p(ut)
                enddo
              endif
            endif
            sum = zero ; L_sum = zero 
            do j = 1, FOGeometry%nfinedivs_p_dn(ut,v)
              dj = path_dn - FOGeometry%xfine_p_dn(ut,j,v) ; tran = exp ( - kn * dj )
              sum  = sum + attenuationsfine_p(ut,j) * tran * FOGeometry%wfine_p_dn(ut,j,v)
              if ( do_profilewfs ) then
                do k = 1, nlayers
                  if ( Qvary(k) .and. k.le.FOGeometry%ntraverse_p_dn(ut,v) ) then
                    do q = 1, Qnums(k)
                      if ( k.eq.np ) then
                        L_tran = - dj * L_extinction(np,q)
                        L_func = ( LP_attenuationsfine_p(ut,j,k,q) + L_tran * attenuationsfine_p(ut,j) ) * tran
                        L_sum(k,q) = L_sum(k,q) + L_func * FOGeometry%wfine_p_dn(ut,j,v)
                      else
                        L_func = LP_attenuationsfine_p(ut,j,k,q)  * tran
                        L_sum(k,q) = L_sum(k,q) + L_func * FOGeometry%wfine_p_dn(ut,j,v)
                      endif
                    enddo
                  endif
                enddo
              endif
            enddo
            multiplier_p(ut) = sum * kn
            if ( do_profilewfs ) then
              do k = 1, nlayers
                if ( Qvary(k).and.k.le.FOGeometry%ntraverse_p_dn(ut,v) ) then
                  do q = 1, Qnums(k)
                    LP_multiplier_p(ut,k,q) = kn * L_sum(k,q) 
                    if ( k.eq.np ) LP_multiplier_p(ut,k,q) = LP_multiplier_p(ut,k,q) + sum * L_extinction(np,q)
                  enddo
                endif
              enddo
            endif
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
          if ( do_profilewfs ) then
            do k = 1, nlayers
              if ( Qvary(k) ) then
                do q = 1, Qnums(k)
                  do o1 = 1, nstokes
                    LP_sources_dn(n,o1,k,q) = shelp(o1) * LP_multiplier(n,k,q)
                    if ( k.eq.n ) then
                      L_Shelp = dot_product(L_exactscat_dn(n,o1,1:ns,q),fluxvec(1:ns))
                      LP_sources_dn(n,o1,k,q) = LP_sources_dn(n,o1,k,q) + L_Shelp * multiplier(n)
                    endif
                  enddo
                enddo
              endif
            enddo
          endif
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
            if ( do_profilewfs ) then
              do k = 1, nlayers
                if ( Qvary(k) ) then
                  do q = 1, Qnums(k)
                    do o1 = 1, nstokes
                      LP_sources_dn_p(ut,o1,k,q) = shelp(o1) * LP_multiplier_p(ut,k,q)
                      if ( k.eq.np ) then
                        L_Shelp = dot_product(L_exactscat_dn(np,o1,1:ns,q),fluxvec(1:ns))
                        LP_sources_dn_p(ut,o1,k,q) = LP_sources_dn_p(ut,o1,k,q) + L_Shelp * multiplier_p(ut)
                      endif
                    enddo
                  enddo
                endif
              enddo
            endif
          endif
        enddo
      endif

!  Source function integration
!  ===========================

!  start recursion

      NC =  0
      CUMSOURCE_DN(NC,:) = zero
      NSTART = 1
      NUT_PREV = NSTART - 1

!  Main Stokes-vector loop over all output optical depths
!     NLEVEL = Layer index for given optical depth
!     Cumulative source terms : Loop over layers working upwards from NSTART to level NUT,
!     Check for updating the recursion. Rob Fix Partials 8/24/16.

      NC = 0
      CUMSOURCE_DN(NC,:) = zero
      NSTART = 1 ; NUT_PREV = NSTART - 1
      DO UTA = 1, N_USER_LEVELS
         NUT    = USER_LEVELS(UTA)
         DO N = NSTART, NUT
            NC = N
            do o1 = 1, nstokes
               CUMSOURCE_DN(NC,O1)  = LOSTRANS_DN(N) * CUMSOURCE_DN(NC-1,O1) + SOURCES_DN(N,O1)
            enddo
         ENDDO
         IF ( Partial_OUTFLAG(UTA) ) THEN
           UT = Partial_OUTINDEX(UTA)
           STOKES_DN(UTA,1:NS,V) = FLUX * ( CUMSOURCE_DN(NC,1:NS) * LOSTRANS_DN_p(UT) + SOURCES_DN_p(UT,1:NS) )
         ELSE
           STOKES_DN(UTA,1:NS,V) = FLUX * CUMSOURCE_DN(NC,1:NS)
         ENDIF
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1 ; NUT_PREV = NUT
      ENDDO

!  Profile Wfs (atmospheric term)

      if ( do_profilewfs ) then
        do k = 1, nlayers
          if ( Qvary(k) ) then
            L_CUMSOURCE = zero
            NSTART = 1 ; NUT_PREV = NSTART - 1
            DO UTA = 1, N_USER_LEVELs
              NUT    = USER_LEVELS(UTA)
              DO N = NSTART, NUT
                NC = N
                if ( k.eq.n ) then
                  do q = 1, Qnums(k)
                    L_cumsource(1:ns,q) = LP_SOURCES_DN(N,1:ns,K,Q)    + &
                            LP_LOSTRANS_DN(N,Q) * CUMSOURCE_DN(NC-1,1:ns) + LOSTRANS_DN(N) * L_CUMSOURCE(1:ns,Q)
                  enddo
                else
                  do q = 1, Qnums(k)
                     L_cumsource(1:ns,q) = LP_SOURCES_DN(N,1:ns,K,Q) +  LOSTRANS_DN(N) * L_CUMSOURCE(1:ns,Q)
                  enddo
                endif
              ENDDO
              IF ( Partial_OUTFLAG(UTA) ) THEN
                UT = Partial_OUTINDEX(UTA)
                if ( k.eq.n ) then
                  do q = 1, Qnums(k)
                    Term1(1:NS) = CUMSOURCE_DN(NC,1:ns) * LP_LOSTRANS_DN_P(UT,q) + &
                                   L_CUMSOURCE(1:ns,Q) * LOSTRANS_DN_p(UT)
                    LP_JACOBIANS_DN(UTA,1:ns,V,K,Q) = FLUX * ( TERM1(1:NS) + LP_SOURCES_DN_p(UT,1:ns,k,Q) )
                  enddo
                else
                  do q = 1, Qnums(k)
                    Term1(1:NS) = L_CUMSOURCE(1:ns,Q) * LOSTRANS_DN_p(UT)
                    LP_JACOBIANS_DN(UTA,1:ns,V,K,Q) = FLUX * ( TERM1(1:NS) + LP_SOURCES_DN_p(UT,1:ns,k,Q) )
                  enddo
                endif
              ELSE
                do q = 1, Qnums(k)
                  LP_JACOBIANS_DN(UTA,1:ns,V,K,Q) = FLUX * L_CUMSOURCE(1:ns,Q)
                enddo
              ENDIF
              IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1  ;  NUT_PREV = NUT
            ENDDO
          endif
        enddo
      endif

!  End geometry loop

   enddo

!  Finish

   return
end subroutine FO_Vector_SSRT_ILPS_DN

!

subroutine FO_Vector_SSRT_ILPS_UPDN &
   ( do_upwelling, do_dnwelling, do_sunlight, do_deltam_scaling, do_Partials, do_lambertian,     & ! Inputs (Flags - General)
     do_surface_leaving, do_water_leaving, do_profilewfs, do_surfacewfs, do_sleavewfs,           & ! Inputs (Flags - Surf/Lin)
     do_PlanPar, do_enhanced_ps, do_sources_up, do_sources_up_p, do_sources_dn, do_sources_dn_p, & ! Inputs (Flags - Geometry)
     n_reflecwfs, n_sleavewfs, n_surfacewfs, Lvaryflags, Lvarynums, LvaryFmat,                   & ! Inputs (Control, Lin)
     nstokes, ngeoms, nlayers, n_user_levels, user_level_mask_up, user_level_mask_dn,            & ! Inputs (Control, Output)
     npartials, partial_outindex, partial_outflag, partial_layeridx, FOGeometry,                 & ! Inputs (Control, Partial)
     flux, fluxvec, extinction, deltaus, omega, truncfac, fmatrix_up, fmatrix_dn,                & ! Inputs (Flux/Optical)
     L_extinction, L_deltaus, L_omega, L_truncfac, L_fmatrix_up, L_fmatrix_dn,                   & ! Inputs (Optical - Lin)
     reflec, slterm, LS_reflec, LSSL_slterm,                                                     & ! Inputs (Surface/Sleave)
     stokes_up, stokes_db, LP_Jacobians_up, LP_Jacobians_db, LS_Jacobians_db,                    & ! Output (Upwelling)
     stokes_dn, LP_Jacobians_dn, cumtrans, LP_cumtrans )                                           ! Output (cumtrans & dnwell.)

!  FO routine for Upwelling and Downwelling Solar-beam Single-scatter (SS)
!    computation of Stokes-vectors and Profile/Surface Jacobians. Inputs: geometry, spherical functions, optical properties.

   Implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ======

!  flags.
!  Version 1.5: --> F-matrix flag added 7/7/16; surface-leaving flag 8/2/16; Partials 8/24/16
!  4/15/20. Version 2.8.2. DO_FMATRIX flag removed; Now we must use FMATRIX (no choice)

   LOGICAL, Intent(in) :: DO_UPWELLING
   LOGICAL, Intent(in) :: DO_DNWELLING
   LOGICAL, Intent(in) :: DO_SUNLIGHT
   LOGICAL, Intent(in) :: DO_DELTAM_SCALING

   LOGICAL, Intent(in) :: DO_LAMBERTIAN
   LOGICAL, Intent(in) :: DO_SURFACE_LEAVING
   LOGICAL, Intent(in) :: DO_WATER_LEAVING      ! 4/9/19 added

   logical, Intent(in) :: DO_Partials
   logical, Intent(in) :: DO_PLANPAR
   logical, Intent(in) :: DO_ENHANCED_PS

!  Jacobian Flags. do_sleavewfs added 8/2/16

   LOGICAL, Intent(in) :: do_profilewfs
   LOGICAL, Intent(in) :: do_surfacewfs
   LOGICAL, Intent(in) :: do_sleavewfs

! Numbers

   integer, Intent(in) :: NSTOKES, NLAYERS, NGEOMS, N_USER_LEVELS
   integer, Intent(in) :: USER_LEVEL_MASK_UP ( MAX_USER_LEVELS )
   integer, Intent(in) :: USER_LEVEL_MASK_DN ( MAX_USER_LEVELS )

!  Numbers for Version 1.5: -->  Partial Control added, 8/22/16

   integer, Intent(in) :: Npartials
   integer, Intent(in) :: partial_layeridx(MAX_PARTLAYERS)
   logical, Intent(in) :: partial_outflag ( MAX_USER_LEVELS )
   integer, Intent(in) :: partial_outindex( MAX_USER_LEVELS )

!  Jacobian control. Reflec and sleave numbers added, 8/2/16
!    Note that n_surfacewfs = n_reflecwfs + n_sleavewfs

   INTEGER, Intent(in) :: n_reflecwfs
   INTEGER, Intent(in) :: n_sleavewfs
   INTEGER, Intent(in) :: n_surfacewfs
   LOGICAL, Intent(in) :: Lvaryflags(maxlayers)
   INTEGER, Intent(in) :: Lvarynums (maxlayers)
   LOGICAL, Intent(in) :: LvaryFmat (maxlayers,max_atmoswfs)

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

!  Surface reflectivity (Could be the albedo) + linearizations
!    Surface leaving input added 8/2/16

   real(ffp), Intent(in) :: REFLEC ( 4, 4, MAX_GEOMETRIES )
   real(ffp), Intent(in) :: SLTERM ( 4,    MAX_GEOMETRIES )
   real(ffp), Intent(in) :: LS_REFLEC   ( 4,4,MAX_GEOMETRIES, max_surfacewfs )
   real(ffp), Intent(in) :: LSSL_SLTERM ( 4,  MAX_GEOMETRIES, max_sleavewfs  )

!  Linearized optical inputs. Fmatrix input added 7/7/16
!  4/15/20. Version 2.8.2. L_GREEKMAT removed; Now we must use L_FMATRIX_UP, L_FMATRIX_DN

   real(ffp), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_OMEGA       ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_TRUNCFAC    ( MAXLAYERS, max_atmoswfs )

   REAL(ffp), Intent(in) :: L_FMATRIX_UP  ( max_atmoswfs, MAXLAYERS, MAX_GEOMETRIES, 6 )
   REAL(ffp), Intent(in) :: L_FMATRIX_DN  ( max_atmoswfs, MAXLAYERS, MAX_GEOMETRIES, 6 )

!  Existence flags. 8/19/16. Criticality enters here

   logical, Intent(in)    :: do_sources_up       (maxlayers,MAX_GEOMETRIES)
   logical, Intent(in)    :: do_sources_up_p     (MAX_PARTLAYERS,MAX_GEOMETRIES)
   logical, Intent(in)    :: do_sources_dn       (maxlayers,MAX_GEOMETRIES)
   logical, Intent(in)    :: do_sources_dn_p     (MAX_PARTLAYERS,MAX_GEOMETRIES)

!  outputs
!  -------

   real(ffp), Intent(Out)  :: stokes_up        ( max_user_levels, 4, MAX_GEOMETRIES )
   real(ffp), Intent(Out)  :: stokes_db        ( max_user_levels, 4, MAX_GEOMETRIES )
   real(ffp), Intent(Out)  :: LP_Jacobians_up  ( max_user_levels, 4, MAX_GEOMETRIES, maxlayers, max_atmoswfs )
   real(ffp), Intent(Out)  :: LP_Jacobians_db  ( max_user_levels, 4, MAX_GEOMETRIES, maxlayers, max_atmoswfs )
   real(ffp), Intent(Out)  :: LS_Jacobians_db  ( max_user_levels, 4, MAX_GEOMETRIES, max_surfacewfs )

   real(ffp), Intent(Out)  :: stokes_dn        ( max_user_levels, 4, MAX_GEOMETRIES )
   real(ffp), Intent(Out)  :: LP_Jacobians_dn  ( max_user_levels, 4, MAX_GEOMETRIES, maxlayers, max_atmoswfs )

!  4/9/19. Additional output for the sleave correction

   real(ffp), Intent(out)  :: CUMTRANS    ( max_user_levels, MAX_GEOMETRIES )
   real(ffp), Intent(out)  :: LP_CUMTRANS ( max_user_levels, MAX_GEOMETRIES, maxlayers, max_atmoswfs )

!  Upwelling
!  4/9/19. Additional output for the sleave correction

   if ( do_upwelling  ) then
      call FO_Vector_SSRT_ILPS_UP &
        ( do_sunlight, do_deltam_scaling, do_Lambertian, do_surface_leaving, do_water_leaving,     & ! Inputs (Flags-General/Surface)
          do_Partials, do_PlanPar, do_enhanced_ps, do_sources_up, do_sources_up_p, do_profilewfs,  & ! Inputs(Flags/sources)
          do_surfacewfs, do_sleavewfs, Lvaryflags, Lvarynums, LvaryFmat, n_reflecwfs, n_sleavewfs, & ! Inputs (Control, Jacobian)
          n_surfacewfs, nstokes, ngeoms, nlayers, n_user_levels, user_level_mask_up,               & ! Inputs (control, output)
          npartials, partial_outindex, partial_outflag, partial_layeridx, FOGeometry,              & ! Inputs (patial/Geometry)
          flux, fluxvec, extinction, deltaus, omega, truncfac, fmatrix_up, reflec, slterm,         & ! Inputs (Optical)
          L_extinction, L_deltaus, L_omega, L_truncfac, L_fmatrix_up, LS_reflec, LSSL_slterm,      & ! Inputs (Linearized)
          Stokes_up, Stokes_db, LP_Jacobians_up, LP_Jacobians_db, LS_Jacobians_db, cumtrans, LP_cumtrans )   ! Output
   endif

   if ( do_dnwelling  ) then
      call FO_Vector_SSRT_ILPS_DN &
        ( do_sunlight, do_deltam_scaling, do_Partials, do_PlanPar, do_enhanced_ps,    & ! Inputs (Flags/flux)
          do_sources_dn, do_sources_dn_p, do_profilewfs, Lvaryflags, Lvarynums,       & ! Inputs (Control, Lin )
          LvaryFmat, nstokes, ngeoms, nlayers, n_user_levels, user_level_mask_dn,     & ! Inputs (Control, Output)
          npartials, partial_outindex, partial_outflag, partial_layeridx, FOGeometry, & ! Inputs (Control, Partial)
          flux, fluxvec, extinction, deltaus, omega, truncfac, fmatrix_dn,            & ! Inputs (Optical)
          L_extinction, L_deltaus, L_omega, L_truncfac, L_fmatrix_dn,                 & ! Inputs (Optical - Lin)
          Stokes_dn, LP_Jacobians_dn )                                                 ! Output
   endif

!  Finish

   return
end subroutine FO_Vector_SSRT_ILPS_UPDN

!  End module

end module FO_Vector_SSRT_ILPS_m

