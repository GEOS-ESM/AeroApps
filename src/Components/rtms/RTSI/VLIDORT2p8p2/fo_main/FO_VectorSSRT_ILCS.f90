
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

module FO_Vector_SSRT_ILCS_m

!  FUNCTION
!  ========

!  For a given wavelength, this Module will calculate First-Order upwelling+downwelling
!         Stokes vectors and any number of LCS Jacobians (column/surface)

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
!       FO_Vector_SSRT_ILCS_UP   (Upwelling only)
!       FO_Vector_SSRT_ILCS_DN   (Downwelling only)
!       FO_Vector_SSRT_ILCS_UPDN (Upwelling and Downwelling)

!  Dependencies
!  ============

   use VLIDORT_PARS_m          , only : fpk, zero, one, Expcutoff, MAX_USER_LEVELS, MAXLAYERS, &
                                        MAX_PARTLAYERS, MAXFINELAYERS, MAX_GEOMETRIES,         &
                                        MAX_ATMOSWFS, MAX_SURFACEWFS, MAX_SLEAVEWFS
   use VLIDORT_Setups_def_m

!  All three subroutines public

public

contains

subroutine FO_Vector_SSRT_ILCS_UP &
        ( do_sunlight, do_deltam_scaling, do_Lambertian, do_surface_leaving, do_water_leaving,   & ! Inputs (Flags-General/Surface)
          do_Partials, do_PlanPar, do_enhanced_ps, do_sources_up, do_sources_up_p,               & ! Inputs(Flags/criticality)
          do_columnwfs, do_surfacewfs, do_sleavewfs, n_columnwfs, LvaryFmat, n_reflecwfs,        & ! Inputs (Control, Jacobian)
          n_sleavewfs, n_surfacewfs, nstokes, ngeoms, nlayers, n_user_levels, user_levels,       & ! Inputs (control, output)
          npartials, partial_outindex, partial_outflag, partial_layeridx, FOGeometry,            & ! Inputs (patial/Geometry)
          flux, fluxvec, extinction, deltaus, omega, truncfac, fmatrix_up, reflec, slterm,       & ! Inputs (Optical)
          L_extinction, L_deltaus, L_omega, L_truncfac, L_fmatrix_up, LS_reflec, LSSL_slterm,    & ! Inputs (Linearized)
          Stokes_up, Stokes_db, LC_Jacobians_up, LC_Jacobians_db, LS_Jacobians_db, cumtrans, LC_cumtrans )  ! Output

!  FO routine for Upwelling Solar-beam Single-scatter (SS)
!    computation of Stokes vectors and LCS Jacobians. Inputs: geometry, spherical functions, optical properties.

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

   real(ffp), Intent(in) :: FLUX, FLUXVEC(4)

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

   real(ffp), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_OMEGA       ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_TRUNCFAC    ( MAXLAYERS, max_atmoswfs )

   REAL(ffp), Intent(in) :: L_FMATRIX_UP  ( max_atmoswfs, MAXLAYERS, MAX_GEOMETRIES, 6 )

!  outputs
!  -------

   REAL(ffp), Intent(Out)  :: stokes_up     ( max_user_levels, 4, MAX_GEOMETRIES )
   REAL(ffp), Intent(Out)  :: stokes_db     ( max_user_levels, 4, MAX_GEOMETRIES )
   real(ffp), Intent(Out)  :: LC_Jacobians_up  ( max_user_levels, 4, MAX_GEOMETRIES, max_atmoswfs )
   real(ffp), Intent(Out)  :: LC_Jacobians_db  ( max_user_levels, 4, MAX_GEOMETRIES, max_atmoswfs )
   real(ffp), Intent(Out)  :: LS_Jacobians_db  ( max_user_levels, 4, MAX_GEOMETRIES, max_surfacewfs )

!  4/9/19. Additional output for the sleave correction

   real(ffp), Intent(out)  :: CUMTRANS    ( max_user_levels, MAX_GEOMETRIES )
   real(ffp), Intent(out)  :: LC_CUMTRANS ( max_user_levels, MAX_GEOMETRIES, max_atmoswfs )

!  LOCAL
!  -----

!  Attenuations. Partials added, 8/22/16

   real(ffp)  :: attenuations      (0:maxlayers)
   real(ffp)  :: LC_attenuations   (0:maxlayers,max_atmoswfs)

   real(ffp)  :: attenuations_p      (MAX_PARTLAYERS)
   real(ffp)  :: LC_attenuations_p   (MAX_PARTLAYERS,max_atmoswfs)

   real(ffp)  :: Attenuationsfine    (maxlayers,maxfinelayers)
   real(ffp)  :: LC_Attenuationsfine (maxlayers,maxfinelayers,max_atmoswfs)

   real(ffp)  :: Attenuationsfine_p    (MAX_PARTLAYERS,maxfinelayers)
   real(ffp)  :: LC_Attenuationsfine_p (MAX_PARTLAYERS,maxfinelayers,max_atmoswfs)

!  Scattering

   real(ffp)  :: tms            (maxlayers)
   real(ffp)  :: exactscat_up   (maxlayers,4,4)
   real(ffp)  :: L_tms          (maxlayers,max_atmoswfs)
   real(ffp)  :: L_exactscat_up (maxlayers,4,4, max_atmoswfs)

!  Source function integration results

   real(ffp)  :: sources_up  (maxlayers,4), sources_up_p  (MAX_PARTLAYERS,4)
   real(ffp)  :: lostrans_up (maxlayers),   lostrans_up_p (MAX_PARTLAYERS)
   real(ffp)  :: LC_sources_up  (maxlayers,4,max_atmoswfs), LC_sources_up_p  (MAX_PARTLAYERS,4,max_atmoswfs)
   real(ffp)  :: LC_lostrans_up (maxlayers,max_atmoswfs)  , LC_lostrans_up_p (MAX_PARTLAYERS,max_atmoswfs)

   real(ffp)  :: multiplier       ( maxlayers )
   real(ffp)  :: LC_multiplier    ( maxlayers, max_atmoswfs )
   real(ffp)  :: multiplier_p     ( MAX_PARTLAYERS )
   real(ffp)  :: LC_multiplier_p  ( MAX_PARTLAYERS, max_atmoswfs )

!  Local cumulative source terms

   real(ffp)  :: cumsource_db      ( 0:maxlayers, 4 )
   real(ffp)  :: cumsource_up      ( 0:maxlayers, 4 )

   real(ffp)  :: L_cumsource       ( 4, max_atmoswfs )
   real(ffp)  :: LS_cumsource      ( 4, max_surfacewfs )

!  Regular_PS or plane-parallel flag

   logical    :: do_RegPSorPP

!  Help

   integer    :: n, ns, j, q, q1, v, o1, uta, nstart, nc, nut, nut_prev, Qnums(maxlayers), nt, np, ut
   logical    :: do_regular_ps, layermask_up(maxlayers), Qvary(maxlayers)

   real(ffp)  :: help, sum, kn, tran, factor1, factor2, m4, m4a, rhelp(4), shelp(4), pi4, term1(4)
   real(ffp)  :: L_help, L_sum(max_atmoswfs), L_tran, L_func, L_factor1, L_factor2, sumd, L_sumd, dj, path_up
   real(ffp)  :: lostau, L_lostau, L_Shelp

   real(ffp)  :: help3c1, help3s1, help4c1, help4s1
   real(ffp)  :: suntau(0:maxlayers), suntau_p(MAX_PARTLAYERS)
   real(ffp)  :: LC_suntau(0:maxlayers,max_atmoswfs),LC_suntau_p(MAX_PARTLAYERS,max_atmoswfs)
   real(ffp)  :: ctrans, LC_ctrans(max_atmoswfs)

!  Number

!mick fix 9/19/2017 - define pi4 as in VLIDORT_PARS
   !pi4 = acos(-one)/4.0_ffp
   pi4 = acos(-one)*4.0_ffp

!  Zero the output. 4/9/19 include CUMTRANS and LC_CUMTRANS

   STOKES_UP       = zero ; STOKES_DB       = zero ; cumtrans = zero
   LC_JACOBIANS_UP = zero ; LC_JACOBIANS_DB = zero ; LS_JACOBIANS_DB = zero ; LC_cumtrans = zero

!  Regular_PS or plane-parallel flag

   do_regular_ps = .false.
   if ( .not.do_Planpar ) do_regular_ps = .not. do_enhanced_ps
   do_RegPSorPP = (do_regular_ps .or. do_PlanPar)

!  Bookkeeping

   ns = nstokes

!mick fix 3/22/2017  - turned on all values of LAYERMASK_UP
   !NUT = USER_LEVELS(1) + 1
   !LAYERMASK_UP = .false.
   !LAYERMASK_UP(NUT:NLAYERS) = .true.
   LAYERMASK_UP = .true.

!  Linearization bookkeeping

   Qvary = .false. ; QNums = 0
   if ( do_columnwfs ) then
      Qvary(1:nlayers) = .true.
      QNums(1:nlayers) =  n_columnwfs
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
! mick mod 3/22/2017 - turned off (not needed)
! mick fix 9/19/2017 - turned "L_exactscat_up" back on

      !lostrans_up    = zero  ; sources_up    = zero ; exactscat_up   = zero ; cumsource_up = zero
      !LC_lostrans_up = zero  ; LC_sources_up = zero
      L_exactscat_up = zero

      !lostrans_up_p    = zero  ; sources_up_p    = zero
      !LC_lostrans_up_p = zero  ; LC_sources_up_p = zero 

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
          endif
        enddo
      endif

!  Attenuations
!  ============

!  Initialize

      Attenuations   = zero ; Attenuationsfine    = zero ; suntau = zero
      Attenuations_p = zero ; Attenuationsfine_p  = zero ; suntau_p = zero

      LC_Attenuations   = zero ; LC_Attenuationsfine    = zero ; LC_suntau = zero
      LC_Attenuations_p = zero ; LC_Attenuationsfine_p  = zero ; LC_suntau_p = zero

!  Attenuations to End points (including TOA). All representations
!    MUST go all the way to NLAYERS (surface term required)

      do n = 0, nlayers
         nt = FOGeometry%ntraverse_up(n,v) ; sumd = dot_product(extinction(1:nt),FOGeometry%sunpaths_up(n,1:nt,v))
         suntau(n) = sumd    ; If (sumd .lt. Expcutoff ) Attenuations(n) = exp( - sumd )
         if ( do_columnwfs ) then
           do q = 1, n_columnwfs
             L_sumd = dot_product(L_extinction(1:nt,q),FOGeometry%sunpaths_up(n,1:nt,v))
             LC_suntau(n,q) = L_sumd ; LC_Attenuations(n,q) = - Attenuations(n) * L_sumd
           enddo
         endif
      enddo

!  RobFix 8/22/16. Attenuations to partial-layer points

      if ( do_Partials ) then
        do ut = 1, npartials
          nt = FOGeometry%ntraverse_p_up(ut,v) ; sumd = dot_product(extinction(1:nt),FOGeometry%sunpaths_p_up(ut,1:nt,v))
          suntau_p(ut) = sumd    ; If (sumd .lt. Expcutoff ) Attenuations_p(ut) = exp( - sumd )
          if ( do_columnwfs ) then
            do q = 1, n_columnwfs
              L_sumd = dot_product(L_extinction(1:nt,q),FOGeometry%sunpaths_p_up(ut,1:nt,v))
              LC_suntau_p(ut,q) = L_sumd ; LC_Attenuations_p(ut,q) = - Attenuations_p(ut) * L_sumd
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
              if ( do_columnwfs ) then
                do q = 1, n_columnwfs
                  L_sumd = dot_product(L_extinction(1:nt,q),FOGeometry%sunpathsfine_up(n,1:nt,j,v))
                  LC_Attenuationsfine(n,j,q) = - Attenuationsfine(n,j) * L_sumd 
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
              if ( do_columnwfs ) then
                do q = 1, n_columnwfs
                  L_sumd = dot_product(L_extinction(1:nt,q),FOGeometry%sunpathsfine_p_up(ut,1:nt,j,v))
!mick fix 9/19/2017 - replaced line below
                  !LC_Attenuationsfine(ut,j,q) = - Attenuationsfine_p(ut,j) * L_sumd
                  LC_Attenuationsfine_p(ut,j,q) = - Attenuationsfine_p(ut,j) * L_sumd
                enddo
              endif
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
              if ( do_columnwfs ) then
                do q = 1, n_columnwfs
                  L_lostau           = L_deltaus(n,q) / FOGeometry%Mu1_up(v)
                  LC_lostrans_up(n,q) = - L_lostau * lostrans_up(n)
                  L_factor1 = LC_Attenuations(n-1,q) - LC_Attenuations(n,q)*lostrans_up(n)
                  L_factor2 = (LC_suntau(n,q) - LC_suntau(n-1,q))/lostau
                  L_factor1 = L_factor1 - Attenuations(n)*LC_lostrans_up(n,q)
                  L_factor2 = L_factor2 - (factor2 - one)*L_lostau/lostau 
                  LC_multiplier(n,q) = ( L_factor1 - multiplier(n)*L_factor2 ) / factor2
                enddo
              endif
            endif
          endif
        enddo
      endif

!  RobFix 8/22/16. Plane/Parallel or Regular-PS, Partial-layer output
!  ------------------------------------------------------------------

      if ( do_RegPSorPP .and. do_Partials ) then
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
              if ( do_columnwfs ) then
                do q = 1, n_columnwfs
                  L_lostau = L_extinction(np,q) * path_up
                  LC_lostrans_up_p(ut,q) = - L_lostau * lostrans_up_p(ut)
                  L_factor1 = LC_Attenuations_p(ut,q) - LC_Attenuations(np,q)*lostrans_up_p(ut)
                  L_factor2 = (LC_suntau(np,q) - LC_suntau_p(ut,q))/lostau
                  L_factor1 = L_factor1 - Attenuations(np)*LC_lostrans_up_p(ut,q)
                  L_factor2 = L_factor2 - (factor2 - one)*L_lostau/lostau 
                  LC_multiplier_p(ut,q) = ( L_factor1 - multiplier_p(ut)*L_factor2 ) / factor2
                enddo
              endif
            endif
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
            if ( do_columnwfs ) then
              do q = 1, Qnums(n)
                L_lostau        = L_extinction(n,q) * path_up
                LC_lostrans_up(n,q) = - L_lostau * lostrans_up(n)
              enddo
            endif
            sum = zero ; L_sum = zero 
            do j = 1, FOGeometry%nfinedivs(n,v)
              dj = FOGeometry%LosW_paths(n,v) - FOGeometry%xfine(n,j,v) ; tran = exp ( - kn * dj )
              sum  = sum + attenuationsfine(n,j) * tran * FOGeometry%wfine(n,j,v)
              if ( do_columnwfs ) then
                do q = 1, n_columnwfs
                  L_tran = - dj * L_extinction(n,q)
                  L_func = ( LC_attenuationsfine(n,j,q) + L_tran * attenuationsfine(n,j) ) * tran
                  L_sum(q) = L_sum(q) + L_func * FOGeometry%wfine(n,j,v)
                enddo
              endif
            enddo
            multiplier(n) = sum * kn
            if ( do_columnwfs ) then
              do q = 1, n_columnwfs
                LC_multiplier(n,q) = kn * L_sum(q) + sum * L_extinction(n,q)
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
            if ( do_columnwfs ) then
              do q = 1, n_columnwfs
                L_lostau = L_extinction(np,q) * path_up
                LC_lostrans_up_p(ut,q) = - L_lostau * lostrans_up_p(ut)
              enddo
            endif
            sum = zero ; L_sum = zero 
            do j = 1, FOGeometry%nfinedivs_p_up(ut,v)
              dj = path_up - FOGeometry%xfine_p_up(ut,j,v) ; tran = exp ( - kn * dj )     ! Correct
              sum  = sum + Attenuationsfine_p(ut,j) * tran * FOGeometry%wfine_p_up(ut,j,v)
              if ( do_columnwfs ) then
                do q = 1, n_columnwfs
                  L_tran = - dj * L_extinction(np,q)
                  L_func = ( LC_attenuationsfine_p(ut,j,q) + L_tran * attenuationsfine_p(ut,j) ) * tran
                  L_sum(q) = L_sum(q) + L_func * FOGeometry%wfine_p_up(ut,j,v)
                enddo
              endif
            enddo
            multiplier_p(ut) = sum * kn
            if ( do_columnwfs ) then
              do q = 1, n_columnwfs
                LC_multiplier_p(ut,q) = kn * L_sum(q) + sum * L_extinction(np,q)
              enddo
            endif
          endif
        enddo
      endif

!  Layer integrated Solar sources
!  ==============================

!  General case, Whole layers
!mick fix 9/19/2017 - defined "sources_up" for all stokes vector elements

      do n = nlayers, 1, -1
        if ( layermask_up(n) .and.do_sources_up(n,v)  ) then
          if ( do_sunlight ) then
             shelp(1:ns) = exactscat_up(n,1:ns,1) * fluxvec(1)
          else
            do o1 = 1, nstokes
              shelp(o1) = dot_product(exactscat_up(n,o1,1:ns),fluxvec(1:ns))
            enddo
          endif
          !sources_up(n,o1) = shelp(o1) * multiplier(n)
          sources_up(n,1:ns) = shelp(1:ns) * multiplier(n)

          if ( do_columnwfs ) then
            do q = 1, n_columnwfs
              do o1 = 1, nstokes
                LC_sources_up(n,o1,q) = shelp(o1) * LC_multiplier(n,q)
                if ( do_sunlight ) then
                  L_Shelp = L_exactscat_up(n,o1,1,q)*fluxvec(1)
                else
                  L_Shelp = dot_product(L_exactscat_up(n,o1,1:ns,q),fluxvec(1:ns))
                endif
                LC_sources_up(n,o1,q) =  LC_sources_up(n,o1,q) + L_Shelp * multiplier(n)
              enddo
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
            if ( do_columnwfs ) then
              do q = 1, n_columnwfs
                do o1 = 1, nstokes !index o2 in SSCOR
                  LC_sources_up_p(ut,o1,q) = shelp(o1) * LC_multiplier_p(ut,q)
                  if ( do_sunlight ) then
                    L_Shelp = L_exactscat_up(np,o1,1,q)*fluxvec(1)
                  else
                    L_Shelp = dot_product(L_exactscat_up(np,o1,1:ns,q),fluxvec(1:ns))
                  endif
                  LC_sources_up_p(ut,o1,q) =  LC_sources_up_p(ut,o1,q) + L_Shelp * multiplier_p(ut)
                enddo
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

      RHELP = zero; M4 = 4.0_ffp * FOGeometry%Mu0_up(v) ; M4A = M4 * attenuations(nlayers)
      if ( DO_LAMBERTIAN ) then
         RHELP(1) = M4 * REFLEC(1,1,V) * fluxvec(1)
         CUMSOURCE_DB(NC,1) = RHELP(1) * attenuations(nlayers)
      else
         do o1 = 1, nstokes
            RHELP(O1) = M4 * dot_product(REFLEC(O1,1:ns,V),fluxvec(1:ns))
            CUMSOURCE_DB(NC,o1) = RHELP(O1) * attenuations(nlayers)
         enddo
      endif

!  surface-leaving term. Added, 8/2/16
!   -- (modeled after the DBCORRECTION code in Version 2.7)
!   -- 4/9/19. Not done for water-leaving, as need to use adjusted values

      IF ( DO_SURFACE_LEAVING .and. .not. DO_WATER_LEAVING ) THEN
         do o1 = 1, nstokes
            CUMSOURCE_DB(NC,o1) = CUMSOURCE_DB(NC,o1) + PI4 * SLTERM(o1,v)
         enddo
      ENDIF

!  Main Stokes vector loop over all output optical depths
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

!  Column Wfs (Atmospheric term)

      if ( do_columnwfs ) then
!mick fix 3/22/2017 - initialize NC
         NC = 0
         L_CUMSOURCE = zero
         NSTART = NLAYERS ; NUT_PREV = NSTART + 1
         DO UTA = N_USER_LEVELS, 1, -1
            NUT = USER_LEVELS(UTA) + 1
            DO N = NSTART, NUT, -1
               NC = NLAYERS + 1 - N
               do q = 1, n_columnwfs
                  L_cumsource(1:ns,q) = LC_SOURCES_UP(N,1:ns,Q)         + &
                                LC_LOSTRANS_UP(N,Q) * CUMSOURCE_UP(NC-1,1:ns) + &
                                   LOSTRANS_UP(N)   * L_CUMSOURCE(1:ns,Q)
               enddo
            ENDDO
            IF ( Partial_OUTFLAG(UTA) ) THEN
              UT = Partial_OUTINDEX(UTA)
              do q = 1, n_columnwfs
                Term1(1:NS) = CUMSOURCE_UP(NC,1:ns) * LC_LOSTRANS_UP_P(UT,q) + &
                              L_CUMSOURCE(1:ns,Q) * LOSTRANS_UP_p(UT)
                LC_JACOBIANS_UP(UTA,1:ns,V,Q) = FLUX * ( TERM1(1:NS) + LC_SOURCES_UP_p(UT,1:ns,Q) )
              enddo
            ELSE
              do q = 1, n_columnwfs
                LC_JACOBIANS_UP(UTA,1:ns,V,Q) = FLUX * L_CUMSOURCE(1:ns,Q)
              enddo
            ENDIF
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
            NUT_PREV = NUT
         ENDDO
      endif

!  Column Wfs (Direct beam term)

      if ( do_columnwfs ) then
!mick fix 3/22/2017 - initialize NC
         NC = 0
         do q = 1, n_columnwfs
            L_CUMSOURCE(1:ns,q) = RHELP(1:ns) * LC_attenuations(nlayers,q)
         enddo
         NSTART = NLAYERS  ;  NUT_PREV = NSTART + 1
         DO UTA = N_USER_LEVELS, 1, -1
            NUT = USER_LEVELS(UTA) + 1
            DO N = NSTART, NUT, -1
               NC = NLAYERS + 1 - N
               do q = 1, n_columnwfs
                  L_cumsource(1:ns,q) =  LC_LOSTRANS_UP(N,Q) * CUMSOURCE_DB(NC-1,1:ns) + &
                                            LOSTRANS_UP(N)   * L_CUMSOURCE(1:ns,Q)
               enddo
            ENDDO
            IF ( Partial_OUTFLAG(UTA) ) THEN
              UT = Partial_OUTINDEX(UTA)
              do q = 1, n_columnwfs
                Term1(1:NS) = CUMSOURCE_DB(NC,1:ns) * LC_LOSTRANS_UP_P(UT,q) + &
                                 L_CUMSOURCE(1:ns,Q) * LOSTRANS_UP_p(UT)
                LC_JACOBIANS_DB(UTA,1:ns,V,Q) = FLUX * TERM1(1:NS) 
              enddo
            ELSE
              do q = 1, n_columnwfs
                LC_JACOBIANS_DB(UTA,1:ns,V,Q) = FLUX * L_CUMSOURCE(1:ns,Q)
              enddo
            ENDIF
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1 ; NUT_PREV = NUT
         ENDDO
      endif

!  4/9/19. WATERLEAVING CASE. Add CUMTRANS calculation
!    Add start of CUMTRANS recursion (CTRANS = 1.0).
!  5/22/20. Version 2.8.2 Upgrades. CUMTRANS calculation was not properly initialized

      if ( do_water_leaving ) then
         NSTART   = NLAYERS              ! 5/22/20 Needed for initializing
         NUT_PREV = NSTART + 1           ! 5/22/20 Needed for initializing
         ctrans   = one
         if ( do_columnwfs ) then
            Q = n_columnwfs ; LC_Ctrans = zero
            DO UTA = N_USER_LEVELS, 1, -1
               NUT    = USER_LEVELS(UTA) + 1
               DO N = NSTART, NUT, -1
                  LC_CTRANS(1:q) = LC_CTRANS(1:q) * LOSTRANS_UP(N) + CTRANS * LC_LOSTRANS_UP(N,1:q)
                  CTRANS = CTRANS * LOSTRANS_UP(N)
               ENDDO
               IF ( Partial_OUTFLAG(UTA) ) THEN
                  UT = Partial_OUTINDEX(UTA)
                  LC_CUMTRANS(uta,v,1:q) = LC_CTRANS(1:q) * LOSTRANS_UP_p(UT) + CTRANS * LC_LOSTRANS_UP_p(UT,1:q)
                  CUMTRANS(UTA,V)     = CTRANS * LOSTRANS_UP_p(UT)
               ELSE
                  CUMTRANS(UTA,V)        = CTRANS
                  LC_CUMTRANS(uta,v,1:q) = LC_CTRANS(1:q)
               ENDIF
               IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
               NUT_PREV = NUT
            ENDDO
         else  
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
      
!  End Geometry Loop

   enddo

!  Finish

   return
end subroutine FO_Vector_SSRT_ILCS_UP

!

subroutine FO_Vector_SSRT_ILCS_DN &
      ( do_sunlight, do_deltam_scaling, do_Partials, do_PlanPar, do_enhanced_ps,    & ! Inputs (Flags/flux)
        do_sources_dn, do_sources_dn_p, do_columnwfs, n_columnwfs, LvaryFmat,       & ! Inputs (control, Jacobian )
        nstokes, ngeoms, nlayers, n_user_levels, user_levels,                       & ! Inputs (control)
        npartials, partial_outindex, partial_outflag, partial_layeridx, FOGeometry, & ! Inputs (control-partial)
        flux, fluxvec, extinction, deltaus, omega, truncfac, fmatrix_dn,            & ! Inputs (Optical)
        L_extinction, L_deltaus, L_omega, L_truncfac, L_fmatrix_dn,                 & ! Inputs (Optical - Linearized)
        Stokes_dn, LC_Jacobians_dn )                                                  ! Output
   
!  FO routine for Downwelling Solar-beam Single-scatter (SS)
!    computation of stokes vectors and LCS Jacobians. Inputs: geometry, spherical functions, optical properties.

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ======

!  flags. F-matrix flag added, 7/7/16
!  4/15/20. Version 2.8.2. DO_FMATRIX flag removed; Now we must use FMATRIX (no choice)

   LOGICAL, Intent(in) :: DO_SUNLIGHT
   LOGICAL, Intent(in) :: DO_DELTAM_SCALING

   logical, Intent(in) :: DO_Partials
   logical, Intent(in) :: DO_PLANPAR
   logical, Intent(in) :: DO_ENHANCED_PS

!  Existence flags. 8/19/16. Criticality enters here

   logical, Intent(in)    :: do_sources_dn       (maxlayers,MAX_GEOMETRIES)
   logical, Intent(in)    :: do_sources_dn_p     (MAX_PARTLAYERS,MAX_GEOMETRIES)

!  Jacobian Flag and control

   LOGICAL, Intent(in) :: do_columnwfs
   INTEGER, Intent(in) :: n_columnwfs
   LOGICAL, Intent(in) :: LvaryFmat (maxlayers,max_atmoswfs)

!  Layer and Level Control Numbers, Number of Moments
!  4/15/20. Version 2.8.2. NGREEK_MOMENTS_INPUT removed; Now we must use FMATRIX (no choice)

   integer, Intent(in) :: NLAYERS, NGEOMS, N_USER_LEVELS, NSTOKES
   integer, Intent(in) :: USER_LEVELS ( MAX_USER_LEVELS )

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

   real(ffp), Intent(in) :: FLUX, FLUXVEC(4)

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

   real(ffp), Intent(in) :: L_FMATRIX_DN  ( max_atmoswfs, MAXLAYERS, MAX_GEOMETRIES, 6 )

!  outputs
!  -------

   real(ffp), Intent(Out)  :: Stokes_dn        ( max_user_levels, 4, MAX_GEOMETRIES )
   real(ffp), Intent(Out)  :: LC_Jacobians_dn  ( max_user_levels, 4, MAX_GEOMETRIES, max_atmoswfs )

!  LOCAL
!  -----

!  Attenuations. Partials added, 8/22/16

   real(ffp)  :: attenuations      (0:maxlayers)
   real(ffp)  :: LC_attenuations   (0:maxlayers,max_atmoswfs)

   real(ffp)  :: attenuations_p     (MAX_PARTLAYERS)
   real(ffp)  :: LC_attenuations_p   (MAX_PARTLAYERS,max_atmoswfs)

   real(ffp)  :: Attenuationsfine    (maxlayers,maxfinelayers)
   real(ffp)  :: LC_Attenuationsfine (maxlayers,maxfinelayers,max_atmoswfs)

   real(ffp)  :: Attenuationsfine_p    (MAX_PARTLAYERS,maxfinelayers)
   real(ffp)  :: LC_Attenuationsfine_p (MAX_PARTLAYERS,maxfinelayers,max_atmoswfs)

!  Scattering

   real(ffp)  :: tms            (maxlayers)
   real(ffp)  :: exactscat_dn   (maxlayers,4,4)
   real(ffp)  :: L_tms          (maxlayers,max_atmoswfs)
   real(ffp)  :: L_exactscat_dn (maxlayers,4,4, max_atmoswfs)

!  Source function integration results

   real(ffp)  :: sources_dn  (maxlayers,4), sources_dn_p  (MAX_PARTLAYERS,4)
   real(ffp)  :: lostrans_dn (maxlayers), lostrans_dn_p (MAX_PARTLAYERS)
   real(ffp)  :: LC_sources_dn  (maxlayers,4,max_atmoswfs), LC_sources_dn_p  (MAX_PARTLAYERS,4,max_atmoswfs)
   real(ffp)  :: LC_lostrans_dn (maxlayers,max_atmoswfs)  , LC_lostrans_dn_p (MAX_PARTLAYERS,max_atmoswfs)

   real(ffp)  :: multiplier       ( maxlayers )
   real(ffp)  :: LC_multiplier    ( maxlayers, max_atmoswfs )
   real(ffp)  :: multiplier_p     ( MAX_PARTLAYERS )
   real(ffp)  :: LC_multiplier_p  ( MAX_PARTLAYERS, max_atmoswfs )

!  Local cumulative source terms

   real(ffp)  :: cumsource_dn      ( 0:maxlayers, 4 )
   real(ffp)  :: L_cumsource       ( 4, max_atmoswfs )

!  Regular_PS or plane-parallel flag

   logical    :: do_RegPSorPP

!  Help

   integer    :: n, ns, j, q, v, o1, uta, nstart, nc, nut, nut_prev, Qnums(maxlayers), nt, np, ut
   logical    :: do_regular_ps, layermask_dn(maxlayers), Qvary(maxlayers)

   real(ffp)  :: help, sum, kn, tran, factor1, factor2, shelp(4), term1(4)
   real(ffp)  :: L_help, L_sum(max_atmoswfs), L_tran, L_func, L_factor1, L_factor2, sumd, L_sumd, dj, path_dn
   real(ffp)  :: lostau, L_lostau, L_Shelp

   real(ffp)  :: help3c1, help3s1, help4c1, help4s1
   real(ffp)  :: suntau(0:maxlayers), suntau_p(MAX_PARTLAYERS)
   real(ffp)  :: LC_suntau(0:maxlayers,max_atmoswfs),LC_suntau_p(MAX_PARTLAYERS,max_atmoswfs)

!  Zero the output 

   STOKES_DN       = zero
   LC_JACOBIANS_DN = zero

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
   if ( do_columnwfs ) then
      Qvary(1:nlayers) = .true.
      QNums(1:nlayers) =  n_columnwfs
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

      lostrans_dn   = zero  ; sources_dn    = zero ; exactscat_dn   = zero ; cumsource_dn = zero
      LC_lostrans_dn = zero  ; LC_sources_dn = zero ; L_exactscat_dn = zero 

      lostrans_dn_p   = zero  ; sources_dn_p   = zero
      LC_lostrans_dn_p = zero  ; LC_sources_dn_p = zero 

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
      LC_Attenuations   = zero ; LC_Attenuationsfine    = zero ; LC_suntau = zero
      LC_Attenuations_p = zero ; LC_Attenuationsfine_p  = zero ; LC_suntau_p = zero
      nstart = nlayers

!  Attenuations to End points (including TOA). All representations
!    MUST go all the way to NLAYERS (surface term required)

      do n = 0, nlayers
         nt   = FOGeometry%ntraverse_dn(n,v)
         sumd = dot_product(extinction(1:nt),FOGeometry%sunpaths_dn(n,1:nt,v))
         suntau(n) = sumd    ; if (sumd .lt. Expcutoff ) Attenuations(n) = exp( - sumd )
         if ( do_columnwfs ) then
           do q = 1, n_columnwfs
             L_sumd = dot_product(L_extinction(1:nt,q),FOGeometry%sunpaths_dn(n,1:nt,v))
             LC_suntau(n,q) = L_sumd ; LC_Attenuations(n,q) = - Attenuations(n) * L_sumd
           enddo
         endif
      enddo

!  RobFix 8/22/16. Attenuations to partial-layer points

      if ( do_Partials ) then
        do ut = 1, npartials
          nt = FOGeometry%ntraverse_p_dn(ut,v) ; sumd = dot_product(extinction(1:nt),FOGeometry%sunpaths_p_dn(ut,1:nt,v))
          suntau_p(ut) = sumd    ; if (sumd .lt. Expcutoff ) Attenuations_p(ut) = exp( - sumd )
          if ( do_columnwfs ) then
            do q = 1, n_columnwfs
              L_sumd = dot_product(L_extinction(1:nt,q),FOGeometry%sunpaths_p_dn(ut,1:nt,v))
              LC_suntau_p(ut,q) = L_sumd ; LC_Attenuations_p(ut,q) = - Attenuations_p(ut) * L_sumd
            enddo
          endif
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
              nt = FOGeometry%ntraversefine_dn(n,j,v) ; sumd = dot_product(extinction(1:nt),FOGeometry%sunpathsfine_dn(n,1:nt,j,v))
              if (sumd .lt. Expcutoff ) Attenuationsfine(n,j) = exp( - sumd )
              if ( do_columnwfs ) then
                do q = 1, n_columnwfs
                  L_sumd = dot_product(L_extinction(1:nt,q),FOGeometry%sunpathsfine_dn(n,1:nt,j,v))
                  LC_Attenuationsfine(n,j,q) = - Attenuationsfine(n,j) * L_sumd 
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
              if (sumd .lt. Expcutoff ) Attenuationsfine_p(ut,j) = exp( - sumd )
              if ( do_columnwfs ) then
                do q = 1, n_columnwfs
                  L_sumd = dot_product(L_extinction(1:nt,q),FOGeometry%sunpathsfine_p_dn(ut,1:nt,j,v))
!mick fix 9/19/2017 - replaced line below
                  !LC_Attenuationsfine(ut,j,q) = - Attenuationsfine_p(ut,j) * L_sumd
                  LC_Attenuationsfine_p(ut,j,q) = - Attenuationsfine_p(ut,j) * L_sumd
                enddo
              endif
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
          if ( layermask_dn(n) .and. do_sources_dn(n,v)  ) then

!  Sources, general case

            if ( FOGeometry%Mu1_dn(v) .gt. zero ) then
              lostau = deltaus(n) / FOGeometry%Mu1_dn(v)
              if ( lostau .lt. Expcutoff ) lostrans_dn(n) = exp( - lostau )
              factor1 = Attenuations(n-1)*lostrans_dn(n) - Attenuations(n)
              factor2 = ((suntau(n) - suntau(n-1))/lostau)  - one
              multiplier(n) = factor1 / factor2

              if ( do_columnwfs ) then
                do q = 1, n_columnwfs
                  L_lostau  = L_deltaus(n,q) / FOGeometry%Mu1_dn(v)
                  LC_lostrans_dn(n,q) = - L_lostau * lostrans_dn(n)
                  L_factor1 = LC_Attenuations(n-1,q)*lostrans_dn(n) - LC_Attenuations(n,q)
                  L_factor2 = ( LC_suntau(n,q) - LC_suntau(n-1,q) ) / lostau
                  L_factor1 = L_factor1 + Attenuations(n-1) * LC_lostrans_dn(n,q)
!mick fix 9/19/2017 - replaced following line
                  !L_factor2 = L_factor2 - factor2 *  L_lostau/lostau
                  L_factor2 = L_factor2 - (factor2 + one) * L_lostau/lostau 
                  LC_multiplier(n,q) = ( L_factor1 - multiplier(n)*L_factor2 ) / factor2
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
          if ( do_sources_dn_p(ut,v) ) then
            np = Partial_layeridx(ut) ; kn = extinction(np)
            path_dn = FOGeometry%LosP_paths(ut,v)
            factor1 = zero ; factor2 = zero

 !  Sources, general case

            if ( FOGeometry%Mu1_dn(v) .gt. zero ) then
              lostau = kn * path_dn
              if ( lostau .lt. Expcutoff ) lostrans_dn_p(ut) = exp( - lostau )
              factor1  = Attenuations(np-1)*lostrans_dn_p(ut) - Attenuations_p(ut)
              factor2 = ((suntau_p(ut) - suntau(np-1))/lostau) - one
              multiplier_p(ut) = factor1 / factor2
              if ( do_columnwfs ) then
                do q = 1, n_columnwfs
                  L_lostau = L_extinction(np,q) * path_dn
                  LC_lostrans_dn_p(ut,q) = - L_lostau * lostrans_dn_p(ut)
                  L_factor1 = LC_Attenuations(np-1,q)*lostrans_dn_p(ut) - LC_Attenuations_p(ut,q)
                  L_factor2 = (LC_suntau_p(ut,q) - LC_suntau(np-1,q))/lostau
                  L_factor1 = L_factor1 + Attenuations(np-1)*LC_lostrans_dn_p(ut,q)
                  L_factor2 = L_factor2 - (factor2 + one) * L_lostau/lostau 
                  L_lostau = L_extinction(np,q) * path_dn
                  LC_multiplier_p(ut,q) = ( L_factor1 - multiplier_p(ut)*L_factor2 ) / factor2
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
          if ( layermask_dn(n) .and. do_sources_dn(n,v)  ) then
!mick fix 3/22/2017 - replaced index "np" with "n" in "LosW_paths"
            kn = extinction(n) ; path_dn = FOGeometry%LosW_paths(n,v)
            lostau = kn * path_dn ; if( lostau.lt.Expcutoff ) lostrans_dn(n) = exp ( - lostau )
            if ( do_columnwfs ) then
              do q = 1, Qnums(n)
                L_lostau        = L_extinction(n,q) * path_dn
                LC_lostrans_dn(n,q) = - L_lostau * lostrans_dn(n)
              enddo
            endif
            sum = zero ; L_sum = zero 
            do j = 1, FOGeometry%nfinedivs(n,v)
              dj = FOGeometry%LosW_paths(n,v) - FOGeometry%xfine(n,j,v) ; tran = exp ( - kn * dj )
              sum  = sum + attenuationsfine(n,j) * tran * FOGeometry%wfine(n,j,v)
              if ( do_columnwfs ) then
                do q = 1, n_columnwfs
                  L_tran = - dj * L_extinction(n,q)
                  L_func = ( LC_attenuationsfine(n,j,q) + L_tran * attenuationsfine(n,j) ) * tran
                  L_sum(q) = L_sum(q) + L_func * FOGeometry%wfine(n,j,v)
                enddo
              endif
            enddo
            multiplier(n) = sum * kn
            if ( do_columnwfs ) then
              do q = 1, n_columnwfs
                LC_multiplier(n,q) = kn * L_sum(q) + sum * L_extinction(n,q)
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
          if ( do_sources_dn_p(ut,v) ) then
            np = partial_layeridx(ut) ; kn = extinction(np)
            path_dn = FOGeometry%LosP_paths(ut,v)
            lostau = kn * path_dn ; if ( lostau.lt.Expcutoff ) lostrans_dn_p(ut) = exp ( - lostau )
            if ( do_columnwfs ) then
              do q = 1, n_columnwfs
                L_lostau        = L_extinction(np,q) * path_dn
                LC_lostrans_dn_p(ut,q) = - L_lostau * lostrans_dn_p(ut)
              enddo
            endif
            sum = zero ; L_sum = zero 
            do j = 1, FOGeometry%nfinedivs_p_dn(ut,v)
              dj = path_dn - FOGeometry%xfine_p_dn(ut,j,v) ; tran = exp ( - kn * dj )
              sum  = sum + attenuationsfine_p(ut,j) * tran * FOGeometry%wfine_p_dn(ut,j,v)
              if ( do_columnwfs ) then
                do q = 1, n_columnwfs
                  L_tran = - dj * L_extinction(np,q)
                  L_func = ( LC_attenuationsfine_p(ut,j,q) + L_tran * attenuationsfine_p(ut,j) ) * tran
                  L_sum(q) = L_sum(q) + L_func * FOGeometry%wfine_p_dn(ut,j,v)
                enddo
              endif
            enddo
            multiplier_p(ut) = sum * kn
            if ( do_columnwfs ) then
              do q = 1, n_columnwfs
                LC_multiplier_p(ut,q) = kn * L_sum(q) + sum * L_extinction(np,q)
              enddo
            endif
          endif
        enddo
      endif

!  Layer sources
!  -------------

!  General case, Whole layers
!mick fix 9/19/2017 - defined "sources_dn" for all stokes vector elements

      do n = 1, nlayers
        if ( layermask_dn(n) .and.do_sources_dn(n,v) ) then
          if ( do_sunlight ) then
             shelp(1:ns) = exactscat_dn(n,1:ns,1) * fluxvec(1)
          else
            do o1 = 1, nstokes
              shelp(o1) = dot_product(exactscat_dn(n,o1,1:ns),fluxvec(1:ns))
            enddo
          endif
          !sources_dn(n,o1) = shelp(o1) * multiplier(n)
          sources_dn(n,1:ns) = shelp(1:ns) * multiplier(n)
          if ( do_columnwfs ) then
            do q = 1, n_columnwfs
              do o1 = 1, nstokes
                LC_sources_dn(n,o1,q) = shelp(o1) * LC_multiplier(n,q)
                if ( do_sunlight ) then
                  L_Shelp = L_exactscat_dn(n,o1,1,q)*fluxvec(1)
                else
                  L_Shelp = dot_product(L_exactscat_dn(n,o1,1:ns,q),fluxvec(1:ns))
                endif
                LC_sources_dn(n,o1,q) =  LC_sources_dn(n,o1,q) + L_Shelp * multiplier(n)
              enddo
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
            if ( do_columnwfs ) then
              do q = 1, n_columnwfs
                do o1 = 1, nstokes
                  LC_sources_dn_p(ut,o1,q) = shelp(o1) * LC_multiplier_p(ut,q)
                  if ( do_sunlight ) then
                    L_Shelp = L_exactscat_dn(np,o1,1,q)*fluxvec(1)
                  else
                    L_Shelp = dot_product(L_exactscat_dn(np,o1,1:ns,q),fluxvec(1:ns))
                  endif
                  LC_sources_dn_p(ut,o1,q) = LC_sources_dn_p(ut,o1,q) + L_Shelp * multiplier_p(ut)
                enddo
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
!     Check for updating the recursion. Rob Fix Partials 8/22/16.

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

!  Column Wfs (atmospheric term)

      if ( do_columnwfs ) then
         L_CUMSOURCE = zero
         NSTART = 1 ; NUT_PREV = NSTART - 1
         DO UTA = 1, N_USER_LEVELS
            NUT = USER_LEVELS(UTA)
            DO N = NSTART, NUT
               NC = N
               do q = 1, n_columnwfs
                  L_cumsource(1:ns,q) = LC_SOURCES_DN(N,1:ns,Q)         + &
                                LC_LOSTRANS_DN(N,Q) * CUMSOURCE_DN(NC-1,1:ns) + &
                                   LOSTRANS_DN(N)   * L_CUMSOURCE(1:ns,Q)
               enddo
            ENDDO

            IF ( Partial_OUTFLAG(UTA) ) THEN
              UT = Partial_OUTINDEX(UTA)
              do q = 1, n_columnwfs
                Term1(1:NS) = CUMSOURCE_DN(NC,1:ns) * LC_LOSTRANS_DN_P(UT,q) + &
                                 L_CUMSOURCE(1:ns,Q) * LOSTRANS_DN_p(UT)
                LC_JACOBIANS_DN(UTA,1:ns,V,Q) = FLUX * ( TERM1(1:NS) + LC_SOURCES_DN_p(UT,1:ns,Q) )
              enddo
            ELSE
              do q = 1, n_columnwfs
                LC_JACOBIANS_DN(UTA,1:ns,V,Q) = FLUX * L_CUMSOURCE(1:ns,Q)
              enddo
            ENDIF
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1 ;  NUT_PREV = NUT
         ENDDO
!stop 'at FO col stopping point'
      endif

!  End geometry Loop

   enddo

!  Finish

   return
end subroutine FO_Vector_SSRT_ILCS_DN

!

subroutine FO_Vector_SSRT_ILCS_UPDN &
        ( do_upwelling, do_dnwelling, do_sunlight, do_deltam_scaling, do_lambertian,          & ! Inputs (Flags)
          do_surface_leaving, do_water_leaving, do_Partials, do_PlanPar, do_enhanced_ps,      & ! Inputs (Flags)
          do_sources_up, do_sources_up_p, do_sources_dn, do_sources_dn_p,                     & ! Inputs (Flags/sources)
          do_columnwfs,  n_columnwfs, LvaryFmat,                                              & ! Inputs (Control Jacobian)
          do_surfacewfs, do_sleavewfs, n_reflecwfs, n_sleavewfs, n_surfacewfs,                & ! Inputs (Control Jacobian)
          nstokes, ngeoms, nlayers, n_user_levels, user_level_mask_up, user_level_mask_dn,    & ! Inputs (control)
          npartials, partial_outindex, partial_outflag, partial_layeridx, FOGeometry,         & ! Inputs (partial/Geometry)
          flux, fluxvec, extinction, deltaus, omega, truncfac, fmatrix_up, fmatrix_dn,        & ! Inputs (Optical)
          L_extinction, L_deltaus, L_omega, L_truncfac, L_fmatrix_up, L_fmatrix_dn,           & ! Inputs (Linearized)
          reflec, slterm, LS_reflec, LSSL_slterm,                                             & ! Inputs (Surface)
          Stokes_up, Stokes_db, LC_Jacobians_up, LC_Jacobians_db, LS_Jacobians_db,            & ! output
          Stokes_dn, LC_Jacobians_dn, cumtrans, LC_cumtrans )                                   ! Output

!  FO routine for Upwelling and Downwelling Solar-beam Single-scatter (SS)
!    computation of Stokes-vectors and Column/Surface Jacobians. Inputs: geometry, spherical functions, optical properties.

   Implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ======

!  flags.
!  Version 1.5: --> F-matrix flag added 7/7/16; surface-leaving flag 8/2/16; Partials 8/22/16
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

   LOGICAL, Intent(in) :: do_columnwfs
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
   INTEGER, Intent(in) :: n_columnwfs
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
!  mick fix 3/2/2020 - moved dimension upper bound "max_atmoswfs" from dim 4 to dim 1 in L_FMATRIX_UP and L_FMATRIX_DN
  
   real(ffp), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_OMEGA       ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_TRUNCFAC    ( MAXLAYERS, max_atmoswfs )

   real(ffp), Intent(in) :: L_FMATRIX_UP  ( max_atmoswfs, MAXLAYERS, MAX_GEOMETRIES, 6 )
   real(ffp), Intent(in) :: L_FMATRIX_DN  ( max_atmoswfs, MAXLAYERS, MAX_GEOMETRIES, 6 )

!  Existence flags. 8/19/16. Criticality enters here

   logical, Intent(in)    :: do_sources_up       (maxlayers,MAX_GEOMETRIES)
   logical, Intent(in)    :: do_sources_up_p     (MAX_PARTLAYERS,MAX_GEOMETRIES)
   logical, Intent(in)    :: do_sources_dn       (maxlayers,MAX_GEOMETRIES)
   logical, Intent(in)    :: do_sources_dn_p     (MAX_PARTLAYERS,MAX_GEOMETRIES)

!  outputs
!  -------

   REAL(ffp), Intent(Out)  :: stokes_up     ( max_user_levels, 4, MAX_GEOMETRIES )
   REAL(ffp), Intent(Out)  :: stokes_db     ( max_user_levels, 4, MAX_GEOMETRIES )
   real(ffp), Intent(Out)  :: LC_Jacobians_up  ( max_user_levels, 4, MAX_GEOMETRIES, max_atmoswfs )
   real(ffp), Intent(Out)  :: LC_Jacobians_db  ( max_user_levels, 4, MAX_GEOMETRIES, max_atmoswfs )
   real(ffp), Intent(Out)  :: LS_Jacobians_db  ( max_user_levels, 4, MAX_GEOMETRIES, max_surfacewfs )

   REAL(ffp), Intent(Out)  :: stokes_dn     ( max_user_levels, 4, MAX_GEOMETRIES )
   real(ffp), Intent(Out)  :: LC_Jacobians_dn  ( max_user_levels, 4, MAX_GEOMETRIES, max_atmoswfs )

!  4/9/19. Additional output for the sleave correction

   real(ffp), Intent(out)  :: CUMTRANS    ( max_user_levels, MAX_GEOMETRIES )
   real(ffp), Intent(out)  :: LC_CUMTRANS ( max_user_levels, MAX_GEOMETRIES, max_atmoswfs )
   
!  upwelling
!  4/9/19. Additional output for the sleave correction

   if ( do_upwelling  ) then
      call FO_Vector_SSRT_ILCS_UP &
        ( do_sunlight, do_deltam_scaling, do_Lambertian, do_surface_leaving, do_water_leaving,    & ! Inputs (Flags-General/Surface)
          do_Partials, do_PlanPar, do_enhanced_ps, do_sources_up, do_sources_up_p,                & ! Inputs(Flags/criticality)
          do_columnwfs, do_surfacewfs, do_sleavewfs, n_columnwfs, LvaryFmat, n_reflecwfs,         & ! Inputs (Control, Jacobian)
          n_sleavewfs, n_surfacewfs, nstokes, ngeoms, nlayers, n_user_levels, user_level_mask_up, & ! Inputs (control, output)
          npartials, partial_outindex, partial_outflag, partial_layeridx, FOGeometry,             & ! Inputs (patial/Geometry)
          flux, fluxvec, extinction, deltaus, omega, truncfac, fmatrix_up, reflec, slterm,        & ! Inputs (Optical)
          L_extinction, L_deltaus, L_omega, L_truncfac, L_fmatrix_up, LS_reflec, LSSL_slterm,     & ! Inputs (Linearized)
          Stokes_up, Stokes_db, LC_Jacobians_up, LC_Jacobians_db, LS_Jacobians_db, cumtrans, LC_cumtrans )  ! Output
   endif

!  Downwelling

   if ( do_dnwelling  ) then
      call FO_Vector_SSRT_ILCS_DN &
      ( do_sunlight, do_deltam_scaling, do_Partials, do_PlanPar, do_enhanced_ps,    & ! Inputs (Flags/flux)
        do_sources_dn, do_sources_dn_p, do_columnwfs, n_columnwfs, LvaryFmat,       & ! Inputs (control, Jacobian )
        nstokes, ngeoms, nlayers, n_user_levels, user_level_mask_dn,                & ! Inputs (control)
        npartials, partial_outindex, partial_outflag, partial_layeridx, FOGeometry, & ! Inputs (control-partial)
        flux, fluxvec, extinction, deltaus, omega, truncfac, fmatrix_dn,            & ! Inputs (Optical)
        L_extinction, L_deltaus, L_omega, L_truncfac, L_fmatrix_dn,                 & ! Inputs (Optical - Linearized)
        Stokes_dn, LC_Jacobians_dn )                                                  ! Output
   endif

!  Finish

   return
end subroutine FO_Vector_SSRT_ILCS_UPDN

!  End module

end module FO_Vector_SSRT_ILCS_m

