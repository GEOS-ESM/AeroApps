
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

module FO_Thermal_DTRT_ILPS_m

!  FUNCTION
!  ========

!  For a given wavelength, this routine will calculate First-Order upwelling+downwelling 
!           Intensities(I), any number of LPS Jacobians (profile/surface

!  1. For the Atmospheric and Surface Direct Thermal Emission (DTE) sources.
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
!  Version     1.5, with optional phase function and partials.

!    Version 1a, 01 December 2011, R. Spurr, RT Solutions Inc.
!    Version 1b, 13 February 2012, R. Spurr, RT Solutions Inc.
!    Version 2,  01 June     2012, R. Spurr, RT Solutions Inc.
!    Version 3,  29 October  2012, Extension to Observational multiple geometries
!    Version 4,  31 July     2013, Lattice Multi-geometry
!    Version 5,  07 July     2016, Optional phase function usage
!    Version 5,  25 August   2016, Partial-layer output
!    Version 5,  11 December 2017, Optional code for polarized surface emissivity

!  For Thermal Emission sources, the subroutines are
!       FO_Thermal_DTRT_ILPS_UP   (Upwelling only)
!       FO_Thermal_DTRT_ILPS_DN   (Downwelling only)
!       FO_Thermal_DTRT_ILPS_UPDN (Upwelling and Downwelling)

!  Dependencies
!  ============

   use VLIDORT_PARS_m          , only : zero, one, Expcutoff, MAX_USER_LEVELS, MAXLAYERS, &
                                        MAXSTOKES, MAX_PARTLAYERS, maxfinelayers, MAX_USER_VZANGLES, &
                                        MAX_ATMOSWFS, MAX_SURFACEWFS
   use VLIDORT_Setups_def_m

public

contains

subroutine FO_Thermal_DTRT_ILPS_UP &
   ( do_deltam_scaling, do_Partials, do_PlanPar, do_enhanced_ps,        & ! Inputs (Flags)
     Do_Polarized_Emissivity, do_sources_up, do_sources_up_p,           & ! Inputs (Flags)
     do_profilewfs, do_surfacewfs, Lvaryflags, Lvarynums, n_surfacewfs, & ! Inputs (Control, Jacobians)
     nstokes, ngeoms, nlayers, n_user_levels, user_levels, npartials,   & ! Inputs (Control output)
     partial_outindex, partial_outflag, partial_layeridx, FOGeometry,   & ! Inputs (Partial/Geometry)
     extinction, deltaus, omega, truncfac, bb_input,                    & ! Inputs (Optical/thermal)
     surfbb, User_Emissivity, User_QUVEmissivity,                       & ! Inputs (Surface)
     L_extinction, L_deltaus, L_omega, L_truncfac,                      & ! Inputs (Optical - Linearized)
     LS_User_Emissivity, LS_User_QUVEmissivity,                         & ! Inputs (Surface - Linearized)
     intensity_dta_up, intensity_dts, LP_Jacobians_dta_up,              & ! Main Outputs
     LP_Jacobians_dts_up, LS_Jacobians_dts, tcom1, L_tcom1,             & ! Main Outputs
     lostrans_up, lostrans_up_p, L_lostrans_up, L_lostrans_up_p,        & ! Other Outputs
     StokesQUV_dts, LP_JacobiansQUV_dts_up, LS_JacobiansQUV_dts )         ! Optional Output. 12/11/17 Rob Add.

!  FO routine for Upwelling Direct-thermal-emission (DTE)
!    computation of Radiances and LPS Jacobians. Inputs: geometry, optical properties, Planck functions, emissivity

!  5/22/20. Version 2.8.2 Upgrades.
!    -  Add hfine/hfine_p inputs for correct DT calculation (Outgoing). These are in FOGeometry
!    -  lostrans_up, lostrans_up_p (and linearizations) are now outputs from this routine

   implicit none         

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ======

!  flags.  Version 1.5:  partials introduced, 8/26/16

   logical, Intent(in) ::  DO_DELTAM_SCALING
   logical, Intent(in) ::  DO_Partials
   logical, Intent(in) ::  DO_PLANPAR
   logical, Intent(in) ::  DO_ENHANCED_PS

!  Optional inputs for polarized emission. 12/11/17 Rob Add.

   logical, intent(in) :: Do_Polarized_Emissivity

!  Existence flags. 8/26/16. Criticality enters here

   logical, Intent(in) :: do_sources_up   (maxlayers,MAX_USER_VZANGLES)
   logical, Intent(in) :: do_sources_up_p (MAX_PARTLAYERS,MAX_USER_VZANGLES)

!  Jacobian Flags and control

   LOGICAL, Intent(in) :: do_surfacewfs
   LOGICAL, Intent(in) :: do_profilewfs
   LOGICAL, Intent(in) :: Lvaryflags(maxlayers)
   INTEGER, Intent(in) :: Lvarynums (maxlayers)
   INTEGER, Intent(in) :: n_surfacewfs

!  Numbers

   integer, Intent(in) ::  NSTOKES, NGEOMS, NLAYERS, N_USER_LEVELS
   integer, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  Numbers for Version 1.5: -->  Partial Control added, 8/26/16

   integer, Intent(in) :: Npartials
   integer, Intent(in) :: partial_layeridx (MAX_PARTLAYERS )
   logical, Intent(in) :: partial_outflag ( MAX_USER_LEVELS )
   integer, Intent(in) :: partial_outindex( MAX_USER_LEVELS )

!  Geometrical inputs
!  ------------------

   Type(VLIDORT_Geometry_FO), Intent(in) :: FOGeometry

!  optical inputs
!  --------------

!  Atmosphere optical

   real(ffp), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(ffp), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(ffp), Intent(in) :: OMEGA       ( MAXLAYERS )
   real(ffp), Intent(in) :: TRUNCFAC    ( MAXLAYERS )

!  Linearized optical inputs

   real(ffp), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_OMEGA       ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_TRUNCFAC    ( MAXLAYERS, max_atmoswfs )

!  Thermal/Surface BB and emissivity

   REAL(ffp), Intent(in) :: SURFBB, BB_INPUT    (0:MAXLAYERS)
   REAL(ffp), Intent(in) :: USER_EMISSIVITY(MAX_USER_VZANGLES),      LS_USER_EMISSIVITY    (MAX_USER_VZANGLES,max_surfacewfs)
   REAL(ffp), Intent(in) :: USER_QUVEmissivity(3,MAX_USER_VZANGLES), LS_USER_QUVEMISSIVITY (3,MAX_USER_VZANGLES,max_surfacewfs)

!  outputs
!  -------

!  Radiances

   real(ffp), Intent(Out)  :: intensity_dta_up ( max_user_levels, MAX_USER_VZANGLES )
   real(ffp), Intent(Out)  :: intensity_dts    ( max_user_levels, MAX_USER_VZANGLES )

!  Jacobians

   real(ffp), Intent(Out)  :: LP_Jacobians_dta_up  ( max_user_levels, MAX_USER_VZANGLES, maxlayers, max_atmoswfs )
   real(ffp), Intent(Out)  :: LP_Jacobians_dts_up  ( max_user_levels, MAX_USER_VZANGLES, maxlayers, max_atmoswfs )
   real(ffp), Intent(Out)  :: LS_Jacobians_dts     ( max_user_levels, MAX_USER_VZANGLES, max_surfacewfs )

!  5/22/20. Version 2.8.2 Upgrades.
!     ==> Add LOSTRANS arrays to output

   real(ffp), Intent(Out)  :: lostrans_up      ( maxlayers,      MAX_USER_VZANGLES )
   real(ffp), Intent(Out)  :: lostrans_up_p    ( MAX_PARTLAYERS, MAX_USER_VZANGLES )

   real(ffp), Intent(Out)  :: L_lostrans_up    ( maxlayers,      MAX_USER_VZANGLES, max_atmoswfs )
   real(ffp), Intent(Out)  :: L_lostrans_up_p  ( MAX_PARTLAYERS, MAX_USER_VZANGLES, max_atmoswfs )

!  Thermal setup

   real(ffp), Intent(InOut)   :: tcom1(maxlayers,2)
   real(ffp), Intent(InOut)   :: L_tcom1(maxlayers,2,max_atmoswfs)

!  Optional outputs for polarized emission. 12/11/17 Rob Add.

   real(ffp), Intent(Out) :: StokesQUV_dts          (max_user_levels,3,MAX_USER_VZANGLES)
   real(ffp), Intent(Out) :: LP_JacobiansQUV_dts_up (max_user_levels,3,MAX_USER_VZANGLES,maxlayers,max_atmoswfs)
   real(ffp), Intent(Out) :: LS_JacobiansQUV_dts    (max_user_levels,3,MAX_USER_VZANGLES,max_surfacewfs)

!  LOCAL
!  -----

!  Source function integration results. Partials added 8/26/16.

   real(ffp)  :: sources_up       ( maxlayers )
   real(ffp)  :: sources_up_p     ( MAX_PARTLAYERS )

   real(ffp)  :: L_sources_up     ( maxlayers,max_atmoswfs )
   real(ffp)  :: L_sources_up_p   ( MAX_PARTLAYERS,max_atmoswfs )

   real(ffp)  :: cumsource_up     ( 0:maxlayers )
   real(ffp)  :: L_cumsource      ( max_atmoswfs )
   real(ffp)  :: LS_cumsource     ( max_surfacewfs )

!  Regular_PS or plane-parallel flag

   logical    :: do_RegPSorPP

!  Help

   integer    :: n, j, k, q, v, uta, nstart, nc, nut, nut_prev, Qnums(maxlayers), np, ut
   logical    :: do_regular_ps, layermask_up(maxlayers), Qvary(maxlayers)


   real(ffp)  :: help, sum, tran, kn, xjkn, dj, zj, zjkn, zjL_kn, path_up, Solutionsfine, Solutionsfine_p
   real(ffp)  :: L_help, L_sum(max_atmoswfs), L_kn, L_Solutionsfine, L_Solutionsfine_p
   real(ffp)  :: t_mult_up(0:2), L_t_mult_up(0:2), thermcoeffs(2)
   real(ffp)  :: tms, L_tms(max_atmoswfs), lostau, LA2, partau, L_Partau
!mick mod 3/22/2017 - switched CUMSOURCE_DSTE from scalar to 1-D array 
   real(ffp)  :: cumsource_dste ( 0:maxlayers )

!  Help variables for Optional polarized emissivity. 12/11/17 Rob add.

   integer    :: ns, nws
   logical    :: do_PolEmiss
   real(ffp)  :: PolEmiss(3,MAX_USER_VZANGLES), LS_PolEmiss(3,MAX_USER_VZANGLES,max_surfacewfs)
   real(ffp)  :: QUV_dts(max_user_levels,3), CUMSOURCE_DSTEQUV( 0:maxlayers, 3) ! Followed Mick change to CUMSOURCE_DSTE
   real(ffp)  :: LS_QUV_dts(max_user_levels,3,max_surfacewfs), LS_CUMSOURCE_DSTEQUV(3,max_surfacewfs)
   real(ffp)  :: LP_QUV_dts(max_user_levels,3,max_atmoswfs),   LP_CUMSOURCE_QUV    (3,max_atmoswfs)

!  Zero the output

   INTENSITY_dta_up    = zero ; INTENSITY_dts       = zero
   LP_JACOBIANS_dta_up = zero ; LP_JACOBIANS_dts_up = zero
   LS_JACOBIANS_dts    = zero

!  Optional code for polarized emissivity, set Proxies and Initialize. 12/11/17 Rob add.
!mick fix 3/2/2020 - trimmed dimensions for defining PolEmiss & LS_PolEmiss

   do_PolEmiss = .false.
   if ( Do_Polarized_Emissivity.and. nstokes.gt.1 ) then
      do_PolEmiss = .true. ; ns = nstokes - 1
      PolEmiss(1:ns,1:ngeoms) = User_QUVEmissivity(1:ns,1:ngeoms) ; StokesQUV_dts = zero
      if ( do_surfacewfs ) then
         LS_PolEmiss(1:ns,1:ngeoms,1:n_surfacewfs) = LS_User_QUVEmissivity(1:ns,1:ngeoms,1:n_surfacewfs)
         LS_JacobiansQUV_dts = zero
      endif
      if ( do_profilewfs ) LP_JacobiansQUV_dts_up = zero
   endif

!  Regular_PS or plane-parallel flag

   do_regular_ps = .false.
   if ( .not.do_Planpar ) do_regular_ps = .not. do_enhanced_ps
   do_RegPSorPP = (do_regular_ps .or. do_PlanPar)

!  Bookkeeping
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

!  Thermal setup factors and linearizations
!     TMS, Initial set of thermal coefficients and TCOM1 variable

   tcom1 = zero ; L_tcom1 = zero
   do n = 1, nlayers
      tms = one - omega(n) ; L_tms = zero
      if ( Qvary(n) ) L_tms(1:Qnums(n)) = - L_omega(n,1:Qnums(n))
      if ( do_deltam_scaling ) then
         help = one - truncfac(n) * omega(n)
         tms = tms / help
         if ( Qvary(n) ) then
            do q = 1, Qnums(n)
               L_help = - L_truncfac(n,q)*omega(n) - truncfac(n) * L_omega(n,q)
               L_tms(q) = ( L_tms(q) - tms * L_help ) / help
            enddo
         endif
      endif
      thermcoeffs(1)  = bb_input(n-1)
      thermcoeffs(2)  = ( bb_input(n)-bb_input(n-1) ) / deltaus(n)
      tcom1(n,1) = thermcoeffs(1) * tms
      tcom1(n,2) = thermcoeffs(2) * tms
      if ( Qvary(n) ) then
         do q = 1, Qnums(n)
            LA2 = L_deltaus(n,q)/deltaus(n)
            L_tcom1(n,1,q) = thermcoeffs(1) * L_tms(q)
            L_tcom1(n,2,q) = thermcoeffs(2) * ( L_tms(q) - LA2  * tms )
         enddo
      endif
   ENDDO

!  5/22/20. Version 2.8.2 Upgrades.
!    ==>  Zero the transmittances

   lostrans_up   = zero ; L_lostrans_up   = zero 
   lostrans_up_p = zero ; L_lostrans_up_p = zero

!  Start Geometry loop
!  ===================

   do v = 1, ngeoms

!  Zero local sources
! mick mod 3/22/2017 - turned off (not needed). 4/24/20 Rob Query

      !sources_up = zero     ; cumsource_up = zero
      !sources_up_p = zero
      !L_sources_up   = zero
      !L_sources_up_p = zero

!  Plane/Parallel and Regular-PS: Layer integrated source terms
!  ============================================================

!  8/26/16. Version 1.5 partials introduced, nadir special case absorbed

!  5/22/20. Version 2.8.2 Upgrades.
!     ==> Add geometry index to Lostrans arrays (now outputs from subroutine)

      if ( do_RegPSorPP ) then
        DO n = 1, nlayers
          if ( layermask_up(n).and.do_sources_up(n,v) ) then
            lostau = deltaus(n) / FOGeometry%Mu1_LOS(v)
            if ( lostau .lt. Expcutoff ) lostrans_up(n,v) = exp( - lostau )
            if ( Qvary(n) ) L_lostrans_up(n,v,1:Qnums(n)) = - lostrans_up(n,v) * L_deltaus(n,1:Qnums(n)) / FOGeometry%Mu1_LOS(v)
            t_mult_up(2) = tcom1(n,2)
            t_mult_up(1) = tcom1(n,1) + t_mult_up(2) * FOGeometry%Mu1_LOS(v)
            sum = t_mult_up(1) + t_mult_up(2) * deltaus(n)
            t_mult_up(0) = - sum
            sources_up(n) = t_mult_up(0) * lostrans_up(n,v) + t_mult_up(1)
            if ( Qvary(n) ) then
              do q = 1, Qnums(n)
                L_t_mult_up(2) = L_tcom1(n,2,q)
                L_t_mult_up(1) = L_tcom1(n,1,q) + L_t_mult_up(2) * FOGeometry%Mu1_LOS(v)
                L_sum(q) = L_t_mult_up(1) + t_mult_up(2) * L_deltaus(n,q) + L_t_mult_up(2) * deltaus(n)
                L_t_mult_up(0) = - L_sum(q)
                L_sources_up(n,q) = L_t_mult_up(0) *   lostrans_up(n,v)   + &
                                      t_mult_up(0) * L_lostrans_up(n,v,q) + L_t_mult_up(1)
              enddo
            endif
          endif
        enddo
      endif

!  Partials. New Code 8/26/16
!mick fix 3/22/2017 - added "do_Partials" to 1st if condition
!                   - moved defining of "np" before 2nd if block
!                   - replaced last three lines with a new set of four in both intensity and
!                     jacobian sections
!  5/22/20. Version 2.8.2 Upgrades.
!     ==> Add geometry index to Lostrans arrays (now outputs from subroutine)

      if ( do_RegPSorPP.and.do_Partials ) then
        DO ut = 1, npartials
          np = partial_layeridx(ut)
          if ( layermask_up(np).and.do_sources_up_p(ut,v) ) then
            !np = partial_layeridx(ut)
            path_up = FOGeometry%losW_paths_LOS(np,v) - FOGeometry%losP_paths_LOS(ut,v) ; kn = extinction(np)
            lostau = kn * path_up ; if ( lostau .lt. Expcutoff ) lostrans_up_p(ut,v) = exp( - lostau )
            partau = lostau * FOGeometry%Mu1_LOS(v)
            if ( Qvary(np) ) L_lostrans_up_p(ut,v,1:Qnums(np)) = - lostrans_up_p(ut,v) * path_up * L_extinction(np,1:Qnums(np))
            t_mult_up(2) = tcom1(np,2)
            t_mult_up(1) = tcom1(np,1) + t_mult_up(2) * FOGeometry%Mu1_LOS(v)

            !sum = t_mult_up(1) + t_mult_up(2) * partau
            !t_mult_up(0) = - sum
            !sources_up_p(ut) = t_mult_up(0) * lostrans_up_p(ut,v) + t_mult_up(1)
            sum = t_mult_up(1) + t_mult_up(2) * deltaus(np)
            t_mult_up(0) = - sum
            sum = t_mult_up(1) + t_mult_up(2) * partau
            sources_up_p(ut) = t_mult_up(0) * lostrans_up_p(ut,v) + sum

            if ( Qvary(np) ) then
              do q = 1, Qnums(np)
                L_partau = L_extinction(np,q) * path_up * FOGeometry%Mu1_LOS(v)
                L_t_mult_up(2) = L_tcom1(np,2,q)
                L_t_mult_up(1) = L_tcom1(np,1,q) + L_t_mult_up(2) * FOGeometry%Mu1_LOS(v)

                !L_sum(q) = L_t_mult_up(1) + t_mult_up(2) * L_partau + L_t_mult_up(2) * partau
                !L_t_mult_up(0) = - L_sum(q)
                !L_sources_up_p(ut,q) = L_t_mult_up(0) *   lostrans_up_p(ut,v)   + &
                !                         t_mult_up(0) * L_lostrans_up_p(ut,q) + L_t_mult_up(1)
                L_sum(q) = L_t_mult_up(1) + t_mult_up(2) * L_deltaus(np,q) + L_t_mult_up(2) * deltaus(np)
                L_t_mult_up(0) = - L_sum(q)
                L_sum(q) = L_t_mult_up(1) + t_mult_up(2) * L_partau + L_t_mult_up(2) * partau
                L_sources_up_p(ut,q) = L_t_mult_up(0) *   lostrans_up_p(ut,v)   + &
                                         t_mult_up(0) * L_lostrans_up_p(ut,v,q) + L_sum(q)
              enddo
            endif
          endif
        enddo
      endif

!  LOS-spherical Layer integrated source terms
!  ===========================================

!  5/22/20. Version 2.8.2 Upgrades.
!     ==> Add geometry index to Lostrans arrays (now outputs from subroutine)
!     ==> Must use vertical distances in Thermal source terms (not path distances), Use zj instead of dj.

      if ( do_enhanced_ps ) then
        do n = nlayers, 1, -1
          if ( layermask_up(n) .and. do_sources_up(n,v) ) then
!mick fix 3/22/2017 - replaced index "np" with "n" in "FOGeometry%losW_paths_LOS"
            kn = extinction(n) ; path_up = FOGeometry%losW_paths_LOS(n,v)
            lostau = kn * path_up ; if( lostau.lt.Expcutoff ) lostrans_up_p(ut,v) = exp ( - lostau )
!mick fix 3/22/2017 - replaced "q" with "1:Qnums(n)" in "L_extinction"
            if ( Qvary(n) ) L_lostrans_up(n,v,1:Qnums(n)) = - lostrans_up(n,v) * path_up * L_extinction(n,1:Qnums(n))
            sum = zero ; L_sum = zero
            do j = 1, FOGeometry%nfinedivs_LOS(n,v)
              dj = FOGeometry%losW_paths_LOS(n,v) - FOGeometry%xfine_LOS(n,j,v) ; xjkn = dj * kn ; tran = exp ( - xjkn )
              zj = FOGeometry%hfine_LOS_up(n,j,v) ; zjkn = zj * kn ; solutionsfine = tcom1(n,1) + zjkn * tcom1(n,2)
              sum  = sum + solutionsfine * tran * FOGeometry%wfine_LOS(n,j,v)
              if ( Qvary(n) ) then
                do q = 1, Qnums(n)
                  L_kn = L_extinction(n,q)  ; zjL_kn = zj * L_kn
                  L_solutionsfine = L_tcom1(n,1,q) + zjkn * L_tcom1(n,2,q) + zjL_kn*tcom1(n,2)
                  L_sum(q)  = L_sum(q) + tran * FOGeometry%wfine_LOS(n,j,v) * &
                                             ( L_solutionsfine - dj * L_kn * solutionsfine )
                enddo
              endif
            enddo
            sources_up(n) = sum * kn
            if ( Qvary(n) ) then
              L_sources_up(n,1:Qnums(n)) = sum * L_extinction(n,1:Qnums(n)) + L_sum(1:Qnums(n)) * kn
            endif       
          endif
        enddo
      endif

!  Partials. 8/26/16.
!  5/22/20. Version 2.8.2 Upgrades.
!     ==> Add geometry index to Lostrans arrays (now outputs from subroutine)
!     ==> Must use vertical distances in Thermal source terms (not path distances), Use zj instead of dj.

      if ( do_enhanced_ps.and.do_Partials ) then
        do ut = 1, npartials
          if ( do_sources_up_p(ut,v) ) then
            np = partial_layeridx(ut) ; kn = extinction(np)
            path_up = FOGeometry%losW_paths_LOS(np,v)- FOGeometry%losP_paths_LOS(ut,v)
            lostau = kn * path_up ; if ( lostau.lt.Expcutoff ) lostrans_up_p(ut,v) = exp ( - lostau )
            if ( Qvary(np) ) L_lostrans_up_p(ut,v,1:Qnums(np)) = - lostrans_up_p(ut,v) * path_up * L_extinction(np,1:Qnums(np))
            sum = zero ; L_sum = zero
            do j = 1, FOGeometry%nfinedivs_p_LOS_up(ut,v)
              dj = path_up - FOGeometry%xfine_p_LOS_up(ut,j,v) ; xjkn = dj * kn ; tran = exp ( - xjkn )     ! Correct
              zj = FOGeometry%hfine_p_LOS_up(ut,j,v) ; zjkn = zj * kn ; solutionsfine_p = tcom1(np,1) + zjkn * tcom1(np,2)
              sum  = sum + solutionsfine_p * tran * FOGeometry%wfine_p_LOS_up(ut,j,v)
              if ( Qvary(np) ) then
                do q = 1, Qnums(np)
                  L_kn = L_extinction(np,q) ; zjL_kn = zj * L_kn
                  L_solutionsfine_p = L_tcom1(np,1,q) + zjkn * L_tcom1(np,2,q) + zjL_kn*tcom1(np,2)
                  L_sum(q)  = L_sum(q) + tran * FOGeometry%wfine_p_LOS_up(ut,j,v) * &
                               ( L_solutionsfine_p - zjL_kn* solutionsfine_p )
                enddo
              endif
            enddo
            sources_up_p(ut) = sum * kn
            if ( Qvary(np) ) then
              L_sources_up_p(ut,1:Qnums(np)) = sum * L_extinction(np,1:Qnums(np)) + L_sum(1:Qnums(np)) * kn
            endif       
          endif
        enddo        
      endif

!  Source function integration
!  ===========================

!  start recursion ( For DSTE term, Use surface emissivity )
!mick mod 3/22/2017 - switched CUMSOURCE_DSTE from scalar to 1-D array
!  5/22/20. Version 2.8.2 Upgrades. ==> Must add geometry index to Lostrans

      NC =  0
      CUMSOURCE_UP(NC)   = zero
      CUMSOURCE_DSTE(NC) = SURFBB * USER_EMISSIVITY(V)
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  Main Intensity loop over all output optical depths
!     NLEVEL = Layer index for given optical depth
!     Cumulative source terms : Loop over layers working upwards from NSTART to level NUT,
!     Check for updating the recursion. Robfix partials, 8/26/16.
!  5/22/20. Version 2.8.2 Upgrades. ==> Must add geometry index to Lostrans

      DO UTA = N_USER_LEVELS, 1, -1
         NUT = USER_LEVELS(UTA) + 1
         DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            CUMSOURCE_UP(NC)   = lostrans_up(n,v) * CUMSOURCE_UP(NC-1) + SOURCES_UP(N)
            CUMSOURCE_DSTE(NC) = lostrans_up(n,v) * CUMSOURCE_DSTE(NC-1)
         ENDDO
         IF ( Partial_OUTFLAG(UTA) ) THEN
           UT = Partial_OUTINDEX(UTA)
           INTENSITY_DTA_UP(UTA,V) = CUMSOURCE_UP(NC) * lostrans_up_p(ut,v) + SOURCES_UP_p(UT)
           INTENSITY_DTS(UTA,V)    = CUMSOURCE_DSTE(NC) * lostrans_up_p(ut,v)
         ELSE
           INTENSITY_DTA_UP(UTA,V) = CUMSOURCE_UP(NC)
           INTENSITY_DTS(UTA,V)    = CUMSOURCE_DSTE(NC)
         ENDIF
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
         NUT_PREV = NUT
      ENDDO

!  Optional code for polarized emissivity. 12/11/17 Rob add.
!  5/22/20. Version 2.8.2 Upgrades. ==> Must add geometry index to Lostrans

      if ( do_PolEmiss ) then
         NC = 0 ; CUMSOURCE_DSTEQUV(NC,1:ns) = SURFBB * PolEmiss(1:ns,v)
         NSTART = NLAYERS ; NUT_PREV = NSTART + 1
         DO UTA = N_USER_LEVELS, 1, -1
            NUT = USER_LEVELS(UTA) + 1
            DO N = NSTART, NUT, -1
              CUMSOURCE_DSTEQUV(NC,1:ns) = lostrans_up(n,v) * CUMSOURCE_DSTEQUV(NC-1,1:ns)
            ENDDO
            IF ( Partial_OUTFLAG(UTA) ) THEN
              UT = Partial_OUTINDEX(UTA)
              QUV_DTS(UTA,1:ns) = CUMSOURCE_DSTEQUV(NC,1:ns)  * lostrans_up_p(ut,v)
            ELSE
              QUV_DTS(UTA,1:ns) = CUMSOURCE_DSTEQUV(NC,1:ns)
            ENDIF
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1; NUT_PREV = NUT
         ENDDO
         StokesQUV_dts(1:N_USER_LEVELS,1:ns,v) = QUV_DTS(1:N_USER_LEVELS,1:ns)
      endif

!  Surface WFs. Partials added, 8/26/16
!mick mod 3/22/2017 - turned off NC (not needed here)
!  5/22/20. Version 2.8.2 Upgrades. ==> Must add geometry index to Lostrans

      if ( do_surfacewfs ) then
        do q = 1, n_surfacewfs
          LS_cumsource(q) = SURFBB * LS_USER_EMISSIVITY(V,q)
        enddo
        NSTART = NLAYERS
        NUT_PREV = NSTART + 1
        DO UTA = N_USER_LEVELS, 1, -1
          NUT    = USER_LEVELS(UTA) + 1
          DO N = NSTART, NUT, -1
            !NC = NLAYERS + 1 - N
            do q = 1, n_surfacewfs
              LS_cumsource(q) = lostrans_up(n,v) * LS_CUMSOURCE(Q)
            enddo
          ENDDO
          IF ( Partial_OUTFLAG(UTA) ) THEN
            UT = Partial_OUTINDEX(UTA)
            LS_Jacobians_DTS(UTA,V,1:n_surfacewfs) = LS_cumsource(1:n_surfacewfs) * lostrans_up_p(ut,v)
          ELSE
            LS_Jacobians_DTS(UTA,V,1:n_surfacewfs) = LS_cumsource(1:n_surfacewfs)
          ENDIF
          IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT
        ENDDO
      endif

!  Optional code for polarized emissivity (Surface Jacobians).12/11/17 Rob add.
!  5/22/20. Version 2.8.2 Upgrades. ==> Must add geometry index to Lostrans

      if ( do_SurfaceWfs .and. do_PolEmiss ) then
         nws = n_SurfaceWfs
         LS_CUMSOURCE_DSTEQUV(1:ns,1:nws) = SURFBB * LS_PolEmiss(1:ns,1:nws,v)
         NSTART = NLAYERS ; NUT_PREV = NSTART + 1
         DO UTA = N_USER_LEVELS, 1, -1
            NUT = USER_LEVELS(UTA) + 1
            DO N = NSTART, NUT, -1
               LS_CUMSOURCE_DSTEQUV(1:ns,1:nws) = lostrans_up(n,v) * LS_CUMSOURCE_DSTEQUV(1:ns,1:nws)
            ENDDO
            IF ( Partial_OUTFLAG(UTA) ) THEN
              UT = Partial_OUTINDEX(UTA)
              LS_QUV_DTS(UTA,1:ns,1:nws) = LS_CUMSOURCE_DSTEQUV(1:ns,1:nws)  * lostrans_up_p(ut,v)
            ELSE
              LS_QUV_DTS(UTA,1:ns,1:nws) = LS_CUMSOURCE_DSTEQUV(1:ns,1:nws)
            ENDIF
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1; NUT_PREV = NUT
         ENDDO
         LS_JacobiansQUV_dts(1:N_USER_LEVELS,1:ns,1:nws,v) = LS_QUV_DTS(1:N_USER_LEVELS,1:ns,1:nws)
      endif

!  Profile Wfs (atmospheric term)
!  5/22/20. Version 2.8.2 Upgrades. ==> Must add geometry index to Lostrans

      if ( do_profilewfs ) then
        do k = 1, nlayers
          if ( Qvary(k) ) then
!mick fix 3/22/2017 - initialize NC
            NC = 0
            L_CUMSOURCE = zero
            NSTART = NLAYERS
            NUT_PREV = NSTART + 1
            DO UTA = N_USER_LEVELS, 1, -1
              NUT    = USER_LEVELS(UTA) + 1
              DO N = NSTART, NUT, -1
                NC = NLAYERS + 1 - N
                if ( k.eq.n ) then
!mick mod 3/22/2017 - turned off 3rd term (not needed)
                  do q = 1, Qnums(k)
                    L_cumsource(q) = L_SOURCES_UP(N,Q) &
                     + L_LOSTRANS_UP(N,v,Q) * CUMSOURCE_UP(NC-1) !+ lostrans_up(n,v) * L_CUMSOURCE(Q)
                  enddo
                else
                  do q = 1, Qnums(k)
                    L_cumsource(q) = lostrans_up(n,v) * L_CUMSOURCE(Q)
                  enddo
                endif
              ENDDO
              IF ( Partial_OUTFLAG(UTA) ) THEN
                UT = Partial_OUTINDEX(UTA) ; np = partial_layeridx(ut)
                if ( k.eq.np ) then
!mick fix 3/22/2017 - added L_SOURCES_UP_P term here
!mick mod 3/22/2017 - turned off 3rd term (not needed)
                  do q = 1, Qnums(k)
                    LP_Jacobians_DTA_UP(UTA,V,K,q) = L_SOURCES_UP_P(UT,q) &
                      + L_LOSTRANS_UP_P(UT,v,Q) * CUMSOURCE_UP(NC) !+ lostrans_up_p(ut,v) * L_CUMSOURCE(Q)
                  enddo
                else
                  do q = 1, Qnums(k)
                    LP_Jacobians_DTA_UP(UTA,V,K,q) = lostrans_up_p(ut,v) * L_CUMSOURCE(q) 
                  enddo
                endif
              ELSE
                do q = 1, Qnums(k)
                  LP_Jacobians_dta_up(UTA,V,K,Q) = L_CUMSOURCE(Q)
                enddo
              endif
              IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
              NUT_PREV = NUT
            ENDDO
          endif
        enddo
      endif

!  Profile Wfs (surface emission term)
!  5/22/20. Version 2.8.2 Upgrades. ==> Must add geometry index to Lostrans

      if ( do_profilewfs ) then
!mick mod 3/22/2017 - defining of CUMSOURCE_DSTE turned off here since CUMSOURCE_DSTE
!                     now an array defined above
        !CUMSOURCE_DSTE = SURFBB * USER_EMISSIVITY(v)
        do k = 1, nlayers
          if ( Qvary(k) ) then
!mick fix 3/22/2017 - initialize NC
!                   - initialize L_CUMSOURCE to zero
            NC = 0
            !L_CUMSOURCE = CUMSOURCE_DSTE
            L_CUMSOURCE = zero
            NSTART = NLAYERS
            NUT_PREV = NSTART + 1
            DO UTA = N_USER_LEVELS, 1, -1
              NUT  = USER_LEVELS(UTA) + 1
              DO N = NSTART, NUT, -1
                NC = NLAYERS + 1 - N
                if ( k.eq.n ) then
!mick fix 3/22/2017 - use CUMSOURCE_DSTE array
                  do q = 1, Qnums(k)
                    !L_cumsource(q) = L_LOSTRANS_UP(N,v,Q) * L_CUMSOURCE(Q)
                    L_cumsource(q) = L_LOSTRANS_UP(N,v,Q) * CUMSOURCE_DSTE(NC-1)
                  enddo
                else
                  do q = 1, Qnums(k)
                    L_cumsource(q) = lostrans_up(n,v) * L_CUMSOURCE(Q)
                  enddo
                endif
              ENDDO
              IF ( Partial_OUTFLAG(UTA) ) THEN
                UT = Partial_OUTINDEX(UTA) ; np = partial_layeridx(ut)
                if ( k.eq.np ) then
!mick fix 3/22/2017 - use CUMSOURCE_DSTE array
                  do q = 1, Qnums(k)
                    !LP_Jacobians_DTS_UP(UTA,V,K,q) = L_LOSTRANS_UP_p(UT,v,Q) * L_CUMSOURCE(Q)
                    LP_Jacobians_DTS_UP(UTA,V,K,q) = L_LOSTRANS_UP_p(UT,v,Q) * CUMSOURCE_DSTE(NC)
                  enddo
                else
                  do q = 1, Qnums(k)
                    LP_Jacobians_DTS_UP(UTA,V,K,q) = lostrans_up_p(ut,v) * L_CUMSOURCE(Q) 
                  enddo
                endif
              ELSE
                do q = 1, Qnums(k)
                  LP_Jacobians_dts_up(UTA,V,K,Q) = L_CUMSOURCE(Q)
                enddo
              endif
              IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
              NUT_PREV = NUT
            ENDDO
          endif
        enddo
      endif

!  Optional code for polarized emissivity (Profile Jacobians).12/11/17 Rob add.
!  5/22/20. Version 2.8.2 Upgrades. ==> Must add geometry index to Lostrans

      if ( do_profilewfs .and. do_PolEmiss ) then
        do k = 1, nlayers
          if ( Qvary(k) ) then
            NC = 0
            LP_CUMSOURCE_QUV = zero
            NSTART = NLAYERS
            NUT_PREV = NSTART + 1
            DO UTA = N_USER_LEVELS, 1, -1
              NUT  = USER_LEVELS(UTA) + 1
              DO N = NSTART, NUT, -1
                NC = NLAYERS + 1 - N
                if ( k.eq.n ) then
                  do q = 1, Qnums(k)
                    LP_CUMSOURCE_QUV(1:ns,q) = L_LOSTRANS_UP(N,v,Q) * CUMSOURCE_DSTEQUV(NC-1,1:ns)
                  enddo
                else
                  do q = 1, Qnums(k)
                    LP_CUMSOURCE_QUV(1:ns,q) = lostrans_up(n,v) * LP_CUMSOURCE_QUV(1:ns,q)
                  enddo
                endif
              ENDDO
              IF ( Partial_OUTFLAG(UTA) ) THEN
                UT = Partial_OUTINDEX(UTA) ; np = partial_layeridx(ut)
                if ( k.eq.np ) then
                  do q = 1, Qnums(k)
                    LP_QUV_DTS(UTA,1:ns,q) = L_LOSTRANS_UP_p(UT,v,Q) * CUMSOURCE_DSTEQUV(NC,1:ns)
                  enddo
                else
                  do q = 1, Qnums(k)
                    LP_QUV_DTS(UTA,1:ns,q) = lostrans_up_p(ut,v) * LP_CUMSOURCE_QUV(1:ns,q) 
                  enddo
                endif
              ELSE
                do q = 1, Qnums(k)
                  LP_QUV_DTS(UTA,1:ns,Q) = LP_CUMSOURCE_QUV(1:ns,q) 
                enddo
              endif
              IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
              NUT_PREV = NUT
            ENDDO
            do q = 1, Qnums(k)
              LP_JacobiansQUV_dts_up(1:N_USER_LEVELS,1:ns,v,k,q) = LP_QUV_DTS(1:N_USER_LEVELS,1:ns,Q)
            enddo
          endif
        enddo
      endif

!  End geometry loop

   enddo

!  Finish

   return
end subroutine FO_Thermal_DTRT_ILPS_UP

!

subroutine FO_Thermal_DTRT_ILPS_DN &
   ( do_deltam_scaling, do_Partials, do_PlanPar, do_enhanced_ps,               & ! Inputs (Flags)
     do_sources_dn, do_sources_dn_p, do_profilewfs, Lvaryflags, Lvarynums,     & ! Inputs (Flags/Jac-control)
     ngeoms, nlayers, n_user_levels, user_levels, npartials,                   & ! Inputs (control output)
     partial_outindex, partial_outflag, partial_layeridx, FOGeometry,          & ! Inputs (control-partial)
     bb_input, extinction, deltaus, omega, truncfac,                           & ! Inputs (Optical - Regular)
     L_extinction, L_deltaus, L_omega, L_truncfac,                             & ! Inputs (Optical - Linearized)
     intensity_dta_dn, LP_Jacobians_dta_dn, tcom1, L_tcom1 )                     ! Output

!  FO routine for Downwelling Direct-thermal-emission (DTE)
!    computation of Radiances and LPS Jacobians. Inputs: geometry, optical properties, Planck functions, emissivity

!  5/22/20. Version 2.8.2 Upgrades.
!    -  Add hfine/hfine_p inputs for correct DT calculation (Outgoing). These are in FOGEOMETRY

   implicit none         

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ======

!  flags. Version 1.5:  Partials 8/26/16

   logical, Intent(in) :: DO_DELTAM_SCALING

   logical, Intent(in) :: DO_Partials
   logical, Intent(in) :: DO_PLANPAR
   logical, Intent(in) :: DO_ENHANCED_PS

!  Existence flags. 8/26/16. Criticality enters here

   logical, Intent(in)    :: do_sources_dn       (maxlayers,MAX_USER_VZANGLES)
   logical, Intent(in)    :: do_sources_dn_p     (MAX_PARTLAYERS,MAX_USER_VZANGLES)

!  Jacobian Flag and control

   LOGICAL, Intent(in) :: do_profilewfs
   LOGICAL, Intent(in) :: Lvaryflags(maxlayers)
   INTEGER, Intent(in) :: Lvarynums (maxlayers)

!  Numbers

   integer, Intent(in) :: NLAYERS, NGEOMS, N_USER_LEVELS
   integer, Intent(in) :: USER_LEVELS ( MAX_USER_LEVELS )

!  Numbers for Version 1.5: -->  Partial Control added, 8/26/16

   integer, Intent(in) :: Npartials
   integer, Intent(in) :: partial_layeridx(MAX_PARTLAYERS)
   logical, Intent(in) :: partial_outflag ( MAX_USER_LEVELS )
   integer, Intent(in) :: partial_outindex( MAX_USER_LEVELS )

!  Geometrical inputs
!  ------------------

   Type(VLIDORT_Geometry_FO), Intent(in) :: FOGeometry

!  optical inputs
!  --------------

!  Atmosphere extinction and deltaus

   real(ffp), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(ffp), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(ffp), Intent(in) :: OMEGA       ( MAXLAYERS )
   real(ffp), Intent(in) :: TRUNCFAC    ( MAXLAYERS )

!  Linearized optical inputs

   real(ffp), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_OMEGA       ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_TRUNCFAC    ( MAXLAYERS, max_atmoswfs )

!  Atmospheric BB functions

   REAL(ffp), Intent(in) :: BB_INPUT (0:MAXLAYERS)

!  outputs
!  -------

!  Radiances

   real(ffp), Intent(Out)  :: intensity_dta_dn ( max_user_levels, MAX_USER_VZANGLES )

!  Jacobians

   real(ffp), Intent(Out)  :: LP_Jacobians_dta_dn  ( max_user_levels, MAX_USER_VZANGLES, maxlayers, max_atmoswfs )

!  Thermal setup

   real(ffp), Intent(InOut)   :: tcom1(maxlayers,2)
   real(ffp), Intent(InOut)   :: L_tcom1(maxlayers,2,max_atmoswfs)

!  LOCAL
!  -----

!  Source function integration results. Partials added 8/26/16.

   real(ffp)  :: sources_dn       ( maxlayers )
   real(ffp)  :: lostrans_dn      ( maxlayers )
   real(ffp)  :: sources_dn_p     ( MAX_PARTLAYERS )
   real(ffp)  :: lostrans_dn_p    ( MAX_PARTLAYERS )

   real(ffp)  :: L_lostrans_dn    ( maxlayers,max_atmoswfs )
   real(ffp)  :: L_sources_dn     ( maxlayers,max_atmoswfs )
   real(ffp)  :: L_lostrans_dn_p  ( MAX_PARTLAYERS,max_atmoswfs )
   real(ffp)  :: L_sources_dn_p   ( MAX_PARTLAYERS,max_atmoswfs )

   real(ffp)  :: cumsource_dn     ( 0:maxlayers )
   real(ffp)  :: L_cumsource      ( max_atmoswfs )

!  Regular_PS or plane-parallel flag

   logical    :: do_RegPSorPP

!  Help

   integer    :: n, j, k, q, v, uta, nstart, nc, nut, nut_prev, Qnums(maxlayers), np, ut
   logical    :: do_regular_PS, layermask_dn(maxlayers), Qvary(maxlayers)

   real(ffp)  :: help, sum, tran, kn, xjkn, dj, zj, zjkn, zjL_kn, lostau, partau, path_dn
   real(ffp)  :: L_help, L_sum(max_atmoswfs), L_kn, L_partau
   real(ffp)  :: t_mult_dn(0:2), L_t_mult_dn(0:2), thermcoeffs(2)
   real(ffp)  :: tms, L_tms(max_atmoswfs), LA2, Solutionsfine, Solutionsfine_p, L_Solutionsfine, L_Solutionsfine_p

   real(ffp), parameter  :: Expcutoff = 88.0d0
   real(ffp), parameter  :: zero   = 0.0_ffp
   real(ffp), parameter  :: one    = 1.0_ffp

!  Zero the output

   INTENSITY_dta_dn    = zero
   LP_JACOBIANS_dta_dn = zero

!  Regular_PS or plane-parallel flag

   do_regular_ps = .false.
   if ( .not.do_Planpar ) do_regular_ps = .not. do_enhanced_ps
   do_RegPSorPP = (do_regular_ps .or. do_PlanPar)

!  Bookkeeping

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

!  Thermal setup factors and linearizations
!     TMS, Initial set of thermal coefficients and TCOM1 variable

   tcom1 = zero ; L_tcom1 = zero
   do n = 1, nlayers
      tms = one - omega(n) ; L_tms = zero
      if ( Qvary(n) ) L_tms(1:Qnums(n)) = - L_omega(n,1:Qnums(n))
      if ( do_deltam_scaling ) then
         help = one - truncfac(n) * omega(n)
         tms = tms / help
         if ( Qvary(n) ) then
            do q = 1, Qnums(n)
               L_help = - L_truncfac(n,q)*omega(n) - truncfac(n) * L_omega(n,q)
               L_tms(q) = ( L_tms(q) - tms * L_help ) / help
            enddo
         endif
      endif
      thermcoeffs(1)  = bb_input(n-1)
      thermcoeffs(2)  = (bb_input(n)-bb_input(n-1)) / deltaus(n)
      tcom1(n,1) = thermcoeffs(1) * tms
      tcom1(n,2) = thermcoeffs(2) * tms
      if ( Qvary(n) ) then
         do q = 1, Qnums(n)
            LA2 = L_deltaus(n,q)/deltaus(n)
            L_tcom1(n,1,q) = thermcoeffs(1) * L_tms(q)
            L_tcom1(n,2,q) = thermcoeffs(2) * ( L_tms(q) - LA2  * tms )
         enddo
      endif
   ENDDO

!  Start Geometry loop
!  ===================

   do v = 1, ngeoms

!  Zero local sources

      lostrans_dn   = zero    ; sources_dn   = zero   ; cumsource_dn = zero
      L_lostrans_dn = zero    ; L_sources_dn = zero

      lostrans_dn_p   = zero  ; sources_dn_p   = zero
      L_lostrans_dn_p = zero  ; L_sources_dn_p = zero

!  Plane/Parallel and Regular-PS: Layer integrated source terms
!  ============================================================

!  Bug Fixed 23 January 2013 (nadir case). Old code commented out and replaced
!  8/26/16. Version 1.5 partials introduced, nadir special case absorbed

      if ( do_RegPSorPP ) then
        DO n = 1, nlayers
          if ( layermask_dn(n).and.do_sources_dn(n,v) ) then
            lostau = deltaus(n) / FOGeometry%Mu1_LOS(v)
            if ( lostau .lt. Expcutoff ) lostrans_dn(n) = exp( - lostau )
            if ( Qvary(n) ) L_lostrans_dn(n,1:Qnums(n)) = - lostrans_dn(n) * L_deltaus(n,1:Qnums(n)) / FOGeometry%Mu1_LOS(v)
            t_mult_dn(2)   = tcom1(n,2)
            t_mult_dn(1)   = tcom1(n,1) - t_mult_dn(2) * FOGeometry%Mu1_LOS(v)
            t_mult_dn(0)   = - t_mult_dn(1)
            sources_dn(n)  = t_mult_dn(0) * lostrans_dn(n)
            sum = t_mult_dn(1) + t_mult_dn(2) * deltaus(n)
            sources_dn(n)  = sources_dn(n) + sum
            if ( Qvary(n) ) then
              do q = 1, Qnums(n)
                 L_t_mult_dn(2) = L_tcom1(n,2,q)
                 L_t_mult_dn(1) = L_tcom1(n,1,q) - L_t_mult_dn(2) * FOGeometry%Mu1_LOS(v)
                 L_t_mult_dn(0) = - L_t_mult_dn(1)
                 L_sources_dn(n,q)  = L_t_mult_dn(0) * lostrans_dn(n) + t_mult_dn(0) * L_lostrans_dn(n,q)
                 L_sum(q) = L_t_mult_dn(1) + t_mult_dn(2) * L_deltaus(n,q) + L_t_mult_dn(2) * deltaus(n)
                 L_sources_dn(n,q) = L_sources_dn(n,q) + L_sum(q)
              enddo
            endif
          endif
        enddo
      endif

!  Partials. New Code 8/26/16
!mick fix 3/22/2017 - added "do_Partials" to 1st if condition
!                   - replaced two lines

      if ( do_RegPSorPP.and.do_Partials ) then
        DO ut = 1, npartials
          np = partial_layeridx(ut)
          if ( layermask_dn(np).and.do_sources_dn_p(ut,v) ) then
            path_dn = FOGeometry%losP_paths_LOS(ut,v) ; kn = extinction(np)
            lostau = kn * path_dn ; if ( lostau .lt. Expcutoff ) lostrans_dn_p(ut) = exp( - lostau )
            if ( Qvary(np) ) L_lostrans_dn_p(ut,1:Qnums(np)) = - lostrans_dn_p(ut) * path_dn * L_extinction(np,1:Qnums(np))
            partau = lostau * FOGeometry%Mu1_LOS(v)
            t_mult_dn(2)   = tcom1(np,2)
            t_mult_dn(1)   = tcom1(np,1) - t_mult_dn(2) * FOGeometry%Mu1_LOS(v)
            t_mult_dn(0)   = - t_mult_dn(1)
            !sources_dn_p(n)  = t_mult_dn(0) * lostrans_dn_p(n)
            sources_dn_p(ut)  = t_mult_dn(0) * lostrans_dn_p(ut)
            sum = t_mult_dn(1) + t_mult_dn(2) * partau
            !sources_dn(n)  = sources_dn(n) + sum
            sources_dn_p(ut) = sources_dn_p(ut) + sum
            if ( Qvary(np) ) then
              do q = 1, Qnums(np)
                 L_partau = L_extinction(np,q) * path_dn * FOGeometry%Mu1_LOS(v)
                 L_t_mult_dn(2) = L_tcom1(np,2,q)
                 L_t_mult_dn(1) = L_tcom1(np,1,q) - L_t_mult_dn(2) * FOGeometry%Mu1_LOS(v)
                 L_t_mult_dn(0) = - L_t_mult_dn(1)
                 L_sources_dn_p(ut,q)  = L_t_mult_dn(0) * lostrans_dn_p(ut) + t_mult_dn(0) * L_lostrans_dn_p(ut,q)
                 L_sum(q) = L_t_mult_dn(1) + t_mult_dn(2) * L_partau + L_t_mult_dn(2) * partau
                 L_sources_dn_p(ut,q) = L_sources_dn_p(ut,q) + L_sum(q)
              enddo
            endif
          endif
        enddo
      endif

!  LOS-spherical Layer integrated source terms
!  ===========================================

!  5/22/20. Version 2.8.2 Upgrades.
!     ==> Must use vertical distances in Thermal source terms (not path distances). ZJ instead of DJ

      if ( do_enhanced_ps ) then
        do n = nlayers, 1, -1
          if ( layermask_dn(n) .and. do_sources_dn(n,v) ) then
!mick fix 3/22/2017 - replaced index "np" with "n" in "FOGeometry%losW_paths_LOS"
            kn = extinction(n) ; path_dn = FOGeometry%losW_paths_LOS(n,v)
            lostau = kn * path_dn ; if( lostau.lt.Expcutoff ) lostrans_dn(n) = exp ( - lostau )
!mick fix 3/22/2017 - replaced "q" with "1:Qnums(n)" in "L_extinction"
            if ( Qvary(n) ) L_lostrans_dn(n,1:Qnums(n)) = - lostrans_dn(n) * path_dn * L_extinction(n,1:Qnums(n))
            sum = zero ; L_sum = zero
            do j = 1, FOGeometry%nfinedivs_LOS(n,v)
              dj = FOGeometry%LosW_paths_LOS(n,v) - FOGeometry%xfine_LOS(n,j,v) ; xjkn = dj * kn ; tran = exp ( - xjkn )
              zj = FOGeometry%hfine_LOS_up(n,j,v) ; zjkn = zj * kn ; solutionsfine = tcom1(n,1) + zjkn * tcom1(n,2)
              sum  = sum + solutionsfine * tran * FOGeometry%wfine_LOS(n,j,v)
              if ( Qvary(n) ) then
                do q = 1, Qnums(n)
                  L_kn = L_extinction(n,q) ; zjL_kn = zj * L_kn
                  L_solutionsfine = L_tcom1(n,1,q) + zjkn * L_tcom1(n,2,q) + zjL_kn * tcom1(n,2)
                  L_sum(q)  = L_sum(q) + tran * FOGeometry%wfine_LOS(n,j,v) * &
                                                   ( L_solutionsfine - zjL_kn  * solutionsfine )
                enddo
              endif
            enddo
            sources_dn(n) = sum * kn
            if ( Qvary(n) ) then
              L_sources_dn(n,1:Qnums(n)) = sum * L_extinction(n,1:Qnums(n)) + L_sum(1:Qnums(n)) * kn
            endif       
          endif
        enddo
      endif

!  Partials. 8/26/16.
!  5/22/20. Version 2.8.2 Upgrades.
!     ==> Must use vertical distances in Thermal source terms (not path distances). ZJ instead of DJ

      if ( do_enhanced_ps.and.do_Partials ) then
        do ut = 1, npartials
          if ( do_sources_dn_p(ut,v) ) then
            np = partial_layeridx(ut) ; kn = extinction(np)
            path_dn = FOGeometry%losP_paths_LOS(ut,v)
            lostau = kn * path_dn ; if ( lostau.lt.Expcutoff ) lostrans_dn_p(ut) = exp ( - lostau )
!mick fix 3/22/2017 - replaced index "n" with "np" in "Qnums"
!                   - replaced "q" with "1:Qnums(np)" in "L_extinction"
            if ( Qvary(np) ) L_lostrans_dn_p(ut,1:Qnums(np)) = - lostrans_dn_p(ut) * path_dn * L_extinction(np,1:Qnums(np))
            sum = zero ; L_sum = zero
            do j = 1, FOGeometry%nfinedivs_p_LOS_dn(ut,v)
              dj = path_dn - FOGeometry%xfine_p_LOS_dn(ut,j,v) ; xjkn = dj * kn ; tran = exp ( - xjkn ) ! Correct
              zj = FOGeometry%hfine_p_LOS_dn(ut,j,v) ; zjkn = zj * kn ; solutionsfine_p = tcom1(np,1) + zjkn * tcom1(np,2)
              sum  = sum + solutionsfine_p * tran * FOGeometry%wfine_p_LOS_dn(ut,j,v)
              if ( Qvary(np) ) then
                do q = 1, Qnums(np)
                  L_kn = L_extinction(np,q) ; zjL_kn = zj * L_kn
                  L_solutionsfine_p = L_tcom1(np,1,q) + zjkn * L_tcom1(np,2,q) + zjL_kn * tcom1(np,2)
                  L_sum(q)  = L_sum(q) + tran * FOGeometry%wfine_p_LOS_dn(ut,j,v) * &
                             ( L_solutionsfine_p - zjL_kn * solutionsfine_p )
                enddo
              endif
            enddo
            sources_dn_p(ut) = sum * kn
            if ( Qvary(np) ) then
              L_sources_dn_p(ut,1:Qnums(np)) = sum * L_extinction(np,1:Qnums(np)) + L_sum(1:Qnums(np)) * kn
            endif       
          endif
        enddo        
      endif

!  Source function integration
!  ===========================

!  start recursion ( For DSTE term, Use surface emissivity )

      NC =  0
      CUMSOURCE_DN(NC) = zero
      NSTART = 1
      NUT_PREV = NSTART - 1

!  Main intensity loop over all output optical depths
!     NLEVEL = Layer index for given optical depth
!     Cumulative source terms : Loop over layers working upwards from NSTART to level NUT,
!     Check for updating the recursion. Rob Fix Partials 8/26/16

      DO UTA = 1, N_USER_LEVELS
         NUT    = USER_LEVELS(UTA)
         DO N = NSTART, NUT
            NC = N
            CUMSOURCE_DN(NC) = SOURCES_DN(N) + LOSTRANS_DN(N) * CUMSOURCE_DN(NC-1)
         ENDDO
         IF ( Partial_OUTFLAG(UTA) ) THEN
            UT = Partial_OUTINDEX(UTA)
            INTENSITY_DTA_DN(UTA,V) = CUMSOURCE_DN(NC) * LOSTRANS_DN_p(UT) + SOURCES_DN_p(UT)
         ELSE
            INTENSITY_DTA_DN(UTA,V) = CUMSOURCE_DN(NC)
         ENDIF
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
         NUT_PREV = NUT
      ENDDO

!  Profile Wfs (atmospheric term)

      if ( do_profilewfs ) then
        do k = 1, nlayers
          if ( Qvary(k) ) then
            L_CUMSOURCE = zero
            NSTART = 1
            NUT_PREV = NSTART - 1
            DO UTA = 1, N_USER_LEVELS
              NUT    = USER_LEVELS(UTA)
              DO N = NSTART, NUT
                NC = N
                if ( k.eq.n ) then
                  do q = 1, Qnums(k)
                    L_cumsource(q) = L_SOURCES_DN(N,Q)  + &
                                L_LOSTRANS_DN(N,Q) * CUMSOURCE_DN(NC-1) + LOSTRANS_DN(N) * L_CUMSOURCE(Q)
                  enddo
                else
                  do q = 1, Qnums(k)
                    L_cumsource(q) = LOSTRANS_DN(N) * L_CUMSOURCE(Q)
                  enddo
                endif
              ENDDO
              IF ( Partial_OUTFLAG(UTA) ) THEN
                UT = Partial_OUTINDEX(UTA) ; np = partial_layeridx(ut)
                if ( k.eq.np ) then
!mick fix 3/22/2017 - added L_SOURCES_DN_P term
                  do q = 1, Qnums(k)
                    LP_Jacobians_DTA_DN(UTA,V,K,q) = L_cumsource(q) * LOSTRANS_DN_p(UT) &
                      +  L_LOSTRANS_DN_P(UT,Q) * CUMSOURCE_DN(NC) + L_SOURCES_DN_P(UT,q)
                  enddo
                else
                  do q = 1, Qnums(k)
                    LP_Jacobians_DTA_DN(UTA,V,K,q) = L_cumsource(q) * LOSTRANS_DN_p(UT)
                  enddo
                endif
              ELSE
                do q = 1, Qnums(k)
                  LP_Jacobians_dta_DN(UTA,V,K,Q) = L_CUMSOURCE(Q)
                enddo
              endif
              IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
              NUT_PREV = NUT
            ENDDO
          endif
        enddo
      endif

!  End geometry loop

   enddo

!  Finish

   return
end subroutine FO_Thermal_DTRT_ILPS_DN

!

subroutine FO_Thermal_DTRT_ILPS_UPDN &
   ( do_upwelling, do_dnwelling, do_deltam_scaling, do_Partials, do_PlanPar, do_enhanced_ps,  & ! Inputs (Flags)
     Do_Polarized_Emissivity, do_sources_up, do_sources_up_p, do_sources_dn, do_sources_dn_p, & ! Inputs (Flags)
     do_profilewfs, do_surfacewfs, Lvaryflags, Lvarynums, n_surfacewfs,          & ! Inputs (Control, Jacobians)
     nstokes, ngeoms, nlayers, n_user_levels, user_levels, npartials,  & ! Inputs (control output)
     partial_outindex, partial_outflag, partial_layeridx, FOGeometry,  & ! Inputs (control-partial)
     extinction, deltaus, omega, truncfac, bb_input,                   & ! Inputs (Thermal/Optical)
     surfbb, User_Emissivity, User_QUVEmissivity,                      & ! Inputs (Surface)
     L_extinction, L_deltaus, L_omega, L_truncfac,                     & ! Inputs (Optical - Linearized)
     LS_User_Emissivity,  LS_User_QUVEmissivity,                       & ! Inputs (Surface - Linearized)
     intensity_dta_up, intensity_dts, intensity_dta_dn,                & ! Main Outputs
     LP_Jacobians_dta_up, LP_Jacobians_dts_up, LS_Jacobians_dts,       & ! Main Outputs
     LP_Jacobians_dta_dn, tcom1, L_tcom1,                              & ! Main Outputs
     lostrans_up, lostrans_up_p, L_lostrans_up, L_lostrans_up_p,       & ! Other Outputs
     StokesQUV_dts, LP_JacobiansQUV_dts_up, LS_JacobiansQUV_dts )        ! Optional Output. 12/11/17 Rob Add.

!  FO routine for Upwelling and downwelling Direct-thermal-emission (DTE)
!    computation of Radiances and LPS Jacobians. Inputs: geometry, optical properties, Planck functions, emissivity

!  5/22/20. Version 2.8.2 Upgrades.
!    -  Add hfine/hfine_p inputs for correct DT calculation (Outgoing). These are in FOGEOMETRY
!    -  lostrans_up, lostrans_up_p (and linearizations) are now outputs from the Upwelling routine

   implicit none         

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ======

!  General flags

   LOGICAL, Intent(in) :: DO_UPWELLING
   LOGICAL, Intent(in) :: DO_DNWELLING

   logical, Intent(in) :: DO_DELTAM_SCALING
   logical, Intent(in) :: DO_Partials
   logical, Intent(in) :: DO_PLANPAR
   logical, Intent(in) :: DO_ENHANCED_PS

!  Optional inputs for polarized emission. 12/11/17 Rob Add.

   logical, intent(in) :: Do_Polarized_Emissivity

!  Existence flags. 8/26/16. Criticality enters here

   logical, Intent(in)    :: do_sources_up       (maxlayers,MAX_USER_VZANGLES)
   logical, Intent(in)    :: do_sources_up_p     (MAX_PARTLAYERS,MAX_USER_VZANGLES)
   logical, Intent(in)    :: do_sources_dn       (maxlayers,MAX_USER_VZANGLES)
   logical, Intent(in)    :: do_sources_dn_p     (MAX_PARTLAYERS,MAX_USER_VZANGLES)

!  Jacobian flags and control

   LOGICAL, Intent(in) :: do_surfacewfs
   LOGICAL, Intent(in) :: do_profilewfs
   LOGICAL, Intent(in) :: Lvaryflags(maxlayers)
   INTEGER, Intent(in) :: Lvarynums (maxlayers)
   INTEGER, Intent(in) :: n_surfacewfs

!  Numbers

   integer, Intent(in) :: NSTOKES, NLAYERS, NGEOMS, N_USER_LEVELS
   integer, Intent(in) :: USER_LEVELS ( MAX_USER_LEVELS )

!  Numbers for Version 1.5: -->  Partial Control added, 8/26/16

   integer, Intent(in) :: Npartials
   integer, Intent(in) :: partial_layeridx(MAX_PARTLAYERS)
   logical, Intent(in) :: partial_outflag ( MAX_USER_LEVELS )
   integer, Intent(in) :: partial_outindex( MAX_USER_LEVELS )

!  Geometrical inputs
!  ------------------

   Type(VLIDORT_Geometry_FO), Intent(in) :: FOGeometry

!  optical inputs
!  --------------

!  Atmosphere extinction and deltaus

   real(ffp), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(ffp), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(ffp), Intent(in) :: OMEGA       ( MAXLAYERS )
   real(ffp), Intent(in) :: TRUNCFAC    ( MAXLAYERS )

!  Linearized optical inputs

   real(ffp), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_OMEGA       ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_TRUNCFAC    ( MAXLAYERS, max_atmoswfs )

!  Thermal/Surface BB and emissivity

   REAL(ffp), Intent(in) :: SURFBB, BB_INPUT    (0:MAXLAYERS)
   REAL(ffp), Intent(in) :: USER_EMISSIVITY(MAX_USER_VZANGLES),      LS_USER_EMISSIVITY    (MAX_USER_VZANGLES,max_surfacewfs)
   REAL(ffp), Intent(in) :: USER_QUVEmissivity(3,MAX_USER_VZANGLES), LS_USER_QUVEMISSIVITY (MAX_USER_VZANGLES,3,max_surfacewfs)

!  outputs
!  -------

!  Upwelling

   real(ffp), Intent(Out)  :: intensity_dta_up     ( max_user_levels, MAX_USER_VZANGLES )
   real(ffp), Intent(Out)  :: intensity_dts        ( max_user_levels, MAX_USER_VZANGLES )
   real(ffp), Intent(Out)  :: LP_Jacobians_dta_up  ( max_user_levels, MAX_USER_VZANGLES, maxlayers, max_atmoswfs )
   real(ffp), Intent(Out)  :: LP_Jacobians_dts_up  ( max_user_levels, MAX_USER_VZANGLES, maxlayers, max_atmoswfs )
   real(ffp), Intent(Out)  :: LS_Jacobians_dts     ( max_user_levels, MAX_USER_VZANGLES, max_surfacewfs )

!  Downwelling

   real(ffp), Intent(Out)  :: intensity_dta_dn     ( max_user_levels, MAX_USER_VZANGLES )
   real(ffp), Intent(Out)  :: LP_Jacobians_dta_dn  ( max_user_levels, MAX_USER_VZANGLES, maxlayers, max_atmoswfs )

!  5/22/20. Version 2.8.2 Upgrades.
!   ==> Add the Lostrans output

   real(ffp), Intent(Out)  :: lostrans_up      ( maxlayers     , MAX_USER_VZANGLES )
   real(ffp), Intent(Out)  :: lostrans_up_p    ( max_partlayers, MAX_USER_VZANGLES )

   real(ffp), Intent(Out)  :: L_lostrans_up    ( maxlayers,      MAX_USER_VZANGLES, max_atmoswfs )
   real(ffp), Intent(Out)  :: L_lostrans_up_p  ( max_partlayers, MAX_USER_VZANGLES, max_atmoswfs )

!  Thermal setup

   real(ffp), Intent(InOut)   :: tcom1(maxlayers,2)
   real(ffp), Intent(InOut)   :: L_tcom1(maxlayers,2,max_atmoswfs)

!  Optional outputs for polarized emission. 12/11/17 Rob Add.

   real(ffp), Intent(Out) :: StokesQUV_dts          (max_user_levels,3,MAX_USER_VZANGLES)
   real(ffp), Intent(Out) :: LP_JacobiansQUV_dts_up (max_user_levels,3,MAX_USER_VZANGLES,maxlayers,max_atmoswfs)
   real(ffp), Intent(Out) :: LS_JacobiansQUV_dts    (max_user_levels,3,MAX_USER_VZANGLES,max_surfacewfs)

!  Upwelling

   if ( do_upwelling ) then
      call FO_Thermal_DTRT_ILPS_UP &
        ( do_deltam_scaling, do_Partials, do_PlanPar, do_enhanced_ps,      & ! Inputs (Flags)
          Do_Polarized_Emissivity, do_sources_up, do_sources_up_p,         & ! Inputs (Flags)
          do_profilewfs, do_surfacewfs, Lvaryflags, Lvarynums, n_surfacewfs,          & ! Inputs (Control, Jacobians)
          nstokes, ngeoms, nlayers, n_user_levels, user_levels, npartials, & ! Inputs (control output)
          partial_outindex, partial_outflag, partial_layeridx, FOGeometry, & ! Inputs (partial/Geometry)
          extinction, deltaus, omega, truncfac, bb_input,                  & ! Inputs (Thermal/Optical)
          surfbb, User_Emissivity, User_QUVEmissivity,                     & ! Inputs (Surface)
          L_extinction, L_deltaus, L_omega, L_truncfac,                    & ! Inputs (Optical - Linearized)
          LS_User_Emissivity,  LS_User_QUVEmissivity,                      & ! Inputs (Surface - Linearized)
          intensity_dta_up, intensity_dts, LP_Jacobians_dta_up,            & ! Main Outputs
          LP_Jacobians_dts_up, LS_Jacobians_dts, tcom1, L_tcom1,           & ! Main Outputs
          lostrans_up, lostrans_up_p, L_lostrans_up, L_lostrans_up_p,      & ! Other Outputs
          StokesQUV_dts, LP_JacobiansQUV_dts_up, LS_JacobiansQUV_dts )       ! Optional Output. 12/11/17 Rob Add.
   endif

!  Downwelling

   if ( do_dnwelling ) then
      call FO_Thermal_DTRT_ILPS_DN &
        ( do_deltam_scaling, do_Partials, do_PlanPar, do_enhanced_ps,          & ! Inputs (Flags)
          do_sources_dn, do_sources_dn_p, do_profilewfs, Lvaryflags, Lvarynums,     & ! Inputs (Flags/Jac-control)
          ngeoms, nlayers, n_user_levels, user_levels, npartials,           & ! Inputs (control output)
          partial_outindex, partial_outflag, partial_layeridx, FOGeometry,  & ! Inputs (control-partial)
          bb_input, extinction, deltaus, omega, truncfac,                   & ! Inputs (Optical - Regular)
          L_extinction, L_deltaus, L_omega, L_truncfac,                     & ! Inputs (Optical - Linearized)
          intensity_dta_dn, LP_Jacobians_dta_dn, tcom1, L_tcom1 )             ! Output
   endif

!  Finish

   return
end subroutine FO_Thermal_DTRT_ILPS_UPDN

!  End module

end module FO_Thermal_DTRT_ILPS_m


