
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

module FO_Thermal_DTRT_m

!  FUNCTION
!  ========

!  For a given wavelength, this routine will calculate First-Order upwelling+downwelling Intensities(I):

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

!    Version 1a, 01 December 2011, R. Spurr, RT Solutions Inc.
!    Version 1b, 13 February 2012, R. Spurr, RT Solutions Inc.
!    Version 2,  01 June     2012, R. Spurr, RT Solutions Inc.
!    Version 3,  29 October  2012, Extension to Observational multiple geometries
!    Version 4,  31 July     2013, Lattice Multi-geometry
!    Version 5,  07 July     2016, Optional calculation using F Matrices directly. NOT RELEVANT HERE for the THERMAL !!!
!    Version 5,  25 August   2016, Partial-layer output
!    Version 5,  11 December 2017, Optional code for polarized surface emissivity

!  For Thermal Emission sources, the subroutines are

!       FO_Thermal_DTRT_UP   (Upwelling only)
!       FO_Thermal_DTRT_DN   (Downwelling only)
!       FO_Thermal_DTRT_UPDN (Upwelling and Downwelling)

!  5/22/20. Version 2.8.2 Upgrades
!    -  Add hfine/hfine_p inputs for correct DT calculation (Outgoing). These are in FOGeometry.
!    -  lostrans_up, lostrans_up_p  are now outputs from this routine

!  Dependencies
!  ============

   use VLIDORT_PARS_m         , only : zero, one, Expcutoff, MAX_USER_LEVELS, MAXLAYERS, &
                                       MAX_PARTLAYERS, MAXFINELAYERS, MAX_USER_VZANGLES
   use VLIDORT_Setups_def_m

!  All three subroutines public

public

contains

subroutine FO_Thermal_DTRT_UP &
   ( do_deltam_scaling, do_Partials, Do_Polarized_Emissivity,            & ! Inputs (Flags)
     do_PlanPar, do_enhanced_ps, do_sources_up, do_sources_up_p,         & ! Inputs (Flags)
     nstokes, ngeoms, nlayers, n_user_levels, user_levels, npartials,    & ! Inputs (control output)
     partial_outindex, partial_outflag, partial_layeridx,                & ! Inputs (control-partial)
     FOGeometry, bb_input, surfbb, user_emissivity, User_QUVEmissivity,  & ! Inputs (Thermal)
     extinction, deltaus, omega, truncfac,                               & ! Inputs (Optical)
     intensity_dta_up, intensity_dts, StokesQUV_dts,                     & ! Main Outputs
     cumsource_up, tcom1, lostrans_up, lostrans_up_p )                     ! Other Outputs

!  FO routine for Upwelling Direct-thermal-emission (DTE) radiation field.
!    Computation of radiance. Inputs: Control, geometry, thermal and optical properties.

   implicit none         

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ======

!  flags.  Version 1.5:  partials introduced, 8/25/16

   logical, Intent(in) ::  DO_DELTAM_SCALING
   logical, Intent(in) ::  DO_Partials
   logical, Intent(in) ::  DO_PLANPAR
   logical, Intent(in) ::  DO_ENHANCED_PS

!  Optional inputs for polarized emission. 12/11/17 Rob Add.

   logical, intent(in) :: Do_Polarized_Emissivity

!  Existence flags. 8/25/16. Criticality enters here

   logical, Intent(in) :: do_sources_up   (MAXLAYERS,     MAX_USER_VZANGLES)
   logical, Intent(in) :: do_sources_up_p (MAX_PARTLAYERS,MAX_USER_VZANGLES)

!  Numbers

   integer, Intent(in) ::  NSTOKES, NGEOMS, NLAYERS, N_USER_LEVELS
   integer, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  Numbers for Version 1.5: -->  Partial Control added, 8/25/16

   integer, Intent(in) :: Npartials
   integer, Intent(in) :: partial_layeridx( MAX_PARTLAYERS )
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

!  Atmospheric BB functions and Surface BB and emissivity

   real(ffp), Intent(in) :: SURFBB, USER_EMISSIVITY(MAX_USER_VZANGLES)
   real(ffp), Intent(in) :: BB_INPUT (0:MAXLAYERS)
   real(ffp), intent(in) :: USER_QUVEmissivity(3,MAX_USER_VZANGLES)

!  outputs
!  -------

   real(ffp), Intent(Out)  :: intensity_dta_up ( MAX_USER_LEVELS, MAX_USER_VZANGLES )
   real(ffp), Intent(Out)  :: intensity_dts    ( MAX_USER_LEVELS, MAX_USER_VZANGLES )
   real(ffp), Intent(Out)  :: cumsource_up     ( 0:MAXLAYERS,     MAX_USER_VZANGLES )

!  5/22/20. Version 2.8.2 Upgrades
!     ==> Add LOSTRANS arrays to output

   real(ffp), Intent(Out)  :: lostrans_up      ( MAXLAYERS     , MAX_USER_VZANGLES )
   real(ffp), Intent(Out)  :: lostrans_up_p    ( MAX_PARTLAYERS, MAX_USER_VZANGLES )

!  Thermal setup

   real(ffp), Intent(InOut) :: tcom1(maxlayers,2)

!  outputs for polarized emission. 12/11/17 Rob Add.

   real(ffp), Intent(Out)   :: StokesQUV_dts    ( MAX_USER_LEVELS, 3, MAX_USER_VZANGLES )

!  LOCAL
!  -----

!  Source function integration results. Partials added 8/25/16.

   real(ffp)  :: sources_up       ( MAXLAYERS )
   real(ffp)  :: sources_up_p     ( MAX_PARTLAYERS  )

!  Regular_PS or plane-parallel flag

   logical    :: do_RegPSorPP

!  Help

   integer    :: n, uta, nstart, nc, nut, nut_prev, j, v, np, ut
   logical    :: do_regular_ps, layermask_up(maxlayers)
   real(ffp)  :: help, sum, tran, kn, xjkn, dj, zjkn, path_up, Solutionsfine, Solutionsfine_p
   real(ffp)  :: cumsource_dste, t_mult_up(0:2), thermcoeffs(2), tms, lostau, partau

!  Help variables for Optional polarized emissivity. 12/11/17 Rob add.

   integer    :: ns
   logical    :: do_PolEmiss
   real(ffp)  :: QUV_dts(max_user_levels,3), PolEmiss(3,max_user_vzangles), CUMSOURCE_DSTEQUV(3)

!  Zero the output

   CUMSOURCE_UP = zero ; INTENSITY_dta_up = zero ; INTENSITY_dts = zero

!  Regular_PS or plane-parallel flag

   do_regular_ps = .false.
   if ( .not.do_Planpar ) do_regular_ps = .not. do_enhanced_ps
   do_RegPSorPP = (do_regular_ps .or. do_PlanPar)

!  Optional code for polarized emissivity, set Proxies and Initialize. 12/11/17 Rob add.
!mick fix 3/2/2020 - trimmed dimensions for defining PolEmiss
 
   do_PolEmiss = .false.
   if ( Do_Polarized_Emissivity.and. nstokes.gt.1 ) then
      do_PolEmiss = .true. ; ns = nstokes - 1
      PolEmiss(1:ns,1:ngeoms) = User_QUVEmissivity(1:ns,1:ngeoms) ; StokesQUV_dts = zero
   endif

!  Bookkeeping

   NUT = USER_LEVELS(1) + 1
   LAYERMASK_UP = .false.
   LAYERMASK_UP(NUT:NLAYERS) = .true.

!  Thermal setup factors
!     TMS, Initial set of thermal coefficients and TCOM1 variable

   tcom1 = zero
   do n = 1, nlayers
      tms = one - omega(n) 
      if ( do_deltam_scaling ) then
         help = one - truncfac(n) * omega(n)
         tms = tms / help
      endif
      thermcoeffs(1)  = bb_input(n-1)
      thermcoeffs(2)  = (bb_input(n)-bb_input(n-1)) / deltaus(n)
      tcom1(n,1) = thermcoeffs(1) * tms
      tcom1(n,2) = thermcoeffs(2) * tms
   ENDDO

!  5/22/20. Version 2.8.2 Upgrades.
!    ==> Zero the LOSTRANS transmittances

   lostrans_up   = zero
   lostrans_up_p = zero 

!  Start Geometry loop
!  ===================

   do v = 1, ngeoms

!  Zero the local sources

      sources_up   = zero
      sources_up_p = zero

!  Plane/Parallel or Regular-PS Layer integrated source terms
!  ==========================================================

!  Bug Fixed 23 January 2013 (nadir case). Old code commented out and replaced
!  8/25/16. Version 1.5 partials introduced, nadir special case absorbed

!  5/2/20. Version 2.8.2 Upgrades.
!     ==> Add geometry index to Lostrans arrays (now outputs from subroutine)

      if ( do_RegPSorPP ) then
        do n = 1, nlayers
          if ( layermask_up(n).and.do_sources_up(n,v) ) then
            lostau = deltaus(n) / FOGeometry%Mu1_LOS(v)
            if ( lostau .lt. Expcutoff ) lostrans_up(n,v) = exp( - lostau )
            t_mult_up(2) = tcom1(n,2)
            t_mult_up(1) = tcom1(n,1) + t_mult_up(2) * FOGeometry%Mu1_LOS(v)
            sum = t_mult_up(1) + t_mult_up(2) * deltaus(n)
            t_mult_up(0) = - sum
            sources_up(n) = t_mult_up(0) * lostrans_up(n,v)  + t_mult_up(1)
          endif
        enddo
      endif

!  Partials. New Code 8/25/16
!mick fix 3/22/2017 - added "do_Partials" to 1st if condition

!  5/2/20. Version 2.8.2 Upgrades.
!     ==> Add geometry index to Lostrans arrays (now outputs from subroutine)

      if ( do_RegPSorPP .and. do_Partials ) then
        do ut = 1, npartials
          np = partial_layeridx(ut)
          if ( layermask_up(np).and.do_sources_up_p(ut,v) ) then
            path_up = FOGeometry%losW_paths_LOS(np,v) - FOGeometry%losP_paths_LOS(ut,v) ; kn = extinction(np)
            lostau = kn * path_up ; if ( lostau .lt. Expcutoff ) lostrans_up_p(ut,v) = exp( - lostau )
            partau = lostau *  FOGeometry%Mu1_LOS(v)
            t_mult_up(2) = tcom1(np,2)
            t_mult_up(1) = tcom1(np,1) + t_mult_up(2) * FOGeometry%Mu1_LOS(v)
            sum = t_mult_up(1) + t_mult_up(2) * partau
            t_mult_up(0) = - sum
            sources_up_p(ut) = t_mult_up(0) * lostrans_up_p(ut,v)  + t_mult_up(1)
          endif
        enddo
      endif

!  LOS-spherical Layer integrated source terms
!  ===========================================

!  5/22/20. Version 2.8.2 Upgrades.
!     ==> Must use vertical distances in Thermal source terms (not path distances, bug corrected)
!     ==> Add geometry index to Lostrans arrays (now outputs from subroutine)

      if ( do_enhanced_ps ) then
         do n = nlayers, 1, -1
           if ( layermask_up(n) .and. do_sources_up(n,v) ) then
!mick fix 3/22/2017 - replaced index "np" with "n" in "LosW_paths"
             kn = extinction(n) ; path_up = FOGeometry%LosW_paths_LOS(n,v)
             lostau = kn * path_up ; if( lostau.lt.Expcutoff ) lostrans_up(n,v) = exp ( - lostau )
             sum = zero
             do j = 1, FOGeometry%nfinedivs_LOS(n,v)
                dj = FOGeometry%LosW_paths_LOS(n,v) - FOGeometry%xfine_LOS(n,j,v) ; xjkn = dj * kn ; tran = exp ( - xjkn )
! Bug            solutionsfine = tcom1(n,1) + xjkn * tcom1(n,2)
                zjkn = FOGeometry%hfine_LOS_up(n,j,v)*kn ; solutionsfine = tcom1(n,1) + zjkn * tcom1(n,2)
                sum  = sum + solutionsfine * tran * FOGeometry%wfine_LOS(n,j,v)
             enddo
             sources_up(n) = sum * kn
           endif
         enddo
      endif

!  Partials. 8/25/16

!  5/22/20. Version 2.8.2 Upgrades.
!     ==> Must use vertical distances in Thermal source terms (not path distances, bug corrected)
!     ==> Add geometry index to Lostrans arrays (now outputs from subroutine)

      if ( do_enhanced_ps.and.do_Partials ) then
         do ut = 1, npartials
          if ( do_sources_up_p(ut,v) ) then
            np = partial_layeridx(ut) ; kn = extinction(np)
            path_up = FOGeometry%LosW_paths_LOS(np,v)- FOGeometry%LosP_paths_LOS(ut,v)
            lostau = kn * path_up ; if ( lostau.lt.Expcutoff ) lostrans_up_p(ut,v) = exp ( - lostau )
            sum = zero
            do j = 1, FOGeometry%nfinedivs_p_LOS_up(ut,v)
              dj = path_up - FOGeometry%xfine_p_LOS_up(ut,j,v) ; xjkn = dj * kn ; tran = exp ( - xjkn )     ! Correct
! Bug          solutionsfine_p = tcom1(np,1) + xjkn * tcom1(np,2)
              zjkn = FOGeometry%hfine_p_LOS_up(ut,j,v)*kn ; solutionsfine_p = tcom1(np,1) + zjkn * tcom1(np,2)
              sum  = sum + solutionsfine_p * tran * FOGeometry%wfine_p_LOS_up(ut,j,v)
            enddo
            sources_up_p(ut) = sum * kn
          endif
        enddo        
      endif

!  Source function integration
!  ===========================

!  start recursion ( For DSTE term, Use surface emissivity )

      NC =  0
      CUMSOURCE_UP(NC,v) = zero
      CUMSOURCE_DSTE = SURFBB * USER_EMISSIVITY(v)
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  Main loop over all output optical depths
!     NLEVEL = Layer index for given optical depth
!     Cumulative source terms : Loop over layers working upwards from NSTART to level NUT,
!     Check for updating the recursion. Robfix partials, 8/25/16.
!  5/22/20. Version 2.8.2 Upgrades. ==> Add geometry index to Lostrans arrays

      DO UTA = N_USER_LEVELS, 1, -1
         NUT = USER_LEVELS(UTA) + 1
         DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            CUMSOURCE_DSTE     = LOSTRANS_UP(N,v) * CUMSOURCE_DSTE
            CUMSOURCE_UP(NC,V) = LOSTRANS_UP(N,v) * CUMSOURCE_UP(NC-1,V) + SOURCES_UP(N)
         ENDDO
         IF ( Partial_OUTFLAG(UTA) ) THEN
           UT = Partial_OUTINDEX(UTA)
           INTENSITY_DTA_UP(UTA,V) = CUMSOURCE_UP(NC,V) * LOSTRANS_UP_p(UT,v) + SOURCES_UP_p(UT)
           INTENSITY_DTS(UTA,V)    = CUMSOURCE_DSTE * LOSTRANS_UP_p(UT,v)
         ELSE
           INTENSITY_DTA_UP(UTA,V) = CUMSOURCE_UP(NC,V)
           INTENSITY_DTS(UTA,V)    = CUMSOURCE_DSTE
         ENDIF
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
         NUT_PREV = NUT
      ENDDO

!  Optional code for polarized emissivity. 12/11/17 Rob add.
!  5/22/20. Version 2.8.2 Upgrades. ==> Add geometry index to Lostrans arrays

      if ( do_PolEmiss ) then
         CUMSOURCE_DSTEQUV(1:ns) = SURFBB * PolEmiss(1:ns,v)
         NSTART = NLAYERS ; NUT_PREV = NSTART + 1
         DO UTA = N_USER_LEVELS, 1, -1
            NUT = USER_LEVELS(UTA) + 1
            DO N = NSTART, NUT, -1
               CUMSOURCE_DSTEQUV(1:ns) = LOSTRANS_UP(N,v) * CUMSOURCE_DSTEQUV(1:ns)
            ENDDO
            IF ( Partial_OUTFLAG(UTA) ) THEN
              UT = Partial_OUTINDEX(UTA)
              QUV_DTS(UTA,1:ns) = CUMSOURCE_DSTEQUV(1:ns)  * LOSTRANS_UP_p(UT,v)
            ELSE
              QUV_DTS(UTA,1:ns) = CUMSOURCE_DSTEQUV(1:ns)
            ENDIF
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1; NUT_PREV = NUT
         ENDDO
         StokesQUV_dts(1:N_USER_LEVELS,1:ns,v) = QUV_DTS(1:N_USER_LEVELS,1:ns)
      endif

!  End geometry loop

   enddo

!  Finish

   return
end subroutine FO_Thermal_DTRT_UP

!

subroutine FO_Thermal_DTRT_DN &
   ( do_deltam_scaling, do_Partials, do_PlanPar,                 & ! Inputs (Flags)
     do_enhanced_ps, do_sources_dn, do_sources_dn_p,             & ! Inputs (Flags)
     ngeoms, nlayers, n_user_levels, user_levels, npartials,     & ! Inputs (control output)
     partial_outindex, partial_outflag, partial_layeridx,        & ! Inputs (control-partial)
     FOGeometry, BB_input, extinction, deltaus, omega, truncfac, & ! Inputs (Optical)
     intensity_dta_dn, cumsource_dn, tcom1 )                       ! Outputs

!  FO routine for Upwelling Direct-thermal-emission (DTE) radiation field.
!    Computation of radiance. Inputs: Control, geometry, optical properties, thermal.

   implicit none         

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ======

!  flags. Version 1.5:  Partials 8/25/16

   logical, Intent(in) :: DO_DELTAM_SCALING

   logical, Intent(in) :: DO_Partials
   logical, Intent(in) :: DO_PLANPAR
   logical, Intent(in) :: DO_ENHANCED_PS

!  Existence flags. 8/25/16. Criticality enters here

   logical, Intent(in)    :: do_sources_dn       (maxlayers,MAX_USER_VZANGLES)
   logical, Intent(in)    :: do_sources_dn_p     (MAX_PARTLAYERS,MAX_USER_VZANGLES)

!  Numbers

   integer, Intent(in) ::  NGEOMS, NLAYERS, N_USER_LEVELS
   integer, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  Numbers for Version 1.5: -->  Partial Control added, 8/25/16

   integer, Intent(in) :: Npartials
   integer, Intent(in) :: partial_layeridx(MAX_PARTLAYERS)
   logical, Intent(in) :: partial_outflag ( MAX_USER_LEVELS )
   integer, Intent(in) :: partial_outindex( MAX_USER_LEVELS )

!  Geometrical inputs
!  ------------------

   Type(VLIDORT_Geometry_FO), Intent(in) :: FOGeometry

!  optical inputs
!  --------------

!  Atmosphere

   real(ffp), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(ffp), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(ffp), Intent(in) :: OMEGA       ( MAXLAYERS )
   real(ffp), Intent(in) :: TRUNCFAC    ( MAXLAYERS )

!  Atmospheric thermal BB functions

   real(ffp), Intent(in) :: BB_INPUT (0:MAXLAYERS)

!  outputs
!  -------

   real(ffp), Intent(Out)  :: intensity_dta_dn ( max_user_levels,MAX_USER_VZANGLES )
   real(ffp), Intent(Out)  :: cumsource_dn     ( 0:maxlayers,MAX_USER_VZANGLES )

!  Thermal setup

   real(ffp), Intent(InOut)   :: tcom1(maxlayers,2)

!  LOCAL
!  -----

!  Source function integration results. Partials added 8/25/16.

   real(ffp)  :: sources_dn       ( maxlayers )
   real(ffp)  :: lostrans_dn      ( maxlayers )
   real(ffp)  :: sources_dn_p     ( MAX_PARTLAYERS )
   real(ffp)  :: lostrans_dn_p    ( MAX_PARTLAYERS )

!  Regular_PS or plane-parallel flag

   logical    :: do_RegPSorPP

!  Help

   integer    :: n, uta, nstart, nc, nut, nut_prev, j, v, np, ut
   logical    :: do_regular_PS, layermask_dn(maxlayers)

   real(ffp)  :: help, sum, tran, kn, xjkn, dj, zjkn, lostau, partau, path_dn
   real(ffp)  :: t_mult_dn(0:2), thermcoeffs(2), tms, Solutionsfine, Solutionsfine_p

!  Zero the output and the local sources

   CUMSOURCE_DN = zero ; INTENSITY_DTA_DN = zero

!  Regular_PS or plane-parallel flag

   do_regular_ps = .false.
   if ( .not.do_Planpar ) do_regular_ps = .not. do_enhanced_ps
   do_RegPSorPP = (do_regular_ps .or. do_PlanPar)

!  Bookkeeping

   NUT = USER_LEVELS(N_USER_LEVELS) + 1
   IF ( NUT > NLAYERS ) NUT = NLAYERS
   LAYERMASK_DN = .false.
   LAYERMASK_DN(1:NUT) = .true.

!  Thermal setup factors
!     TMS, Initial set of thermal coefficients and TCOM1 variable

   tcom1 = zero
   do n = 1, nlayers
      tms = one - omega(n) 
      if ( do_deltam_scaling ) then
         help = one - truncfac(n) * omega(n)
         tms = tms / help
      endif
      thermcoeffs(1)  = bb_input(n-1)
      thermcoeffs(2)  = (bb_input(n)-bb_input(n-1)) / deltaus(n)
      tcom1(n,1) = thermcoeffs(1) * tms
      tcom1(n,2) = thermcoeffs(2) * tms
   ENDDO

!  Start geometry loop
!  ===================

   do v = 1, ngeoms

!  Zero the local sources

      lostrans_dn = zero    ; sources_dn = zero
      lostrans_dn_p = zero  ; sources_dn_p = zero

!  Plane/Parallel or Regular-PS Layer integrated source terms
!  ==========================================================

!  Bug Fixed 23 January 2013 (nadir case). Old code commented out and replaced
!  8/25/16. Version 1.5 partials introduced, nadir special case absorbed

      if ( do_RegPSorPP ) then
        DO n = 1, nlayers
          if ( layermask_dn(n).and.do_sources_dn(n,v) ) then
            lostau = deltaus(n) / FOGeometry%Mu1_LOS(v)
            if ( lostau .lt. Expcutoff ) lostrans_dn(n) = exp( - lostau )
            t_mult_dn(2)   = tcom1(n,2)
            t_mult_dn(1)   = tcom1(n,1) - t_mult_dn(2) * FOGeometry%Mu1_LOS(v)
            t_mult_dn(0)   = - t_mult_dn(1)
            sources_dn(n)  = t_mult_dn(0) * lostrans_dn(n)
            sum = t_mult_dn(1) + t_mult_dn(2) * deltaus(n)
            sources_dn(n)  = sources_dn(n) + sum
          endif
        enddo
      endif

!  Partials. New Code 8/25/16
!mick fix 3/22/2017 - added "do_Partials" to 1st if condition
!                   - replaced two lines

      if ( do_RegPSorPP .and. do_Partials ) then
        DO ut = 1, npartials
          np = partial_layeridx(ut)
          if ( layermask_dn(np).and.do_sources_dn_p(ut,v) ) then
            path_dn = FOGeometry%losP_paths_LOS(ut,v) ; kn = extinction(np)
            lostau = kn * path_dn ; if ( lostau .lt. Expcutoff ) lostrans_dn_p(ut) = exp( - lostau )
            partau = lostau * FOGeometry%Mu1_LOS(v)
            t_mult_dn(2)   = tcom1(np,2)
            t_mult_dn(1)   = tcom1(np,1) - t_mult_dn(2) * FOGeometry%Mu1_LOS(v)
            t_mult_dn(0)   = - t_mult_dn(1)
            !sources_dn_p(n)  = t_mult_dn(0) * lostrans_dn_p(n)
            sources_dn_p(ut)  = t_mult_dn(0) * lostrans_dn_p(ut)
            sum = t_mult_dn(1) + t_mult_dn(2) * partau
            !sources_dn(n)  = sources_dn(n) + sum
            sources_dn_p(ut) = sources_dn_p(ut) + sum
          endif
        enddo
      endif

!  LOS-spherical Layer integrated source terms
!  ===========================================

!  5/22/20. Version 2.8.2 Upgrades.
!     ==> Must use vertical distances in Thermal source terms (not path distances, bug corrected)

      if ( do_enhanced_ps ) then
         do n = nlayers, 1, -1
           if ( layermask_dn(n) .and. do_sources_dn(n,v) ) then
!mick fix 3/22/2017 - replaced index "np" with "n" in "LosW_paths"
             kn = extinction(n) ; path_dn = FOGeometry%LosW_paths_LOS(n,v)
             lostau = kn * path_dn ; if( lostau.lt.Expcutoff ) lostrans_dn(n) = exp ( - lostau )
             sum = zero
             do j = 1, FOGeometry%nfinedivs_LOS(n,v)
                dj = FOGeometry%LosW_paths_LOS(n,v) - FOGeometry%xfine_LOS(n,j,v) ; xjkn = dj * kn ; tran = exp ( - xjkn )
!                solutionsfine = tcom1(n,1) + xjkn * tcom1(n,2)
                zjkn = FOGeometry%hfine_LOS_dn(n,j,v) * kn ; solutionsfine = tcom1(n,1) + zjkn * tcom1(n,2)
                sum  = sum + solutionsfine * tran * FOGeometry%wfine_LOS(n,j,v)
             enddo
             sources_dn(n) = sum * kn
           endif
         enddo
      endif

!  Partials
!mick fix 3/2/2020 - apparent cut & paste error: replaced "whole layer" do loop material
!                    with "partial layer" do loop material
 
!  5/22/20. Version 2.8.2 Upgrades.
!     ==> Must use vertical distances in Thermal source terms (not path distances, bug corrected)

      if ( do_enhanced_ps.and.do_Partials ) then
         do ut = 1, npartials
           if ( do_sources_dn_p(ut,v) ) then
             np = partial_layeridx(ut) ; kn = extinction(np)
             path_dn = FOGeometry%losP_paths_LOS(ut,v)
             lostau = kn * path_dn ; if ( lostau.lt.Expcutoff ) lostrans_dn_p(ut) = exp ( - lostau )
             sum = zero
             do j = 1, FOGeometry%nfinedivs_p_LOS_dn(ut,v)
               dj = path_dn - FOGeometry%xfine_p_LOS_dn(ut,j,v) ; xjkn = dj * kn ; tran = exp ( - xjkn )  ! Correct
!              solutionsfine_p = tcom1(np,1) + xjkn * tcom1(np,2)
               zjkn = FOGeometry%hfine_p_LOS_dn(ut,j,v)*kn ; solutionsfine_p = tcom1(np,1) + zjkn * tcom1(np,2)
               sum  = sum + solutionsfine_p * tran * FOGeometry%wfine_p_LOS_dn(ut,j,v)
             enddo
             sources_dn_p(ut) = sum * kn
           endif
         enddo        
      endif

!  Source function integration
!  ===========================

!  start recursion

      NC =  0
      CUMSOURCE_DN(NC,V) = zero
      NSTART = 1
      NUT_PREV = NSTART - 1

!  Main loop over all output optical depths
!     NLEVEL = Layer index for given optical depth
!     Cumulative source terms : Loop over layers working Downn from NSTART to NUT
!     Check for dndating the recursion. Rob Fix Partials 8/25/16.

      DO UTA = 1, N_USER_LEVELS
         NUT    = USER_LEVELS(UTA)
         DO N = NSTART, NUT
            NC = N
            CUMSOURCE_DN(NC,V) = SOURCES_DN(N) + LOSTRANS_DN(N) * CUMSOURCE_DN(NC-1,V)
         ENDDO
         IF ( Partial_OUTFLAG(UTA) ) THEN
            UT = Partial_OUTINDEX(UTA)
            INTENSITY_DTA_DN(UTA,V) = CUMSOURCE_DN(NC,V) * LOSTRANS_DN_p(UT) + SOURCES_DN_p(UT)
         ELSE
            INTENSITY_DTA_DN(UTA,V) = CUMSOURCE_DN(NC,v)
         ENDIF
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
         NUT_PREV = NUT
      ENDDO

!  End geometry loop

   enddo

!  Finish

   return

end subroutine FO_Thermal_DTRT_DN

!

subroutine FO_Thermal_DTRT_UPDN   &
   ( do_upwelling, do_dnwelling, do_deltam_scaling, do_Partials, do_PlanPar, do_enhanced_ps,  & ! Inputs (Flags)
     Do_Polarized_Emissivity, do_sources_up, do_sources_up_p, do_sources_dn, do_sources_dn_p,  & ! Inputs (Flags)
     nstokes, ngeoms, nlayers, n_user_levels, user_levels, npartials,   & ! Inputs (control output)
     partial_outindex, partial_outflag, partial_layeridx,               & ! Inputs (control-partial)
     FOGeometry, bb_input, surfbb, user_emissivity, User_QUVEmissivity, & ! Inputs (Geometry/thermal)
     extinction, deltaus, omega, truncfac,                              & ! Inputs (Optical)
     intensity_dta_up, intensity_dts, cumsource_up, lostrans_up, lostrans_up_p, & ! Main Outputs
     intensity_dta_dn, cumsource_dn, tcom1, StokesQUV_dts )                       ! Main outputs

!  FO routine for Upwelling and Downwelling Direct-thermal-emission (DTE) radiation field.
!    Computation of radiance. Inputs: Control, geometry, optical properties, thermal.

!  5/22/20. Version 2.8.2 Upgrades.
!    -  Add hfine/hfine_p inputs for correct DT calculation (Outgoing). These are in FOGeometry
!    -  lostrans_up, lostrans_up_p  are now outputs from the Upwelling routine

   implicit none         

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ======

!  flags

   logical, Intent(in) :: DO_UPWELLING
   logical, Intent(in) :: DO_DNWELLING
   logical, Intent(in) :: DO_DELTAM_SCALING

   logical, Intent(in) :: DO_Partials
   logical, Intent(in) :: DO_PLANPAR
   logical, Intent(in) :: DO_ENHANCED_PS

!  Optional inputs for polarized emission. 12/11/17 Rob Add.

   logical, intent(in) :: Do_Polarized_Emissivity

!  Existence flags. 8/25/16. Criticality enters here

   logical, Intent(in)    :: do_sources_up       (maxlayers,MAX_USER_VZANGLES)
   logical, Intent(in)    :: do_sources_up_p     (MAX_PARTLAYERS,MAX_USER_VZANGLES)
   logical, Intent(in)    :: do_sources_dn       (maxlayers,MAX_USER_VZANGLES)
   logical, Intent(in)    :: do_sources_dn_p     (MAX_PARTLAYERS,MAX_USER_VZANGLES)

!  Numbers

   integer, Intent(in) ::  NSTOKES, NGEOMS, NLAYERS, N_USER_LEVELS
   integer, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  Numbers for Version 1.5: -->  Partial Control added, 8/25/16

   integer, Intent(in) :: Npartials
   integer, Intent(in) :: partial_layeridx(MAX_PARTLAYERS)
   logical, Intent(in) :: partial_outflag ( MAX_USER_LEVELS )
   integer, Intent(in) :: partial_outindex( MAX_USER_LEVELS )

!  Geometrical inputs
!  ------------------

   Type(VLIDORT_Geometry_FO), Intent(in) :: FOGeometry

!  optical inputs
!  --------------

!  Atmosphere

   real(ffp), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(ffp), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(ffp), Intent(in) :: OMEGA       ( MAXLAYERS )
   real(ffp), Intent(in) :: TRUNCFAC    ( MAXLAYERS )

!  Atmospheric BB functions and Surface BB and emissivity

   real(ffp), Intent(in) :: SURFBB, USER_EMISSIVITY(MAX_USER_VZANGLES)
   real(ffp), Intent(in) :: BB_INPUT (0:MAXLAYERS)
   real(ffp), intent(in) :: USER_QUVEmissivity(3,MAX_USER_VZANGLES)

!  outputs
!  -------

   real(ffp), Intent(Out)  :: intensity_dta_up     ( max_user_levels,MAX_USER_VZANGLES )
   real(ffp), Intent(Out)  :: intensity_dts        ( max_user_levels,MAX_USER_VZANGLES )
   real(ffp), Intent(Out)  :: cumsource_up         ( 0:maxlayers,MAX_USER_VZANGLES )
   real(ffp), Intent(Out)  :: intensity_dta_dn     ( max_user_levels,MAX_USER_VZANGLES )
   real(ffp), Intent(Out)  :: cumsource_dn         ( 0:maxlayers,MAX_USER_VZANGLES )

!  Thermal setup

   real(ffp), Intent(InOut)   :: tcom1(maxlayers,2)

!  5/22/20. Version 2.8.2 Upgrades.
!   ===>  Add the Lostrans output

   real(ffp), Intent(Out)  :: lostrans_up      ( MAXLAYERS     ,MAX_USER_VZANGLES )
   real(ffp), Intent(Out)  :: lostrans_up_p    ( MAX_PARTLAYERS,MAX_USER_VZANGLES )

!  outputs for polarized emission. 12/11/17 Rob Add.

   real(ffp), Intent(Out)   :: StokesQUV_dts    ( MAX_USER_LEVELS, 3, MAX_USER_VZANGLES )

!  Upwelling
!  ---------

   if ( do_upwelling ) then
      Call FO_Thermal_DTRT_UP &
        ( do_deltam_scaling, do_Partials, Do_Polarized_Emissivity,            & ! Inputs (Flags)
          do_PlanPar, do_enhanced_ps, do_sources_up, do_sources_up_p,         & ! Inputs (Flags)
          nstokes, ngeoms, nlayers, n_user_levels, user_levels, npartials,    & ! Inputs (control output)
          partial_outindex, partial_outflag, partial_layeridx,                & ! Inputs (control-partial)
          FOGeometry, bb_input, surfbb, user_emissivity, User_QUVEmissivity,  & ! Inputs (Thermal)
          extinction, deltaus, omega, truncfac,                               & ! Inputs (Optical)
          intensity_dta_up, intensity_dts, StokesQUV_dts,                     & ! Main Outputs
          cumsource_up, tcom1, lostrans_up, lostrans_up_p )                     ! Other Outputs
   endif

   if ( do_dnwelling ) then
      call FO_Thermal_DTRT_DN &
        ( do_deltam_scaling, do_Partials, do_PlanPar,                 & ! Inputs (Flags)
          do_enhanced_ps, do_sources_dn, do_sources_dn_p,             & ! Inputs (Flags)
          ngeoms, nlayers, n_user_levels, user_levels, npartials,     & ! Inputs (control output)
          partial_outindex, partial_outflag, partial_layeridx,        & ! Inputs (control-partial)
          FOGeometry, BB_input, extinction, deltaus, omega, truncfac, & ! Inputs (Optical)
          intensity_dta_dn, cumsource_dn, tcom1 )                       ! Outputs
   endif

!  Finish

   return
end subroutine FO_Thermal_DTRT_UPDN

!  End module

end module FO_Thermal_DTRT_m

