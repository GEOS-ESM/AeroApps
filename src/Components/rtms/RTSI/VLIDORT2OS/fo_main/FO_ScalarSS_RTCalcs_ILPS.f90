! ###########################################################
! #                                                         #
! #                    THE LIDORT FAMILY                    #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #       --         -        -        -         -          #
! #                                                         #
! ###########################################################

! ###########################################################
! #                                                         #
! #  Author :      Robert. J. D. Spurr                      #
! #                                                         #
! #  Address :     RT Solutions, Inc.                       #
! #                9 Channing Street                        #
! #                Cambridge, MA 02138, USA                 #
! #                                                         #
! #  Tel:          (617) 492 1183                           #
! #  Email :        rtsolutions@verizon.net                 #
! #                                                         #
! ###########################################################

! ###########################################################
! #                                                         #
! #     FIRST-ORDER SCALAR MODEL (EXACT SINGLE SCATTERING)  #
! #                                                         #
! #  This Version :   1.3 F90                               #
! #  Release Date :   March 2013                            #
! #                                                         #
! #   Version 1.1,  13 February 2012, First Code            #
! #   Version 1.2,  01 June     2012, Modularization        #
! #   Version 1.3,  29 October  2012, Multiple geometries   #
! #                                                         #
! ###########################################################

! ##########################################################
! #                                                        #
! #   This Version of FIRST_ORDER comes with a GNU-style   #
! #   license. Please read the license carefully.          #
! #                                                        #
! ##########################################################

module FO_ScalarSS_RTCalcs_ILPS_m

!  For a given wavelength, this routine will calculate upwelling and downwelling
!  First Order Intensities(I), and any number of LPS Jacobians (profile/surface)

!     (1) For the Atmospheric Solar Single-scatter and Surface Direct-Beam (SS) sources.
!  This is based on Precalculated Geometrical quantities and appropriate Optical properties.

!  This will perform Enhanced-PS calculations (incoming solar and outgoing LOS-path sphericity) 
!  This will perform Regular-PS  calculations (plane-parallel or incoming solar pseudo-spherical)

!  This is Versions 1-3, without Partials. Code is stand alone with no dependencies.
!    Version 1a, 01 December 2011, R. Spurr, RT Solutions Inc.
!    Version 1b, 02 February 2012, R. Spurr, RT Solutions Inc.
!    Version 2 , 01 June     2012, R. Spurr, RT Solutions Inc.
!    Version 3 , 19 December 2012, Extension to Multiple Geometries, LCS separation

!  For Solar sources, the subroutines are
!       SS_Integral_ILPS_UP   (Upwelling only)
!       SS_Integral_ILPS_DN   (Downwelling only)
!       SS_Integral_ILPS_UPDN (Upwelling and Downwelling)

!  All subroutines public

public

contains

subroutine SS_Integral_ILPS_UP &
   ( maxgeoms, maxlayers, maxfinelayers, maxmoments_input,                             & ! Inputs (dimensioning)
     max_user_levels, max_atmoswfs, max_surfacewfs,                                    & ! Inputs (dimensioning)
     do_deltam_scaling, do_PlanPar, do_regular_ps, do_enhanced_ps, doNadir,            & ! Inputs (Flags - General  )
     do_profilewfs, do_surfacewfs, n_surfacewfs, Lvaryflags, Lvarynums, Lvarymoms,     & ! Inputs (Control - Jacobian )
     ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,           & ! Inputs (control - output)
     reflec, extinction, deltaus, omega, truncfac, phasmoms, flux,                     & ! Inputs (Optical - Regular)
     LS_reflec, L_extinction, L_deltaus, L_omega, L_truncfac, L_phasmoms,              & ! Inputs (Optical - Linearized)
     Mu0, Mu1, LegPoly_up, NCrit, xfine, wfine, csqfine, cotfine,                      & ! Inputs (Geometry)
     Raycon, cota, sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,                 & ! Inputs (Geometry)
     intensity_up, intensity_db, LP_Jacobians_up, LP_Jacobians_db, LS_Jacobians_db )    ! Output

!  Stand-alone routine for Upwelling Solar-beam Single-scatter (SS)
!    computation of Radiances and LPS Jacobians. Inputs: geometry, spherical functions, optical properties.

!  This version, revised by R. Spurr, 01 June 2012
!   Extension to multiple geometries, 19 December 2012

   implicit none         

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs
!  ======

!  Dimensions

   integer, Intent(in) :: maxgeoms
   integer, Intent(in) :: maxlayers
   integer, Intent(in) :: maxfinelayers
   integer, Intent(in) :: maxmoments_input
   integer, Intent(in) :: max_user_levels
   INTEGER, Intent(in) :: max_atmoswfs
   INTEGER, Intent(in) :: max_surfacewfs

!  flags

   logical, Intent(in) :: DO_DELTAM_SCALING
   logical, Intent(in) :: DO_PLANPAR
   logical, Intent(in) :: DO_REGULAR_PS
   logical, Intent(in) :: DO_ENHANCED_PS
   logical, Intent(in) :: DONADIR(MAXGEOMS)

!  Jacobian Flags

   LOGICAL, Intent(in) :: do_surfacewfs
   LOGICAL, Intent(in) :: do_profilewfs

!  Layer and Level Control Numbers, Number of Moments

   integer, Intent(in) :: NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   integer, Intent(in) :: NGEOMS, NMOMENTS_INPUT

   integer, Intent(in) :: N_USER_LEVELS
   integer, Intent(in) :: USER_LEVELS ( MAX_USER_LEVELS )

!  Jacobian control

   INTEGER, Intent(in) :: n_surfacewfs
   LOGICAL, Intent(in) :: Lvaryflags(maxlayers)
   INTEGER, Intent(in) :: Lvarynums (maxlayers)
   LOGICAL, Intent(in) :: Lvarymoms (maxlayers,max_atmoswfs)

!  optical inputs
!  --------------

!  Atmosphere

   real(fpk), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(fpk), Intent(in) :: OMEGA       ( MAXLAYERS )
   real(fpk), Intent(in) :: TRUNCFAC    ( MAXLAYERS )
   real(fpk), Intent(in) :: PHASMOMS    ( MAXLAYERS,0:MAXMOMENTS_INPUT )

!  Solar Flux and Surface reflectivity (Could be the albedo)

   real(fpk), Intent(in) :: REFLEC(MAXGEOMS), FLUX

!  Linearized optical inputs

   real(fpk), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
   real(fpk), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
   real(fpk), Intent(in) :: L_OMEGA       ( MAXLAYERS, max_atmoswfs )
   real(fpk), Intent(in) :: L_TRUNCFAC    ( MAXLAYERS, max_atmoswfs )
   real(fpk), Intent(in) :: L_PHASMOMS    ( MAXLAYERS,0:MAXMOMENTS_INPUT, max_atmoswfs )
   real(fpk), Intent(in) :: LS_REFLEC     ( MAXGEOMS, max_surfacewfs)

!  Geometrical inputs
!  ------------------

!       Ray constant, Cotangents
!       Mu0 = cos(theta_boa), required for the surface term (both regular and enhanced)
!       Mu1 = cos(alpha_boa), required for the Regular PS only
!       solar paths, Legendre Polynomials

   integer  , Intent(in)  :: NCrit(maxgeoms)
   real(fpk), Intent(in)  :: Raycon(maxgeoms), cota(0:maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: Mu0(maxgeoms), Mu1(maxgeoms)

!  solar paths 

   integer  , Intent(in)  :: ntraverse  (0:maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: sunpaths   (0:maxlayers,maxlayers,maxgeoms)
   integer  , Intent(in)  :: ntraverse_fine(maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: sunpaths_fine (maxlayers,maxlayers,maxfinelayers,maxgeoms)

!  Legendres

   real(fpk), Intent(in)  :: LegPoly_up(0:maxmoments_input,maxgeoms)

!  LOS Quadratures for Enhanced PS

   real(fpk), Intent(in)  :: xfine   (maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: wfine   (maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: csqfine (maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: cotfine (maxlayers,maxfinelayers,maxgeoms)

!  outputs
!  -------

   real(fpk), Intent(Out)  :: intensity_up     ( max_user_levels, maxgeoms )
   real(fpk), Intent(Out)  :: intensity_db     ( max_user_levels, maxgeoms )
   real(fpk), Intent(Out)  :: LP_Jacobians_up  ( max_user_levels, maxgeoms, maxlayers, max_atmoswfs )
   real(fpk), Intent(Out)  :: LP_Jacobians_db  ( max_user_levels, maxgeoms, maxlayers, max_atmoswfs )
   real(fpk), Intent(Out)  :: LS_Jacobians_db  ( max_user_levels, maxgeoms, max_surfacewfs )

!  LOCAL
!  -----

!  Attenuations

   real(fpk)  :: attenuations     (0:maxlayers)
   real(fpk)  :: LP_attenuations   (0:maxlayers,maxlayers,max_atmoswfs)

!  Solutions for Enhaced-PS

   real(fpk)  :: Solutions_fine   (maxlayers,maxfinelayers)
   real(fpk)  :: LP_solutions_fine (maxlayers,maxfinelayers,maxlayers,max_atmoswfs)

!  Scattering

   real(fpk)  :: tms            (maxlayers)
   real(fpk)  :: exactscat_up   (maxlayers)
   real(fpk)  :: L_tms          (maxlayers,max_atmoswfs)
   real(fpk)  :: L_exactscat_up (maxlayers,max_atmoswfs)

!  Source function integration results

   real(fpk) :: sources_up         ( maxlayers )
   real(fpk) :: LP_sources_up       ( maxlayers,maxlayers,max_atmoswfs )

   real(fpk) :: lostrans_up        ( maxlayers )
   real(fpk) :: L_lostrans_up      ( maxlayers,max_atmoswfs )

   real(fpk) :: cumsource_db      ( 0:maxlayers )
   real(fpk) :: cumsource_up      ( 0:maxlayers )
   real(fpk) :: L_cumsource       ( max_atmoswfs )
   real(fpk) :: LS_cumsource      ( max_surfacewfs )

!  Regular_PS or plane-parallel flag

   logical    :: do_RegPSorPP

!  Help

   integer    :: n, k, j, q, L, v, uta, nstart, nc, nut, nut_prev, Qnums(maxlayers)
   logical    :: layermask_up(maxlayers), Qvary(maxlayers)

   real(fpk)  :: argum(maxfinelayers), tran(maxfinelayers), func(maxfinelayers)
   real(fpk)  :: cons, help, sum, tran_1, kn, ke, factor1, factor2, m4, term1
   real(fpk)  :: L_help, L_sum, L_tran, L_func, L_factor1, L_factor2
   real(fpk)  :: cot_1, cot_2, multiplier, suntau(0:maxlayers), lostau
   real(fpk)  :: L_multiplier, LP_suntau(0:maxlayers,maxlayers,max_atmoswfs), L_lostau(max_atmoswfs)
   real(fpk)  :: attenuations_fine, L_attenuations_fine, sumd, L_sumd, L_exactscat

   real(fpk), parameter  :: cutoff = 88.0_fpk
   real(fpk), parameter  :: zero   = 0.0_fpk
   real(fpk), parameter  :: one    = 1.0_fpk

!  Zero the output

   INTENSITY_UP    = zero
   LP_JACOBIANS_UP = zero

   INTENSITY_DB    = zero
   LP_JACOBIANS_DB = zero
   LS_JACOBIANS_DB = zero

!  Regular_PS or plane-parallel flag

   do_RegPSorPP = (do_regular_ps .or. do_PlanPar)

!  Bookkeeping

   NUT = USER_LEVELS(1) + 1
   LAYERMASK_UP = .false.
   LAYERMASK_UP(NUT:NLAYERS) = .true.

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

      lostrans_up   = zero  ; sources_up    = zero ; cumsource_up = zero
      L_lostrans_up = zero  ; LP_sources_up = zero

!  Scattering functions and Linearization

      do n = 1, nlayers
         if ( layermask_up(n) ) then
            sum = zero
            do L = 0, nmoments_input
               sum = sum + LegPoly_Up(L,V) * phasmoms(n,L)
            enddo
            exactscat_up(n) = sum * tms(n)
            if ( Qvary(n) ) then
               do q = 1, Qnums(n)
                  if ( Lvarymoms(n,q) ) then
                     L_sum = zero
                     do L = 0, nmoments_input
                        L_sum = L_sum + LegPoly_Up(L,V) * L_phasmoms(n,L,q)
                     enddo
                     L_exactscat_up(n,q) = L_sum * tms(n) + sum * L_tms(n,q)
                  else
                     L_exactscat_up(n,q) = sum * L_tms(n,q)
                  endif
               enddo
            endif               
         endif
      enddo

!  Attenuations and Solar solutions
!  ================================

!  Initialize, only to layer Ncrit if applicable

      Attenuations = ZERO   ; LP_Attenuations = ZERO
      Suntau       = ZERO   ; LP_suntau       = ZERO
      Solutions_fine = zero ; LP_Solutions_fine = zero

      nstart = nlayers ; if (Ncrit(v).ne.0) nstart = nCrit(v)

!  Attenuations to End points (including TOA). All representations
!  ===============================================================

      do n = 0, nlayers
         sumd = ZERO
         do k = 1, ntraverse(n,v)
            sumd = sumd + extinction(k) * sunpaths(n,k,v)
         enddo
         suntau(n) = sumd
         If (sumd .lt. cutoff ) Attenuations(n) = exp( - sumd )
         if ( do_profilewfs ) then
            do k = 1, nlayers
               if ( Qvary(k) .and. k.le.ntraverse(n,v) ) then
                  do q = 1, Qnums(k)
                     LP_suntau(n,k,q) = L_extinction(k,q) * sunpaths(n,k,v)
                     LP_Attenuations(n,k,q) = - Attenuations(n) * LP_suntau(n,k,q)
                  enddo
               endif
            enddo
         endif
      enddo

!  Enhanced-spherical, fine-layer attenuations
!  ===========================================

      if ( do_enhanced_ps ) then
         do n = 1, nstart
            if ( layermask_up(n) ) then
               do j = 1, nfinedivs(n,v)
                  sumd = ZERO
                  do k = 1, ntraverse_fine(n,j,v)
                     sumd = sumd + extinction(k) * sunpaths_fine(n,k,j,v)
                  enddo
                  if (sumd .lt. cutoff ) Attenuations_fine = exp( - sumd )
                  Solutions_fine(n,j) = exactscat_up(n) * Attenuations_fine
                  if ( do_profilewfs ) then
                     do k = 1, nlayers
                        if ( Qvary(k) .and. k.le.ntraverse_fine(n,j,v) ) then
                           do q = 1, Qnums(k)
                               L_exactscat = zero ; if ( k.eq.n) L_exactscat = L_exactscat_up(n,q)
                               L_sumd = L_extinction(k,q) * sunpaths_fine(n,k,j,v)
                               L_Attenuations_fine = - Attenuations_fine * L_sumd 
                               LP_Solutions_fine(n,j,k,q) = L_exactscat       *   Attenuations_fine   + &
                                                              exactscat_up(n) * L_Attenuations_fine
                           enddo
                        endif
                     enddo
                  endif
               enddo
            endif
         enddo
      endif

!  Plane-Parallel or Regular-PS (Average secant formulation): Layer integrated Solar sources
!  =========================================================================================

      if ( do_RegPSorPP ) then
         do n = nlayers, 1, -1

!  LOS transmittance (not for the Horizontal case)

            if ( Mu1(v) .gt. zero ) then
               lostau = deltaus(n)  / Mu1(v)
               if ( lostau .lt. cutoff ) lostrans_up(n) = exp( - lostau )
               if ( do_profilewfs ) then
                  if ( Qvary(n) ) then
                     do q = 1, Qnums(n)
                        L_lostau(q)        = L_deltaus(n,q) / Mu1(v)
                        L_lostrans_up(n,q) = - L_lostau(q) * lostrans_up(n)
                     enddo
                  endif
               endif
            endif

!  Sources, general case

            if ( layermask_up(n) .and. n.le.nstart  ) then
              if ( Mu1(v) .gt. zero ) then
                factor1 = Attenuations(n-1) - Attenuations(n)*lostrans_up(n)
                factor2 = one + (suntau(n) - suntau(n-1))/lostau
                multiplier = factor1 / factor2
                sources_up(n) = exactscat_up(n) * multiplier
                if ( do_profilewfs ) then
                  do k = 1, nlayers
                     if ( Qvary(k).and.k.le.ntraverse(n,v) ) then
                        do q = 1, Qnums(k)
                           L_factor1 = LP_Attenuations(n-1,k,q) - LP_Attenuations(n,k,q)*lostrans_up(n)
                           L_factor2 = (LP_suntau(n,k,q) - LP_suntau(n-1,k,q))/lostau
                           if ( k.eq.n) then
                              L_factor1 = L_factor1 - Attenuations(n)*L_lostrans_up(n,q)
                              L_factor2 = L_factor2 - (factor2 - one)*L_lostau(q)/lostau 
                              L_multiplier = ( L_factor1 - multiplier*L_factor2 ) / factor2
                              LP_sources_up(n,k,q) = L_exactscat_up(n,q) *   multiplier + &
                                                      exactscat_up(n)   * L_multiplier
                           else
                              L_multiplier = ( L_factor1 - multiplier * L_factor2 ) / factor2
                              LP_sources_up(n,k,q) = exactscat_up(n)   * L_multiplier
                           endif
                        enddo
                     endif
                  enddo
                endif
              endif
            endif

!  Sources, special case (horizonal view)

            if ( layermask_up(n) .and. n.le.nstart  ) then
              if ( Mu1(v) .eq. zero ) then
                factor1 = Attenuations(n-1) - Attenuations(n)*lostrans_up(n)
                sources_up(n) = exactscat_up(n) * factor1
                if ( do_profilewfs ) then
                  do k = 1, nlayers
                     if ( Qvary(k).and.k.le.ntraverse(n,v) ) then
                        do q = 1, Qnums(k)
                           L_factor1 = LP_Attenuations(n-1,k,q)
                           if ( k.eq.n) then
                              LP_sources_up(n,k,q) = L_exactscat_up(n,q) *   factor1  + &
                                                       exactscat_up(n)   * L_factor1
                           else
                              LP_sources_up(n,k,q) = exactscat_up(n)   * L_factor1
                           endif
                        enddo
                     endif
                  enddo
                endif
              endif
            endif

!  End layers and regular-PS formulation

         enddo
      endif

!  Enhanced PS: special case (nadir viewing). Layer integrated Solar sources
!  =========================================================================

      if ( do_enhanced_ps .and. doNadir(v) ) then
         do n = nlayers, 1, -1

!  LOS transmittance

            kn     = extinction(n)
            lostau = deltaus(n)
            if ( lostau .lt. cutoff ) lostrans_up(n) = exp( - lostau )
            if ( do_profilewfs ) then
               if ( Qvary(n) ) then
                  do q = 1, Qnums(n)
                     L_lostau(q)        = L_deltaus(n,q)
                     L_lostrans_up(n,q) = - L_lostau(q) * lostrans_up(n)
                  enddo
               endif
            endif

!  Sources

            if ( layermask_up(n) .and. n.le.nstart  ) then
               sum = zero
               do j = 1, nfinedivs(n,v)
                  argum(j) = xfine(n,j,v)
                  tran(j)  = exp ( - argum(j) * kn )
                  func(j)  = solutions_fine(n,j) * tran(j)
                  sum = sum + func(j) * wfine(n,j,v)
               enddo
               sources_up(n) = sum * kn
               if ( do_profilewfs ) then
                  do k = 1, nlayers
                     if ( Qvary(k).and.k.le.ntraverse(n,v) ) then
                        do q = 1, Qnums(k)
                           if ( k.eq.n ) then
                              L_sum = zero
                              do j = 1, nfinedivs(n,v)
                                 L_tran = - argum(j) * L_extinction(n,q)
                                 L_func = LP_solutions_fine(n,j,k,q) * tran(j) + L_tran * func(j)
                                 L_sum  = L_sum + L_func * wfine(n,j,v)
                              enddo
                              LP_sources_up(n,k,q)  = L_sum * kn + L_extinction(N,q) * sum
                           else
                              L_sum = zero
                              do j = 1, nfinedivs(n,v)
                                 L_func = LP_solutions_fine(n,j,k,q)  * tran(j)
                                 L_sum  = L_sum + L_func * wfine(n,j,v)
                              enddo
                              LP_sources_up(n,k,q)  = L_sum * kn
                           endif
                        enddo
                     endif
                  enddo
               endif
            endif

!  End layer loop and Nadir Enhanced PS case

         enddo
      endif

!  Enhanced PS: General case. Layer integrated Solar sources
!  =========================================================

      if ( do_enhanced_ps .and. .not. doNadir(v) ) then
         do n = nlayers, 1, -1

!  LOS transmittance

            cot_2 = cota(n-1,v) ; cot_1 = cota(n,v)
            kn = extinction(n) ;  ke = raycon(v) * kn ; cons = raycon(v) * ( cot_2 - cot_1 )
            tran_1 = kn * cons
            if ( tran_1 .lt. cutoff ) lostrans_up(n) = exp ( - tran_1 )
            if ( do_profilewfs ) then
               if ( Qvary(n) ) then
                  do q = 1, Qnums(n)
                     L_lostau(q)        = L_extinction(n,q) * cons
                     L_lostrans_up(n,q) = - L_lostau(q) * lostrans_up(n)
                  enddo
               endif
            endif

!  Sources

            if ( layermask_up(n) .and. n.le.nstart  ) then
               sum = zero
               do j = 1, nfinedivs(n,v)
                  argum(j) = Raycon(v) * ( cot_2 - cotfine(n,j,v) )
                  tran(j)  = exp ( - kn * argum(j) )
                  func(j)  = solutions_fine(n,j) * csqfine(n,j,v) * tran(j)
                  sum      = sum + func(j) * wfine(n,j,v)
               enddo
               sources_up(n) = sum * ke 
               if ( do_profilewfs ) then
                  do k = 1, nlayers
                     if ( Qvary(k).and.k.le.ntraverse(n,v) ) then
                        do q = 1, Qnums(k)
                           if ( k.eq.n ) then
                              L_sum = zero
                              do j = 1, nfinedivs(n,v)
                                 L_tran = - argum(j) * L_extinction(n,q)
                                 L_func = LP_solutions_fine(n,j,k,q) * csqfine(n,j,v) * tran(j) + L_tran * func(j)
                                 L_sum  = L_sum + L_func * wfine(n,j,v)
                              enddo
                              LP_sources_up(n,k,q)  = L_sum * ke + L_extinction(N,q) * Raycon(v) * sum
                           else
                              L_sum = zero
                              do j = 1, nfinedivs(n,v)
                                 L_func = LP_solutions_fine(n,j,k,q) * csqfine(n,j,v) * tran(j)
                                 L_sum  = L_sum + L_func * wfine(n,j,v)
                              enddo
                              LP_sources_up(n,k,q)  = L_sum * ke
                           endif
                        enddo
                     endif
                  enddo
               endif
            endif

!  End layer loop and general Enhanced PS case

         enddo
      endif

!  Source function integration
!  ===========================

!  NLEVEL = Layer index for given optical depth
!  Cumulative source terms : Loop over layers working upwards from NSTART to level NUT,
!  Check for updating the recursion

!  INTENSITY Main loop over all output optical depths
!          Cumulative source term will be saved

      NC = 0 ; M4 = 4.0_fpk * Mu0(v)
      CUMSOURCE_UP(NC) = zero
      CUMSOURCE_DB(NC) = M4 * REFLEC(v) * attenuations(nlayers)
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1
      DO UTA = N_USER_LEVELS, 1, -1
         NUT    = USER_LEVELS(UTA) + 1
         DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            CUMSOURCE_DB(NC) = LOSTRANS_UP(N) * CUMSOURCE_DB(NC-1)
            CUMSOURCE_UP(NC) = SOURCES_UP(N) + LOSTRANS_UP(N) * CUMSOURCE_UP(NC-1)
         ENDDO
         INTENSITY_DB(UTA,V) = FLUX * CUMSOURCE_DB(NC)
         INTENSITY_UP(UTA,V) = FLUX * CUMSOURCE_UP(NC)
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
         NUT_PREV = NUT
      ENDDO

!  Surface WFs

      if ( do_surfacewfs ) then
         Term1 = M4  * attenuations(nlayers)
         do q = 1, n_surfacewfs
            LS_cumsource(q) = term1 * LS_reflec(v,q)
         enddo
         NSTART = NLAYERS
         NUT_PREV = NSTART + 1
         DO UTA = N_USER_LEVELS, 1, -1
            NUT    = USER_LEVELS(UTA) + 1
            DO N = NSTART, NUT, -1
               NC = NLAYERS + 1 - N
               do q = 1, n_surfacewfs
                  LS_cumsource(q) = LOSTRANS_UP(N) * LS_CUMSOURCE(Q)
               enddo
            ENDDO
            do q = 1, n_surfacewfs
               LS_JACOBIANS_DB(UTA,V,Q) = FLUX * LS_CUMSOURCE(Q)
            enddo
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
            NUT_PREV = NUT
         ENDDO
      endif

!  Profile Wfs (atmospheric term)

      if ( do_profilewfs ) then
         do k = 1, nlayers
            if ( Qvary(k) ) then
               L_CUMSOURCE = zero
               NSTART = NLAYERS
               NUT_PREV = NSTART + 1
               DO UTA = N_USER_LEVELS, 1, -1
                  NUT    = USER_LEVELS(UTA) + 1
                  DO N = NSTART, NUT, -1
                     NC = NLAYERS + 1 - N
                     if ( k.eq.n ) then
                        do q = 1, Qnums(k)
                           L_cumsource(q) = LP_SOURCES_UP(N,K,Q)         + &
                                L_LOSTRANS_UP(N,Q) * CUMSOURCE_UP(NC-1) + &
                                  LOSTRANS_UP(N)   * L_CUMSOURCE(Q)
                        enddo
                     else
                        do q = 1, Qnums(k)
                           L_cumsource(q) = LP_SOURCES_UP(N,K,Q)         + &
                                  LOSTRANS_UP(N)   * L_CUMSOURCE(Q)
                        enddo
                     endif
                  ENDDO
                  do q = 1, Qnums(k)
                     LP_JACOBIANS_UP(UTA,V,K,Q) = FLUX * L_CUMSOURCE(Q)
                  enddo
                  IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
                  NUT_PREV = NUT
               ENDDO
            endif
         enddo
      endif

!  Profile Wfs (direct beam term)

      if ( do_profilewfs ) then
         term1 = M4 * reflec(v)
         do k = 1, nlayers
            if ( Qvary(k) ) then
               do q = 1, Qnums(k)
                  L_CUMSOURCE(q) = Term1 * LP_attenuations(nlayers,k,q)
               enddo
               NSTART = NLAYERS
               NUT_PREV = NSTART + 1
               DO UTA = N_USER_LEVELS, 1, -1
                  NUT    = USER_LEVELS(UTA) + 1
                  DO N = NSTART, NUT, -1
                     NC = NLAYERS + 1 - N
                     if ( k.eq.n ) then
                        do q = 1, Qnums(k)
                           L_cumsource(q) =  L_LOSTRANS_UP(N,Q) * CUMSOURCE_DB(NC-1) + &
                                               LOSTRANS_UP(N)   * L_CUMSOURCE(Q)
                        enddo
                     else
                        do q = 1, Qnums(k)
                           L_cumsource(q) = LOSTRANS_UP(N) * L_CUMSOURCE(Q)
                        enddo
                     endif
                  ENDDO
                  do q = 1, Qnums(k)
                     LP_JACOBIANS_DB(UTA,V,K,Q) = FLUX * L_CUMSOURCE(Q)
                  enddo
                  IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
                  NUT_PREV = NUT
               ENDDO
            endif
         enddo
      endif

!  End GEometry Loop

   enddo

!  Finish

   return
end subroutine SS_Integral_ILPS_UP

!

subroutine SS_Integral_ILPS_DN &
   ( maxgeoms, maxlayers, maxfinelayers, maxmoments_input,                              & ! Inputs (dimensioning)
     max_user_levels, max_atmoswfs,                                                     & ! Inputs (dimensioning)
     do_deltam_scaling, do_PlanPar, do_regular_ps, do_enhanced_ps, doNadir,             & ! Inputs (Flags)
     do_profilewfs, Lvaryflags, Lvarynums, Lvarymoms,                                   & ! Inputs (control, Jacobian )
     ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,            & ! Inputs (control,  output)
     extinction, deltaus, omega, truncfac, phasmoms, flux,                              & ! Inputs (Optical - Regular)
     L_extinction, L_deltaus, L_omega, L_truncfac, L_phasmoms,                          & ! Inputs (Optical - Linearized)
     Mu1, LegPoly_dn, radii, NCrit, RadCrit, CotCrit, xfine, wfine, csqfine, cotfine,   & ! Inputs (Geometry)
     Raycon, cota, sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,                  & ! Inputs (Geometry)
     intensity_dn, LP_Jacobians_dn )                                                      ! Output

!  Stand-alone routine for Downwelling Solar-beam Single-scatter (SS)
!    computation of Radiances and LPS Jacobians. Inputs: geometry, spherical functions, optical properties.

!  This version, revised by R. Spurr, 01 June 2012
!   Extension to multiple geometries, 19 December 2012

   implicit none         

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs
!  ======

!  Dimensions

   integer, Intent(in) :: maxgeoms
   integer, Intent(in) :: maxlayers
   integer, Intent(in) :: maxfinelayers
   integer, Intent(in) :: maxmoments_input
   integer, Intent(in) :: max_user_levels
   INTEGER, Intent(in) :: max_atmoswfs

!  flags

   logical, Intent(in) :: DO_DELTAM_SCALING
   logical, Intent(in) :: DO_PLANPAR
   logical, Intent(in) :: DO_REGULAR_PS
   logical, Intent(in) :: DO_ENHANCED_PS
   logical, Intent(in) :: DONADIR(MAXGEOMS)

!  Jacobian Flag

   LOGICAL, Intent(in) :: do_profilewfs

!  Layer and Level Control Numbers, Number of Moments

   integer, Intent(in) :: NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   integer, Intent(in) :: NGEOMS, NMOMENTS_INPUT

   integer, Intent(in) :: N_USER_LEVELS
   integer, Intent(in) :: USER_LEVELS ( MAX_USER_LEVELS )

!  Jacobian control

   LOGICAL, Intent(in) :: Lvaryflags(maxlayers)
   INTEGER, Intent(in) :: Lvarynums (maxlayers)
   LOGICAL, Intent(in) :: Lvarymoms (maxlayers,max_atmoswfs)

!  optical inputs
!  --------------

!  Atmosphere

   real(fpk), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(fpk), Intent(in) :: OMEGA       ( MAXLAYERS )
   real(fpk), Intent(in) :: TRUNCFAC    ( MAXLAYERS )
   real(fpk), Intent(in) :: PHASMOMS    ( MAXLAYERS,0:MAXMOMENTS_INPUT )

!  Solar Flux

   real(fpk), Intent(in) :: FLUX

!  Linearized optical inputs

   real(fpk), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
   real(fpk), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
   real(fpk), Intent(in) :: L_OMEGA       ( MAXLAYERS, max_atmoswfs )
   real(fpk), Intent(in) :: L_TRUNCFAC    ( MAXLAYERS, max_atmoswfs )
   real(fpk), Intent(in) :: L_PHASMOMS    ( MAXLAYERS,0:MAXMOMENTS_INPUT, max_atmoswfs )

!  Geometrical inputs
!  ------------------

!       Ray constant, Cotangents
!       Mu1 = cos(alpha_boa), required for the Regular PS only
!       solar paths, Legendre Polynomials

   integer  , Intent(in)  :: NCrit(maxgeoms)
   real(fpk), Intent(in)  :: RadCrit(maxgeoms), CotCrit(maxgeoms)
   real(fpk), Intent(in)  :: Raycon(maxgeoms), cota(0:maxlayers,maxgeoms), radii(0:maxlayers)
   real(fpk), Intent(in)  :: Mu1(maxgeoms)

!  solar paths 

   integer  , Intent(in)  :: ntraverse  (0:maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: sunpaths   (0:maxlayers,maxlayers,maxgeoms)
   integer  , Intent(in)  :: ntraverse_fine(maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: sunpaths_fine (maxlayers,maxlayers,maxfinelayers,maxgeoms)

!  Legendres

   real(fpk), Intent(in)  :: LegPoly_dn(0:maxmoments_input,maxgeoms)

!  LOS Quadratures for Enhanced PS

   real(fpk), Intent(in)  :: xfine   (maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: wfine   (maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: csqfine (maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: cotfine (maxlayers,maxfinelayers,maxgeoms)

!  outputs
!  -------

   real(fpk), Intent(Out)  :: intensity_dn     ( max_user_levels, maxgeoms )
   real(fpk), Intent(Out)  :: LP_Jacobians_dn  ( max_user_levels, maxgeoms, maxlayers, max_atmoswfs )

!  LOCAL
!  -----

!  Attenuations

   real(fpk)  :: attenuations      (0:maxlayers)
   real(fpk)  :: LP_attenuations   (0:maxlayers,maxlayers,max_atmoswfs)

!  Solutions for Enhaced-PS

   real(fpk)  :: Solutions_fine    (maxlayers,maxfinelayers)
   real(fpk)  :: LP_solutions_fine (maxlayers,maxfinelayers,maxlayers,max_atmoswfs)

!  Scattering

   real(fpk)  :: tms            (maxlayers)
   real(fpk)  :: exactscat_dn   (maxlayers)
   real(fpk)  :: L_tms          (maxlayers,max_atmoswfs)
   real(fpk)  :: L_exactscat_dn (maxlayers,max_atmoswfs)

!  Source function integration results

   real(fpk) :: sources_dn         ( maxlayers )
   real(fpk) :: LP_sources_dn      ( maxlayers,maxlayers,max_atmoswfs )

   real(fpk) :: lostrans_dn        ( maxlayers )
   real(fpk) :: L_lostrans_dn      ( maxlayers,max_atmoswfs )

   real(fpk) :: cumsource_dn      ( 0:maxlayers )
   real(fpk) :: L_cumsource       ( max_atmoswfs )

!  Regular_PS or plane-parallel flag

   logical    :: do_RegPSorPP

!  Help

   integer    :: n, k, j, q, L, v, uta, nstart, nc, nut, nut_prev, Qnums(maxlayers)
   logical    :: layermask_dn(maxlayers), Qvary(maxlayers)

   real(fpk)  :: argum(maxfinelayers), tran(maxfinelayers), func(maxfinelayers)
   real(fpk)  :: cons, help, sum, tran_1, kn, ke, factor1, factor2
   real(fpk)  :: L_help, L_sum, L_tran, L_func, L_factor1, L_factor2
   real(fpk)  :: cot_1, cot_2, multiplier, suntau(0:maxlayers), lostau, rdiff, cot_c
   real(fpk)  :: L_multiplier, LP_suntau(0:maxlayers,maxlayers,max_atmoswfs), L_lostau(max_atmoswfs)
   real(fpk)  :: attenuations_fine, L_attenuations_fine, sumd, L_sumd, L_exactscat

   real(fpk)  :: consc, trand

      real(fpk), parameter  :: cutoff = 88.0_fpk
      real(fpk), parameter  :: zero   = 0.0_fpk
      real(fpk), parameter  :: one    = 1.0_fpk

!  Zero the output

   INTENSITY_DN    = zero
   LP_JACOBIANS_DN = zero

!  Regular_PS or plane-parallel flag

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

      lostrans_dn   = zero  ; sources_dn    = zero ; cumsource_dn = zero
      L_lostrans_dn = zero  ; LP_sources_dn = zero

!  Scattering functions and Linearization
!  ======================================

      do n = 1, nlayers
         if ( layermask_dn(n) ) then
            sum = zero
            do L = 0, nmoments_input
               sum = sum + LegPoly_dn(L,v) * phasmoms(n,L)
            enddo
            exactscat_dn(n) = sum * tms(n)
            if ( Qvary(n) ) then
               do q = 1, Qnums(n)
                  if ( Lvarymoms(n,q) ) then
                     L_sum = zero
                     do L = 0, nmoments_input
                        L_sum = L_sum + LegPoly_dn(L,v) * L_phasmoms(n,L,q)
                     enddo
                     L_exactscat_dn(n,q) = L_sum * tms(n) + sum * L_tms(n,q)
                  else
                     L_exactscat_dn(n,q) = sum * L_tms(n,q)
                  endif
               enddo
            endif               
         endif
      enddo

!  Attenuations and Solar solutions
!  ================================

!  Initialize, only to layer Ncrit if applicable

      Attenuations = ZERO   ; LP_Attenuations = ZERO
      Suntau       = ZERO   ; LP_suntau       = ZERO
      Solutions_fine = zero ; LP_Solutions_fine = zero

      nstart = nlayers ; if (Ncrit(v).ne.0) nstart = nCrit(v)
      consc = zero ; trand = one

!  Attenuations to End points (including TOA). All representations
!  ===============================================================

      do n = 0, nlayers
         sumd = ZERO
         do k = 1, ntraverse(n,v)
            sumd = sumd + extinction(k) * sunpaths(n,k,v)
         enddo
         suntau(n) = sumd
         If (sumd .lt. cutoff ) Attenuations(n) = exp( - sumd )
         if ( do_profilewfs ) then
            do k = 1, nlayers
               if ( Qvary(k) .and. k.le.ntraverse(n,v) ) then
                  do q = 1, Qnums(k)
                     LP_suntau(n,k,q) = L_extinction(k,q) * sunpaths(n,k,v)
                     LP_Attenuations(n,k,q) = - Attenuations(n) * LP_suntau(n,k,q)
                  enddo
               endif
            enddo
         endif
      enddo

!  Enhanced-spherical, fine-layer attenuations
!  ===========================================

      if ( do_enhanced_ps ) then
         do n = 1, nstart
            if ( layermask_dn(n) ) then
               do j = 1, nfinedivs(n,v)
                  sumd = ZERO
                  do k = 1, ntraverse_fine(n,j,v)
                     sumd = sumd + extinction(k) * sunpaths_fine(n,k,j,v)
                  enddo
                  if (sumd .lt. cutoff ) Attenuations_fine = exp( - sumd )
                  Solutions_fine(n,j) = exactscat_dn(n) * Attenuations_fine
                  if ( do_profilewfs ) then
                     do k = 1, nlayers
                        if ( Qvary(k) .and. k.le.ntraverse_fine(n,j,v) ) then
                           do q = 1, Qnums(k)
                               L_exactscat = zero ; if ( k.eq.n) L_exactscat = L_exactscat_dn(n,q)
                               L_sumd = L_extinction(k,q) * sunpaths_fine(n,k,j,v)
                               L_Attenuations_fine = - Attenuations_fine * L_sumd 
                               LP_Solutions_fine(n,j,k,q) = L_exactscat       *   Attenuations_fine   + &
                                                              exactscat_dn(n) * L_Attenuations_fine
                           enddo
                        endif
                     enddo
                  endif
               enddo
            endif
         enddo
      endif

!  Plane-Parallel or Regular-PS (Average secant formulation): Layer integrated Solar sources
!  =========================================================================================

      if ( do_RegPSorPP ) then
         do n = nlayers, 1, -1

!  LOS transmittance (Not for the Horizontal View)

            if ( Mu1(v) .gt. zero ) then
               lostau = deltaus(n) / Mu1(v)
               if ( lostau .lt. cutoff ) lostrans_dn(n) = exp( - lostau )
               if ( do_profilewfs ) then
                  if ( Qvary(n) ) then
                     do q = 1, Qnums(n)
                        L_lostau(q)        = L_deltaus(n,q)  / Mu1(v)
                        L_lostrans_dn(n,q) = - L_lostau(q) * lostrans_dn(n)
                     enddo
                  endif
               endif
            endif

!  Sources, general case

            if ( layermask_dn(n) .and. n.le.nstart  ) then
              if ( Mu1(v) .gt. zero ) then
                factor1 = Attenuations(n-1)*lostrans_dn(n) - Attenuations(n)
                factor2 = ((suntau(n) - suntau(n-1))/lostau) - one
                multiplier = factor1 / factor2
                sources_dn(n) = exactscat_dn(n) * multiplier
                if ( do_profilewfs ) then
                  do k = 1, nlayers
                     if ( Qvary(k).and.k.le.ntraverse(n,v) ) then
                        do q = 1, Qnums(k)
                           L_factor1 = LP_Attenuations(n-1,k,q)*lostrans_dn(n) - LP_Attenuations(n,k,q)
                           L_factor2 = (LP_suntau(n,k,q) - LP_suntau(n-1,k,q))/lostau
                           if ( k.eq.n) then
                              L_factor1 = L_factor1 + Attenuations(n-1)*L_lostrans_dn(n,q)
                              L_factor2 = L_factor2 - (factor2 + one)*L_lostau(q)/lostau 
                              L_multiplier = ( L_factor1 - multiplier*L_factor2 ) / factor2
                              LP_sources_dn(n,k,q) = L_exactscat_dn(n,q) *   multiplier + &
                                                       exactscat_dn(n)   * L_multiplier
                           else
                              L_multiplier = ( L_factor1 - multiplier*L_factor2 ) / factor2
                              LP_sources_dn(n,k,q) = exactscat_dn(n)   * L_multiplier
                           endif
                        enddo
                     endif
                  enddo
                endif
              endif
            endif

!  Sources, Special case

            if ( layermask_dn(n) .and. n.le.nstart  ) then
              if ( Mu1(v) .eq. zero ) then
                factor1 = Attenuations(n)
                sources_dn(n) = exactscat_dn(n) * factor1
                if ( do_profilewfs ) then
                  do k = 1, nlayers
                     if ( Qvary(k).and.k.le.ntraverse(n,v) ) then
                        do q = 1, Qnums(k)
                           L_factor1 = LP_Attenuations(n,k,q)
                           if ( k.eq.n) then
                              LP_sources_dn(n,k,q) = L_exactscat_dn(n,q) *   factor1 + &
                                                       exactscat_dn(n)   * L_factor1
                           else
                              LP_sources_dn(n,k,q) = exactscat_dn(n)   * L_factor1
                           endif
                        enddo
                     endif
                  enddo
                endif
              endif
            endif

!  End layers and regular-PS formulation

         enddo
      endif

!  Enhanced PS: special case (nadir viewing). Layer integrated Solar sources
!  =========================================================================

      if ( do_enhanced_ps .and. doNadir(v) ) then
         do n = nlayers, 1, -1

!  LOS transmittance

            kn     = extinction(n)
            lostau = deltaus(n)
            if ( lostau .lt. cutoff ) lostrans_dn(n) = exp( - lostau )
            if ( do_profilewfs ) then
               if ( Qvary(n) ) then
                  do q = 1, Qnums(n)
                     L_lostau(q)        = L_deltaus(n,q)
                     L_lostrans_dn(n,q) = - L_lostau(q) * lostrans_dn(n)
                  enddo
               endif
            endif

!  Sources

            if ( layermask_dn(n) .and. n.le.nstart  ) then
               rdiff = radii(n-1) - radii(n)
               if ( n.eq.NCrit(v)) rdiff = radii(n-1) - RadCrit(v)
               sum = zero
               do j = 1, nfinedivs(n,v)
                  argum(j) = rdiff - xfine(n,j,v)
                  tran(j)  = exp ( - argum(j) * kn )
                  func(j)  = solutions_fine(n,j) * tran(j)
                  sum = sum + func(j) * wfine(n,j,v)
               enddo
               sources_dn(n) = sum * kn
               if ( n.eq.NCrit(v) ) then
                  trand = exp ( - kn * rdiff )
                  sources_dn(n) = sources_dn(n) * trand
               endif
               if ( do_profilewfs ) then
                  do k = 1, nlayers
                     if ( Qvary(k).and.k.le.ntraverse(n,v) ) then
                        do q = 1, Qnums(k)
                           if ( k.eq.n ) then
                              L_sum = zero
                              do j = 1, nfinedivs(n,v)
                                 L_tran = - argum(j) * L_extinction(n,q)
                                 L_func = LP_solutions_fine(n,j,k,q) * tran(j) + L_tran * func(j)
                                 L_sum  = L_sum + L_func * wfine(n,j,v)
                              enddo
                              LP_sources_dn(n,k,q)  = L_sum * kn + L_extinction(N,q) * sum
                              if ( n.eq.NCrit(v) ) then
                                 LP_sources_dn(n,k,q) =  LP_sources_dn(n,k,q) * trand - &
                                           sources_dn(n)* L_extinction(N,q) * rdiff
                              endif
                           else
                              L_sum = zero
                              do j = 1, nfinedivs(n,v)
                                 L_func = LP_solutions_fine(n,j,k,q)  * tran(j)
                                 L_sum  = L_sum + L_func * wfine(n,j,v)
                              enddo
                              LP_sources_dn(n,k,q)  = L_sum * kn
                              if ( n.eq.NCrit(v) ) LP_sources_dn(n,k,q) = LP_sources_dn(n,k,q) * trand
                           endif
                        enddo
                     endif
                  enddo
               endif
            endif

!  End layer loop and Nadir Enhanced PS case

         enddo
      endif

!  Enhanced PS: General case. Layer integrated Solar sources
!  =========================================================

      if ( do_enhanced_ps .and. .not. doNadir(v) ) then
         do n = nlayers, 1, -1

!  LOS transmittance

            cot_2 = cota(n-1,v) ; cot_1 = cota(n,v)
            cot_c = cot_1  ; if ( n.eq.NCrit(v) ) cot_c = CotCrit(v)
            kn = extinction(n) ;  ke = raycon(v) * kn ; cons = raycon(v) * ( cot_2 - cot_1 )
            tran_1 = kn * cons
            if ( tran_1 .lt. cutoff ) lostrans_dn(n) = exp ( - tran_1 )
            if ( n.eq.NCrit(v) ) then
               consc = raycon(v) * ( CotCrit(v) - cot_1 )
               trand = exp ( - kn * consc )
            endif
            if ( do_profilewfs ) then
               if ( Qvary(n) ) then
                  do q = 1, Qnums(n)
                     L_lostau(q)        = L_extinction(n,q) * cons
                     L_lostrans_dn(n,q) = - L_lostau(q) * lostrans_dn(n)
                  enddo
               endif
            endif

!  Sources

            if ( layermask_dn(n) .and. n.le.nstart  ) then
               sum = zero
               do j = 1, nfinedivs(n,v)
                  argum(j) = Raycon(v) * ( cotfine(n,j,v) - cot_c )
                  tran(j)  = exp ( - kn * argum(j) )
                  func(j)  = solutions_fine(n,j) * csqfine(n,j,v) * tran(j)
                  sum      = sum + func(j) * wfine(n,j,v)
               enddo
               sources_dn(n) = sum * ke
               if ( n.eq.NCrit(v) ) sources_dn(n) = sources_dn(n) * trand
               if ( do_profilewfs ) then
                  do k = 1, nlayers
                     if ( Qvary(k).and.k.le.ntraverse(n,v) ) then
                        do q = 1, Qnums(k)
                           if ( k.eq.n ) then
                              L_sum = zero
                              do j = 1, nfinedivs(n,v)
                                 L_tran = - argum(j) * L_extinction(n,q)
                                 L_func = LP_solutions_fine(n,j,k,q) * csqfine(n,j,v) * tran(j) + L_tran * func(j)
                                 L_sum  = L_sum + L_func * wfine(n,j,v)
                              enddo
                              LP_sources_dn(n,k,q)  = L_sum * ke + L_extinction(N,q) * Raycon(v) * sum
                              if ( n.eq.NCrit(v) ) then
                                 LP_sources_dn(n,k,q) =  LP_sources_dn(n,k,q) * trand - &
                                           sources_dn(n)* L_extinction(N,q) * consc
                              endif
                           else
                              L_sum = zero
                              do j = 1, nfinedivs(n,v)
                                 L_func = LP_solutions_fine(n,j,k,q) * csqfine(n,j,v) * tran(j)
                                 L_sum  = L_sum + L_func * wfine(n,j,v)
                              enddo
                              LP_sources_dn(n,k,q)  = L_sum * ke
                              if ( n.eq.NCrit(v) ) LP_sources_dn(n,k,q) = LP_sources_dn(n,k,q) * trand
                           endif
                        enddo
                     endif
                  enddo
               endif
            endif

!  End layer loop and general Enhanced PS case

         enddo
      endif

!  Source function integration
!  ===========================

!  NLEVEL = Layer index for given optical depth
!  Cumulative source terms : Loop over layers working upwards from NSTART to level NUT,
!  Check for updating the recursion

!  INTENSITY Main loop over all output optical depths
!          Cumulative source term will be saved

      NC = 0
      CUMSOURCE_DN(NC) = zero
      NSTART = 1
      NUT_PREV = NSTART - 1
      DO UTA = 1, N_USER_LEVELS
         NUT    = USER_LEVELS(UTA)
         DO N = NSTART, NUT
            NC = N
            CUMSOURCE_DN(NC) = SOURCES_DN(N) + LOSTRANS_DN(N) * CUMSOURCE_DN(NC-1)
         ENDDO
         INTENSITY_DN(UTA,V) = FLUX * CUMSOURCE_DN(NC)
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
         NUT_PREV = NUT
      ENDDO

!  Profile Wfs

      if ( do_profilewfs ) then
         do k = 1, nlayers
            if ( Qvary(k) ) then
               do q = 1, Qnums(k)
                  L_CUMSOURCE(q) = ZERO
               enddo
               NSTART = 1
               NUT_PREV = NSTART - 1
               DO UTA = 1, N_USER_LEVELs
                  NUT    = USER_LEVELS(UTA)
                  DO N = NSTART, NUT
                     NC = N
                     if ( k.eq.n ) then
                        do q = 1, Qnums(k)
                           L_cumsource(q) = LP_SOURCES_DN(N,K,Q)         + &
                                L_LOSTRANS_DN(N,Q) * CUMSOURCE_DN(NC-1) + &
                                  LOSTRANS_DN(N)   * L_CUMSOURCE(Q)
                        enddo
                     else
                        do q = 1, Qnums(k)
                           L_cumsource(q) = LP_SOURCES_DN(N,K,Q)         + &
                                  LOSTRANS_DN(N)   * L_CUMSOURCE(Q)
                        enddo
                     endif
                  ENDDO
                  do q = 1, Qnums(k)
                     LP_JACOBIANS_DN(UTA,V,K,Q) = FLUX * L_CUMSOURCE(Q)
                  enddo
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
end subroutine SS_Integral_ILPS_DN

!

subroutine SS_Integral_ILPS_UPDN &
   ( maxgeoms, maxlayers, maxfinelayers, maxmoments_input,                             & ! Inputs (dimensioning)
     max_user_levels, max_atmoswfs, max_surfacewfs,                                    & ! Inputs (dimensioning)
     do_upwelling, do_dnwelling,                                                       & ! Inputs (Flags)
     do_deltam_scaling, do_PlanPar, do_regular_ps, do_enhanced_ps, doNadir,            & ! Inputs (Flags)
     do_profilewfs, do_surfacewfs, n_surfacewfs, Lvaryflags, Lvarynums, Lvarymoms,     & ! Inputs (control, Jacobian )
     ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,           & ! Inputs (control,  output)
     reflec, extinction, deltaus, omega, truncfac, phasmoms, flux,                     & ! Inputs (Optical - Regular)
     LS_reflec, L_extinction, L_deltaus, L_omega, L_truncfac, L_phasmoms,              & ! Inputs (Optical - Linearized)
     Mu0, Mu1, LegPoly_up, LegPoly_dn, radii, NCrit, RadCrit, CotCrit,                 & ! Inputs (Geometry)
     xfine, wfine, csqfine, cotfine, Raycon, cota,                                     & ! Inputs (Geometry)
     sunpaths_up, ntraverse_up, sunpaths_dn, ntraverse_dn,                             & ! Inputs (Whole layer Geometry)
     sunpaths_fine_up, ntraverse_fine_up, sunpaths_fine_dn, ntraverse_fine_dn,         & ! Inputs (Fine  layer Geometry)
     intensity_up, intensity_db, LP_Jacobians_up, LP_Jacobians_db, LS_Jacobians_db,    & ! Output (upwelling)
     intensity_dn, LP_Jacobians_dn )                                                     ! Output (downwelling)

   implicit none

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  Dimensions

   integer, Intent(in) :: maxgeoms
   integer, Intent(in) :: maxlayers
   integer, Intent(in) :: maxfinelayers
   integer, Intent(in) :: maxmoments_input
   integer, Intent(in) :: max_user_levels
   INTEGER, Intent(in) :: max_atmoswfs
   INTEGER, Intent(in) :: max_surfacewfs

!  General flags

   LOGICAL, Intent(in) :: DO_UPWELLING
   LOGICAL, Intent(in) :: DO_DNWELLING

!  flags

   logical, Intent(in) :: DO_DELTAM_SCALING
   logical, Intent(in) :: DO_PLANPAR
   logical, Intent(in) :: DO_REGULAR_PS
   logical, Intent(in) :: DO_ENHANCED_PS
   logical, Intent(in) :: DONADIR(MAXGEOMS)

!  Jacobian Flags

   LOGICAL, Intent(in) :: do_surfacewfs
   LOGICAL, Intent(in) :: do_profilewfs

!  Layer and Level Control Numbers, Number of Moments

   integer, Intent(in) :: NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   integer, Intent(in) :: NGEOMS, NMOMENTS_INPUT

   integer, Intent(in) :: N_USER_LEVELS
   integer, Intent(in) :: USER_LEVELS ( MAX_USER_LEVELS )

!  Jacobian control

   INTEGER, Intent(in) :: n_surfacewfs
   LOGICAL, Intent(in) :: Lvaryflags(maxlayers)
   INTEGER, Intent(in) :: Lvarynums (maxlayers)
   LOGICAL, Intent(in) :: Lvarymoms (maxlayers,max_atmoswfs)

!  optical inputs
!  --------------

!  Atmosphere

   real(fpk), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(fpk), Intent(in) :: OMEGA       ( MAXLAYERS )
   real(fpk), Intent(in) :: TRUNCFAC    ( MAXLAYERS )
   real(fpk), Intent(in) :: PHASMOMS    ( MAXLAYERS,0:MAXMOMENTS_INPUT )

!  Solar Flux and Surface reflectivity (Could be the albedo)

   real(fpk), Intent(in) :: REFLEC(MAXGEOMS), FLUX

!  Linearized optical inputs

   real(fpk), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
   real(fpk), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
   real(fpk), Intent(in) :: L_OMEGA       ( MAXLAYERS, max_atmoswfs )
   real(fpk), Intent(in) :: L_TRUNCFAC    ( MAXLAYERS, max_atmoswfs )
   real(fpk), Intent(in) :: L_PHASMOMS    ( MAXLAYERS,0:MAXMOMENTS_INPUT, max_atmoswfs )
   real(fpk), Intent(in) :: LS_REFLEC     ( MAXGEOMS, max_surfacewfs)

!  Geometrical inputs
!  ------------------

!       Ray constant, Cotangents
!       Mu0 = cos(theta_boa), required for the surface term (both regular and enhanced)
!       Mu1 = cos(alpha_boa), required for the Regular PS only
!       solar paths, Legendre Polynomials

   integer  , Intent(in)  :: NCrit(maxgeoms)
   real(fpk), Intent(in)  :: RadCrit(maxgeoms), CotCrit(maxgeoms)
   real(fpk), Intent(in)  :: Raycon(maxgeoms), cota(0:maxlayers,maxgeoms), radii(0:maxlayers)
   real(fpk), Intent(in)  :: Mu0(maxgeoms), Mu1(maxgeoms)

!  solar paths 

   integer  , Intent(in)  :: ntraverse_up  (0:maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: sunpaths_up   (0:maxlayers,maxlayers,maxgeoms)
   integer  , Intent(in)  :: ntraverse_fine_up(maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: sunpaths_fine_up (maxlayers,maxlayers,maxfinelayers,maxgeoms)

   integer  , Intent(in)  :: ntraverse_dn  (0:maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: sunpaths_dn   (0:maxlayers,maxlayers,maxgeoms)
   integer  , Intent(in)  :: ntraverse_fine_dn(maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: sunpaths_fine_dn (maxlayers,maxlayers,maxfinelayers,maxgeoms)

!  Legendres

   real(fpk), Intent(in)  :: LegPoly_up(0:maxmoments_input,maxgeoms)
   real(fpk), Intent(in)  :: LegPoly_dn(0:maxmoments_input,maxgeoms)

!  LOS Quadratures for Enhanced PS

   real(fpk), Intent(in)  :: xfine   (maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: wfine   (maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: csqfine (maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: cotfine (maxlayers,maxfinelayers,maxgeoms)

!  outputs
!  -------

   real(fpk), Intent(Out)  :: intensity_up     ( max_user_levels, maxgeoms )
   real(fpk), Intent(Out)  :: intensity_db     ( max_user_levels, maxgeoms )
   real(fpk), Intent(Out)  :: LP_Jacobians_up  ( max_user_levels, maxgeoms, maxlayers, max_atmoswfs )
   real(fpk), Intent(Out)  :: LP_Jacobians_db  ( max_user_levels, maxgeoms, maxlayers, max_atmoswfs )
   real(fpk), Intent(Out)  :: LS_Jacobians_db  ( max_user_levels, maxgeoms, max_surfacewfs )

   real(fpk), Intent(Out)  :: intensity_dn     ( max_user_levels, maxgeoms )
   real(fpk), Intent(Out)  :: LP_Jacobians_dn  ( max_user_levels, maxgeoms, maxlayers, max_atmoswfs )

!  Upwelling

   if ( do_upwelling  ) then
      call SS_Integral_ILPS_UP &
   ( maxgeoms, maxlayers, maxfinelayers, maxmoments_input,                             & ! Inputs (dimensioning)
     max_user_levels, max_atmoswfs, max_surfacewfs,                                    & ! Inputs (dimensioning)
     do_deltam_scaling, do_PlanPar, do_regular_ps, do_enhanced_ps, doNadir,            & ! Inputs (Flags - General  )
     do_profilewfs, do_surfacewfs, n_surfacewfs, Lvaryflags, Lvarynums, Lvarymoms,     & ! Inputs (Control - Jacobian )
     ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,           & ! Inputs (control - output)
     reflec, extinction, deltaus, omega, truncfac, phasmoms, flux,                     & ! Inputs (Optical - Regular)
     LS_reflec, L_extinction, L_deltaus, L_omega, L_truncfac, L_phasmoms,              & ! Inputs (Optical - Linearized)
     Mu0, Mu1, LegPoly_up, NCrit, xfine, wfine, csqfine, cotfine,                      & ! Inputs (Geometry)
     Raycon, cota, sunpaths_up, ntraverse_up, sunpaths_fine_up, ntraverse_fine_up,     & ! Inputs (Geometry)
     intensity_up, intensity_db, LP_Jacobians_up, LP_Jacobians_db, LS_Jacobians_db )     ! Output

   endif

!  Downwelling

   if ( do_dnwelling  ) then
      call SS_Integral_ILPS_DN &
   ( maxgeoms, maxlayers, maxfinelayers, maxmoments_input,                              & ! Inputs (dimensioning)
     max_user_levels, max_atmoswfs,                                                     & ! Inputs (dimensioning)
     do_deltam_scaling, do_PlanPar, do_regular_ps, do_enhanced_ps, doNadir,             & ! Inputs (Flags)
     do_profilewfs, Lvaryflags, Lvarynums, Lvarymoms,                                   & ! Inputs (control, Jacobian )
     ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,            & ! Inputs (control,  output)
     extinction, deltaus, omega, truncfac, phasmoms, flux,                              & ! Inputs (Optical - Regular)
     L_extinction, L_deltaus, L_omega, L_truncfac, L_phasmoms,                          & ! Inputs (Optical - Linearized)
     Mu1, LegPoly_dn, radii, NCrit, RadCrit, CotCrit, xfine, wfine, csqfine, cotfine,   & ! Inputs (Geometry)
     Raycon, cota, sunpaths_dn, ntraverse_dn, sunpaths_fine_dn, ntraverse_fine_dn,      & ! Inputs (Geometry)
     intensity_dn, LP_Jacobians_dn )                                                      ! Output
   endif

!  Finish

   return
end subroutine SS_Integral_ILPS_UPDN

!  End module

end module FO_ScalarSS_RTCalcs_ILPS_m


