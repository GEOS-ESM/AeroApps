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
! #     FIRST-ORDER VECTOR MODEL (EXACT SINGLE SCATTERING)  #
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

module FO_VectorSS_RTCalcs_ILCS_m

!  For a given wavelength, this routine will calculate upwelling and downwelling
!  First Order Stokes vectors, and any number of LCS Jacobians (column/surface)

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
!       SSV_Integral_ILCS_UP   (Upwelling only)
!       SSV_Integral_ILCS_DN   (Downwelling only)
!       SSV_Integral_ILCS_UPDN (Upwelling and Downwelling)

!  All subroutines public

public

contains

subroutine SSV_Integral_ILCS_UP &
   ( maxgeoms, maxlayers, maxfinelayers, maxmoments_input, & ! Inputs (dimensioning)
     max_user_levels, max_atmoswfs, max_surfacewfs,        & ! Inputs (dimensioning)
     do_sunlight, do_deltam_scaling, do_lambertian,        & ! Inputs (Flags - General)
     do_PlanPar, do_regular_ps, do_enhanced_ps, doNadir, nstokes,               & ! Inputs (Flags - General)
     do_columnwfs, do_surfacewfs, n_columnwfs, n_surfacewfs, Lvarymoms,         & ! Inputs (control, Jacobian )
     ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,    & ! Inputs (control,  output)
     reflec, extinction, deltaus, omega, truncfac, greekmat, flux, fluxvec,     & ! Inputs (Optical)
     LS_reflec, L_extinction, L_deltaus, L_omega, L_truncfac, L_greekmat,       & ! Inputs (Optical - Linearized)
     Mu0, Mu1, GenSpher, Rotations, NCrit, xfine, wfine, csqfine, cotfine,      & ! Inputs (Geometry)
     Raycon, cota, sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,          & ! Inputs (Geometry)
     Stokes_up, Stokes_db, LC_Jacobians_up, LC_Jacobians_db, LS_Jacobians_db )    ! Output

!  Stand-alone routine for Upwelling Solar-beam Single-scatter (SS)
!    computation of Stokes vectors and LCS Jacobians. Inputs: geometry, spherical functions, optical properties.

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

   LOGICAL, Intent(in) :: DO_SUNLIGHT
   LOGICAL, Intent(in) :: DO_DELTAM_SCALING
   LOGICAL, Intent(in) :: DO_LAMBERTIAN

   logical, Intent(in) :: DO_PLANPAR
   logical, Intent(in) :: DO_REGULAR_PS
   logical, Intent(in) :: DO_ENHANCED_PS
   logical, Intent(in) :: DONADIR(MAXGEOMS)

!  Jacobian Flags

   LOGICAL, Intent(in) :: do_surfacewfs
   LOGICAL, Intent(in) :: do_columnwfs

!  Layer and Level Control Numbers, Number of Moments

   integer, Intent(in) :: NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   integer, Intent(in) :: NGEOMS, NMOMENTS_INPUT, NSTOKES

   integer, Intent(in) :: N_USER_LEVELS
   integer, Intent(in) :: USER_LEVELS ( MAX_USER_LEVELS )

!  Jacobian control

   INTEGER, Intent(in) :: n_columnwfs
   INTEGER, Intent(in) :: n_surfacewfs
   LOGICAL, Intent(in) :: Lvarymoms (maxlayers,max_atmoswfs)

!  optical inputs
!  --------------

!  Atmosphere

   real(fpk), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(fpk), Intent(in) :: OMEGA       ( MAXLAYERS )
   real(fpk), Intent(in) :: TRUNCFAC    ( MAXLAYERS )
   real(fpk), Intent(in) :: GREEKMAT    ( MAXLAYERS,0:MAXMOMENTS_INPUT, 4, 4 )

!  Solar Flux and Surface reflectivity (Could be the albedo)

   real(fpk), Intent(in) :: REFLEC(4,4,MAXGEOMS), FLUX, FLUXVEC(4)

!  Linearized optical inputs

   real(fpk), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
   real(fpk), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
   real(fpk), Intent(in) :: L_OMEGA       ( MAXLAYERS, max_atmoswfs )
   real(fpk), Intent(in) :: L_TRUNCFAC    ( MAXLAYERS, max_atmoswfs )
   real(fpk), Intent(in) :: L_GREEKMAT    ( MAXLAYERS,0:MAXMOMENTS_INPUT, 4, 4, max_atmoswfs )
   real(fpk), Intent(in) :: LS_REFLEC     ( 4,4,MAXGEOMS, max_surfacewfs)

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

!  Generalized spherical functions.
!    Rotations(1-4)    = C1, S1, C2, S2

   REAL(fpk), Intent(in)  :: GenSpher(0:maxmoments_input,4,maxgeoms)
   REAL(fpk), Intent(in)  :: Rotations(4,maxgeoms)

!  LOS Quadratures for Enhanced PS

   real(fpk), Intent(in)  :: xfine   (maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: wfine   (maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: csqfine (maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: cotfine (maxlayers,maxfinelayers,maxgeoms)

!  outputs
!  -------

   REAL(fpk), Intent(Out)  :: stokes_up     ( max_user_levels, 4, maxgeoms )
   REAL(fpk), Intent(Out)  :: stokes_db     ( max_user_levels, 4, maxgeoms )
   real(fpk), Intent(Out)  :: LC_Jacobians_up  ( max_user_levels, 4, maxgeoms, max_atmoswfs )
   real(fpk), Intent(Out)  :: LC_Jacobians_db  ( max_user_levels, 4, maxgeoms, max_atmoswfs )
   real(fpk), Intent(Out)  :: LS_Jacobians_db  ( max_user_levels, 4, maxgeoms, max_surfacewfs )

!  LOCAL
!  -----

!  Attenuations

   real(fpk)  :: suntau              (0:maxlayers)
   real(fpk)  :: attenuations        (0:maxlayers)
   real(fpk)  :: attenuations_fine   (maxlayers,maxfinelayers)

   real(fpk)  :: LC_suntau            (0:maxlayers,max_atmoswfs)
   real(fpk)  :: LC_attenuations      (0:maxlayers,max_atmoswfs)
   real(fpk)  :: LC_attenuations_fine (maxlayers,maxfinelayers,max_atmoswfs)

!  Scattering

   real(fpk)  :: tms            (maxlayers)
   real(fpk)  :: exactscat_up   (maxlayers,4,4)

   real(fpk)  :: L_tms          (maxlayers,max_atmoswfs)
   real(fpk)  :: L_exactscat_up (maxlayers,4,4, max_atmoswfs)

!  Source function integration results

   real(fpk)  :: sources_up       ( maxlayers, 4 )
   real(fpk)  :: lostrans_up      ( maxlayers )
   real(fpk)  :: multiplier_up    ( maxlayers )

   real(fpk)  :: LC_sources_up    ( maxlayers, 4, max_atmoswfs )
   real(fpk)  :: L_lostrans_up    ( maxlayers, max_atmoswfs )
   real(fpk)  :: L_multiplier_up  ( maxlayers, max_atmoswfs )

!  Local cumulative source terms

   real(fpk)  :: cumsource_db      ( 0:maxlayers, 4 )
   real(fpk)  :: cumsource_up      ( 0:maxlayers, 4 )

   real(fpk)  :: L_cumsource       ( 4, max_atmoswfs )
   real(fpk)  :: LS_cumsource      ( 4, max_surfacewfs )

!  Regular_PS or plane-parallel flag

   logical    :: do_RegPSorPP

!  Help

   integer    :: n, ns, k, j, q, L, v, o1, uta, nstart, nc, nut, nut_prev, Qnums(maxlayers)
   logical    :: layermask_up(maxlayers), Qvary(maxlayers)

   real(fpk)  :: argum(maxfinelayers), tran(maxfinelayers), func(maxfinelayers)
   real(fpk)  :: cons, help, sum, tran_1, kn, ke, factor1, factor2, factor3, m4, m4a, rhelp(4), shelp(4)
   real(fpk)  :: cot_1, cot_2, L_help, L_sum, L_tran, L_func, L_factor1, L_factor2, sumd, L_sumd
   real(fpk)  :: lostau, L_lostau(max_atmoswfs), fmat(6), L_fmat(6), sum23, dif23, L_Shelp
   real(fpk)  :: help3c1, help3s1, help4c1, help4s1

   real(fpk), parameter  :: cutoff = 88.0_fpk
   real(fpk), parameter  :: zero   = 0.0_fpk
   real(fpk), parameter  :: one    = 1.0_fpk

!  Zero the output

   STOKES_UP       = zero
   LC_JACOBIANS_UP = zero

   STOKES_DB       = zero
   LC_JACOBIANS_DB = zero
   LS_JACOBIANS_DB = zero

!  Regular_PS or plane-parallel flag

   do_RegPSorPP = (do_regular_ps .or. do_PlanPar)

!  Bookkeeping

   ns = nstokes ; fmat = zero
   NUT = USER_LEVELS(1) + 1
   LAYERMASK_UP = .false.
   LAYERMASK_UP(NUT:NLAYERS) = .true.

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

      lostrans_up   = zero  ; sources_up   = zero  ; exactscat_up   = zero ; multiplier_up   = zero
      L_lostrans_up = zero  ; LC_sources_up = zero ; L_exactscat_up = zero ; L_multiplier_up = zero

!  Scattering functions and Linearization
!  --------------------------------------

!  Scalar only

      if ( nstokes .eq. 1 ) then
         do n = 1, nlayers
            if ( layermask_up(n) ) then
              sum = zero
               do L = 0, nmoments_input
                  sum = sum + GenSpher(L,1,V) * Greekmat(n,L,1,1)
               enddo
               exactscat_up(n,1,1) = sum * tms(n)
               if ( Qvary(n) ) then
                  do q = 1, Qnums(n)
                     if ( Lvarymoms(n,q) ) then
                        L_sum = zero
                        do L = 0, nmoments_input
                           L_sum = L_sum + GenSpher(L,1,V)  * L_Greekmat(n,L,1,1,q)
                        enddo
                        L_exactscat_up(n,1,1,q) = L_sum * tms(n) + sum * L_tms(n,q)
                     else
                        L_exactscat_up(n,1,1,q) = sum * L_tms(n,q)
                     endif
                  enddo
               endif
            endif
         enddo
      endif

!  Vector with Sunlight

      if ( nstokes .gt. 1 .and. do_sunlight ) then
         do n = 1, nlayers
            if ( layermask_up(n) ) then
               fmat(1:2) = zero
               do L = 0, nmoments_input
                  fmat(1) = fmat(1) + Genspher(L,1,V) * greekmat(n,L,1,1)
                  fmat(2) = fmat(2) + Genspher(L,2,V) * greekmat(n,L,1,2)
               enddo
               exactscat_up(n,1,1) = + fmat(1)
               exactscat_up(n,2,1) = - fmat(2) * Rotations(3,V)
               exactscat_up(n,3,1) = + fmat(2) * Rotations(4,V)
               if ( Qvary(n) ) then
                  do q = 1, Qnums(n)
                     if ( Lvarymoms(n,q) ) then
                        L_fmat(1:2) = zero
                        do L = 0, nmoments_input
                           L_fmat(1) = L_fmat(1) + GenSpher(L,1,V)  * L_Greekmat(n,L,1,1,q)
                           L_fmat(2) = L_fmat(2) + GenSpher(L,2,V)  * L_Greekmat(n,L,1,2,q)
                        enddo
                        L_exactscat_up(n,1,1,q) = + L_fmat(1)
                        L_exactscat_up(n,2,1,q) = - L_fmat(2) * Rotations(3,V)
                        L_exactscat_up(n,3,1,q) = + L_fmat(2) * Rotations(4,V)
                     endif
                     L_exactscat_up(n,1:ns,1,q) = L_exactscat_up(n,1:ns,1,q) *   tms(n) &
                                                  + exactscat_up(n,1:ns,1)   * L_tms(n,q)
                  enddo
               endif
               exactscat_up(n,1:ns,1) = tms(n) * exactscat_up(n,1:ns,1) 
            endif
         enddo
      endif

!  Vector General case, USE FULL 4X4 MATRIX
!    CODE INTRODUCED BUT NOT TESTED, 05 OCTOBER 2010

      if ( nstokes .gt. 1 .and. .not. do_sunlight ) then
         do n = 1, nlayers
            if ( layermask_up(n) ) then
               fmat = zero
               do L = 0, nmoments_input
                  fmat(1) = fmat(1) + Genspher(L,1,V) * greekmat(n,L,1,1)
                  fmat(2) = fmat(2) + Genspher(L,2,V) * greekmat(n,L,1,2)
                  sum23 = greekmat(n,L,2,2) + greekmat(n,L,3,3)
                  dif23 = greekmat(n,L,2,2) - greekmat(n,L,3,3)
                  fmat(3) = fmat(3) + Genspher(L,3,V) * sum23
                  fmat(4) = fmat(4) + Genspher(L,4,V) * dif23
               enddo
               fmat(3) = ( fmat(3) + fmat(4) ) * 0.5_fpk
               fmat(4) = ( fmat(3) - fmat(4) )
               if ( nstokes.eq.4) then
                  do L = 0, nmoments_input
                     fmat(5) = fmat(5) + Genspher(L,2,V) * greekmat(n,L,3,4)
                     fmat(6) = fmat(6) + Genspher(L,1,V) * greekmat(n,L,4,4)
                  enddo
               endif
               help3c1 = fmat(3) * Rotations(1,V)
               help3s1 = fmat(3) * Rotations(2,V)
               help4c1 = fmat(4) * Rotations(1,V)
               help4s1 = fmat(4) * Rotations(2,V)
               exactscat_up(n,1,1) = + fmat(1)
               exactscat_up(n,2,1) = - fmat(2) * Rotations(3,V)
               exactscat_up(n,3,1) = + fmat(2) * Rotations(4,V)
               exactscat_up(n,2,2) = + help3c1 * Rotations(3,V) - help4s1 * Rotations(4,V)
               exactscat_up(n,2,3) = - help3s1 * Rotations(3,V) - help4c1 * Rotations(4,V)
               exactscat_up(n,3,2) = + help3c1 * Rotations(4,V) + help4s1 * Rotations(3,V)
               exactscat_up(n,3,3) = - help3s1 * Rotations(4,V) + help4c1 * Rotations(3,V)
               if ( nstokes .eq. 4 ) then
                  exactscat_up(n,2,4) = - fmat(5) * Rotations(4,V) 
                  exactscat_up(n,4,2) = - fmat(5) * Rotations(2,V) 
                  exactscat_up(n,3,4) = + fmat(5) * Rotations(3,V) 
                  exactscat_up(n,4,3) = - fmat(5) * Rotations(1,V) 
                  exactscat_up(n,4,4) = + fmat(6)
               endif
               if ( Qvary(n) ) then
                  do q = 1, Qnums(n)
                     if ( Lvarymoms(n,q) ) then
                        L_fmat = zero
                        do L = 0, nmoments_input
                           L_fmat(1) = L_fmat(1) + GenSpher(L,1,V)  * L_Greekmat(n,L,1,1,q)
                           L_fmat(2) = L_fmat(2) + GenSpher(L,2,V)  * L_Greekmat(n,L,1,2,q)
                           sum23 = L_greekmat(n,L,2,2,q) + L_greekmat(n,L,3,3,q)
                           dif23 = L_greekmat(n,L,2,2,q) - L_greekmat(n,L,3,3,q)
                           L_fmat(3) = L_fmat(3) + Genspher(L,3,V) * sum23
                           L_fmat(4) = L_fmat(4) + Genspher(L,4,V) * dif23
                       enddo
                       L_fmat(3) = ( L_fmat(3) + L_fmat(4) ) * 0.5_fpk
                       L_fmat(4) = ( L_fmat(3) - L_fmat(4) )
                       if ( nstokes.eq.4) then
                          do L = 0, nmoments_input
                             L_fmat(5) = L_fmat(5) + Genspher(L,2,V) * L_greekmat(n,L,3,4,q)
                             L_fmat(6) = L_fmat(6) + Genspher(L,1,V) * L_greekmat(n,L,4,4,q)
                          enddo
                       endif
                       help3c1 = L_fmat(3) * Rotations(1,V)
                       help3s1 = L_fmat(3) * Rotations(2,V)
                       help4c1 = L_fmat(4) * Rotations(1,V)
                       help4s1 = L_fmat(4) * Rotations(2,V)
                       L_exactscat_up(n,1,1,q) = + L_fmat(1)
                       L_exactscat_up(n,2,1,q) = - L_fmat(2) * Rotations(3,V)
                       L_exactscat_up(n,3,1,q) = + L_fmat(2) * Rotations(4,V)
                       L_exactscat_up(n,2,2,q) = + help3c1 * Rotations(3,V) - help4s1 * Rotations(4,V)
                       L_exactscat_up(n,2,3,q) = - help3s1 * Rotations(3,V) - help4c1 * Rotations(4,V)
                       L_exactscat_up(n,3,2,q) = + help3c1 * Rotations(4,V) + help4s1 * Rotations(3,V)
                       L_exactscat_up(n,3,3,q) = - help3s1 * Rotations(4,V) + help4c1 * Rotations(3,V)
                       if ( nstokes .eq. 4 ) then
                          L_exactscat_up(n,2,4,q) = - L_fmat(5) * Rotations(4,V) 
                          L_exactscat_up(n,4,2,q) = - L_fmat(5) * Rotations(2,V) 
                          L_exactscat_up(n,3,4,q) = + L_fmat(5) * Rotations(3,V) 
                          L_exactscat_up(n,4,3,q) = - L_fmat(5) * Rotations(1,V) 
                          L_exactscat_up(n,4,4,q) = + L_fmat(6)
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

!  Initialize, only to layer Ncrit if applicable

      Attenuations   = zero  ; Attenuations_fine   = zero
      LC_Attenuations = zero ; LC_Attenuations_fine = zero
      Suntau         = zero  ; LC_suntau            = zero
      nstart = nlayers ; if (Ncrit(v).ne.0) nstart = nCrit(v)

!  Attenuations to End points (including TOA). All representations
!    MUST go all the way to NLAYERS (surface term required)

      do n = 0, nlayers
         sumd = ZERO
         do k = 1, ntraverse(n,v)
            sumd = sumd + extinction(k) * sunpaths(n,k,v)
         enddo
         suntau(n) = sumd
         If (sumd .lt. cutoff ) Attenuations(n) = exp( - sumd )
         if ( do_columnwfs ) then
            do q = 1, n_columnwfs
               L_sumd = ZERO
               do k = 1, ntraverse(n,v)
                  L_sumd = L_sumd + L_extinction(k,q) * sunpaths(n,k,v)
               enddo
               LC_suntau(n,q) = L_sumd
               LC_Attenuations(n,q) = - Attenuations(n) * L_sumd
            enddo
         endif
      enddo

!  Adjust nstart

      do n = 1, nlayers
         if ( layermask_up(n) .and. attenuations(n-1).ne.zero )  nstart = n
      enddo

!  Enhanced-spherical, fine-layer attenuations

      if ( do_enhanced_ps ) then
         do n = 1, nstart
            if ( layermask_up(n) ) then
               do j = 1, nfinedivs(n,v)
                  sumd = zero
                  do k = 1, ntraverse_fine(n,j,v)
                     sumd = sumd + extinction(k) * sunpaths_fine(n,k,j,v)
                  enddo
                  if (sumd .lt. cutoff ) Attenuations_fine(n,j) = exp( - sumd )
                  if ( do_columnwfs ) then
                     do q = 1, n_columnwfs
                        L_sumd = zero
                        do k = 1, ntraverse_fine(n,j,v)
                           L_sumd = L_sumd + L_extinction(k,q) * sunpaths_fine(n,k,j,v)
                        enddo
                        LC_Attenuations_fine(n,j,q) = - Attenuations_fine(n,j) * L_sumd 
                     enddo
                  endif
               enddo
            endif
         enddo
      endif

!  Layer integrated Solar sources
!  ==============================

!  Plane-Parallel or Regular-PS (Average secant formulation)

      if ( do_RegPSorPP ) then
         do n = nlayers, 1, -1

!  LOS transmittance (not for the Horizontal case)

            if ( Mu1(v) .gt. zero ) then
               lostau = deltaus(n)  / Mu1(v)
               if ( lostau .lt. cutoff ) lostrans_up(n) = exp( - lostau )
               if ( do_columnwfs ) then
                  if ( Qvary(n) ) then
                     do q = 1, Qnums(n)
                        L_lostau(q)        = L_deltaus(n,q) / Mu1(v)
                        L_lostrans_up(n,q) = - L_lostau(q) * lostrans_up(n)
                     enddo
                  endif
               endif
            endif

!  Multipliers

            if ( layermask_up(n) .and. n.le.nstart  ) then
              if ( Mu1(v) .gt. zero ) then
                if ( attenuations(n-1).ne.zero ) then
                  factor1 = Attenuations(n-1) - Attenuations(n)*lostrans_up(n)
                  factor2 = (suntau(n) - suntau(n-1))/lostau
                  factor3 = one + factor2
                  multiplier_up(n) = factor1 / factor3
                  if ( do_columnwfs ) then
                     do q = 1, n_columnwfs
                        L_factor1 = LC_Attenuations(n-1,q) - LC_Attenuations(n,q)*lostrans_up(n)
                        L_factor2 = ( LC_suntau(n,q) - LC_suntau(n-1,q) ) / lostau
                        L_factor1 = L_factor1 - Attenuations(n) * L_lostrans_up(n,q)
                        L_factor2 = L_factor2 - factor2 *  L_lostau(q) / lostau
                        L_multiplier_up(n,q) = ( L_factor1 - multiplier_up(n) * L_factor2 ) / factor3
                     enddo
                  endif
                endif
              else
                if ( attenuations(n-1).ne.zero ) then
                  multiplier_up(n) = Attenuations(n-1)
                  if ( do_columnwfs ) then
                     do q = 1, n_columnwfs
                        L_multiplier_up(n,q) = LC_Attenuations(n-1,q)
                     enddo
                  endif
                endif
              endif
            endif

!  End layers and regular-PS formulation

         enddo
      endif

!  Enhanced PS multipliers and LOSTRANS: special case (nadir viewing)
!  ------------------------------------------------------------------

      if ( do_enhanced_ps .and. doNadir(v) ) then
         do n = nlayers, 1, -1

!  LOS transmittance

            kn     = extinction(n)
            lostau = deltaus(n)
            if ( lostau .lt. cutoff ) lostrans_up(n) = exp( - lostau )
            if ( do_columnwfs ) then
               if ( Qvary(n) ) then
                  do q = 1, Qnums(n)
                     L_lostau(q)        = L_deltaus(n,q)
                     L_lostrans_up(n,q) = - L_lostau(q) * lostrans_up(n)
                  enddo
               endif
            endif

!  Multipliers

            if ( layermask_up(n) .and. n.le.nstart  ) then
               sum = zero
               do j = 1, nfinedivs(n,v)
                  argum(j) = xfine(n,j,v)
                  tran(j)  = exp ( - argum(j) * kn )
                  func(j)  = attenuations_fine(n,j) * tran(j)
                  sum = sum + func(j) * wfine(n,j,v)
               enddo
               multiplier_up(n) = sum * kn
               if ( do_columnwfs ) then
                  do q = 1, n_columnwfs
                     L_sum = zero
                     do j = 1, nfinedivs(n,v)
                        L_tran = - argum(j) * L_extinction(n,q)
                        L_func = LC_attenuations_fine(n,j,q) * tran(j) + L_tran * func(j)
                        L_sum  = L_sum + L_func * wfine(n,j,v)
                     enddo
                     L_multiplier_up(n,q)  = L_sum * kn + L_extinction(n,q) * sum
                  enddo
               endif
            endif

!  End layer loop and Nadir Enhanced PS case

         enddo
      endif

!  Enhanced PS multipliers and LOSTRANS: General case
!  --------------------------------------------------

      if ( do_enhanced_ps .and. .not. doNadir(v) ) then
         do n = nlayers, 1, -1

!  LOS transmittance

            cot_2 = cota(n-1,v) ; cot_1 = cota(n,v)
            kn = extinction(n) ;  ke = raycon(v) * kn ; cons = raycon(v) * ( cot_2 - cot_1 )
            tran_1 = kn * cons
            if ( tran_1 .lt. cutoff ) lostrans_up(n) = exp ( - tran_1 )
            if ( do_columnwfs ) then
               if ( Qvary(n) ) then
                  do q = 1, Qnums(n)
                     L_lostau(q)        = L_extinction(n,q) * cons
                     L_lostrans_up(n,q) = - L_lostau(q) * lostrans_up(n)
                  enddo
               endif
            endif

!  Multipliers

            if ( layermask_up(n) .and. n.le.nstart  ) then
               sum = zero
               do j = 1, nfinedivs(n,v)
                  argum(j) = Raycon(v) * ( cot_2 - cotfine(n,j,v) )
                  tran(j)  = exp ( - kn * argum(j) )
                  func(j)  = attenuations_fine(n,j) * csqfine(n,j,v) * tran(j)
                  sum      = sum + func(j) * wfine(n,j,v)
               enddo
               Multiplier_up(n) = sum * ke 
               if ( do_columnwfs ) then
                  do q = 1, n_columnwfs
                     L_sum = zero
                     do j = 1, nfinedivs(n,v)
                        L_tran = - argum(j) * L_extinction(n,q)
                        L_func = LC_attenuations_fine(n,j,q) * csqfine(n,j,v) * tran(j) + L_tran * func(j)
                        L_sum  = L_sum + L_func * wfine(n,j,v)
                     enddo
                     L_multiplier_up(n,q)  = L_sum * ke + L_extinction(N,q) * Raycon(v) * sum
                  enddo
               endif
            endif

!  End layer loop and general Enhanced PS case

         enddo
      endif

!  Layer sources
!  -------------

      do n = nlayers, 1, -1
         if ( layermask_up(n) .and. n.le.nstart  ) then
            do o1 = 1, nstokes
               shelp(o1) = dot_product(exactscat_up(n,o1,1:ns),fluxvec(1:ns))
               sources_up(n,o1) = shelp(o1) * multiplier_up(n)
            enddo
            if ( do_columnwfs ) then
               do q = 1, n_columnwfs
                  do o1 = 1, nstokes
                     LC_sources_up(n,o1,q) = shelp(o1) * L_multiplier_up(n,q)
                     L_Shelp = dot_product(L_exactscat_up(n,o1,1:ns,q),fluxvec(1:ns))
                     LC_sources_up(n,o1,q) =  LC_sources_up(n,o1,q) + L_Shelp * multiplier_up(n)
                  enddo
               enddo
            endif
         endif
      enddo

!  Source function integration
!  ===========================

!  NLEVEL = Layer index for given optical depth
!  Cumulative source terms : Loop over layers working upwards from NSTART to level NUT,
!  Check for updating the recursion

!  INTENSITY Main loop over all output optical depths
!          Cumulative source term will be saved

      NC = 0 
      CUMSOURCE_UP(NC,:) = zero
      CUMSOURCE_DB(NC,:) = zero
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  Surface term

      RHELP = zero; M4 = 4.0_fpk * MU0(v) ; M4A = M4 * attenuations(nlayers)
      if ( DO_LAMBERTIAN ) then
         RHELP(1) = M4 * REFLEC(1,1,V) * Fluxvec(1)
         CUMSOURCE_DB(NC,1) = RHELP(1) * attenuations(nlayers)
      else
         do o1 = 1, nstokes
            RHELP(O1) = M4 * dot_product(REFLEC(O1,1:ns,V),fluxvec(1:ns))
            CUMSOURCE_DB(NC,o1) = RHELP(O1) * attenuations(nlayers)
         enddo
      endif

!  Main loop over all output optical depths
!     NLEVEL = Layer index for given optical depth
!     Cumulative source terms : Loop over layers working upwards from NSTART to level NUT,
!     Check for updating the recursion

      DO UTA = N_USER_LEVELS, 1, -1
         NUT    = USER_LEVELS(UTA) + 1
         DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            do o1 = 1, nstokes
               CUMSOURCE_DB(NC,O1) = LOSTRANS_UP(N) * CUMSOURCE_DB(NC-1,O1)
               CUMSOURCE_UP(NC,O1) = LOSTRANS_UP(N) * CUMSOURCE_UP(NC-1,O1) + SOURCES_UP(N,O1)
            enddo
         ENDDO
         do o1 = 1, nstokes
            STOKES_UP(UTA,O1,V) = FLUX * CUMSOURCE_UP(NC,O1)
            STOKES_DB(UTA,O1,V) = FLUX * CUMSOURCE_DB(NC,O1)
         enddo
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1 ; NUT_PREV = NUT
      ENDDO

!  Surface WFs

      if ( do_surfacewfs ) then
         LS_CUMSOURCE = zero
         do q = 1, n_surfacewfs
            if ( DO_LAMBERTIAN ) then
               LS_cumsource(1,q) = M4A * LS_REFLEC(1,1,v,q)
            else
               do o1 = 1, nstokes
                  LS_cumsource(o1,q) = M4A * dot_product(LS_REFLEC(O1,1:ns,v,q),fluxvec(1:ns))
               enddo
            endif
         enddo
         NSTART = NLAYERS ; NUT_PREV = NSTART + 1
         DO UTA = N_USER_LEVELS, 1, -1
            NUT    = USER_LEVELS(UTA) + 1
            DO N = NSTART, NUT, -1
               do q = 1, n_surfacewfs
                  LS_cumsource(1:ns,q) = LOSTRANS_UP(N) * LS_CUMSOURCE(1:ns,q)
               enddo
            ENDDO
            do q = 1, n_surfacewfs
               LS_JACOBIANS_DB(UTA,1:ns,V,Q) = FLUX * LS_CUMSOURCE(1:ns,Q)
            enddo
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1 ; NUT_PREV = NUT
         ENDDO
      endif

!  Column Wfs (Atmospheric term)

      if ( do_columnwfs ) then
         L_CUMSOURCE = zero
         NSTART = NLAYERS ; NUT_PREV = NSTART + 1
         DO UTA = N_USER_LEVELS, 1, -1
            NUT    = USER_LEVELS(UTA) + 1
            DO N = NSTART, NUT, -1
               NC = NLAYERS + 1 - N
               do q = 1, n_columnwfs
                  L_cumsource(1:ns,q) = LC_SOURCES_UP(N,1:ns,Q)         + &
                                L_LOSTRANS_UP(N,Q) * CUMSOURCE_UP(NC-1,1:ns) + &
                                  LOSTRANS_UP(N)   * L_CUMSOURCE(1:ns,Q)
               enddo
            ENDDO
            do q = 1, n_columnwfs
               LC_JACOBIANS_UP(UTA,1:ns,V,Q) = FLUX * L_CUMSOURCE(1:ns,Q)
            enddo
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
            NUT_PREV = NUT
         ENDDO
      endif

!  Column Wfs (Surface term)

      if ( do_columnwfs ) then
         do q = 1, n_columnwfs
            L_CUMSOURCE(1:ns,q) = RHELP(1:ns) * LC_attenuations(nlayers,q)
         enddo
         NSTART = NLAYERS  ;  NUT_PREV = NSTART + 1
         DO UTA = N_USER_LEVELS, 1, -1
            NUT    = USER_LEVELS(UTA) + 1
            DO N = NSTART, NUT, -1
               NC = NLAYERS + 1 - N
               do q = 1, n_columnwfs
                  L_cumsource(1:ns,q) =  L_LOSTRANS_UP(N,Q) *   CUMSOURCE_DB(NC-1,1:ns) + &
                                           LOSTRANS_UP(N)   * L_CUMSOURCE(1:ns,Q)
               enddo
            ENDDO
            do q = 1, n_columnwfs
               LC_JACOBIANS_DB(UTA,1:ns,V,Q) = FLUX * L_CUMSOURCE(1:ns,Q)
            enddo
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1 ; NUT_PREV = NUT
         ENDDO
      endif

!  End GEometry Loop

   enddo

!  Finish

   return
end subroutine SSV_Integral_ILCS_UP

!

subroutine SSV_Integral_ILCS_DN &
   ( maxgeoms, maxlayers, maxfinelayers, maxmoments_input,                      & ! Inputs (dimensioning)
     max_user_levels, max_atmoswfs, do_sunlight, do_deltam_scaling,             & ! Inputs (Flags - General)
     do_PlanPar, do_regular_ps, do_enhanced_ps, doNadir, nstokes,               & ! Inputs (Flags - General)
     do_columnwfs, n_columnwfs, Lvarymoms,                                      & ! Inputs (control, Jacobian )
     ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,    & ! Inputs (control,  output)
     extinction, deltaus, omega, truncfac, greekmat, flux, fluxvec,             & ! Inputs (Optical)
     L_extinction, L_deltaus, L_omega, L_truncfac, L_greekmat,                  & ! Inputs (Optical - Linearized)
     Mu1, GenSpher, Rotations, NCrit, RadCrit, CotCrit,                         & ! Inputs (Geometry)
     xfine, wfine, csqfine, cotfine, Raycon, radii, cota,                       & ! Inputs (Geometry)
     sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,                        & ! Inputs (Geometry)
     Stokes_dn, LC_Jacobians_dn )                                                 ! Output

!  Stand-alone routine for Downwelling Solar-beam Single-scatter (SS)
!    computation of stokes vectors and LCS Jacobians. Inputs: geometry, spherical functions, optical properties.

!  This version, revised by R. Spurr, 01 June 2012
!   Extension to multiple geometries, 19 December 2012

   use FO_Taylor_m

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

   LOGICAL, Intent(in) :: DO_SUNLIGHT
   LOGICAL, Intent(in) :: DO_DELTAM_SCALING

   logical, Intent(in) :: DO_PLANPAR
   logical, Intent(in) :: DO_REGULAR_PS
   logical, Intent(in) :: DO_ENHANCED_PS
   logical, Intent(in) :: DONADIR(MAXGEOMS)

!  Jacobian Flag

   LOGICAL, Intent(in) :: do_columnwfs

!  Layer and Level Control Numbers, Number of Moments

   integer, Intent(in) :: NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   integer, Intent(in) :: NGEOMS, NMOMENTS_INPUT, NSTOKES

   integer, Intent(in) :: N_USER_LEVELS
   integer, Intent(in) :: USER_LEVELS ( MAX_USER_LEVELS )

!  Jacobian control

   INTEGER, Intent(in) :: n_columnwfs
   LOGICAL, Intent(in) :: Lvarymoms (maxlayers,max_atmoswfs)

!  optical inputs
!  --------------

!  Atmosphere

   real(fpk), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(fpk), Intent(in) :: OMEGA       ( MAXLAYERS )
   real(fpk), Intent(in) :: TRUNCFAC    ( MAXLAYERS )
   real(fpk), Intent(in) :: GREEKMAT    ( MAXLAYERS,0:MAXMOMENTS_INPUT, 4, 4 )

!  Solar Flux

   real(fpk), Intent(in) :: FLUX, FLUXVEC(4)

!  Linearized optical inputs

   real(fpk), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
   real(fpk), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
   real(fpk), Intent(in) :: L_OMEGA       ( MAXLAYERS, max_atmoswfs )
   real(fpk), Intent(in) :: L_TRUNCFAC    ( MAXLAYERS, max_atmoswfs )
   real(fpk), Intent(in) :: L_GREEKMAT    ( MAXLAYERS,0:MAXMOMENTS_INPUT, 4, 4, max_atmoswfs )

!  Geometrical inputs
!  ------------------

!  Ray constant, Cotangents, Critical layer
!    Mu1 = cos(alpha_boa), required for the Regular PS only

   INTEGER  , Intent(in)  :: NCrit(maxgeoms)
   REAL(fpk), Intent(in)  :: RadCrit(maxgeoms), CotCrit(maxgeoms)
   REAL(fpk), Intent(in)  :: Raycon(maxgeoms), cota(0:maxlayers,maxgeoms), radii(0:maxlayers)
   REAL(fpk), Intent(in)  :: Mu1(maxgeoms)

!  solar paths 

   INTEGER  , Intent(in)  :: ntraverse  (0:maxlayers,maxgeoms)
   REAL(fpk), Intent(in)  :: sunpaths   (0:maxlayers,maxlayers,maxgeoms)
   INTEGER  , Intent(in)  :: ntraverse_fine(maxlayers,maxfinelayers,maxgeoms)
   REAL(fpk), Intent(in)  :: sunpaths_fine (maxlayers,maxlayers,maxfinelayers,maxgeoms)

!  Generalized spherical functions.
!    Rotations(1-4)    = C1, S1, C2, S2

   REAL(fpk), Intent(in)  :: GenSpher(0:maxmoments_input,4,maxgeoms)
   REAL(fpk), Intent(in)  :: Rotations(4,maxgeoms)

!  LOS Quadratures for Enhanced PS

   REAL(fpk), Intent(in)  :: xfine   (maxlayers,maxfinelayers,maxgeoms)
   REAL(fpk), Intent(in)  :: wfine   (maxlayers,maxfinelayers,maxgeoms)
   REAL(fpk), Intent(in)  :: csqfine (maxlayers,maxfinelayers,maxgeoms)
   REAL(fpk), Intent(in)  :: cotfine (maxlayers,maxfinelayers,maxgeoms)

!  outputs
!  -------

   real(fpk), Intent(Out)  :: Stokes_dn        ( max_user_levels, 4, maxgeoms )
   real(fpk), Intent(Out)  :: LC_Jacobians_dn  ( max_user_levels, 4, maxgeoms, max_atmoswfs )

!  LOCAL
!  -----

!  Attenuations

   real(fpk)  :: suntau              (0:maxlayers)
   real(fpk)  :: attenuations        (0:maxlayers)
   real(fpk)  :: attenuations_fine   (maxlayers,maxfinelayers)

   real(fpk)  :: LC_suntau            (0:maxlayers,max_atmoswfs)
   real(fpk)  :: LC_attenuations      (0:maxlayers,max_atmoswfs)
   real(fpk)  :: LC_attenuations_fine (maxlayers,maxfinelayers,max_atmoswfs)

!  Scattering

   real(fpk)  :: tms            (maxlayers)
   real(fpk)  :: exactscat_dn   (maxlayers,4,4)

   real(fpk)  :: L_tms          (maxlayers,max_atmoswfs)
   real(fpk)  :: L_exactscat_dn (maxlayers,4,4, max_atmoswfs)

!  Source function integration results

   real(fpk)  :: sources_dn       ( maxlayers, 4 )
   real(fpk)  :: lostrans_dn      ( maxlayers )
   real(fpk)  :: multiplier_dn    ( maxlayers )

   real(fpk)  :: LC_sources_dn    ( maxlayers,4,max_atmoswfs )
   real(fpk)  :: L_lostrans_dn    ( maxlayers,  max_atmoswfs )
   real(fpk)  :: L_multiplier_dn  ( maxlayers,  max_atmoswfs )

!  Local cumulative source terms

   real(fpk)  :: cumsource_dn      ( 0:maxlayers, 4 )
   real(fpk)  :: L_cumsource       ( 4, max_atmoswfs )

!  Regular_PS or plane-parallel flag

   logical    :: do_RegPSorPP

!  Help

   integer    :: n, ns, k, j, q, L, v, o1, uta, nstart, nc, nut, nut_prev, Qnums(maxlayers)
   logical    :: layermask_dn(maxlayers), Qvary(maxlayers)

   real(fpk)  :: argum(maxfinelayers), tran(maxfinelayers), func(maxfinelayers)
   real(fpk)  :: cons, help, sum, tran_1, kn, ke, factor1, factor2, factor3, shelp(4)
   real(fpk)  :: cot_1, cot_2, L_help, L_sum, L_tran, L_func, L_factor1, L_factor2, sumd, L_sumd
   real(fpk)  :: lostau, L_lostau(max_atmoswfs), fmat(6), L_fmat(6), sum23, dif23, L_Shelp
   real(fpk)  :: help3c1, help3s1, help4c1, help4s1, cot_c, rdiff

   real(fpk)  :: eps, lospath, term2, mult, L_term2, L_mult

   real(fpk), parameter  :: cutoff = 88.0_fpk
   real(fpk), parameter  :: zero   = 0.0_fpk
   real(fpk), parameter  :: one    = 1.0_fpk

!  Zero the output 

   STOKES_DN       = zero
   LC_JACOBIANS_DN = zero

!  Regular_PS or plane-parallel flag

   do_RegPSorPP = (do_regular_ps .or. do_PlanPar)

!  Bookkeeping

   ns = nstokes ; fmat = zero
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

      lostrans_dn   = zero  ; sources_dn    = zero ; exactscat_dn   = zero ; multiplier_dn   = zero
      L_lostrans_dn = zero  ; LC_sources_dn = zero ; L_exactscat_dn = zero ; L_multiplier_dn = zero
 
!  Scattering functions and Linearization
!  --------------------------------------

!  Scalar only

      if ( nstokes .eq. 1 ) then
         do n = 1, nlayers
            if ( layermask_dn(n) ) then
              sum = zero
               do L = 0, nmoments_input
                  sum = sum + GenSpher(L,1,V) * Greekmat(n,L,1,1)
               enddo
               exactscat_dn(n,1,1) = sum * tms(n)
               if ( Qvary(n) ) then
                  do q = 1, Qnums(n)
                     if ( Lvarymoms(n,q) ) then
                        L_sum = zero
                        do L = 0, nmoments_input
                           L_sum = L_sum + GenSpher(L,1,V)  * L_Greekmat(n,L,1,1,q)
                        enddo
                        L_exactscat_dn(n,1,1,q) = L_sum * tms(n) + sum * L_tms(n,q)
                     else
                        L_exactscat_dn(n,1,1,q) = sum * L_tms(n,q)
                     endif
                  enddo
               endif               
            endif
         enddo
      endif

!  Vector with Sunlight

      if ( nstokes .gt. 1 .and. do_sunlight ) then
         do n = 1, nlayers
            if ( layermask_dn(n) ) then
               fmat(1:2) = zero
               do L = 0, nmoments_input
                  fmat(1) = fmat(1) + Genspher(L,1,V) * greekmat(n,L,1,1)
                  fmat(2) = fmat(2) + Genspher(L,2,V) * greekmat(n,L,1,2)
               enddo
               exactscat_dn(n,1,1) = + fmat(1)
               exactscat_dn(n,2,1) = - fmat(2) * Rotations(3,V)
               exactscat_dn(n,3,1) = + fmat(2) * Rotations(4,V)
               if ( Qvary(n) ) then
                  do q = 1, Qnums(n)
                     if ( Lvarymoms(n,q) ) then
                        L_fmat(1:2) = zero
                        do L = 0, nmoments_input
                           L_fmat(1) = L_fmat(1) + GenSpher(L,1,V)  * L_Greekmat(n,L,1,1,q)
                           L_fmat(2) = L_fmat(2) + GenSpher(L,2,V)  * L_Greekmat(n,L,1,2,q)
                        enddo
                        L_exactscat_dn(n,1,1,q) = + L_fmat(1)
                        L_exactscat_dn(n,2,1,q) = - L_fmat(2) * Rotations(3,V)
                        L_exactscat_dn(n,3,1,q) = + L_fmat(2) * Rotations(4,V)
                     endif
                     L_exactscat_dn(n,1:ns,1,q) = L_exactscat_dn(n,1:ns,1,q) *   tms(n) &
                                                  + exactscat_dn(n,1:ns,1)   * L_tms(n,q)
                  enddo
               endif               
               exactscat_dn(n,1:ns,1) = tms(n) * exactscat_dn(n,1:ns,1) 
            endif
         enddo
      endif

!  Vector General case, USE FULL 4X4 MATRIX
!    CODE INTRODUCED BUT NOT TESTED, 05 OCTOBER 2010

      if ( nstokes .gt. 1 .and. .not. do_sunlight ) then
         do n = 1, nlayers
            if ( layermask_dn(n) ) then
               fmat = zero
               do L = 0, nmoments_input
                  fmat(1) = fmat(1) + Genspher(L,1,V) * greekmat(n,L,1,1)
                  fmat(2) = fmat(2) + Genspher(L,2,V) * greekmat(n,L,1,2)
                  sum23 = greekmat(n,L,2,2) + greekmat(n,L,3,3)
                  dif23 = greekmat(n,L,2,2) - greekmat(n,L,3,3)
                  fmat(3) = fmat(3) + Genspher(L,3,V) * sum23
                  fmat(4) = fmat(4) + Genspher(L,4,V) * dif23
               enddo
               fmat(3) = ( fmat(3) + fmat(4) ) * 0.5_fpk
               fmat(4) = ( fmat(3) - fmat(4) )
               if ( nstokes.eq.4) then
                  do L = 0, nmoments_input
                     fmat(5) = fmat(5) + Genspher(L,2,V) * greekmat(n,L,3,4)
                     fmat(6) = fmat(6) + Genspher(L,1,V) * greekmat(n,L,4,4)
                  enddo
               endif
               help3c1 = fmat(3) * Rotations(1,V)
               help3s1 = fmat(3) * Rotations(2,V)
               help4c1 = fmat(4) * Rotations(1,V)
               help4s1 = fmat(4) * Rotations(2,V)
               exactscat_dn(n,1,1) = + fmat(1)
               exactscat_dn(n,2,1) = - fmat(2) * Rotations(3,V)
               exactscat_dn(n,3,1) = + fmat(2) * Rotations(4,V)
               exactscat_dn(n,2,2) = + help3c1 * Rotations(3,V) - help4s1 * Rotations(4,V)
               exactscat_dn(n,2,3) = - help3s1 * Rotations(3,V) - help4c1 * Rotations(4,V)
               exactscat_dn(n,3,2) = + help3c1 * Rotations(4,V) + help4s1 * Rotations(3,V)
               exactscat_dn(n,3,3) = - help3s1 * Rotations(4,V) + help4c1 * Rotations(3,V)
               if ( nstokes .eq. 4 ) then
                  exactscat_dn(n,2,4) = - fmat(5) * Rotations(4,V) 
                  exactscat_dn(n,4,2) = - fmat(5) * Rotations(2,V) 
                  exactscat_dn(n,3,4) = + fmat(5) * Rotations(3,V) 
                  exactscat_dn(n,4,3) = - fmat(5) * Rotations(1,V) 
                  exactscat_dn(n,4,4) = + fmat(6)
               endif
               if ( Qvary(n) ) then
                  do q = 1, Qnums(n)
                     if ( Lvarymoms(n,q) ) then
                        L_fmat = zero
                        do L = 0, nmoments_input
                           L_fmat(1) = L_fmat(1) + GenSpher(L,1,V)  * L_Greekmat(n,L,1,1,q)
                           L_fmat(2) = L_fmat(2) + GenSpher(L,2,V)  * L_Greekmat(n,L,1,2,q)
                           sum23 = L_greekmat(n,L,2,2,q) + L_greekmat(n,L,3,3,q)
                           dif23 = L_greekmat(n,L,2,2,q) - L_greekmat(n,L,3,3,q)
                           L_fmat(3) = L_fmat(3) + Genspher(L,3,V) * sum23
                           L_fmat(4) = L_fmat(4) + Genspher(L,4,V) * dif23
                       enddo
                       L_fmat(3) = ( L_fmat(3) + L_fmat(4) ) * 0.5_fpk
                       L_fmat(4) = ( L_fmat(3) - L_fmat(4) )
                       if ( nstokes.eq.4) then
                          do L = 0, nmoments_input
                             L_fmat(5) = L_fmat(5) + GenSpher(L,2,V) * L_greekmat(n,L,3,4,q)
                             L_fmat(6) = L_fmat(6) + GenSpher(L,1,V) * L_greekmat(n,L,4,4,q)
                          enddo
                       endif
                       help3c1 = L_fmat(3) * Rotations(1,V)
                       help3s1 = L_fmat(3) * Rotations(2,V)
                       help4c1 = L_fmat(4) * Rotations(1,V)
                       help4s1 = L_fmat(4) * Rotations(2,V)
                       L_exactscat_dn(n,1,1,q) = + L_fmat(1)
                       L_exactscat_dn(n,2,1,q) = - L_fmat(2) * Rotations(3,V)
                       L_exactscat_dn(n,3,1,q) = + L_fmat(2) * Rotations(4,V)
                       L_exactscat_dn(n,2,2,q) = + help3c1 * Rotations(3,V) - help4s1 * Rotations(4,V)
                       L_exactscat_dn(n,2,3,q) = - help3s1 * Rotations(3,V) - help4c1 * Rotations(4,V)
                       L_exactscat_dn(n,3,2,q) = + help3c1 * Rotations(4,V) + help4s1 * Rotations(3,V)
                       L_exactscat_dn(n,3,3,q) = - help3s1 * Rotations(4,V) + help4c1 * Rotations(3,V)
                       if ( nstokes .eq. 4 ) then
                          L_exactscat_dn(n,2,4,q) = - L_fmat(5) * Rotations(4,V) 
                          L_exactscat_dn(n,4,2,q) = - L_fmat(5) * Rotations(2,V) 
                          L_exactscat_dn(n,3,4,q) = + L_fmat(5) * Rotations(3,V) 
                          L_exactscat_dn(n,4,3,q) = - L_fmat(5) * Rotations(1,V) 
                          L_exactscat_dn(n,4,4,q) = + L_fmat(6)
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

!  Initialize, only to layer Ncrit if applicable

      Attenuations   = zero  ; Attenuations_fine    = zero
      LC_Attenuations = zero ; LC_Attenuations_fine = zero
      Suntau         = zero  ; LC_suntau            = zero
      nstart = nlayers ; if (Ncrit(v).ne.0) nstart = nCrit(v)

!  Attenuations to End points (including TOA). All representations
!    MUST go all the way to NLAYERS (surface term required)

      do n = 0, nlayers
         sumd = ZERO
         do k = 1, ntraverse(n,v)
            sumd = sumd + extinction(k) * sunpaths(n,k,v)
         enddo
         suntau(n) = sumd
         If (sumd .lt. cutoff ) Attenuations(n) = exp( - sumd )
         if ( do_columnwfs ) then
            do q = 1, n_columnwfs
               L_sumd = ZERO
               do k = 1, ntraverse(n,v)
                  L_sumd = L_sumd + L_extinction(k,q) * sunpaths(n,k,v)
               enddo
               LC_suntau(n,q) = L_sumd
               LC_Attenuations(n,q) = - Attenuations(n) * L_sumd
            enddo
         endif
      enddo

!  Adjust nstart

      do n = 1, nlayers
         if ( layermask_dn(n) .and. attenuations(n-1).ne.zero )  nstart = n
      enddo

!  Enhanced-spherical, fine-layer attenuations

      if ( do_enhanced_ps ) then
         do n = 1, nstart
            if ( layermask_dn(n) ) then
               do j = 1, nfinedivs(n,v)
                  sumd = zero
                  do k = 1, ntraverse_fine(n,j,v)
                     sumd = sumd + extinction(k) * sunpaths_fine(n,k,j,v)
                  enddo
                  if (sumd .lt. cutoff ) Attenuations_fine(n,j) = exp( - sumd )
                  if ( do_columnwfs ) then
                     do q = 1, n_columnwfs
                        L_sumd = zero
                        do k = 1, ntraverse_fine(n,j,v)
                           L_sumd = L_sumd + L_extinction(k,q) * sunpaths_fine(n,k,j,v)
                        enddo
                        LC_Attenuations_fine(n,j,q) = - Attenuations_fine(n,j) * L_sumd 
                     enddo
                  endif
               enddo
            endif
         enddo
      endif

!  Layer integrated Solar sources
!  ==============================

!  Plane-parallel or Regular-PS multipliers and LOSTRANS

      if ( do_RegPSorPP ) then
         do n = nlayers, 1, -1

!  LOS transmittance (Not for the Horizontal View)

            if ( Mu1(v) .gt. zero ) then
               !Case: UZA > 0.0
               lostau = deltaus(n) / Mu1(v)
               if ( lostau .lt. cutoff ) lostrans_dn(n) = exp( - lostau )
               if ( do_columnwfs ) then
                  if ( Qvary(n) ) then
                     do q = 1, Qnums(n)
                        L_lostau(q)        = L_deltaus(n,q)  / Mu1(v)
                        L_lostrans_dn(n,q) = - L_lostau(q) * lostrans_dn(n)
                     enddo
                  endif
               endif
            endif

!  Multipliers

            if ( layermask_dn(n) .and. n.le.nstart  ) then
              if ( Mu1(v) .gt. zero ) then
                !Case: UZA > 0.0
                if ( attenuations(n-1).ne.zero ) then
!mick fix 8/11/2013 - added "taylor small" section
                  lospath = lostau/extinction(n)
                  eps     = sunpaths(n,n,v) - lospath
                  if ( abs(eps) .lt. TAYLOR_SMALL ) then
                    term2 = exp( - extinction(n) * sunpaths(n,n,v) )
                    call taylor_series_1 (eps, extinction(n), term2, mult)
                    multiplier_dn(n) = Attenuations(n-1) * lospath * mult
                  else
                    factor1 = Attenuations(n-1)*lostrans_dn(n) - Attenuations(n)
                    factor2 = (suntau(n) - suntau(n-1))/lostau
                    factor3 = factor2 - one
                    multiplier_dn(n) = factor1 / factor3
                  end if

                  if ( do_columnwfs ) then
                     do q = 1, n_columnwfs
                       if ( abs(eps) .lt. TAYLOR_SMALL ) then
!mick test 8/12/2013 - added "taylor small" section
                         L_term2 = -L_extinction(n,q)*sunpaths(n,n,v)*term2
                         call taylor_series_L_1 (eps, extinction(n), term2, &
                                                 L_extinction(n,q), L_term2, L_mult)
                         L_multiplier_dn(n,q) = lospath * &
                           (LC_Attenuations(n-1,q) * mult + Attenuations(n-1) * L_mult)
                       else
                         L_factor1 = LC_Attenuations(n-1,q)*lostrans_dn(n) - LC_Attenuations(n,q)
                         L_factor2 = ( LC_suntau(n,q) - LC_suntau(n-1,q) ) / lostau
                         L_factor1 = L_factor1 + Attenuations(n-1) * L_lostrans_dn(n,q)
                         L_factor2 = L_factor2 - factor2 *  L_lostau(q) / lostau
                         L_multiplier_dn(n,q) = ( L_factor1 - multiplier_dn(n) * L_factor2 ) / factor3
                       end if
                     enddo
                  endif
                endif
              else
                !Case: UZA = 0.0
                if ( attenuations(n-1).ne.zero ) then
                  multiplier_dn(n) = Attenuations(n)
                  if ( do_columnwfs ) then
                     do q = 1, n_columnwfs
                        L_multiplier_dn(n,q) = LC_Attenuations(n-1,q)
                     enddo
                  endif
                endif
              endif
            endif

!  End layers and regular-PS formulation

         enddo
      endif

!  Enhanced PS multipliers and LOSTRANS: special case (nadir viewing)

      if ( do_enhanced_ps .and. doNadir(v) ) then
         do n = nlayers, 1, -1

!  LOS transmittance

            kn     = extinction(n)
            lostau = deltaus(n)
            if ( lostau .lt. cutoff ) lostrans_dn(n) = exp( - lostau )
            if ( do_columnwfs ) then
               if ( Qvary(n) ) then
                  do q = 1, Qnums(n)
                     L_lostau(q)        = L_deltaus(n,q)
                     L_lostrans_dn(n,q) = - L_lostau(q) * lostrans_dn(n)
                  enddo
               endif
            endif

!  Multipliers

            if ( layermask_dn(n) .and. n.le.nstart  ) then
               rdiff = radii(n-1) - radii(n)
               if ( n.eq.NCrit(v)) rdiff = radii(n-1) - RadCrit(v)
               sum = zero
               do j = 1, nfinedivs(n,v)
                  argum(j) = rdiff - xfine(n,j,v)
                  tran(j)  = exp ( - argum(j) * kn )
                  func(j)  = attenuations_fine(n,j) * tran(j)
                  sum = sum + func(j) * wfine(n,j,v)
               enddo
               multiplier_dn(n) = sum * kn
               if ( do_columnwfs ) then
                  do q = 1, n_columnwfs
                     L_sum = zero
                     do j = 1, nfinedivs(n,v)
                        L_tran = - argum(j) * L_extinction(n,q)
                        L_func = LC_attenuations_fine(n,j,q) * tran(j) + L_tran * func(j)
                        L_sum  = L_sum + L_func * wfine(n,j,v)
                     enddo
                     L_multiplier_dn(n,q)  = L_sum * kn + L_extinction(N,q) * sum
                  enddo
               endif
            endif

!  End layer loop and Nadir Enhanced PS case

         enddo
      endif

!  Enhanced PS multipliers and LOSTRANS: General case

      if ( do_enhanced_ps .and. .not. doNadir(v) ) then
         do n = nlayers, 1, -1

!  LOS transmittance

            cot_2 = cota(n-1,v) ; cot_1 = cota(n,v)
            cot_c = cot_1  ; if ( n.eq.NCrit(v) ) cot_c = CotCrit(v)
            kn = extinction(n) ;  ke = raycon(v) * kn ; cons = raycon(v) * ( cot_2 - cot_1 )
            tran_1 = kn * cons
            if ( tran_1 .lt. cutoff ) lostrans_dn(n) = exp ( - tran_1 )
            if ( do_columnwfs ) then
               if ( Qvary(n) ) then
                  do q = 1, Qnums(n)
                     L_lostau(q)        = L_extinction(n,q) * cons
                     L_lostrans_dn(n,q) = - L_lostau(q) * lostrans_dn(n)
                  enddo
               endif
            endif

!  multiplier

            if ( layermask_dn(n) .and. n.le.nstart  ) then
               sum = zero
               do j = 1, nfinedivs(n,v)
                  argum(j) = Raycon(v) * ( cotfine(n,j,v) - cot_c )
                  tran(j)  = exp ( - kn * argum(j) )
                  func(j)  = attenuations_fine(n,j) * csqfine(n,j,v) * tran(j)
                  sum      = sum + func(j) * wfine(n,j,v)
               enddo
               multiplier_dn(n) = sum * ke 
               if ( do_columnwfs ) then
                  do q = 1, n_columnwfs
                     L_sum = zero
                     do j = 1, nfinedivs(n,v)
                        L_tran = - argum(j) * L_extinction(n,q)
                        L_func = LC_attenuations_fine(n,j,q) * csqfine(n,j,v) * tran(j) + L_tran * func(j)
                        L_sum  = L_sum + L_func * wfine(n,j,v)
                     enddo
                     L_multiplier_dn(n,q)  = L_sum * ke + L_extinction(N,q) * Raycon(v) * sum
                  enddo
               endif
            endif

!  End layer loop and general Enhanced PS case

         enddo
      endif

!  Layer sources
!  -------------

      do n = nlayers, 1, -1
         if ( layermask_dn(n) .and. n.le.nstart  ) then
            do o1 = 1, nstokes
               shelp(o1) = dot_product(exactscat_dn(n,o1,1:ns),fluxvec(1:ns))
               sources_dn(n,o1) = shelp(o1) * multiplier_dn(n)
            enddo
            if ( do_columnwfs ) then
               do q = 1, n_columnwfs
                  do o1 = 1, nstokes
                     LC_sources_dn(n,o1,q) = shelp(o1) * L_multiplier_dn(n,q)
                     L_Shelp = dot_product(L_exactscat_dn(n,o1,1:ns,q),fluxvec(1:ns))
                     LC_sources_dn(n,o1,q) =  LC_sources_dn(n,o1,q) + L_Shelp * multiplier_dn(n)
                  enddo
               enddo
            endif
         endif
      enddo

!  Source function integration
!  ===========================

!  start recursion

      NC =  0
      CUMSOURCE_DN(NC,:) = zero
      NSTART = 1
      NUT_PREV = NSTART - 1

!  Main loop over all output optical depths
!     NLEVEL = Layer index for given optical depth
!     Cumulative source terms : Loop over layers working upwards from NSTART to level NUT,
!     Check for updating the recursion

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
         do o1 = 1, nstokes
            STOKES_DN(UTA,O1,V) = FLUX * CUMSOURCE_DN(NC,O1)
         enddo
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1 ; NUT_PREV = NUT
      ENDDO

!  Column Wfs

      if ( do_columnwfs ) then
         L_CUMSOURCE = zero
         NSTART = 1 ; NUT_PREV = NSTART - 1
         DO UTA = 1, N_USER_LEVELS
            NUT    = USER_LEVELS(UTA)
            DO N = NSTART, NUT
               NC = N
               do q = 1, n_columnwfs
                  L_cumsource(1:ns,q) = LC_SOURCES_DN(N,1:ns,Q)         + &
                                L_LOSTRANS_DN(N,Q) * CUMSOURCE_DN(NC-1,1:ns) + &
                                  LOSTRANS_DN(N)   * L_CUMSOURCE(1:ns,Q)
               enddo
            ENDDO
            do q = 1, n_columnwfs
               LC_JACOBIANS_DN(UTA,1:ns,V,Q) = FLUX * L_CUMSOURCE(1:ns,Q)
            enddo
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1 ;  NUT_PREV = NUT
         ENDDO
      endif

!  End GEometry Loop

   enddo

!  Finish

   return
end subroutine SSV_Integral_ILCS_DN


!

subroutine SSV_Integral_ILCS_UPDN &
   ( maxgeoms, maxlayers, maxfinelayers, maxmoments_input,               & ! Inputs (dimensioning)
     max_user_levels, max_atmoswfs, max_surfacewfs,                      & ! Inputs (dimensioning)
     do_upwelling, do_dnwelling, do_sunlight, do_deltam_scaling,         & ! Inputs (Flags - General)
     do_lambertian, do_PlanPar, do_regular_ps, do_enhanced_ps, doNadir, nstokes,       & ! Inputs (Flags - General)
     do_columnwfs, do_surfacewfs, n_columnwfs, n_surfacewfs, Lvarymoms,                & ! Inputs (control, Jacobian )
     ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,           & ! Inputs (Layers/Levels control)
     reflec, extinction, deltaus, omega, truncfac, Greekmat, flux, fluxvec,            & ! Inputs (Optical)
     LS_reflec, L_extinction, L_deltaus, L_omega, L_truncfac, L_greekmat,              & ! Inputs (Optical - Linearized)
     Mu0, Mu1, GenSpher_up, GenSpher_dn, Rotations_up, Rotations_dn, NCrit,            & ! Inputs (Geometry)
     RadCrit, CotCrit, xfine, wfine, csqfine, cotfine, Raycon, radii, cota,            & ! Inputs (Geometry)
     sunpaths_up, ntraverse_up, sunpaths_fine_up, ntraverse_fine_up,                   & ! Inputs (Geometry)
     sunpaths_dn, ntraverse_dn, sunpaths_fine_dn, ntraverse_fine_dn,                   & ! Inputs (Geometry)
     Stokes_up, Stokes_db, LC_Jacobians_up, LC_Jacobians_db, LS_Jacobians_db,          & ! Output
     Stokes_dn, LC_Jacobians_dn)                                                         ! Output

   Implicit none

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

!  General flags

   LOGICAL, Intent(in) :: DO_UPWELLING
   LOGICAL, Intent(in) :: DO_DNWELLING

!  flags

   LOGICAL, Intent(in) :: DO_SUNLIGHT
   LOGICAL, Intent(in) :: DO_DELTAM_SCALING
   LOGICAL, Intent(in) :: DO_LAMBERTIAN

   logical, Intent(in) :: DO_PLANPAR
   logical, Intent(in) :: DO_REGULAR_PS
   logical, Intent(in) :: DO_ENHANCED_PS
   logical, Intent(in) :: DONADIR(MAXGEOMS)

!  Jacobian Flags

   LOGICAL, Intent(in) :: do_surfacewfs
   LOGICAL, Intent(in) :: do_columnwfs

!  Layer and Level Control Numbers, Number of Moments

   integer, Intent(in) :: NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   integer, Intent(in) :: NGEOMS, NMOMENTS_INPUT, NSTOKES

   integer, Intent(in) :: N_USER_LEVELS
   integer, Intent(in) :: USER_LEVELS ( MAX_USER_LEVELS )

!  Jacobian control

   INTEGER, Intent(in) :: n_columnwfs
   INTEGER, Intent(in) :: n_surfacewfs
   LOGICAL, Intent(in) :: Lvarymoms (maxlayers,max_atmoswfs)

!  optical inputs
!  --------------

!  Atmosphere

   real(fpk), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(fpk), Intent(in) :: OMEGA       ( MAXLAYERS )
   real(fpk), Intent(in) :: TRUNCFAC    ( MAXLAYERS )
   real(fpk), Intent(in) :: GREEKMAT    ( MAXLAYERS,0:MAXMOMENTS_INPUT, 4, 4 )

!  Solar Flux and Surface reflectivity (Could be the albedo)

   real(fpk), Intent(in) :: REFLEC(4,4,MAXGEOMS), FLUX, FLUXVEC(4)

!  Linearized optical inputs

   real(fpk), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
   real(fpk), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
   real(fpk), Intent(in) :: L_OMEGA       ( MAXLAYERS, max_atmoswfs )
   real(fpk), Intent(in) :: L_TRUNCFAC    ( MAXLAYERS, max_atmoswfs )
   real(fpk), Intent(in) :: L_GREEKMAT    ( MAXLAYERS,0:MAXMOMENTS_INPUT, 4, 4, max_atmoswfs )
   real(fpk), Intent(in) :: LS_REFLEC     ( 4,4,MAXGEOMS, max_surfacewfs)

!  Ray constant, Cotangents, Critical layer
!    Mu0 = cos(theta_boa), required for surface term (both regular & enhanced)
!    Mu1 = cos(alpha_boa), required for the Regular PS only

   INTEGER  , Intent(in)  :: NCrit(maxgeoms)
   REAL(fpk), Intent(in)  :: Raycon(maxgeoms), cota(0:maxlayers,maxgeoms), radii(0:maxlayers)
   REAL(fpk), Intent(in)  :: Mu1(maxgeoms), Mu0(maxgeoms), RadCrit(maxgeoms), CotCrit(maxgeoms)

!  solar paths 

   INTEGER  , Intent(in)  :: ntraverse_up  (0:maxlayers,maxgeoms)
   REAL(fpk), Intent(in)  :: sunpaths_up   (0:maxlayers,maxlayers,maxgeoms)
   INTEGER  , Intent(in)  :: ntraverse_fine_up(maxlayers,maxfinelayers,maxgeoms)
   REAL(fpk), Intent(in)  :: sunpaths_fine_up (maxlayers,maxlayers,maxfinelayers,maxgeoms)

   INTEGER  , Intent(in)  :: ntraverse_dn  (0:maxlayers,maxgeoms)
   REAL(fpk), Intent(in)  :: sunpaths_dn   (0:maxlayers,maxlayers,maxgeoms)
   INTEGER  , Intent(in)  :: ntraverse_fine_dn(maxlayers,maxfinelayers,maxgeoms)
   REAL(fpk), Intent(in)  :: sunpaths_fine_dn (maxlayers,maxlayers,maxfinelayers,maxgeoms)

!  Generalized spherical functions.
!    Rotations(1-4)    = C1, S1, C2, S2

   REAL(fpk), Intent(in)  :: GenSpher_up(0:maxmoments_input,4,maxgeoms)
   REAL(fpk), Intent(in)  :: GenSpher_dn(0:maxmoments_input,4,maxgeoms)
   REAL(fpk), Intent(in)  :: Rotations_up(4,maxgeoms)
   REAL(fpk), Intent(in)  :: Rotations_dn(4,maxgeoms)

!  LOS Quadratures for Enhanced PS

   real(fpk), Intent(in)  :: xfine   (maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: wfine   (maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: csqfine (maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: cotfine (maxlayers,maxfinelayers,maxgeoms)

!  outputs
!  -------

   REAL(fpk), Intent(Out)  :: stokes_up     ( max_user_levels, 4, maxgeoms )
   REAL(fpk), Intent(Out)  :: stokes_db     ( max_user_levels, 4, maxgeoms )
   real(fpk), Intent(Out)  :: LC_Jacobians_up  ( max_user_levels, 4, maxgeoms, max_atmoswfs )
   real(fpk), Intent(Out)  :: LC_Jacobians_db  ( max_user_levels, 4, maxgeoms, max_atmoswfs )
   real(fpk), Intent(Out)  :: LS_Jacobians_db  ( max_user_levels, 4, maxgeoms, max_surfacewfs )

   REAL(fpk), Intent(Out)  :: stokes_dn     ( max_user_levels, 4, maxgeoms )
   real(fpk), Intent(Out)  :: LC_Jacobians_dn  ( max_user_levels, 4, maxgeoms, max_atmoswfs )

!  upwelling

   if ( do_upwelling  ) then
      call SSV_Integral_ILCS_UP &
   ( maxgeoms, maxlayers, maxfinelayers, maxmoments_input, & ! Inputs (dimensioning)
     max_user_levels, max_atmoswfs, max_surfacewfs,        & ! Inputs (dimensioning)
     do_sunlight, do_deltam_scaling, do_lambertian,        & ! Inputs (Flags - General)
     do_PlanPar, do_regular_ps, do_enhanced_ps, doNadir, nstokes,                  & ! Inputs (Flags - General)
     do_columnwfs, do_surfacewfs, n_columnwfs, n_surfacewfs, Lvarymoms,            & ! Inputs (control, Jacobian )
     ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,       & ! Inputs (control,  output)
     reflec, extinction, deltaus, omega, truncfac, greekmat, flux, fluxvec,        & ! Inputs (Optical)
     LS_reflec, L_extinction, L_deltaus, L_omega, L_truncfac, L_greekmat,          & ! Inputs (Optical - Linearized)
     Mu0, Mu1, GenSpher_up, Rotations_up, NCrit, xfine, wfine, csqfine, cotfine,   & ! Inputs (Geometry)
     Raycon, cota, sunpaths_up, ntraverse_up, sunpaths_fine_up, ntraverse_fine_up, & ! Inputs (Geometry)
     Stokes_up, Stokes_db, LC_Jacobians_up, LC_Jacobians_db, LS_Jacobians_db )       ! Output
   endif

!  Downwelling

   if ( do_dnwelling  ) then
      call SSV_Integral_ILCS_DN &
   ( maxgeoms, maxlayers, maxfinelayers, maxmoments_input,                      & ! Inputs (dimensioning)
     max_user_levels, max_atmoswfs, do_sunlight, do_deltam_scaling,             & ! Inputs (Flags - General)
     do_PlanPar, do_regular_ps, do_enhanced_ps, doNadir, nstokes,               & ! Inputs (Flags - General)
     do_columnwfs, n_columnwfs, Lvarymoms,                                      & ! Inputs (control, Jacobian )
     ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,    & ! Inputs (control,  output)
     extinction, deltaus, omega, truncfac, greekmat, flux, fluxvec,             & ! Inputs (Optical)
     L_extinction, L_deltaus, L_omega, L_truncfac, L_greekmat,                  & ! Inputs (Optical - Linearized)
     Mu1, GenSpher_dn, Rotations_dn, NCrit, RadCrit, CotCrit,                   & ! Inputs (Geometry)
     xfine, wfine, csqfine, cotfine, Raycon, radii, cota,                       & ! Inputs (Geometry)
     sunpaths_dn, ntraverse_dn, sunpaths_fine_dn, ntraverse_fine_dn,            & ! Inputs (Geometry)
     Stokes_dn, LC_Jacobians_dn )                                                 ! Output
   endif

!  Finish

   return
end subroutine SSV_Integral_ILCS_UPDN

!  End module

end module FO_VectorSS_RTCalcs_ILCS_m

