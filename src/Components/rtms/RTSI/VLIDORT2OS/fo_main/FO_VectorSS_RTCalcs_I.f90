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

module FO_VectorSS_RTCalcs_I_m

!  For given wavelength, routine calculates First-Order upwelling+downwelling Stokes vectors

!     (1) For the Atmospheric Solar Single-scatter and Surface Direct-Beam solar sources,

!  This is based on Precalculated Geometrical quantities and appropriate Optical properties.

!  This will perform Enhanced-PS calculations (incoming solar and outgoing LOS-path sphericity) 
!  This will perform Regular-PS  calculations (plane-parallel or incoming solar pseudo-spherical)

!  This is Versions 1-3, without Partials. Code is stand alone with no dependencies.
!    Version 1a, 01 December 2011, R. Spurr, RT Solutions Inc.
!    Version 1b, 13 February 2012, R. Spurr, RT Solutions Inc.
!    Version 2,  01 June     2012, R. Spurr, RT Solutions Inc.
!    Version 3,  29 October  2012, Extension to multiple geometries

!  For Solar sources, the subroutines are
!       SSV_Integral_I_UP   (Upwelling only)
!       SSV_Integral_I_DN   (Downwelling only)
!       SSV_Integral_I_UPDN (Upwelling and Downwelling)

!  All subroutines public

public

contains

subroutine SSV_Integral_I_UP &
   ( maxgeoms, maxlayers, maxfinelayers, maxmoments_input, max_user_levels,   & ! Inputs (dimension)
     do_sunlight, do_deltam_scaling, do_lambertian,                           & ! Inputs (Flags)
     do_PlanPar, do_regular_ps, do_enhanced_ps, doNadir, nstokes,             & ! Inputs (Flags)
     ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,  & ! Inputs (control)
     reflec, extinction, deltaus, omega, truncfac, greekmat, flux, fluxvec,   & ! Inputs (Optical)
     Mu0, Mu1, GenSpher, Rotations, NCrit, xfine, wfine, csqfine, cotfine,    & ! Inputs (Geometry)
     Raycon, cota, sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,        & ! Inputs (Geometry)
     stokes_up, stokes_db, cumsource_up )                                       ! Outputs

!  Stand-alone routine for Upwelling Solar-beam Single-scatter (SS)
!    computation of Stokes-vector. Inputs: geometry, spherical functions, optical properties.

!  This version, revised by R. Spurr, 01 June 2012
!   Extension to multiple geometries, 29 October 2012

!    No partials

   implicit none

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs
!  ======

!  Dimensions

   INTEGER, Intent(in) :: maxgeoms
   INTEGER, Intent(in) :: maxlayers
   INTEGER, Intent(in) :: maxfinelayers
   INTEGER, Intent(in) :: maxmoments_input
   INTEGER, Intent(in) :: max_user_levels

!  flags

   LOGICAL, Intent(in) ::  DO_SUNLIGHT
   LOGICAL, Intent(in) ::  DO_DELTAM_SCALING
   LOGICAL, Intent(in) ::  DO_LAMBERTIAN

   LOGICAL, Intent(in) ::  DO_PLANPAR
   LOGICAL, Intent(in) ::  DO_REGULAR_PS
   LOGICAL, Intent(in) ::  DO_ENHANCED_PS
   LOGICAL, Intent(in) ::  DONADIR(maxgeoms)

!  Numbers

   INTEGER, Intent(in) ::  NSTOKES
   INTEGER, Intent(in) ::  NGEOMS, NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   INTEGER, Intent(in) ::  NMOMENTS_INPUT

   INTEGER, Intent(in) ::  N_USER_LEVELS
   INTEGER, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  optical inputs
!  --------------

!  Atmosphere

   REAL(fpk), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   REAL(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
   REAL(fpk), Intent(in) :: OMEGA       ( MAXLAYERS )
   REAL(fpk), Intent(in) :: TRUNCFAC    ( MAXLAYERS )
   REAL(fpk), Intent(in) :: GREEKMAT    ( MAXLAYERS, 0:MAXMOMENTS_INPUT, 4, 4 )

!  Solar Flux and Surface reflectivity (Could be the albedo)

   REAL(fpk), Intent(in) :: REFLEC(4,4,maxgeoms), FLUX, FLUXVEC(4)

!  Geometrical inputs
!  ------------------

!  Ray constant, Cotangents, Critical layer
!    Mu0 = cos(theta_boa), required for surface term (both regular & enhanced)
!    Mu1 = cos(alpha_boa), required for the Regular PS only

   INTEGER  , Intent(in)  :: NCrit(maxgeoms)
   REAL(fpk), Intent(in)  :: Raycon(maxgeoms), cota(0:maxlayers,maxgeoms)
   REAL(fpk), Intent(in)  :: Mu0(maxgeoms), Mu1(maxgeoms)

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

   REAL(fpk), Intent(Out)  :: stokes_up     ( max_user_levels, 4, maxgeoms )
   REAL(fpk), Intent(Out)  :: stokes_db     ( max_user_levels, 4, maxgeoms )
   REAL(fpk), Intent(Out)  :: cumsource_up  ( 0:maxlayers,     4, maxgeoms )

!  LOCAL
!  -----

!  Attenuations

   REAL(fpk)  :: suntau            (0:maxlayers)
   REAL(fpk)  :: attenuations      (0:maxlayers)
   REAL(fpk)  :: attenuations_fine (maxlayers,maxfinelayers)

!  Scattering

   REAL(fpk)  :: tms (maxlayers)
   REAL(fpk)  :: exactscat_up (maxlayers,4,4)

!  Source function integration results

   REAL(fpk)  :: sources_up       ( maxlayers, 4 )
   REAL(fpk)  :: lostrans_up      ( maxlayers )
   REAL(fpk)  :: multiplier_up    ( maxlayers )

!  Regular_PS or plane-parallel flag

   LOGICAL    :: do_RegPSorPP

!  Help

   INTEGER    :: n, ns, uta, nstart, nc, nut, nut_prev, j, k, L, o1, o2, v
   LOGICAL    :: layermask_up(maxlayers)

   REAL(fpk)  :: sumd, help, sum, tran_1, tran, func, kn, ke, xjkn
   REAL(fpk)  :: cot_1, cot_2, lostau, factor1, factor2
   REAL(fpk)  :: cumsource_db(4), fmat(6), sum23, dif23
   REAL(fpk)  :: help3c1, help3s1, help4c1, help4s1, m4

   REAL(fpk), parameter  :: cutoff = 88.0_fpk
   REAL(fpk), parameter  :: zero   = 0.0_fpk
   REAL(fpk), parameter  :: one    = 1.0_fpk

!  Zero the output

   CUMSOURCE_UP  = zero ; STOKES_UP  = zero ; STOKES_DB    = zero

!  Regular_PS or plane-parallel flag

   do_RegPSorPP = (do_regular_ps .or. do_PlanPar)

!  Bookkeeping

   ns = nstokes ; fmat = zero
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

      lostrans_up = zero ; sources_up = zero ; exactscat_up = zero ; multiplier_up = zero

!  Scattering functions
!  --------------------

!  Scalar only

      if ( nstokes .eq. 1 ) then
         do n = 1, nlayers
            if ( layermask_up(n) ) then
              sum = zero
               do L = 0, nmoments_input
                  sum = sum + GenSpher(L,1,V) * Greekmat(n,L,1,1)
               enddo
               exactscat_up(n,1,1) = sum * tms(n)
            endif
         enddo
      endif

!  Vector with Sunlight

      if ( nstokes .gt. 1 .and. do_sunlight ) then
         do n = 1, nlayers
            if ( layermask_up(n) ) then
               fmat(1:2) = zero
               do L = 0, nmoments_input
                  fmat(1) = fmat(1) + Genspher(L,1,v) * greekmat(n,L,1,1)
                  fmat(2) = fmat(2) + Genspher(L,2,v) * greekmat(n,L,1,2)
               enddo
               exactscat_up(n,1,1) = + fmat(1)
               exactscat_up(n,2,1) = - fmat(2) * Rotations(3,v)
               exactscat_up(n,3,1) = + fmat(2) * Rotations(4,v)
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
                  fmat(1) = fmat(1) + Genspher(L,1,v) * greekmat(n,L,1,1)
                  fmat(2) = fmat(2) + Genspher(L,2,v) * greekmat(n,L,1,2)
                  sum23 = greekmat(n,L,2,2) + greekmat(n,L,3,3)
                  dif23 = greekmat(n,L,2,2) - greekmat(n,L,3,3)
                  fmat(3) = fmat(3) + Genspher(L,3,v) * sum23
                  fmat(4) = fmat(4) + Genspher(L,4,v) * dif23
               enddo
               fmat(3) = ( fmat(3) + fmat(4) ) * 0.5_fpk
               fmat(4) = ( fmat(3) - fmat(4) )
               if ( nstokes.eq.4) then
                  do L = 0, nmoments_input
                     fmat(5) = fmat(5) + Genspher(L,2,v) * greekmat(n,L,3,4)
                     fmat(6) = fmat(6) + Genspher(L,1,v) * greekmat(n,L,4,4)
                  enddo
               endif
               help3c1 = fmat(3) * Rotations(1,v)
               help3s1 = fmat(3) * Rotations(2,v)
               help4c1 = fmat(4) * Rotations(1,v)
               help4s1 = fmat(4) * Rotations(2,v)
               exactscat_up(n,1,1) = + fmat(1)
               exactscat_up(n,2,1) = - fmat(2) * Rotations(3,v)
               exactscat_up(n,3,1) = + fmat(2) * Rotations(4,v)
               exactscat_up(n,1,2) = + fmat(2) * Rotations(1,v)
               exactscat_up(n,1,3) = - fmat(2) * Rotations(2,v)
               exactscat_up(n,2,2) = + help3c1 * Rotations(3,v) - help4s1 * Rotations(4,v)
               exactscat_up(n,2,3) = - help3s1 * Rotations(3,v) - help4c1 * Rotations(4,v)
               exactscat_up(n,3,2) = + help3c1 * Rotations(4,v) + help4s1 * Rotations(3,v)
               exactscat_up(n,3,3) = - help3s1 * Rotations(4,v) + help4c1 * Rotations(3,v)
               if ( nstokes .eq. 4 ) then
                  exactscat_up(n,2,4) = - fmat(5) * Rotations(4,v) 
                  exactscat_up(n,4,2) = - fmat(5) * Rotations(2,v) 
                  exactscat_up(n,3,4) = + fmat(5) * Rotations(3,v) 
                  exactscat_up(n,4,3) = - fmat(5) * Rotations(1,v) 
                  exactscat_up(n,4,4) = + fmat(6)
               endif
               exactscat_up(n,1:ns,1:ns) = tms(n)*exactscat_up(n,1:ns,1:ns)
            endif
         enddo
      endif

!  Attenuations
!  ============

!  Initialize, only to layer Ncrit if applicable

      Attenuations = zero ; Attenuations_fine = zero ; Suntau = zero 
      nstart = nlayers ; if (Ncrit(v).ne.0) nstart = nCrit(v)

!  Attenuations to End points (including TOA). All representations
!    MUST go all the way to NLAYERS (surface term required)

      do n = 0, nlayers
         sumd = zero
         do k = 1, ntraverse(n,v)
            sumd = sumd + extinction(k) * sunpaths(n,k,v)
         enddo
         suntau(n) = sumd
         If (sumd .lt. cutoff ) Attenuations(n) = exp( - sumd )
!         if ( n.gt.87 ) write(*,*)n,sumd,Attenuations(n)
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
                  sumd = ZERO
                  do k = 1, ntraverse_fine(n,j,v)
                     sumd = sumd + extinction(k) * sunpaths_fine(n,k,j,v)
                  enddo
                  if (sumd .lt. cutoff ) Attenuations_fine(n,j) = exp( - sumd )
               enddo
            endif
         enddo
      endif

!  Layer integrated Solar sources
!  ==============================

!  Plane-parallel or Regular-PS, multipliers and LOSTRANS

      if ( do_RegPSorPP ) then
         if ( Mu1(v) .eq. zero ) then
            do n = nlayers, 1, -1
               if ( layermask_up(n) .and. n.le.nstart ) then
                  factor1 = Attenuations(n-1) - Attenuations(n)*lostrans_up(n)
                  multiplier_up(n) = factor1 
               endif
            enddo 
         else
            do n = nlayers, 1, -1
               lostau = deltaus(n) / Mu1(v)
               if ( lostau .lt. cutoff ) lostrans_up(n) = exp( - lostau )
               if ( layermask_up(n) .and. n.le.nstart ) then
                  factor1 = Attenuations(n-1) - Attenuations(n)*lostrans_up(n)
                  factor2 = (suntau(n) - suntau(n-1))/lostau
                  multiplier_up(n) = factor1 / (factor2 + one)
               endif
!               write(76,*)n,lostau, lostrans_up(n),multiplier_up(n)
            enddo
         endif
      endif

!  Enhanced PS multiplier: special case (nadir viewing)

      if ( do_enhanced_ps .and. doNadir(v) ) then
         do n = nlayers, 1, -1
            kn = extinction(n)
            lostrans_up(n)  = exp ( - deltaus(n))
            if ( layermask_up(n) .and. n.le.nstart ) then
               sum = zero
               do j = 1, nfinedivs(n,v)
                  xjkn = xfine(n,j,v) * kn
                  func = attenuations_fine(n,j) * exp ( - xjkn )
                  sum = sum + func * wfine(n,j,v)
               enddo
               multiplier_up(n) = sum * kn
            endif
         enddo
      endif

!  Enhanced PS multiplier: General case

      if ( do_enhanced_ps .and. .not. doNadir(v) ) then
         do n = nlayers, 1, -1
            cot_2 = cota(n-1,v) ; cot_1 = cota(n,v)
            kn = extinction(n) ;  ke = raycon(v) * kn
            tran_1 = exp ( - ke * ( cot_2 - cot_1 ) )
            lostrans_up(n) = tran_1
            if ( n.le.nstart ) then
               sum = zero
               do j = 1, nfinedivs(n,v)
                  tran = exp ( - ke * ( cot_2 - cotfine(n,j,v) ) )
                  func = attenuations_fine(n,j) * csqfine(n,j,v) * tran
                  sum  = sum + func * wfine(n,j,v)
               enddo
               multiplier_up(n) = sum * ke
            endif
         enddo
      endif

!  Layer sources

      do n = nlayers, 1, -1
         if ( layermask_up(n) .and. n.le.nstart  ) then
            if ( do_sunlight ) then
               do o1 = 1, nstokes
                  sources_up(n,o1) = exactscat_up(n,o1,1) * multiplier_up(n) * fluxvec(1)
               enddo
            else
               do o1 = 1, nstokes
                  sources_up(n,o1) = dot_product(exactscat_up(n,o1,1:ns),fluxvec(1:ns)) * multiplier_up(n)
               enddo
            endif
         endif
      enddo

!  Source function integration
!  ===========================

!  initialize recursion ( For Direct Beam, use PI.mu0.R.Atten )

      NC =  0
      CUMSOURCE_UP(NC,:,V) = zero
      CUMSOURCE_DB         = zero
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  Surface term

      M4 = 4.0d0 * MU0(V)

      if ( DO_LAMBERTIAN ) then
         CUMSOURCE_DB(1) = M4 * REFLEC(1,1,V) * attenuations(nlayers) * fluxvec(1)
      else
         do o1 = 1, nstokes
            sum = zero
            do o2 = 1, nstokes
               sum = sum + reflec(o1,o2,v) * fluxvec(o2)
            enddo
!mick fix 7/22/2013 - multiply by attenuation
            !CUMSOURCE_DB(o1) = M4 * sum
            CUMSOURCE_DB(o1) = M4 * attenuations(nlayers) * sum
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
               CUMSOURCE_DB(O1)       = LOSTRANS_UP(N) * CUMSOURCE_DB(O1)
               CUMSOURCE_UP(NC,O1,V)  = LOSTRANS_UP(N) * CUMSOURCE_UP(NC-1,O1,V) + SOURCES_UP(N,O1)
            enddo
         ENDDO
         do o1 = 1, nstokes
            STOKES_UP(UTA,O1,V) = FLUX * CUMSOURCE_UP(NC,O1,V)
            STOKES_DB(UTA,O1,V) = FLUX * CUMSOURCE_DB(O1)
         enddo
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
         NUT_PREV = NUT
      ENDDO

!  Finish geometry loop

   enddo

!  Finish

   return
end subroutine SSV_Integral_I_UP

!

subroutine SSV_Integral_I_DN &
   ( maxgeoms, maxlayers, maxfinelayers, maxmoments_input, max_user_levels, do_sunlight, & ! Inputs (dimension)
     do_deltam_scaling, do_PlanPar, do_regular_ps, do_enhanced_ps, doNadir,              & ! Inputs (Flags)
     ngeoms, nstokes, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,    & ! Inputs (control)
     extinction, deltaus, omega, truncfac, greekmat, flux, fluxvec,                      & ! Inputs (Optical)
     Mu1, GenSpher, Rotations, NCrit, RadCrit, CotCrit,                                  & ! Inputs (Geometry)
     xfine, wfine, csqfine, cotfine, Raycon, radii, cota,                                & ! Inputs (Geometry)
     sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,                                 & ! Inputs (Geometry)
     stokes_dn, cumsource_dn )                                                             ! Outputs

!  Stand alone routine for SS field with Solar sources alone
!    No partials

   use FO_Taylor_m

   implicit none

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs
!  ======

!  Dimensions

   INTEGER, Intent(in) :: maxgeoms
   INTEGER, Intent(in) :: maxlayers
   INTEGER, Intent(in) :: maxfinelayers
   INTEGER, Intent(in) :: maxmoments_input
   INTEGER, Intent(in) :: max_user_levels

!  flags

   LOGICAL, Intent(in) ::  DO_SUNLIGHT
   LOGICAL, Intent(in) ::  DO_DELTAM_SCALING

   LOGICAL, Intent(in) ::  DO_PLANPAR
   LOGICAL, Intent(in) ::  DO_REGULAR_PS
   LOGICAL, Intent(in) ::  DO_ENHANCED_PS
   LOGICAL, Intent(in) ::  DONADIR(MAXGEOMS)

!  Numbers

   INTEGER, Intent(in) ::  NSTOKES
   INTEGER, Intent(in) ::  NGEOMS, NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   INTEGER, Intent(in) ::  NMOMENTS_INPUT

   INTEGER, Intent(in) ::  N_USER_LEVELS
   INTEGER, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  optical inputs
!  --------------

!  Atmosphere

   REAL(fpk), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   REAL(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
   REAL(fpk), Intent(in) :: OMEGA       ( MAXLAYERS )
   REAL(fpk), Intent(in) :: TRUNCFAC    ( MAXLAYERS )
   REAL(fpk), Intent(in) :: GREEKMAT    ( MAXLAYERS,0:MAXMOMENTS_INPUT, 4, 4 )

!  Solar Flux 

   REAL(fpk), Intent(in) :: FLUX, FLUXVEC(4)

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

   REAL(fpk), Intent(Out)  :: stokes_dn     ( max_user_levels, 4, maxgeoms )
   REAL(fpk), Intent(Out)  :: cumsource_dn  ( 0:maxlayers,     4, maxgeoms )

!  LOCAL
!  -----

!  Attenuations

   REAL(fpk)  :: suntau            (0:maxlayers)
   REAL(fpk)  :: attenuations      (0:maxlayers)
   REAL(fpk)  :: attenuations_fine (maxlayers,maxfinelayers)

!  Scattering

   REAL(fpk)  :: tms (maxlayers)
   REAL(fpk)  :: exactscat_dn (maxlayers,4,4)

!  Source function integration results

   REAL(fpk)  :: sources_dn       ( maxlayers, 4 )
   REAL(fpk)  :: lostrans_dn      ( maxlayers )
   REAL(fpk)  :: multiplier_dn    ( maxlayers )

!  Regular_PS or plane-parallel flag

   LOGICAL    :: do_RegPSorPP

!  Help

   INTEGER    :: n, ns, uta, nstart, nc, nut, nut_prev, j, k, L, O1, O2, v
   LOGICAL    :: layermask_dn(maxlayers)
   REAL(fpk)  :: sumd, help, sum, tran_1, tran, func, kn, ke, xjkn, trand
   REAL(fpk)  :: cot_1, cot_2, factor1, factor2, factor3, lostau, rdiff, cot_c
   REAL(fpk)  :: fmat(6), sum23, dif23, help3c1, help3s1, help4c1, help4s1

   REAL(fpk)  :: eps, lospath, mult, term2

   REAL(fpk), parameter  :: cutoff = 88.0_fpk
   REAL(fpk), parameter  :: zero   = 0.0_fpk
   REAL(fpk), parameter  :: one    = 1.0_fpk

!  Zero the output

   CUMSOURCE_DN = zero ; STOKES_DN  = zero

!  Regular_PS or plane-parallel flag

   do_RegPSorPP = (do_regular_ps .or. do_PlanPar)

!  Bookkeeping

   ns = nstokes ; fmat = zero
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

      lostrans_dn  = zero ; sources_dn = zero ; exactscat_dn = zero ; multiplier_dn = zero

!  Scattering functions
!  --------------------

!  Scalar only

      if ( nstokes .eq. 1 ) then
         do n = 1, nlayers
            if ( layermask_dn(n) ) then
              sum = zero
               do L = 0, nmoments_input
                  sum = sum + GenSpher(L,1,v) * Greekmat(n,L,1,1)
               enddo
               exactscat_dn(n,1,1) = sum * tms(n)
            endif
         enddo
      endif

!  Vector with Sunlight

      if ( nstokes .gt. 1 .and. do_sunlight ) then
         do n = 1, nlayers
            if ( layermask_dn(n) ) then
               fmat(1:2) = zero
               do L = 0, nmoments_input
                  fmat(1) = fmat(1) + Genspher(L,1,v) * greekmat(n,L,1,1)
                  fmat(2) = fmat(2) + Genspher(L,2,v) * greekmat(n,L,1,2)
               enddo
               exactscat_dn(n,1,1) = + fmat(1)
               exactscat_dn(n,2,1) = - fmat(2) * Rotations(3,v)
               exactscat_dn(n,3,1) = + fmat(2) * Rotations(4,v)
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
                  fmat(1) = fmat(1) + Genspher(L,1,v) * greekmat(n,L,1,1)
                  fmat(2) = fmat(2) + Genspher(L,2,v) * greekmat(n,L,1,2)
                  sum23 = greekmat(n,L,2,2) + greekmat(n,L,3,3)
                  dif23 = greekmat(n,L,2,2) - greekmat(n,L,3,3)
                  fmat(3) = fmat(3) + Genspher(L,3,v) * sum23
                  fmat(4) = fmat(4) + Genspher(L,4,v) * dif23
               enddo
               fmat(3) = ( fmat(3) + fmat(4) ) * 0.5_fpk
               fmat(4) = ( fmat(3) - fmat(4) )
               if ( nstokes.eq.4) then
                  do L = 0, nmoments_input
                     fmat(5) = fmat(5) + Genspher(L,2,v) * greekmat(n,L,3,4)
                     fmat(6) = fmat(6) + Genspher(L,1,v) * greekmat(n,L,4,4)
                  enddo
               endif
               help3c1 = fmat(3) * Rotations(1,v)
               help3s1 = fmat(3) * Rotations(2,v)
               help4c1 = fmat(4) * Rotations(1,v)
               help4s1 = fmat(4) * Rotations(2,v)
               exactscat_dn(n,1,1) = + fmat(1)
               exactscat_dn(n,2,1) = - fmat(2) * Rotations(3,v)
               exactscat_dn(n,3,1) = + fmat(2) * Rotations(4,v)
               exactscat_dn(n,1,2) = + fmat(2) * Rotations(1,v)
               exactscat_dn(n,1,3) = - fmat(2) * Rotations(2,v)
               exactscat_dn(n,2,2) = + help3c1 * Rotations(3,v) - help4s1 * Rotations(4,v)
               exactscat_dn(n,2,3) = - help3s1 * Rotations(3,v) - help4c1 * Rotations(4,v)
               exactscat_dn(n,3,2) = + help3c1 * Rotations(4,v) + help4s1 * Rotations(3,v)
               exactscat_dn(n,3,3) = - help3s1 * Rotations(4,v) + help4c1 * Rotations(3,v)
               if ( nstokes .eq. 4 ) then
                  exactscat_dn(n,2,4) = - fmat(5) * Rotations(4,v) 
                  exactscat_dn(n,4,2) = - fmat(5) * Rotations(2,v) 
                  exactscat_dn(n,3,4) = + fmat(5) * Rotations(3,v) 
                  exactscat_dn(n,4,3) = - fmat(5) * Rotations(1,v) 
                  exactscat_dn(n,4,4) = + fmat(6)
               endif
               exactscat_dn(n,1:ns,1:ns) = tms(n)*exactscat_dn(n,1:ns,1:ns)
            endif
         enddo
      endif

!  Attenuations and Solar solutions
!  ================================

!  Initialize, only to layer Ncrit if applicable

      Attenuations = ZERO ; Attenuations_fine = ZERO
      nstart = nlayers ; if (Ncrit(v).ne.0) nstart = nCrit(v)

!  Attenuations to End points (including TOA). All representations
!    MUST go all the way to NLAYERS (surface term required)

      do n = 0, nlayers
         sumd = zero
         do k = 1, ntraverse(n,v)
            sumd = sumd + extinction(k) * sunpaths(n,k,v)
         enddo
         suntau(n) = sumd
         If (sumd .lt. cutoff ) Attenuations(n) = exp( - sumd )
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
               enddo
            endif
         enddo
      endif

!  Layer integrated Solar sources
!  ==============================

!  Plane-parallel or Regular-PS, multipliers and LOSTRANS

      if ( do_RegPSorPP ) then
         if( Mu1(v) .eq. zero ) then
            !Case: UZA = 0.0
            do n = 1, nlayers
               if ( layermask_dn(n).and. n.le.nstart ) then
                  if ( attenuations(n-1).ne.zero ) then
                     factor1 = Attenuations(n)
                     multiplier_dn(n) = factor1
                  endif
               endif
            enddo
         else
            !Case: UZA /= 0.0
            do n = 1, nlayers
               lostau = deltaus(n) / Mu1(v)
               if ( lostau .lt. cutoff ) lostrans_dn(n) = exp( - lostau )
               if ( layermask_dn(n) .and. n.le.nstart ) then
!mick fix 7/18/2013 - added "taylor small" section
                  lospath = lostau/extinction(n)
                  eps     = sunpaths(n,n,v) - lospath
                  if ( abs(eps) .lt. TAYLOR_SMALL ) then
                    term2 = exp( - extinction(n) * sunpaths(n,n,v) )
                    call taylor_series_1 (eps, extinction(n), term2, mult)
                    multiplier_dn(n) = Attenuations(n-1) * lospath * mult
                  else
                    factor1 = lostrans_dn(n) * Attenuations(n-1) - Attenuations(n)
                    factor2 = (suntau(n) - suntau(n-1))/lostau
                    factor3 = factor2 - one
                    multiplier_dn(n) = factor1 / factor3
                  end if
               endif
            enddo
         endif
      endif

!  Enhanced PS multipliers and LOSTRANS: special case (nadir viewing)

      if ( do_enhanced_ps .and. doNadir(v) ) then
         do n = nlayers, 1, -1
            kn = extinction(n)
            lostrans_dn(n)  = exp ( - deltaus(n))
            rdiff = radii(n-1) - radii(n) ; if ( n.eq.NCrit(v)) rdiff = radii(n-1) - RadCrit(v)
            trand = one ; if ( n.eq.NCrit(v)) trand = exp ( -kn * (RadCrit(v) -radii(n) ) )
            if ( layermask_dn(n) .and. n.le.nstart  ) then
               sum = zero
               do j = 1, nfinedivs(n,v)
                  xjkn = ( rdiff - xfine(n,j,v) ) * kn
                  func = attenuations_fine(n,j) * exp ( - xjkn )
                  sum = sum + func * wfine(n,j,v)
               enddo
               multiplier_dn(n) = sum * kn
            endif
         enddo
      endif

!  Enhanced PS multipliers and LOSTRANS: General case

      if ( do_enhanced_ps .and. .not. doNadir(v) ) then
         do n = nlayers, 1, -1
            cot_2 = cota(n-1,v) ; cot_1 = cota(n,v)
            cot_c = cot_1  ; if ( n.eq.NCrit(v) ) cot_c = CotCrit(v)
            kn = extinction(n) ;  ke = raycon(v) * kn
            trand = one  ; if ( n.eq.NCrit(v) ) trand = exp ( - ke * ( CotCrit(v) - cot_1 ) )
            tran_1 = exp ( - ke * ( cot_2 - cot_1 ) )
            lostrans_dn(n) = tran_1
            if ( layermask_dn(n) .and. n.le.nstart  ) then
               sum = zero
               do j = 1, nfinedivs(n,v)
                  tran = exp ( - ke * ( cotfine(n,j,v) - cot_c ) )   !  Down
!                 tran = exp ( - ke * ( cot_2 - cotfine(n,j) ) )   !  Up
                  func = Attenuations_fine(n,j) * csqfine(n,j,v) * tran
                  sum  = sum + func * wfine(n,j,v)
               enddo
               multiplier_dn(n) = sum * ke * trand
            endif
         enddo
      endif

!  Layer sources

      do n = nlayers, 1, -1
         if ( layermask_dn(n) .and. n.le.nstart ) then
            do o1 = 1, nstokes
               sum = zero
               do o2 = 1, nstokes
                  sum = sum + exactscat_dn(n,o1,o2) * fluxvec(o2)
               enddo
               sources_dn(n,o1) = sum * multiplier_dn(n)
            enddo
         endif
      enddo

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
         do o1 = 1, nstokes
            STOKES_DN(UTA,O1,v) = FLUX * CUMSOURCE_DN(NC,O1,v)
         enddo
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
         NUT_PREV = NUT
      ENDDO

!  Finish geometry loop

   enddo

!  Finish

   return
end subroutine SSV_Integral_I_DN

!

subroutine SSV_Integral_I_UPDN   &
   ( maxgeoms, maxlayers, maxfinelayers, maxmoments_input, max_user_levels,            & ! Inputs (dimensioning)
     do_upwelling, do_dnwelling, do_deltam_scaling, do_sunlight,                       & ! Inputs (Flags)
     do_lambertian, do_Planpar, do_regular_ps, do_enhanced_ps, doNadir,                & ! Inputs (Flags)
     ngeoms, nstokes, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,  & ! Inputs (control output)
     reflec, extinction, deltaus, omega, truncfac, Greekmat, flux, fluxvec,            & ! Inputs (Optical)
     Mu0, Mu1, GenSpher_up, GenSpher_dn, Rotations_up, Rotations_dn, NCrit,            & ! Inputs (Geometry)
     RadCrit, CotCrit, xfine, wfine, csqfine, cotfine, Raycon, radii, cota,            & ! Inputs (Geometry)
     sunpaths_up, ntraverse_up, sunpaths_fine_up, ntraverse_fine_up,                   & ! Inputs (Geometry)
     sunpaths_dn, ntraverse_dn, sunpaths_fine_dn, ntraverse_fine_dn,                   & ! Inputs (Geometry)
     stokes_up, stokes_db, cumsource_up, stokes_dn, cumsource_dn )                       ! Outputs

!  Stand alone routine for SS field with Solar sources alone
!    No partials

   implicit none

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs
!  ======

!  Dimensions

   INTEGER, Intent(in) :: maxgeoms
   INTEGER, Intent(in) :: maxlayers
   INTEGER, Intent(in) :: maxfinelayers
   INTEGER, Intent(in) :: maxmoments_input
   INTEGER, Intent(in) :: max_user_levels

!  flags

   LOGICAL, Intent(in) ::  DO_UPWELLING
   LOGICAL, Intent(in) ::  DO_DNWELLING
   LOGICAL, Intent(in) ::  DO_SUNLIGHT
   LOGICAL, Intent(in) ::  DO_DELTAM_SCALING
   LOGICAL, Intent(in) ::  DO_LAMBERTIAN

   LOGICAL, Intent(in) ::  DO_PLANPAR
   LOGICAL, Intent(in) ::  DO_REGULAR_PS
   LOGICAL, Intent(in) ::  DO_ENHANCED_PS
   LOGICAL, Intent(in) ::  DONADIR(maxgeoms)

!  Numbers

   INTEGER, Intent(in) ::  NGEOMS, NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   INTEGER, Intent(in) ::  NMOMENTS_INPUT, NSTOKES

   INTEGER, Intent(in) ::  N_USER_LEVELS
   INTEGER, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  optical inputs
!  --------------

!  Atmosphere

   REAL(fpk), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   REAL(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
   REAL(fpk), Intent(in) :: OMEGA       ( MAXLAYERS )
   REAL(fpk), Intent(in) :: TRUNCFAC    ( MAXLAYERS )
   REAL(fpk), Intent(in) :: GREEKMAT    ( MAXLAYERS, 0:MAXMOMENTS_INPUT, 4, 4 )

!  Solar Flux and Surface reflectivity (Could be the albedo)

   REAL(fpk), Intent(in) :: REFLEC(4,4,MAXGEOMS), FLUX, FLUXVEC(4)

!  Geometrical inputs
!  ------------------

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

   REAL(fpk), Intent(in)  :: xfine   (maxlayers,maxfinelayers,maxgeoms)
   REAL(fpk), Intent(in)  :: wfine   (maxlayers,maxfinelayers,maxgeoms)
   REAL(fpk), Intent(in)  :: csqfine (maxlayers,maxfinelayers,maxgeoms)
   REAL(fpk), Intent(in)  :: cotfine (maxlayers,maxfinelayers,maxgeoms)

!  outputs
!  -------

   REAL(fpk), Intent(Out)  :: stokes_up     ( max_user_levels, 4, maxgeoms )
   REAL(fpk), Intent(Out)  :: stokes_db     ( max_user_levels, 4, maxgeoms )
   REAL(fpk), Intent(Out)  :: cumsource_up  ( 0:maxlayers, 4, maxgeoms )
   REAL(fpk), Intent(Out)  :: stokes_dn     ( max_user_levels, 4, maxgeoms )
   REAL(fpk), Intent(Out)  :: cumsource_dn  ( 0:maxlayers, 4, maxgeoms )

!  Upwelling
!  ---------

   if ( do_upwelling ) then
       call SSV_Integral_I_UP &
   ( maxgeoms, maxlayers, maxfinelayers, maxmoments_input, max_user_levels,           & ! Inputs (dimension)
     do_sunlight, do_deltam_scaling, do_lambertian,                                   & ! Inputs (Flags)
     do_PlanPar, do_regular_ps, do_enhanced_ps, doNadir,                              & ! Inputs (Flags)
     ngeoms, nstokes, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels, & ! Inputs (control)
     reflec, extinction, deltaus, omega, truncfac, greekmat, flux, fluxvec,           & ! Inputs (Optical)
     Mu0, Mu1, GenSpher_up, Rotations_up, NCrit, xfine, wfine, csqfine, cotfine,      & ! Inputs (Geometry)
     Raycon, cota, sunpaths_up, ntraverse_up, sunpaths_fine_up, ntraverse_fine_up,    & ! Inputs (Geometry)
     stokes_up, stokes_db, cumsource_up )                                               ! Outputs
   endif

   if ( do_dnwelling ) then
       call SSV_Integral_I_DN &
   ( maxgeoms, maxlayers, maxfinelayers, maxmoments_input, max_user_levels, do_sunlight, & ! Inputs (dimension)
     do_deltam_scaling, do_PlanPar, do_regular_ps, do_enhanced_ps, doNadir,              & ! Inputs (Flags)
     ngeoms, nstokes, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,    & ! Inputs (control)
     extinction, deltaus, omega, truncfac, greekmat, flux, fluxvec,                      & ! Inputs (Optical)
     Mu1, GenSpher_dn, Rotations_dn, NCrit, RadCrit, CotCrit,                            & ! Inputs (Geometry)
     xfine, wfine, csqfine, cotfine, Raycon, radii, cota,                                & ! Inputs (Geometry)
     sunpaths_dn, ntraverse_dn, sunpaths_fine_dn, ntraverse_fine_dn,                     & ! Inputs (Geometry)
     stokes_dn, cumsource_dn )                                                             ! Outputs
   endif

!  Finish

   return
end subroutine SSV_Integral_I_UPDN


!  End module

end module FO_VectorSS_RTCalcs_I_m
