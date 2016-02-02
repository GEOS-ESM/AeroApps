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
! #     FIRST-ORDER THERMAL MODEL (Direct term)             #
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

module FO_Thermal_RTCalcs_ILPS_m

!  For a given wavelength, this routine will calculate upwelling and downwelling
!  First Order Intensities(I), and any number of LPS Jacobians (profile/surface)

!     (1) For the Atmospheric and Surface Direct Thermal Emission (DTE) sources.

!  This is based on Precalculated Geometrical quantities and appropriate Optical properties.

!  This will perform Enhanced-PS calculations (outgoing LOS-path sphericity) 
!  This will perform Regular-PS  calculations (plane-parallel LOS-path)

!  This is Versions 1-3, without Partials. Code is stand alone with no dependencies.
!    Version 1a, 01 December 2011, R. Spurr, RT Solutions Inc.
!    Version 1b, 02 February 2012, R. Spurr, RT Solutions Inc.
!    Version 2 , 01 June     2012, R. Spurr, RT Solutions Inc.
!    Version 3 , 19 December 2012, Extension to Multiple Geometries, LCS separation

!  For Thermal Emission sources, the subroutines are
!       DTE_Integral_ILPS_UP   (Upwelling only)
!       DTE_Integral_ILPS_DN   (Downwelling only)
!       DTE_Integral_ILPS_UPDN (Upwelling and Downwelling)

!  All subroutines public

public

contains

subroutine DTE_Integral_ILPS_UP &
   ( maxgeoms, maxlayers, maxfinelayers, max_user_levels,               & ! Inputs (dimensioning)
     max_atmoswfs, max_surfacewfs,                                      & ! Inputs (dimensioning)
     Do_Thermset, do_deltam_scaling, do_PlanPar,                        & ! Inputs (Flags)
     do_regular_ps, do_enhanced_ps, doNadir,                            & ! Inputs (Flags)
     do_profilewfs, do_surfacewfs, Lvaryflags, Lvarynums, n_surfacewfs, & ! Inputs (Control, Jacobians)
     ngeoms, nlayers, nfinedivs, n_user_levels, user_levels,            & ! Inputs (control, output)
     bb_input, surfbb, user_emissivity, LS_user_emissivity,             & ! Inputs (Thermal)
     extinction, deltaus, omega, truncfac,                              & ! Inputs (Optical - Regular)
     L_extinction, L_deltaus, L_omega, L_truncfac,                      & ! Inputs (Optical - Linearized)
     Mu1, NCrit, Raycon, cota, xfine, wfine, csqfine, cotfine,          & ! Inputs (Geometry)
     intensity_dta_up, intensity_dts, LP_Jacobians_dta_up,              & ! Outputs
     LP_Jacobians_dts_up, LS_Jacobians_dts, tcom1, L_tcom1  )             ! Output

!  Stand alone routine for Upwelling Direct-thermal-emission (DTE)
!    computation of Radiances and LPS Jacobians. Inputs: geometry, optical properties, Planck functions, emissivity

!  This version, revised by R. Spurr, 01 June 2012
!   Extension to multiple geometries, 19 December 2012

!    No partials

   implicit none         

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs
!  ======

!  Dimensions

   integer, Intent(in) :: maxgeoms
   integer, Intent(in) :: maxlayers
   integer, Intent(in) :: maxfinelayers
   integer, Intent(in) :: max_user_levels
   INTEGER, Intent(in) :: max_atmoswfs
   INTEGER, Intent(in) :: max_surfacewfs

!  Thermal setup flag (for TCOM1)

   LOGICAL, Intent(inout) :: Do_Thermset

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
   integer, Intent(in) :: NGEOMS

   integer, Intent(in) :: N_USER_LEVELS
   integer, Intent(in) :: USER_LEVELS ( MAX_USER_LEVELS )

!  Jacobian control

   LOGICAL, Intent(in) :: Lvaryflags(maxlayers)
   INTEGER, Intent(in) :: Lvarynums (maxlayers)
   INTEGER, Intent(in) :: n_surfacewfs

!  optical inputs
!  --------------

!  Atmosphere extinction and deltaus

   real(fpk), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(fpk), Intent(in) :: OMEGA       ( MAXLAYERS )
   real(fpk), Intent(in) :: TRUNCFAC    ( MAXLAYERS )

!  Linearized optical inputs

   real(fpk), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
   real(fpk), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
   real(fpk), Intent(in) :: L_OMEGA       ( MAXLAYERS, max_atmoswfs )
   real(fpk), Intent(in) :: L_TRUNCFAC    ( MAXLAYERS, max_atmoswfs )

!  Atmospheric BB functions and Surface BB and emissivity

   REAL(fpk), Intent(in) :: SURFBB, USER_EMISSIVITY(MAXGEOMS)
   REAL(fpk), Intent(in) :: BB_INPUT (0:MAXLAYERS)
   REAL(fpk), Intent(in) :: LS_USER_EMISSIVITY  (maxgeoms, max_surfacewfs)

!  Geometrical inputs
!  ------------------

!  Ray constant, Cotangents, Critical layer
!    Mu1 = cos(alpha_boa), required for the Regular PS only

   integer  , Intent(in)  :: NCrit(maxgeoms)
   real(fpk), Intent(in)  :: Raycon(maxgeoms), cota(0:maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: Mu1(maxgeoms)

!  LOS Quadratures for Enhanced PS

   real(fpk), Intent(in)  :: xfine   (maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: wfine   (maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: csqfine (maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: cotfine (maxlayers,maxfinelayers,maxgeoms)

!  outputs
!  -------

!  Radiances

   real(fpk), Intent(Out)  :: intensity_dta_up ( max_user_levels, maxgeoms )
   real(fpk), Intent(Out)  :: intensity_dts    ( max_user_levels, maxgeoms )

   real(fpk), Intent(Out)  :: LP_Jacobians_dta_up  ( max_user_levels, maxgeoms, maxlayers, max_atmoswfs )
   real(fpk), Intent(Out)  :: LP_Jacobians_dts_up  ( max_user_levels, maxgeoms, maxlayers, max_atmoswfs )
   real(fpk), Intent(Out)  :: LS_Jacobians_dts     ( max_user_levels, maxgeoms, max_surfacewfs )

!  Thermal setup

   real(fpk), Intent(InOut)   :: tcom1(maxlayers,2)
   real(fpk), Intent(InOut)   :: L_tcom1(maxlayers,2,max_atmoswfs)

!  LOCAL
!  -----

!  Local solutions (enhanced_ps case)

   real(fpk)  :: Solutions_fine (maxlayers,maxfinelayers)
   real(fpk)  :: Wtrans_fine    (maxlayers,maxfinelayers)
   real(fpk)  :: L_Solutions_fine (maxlayers,maxfinelayers,max_atmoswfs)
   real(fpk)  :: L_Wtrans_fine    (maxlayers,maxfinelayers,max_atmoswfs)

!  Source function integration results

   real(fpk)  :: sources_up       ( maxlayers )
   real(fpk)  :: lostrans_up      ( maxlayers )
   real(fpk)  :: cumsource_up     ( 0:maxlayers )

   real(fpk)  :: L_lostrans_up    ( maxlayers,max_atmoswfs )
   real(fpk)  :: L_cumsource      ( max_atmoswfs )
   real(fpk)  :: L_sources_up     ( maxlayers,max_atmoswfs )
   real(fpk)  :: LS_cumsource     ( max_surfacewfs )

!  Regular_PS or plane-parallel flag

   logical    :: do_RegPSorPP

!  Help

   integer    :: n, j, k, nj, q, v, uta, nstart, nc, nut, nut_prev, Qnums(maxlayers)
   logical    :: layermask_up(maxlayers), Qvary(maxlayers)

   real(fpk)  :: argum(maxfinelayers), tran(maxfinelayers)
   real(fpk)  :: help, sum, kn, ke, xjkn, cot_1, cot_2, arg, tcw
   real(fpk)  :: L_help, L_sum, L_tran, L_kn
   real(fpk)  :: cumsource_dste, t_mult_up(0:2), L_t_mult_up(0:2), thermcoeffs(2)
   real(fpk)  :: tms, L_tms(max_atmoswfs), lostau, LA2

   real(fpk), parameter  :: cutoff = 88.0d0
   real(fpk), parameter  :: zero   = 0.0_fpk
   real(fpk), parameter  :: one    = 1.0_fpk

!  Zero the output

   INTENSITY_dta_up = zero
   INTENSITY_dts    = zero

   LP_JACOBIANS_dta_up = zero
   LP_JACOBIANS_dts_up = zero
   LS_JACOBIANS_dts    = zero

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

!  Thermal setup factors and linearizations
!     TMS, Initial set of thermal coefficients and TCOM1 variable

   if ( do_Thermset ) then
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
      do_Thermset = .false.
   endif

!  Start Geometry loop
!  ===================

   do v = 1, ngeoms

!  Zero local sources

      lostrans_up   = zero  ; sources_up   = zero ; cumsource_up = zero
      L_lostrans_up = zero  ; L_sources_up = zero

!  Bookkeeping

      nstart = nlayers ; if (Ncrit(v).ne.0) nstart = nCrit(v)
      if (nstart.lt.nlayers) LAYERMASK_UP(nstart+1:nlayers) = .false.

!  Plane/Parallel and Regular-PS: Layer integrated source terms
!  ============================================================

!  Bug Fixed 23 January 2013 (nadir case). Old code commented out.

      if ( do_RegPSorPP ) then
        if ( doNadir(v) ) then
          DO n = 1, nlayers
            lostau = deltaus(n)
            if ( lostau .lt. cutoff ) lostrans_up(n) = exp( - lostau )
            if ( Qvary(n) ) L_lostrans_up(n,1:Qnums(n)) = - lostrans_up(n) * L_deltaus(n,1:Qnums(n))
            if ( layermask_up(n) ) then
              t_mult_up(2) = tcom1(n,2)
              t_mult_up(1) = tcom1(n,1) + t_mult_up(2) 
              sum = t_mult_up(1) + t_mult_up(2) * deltaus(n)
              t_mult_up(0) = - sum
              sources_up(n) = t_mult_up(0) * lostrans_up(n)  + t_mult_up(1)
!              t_mult_up(1:2) = tcom1(n,1:2)
!              sum = t_mult_up(1) + t_mult_up(2) * deltaus(n)
!              t_mult_up(0) = - sum
!              sources_up(n) = t_mult_up(1)           ! Check this formula
              if ( Qvary(n) ) then
                do q = 1, Qnums(n)
                   L_t_mult_up(2) = L_tcom1(n,2,q)
                   L_t_mult_up(1) = L_tcom1(n,1,q) + L_t_mult_up(2) 
                   L_sum = L_t_mult_up(1) + t_mult_up(2) * L_deltaus(n,q) + L_t_mult_up(2) * deltaus(n)
                   L_t_mult_up(0) = - L_sum
                   L_sources_up(n,q) = L_t_mult_up(0) * lostrans_up(n) + t_mult_up(0) * L_lostrans_up(n,q) + L_t_mult_up(1) 
!                   L_t_mult_up(1:2) = L_tcom1(n,1:2,q)
!                   L_sum = L_t_mult_up(1) + t_mult_up(2) * L_deltaus(n,q) + L_t_mult_up(2) * deltaus(n)
!                   L_t_mult_up(0) = - sum
!                   L_sources_up(n,q) = L_t_mult_up(1)           ! Check this formula
                enddo
              endif
            endif
          enddo
        else
          DO n = 1, nlayers
            lostau = deltaus(n) / Mu1(v)
            if ( lostau .lt. cutoff ) lostrans_up(n) = exp( - lostau )
            if ( Qvary(n) ) L_lostrans_up(n,1:Qnums(n)) = - lostrans_up(n) * L_deltaus(n,1:Qnums(n)) / Mu1(v)
            if ( layermask_up(n) ) then
              t_mult_up(2) = tcom1(n,2)
              t_mult_up(1) = tcom1(n,1) + t_mult_up(2) * Mu1(v)
              sum = t_mult_up(1) + t_mult_up(2) * deltaus(n)
              t_mult_up(0) = - sum
              sources_up(n) = t_mult_up(0) * lostrans_up(n)  + t_mult_up(1)
              if ( Qvary(n) ) then
                do q = 1, Qnums(n)
                   L_t_mult_up(2) = L_tcom1(n,2,q)
                   L_t_mult_up(1) = L_tcom1(n,1,q) + L_t_mult_up(2) * Mu1(v)
                   L_sum = L_t_mult_up(1) + t_mult_up(2) * L_deltaus(n,q) + L_t_mult_up(2) * deltaus(n)
                   L_t_mult_up(0) = - L_sum
                   L_sources_up(n,q) = L_t_mult_up(0) *   lostrans_up(n)   + &
                                         t_mult_up(0) * L_lostrans_up(n,q) + L_t_mult_up(1)
                enddo
              endif
            endif
          enddo
        endif
      endif

!  LOS-spherical Layer integrated source terms
!  ===========================================

      if ( do_enhanced_ps ) then
         do n = nlayers, 1, -1
            kn = extinction(n) ; nj = nfinedivs(n,v)
            if (  doNadir(v)  .and. layermask_up(n) ) then
               lostau = deltaus(n)
               if ( lostau .lt. cutoff ) lostrans_up(n) = exp( - lostau )
               if ( Qvary(n) ) L_lostrans_up(n,1:Qnums(n)) = -lostrans_up(n) * L_deltaus(n,1:Qnums(n))
               do j = 1, nj
                  argum(j) = xfine(n,j,v) ; xjkn = argum(j) * kn
                  tran(j)  = exp ( -xjkn )
                  solutions_fine(n,j) = tcom1(n,1) + xjkn * tcom1(n,2)
                  wtrans_fine(n,j)    = kn * tran(j) * wfine(n,j,v)
                  if ( Qvary(n) ) then
                     do q = 1, Qnums(n)
                        L_kn = L_extinction(n,q) ; L_tran = argum(j) * L_kn
                        L_solutions_fine(n,j,q) = L_tcom1(n,1,q) + xjkn * L_tcom1(n,2,q) + L_tran * tcom1(n,2)
                        L_wtrans_fine(n,j,q)    = ( L_kn - kn * L_tran ) * tran(j) * wfine(n,j,v)
                     enddo
                  endif
               enddo
            else if ( .not. doNadir(v) .and. layermask_up(n) ) then
               cot_2 = cota(n-1,v) ; cot_1 = cota(n,v)
               ke = raycon(v) * kn  ; arg = raycon(v) * ( cot_2 - cot_1 ) ; lostau = kn * arg
               if ( lostau .lt. cutoff ) lostrans_up(n) = exp( - lostau )
               if ( Qvary(n) ) L_lostrans_up(n,1:Qnums(n)) = - arg * lostrans_up(n) * L_extinction(n,1:Qnums(n))
               do j = 1, nj
                  argum(j) = Raycon(v) * ( cot_2 - cotfine(n,j,v) )  ; xjkn = xfine(n,j,v) * kn
                  tran(j)  = exp ( - kn * argum(j) ) ; tcw = tran(j) * csqfine(n,j,v) * wfine(n,j,v)
                  solutions_fine(n,j) = tcom1(n,1) + xjkn * tcom1(n,2)
                  wtrans_fine(n,j)    = ke * tcw
                  if ( Qvary(n) ) then
                     do q = 1, Qnums(n)
                        L_kn = L_extinction(n,q) ; L_tran = argum(j) * L_kn
                        L_solutions_fine(n,j,q) = L_tcom1(n,1,q) + xfine(n,j,v)*(kn*L_tcom1(n,2,q) + L_kn*tcom1(n,2))
                        L_wtrans_fine(n,j,q)    = raycon(v) * ( L_kn - kn * L_tran ) * tcw
                     enddo
                  endif
               enddo
            endif
            sources_up(n) = dot_product(solutions_fine(n,1:nj),wtrans_fine(n,1:nj))
            if ( Qvary(n) ) then
               do q = 1, Qnums(n)
                  L_sources_up(n,q) = dot_product(L_solutions_fine(n,1:nj,q),  wtrans_fine(n,1:nj)) + &
                                      dot_product(  solutions_fine(n,1:nj),  L_wtrans_fine(n,1:nj,q))
               enddo
            endif
         enddo
      endif

!  Source function integration
!  ===========================

!  start recursion ( For DSTE term, Use surface emissivity )

      NC =  0
      CUMSOURCE_UP(NC) = zero
      CUMSOURCE_DSTE   = SURFBB * USER_EMISSIVITY(V)
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  Main loop over all output optical depths
!     NLEVEL = Layer index for given optical depth
!     Cumulative source terms : Loop over layers working upwards from NSTART to level NUT,
!     Check for updating the recursion

      DO UTA = N_USER_LEVELS, 1, -1
         NUT    = USER_LEVELS(UTA) + 1
         DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            CUMSOURCE_DSTE   = LOSTRANS_UP(N) * CUMSOURCE_DSTE
            CUMSOURCE_UP(NC) = LOSTRANS_UP(N) * CUMSOURCE_UP(NC-1) + SOURCES_UP(N)
        ENDDO
         INTENSITY_DTA_UP(UTA,V) = CUMSOURCE_UP(NC)
         INTENSITY_DTS(UTA,V)    = CUMSOURCE_DSTE
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
         NUT_PREV = NUT
      ENDDO

!  Surface WFs

      if ( do_surfacewfs ) then
         do q = 1, n_surfacewfs
            LS_cumsource(q) = SURFBB * LS_USER_EMISSIVITY(V,q)
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
               LS_Jacobians_dts(UTA,V,Q) = LS_CUMSOURCE(Q)
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
                           L_cumsource(q) = L_SOURCES_UP(N,Q)           + &
                                L_LOSTRANS_UP(N,Q) * CUMSOURCE_UP(NC-1) + &
                                  LOSTRANS_UP(N)   * L_CUMSOURCE(Q)
                        enddo
                     else
                        do q = 1, Qnums(k)
                           L_cumsource(q) = LOSTRANS_UP(N) * L_CUMSOURCE(Q)
                        enddo
                     endif
                  ENDDO
                  do q = 1, Qnums(k)
                     LP_Jacobians_dta_up(UTA,V,K,Q) = L_CUMSOURCE(Q)
                  enddo
                  IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
                  NUT_PREV = NUT
               ENDDO
            endif
         enddo
      endif

!  Profile Wfs (surface emission term)

      if ( do_profilewfs ) then
         CUMSOURCE_DSTE   = SURFBB * USER_EMISSIVITY(v)
         do k = 1, nlayers
            if ( Qvary(k) ) then
               L_CUMSOURCE = CUMSOURCE_DSTE
               NSTART = NLAYERS
               NUT_PREV = NSTART + 1
               DO UTA = N_USER_LEVELS, 1, -1
                  NUT    = USER_LEVELS(UTA) + 1
                  DO N = NSTART, NUT, -1
                     NC = NLAYERS + 1 - N
                     if ( k.eq.n ) then
                        do q = 1, Qnums(k)
                           L_cumsource(q) =  L_LOSTRANS_UP(N,Q) * L_CUMSOURCE(Q)
                        enddo
                     else
                        do q = 1, Qnums(k)
                           L_cumsource(q) = LOSTRANS_UP(N) * L_CUMSOURCE(Q)
                        enddo
                     endif
                  ENDDO
                  do q = 1, Qnums(k)
                     LP_Jacobians_dts_up(UTA,V,K,Q) = L_CUMSOURCE(Q)
                  enddo
                  IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
                  NUT_PREV = NUT
               ENDDO
            endif
         enddo
      endif

!  End geometry loop

   enddo

!  Finish

   return
end subroutine DTE_Integral_ILPS_UP

!

subroutine DTE_Integral_ILPS_DN &
   ( maxgeoms, maxlayers, maxfinelayers, max_user_levels,         & ! Inputs (dimensioning)
     max_atmoswfs,                                                & ! Inputs (dimensioning)
     Do_Thermset, do_deltam_scaling, do_PlanPar,                  & ! Inputs (Flags)
     do_regular_ps, do_enhanced_ps, doNadir,                      & ! Inputs (Flags)
     do_profilewfs, Lvaryflags, Lvarynums,                        & ! Inputs (Control, Jacobians)
     ngeoms, nlayers, nfinedivs, n_user_levels, user_levels,      & ! Inputs (control, output)
     bb_input, extinction, deltaus, omega, truncfac,              & ! Inputs (Optical - Regular)
     L_extinction, L_deltaus, L_omega, L_truncfac,                & ! Inputs (Optical - Linearized)
     Mu1, NCrit, Radcrit, CotCrit, Raycon, cota,                  & ! Inputs (Geometry)
     radii, xfine, wfine, csqfine, cotfine,                       & ! Inputs (Geometry)
     intensity_dta_dn, LP_Jacobians_dta_dn, tcom1, L_tcom1 )        ! Output

!  Stand alone routine for Downwelling Direct-thermal-emission (DTE)
!    computation of Radiances and LPS Jacobians. Inputs: geometry, optical properties, Planck functions, emissivity

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
   integer, Intent(in) :: max_user_levels
   INTEGER, Intent(in) :: max_atmoswfs

!  Thermal setup flag (for TCOM1)

   LOGICAL, Intent(inout) :: Do_Thermset

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
   integer, Intent(in) :: NGEOMS

   integer, Intent(in) :: N_USER_LEVELS
   integer, Intent(in) :: USER_LEVELS ( MAX_USER_LEVELS )

!  Jacobian control

   LOGICAL, Intent(in) :: Lvaryflags(maxlayers)
   INTEGER, Intent(in) :: Lvarynums (maxlayers)

!  optical inputs
!  --------------

!  Atmosphere extinction and deltaus

   real(fpk), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(fpk), Intent(in) :: OMEGA       ( MAXLAYERS )
   real(fpk), Intent(in) :: TRUNCFAC    ( MAXLAYERS )

!  Linearized optical inputs

   real(fpk), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
   real(fpk), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
   real(fpk), Intent(in) :: L_OMEGA       ( MAXLAYERS, max_atmoswfs )
   real(fpk), Intent(in) :: L_TRUNCFAC    ( MAXLAYERS, max_atmoswfs )

!  Atmospheric BB functions

   REAL(fpk), Intent(in) :: BB_INPUT (0:MAXLAYERS)

!  Geometrical inputs
!  ------------------

!  Ray constant, Cotangents, Critical layer
!    Mu1 = cos(alpha_boa), required for the Regular PS only

   integer  , Intent(in)  :: NCrit(maxgeoms)
   real(fpk), Intent(in)  :: RadCrit(maxgeoms), CotCrit(maxgeoms)
   real(fpk), Intent(in)  :: Raycon(maxgeoms), cota(0:maxlayers,maxgeoms), radii(0:maxlayers)
   real(fpk), Intent(in)  :: Mu1(maxgeoms)

!  LOS Quadratures for Enhanced PS

   real(fpk), Intent(in)  :: xfine   (maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: wfine   (maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: csqfine (maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: cotfine (maxlayers,maxfinelayers,maxgeoms)

!  outputs
!  -------

!  Radiances

   real(fpk), Intent(Out)  :: intensity_dta_dn ( max_user_levels, maxgeoms )

!  Jacobians

   real(fpk), Intent(Out)  :: LP_Jacobians_dta_dn  ( max_user_levels, maxgeoms, maxlayers, max_atmoswfs )

!  Thermal setup

   real(fpk), Intent(InOut)   :: tcom1(maxlayers,2)
   real(fpk), Intent(InOut)   :: L_tcom1(maxlayers,2,max_atmoswfs)

!  LOCAL
!  -----

!  Local solutions (enhanced_ps case)

   real(fpk)  :: Solutions_fine (maxlayers,maxfinelayers)
   real(fpk)  :: Wtrans_fine    (maxlayers,maxfinelayers)
   real(fpk)  :: L_Solutions_fine (maxlayers,maxfinelayers,max_atmoswfs)
   real(fpk)  :: L_Wtrans_fine    (maxlayers,maxfinelayers,max_atmoswfs)

!  Source function integration results

   real(fpk)  :: sources_dn       ( maxlayers )
   real(fpk)  :: lostrans_dn      ( maxlayers )
   real(fpk)  :: cumsource_dn     ( 0:maxlayers )

   real(fpk)  :: L_lostrans_dn    ( maxlayers,max_atmoswfs )
   real(fpk)  :: L_cumsource      ( max_atmoswfs )
   real(fpk)  :: L_sources_dn     ( maxlayers,max_atmoswfs )

!  Regular_PS or plane-parallel flag

   logical    :: do_RegPSorPP

!  Help

   integer    :: n, j, k, nj, q, v, uta, nstart, nc, nut, nut_prev, Qnums(maxlayers)
   logical    :: layermask_dn(maxlayers), Qvary(maxlayers)

   real(fpk)  :: argum(maxfinelayers), tran(maxfinelayers)
   real(fpk)  :: help, sum, kn, ke, xjkn, cot_1, cot_2, cot_c, rdiff, arg0, arg, trand, tcw
   real(fpk)  :: L_help, L_sum, L_tran, L_kn, L_trand
   real(fpk)  :: t_mult_dn(0:2), L_t_mult_dn(0:2), thermcoeffs(2)
   real(fpk)  :: tms, L_tms(max_atmoswfs), lostau, LA2

   real(fpk), parameter  :: cutoff = 88.0d0
   real(fpk), parameter  :: zero   = 0.0_fpk
   real(fpk), parameter  :: one    = 1.0_fpk

!  Zero the output

   INTENSITY_dta_dn    = zero
   LP_JACOBIANS_dta_dn = zero

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

!  Thermal setup factors and linearizations
!     TMS, Initial set of thermal coefficients and TCOM1 variable

   if ( do_Thermset ) then
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
      do_Thermset = .false.
   endif

!  Start Geometry loop
!  ===================

   do v = 1, ngeoms

!  Zero local sources

      lostrans_dn   = zero  ; sources_dn   = zero ; cumsource_dn = zero
      L_lostrans_dn = zero  ; L_sources_dn = zero
      trand = one ; L_trand = zero

!  Criticality

      nstart = nlayers ; if (Ncrit(v).ne.0) nstart = nCrit(v)
      if (nstart.lt.nlayers) LAYERMASK_DN(nstart+1:nlayers) = .false.

!  Plane/Parallel and Regular-PS: Layer integrated source terms
!  ============================================================

!  Bug Fixed 23 January 2013 (nadir case). Old code commented out and replaced

      if ( do_RegPSorPP ) then
        if ( doNadir(v) ) then
          DO n = 1, nlayers
            lostau = deltaus(n)
            if ( lostau .lt. cutoff ) lostrans_dn(n) = exp( - lostau )
            if ( Qvary(n) ) L_lostrans_dn(n,1:Qnums(n)) = -lostrans_dn(n) * L_deltaus(n,1:Qnums(n))
            if ( layermask_dn(n) ) then
              t_mult_dn(2)   = tcom1(n,2)
              t_mult_dn(1)   = tcom1(n,1) - t_mult_dn(2)
              t_mult_dn(0)   = - t_mult_dn(1)
              sources_dn(n)  = t_mult_dn(0) * lostrans_dn(n)
              sum = t_mult_dn(1) + t_mult_dn(2) * deltaus(n)
              sources_dn(n)  = sources_dn(n) + sum
!              t_mult_dn(1:2) = tcom1(n,1:2)
!              t_mult_dn(0)   = - t_mult_dn(1)
!              sum = t_mult_dn(1) + t_mult_dn(2) * deltaus(n)
!              sources_dn(n)  = sum
              if ( Qvary(n) ) then
                do q = 1, Qnums(n)
                   L_t_mult_dn(2) = L_tcom1(n,2,q)
                   L_t_mult_dn(1) = L_tcom1(n,1,q) - L_t_mult_dn(2)
                   L_t_mult_dn(0)   = - L_t_mult_dn(1)
                   L_sources_dn(n,q)  = L_t_mult_dn(0) * lostrans_dn(n) + t_mult_dn(0) * L_lostrans_dn(n,q)
                   L_sum = L_t_mult_dn(1) + L_t_mult_dn(2) * deltaus(n) + t_mult_dn(2) * L_deltaus(n,q)
                   L_sources_dn(n,q) = L_sources_dn(n,q) + L_sum
!                   L_t_mult_dn(1:2) = L_tcom1(n,1:2,q)
!                   L_t_mult_dn(0)   = - L_t_mult_dn(1)
!                   L_sum = L_t_mult_dn(1) + t_mult_dn(2) * L_deltaus(n,q) + L_t_mult_dn(2) * deltaus(n)
!                   L_sources_dn(n,q) = L_sum          ! Check this formula
                enddo
              endif
            endif
          enddo
        else
          DO n = 1, nlayers
            lostau = deltaus(n) / Mu1(v)
            if ( lostau .lt. cutoff ) lostrans_dn(n) = exp( - lostau )
            if ( Qvary(n) ) L_lostrans_dn(n,1:Qnums(n)) = - lostrans_dn(n) * L_deltaus(n,1:Qnums(n)) / Mu1(v)
            if ( layermask_dn(n) ) then
              t_mult_dn(2)   = tcom1(n,2)
              t_mult_dn(1)   = tcom1(n,1) - t_mult_dn(2) * Mu1(v)
              t_mult_dn(0)   = - t_mult_dn(1)
              sources_dn(n)  = t_mult_dn(0) * lostrans_dn(n)
              sum = t_mult_dn(1) + t_mult_dn(2) * deltaus(n)
              sources_dn(n)  = sources_dn(n) + sum
              if ( Qvary(n) ) then
                do q = 1, Qnums(n)
                   L_t_mult_dn(2) = L_tcom1(n,2,q)
                   L_t_mult_dn(1) = L_tcom1(n,1,q) - L_t_mult_dn(2) * Mu1(v)
                   L_t_mult_dn(0) = - L_t_mult_dn(1)
                   L_sources_dn(n,q)  = L_t_mult_dn(0) * lostrans_dn(n) + t_mult_dn(0) * L_lostrans_dn(n,q)
                   L_sum = L_t_mult_dn(1) + t_mult_dn(2) * L_deltaus(n,q) + L_t_mult_dn(2) * deltaus(n)
                   L_sources_dn(n,q) = L_sources_dn(n,q) + L_sum
                enddo
              endif
            endif
          enddo
        endif
      endif

!  LOS-spherical Layer integrated source terms
!  ===========================================

!  Bug Fixed 23 January 2013 (nadir case). Old code commented out and replaced

      if ( do_enhanced_ps ) then
         do n = nlayers, 1, -1
            kn = extinction(n) ; nj = nfinedivs(n,v)
            if (  doNadir(v)  .and. layermask_dn(n) ) then
               rdiff = radii(n-1) - radii(n) ; if ( n.eq.NCrit(v)) rdiff = radii(n-1) - RadCrit(v)
               trand = one ; if ( n.eq.NCrit(v)) trand = exp ( -kn * (RadCrit(v) -radii(n) ) )
               lostau = deltaus(n)
               if ( lostau .lt. cutoff ) lostrans_dn(n) = exp( - lostau )
               if ( Qvary(n) ) L_lostrans_dn(n,1:Qnums(n)) = -lostrans_dn(n) * L_deltaus(n,1:Qnums(n))
               do j = 1, nj
                  argum(j) = rdiff - xfine(n,j,v) ; xjkn = xfine(n,j,v) * kn
                  tran(j)  = exp ( - kn * argum(j) )
!                  tran(j)  = exp ( -xjkn )
                  solutions_fine(n,j) = tcom1(n,1) + xjkn * tcom1(n,2)
                  wtrans_fine(n,j)    = kn * tran(j) * wfine(n,j,v)
                  if ( Qvary(n) ) then
                     do q = 1, Qnums(n)
                        L_kn = L_extinction(n,q) ; L_tran = argum(j) * L_kn
                        L_solutions_fine(n,j,q) = L_tcom1(n,1,q) + xjkn * L_tcom1(n,2,q) + L_tran * tcom1(n,2)
                        L_wtrans_fine(n,j,q)    = ( L_kn - kn * L_tran ) * tran(j) * wfine(n,j,v)
                     enddo
                  endif
               enddo
            else if ( .not. doNadir(v) .and. layermask_dn(n) ) then
               cot_2 = cota(n-1,v) ; cot_1 = cota(n,v)
               cot_c = cot_1     ; if ( n.eq.NCrit(v) ) cot_c = CotCrit(v)
               ke = raycon(v) * kn  ; arg = raycon(v) * ( cot_2 - cot_1 ) ; lostau = kn * arg
               if ( n.eq.NCrit(v) ) then
                  arg0 = Raycon(v) * ( CotCrit(v) - cot_1 ) ; trand = exp ( - kn * arg0 )
               endif
               if ( lostau .lt. cutoff ) lostrans_dn(n) = exp( - lostau )
               if ( Qvary(n) ) L_lostrans_dn(n,1:Qnums(n)) = - arg * lostrans_dn(n) * L_extinction(n,1:Qnums(n))
               do j = 1, nj
                  argum(j) = raycon(v) * ( cotfine(n,j,v) - cot_1 )  ; xjkn = xfine(n,j,v) * kn
                  tran(j)  = exp ( - kn * argum(j) ) ; tcw = tran(j) * csqfine(n,j,v) * wfine(n,j,v)
                  solutions_fine(n,j) = tcom1(n,1) + xjkn * tcom1(n,2)
                  wtrans_fine(n,j)    = ke * tcw
                  if ( Qvary(n) ) then
                     do q = 1, Qnums(n)
                        L_kn = L_extinction(n,q) ; L_tran = argum(j) * L_kn
                        L_solutions_fine(n,j,q) = L_tcom1(n,1,q) + xfine(n,j,v)*(kn*L_tcom1(n,2,q) + L_kn*tcom1(n,2))
                        L_wtrans_fine(n,j,q)    = Raycon(v) * ( L_kn - kn * L_tran ) * tcw
                     enddo
                  endif
               enddo
            endif
            sources_dn(n) = dot_product(solutions_fine(n,1:nj),wtrans_fine(n,1:nj))
            if ( Qvary(n) ) then
               do q = 1, Qnums(n)
                  L_sources_dn(n,q) = dot_product(L_solutions_fine(n,1:nj,q),  wtrans_fine(n,1:nj)) + &
                                      dot_product(  solutions_fine(n,1:nj),  L_wtrans_fine(n,1:nj,q) )
                  if ( n.eq.NCrit(v) ) then
                     L_trand =  - L_extinction(n,q) * trand * arg0
                     L_sources_dn(n,q) = trand  * L_sources_dn(n,q) + L_trand * sources_dn(n)    !@@ Robfix
                  endif
               enddo
            endif
            if ( n.eq.NCrit(v) ) sources_dn(n) = sources_dn(n) * trand         !@@ Robfix
         enddo
      endif

!  Source function integration
!  ===========================

!  start recursion ( For DSTE term, Use surface emissivity )

      NC =  0
      CUMSOURCE_DN(NC) = zero
      NSTART = 1
      NUT_PREV = NSTART - 1

!  Main loop over all output optical depths
!     NLEVEL = Layer index for given optical depth
!     Cumulative source terms : Loop over layers working upwards from NSTART to level NUT,
!     Check for updating the recursion

      DO UTA = 1, N_USER_LEVELS
         NUT    = USER_LEVELS(UTA)
         DO N = NSTART, NUT
            NC = N
            CUMSOURCE_DN(NC) = SOURCES_DN(N) + LOSTRANS_DN(N) * CUMSOURCE_DN(NC-1)
         ENDDO
         INTENSITY_DTA_DN(UTA,V) = CUMSOURCE_DN(NC)
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
                           L_cumsource(q) = L_SOURCES_DN(N,Q)           + &
                                L_LOSTRANS_DN(N,Q) * CUMSOURCE_DN(NC-1) + &
                                  LOSTRANS_DN(N)   * L_CUMSOURCE(Q)
                        enddo
                     else
                        do q = 1, Qnums(k)
                           L_cumsource(q) = LOSTRANS_DN(N) * L_CUMSOURCE(Q)
                        enddo
                     endif
                  ENDDO
                  do q = 1, Qnums(k)
                     LP_Jacobians_dta_DN(UTA,V,K,Q) = L_CUMSOURCE(Q)
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
end subroutine DTE_Integral_ILPS_DN

!

subroutine DTE_Integral_ILPS_UPDN &
   ( maxgeoms, maxlayers, maxfinelayers, max_user_levels,               & ! Inputs (dimensioning)
     max_atmoswfs, max_surfacewfs,                                      & ! Inputs (dimensioning)
     do_upwelling, do_dnwelling,                                        & ! Inputs (Flags)
     Do_Thermset, do_deltam_scaling, do_PlanPar,                        & ! Inputs (Flags)
     do_regular_ps, do_enhanced_ps, doNadir,                            & ! Inputs (Flags)
     do_profilewfs, do_surfacewfs, Lvaryflags, Lvarynums, n_surfacewfs, & ! Inputs (Control, Jacobians)
     ngeoms, nlayers, nfinedivs, n_user_levels, user_levels,            & ! Inputs (control, output)
     bb_input, surfbb, user_emissivity, LS_user_emissivity,             & ! Inputs (Thermal)
     extinction, deltaus, omega, truncfac,                              & ! Inputs (Optical - Regular)
     L_extinction, L_deltaus, L_omega, L_truncfac,                      & ! Inputs (Optical - Linearized)
     Mu1, NCrit, Radcrit, CotCrit, Raycon, cota,                        & ! Inputs (Geometry)
     radii, xfine, wfine, csqfine, cotfine,                             & ! Inputs (Geometry)
     intensity_dta_up, intensity_dts, intensity_dta_dn,                 & ! Outputs
     LP_Jacobians_dta_up, LP_Jacobians_dts_up, LS_Jacobians_dts,        & ! Output
     LP_Jacobians_dta_dn, tcom1, L_tcom1 )                                ! Output

!  Stand alone routine for Upwelling and downwelling Direct-thermal-emission (DTE)
!    computation of Radiances and LPS Jacobians. Inputs: geometry, optical properties, Planck functions, emissivity

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
   integer, Intent(in) :: max_user_levels
   INTEGER, Intent(in) :: max_atmoswfs
   INTEGER, Intent(in) :: max_surfacewfs

!  General flags

   LOGICAL, Intent(in) :: DO_UPWELLING
   LOGICAL, Intent(in) :: DO_DNWELLING

!  Thermal setup flag (for TCOM1)

   LOGICAL, Intent(inout) :: Do_Thermset

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
   integer, Intent(in) :: NGEOMS

   integer, Intent(in) :: N_USER_LEVELS
   integer, Intent(in) :: USER_LEVELS ( MAX_USER_LEVELS )

!  Jacobian control

   LOGICAL, Intent(in) :: Lvaryflags(maxlayers)
   INTEGER, Intent(in) :: Lvarynums (maxlayers)
   INTEGER, Intent(in) :: n_surfacewfs

!  optical inputs
!  --------------

!  Atmosphere extinction and deltaus

   real(fpk), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(fpk), Intent(in) :: OMEGA       ( MAXLAYERS )
   real(fpk), Intent(in) :: TRUNCFAC    ( MAXLAYERS )

!  Linearized optical inputs

   real(fpk), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
   real(fpk), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
   real(fpk), Intent(in) :: L_OMEGA       ( MAXLAYERS, max_atmoswfs )
   real(fpk), Intent(in) :: L_TRUNCFAC    ( MAXLAYERS, max_atmoswfs )

!  Atmospheric BB functions and Surface BB and emissivity

   REAL(fpk), Intent(in) :: SURFBB, USER_EMISSIVITY(MAXGEOMS)
   REAL(fpk), Intent(in) :: BB_INPUT (0:MAXLAYERS)
   REAL(fpk), Intent(in) :: LS_USER_EMISSIVITY  (maxgeoms, max_surfacewfs)

!  Geometrical inputs
!  ------------------

!  Ray constant, Cotangents, Critical layer
!    Mu1 = cos(alpha_boa), required for the Regular PS only

   integer  , Intent(in)  :: NCrit(maxgeoms)
   real(fpk), Intent(in)  :: RadCrit(maxgeoms), CotCrit(maxgeoms)
   real(fpk), Intent(in)  :: Raycon(maxgeoms), cota(0:maxlayers,maxgeoms), radii(0:maxlayers)
   real(fpk), Intent(in)  :: Mu1(maxgeoms)

!  LOS Quadratures for Enhanced PS

   real(fpk), Intent(in)  :: xfine   (maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: wfine   (maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: csqfine (maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: cotfine (maxlayers,maxfinelayers,maxgeoms)

!  outputs
!  -------

!  Upwelling

   real(fpk), Intent(Out)  :: intensity_dta_up     ( max_user_levels, maxgeoms )
   real(fpk), Intent(Out)  :: intensity_dts        ( max_user_levels, maxgeoms )
   real(fpk), Intent(Out)  :: LP_Jacobians_dta_up  ( max_user_levels, maxgeoms, maxlayers, max_atmoswfs )
   real(fpk), Intent(Out)  :: LP_Jacobians_dts_up  ( max_user_levels, maxgeoms, maxlayers, max_atmoswfs )
   real(fpk), Intent(Out)  :: LS_Jacobians_dts     ( max_user_levels, maxgeoms, max_surfacewfs )

!  Downwelling

   real(fpk), Intent(Out)  :: intensity_dta_dn     ( max_user_levels, maxgeoms )
   real(fpk), Intent(Out)  :: LP_Jacobians_dta_dn  ( max_user_levels, maxgeoms, max_atmoswfs )

!  Thermal setup

   real(fpk), Intent(InOut)   :: tcom1(maxlayers,2)
   real(fpk), Intent(InOut)   :: L_tcom1(maxlayers,2,max_atmoswfs)


!  Upwelling

   if ( do_upwelling ) then
      call DTE_Integral_ILPS_UP &
   ( maxgeoms, maxlayers, maxfinelayers, max_user_levels,               & ! Inputs (dimensioning)
     max_atmoswfs, max_surfacewfs,                                      & ! Inputs (dimensioning)
     Do_Thermset, do_deltam_scaling, do_PlanPar,                        & ! Inputs (Flags)
     do_regular_ps, do_enhanced_ps, doNadir,                            & ! Inputs (Flags)
     do_profilewfs, do_surfacewfs, Lvaryflags, Lvarynums, n_surfacewfs, & ! Inputs (Control, Jacobians)
     ngeoms, nlayers, nfinedivs, n_user_levels, user_levels,            & ! Inputs (control, output)
     bb_input, surfbb, user_emissivity, LS_user_emissivity,             & ! Inputs (Thermal)
     extinction, deltaus, omega, truncfac,                              & ! Inputs (Optical - Regular)
     L_extinction, L_deltaus, L_omega, L_truncfac,                      & ! Inputs (Optical - Linearized)
     Mu1, NCrit, Raycon, cota, xfine, wfine, csqfine, cotfine,          & ! Inputs (Geometry)
     intensity_dta_up, intensity_dts, LP_Jacobians_dta_up,              & ! Outputs
     LP_Jacobians_dts_up, LS_Jacobians_dts, tcom1, L_tcom1  )             ! Output
   endif

!  Downwelling

   if ( do_dnwelling ) then
      call DTE_Integral_ILPS_DN &
   ( maxgeoms, maxlayers, maxfinelayers, max_user_levels,         & ! Inputs (dimensioning)
     max_atmoswfs,                                                & ! Inputs (dimensioning)
     Do_Thermset, do_deltam_scaling, do_PlanPar,                  & ! Inputs (Flags)
     do_regular_ps, do_enhanced_ps, doNadir,                      & ! Inputs (Flags)
     do_profilewfs, Lvaryflags, Lvarynums,                        & ! Inputs (Control, Jacobians)
     ngeoms, nlayers, nfinedivs, n_user_levels, user_levels,      & ! Inputs (control, output)
     bb_input, extinction, deltaus, omega, truncfac,              & ! Inputs (Optical - Regular)
     L_extinction, L_deltaus, L_omega, L_truncfac,                & ! Inputs (Optical - Linearized)
     Mu1, NCrit, Radcrit, CotCrit, Raycon, cota,                  & ! Inputs (Geometry)
     radii, xfine, wfine, csqfine, cotfine,                       & ! Inputs (Geometry)
     intensity_dta_dn, LP_Jacobians_dta_dn, tcom1, L_tcom1 )        ! Output
   endif

!  Finish

   return
end subroutine DTE_Integral_ILPS_UPDN

!  End module

end module FO_Thermal_RTCalcs_ILPS_m


