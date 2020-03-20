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

module FO_ScalarSS_Masters

!  All subroutines public

public

contains

subroutine FO_MASTER &
       ( maxgeoms, maxlayers, maxfine, maxmoments_input,               & ! Inputs (dimensioning)
         max_user_levels,                                              & ! Inputs (dimensioning)
         do_solar_sources, do_thermal_emission, do_surface_emission,   & ! Input flags
         do_planpar, do_regular_ps, do_enhanced_ps, do_deltam_scaling, & ! Input flags
         do_upwelling, do_dnwelling,                                   & ! Input flags
         ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels,    & ! Input numbers
         user_levels, dtr, Pie, doCrit, Acrit, vsign, eradius,         & ! Input geometry
         heights, doNadir, alpha_boa, theta_boa, phi_boa, Mu0, Mu1,    & ! Input geometry
         reflec, surfbb, bb_input, emiss,                              & ! Inputs (Optical - Surface)
         extinction, deltaus, omega, truncfac, phasmoms, flux,         & ! Inputs (Optical - Regular)
         intensity, fail, message, trace )                                                     ! Output

   USE FO_geometry_DTonly_m
   USE FO_geometry_SSonly_m

   USE FO_ScalarSS_spherfuncs_m
   USE FO_ScalarSS_RTCalcs_I_m

   USE FO_Thermal_RTCalcs_I_m

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs for testing
!  ==================

!  Dimensions
!  ----------

   integer  :: maxgeoms, maxlayers, maxfine, maxmoments_input
   integer  :: max_user_levels

!  Configuration inputs
!  --------------------

!  Sources control, including thermal

   logical  :: DO_THERMAL_EMISSION
   logical  :: DO_SURFACE_EMISSION
   logical  :: DO_SOLAR_SOURCES

!  Flags (sphericity flags should be mutually exclusive)

   logical  :: do_regular_ps
   logical  :: do_enhanced_ps
   logical  :: do_planpar

!  deltam scaling flag

   logical  :: DO_DELTAM_SCALING

!  Directional Flags

   logical  :: DO_UPWELLING, DO_DNWELLING

!  General and Geoemtry inputs
!  ---------------------------

!  Layer and geometry control. Finelayer divisions may be changed

   integer  :: ngeoms, nlayers
   integer  :: nfinedivs(maxlayers,maxgeoms)

!  output levels

   integer  :: N_USER_LEVELS
   integer  :: USER_LEVELS ( MAX_USER_LEVELS )

!  Earth radius + heights

   real(ffp) :: eradius, heights (0:maxlayers)

!  dtr = degrees-to-Radians. Pie = 3.14159..., VSIGN = +1 (Up); -1(Down)

   real(ffp) :: dtr, Pie, vsign

!  Flag for the Nadir case. Intent(inout), input if DO_LOSpaths set)

   logical  :: doNadir(maxgeoms)

!  input angles (Degrees)
!    Mu0 = cos(theta_boa), required for surface term (both regular & enhanced)
!    Mu1 = cos(alpha_boa), required for the Regular PS only

   real(ffp) :: Mu0(maxgeoms), Mu1(maxgeoms)
   real(ffp) :: alpha_boa(maxgeoms)
   real(ffp) :: theta_boa(maxgeoms)
   real(ffp) :: phi_boa(maxgeoms)

!  Critical adjustment for cloud layers

   logical  :: doCrit
   real(ffp) :: Acrit

!  optical inputs
!  --------------

!  Atmosphere

   integer   :: NMOMENTS_INPUT
   real(ffp) :: EXTINCTION  ( MAXLAYERS )
   real(ffp) :: DELTAUS     ( MAXLAYERS )
   real(ffp) :: OMEGA       ( MAXLAYERS )
   real(ffp) :: TRUNCFAC    ( MAXLAYERS )
   real(ffp) :: PHASMOMS    ( MAXLAYERS,0:MAXMOMENTS_INPUT )

!  Solar Flux and Surface reflectivity (Could be the albedo)

   real(ffp) :: REFLEC(MAXGEOMS), FLUX

!  Thermal inputs

   real(ffp) :: SURFBB
   real(ffp) :: BB_INPUT ( 0:MAXLAYERS )

!  Emissivity

   real(ffp) :: EMISS ( maxgeoms )

!  Subroutine outputs
!  ==================

!  Output values

   real(ffp) :: intensity(maxgeoms)

!  Exception handling

   logical           :: fail
   character (len=*) :: message
   character (len=*) :: trace

!  Other variables
!  ===============

!  Geometry routine outputs
!  ------------------------

!  Alphas,  Cotangents, Radii, Ray constant. 

   real(ffp)  :: radii    (0:maxlayers)
   real(ffp)  :: Raycon   (maxgeoms)
   real(ffp)  :: alpha    (0:maxlayers,maxgeoms)
   real(ffp)  :: cota     (0:maxlayers,maxgeoms)

!  Critical layer

   integer    :: Ncrit(maxgeoms)
   real(ffp)  :: RadCrit(maxgeoms), CotCrit(maxgeoms)

!  solar paths 

   integer    :: ntraverse  (0:maxlayers,maxgeoms)
   real(ffp)  :: sunpaths   (0:maxlayers,maxlayers,maxgeoms)
   integer    :: ntraverse_fine(maxlayers,maxfine,maxgeoms)
   real(ffp)  :: sunpaths_fine (maxlayers,maxlayers,maxfine,maxgeoms)
   real(ffp)  :: Chapfacs   (maxlayers,maxlayers,maxgeoms)

!  Cosine scattering angle

   real(ffp)  :: cosscat(maxgeoms)

!  LOS Quadratures for Enhanced PS

   real(ffp)  :: xfine   (maxlayers,maxfine,maxgeoms)
   real(ffp)  :: wfine   (maxlayers,maxfine,maxgeoms)
   real(ffp)  :: csqfine (maxlayers,maxfine,maxgeoms)
   real(ffp)  :: cotfine (maxlayers,maxfine,maxgeoms)

!  Fine layering output

   real(ffp)  :: alphafine (maxlayers,maxfine,maxgeoms)
   real(ffp)  :: radiifine (maxlayers,maxfine,maxgeoms)

!  Spherfunc routine outputs
!  -------------------------

!  Help variables

   real(ffp) :: DF1(MAXMOMENTS_INPUT)
   real(ffp) :: DF2(MAXMOMENTS_INPUT)

!  Legendre Polynomials

   real(ffp)  :: LegPoly_up(0:maxmoments_input,maxgeoms)
   real(ffp)  :: LegPoly_dn(0:maxmoments_input,maxgeoms)

!  RT Calculation outputs
!  ----------------------

!  SS routines output

   real(ffp)  :: intensity_up     ( max_user_levels,maxgeoms )
   real(ffp)  :: intensity_dn     ( max_user_levels,maxgeoms )
   real(ffp)  :: intensity_db     ( max_user_levels,maxgeoms )

!  Thermal routines output

   real(ffp)  :: intensity_dta_up ( max_user_levels,maxgeoms )
   real(ffp)  :: intensity_dta_dn ( max_user_levels,maxgeoms )
   real(ffp)  :: intensity_dts    ( max_user_levels,maxgeoms )

!  Thermal setup and linearization

   real(ffp)  :: tcom1(maxlayers,2)

!  Dummies
   real(ffp)  :: cumsource_up     ( 0:maxlayers,maxgeoms )
   real(ffp)  :: cumsource_dn     ( 0:maxlayers,maxgeoms )

!  LOCAL
!  -----

!  numbers

   real(ffp), parameter :: zero = 0.0_ffp

!  help variables

   integer   :: n
   logical   :: STARTER, do_Thermset, do_lospaths, do_Chapman

!  Initialize

   intensity   = zero

   fail    = .false.
   message = ' '
   trace   = ' '

!  Flags to be set foreach calculation (safety)

   do_Chapman  = .false.
   do_Thermset = .true.
   do_lospaths = .false.
   starter     = .true.

!  Solar sources run (NO THERMAL)
!  ------------------------------

        if ( do_solar_sources.and..not. do_thermal_emission ) then

!  Geometry call

           call SS_Geometry_1 &
       ( maxgeoms, maxlayers, maxfine, ngeoms, nlayers, nfinedivs, dtr, Pie, & ! Input
         vsign, eradius, heights, alpha_boa, theta_boa, phi_boa, do_Chapman, & ! Input
         do_LOSpaths, do_planpar, do_regular_ps, do_enhanced_ps, doNadir,    & ! Input flags
         doCrit, Acrit, extinction, Raycon, radii, alpha, cota,              & ! Input/Output
         xfine, wfine, csqfine, cotfine, alphafine, radiifine,               & ! Input/Output
         NCrit, RadCrit, CotCrit, cosscat, chapfacs,                         & ! Output
         sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,                 & ! Output
         fail, message, trace )                                                ! Output

           !if ( fail ) then
           !   write(*,'(a)')message
           !   write(*,'(a)')trace
           !   stop'SS GEOMETRY failed'
           !endif
           if ( fail ) return

!  Spherical functions call

           if ( do_upwelling ) then
              Call FO_ScalarSS_spherfuncs ( STARTER, MAXMOMENTS_INPUT, MAXGEOMS, &
                 NMOMENTS_INPUT, NGEOMS, DF1, DF2, COSSCAT, LEGPOLY_UP )
           else if ( do_dnwelling ) then
              Call FO_ScalarSS_spherfuncs ( STARTER, MAXMOMENTS_INPUT, MAXGEOMS, &
                 NMOMENTS_INPUT, NGEOMS, DF1, DF2, COSSCAT, LEGPOLY_DN )
           endif

!  RT Call

           if ( do_upwelling ) then
              call SS_Integral_I_UP &
      ( maxgeoms, maxlayers, maxfine, maxmoments_input, max_user_levels,        & ! Inputs (dimensioning)
        do_deltam_scaling, do_PlanPar, do_regular_ps, do_enhanced_ps, doNadir,  & ! Inputs (Flags)
        ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels, & ! Inputs (control output)
        reflec, extinction, deltaus, omega, truncfac, phasmoms, flux,           & ! Inputs (Optical)
        Mu0, Mu1, LegPoly_up, NCrit, xfine, wfine, csqfine, cotfine,            & ! Inputs (Geometry)
        Raycon, cota, sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,       & ! Inputs (Geometry)
        intensity_up, intensity_db, cumsource_up )                                ! Outputs
           else if ( do_dnwelling ) then
              call SS_Integral_I_DN &
      ( maxgeoms, maxlayers, maxfine, maxmoments_input, max_user_levels,          & ! Inputs (dimensioning)
        do_deltam_scaling, do_PlanPar, do_regular_ps, do_enhanced_ps, doNadir,    & ! Inputs (Flags)
        ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,   & ! Inputs (control output)
        extinction, deltaus, omega, truncfac, phasmoms, flux,                     & ! Inputs (Optical)
        Mu1, LegPoly_dn, NCrit, RadCrit, CotCrit, xfine, wfine, csqfine, cotfine, & ! Inputs (Geometry)
        Raycon, radii, cota, sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,  & ! Inputs (Geometry)
        intensity_dn, cumsource_dn )                                                ! Outputs
           endif

!  save results (TOA or BOA)

           if ( do_upwelling ) then
              intensity(1:ngeoms) = intensity_up(1,1:ngeoms) + intensity_db(1,1:ngeoms)
           else
              intensity(1:ngeoms) = intensity_dn(n_user_levels,1:ngeoms)
           endif

!  End solar only

        endif

!  Thermal sources run (NO SOLAR)
!  ------------------------------

        if ( do_thermal_emission.and..not. do_solar_sources ) then

!  Geometry call

           call DT_Geometry &
       ( maxgeoms, maxlayers, maxfine, ngeoms, nlayers, nfinedivs,            & ! Inputs
         dtr, eradius, heights, alpha_boa,                                    & ! Input
         do_LOSpaths, do_planpar, do_regular_ps, do_enhanced_ps, doNadir,     & ! Input
         doCrit, Acrit, extinction, Raycon, radii, alpha, cota, xfine, wfine, & ! Input/Output
         csqfine, cotfine, alphafine, radiifine, NCrit, RadCrit, CotCrit,     & ! Output
         fail, message, trace )                                                 ! Output

           !if ( fail ) then
           !   write(*,'(a)')message
           !   write(*,'(a)')trace
           !   stop'DT GEOMETRY failed'
           !endif
           if ( fail ) return

!  RT Call

           if ( do_upwelling ) then
              call DTE_Integral_I_UP &
      ( maxgeoms, maxlayers, maxfine, max_user_levels,               & ! Inputs (dimensioning)
        Do_Thermset, do_deltam_scaling, do_PlanPar,                  & ! Inputs (Flags)
        do_regular_ps, do_enhanced_ps, doNadir,                      & ! Inputs (Flags)
        ngeoms, nlayers, nfinedivs, n_user_levels, user_levels,      & ! Inputs (control output)
        bb_input, surfbb, emiss,                                     & ! Inputs (Thermal)
        extinction, deltaus, omega, truncfac,                        & ! Inputs (Optical)
        Mu1, NCrit, Raycon, cota, xfine, wfine, csqfine, cotfine,    & ! Inputs (Geometry)
        intensity_dta_up, intensity_dts, cumsource_up, tcom1 )         ! Outputs
           else
              call DTE_Integral_I_DN &
      ( maxgeoms, maxlayers, maxfine, max_user_levels,           & ! Inputs (dimensioning)
        Do_Thermset, do_deltam_scaling, do_PlanPar,              & ! Inputs (Flags)
        do_regular_ps, do_enhanced_ps, doNadir,                  & ! Inputs (Flags)
        ngeoms, nlayers, nfinedivs, n_user_levels, user_levels,  & ! Inputs (control output)
        BB_input, extinction, deltaus, omega, truncfac,          & ! Inputs (Optical)
        Mu1, NCrit, RadCrit, CotCrit, Raycon, radii,             & ! Inputs (Geometry)
        cota, xfine, wfine, csqfine, cotfine,                    & ! Inputs (Geometry)
        intensity_dta_dn, cumsource_dn, tcom1 )                    ! Outputs
           endif

!  save results (TOA or BOA)

           if ( do_upwelling ) then
              intensity(1:ngeoms) = intensity_dta_up(1,1:ngeoms) + intensity_dts(1,1:ngeoms)
           else
              intensity(1:ngeoms) = intensity_dta_dn(n_user_levels,1:ngeoms)
           endif

!  End thermal

        endif

!  Solar and Thermal sources run
!  -----------------------------

        if ( do_solar_sources.and.do_thermal_emission ) then

!  SS Geometry call

           call SS_Geometry_1 &
       ( maxgeoms, maxlayers, maxfine, ngeoms, nlayers, nfinedivs, dtr, Pie, & ! Input
         vsign, eradius, heights, alpha_boa, theta_boa, phi_boa, do_Chapman, & ! Input
         do_LOSpaths, do_planpar, do_regular_ps, do_enhanced_ps, doNadir,    & ! Input flags
         doCrit, Acrit, extinction, Raycon, radii, alpha, cota,              & ! Input/Output
         xfine, wfine, csqfine, cotfine, alphafine, radiifine,               & ! Input/Output
         NCrit, RadCrit, CotCrit, cosscat, chapfacs,                         & ! Output
         sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,                 & ! Output
         fail, message, trace )                                                ! Output

           !if ( fail ) then
           !   write(*,'(a)')message
           !   write(*,'(a)')trace
           !   stop'SS GEOMETRY failed'
           !endif
           if ( fail ) return

!  Spherical functions call

           if ( do_upwelling ) then
              Call FO_ScalarSS_spherfuncs ( STARTER, MAXMOMENTS_INPUT, MAXGEOMS, &
                 NMOMENTS_INPUT, NGEOMS, DF1, DF2, COSSCAT, LEGPOLY_UP )
           else if ( do_dnwelling ) then
              Call FO_ScalarSS_spherfuncs ( STARTER, MAXMOMENTS_INPUT, MAXGEOMS, &
                 NMOMENTS_INPUT, NGEOMS, DF1, DF2, COSSCAT, LEGPOLY_DN )
           endif

!  SOLAR RT Call

           if ( do_upwelling ) then
              call SS_Integral_I_UP &
      ( maxgeoms, maxlayers, maxfine, maxmoments_input, max_user_levels,        & ! Inputs (dimensioning)
        do_deltam_scaling, do_PlanPar, do_regular_ps, do_enhanced_ps, doNadir,  & ! Inputs (Flags)
        ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels, & ! Inputs (control output)
        reflec, extinction, deltaus, omega, truncfac, phasmoms, flux,           & ! Inputs (Optical)
        Mu0, Mu1, LegPoly_up, NCrit, xfine, wfine, csqfine, cotfine,            & ! Inputs (Geometry)
        Raycon, cota, sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,       & ! Inputs (Geometry)
        intensity_up, intensity_db, cumsource_up )                                ! Outputs
           else if ( do_dnwelling ) then
              call SS_Integral_I_DN &
      ( maxgeoms, maxlayers, maxfine, maxmoments_input, max_user_levels,          & ! Inputs (dimensioning)
        do_deltam_scaling, do_PlanPar, do_regular_ps, do_enhanced_ps, doNadir,    & ! Inputs (Flags)
        ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,   & ! Inputs (control output)
        extinction, deltaus, omega, truncfac, phasmoms, flux,                     & ! Inputs (Optical)
        Mu1, LegPoly_dn, NCrit, RadCrit, CotCrit, xfine, wfine, csqfine, cotfine, & ! Inputs (Geometry)
        Raycon, radii, cota, sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,  & ! Inputs (Geometry)
        intensity_dn, cumsource_dn )                                                ! Outputs
           endif

!  save results (TOA or BOA)

           if ( do_upwelling ) then
              intensity(1:ngeoms) = intensity_up(1,1:ngeoms) + intensity_db(1,1:ngeoms)
           else
              intensity(1:ngeoms) = intensity_dn(n_user_levels,1:ngeoms)
           endif

!  DT Geometry call

           call DT_Geometry &
       ( maxgeoms, maxlayers, maxfine, ngeoms, nlayers, nfinedivs,            & ! Inputs
         dtr, eradius, heights, alpha_boa,                                    & ! Input
         do_LOSpaths, do_planpar, do_regular_ps, do_enhanced_ps, doNadir,     & ! Input
         doCrit, Acrit, extinction, Raycon, radii, alpha, cota, xfine, wfine, & ! Input/Output
         csqfine, cotfine, alphafine, radiifine, NCrit, RadCrit, CotCrit,     & ! Output
         fail, message, trace )                                                 ! Output

           !if ( fail ) then
           !   write(*,'(a)')message
           !   write(*,'(a)')trace
           !   stop'DT GEOMETRY failed'
           !endif
           if ( fail ) return

!  DT RT Call

           if ( do_upwelling ) then
              call DTE_Integral_I_UP &
      ( maxgeoms, maxlayers, maxfine, max_user_levels,               & ! Inputs (dimensioning)
        Do_Thermset, do_deltam_scaling, do_PlanPar,                  & ! Inputs (Flags)
        do_regular_ps, do_enhanced_ps, doNadir,                      & ! Inputs (Flags)
        ngeoms, nlayers, nfinedivs, n_user_levels, user_levels,      & ! Inputs (control output)
        bb_input, surfbb, emiss,                                     & ! Inputs (Thermal)
        extinction, deltaus, omega, truncfac,                        & ! Inputs (Optical)
        Mu1, NCrit, Raycon, cota, xfine, wfine, csqfine, cotfine,    & ! Inputs (Geometry)
        intensity_dta_up, intensity_dts, cumsource_up, tcom1 )         ! Outputs
           else
              call DTE_Integral_I_DN &
      ( maxgeoms, maxlayers, maxfine, max_user_levels,           & ! Inputs (dimensioning)
        Do_Thermset, do_deltam_scaling, do_PlanPar,              & ! Inputs (Flags)
        do_regular_ps, do_enhanced_ps, doNadir,                  & ! Inputs (Flags)
        ngeoms, nlayers, nfinedivs, n_user_levels, user_levels,  & ! Inputs (control output)
        BB_input, extinction, deltaus, omega, truncfac,          & ! Inputs (Optical)
        Mu1, NCrit, RadCrit, CotCrit, Raycon, radii,             & ! Inputs (Geometry)
        cota, xfine, wfine, csqfine, cotfine,                    & ! Inputs (Geometry)
        intensity_dta_dn, cumsource_dn, tcom1 )                    ! Outputs
           endif

!  save results (TOA or BOA). Add to existing

           if ( do_upwelling ) then
              intensity(1:ngeoms) = intensity(1:ngeoms) + &
                                    intensity_dta_up(1,1:ngeoms) + intensity_dts(1,1:ngeoms)
           else
              intensity(1:ngeoms) = intensity(1:ngeoms) + &
                                    intensity_dta_dn(n_user_levels,1:ngeoms)
           endif

!  End solar+thermal

        endif

!  Finish

   return
end subroutine FO_MASTER

end module FO_ScalarSS_Masters
