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

module FO_VectorSS_masters

!  All subroutines public

public

contains

subroutine VFO_MASTER &
       ( maxgeoms, maxlayers, maxfine, maxmoments_input, max_user_levels,    & ! Input max dims
         ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, nstokes, & ! Input dims
         do_solar_sources, do_thermal_emission, do_surface_emission,         & ! Input flags
         do_planpar, do_regular_ps, do_enhanced_ps, do_deltam_scaling,       & ! Input flags
         do_upwelling, do_dnwelling, do_lambertian, do_sunlight,             & ! Input flags
         dtr, Pie, doCrit, Acrit, eradius, heights, user_levels,             & ! Input general
         theta_boa, alpha_boa, phi_boa, Mu0, Mu1, doNadir,                   & ! Input geometry
         flux, fluxvec, extinction, deltaus, omega, greekmat, truncfac,      & ! Input atmos optical
         bb_input, reflec, surfbb, emiss,                                    & ! Input surf optical
         fo_stokes_ss, fo_stokes_db, fo_stokes_dta, fo_stokes_dts,           & ! Output
         fo_stokes, fail, message, trace )                                     ! Output

   USE FO_geometry_DTonly_m
   USE FO_geometry_SSonly_m

   USE FO_VectorSS_spherfuncs_m
   USE FO_VectorSS_RTCalcs_I_m

   USE FO_Thermal_RTCalcs_I_m

   implicit none

!  parameter arguments

   integer, parameter :: ffp = selected_real_kind(15),&
                         maxstokes = 4,&
                         max_directions = 2,&
                         upidx = 1, dnidx = 2

!  Subroutine inputs
!  =================

!  Max dimensions
!  --------------

   integer, intent(in)  :: maxgeoms, maxlayers, maxfine, maxmoments_input,&
                           max_user_levels

!  Dimensions
!  ----------

!  Layer and geometry control. Finelayer divisions may be changed

   integer, intent(in)    :: ngeoms, nlayers
   integer, intent(inout) :: nfinedivs(maxlayers,maxgeoms)
   integer, intent(in)    :: nmoments_input
   integer, intent(in)    :: n_user_levels

!  Number of Stokes components

   integer, intent(in)    :: nstokes

!  Configuration inputs
!  --------------------

!  Sources control, including thermal

   logical, intent(in)  :: do_solar_sources
   logical, intent(in)  :: do_thermal_emission
   logical, intent(in)  :: do_surface_emission

!  Flags (sphericity flags should be mutually exclusive)

   logical, intent(in)  :: do_planpar
   logical, intent(in)  :: do_regular_ps
   logical, intent(in)  :: do_enhanced_ps

!  deltam scaling flag

   logical, intent(in)  :: do_deltam_scaling

!  Directional Flags

   logical, intent(in)  :: do_upwelling, do_dnwelling

!  Lambertian surface flag

   logical, intent(in)  :: do_lambertian

!  Vector sunlight flag

   logical, intent(in)  :: do_sunlight

!  General inputs
!  --------------

!  DTR = degrees-to-Radians. Pie = 3.14159...

   real(ffp), intent(in) :: dtr, Pie

!  Critical adjustment for cloud layers

   logical, intent(inout) :: doCrit
   real(ffp), intent(in)  :: Acrit

!  Earth radius + heights

   real(ffp), intent(in) :: eradius, heights (0:maxlayers)

!  Output levels

   integer, intent(in)   :: user_levels ( max_user_levels )

!  Geometry inputs
!  ---------------

!  Input angles (Degrees)

   real(ffp), intent(inout) :: theta_boa(maxgeoms) !SZA
   real(ffp), intent(inout) :: alpha_boa(maxgeoms) !UZA
   real(ffp), intent(inout) :: phi_boa(maxgeoms)   !RAA

!  Mu0 = cos(theta_boa), required for surface term (both regular & enhanced)
!  Mu1 = cos(alpha_boa), required for the Regular PS only

   real(ffp), intent(in)    :: Mu0(maxgeoms), Mu1(maxgeoms)

!  Flag for the Nadir case. Intent(inout), input if DO_LOSpaths set)

   logical, intent(inout)   :: doNadir(maxgeoms)

!  Optical inputs
!  --------------

!  Solar flux

   real(ffp), intent(in) :: flux, fluxvec(maxstokes)

!  Atmosphere

   real(ffp), intent(in) :: extinction  ( maxlayers )
   real(ffp), intent(in) :: deltaus     ( maxlayers )
   real(ffp), intent(in) :: omega       ( maxlayers )
   real(ffp), intent(in) :: greekmat    ( maxlayers, 0:maxmoments_input, &
                                          maxstokes, maxstokes )

!  For TMS correction

   real(ffp), intent(in) :: truncfac    ( maxlayers )

!  Thermal inputs

   real(ffp), intent(in) :: bb_input ( 0:maxlayers )

!  Surface properties - reflective (could be the albedo)

   real(ffp), intent(in) :: reflec(maxstokes,maxstokes,maxgeoms)

!  Surface properties - emissive

   real(ffp), intent(in) :: surfbb
   real(ffp), intent(in) :: emiss ( maxgeoms )

!  Subroutine outputs
!  ==================

!  Solar

   real(ffp), intent(out) :: fo_stokes_ss ( max_user_levels,maxgeoms,maxstokes,max_directions )
   real(ffp), intent(out) :: fo_stokes_db ( max_user_levels,maxgeoms,maxstokes )

!  Thermal

   real(ffp), intent(out) :: fo_stokes_dta ( max_user_levels,maxgeoms,maxstokes,max_directions )
   real(ffp), intent(out) :: fo_stokes_dts ( max_user_levels,maxgeoms,maxstokes )

!  Composite

   !real(ffp), intent(out) :: fo_stokes_atmos ( max_user_levels,maxgeoms,maxstokes,max_directions )
   !real(ffp), intent(out) :: fo_stokes_surf  ( max_user_levels,maxgeoms,maxstokes )
   real(ffp), intent(out) :: fo_stokes       ( max_user_levels,maxgeoms,maxstokes,max_directions )

!  Exception handling

   logical, intent(out)           :: fail
   character (len=*), intent(out) :: message
   character (len=*), intent(out) :: trace

!  Other variables
!  ===============

!  Geometry routine outputs
!  ------------------------

!  VSIGN = +1 (Up); -1(Down)

   real(ffp)  :: vsign

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

   real(ffp)  :: DF1(MAXMOMENTS_INPUT)
   real(ffp)  :: DF2(MAXMOMENTS_INPUT)

!  Spherical functions, rotation angles

   real(ffp)  :: rotations_up(maxstokes,maxgeoms)
   real(ffp)  :: genspher_up(0:maxmoments_input,maxstokes,maxgeoms)
   real(ffp)  :: rotations_dn(maxstokes,maxgeoms)
   real(ffp)  :: genspher_dn(0:maxmoments_input,maxstokes,maxgeoms)
   real(ffp)  :: gshelp(7,0:maxmoments_input)

!  RT calculation outputs
!  ----------------------

!  Solar routines

   real(ffp)  :: stokes_up    ( max_user_levels,maxstokes,maxgeoms )
   real(ffp)  :: stokes_dn    ( max_user_levels,maxstokes,maxgeoms )
   real(ffp)  :: stokes_db    ( max_user_levels,maxstokes,maxgeoms )

!  Thermal routines (Scalar, no polarization)

   real(ffp)  :: intensity_dta_up ( max_user_levels,maxgeoms )
   real(ffp)  :: intensity_dta_dn ( max_user_levels,maxgeoms )
   real(ffp)  :: intensity_dts    ( max_user_levels,maxgeoms )

!  Intermediate RT products
!  ------------------------

!  Solar

   !real(ffp)  :: fo_stokes_ss ( max_user_levels,maxgeoms,maxstokes,max_directions )
   !real(ffp)  :: fo_stokes_db ( max_user_levels,maxgeoms,maxstokes )

!  Thermal

   !real(ffp)  :: fo_stokes_dta ( max_user_levels,maxgeoms,maxstokes,max_directions )
   !real(ffp)  :: fo_stokes_dts ( max_user_levels,maxgeoms,maxstokes )

!  Composite

   real(ffp)  :: fo_stokes_atmos ( max_user_levels,maxgeoms,maxstokes,max_directions )
   real(ffp)  :: fo_stokes_surf  ( max_user_levels,maxgeoms,maxstokes )

!  Other products
!  --------------

!  Thermal setup

   real(ffp)  :: tcom1(maxlayers,2)

!  Dummies

   real(ffp)  :: SScumsource_up ( 0:maxlayers,maxstokes,maxgeoms )
   real(ffp)  :: SScumsource_dn ( 0:maxlayers,maxstokes,maxgeoms )
   real(ffp)  :: DTcumsource_up ( 0:maxlayers,maxgeoms )
   real(ffp)  :: DTcumsource_dn ( 0:maxlayers,maxgeoms )

!  LOCAL
!  -----

!  numbers

   real(ffp), parameter :: zero = 0.0_ffp, one = 1.0_ffp

!  help variables

   integer   :: g,lev,n,o1
   logical   :: do_Chapman, do_lospaths, do_Thermset, starter

!  Initialize output

   fo_stokes_ss  = zero
   fo_stokes_db  = zero
   fo_stokes_dta = zero
   fo_stokes_dts = zero
   fo_stokes     = zero

   fail    = .false.
   message = ' '
   trace   = ' '

!  Flags to be set for each calculation (safety)

   do_Chapman  = .false.
   do_lospaths = .false.
   do_Thermset = .true.
   starter     = .true.

!  Solar sources run (NO THERMAL)
!  ------------------------------

   if ( do_solar_sources ) then

     if ( do_upwelling ) then
       vsign =  one

       !Geometry call
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
       !   stop 'SS GEOMETRY failed'
       !endif
       if ( fail ) return

       !Spherical functions call
       Call FO_VectorSS_spherfuncs &
       ( starter, maxmoments_input, maxgeoms, nmoments_input, ngeoms, nstokes, & ! inputs
         dtr, do_sunlight, theta_boa, alpha_boa, phi_boa, cosscat, vsign,      & ! inputs
         rotations_up, gshelp, genspher_up )                                     ! Outputs

       !RT call - solar only
       call SSV_Integral_I_UP &
       ( maxgeoms, maxlayers, maxfine, maxmoments_input, max_user_levels,            & ! Inputs (dimensioning)
         do_sunlight, do_deltam_scaling, do_lambertian,                              & ! Inputs (Flags)
         do_PlanPar, do_regular_ps, do_enhanced_ps, doNadir, nstokes,                & ! Inputs (Flags)
         ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,     & ! Inputs (control output)
         reflec, extinction, deltaus, omega, truncfac, greekmat, flux, fluxvec,      & ! Inputs (Optical)
         Mu0, Mu1, GenSpher_up, Rotations_up, NCrit, xfine, wfine, csqfine, cotfine, & ! Inputs (Geometry)
         Raycon, cota, sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,           & ! Inputs (Geometry)
         stokes_up, stokes_db, SScumsource_up )                                        ! Outputs

       !Save results
       do o1 = 1, nstokes
         do g = 1, ngeoms
           do lev=1,n_user_levels
             fo_stokes_ss(lev,g,o1,upidx) = stokes_up(lev,o1,g)
             fo_stokes_db(lev,g,o1)       = stokes_db(lev,o1,g)
             fo_stokes(lev,g,o1,upidx)    = fo_stokes_ss(lev,g,o1,upidx) &
                                          + fo_stokes_db(lev,g,o1)
           enddo
         enddo
       enddo
     endif

     if ( do_dnwelling ) then
       vsign =  -one

       !Geometry call
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
       !   stop 'SS GEOMETRY failed'
       !endif
       if ( fail ) return

       !Spherical functions call
       call FO_VectorSS_spherfuncs &
       ( starter, maxmoments_input, maxgeoms, nmoments_input, ngeoms, nstokes, & ! inputs
         dtr, do_sunlight, theta_boa, alpha_boa, phi_boa, cosscat, vsign,      & ! inputs
         rotations_dn, gshelp, genspher_dn )                                     ! Outputs

       !RT call - solar only
       call SSV_Integral_I_DN &
       ( maxgeoms, maxlayers, maxfine, maxmoments_input, max_user_levels, do_sunlight,    & ! Inputs (dimension)
         do_deltam_scaling, do_PlanPar, do_regular_ps, do_enhanced_ps, doNadir,           & ! Inputs (Flags)
         ngeoms, nstokes, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels, & ! Inputs (control)
         extinction, deltaus, omega, truncfac, greekmat, flux, fluxvec,                   & ! Inputs (Optical)
         Mu1, GenSpher_dn, Rotations_dn, NCrit, RadCrit, CotCrit,                         & ! Inputs (Geometry)
         xfine, wfine, csqfine, cotfine, Raycon, radii, cota,                             & ! Inputs (Geometry)
         sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,                              & ! Inputs (Geometry)
         stokes_dn, SScumsource_dn )                                                        ! Outputs

       !Save results
       do o1 = 1, nstokes
         do g = 1, ngeoms
           do lev=1,n_user_levels
             fo_stokes_ss(lev,g,o1,dnidx) = stokes_dn(lev,o1,g)
             fo_stokes(lev,g,o1,dnidx)    = fo_stokes_ss(lev,g,o1,dnidx)
           enddo
         enddo
       enddo
     endif

!  End solar only
   endif

!  Thermal sources run (NO SOLAR)
!  ------------------------------

   if ( do_thermal_emission ) then

     !Geometry call
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
     !   stop 'DT GEOMETRY failed'
     !endif
     if ( fail ) return

     if ( do_upwelling ) then
       !RT call - thermal only
       call DTE_Integral_I_UP &
       ( maxgeoms, maxlayers, maxfine, max_user_levels,               & ! Inputs (dimensioning)
         Do_Thermset, do_deltam_scaling, do_PlanPar,                  & ! Inputs (Flags)
         do_regular_ps, do_enhanced_ps, doNadir,                      & ! Inputs (Flags)
         ngeoms, nlayers, nfinedivs, n_user_levels, user_levels,      & ! Inputs (control output)
         bb_input, surfbb, emiss,                                     & ! Inputs (Thermal)
         extinction, deltaus, omega, truncfac,                        & ! Inputs (Optical)
         Mu1, NCrit, Raycon, cota, xfine, wfine, csqfine, cotfine,    & ! Inputs (Geometry)
         intensity_dta_up, intensity_dts, DTcumsource_up, tcom1 )       ! Outputs

       !Save results
       o1=1
       do g = 1, ngeoms
         do lev=1,n_user_levels
           fo_stokes_dta(lev,g,o1,upidx) = intensity_dta_up(lev,g)
           fo_stokes_dts(lev,g,o1)       = intensity_dts(lev,g)

           fo_stokes(lev,g,o1,upidx)     = fo_stokes_dta(lev,g,o1,upidx) &
                                         + fo_stokes_dts(lev,g,o1)
         enddo
       enddo

     !End upwelling
     endif

     if ( do_dnwelling ) then
       !RT call - thermal only
       call DTE_Integral_I_DN &
       ( maxgeoms, maxlayers, maxfine, max_user_levels,           & ! Inputs (dimensioning)
         Do_Thermset, do_deltam_scaling, do_PlanPar,              & ! Inputs (Flags)
         do_regular_ps, do_enhanced_ps, doNadir,                  & ! Inputs (Flags)
         ngeoms, nlayers, nfinedivs, n_user_levels, user_levels,  & ! Inputs (control output)
         BB_input, extinction, deltaus, omega, truncfac,          & ! Inputs (Optical)
         Mu1, NCrit, RadCrit, CotCrit, Raycon, radii,             & ! Inputs (Geometry)
         cota, xfine, wfine, csqfine, cotfine,                    & ! Inputs (Geometry)
         intensity_dta_dn, DTcumsource_dn, tcom1 )                  ! Outputs

       !Save results
       o1=1
       do g = 1, ngeoms
         do lev=1,n_user_levels
           fo_stokes_dta(lev,g,o1,dnidx) = intensity_dta_dn(lev,g)
           fo_stokes(lev,g,o1,dnidx)     = fo_stokes_dta(lev,g,o1,dnidx)
         enddo
       enddo
     endif

!  End thermal

   endif

!  Solar and Thermal sources run
!  -----------------------------

   if ( do_solar_sources .and. do_thermal_emission ) then
     !Add solar and thermal components
     if ( do_upwelling ) then
       o1=1
       do g = 1, ngeoms
         do lev=1,n_user_levels
           fo_stokes_atmos(lev,g,o1,upidx) = fo_stokes_ss(lev,g,o1,upidx) &
                                           + fo_stokes_dta(lev,g,o1,upidx)
           fo_stokes_surf(lev,g,o1)        = fo_stokes_db(lev,g,o1)       &
                                           + fo_stokes_dts(lev,g,o1)

           fo_stokes(lev,g,o1,upidx)       = fo_stokes_atmos(lev,g,o1,upidx) &
                                           + fo_stokes_surf(lev,g,o1)
         enddo
       enddo
     endif

     if ( do_dnwelling ) then
       o1=1
       do g = 1, ngeoms
         do lev=1,n_user_levels
           fo_stokes_atmos(lev,g,o1,dnidx) = fo_stokes_ss(lev,g,o1,dnidx) &
                                           + fo_stokes_dta(lev,g,o1,dnidx)

           fo_stokes(lev,g,o1,dnidx)       = fo_stokes_atmos(lev,g,o1,dnidx)
         enddo
       enddo
     endif
   endif

!  Finish

end subroutine VFO_MASTER

end module FO_VectorSS_masters
