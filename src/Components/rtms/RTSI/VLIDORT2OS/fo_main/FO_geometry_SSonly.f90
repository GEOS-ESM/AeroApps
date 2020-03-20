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

module FO_geometry_SSonly_m

!  Stand alone geometry for solar scattering only

!   Plane-parallel option, September 2012               (Version 1.2)
!   Extension to multiple geometries, 29 October 2012   (Version 1.3)

   use FO_geometry_Pool_m

private
public SS_Geometry_1, SS_Geometry_2, RegularPS_sphergeom, PlanePargeom

contains

subroutine SS_Geometry_1 &
       ( maxgeoms, maxlayers, maxfine, ngeoms, nlayers, nfinedivs, dtr, Pie, & ! Input
         vsign, eradius, heights, alpha_boa, theta_boa, phi_boa, do_Chapman, & ! Input
         do_LOSpaths, do_planpar, do_regular_ps, do_enhanced_ps, doNadir,    & ! Input flags
         doCrit, Acrit, extinc, Raycon, radii, alpha, cota,                  & ! Input/Output
         xfine, wfine, csqfine, cotfine, alphafine, radiifine,               & ! Input/Output
         NCrit, RadCrit, CotCrit, cosscat, chapfacs,                         & ! Output
         sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,                 & ! Output
         fail, message, trace )                                                ! Output

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Input arguments
!  ===============

!  Dimensions

   integer  , intent(in)     :: maxgeoms, maxlayers, maxfine

!  Layer and geometry control. Finelayer divisions may be changed

   integer  , intent(in)     :: ngeoms, nlayers
   integer  , intent(inout)  :: nfinedivs(maxlayers,maxgeoms)

!  Radius + heights

   real(ffp), intent(in)     :: eradius, heights (0:maxlayers)

!  input angles (Degrees), dtr = degrees-to-Radians. VSIGN = +1 (Up); -1(Down)

   real(ffp), intent(in)     :: dtr, Pie, vsign
   real(ffp), intent(InOut)  :: alpha_boa(maxgeoms), theta_boa(maxgeoms), phi_boa(maxgeoms)

!  LOS paths flag is new, 20 January 2012

   logical  , intent(in)     :: do_Chapman
   logical  , intent(in)     :: do_LOSpaths

!  Flags (sphericity flags should be mutually exclusive)

   logical  , intent(in)     :: do_regular_ps
   logical  , intent(in)     :: do_enhanced_ps
   logical  , intent(in)     :: do_planpar

!  Critical adjustment for cloud layers

   logical  , intent(inout)  :: doCrit
   real(ffp), intent(in)     :: extinc(maxlayers)
   real(ffp), intent(in)     :: Acrit

!  Input/Output Arguments (These may already be set if do_LOSpaths )
!  ======================

!  Flag for the Nadir case. Intent(inout), input if DO_LOSpaths set)

   logical  , intent(inout)  :: doNadir(maxgeoms)
  
!  Alphas,  Cotangents, Radii, Ray constant. Intent(inout), input if DO_LOSpaths set)

   real(ffp), intent(inout)  :: radii    (0:maxlayers)
   real(ffp), intent(inout)  :: Raycon   (maxgeoms)
   real(ffp), intent(inout)  :: alpha    (0:maxlayers,maxgeoms)
   real(ffp), intent(inout)  :: cota     (0:maxlayers,maxgeoms)

!  LOS Quadratures for Enhanced PS

   real(ffp), intent(inout)  :: xfine    (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: wfine    (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: csqfine  (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: cotfine  (maxlayers,maxfine,maxgeoms)

!  Fine layering output

   real(ffp), intent(inout)  :: alphafine (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: radiifine (maxlayers,maxfine,maxgeoms)

!  Output arguments
!  ================

!  Critical layer

   integer  , intent(out)  :: Ncrit(maxgeoms)
   real(ffp), intent(out)  :: RadCrit(maxgeoms), CotCrit(maxgeoms)

!  solar paths 

   integer  , Intent(out)  :: ntraverse  (0:maxlayers,maxgeoms)
   real(ffp), Intent(out)  :: sunpaths   (0:maxlayers,maxlayers,maxgeoms)
   integer  , Intent(out)  :: ntraverse_fine(maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(out)  :: sunpaths_fine (maxlayers,maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(out)  :: Chapfacs   (maxlayers,maxlayers,maxgeoms)

!  Cosine scattering angle

   real(ffp), Intent(out)  :: cosscat(maxgeoms)

!  Exception handling

   logical      , intent(out)    :: fail
   character*(*), intent(out)    :: message
   character*(*), intent(out)    :: trace

!  Local arguments
!  ===============

!  LOS path lengths

   real(ffp)  :: Lospaths (maxlayers,maxgeoms)

!  Other angles

   real(ffp)  :: theta_all  (0:maxlayers,maxgeoms)
   real(ffp)  :: phi_all    (0:maxlayers,maxgeoms)
   real(ffp)  :: cosa       (0:maxlayers,maxgeoms)
   real(ffp)  :: sina       (0:maxlayers,maxgeoms)

!  Critical values

   real(ffp)  :: AlphaCrit(maxgeoms)

!  Help variables

   real(ffp), parameter   :: zero = 0.0_ffp
   real(ffp), parameter   :: one  = 1.0_ffp
   real(ffp)              :: term1(maxgeoms), term2(maxgeoms), cutoff
   character*2            :: c2
   integer                :: v

!  Initialize output

   fail = .false. ; message = ' ' ; trace = ' '
   NCrit = 0 ; AlphaCrit = zero ; RadCrit = zero ; CotCrit = zero
   cutoff = zero; if (ACrit.gt.zero) cutoff = -log(ACrit)

!  check range of inputs
!  ---------------------

!  Cannot have Plane-parallel and Regular PS

   if ( do_planpar .and. do_regular_ps ) then
      message = 'Cannot have BOTH Plane-parallel and Regular PS options'
      trace   = 'Initial Flag Check in SS_Geometry_1'
      fail    = .true. ;  return
   endif

!  Cannot have Plane-parallel and Enhanced PS

   if ( do_planpar .and. do_enhanced_ps ) then
      message = 'Cannot have BOTH Plane-parallel and Enhanced PS options'
      trace   = 'Initial Flag Check in SS_Geometry_1'
      fail    = .true. ;  return
   endif

!  Cannot have Regular PS and Enhanced PS

   if ( do_regular_ps .and. do_enhanced_ps ) then
      message = 'Cannot have BOTH Regular and Enhanced PS options'
      trace   = 'Initial Flag Check in SS_Geometry_1'
      fail    = .true. ;  return
   endif

!  Check geometry angles
!  ---------------------

   do v = 1, ngeoms

!  VZA can be 0-90 degrees inclusive, but not outside this range

      if ( alpha_boa(v).gt.90.0_ffp.or.alpha_boa(v).lt.zero ) then
         write(c2,'(I2)')v
         message = 'Boa LOS angle outside range [0,90]); Check it!'
         trace   = 'Geometry # '//c2//'; Initial Angle Check in SS_Geometry_1'
         fail    = .true. ;  return
      endif

!  PHI is not limited to <= 180 degs. Also, not negative.
!     Old Code :     if ( phi_boa.gt.180.0_ffp ) phi_boa = 360.0_ffp - phi_boa

      if ( phi_boa(v).lt.zero )   phi_boa(v) = - phi_boa(v)
      if ( phi_boa(v).gt.360.0_ffp ) phi_boa(v) = 360.0_ffp - phi_boa(v) - 360.0_ffp

!  SZA can be 0-90 degrees inclusive, but not outside this range
!    For plane-parallel, 90 degrees is not allowed

      if ( do_planpar ) then
         if ( theta_boa(v).ge.90.0_ffp.or.theta_boa(v).lt.zero ) then
            write(c2,'(I2)')v
            message = 'Plane-parallel: Boa SZA angle outside range [0,90)); Check it!'
            trace   = 'Geometry # '//c2//'; Initial Angle Check in SS_Geometry_1'
            fail    = .true. ;  return
         endif
      else
         if ( theta_boa(v).gt.90.0_ffp.or.theta_boa(v).lt.zero ) then
            write(c2,'(I2)')v
            message = 'Pseudo-spherical : Boa SZA angle outside range [0,90]); Check it!'
            trace   = 'Geometry # '//c2//'; Initial Angle Check in SS_Geometry_1'
            fail    = .true. ;  return
         endif
      endif

   enddo

!  Plane-parallel, One routine only
!  --------------------------------

   if ( do_planpar ) then
      CALL PlanePargeom &
       ( maxgeoms, maxlayers, ngeoms, nlayers, heights,          & ! Inputs
         do_Chapman, alpha_boa, theta_boa, phi_boa, vsign, dtr,  & ! Inputs
         alpha, sunpaths, ntraverse,                             & ! Outputs
         chapfacs, cosscat, term1, term2 )                         ! Outputs
      return
   endif

!  Regular PS, One routine only
!  ----------------------------

   if ( do_regular_ps ) then
      CALL RegularPS_sphergeom &
       ( maxgeoms, maxlayers, ngeoms, nlayers, heights, eradius,        & ! Inputs
         do_Chapman, alpha_boa, theta_boa, phi_boa, vsign, dtr,         & ! Inputs
         Raycon, radii, alpha, sunpaths, ntraverse,                     & ! Outputs
         chapfacs, cosscat, term1, term2 )                                ! Outputs
      return
   endif

!  Enhanced PS; proceed in 4 Steps
!  ===============================

!  Step 1; Initial LOS-path quantities
!  -----------------------------------

!    Given heights and BOA LOS angle, compute path angles and radii
!   Only need to do this if LOSpaths flag is not set.

   if ( .not. do_LOSpaths ) then
      CALL STD_outgoing_sphergeom_Initial &
       ( maxgeoms, maxlayers, ngeoms, nlayers, heights, eradius, alpha_boa, dtr, & ! Input
         doNadir, radii, Raycon, Lospaths, alpha, sina, cosa, cota )               ! Output
   endif

!  Step 2, Adjust for Criticality
!  ------------------------------
   
!  Step 2a; Outgoing, Find Critical-layer adjustments (Optional)

   if ( doCrit) then
      CALL STD_outgoing_FindCritlayer &
       ( maxgeoms, maxlayers, ngeoms, nlayers, Acrit, Cutoff, doNadir, & ! Inputs
         extinc, Lospaths, sina, cosa, radii, nfinedivs,               & ! Input
         Ncrit, AlphaCrit, RadCrit, CotCrit, fail, message )             ! Outputs
      if ( Fail ) then
         trace = 'Error from STD_Outgoing_FindCritLayer in SS_Geometry_1' ; return
      endif
   endif

!  Step 2b; Incoming, Find Critical-layer adjustments (Optional)

   if ( doCrit) then
      call SD_incoming_FindCritLayer &
       ( maxgeoms, maxlayers, ngeoms, nlayers, doNadir, dtr, Acrit, & ! Input
         cutoff, alpha, radii, extinc, Raycon, theta_boa,           & ! Inputs
         doCrit, Ncrit, nfinedivs, AlphaCrit, RadCrit, CotCrit,     & ! Outputs
         fail, message )                                              ! Outputs
      if ( Fail ) then
         trace = 'Error from SD_incoming_FindCritLayer in SS_Geometry_1' ; return
      endif
   endif

!  Step 3. Set Quadratures
!  -----------------------

!  Step 3a. LOS fine-layer quadratures (Regular, non-adjusted, no Criticality)
!           Regular quadratures may be avoided if do_LOSpaths is set

   if ( .not. doCrit) then
     CALL STD_outgoing_sphergeom_Qbasic &
       ( maxgeoms, maxlayers, maxfine, ngeoms, nlayers,     & ! Input
         nfinedivs, doNadir, radii, alpha, Raycon,          & ! Input
         radiifine, alphafine, xfine, wfine,                & ! Output
         csqfine, cotfine )                                   ! Output
   endif

!  Step 3b. LOS fine-layer quadratures
!           Critical-layer adjustment of Quadrature done here.
!           Regular quadratures may be avoided if do_LOSpaths is set

   if ( doCrit) then
      CALL STD_outgoing_sphergeom_Qadjusted &
       ( maxgeoms, maxlayers, ngeoms, maxfine, nlayers,         & ! Input
         nfinedivs, do_LOSpaths, doNadir, radii, alpha, Raycon, & ! Input
         doCrit, Ncrit, AlphaCrit, RadCrit,                     & ! Input
         radiifine, alphafine, xfine, wfine,                    & ! Output
         csqfine, cotfine )                                       ! Output
   endif

!  Step 4. Solar path lengths
!  --------------------------

   CALL SD_incoming_sphergeom &
       ( maxgeoms, maxlayers, maxfine, ngeoms, nlayers, nfinedivs, do_Chapman, & ! Input
         doNadir, DoCrit, NCrit, alpha_boa, theta_boa, phi_boa, radii, alpha,  & ! Input
         vsign, dtr, Pie, RadCrit, AlphaCrit, radiifine, alphafine,            & ! Input
         sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,                   & ! Output
         chapfacs, cosscat, theta_all, phi_all )                                 ! Output

!  Finish

   return
end subroutine SS_Geometry_1

!

subroutine SS_Geometry_2  &
       ( maxgeoms, maxlayers, maxfine, ngeoms, nlayers, nfinedivs, dtr, Pie, & ! Input
         vsign, eradius, heights, alpha_boa, theta_boa, phi_boa, do_Chapman, & ! Input
         do_LOSpaths, do_planpar, do_regular_ps, do_enhanced_ps, doNadir,    & ! Input
         doCrit, Acrit, extinc, Raycon, radii, alpha, cota,                  & ! Input/Output
         xfine, wfine, csqfine, cotfine, alphafine, radiifine,               & ! Input/Output
         NCrit, RadCrit, CotCrit,                                            & ! Output
         cosscat_up, chapfacs_up, cosscat_dn, chapfacs_dn,                   & ! Output
         sunpaths_up, ntraverse_up, sunpaths_fine_up, ntraverse_fine_up,     & ! Output
         sunpaths_dn, ntraverse_dn, sunpaths_fine_dn, ntraverse_fine_dn,     & ! Output
         fail, message, trace )                                                ! Output

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Input arguments
!  ===============

!  Dimensions

   integer  , intent(in)     :: maxgeoms, maxlayers, maxfine

!  Layer and geometry control. Finelayer divisions may be changed

   integer  , intent(in)     :: ngeoms, nlayers
   integer  , intent(inout)  :: nfinedivs(maxlayers,maxgeoms)

!  Radius + heights

   real(ffp), intent(in)     :: eradius, heights (0:maxlayers)

!  input angles (Degrees), dtr = degrees-to-Radians

   real(ffp), intent(in)     :: dtr, Pie
   real(ffp), intent(InOut)  :: alpha_boa(maxgeoms), theta_boa(maxgeoms), phi_boa(maxgeoms)

!  LOS paths flag is new, 20 January 2012

   logical  , intent(in)     :: do_Chapman
   logical  , intent(in)     :: do_LOSpaths

!  Flags (sphericity flags should be mutually exclusive)

   logical  , intent(in)     :: do_regular_ps
   logical  , intent(in)     :: do_enhanced_ps
   logical  , intent(in)     :: do_planpar

!  Critical adjustment for cloud layers

   logical  , intent(inout)  :: doCrit
   real(ffp), intent(in)     :: extinc(maxlayers)
   real(ffp), intent(in)     :: Acrit

!  Input/Output Arguments (These may already be set if do_LOSpaths )
!  ======================

!  Flag for the Nadir case. Intent(inout), input if DO_LOSpaths set)

   logical  , intent(inout)  :: doNadir(maxgeoms)
  
!  Alphas,  Cotangents, Radii, Ray constant. Intent(inout), input if DO_LOSpaths set)

   real(ffp), intent(inout)  :: radii    (0:maxlayers)
   real(ffp), intent(inout)  :: Raycon   (maxgeoms)
   real(ffp), intent(inout)  :: alpha    (0:maxlayers,maxgeoms)
   real(ffp), intent(inout)  :: cota     (0:maxlayers,maxgeoms)

!  LOS Quadratures for Enhanced PS

   real(ffp), intent(inout)  :: xfine    (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: wfine    (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: csqfine  (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: cotfine  (maxlayers,maxfine,maxgeoms)

!  Fine layering output

   real(ffp), intent(inout)  :: alphafine (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: radiifine (maxlayers,maxfine,maxgeoms)

!  Output arguments
!  ================

!  Critical layer

   integer  , intent(out)  :: Ncrit(maxgeoms)
   real(ffp), intent(out)  :: RadCrit(maxgeoms), CotCrit(maxgeoms)

!  solar paths upwelling

   integer  , Intent(out)  :: ntraverse_up  (0:maxlayers,maxgeoms)
   real(ffp), Intent(out)  :: sunpaths_up   (0:maxlayers,maxlayers,maxgeoms)
   integer  , Intent(out)  :: ntraverse_fine_up(maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(out)  :: sunpaths_fine_up (maxlayers,maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(out)  :: Chapfacs_up  (maxlayers,maxlayers,maxgeoms)

!  solar paths downwelling

   integer  , Intent(out)  :: ntraverse_dn  (0:maxlayers,maxgeoms)
   real(ffp), Intent(out)  :: sunpaths_dn   (0:maxlayers,maxlayers,maxgeoms)
   integer  , Intent(out)  :: ntraverse_fine_dn(maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(out)  :: sunpaths_fine_dn (maxlayers,maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(out)  :: Chapfacs_dn  (maxlayers,maxlayers,maxgeoms)

!  Cosine scattering angles

   real(ffp), Intent(out)  :: cosscat_up(maxgeoms), cosscat_dn(maxgeoms)

!  Exception handling

   logical      , intent(out)    :: fail
   character*(*), intent(out)    :: message
   character*(*), intent(out)    :: trace

!  Local arguments
!  ===============

!  LOS path lengths

   real(ffp)  :: Lospaths (maxlayers,maxgeoms)

!  Other angles

   real(ffp)  :: theta_all  (0:maxlayers,maxgeoms)
   real(ffp)  :: phi_all    (0:maxlayers,maxgeoms)
   real(ffp)  :: cosa       (0:maxlayers,maxgeoms)
   real(ffp)  :: sina       (0:maxlayers,maxgeoms)

!  Critical values

   real(ffp)  :: AlphaCrit(maxgeoms)

!  Help variables
!    VSIGN = +1 for Upwelling, -1 for Downwelling

   real(ffp), parameter   :: zero = 0.0_ffp
   real(ffp), parameter   :: one  = 1.0_ffp
   real(ffp)              :: vsign, term1(maxgeoms), term2(maxgeoms), cutoff
   character*2            :: c2
   integer                :: v

!  Initialize output

   fail = .false. ; message = ' ' ; trace = ' '
   NCrit = 0 ; AlphaCrit = zero ; RadCrit = zero ; CotCrit = zero
   cutoff = -log(ACrit)

!  check range of inputs
!  ---------------------

!  Cannot have Plane-parallel and Regular PS

   if ( do_planpar .and. do_regular_ps ) then
      message = 'Cannot have BOTH Plane-parallel and Regular PS options'
      trace   = 'Initial Flag Check in SS_Geometry_2'
      fail    = .true. ;  return
   endif

!  Cannot have Plane-parallel and Enhanced PS

   if ( do_planpar .and. do_enhanced_ps ) then
      message = 'Cannot have BOTH Plane-parallel and Enhanced PS options'
      trace   = 'Initial Flag Check in SS_Geometry_1'
      fail    = .true. ;  return
   endif

!  Cannot have Regular PS and Enhanced PS

   if ( do_regular_ps .and. do_enhanced_ps ) then
      message = 'Cannot have BOTH Regular and Enhanced PS options'
      trace   = 'Initial Flag Check in SS_Geometry_1'
      fail    = .true. ;  return
   endif

!  Check geometry angles
!  ---------------------

   do v = 1, ngeoms

!  VZA can be 0-90 degrees inclusive, but not outside this range

      if ( alpha_boa(v).gt.90.0_ffp.or.alpha_boa(v).lt.zero ) then
         write(c2,'(I2)')v
         message = 'Boa LOS angle outside range [0,90]); Check it!'
         trace   = 'Geometry # '//c2//'; Initial Angle Check in SS_Geometry_1'
         fail    = .true. ;  return
      endif

!  PHI is not limited to <= 180 degs. Also, not negative.
!     Old Code :     if ( phi_boa.gt.180.0_ffp ) phi_boa = 360.0_ffp - phi_boa

      if ( phi_boa(v).lt.zero )   phi_boa(v) = - phi_boa(v)
      if ( phi_boa(v).gt.360.0_ffp ) phi_boa(v) = 360.0_ffp - phi_boa(v) - 360.0_ffp

!  SZA can be 0-90 degrees inclusive, but not outside this range
!    For plane-parallel, 90 degrees is not allowed

      if ( do_planpar ) then
         if ( theta_boa(v).ge.90.0_ffp.or.theta_boa(v).lt.zero ) then
            write(c2,'(I2)')v
            message = 'Plane-parallel: Boa SZA angle outside range [0,90)); Check it!'
            trace   = 'Geometry # '//c2//'; Initial Angle Check in SS_Geometry_1'
            fail    = .true. ;  return
         endif
      else
         if ( theta_boa(v).gt.90.0_ffp.or.theta_boa(v).lt.zero ) then
            write(c2,'(I2)')v
            message = 'Pseudo-spherical : Boa SZA angle outside range [0,90]); Check it!'
            trace   = 'Geometry # '//c2//'; Initial Angle Check in SS_Geometry_1'
            fail    = .true. ;  return
         endif
      endif

   enddo

!  Plane-parallel, One routine only
!  --------------------------------

!  Single call to upwelling (VSIGN = +1)
!  Downwelling only requires new COSSCAT, rest is copied

   if ( do_planpar ) then
      vsign = +one
      CALL PlanePargeom &
       ( maxgeoms, maxlayers, ngeoms, nlayers, heights,         & ! Inputs
         do_Chapman, alpha_boa, theta_boa, phi_boa, vsign, dtr, & ! Inputs
         alpha, sunpaths_up, ntraverse_up,                      & ! Outputs
         chapfacs_up, cosscat_up, term1, term2 )                  ! Outputs
      vsign = -one
      do v = 1, ngeoms
         cosscat_dn(v) = -vsign * term1(v) + term2(v)
         ntraverse_dn(:,v) = ntraverse_up(:,v)
         sunpaths_dn(:,:,v)  = sunpaths_up(:,:,v)
         chapfacs_dn(:,:,v)  = chapfacs_up(:,:,v)
      enddo
      return
   endif

!  Regular PS, One routine only
!  ----------------------------

!  Single call to upwelling (VSIGN = +1)
!  Downwelling only requires new COSSCAT, rest is copied

   if ( do_regular_ps ) then
      vsign = +one
      CALL RegularPS_sphergeom &
       ( maxgeoms, maxlayers, ngeoms, nlayers, heights, eradius,  & ! Inputs
         do_Chapman, alpha_boa, theta_boa, phi_boa, vsign, dtr,   & ! Inputs
         Raycon, radii, alpha, sunpaths_up, ntraverse_up,         & ! Outputs
         chapfacs_up, cosscat_up, term1, term2 )                    ! Outputs
      vsign = -one
      do v = 1, ngeoms
         cosscat_dn(v) = -vsign * term1(v) + term2(v)
         ntraverse_dn(:,v) = ntraverse_up(:,v)
         sunpaths_dn(:,:,v)  = sunpaths_up(:,:,v)
         chapfacs_dn(:,:,v)  = chapfacs_up(:,:,v)
      enddo
      return
   endif

!  Enhanced PS; proceed in 4 Steps
!  ===============================

!  Step 1; Initial LOS-path quantities
!  -----------------------------------

!    Given heights and BOA LOS angle, compute path angles and radii
!   Only need to do this if LOSpaths flag is not set.

   if ( .not. do_LOSpaths ) then
      CALL STD_outgoing_sphergeom_Initial &
       ( maxgeoms, maxlayers, ngeoms, nlayers, heights, eradius, alpha_boa, dtr, & ! Input
         doNadir, radii, Raycon, Lospaths, alpha, sina, cosa, cota )               ! Output
   endif

!  Step 2, Adjust for Criticality
!  ------------------------------
   
!  Step 2a; Outgoing, Find Critical-layer adjustments (Optional)

   if ( doCrit) then
      CALL STD_outgoing_FindCritlayer &
       ( maxgeoms, maxlayers, ngeoms, nlayers, Acrit, Cutoff, doNadir, & ! Inputs
         extinc, Lospaths, sina, cosa, radii, nfinedivs,               & ! Input
         Ncrit, AlphaCrit, RadCrit, CotCrit, fail, message )             ! Outputs
      if ( Fail ) then
         trace = 'Error from STD_Outgoing_FindCritLayer in SS_Geometry_1' ; return
      endif
   endif

!  Step 2b; Incoming, Find Critical-layer adjustments (Optional)

   if ( doCrit) then
      call SD_incoming_FindCritLayer &
       ( maxgeoms, maxlayers, ngeoms, nlayers, doNadir, dtr, Acrit, & ! Input
         cutoff, alpha, radii, extinc, Raycon, theta_boa,           & ! Inputs
         doCrit, Ncrit, nfinedivs, AlphaCrit, RadCrit, CotCrit,     & ! Outputs
         fail, message )                                              ! Outputs
      if ( Fail ) then
         trace = 'Error from SD_incoming_FindCritLayer in SS_Geometry_1' ; return
      endif
   endif

!  Step 3. Set Quadratures
!  -----------------------

!  Step 3a. LOS fine-layer quadratures (Regular, non-adjusted, no Criticality)
!           Regular quadratures may be avoided if do_LOSpaths is set

   if ( .not. doCrit) then
      CALL STD_outgoing_sphergeom_Qbasic &
       ( maxgeoms, maxlayers, maxfine, ngeoms, nlayers,     & ! Input
         nfinedivs, doNadir, radii, alpha, Raycon,          & ! Input
         radiifine, alphafine, xfine, wfine,                & ! Output
         csqfine, cotfine )                                   ! Output
   endif

!  Step 3b. LOS fine-layer quadratures
!           Critical-layer adjustment of Qudarature done here.
!           Regular quadratures may be avoided if do_LOSpaths is set

   if ( doCrit) then
      CALL STD_outgoing_sphergeom_Qadjusted &
       ( maxgeoms, maxlayers, ngeoms, maxfine, nlayers,         & ! Input
         nfinedivs, do_LOSpaths, doNadir, radii, alpha, Raycon, & ! Input
         doCrit, Ncrit, AlphaCrit, RadCrit,                     & ! Input
         radiifine, alphafine, xfine, wfine,                    & ! Output
         csqfine, cotfine )                                       ! Output
   endif

!  Step 4. Solar path lengths
!  --------------------------

!  Step 4a. Solar path lengths, Upwelling

   vsign = + one
   CALL SD_incoming_sphergeom &
       ( maxgeoms, maxlayers, ngeoms, maxfine, nlayers, nfinedivs, do_Chapman, & ! Input
         doNadir, DoCrit, NCrit, alpha_boa, theta_boa, phi_boa, radii, alpha,  & ! Input
         vsign, dtr, Pie, RadCrit, AlphaCrit, radiifine, alphafine,            & ! Input
         sunpaths_up, ntraverse_up, sunpaths_fine_up, ntraverse_fine_up,       & ! Output
         chapfacs_up, cosscat_up, theta_all, phi_all )                           ! Output

!  Step 4b. Solar path lengths, Downwelling

   vsign = - one
   CALL SD_incoming_sphergeom &
       ( maxgeoms, maxlayers, ngeoms, maxfine, nlayers, nfinedivs, do_Chapman, & ! Input
         doNadir, DoCrit, NCrit, alpha_boa, theta_boa, phi_boa, radii, alpha,  & ! Input
         vsign, dtr, Pie, RadCrit, AlphaCrit, radiifine, alphafine,            & ! Input
         sunpaths_dn, ntraverse_dn, sunpaths_fine_dn, ntraverse_fine_dn,       & ! Output
         chapfacs_dn, cosscat_dn, theta_all, phi_all )                           ! Output

!  Finish

   return
end subroutine SS_Geometry_2

!

subroutine Planepargeom &
       ( maxgeoms, maxlayers, ngeoms, nlayers, heights,         & ! Inputs
         do_Chapman, alpha_boa, theta_boa, phi_boa, vsign, dtr, & ! Inputs
         alpha, sunpaths, ntraverse,                            & ! Outputs
         chapfacs, cosscat, term1, term2 )                        ! Outputs

!  Completely stand-alone geometry routine for Accurate SS
!     This is for the Plane-parallel choice
!     This is applicable to the Upwelling and/or/Downwelling LOS-path geometries
!     No partials, this routine

!    starting inputs are the BOA values of SZA, VZA and PHI
!    need also the height grids, earth radius and control

      implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs

   integer  , intent(In)    :: maxlayers, maxgeoms
   integer  , intent(In)    :: nlayers, ngeoms
   real(ffp), intent(InOut) :: alpha_boa(maxgeoms), theta_boa(maxgeoms), phi_boa(maxgeoms)
   logical  , intent(In)    :: do_Chapman
   real(ffp), intent(In)    :: dtr, heights (0:maxlayers), vsign

!  Los geometry

   real(ffp), intent(Out)  :: alpha         (0:maxlayers,maxgeoms)

!  main outputs (geometry)

   integer  , intent(Out)  :: ntraverse  (0:maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: sunpaths   (0:maxlayers,maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: chapfacs   (maxlayers,maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: cosscat(maxgeoms), term1(maxgeoms), term2(maxgeoms)

!  Local

   logical       :: Do_OverheadSun
   integer       :: n, v
   real(ffp)     :: alpha_boa_R, theta_boa_R
   real(ffp)     :: salpha_boa, calpha_boa
   real(ffp)     :: stheta_boa, ctheta_boa, utheta_boa, cphi_boa, diffhts(maxlayers)

   real(ffp), parameter  :: zero = 0.0_ffp
   real(ffp), parameter  :: one  = 1.0_ffp

!  Initialise output

   alpha = zero
   ntraverse = 0
   sunpaths  = zero ; chapfacs  = zero
   cosscat   = zero ; term1 = zero ; term2 = zero

!  start geometry loop

   do v = 1, ngeoms

!  BOA angles

      alpha_boa_R    = alpha_boa(v) * DTR
      if ( alpha_boa(v).eq.90.0_ffp ) then
         calpha_boa     = zero
         salpha_boa     = one
      else
         calpha_boa     = cos(alpha_boa_R)
         salpha_boa     = sin(alpha_boa_R)
      endif

      theta_boa_R    = theta_boa(v) * DTR
      stheta_boa     = sin(theta_boa_R)
      ctheta_boa     = cos(theta_boa_R)

      cphi_boa       = cos(phi_boa(v) * dtr)

!  Nominal traverse paths for Full illumination. Difference heights

      do n = 1, nlayers
         diffhts(n) = heights(n-1) - heights(n)
         ntraverse(n,v) = n
      enddo

!  Overhead Sun

      Do_OverheadSun = (theta_boa(v).eq.zero) 

!  Set Alpha, scattering angle

      alpha(1:nlayers,v) = alpha_boa_R
      if ( Do_OverheadSun ) then
         term1(v) = zero
         term2(v) = calpha_boa
         cosscat(v) = - vsign * term2(v) ; if (term2(v).eq.zero) cosscat(v) = term2(v)
      else
         term1(v) = salpha_boa * stheta_boa * cphi_boa
         term2(v) = calpha_boa * ctheta_boa
         cosscat(v) = - vsign * term2(v) + term1(v) 
      endif

!  Sunpath/Chapman factor calculations

      utheta_boa     = one / ctheta_boa
      do n = 1, nlayers
         sunpaths(n,1:n,v) = diffhts(1:n) * utheta_boa
         if ( do_Chapman ) chapfacs(n,1:n,v) = utheta_boa
      enddo

!  End geometry routine

   enddo

!  Finish

   return
end subroutine Planepargeom

!

subroutine RegularPS_sphergeom &
       ( maxgeoms, maxlayers, ngeoms, nlayers, heights, eradius, & ! Inputs
         do_Chapman, alpha_boa, theta_boa, phi_boa, vsign, dtr,  & ! Inputs
         Raycon, radii, alpha, sunpaths, ntraverse,              & ! Outputs
         chapfacs, cosscat, term1, term2 )                         ! Outputs

!  Completely stand-alone geometry routine for Accurate SS
!     This is for the Regular PS choice
!     This is applicable to the Upwelling and/or/Downwelling LOS-path geometries
!     No partials, this routine

!    starting inputs are the BOA values of SZA, VZA and PHI
!    need also the height grids, earth radius and control

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs

   integer  , intent(In)    :: maxlayers, maxgeoms
   integer  , intent(In)    :: nlayers, ngeoms
   real(ffp), intent(InOut) :: alpha_boa(maxgeoms), theta_boa(maxgeoms), phi_boa(maxgeoms)
   logical  , intent(In)    :: do_Chapman
   real(ffp), intent(In)    :: dtr, eradius, heights (0:maxlayers), vsign

!  Los geometry

   real(ffp), intent(Out)  :: alpha         (0:maxlayers, maxgeoms)
   real(ffp), intent(Out)  :: radii         (0:maxlayers)
   real(ffp), intent(Out)  :: Raycon (maxgeoms)

!  main outputs (geometry)

   integer  , intent(Out)  :: ntraverse  (0:maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: sunpaths   (0:maxlayers,maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: chapfacs   (maxlayers,maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: cosscat(maxgeoms), term1(maxgeoms), term2(maxgeoms)

!  Local

   logical       :: Do_OverheadSun
   integer       :: n, k, v
   real(ffp)     :: alpha_boa_R, theta_boa_R
   real(ffp)     :: salpha_boa, calpha_boa
   real(ffp)     :: stheta_boa, ctheta_boa, cphi_boa, sunpaths_local(maxlayers)
   real(ffp), parameter   :: zero = 0.0_ffp
   real(ffp), parameter   :: one  = 1.0_ffp

!  Initialise output

   radii = zero ; alpha = zero ; Raycon = zero
   ntraverse = 0
   sunpaths  = zero ; chapfacs  = zero
   cosscat   = zero ; term1 = zero ; term2 = zero

!  Radii

   do n = 0, nlayers
      radii(n) = eradius + heights(n)
   enddo

!  Start geometry loop

   do v = 1, ngeoms

!  BOA angles

      alpha_boa_R    = alpha_boa(v) * DTR
      if ( alpha_boa(v).eq.90.0_ffp ) then
         calpha_boa     = zero
         salpha_boa     = one
      else
         calpha_boa     = cos(alpha_boa_R)
         salpha_boa     = sin(alpha_boa_R)
      endif

      theta_boa_R    = theta_boa(v) * DTR
      if ( theta_boa(v).eq.90.0_ffp ) then
         ctheta_boa     = zero
         stheta_boa     = one
      else
         stheta_boa     = sin(theta_boa_R)
         ctheta_boa     = cos(theta_boa_R)
      endif

      cphi_boa       = cos(phi_boa(v) * dtr)

!  Nominal traverse paths for Full illumination

      do n = 1, nlayers
         ntraverse(n,v) = n
      enddo

!  Overhead Sun

      Do_OverheadSun = (theta_boa(v).eq.zero) 

!  Set Alpha, ray constant, scattering angle

      alpha(1:nlayers,v) = alpha_boa_R
      Raycon           = salpha_boa * radii(nlayers)
      if ( Do_OverheadSun ) then
         term1(v) = zero
         term2(v) = calpha_boa
         cosscat(v) = - vsign * term2(v) ; if (term2(v).eq.zero) cosscat(v) = term2(v)
      else
         term1(v) = salpha_boa * stheta_boa * cphi_boa
         term2(v) = calpha_boa * ctheta_boa
         cosscat(v) = - vsign * term2(v) + term1(v)
      endif

!  Sunpath/Chapman factor calculations

!mick fix 4/12/12 - adjusted dimension of array "sunpaths" assignments to be
!                   compatible with subroutine "FindSunPaths_D" computations
      !sunpaths(0,1:maxlayers) = zero
      sunpaths(0,1:nlayers,v) = zero
      do n = 1, nlayers
         call FindSunPaths_D (Do_OverheadSun,Maxlayers,radii(n),Radii,&
           theta_boa_R,stheta_boa,N,sunpaths_local)
         !sunpaths(n,1:maxlayers) = sunpaths_local(1:maxlayers)
         sunpaths(n,1:n,v) = sunpaths_local(1:n)
         if ( do_Chapman ) then
            do k = 1, n
               chapfacs(n,k,v) = sunpaths(n,k,v)/(radii(k-1)-radii(k))
            enddo
         endif
      enddo

!  End geometry loop

   enddo
!  Finish

   return
end subroutine RegularPS_sphergeom

!  Finish

end module FO_geometry_SSonly_m

