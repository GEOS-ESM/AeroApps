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

module FO_geometry_DTonly_m

   use FO_geometry_Pool_m, only : STD_outgoing_sphergeom_initial, STD_outgoing_FindCritlayer, &
                                  STD_outgoing_sphergeom_Qbasic , STD_outgoing_sphergeom_Qadjusted

!  Routines for Geometry with Thermal only

!   Extension to multiple geometries, 29 October 2012   (Version 1.3)

public DT_Geometry

contains

subroutine DT_Geometry  &
       ( maxgeoms, maxlayers, maxfine, ngeoms, nlayers, nfinedivs, dtr, eradius, heights, alpha_boa, & ! Input
         do_LOSpaths, do_planpar, do_regular_ps, do_enhanced_ps, doNadir, doCrit,  & ! Inpu
         Acrit, extinc, Raycon, radii, alpha, cota, xfine, wfine,                  & ! Input/Output
         csqfine, cotfine, alphafine, radiifine, NCrit, RadCrit, CotCrit,          & ! Output
         fail, message, trace )                                                      ! Output

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Input arguments
!  ===============

!  Dimensions

   integer  , intent(in)     :: maxlayers, maxfine, maxgeoms

!  Layer and geometry control. Finelayer divisions may be changed

   integer  , intent(in)     :: ngeoms, nlayers
   integer  , intent(inout)  :: nfinedivs(maxlayers,maxgeoms)

!  Radius + heights

   real(ffp), intent(in)     :: eradius, heights (0:maxlayers)

!  input angle (Degrees), dtr = degrees-to-Radians

   real(ffp), intent(in)     :: dtr
   real(ffp), intent(InOut)  :: alpha_boa(maxgeoms)

!    LOS paths flag is new, 20 January 2012

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

!  Exception handling

   logical      , intent(out)    :: fail
   character*(*), intent(out)    :: message
   character*(*), intent(out)    :: trace

!  Local arguments
!  ===============

!  LOS path lengths

   real(ffp)  :: Lospaths (maxlayers,maxgeoms)

!  Other angles

   real(ffp)  :: cosa       (0:maxlayers,maxgeoms)
   real(ffp)  :: sina       (0:maxlayers,maxgeoms)

!  Critical values

   real(ffp)  :: AlphaCrit(maxgeoms)

!  Help variables

   real(ffp), parameter   :: zero = 0.0_ffp
   real(ffp)              :: cutoff, alpha_boa_R
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
      trace   = 'Initial Flag Check in DT_Geometry'
      fail    = .true. ;  return
   endif

!  Cannot have Plane-parallel and Enhanced PS

   if ( do_planpar .and. do_enhanced_ps ) then
      message = 'Cannot have BOTH Plane-parallel and Enhanced PS options'
      trace   = 'Initial Flag Check in DT_Geometry'
      fail    = .true. ;  return
   endif

!  Cannot have Regular PS and Enhanced PS

   if ( do_regular_ps .and. do_enhanced_ps ) then
      message = 'Cannot have BOTH Regular and Enhanced PS options'
      trace   = 'Initial Flag Check in DT_Geometry'
      fail    = .true. ;  return
   endif

!  VZA can be 0-90 degrees inclusive, but not outside this range

   do v = 1, ngeoms
      if ( alpha_boa(v).gt.90.0_ffp.or.alpha_boa(v).lt.zero ) then
         write(c2,'(I2)')v
         message = 'Boa LOS angle outside range [0,90]); Check it!'
         trace   = 'Geometry # '//c2//'; Initial Angle Check in DT_Geometry_1'
         fail    = .true. ;  return
      endif
    enddo

!  Plane-parallel or regular-PS, just set alpha and return
!  -------------------------------------------------------

   if ( do_planpar .or. do_regular_ps ) then
      alpha = zero
      do v = 1, ngeoms
         alpha_boa_R      = alpha_boa(v) * DTR
         alpha(1:nlayers,v) = alpha_boa_R
      enddo
      return
   endif

!  Enhanced PS; proceed in 3 Steps
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
      CALL STD_outgoing_FindCritlayer  &
       ( maxgeoms, maxlayers, ngeoms, nlayers, Acrit, Cutoff, doNadir, & ! Inputs
         extinc, Lospaths, sina, cosa, radii, nfinedivs,               & ! Input
         Ncrit, AlphaCrit, RadCrit, CotCrit, fail, message )             ! Outputs
      if ( Fail ) then
         trace = 'Error from STD_Outgoing_FindCritLayer in DT_Geometry' ; return
      endif
   endif

!  Step 3. Set Quadratures
!  -----------------------

!  Step 3a. LOS fine-layer quadratures (Regular, non-adjusted, no Criticality)
!           Regular quadratures may be avoided if do_LOSpaths is set

   if ( .not.doCrit) then
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
       ( maxgeoms, maxlayers, maxfine, ngeoms, nlayers, nfinedivs, & ! Input
         do_LOSpaths, doNadir, radii, alpha, Raycon,               & ! Input
         doCrit, Ncrit, AlphaCrit, RadCrit,                        & ! Input
         radiifine, alphafine, xfine, wfine,                       & ! Output
         csqfine, cotfine )                                          ! Output
   endif

!  Finish

   return
end subroutine DT_Geometry

!  Finish

end module FO_geometry_DTonly_m

