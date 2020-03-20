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

module FO_geometry_Pool_m

!  Following routines are for Outgoing beams
!    They apply equally to Thermal and Solar-scatter Cases

! subroutine STD_outgoing_sphergeom_initial
! subroutine STD_outgoing_FindCritlayer
! subroutine STD_outgoing_sphergeom_Qbasic 
! subroutine STD_outgoing_sphergeom_Qadjusted 

!  Following routines are for Incoming Solar beams

! subroutine SD_incoming_FindCritLayer 
! subroutine SD_incoming_sphergeom 
! subroutine FindSun
! subroutine FindSunPaths_D
! subroutine FindSunPaths_T
! subroutine FindSunPath

!  Following is general purpose

!  Subroutine GAULEG_NG

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  EVERYTHING PUBLIC HERE
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

public

contains

subroutine STD_outgoing_sphergeom_Initial  &
       ( maxgeoms, maxlayers, ngeoms, nlayers, heights, eradius, alpha_boa, dtr, & ! Input
         doNadir, radii, Raycon, Lospaths, alpha, sina, cosa, cota )               ! Output

!  Completely stand-alone geometry routine for the outgoing STD correction
!     This is applicable to Both path geometries (up and down)
!     No Partial layer stuff here

!  Extension to Multiple Geometries, 29 October 2012

!  This routine: Initial LOS path setup

!    starting inputs are - BOA values of VZA (alpha_boa), in degrees
!                        - height grid, earth radius

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs
!  ------

!  Dimensions

   integer  , intent(in)   :: maxlayers, maxgeoms

!  number of geometries

   integer  , intent(in)   :: ngeoms

!  Layer control

   integer  , intent(in)   :: nlayers
   real(ffp), intent(in)   :: eradius, heights (0:maxlayers)

!  input angles

   real(ffp), intent(in)   :: alpha_boa(maxgeoms), dtr

!  Flag for the Nadir case

   logical  , intent(out)  :: doNadir(maxgeoms)
  
!  Alphas, Radii, Ray constant, Lospaths

   real(ffp), intent(out)  :: radii    (0:maxlayers)
   real(ffp), intent(out)  :: Raycon(maxgeoms)
   real(ffp), intent(out)  :: Lospaths (maxlayers,maxgeoms)
   real(ffp), intent(out)  :: alpha    (0:maxlayers,maxgeoms)
   real(ffp), intent(out)  :: cosa     (0:maxlayers,maxgeoms)
   real(ffp), intent(out)  :: sina     (0:maxlayers,maxgeoms)
   real(ffp), intent(out)  :: cota     (0:maxlayers,maxgeoms)

!  Local
!  -----

   integer      :: n, n1, v
   real(ffp)    :: salpha_boa, difh, alpha_boa_R
   real(ffp)    :: calpha, calpha1

   real(ffp), parameter :: zero = 0.0_ffp
   real(ffp), parameter :: one  = 1.0_ffp

!  Zero output

   cota = zero ; cosa = zero ; sina = zero ; alpha = zero
   donadir = .false. ; Raycon = zero ; Lospaths = zero

!  Radii

   do n = 0, nlayers
     radii(n) = eradius + heights(n)
   enddo

!  START LOOP
!  ==========

   do v = 1, ngeoms

!  Special case. Direct nadir viewing. Compute everything and Exit.

      if ( alpha_boa(v).eq.zero ) then
         doNadir(v) = .true.
         do n = nlayers,1,-1
            difh = radii(n-1) - radii(n) ; Lospaths(n,v) = difh
         enddo
         go to 67
      endif

!  Outgoing sphericity geometry (General case)
!  ===========================================

!  start at BOA

      alpha_boa_R    = alpha_boa(v) * DTR
      if ( alpha_boa(v) .eq. 90.0_ffp ) then
         salpha_boa     = one
         calpha1        = zero
      else
         salpha_boa     = sin(alpha_boa_R)
         calpha1        = cos(alpha_boa_R)
      endif

      cosa(nlayers,v)  = calpha1
      sina(nlayers,v)  = salpha_boa
      cota(nlayers,v)  = calpha1 / salpha_boa
      alpha(nlayers,v) = alpha_boa_R

!  Ray constant

      Raycon(v) = salpha_boa * radii(nlayers)

!  Whole layer values

      do n = nlayers - 1, 0, -1
         n1 = n + 1
         sina(n,v) = Raycon(v) / radii(n) ; alpha(n,v) = asin(sina(n,v))
         calpha  = cos(alpha(n,v)) ; cosa(n,v) = calpha 
         cota(n,v) = cosa(n,v) / sina(n,v)
         Lospaths(n1,v) = radii(n)*calpha - radii(n1)*calpha1
         calpha1 = calpha
      enddo

!  End loop

67    continue
   enddo

!  Finish

   return
end subroutine STD_outgoing_sphergeom_initial

!

SUBROUTINE STD_outgoing_FindCritlayer &
       ( maxgeoms, maxlayers, ngeoms, nlayers, Acrit, Cutoff, doNadir, & ! Inputs
         extinc, Lospaths, sina, cosa, radii, nfinedivs,               & ! Input
         Ncrit, AlphaCrit, RadCrit, CotCrit, fail, message )             ! Outputs

!  Purpose: Given a list of Maximum extinctions and LOS angles
!     Then find Critical layers (NCrit) and point where LOS attenuation wipe-outs (Acrit) are achieved
!     Then find the LOS angles and Radii (AlphaCrit,RadCrit) for these Critical Points

!  Extension to Multiple Geometries, 29 October 2012

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ------

!  Dimensioning

   integer, intent(in) :: maxlayers, maxgeoms

!  Layer and geometry numbers control

   integer, intent(in) :: nlayers, ngeoms

!  Attenuation and other parameters

   real(ffp), intent(in)  :: Acrit, Cutoff

!  Special case, Nadir viewing

   logical, intent(in)  :: doNadir(maxgeoms)

!  View angles and Radii at layer boundaries

   real(ffp), intent(in)  :: sina (0:maxlayers,maxgeoms)
   real(ffp), intent(in)  :: cosa (0:maxlayers,maxgeoms)
   real(ffp), intent(in)  :: radii(0:maxlayers)

!  Extinctions

   real(ffp), intent(in)  :: Lospaths(maxlayers,maxgeoms)
   real(ffp), intent(in)  :: extinc(maxlayers)

!  Modified inputs
!  ---------------

!  Number of Fine divisions

   integer, intent(inout) :: nfinedivs(maxlayers,maxgeoms)

!  outputs
!  -------

!  Critical layer, Number of Fine divisions for this layer

   integer  , intent(out)  :: Ncrit(maxgeoms)

!  Critical angle and radius and cotangent

   real(ffp), intent(out)  :: AlphaCrit(maxgeoms)
   real(ffp), intent(out)  :: RadCrit(maxgeoms), CotCrit(maxgeoms)

!  Exception handling

   logical      , intent(out) :: fail
   character*(*), intent(out) :: message

!  Local variables
!  ---------------

   real(ffp), parameter :: zero = 0.0_ffp
   real(ffp), parameter :: one  = 1.0_ffp

!  Other variables

   logical    ::  trawl
   integer    ::  N, N1, ntrans, v
   real(ffp)  ::  opdep, cumtrans_0, cumtrans_1,trans,transcrit,dcrit,TanCrit

!  Initialize

   Ncrit     = 0
   AlphaCrit = zero
   RadCrit   = zero
   CotCrit   = zero

   fail    = .false.
   message = ' '

!  Start Geometry loop

   do v = 1, ngeoms

!  Set trawl
!   Tested March 17th

      trawl = .true. ; n = 0 ; cumtrans_0 = one
      do while (trawl .and.n.lt.nlayers) 
         n = n + 1
         opdep = Lospaths(n,v) * extinc(n) ; trans = zero
         if ( opdep .lt. cutoff ) trans = exp ( - opdep )
         cumtrans_1 = cumtrans_0 * trans
         if ( cumtrans_1 .gt. Acrit ) then
            ntrans = int(-log10(trans) + 1)
            nfinedivs(n,v) = max(nfinedivs(n,v),ntrans)
            cumtrans_0 = cumtrans_1
         else
            NCrit(v) = n ; trawl = .false.
            transcrit    = Acrit / cumtrans_0
            ntrans = int(-log10(transcrit) + 1)
            nfinedivs(n,v) = max(nfinedivs(n,v),ntrans)
            dcrit        = - log(transcrit) / extinc(n)
            if ( doNadir(v) ) then
               Radcrit(v)    = radii(n-1) - dcrit      
            else
               n1 = n-1 ; TanCrit = radii(n1)*sina(n1,v)/(radii(n1)*cosa(n1,v)-dcrit)
               Cotcrit(v)    = one / Tancrit
               alphacrit(v)  = atan( TanCrit)
               radcrit(v)    = sina(n,v) * radii(n) / sin(alphacrit(v))    
            endif
         endif
      enddo

!  Zero the rest

      if ( NCrit(v) .ne. 0 ) nfinedivs(NCrit(v)+1:nlayers,v) = 0

!  End geometry loop

   enddo

!  Finish

   return
end subroutine STD_outgoing_FindCritlayer


subroutine SD_incoming_FindCritLayer &
       ( maxgeoms, maxlayers, ngeoms, nlayers, doNadir, dtr, Acrit, & ! Input
         cutoff, alpha, radii, extinc, Raycon, theta_boa,           & ! Inputs
         doCrit, Ncrit, nfinedivs, AlphaCrit, RadCrit, CotCrit,     & ! Outputs
         fail, message )                                              ! Outputs

!  Purpose: Given a list of Maximum extinctions and solar angles at BOA
!           Then find Critical layers (NCrit) and points where TOA attenuation wipe-outs (Acrit) are achieved
!           Then find the LOS angles and Radii (AlphaCrit,RadCrit) for these Critical Points
!           Nadir case, Alpha = 0.0, find only the radius (RadCrit)

!  Extension to Multiple Geometries, 29 October 2012

!  Find the Critical Radius (or angle) in layer Ncrit_i, Use Bisection 
!      based on the Function F(x) = opdep(x) - Crit_opdep

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ------

!  Dimensioning

   integer, intent(in) :: maxlayers, maxgeoms

!  Layer and geometry numbers control

   integer, intent(in) :: nlayers, ngeoms

!  Special case, Nadir viewing

   logical, intent(in)  :: doNadir(maxgeoms)

!  View angles and Radii at layer boundaries, Ray constant

   real(ffp), intent(in)  :: alpha(0:maxlayers,maxgeoms)
   real(ffp), intent(in)  :: radii(0:maxlayers)
   real(ffp), intent(in)  :: Raycon(maxgeoms)

!  Extinctions

   real(ffp), intent(in)  :: extinc(maxlayers)

!  Solar control and other parameters

   real(ffp), intent(in)  :: Acrit, theta_boa(maxgeoms), dtr, cutoff

!  Modified inputs (outputs)
!  -------------------------

!  Overall control (May be switched off if Critical test is negative for all Geometries)

   logical  , intent(inout)  :: DoCrit

!  Critical layer

   integer  , intent(inout)  :: Ncrit(maxgeoms)

!  Number of Fine divisions. This is updated according to Criticality

   integer  , intent(inout)  :: nfinedivs(maxlayers,maxgeoms)

!  Critical angle and radius and cotangent

   real(ffp), intent(inout)  :: AlphaCrit(maxgeoms)
   real(ffp), intent(inout)  :: RadCrit(maxgeoms), CotCrit(maxgeoms)

!  Outputs
!  -------

!  Exception handling

   logical      , intent(out) :: fail
   character*(*), intent(out) :: message

!  Local variables
!  ---------------

!  Bisection accuracy (lower for the Nadir case, as using distances)

   real(ffp), parameter  :: BisectionAccuracy_Nadir   = 1.0d-6
   real(ffp), parameter  :: BisectionAccuracy_General = 1.0d-9
   integer  , parameter  :: jmax = 50

   real(ffp), parameter  :: zero = 0.0_ffp
   real(ffp), parameter  :: one  = 1.0_ffp

!  Other variables

   character*2::  c2
   logical    ::  Finding, trawl, do_ZeroSunBOA, doCrit_local(maxgeoms)
   integer    ::  j, n, ntrans, ncrit_i, k, v
   real(ffp)  ::  s0, x1, x2, xmid, rtbis, dx, f, fmid, suncon, radii_x, sx, accuracy
   real(ffp)  ::  dist, theta_0, theta_1, stheta_0, stheta_1, ground_R, sunpaths(maxlayers)
   real(ffp)  ::  opdep, trans, atten_0, atten_1, transcrit, theta_n, stheta_n, theta_boa_R

!  Initialize

   fail    = .false.
   message = ' '

!  Start Geometry loop

   do v = 1, ngeoms

!  Initial setups

      doCrit_local(v) = .true.
      atten_0 = one ; theta_boa_R = theta_boa(v) * dtr
      if ( doNadir(v) ) then
         s0 = sin(theta_boa_R) ; ground_R = zero ; accuracy = BisectionAccuracy_Nadir
      else
         s0 = zero ; ground_R = alpha(nlayers,v) + theta_boa_R ; accuracy = BisectionAccuracy_General
      endif

!  Trawl through layers until Critical layer is reached. Nfinedivs is updated.
!    Only go down to Initial (LOS-path) Critical layer 
!       Condition on ZerosunBOA changed 27 March 2012

      NCrit_i = 0 ; trawl = .true. ; n = 0
      do while ( trawl .and.( (Ncrit(v).eq.0.and.n.lt.nlayers).or.(NCrit(v).ne.0.and.n.lt.NCrit(v)) ) )
         n = n + 1
         do_ZeroSunBOA = (n.eq.nlayers.and.theta_boa(v).eq.zero).or.(doNadir(v).and.theta_boa(v).eq.zero)
         if ( doNadir(v) ) then
            theta_n = theta_boa_R ; stheta_n = s0
         else
            theta_n = ground_R - alpha(n,v) ; stheta_n = sin(theta_n)
         endif
         call FindSunPaths_D (do_ZeroSunBOA,Maxlayers,radii(n),Radii,&
                              theta_n,stheta_n,n,sunpaths(:))
         opdep = sum(extinc(1:n)*sunpaths(1:n)) ; trans = zero
         atten_1 = zero; if ( opdep .lt. cutoff ) atten_1 = exp ( - opdep )
         if ( atten_1 .gt. Acrit ) then
            trans = atten_1 / atten_0
            ntrans = int(-log10(trans) + 1)
            nfinedivs(n,v) = max(nfinedivs(n,v),ntrans)
            atten_0 = atten_1
         else
            NCrit_i = n ; trawl = .false.
            transcrit    = Acrit / atten_0
            ntrans = int(-log10(transcrit) + 1)
            nfinedivs(n,v) = max(nfinedivs(n,v),ntrans)
         endif
      enddo

!  Nothing to do if No criticality (previous Critical values are unchanged)

      if ( trawl .and. NCrit_i.eq. 0 ) then
         if ( trawl .and. NCrit(v) .eq. 0 ) DoCrit_local(v) = .false.
         go to 67
      endif

!  Bisection: set Highest/Lowest value of Function (layer bottom/top). 

      if ( doNadir(v) ) then
         x1 = zero             ; x2 = radii(NCrit_i-1) - radii(NCrit_i)
      else
         x1 = alpha(NCrit_i-1,v) ; x2 = alpha(NCrit_i,v) 
      endif
      F     = log(Atten_0) - cutoff
      FMID  = opdep        - cutoff

!  Bisection: Check bracketing, if OK, perform bisection

      IF(F*FMID.GE.zero) then
         write(c2,'(I2)')v
         fail = .true. ; message = 'Root must be bracketed for bisection, geometry #'//c2 ; return
      ENDIF
      IF(F.LT.zero)THEN
         RTBIS=X1 ; DX=X2-X1
      ELSE
         RTBIS=X2 ; DX=X1-X2
      ENDIF

!  Bisection: Iterate to find the answer

      Finding = .true. ; J = 0
      DO While (Finding .and. j .lt. JMAX)
         J = J + 1 ; dx = 0.5_ffp * dx ; XMID = RTBIS + DX
         if ( doNadir(v) ) then
            theta_0 = theta_boa_R ; stheta_0 = s0
         radii_x = radii(NCrit_i-1) - xmid 
         else
            theta_0 = ground_R - xmid ; stheta_0 = sin(theta_0)
            radii_x = Raycon(v) / sin(xmid)
         endif
         suncon = radii_x * stheta_0
         stheta_1 = suncon / radii(NCrit_i-1) ;  theta_1 = asin(stheta_1)
         dist = radii(NCrit_i-1) * sin(theta_0-theta_1) / stheta_0
         opdep = dist * extinc(NCrit_i)
         theta_0 = theta_1 ; stheta_1 = stheta_0
         do k = n - 1, 1, -1
            stheta_1 = suncon / radii(k-1) ; theta_1 = asin(stheta_1)
            dist = radii(k-1) * sin(theta_0-theta_1) / stheta_0
            opdep = opdep + dist * extinc(k)
            theta_0 = theta_1 ; stheta_0 = stheta_1
         enddo
         fmid = opdep - cutoff
         IF ( FMID.LE.zero ) RTBIS = XMID
         IF(ABS(DX).LT.Accuracy .OR. FMID.EQ.0.) Finding = .false.
      ENDDO

!  Exception (too many bisections)

      if ( Finding ) Then
         write(c2,'(I2)')v
         fail = .true. ; message = 'Too many Bisections (540); Root not found, geometry #'//c2 ; return
      endif

!  Set final output if successful

      if ( doNadir(v) ) then
         RadCrit(v)   =  radii(NCrit_i-1) - RTBIS
      else
         AlphaCrit(v) = RTBIS ;  SX = sin(AlphaCrit(v))
         RadCrit(v)   = Raycon(v) / SX
         CotCrit(v)   = sqrt(one-SX*SX)/SX
      endif
      NCrit(v)     = NCrit_i

!  Continuation point

67    continue

!  Zero the rest

      if ( NCrit(v) .ne. 0 ) nfinedivs(NCrit(v)+1:nlayers,v) = 0

!  End geometry loop

   enddo

!  Reset criticality

   J = 0
   do v = 1, ngeoms
      if (doCrit_local(v) ) J = J + 1
   enddo
   doCrit = ( J.gt.0 )
 
!  Finish

   return
end subroutine SD_incoming_FindCritLayer

!

subroutine STD_outgoing_sphergeom_Qbasic                 &
       ( maxgeoms, maxlayers, maxfine, ngeoms, nlayers,  & ! Input
         nfinedivs, doNadir, radii, alpha, Raycon,       & ! Input
         radii_fine, alpha_fine, xfine, wfine,           & ! Output
         csqfine, cotfine )                                ! Output

!  Completely stand-alone geometry routine for the outgoing STD correction
!     This is applicable to Both path geometries (up and down)
!     No Partial layer stuff here

!  Extension to Multiple Geometries, 29 October 2012

!    starting inputs are - BOA values of VZA (alpha_boa), in degrees
!                        - height grid, earth radius, Layer control

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs
!  ------

!  Dimensions

   integer, intent(in)          :: maxlayers, maxfine, maxgeoms

!  Layer and geometry numbers control

   integer, intent(in) :: nlayers, ngeoms

!  Finelayer divisions is Strictly input

   integer, intent(in)       :: nfinedivs(maxlayers,maxgeoms)

!  Flag for the Nadir case

   logical  , intent(in)     :: doNadir(maxgeoms)
  
!  Alphas, Radii, Ray constant

   real(ffp), intent(in)  :: alpha      (0:maxlayers,maxgeoms)
   real(ffp), intent(in)  :: radii      (0:maxlayers)
   real(ffp), intent(in)  :: Raycon      (maxgeoms)

!  Outputs
!  =======

!  Fine layering

   real(ffp), intent(out)  :: alpha_fine (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(out)  :: radii_fine (maxlayers,maxfine,maxgeoms)

!  Quadratures

   real(ffp), intent(out)  :: xfine    (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(out)  :: wfine    (maxlayers,maxfine,maxgeoms)

!  Local geoemetry arrays

   real(ffp), intent(out)  :: csqfine  (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(out)  :: cotfine  (maxlayers,maxfine,maxgeoms)

!  Local
!  -----

   integer            :: n, n1, j, nfine, v
   real(ffp)          :: difh, csfine
   real(ffp)          :: tfine(maxfine), afine(maxfine)

   real(ffp), parameter :: zero = 0.0_ffp
   real(ffp), parameter :: one  = 1.0_ffp

!  Zero output

   alpha_fine = zero    ; radii_fine = zero
   xfine      = zero    ; wfine      = zero
   cotfine    = zero    ; csqfine    = zero

!  Start geometry loop
!  ===================

   do v = 1, ngeoms

!  Special case. Direct nadir viewing
!  ==================================

!  Compute everything and Exit. Qudratures are height-oriented
!    (This should be the same as the regular pseudo-spherical )

      if ( doNadir(v) ) then
         do n = nlayers,1,-1
            difh  = radii(n-1) - radii(n)
            nfine = nfinedivs(n,v)
            call gauleg_ng (zero,difh,tfine,afine,nfine,maxfine)
            do j = 1, nfine
               radii_fine(n,j,v) = radii(n-1) - tfine(j)
               xfine(n,j,v) = tfine(j)
               wfine(n,j,v) = afine(j)
            enddo
         enddo
         go to 67
      endif

!  Outgoing sphericity geometry (General case)
!  ===========================================

!  Whole layer values

      do n = nlayers, 1, -1
         n1 = n - 1
         nfine = nfinedivs(n,v)
         call gauleg_ng (alpha(n1,v),alpha(n,v),tfine,afine,nfine,maxfine)
         do j = 1,  nfine
            csfine = one / sin(tfine(j))
            radii_fine(n,j,v) = raycon(v) * csfine
            alpha_fine(n,j,v) = tfine(j)
            xfine(n,j,v)   = radii(n1) - radii_fine(n,j,v)
            wfine(n,j,v)   = afine(j)
            cotfine(n,j,v) = cos(tfine(j)) * csfine
            csqfine(n,j,v) = csfine * csfine
         enddo
      enddo

!  Continuation point

67    continue

!  End geometry loop

   enddo

!  Finish

   return
end subroutine STD_outgoing_sphergeom_Qbasic

!

subroutine STD_outgoing_sphergeom_Qadjusted &
       ( maxgeoms, maxlayers, maxfine, ngeoms, nlayers, nfinedivs,         & ! Input
         do_LOSpaths, doNadir, radii, alpha, Raycon,     & ! Input
         doCrit, Ncrit, AlphaCrit, RadCrit,              & ! Input
         radii_fine, alpha_fine, xfine, wfine,           & ! Output/Input
         csqfine, cotfine )                                ! Output/Input

!  Completely stand-alone geometry routine for the outgoing STD correction
!     This is applicable to Both path geometries (up and down)
!     No Partial layer stuff here

!    starting inputs are - BOA values of VZA (alpha_boa), in degrees
!                        - height grid, earth radius, Layer control
!                        - Critical Layer control

!  Regular Quadrature need not be done if LOSPATHS is set

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs
!  ------

!  Dimensions

   integer, intent(in)          :: maxlayers, maxfine, maxgeoms

!  Layer and geometry numbers control

   integer, intent(in) :: nlayers, ngeoms

!  Finelayer divisions may be changed

   integer, intent(inout)       :: nfinedivs(maxlayers,maxgeoms)

!  Flag for the Nadir case

   logical  , intent(in)     :: doNadir(maxgeoms)
  
!  Flag for pre-existing calculations

   logical  , intent(in)     :: do_LOSpaths
  
!  Alphas, Radii, Ray constant

   real(ffp), intent(in)  :: alpha      (0:maxlayers,maxgeoms)
   real(ffp), intent(in)  :: radii      (0:maxlayers)
   real(ffp), intent(in)  :: Raycon(maxgeoms)

!  Critical stuff

   logical  , intent(in)  :: doCrit
   integer  , intent(in)  :: Ncrit(maxgeoms)
   real(ffp), intent(in)  :: AlphaCrit(maxgeoms)
   real(ffp), intent(in)  :: RadCrit(maxgeoms)

!  Outputs
!  =======

!  Fine layering

   real(ffp), intent(out)  :: alpha_fine (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(out)  :: radii_fine (maxlayers,maxfine,maxgeoms)

!  Quadratures

   real(ffp), intent(out)  :: xfine    (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(out)  :: wfine    (maxlayers,maxfine,maxgeoms)

!  Local geoemetry arrays

   real(ffp), intent(out)  :: csqfine  (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(out)  :: cotfine  (maxlayers,maxfine,maxgeoms)

!  Local
!  -----

   integer            :: n, n1, j, nfine, v
   real(ffp)          :: difh, csfine
   real(ffp)          :: tfine(maxfine), afine(maxfine)

   real(ffp), parameter  :: zero = 0.0_ffp
   real(ffp), parameter  :: one  = 1.0_ffp

!  Zero output

   if ( .not. do_LOSpaths ) then
      alpha_fine = zero    ; radii_fine = zero
      xfine      = zero    ; wfine      = zero
      cotfine    = zero    ; csqfine    = zero
   endif

!  Start geometry loop
!  ===================

   do v = 1, ngeoms

!  Special case. Direct nadir viewing
!  ==================================

!  Compute everything and Exit. Qudratures are height-oriented
!    (This should be the same as the regular pseudo-spherical )

      if ( doNadir(v) ) then

!  For normal atmosphere, Regular quadrature only if flagged

         if ( .not. doCrit .or. NCrit(v) .eq. 0 ) then
            if ( .not. do_LOSpaths ) then
               do n = nlayers,1,-1
                  difh  = radii(n-1) - radii(n)
                  nfine = nfinedivs(n,v)
                  call gauleg_ng (zero,difh,tfine,afine,nfine,maxfine)
                  do j = 1, nfine
                     radii_fine(n,j,v) = radii(n-1) - tfine(j)
                     xfine(n,j,v) = tfine(j)
                     wfine(n,j,v) = afine(j)
                  enddo
               enddo
            endif

!  Otherwise.....

         else

!    -- Adjust quadrature for the Critical layer
            
            n = NCrit(v) ; difh = radii(n-1) - Radcrit(v) ; nfine = nfinedivs(n,v)
            call gauleg_ng (zero,difh,tfine,afine,nfine,maxfine)
            do j = 1, nfine
               radii_fine(n,j,v) = radii(n-1) - tfine(j)
               xfine(n,j,v) = tfine(j)
               wfine(n,j,v) = afine(j)
            enddo

!    -- For all layers above Critical layer, Regular quadrature only if flagged

            if ( .not. do_LOSpaths ) then
               do n = NCrit(v)-1,1,-1
                  difh  = radii(n-1) - radii(n) ; nfine = nfinedivs(n,v)
                  call gauleg_ng (zero,difh,tfine,afine,nfine,maxfine)
                  do j = 1, nfine
                     radii_fine(n,j,v) = radii(n-1) - tfine(j)
                     xfine(n,j,v) = tfine(j)
                     wfine(n,j,v) = afine(j)
                  enddo
               enddo
            endif
         endif

!  Done Nadir case so continue

         go to 67
      endif

!  Outgoing sphericity geometry (General case)
!  ===========================================

!  For normal atmosphere, Regular quadrature only if flagged

      if ( .not. doCrit .or. NCrit(v) .eq. 0 ) then
         if ( .not. do_LOSpaths ) then
            do n = nlayers, 1, -1
               n1 = n - 1
               nfine = nfinedivs(n,v)
               call gauleg_ng (alpha(n1,v),alpha(n,v),tfine,afine,nfine,maxfine)
               do j = 1,  nfine
                  csfine = one / sin(tfine(j))
                  radii_fine(n,j,v) = raycon(v) * csfine
                  alpha_fine(n,j,v) = tfine(j)
                  xfine(n,j,v)   = radii(n1) - radii_fine(n,j,v)
                  wfine(n,j,v)   = afine(j)
                  cotfine(n,j,v) = cos(tfine(j)) * csfine
                  csqfine(n,j,v) = csfine * csfine
               enddo
            enddo
         endif

!  Otherwise

      else

!    -- Adjust quadrature for the Critical layer

         n = NCrit(v) ; n1 = n - 1 ; nfine = nfinedivs(n,v)
         call gauleg_ng (alpha(n1,v),AlphaCrit(v),tfine,afine,nfine,maxfine)
         do j = 1,  nfine
            csfine = one / sin(tfine(j))
            radii_fine(n,j,v) = raycon(v) * csfine
            alpha_fine(n,j,v) = tfine(j)
            xfine(n,j,v)   = radii(n1) - radii_fine(n,j,v)
            wfine(n,j,v)   = afine(j)
            cotfine(n,j,v) = cos(tfine(j)) * csfine
            csqfine(n,j,v) = csfine * csfine
         enddo

!    -- For all layers above Critical layer, Regular quadrature only if flagged

         if ( .not. do_LOSpaths ) then
            do n = NCrit(v)-1,1,-1
               n1 = n - 1
               nfine = nfinedivs(n,v)
               call gauleg_ng (alpha(n1,v),alpha(n,v),tfine,afine,nfine,maxfine)
               do j = 1,  nfine
                  csfine = one / sin(tfine(j))
                  radii_fine(n,j,v) = raycon(v) * csfine
                  alpha_fine(n,j,v) = tfine(j)
                  xfine(n,j,v)   = radii(n1) - radii_fine(n,j,v)
                  wfine(n,j,v)   = afine(j)
                  cotfine(n,j,v) = cos(tfine(j)) * csfine
                  csqfine(n,j,v) = csfine * csfine
               enddo
            enddo
         endif

!  Done

      endif

!  Continuation point

67    continue

!  End geometry loop

   enddo

!  Finish

   return
end subroutine STD_outgoing_sphergeom_Qadjusted

!

subroutine SD_incoming_sphergeom &
       ( maxgeoms, maxlayers, maxfine, ngeoms, nlayers, nfinedivs, do_Chapman, doNadir, & ! Input
         DoCrit, NCrit, alpha_boa, theta_boa, phi_boa, radii, alpha,  & ! Input
         vsign, dtr, Pie, RadCrit, AlphaCrit, radii_fine, alpha_fine, & ! Input
         sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,          & ! Output
         chapfacs, cosscat, theta_all, phi_all )                        ! Output

!  Completely stand-alone geometry routine for Accurate SS
!     This is for the incoming Solar Beams
!     This is applicable to Both Upwelling and Downwelling LOS-path geometries
!     No partials, this routine

!    starting inputs are the BOA values of SZA, VZA and PHI
!    need also the height grids, earth radius and control
!    need also the complete values of all VZAs along outgoing paths

!  This routine has the fine gridding treatment

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs

   integer  , intent(In)    :: maxlayers, maxfine, maxgeoms
   integer  , intent(In)    :: ngeoms, nlayers, nfinedivs(maxlayers,maxgeoms), NCrit(maxgeoms)
   real(ffp), intent(InOut) :: alpha_boa(maxgeoms), theta_boa(maxgeoms), phi_boa(maxgeoms)
   real(ffp), intent(In)    :: vsign
   logical  , intent(In)    :: do_Chapman, DoCrit, doNadir(maxgeoms)

!  Los geometry

   real(ffp), intent(In)   :: alpha         (0:maxlayers,maxgeoms)
   real(ffp), intent(In)   :: alpha_fine    (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(In)   :: radii         (0:maxlayers)
   real(ffp), intent(In)   :: radii_fine    (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(In)   :: AlphaCrit(maxgeoms), RadCrit(maxgeoms), dtr, pie

!  main outputs (geometry)

   integer  , intent(Out)  :: ntraverse  (0:maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: sunpaths   (0:maxlayers,maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: chapfacs   (maxlayers,maxlayers,maxgeoms)

!  Fine level output (geometry)

   integer  , intent(Out)  :: ntraverse_fine(maxlayers,maxfine,maxgeoms)
   real(ffp), intent(Out)  :: sunpaths_fine (maxlayers,maxlayers,maxfine,maxgeoms)

!  scattering angle and associated angles

   real(ffp), intent(Out)  :: cosscat     (maxgeoms)
   real(ffp), intent(Out)  :: theta_all  (0:maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: phi_all    (0:maxlayers,maxgeoms)

!  Local

   logical       :: DirectSun, Do_OverheadSun, Do_ZeroSunBOA, Do_Normal
   integer       :: n, j, k, v
   real(ffp)     :: SolarDirection(3), Radstart, term1, term2
   real(ffp)     :: salpha_boa, calpha_boa, phi_boa_R, sphi_boa
   real(ffp)     :: theta_boa_R, stheta_boa, ctheta_boa, cphi_boa
   real(ffp)     :: ctheta, stheta, calpha, salpha, cphi, CumAngle

   real(ffp), parameter  :: zero = 0.0_ffp
   real(ffp), parameter  :: one  = 1.0_ffp

!  Check
!   real(ffp)  :: sumd, sume, sth1

!  Local arrays associated with fine grid output

   logical         :: DirectSunf(maxfine)
   real(ffp)       :: thetaf(maxfine)
   real(ffp)       :: sthetaf(maxfine)
   real(ffp)       :: cthetaf(maxfine)

!  Initialise output

   ntraverse = 0     ; ntraverse_fine = 0
   sunpaths  = zero ; sunpaths_fine  = zero
   chapfacs  = zero

   phi_all = zero   ; theta_all = zero ; cosscat = zero

!  Start geometry loop

   do v = 1, ngeoms

!  Nominal traverse paths for Full illumination

      ntraverse(0,v) = 0
      do n = 1, nlayers
         ntraverse(n,v) = n
         do j = 1, nfinedivs(n,v)
            ntraverse_fine(n,j,v) = n
         enddo
      enddo

!  check range of inputs, already done.....

!  Special case

      Do_OverheadSun = theta_boa(v).eq.zero

!  BOA angles

      if ( alpha_boa(v).eq.90.0_ffp ) then
         calpha_boa     = zero
         salpha_boa     = one
      else
         salpha_boa  = sin(alpha(nlayers,v))
         calpha_boa  = cos(alpha(nlayers,v))
      endif

      theta_boa_R    = theta_boa(v) * DTR
      if ( theta_boa(v).eq.90.0_ffp ) then
         ctheta_boa     = zero
         stheta_boa     = one
      else
         stheta_boa     = sin(theta_boa_R)
         ctheta_boa     = cos(theta_boa_R)
      endif

      phi_boa_R   = phi_boa(v) * dtr
      cphi_boa    = cos(phi_boa_R)
      sphi_boa    = sin(phi_boa_R)

!  define Unit solar vector at BOA

      if ( Do_OverheadSun ) then
         SolarDirection = 0.0_ffp
      else
         SolarDirection(1) = - stheta_boa * cphi_boa * vsign
         SolarDirection(2) = - stheta_boa * sphi_boa
         SolarDirection(3) = - ctheta_boa
      endif

!  Cosine of scattering angle at boa

      if ( Do_OverheadSun ) then
         term1 = zero
         term2 = calpha_boa
         cosscat(v) = - vsign * term2 ; if (term2.eq.zero) cosscat(v) = term2
      else
         term1 = salpha_boa * stheta_boa * cphi_boa
         term2 = calpha_boa * ctheta_boa
         cosscat(V) = - vsign * term2 + term1 
      endif

!  General case: LOS path in spherical geometry
!  ============================================

!  Start loop over all layers

      do n = nlayers, 1, -1

!  Special cases

        DO_ZeroSunBOA  = Do_OverheadSun.and.(n.eq.nlayers.or.doNadir(v))
        DO_Normal      = .not. doCrit .or. ( doCrit .and. n.le. NCrit(v) )

!  Layer boundary Sun position
!     * Local save of angles, cosines, sines and  illumination flags
!     * Use critical ALPHA and RADIUS if N = NCrit
!     * Use Bottom-of-layer values if N < NCrit (BOA values if illuminated)

!        if ( do_Normal ) then                               !   @@RTSFix 9/5/12 (Comment out line)
           if ( doCrit .and. n .eq. NCrit(v) ) then
              CumAngle = alpha(nlayers,v) - AlphaCrit(v) ; Radstart = RadCrit(v)
              call FindSun(DoNadir(v),Do_OverheadSun,Radstart,SolarDirection,CumAngle,theta_boa_R,&
                           theta_all(n,v),stheta,ctheta,DirectSun)
           else
              Radstart = radii(n)
              if ( n.eq. nlayers ) then
                 theta_all(n,v) = theta_boa(v)*dtr ; stheta = stheta_boa ; ctheta = ctheta_boa ; DirectSun = .true.
              else
                 CumAngle = alpha(nlayers,v) - alpha(n,v)
                 call FindSun(DoNadir(v),Do_OverheadSun,radii(n),SolarDirection,CumAngle,theta_boa_R,&
                              theta_all(n,v),stheta,ctheta,DirectSun)
              endif
           endif
!        endif                                               !   @@RTSFix 9/5/12 (Comment out line)

!  Fine-layer sun positions

        if ( Do_Normal ) then
           do j = 1, nfinedivs(n,v)
              CumAngle = alpha(nlayers,v) - alpha_fine(n,j,v)
              call FindSun(DoNadir(v),Do_OverheadSun,radii_fine(n,j,v),SolarDirection,CumAngle,theta_boa_R,&
                           thetaf(j),sthetaf(j),cthetaf(j),DirectSunf(j))
           enddo
        endif

!  Sun paths in layer

!        if ( do_Normal ) then                               !   @@RTSFix 9/5/12 (Comment out line)
           if ( DirectSun ) then
              call FindSunPaths_D (Do_ZeroSunBOA,Maxlayers,Radstart,Radii,&
                theta_all(n,v),stheta,N,sunpaths(n,:,v))
           else
              call FindSunPaths_T (Maxlayers,Pie,Radstart,Radii,theta_all(n,v),stheta,N,sunpaths(n,:,v),ntraverse(n,v))
           endif
        if ( Do_Normal ) then                                !   @@RTSFix 9/5/12 (Addline)
           do j = 1, nfinedivs(n,v) 
              if ( DirectSunf(j) ) then
                 call FindSunPaths_D &
                  (Do_ZeroSunBOA,Maxlayers,Radii_fine(n,j,v),Radii,&
                   thetaf(j),sthetaf(j),N,sunpaths_fine(n,:,J,v))
              else
                 call FindSunPaths_T &
                  (Maxlayers,Pie,Radii_fine(n,j,v),Radii,thetaf(j),sthetaf(j),N,sunpaths_fine(n,:,J,v),ntraverse_fine(n,J,v))
              endif
!             if ( n.eq.14 ) write(*,*)j,n,Radii_fine(n,j)-radii(n)
           enddo
        endif

!  debugging

!        if ( n.eq.14) then
!       sumd = SUM(sunpaths(n,1:ntraverse(n),v))
!       sth1 = stheta*RadCrit(v)/radii(0)
!       sume = sin(theta_all(n,v) - asin(sth1))*radii(0)/stheta
!       write(*,*)n,sumd,sume
!       do j = 1, nfinedivs(n)
!         sumd = SUM(sunpaths_fine(n,1:ntraverse_fine(n,j,v),j,v))
!         sth1 = sthetaf(j)*radii_fine(n,j,v)/radii(0)
!         sume = sin(thetaf(j) - asin(sth1))*radii(0)/sthetaf(j)
!         write(*,*)j,sumd,sume
!       enddo
!       pause
!      endif

!  Fix phi by using constancy of scatter angle
!     If AZM > 180, Subtract from 360 for consistency. (VLIDORT code, 10 October 2011)

        if (Do_OverheadSun.or.doNadir(v) ) then
           phi_all(n,v)     = phi_boa(v) * dtr
        else
!           if ( do_Normal ) then                               !   @@RTSFix 9/5/12 (Comment out line)
              if ( doCrit .and. n .eq. NCrit(v) ) then
                 salpha = sin(AlphaCrit(v))
                 calpha = cos(AlphaCrit(v))
              else
                 salpha = sin(alpha(n,v))
                 calpha = cos(alpha(n,v))
              endif
              cphi = (cosscat(v)+vsign*calpha*ctheta)/stheta/salpha
              if ( cphi.gt.one)  cphi = one
              if ( cphi.lt.-one) cphi = -one
              phi_all(n,v)     = acos(cphi)
              if ( phi_boa(v).gt.180.0_ffp) phi_all(n,v) = 2.0_ffp * Pie - phi_all(n,v)
!           endif                                               !   @@RTSFix 9/5/12 (Comment out line)
        endif

!  End layer loop

      enddo

!  TOA Sun angle sunpaths and PHI.
!    (No sunpaths if directly illuminated)

      DO_ZeroSunBOA  = Do_OverheadSun.and.doNadir(v)
      CumAngle = alpha(nlayers,v) - alpha(0,v) ; Radstart = radii(0)
      call FindSun(DoNadir(v),Do_OverheadSun,Radstart,SolarDirection,CumAngle,theta_boa_R,&
                   theta_all(0,v),stheta,ctheta,DirectSun)
      if (.not.DirectSun ) then
          call FindSunPaths_T (Maxlayers,Pie,Radii(0),Radii,theta_all(0,v),stheta,1,sunpaths(0,:,v),ntraverse(0,v))
      endif
      if ( Do_OverheadSun .or. doNadir(v) ) then
         phi_all(0,v)     = phi_boa(v) * dtr
      else
         cphi = (cosscat(v)+vsign*calpha*ctheta)/stheta/salpha
         if ( cphi.gt.one)  cphi = one ; if ( cphi.lt.-one) cphi = -one
         phi_all(0,v)     = acos(cphi)
         if ( phi_boa(v).gt.180.0_ffp) phi_all(0,v) = 2.0_ffp * Pie - phi_all(0,v)
      endif

!  Chapman factor calculations
!  ---------------------------

      if ( do_Chapman ) then
         do n = 1, nlayers
            call FindSunPaths_D (Do_OverheadSun,Maxlayers,radii(n),Radii,&
              theta_boa_R,stheta_boa,N,chapfacs(n,:,v))
            do k = 1, n
               chapfacs(n,k,v) = chapfacs(n,k,v)/(radii(k-1)-radii(k))
            enddo
         enddo
      endif

!  End geometry loop

   enddo

!  Finish

   return
end subroutine SD_incoming_sphergeom

!

!  General Routines for Sun positioning

subroutine FindSun(DoNadir,Do_OverheadSun,Radius,SolarDirection,CumAngle,theta_boa_R,theta,stheta,ctheta,DirSun)

!  Find the solar anlge along the LOS path, for given radius and cumulative angle from BOA
!    SolarDirection is defined at BOA, with azimuth relative to the LOS direction.

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs

   logical   , Intent(In)    :: DoNadir,Do_OverheadSun
   real(ffp) , Intent(in)    :: Radius,SolarDirection(3),CumAngle,theta_boa_R

!  Outputs

   real(ffp) , Intent(out)   :: theta,stheta,ctheta
   logical   , Intent(InOut) :: DirSun

!  Local

   real(ffp) :: px(3),b
   real(ffp), parameter :: zero = 0.0_ffp
   real(ffp), parameter :: one  = 1.0_ffp

!  Calculation (Nadir view scenario)

   if ( doNadir ) then
      DirSun = .true.
      theta = theta_boa_R
      ctheta = cos(theta_boa_R)
      stheta = sin(theta_boa_R)
      return
   endif

!  Calculation (overhead sun)

   if ( Do_OverheadSun ) then
      DirSun = .true.
      ctheta = cos(CumAngle)
      stheta = sin(CumAngle)
      theta  = CumAngle
      return
   endif

!  Calculation (General)

   px(1) = - Radius * sin(CumAngle)
   px(2) = zero
   px(3) =   Radius * cos(CumAngle)
   b = DOT_PRODUCT(px,SolarDirection)
   ctheta = -b/Radius
   DirSun = ( ctheta.ge.zero )
   stheta = sqrt(one-ctheta*ctheta)
   theta  = acos(ctheta)

!  Done

   return
end subroutine FindSun


subroutine FindSunPaths_D (Do_ZeroSunBOA,Maxlayers,Radstart,Radii,&
                           thstart,sthstart,N,sunpaths)

!  Sunpaths for the Direct-sun illumination
!  Starting point is Radstart on the LOS path, with solar angle thstart, in layer N
!  Special case = Overhead sun at BOA

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs

   LOGICAL   , Intent(In)   :: Do_ZeroSunBOA
   INTEGER   , Intent(In)   :: maxlayers, N
   real(ffp) , Intent(In)   :: Radstart,Radii(0:maxlayers)
   real(ffp) , Intent(In)   :: thstart,sthstart

!  Output

   real(ffp), Intent(InOut) :: Sunpaths(maxlayers)

!  Local

   integer    :: n1, k
   real(ffp)  :: sth0, th0, sth1, th1, ks1
   real(ffp), parameter :: zero = 0.0_ffp
   real(ffp), parameter :: one  = 1.0_ffp

!  Layer boundary upper

   N1 = N - 1

!  SBOA condition

   if ( Do_ZeroSunBOA ) then
      sunpaths(n) = radii(n1) - Radstart
      do k = n1, 1, -1
         sunpaths(k) = radii(k-1) - radii(k)
      enddo
      return
   endif

!  First layer

   sth0 = sthstart
   th0  = thstart
   sth1 = sth0*Radstart/radii(N1)
   th1  = asin(sth1)
   ks1  = th0-th1
   sunpaths(n) = sin(ks1)*Radstart/sth1

!  Other layers to TOA

   sth0 = sth1
   th0  = th1
   do k = n1, 1, -1
      sth1 = sth0*radii(k)/radii(k-1)
      th1  = asin(sth1)
      ks1  = th0-th1
      sunpaths(k) = sin(ks1)*radii(k)/sth1
      sth0 = sth1
      th0  = th1
   enddo

!  Done

   return
end subroutine FindSunPaths_D

subroutine FindSunPaths_T (Maxlayers,Pie,Radstart,Radii,thstart,sthstart,N,sunpaths,NT)

!  Sunpaths for the Tangent-height illumination
!  Starting point is Radstart on the LOS path, with solar angle thstart, in layer N

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs

   INTEGER   , Intent(In)   :: maxlayers, N
   real(ffp) , Intent(In)   :: Radstart,Radii(0:maxlayers)
   real(ffp) , Intent(In)   :: thstart,sthstart,Pie

!  Output

   INTEGER   , Intent(InOut)   :: NT
   real(ffp), Intent(InOut)    :: Sunpaths(maxlayers)

!  Local

   logical    :: trawl
   integer    :: n1, k
   real(ffp)  :: sth0, th0, sth1, th1, ks1, tanr

   real(ffp), parameter :: zero = 0.0_ffp
   real(ffp), parameter :: one  = 1.0_ffp
   real(ffp), parameter :: two  = 2.0_ffp

!  Layer boundary upper

   N1 = N - 1

!  tangent height, Find which layer NT

   NT = N
   tanr = sthstart * Radstart
   k = n1 ; trawl = .true.
   do while (k.ge.n1.and.trawl)
      trawl = (radii(k).gt.tanr) ; k = k + 1
   enddo
   nt = k-1 !; write(*,*)n,nt

!  Distances for layers N and below to NT

   if ( nt.gt.n ) then
      th0  = pie - thstart ; sth0 = sthstart
      sth1 = sth0*Radstart/radii(n)
      th1  = asin(sth1) ; ks1  = th0-th1
      sunpaths(n) = two * sin(ks1)*Radstart/sth1
      sth0 = sth1 ; th0 = th1
      do k = n+1,nt-1
        sth1 = sth0*radii(k-1)/radii(k)
        th1  = asin(sth1) ; ks1  = th0-th1
        sunpaths(k) = two * sin(ks1)*radii(k)/sth0
        sth0 = sth1 ; th0 = th1
      enddo
      sth1 = one ; ks1 = 0.5_ffp * pie - th0
      sunpaths(nt) = two * sin(ks1)*radii(nt-1)
   else if ( nt.eq.n ) then
      sunpaths(n) = - two * Radstart * cos(thstart)
   endif

!  Rest of layer n up to the upper boundary

   th0 = pie - thstart ; sth0 = sthstart
   sth1 = sth0*Radstart/radii(N1)
   th1  = asin(sth1) ; ks1  = th0-th1
   sunpaths(n) = sunpaths(n) + sin(ks1)*Radstart/sth1
   sth0 = sth1 ; th0 = th1

!  Trawl up from layers above n, to TOA

   do k = n1, 1, -1
      sth1 = sth0*radii(k)/radii(k-1)
      th1  = asin(sth1)
      ks1  = th0-th1
      sunpaths(k) = sin(ks1)*radii(k)/sth1 
      sth0 = sth1
      th0  = th1
   enddo

!  Done

   return
end subroutine FindSunPaths_T

SUBROUTINE GAULEG_NG(X1,X2,X,W,N,NMAX)

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Input/Output

      INTEGER  , intent(in)  :: N,NMAX
      REAL(ffp), intent(in)  :: X1, X2
      REAL(ffp), intent(out) :: X(NMAX),W(NMAX)

      INTEGER     :: I, M, J
      REAL(ffp)   :: XM,XL,P1,P2,P3,PP,Z,Z1
      REAL(ffp)   :: QUANT, RN, RJ, PIE, ARG

      REAL(ffp), PARAMETER :: EPS = 3.0D-14
      real(ffp), parameter :: zero = 0.0_ffp
      real(ffp), parameter :: one  = 1.0_ffp
      real(ffp), parameter :: two  = 2.0_ffp
      real(ffp), parameter :: half  = 0.5_ffp
      real(ffp), parameter :: qtr   = 0.25_ffp

      M=(N+1)/2
      XM = half * (X2+X1)
      XL = half * (X2-X1)
      RN = real(N,ffp)
      Z1 = zero
      pie = acos(-one)

      DO I=1,M
            arg = ( real(i,ffp) - qtr ) / ( rn + half )
!            Z=COS(3.141592654D0*(I-.25D0)/(N+.5D0))
            Z  = COS ( pie * arg )
            DO WHILE (ABS(Z-Z1).GT.EPS)
                  P1=one
                  P2=zero
                  DO J=1,N
                        RJ = real(J,ffp)
                        P3=P2
                        P2=P1
                        P1= ( ( two*RJ-one)*Z*P2-(RJ-one)*P3 ) / RJ
                        P1= ( ( two*RJ-one)*Z*P2-(RJ-one)*P3 ) / RJ
                  ENDDO
                  PP=RN*(Z*P1-P2)/(Z*Z-one)
                  Z1=Z
                  Z=Z1-P1/PP
            ENDDO
            QUANT = two / ( ( one-Z*Z) * PP *PP )
            X(I)     = XM - XL*Z
            X(N+1-I) = XM + XL*Z
            W(I)     = QUANT * XL
            W(N+1-I) = W(I)
      ENDDO
      RETURN
END SUBROUTINE GAULEG_NG


subroutine FindSunPath  ( x, xboa, rtop, raycon, sundir, sundist, theta0 )

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  I/O

   real(ffp), intent(in)  ::  x, xboa, rtop, raycon, sundir(3)
   real(ffp), intent(out) ::  sundist, theta0

!  Local

   real(ffp) :: xicum, sinx, rad, c0, s0, s1, c1
   real(ffp), parameter  :: zero = 0.0_ffp
   real(ffp), parameter  :: one  = 1.0_ffp

!  Subroutine for the quick calculation of sunpath from point X on View Path to point at Top of layer

   xicum = xboa - x
   sinx  = sin(x)
   rad   = Raycon / sinx
   c0 = sundir(1) * sin(xicum) - sundir(3) * cos(xicum)
   theta0 = - acos(c0)
   s0 = sqrt(one-c0*c0)
   s1 = s0 * rad / rtop
   c1 = sqrt(one-s1*s1)
   sundist = -rtop * (s1*c0-s0*c1)/s0

!  finish

   return
end subroutine FindSunPath

!  Finish

end module FO_geometry_Pool_m

