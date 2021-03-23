
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
! #  Authors :     Robert. J. D. Spurr (1)                      #
! #                Matt Christi                                 #
! #                                                             #
! #  Address (1) : RT Solutions, inc.                           #
! #                9 Channing Street                            #
! #                Cambridge, MA 02138, USA                     #
! #                                                             #
! #  Tel:          (617) 492 1183                               #
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  This Version :   VLIDORT_2p8p2                             #
! #  Release Date :   15 April 2020                             #
! #                                                             #
! #  Previous VLIDORT Versions under Standard GPL 3.0:          #
! #  ------------------------------------------------           #
! #                                                             #
! #      2.7   F90, released August 2014                        #
! #      2.8   F90, released May    2017                        #
! #      2.8.1 F90, released August 2019                        # 
! #                                                             #
! #  Features Summary of Recent VLIDORT Versions:               #
! #  -------------------------------------------                #
! #                                                             #
! #      NEW: TOTAL COLUMN JACOBIANS         (2.4)              #
! #      NEW: BPDF Land-surface KERNELS      (2.4R)             #
! #      NEW: Thermal Emission Treatment     (2.4RT)            #
! #      Consolidated BRDF treatment         (2.4RTC)           #
! #      f77/f90 Release                     (2.5)              #
! #      External SS / New I/O Structures    (2.6)              #
! #                                                             #
! #      SURFACE-LEAVING / BRDF-SCALING      (2.7)              #
! #      TAYLOR Series / OMP THREADSAFE      (2.7)              #
! #      New Water-Leaving Treatment         (2.8)              #
! #      LBBF & BRDF-Telescoping, enabled    (2.8)              #
! #      Several Performance Enhancements    (2.8)              #
! #      Water-leaving coupled code          (2.8.1)            #
! #      Planetary problem, media properties (2.8.1)            #
! #                                                             #
! #  Features Summary of This VLIDORT Version                   #
! #  ----------------------------------------                   #
! #                                                             #
! #   2.8.2, released 15 April 2020.                            #
! #     ==> Geometry (FO/MS), check/derive separation           #
! #     ==> New setup_master for Geometry/Check/Derive          #
! #     ==> Reduction of zeroing, some dynamic memory           #
! #     ==> Use of F-matrixes only in FO code                   #
! #     ==> Use I/O type structures directly                    #
! #     ==> Doublet geometry post-processing option             #
! #                                                             #
! ###############################################################

! ###################################################################
! #                                                                 #
! # This is Version 2.8.2 of the VLIDORT_2p8 software library.      #
! # This library comes with the Standard GNU General Public License,#
! # Version 3.0, 29 June 2007. Please read this license carefully.  #
! #                                                                 #
! #      VLIDORT Copyright (c) 2003-2020.                           #
! #          Robert Spurr, RT Solutions, Inc.                       #
! #          9 Channing Street, Cambridge, MA 02138, USA.           #
! #                                                                 #
! # This file is part of VLIDORT_2p8p2 ( Version 2.8.2 )            #
! #                                                                 #
! # VLIDORT_2p8p2 is free software: you can redistribute it         #
! # and/or modify it under the terms of the Standard GNU GPL        #
! # (General Public License) as published by the Free Software      #
! # Foundation, either version 3.0 of the License, or any           #
! # later version.                                                  #
! #                                                                 #
! # VLIDORT_2p8p2 is distributed in the hope that it will be        #
! # useful, but WITHOUT ANY WARRANTY; without even the implied      #
! # warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR         #
! # PURPOSE. See the Standard GNU General Public License (GPL)      #
! # for more details.                                               #
! #                                                                 #
! # You should have received a copy of the Standard GNU General     #
! # Public License (GPL) Version 3.0, along with the VLIDORT_2p8p2  #
! # code package. If not, see <http://www.gnu.org/licenses/>.       #
! #                                                                 #
! ###################################################################

! ###############################################################
! #                                                             #
! # WaterLeaving Subroutines in this Module                     #
! #                                                             #
! #         Linearized_WaterLeaving (Top-level)                 #
! #           - Lin_Ocean_Reflectance (ex. MORCASIWAT)          #
! #           - Interpolate_Lin_fOQ_BS1                         #
! #           - Interpolate_Lin_fOQ_BS2                         #
! #           - Lin_WhiteCap_Reflectance                        #
! #           - Lin_Water_Transmittance                         #
! #               * Lin_GENERAL_SUNGLINT                        #
! #                                                             #
! ###############################################################

      MODULE vsleave_lin_sup_routines_2_m

!  Version 2.6 Notes
!  -----------------

!  INDWAT, MORCASIWAT Routines taken straight from Clark Weaver code
!      Compiled here by R. Spurr, 11 July 2012
!  get_fluorescence_755 Routines taken straight from Chris O'dell code
!      Compiled here by R. Spurr, 12 July 2012

!  Upgrade for Version 2.7
!  -----------------------

!  New Water-leaving code developed/tested by R. Spurr, April 2014
!  Based in part on Modified-6S code by A. Sayer.
!  Validated against Modified-6 OCEABRDF.F code, 24  April 2014.

!  Water-leaving upgrade according to the modified 6S specification
!  New code calculates transmittances into and out of ocean, using
!  usual sun-glint rough surface approximations. In addition, the
!  water-leaving term itself is now SZA-dependent (A. Sayer), and
!  there is now a correction for Whitecaps (again, from 6S)

!  This upgrade gives the  water-leaving terms some directionality,
!  but they are still azimuth-independent

!  Earlier version inputs were just Wavelength/Salinity/PigmentConc
!  This was enough for the isotropic case (Fast Option) in Version 2.6.
!  For Version 2.7, we require additional inputs, including:
!    - Wind-speed and Direction (direction was not used in earlier version) 
!    - flags to control use sunglint shadowing and foam (whitecaps) correction.

!  This water-leaving option is designed to work alongside the "NewCM" glint
!  reflectance option in the VBRDF code. The glint and whitecap calculations
!  in the two supplements are the same.

!  You need to make sure that the wind input information is the same as that
!  used for the "NewCM" glint option in the VBRDF supplement. Also, the Foam
!  correction applied here in the surface-leaving code should also be 
!  applied for "NewCM" glint.

!  R. Spurr, 06 October 0215. 
!    -- Expanded list of Private routines, including foQ and alternative Ocean reflectance
!    -- Introduced the Rough_Surface parameter into main argument list.

!  Upgrade for Version 2.8
!  -----------------------

!  R. Spurr 4/12/18.
!   --  Fixed the "MU0 factor" bug in the Water Leaving calculation.
!   --  Fixed the indexing bug in Interpolate_foQ_2 (Anisotropic case)

      use vsleave_sup_aux_m
      use vsleave_sup_routines_2_m, ONLY :  Water_RefracIndex,         &
                                            Water_Transmittance_Quads, &
                                            Fresnel_Complex,           &
                                            TaSav_Interpolate

      PRIVATE :: Lin_WhiteCap_Reflectance,    &
                 Lin_Water_Transmittance,     &
                 Interpolate_Lin_fOQ_BS1,     &
                 Interpolate_Lin_fOQ_BS2,     &
                 Lin_Ocean_Reflectance_Basic, &
                 Lin_GENERAL_SUNGLINT

      PUBLIC  :: Linearized_WaterLeaving_2

      CONTAINS

subroutine Linearized_WaterLeaving_2 &
     ( Maxszas, Maxvzas, Maxstreams,                                       &
       FoQFile, TaRayFile, do_Approximate_Ta, do_Isotropy,                 &
       Do_Rough_Surface, Do_FoamOption, Do_GlintShadow, Do_FacetIsotropy,  &
       Wavelength, Salinity, PigmentConc, Windspeed, WindSunAngles,        &
       nszas, nvzas, nstreams, szas, vzas, streams,                        &
       WLeaving_ISO, WLeaving_SD, WLeaving_SV,                             &
       dC_WLeaving_ISO, dC_WLeaving_SD, dC_WLeaving_SV,                    &
       dW_WLeaving_ISO, dW_WLeaving_SD, dW_WLeaving_SV,                    & 
       TaSav, fail, message )                                           

! mick fix 12/28/2014 - Using normalization correction by A Sayer, 04 Nov 2014
! Apply a normalisation factor of 1/(pi/mu0) to output water-leaving reflectance, to
! bring things in line with results from e.g. 6S simulations and expected behaviour.
! Think this is a subtlety related to reflectance vs. normalised radiance treatment,
! although it is very obvious if you don't do it. Correction applied at end of the
! subroutine.

! Rob  Fix 10/06/15  Applying comments by ! A Sayer, 22 Sep 2015

! Above normalisation factor was not correct. Instead of being 1/(pi/mu0), it should
! have been Ta/Q, where Ta is downwelling atmospheric transmittance and Q is a
! correction term. Note the previous normalisation worked in most cases since Ta can
! be close to mu0 and Q can be close to pi. However this updated version is better.
! R. Spurr will provide exact Ta for use in the future.

! Rob  Fix 10/06/15  Applying comments by ! A Sayer 30 Sep 2015

! Q factor has been moved into OceanReflecs part. Rather than have f and Q separately,
! as before, we will now use a lookup table of the f/Q ratio provided by David Antoine.

! Rob  Fix 10/06/15  Applying comments by ! R. Spurr, 01 October 2015

! Use of approximate formula for Ta (downwelling atmospheric transmittance) is now
! controlled by an internal parameter flag (do_Approximate_Ta). Also, the f/q lookup
! table and interpolation thereby has been implemented in a first version using 
! averaged fOQ tables for each wavelength, pigment value and solar zenith angle.
! With this first attempt, there is no need to change the routine I/O.
! The introduction of a more-sophisticated interpolation scheme for the whole
! f/Q LUT requires a major overhaul of the entire subroutine dependencies, as
! we should then have to introduce VZA and AZM inputs....

   IMPLICIT NONE
   INTEGER  , parameter:: fpk = SELECTED_REAL_KIND(15)

!  Newly constructed Water-Leaving supplement routine. Version 2.7 VLIDORT code.
!    23-24 April 2014. R. Spurr, RT Solutions Inc.
!    Validated against the Modified 6S code by A. Sayer.. 

!  This is a Stand-alone subroutine.

!  inputs
!  ======

!  Dimensioning

   integer  , intent(in)   :: Maxszas, Maxvzas, Maxstreams

!  Files

   character*(*), intent(in) :: FoQFile, TaRayFile

!  Approximate Ta flag (New 10/5/15)

   logical, intent(in) :: do_Approximate_Ta

!  Logical flags
!  -------------

!  Isotropic (Fast Calculation) option, assumes all transmittances = 1

   Logical  , intent(in)   :: do_Isotropy

!  Rough surface option is now separate

   Logical  , intent(in)   :: do_Rough_Surface

!  These 3 flags all apply to the Rough Surface option
!     Optional inclusion of Foam term
!     Optional inclusion of Shadow term for Glitter (Air-water only?)
!     Flag for using Isotropic Facet distribution

   Logical  , intent(in)   :: Do_FoamOption
   Logical  , intent(in)   :: Do_GlintShadow
   LOGICAL  , intent(in)   :: Do_FacetIsotropy

!  Physical
!  --------

!  Wavelength in Micrometers

   real(fpk), intent(in)   :: Wavelength

!  Salinity

   real(fpk), intent(in)   :: Salinity

!  Pigment concentration

   real(fpk), intent(in)   :: PigmentConc

!  Windspeed m/s, wind-Sun Azimuth angles in Radians
!    -- Rough Surface options only

   REAL(fpk), intent(in)   :: WINDSPEED
   REAL(fpk), intent(in)   :: WindSunAngles ( Maxszas )

!  Sun, viewing and stream angles
!  ------------------------------

   integer  , intent(in) :: nszas, nvzas, nstreams
   real(fpk), intent(in) :: szas   (Maxszas)
   real(fpk), intent(in) :: vzas   (Maxvzas)
   real(fpk), intent(in) :: streams(Maxstreams)

!  OUTPUT
!  ======

!  Isotropic value. Fast calculation

   REAL(fpk), intent(out)    :: WLeaving_ISO    ( Maxszas )
   REAL(fpk), intent(out)    :: dC_WLeaving_ISO ( Maxszas )
   REAL(fpk), intent(out)    :: dW_WLeaving_ISO ( Maxszas )

!  Input solar, output stream angles

   REAL(fpk), intent(out)    :: WLeaving_SD    ( Maxszas, Maxstreams )
   REAL(fpk), intent(out)    :: dC_WLeaving_SD ( Maxszas, Maxstreams )
   REAL(fpk), intent(out)    :: dW_WLeaving_SD ( Maxszas, Maxstreams )

!  input solar, output view angles

   REAL(fpk), intent(out)    :: WLeaving_SV    ( Maxszas, Maxvzas )
!mick fix 9/19/2017 - changed extent of 2nd dim
   !REAL(fpk), intent(out)    :: dC_WLeaving_SV ( Maxszas, Maxstreams )
   !REAL(fpk), intent(out)    :: dW_WLeaving_SV ( Maxszas, Maxstreams )
   REAL(fpk), intent(out)    :: dC_WLeaving_SV ( Maxszas, Maxvzas )
   REAL(fpk), intent(out)    :: dW_WLeaving_SV ( Maxszas, Maxvzas )

!  Atmospheric Transmittance (Diagnostic output)

   REAL(fpk), intent(out)    :: TaSav ( Maxszas, 4 )

!  Exception handling

   logical      , intent(out)  :: fail
   character*(*), intent(out)  :: message
 
!  remark: Still no Azimuth dependence here.....

!  LOCAL
!  =====

!  Variables for the Rough Surface Option
!  --------------------------------------

!  Transmittance Quadratures

   integer, parameter   :: Max_PolarQuads=24, Max_AzimQuads=48
   real(fpk) :: PolarQuads    (Max_PolarQuads)   ! Radians
   real(fpk) :: CosPolarQuads (Max_PolarQuads)
   real(fpk) :: SinPolarQuads (Max_PolarQuads)
   real(fpk) :: PolarWeights  (Max_PolarQuads)

   real(fpk) :: AzimQuads    (Max_AzimQuads)     ! Radians
   real(fpk) :: CosAzimQuads (Max_AzimQuads)
   real(fpk) :: SinAzimQuads (Max_AzimQuads)
   real(fpk) :: AzimWeights  (Max_AzimQuads)

!  Glitter control
!  ---------------

   logical   :: do_coeffs, local_do_shadow
   REAL(fpk) :: SUNGLINT_COEFFS(7), dSUNGLINT_COEFFS(7), Refrac_R, Refrac_I, Refrac_sq
   real(fpk) :: phi_w, cphi_w, sphi_w, local_refrac_R, local_refrac_I

!  other variables
!  ---------------

!  Help

   logical   :: noWL
   integer   :: J, I, k
   real(fpk) :: dtr, pi, Albedo, Const0, Const1, Const_Iso, denom, Rwprime, eta, f
   real(fpk) :: incident_angle, Local_Sine, SZA_cosines(Maxszas), mu0
   real(fpk) :: dW_Const0, dW_Const1, dC_Rwprime, dW_Const_Iso, dC_Const_Iso, dC_eta, dC_f
   real(fpk) :: Foam_correction, WC_Reflectance, WC_Lambertian
   real(fpk) :: dW_Foam_correction, dWC_Reflectance, dWC_Lambertian
   real(fpk), parameter :: zero = 0.0_fpk
   real(fpk), parameter :: one  = 1.0_fpk

!  Ocean reflectances
!    -- New directional arrays, R. Spurr 06 October 2015

   real(fpk) :: Ocean_Reflec_Basic                      ! new
   real(fpk) :: Ocean_Iso_Reflecs(Maxszas)              ! renamed
   real(fpk) :: Ocean_SV_Reflecs (Maxszas, Maxvzas)     ! New
   real(fpk) :: Ocean_SD_Reflecs (Maxszas, Maxstreams)  ! New
   real(fpk) :: dC_Ocean_Reflec_Basic                      ! new
   real(fpk) :: dC_Ocean_Iso_Reflecs(Maxszas)              ! renamed
   real(fpk) :: dC_Ocean_SV_Reflecs (Maxszas, Maxvzas)     ! New
   real(fpk) :: dC_Ocean_SD_Reflecs (Maxszas, Maxstreams)  ! New
   real(fpk), parameter :: Minimum_Reflectance = 0.0001_fpk

! f/Q calculation.  R. Spurr 03-06 October 2015 (with lienarizations)

   real(fpk) :: foQ_Int_1  ( maxszas )
   real(fpk) :: foQ_Int_SV ( maxszas, Maxvzas )
   real(fpk) :: foQ_Int_SD ( maxszas, Maxstreams )
   real(fpk) :: dC_foQ_Int_1  ( maxszas )
   real(fpk) :: dC_foQ_Int_SV ( maxszas, Maxvzas )
   real(fpk) :: dC_foQ_Int_SD ( maxszas, Maxstreams )

!  Atmospheric downwelling transmittance

!  Rob Fic 10/06/15, Appllying earlier work to this module

!  A Sayer 22 Sep 2015. Downwelling transmittance, Q factor, approximate Rayleigh optical depth.

!  R. Spurr, 10/01/15. Use of an approximate form of diffuse-field transmittance "Ta"
!                      Based on Gordon and Wang formula (1994). Thanks to A. Sayer.
!                      If Not set here, then Ta will default to 1.0, and should then
!                      be calculated inside of VLIDORT by a dedicated routine.
!  R. Spurr, 10/05/15. Instead of Ta defaulting to 1.0, interpolate from a dedicated set
!                      of Rayleigh-atmosphere transmittances (diagnostic output).

   integer      :: nTa_szas, nTa_wavs
   real(fpk)    :: Ta, tau_rayleigh, dum
   real(fpk)    :: Ta_szas(14), Ta_wavs(64), TaSavData(64,14,4)

!  rough surface transmittance help variables

   logical   :: do_RS_transmittances
   real(fpk) :: Trans_Norm, Trans_solar(Maxszas), dW_Trans_solar(Maxszas)
   real(fpk) :: Trans_stream(Maxstreams), Trans_Viewing(Maxvzas)
   real(fpk) :: dW_Trans_stream(Maxstreams), dW_Trans_Viewing(Maxvzas)

!  Initial Setup
!  -------------

!  conversions

   pi = acos(-1.0_fpk)
   dtr = pi / 180.0_fpk
   SZA_cosines = zero
   do J = 1, nszas
      SZA_cosines(J) = cos(szas(J)*dtr)
   enddo

!  Zero the output

   fail = .false. ; message = ' '
   WLeaving_ISO = zero ; dC_WLeaving_ISO = zero ; dW_WLeaving_ISO = zero
   WLeaving_SD  = zero ; dC_WLeaving_SD  = zero ; dW_WLeaving_SD  = zero
   WLeaving_SV  = zero ; dC_WLeaving_SV  = zero ; dW_WLeaving_SV  = zero
   TaSav        = zero

!  Refractive index. Formerly INDWAT

   Call  Water_RefracIndex  ( Wavelength, Salinity, Refrac_R, Refrac_I )
   Refrac_sq = Refrac_R * Refrac_R + Refrac_I * Refrac_I

!  Ocean Leaving
!  -------------

! R. Spurr 01 October 2015

!  First Version of the interpolation of David Antoine's lookup table of foQ
!  Interpolates averaged values (average over all 221 VZA/AZM combinations)
!  Methodology --
!    (1) Linear with wavelength, Cos(SZA), Log(Pigment)
!    (2) No Extrapolation. If out-of-range, use end-values

! R. Spurr 03 October 2015. Revised to include NON-ISOTROPIC option
!    -- This arises because of the directionality of foQ.
!    -- Original Isotropic routine has been preserved ( Formerly MORCASIWAT )
!    -- New directional routine has additiona stream/viewing angles input, + directional output
!    -- re-engineered to privde first the basic reflectance, then FoQ factoring

!  First make the basic call. This will then need to be multiplied by an foQ factor

   call  Lin_Ocean_Reflectance_Basic &
     ( Wavelength, PigmentConc, noWL, Ocean_Reflec_Basic, eta, dC_Ocean_Reflec_Basic, dC_eta )

!  Ocean Leaving ; Exit if no contribution (outside range 200-900 nm)

   if ( noWL ) return

!  Now Derive the FoQ factor from database and interpolation
!  ---------------------------------------------------------

!  Isotropy case, FoQ factor from database + 221-averaged interpolation

   if ( do_Isotropy ) then
      call Interpolate_Lin_fOQ_BS1 &
           ( FoQFile, Maxszas, nszas, szas, Wavelength, PigmentConc, &
             foQ_Int_1, dC_foQ_Int_1, fail, message )
      if ( fail ) return
      do j = 1, nszas
         Ocean_Iso_Reflecs(J)    = Ocean_Reflec_Basic * foQ_Int_1(J)
         dC_Ocean_Iso_Reflecs(J) = dC_Ocean_Reflec_Basic *    foQ_Int_1(J) &
                                    + Ocean_Reflec_Basic * dC_foQ_Int_1(J)
      enddo
   endif

!  Directional case, FoQ factor from database + 13-averaged interpolation

   if ( .not. do_Isotropy ) then

      call Interpolate_Lin_fOQ_BS2 &
           ( FoQFile, Maxszas, Maxstreams, Maxvzas, nszas, nvzas, nstreams, &
             szas, vzas, streams, refrac_R, Wavelength, PigmentConc,        &
             foQ_Int_1, foQ_Int_SV, foQ_Int_SD,                             &
             dC_foQ_Int_1, dC_foQ_Int_SV, dC_foQ_Int_SD, fail, message ) 
      if ( fail ) return
      do j = 1, nszas
         Ocean_Iso_Reflecs(J)    = Ocean_Reflec_Basic * foQ_Int_1(J)
         dC_Ocean_Iso_Reflecs(J) = dC_Ocean_Reflec_Basic *    foQ_Int_1(J) &
                                    + Ocean_Reflec_Basic * dC_foQ_Int_1(J)
         do i = 1, nstreams 
            Ocean_SD_Reflecs(j,i)    = Ocean_Reflec_Basic * foQ_Int_SD(J,I)
            dC_Ocean_SD_Reflecs(j,i) = dC_Ocean_Reflec_Basic *    foQ_Int_SD(J,I) &
                                        + Ocean_Reflec_Basic * dC_foQ_Int_SD(J,I)
         enddo
         do i = 1, nvzas
            Ocean_SV_Reflecs(j,i)    = Ocean_Reflec_Basic * foQ_Int_SV(J,I)
            dC_Ocean_SV_Reflecs(j,i) = dC_Ocean_Reflec_Basic *    foQ_Int_SV(J,I) &
                                        + Ocean_Reflec_Basic * dC_foQ_Int_SV(J,I)
         enddo
      enddo
   endif

!  Here is the original code again, just for reminder !!!!!!!!!!!!!
!! Morel and Gentili, 1991.
!     do j = 1, nszas
!! Now we have taken the Q factor out of main WaterLeaving routine,
!! replace f calculation here with the f/Q ratio from Morel et al (AO, 2002).
!!        mu_sol = real(SZA_cosines(J))
!!        f=0.6279 - (0.2227*eta) - (0.0513*eta*eta) + (0.2465*eta -0.3119 )*mu_sol
!!        R2=f*(b_tot/a_tot)
!! R. Spurr 01 October 2015
!!        foQ = 0.09 ! Placeholder until we have implementation of David Antoine's lookup table.
!!        Removed placeholder and inserted interpolated values
!        foQ=foQ_Int_1(J) 
!!        write(77,'(2f8.3,f10.6)')  wl*1000.0d0, Log(C),foQ
!        R2=foQ*(b_tot/a_tot)
!        Ocean_Reflecs(J) = real(R2,fpk)
!     enddo
!   write(*,*) trans_norm ; pause'after quads'

!  Prepare for Trans_Atmos output
!  ------------------------------

! R Spurr 01-05 October 2015

      ! Introduced the flag to control computation of this approximation
      ! More Exact calculation inside VLIDORT now programmed and tested
      ! Read a pre-calculated Rayleigh-atmopshere data set. Somewhat of a fiddle

!  If approximate formula, use estimated Tau_Rayleigh
!        tau_rayleigh approximate calculation is from Hansen and Travis (1984).
!  If use database, read the table for interpolation

      if ( Do_Approximate_Ta ) then
         tau_rayleigh=0.008569*(Wavelength**(-4))*(1+0.0113*(Wavelength**(-2)) +0.00013*(Wavelength**(-4)))
      else
!mick fix 9/19/2017 - included error handling 'err=' clause
         open(76,file=Trim(TaRayFile),err=88,status='old',action='read')
         read(76,*) ; read(76,*) nTa_szas, nTa_wavs
         if ( nTa_szas.ne.14 .or. nTa_wavs .ne. 64 ) stop 'TaRaydata not valid for file-read'
         read(76,*) ; read(76,*) Ta_szas(1:nTa_szas)
         read(76,*) ; read(76,*) Ta_wavs(1:nTa_wavs) ; read(76,*)
         do i = 1, nTa_wavs
            read(76,'(f8.2,4(2x,14e14.6))')dum,((TaSavData(i,j,k),j = 1,nTa_szas),k=1,4)
         enddo
         close(76)
      endif

!  Introduce Rough surface control
!  -------------------------------

!  R Spurr 03 Oct 2015 Introduce Rough surface control
!     -- rename do_transmittances --> do_RS_Transmittances

!  Initialize Transmittances to 1.0.
!   These will be the flat surface defaults.

   trans_norm = one
   Trans_solar   = one ; dW_Trans_solar   = zero
   Trans_viewing = one ; dW_Trans_viewing = zero
   Trans_stream  = one ; dW_Trans_stream  = zero

!  Foam-reflectance correction. Rough Surface Only

   Foam_correction = one ; dW_Foam_correction = zero
   if ( do_rough_surface .and. Do_FoamOption ) then
      call Lin_WhiteCap_Reflectance &
         ( WindSpeed, Wavelength, WC_Reflectance, WC_Lambertian, &
                        DWC_Reflectance, DWC_Lambertian )
      Foam_correction    = one - WC_Reflectance
      dW_Foam_correction =     - dWC_Reflectance
   endif

!  set  Coeffs flag, initialize local shadow flag

   do_coeffs       = .true.
   local_do_shadow = do_GlintShadow

!  Set Rough Surface Transmittances flag. Not required for the Fast calculation option

   do_RS_transmittances = .false.
   if ( .not. Do_Isotropy.and.do_Rough_Surface ) then
!mick mod 9/19/2017 - replaced hardcoded min value with "Minimum_Reflectance"
      do J = 1, nszas
         !if ( Ocean_Iso_Reflecs(J).gt.0.0001 ) do_RS_transmittances =.true.
         if ( Ocean_Iso_Reflecs(J).gt.Minimum_Reflectance ) do_RS_transmittances =.true.
      enddo
   endif

!  Get quadratures

   if ( do_RS_transmittances ) then
      call Water_Transmittance_Quads &
       ( Max_PolarQuads, Max_AzimQuads,                          & ! Input
         PolarQuads, CosPolarQuads, SinPolarQuads, PolarWeights, & ! Output
         AzimQuads,  CosAzimQuads,  SinAzimQuads,  AzimWeights,  & ! Output
         TRANS_NORM )
   endif

!  Solar angle incident loop
!  =========================

   do J = 1, nszas

!  Rob Fix 10/03/15. Proper control for the transmittance calculations
!        (All values are initialized to 1.0)

!  Rough surface transmittances
!  ----------------------------

      if ( do_RS_transmittances ) then

!  Downward Solar transmittances. Only perform if the Ocean_Iso_Reflecs term is non-trivial
!       Passing from Air to Water, Local shadow flag is required.
!       Minimum_Reflectance value set as a parameter, formerly 0.0001

         if ( Ocean_Iso_Reflecs(J).gt.Minimum_Reflectance ) then
            phi_w = WindSunAngles(J)
            cphi_w = cos(phi_w*dtr)
            sphi_w = sin(phi_w*dtr)
            local_do_shadow = do_GlintShadow
            call Lin_Water_Transmittance &
               ( Max_PolarQuads, Max_AzimQuads,                          & ! Input
                 PolarQuads, CosPolarQuads, SinPolarQuads, PolarWeights, & ! Input
                 AzimQuads,  CosAzimQuads,  SinAzimQuads,  AzimWeights,  & ! Input
                 do_FacetIsotropy, LOCAL_DO_SHADOW, DO_COEFFS,           & ! Input
                 szas(j), REFRAC_R, REFRAC_I,                            & ! Input
                 WINDSPEED, SUNGLINT_COEFFS, dSUNGLINT_COEFFS,           & ! Input
                 PHI_W, CPHI_W, SPHI_W, TRANS_NORM,                      & ! Inpu
                 TRANS_SOLAR(J), dW_TRANS_SOLAR(J) )
         endif

!  Upward transmittances into View directions
!     Passing from water to air, use Snell's Law.  no absorption
!     Local shadow flag turned off here.

         local_do_shadow  = .false.
         local_refrac_R   = one / refrac_R
         local_refrac_I   = zero
         if ( Ocean_Iso_Reflecs(J).gt.Minimum_Reflectance ) then
            do i = 1, nvzas
               incident_angle = asin(sin(vzas(i) * dtr)/refrac_R)/dtr
               call Lin_Water_Transmittance &
                  ( Max_PolarQuads, Max_AzimQuads,                          & ! Input
                    PolarQuads, CosPolarQuads, SinPolarQuads, PolarWeights, & ! Input
                    AzimQuads,  CosAzimQuads,  SinAzimQuads,  AzimWeights,  & ! Input
                    do_FacetIsotropy, LOCAL_DO_SHADOW, DO_COEFFS,           & ! Input
                    incident_angle, local_refrac_R, local_refrac_I,         & ! Input
                    WINDSPEED, SUNGLINT_COEFFS, dSUNGLINT_COEFFS,           & ! Input
                    PHI_W, CPHI_W, SPHI_W, TRANS_NORM,                      & ! Input
                    TRANS_VIEWING(I), dW_TRANS_VIEWING(I) )
            enddo
         endif

!  Upward transmittances into stream directions
!     Passing from water to air, use Snell's Law.  no absorption

         local_do_shadow  = .false.
         local_refrac_R   = one / refrac_R
         local_refrac_I   = zero
         if ( Ocean_Iso_Reflecs(J).gt.Minimum_Reflectance ) then
            do i = 1, nstreams
               local_sine = sqrt ( one - streams(i) * streams(i) )
               incident_angle = asin(local_sine/refrac_R)/dtr
               call Lin_Water_Transmittance &
                  ( Max_PolarQuads, Max_AzimQuads,                          & ! Input
                    PolarQuads, CosPolarQuads, SinPolarQuads, PolarWeights, & ! Input
                    AzimQuads,  CosAzimQuads,  SinAzimQuads,  AzimWeights,  & ! Input
                    do_FacetIsotropy, LOCAL_DO_SHADOW, DO_COEFFS,           & ! Input
                    incident_angle, local_refrac_R, local_refrac_I,         & ! Input
                    WINDSPEED, SUNGLINT_COEFFS, dSUNGLINT_COEFFS,           & ! Input
                    PHI_W, CPHI_W, SPHI_W, TRANS_NORM,                      & ! Input
                    TRANS_STREAM(I), dW_TRANS_STREAM(I) )
            enddo
         endif

!  End clause for RS transmittances

      endif

!  R Spurr 03 Oct 2015. Revised Final Computation
!  ----------------------------------------------

!     Trans_solar, Trans_Stream, Trans_Viewing are all 1.0 for the flat surface (default)

      Albedo  = 0.485_fpk

!  Isotropy term
!   -- Constructed as before, using the Ocean_Iso_Reflecs
!      Rwprime = Ocean_Iso_Reflecs(J) / ( one - Albedo * Ocean_Iso_Reflecs(J)  )
!      Const_Iso = ( one / refrac_sq ) * trans_solar(J) * Rwprime
!      Const_Iso = Const_Iso * Foam_correction
!      WLeaving_Iso(J) = Const_Iso
!  Directional terms
!    -- Must use directional Ocean reflectance terms for the NON-ISOTROPY case
!      if ( do_isotropy ) then
!         WLeaving_SD(J,1:nstreams) = Const_Iso * Trans_Stream(1:nstreams)
!         WLeaving_SV(J,1:nvzas)    = Const_Iso * Trans_Viewing(1:nvzas)
!      else
!         do I = 1, nstreams
!            Rwprime = Ocean_SD_Reflecs(J,I) / ( one - Albedo * Ocean_SD_Reflecs(J,I)  )
!            Const = ( one / refrac_sq ) * trans_solar(J) * Rwprime
!            Const = Const * Foam_correction
!            WLeaving_SD(J,I) = Const * Trans_Stream(I)
!         enddo
!         do I = 1, nvzas
!            Rwprime = Ocean_SV_Reflecs(J,I) / ( one - Albedo * Ocean_SV_Reflecs(J,I)  )
!            Const = ( one / refrac_sq ) * trans_solar(J) * Rwprime
!            Const = Const * Foam_correction
!            WLeaving_SV(J,I) = Const * Trans_Viewing(I)
!         enddo
!      endif

!  Alternative. Should really use the F-term in ( one - Albedo * f * Ocean_Reflec_Basic  )
!               slightly larger values

!  Basic quantities

      mu0 = SZA_cosines(J)
      f    = 0.6279 - (0.2227*eta) - (0.0513*eta*eta) + (0.2465*eta -0.3119 ) * mu0
      denom =  ( one - Albedo * f * Ocean_Reflec_Basic )
      Rwprime = one / denom
      Const0       = ( one / refrac_sq ) * trans_solar(J)
      Const1       = Const0 * Foam_correction
      Const_Iso    = Const1 * Rwprime

!  C-derivatives
  
      dC_f = ( - 0.2227 - 0.1026*eta + 0.2465* mu0 ) * dC_eta
      dC_Rwprime = Rwprime * Rwprime * Albedo &
                  * ( dC_f * Ocean_Reflec_Basic + f * dC_Ocean_Reflec_Basic )
      dC_Const_Iso = Const1 * dC_Rwprime

! @@@@@@@@@@@@@@@@@@@@@ PATCH, 12 APRIL 2018. Add Mu0 factor @@@@@@@@@@@@@@@@@@@@@@@@@@@

!  R. Spurr and A. Sayer, 11-12 April 2018

!  Assign Isotropic Terms
!    -- Multiplied by mu0

      WLeaving_Iso(J)    =      Const_Iso *    Ocean_Iso_Reflecs(J) * mu0
      dC_WLeaving_Iso(J) = ( dC_Const_Iso *    Ocean_Iso_Reflecs(J) &
                              + Const_Iso * dC_Ocean_Iso_Reflecs(J) ) * mu0

!  Assign Directional terms
!    -- Must use directional Ocean reflectance terms for the NON-ISOTROPY case
!    -- Multiplied by mu0 for non-isotropic case

      if ( do_isotropy ) then
         WLeaving_SD(J,1:nstreams)    =    WLeaving_Iso(J) * Trans_Stream(1:nstreams)
         WLeaving_SV(J,1:nvzas)       =    WLeaving_Iso(J) * Trans_Viewing(1:nvzas)
         dC_WLeaving_SD(J,1:nstreams) = dC_WLeaving_Iso(J) * Trans_Stream(1:nstreams)
         dC_WLeaving_SV(J,1:nvzas)    = dC_WLeaving_Iso(J) * Trans_Viewing(1:nvzas)
      else
         do I = 1, nstreams
            WLeaving_SD(J,I)    = Const_Iso * Ocean_SD_Reflecs(J,I) * Trans_Stream(I) * mu0
            dC_WLeaving_SD(J,I) = ( dC_Const_Iso *    Ocean_SD_Reflecs(J,I) &
                                     + Const_Iso * dC_Ocean_SD_Reflecs(J,I) ) * Trans_Stream(I) * mu0
         enddo
         do I = 1, nvzas
            WLeaving_SV(J,I)    = Const_Iso * Ocean_SV_Reflecs(J,I) * Trans_Viewing(I) * mu0
            dC_WLeaving_SV(J,I) = ( dC_Const_Iso *    Ocean_SV_Reflecs(J,I) &
                                     + Const_Iso * dC_Ocean_SV_Reflecs(J,I) ) * Trans_Viewing(I) * mu0
         enddo
      endif

!  Rough Surface wind-speed derivatives
!    -- Multiplied by mu0.

      if ( do_Rough_Surface ) then
         dW_const0  = Const0 * dW_trans_solar(J) / trans_solar(J)
         dW_const1  = Const0 * dW_Foam_correction + dW_Const0 * Foam_correction
         dW_const_Iso = dW_Const1 * Rwprime
         dW_WLeaving_Iso(J) = dW_const_Iso * Ocean_Iso_Reflecs(J) * mu0
         if ( do_isotropy ) then
            dW_WLeaving_SD(J,1:nstreams) = dW_WLeaving_Iso(J) * Trans_Stream(1:nstreams) &
                                           +  WLeaving_Iso(J) * dW_Trans_Stream(1:nstreams)
            dW_WLeaving_SV(J,1:nvzas)    = dW_WLeaving_Iso(J) * Trans_Viewing(1:nvzas) &
                                           +  WLeaving_Iso(J) * dW_Trans_Viewing(1:nvzas)
         else
            do I = 1, nstreams
               dW_WLeaving_SD(J,I) = Ocean_SD_Reflecs(J,I) * &
                                ( dW_Const_Iso*Trans_Stream(I) + Const_Iso*dW_Trans_Stream(I) ) * mu0
            enddo
            do I = 1, nvzas
               dW_WLeaving_SV(J,I) = Ocean_SV_Reflecs(J,I) * &
                                ( dW_Const_Iso*Trans_Viewing(I) + Const_Iso*dW_Trans_Viewing(I) ) * mu0
            enddo
         endif
      endif

! @@@@@@@@@@@@@@@@@@@@@ END PATCH, 12 APRIL 2018. Add Mu0 factor @@@@@@@@@@@@@@@@@@@@@@@@@@@

! A Sayer 22 Sep 2015

      ! Normalisation factor Ta/Q. Note the prior version used mu0/pi, which gave
      ! similar numbers often but was not quite correct.

      ! This is from Gordon and Wang, AO, 1994, and is an approximate calculation.
      ! t = exp(-0.5*tau_rayleigh/mu0)

! R Spurr 01-05 October 2015

      ! Introduced the flag to control computation of this approximation (Now in Masater routine)
      ! More Exact calculation inside VLIDORT now programmed and tested
      ! Read a pre-calculated Rayleigh-atmopshere data set. Somewhat of a fiddle

      if ( Do_Approximate_Ta ) then
         Ta=exp(-0.5*tau_rayleigh/SZA_cosines(J))
      else
         Call TaSav_Interpolate &
            ( nTa_szas, nTa_wavs, Ta_szas, Ta_wavs, TaSavData, Wavelength, mu0, Ta)
      endif

!  Save the output

      TaSav(j,1) = Ta

! 25 September 2015 A Sayer

      ! I think that these terms may also need to be multiplied by the TOA solar flux. Presently we are taking
      ! this as equal to 1, so there is no issue. But the flux can be changed in VLIDORT so this parameter
      ! should really be propagated down here. Rob, can you check and, if necessary, do that? It might not be
      ! needed, if this value is scaled by the solar flux somewhere else.

! R Spurr 01 October 2015
!   17 December 2015. This is not required if we use the TRANSFLUX routine IN VLIDORT masters

      ! Checked. Multiplication by F0 is not necessary. It is the job of the water-leaving supplement
      ! to provide flux-normalized radiance sources.

!      WLeaving_Iso(J) = WLeaving_Iso(J)*Ta
!      do i = 1, nstreams
!         WLeaving_SD(J,I) = WLeaving_SD(J,I)*Ta
!      enddo
!      do i = 1, nvzas
!         WLeaving_SV(J,I) = WLeaving_SV(J,I)*Ta
!      enddo
! And I think we need to do the same for the gradients, too.
!      dC_WLeaving_Iso(J)=dC_WLeaving_Iso(J)*Ta
!      dW_WLeaving_Iso(J)=dW_WLeaving_Iso(J)*Ta
!      do i = 1, nstreams
!         dC_WLeaving_SD(J,I) = dC_WLeaving_SD(J,I)*Ta
!         dW_WLeaving_SD(J,I) = dW_WLeaving_SD(J,I)*Ta
!      enddo
!      do i = 1, nvzas
!         dC_WLeaving_SV(J,I) = dC_WLeaving_SV(J,I)*Ta
!         dW_WLeaving_SV(J,I) = dW_WLeaving_SV(J,I)*Ta
!      enddo

!  End solar loop

   enddo

!  Normal return

   return

!mick fix 9/19/2017 - added error return section
!  Error return

88 continue
   fail = .true.
   message = 'Openfile error in Linearized_WaterLeaving_2; file not found: ' // Trim(TaRayFile)

end subroutine Linearized_WaterLeaving_2

!

subroutine Lin_WhiteCap_Reflectance &
    ( WindSpeed, Wavelength, WC_Reflectance, WC_Lambertian, DWC_Reflectance, DWC_Lambertian )

!  Stand-alone routine for computing the WhiteCap Reflectance
!   Based on 6S code, as updated by A. Sayer (2011)

!  Linearization with respect to Wind-speed

   implicit none
   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Inputs
!    (Wind speed in [m/s], Wavelength in Microns)

   real(fpk), intent(in)  :: WindSpeed
   real(fpk), intent(in)  :: Wavelength

!  output

   real(fpk), intent(out) :: WC_Reflectance
   real(fpk), intent(out) :: WC_Lambertian
   real(fpk), intent(out) :: DWC_Reflectance
   real(fpk), intent(out) :: DWC_Lambertian

!  Data
!  ----

!  Single precision

   real :: Effective_WCRef(39)

! effective reflectance of the whitecaps (Koepke, 1984)
! These are the original values - superseded, A Sayer 05 Jul 2011.
!      data Effective_WCRef/ &
!     0.220,0.220,0.220,0.220,0.220,0.220,0.215,0.210,0.200,0.190,&
!     0.175,0.155,0.130,0.080,0.100,0.105,0.100,0.080,0.045,0.055,&
!     0.065,0.060,0.055,0.040,0.000,0.000,0.000,0.000,0.000,0.000,&
!     0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000/

! effective reflectance of the whitecaps (Frouin et al, 1996)
! Assume linear trends between the node points they give
! This is the spectral shape

      data Effective_WCRef/ &
     1.000,1.000,1.000,1.000,0.950,0.900,0.700,0.550,0.500,0.450,&
     0.400,0.350,0.300,0.250,0.200,0.150,0.100,0.050,0.000,0.000,&
     0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,&
     0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000/

!  Local variables
!  ---------------

!  Single precision in the original code

   integer :: iwl, iref
   real    :: Wlb, DWlb, WLP, Ref(39), wspd, wl, Ref_i, Rwc, DRwc

!  Initialize

   WC_Reflectance = 0.0_fpk
   WC_Lambertian  = 0.0_fpk
   DWC_Lambertian = 0.0_fpk

!  Single precision inputs in the original

   wspd = real(WindSpeed)
   wl   = real(Wavelength)

!  Scale data for value of 0.22 in the midvisible.

   DO iref = 1,39
      Ref(iref) = 0.22 * Effective_WCRef(iref)
   ENDDO

!  COMPUTE WHITECAPS REFLECTANCE (LAMBERTIAN)

! Bugfixed to make sure whitecap fraction never goes negative
! (i.e. check for ws < 3.7) A Sayer 21 Feb 2017

   Wlb    = 0.0 ; DWlb = 0.0
   IF (wspd .le. 9.25) THEN
      Wlb  = 0.01*((3.18e-03)*((wspd-3.7)**3.0))
      DWlb = 0.03*((3.18e-03)*((wspd-3.7)**2.0))
   ELSE IF (wspd .gt. 9.25) THEN
      Wlb  = 0.01*((4.82e-04)*((wspd+1.8)**3.0))
      DWlb = 0.03*((4.82e-04)*((wspd+1.8)**2.0))
   END IF
   IF (wspd .le. 3.7) THEN
      Wlb = 0.0 ; DWlb = 0.0
   END IF

! Original whitecap calculation - superseded, A. Sayer 05 Jul 2011.
!      W=2.95e-06*(wspd**3.52)

!  Find data point, Linearly interpolate

   iwl   = 1+int((wl-0.2)/0.1)
   wlp   = 0.5+(iwl-1)*0.1
   Ref_i = Ref(iwl+1) + ( wl-wlp)/0.1*(Ref(iwl)-Ref(iwl+1))
   Rwc   = Wlb*Ref_i
   DRwc  = DWlb*Ref_i

!  Final values

   WC_Lambertian   = real(Wlb,fpk)
   DWC_Lambertian  = real(DWlb,fpk)
   WC_Reflectance  = real(Rwc,fpk)
   DWC_Reflectance = real(DRwc,fpk)

!  Finish

   return
end subroutine Lin_WhiteCap_Reflectance

!

subroutine Lin_Ocean_Reflectance_Basic &
       ( Wavelength, PigmentConc, noWL, Ocean_Reflec, eta_out, dC_Ocean_Reflec, dC_eta_out )

!  THIS IS FORMERLY CALLED "MORCASIWAT", as modified by A. Sayer for 6S

! mick fix 12/28/2014 - Using updates by A Sayer November 03 2014:
! Extended functionality down to 200 nm. Achieved by:
! - Extended data arrays down to 200 nm (another 40 elements).
! - Changed logic check for contribution to 0.2-0.9 microns from 0.4-0.9 microns, and started table
!   lookup calculation from 0.2 microns instead of 0.4 microns.
! Note, this is based on a simple extension of the published optical model for vis wavelengths.
!   Possible that other scatterers/absorbers.
! which are neglected in this model may be important at UV wavelengths.
! Do linear interpolation of optical property LUTs, rather than nearest neighbour, to remove discontinuities. Achieved by:
! - Replicated final element of LUTs to avoid potential for extrapolation errors.
! - Replace nint() call with floor() call to correctly get lower bound
! - Define variable dwl, fractional distance along the 5 nm LUT grid
! - Implement the interpolation using dwl rather than direct lookup of nearest value.
! Also:
! - Corrected Prieur and Sathyendranath, Limnol. Oceanogr. reference year to 1981 instead of 1983.
! - Corrected typo in water scattering coefficient at 410 nm: 0.0068 was written instead of 0.0061.
!   Removes artificial spike at 410 nm in calculated reflectance.

! Updated A Sayer August 07 2015:

! - Updated pure water absorption coefficient from Smith and Baker (1981) to Lee et al (2015), up to 550 nm. This has the
!   effect of decreasing water absorption, particularly in the blue and UV.
! - Updated pure water absorption coefficient below 350 nm using Quickenden and Irvin, 1980:
!   http://scitation.aip.org/content/aip/journal/jcp/72/8/10.1063/1.439733
! - Updated pure water absorption coefficient between 555 nm and 725 nm to Pope and Fry (1997).
! - Updated pure water absorption coefficient above 725 nm to Hale and Querry (1973).
! - Updated water scattering from Smith and Baker (2009) to Zhongping Lee's analytical summary from Zhang et al (2009).
! Updated A Sayer September 22 2015:
! - Bugfix of R=f*(b/a) instead of R=f*(b/[a+b]). Thanks to A. Vasilkov for pointing this out. Note effect is minor at midvisible
!   and longer wavelengths.

! Updated A Sayer September 28 2015:

! - Use Vasilkov et al (2005), Applied optics 44(14), 2863-2869, to calculate Chl absorption for 400 nm or lower. This uses
!   a different set of coefficients to the current Chl calculation. Note source data are at 2 nm intervals but as it is fairly
!   linear, and there is some scatter about it anyway, subsample to 10 nm intervals instead.
! Updated A Sayer September 29 2015:
! - Updated Chl absorption from 400-720 nm using empirical model from Lee et al (AO 37(27), 1998. This accounts for Chl-dependence
!   of spectral shape of pigement absorption, replacing Prieur and Sathyendranath (1981) spectrum.
! Updated A Sayer September 30 2015:
! - Instead of having f in here and Q in the main WaterLeaving part, we now use the ratio f/Q ('foQ') in here. These data were
!   provided by David Antoine, described in Morel et al (AO, 2002), doi: 10.1364/AO.41.006289, and are still considered current
!   and valid.
!   This means that we now account more fully for the bidirectional nature of the underlight term.
!   Note that for now I have put in a placeholder for foQ=0.09; Rob Spurr will implement the main formulation.

! Updated R. Spurr October 01-06 2015:

!mick fix 9/19/2017 - upgraded these calculations to double precision

   implicit none
   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Input/Output
!  ------------

!  Angle dependence has been removed. R. Spurr, 03 October 2015
!   integer  , intent(in)   :: Maxszas, nszas
!   real(fpk), intent(in)   :: szas  (Maxszas)

   real(fpk), intent(in)   :: Wavelength
   real(fpk), intent(in)   :: PigmentConc

   logical  , intent(out)  :: noWL
   real(fpk), intent(out)  :: Ocean_Reflec, eta_out
   real(fpk), intent(out)  :: dC_Ocean_Reflec, dC_eta_out      ! Derivatives w.r.t PigmentConc

!  Local
!  -----

!      subroutine morcasiwat(wl,C,R2,mu_sol)
! Rewritten, beginning 07 July 2011, Andrew Sayer
! Now extends underlight calculations out to 400-900 nm.
! and more modern formulations,
! but still based on Case 1 principle.

! Spectral diffuse attenuation coefficient of Case I Waters as Predicted
! by MOREL within the spectral range 400-700nm (1988, Journal of Geophysical
! Research, Vol.93, No C9, pp 10749-10768)
!
! input parameters:	wl wavelength (IN MICROMETERS)
!			C  pigment concentration
!                       mu_sol : cosine of solar zenith angle
! output parameter:	R2  reflectance of water below the surface

! Tabulated absorption/scattering coefficients. See comments by data arrays for info.
      real(fpk) water_abs_coeff(142)
      real(fpk) vasilkov_chl_coeff_a(11),vasilkov_chl_coeff_b(11),lee_chl_coeff_a(33),lee_chl_coeff_b(33)
! Input/output parameters
      real(fpk) r2,C,wl
! Absorption/scattering terms, and parameters for calculation of f
      real(fpk) a_wat,b_wat,a_chl,a_tot,b_tot,a_ph,a_cdom,v,bp,bbp,eta,dwl
      real(fpk) vas_chl_a,vas_chl_b,lee_chl_a,lee_chl_b,a_440,helpw
! Wavelength index for tables
      integer iwl
! f/Q calculation. Now Outside the routine. R. Spurr 03 October 2015
!      real foQ, foQ_Int_1 ( maxszas)
!      character(Len=100) FoQPath

!  derivative variables

      real(fpk) da_tot,db_tot,da_ph,da_cdom,dv,dbp,dbbp,deta,dx,dz,x,z,dR2
   
! 07 AUG 2015 Updated A Sayer

! Quickenden and Irvin (1980) from 200-320 nm (http://scitation.aip.org/content/aip/journal/jcp/72/8/10.1063/1.439733),
! interpolate with Lee et al. (2015) between 325 and 345 nm
! Lee et al. (2015) for 350-550 nm. (https://www.osapublishing.org/ao/abstract.cfm?uri=ao-54-3-546)
! Pope and Fry (1997) for 555-725 nm. These are pretty similar to Lee et al in the 520-550 nm range.
! Hale & Qurry, AO(12) 1973, table 1, for 725-900 nm. This has 25 nm increments so
! linearly interpolate between these. Provided as extinction coefficient
! so convert using a=4*pi*k/lambda (note lambda in m for units m^{-1})
! From 200 - 900 nm in 5 nm increments.
       data water_abs_coeff/ &
       0.3240_fpk,0.1560_fpk,0.1260_fpk,0.1010_fpk,0.0850_fpk,&  ! 200-220
       0.0650_fpk,0.0595_fpk,0.0542_fpk,0.0483_fpk,0.0392_fpk,&  ! 225-245
       0.0376_fpk,0.0326_fpk,0.0308_fpk,0.0251_fpk,0.0236_fpk,&  ! 250-270
       0.0216_fpk,0.0222_fpk,0.0178_fpk,0.0163_fpk,0.0145_fpk,&  ! 275-295
       0.0124_fpk,0.0112_fpk,0.0112_fpk,0.0105_fpk,0.0100_fpk,&  ! 300-320
       0.0095_fpk,0.0091_fpk,0.0085_fpk,0.0080_fpk,0.0075_fpk,&  ! 325-345
       0.0071_fpk,0.0062_fpk,0.0056_fpk,0.0050_fpk,0.0046_fpk,&  ! 350-370
       0.0041_fpk,0.0037_fpk,0.0035_fpk,0.0032_fpk,0.0032_fpk,&  ! 375-395
       0.0032_fpk,0.0032_fpk,0.0031_fpk,0.0031_fpk,0.0032_fpk,&  ! 400-420
       0.0033_fpk,0.0036_fpk,0.0038_fpk,0.0044_fpk,0.0054_fpk,&  ! 425-445
       0.0068_fpk,0.0073_fpk,0.0076_fpk,0.0081_fpk,0.0089_fpk,&  ! 450-470
       0.0099_fpk,0.0109_fpk,0.0118_fpk,0.0132_fpk,0.0154_fpk,&  ! 475-495
       0.0187_fpk,0.0230_fpk,0.0302_fpk,0.0368_fpk,0.0387_fpk,&  ! 500-520
       0.0400_fpk,0.0418_fpk,0.0443_fpk,0.0470_fpk,0.0507_fpk,&  ! 525-545
       0.0562_fpk,0.0596_fpk,0.0619_fpk,0.0642_fpk,0.0695_fpk,&  ! 550-570
       0.0772_fpk,0.0896_fpk,0.1100_fpk,0.1351_fpk,0.1672_fpk,&  ! 575-595
       0.2224_fpk,0.2577_fpk,0.2644_fpk,0.2678_fpk,0.2755_fpk,&  ! 600-620
       0.2834_fpk,0.2916_fpk,0.3012_fpk,0.3108_fpk,0.3250_fpk,&  ! 625-645
       0.3400_fpk,0.3710_fpk,0.4100_fpk,0.4290_fpk,0.4390_fpk,&  ! 650-670
       0.4480_fpk,0.4650_fpk,0.4860_fpk,0.5160_fpk,0.5590_fpk,&  ! 675-695
       0.6240_fpk,0.7040_fpk,0.8270_fpk,1.0070_fpk,1.2310_fpk,&  ! 700-720
       1.4890_fpk,1.7990_fpk,2.0895_fpk,2.3800_fpk,2.4250_fpk,&  ! 725-745
       2.4700_fpk,2.5100_fpk,2.5500_fpk,2.5300_fpk,2.5100_fpk,&  ! 750-770
       2.4350_fpk,2.3600_fpk,2.2600_fpk,2.1600_fpk,2.1150_fpk,&  ! 775-795
       2.0700_fpk,2.2104_fpk,2.3509_fpk,2.4913_fpk,2.6318_fpk,&  ! 800-820
       2.7722_fpk,3.0841_fpk,3.3960_fpk,3.7079_fpk,4.0198_fpk,&  ! 825-845
       4.3317_fpk,4.5884_fpk,4.8451_fpk,5.1019_fpk,5.3586_fpk,&  ! 850-870
       5.6153_fpk,5.8495_fpk,6.0836_fpk,6.3177_fpk,6.5518_fpk,&  ! 875-895
       6.7858_fpk,6.7858_fpk/                                    ! 900-905

! Coefficients to calculate Chl absorption for 300-400 nm. Use these as 400 nm and below. Use 300 nm data for wavelengths <300 nm.
! From Vasilkov et al (2005), Applied optics 44(14), 2863-2869, table 1.
! From 300 - 400 nm in 10 nm increments.
! Note that this does lead to a discontinuity in predicted reflectance at 400 nm for some geometries and Chl concentrations.

       data vasilkov_chl_coeff_a/&
       0.1023_fpk,0.0986_fpk,0.0958_fpk,0.0924_fpk,0.0841_fpk,&
       0.0709_fpk,0.0604_fpk,0.0580_fpk,0.0556_fpk,0.0532_fpk,&
       0.0520_fpk/

       data vasilkov_chl_coeff_b/&
       0.0983_fpk,0.103_fpk,0.120_fpk,0.137_fpk,0.147_fpk,&
       0.153_fpk,0.160_fpk,0.172_fpk,0.182_fpk,0.191_fpk,&
       0.198_fpk/

! Coefficient to calculate Chl absorption for 400-720 nm.
! See Table 2 and Equation 12 of Lee et al (1998), Applied Optics 37(27), 6329-6338.
       data lee_chl_coeff_a/&
       0.6843_fpk,0.7782_fpk,0.8637_fpk,0.9603_fpk,1.0000_fpk,&
       0.9634_fpk,0.9311_fpk,0.8697_fpk,0.7890_fpk,0.7558_fpk,&
       0.7333_fpk,0.6911_fpk,0.6327_fpk,0.5681_fpk,0.5046_fpk,&
       0.4262_fpk,0.3433_fpk,0.2950_fpk,0.2784_fpk,0.2595_fpk,&
       0.2389_fpk,0.2745_fpk,0.3197_fpk,0.3421_fpk,0.3331_fpk,&
       0.3502_fpk,0.5610_fpk,0.8435_fpk,0.7485_fpk,0.3890_fpk,&
       0.1360_fpk,0.0545_fpk,0.0250_fpk/
       
       data lee_chl_coeff_b/&
       0.0205_fpk,0.0129_fpk,0.0060_fpk,0.0020_fpk,0.0000_fpk,&
       0.0060_fpk,0.0109_fpk,0.0157_fpk,0.0152_fpk,0.0256_fpk,&
       0.0559_fpk,0.0865_fpk,0.0981_fpk,0.0969_fpk,0.0900_fpk,&
       0.0781_fpk,0.0659_fpk,0.0600_fpk,0.0581_fpk,0.0540_fpk,&
       0.0495_fpk,0.0578_fpk,0.0674_fpk,0.0718_fpk,0.0685_fpk,&
       0.0713_fpk,0.1128_fpk,0.1595_fpk,0.1388_fpk,0.0812_fpk,&
       0.0317_fpk,0.0128_fpk,0.0050_fpk/

!  Initialize

      wl     = Wavelength
      C      = PigmentConc
      noWL = .false.  ; Ocean_Reflec = 0.0_fpk ; dC_Ocean_Reflec = 0.0_fpk

! If wavelength out of range, no need to calculate underlight

      if (wl.lt.0.200_fpk.or.wl.gt.0.900_fpk)then
        noWL = .true. ; goto 60
      endif

      ! Get water absorption.
      ! Extract tabulated values for this wavelength
      iwl=1+floor((wl-0.200_fpk)/0.005_fpk)
      dwl=(wl-0.200_fpk)/0.005_fpk-floor((wl-0.200_fpk)/0.005_fpk)
      a_wat=water_abs_coeff(iwl)+dwl*(water_abs_coeff(iwl+1)-water_abs_coeff(iwl))

      ! Get Chl absorption. Updated A. Sayer, September 2015.
      ! If the wavelength is <400 nm, then use Vasilkov et al (AO, 2005) to get a_ph.
      if (wl.lt.0.400_fpk)then
         ! Extract tabulated values for this wavelength
         iwl=1+floor((wl-0.300_fpk)/0.01_fpk)
         dwl=(wl-0.300_fpk)/0.01_fpk-floor((wl-0.300_fpk)/0.01_fpk)
         if (wl.lt.0.300_fpk) then ! Use 300 nm values for wavelengths below 300 nm
            iwl=1
            dwl=0
         endif
         vas_chl_a=vasilkov_chl_coeff_a(iwl)+dwl*(vasilkov_chl_coeff_a(iwl+1)-vasilkov_chl_coeff_a(iwl))
         vas_chl_b=vasilkov_chl_coeff_b(iwl)+dwl*(vasilkov_chl_coeff_b(iwl+1)-vasilkov_chl_coeff_b(iwl))
         a_chl=vas_chl_a*(C**(-vas_chl_b))
         a_ph=C*a_chl
         da_ph = a_chl*(1.0_fpk-vas_chl_b)
      endif
      ! If the wavelength is 400 nm or above, use Lee et al (1998) data. Use 720 nm data for wavelengths > 720 nm.
      ! This has a minor influence because water absorption dominates at 720 nm and longer.
      if (wl.ge.0.400_fpk) then
         ! Extract tabulated values for this wavelength
         iwl=1+floor((wl-0.400_fpk)/0.01_fpk)
         dwl=(wl-0.400_fpk)/0.01_fpk-floor((wl-0.400_fpk)/0.01_fpk)
         if (wl.gt.0.720_fpk) then ! Use 720 nm values for wavelengths above 720 nm
            iwl=33
            dwl=0
         endif
         a_440=0.06_fpk*(C**0.65_fpk) ! e.g. Morel and Maritorena, 2001; also used by Lee et al papers
         lee_chl_a=lee_chl_coeff_a(iwl)+dwl*(lee_chl_coeff_a(iwl+1)-lee_chl_coeff_a(iwl))
         lee_chl_b=lee_chl_coeff_b(iwl)+dwl*(lee_chl_coeff_b(iwl+1)-lee_chl_coeff_b(iwl))
         a_ph=(lee_chl_a+lee_chl_b*log(a_440))*a_440
         da_ph = 0.65_fpk * (a_ph + lee_chl_b*a_440 ) / C
      endif

      ! Get CDOM absorption.
      ! Equations 2 and 4a from Morel and Gentili, RSE 113, 2009.
      ! This is assuming that CDOM absorption at reference wavelength is coupled to Chl concentration.
      ! Revision. 11/28/18, upgrade for VLIDORT 2.8a, 3/18/19.

!  Here is the Old code
!       a_cdom=0.0524_fpk*(C**0.63_fpk)*exp(-0.018_fpk*(wl*1000.0_fpk-412.0_fpk))
!       da_cdom = 0.63_fpk * a_cdom / C

      ! Revision. 11/28/18, upgrade for VLIDORT 2.8a, 3/18/19.
!  changed to the commonly accepted equation from Morel&Maritorena, JGR, 2001;
!  value of pure water absorption at 440 nm aw=0.0065 is taken from Morel et
!  al, Limnol. Oceanogr, 2007
!     a_cdom=0.2*(0.00635+0.06*(C**0.65))*exp(-0.014*(wl-0.440))
      !  the correct form if wl in microns, 11/28/18
      helpw   = 0.2_fpk*exp(-0.014_fpk*(wl*1000.0_fpk-440.0_fpk))
      a_cdom  = helpw * (0.00635_fpk+0.06_fpk*(C**0.65_fpk))
      da_cdom = helpw * 0.06_fpk*0.65_fpk*(C**(-0.35_fpk))

!  Total

      a_tot=a_wat + a_ph + a_cdom
      da_tot  = da_ph + da_cdom

! 07 AUG 2015 Updated b_wat - Zhongping Lee's quick form of Zhang et al (2009)
!    https://www.osapublishing.org/oe/abstract.cfm?uri=oe-17-7-5698
!      b_wat=0.0026_fpk*(0.42_fpk/wl)**4.3_fpk
      
! Revision. 11/28/18, upgrade for VLIDORT 2.8a, 3/18/19.
!     b_wat value from Morel et al., Limonol. Oceanogr, 2007
!     Simply change the coefficient of 0.0026 to 0.0028 so that you have

      b_wat=0.0028_fpk*(0.42_fpk/wl)**4.3_fpk

! Morel and Maritorena, 2001 (also earlier work by Loisel and Morel)
! exponent for spectral dependence of scattering
      x=log10(C) ; dx = 1.0_fpk/C/log(10.0_fpk) 
      if (C .le. 2.0_fpk) then
        v=0.5_fpk*(x-0.3_fpk) ; dv = 0.5_fpk*dx
      endif
      if (C .gt. 2.0_fpk) then
        v=0.0_fpk ; dv = 0.0_fpk
      endif

      bp=0.416_fpk*(C**0.766_fpk)
      dbp = 0.766_fpk*bp/C

      z = (wl/0.55_fpk)**v
      dz = z * dv* log(wl/0.55_fpk)

      bbp  = 0.002_fpk+0.01_fpk*(0.5_fpk-0.25_fpk*x)*z
      dbbp = -0.0025_fpk*dx*z+0.01_fpk*(0.5_fpk-0.25_fpk*x)*dz

      b_tot=b_wat + bbp*bp ; db_tot = bbp * dbp + dbbp * bp
      eta=b_wat/b_tot      ; deta = -eta * db_tot / b_tot
      eta_out = real(eta,fpk) ; dC_eta_out = real(deta,fpk)

!  Basic reflectance, before foQ factor
!     R2 = b_tot/(b_tot+a_tot)  ! Alternative. Better choice after talking w/ Sasha on 4/12/18 wqin

      R2 = (b_tot/a_tot)                 ! Sasha's suggestion on 11/28/18, Implemented 3/18/19 Version 2.8a
      dR2 = ( db_tot - R2*da_tot ) /a_tot
      Ocean_Reflec    = real(R2,fpk)
      dC_Ocean_Reflec = real(dR2,fpk)

!  Former code from earlier versions, with comments

! Morel and Gentili, 1991.
!     do j = 1, nszas
! Now we have taken the Q factor out of main WaterLeaving routine,
! replace f calculation here with the f/Q ratio from Morel et al (AO, 2002).
!        mu_sol = real(SZA_cosines(J))
!        f=0.6279 - (0.2227*eta) - (0.0513*eta*eta) + (0.2465*eta -0.3119 )*mu_sol
!        R2=f*(b_tot/a_tot)
! R. Spurr 01 October 2015
!        foQ = 0.09 ! Placeholder until we have implementation of David Antoine's lookup table.
!        Removed placeholder and inserted interpolated values
!        foQ=foQ_Int_1(J) 
!        write(77,'(2f8.3,f10.6)')  wl*1000.0d0, Log(C),foQ
!        R2=foQ*(b_tot/a_tot)
!        Ocean_Reflecs(J) = real(R2,fpk)
!     enddo

!  continuation point

 60  continue

!  Finish

      return
end subroutine Lin_Ocean_Reflectance_Basic

!

subroutine Interpolate_Lin_fOQ_BS1 &
      ( FoQFile, Maxszas, nszas, szas, Wavelength, PigmentConc, &
        foQ_Int_1, dC_foQ_Int_1, fail, message )

!  I/O double precision. Local computations are all single precision

   implicit none
   integer, parameter :: fpk = SELECTED_REAL_KIND(15)
   integer, parameter :: spk = SELECTED_REAL_KIND(6)

!  Input/Output
!  ------------

!  input is double precision from the main routine

   character*(*), intent(in)   :: FoQFile
   integer      , intent(in)   :: Maxszas, nszas
   real(fpk)    , intent(in)   :: szas  (Maxszas)
   real(fpk)    , intent(in)   :: Wavelength
   real(fpk)    , intent(in)   :: PigmentConc

!  output is double precision

   real(fpk)    , intent(out)  :: foQ_Int_1(Maxszas)
   real(fpk)    , intent(out)  :: dC_foQ_Int_1(Maxszas)
   logical      , intent(out)  :: fail
   character*(*), intent(out)  :: message

!  Local
!  -----

!  Table information (single precision)

   real(spk) :: Logpigs(6), Lambdas(7), cossuns(6), lams(7), pigs(6), suns(6)
   real(spk) :: FoQ_Table(17,13), fOQ_averaged(7,6,6)

   data suns / 0.0, 15.0, 30.0, 45.0, 60.0, 75.0 /
   data pigs / 0.03, 0.1, 0.3, 1.0, 3.0, 10.0    /
   data lams / 412.5, 442.5, 490.0, 510.0, 560.0, 620.0, 660.0 /

!  help variables (all single precision)

   real(spk) :: dtr, fw1, fw2, fs1, fs2, cszas(Maxszas), yspline(6), bbas(6), cbas(6), dbas(6)
   real(spk) :: Wave, csza, C, LogC, lam, sun, pigc, fval, Logder, foQ_suns(6), dC_foQ_suns(6)
   integer   :: i, j, k, m, n, w1, w2, s1, s2, j1

!  initialize

   FoQ_Int_1    = 0.0_fpk
   dC_FoQ_Int_1 = 0.0_fpk
   fail = .false. ; message = ' '

!  Basic interpolation quantities

   dtr = acos(-1.0)/180.0
   do j = 1, 6
      cossuns(7-j) = cos(suns(j)*dtr)
   enddo
   do k = 1, 6
      LogPigs(k) = log(pigs(k))
   enddo
   LAMBDAS(1:7) = lams(1:7)

!  Develop cosine streams

   do j = 1, nszas
      cszas(j) = cos(real(szas(j),spk)*dtr)
   enddo

!  Obtain table averages

   open(45,file=Trim(FoQFile),err=88,status='old',action='read')
   do i = 1, 7
     do j = 1, 6
       do k = 1, 6
         read(45,*)lam, sun, pigc, fval ;j1 = 7-j
         do m = 1, 17
           read(45,*)FoQ_Table(m,1:13)
         enddo
         fOQ_averaged(i,j1,k) = sum(FoQ_Table(1:17,1:13))/221.0
       enddo
     enddo
   enddo
   close(45)

!  Find boundaries

   Wave = real(Wavelength) * 1000.0 ! Convert to [nm]
   if ( Wave.le.Lambdas(1) ) then
      w1 = 1 ; w2 = w1 + 1 ; fw1 = 1.0 ; fw2 = 0.0
   else if ( Wave.ge.Lambdas(7) ) then
      w1 = 6 ; w2 = w1 + 1 ; fw1 = 0.0 ; fw2 = 1.0
   else
      do i = 1, 6
         if ( Wave.gt.Lambdas(i).and.Wave.le.Lambdas(i+1)) w1 = i
      enddo
      w2 = w1 + 1 ; fw1 = (Lambdas(w2)-Wave)/(Lambdas(w2)-Lambdas(w1)) ; fw2 = 1.0 - fw1
   endif

!  Carry out Linear wavelength and Splined Pigment interpolations for all solar angles
!      SUBROUTINE splint(xa,ya,y2a,n,x,y)
!      SUBROUTINE spline(x,y,n,yp1,ypn,y2)

   C = Real(PigmentConc) ; LogC = log(C)
   do j = 1, 6
      do k = 1, 6
         yspline(k) = fw1*foQ_averaged(w1,j,k) + fw2*foQ_averaged(w2,j,k) 
      enddo
      Call bspline(6,6,LogPigs,yspline,bbas,cbas,dbas)
      Call dSeval(6,6,LogC,LogPigs,yspline,bbas,cbas,dbas,foQ_suns(j),Logder) ; dC_foQ_suns(j) = LogDer / C
   enddo

!  Solar angles

   do n = 1, nszas
      csza = cos(real(szas(n))*dtr)
      if ( csza.le.cossuns(1) ) then
         s1 = 1 ; s2 = s1 + 1 ; fs1 = 1.0 ; fs2 = 0.0
      else if ( csza.ge.cossuns(6) ) then
         s1 = 5 ; s2 = s1 + 1 ; fs1 = 0.0 ; fs2 = 1.0
      else
         do j = 1, 5
            if ( csza .gt.cossuns(j).and.csza.le.cossuns(j+1)) s1 = j
         enddo
         s2 = s1 + 1 ; fs1 = (cossuns(s2)-csza)/(cossuns(s2)-cossuns(s1)) ; fs2 = 1.0 - fs1
      endif
      FoQ_Int_1(n)    = dble ( fs1 * foQ_suns(s1)    + fs2 * foQ_suns(s2)    )
      dC_FoQ_Int_1(n) = dble ( fs1 * dC_foQ_suns(s1) + fs2 * dC_foQ_suns(s2) )
!      write(*,*)'Sun',s1,s2,fs1,fs2

  enddo

!  Normal return

   return

!  Error return

88 continue
   fail = .true.
!mick fix 3/22/2017 - upgraded error msg to include file name
   message = 'Openfile error in Interpolate_Lin_fOQ_BS1; file not found: ' // Trim(FoQFile)

   return
end subroutine Interpolate_Lin_fOQ_BS1

!

subroutine Interpolate_Lin_fOQ_BS2 &
      ( FoQFile, Maxszas, Maxstreams, Maxvzas, nszas, nvzas, nstreams, &
        szas, vzas, streams, refrac_R, Wavelength, PigmentConc,        &
        foQ_Int_1, foQ_Int_SV, foQ_Int_SD,                             &
        dC_foQ_Int_1, dC_foQ_Int_SV, dC_foQ_Int_SD, fail, message ) 

!  I/O double precision. Local computations are all single precision

   implicit none
   integer, parameter :: fpk = SELECTED_REAL_KIND(15)
   integer, parameter :: spk = SELECTED_REAL_KIND(6)

!  Input/Output
!  ------------

!  input is double precision from the main routine

   character*(*), intent(in)   :: FoQFile

   integer      , intent(in)   :: Maxszas, Maxvzas, Maxstreams
   integer      , intent(in)   :: nszas, nvzas, nstreams
   real(fpk)    , intent(in)   :: szas   (Maxszas)
   real(fpk)    , intent(in)   :: vzas   (Maxvzas)
   real(fpk)    , intent(in)   :: streams(Maxstreams)

   real(fpk)    , intent(in)   :: Wavelength
   real(fpk)    , intent(in)   :: PigmentConc
   real(fpk)    , intent(in)   :: Refrac_R

!  output is double precision

   real(fpk)    , intent(out)  :: foQ_Int_1  ( Maxszas )
   real(fpk)    , intent(out)  :: foQ_Int_SD ( maxszas, Maxstreams )
   real(fpk)    , intent(out)  :: foQ_Int_SV ( maxszas, Maxvzas )
   real(fpk)    , intent(out)  :: dC_foQ_Int_1  ( Maxszas )
   real(fpk)    , intent(out)  :: dC_foQ_Int_SD ( maxszas, Maxstreams )
   real(fpk)    , intent(out)  :: dC_foQ_Int_SV ( maxszas, Maxvzas )
   logical      , intent(out)  :: fail
   character*(*), intent(out)  :: message

!  Local
!  -----

!  Table information (single precision)

   real(spk) :: Logpigs(6), Lambdas(7), cossuns(6), cosnads(17), lams(7), pigs(6), suns(6), nads(17)
   real(spk) :: FoQ_Table(17,13), fOQ_averaged(7,6,6), fOQ_averaged_2(7,6,6,17)

   data suns / 0.0, 15.0, 30.0, 45.0, 60.0, 75.0 /
   data pigs / 0.03, 0.1, 0.3, 1.0, 3.0, 10.0    /
   data lams / 412.5, 442.5, 490.0, 510.0, 560.0, 620.0, 660.0 /
   data nads /  1.078,  3.411,  6.289,  9.278, 12.300, 15.330, 18.370, 21.410, 24.450, &
               27.500, 30.540, 33.590, 36.640, 39.690, 42.730, 45.780, 48.830 /

!  help variables (all single precision)

   real(spk) :: fw1, fw2, fs1, fs2, fd1, fd2, yspline(6), bbas(6), cbas(6), dbas(6)
   real(spk) :: Wave, C, LogC, cs, cd, lam, sun, pigc, fval, foQ_1, foQ_2, Logder
   real(spk) :: foQ_suns(6), foQ_nads(6,17), dC_foQ_suns(6), dC_foQ_nads(6,17)
   real(spk) :: dtr, local_sine, local_cos, incident_angle, refrac_R_sp
   real(spk) :: cstreams(Maxstreams), cszas(Maxszas), cvzas(Maxvzas)
   integer   :: i, j, k, m, n, w1, w2, s1, s2, d1, d2, j1, m1

!  initialize

   FoQ_Int_1  = 0.0_fpk
   FoQ_Int_SV = 0.0_fpk
   FoQ_Int_SD = 0.0_fpk
   dC_FoQ_Int_1  = 0.0_fpk
   dC_FoQ_Int_SV = 0.0_fpk
   dC_FoQ_Int_SD = 0.0_fpk
   fail = .false. ; message = ' '

!  Basic interpolation quantities

   dtr = acos(-1.0)/180.0
   do m = 1, 17
     cosnads(18-m) = cos(nads(m)*dtr)
   enddo
   do j = 1, 6
      cossuns(7-j) = cos(suns(j)*dtr)
   enddo
   do k = 1, 6
      LogPigs(k) = log(pigs(k))
   enddo
   LAMBDAS(1:7) = lams(1:7)

!  Develop cosine streams
!    -- refract the viewing angles

   refrac_R_sp = real(refrac_R,spk)
   do i = 1, nstreams
      local_cos  = real (streams(i),spk)
      local_sine = sqrt ( 1.0 - local_cos * local_cos )
      incident_angle = asin(local_sine/refrac_R_sp)
      cstreams(i) = cos(incident_angle)
   enddo
   do i = 1, nvzas
      local_sine = sin(real(vzas(i),spk)*dtr)
      incident_angle = asin(local_sine/refrac_R_sp)
      cvzas(i) = cos(incident_angle)
   enddo
   do j = 1, nszas
      cszas(j) = cos(real(szas(j),spk)*dtr)
   enddo

!  Read table and obtain averages

   open(45,file=Trim(FoQFile),err=88,status='old',action='read')
   do i = 1, 7
     do j = 1, 6
       do k = 1, 6
         read(45,*)lam, sun, pigc, fval ; j1 = 7-j
         do m = 1, 17
           m1 = 18-m ; read(45,*)FoQ_Table(m,1:13)
           fOQ_averaged_2(i,j1,k,m1) = sum(FoQ_Table(m,1:13))/13.0
         enddo
         fOQ_averaged(i,j1,k) = sum(FoQ_Table(1:17,1:13))/221.0
       enddo
     enddo
   enddo
   close(45)

!  Find boundaries

   Wave = real(Wavelength) * 1000.0 ! Convert to [nm]
   if ( Wave.le.Lambdas(1) ) then
      w1 = 1 ; w2 = w1 + 1 ; fw1 = 1.0 ; fw2 = 0.0
   else if ( Wave.ge.Lambdas(7) ) then
      w1 = 6 ; w2 = w1 + 1 ; fw1 = 0.0 ; fw2 = 1.0
   else
      do i = 1, 6
         if ( Wave.gt.Lambdas(i).and.Wave.le.Lambdas(i+1)) w1 = i
      enddo
      w2 = w1 + 1 ; fw1 = (Lambdas(w2)-Wave)/(Lambdas(w2)-Lambdas(w1)) ; fw2 = 1.0 - fw1
   endif
!   write(*,*)'Wav',w1,w2,fw1,fw2

!  Carry out Linear wavelength and Splined Pigment interpolations for all solar angles
!      SUBROUTINE splint(xa,ya,y2a,n,x,y)
!      SUBROUTINE spline(x,y,n,yp1,ypn,y2)

   C = Real(PigmentConc) ; LogC = log(C) 
   do j = 1, 6
      do k = 1, 6
         yspline(k) = fw1*foQ_averaged(w1,j,k) + fw2*foQ_averaged(w2,j,k) 
      enddo
      Call bspline(6,6,LogPigs,yspline,bbas,cbas,dbas)
      Call dSeval(6,6,LogC,LogPigs,yspline,bbas,cbas,dbas,foQ_suns(j),Logder) ; dC_foQ_suns(j) = LogDer / C
      do m = 1, 17
         do k = 1, 6
            yspline(k) = fw1*foQ_averaged_2(w1,j,k,m) + fw2*foQ_averaged_2(w2,j,k,m) 
         enddo
         Call bspline(6,6,LogPigs,yspline,bbas,cbas,dbas)
         Call dSeval(6,6,LogC,LogPigs,yspline,bbas,cbas,dbas,foQ_nads(j,m),Logder) ; dC_foQ_nads(j,m) = LogDer / C
      enddo
   enddo

!  Solar angle loop

   do n = 1, nszas
      cs = cszas(n)
      if ( cs.le.cossuns(1) ) then
         s1 = 1 ; s2 = s1 + 1 ; fs1 = 1.0 ; fs2 = 0.0
      else if ( cs.ge.cossuns(6) ) then
         s1 = 5 ; s2 = s1 + 1 ; fs1 = 0.0 ; fs2 = 1.0
      else
         do j = 1, 5
            if ( cs .gt.cossuns(j).and.cs.le.cossuns(j+1)) s1 = j
         enddo
         s2 = s1 + 1 ; fs1 = (cossuns(s2)-cs)/(cossuns(s2)-cossuns(s1)) ; fs2 = 1.0 - fs1
      endif
      FoQ_Int_1(n)    = dble ( fs1 * foQ_suns(s1)    + fs2 * foQ_suns(s2)    )
      dC_FoQ_Int_1(n) = dble ( fs1 * dC_foQ_suns(s1) + fs2 * dC_foQ_suns(s2) )

!  Stream angles
!   -- R. Spurr, 12  April 2018.  Fixed glitch in f/Q interpolation

      do i = 1, nstreams
         cd = cstreams(i)
         if ( cd.le.cosnads(1) ) then
            d1 = 1 ; d2 = d1 + 1 ; fd1 = 1.0 ; fd2 = 0.0
         else if ( cd.ge.cosnads(17) ) then
!            d1 = 5 ; d2 = d1 + 1 ; fd1 = 0.0 ; fd2 = 1.0   !  Wrong offset, produced glitches
            d1 = 16 ; d2 = d1 + 1 ; fd1 = 0.0 ; fd2 = 1.0
         else
            do j = 1, 16
               if ( cd .gt.cosnads(j).and.cd.le.cosnads(j+1)) d1 = j
            enddo
            d2 = d1 + 1 ; fd1 = (cosnads(d2)-cd)/(cosnads(d2)-cosnads(d1)) ; fd2 = 1.0 - fd1
         endif
         foQ_1 = fd1*foQ_nads(s1,d1) + fd2*foQ_nads(s1,d2) 
         foQ_2 = fd1*foQ_nads(s2,d1) + fd2*foQ_nads(s2,d2) 
         FoQ_Int_SD(n,i) = dble ( fs1 * foQ_1 + fs2 * foQ_2 )
         foQ_1 = fd1*dC_foQ_nads(s1,d1) + fd2*dC_foQ_nads(s1,d2) 
         foQ_2 = fd1*dC_foQ_nads(s2,d1) + fd2*dC_foQ_nads(s2,d2) 
         dC_FoQ_Int_SD(n,i) = dble ( fs1 * foQ_1 + fs2 * foQ_2 )
      enddo

!  Viewing angles
!   -- R. Spurr, 12  April 2018.  Fixed glitch in f/Q interpolation

      do i = 1, nvzas
         cd = cvzas(i)
         if ( cd.le.cosnads(1) ) then
            d1 = 1 ; d2 = d1 + 1 ; fd1 = 1.0 ; fd2 = 0.0
         else if ( cd.ge.cosnads(17) ) then
!            d1 = 5 ; d2 = d1 + 1 ; fd1 = 0.0 ; fd2 = 1.0   !  Wrong offset, produced glitches
            d1 = 16 ; d2 = d1 + 1 ; fd1 = 0.0 ; fd2 = 1.0
         else
            do j = 1, 16
               if ( cd .gt.cosnads(j).and.cd.le.cosnads(j+1)) d1 = j
            enddo
            d2 = d1 + 1 ; fd1 = (cosnads(d2)-cd)/(cosnads(d2)-cosnads(d1)) ; fd2 = 1.0 - fd1
         endif
         foQ_1 = fd1*foQ_nads(s1,d1) + fd2*foQ_nads(s1,d2) 
         foQ_2 = fd1*foQ_nads(s2,d1) + fd2*foQ_nads(s2,d2) 
         FoQ_Int_SV(n,i) = dble ( fs1 * foQ_1 + fs2 * foQ_2 )
         foQ_1 = fd1*dC_foQ_nads(s1,d1) + fd2*dC_foQ_nads(s1,d2) 
         foQ_2 = fd1*dC_foQ_nads(s2,d1) + fd2*dC_foQ_nads(s2,d2) 
         dC_FoQ_Int_SV(n,i) = dble ( fs1 * foQ_1 + fs2 * foQ_2 )
      enddo

!  end solar loop

   enddo

!  Normal return

   return

!  Error return

88 continue
   fail = .true.
!mick fix 3/22/2017 - upgraded error msg to include file name
   message = 'Openfile error in Interpolate_Lin_fOQ_BS2; file not found: ' // Trim(FoQFile)

   return
end subroutine Interpolate_Lin_fOQ_BS2

!

subroutine Lin_Water_Transmittance &
    ( Max_PolarQuads, Max_AzimQuads,                          & ! Input
      PolarQuads, CosPolarQuads, SinPolarQuads, PolarWeights, & ! Input
      AzimQuads,  CosAzimQuads,  SinAzimQuads,  AzimWeights,  & ! Input
      DO_ISOTROPIC, DO_SHADOW, DO_COEFFS,                     & ! Input
      INCIDENT_ANGLE, REFRAC_R, REFRAC_I,                     & ! Input
      WINDSPEED, SUNGLINT_COEFFS, dSUNGLINT_COEFFS,           & ! Input
      PHI_W, CPHI_W, SPHI_W, TRANS_NORM,                      & ! Input
      TRANS, dTRANS )

      IMPLICIT NONE
      INTEGER  , parameter:: fpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ------

   integer, intent(in)   :: Max_PolarQuads, Max_AzimQuads

!  Quadratures

   real(fpk), intent(in) :: PolarQuads    (Max_PolarQuads)   ! Radians
   real(fpk), intent(in) :: CosPolarQuads (Max_PolarQuads)
   real(fpk), intent(in) :: SinPolarQuads (Max_PolarQuads)
   real(fpk), intent(in) :: PolarWeights  (Max_PolarQuads)

   real(fpk), intent(in) :: AzimQuads    (Max_AzimQuads)     ! Radians
   real(fpk), intent(in) :: CosAzimQuads (Max_AzimQuads)
   real(fpk), intent(in) :: SinAzimQuads (Max_AzimQuads)
   real(fpk), intent(in) :: AzimWeights  (Max_AzimQuads)

!  Windspeed, coefficients
!  -----------------------

!  Flag for Calculating Cox-Munk Coefficients
!     Only needs to be done once, so intent(inout)

      LOGICAL  , intent(inout) :: DO_COEFFS

!  Windspeed m/s

      REAL(fpk), intent(in)    :: WINDSPEED

!  Azimuth between Sun and Wind directions. angle in Radians + Cosine/sine

      REAL(fpk), intent(in)    :: PHI_W, CPHI_W, SPHI_W

!  Cox-Munk Coefficients. Intent(inout).

      REAL(fpk), intent(inout) :: SUNGLINT_COEFFS(7)
      REAL(fpk), intent(inout) :: dSUNGLINT_COEFFS(7)

!  Other inputs
!  ------------

!  Flag for using Isotropic Facet distribution

      LOGICAL  , intent(in)    :: DO_ISOTROPIC

!  Flag for including Shadow effect

      LOGICAL  , intent(in)    :: DO_SHADOW

!  incident angle in degrees

      REAL(fpk), intent(in)    :: INCIDENT_ANGLE

!  Real and imaginary parts of refractive index

      REAL(fpk), intent(in)    :: REFRAC_R
      REAL(fpk), intent(in)    :: REFRAC_I

!  Pre-computed Norm

      REAL(fpk), intent(in)    :: TRANS_NORM

!  Output
!  ======

      REAL(fpk), intent(out)   :: TRANS
      REAL(fpk), intent(out)   :: dTRANS

!  Local
!  =====

      integer   :: i, k
      real(fpk) :: dtr, xj, sxj, xi, sxi, phi, cphi, sphi, weight
      real(fpk) :: SUNGLINT_REFLEC, dSUNGLINT_REFLEC

!  Computation setup

      TRANS  = 0.0_fpk
      dTRANS = 0.0_fpk
      DTR   = ACOS(-1.0d0) / 180.0_fpk
      XJ  = COS ( INCIDENT_ANGLE * DTR )
      SXJ = SQRT ( 1.0_fpk - XJ * XJ )

!  Loops

      do k = 1, Max_AzimQuads
         PHI  = AzimQuads(K)/dtr
         CPHI = CosAzimQuads(K)
         SPHI = SinAzimQuads(K)
         do i = 1, Max_PolarQuads
            XI  = CosPolarQuads(I)
            SXI = SinPolarQuads(I)
            Call Lin_GENERAL_SUNGLINT &
             ( DO_ISOTROPIC, DO_SHADOW, DO_COEFFS,    &
               REFRAC_R, REFRAC_I, WINDSPEED,         &
               PHI_W, CPHI_W, SPHI_W,                 &
               XJ, SXJ, XI, SXI, PHI, CPHI, SPHI,     &
               SUNGLINT_COEFFS, DSUNGLINT_COEFFS,     &
               SUNGLINT_REFLEC, DSUNGLINT_REFLEC )
            WEIGHT = PolarWeights(I) * AzimWeights(k)
            TRANS  = TRANS  + SUNGLINT_REFLEC  * WEIGHT
            dTRANS = dTRANS + dSUNGLINT_REFLEC * WEIGHT
         enddo
      enddo
      dTRANS = -dTRANS/TRANS_NORM
      TRANS = 1.0_fpk - (TRANS/TRANS_NORM)

!  done

      RETURN
end subroutine Lin_Water_Transmittance

!

subroutine Lin_GENERAL_SUNGLINT &
         ( DO_ISOTROPIC, DO_SHADOW, DO_COEFFS,    &
           REFRAC_R, REFRAC_I, WINDSPEED,         &
           PHI_W, CPHI_W, SPHI_W,                 &
           XJ, SXJ, XI, SXI, PHI, CPHI, SPHI,     &
           SUNGLINT_COEFFS, DSUNGLINT_COEFFS,     &
           SUNGLINT_REFLEC, DSUNGLINT_REFLEC )

      implicit none

      integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  subroutine Input arguments
!  --------------------------

!  Flag for using Isotropic Facet distribution

      LOGICAL  , intent(in)    :: DO_ISOTROPIC

!  Flag for including Shadow effect

      LOGICAL  , intent(in)    :: DO_SHADOW

!  Flag for Calculating Cox-Munk Coefficients
!     Only needs to be done once, so intent(inout)

      LOGICAL  , intent(inout) :: DO_COEFFS

!  Real and imaginary parts of refractive index

      REAL(fpk), intent(in)    :: REFRAC_R
      REAL(fpk), intent(in)    :: REFRAC_I

!  Windspeed m/s

      REAL(fpk), intent(in)    :: WINDSPEED

!  Azimuth between Sun and Wind directions. angle in Radians + Cosine/sine

      REAL(fpk), intent(in)    :: PHI_W, CPHI_W, SPHI_W

!  Incident and reflected ddirections: sines/cosines. Relative azimuth (angle in radians)

      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SPHI

!  subroutine output arguments
!  ---------------------------

!   Glitter reflectance

      REAL(fpk), intent(out)   :: SUNGLINT_REFLEC
      REAL(fpk), intent(out)   :: dSUNGLINT_REFLEC

!  Cox-Munk Coefficients. Intent(inout).

      REAL(fpk), intent(inout) :: SUNGLINT_COEFFS(7)
      REAL(fpk), intent(inout) :: DSUNGLINT_COEFFS(7)

!  Local arguments
!  ---------------

!  parameters from LIDORT

   real(fpk), PARAMETER :: ONE = 1.0_fpk, ZERO  = 0.0_fpk, ONEP5 = 1.5_fpk
   real(fpk), PARAMETER :: TWO = 2.0_fpk, THREE = 3.0_fpk, FOUR  = 4.0_fpk
   real(fpk), PARAMETER :: six = two * three, twentyfour = six * four
   real(fpk), PARAMETER :: QUARTER = 0.25_fpk, HALF = 0.5_fpk

   real(fpk), PARAMETER :: MINUS_ONE = - ONE
   real(fpk), PARAMETER :: MINUS_TWO = - TWO

   real(fpk), PARAMETER :: PIE = ACOS(MINUS_ONE)
   real(fpk), PARAMETER :: DEG_TO_RAD = PIE/180.0_fpk

   real(fpk), PARAMETER :: PI2  = TWO  * PIE
   real(fpk), PARAMETER :: PI4  = FOUR * PIE
   real(fpk), PARAMETER :: PIO2 = HALF * PIE
   real(fpk), PARAMETER :: PIO4 = QUARTER * PIE

   REAL(fpk), PARAMETER   ::  CRITEXP = 88.0D0

!  Local variables

      REAL(fpk)  :: B, ZX, ZY, Z, Z1, Z2, XMP
      REAL(fpk)  :: TILT, TANTILT, TANTILT_SQ, COSTILT
      REAL(fpk)  :: ARGUMENT, PROB, FAC2, COEFF, VAR, WSigC, WSigU
      REAL(fpk)  :: XE, XN, XE_sq, XN_sq, XE_sq_1, XN_sq_1
      REAL(fpk)  :: XPHI, CKPHI, SKPHI, XPHI_W, CKPHI_W, SKPHI_W
      REAL(fpk)  :: S1, S2, S3, XXI, XXJ, T1, T2, DCOT
      REAL(fpk)  :: SHADOWI, SHADOWR, SHADOW

      REAL(fpk)  :: dARGUMENT, dPROB, dCOEFF, dVAR, dWSigC, dWSigU
      REAL(fpk)  :: AN, AE, dXE, dXN, dXE_sq, dXN_sq, EXPO, dEXPO
      REAL(fpk)  :: T3, T4, T5, T6, T7, dT3, dT4, dT5, dT6, dT7

!  Initialise output

      SUNGLINT_REFLEC = ZERO
      dSUNGLINT_REFLEC = ZERO

!  COmpute coefficients, according to 6S formulation

      IF ( DO_COEFFS ) THEN
         SUNGLINT_COEFFS = zero ; DSUNGLINT_COEFFS = zero
         IF ( DO_ISOTROPIC ) THEN
            SUNGLINT_COEFFS(1)  = 0.003_fpk + 0.00512_fpk * WINDSPEED
            DSUNGLINT_COEFFS(1) = 0.00512_fpk
         ELSE
            SUNGLINT_COEFFS(1) = 0.003_fpk + 0.00192_fpk * WINDSPEED ! sigmaC
            SUNGLINT_COEFFS(2) =             0.00316_fpk * WINDSPEED ! sigmaU
            SUNGLINT_COEFFS(3) = 0.010_fpk - 0.00860_fpk * WINDSPEED ! C21
            SUNGLINT_COEFFS(4) = 0.040_fpk - 0.03300_fpk * WINDSPEED ! C03
            SUNGLINT_COEFFS(5) = 0.400_fpk                           ! C40
            SUNGLINT_COEFFS(6) = 0.230_fpk                           ! C04
            SUNGLINT_COEFFS(7) = 0.120_fpk                           ! C22
            DSUNGLINT_COEFFS(1) = 0.00192_fpk
            DSUNGLINT_COEFFS(2) = 0.00316_fpk 
            DSUNGLINT_COEFFS(3) = - 0.00860_fpk
            DSUNGLINT_COEFFS(4) = - 0.03300_fpk
         ENDIF
         DO_COEFFS = .false.
      ENDIF

!  Local angles

      XPHI   = PIE - PHI       ! Not used
!     CKPHI  = - CPHI          ! Original, not correct.

      CKPHI  = + CPHI
      SKPHI  = + SPHI

      XPHI_W  = PHI_W
      CKPHI_W = CPHI_W
      SKPHI_W = SPHI_W

!  Tilt angle

      B  = ONE / ( XI + XJ )
      ZX = - SXI * SKPHI * B
      ZY = ( SXJ + SXI * CKPHI ) * B
      TANTILT_SQ = ZX * ZX + ZY * ZY
      TANTILT    = SQRT ( TANTILT_SQ )
      TILT       = ATAN(TANTILT)
      COSTILT    = COS(TILT)

!  Scatter angle

      Z = XI * XJ + SXI * SXJ * CKPHI
      IF ( Z .GT. ONE) Z = ONE
      Z1 = ACOS(Z)
      Z2 = COS(Z1*HALF)

!  Fresnel
!  -------

       CALL Fresnel_Complex ( REFRAC_R, REFRAC_I, Z2, XMP )

!  Anisotropic
!  -----------

      IF ( .not. DO_ISOTROPIC ) THEN

!  Variance

         WSigC = ONE / Sqrt(SUNGLINT_COEFFS(1)) ; dWSigC = - half * WSigC * WSigC * WSigC * DSUNGLINT_COEFFS(1)
         WSigU = ONE / Sqrt(SUNGLINT_COEFFS(2)) ; dWSigU = - half * WSigU * WSigU * WSigU * DSUNGLINT_COEFFS(2)
         VAR   = WSigC * WSigU * HALF ; dVAR = half * ( dWSigC * WSigU + WSigC * dWSigU )
         VAR = ONE / VAR ; dVAR =  - VAR * VAR * dVAR
  
!  angles

         AE = (  CKPHI_W * ZX + SKPHI_W * ZY )
         AN = ( -SKPHI_W * ZX + CKPHI_W * ZY )
         XE = AE * WSigC ; XE_sq = XE * XE ; XE_sq_1 = xe_sq - one ; dXE = AE * dWSigC ; dXE_sq = two * dXE * XE
         XN = AN * WSigU ; XN_sq = XN * XN ; XN_sq_1 = xn_sq - one ; dXN = AN * dWSigU ; dXN_sq = two * dXN * XN

!  GC Coefficient

         T3  = XE_sq_1 * XN * half
         dT3 = ( XE_sq_1 * dXN + dXE_sq * XN ) * half
         T4  = ( XN_sq - three ) * XN / six
         dT4 = ( ( XN_sq - three ) * dXN + dXN_sq * XN ) / six
         T5  = ( XE_sq * XE_sq - six * XE_sq + three ) / twentyfour
         dT5 = ( two * dXE_sq * XE_sq - six * dXE_sq ) / twentyfour
         T6  = ( XN_sq * XN_sq - six * XN_sq + three ) / twentyfour
         dT6 = ( two * dXN_sq * XN_sq - six * dXN_sq ) / twentyfour
         T7  = XE_sq_1 * XN_sq_1 / four
         dT7 = ( dXE_sq * XN_sq_1 + XE_sq_1 * dXN_sq ) / four

         Coeff  = ONE - SUNGLINT_COEFFS(3) * T3 &
                      - SUNGLINT_COEFFS(4) * T4 &
                      + SUNGLINT_COEFFS(5) * T5 &
                      + SUNGLINT_COEFFS(6) * T6 &
                      + SUNGLINT_COEFFS(7) * T7
         dCoeff =  - dSUNGLINT_COEFFS(3) * T3 - SUNGLINT_COEFFS(3) * dT3 &
                   - dSUNGLINT_COEFFS(4) * T4 - SUNGLINT_COEFFS(4) * dT4 &
                                              + SUNGLINT_COEFFS(5) * dT5 &
                                              + SUNGLINT_COEFFS(6) * dT6 &
                                              + SUNGLINT_COEFFS(7) * dT7

!  Probability and finish

         ARGUMENT  = (  XE_sq  +  XN_sq ) * HALF
         dARGUMENT = ( dXE_sq  + dXN_sq ) * HALF
         IF ( ARGUMENT .LT. CRITEXP ) THEN
            EXPO = EXP ( - ARGUMENT ) ; dEXPO = - dARGUMENT * EXPO
            PROB = COEFF * EXPO / VAR ; dPROB =  ( dCOEFF * EXPO + COEFF * dEXPO - PROB * dVAR ) / VAR
            FAC2 = QUARTER / XI / XJ / ( COSTILT ** FOUR )
            SUNGLINT_REFLEC  = XMP * PROB  * FAC2
            dSUNGLINT_REFLEC = XMP * dPROB * FAC2
         ENDIF

      ENDIF

!  Isotropic
!  ---------

      IF ( DO_ISOTROPIC ) THEN

!  Compute Probability and finish

         VAR   = SUNGLINT_COEFFS(1) ; dVAR = dSUNGLINT_COEFFS(1)
         ARGUMENT = TANTILT_SQ / VAR
         dARGUMENT = - ARGUMENT * dVAR / VAR
         IF ( ARGUMENT .LT. CRITEXP ) THEN
            EXPO = EXP ( - ARGUMENT ) ; dEXPO = - dARGUMENT * EXPO
            PROB = EXPO / VAR ; dPROB =  ( dEXPO - PROB * dVAR ) / VAR
            FAC2 = QUARTER / XI / XJ / ( COSTILT ** FOUR )
            SUNGLINT_REFLEC  = XMP * PROB  * FAC2
            dSUNGLINT_REFLEC = XMP * dPROB * FAC2
         ENDIF

      ENDIF

!  No Shadow code if not flagged

      IF ( .not. DO_SHADOW  ) RETURN

!  Shadow code

      S1 = SQRT ( VAR / PIE )
      S3 = ONE / ( SQRT(VAR) )
      S2 = S3 * S3

      XXI  = XI*XI
      DCOT = XI / SQRT ( ONE - XXI )
      T1   = EXP ( - DCOT * DCOT * S2 )
      T2   = DCOT * S3 ; CALL HOMEGROWN_ERRFUNC ( T2 )     !  Error function
      SHADOWI = HALF * ( S1 * T1 / DCOT - T2 )

      XXJ  = XJ*XJ
      DCOT = XJ / SQRT ( ONE - XXJ )
      T1   = EXP ( - DCOT * DCOT * S2 )
      T2   = DCOT * S3 ; CALL HOMEGROWN_ERRFUNC ( T2 )     !  Error function
      SHADOWR = HALF * ( S1 * T1 / DCOT - T2 )

      SHADOW = ONE / ( ONE + SHADOWI + SHADOWR )
      SUNGLINT_REFLEC = SUNGLINT_REFLEC * SHADOW

!     Finish

      RETURN
END subroutine Lin_GENERAL_SUNGLINT

!  End module

END MODULE vsleave_lin_sup_routines_2_m

