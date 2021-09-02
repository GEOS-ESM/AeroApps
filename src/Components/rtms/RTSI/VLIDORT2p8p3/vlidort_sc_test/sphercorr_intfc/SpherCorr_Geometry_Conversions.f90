
! ###############################################################
! #                                                             #
! #                       VLIDORT_2p8p3                         #
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
! #  This Version :   VLIDORT_2p8p3                             #
! #  Release Date :   31 March 2021                             #
! #                                                             #
! #  Previous VLIDORT Versions under Standard GPL 3.0:          #
! #  ------------------------------------------------           #
! #                                                             #
! #      2.7   F90, released        August 2014                 #
! #      2.8   F90, released        May    2017                 #
! #      2.8.1 F90, released        August 2019                 # 
! #      2.8.2 F90, limited release May    2020                 # 
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
! #      Doublet geometry post-processing    (2.8.2)            #
! #      Reduction zeroing, dynamic memory   (2.8.2)            #
! #                                                             #
! #  Features Summary of This VLIDORT Version                   #
! #  ----------------------------------------                   #
! #                                                             #
! #   2.8.3, released 31 March 2021.                            #
! #     ==> Green's function RT solutions (Nstokes = 1 or 3)    #
! #     ==> Sphericity Corrections using MS source terms        #
! #     ==> BRDF upgrades, including new snow reflectance       #
! #     ==> SLEAVE Upgrades, extended water-leaving treatment   #
! #                                                             #
! ###############################################################

! ###################################################################
! #                                                                 #
! # This is Version 2.8.3 of the VLIDORT_2p8 software library.      #
! # This library comes with the Standard GNU General Public License,#
! # Version 3.0, 29 June 2007. Please read this license carefully.  #
! #                                                                 #
! #      VLIDORT Copyright (c) 2003-2021.                           #
! #          Robert Spurr, RT Solutions, Inc.                       #
! #          9 Channing Street, Cambridge, MA 02138, USA.           #
! #                                                                 #
! # This file is part of VLIDORT_2p8p3 ( Version 2.8.3 )            #
! #                                                                 #
! # VLIDORT_2p8p3 is free software: you can redistribute it         #
! # and/or modify it under the terms of the Standard GNU GPL        #
! # (General Public License) as published by the Free Software      #
! # Foundation, either version 3.0 of the License, or any           #
! # later version.                                                  #
! #                                                                 #
! # VLIDORT_2p8p3 is distributed in the hope that it will be        #
! # useful, but WITHOUT ANY WARRANTY; without even the implied      #
! # warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR         #
! # PURPOSE. See the Standard GNU General Public License (GPL)      #
! # for more details.                                               #
! #                                                                 #
! # You should have received a copy of the Standard GNU General     #
! # Public License (GPL) Version 3.0, along with the VLIDORT_2p8p3  #
! # code package. If not, see <http://www.gnu.org/licenses/>.       #
! #                                                                 #
! ###################################################################

MODULE SPHERCORR_GEOMETRY_CONVERSIONS_m

!  HERE IS THE HISTORY
!  ===================

!  2-point correction history --
!    V1: 11/20/19. Originally Coded 20-29 November 2019 for VLIDORT Version 2.8.1
!    V2: 12/18/19. Extension to include BOA downwelling situation as well as TOA upwelling.
!    V3: 03/01/20. Renamed, along with addition of 3-point and multi-point Corrections

!  3-point correction history --
!    V1: 03/01/20. New 3pt Correction

!  Multi-point correction history --
!    V1: 03/01/20. New  Correction

!  1/31/21. Developed for Official VLIDORT Version 2.8.3 Package
!     R. Spurr. RT Solutions Inc.

!  2/28/21. Developed for Official LIDORT Version 3.8.3 Package
!     R. Spurr. RT Solutions Inc.

!  THIS MODULE - geometry conversion routines (all stand-alone)
!       BOATOA_2point_conversion
!       BOATOA_3point_conversion
!       BOATOA_Mpoint_conversion

! Everything public

public

contains


subroutine BOATOA_2point_conversion ( DirIdx, &
     Rearth, Hatmos, theta_boa, alpha_boa, phi_boa, & ! input TOA/BOA heights, BOA Geometry
     cosscat, theta_toa, alpha_toa, phi_toa )         ! output TOA geometry, scattering angle

!  Completely stand-alone 2-point conversion routine.
!   12/18/19. Add direction index for downwelling case

   implicit none

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs/Outputs

   integer  , intent(in)  :: DirIdx
   real(fpk), intent(in)  :: Rearth, Hatmos, theta_boa, alpha_boa, phi_boa
   real(fpk), intent(out) :: theta_toa, alpha_toa, phi_toa, cosscat

!  Local

   real(fpk) :: alpha_boa_R, theta_boa_R, phi_boa_R, dtr, pie, zero, one
   real(fpk) :: SolarDirection(3), vsign, vsign_solar, b, px(3)
   real(fpk) :: salpha_boa, calpha_boa, stheta_boa, ctheta_boa, cphi_boa, sphi_boa
   real(fpk) :: salpha_toa, calpha_toa, stheta_toa, ctheta_toa, cphi_toa
   real(fpk) :: term1, term2, rtoa, cumangle
   real(fpk) :: alpha_toa_R, theta_toa_R, phi_toa_R

!  constants

   zero = 0.0_fpk ; one = 1.0_fpk
   Pie = acos(-1.0_fpk)
   dtr = Pie / 180.0_fpk

!  BOA angles in radians (_R) and cos/sin

   alpha_boa_R = alpha_boa * dtr
   salpha_boa  = sin(alpha_boa_R)
   calpha_boa  = cos(alpha_boa_R)

   theta_boa_R    = theta_boa * dtr
   stheta_boa     = sin(theta_boa_R)
   ctheta_boa     = cos(theta_boa_R)
      
   phi_boa_R   = phi_boa * dtr
   cphi_boa    = cos(phi_boa_R)
   sphi_boa    = sin(phi_boa_R)

!  Solar direction (unit vector)

   vsign       = 1.0_fpk ; if ( DirIdx .eq. 2 ) vsign = -1.0_fpk
   vsign_solar = vsign ! vsign_solar = 1.0_fpk
   SolarDirection(1) = - stheta_boa * cphi_boa * vsign_solar
   SolarDirection(2) = - stheta_boa * sphi_boa
   SolarDirection(3) = - ctheta_boa

!  scattering angle

   term1 = salpha_boa * stheta_boa * cphi_boa
   term2 = calpha_boa * ctheta_boa
   cosscat = -vsign * term2 + term1 

!  TOA calculation
!  ===============

!  Vza at TOA (sine rule)

   rtoa = rearth+hatmos
   salpha_toa = rearth * salpha_boa / rtoa
   calpha_toa = sqrt(one-salpha_toa*salpha_toa)
   alpha_toa_R = asin(salpha_toa)
   cumangle = alpha_boa_R - alpha_toa_R

!  Calculate SZA at TOA, using vector geometry

   px(1) = - rtoa * sin(CumAngle)
   px(2) = 0.0_fpk
   px(3) =   rtoa * cos(CumAngle)
   b = DOT_PRODUCT(px,SolarDirection)
   ctheta_toa  = - b/rtoa
   stheta_toa  = sqrt(one-ctheta_toa*ctheta_toa)
   theta_toa_R = acos(ctheta_toa)

!  Calculate PHI at TOA using constancy of scattering angle

   cphi_toa = (cosscat+vsign*calpha_toa*ctheta_toa)/stheta_toa/salpha_toa
   if ( cphi_toa.gt.one)  cphi_toa = one
   if ( cphi_toa.lt.-one) cphi_toa = -one
   phi_toa_R = acos(cphi_toa)
   if ( phi_boa.gt.180.0_fpk) phi_toa_R = 2.0_fpk * Pie - phi_toa_R

!  Final answers for TOA geometry (in degrees)

   alpha_toa = alpha_toa_R / dtr
   theta_toa = theta_toa_R / dtr
   phi_toa   = phi_toa_R   / dtr

!  finish

   return
end subroutine BOATOA_2point_conversion

!

subroutine BOATOA_3point_conversion ( DirIdx, &
     Rearth, Hatmos, BOA_geometry,                & ! input TOA height, BOA Geometry
     Hmid, cosscat, TOA_geometry, MID_geometry )    ! output TOA/MID geometries, scattering angle, Midheight

!  3/1/20. Completely new, stand-alone 3-point geometry conversion routine.

   implicit none

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs/Outputs

   integer  , intent(in)  :: DirIdx
   real(fpk), intent(in)  :: Rearth, Hatmos, BOA_geometry(3)
   real(fpk), intent(out) :: TOA_geometry(3), MID_geometry(3), cosscat, Hmid

!  Local

   real(fpk) :: alpha_boa, theta_boa, phi_boa
   real(fpk) :: alpha_boa_R, theta_boa_R, phi_boa_R, dtr, pie, zero, one
   real(fpk) :: SolarDirection(3), vsign, vsign_solar, b, px(3), term1, term2, cumangle
   real(fpk) :: salpha_boa, calpha_boa, stheta_boa, ctheta_boa, cphi_boa, sphi_boa
   real(fpk) :: salpha_toa, calpha_toa, stheta_toa, ctheta_toa, cphi_toa
   real(fpk) :: rtoa, alpha_toa_R, theta_toa_R, phi_toa_R
   real(fpk) :: salpha_mid, calpha_mid, stheta_mid, ctheta_mid, cphi_mid
   real(fpk) :: rmid, alpha_mid_R, theta_mid_R, phi_mid_R

!  constants

   zero = 0.0_fpk ; one = 1.0_fpk
   Pie = acos(-1.0_fpk)
   dtr = Pie / 180.0_fpk

!  BOA angles in radians (_R) and cos/sin

   theta_boa = BOA_GEOMETRY(1)
   alpha_boa = BOA_GEOMETRY(2)
   phi_boa   = BOA_GEOMETRY(3)

   alpha_boa_R = alpha_boa * dtr
   salpha_boa  = sin(alpha_boa_R)
   calpha_boa  = cos(alpha_boa_R)

   theta_boa_R    = theta_boa * dtr
   stheta_boa     = sin(theta_boa_R)
   ctheta_boa     = cos(theta_boa_R)
      
   phi_boa_R   = phi_boa * dtr
   cphi_boa    = cos(phi_boa_R)
   sphi_boa    = sin(phi_boa_R)

!  Solar direction (unit vector)

   vsign       = 1.0_fpk ; if ( DirIdx .eq. 2 ) vsign = -1.0_fpk
   vsign_solar = vsign ! vsign_solar = 1.0_fpk
   SolarDirection(1) = - stheta_boa * cphi_boa * vsign_solar
   SolarDirection(2) = - stheta_boa * sphi_boa
   SolarDirection(3) = - ctheta_boa

!  scattering angle

   term1 = salpha_boa * stheta_boa * cphi_boa
   term2 = calpha_boa * ctheta_boa
   cosscat = -vsign * term2 + term1 

!  TOA calculation
!  ===============

!  Vza at TOA (sine rule)

   rtoa = rearth+hatmos
   salpha_toa = rearth * salpha_boa / rtoa
   calpha_toa = sqrt(one-salpha_toa*salpha_toa)
   alpha_toa_R = asin(salpha_toa)
   cumangle = alpha_boa_R - alpha_toa_R

!  Calculate SZA at TOA, using vector geometry

   px(1) = - rtoa * sin(CumAngle)
   px(2) = 0.0_fpk
   px(3) =   rtoa * cos(CumAngle)
   b = DOT_PRODUCT(px,SolarDirection)
   ctheta_toa  = - b/rtoa
   stheta_toa  = sqrt(one-ctheta_toa*ctheta_toa)
   theta_toa_R = acos(ctheta_toa)

!  Calculate PHI at TOA using constancy of scattering angle

   cphi_toa = (cosscat+vsign*calpha_toa*ctheta_toa)/stheta_toa/salpha_toa
   if ( cphi_toa.gt.one)  cphi_toa = one
   if ( cphi_toa.lt.-one) cphi_toa = -one
   phi_toa_R = acos(cphi_toa)
   if ( phi_boa.gt.180.0_fpk) phi_toa_R = 2.0_fpk * Pie - phi_toa_R

!  Final answers for TOA geometry (in degrees)

   TOA_GEOMETRY(1) = theta_toa_R / dtr
   TOA_GEOMETRY(2) = alpha_toa_R / dtr
   TOA_GEOMETRY(3) = phi_toa_R   / dtr

!  MID_HEIGHT calculation
!  ======================

!  determine mid height by averaging cos(VZA) of the TOA and BOA values.

   calpha_mid = ( calpha_toa + calpha_boa ) * 0.5_fpk
   salpha_mid = sqrt ( 1.0_fpk - calpha_mid * calpha_mid )

!  Use sine rule to get the mid-height value

   rmid  = rearth * salpha_boa / salpha_mid
   hmid  = rmid - rearth

!  VZA and earth-centered angles at Mid-height

   alpha_mid_R = asin(salpha_mid)
   cumangle = alpha_boa_R - alpha_mid_R

!  Calculate SZA at mid-height, using via vector geometry

   px(1) = - rmid * sin(CumAngle)
   px(2) = 0.0_fpk
   px(3) =   rmid * cos(CumAngle)
   b = DOT_PRODUCT(px,SolarDirection)
   ctheta_mid  = - b/rmid
   stheta_mid  = sqrt(one-ctheta_mid*ctheta_mid)
   theta_mid_R = acos(ctheta_mid)

!  Calculate PHI at mid-height using constancy of scattering angle

   cphi_mid = (cosscat+vsign*calpha_mid*ctheta_mid)/stheta_mid/salpha_mid
   if ( cphi_mid.gt.one)  cphi_mid = one
   if ( cphi_mid.lt.-one) cphi_mid = -one
   phi_mid_R = acos(cphi_mid)
   if ( phi_boa.gt.180.0_fpk) phi_mid_R = 2.0_fpk * Pie - phi_mid_R

!  Final answers for mid-height geometry (in degrees)

   MID_GEOMETRY(1) = theta_mid_R / dtr
   MID_GEOMETRY(2) = alpha_mid_R / dtr
   MID_GEOMETRY(3) = phi_mid_R   / dtr

!  finish

   return
end subroutine BOATOA_3point_conversion

!

subroutine BOATOA_Mpoint_conversion ( &
     maxgeoms, maxlayers, nlayers, ngeoms, DirIdx, & ! Input numbers
     Rearth, Heights, BOA_geometry,                & ! input Heights, BOA Geometry
     cosscat, ALL_geometry )                         ! output all geometries, Cosine scattering angle

!  3/2/20. Completely new, stand-alone Multi-point geometry conversion routine.

   implicit none

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs/Outputs

   integer  , intent(in)  :: DirIdx, maxgeoms, maxlayers, nlayers, ngeoms
   real(fpk), intent(in)  :: Rearth, Heights(0:maxlayers), BOA_geometry(3)
   real(fpk), intent(out) :: ALL_geometry(maxgeoms,3), cosscat

!  Local

   integer   :: n, ng
   real(fpk) :: alpha_boa, theta_boa, phi_boa
   real(fpk) :: alpha_boa_R, theta_boa_R, phi_boa_R, dtr, pie, zero, one
   real(fpk) :: SolarDirection(3), vsign, vsign_solar, b, px(3), term1, term2, cumangle
   real(fpk) :: salpha_boa, calpha_boa, stheta_boa, ctheta_boa, cphi_boa, sphi_boa
   real(fpk) :: salpha, calpha, stheta, ctheta, cphi
   real(fpk) :: rad, alpha_R, theta_R, phi_R

!   real(fpk) :: stheta_mid, ctheta_mid, cphi_mid
!   real(fpk) :: theta_mid_R, phi_mid_R
!   real(fpk) :: alpha_mid_R, calpha_mid, salpha_mid, rmid

!  constants

   zero = 0.0_fpk ; one = 1.0_fpk
   Pie = acos(-1.0_fpk)
   dtr = Pie / 180.0_fpk

!  BOA angles in radians (_R) and cos/sin

   All_Geometry(1,1:3) = BOA_Geometry(1:3)
   theta_boa = BOA_GEOMETRY(1)
   alpha_boa = BOA_GEOMETRY(2)
   phi_boa   = BOA_GEOMETRY(3)

   alpha_boa_R = alpha_boa * dtr
   salpha_boa  = sin(alpha_boa_R)
   calpha_boa  = cos(alpha_boa_R)

   theta_boa_R    = theta_boa * dtr
   stheta_boa     = sin(theta_boa_R)
   ctheta_boa     = cos(theta_boa_R)
      
   phi_boa_R   = phi_boa * dtr
   cphi_boa    = cos(phi_boa_R)
   sphi_boa    = sin(phi_boa_R)

!  Solar direction (unit vector)

   vsign       = 1.0_fpk ; if ( DirIdx .eq. 2 ) vsign = -1.0_fpk
   vsign_solar = vsign ! vsign_solar = 1.0_fpk
   SolarDirection(1) = - stheta_boa * cphi_boa * vsign_solar
   SolarDirection(2) = - stheta_boa * sphi_boa
   SolarDirection(3) = - ctheta_boa

!  scattering angle

   term1 = salpha_boa * stheta_boa * cphi_boa
   term2 = calpha_boa * ctheta_boa
   cosscat = -vsign * term2 + term1 

!  calculation
!  ===========

   do ng = 2, ngeoms
     n = ngeoms - ng

!  Vza at boundary-level height (sine rule)

     rad = rearth+heights(n)
     salpha = rearth * salpha_boa / rad
     calpha = sqrt(one-salpha*salpha)
     alpha_R = asin(salpha)
     cumangle = alpha_boa_R - alpha_R

!  Calculate SZA at this boundary, using vector geometry

     px(1) = - rad * sin(CumAngle)
     px(2) = 0.0_fpk
     px(3) =   rad * cos(CumAngle)
     b = DOT_PRODUCT(px,SolarDirection)
     ctheta  = - b/rad
     stheta  = sqrt(one-ctheta*ctheta)
     theta_R = acos(ctheta)

!  Calculate PHI using constancy of scattering angle

     cphi = (cosscat+vsign*calpha*ctheta)/stheta/salpha
     if ( cphi.gt.one)  cphi = one
     if ( cphi.lt.-one) cphi = -one
     phi_R = acos(cphi)
     if ( phi_boa.gt.180.0_fpk) phi_R = 2.0_fpk * Pie - phi_R

!  Final answers for TOA geometry (in degrees)

     ALL_GEOMETRY(ng,1) = theta_R / dtr
     ALL_GEOMETRY(ng,2) = alpha_R / dtr
     ALL_GEOMETRY(ng,3) = phi_R   / dtr

   enddo

!  finish

   return
end subroutine BOATOA_Mpoint_conversion

!  End module

END MODULE SPHERCORR_GEOMETRY_CONVERSIONS_m

