
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
! #    Taylor series, small number expansions:                  #
! #             TAYLOR_SERIES_1                                 #
! #             TAYLOR_SERIES_L_1                               #
! #                                                             #
! ###############################################################

!  4/15/20. Version 2.8.2. No changes

MODULE vlidort_Taylor_m

   USE VLIDORT_PARS_m, only : max_Taylor_terms, Taylor_small, zero, one, two


PUBLIC

CONTAINS

subroutine Taylor_series_1 &
          ( order, eps, delta, udel, sm, mult )

!  Good for the particular and homogeneous-solution multipliers.
!      Small number expansion to any order up to 4.

!  Note: In LIDORT, this subroutine is applied to quantities of the form
!                [exp(-b*DELTA) - exp(-a*DELTA)]
!            q = -------------------------------
!                            (a-b)
!        where (a-b) has become small to the point of causing instability.
!        Note the positions of "a" and "b" in the numerator are swapped
!        from their positions in the denominator.
!
!        Using the above form for "q", the I/O for the subroutine is as
!        follows:
!        * EPS   = (a-b)
!        * DELTA = optical thickness (whole or partial layer)
!        * TERM2 = exp(-a*DELTA) (usually UDEL or WDEL)
!        * MULT  = q

   IMPLICIT NONE

!  arguments

   INTEGER         , INTENT(IN)  :: order
   double precision, INTENT(IN)  :: eps, delta, udel, sm
   double precision, INTENT(OUT) :: mult

!  local declarations

   integer          :: mterms, m
   double precision :: power, d(0:max_Taylor_terms)

!  exp(De) expansion coefficients

   mterms = order + 1 ; d(0) = one
   do m = 1, mterms
      d(m) = delta * d(m-1) / dble(m)
   enddo

!  evaluate multiplier

   mult = d(1) ; power = one
   do m = 2, mterms
      power = power * eps ; mult = mult + power * d(m)
   enddo
   mult = mult * udel * sm

!  Equivalent to the following, for order = 3 (highest power of eps)
!      power   = delta*eps ;  power2  = power*power ; power3 = power2 * power
!      mult = udel * sm * delta *(one + half*power + power2/6.0_fpk + power3/24.0_fpk)

!  Finish

   return
end subroutine Taylor_Series_1

!

subroutine Taylor_series_L_1 ( order, eps, delta, ddot, kdot, Ldot, uterm, sm, L_mult )

!  Small number expansion for derivatives of Series 1 quantities
!     sm is required, but result is NOT SCALED

!   L_HMULT --> sm = user-secant   ,    Ldot = zero
!   L_EMULT --> sm = user_secant   ,    Ldot = zero

   implicit none

!  arguments

   INTEGER         , INTENT(IN)  :: order
   double precision, INTENT(IN)  :: eps, delta, ddot, kdot, Ldot, uterm, sm
   double precision, INTENT(OUT) :: L_mult

!  local declarations

   integer          :: mterms, m, m1, m2
   double precision :: power, d(0:max_Taylor_terms), series1, series2, series3

!  exp(De) expansion coefficients

   mterms = order + 2 ; d(0) = one
   do m = 1, mterms
      d(m) = delta * d(m-1) / dble(m)
   enddo

!  Develop  series
!  Series 3 absent for HMULT/EMULT, only present for GMULT

   power = one
   series1 = d(0) - d(1)*sm
   series2 = d(2) - d(1)*delta
   series3 = d(2)
   do m = 1, order
      m1 = m + 1 ; m2 = m1 + 1
      power = power * eps
      series1 = series1 + power * (d(m)  - d(m1)*sm)
      series2 = series2 + power * (d(m2) - d(m1)*delta)
      series3 = series3 + power * d(m2)
   enddo

!  final

   L_mult = ( ddot*series1 - Ldot*series3 + kdot*series2 ) * uterm

   return
end subroutine Taylor_series_L_1


!  Finish module

END MODULE vlidort_Taylor_m
