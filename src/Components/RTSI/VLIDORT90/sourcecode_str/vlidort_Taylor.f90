! ###############################################################
! #                                                             #
! #                    THE VECTOR LIDORT MODEL                  #
! #                                                             #
! #  (Vector LInearized Discrete Ordinate Radiative Transfer)   #
! #   -      --         -        -        -         -           #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Author :      Robert. J. D. Spurr                          #
! #                                                             #
! #  Address :     RT Solutions, inc.                           #
! #                9 Channing Street                            #
! #                Cambridge, MA 02138, USA                     #
! #                Tel: (617) 492 1183                          #
! #                                                             #
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  Versions     :   2.0, 2.2, 2.3, 2.4, 2.4R, 2.4RT, 2.4RTC,  #
! #                   2.5, 2.6, 2.7                             #
! #  Release Date :   December 2005  (2.0)                      #
! #  Release Date :   March 2007     (2.2)                      #
! #  Release Date :   October 2007   (2.3)                      #
! #  Release Date :   December 2008  (2.4)                      #
! #  Release Date :   April 2009     (2.4R)                     #
! #  Release Date :   July 2009      (2.4RT)                    #
! #  Release Date :   October 2010   (2.4RTC)                   #
! #  Release Date :   March 2011     (2.5)                      #
! #  Release Date :   May 2012       (2.6)                      #
! #  Release Date :   August 2014    (2.7)                      #
! #                                                             #
! #       NEW: TOTAL COLUMN JACOBIANS         (2.4)             #
! #       NEW: BPDF Land-surface KERNELS      (2.4R)            #
! #       NEW: Thermal Emission Treatment     (2.4RT)           #
! #       Consolidated BRDF treatment         (2.4RTC)          #
! #       f77/f90 Release                     (2.5)             #
! #       External SS / New I/O Structures    (2.6)             #
! #                                                             #
! #       SURFACE-LEAVING / BRDF-SCALING      (2.7)             #
! #       TAYLOR Series / OMP THREADSAFE      (2.7)             #
! #                                                             #
! ###############################################################

!    #####################################################
!    #                                                   #
!    #   This Version of VLIDORT comes with a GNU-style  #
!    #   license. Please read the license carefully.     #
!    #                                                   #
!    #####################################################

! ###############################################################
! #                                                             #
! #    Taylor series, small number expansions:                  #
! #             TAYLOR_SERIES_1                                 #
! #             TAYLOR_SERIES_L_1                               #
! #                                                             #
! ###############################################################

MODULE vlidort_Taylor_m

   USE VLIDORT_PARS, only : max_Taylor_terms, Taylor_small, zero, one, two


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
