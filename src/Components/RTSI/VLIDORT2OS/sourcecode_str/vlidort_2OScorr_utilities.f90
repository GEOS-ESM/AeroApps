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

! ###########################################################
! #                                                         #
! #  2OS History :                                          #
! #                                                         #
! #     Mark 1: 2007 for OCO-Mk1, publication               #
! #     Mark 2: 2009, with BRDFs                            #
! #     Mark 3: 2013, Re-engineered Model:                  #
! #              * Including Multiple Geometry              #
! #              * Merging with Stand-alone FO code         #
! #              * New linearization for Bulk properties    #
! #     Mark 4: 2014, integration with VLIDORT              #
! #                                                         #
! #                                                         #
! ###########################################################

module vlidort_2OScorr_utilities

!  Purpose: Various utilities:-

!    a. ExpTrans         : computes Transmittance exp(-x)
!    b. Taylor_1, 2, 3   : 3 Taylor series expansions

!    c. L_Taylor_1       : L_Taylor series expansion
!    d. ExpTrans_L_1     : Transmittance exp(-x) + Linearization L[exp(-x)]
!    e. ExpTrans_L_2     : Transmittance exp(-bx) + Linearization L[exp(-bx)]

!    f. Make_Trans23     : Transposes Matrix S(MaxDims,4,4)
!    g. Make_Trans23_P   : Transposes Matrix L_S(MaxDims,4,4,MaxPars)

!    h. Gaussquad_2os    : Half-space Gaussian Quadrature.

   use VLIDORT_PARS, only : fpk, zero, one, two, three,   &
                            quarter, half, third, pie, Gaussian_eps

!  Everything public

public

contains

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!      TAYLOR-SERIES and EXPONENTIAL ROUTINES
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  The Exponential Transmittance function

Subroutine ExpTrans(BigValue,x,yfunc)
   implicit none
   real(fpk), intent(in)    :: BigValue,x
   real(fpk), intent(inout) :: yfunc
   yfunc = zero ; if ( x .lt. BigValue) yfunc = exp(-x)
   return
end Subroutine ExpTrans

!  Three Taylor series subroutine

Subroutine Taylor_1(eps,x,yfunc)
   implicit none
   real(fpk), intent(in)    :: eps, x
   real(fpk), intent(inout) :: yfunc
   real(fpk) :: xx
   xx = eps * x ; yfunc = one - half * xx * ( one - third * xx )
   return
end Subroutine Taylor_1

Subroutine Taylor_3(eps,x1,x2,x3,yfunc)
   implicit none
   real(fpk), intent(in)    :: eps, x1, x2, x3
   real(fpk), intent(inout) :: yfunc
   real(fpk) :: xx1, xx2, t1, t2, t3
   xx1 = eps * x1 ; xx2 = eps * x2
   t1 =    x3    * x1 * ( one - half * xx1 * ( one - third * xx1 ) )
   t2 = (one-x3) * x2 *  ( one - xx2 * ( one - xx2 ) )
   t3 =    x3    * x2 * xx1 * ( one - half * xx1 - xx2 )
   yfunc = - t1 + t2 + t3
   return
end Subroutine Taylor_3

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!            LINEARIZED ROUTINES
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

Subroutine LS_Taylor_1(x1,npars,L_x1,L_yfunc)
   implicit none
   integer  , intent(in)    :: npars
   real(fpk), intent(in)    :: x1,L_x1(npars)
   real(fpk), intent(inout) :: L_yfunc(npars)
   real(fpk)                :: xhelp
   integer                  :: p
   xhelp  = one - x1 + half * x1 * x1
   do p = 1, npars
      L_yfunc(p) = L_x1(p) * xhelp
   enddo
end Subroutine LS_Taylor_1

Subroutine L_Taylor_1(x1,x2,npars,L_x1,L_x2,L_yfunc)
   implicit none
   integer  , intent(in)    :: npars
   real(fpk), intent(in)    :: x1,x2,L_x1(npars),L_x2(npars)
   real(fpk), intent(inout) :: L_yfunc(npars)
   real(fpk)                :: xhelp, xhelp3, L_xhelp, L_xhelp3
   integer                  :: p
   xhelp3 = one - third * x2
   xhelp  = one - half  * x2 * xhelp3
   do p = 1, npars
      L_xhelp3   = - third * L_x2(p)
      L_xhelp    = - ( L_x2(p) * xhelp3 + x2 * L_xhelp3 ) * half
      L_yfunc(p) = L_x1(p) * xhelp + x1 * L_xhelp
   enddo
end Subroutine L_Taylor_1

Subroutine L_Taylor_2(d,x1,x2,x3,x4,npars,L_d,L_x1,L_x2,L_x3,L_x4,L_yfunc)
   implicit none
   integer  , intent(in)    :: npars
   real(fpk), intent(in)    :: d,x1,x2,x3,x4, L_d(npars)
   real(fpk), intent(in)    :: L_x1(npars),L_x2(npars),L_x3(npars),L_x4(npars)
   real(fpk), intent(inout) :: L_yfunc(npars)
   integer                  :: p
   real(fpk)                :: H12,H34,HAD,HMD,HND,L34,L12,L34d,LMD
   H12 = ( x1 + two * x2 ) * half
   H34 = ( x3 * x3 + x3 * x4 + x4 * x4 ) * third
   HAD = H12 - d * H34 ; HMD = d * HAD ; HND = one - HMD
   do p = 1, npars
      L34 = third*(two*L_x3(p)*x3+L_x3(p)*x4+x3*L_x4(p)+two*L_x4(p)*x4)
      L12 = half *( L_x1(p) + two * L_x2(p) )
      L34d =  L_d(p) * H34 * d * L34
      LMD  = d * (L12 - L34d) + L_d(p) * HAD 
      L_yfunc(p)  = - d * LMD  + L_d(p) * HND 
   enddo
   return
end Subroutine L_Taylor_2

Subroutine L_Taylor_3(z,x1,x2,x3,npars,L_z,L_x1,L_x2,L_x3,L_yfunc)
   implicit none
   integer  , intent(in)    :: npars
   real(fpk), intent(in)    :: z,x1,x2,x3,L_z(npars)
   real(fpk), intent(in)    :: L_x1(npars),L_x2(npars),L_x3(npars)
   real(fpk), intent(inout) :: L_yfunc(npars)
   integer                  :: p
   real(fpk)                :: H123,L123
   H123 = ( x1 + x2 + two * x3 ) * third
   do p = 1, npars
     L123 = ( L_x1(p) + L_x2(p) + two * L_x3(p) ) * third
     L_yfunc(p) = L_z(p) * ( one - H123 ) - z * L123
   enddo
   return
end Subroutine L_Taylor_3

subroutine ExpTrans_L(BigValue,npars,x,L_x,yfunc,L_yfunc)
!  yfunc assumed known
   implicit none
   integer  , intent(in)    :: npars
   real(fpk), intent(in)    :: BigValue,x,L_x(npars)
   real(fpk), intent(inout) :: yfunc,L_yfunc(npars)
   L_yfunc = zero
   if ( x .lt. BigValue) then
      L_yfunc(1:npars) = - L_x(1:npars) * yfunc
   endif
   return
end subroutine ExpTrans_L

Subroutine ExpTrans_L_2(BigValue,b,x,L_b,L_x,yfunc,L_yfunc)
   implicit none
   real(fpk), intent(in)    :: BigValue,b,x,L_b,L_x
   real(fpk), intent(inout) :: yfunc,L_yfunc
   real(fpk)                :: bx
   L_yfunc = zero ; bx = b * x
   if ( bx .lt. BigValue) then
      L_yfunc = - ( L_b*x + L_x*b ) * yfunc
   endif
   return
end Subroutine ExpTrans_L_2

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!         Transformation (index switching)
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine Make_Trans23(MaxDim,n,S)
   implicit none
   integer  , intent(in)    :: MaxDim,n
   real(fpk), intent(inout) :: S(MaxDim,4,4)
   integer   :: i,j
   real(fpk) :: s1,s2

   do i = 2,2
      do j = 1,4
         s1 = S(n,i,j)
         s2 = S(n,i+1,j)
         S(n,i,j) = s1+s2
         S(n,i+1,j) = s1-s2
      enddo
      do j = 1,4
         s1 = S(n,j,i)
         s2 = S(n,j,i+1)
         S(n,j,i) = s1+s2
         S(n,j,i+1) = s1-s2
      enddo
    enddo

   return
end subroutine Make_Trans23

subroutine Make_Trans23_P(MaxDim,MaxPars,n,p,S)
   implicit none
   integer  , intent(in)    :: MaxDim,MaxPars,n,p
   real(fpk), intent(inout) :: S(MaxDim,4,4,MaxPars)
   integer   :: i,j
   real(fpk) :: s1,s2

   do i = 2,2
      do j = 1,4
         s1 = S(n,i,j,p)
         s2 = S(n,i+1,j,p)
         S(n,i,j,p) = s1+s2
         S(n,i+1,j,p) = s1-s2
      enddo
      do j = 1,4
         s1 = S(n,j,i,p)
         s2 = S(n,j,i+1,p)
         S(n,j,i,p) = s1+s2
         S(n,j,i+1,p) = s1-s2
      enddo
    enddo

   return
end subroutine Make_Trans23_P

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!         Gaussian Quadrature routine
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine Gaussquad_2os (Maxstreams,nstreams,a,b,x,w)

   implicit none

!  Arguments

   integer  , intent(in)  :: Maxstreams,nstreams
   real(fpk), intent(in)  :: a,b
   real(fpk), intent(out) :: x(Maxstreams),w(Maxstreams)

!  local variables

   integer   :: m,i,j,ii
   logical   :: cond
   real(fpk) :: xm,xl,z,P1,P2,P3,Pa,z1,ri,rj,rj1,rj21,rs

   x = zero ; w = zero

   m = (nstreams+1)/2
   rs = real(nstreams,fpk)

   xm = half*(a+b)
   xl = half*(b-a)
   do i = 1, m
      ii = nstreams + 1 - i
      ri = real(i,fpk)
      z = cos(pie*(ri-quarter)/(rs+half))
      cond = .true.
      do while (cond)
         P1 = one ; P2 = zero
         do j = 1, nstreams
            P3 = P2 ; P2 = P1
            rj = real(j,fpk) ; rj1 = rj-one ;rj21 = rj + rj1
            P1 = ( rj21 * z * P2 - rj1 * P3 ) / rj
         enddo
         Pa = rs*(z*P1-P2)/(z*z-one)
         z1 = z ; z  = z1-P1/Pa
         if ( abs(z-z1) .le. Gaussian_eps) cond = .false.
      enddo
      x(i) = xm-xl*z                    ;  x(ii) = xm+xl*z
      w(i) = two*xl / ((one-z*z)*Pa*Pa) ;  w(ii) = w(i)
   enddo

   return
end subroutine Gaussquad_2os

!  FInish module

end module vlidort_2OScorr_utilities

