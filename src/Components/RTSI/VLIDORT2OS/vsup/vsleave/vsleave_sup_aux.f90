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
! # Auxiliary subroutines for Surface-Leaving                   #
! #                                                             #
! #               * VSleave_GAULEG                              #
! #               * HOMEGROWN_ERRFUNC                           #
! #                                                             #
! ###############################################################

      MODULE vsleave_sup_aux_m


      PUBLIC
      CONTAINS


SUBROUTINE VSleave_GAULEG(X1,X2,X,W,N)

      IMPLICIT NONE
      INTEGER  , parameter:: fpk = SELECTED_REAL_KIND(15)

      INTEGER  , intent(in)  :: N
      REAL(fpk), intent(in)  :: X1, X2
      REAL(fpk), intent(out) :: X(N),W(N)

      INTEGER     :: I, M, J
      REAL(fpk)   :: XM,XL,P1,P2,P3,PP,Z,Z1
      REAL(fpk), PARAMETER :: EPS = 3.0D-14

      M=(N+1)/2
      XM=0.5D0*(X2+X1)
      XL=0.5D0*(X2-X1)

      DO I=1,M
            Z=DCOS(3.141592654D0*(I-.25D0)/(N+.5D0))
            Z1 = 0.0d0
            DO WHILE (DABS(Z-Z1).GT.EPS)
                  P1=1.D0
                  P2=0.D0
                  DO J=1,N
                        P3=P2
                        P2=P1
                        P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
                  ENDDO
                  PP=N*(Z*P1-P2)/(Z*Z-1.D0)
                  Z1=Z
                  Z=Z1-P1/PP
            ENDDO
            X(I)=XM-XL*Z
            X(N+1-I)=XM+XL*Z
            W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
            W(N+1-I)=W(I)
      ENDDO
      RETURN
      END SUBROUTINE VSleave_GAULEG

!

      SUBROUTINE HOMEGROWN_ERRFUNC ( X )

      IMPLICIT NONE

      integer, parameter :: fpk = SELECTED_REAL_KIND(15)

      REAL(fpk), intent(inout) :: x

! Returns the complementary error function erfc(x) with fractional error
! everywhere less than 1.2 * 10^7. 

      REAL(fpk) :: t,z,xin

      z = abs(x) ; xin = x
      t = 1.d0/(1.d0+0.5d0*z)
      x = t*dexp(-z*z-1.26551223d0+t*(1.00002368d0+t*(.37409196d0+ &
              t*(.09678418d0+t*(-.18628806d0+t*(.27886807d0+t* &
              (-1.13520398d0+t*(1.48851587d0+t*(-.82215223d0+t* &
              .17087277d0)))))))))
      if ( xin .lt. 0.d0 ) x = 2.d0 - x

      return
      END SUBROUTINE HOMEGROWN_ERRFUNC

!  End module

      END MODULE vsleave_sup_aux_m

