
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
! #       NEW: TOTAL COLUMN JACOBIANS          (2.4)            #
! #       NEW: BPDF Land-surface KERNELS       (2.4R)           #
! #       NEW: Thermal Emission Treatment      (2.4RT)          #
! #       Consolidated BRDF treatment          (2.4RTC)         #
! #       f77/f90 Release                      (2.5)            #
! #       External SS / New I/O Structures     (2.6)            #
! #                                                             #
! #       SURFACE-LEAVING / BRDF-SCALING       (2.7)            #
! #       TAYLOR Series / OMP THREADSAFE       (2.7)            #
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
! #  This is vbrdf_aux.f. Utility routines                      #
! #  The subroutines in vbrdf_aux are listed below with their   #
! #      source of origin (order of appearance).                #
! #                                                             #
! #     brdf_gauleg:               Numerical Recipes, 1992      #
! #     derfc_e:                   V. Natraj, 2005              #
! #     VBRDF_Fresnel_Complex      R. Spurr, 2014 (Version 2.7) #
! #     BRDF_QUADRATURE_Gaussian   R. Spurr, 2004               #
! #     BRDF_QUADRATURE_Trapezoid  R. Spurr, 2004 (not used)    #
! #                                                             #
! ###############################################################

      MODULE vbrdf_sup_aux_m

      USE VLIDORT_PARS, only : zero, one, two, half, pie
 
!  Everything public here

      public

      CONTAINS


      SUBROUTINE BRDF_GAULEG(X1,X2,X,W,N)

      implicit none

      INTEGER             N
      DOUBLE PRECISION X1,X2,X(N),W(N)
      INTEGER             I, M, J
      DOUBLE PRECISION EPS,XM,XL,P1,P2,P3,PP,Z,Z1
      PARAMETER (EPS=3.D-14)
      M=(N+1)/2
      XM=0.5D0*(X2+X1)
      XL=0.5D0*(X2-X1)
      DO I=1,M
            Z=DCOS(3.141592654D0*(I-.25D0)/(N+.5D0))
1            CONTINUE
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
            IF(DABS(Z-Z1).GT.EPS)GO TO 1
            X(I)=XM-XL*Z
            X(N+1-I)=XM+XL*Z
            W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
            W(N+1-I)=W(I)
      ENDDO
      RETURN
      END SUBROUTINE BRDF_GAULEG

!

      double precision function derfc_e(x)
      implicit none
      double precision :: x

! Returns the complementary error function erfc(x) with fractional error
! everywhere less than 1.2 * 10^7.

      double precision :: t,z

      z = dabs(x)
      t = 1.d0/(1.d0+0.5d0*z)
      derfc_e = t*dexp(-z*z-1.26551223d0+t*(1.00002368d0+t*(.37409196d0+ &
              t*(.09678418d0+t*(-.18628806d0+t*(.27886807d0+t* &
              (-1.13520398d0+t*(1.48851587d0+t*(-.82215223d0+t* &
              .17087277d0)))))))))
      if (x .lt. 0.d0) derfc_e = 2.d0-derfc_e

      return
      END function derfc_e

!

Subroutine VBRDF_Fresnel_Complex ( MR, MI, COSCHI, FP )

!  Renamed for the BRDF supplement.
!    (Same routine occurs also in the VSLEAVE suite)
!    Adapted from SixS code, this is essentially Born/Wolf computation

   implicit none
  
!  Arguments

   double precision, intent(in)  :: MR, MI, COSCHI
   double precision, intent(out) :: FP

!  Local

   double precision :: MRSQ, MISQ, MSQ, MRMI2, SINCHI_SQ, AA, A1, A2, B1, B2
   double precision :: U, V, VSQ, CMU, CPU, RR2
   double precision :: B1MU, B1PU, B2MV, B2PV, RL2

!  Calculation of FP, Complex RI

   IF ( MI.eq.zero) goto 67

   MRSQ = MR * MR ; MISQ = MI * MI
   MSQ   = MRSQ - MISQ
   MRMI2 = two * MR * MI

   SINCHI_SQ = one - COSCHI * COSCHI 
   AA = MSQ - SINCHI_SQ
   A1 = abs(AA)
   A2 = SQRT ( AA*AA + MRMI2 * MRMI2 )

   U = sqrt(half*abs(A1+A2))
   V = sqrt(half*abs(-A1+A2))
   VSQ = V * V
   CMU = ( COSCHI - U ) ; CPU = ( COSCHI + U )
   RR2 = ( CMU*CMU + VSQ ) / ( CPU*CPU + VSQ )

   B1 = MSQ   * COSCHI
   B2 = MRMI2 * COSCHI
   B1MU = B1 - U ; B1PU = B1 + U 
   B2PV = B2 + V ; B2MV = B2 - V 

   RL2 = ( B1MU*B1MU + B2PV*B2PV ) / ( B1PU*B1PU + B2MV*B2MV )
   FP = half * ( RR2 + RL2 )
   return

!  Calculation of FP. Real RI

67 continue
   MSQ = MR * MR
   SINCHI_SQ = one - COSCHI * COSCHI 
   U = sqrt(abs(MSQ - SINCHI_SQ))
   CMU = ( COSCHI - U ) ; CPU = ( COSCHI + U )
   RR2  = CMU*CMU / ( CPU*CPU )
   B1   = MSQ * COSCHI
   B1MU = B1 - U ; B1PU = B1 + U 
   RL2  = B1MU*B1MU / ( B1PU*B1PU )
   FP   = half * ( RR2 + RL2 )

!  Finish

   return
end subroutine VBRDF_Fresnel_Complex


      SUBROUTINE BRDF_QUADRATURE_GAUSSIAN &
         ( DO_BRDF_SURFEMISSION, NSTREAMS_BRDF, NBRDF_HALF, &
           X_BRDF, CX_BRDF, SX_BRDF, A_BRDF, &
           BAX_BRDF, CXE_BRDF, SXE_BRDF )

!  include file of dimensions and numbers

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Input
!  =====

!  Emission flag

      LOGICAL ::          DO_BRDF_SURFEMISSION

!  Number of streams

      INTEGER ::          NSTREAMS_BRDF, NBRDF_HALF

!  OUTPUT
!  ======

!  azimuth quadrature streams for BRDF

      DOUBLE PRECISION :: X_BRDF  ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: CX_BRDF ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: SX_BRDF ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: A_BRDF  ( MAXSTREAMS_BRDF )

!  For emission calculations

      DOUBLE PRECISION :: BAX_BRDF ( MAXSTHALF_BRDF )
      DOUBLE PRECISION :: CXE_BRDF ( MAXSTHALF_BRDF )
      DOUBLE PRECISION :: SXE_BRDF ( MAXSTHALF_BRDF )

!  local variables
!  ---------------

      INTEGER ::          I, I1, K

!  BRDF quadrature (Gauss-Legendre)
!  ---------------

!  Save these quantities for efficient coding

      CALL BRDF_GAULEG ( ZERO, ONE, X_BRDF, A_BRDF, NBRDF_HALF )
      DO I = 1, NBRDF_HALF
        I1 = I + NBRDF_HALF
        X_BRDF(I1) = - X_BRDF(I)
        A_BRDF(I1) =   A_BRDF(I)
        CXE_BRDF(I) = X_BRDF(I)
        SXE_BRDF(I) = DSQRT(ONE-X_BRDF(I)*X_BRDF(I))
      ENDDO
      DO I = 1, NSTREAMS_BRDF
        X_BRDF(I)  = PIE * X_BRDF(I)
        CX_BRDF(I) = DCOS ( X_BRDF(I) )
        SX_BRDF(I) = DSIN ( X_BRDF(I) )
      ENDDO

!  Half space cosine-weight arrays (emission only, non-Lambertian)

      IF ( DO_BRDF_SURFEMISSION ) THEN
        DO K = 1, NBRDF_HALF
          BAX_BRDF(K) = X_BRDF(K) * A_BRDF(K) / PIE
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BRDF_QUADRATURE_GAUSSIAN

! 

      SUBROUTINE BRDF_QUADRATURE_TRAPEZOID &
         ( DO_BRDF_SURFEMISSION, NSTREAMS_BRDF, NBRDF_HALF, &
           X_BRDF, CX_BRDF, SX_BRDF, A_BRDF, &
           BAX_BRDF, CXE_BRDF, SXE_BRDF )

!  include file of dimensions and numbers

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Input
!  =====

!  Emission flag

      LOGICAL ::          DO_BRDF_SURFEMISSION

!  Number of streams

      INTEGER ::          NSTREAMS_BRDF, NBRDF_HALF

!  OUTPUT
!  ======

!  azimuth quadrature streams for BRDF

      DOUBLE PRECISION :: X_BRDF  ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: CX_BRDF ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: SX_BRDF ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: A_BRDF  ( MAXSTREAMS_BRDF )

!  For emission calculations

      DOUBLE PRECISION :: BAX_BRDF ( MAXSTHALF_BRDF )
      DOUBLE PRECISION :: CXE_BRDF ( MAXSTHALF_BRDF )
      DOUBLE PRECISION :: SXE_BRDF ( MAXSTHALF_BRDF )

!  local variables
!  ---------------

      INTEGER ::          I, I1, K
      DOUBLE PRECISION :: DF1, DEL

!  BRDF quadrature (Trapezium)
!  ---------------

!  Save these quantities for efficient coding

      DF1 = DBLE(NSTREAMS_BRDF - 1 )
      DEL = TWO * PIE / DF1
        DO I = 1, NSTREAMS_BRDF
        I1 = I - 1
        X_BRDF(I) = DBLE(I1) * DEL - PIE
        X_BRDF(I) = DBLE(I1) * DEL
        CX_BRDF(I) = DCOS ( X_BRDF(I) )
        SX_BRDF(I) = DSIN ( X_BRDF(I) )
        CXE_BRDF(I) = CX_BRDF(I)
        SXE_BRDF(I) = SX_BRDF(I)
      ENDDO
        DO I = 2, NSTREAMS_BRDF - 1
        A_BRDF(I)  = DEL / PIE
      ENDDO
      A_BRDF(1)              = DEL * HALF / PIE
      A_BRDF(NSTREAMS_BRDF)  = DEL * HALF / PIE

!  Half space cosine-weight arrays (emission only, non-Lambertian)

      IF ( DO_BRDF_SURFEMISSION ) THEN
        DO K = 1, NBRDF_HALF
          BAX_BRDF(K) = X_BRDF(K) * A_BRDF(K) / PIE
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BRDF_QUADRATURE_TRAPEZOID

!  End module

      END MODULE vbrdf_sup_aux_m

