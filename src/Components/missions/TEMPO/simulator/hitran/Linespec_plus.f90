
module Linespec_plus_m

!  Contains the following routines:

!      SUBROUTINE q_lookup_lin
!      SUBROUTINE Lin_humlik

!      SUBROUTINE zhagloul_voigt    ! NOT INCLUDED

!   - All Subroutines are stand-alone, no dependencies
!   - All Subroutines are Implicit none, and public.

INTEGER, PARAMETER :: dp = KIND(1.0D0)

public
private :: dp

contains

SUBROUTINE q_lookup_lin (maxmols, maxiso, q_input, molnum, temp, if_q, q, dqdt)

! 3rd order lagrange interpolation, for evenly-spaced x-values. the known values
! of x = -1, 0, 1, 2. the unknown x normally is 0 < x < 1, although that is not
! strictly necessary (the intent is to interpolate more evenly).

IMPLICIT none
INTEGER, INTENT(IN)      :: maxmols, maxiso, molnum
REAL(KIND=dp), INTENT(IN) :: temp
LOGICAL, DIMENSION(maxmols, maxiso), INTENT(IN)               :: if_q
REAL(KIND=dp), DIMENSION(maxmols, maxiso, 148:342), INTENT(IN) :: q_input
REAL (KIND=dp), DIMENSION(maxmols, maxiso), INTENT(OUT)        :: q
REAL (KIND=dp), DIMENSION(maxmols, maxiso), INTENT(OUT)        :: dqdt

!  Linearization by R. Spurr

INTEGER       :: ia, ib, ic, id, i, j
REAL(KIND=dp)  :: tcalc, ya, yb, yc, yd, a, b, c, d,ab,ac,ad,bc,bd,cd

!LOGICAL, SAVE :: first_call = .true.
LOGICAL       :: temp_int

q    = 0.d0
dqdt = 0.d0
temp_int = .FALSE.

! temperature limits (148 - 342 degrees included here for interpolation of
! partition functions for 150 - 340 degrees).
IF (temp .LT. 150.d0 .OR. temp .GT. 340.d0) THEN
  WRITE (*, *) 'temperature out of bounds'
  STOP
END IF

!IF (first_call) THEN
!  first_call = .FALSE.
!  CALL q_load (maxmols, maxiso, q_input)
!  WRITE (*, *) 'partition function database loaded'
!END IF

! find q for the input temperature (= temp) interpolating IF necessary.
! ib = int (temp) is the location of temperature just below temp in the evenly-
! spaced array of q values. interpolate among ib - 1, ib, ib + 1, ib + 2.
ib = INT (temp)
ia = ib - 1
ic = ib + 1
id = ib + 2
! integral value of temperature? IF so, don't interpolate.
IF (temp .EQ. AINT (temp)) THEN
  temp_int = .TRUE.
ELSE
! non-integral values.
  tcalc = temp - ib
END IF

! calculate q (i, j) for existing mol, iso pairs
DO i = molnum, molnum
  DO j = 1, maxiso
    IF (if_q (i, j)) THEN
      ya = q_input (i, j, ia)
      yb = q_input (i, j, ib)
      yc = q_input (i, j, ic)
      yd = q_input (i, j, id)
      IF (temp_int) THEN
!       don't interpolate
        q (i, j)  = yb
        dqdt(i,j) = - (ya/3.d0)-(yb/2.d0) + yc - (yd/6.d0) 
      ELSE
!       interpolate
        a = tcalc + 1.d0
        b = tcalc
        c = tcalc - 1.d0
        d = tcalc - 2.d0
        ab = a*b ; ac = a*c ; ad = a*d
        bc = b*c ; bd = b*d ; cd = c*d
        q   (i,j) = - (ya * b * cd / 6.d0) + (yb * a * cd / 2.d0) - &
                      (yc * a * bd / 2.d0) + (yd * a * bc / 6.d0)
        dqdt(i,j) = - (ya*(cd+bc+bd)/6.d0) + (yb*(cd+ac+ad)/2.d0) - &
                      (yc*(ad+ab+bd)/2.d0) + (yd*(ab+bc+ac)/6.d0)
      END IF
    END IF
  END DO
END DO

RETURN
END SUBROUTINE q_lookup_lin


SUBROUTINE Lin_HUMLIK ( N, m, X, Y, dX, dY, K, dK )

!     To calculate the Faddeeva function with relative error less than 10^(-4).

implicit none

! Arguments
INTEGER N, M                                                          ! IN   Number of points
REAL(KIND=dp)    X(0:N-1)                                                  ! IN   Input x array
REAL(KIND=dp)    Y                                                         ! IN   Input y value >=0.0
REAL(KIND=dp)    K(0:N-1)                                                  ! OUT  Real (Voigt) array

REAL(KIND=dp)    dX(0:N-1,2)                                                  ! IN   Input x array
REAL(KIND=dp)    dY(2)                                                        ! IN   Input y value >=0.0
REAL(KIND=dp)    dK(0:N-1,2)                                                  ! OUT  Real (Voigt) array

! Constants
REAL        RRTPI                                                 ! 1/SQRT(pi)
PARAMETER ( RRTPI = 0.56418958 )
REAL        Y0,       Y0PY0,         Y0Q                          ! for CPF12 algorithm
PARAMETER ( Y0 = 1.5, Y0PY0 = Y0+Y0, Y0Q = Y0*Y0  )
REAL  C(0:5), S(0:5), T(0:5)
SAVE  C,      S,      T
!     SAVE preserves values of C, S and T (static) arrays between procedure calls
DATA C / 1.0117281,     -0.75197147,        0.012557727, &
     0.010022008,   -0.00024206814,     0.00000050084806 /
DATA S / 1.393237,       0.23115241,       -0.15535147,  &
     0.0062183662,   0.000091908299,   -0.00000062752596 /
DATA T / 0.31424038,     0.94778839,        1.5976826,  &
     2.2795071,      3.0206370,         3.8897249 /

! Local variables
INTEGER I, J                                                      ! Loop variables
INTEGER RG1, RG2, RG3                                             ! y polynomial flags
REAL ABX, XQ, YQ, XPQ, YRRTPI                                     ! |x|, x^2, y^2, y/SQRT(pi)
REAL XLIM0, XLIM1, XLIM2, XLIM3, XLIM4                            ! |x| on region boundaries
REAL A0, D0, D2, E0, E2, E4, H0, H2, H4, H6                       ! W4 temporary variables
REAL P0, P2, P4, P6, P8, Z0, Z2, Z4, Z6, Z8
REAL XP(0:5), XM(0:5), YP(0:5), YM(0:5)                           ! CPF12 temporary values
REAL MQ(0:5), PQ(0:5), MF(0:5), PF(0:5)
REAL D, YF, YPY0, YPY0Q
 
REAL Den, Num, MDen, Mnum, PDen, PNum, KEXP
REAL A1, A2, A3, A4, A5, A6, A7, A8, A9
REAL S0, S2, S4, S6, T0, T2, T4, T6, T8

REAL dXQ(2), dYQ(2), dXPQ(2), dYRRTPI(2)                                            ! |x|, x^2, y^2, y/SQRT(pi)
REAL dD0(2), dD2(2), dE0(2), dE2(2), dE4(2), dH0(2), dH2(2), dH4(2), dH6(2)         ! W4 temporary variables
REAL dP0(2), dP2(2), dP4(2), dP6(2), dP8(2), dZ0(2), dZ2(2), dZ4(2), dZ6(2), dZ8(2)
REAL dXP(0:5,2), dXM(0:5,2), dYP(0:5,2), dYM(0:5,2)                                 ! CPF12 temporary values
REAL dMQ(0:5,2), dPQ(0:5,2), dMF(0:5,2), dPF(0:5,2)
REAL dYF(2), dYPY0(2), dYPY0Q(2) 
 
REAL dDen(2), dNum(2), dMDen(2), dMnum(2), dPDen(2), dPNum(2)
REAL dA1(2), dA2(2), dA3(2), dA4(2), dA5(2), dA6(2), dA7(2), dA8(2), dA9(2)
REAL dS0(2), dS2(2), dS4(2), dS6(2), dT0(2), dT2(2), dT4(2), dT6(2), dT8(2)

!logical, parameter :: verbose = .true.
logical, parameter :: verbose = .false.

!***** Start of executable code *****************************************

YQ  = Y*Y                    ; dYQ(1:m)     = 2.0d0*Y*dY(1:m)   ! y^2
YRRTPI = Y*RRTPI             ; dYRRTPI(1:m) = dY(1:m)*RRTPI     ! y/SQRT(pi)

IF ( Y .GE. 70.55 ) THEN                                          ! All points
   DO I = 0, N-1                                                   ! in Region 0
      XQ   = X(I)*X(I)       ; dxQ(1:m)  = 2.0*X(I)*dX(I,1:m)
      XPQ = XQ + YQ          ; dXPQ(1:m) = dXQ(1:m) + dYQ(1:m)
      K(I) = YRRTPI / XPQ    ; dK(I,1:m) = K(I) * ( (dY(1:m)/Y)- (dXPQ(1:m)/XPQ) )
   ENDDO
   RETURN
ENDIF

! ########################################


RG1 = 1                                                           ! Set flags
RG2 = 1
RG3 = 1

if ( verbose) write(*,*)'Y < 70.55',Y,X(0)
XLIM0 = SQRT ( 15100.0 + Y*(40.0 - Y*3.6) )                       ! y<70.55
if ( verbose) write(*,*)'XLIM0 = ',XLIM0
IF ( Y .GE. 8.425 ) THEN
   XLIM1 = 0.0
ELSE
   XLIM1 = SQRT ( 164.0 - Y*(4.3 + Y*1.8) )
ENDIF
XLIM2 = 6.8 - Y
XLIM3 = 2.4*Y
XLIM4 = 18.1*Y + 1.65
IF ( Y .LE. 0.000001 ) THEN                                       ! When y<10^-6
   XLIM1 = XLIM0                                                    ! avoid W4 algorithm
   XLIM2 = XLIM0
ENDIF
if ( verbose) write(*,*)'XLIM1 = ',XLIM1
if ( verbose) write(*,*)'XLIM2 = ',XLIM2
if ( verbose) write(*,*)'XLIM3 = ',XLIM3
if ( verbose) write(*,*)'XLIM4 = ',XLIM4


!..... MAIN LOOP &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


DO I = 0, N-1                                                     ! Loop over all points

   ABX = ABS ( X(I) )                                             ! |x|
   XQ  = ABX*ABX          ; dXQ(1:m) = 2.0 * X(I) * dX(I,1:m)     ! x^2

   IF  ( ABX .GE. XLIM0 ) THEN                                   ! Region 0 algorithm
      XPQ = XQ + YQ            ; dXPQ(1:m) = dXQ(1:m) + dYQ(1:m)
      K(I) = YRRTPI / XPQ      ; dK(I,1:m) = K(I) * ( (dY(1:m)/Y)- (dXPQ(1:m)/XPQ) )

   ELSE IF ( ABX .GE. XLIM1 ) THEN                                   ! Humlicek W4 Region 1
      if ( verbose) write(*,*)'W4 Region 1'
      IF ( RG1 .NE. 0 ) THEN                                          ! First point in Region 1
         RG1 = 0
         A0 = YQ + 0.5                                                  ! Region 1 y-dependents
         D0 = A0*A0           ; dD0(1:m) = 2.0d0 * A0 * dYQ(1:m)
         D2 = YQ + YQ - 1.0   ; dD2(1:m) = 2.0d0 * dYQ(1:m)
      ENDIF
      S0 = A0 + XQ            ; dS0(1:m)  = dYQ(1:m) + dXQ(1:m)
      T2 = D2 + XQ            ; dT2(1:m)  = dD2(1:m) + dXQ(1:m)
      T0 = D0 + XQ * T2       ; dT0(1:m)  = dD0(1:m) + dXQ(1:m) * T2 + XQ * dT2(1:m)
      Den  =  RRTPI / T0      ; dDen(1:m) = - dT0(1:m) * Den / T0
      Num  =  Y * S0          ; dNum(1:m) = dY(1:m) * S0 + Y * dS0(1:m)
      K(I) = Den * Num        ; dK(I,1:m) = dDen(1:m) * Num + Den * dNum(1:m)
   
   ELSE IF ( ABX .GT. XLIM2 ) THEN                                   ! Humlicek W4 Region 2 
      if ( verbose) write(*,*)'W4 Region 2'
      IF ( RG2 .NE. 0 ) THEN                                          ! First point in Region 2
         RG2 = 0       ! Region 2 y-dependents

         a6 = 6.0  + yq        ; da6(1:m) = dyq(1:m)
         a4 = 10.5 + yq*a6     ; da4(1:m) = a6 * dyq(1:m) + yq * da6(1:m)
         a2 = 4.5  + yq*a4     ; da2(1:m) = a4 * dyq(1:m) + yq * da4(1:m)
         H0 = 0.5625 + YQ*A2   ; dH0(1:m) = a2 * dyq(1:m) + yq * da2(1:m)

         a4 = 6.0  + 4.0 * yq  ; da4(1:m) = 4.0 * dyq(1:m)
         a2 = 9.0  + yq * a4   ; da2(1:m) = a4 * dyq(1:m) + yq * da4(1:m)
         H2 = -4.5 + YQ * A2   ; dH2(1:m) = a2 * dyq(1:m) + yq * da2(1:m)

         a2 = 6.0  - 6.0 * yq  ; da2(1:m) = - 6.0 * dyq(1:m)
         H4 = 10.5 - YQ * A2   ; dH4(1:m) = - a2 * dyq(1:m) -  yq * da2(1:m)

         H6 = -6.0 + YQ* 4.0   ; dH6(1:m) = 4.0 * dyq(1:m) 


         a4 = 5.5  + yq        ; da4(1:m) = dyq(1:m)
         a2 = 8.25  + yq * a4  ; da2(1:m) = a4 * dyq(1:m) + yq * da4(1:m)
         E0 = 1.875 + YQ * A2  ; dE0(1:m) = a2 * dyq(1:m) + yq * da2(1:m)

         a2 = 1.0  + 3.0 * yq  ; da2(1:m) = 3.0 * dyq(1:m)
         E2 = 5.25 + YQ * A2   ; dE2(1:m) = a2 * dyq(1:m) +  yq * da2(1:m)

         E4 =  0.75*H6         ; dE4(1:m) = 0.75 * dH6(1:m)  

      ENDIF

      S4 = E4 + XQ            ; dS4(1:m)  = dE4(1:m) + dXQ(1:m)
      S2 = E2 + XQ * S4       ; dS2(1:m)  = dE2(1:m) + dXQ(1:m) * S4 + XQ * dS4(1:m)
      S0 = E0 + XQ * S2       ; dS0(1:m)  = dE0(1:m) + dXQ(1:m) * S2 + XQ * dS2(1:m)

      T6 = H6 + XQ            ; dT6(1:m)  = dH6(1:m) + dXQ(1:m)
      T4 = H4 + XQ * T6       ; dT4(1:m)  = dH4(1:m) + dXQ(1:m) * T6 + XQ * dT6(1:m)
      T2 = H2 + XQ * T4       ; dT2(1:m)  = dH2(1:m) + dXQ(1:m) * T4 + XQ * dT4(1:m)
      T0 = H0 + XQ * T2       ; dT0(1:m)  = dH0(1:m) + dXQ(1:m) * T2 + XQ * dT2(1:m)

      Den  =  RRTPI / T0      ; dDen(1:m) = - dT0(1:m) * Den / T0
      Num  =  Y * S0          ; dNum(1:m) = dY(1:m) * S0 + Y * dS0(1:m)
      K(I) = Den * Num        ; dK(I,1:m) = dDen(1:m) * Num + Den * dNum(1:m)

   ELSE IF ( ABX .LT. XLIM3 ) THEN                                   ! Humlicek W4 Region 3
      IF ( RG3 .NE. 0 ) THEN                                          ! First point in Region 3
      if ( verbose) write(*,*)'W4 Region 3'

         RG3 = 0

         a9 = 13.39880 + Y      ; da9(1:m) = dY(1:m)
         a8 = 88.26741 + Y*a9   ; da8(1:m) = dY(1:m) * a9 + Y*da9(1:m)
         a7 = 369.1989 + Y*a8   ; da7(1:m) = dY(1:m) * a8 + Y*da8(1:m)
         a6 = 1074.409 + Y*a7   ; da6(1:m) = dY(1:m) * a7 + Y*da7(1:m)
         a5 = 2256.981 + Y*a6   ; da5(1:m) = dY(1:m) * a6 + Y*da6(1:m)
         a4 = 3447.629 + Y*a5   ; da4(1:m) = dY(1:m) * a5 + Y*da5(1:m)
         a3 = 3764.966 + Y*a4   ; da3(1:m) = dY(1:m) * a4 + Y*da4(1:m)
         a2 = 2802.870 + Y*a3   ; da2(1:m) = dY(1:m) * a3 + Y*da3(1:m)
         a1 = 1280.829 + Y*a2   ; da1(1:m) = dY(1:m) * a2 + Y*da2(1:m)
         Z0 = 272.1014 + Y*a1   ; dZ0(1:m) = dY(1:m) * a1 + Y*da1(1:m)

         a7 = 53.59518 + Y*5.0  ; da7(1:m) = dY(1:m) * 5.0
         a6 = 266.2987 + Y*a7   ; da6(1:m) = dY(1:m) * a7 + Y*da7(1:m)
         a5 = 793.4273 + Y*a6   ; da5(1:m) = dY(1:m) * a6 + Y*da6(1:m)
         a4 = 1549.675 + Y*a5   ; da4(1:m) = dY(1:m) * a5 + Y*da5(1:m)
         a3 = 2037.310 + Y*a4   ; da3(1:m) = dY(1:m) * a4 + Y*da4(1:m)
         a2 = 1758.336 + Y*a3   ; da2(1:m) = dY(1:m) * a3 + Y*da3(1:m)
         a1 = 902.3066 + Y*a2   ; da1(1:m) = dY(1:m) * a2 + Y*da2(1:m)
         Z2 = 211.678  + Y*a1   ; dZ2(1:m) = dY(1:m) * a1 + Y*da1(1:m)

         a5 = 80.39278 + Y*10.0 ; da5(1:m) = dY(1:m) * 10.0
         a4 = 269.2916 + Y*a5   ; da4(1:m) = dY(1:m) * a5 + Y*da5(1:m)
         a3 = 479.2576 + Y*a4   ; da3(1:m) = dY(1:m) * a4 + Y*da4(1:m)
         a2 = 497.3014 + Y*a3   ; da2(1:m) = dY(1:m) * a3 + Y*da3(1:m)
         a1 = 308.1852 + Y*a2   ; da1(1:m) = dY(1:m) * a2 + Y*da2(1:m)
         Z4 = 78.86585 + Y*a1   ; dZ4(1:m) = dY(1:m) * a1 + Y*da1(1:m)

         a3 = 53.59518 + Y*10.0 ; da3(1:m) = dY(1:m) * 10.0
         a2 = 92.75679 + Y*a3   ; da2(1:m) = dY(1:m) * a3 + Y*da3(1:m)
         a1 = 55.02933 + Y*a2   ; da1(1:m) = dY(1:m) * a2 + Y*da2(1:m)
         Z6 = 22.03523 + Y*a1   ; dZ6(1:m) = dY(1:m) * a1 + Y*da1(1:m)

         a1 = 13.39880 + Y*5.0  ; da1(1:m) = dY(1:m) * 5.0
         Z8 = 1.496460 + Y*a1   ; dZ8(1:m) = dY(1:m) * a1 + Y*da1(1:m)

         a8 = 4.264678 + Y*0.3183291   ; da8(1:m) = dY(1:m)*0.3183291
         a7 = 27.93941 + Y*a8          ; da7(1:m) = dY(1:m) * a8 + Y*da8(1:m)
         a6 = 115.3772 + Y*a7          ; da6(1:m) = dY(1:m) * a7 + Y*da7(1:m)
         a5 = 328.2151 + Y*a6          ; da5(1:m) = dY(1:m) * a6 + Y*da6(1:m)
         a4 = 662.8097 + Y*a5          ; da4(1:m) = dY(1:m) * a5 + Y*da5(1:m)
         a3 = 946.8970 + Y*a4          ; da3(1:m) = dY(1:m) * a4 + Y*da4(1:m)
         a2 = 919.4955 + Y*a3          ; da2(1:m) = dY(1:m) * a3 + Y*da3(1:m)
         a1 = 549.3954 + Y*a2          ; da1(1:m) = dY(1:m) * a2 + Y*da2(1:m)
         P0 = 153.5168 + Y*a1          ; dP0(1:m) = dY(1:m) * a1 + Y*da1(1:m)

         a6 = 12.79458 + Y*1.2733163   ; da6(1:m) = dY(1:m)*1.2733163
         a5 = 56.81652 + Y*a6          ; da5(1:m) = dY(1:m) * a6 + Y*da6(1:m)
         a4 = 139.4665 + Y*a5          ; da4(1:m) = dY(1:m) * a5 + Y*da5(1:m)
         a3 = 189.7730 + Y*a4          ; da3(1:m) = dY(1:m) * a4 + Y*da4(1:m)
         a2 = 124.5975 + Y*a3          ; da2(1:m) = dY(1:m) * a3 + Y*da3(1:m)
         a1 =-1.322256 + Y*a2          ; da1(1:m) = dY(1:m) * a2 + Y*da2(1:m)
         P2 =-34.16955 + Y*a1          ; dP2(1:m) = dY(1:m) * a1 + Y*da1(1:m)

         a4 = 12.79568 + Y*1.9099744   ; da4(1:m) = dY(1:m)*1.9099744
         a3 = 29.81482 + Y*a4          ; da3(1:m) = dY(1:m) * a4 + Y*da4(1:m)
         a2 = 24.01655 + Y*a3          ; da2(1:m) = dY(1:m) * a3 + Y*da3(1:m)
         a1 = 10.46332 + Y*a2          ; da1(1:m) = dY(1:m) * a2 + Y*da2(1:m)
         P4 = 2.584042 + Y*a1          ; dP4(1:m) = dY(1:m) * a1 + Y*da1(1:m)

         a2 = 4.266322   + Y*1.273316  ; da2(1:m) = dY(1:m)*1.273316
         a1 = 0.9377051  + Y*a2        ; da1(1:m) = dY(1:m) * a2 + Y*da2(1:m)
         P6 =-0.07272979 + Y*a1        ; dP6(1:m) = dY(1:m) * a1 + Y*da1(1:m)

         P8 = 0.0005480304 + Y*0.3183291 ; dP8(1:m) = dY(1:m) * 0.3183291

      ENDIF

      S6 = P6 + XQ * P8       ; dS6(1:m)  = dP6(1:m) + dXQ(1:m) * P8 + XQ * dP8(1:m)
      S4 = P4 + XQ * S6       ; dS4(1:m)  = dP4(1:m) + dXQ(1:m) * S6 + XQ * dS6(1:m)
      S2 = P2 + XQ * S4       ; dS2(1:m)  = dP2(1:m) + dXQ(1:m) * S4 + XQ * dS4(1:m)
      S0 = P0 + XQ * S2       ; dS0(1:m)  = dP0(1:m) + dXQ(1:m) * S2 + XQ * dS2(1:m)

      T8 = Z8 + XQ            ; dT8(1:m)  = dZ8(1:m) + dXQ(1:m)
      T6 = Z6 + XQ * T8       ; dT6(1:m)  = dZ6(1:m) + dXQ(1:m) * T8 + XQ * dT8(1:m)
      T4 = Z4 + XQ * T6       ; dT4(1:m)  = dZ4(1:m) + dXQ(1:m) * T6 + XQ * dT6(1:m)
      T2 = Z2 + XQ * T4       ; dT2(1:m)  = dZ2(1:m) + dXQ(1:m) * T4 + XQ * dT4(1:m)
      T0 = Z0 + XQ * T2       ; dT0(1:m)  = dZ0(1:m) + dXQ(1:m) * T2 + XQ * dT2(1:m)

      Den  =  1.7724538 / T0  ; dDen(1:m) = - dT0(1:m) * Den / T0
      Num  =  S0              ; dNum(1:m) = dS0(1:m)
      K(I) = Den * Num        ; dK(I,1:m) = dDen(1:m) * Num + Den * dNum(1:m)

     
   ELSE                                                             ! Humlicek CPF12 algorithm

      YPY0 = Y + Y0          ; dYPY0(1:m)  = dY(1:m)
      YPY0Q = YPY0*YPY0      ; dYPY0Q(1:m) = 2.0 * YPY0 *dY(1:m)
      K(I) = 0.0             ; dK(I,1:m)   = 0.0
      DO J = 0, 5
         D = X(I) - T(J)
         MQ(J) = D*D                       ; dMQ(J,1:m) = 2.0 * D * dX(I,1:m)
         MF(J) = 1.0 / (MQ(J) + YPY0Q)     ; dMF(J,1:m) = - MF(J) * MF(J) * ( dMQ(J,1:m) + dYPY0Q(1:m) )
         XM(J) = MF(J)*D                   ; dXM(J,1:m) = dMF(J,1:m) * D    + MF(J) * dX(I,1:m)
         YM(J) = MF(J)*YPY0                ; dYM(J,1:m) = dMF(J,1:m) * YPY0 + MF(J) * dYPY0(1:m)
         D = X(I) + T(J)
         PQ(J) = D*D                       ; dPQ(J,1:m) = 2.0 * D * dX(I,1:m)
         PF(J) = 1.0 / (PQ(J) + YPY0Q)     ; dPF(J,1:m) = - PF(J) * PF(J) * ( dPQ(J,1:m) + dYPY0Q(1:m) )
         XP(J) = PF(J)*D                   ; dXP(J,1:m) = dPF(J,1:m) * D    + PF(J) * dX(I,1:m)
         YP(J) = PF(J)*YPY0                ; dYP(J,1:m) = dPF(J,1:m) * YPY0 + PF(J) * dYPY0(1:m)
      ENDDO
      
      IF ( ABX .LE. XLIM4 ) THEN                                      ! Humlicek CPF12 Region I
      if ( verbose) write(*,*)'CPF12 Region I'
         DO J = 0, 5
            K(I) = K(I) + C(J)*(YM(J)+YP(J)) - S(J)*(XM(J)-XP(J))
            dK(I,1:m) = dK(I,1:m) + C(J)*(dYM(J,1:m)+dYP(J,1:m)) - S(J)*(dXM(J,1:m)-dXP(J,1:m))
         ENDDO
      ELSE                                                            ! Humlicek CPF12 Region II
      if ( verbose) write(*,*)'CPF12 Region II'
         YF   = Y + Y0PY0                 ; dYF(1:m) = dY(1:m)
         DO J = 0, 5
            MNum = C(J)*(MQ(J)*MF(J)-Y0*YM(J)) + S(J)*YF*XM(J)  ; MDen = 1.0 / (MQ(J)+Y0Q)
            PNum = C(J)*(PQ(J)*PF(J)-Y0*YP(J)) - S(J)*YF*XP(J)  ; PDen = 1.0 / (PQ(J)+Y0Q)
            dMnum(1:m) = C(J) * ( dMQ(J,1:m)*MF(J) + MQ(J)*dMF(J,1:m) - Y0*dYM(J,1:m) ) &
                       + S(J) * ( dYF(1:m)*XM(J) + YF*dXM(J,1:m) )
            dPnum(1:m) = C(J) * ( dPQ(J,1:m)*PF(J) + PQ(J)*dPF(J,1:m) - Y0*dYP(J,1:m) ) &
                       - S(J) * ( dYF(1:m)*XP(J) + YF*dXP(J,1:m) )
            dMDen(1:m) = - dMQ(J,1:m) * Mden * MDen
            dPDen(1:m) = - dPQ(J,1:m) * Pden * PDen
            K(I) = K(I) + Mnum * MDen + Pnum * PDen
            dK(I,1:m) = dK(I,1:m) +  dMnum(1:m) * MDen + Mnum * dMDen(1:m) + dPnum(1:m) * PDen + Pnum * dPDen(1:m)
         ENDDO
         KEXP = EXP ( -XQ )
         dK(I,1:m) = Y * dK(I,1:m) + dY(1:m)*K(I) - KEXP * dXQ(1:m)
         K(I) = Y*K(I) + KEXP
      ENDIF

! End Region W4, CPF12

   ENDIF

!  End points

ENDDO
!.....
END SUBROUTINE Lin_HUMLIK

!  End module

end module Linespec_plus_m
