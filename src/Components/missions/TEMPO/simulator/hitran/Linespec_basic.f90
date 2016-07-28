
module Linespec_basic_m

!  Contains the following routines:

!      SUBROUTINE q_load
!      SUBROUTINE q_lookup
!      SUBROUTINE hitran_setup
!      SUBROUTINE voigt
!      SUBROUTINE humlik

!      SUBROUTINE zhagloul_voigt    ! NOT INCLUDED

!   - All Subroutines are stand-alone, no dependencies
!   - All Subroutines are Implicit none, and public.

INTEGER, PARAMETER :: dp = KIND(1.0D0)

public
private :: dp

contains

!

!SUBROUTINE hitran_setup (maxmols, maxiso, qpower, amu)
!
!IMPLICIT NONE
!INTEGER, INTENT(IN) :: maxmols, maxiso
!REAL(KIND=dp), DIMENSION(maxmols), INTENT(OUT)         :: qpower
!REAL(KIND=dp), DIMENSION(maxmols, maxiso), INTENT(OUT) :: amu
!
!INTEGER :: i
!
!qpower = 1.5d0
!amu = 0.d0
!
!! hitran numbers:
!!   h2o (1)
!!   co2 (2)
!!    o3 (3)
!!   n2o (4)
!!    co (5)
!!   ch4 (6)
!!    o2 (7)
!!    no (8)
!!   so2 (9)
!!   no2 (10)
!!   nh3 (11)
!!  hno3 (12)
!!    oh (13)
!!    hf (14)
!!   hcl (15)
!!   hbr (16)
!!    hi (17)
!!   clo (18)
!!   ocs (19)
!!  h2co (20)
!!  hocl (21)
!!    n2 (22)
!!   hcn (23)
!! ch3cl (24)
!!  h2o2 (25)
!!  c2h2 (26)
!!  c2h6 (27)
!!   ph3 (28)
!!  cof2 (29)
!!   sf6 (30)
!!   h2s (31)
!! hcooh (32)
!!   ho2 (33)
!!     o (34)
!!clono2 (35)
!!   no+ (36)
!!  hobr (37)
!!  c2h4 (38)
!! ch3oh (39)
!! ch3br (40)
!! ch3cn (41)
!!   cf4 (42)
!
!do i = 1, maxmols
!  if (i .eq. 34) qpower (i) = 0.d0
!  if (i .eq. 2 .or. i .eq. 4 .or. i .eq. 5 .or. i .eq. 7 .or. i .eq. 8 &
!  .or. i .eq. 13 .or. i .eq. 14 .or. i .eq. 15 .or. i .eq. 16 .or. i &
!  .eq. 17 .or. i .eq. 18 .or. i .eq. 19 .or. i .eq. 22 .or. i .eq. 23 &
!  .or. i .eq. 26 .or. i .eq. 36) qpower (i) = 1.d0
!end do
!amu (1, 1) = 18.010565  
!amu (1, 2) = 20.014811  
!amu (1, 3) = 19.014780
!amu (1, 4) = 19.016740  
!amu (1, 5) = 21.020985  
!amu (1, 6) = 20.020956  
!amu (2, 1) = 43.989830  
!amu (2, 2) = 44.993185  
!amu (2, 3) = 45.994076  
!amu (2, 4) = 44.994045  
!amu (2, 5) = 46.997431  
!amu (2, 6) = 45.997400  
!amu (2, 7) = 47.998322  
!amu (2, 8) = 46.998291  
!amu (3, 1) = 47.984745  
!amu (3, 2) = 49.988991  
!amu (3, 3) = 49.988991  
!amu (3, 4) = 48.988960  
!amu (3, 5) = 48.988960  
!amu (4, 1) = 44.001062  
!amu (4, 2) = 44.998096  
!amu (4, 3) = 44.998096  
!amu (4, 4) = 46.005308  
!amu (4, 5) = 45.005278  
!amu (5, 1) = 27.994915  
!amu (5, 2) = 28.998270  
!amu (5, 3) = 29.999161  
!amu (5, 4) = 28.999130  
!amu (5, 5) = 31.002516  
!amu (5, 6) = 30.002485  
!amu (6, 1) = 16.031300  
!amu (6, 2) = 17.034655  
!amu (6, 3) = 17.037475  
!amu (7, 1) = 31.989830  
!amu (7, 2) = 33.994076  
!amu (7, 3) = 32.994045  
!amu (8, 1) = 29.997989  
!amu (8, 2) = 30.995023  
!amu (8, 3) = 32.002234  
!amu (9, 1) = 63.961901  
!amu (9, 2) = 65.957695  
!amu (10, 1) = 45.992904  
!amu (11, 1) = 17.026549  
!amu (11, 2) = 18.023583  
!amu (12, 1) = 62.995644  
!amu (13, 1) = 17.002740  
!amu (13, 2) = 19.006986  
!amu (13, 3) = 18.008915  
!amu (14, 1) = 20.006229  
!amu (15, 1) = 35.976678  
!amu (15, 2) = 37.973729  
!amu (16, 1) = 79.926160  
!amu (16, 2) = 81.924115  
!amu (17, 1) = 127.912297  
!amu (18, 1) = 50.963768  
!amu (18, 2) = 52.960819  
!amu (19, 1) = 59.966986  
!amu (19, 2) = 61.962780  
!amu (19, 3) = 60.970341  
!amu (19, 4) = 60.966371  
!amu (19, 5) = 61.971231  
!amu (20, 1) = 30.010565  
!amu (20, 2) = 31.013920  
!amu (20, 3) = 32.014811  
!amu (21, 1) = 51.971593  
!amu (21, 2) = 53.968644  
!amu (22, 1) = 28.006147  
!amu (23, 1) = 27.010899  
!amu (23, 2) = 28.014254  
!amu (23, 3) = 28.007933  
!amu (24, 1) = 49.992328  
!amu (24, 2) = 51.989379  
!amu (25, 1) = 34.005480  
!amu (26, 1) = 26.015650  
!amu (26, 2) = 27.019005  
!amu (27, 1) = 30.046950  
!amu (28, 1) = 33.997238  
!amu (29, 1) = 65.991722  
!amu (30, 1) = 145.962492
!amu (31, 1) = 33.987721
!amu (31, 2) = 35.983515
!amu (31, 3) = 34.987105
!amu (32, 1) = 46.005480
!amu (33, 1) = 32.997655
!amu (34, 1) = 15.994915
!amu (35, 1) = 96.956672
!amu (35, 2) = 98.953723
!amu (26, 1) = 29.997989
!amu (37, 1) = 95.921076
!amu (37, 2) = 97.919027
!amu (38, 1) = 28.031300
!amu (38, 2) = 29.034655
!amu (39, 1) = 32.026215
!
!RETURN
!END SUBROUTINE hitran_setup

!SUBROUTINE hitran_setup (maxmols, maxiso, qpower, amu)
!
!IMPLICIT NONE
!INTEGER, INTENT(IN) :: maxmols, maxiso
!REAL(KIND=dp), DIMENSION(maxmols), INTENT(OUT)         :: qpower
!REAL(KIND=dp), DIMENSION(maxmols, maxiso), INTENT(OUT) :: amu
!
!INTEGER :: i
!
!qpower = 1.5d0
!amu = 0.d0
!

SUBROUTINE hitran_setup (maxmols, maxiso, amu)

IMPLICIT NONE
INTEGER, INTENT(IN)                                   :: maxmols, maxiso
REAL(KIND=dp), DIMENSION(maxmols, maxiso), INTENT(OUT) :: amu

! hitran numbers: h2o (1), co2 (2), o3 (3), n2o (4), co (5), ch4 (6), o2 (7),
! no (8), so2 (9), no2 (10), nh3 (11), hno3 (12), oh (13), hf (14), hcl (15),
! hbr (16), hi (17), clo (18), ocs (19), h2co (20), hocl (21), n2 (22),
! hcn (23), ch3cl (24), h2o2 (25), c2h2 (26), c2h6 (27), ph3 (28), cof2 (29),
! sf6 (30), h2s (31), hcooh (32), ho2 (33), o (34), clono2 (35), no+ (36),
! hobr (37), c2h4 (38), ch3oh (39), ch3br (40; not included in parsum.dat),
! ch3cn (41), cf4 (42)

amu = 0.d0
amu (1, 1) = 18.010565
amu (1, 2) = 20.014811
amu (1, 3) = 19.014780
amu (1, 4) = 19.016740
amu (1, 5) = 21.020985
amu (1, 6) = 20.020956
amu (2, 1) = 43.989830
amu (2, 2) = 44.993185
amu (2, 3) = 45.994076
amu (2, 4) = 44.994045
amu (2, 5) = 46.997431
amu (2, 6) = 45.997400
amu (2, 7) = 47.998322
amu (2, 8) = 46.998291
amu (3, 1) = 47.984745
amu (3, 2) = 49.988991
amu (3, 3) = 49.988991
amu (3, 4) = 48.988960
amu (3, 5) = 48.988960
amu (4, 1) = 44.001062
amu (4, 2) = 44.998096
amu (4, 3) = 44.998096
amu (4, 4) = 46.005308
amu (4, 5) = 45.005278
amu (5, 1) = 27.994915
amu (5, 2) = 28.998270
amu (5, 3) = 29.999161
amu (5, 4) = 28.999130
amu (5, 5) = 31.002516
amu (5, 6) = 30.002485
amu (6, 1) = 16.031300
amu (6, 2) = 17.034655
amu (6, 3) = 17.037475
amu (7, 1) = 31.989830
amu (7, 2) = 33.994076
amu (7, 3) = 32.994045
amu (8, 1) = 29.997989
amu (8, 2) = 30.995023
amu (8, 3) = 32.002234
amu (9, 1) = 63.961901
amu (9, 2) = 65.957695
amu (10, 1) = 45.992904
amu (11, 1) = 17.026549
amu (11, 2) = 18.023583
amu (12, 1) = 62.995644
amu (13, 1) = 17.002740
amu (13, 2) = 19.006986
amu (13, 3) = 18.008915
amu (14, 1) = 20.006229
amu (15, 1) = 35.976678
amu (15, 2) = 37.973729
amu (16, 1) = 79.926160
amu (16, 2) = 81.924115
amu (17, 1) = 127.912297
amu (18, 1) = 50.963768
amu (18, 2) = 52.960819
amu (19, 1) = 59.966986
amu (19, 2) = 61.962780
amu (19, 3) = 60.970341
amu (19, 4) = 60.966371
amu (19, 5) = 61.971231
amu (20, 1) = 30.010565
amu (20, 2) = 31.013920
amu (20, 3) = 32.014811
amu (21, 1) = 51.971593
amu (21, 2) = 53.968644
amu (22, 1) = 28.006147
amu (23, 1) = 27.010899
amu (23, 2) = 28.014254
amu (23, 3) = 28.007933
amu (24, 1) = 49.992328
amu (24, 2) = 51.989379
amu (25, 1) = 34.005480
amu (26, 1) = 26.015650
amu (26, 2) = 27.019005
amu (27, 1) = 30.046950
amu (28, 1) = 33.997238
amu (29, 1) = 65.991722
amu (30, 1) = 145.962492
amu (31, 1) = 33.987721
amu (31, 2) = 35.983515
amu (31, 3) = 34.987105
amu (32, 1) = 46.005480
amu (33, 1) = 32.997655
amu (34, 1) = 15.994915
amu (35, 1) = 96.956672
amu (35, 2) = 98.953723
amu (36, 1) = 29.997989
amu (37, 1) = 95.921076
amu (37, 2) = 97.919027
amu (38, 1) = 28.031300
amu (38, 2) = 29.034655
amu (39, 1) = 32.026215
amu (41, 1) = 41.026549
amu (42, 1) = 87.993616

RETURN
END SUBROUTINE hitran_setup

SUBROUTINE q_lookup (maxmols, maxiso, q_input, molnum, temp, if_q, q)

! 3rd order lagrange interpolation, for evenly-spaced x-values. the known values
! of x = -1, 0, 1, 2. the unknown x normally is 0 < x < 1, although that is not
! strictly necessary (the intent is to interpolate more evenly).

IMPLICIT none
INTEGER, INTENT(IN)      :: maxmols, maxiso, molnum
REAL(KIND=dp), INTENT(IN) :: temp
LOGICAL, DIMENSION(maxmols, maxiso), INTENT(IN)               :: if_q
REAL(KIND=dp), DIMENSION(maxmols, maxiso, 148:342), INTENT(IN) :: q_input
REAL (KIND=dp), DIMENSION(maxmols, maxiso), INTENT(OUT)        :: q

INTEGER       :: ia, ib, ic, id, i, j
REAL(KIND=dp)  :: tcalc, ya, yb, yc, yd, a, b, c, d

!LOGICAL, SAVE :: first_call = .true.
LOGICAL       :: temp_int

q = 0.d0
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
! integral value of temperature? IF so, don't interpolate.
IF (temp .EQ. AINT (temp)) THEN
  temp_int = .TRUE.
ELSE
! non-integral values.
  ia = ib - 1
  ic = ib + 1
  id = ib + 2
  tcalc = temp - ib
END IF

! calculate q (i, j) for existing mol, iso pairs
DO i = molnum, molnum
  DO j = 1, maxiso
    IF (if_q (i, j)) THEN
      IF (temp_int) THEN
!       don't interpolate
        q (i, j) = q_input (i, j, ib)
      ELSE
!       interpolate
        ya = q_input (i, j, ia)
        yb = q_input (i, j, ib)
        yc = q_input (i, j, ic)
        yd = q_input (i, j, id)
        a = tcalc + 1.d0
        b = tcalc
        c = tcalc - 1.d0
        d = tcalc - 2.d0
        q (i, j) = - (ya * b * c * d / 6.d0) + (yb * a * c * d / 2.d0) - &
        (yc * a * b * d / 2.d0) + (yd * a * b * c / 6.d0)
      END IF
    END IF
  END DO
END DO

RETURN
END SUBROUTINE q_lookup


SUBROUTINE q_load (maxmols, maxiso, q_input, q_file)

IMPLICIT NONE
INTEGER, INTENT(IN) :: maxmols, maxiso
REAL (KIND=dp), DIMENSION(maxmols, maxiso, 148:342), INTENT(out) :: q_input
CHARACTER (len=*), INTENT(IN)   :: q_file 

INTEGER             :: i, j, k, l
!CHARACTER (len=132) :: q_file 

q_input = 0.d0
!q_file  = '../geocape_data/HITRAN/hitran08-parsum.resorted'
!OPEN (unit = 22, file = q_file, status = 'old')

OPEN (unit = 22, file = adjustl(trim(q_file)), status = 'old')

DO i = 1, 95
   READ (22, *) j, k
   DO l = 148, 342
      READ (22, *) q_input (j, k, l)
   ENDDO
ENDDO

CLOSE (unit = 22)
RETURN
END SUBROUTINE q_load


SUBROUTINE voigt (x, a, v, ndim, nvlo, nvhi)

! the following calculated voigt values at all grid values for each point
! subroutine voigt (x, a, v, nx)
! voigt first and second derivatives commented out
! subroutine voigt (x, a, v, dv, d2v, nx)

IMPLICIT NONE
INTEGER, INTENT(IN)      :: ndim, nvlo, nvhi
REAL(KIND=dp), INTENT(IN) :: a
REAL(KIND=dp), DIMENSION(ndim), INTENT(IN)  :: x
REAL(KIND=dp), DIMENSION(ndim), INTENT(OUT) :: v ! , dv, d2v 

REAL (KIND=dp), PARAMETER  :: sqrtpi =1.77245385090551d0, &
     twooverpi = 0.63661977236758d0,  fouroverpi = 1.27323954473516d0

INTEGER      :: capn, i, nu, n, in, np1
REAL(KIND=dp) :: lamda, sfac, absx, s, h, h2, r1, r2, s1, s2, t1, t2, c !, c2v
LOGICAL      :: b

! a = 0.
IF (a < 1.0d-8) THEN
   v(nvlo:nvhi) = dexp(-x(nvlo:nvhi)**2) / sqrtpi
   !dv(nvlo:nvhi) = -2.0d0 * x(nvlo:nvhi) * v(nvlo:nvhi)
   !d2v(nvlo:nvhi) = (4.0d0 * x(nvlo:nvhi) ** 2 - 2.d0) * v(nvlo:nvhi)

   ! add lorentzian check here, for speed
ELSE
  ! coefficient for second derivative
  ! c2v = 4.0d0 * a * a + 2.d0

  sfac = 1.0d0 - a / 4.29d0
  DO i = nvlo, nvhi
  ! do i = 1, nx
     absx = dabs (x (i))
     IF ((a < 4.29d0) .AND. (absx < 5.33d0)) THEN
        s = sfac * dsqrt (1.d0 - (x (i) / 5.33d0)**2)
        h = 1.6d0 * s
        h2 = 2.0d0 * h
        capn = 6.d0 + 23.0d0 * s
        lamda = h2**capn
        nu = 9.0d0 + 21.0d0 * s
     ELSE
        h = 0.0d0
        capn = 0
        nu = 8
     ENDIF
     b = (h == 0.0d0) .or. (lamda == 0.0d0)
     r1 = 0.0d0
     r2 = 0.0d0
     s1 = 0.0d0
     s2 = 0.0d0
     n = nu
     DO in = 1, nu + 1
        np1 = n + 1
!        t1 = a + h + dfloat (np1) * r1
!        t2 = absx - dfloat (np1) * r2
        t1 = a + h + real(np1,dp) * r1
        t2 = absx - real(np1,dp) * r2
        c = .5d0 / (t1 * t1 + t2 * t2)
        r1 = c * t1
        r2 = c * t2
        IF ((h > 0.0d0) .AND. (n <= capn)) THEN
           t1 = lamda + s1
           s1 = r1 * t1 - r2 * s2
           s2 = r2 * t1 + r1 * s2
           lamda = lamda / h2
        ENDIF
        n = n - 1
     ENDDO
     IF (b) THEN
        v (i) = twooverpi * r1
        !     dv (i) = fouroverpi * (a * r2 - absx * r1)
     ELSE
        v (i) = twooverpi * s1
        !     dv (i) = fouroverpi * (a * s2 - absx * s1)
     ENDIF
     !   dv (i) = -dsign (dv (i), x (i))
     !   d2v (i) = fouroverpi * a - (c2v + 4.d0 * x (i) * x (i)) * &
     !   v (i) - 4.d0 * x (i) * dv (i)
  ENDDO
ENDIF

RETURN
END SUBROUTINE voigt

SUBROUTINE HUMLIK ( N, X, Y, K )

!     To calculate the Faddeeva function with relative error less than 10^(-4).

implicit none

! Arguments
INTEGER N                                                         ! IN   Number of points
REAL(KIND=dp) X(0:N-1)                                            ! IN   Input x array
REAL(KIND=dp) Y                                                   ! IN   Input y value >=0.0
REAL(KIND=dp) K(0:N-1)                                            ! OUT  Real (Voigt) array

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
REAL ABX, XQ, YQ, YRRTPI                                          ! |x|, x^2, y^2, y/SQRT(pi)
REAL XLIM0, XLIM1, XLIM2, XLIM3, XLIM4                            ! |x| on region boundaries
REAL A0, D0, D2, E0, E2, E4, H0, H2, H4, H6                       ! W4 temporary variables
REAL P0, P2, P4, P6, P8, Z0, Z2, Z4, Z6, Z8
REAL XP(0:5), XM(0:5), YP(0:5), YM(0:5)                           ! CPF12 temporary values
REAL MQ(0:5), PQ(0:5), MF(0:5), PF(0:5)
REAL D, YF, YPY0, YPY0Q  
!logical, parameter :: verbose = .true.
logical, parameter :: verbose = .false.

!***** Start of executable code *****************************************

YQ  = Y*Y                                                         ! y^2
YRRTPI = Y*RRTPI                                                  ! y/SQRT(pi)

IF ( Y .GE. 70.55 ) THEN                                          ! All points
   if ( verbose) write(*,*)'Y.ge.70.55',Y
   DO I = 0, N-1                                                   ! in Region 0
      XQ   = X(I)*X(I)
      K(I) = YRRTPI / (XQ + YQ)
   ENDDO
   RETURN
ENDIF

RG1 = 1                                                           ! Set flags
RG2 = 1
RG3 = 1

if ( verbose) write(*,*)'Y < 70.55',Y
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

!.....
DO I = 0, N-1                                                     ! Loop over all points
   ABX = ABS ( X(I) )                                               ! |x|
   XQ  = ABX*ABX                                                    ! x^2
   IF     ( ABX .GE. XLIM0 ) THEN                                   ! Region 0 algorithm
      K(I) = YRRTPI / (XQ + YQ)
       if ( verbose) write(*,*)'W4 Region 0'
     
   ELSEIF ( ABX .GE. XLIM1 ) THEN                                   ! Humlicek W4 Region 1
      IF ( RG1 .NE. 0 ) THEN                                          ! First point in Region 1
         RG1 = 0
         A0 = YQ + 0.5                                                  ! Region 1 y-dependents
         D0 = A0*A0
         D2 = YQ + YQ - 1.0
      if ( verbose) write(*,*)'W4 Region 1'
      ENDIF
      D = RRTPI / (D0 + XQ*(D2 + XQ))
      K(I) = D*Y   *(A0 + XQ)
      
   ELSEIF ( ABX .GT. XLIM2 ) THEN                                   ! Humlicek W4 Region 2 
      IF ( RG2 .NE. 0 ) THEN                                          ! First point in Region 2
         RG2 = 0
         H0 =  0.5625 + YQ*(4.5 + YQ*(10.5 + YQ*(6.0 + YQ)))            ! Region 2 y-dependents
         H2 = -4.5    + YQ*(9.0 + YQ*( 6.0 + YQ* 4.0))
         H4 = 10.5    - YQ*(6.0 - YQ*  6.0)
         H6 = -6.0    + YQ* 4.0
         E0 =  1.875  + YQ*(8.25 + YQ*(5.5 + YQ))
         E2 =  5.25   + YQ*(1.0  + YQ* 3.0)
         E4 =  0.75*H6
      ENDIF
      if ( verbose) write(*,*)'W4 Region 2'
      D = RRTPI / (H0 + XQ*(H2 + XQ*(H4 + XQ*(H6 + XQ))))
      K(I) = D*Y   *(E0 + XQ*(E2 + XQ*(E4 + XQ)))
      
   ELSEIF ( ABX .LT. XLIM3 ) THEN                                   ! Humlicek W4 Region 3
      if ( verbose) write(*,*)'W4 Region 3'
      IF ( RG3 .NE. 0 ) THEN                                          ! First point in Region 3
         RG3 = 0
         Z0 = 272.1014     + Y*(1280.829 + Y*(2802.870 + Y*(3764.966  &    ! Region 3 y-dependents
              + Y*(3447.629 + Y*(2256.981 + Y*(1074.409 + Y*(369.1989 &
              + Y*(88.26741 + Y*(13.39880 + Y)))))))))
         Z2 = 211.678      + Y*(902.3066 + Y*(1758.336 + Y*(2037.310  &
              + Y*(1549.675 + Y*(793.4273 + Y*(266.2987 &
              + Y*(53.59518 + Y*5.0)))))))
         Z4 = 78.86585     + Y*(308.1852 + Y*(497.3014 + Y*(479.2576  &
              + Y*(269.2916 + Y*(80.39278 + Y*10.0)))))
         Z6 = 22.03523     + Y*(55.02933 + Y*(92.75679 + Y*(53.59518  &
              + Y*10.0)))
         Z8 = 1.496460     + Y*(13.39880 + Y*5.0)
         P0 = 153.5168     + Y*(549.3954 + Y*(919.4955 + Y*(946.8970  &
              + Y*(662.8097 + Y*(328.2151 + Y*(115.3772 + Y*(27.93941 &
              + Y*(4.264678 + Y*0.3183291))))))))
         P2 = -34.16955    + Y*(-1.322256+ Y*(124.5975 + Y*(189.7730  &
              + Y*(139.4665 + Y*(56.81652 + Y*(12.79458               &
              + Y*1.2733163))))))
         P4 = 2.584042     + Y*(10.46332 + Y*(24.01655 + Y*(29.81482  &
              + Y*(12.79568 + Y*1.9099744))))
         P6 = -0.07272979  + Y*(0.9377051+ Y*(4.266322 + Y*1.273316))
         P8 = 0.0005480304 + Y*0.3183291
      ENDIF
      D = 1.7724538 / (Z0 + XQ*(Z2 + XQ*(Z4 + XQ*(Z6 + XQ*(Z8+XQ)))))
      K(I) = D*(P0 + XQ*(P2 + XQ*(P4 + XQ*(P6 + XQ*P8))))
      
   ELSE                                                             ! Humlicek CPF12 algorithm
      YPY0 = Y + Y0
      YPY0Q = YPY0*YPY0
      K(I) = 0.0
      DO J = 0, 5
         D = X(I) - T(J)
         MQ(J) = D*D
         MF(J) = 1.0 / (MQ(J) + YPY0Q)
         XM(J) = MF(J)*D
         YM(J) = MF(J)*YPY0
         D = X(I) + T(J)
         PQ(J) = D*D
         PF(J) = 1.0 / (PQ(J) + YPY0Q)
         XP(J) = PF(J)*D
         YP(J) = PF(J)*YPY0
      ENDDO
      
      IF ( ABX .LE. XLIM4 ) THEN                                      ! Humlicek CPF12 Region I
         if ( verbose) write(*,*)'CPF12 Region I'
         DO J = 0, 5
            K(I) = K(I) + C(J)*(YM(J)+YP(J)) - S(J)*(XM(J)-XP(J))
         ENDDO
 
      ELSE                                                            ! Humlicek CPF12 Region II
         if ( verbose) write(*,*)'CPF12 Region II'
         YF   = Y + Y0PY0
         DO J = 0, 5
            K(I) = K(I)   &
                 + (C(J)*(MQ(J)*MF(J)-Y0*YM(J)) + S(J)*YF*XM(J)) / (MQ(J)+Y0Q) &
                 + (C(J)*(PQ(J)*PF(J)-Y0*YP(J)) - S(J)*YF*XP(J)) / (PQ(J)+Y0Q)
         ENDDO
         K(I) = Y*K(I) + EXP ( -XQ )
      ENDIF
   ENDIF
ENDDO
!.....
END SUBROUTINE HUMLIK

!  End module

end module Linespec_basic_m
