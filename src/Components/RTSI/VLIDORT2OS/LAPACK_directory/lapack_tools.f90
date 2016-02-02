
      MODULE lapack_tools

      CONTAINS

      SUBROUTINE DGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            INFO, KL, KU, LDAB, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   AB( LDAB, * )
!     ..
!
!  Purpose
!  =======
!
!  DGBTF2 computes an LU factorization of a real m-by-n band matrix A
!  using partial pivoting with row interchanges.
!
!  This is the unblocked version of the algorithm, calling Level 2 BLAS.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  KL      (input) INTEGER
!          The number of subdiagonals within the band of A.  KL >= 0.
!
!  KU      (input) INTEGER
!          The number of superdiagonals within the band of A.  KU >= 0.
!
!  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
!          On entry, the matrix A in band storage, in rows KL+1 to
!          2*KL+KU+1; rows 1 to KL of the array need not be set.
!          The j-th column of A is stored in the j-th column of the
!          array AB as follows:
!          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
!
!          On exit, details of the factorization: U is stored as an
!          upper triangular band matrix with KL+KU superdiagonals in
!          rows 1 to KL+KU+1, and the multipliers used during the
!          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
!          See below for further details.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  Further Details
!  ===============
!
!  The band storage scheme is illustrated by the following example, when
!  M = N = 6, KL = 2, KU = 1:
!
!  On entry:                       On exit:
!
!      *    *    *    +    +    +       *    *    *   u14  u25  u36
!      *    *    +    +    +    +       *    *   u13  u24  u35  u46
!      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
!     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
!
!  Array elements marked * are not used by the routine; elements marked
!  + need not be set on entry, but are required by the routine to store
!  elements of U, because of fill-in resulting from the row
!  interchanges.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J, JP, JU, KM, KV
!     ..
!     .. External Functions ..
!      INTEGER            IDAMAX
!      EXTERNAL           IDAMAX
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DGER, DSCAL, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     KV is the number of superdiagonals in the factor U, allowing for
!     fill-in.
!
      KV = KU + KL
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KL+KV+1 ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGBTF2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 ) &
         RETURN
!
!     Gaussian elimination with partial pivoting
!
!     Set fill-in elements in columns KU+2 to KV to zero.
!
      DO 20 J = KU + 2, MIN( KV, N )
         DO 10 I = KV - J + 2, KL
            AB( I, J ) = ZERO
   10    CONTINUE
   20 CONTINUE
!
!     JU is the index of the last column affected by the current stage
!     of the factorization.
!
      JU = 1
!
      DO 40 J = 1, MIN( M, N )
!
!        Set fill-in elements in column J+KV to zero.
!
         IF( J+KV.LE.N ) THEN
            DO 30 I = 1, KL
               AB( I, J+KV ) = ZERO
   30       CONTINUE
         END IF
!
!        Find pivot and test for singularity. KM is the number of
!        subdiagonal elements in the current column.
!
         KM = MIN( KL, M-J )
         JP = IDAMAX( KM+1, AB( KV+1, J ), 1 )
         IPIV( J ) = JP + J - 1
         IF( AB( KV+JP, J ).NE.ZERO ) THEN
            JU = MAX( JU, MIN( J+KU+JP-1, N ) )
!
!           Apply interchange to columns J to JU.
!
            IF( JP.NE.1 ) &
               CALL DSWAP( JU-J+1, AB( KV+JP, J ), LDAB-1, &
                           AB( KV+1, J ), LDAB-1 )
!
            IF( KM.GT.0 ) THEN
!
!              Compute multipliers.
!
               CALL DSCAL( KM, ONE / AB( KV+1, J ), AB( KV+2, J ), 1 )
!
!              Update trailing submatrix within the band.
!
               IF( JU.GT.J ) &
                  CALL DGER( KM, JU-J, -ONE, AB( KV+2, J ), 1, &
                             AB( KV, J+1 ), LDAB-1, AB( KV+1, J+1 ), &
                             LDAB-1 )
            END IF
         ELSE
!
!           If pivot is zero, set INFO to the index of the pivot
!           unless a zero pivot has already been found.
!
            IF( INFO.EQ.0 ) &
               INFO = J
         END IF
   40 CONTINUE
      RETURN
!
!     End of DGBTF2
!
      END SUBROUTINE DGBTF2


      SUBROUTINE DLABAD( SMALL, LARGE )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   LARGE, SMALL
!     ..
!
!  Purpose
!  =======
!
!  DLABAD takes as input the values computed by DLAMCH for underflow and
!  overflow, and returns the square root of each of these values if the
!  log of LARGE is sufficiently large.  This subroutine is intended to
!  identify machines with a large exponent range, such as the Crays, and
!  redefine the underflow and overflow limits to be the square roots of
!  the values computed by DLAMCH.  This subroutine is needed because
!  DLAMCH does not compensate for poor arithmetic in the upper half of
!  the exponent range, as is found on a Cray.
!
!  Arguments
!  =========
!
!  SMALL   (input/output) DOUBLE PRECISION
!          On entry, the underflow threshold as computed by DLAMCH.
!          On exit, if LOG10(LARGE) is sufficiently large, the square
!          root of SMALL, otherwise unchanged.
!
!  LARGE   (input/output) DOUBLE PRECISION
!          On entry, the overflow threshold as computed by DLAMCH.
!          On exit, if LOG10(LARGE) is sufficiently large, the square
!          root of LARGE, otherwise unchanged.
!
!  =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC          LOG10, SQRT
!     ..
!     .. Executable Statements ..
!
!     If it looks like we're on a Cray, take the square root of
!     SMALL and LARGE to avoid overflow and underflow problems.
!
      IF( LOG10( LARGE ).GT.2000.D0 ) THEN
         SMALL = SQRT( SMALL )
         LARGE = SQRT( LARGE )
      END IF
!
      RETURN
!
!     End of DLABAD
!
      END SUBROUTINE DLABAD


      SUBROUTINE DLACON( N, V, X, ISGN, EST, KASE )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            KASE, N
      DOUBLE PRECISION   EST
!     ..
!     .. Array Arguments ..
      INTEGER            ISGN( * )
      DOUBLE PRECISION   V( * ), X( * )
!     ..
!
!  Purpose
!  =======
!
!  DLACON estimates the 1-norm of a square, real matrix A.
!  Reverse communication is used for evaluating matrix-vector products.
!
!  Arguments
!  =========
!
!  N      (input) INTEGER
!         The order of the matrix.  N >= 1.
!
!  V      (workspace) DOUBLE PRECISION array, dimension (N)
!         On the final return, V = A*W,  where  EST = norm(V)/norm(W)
!         (W is not returned).
!
!  X      (input/output) DOUBLE PRECISION array, dimension (N)
!         On an intermediate return, X should be overwritten by
!               A * X,   if KASE=1,
!               A' * X,  if KASE=2,
!         and DLACON must be re-called with all the other parameters
!         unchanged.
!
!  ISGN   (workspace) INTEGER array, dimension (N)
!
!  EST    (input/output) DOUBLE PRECISION
!         On entry with KASE = 1 or 2 and JUMP = 3, EST should be
!         unchanged from the previous call to DLACON.
!         On exit, EST is an estimate (a lower bound) for norm(A).
!
!  KASE   (input/output) INTEGER
!         On the initial call to DLACON, KASE should be 0.
!         On an intermediate return, KASE will be 1 or 2, indicating
!         whether X should be overwritten by A * X  or A' * X.
!         On the final return from DLACON, KASE will again be 0.
!
!  Further Details
!  ======= =======
!
!  Contributed by Nick Higham, University of Manchester.
!  Originally named SONEST, dated March 16, 1988.
!
!  Reference: N.J. Higham, "FORTRAN codes for estimating the one-norm of
!  a real or complex matrix, with applications to condition estimation",
!  ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            ITMAX
      PARAMETER          ( ITMAX = 5 )
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, ITER, J, JLAST, JUMP
      DOUBLE PRECISION   ALTSGN, ESTOLD, TEMP
!     ..
!     .. External Functions ..
!      INTEGER            IDAMAX
!      DOUBLE PRECISION   DASUM
!      EXTERNAL           IDAMAX, DASUM
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DCOPY
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, NINT, SIGN
!     ..
!     .. Save statement ..
      SAVE
!     ..
!     .. Executable Statements ..
!
      IF( KASE.EQ.0 ) THEN
         DO 10 I = 1, N
            X( I ) = ONE / DBLE( N )
   10    CONTINUE
         KASE = 1
         JUMP = 1
         RETURN
      END IF
!
      GO TO ( 20, 40, 70, 110, 140 )JUMP
!
!     ................ ENTRY   (JUMP = 1)
!     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
!
   20 CONTINUE
      IF( N.EQ.1 ) THEN
         V( 1 ) = X( 1 )
         EST = ABS( V( 1 ) )
!        ... QUIT
         GO TO 150
      END IF
      EST = DASUM( N, X, 1 )
!
      DO 30 I = 1, N
         X( I ) = SIGN( ONE, X( I ) )
         ISGN( I ) = NINT( X( I ) )
   30 CONTINUE
      KASE = 2
      JUMP = 2
      RETURN
!
!     ................ ENTRY   (JUMP = 2)
!     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
!
   40 CONTINUE
      J = IDAMAX( N, X, 1 )
      ITER = 2
!
!     MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
!
   50 CONTINUE
      DO 60 I = 1, N
         X( I ) = ZERO
   60 CONTINUE
      X( J ) = ONE
      KASE = 1
      JUMP = 3
      RETURN
!
!     ................ ENTRY   (JUMP = 3)
!     X HAS BEEN OVERWRITTEN BY A*X.
!
   70 CONTINUE
      CALL DCOPY( N, X, 1, V, 1 )
      ESTOLD = EST
      EST = DASUM( N, V, 1 )
      DO 80 I = 1, N
         IF( NINT( SIGN( ONE, X( I ) ) ).NE.ISGN( I ) ) &
            GO TO 90
   80 CONTINUE
!     REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
      GO TO 120
!
   90 CONTINUE
!     TEST FOR CYCLING.
      IF( EST.LE.ESTOLD ) &
         GO TO 120
!
      DO 100 I = 1, N
         X( I ) = SIGN( ONE, X( I ) )
         ISGN( I ) = NINT( X( I ) )
  100 CONTINUE
      KASE = 2
      JUMP = 4
      RETURN
!
!     ................ ENTRY   (JUMP = 4)
!     X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
!
  110 CONTINUE
      JLAST = J
      J = IDAMAX( N, X, 1 )
      IF( ( X( JLAST ).NE.ABS( X( J ) ) ) .AND. ( ITER.LT.ITMAX ) ) THEN
         ITER = ITER + 1
         GO TO 50
      END IF
!
!     ITERATION COMPLETE.  FINAL STAGE.
!
  120 CONTINUE
      ALTSGN = ONE
      DO 130 I = 1, N
         X( I ) = ALTSGN*( ONE+DBLE( I-1 ) / DBLE( N-1 ) )
         ALTSGN = -ALTSGN
  130 CONTINUE
      KASE = 1
      JUMP = 5
      RETURN
!
!     ................ ENTRY   (JUMP = 5)
!     X HAS BEEN OVERWRITTEN BY A*X.
!
  140 CONTINUE
      TEMP = TWO*( DASUM( N, X, 1 ) / DBLE( 3*N ) )
      IF( TEMP.GT.EST ) THEN
         CALL DCOPY( N, X, 1, V, 1 )
         EST = TEMP
      END IF
!
  150 CONTINUE
      KASE = 0
      RETURN
!
!     End of DLACON
!
      END SUBROUTINE DLACON


      SUBROUTINE DLASSQ( N, X, INCX, SCALE, SUMSQ )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   SCALE, SUMSQ
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
!     ..
!
!  Purpose
!  =======
!
!  DLASSQ  returns the values  scl  and  smsq  such that
!
!     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
!
!  where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
!  assumed to be non-negative and  scl  returns the value
!
!     scl = max( scale, abs( x( i ) ) ).
!
!  scale and sumsq must be supplied in SCALE and SUMSQ and
!  scl and smsq are overwritten on SCALE and SUMSQ respectively.
!
!  The routine makes only one pass through the vector x.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of elements to be used from the vector X.
!
!  X       (input) DOUBLE PRECISION array, dimension (N)
!          The vector for which a scaled sum of squares is computed.
!             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
!
!  INCX    (input) INTEGER
!          The increment between successive values of the vector X.
!          INCX > 0.
!
!  SCALE   (input/output) DOUBLE PRECISION
!          On entry, the value  scale  in the equation above.
!          On exit, SCALE is overwritten with  scl , the scaling factor
!          for the sum of squares.
!
!  SUMSQ   (input/output) DOUBLE PRECISION
!          On entry, the value  sumsq  in the equation above.
!          On exit, SUMSQ is overwritten with  smsq , the basic sum of
!          squares from which  scl  has been factored out.
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            IX
      DOUBLE PRECISION   ABSXI
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS
!     ..
!     .. Executable Statements ..
!
      IF( N.GT.0 ) THEN
         DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
            IF( X( IX ).NE.ZERO ) THEN
               ABSXI = ABS( X( IX ) )
               IF( SCALE.LT.ABSXI ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SUMSQ = SUMSQ + ( ABSXI / SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
      RETURN
!
!     End of DLASSQ
!
      END SUBROUTINE DLASSQ


      SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DLASWP performs a series of row interchanges on the matrix A.
!  One row interchange is initiated for each of rows K1 through K2 of A.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the matrix of column dimension N to which the row
!          interchanges will be applied.
!          On exit, the permuted matrix.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!
!  K1      (input) INTEGER
!          The first element of IPIV for which a row interchange will
!          be done.
!
!  K2      (input) INTEGER
!          The last element of IPIV for which a row interchange will
!          be done.
!
!  IPIV    (input) INTEGER array, dimension (K2*abs(INCX))
!          The vector of pivot indices.  Only the elements in positions
!          K1 through K2 of IPIV are accessed.
!          IPIV(K) = L implies rows K and L are to be interchanged.
!
!  INCX    (input) INTEGER
!          The increment between successive values of IPIV.  If IPIV
!          is negative, the pivots are applied in reverse order.
!
!  Further Details
!  ===============
!
!  Modified by
!   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
      DOUBLE PRECISION   TEMP
!     ..
!     .. Executable Statements ..
!
!     Interchange row I with row IPIV(I) for each of rows K1 through K2.
!
      IF( INCX.GT.0 ) THEN
         IX0 = K1
         I1 = K1
         I2 = K2
         INC = 1
      ELSE IF( INCX.LT.0 ) THEN
         IX0 = 1 + ( 1-K2 )*INCX
         I1 = K2
         I2 = K1
         INC = -1
      ELSE
         RETURN
      END IF
!
      N32 = ( N / 32 )*32
      IF( N32.NE.0 ) THEN
         DO 30 J = 1, N32, 32
            IX = IX0
            DO 20 I = I1, I2, INC
               IP = IPIV( IX )
               IF( IP.NE.I ) THEN
                  DO 10 K = J, J + 31
                     TEMP = A( I, K )
                     A( I, K ) = A( IP, K )
                     A( IP, K ) = TEMP
   10             CONTINUE
               END IF
               IX = IX + INCX
   20       CONTINUE
   30    CONTINUE
      END IF
      IF( N32.NE.N ) THEN
         N32 = N32 + 1
         IX = IX0
         DO 50 I = I1, I2, INC
            IP = IPIV( IX )
            IF( IP.NE.I ) THEN
               DO 40 K = N32, N
                  TEMP = A( I, K )
                  A( I, K ) = A( IP, K )
                  A( IP, K ) = TEMP
   40          CONTINUE
            END IF
            IX = IX + INCX
   50    CONTINUE
      END IF
!
      RETURN
!
!     End of DLASWP
!
      END SUBROUTINE DLASWP


      SUBROUTINE DLATBS( UPLO, TRANS, DIAG, NORMIN, N, KD, AB, LDAB, X, &
                         SCALE, CNORM, INFO )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          DIAG, NORMIN, TRANS, UPLO
      INTEGER            INFO, KD, LDAB, N
      DOUBLE PRECISION   SCALE
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   AB( LDAB, * ), CNORM( * ), X( * )
!     ..
!
!  Purpose
!  =======
!
!  DLATBS solves one of the triangular systems
!
!     A *x = s*b  or  A'*x = s*b
!
!  with scaling to prevent overflow, where A is an upper or lower
!  triangular band matrix.  Here A' denotes the transpose of A, x and b
!  are n-element vectors, and s is a scaling factor, usually less than
!  or equal to 1, chosen so that the components of x will be less than
!  the overflow threshold.  If the unscaled problem will not cause
!  overflow, the Level 2 BLAS routine DTBSV is called.  If the matrix A
!  is singular (A(j,j) = 0 for some j), then s is set to 0 and a
!  non-trivial solution to A*x = 0 is returned.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the matrix A is upper or lower triangular.
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  TRANS   (input) CHARACTER*1
!          Specifies the operation applied to A.
!          = 'N':  Solve A * x = s*b  (No transpose)
!          = 'T':  Solve A'* x = s*b  (Transpose)
!          = 'C':  Solve A'* x = s*b  (Conjugate transpose = Transpose)
!
!  DIAG    (input) CHARACTER*1
!          Specifies whether or not the matrix A is unit triangular.
!          = 'N':  Non-unit triangular
!          = 'U':  Unit triangular
!
!  NORMIN  (input) CHARACTER*1
!          Specifies whether CNORM has been set or not.
!          = 'Y':  CNORM contains the column norms on entry
!          = 'N':  CNORM is not set on entry.  On exit, the norms will
!                  be computed and stored in CNORM.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  KD      (input) INTEGER
!          The number of subdiagonals or superdiagonals in the
!          triangular matrix A.  KD >= 0.
!
!  AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)
!          The upper or lower triangular band matrix A, stored in the
!          first KD+1 rows of the array. The j-th column of A is stored
!          in the j-th column of the array AB as follows:
!          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= KD+1.
!
!  X       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the right hand side b of the triangular system.
!          On exit, X is overwritten by the solution vector x.
!
!  SCALE   (output) DOUBLE PRECISION
!          The scaling factor s for the triangular system
!             A * x = s*b  or  A'* x = s*b.
!          If SCALE = 0, the matrix A is singular or badly scaled, and
!          the vector x is an exact or approximate solution to A*x = 0.
!
!  CNORM   (input or output) DOUBLE PRECISION array, dimension (N)
!
!          If NORMIN = 'Y', CNORM is an input argument and CNORM(j)
!          contains the norm of the off-diagonal part of the j-th column
!          of A.  If TRANS = 'N', CNORM(j) must be greater than or equal
!          to the infinity-norm, and if TRANS = 'T' or 'C', CNORM(j)
!          must be greater than or equal to the 1-norm.
!
!          If NORMIN = 'N', CNORM is an output argument and CNORM(j)
!          returns the 1-norm of the offdiagonal part of the j-th column
!          of A.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -k, the k-th argument had an illegal value
!
!  Further Details
!  ======= =======
!
!  A rough bound on x is computed; if that is less than overflow, DTBSV
!  is called, otherwise, specific code is used which checks for possible
!  overflow or divide-by-zero at every operation.
!
!  A columnwise scheme is used for solving A*x = b.  The basic algorithm
!  if A is lower triangular is
!
!       x[1:n] := b[1:n]
!       for j = 1, ..., n
!            x(j) := x(j) / A(j,j)
!            x[j+1:n] := x[j+1:n] - x(j) * A[j+1:n,j]
!       end
!
!  Define bounds on the components of x after j iterations of the loop:
!     M(j) = bound on x[1:j]
!     G(j) = bound on x[j+1:n]
!  Initially, let M(0) = 0 and G(0) = max{x(i), i=1,...,n}.
!
!  Then for iteration j+1 we have
!     M(j+1) <= G(j) / | A(j+1,j+1) |
!     G(j+1) <= G(j) + M(j+1) * | A[j+2:n,j+1] |
!            <= G(j) ( 1 + CNORM(j+1) / | A(j+1,j+1) | )
!
!  where CNORM(j+1) is greater than or equal to the infinity-norm of
!  column j+1 of A, not counting the diagonal.  Hence
!
!     G(j) <= G(0) product ( 1 + CNORM(i) / | A(i,i) | )
!                  1<=i<=j
!  and
!
!     |x(j)| <= ( G(0) / |A(j,j)| ) product ( 1 + CNORM(i) / |A(i,i)| )
!                                   1<=i< j
!
!  Since |x(j)| <= M(j), we use the Level 2 BLAS routine DTBSV if the
!  reciprocal of the largest M(j), j=1,..,n, is larger than
!  max(underflow, 1/overflow).
!
!  The bound on x(j) is also used to determine when a step in the
!  columnwise method can be performed without fear of overflow.  If
!  the computed bound is greater than a large constant, x is scaled to
!  prevent overflow, but if the bound overflows, x is set to 0, x(j) to
!  1, and scale to 0, and a non-trivial solution to A*x = 0 is found.
!
!  Similarly, a row-wise scheme is used to solve A'*x = b.  The basic
!  algorithm for A upper triangular is
!
!       for j = 1, ..., n
!            x(j) := ( b(j) - A[1:j-1,j]' * x[1:j-1] ) / A(j,j)
!       end
!
!  We simultaneously compute two bounds
!       G(j) = bound on ( b(i) - A[1:i-1,i]' * x[1:i-1] ), 1<=i<=j
!       M(j) = bound on x(i), 1<=i<=j
!
!  The initial values are G(0) = 0, M(0) = max{b(i), i=1,..,n}, and we
!  add the constraint G(j) >= G(j-1) and M(j) >= M(j-1) for j >= 1.
!  Then the bound on x(j) is
!
!       M(j) <= M(j-1) * ( 1 + CNORM(j) ) / | A(j,j) |
!
!            <= M(0) * product ( ( 1 + CNORM(i) ) / |A(i,i)| )
!                      1<=i<=j
!
!  and we can safely call DTBSV if 1/M(n) and 1/G(n) are both greater
!  than max(underflow, 1/overflow).
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, HALF, ONE
      PARAMETER          ( ZERO = 0.0D+0, HALF = 0.5D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOTRAN, NOUNIT, UPPER
      INTEGER            I, IMAX, J, JFIRST, JINC, JLAST, JLEN, MAIND
      DOUBLE PRECISION   BIGNUM, GROW, REC, SMLNUM, SUMJ, TJJ, TJJS, &
                         TMAX, TSCAL, USCAL, XBND, XJ, XMAX
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      INTEGER            IDAMAX
!      DOUBLE PRECISION   DASUM, DDOT, DLAMCH
!      EXTERNAL           LSAME, IDAMAX, DASUM, DDOT, DLAMCH
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DAXPY, DSCAL, DTBSV, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
!     ..
!     .. Executable Statements ..
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOTRAN = LSAME( TRANS, 'N' )
      NOUNIT = LSAME( DIAG, 'N' )
!
!     Test the input parameters.
!
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. &
               LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -3
      ELSE IF( .NOT.LSAME( NORMIN, 'Y' ) .AND. .NOT. &
               LSAME( NORMIN, 'N' ) ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( KD.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLATBS', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
!
!     Determine machine dependent parameters to control overflow.
!
      SMLNUM = DLAMCH( 'Safe minimum' ) / DLAMCH( 'Precision' )
      BIGNUM = ONE / SMLNUM
      SCALE = ONE
!
      IF( LSAME( NORMIN, 'N' ) ) THEN
!
!        Compute the 1-norm of each column, not including the diagonal.
!
         IF( UPPER ) THEN
!
!           A is upper triangular.
!
            DO 10 J = 1, N
               JLEN = MIN( KD, J-1 )
               CNORM( J ) = DASUM( JLEN, AB( KD+1-JLEN, J ), 1 )
   10       CONTINUE
         ELSE
!
!           A is lower triangular.
!
            DO 20 J = 1, N
               JLEN = MIN( KD, N-J )
               IF( JLEN.GT.0 ) THEN
                  CNORM( J ) = DASUM( JLEN, AB( 2, J ), 1 )
               ELSE
                  CNORM( J ) = ZERO
               END IF
   20       CONTINUE
         END IF
      END IF
!
!     Scale the column norms by TSCAL if the maximum element in CNORM is
!     greater than BIGNUM.
!
      IMAX = IDAMAX( N, CNORM, 1 )
      TMAX = CNORM( IMAX )
      IF( TMAX.LE.BIGNUM ) THEN
         TSCAL = ONE
      ELSE
         TSCAL = ONE / ( SMLNUM*TMAX )
         CALL DSCAL( N, TSCAL, CNORM, 1 )
      END IF
!
!     Compute a bound on the computed solution vector to see if the
!     Level 2 BLAS routine DTBSV can be used.
!
      J = IDAMAX( N, X, 1 )
      XMAX = ABS( X( J ) )
      XBND = XMAX
      IF( NOTRAN ) THEN
!
!        Compute the growth in A * x = b.
!
         IF( UPPER ) THEN
            JFIRST = N
            JLAST = 1
            JINC = -1
            MAIND = KD + 1
         ELSE
            JFIRST = 1
            JLAST = N
            JINC = 1
            MAIND = 1
         END IF
!
         IF( TSCAL.NE.ONE ) THEN
            GROW = ZERO
            GO TO 50
         END IF
!
         IF( NOUNIT ) THEN
!
!           A is non-unit triangular.
!
!           Compute GROW = 1/G(j) and XBND = 1/M(j).
!           Initially, G(0) = max{x(i), i=1,...,n}.
!
            GROW = ONE / MAX( XBND, SMLNUM )
            XBND = GROW
            DO 30 J = JFIRST, JLAST, JINC
!
!              Exit the loop if the growth factor is too small.
!
               IF( GROW.LE.SMLNUM ) &
                  GO TO 50
!
!              M(j) = G(j-1) / abs(A(j,j))
!
               TJJ = ABS( AB( MAIND, J ) )
               XBND = MIN( XBND, MIN( ONE, TJJ )*GROW )
               IF( TJJ+CNORM( J ).GE.SMLNUM ) THEN
!
!                 G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) )
!
                  GROW = GROW*( TJJ / ( TJJ+CNORM( J ) ) )
               ELSE
!
!                 G(j) could overflow, set GROW to 0.
!
                  GROW = ZERO
               END IF
   30       CONTINUE
            GROW = XBND
         ELSE
!
!           A is unit triangular.
!
!           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
!
            GROW = MIN( ONE, ONE / MAX( XBND, SMLNUM ) )
            DO 40 J = JFIRST, JLAST, JINC
!
!              Exit the loop if the growth factor is too small.
!
               IF( GROW.LE.SMLNUM ) &
                  GO TO 50
!
!              G(j) = G(j-1)*( 1 + CNORM(j) )
!
               GROW = GROW*( ONE / ( ONE+CNORM( J ) ) )
   40       CONTINUE
         END IF
   50    CONTINUE
!
      ELSE
!
!        Compute the growth in A' * x = b.
!
         IF( UPPER ) THEN
            JFIRST = 1
            JLAST = N
            JINC = 1
            MAIND = KD + 1
         ELSE
            JFIRST = N
            JLAST = 1
            JINC = -1
            MAIND = 1
         END IF
!
         IF( TSCAL.NE.ONE ) THEN
            GROW = ZERO
            GO TO 80
         END IF
!
         IF( NOUNIT ) THEN
!
!           A is non-unit triangular.
!
!           Compute GROW = 1/G(j) and XBND = 1/M(j).
!           Initially, M(0) = max{x(i), i=1,...,n}.
!
            GROW = ONE / MAX( XBND, SMLNUM )
            XBND = GROW
            DO 60 J = JFIRST, JLAST, JINC
!
!              Exit the loop if the growth factor is too small.
!
               IF( GROW.LE.SMLNUM ) &
                  GO TO 80
!
!              G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )
!
               XJ = ONE + CNORM( J )
               GROW = MIN( GROW, XBND / XJ )
!
!              M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))
!
               TJJ = ABS( AB( MAIND, J ) )
               IF( XJ.GT.TJJ ) &
                  XBND = XBND*( TJJ / XJ )
   60       CONTINUE
            GROW = MIN( GROW, XBND )
         ELSE
!
!           A is unit triangular.
!
!           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
!
            GROW = MIN( ONE, ONE / MAX( XBND, SMLNUM ) )
            DO 70 J = JFIRST, JLAST, JINC
!
!              Exit the loop if the growth factor is too small.
!
               IF( GROW.LE.SMLNUM ) &
                  GO TO 80
!
!              G(j) = ( 1 + CNORM(j) )*G(j-1)
!
               XJ = ONE + CNORM( J )
               GROW = GROW / XJ
   70       CONTINUE
         END IF
   80    CONTINUE
      END IF
!
      IF( ( GROW*TSCAL ).GT.SMLNUM ) THEN
!
!        Use the Level 2 BLAS solve if the reciprocal of the bound on
!        elements of X is not too small.
!
         CALL DTBSV( UPLO, TRANS, DIAG, N, KD, AB, LDAB, X, 1 )
      ELSE
!
!        Use a Level 1 BLAS solve, scaling intermediate results.
!
         IF( XMAX.GT.BIGNUM ) THEN
!
!           Scale X so that its components are less than or equal to
!           BIGNUM in absolute value.
!
            SCALE = BIGNUM / XMAX
            CALL DSCAL( N, SCALE, X, 1 )
            XMAX = BIGNUM
         END IF
!
         IF( NOTRAN ) THEN
!
!           Solve A * x = b
!
            DO 110 J = JFIRST, JLAST, JINC
!
!              Compute x(j) = b(j) / A(j,j), scaling x if necessary.
!
               XJ = ABS( X( J ) )
               IF( NOUNIT ) THEN
                  TJJS = AB( MAIND, J )*TSCAL
               ELSE
                  TJJS = TSCAL
                  IF( TSCAL.EQ.ONE ) &
                     GO TO 100
               END IF
               TJJ = ABS( TJJS )
               IF( TJJ.GT.SMLNUM ) THEN
!
!                    abs(A(j,j)) > SMLNUM:
!
                  IF( TJJ.LT.ONE ) THEN
                     IF( XJ.GT.TJJ*BIGNUM ) THEN
!
!                          Scale x by 1/b(j).
!
                        REC = ONE / XJ
                        CALL DSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
                  X( J ) = X( J ) / TJJS
                  XJ = ABS( X( J ) )
               ELSE IF( TJJ.GT.ZERO ) THEN
!
!                    0 < abs(A(j,j)) <= SMLNUM:
!
                  IF( XJ.GT.TJJ*BIGNUM ) THEN
!
!                       Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM
!                       to avoid overflow when dividing by A(j,j).
!
                     REC = ( TJJ*BIGNUM ) / XJ
                     IF( CNORM( J ).GT.ONE ) THEN
!
!                          Scale by 1/CNORM(j) to avoid overflow when
!                          multiplying x(j) times column j.
!
                        REC = REC / CNORM( J )
                     END IF
                     CALL DSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  END IF
                  X( J ) = X( J ) / TJJS
                  XJ = ABS( X( J ) )
               ELSE
!
!                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
!                    scale = 0, and compute a solution to A*x = 0.
!
                  DO 90 I = 1, N
                     X( I ) = ZERO
   90             CONTINUE
                  X( J ) = ONE
                  XJ = ONE
                  SCALE = ZERO
                  XMAX = ZERO
               END IF
  100          CONTINUE
!
!              Scale x if necessary to avoid overflow when adding a
!              multiple of column j of A.
!
               IF( XJ.GT.ONE ) THEN
                  REC = ONE / XJ
                  IF( CNORM( J ).GT.( BIGNUM-XMAX )*REC ) THEN
!
!                    Scale x by 1/(2*abs(x(j))).
!
                     REC = REC*HALF
                     CALL DSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                  END IF
               ELSE IF( XJ*CNORM( J ).GT.( BIGNUM-XMAX ) ) THEN
!
!                 Scale x by 1/2.
!
                  CALL DSCAL( N, HALF, X, 1 )
                  SCALE = SCALE*HALF
               END IF
!
               IF( UPPER ) THEN
                  IF( J.GT.1 ) THEN
!
!                    Compute the update
!                       x(max(1,j-kd):j-1) := x(max(1,j-kd):j-1) -
!                                             x(j)* A(max(1,j-kd):j-1,j)
!
                     JLEN = MIN( KD, J-1 )
                     CALL DAXPY( JLEN, -X( J )*TSCAL, &
                                 AB( KD+1-JLEN, J ), 1, X( J-JLEN ), 1 )
                     I = IDAMAX( J-1, X, 1 )
                     XMAX = ABS( X( I ) )
                  END IF
               ELSE IF( J.LT.N ) THEN
!
!                 Compute the update
!                    x(j+1:min(j+kd,n)) := x(j+1:min(j+kd,n)) -
!                                          x(j) * A(j+1:min(j+kd,n),j)
!
                  JLEN = MIN( KD, N-J )
                  IF( JLEN.GT.0 ) &
                     CALL DAXPY( JLEN, -X( J )*TSCAL, AB( 2, J ), 1, &
                                 X( J+1 ), 1 )
                  I = J + IDAMAX( N-J, X( J+1 ), 1 )
                  XMAX = ABS( X( I ) )
               END IF
  110       CONTINUE
!
         ELSE
!
!           Solve A' * x = b
!
            DO 160 J = JFIRST, JLAST, JINC
!
!              Compute x(j) = b(j) - sum A(k,j)*x(k).
!                                    k<>j
!
               XJ = ABS( X( J ) )
               USCAL = TSCAL
               REC = ONE / MAX( XMAX, ONE )
               IF( CNORM( J ).GT.( BIGNUM-XJ )*REC ) THEN
!
!                 If x(j) could overflow, scale x by 1/(2*XMAX).
!
                  REC = REC*HALF
                  IF( NOUNIT ) THEN
                     TJJS = AB( MAIND, J )*TSCAL
                  ELSE
                     TJJS = TSCAL
                  END IF
                  TJJ = ABS( TJJS )
                  IF( TJJ.GT.ONE ) THEN
!
!                       Divide by A(j,j) when scaling x if A(j,j) > 1.
!
                     REC = MIN( ONE, REC*TJJ )
                     USCAL = USCAL / TJJS
                  END IF
                  IF( REC.LT.ONE ) THEN
                     CALL DSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  END IF
               END IF
!
               SUMJ = ZERO
               IF( USCAL.EQ.ONE ) THEN
!
!                 If the scaling needed for A in the dot product is 1,
!                 call DDOT to perform the dot product.
!
                  IF( UPPER ) THEN
                     JLEN = MIN( KD, J-1 )
                     SUMJ = DDOT( JLEN, AB( KD+1-JLEN, J ), 1, &
                            X( J-JLEN ), 1 )
                  ELSE
                     JLEN = MIN( KD, N-J )
                     IF( JLEN.GT.0 ) &
                        SUMJ = DDOT( JLEN, AB( 2, J ), 1, X( J+1 ), 1 )
                  END IF
               ELSE
!
!                 Otherwise, use in-line code for the dot product.
!
                  IF( UPPER ) THEN
                     JLEN = MIN( KD, J-1 )
                     DO 120 I = 1, JLEN
                        SUMJ = SUMJ + ( AB( KD+I-JLEN, J )*USCAL )* &
                               X( J-JLEN-1+I )
  120                CONTINUE
                  ELSE
                     JLEN = MIN( KD, N-J )
                     DO 130 I = 1, JLEN
                        SUMJ = SUMJ + ( AB( I+1, J )*USCAL )*X( J+I )
  130                CONTINUE
                  END IF
               END IF
!
               IF( USCAL.EQ.TSCAL ) THEN
!
!                 Compute x(j) := ( x(j) - sumj ) / A(j,j) if 1/A(j,j)
!                 was not used to scale the dotproduct.
!
                  X( J ) = X( J ) - SUMJ
                  XJ = ABS( X( J ) )
                  IF( NOUNIT ) THEN
!
!                    Compute x(j) = x(j) / A(j,j), scaling if necessary.
!
                     TJJS = AB( MAIND, J )*TSCAL
                  ELSE
                     TJJS = TSCAL
                     IF( TSCAL.EQ.ONE ) &
                        GO TO 150
                  END IF
                  TJJ = ABS( TJJS )
                  IF( TJJ.GT.SMLNUM ) THEN
!
!                       abs(A(j,j)) > SMLNUM:
!
                     IF( TJJ.LT.ONE ) THEN
                        IF( XJ.GT.TJJ*BIGNUM ) THEN
!
!                             Scale X by 1/abs(x(j)).
!
                           REC = ONE / XJ
                           CALL DSCAL( N, REC, X, 1 )
                           SCALE = SCALE*REC
                           XMAX = XMAX*REC
                        END IF
                     END IF
                     X( J ) = X( J ) / TJJS
                  ELSE IF( TJJ.GT.ZERO ) THEN
!
!                       0 < abs(A(j,j)) <= SMLNUM:
!
                     IF( XJ.GT.TJJ*BIGNUM ) THEN
!
!                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
!
                        REC = ( TJJ*BIGNUM ) / XJ
                        CALL DSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                     X( J ) = X( J ) / TJJS
                  ELSE
!
!                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
!                       scale = 0, and compute a solution to A'*x = 0.
!
                     DO 140 I = 1, N
                        X( I ) = ZERO
  140                CONTINUE
                     X( J ) = ONE
                     SCALE = ZERO
                     XMAX = ZERO
                  END IF
  150             CONTINUE
               ELSE
!
!                 Compute x(j) := x(j) / A(j,j) - sumj if the dot
!                 product has already been divided by 1/A(j,j).
!
                  X( J ) = X( J ) / TJJS - SUMJ
               END IF
               XMAX = MAX( XMAX, ABS( X( J ) ) )
  160       CONTINUE
         END IF
         SCALE = SCALE / TSCAL
      END IF
!
!     Scale the column norms by 1/TSCAL for return.
!
      IF( TSCAL.NE.ONE ) THEN
         CALL DSCAL( N, ONE / TSCAL, CNORM, 1 )
      END IF
!
      RETURN
!
!     End of DLATBS
!
      END SUBROUTINE DLATBS


      SUBROUTINE DRSCL( N, SA, SX, INCX )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   SA
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   SX( * )
!     ..
!
!  Purpose
!  =======
!
!  DRSCL multiplies an n-element real vector x by the real scalar 1/a.
!  This is done without overflow or underflow as long as
!  the final result x/a does not overflow or underflow.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of components of the vector x.
!
!  SA      (input) DOUBLE PRECISION
!          The scalar a which is used to divide each component of x.
!          SA must be >= 0, or the subroutine will divide by zero.
!
!  SX      (input/output) DOUBLE PRECISION array, dimension
!                         (1+(N-1)*abs(INCX))
!          The n-element vector x.
!
!  INCX    (input) INTEGER
!          The increment between successive values of the vector SX.
!          > 0:  SX(1) = X(1) and SX(1+(i-1)*INCX) = x(i),     1< i<= n
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            DONE
      DOUBLE PRECISION   BIGNUM, CDEN, CDEN1, CNUM, CNUM1, MUL, SMLNUM
!     ..
!     .. External Functions ..
!      DOUBLE PRECISION   DLAMCH
!      EXTERNAL           DLAMCH
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DSCAL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( N.LE.0 ) &
         RETURN
!
!     Get machine parameters
!
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
!
!     Initialize the denominator to SA and the numerator to 1.
!
      CDEN = SA
      CNUM = ONE
!
   10 CONTINUE
      CDEN1 = CDEN*SMLNUM
      CNUM1 = CNUM / BIGNUM
      IF( ABS( CDEN1 ).GT.ABS( CNUM ) .AND. CNUM.NE.ZERO ) THEN
!
!        Pre-multiply X by SMLNUM if CDEN is large compared to CNUM.
!
         MUL = SMLNUM
         DONE = .FALSE.
         CDEN = CDEN1
      ELSE IF( ABS( CNUM1 ).GT.ABS( CDEN ) ) THEN
!
!        Pre-multiply X by BIGNUM if CDEN is small compared to CNUM.
!
         MUL = BIGNUM
         DONE = .FALSE.
         CNUM = CNUM1
      ELSE
!
!        Multiply X by CNUM / CDEN and return.
!
         MUL = CNUM / CDEN
         DONE = .TRUE.
      END IF
!
!     Scale the vector X by MUL
!
      CALL DSCAL( N, MUL, SX, INCX )
!
      IF( .NOT.DONE ) &
         GO TO 10
!
      RETURN
!
!     End of DRSCL
!
      END SUBROUTINE DRSCL


      SUBROUTINE XERBLA( SRNAME, INFO )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER*6        SRNAME
      INTEGER            INFO
!     ..
!
!  Purpose
!  =======
!
!  XERBLA  is an error handler for the LAPACK routines.
!  It is called by an LAPACK routine if an input parameter has an
!  invalid value.  A message is printed and execution stops.
!
!  Installers may consider modifying the STOP statement in order to
!  call system-specific exception-handling facilities.
!
!  Arguments
!  =========
!
!  SRNAME  (input) CHARACTER*6
!          The name of the routine which called XERBLA.
!
!  INFO    (input) INTEGER
!          The position of the invalid parameter in the parameter list
!          of the calling routine.
!
! =====================================================================
!
!     .. Executable Statements ..
!
      WRITE( *, FMT = 9999 )SRNAME, INFO
!
      STOP
!
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ', &
            'an illegal value' )
!
!     End of XERBLA
!
      END SUBROUTINE XERBLA


      SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DGETF2 computes an LU factorization of a general m-by-n matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 2 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the m by n matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   SFMIN
      INTEGER            I, J, JP
!     ..
!     .. External Functions ..
!      DOUBLE PRECISION   DLAMCH
!      INTEGER            IDAMAX
!      EXTERNAL           DLAMCH, IDAMAX
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DGER, DSCAL, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETF2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 ) &
         RETURN
!
!     Compute machine safe minimum
!
      SFMIN = DLAMCH('S')
!
      DO 10 J = 1, MIN( M, N )
!
!        Find pivot and test for singularity.
!
         JP = J - 1 + IDAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN
!
!           Apply the interchange to columns 1:N.
!
            IF( JP.NE.J ) &
               CALL DSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
!
!           Compute elements J+1:M of J-th column.
!
            IF( J.LT.M ) THEN
               IF( ABS(A( J, J )) .GE. SFMIN ) THEN
                  CALL DSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
               ELSE
                 DO 20 I = 1, M-J
                    A( J+I, J ) = A( J+I, J ) / A( J, J )
   20            CONTINUE
               END IF
            END IF
!
         ELSE IF( INFO.EQ.0 ) THEN
!
            INFO = J
         END IF
!
         IF( J.LT.MIN( M, N ) ) THEN
!
!           Update trailing submatrix.
!
            CALL DGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA, &
                       A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN
!
!     End of DGETF2
!
      END SUBROUTINE DGETF2


      SUBROUTINE DLATRS( UPLO, TRANS, DIAG, NORMIN, N, A, LDA, X, SCALE, &
                         CNORM, INFO )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          DIAG, NORMIN, TRANS, UPLO
      INTEGER            INFO, LDA, N
      DOUBLE PRECISION   SCALE
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), CNORM( * ), X( * )
!     ..
!
!  Purpose
!  =======
!
!  DLATRS solves one of the triangular systems
!
!     A *x = s*b  or  A'*x = s*b
!
!  with scaling to prevent overflow.  Here A is an upper or lower
!  triangular matrix, A' denotes the transpose of A, x and b are
!  n-element vectors, and s is a scaling factor, usually less than
!  or equal to 1, chosen so that the components of x will be less than
!  the overflow threshold.  If the unscaled problem will not cause
!  overflow, the Level 2 BLAS routine DTRSV is called.  If the matrix A
!  is singular (A(j,j) = 0 for some j), then s is set to 0 and a
!  non-trivial solution to A*x = 0 is returned.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the matrix A is upper or lower triangular.
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  TRANS   (input) CHARACTER*1
!          Specifies the operation applied to A.
!          = 'N':  Solve A * x = s*b  (No transpose)
!          = 'T':  Solve A'* x = s*b  (Transpose)
!          = 'C':  Solve A'* x = s*b  (Conjugate transpose = Transpose)
!
!  DIAG    (input) CHARACTER*1
!          Specifies whether or not the matrix A is unit triangular.
!          = 'N':  Non-unit triangular
!          = 'U':  Unit triangular
!
!  NORMIN  (input) CHARACTER*1
!          Specifies whether CNORM has been set or not.
!          = 'Y':  CNORM contains the column norms on entry
!          = 'N':  CNORM is not set on entry.  On exit, the norms will
!                  be computed and stored in CNORM.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!          The triangular matrix A.  If UPLO = 'U', the leading n by n
!          upper triangular part of the array A contains the upper
!          triangular matrix, and the strictly lower triangular part of
!          A is not referenced.  If UPLO = 'L', the leading n by n lower
!          triangular part of the array A contains the lower triangular
!          matrix, and the strictly upper triangular part of A is not
!          referenced.  If DIAG = 'U', the diagonal elements of A are
!          also not referenced and are assumed to be 1.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max (1,N).
!
!  X       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the right hand side b of the triangular system.
!          On exit, X is overwritten by the solution vector x.
!
!  SCALE   (output) DOUBLE PRECISION
!          The scaling factor s for the triangular system
!             A * x = s*b  or  A'* x = s*b.
!          If SCALE = 0, the matrix A is singular or badly scaled, and
!          the vector x is an exact or approximate solution to A*x = 0.
!
!  CNORM   (input or output) DOUBLE PRECISION array, dimension (N)
!
!          If NORMIN = 'Y', CNORM is an input argument and CNORM(j)
!          contains the norm of the off-diagonal part of the j-th column
!          of A.  If TRANS = 'N', CNORM(j) must be greater than or equal
!          to the infinity-norm, and if TRANS = 'T' or 'C', CNORM(j)
!          must be greater than or equal to the 1-norm.
!
!          If NORMIN = 'N', CNORM is an output argument and CNORM(j)
!          returns the 1-norm of the offdiagonal part of the j-th column
!          of A.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -k, the k-th argument had an illegal value
!
!  Further Details
!  ======= =======
!
!  A rough bound on x is computed; if that is less than overflow, DTRSV
!  is called, otherwise, specific code is used which checks for possible
!  overflow or divide-by-zero at every operation.
!
!  A columnwise scheme is used for solving A*x = b.  The basic algorithm
!  if A is lower triangular is
!
!       x[1:n] := b[1:n]
!       for j = 1, ..., n
!            x(j) := x(j) / A(j,j)
!            x[j+1:n] := x[j+1:n] - x(j) * A[j+1:n,j]
!       end
!
!  Define bounds on the components of x after j iterations of the loop:
!     M(j) = bound on x[1:j]
!     G(j) = bound on x[j+1:n]
!  Initially, let M(0) = 0 and G(0) = max{x(i), i=1,...,n}.
!
!  Then for iteration j+1 we have
!     M(j+1) <= G(j) / | A(j+1,j+1) |
!     G(j+1) <= G(j) + M(j+1) * | A[j+2:n,j+1] |
!            <= G(j) ( 1 + CNORM(j+1) / | A(j+1,j+1) | )
!
!  where CNORM(j+1) is greater than or equal to the infinity-norm of
!  column j+1 of A, not counting the diagonal.  Hence
!
!     G(j) <= G(0) product ( 1 + CNORM(i) / | A(i,i) | )
!                  1<=i<=j
!  and
!
!     |x(j)| <= ( G(0) / |A(j,j)| ) product ( 1 + CNORM(i) / |A(i,i)| )
!                                   1<=i< j
!
!  Since |x(j)| <= M(j), we use the Level 2 BLAS routine DTRSV if the
!  reciprocal of the largest M(j), j=1,..,n, is larger than
!  max(underflow, 1/overflow).
!
!  The bound on x(j) is also used to determine when a step in the
!  columnwise method can be performed without fear of overflow.  If
!  the computed bound is greater than a large constant, x is scaled to
!  prevent overflow, but if the bound overflows, x is set to 0, x(j) to
!  1, and scale to 0, and a non-trivial solution to A*x = 0 is found.
!
!  Similarly, a row-wise scheme is used to solve A'*x = b.  The basic
!  algorithm for A upper triangular is
!
!       for j = 1, ..., n
!            x(j) := ( b(j) - A[1:j-1,j]' * x[1:j-1] ) / A(j,j)
!       end
!
!  We simultaneously compute two bounds
!       G(j) = bound on ( b(i) - A[1:i-1,i]' * x[1:i-1] ), 1<=i<=j
!       M(j) = bound on x(i), 1<=i<=j
!
!  The initial values are G(0) = 0, M(0) = max{b(i), i=1,..,n}, and we
!  add the constraint G(j) >= G(j-1) and M(j) >= M(j-1) for j >= 1.
!  Then the bound on x(j) is
!
!       M(j) <= M(j-1) * ( 1 + CNORM(j) ) / | A(j,j) |
!
!            <= M(0) * product ( ( 1 + CNORM(i) ) / |A(i,i)| )
!                      1<=i<=j
!
!  and we can safely call DTRSV if 1/M(n) and 1/G(n) are both greater
!  than max(underflow, 1/overflow).
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, HALF, ONE
      PARAMETER          ( ZERO = 0.0D+0, HALF = 0.5D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOTRAN, NOUNIT, UPPER
      INTEGER            I, IMAX, J, JFIRST, JINC, JLAST
      DOUBLE PRECISION   BIGNUM, GROW, REC, SMLNUM, SUMJ, TJJ, TJJS, &
                         TMAX, TSCAL, USCAL, XBND, XJ, XMAX
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      INTEGER            IDAMAX
!      DOUBLE PRECISION   DASUM, DDOT, DLAMCH
!      EXTERNAL           LSAME, IDAMAX, DASUM, DDOT, DLAMCH
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DAXPY, DSCAL, DTRSV, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
!     ..
!     .. Executable Statements ..
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOTRAN = LSAME( TRANS, 'N' )
      NOUNIT = LSAME( DIAG, 'N' )
!
!     Test the input parameters.
!
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. &
               LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -3
      ELSE IF( .NOT.LSAME( NORMIN, 'Y' ) .AND. .NOT. &
               LSAME( NORMIN, 'N' ) ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLATRS', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
!
!     Determine machine dependent parameters to control overflow.
!
      SMLNUM = DLAMCH( 'Safe minimum' ) / DLAMCH( 'Precision' )
      BIGNUM = ONE / SMLNUM
      SCALE = ONE
!
      IF( LSAME( NORMIN, 'N' ) ) THEN
!
!        Compute the 1-norm of each column, not including the diagonal.
!
         IF( UPPER ) THEN
!
!           A is upper triangular.
!
            DO 10 J = 1, N
               CNORM( J ) = DASUM( J-1, A( 1, J ), 1 )
   10       CONTINUE
         ELSE
!
!           A is lower triangular.
!
            DO 20 J = 1, N - 1
               CNORM( J ) = DASUM( N-J, A( J+1, J ), 1 )
   20       CONTINUE
            CNORM( N ) = ZERO
         END IF
      END IF
!
!     Scale the column norms by TSCAL if the maximum element in CNORM is
!     greater than BIGNUM.
!
      IMAX = IDAMAX( N, CNORM, 1 )
      TMAX = CNORM( IMAX )
      IF( TMAX.LE.BIGNUM ) THEN
         TSCAL = ONE
      ELSE
         TSCAL = ONE / ( SMLNUM*TMAX )
         CALL DSCAL( N, TSCAL, CNORM, 1 )
      END IF
!
!     Compute a bound on the computed solution vector to see if the
!     Level 2 BLAS routine DTRSV can be used.
!
      J = IDAMAX( N, X, 1 )
      XMAX = ABS( X( J ) )
      XBND = XMAX
      IF( NOTRAN ) THEN
!
!        Compute the growth in A * x = b.
!
         IF( UPPER ) THEN
            JFIRST = N
            JLAST = 1
            JINC = -1
         ELSE
            JFIRST = 1
            JLAST = N
            JINC = 1
         END IF
!
         IF( TSCAL.NE.ONE ) THEN
            GROW = ZERO
            GO TO 50
         END IF
!
         IF( NOUNIT ) THEN
!
!           A is non-unit triangular.
!
!           Compute GROW = 1/G(j) and XBND = 1/M(j).
!           Initially, G(0) = max{x(i), i=1,...,n}.
!
            GROW = ONE / MAX( XBND, SMLNUM )
            XBND = GROW
            DO 30 J = JFIRST, JLAST, JINC
!
!              Exit the loop if the growth factor is too small.
!
               IF( GROW.LE.SMLNUM ) &
                  GO TO 50
!
!              M(j) = G(j-1) / abs(A(j,j))
!
               TJJ = ABS( A( J, J ) )
               XBND = MIN( XBND, MIN( ONE, TJJ )*GROW )
               IF( TJJ+CNORM( J ).GE.SMLNUM ) THEN
!
!                 G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) )
!
                  GROW = GROW*( TJJ / ( TJJ+CNORM( J ) ) )
               ELSE
!
!                 G(j) could overflow, set GROW to 0.
!
                  GROW = ZERO
               END IF
   30       CONTINUE
            GROW = XBND
         ELSE
!
!           A is unit triangular.
!
!           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
!
            GROW = MIN( ONE, ONE / MAX( XBND, SMLNUM ) )
            DO 40 J = JFIRST, JLAST, JINC
!
!              Exit the loop if the growth factor is too small.
!
               IF( GROW.LE.SMLNUM ) &
                  GO TO 50
!
!              G(j) = G(j-1)*( 1 + CNORM(j) )
!
               GROW = GROW*( ONE / ( ONE+CNORM( J ) ) )
   40       CONTINUE
         END IF
   50    CONTINUE
!
      ELSE
!
!        Compute the growth in A' * x = b.
!
         IF( UPPER ) THEN
            JFIRST = 1
            JLAST = N
            JINC = 1
         ELSE
            JFIRST = N
            JLAST = 1
            JINC = -1
         END IF
!
         IF( TSCAL.NE.ONE ) THEN
            GROW = ZERO
            GO TO 80
         END IF
!
         IF( NOUNIT ) THEN
!
!           A is non-unit triangular.
!
!           Compute GROW = 1/G(j) and XBND = 1/M(j).
!           Initially, M(0) = max{x(i), i=1,...,n}.
!
            GROW = ONE / MAX( XBND, SMLNUM )
            XBND = GROW
            DO 60 J = JFIRST, JLAST, JINC
!
!              Exit the loop if the growth factor is too small.
!
               IF( GROW.LE.SMLNUM ) &
                  GO TO 80
!
!              G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )
!
               XJ = ONE + CNORM( J )
               GROW = MIN( GROW, XBND / XJ )
!
!              M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))
!
               TJJ = ABS( A( J, J ) )
               IF( XJ.GT.TJJ ) &
                  XBND = XBND*( TJJ / XJ )
   60       CONTINUE
            GROW = MIN( GROW, XBND )
         ELSE
!
!           A is unit triangular.
!
!           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
!
            GROW = MIN( ONE, ONE / MAX( XBND, SMLNUM ) )
            DO 70 J = JFIRST, JLAST, JINC
!
!              Exit the loop if the growth factor is too small.
!
               IF( GROW.LE.SMLNUM ) &
                  GO TO 80
!
!              G(j) = ( 1 + CNORM(j) )*G(j-1)
!
               XJ = ONE + CNORM( J )
               GROW = GROW / XJ
   70       CONTINUE
         END IF
   80    CONTINUE
      END IF
!
      IF( ( GROW*TSCAL ).GT.SMLNUM ) THEN
!
!        Use the Level 2 BLAS solve if the reciprocal of the bound on
!        elements of X is not too small.
!
         CALL DTRSV( UPLO, TRANS, DIAG, N, A, LDA, X, 1 )
      ELSE
!
!        Use a Level 1 BLAS solve, scaling intermediate results.
!
         IF( XMAX.GT.BIGNUM ) THEN
!
!           Scale X so that its components are less than or equal to
!           BIGNUM in absolute value.
!
            SCALE = BIGNUM / XMAX
            CALL DSCAL( N, SCALE, X, 1 )
            XMAX = BIGNUM
         END IF
!
         IF( NOTRAN ) THEN
!
!           Solve A * x = b
!
            DO 110 J = JFIRST, JLAST, JINC
!
!              Compute x(j) = b(j) / A(j,j), scaling x if necessary.
!
               XJ = ABS( X( J ) )
               IF( NOUNIT ) THEN
                  TJJS = A( J, J )*TSCAL
               ELSE
                  TJJS = TSCAL
                  IF( TSCAL.EQ.ONE ) &
                     GO TO 100
               END IF
               TJJ = ABS( TJJS )
               IF( TJJ.GT.SMLNUM ) THEN
!
!                    abs(A(j,j)) > SMLNUM:
!
                  IF( TJJ.LT.ONE ) THEN
                     IF( XJ.GT.TJJ*BIGNUM ) THEN
!
!                          Scale x by 1/b(j).
!
                        REC = ONE / XJ
                        CALL DSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
                  X( J ) = X( J ) / TJJS
                  XJ = ABS( X( J ) )
               ELSE IF( TJJ.GT.ZERO ) THEN
!
!                    0 < abs(A(j,j)) <= SMLNUM:
!
                  IF( XJ.GT.TJJ*BIGNUM ) THEN
!
!                       Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM
!                       to avoid overflow when dividing by A(j,j).
!
                     REC = ( TJJ*BIGNUM ) / XJ
                     IF( CNORM( J ).GT.ONE ) THEN
!
!                          Scale by 1/CNORM(j) to avoid overflow when
!                          multiplying x(j) times column j.
!
                        REC = REC / CNORM( J )
                     END IF
                     CALL DSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  END IF
                  X( J ) = X( J ) / TJJS
                  XJ = ABS( X( J ) )
               ELSE
!
!                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
!                    scale = 0, and compute a solution to A*x = 0.
!
                  DO 90 I = 1, N
                     X( I ) = ZERO
   90             CONTINUE
                  X( J ) = ONE
                  XJ = ONE
                  SCALE = ZERO
                  XMAX = ZERO
               END IF
  100          CONTINUE
!
!              Scale x if necessary to avoid overflow when adding a
!              multiple of column j of A.
!
               IF( XJ.GT.ONE ) THEN
                  REC = ONE / XJ
                  IF( CNORM( J ).GT.( BIGNUM-XMAX )*REC ) THEN
!
!                    Scale x by 1/(2*abs(x(j))).
!
                     REC = REC*HALF
                     CALL DSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                  END IF
               ELSE IF( XJ*CNORM( J ).GT.( BIGNUM-XMAX ) ) THEN
!
!                 Scale x by 1/2.
!
                  CALL DSCAL( N, HALF, X, 1 )
                  SCALE = SCALE*HALF
               END IF
!
               IF( UPPER ) THEN
                  IF( J.GT.1 ) THEN
!
!                    Compute the update
!                       x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j)
!
                     CALL DAXPY( J-1, -X( J )*TSCAL, A( 1, J ), 1, X, &
                                 1 )
                     I = IDAMAX( J-1, X, 1 )
                     XMAX = ABS( X( I ) )
                  END IF
               ELSE
                  IF( J.LT.N ) THEN
!
!                    Compute the update
!                       x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j)
!
                     CALL DAXPY( N-J, -X( J )*TSCAL, A( J+1, J ), 1, &
                                 X( J+1 ), 1 )
                     I = J + IDAMAX( N-J, X( J+1 ), 1 )
                     XMAX = ABS( X( I ) )
                  END IF
               END IF
  110       CONTINUE
!
         ELSE
!
!           Solve A' * x = b
!
            DO 160 J = JFIRST, JLAST, JINC
!
!              Compute x(j) = b(j) - sum A(k,j)*x(k).
!                                    k<>j
!
               XJ = ABS( X( J ) )
               USCAL = TSCAL
               REC = ONE / MAX( XMAX, ONE )
               IF( CNORM( J ).GT.( BIGNUM-XJ )*REC ) THEN
!
!                 If x(j) could overflow, scale x by 1/(2*XMAX).
!
                  REC = REC*HALF
                  IF( NOUNIT ) THEN
                     TJJS = A( J, J )*TSCAL
                  ELSE
                     TJJS = TSCAL
                  END IF
                  TJJ = ABS( TJJS )
                  IF( TJJ.GT.ONE ) THEN
!
!                       Divide by A(j,j) when scaling x if A(j,j) > 1.
!
                     REC = MIN( ONE, REC*TJJ )
                     USCAL = USCAL / TJJS
                  END IF
                  IF( REC.LT.ONE ) THEN
                     CALL DSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  END IF
               END IF
!
               SUMJ = ZERO
               IF( USCAL.EQ.ONE ) THEN
!
!                 If the scaling needed for A in the dot product is 1,
!                 call DDOT to perform the dot product.
!
                  IF( UPPER ) THEN
                     SUMJ = DDOT( J-1, A( 1, J ), 1, X, 1 )
                  ELSE IF( J.LT.N ) THEN
                     SUMJ = DDOT( N-J, A( J+1, J ), 1, X( J+1 ), 1 )
                  END IF
               ELSE
!
!                 Otherwise, use in-line code for the dot product.
!
                  IF( UPPER ) THEN
                     DO 120 I = 1, J - 1
                        SUMJ = SUMJ + ( A( I, J )*USCAL )*X( I )
  120                CONTINUE
                  ELSE IF( J.LT.N ) THEN
                     DO 130 I = J + 1, N
                        SUMJ = SUMJ + ( A( I, J )*USCAL )*X( I )
  130                CONTINUE
                  END IF
               END IF
!
               IF( USCAL.EQ.TSCAL ) THEN
!
!                 Compute x(j) := ( x(j) - sumj ) / A(j,j) if 1/A(j,j)
!                 was not used to scale the dotproduct.
!
                  X( J ) = X( J ) - SUMJ
                  XJ = ABS( X( J ) )
                  IF( NOUNIT ) THEN
                     TJJS = A( J, J )*TSCAL
                  ELSE
                     TJJS = TSCAL
                     IF( TSCAL.EQ.ONE ) &
                        GO TO 150
                  END IF
!
!                    Compute x(j) = x(j) / A(j,j), scaling if necessary.
!
                  TJJ = ABS( TJJS )
                  IF( TJJ.GT.SMLNUM ) THEN
!
!                       abs(A(j,j)) > SMLNUM:
!
                     IF( TJJ.LT.ONE ) THEN
                        IF( XJ.GT.TJJ*BIGNUM ) THEN
!
!                             Scale X by 1/abs(x(j)).
!
                           REC = ONE / XJ
                           CALL DSCAL( N, REC, X, 1 )
                           SCALE = SCALE*REC
                           XMAX = XMAX*REC
                        END IF
                     END IF
                     X( J ) = X( J ) / TJJS
                  ELSE IF( TJJ.GT.ZERO ) THEN
!
!                       0 < abs(A(j,j)) <= SMLNUM:
!
                     IF( XJ.GT.TJJ*BIGNUM ) THEN
!
!                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
!
                        REC = ( TJJ*BIGNUM ) / XJ
                        CALL DSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                     X( J ) = X( J ) / TJJS
                  ELSE
!
!                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
!                       scale = 0, and compute a solution to A'*x = 0.
!
                     DO 140 I = 1, N
                        X( I ) = ZERO
  140                CONTINUE
                     X( J ) = ONE
                     SCALE = ZERO
                     XMAX = ZERO
                  END IF
  150             CONTINUE
               ELSE
!
!                 Compute x(j) := x(j) / A(j,j)  - sumj if the dot
!                 product has already been divided by 1/A(j,j).
!
                  X( J ) = X( J ) / TJJS - SUMJ
               END IF
               XMAX = MAX( XMAX, ABS( X( J ) ) )
  160       CONTINUE
         END IF
         SCALE = SCALE / TSCAL
      END IF
!
!     Scale the column norms by 1/TSCAL for return.
!
      IF( TSCAL.NE.ONE ) THEN
         CALL DSCAL( N, ONE / TSCAL, CNORM, 1 )
      END IF
!
      RETURN
!
!     End of DLATRS
!
      END SUBROUTINE DLATRS


      INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
!
!  -- LAPACK auxiliary routine (version 3.1.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     January 2007
!
!     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
!     ..
!
!  Purpose
!  =======
!
!  ILAENV is called from the LAPACK routines to choose problem-dependent
!  parameters for the local environment.  See ISPEC for a description of
!  the parameters.
!
!  ILAENV returns an INTEGER
!  if ILAENV >= 0: ILAENV returns the value of the parameter specified b
!  if ILAENV < 0:  if ILAENV = -k, the k-th argument had an illegal valu
!
!  This version provides a set of parameters which should give good,
!  but not optimal, performance on many of the currently available
!  computers.  Users are encouraged to modify this subroutine to set
!  the tuning parameters for their particular machine using the option
!  and problem size information in the arguments.
!
!  This routine will not function correctly if it is converted to all
!  lower case.  Converting it to all upper case is allowed.
!
!  Arguments
!  =========
!
!  ISPEC   (input) INTEGER
!          Specifies the parameter to be returned as the value of
!          ILAENV.
!          = 1: the optimal blocksize; if this value is 1, an unblocked
!               algorithm will give the best performance.
!          = 2: the minimum block size for which the block routine
!               should be used; if the usable block size is less than
!               this value, an unblocked routine should be used.
!          = 3: the crossover point (in a block routine, for N less
!               than this value, an unblocked routine should be used)
!          = 4: the number of shifts, used in the nonsymmetric
!               eigenvalue routines (DEPRECATED)
!          = 5: the minimum column dimension for blocking to be used;
!               rectangular blocks must have dimension at least k by m,
!               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!          = 6: the crossover point for the SVD (when reducing an m by n
!               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!               this value, a QR factorization is used first to reduce
!               the matrix to a triangular form.)
!          = 7: the number of processors
!          = 8: the crossover point for the multishift QR method
!               for nonsymmetric eigenvalue problems (DEPRECATED)
!          = 9: maximum size of the subproblems at the bottom of the
!               computation tree in the divide-and-conquer algorithm
!               (used by xGELSD and xGESDD)
!          =10: ieee NaN arithmetic can be trusted not to trap
!          =11: infinity arithmetic can be trusted not to trap
!          12 <= ISPEC <= 16:
!               xHSEQR or one of its subroutines,
!               see IPARMQ for detailed explanation
!
!  NAME    (input) CHARACTER*(*)
!          The name of the calling subroutine, in either upper case or
!          lower case.
!
!  OPTS    (input) CHARACTER*(*)
!          The character options to the subroutine NAME, concatenated
!          into a single character string.  For example, UPLO = 'U',
!          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!          be specified as OPTS = 'UTN'.
!
!  N1      (input) INTEGER
!  N2      (input) INTEGER
!  N3      (input) INTEGER
!  N4      (input) INTEGER
!          Problem dimensions for the subroutine NAME; these may not all
!          be required.
!
!  Further Details
!  ===============
!
!  The following conventions have been used when calling ILAENV from the
!  LAPACK routines:
!  1)  OPTS is a concatenation of all of the character options to
!      subroutine NAME, in the same order that they appear in the
!      argument list for NAME, even if they are not used in determining
!      the value of the parameter specified by ISPEC.
!  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
!      that they appear in the argument list for NAME.  N1 is used
!      first, N2 second, and so on, and unused problem dimensions are
!      passed a value of -1.
!  3)  The parameter value returned by ILAENV is checked for validity in
!      the calling subroutine.  For example, ILAENV is used to retrieve
!      the optimal blocksize for STRTRI as follows:
!
!      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!      IF( NB.LE.1 ) NB = MAX( 1, N )
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, IC, IZ, NB, NBMIN, NX
      LOGICAL            CNAME, SNAME
      CHARACTER          C1*1, C2*2, C4*2, C3*3, SUBNAM*6
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
!     ..
!     .. External Functions ..
!      INTEGER            IEEECK, IPARMQ
!      EXTERNAL           IEEECK, IPARMQ
!     ..
!     .. Executable Statements ..
!
      GO TO ( 10, 10, 10, 80, 90, 100, 110, 120, &
              130, 140, 150, 160, 160, 160, 160, 160 )ISPEC
!
!     Invalid value for ISPEC
!
      ILAENV = -1
      RETURN
!
   10 CONTINUE
!
!     Convert NAME to upper case if the first character is lower case.
!
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1: 1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
!
!        ASCII character set
!
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( IC.GE.97 .AND. IC.LE.122 ) &
                  SUBNAM( I: I ) = CHAR( IC-32 )
   20       CONTINUE
         END IF
!
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
!
!        EBCDIC character set
!
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
             ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
             ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC+64 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
                   ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
                   ( IC.GE.162 .AND. IC.LE.169 ) )SUBNAM( I: &
                   I ) = CHAR( IC+64 )
   30       CONTINUE
         END IF
!
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
!
!        Prime machines:  ASCII+128
!
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 40 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( IC.GE.225 .AND. IC.LE.250 ) &
                  SUBNAM( I: I ) = CHAR( IC-32 )
   40       CONTINUE
         END IF
      END IF
!
      C1 = SUBNAM( 1: 1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) ) &
         RETURN
      C2 = SUBNAM( 2: 3 )
      C3 = SUBNAM( 4: 6 )
      C4 = C3( 2: 3 )
!
      GO TO ( 50, 60, 70 )ISPEC
!
   50 CONTINUE
!
!     ISPEC = 1:  block size
!
!     In these examples, separate code is provided for setting NB for
!     real and complex.  We assume that NB will take the same value in
!     single or double precision.
!
      NB = 1
!
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
                  C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
                 THEN
               NB = 32
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
                 THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
                 THEN
               NB = 32
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
                 THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
!
   60 CONTINUE
!
!     ISPEC = 2:  minimum block size
!
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ. &
             'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
                 THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
                 THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
                 THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
                 THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
!
   70 CONTINUE
!
!     ISPEC = 3:  crossover point
!
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ. &
             'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
                 THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
                 THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
!
   80 CONTINUE
!
!     ISPEC = 4:  number of shifts (used by xHSEQR)
!
      ILAENV = 6
      RETURN
!
   90 CONTINUE
!
!     ISPEC = 5:  minimum column dimension (not used)
!
      ILAENV = 2
      RETURN
!
  100 CONTINUE
!
!     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
!
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
!
  110 CONTINUE
!
!     ISPEC = 7:  number of processors (not used)
!
      ILAENV = 1
      RETURN
!
  120 CONTINUE
!
!     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
!
      ILAENV = 50
      RETURN
!
  130 CONTINUE
!
!     ISPEC = 9:  maximum size of the subproblems at the bottom of the
!                 computation tree in the divide-and-conquer algorithm
!                 (used by xGELSD and xGESDD)
!
      ILAENV = 25
      RETURN
!
  140 CONTINUE
!
!     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
!
!     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 0, 0.0, 1.0 )
      END IF
      RETURN
!
  150 CONTINUE
!
!     ISPEC = 11: infinity arithmetic can be trusted not to trap
!
!     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 1, 0.0, 1.0 )
      END IF
      RETURN
!
  160 CONTINUE
!
!     12 <= ISPEC <= 16: xHSEQR or one of its subroutines.
!
      ILAENV = IPARMQ( ISPEC, N2, N3 )
      RETURN
!
!     End of ILAENV
!
      END FUNCTION ILAENV


      LOGICAL FUNCTION LSAME( CA, CB )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          CA, CB
!     ..
!
!  Purpose
!  =======
!
!  LSAME returns .TRUE. if CA is the same letter as CB regardless of
!  case.
!
!  Arguments
!  =========
!
!  CA      (input) CHARACTER*1
!  CB      (input) CHARACTER*1
!          CA and CB specify the single characters to be compared.
!
! =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
!     ..
!     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
!     ..
!     .. Executable Statements ..
!
!     Test if the characters are equal
!
      LSAME = CA.EQ.CB
      IF( LSAME ) &
         RETURN
!
!     Now test for equivalence if both characters are alphabetic.
!
      ZCODE = ICHAR( 'Z' )
!
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDIC machine.
!
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
!
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
!
!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.
!
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
!
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
!
!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!        upper case 'Z'.
!
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR. &
             INTA.GE.145 .AND. INTA.LE.153 .OR. &
             INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR. &
             INTB.GE.145 .AND. INTB.LE.153 .OR. &
             INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
!
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
!
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
!
!     RETURN
!
!     End of LSAME
!
      END FUNCTION LSAME

      INTEGER FUNCTION IEEECK( ISPEC, ZERO, ONE )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            ISPEC
      REAL               ONE, ZERO
!     ..
!
!  Purpose
!  =======
!
!  IEEECK is called from the ILAENV to verify that Infinity and
!  possibly NaN arithmetic is safe (i.e. will not trap).
!
!  Arguments
!  =========
!
!  ISPEC   (input) INTEGER
!          Specifies whether to test just for inifinity arithmetic
!          or whether to test for infinity and NaN arithmetic.
!          = 0: Verify infinity arithmetic only.
!          = 1: Verify infinity and NaN arithmetic.
!
!  ZERO    (input) REAL
!          Must contain the value 0.0
!          This is passed to prevent the compiler from optimizing
!          away this code.
!
!  ONE     (input) REAL
!          Must contain the value 1.0
!          This is passed to prevent the compiler from optimizing
!          away this code.
!
!  RETURN VALUE:  INTEGER
!          = 0:  Arithmetic failed to produce the correct answers
!          = 1:  Arithmetic produced the correct answers
!
!     .. Local Scalars ..
      REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF, &
                         NEGZRO, NEWZRO, POSINF
!     ..
!     .. Executable Statements ..
      IEEECK = 1
!
      POSINF = ONE / ZERO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEGINF = -ONE / ZERO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEGZRO = ONE / ( NEGINF+ONE )
      IF( NEGZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEGINF = ONE / NEGZRO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEWZRO = NEGZRO + ZERO
      IF( NEWZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      POSINF = ONE / NEWZRO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEGINF = NEGINF*POSINF
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      POSINF = POSINF*POSINF
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
!
!
!
!
!     Return if we were only asked to check infinity arithmetic
!
      IF( ISPEC.EQ.0 ) &
         RETURN
!
      NAN1 = POSINF + NEGINF
!
      NAN2 = POSINF / NEGINF
!
      NAN3 = POSINF / POSINF
!
      NAN4 = POSINF*ZERO
!
      NAN5 = NEGINF*NEGZRO
!
      NAN6 = NAN5*0.0
!
      IF( NAN1.EQ.NAN1 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN2.EQ.NAN2 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN3.EQ.NAN3 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN4.EQ.NAN4 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN5.EQ.NAN5 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN6.EQ.NAN6 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      RETURN
      END FUNCTION IEEECK


      INTEGER FUNCTION IPARMQ( ISPEC, ILO, IHI )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            IHI, ILO, ISPEC
!
!  Purpose
!  =======
!
!       This program sets problem and machine dependent parameters
!       useful for xHSEQR and its subroutines. It is called whenever
!       ILAENV is called with 12 <= ISPEC <= 16
!
!  Arguments
!  =========
!
!       ISPEC  (input) integer scalar
!              ISPEC specifies which tunable parameter IPARMQ should
!              return.
!
!              ISPEC=12: (INMIN)  Matrices of order nmin or less
!                        are sent directly to xLAHQR, the implicit
!                        double shift QR algorithm.  NMIN must be
!                        at least 11.
!
!              ISPEC=13: (INWIN)  Size of the deflation window.
!                        This is best set greater than or equal to
!                        the number of simultaneous shifts NS.
!                        Larger matrices benefit from larger deflation
!                        windows.
!
!              ISPEC=14: (INIBL) Determines when to stop nibbling and
!                        invest in an (expensive) multi-shift QR sweep.
!                        If the aggressive early deflation subroutine
!                        finds LD converged eigenvalues from an order
!                        NW deflation window and LD.GT.(NW*NIBBLE)/100,
!                        then the next QR sweep is skipped and early
!                        deflation is applied immediately to the
!                        remaining active diagonal block.  Setting
!                        IPARMQ(ISPEC=14) = 0 causes TTQRE to skip a
!                        multi-shift QR sweep whenever early deflation
!                        finds a converged eigenvalue.  Setting
!                        IPARMQ(ISPEC=14) greater than or equal to 100
!                        prevents TTQRE from skipping a multi-shift
!                        QR sweep.
!
!              ISPEC=15: (NSHFTS) The number of simultaneous shifts in
!                        a multi-shift QR iteration.
!
!              ISPEC=16: (IACC22) IPARMQ is set to 0, 1 or 2 with the
!                        following meanings.
!                        0:  During the multi-shift QR sweep,
!                            xLAQR5 does not accumulate reflections and
!                            does not use matrix-matrix multiply to
!                            update the far-from-diagonal matrix
!                            entries.
!                        1:  During the multi-shift QR sweep,
!                            xLAQR5 and/or xLAQRaccumulates reflections
!                            matrix-matrix multiply to update the
!                            far-from-diagonal matrix entries.
!                        2:  During the multi-shift QR sweep.
!                            xLAQR5 accumulates reflections and takes
!                            advantage of 2-by-2 block structure during
!                            matrix-matrix multiplies.
!                        (If xTRMM is slower than xGEMM, then
!                        IPARMQ(ISPEC=16)=1 may be more efficient than
!                        IPARMQ(ISPEC=16)=2 despite the greater level of
!                        arithmetic work implied by the latter choice.)
!
!       ILO     (input) INTEGER
!       IHI     (input) INTEGER
!               It is assumed that H is already upper triangular
!               in rows and columns 1:ILO-1 and IHI+1:N.
!
!  Further Details
!  ===============
!
!       Little is known about how best to choose these parameters.
!       It is possible to use different values of the parameters
!       for each of CHSEQR, DHSEQR, SHSEQR and ZHSEQR.
!
!       It is probably best to choose different parameters for
!       different matrices and different parameters at different
!       times during the iteration, but this has not been
!       implemented --- yet.
!
!
!       The best choices of most of the parameters depend
!       in an ill-understood way on the relative execution
!       rate of xLAQR3 and xLAQR5 and on the nature of each
!       particular eigenvalue problem.  Experiment may be the
!       only practical way to determine which choices are most
!       effective.
!
!       Following is a list of default values supplied by IPARMQ.
!       These defaults may be adjusted in order to attain better
!       performance in any particular computational environment.
!
!       IPARMQ(ISPEC=12) The xLAHQR vs xLAQR0 crossover point.
!                        Default: 75. (Must be at least 11.)
!
!       IPARMQ(ISPEC=13) Recommended deflation window size.
!                        This depends on ILO, IHI and NS, the
!                        number of simultaneous shifts returned
!                        by IPARMQ(ISPEC=15).  The default for
!                        (IHI-ILO+1).LE.500 is NS.  The default
!                        for (IHI-ILO+1).GT.500 is 3*NS/2.
!
!       IPARMQ(ISPEC=14) Nibble crossover point.  Default: 14.
!
!       IPARMQ(ISPEC=15) Number of simultaneous shifts, NS.
!                        a multi-shift QR iteration.
!
!                        If IHI-ILO+1 is ...
!
!                        greater than      ...but less    ... the
!                        or equal to ...      than        default is
!
!                                0               30       NS =   2+
!                               30               60       NS =   4+
!                               60              150       NS =  10
!                              150              590       NS =  **
!                              590             3000       NS =  64
!                             3000             6000       NS = 128
!                             6000             infinity   NS = 256
!
!                    (+)  By default matrices of this order are
!                         passed to the implicit double shift routine
!                         xLAHQR.  See IPARMQ(ISPEC=12) above.   These
!                         values of NS are used only in case of a rare
!                         xLAHQR failure.
!
!                    (**) The asterisks (**) indicate an ad-hoc
!                         function increasing from 10 to 64.
!
!       IPARMQ(ISPEC=16) Select structured matrix multiply.
!                        (See ISPEC=16 above for details.)
!                        Default: 3.
!
!     ================================================================
!     .. Parameters ..
      INTEGER            INMIN, INWIN, INIBL, ISHFTS, IACC22
      PARAMETER          ( INMIN = 12, INWIN = 13, INIBL = 14, &
                         ISHFTS = 15, IACC22 = 16 )
      INTEGER            NMIN, K22MIN, KACMIN, NIBBLE, KNWSWP
      PARAMETER          ( NMIN = 75, K22MIN = 14, KACMIN = 14, &
                         NIBBLE = 14, KNWSWP = 500 )
      REAL               TWO
      PARAMETER          ( TWO = 2.0 )
!     ..
!     .. Local Scalars ..
      INTEGER            NH, NS
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          LOG, MAX, MOD, NINT, REAL
!     ..
!     .. Executable Statements ..
      IF( ( ISPEC.EQ.ISHFTS ) .OR. ( ISPEC.EQ.INWIN ) .OR. &
          ( ISPEC.EQ.IACC22 ) ) THEN
!
!        ==== Set the number simultaneous shifts ====
!
         NH = IHI - ILO + 1
         NS = 2
         IF( NH.GE.30 ) &
            NS = 4
         IF( NH.GE.60 ) &
            NS = 10
         IF( NH.GE.150 ) &
            NS = MAX( 10, NH / NINT( LOG( REAL( NH ) ) / LOG( TWO ) ) )
         IF( NH.GE.590 ) &
            NS = 64
         IF( NH.GE.3000 ) &
            NS = 128
         IF( NH.GE.6000 ) &
            NS = 256
         NS = MAX( 2, NS-MOD( NS, 2 ) )
      END IF
!
      IF( ISPEC.EQ.INMIN ) THEN
!
!
!        ===== Matrices of order smaller than NMIN get sent
!        .     to xLAHQR, the classic double shift algorithm.
!        .     This must be at least 11. ====
!
         IPARMQ = NMIN
!
      ELSE IF( ISPEC.EQ.INIBL ) THEN
!
!        ==== INIBL: skip a multi-shift qr iteration and
!        .    whenever aggressive early deflation finds
!        .    at least (NIBBLE*(window size)/100) deflations. ====
!
         IPARMQ = NIBBLE
!
      ELSE IF( ISPEC.EQ.ISHFTS ) THEN
!
!        ==== NSHFTS: The number of simultaneous shifts =====
!
         IPARMQ = NS
!
      ELSE IF( ISPEC.EQ.INWIN ) THEN
!
!        ==== NW: deflation window size.  ====
!
         IF( NH.LE.KNWSWP ) THEN
            IPARMQ = NS
         ELSE
            IPARMQ = 3*NS / 2
         END IF
!
      ELSE IF( ISPEC.EQ.IACC22 ) THEN
!
!        ==== IACC22: Whether to accumulate reflections
!        .     before updating the far-from-diagonal elements
!        .     and whether to use 2-by-2 block structure while
!        .     doing it.  A small amount of work could be saved
!        .     by making this choice dependent also upon the
!        .     NH=IHI-ILO+1.
!
         IPARMQ = 0
         IF( NS.GE.KACMIN ) &
            IPARMQ = 1
         IF( NS.GE.K22MIN ) &
            IPARMQ = 2
!
      ELSE
!        ===== invalid value of ispec =====
         IPARMQ = -1
!
      END IF
!
!     ==== End of IPARMQ ====
!
      END FUNCTION IPARMQ


      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
!     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
!     ..
!
!  Purpose
!  =======
!
!     constant times a vector plus a vector.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!
!     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MOD
!     ..
      IF (N.LE.0) RETURN
      IF (DA.EQ.0.0d0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          DY(IY) = DY(IY) + DA*DX(IX)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 M = MOD(N,4)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF (N.LT.4) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
          DY(I) = DY(I) + DA*DX(I)
          DY(I+1) = DY(I+1) + DA*DX(I+1)
          DY(I+2) = DY(I+2) + DA*DX(I+2)
          DY(I+3) = DY(I+3) + DA*DX(I+3)
   50 CONTINUE
      RETURN
      END SUBROUTINE DAXPY


      SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
!     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
!     ..
!
!  Purpose
!  =======
!
!     copies a vector, x, to a vector, y.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!
!     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MOD
!     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          DY(IY) = DX(IX)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 M = MOD(N,7)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DY(I) = DX(I)
   30 CONTINUE
      IF (N.LT.7) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
          DY(I) = DX(I)
          DY(I+1) = DX(I+1)
          DY(I+2) = DX(I+2)
          DY(I+3) = DX(I+3)
          DY(I+4) = DX(I+4)
          DY(I+5) = DX(I+5)
          DY(I+6) = DX(I+6)
   50 CONTINUE
      RETURN
      END SUBROUTINE DCOPY


      SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
!     ..
!
!  Purpose
!  =======
!
!  DGEMM  performs one of the matrix-matrix operations
!
!     C := alpha*op( A )*op( B ) + beta*C,
!
!  where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X',
!
!  alpha and beta are scalars, and A, B and C are matrices, with op( A )
!  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!
!  Arguments
!  ==========
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!
!              TRANSA = 'T' or 't',  op( A ) = A'.
!
!              TRANSA = 'C' or 'c',  op( A ) = A'.
!
!           Unchanged on exit.
!
!  TRANSB - CHARACTER*1.
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!
!              TRANSB = 'T' or 't',  op( B ) = B'.
!
!              TRANSB = 'C' or 'c',  op( B ) = B'.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
!           Unchanged on exit.
!
!  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
!           Unchanged on exit.
!
!  BETA   - DOUBLE PRECISION.
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!           Unchanged on exit.
!
!  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).
!
!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
!      LOGICAL LSAME
!      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,L,NCOLA,NROWA,NROWB
      LOGICAL NOTA,NOTB
!     ..
!     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
!     and  columns of  A  and the  number of  rows  of  B  respectively.
!
      NOTA = LSAME(TRANSA,'N')
      NOTB = LSAME(TRANSB,'N')
      IF (NOTA) THEN
          NROWA = M
          NCOLA = K
      ELSE
          NROWA = K
          NCOLA = M
      END IF
      IF (NOTB) THEN
          NROWB = K
      ELSE
          NROWB = N
      END IF
!
!     Test the input parameters.
!
      INFO = 0
      IF ((.NOT.NOTA) .AND. (.NOT.LSAME(TRANSA,'C')) .AND. &
          (.NOT.LSAME(TRANSA,'T'))) THEN
          INFO = 1
      ELSE IF ((.NOT.NOTB) .AND. (.NOT.LSAME(TRANSB,'C')) .AND. &
               (.NOT.LSAME(TRANSB,'T'))) THEN
          INFO = 2
      ELSE IF (M.LT.0) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (K.LT.0) THEN
          INFO = 5
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 8
      ELSE IF (LDB.LT.MAX(1,NROWB)) THEN
          INFO = 10
      ELSE IF (LDC.LT.MAX(1,M)) THEN
          INFO = 13
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DGEMM ',INFO)
          RETURN
      END IF
!
!     Quick return if possible.
!
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. &
          (((ALPHA.EQ.ZERO).OR. (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
!
!     And if  alpha.eq.zero.
!
      IF (ALPHA.EQ.ZERO) THEN
          IF (BETA.EQ.ZERO) THEN
              DO 20 J = 1,N
                  DO 10 I = 1,M
                      C(I,J) = ZERO
   10             CONTINUE
   20         CONTINUE
          ELSE
              DO 40 J = 1,N
                  DO 30 I = 1,M
                      C(I,J) = BETA*C(I,J)
   30             CONTINUE
   40         CONTINUE
          END IF
          RETURN
      END IF
!
!     Start the operations.
!
      IF (NOTB) THEN
          IF (NOTA) THEN
!
!           Form  C := alpha*A*B + beta*C.
!
              DO 90 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 50 I = 1,M
                          C(I,J) = ZERO
   50                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 60 I = 1,M
                          C(I,J) = BETA*C(I,J)
   60                 CONTINUE
                  END IF
                  DO 80 L = 1,K
                      IF (B(L,J).NE.ZERO) THEN
                          TEMP = ALPHA*B(L,J)
                          DO 70 I = 1,M
                              C(I,J) = C(I,J) + TEMP*A(I,L)
   70                     CONTINUE
                      END IF
   80             CONTINUE
   90         CONTINUE
          ELSE
!
!           Form  C := alpha*A'*B + beta*C
!
              DO 120 J = 1,N
                  DO 110 I = 1,M
                      TEMP = ZERO
                      DO 100 L = 1,K
                          TEMP = TEMP + A(L,I)*B(L,J)
  100                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  110             CONTINUE
  120         CONTINUE
          END IF
      ELSE
          IF (NOTA) THEN
!
!           Form  C := alpha*A*B' + beta*C
!
              DO 170 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 130 I = 1,M
                          C(I,J) = ZERO
  130                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 140 I = 1,M
                          C(I,J) = BETA*C(I,J)
  140                 CONTINUE
                  END IF
                  DO 160 L = 1,K
                      IF (B(J,L).NE.ZERO) THEN
                          TEMP = ALPHA*B(J,L)
                          DO 150 I = 1,M
                              C(I,J) = C(I,J) + TEMP*A(I,L)
  150                     CONTINUE
                      END IF
  160             CONTINUE
  170         CONTINUE
          ELSE
!
!           Form  C := alpha*A'*B' + beta*C
!
              DO 200 J = 1,N
                  DO 190 I = 1,M
                      TEMP = ZERO
                      DO 180 L = 1,K
                          TEMP = TEMP + A(L,I)*B(J,L)
  180                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  190             CONTINUE
  200         CONTINUE
          END IF
      END IF
!
      RETURN
!
!     End of DGEMM .
!
      END SUBROUTINE DGEMM


      SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER INCX,INCY,LDA,M,N
      CHARACTER TRANS
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*),Y(*)
!     ..
!
!  Purpose
!  =======
!
!  DGEMV  performs one of the matrix-vector operations
!
!     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
!
!  where alpha and beta are scalars, x and y are vectors and A is an
!  m by n matrix.
!
!  Arguments
!  ==========
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!
!              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
!
!              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - DOUBLE PRECISION.
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - DOUBLE PRECISION array of DIMENSION at least
!           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!           Before entry with BETA non-zero, the incremented array Y
!           must contain the vector y. On exit, Y is overwritten by the
!           updated vector y.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
!     ..
!     .. External Functions ..
!      LOGICAL LSAME
!      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!
!     Test the input parameters.
!
      INFO = 0
      IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND. &
          .NOT.LSAME(TRANS,'C')) THEN
          INFO = 1
      ELSE IF (M.LT.0) THEN
          INFO = 2
      ELSE IF (N.LT.0) THEN
          INFO = 3
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 6
      ELSE IF (INCX.EQ.0) THEN
          INFO = 8
      ELSE IF (INCY.EQ.0) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DGEMV ',INFO)
          RETURN
      END IF
!
!     Quick return if possible.
!
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. &
          ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
!
!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.
!
      IF (LSAME(TRANS,'N')) THEN
          LENX = N
          LENY = M
      ELSE
          LENX = M
          LENY = N
      END IF
      IF (INCX.GT.0) THEN
          KX = 1
      ELSE
          KX = 1 - (LENX-1)*INCX
      END IF
      IF (INCY.GT.0) THEN
          KY = 1
      ELSE
          KY = 1 - (LENY-1)*INCY
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
!     First form  y := beta*y.
!
      IF (BETA.NE.ONE) THEN
          IF (INCY.EQ.1) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 10 I = 1,LENY
                      Y(I) = ZERO
   10             CONTINUE
              ELSE
                  DO 20 I = 1,LENY
                      Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
          ELSE
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO 30 I = 1,LENY
                      Y(IY) = ZERO
                      IY = IY + INCY
   30             CONTINUE
              ELSE
                  DO 40 I = 1,LENY
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
   40             CONTINUE
              END IF
          END IF
      END IF
      IF (ALPHA.EQ.ZERO) RETURN
      IF (LSAME(TRANS,'N')) THEN
!
!        Form  y := alpha*A*x + y.
!
          JX = KX
          IF (INCY.EQ.1) THEN
              DO 60 J = 1,N
                  IF (X(JX).NE.ZERO) THEN
                      TEMP = ALPHA*X(JX)
                      DO 50 I = 1,M
                          Y(I) = Y(I) + TEMP*A(I,J)
   50                 CONTINUE
                  END IF
                  JX = JX + INCX
   60         CONTINUE
          ELSE
              DO 80 J = 1,N
                  IF (X(JX).NE.ZERO) THEN
                      TEMP = ALPHA*X(JX)
                      IY = KY
                      DO 70 I = 1,M
                          Y(IY) = Y(IY) + TEMP*A(I,J)
                          IY = IY + INCY
   70                 CONTINUE
                  END IF
                  JX = JX + INCX
   80         CONTINUE
          END IF
      ELSE
!
!        Form  y := alpha*A'*x + y.
!
          JY = KY
          IF (INCX.EQ.1) THEN
              DO 100 J = 1,N
                  TEMP = ZERO
                  DO 90 I = 1,M
                      TEMP = TEMP + A(I,J)*X(I)
   90             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  100         CONTINUE
          ELSE
              DO 120 J = 1,N
                  TEMP = ZERO
                  IX = KX
                  DO 110 I = 1,M
                      TEMP = TEMP + A(I,J)*X(IX)
                      IX = IX + INCX
  110             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  120         CONTINUE
          END IF
      END IF
!
      RETURN
!
!     End of DGEMV .
!
      END SUBROUTINE DGEMV


      SUBROUTINE DGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
!     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER INCX,INCY,LDA,M,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*),Y(*)
!     ..
!
!  Purpose
!  =======
!
!  DGER   performs the rank 1 operation
!
!     A := alpha*x*y' + A,
!
!  where alpha is a scalar, x is an m element vector, y is an n element
!  vector and A is an m by n matrix.
!
!  Arguments
!  ==========
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( m - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Y      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,J,JY,KX
!     ..
!     .. External Subroutines ..
!      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!
!     Test the input parameters.
!
      INFO = 0
      IF (M.LT.0) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 5
      ELSE IF (INCY.EQ.0) THEN
          INFO = 7
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 9
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DGER  ',INFO)
          RETURN
      END IF
!
!     Quick return if possible.
!
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      IF (INCY.GT.0) THEN
          JY = 1
      ELSE
          JY = 1 - (N-1)*INCY
      END IF
      IF (INCX.EQ.1) THEN
          DO 20 J = 1,N
              IF (Y(JY).NE.ZERO) THEN
                  TEMP = ALPHA*Y(JY)
                  DO 10 I = 1,M
                      A(I,J) = A(I,J) + X(I)*TEMP
   10             CONTINUE
              END IF
              JY = JY + INCY
   20     CONTINUE
      ELSE
          IF (INCX.GT.0) THEN
              KX = 1
          ELSE
              KX = 1 - (M-1)*INCX
          END IF
          DO 40 J = 1,N
              IF (Y(JY).NE.ZERO) THEN
                  TEMP = ALPHA*Y(JY)
                  IX = KX
                  DO 30 I = 1,M
                      A(I,J) = A(I,J) + X(IX)*TEMP
                      IX = IX + INCX
   30             CONTINUE
              END IF
              JY = JY + INCY
   40     CONTINUE
      END IF
!
      RETURN
!
!     End of DGER  .
!
      END SUBROUTINE DGER


      SUBROUTINE DSCAL(N,DA,DX,INCX)
!     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
!     ..
!
!  Purpose
!  =======
!*
!     scales a vector by a constant.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!
!     .. Local Scalars ..
      INTEGER I,M,MP1,NINCX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MOD
!     ..
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) GO TO 20
!
!        code for increment not equal to 1
!
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
          DX(I) = DA*DX(I)
   10 CONTINUE
      RETURN
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
   20 M = MOD(N,5)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DX(I) = DA*DX(I)
   30 CONTINUE
      IF (N.LT.5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
          DX(I) = DA*DX(I)
          DX(I+1) = DA*DX(I+1)
          DX(I+2) = DA*DX(I+2)
          DX(I+3) = DA*DX(I+3)
          DX(I+4) = DA*DX(I+4)
   50 CONTINUE
      RETURN
      END SUBROUTINE DSCAL


      SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)
!     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
!     ..
!
!  Purpose
!  =======
!
!     interchanges two vectors.
!     uses unrolled loops for increments equal one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!
!     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY,M,MP1
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MOD
!     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
!
!       code for unequal increments or equal increments not equal
!         to 1
!
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          DTEMP = DX(IX)
          DX(IX) = DY(IY)
          DY(IY) = DTEMP
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
!
!       code for both increments equal to 1
!
!
!       clean-up loop
!
   20 M = MOD(N,3)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DTEMP = DX(I)
          DX(I) = DY(I)
          DY(I) = DTEMP
   30 CONTINUE
      IF (N.LT.3) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
          DTEMP = DX(I)
          DX(I) = DY(I)
          DY(I) = DTEMP
          DTEMP = DX(I+1)
          DX(I+1) = DY(I+1)
          DY(I+1) = DTEMP
          DTEMP = DX(I+2)
          DX(I+2) = DY(I+2)
          DY(I+2) = DTEMP
   50 CONTINUE
      RETURN
      END SUBROUTINE DSWAP


      SUBROUTINE DTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
!     .. Scalar Arguments ..
      INTEGER INCX,K,LDA,N
      CHARACTER DIAG,TRANS,UPLO
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*)
!     ..
!
!  Purpose
!  =======
!
!  DTBSV  solves one of the systems of equations
!
!     A*x = b,   or   A'*x = b,
!
!  where b and x are n element vectors and A is an n by n unit, or
!  non-unit, upper or lower triangular band matrix, with ( k + 1 )
!  diagonals.
!
!  No test for singularity or near-singularity is included in this
!  routine. Such tests must be performed before calling this routine.
!
!  Arguments
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the equations to be solved as
!           follows:
!
!              TRANS = 'N' or 'n'   A*x = b.
!
!              TRANS = 'T' or 't'   A'*x = b.
!
!              TRANS = 'C' or 'c'   A'*x = b.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry with UPLO = 'U' or 'u', K specifies the number of
!           super-diagonals of the matrix A.
!           On entry with UPLO = 'L' or 'l', K specifies the number of
!           sub-diagonals of the matrix A.
!           K must satisfy  0 .le. K.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
!           by n part of the array A must contain the upper triangular
!           band part of the matrix of coefficients, supplied column by
!           column, with the leading diagonal of the matrix in row
!           ( k + 1 ) of the array, the first super-diagonal starting at
!           position 2 in row k, and so on. The top left k by k triangle
!           of the array A is not referenced.
!           The following program segment will transfer an upper
!           triangular band matrix from conventional full matrix storage
!           to band storage:
!
!                 DO 20, J = 1, N
!                    M = K + 1 - J
!                    DO 10, I = MAX( 1, J - K ), J
!                       A( M + I, J ) = matrix( I, J )
!              10    CONTINUE
!              20 CONTINUE
!
!           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
!           by n part of the array A must contain the lower triangular
!           band part of the matrix of coefficients, supplied column by
!           column, with the leading diagonal of the matrix in row 1 of
!           the array, the first sub-diagonal starting at position 1 in
!           row 2, and so on. The bottom right k by k triangle of the
!           array A is not referenced.
!           The following program segment will transfer a lower
!           triangular band matrix from conventional full matrix storage
!           to band storage:
!
!                 DO 20, J = 1, N
!                    M = 1 - J
!                    DO 10, I = J, MIN( N, J + K )
!                       A( M + I, J ) = matrix( I, J )
!              10    CONTINUE
!              20 CONTINUE
!
!           Note that when DIAG = 'U' or 'u' the elements of the array A
!           corresponding to the diagonal elements of the matrix are not
!           referenced, but are assumed to be unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           ( k + 1 ).
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element right-hand side vector b. On exit, X is overwritten
!           with the solution vector x.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,J,JX,KPLUS1,KX,L
      LOGICAL NOUNIT
!     ..
!     .. External Functions ..
!      LOGICAL LSAME
!      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX,MIN
!     ..
!
!     Test the input parameters.
!
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND. &
               .NOT.LSAME(TRANS,'C')) THEN
          INFO = 2
      ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (K.LT.0) THEN
          INFO = 5
      ELSE IF (LDA.LT. (K+1)) THEN
          INFO = 7
      ELSE IF (INCX.EQ.0) THEN
          INFO = 9
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DTBSV ',INFO)
          RETURN
      END IF
!
!     Quick return if possible.
!
      IF (N.EQ.0) RETURN
!
      NOUNIT = LSAME(DIAG,'N')
!
!     Set up the start point in X if the increment is not unity. This
!     will be  ( N - 1 )*INCX  too small for descending loops.
!
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed by sequentially with one pass through A.
!
      IF (LSAME(TRANS,'N')) THEN
!
!        Form  x := inv( A )*x.
!
          IF (LSAME(UPLO,'U')) THEN
              KPLUS1 = K + 1
              IF (INCX.EQ.1) THEN
                  DO 20 J = N,1,-1
                      IF (X(J).NE.ZERO) THEN
                          L = KPLUS1 - J
                          IF (NOUNIT) X(J) = X(J)/A(KPLUS1,J)
                          TEMP = X(J)
                          DO 10 I = J - 1,MAX(1,J-K),-1
                              X(I) = X(I) - TEMP*A(L+I,J)
   10                     CONTINUE
                      END IF
   20             CONTINUE
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 40 J = N,1,-1
                      KX = KX - INCX
                      IF (X(JX).NE.ZERO) THEN
                          IX = KX
                          L = KPLUS1 - J
                          IF (NOUNIT) X(JX) = X(JX)/A(KPLUS1,J)
                          TEMP = X(JX)
                          DO 30 I = J - 1,MAX(1,J-K),-1
                              X(IX) = X(IX) - TEMP*A(L+I,J)
                              IX = IX - INCX
   30                     CONTINUE
                      END IF
                      JX = JX - INCX
   40             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 60 J = 1,N
                      IF (X(J).NE.ZERO) THEN
                          L = 1 - J
                          IF (NOUNIT) X(J) = X(J)/A(1,J)
                          TEMP = X(J)
                          DO 50 I = J + 1,MIN(N,J+K)
                              X(I) = X(I) - TEMP*A(L+I,J)
   50                     CONTINUE
                      END IF
   60             CONTINUE
              ELSE
                  JX = KX
                  DO 80 J = 1,N
                      KX = KX + INCX
                      IF (X(JX).NE.ZERO) THEN
                          IX = KX
                          L = 1 - J
                          IF (NOUNIT) X(JX) = X(JX)/A(1,J)
                          TEMP = X(JX)
                          DO 70 I = J + 1,MIN(N,J+K)
                              X(IX) = X(IX) - TEMP*A(L+I,J)
                              IX = IX + INCX
   70                     CONTINUE
                      END IF
                      JX = JX + INCX
   80             CONTINUE
              END IF
          END IF
      ELSE
!
!        Form  x := inv( A')*x.
!
          IF (LSAME(UPLO,'U')) THEN
              KPLUS1 = K + 1
              IF (INCX.EQ.1) THEN
                  DO 100 J = 1,N
                      TEMP = X(J)
                      L = KPLUS1 - J
                      DO 90 I = MAX(1,J-K),J - 1
                          TEMP = TEMP - A(L+I,J)*X(I)
   90                 CONTINUE
                      IF (NOUNIT) TEMP = TEMP/A(KPLUS1,J)
                      X(J) = TEMP
  100             CONTINUE
              ELSE
                  JX = KX
                  DO 120 J = 1,N
                      TEMP = X(JX)
                      IX = KX
                      L = KPLUS1 - J
                      DO 110 I = MAX(1,J-K),J - 1
                          TEMP = TEMP - A(L+I,J)*X(IX)
                          IX = IX + INCX
  110                 CONTINUE
                      IF (NOUNIT) TEMP = TEMP/A(KPLUS1,J)
                      X(JX) = TEMP
                      JX = JX + INCX
                      IF (J.GT.K) KX = KX + INCX
  120             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 140 J = N,1,-1
                      TEMP = X(J)
                      L = 1 - J
                      DO 130 I = MIN(N,J+K),J + 1,-1
                          TEMP = TEMP - A(L+I,J)*X(I)
  130                 CONTINUE
                      IF (NOUNIT) TEMP = TEMP/A(1,J)
                      X(J) = TEMP
  140             CONTINUE
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 160 J = N,1,-1
                      TEMP = X(JX)
                      IX = KX
                      L = 1 - J
                      DO 150 I = MIN(N,J+K),J + 1,-1
                          TEMP = TEMP - A(L+I,J)*X(IX)
                          IX = IX - INCX
  150                 CONTINUE
                      IF (NOUNIT) TEMP = TEMP/A(1,J)
                      X(JX) = TEMP
                      JX = JX - INCX
                      IF ((N-J).GE.K) KX = KX - INCX
  160             CONTINUE
              END IF
          END IF
      END IF
!
      RETURN
!
!     End of DTBSV .
!
      END SUBROUTINE DTBSV


      SUBROUTINE DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
!     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*)
!     ..
!
!  Purpose
!  =======
!
!  DTRSM  solves one of the matrix equations
!
!     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
!
!  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
!  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!
!     op( A ) = A   or   op( A ) = A'.
!
!  The matrix X is overwritten on B.
!
!  Arguments
!  ==========
!
!  SIDE   - CHARACTER*1.
!           On entry, SIDE specifies whether op( A ) appears on the left
!           or right of X as follows:
!
!              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
!
!              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
!
!           Unchanged on exit.
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix A is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n'   op( A ) = A.
!
!              TRANSA = 'T' or 't'   op( A ) = A'.
!
!              TRANSA = 'C' or 'c'   op( A ) = A'.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit triangular
!           as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of B. M must be at
!           least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of B.  N must be
!           at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
!           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
!           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!           upper triangular part of the array  A must contain the upper
!           triangular matrix  and the strictly lower triangular part of
!           A is not referenced.
!           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!           lower triangular part of the array  A must contain the lower
!           triangular matrix  and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!           A  are not referenced either,  but are assumed to be  unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!           Unchanged on exit.
!
!  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain  the  right-hand  side  matrix  B,  and  on exit  is
!           overwritten by the solution matrix  X.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
!      LOGICAL LSAME
!      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOUNIT,UPPER
!     ..
!     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!
!     Test the input parameters.
!
      LSIDE = LSAME(SIDE,'L')
      IF (LSIDE) THEN
          NROWA = M
      ELSE
          NROWA = N
      END IF
      NOUNIT = LSAME(DIAG,'N')
      UPPER = LSAME(UPLO,'U')
!
      INFO = 0
      IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
          INFO = 1
      ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 2
      ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND. &
               (.NOT.LSAME(TRANSA,'T')) .AND. &
               (.NOT.LSAME(TRANSA,'C'))) THEN
          INFO = 3
      ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) THEN
          INFO = 4
      ELSE IF (M.LT.0) THEN
          INFO = 5
      ELSE IF (N.LT.0) THEN
          INFO = 6
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 9
      ELSE IF (LDB.LT.MAX(1,M)) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DTRSM ',INFO)
          RETURN
      END IF
!
!     Quick return if possible.
!
      IF (N.EQ.0) RETURN
!
!     And when  alpha.eq.zero.
!
      IF (ALPHA.EQ.ZERO) THEN
          DO 20 J = 1,N
              DO 10 I = 1,M
                  B(I,J) = ZERO
   10         CONTINUE
   20     CONTINUE
          RETURN
      END IF
!
!     Start the operations.
!
      IF (LSIDE) THEN
          IF (LSAME(TRANSA,'N')) THEN
!
!           Form  B := alpha*inv( A )*B.
!
              IF (UPPER) THEN
                  DO 60 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 30 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
   30                     CONTINUE
                      END IF
                      DO 50 K = M,1,-1
                          IF (B(K,J).NE.ZERO) THEN
                              IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                              DO 40 I = 1,K - 1
                                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
   40                         CONTINUE
                          END IF
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 100 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 70 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
   70                     CONTINUE
                      END IF
                      DO 90 K = 1,M
                          IF (B(K,J).NE.ZERO) THEN
                              IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                              DO 80 I = K + 1,M
                                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
   80                         CONTINUE
                          END IF
   90                 CONTINUE
  100             CONTINUE
              END IF
          ELSE
!
!           Form  B := alpha*inv( A' )*B.
!
              IF (UPPER) THEN
                  DO 130 J = 1,N
                      DO 120 I = 1,M
                          TEMP = ALPHA*B(I,J)
                          DO 110 K = 1,I - 1
                              TEMP = TEMP - A(K,I)*B(K,J)
  110                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(I,I)
                          B(I,J) = TEMP
  120                 CONTINUE
  130             CONTINUE
              ELSE
                  DO 160 J = 1,N
                      DO 150 I = M,1,-1
                          TEMP = ALPHA*B(I,J)
                          DO 140 K = I + 1,M
                              TEMP = TEMP - A(K,I)*B(K,J)
  140                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(I,I)
                          B(I,J) = TEMP
  150                 CONTINUE
  160             CONTINUE
              END IF
          END IF
      ELSE
          IF (LSAME(TRANSA,'N')) THEN
!
!           Form  B := alpha*B*inv( A ).
!
              IF (UPPER) THEN
                  DO 210 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 170 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
  170                     CONTINUE
                      END IF
                      DO 190 K = 1,J - 1
                          IF (A(K,J).NE.ZERO) THEN
                              DO 180 I = 1,M
                                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
  180                         CONTINUE
                          END IF
  190                 CONTINUE
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(J,J)
                          DO 200 I = 1,M
                              B(I,J) = TEMP*B(I,J)
  200                     CONTINUE
                      END IF
  210             CONTINUE
              ELSE
                  DO 260 J = N,1,-1
                      IF (ALPHA.NE.ONE) THEN
                          DO 220 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
  220                     CONTINUE
                      END IF
                      DO 240 K = J + 1,N
                          IF (A(K,J).NE.ZERO) THEN
                              DO 230 I = 1,M
                                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
  230                         CONTINUE
                          END IF
  240                 CONTINUE
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(J,J)
                          DO 250 I = 1,M
                              B(I,J) = TEMP*B(I,J)
  250                     CONTINUE
                      END IF
  260             CONTINUE
              END IF
          ELSE
!
!           Form  B := alpha*B*inv( A' ).
!
              IF (UPPER) THEN
                  DO 310 K = N,1,-1
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(K,K)
                          DO 270 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  270                     CONTINUE
                      END IF
                      DO 290 J = 1,K - 1
                          IF (A(J,K).NE.ZERO) THEN
                              TEMP = A(J,K)
                              DO 280 I = 1,M
                                  B(I,J) = B(I,J) - TEMP*B(I,K)
  280                         CONTINUE
                          END IF
  290                 CONTINUE
                      IF (ALPHA.NE.ONE) THEN
                          DO 300 I = 1,M
                              B(I,K) = ALPHA*B(I,K)
  300                     CONTINUE
                      END IF
  310             CONTINUE
              ELSE
                  DO 360 K = 1,N
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(K,K)
                          DO 320 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  320                     CONTINUE
                      END IF
                      DO 340 J = K + 1,N
                          IF (A(J,K).NE.ZERO) THEN
                              TEMP = A(J,K)
                              DO 330 I = 1,M
                                  B(I,J) = B(I,J) - TEMP*B(I,K)
  330                         CONTINUE
                          END IF
  340                 CONTINUE
                      IF (ALPHA.NE.ONE) THEN
                          DO 350 I = 1,M
                              B(I,K) = ALPHA*B(I,K)
  350                     CONTINUE
                      END IF
  360             CONTINUE
              END IF
          END IF
      END IF
!
      RETURN
!
!     End of DTRSM .
!
      END SUBROUTINE DTRSM


      SUBROUTINE DTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
!     .. Scalar Arguments ..
      INTEGER INCX,LDA,N
      CHARACTER DIAG,TRANS,UPLO
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*)
!     ..
!
!  Purpose
!  =======
!
!  DTRSV  solves one of the systems of equations
!
!     A*x = b,   or   A'*x = b,
!
!  where b and x are n element vectors and A is an n by n unit, or
!  non-unit, upper or lower triangular matrix.
!
!  No test for singularity or near-singularity is included in this
!  routine. Such tests must be performed before calling this routine.
!
!  Arguments
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the equations to be solved as
!           follows:
!
!              TRANS = 'N' or 'n'   A*x = b.
!
!              TRANS = 'T' or 't'   A'*x = b.
!
!              TRANS = 'C' or 'c'   A'*x = b.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular matrix and the strictly lower triangular part of
!           A is not referenced.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular matrix and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u', the diagonal elements of
!           A are not referenced either, but are assumed to be unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element right-hand side vector b. On exit, X is overwritten
!           with the solution vector x.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,J,JX,KX
      LOGICAL NOUNIT
!     ..
!     .. External Functions ..
!      LOGICAL LSAME
!      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!
!     Test the input parameters.
!
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND. &
               .NOT.LSAME(TRANS,'C')) THEN
          INFO = 2
      ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (LDA.LT.MAX(1,N)) THEN
          INFO = 6
      ELSE IF (INCX.EQ.0) THEN
          INFO = 8
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DTRSV ',INFO)
          RETURN
      END IF
!
!     Quick return if possible.
!
      IF (N.EQ.0) RETURN
!
      NOUNIT = LSAME(DIAG,'N')
!
!     Set up the start point in X if the increment is not unity. This
!     will be  ( N - 1 )*INCX  too small for descending loops.
!
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      IF (LSAME(TRANS,'N')) THEN
!
!        Form  x := inv( A )*x.
!
          IF (LSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 20 J = N,1,-1
                      IF (X(J).NE.ZERO) THEN
                          IF (NOUNIT) X(J) = X(J)/A(J,J)
                          TEMP = X(J)
                          DO 10 I = J - 1,1,-1
                              X(I) = X(I) - TEMP*A(I,J)
   10                     CONTINUE
                      END IF
   20             CONTINUE
              ELSE
                  JX = KX + (N-1)*INCX
                  DO 40 J = N,1,-1
                      IF (X(JX).NE.ZERO) THEN
                          IF (NOUNIT) X(JX) = X(JX)/A(J,J)
                          TEMP = X(JX)
                          IX = JX
                          DO 30 I = J - 1,1,-1
                              IX = IX - INCX
                              X(IX) = X(IX) - TEMP*A(I,J)
   30                     CONTINUE
                      END IF
                      JX = JX - INCX
   40             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 60 J = 1,N
                      IF (X(J).NE.ZERO) THEN
                          IF (NOUNIT) X(J) = X(J)/A(J,J)
                          TEMP = X(J)
                          DO 50 I = J + 1,N
                              X(I) = X(I) - TEMP*A(I,J)
   50                     CONTINUE
                      END IF
   60             CONTINUE
              ELSE
                  JX = KX
                  DO 80 J = 1,N
                      IF (X(JX).NE.ZERO) THEN
                          IF (NOUNIT) X(JX) = X(JX)/A(J,J)
                          TEMP = X(JX)
                          IX = JX
                          DO 70 I = J + 1,N
                              IX = IX + INCX
                              X(IX) = X(IX) - TEMP*A(I,J)
   70                     CONTINUE
                      END IF
                      JX = JX + INCX
   80             CONTINUE
              END IF
          END IF
      ELSE
!
!        Form  x := inv( A' )*x.
!
          IF (LSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 100 J = 1,N
                      TEMP = X(J)
                      DO 90 I = 1,J - 1
                          TEMP = TEMP - A(I,J)*X(I)
   90                 CONTINUE
                      IF (NOUNIT) TEMP = TEMP/A(J,J)
                      X(J) = TEMP
  100             CONTINUE
              ELSE
                  JX = KX
                  DO 120 J = 1,N
                      TEMP = X(JX)
                      IX = KX
                      DO 110 I = 1,J - 1
                          TEMP = TEMP - A(I,J)*X(IX)
                          IX = IX + INCX
  110                 CONTINUE
                      IF (NOUNIT) TEMP = TEMP/A(J,J)
                      X(JX) = TEMP
                      JX = JX + INCX
  120             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 140 J = N,1,-1
                      TEMP = X(J)
                      DO 130 I = N,J + 1,-1
                          TEMP = TEMP - A(I,J)*X(I)
  130                 CONTINUE
                      IF (NOUNIT) TEMP = TEMP/A(J,J)
                      X(J) = TEMP
  140             CONTINUE
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 160 J = N,1,-1
                      TEMP = X(JX)
                      IX = KX
                      DO 150 I = N,J + 1,-1
                          TEMP = TEMP - A(I,J)*X(IX)
                          IX = IX - INCX
  150                 CONTINUE
                      IF (NOUNIT) TEMP = TEMP/A(J,J)
                      X(JX) = TEMP
                      JX = JX - INCX
  160             CONTINUE
              END IF
          END IF
      END IF
!
      RETURN
!
!     End of DTRSV .
!
      END SUBROUTINE DTRSV


      DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX)
!     .. Scalar Arguments ..
      INTEGER INCX,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
!     ..
!
!  Purpose
!  =======
!
!     takes the sum of the absolute values.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!
!     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,M,MP1,NINCX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DABS,MOD
!     ..
      DASUM = 0.0d0
      DTEMP = 0.0d0
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) GO TO 20
!
!        code for increment not equal to 1
!
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
          DTEMP = DTEMP + DABS(DX(I))
   10 CONTINUE
      DASUM = DTEMP
      RETURN
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
   20 M = MOD(N,6)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DTEMP = DTEMP + DABS(DX(I))
   30 CONTINUE
      IF (N.LT.6) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
          DTEMP = DTEMP + DABS(DX(I)) + DABS(DX(I+1)) + DABS(DX(I+2)) + &
                  DABS(DX(I+3)) + DABS(DX(I+4)) + DABS(DX(I+5))
   50 CONTINUE
   60 DASUM = DTEMP
      RETURN
      END FUNCTION DASUM


      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
!     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
!     ..
!
!  Purpose
!  =======
!
!     forms the dot product of two vectors.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!
!     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY,M,MP1
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MOD
!     ..
      DDOT = 0.0d0
      DTEMP = 0.0d0
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          DTEMP = DTEMP + DX(IX)*DY(IY)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      DDOT = DTEMP
      RETURN
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 M = MOD(N,5)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DTEMP = DTEMP + DX(I)*DY(I)
   30 CONTINUE
      IF (N.LT.5) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
          DTEMP = DTEMP + DX(I)*DY(I) + DX(I+1)*DY(I+1) + &
                  DX(I+2)*DY(I+2) + DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
   50 CONTINUE
   60 DDOT = DTEMP
      RETURN
      END FUNCTION DDOT


      INTEGER FUNCTION IDAMAX(N,DX,INCX)
!     .. Scalar Arguments ..
      INTEGER INCX,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
!     ..
!
!  Purpose
!  =======
!
!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!
!     .. Local Scalars ..
      DOUBLE PRECISION DMAX
      INTEGER I,IX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DABS
!     ..
      IDAMAX = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      IDAMAX = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) GO TO 20
!
!        code for increment not equal to 1
!
      IX = 1
      DMAX = DABS(DX(1))
      IX = IX + INCX
      DO 10 I = 2,N
          IF (DABS(DX(IX)).LE.DMAX) GO TO 5
          IDAMAX = I
          DMAX = DABS(DX(IX))
    5     IX = IX + INCX
   10 CONTINUE
      RETURN
!
!        code for increment equal to 1
!
   20 DMAX = DABS(DX(1))
      DO 30 I = 2,N
          IF (DABS(DX(I)).LE.DMAX) GO TO 30
          IDAMAX = I
          DMAX = DABS(DX(I))
   30 CONTINUE
      RETURN
      END FUNCTION IDAMAX


      SUBROUTINE DROT(N,DX,INCX,DY,INCY,C,S)
!     .. Scalar Arguments ..
      DOUBLE PRECISION C,S
      INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
!     ..
!
!  Purpose
!  =======
!
!     applies a plane rotation.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!
!     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY
!     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
!
!       code for unequal increments or equal increments not equal
!         to 1
!
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          DTEMP = C*DX(IX) + S*DY(IY)
          DY(IY) = C*DY(IY) - S*DX(IX)
          DX(IX) = DTEMP
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
!
!       code for both increments equal to 1
!
   20 DO 30 I = 1,N
          DTEMP = C*DX(I) + S*DY(I)
          DY(I) = C*DY(I) - S*DX(I)
          DX(I) = DTEMP
   30 CONTINUE
      RETURN
      END SUBROUTINE DROT

      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          CMACH
!     ..
!
!  Purpose
!  =======
!
!  DLAMCH determines double precision machine parameters.
!
!  Arguments
!  =========
!
!  CMACH   (input) CHARACTER*1
!          Specifies the value to be returned by DLAMCH:
!          = 'E' or 'e',   DLAMCH := eps
!          = 'S' or 's ,   DLAMCH := sfmin
!          = 'B' or 'b',   DLAMCH := base
!          = 'P' or 'p',   DLAMCH := eps*base
!          = 'N' or 'n',   DLAMCH := t
!          = 'R' or 'r',   DLAMCH := rnd
!          = 'M' or 'm',   DLAMCH := emin
!          = 'U' or 'u',   DLAMCH := rmin
!          = 'L' or 'l',   DLAMCH := emax
!          = 'O' or 'o',   DLAMCH := rmax
!
!          where
!
!          eps   = relative machine precision
!          sfmin = safe minimum, such that 1/sfmin does not overflow
!          base  = base of the machine
!          prec  = eps*base
!          t     = number of (base) digits in the mantissa
!          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
!          emin  = minimum exponent before (gradual) underflow
!          rmin  = underflow threshold - base**(emin-1)
!          emax  = largest exponent before overflow
!          rmax  = overflow threshold  - (base**emax)*(1-eps)
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            FIRST, LRND
      INTEGER            BETA, IMAX, IMIN, IT
      DOUBLE PRECISION   BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN, &
                         RND, SFMIN, SMALL, T
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DLAMC2
!     ..
!     .. Save statement ..
      SAVE               FIRST, EPS, SFMIN, BASE, T, RND, EMIN, RMIN, &
                         EMAX, RMAX, PREC
!     ..
!     .. Data statements ..
      DATA               FIRST / .TRUE. /
!     ..
!     .. Executable Statements ..
!
      IF( FIRST ) THEN
         CALL DLAMC2( BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX )
         BASE = BETA
         T = IT
         IF( LRND ) THEN
            RND = ONE
            EPS = ( BASE**( 1-IT ) ) / 2
         ELSE
            RND = ZERO
            EPS = BASE**( 1-IT )
         END IF
         PREC = EPS*BASE
         EMIN = IMIN
         EMAX = IMAX
         SFMIN = RMIN
         SMALL = ONE / RMAX
         IF( SMALL.GE.SFMIN ) THEN
!
!           Use SMALL plus a bit, to avoid the possibility of rounding
!           causing overflow when computing  1/sfmin.
!
            SFMIN = SMALL*( ONE+EPS )
         END IF
      END IF
!
      IF( LSAME( CMACH, 'E' ) ) THEN
         RMACH = EPS
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
         RMACH = SFMIN
      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
         RMACH = BASE
      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
         RMACH = PREC
      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
         RMACH = T
      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
         RMACH = RND
      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
         RMACH = EMIN
      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
         RMACH = RMIN
      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
         RMACH = EMAX
      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
         RMACH = RMAX
      END IF
!
      DLAMCH = RMACH
      FIRST  = .FALSE.
      RETURN
!
!     End of DLAMCH
!
      END FUNCTION DLAMCH
!
!***********************************************************************
!
      SUBROUTINE DLAMC1( BETA, T, RND, IEEE1 )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      LOGICAL            IEEE1, RND
      INTEGER            BETA, T
!     ..
!
!  Purpose
!  =======
!
!  DLAMC1 determines the machine parameters given by BETA, T, RND, and
!  IEEE1.
!
!  Arguments
!  =========
!
!  BETA    (output) INTEGER
!          The base of the machine.
!
!  T       (output) INTEGER
!          The number of ( BETA ) digits in the mantissa.
!
!  RND     (output) LOGICAL
!          Specifies whether proper rounding  ( RND = .TRUE. )  or
!          chopping  ( RND = .FALSE. )  occurs in addition. This may not
!          be a reliable guide to the way in which the machine performs
!          its arithmetic.
!
!  IEEE1   (output) LOGICAL
!          Specifies whether rounding appears to be done in the IEEE
!          'round to nearest' style.
!
!  Further Details
!  ===============
!
!  The routine is based on the routine  ENVRON  by Malcolm and
!  incorporates suggestions by Gentleman and Marovich. See
!
!     Malcolm M. A. (1972) Algorithms to reveal properties of
!        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
!
!     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
!        that reveal properties of floating point arithmetic units.
!        Comms. of the ACM, 17, 276-277.
!
! =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            FIRST, LIEEE1, LRND
      INTEGER            LBETA, LT
      DOUBLE PRECISION   A, B, C, F, ONE, QTR, SAVEC, T1, T2
!     ..
!     .. External Functions ..
!      DOUBLE PRECISION   DLAMC3
!      EXTERNAL           DLAMC3
!     ..
!     .. Save statement ..
      SAVE               FIRST, LIEEE1, LBETA, LRND, LT
!     ..
!     .. Data statements ..
      DATA               FIRST / .TRUE. /
!     ..
!     .. Executable Statements ..
!
      IF( FIRST ) THEN
         ONE = 1
!
!        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
!        IEEE1, T and RND.
!
!        Throughout this routine  we use the function  DLAMC3  to ensure
!        that relevant values are  stored and not held in registers,  or
!        are not affected by optimizers.
!
!        Compute  a = 2.0**m  with the  smallest positive integer m such
!        that
!
!           fl( a + 1.0 ) = a.
!
         A = 1
         C = 1
!
!+       WHILE( C.EQ.ONE )LOOP
   10    CONTINUE
         IF( C.EQ.ONE ) THEN
            A = 2*A
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 10
         END IF
!+       END WHILE
!
!        Now compute  b = 2.0**m  with the smallest positive integer m
!        such that
!
!           fl( a + b ) .gt. a.
!
         B = 1
         C = DLAMC3( A, B )
!
!+       WHILE( C.EQ.A )LOOP
   20    CONTINUE
         IF( C.EQ.A ) THEN
            B = 2*B
            C = DLAMC3( A, B )
            GO TO 20
         END IF
!+       END WHILE
!
!        Now compute the base.  a and c  are neighbouring floating point
!        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
!        their difference is beta. Adding 0.25 to c is to ensure that it
!        is truncated to beta and not ( beta - 1 ).
!
         QTR = ONE / 4
         SAVEC = C
         C = DLAMC3( C, -A )
         LBETA = C + QTR
!
!        Now determine whether rounding or chopping occurs,  by adding a
!        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
!
         B = LBETA
         F = DLAMC3( B / 2, -B / 100 )
         C = DLAMC3( F, A )
         IF( C.EQ.A ) THEN
            LRND = .TRUE.
         ELSE
            LRND = .FALSE.
         END IF
         F = DLAMC3( B / 2, B / 100 )
         C = DLAMC3( F, A )
         IF( ( LRND ) .AND. ( C.EQ.A ) ) &
            LRND = .FALSE.
!
!        Try and decide whether rounding is done in the  IEEE  'round to
!        nearest' style. B/2 is half a unit in the last place of the two
!        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
!        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
!        A, but adding B/2 to SAVEC should change SAVEC.
!
         T1 = DLAMC3( B / 2, A )
         T2 = DLAMC3( B / 2, SAVEC )
         LIEEE1 = ( T1.EQ.A ) .AND. ( T2.GT.SAVEC ) .AND. LRND
!
!        Now find  the  mantissa, t.  It should  be the  integer part of
!        log to the base beta of a,  however it is safer to determine  t
!        by powering.  So we find t as the smallest positive integer for
!        which
!
!           fl( beta**t + 1.0 ) = 1.0.
!
         LT = 0
         A = 1
         C = 1
!
!+       WHILE( C.EQ.ONE )LOOP
   30    CONTINUE
         IF( C.EQ.ONE ) THEN
            LT = LT + 1
            A = A*LBETA
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 30
         END IF
!+       END WHILE
!
      END IF
!
      BETA = LBETA
      T = LT
      RND = LRND
      IEEE1 = LIEEE1
      FIRST = .FALSE.
      RETURN
!
!     End of DLAMC1
!
      END SUBROUTINE DLAMC1
!
!***********************************************************************
!
      SUBROUTINE DLAMC2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      LOGICAL            RND
      INTEGER            BETA, EMAX, EMIN, T
      DOUBLE PRECISION   EPS, RMAX, RMIN
!     ..
!
!  Purpose
!  =======
!
!  DLAMC2 determines the machine parameters specified in its argument
!  list.
!
!  Arguments
!  =========
!
!  BETA    (output) INTEGER
!          The base of the machine.
!
!  T       (output) INTEGER
!          The number of ( BETA ) digits in the mantissa.
!
!  RND     (output) LOGICAL
!          Specifies whether proper rounding  ( RND = .TRUE. )  or
!          chopping  ( RND = .FALSE. )  occurs in addition. This may not
!          be a reliable guide to the way in which the machine performs
!          its arithmetic.
!
!  EPS     (output) DOUBLE PRECISION
!          The smallest positive number such that
!
!             fl( 1.0 - EPS ) .LT. 1.0,
!
!          where fl denotes the computed value.
!
!  EMIN    (output) INTEGER
!          The minimum exponent before (gradual) underflow occurs.
!
!  RMIN    (output) DOUBLE PRECISION
!          The smallest normalized number for the machine, given by
!          BASE**( EMIN - 1 ), where  BASE  is the floating point value
!          of BETA.
!
!  EMAX    (output) INTEGER
!          The maximum exponent before overflow occurs.
!
!  RMAX    (output) DOUBLE PRECISION
!          The largest positive number for the machine, given by
!          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
!          value of BETA.
!
!  Further Details
!  ===============
!
!  The computation of  EPS  is based on a routine PARANOIA by
!  W. Kahan of the University of California at Berkeley.
!
! =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            FIRST, IEEE, IWARN, LIEEE1, LRND
      INTEGER            GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT, &
                         NGNMIN, NGPMIN
      DOUBLE PRECISION   A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE, &
                         SIXTH, SMALL, THIRD, TWO, ZERO
!     ..
!     .. External Functions ..
!      DOUBLE PRECISION   DLAMC3
!      EXTERNAL           DLAMC3
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DLAMC1, DLAMC4, DLAMC5
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
!     ..
!     .. Save statement ..
      SAVE               FIRST, IWARN, LBETA, LEMAX, LEMIN, LEPS, LRMAX, &
                         LRMIN, LT
!     ..
!     .. Data statements ..
      DATA               FIRST / .TRUE. / , IWARN / .FALSE. /
!     ..
!     .. Executable Statements ..
!
      IF( FIRST ) THEN
         ZERO = 0
         ONE = 1
         TWO = 2
!
!        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
!        BETA, T, RND, EPS, EMIN and RMIN.
!
!        Throughout this routine  we use the function  DLAMC3  to ensure
!        that relevant values are stored  and not held in registers,  or
!        are not affected by optimizers.
!
!        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
!
         CALL DLAMC1( LBETA, LT, LRND, LIEEE1 )
!
!        Start to find EPS.
!
         B = LBETA
         A = B**( -LT )
         LEPS = A
!
!        Try some tricks to see whether or not this is the correct  EPS.
!
         B = TWO / 3
         HALF = ONE / 2
         SIXTH = DLAMC3( B, -HALF )
         THIRD = DLAMC3( SIXTH, SIXTH )
         B = DLAMC3( THIRD, -HALF )
         B = DLAMC3( B, SIXTH )
         B = ABS( B )
         IF( B.LT.LEPS ) &
            B = LEPS
!
         LEPS = 1
!
!+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
   10    CONTINUE
         IF( ( LEPS.GT.B ) .AND. ( B.GT.ZERO ) ) THEN
            LEPS = B
            C = DLAMC3( HALF*LEPS, ( TWO**5 )*( LEPS**2 ) )
            C = DLAMC3( HALF, -C )
            B = DLAMC3( HALF, C )
            C = DLAMC3( HALF, -B )
            B = DLAMC3( HALF, C )
            GO TO 10
         END IF
!+       END WHILE
!
         IF( A.LT.LEPS ) &
            LEPS = A
!
!        Computation of EPS complete.
!
!        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
!        Keep dividing  A by BETA until (gradual) underflow occurs. This
!        is detected when we cannot recover the previous A.
!
         RBASE = ONE / LBETA
         SMALL = ONE
         DO 20 I = 1, 3
            SMALL = DLAMC3( SMALL*RBASE, ZERO )
   20    CONTINUE
         A = DLAMC3( ONE, SMALL )
         CALL DLAMC4( NGPMIN, ONE, LBETA )
         CALL DLAMC4( NGNMIN, -ONE, LBETA )
         CALL DLAMC4( GPMIN, A, LBETA )
         CALL DLAMC4( GNMIN, -A, LBETA )
         IEEE = .FALSE.
!
         IF( ( NGPMIN.EQ.NGNMIN ) .AND. ( GPMIN.EQ.GNMIN ) ) THEN
            IF( NGPMIN.EQ.GPMIN ) THEN
               LEMIN = NGPMIN
!            ( Non twos-complement machines, no gradual underflow;
!              e.g.,  VAX )
            ELSE IF( ( GPMIN-NGPMIN ).EQ.3 ) THEN
               LEMIN = NGPMIN - 1 + LT
               IEEE = .TRUE.
!            ( Non twos-complement machines, with gradual underflow;
!              e.g., IEEE standard followers )
            ELSE
               LEMIN = MIN( NGPMIN, GPMIN )
!            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
!
         ELSE IF( ( NGPMIN.EQ.GPMIN ) .AND. ( NGNMIN.EQ.GNMIN ) ) THEN
            IF( ABS( NGPMIN-NGNMIN ).EQ.1 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN )
!            ( Twos-complement machines, no gradual underflow;
!              e.g., CYBER 205 )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
!            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
!
         ELSE IF( ( ABS( NGPMIN-NGNMIN ).EQ.1 ) .AND. &
                  ( GPMIN.EQ.GNMIN ) ) THEN
            IF( ( GPMIN-MIN( NGPMIN, NGNMIN ) ).EQ.3 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN ) - 1 + LT
!            ( Twos-complement machines with gradual underflow;
!              no known machine )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
!            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
!
         ELSE
            LEMIN = MIN( NGPMIN, NGNMIN, GPMIN, GNMIN )
!         ( A guess; no known machine )
            IWARN = .TRUE.
         END IF
         FIRST = .FALSE.
!**
! Comment out this if block if EMIN is ok
         IF( IWARN ) THEN
            FIRST = .TRUE.
            WRITE( 6, FMT = 9999 )LEMIN
         END IF
!**
!
!        Assume IEEE arithmetic if we found denormalised  numbers above,
!        or if arithmetic seems to round in the  IEEE style,  determined
!        in routine DLAMC1. A true IEEE machine should have both  things
!        true; however, faulty machines may have one or the other.
!
         IEEE = IEEE .OR. LIEEE1
!
!        Compute  RMIN by successive division by  BETA. We could compute
!        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
!        this computation.
!
         LRMIN = 1
         DO 30 I = 1, 1 - LEMIN
            LRMIN = DLAMC3( LRMIN*RBASE, ZERO )
   30    CONTINUE
!
!        Finally, call DLAMC5 to compute EMAX and RMAX.
!
         CALL DLAMC5( LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX )
      END IF
!
      BETA = LBETA
      T = LT
      RND = LRND
      EPS = LEPS
      EMIN = LEMIN
      RMIN = LRMIN
      EMAX = LEMAX
      RMAX = LRMAX
!
      RETURN
!
 9999 FORMAT( / / ' WARNING. The value EMIN may be incorrect:-', &
            '  EMIN = ', I8, / &
            ' If, after inspection, the value EMIN looks', &
            ' acceptable please comment out ', &
            / ' the IF block as marked within the code of routine', &
            ' DLAMC2,', / ' otherwise supply EMIN explicitly.', / )
!
!     End of DLAMC2
!
      END SUBROUTINE DLAMC2
!
!***********************************************************************
!
      DOUBLE PRECISION FUNCTION DLAMC3( A, B )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B
!     ..
!
!  Purpose
!  =======
!
!  DLAMC3  is intended to force  A  and  B  to be stored prior to doing
!  the addition of  A  and  B ,  for use in situations where optimizers
!  might hold one of these in a register.
!
!  Arguments
!  =========
!
!  A       (input) DOUBLE PRECISION
!  B       (input) DOUBLE PRECISION
!          The values A and B.
!
! =====================================================================
!
!     .. Executable Statements ..
!
      DLAMC3 = A + B
!
      RETURN
!
!     End of DLAMC3
!
      END FUNCTION DLAMC3
!
!***********************************************************************
!
      SUBROUTINE DLAMC4( EMIN, START, BASE )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            BASE, EMIN
      DOUBLE PRECISION   START
!     ..
!
!  Purpose
!  =======
!
!  DLAMC4 is a service routine for DLAMC2.
!
!  Arguments
!  =========
!
!  EMIN    (output) INTEGER
!          The minimum exponent before (gradual) underflow, computed by
!          setting A = START and dividing by BASE until the previous A
!          can not be recovered.
!
!  START   (input) DOUBLE PRECISION
!          The starting point for determining EMIN.
!
!  BASE    (input) INTEGER
!          The base of the machine.
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
!     ..
!     .. External Functions ..
!      DOUBLE PRECISION   DLAMC3
!      EXTERNAL           DLAMC3
!     ..
!     .. Executable Statements ..
!
      A = START
      ONE = 1
      RBASE = ONE / BASE
      ZERO = 0
      EMIN = 1
      B1 = DLAMC3( A*RBASE, ZERO )
      C1 = A
      C2 = A
      D1 = A
      D2 = A
!+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
!    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
   10 CONTINUE
      IF( ( C1.EQ.A ) .AND. ( C2.EQ.A ) .AND. ( D1.EQ.A ) .AND. &
          ( D2.EQ.A ) ) THEN
         EMIN = EMIN - 1
         A = B1
         B1 = DLAMC3( A / BASE, ZERO )
         C1 = DLAMC3( B1*BASE, ZERO )
         D1 = ZERO
         DO 20 I = 1, BASE
            D1 = D1 + B1
   20    CONTINUE
         B2 = DLAMC3( A*RBASE, ZERO )
         C2 = DLAMC3( B2 / RBASE, ZERO )
         D2 = ZERO
         DO 30 I = 1, BASE
            D2 = D2 + B2
   30    CONTINUE
         GO TO 10
      END IF
!+    END WHILE
!
      RETURN
!
!     End of DLAMC4
!
      END SUBROUTINE DLAMC4
!
!***********************************************************************
!
      SUBROUTINE DLAMC5( BETA, P, EMIN, IEEE, EMAX, RMAX )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            BETA, EMAX, EMIN, P
      DOUBLE PRECISION   RMAX
!     ..
!
!  Purpose
!  =======
!
!  DLAMC5 attempts to compute RMAX, the largest machine floating-point
!  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
!  approximately to a power of 2.  It will fail on machines where this
!  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
!  EMAX = 28718).  It will also fail if the value supplied for EMIN is
!  too large (i.e. too close to zero), probably with overflow.
!
!  Arguments
!  =========
!
!  BETA    (input) INTEGER
!          The base of floating-point arithmetic.
!
!  P       (input) INTEGER
!          The number of base BETA digits in the mantissa of a
!          floating-point value.
!
!  EMIN    (input) INTEGER
!          The minimum exponent before (gradual) underflow.
!
!  IEEE    (input) LOGICAL
!          A logical flag specifying whether or not the arithmetic
!          system is thought to comply with the IEEE standard.
!
!  EMAX    (output) INTEGER
!          The largest exponent before overflow
!
!  RMAX    (output) DOUBLE PRECISION
!          The largest machine floating-point number.
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      INTEGER            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
      DOUBLE PRECISION   OLDY, RECBAS, Y, Z
!     ..
!     .. External Functions ..
!      DOUBLE PRECISION   DLAMC3
!      EXTERNAL           DLAMC3
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MOD
!     ..
!     .. Executable Statements ..
!
!     First compute LEXP and UEXP, two powers of 2 that bound
!     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
!     approximately to the bound that is closest to abs(EMIN).
!     (EMAX is the exponent of the required number RMAX).
!
      LEXP = 1
      EXBITS = 1
   10 CONTINUE
      TRY = LEXP*2
      IF( TRY.LE.( -EMIN ) ) THEN
         LEXP = TRY
         EXBITS = EXBITS + 1
         GO TO 10
      END IF
      IF( LEXP.EQ.-EMIN ) THEN
         UEXP = LEXP
      ELSE
         UEXP = TRY
         EXBITS = EXBITS + 1
      END IF
!
!     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
!     than or equal to EMIN. EXBITS is the number of bits needed to
!     store the exponent.
!
      IF( ( UEXP+EMIN ).GT.( -LEXP-EMIN ) ) THEN
         EXPSUM = 2*LEXP
      ELSE
         EXPSUM = 2*UEXP
      END IF
!
!     EXPSUM is the exponent range, approximately equal to
!     EMAX - EMIN + 1 .
!
      EMAX = EXPSUM + EMIN - 1
      NBITS = 1 + EXBITS + P
!
!     NBITS is the total number of bits needed to store a
!     floating-point number.
!
      IF( ( MOD( NBITS, 2 ).EQ.1 ) .AND. ( BETA.EQ.2 ) ) THEN
!
!        Either there are an odd number of bits used to store a
!        floating-point number, which is unlikely, or some bits are
!        not used in the representation of numbers, which is possible,
!        (e.g. Cray machines) or the mantissa has an implicit bit,
!        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
!        most likely. We have to assume the last alternative.
!        If this is true, then we need to reduce EMAX by one because
!        there must be some way of representing zero in an implicit-bit
!        system. On machines like Cray, we are reducing EMAX by one
!        unnecessarily.
!
         EMAX = EMAX - 1
      END IF
!
      IF( IEEE ) THEN
!
!        Assume we are on an IEEE machine which reserves one exponent
!        for infinity and NaN.
!
         EMAX = EMAX - 1
      END IF
!
!     Now create RMAX, the largest machine number, which should
!     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
!
!     First compute 1.0 - BETA**(-P), being careful that the
!     result is less than 1.0 .
!
      RECBAS = ONE / BETA
      Z = BETA - ONE
      Y = ZERO
      DO 20 I = 1, P
         Z = Z*RECBAS
         IF( Y.LT.ONE ) &
            OLDY = Y
         Y = DLAMC3( Y, Z )
   20 CONTINUE
      IF( Y.GE.ONE ) &
         Y = OLDY
!
!     Now multiply by BETA**EMAX to get RMAX.
!
      DO 30 I = 1, EMAX
         Y = DLAMC3( Y*BETA, ZERO )
   30 CONTINUE
!
      RMAX = Y
      RETURN
!
!     End of DLAMC5
!
      END SUBROUTINE DLAMC5

      SUBROUTINE DGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            INFO, KL, KU, LDAB, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   AB( LDAB, * )
!     ..
!
!  Purpose
!  =======
!
!  DGBTRF computes an LU factorization of a real m-by-n band matrix A
!  using partial pivoting with row interchanges.
!
!  This is the blocked version of the algorithm, calling Level 3 BLAS.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  KL      (input) INTEGER
!          The number of subdiagonals within the band of A.  KL >= 0.
!
!  KU      (input) INTEGER
!          The number of superdiagonals within the band of A.  KU >= 0.
!
!  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
!          On entry, the matrix A in band storage, in rows KL+1 to
!          2*KL+KU+1; rows 1 to KL of the array need not be set.
!          The j-th column of A is stored in the j-th column of the
!          array AB as follows:
!          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
!
!          On exit, details of the factorization: U is stored as an
!          upper triangular band matrix with KL+KU superdiagonals in
!          rows 1 to KL+KU+1, and the multipliers used during the
!          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
!          See below for further details.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  Further Details
!  ===============
!
!  The band storage scheme is illustrated by the following example, when
!  M = N = 6, KL = 2, KU = 1:
!
!  On entry:                       On exit:
!
!      *    *    *    +    +    +       *    *    *   u14  u25  u36
!      *    *    +    +    +    +       *    *   u13  u24  u35  u46
!      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
!     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
!
!  Array elements marked * are not used by the routine; elements marked
!  + need not be set on entry, but are required by the routine to store
!  elements of U because of fill-in resulting from the row interchanges.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      INTEGER            NBMAX, LDWORK
      PARAMETER          ( NBMAX = 64, LDWORK = NBMAX+1 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, I2, I3, II, IP, J, J2, J3, JB, JJ, JM, JP, &
                         JU, K2, KM, KV, NB, NW
      DOUBLE PRECISION   TEMP
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   WORK13( LDWORK, NBMAX ), &
                         WORK31( LDWORK, NBMAX )
!     ..
!     .. External Functions ..
!      INTEGER            IDAMAX, ILAENV
!      EXTERNAL           IDAMAX, ILAENV
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DCOPY, DGBTF2, DGEMM, DGER, DLASWP, DSCAL, &
!                         DSWAP, DTRSM, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     KV is the number of superdiagonals in the factor U, allowing for
!     fill-in
!
      KV = KU + KL
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KL+KV+1 ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGBTRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 ) &
         RETURN
!
!     Determine the block size for this environment
!
      NB = ILAENV( 1, 'DGBTRF', ' ', M, N, KL, KU )
!
!     The block size must not exceed the limit set by the size of the
!     local arrays WORK13 and WORK31.
!
      NB = MIN( NB, NBMAX )
!
      IF( NB.LE.1 .OR. NB.GT.KL ) THEN
!
!        Use unblocked code
!
         CALL DGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
      ELSE
!
!        Use blocked code
!
!        Zero the superdiagonal elements of the work array WORK13
!
         DO 20 J = 1, NB
            DO 10 I = 1, J - 1
               WORK13( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
!
!        Zero the subdiagonal elements of the work array WORK31
!
         DO 40 J = 1, NB
            DO 30 I = J + 1, NB
               WORK31( I, J ) = ZERO
   30       CONTINUE
   40    CONTINUE
!
!        Gaussian elimination with partial pivoting
!
!        Set fill-in elements in columns KU+2 to KV to zero
!
         DO 60 J = KU + 2, MIN( KV, N )
            DO 50 I = KV - J + 2, KL
               AB( I, J ) = ZERO
   50       CONTINUE
   60    CONTINUE
!
!        JU is the index of the last column affected by the current
!        stage of the factorization
!
         JU = 1
!
         DO 180 J = 1, MIN( M, N ), NB
            JB = MIN( NB, MIN( M, N )-J+1 )
!
!           The active part of the matrix is partitioned
!
!              A11   A12   A13
!              A21   A22   A23
!              A31   A32   A33
!
!           Here A11, A21 and A31 denote the current block of JB columns
!           which is about to be factorized. The number of rows in the
!           partitioning are JB, I2, I3 respectively, and the numbers
!           of columns are JB, J2, J3. The superdiagonal elements of A13
!           and the subdiagonal elements of A31 lie outside the band.
!
            I2 = MIN( KL-JB, M-J-JB+1 )
            I3 = MIN( JB, M-J-KL+1 )
!
!           J2 and J3 are computed after JU has been updated.
!
!           Factorize the current block of JB columns
!
            DO 80 JJ = J, J + JB - 1
!
!              Set fill-in elements in column JJ+KV to zero
!
               IF( JJ+KV.LE.N ) THEN
                  DO 70 I = 1, KL
                     AB( I, JJ+KV ) = ZERO
   70             CONTINUE
               END IF
!
!              Find pivot and test for singularity. KM is the number of
!              subdiagonal elements in the current column.
!
               KM = MIN( KL, M-JJ )
               JP = IDAMAX( KM+1, AB( KV+1, JJ ), 1 )
               IPIV( JJ ) = JP + JJ - J
               IF( AB( KV+JP, JJ ).NE.ZERO ) THEN
                  JU = MAX( JU, MIN( JJ+KU+JP-1, N ) )
                  IF( JP.NE.1 ) THEN
!
!                    Apply interchange to columns J to J+JB-1
!
                     IF( JP+JJ-1.LT.J+KL ) THEN
!
                        CALL DSWAP( JB, AB( KV+1+JJ-J, J ), LDAB-1, &
                                    AB( KV+JP+JJ-J, J ), LDAB-1 )
                     ELSE
!
!                       The interchange affects columns J to JJ-1 of A31
!                       which are stored in the work array WORK31
!
                        CALL DSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1, &
                                    WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                        CALL DSWAP( J+JB-JJ, AB( KV+1, JJ ), LDAB-1, &
                                    AB( KV+JP, JJ ), LDAB-1 )
                     END IF
                  END IF
!
!                 Compute multipliers
!
                  CALL DSCAL( KM, ONE / AB( KV+1, JJ ), AB( KV+2, JJ ), &
                              1 )
!
!                 Update trailing submatrix within the band and within
!                 the current block. JM is the index of the last column
!                 which needs to be updated.
!
                  JM = MIN( JU, J+JB-1 )
                  IF( JM.GT.JJ ) &
                     CALL DGER( KM, JM-JJ, -ONE, AB( KV+2, JJ ), 1, &
                                AB( KV, JJ+1 ), LDAB-1, &
                                AB( KV+1, JJ+1 ), LDAB-1 )
               ELSE
!
!                 If pivot is zero, set INFO to the index of the pivot
!                 unless a zero pivot has already been found.
!
                  IF( INFO.EQ.0 ) &
                     INFO = JJ
               END IF
!
!              Copy current column of A31 into the work array WORK31
!
               NW = MIN( JJ-J+1, I3 )
               IF( NW.GT.0 ) &
                  CALL DCOPY( NW, AB( KV+KL+1-JJ+J, JJ ), 1, &
                              WORK31( 1, JJ-J+1 ), 1 )
   80       CONTINUE
            IF( J+JB.LE.N ) THEN
!
!              Apply the row interchanges to the other blocks.
!
               J2 = MIN( JU-J+1, KV ) - JB
               J3 = MAX( 0, JU-J-KV+1 )
!
!              Use DLASWP to apply the row interchanges to A12, A22, and
!              A32.
!
               CALL DLASWP( J2, AB( KV+1-JB, J+JB ), LDAB-1, 1, JB, &
                            IPIV( J ), 1 )
!
!              Adjust the pivot indices.
!
               DO 90 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
   90          CONTINUE
!
!              Apply the row interchanges to A13, A23, and A33
!              columnwise.
!
               K2 = J - 1 + JB + J2
               DO 110 I = 1, J3
                  JJ = K2 + I
                  DO 100 II = J + I - 1, J + JB - 1
                     IP = IPIV( II )
                     IF( IP.NE.II ) THEN
                        TEMP = AB( KV+1+II-JJ, JJ )
                        AB( KV+1+II-JJ, JJ ) = AB( KV+1+IP-JJ, JJ )
                        AB( KV+1+IP-JJ, JJ ) = TEMP
                     END IF
  100             CONTINUE
  110          CONTINUE
!
!              Update the relevant part of the trailing submatrix
!
               IF( J2.GT.0 ) THEN
!
!                 Update A12
!
                  CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', &
                              JB, J2, ONE, AB( KV+1, J ), LDAB-1, &
                              AB( KV+1-JB, J+JB ), LDAB-1 )
!
                  IF( I2.GT.0 ) THEN
!
!                    Update A22
!
                     CALL DGEMM( 'No transpose', 'No transpose', I2, J2, &
                                 JB, -ONE, AB( KV+1+JB, J ), LDAB-1, &
                                 AB( KV+1-JB, J+JB ), LDAB-1, ONE, &
                                 AB( KV+1, J+JB ), LDAB-1 )
                  END IF
!
                  IF( I3.GT.0 ) THEN
!
!                    Update A32
!
                     CALL DGEMM( 'No transpose', 'No transpose', I3, J2, &
                                 JB, -ONE, WORK31, LDWORK, &
                                 AB( KV+1-JB, J+JB ), LDAB-1, ONE, &
                                 AB( KV+KL+1-JB, J+JB ), LDAB-1 )
                  END IF
               END IF
!
               IF( J3.GT.0 ) THEN
!
!                 Copy the lower triangle of A13 into the work array
!                 WORK13
!
                  DO 130 JJ = 1, J3
                     DO 120 II = JJ, JB
                        WORK13( II, JJ ) = AB( II-JJ+1, JJ+J+KV-1 )
  120                CONTINUE
  130             CONTINUE
!
!                 Update A13 in the work array
!
                  CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', &
                              JB, J3, ONE, AB( KV+1, J ), LDAB-1, &
                              WORK13, LDWORK )
!
                  IF( I2.GT.0 ) THEN
!
!                    Update A23
!
                     CALL DGEMM( 'No transpose', 'No transpose', I2, J3, &
                                 JB, -ONE, AB( KV+1+JB, J ), LDAB-1, &
                                 WORK13, LDWORK, ONE, AB( 1+JB, J+KV ), &
                                 LDAB-1 )
                  END IF
!
                  IF( I3.GT.0 ) THEN
!
!                    Update A33
!
                     CALL DGEMM( 'No transpose', 'No transpose', I3, J3, &
                                 JB, -ONE, WORK31, LDWORK, WORK13, &
                                 LDWORK, ONE, AB( 1+KL, J+KV ), LDAB-1 )
                  END IF
!
!                 Copy the lower triangle of A13 back into place
!
                  DO 150 JJ = 1, J3
                     DO 140 II = JJ, JB
                        AB( II-JJ+1, JJ+J+KV-1 ) = WORK13( II, JJ )
  140                CONTINUE
  150             CONTINUE
               END IF
            ELSE
!
!              Adjust the pivot indices.
!
               DO 160 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
  160          CONTINUE
            END IF
!
!           Partially undo the interchanges in the current block to
!           restore the upper triangular form of A31 and copy the upper
!           triangle of A31 back into place
!
            DO 170 JJ = J + JB - 1, J, -1
               JP = IPIV( JJ ) - JJ + 1
               IF( JP.NE.1 ) THEN
!
!                 Apply interchange to columns J to JJ-1
!
                  IF( JP+JJ-1.LT.J+KL ) THEN
!
!                    The interchange does not affect A31
!
                     CALL DSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1, &
                                 AB( KV+JP+JJ-J, J ), LDAB-1 )
                  ELSE
!
!                    The interchange does affect A31
!
                     CALL DSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1, &
                                 WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                  END IF
               END IF
!
!              Copy the current column of A31 back into place
!
               NW = MIN( I3, JJ-J+1 )
               IF( NW.GT.0 ) &
                  CALL DCOPY( NW, WORK31( 1, JJ-J+1 ), 1, &
                              AB( KV+KL+1-JJ+J, JJ ), 1 )
  170       CONTINUE
  180    CONTINUE
      END IF
!
      RETURN
!
!     End of DGBTRF
!
      END SUBROUTINE DGBTRF


      SUBROUTINE DGBCON( NORM, N, KL, KU, AB, LDAB, IPIV, ANORM, RCOND, &
                         WORK, IWORK, INFO )
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     Modified to call DLACN2 in place of DLACON, 5 Feb 03, SJH.
!
!     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            INFO, KL, KU, LDAB, N
      DOUBLE PRECISION   ANORM, RCOND
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * ), IWORK( * )
      DOUBLE PRECISION   AB( LDAB, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DGBCON estimates the reciprocal of the condition number of a real
!  general band matrix A, in either the 1-norm or the infinity-norm,
!  using the LU factorization computed by DGBTRF.
!
!  An estimate is obtained for norm(inv(A)), and the reciprocal of the
!  condition number is computed as
!     RCOND = 1 / ( norm(A) * norm(inv(A)) ).
!
!  Arguments
!  =========
!
!  NORM    (input) CHARACTER*1
!          Specifies whether the 1-norm condition number or the
!          infinity-norm condition number is required:
!          = '1' or 'O':  1-norm;
!          = 'I':         Infinity-norm.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  KL      (input) INTEGER
!          The number of subdiagonals within the band of A.  KL >= 0.
!
!  KU      (input) INTEGER
!          The number of superdiagonals within the band of A.  KU >= 0.
!
!  AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)
!          Details of the LU factorization of the band matrix A, as
!          computed by DGBTRF.  U is stored as an upper triangular band
!          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
!          the multipliers used during the factorization are stored in
!          rows KL+KU+2 to 2*KL+KU+1.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices; for 1 <= i <= N, row i of the matrix was
!          interchanged with row IPIV(i).
!
!  ANORM   (input) DOUBLE PRECISION
!          If NORM = '1' or 'O', the 1-norm of the original matrix A.
!          If NORM = 'I', the infinity-norm of the original matrix A.
!
!  RCOND   (output) DOUBLE PRECISION
!          The reciprocal of the condition number of the matrix A,
!          computed as RCOND = 1/(norm(A) * norm(inv(A))).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (3*N)
!
!  IWORK   (workspace) INTEGER array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LNOTI, ONENRM
      CHARACTER          NORMIN
      INTEGER            IX, J, JP, KASE, KASE1, KD, LM
      DOUBLE PRECISION   AINVNM, SCALE, SMLNUM, T
!     ..
!     .. Local Arrays ..
      INTEGER            ISAVE( 3 )
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      INTEGER            IDAMAX
!      DOUBLE PRECISION   DDOT, DLAMCH
!      EXTERNAL           LSAME, IDAMAX, DDOT, DLAMCH
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DAXPY, DLACN2, DLATBS, DRSCL, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      ONENRM = NORM.EQ.'1' .OR. LSAME( NORM, 'O' )
      IF( .NOT.ONENRM .AND. .NOT.LSAME( NORM, 'I' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.2*KL+KU+1 ) THEN
         INFO = -6
      ELSE IF( ANORM.LT.ZERO ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGBCON', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      RCOND = ZERO
      IF( N.EQ.0 ) THEN
         RCOND = ONE
         RETURN
      ELSE IF( ANORM.EQ.ZERO ) THEN
         RETURN
      END IF
!
      SMLNUM = DLAMCH( 'Safe minimum' )
!
!     Estimate the norm of inv(A).
!
      AINVNM = ZERO
      NORMIN = 'N'
      IF( ONENRM ) THEN
         KASE1 = 1
      ELSE
         KASE1 = 2
      END IF
      KD = KL + KU + 1
      LNOTI = KL.GT.0
      KASE = 0
   10 CONTINUE
      CALL DLACN2( N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE )
      IF( KASE.NE.0 ) THEN
         IF( KASE.EQ.KASE1 ) THEN
!
!           Multiply by inv(L).
!
            IF( LNOTI ) THEN
               DO 20 J = 1, N - 1
                  LM = MIN( KL, N-J )
                  JP = IPIV( J )
                  T = WORK( JP )
                  IF( JP.NE.J ) THEN
                     WORK( JP ) = WORK( J )
                     WORK( J ) = T
                  END IF
                  CALL DAXPY( LM, -T, AB( KD+1, J ), 1, WORK( J+1 ), 1 )
   20          CONTINUE
            END IF
!
!           Multiply by inv(U).
!
            CALL DLATBS( 'Upper', 'No transpose', 'Non-unit', NORMIN, N, &
                         KL+KU, AB, LDAB, WORK, SCALE, WORK( 2*N+1 ), &
                         INFO )
         ELSE
!
!           Multiply by inv(U').
!
            CALL DLATBS( 'Upper', 'Transpose', 'Non-unit', NORMIN, N, &
                         KL+KU, AB, LDAB, WORK, SCALE, WORK( 2*N+1 ), &
                         INFO )
!
!           Multiply by inv(L').
!
            IF( LNOTI ) THEN
               DO 30 J = N - 1, 1, -1
                  LM = MIN( KL, N-J )
                  WORK( J ) = WORK( J ) - DDOT( LM, AB( KD+1, J ), 1, &
                              WORK( J+1 ), 1 )
                  JP = IPIV( J )
                  IF( JP.NE.J ) THEN
                     T = WORK( JP )
                     WORK( JP ) = WORK( J )
                     WORK( J ) = T
                  END IF
   30          CONTINUE
            END IF
         END IF
!
!        Divide X by 1/SCALE if doing so will not cause overflow.
!
         NORMIN = 'Y'
         IF( SCALE.NE.ONE ) THEN
            IX = IDAMAX( N, WORK, 1 )
            IF( SCALE.LT.ABS( WORK( IX ) )*SMLNUM .OR. SCALE.EQ.ZERO ) &
               GO TO 40
            CALL DRSCL( N, SCALE, WORK, 1 )
         END IF
         GO TO 10
      END IF
!
!     Compute the estimate of the reciprocal condition number.
!
      IF( AINVNM.NE.ZERO ) &
         RCOND = ( ONE / AINVNM ) / ANORM
!
   40 CONTINUE
      RETURN
!
!     End of DGBCON
!
      END SUBROUTINE DGBCON


      SUBROUTINE DGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, &
                         INFO )
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   AB( LDAB, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DGBTRS solves a system of linear equations
!     A * X = B  or  A' * X = B
!  with a general band matrix A using the LU factorization computed
!  by DGBTRF.
!
!  Arguments
!  =========
!
!  TRANS   (input) CHARACTER*1
!          Specifies the form of the system of equations.
!          = 'N':  A * X = B  (No transpose)
!          = 'T':  A'* X = B  (Transpose)
!          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  KL      (input) INTEGER
!          The number of subdiagonals within the band of A.  KL >= 0.
!
!  KU      (input) INTEGER
!          The number of superdiagonals within the band of A.  KU >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)
!          Details of the LU factorization of the band matrix A, as
!          computed by DGBTRF.  U is stored as an upper triangular band
!          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
!          the multipliers used during the factorization are stored in
!          rows KL+KU+2 to 2*KL+KU+1.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices; for 1 <= i <= N, row i of the matrix was
!          interchanged with row IPIV(i).
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LNOTI, NOTRAN
      INTEGER            I, J, KD, L, LM
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DGEMV, DGER, DSWAP, DTBSV, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. &
          LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDAB.LT.( 2*KL+KU+1 ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGBTRS', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) &
         RETURN
!
      KD = KU + KL + 1
      LNOTI = KL.GT.0
!
      IF( NOTRAN ) THEN
!
!        Solve  A*X = B.
!
!        Solve L*X = B, overwriting B with X.
!
!        L is represented as a product of permutations and unit lower
!        triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
!        where each transformation L(i) is a rank-one modification of
!        the identity matrix.
!
         IF( LNOTI ) THEN
            DO 10 J = 1, N - 1
               LM = MIN( KL, N-J )
               L = IPIV( J )
               IF( L.NE.J ) &
                  CALL DSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB )
               CALL DGER( LM, NRHS, -ONE, AB( KD+1, J ), 1, B( J, 1 ), &
                          LDB, B( J+1, 1 ), LDB )
   10       CONTINUE
         END IF
!
         DO 20 I = 1, NRHS
!
!           Solve U*X = B, overwriting B with X.
!
            CALL DTBSV( 'Upper', 'No transpose', 'Non-unit', N, KL+KU, &
                        AB, LDAB, B( 1, I ), 1 )
   20    CONTINUE
!
      ELSE
!
!        Solve A'*X = B.
!
         DO 30 I = 1, NRHS
!
!           Solve U'*X = B, overwriting B with X.
!
            CALL DTBSV( 'Upper', 'Transpose', 'Non-unit', N, KL+KU, AB, &
                        LDAB, B( 1, I ), 1 )
   30    CONTINUE
!
!        Solve L'*X = B, overwriting B with X.
!
         IF( LNOTI ) THEN
            DO 40 J = N - 1, 1, -1
               LM = MIN( KL, N-J )
               CALL DGEMV( 'Transpose', LM, NRHS, -ONE, B( J+1, 1 ), &
                           LDB, AB( KD+1, J ), 1, ONE, B( J, 1 ), LDB )
               L = IPIV( J )
               IF( L.NE.J ) &
                  CALL DSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB )
   40       CONTINUE
         END IF
      END IF
      RETURN
!
!     End of DGBTRS
!
      END SUBROUTINE DGBTRS


      SUBROUTINE DGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, IWORK, &
                         INFO )
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     Modified to call DLACN2 in place of DLACON, 5 Feb 03, SJH.
!
!     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            INFO, LDA, N
      DOUBLE PRECISION   ANORM, RCOND
!     ..
!     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DGECON estimates the reciprocal of the condition number of a general
!  real matrix A, in either the 1-norm or the infinity-norm, using
!  the LU factorization computed by DGETRF.
!
!  An estimate is obtained for norm(inv(A)), and the reciprocal of the
!  condition number is computed as
!     RCOND = 1 / ( norm(A) * norm(inv(A)) ).
!
!  Arguments
!  =========
!
!  NORM    (input) CHARACTER*1
!          Specifies whether the 1-norm condition number or the
!          infinity-norm condition number is required:
!          = '1' or 'O':  1-norm;
!          = 'I':         Infinity-norm.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!          The factors L and U from the factorization A = P*L*U
!          as computed by DGETRF.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  ANORM   (input) DOUBLE PRECISION
!          If NORM = '1' or 'O', the 1-norm of the original matrix A.
!          If NORM = 'I', the infinity-norm of the original matrix A.
!
!  RCOND   (output) DOUBLE PRECISION
!          The reciprocal of the condition number of the matrix A,
!          computed as RCOND = 1/(norm(A) * norm(inv(A))).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (4*N)
!
!  IWORK   (workspace) INTEGER array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            ONENRM
      CHARACTER          NORMIN
      INTEGER            IX, KASE, KASE1
      DOUBLE PRECISION   AINVNM, SCALE, SL, SMLNUM, SU
!     ..
!     .. Local Arrays ..
      INTEGER            ISAVE( 3 )
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      INTEGER            IDAMAX
!      DOUBLE PRECISION   DLAMCH
!     EXTERNAL           LSAME, IDAMAX, DLAMCH
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DLACN2, DLATRS, DRSCL, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      ONENRM = NORM.EQ.'1' .OR. LSAME( NORM, 'O' )
      IF( .NOT.ONENRM .AND. .NOT.LSAME( NORM, 'I' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( ANORM.LT.ZERO ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGECON', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      RCOND = ZERO
      IF( N.EQ.0 ) THEN
         RCOND = ONE
         RETURN
      ELSE IF( ANORM.EQ.ZERO ) THEN
         RETURN
      END IF
!
      SMLNUM = DLAMCH( 'Safe minimum' )
!
!     Estimate the norm of inv(A).
!
      AINVNM = ZERO
      NORMIN = 'N'
      IF( ONENRM ) THEN
         KASE1 = 1
      ELSE
         KASE1 = 2
      END IF
      KASE = 0
   10 CONTINUE
      CALL DLACN2( N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE )
      IF( KASE.NE.0 ) THEN
         IF( KASE.EQ.KASE1 ) THEN
!
!           Multiply by inv(L).
!
            CALL DLATRS( 'Lower', 'No transpose', 'Unit', NORMIN, N, A, &
                         LDA, WORK, SL, WORK( 2*N+1 ), INFO )
!
!           Multiply by inv(U).
!
            CALL DLATRS( 'Upper', 'No transpose', 'Non-unit', NORMIN, N, &
                         A, LDA, WORK, SU, WORK( 3*N+1 ), INFO )
         ELSE
!
!           Multiply by inv(U').
!
            CALL DLATRS( 'Upper', 'Transpose', 'Non-unit', NORMIN, N, A, &
                         LDA, WORK, SU, WORK( 3*N+1 ), INFO )
!
!           Multiply by inv(L').
!
            CALL DLATRS( 'Lower', 'Transpose', 'Unit', NORMIN, N, A, &
                         LDA, WORK, SL, WORK( 2*N+1 ), INFO )
         END IF
!
!        Divide X by 1/(SL*SU) if doing so will not cause overflow.
!
         SCALE = SL*SU
         NORMIN = 'Y'
         IF( SCALE.NE.ONE ) THEN
            IX = IDAMAX( N, WORK, 1 )
            IF( SCALE.LT.ABS( WORK( IX ) )*SMLNUM .OR. SCALE.EQ.ZERO ) &
               GO TO 20
            CALL DRSCL( N, SCALE, WORK, 1 )
         END IF
         GO TO 10
      END IF
!
!     Compute the estimate of the reciprocal condition number.
!
      IF( AINVNM.NE.ZERO ) &
         RCOND = ( ONE / AINVNM ) / ANORM
!
   20 CONTINUE
      RETURN
!
!     End of DGECON
!
      END SUBROUTINE DGECON


      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DGETRF computes an LU factorization of a general M-by-N matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 3 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the M-by-N matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
!                has been completed, but the factor U is exactly
!                singular, and division by zero will occur if it is used
!                to solve a system of equations.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DGEMM, DGETF2, DLASWP, DTRSM, XERBLA
!     ..
!     .. External Functions ..
!      INTEGER            ILAENV
!      EXTERNAL           ILAENV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 ) &
         RETURN
!
!     Determine the block size for this environment.
!
      NB = ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
!
!        Use unblocked code.
!
         CALL DGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE
!
!        Use blocked code.
!
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
!
!           Factor diagonal and subdiagonal blocks and test for exact
!           singularity.
!
            CALL DGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
!
!           Adjust INFO and the pivot indices.
!
            IF( INFO.EQ.0 .AND. IINFO.GT.0 ) &
               INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
!
!           Apply interchanges to columns 1:J-1.
!
            CALL DLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
!
            IF( J+JB.LE.N ) THEN
!
!              Apply interchanges to columns J+JB:N.
!
               CALL DLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1, &
                            IPIV, 1 )
!
!              Compute block row of U.
!
               CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB, &
                           N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ), &
                           LDA )
               IF( J+JB.LE.M ) THEN
!
!                 Update trailing submatrix.
!
                  CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1, &
                              N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA, &
                              A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ), &
                              LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
!
!     End of DGETRF
!
      END SUBROUTINE DGETRF


      SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DGETRS solves a system of linear equations
!     A * X = B  or  A' * X = B
!  with a general N-by-N matrix A using the LU factorization computed
!  by DGETRF.
!
!  Arguments
!  =========
!
!  TRANS   (input) CHARACTER*1
!          Specifies the form of the system of equations:
!          = 'N':  A * X = B  (No transpose)
!          = 'T':  A'* X = B  (Transpose)
!          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!          The factors L and U from the factorization A = P*L*U
!          as computed by DGETRF.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices from DGETRF; for 1<=i<=N, row i of the
!          matrix was interchanged with row IPIV(i).
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOTRAN
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DLASWP, DTRSM, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. &
          LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRS', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) &
         RETURN
!
      IF( NOTRAN ) THEN
!
!        Solve A * X = B.
!
!        Apply row interchanges to the right hand sides.
!
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
!
!        Solve L*X = B, overwriting B with X.
!
         CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS, &
                     ONE, A, LDA, B, LDB )
!
!        Solve U*X = B, overwriting B with X.
!
         CALL DTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N, &
                     NRHS, ONE, A, LDA, B, LDB )
      ELSE
!
!        Solve A' * X = B.
!
!        Solve U'*X = B, overwriting B with X.
!
         CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS, &
                     ONE, A, LDA, B, LDB )
!
!        Solve L'*X = B, overwriting B with X.
!
         CALL DTRSM( 'Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE, &
                     A, LDA, B, LDB )
!
!        Apply row interchanges to the solution vectors.
!
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF
!
      RETURN
!
!     End of DGETRS
!
      END SUBROUTINE DGETRS


      DOUBLE PRECISION FUNCTION DLANGB( NORM, N, KL, KU, AB, LDAB, &
                       WORK )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            KL, KU, LDAB, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   AB( LDAB, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DLANGB  returns the value of the one norm,  or the Frobenius norm, or
!  the  infinity norm,  or the element of  largest absolute value  of an
!  n by n band matrix  A,  with kl sub-diagonals and ku super-diagonals.
!
!  Description
!  ===========
!
!  DLANGB returns the value
!
!     DLANGB = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!              (
!              ( norm1(A),         NORM = '1', 'O' or 'o'
!              (
!              ( normI(A),         NORM = 'I' or 'i'
!              (
!              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!
!  where  norm1  denotes the  one norm of a matrix (maximum column sum),
!  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!  normF  denotes the  Frobenius norm of a matrix (square root of sum of
!  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix no
!
!  Arguments
!  =========
!
!  NORM    (input) CHARACTER*1
!          Specifies the value to be returned in DLANGB as described
!          above.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.  When N = 0, DLANGB is
!          set to zero.
!
!  KL      (input) INTEGER
!          The number of sub-diagonals of the matrix A.  KL >= 0.
!
!  KU      (input) INTEGER
!          The number of super-diagonals of the matrix A.  KU >= 0.
!
!  AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)
!          The band matrix A, stored in rows 1 to KL+KU+1.  The j-th
!          column of A is stored in the j-th column of the array AB as
!          follows:
!          AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(n,j+kl).
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= KL+KU+1.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)),
!          where LWORK >= N when NORM = 'I'; otherwise, WORK is not
!          referenced.
!
! =====================================================================
!
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J, K, L
      DOUBLE PRECISION   SCALE, SUM, VALUE
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DLASSQ
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
      IF( N.EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
!
!        Find max(abs(A(i,j))).
!
         VALUE = ZERO
         DO 20 J = 1, N
            DO 10 I = MAX( KU+2-J, 1 ), MIN( N+KU+1-J, KL+KU+1 )
               VALUE = MAX( VALUE, ABS( AB( I, J ) ) )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) THEN
!
!        Find norm1(A).
!
         VALUE = ZERO
         DO 40 J = 1, N
            SUM = ZERO
            DO 30 I = MAX( KU+2-J, 1 ), MIN( N+KU+1-J, KL+KU+1 )
               SUM = SUM + ABS( AB( I, J ) )
   30       CONTINUE
            VALUE = MAX( VALUE, SUM )
   40    CONTINUE
      ELSE IF( LSAME( NORM, 'I' ) ) THEN
!
!        Find normI(A).
!
         DO 50 I = 1, N
            WORK( I ) = ZERO
   50    CONTINUE
         DO 70 J = 1, N
            K = KU + 1 - J
            DO 60 I = MAX( 1, J-KU ), MIN( N, J+KL )
               WORK( I ) = WORK( I ) + ABS( AB( K+I, J ) )
   60       CONTINUE
   70    CONTINUE
         VALUE = ZERO
         DO 80 I = 1, N
            VALUE = MAX( VALUE, WORK( I ) )
   80    CONTINUE
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
         SCALE = ZERO
         SUM = ONE
         DO 90 J = 1, N
            L = MAX( 1, J-KU )
            K = KU + 1 - J + L
            CALL DLASSQ( MIN( N, J+KL )-L+1, AB( K, J ), 1, SCALE, SUM )
   90    CONTINUE
         VALUE = SCALE*SQRT( SUM )
      END IF
!
      DLANGB = VALUE
      RETURN
!
!     End of DLANGB
!
      END FUNCTION DLANGB


      DOUBLE PRECISION FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            LDA, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DLANGE  returns the value of the one norm,  or the Frobenius norm, or
!  the  infinity norm,  or the  element of  largest absolute value  of a
!  real matrix A.
!
!  Description
!  ===========
!
!  DLANGE returns the value
!
!     DLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!              (
!              ( norm1(A),         NORM = '1', 'O' or 'o'
!              (
!              ( normI(A),         NORM = 'I' or 'i'
!              (
!              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!
!  where  norm1  denotes the  one norm of a matrix (maximum column sum),
!  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!  normF  denotes the  Frobenius norm of a matrix (square root of sum of
!  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix no
!
!  Arguments
!  =========
!
!  NORM    (input) CHARACTER*1
!          Specifies the value to be returned in DLANGE as described
!          above.
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.  When M = 0,
!          DLANGE is set to zero.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.  When N = 0,
!          DLANGE is set to zero.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!          The m by n matrix A.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(M,1).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)),
!          where LWORK >= M when NORM = 'I'; otherwise, WORK is not
!          referenced.
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   SCALE, SUM, VALUE
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DLASSQ
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
      IF( MIN( M, N ).EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
!
!        Find max(abs(A(i,j))).
!
         VALUE = ZERO
         DO 20 J = 1, N
            DO 10 I = 1, M
               VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) THEN
!
!        Find norm1(A).
!
         VALUE = ZERO
         DO 40 J = 1, N
            SUM = ZERO
            DO 30 I = 1, M
               SUM = SUM + ABS( A( I, J ) )
   30       CONTINUE
            VALUE = MAX( VALUE, SUM )
   40    CONTINUE
      ELSE IF( LSAME( NORM, 'I' ) ) THEN
!
!        Find normI(A).
!
         DO 50 I = 1, M
            WORK( I ) = ZERO
   50    CONTINUE
         DO 70 J = 1, N
            DO 60 I = 1, M
               WORK( I ) = WORK( I ) + ABS( A( I, J ) )
   60       CONTINUE
   70    CONTINUE
         VALUE = ZERO
         DO 80 I = 1, M
            VALUE = MAX( VALUE, WORK( I ) )
   80    CONTINUE
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
         SCALE = ZERO
         SUM = ONE
         DO 90 J = 1, N
            CALL DLASSQ( M, A( 1, J ), 1, SCALE, SUM )
   90    CONTINUE
         VALUE = SCALE*SQRT( SUM )
      END IF
!
      DLANGE = VALUE
      RETURN
!
!     End of DLANGE
!
      END FUNCTION DLANGE


      SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, &
                        LDVR, WORK, LWORK, INFO )
!
!  -- LAPACK driver routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          JOBVL, JOBVR
      INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), &
                         WI( * ), WORK( * ), WR( * )
!     ..
!
!  Purpose
!  =======
!
!  DGEEV computes for an N-by-N real nonsymmetric matrix A, the
!  eigenvalues and, optionally, the left and/or right eigenvectors.
!
!  The right eigenvector v(j) of A satisfies
!                   A * v(j) = lambda(j) * v(j)
!  where lambda(j) is its eigenvalue.
!  The left eigenvector u(j) of A satisfies
!                u(j)**H * A = lambda(j) * u(j)**H
!  where u(j)**H denotes the conjugate transpose of u(j).
!
!  The computed eigenvectors are normalized to have Euclidean norm
!  equal to 1 and largest component real.
!
!  Arguments
!  =========
!
!  JOBVL   (input) CHARACTER*1
!          = 'N': left eigenvectors of A are not computed;
!          = 'V': left eigenvectors of A are computed.
!
!  JOBVR   (input) CHARACTER*1
!          = 'N': right eigenvectors of A are not computed;
!          = 'V': right eigenvectors of A are computed.
!
!  N       (input) INTEGER
!          The order of the matrix A. N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the N-by-N matrix A.
!          On exit, A has been overwritten.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  WR      (output) DOUBLE PRECISION array, dimension (N)
!  WI      (output) DOUBLE PRECISION array, dimension (N)
!          WR and WI contain the real and imaginary parts,
!          respectively, of the computed eigenvalues.  Complex
!          conjugate pairs of eigenvalues appear consecutively
!          with the eigenvalue having the positive imaginary part
!          first.
!
!  VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)
!          If JOBVL = 'V', the left eigenvectors u(j) are stored one
!          after another in the columns of VL, in the same order
!          as their eigenvalues.
!          If JOBVL = 'N', VL is not referenced.
!          If the j-th eigenvalue is real, then u(j) = VL(:,j),
!          the j-th column of VL.
!          If the j-th and (j+1)-st eigenvalues form a complex
!          conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and
!          u(j+1) = VL(:,j) - i*VL(:,j+1).
!
!  LDVL    (input) INTEGER
!          The leading dimension of the array VL.  LDVL >= 1; if
!          JOBVL = 'V', LDVL >= N.
!
!  VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)
!          If JOBVR = 'V', the right eigenvectors v(j) are stored one
!          after another in the columns of VR, in the same order
!          as their eigenvalues.
!          If JOBVR = 'N', VR is not referenced.
!          If the j-th eigenvalue is real, then v(j) = VR(:,j),
!          the j-th column of VR.
!          If the j-th and (j+1)-st eigenvalues form a complex
!          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
!          v(j+1) = VR(:,j) - i*VR(:,j+1).
!
!  LDVR    (input) INTEGER
!          The leading dimension of the array VR.  LDVR >= 1; if
!          JOBVR = 'V', LDVR >= N.
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,L
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.  LWORK >= max(1,3*N), and
!          if JOBVL = 'V' or JOBVR = 'V', LWORK >= 4*N.  For good
!          performance, LWORK must generally be larger.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  if INFO = i, the QR algorithm failed to compute all the
!                eigenvalues, and no eigenvectors have been computed;
!                elements i+1:N of WR and WI contain eigenvalues which
!                have converged.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY, SCALEA, WANTVL, WANTVR
      CHARACTER          SIDE
      INTEGER            HSWORK, I, IBAL, IERR, IHI, ILO, ITAU, IWRK, K, &
                         MAXWRK, MINWRK, NOUT
      DOUBLE PRECISION   ANRM, BIGNUM, CS, CSCALE, EPS, R, SCL, SMLNUM, &
                         SN
!     ..
!     .. Local Arrays ..
      LOGICAL            SELECT( 1 )
      DOUBLE PRECISION   DUM( 1 )
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DGEBAK, DGEBAL, DGEHRD, DHSEQR, DLABAD, DLACPY, &
!                         DLARTG, DLASCL, DORGHR, DROT, DSCAL, DTREVC, &
!                         XERBLA
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      INTEGER            IDAMAX, ILAENV
!      DOUBLE PRECISION   DLAMCH, DLANGE, DLAPY2, DNRM2
!      EXTERNAL           LSAME, IDAMAX, ILAENV, DLAMCH, DLANGE, DLAPY2, &
!                         DNRM2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      WANTVL = LSAME( JOBVL, 'V' )
      WANTVR = LSAME( JOBVR, 'V' )
      IF( ( .NOT.WANTVL ) .AND. ( .NOT.LSAME( JOBVL, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( ( .NOT.WANTVR ) .AND. ( .NOT.LSAME( JOBVR, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDVL.LT.1 .OR. ( WANTVL .AND. LDVL.LT.N ) ) THEN
         INFO = -9
      ELSE IF( LDVR.LT.1 .OR. ( WANTVR .AND. LDVR.LT.N ) ) THEN
         INFO = -11
      END IF
!
!     Compute workspace
!      (Note: Comments in the code beginning "Workspace:" describe the
!       minimal amount of workspace needed at that point in the code,
!       as well as the preferred amount for good performance.
!       NB refers to the optimal block size for the immediately
!       following subroutine, as returned by ILAENV.
!       HSWORK refers to the workspace preferred by DHSEQR, as
!       calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
!       the worst case.)
!
      IF( INFO.EQ.0 ) THEN
         IF( N.EQ.0 ) THEN
            MINWRK = 1
            MAXWRK = 1
         ELSE
            MAXWRK = 2*N + N*ILAENV( 1, 'DGEHRD', ' ', N, 1, N, 0 )
            IF( WANTVL ) THEN
               MINWRK = 4*N
               MAXWRK = MAX( MAXWRK, 2*N + ( N - 1 )*ILAENV( 1, &
                             'DORGHR', ' ', N, 1, N, -1 ) )
               CALL DHSEQR( 'S', 'V', N, 1, N, A, LDA, WR, WI, VL, LDVL, &
                      WORK, -1, INFO )
               HSWORK = WORK( 1 )
               MAXWRK = MAX( MAXWRK, N + 1, N + HSWORK )
               MAXWRK = MAX( MAXWRK, 4*N )
            ELSE IF( WANTVR ) THEN
               MINWRK = 4*N
               MAXWRK = MAX( MAXWRK, 2*N + ( N - 1 )*ILAENV( 1, &
                             'DORGHR', ' ', N, 1, N, -1 ) )
               CALL DHSEQR( 'S', 'V', N, 1, N, A, LDA, WR, WI, VR, LDVR, &
                      WORK, -1, INFO )
               HSWORK = WORK( 1 )
               MAXWRK = MAX( MAXWRK, N + 1, N + HSWORK )
               MAXWRK = MAX( MAXWRK, 4*N )
            ELSE
               MINWRK = 3*N
               CALL DHSEQR( 'E', 'N', N, 1, N, A, LDA, WR, WI, VR, LDVR, &
                      WORK, -1, INFO )
               HSWORK = WORK( 1 )
               MAXWRK = MAX( MAXWRK, N + 1, N + HSWORK )
            END IF
            MAXWRK = MAX( MAXWRK, MINWRK )
         END IF
         WORK( 1 ) = MAXWRK
!
         IF( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) THEN
            INFO = -13
         END IF
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEEV ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
!
!     Get machine constants
!
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM
!
!     Scale A if max element outside range [SMLNUM,BIGNUM]
!
      ANRM = DLANGE( 'M', N, N, A, LDA, DUM )
      SCALEA = .FALSE.
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         SCALEA = .TRUE.
         CSCALE = SMLNUM
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         SCALEA = .TRUE.
         CSCALE = BIGNUM
      END IF
      IF( SCALEA ) &
         CALL DLASCL( 'G', 0, 0, ANRM, CSCALE, N, N, A, LDA, IERR )
!
!     Balance the matrix
!     (Workspace: need N)
!
      IBAL = 1
      CALL DGEBAL( 'B', N, A, LDA, ILO, IHI, WORK( IBAL ), IERR )
!
!     Reduce to upper Hessenberg form
!     (Workspace: need 3*N, prefer 2*N+N*NB)
!
      ITAU = IBAL + N
      IWRK = ITAU + N
      CALL DGEHRD( N, ILO, IHI, A, LDA, WORK( ITAU ), WORK( IWRK ), &
                   LWORK-IWRK+1, IERR )
!
      IF( WANTVL ) THEN
!
!        Want left eigenvectors
!        Copy Householder vectors to VL
!
         SIDE = 'L'
         CALL DLACPY( 'L', N, N, A, LDA, VL, LDVL )
!
!        Generate orthogonal matrix in VL
!        (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)
!
         CALL DORGHR( N, ILO, IHI, VL, LDVL, WORK( ITAU ), WORK( IWRK ), &
                      LWORK-IWRK+1, IERR )
!
!        Perform QR iteration, accumulating Schur vectors in VL
!        (Workspace: need N+1, prefer N+HSWORK (see comments) )
!
         IWRK = ITAU
         CALL DHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, WR, WI, VL, LDVL, &
                      WORK( IWRK ), LWORK-IWRK+1, INFO )
!
         IF( WANTVR ) THEN
!
!           Want left and right eigenvectors
!           Copy Schur vectors to VR
!
            SIDE = 'B'
            CALL DLACPY( 'F', N, N, VL, LDVL, VR, LDVR )
         END IF
!
      ELSE IF( WANTVR ) THEN
!
!        Want right eigenvectors
!        Copy Householder vectors to VR
!
         SIDE = 'R'
         CALL DLACPY( 'L', N, N, A, LDA, VR, LDVR )
!
!        Generate orthogonal matrix in VR
!        (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)
!
         CALL DORGHR( N, ILO, IHI, VR, LDVR, WORK( ITAU ), WORK( IWRK ), &
                      LWORK-IWRK+1, IERR )
!
!        Perform QR iteration, accumulating Schur vectors in VR
!        (Workspace: need N+1, prefer N+HSWORK (see comments) )
!
         IWRK = ITAU
         CALL DHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, WR, WI, VR, LDVR, &
                      WORK( IWRK ), LWORK-IWRK+1, INFO )
!
      ELSE
!
!        Compute eigenvalues only
!        (Workspace: need N+1, prefer N+HSWORK (see comments) )
!
         IWRK = ITAU
         CALL DHSEQR( 'E', 'N', N, ILO, IHI, A, LDA, WR, WI, VR, LDVR, &
                      WORK( IWRK ), LWORK-IWRK+1, INFO )
      END IF
!
!     If INFO > 0 from DHSEQR, then quit
!
      IF( INFO.GT.0 ) &
         GO TO 50
!
      IF( WANTVL .OR. WANTVR ) THEN
!
!        Compute left and/or right eigenvectors
!        (Workspace: need 4*N)
!
         CALL DTREVC( SIDE, 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, &
                      N, NOUT, WORK( IWRK ), IERR )
      END IF
!
      IF( WANTVL ) THEN
!
!        Undo balancing of left eigenvectors
!        (Workspace: need N)
!
         CALL DGEBAK( 'B', 'L', N, ILO, IHI, WORK( IBAL ), N, VL, LDVL, &
                      IERR )
!
!        Normalize left eigenvectors and make largest component real
!
         DO 20 I = 1, N
            IF( WI( I ).EQ.ZERO ) THEN
               SCL = ONE / DNRM2( N, VL( 1, I ), 1 )
               CALL DSCAL( N, SCL, VL( 1, I ), 1 )
            ELSE IF( WI( I ).GT.ZERO ) THEN
               SCL = ONE / DLAPY2( DNRM2( N, VL( 1, I ), 1 ), &
                     DNRM2( N, VL( 1, I+1 ), 1 ) )
               CALL DSCAL( N, SCL, VL( 1, I ), 1 )
               CALL DSCAL( N, SCL, VL( 1, I+1 ), 1 )
               DO 10 K = 1, N
                  WORK( IWRK+K-1 ) = VL( K, I )**2 + VL( K, I+1 )**2
   10          CONTINUE
               K = IDAMAX( N, WORK( IWRK ), 1 )
               CALL DLARTG( VL( K, I ), VL( K, I+1 ), CS, SN, R )
               CALL DROT( N, VL( 1, I ), 1, VL( 1, I+1 ), 1, CS, SN )
               VL( K, I+1 ) = ZERO
            END IF
   20    CONTINUE
      END IF
!
      IF( WANTVR ) THEN
!
!        Undo balancing of right eigenvectors
!        (Workspace: need N)
!
         CALL DGEBAK( 'B', 'R', N, ILO, IHI, WORK( IBAL ), N, VR, LDVR, &
                      IERR )
!
!        Normalize right eigenvectors and make largest component real
!
         DO 40 I = 1, N
            IF( WI( I ).EQ.ZERO ) THEN
               SCL = ONE / DNRM2( N, VR( 1, I ), 1 )
               CALL DSCAL( N, SCL, VR( 1, I ), 1 )
            ELSE IF( WI( I ).GT.ZERO ) THEN
               SCL = ONE / DLAPY2( DNRM2( N, VR( 1, I ), 1 ), &
                     DNRM2( N, VR( 1, I+1 ), 1 ) )
               CALL DSCAL( N, SCL, VR( 1, I ), 1 )
               CALL DSCAL( N, SCL, VR( 1, I+1 ), 1 )
               DO 30 K = 1, N
                  WORK( IWRK+K-1 ) = VR( K, I )**2 + VR( K, I+1 )**2
   30          CONTINUE
               K = IDAMAX( N, WORK( IWRK ), 1 )
               CALL DLARTG( VR( K, I ), VR( K, I+1 ), CS, SN, R )
               CALL DROT( N, VR( 1, I ), 1, VR( 1, I+1 ), 1, CS, SN )
               VR( K, I+1 ) = ZERO
            END IF
   40    CONTINUE
      END IF
!
!     Undo scaling if necessary
!
   50 CONTINUE
      IF( SCALEA ) THEN
         CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, N-INFO, 1, WR( INFO+1 ), &
                      MAX( N-INFO, 1 ), IERR )
         CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, N-INFO, 1, WI( INFO+1 ), &
                      MAX( N-INFO, 1 ), IERR )
         IF( INFO.GT.0 ) THEN
            CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, ILO-1, 1, WR, N, &
                         IERR )
            CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, ILO-1, 1, WI, N, &
                         IERR )
         END IF
      END IF
!
      WORK( 1 ) = MAXWRK
      RETURN
!
!     End of DGEEV
!
      END SUBROUTINE DGEEV


      SUBROUTINE DGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, WR, WI, &
                         VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, &
                         RCONDE, RCONDV, WORK, LWORK, IWORK, INFO )
!
!  -- LAPACK driver routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          BALANC, JOBVL, JOBVR, SENSE
      INTEGER            IHI, ILO, INFO, LDA, LDVL, LDVR, LWORK, N
      DOUBLE PRECISION   ABNRM
!     ..
!     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), RCONDE( * ), RCONDV( * ), &
                         SCALE( * ), VL( LDVL, * ), VR( LDVR, * ), &
                         WI( * ), WORK( * ), WR( * )
!     ..
!
!  Purpose
!  =======
!
!  DGEEVX computes for an N-by-N real nonsymmetric matrix A, the
!  eigenvalues and, optionally, the left and/or right eigenvectors.
!
!  Optionally also, it computes a balancing transformation to improve
!  the conditioning of the eigenvalues and eigenvectors (ILO, IHI,
!  SCALE, and ABNRM), reciprocal condition numbers for the eigenvalues
!  (RCONDE), and reciprocal condition numbers for the right
!  eigenvectors (RCONDV).
!
!  The right eigenvector v(j) of A satisfies
!                   A * v(j) = lambda(j) * v(j)
!  where lambda(j) is its eigenvalue.
!  The left eigenvector u(j) of A satisfies
!                u(j)**H * A = lambda(j) * u(j)**H
!  where u(j)**H denotes the conjugate transpose of u(j).
!
!  The computed eigenvectors are normalized to have Euclidean norm
!  equal to 1 and largest component real.
!
!  Balancing a matrix means permuting the rows and columns to make it
!  more nearly upper triangular, and applying a diagonal similarity
!  transformation D * A * D**(-1), where D is a diagonal matrix, to
!  make its rows and columns closer in norm and the condition numbers
!  of its eigenvalues and eigenvectors smaller.  The computed
!  reciprocal condition numbers correspond to the balanced matrix.
!  Permuting rows and columns will not change the condition numbers
!  (in exact arithmetic) but diagonal scaling will.  For further
!  explanation of balancing, see section 4.10.2 of the LAPACK
!  Users' Guide.
!
!  Arguments
!  =========
!
!  BALANC  (input) CHARACTER*1
!          Indicates how the input matrix should be diagonally scaled
!          and/or permuted to improve the conditioning of its
!          eigenvalues.
!          = 'N': Do not diagonally scale or permute;
!          = 'P': Perform permutations to make the matrix more nearly
!                 upper triangular. Do not diagonally scale;
!          = 'S': Diagonally scale the matrix, i.e. replace A by
!                 D*A*D**(-1), where D is a diagonal matrix chosen
!                 to make the rows and columns of A more equal in
!                 norm. Do not permute;
!          = 'B': Both diagonally scale and permute A.
!
!          Computed reciprocal condition numbers will be for the matrix
!          after balancing and/or permuting. Permuting does not change
!          condition numbers (in exact arithmetic), but balancing does.
!
!  JOBVL   (input) CHARACTER*1
!          = 'N': left eigenvectors of A are not computed;
!          = 'V': left eigenvectors of A are computed.
!          If SENSE = 'E' or 'B', JOBVL must = 'V'.
!
!  JOBVR   (input) CHARACTER*1
!          = 'N': right eigenvectors of A are not computed;
!          = 'V': right eigenvectors of A are computed.
!          If SENSE = 'E' or 'B', JOBVR must = 'V'.
!
!  SENSE   (input) CHARACTER*1
!          Determines which reciprocal condition numbers are computed.
!          = 'N': None are computed;
!          = 'E': Computed for eigenvalues only;
!          = 'V': Computed for right eigenvectors only;
!          = 'B': Computed for eigenvalues and right eigenvectors.
!
!          If SENSE = 'E' or 'B', both left and right eigenvectors
!          must also be computed (JOBVL = 'V' and JOBVR = 'V').
!
!  N       (input) INTEGER
!          The order of the matrix A. N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the N-by-N matrix A.
!          On exit, A has been overwritten.  If JOBVL = 'V' or
!          JOBVR = 'V', A contains the real Schur form of the balanced
!          version of the input matrix A.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  WR      (output) DOUBLE PRECISION array, dimension (N)
!  WI      (output) DOUBLE PRECISION array, dimension (N)
!          WR and WI contain the real and imaginary parts,
!          respectively, of the computed eigenvalues.  Complex
!          conjugate pairs of eigenvalues will appear consecutively
!          with the eigenvalue having the positive imaginary part
!          first.
!
!  VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)
!          If JOBVL = 'V', the left eigenvectors u(j) are stored one
!          after another in the columns of VL, in the same order
!          as their eigenvalues.
!          If JOBVL = 'N', VL is not referenced.
!          If the j-th eigenvalue is real, then u(j) = VL(:,j),
!          the j-th column of VL.
!          If the j-th and (j+1)-st eigenvalues form a complex
!          conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and
!          u(j+1) = VL(:,j) - i*VL(:,j+1).
!
!  LDVL    (input) INTEGER
!          The leading dimension of the array VL.  LDVL >= 1; if
!          JOBVL = 'V', LDVL >= N.
!
!  VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)
!          If JOBVR = 'V', the right eigenvectors v(j) are stored one
!          after another in the columns of VR, in the same order
!          as their eigenvalues.
!          If JOBVR = 'N', VR is not referenced.
!          If the j-th eigenvalue is real, then v(j) = VR(:,j),
!          the j-th column of VR.
!          If the j-th and (j+1)-st eigenvalues form a complex
!          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
!          v(j+1) = VR(:,j) - i*VR(:,j+1).
!
!  LDVR    (input) INTEGER
!          The leading dimension of the array VR.  LDVR >= 1, and if
!          JOBVR = 'V', LDVR >= N.
!
!  ILO     (output) INTEGER
!  IHI     (output) INTEGER
!          ILO and IHI are integer values determined when A was
!          balanced.  The balanced A(i,j) = 0 if I > J and
!          J = 1,...,ILO-1 or I = IHI+1,...,N.
!
!  SCALE   (output) DOUBLE PRECISION array, dimension (N)
!          Details of the permutations and scaling factors applied
!          when balancing A.  If P(j) is the index of the row and column
!          interchanged with row and column j, and D(j) is the scaling
!          factor applied to row and column j, then
!          SCALE(J) = P(J),    for J = 1,...,ILO-1
!                   = D(J),    for J = ILO,...,IHI
!                   = P(J)     for J = IHI+1,...,N.
!          The order in which the interchanges are made is N to IHI+1,
!          then 1 to ILO-1.
!
!  ABNRM   (output) DOUBLE PRECISION
!          The one-norm of the balanced matrix (the maximum
!          of the sum of absolute values of elements of any column).
!
!  RCONDE  (output) DOUBLE PRECISION array, dimension (N)
!          RCONDE(j) is the reciprocal condition number of the j-th
!          eigenvalue.
!
!  RCONDV  (output) DOUBLE PRECISION array, dimension (N)
!          RCONDV(j) is the reciprocal condition number of the j-th
!          right eigenvector.
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,L
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.   If SENSE = 'N' or 'E',
!          LWORK >= max(1,2*N), and if JOBVL = 'V' or JOBVR = 'V',
!          LWORK >= 3*N.  If SENSE = 'V' or 'B', LWORK >= N*(N+6).
!          For good performance, LWORK must generally be larger.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  IWORK   (workspace) INTEGER array, dimension (2*N-2)
!          If SENSE = 'N' or 'E', not referenced.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  if INFO = i, the QR algorithm failed to compute all the
!                eigenvalues, and no eigenvectors or condition numbers
!                have been computed; elements 1:ILO-1 and i+1:N of WR
!                and WI contain eigenvalues which have converged.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY, SCALEA, WANTVL, WANTVR, WNTSNB, WNTSNE, &
                         WNTSNN, WNTSNV
      CHARACTER          JOB, SIDE
      INTEGER            HSWORK, I, ICOND, IERR, ITAU, IWRK, K, MAXWRK, &
                         MINWRK, NOUT
      DOUBLE PRECISION   ANRM, BIGNUM, CS, CSCALE, EPS, R, SCL, SMLNUM, &
                         SN
!     ..
!     .. Local Arrays ..
      LOGICAL            SELECT( 1 )
      DOUBLE PRECISION   DUM( 1 )
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DGEBAK, DGEBAL, DGEHRD, DHSEQR, DLABAD, DLACPY, &
!                         DLARTG, DLASCL, DORGHR, DROT, DSCAL, DTREVC, &
!                         DTRSNA, XERBLA
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      INTEGER            IDAMAX, ILAENV
!      DOUBLE PRECISION   DLAMCH, DLANGE, DLAPY2, DNRM2
!      EXTERNAL           LSAME, IDAMAX, ILAENV, DLAMCH, DLANGE, DLAPY2, &
!                         DNRM2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      WANTVL = LSAME( JOBVL, 'V' )
      WANTVR = LSAME( JOBVR, 'V' )
      WNTSNN = LSAME( SENSE, 'N' )
      WNTSNE = LSAME( SENSE, 'E' )
      WNTSNV = LSAME( SENSE, 'V' )
      WNTSNB = LSAME( SENSE, 'B' )
      IF( .NOT.( LSAME( BALANC, 'N' ) .OR. LSAME( BALANC, &
          'S' ) .OR. LSAME( BALANC, 'P' ) .OR. LSAME( BALANC, 'B' ) ) ) &
           THEN
         INFO = -1
      ELSE IF( ( .NOT.WANTVL ) .AND. ( .NOT.LSAME( JOBVL, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( ( .NOT.WANTVR ) .AND. ( .NOT.LSAME( JOBVR, 'N' ) ) ) THEN
         INFO = -3
      ELSE IF( .NOT.( WNTSNN .OR. WNTSNE .OR. WNTSNB .OR. WNTSNV ) .OR. &
               ( ( WNTSNE .OR. WNTSNB ) .AND. .NOT.( WANTVL .AND. &
               WANTVR ) ) ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDVL.LT.1 .OR. ( WANTVL .AND. LDVL.LT.N ) ) THEN
         INFO = -11
      ELSE IF( LDVR.LT.1 .OR. ( WANTVR .AND. LDVR.LT.N ) ) THEN
         INFO = -13
      END IF
!
!     Compute workspace
!      (Note: Comments in the code beginning "Workspace:" describe the
!       minimal amount of workspace needed at that point in the code,
!       as well as the preferred amount for good performance.
!       NB refers to the optimal block size for the immediately
!       following subroutine, as returned by ILAENV.
!       HSWORK refers to the workspace preferred by DHSEQR, as
!       calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
!       the worst case.)
!
      IF( INFO.EQ.0 ) THEN
         IF( N.EQ.0 ) THEN
            MINWRK = 1
            MAXWRK = 1
         ELSE
            MAXWRK = N + N*ILAENV( 1, 'DGEHRD', ' ', N, 1, N, 0 )
!
            IF( WANTVL ) THEN
               CALL DHSEQR( 'S', 'V', N, 1, N, A, LDA, WR, WI, VL, LDVL, &
                      WORK, -1, INFO )
            ELSE IF( WANTVR ) THEN
               CALL DHSEQR( 'S', 'V', N, 1, N, A, LDA, WR, WI, VR, LDVR, &
                      WORK, -1, INFO )
            ELSE
               IF( WNTSNN ) THEN
                  CALL DHSEQR( 'E', 'N', N, 1, N, A, LDA, WR, WI, VR, &
                      LDVR, WORK, -1, INFO )
               ELSE
                  CALL DHSEQR( 'S', 'N', N, 1, N, A, LDA, WR, WI, VR, &
                      LDVR, WORK, -1, INFO )
               END IF
            END IF
            HSWORK = WORK( 1 )
!
            IF( ( .NOT.WANTVL ) .AND. ( .NOT.WANTVR ) ) THEN
               MINWRK = 2*N
               IF( .NOT.WNTSNN ) &
                  MINWRK = MAX( MINWRK, N*N+6*N )
               MAXWRK = MAX( MAXWRK, HSWORK )
               IF( .NOT.WNTSNN ) &
                  MAXWRK = MAX( MAXWRK, N*N + 6*N )
            ELSE
               MINWRK = 3*N
               IF( ( .NOT.WNTSNN ) .AND. ( .NOT.WNTSNE ) ) &
                  MINWRK = MAX( MINWRK, N*N + 6*N )
               MAXWRK = MAX( MAXWRK, HSWORK )
               MAXWRK = MAX( MAXWRK, N + ( N - 1 )*ILAENV( 1, 'DORGHR', &
                             ' ', N, 1, N, -1 ) )
               IF( ( .NOT.WNTSNN ) .AND. ( .NOT.WNTSNE ) ) &
                  MAXWRK = MAX( MAXWRK, N*N + 6*N )
               MAXWRK = MAX( MAXWRK, 3*N )
            END IF
            MAXWRK = MAX( MAXWRK, MINWRK )
         END IF
         WORK( 1 ) = MAXWRK
!
         IF( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) THEN
            INFO = -21
         END IF
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEEVX', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
!
!     Get machine constants
!
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM
!
!     Scale A if max element outside range [SMLNUM,BIGNUM]
!
      ICOND = 0
      ANRM = DLANGE( 'M', N, N, A, LDA, DUM )
      SCALEA = .FALSE.
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         SCALEA = .TRUE.
         CSCALE = SMLNUM
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         SCALEA = .TRUE.
         CSCALE = BIGNUM
      END IF
      IF( SCALEA ) &
         CALL DLASCL( 'G', 0, 0, ANRM, CSCALE, N, N, A, LDA, IERR )
!
!     Balance the matrix and compute ABNRM
!
      CALL DGEBAL( BALANC, N, A, LDA, ILO, IHI, SCALE, IERR )
      ABNRM = DLANGE( '1', N, N, A, LDA, DUM )
      IF( SCALEA ) THEN
         DUM( 1 ) = ABNRM
         CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, 1, 1, DUM, 1, IERR )
         ABNRM = DUM( 1 )
      END IF
!
!     Reduce to upper Hessenberg form
!     (Workspace: need 2*N, prefer N+N*NB)
!
      ITAU = 1
      IWRK = ITAU + N
      CALL DGEHRD( N, ILO, IHI, A, LDA, WORK( ITAU ), WORK( IWRK ), &
                   LWORK-IWRK+1, IERR )
!
      IF( WANTVL ) THEN
!
!        Want left eigenvectors
!        Copy Householder vectors to VL
!
         SIDE = 'L'
         CALL DLACPY( 'L', N, N, A, LDA, VL, LDVL )
!
!        Generate orthogonal matrix in VL
!        (Workspace: need 2*N-1, prefer N+(N-1)*NB)
!
         CALL DORGHR( N, ILO, IHI, VL, LDVL, WORK( ITAU ), WORK( IWRK ), &
                      LWORK-IWRK+1, IERR )
!
!        Perform QR iteration, accumulating Schur vectors in VL
!        (Workspace: need 1, prefer HSWORK (see comments) )
!
         IWRK = ITAU
         CALL DHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, WR, WI, VL, LDVL, &
                      WORK( IWRK ), LWORK-IWRK+1, INFO )
!
         IF( WANTVR ) THEN
!
!           Want left and right eigenvectors
!           Copy Schur vectors to VR
!
            SIDE = 'B'
            CALL DLACPY( 'F', N, N, VL, LDVL, VR, LDVR )
         END IF
!
      ELSE IF( WANTVR ) THEN
!
!        Want right eigenvectors
!        Copy Householder vectors to VR
!
         SIDE = 'R'
         CALL DLACPY( 'L', N, N, A, LDA, VR, LDVR )
!
!        Generate orthogonal matrix in VR
!        (Workspace: need 2*N-1, prefer N+(N-1)*NB)
!
         CALL DORGHR( N, ILO, IHI, VR, LDVR, WORK( ITAU ), WORK( IWRK ), &
                      LWORK-IWRK+1, IERR )
!
!        Perform QR iteration, accumulating Schur vectors in VR
!        (Workspace: need 1, prefer HSWORK (see comments) )
!
         IWRK = ITAU
         CALL DHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, WR, WI, VR, LDVR, &
                      WORK( IWRK ), LWORK-IWRK+1, INFO )
!
      ELSE
!
!        Compute eigenvalues only
!        If condition numbers desired, compute Schur form
!
         IF( WNTSNN ) THEN
            JOB = 'E'
         ELSE
            JOB = 'S'
         END IF
!
!        (Workspace: need 1, prefer HSWORK (see comments) )
!
         IWRK = ITAU
         CALL DHSEQR( JOB, 'N', N, ILO, IHI, A, LDA, WR, WI, VR, LDVR, &
                      WORK( IWRK ), LWORK-IWRK+1, INFO )
      END IF
!
!     If INFO > 0 from DHSEQR, then quit
!
      IF( INFO.GT.0 ) &
         GO TO 50
!
      IF( WANTVL .OR. WANTVR ) THEN
!
!        Compute left and/or right eigenvectors
!        (Workspace: need 3*N)
!
         CALL DTREVC( SIDE, 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, &
                      N, NOUT, WORK( IWRK ), IERR )
      END IF
!
!     Compute condition numbers if desired
!     (Workspace: need N*N+6*N unless SENSE = 'E')
!
      IF( .NOT.WNTSNN ) THEN
         CALL DTRSNA( SENSE, 'A', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, &
                      RCONDE, RCONDV, N, NOUT, WORK( IWRK ), N, IWORK, &
                      ICOND )
      END IF
!
      IF( WANTVL ) THEN
!
!        Undo balancing of left eigenvectors
!
         CALL DGEBAK( BALANC, 'L', N, ILO, IHI, SCALE, N, VL, LDVL, &
                      IERR )
!
!        Normalize left eigenvectors and make largest component real
!
         DO 20 I = 1, N
            IF( WI( I ).EQ.ZERO ) THEN
               SCL = ONE / DNRM2( N, VL( 1, I ), 1 )
               CALL DSCAL( N, SCL, VL( 1, I ), 1 )
            ELSE IF( WI( I ).GT.ZERO ) THEN
               SCL = ONE / DLAPY2( DNRM2( N, VL( 1, I ), 1 ), &
                     DNRM2( N, VL( 1, I+1 ), 1 ) )
               CALL DSCAL( N, SCL, VL( 1, I ), 1 )
               CALL DSCAL( N, SCL, VL( 1, I+1 ), 1 )
               DO 10 K = 1, N
                  WORK( K ) = VL( K, I )**2 + VL( K, I+1 )**2
   10          CONTINUE
               K = IDAMAX( N, WORK, 1 )
               CALL DLARTG( VL( K, I ), VL( K, I+1 ), CS, SN, R )
               CALL DROT( N, VL( 1, I ), 1, VL( 1, I+1 ), 1, CS, SN )
               VL( K, I+1 ) = ZERO
            END IF
   20    CONTINUE
      END IF
!
      IF( WANTVR ) THEN
!
!        Undo balancing of right eigenvectors
!
         CALL DGEBAK( BALANC, 'R', N, ILO, IHI, SCALE, N, VR, LDVR, &
                      IERR )
!
!        Normalize right eigenvectors and make largest component real
!
         DO 40 I = 1, N
            IF( WI( I ).EQ.ZERO ) THEN
               SCL = ONE / DNRM2( N, VR( 1, I ), 1 )
               CALL DSCAL( N, SCL, VR( 1, I ), 1 )
            ELSE IF( WI( I ).GT.ZERO ) THEN
               SCL = ONE / DLAPY2( DNRM2( N, VR( 1, I ), 1 ), &
                     DNRM2( N, VR( 1, I+1 ), 1 ) )
               CALL DSCAL( N, SCL, VR( 1, I ), 1 )
               CALL DSCAL( N, SCL, VR( 1, I+1 ), 1 )
               DO 30 K = 1, N
                  WORK( K ) = VR( K, I )**2 + VR( K, I+1 )**2
   30          CONTINUE
               K = IDAMAX( N, WORK, 1 )
               CALL DLARTG( VR( K, I ), VR( K, I+1 ), CS, SN, R )
               CALL DROT( N, VR( 1, I ), 1, VR( 1, I+1 ), 1, CS, SN )
               VR( K, I+1 ) = ZERO
            END IF
   40    CONTINUE
      END IF
!
!     Undo scaling if necessary
!
   50 CONTINUE
      IF( SCALEA ) THEN
         CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, N-INFO, 1, WR( INFO+1 ), &
                      MAX( N-INFO, 1 ), IERR )
         CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, N-INFO, 1, WI( INFO+1 ), &
                      MAX( N-INFO, 1 ), IERR )
         IF( INFO.EQ.0 ) THEN
            IF( ( WNTSNV .OR. WNTSNB ) .AND. ICOND.EQ.0 ) &
               CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, N, 1, RCONDV, N, &
                            IERR )
         ELSE
            CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, ILO-1, 1, WR, N, &
                         IERR )
            CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, ILO-1, 1, WI, N, &
                         IERR )
         END IF
      END IF
!
      WORK( 1 ) = MAXWRK
      RETURN
!
!     End of DGEEVX
!
      END SUBROUTINE DGEEVX


      SUBROUTINE DGEHRD( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION  A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DGEHRD reduces a real general matrix A to upper Hessenberg form H by
!  an orthogonal similarity transformation:  Q' * A * Q = H .
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  ILO     (input) INTEGER
!  IHI     (input) INTEGER
!          It is assumed that A is already upper triangular in rows
!          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
!          set by a previous call to DGEBAL; otherwise they should be
!          set to 1 and N respectively. See Further Details.
!          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the N-by-N general matrix to be reduced.
!          On exit, the upper triangle and the first subdiagonal of A
!          are overwritten with the upper Hessenberg matrix H, and the
!          elements below the first subdiagonal, with the array TAU,
!          represent the orthogonal matrix Q as a product of elementary
!          reflectors. See Further Details.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  TAU     (output) DOUBLE PRECISION array, dimension (N-1)
!          The scalar factors of the elementary reflectors (see Further
!          Details). Elements 1:ILO-1 and IHI:N-1 of TAU are set to
!          zero.
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The length of the array WORK.  LWORK >= max(1,N).
!          For optimum performance LWORK >= N*NB, where NB is the
!          optimal blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  The matrix Q is represented as a product of (ihi-ilo) elementary
!  reflectors
!
!     Q = H(ilo) H(ilo+1) . . . H(ihi-1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on
!  exit in A(i+2:ihi,i), and tau in TAU(i).
!
!  The contents of A are illustrated by the following example, with
!  n = 7, ilo = 2 and ihi = 6:
!
!  on entry,                        on exit,
!
!  ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )
!  (     a   a   a   a   a   a )    (      a   h   h   h   h   a )
!  (     a   a   a   a   a   a )    (      h   h   h   h   h   h )
!  (     a   a   a   a   a   a )    (      v2  h   h   h   h   h )
!  (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )
!  (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )
!  (                         a )    (                          a )
!
!  where a denotes an element of the original matrix A, h denotes a
!  modified element of the upper Hessenberg matrix H, and vi denotes an
!  element of the vector defining H(i).
!
!  This file is a slight modification of LAPACK-3.0's DGEHRD
!  subroutine incorporating improvements proposed by Quintana-Orti and
!  Van de Geijn (2005).
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            NBMAX, LDT
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1 )
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, &
                           ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, LDWORK, LWKOPT, NB, &
                         NBMIN, NH, NX
      DOUBLE PRECISION  EI
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION  T( LDT, NBMAX )
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DAXPY, DGEHD2, DGEMM, DLAHR2, DLARFB, DTRMM, &
!                         XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. External Functions ..
!      INTEGER            ILAENV
!      EXTERNAL           ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0
      NB = MIN( NBMAX, ILAENV( 1, 'DGEHRD', ' ', N, ILO, IHI, -1 ) )
      LWKOPT = N*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEHRD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Set elements 1:ILO-1 and IHI:N-1 of TAU to zero
!
      DO 10 I = 1, ILO - 1
         TAU( I ) = ZERO
   10 CONTINUE
      DO 20 I = MAX( 1, IHI ), N - 1
         TAU( I ) = ZERO
   20 CONTINUE
!
!     Quick return if possible
!
      NH = IHI - ILO + 1
      IF( NH.LE.1 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
!     Determine the block size
!
      NB = MIN( NBMAX, ILAENV( 1, 'DGEHRD', ' ', N, ILO, IHI, -1 ) )
      NBMIN = 2
      IWS = 1
      IF( NB.GT.1 .AND. NB.LT.NH ) THEN
!
!        Determine when to cross over from blocked to unblocked code
!        (last block is always handled by unblocked code)
!
         NX = MAX( NB, ILAENV( 3, 'DGEHRD', ' ', N, ILO, IHI, -1 ) )
         IF( NX.LT.NH ) THEN
!
!           Determine if workspace is large enough for blocked code
!
            IWS = N*NB
            IF( LWORK.LT.IWS ) THEN
!
!              Not enough workspace to use optimal NB:  determine the
!              minimum value of NB, and reduce NB or force use of
!              unblocked code
!
               NBMIN = MAX( 2, ILAENV( 2, 'DGEHRD', ' ', N, ILO, IHI, &
                       -1 ) )
               IF( LWORK.GE.N*NBMIN ) THEN
                  NB = LWORK / N
               ELSE
                  NB = 1
               END IF
            END IF
         END IF
      END IF
      LDWORK = N
!
      IF( NB.LT.NBMIN .OR. NB.GE.NH ) THEN
!
!        Use unblocked code below
!
         I = ILO
!
      ELSE
!
!        Use blocked code
!
         DO 40 I = ILO, IHI - 1 - NX, NB
            IB = MIN( NB, IHI-I )
!
!           Reduce columns i:i+ib-1 to Hessenberg form, returning the
!           matrices V and T of the block reflector H = I - V*T*V'
!           which performs the reduction, and also the matrix Y = A*V*T
!
            CALL DLAHR2( IHI, I, IB, A( 1, I ), LDA, TAU( I ), T, LDT, &
                         WORK, LDWORK )
!
!           Apply the block reflector H to A(1:ihi,i+ib:ihi) from the
!           right, computing  A := A - Y * V'. V(i+ib,ib-1) must be set
!           to 1
!
            EI = A( I+IB, I+IB-1 )
            A( I+IB, I+IB-1 ) = ONE
            CALL DGEMM( 'No transpose', 'Transpose', &
                        IHI, IHI-I-IB+1, &
                        IB, -ONE, WORK, LDWORK, A( I+IB, I ), LDA, ONE, &
                        A( 1, I+IB ), LDA )
            A( I+IB, I+IB-1 ) = EI
!
!           Apply the block reflector H to A(1:i,i+1:i+ib-1) from the
!           right
!
            CALL DTRMM( 'Right', 'Lower', 'Transpose', &
                        'Unit', I, IB-1, &
                        ONE, A( I+1, I ), LDA, WORK, LDWORK )
            DO 30 J = 0, IB-2
               CALL DAXPY( I, -ONE, WORK( LDWORK*J+1 ), 1, &
                           A( 1, I+J+1 ), 1 )
   30       CONTINUE
!
!           Apply the block reflector H to A(i+1:ihi,i+ib:n) from the
!           left
!
            CALL DLARFB( 'Left', 'Transpose', 'Forward', &
                         'Columnwise', &
                         IHI-I, N-I-IB+1, IB, A( I+1, I ), LDA, T, LDT, &
                         A( I+1, I+IB ), LDA, WORK, LDWORK )
   40    CONTINUE
      END IF
!
!     Use unblocked code to reduce the rest of the matrix
!
      CALL DGEHD2( N, I, IHI, A, LDA, TAU, WORK, IINFO )
      WORK( 1 ) = IWS
!
      RETURN
!
!     End of DGEHRD
!
      END SUBROUTINE DGEHRD


      SUBROUTINE DTREVC( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, &
                         LDVR, MM, M, WORK, INFO )
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          HOWMNY, SIDE
      INTEGER            INFO, LDT, LDVL, LDVR, M, MM, N
!     ..
!     .. Array Arguments ..
      LOGICAL            SELECT( * )
      DOUBLE PRECISION   T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ), &
                         WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DTREVC computes some or all of the right and/or left eigenvectors of
!  a real upper quasi-triangular matrix T.
!  Matrices of this type are produced by the Schur factorization of
!  a real general matrix:  A = Q*T*Q**T, as computed by DHSEQR.
!
!  The right eigenvector x and the left eigenvector y of T corresponding
!  to an eigenvalue w are defined by:
!
!     T*x = w*x,     (y**H)*T = w*(y**H)
!
!  where y**H denotes the conjugate transpose of y.
!  The eigenvalues are not input to this routine, but are read directly
!  from the diagonal blocks of T.
!
!  This routine returns the matrices X and/or Y of right and left
!  eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an
!  input matrix.  If Q is the orthogonal factor that reduces a matrix
!  A to Schur form T, then Q*X and Q*Y are the matrices of right and
!  left eigenvectors of A.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'R':  compute right eigenvectors only;
!          = 'L':  compute left eigenvectors only;
!          = 'B':  compute both right and left eigenvectors.
!
!  HOWMNY  (input) CHARACTER*1
!          = 'A':  compute all right and/or left eigenvectors;
!          = 'B':  compute all right and/or left eigenvectors,
!                  backtransformed by the matrices in VR and/or VL;
!          = 'S':  compute selected right and/or left eigenvectors,
!                  as indicated by the logical array SELECT.
!
!  SELECT  (input/output) LOGICAL array, dimension (N)
!          If HOWMNY = 'S', SELECT specifies the eigenvectors to be
!          computed.
!          If w(j) is a real eigenvalue, the corresponding real
!          eigenvector is computed if SELECT(j) is .TRUE..
!          If w(j) and w(j+1) are the real and imaginary parts of a
!          complex eigenvalue, the corresponding complex eigenvector is
!          computed if either SELECT(j) or SELECT(j+1) is .TRUE., and
!          on exit SELECT(j) is set to .TRUE. and SELECT(j+1) is set to
!          .FALSE..
!          Not referenced if HOWMNY = 'A' or 'B'.
!
!  N       (input) INTEGER
!          The order of the matrix T. N >= 0.
!
!  T       (input) DOUBLE PRECISION array, dimension (LDT,N)
!          The upper quasi-triangular matrix T in Schur canonical form.
!
!  LDT     (input) INTEGER
!          The leading dimension of the array T. LDT >= max(1,N).
!
!  VL      (input/output) DOUBLE PRECISION array, dimension (LDVL,MM)
!          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must
!          contain an N-by-N matrix Q (usually the orthogonal matrix Q
!          of Schur vectors returned by DHSEQR).
!          On exit, if SIDE = 'L' or 'B', VL contains:
!          if HOWMNY = 'A', the matrix Y of left eigenvectors of T;
!          if HOWMNY = 'B', the matrix Q*Y;
!          if HOWMNY = 'S', the left eigenvectors of T specified by
!                           SELECT, stored consecutively in the columns
!                           of VL, in the same order as their
!                           eigenvalues.
!          A complex eigenvector corresponding to a complex eigenvalue
!          is stored in two consecutive columns, the first holding the
!          real part, and the second the imaginary part.
!          Not referenced if SIDE = 'R'.
!
!  LDVL    (input) INTEGER
!          The leading dimension of the array VL.  LDVL >= 1, and if
!          SIDE = 'L' or 'B', LDVL >= N.
!
!  VR      (input/output) DOUBLE PRECISION array, dimension (LDVR,MM)
!          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must
!          contain an N-by-N matrix Q (usually the orthogonal matrix Q
!          of Schur vectors returned by DHSEQR).
!          On exit, if SIDE = 'R' or 'B', VR contains:
!          if HOWMNY = 'A', the matrix X of right eigenvectors of T;
!          if HOWMNY = 'B', the matrix Q*X;
!          if HOWMNY = 'S', the right eigenvectors of T specified by
!                           SELECT, stored consecutively in the columns
!                           of VR, in the same order as their
!                           eigenvalues.
!          A complex eigenvector corresponding to a complex eigenvalue
!          is stored in two consecutive columns, the first holding the
!          real part and the second the imaginary part.
!          Not referenced if SIDE = 'L'.
!
!  LDVR    (input) INTEGER
!          The leading dimension of the array VR.  LDVR >= 1, and if
!          SIDE = 'R' or 'B', LDVR >= N.
!
!  MM      (input) INTEGER
!          The number of columns in the arrays VL and/or VR. MM >= M.
!
!  M       (output) INTEGER
!          The number of columns in the arrays VL and/or VR actually
!          used to store the eigenvectors.
!          If HOWMNY = 'A' or 'B', M is set to N.
!          Each selected real eigenvector occupies one column and each
!          selected complex eigenvector occupies two columns.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (3*N)
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  Further Details
!  ===============
!
!  The algorithm used in this program is basically backward (forward)
!  substitution, with scaling to make the the code robust against
!  possible overflow.
!
!  Each eigenvector is normalized so that the element of largest
!  magnitude has magnitude 1; here the magnitude of a complex number
!  (x,y) is taken to be |x| + |y|.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            ALLV, BOTHV, LEFTV, OVER, PAIR, RIGHTV, SOMEV
      INTEGER            I, IERR, II, IP, IS, J, J1, J2, JNXT, K, KI, N2
      DOUBLE PRECISION   BETA, BIGNUM, EMAX, OVFL, REC, REMAX, SCALE, &
                         SMIN, SMLNUM, ULP, UNFL, VCRIT, VMAX, WI, WR, &
                         XNORM
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      INTEGER            IDAMAX
!      DOUBLE PRECISION   DDOT, DLAMCH
!      EXTERNAL           LSAME, IDAMAX, DDOT, DLAMCH
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DAXPY, DCOPY, DGEMV, DLALN2, DSCAL, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   X( 2, 2 )
!     ..
!     .. Executable Statements ..
!
!     Decode and test the input parameters
!
      BOTHV = LSAME( SIDE, 'B' )
      RIGHTV = LSAME( SIDE, 'R' ) .OR. BOTHV
      LEFTV = LSAME( SIDE, 'L' ) .OR. BOTHV
!
      ALLV = LSAME( HOWMNY, 'A' )
      OVER = LSAME( HOWMNY, 'B' )
      SOMEV = LSAME( HOWMNY, 'S' )
!
      INFO = 0
      IF( .NOT.RIGHTV .AND. .NOT.LEFTV ) THEN
         INFO = -1
      ELSE IF( .NOT.ALLV .AND. .NOT.OVER .AND. .NOT.SOMEV ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDVL.LT.1 .OR. ( LEFTV .AND. LDVL.LT.N ) ) THEN
         INFO = -8
      ELSE IF( LDVR.LT.1 .OR. ( RIGHTV .AND. LDVR.LT.N ) ) THEN
         INFO = -10
      ELSE
!
!        Set M to the number of columns required to store the selected
!        eigenvectors, standardize the array SELECT if necessary, and
!        test MM.
!
         IF( SOMEV ) THEN
            M = 0
            PAIR = .FALSE.
            DO 10 J = 1, N
               IF( PAIR ) THEN
                  PAIR = .FALSE.
                  SELECT( J ) = .FALSE.
               ELSE
                  IF( J.LT.N ) THEN
                     IF( T( J+1, J ).EQ.ZERO ) THEN
                        IF( SELECT( J ) ) &
                           M = M + 1
                     ELSE
                        PAIR = .TRUE.
                        IF( SELECT( J ) .OR. SELECT( J+1 ) ) THEN
                           SELECT( J ) = .TRUE.
                           M = M + 2
                        END IF
                     END IF
                  ELSE
                     IF( SELECT( N ) ) &
                        M = M + 1
                  END IF
               END IF
   10       CONTINUE
         ELSE
            M = N
         END IF
!
         IF( MM.LT.M ) THEN
            INFO = -11
         END IF
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTREVC', -INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( N.EQ.0 ) &
         RETURN
!
!     Set the constants to control overflow.
!
      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      CALL DLABAD( UNFL, OVFL )
      ULP = DLAMCH( 'Precision' )
      SMLNUM = UNFL*( N / ULP )
      BIGNUM = ( ONE-ULP ) / SMLNUM
!
!     Compute 1-norm of each column of strictly upper triangular
!     part of T to control overflow in triangular solver.
!
      WORK( 1 ) = ZERO
      DO 30 J = 2, N
         WORK( J ) = ZERO
         DO 20 I = 1, J - 1
            WORK( J ) = WORK( J ) + ABS( T( I, J ) )
   20    CONTINUE
   30 CONTINUE
!
!     Index IP is used to specify the real or complex eigenvalue:
!       IP = 0, real eigenvalue,
!            1, first of conjugate complex pair: (wr,wi)
!           -1, second of conjugate complex pair: (wr,wi)
!
      N2 = 2*N
!
      IF( RIGHTV ) THEN
!
!        Compute right eigenvectors.
!
         IP = 0
         IS = M
         DO 140 KI = N, 1, -1
!
            IF( IP.EQ.1 ) &
               GO TO 130
            IF( KI.EQ.1 ) &
               GO TO 40
            IF( T( KI, KI-1 ).EQ.ZERO ) &
               GO TO 40
            IP = -1
!
   40       CONTINUE
            IF( SOMEV ) THEN
               IF( IP.EQ.0 ) THEN
                  IF( .NOT.SELECT( KI ) ) &
                     GO TO 130
               ELSE
                  IF( .NOT.SELECT( KI-1 ) ) &
                     GO TO 130
               END IF
            END IF
!
!           Compute the KI-th eigenvalue (WR,WI).
!
            WR = T( KI, KI )
            WI = ZERO
            IF( IP.NE.0 ) &
               WI = SQRT( ABS( T( KI, KI-1 ) ) )* &
                    SQRT( ABS( T( KI-1, KI ) ) )
            SMIN = MAX( ULP*( ABS( WR )+ABS( WI ) ), SMLNUM )
!
            IF( IP.EQ.0 ) THEN
!
!              Real right eigenvector
!
               WORK( KI+N ) = ONE
!
!              Form right-hand side
!
               DO 50 K = 1, KI - 1
                  WORK( K+N ) = -T( K, KI )
   50          CONTINUE
!
!              Solve the upper quasi-triangular system:
!                 (T(1:KI-1,1:KI-1) - WR)*X = SCALE*WORK.
!
               JNXT = KI - 1
               DO 60 J = KI - 1, 1, -1
                  IF( J.GT.JNXT ) &
                     GO TO 60
                  J1 = J
                  J2 = J
                  JNXT = J - 1
                  IF( J.GT.1 ) THEN
                     IF( T( J, J-1 ).NE.ZERO ) THEN
                        J1 = J - 1
                        JNXT = J - 2
                     END IF
                  END IF
!
                  IF( J1.EQ.J2 ) THEN
!
!                    1-by-1 diagonal block
!
                     CALL DLALN2( .FALSE., 1, 1, SMIN, ONE, T( J, J ), &
                                  LDT, ONE, ONE, WORK( J+N ), N, WR, &
                                  ZERO, X, 2, SCALE, XNORM, IERR )
!
!                    Scale X(1,1) to avoid overflow when updating
!                    the right-hand side.
!
                     IF( XNORM.GT.ONE ) THEN
                        IF( WORK( J ).GT.BIGNUM / XNORM ) THEN
                           X( 1, 1 ) = X( 1, 1 ) / XNORM
                           SCALE = SCALE / XNORM
                        END IF
                     END IF
!
!                    Scale if necessary
!
                     IF( SCALE.NE.ONE ) &
                        CALL DSCAL( KI, SCALE, WORK( 1+N ), 1 )
                     WORK( J+N ) = X( 1, 1 )
!
!                    Update right-hand side
!
                     CALL DAXPY( J-1, -X( 1, 1 ), T( 1, J ), 1, &
                                 WORK( 1+N ), 1 )
!
                  ELSE
!
!                    2-by-2 diagonal block
!
                     CALL DLALN2( .FALSE., 2, 1, SMIN, ONE, &
                                  T( J-1, J-1 ), LDT, ONE, ONE, &
                                  WORK( J-1+N ), N, WR, ZERO, X, 2, &
                                  SCALE, XNORM, IERR )
!
!                    Scale X(1,1) and X(2,1) to avoid overflow when
!                    updating the right-hand side.
!
                     IF( XNORM.GT.ONE ) THEN
                        BETA = MAX( WORK( J-1 ), WORK( J ) )
                        IF( BETA.GT.BIGNUM / XNORM ) THEN
                           X( 1, 1 ) = X( 1, 1 ) / XNORM
                           X( 2, 1 ) = X( 2, 1 ) / XNORM
                           SCALE = SCALE / XNORM
                        END IF
                     END IF
!
!                    Scale if necessary
!
                     IF( SCALE.NE.ONE ) &
                        CALL DSCAL( KI, SCALE, WORK( 1+N ), 1 )
                     WORK( J-1+N ) = X( 1, 1 )
                     WORK( J+N ) = X( 2, 1 )
!
!                    Update right-hand side
!
                     CALL DAXPY( J-2, -X( 1, 1 ), T( 1, J-1 ), 1, &
                                 WORK( 1+N ), 1 )
                     CALL DAXPY( J-2, -X( 2, 1 ), T( 1, J ), 1, &
                                 WORK( 1+N ), 1 )
                  END IF
   60          CONTINUE
!
!              Copy the vector x or Q*x to VR and normalize.
!
               IF( .NOT.OVER ) THEN
                  CALL DCOPY( KI, WORK( 1+N ), 1, VR( 1, IS ), 1 )
!
                  II = IDAMAX( KI, VR( 1, IS ), 1 )
                  REMAX = ONE / ABS( VR( II, IS ) )
                  CALL DSCAL( KI, REMAX, VR( 1, IS ), 1 )
!
                  DO 70 K = KI + 1, N
                     VR( K, IS ) = ZERO
   70             CONTINUE
               ELSE
                  IF( KI.GT.1 ) &
                     CALL DGEMV( 'N', N, KI-1, ONE, VR, LDVR, &
                                 WORK( 1+N ), 1, WORK( KI+N ), &
                                 VR( 1, KI ), 1 )
!
                  II = IDAMAX( N, VR( 1, KI ), 1 )
                  REMAX = ONE / ABS( VR( II, KI ) )
                  CALL DSCAL( N, REMAX, VR( 1, KI ), 1 )
               END IF
!
            ELSE
!
!              Complex right eigenvector.
!
!              Initial solve
!                [ (T(KI-1,KI-1) T(KI-1,KI) ) - (WR + I* WI)]*X = 0.
!                [ (T(KI,KI-1)   T(KI,KI)   )               ]
!
               IF( ABS( T( KI-1, KI ) ).GE.ABS( T( KI, KI-1 ) ) ) THEN
                  WORK( KI-1+N ) = ONE
                  WORK( KI+N2 ) = WI / T( KI-1, KI )
               ELSE
                  WORK( KI-1+N ) = -WI / T( KI, KI-1 )
                  WORK( KI+N2 ) = ONE
               END IF
               WORK( KI+N ) = ZERO
               WORK( KI-1+N2 ) = ZERO
!
!              Form right-hand side
!
               DO 80 K = 1, KI - 2
                  WORK( K+N ) = -WORK( KI-1+N )*T( K, KI-1 )
                  WORK( K+N2 ) = -WORK( KI+N2 )*T( K, KI )
   80          CONTINUE
!
!              Solve upper quasi-triangular system:
!              (T(1:KI-2,1:KI-2) - (WR+i*WI))*X = SCALE*(WORK+i*WORK2)
!
               JNXT = KI - 2
               DO 90 J = KI - 2, 1, -1
                  IF( J.GT.JNXT ) &
                     GO TO 90
                  J1 = J
                  J2 = J
                  JNXT = J - 1
                  IF( J.GT.1 ) THEN
                     IF( T( J, J-1 ).NE.ZERO ) THEN
                        J1 = J - 1
                        JNXT = J - 2
                     END IF
                  END IF
!
                  IF( J1.EQ.J2 ) THEN
!
!                    1-by-1 diagonal block
!
                     CALL DLALN2( .FALSE., 1, 2, SMIN, ONE, T( J, J ), &
                                  LDT, ONE, ONE, WORK( J+N ), N, WR, WI, &
                                  X, 2, SCALE, XNORM, IERR )
!
!                    Scale X(1,1) and X(1,2) to avoid overflow when
!                    updating the right-hand side.
!
                     IF( XNORM.GT.ONE ) THEN
                        IF( WORK( J ).GT.BIGNUM / XNORM ) THEN
                           X( 1, 1 ) = X( 1, 1 ) / XNORM
                           X( 1, 2 ) = X( 1, 2 ) / XNORM
                           SCALE = SCALE / XNORM
                        END IF
                     END IF
!
!                    Scale if necessary
!
                     IF( SCALE.NE.ONE ) THEN
                        CALL DSCAL( KI, SCALE, WORK( 1+N ), 1 )
                        CALL DSCAL( KI, SCALE, WORK( 1+N2 ), 1 )
                     END IF
                     WORK( J+N ) = X( 1, 1 )
                     WORK( J+N2 ) = X( 1, 2 )
!
!                    Update the right-hand side
!
                     CALL DAXPY( J-1, -X( 1, 1 ), T( 1, J ), 1, &
                                 WORK( 1+N ), 1 )
                     CALL DAXPY( J-1, -X( 1, 2 ), T( 1, J ), 1, &
                                 WORK( 1+N2 ), 1 )
!
                  ELSE
!
!                    2-by-2 diagonal block
!
                     CALL DLALN2( .FALSE., 2, 2, SMIN, ONE, &
                                  T( J-1, J-1 ), LDT, ONE, ONE, &
                                  WORK( J-1+N ), N, WR, WI, X, 2, SCALE, &
                                  XNORM, IERR )
!
!                    Scale X to avoid overflow when updating
!                    the right-hand side.
!
                     IF( XNORM.GT.ONE ) THEN
                        BETA = MAX( WORK( J-1 ), WORK( J ) )
                        IF( BETA.GT.BIGNUM / XNORM ) THEN
                           REC = ONE / XNORM
                           X( 1, 1 ) = X( 1, 1 )*REC
                           X( 1, 2 ) = X( 1, 2 )*REC
                           X( 2, 1 ) = X( 2, 1 )*REC
                           X( 2, 2 ) = X( 2, 2 )*REC
                           SCALE = SCALE*REC
                        END IF
                     END IF
!
!                    Scale if necessary
!
                     IF( SCALE.NE.ONE ) THEN
                        CALL DSCAL( KI, SCALE, WORK( 1+N ), 1 )
                        CALL DSCAL( KI, SCALE, WORK( 1+N2 ), 1 )
                     END IF
                     WORK( J-1+N ) = X( 1, 1 )
                     WORK( J+N ) = X( 2, 1 )
                     WORK( J-1+N2 ) = X( 1, 2 )
                     WORK( J+N2 ) = X( 2, 2 )
!
!                    Update the right-hand side
!
                     CALL DAXPY( J-2, -X( 1, 1 ), T( 1, J-1 ), 1, &
                                 WORK( 1+N ), 1 )
                     CALL DAXPY( J-2, -X( 2, 1 ), T( 1, J ), 1, &
                                 WORK( 1+N ), 1 )
                     CALL DAXPY( J-2, -X( 1, 2 ), T( 1, J-1 ), 1, &
                                 WORK( 1+N2 ), 1 )
                     CALL DAXPY( J-2, -X( 2, 2 ), T( 1, J ), 1, &
                                 WORK( 1+N2 ), 1 )
                  END IF
   90          CONTINUE
!
!              Copy the vector x or Q*x to VR and normalize.
!
               IF( .NOT.OVER ) THEN
                  CALL DCOPY( KI, WORK( 1+N ), 1, VR( 1, IS-1 ), 1 )
                  CALL DCOPY( KI, WORK( 1+N2 ), 1, VR( 1, IS ), 1 )
!
                  EMAX = ZERO
                  DO 100 K = 1, KI
                     EMAX = MAX( EMAX, ABS( VR( K, IS-1 ) )+ &
                            ABS( VR( K, IS ) ) )
  100             CONTINUE
!
                  REMAX = ONE / EMAX
                  CALL DSCAL( KI, REMAX, VR( 1, IS-1 ), 1 )
                  CALL DSCAL( KI, REMAX, VR( 1, IS ), 1 )
!
                  DO 110 K = KI + 1, N
                     VR( K, IS-1 ) = ZERO
                     VR( K, IS ) = ZERO
  110             CONTINUE
!
               ELSE
!
                  IF( KI.GT.2 ) THEN
                     CALL DGEMV( 'N', N, KI-2, ONE, VR, LDVR, &
                                 WORK( 1+N ), 1, WORK( KI-1+N ), &
                                 VR( 1, KI-1 ), 1 )
                     CALL DGEMV( 'N', N, KI-2, ONE, VR, LDVR, &
                                 WORK( 1+N2 ), 1, WORK( KI+N2 ), &
                                 VR( 1, KI ), 1 )
                  ELSE
                     CALL DSCAL( N, WORK( KI-1+N ), VR( 1, KI-1 ), 1 )
                     CALL DSCAL( N, WORK( KI+N2 ), VR( 1, KI ), 1 )
                  END IF
!
                  EMAX = ZERO
                  DO 120 K = 1, N
                     EMAX = MAX( EMAX, ABS( VR( K, KI-1 ) )+ &
                            ABS( VR( K, KI ) ) )
  120             CONTINUE
                  REMAX = ONE / EMAX
                  CALL DSCAL( N, REMAX, VR( 1, KI-1 ), 1 )
                  CALL DSCAL( N, REMAX, VR( 1, KI ), 1 )
               END IF
            END IF
!
            IS = IS - 1
            IF( IP.NE.0 ) &
               IS = IS - 1
  130       CONTINUE
            IF( IP.EQ.1 ) &
               IP = 0
            IF( IP.EQ.-1 ) &
               IP = 1
  140    CONTINUE
      END IF
!
      IF( LEFTV ) THEN
!
!        Compute left eigenvectors.
!
         IP = 0
         IS = 1
         DO 260 KI = 1, N
!
            IF( IP.EQ.-1 ) &
               GO TO 250
            IF( KI.EQ.N ) &
               GO TO 150
            IF( T( KI+1, KI ).EQ.ZERO ) &
               GO TO 150
            IP = 1
!
  150       CONTINUE
            IF( SOMEV ) THEN
               IF( .NOT.SELECT( KI ) ) &
                  GO TO 250
            END IF
!
!           Compute the KI-th eigenvalue (WR,WI).
!
            WR = T( KI, KI )
            WI = ZERO
            IF( IP.NE.0 ) &
               WI = SQRT( ABS( T( KI, KI+1 ) ) )* &
                    SQRT( ABS( T( KI+1, KI ) ) )
            SMIN = MAX( ULP*( ABS( WR )+ABS( WI ) ), SMLNUM )
!
            IF( IP.EQ.0 ) THEN
!
!              Real left eigenvector.
!
               WORK( KI+N ) = ONE
!
!              Form right-hand side
!
               DO 160 K = KI + 1, N
                  WORK( K+N ) = -T( KI, K )
  160          CONTINUE
!
!              Solve the quasi-triangular system:
!                 (T(KI+1:N,KI+1:N) - WR)'*X = SCALE*WORK
!
               VMAX = ONE
               VCRIT = BIGNUM
!
               JNXT = KI + 1
               DO 170 J = KI + 1, N
                  IF( J.LT.JNXT ) &
                     GO TO 170
                  J1 = J
                  J2 = J
                  JNXT = J + 1
                  IF( J.LT.N ) THEN
                     IF( T( J+1, J ).NE.ZERO ) THEN
                        J2 = J + 1
                        JNXT = J + 2
                     END IF
                  END IF
!
                  IF( J1.EQ.J2 ) THEN
!
!                    1-by-1 diagonal block
!
!                    Scale if necessary to avoid overflow when forming
!                    the right-hand side.
!
                     IF( WORK( J ).GT.VCRIT ) THEN
                        REC = ONE / VMAX
                        CALL DSCAL( N-KI+1, REC, WORK( KI+N ), 1 )
                        VMAX = ONE
                        VCRIT = BIGNUM
                     END IF
!
                     WORK( J+N ) = WORK( J+N ) - &
                                   DDOT( J-KI-1, T( KI+1, J ), 1, &
                                   WORK( KI+1+N ), 1 )
!
!                    Solve (T(J,J)-WR)'*X = WORK
!
                     CALL DLALN2( .FALSE., 1, 1, SMIN, ONE, T( J, J ), &
                                  LDT, ONE, ONE, WORK( J+N ), N, WR, &
                                  ZERO, X, 2, SCALE, XNORM, IERR )
!
!                    Scale if necessary
!
                     IF( SCALE.NE.ONE ) &
                        CALL DSCAL( N-KI+1, SCALE, WORK( KI+N ), 1 )
                     WORK( J+N ) = X( 1, 1 )
                     VMAX = MAX( ABS( WORK( J+N ) ), VMAX )
                     VCRIT = BIGNUM / VMAX
!
                  ELSE
!
!                    2-by-2 diagonal block
!
!                    Scale if necessary to avoid overflow when forming
!                    the right-hand side.
!
                     BETA = MAX( WORK( J ), WORK( J+1 ) )
                     IF( BETA.GT.VCRIT ) THEN
                        REC = ONE / VMAX
                        CALL DSCAL( N-KI+1, REC, WORK( KI+N ), 1 )
                        VMAX = ONE
                        VCRIT = BIGNUM
                     END IF
!
                     WORK( J+N ) = WORK( J+N ) - &
                                   DDOT( J-KI-1, T( KI+1, J ), 1, &
                                   WORK( KI+1+N ), 1 )
!
                     WORK( J+1+N ) = WORK( J+1+N ) - &
                                     DDOT( J-KI-1, T( KI+1, J+1 ), 1, &
                                     WORK( KI+1+N ), 1 )
!
!                    Solve
!                      [T(J,J)-WR   T(J,J+1)     ]'* X = SCALE*( WORK1 )
!                      [T(J+1,J)    T(J+1,J+1)-WR]             ( WORK2 )
!
                     CALL DLALN2( .TRUE., 2, 1, SMIN, ONE, T( J, J ), &
                                  LDT, ONE, ONE, WORK( J+N ), N, WR, &
                                  ZERO, X, 2, SCALE, XNORM, IERR )
!
!                    Scale if necessary
!
                     IF( SCALE.NE.ONE ) &
                        CALL DSCAL( N-KI+1, SCALE, WORK( KI+N ), 1 )
                     WORK( J+N ) = X( 1, 1 )
                     WORK( J+1+N ) = X( 2, 1 )
!
                     VMAX = MAX( ABS( WORK( J+N ) ), &
                            ABS( WORK( J+1+N ) ), VMAX )
                     VCRIT = BIGNUM / VMAX
!
                  END IF
  170          CONTINUE
!
!              Copy the vector x or Q*x to VL and normalize.
!
               IF( .NOT.OVER ) THEN
                  CALL DCOPY( N-KI+1, WORK( KI+N ), 1, VL( KI, IS ), 1 )
!
                  II = IDAMAX( N-KI+1, VL( KI, IS ), 1 ) + KI - 1
                  REMAX = ONE / ABS( VL( II, IS ) )
                  CALL DSCAL( N-KI+1, REMAX, VL( KI, IS ), 1 )
!
                  DO 180 K = 1, KI - 1
                     VL( K, IS ) = ZERO
  180             CONTINUE
!
               ELSE
!
                  IF( KI.LT.N ) &
                     CALL DGEMV( 'N', N, N-KI, ONE, VL( 1, KI+1 ), LDVL, &
                                 WORK( KI+1+N ), 1, WORK( KI+N ), &
                                 VL( 1, KI ), 1 )
!
                  II = IDAMAX( N, VL( 1, KI ), 1 )
                  REMAX = ONE / ABS( VL( II, KI ) )
                  CALL DSCAL( N, REMAX, VL( 1, KI ), 1 )
!
               END IF
!
            ELSE
!
!              Complex left eigenvector.
!
!               Initial solve:
!                 ((T(KI,KI)    T(KI,KI+1) )' - (WR - I* WI))*X = 0.
!                 ((T(KI+1,KI) T(KI+1,KI+1))                )
!
               IF( ABS( T( KI, KI+1 ) ).GE.ABS( T( KI+1, KI ) ) ) THEN
                  WORK( KI+N ) = WI / T( KI, KI+1 )
                  WORK( KI+1+N2 ) = ONE
               ELSE
                  WORK( KI+N ) = ONE
                  WORK( KI+1+N2 ) = -WI / T( KI+1, KI )
               END IF
               WORK( KI+1+N ) = ZERO
               WORK( KI+N2 ) = ZERO
!
!              Form right-hand side
!
               DO 190 K = KI + 2, N
                  WORK( K+N ) = -WORK( KI+N )*T( KI, K )
                  WORK( K+N2 ) = -WORK( KI+1+N2 )*T( KI+1, K )
  190          CONTINUE
!
!              Solve complex quasi-triangular system:
!              ( T(KI+2,N:KI+2,N) - (WR-i*WI) )*X = WORK1+i*WORK2
!
               VMAX = ONE
               VCRIT = BIGNUM
!
               JNXT = KI + 2
               DO 200 J = KI + 2, N
                  IF( J.LT.JNXT ) &
                     GO TO 200
                  J1 = J
                  J2 = J
                  JNXT = J + 1
                  IF( J.LT.N ) THEN
                     IF( T( J+1, J ).NE.ZERO ) THEN
                        J2 = J + 1
                        JNXT = J + 2
                     END IF
                  END IF
!
                  IF( J1.EQ.J2 ) THEN
!
!                    1-by-1 diagonal block
!
!                    Scale if necessary to avoid overflow when
!                    forming the right-hand side elements.
!
                     IF( WORK( J ).GT.VCRIT ) THEN
                        REC = ONE / VMAX
                        CALL DSCAL( N-KI+1, REC, WORK( KI+N ), 1 )
                        CALL DSCAL( N-KI+1, REC, WORK( KI+N2 ), 1 )
                        VMAX = ONE
                        VCRIT = BIGNUM
                     END IF
!
                     WORK( J+N ) = WORK( J+N ) - &
                                   DDOT( J-KI-2, T( KI+2, J ), 1, &
                                   WORK( KI+2+N ), 1 )
                     WORK( J+N2 ) = WORK( J+N2 ) - &
                                    DDOT( J-KI-2, T( KI+2, J ), 1, &
                                    WORK( KI+2+N2 ), 1 )
!
!                    Solve (T(J,J)-(WR-i*WI))*(X11+i*X12)= WK+I*WK2
!
                     CALL DLALN2( .FALSE., 1, 2, SMIN, ONE, T( J, J ), &
                                  LDT, ONE, ONE, WORK( J+N ), N, WR, &
                                  -WI, X, 2, SCALE, XNORM, IERR )
!
!                    Scale if necessary
!
                     IF( SCALE.NE.ONE ) THEN
                        CALL DSCAL( N-KI+1, SCALE, WORK( KI+N ), 1 )
                        CALL DSCAL( N-KI+1, SCALE, WORK( KI+N2 ), 1 )
                     END IF
                     WORK( J+N ) = X( 1, 1 )
                     WORK( J+N2 ) = X( 1, 2 )
                     VMAX = MAX( ABS( WORK( J+N ) ), &
                            ABS( WORK( J+N2 ) ), VMAX )
                     VCRIT = BIGNUM / VMAX
!
                  ELSE
!
!                    2-by-2 diagonal block
!
!                    Scale if necessary to avoid overflow when forming
!                    the right-hand side elements.
!
                     BETA = MAX( WORK( J ), WORK( J+1 ) )
                     IF( BETA.GT.VCRIT ) THEN
                        REC = ONE / VMAX
                        CALL DSCAL( N-KI+1, REC, WORK( KI+N ), 1 )
                        CALL DSCAL( N-KI+1, REC, WORK( KI+N2 ), 1 )
                        VMAX = ONE
                        VCRIT = BIGNUM
                     END IF
!
                     WORK( J+N ) = WORK( J+N ) - &
                                   DDOT( J-KI-2, T( KI+2, J ), 1, &
                                   WORK( KI+2+N ), 1 )
!
                     WORK( J+N2 ) = WORK( J+N2 ) - &
                                    DDOT( J-KI-2, T( KI+2, J ), 1, &
                                    WORK( KI+2+N2 ), 1 )
!
                     WORK( J+1+N ) = WORK( J+1+N ) - &
                                     DDOT( J-KI-2, T( KI+2, J+1 ), 1, &
                                     WORK( KI+2+N ), 1 )
!
                     WORK( J+1+N2 ) = WORK( J+1+N2 ) - &
                                      DDOT( J-KI-2, T( KI+2, J+1 ), 1, &
                                      WORK( KI+2+N2 ), 1 )
!
!                    Solve 2-by-2 complex linear equation
!                      ([T(j,j)   T(j,j+1)  ]'-(wr-i*wi)*I)*X = SCALE*B
!                      ([T(j+1,j) T(j+1,j+1)]             )
!
                     CALL DLALN2( .TRUE., 2, 2, SMIN, ONE, T( J, J ), &
                                  LDT, ONE, ONE, WORK( J+N ), N, WR, &
                                  -WI, X, 2, SCALE, XNORM, IERR )
!
!                    Scale if necessary
!
                     IF( SCALE.NE.ONE ) THEN
                        CALL DSCAL( N-KI+1, SCALE, WORK( KI+N ), 1 )
                        CALL DSCAL( N-KI+1, SCALE, WORK( KI+N2 ), 1 )
                     END IF
                     WORK( J+N ) = X( 1, 1 )
                     WORK( J+N2 ) = X( 1, 2 )
                     WORK( J+1+N ) = X( 2, 1 )
                     WORK( J+1+N2 ) = X( 2, 2 )
                     VMAX = MAX( ABS( X( 1, 1 ) ), ABS( X( 1, 2 ) ), &
                            ABS( X( 2, 1 ) ), ABS( X( 2, 2 ) ), VMAX )
                     VCRIT = BIGNUM / VMAX
!
                  END IF
  200          CONTINUE
!
!              Copy the vector x or Q*x to VL and normalize.
!
               IF( .NOT.OVER ) THEN
                  CALL DCOPY( N-KI+1, WORK( KI+N ), 1, VL( KI, IS ), 1 )
                  CALL DCOPY( N-KI+1, WORK( KI+N2 ), 1, VL( KI, IS+1 ), &
                              1 )
!
                  EMAX = ZERO
                  DO 220 K = KI, N
                     EMAX = MAX( EMAX, ABS( VL( K, IS ) )+ &
                            ABS( VL( K, IS+1 ) ) )
  220             CONTINUE
                  REMAX = ONE / EMAX
                  CALL DSCAL( N-KI+1, REMAX, VL( KI, IS ), 1 )
                  CALL DSCAL( N-KI+1, REMAX, VL( KI, IS+1 ), 1 )
!
                  DO 230 K = 1, KI - 1
                     VL( K, IS ) = ZERO
                     VL( K, IS+1 ) = ZERO
  230             CONTINUE
               ELSE
                  IF( KI.LT.N-1 ) THEN
                     CALL DGEMV( 'N', N, N-KI-1, ONE, VL( 1, KI+2 ), &
                                 LDVL, WORK( KI+2+N ), 1, WORK( KI+N ), &
                                 VL( 1, KI ), 1 )
                     CALL DGEMV( 'N', N, N-KI-1, ONE, VL( 1, KI+2 ), &
                                 LDVL, WORK( KI+2+N2 ), 1, &
                                 WORK( KI+1+N2 ), VL( 1, KI+1 ), 1 )
                  ELSE
                     CALL DSCAL( N, WORK( KI+N ), VL( 1, KI ), 1 )
                     CALL DSCAL( N, WORK( KI+1+N2 ), VL( 1, KI+1 ), 1 )
                  END IF
!
                  EMAX = ZERO
                  DO 240 K = 1, N
                     EMAX = MAX( EMAX, ABS( VL( K, KI ) )+ &
                            ABS( VL( K, KI+1 ) ) )
  240             CONTINUE
                  REMAX = ONE / EMAX
                  CALL DSCAL( N, REMAX, VL( 1, KI ), 1 )
                  CALL DSCAL( N, REMAX, VL( 1, KI+1 ), 1 )
!
               END IF
!
            END IF
!
            IS = IS + 1
            IF( IP.NE.0 ) &
               IS = IS + 1
  250       CONTINUE
            IF( IP.EQ.-1 ) &
               IP = 0
            IF( IP.EQ.1 ) &
               IP = -1
!
  260    CONTINUE
!
      END IF
!
      RETURN
!
!     End of DTREVC
!
      END SUBROUTINE DTREVC


      SUBROUTINE DHSEQR( JOB, COMPZ, N, ILO, IHI, H, LDH, WR, WI, Z, &
                         LDZ, WORK, LWORK, INFO )
!
!  -- LAPACK driver routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDH, LDZ, LWORK, N
      CHARACTER          COMPZ, JOB
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), WI( * ), WORK( * ), WR( * ), &
                         Z( LDZ, * )
!     ..
!     Purpose
!     =======
!
!     DHSEQR computes the eigenvalues of a Hessenberg matrix H
!     and, optionally, the matrices T and Z from the Schur decomposition
!     H = Z T Z**T, where T is an upper quasi-triangular matrix (the
!     Schur form), and Z is the orthogonal matrix of Schur vectors.
!
!     Optionally Z may be postmultiplied into an input orthogonal
!     matrix Q so that this routine can give the Schur factorization
!     of a matrix A which has been reduced to the Hessenberg form H
!     by the orthogonal matrix Q:  A = Q*H*Q**T = (QZ)*T*(QZ)**T.
!
!     Arguments
!     =========
!
!     JOB   (input) CHARACTER*1
!           = 'E':  compute eigenvalues only;
!           = 'S':  compute eigenvalues and the Schur form T.
!
!     COMPZ (input) CHARACTER*1
!           = 'N':  no Schur vectors are computed;
!           = 'I':  Z is initialized to the unit matrix and the matrix Z
!                   of Schur vectors of H is returned;
!           = 'V':  Z must contain an orthogonal matrix Q on entry, and
!                   the product Q*Z is returned.
!
!     N     (input) INTEGER
!           The order of the matrix H.  N .GE. 0.
!
!     ILO   (input) INTEGER
!     IHI   (input) INTEGER
!           It is assumed that H is already upper triangular in rows
!           and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
!           set by a previous call to DGEBAL, and then passed to DGEHRD
!           when the matrix output by DGEBAL is reduced to Hessenberg
!           form. Otherwise ILO and IHI should be set to 1 and N
!           respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N.
!           If N = 0, then ILO = 1 and IHI = 0.
!
!     H     (input/output) DOUBLE PRECISION array, dimension (LDH,N)
!           On entry, the upper Hessenberg matrix H.
!           On exit, if INFO = 0 and JOB = 'S', then H contains the
!           upper quasi-triangular matrix T from the Schur decomposition
!           (the Schur form); 2-by-2 diagonal blocks (corresponding to
!           complex conjugate pairs of eigenvalues) are returned in
!           standard form, with H(i,i) = H(i+1,i+1) and
!           H(i+1,i)*H(i,i+1).LT.0. If INFO = 0 and JOB = 'E', the
!           contents of H are unspecified on exit.  (The output value of
!           H when INFO.GT.0 is given under the description of INFO
!           below.)
!
!           Unlike earlier versions of DHSEQR, this subroutine may
!           explicitly H(i,j) = 0 for i.GT.j and j = 1, 2, ... ILO-1
!           or j = IHI+1, IHI+2, ... N.
!
!     LDH   (input) INTEGER
!           The leading dimension of the array H. LDH .GE. max(1,N).
!
!     WR    (output) DOUBLE PRECISION array, dimension (N)
!     WI    (output) DOUBLE PRECISION array, dimension (N)
!           The real and imaginary parts, respectively, of the computed
!           eigenvalues. If two eigenvalues are computed as a complex
!           conjugate pair, they are stored in consecutive elements of
!           WR and WI, say the i-th and (i+1)th, with WI(i) .GT. 0 and
!           WI(i+1) .LT. 0. If JOB = 'S', the eigenvalues are stored in
!           the same order as on the diagonal of the Schur form returned
!           in H, with WR(i) = H(i,i) and, if H(i:i+1,i:i+1) is a 2-by-2
!           diagonal block, WI(i) = sqrt(-H(i+1,i)*H(i,i+1)) and
!           WI(i+1) = -WI(i).
!
!     Z     (input/output) DOUBLE PRECISION array, dimension (LDZ,N)
!           If COMPZ = 'N', Z is not referenced.
!           If COMPZ = 'I', on entry Z need not be set and on exit,
!           if INFO = 0, Z contains the orthogonal matrix Z of the Schur
!           vectors of H.  If COMPZ = 'V', on entry Z must contain an
!           N-by-N matrix Q, which is assumed to be equal to the unit
!           matrix except for the submatrix Z(ILO:IHI,ILO:IHI). On exit,
!           if INFO = 0, Z contains Q*Z.
!           Normally Q is the orthogonal matrix generated by DORGHR
!           after the call to DGEHRD which formed the Hessenberg matrix
!           H. (The output value of Z when INFO.GT.0 is given under
!           the description of INFO below.)
!
!     LDZ   (input) INTEGER
!           The leading dimension of the array Z.  if COMPZ = 'I' or
!           COMPZ = 'V', then LDZ.GE.MAX(1,N).  Otherwize, LDZ.GE.1.
!
!     WORK  (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!           On exit, if INFO = 0, WORK(1) returns an estimate of
!           the optimal value for LWORK.
!
!     LWORK (input) INTEGER
!           The dimension of the array WORK.  LWORK .GE. max(1,N)
!           is sufficient, but LWORK typically as large as 6*N may
!           be required for optimal performance.  A workspace query
!           to determine the optimal workspace size is recommended.
!
!           If LWORK = -1, then DHSEQR does a workspace query.
!           In this case, DHSEQR checks the input parameters and
!           estimates the optimal workspace size for the given
!           values of N, ILO and IHI.  The estimate is returned
!           in WORK(1).  No error message related to LWORK is
!           issued by XERBLA.  Neither H nor Z are accessed.
!
!
!     INFO  (output) INTEGER
!             =  0:  successful exit
!           .LT. 0:  if INFO = -i, the i-th argument had an illegal
!                    value
!           .GT. 0:  if INFO = i, DHSEQR failed to compute all of
!                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR
!                and WI contain those eigenvalues which have been
!                successfully computed.  (Failures are rare.)
!
!                If INFO .GT. 0 and JOB = 'E', then on exit, the
!                remaining unconverged eigenvalues are the eigen-
!                values of the upper Hessenberg matrix rows and
!                columns ILO through INFO of the final, output
!                value of H.
!
!                If INFO .GT. 0 and JOB   = 'S', then on exit
!
!           (*)  (initial value of H)*U  = U*(final value of H)
!
!                where U is an orthogonal matrix.  The final
!                value of H is upper Hessenberg and quasi-triangular
!                in rows and columns INFO+1 through IHI.
!
!                If INFO .GT. 0 and COMPZ = 'V', then on exit
!
!                  (final value of Z)  =  (initial value of Z)*U
!
!                where U is the orthogonal matrix in (*) (regard-
!                less of the value of JOB.)
!
!                If INFO .GT. 0 and COMPZ = 'I', then on exit
!                      (final value of Z)  = U
!                where U is the orthogonal matrix in (*) (regard-
!                less of the value of JOB.)
!
!                If INFO .GT. 0 and COMPZ = 'N', then Z is not
!                accessed.
!
!     ================================================================
!             Default values supplied by
!             ILAENV(ISPEC,'DHSEQR',JOB(:1)//COMPZ(:1),N,ILO,IHI,LWORK).
!             It is suggested that these defaults be adjusted in order
!             to attain best performance in each particular
!             computational environment.
!
!            ISPEC=1:  The DLAHQR vs DLAQR0 crossover point.
!                      Default: 75. (Must be at least 11.)
!
!            ISPEC=2:  Recommended deflation window size.
!                      This depends on ILO, IHI and NS.  NS is the
!                      number of simultaneous shifts returned
!                      by ILAENV(ISPEC=4).  (See ISPEC=4 below.)
!                      The default for (IHI-ILO+1).LE.500 is NS.
!                      The default for (IHI-ILO+1).GT.500 is 3*NS/2.
!
!            ISPEC=3:  Nibble crossover point. (See ILAENV for
!                      details.)  Default: 14% of deflation window
!                      size.
!
!            ISPEC=4:  Number of simultaneous shifts, NS, in
!                      a multi-shift QR iteration.
!
!                      If IHI-ILO+1 is ...
!
!                      greater than      ...but less    ... the
!                      or equal to ...      than        default is
!
!                           1               30          NS -   2(+)
!                          30               60          NS -   4(+)
!                          60              150          NS =  10(+)
!                         150              590          NS =  **
!                         590             3000          NS =  64
!                        3000             6000          NS = 128
!                        6000             infinity      NS = 256
!
!                  (+)  By default some or all matrices of this order
!                       are passed to the implicit double shift routine
!                       DLAHQR and NS is ignored.  See ISPEC=1 above
!                       and comments in IPARM for details.
!
!                       The asterisks (**) indicate an ad-hoc
!                       function of N increasing from 10 to 64.
!
!            ISPEC=5:  Select structured matrix multiply.
!                      (See ILAENV for details.) Default: 3.
!
!     ================================================================
!     Based on contributions by
!        Karen Braman and Ralph Byers, Department of Mathematics,
!        University of Kansas, USA
!
!     ================================================================
!     References:
!       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
!       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3
!       Performance, SIAM Journal of Matrix Analysis, volume 23, pages
!       929--947, 2002.
!
!       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
!       Algorithm Part II: Aggressive Early Deflation, SIAM Journal
!       of Matrix Analysis, volume 23, pages 948--973, 2002.
!
!     ================================================================
!     .. Parameters ..
!
!     ==== Matrices of order NTINY or smaller must be processed by
!     .    DLAHQR because of insufficient subdiagonal scratch space.
!     .    (This is a hard limit.) ====
!
!     ==== NL allocates some local workspace to help small matrices
!     .    through a rare DLAHQR failure.  NL .GT. NTINY = 11 is
!     .    required and NL .LE. NMIN = ILAENV(ISPEC=1,...) is recom-
!     .    mended.  (The default value of NMIN is 75.)  Using NL = 49
!     .    allows up to six simultaneous shifts and a 16-by-16
!     .    deflation window.  ====
!
      INTEGER            NTINY
      PARAMETER          ( NTINY = 11 )
      INTEGER            NL
      PARAMETER          ( NL = 49 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0d0, ONE = 1.0d0 )
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   HL( NL, NL ), WORKL( NL )
!     ..
!     .. Local Scalars ..
      INTEGER            I, KBOT, NMIN
      LOGICAL            INITZ, LQUERY, WANTT, WANTZ
!     ..
!     .. External Functions ..
!      INTEGER            ILAENV
!      LOGICAL            LSAME
!      EXTERNAL           ILAENV, LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DLACPY, DLAHQR, DLAQR0, DLASET, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     ==== Decode and check the input parameters. ====
!
      WANTT = LSAME( JOB, 'S' )
      INITZ = LSAME( COMPZ, 'I' )
      WANTZ = INITZ .OR. LSAME( COMPZ, 'V' )
      WORK( 1 ) = DBLE( MAX( 1, N ) )
      LQUERY = LWORK.EQ.-1
!
      INFO = 0
      IF( .NOT.LSAME( JOB, 'E' ) .AND. .NOT.WANTT ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( COMPZ, 'N' ) .AND. .NOT.WANTZ ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -5
      ELSE IF( LDH.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.MAX( 1, N ) ) ) THEN
         INFO = -11
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -13
      END IF
!
      IF( INFO.NE.0 ) THEN
!
!        ==== Quick return in case of invalid argument. ====
!
         CALL XERBLA( 'DHSEQR', -INFO )
         RETURN
!
      ELSE IF( N.EQ.0 ) THEN
!
!        ==== Quick return in case N = 0; nothing to do. ====
!
         RETURN
!
      ELSE IF( LQUERY ) THEN
!
!        ==== Quick return in case of a workspace query ====
!
         CALL DLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILO, &
                      IHI, Z, LDZ, WORK, LWORK, INFO )
!        ==== Ensure reported workspace size is backward-compatible with
!        .    previous LAPACK versions. ====
         WORK( 1 ) = MAX( DBLE( MAX( 1, N ) ), WORK( 1 ) )
         RETURN
!
      ELSE
!
!        ==== copy eigenvalues isolated by DGEBAL ====
!
         DO 10 I = 1, ILO - 1
            WR( I ) = H( I, I )
            WI( I ) = ZERO
   10    CONTINUE
         DO 20 I = IHI + 1, N
            WR( I ) = H( I, I )
            WI( I ) = ZERO
   20    CONTINUE
!
!        ==== Initialize Z, if requested ====
!
         IF( INITZ ) &
            CALL DLASET( 'A', N, N, ZERO, ONE, Z, LDZ )
!
!        ==== Quick return if possible ====
!
         IF( ILO.EQ.IHI ) THEN
            WR( ILO ) = H( ILO, ILO )
            WI( ILO ) = ZERO
            RETURN
         END IF
!
!        ==== DLAHQR/DLAQR0 crossover point ====
!
         NMIN = ILAENV( 12, 'DHSEQR', JOB( : 1 ) // COMPZ( : 1 ), N, &
                ILO, IHI, LWORK )
         NMIN = MAX( NTINY, NMIN )
!
!        ==== DLAQR0 for big matrices; DLAHQR for small ones ====
!
         IF( N.GT.NMIN ) THEN
            CALL DLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILO, &
                         IHI, Z, LDZ, WORK, LWORK, INFO )
         ELSE
!
!           ==== Small matrix ====
!
            CALL DLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILO, &
                         IHI, Z, LDZ, INFO )
!
            IF( INFO.GT.0 ) THEN
!
!              ==== A rare DLAHQR failure!  DLAQR0 sometimes succeeds
!              .    when DLAHQR fails. ====
!
               KBOT = INFO
!
               IF( N.GE.NL ) THEN
!
!                 ==== Larger matrices have enough subdiagonal scratch
!                 .    space to call DLAQR0 directly. ====
!
                  CALL DLAQR0( WANTT, WANTZ, N, ILO, KBOT, H, LDH, WR, &
                               WI, ILO, IHI, Z, LDZ, WORK, LWORK, INFO )
!
               ELSE
!
!                 ==== Tiny matrices don't have enough subdiagonal
!                 .    scratch space to benefit from DLAQR0.  Hence,
!                 .    tiny matrices must be copied into a larger
!                 .    array before calling DLAQR0. ====
!
                  CALL DLACPY( 'A', N, N, H, LDH, HL, NL )
                  HL( N+1, N ) = ZERO
                  CALL DLASET( 'A', NL, NL-N, ZERO, ZERO, HL( 1, N+1 ), &
                               NL )
                  CALL DLAQR0( WANTT, WANTZ, NL, ILO, KBOT, HL, NL, WR, &
                               WI, ILO, IHI, Z, LDZ, WORKL, NL, INFO )
                  IF( WANTT .OR. INFO.NE.0 ) &
                     CALL DLACPY( 'A', N, N, HL, NL, H, LDH )
               END IF
            END IF
         END IF
!
!        ==== Clear out the trash, if necessary. ====
!
         IF( ( WANTT .OR. INFO.NE.0 ) .AND. N.GT.2 ) &
            CALL DLASET( 'L', N-2, N-2, ZERO, ZERO, H( 3, 1 ), LDH )
!
!        ==== Ensure reported workspace size is backward-compatible with
!        .    previous LAPACK versions. ====
!
         WORK( 1 ) = MAX( DBLE( MAX( 1, N ) ), WORK( 1 ) )
      END IF
!
!     ==== End of DHSEQR ====
!
      END SUBROUTINE DHSEQR


      SUBROUTINE DLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          TYPE
      INTEGER            INFO, KL, KU, LDA, M, N
      DOUBLE PRECISION   CFROM, CTO
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DLASCL multiplies the M by N real matrix A by the real scalar
!  CTO/CFROM.  This is done without over/underflow as long as the final
!  result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
!  A may be full, upper triangular, lower triangular, upper Hessenberg,
!  or banded.
!
!  Arguments
!  =========
!
!  TYPE    (input) CHARACTER*1
!          TYPE indices the storage type of the input matrix.
!          = 'G':  A is a full matrix.
!          = 'L':  A is a lower triangular matrix.
!          = 'U':  A is an upper triangular matrix.
!          = 'H':  A is an upper Hessenberg matrix.
!          = 'B':  A is a symmetric band matrix with lower bandwidth KL
!                  and upper bandwidth KU and with the only the lower
!                  half stored.
!          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
!                  and upper bandwidth KU and with the only the upper
!                  half stored.
!          = 'Z':  A is a band matrix with lower bandwidth KL and upper
!                  bandwidth KU.
!
!  KL      (input) INTEGER
!          The lower bandwidth of A.  Referenced only if TYPE = 'B',
!          'Q' or 'Z'.
!
!  KU      (input) INTEGER
!          The upper bandwidth of A.  Referenced only if TYPE = 'B',
!          'Q' or 'Z'.
!
!  CFROM   (input) DOUBLE PRECISION
!  CTO     (input) DOUBLE PRECISION
!          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
!          without over/underflow if the final result CTO*A(I,J)/CFROM
!          can be represented without over/underflow.  CFROM must be
!          nonzero.
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
!          storage type.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  INFO    (output) INTEGER
!          0  - successful exit
!          <0 - if INFO = -i, the i-th argument had an illegal value.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            DONE
      INTEGER            I, ITYPE, J, K1, K2, K3, K4
      DOUBLE PRECISION   BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      DOUBLE PRECISION   DLAMCH
!      EXTERNAL           LSAME, DLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
!     ..
!     .. External Subroutines ..
!      EXTERNAL           XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
!
      IF( LSAME( TYPE, 'G' ) ) THEN
         ITYPE = 0
      ELSE IF( LSAME( TYPE, 'L' ) ) THEN
         ITYPE = 1
      ELSE IF( LSAME( TYPE, 'U' ) ) THEN
         ITYPE = 2
      ELSE IF( LSAME( TYPE, 'H' ) ) THEN
         ITYPE = 3
      ELSE IF( LSAME( TYPE, 'B' ) ) THEN
         ITYPE = 4
      ELSE IF( LSAME( TYPE, 'Q' ) ) THEN
         ITYPE = 5
      ELSE IF( LSAME( TYPE, 'Z' ) ) THEN
         ITYPE = 6
      ELSE
         ITYPE = -1
      END IF
!
      IF( ITYPE.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( CFROM.EQ.ZERO ) THEN
         INFO = -4
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( N.LT.0 .OR. ( ITYPE.EQ.4 .AND. N.NE.M ) .OR. &
               ( ITYPE.EQ.5 .AND. N.NE.M ) ) THEN
         INFO = -7
      ELSE IF( ITYPE.LE.3 .AND. LDA.LT.MAX( 1, M ) ) THEN
         INFO = -9
      ELSE IF( ITYPE.GE.4 ) THEN
         IF( KL.LT.0 .OR. KL.GT.MAX( M-1, 0 ) ) THEN
            INFO = -2
         ELSE IF( KU.LT.0 .OR. KU.GT.MAX( N-1, 0 ) .OR. &
                  ( ( ITYPE.EQ.4 .OR. ITYPE.EQ.5 ) .AND. KL.NE.KU ) ) &
                   THEN
            INFO = -3
         ELSE IF( ( ITYPE.EQ.4 .AND. LDA.LT.KL+1 ) .OR. &
                  ( ITYPE.EQ.5 .AND. LDA.LT.KU+1 ) .OR. &
                  ( ITYPE.EQ.6 .AND. LDA.LT.2*KL+KU+1 ) ) THEN
            INFO = -9
         END IF
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASCL', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. M.EQ.0 ) &
         RETURN
!
!     Get machine parameters
!
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
!
      CFROMC = CFROM
      CTOC = CTO
!
   10 CONTINUE
      CFROM1 = CFROMC*SMLNUM
      CTO1 = CTOC / BIGNUM
      IF( ABS( CFROM1 ).GT.ABS( CTOC ) .AND. CTOC.NE.ZERO ) THEN
         MUL = SMLNUM
         DONE = .FALSE.
         CFROMC = CFROM1
      ELSE IF( ABS( CTO1 ).GT.ABS( CFROMC ) ) THEN
         MUL = BIGNUM
         DONE = .FALSE.
         CTOC = CTO1
      ELSE
         MUL = CTOC / CFROMC
         DONE = .TRUE.
      END IF
!
      IF( ITYPE.EQ.0 ) THEN
!
!        Full matrix
!
         DO 30 J = 1, N
            DO 20 I = 1, M
               A( I, J ) = A( I, J )*MUL
   20       CONTINUE
   30    CONTINUE
!
      ELSE IF( ITYPE.EQ.1 ) THEN
!
!        Lower triangular matrix
!
         DO 50 J = 1, N
            DO 40 I = J, M
               A( I, J ) = A( I, J )*MUL
   40       CONTINUE
   50    CONTINUE
!
      ELSE IF( ITYPE.EQ.2 ) THEN
!
!        Upper triangular matrix
!
         DO 70 J = 1, N
            DO 60 I = 1, MIN( J, M )
               A( I, J ) = A( I, J )*MUL
   60       CONTINUE
   70    CONTINUE
!
      ELSE IF( ITYPE.EQ.3 ) THEN
!
!        Upper Hessenberg matrix
!
         DO 90 J = 1, N
            DO 80 I = 1, MIN( J+1, M )
               A( I, J ) = A( I, J )*MUL
   80       CONTINUE
   90    CONTINUE
!
      ELSE IF( ITYPE.EQ.4 ) THEN
!
!        Lower half of a symmetric band matrix
!
         K3 = KL + 1
         K4 = N + 1
         DO 110 J = 1, N
            DO 100 I = 1, MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  100       CONTINUE
  110    CONTINUE
!
      ELSE IF( ITYPE.EQ.5 ) THEN
!
!        Upper half of a symmetric band matrix
!
         K1 = KU + 2
         K3 = KU + 1
         DO 130 J = 1, N
            DO 120 I = MAX( K1-J, 1 ), K3
               A( I, J ) = A( I, J )*MUL
  120       CONTINUE
  130    CONTINUE
!
      ELSE IF( ITYPE.EQ.6 ) THEN
!
!        Band matrix
!
         K1 = KL + KU + 2
         K2 = KL + 1
         K3 = 2*KL + KU + 1
         K4 = KL + KU + 1 + M
         DO 150 J = 1, N
            DO 140 I = MAX( K1-J, K2 ), MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  140       CONTINUE
  150    CONTINUE
!
      END IF
!
      IF( .NOT.DONE ) &
         GO TO 10
!
      RETURN
!
!     End of DLASCL
!
      END SUBROUTINE DLASCL


      SUBROUTINE DLARTG( F, G, CS, SN, R )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   CS, F, G, R, SN
!     ..
!
!  Purpose
!  =======
!
!  DLARTG generate a plane rotation so that
!
!     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
!     [ -SN  CS  ]     [ G ]     [ 0 ]
!
!  This is a slower, more accurate version of the BLAS1 routine DROTG,
!  with the following other differences:
!     F and G are unchanged on return.
!     If G=0, then CS=1 and SN=0.
!     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
!        floating point operations (saves work in DBDSQR when
!        there are zeros on the diagonal).
!
!  If F exceeds G in magnitude, CS will be positive.
!
!  Arguments
!  =========
!
!  F       (input) DOUBLE PRECISION
!          The first component of vector to be rotated.
!
!  G       (input) DOUBLE PRECISION
!          The second component of vector to be rotated.
!
!  CS      (output) DOUBLE PRECISION
!          The cosine of the rotation.
!
!  SN      (output) DOUBLE PRECISION
!          The sine of the rotation.
!
!  R       (output) DOUBLE PRECISION
!          The nonzero component of the rotated vector.
!
!  This version has a few statements commented out for thread safety
!  (machine parameters are computed on each entry). 10 feb 03, SJH.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
!     ..
!     .. Local Scalars ..
!     LOGICAL            FIRST
      INTEGER            COUNT, I
      DOUBLE PRECISION   EPS, F1, G1, SAFMIN, SAFMN2, SAFMX2, SCALE
!     ..
!     .. External Functions ..
!      DOUBLE PRECISION   DLAMCH
!      EXTERNAL           DLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, SQRT
!     ..
!     .. Save statement ..
!     SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
!     ..
!     .. Data statements ..
!     DATA               FIRST / .TRUE. /
!     ..
!     .. Executable Statements ..
!
!     IF( FIRST ) THEN
         SAFMIN = DLAMCH( 'S' )
         EPS = DLAMCH( 'E' )
         SAFMN2 = DLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) / &
                  LOG( DLAMCH( 'B' ) ) / TWO )
         SAFMX2 = ONE / SAFMN2
!        FIRST = .FALSE.
!     END IF
      IF( G.EQ.ZERO ) THEN
         CS = ONE
         SN = ZERO
         R = F
      ELSE IF( F.EQ.ZERO ) THEN
         CS = ZERO
         SN = ONE
         R = G
      ELSE
         F1 = F
         G1 = G
         SCALE = MAX( ABS( F1 ), ABS( G1 ) )
         IF( SCALE.GE.SAFMX2 ) THEN
            COUNT = 0
   10       CONTINUE
            COUNT = COUNT + 1
            F1 = F1*SAFMN2
            G1 = G1*SAFMN2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE.GE.SAFMX2 ) &
               GO TO 10
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
            DO 20 I = 1, COUNT
               R = R*SAFMX2
   20       CONTINUE
         ELSE IF( SCALE.LE.SAFMN2 ) THEN
            COUNT = 0
   30       CONTINUE
            COUNT = COUNT + 1
            F1 = F1*SAFMX2
            G1 = G1*SAFMX2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE.LE.SAFMN2 ) &
               GO TO 30
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
            DO 40 I = 1, COUNT
               R = R*SAFMN2
   40       CONTINUE
         ELSE
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
         END IF
         IF( ABS( F ).GT.ABS( G ) .AND. CS.LT.ZERO ) THEN
            CS = -CS
            SN = -SN
            R = -R
         END IF
      END IF
      RETURN
!
!     End of DLARTG
!
      END SUBROUTINE DLARTG


      SUBROUTINE DGEBAK( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV, &
                         INFO )
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          JOB, SIDE
      INTEGER            IHI, ILO, INFO, LDV, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   SCALE( * ), V( LDV, * )
!     ..
!
!  Purpose
!  =======
!
!  DGEBAK forms the right or left eigenvectors of a real general matrix
!  by backward transformation on the computed eigenvectors of the
!  balanced matrix output by DGEBAL.
!
!  Arguments
!  =========
!
!  JOB     (input) CHARACTER*1
!          Specifies the type of backward transformation required:
!          = 'N', do nothing, return immediately;
!          = 'P', do backward transformation for permutation only;
!          = 'S', do backward transformation for scaling only;
!          = 'B', do backward transformations for both permutation and
!                 scaling.
!          JOB must be the same as the argument JOB supplied to DGEBAL.
!
!  SIDE    (input) CHARACTER*1
!          = 'R':  V contains right eigenvectors;
!          = 'L':  V contains left eigenvectors.
!
!  N       (input) INTEGER
!          The number of rows of the matrix V.  N >= 0.
!
!  ILO     (input) INTEGER
!  IHI     (input) INTEGER
!          The integers ILO and IHI determined by DGEBAL.
!          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
!
!  SCALE   (input) DOUBLE PRECISION array, dimension (N)
!          Details of the permutation and scaling factors, as returned
!          by DGEBAL.
!
!  M       (input) INTEGER
!          The number of columns of the matrix V.  M >= 0.
!
!  V       (input/output) DOUBLE PRECISION array, dimension (LDV,M)
!          On entry, the matrix of right or left eigenvectors to be
!          transformed, as returned by DHSEIN or DTREVC.
!          On exit, V is overwritten by the transformed eigenvectors.
!
!  LDV     (input) INTEGER
!          The leading dimension of the array V. LDV >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFTV, RIGHTV
      INTEGER            I, II, K
      DOUBLE PRECISION   S
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DSCAL, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Decode and Test the input parameters
!
      RIGHTV = LSAME( SIDE, 'R' )
      LEFTV = LSAME( SIDE, 'L' )
!
      INFO = 0
      IF( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.LSAME( JOB, 'P' ) .AND. &
          .NOT.LSAME( JOB, 'S' ) .AND. .NOT.LSAME( JOB, 'B' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.RIGHTV .AND. .NOT.LEFTV ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -5
      ELSE IF( M.LT.0 ) THEN
         INFO = -7
      ELSE IF( LDV.LT.MAX( 1, N ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEBAK', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
      IF( M.EQ.0 ) &
         RETURN
      IF( LSAME( JOB, 'N' ) ) &
         RETURN
!
      IF( ILO.EQ.IHI ) &
         GO TO 30
!
!     Backward balance
!
      IF( LSAME( JOB, 'S' ) .OR. LSAME( JOB, 'B' ) ) THEN
!
         IF( RIGHTV ) THEN
            DO 10 I = ILO, IHI
               S = SCALE( I )
               CALL DSCAL( M, S, V( I, 1 ), LDV )
   10       CONTINUE
         END IF
!
         IF( LEFTV ) THEN
            DO 20 I = ILO, IHI
               S = ONE / SCALE( I )
               CALL DSCAL( M, S, V( I, 1 ), LDV )
   20       CONTINUE
         END IF
!
      END IF
!
!     Backward permutation
!
!     For  I = ILO-1 step -1 until 1,
!              IHI+1 step 1 until N do --
!
   30 CONTINUE
      IF( LSAME( JOB, 'P' ) .OR. LSAME( JOB, 'B' ) ) THEN
         IF( RIGHTV ) THEN
            DO 40 II = 1, N
               I = II
               IF( I.GE.ILO .AND. I.LE.IHI ) &
                  GO TO 40
               IF( I.LT.ILO ) &
                  I = ILO - II
               K = SCALE( I )
               IF( K.EQ.I ) &
                  GO TO 40
               CALL DSWAP( M, V( I, 1 ), LDV, V( K, 1 ), LDV )
   40       CONTINUE
         END IF
!
         IF( LEFTV ) THEN
            DO 50 II = 1, N
               I = II
               IF( I.GE.ILO .AND. I.LE.IHI ) &
                  GO TO 50
               IF( I.LT.ILO ) &
                  I = ILO - II
               K = SCALE( I )
               IF( K.EQ.I ) &
                  GO TO 50
               CALL DSWAP( M, V( I, 1 ), LDV, V( K, 1 ), LDV )
   50       CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of DGEBAK
!
      END SUBROUTINE DGEBAK


      SUBROUTINE DGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO )
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          JOB
      INTEGER            IHI, ILO, INFO, LDA, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), SCALE( * )
!     ..
!
!  Purpose
!  =======
!
!  DGEBAL balances a general real matrix A.  This involves, first,
!  permuting A by a similarity transformation to isolate eigenvalues
!  in the first 1 to ILO-1 and last IHI+1 to N elements on the
!  diagonal; and second, applying a diagonal similarity transformation
!  to rows and columns ILO to IHI to make the rows and columns as
!  close in norm as possible.  Both steps are optional.
!
!  Balancing may reduce the 1-norm of the matrix, and improve the
!  accuracy of the computed eigenvalues and/or eigenvectors.
!
!  Arguments
!  =========
!
!  JOB     (input) CHARACTER*1
!          Specifies the operations to be performed on A:
!          = 'N':  none:  simply set ILO = 1, IHI = N, SCALE(I) = 1.0
!                  for i = 1,...,N;
!          = 'P':  permute only;
!          = 'S':  scale only;
!          = 'B':  both permute and scale.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the input matrix A.
!          On exit,  A is overwritten by the balanced matrix.
!          If JOB = 'N', A is not referenced.
!          See Further Details.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  ILO     (output) INTEGER
!  IHI     (output) INTEGER
!          ILO and IHI are set to integers such that on exit
!          A(i,j) = 0 if i > j and j = 1,...,ILO-1 or I = IHI+1,...,N.
!          If JOB = 'N' or 'S', ILO = 1 and IHI = N.
!
!  SCALE   (output) DOUBLE PRECISION array, dimension (N)
!          Details of the permutations and scaling factors applied to
!          A.  If P(j) is the index of the row and column interchanged
!          with row and column j and D(j) is the scaling factor
!          applied to row and column j, then
!          SCALE(j) = P(j)    for j = 1,...,ILO-1
!                   = D(j)    for j = ILO,...,IHI
!                   = P(j)    for j = IHI+1,...,N.
!          The order in which the interchanges are made is N to IHI+1,
!          then 1 to ILO-1.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  The permutations consist of row and column interchanges which put
!  the matrix in the form
!
!             ( T1   X   Y  )
!     P A P = (  0   B   Z  )
!             (  0   0   T2 )
!
!  where T1 and T2 are upper triangular matrices whose eigenvalues lie
!  along the diagonal.  The column indices ILO and IHI mark the starting
!  and ending columns of the submatrix B. Balancing consists of applying
!  a diagonal similarity transformation inv(D) * B * D to make the
!  1-norms of each row of B and its corresponding column nearly equal.
!  The output matrix is
!
!     ( T1     X*D          Y    )
!     (  0  inv(D)*B*D  inv(D)*Z ).
!     (  0      0           T2   )
!
!  Information about the permutations P and the diagonal matrix D is
!  returned in the vector SCALE.
!
!  This subroutine is based on the EISPACK routine BALANC.
!
!  Modified by Tzu-Yi Chen, Computer Science Division, University of
!    California at Berkeley, USA
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   SCLFAC
      PARAMETER          ( SCLFAC = 2.0D+0 )
      DOUBLE PRECISION   FACTOR
      PARAMETER          ( FACTOR = 0.95D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOCONV
      INTEGER            I, ICA, IEXC, IRA, J, K, L, M
      DOUBLE PRECISION   C, CA, F, G, R, RA, S, SFMAX1, SFMAX2, SFMIN1, &
                         SFMIN2
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      INTEGER            IDAMAX
!      DOUBLE PRECISION   DLAMCH
!      EXTERNAL           LSAME, IDAMAX, DLAMCH
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DSCAL, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0
      IF( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.LSAME( JOB, 'P' ) .AND. &
          .NOT.LSAME( JOB, 'S' ) .AND. .NOT.LSAME( JOB, 'B' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEBAL', -INFO )
         RETURN
      END IF
!
      K = 1
      L = N
!
      IF( N.EQ.0 ) &
         GO TO 210
!
      IF( LSAME( JOB, 'N' ) ) THEN
         DO 10 I = 1, N
            SCALE( I ) = ONE
   10    CONTINUE
         GO TO 210
      END IF
!
      IF( LSAME( JOB, 'S' ) ) &
         GO TO 120
!
!     Permutation to isolate eigenvalues if possible
!
      GO TO 50
!
!     Row and column exchange.
!
   20 CONTINUE
      SCALE( M ) = J
      IF( J.EQ.M ) &
         GO TO 30
!
      CALL DSWAP( L, A( 1, J ), 1, A( 1, M ), 1 )
      CALL DSWAP( N-K+1, A( J, K ), LDA, A( M, K ), LDA )
!
   30 CONTINUE
      GO TO ( 40, 80 )IEXC
!
!     Search for rows isolating an eigenvalue and push them down.
!
   40 CONTINUE
      IF( L.EQ.1 ) &
         GO TO 210
      L = L - 1
!
   50 CONTINUE
      DO 70 J = L, 1, -1
!
         DO 60 I = 1, L
            IF( I.EQ.J ) &
               GO TO 60
            IF( A( J, I ).NE.ZERO ) &
               GO TO 70
   60    CONTINUE
!
         M = L
         IEXC = 1
         GO TO 20
   70 CONTINUE
!
      GO TO 90
!
!     Search for columns isolating an eigenvalue and push them left.
!
   80 CONTINUE
      K = K + 1
!
   90 CONTINUE
      DO 110 J = K, L
!
         DO 100 I = K, L
            IF( I.EQ.J ) &
               GO TO 100
            IF( A( I, J ).NE.ZERO ) &
               GO TO 110
  100    CONTINUE
!
         M = K
         IEXC = 2
         GO TO 20
  110 CONTINUE
!
  120 CONTINUE
      DO 130 I = K, L
         SCALE( I ) = ONE
  130 CONTINUE
!
      IF( LSAME( JOB, 'P' ) ) &
         GO TO 210
!
!     Balance the submatrix in rows K to L.
!
!     Iterative loop for norm reduction
!
      SFMIN1 = DLAMCH( 'S' ) / DLAMCH( 'P' )
      SFMAX1 = ONE / SFMIN1
      SFMIN2 = SFMIN1*SCLFAC
      SFMAX2 = ONE / SFMIN2
  140 CONTINUE
      NOCONV = .FALSE.
!
      DO 200 I = K, L
         C = ZERO
         R = ZERO
!
         DO 150 J = K, L
            IF( J.EQ.I ) &
               GO TO 150
            C = C + ABS( A( J, I ) )
            R = R + ABS( A( I, J ) )
  150    CONTINUE
         ICA = IDAMAX( L, A( 1, I ), 1 )
         CA = ABS( A( ICA, I ) )
         IRA = IDAMAX( N-K+1, A( I, K ), LDA )
         RA = ABS( A( I, IRA+K-1 ) )
!
!        Guard against zero C or R due to underflow.
!
         IF( C.EQ.ZERO .OR. R.EQ.ZERO ) &
            GO TO 200
         G = R / SCLFAC
         F = ONE
         S = C + R
  160    CONTINUE
         IF( C.GE.G .OR. MAX( F, C, CA ).GE.SFMAX2 .OR. &
             MIN( R, G, RA ).LE.SFMIN2 )GO TO 170
         F = F*SCLFAC
         C = C*SCLFAC
         CA = CA*SCLFAC
         R = R / SCLFAC
         G = G / SCLFAC
         RA = RA / SCLFAC
         GO TO 160
!
  170    CONTINUE
         G = C / SCLFAC
  180    CONTINUE
         IF( G.LT.R .OR. MAX( R, RA ).GE.SFMAX2 .OR. &
             MIN( F, C, G, CA ).LE.SFMIN2 )GO TO 190
         F = F / SCLFAC
         C = C / SCLFAC
         G = G / SCLFAC
         CA = CA / SCLFAC
         R = R*SCLFAC
         RA = RA*SCLFAC
         GO TO 180
!
!        Now balance.
!
  190    CONTINUE
         IF( ( C+R ).GE.FACTOR*S ) &
            GO TO 200
         IF( F.LT.ONE .AND. SCALE( I ).LT.ONE ) THEN
            IF( F*SCALE( I ).LE.SFMIN1 ) &
               GO TO 200
         END IF
         IF( F.GT.ONE .AND. SCALE( I ).GT.ONE ) THEN
            IF( SCALE( I ).GE.SFMAX1 / F ) &
               GO TO 200
         END IF
         G = ONE / F
         SCALE( I ) = SCALE( I )*F
         NOCONV = .TRUE.
!
         CALL DSCAL( N-K+1, G, A( I, K ), LDA )
         CALL DSCAL( L, F, A( 1, I ), 1 )
!
  200 CONTINUE
!
      IF( NOCONV ) &
         GO TO 140
!
  210 CONTINUE
      ILO = K
      IHI = L
!
      RETURN
!
!     End of DGEBAL
!
      END SUBROUTINE DGEBAL


      SUBROUTINE DORGHR( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DORGHR generates a real orthogonal matrix Q which is defined as the
!  product of IHI-ILO elementary reflectors of order N, as returned by
!  DGEHRD:
!
!  Q = H(ilo) H(ilo+1) . . . H(ihi-1).
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix Q. N >= 0.
!
!  ILO     (input) INTEGER
!  IHI     (input) INTEGER
!          ILO and IHI must have the same values as in the previous call
!          of DGEHRD. Q is equal to the unit matrix except in the
!          submatrix Q(ilo+1:ihi,ilo+1:ihi).
!          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the vectors which define the elementary reflectors,
!          as returned by DGEHRD.
!          On exit, the N-by-N orthogonal matrix Q.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A. LDA >= max(1,N).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (N-1)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEHRD.
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,L
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK. LWORK >= IHI-ILO.
!          For optimum performance LWORK >= (IHI-ILO)*NB, where NB is
!          the optimal blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IINFO, J, LWKOPT, NB, NH
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DORGQR, XERBLA
!     ..
!     .. External Functions ..
!      INTEGER            ILAENV
!      EXTERNAL           ILAENV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      NH = IHI - ILO
      LQUERY = ( LWORK.EQ.-1 )
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, NH ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
!
      IF( INFO.EQ.0 ) THEN
         NB = ILAENV( 1, 'DORGQR', ' ', NH, NH, NH, -1 )
         LWKOPT = MAX( 1, NH )*NB
         WORK( 1 ) = LWKOPT
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGHR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
!     Shift the vectors which define the elementary reflectors one
!     column to the right, and set the first ilo and the last n-ihi
!     rows and columns to those of the unit matrix
!
      DO 40 J = IHI, ILO + 1, -1
         DO 10 I = 1, J - 1
            A( I, J ) = ZERO
   10    CONTINUE
         DO 20 I = J + 1, IHI
            A( I, J ) = A( I, J-1 )
   20    CONTINUE
         DO 30 I = IHI + 1, N
            A( I, J ) = ZERO
   30    CONTINUE
   40 CONTINUE
      DO 60 J = 1, ILO
         DO 50 I = 1, N
            A( I, J ) = ZERO
   50    CONTINUE
         A( J, J ) = ONE
   60 CONTINUE
      DO 80 J = IHI + 1, N
         DO 70 I = 1, N
            A( I, J ) = ZERO
   70    CONTINUE
         A( J, J ) = ONE
   80 CONTINUE
!
      IF( NH.GT.0 ) THEN
!
!        Generate Q(ilo+1:ihi,ilo+1:ihi)
!
         CALL DORGQR( NH, NH, NH, A( ILO+1, ILO+1 ), LDA, TAU( ILO ), &
                      WORK, LWORK, IINFO )
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
!
!     End of DORGHR
!
      END SUBROUTINE DORGHR


      SUBROUTINE DLACPY( UPLO, M, N, A, LDA, B, LDB )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDB, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DLACPY copies all or part of a two-dimensional matrix A to another
!  matrix B.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies the part of the matrix A to be copied to B.
!          = 'U':      Upper triangular part
!          = 'L':      Lower triangular part
!          Otherwise:  All of the matrix A
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!          The m by n matrix A.  If UPLO = 'U', only the upper triangle
!          or trapezoid is accessed; if UPLO = 'L', only the lower
!          triangle or trapezoid is accessed.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  B       (output) DOUBLE PRECISION array, dimension (LDB,N)
!          On exit, B = A in the locations specified by UPLO.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,M).
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, J
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 20 J = 1, N
            DO 10 I = 1, MIN( J, M )
               B( I, J ) = A( I, J )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
         DO 40 J = 1, N
            DO 30 I = J, M
               B( I, J ) = A( I, J )
   30       CONTINUE
   40    CONTINUE
      ELSE
         DO 60 J = 1, N
            DO 50 I = 1, M
               B( I, J ) = A( I, J )
   50       CONTINUE
   60    CONTINUE
      END IF
      RETURN
!
!     End of DLACPY
!
      END SUBROUTINE DLACPY


      SUBROUTINE DTRSNA( JOB, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, &
                         LDVR, S, SEP, MM, M, WORK, LDWORK, IWORK, &
                         INFO )
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     Modified to call DLACN2 in place of DLACON, 5 Feb 03, SJH.
!
!     .. Scalar Arguments ..
      CHARACTER          HOWMNY, JOB
      INTEGER            INFO, LDT, LDVL, LDVR, LDWORK, M, MM, N
!     ..
!     .. Array Arguments ..
      LOGICAL            SELECT( * )
      INTEGER            IWORK( * )
      DOUBLE PRECISION   S( * ), SEP( * ), T( LDT, * ), VL( LDVL, * ), &
                         VR( LDVR, * ), WORK( LDWORK, * )
!     ..
!
!  Purpose
!  =======
!
!  DTRSNA estimates reciprocal condition numbers for specified
!  eigenvalues and/or right eigenvectors of a real upper
!  quasi-triangular matrix T (or of any matrix Q*T*Q**T with Q
!  orthogonal).
!
!  T must be in Schur canonical form (as returned by DHSEQR), that is,
!  block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
!  2-by-2 diagonal block has its diagonal elements equal and its
!  off-diagonal elements of opposite sign.
!
!  Arguments
!  =========
!
!  JOB     (input) CHARACTER*1
!          Specifies whether condition numbers are required for
!          eigenvalues (S) or eigenvectors (SEP):
!          = 'E': for eigenvalues only (S);
!          = 'V': for eigenvectors only (SEP);
!          = 'B': for both eigenvalues and eigenvectors (S and SEP).
!
!  HOWMNY  (input) CHARACTER*1
!          = 'A': compute condition numbers for all eigenpairs;
!          = 'S': compute condition numbers for selected eigenpairs
!                 specified by the array SELECT.
!
!  SELECT  (input) LOGICAL array, dimension (N)
!          If HOWMNY = 'S', SELECT specifies the eigenpairs for which
!          condition numbers are required. To select condition numbers
!          for the eigenpair corresponding to a real eigenvalue w(j),
!          SELECT(j) must be set to .TRUE.. To select condition numbers
!          corresponding to a complex conjugate pair of eigenvalues w(j)
!          and w(j+1), either SELECT(j) or SELECT(j+1) or both, must be
!          set to .TRUE..
!          If HOWMNY = 'A', SELECT is not referenced.
!
!  N       (input) INTEGER
!          The order of the matrix T. N >= 0.
!
!  T       (input) DOUBLE PRECISION array, dimension (LDT,N)
!          The upper quasi-triangular matrix T, in Schur canonical form.
!
!  LDT     (input) INTEGER
!          The leading dimension of the array T. LDT >= max(1,N).
!
!  VL      (input) DOUBLE PRECISION array, dimension (LDVL,M)
!          If JOB = 'E' or 'B', VL must contain left eigenvectors of T
!          (or of any Q*T*Q**T with Q orthogonal), corresponding to the
!          eigenpairs specified by HOWMNY and SELECT. The eigenvectors
!          must be stored in consecutive columns of VL, as returned by
!          DHSEIN or DTREVC.
!          If JOB = 'V', VL is not referenced.
!
!  LDVL    (input) INTEGER
!          The leading dimension of the array VL.
!          LDVL >= 1; and if JOB = 'E' or 'B', LDVL >= N.
!
!  VR      (input) DOUBLE PRECISION array, dimension (LDVR,M)
!          If JOB = 'E' or 'B', VR must contain right eigenvectors of T
!          (or of any Q*T*Q**T with Q orthogonal), corresponding to the
!          eigenpairs specified by HOWMNY and SELECT. The eigenvectors
!          must be stored in consecutive columns of VR, as returned by
!          DHSEIN or DTREVC.
!          If JOB = 'V', VR is not referenced.
!
!  LDVR    (input) INTEGER
!          The leading dimension of the array VR.
!          LDVR >= 1; and if JOB = 'E' or 'B', LDVR >= N.
!
!  S       (output) DOUBLE PRECISION array, dimension (MM)
!          If JOB = 'E' or 'B', the reciprocal condition numbers of the
!          selected eigenvalues, stored in consecutive elements of the
!          array. For a complex conjugate pair of eigenvalues two
!          consecutive elements of S are set to the same value. Thus
!          S(j), SEP(j), and the j-th columns of VL and VR all
!          correspond to the same eigenpair (but not in general the
!          j-th eigenpair, unless all eigenpairs are selected).
!          If JOB = 'V', S is not referenced.
!
!  SEP     (output) DOUBLE PRECISION array, dimension (MM)
!          If JOB = 'V' or 'B', the estimated reciprocal condition
!          numbers of the selected eigenvectors, stored in consecutive
!          elements of the array. For a complex eigenvector two
!          consecutive elements of SEP are set to the same value. If
!          the eigenvalues cannot be reordered to compute SEP(j), SEP(j)
!          is set to 0; this can only occur when the true value would be
!          very small anyway.
!          If JOB = 'E', SEP is not referenced.
!
!  MM      (input) INTEGER
!          The number of elements in the arrays S (if JOB = 'E' or 'B')
!           and/or SEP (if JOB = 'V' or 'B'). MM >= M.
!
!  M       (output) INTEGER
!          The number of elements of the arrays S and/or SEP actually
!          used to store the estimated condition numbers.
!          If HOWMNY = 'A', M is set to N.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (LDWORK,N+6)
!          If JOB = 'E', WORK is not referenced.
!
!  LDWORK  (input) INTEGER
!          The leading dimension of the array WORK.
!          LDWORK >= 1; and if JOB = 'V' or 'B', LDWORK >= N.
!
!  IWORK   (workspace) INTEGER array, dimension (2*(N-1))
!          If JOB = 'E', IWORK is not referenced.
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!
!  Further Details
!  ===============
!
!  The reciprocal of the condition number of an eigenvalue lambda is
!  defined as
!
!          S(lambda) = |v'*u| / (norm(u)*norm(v))
!
!  where u and v are the right and left eigenvectors of T corresponding
!  to lambda; v' denotes the conjugate-transpose of v, and norm(u)
!  denotes the Euclidean norm. These reciprocal condition numbers always
!  lie between zero (very badly conditioned) and one (very well
!  conditioned). If n = 1, S(lambda) is defined to be 1.
!
!  An approximate error bound for a computed eigenvalue W(i) is given by
!
!                      EPS * norm(T) / S(i)
!
!  where EPS is the machine precision.
!
!  The reciprocal of the condition number of the right eigenvector u
!  corresponding to lambda is defined as follows. Suppose
!
!              T = ( lambda  c  )
!                  (   0    T22 )
!
!  Then the reciprocal condition number is
!
!          SEP( lambda, T22 ) = sigma-min( T22 - lambda*I )
!
!  where sigma-min denotes the smallest singular value. We approximate
!  the smallest singular value by the reciprocal of an estimate of the
!  one-norm of the inverse of T22 - lambda*I. If n = 1, SEP(1) is
!  defined to be abs(T(1,1)).
!
!  An approximate error bound for a computed right eigenvector VR(i)
!  is given by
!
!                      EPS * norm(T) / SEP(i)
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            PAIR, SOMCON, WANTBH, WANTS, WANTSP
      INTEGER            I, IERR, IFST, ILST, J, K, KASE, KS, N2, NN
      DOUBLE PRECISION   BIGNUM, COND, CS, DELTA, DUMM, EPS, EST, LNRM, &
                         MU, PROD, PROD1, PROD2, RNRM, SCALE, SMLNUM, SN
!     ..
!     .. Local Arrays ..
      INTEGER            ISAVE( 3 )
      DOUBLE PRECISION   DUMMY( 1 )
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      DOUBLE PRECISION   DDOT, DLAMCH, DLAPY2, DNRM2
!      EXTERNAL           LSAME, DDOT, DLAMCH, DLAPY2, DNRM2
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DLACN2, DLACPY, DLAQTR, DTREXC, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
!     ..
!     .. Executable Statements ..
!
!     Decode and test the input parameters
!
      WANTBH = LSAME( JOB, 'B' )
      WANTS = LSAME( JOB, 'E' ) .OR. WANTBH
      WANTSP = LSAME( JOB, 'V' ) .OR. WANTBH
!
      SOMCON = LSAME( HOWMNY, 'S' )
!
      INFO = 0
      IF( .NOT.WANTS .AND. .NOT.WANTSP ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( HOWMNY, 'A' ) .AND. .NOT.SOMCON ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDVL.LT.1 .OR. ( WANTS .AND. LDVL.LT.N ) ) THEN
         INFO = -8
      ELSE IF( LDVR.LT.1 .OR. ( WANTS .AND. LDVR.LT.N ) ) THEN
         INFO = -10
      ELSE
!
!        Set M to the number of eigenpairs for which condition numbers
!        are required, and test MM.
!
         IF( SOMCON ) THEN
            M = 0
            PAIR = .FALSE.
            DO 10 K = 1, N
               IF( PAIR ) THEN
                  PAIR = .FALSE.
               ELSE
                  IF( K.LT.N ) THEN
                     IF( T( K+1, K ).EQ.ZERO ) THEN
                        IF( SELECT( K ) ) &
                           M = M + 1
                     ELSE
                        PAIR = .TRUE.
                        IF( SELECT( K ) .OR. SELECT( K+1 ) ) &
                           M = M + 2
                     END IF
                  ELSE
                     IF( SELECT( N ) ) &
                        M = M + 1
                  END IF
               END IF
   10       CONTINUE
         ELSE
            M = N
         END IF
!
         IF( MM.LT.M ) THEN
            INFO = -13
         ELSE IF( LDWORK.LT.1 .OR. ( WANTSP .AND. LDWORK.LT.N ) ) THEN
            INFO = -16
         END IF
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTRSNA', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
!
      IF( N.EQ.1 ) THEN
         IF( SOMCON ) THEN
            IF( .NOT.SELECT( 1 ) ) &
               RETURN
         END IF
         IF( WANTS ) &
            S( 1 ) = ONE
         IF( WANTSP ) &
            SEP( 1 ) = ABS( T( 1, 1 ) )
         RETURN
      END IF
!
!     Get machine constants
!
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' ) / EPS
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
!
      KS = 0
      PAIR = .FALSE.
      DO 60 K = 1, N
!
!        Determine whether T(k,k) begins a 1-by-1 or 2-by-2 block.
!
         IF( PAIR ) THEN
            PAIR = .FALSE.
            GO TO 60
         ELSE
            IF( K.LT.N ) &
               PAIR = T( K+1, K ).NE.ZERO
         END IF
!
!        Determine whether condition numbers are required for the k-th
!        eigenpair.
!
         IF( SOMCON ) THEN
            IF( PAIR ) THEN
               IF( .NOT.SELECT( K ) .AND. .NOT.SELECT( K+1 ) ) &
                  GO TO 60
            ELSE
               IF( .NOT.SELECT( K ) ) &
                  GO TO 60
            END IF
         END IF
!
         KS = KS + 1
!
         IF( WANTS ) THEN
!
!           Compute the reciprocal condition number of the k-th
!           eigenvalue.
!
            IF( .NOT.PAIR ) THEN
!
!              Real eigenvalue.
!
               PROD = DDOT( N, VR( 1, KS ), 1, VL( 1, KS ), 1 )
               RNRM = DNRM2( N, VR( 1, KS ), 1 )
               LNRM = DNRM2( N, VL( 1, KS ), 1 )
               S( KS ) = ABS( PROD ) / ( RNRM*LNRM )
            ELSE
!
!              Complex eigenvalue.
!
               PROD1 = DDOT( N, VR( 1, KS ), 1, VL( 1, KS ), 1 )
               PROD1 = PROD1 + DDOT( N, VR( 1, KS+1 ), 1, VL( 1, KS+1 ), &
                       1 )
               PROD2 = DDOT( N, VL( 1, KS ), 1, VR( 1, KS+1 ), 1 )
               PROD2 = PROD2 - DDOT( N, VL( 1, KS+1 ), 1, VR( 1, KS ), &
                       1 )
               RNRM = DLAPY2( DNRM2( N, VR( 1, KS ), 1 ), &
                      DNRM2( N, VR( 1, KS+1 ), 1 ) )
               LNRM = DLAPY2( DNRM2( N, VL( 1, KS ), 1 ), &
                      DNRM2( N, VL( 1, KS+1 ), 1 ) )
               COND = DLAPY2( PROD1, PROD2 ) / ( RNRM*LNRM )
               S( KS ) = COND
               S( KS+1 ) = COND
            END IF
         END IF
!
         IF( WANTSP ) THEN
!
!           Estimate the reciprocal condition number of the k-th
!           eigenvector.
!
!           Copy the matrix T to the array WORK and swap the diagonal
!           block beginning at T(k,k) to the (1,1) position.
!
            CALL DLACPY( 'Full', N, N, T, LDT, WORK, LDWORK )
            IFST = K
            ILST = 1
            CALL DTREXC( 'No Q', N, WORK, LDWORK, DUMMY, 1, IFST, ILST, &
                         WORK( 1, N+1 ), IERR )
!
            IF( IERR.EQ.1 .OR. IERR.EQ.2 ) THEN
!
!              Could not swap because blocks not well separated
!
               SCALE = ONE
               EST = BIGNUM
            ELSE
!
!              Reordering successful
!
               IF( WORK( 2, 1 ).EQ.ZERO ) THEN
!
!                 Form C = T22 - lambda*I in WORK(2:N,2:N).
!
                  DO 20 I = 2, N
                     WORK( I, I ) = WORK( I, I ) - WORK( 1, 1 )
   20             CONTINUE
                  N2 = 1
                  NN = N - 1
               ELSE
!
!                 Triangularize the 2 by 2 block by unitary
!                 transformation U = [  cs   i*ss ]
!                                    [ i*ss   cs  ].
!                 such that the (1,1) position of WORK is complex
!                 eigenvalue lambda with positive imaginary part. (2,2)
!                 position of WORK is the complex eigenvalue lambda
!                 with negative imaginary  part.
!
                  MU = SQRT( ABS( WORK( 1, 2 ) ) )* &
                       SQRT( ABS( WORK( 2, 1 ) ) )
                  DELTA = DLAPY2( MU, WORK( 2, 1 ) )
                  CS = MU / DELTA
                  SN = -WORK( 2, 1 ) / DELTA
!
!                 Form
!
!                 C' = WORK(2:N,2:N) + i*[rwork(1) ..... rwork(n-1) ]
!                                        [   mu                     ]
!                                        [         ..               ]
!                                        [             ..           ]
!                                        [                  mu      ]
!                 where C' is conjugate transpose of complex matrix C,
!                 and RWORK is stored starting in the N+1-st column of
!                 WORK.
!
                  DO 30 J = 3, N
                     WORK( 2, J ) = CS*WORK( 2, J )
                     WORK( J, J ) = WORK( J, J ) - WORK( 1, 1 )
   30             CONTINUE
                  WORK( 2, 2 ) = ZERO
!
                  WORK( 1, N+1 ) = TWO*MU
                  DO 40 I = 2, N - 1
                     WORK( I, N+1 ) = SN*WORK( 1, I+1 )
   40             CONTINUE
                  N2 = 2
                  NN = 2*( N-1 )
               END IF
!
!              Estimate norm(inv(C'))
!
               EST = ZERO
               KASE = 0
   50          CONTINUE
               CALL DLACN2( NN, WORK( 1, N+2 ), WORK( 1, N+4 ), IWORK, &
                            EST, KASE, ISAVE )
               IF( KASE.NE.0 ) THEN
                  IF( KASE.EQ.1 ) THEN
                     IF( N2.EQ.1 ) THEN
!
!                       Real eigenvalue: solve C'*x = scale*c.
!
                        CALL DLAQTR( .TRUE., .TRUE., N-1, WORK( 2, 2 ), &
                                     LDWORK, DUMMY, DUMM, SCALE, &
                                     WORK( 1, N+4 ), WORK( 1, N+6 ), &
                                     IERR )
                     ELSE
!
!                       Complex eigenvalue: solve
!                       C'*(p+iq) = scale*(c+id) in real arithmetic.
!
                        CALL DLAQTR( .TRUE., .FALSE., N-1, WORK( 2, 2 ), &
                                     LDWORK, WORK( 1, N+1 ), MU, SCALE, &
                                     WORK( 1, N+4 ), WORK( 1, N+6 ), &
                                     IERR )
                     END IF
                  ELSE
                     IF( N2.EQ.1 ) THEN
!
!                       Real eigenvalue: solve C*x = scale*c.
!
                        CALL DLAQTR( .FALSE., .TRUE., N-1, WORK( 2, 2 ), &
                                     LDWORK, DUMMY, DUMM, SCALE, &
                                     WORK( 1, N+4 ), WORK( 1, N+6 ), &
                                     IERR )
                     ELSE
!
!                       Complex eigenvalue: solve
!                       C*(p+iq) = scale*(c+id) in real arithmetic.
!
                        CALL DLAQTR( .FALSE., .FALSE., N-1, &
                                     WORK( 2, 2 ), LDWORK, &
                                     WORK( 1, N+1 ), MU, SCALE, &
                                     WORK( 1, N+4 ), WORK( 1, N+6 ), &
                                     IERR )
!
                     END IF
                  END IF
!
                  GO TO 50
               END IF
            END IF
!
            SEP( KS ) = SCALE / MAX( EST, SMLNUM )
            IF( PAIR ) &
               SEP( KS+1 ) = SEP( KS )
         END IF
!
         IF( PAIR ) &
            KS = KS + 1
!
   60 CONTINUE
      RETURN
!
!     End of DTRSNA
!
      END SUBROUTINE DTRSNA


      SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DORGQR generates an M-by-N real matrix Q with orthonormal columns,
!  which is defined as the first N columns of a product of K elementary
!  reflectors of order M
!
!        Q  =  H(1) H(2) . . . H(k)
!
!  as returned by DGEQRF.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix Q. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix Q. M >= N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines the
!          matrix Q. N >= K >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the i-th column must contain the vector which
!          defines the elementary reflector H(i), for i = 1,2,...,k, as
!          returned by DGEQRF in the first k columns of its array
!          argument A.
!          On exit, the M-by-N matrix Q.
!
!  LDA     (input) INTEGER
!          The first dimension of the array A. LDA >= max(1,M).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQRF.
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,L
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK. LWORK >= max(1,N).
!          For optimum performance LWORK >= N*NB, where NB is the
!          optimal blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument has an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, KI, KK, L, LDWORK, &
                         LWKOPT, NB, NBMIN, NX
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DLARFB, DLARFT, DORG2R, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. External Functions ..
!      INTEGER            ILAENV
!      EXTERNAL           ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      NB = ILAENV( 1, 'DORGQR', ' ', M, N, K, -1 )
      LWKOPT = MAX( 1, N )*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGQR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
!
!        Determine when to cross over from blocked to unblocked code.
!
         NX = MAX( 0, ILAENV( 3, 'DORGQR', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
!
!              Not enough workspace to use optimal NB:  reduce NB and
!              determine the minimum value of NB.
!
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DORGQR', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
!
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
!
!        Use blocked code after the last block.
!        The first kk columns are handled by the block method.
!
         KI = ( ( K-NX-1 ) / NB )*NB
         KK = MIN( K, KI+NB )
!
!        Set A(1:kk,kk+1:n) to zero.
!
         DO 20 J = KK + 1, N
            DO 10 I = 1, KK
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
!
!     Use unblocked code for the last or only block.
!
      IF( KK.LT.N ) &
         CALL DORG2R( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA, &
                      TAU( KK+1 ), WORK, IINFO )
!
      IF( KK.GT.0 ) THEN
!
!        Use blocked code
!
         DO 50 I = KI + 1, 1, -NB
            IB = MIN( NB, K-I+1 )
            IF( I+IB.LE.N ) THEN
!
!              Form the triangular factor of the block reflector
!              H = H(i) H(i+1) . . . H(i+ib-1)
!
               CALL DLARFT( 'Forward', 'Columnwise', M-I+1, IB, &
                            A( I, I ), LDA, TAU( I ), WORK, LDWORK )
!
!              Apply H to A(i:m,i+ib:n) from the left
!
               CALL DLARFB( 'Left', 'No transpose', 'Forward', &
                            'Columnwise', M-I+1, N-I-IB+1, IB, &
                            A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ), &
                            LDA, WORK( IB+1 ), LDWORK )
            END IF
!
!           Apply H to rows i:m of current block
!
            CALL DORG2R( M-I+1, IB, IB, A( I, I ), LDA, TAU( I ), WORK, &
                         IINFO )
!
!           Set rows 1:i-1 of current block to zero
!
            DO 40 J = I, I + IB - 1
               DO 30 L = 1, I - 1
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
!
      WORK( 1 ) = IWS
      RETURN
!
!     End of DORGQR
!
      END SUBROUTINE DORGQR


      SUBROUTINE DLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, M, N
      DOUBLE PRECISION   ALPHA, BETA
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DLASET initializes an m-by-n matrix A to BETA on the diagonal and
!  ALPHA on the offdiagonals.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies the part of the matrix A to be set.
!          = 'U':      Upper triangular part is set; the strictly lower
!                      triangular part of A is not changed.
!          = 'L':      Lower triangular part is set; the strictly upper
!                      triangular part of A is not changed.
!          Otherwise:  All of the matrix A is set.
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  ALPHA   (input) DOUBLE PRECISION
!          The constant to which the offdiagonal elements are to be set.
!
!  BETA    (input) DOUBLE PRECISION
!          The constant to which the diagonal elements are to be set.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On exit, the leading m-by-n submatrix of A is set as follows:
!
!          if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,
!          if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,
!          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,
!
!          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, J
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
      IF( LSAME( UPLO, 'U' ) ) THEN
!
!        Set the strictly upper triangular or trapezoidal part of the
!        array to ALPHA.
!
         DO 20 J = 2, N
            DO 10 I = 1, MIN( J-1, M )
               A( I, J ) = ALPHA
   10       CONTINUE
   20    CONTINUE
!
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
!
!        Set the strictly lower triangular or trapezoidal part of the
!        array to ALPHA.
!
         DO 40 J = 1, MIN( M, N )
            DO 30 I = J + 1, M
               A( I, J ) = ALPHA
   30       CONTINUE
   40    CONTINUE
!
      ELSE
!
!        Set the leading m-by-n submatrix to ALPHA.
!
         DO 60 J = 1, N
            DO 50 I = 1, M
               A( I, J ) = ALPHA
   50       CONTINUE
   60    CONTINUE
      END IF
!
!     Set the first min(M,N) diagonal elements to BETA.
!
      DO 70 I = 1, MIN( M, N )
         A( I, I ) = BETA
   70 CONTINUE
!
      RETURN
!
!     End of DLASET
!
      END SUBROUTINE DLASET


      SUBROUTINE DLARFX( SIDE, M, N, V, TAU, C, LDC, WORK )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            LDC, M, N
      DOUBLE PRECISION   TAU
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DLARFX applies a real elementary reflector H to a real m by n
!  matrix C, from either the left or the right. H is represented in the
!  form
!
!        H = I - tau * v * v'
!
!  where tau is a real scalar and v is a real vector.
!
!  If tau = 0, then H is taken to be the unit matrix
!
!  This version uses inline code if H has order < 11.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': form  H * C
!          = 'R': form  C * H
!
!  M       (input) INTEGER
!          The number of rows of the matrix C.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C.
!
!  V       (input) DOUBLE PRECISION array, dimension (M) if SIDE = 'L'
!                                     or (N) if SIDE = 'R'
!          The vector v in the representation of H.
!
!  TAU     (input) DOUBLE PRECISION
!          The value tau in the representation of H.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
!          or C * H if SIDE = 'R'.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDA >= (1,M).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension
!                      (N) if SIDE = 'L'
!                      or (M) if SIDE = 'R'
!          WORK is not referenced if H has order < 11.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            J
      DOUBLE PRECISION   SUM, T1, T10, T2, T3, T4, T5, T6, T7, T8, T9, &
                         V1, V10, V2, V3, V4, V5, V6, V7, V8, V9
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DGEMV, DGER
!     ..
!     .. Executable Statements ..
!
      IF( TAU.EQ.ZERO ) &
         RETURN
      IF( LSAME( SIDE, 'L' ) ) THEN
!
!        Form  H * C, where H has order m.
!
         GO TO ( 10, 30, 50, 70, 90, 110, 130, 150, &
                 170, 190 )M
!
!        Code for general M
!
!        w := C'*v
!
         CALL DGEMV( 'Transpose', M, N, ONE, C, LDC, V, 1, ZERO, WORK, &
                     1 )
!
!        C := C - tau * v * w'
!
         CALL DGER( M, N, -TAU, V, 1, WORK, 1, C, LDC )
         GO TO 410
   10    CONTINUE
!
!        Special code for 1 x 1 Householder
!
         T1 = ONE - TAU*V( 1 )*V( 1 )
         DO 20 J = 1, N
            C( 1, J ) = T1*C( 1, J )
   20    CONTINUE
         GO TO 410
   30    CONTINUE
!
!        Special code for 2 x 2 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         DO 40 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
   40    CONTINUE
         GO TO 410
   50    CONTINUE
!
!        Special code for 3 x 3 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         DO 60 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
   60    CONTINUE
         GO TO 410
   70    CONTINUE
!
!        Special code for 4 x 4 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         DO 80 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + &
                  V4*C( 4, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
            C( 4, J ) = C( 4, J ) - SUM*T4
   80    CONTINUE
         GO TO 410
   90    CONTINUE
!
!        Special code for 5 x 5 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         DO 100 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + &
                  V4*C( 4, J ) + V5*C( 5, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
            C( 4, J ) = C( 4, J ) - SUM*T4
            C( 5, J ) = C( 5, J ) - SUM*T5
  100    CONTINUE
         GO TO 410
  110    CONTINUE
!
!        Special code for 6 x 6 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         DO 120 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + &
                  V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
            C( 4, J ) = C( 4, J ) - SUM*T4
            C( 5, J ) = C( 5, J ) - SUM*T5
            C( 6, J ) = C( 6, J ) - SUM*T6
  120    CONTINUE
         GO TO 410
  130    CONTINUE
!
!        Special code for 7 x 7 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         DO 140 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + &
                  V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J ) + &
                  V7*C( 7, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
            C( 4, J ) = C( 4, J ) - SUM*T4
            C( 5, J ) = C( 5, J ) - SUM*T5
            C( 6, J ) = C( 6, J ) - SUM*T6
            C( 7, J ) = C( 7, J ) - SUM*T7
  140    CONTINUE
         GO TO 410
  150    CONTINUE
!
!        Special code for 8 x 8 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         V8 = V( 8 )
         T8 = TAU*V8
         DO 160 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + &
                  V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J ) + &
                  V7*C( 7, J ) + V8*C( 8, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
            C( 4, J ) = C( 4, J ) - SUM*T4
            C( 5, J ) = C( 5, J ) - SUM*T5
            C( 6, J ) = C( 6, J ) - SUM*T6
            C( 7, J ) = C( 7, J ) - SUM*T7
            C( 8, J ) = C( 8, J ) - SUM*T8
  160    CONTINUE
         GO TO 410
  170    CONTINUE
!
!        Special code for 9 x 9 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         V8 = V( 8 )
         T8 = TAU*V8
         V9 = V( 9 )
         T9 = TAU*V9
         DO 180 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + &
                  V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J ) + &
                  V7*C( 7, J ) + V8*C( 8, J ) + V9*C( 9, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
            C( 4, J ) = C( 4, J ) - SUM*T4
            C( 5, J ) = C( 5, J ) - SUM*T5
            C( 6, J ) = C( 6, J ) - SUM*T6
            C( 7, J ) = C( 7, J ) - SUM*T7
            C( 8, J ) = C( 8, J ) - SUM*T8
            C( 9, J ) = C( 9, J ) - SUM*T9
  180    CONTINUE
         GO TO 410
  190    CONTINUE
!
!        Special code for 10 x 10 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         V8 = V( 8 )
         T8 = TAU*V8
         V9 = V( 9 )
         T9 = TAU*V9
         V10 = V( 10 )
         T10 = TAU*V10
         DO 200 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + &
                  V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J ) + &
                  V7*C( 7, J ) + V8*C( 8, J ) + V9*C( 9, J ) + &
                  V10*C( 10, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
            C( 4, J ) = C( 4, J ) - SUM*T4
            C( 5, J ) = C( 5, J ) - SUM*T5
            C( 6, J ) = C( 6, J ) - SUM*T6
            C( 7, J ) = C( 7, J ) - SUM*T7
            C( 8, J ) = C( 8, J ) - SUM*T8
            C( 9, J ) = C( 9, J ) - SUM*T9
            C( 10, J ) = C( 10, J ) - SUM*T10
  200    CONTINUE
         GO TO 410
      ELSE
!
!        Form  C * H, where H has order n.
!
         GO TO ( 210, 230, 250, 270, 290, 310, 330, 350, &
                 370, 390 )N
!
!        Code for general N
!
!        w := C * v
!
         CALL DGEMV( 'No transpose', M, N, ONE, C, LDC, V, 1, ZERO, &
                     WORK, 1 )
!
!        C := C - tau * w * v'
!
         CALL DGER( M, N, -TAU, WORK, 1, V, 1, C, LDC )
         GO TO 410
  210    CONTINUE
!
!        Special code for 1 x 1 Householder
!
         T1 = ONE - TAU*V( 1 )*V( 1 )
         DO 220 J = 1, M
            C( J, 1 ) = T1*C( J, 1 )
  220    CONTINUE
         GO TO 410
  230    CONTINUE
!
!        Special code for 2 x 2 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         DO 240 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
  240    CONTINUE
         GO TO 410
  250    CONTINUE
!
!        Special code for 3 x 3 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         DO 260 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
  260    CONTINUE
         GO TO 410
  270    CONTINUE
!
!        Special code for 4 x 4 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         DO 280 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + &
                  V4*C( J, 4 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
            C( J, 4 ) = C( J, 4 ) - SUM*T4
  280    CONTINUE
         GO TO 410
  290    CONTINUE
!
!        Special code for 5 x 5 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         DO 300 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + &
                  V4*C( J, 4 ) + V5*C( J, 5 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
            C( J, 4 ) = C( J, 4 ) - SUM*T4
            C( J, 5 ) = C( J, 5 ) - SUM*T5
  300    CONTINUE
         GO TO 410
  310    CONTINUE
!
!        Special code for 6 x 6 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         DO 320 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + &
                  V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
            C( J, 4 ) = C( J, 4 ) - SUM*T4
            C( J, 5 ) = C( J, 5 ) - SUM*T5
            C( J, 6 ) = C( J, 6 ) - SUM*T6
  320    CONTINUE
         GO TO 410
  330    CONTINUE
!
!        Special code for 7 x 7 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         DO 340 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + &
                  V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 ) + &
                  V7*C( J, 7 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
            C( J, 4 ) = C( J, 4 ) - SUM*T4
            C( J, 5 ) = C( J, 5 ) - SUM*T5
            C( J, 6 ) = C( J, 6 ) - SUM*T6
            C( J, 7 ) = C( J, 7 ) - SUM*T7
  340    CONTINUE
         GO TO 410
  350    CONTINUE
!
!        Special code for 8 x 8 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         V8 = V( 8 )
         T8 = TAU*V8
         DO 360 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + &
                  V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 ) + &
                  V7*C( J, 7 ) + V8*C( J, 8 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
            C( J, 4 ) = C( J, 4 ) - SUM*T4
            C( J, 5 ) = C( J, 5 ) - SUM*T5
            C( J, 6 ) = C( J, 6 ) - SUM*T6
            C( J, 7 ) = C( J, 7 ) - SUM*T7
            C( J, 8 ) = C( J, 8 ) - SUM*T8
  360    CONTINUE
         GO TO 410
  370    CONTINUE
!
!        Special code for 9 x 9 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         V8 = V( 8 )
         T8 = TAU*V8
         V9 = V( 9 )
         T9 = TAU*V9
         DO 380 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + &
                  V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 ) + &
                  V7*C( J, 7 ) + V8*C( J, 8 ) + V9*C( J, 9 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
            C( J, 4 ) = C( J, 4 ) - SUM*T4
            C( J, 5 ) = C( J, 5 ) - SUM*T5
            C( J, 6 ) = C( J, 6 ) - SUM*T6
            C( J, 7 ) = C( J, 7 ) - SUM*T7
            C( J, 8 ) = C( J, 8 ) - SUM*T8
            C( J, 9 ) = C( J, 9 ) - SUM*T9
  380    CONTINUE
         GO TO 410
  390    CONTINUE
!
!        Special code for 10 x 10 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         V8 = V( 8 )
         T8 = TAU*V8
         V9 = V( 9 )
         T9 = TAU*V9
         V10 = V( 10 )
         T10 = TAU*V10
         DO 400 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + &
                  V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 ) + &
                  V7*C( J, 7 ) + V8*C( J, 8 ) + V9*C( J, 9 ) + &
                  V10*C( J, 10 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
            C( J, 4 ) = C( J, 4 ) - SUM*T4
            C( J, 5 ) = C( J, 5 ) - SUM*T5
            C( J, 6 ) = C( J, 6 ) - SUM*T6
            C( J, 7 ) = C( J, 7 ) - SUM*T7
            C( J, 8 ) = C( J, 8 ) - SUM*T8
            C( J, 9 ) = C( J, 9 ) - SUM*T9
            C( J, 10 ) = C( J, 10 ) - SUM*T10
  400    CONTINUE
         GO TO 410
      END IF
  410 CONTINUE
      RETURN
!
!     End of DLARFX
!
      END SUBROUTINE DLARFX


      SUBROUTINE DLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, &
                         ILOZ, IHIZ, Z, LDZ, INFO )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
      LOGICAL            WANTT, WANTZ
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), WI( * ), WR( * ), Z( LDZ, * )
!     ..
!
!     Purpose
!     =======
!
!     DLAHQR is an auxiliary routine called by DHSEQR to update the
!     eigenvalues and Schur decomposition already computed by DHSEQR, by
!     dealing with the Hessenberg submatrix in rows and columns ILO to
!     IHI.
!
!     Arguments
!     =========
!
!     WANTT   (input) LOGICAL
!          = .TRUE. : the full Schur form T is required;
!          = .FALSE.: only eigenvalues are required.
!
!     WANTZ   (input) LOGICAL
!          = .TRUE. : the matrix of Schur vectors Z is required;
!          = .FALSE.: Schur vectors are not required.
!
!     N       (input) INTEGER
!          The order of the matrix H.  N >= 0.
!
!     ILO     (input) INTEGER
!     IHI     (input) INTEGER
!          It is assumed that H is already upper quasi-triangular in
!          rows and columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless
!          ILO = 1). DLAHQR works primarily with the Hessenberg
!          submatrix in rows and columns ILO to IHI, but applies
!          transformations to all of H if WANTT is .TRUE..
!          1 <= ILO <= max(1,IHI); IHI <= N.
!
!     H       (input/output) DOUBLE PRECISION array, dimension (LDH,N)
!          On entry, the upper Hessenberg matrix H.
!          On exit, if INFO is zero and if WANTT is .TRUE., H is upper
!          quasi-triangular in rows and columns ILO:IHI, with any
!          2-by-2 diagonal blocks in standard form. If INFO is zero
!          and WANTT is .FALSE., the contents of H are unspecified on
!          exit.  The output state of H if INFO is nonzero is given
!          below under the description of INFO.
!
!     LDH     (input) INTEGER
!          The leading dimension of the array H. LDH >= max(1,N).
!
!     WR      (output) DOUBLE PRECISION array, dimension (N)
!     WI      (output) DOUBLE PRECISION array, dimension (N)
!          The real and imaginary parts, respectively, of the computed
!          eigenvalues ILO to IHI are stored in the corresponding
!          elements of WR and WI. If two eigenvalues are computed as a
!          complex conjugate pair, they are stored in consecutive
!          elements of WR and WI, say the i-th and (i+1)th, with
!          WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., the
!          eigenvalues are stored in the same order as on the diagonal
!          of the Schur form returned in H, with WR(i) = H(i,i), and, if
!          H(i:i+1,i:i+1) is a 2-by-2 diagonal block,
!          WI(i) = sqrt(H(i+1,i)*H(i,i+1)) and WI(i+1) = -WI(i).
!
!     ILOZ    (input) INTEGER
!     IHIZ    (input) INTEGER
!          Specify the rows of Z to which transformations must be
!          applied if WANTZ is .TRUE..
!          1 <= ILOZ <= ILO; IHI <= IHIZ <= N.
!
!     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N)
!          If WANTZ is .TRUE., on entry Z must contain the current
!          matrix Z of transformations accumulated by DHSEQR, and on
!          exit Z has been updated; transformations are applied only to
!          the submatrix Z(ILOZ:IHIZ,ILO:IHI).
!          If WANTZ is .FALSE., Z is not referenced.
!
!     LDZ     (input) INTEGER
!          The leading dimension of the array Z. LDZ >= max(1,N).
!
!     INFO    (output) INTEGER
!           =   0: successful exit
!          .GT. 0: If INFO = i, DLAHQR failed to compute all the
!                  eigenvalues ILO to IHI in a total of 30 iterations
!                  per eigenvalue; elements i+1:ihi of WR and WI
!                  contain those eigenvalues which have been
!                  successfully computed.
!
!                  If INFO .GT. 0 and WANTT is .FALSE., then on exit,
!                  the remaining unconverged eigenvalues are the
!                  eigenvalues of the upper Hessenberg matrix rows
!                  and columns ILO thorugh INFO of the final, output
!                  value of H.
!
!                  If INFO .GT. 0 and WANTT is .TRUE., then on exit
!          (*)       (initial value of H)*U  = U*(final value of H)
!                  where U is an orthognal matrix.    The final
!                  value of H is upper Hessenberg and triangular in
!                  rows and columns INFO+1 through IHI.
!
!                  If INFO .GT. 0 and WANTZ is .TRUE., then on exit
!                      (final value of Z)  = (initial value of Z)*U
!                  where U is the orthogonal matrix in (*)
!                  (regardless of the value of WANTT.)
!
!     Further Details
!     ===============
!
!     02-96 Based on modifications by
!     David Day, Sandia National Laboratory, USA
!
!     12-04 Further modifications by
!     Ralph Byers, University of Kansas, USA
!
!       This is a modified version of DLAHQR from LAPACK version 3.0.
!       It is (1) more robust against overflow and underflow and
!       (2) adopts the more conservative Ahues & Tisseur stopping
!       criterion (LAWN 122, 1997).
!
!     =========================================================
!
!     .. Parameters ..
      INTEGER            ITMAX
      PARAMETER          ( ITMAX = 30 )
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0d0, ONE = 1.0d0, TWO = 2.0d0 )
      DOUBLE PRECISION   DAT1, DAT2
      PARAMETER          ( DAT1 = 3.0d0 / 4.0d0, DAT2 = -0.4375d0 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   AA, AB, BA, BB, CS, DET, H11, H12, H21, H21S, &
                         H22, RT1I, RT1R, RT2I, RT2R, RTDISC, S, SAFMAX, &
                         SAFMIN, SMLNUM, SN, SUM, T1, T2, T3, TR, TST, &
                         ULP, V2, V3
      INTEGER            I, I1, I2, ITS, J, K, L, M, NH, NR, NZ
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   V( 3 )
!     ..
!     .. External Functions ..
!      DOUBLE PRECISION   DLAMCH
!      EXTERNAL           DLAMCH
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DCOPY, DLABAD, DLANV2, DLARFG, DROT
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
      INFO = 0
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
      IF( ILO.EQ.IHI ) THEN
         WR( ILO ) = H( ILO, ILO )
         WI( ILO ) = ZERO
         RETURN
      END IF
!
!     ==== clear out the trash ====
      DO 10 J = ILO, IHI - 3
         H( J+2, J ) = ZERO
         H( J+3, J ) = ZERO
   10 CONTINUE
      IF( ILO.LE.IHI-2 ) &
         H( IHI, IHI-2 ) = ZERO
!
      NH = IHI - ILO + 1
      NZ = IHIZ - ILOZ + 1
!
!     Set machine-dependent constants for the stopping criterion.
!
      SAFMIN = DLAMCH( 'SAFE MINIMUM' )
      SAFMAX = ONE / SAFMIN
      CALL DLABAD( SAFMIN, SAFMAX )
      ULP = DLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( DBLE( NH ) / ULP )
!
!     I1 and I2 are the indices of the first row and last column of H
!     to which transformations must be applied. If eigenvalues only are
!     being computed, I1 and I2 are set inside the main loop.
!
      IF( WANTT ) THEN
         I1 = 1
         I2 = N
      END IF
!
!     The main loop begins here. I is the loop index and decreases from
!     IHI to ILO in steps of 1 or 2. Each iteration of the loop works
!     with the active submatrix in rows and columns L to I.
!     Eigenvalues I+1 to IHI have already converged. Either L = ILO or
!     H(L,L-1) is negligible so that the matrix splits.
!
      I = IHI
   20 CONTINUE
      L = ILO
      IF( I.LT.ILO ) &
         GO TO 160
!
!     Perform QR iterations on rows and columns ILO to I until a
!     submatrix of order 1 or 2 splits off at the bottom because a
!     subdiagonal element has become negligible.
!
      DO 140 ITS = 0, ITMAX
!
!        Look for a single small subdiagonal element.
!
         DO 30 K = I, L + 1, -1
            IF( ABS( H( K, K-1 ) ).LE.SMLNUM ) &
               GO TO 40
            TST = ABS( H( K-1, K-1 ) ) + ABS( H( K, K ) )
            IF( TST.EQ.ZERO ) THEN
               IF( K-2.GE.ILO ) &
                  TST = TST + ABS( H( K-1, K-2 ) )
               IF( K+1.LE.IHI ) &
                  TST = TST + ABS( H( K+1, K ) )
            END IF
!           ==== The following is a conservative small subdiagonal
!           .    deflation  criterion due to Ahues & Tisseur (LAWN 122,
!           .    1997). It has better mathematical foundation and
!           .    improves accuracy in some cases.  ====
            IF( ABS( H( K, K-1 ) ).LE.ULP*TST ) THEN
               AB = MAX( ABS( H( K, K-1 ) ), ABS( H( K-1, K ) ) )
               BA = MIN( ABS( H( K, K-1 ) ), ABS( H( K-1, K ) ) )
               AA = MAX( ABS( H( K, K ) ), &
                    ABS( H( K-1, K-1 )-H( K, K ) ) )
               BB = MIN( ABS( H( K, K ) ), &
                    ABS( H( K-1, K-1 )-H( K, K ) ) )
               S = AA + AB
               IF( BA*( AB / S ).LE.MAX( SMLNUM, &
                   ULP*( BB*( AA / S ) ) ) )GO TO 40
            END IF
   30    CONTINUE
   40    CONTINUE
         L = K
         IF( L.GT.ILO ) THEN
!
!           H(L,L-1) is negligible
!
            H( L, L-1 ) = ZERO
         END IF
!
!        Exit from loop if a submatrix of order 1 or 2 has split off.
!
         IF( L.GE.I-1 ) &
            GO TO 150
!
!        Now the active submatrix is in rows and columns L to I. If
!        eigenvalues only are being computed, only the active submatrix
!        need be transformed.
!
         IF( .NOT.WANTT ) THEN
            I1 = L
            I2 = I
         END IF
!
         IF( ITS.EQ.10 .OR. ITS.EQ.20 ) THEN
!
!           Exceptional shift.
!
            H11 = DAT1*S + H( I, I )
            H12 = DAT2*S
            H21 = S
            H22 = H11
         ELSE
!
!           Prepare to use Francis' double shift
!           (i.e. 2nd degree generalized Rayleigh quotient)
!
            H11 = H( I-1, I-1 )
            H21 = H( I, I-1 )
            H12 = H( I-1, I )
            H22 = H( I, I )
         END IF
         S = ABS( H11 ) + ABS( H12 ) + ABS( H21 ) + ABS( H22 )
         IF( S.EQ.ZERO ) THEN
            RT1R = ZERO
            RT1I = ZERO
            RT2R = ZERO
            RT2I = ZERO
         ELSE
            H11 = H11 / S
            H21 = H21 / S
            H12 = H12 / S
            H22 = H22 / S
            TR = ( H11+H22 ) / TWO
            DET = ( H11-TR )*( H22-TR ) - H12*H21
            RTDISC = SQRT( ABS( DET ) )
            IF( DET.GE.ZERO ) THEN
!
!              ==== complex conjugate shifts ====
!
               RT1R = TR*S
               RT2R = RT1R
               RT1I = RTDISC*S
               RT2I = -RT1I
            ELSE
!
!              ==== real shifts (use only one of them)  ====
!
               RT1R = TR + RTDISC
               RT2R = TR - RTDISC
               IF( ABS( RT1R-H22 ).LE.ABS( RT2R-H22 ) ) THEN
                  RT1R = RT1R*S
                  RT2R = RT1R
               ELSE
                  RT2R = RT2R*S
                  RT1R = RT2R
               END IF
               RT1I = ZERO
               RT2I = ZERO
            END IF
         END IF
!
!        Look for two consecutive small subdiagonal elements.
!
         DO 50 M = I - 2, L, -1
!           Determine the effect of starting the double-shift QR
!           iteration at row M, and see if this would make H(M,M-1)
!           negligible.  (The following uses scaling to avoid
!           overflows and most underflows.)
!
            H21S = H( M+1, M )
            S = ABS( H( M, M )-RT2R ) + ABS( RT2I ) + ABS( H21S )
            H21S = H( M+1, M ) / S
            V( 1 ) = H21S*H( M, M+1 ) + ( H( M, M )-RT1R )* &
                     ( ( H( M, M )-RT2R ) / S ) - RT1I*( RT2I / S )
            V( 2 ) = H21S*( H( M, M )+H( M+1, M+1 )-RT1R-RT2R )
            V( 3 ) = H21S*H( M+2, M+1 )
            S = ABS( V( 1 ) ) + ABS( V( 2 ) ) + ABS( V( 3 ) )
            V( 1 ) = V( 1 ) / S
            V( 2 ) = V( 2 ) / S
            V( 3 ) = V( 3 ) / S
            IF( M.EQ.L ) &
               GO TO 60
            IF( ABS( H( M, M-1 ) )*( ABS( V( 2 ) )+ABS( V( 3 ) ) ).LE. &
                ULP*ABS( V( 1 ) )*( ABS( H( M-1, M-1 ) )+ABS( H( M, &
                M ) )+ABS( H( M+1, M+1 ) ) ) )GO TO 60
   50    CONTINUE
   60    CONTINUE
!
!        Double-shift QR step
!
         DO 130 K = M, I - 1
!
!           The first iteration of this loop determines a reflection G
!           from the vector V and applies it from left and right to H,
!           thus creating a nonzero bulge below the subdiagonal.
!
!           Each subsequent iteration determines a reflection G to
!           restore the Hessenberg form in the (K-1)th column, and thus
!           chases the bulge one step toward the bottom of the active
!           submatrix. NR is the order of G.
!
            NR = MIN( 3, I-K+1 )
            IF( K.GT.M ) &
               CALL DCOPY( NR, H( K, K-1 ), 1, V, 1 )
            CALL DLARFG( NR, V( 1 ), V( 2 ), 1, T1 )
            IF( K.GT.M ) THEN
               H( K, K-1 ) = V( 1 )
               H( K+1, K-1 ) = ZERO
               IF( K.LT.I-1 ) &
                  H( K+2, K-1 ) = ZERO
            ELSE IF( M.GT.L ) THEN
               H( K, K-1 ) = -H( K, K-1 )
            END IF
            V2 = V( 2 )
            T2 = T1*V2
            IF( NR.EQ.3 ) THEN
               V3 = V( 3 )
               T3 = T1*V3
!
!              Apply G from the left to transform the rows of the matrix
!              in columns K to I2.
!
               DO 70 J = K, I2
                  SUM = H( K, J ) + V2*H( K+1, J ) + V3*H( K+2, J )
                  H( K, J ) = H( K, J ) - SUM*T1
                  H( K+1, J ) = H( K+1, J ) - SUM*T2
                  H( K+2, J ) = H( K+2, J ) - SUM*T3
   70          CONTINUE
!
!              Apply G from the right to transform the columns of the
!              matrix in rows I1 to min(K+3,I).
!
               DO 80 J = I1, MIN( K+3, I )
                  SUM = H( J, K ) + V2*H( J, K+1 ) + V3*H( J, K+2 )
                  H( J, K ) = H( J, K ) - SUM*T1
                  H( J, K+1 ) = H( J, K+1 ) - SUM*T2
                  H( J, K+2 ) = H( J, K+2 ) - SUM*T3
   80          CONTINUE
!
               IF( WANTZ ) THEN
!
!                 Accumulate transformations in the matrix Z
!
                  DO 90 J = ILOZ, IHIZ
                     SUM = Z( J, K ) + V2*Z( J, K+1 ) + V3*Z( J, K+2 )
                     Z( J, K ) = Z( J, K ) - SUM*T1
                     Z( J, K+1 ) = Z( J, K+1 ) - SUM*T2
                     Z( J, K+2 ) = Z( J, K+2 ) - SUM*T3
   90             CONTINUE
               END IF
            ELSE IF( NR.EQ.2 ) THEN
!
!              Apply G from the left to transform the rows of the matrix
!              in columns K to I2.
!
               DO 100 J = K, I2
                  SUM = H( K, J ) + V2*H( K+1, J )
                  H( K, J ) = H( K, J ) - SUM*T1
                  H( K+1, J ) = H( K+1, J ) - SUM*T2
  100          CONTINUE
!
!              Apply G from the right to transform the columns of the
!              matrix in rows I1 to min(K+3,I).
!
               DO 110 J = I1, I
                  SUM = H( J, K ) + V2*H( J, K+1 )
                  H( J, K ) = H( J, K ) - SUM*T1
                  H( J, K+1 ) = H( J, K+1 ) - SUM*T2
  110          CONTINUE
!
               IF( WANTZ ) THEN
!
!                 Accumulate transformations in the matrix Z
!
                  DO 120 J = ILOZ, IHIZ
                     SUM = Z( J, K ) + V2*Z( J, K+1 )
                     Z( J, K ) = Z( J, K ) - SUM*T1
                     Z( J, K+1 ) = Z( J, K+1 ) - SUM*T2
  120             CONTINUE
               END IF
            END IF
  130    CONTINUE
!
  140 CONTINUE
!
!     Failure to converge in remaining number of iterations
!
      INFO = I
      RETURN
!
  150 CONTINUE
!
      IF( L.EQ.I ) THEN
!
!        H(I,I-1) is negligible: one eigenvalue has converged.
!
         WR( I ) = H( I, I )
         WI( I ) = ZERO
      ELSE IF( L.EQ.I-1 ) THEN
!
!        H(I-1,I-2) is negligible: a pair of eigenvalues have converged.
!
!        Transform the 2-by-2 submatrix to standard Schur form,
!        and compute and store the eigenvalues.
!
         CALL DLANV2( H( I-1, I-1 ), H( I-1, I ), H( I, I-1 ), &
                      H( I, I ), WR( I-1 ), WI( I-1 ), WR( I ), WI( I ), &
                      CS, SN )
!
         IF( WANTT ) THEN
!
!           Apply the transformation to the rest of H.
!
            IF( I2.GT.I ) &
               CALL DROT( I2-I, H( I-1, I+1 ), LDH, H( I, I+1 ), LDH, &
                          CS, SN )
            CALL DROT( I-I1-1, H( I1, I-1 ), 1, H( I1, I ), 1, CS, SN )
         END IF
         IF( WANTZ ) THEN
!
!           Apply the transformation to Z.
!
            CALL DROT( NZ, Z( ILOZ, I-1 ), 1, Z( ILOZ, I ), 1, CS, SN )
         END IF
      END IF
!
!     return to start of the main loop with new value of I.
!
      I = L - 1
      GO TO 20
!
  160 CONTINUE
      RETURN
!
!     End of DLAHQR
!
      END SUBROUTINE DLAHQR


      SUBROUTINE DLAHRD( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            K, LDA, LDT, LDY, N, NB
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), T( LDT, NB ), TAU( NB ), &
                         Y( LDY, NB )
!     ..
!
!  Purpose
!  =======
!
!  DLAHRD reduces the first NB columns of a real general n-by-(n-k+1)
!  matrix A so that elements below the k-th subdiagonal are zero. The
!  reduction is performed by an orthogonal similarity transformation
!  Q' * A * Q. The routine returns the matrices V and T which determine
!  Q as a block reflector I - V*T*V', and also the matrix Y = A * V * T.
!
!  This is an OBSOLETE auxiliary routine.
!  This routine will be 'deprecated' in a  future release.
!  Please use the new routine DLAHR2 instead.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix A.
!
!  K       (input) INTEGER
!          The offset for the reduction. Elements below the k-th
!          subdiagonal in the first NB columns are reduced to zero.
!
!  NB      (input) INTEGER
!          The number of columns to be reduced.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N-K+1)
!          On entry, the n-by-(n-k+1) general matrix A.
!          On exit, the elements on and above the k-th subdiagonal in
!          the first NB columns are overwritten with the corresponding
!          elements of the reduced matrix; the elements below the k-th
!          subdiagonal, with the array TAU, represent the matrix Q as a
!          product of elementary reflectors. The other columns of A are
!          unchanged. See Further Details.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  TAU     (output) DOUBLE PRECISION array, dimension (NB)
!          The scalar factors of the elementary reflectors. See Further
!          Details.
!
!  T       (output) DOUBLE PRECISION array, dimension (LDT,NB)
!          The upper triangular matrix T.
!
!  LDT     (input) INTEGER
!          The leading dimension of the array T.  LDT >= NB.
!
!  Y       (output) DOUBLE PRECISION array, dimension (LDY,NB)
!          The n-by-nb matrix Y.
!
!  LDY     (input) INTEGER
!          The leading dimension of the array Y. LDY >= N.
!
!  Further Details
!  ===============
!
!  The matrix Q is represented as a product of nb elementary reflectors
!
!     Q = H(1) H(2) . . . H(nb).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i+k-1) = 0, v(i+k) = 1; v(i+k+1:n) is stored on exit in
!  A(i+k+1:n,i), and tau in TAU(i).
!
!  The elements of the vectors v together form the (n-k+1)-by-nb matrix
!  V which is needed, with T and Y, to apply the transformation to the
!  unreduced part of the matrix, using an update of the form:
!  A := (I - V*T*V') * (A - Y*V').
!
!  The contents of A on exit are illustrated by the following example
!  with n = 7, k = 3 and nb = 2:
!
!     ( a   h   a   a   a )
!     ( a   h   a   a   a )
!     ( a   h   a   a   a )
!     ( h   h   a   a   a )
!     ( v1  h   a   a   a )
!     ( v1  v2  a   a   a )
!     ( v1  v2  a   a   a )
!
!  where a denotes an element of the original matrix A, h denotes a
!  modified element of the upper Hessenberg matrix H, and vi denotes an
!  element of the vector defining H(i).
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   EI
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DAXPY, DCOPY, DGEMV, DLARFG, DSCAL, DTRMV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( N.LE.1 ) &
         RETURN
!
      DO 10 I = 1, NB
         IF( I.GT.1 ) THEN
!
!           Update A(1:n,i)
!
!           Compute i-th column of A - Y * V'
!
            CALL DGEMV( 'No transpose', N, I-1, -ONE, Y, LDY, &
                        A( K+I-1, 1 ), LDA, ONE, A( 1, I ), 1 )
!
!           Apply I - V * T' * V' to this column (call it b) from the
!           left, using the last column of T as workspace
!
!           Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows)
!                    ( V2 )             ( b2 )
!
!           where V1 is unit lower triangular
!
!           w := V1' * b1
!
            CALL DCOPY( I-1, A( K+1, I ), 1, T( 1, NB ), 1 )
            CALL DTRMV( 'Lower', 'Transpose', 'Unit', I-1, A( K+1, 1 ), &
                        LDA, T( 1, NB ), 1 )
!
!           w := w + V2'*b2
!
            CALL DGEMV( 'Transpose', N-K-I+1, I-1, ONE, A( K+I, 1 ), &
                        LDA, A( K+I, I ), 1, ONE, T( 1, NB ), 1 )
!
!           w := T'*w
!
            CALL DTRMV( 'Upper', 'Transpose', 'Non-unit', I-1, T, LDT, &
                        T( 1, NB ), 1 )
!
!           b2 := b2 - V2*w
!
            CALL DGEMV( 'No transpose', N-K-I+1, I-1, -ONE, A( K+I, 1 ), &
                        LDA, T( 1, NB ), 1, ONE, A( K+I, I ), 1 )
!
!           b1 := b1 - V1*w
!
            CALL DTRMV( 'Lower', 'No transpose', 'Unit', I-1, &
                        A( K+1, 1 ), LDA, T( 1, NB ), 1 )
            CALL DAXPY( I-1, -ONE, T( 1, NB ), 1, A( K+1, I ), 1 )
!
            A( K+I-1, I-1 ) = EI
         END IF
!
!        Generate the elementary reflector H(i) to annihilate
!        A(k+i+1:n,i)
!
         CALL DLARFG( N-K-I+1, A( K+I, I ), A( MIN( K+I+1, N ), I ), 1, &
                      TAU( I ) )
         EI = A( K+I, I )
         A( K+I, I ) = ONE
!
!        Compute  Y(1:n,i)
!
         CALL DGEMV( 'No transpose', N, N-K-I+1, ONE, A( 1, I+1 ), LDA, &
                     A( K+I, I ), 1, ZERO, Y( 1, I ), 1 )
         CALL DGEMV( 'Transpose', N-K-I+1, I-1, ONE, A( K+I, 1 ), LDA, &
                     A( K+I, I ), 1, ZERO, T( 1, I ), 1 )
         CALL DGEMV( 'No transpose', N, I-1, -ONE, Y, LDY, T( 1, I ), 1, &
                     ONE, Y( 1, I ), 1 )
         CALL DSCAL( N, TAU( I ), Y( 1, I ), 1 )
!
!        Compute T(1:i,i)
!
         CALL DSCAL( I-1, -TAU( I ), T( 1, I ), 1 )
         CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1, T, LDT, &
                     T( 1, I ), 1 )
         T( I, I ) = TAU( I )
!
   10 CONTINUE
      A( K+NB, NB ) = EI
!
      RETURN
!
!     End of DLAHRD
!
      END SUBROUTINE DLAHRD


      SUBROUTINE DLAQTR( LTRAN, LREAL, N, T, LDT, B, W, SCALE, X, WORK, &
                         INFO )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      LOGICAL            LREAL, LTRAN
      INTEGER            INFO, LDT, N
      DOUBLE PRECISION   SCALE, W
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   B( * ), T( LDT, * ), WORK( * ), X( * )
!     ..
!
!  Purpose
!  =======
!
!  DLAQTR solves the real quasi-triangular system
!
!               op(T)*p = scale*c,               if LREAL = .TRUE.
!
!  or the complex quasi-triangular systems
!
!             op(T + iB)*(p+iq) = scale*(c+id),  if LREAL = .FALSE.
!
!  in real arithmetic, where T is upper quasi-triangular.
!  If LREAL = .FALSE., then the first diagonal block of T must be
!  1 by 1, B is the specially structured matrix
!
!                 B = [ b(1) b(2) ... b(n) ]
!                     [       w            ]
!                     [           w        ]
!                     [              .     ]
!                     [                 w  ]
!
!  op(A) = A or A', A' denotes the conjugate transpose of
!  matrix A.
!
!  On input, X = [ c ].  On output, X = [ p ].
!                [ d ]                  [ q ]
!
!  This subroutine is designed for the condition number estimation
!  in routine DTRSNA.
!
!  Arguments
!  =========
!
!  LTRAN   (input) LOGICAL
!          On entry, LTRAN specifies the option of conjugate transpose:
!             = .FALSE.,    op(T+i*B) = T+i*B,
!             = .TRUE.,     op(T+i*B) = (T+i*B)'.
!
!  LREAL   (input) LOGICAL
!          On entry, LREAL specifies the input matrix structure:
!             = .FALSE.,    the input is complex
!             = .TRUE.,     the input is real
!
!  N       (input) INTEGER
!          On entry, N specifies the order of T+i*B. N >= 0.
!
!  T       (input) DOUBLE PRECISION array, dimension (LDT,N)
!          On entry, T contains a matrix in Schur canonical form.
!          If LREAL = .FALSE., then the first diagonal block of T mu
!          be 1 by 1.
!
!  LDT     (input) INTEGER
!          The leading dimension of the matrix T. LDT >= max(1,N).
!
!  B       (input) DOUBLE PRECISION array, dimension (N)
!          On entry, B contains the elements to form the matrix
!          B as described above.
!          If LREAL = .TRUE., B is not referenced.
!
!  W       (input) DOUBLE PRECISION
!          On entry, W is the diagonal element of the matrix B.
!          If LREAL = .TRUE., W is not referenced.
!
!  SCALE   (output) DOUBLE PRECISION
!          On exit, SCALE is the scale factor.
!
!  X       (input/output) DOUBLE PRECISION array, dimension (2*N)
!          On entry, X contains the right hand side of the system.
!          On exit, X is overwritten by the solution.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
!
!  INFO    (output) INTEGER
!          On exit, INFO is set to
!             0: successful exit.
!               1: the some diagonal 1 by 1 block has been perturbed by
!                  a small number SMIN to keep nonsingularity.
!               2: the some diagonal 2 by 2 block has been perturbed by
!                  a small number in DLALN2 to keep nonsingularity.
!          NOTE: In the interests of speed, this routine does not
!                check the inputs for errors.
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOTRAN
      INTEGER            I, IERR, J, J1, J2, JNEXT, K, N1, N2
      DOUBLE PRECISION   BIGNUM, EPS, REC, SCALOC, SI, SMIN, SMINW, &
                         SMLNUM, SR, TJJ, TMP, XJ, XMAX, XNORM, Z
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   D( 2, 2 ), V( 2, 2 )
!     ..
!     .. External Functions ..
!      INTEGER            IDAMAX
!      DOUBLE PRECISION   DASUM, DDOT, DLAMCH, DLANGE
!      EXTERNAL           IDAMAX, DASUM, DDOT, DLAMCH, DLANGE
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DAXPY, DLADIV, DLALN2, DSCAL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
!     ..
!     .. Executable Statements ..
!
!     Do not test the input parameters for errors
!
      NOTRAN = .NOT.LTRAN
      INFO = 0
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
!
!     Set constants to control overflow
!
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' ) / EPS
      BIGNUM = ONE / SMLNUM
!
      XNORM = DLANGE( 'M', N, N, T, LDT, D )
      IF( .NOT.LREAL ) &
         XNORM = MAX( XNORM, ABS( W ), DLANGE( 'M', N, 1, B, N, D ) )
      SMIN = MAX( SMLNUM, EPS*XNORM )
!
!     Compute 1-norm of each column of strictly upper triangular
!     part of T to control overflow in triangular solver.
!
      WORK( 1 ) = ZERO
      DO 10 J = 2, N
         WORK( J ) = DASUM( J-1, T( 1, J ), 1 )
   10 CONTINUE
!
      IF( .NOT.LREAL ) THEN
         DO 20 I = 2, N
            WORK( I ) = WORK( I ) + ABS( B( I ) )
   20    CONTINUE
      END IF
!
      N2 = 2*N
      N1 = N
      IF( .NOT.LREAL ) &
         N1 = N2
      K = IDAMAX( N1, X, 1 )
      XMAX = ABS( X( K ) )
      SCALE = ONE
!
      IF( XMAX.GT.BIGNUM ) THEN
         SCALE = BIGNUM / XMAX
         CALL DSCAL( N1, SCALE, X, 1 )
         XMAX = BIGNUM
      END IF
!
      IF( LREAL ) THEN
!
         IF( NOTRAN ) THEN
!
!           Solve T*p = scale*c
!
            JNEXT = N
            DO 30 J = N, 1, -1
               IF( J.GT.JNEXT ) &
                  GO TO 30
               J1 = J
               J2 = J
               JNEXT = J - 1
               IF( J.GT.1 ) THEN
                  IF( T( J, J-1 ).NE.ZERO ) THEN
                     J1 = J - 1
                     JNEXT = J - 2
                  END IF
               END IF
!
               IF( J1.EQ.J2 ) THEN
!
!                 Meet 1 by 1 diagonal block
!
!                 Scale to avoid overflow when computing
!                     x(j) = b(j)/T(j,j)
!
                  XJ = ABS( X( J1 ) )
                  TJJ = ABS( T( J1, J1 ) )
                  TMP = T( J1, J1 )
                  IF( TJJ.LT.SMIN ) THEN
                     TMP = SMIN
                     TJJ = SMIN
                     INFO = 1
                  END IF
!
                  IF( XJ.EQ.ZERO ) &
                     GO TO 30
!
                  IF( TJJ.LT.ONE ) THEN
                     IF( XJ.GT.BIGNUM*TJJ ) THEN
                        REC = ONE / XJ
                        CALL DSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
                  X( J1 ) = X( J1 ) / TMP
                  XJ = ABS( X( J1 ) )
!
!                 Scale x if necessary to avoid overflow when adding a
!                 multiple of column j1 of T.
!
                  IF( XJ.GT.ONE ) THEN
                     REC = ONE / XJ
                     IF( WORK( J1 ).GT.( BIGNUM-XMAX )*REC ) THEN
                        CALL DSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                     END IF
                  END IF
                  IF( J1.GT.1 ) THEN
                     CALL DAXPY( J1-1, -X( J1 ), T( 1, J1 ), 1, X, 1 )
                     K = IDAMAX( J1-1, X, 1 )
                     XMAX = ABS( X( K ) )
                  END IF
!
               ELSE
!
!                 Meet 2 by 2 diagonal block
!
!                 Call 2 by 2 linear system solve, to take
!                 care of possible overflow by scaling factor.
!
                  D( 1, 1 ) = X( J1 )
                  D( 2, 1 ) = X( J2 )
                  CALL DLALN2( .FALSE., 2, 1, SMIN, ONE, T( J1, J1 ), &
                               LDT, ONE, ONE, D, 2, ZERO, ZERO, V, 2, &
                               SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 ) &
                     INFO = 2
!
                  IF( SCALOC.NE.ONE ) THEN
                     CALL DSCAL( N, SCALOC, X, 1 )
                     SCALE = SCALE*SCALOC
                  END IF
                  X( J1 ) = V( 1, 1 )
                  X( J2 ) = V( 2, 1 )
!
!                 Scale V(1,1) (= X(J1)) and/or V(2,1) (=X(J2))
!                 to avoid overflow in updating right-hand side.
!
                  XJ = MAX( ABS( V( 1, 1 ) ), ABS( V( 2, 1 ) ) )
                  IF( XJ.GT.ONE ) THEN
                     REC = ONE / XJ
                     IF( MAX( WORK( J1 ), WORK( J2 ) ).GT. &
                         ( BIGNUM-XMAX )*REC ) THEN
                        CALL DSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                     END IF
                  END IF
!
!                 Update right-hand side
!
                  IF( J1.GT.1 ) THEN
                     CALL DAXPY( J1-1, -X( J1 ), T( 1, J1 ), 1, X, 1 )
                     CALL DAXPY( J1-1, -X( J2 ), T( 1, J2 ), 1, X, 1 )
                     K = IDAMAX( J1-1, X, 1 )
                     XMAX = ABS( X( K ) )
                  END IF
!
               END IF
!
   30       CONTINUE
!
         ELSE
!
!           Solve T'*p = scale*c
!
            JNEXT = 1
            DO 40 J = 1, N
               IF( J.LT.JNEXT ) &
                  GO TO 40
               J1 = J
               J2 = J
               JNEXT = J + 1
               IF( J.LT.N ) THEN
                  IF( T( J+1, J ).NE.ZERO ) THEN
                     J2 = J + 1
                     JNEXT = J + 2
                  END IF
               END IF
!
               IF( J1.EQ.J2 ) THEN
!
!                 1 by 1 diagonal block
!
!                 Scale if necessary to avoid overflow in forming the
!                 right-hand side element by inner product.
!
                  XJ = ABS( X( J1 ) )
                  IF( XMAX.GT.ONE ) THEN
                     REC = ONE / XMAX
                     IF( WORK( J1 ).GT.( BIGNUM-XJ )*REC ) THEN
                        CALL DSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
!
                  X( J1 ) = X( J1 ) - DDOT( J1-1, T( 1, J1 ), 1, X, 1 )
!
                  XJ = ABS( X( J1 ) )
                  TJJ = ABS( T( J1, J1 ) )
                  TMP = T( J1, J1 )
                  IF( TJJ.LT.SMIN ) THEN
                     TMP = SMIN
                     TJJ = SMIN
                     INFO = 1
                  END IF
!
                  IF( TJJ.LT.ONE ) THEN
                     IF( XJ.GT.BIGNUM*TJJ ) THEN
                        REC = ONE / XJ
                        CALL DSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
                  X( J1 ) = X( J1 ) / TMP
                  XMAX = MAX( XMAX, ABS( X( J1 ) ) )
!
               ELSE
!
!                 2 by 2 diagonal block
!
!                 Scale if necessary to avoid overflow in forming the
!                 right-hand side elements by inner product.
!
                  XJ = MAX( ABS( X( J1 ) ), ABS( X( J2 ) ) )
                  IF( XMAX.GT.ONE ) THEN
                     REC = ONE / XMAX
                     IF( MAX( WORK( J2 ), WORK( J1 ) ).GT.( BIGNUM-XJ )* &
                         REC ) THEN
                        CALL DSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
!
                  D( 1, 1 ) = X( J1 ) - DDOT( J1-1, T( 1, J1 ), 1, X, &
                              1 )
                  D( 2, 1 ) = X( J2 ) - DDOT( J1-1, T( 1, J2 ), 1, X, &
                              1 )
!
                  CALL DLALN2( .TRUE., 2, 1, SMIN, ONE, T( J1, J1 ), &
                               LDT, ONE, ONE, D, 2, ZERO, ZERO, V, 2, &
                               SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 ) &
                     INFO = 2
!
                  IF( SCALOC.NE.ONE ) THEN
                     CALL DSCAL( N, SCALOC, X, 1 )
                     SCALE = SCALE*SCALOC
                  END IF
                  X( J1 ) = V( 1, 1 )
                  X( J2 ) = V( 2, 1 )
                  XMAX = MAX( ABS( X( J1 ) ), ABS( X( J2 ) ), XMAX )
!
               END IF
   40       CONTINUE
         END IF
!
      ELSE
!
         SMINW = MAX( EPS*ABS( W ), SMIN )
         IF( NOTRAN ) THEN
!
!           Solve (T + iB)*(p+iq) = c+id
!
            JNEXT = N
            DO 70 J = N, 1, -1
               IF( J.GT.JNEXT ) &
                  GO TO 70
               J1 = J
               J2 = J
               JNEXT = J - 1
               IF( J.GT.1 ) THEN
                  IF( T( J, J-1 ).NE.ZERO ) THEN
                     J1 = J - 1
                     JNEXT = J - 2
                  END IF
               END IF
!
               IF( J1.EQ.J2 ) THEN
!
!                 1 by 1 diagonal block
!
!                 Scale if necessary to avoid overflow in division
!
                  Z = W
                  IF( J1.EQ.1 ) &
                     Z = B( 1 )
                  XJ = ABS( X( J1 ) ) + ABS( X( N+J1 ) )
                  TJJ = ABS( T( J1, J1 ) ) + ABS( Z )
                  TMP = T( J1, J1 )
                  IF( TJJ.LT.SMINW ) THEN
                     TMP = SMINW
                     TJJ = SMINW
                     INFO = 1
                  END IF
!
                  IF( XJ.EQ.ZERO ) &
                     GO TO 70
!
                  IF( TJJ.LT.ONE ) THEN
                     IF( XJ.GT.BIGNUM*TJJ ) THEN
                        REC = ONE / XJ
                        CALL DSCAL( N2, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
                  CALL DLADIV( X( J1 ), X( N+J1 ), TMP, Z, SR, SI )
                  X( J1 ) = SR
                  X( N+J1 ) = SI
                  XJ = ABS( X( J1 ) ) + ABS( X( N+J1 ) )
!
!                 Scale x if necessary to avoid overflow when adding a
!                 multiple of column j1 of T.
!
                  IF( XJ.GT.ONE ) THEN
                     REC = ONE / XJ
                     IF( WORK( J1 ).GT.( BIGNUM-XMAX )*REC ) THEN
                        CALL DSCAL( N2, REC, X, 1 )
                        SCALE = SCALE*REC
                     END IF
                  END IF
!
                  IF( J1.GT.1 ) THEN
                     CALL DAXPY( J1-1, -X( J1 ), T( 1, J1 ), 1, X, 1 )
                     CALL DAXPY( J1-1, -X( N+J1 ), T( 1, J1 ), 1, &
                                 X( N+1 ), 1 )
!
                     X( 1 ) = X( 1 ) + B( J1 )*X( N+J1 )
                     X( N+1 ) = X( N+1 ) - B( J1 )*X( J1 )
!
                     XMAX = ZERO
                     DO 50 K = 1, J1 - 1
                        XMAX = MAX( XMAX, ABS( X( K ) )+ &
                               ABS( X( K+N ) ) )
   50                CONTINUE
                  END IF
!
               ELSE
!
!                 Meet 2 by 2 diagonal block
!
                  D( 1, 1 ) = X( J1 )
                  D( 2, 1 ) = X( J2 )
                  D( 1, 2 ) = X( N+J1 )
                  D( 2, 2 ) = X( N+J2 )
                  CALL DLALN2( .FALSE., 2, 2, SMINW, ONE, T( J1, J1 ), &
                               LDT, ONE, ONE, D, 2, ZERO, -W, V, 2, &
                               SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 ) &
                     INFO = 2
!
                  IF( SCALOC.NE.ONE ) THEN
                     CALL DSCAL( 2*N, SCALOC, X, 1 )
                     SCALE = SCALOC*SCALE
                  END IF
                  X( J1 ) = V( 1, 1 )
                  X( J2 ) = V( 2, 1 )
                  X( N+J1 ) = V( 1, 2 )
                  X( N+J2 ) = V( 2, 2 )
!
!                 Scale X(J1), .... to avoid overflow in
!                 updating right hand side.
!
                  XJ = MAX( ABS( V( 1, 1 ) )+ABS( V( 1, 2 ) ), &
                       ABS( V( 2, 1 ) )+ABS( V( 2, 2 ) ) )
                  IF( XJ.GT.ONE ) THEN
                     REC = ONE / XJ
                     IF( MAX( WORK( J1 ), WORK( J2 ) ).GT. &
                         ( BIGNUM-XMAX )*REC ) THEN
                        CALL DSCAL( N2, REC, X, 1 )
                        SCALE = SCALE*REC
                     END IF
                  END IF
!
!                 Update the right-hand side.
!
                  IF( J1.GT.1 ) THEN
                     CALL DAXPY( J1-1, -X( J1 ), T( 1, J1 ), 1, X, 1 )
                     CALL DAXPY( J1-1, -X( J2 ), T( 1, J2 ), 1, X, 1 )
!
                     CALL DAXPY( J1-1, -X( N+J1 ), T( 1, J1 ), 1, &
                                 X( N+1 ), 1 )
                     CALL DAXPY( J1-1, -X( N+J2 ), T( 1, J2 ), 1, &
                                 X( N+1 ), 1 )
!
                     X( 1 ) = X( 1 ) + B( J1 )*X( N+J1 ) + &
                              B( J2 )*X( N+J2 )
                     X( N+1 ) = X( N+1 ) - B( J1 )*X( J1 ) - &
                                B( J2 )*X( J2 )
!
                     XMAX = ZERO
                     DO 60 K = 1, J1 - 1
                        XMAX = MAX( ABS( X( K ) )+ABS( X( K+N ) ), &
                               XMAX )
   60                CONTINUE
                  END IF
!
               END IF
   70       CONTINUE
!
         ELSE
!
!           Solve (T + iB)'*(p+iq) = c+id
!
            JNEXT = 1
            DO 80 J = 1, N
               IF( J.LT.JNEXT ) &
                  GO TO 80
               J1 = J
               J2 = J
               JNEXT = J + 1
               IF( J.LT.N ) THEN
                  IF( T( J+1, J ).NE.ZERO ) THEN
                     J2 = J + 1
                     JNEXT = J + 2
                  END IF
               END IF
!
               IF( J1.EQ.J2 ) THEN
!
!                 1 by 1 diagonal block
!
!                 Scale if necessary to avoid overflow in forming the
!                 right-hand side element by inner product.
!
                  XJ = ABS( X( J1 ) ) + ABS( X( J1+N ) )
                  IF( XMAX.GT.ONE ) THEN
                     REC = ONE / XMAX
                     IF( WORK( J1 ).GT.( BIGNUM-XJ )*REC ) THEN
                        CALL DSCAL( N2, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
!
                  X( J1 ) = X( J1 ) - DDOT( J1-1, T( 1, J1 ), 1, X, 1 )
                  X( N+J1 ) = X( N+J1 ) - DDOT( J1-1, T( 1, J1 ), 1, &
                              X( N+1 ), 1 )
                  IF( J1.GT.1 ) THEN
                     X( J1 ) = X( J1 ) - B( J1 )*X( N+1 )
                     X( N+J1 ) = X( N+J1 ) + B( J1 )*X( 1 )
                  END IF
                  XJ = ABS( X( J1 ) ) + ABS( X( J1+N ) )
!
                  Z = W
                  IF( J1.EQ.1 ) &
                     Z = B( 1 )
!
!                 Scale if necessary to avoid overflow in
!                 complex division
!
                  TJJ = ABS( T( J1, J1 ) ) + ABS( Z )
                  TMP = T( J1, J1 )
                  IF( TJJ.LT.SMINW ) THEN
                     TMP = SMINW
                     TJJ = SMINW
                     INFO = 1
                  END IF
!
                  IF( TJJ.LT.ONE ) THEN
                     IF( XJ.GT.BIGNUM*TJJ ) THEN
                        REC = ONE / XJ
                        CALL DSCAL( N2, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
                  CALL DLADIV( X( J1 ), X( N+J1 ), TMP, -Z, SR, SI )
                  X( J1 ) = SR
                  X( J1+N ) = SI
                  XMAX = MAX( ABS( X( J1 ) )+ABS( X( J1+N ) ), XMAX )
!
               ELSE
!
!                 2 by 2 diagonal block
!
!                 Scale if necessary to avoid overflow in forming the
!                 right-hand side element by inner product.
!
                  XJ = MAX( ABS( X( J1 ) )+ABS( X( N+J1 ) ), &
                       ABS( X( J2 ) )+ABS( X( N+J2 ) ) )
                  IF( XMAX.GT.ONE ) THEN
                     REC = ONE / XMAX
                     IF( MAX( WORK( J1 ), WORK( J2 ) ).GT. &
                         ( BIGNUM-XJ ) / XMAX ) THEN
                        CALL DSCAL( N2, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
!
                  D( 1, 1 ) = X( J1 ) - DDOT( J1-1, T( 1, J1 ), 1, X, &
                              1 )
                  D( 2, 1 ) = X( J2 ) - DDOT( J1-1, T( 1, J2 ), 1, X, &
                              1 )
                  D( 1, 2 ) = X( N+J1 ) - DDOT( J1-1, T( 1, J1 ), 1, &
                              X( N+1 ), 1 )
                  D( 2, 2 ) = X( N+J2 ) - DDOT( J1-1, T( 1, J2 ), 1, &
                              X( N+1 ), 1 )
                  D( 1, 1 ) = D( 1, 1 ) - B( J1 )*X( N+1 )
                  D( 2, 1 ) = D( 2, 1 ) - B( J2 )*X( N+1 )
                  D( 1, 2 ) = D( 1, 2 ) + B( J1 )*X( 1 )
                  D( 2, 2 ) = D( 2, 2 ) + B( J2 )*X( 1 )
!
                  CALL DLALN2( .TRUE., 2, 2, SMINW, ONE, T( J1, J1 ), &
                               LDT, ONE, ONE, D, 2, ZERO, W, V, 2, &
                               SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 ) &
                     INFO = 2
!
                  IF( SCALOC.NE.ONE ) THEN
                     CALL DSCAL( N2, SCALOC, X, 1 )
                     SCALE = SCALOC*SCALE
                  END IF
                  X( J1 ) = V( 1, 1 )
                  X( J2 ) = V( 2, 1 )
                  X( N+J1 ) = V( 1, 2 )
                  X( N+J2 ) = V( 2, 2 )
                  XMAX = MAX( ABS( X( J1 ) )+ABS( X( N+J1 ) ), &
                         ABS( X( J2 ) )+ABS( X( N+J2 ) ), XMAX )
!
               END IF
!
   80       CONTINUE
!
         END IF
!
      END IF
!
      RETURN
!
!     End of DLAQTR
!
      END SUBROUTINE DLAQTR


      SUBROUTINE DGEHD2( N, ILO, IHI, A, LDA, TAU, WORK, INFO )
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDA, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DGEHD2 reduces a real general matrix A to upper Hessenberg form H by
!  an orthogonal similarity transformation:  Q' * A * Q = H .
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  ILO     (input) INTEGER
!  IHI     (input) INTEGER
!          It is assumed that A is already upper triangular in rows
!          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
!          set by a previous call to DGEBAL; otherwise they should be
!          set to 1 and N respectively. See Further Details.
!          1 <= ILO <= IHI <= max(1,N).
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the n by n general matrix to be reduced.
!          On exit, the upper triangle and the first subdiagonal of A
!          are overwritten with the upper Hessenberg matrix H, and the
!          elements below the first subdiagonal, with the array TAU,
!          represent the orthogonal matrix Q as a product of elementary
!          reflectors. See Further Details.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  TAU     (output) DOUBLE PRECISION array, dimension (N-1)
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  The matrix Q is represented as a product of (ihi-ilo) elementary
!  reflectors
!
!     Q = H(ilo) H(ilo+1) . . . H(ihi-1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on
!  exit in A(i+2:ihi,i), and tau in TAU(i).
!
!  The contents of A are illustrated by the following example, with
!  n = 7, ilo = 2 and ihi = 6:
!
!  on entry,                        on exit,
!
!  ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )
!  (     a   a   a   a   a   a )    (      a   h   h   h   h   a )
!  (     a   a   a   a   a   a )    (      h   h   h   h   h   h )
!  (     a   a   a   a   a   a )    (      v2  h   h   h   h   h )
!  (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )
!  (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )
!  (                         a )    (                          a )
!
!  where a denotes an element of the original matrix A, h denotes a
!  modified element of the upper Hessenberg matrix H, and vi denotes an
!  element of the vector defining H(i).
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   AII
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DLARF, DLARFG, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEHD2', -INFO )
         RETURN
      END IF
!
      DO 10 I = ILO, IHI - 1
!
!        Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)
!
         CALL DLARFG( IHI-I, A( I+1, I ), A( MIN( I+2, N ), I ), 1, &
                      TAU( I ) )
         AII = A( I+1, I )
         A( I+1, I ) = ONE
!
!        Apply H(i) to A(1:ihi,i+1:ihi) from the right
!
         CALL DLARF( 'Right', IHI, IHI-I, A( I+1, I ), 1, TAU( I ), &
                     A( 1, I+1 ), LDA, WORK )
!
!        Apply H(i) to A(i+1:ihi,i+1:n) from the left
!
         CALL DLARF( 'Left', IHI-I, N-I, A( I+1, I ), 1, TAU( I ), &
                     A( I+1, I+1 ), LDA, WORK )
!
         A( I+1, I ) = AII
   10 CONTINUE
!
      RETURN
!
!     End of DGEHD2
!
      END SUBROUTINE DGEHD2


      SUBROUTINE DLALN2( LTRANS, NA, NW, SMIN, CA, A, LDA, D1, D2, B, &
                         LDB, WR, WI, X, LDX, SCALE, XNORM, INFO )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      LOGICAL            LTRANS
      INTEGER            INFO, LDA, LDB, LDX, NA, NW
      DOUBLE PRECISION   CA, D1, D2, SCALE, SMIN, WI, WR, XNORM
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), X( LDX, * )
!     ..
!
!  Purpose
!  =======
!
!  DLALN2 solves a system of the form  (ca A - w D ) X = s B
!  or (ca A' - w D) X = s B   with possible scaling ("s") and
!  perturbation of A.  (A' means A-transpose.)
!
!  A is an NA x NA real matrix, ca is a real scalar, D is an NA x NA
!  real diagonal matrix, w is a real or complex value, and X and B are
!  NA x 1 matrices -- real if w is real, complex if w is complex.  NA
!  may be 1 or 2.
!
!  If w is complex, X and B are represented as NA x 2 matrices,
!  the first column of each being the real part and the second
!  being the imaginary part.
!
!  "s" is a scaling factor (.LE. 1), computed by DLALN2, which is
!  so chosen that X can be computed without overflow.  X is further
!  scaled if necessary to assure that norm(ca A - w D)*norm(X) is less
!  than overflow.
!
!  If both singular values of (ca A - w D) are less than SMIN,
!  SMIN*identity will be used instead of (ca A - w D).  If only one
!  singular value is less than SMIN, one element of (ca A - w D) will be
!  perturbed enough to make the smallest singular value roughly SMIN.
!  If both singular values are at least SMIN, (ca A - w D) will not be
!  perturbed.  In any case, the perturbation will be at most some small
!  multiple of max( SMIN, ulp*norm(ca A - w D) ).  The singular values
!  are computed by infinity-norm approximations, and thus will only be
!  correct to a factor of 2 or so.
!
!  Note: all input quantities are assumed to be smaller than overflow
!  by a reasonable factor.  (See BIGNUM.)
!
!  Arguments
!  ==========
!
!  LTRANS  (input) LOGICAL
!          =.TRUE.:  A-transpose will be used.
!          =.FALSE.: A will be used (not transposed.)
!
!  NA      (input) INTEGER
!          The size of the matrix A.  It may (only) be 1 or 2.
!
!  NW      (input) INTEGER
!          1 if "w" is real, 2 if "w" is complex.  It may only be 1
!          or 2.
!
!  SMIN    (input) DOUBLE PRECISION
!          The desired lower bound on the singular values of A.  This
!          should be a safe distance away from underflow or overflow,
!          say, between (underflow/machine precision) and  (machine
!          precision * overflow ).  (See BIGNUM and ULP.)
!
!  CA      (input) DOUBLE PRECISION
!          The coefficient c, which A is multiplied by.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,NA)
!          The NA x NA matrix A.
!
!  LDA     (input) INTEGER
!          The leading dimension of A.  It must be at least NA.
!
!  D1      (input) DOUBLE PRECISION
!          The 1,1 element in the diagonal matrix D.
!
!  D2      (input) DOUBLE PRECISION
!          The 2,2 element in the diagonal matrix D.  Not used if NW=1.
!
!  B       (input) DOUBLE PRECISION array, dimension (LDB,NW)
!          The NA x NW matrix B (right-hand side).  If NW=2 ("w" is
!          complex), column 1 contains the real part of B and column 2
!          contains the imaginary part.
!
!  LDB     (input) INTEGER
!          The leading dimension of B.  It must be at least NA.
!
!  WR      (input) DOUBLE PRECISION
!          The real part of the scalar "w".
!
!  WI      (input) DOUBLE PRECISION
!          The imaginary part of the scalar "w".  Not used if NW=1.
!
!  X       (output) DOUBLE PRECISION array, dimension (LDX,NW)
!          The NA x NW matrix X (unknowns), as computed by DLALN2.
!          If NW=2 ("w" is complex), on exit, column 1 will contain
!          the real part of X and column 2 will contain the imaginary
!          part.
!
!  LDX     (input) INTEGER
!          The leading dimension of X.  It must be at least NA.
!
!  SCALE   (output) DOUBLE PRECISION
!          The scale factor that B must be multiplied by to insure
!          that overflow does not occur when computing X.  Thus,
!          (ca A - w D) X  will be SCALE*B, not B (ignoring
!          perturbations of A.)  It will be at most 1.
!
!  XNORM   (output) DOUBLE PRECISION
!          The infinity-norm of X, when X is regarded as an NA x NW
!          real matrix.
!
!  INFO    (output) INTEGER
!          An error flag.  It will be set to zero if no error occurs,
!          a negative number if an argument is in error, or a positive
!          number if  ca A - w D  had to be perturbed.
!          The possible values are:
!          = 0: No error occurred, and (ca A - w D) did not have to be
!                 perturbed.
!          = 1: (ca A - w D) had to be perturbed to make its smallest
!               (or only) singular value greater than SMIN.
!          NOTE: In the interests of speed, this routine does not
!                check the inputs for errors.
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
!     ..
!     .. Local Scalars ..
      INTEGER            ICMAX, J
      DOUBLE PRECISION   BBND, BI1, BI2, BIGNUM, BNORM, BR1, BR2, CI21, &
                         CI22, CMAX, CNORM, CR21, CR22, CSI, CSR, LI21, &
                         LR21, SMINI, SMLNUM, TEMP, U22ABS, UI11, UI11R, &
                         UI12, UI12S, UI22, UR11, UR11R, UR12, UR12S, &
                         UR22, XI1, XI2, XR1, XR2
!     ..
!     .. Local Arrays ..
      LOGICAL            RSWAP( 4 ), ZSWAP( 4 )
      INTEGER            IPIVOT( 4, 4 )
      DOUBLE PRECISION   CI( 2, 2 ), CIV( 4 ), CR( 2, 2 ), CRV( 4 )
!     ..
!     .. External Functions ..
!      DOUBLE PRECISION   DLAMCH
!      EXTERNAL           DLAMCH
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DLADIV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
!     ..
!     .. Equivalences ..
      EQUIVALENCE        ( CI( 1, 1 ), CIV( 1 ) ), &
                         ( CR( 1, 1 ), CRV( 1 ) )
!     ..
!     .. Data statements ..
      DATA               ZSWAP / .FALSE., .FALSE., .TRUE., .TRUE. /
      DATA               RSWAP / .FALSE., .TRUE., .FALSE., .TRUE. /
      DATA               IPIVOT / 1, 2, 3, 4, 2, 1, 4, 3, 3, 4, 1, 2, 4, &
                         3, 2, 1 /
!     ..
!     .. Executable Statements ..
!
!     Compute BIGNUM
!
      SMLNUM = TWO*DLAMCH( 'Safe minimum' )
      BIGNUM = ONE / SMLNUM
      SMINI = MAX( SMIN, SMLNUM )
!
!     Don't check for input errors
!
      INFO = 0
!
!     Standard Initializations
!
      SCALE = ONE
!
      IF( NA.EQ.1 ) THEN
!
!        1 x 1  (i.e., scalar) system   C X = B
!
         IF( NW.EQ.1 ) THEN
!
!           Real 1x1 system.
!
!           C = ca A - w D
!
            CSR = CA*A( 1, 1 ) - WR*D1
            CNORM = ABS( CSR )
!
!           If | C | < SMINI, use C = SMINI
!
            IF( CNORM.LT.SMINI ) THEN
               CSR = SMINI
               CNORM = SMINI
               INFO = 1
            END IF
!
!           Check scaling for  X = B / C
!
            BNORM = ABS( B( 1, 1 ) )
            IF( CNORM.LT.ONE .AND. BNORM.GT.ONE ) THEN
               IF( BNORM.GT.BIGNUM*CNORM ) &
                  SCALE = ONE / BNORM
            END IF
!
!           Compute X
!
            X( 1, 1 ) = ( B( 1, 1 )*SCALE ) / CSR
            XNORM = ABS( X( 1, 1 ) )
         ELSE
!
!           Complex 1x1 system (w is complex)
!
!           C = ca A - w D
!
            CSR = CA*A( 1, 1 ) - WR*D1
            CSI = -WI*D1
            CNORM = ABS( CSR ) + ABS( CSI )
!
!           If | C | < SMINI, use C = SMINI
!
            IF( CNORM.LT.SMINI ) THEN
               CSR = SMINI
               CSI = ZERO
               CNORM = SMINI
               INFO = 1
            END IF
!
!           Check scaling for  X = B / C
!
            BNORM = ABS( B( 1, 1 ) ) + ABS( B( 1, 2 ) )
            IF( CNORM.LT.ONE .AND. BNORM.GT.ONE ) THEN
               IF( BNORM.GT.BIGNUM*CNORM ) &
                  SCALE = ONE / BNORM
            END IF
!
!           Compute X
!
            CALL DLADIV( SCALE*B( 1, 1 ), SCALE*B( 1, 2 ), CSR, CSI, &
                         X( 1, 1 ), X( 1, 2 ) )
            XNORM = ABS( X( 1, 1 ) ) + ABS( X( 1, 2 ) )
         END IF
!
      ELSE
!
!        2x2 System
!
!        Compute the real part of  C = ca A - w D  (or  ca A' - w D )
!
         CR( 1, 1 ) = CA*A( 1, 1 ) - WR*D1
         CR( 2, 2 ) = CA*A( 2, 2 ) - WR*D2
         IF( LTRANS ) THEN
            CR( 1, 2 ) = CA*A( 2, 1 )
            CR( 2, 1 ) = CA*A( 1, 2 )
         ELSE
            CR( 2, 1 ) = CA*A( 2, 1 )
            CR( 1, 2 ) = CA*A( 1, 2 )
         END IF
!
         IF( NW.EQ.1 ) THEN
!
!           Real 2x2 system  (w is real)
!
!           Find the largest element in C
!
            CMAX = ZERO
            ICMAX = 0
!
            DO 10 J = 1, 4
               IF( ABS( CRV( J ) ).GT.CMAX ) THEN
                  CMAX = ABS( CRV( J ) )
                  ICMAX = J
               END IF
   10       CONTINUE
!
!           If norm(C) < SMINI, use SMINI*identity.
!
            IF( CMAX.LT.SMINI ) THEN
               BNORM = MAX( ABS( B( 1, 1 ) ), ABS( B( 2, 1 ) ) )
               IF( SMINI.LT.ONE .AND. BNORM.GT.ONE ) THEN
                  IF( BNORM.GT.BIGNUM*SMINI ) &
                     SCALE = ONE / BNORM
               END IF
               TEMP = SCALE / SMINI
               X( 1, 1 ) = TEMP*B( 1, 1 )
               X( 2, 1 ) = TEMP*B( 2, 1 )
               XNORM = TEMP*BNORM
               INFO = 1
               RETURN
            END IF
!
!           Gaussian elimination with complete pivoting.
!
            UR11 = CRV( ICMAX )
            CR21 = CRV( IPIVOT( 2, ICMAX ) )
            UR12 = CRV( IPIVOT( 3, ICMAX ) )
            CR22 = CRV( IPIVOT( 4, ICMAX ) )
            UR11R = ONE / UR11
            LR21 = UR11R*CR21
            UR22 = CR22 - UR12*LR21
!
!           If smaller pivot < SMINI, use SMINI
!
            IF( ABS( UR22 ).LT.SMINI ) THEN
               UR22 = SMINI
               INFO = 1
            END IF
            IF( RSWAP( ICMAX ) ) THEN
               BR1 = B( 2, 1 )
               BR2 = B( 1, 1 )
            ELSE
               BR1 = B( 1, 1 )
               BR2 = B( 2, 1 )
            END IF
            BR2 = BR2 - LR21*BR1
            BBND = MAX( ABS( BR1*( UR22*UR11R ) ), ABS( BR2 ) )
            IF( BBND.GT.ONE .AND. ABS( UR22 ).LT.ONE ) THEN
               IF( BBND.GE.BIGNUM*ABS( UR22 ) ) &
                  SCALE = ONE / BBND
            END IF
!
            XR2 = ( BR2*SCALE ) / UR22
            XR1 = ( SCALE*BR1 )*UR11R - XR2*( UR11R*UR12 )
            IF( ZSWAP( ICMAX ) ) THEN
               X( 1, 1 ) = XR2
               X( 2, 1 ) = XR1
            ELSE
               X( 1, 1 ) = XR1
               X( 2, 1 ) = XR2
            END IF
            XNORM = MAX( ABS( XR1 ), ABS( XR2 ) )
!
!           Further scaling if  norm(A) norm(X) > overflow
!
            IF( XNORM.GT.ONE .AND. CMAX.GT.ONE ) THEN
               IF( XNORM.GT.BIGNUM / CMAX ) THEN
                  TEMP = CMAX / BIGNUM
                  X( 1, 1 ) = TEMP*X( 1, 1 )
                  X( 2, 1 ) = TEMP*X( 2, 1 )
                  XNORM = TEMP*XNORM
                  SCALE = TEMP*SCALE
               END IF
            END IF
         ELSE
!
!           Complex 2x2 system  (w is complex)
!
!           Find the largest element in C
!
            CI( 1, 1 ) = -WI*D1
            CI( 2, 1 ) = ZERO
            CI( 1, 2 ) = ZERO
            CI( 2, 2 ) = -WI*D2
            CMAX = ZERO
            ICMAX = 0
!
            DO 20 J = 1, 4
               IF( ABS( CRV( J ) )+ABS( CIV( J ) ).GT.CMAX ) THEN
                  CMAX = ABS( CRV( J ) ) + ABS( CIV( J ) )
                  ICMAX = J
               END IF
   20       CONTINUE
!
!           If norm(C) < SMINI, use SMINI*identity.
!
            IF( CMAX.LT.SMINI ) THEN
               BNORM = MAX( ABS( B( 1, 1 ) )+ABS( B( 1, 2 ) ), &
                       ABS( B( 2, 1 ) )+ABS( B( 2, 2 ) ) )
               IF( SMINI.LT.ONE .AND. BNORM.GT.ONE ) THEN
                  IF( BNORM.GT.BIGNUM*SMINI ) &
                     SCALE = ONE / BNORM
               END IF
               TEMP = SCALE / SMINI
               X( 1, 1 ) = TEMP*B( 1, 1 )
               X( 2, 1 ) = TEMP*B( 2, 1 )
               X( 1, 2 ) = TEMP*B( 1, 2 )
               X( 2, 2 ) = TEMP*B( 2, 2 )
               XNORM = TEMP*BNORM
               INFO = 1
               RETURN
            END IF
!
!           Gaussian elimination with complete pivoting.
!
            UR11 = CRV( ICMAX )
            UI11 = CIV( ICMAX )
            CR21 = CRV( IPIVOT( 2, ICMAX ) )
            CI21 = CIV( IPIVOT( 2, ICMAX ) )
            UR12 = CRV( IPIVOT( 3, ICMAX ) )
            UI12 = CIV( IPIVOT( 3, ICMAX ) )
            CR22 = CRV( IPIVOT( 4, ICMAX ) )
            CI22 = CIV( IPIVOT( 4, ICMAX ) )
            IF( ICMAX.EQ.1 .OR. ICMAX.EQ.4 ) THEN
!
!              Code when off-diagonals of pivoted C are real
!
               IF( ABS( UR11 ).GT.ABS( UI11 ) ) THEN
                  TEMP = UI11 / UR11
                  UR11R = ONE / ( UR11*( ONE+TEMP**2 ) )
                  UI11R = -TEMP*UR11R
               ELSE
                  TEMP = UR11 / UI11
                  UI11R = -ONE / ( UI11*( ONE+TEMP**2 ) )
                  UR11R = -TEMP*UI11R
               END IF
               LR21 = CR21*UR11R
               LI21 = CR21*UI11R
               UR12S = UR12*UR11R
               UI12S = UR12*UI11R
               UR22 = CR22 - UR12*LR21
               UI22 = CI22 - UR12*LI21
            ELSE
!
!              Code when diagonals of pivoted C are real
!
               UR11R = ONE / UR11
               UI11R = ZERO
               LR21 = CR21*UR11R
               LI21 = CI21*UR11R
               UR12S = UR12*UR11R
               UI12S = UI12*UR11R
               UR22 = CR22 - UR12*LR21 + UI12*LI21
               UI22 = -UR12*LI21 - UI12*LR21
            END IF
            U22ABS = ABS( UR22 ) + ABS( UI22 )
!
!           If smaller pivot < SMINI, use SMINI
!
            IF( U22ABS.LT.SMINI ) THEN
               UR22 = SMINI
               UI22 = ZERO
               INFO = 1
            END IF
            IF( RSWAP( ICMAX ) ) THEN
               BR2 = B( 1, 1 )
               BR1 = B( 2, 1 )
               BI2 = B( 1, 2 )
               BI1 = B( 2, 2 )
            ELSE
               BR1 = B( 1, 1 )
               BR2 = B( 2, 1 )
               BI1 = B( 1, 2 )
               BI2 = B( 2, 2 )
            END IF
            BR2 = BR2 - LR21*BR1 + LI21*BI1
            BI2 = BI2 - LI21*BR1 - LR21*BI1
            BBND = MAX( ( ABS( BR1 )+ABS( BI1 ) )* &
                   ( U22ABS*( ABS( UR11R )+ABS( UI11R ) ) ), &
                   ABS( BR2 )+ABS( BI2 ) )
            IF( BBND.GT.ONE .AND. U22ABS.LT.ONE ) THEN
               IF( BBND.GE.BIGNUM*U22ABS ) THEN
                  SCALE = ONE / BBND
                  BR1 = SCALE*BR1
                  BI1 = SCALE*BI1
                  BR2 = SCALE*BR2
                  BI2 = SCALE*BI2
               END IF
            END IF
!
            CALL DLADIV( BR2, BI2, UR22, UI22, XR2, XI2 )
            XR1 = UR11R*BR1 - UI11R*BI1 - UR12S*XR2 + UI12S*XI2
            XI1 = UI11R*BR1 + UR11R*BI1 - UI12S*XR2 - UR12S*XI2
            IF( ZSWAP( ICMAX ) ) THEN
               X( 1, 1 ) = XR2
               X( 2, 1 ) = XR1
               X( 1, 2 ) = XI2
               X( 2, 2 ) = XI1
            ELSE
               X( 1, 1 ) = XR1
               X( 2, 1 ) = XR2
               X( 1, 2 ) = XI1
               X( 2, 2 ) = XI2
            END IF
            XNORM = MAX( ABS( XR1 )+ABS( XI1 ), ABS( XR2 )+ABS( XI2 ) )
!
!           Further scaling if  norm(A) norm(X) > overflow
!
            IF( XNORM.GT.ONE .AND. CMAX.GT.ONE ) THEN
               IF( XNORM.GT.BIGNUM / CMAX ) THEN
                  TEMP = CMAX / BIGNUM
                  X( 1, 1 ) = TEMP*X( 1, 1 )
                  X( 2, 1 ) = TEMP*X( 2, 1 )
                  X( 1, 2 ) = TEMP*X( 1, 2 )
                  X( 2, 2 ) = TEMP*X( 2, 2 )
                  XNORM = TEMP*XNORM
                  SCALE = TEMP*SCALE
               END IF
            END IF
         END IF
      END IF
!
      RETURN
!
!     End of DLALN2
!
      END SUBROUTINE DLALN2


      SUBROUTINE DTREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, WORK, &
                         INFO )
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          COMPQ
      INTEGER            IFST, ILST, INFO, LDQ, LDT, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   Q( LDQ, * ), T( LDT, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DTREXC reorders the real Schur factorization of a real matrix
!  A = Q*T*Q**T, so that the diagonal block of T with row index IFST is
!  moved to row ILST.
!
!  The real Schur form T is reordered by an orthogonal similarity
!  transformation Z**T*T*Z, and optionally the matrix Q of Schur vectors
!  is updated by postmultiplying it with Z.
!
!  T must be in Schur canonical form (as returned by DHSEQR), that is,
!  block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
!  2-by-2 diagonal block has its diagonal elements equal and its
!  off-diagonal elements of opposite sign.
!
!  Arguments
!  =========
!
!  COMPQ   (input) CHARACTER*1
!          = 'V':  update the matrix Q of Schur vectors;
!          = 'N':  do not update Q.
!
!  N       (input) INTEGER
!          The order of the matrix T. N >= 0.
!
!  T       (input/output) DOUBLE PRECISION array, dimension (LDT,N)
!          On entry, the upper quasi-triangular matrix T, in Schur
!          Schur canonical form.
!          On exit, the reordered upper quasi-triangular matrix, again
!          in Schur canonical form.
!
!  LDT     (input) INTEGER
!          The leading dimension of the array T. LDT >= max(1,N).
!
!  Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
!          On entry, if COMPQ = 'V', the matrix Q of Schur vectors.
!          On exit, if COMPQ = 'V', Q has been postmultiplied by the
!          orthogonal transformation matrix Z which reorders T.
!          If COMPQ = 'N', Q is not referenced.
!
!  LDQ     (input) INTEGER
!          The leading dimension of the array Q.  LDQ >= max(1,N).
!
!  IFST    (input/output) INTEGER
!  ILST    (input/output) INTEGER
!          Specify the reordering of the diagonal blocks of T.
!          The block with row index IFST is moved to row ILST, by a
!          sequence of transpositions between adjacent blocks.
!          On exit, if IFST pointed on entry to the second row of a
!          2-by-2 block, it is changed to point to the first row; ILST
!          always points to the first row of the block in its final
!          position (which may differ from its input value by +1 or -1).
!          1 <= IFST <= N; 1 <= ILST <= N.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          = 1:  two adjacent blocks were too close to swap (the problem
!                is very ill-conditioned); T may have been partially
!                reordered, and ILST points to the first row of the
!                current position of the block being moved.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            WANTQ
      INTEGER            HERE, NBF, NBL, NBNEXT
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DLAEXC, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Decode and test the input arguments.
!
      INFO = 0
      WANTQ = LSAME( COMPQ, 'V' )
      IF( .NOT.WANTQ .AND. .NOT.LSAME( COMPQ, 'N' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDQ.LT.1 .OR. ( WANTQ .AND. LDQ.LT.MAX( 1, N ) ) ) THEN
         INFO = -6
      ELSE IF( IFST.LT.1 .OR. IFST.GT.N ) THEN
         INFO = -7
      ELSE IF( ILST.LT.1 .OR. ILST.GT.N ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTREXC', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.LE.1 ) &
         RETURN
!
!     Determine the first row of specified block
!     and find out it is 1 by 1 or 2 by 2.
!
      IF( IFST.GT.1 ) THEN
         IF( T( IFST, IFST-1 ).NE.ZERO ) &
            IFST = IFST - 1
      END IF
      NBF = 1
      IF( IFST.LT.N ) THEN
         IF( T( IFST+1, IFST ).NE.ZERO ) &
            NBF = 2
      END IF
!
!     Determine the first row of the final block
!     and find out it is 1 by 1 or 2 by 2.
!
      IF( ILST.GT.1 ) THEN
         IF( T( ILST, ILST-1 ).NE.ZERO ) &
            ILST = ILST - 1
      END IF
      NBL = 1
      IF( ILST.LT.N ) THEN
         IF( T( ILST+1, ILST ).NE.ZERO ) &
            NBL = 2
      END IF
!
      IF( IFST.EQ.ILST ) &
         RETURN
!
      IF( IFST.LT.ILST ) THEN
!
!        Update ILST
!
         IF( NBF.EQ.2 .AND. NBL.EQ.1 ) &
            ILST = ILST - 1
         IF( NBF.EQ.1 .AND. NBL.EQ.2 ) &
            ILST = ILST + 1
!
         HERE = IFST
!
   10    CONTINUE
!
!        Swap block with next one below
!
         IF( NBF.EQ.1 .OR. NBF.EQ.2 ) THEN
!
!           Current block either 1 by 1 or 2 by 2
!
            NBNEXT = 1
            IF( HERE+NBF+1.LE.N ) THEN
               IF( T( HERE+NBF+1, HERE+NBF ).NE.ZERO ) &
                  NBNEXT = 2
            END IF
            CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, NBF, NBNEXT, &
                         WORK, INFO )
            IF( INFO.NE.0 ) THEN
               ILST = HERE
               RETURN
            END IF
            HERE = HERE + NBNEXT
!
!           Test if 2 by 2 block breaks into two 1 by 1 blocks
!
            IF( NBF.EQ.2 ) THEN
               IF( T( HERE+1, HERE ).EQ.ZERO ) &
                  NBF = 3
            END IF
!
         ELSE
!
!           Current block consists of two 1 by 1 blocks each of which
!           must be swapped individually
!
            NBNEXT = 1
            IF( HERE+3.LE.N ) THEN
               IF( T( HERE+3, HERE+2 ).NE.ZERO ) &
                  NBNEXT = 2
            END IF
            CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE+1, 1, NBNEXT, &
                         WORK, INFO )
            IF( INFO.NE.0 ) THEN
               ILST = HERE
               RETURN
            END IF
            IF( NBNEXT.EQ.1 ) THEN
!
!              Swap two 1 by 1 blocks, no problems possible
!
               CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, 1, NBNEXT, &
                            WORK, INFO )
               HERE = HERE + 1
            ELSE
!
!              Recompute NBNEXT in case 2 by 2 split
!
               IF( T( HERE+2, HERE+1 ).EQ.ZERO ) &
                  NBNEXT = 1
               IF( NBNEXT.EQ.2 ) THEN
!
!                 2 by 2 Block did not split
!
                  CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, 1, &
                               NBNEXT, WORK, INFO )
                  IF( INFO.NE.0 ) THEN
                     ILST = HERE
                     RETURN
                  END IF
                  HERE = HERE + 2
               ELSE
!
!                 2 by 2 Block did split
!
                  CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, 1, 1, &
                               WORK, INFO )
                  CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE+1, 1, 1, &
                               WORK, INFO )
                  HERE = HERE + 2
               END IF
            END IF
         END IF
         IF( HERE.LT.ILST ) &
            GO TO 10
!
      ELSE
!
         HERE = IFST
   20    CONTINUE
!
!        Swap block with next one above
!
         IF( NBF.EQ.1 .OR. NBF.EQ.2 ) THEN
!
!           Current block either 1 by 1 or 2 by 2
!
            NBNEXT = 1
            IF( HERE.GE.3 ) THEN
               IF( T( HERE-1, HERE-2 ).NE.ZERO ) &
                  NBNEXT = 2
            END IF
            CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE-NBNEXT, NBNEXT, &
                         NBF, WORK, INFO )
            IF( INFO.NE.0 ) THEN
               ILST = HERE
               RETURN
            END IF
            HERE = HERE - NBNEXT
!
!           Test if 2 by 2 block breaks into two 1 by 1 blocks
!
            IF( NBF.EQ.2 ) THEN
               IF( T( HERE+1, HERE ).EQ.ZERO ) &
                  NBF = 3
            END IF
!
         ELSE
!
!           Current block consists of two 1 by 1 blocks each of which
!           must be swapped individually
!
            NBNEXT = 1
            IF( HERE.GE.3 ) THEN
               IF( T( HERE-1, HERE-2 ).NE.ZERO ) &
                  NBNEXT = 2
            END IF
            CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE-NBNEXT, NBNEXT, &
                         1, WORK, INFO )
            IF( INFO.NE.0 ) THEN
               ILST = HERE
               RETURN
            END IF
            IF( NBNEXT.EQ.1 ) THEN
!
!              Swap two 1 by 1 blocks, no problems possible
!
               CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, NBNEXT, 1, &
                            WORK, INFO )
               HERE = HERE - 1
            ELSE
!
!              Recompute NBNEXT in case 2 by 2 split
!
               IF( T( HERE, HERE-1 ).EQ.ZERO ) &
                  NBNEXT = 1
               IF( NBNEXT.EQ.2 ) THEN
!
!                 2 by 2 Block did not split
!
                  CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE-1, 2, 1, &
                               WORK, INFO )
                  IF( INFO.NE.0 ) THEN
                     ILST = HERE
                     RETURN
                  END IF
                  HERE = HERE - 2
               ELSE
!
!                 2 by 2 Block did split
!
                  CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, 1, 1, &
                               WORK, INFO )
                  CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE-1, 1, 1, &
                               WORK, INFO )
                  HERE = HERE - 2
               END IF
            END IF
         END IF
         IF( HERE.GT.ILST ) &
            GO TO 20
      END IF
      ILST = HERE
!
      RETURN
!
!     End of DTREXC
!
      END SUBROUTINE DTREXC


      SUBROUTINE DLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV, &
                         T, LDT, C, LDC, WORK, LDWORK )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          DIRECT, SIDE, STOREV, TRANS
      INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), T( LDT, * ), V( LDV, * ), &
                         WORK( LDWORK, * )
!     ..
!
!  Purpose
!  =======
!
!  DLARFB applies a real block reflector H or its transpose H' to a
!  real m by n matrix C, from either the left or the right.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply H or H' from the Left
!          = 'R': apply H or H' from the Right
!
!  TRANS   (input) CHARACTER*1
!          = 'N': apply H (No transpose)
!          = 'T': apply H' (Transpose)
!
!  DIRECT  (input) CHARACTER*1
!          Indicates how H is formed from a product of elementary
!          reflectors
!          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!
!  STOREV  (input) CHARACTER*1
!          Indicates how the vectors which define the elementary
!          reflectors are stored:
!          = 'C': Columnwise
!          = 'R': Rowwise
!
!  M       (input) INTEGER
!          The number of rows of the matrix C.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C.
!
!  K       (input) INTEGER
!          The order of the matrix T (= the number of elementary
!          reflectors whose product defines the block reflector).
!
!  V       (input) DOUBLE PRECISION array, dimension
!                                (LDV,K) if STOREV = 'C'
!                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
!                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
!          The matrix V. See further details.
!
!  LDV     (input) INTEGER
!          The leading dimension of the array V.
!          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
!          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
!          if STOREV = 'R', LDV >= K.
!
!  T       (input) DOUBLE PRECISION array, dimension (LDT,K)
!          The triangular k by k matrix T in the representation of the
!          block reflector.
!
!  LDT     (input) INTEGER
!          The leading dimension of the array T. LDT >= K.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by H*C or H'*C or C*H or C*H'.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDA >= max(1,M).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (LDWORK,K)
!
!  LDWORK  (input) INTEGER
!          The leading dimension of the array WORK.
!          If SIDE = 'L', LDWORK >= max(1,N);
!          if SIDE = 'R', LDWORK >= max(1,M).
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      CHARACTER          TRANST
      INTEGER            I, J
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DCOPY, DGEMM, DTRMM
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( M.LE.0 .OR. N.LE.0 ) &
         RETURN
!
      IF( LSAME( TRANS, 'N' ) ) THEN
         TRANST = 'T'
      ELSE
         TRANST = 'N'
      END IF
!
      IF( LSAME( STOREV, 'C' ) ) THEN
!
         IF( LSAME( DIRECT, 'F' ) ) THEN
!
!           Let  V =  ( V1 )    (first K rows)
!                     ( V2 )
!           where  V1  is unit lower triangular.
!
            IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
!
!              W := C1'
!
               DO 10 J = 1, K
                  CALL DCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
   10          CONTINUE
!
!              W := W * V1
!
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N, &
                           K, ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
!
!                 W := W + C2'*V2
!
                  CALL DGEMM( 'Transpose', 'No transpose', N, K, M-K, &
                              ONE, C( K+1, 1 ), LDC, V( K+1, 1 ), LDV, &
                              ONE, WORK, LDWORK )
               END IF
!
!              W := W * T'  or  W * T
!
               CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V * W'
!
               IF( M.GT.K ) THEN
!
!                 C2 := C2 - V2 * W'
!
                  CALL DGEMM( 'No transpose', 'Transpose', M-K, N, K, &
                              -ONE, V( K+1, 1 ), LDV, WORK, LDWORK, ONE, &
                              C( K+1, 1 ), LDC )
               END IF
!
!              W := W * V1'
!
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K, &
                           ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W'
!
               DO 30 J = 1, K
                  DO 20 I = 1, N
                     C( J, I ) = C( J, I ) - WORK( I, J )
   20             CONTINUE
   30          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!
!              W := C1
!
               DO 40 J = 1, K
                  CALL DCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
   40          CONTINUE
!
!              W := W * V1
!
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M, &
                           K, ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
!
!                 W := W + C2 * V2
!
                  CALL DGEMM( 'No transpose', 'No transpose', M, K, N-K, &
                              ONE, C( 1, K+1 ), LDC, V( K+1, 1 ), LDV, &
                              ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T'
!
               CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V'
!
               IF( N.GT.K ) THEN
!
!                 C2 := C2 - W * V2'
!
                  CALL DGEMM( 'No transpose', 'Transpose', M, N-K, K, &
                              -ONE, WORK, LDWORK, V( K+1, 1 ), LDV, ONE, &
                              C( 1, K+1 ), LDC )
               END IF
!
!              W := W * V1'
!
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K, &
                           ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
               DO 60 J = 1, K
                  DO 50 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
   50             CONTINUE
   60          CONTINUE
            END IF
!
         ELSE
!
!           Let  V =  ( V1 )
!                     ( V2 )    (last K rows)
!           where  V2  is unit upper triangular.
!
            IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
!
!              W := C2'
!
               DO 70 J = 1, K
                  CALL DCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
   70          CONTINUE
!
!              W := W * V2
!
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N, &
                           K, ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
!
!                 W := W + C1'*V1
!
                  CALL DGEMM( 'Transpose', 'No transpose', N, K, M-K, &
                              ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T'  or  W * T
!
               CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V * W'
!
               IF( M.GT.K ) THEN
!
!                 C1 := C1 - V1 * W'
!
                  CALL DGEMM( 'No transpose', 'Transpose', M-K, N, K, &
                              -ONE, V, LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
!
!              W := W * V2'
!
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K, &
                           ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
!
!              C2 := C2 - W'
!
               DO 90 J = 1, K
                  DO 80 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
   80             CONTINUE
   90          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!
!              W := C2
!
               DO 100 J = 1, K
                  CALL DCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  100          CONTINUE
!
!              W := W * V2
!
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M, &
                           K, ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
!
!                 W := W + C1 * V1
!
                  CALL DGEMM( 'No transpose', 'No transpose', M, K, N-K, &
                              ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T'
!
               CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V'
!
               IF( N.GT.K ) THEN
!
!                 C1 := C1 - W * V1'
!
                  CALL DGEMM( 'No transpose', 'Transpose', M, N-K, K, &
                              -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
               END IF
!
!              W := W * V2'
!
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K, &
                           ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
!
!              C2 := C2 - W
!
               DO 120 J = 1, K
                  DO 110 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  110             CONTINUE
  120          CONTINUE
            END IF
         END IF
!
      ELSE IF( LSAME( STOREV, 'R' ) ) THEN
!
         IF( LSAME( DIRECT, 'F' ) ) THEN
!
!           Let  V =  ( V1  V2 )    (V1: first K columns)
!           where  V1  is unit upper triangular.
!
            IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
!
!              W := C1'
!
               DO 130 J = 1, K
                  CALL DCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
  130          CONTINUE
!
!              W := W * V1'
!
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K, &
                           ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
!
!                 W := W + C2'*V2'
!
                  CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE, &
                              C( K+1, 1 ), LDC, V( 1, K+1 ), LDV, ONE, &
                              WORK, LDWORK )
               END IF
!
!              W := W * T'  or  W * T
!
               CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V' * W'
!
               IF( M.GT.K ) THEN
!
!                 C2 := C2 - V2' * W'
!
                  CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE, &
                              V( 1, K+1 ), LDV, WORK, LDWORK, ONE, &
                              C( K+1, 1 ), LDC )
               END IF
!
!              W := W * V1
!
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N, &
                           K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W'
!
               DO 150 J = 1, K
                  DO 140 I = 1, N
                     C( J, I ) = C( J, I ) - WORK( I, J )
  140             CONTINUE
  150          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
!
!              W := C1
!
               DO 160 J = 1, K
                  CALL DCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
  160          CONTINUE
!
!              W := W * V1'
!
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K, &
                           ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
!
!                 W := W + C2 * V2'
!
                  CALL DGEMM( 'No transpose', 'Transpose', M, K, N-K, &
                              ONE, C( 1, K+1 ), LDC, V( 1, K+1 ), LDV, &
                              ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T'
!
               CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V
!
               IF( N.GT.K ) THEN
!
!                 C2 := C2 - W * V2
!
                  CALL DGEMM( 'No transpose', 'No transpose', M, N-K, K, &
                              -ONE, WORK, LDWORK, V( 1, K+1 ), LDV, ONE, &
                              C( 1, K+1 ), LDC )
               END IF
!
!              W := W * V1
!
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M, &
                           K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
               DO 180 J = 1, K
                  DO 170 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
  170             CONTINUE
  180          CONTINUE
!
            END IF
!
         ELSE
!
!           Let  V =  ( V1  V2 )    (V2: last K columns)
!           where  V2  is unit lower triangular.
!
            IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
!
!              W := C2'
!
               DO 190 J = 1, K
                  CALL DCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
  190          CONTINUE
!
!              W := W * V2'
!
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K, &
                           ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
!
!                 W := W + C1'*V1'
!
                  CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE, &
                              C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T'  or  W * T
!
               CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V' * W'
!
               IF( M.GT.K ) THEN
!
!                 C1 := C1 - V1' * W'
!
                  CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE, &
                              V, LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
!
!              W := W * V2
!
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N, &
                           K, ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
!
!              C2 := C2 - W'
!
               DO 210 J = 1, K
                  DO 200 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
  200             CONTINUE
  210          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
!
!              W := C2
!
               DO 220 J = 1, K
                  CALL DCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  220          CONTINUE
!
!              W := W * V2'
!
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K, &
                           ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
!
!                 W := W + C1 * V1'
!
                  CALL DGEMM( 'No transpose', 'Transpose', M, K, N-K, &
                              ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T'
!
               CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V
!
               IF( N.GT.K ) THEN
!
!                 C1 := C1 - W * V1
!
                  CALL DGEMM( 'No transpose', 'No transpose', M, N-K, K, &
                              -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
               END IF
!
!              W := W * V2
!
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M, &
                           K, ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
               DO 240 J = 1, K
                  DO 230 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  230             CONTINUE
  240          CONTINUE
!
            END IF
!
         END IF
      END IF
!
      RETURN
!
!     End of DLARFB
!
      END SUBROUTINE DLARFB


      SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   ALPHA, TAU
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
!     ..
!
!  Purpose
!  =======
!
!  DLARFG generates a real elementary reflector H of order n, such
!  that
!
!        H * ( alpha ) = ( beta ),   H' * H = I.
!            (   x   )   (   0  )
!
!  where alpha and beta are scalars, and x is an (n-1)-element real
!  vector. H is represented in the form
!
!        H = I - tau * ( 1 ) * ( 1 v' ) ,
!                      ( v )
!
!  where tau is a real scalar and v is a real (n-1)-element
!  vector.
!
!  If the elements of x are all zero, then tau = 0 and H is taken to be
!  the unit matrix.
!
!  Otherwise  1 <= tau <= 2.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the elementary reflector.
!
!  ALPHA   (input/output) DOUBLE PRECISION
!          On entry, the value alpha.
!          On exit, it is overwritten with the value beta.
!
!  X       (input/output) DOUBLE PRECISION array, dimension
!                         (1+(N-2)*abs(INCX))
!          On entry, the vector x.
!          On exit, it is overwritten with the vector v.
!
!  INCX    (input) INTEGER
!          The increment between elements of X. INCX > 0.
!
!  TAU     (output) DOUBLE PRECISION
!          The value tau.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            J, KNT
      DOUBLE PRECISION   BETA, RSAFMN, SAFMIN, XNORM
!     ..
!     .. External Functions ..
!      DOUBLE PRECISION   DLAMCH, DLAPY2, DNRM2
!      EXTERNAL           DLAMCH, DLAPY2, DNRM2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DSCAL
!     ..
!     .. Executable Statements ..
!
      IF( N.LE.1 ) THEN
         TAU = ZERO
         RETURN
      END IF
!
      XNORM = DNRM2( N-1, X, INCX )
!
      IF( XNORM.EQ.ZERO ) THEN
!
!        H  =  I
!
         TAU = ZERO
      ELSE
!
!        general case
!
         BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
         SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )
         IF( ABS( BETA ).LT.SAFMIN ) THEN
!
!           XNORM, BETA may be inaccurate; scale X and recompute them
!
            RSAFMN = ONE / SAFMIN
            KNT = 0
   10       CONTINUE
            KNT = KNT + 1
            CALL DSCAL( N-1, RSAFMN, X, INCX )
            BETA = BETA*RSAFMN
            ALPHA = ALPHA*RSAFMN
            IF( ABS( BETA ).LT.SAFMIN ) &
               GO TO 10
!
!           New BETA is at most 1, at least SAFMIN
!
            XNORM = DNRM2( N-1, X, INCX )
            BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
            TAU = ( BETA-ALPHA ) / BETA
            CALL DSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
!
!           If ALPHA is subnormal, it may lose relative accuracy
!
            ALPHA = BETA
            DO 20 J = 1, KNT
               ALPHA = ALPHA*SAFMIN
   20       CONTINUE
         ELSE
            TAU = ( BETA-ALPHA ) / BETA
            CALL DSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
            ALPHA = BETA
         END IF
      END IF
!
      RETURN
!
!     End of DLARFG
!
      END SUBROUTINE DLARFG

      SUBROUTINE DORG2R( M, N, K, A, LDA, TAU, WORK, INFO )
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DORG2R generates an m by n real matrix Q with orthonormal columns,
!  which is defined as the first n columns of a product of k elementary
!  reflectors of order m
!
!        Q  =  H(1) H(2) . . . H(k)
!
!  as returned by DGEQRF.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix Q. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix Q. M >= N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines the
!          matrix Q. N >= K >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the i-th column must contain the vector which
!          defines the elementary reflector H(i), for i = 1,2,...,k, as
!          returned by DGEQRF in the first k columns of its array
!          argument A.
!          On exit, the m-by-n matrix Q.
!
!  LDA     (input) INTEGER
!          The first dimension of the array A. LDA >= max(1,M).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQRF.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument has an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J, L
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DLARF, DSCAL, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORG2R', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.LE.0 ) &
         RETURN
!
!     Initialise columns k+1:n to columns of the unit matrix
!
      DO 20 J = K + 1, N
         DO 10 L = 1, M
            A( L, J ) = ZERO
   10    CONTINUE
         A( J, J ) = ONE
   20 CONTINUE
!
      DO 40 I = K, 1, -1
!
!        Apply H(i) to A(i:m,i:n) from the left
!
         IF( I.LT.N ) THEN
            A( I, I ) = ONE
            CALL DLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ), &
                        A( I, I+1 ), LDA, WORK )
         END IF
         IF( I.LT.M ) &
            CALL DSCAL( M-I, -TAU( I ), A( I+1, I ), 1 )
         A( I, I ) = ONE - TAU( I )
!
!        Set A(1:i-1,i) to zero
!
         DO 30 L = 1, I - 1
            A( L, I ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
!
!     End of DORG2R
!
      END SUBROUTINE DORG2R


      SUBROUTINE DLADIV( A, B, C, D, P, Q )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, D, P, Q
!     ..
!
!  Purpose
!  =======
!
!  DLADIV performs complex division in  real arithmetic
!
!                        a + i*b
!             p + i*q = ---------
!                        c + i*d
!
!  The algorithm is due to Robert L. Smith and can be found
!  in D. Knuth, The art of Computer Programming, Vol.2, p.195
!
!  Arguments
!  =========
!
!  A       (input) DOUBLE PRECISION
!  B       (input) DOUBLE PRECISION
!  C       (input) DOUBLE PRECISION
!  D       (input) DOUBLE PRECISION
!          The scalars a, b, c, and d in the above expression.
!
!  P       (output) DOUBLE PRECISION
!  Q       (output) DOUBLE PRECISION
!          The scalars p and q in the above expression.
!
!  =====================================================================
!
!     .. Local Scalars ..
      DOUBLE PRECISION   E, F
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS
!     ..
!     .. Executable Statements ..
!
      IF( ABS( D ).LT.ABS( C ) ) THEN
         E = D / C
         F = C + D*E
         P = ( A+B*E ) / F
         Q = ( B-A*E ) / F
      ELSE
         E = C / D
         F = D + C*E
         P = ( B+A*E ) / F
         Q = ( -A+B*E ) / F
      END IF
!
      RETURN
!
!     End of DLADIV
!
      END SUBROUTINE DLADIV


      SUBROUTINE DLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          DIRECT, STOREV
      INTEGER            K, LDT, LDV, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   T( LDT, * ), TAU( * ), V( LDV, * )
!     ..
!
!  Purpose
!  =======
!
!  DLARFT forms the triangular factor T of a real block reflector H
!  of order n, which is defined as a product of k elementary reflectors.
!
!  If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
!
!  If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
!
!  If STOREV = 'C', the vector which defines the elementary reflector
!  H(i) is stored in the i-th column of the array V, and
!
!     H  =  I - V * T * V'
!
!  If STOREV = 'R', the vector which defines the elementary reflector
!  H(i) is stored in the i-th row of the array V, and
!
!     H  =  I - V' * T * V
!
!  Arguments
!  =========
!
!  DIRECT  (input) CHARACTER*1
!          Specifies the order in which the elementary reflectors are
!          multiplied to form the block reflector:
!          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!
!  STOREV  (input) CHARACTER*1
!          Specifies how the vectors which define the elementary
!          reflectors are stored (see also Further Details):
!          = 'C': columnwise
!          = 'R': rowwise
!
!  N       (input) INTEGER
!          The order of the block reflector H. N >= 0.
!
!  K       (input) INTEGER
!          The order of the triangular factor T (= the number of
!          elementary reflectors). K >= 1.
!
!  V       (input/output) DOUBLE PRECISION array, dimension
!                               (LDV,K) if STOREV = 'C'
!                               (LDV,N) if STOREV = 'R'
!          The matrix V. See further details.
!
!  LDV     (input) INTEGER
!          The leading dimension of the array V.
!          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i).
!
!  T       (output) DOUBLE PRECISION array, dimension (LDT,K)
!          The k by k triangular factor T of the block reflector.
!          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
!          lower triangular. The rest of the array is not used.
!
!  LDT     (input) INTEGER
!          The leading dimension of the array T. LDT >= K.
!
!  Further Details
!  ===============
!
!  The shape of the matrix V and the storage of the vectors which define
!  the H(i) is best illustrated by the following example with n = 5 and
!  k = 3. The elements equal to 1 are not stored; the corresponding
!  array elements are modified but restored on exit. The rest of the
!  array is not used.
!
!  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
!
!               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
!                   ( v1  1    )                     (     1 v2 v2 v2 )
!                   ( v1 v2  1 )                     (        1 v3 v3 )
!                   ( v1 v2 v3 )
!                   ( v1 v2 v3 )
!
!  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
!
!               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
!                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
!                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
!                   (     1 v3 )
!                   (        1 )
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   VII
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DGEMV, DTRMV
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           LSAME
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
!
      IF( LSAME( DIRECT, 'F' ) ) THEN
         DO 20 I = 1, K
            IF( TAU( I ).EQ.ZERO ) THEN
!
!              H(i)  =  I
!
               DO 10 J = 1, I
                  T( J, I ) = ZERO
   10          CONTINUE
            ELSE
!
!              general case
!
               VII = V( I, I )
               V( I, I ) = ONE
               IF( LSAME( STOREV, 'C' ) ) THEN
!
!                 T(1:i-1,i) := - tau(i) * V(i:n,1:i-1)' * V(i:n,i)
!
                  CALL DGEMV( 'Transpose', N-I+1, I-1, -TAU( I ), &
                              V( I, 1 ), LDV, V( I, I ), 1, ZERO, &
                              T( 1, I ), 1 )
               ELSE
!
!                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:n) * V(i,i:n)'
!
                  CALL DGEMV( 'No transpose', I-1, N-I+1, -TAU( I ), &
                              V( 1, I ), LDV, V( I, I ), LDV, ZERO, &
                              T( 1, I ), 1 )
               END IF
               V( I, I ) = VII
!
!              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
!
               CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1, T, &
                           LDT, T( 1, I ), 1 )
               T( I, I ) = TAU( I )
            END IF
   20    CONTINUE
      ELSE
         DO 40 I = K, 1, -1
            IF( TAU( I ).EQ.ZERO ) THEN
!
!              H(i)  =  I
!
               DO 30 J = I, K
                  T( J, I ) = ZERO
   30          CONTINUE
            ELSE
!
!              general case
!
               IF( I.LT.K ) THEN
                  IF( LSAME( STOREV, 'C' ) ) THEN
                     VII = V( N-K+I, I )
                     V( N-K+I, I ) = ONE
!
!                    T(i+1:k,i) :=
!                            - tau(i) * V(1:n-k+i,i+1:k)' * V(1:n-k+i,i)
!
                     CALL DGEMV( 'Transpose', N-K+I, K-I, -TAU( I ), &
                                 V( 1, I+1 ), LDV, V( 1, I ), 1, ZERO, &
                                 T( I+1, I ), 1 )
                     V( N-K+I, I ) = VII
                  ELSE
                     VII = V( I, N-K+I )
                     V( I, N-K+I ) = ONE
!
!                    T(i+1:k,i) :=
!                            - tau(i) * V(i+1:k,1:n-k+i) * V(i,1:n-k+i)'
!
                     CALL DGEMV( 'No transpose', K-I, N-K+I, -TAU( I ), &
                                 V( I+1, 1 ), LDV, V( I, 1 ), LDV, ZERO, &
                                 T( I+1, I ), 1 )
                     V( I, N-K+I ) = VII
                  END IF
!
!                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
!
                  CALL DTRMV( 'Lower', 'No transpose', 'Non-unit', K-I, &
                              T( I+1, I+1 ), LDT, T( I+1, I ), 1 )
               END IF
               T( I, I ) = TAU( I )
            END IF
   40    CONTINUE
      END IF
      RETURN
!
!     End of DLARFT
!
      END SUBROUTINE DLARFT


      SUBROUTINE DLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      DOUBLE PRECISION   TAU
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DLARF applies a real elementary reflector H to a real m by n matrix
!  C, from either the left or the right. H is represented in the form
!
!        H = I - tau * v * v'
!
!  where tau is a real scalar and v is a real vector.
!
!  If tau = 0, then H is taken to be the unit matrix.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': form  H * C
!          = 'R': form  C * H
!
!  M       (input) INTEGER
!          The number of rows of the matrix C.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C.
!
!  V       (input) DOUBLE PRECISION array, dimension
!                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
!                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
!          The vector v in the representation of H. V is not used if
!          TAU = 0.
!
!  INCV    (input) INTEGER
!          The increment between elements of v. INCV <> 0.
!
!  TAU     (input) DOUBLE PRECISION
!          The value tau in the representation of H.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
!          or C * H if SIDE = 'R'.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension
!                         (N) if SIDE = 'L'
!                      or (M) if SIDE = 'R'
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DGEMV, DGER
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           LSAME
!     ..
!     .. Executable Statements ..
!
      IF( LSAME( SIDE, 'L' ) ) THEN
!
!        Form  H * C
!
         IF( TAU.NE.ZERO ) THEN
!
!           w := C' * v
!
            CALL DGEMV( 'Transpose', M, N, ONE, C, LDC, V, INCV, ZERO, &
                        WORK, 1 )
!
!           C := C - v * w'
!
            CALL DGER( M, N, -TAU, V, INCV, WORK, 1, C, LDC )
         END IF
      ELSE
!
!        Form  C * H
!
         IF( TAU.NE.ZERO ) THEN
!
!           w := C * v
!
            CALL DGEMV( 'No transpose', M, N, ONE, C, LDC, V, INCV, &
                        ZERO, WORK, 1 )
!
!           C := C - w * v'
!
            CALL DGER( M, N, -TAU, WORK, 1, V, INCV, C, LDC )
         END IF
      END IF
      RETURN
!
!     End of DLARF
!
      END SUBROUTINE DLARF


      SUBROUTINE DLAEXC( WANTQ, N, T, LDT, Q, LDQ, J1, N1, N2, WORK, &
                         INFO )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      LOGICAL            WANTQ
      INTEGER            INFO, J1, LDQ, LDT, N, N1, N2
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   Q( LDQ, * ), T( LDT, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DLAEXC swaps adjacent diagonal blocks T11 and T22 of order 1 or 2 in
!  an upper quasi-triangular matrix T by an orthogonal similarity
!  transformation.
!
!  T must be in Schur canonical form, that is, block upper triangular
!  with 1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block
!  has its diagonal elemnts equal and its off-diagonal elements of
!  opposite sign.
!
!  Arguments
!  =========
!
!  WANTQ   (input) LOGICAL
!          = .TRUE. : accumulate the transformation in the matrix Q;
!          = .FALSE.: do not accumulate the transformation.
!
!  N       (input) INTEGER
!          The order of the matrix T. N >= 0.
!
!  T       (input/output) DOUBLE PRECISION array, dimension (LDT,N)
!          On entry, the upper quasi-triangular matrix T, in Schur
!          canonical form.
!          On exit, the updated matrix T, again in Schur canonical form.
!
!  LDT     (input)  INTEGER
!          The leading dimension of the array T. LDT >= max(1,N).
!
!  Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
!          On entry, if WANTQ is .TRUE., the orthogonal matrix Q.
!          On exit, if WANTQ is .TRUE., the updated matrix Q.
!          If WANTQ is .FALSE., Q is not referenced.
!
!  LDQ     (input) INTEGER
!          The leading dimension of the array Q.
!          LDQ >= 1; and if WANTQ is .TRUE., LDQ >= N.
!
!  J1      (input) INTEGER
!          The index of the first row of the first block T11.
!
!  N1      (input) INTEGER
!          The order of the first block T11. N1 = 0, 1 or 2.
!
!  N2      (input) INTEGER
!          The order of the second block T22. N2 = 0, 1 or 2.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          = 1: the transformed matrix T would be too far from Schur
!               form; the blocks are not swapped and T and Q are
!               unchanged.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   TEN
      PARAMETER          ( TEN = 1.0D+1 )
      INTEGER            LDD, LDX
      PARAMETER          ( LDD = 4, LDX = 2 )
!     ..
!     .. Local Scalars ..
      INTEGER            IERR, J2, J3, J4, K, ND
      DOUBLE PRECISION   CS, DNORM, EPS, SCALE, SMLNUM, SN, T11, T22, &
                         T33, TAU, TAU1, TAU2, TEMP, THRESH, WI1, WI2, &
                         WR1, WR2, XNORM
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   D( LDD, 4 ), U( 3 ), U1( 3 ), U2( 3 ), &
                         X( LDX, 2 )
!     ..
!     .. External Functions ..
!      DOUBLE PRECISION   DLAMCH, DLANGE
!      EXTERNAL           DLAMCH, DLANGE
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DLACPY, DLANV2, DLARFG, DLARFX, DLARTG, DLASY2, &
!                         DROT
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
!     ..
!     .. Executable Statements ..
!
      INFO = 0
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. N1.EQ.0 .OR. N2.EQ.0 ) &
         RETURN
      IF( J1+N1.GT.N ) &
         RETURN
!
      J2 = J1 + 1
      J3 = J1 + 2
      J4 = J1 + 3
!
      IF( N1.EQ.1 .AND. N2.EQ.1 ) THEN
!
!        Swap two 1-by-1 blocks.
!
         T11 = T( J1, J1 )
         T22 = T( J2, J2 )
!
!        Determine the transformation to perform the interchange.
!
         CALL DLARTG( T( J1, J2 ), T22-T11, CS, SN, TEMP )
!
!        Apply transformation to the matrix T.
!
         IF( J3.LE.N ) &
            CALL DROT( N-J1-1, T( J1, J3 ), LDT, T( J2, J3 ), LDT, CS, &
                       SN )
         CALL DROT( J1-1, T( 1, J1 ), 1, T( 1, J2 ), 1, CS, SN )
!
         T( J1, J1 ) = T22
         T( J2, J2 ) = T11
!
         IF( WANTQ ) THEN
!
!           Accumulate transformation in the matrix Q.
!
            CALL DROT( N, Q( 1, J1 ), 1, Q( 1, J2 ), 1, CS, SN )
         END IF
!
      ELSE
!
!        Swapping involves at least one 2-by-2 block.
!
!        Copy the diagonal block of order N1+N2 to the local array D
!        and compute its norm.
!
         ND = N1 + N2
         CALL DLACPY( 'Full', ND, ND, T( J1, J1 ), LDT, D, LDD )
         DNORM = DLANGE( 'Max', ND, ND, D, LDD, WORK )
!
!        Compute machine-dependent threshold for test for accepting
!        swap.
!
         EPS = DLAMCH( 'P' )
         SMLNUM = DLAMCH( 'S' ) / EPS
         THRESH = MAX( TEN*EPS*DNORM, SMLNUM )
!
!        Solve T11*X - X*T22 = scale*T12 for X.
!
         CALL DLASY2( .FALSE., .FALSE., -1, N1, N2, D, LDD, &
                      D( N1+1, N1+1 ), LDD, D( 1, N1+1 ), LDD, SCALE, X, &
                      LDX, XNORM, IERR )
!
!        Swap the adjacent diagonal blocks.
!
         K = N1 + N1 + N2 - 3
         GO TO ( 10, 20, 30 )K
!
   10    CONTINUE
!
!        N1 = 1, N2 = 2: generate elementary reflector H so that:
!
!        ( scale, X11, X12 ) H = ( 0, 0, * )
!
         U( 1 ) = SCALE
         U( 2 ) = X( 1, 1 )
         U( 3 ) = X( 1, 2 )
         CALL DLARFG( 3, U( 3 ), U, 1, TAU )
         U( 3 ) = ONE
         T11 = T( J1, J1 )
!
!        Perform swap provisionally on diagonal block in D.
!
         CALL DLARFX( 'L', 3, 3, U, TAU, D, LDD, WORK )
         CALL DLARFX( 'R', 3, 3, U, TAU, D, LDD, WORK )
!
!        Test whether to reject swap.
!
         IF( MAX( ABS( D( 3, 1 ) ), ABS( D( 3, 2 ) ), ABS( D( 3, &
             3 )-T11 ) ).GT.THRESH )GO TO 50
!
!        Accept swap: apply transformation to the entire matrix T.
!
         CALL DLARFX( 'L', 3, N-J1+1, U, TAU, T( J1, J1 ), LDT, WORK )
         CALL DLARFX( 'R', J2, 3, U, TAU, T( 1, J1 ), LDT, WORK )
!
         T( J3, J1 ) = ZERO
         T( J3, J2 ) = ZERO
         T( J3, J3 ) = T11
!
         IF( WANTQ ) THEN
!
!           Accumulate transformation in the matrix Q.
!
            CALL DLARFX( 'R', N, 3, U, TAU, Q( 1, J1 ), LDQ, WORK )
         END IF
         GO TO 40
!
   20    CONTINUE
!
!        N1 = 2, N2 = 1: generate elementary reflector H so that:
!
!        H (  -X11 ) = ( * )
!          (  -X21 ) = ( 0 )
!          ( scale ) = ( 0 )
!
         U( 1 ) = -X( 1, 1 )
         U( 2 ) = -X( 2, 1 )
         U( 3 ) = SCALE
         CALL DLARFG( 3, U( 1 ), U( 2 ), 1, TAU )
         U( 1 ) = ONE
         T33 = T( J3, J3 )
!
!        Perform swap provisionally on diagonal block in D.
!
         CALL DLARFX( 'L', 3, 3, U, TAU, D, LDD, WORK )
         CALL DLARFX( 'R', 3, 3, U, TAU, D, LDD, WORK )
!
!        Test whether to reject swap.
!
         IF( MAX( ABS( D( 2, 1 ) ), ABS( D( 3, 1 ) ), ABS( D( 1, &
             1 )-T33 ) ).GT.THRESH )GO TO 50
!
!        Accept swap: apply transformation to the entire matrix T.
!
         CALL DLARFX( 'R', J3, 3, U, TAU, T( 1, J1 ), LDT, WORK )
         CALL DLARFX( 'L', 3, N-J1, U, TAU, T( J1, J2 ), LDT, WORK )
!
         T( J1, J1 ) = T33
         T( J2, J1 ) = ZERO
         T( J3, J1 ) = ZERO
!
         IF( WANTQ ) THEN
!
!           Accumulate transformation in the matrix Q.
!
            CALL DLARFX( 'R', N, 3, U, TAU, Q( 1, J1 ), LDQ, WORK )
         END IF
         GO TO 40
!
   30    CONTINUE
!
!        N1 = 2, N2 = 2: generate elementary reflectors H(1) and H(2) so
!        that:
!
!        H(2) H(1) (  -X11  -X12 ) = (  *  * )
!                  (  -X21  -X22 )   (  0  * )
!                  ( scale    0  )   (  0  0 )
!                  (    0  scale )   (  0  0 )
!
         U1( 1 ) = -X( 1, 1 )
         U1( 2 ) = -X( 2, 1 )
         U1( 3 ) = SCALE
         CALL DLARFG( 3, U1( 1 ), U1( 2 ), 1, TAU1 )
         U1( 1 ) = ONE
!
         TEMP = -TAU1*( X( 1, 2 )+U1( 2 )*X( 2, 2 ) )
         U2( 1 ) = -TEMP*U1( 2 ) - X( 2, 2 )
         U2( 2 ) = -TEMP*U1( 3 )
         U2( 3 ) = SCALE
         CALL DLARFG( 3, U2( 1 ), U2( 2 ), 1, TAU2 )
         U2( 1 ) = ONE
!
!        Perform swap provisionally on diagonal block in D.
!
         CALL DLARFX( 'L', 3, 4, U1, TAU1, D, LDD, WORK )
         CALL DLARFX( 'R', 4, 3, U1, TAU1, D, LDD, WORK )
         CALL DLARFX( 'L', 3, 4, U2, TAU2, D( 2, 1 ), LDD, WORK )
         CALL DLARFX( 'R', 4, 3, U2, TAU2, D( 1, 2 ), LDD, WORK )
!
!        Test whether to reject swap.
!
         IF( MAX( ABS( D( 3, 1 ) ), ABS( D( 3, 2 ) ), ABS( D( 4, 1 ) ), &
             ABS( D( 4, 2 ) ) ).GT.THRESH )GO TO 50
!
!        Accept swap: apply transformation to the entire matrix T.
!
         CALL DLARFX( 'L', 3, N-J1+1, U1, TAU1, T( J1, J1 ), LDT, WORK )
         CALL DLARFX( 'R', J4, 3, U1, TAU1, T( 1, J1 ), LDT, WORK )
         CALL DLARFX( 'L', 3, N-J1+1, U2, TAU2, T( J2, J1 ), LDT, WORK )
         CALL DLARFX( 'R', J4, 3, U2, TAU2, T( 1, J2 ), LDT, WORK )
!
         T( J3, J1 ) = ZERO
         T( J3, J2 ) = ZERO
         T( J4, J1 ) = ZERO
         T( J4, J2 ) = ZERO
!
         IF( WANTQ ) THEN
!
!           Accumulate transformation in the matrix Q.
!
            CALL DLARFX( 'R', N, 3, U1, TAU1, Q( 1, J1 ), LDQ, WORK )
            CALL DLARFX( 'R', N, 3, U2, TAU2, Q( 1, J2 ), LDQ, WORK )
         END IF
!
   40    CONTINUE
!
         IF( N2.EQ.2 ) THEN
!
!           Standardize new 2-by-2 block T11
!
            CALL DLANV2( T( J1, J1 ), T( J1, J2 ), T( J2, J1 ), &
                         T( J2, J2 ), WR1, WI1, WR2, WI2, CS, SN )
            CALL DROT( N-J1-1, T( J1, J1+2 ), LDT, T( J2, J1+2 ), LDT, &
                       CS, SN )
            CALL DROT( J1-1, T( 1, J1 ), 1, T( 1, J2 ), 1, CS, SN )
            IF( WANTQ ) &
               CALL DROT( N, Q( 1, J1 ), 1, Q( 1, J2 ), 1, CS, SN )
         END IF
!
         IF( N1.EQ.2 ) THEN
!
!           Standardize new 2-by-2 block T22
!
            J3 = J1 + N2
            J4 = J3 + 1
            CALL DLANV2( T( J3, J3 ), T( J3, J4 ), T( J4, J3 ), &
                         T( J4, J4 ), WR1, WI1, WR2, WI2, CS, SN )
            IF( J3+2.LE.N ) &
               CALL DROT( N-J3-1, T( J3, J3+2 ), LDT, T( J4, J3+2 ), &
                          LDT, CS, SN )
            CALL DROT( J3-1, T( 1, J3 ), 1, T( 1, J4 ), 1, CS, SN )
            IF( WANTQ ) &
               CALL DROT( N, Q( 1, J3 ), 1, Q( 1, J4 ), 1, CS, SN )
         END IF
!
      END IF
      RETURN
!
!     Exit with INFO = 1 if swap was rejected.
!
   50 CONTINUE
      INFO = 1
      RETURN
!
!     End of DLAEXC
!
      END SUBROUTINE DLAEXC


      SUBROUTINE DLANV2( A, B, C, D, RT1R, RT1I, RT2R, RT2I, CS, SN )
!
!  -- LAPACK driver routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, CS, D, RT1I, RT1R, RT2I, RT2R, SN
!     ..
!
!  Purpose
!  =======
!
!  DLANV2 computes the Schur factorization of a real 2-by-2 nonsymmetric
!  matrix in standard form:
!
!       [ A  B ] = [ CS -SN ] [ AA  BB ] [ CS  SN ]
!       [ C  D ]   [ SN  CS ] [ CC  DD ] [-SN  CS ]
!
!  where either
!  1) CC = 0 so that AA and DD are real eigenvalues of the matrix, or
!  2) AA = DD and BB*CC < 0, so that AA + or - sqrt(BB*CC) are complex
!  conjugate eigenvalues.
!
!  Arguments
!  =========
!
!  A       (input/output) DOUBLE PRECISION
!  B       (input/output) DOUBLE PRECISION
!  C       (input/output) DOUBLE PRECISION
!  D       (input/output) DOUBLE PRECISION
!          On entry, the elements of the input matrix.
!          On exit, they are overwritten by the elements of the
!          standardised Schur form.
!
!  RT1R    (output) DOUBLE PRECISION
!  RT1I    (output) DOUBLE PRECISION
!  RT2R    (output) DOUBLE PRECISION
!  RT2I    (output) DOUBLE PRECISION
!          The real and imaginary parts of the eigenvalues. If the
!          eigenvalues are a complex conjugate pair, RT1I > 0.
!
!  CS      (output) DOUBLE PRECISION
!  SN      (output) DOUBLE PRECISION
!          Parameters of the rotation matrix.
!
!  Further Details
!  ===============
!
!  Modified by V. Sima, Research Institute for Informatics, Bucharest,
!  Romania, to reduce the risk of cancellation errors,
!  when computing real eigenvalues, and to ensure, if possible, that
!  abs(RT1R) >= abs(RT2R).
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, HALF, ONE
      PARAMETER          ( ZERO = 0.0D+0, HALF = 0.5D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   MULTPL
      PARAMETER          ( MULTPL = 4.0D+0 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   AA, BB, BCMAX, BCMIS, CC, CS1, DD, EPS, P, SAB, &
                         SAC, SCALE, SIGMA, SN1, TAU, TEMP, Z
!     ..
!     .. External Functions ..
!      DOUBLE PRECISION   DLAMCH, DLAPY2
!      EXTERNAL           DLAMCH, DLAPY2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
      EPS = DLAMCH( 'P' )
      IF( C.EQ.ZERO ) THEN
         CS = ONE
         SN = ZERO
         GO TO 10
!
      ELSE IF( B.EQ.ZERO ) THEN
!
!        Swap rows and columns
!
         CS = ZERO
         SN = ONE
         TEMP = D
         D = A
         A = TEMP
         B = -C
         C = ZERO
         GO TO 10
      ELSE IF( ( A-D ).EQ.ZERO .AND. SIGN( ONE, B ).NE.SIGN( ONE, C ) ) &
                THEN
         CS = ONE
         SN = ZERO
         GO TO 10
      ELSE
!
         TEMP = A - D
         P = HALF*TEMP
         BCMAX = MAX( ABS( B ), ABS( C ) )
         BCMIS = MIN( ABS( B ), ABS( C ) )*SIGN( ONE, B )*SIGN( ONE, C )
         SCALE = MAX( ABS( P ), BCMAX )
         Z = ( P / SCALE )*P + ( BCMAX / SCALE )*BCMIS
!
!        If Z is of the order of the machine accuracy, postpone the
!        decision on the nature of eigenvalues
!
         IF( Z.GE.MULTPL*EPS ) THEN
!
!           Real eigenvalues. Compute A and D.
!
            Z = P + SIGN( SQRT( SCALE )*SQRT( Z ), P )
            A = D + Z
            D = D - ( BCMAX / Z )*BCMIS
!
!           Compute B and the rotation matrix
!
            TAU = DLAPY2( C, Z )
            CS = Z / TAU
            SN = C / TAU
            B = B - C
            C = ZERO
         ELSE
!
!           Complex eigenvalues, or real (almost) equal eigenvalues.
!           Make diagonal elements equal.
!
            SIGMA = B + C
            TAU = DLAPY2( SIGMA, TEMP )
            CS = SQRT( HALF*( ONE+ABS( SIGMA ) / TAU ) )
            SN = -( P / ( TAU*CS ) )*SIGN( ONE, SIGMA )
!
!           Compute [ AA  BB ] = [ A  B ] [ CS -SN ]
!                   [ CC  DD ]   [ C  D ] [ SN  CS ]
!
            AA = A*CS + B*SN
            BB = -A*SN + B*CS
            CC = C*CS + D*SN
            DD = -C*SN + D*CS
!
!           Compute [ A  B ] = [ CS  SN ] [ AA  BB ]
!                   [ C  D ]   [-SN  CS ] [ CC  DD ]
!
            A = AA*CS + CC*SN
            B = BB*CS + DD*SN
            C = -AA*SN + CC*CS
            D = -BB*SN + DD*CS
!
            TEMP = HALF*( A+D )
            A = TEMP
            D = TEMP
!
            IF( C.NE.ZERO ) THEN
               IF( B.NE.ZERO ) THEN
                  IF( SIGN( ONE, B ).EQ.SIGN( ONE, C ) ) THEN
!
!                    Real eigenvalues: reduce to upper triangular form
!
                     SAB = SQRT( ABS( B ) )
                     SAC = SQRT( ABS( C ) )
                     P = SIGN( SAB*SAC, C )
                     TAU = ONE / SQRT( ABS( B+C ) )
                     A = TEMP + P
                     D = TEMP - P
                     B = B - C
                     C = ZERO
                     CS1 = SAB*TAU
                     SN1 = SAC*TAU
                     TEMP = CS*CS1 - SN*SN1
                     SN = CS*SN1 + SN*CS1
                     CS = TEMP
                  END IF
               ELSE
                  B = -C
                  C = ZERO
                  TEMP = CS
                  CS = -SN
                  SN = TEMP
               END IF
            END IF
         END IF
!
      END IF
!
   10 CONTINUE
!
!     Store eigenvalues in (RT1R,RT1I) and (RT2R,RT2I).
!
      RT1R = A
      RT2R = D
      IF( C.EQ.ZERO ) THEN
         RT1I = ZERO
         RT2I = ZERO
      ELSE
         RT1I = SQRT( ABS( B ) )*SQRT( ABS( C ) )
         RT2I = -RT1I
      END IF
      RETURN
!
!     End of DLANV2
!
      END SUBROUTINE DLANV2


      SUBROUTINE DLASY2( LTRANL, LTRANR, ISGN, N1, N2, TL, LDTL, TR, &
                         LDTR, B, LDB, SCALE, X, LDX, XNORM, INFO )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      LOGICAL            LTRANL, LTRANR
      INTEGER            INFO, ISGN, LDB, LDTL, LDTR, LDX, N1, N2
      DOUBLE PRECISION   SCALE, XNORM
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   B( LDB, * ), TL( LDTL, * ), TR( LDTR, * ), &
                         X( LDX, * )
!     ..
!
!  Purpose
!  =======
!
!  DLASY2 solves for the N1 by N2 matrix X, 1 <= N1,N2 <= 2, in
!
!         op(TL)*X + ISGN*X*op(TR) = SCALE*B,
!
!  where TL is N1 by N1, TR is N2 by N2, B is N1 by N2, and ISGN = 1 or
!  -1.  op(T) = T or T', where T' denotes the transpose of T.
!
!  Arguments
!  =========
!
!  LTRANL  (input) LOGICAL
!          On entry, LTRANL specifies the op(TL):
!             = .FALSE., op(TL) = TL,
!             = .TRUE., op(TL) = TL'.
!
!  LTRANR  (input) LOGICAL
!          On entry, LTRANR specifies the op(TR):
!            = .FALSE., op(TR) = TR,
!            = .TRUE., op(TR) = TR'.
!
!  ISGN    (input) INTEGER
!          On entry, ISGN specifies the sign of the equation
!          as described before. ISGN may only be 1 or -1.
!
!  N1      (input) INTEGER
!          On entry, N1 specifies the order of matrix TL.
!          N1 may only be 0, 1 or 2.
!
!  N2      (input) INTEGER
!          On entry, N2 specifies the order of matrix TR.
!          N2 may only be 0, 1 or 2.
!
!  TL      (input) DOUBLE PRECISION array, dimension (LDTL,2)
!          On entry, TL contains an N1 by N1 matrix.
!
!  LDTL    (input) INTEGER
!          The leading dimension of the matrix TL. LDTL >= max(1,N1).
!
!  TR      (input) DOUBLE PRECISION array, dimension (LDTR,2)
!          On entry, TR contains an N2 by N2 matrix.
!
!  LDTR    (input) INTEGER
!          The leading dimension of the matrix TR. LDTR >= max(1,N2).
!
!  B       (input) DOUBLE PRECISION array, dimension (LDB,2)
!          On entry, the N1 by N2 matrix B contains the right-hand
!          side of the equation.
!
!  LDB     (input) INTEGER
!          The leading dimension of the matrix B. LDB >= max(1,N1).
!
!  SCALE   (output) DOUBLE PRECISION
!          On exit, SCALE contains the scale factor. SCALE is chosen
!          less than or equal to 1 to prevent the solution overflowing.
!
!  X       (output) DOUBLE PRECISION array, dimension (LDX,2)
!          On exit, X contains the N1 by N2 solution.
!
!  LDX     (input) INTEGER
!          The leading dimension of the matrix X. LDX >= max(1,N1).
!
!  XNORM   (output) DOUBLE PRECISION
!          On exit, XNORM is the infinity-norm of the solution.
!
!  INFO    (output) INTEGER
!          On exit, INFO is set to
!             0: successful exit.
!             1: TL and TR have too close eigenvalues, so TL or
!                TR is perturbed to get a nonsingular equation.
!          NOTE: In the interests of speed, this routine does not
!                check the inputs for errors.
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   TWO, HALF, EIGHT
      PARAMETER          ( TWO = 2.0D+0, HALF = 0.5D+0, EIGHT = 8.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            BSWAP, XSWAP
      INTEGER            I, IP, IPIV, IPSV, J, JP, JPSV, K
      DOUBLE PRECISION   BET, EPS, GAM, L21, SGN, SMIN, SMLNUM, TAU1, &
                         TEMP, U11, U12, U22, XMAX
!     ..
!     .. Local Arrays ..
      LOGICAL            BSWPIV( 4 ), XSWPIV( 4 )
      INTEGER            JPIV( 4 ), LOCL21( 4 ), LOCU12( 4 ), &
                         LOCU22( 4 )
      DOUBLE PRECISION   BTMP( 4 ), T16( 4, 4 ), TMP( 4 ), X2( 2 )
!     ..
!     .. External Functions ..
!      INTEGER            IDAMAX
!      DOUBLE PRECISION   DLAMCH
!      EXTERNAL           IDAMAX, DLAMCH
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DCOPY, DSWAP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
!     ..
!     .. Data statements ..
      DATA               LOCU12 / 3, 4, 1, 2 / , LOCL21 / 2, 1, 4, 3 / , &
                         LOCU22 / 4, 3, 2, 1 /
      DATA               XSWPIV / .FALSE., .FALSE., .TRUE., .TRUE. /
      DATA               BSWPIV / .FALSE., .TRUE., .FALSE., .TRUE. /
!     ..
!     .. Executable Statements ..
!
!     Do not check the input parameters for errors
!
      INFO = 0
!
!     Quick return if possible
!
      IF( N1.EQ.0 .OR. N2.EQ.0 ) &
         RETURN
!
!     Set constants to control overflow
!
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' ) / EPS
      SGN = ISGN
!
      K = N1 + N1 + N2 - 2
      GO TO ( 10, 20, 30, 50 )K
!
!     1 by 1: TL11*X + SGN*X*TR11 = B11
!
   10 CONTINUE
      TAU1 = TL( 1, 1 ) + SGN*TR( 1, 1 )
      BET = ABS( TAU1 )
      IF( BET.LE.SMLNUM ) THEN
         TAU1 = SMLNUM
         BET = SMLNUM
         INFO = 1
      END IF
!
      SCALE = ONE
      GAM = ABS( B( 1, 1 ) )
      IF( SMLNUM*GAM.GT.BET ) &
         SCALE = ONE / GAM
!
      X( 1, 1 ) = ( B( 1, 1 )*SCALE ) / TAU1
      XNORM = ABS( X( 1, 1 ) )
      RETURN
!
!     1 by 2:
!     TL11*[X11 X12] + ISGN*[X11 X12]*op[TR11 TR12]  = [B11 B12]
!                                       [TR21 TR22]
!
   20 CONTINUE
!
      SMIN = MAX( EPS*MAX( ABS( TL( 1, 1 ) ), ABS( TR( 1, 1 ) ), &
             ABS( TR( 1, 2 ) ), ABS( TR( 2, 1 ) ), ABS( TR( 2, 2 ) ) ), &
             SMLNUM )
      TMP( 1 ) = TL( 1, 1 ) + SGN*TR( 1, 1 )
      TMP( 4 ) = TL( 1, 1 ) + SGN*TR( 2, 2 )
      IF( LTRANR ) THEN
         TMP( 2 ) = SGN*TR( 2, 1 )
         TMP( 3 ) = SGN*TR( 1, 2 )
      ELSE
         TMP( 2 ) = SGN*TR( 1, 2 )
         TMP( 3 ) = SGN*TR( 2, 1 )
      END IF
      BTMP( 1 ) = B( 1, 1 )
      BTMP( 2 ) = B( 1, 2 )
      GO TO 40
!
!     2 by 1:
!          op[TL11 TL12]*[X11] + ISGN* [X11]*TR11  = [B11]
!            [TL21 TL22] [X21]         [X21]         [B21]
!
   30 CONTINUE
      SMIN = MAX( EPS*MAX( ABS( TR( 1, 1 ) ), ABS( TL( 1, 1 ) ), &
             ABS( TL( 1, 2 ) ), ABS( TL( 2, 1 ) ), ABS( TL( 2, 2 ) ) ), &
             SMLNUM )
      TMP( 1 ) = TL( 1, 1 ) + SGN*TR( 1, 1 )
      TMP( 4 ) = TL( 2, 2 ) + SGN*TR( 1, 1 )
      IF( LTRANL ) THEN
         TMP( 2 ) = TL( 1, 2 )
         TMP( 3 ) = TL( 2, 1 )
      ELSE
         TMP( 2 ) = TL( 2, 1 )
         TMP( 3 ) = TL( 1, 2 )
      END IF
      BTMP( 1 ) = B( 1, 1 )
      BTMP( 2 ) = B( 2, 1 )
   40 CONTINUE
!
!     Solve 2 by 2 system using complete pivoting.
!     Set pivots less than SMIN to SMIN.
!
      IPIV = IDAMAX( 4, TMP, 1 )
      U11 = TMP( IPIV )
      IF( ABS( U11 ).LE.SMIN ) THEN
         INFO = 1
         U11 = SMIN
      END IF
      U12 = TMP( LOCU12( IPIV ) )
      L21 = TMP( LOCL21( IPIV ) ) / U11
      U22 = TMP( LOCU22( IPIV ) ) - U12*L21
      XSWAP = XSWPIV( IPIV )
      BSWAP = BSWPIV( IPIV )
      IF( ABS( U22 ).LE.SMIN ) THEN
         INFO = 1
         U22 = SMIN
      END IF
      IF( BSWAP ) THEN
         TEMP = BTMP( 2 )
         BTMP( 2 ) = BTMP( 1 ) - L21*TEMP
         BTMP( 1 ) = TEMP
      ELSE
         BTMP( 2 ) = BTMP( 2 ) - L21*BTMP( 1 )
      END IF
      SCALE = ONE
      IF( ( TWO*SMLNUM )*ABS( BTMP( 2 ) ).GT.ABS( U22 ) .OR. &
          ( TWO*SMLNUM )*ABS( BTMP( 1 ) ).GT.ABS( U11 ) ) THEN
         SCALE = HALF / MAX( ABS( BTMP( 1 ) ), ABS( BTMP( 2 ) ) )
         BTMP( 1 ) = BTMP( 1 )*SCALE
         BTMP( 2 ) = BTMP( 2 )*SCALE
      END IF
      X2( 2 ) = BTMP( 2 ) / U22
      X2( 1 ) = BTMP( 1 ) / U11 - ( U12 / U11 )*X2( 2 )
      IF( XSWAP ) THEN
         TEMP = X2( 2 )
         X2( 2 ) = X2( 1 )
         X2( 1 ) = TEMP
      END IF
      X( 1, 1 ) = X2( 1 )
      IF( N1.EQ.1 ) THEN
         X( 1, 2 ) = X2( 2 )
         XNORM = ABS( X( 1, 1 ) ) + ABS( X( 1, 2 ) )
      ELSE
         X( 2, 1 ) = X2( 2 )
         XNORM = MAX( ABS( X( 1, 1 ) ), ABS( X( 2, 1 ) ) )
      END IF
      RETURN
!
!     2 by 2:
!     op[TL11 TL12]*[X11 X12] +ISGN* [X11 X12]*op[TR11 TR12] = [B11 B12]
!       [TL21 TL22] [X21 X22]        [X21 X22]   [TR21 TR22]   [B21 B22]
!
!     Solve equivalent 4 by 4 system using complete pivoting.
!     Set pivots less than SMIN to SMIN.
!
   50 CONTINUE
      SMIN = MAX( ABS( TR( 1, 1 ) ), ABS( TR( 1, 2 ) ), &
             ABS( TR( 2, 1 ) ), ABS( TR( 2, 2 ) ) )
      SMIN = MAX( SMIN, ABS( TL( 1, 1 ) ), ABS( TL( 1, 2 ) ), &
             ABS( TL( 2, 1 ) ), ABS( TL( 2, 2 ) ) )
      SMIN = MAX( EPS*SMIN, SMLNUM )
      BTMP( 1 ) = ZERO
      CALL DCOPY( 16, BTMP, 0, T16, 1 )
      T16( 1, 1 ) = TL( 1, 1 ) + SGN*TR( 1, 1 )
      T16( 2, 2 ) = TL( 2, 2 ) + SGN*TR( 1, 1 )
      T16( 3, 3 ) = TL( 1, 1 ) + SGN*TR( 2, 2 )
      T16( 4, 4 ) = TL( 2, 2 ) + SGN*TR( 2, 2 )
      IF( LTRANL ) THEN
         T16( 1, 2 ) = TL( 2, 1 )
         T16( 2, 1 ) = TL( 1, 2 )
         T16( 3, 4 ) = TL( 2, 1 )
         T16( 4, 3 ) = TL( 1, 2 )
      ELSE
         T16( 1, 2 ) = TL( 1, 2 )
         T16( 2, 1 ) = TL( 2, 1 )
         T16( 3, 4 ) = TL( 1, 2 )
         T16( 4, 3 ) = TL( 2, 1 )
      END IF
      IF( LTRANR ) THEN
         T16( 1, 3 ) = SGN*TR( 1, 2 )
         T16( 2, 4 ) = SGN*TR( 1, 2 )
         T16( 3, 1 ) = SGN*TR( 2, 1 )
         T16( 4, 2 ) = SGN*TR( 2, 1 )
      ELSE
         T16( 1, 3 ) = SGN*TR( 2, 1 )
         T16( 2, 4 ) = SGN*TR( 2, 1 )
         T16( 3, 1 ) = SGN*TR( 1, 2 )
         T16( 4, 2 ) = SGN*TR( 1, 2 )
      END IF
      BTMP( 1 ) = B( 1, 1 )
      BTMP( 2 ) = B( 2, 1 )
      BTMP( 3 ) = B( 1, 2 )
      BTMP( 4 ) = B( 2, 2 )
!
!     Perform elimination
!
      DO 100 I = 1, 3
         XMAX = ZERO
         DO 70 IP = I, 4
            DO 60 JP = I, 4
               IF( ABS( T16( IP, JP ) ).GE.XMAX ) THEN
                  XMAX = ABS( T16( IP, JP ) )
                  IPSV = IP
                  JPSV = JP
               END IF
   60       CONTINUE
   70    CONTINUE
         IF( IPSV.NE.I ) THEN
            CALL DSWAP( 4, T16( IPSV, 1 ), 4, T16( I, 1 ), 4 )
            TEMP = BTMP( I )
            BTMP( I ) = BTMP( IPSV )
            BTMP( IPSV ) = TEMP
         END IF
         IF( JPSV.NE.I ) &
            CALL DSWAP( 4, T16( 1, JPSV ), 1, T16( 1, I ), 1 )
         JPIV( I ) = JPSV
         IF( ABS( T16( I, I ) ).LT.SMIN ) THEN
            INFO = 1
            T16( I, I ) = SMIN
         END IF
         DO 90 J = I + 1, 4
            T16( J, I ) = T16( J, I ) / T16( I, I )
            BTMP( J ) = BTMP( J ) - T16( J, I )*BTMP( I )
            DO 80 K = I + 1, 4
               T16( J, K ) = T16( J, K ) - T16( J, I )*T16( I, K )
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE
      IF( ABS( T16( 4, 4 ) ).LT.SMIN ) &
         T16( 4, 4 ) = SMIN
      SCALE = ONE
      IF( ( EIGHT*SMLNUM )*ABS( BTMP( 1 ) ).GT.ABS( T16( 1, 1 ) ) .OR. &
          ( EIGHT*SMLNUM )*ABS( BTMP( 2 ) ).GT.ABS( T16( 2, 2 ) ) .OR. &
          ( EIGHT*SMLNUM )*ABS( BTMP( 3 ) ).GT.ABS( T16( 3, 3 ) ) .OR. &
          ( EIGHT*SMLNUM )*ABS( BTMP( 4 ) ).GT.ABS( T16( 4, 4 ) ) ) THEN
         SCALE = ( ONE / EIGHT ) / MAX( ABS( BTMP( 1 ) ), &
                 ABS( BTMP( 2 ) ), ABS( BTMP( 3 ) ), ABS( BTMP( 4 ) ) )
         BTMP( 1 ) = BTMP( 1 )*SCALE
         BTMP( 2 ) = BTMP( 2 )*SCALE
         BTMP( 3 ) = BTMP( 3 )*SCALE
         BTMP( 4 ) = BTMP( 4 )*SCALE
      END IF
      DO 120 I = 1, 4
         K = 5 - I
         TEMP = ONE / T16( K, K )
         TMP( K ) = BTMP( K )*TEMP
         DO 110 J = K + 1, 4
            TMP( K ) = TMP( K ) - ( TEMP*T16( K, J ) )*TMP( J )
  110    CONTINUE
  120 CONTINUE
      DO 130 I = 1, 3
         IF( JPIV( 4-I ).NE.4-I ) THEN
            TEMP = TMP( 4-I )
            TMP( 4-I ) = TMP( JPIV( 4-I ) )
            TMP( JPIV( 4-I ) ) = TEMP
         END IF
  130 CONTINUE
      X( 1, 1 ) = TMP( 1 )
      X( 2, 1 ) = TMP( 2 )
      X( 1, 2 ) = TMP( 3 )
      X( 2, 2 ) = TMP( 4 )
      XNORM = MAX( ABS( TMP( 1 ) )+ABS( TMP( 3 ) ), &
              ABS( TMP( 2 ) )+ABS( TMP( 4 ) ) )
      RETURN
!
!     End of DLASY2
!
      END SUBROUTINE DLASY2


      DOUBLE PRECISION FUNCTION DLAPY2( X, Y )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   X, Y
!     ..
!
!  Purpose
!  =======
!
!  DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
!  overflow.
!
!  Arguments
!  =========
!
!  X       (input) DOUBLE PRECISION
!  Y       (input) DOUBLE PRECISION
!          X and Y specify the values x and y.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   W, XABS, YABS, Z
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
      XABS = ABS( X )
      YABS = ABS( Y )
      W = MAX( XABS, YABS )
      Z = MIN( XABS, YABS )
      IF( Z.EQ.ZERO ) THEN
         DLAPY2 = W
      ELSE
         DLAPY2 = W*SQRT( ONE+( Z / W )**2 )
      END IF
      RETURN
!
!     End of DLAPY2
!
      END FUNCTION DLAPY2


      DOUBLE PRECISION FUNCTION DLANHS( NORM, N, A, LDA, WORK )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            LDA, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DLANHS  returns the value of the one norm,  or the Frobenius norm, or
!  the  infinity norm,  or the  element of  largest absolute value  of a
!  Hessenberg matrix A.
!
!  Description
!  ===========
!
!  DLANHS returns the value
!
!     DLANHS = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!              (
!              ( norm1(A),         NORM = '1', 'O' or 'o'
!              (
!              ( normI(A),         NORM = 'I' or 'i'
!              (
!              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!
!  where  norm1  denotes the  one norm of a matrix (maximum column sum),
!  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!  normF  denotes the  Frobenius norm of a matrix (square root of sum of
!  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix no
!
!  Arguments
!  =========
!
!  NORM    (input) CHARACTER*1
!          Specifies the value to be returned in DLANHS as described
!          above.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.  When N = 0, DLANHS is
!          set to zero.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!          The n by n upper Hessenberg matrix A; the part of A below the
!          first sub-diagonal is not referenced.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(N,1).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)),
!          where LWORK >= N when NORM = 'I'; otherwise, WORK is not
!          referenced.
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   SCALE, SUM, VALUE
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DLASSQ
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
      IF( N.EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
!
!        Find max(abs(A(i,j))).
!
         VALUE = ZERO
         DO 20 J = 1, N
            DO 10 I = 1, MIN( N, J+1 )
               VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) THEN
!
!        Find norm1(A).
!
         VALUE = ZERO
         DO 40 J = 1, N
            SUM = ZERO
            DO 30 I = 1, MIN( N, J+1 )
               SUM = SUM + ABS( A( I, J ) )
   30       CONTINUE
            VALUE = MAX( VALUE, SUM )
   40    CONTINUE
      ELSE IF( LSAME( NORM, 'I' ) ) THEN
!
!        Find normI(A).
!
         DO 50 I = 1, N
            WORK( I ) = ZERO
   50    CONTINUE
         DO 70 J = 1, N
            DO 60 I = 1, MIN( N, J+1 )
               WORK( I ) = WORK( I ) + ABS( A( I, J ) )
   60       CONTINUE
   70    CONTINUE
         VALUE = ZERO
         DO 80 I = 1, N
            VALUE = MAX( VALUE, WORK( I ) )
   80    CONTINUE
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
         SCALE = ZERO
         SUM = ONE
         DO 90 J = 1, N
            CALL DLASSQ( MIN( N, J+1 ), A( 1, J ), 1, SCALE, SUM )
   90    CONTINUE
         VALUE = SCALE*SQRT( SUM )
      END IF
!
      DLANHS = VALUE
      RETURN
!
!     End of DLANHS
!
      END FUNCTION DLANHS


      SUBROUTINE DTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
!     .. Scalar Arguments ..
      INTEGER INCX,LDA,N
      CHARACTER DIAG,TRANS,UPLO
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*)
!     ..
!
!  Purpose
!  =======
!
!  DTRMV  performs one of the matrix-vector operations
!
!     x := A*x,   or   x := A'*x,
!
!  where x is an n element vector and  A is an n by n unit, or non-unit,
!  upper or lower triangular matrix.
!
!  Arguments
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   x := A*x.
!
!              TRANS = 'T' or 't'   x := A'*x.
!
!              TRANS = 'C' or 'c'   x := A'*x.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular matrix and the strictly lower triangular part of
!           A is not referenced.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular matrix and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u', the diagonal elements of
!           A are not referenced either, but are assumed to be unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x. On exit, X is overwritten with the
!           tranformed vector x.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,J,JX,KX
      LOGICAL NOUNIT
!     ..
!     .. External Functions ..
!      LOGICAL LSAME
!      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!
!     Test the input parameters.
!
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND. &
               .NOT.LSAME(TRANS,'C')) THEN
          INFO = 2
      ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (LDA.LT.MAX(1,N)) THEN
          INFO = 6
      ELSE IF (INCX.EQ.0) THEN
          INFO = 8
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DTRMV ',INFO)
          RETURN
      END IF
!
!     Quick return if possible.
!
      IF (N.EQ.0) RETURN
!
      NOUNIT = LSAME(DIAG,'N')
!
!     Set up the start point in X if the increment is not unity. This
!     will be  ( N - 1 )*INCX  too small for descending loops.
!
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      IF (LSAME(TRANS,'N')) THEN
!
!        Form  x := A*x.
!
          IF (LSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 20 J = 1,N
                      IF (X(J).NE.ZERO) THEN
                          TEMP = X(J)
                          DO 10 I = 1,J - 1
                              X(I) = X(I) + TEMP*A(I,J)
   10                     CONTINUE
                          IF (NOUNIT) X(J) = X(J)*A(J,J)
                      END IF
   20             CONTINUE
              ELSE
                  JX = KX
                  DO 40 J = 1,N
                      IF (X(JX).NE.ZERO) THEN
                          TEMP = X(JX)
                          IX = KX
                          DO 30 I = 1,J - 1
                              X(IX) = X(IX) + TEMP*A(I,J)
                              IX = IX + INCX
   30                     CONTINUE
                          IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                      END IF
                      JX = JX + INCX
   40             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 60 J = N,1,-1
                      IF (X(J).NE.ZERO) THEN
                          TEMP = X(J)
                          DO 50 I = N,J + 1,-1
                              X(I) = X(I) + TEMP*A(I,J)
   50                     CONTINUE
                          IF (NOUNIT) X(J) = X(J)*A(J,J)
                      END IF
   60             CONTINUE
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 80 J = N,1,-1
                      IF (X(JX).NE.ZERO) THEN
                          TEMP = X(JX)
                          IX = KX
                          DO 70 I = N,J + 1,-1
                              X(IX) = X(IX) + TEMP*A(I,J)
                              IX = IX - INCX
   70                     CONTINUE
                          IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                      END IF
                      JX = JX - INCX
   80             CONTINUE
              END IF
          END IF
      ELSE
!
!        Form  x := A'*x.
!
          IF (LSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 100 J = N,1,-1
                      TEMP = X(J)
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      DO 90 I = J - 1,1,-1
                          TEMP = TEMP + A(I,J)*X(I)
   90                 CONTINUE
                      X(J) = TEMP
  100             CONTINUE
              ELSE
                  JX = KX + (N-1)*INCX
                  DO 120 J = N,1,-1
                      TEMP = X(JX)
                      IX = JX
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      DO 110 I = J - 1,1,-1
                          IX = IX - INCX
                          TEMP = TEMP + A(I,J)*X(IX)
  110                 CONTINUE
                      X(JX) = TEMP
                      JX = JX - INCX
  120             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 140 J = 1,N
                      TEMP = X(J)
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      DO 130 I = J + 1,N
                          TEMP = TEMP + A(I,J)*X(I)
  130                 CONTINUE
                      X(J) = TEMP
  140             CONTINUE
              ELSE
                  JX = KX
                  DO 160 J = 1,N
                      TEMP = X(JX)
                      IX = JX
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      DO 150 I = J + 1,N
                          IX = IX + INCX
                          TEMP = TEMP + A(I,J)*X(IX)
  150                 CONTINUE
                      X(JX) = TEMP
                      JX = JX + INCX
  160             CONTINUE
              END IF
          END IF
      END IF
!
      RETURN
!
!     End of DTRMV .
!
      END SUBROUTINE DTRMV


      SUBROUTINE DTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
!     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*)
!     ..
!
!  Purpose
!  =======
!
!  DTRMM  performs one of the matrix-matrix operations
!
!     B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
!
!  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
!  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!
!     op( A ) = A   or   op( A ) = A'.
!
!  Arguments
!  ==========
!
!  SIDE   - CHARACTER*1.
!           On entry,  SIDE specifies whether  op( A ) multiplies B from
!           the left or right as follows:
!
!              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
!
!              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
!
!           Unchanged on exit.
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix A is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n'   op( A ) = A.
!
!              TRANSA = 'T' or 't'   op( A ) = A'.
!
!              TRANSA = 'C' or 'c'   op( A ) = A'.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit triangular
!           as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of B. M must be at
!           least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of B.  N must be
!           at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
!           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
!           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!           upper triangular part of the array  A must contain the upper
!           triangular matrix  and the strictly lower triangular part of
!           A is not referenced.
!           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!           lower triangular part of the array  A must contain the lower
!           triangular matrix  and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!           A  are not referenced either,  but are assumed to be  unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!           Unchanged on exit.
!
!  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain the matrix  B,  and  on exit  is overwritten  by the
!           transformed matrix.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
!      LOGICAL LSAME
!      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOUNIT,UPPER
!     ..
!     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!
!     Test the input parameters.
!
      LSIDE = LSAME(SIDE,'L')
      IF (LSIDE) THEN
          NROWA = M
      ELSE
          NROWA = N
      END IF
      NOUNIT = LSAME(DIAG,'N')
      UPPER = LSAME(UPLO,'U')
!
      INFO = 0
      IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
          INFO = 1
      ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 2
      ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND. &
               (.NOT.LSAME(TRANSA,'T')) .AND. &
               (.NOT.LSAME(TRANSA,'C'))) THEN
          INFO = 3
      ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) THEN
          INFO = 4
      ELSE IF (M.LT.0) THEN
          INFO = 5
      ELSE IF (N.LT.0) THEN
          INFO = 6
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 9
      ELSE IF (LDB.LT.MAX(1,M)) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DTRMM ',INFO)
          RETURN
      END IF
!
!     Quick return if possible.
!
      IF (N.EQ.0) RETURN
!
!     And when  alpha.eq.zero.
!
      IF (ALPHA.EQ.ZERO) THEN
          DO 20 J = 1,N
              DO 10 I = 1,M
                  B(I,J) = ZERO
   10         CONTINUE
   20     CONTINUE
          RETURN
      END IF
!
!     Start the operations.
!
      IF (LSIDE) THEN
          IF (LSAME(TRANSA,'N')) THEN
!
!           Form  B := alpha*A*B.
!
              IF (UPPER) THEN
                  DO 50 J = 1,N
                      DO 40 K = 1,M
                          IF (B(K,J).NE.ZERO) THEN
                              TEMP = ALPHA*B(K,J)
                              DO 30 I = 1,K - 1
                                  B(I,J) = B(I,J) + TEMP*A(I,K)
   30                         CONTINUE
                              IF (NOUNIT) TEMP = TEMP*A(K,K)
                              B(K,J) = TEMP
                          END IF
   40                 CONTINUE
   50             CONTINUE
              ELSE
                  DO 80 J = 1,N
                      DO 70 K = M,1,-1
                          IF (B(K,J).NE.ZERO) THEN
                              TEMP = ALPHA*B(K,J)
                              B(K,J) = TEMP
                              IF (NOUNIT) B(K,J) = B(K,J)*A(K,K)
                              DO 60 I = K + 1,M
                                  B(I,J) = B(I,J) + TEMP*A(I,K)
   60                         CONTINUE
                          END IF
   70                 CONTINUE
   80             CONTINUE
              END IF
          ELSE
!
!           Form  B := alpha*A'*B.
!
              IF (UPPER) THEN
                  DO 110 J = 1,N
                      DO 100 I = M,1,-1
                          TEMP = B(I,J)
                          IF (NOUNIT) TEMP = TEMP*A(I,I)
                          DO 90 K = 1,I - 1
                              TEMP = TEMP + A(K,I)*B(K,J)
   90                     CONTINUE
                          B(I,J) = ALPHA*TEMP
  100                 CONTINUE
  110             CONTINUE
              ELSE
                  DO 140 J = 1,N
                      DO 130 I = 1,M
                          TEMP = B(I,J)
                          IF (NOUNIT) TEMP = TEMP*A(I,I)
                          DO 120 K = I + 1,M
                              TEMP = TEMP + A(K,I)*B(K,J)
  120                     CONTINUE
                          B(I,J) = ALPHA*TEMP
  130                 CONTINUE
  140             CONTINUE
              END IF
          END IF
      ELSE
          IF (LSAME(TRANSA,'N')) THEN
!
!           Form  B := alpha*B*A.
!
              IF (UPPER) THEN
                  DO 180 J = N,1,-1
                      TEMP = ALPHA
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      DO 150 I = 1,M
                          B(I,J) = TEMP*B(I,J)
  150                 CONTINUE
                      DO 170 K = 1,J - 1
                          IF (A(K,J).NE.ZERO) THEN
                              TEMP = ALPHA*A(K,J)
                              DO 160 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  160                         CONTINUE
                          END IF
  170                 CONTINUE
  180             CONTINUE
              ELSE
                  DO 220 J = 1,N
                      TEMP = ALPHA
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      DO 190 I = 1,M
                          B(I,J) = TEMP*B(I,J)
  190                 CONTINUE
                      DO 210 K = J + 1,N
                          IF (A(K,J).NE.ZERO) THEN
                              TEMP = ALPHA*A(K,J)
                              DO 200 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  200                         CONTINUE
                          END IF
  210                 CONTINUE
  220             CONTINUE
              END IF
          ELSE
!
!           Form  B := alpha*B*A'.
!
              IF (UPPER) THEN
                  DO 260 K = 1,N
                      DO 240 J = 1,K - 1
                          IF (A(J,K).NE.ZERO) THEN
                              TEMP = ALPHA*A(J,K)
                              DO 230 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  230                         CONTINUE
                          END IF
  240                 CONTINUE
                      TEMP = ALPHA
                      IF (NOUNIT) TEMP = TEMP*A(K,K)
                      IF (TEMP.NE.ONE) THEN
                          DO 250 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  250                     CONTINUE
                      END IF
  260             CONTINUE
              ELSE
                  DO 300 K = N,1,-1
                      DO 280 J = K + 1,N
                          IF (A(J,K).NE.ZERO) THEN
                              TEMP = ALPHA*A(J,K)
                              DO 270 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  270                         CONTINUE
                          END IF
  280                 CONTINUE
                      TEMP = ALPHA
                      IF (NOUNIT) TEMP = TEMP*A(K,K)
                      IF (TEMP.NE.ONE) THEN
                          DO 290 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  290                     CONTINUE
                      END IF
  300             CONTINUE
              END IF
          END IF
      END IF
!
      RETURN
!
!     End of DTRMM .
!
      END SUBROUTINE DTRMM


      DOUBLE PRECISION FUNCTION DNRM2(N,X,INCX)
!     .. Scalar Arguments ..
      INTEGER INCX,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION X(*)
!     ..
!
!  Purpose
!  =======
!
!  DNRM2 returns the euclidean norm of a vector via the function
!  name, so that
!
!     DNRM2 := sqrt( x'*x )
!
!
!  -- This version written on 25-October-1982.
!     Modified on 14-October-1993 to inline the call to DLASSQ.
!     Sven Hammarling, Nag Ltd.
!
!
!     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION ABSXI,NORM,SCALE,SSQ
      INTEGER IX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS,SQRT
!     ..
      IF (N.LT.1 .OR. INCX.LT.1) THEN
          NORM = ZERO
      ELSE IF (N.EQ.1) THEN
          NORM = ABS(X(1))
      ELSE
          SCALE = ZERO
          SSQ = ONE
!        The following loop is equivalent to this call to the LAPACK
!        auxiliary routine:
!        CALL DLASSQ( N, X, INCX, SCALE, SSQ )
!
          DO 10 IX = 1,1 + (N-1)*INCX,INCX
              IF (X(IX).NE.ZERO) THEN
                  ABSXI = ABS(X(IX))
                  IF (SCALE.LT.ABSXI) THEN
                      SSQ = ONE + SSQ* (SCALE/ABSXI)**2
                      SCALE = ABSXI
                  ELSE
                      SSQ = SSQ + (ABSXI/SCALE)**2
                  END IF
              END IF
   10     CONTINUE
          NORM = SCALE*SQRT(SSQ)
      END IF
!
      DNRM2 = NORM
      RETURN
!
!     End of DNRM2.
!
      END FUNCTION DNRM2


      SUBROUTINE DLACN2( N, V, X, ISGN, EST, KASE, ISAVE )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            KASE, N
      DOUBLE PRECISION   EST
!     ..
!     .. Array Arguments ..
      INTEGER            ISGN( * ), ISAVE( 3 )
      DOUBLE PRECISION   V( * ), X( * )
!     ..
!
!  Purpose
!  =======
!
!  DLACN2 estimates the 1-norm of a square, real matrix A.
!  Reverse communication is used for evaluating matrix-vector products.
!
!  Arguments
!  =========
!
!  N      (input) INTEGER
!         The order of the matrix.  N >= 1.
!
!  V      (workspace) DOUBLE PRECISION array, dimension (N)
!         On the final return, V = A*W,  where  EST = norm(V)/norm(W)
!         (W is not returned).
!
!  X      (input/output) DOUBLE PRECISION array, dimension (N)
!         On an intermediate return, X should be overwritten by
!               A * X,   if KASE=1,
!               A' * X,  if KASE=2,
!         and DLACN2 must be re-called with all the other parameters
!         unchanged.
!
!  ISGN   (workspace) INTEGER array, dimension (N)
!
!  EST    (input/output) DOUBLE PRECISION
!         On entry with KASE = 1 or 2 and ISAVE(1) = 3, EST should be
!         unchanged from the previous call to DLACN2.
!         On exit, EST is an estimate (a lower bound) for norm(A).
!
!  KASE   (input/output) INTEGER
!         On the initial call to DLACN2, KASE should be 0.
!         On an intermediate return, KASE will be 1 or 2, indicating
!         whether X should be overwritten by A * X  or A' * X.
!         On the final return from DLACN2, KASE will again be 0.
!
!  ISAVE  (input/output) INTEGER array, dimension (3)
!         ISAVE is used to save variables between calls to DLACN2
!
!  Further Details
!  ======= =======
!
!  Contributed by Nick Higham, University of Manchester.
!  Originally named SONEST, dated March 16, 1988.
!
!  Reference: N.J. Higham, "FORTRAN codes for estimating the one-norm of
!  a real or complex matrix, with applications to condition estimation",
!  ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.
!
!  This is a thread safe version of DLACON, which uses the array ISAVE
!  in place of a SAVE statement, as follows:
!
!     DLACON     DLACN2
!      JUMP     ISAVE(1)
!      J        ISAVE(2)
!      ITER     ISAVE(3)
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            ITMAX
      PARAMETER          ( ITMAX = 5 )
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, JLAST
      DOUBLE PRECISION   ALTSGN, ESTOLD, TEMP
!     ..
!     .. External Functions ..
!      INTEGER            IDAMAX
!      DOUBLE PRECISION   DASUM
!      EXTERNAL           IDAMAX, DASUM
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DCOPY
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, NINT, SIGN
!     ..
!     .. Executable Statements ..
!
      IF( KASE.EQ.0 ) THEN
         DO 10 I = 1, N
            X( I ) = ONE / DBLE( N )
   10    CONTINUE
         KASE = 1
         ISAVE( 1 ) = 1
         RETURN
      END IF
!
      GO TO ( 20, 40, 70, 110, 140 )ISAVE( 1 )
!
!     ................ ENTRY   (ISAVE( 1 ) = 1)
!     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
!
   20 CONTINUE
      IF( N.EQ.1 ) THEN
         V( 1 ) = X( 1 )
         EST = ABS( V( 1 ) )
!        ... QUIT
         GO TO 150
      END IF
      EST = DASUM( N, X, 1 )
!
      DO 30 I = 1, N
         X( I ) = SIGN( ONE, X( I ) )
         ISGN( I ) = NINT( X( I ) )
   30 CONTINUE
      KASE = 2
      ISAVE( 1 ) = 2
      RETURN
!
!     ................ ENTRY   (ISAVE( 1 ) = 2)
!     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
!
   40 CONTINUE
      ISAVE( 2 ) = IDAMAX( N, X, 1 )
      ISAVE( 3 ) = 2
!
!     MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
!
   50 CONTINUE
      DO 60 I = 1, N
         X( I ) = ZERO
   60 CONTINUE
      X( ISAVE( 2 ) ) = ONE
      KASE = 1
      ISAVE( 1 ) = 3
      RETURN
!
!     ................ ENTRY   (ISAVE( 1 ) = 3)
!     X HAS BEEN OVERWRITTEN BY A*X.
!
   70 CONTINUE
      CALL DCOPY( N, X, 1, V, 1 )
      ESTOLD = EST
      EST = DASUM( N, V, 1 )
      DO 80 I = 1, N
         IF( NINT( SIGN( ONE, X( I ) ) ).NE.ISGN( I ) ) &
            GO TO 90
   80 CONTINUE
!     REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
      GO TO 120
!
   90 CONTINUE
!     TEST FOR CYCLING.
      IF( EST.LE.ESTOLD ) &
         GO TO 120
!
      DO 100 I = 1, N
         X( I ) = SIGN( ONE, X( I ) )
         ISGN( I ) = NINT( X( I ) )
  100 CONTINUE
      KASE = 2
      ISAVE( 1 ) = 4
      RETURN
!
!     ................ ENTRY   (ISAVE( 1 ) = 4)
!     X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
!
  110 CONTINUE
      JLAST = ISAVE( 2 )
      ISAVE( 2 ) = IDAMAX( N, X, 1 )
      IF( ( X( JLAST ).NE.ABS( X( ISAVE( 2 ) ) ) ) .AND. &
          ( ISAVE( 3 ).LT.ITMAX ) ) THEN
         ISAVE( 3 ) = ISAVE( 3 ) + 1
         GO TO 50
      END IF
!
!     ITERATION COMPLETE.  FINAL STAGE.
!
  120 CONTINUE
      ALTSGN = ONE
      DO 130 I = 1, N
         X( I ) = ALTSGN*( ONE+DBLE( I-1 ) / DBLE( N-1 ) )
         ALTSGN = -ALTSGN
  130 CONTINUE
      KASE = 1
      ISAVE( 1 ) = 5
      RETURN
!
!     ................ ENTRY   (ISAVE( 1 ) = 5)
!     X HAS BEEN OVERWRITTEN BY A*X.
!
  140 CONTINUE
      TEMP = TWO*( DASUM( N, X, 1 ) / DBLE( 3*N ) )
      IF( TEMP.GT.EST ) THEN
         CALL DCOPY( N, X, 1, V, 1 )
         EST = TEMP
      END IF
!
  150 CONTINUE
      KASE = 0
      RETURN
!
!     End of DLACN2
!
      END SUBROUTINE DLACN2


      SUBROUTINE DLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, &
                         ILOZ, IHIZ, Z, LDZ, WORK, LWORK, INFO )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N
      LOGICAL            WANTT, WANTZ
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), WI( * ), WORK( * ), WR( * ), &
                         Z( LDZ, * )
!     ..
!
!     Purpose
!     =======
!
!     DLAQR0 computes the eigenvalues of a Hessenberg matrix H
!     and, optionally, the matrices T and Z from the Schur decomposition
!     H = Z T Z**T, where T is an upper quasi-triangular matrix (the
!     Schur form), and Z is the orthogonal matrix of Schur vectors.
!
!     Optionally Z may be postmultiplied into an input orthogonal
!     matrix Q so that this routine can give the Schur factorization
!     of a matrix A which has been reduced to the Hessenberg form H
!     by the orthogonal matrix Q:  A = Q*H*Q**T = (QZ)*T*(QZ)**T.
!
!     Arguments
!     =========
!
!     WANTT   (input) LOGICAL
!          = .TRUE. : the full Schur form T is required;
!          = .FALSE.: only eigenvalues are required.
!
!     WANTZ   (input) LOGICAL
!          = .TRUE. : the matrix of Schur vectors Z is required;
!          = .FALSE.: Schur vectors are not required.
!
!     N     (input) INTEGER
!           The order of the matrix H.  N .GE. 0.
!
!     ILO   (input) INTEGER
!     IHI   (input) INTEGER
!           It is assumed that H is already upper triangular in rows
!           and columns 1:ILO-1 and IHI+1:N and, if ILO.GT.1,
!           H(ILO,ILO-1) is zero. ILO and IHI are normally set by a
!           previous call to DGEBAL, and then passed to DGEHRD when the
!           matrix output by DGEBAL is reduced to Hessenberg form.
!           Otherwise, ILO and IHI should be set to 1 and N,
!           respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N.
!           If N = 0, then ILO = 1 and IHI = 0.
!
!     H     (input/output) DOUBLE PRECISION array, dimension (LDH,N)
!           On entry, the upper Hessenberg matrix H.
!           On exit, if INFO = 0 and WANTT is .TRUE., then H contains
!           the upper quasi-triangular matrix T from the Schur
!           decomposition (the Schur form); 2-by-2 diagonal blocks
!           (corresponding to complex conjugate pairs of eigenvalues)
!           are returned in standard form, with H(i,i) = H(i+1,i+1)
!           and H(i+1,i)*H(i,i+1).LT.0. If INFO = 0 and WANTT is
!           .FALSE., then the contents of H are unspecified on exit.
!           (The output value of H when INFO.GT.0 is given under the
!           description of INFO below.)
!
!           This subroutine may explicitly set H(i,j) = 0 for i.GT.j and
!           j = 1, 2, ... ILO-1 or j = IHI+1, IHI+2, ... N.
!
!     LDH   (input) INTEGER
!           The leading dimension of the array H. LDH .GE. max(1,N).
!
!     WR    (output) DOUBLE PRECISION array, dimension (IHI)
!     WI    (output) DOUBLE PRECISION array, dimension (IHI)
!           The real and imaginary parts, respectively, of the computed
!           eigenvalues of H(ILO:IHI,ILO:IHI) are stored WR(ILO:IHI)
!           and WI(ILO:IHI). If two eigenvalues are computed as a
!           complex conjugate pair, they are stored in consecutive
!           elements of WR and WI, say the i-th and (i+1)th, with
!           WI(i) .GT. 0 and WI(i+1) .LT. 0. If WANTT is .TRUE., then
!           the eigenvalues are stored in the same order as on the
!           diagonal of the Schur form returned in H, with
!           WR(i) = H(i,i) and, if H(i:i+1,i:i+1) is a 2-by-2 diagonal
!           block, WI(i) = sqrt(-H(i+1,i)*H(i,i+1)) and
!           WI(i+1) = -WI(i).
!
!     ILOZ     (input) INTEGER
!     IHIZ     (input) INTEGER
!           Specify the rows of Z to which transformations must be
!           applied if WANTZ is .TRUE..
!           1 .LE. ILOZ .LE. ILO; IHI .LE. IHIZ .LE. N.
!
!     Z     (input/output) DOUBLE PRECISION array, dimension (LDZ,IHI)
!           If WANTZ is .FALSE., then Z is not referenced.
!           If WANTZ is .TRUE., then Z(ILO:IHI,ILOZ:IHIZ) is
!           replaced by Z(ILO:IHI,ILOZ:IHIZ)*U where U is the
!           orthogonal Schur factor of H(ILO:IHI,ILO:IHI).
!           (The output value of Z when INFO.GT.0 is given under
!           the description of INFO below.)
!
!     LDZ   (input) INTEGER
!           The leading dimension of the array Z.  if WANTZ is .TRUE.
!           then LDZ.GE.MAX(1,IHIZ).  Otherwize, LDZ.GE.1.
!
!     WORK  (workspace/output) DOUBLE PRECISION array, dimension LWORK
!           On exit, if LWORK = -1, WORK(1) returns an estimate of
!           the optimal value for LWORK.
!
!     LWORK (input) INTEGER
!           The dimension of the array WORK.  LWORK .GE. max(1,N)
!           is sufficient, but LWORK typically as large as 6*N may
!           be required for optimal performance.  A workspace query
!           to determine the optimal workspace size is recommended.
!
!           If LWORK = -1, then DLAQR0 does a workspace query.
!           In this case, DLAQR0 checks the input parameters and
!           estimates the optimal workspace size for the given
!           values of N, ILO and IHI.  The estimate is returned
!           in WORK(1).  No error message related to LWORK is
!           issued by XERBLA.  Neither H nor Z are accessed.
!
!
!     INFO  (output) INTEGER
!             =  0:  successful exit
!           .GT. 0:  if INFO = i, DLAQR0 failed to compute all of
!                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR
!                and WI contain those eigenvalues which have been
!                successfully computed.  (Failures are rare.)
!
!                If INFO .GT. 0 and WANT is .FALSE., then on exit,
!                the remaining unconverged eigenvalues are the eigen-
!                values of the upper Hessenberg matrix rows and
!                columns ILO through INFO of the final, output
!                value of H.
!
!                If INFO .GT. 0 and WANTT is .TRUE., then on exit
!
!           (*)  (initial value of H)*U  = U*(final value of H)
!
!                where U is an orthogonal matrix.  The final
!                value of H is upper Hessenberg and quasi-triangular
!                in rows and columns INFO+1 through IHI.
!
!                If INFO .GT. 0 and WANTZ is .TRUE., then on exit
!
!                  (final value of Z(ILO:IHI,ILOZ:IHIZ)
!                   =  (initial value of Z(ILO:IHI,ILOZ:IHIZ)*U
!
!                where U is the orthogonal matrix in (*) (regard-
!                less of the value of WANTT.)
!
!                If INFO .GT. 0 and WANTZ is .FALSE., then Z is not
!                accessed.
!
!
!     ================================================================
!     Based on contributions by
!        Karen Braman and Ralph Byers, Department of Mathematics,
!        University of Kansas, USA
!
!     ================================================================
!
!     References:
!       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
!       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3
!       Performance, SIAM Journal of Matrix Analysis, volume 23, pages
!       929--947, 2002.
!
!       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
!       Algorithm Part II: Aggressive Early Deflation, SIAM Journal
!       of Matrix Analysis, volume 23, pages 948--973, 2002.
!
!     ================================================================
!     .. Parameters ..
!
!     ==== Matrices of order NTINY or smaller must be processed by
!     .    DLAHQR because of insufficient subdiagonal scratch space.
!     .    (This is a hard limit.) ====
!
!     ==== Exceptional deflation windows:  try to cure rare
!     .    slow convergence by increasing the size of the
!     .    deflation window after KEXNW iterations. =====
!
!     ==== Exceptional shifts: try to cure rare slow convergence
!     .    with ad-hoc exceptional shifts every KEXSH iterations.
!     .    The constants WILK1 and WILK2 are used to form the
!     .    exceptional shifts. ====
!
      INTEGER            NTINY
      PARAMETER          ( NTINY = 11 )
      INTEGER            KEXNW, KEXSH
      PARAMETER          ( KEXNW = 5, KEXSH = 6 )
      DOUBLE PRECISION   WILK1, WILK2
      PARAMETER          ( WILK1 = 0.75d0, WILK2 = -0.4375d0 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0d0, ONE = 1.0d0 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   AA, BB, CC, CS, DD, SN, SS, SWAP
      INTEGER            I, INF, IT, ITMAX, K, KACC22, KBOT, KDU, KS, &
                         KT, KTOP, KU, KV, KWH, KWTOP, KWV, LD, LS, &
                         LWKOPT, NDFL, NH, NHO, NIBBLE, NMIN, NS, NSMAX, &
                         NSR, NVE, NW, NWMAX, NWR
      LOGICAL            NWINC, SORTED
      CHARACTER          JBCMPZ*2
!     ..
!     .. External Functions ..
!      INTEGER            ILAENV
!      EXTERNAL           ILAENV
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   ZDUM( 1, 1 )
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DLACPY, DLAHQR, DLANV2, DLAQR3, DLAQR4, DLAQR5
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, MAX, MIN, MOD
!     ..
!     .. Executable Statements ..
      INFO = 0
!
!     ==== Quick return for N = 0: nothing to do. ====
!
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = ONE
         RETURN
      END IF
!
!     ==== Set up job flags for ILAENV. ====
!
      IF( WANTT ) THEN
         JBCMPZ( 1: 1 ) = 'S'
      ELSE
         JBCMPZ( 1: 1 ) = 'E'
      END IF
      IF( WANTZ ) THEN
         JBCMPZ( 2: 2 ) = 'V'
      ELSE
         JBCMPZ( 2: 2 ) = 'N'
      END IF
!
!     ==== Tiny matrices must use DLAHQR. ====
!
      IF( N.LE.NTINY ) THEN
!
!        ==== Estimate optimal workspace. ====
!
         LWKOPT = 1
         IF( LWORK.NE.-1 ) &
            CALL DLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, &
                         ILOZ, IHIZ, Z, LDZ, INFO )
      ELSE
!
!        ==== Use small bulge multi-shift QR with aggressive early
!        .    deflation on larger-than-tiny matrices. ====
!
!        ==== Hope for the best. ====
!
         INFO = 0
!
!        ==== NWR = recommended deflation window size.  At this
!        .    point,  N .GT. NTINY = 11, so there is enough
!        .    subdiagonal workspace for NWR.GE.2 as required.
!        .    (In fact, there is enough subdiagonal space for
!        .    NWR.GE.3.) ====
!
         NWR = ILAENV( 13, 'DLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
         NWR = MAX( 2, NWR )
         NWR = MIN( IHI-ILO+1, ( N-1 ) / 3, NWR )
         NW = NWR
!
!        ==== NSR = recommended number of simultaneous shifts.
!        .    At this point N .GT. NTINY = 11, so there is at
!        .    enough subdiagonal workspace for NSR to be even
!        .    and greater than or equal to two as required. ====
!
         NSR = ILAENV( 15, 'DLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
         NSR = MIN( NSR, ( N+6 ) / 9, IHI-ILO )
         NSR = MAX( 2, NSR-MOD( NSR, 2 ) )
!
!        ==== Estimate optimal workspace ====
!
!        ==== Workspace query call to DLAQR3 ====
!
         CALL DLAQR3( WANTT, WANTZ, N, ILO, IHI, NWR+1, H, LDH, ILOZ, &
                      IHIZ, Z, LDZ, LS, LD, WR, WI, H, LDH, N, H, LDH, &
                      N, H, LDH, WORK, -1 )
!
!        ==== Optimal workspace = MAX(DLAQR5, DLAQR3) ====
!
         LWKOPT = MAX( 3*NSR / 2, INT( WORK( 1 ) ) )
!
!        ==== Quick return in case of workspace query. ====
!
         IF( LWORK.EQ.-1 ) THEN
            WORK( 1 ) = DBLE( LWKOPT )
            RETURN
         END IF
!
!        ==== DLAHQR/DLAQR0 crossover point ====
!
         NMIN = ILAENV( 12, 'DLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
         NMIN = MAX( NTINY, NMIN )
!
!        ==== Nibble crossover point ====
!
         NIBBLE = ILAENV( 14, 'DLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
         NIBBLE = MAX( 0, NIBBLE )
!
!        ==== Accumulate reflections during ttswp?  Use block
!        .    2-by-2 structure during matrix-matrix multiply? ====
!
         KACC22 = ILAENV( 16, 'DLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
         KACC22 = MAX( 0, KACC22 )
         KACC22 = MIN( 2, KACC22 )
!
!        ==== NWMAX = the largest possible deflation window for
!        .    which there is sufficient workspace. ====
!
         NWMAX = MIN( ( N-1 ) / 3, LWORK / 2 )
!
!        ==== NSMAX = the Largest number of simultaneous shifts
!        .    for which there is sufficient workspace. ====
!
         NSMAX = MIN( ( N+6 ) / 9, 2*LWORK / 3 )
         NSMAX = NSMAX - MOD( NSMAX, 2 )
!
!        ==== NDFL: an iteration count restarted at deflation. ====
!
         NDFL = 1
!
!        ==== ITMAX = iteration limit ====
!
         ITMAX = MAX( 30, 2*KEXSH )*MAX( 10, ( IHI-ILO+1 ) )
!
!        ==== Last row and column in the active block ====
!
         KBOT = IHI
!
!        ==== Main Loop ====
!
         DO 80 IT = 1, ITMAX
!
!           ==== Done when KBOT falls below ILO ====
!
            IF( KBOT.LT.ILO ) &
               GO TO 90
!
!           ==== Locate active block ====
!
            DO 10 K = KBOT, ILO + 1, -1
               IF( H( K, K-1 ).EQ.ZERO ) &
                  GO TO 20
   10       CONTINUE
            K = ILO
   20       CONTINUE
            KTOP = K
!
!           ==== Select deflation window size ====
!
            NH = KBOT - KTOP + 1
            IF( NDFL.LT.KEXNW .OR. NH.LT.NW ) THEN
!
!              ==== Typical deflation window.  If possible and
!              .    advisable, nibble the entire active block.
!              .    If not, use size NWR or NWR+1 depending upon
!              .    which has the smaller corresponding subdiagonal
!              .    entry (a heuristic). ====
!
               NWINC = .TRUE.
               IF( NH.LE.MIN( NMIN, NWMAX ) ) THEN
                  NW = NH
               ELSE
                  NW = MIN( NWR, NH, NWMAX )
                  IF( NW.LT.NWMAX ) THEN
                     IF( NW.GE.NH-1 ) THEN
                        NW = NH
                     ELSE
                        KWTOP = KBOT - NW + 1
                        IF( ABS( H( KWTOP, KWTOP-1 ) ).GT. &
                            ABS( H( KWTOP-1, KWTOP-2 ) ) )NW = NW + 1
                     END IF
                  END IF
               END IF
            ELSE
!
!              ==== Exceptional deflation window.  If there have
!              .    been no deflations in KEXNW or more iterations,
!              .    then vary the deflation window size.   At first,
!              .    because, larger windows are, in general, more
!              .    powerful than smaller ones, rapidly increase the
!              .    window up to the maximum reasonable and possible.
!              .    Then maybe try a slightly smaller window.  ====
!
               IF( NWINC .AND. NW.LT.MIN( NWMAX, NH ) ) THEN
                  NW = MIN( NWMAX, NH, 2*NW )
               ELSE
                  NWINC = .FALSE.
                  IF( NW.EQ.NH .AND. NH.GT.2 ) &
                     NW = NH - 1
               END IF
            END IF
!
!           ==== Aggressive early deflation:
!           .    split workspace under the subdiagonal into
!           .      - an nw-by-nw work array V in the lower
!           .        left-hand-corner,
!           .      - an NW-by-at-least-NW-but-more-is-better
!           .        (NW-by-NHO) horizontal work array along
!           .        the bottom edge,
!           .      - an at-least-NW-but-more-is-better (NHV-by-NW)
!           .        vertical work array along the left-hand-edge.
!           .        ====
!
            KV = N - NW + 1
            KT = NW + 1
            NHO = ( N-NW-1 ) - KT + 1
            KWV = NW + 2
            NVE = ( N-NW ) - KWV + 1
!
!           ==== Aggressive early deflation ====
!
            CALL DLAQR3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ, &
                         IHIZ, Z, LDZ, LS, LD, WR, WI, H( KV, 1 ), LDH, &
                         NHO, H( KV, KT ), LDH, NVE, H( KWV, 1 ), LDH, &
                         WORK, LWORK )
!
!           ==== Adjust KBOT accounting for new deflations. ====
!
            KBOT = KBOT - LD
!
!           ==== KS points to the shifts. ====
!
            KS = KBOT - LS + 1
!
!           ==== Skip an expensive QR sweep if there is a (partly
!           .    heuristic) reason to expect that many eigenvalues
!           .    will deflate without it.  Here, the QR sweep is
!           .    skipped if many eigenvalues have just been deflated
!           .    or if the remaining active block is small.
!
            IF( ( LD.EQ.0 ) .OR. ( ( 100*LD.LE.NW*NIBBLE ) .AND. ( KBOT- &
                KTOP+1.GT.MIN( NMIN, NWMAX ) ) ) ) THEN
!
!              ==== NS = nominal number of simultaneous shifts.
!              .    This may be lowered (slightly) if DLAQR3
!              .    did not provide that many shifts. ====
!
               NS = MIN( NSMAX, NSR, MAX( 2, KBOT-KTOP ) )
               NS = NS - MOD( NS, 2 )
!
!              ==== If there have been no deflations
!              .    in a multiple of KEXSH iterations,
!              .    then try exceptional shifts.
!              .    Otherwise use shifts provided by
!              .    DLAQR3 above or from the eigenvalues
!              .    of a trailing principal submatrix. ====
!
               IF( MOD( NDFL, KEXSH ).EQ.0 ) THEN
                  KS = KBOT - NS + 1
                  DO 30 I = KBOT, MAX( KS+1, KTOP+2 ), -2
                     SS = ABS( H( I, I-1 ) ) + ABS( H( I-1, I-2 ) )
                     AA = WILK1*SS + H( I, I )
                     BB = SS
                     CC = WILK2*SS
                     DD = AA
                     CALL DLANV2( AA, BB, CC, DD, WR( I-1 ), WI( I-1 ), &
                                  WR( I ), WI( I ), CS, SN )
   30             CONTINUE
                  IF( KS.EQ.KTOP ) THEN
                     WR( KS+1 ) = H( KS+1, KS+1 )
                     WI( KS+1 ) = ZERO
                     WR( KS ) = WR( KS+1 )
                     WI( KS ) = WI( KS+1 )
                  END IF
               ELSE
!
!                 ==== Got NS/2 or fewer shifts? Use DLAQR4 or
!                 .    DLAHQR on a trailing principal submatrix to
!                 .    get more. (Since NS.LE.NSMAX.LE.(N+6)/9,
!                 .    there is enough space below the subdiagonal
!                 .    to fit an NS-by-NS scratch array.) ====
!
                  IF( KBOT-KS+1.LE.NS / 2 ) THEN
                     KS = KBOT - NS + 1
                     KT = N - NS + 1
                     CALL DLACPY( 'A', NS, NS, H( KS, KS ), LDH, &
                                  H( KT, 1 ), LDH )
                     IF( NS.GT.NMIN ) THEN
                        CALL DLAQR4( .false., .false., NS, 1, NS, &
                                     H( KT, 1 ), LDH, WR( KS ), &
                                     WI( KS ), 1, 1, ZDUM, 1, WORK, &
                                     LWORK, INF )
                     ELSE
                        CALL DLAHQR( .false., .false., NS, 1, NS, &
                                     H( KT, 1 ), LDH, WR( KS ), &
                                     WI( KS ), 1, 1, ZDUM, 1, INF )
                     END IF
                     KS = KS + INF
!
!                    ==== In case of a rare QR failure use
!                    .    eigenvalues of the trailing 2-by-2
!                    .    principal submatrix.  ====
!
                     IF( KS.GE.KBOT ) THEN
                        AA = H( KBOT-1, KBOT-1 )
                        CC = H( KBOT, KBOT-1 )
                        BB = H( KBOT-1, KBOT )
                        DD = H( KBOT, KBOT )
                        CALL DLANV2( AA, BB, CC, DD, WR( KBOT-1 ), &
                                     WI( KBOT-1 ), WR( KBOT ), &
                                     WI( KBOT ), CS, SN )
                        KS = KBOT - 1
                     END IF
                  END IF
!
                  IF( KBOT-KS+1.GT.NS ) THEN
!
!                    ==== Sort the shifts (Helps a little)
!                    .    Bubble sort keeps complex conjugate
!                    .    pairs together. ====
!
                     SORTED = .false.
                     DO 50 K = KBOT, KS + 1, -1
                        IF( SORTED ) &
                           GO TO 60
                        SORTED = .true.
                        DO 40 I = KS, K - 1
                           IF( ABS( WR( I ) )+ABS( WI( I ) ).LT. &
                               ABS( WR( I+1 ) )+ABS( WI( I+1 ) ) ) THEN
                              SORTED = .false.
!
                              SWAP = WR( I )
                              WR( I ) = WR( I+1 )
                              WR( I+1 ) = SWAP
!
                              SWAP = WI( I )
                              WI( I ) = WI( I+1 )
                              WI( I+1 ) = SWAP
                           END IF
   40                   CONTINUE
   50                CONTINUE
   60                CONTINUE
                  END IF
!
!                 ==== Shuffle shifts into pairs of real shifts
!                 .    and pairs of complex conjugate shifts
!                 .    assuming complex conjugate shifts are
!                 .    already adjacent to one another. (Yes,
!                 .    they are.)  ====
!
                  DO 70 I = KBOT, KS + 2, -2
                     IF( WI( I ).NE.-WI( I-1 ) ) THEN
!
                        SWAP = WR( I )
                        WR( I ) = WR( I-1 )
                        WR( I-1 ) = WR( I-2 )
                        WR( I-2 ) = SWAP
!
                        SWAP = WI( I )
                        WI( I ) = WI( I-1 )
                        WI( I-1 ) = WI( I-2 )
                        WI( I-2 ) = SWAP
                     END IF
   70             CONTINUE
               END IF
!
!              ==== If there are only two shifts and both are
!              .    real, then use only one.  ====
!
               IF( KBOT-KS+1.EQ.2 ) THEN
                  IF( WI( KBOT ).EQ.ZERO ) THEN
                     IF( ABS( WR( KBOT )-H( KBOT, KBOT ) ).LT. &
                         ABS( WR( KBOT-1 )-H( KBOT, KBOT ) ) ) THEN
                        WR( KBOT-1 ) = WR( KBOT )
                     ELSE
                        WR( KBOT ) = WR( KBOT-1 )
                     END IF
                  END IF
               END IF
!
!              ==== Use up to NS of the the smallest magnatiude
!              .    shifts.  If there aren't NS shifts available,
!              .    then use them all, possibly dropping one to
!              .    make the number of shifts even. ====
!
               NS = MIN( NS, KBOT-KS+1 )
               NS = NS - MOD( NS, 2 )
               KS = KBOT - NS + 1
!
!              ==== Small-bulge multi-shift QR sweep:
!              .    split workspace under the subdiagonal into
!              .    - a KDU-by-KDU work array U in the lower
!              .      left-hand-corner,
!              .    - a KDU-by-at-least-KDU-but-more-is-better
!              .      (KDU-by-NHo) horizontal work array WH along
!              .      the bottom edge,
!              .    - and an at-least-KDU-but-more-is-better-by-KDU
!              .      (NVE-by-KDU) vertical work WV arrow along
!              .      the left-hand-edge. ====
!
               KDU = 3*NS - 3
               KU = N - KDU + 1
               KWH = KDU + 1
               NHO = ( N-KDU+1-4 ) - ( KDU+1 ) + 1
               KWV = KDU + 4
               NVE = N - KDU - KWV + 1
!
!              ==== Small-bulge multi-shift QR sweep ====
!
               CALL DLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NS, &
                            WR( KS ), WI( KS ), H, LDH, ILOZ, IHIZ, Z, &
                            LDZ, WORK, 3, H( KU, 1 ), LDH, NVE, &
                            H( KWV, 1 ), LDH, NHO, H( KU, KWH ), LDH )
            END IF
!
!           ==== Note progress (or the lack of it). ====
!
            IF( LD.GT.0 ) THEN
               NDFL = 1
            ELSE
               NDFL = NDFL + 1
            END IF
!
!           ==== End of main loop ====
   80    CONTINUE
!
!        ==== Iteration limit exceeded.  Set INFO to show where
!        .    the problem occurred and exit. ====
!
         INFO = KBOT
   90    CONTINUE
      END IF
!
!     ==== Return the optimal value of LWORK. ====
!
      WORK( 1 ) = DBLE( LWKOPT )
!
!     ==== End of DLAQR0 ====
!
      END SUBROUTINE DLAQR0


      SUBROUTINE DLAHR2( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            K, LDA, LDT, LDY, N, NB
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION  A( LDA, * ), T( LDT, NB ), TAU( NB ), &
                         Y( LDY, NB )
!     ..
!
!  Purpose
!  =======
!
!  DLAHR2 reduces the first NB columns of A real general n-BY-(n-k+1)
!  matrix A so that elements below the k-th subdiagonal are zero. The
!  reduction is performed by an orthogonal similarity transformation
!  Q' * A * Q. The routine returns the matrices V and T which determine
!  Q as a block reflector I - V*T*V', and also the matrix Y = A * V * T.
!
!  This is an auxiliary routine called by DGEHRD.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix A.
!
!  K       (input) INTEGER
!          The offset for the reduction. Elements below the k-th
!          subdiagonal in the first NB columns are reduced to zero.
!          K < N.
!
!  NB      (input) INTEGER
!          The number of columns to be reduced.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N-K+1)
!          On entry, the n-by-(n-k+1) general matrix A.
!          On exit, the elements on and above the k-th subdiagonal in
!          the first NB columns are overwritten with the corresponding
!          elements of the reduced matrix; the elements below the k-th
!          subdiagonal, with the array TAU, represent the matrix Q as a
!          product of elementary reflectors. The other columns of A are
!          unchanged. See Further Details.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  TAU     (output) DOUBLE PRECISION array, dimension (NB)
!          The scalar factors of the elementary reflectors. See Further
!          Details.
!
!  T       (output) DOUBLE PRECISION array, dimension (LDT,NB)
!          The upper triangular matrix T.
!
!  LDT     (input) INTEGER
!          The leading dimension of the array T.  LDT >= NB.
!
!  Y       (output) DOUBLE PRECISION array, dimension (LDY,NB)
!          The n-by-nb matrix Y.
!
!  LDY     (input) INTEGER
!          The leading dimension of the array Y. LDY >= N.
!
!  Further Details
!  ===============
!
!  The matrix Q is represented as a product of nb elementary reflectors
!
!     Q = H(1) H(2) . . . H(nb).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i+k-1) = 0, v(i+k) = 1; v(i+k+1:n) is stored on exit in
!  A(i+k+1:n,i), and tau in TAU(i).
!
!  The elements of the vectors v together form the (n-k+1)-by-nb matrix
!  V which is needed, with T and Y, to apply the transformation to the
!  unreduced part of the matrix, using an update of the form:
!  A := (I - V*T*V') * (A - Y*V').
!
!  The contents of A on exit are illustrated by the following example
!  with n = 7, k = 3 and nb = 2:
!
!     ( a   a   a   a   a )
!     ( a   a   a   a   a )
!     ( a   a   a   a   a )
!     ( h   h   a   a   a )
!     ( v1  h   a   a   a )
!     ( v1  v2  a   a   a )
!     ( v1  v2  a   a   a )
!
!  where a denotes an element of the original matrix A, h denotes a
!  modified element of the upper Hessenberg matrix H, and vi denotes an
!  element of the vector defining H(i).
!
!  This file is a slight modification of LAPACK-3.0's DLAHRD
!  incorporating improvements proposed by Quintana-Orti and Van de
!  Gejin. Note that the entries of A(1:K,2:NB) differ from those
!  returned by the original LAPACK routine. This function is
!  not backward compatible with LAPACK3.0.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, &
                           ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION  EI
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DAXPY, DCOPY, DGEMM, DGEMV, DLACPY, &
!                         DLARFG, DSCAL, DTRMM, DTRMV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( N.LE.1 ) &
         RETURN
!
      DO 10 I = 1, NB
         IF( I.GT.1 ) THEN
!
!           Update A(K+1:N,I)
!
!           Update I-th column of A - Y * V'
!
            CALL DGEMV( 'NO TRANSPOSE', N-K, I-1, -ONE, Y(K+1,1), LDY, &
                        A( K+I-1, 1 ), LDA, ONE, A( K+1, I ), 1 )
!
!           Apply I - V * T' * V' to this column (call it b) from the
!           left, using the last column of T as workspace
!
!           Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows)
!                    ( V2 )             ( b2 )
!
!           where V1 is unit lower triangular
!
!           w := V1' * b1
!
            CALL DCOPY( I-1, A( K+1, I ), 1, T( 1, NB ), 1 )
            CALL DTRMV( 'Lower', 'Transpose', 'UNIT', &
                        I-1, A( K+1, 1 ), &
                        LDA, T( 1, NB ), 1 )
!
!           w := w + V2'*b2
!
            CALL DGEMV( 'Transpose', N-K-I+1, I-1, &
                        ONE, A( K+I, 1 ), &
                        LDA, A( K+I, I ), 1, ONE, T( 1, NB ), 1 )
!
!           w := T'*w
!
            CALL DTRMV( 'Upper', 'Transpose', 'NON-UNIT', &
                        I-1, T, LDT, &
                        T( 1, NB ), 1 )
!
!           b2 := b2 - V2*w
!
            CALL DGEMV( 'NO TRANSPOSE', N-K-I+1, I-1, -ONE, &
                        A( K+I, 1 ), &
                        LDA, T( 1, NB ), 1, ONE, A( K+I, I ), 1 )
!
!           b1 := b1 - V1*w
!
            CALL DTRMV( 'Lower', 'NO TRANSPOSE', &
                        'UNIT', I-1, &
                        A( K+1, 1 ), LDA, T( 1, NB ), 1 )
            CALL DAXPY( I-1, -ONE, T( 1, NB ), 1, A( K+1, I ), 1 )
!
            A( K+I-1, I-1 ) = EI
         END IF
!
!        Generate the elementary reflector H(I) to annihilate
!        A(K+I+1:N,I)
!
         CALL DLARFG( N-K-I+1, A( K+I, I ), A( MIN( K+I+1, N ), I ), 1, &
                      TAU( I ) )
         EI = A( K+I, I )
         A( K+I, I ) = ONE
!
!        Compute  Y(K+1:N,I)
!
         CALL DGEMV( 'NO TRANSPOSE', N-K, N-K-I+1, &
                     ONE, A( K+1, I+1 ), &
                     LDA, A( K+I, I ), 1, ZERO, Y( K+1, I ), 1 )
         CALL DGEMV( 'Transpose', N-K-I+1, I-1, &
                     ONE, A( K+I, 1 ), LDA, &
                     A( K+I, I ), 1, ZERO, T( 1, I ), 1 )
         CALL DGEMV( 'NO TRANSPOSE', N-K, I-1, -ONE, &
                     Y( K+1, 1 ), LDY, &
                     T( 1, I ), 1, ONE, Y( K+1, I ), 1 )
         CALL DSCAL( N-K, TAU( I ), Y( K+1, I ), 1 )
!
!        Compute T(1:I,I)
!
         CALL DSCAL( I-1, -TAU( I ), T( 1, I ), 1 )
         CALL DTRMV( 'Upper', 'No Transpose', 'NON-UNIT', &
                     I-1, T, LDT, &
                     T( 1, I ), 1 )
         T( I, I ) = TAU( I )
!
   10 CONTINUE
      A( K+NB, NB ) = EI
!
!     Compute Y(1:K,1:NB)
!
      CALL DLACPY( 'ALL', K, NB, A( 1, 2 ), LDA, Y, LDY )
      CALL DTRMM( 'RIGHT', 'Lower', 'NO TRANSPOSE', &
                  'UNIT', K, NB, &
                  ONE, A( K+1, 1 ), LDA, Y, LDY )
      IF( N.GT.K+NB ) &
         CALL DGEMM( 'NO TRANSPOSE', 'NO TRANSPOSE', K, &
                     NB, N-K-NB, ONE, &
                     A( 1, 2+NB ), LDA, A( K+1+NB, 1 ), LDA, ONE, Y, &
                     LDY )
      CALL DTRMM( 'RIGHT', 'Upper', 'NO TRANSPOSE', &
                  'NON-UNIT', K, NB, &
                  ONE, T, LDT, Y, LDY )
!
      RETURN
!
!     End of DLAHR2
!
      END SUBROUTINE DLAHR2


      SUBROUTINE DLAQR3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ, &
                         IHIZ, Z, LDZ, NS, ND, SR, SI, V, LDV, NH, T, &
                         LDT, NV, WV, LDWV, WORK, LWORK )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV, &
                         LDZ, LWORK, N, ND, NH, NS, NV, NW
      LOGICAL            WANTT, WANTZ
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), SI( * ), SR( * ), T( LDT, * ), &
                         V( LDV, * ), WORK( * ), WV( LDWV, * ), &
                         Z( LDZ, * )
!     ..
!
!     ******************************************************************
!     Aggressive early deflation:
!
!     This subroutine accepts as input an upper Hessenberg matrix
!     H and performs an orthogonal similarity transformation
!     designed to detect and deflate fully converged eigenvalues from
!     a trailing principal submatrix.  On output H has been over-
!     written by a new Hessenberg matrix that is a perturbation of
!     an orthogonal similarity transformation of H.  It is to be
!     hoped that the final version of H has many zero subdiagonal
!     entries.
!
!     ******************************************************************
!     WANTT   (input) LOGICAL
!          If .TRUE., then the Hessenberg matrix H is fully updated
!          so that the quasi-triangular Schur factor may be
!          computed (in cooperation with the calling subroutine).
!          If .FALSE., then only enough of H is updated to preserve
!          the eigenvalues.
!
!     WANTZ   (input) LOGICAL
!          If .TRUE., then the orthogonal matrix Z is updated so
!          so that the orthogonal Schur factor may be computed
!          (in cooperation with the calling subroutine).
!          If .FALSE., then Z is not referenced.
!
!     N       (input) INTEGER
!          The order of the matrix H and (if WANTZ is .TRUE.) the
!          order of the orthogonal matrix Z.
!
!     KTOP    (input) INTEGER
!          It is assumed that either KTOP = 1 or H(KTOP,KTOP-1)=0.
!          KBOT and KTOP together determine an isolated block
!          along the diagonal of the Hessenberg matrix.
!
!     KBOT    (input) INTEGER
!          It is assumed without a check that either
!          KBOT = N or H(KBOT+1,KBOT)=0.  KBOT and KTOP together
!          determine an isolated block along the diagonal of the
!          Hessenberg matrix.
!
!     NW      (input) INTEGER
!          Deflation window size.  1 .LE. NW .LE. (KBOT-KTOP+1).
!
!     H       (input/output) DOUBLE PRECISION array, dimension (LDH,N)
!          On input the initial N-by-N section of H stores the
!          Hessenberg matrix undergoing aggressive early deflation.
!          On output H has been transformed by an orthogonal
!          similarity transformation, perturbed, and the returned
!          to Hessenberg form that (it is to be hoped) has some
!          zero subdiagonal entries.
!
!     LDH     (input) integer
!          Leading dimension of H just as declared in the calling
!          subroutine.  N .LE. LDH
!
!     ILOZ    (input) INTEGER
!     IHIZ    (input) INTEGER
!          Specify the rows of Z to which transformations must be
!          applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N.
!
!     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,IHI)
!          IF WANTZ is .TRUE., then on output, the orthogonal
!          similarity transformation mentioned above has been
!          accumulated into Z(ILOZ:IHIZ,ILO:IHI) from the right.
!          If WANTZ is .FALSE., then Z is unreferenced.
!
!     LDZ     (input) integer
!          The leading dimension of Z just as declared in the
!          calling subroutine.  1 .LE. LDZ.
!
!     NS      (output) integer
!          The number of unconverged (ie approximate) eigenvalues
!          returned in SR and SI that may be used as shifts by the
!          calling subroutine.
!
!     ND      (output) integer
!          The number of converged eigenvalues uncovered by this
!          subroutine.
!
!     SR      (output) DOUBLE PRECISION array, dimension KBOT
!     SI      (output) DOUBLE PRECISION array, dimension KBOT
!          On output, the real and imaginary parts of approximate
!          eigenvalues that may be used for shifts are stored in
!          SR(KBOT-ND-NS+1) through SR(KBOT-ND) and
!          SI(KBOT-ND-NS+1) through SI(KBOT-ND), respectively.
!          The real and imaginary parts of converged eigenvalues
!          are stored in SR(KBOT-ND+1) through SR(KBOT) and
!          SI(KBOT-ND+1) through SI(KBOT), respectively.
!
!     V       (workspace) DOUBLE PRECISION array, dimension (LDV,NW)
!          An NW-by-NW work array.
!
!     LDV     (input) integer scalar
!          The leading dimension of V just as declared in the
!          calling subroutine.  NW .LE. LDV
!
!     NH      (input) integer scalar
!          The number of columns of T.  NH.GE.NW.
!
!     T       (workspace) DOUBLE PRECISION array, dimension (LDT,NW)
!
!     LDT     (input) integer
!          The leading dimension of T just as declared in the
!          calling subroutine.  NW .LE. LDT
!
!     NV      (input) integer
!          The number of rows of work array WV available for
!          workspace.  NV.GE.NW.
!
!     WV      (workspace) DOUBLE PRECISION array, dimension (LDWV,NW)
!
!     LDWV    (input) integer
!          The leading dimension of W just as declared in the
!          calling subroutine.  NW .LE. LDV
!
!     WORK    (workspace) DOUBLE PRECISION array, dimension LWORK.
!          On exit, WORK(1) is set to an estimate of the optimal value
!          of LWORK for the given values of N, NW, KTOP and KBOT.
!
!     LWORK   (input) integer
!          The dimension of the work array WORK.  LWORK = 2*NW
!          suffices, but greater efficiency may result from larger
!          values of LWORK.
!
!          If LWORK = -1, then a workspace query is assumed; DLAQR3
!          only estimates the optimal workspace size for the given
!          values of N, NW, KTOP and KBOT.  The estimate is returned
!          in WORK(1).  No error message related to LWORK is issued
!          by XERBLA.  Neither H nor Z are accessed.
!
!     ================================================================
!     Based on contributions by
!        Karen Braman and Ralph Byers, Department of Mathematics,
!        University of Kansas, USA
!
!     ==================================================================
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0d0, ONE = 1.0d0 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   AA, BB, BETA, CC, CS, DD, EVI, EVK, FOO, S, &
                         SAFMAX, SAFMIN, SMLNUM, SN, TAU, ULP
      INTEGER            I, IFST, ILST, INFO, INFQR, J, JW, K, KCOL, &
                         KEND, KLN, KROW, KWTOP, LTOP, LWK1, LWK2, LWK3, &
                         LWKOPT, NMIN
      LOGICAL            BULGE, SORTED
!     ..
!     .. External Functions ..
!      DOUBLE PRECISION   DLAMCH
!      INTEGER            ILAENV
!      EXTERNAL           DLAMCH, ILAENV
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DCOPY, DGEHRD, DGEMM, DLABAD, DLACPY, DLAHQR, &
!                         DLANV2, DLAQR4, DLARF, DLARFG, DLASET, DORGHR, &
!                         DTREXC
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
!     ==== Estimate optimal workspace. ====
!
      JW = MIN( NW, KBOT-KTOP+1 )
      IF( JW.LE.2 ) THEN
         LWKOPT = 1
      ELSE
!
!        ==== Workspace query call to DGEHRD ====
!
         CALL DGEHRD( JW, 1, JW-1, T, LDT, WORK, WORK, -1, INFO )
         LWK1 = INT( WORK( 1 ) )
!
!        ==== Workspace query call to DORGHR ====
!
         CALL DORGHR( JW, 1, JW-1, T, LDT, WORK, WORK, -1, INFO )
         LWK2 = INT( WORK( 1 ) )
!
!        ==== Workspace query call to DLAQR4 ====
!
         CALL DLAQR4( .true., .true., JW, 1, JW, T, LDT, SR, SI, 1, JW, &
                      V, LDV, WORK, -1, INFQR )
         LWK3 = INT( WORK( 1 ) )
!
!        ==== Optimal workspace ====
!
         LWKOPT = MAX( JW+MAX( LWK1, LWK2 ), LWK3 )
      END IF
!
!     ==== Quick return in case of workspace query. ====
!
      IF( LWORK.EQ.-1 ) THEN
         WORK( 1 ) = DBLE( LWKOPT )
         RETURN
      END IF
!
!     ==== Nothing to do ...
!     ... for an empty active block ... ====
      NS = 0
      ND = 0
      IF( KTOP.GT.KBOT ) &
         RETURN
!     ... nor for an empty deflation window. ====
      IF( NW.LT.1 ) &
         RETURN
!
!     ==== Machine constants ====
!
      SAFMIN = DLAMCH( 'SAFE MINIMUM' )
      SAFMAX = ONE / SAFMIN
      CALL DLABAD( SAFMIN, SAFMAX )
      ULP = DLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( DBLE( N ) / ULP )
!
!     ==== Setup deflation window ====
!
      JW = MIN( NW, KBOT-KTOP+1 )
      KWTOP = KBOT - JW + 1
      IF( KWTOP.EQ.KTOP ) THEN
         S = ZERO
      ELSE
         S = H( KWTOP, KWTOP-1 )
      END IF
!
      IF( KBOT.EQ.KWTOP ) THEN
!
!        ==== 1-by-1 deflation window: not much to do ====
!
         SR( KWTOP ) = H( KWTOP, KWTOP )
         SI( KWTOP ) = ZERO
         NS = 1
         ND = 0
         IF( ABS( S ).LE.MAX( SMLNUM, ULP*ABS( H( KWTOP, KWTOP ) ) ) ) &
              THEN
            NS = 0
            ND = 1
            IF( KWTOP.GT.KTOP ) &
               H( KWTOP, KWTOP-1 ) = ZERO
         END IF
         RETURN
      END IF
!
!     ==== Convert to spike-triangular form.  (In case of a
!     .    rare QR failure, this routine continues to do
!     .    aggressive early deflation using that part of
!     .    the deflation window that converged using INFQR
!     .    here and there to keep track.) ====
!
      CALL DLACPY( 'U', JW, JW, H( KWTOP, KWTOP ), LDH, T, LDT )
      CALL DCOPY( JW-1, H( KWTOP+1, KWTOP ), LDH+1, T( 2, 1 ), LDT+1 )
!
      CALL DLASET( 'A', JW, JW, ZERO, ONE, V, LDV )
      NMIN = ILAENV( 12, 'DLAQR3', 'SV', JW, 1, JW, LWORK )
      IF( JW.GT.NMIN ) THEN
         CALL DLAQR4( .true., .true., JW, 1, JW, T, LDT, SR( KWTOP ), &
                      SI( KWTOP ), 1, JW, V, LDV, WORK, LWORK, INFQR )
      ELSE
         CALL DLAHQR( .true., .true., JW, 1, JW, T, LDT, SR( KWTOP ), &
                      SI( KWTOP ), 1, JW, V, LDV, INFQR )
      END IF
!
!     ==== DTREXC needs a clean margin near the diagonal ====
!
      DO 10 J = 1, JW - 3
         T( J+2, J ) = ZERO
         T( J+3, J ) = ZERO
   10 CONTINUE
      IF( JW.GT.2 ) &
         T( JW, JW-2 ) = ZERO
!
!     ==== Deflation detection loop ====
!
      NS = JW
      ILST = INFQR + 1
   20 CONTINUE
      IF( ILST.LE.NS ) THEN
         IF( NS.EQ.1 ) THEN
            BULGE = .FALSE.
         ELSE
            BULGE = T( NS, NS-1 ).NE.ZERO
         END IF
!
!        ==== Small spike tip test for deflation ====
!
         IF( .NOT.BULGE ) THEN
!
!           ==== Real eigenvalue ====
!
            FOO = ABS( T( NS, NS ) )
            IF( FOO.EQ.ZERO ) &
               FOO = ABS( S )
            IF( ABS( S*V( 1, NS ) ).LE.MAX( SMLNUM, ULP*FOO ) ) THEN
!
!              ==== Deflatable ====
!
               NS = NS - 1
            ELSE
!
!              ==== Undeflatable.   Move it up out of the way.
!              .    (DTREXC can not fail in this case.) ====
!
               IFST = NS
               CALL DTREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST, WORK, &
                            INFO )
               ILST = ILST + 1
            END IF
         ELSE
!
!           ==== Complex conjugate pair ====
!
            FOO = ABS( T( NS, NS ) ) + SQRT( ABS( T( NS, NS-1 ) ) )* &
                  SQRT( ABS( T( NS-1, NS ) ) )
            IF( FOO.EQ.ZERO ) &
               FOO = ABS( S )
            IF( MAX( ABS( S*V( 1, NS ) ), ABS( S*V( 1, NS-1 ) ) ).LE. &
                MAX( SMLNUM, ULP*FOO ) ) THEN
!
!              ==== Deflatable ====
!
               NS = NS - 2
            ELSE
!
!              ==== Undflatable. Move them up out of the way.
!              .    Fortunately, DTREXC does the right thing with
!              .    ILST in case of a rare exchange failure. ====
!
               IFST = NS
               CALL DTREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST, WORK, &
                            INFO )
               ILST = ILST + 2
            END IF
         END IF
!
!        ==== End deflation detection loop ====
!
         GO TO 20
      END IF
!
!        ==== Return to Hessenberg form ====
!
      IF( NS.EQ.0 ) &
         S = ZERO
!
      IF( NS.LT.JW ) THEN
!
!        ==== sorting diagonal blocks of T improves accuracy for
!        .    graded matrices.  Bubble sort deals well with
!        .    exchange failures. ====
!
         SORTED = .false.
         I = NS + 1
   30    CONTINUE
         IF( SORTED ) &
            GO TO 50
         SORTED = .true.
!
         KEND = I - 1
         I = INFQR + 1
         IF( I.EQ.NS ) THEN
            K = I + 1
         ELSE IF( T( I+1, I ).EQ.ZERO ) THEN
            K = I + 1
         ELSE
            K = I + 2
         END IF
   40    CONTINUE
         IF( K.LE.KEND ) THEN
            IF( K.EQ.I+1 ) THEN
               EVI = ABS( T( I, I ) )
            ELSE
               EVI = ABS( T( I, I ) ) + SQRT( ABS( T( I+1, I ) ) )* &
                     SQRT( ABS( T( I, I+1 ) ) )
            END IF
!
            IF( K.EQ.KEND ) THEN
               EVK = ABS( T( K, K ) )
            ELSE IF( T( K+1, K ).EQ.ZERO ) THEN
               EVK = ABS( T( K, K ) )
            ELSE
               EVK = ABS( T( K, K ) ) + SQRT( ABS( T( K+1, K ) ) )* &
                     SQRT( ABS( T( K, K+1 ) ) )
            END IF
!
            IF( EVI.GE.EVK ) THEN
               I = K
            ELSE
               SORTED = .false.
               IFST = I
               ILST = K
               CALL DTREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST, WORK, &
                            INFO )
               IF( INFO.EQ.0 ) THEN
                  I = ILST
               ELSE
                  I = K
               END IF
            END IF
            IF( I.EQ.KEND ) THEN
               K = I + 1
            ELSE IF( T( I+1, I ).EQ.ZERO ) THEN
               K = I + 1
            ELSE
               K = I + 2
            END IF
            GO TO 40
         END IF
         GO TO 30
   50    CONTINUE
      END IF
!
!     ==== Restore shift/eigenvalue array from T ====
!
      I = JW
   60 CONTINUE
      IF( I.GE.INFQR+1 ) THEN
         IF( I.EQ.INFQR+1 ) THEN
            SR( KWTOP+I-1 ) = T( I, I )
            SI( KWTOP+I-1 ) = ZERO
            I = I - 1
         ELSE IF( T( I, I-1 ).EQ.ZERO ) THEN
            SR( KWTOP+I-1 ) = T( I, I )
            SI( KWTOP+I-1 ) = ZERO
            I = I - 1
         ELSE
            AA = T( I-1, I-1 )
            CC = T( I, I-1 )
            BB = T( I-1, I )
            DD = T( I, I )
            CALL DLANV2( AA, BB, CC, DD, SR( KWTOP+I-2 ), &
                         SI( KWTOP+I-2 ), SR( KWTOP+I-1 ), &
                         SI( KWTOP+I-1 ), CS, SN )
            I = I - 2
         END IF
         GO TO 60
      END IF
!
      IF( NS.LT.JW .OR. S.EQ.ZERO ) THEN
         IF( NS.GT.1 .AND. S.NE.ZERO ) THEN
!
!           ==== Reflect spike back into lower triangle ====
!
            CALL DCOPY( NS, V, LDV, WORK, 1 )
            BETA = WORK( 1 )
            CALL DLARFG( NS, BETA, WORK( 2 ), 1, TAU )
            WORK( 1 ) = ONE
!
            CALL DLASET( 'L', JW-2, JW-2, ZERO, ZERO, T( 3, 1 ), LDT )
!
            CALL DLARF( 'L', NS, JW, WORK, 1, TAU, T, LDT, &
                        WORK( JW+1 ) )
            CALL DLARF( 'R', NS, NS, WORK, 1, TAU, T, LDT, &
                        WORK( JW+1 ) )
            CALL DLARF( 'R', JW, NS, WORK, 1, TAU, V, LDV, &
                        WORK( JW+1 ) )
!
            CALL DGEHRD( JW, 1, NS, T, LDT, WORK, WORK( JW+1 ), &
                         LWORK-JW, INFO )
         END IF
!
!        ==== Copy updated reduced window into place ====
!
         IF( KWTOP.GT.1 ) &
            H( KWTOP, KWTOP-1 ) = S*V( 1, 1 )
         CALL DLACPY( 'U', JW, JW, T, LDT, H( KWTOP, KWTOP ), LDH )
         CALL DCOPY( JW-1, T( 2, 1 ), LDT+1, H( KWTOP+1, KWTOP ), &
                     LDH+1 )
!
!        ==== Accumulate orthogonal matrix in order update
!        .    H and Z, if requested.  (A modified version
!        .    of  DORGHR that accumulates block Householder
!        .    transformations into V directly might be
!        .    marginally more efficient than the following.) ====
!
         IF( NS.GT.1 .AND. S.NE.ZERO ) THEN
            CALL DORGHR( JW, 1, NS, T, LDT, WORK, WORK( JW+1 ), &
                         LWORK-JW, INFO )
            CALL DGEMM( 'N', 'N', JW, NS, NS, ONE, V, LDV, T, LDT, ZERO, &
                        WV, LDWV )
            CALL DLACPY( 'A', JW, NS, WV, LDWV, V, LDV )
         END IF
!
!        ==== Update vertical slab in H ====
!
         IF( WANTT ) THEN
            LTOP = 1
         ELSE
            LTOP = KTOP
         END IF
         DO 70 KROW = LTOP, KWTOP - 1, NV
            KLN = MIN( NV, KWTOP-KROW )
            CALL DGEMM( 'N', 'N', KLN, JW, JW, ONE, H( KROW, KWTOP ), &
                        LDH, V, LDV, ZERO, WV, LDWV )
            CALL DLACPY( 'A', KLN, JW, WV, LDWV, H( KROW, KWTOP ), LDH )
   70    CONTINUE
!
!        ==== Update horizontal slab in H ====
!
         IF( WANTT ) THEN
            DO 80 KCOL = KBOT + 1, N, NH
               KLN = MIN( NH, N-KCOL+1 )
               CALL DGEMM( 'C', 'N', JW, KLN, JW, ONE, V, LDV, &
                           H( KWTOP, KCOL ), LDH, ZERO, T, LDT )
               CALL DLACPY( 'A', JW, KLN, T, LDT, H( KWTOP, KCOL ), &
                            LDH )
   80       CONTINUE
         END IF
!
!        ==== Update vertical slab in Z ====
!
         IF( WANTZ ) THEN
            DO 90 KROW = ILOZ, IHIZ, NV
               KLN = MIN( NV, IHIZ-KROW+1 )
               CALL DGEMM( 'N', 'N', KLN, JW, JW, ONE, Z( KROW, KWTOP ), &
                           LDZ, V, LDV, ZERO, WV, LDWV )
               CALL DLACPY( 'A', KLN, JW, WV, LDWV, Z( KROW, KWTOP ), &
                            LDZ )
   90       CONTINUE
         END IF
      END IF
!
!     ==== Return the number of deflations ... ====
!
      ND = JW - NS
!
!     ==== ... and the number of shifts. (Subtracting
!     .    INFQR from the spike length takes care
!     .    of the case of a rare QR failure while
!     .    calculating eigenvalues of the deflation
!     .    window.)  ====
!
      NS = NS - INFQR
!
!      ==== Return optimal workspace. ====
!
      WORK( 1 ) = DBLE( LWKOPT )
!
!     ==== End of DLAQR3 ====
!
      END SUBROUTINE DLAQR3


      SUBROUTINE DLAQR4( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, &
                         ILOZ, IHIZ, Z, LDZ, WORK, LWORK, INFO )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N
      LOGICAL            WANTT, WANTZ
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), WI( * ), WORK( * ), WR( * ), &
                         Z( LDZ, * )
!     ..
!
!     This subroutine implements one level of recursion for DLAQR0.
!     It is a complete implementation of the small bulge multi-shift
!     QR algorithm.  It may be called by DLAQR0 and, for large enough
!     deflation window size, it may be called by DLAQR3.  This
!     subroutine is identical to DLAQR0 except that it calls DLAQR2
!     instead of DLAQR3.
!
!     Purpose
!     =======
!
!     DLAQR4 computes the eigenvalues of a Hessenberg matrix H
!     and, optionally, the matrices T and Z from the Schur decomposition
!     H = Z T Z**T, where T is an upper quasi-triangular matrix (the
!     Schur form), and Z is the orthogonal matrix of Schur vectors.
!
!     Optionally Z may be postmultiplied into an input orthogonal
!     matrix Q so that this routine can give the Schur factorization
!     of a matrix A which has been reduced to the Hessenberg form H
!     by the orthogonal matrix Q:  A = Q*H*Q**T = (QZ)*T*(QZ)**T.
!
!     Arguments
!     =========
!
!     WANTT   (input) LOGICAL
!          = .TRUE. : the full Schur form T is required;
!          = .FALSE.: only eigenvalues are required.
!
!     WANTZ   (input) LOGICAL
!          = .TRUE. : the matrix of Schur vectors Z is required;
!          = .FALSE.: Schur vectors are not required.
!
!     N     (input) INTEGER
!           The order of the matrix H.  N .GE. 0.
!
!     ILO   (input) INTEGER
!     IHI   (input) INTEGER
!           It is assumed that H is already upper triangular in rows
!           and columns 1:ILO-1 and IHI+1:N and, if ILO.GT.1,
!           H(ILO,ILO-1) is zero. ILO and IHI are normally set by a
!           previous call to DGEBAL, and then passed to DGEHRD when the
!           matrix output by DGEBAL is reduced to Hessenberg form.
!           Otherwise, ILO and IHI should be set to 1 and N,
!           respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N.
!           If N = 0, then ILO = 1 and IHI = 0.
!
!     H     (input/output) DOUBLE PRECISION array, dimension (LDH,N)
!           On entry, the upper Hessenberg matrix H.
!           On exit, if INFO = 0 and WANTT is .TRUE., then H contains
!           the upper quasi-triangular matrix T from the Schur
!           decomposition (the Schur form); 2-by-2 diagonal blocks
!           (corresponding to complex conjugate pairs of eigenvalues)
!           are returned in standard form, with H(i,i) = H(i+1,i+1)
!           and H(i+1,i)*H(i,i+1).LT.0. If INFO = 0 and WANTT is
!           .FALSE., then the contents of H are unspecified on exit.
!           (The output value of H when INFO.GT.0 is given under the
!           description of INFO below.)
!
!           This subroutine may explicitly set H(i,j) = 0 for i.GT.j and
!           j = 1, 2, ... ILO-1 or j = IHI+1, IHI+2, ... N.
!
!     LDH   (input) INTEGER
!           The leading dimension of the array H. LDH .GE. max(1,N).
!
!     WR    (output) DOUBLE PRECISION array, dimension (IHI)
!     WI    (output) DOUBLE PRECISION array, dimension (IHI)
!           The real and imaginary parts, respectively, of the computed
!           eigenvalues of H(ILO:IHI,ILO:IHI) are stored WR(ILO:IHI)
!           and WI(ILO:IHI). If two eigenvalues are computed as a
!           complex conjugate pair, they are stored in consecutive
!           elements of WR and WI, say the i-th and (i+1)th, with
!           WI(i) .GT. 0 and WI(i+1) .LT. 0. If WANTT is .TRUE., then
!           the eigenvalues are stored in the same order as on the
!           diagonal of the Schur form returned in H, with
!           WR(i) = H(i,i) and, if H(i:i+1,i:i+1) is a 2-by-2 diagonal
!           block, WI(i) = sqrt(-H(i+1,i)*H(i,i+1)) and
!           WI(i+1) = -WI(i).
!
!     ILOZ     (input) INTEGER
!     IHIZ     (input) INTEGER
!           Specify the rows of Z to which transformations must be
!           applied if WANTZ is .TRUE..
!           1 .LE. ILOZ .LE. ILO; IHI .LE. IHIZ .LE. N.
!
!     Z     (input/output) DOUBLE PRECISION array, dimension (LDZ,IHI)
!           If WANTZ is .FALSE., then Z is not referenced.
!           If WANTZ is .TRUE., then Z(ILO:IHI,ILOZ:IHIZ) is
!           replaced by Z(ILO:IHI,ILOZ:IHIZ)*U where U is the
!           orthogonal Schur factor of H(ILO:IHI,ILO:IHI).
!           (The output value of Z when INFO.GT.0 is given under
!           the description of INFO below.)
!
!     LDZ   (input) INTEGER
!           The leading dimension of the array Z.  if WANTZ is .TRUE.
!           then LDZ.GE.MAX(1,IHIZ).  Otherwize, LDZ.GE.1.
!
!     WORK  (workspace/output) DOUBLE PRECISION array, dimension LWORK
!           On exit, if LWORK = -1, WORK(1) returns an estimate of
!           the optimal value for LWORK.
!
!     LWORK (input) INTEGER
!           The dimension of the array WORK.  LWORK .GE. max(1,N)
!           is sufficient, but LWORK typically as large as 6*N may
!           be required for optimal performance.  A workspace query
!           to determine the optimal workspace size is recommended.
!
!           If LWORK = -1, then DLAQR4 does a workspace query.
!           In this case, DLAQR4 checks the input parameters and
!           estimates the optimal workspace size for the given
!           values of N, ILO and IHI.  The estimate is returned
!           in WORK(1).  No error message related to LWORK is
!           issued by XERBLA.  Neither H nor Z are accessed.
!
!
!     INFO  (output) INTEGER
!             =  0:  successful exit
!           .GT. 0:  if INFO = i, DLAQR4 failed to compute all of
!                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR
!                and WI contain those eigenvalues which have been
!                successfully computed.  (Failures are rare.)
!
!                If INFO .GT. 0 and WANT is .FALSE., then on exit,
!                the remaining unconverged eigenvalues are the eigen-
!                values of the upper Hessenberg matrix rows and
!                columns ILO through INFO of the final, output
!                value of H.
!
!                If INFO .GT. 0 and WANTT is .TRUE., then on exit
!
!           (*)  (initial value of H)*U  = U*(final value of H)
!
!                where U is an orthogonal matrix.  The final
!                value of H is upper Hessenberg and quasi-triangular
!                in rows and columns INFO+1 through IHI.
!
!                If INFO .GT. 0 and WANTZ is .TRUE., then on exit
!
!                  (final value of Z(ILO:IHI,ILOZ:IHIZ)
!                   =  (initial value of Z(ILO:IHI,ILOZ:IHIZ)*U
!
!                where U is the orthogonal matrix in (*) (regard-
!                less of the value of WANTT.)
!
!                If INFO .GT. 0 and WANTZ is .FALSE., then Z is not
!                accessed.
!
!     ================================================================
!     Based on contributions by
!        Karen Braman and Ralph Byers, Department of Mathematics,
!        University of Kansas, USA
!
!     ================================================================
!     References:
!       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
!       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3
!       Performance, SIAM Journal of Matrix Analysis, volume 23, pages
!       929--947, 2002.
!
!       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
!       Algorithm Part II: Aggressive Early Deflation, SIAM Journal
!       of Matrix Analysis, volume 23, pages 948--973, 2002.
!
!     ================================================================
!     .. Parameters ..
!
!     ==== Matrices of order NTINY or smaller must be processed by
!     .    DLAHQR because of insufficient subdiagonal scratch space.
!     .    (This is a hard limit.) ====
!
!     ==== Exceptional deflation windows:  try to cure rare
!     .    slow convergence by increasing the size of the
!     .    deflation window after KEXNW iterations. =====
!
!     ==== Exceptional shifts: try to cure rare slow convergence
!     .    with ad-hoc exceptional shifts every KEXSH iterations.
!     .    The constants WILK1 and WILK2 are used to form the
!     .    exceptional shifts. ====
!
      INTEGER            NTINY
      PARAMETER          ( NTINY = 11 )
      INTEGER            KEXNW, KEXSH
      PARAMETER          ( KEXNW = 5, KEXSH = 6 )
      DOUBLE PRECISION   WILK1, WILK2
      PARAMETER          ( WILK1 = 0.75d0, WILK2 = -0.4375d0 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0d0, ONE = 1.0d0 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   AA, BB, CC, CS, DD, SN, SS, SWAP
      INTEGER            I, INF, IT, ITMAX, K, KACC22, KBOT, KDU, KS, &
                         KT, KTOP, KU, KV, KWH, KWTOP, KWV, LD, LS, &
                         LWKOPT, NDFL, NH, NHO, NIBBLE, NMIN, NS, NSMAX, &
                         NSR, NVE, NW, NWMAX, NWR
      LOGICAL            NWINC, SORTED
      CHARACTER          JBCMPZ*2
!     ..
!     .. External Functions ..
!      INTEGER            ILAENV
!      EXTERNAL           ILAENV
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   ZDUM( 1, 1 )
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DLACPY, DLAHQR, DLANV2, DLAQR2, DLAQR5
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, MAX, MIN, MOD
!     ..
!     .. Executable Statements ..
      INFO = 0
!
!     ==== Quick return for N = 0: nothing to do. ====
!
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = ONE
         RETURN
      END IF
!
!     ==== Set up job flags for ILAENV. ====
!
      IF( WANTT ) THEN
         JBCMPZ( 1: 1 ) = 'S'
      ELSE
         JBCMPZ( 1: 1 ) = 'E'
      END IF
      IF( WANTZ ) THEN
         JBCMPZ( 2: 2 ) = 'V'
      ELSE
         JBCMPZ( 2: 2 ) = 'N'
      END IF
!
!     ==== Tiny matrices must use DLAHQR. ====
!
      IF( N.LE.NTINY ) THEN
!
!        ==== Estimate optimal workspace. ====
!
         LWKOPT = 1
         IF( LWORK.NE.-1 ) &
            CALL DLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, &
                         ILOZ, IHIZ, Z, LDZ, INFO )
      ELSE
!
!        ==== Use small bulge multi-shift QR with aggressive early
!        .    deflation on larger-than-tiny matrices. ====
!
!        ==== Hope for the best. ====
!
         INFO = 0
!
!        ==== NWR = recommended deflation window size.  At this
!        .    point,  N .GT. NTINY = 11, so there is enough
!        .    subdiagonal workspace for NWR.GE.2 as required.
!        .    (In fact, there is enough subdiagonal space for
!        .    NWR.GE.3.) ====
!
         NWR = ILAENV( 13, 'DLAQR4', JBCMPZ, N, ILO, IHI, LWORK )
         NWR = MAX( 2, NWR )
         NWR = MIN( IHI-ILO+1, ( N-1 ) / 3, NWR )
         NW = NWR
!
!        ==== NSR = recommended number of simultaneous shifts.
!        .    At this point N .GT. NTINY = 11, so there is at
!        .    enough subdiagonal workspace for NSR to be even
!        .    and greater than or equal to two as required. ====
!
         NSR = ILAENV( 15, 'DLAQR4', JBCMPZ, N, ILO, IHI, LWORK )
         NSR = MIN( NSR, ( N+6 ) / 9, IHI-ILO )
         NSR = MAX( 2, NSR-MOD( NSR, 2 ) )
!
!        ==== Estimate optimal workspace ====
!
!        ==== Workspace query call to DLAQR2 ====
!
         CALL DLAQR2( WANTT, WANTZ, N, ILO, IHI, NWR+1, H, LDH, ILOZ, &
                      IHIZ, Z, LDZ, LS, LD, WR, WI, H, LDH, N, H, LDH, &
                      N, H, LDH, WORK, -1 )
!
!        ==== Optimal workspace = MAX(DLAQR5, DLAQR2) ====
!
         LWKOPT = MAX( 3*NSR / 2, INT( WORK( 1 ) ) )
!
!        ==== Quick return in case of workspace query. ====
!
         IF( LWORK.EQ.-1 ) THEN
            WORK( 1 ) = DBLE( LWKOPT )
            RETURN
         END IF
!
!        ==== DLAHQR/DLAQR0 crossover point ====
!
         NMIN = ILAENV( 12, 'DLAQR4', JBCMPZ, N, ILO, IHI, LWORK )
         NMIN = MAX( NTINY, NMIN )
!
!        ==== Nibble crossover point ====
!
         NIBBLE = ILAENV( 14, 'DLAQR4', JBCMPZ, N, ILO, IHI, LWORK )
         NIBBLE = MAX( 0, NIBBLE )
!
!        ==== Accumulate reflections during ttswp?  Use block
!        .    2-by-2 structure during matrix-matrix multiply? ====
!
         KACC22 = ILAENV( 16, 'DLAQR4', JBCMPZ, N, ILO, IHI, LWORK )
         KACC22 = MAX( 0, KACC22 )
         KACC22 = MIN( 2, KACC22 )
!
!        ==== NWMAX = the largest possible deflation window for
!        .    which there is sufficient workspace. ====
!
         NWMAX = MIN( ( N-1 ) / 3, LWORK / 2 )
!
!        ==== NSMAX = the Largest number of simultaneous shifts
!        .    for which there is sufficient workspace. ====
!
         NSMAX = MIN( ( N+6 ) / 9, 2*LWORK / 3 )
         NSMAX = NSMAX - MOD( NSMAX, 2 )
!
!        ==== NDFL: an iteration count restarted at deflation. ====
!
         NDFL = 1
!
!        ==== ITMAX = iteration limit ====
!
         ITMAX = MAX( 30, 2*KEXSH )*MAX( 10, ( IHI-ILO+1 ) )
!
!        ==== Last row and column in the active block ====
!
         KBOT = IHI
!
!        ==== Main Loop ====
!
         DO 80 IT = 1, ITMAX
!
!           ==== Done when KBOT falls below ILO ====
!
            IF( KBOT.LT.ILO ) &
               GO TO 90
!
!           ==== Locate active block ====
!
            DO 10 K = KBOT, ILO + 1, -1
               IF( H( K, K-1 ).EQ.ZERO ) &
                  GO TO 20
   10       CONTINUE
            K = ILO
   20       CONTINUE
            KTOP = K
!
!           ==== Select deflation window size ====
!
            NH = KBOT - KTOP + 1
            IF( NDFL.LT.KEXNW .OR. NH.LT.NW ) THEN
!
!              ==== Typical deflation window.  If possible and
!              .    advisable, nibble the entire active block.
!              .    If not, use size NWR or NWR+1 depending upon
!              .    which has the smaller corresponding subdiagonal
!              .    entry (a heuristic). ====
!
               NWINC = .TRUE.
               IF( NH.LE.MIN( NMIN, NWMAX ) ) THEN
                  NW = NH
               ELSE
                  NW = MIN( NWR, NH, NWMAX )
                  IF( NW.LT.NWMAX ) THEN
                     IF( NW.GE.NH-1 ) THEN
                        NW = NH
                     ELSE
                        KWTOP = KBOT - NW + 1
                        IF( ABS( H( KWTOP, KWTOP-1 ) ).GT. &
                            ABS( H( KWTOP-1, KWTOP-2 ) ) )NW = NW + 1
                     END IF
                  END IF
               END IF
            ELSE
!
!              ==== Exceptional deflation window.  If there have
!              .    been no deflations in KEXNW or more iterations,
!              .    then vary the deflation window size.   At first,
!              .    because, larger windows are, in general, more
!              .    powerful than smaller ones, rapidly increase the
!              .    window up to the maximum reasonable and possible.
!              .    Then maybe try a slightly smaller window.  ====
!
               IF( NWINC .AND. NW.LT.MIN( NWMAX, NH ) ) THEN
                  NW = MIN( NWMAX, NH, 2*NW )
               ELSE
                  NWINC = .FALSE.
                  IF( NW.EQ.NH .AND. NH.GT.2 ) &
                     NW = NH - 1
               END IF
            END IF
!
!           ==== Aggressive early deflation:
!           .    split workspace under the subdiagonal into
!           .      - an nw-by-nw work array V in the lower
!           .        left-hand-corner,
!           .      - an NW-by-at-least-NW-but-more-is-better
!           .        (NW-by-NHO) horizontal work array along
!           .        the bottom edge,
!           .      - an at-least-NW-but-more-is-better (NHV-by-NW)
!           .        vertical work array along the left-hand-edge.
!           .        ====
!
            KV = N - NW + 1
            KT = NW + 1
            NHO = ( N-NW-1 ) - KT + 1
            KWV = NW + 2
            NVE = ( N-NW ) - KWV + 1
!
!           ==== Aggressive early deflation ====
!
            CALL DLAQR2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ, &
                         IHIZ, Z, LDZ, LS, LD, WR, WI, H( KV, 1 ), LDH, &
                         NHO, H( KV, KT ), LDH, NVE, H( KWV, 1 ), LDH, &
                         WORK, LWORK )
!
!           ==== Adjust KBOT accounting for new deflations. ====
!
            KBOT = KBOT - LD
!
!           ==== KS points to the shifts. ====
!
            KS = KBOT - LS + 1
!
!           ==== Skip an expensive QR sweep if there is a (partly
!           .    heuristic) reason to expect that many eigenvalues
!           .    will deflate without it.  Here, the QR sweep is
!           .    skipped if many eigenvalues have just been deflated
!           .    or if the remaining active block is small.
!
            IF( ( LD.EQ.0 ) .OR. ( ( 100*LD.LE.NW*NIBBLE ) .AND. ( KBOT- &
                KTOP+1.GT.MIN( NMIN, NWMAX ) ) ) ) THEN
!
!              ==== NS = nominal number of simultaneous shifts.
!              .    This may be lowered (slightly) if DLAQR2
!              .    did not provide that many shifts. ====
!
               NS = MIN( NSMAX, NSR, MAX( 2, KBOT-KTOP ) )
               NS = NS - MOD( NS, 2 )
!
!              ==== If there have been no deflations
!              .    in a multiple of KEXSH iterations,
!              .    then try exceptional shifts.
!              .    Otherwise use shifts provided by
!              .    DLAQR2 above or from the eigenvalues
!              .    of a trailing principal submatrix. ====
!
               IF( MOD( NDFL, KEXSH ).EQ.0 ) THEN
                  KS = KBOT - NS + 1
                  DO 30 I = KBOT, MAX( KS+1, KTOP+2 ), -2
                     SS = ABS( H( I, I-1 ) ) + ABS( H( I-1, I-2 ) )
                     AA = WILK1*SS + H( I, I )
                     BB = SS
                     CC = WILK2*SS
                     DD = AA
                     CALL DLANV2( AA, BB, CC, DD, WR( I-1 ), WI( I-1 ), &
                                  WR( I ), WI( I ), CS, SN )
   30             CONTINUE
                  IF( KS.EQ.KTOP ) THEN
                     WR( KS+1 ) = H( KS+1, KS+1 )
                     WI( KS+1 ) = ZERO
                     WR( KS ) = WR( KS+1 )
                     WI( KS ) = WI( KS+1 )
                  END IF
               ELSE
!
!                 ==== Got NS/2 or fewer shifts? Use DLAHQR
!                 .    on a trailing principal submatrix to
!                 .    get more. (Since NS.LE.NSMAX.LE.(N+6)/9,
!                 .    there is enough space below the subdiagonal
!                 .    to fit an NS-by-NS scratch array.) ====
!
                  IF( KBOT-KS+1.LE.NS / 2 ) THEN
                     KS = KBOT - NS + 1
                     KT = N - NS + 1
                     CALL DLACPY( 'A', NS, NS, H( KS, KS ), LDH, &
                                  H( KT, 1 ), LDH )
                     CALL DLAHQR( .false., .false., NS, 1, NS, &
                                  H( KT, 1 ), LDH, WR( KS ), WI( KS ), &
                                  1, 1, ZDUM, 1, INF )
                     KS = KS + INF
!
!                    ==== In case of a rare QR failure use
!                    .    eigenvalues of the trailing 2-by-2
!                    .    principal submatrix.  ====
!
                     IF( KS.GE.KBOT ) THEN
                        AA = H( KBOT-1, KBOT-1 )
                        CC = H( KBOT, KBOT-1 )
                        BB = H( KBOT-1, KBOT )
                        DD = H( KBOT, KBOT )
                        CALL DLANV2( AA, BB, CC, DD, WR( KBOT-1 ), &
                                     WI( KBOT-1 ), WR( KBOT ), &
                                     WI( KBOT ), CS, SN )
                        KS = KBOT - 1
                     END IF
                  END IF
!
                  IF( KBOT-KS+1.GT.NS ) THEN
!
!                    ==== Sort the shifts (Helps a little)
!                    .    Bubble sort keeps complex conjugate
!                    .    pairs together. ====
!
                     SORTED = .false.
                     DO 50 K = KBOT, KS + 1, -1
                        IF( SORTED ) &
                           GO TO 60
                        SORTED = .true.
                        DO 40 I = KS, K - 1
                           IF( ABS( WR( I ) )+ABS( WI( I ) ).LT. &
                               ABS( WR( I+1 ) )+ABS( WI( I+1 ) ) ) THEN
                              SORTED = .false.
!
                              SWAP = WR( I )
                              WR( I ) = WR( I+1 )
                              WR( I+1 ) = SWAP
!
                              SWAP = WI( I )
                              WI( I ) = WI( I+1 )
                              WI( I+1 ) = SWAP
                           END IF
   40                   CONTINUE
   50                CONTINUE
   60                CONTINUE
                  END IF
!
!                 ==== Shuffle shifts into pairs of real shifts
!                 .    and pairs of complex conjugate shifts
!                 .    assuming complex conjugate shifts are
!                 .    already adjacent to one another. (Yes,
!                 .    they are.)  ====
!
                  DO 70 I = KBOT, KS + 2, -2
                     IF( WI( I ).NE.-WI( I-1 ) ) THEN
!
                        SWAP = WR( I )
                        WR( I ) = WR( I-1 )
                        WR( I-1 ) = WR( I-2 )
                        WR( I-2 ) = SWAP
!
                        SWAP = WI( I )
                        WI( I ) = WI( I-1 )
                        WI( I-1 ) = WI( I-2 )
                        WI( I-2 ) = SWAP
                     END IF
   70             CONTINUE
               END IF
!
!              ==== If there are only two shifts and both are
!              .    real, then use only one.  ====
!
               IF( KBOT-KS+1.EQ.2 ) THEN
                  IF( WI( KBOT ).EQ.ZERO ) THEN
                     IF( ABS( WR( KBOT )-H( KBOT, KBOT ) ).LT. &
                         ABS( WR( KBOT-1 )-H( KBOT, KBOT ) ) ) THEN
                        WR( KBOT-1 ) = WR( KBOT )
                     ELSE
                        WR( KBOT ) = WR( KBOT-1 )
                     END IF
                  END IF
               END IF
!
!              ==== Use up to NS of the the smallest magnatiude
!              .    shifts.  If there aren't NS shifts available,
!              .    then use them all, possibly dropping one to
!              .    make the number of shifts even. ====
!
               NS = MIN( NS, KBOT-KS+1 )
               NS = NS - MOD( NS, 2 )
               KS = KBOT - NS + 1
!
!              ==== Small-bulge multi-shift QR sweep:
!              .    split workspace under the subdiagonal into
!              .    - a KDU-by-KDU work array U in the lower
!              .      left-hand-corner,
!              .    - a KDU-by-at-least-KDU-but-more-is-better
!              .      (KDU-by-NHo) horizontal work array WH along
!              .      the bottom edge,
!              .    - and an at-least-KDU-but-more-is-better-by-KDU
!              .      (NVE-by-KDU) vertical work WV arrow along
!              .      the left-hand-edge. ====
!
               KDU = 3*NS - 3
               KU = N - KDU + 1
               KWH = KDU + 1
               NHO = ( N-KDU+1-4 ) - ( KDU+1 ) + 1
               KWV = KDU + 4
               NVE = N - KDU - KWV + 1
!
!              ==== Small-bulge multi-shift QR sweep ====
!
               CALL DLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NS, &
                            WR( KS ), WI( KS ), H, LDH, ILOZ, IHIZ, Z, &
                            LDZ, WORK, 3, H( KU, 1 ), LDH, NVE, &
                            H( KWV, 1 ), LDH, NHO, H( KU, KWH ), LDH )
            END IF
!
!           ==== Note progress (or the lack of it). ====
!
            IF( LD.GT.0 ) THEN
               NDFL = 1
            ELSE
               NDFL = NDFL + 1
            END IF
!
!           ==== End of main loop ====
   80    CONTINUE
!
!        ==== Iteration limit exceeded.  Set INFO to show where
!        .    the problem occurred and exit. ====
!
         INFO = KBOT
   90    CONTINUE
      END IF
!
!     ==== Return the optimal value of LWORK. ====
!
      WORK( 1 ) = DBLE( LWKOPT )
!
!     ==== End of DLAQR4 ====
!
      END SUBROUTINE DLAQR4


      SUBROUTINE DLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NSHFTS, &
                         SR, SI, H, LDH, ILOZ, IHIZ, Z, LDZ, V, LDV, U, &
                         LDU, NV, WV, LDWV, NH, WH, LDWH )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            IHIZ, ILOZ, KACC22, KBOT, KTOP, LDH, LDU, LDV, &
                         LDWH, LDWV, LDZ, N, NH, NSHFTS, NV
      LOGICAL            WANTT, WANTZ
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), SI( * ), SR( * ), U( LDU, * ), &
                         V( LDV, * ), WH( LDWH, * ), WV( LDWV, * ), &
                         Z( LDZ, * )
!     ..
!
!     This auxiliary subroutine called by DLAQR0 performs a
!     single small-bulge multi-shift QR sweep.
!
!      WANTT  (input) logical scalar
!             WANTT = .true. if the quasi-triangular Schur factor
!             is being computed.  WANTT is set to .false. otherwise.
!
!      WANTZ  (input) logical scalar
!             WANTZ = .true. if the orthogonal Schur factor is being
!             computed.  WANTZ is set to .false. otherwise.
!
!      KACC22 (input) integer with value 0, 1, or 2.
!             Specifies the computation mode of far-from-diagonal
!             orthogonal updates.
!        = 0: DLAQR5 does not accumulate reflections and does not
!             use matrix-matrix multiply to update far-from-diagonal
!             matrix entries.
!        = 1: DLAQR5 accumulates reflections and uses matrix-matrix
!             multiply to update the far-from-diagonal matrix entries.
!        = 2: DLAQR5 accumulates reflections, uses matrix-matrix
!             multiply to update the far-from-diagonal matrix entries,
!             and takes advantage of 2-by-2 block structure during
!             matrix multiplies.
!
!      N      (input) integer scalar
!             N is the order of the Hessenberg matrix H upon which this
!             subroutine operates.
!
!      KTOP   (input) integer scalar
!      KBOT   (input) integer scalar
!             These are the first and last rows and columns of an
!             isolated diagonal block upon which the QR sweep is to be
!             applied. It is assumed without a check that
!                       either KTOP = 1  or   H(KTOP,KTOP-1) = 0
!             and
!                       either KBOT = N  or   H(KBOT+1,KBOT) = 0.
!
!      NSHFTS (input) integer scalar
!             NSHFTS gives the number of simultaneous shifts.  NSHFTS
!             must be positive and even.
!
!      SR     (input) DOUBLE PRECISION array of size (NSHFTS)
!      SI     (input) DOUBLE PRECISION array of size (NSHFTS)
!             SR contains the real parts and SI contains the imaginary
!             parts of the NSHFTS shifts of origin that define the
!             multi-shift QR sweep.
!
!      H      (input/output) DOUBLE PRECISION array of size (LDH,N)
!             On input H contains a Hessenberg matrix.  On output a
!             multi-shift QR sweep with shifts SR(J)+i*SI(J) is applied
!             to the isolated diagonal block in rows and columns KTOP
!             through KBOT.
!
!      LDH    (input) integer scalar
!             LDH is the leading dimension of H just as declared in the
!             calling procedure.  LDH.GE.MAX(1,N).
!
!      ILOZ   (input) INTEGER
!      IHIZ   (input) INTEGER
!             Specify the rows of Z to which transformations must be
!             applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N
!
!      Z      (input/output) DOUBLE PRECISION array of size (LDZ,IHI)
!             If WANTZ = .TRUE., then the QR Sweep orthogonal
!             similarity transformation is accumulated into
!             Z(ILOZ:IHIZ,ILO:IHI) from the right.
!             If WANTZ = .FALSE., then Z is unreferenced.
!
!      LDZ    (input) integer scalar
!             LDA is the leading dimension of Z just as declared in
!             the calling procedure. LDZ.GE.N.
!
!      V      (workspace) DOUBLE PRECISION array of size (LDV,NSHFTS/2)
!
!      LDV    (input) integer scalar
!             LDV is the leading dimension of V as declared in the
!             calling procedure.  LDV.GE.3.
!
!      U      (workspace) DOUBLE PRECISION array of size
!             (LDU,3*NSHFTS-3)
!
!      LDU    (input) integer scalar
!             LDU is the leading dimension of U just as declared in the
!             in the calling subroutine.  LDU.GE.3*NSHFTS-3.
!
!      NH     (input) integer scalar
!             NH is the number of columns in array WH available for
!             workspace. NH.GE.1.
!
!      WH     (workspace) DOUBLE PRECISION array of size (LDWH,NH)
!
!      LDWH   (input) integer scalar
!             Leading dimension of WH just as declared in the
!             calling procedure.  LDWH.GE.3*NSHFTS-3.
!
!      NV     (input) integer scalar
!             NV is the number of rows in WV agailable for workspace.
!             NV.GE.1.
!
!      WV     (workspace) DOUBLE PRECISION array of size
!             (LDWV,3*NSHFTS-3)
!
!      LDWV   (input) integer scalar
!             LDWV is the leading dimension of WV as declared in the
!             in the calling subroutine.  LDWV.GE.NV.
!
!
!     ================================================================
!     Based on contributions by
!        Karen Braman and Ralph Byers, Department of Mathematics,
!        University of Kansas, USA
!
!     ============================================================
!     Reference:
!
!     K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
!     Algorithm Part I: Maintaining Well Focused Shifts, and
!     Level 3 Performance, SIAM Journal of Matrix Analysis,
!     volume 23, pages 929--947, 2002.
!
!     ============================================================
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0d0, ONE = 1.0d0 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   ALPHA, BETA, H11, H12, H21, H22, REFSUM, &
                         SAFMAX, SAFMIN, SCL, SMLNUM, SWAP, TST1, TST2, &
                         ULP
      INTEGER            I, I2, I4, INCOL, J, J2, J4, JBOT, JCOL, JLEN, &
                         JROW, JTOP, K, K1, KDU, KMS, KNZ, KRCOL, KZS, &
                         M, M22, MBOT, MEND, MSTART, MTOP, NBMPS, NDCOL, &
                         NS, NU
      LOGICAL            ACCUM, BLK22, BMP22
!     ..
!     .. External Functions ..
!      DOUBLE PRECISION   DLAMCH
!      EXTERNAL           DLAMCH
!     ..
!     .. Intrinsic Functions ..
!
      INTRINSIC          ABS, DBLE, MAX, MIN, MOD
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   VT( 3 )
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DGEMM, DLABAD, DLACPY, DLAQR1, DLARFG, DLASET, &
!                         DTRMM
!     ..
!     .. Executable Statements ..
!
!     ==== If there are no shifts, then there is nothing to do. ====
!
      IF( NSHFTS.LT.2 ) &
         RETURN
!
!     ==== If the active block is empty or 1-by-1, then there
!     .    is nothing to do. ====
!
      IF( KTOP.GE.KBOT ) &
         RETURN
!
!     ==== Shuffle shifts into pairs of real shifts and pairs
!     .    of complex conjugate shifts assuming complex
!     .    conjugate shifts are already adjacent to one
!     .    another. ====
!
      DO 10 I = 1, NSHFTS - 2, 2
         IF( SI( I ).NE.-SI( I+1 ) ) THEN
!
            SWAP = SR( I )
            SR( I ) = SR( I+1 )
            SR( I+1 ) = SR( I+2 )
            SR( I+2 ) = SWAP
!
            SWAP = SI( I )
            SI( I ) = SI( I+1 )
            SI( I+1 ) = SI( I+2 )
            SI( I+2 ) = SWAP
         END IF
   10 CONTINUE
!
!     ==== NSHFTS is supposed to be even, but if is odd,
!     .    then simply reduce it by one.  The shuffle above
!     .    ensures that the dropped shift is real and that
!     .    the remaining shifts are paired. ====
!
      NS = NSHFTS - MOD( NSHFTS, 2 )
!
!     ==== Machine constants for deflation ====
!
      SAFMIN = DLAMCH( 'SAFE MINIMUM' )
      SAFMAX = ONE / SAFMIN
      CALL DLABAD( SAFMIN, SAFMAX )
      ULP = DLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( DBLE( N ) / ULP )
!
!     ==== Use accumulated reflections to update far-from-diagonal
!     .    entries ? ====
!
      ACCUM = ( KACC22.EQ.1 ) .OR. ( KACC22.EQ.2 )
!
!     ==== If so, exploit the 2-by-2 block structure? ====
!
      BLK22 = ( NS.GT.2 ) .AND. ( KACC22.EQ.2 )
!
!     ==== clear trash ====
!
      IF( KTOP+2.LE.KBOT ) &
         H( KTOP+2, KTOP ) = ZERO
!
!     ==== NBMPS = number of 2-shift bulges in the chain ====
!
      NBMPS = NS / 2
!
!     ==== KDU = width of slab ====
!
      KDU = 6*NBMPS - 3
!
!     ==== Create and chase chains of NBMPS bulges ====
!
      DO 220 INCOL = 3*( 1-NBMPS ) + KTOP - 1, KBOT - 2, 3*NBMPS - 2
         NDCOL = INCOL + KDU
         IF( ACCUM ) &
            CALL DLASET( 'ALL', KDU, KDU, ZERO, ONE, U, LDU )
!
!        ==== Near-the-diagonal bulge chase.  The following loop
!        .    performs the near-the-diagonal part of a small bulge
!        .    multi-shift QR sweep.  Each 6*NBMPS-2 column diagonal
!        .    chunk extends from column INCOL to column NDCOL
!        .    (including both column INCOL and column NDCOL). The
!        .    following loop chases a 3*NBMPS column long chain of
!        .    NBMPS bulges 3*NBMPS-2 columns to the right.  (INCOL
!        .    may be less than KTOP and and NDCOL may be greater than
!        .    KBOT indicating phantom columns from which to chase
!        .    bulges before they are actually introduced or to which
!        .    to chase bulges beyond column KBOT.)  ====
!
         DO 150 KRCOL = INCOL, MIN( INCOL+3*NBMPS-3, KBOT-2 )
!
!           ==== Bulges number MTOP to MBOT are active double implicit
!           .    shift bulges.  There may or may not also be small
!           .    2-by-2 bulge, if there is room.  The inactive bulges
!           .    (if any) must wait until the active bulges have moved
!           .    down the diagonal to make room.  The phantom matrix
!           .    paradigm described above helps keep track.  ====
!
            MTOP = MAX( 1, ( ( KTOP-1 )-KRCOL+2 ) / 3+1 )
            MBOT = MIN( NBMPS, ( KBOT-KRCOL ) / 3 )
            M22 = MBOT + 1
            BMP22 = ( MBOT.LT.NBMPS ) .AND. ( KRCOL+3*( M22-1 ) ).EQ. &
                    ( KBOT-2 )
!
!           ==== Generate reflections to chase the chain right
!           .    one column.  (The minimum value of K is KTOP-1.) ====
!
            DO 20 M = MTOP, MBOT
               K = KRCOL + 3*( M-1 )
               IF( K.EQ.KTOP-1 ) THEN
                  CALL DLAQR1( 3, H( KTOP, KTOP ), LDH, SR( 2*M-1 ), &
                               SI( 2*M-1 ), SR( 2*M ), SI( 2*M ), &
                               V( 1, M ) )
                  ALPHA = V( 1, M )
                  CALL DLARFG( 3, ALPHA, V( 2, M ), 1, V( 1, M ) )
               ELSE
                  BETA = H( K+1, K )
                  V( 2, M ) = H( K+2, K )
                  V( 3, M ) = H( K+3, K )
                  CALL DLARFG( 3, BETA, V( 2, M ), 1, V( 1, M ) )
!
!                 ==== A Bulge may collapse because of vigilant
!                 .    deflation or destructive underflow.  (The
!                 .    initial bulge is always collapsed.) Use
!                 .    the two-small-subdiagonals trick to try
!                 .    to get it started again. If V(2,M).NE.0 and
!                 .    V(3,M) = H(K+3,K+1) = H(K+3,K+2) = 0, then
!                 .    this bulge is collapsing into a zero
!                 .    subdiagonal.  It will be restarted next
!                 .    trip through the loop.)
!
                  IF( V( 1, M ).NE.ZERO .AND. &
                      ( V( 3, M ).NE.ZERO .OR. ( H( K+3, &
                      K+1 ).EQ.ZERO .AND. H( K+3, K+2 ).EQ.ZERO ) ) ) &
                       THEN
!
!                    ==== Typical case: not collapsed (yet). ====
!
                     H( K+1, K ) = BETA
                     H( K+2, K ) = ZERO
                     H( K+3, K ) = ZERO
                  ELSE
!
!                    ==== Atypical case: collapsed.  Attempt to
!                    .    reintroduce ignoring H(K+1,K).  If the
!                    .    fill resulting from the new reflector
!                    .    is too large, then abandon it.
!                    .    Otherwise, use the new one. ====
!
                     CALL DLAQR1( 3, H( K+1, K+1 ), LDH, SR( 2*M-1 ), &
                                  SI( 2*M-1 ), SR( 2*M ), SI( 2*M ), &
                                  VT )
                     SCL = ABS( VT( 1 ) ) + ABS( VT( 2 ) ) + &
                           ABS( VT( 3 ) )
                     IF( SCL.NE.ZERO ) THEN
                        VT( 1 ) = VT( 1 ) / SCL
                        VT( 2 ) = VT( 2 ) / SCL
                        VT( 3 ) = VT( 3 ) / SCL
                     END IF
!
!                    ==== The following is the traditional and
!                    .    conservative two-small-subdiagonals
!                    .    test.  ====
!                    .
                     IF( ABS( H( K+1, K ) )*( ABS( VT( 2 ) )+ &
                         ABS( VT( 3 ) ) ).GT.ULP*ABS( VT( 1 ) )* &
                         ( ABS( H( K, K ) )+ABS( H( K+1, &
                         K+1 ) )+ABS( H( K+2, K+2 ) ) ) ) THEN
!
!                       ==== Starting a new bulge here would
!                       .    create non-negligible fill.   If
!                       .    the old reflector is diagonal (only
!                       .    possible with underflows), then
!                       .    change it to I.  Otherwise, use
!                       .    it with trepidation. ====
!
                        IF( V( 2, M ).EQ.ZERO .AND. V( 3, M ).EQ.ZERO ) &
                             THEN
                           V( 1, M ) = ZERO
                        ELSE
                           H( K+1, K ) = BETA
                           H( K+2, K ) = ZERO
                           H( K+3, K ) = ZERO
                        END IF
                     ELSE
!
!                       ==== Stating a new bulge here would
!                       .    create only negligible fill.
!                       .    Replace the old reflector with
!                       .    the new one. ====
!
                        ALPHA = VT( 1 )
                        CALL DLARFG( 3, ALPHA, VT( 2 ), 1, VT( 1 ) )
                        REFSUM = H( K+1, K ) + H( K+2, K )*VT( 2 ) + &
                                 H( K+3, K )*VT( 3 )
                        H( K+1, K ) = H( K+1, K ) - VT( 1 )*REFSUM
                        H( K+2, K ) = ZERO
                        H( K+3, K ) = ZERO
                        V( 1, M ) = VT( 1 )
                        V( 2, M ) = VT( 2 )
                        V( 3, M ) = VT( 3 )
                     END IF
                  END IF
               END IF
   20       CONTINUE
!
!           ==== Generate a 2-by-2 reflection, if needed. ====
!
            K = KRCOL + 3*( M22-1 )
            IF( BMP22 ) THEN
               IF( K.EQ.KTOP-1 ) THEN
                  CALL DLAQR1( 2, H( K+1, K+1 ), LDH, SR( 2*M22-1 ), &
                               SI( 2*M22-1 ), SR( 2*M22 ), SI( 2*M22 ), &
                               V( 1, M22 ) )
                  BETA = V( 1, M22 )
                  CALL DLARFG( 2, BETA, V( 2, M22 ), 1, V( 1, M22 ) )
               ELSE
                  BETA = H( K+1, K )
                  V( 2, M22 ) = H( K+2, K )
                  CALL DLARFG( 2, BETA, V( 2, M22 ), 1, V( 1, M22 ) )
                  H( K+1, K ) = BETA
                  H( K+2, K ) = ZERO
               END IF
            ELSE
!
!              ==== Initialize V(1,M22) here to avoid possible undefined
!              .    variable problems later. ====
!
               V( 1, M22 ) = ZERO
            END IF
!
!           ==== Multiply H by reflections from the left ====
!
            IF( ACCUM ) THEN
               JBOT = MIN( NDCOL, KBOT )
            ELSE IF( WANTT ) THEN
               JBOT = N
            ELSE
               JBOT = KBOT
            END IF
            DO 40 J = MAX( KTOP, KRCOL ), JBOT
               MEND = MIN( MBOT, ( J-KRCOL+2 ) / 3 )
               DO 30 M = MTOP, MEND
                  K = KRCOL + 3*( M-1 )
                  REFSUM = V( 1, M )*( H( K+1, J )+V( 2, M )* &
                           H( K+2, J )+V( 3, M )*H( K+3, J ) )
                  H( K+1, J ) = H( K+1, J ) - REFSUM
                  H( K+2, J ) = H( K+2, J ) - REFSUM*V( 2, M )
                  H( K+3, J ) = H( K+3, J ) - REFSUM*V( 3, M )
   30          CONTINUE
   40       CONTINUE
            IF( BMP22 ) THEN
               K = KRCOL + 3*( M22-1 )
               DO 50 J = MAX( K+1, KTOP ), JBOT
                  REFSUM = V( 1, M22 )*( H( K+1, J )+V( 2, M22 )* &
                           H( K+2, J ) )
                  H( K+1, J ) = H( K+1, J ) - REFSUM
                  H( K+2, J ) = H( K+2, J ) - REFSUM*V( 2, M22 )
   50          CONTINUE
            END IF
!
!           ==== Multiply H by reflections from the right.
!           .    Delay filling in the last row until the
!           .    vigilant deflation check is complete. ====
!
            IF( ACCUM ) THEN
               JTOP = MAX( KTOP, INCOL )
            ELSE IF( WANTT ) THEN
               JTOP = 1
            ELSE
               JTOP = KTOP
            END IF
            DO 90 M = MTOP, MBOT
               IF( V( 1, M ).NE.ZERO ) THEN
                  K = KRCOL + 3*( M-1 )
                  DO 60 J = JTOP, MIN( KBOT, K+3 )
                     REFSUM = V( 1, M )*( H( J, K+1 )+V( 2, M )* &
                              H( J, K+2 )+V( 3, M )*H( J, K+3 ) )
                     H( J, K+1 ) = H( J, K+1 ) - REFSUM
                     H( J, K+2 ) = H( J, K+2 ) - REFSUM*V( 2, M )
                     H( J, K+3 ) = H( J, K+3 ) - REFSUM*V( 3, M )
   60             CONTINUE
!
                  IF( ACCUM ) THEN
!
!                    ==== Accumulate U. (If necessary, update Z later
!                    .    with with an efficient matrix-matrix
!                    .    multiply.) ====
!
                     KMS = K - INCOL
                     DO 70 J = MAX( 1, KTOP-INCOL ), KDU
                        REFSUM = V( 1, M )*( U( J, KMS+1 )+V( 2, M )* &
                                 U( J, KMS+2 )+V( 3, M )*U( J, KMS+3 ) )
                        U( J, KMS+1 ) = U( J, KMS+1 ) - REFSUM
                        U( J, KMS+2 ) = U( J, KMS+2 ) - REFSUM*V( 2, M )
                        U( J, KMS+3 ) = U( J, KMS+3 ) - REFSUM*V( 3, M )
   70                CONTINUE
                  ELSE IF( WANTZ ) THEN
!
!                    ==== U is not accumulated, so update Z
!                    .    now by multiplying by reflections
!                    .    from the right. ====
!
                     DO 80 J = ILOZ, IHIZ
                        REFSUM = V( 1, M )*( Z( J, K+1 )+V( 2, M )* &
                                 Z( J, K+2 )+V( 3, M )*Z( J, K+3 ) )
                        Z( J, K+1 ) = Z( J, K+1 ) - REFSUM
                        Z( J, K+2 ) = Z( J, K+2 ) - REFSUM*V( 2, M )
                        Z( J, K+3 ) = Z( J, K+3 ) - REFSUM*V( 3, M )
   80                CONTINUE
                  END IF
               END IF
   90       CONTINUE
!
!           ==== Special case: 2-by-2 reflection (if needed) ====
!
            K = KRCOL + 3*( M22-1 )
            IF( BMP22 .AND. ( V( 1, M22 ).NE.ZERO ) ) THEN
               DO 100 J = JTOP, MIN( KBOT, K+3 )
                  REFSUM = V( 1, M22 )*( H( J, K+1 )+V( 2, M22 )* &
                           H( J, K+2 ) )
                  H( J, K+1 ) = H( J, K+1 ) - REFSUM
                  H( J, K+2 ) = H( J, K+2 ) - REFSUM*V( 2, M22 )
  100          CONTINUE
!
               IF( ACCUM ) THEN
                  KMS = K - INCOL
                  DO 110 J = MAX( 1, KTOP-INCOL ), KDU
                     REFSUM = V( 1, M22 )*( U( J, KMS+1 )+V( 2, M22 )* &
                              U( J, KMS+2 ) )
                     U( J, KMS+1 ) = U( J, KMS+1 ) - REFSUM
                     U( J, KMS+2 ) = U( J, KMS+2 ) - REFSUM*V( 2, M22 )
  110             CONTINUE
               ELSE IF( WANTZ ) THEN
                  DO 120 J = ILOZ, IHIZ
                     REFSUM = V( 1, M22 )*( Z( J, K+1 )+V( 2, M22 )* &
                              Z( J, K+2 ) )
                     Z( J, K+1 ) = Z( J, K+1 ) - REFSUM
                     Z( J, K+2 ) = Z( J, K+2 ) - REFSUM*V( 2, M22 )
  120             CONTINUE
               END IF
            END IF
!
!           ==== Vigilant deflation check ====
!
            MSTART = MTOP
            IF( KRCOL+3*( MSTART-1 ).LT.KTOP ) &
               MSTART = MSTART + 1
            MEND = MBOT
            IF( BMP22 ) &
               MEND = MEND + 1
            IF( KRCOL.EQ.KBOT-2 ) &
               MEND = MEND + 1
            DO 130 M = MSTART, MEND
               K = MIN( KBOT-1, KRCOL+3*( M-1 ) )
!
!              ==== The following convergence test requires that
!              .    the tradition small-compared-to-nearby-diagonals
!              .    criterion and the Ahues & Tisseur (LAWN 122, 1997)
!              .    criteria both be satisfied.  The latter improves
!              .    accuracy in some examples. Falling back on an
!              .    alternate convergence criterion when TST1 or TST2
!              .    is zero (as done here) is traditional but probably
!              .    unnecessary. ====
!
               IF( H( K+1, K ).NE.ZERO ) THEN
                  TST1 = ABS( H( K, K ) ) + ABS( H( K+1, K+1 ) )
                  IF( TST1.EQ.ZERO ) THEN
                     IF( K.GE.KTOP+1 ) &
                        TST1 = TST1 + ABS( H( K, K-1 ) )
                     IF( K.GE.KTOP+2 ) &
                        TST1 = TST1 + ABS( H( K, K-2 ) )
                     IF( K.GE.KTOP+3 ) &
                        TST1 = TST1 + ABS( H( K, K-3 ) )
                     IF( K.LE.KBOT-2 ) &
                        TST1 = TST1 + ABS( H( K+2, K+1 ) )
                     IF( K.LE.KBOT-3 ) &
                        TST1 = TST1 + ABS( H( K+3, K+1 ) )
                     IF( K.LE.KBOT-4 ) &
                        TST1 = TST1 + ABS( H( K+4, K+1 ) )
                  END IF
                  IF( ABS( H( K+1, K ) ).LE.MAX( SMLNUM, ULP*TST1 ) ) &
                       THEN
                     H12 = MAX( ABS( H( K+1, K ) ), ABS( H( K, K+1 ) ) )
                     H21 = MIN( ABS( H( K+1, K ) ), ABS( H( K, K+1 ) ) )
                     H11 = MAX( ABS( H( K+1, K+1 ) ), &
                           ABS( H( K, K )-H( K+1, K+1 ) ) )
                     H22 = MIN( ABS( H( K+1, K+1 ) ), &
                           ABS( H( K, K )-H( K+1, K+1 ) ) )
                     SCL = H11 + H12
                     TST2 = H22*( H11 / SCL )
!
                     IF( TST2.EQ.ZERO .OR. H21*( H12 / SCL ).LE. &
                         MAX( SMLNUM, ULP*TST2 ) )H( K+1, K ) = ZERO
                  END IF
               END IF
  130       CONTINUE
!
!           ==== Fill in the last row of each bulge. ====
!
            MEND = MIN( NBMPS, ( KBOT-KRCOL-1 ) / 3 )
            DO 140 M = MTOP, MEND
               K = KRCOL + 3*( M-1 )
               REFSUM = V( 1, M )*V( 3, M )*H( K+4, K+3 )
               H( K+4, K+1 ) = -REFSUM
               H( K+4, K+2 ) = -REFSUM*V( 2, M )
               H( K+4, K+3 ) = H( K+4, K+3 ) - REFSUM*V( 3, M )
  140       CONTINUE
!
!           ==== End of near-the-diagonal bulge chase. ====
!
  150    CONTINUE
!
!        ==== Use U (if accumulated) to update far-from-diagonal
!        .    entries in H.  If required, use U to update Z as
!        .    well. ====
!
         IF( ACCUM ) THEN
            IF( WANTT ) THEN
               JTOP = 1
               JBOT = N
            ELSE
               JTOP = KTOP
               JBOT = KBOT
            END IF
            IF( ( .NOT.BLK22 ) .OR. ( INCOL.LT.KTOP ) .OR. &
                ( NDCOL.GT.KBOT ) .OR. ( NS.LE.2 ) ) THEN
!
!              ==== Updates not exploiting the 2-by-2 block
!              .    structure of U.  K1 and NU keep track of
!              .    the location and size of U in the special
!              .    cases of introducing bulges and chasing
!              .    bulges off the bottom.  In these special
!              .    cases and in case the number of shifts
!              .    is NS = 2, there is no 2-by-2 block
!              .    structure to exploit.  ====
!
               K1 = MAX( 1, KTOP-INCOL )
               NU = ( KDU-MAX( 0, NDCOL-KBOT ) ) - K1 + 1
!
!              ==== Horizontal Multiply ====
!
               DO 160 JCOL = MIN( NDCOL, KBOT ) + 1, JBOT, NH
                  JLEN = MIN( NH, JBOT-JCOL+1 )
                  CALL DGEMM( 'C', 'N', NU, JLEN, NU, ONE, U( K1, K1 ), &
                              LDU, H( INCOL+K1, JCOL ), LDH, ZERO, WH, &
                              LDWH )
                  CALL DLACPY( 'ALL', NU, JLEN, WH, LDWH, &
                               H( INCOL+K1, JCOL ), LDH )
  160          CONTINUE
!
!              ==== Vertical multiply ====
!
               DO 170 JROW = JTOP, MAX( KTOP, INCOL ) - 1, NV
                  JLEN = MIN( NV, MAX( KTOP, INCOL )-JROW )
                  CALL DGEMM( 'N', 'N', JLEN, NU, NU, ONE, &
                              H( JROW, INCOL+K1 ), LDH, U( K1, K1 ), &
                              LDU, ZERO, WV, LDWV )
                  CALL DLACPY( 'ALL', JLEN, NU, WV, LDWV, &
                               H( JROW, INCOL+K1 ), LDH )
  170          CONTINUE
!
!              ==== Z multiply (also vertical) ====
!
               IF( WANTZ ) THEN
                  DO 180 JROW = ILOZ, IHIZ, NV
                     JLEN = MIN( NV, IHIZ-JROW+1 )
                     CALL DGEMM( 'N', 'N', JLEN, NU, NU, ONE, &
                                 Z( JROW, INCOL+K1 ), LDZ, U( K1, K1 ), &
                                 LDU, ZERO, WV, LDWV )
                     CALL DLACPY( 'ALL', JLEN, NU, WV, LDWV, &
                                  Z( JROW, INCOL+K1 ), LDZ )
  180             CONTINUE
               END IF
            ELSE
!
!              ==== Updates exploiting U's 2-by-2 block structure.
!              .    (I2, I4, J2, J4 are the last rows and columns
!              .    of the blocks.) ====
!
               I2 = ( KDU+1 ) / 2
               I4 = KDU
               J2 = I4 - I2
               J4 = KDU
!
!              ==== KZS and KNZ deal with the band of zeros
!              .    along the diagonal of one of the triangular
!              .    blocks. ====
!
               KZS = ( J4-J2 ) - ( NS+1 )
               KNZ = NS + 1
!
!              ==== Horizontal multiply ====
!
               DO 190 JCOL = MIN( NDCOL, KBOT ) + 1, JBOT, NH
                  JLEN = MIN( NH, JBOT-JCOL+1 )
!
!                 ==== Copy bottom of H to top+KZS of scratch ====
!                  (The first KZS rows get multiplied by zero.) ====
!
                  CALL DLACPY( 'ALL', KNZ, JLEN, H( INCOL+1+J2, JCOL ), &
                               LDH, WH( KZS+1, 1 ), LDWH )
!
!                 ==== Multiply by U21' ====
!
                  CALL DLASET( 'ALL', KZS, JLEN, ZERO, ZERO, WH, LDWH )
                  CALL DTRMM( 'L', 'U', 'C', 'N', KNZ, JLEN, ONE, &
                              U( J2+1, 1+KZS ), LDU, WH( KZS+1, 1 ), &
                              LDWH )
!
!                 ==== Multiply top of H by U11' ====
!
                  CALL DGEMM( 'C', 'N', I2, JLEN, J2, ONE, U, LDU, &
                              H( INCOL+1, JCOL ), LDH, ONE, WH, LDWH )
!
!                 ==== Copy top of H bottom of WH ====
!
                  CALL DLACPY( 'ALL', J2, JLEN, H( INCOL+1, JCOL ), LDH, &
                               WH( I2+1, 1 ), LDWH )
!
!                 ==== Multiply by U21' ====
!
                  CALL DTRMM( 'L', 'L', 'C', 'N', J2, JLEN, ONE, &
                              U( 1, I2+1 ), LDU, WH( I2+1, 1 ), LDWH )
!
!                 ==== Multiply by U22 ====
!
                  CALL DGEMM( 'C', 'N', I4-I2, JLEN, J4-J2, ONE, &
                              U( J2+1, I2+1 ), LDU, &
                              H( INCOL+1+J2, JCOL ), LDH, ONE, &
                              WH( I2+1, 1 ), LDWH )
!
!                 ==== Copy it back ====
!
                  CALL DLACPY( 'ALL', KDU, JLEN, WH, LDWH, &
                               H( INCOL+1, JCOL ), LDH )
  190          CONTINUE
!
!              ==== Vertical multiply ====
!
               DO 200 JROW = JTOP, MAX( INCOL, KTOP ) - 1, NV
                  JLEN = MIN( NV, MAX( INCOL, KTOP )-JROW )
!
!                 ==== Copy right of H to scratch (the first KZS
!                 .    columns get multiplied by zero) ====
!
                  CALL DLACPY( 'ALL', JLEN, KNZ, H( JROW, INCOL+1+J2 ), &
                               LDH, WV( 1, 1+KZS ), LDWV )
!
!                 ==== Multiply by U21 ====
!
                  CALL DLASET( 'ALL', JLEN, KZS, ZERO, ZERO, WV, LDWV )
                  CALL DTRMM( 'R', 'U', 'N', 'N', JLEN, KNZ, ONE, &
                              U( J2+1, 1+KZS ), LDU, WV( 1, 1+KZS ), &
                              LDWV )
!
!                 ==== Multiply by U11 ====
!
                  CALL DGEMM( 'N', 'N', JLEN, I2, J2, ONE, &
                              H( JROW, INCOL+1 ), LDH, U, LDU, ONE, WV, &
                              LDWV )
!
!                 ==== Copy left of H to right of scratch ====
!
                  CALL DLACPY( 'ALL', JLEN, J2, H( JROW, INCOL+1 ), LDH, &
                               WV( 1, 1+I2 ), LDWV )
!
!                 ==== Multiply by U21 ====
!
                  CALL DTRMM( 'R', 'L', 'N', 'N', JLEN, I4-I2, ONE, &
                              U( 1, I2+1 ), LDU, WV( 1, 1+I2 ), LDWV )
!
!                 ==== Multiply by U22 ====
!
                  CALL DGEMM( 'N', 'N', JLEN, I4-I2, J4-J2, ONE, &
                              H( JROW, INCOL+1+J2 ), LDH, &
                              U( J2+1, I2+1 ), LDU, ONE, WV( 1, 1+I2 ), &
                              LDWV )
!
!                 ==== Copy it back ====
!
                  CALL DLACPY( 'ALL', JLEN, KDU, WV, LDWV, &
                               H( JROW, INCOL+1 ), LDH )
  200          CONTINUE
!
!              ==== Multiply Z (also vertical) ====
!
               IF( WANTZ ) THEN
                  DO 210 JROW = ILOZ, IHIZ, NV
                     JLEN = MIN( NV, IHIZ-JROW+1 )
!
!                    ==== Copy right of Z to left of scratch (first
!                    .     KZS columns get multiplied by zero) ====
!
                     CALL DLACPY( 'ALL', JLEN, KNZ, &
                                  Z( JROW, INCOL+1+J2 ), LDZ, &
                                  WV( 1, 1+KZS ), LDWV )
!
!                    ==== Multiply by U12 ====
!
                     CALL DLASET( 'ALL', JLEN, KZS, ZERO, ZERO, WV, &
                                  LDWV )
                     CALL DTRMM( 'R', 'U', 'N', 'N', JLEN, KNZ, ONE, &
                                 U( J2+1, 1+KZS ), LDU, WV( 1, 1+KZS ), &
                                 LDWV )
!
!                    ==== Multiply by U11 ====
!
                     CALL DGEMM( 'N', 'N', JLEN, I2, J2, ONE, &
                                 Z( JROW, INCOL+1 ), LDZ, U, LDU, ONE, &
                                 WV, LDWV )
!
!                    ==== Copy left of Z to right of scratch ====
!
                     CALL DLACPY( 'ALL', JLEN, J2, Z( JROW, INCOL+1 ), &
                                  LDZ, WV( 1, 1+I2 ), LDWV )
!
!                    ==== Multiply by U21 ====
!
                     CALL DTRMM( 'R', 'L', 'N', 'N', JLEN, I4-I2, ONE, &
                                 U( 1, I2+1 ), LDU, WV( 1, 1+I2 ), &
                                 LDWV )
!
!                    ==== Multiply by U22 ====
!
                     CALL DGEMM( 'N', 'N', JLEN, I4-I2, J4-J2, ONE, &
                                 Z( JROW, INCOL+1+J2 ), LDZ, &
                                 U( J2+1, I2+1 ), LDU, ONE, &
                                 WV( 1, 1+I2 ), LDWV )
!
!                    ==== Copy the result back to Z ====
!
                     CALL DLACPY( 'ALL', JLEN, KDU, WV, LDWV, &
                                  Z( JROW, INCOL+1 ), LDZ )
  210             CONTINUE
               END IF
            END IF
         END IF
  220 CONTINUE
!
!     ==== End of DLAQR5 ====
!
      END SUBROUTINE DLAQR5

      SUBROUTINE DLAQR1( N, H, LDH, SR1, SI1, SR2, SI2, V )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   SI1, SI2, SR1, SR2
      INTEGER            LDH, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), V( * )
!     ..
!
!       Given a 2-by-2 or 3-by-3 matrix H, DLAQR1 sets v to a
!       scalar multiple of the first column of the product
!
!       (*)  K = (H - (sr1 + i*si1)*I)*(H - (sr2 + i*si2)*I)
!
!       scaling to avoid overflows and most underflows. It
!       is assumed that either
!
!               1) sr1 = sr2 and si1 = -si2
!           or
!               2) si1 = si2 = 0.
!
!       This is useful for starting double implicit shift bulges
!       in the QR algorithm.
!
!
!       N      (input) integer
!              Order of the matrix H. N must be either 2 or 3.
!
!       H      (input) DOUBLE PRECISION array of dimension (LDH,N)
!              The 2-by-2 or 3-by-3 matrix H in (*).
!
!       LDH    (input) integer
!              The leading dimension of H as declared in
!              the calling procedure.  LDH.GE.N
!
!       SR1    (input) DOUBLE PRECISION
!       SI1    The shifts in (*).
!       SR2
!       SI2
!
!       V      (output) DOUBLE PRECISION array of dimension N
!              A scalar multiple of the first column of the
!              matrix K in (*).
!
!     ================================================================
!     Based on contributions by
!        Karen Braman and Ralph Byers, Department of Mathematics,
!        University of Kansas, USA
!
!     ================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0d0 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   H21S, H31S, S
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS
!     ..
!     .. Executable Statements ..
      IF( N.EQ.2 ) THEN
         S = ABS( H( 1, 1 )-SR2 ) + ABS( SI2 ) + ABS( H( 2, 1 ) )
         IF( S.EQ.ZERO ) THEN
            V( 1 ) = ZERO
            V( 2 ) = ZERO
         ELSE
            H21S = H( 2, 1 ) / S
            V( 1 ) = H21S*H( 1, 2 ) + ( H( 1, 1 )-SR1 )* &
                     ( ( H( 1, 1 )-SR2 ) / S ) - SI1*( SI2 / S )
            V( 2 ) = H21S*( H( 1, 1 )+H( 2, 2 )-SR1-SR2 )
         END IF
      ELSE
         S = ABS( H( 1, 1 )-SR2 ) + ABS( SI2 ) + ABS( H( 2, 1 ) ) + &
             ABS( H( 3, 1 ) )
         IF( S.EQ.ZERO ) THEN
            V( 1 ) = ZERO
            V( 2 ) = ZERO
            V( 3 ) = ZERO
         ELSE
            H21S = H( 2, 1 ) / S
            H31S = H( 3, 1 ) / S
            V( 1 ) = ( H( 1, 1 )-SR1 )*( ( H( 1, 1 )-SR2 ) / S ) - &
                     SI1*( SI2 / S ) + H( 1, 2 )*H21S + H( 1, 3 )*H31S
            V( 2 ) = H21S*( H( 1, 1 )+H( 2, 2 )-SR1-SR2 ) + &
                     H( 2, 3 )*H31S
            V( 3 ) = H31S*( H( 1, 1 )+H( 3, 3 )-SR1-SR2 ) + &
                     H21S*H( 3, 2 )
         END IF
      END IF
      END SUBROUTINE DLAQR1

      SUBROUTINE DLAQR2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ, &
                         IHIZ, Z, LDZ, NS, ND, SR, SI, V, LDV, NH, T, &
                         LDT, NV, WV, LDWV, WORK, LWORK )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV, &
                         LDZ, LWORK, N, ND, NH, NS, NV, NW
      LOGICAL            WANTT, WANTZ
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), SI( * ), SR( * ), T( LDT, * ), &
                         V( LDV, * ), WORK( * ), WV( LDWV, * ), &
                         Z( LDZ, * )
!     ..
!
!     This subroutine is identical to DLAQR3 except that it avoids
!     recursion by calling DLAHQR instead of DLAQR4.
!
!
!     ******************************************************************
!     Aggressive early deflation:
!
!     This subroutine accepts as input an upper Hessenberg matrix
!     H and performs an orthogonal similarity transformation
!     designed to detect and deflate fully converged eigenvalues from
!     a trailing principal submatrix.  On output H has been over-
!     written by a new Hessenberg matrix that is a perturbation of
!     an orthogonal similarity transformation of H.  It is to be
!     hoped that the final version of H has many zero subdiagonal
!     entries.
!
!     ******************************************************************
!     WANTT   (input) LOGICAL
!          If .TRUE., then the Hessenberg matrix H is fully updated
!          so that the quasi-triangular Schur factor may be
!          computed (in cooperation with the calling subroutine).
!          If .FALSE., then only enough of H is updated to preserve
!          the eigenvalues.
!
!     WANTZ   (input) LOGICAL
!          If .TRUE., then the orthogonal matrix Z is updated so
!          so that the orthogonal Schur factor may be computed
!          (in cooperation with the calling subroutine).
!          If .FALSE., then Z is not referenced.
!
!     N       (input) INTEGER
!          The order of the matrix H and (if WANTZ is .TRUE.) the
!          order of the orthogonal matrix Z.
!
!     KTOP    (input) INTEGER
!          It is assumed that either KTOP = 1 or H(KTOP,KTOP-1)=0.
!          KBOT and KTOP together determine an isolated block
!          along the diagonal of the Hessenberg matrix.
!
!     KBOT    (input) INTEGER
!          It is assumed without a check that either
!          KBOT = N or H(KBOT+1,KBOT)=0.  KBOT and KTOP together
!          determine an isolated block along the diagonal of the
!          Hessenberg matrix.
!
!     NW      (input) INTEGER
!          Deflation window size.  1 .LE. NW .LE. (KBOT-KTOP+1).
!
!     H       (input/output) DOUBLE PRECISION array, dimension (LDH,N)
!          On input the initial N-by-N section of H stores the
!          Hessenberg matrix undergoing aggressive early deflation.
!          On output H has been transformed by an orthogonal
!          similarity transformation, perturbed, and the returned
!          to Hessenberg form that (it is to be hoped) has some
!          zero subdiagonal entries.
!
!     LDH     (input) integer
!          Leading dimension of H just as declared in the calling
!          subroutine.  N .LE. LDH
!
!     ILOZ    (input) INTEGER
!     IHIZ    (input) INTEGER
!          Specify the rows of Z to which transformations must be
!          applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N.
!
!     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,IHI)
!          IF WANTZ is .TRUE., then on output, the orthogonal
!          similarity transformation mentioned above has been
!          accumulated into Z(ILOZ:IHIZ,ILO:IHI) from the right.
!          If WANTZ is .FALSE., then Z is unreferenced.
!
!     LDZ     (input) integer
!          The leading dimension of Z just as declared in the
!          calling subroutine.  1 .LE. LDZ.
!
!     NS      (output) integer
!          The number of unconverged (ie approximate) eigenvalues
!          returned in SR and SI that may be used as shifts by the
!          calling subroutine.
!
!     ND      (output) integer
!          The number of converged eigenvalues uncovered by this
!          subroutine.
!
!     SR      (output) DOUBLE PRECISION array, dimension KBOT
!     SI      (output) DOUBLE PRECISION array, dimension KBOT
!          On output, the real and imaginary parts of approximate
!          eigenvalues that may be used for shifts are stored in
!          SR(KBOT-ND-NS+1) through SR(KBOT-ND) and
!          SI(KBOT-ND-NS+1) through SI(KBOT-ND), respectively.
!          The real and imaginary parts of converged eigenvalues
!          are stored in SR(KBOT-ND+1) through SR(KBOT) and
!          SI(KBOT-ND+1) through SI(KBOT), respectively.
!
!     V       (workspace) DOUBLE PRECISION array, dimension (LDV,NW)
!          An NW-by-NW work array.
!
!     LDV     (input) integer scalar
!          The leading dimension of V just as declared in the
!          calling subroutine.  NW .LE. LDV
!
!     NH      (input) integer scalar
!          The number of columns of T.  NH.GE.NW.
!
!     T       (workspace) DOUBLE PRECISION array, dimension (LDT,NW)
!
!     LDT     (input) integer
!          The leading dimension of T just as declared in the
!          calling subroutine.  NW .LE. LDT
!
!     NV      (input) integer
!          The number of rows of work array WV available for
!          workspace.  NV.GE.NW.
!
!     WV      (workspace) DOUBLE PRECISION array, dimension (LDWV,NW)
!
!     LDWV    (input) integer
!          The leading dimension of W just as declared in the
!          calling subroutine.  NW .LE. LDV
!
!     WORK    (workspace) DOUBLE PRECISION array, dimension LWORK.
!          On exit, WORK(1) is set to an estimate of the optimal value
!          of LWORK for the given values of N, NW, KTOP and KBOT.
!
!     LWORK   (input) integer
!          The dimension of the work array WORK.  LWORK = 2*NW
!          suffices, but greater efficiency may result from larger
!          values of LWORK.
!
!          If LWORK = -1, then a workspace query is assumed; DLAQR2
!          only estimates the optimal workspace size for the given
!          values of N, NW, KTOP and KBOT.  The estimate is returned
!          in WORK(1).  No error message related to LWORK is issued
!          by XERBLA.  Neither H nor Z are accessed.
!
!     ================================================================
!     Based on contributions by
!        Karen Braman and Ralph Byers, Department of Mathematics,
!        University of Kansas, USA
!
!     ================================================================
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0d0, ONE = 1.0d0 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   AA, BB, BETA, CC, CS, DD, EVI, EVK, FOO, S, &
                         SAFMAX, SAFMIN, SMLNUM, SN, TAU, ULP
      INTEGER            I, IFST, ILST, INFO, INFQR, J, JW, K, KCOL, &
                         KEND, KLN, KROW, KWTOP, LTOP, LWK1, LWK2, &
                         LWKOPT
      LOGICAL            BULGE, SORTED
!     ..
!     .. External Functions ..
!      DOUBLE PRECISION   DLAMCH
!      EXTERNAL           DLAMCH
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DCOPY, DGEHRD, DGEMM, DLABAD, DLACPY, DLAHQR, &
!                         DLANV2, DLARF, DLARFG, DLASET, DORGHR, DTREXC
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
!     ==== Estimate optimal workspace. ====
!
      JW = MIN( NW, KBOT-KTOP+1 )
      IF( JW.LE.2 ) THEN
         LWKOPT = 1
      ELSE
!
!        ==== Workspace query call to DGEHRD ====
!
         CALL DGEHRD( JW, 1, JW-1, T, LDT, WORK, WORK, -1, INFO )
         LWK1 = INT( WORK( 1 ) )
!
!        ==== Workspace query call to DORGHR ====
!
         CALL DORGHR( JW, 1, JW-1, T, LDT, WORK, WORK, -1, INFO )
         LWK2 = INT( WORK( 1 ) )
!
!        ==== Optimal workspace ====
!
         LWKOPT = JW + MAX( LWK1, LWK2 )
      END IF
!
!     ==== Quick return in case of workspace query. ====
!
      IF( LWORK.EQ.-1 ) THEN
         WORK( 1 ) = DBLE( LWKOPT )
         RETURN
      END IF
!
!     ==== Nothing to do ...
!     ... for an empty active block ... ====
      NS = 0
      ND = 0
      IF( KTOP.GT.KBOT ) &
         RETURN
!     ... nor for an empty deflation window. ====
      IF( NW.LT.1 ) &
         RETURN
!
!     ==== Machine constants ====
!
      SAFMIN = DLAMCH( 'SAFE MINIMUM' )
      SAFMAX = ONE / SAFMIN
      CALL DLABAD( SAFMIN, SAFMAX )
      ULP = DLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( DBLE( N ) / ULP )
!
!     ==== Setup deflation window ====
!
      JW = MIN( NW, KBOT-KTOP+1 )
      KWTOP = KBOT - JW + 1
      IF( KWTOP.EQ.KTOP ) THEN
         S = ZERO
      ELSE
         S = H( KWTOP, KWTOP-1 )
      END IF
!
      IF( KBOT.EQ.KWTOP ) THEN
!
!        ==== 1-by-1 deflation window: not much to do ====
!
         SR( KWTOP ) = H( KWTOP, KWTOP )
         SI( KWTOP ) = ZERO
         NS = 1
         ND = 0
         IF( ABS( S ).LE.MAX( SMLNUM, ULP*ABS( H( KWTOP, KWTOP ) ) ) ) &
              THEN
            NS = 0
            ND = 1
            IF( KWTOP.GT.KTOP ) &
               H( KWTOP, KWTOP-1 ) = ZERO
         END IF
         RETURN
      END IF
!
!     ==== Convert to spike-triangular form.  (In case of a
!     .    rare QR failure, this routine continues to do
!     .    aggressive early deflation using that part of
!     .    the deflation window that converged using INFQR
!     .    here and there to keep track.) ====
!
      CALL DLACPY( 'U', JW, JW, H( KWTOP, KWTOP ), LDH, T, LDT )
      CALL DCOPY( JW-1, H( KWTOP+1, KWTOP ), LDH+1, T( 2, 1 ), LDT+1 )
!
      CALL DLASET( 'A', JW, JW, ZERO, ONE, V, LDV )
      CALL DLAHQR( .true., .true., JW, 1, JW, T, LDT, SR( KWTOP ), &
                   SI( KWTOP ), 1, JW, V, LDV, INFQR )
!
!     ==== DTREXC needs a clean margin near the diagonal ====
!
      DO 10 J = 1, JW - 3
         T( J+2, J ) = ZERO
         T( J+3, J ) = ZERO
   10 CONTINUE
      IF( JW.GT.2 ) &
         T( JW, JW-2 ) = ZERO
!
!     ==== Deflation detection loop ====
!
      NS = JW
      ILST = INFQR + 1
   20 CONTINUE
      IF( ILST.LE.NS ) THEN
         IF( NS.EQ.1 ) THEN
            BULGE = .FALSE.
         ELSE
            BULGE = T( NS, NS-1 ).NE.ZERO
         END IF
!
!        ==== Small spike tip test for deflation ====
!
         IF( .NOT.BULGE ) THEN
!
!           ==== Real eigenvalue ====
!
            FOO = ABS( T( NS, NS ) )
            IF( FOO.EQ.ZERO ) &
               FOO = ABS( S )
            IF( ABS( S*V( 1, NS ) ).LE.MAX( SMLNUM, ULP*FOO ) ) THEN
!
!              ==== Deflatable ====
!
               NS = NS - 1
            ELSE
!
!              ==== Undeflatable.   Move it up out of the way.
!              .    (DTREXC can not fail in this case.) ====
!
               IFST = NS
               CALL DTREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST, WORK, &
                            INFO )
               ILST = ILST + 1
            END IF
         ELSE
!
!           ==== Complex conjugate pair ====
!
            FOO = ABS( T( NS, NS ) ) + SQRT( ABS( T( NS, NS-1 ) ) )* &
                  SQRT( ABS( T( NS-1, NS ) ) )
            IF( FOO.EQ.ZERO ) &
               FOO = ABS( S )
            IF( MAX( ABS( S*V( 1, NS ) ), ABS( S*V( 1, NS-1 ) ) ).LE. &
                MAX( SMLNUM, ULP*FOO ) ) THEN
!
!              ==== Deflatable ====
!
               NS = NS - 2
            ELSE
!
!              ==== Undflatable. Move them up out of the way.
!              .    Fortunately, DTREXC does the right thing with
!              .    ILST in case of a rare exchange failure. ====
!
               IFST = NS
               CALL DTREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST, WORK, &
                            INFO )
               ILST = ILST + 2
            END IF
         END IF
!
!        ==== End deflation detection loop ====
!
         GO TO 20
      END IF
!
!        ==== Return to Hessenberg form ====
!
      IF( NS.EQ.0 ) &
         S = ZERO
!
      IF( NS.LT.JW ) THEN
!
!        ==== sorting diagonal blocks of T improves accuracy for
!        .    graded matrices.  Bubble sort deals well with
!        .    exchange failures. ====
!
         SORTED = .false.
         I = NS + 1
   30    CONTINUE
         IF( SORTED ) &
            GO TO 50
         SORTED = .true.
!
         KEND = I - 1
         I = INFQR + 1
         IF( I.EQ.NS ) THEN
            K = I + 1
         ELSE IF( T( I+1, I ).EQ.ZERO ) THEN
            K = I + 1
         ELSE
            K = I + 2
         END IF
   40    CONTINUE
         IF( K.LE.KEND ) THEN
            IF( K.EQ.I+1 ) THEN
               EVI = ABS( T( I, I ) )
            ELSE
               EVI = ABS( T( I, I ) ) + SQRT( ABS( T( I+1, I ) ) )* &
                     SQRT( ABS( T( I, I+1 ) ) )
            END IF
!
            IF( K.EQ.KEND ) THEN
               EVK = ABS( T( K, K ) )
            ELSE IF( T( K+1, K ).EQ.ZERO ) THEN
               EVK = ABS( T( K, K ) )
            ELSE
               EVK = ABS( T( K, K ) ) + SQRT( ABS( T( K+1, K ) ) )* &
                     SQRT( ABS( T( K, K+1 ) ) )
            END IF
!
            IF( EVI.GE.EVK ) THEN
               I = K
            ELSE
               SORTED = .false.
               IFST = I
               ILST = K
               CALL DTREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST, WORK, &
                            INFO )
               IF( INFO.EQ.0 ) THEN
                  I = ILST
               ELSE
                  I = K
               END IF
            END IF
            IF( I.EQ.KEND ) THEN
               K = I + 1
            ELSE IF( T( I+1, I ).EQ.ZERO ) THEN
               K = I + 1
            ELSE
               K = I + 2
            END IF
            GO TO 40
         END IF
         GO TO 30
   50    CONTINUE
      END IF
!
!     ==== Restore shift/eigenvalue array from T ====
!
      I = JW
   60 CONTINUE
      IF( I.GE.INFQR+1 ) THEN
         IF( I.EQ.INFQR+1 ) THEN
            SR( KWTOP+I-1 ) = T( I, I )
            SI( KWTOP+I-1 ) = ZERO
            I = I - 1
         ELSE IF( T( I, I-1 ).EQ.ZERO ) THEN
            SR( KWTOP+I-1 ) = T( I, I )
            SI( KWTOP+I-1 ) = ZERO
            I = I - 1
         ELSE
            AA = T( I-1, I-1 )
            CC = T( I, I-1 )
            BB = T( I-1, I )
            DD = T( I, I )
            CALL DLANV2( AA, BB, CC, DD, SR( KWTOP+I-2 ), &
                         SI( KWTOP+I-2 ), SR( KWTOP+I-1 ), &
                         SI( KWTOP+I-1 ), CS, SN )
            I = I - 2
         END IF
         GO TO 60
      END IF
!
      IF( NS.LT.JW .OR. S.EQ.ZERO ) THEN
         IF( NS.GT.1 .AND. S.NE.ZERO ) THEN
!
!           ==== Reflect spike back into lower triangle ====
!
            CALL DCOPY( NS, V, LDV, WORK, 1 )
            BETA = WORK( 1 )
            CALL DLARFG( NS, BETA, WORK( 2 ), 1, TAU )
            WORK( 1 ) = ONE
!
            CALL DLASET( 'L', JW-2, JW-2, ZERO, ZERO, T( 3, 1 ), LDT )
!
            CALL DLARF( 'L', NS, JW, WORK, 1, TAU, T, LDT, &
                        WORK( JW+1 ) )
            CALL DLARF( 'R', NS, NS, WORK, 1, TAU, T, LDT, &
                        WORK( JW+1 ) )
            CALL DLARF( 'R', JW, NS, WORK, 1, TAU, V, LDV, &
                        WORK( JW+1 ) )
!
            CALL DGEHRD( JW, 1, NS, T, LDT, WORK, WORK( JW+1 ), &
                         LWORK-JW, INFO )
         END IF
!
!        ==== Copy updated reduced window into place ====
!
         IF( KWTOP.GT.1 ) &
            H( KWTOP, KWTOP-1 ) = S*V( 1, 1 )
         CALL DLACPY( 'U', JW, JW, T, LDT, H( KWTOP, KWTOP ), LDH )
         CALL DCOPY( JW-1, T( 2, 1 ), LDT+1, H( KWTOP+1, KWTOP ), &
                     LDH+1 )
!
!        ==== Accumulate orthogonal matrix in order update
!        .    H and Z, if requested.  (A modified version
!        .    of  DORGHR that accumulates block Householder
!        .    transformations into V directly might be
!        .    marginally more efficient than the following.) ====
!
         IF( NS.GT.1 .AND. S.NE.ZERO ) THEN
            CALL DORGHR( JW, 1, NS, T, LDT, WORK, WORK( JW+1 ), &
                         LWORK-JW, INFO )
            CALL DGEMM( 'N', 'N', JW, NS, NS, ONE, V, LDV, T, LDT, ZERO, &
                        WV, LDWV )
            CALL DLACPY( 'A', JW, NS, WV, LDWV, V, LDV )
         END IF
!
!        ==== Update vertical slab in H ====
!
         IF( WANTT ) THEN
            LTOP = 1
         ELSE
            LTOP = KTOP
         END IF
         DO 70 KROW = LTOP, KWTOP - 1, NV
            KLN = MIN( NV, KWTOP-KROW )
            CALL DGEMM( 'N', 'N', KLN, JW, JW, ONE, H( KROW, KWTOP ), &
                        LDH, V, LDV, ZERO, WV, LDWV )
            CALL DLACPY( 'A', KLN, JW, WV, LDWV, H( KROW, KWTOP ), LDH )
   70    CONTINUE
!
!        ==== Update horizontal slab in H ====
!
         IF( WANTT ) THEN
            DO 80 KCOL = KBOT + 1, N, NH
               KLN = MIN( NH, N-KCOL+1 )
               CALL DGEMM( 'C', 'N', JW, KLN, JW, ONE, V, LDV, &
                           H( KWTOP, KCOL ), LDH, ZERO, T, LDT )
               CALL DLACPY( 'A', JW, KLN, T, LDT, H( KWTOP, KCOL ), &
                            LDH )
   80       CONTINUE
         END IF
!
!        ==== Update vertical slab in Z ====
!
         IF( WANTZ ) THEN
            DO 90 KROW = ILOZ, IHIZ, NV
               KLN = MIN( NV, IHIZ-KROW+1 )
               CALL DGEMM( 'N', 'N', KLN, JW, JW, ONE, Z( KROW, KWTOP ), &
                           LDZ, V, LDV, ZERO, WV, LDWV )
               CALL DLACPY( 'A', KLN, JW, WV, LDWV, Z( KROW, KWTOP ), &
                            LDZ )
   90       CONTINUE
         END IF
      END IF
!
!     ==== Return the number of deflations ... ====
!
      ND = JW - NS
!
!     ==== ... and the number of shifts. (Subtracting
!     .    INFQR from the spike length takes care
!     .    of the case of a rare QR failure while
!     .    calculating eigenvalues of the deflation
!     .    window.)  ====
!
      NS = NS - INFQR
!
!      ==== Return optimal workspace. ====
!
      WORK( 1 ) = DBLE( LWKOPT )
!
!     ==== End of DLAQR2 ====
!
      END SUBROUTINE DLAQR2


      END MODULE lapack_tools
