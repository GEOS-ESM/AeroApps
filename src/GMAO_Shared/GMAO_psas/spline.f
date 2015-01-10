**********************************************************************
*  FILE SPLINE - LAST CHANGE: 06/06/88 (AMS)                         *
*--------------------------------------------------------------------*
*                                                                    *
*      THIS FILE CONTAINS ROUTINES TO COMPUTE CUBIC SPLINES. THE     *
*  ROUTINES ARE BASED ON:                                            *
*                                                                    *
*  FORSYTHE, G. E., M. A. MALCOLN AND C. B. MOLER (1977):            *
*  "COMPUTER METHODS FOR MATHEMATICAL COMPUTATIONS", PrentiCe-Hall.  *
*                                                                    *
**********************************************************************
 
      subroutine SPLINE ( b, c, d, n, x, y )
 
      integer  n
      real     b(n), c(n), d(n), x(n), y(n)
 
**********************************************************************
*   FIRST VERSION: 02/13/86 (AMS)   CURRENT VERSION: 02/13/86 (AMS)  *
*--------------------------------------------------------------------*
*                                                                    *
*     THIS ROUTINE COMPUTES THE COEFFICIENTS  B(I), C(I), D(I),      *
*   I = 1, ..., N FOR A CUBIC INTERPOLATING SPLINE.                  *
*                                                                    *
*    S(X) = Y(I) + B(I) * ( X - X(I) ) + C(I) * ( X - X(I) )**2      *
*                + D(I) * ( X - X(I) )**3                            *
*                                                                    *
*--------------------------------------------------------------------*
*                                                                    *
*   ON INPUT                                                         *
*   --------                                                         *
*   N   ---   NUMBER OF DATA POINTS OR KNOTS ( N .GE. 2 )            *
*   X   ---   THE ABSCISSAS OF THE KNOTS IN STRICTLY INCREASING      *
*             ORDER                                                  *
*   Y   ---   THE ORDINATES OF THE KNOTS                             *
*                                                                    *
*   ON OUTPUT                                                        *
*   ---------                                                        *
*   B, C, D   ---   ARRAYS OF SPLINE COEFFICIENTS AS DEFINED ABOVE.  *
*                                                                    *
*--------------------------------------------------------------------*
*                                                                    *
*   INTERPRETAION:                                                   *
*   -------------                                                    *
*                                                                    *
*   Y(I)   ---   S    ( X(I) )                                       *
*   B(I)   ---   S'   ( X(I) )                                       *
*   C(I)   ---   S''  ( X(I) ) / 2                                   *
*   D(I)   ---   S''' ( X(I) ) / 6  ( DERIVATIVE FROM THE RIGHT )    *
*                                                                    *
*   WHERE ' DENOTES DIFFERENTIATION. THE ACCOMPANYING SUBPROGRAM     *
*   FUNCTION  SEVAL  CAN BE USED TO EVALUATE THE SPLINE.             *
*                                                                    *
**********************************************************************
*
      integer  nm1, ib, i
      real     t
*
      nm1 = n - 1
      if ( n .lt. 2 )  return
      if ( n .lt. 3 )  go to 50
*
*     SET UP TRIDIAGONAL SYSTEM
*     -------------------------
*
*     B = DIAGONAL, D = DIAGONAL, C = RIGHT HAND SIDE
*     -----------------------------------------------
      d(1) = x(2) - x(1)
      c(2) = ( y(2) - y(1) ) / d(1)
      do 10 i = 2, nm1
         d(i)   = x(i+1) - x(i)
         b(i)   = 2.0 * ( d(i-1) + d(i) )
         c(i+1) = ( y(i+1) - y(i) ) / d(i)
         c(i)   = c(i+1) - c(i)
 10   continue
*
*     END CONDITIONS. THIRD DERIVATIVES AT  X(1)  AND  X(N)
*     OBTAINED FROM DIVIDED DIFFERENCES
*     -----------------------------------------------------
      b(1) = - d(1)
      b(n) = - d(n-1)
      c(1) = 0.0
      c(n) = 0.0
      if ( n .eq. 3 )  go to 15
      c(1) = c(3)   / ( x(4) - x(2)   )
     1     - c(2)   / ( x(3)   - x(1)   )
      c(n) = c(n-1) / ( x(n) - x(n-2) )
     1     - c(n-2) / ( x(n-1) - x(n-3) )
      c(1) = c(1) * d(1)**2 / ( x(4) - x(1) )
      c(n) = -c(n) * d(n-1)**2 / ( x(n) - x(n-3) )
*
*     FORWARD ELIMINATION
*     -------------------
 15   do 20 i = 2, n
          t    = d(i-1) / b(i-1)
          b(i) = b(i)  -  t * d(i-1)
          c(i) = c(i)  -  t * c(i-1)
 20   continue
*
*     BACK SUBSTITUTION
*     -----------------
      c(n) = c(n) / b(n)
      do 30 ib = 1, nm1
         i = n - ib
         c(i) = ( c(i) - d(i) * c(i+1) ) / b(i)
 30   continue
*
*     C(I) NOW CONTAINS SIGMA(I)
*     --------------------------
*
*     COMPUTE POLYNOMIAL COEFFICIENTS
*     -------------------------------
      b(n) = ( y(n) - y(nm1) ) / d(nm1)
     1     + d(nm1) * ( c(nm1) + 2.0 * c(n) )
      do 40 i = 1, nm1
         b(i) = ( y(i+1) - y(i) ) / d(i)
     1        - d(i) * ( c(i+1) + 2.0 * c(i) )
         d(i) = ( c(i+1) - c(i) ) / d(i)
         c(i) = 3.0 * c(i)
 40   continue
      c(n) = 3.0 * c(n)
      d(n) = d(n-1)
*
      return
*
 50   b(1) = ( y(2) - y(1) ) / ( x(2) - x(1) )
      c(1) = 0.0
      d(1) = 0.0
      b(2) = b(1)
      c(2) = 0.0
      d(2) = 0.0
*
      return
      end

*............................................................................

      real function SEVAL ( u, x, y, n, b, c, d )

      integer  n
      real     u, x(n), y(n), b(n), c(n), d(n)

**********************************************************************
*   FIRST VERSION: 02/13/86 (AMS)   CURRENT VERSION: 02/13/86 (AMS)  *
*--------------------------------------------------------------------*
*                                                                    *
*       THIS ROUTINE  RETURNS THE CUBIC SPLINE FUNCTION              *
*                                                                    *
*    SEVAL = Y(I) + B(I) * ( U - X(I) ) + C(I) * ( U - X(I) )**2     *
*                 + D(I) * ( U - X(I) )**3                           *
*                                                                    *
*    WHERE  X(I) .LT. U .LT. X(I+1), USING HORNER'S RULE.            *
*       IF U .LT. X(1) THEN I=1 IS USED, IF U .GE. X(N) THEN I=N     *
*    IS USED.                                                        *
*                                                                    *
*--------------------------------------------------------------------*
*                                                                    *
*   ON INPUT                                                         *
*   --------                                                         *
*   N        ---   NUMBER OF DATA POINTS OR KNOTS ( N .GE. 2 )       *
*   U        ---   ABSCISSA AT WHICH THE SPLINE IS TO BE EVALUATED   *
*   X, Y     ---   THE ARRAYS OF ABSCISSAS AND ORDINATES.            *
*   B, C, D  ---   ARRAYS OF SPLINE COEFFICIENTS AS DEFINED ABOVE.   *
*                                                                    *
*      IF  U  IS NOT IN THE SAME INTERVAL AS THE PREVIOUS CALL,      *
*   THEN A BINARY SEARCH IS PERFORMED TO DETERMINE THE PROPER        *
*   INTERVAL.                                                        *
*                                                                    *
**********************************************************************

      save i
      integer  i, j, k
      real     dx

      data i / 1 /

      if ( i .ge. n      ) i = 1
      if ( u .lt. x(i)   ) go to 10
      if ( u .le. x(i+1) ) go to 30

*     BINARY SEARCH
*     -------------
 10   i = 1
      j = n + 1
 20   k = ( i + j ) / 2
      if ( u .lt. x(k) )  j = k
      if ( u .ge. x(k) )  i = k
      if ( j .gt. i+1  )  go to 20

*     EVALUATE SPLINE
*     ---------------
 30   dx    = u - x(i)
      seval = y(i)
     1      + dx * ( b(i) + dx * ( c(i) + dx * d(i) ) )

      return
      end

*............................................................................

      subroutine SEVAL1 ( s, sp, spp, u, x, y, n, b, c, d )

      integer  n
      real     s, sp, spp, u, x(n), y(n), b(n), c(n), d(n)


**********************************************************************
*   FIRST VERSION: 02/17/87 (AMS)   CURRENT VERSION: 02/17/87 (AMS)  *
*--------------------------------------------------------------------*
*                                                                    *
*       THIS ROUTINE  RETURNS THE CUBIC SPLINE FUNCTION              *
*                                                                    *
*    S(U)  = Y(I) + B(I) * ( U - X(I) ) + C(I) * ( U - X(I) )**2     *
*                 + D(I) * ( U - X(I) )**3                           *
*                                                                    *
*    AND ITS RESPECTIVE FIRST AND SECOND DERIVATIVES:                *
*                                                                    *
*                 SP = S'(U)    AND    SPP = S''(U)                  *
*                                                                    *
*    WHERE  X(I) .LT. U .LT. X(I+1), USING HORNER'S RULE.            *
*       IF U .LT. X(1) THEN I=1 IS USED, IF U .GE. X(N) THEN I=N     *
*    IS USED.                                                        *
*                                                                    *
*--------------------------------------------------------------------*
*                                                                    *
*   ON INPUT                                                         *
*   --------                                                         *
*   N        ---   NUMBER OF DATA POINTS OR KNOTS ( N .GE. 2 )       *
*   U        ---   ABSCISSA AT WHICH THE SPLINE IS TO BE EVALUATED   *
*   X, Y     ---   THE ARRAYS OF ABSCISSAS AND ORDINATES.            *
*   B, C, D  ---   ARRAYS OF SPLINE COEFFICIENTS AS DEFINED ABOVE.   *
*                                                                    *
*      IF  U  IS NOT IN THE SAME INTERVAL AS THE PREVIOUS CALL,      *
*   THEN A BINARY SEARCH IS PERFORMED TO DETERMINE THE PROPER        *
*   INTERVAL.                                                        *
*                                                                    *
**********************************************************************

      save     i
      integer  i, j, k
      real     dx
   
      data i / 1 /

      if ( i .ge. n      ) i = 1
      if ( u .lt. x(i)   ) go to 10
      if ( u .le. x(i+1) ) go to 30

*     BINARY SEARCH
*     -------------
 10   i = 1
      j = n + 1
 20   k = ( i + j ) / 2
      if ( u .lt. x(k) )  j = k
      if ( u .ge. x(k) )  i = k
      if ( j .gt. i+1  )  go to 20

*     EVALUATE SPLINE
*     ---------------
 30   dx = u - x(i)
      s   = y(i) + dx * ( b(i) + dx * ( c(i) + dx * d(i) ) )
      sp  = b(i) + dx * ( 2.0 * c(i) + 3.0 * d(i) * dx )
      spp = 2.0 * c(i) + 6.0 * d(i) * dx

      return
      end
 
 
*............................................................................


      subroutine BICSPL ( a1, x1, y1, nx1, ny1,
     1                    a2, x2, y2, nx2, ny2,
     2                    b, c, d )
 
*
*                           A2 contains input array; A1 receives output
*                             dim ( A2 ):   max(NX1,NX2) x max(NY1,NY2)
*                             dim ( A1 ):   max(NX1,NX2) x max(NY1,NY2)
 
      real   a1(*), a2(*)
 
      real   x1(nx1), x2(nx2), y1(ny1), y2(ny2)
*
*                                  these vectors must have dimension
*                                   >= max ( nx2, ny2 )
      real   b(*), c(*), d(*)
 
**********************************************************************
*   FIRST VERSION: 06/06/88         CURRENT VERSION: 06/06/88        *
*--------------------------------------------------------------------*
*                                                                    *
*     THIS ROUTINE COMPUTES INTRPOLATES THE BI-DIMENSIONAL FIELD  A2 *
*   GIVEN ON A GRID  X2 x Y2  TO A FIELD  A1  DEFINED ON A GRID      *
*   X1 x Y1.                                                         *
*                                                                    *
*         A1, A2 should not share storage.   No attempt is made to   *
*   save storage.  A2 is overwritten.                                *
*                                                                    *
**********************************************************************
 
 
*     X-INTERPOLATION
*     ---------------
      do 10 j = 1, ny2
         l1j = ( j - 1 ) * nx2 + 1
         call spline ( b, c, d, nx2, x2, a2(l1j) )
         do 20 i = 1, nx1
            lij = ( j - 1 ) * nx1 + i
            a1(lij) = seval ( x1(i), x2, a2(l1j), nx2, b, c, d )
 20      continue
 10   continue
 
*     TRANSPOSITION
*     -------------
      do 30 i = 1, nx1
         do 30 j = 1, ny2
            lij = ( j - 1 ) * nx1 + i
            lji = ( i - 1 ) * ny2 + j
            a2(lji) = a1(lij)
 30   continue
 
*     Y-INTERPOLATION
*     ---------------
      do 40 i = 1, nx1
         l1i = ( i - 1 ) * ny2 + 1
         call SPLINE ( b, c, d, ny2, y2, a2(l1i) )
         do 50 j = 1, ny1
            lij = ( j - 1 ) * nx1 + i
            a1(lij) = seval ( y1(j), y2, a2(l1i), ny2, b, c, d )
 50      continue
 40   continue
 
 
      return
      end





