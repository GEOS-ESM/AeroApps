
C ###############################################################
C #                                                             #
C #                    THE VECTOR LIDORT MODEL                  #
C #                                                             #
C #  (Vector LInearized Discrete Ordinate Radiative Transfer)   #
C #   -      --         -        -        -         -           #
C #                                                             #
C ###############################################################

C ###############################################################
C #                                                             #
C #  Author :      Robert. J. D. Spurr                          #
C #                                                             #
C #  Address :      RT Solutions, inc.                          #
C #            9 Channing Street                                #
C #             Cambridge, MA 02138, USA                        #
C #            Tel: (617) 492 1183                              #
C #                                                             #
C #  Email :      rtsolutions@verizon.net                       #
C #                                                             #
C #  Versions     :   2.0, 2.2, 2.3, 2.4, 2.4R, 2.4RT           #
C #  Release Date :   December 2005  (2.0)                      #
C #  Release Date :   March 2007     (2.2)                      #
C #  Release Date :   October 2007   (2.3)                      #
C #  Release Date :   December 2008  (2.4)                      #
C #  Release Date :   April/May 2009 (2.4R)                     #
C #  Release Date :   July 2009      (2.4RT)                    #
C #                                                             #
C #       NEW: TOTAL COLUMN JACOBIANS         (2.4)             #
C #       NEW: BPDF Land-surface KERNELS      (2.4R)            #
C #       NEW: Thermal Emission Treatment     (2.4RT)           #
C #                                                             #
C ###############################################################

C    #####################################################
C    #                                                   #
C    #   This Version of VLIDORT comes with a GNU-style  #
C    #   license. Please read the license carefully.     #
C    #                                                   #
C    #####################################################

C ###############################################################
C #                                                             #
C #  This is VLIDORT_AUX.f. This is a collection of utility     #
C #  modules compiled from a variety of sources. This is the    #
C #  only part of the VLIDORT package not written by the Author.#
C #  The subroutines in VLIDORT_AUX are listed below with their #
C #  source of origin (order of appearance).                    #
C #                                                             #
C #      lcstring:              E. Mikusch, DLR 1994            #
C #      len_string:            E. Mikusch, DLR 1994            #
C #      vlidort_status:        R.J.D. Spurr, RTS, 2005         #
C #      vlidort_lapack_error:  R.J.D. Spurr, CFA, 2004         #
C #      vlidort_error_trace:   R.J.D. Spurr, CFA, 2004         #
C #      asymtx:                adapted by RJDS from the        #
C #                             DISORT module ASYMTX, 1988      #
C #      gauleg:                Numerical Recipes, 1992         #
C #      hpsort:                Numerical Recipes, 1992         #
C #      indexx:                Numerical Recipes, 1992         #
C #      cfplgarr:              T.P. Kurosu, 1997 (adapted      #
C #                               Numerical Recipes).           #
C #      xfindpar:              J. Lavagnino, SAO, 1991         #
C #      gfindpar:              J. Lavagnino, SAO, 1991         #
C #                                                             #
C ###############################################################

C  Convert a string to lower-case.  (ASCII character set assumed.)

      SUBROUTINE LCSTRING ( STRING )

C  Modified argument.

      CHARACTER * (*)     STRING

C  Local variables.

      INTEGER         CCC
      INTEGER         III
      INTEGER         OFFSET

C***********************************************************************

      OFFSET = ICHAR ( 'a' ) - ICHAR ( 'A' )
      DO 10 III = 1, LEN ( STRING )
        CCC = ICHAR ( STRING ( III : III ) )
        IF ( CCC .GE. ICHAR ( 'A' ) .AND.
     $          CCC .LE. ICHAR ( 'Z' ) ) THEN
          STRING ( III : III ) = CHAR ( CCC + OFFSET )
        END IF
10    CONTINUE

      END

C

      INTEGER FUNCTION LEN_STRING(STR)
C       ----------------------------------------------------------------
C       Length of string.
C       This function determines the length of the string contained in 
C       the character variable. Therefore it searches the last non-blank
C       character and returnes the position of this value as the length.
C       If a character variable contains blanks only, 0 is returned.
C        2 Feb 1994, Eberhard Mikusch 
C       ----------------------------------------------------------------

C     ----------------------------------------------------------------
C     input argument
C     ----------------------------------------------------------------
      CHARACTER * (*) STR

C     ----------------------------------------------------------------
C     get last non-blank character
C     ----------------------------------------------------------------
      LEN_STRING = LEN(STR)

      DO WHILE ((LEN_STRING .GT.0) .AND. 
     &            (STR(LEN_STRING:LEN_STRING) .EQ. ' '))
        LEN_STRING = LEN_STRING - 1
      END DO 

C     ----------------------------------------------------------------
C     return
C     ----------------------------------------------------------------
      RETURN
      END
C     ----------------------------------------------------------------

C  status check

      SUBROUTINE VLIDORT_STATUS
     &     ( STATUS_INPUTCHECK, STATUS_CALCULATION )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  subroutine arguments
C  --------------------

      INTEGER          STATUS_INPUTCHECK
      INTEGER          STATUS_CALCULATION

      IF ( STATUS_INPUTCHECK. EQ. VLIDORT_WARNING ) THEN
        WRITE(*,*)' WARNING: Bad input from VLIDORT input check'
        WRITE(*,*)'  ------ VLIDORT has executed with internal defaults'
        WRITE(*,*)'  ------ Check the Logfile for details'        
      ELSE IF ( STATUS_INPUTCHECK. EQ. VLIDORT_SERIOUS ) THEN
        WRITE(*,*)' FATAL: Bad input from VLIDORT input check'
        WRITE(*,*)'  ------ VLIDORT will not execute'
        WRITE(*,*)'  ------ Check the Logfile for details'
      ELSE
c        WRITE(*,*)' SUCCESS with VLIDORT input check'
      ENDIF

      IF ( STATUS_INPUTCHECK. NE. VLIDORT_SERIOUS ) THEN
       IF ( STATUS_CALCULATION. NE. VLIDORT_SUCCESS ) THEN
        WRITE(*,*)' FATAL: VLIDORT execution failure'
        WRITE(*,*)'  ------ Check the Logfile for details'
       ELSE
c        WRITE(*,*)' SUCCESS with VLIDORT execution'
       ENDIF
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE VLIDORT_LAPACK_ERROR
     I  ( INFO, LAYER, INPUT_TRACE, STATUS )

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Input arguments
C  ---------------

C  layer and information counter

      INTEGER             INFO, LAYER

C  error tracing

      CHARACTER*(*)       INPUT_TRACE

C  Output
C  ------

      INTEGER             STATUS

C  local variables for error tracing

      CHARACTER*(70)      MAIL, TRACE
      CHARACTER*3         CN, CI
      INTEGER             LENGTH, LEN_STRING
      EXTERNAL            LEN_STRING

C  Initialise

      STATUS = VLIDORT_SUCCESS

C  Check info

      WRITE(CI, '(I3)' ) INFO
      WRITE(CN,'(I3)')LAYER
      LENGTH = LEN_STRING(INPUT_TRACE)
      TRACE = INPUT_TRACE(1:LENGTH)//' call for Layer '//CN

      IF ( INFO .LT. 0 ) THEN
        MAIL  = 'argument i illegal value, for i = '//CI
        STATUS = VLIDORT_SERIOUS
        CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
      ELSE IF ( INFO .GT. 0 ) THEN
        MAIL  = 'Singular matrix, u(i,i)=0, for i = '//CI
        STATUS = VLIDORT_SERIOUS
        CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE VLIDORT_ERROR_TRACE ( MESSAGE, TRACE, STATUS )

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include file of control variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'

C  subroutine arguments

      INTEGER             STATUS
      CHARACTER*(*)       TRACE, MESSAGE

C  first call - open file

      IF ( VLIDORT_ERROR_INIT ) THEN
        OPEN ( VLIDORT_ERRUNIT,
     &           FILE = VLIDORT_ERROR_FILENAME, STATUS = 'UNKNOWN' )
      ENDIF

C  write messages

      IF ( STATUS .EQ. VLIDORT_SERIOUS .OR.
     &       STATUS .EQ. VLIDORT_WARNING ) THEN
        IF ( TRACE .EQ. ' ' ) THEN
          WRITE(VLIDORT_ERRUNIT,'(a)')MESSAGE
        ELSE
          WRITE(VLIDORT_ERRUNIT,'(a)')MESSAGE
          WRITE(VLIDORT_ERRUNIT,'(a)')TRACE
        ENDIF
      ENDIF

C  turn off error initialization flag

      VLIDORT_ERROR_INIT = .FALSE.

C  finish

      RETURN
      END

C

      SUBROUTINE  ASYMTX( AAD, M, IA, IEVEC,
     $                      EVECD, EVALD, IER, WKD,
     $                      MESSAGE, BAD_STATUS )

C    =======  D O U B L E    P R E C I S I O N    V E R S I O N  ======

c       Solves eigenfunction problem for real asymmetric matrix
c       for which it is known a priori that the eigenvalues are real.

c       This is an adaptation of a subroutine EIGRF in the IMSL
c       library to use real instead of complex arithmetic, accounting
c       for the known fact that the eigenvalues and eigenvectors in
c       the discrete ordinate solution are real.  Other changes include
c       putting all the called subroutines in-line, deleting the
c       performance index calculation, updating many DO-loops
c       to Fortran77, and in calculating the machine precision
c       TOL instead of specifying it in a data statement.

c       EIGRF is based primarily on EISPACK routines.  The matrix is
c       first balanced using the parlett-reinsch algorithm.  Then
c       the Martin-Wilkinson algorithm is applied.

C       References:
C          Dongarra, J. and C. Moler, EISPACK -- A Package for Solving
C             Matrix Eigenvalue Problems, in Cowell, ed., 1984:
C             Sources and Development of Mathematical Software,
C             Prentice-Hall, Englewood Cliffs, NJ
C         Parlett and Reinsch, 1969: Balancing a Matrix for Calculation
C             of Eigenvalues and Eigenvectors, Num. Math. 13, 293-304
C         Wilkinson, J., 1965: The Algebraic Eigenvalue Problem,
C             Clarendon Press, Oxford

C   I N P U T    V A R I A B L E S:

C        AAD  :  input asymmetric matrix, destroyed after solved
C        M    :  order of  A
C       IA    :  first dimension of  A
C    IEVEC    :  first dimension of  EVECD

C   O U T P U T    V A R I A B L E S:

C       EVECD :  (unnormalized) eigenvectors of  A 
C                   ( column J corresponds to EVALD(J) )

C       EVALD :  (unordered) eigenvalues of  A ( dimension at least M )

C       IER   :  if .NE. 0, signals that EVALD(IER) failed to converge;
C                   in that case eigenvalues IER+1,IER+2,...,M  are
C                   correct but eigenvalues 1,...,IER are set to zero.

      LOGICAL            BAD_STATUS
      CHARACTER*(*)      MESSAGE

C   S C R A T C H   V A R I A B L E S:

C       WKD    :  WORK AREA ( DIMENSION AT LEAST 2*M )
C+---------------------------------------------------------------------+

C  include file for dimensioning

      INCLUDE '../includes/VLIDORT.PARS'

C  input/output arguments

      INTEGER              M, IA, IEVEC, IER
      DOUBLE PRECISION 
     &         AAD(MAXSTRMSTKS,MAXSTRMSTKS), WKD(4*MAXSTRMSTKS),
     &         EVALD(MAXSTRMSTKS), EVECD(MAXSTRMSTKS,MAXSTRMSTKS)

C  local variables (explicit declaration

      LOGICAL           NOCONV, NOTLAS
      INTEGER              I, J, L, K, KKK, LLL
      INTEGER              N, N1, N2, IN, LB, KA, II
      DOUBLE PRECISION  C1, C2, C3, C4, C5, C6
      DATA             C1 / 0.4375D0 /
      DATA               C2 / 0.5D0 /
      DATA               C3 / 0.75D0 /
      DATA              C4 / 0.95D0 /
      DATA              C5 / 16.0D0 /
      DATA              C6 / 256.0D0 /
      DOUBLE PRECISION  TOL, DISCRI, SGN, RNORM, W, F, G, H, P, Q, R
      DOUBLE PRECISION  REPL, COL, ROW, SCALE, T, X, Z, S, Y, UU, VV

      IER = 0
      BAD_STATUS = .FALSE.
      MESSAGE = ' '

C---4/18/08. Use 1.0d-12 (1.0d-24 has convergence issues in Rayleigh layers)
c      TOL = 0.0000001
      TOL = 1.0D-12
c      TOL = 1.0D-24
c      TOL = 1.0D-20

C       Here change to bypass D1MACH:
C        TOL = D1MACH(4)

      IF ( M.LT.1 .OR. IA.LT.M .OR. IEVEC.LT.M ) THEN
        MESSAGE = 'ASYMTX--bad input variable(s)'
        BAD_STATUS = .TRUE.
        RETURN
      ENDIF
C                           ** HANDLE 1X1 AND 2X2 SPECIAL CASES
      IF ( M.EQ.1 )  THEN
         EVALD(1) = AAD(1,1)
         EVECD(1,1) = 1.0D0
         RETURN
      ELSE IF ( M.EQ.2 )  THEN
         DISCRI = ( AAD(1,1) - AAD(2,2) )**2 + 4.0D0*AAD(1,2)*AAD(2,1)
         IF ( DISCRI.LT.ZERO ) THEN
           MESSAGE = 'ASYMTX--COMPLEX EVALS IN 2X2 CASE'
           BAD_STATUS = .TRUE.
           RETURN
         ENDIF
         SGN = ONE
         IF ( AAD(1,1).LT.AAD(2,2) )  SGN = - ONE
         EVALD(1) = 0.5D0*( AAD(1,1) + AAD(2,2) + SGN*DSQRT(DISCRI) )
         EVALD(2) = 0.5D0*( AAD(1,1) + AAD(2,2) - SGN*DSQRT(DISCRI) )
         EVECD(1,1) = ONE
         EVECD(2,2) = ONE
         IF ( AAD(1,1).EQ.AAD(2,2) .AND.
     $         (AAD(2,1).EQ.ZERO.OR.AAD(1,2).EQ.ZERO) ) THEN
            RNORM = DABS(AAD(1,1))+DABS(AAD(1,2))+
     $                DABS(AAD(2,1))+DABS(AAD(2,2))
            W = TOL * RNORM
            EVECD(2,1) = AAD(2,1) / W
            EVECD(1,2) = - AAD(1,2) / W
         ELSE
            EVECD(2,1) = AAD(2,1) / ( EVALD(1) - AAD(2,2) )
            EVECD(1,2) = AAD(1,2) / ( EVALD(2) - AAD(1,1) )
         ENDIF
         RETURN
      END IF
C                                        ** INITIALIZE OUTPUT VARIABLES
      DO 20 I = 1, M
         EVALD(I) = ZERO
         DO 10 J = 1, M
            EVECD(I,J) = ZERO
10       CONTINUE
         EVECD(I,I) = ONE
20    CONTINUE
C                  ** BALANCE THE INPUT MATRIX AND REDUCE ITS NORM BY
C                  ** DIAGONAL SIMILARITY TRANSFORMATION STORED IN WK;
C                  ** THEN SEARCH FOR ROWS ISOLATING AN EIGENVALUE
C                  ** AND PUSH THEM DOWN
      RNORM = ZERO
      L  = 1
      K  = M

30    KKK = K
         DO 70  J = KKK, 1, -1
            ROW = ZERO
            DO 40 I = 1, K
               IF ( I.NE.J ) ROW = ROW + DABS( AAD(J,I) )
40          CONTINUE
            IF ( ROW.EQ.ZERO ) THEN
               WKD(K) = J
               IF ( J.NE.K ) THEN
                  DO 50 I = 1, K
                     REPL   = AAD(I,J)
                     AAD(I,J) = AAD(I,K)
                     AAD(I,K) = REPL
50                CONTINUE
                  DO 60 I = L, M
                     REPL   = AAD(J,I)
                     AAD(J,I) = AAD(K,I)
                     AAD(K,I) = REPL
60                CONTINUE
               END IF
               K = K - 1
               GO TO 30
            END IF
70       CONTINUE
C                                     ** SEARCH FOR COLUMNS ISOLATING AN
C                                       ** EIGENVALUE AND PUSH THEM LEFT
80    LLL = L
         DO 120 J = LLL, K
            COL = ZERO
            DO 90 I = L, K
               IF ( I.NE.J ) COL = COL + DABS( AAD(I,J) )
90          CONTINUE
            IF ( COL.EQ.ZERO ) THEN
               WKD(L) = J
               IF ( J.NE.L ) THEN
                  DO 100 I = 1, K
                     REPL   = AAD(I,J)
                     AAD(I,J) = AAD(I,L)
                     AAD(I,L) = REPL
100               CONTINUE
                  DO 110 I = L, M
                     REPL   = AAD(J,I)
                     AAD(J,I) = AAD(L,I)
                     AAD(L,I) = REPL
110               CONTINUE
               END IF
               L = L + 1
               GO TO 80
            END IF
120      CONTINUE
C                           ** BALANCE THE SUBMATRIX IN ROWS L THROUGH K
      DO 130 I = L, K
         WKD(I) = ONE
130   CONTINUE

140   NOCONV = .FALSE.
         DO 200 I = L, K
            COL = ZERO
            ROW = ZERO
            DO 150 J = L, K
               IF ( J.NE.I ) THEN
                  COL = COL + DABS( AAD(J,I) )
                  ROW = ROW + DABS( AAD(I,J) )
               END IF
150         CONTINUE
            F = ONE
            G = ROW / C5
            H = COL + ROW
160         IF ( COL.LT.G ) THEN
               F   = F * C5
               COL = COL * C6
               GO TO 160
            END IF
            G = ROW * C5
170         IF ( COL.GE.G ) THEN
               F   = F / C5
               COL = COL / C6
               GO TO 170
            END IF
C                                                         ** NOW BALANCE
            IF ( (COL+ROW)/F .LT. C4*H ) THEN
               WKD(I)  = WKD(I) * F
               NOCONV = .TRUE.
               DO 180 J = L, M
                  AAD(I,J) = AAD(I,J) / F
180            CONTINUE
               DO 190 J = 1, K
                  AAD(J,I) = AAD(J,I) * F
190            CONTINUE
            END IF
200      CONTINUE

      IF ( NOCONV ) GO TO 140
C                                  ** IS -A- ALREADY IN HESSENBERG FORM?
      IF ( K-1 .LT. L+1 ) GO TO 350
C                                   ** TRANSFER -A- TO A HESSENBERG FORM
      DO 290 N = L+1, K-1
         H        = ZERO
         WKD(N+M) = ZERO
         SCALE    = ZERO
C                                                        ** SCALE COLUMN
         DO 210 I = N, K
            SCALE = SCALE + DABS(AAD(I,N-1))
210      CONTINUE
         IF ( SCALE.NE.ZERO ) THEN
            DO 220 I = K, N, -1
               WKD(I+M) = AAD(I,N-1) / SCALE
               H = H + WKD(I+M)**2
220         CONTINUE
            G = - SIGN( DSQRT(H), WKD(N+M) )
            H = H - WKD(N+M) * G
            WKD(N+M) = WKD(N+M) - G
C                                                 ** FORM (I-(U*UT)/H)*A
            DO 250 J = N, M
               F = ZERO
               DO 230  I = K, N, -1
                  F = F + WKD(I+M) * AAD(I,J)
230            CONTINUE
               DO 240 I = N, K
                  AAD(I,J) = AAD(I,J) - WKD(I+M) * F / H
240            CONTINUE
250         CONTINUE
C                                    ** FORM (I-(U*UT)/H)*A*(I-(U*UT)/H)
            DO 280 I = 1, K
               F = ZERO
               DO 260  J = K, N, -1
                  F = F + WKD(J+M) * AAD(I,J)
260            CONTINUE
               DO 270 J = N, K
                  AAD(I,J) = AAD(I,J) - WKD(J+M) * F / H
270            CONTINUE
280         CONTINUE
            WKD(N+M)  = SCALE * WKD(N+M)
            AAD(N,N-1) = SCALE * G
         END IF
290   CONTINUE

      DO 340  N = K-2, L, -1
         N1 = N + 1
         N2 = N + 2
         F  = AAD(N1,N)
         IF ( F.NE.ZERO ) THEN
            F  = F * WKD(N1+M)
            DO 300 I = N2, K
               WKD(I+M) = AAD(I,N)
300         CONTINUE
            IF ( N1.LE.K ) THEN
               DO 330 J = 1, M
                  G = ZERO
                  DO 310 I = N1, K
                     G = G + WKD(I+M) * EVECD(I,J)
310               CONTINUE
                  G = G / F
                  DO 320 I = N1, K
                     EVECD(I,J) = EVECD(I,J) + G * WKD(I+M)
320               CONTINUE
330            CONTINUE
            END IF
         END IF
340   CONTINUE

350   CONTINUE
      N = 1
      DO 370 I = 1, M
         DO 360 J = N, M
            RNORM = RNORM + DABS(AAD(I,J))
360      CONTINUE
         N = I
         IF ( I.LT.L .OR. I.GT.K ) EVALD(I) = AAD(I,I)
370   CONTINUE
      N = K
      T = ZERO
C                                         ** SEARCH FOR NEXT EIGENVALUES
380   IF ( N.LT.L ) GO TO 530
      IN = 0
      N1 = N - 1
      N2 = N - 2
C                          ** LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
390   CONTINUE
      DO 400 I = L, N
         LB = N+L - I
         IF ( LB.EQ.L ) GO TO 410
         S = DABS( AAD(LB-1,LB-1) ) + DABS( AAD(LB,LB) )
         IF ( S.EQ.ZERO ) S = RNORM
         IF ( DABS(AAD(LB,LB-1)) .LE. TOL*S ) GO TO 410
400   CONTINUE

410   X = AAD(N,N)
      IF ( LB.EQ.N ) THEN
C                                        ** ONE EIGENVALUE FOUND
         AAD(N,N)  = X + T
         EVALD(N) = AAD(N,N)
         N = N1
         GO TO 380
      END IF

      Y = AAD(N1,N1)
      W = AAD(N,N1) * AAD(N1,N)
      IF ( LB.EQ.N1 ) THEN
C                                        ** TWO EIGENVALUES FOUND
         P = (Y-X) * C2
         Q = P**2 + W
         Z = DSQRT( DABS(Q) )
         AAD(N,N) = X + T
         X = AAD(N,N)
         AAD(N1,N1) = Y + T
C                                        ** REAL PAIR
         Z = P + SIGN(Z,P)
         EVALD(N1) = X + Z
         EVALD(N)  = EVALD(N1)
         IF ( Z.NE.ZERO ) EVALD(N) = X - W / Z
         X = AAD(N,N1)
C                                  ** EMPLOY SCALE FACTOR IN CASE
C                                  ** X AND Z ARE VERY SMALL
         R = SQRT( X*X + Z*Z )
         P = X / R
         Q = Z / R
C                                             ** ROW MODIFICATION
         DO 420 J = N1, M
            Z = AAD(N1,J)
            AAD(N1,J) = Q * Z + P * AAD(N,J)
            AAD(N,J)  = Q * AAD(N,J) - P * Z
420      CONTINUE
C                                             ** COLUMN MODIFICATION
         DO 430 I = 1, N
            Z = AAD(I,N1)
            AAD(I,N1) = Q * Z + P * AAD(I,N)
            AAD(I,N)  = Q * AAD(I,N) - P * Z
430      CONTINUE
C                                          ** ACCUMULATE TRANSFORMATIONS
         DO 440 I = L, K
            Z = EVECD(I,N1)
            EVECD(I,N1) = Q * Z + P * EVECD(I,N)
            EVECD(I,N)  = Q * EVECD(I,N) - P * Z
440      CONTINUE

         N = N2
         GO TO 380
      END IF

      IF ( IN.EQ.30 ) THEN
C                    ** NO CONVERGENCE AFTER 30 ITERATIONS; SET ERROR
C                    ** INDICATOR TO THE INDEX OF THE CURRENT EIGENVALUE
         IER = N
         GO TO 670
      END IF
C                                                          ** FORM SHIFT
      IF ( IN.EQ.10 .OR. IN.EQ.20 ) THEN
         T = T + X
         DO 450 I = L, N
            AAD(I,I) = AAD(I,I) - X
450      CONTINUE
         S = DABS(AAD(N,N1)) + DABS(AAD(N1,N2))
         X = C3 * S
         Y = X
         W = - C1 * S**2
      END IF

      IN = IN + 1
C                ** LOOK FOR TWO CONSECUTIVE SMALL SUB-DIAGONAL ELEMENTS

      DO 460 J = LB, N2
         I = N2+LB - J
         Z = AAD(I,I)
         R = X - Z
         S = Y - Z
         P = ( R * S - W ) / AAD(I+1,I) + AAD(I,I+1)
         Q = AAD(I+1,I+1) - Z - R - S
         R = AAD(I+2,I+1)
         S = DABS(P) + DABS(Q) + DABS(R)
         P = P / S
         Q = Q / S
         R = R / S
         IF ( I.EQ.LB ) GO TO 470
         UU = DABS( AAD(I,I-1) ) * ( DABS(Q) + DABS(R) )
         VV = DABS(P)*(DABS(AAD(I-1,I-1))+DABS(Z)+DABS(AAD(I+1,I+1)))
         IF ( UU .LE. TOL*VV ) GO TO 470
460   CONTINUE

470   CONTINUE
      AAD(I+2,I) = ZERO
      DO 480 J = I+3, N
         AAD(J,J-2) = ZERO
         AAD(J,J-3) = ZERO
480   CONTINUE

C             ** DOUBLE QR STEP INVOLVING ROWS K TO N AND COLUMNS M TO N

      DO 520 KA = I, N1
         NOTLAS = KA.NE.N1
         IF ( KA.EQ.I ) THEN
            S = SIGN( DSQRT( P*P + Q*Q + R*R ), P )
            IF ( LB.NE.I ) AAD(KA,KA-1) = - AAD(KA,KA-1)
         ELSE
            P = AAD(KA,KA-1)
            Q = AAD(KA+1,KA-1)
            R = ZERO
            IF ( NOTLAS ) R = AAD(KA+2,KA-1)
            X = DABS(P) + DABS(Q) + DABS(R)
            IF ( X.EQ.ZERO ) GO TO 520
            P = P / X
            Q = Q / X
            R = R / X
            S = SIGN( DSQRT( P*P + Q*Q + R*R ), P )
            AAD(KA,KA-1) = - S * X
         END IF
         P = P + S
         X = P / S
         Y = Q / S
         Z = R / S
         Q = Q / P
         R = R / P
C                                                    ** ROW MODIFICATION
         DO 490 J = KA, M
            P = AAD(KA,J) + Q * AAD(KA+1,J)
            IF ( NOTLAS ) THEN
               P = P + R * AAD(KA+2,J)
               AAD(KA+2,J) = AAD(KA+2,J) - P * Z
            END IF
            AAD(KA+1,J) = AAD(KA+1,J) - P * Y
            AAD(KA,J)   = AAD(KA,J)   - P * X
490      CONTINUE
C                                                 ** COLUMN MODIFICATION
         DO 500 II = 1, MIN0(N,KA+3)
            P = X * AAD(II,KA) + Y * AAD(II,KA+1)
            IF ( NOTLAS ) THEN
               P = P + Z * AAD(II,KA+2)
               AAD(II,KA+2) = AAD(II,KA+2) - P * R
            END IF
            AAD(II,KA+1) = AAD(II,KA+1) - P * Q
            AAD(II,KA)   = AAD(II,KA) - P
500      CONTINUE
C                                          ** ACCUMULATE TRANSFORMATIONS
         DO 510 II = L, K
            P = X * EVECD(II,KA) + Y * EVECD(II,KA+1)
            IF ( NOTLAS ) THEN
               P = P + Z * EVECD(II,KA+2)
               EVECD(II,KA+2) = EVECD(II,KA+2) - P * R
            END IF
            EVECD(II,KA+1) = EVECD(II,KA+1) - P * Q
            EVECD(II,KA)   = EVECD(II,KA) - P
510      CONTINUE

520   CONTINUE
      GO TO 390
C                     ** ALL EVALS FOUND, NOW BACKSUBSTITUTE REAL VECTOR
530   CONTINUE
      IF ( RNORM.NE.ZERO ) THEN
         DO 560  N = M, 1, -1
            N2 = N
            AAD(N,N) = ONE
            DO 550  I = N-1, 1, -1
               W = AAD(I,I) - EVALD(N)
               IF ( W.EQ.ZERO ) W = TOL * RNORM
               R = AAD(I,N)
               DO 540 J = N2, N-1
                  R = R + AAD(I,J) * AAD(J,N)
540            CONTINUE
               AAD(I,N) = - R / W
               N2 = I
550         CONTINUE
560      CONTINUE
C                      ** END BACKSUBSTITUTION VECTORS OF ISOLATED EVALS

         DO 580 I = 1, M
            IF ( I.LT.L .OR. I.GT.K ) THEN
               DO 570 J = I, M
                  EVECD(I,J) = AAD(I,J)
570            CONTINUE
            END IF
580      CONTINUE
C                                   ** MULTIPLY BY TRANSFORMATION MATRIX
         IF ( K.NE.0 ) THEN
            DO 600  J = M, L, -1
               DO 600 I = L, K
                  Z = ZERO
                  DO 590 N = L, MIN0(J,K)
                     Z = Z + EVECD(I,N) * AAD(N,J)
590               CONTINUE
                  EVECD(I,J) = Z
600         CONTINUE
         END IF

      END IF

      DO 620 I = L, K
         DO 620 J = 1, M
            EVECD(I,J) = EVECD(I,J) * WKD(I)
620   CONTINUE
C                           ** INTERCHANGE ROWS IF PERMUTATIONS OCCURRED
      DO 640  I = L-1, 1, -1
         J = INT(WKD(I))
         IF ( I.NE.J ) THEN
            DO 630 N = 1, M
               REPL       = EVECD(I,N)
               EVECD(I,N) = EVECD(J,N)
               EVECD(J,N) = REPL
630         CONTINUE
         END IF
640   CONTINUE

      DO 660 I = K+1, M
         J = INT(WKD(I))
         IF ( I.NE.J ) THEN
            DO 650 N = 1, M
               REPL       = EVECD(I,N)
               EVECD(I,N) = EVECD(J,N)
               EVECD(J,N) = REPL
650         CONTINUE
         END IF
660   CONTINUE
C                    
  670 CONTINUE

      RETURN
      END

C

      SUBROUTINE GAULEG(X1,X2,X,W,N)
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
      END

C

      SUBROUTINE HPSORT(N,RA)
      INTEGER             N, L, IR, I, J
      DOUBLE PRECISION RA(N), RRA
      L=N/2+1
      IR=N
10      CONTINUE
            IF(L.GT.1)THEN
                  L=L-1
                  RRA=RA(L)
            ELSE
                  RRA=RA(IR)
                  RA(IR)=RA(1)
                  IR=IR-1
                  IF(IR.EQ.1)THEN
                        RA(1)=RRA
                        RETURN
                  ENDIF
            ENDIF
            I=L
            J=L+L
20            IF(J.LE.IR)THEN
                  IF(J.LT.IR)THEN
                        IF(RA(J).LT.RA(J+1))J=J+1
                  ENDIF
                  IF(RRA.LT.RA(J))THEN
                        RA(I)=RA(J)
                        I=J
                        J=J+J
                  ELSE
                        J=IR+1
                  ENDIF
            GO TO 20
            ENDIF
            RA(I)=RRA
      GO TO 10
      END

C

      SUBROUTINE INDEXX(N,ARRIN,INDX)
      INTEGER             N
      INTEGER             INDX(N)
      DOUBLE PRECISION ARRIN(N)
      INTEGER             I, J, L, IR, INDXT
      DOUBLE PRECISION Q

      DO J=1,N
            INDX(J)=J
      ENDDO
      L=N/2+1
      IR=N
10      CONTINUE
            IF(L.GT.1)THEN
                  L=L-1
                  INDXT=INDX(L)
                  Q=ARRIN(INDXT)
            ELSE
                  INDXT=INDX(IR)
                  Q=ARRIN(INDXT)
                  INDX(IR)=INDX(1)
                  IR=IR-1
                  IF(IR.EQ.1)THEN
                        INDX(1)=INDXT
                        RETURN
                  ENDIF
            ENDIF
            I=L
            J=L+L
20            IF(J.LE.IR)THEN
                  IF(J.LT.IR)THEN
                    IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
                  ENDIF
                  IF(Q.LT.ARRIN(INDX(J)))THEN
                        INDX(I)=INDX(J)
                        I=J
                        J=J+J
                  ELSE
                        J=IR+1
                  ENDIF
            GO TO 20
            ENDIF
            INDX(I)=INDXT
      GO TO 10
      END

C

      Subroutine CFPLGARR (MAXMOM, L, M, X, CFPLG)

C-------- This Routine Calculates  -----------

C ###########################################################
C #                                                         #
C #          [(l - m)! / (l + m)!]^{1/2} P(l,m)(x)          #
C #                                                         #
C ###########################################################

c     Based on the PLG routine from the numerical receipes,
c     this revised routine calculates all legendre polynomials
c     from Plg(m,m) to Plg(l,m).
c     The maximum order for the Plg is OMAX


c     ---------------
c     Input Variables
c     ---------------

      Integer            maxmom
      Integer            l, m       ! order of legendre polynomial
      Double Precision   x          ! abscissa


c     ---------------
c     Local Variables
c     ---------------

      Integer            i, ll

      Double Precision   fact, fll, fmm, fmmp1, somx2

c     ---------------
c     Output Variables
c     ---------------

      Double Precision   cfplg (0:maxmom) ! array of Plgs

C     --- Compute F(M,M)  ---

      fmm = 1.0d0

      If (m. gt. 0) Then
         somx2 = Dsqrt( (1.0d0-x) * (1.0d0+x) )
         fact = 1.0d0
         Do i = 1, m
            fmm = -fmm*fact*somx2 / Dsqrt( DBLE( (2*i)*(2*i-1) ) )
            fact = fact + 2.0d0
         End Do
      End If

      cfplg(m) = fmm

C     --- Compute F(M+1,M)  ---

      If (l .gt. m) Then
         fmmp1 = x * Dsqrt( DBLE( 2*m+1 ) ) * fmm
         cfplg(m+1) = fmmp1

C     --- Compute F(L,M) ---

         If ( l .gt. (m+1) ) Then
         
            Do ll = m+2, l
               fll = (fmmp1 * x * DBLE(2*ll-1) - 
     &              fmm * DBLE(ll+m-1) 
     &              * Dsqrt( (DBLE( ll-m-1 ) ) / DBLE( ll+m-1) ) )
     &              / Dsqrt( DBLE( (ll+m)*(ll-m) ) )
               fmm   = fmmp1
               fmmp1 = fll
               cfplg(ll) = fmmp1
            End Do

         End If

      End If

      Return
      End  



C  FINDPAR input subroutines.

C----- The basic idea

C  FINDPAR is a poor man's RDPAR.  It implements some of the
C central features of RDPAR --- the ability to find a named set
C of data values in an input file, the relaxation of the need for
C data in the input file to appear in the order in which it's
C used, the provision for comments in the input file --- in
C (almost) standard FORTRAN 77.

C  Among the features of RDPAR that are not implemented in FINDPAR:
C all redirection of input (from one parameter name to another, from
C file to terminal, from file to command line); freer-format data
C specification on input; tracing options.

C----- Using FINDPAR

C  To use FINDPAR, you open your input file using ordinary FORTRAN
C procedures.  You will also use ordinary FORTRAN to read data
C values.  The only difference is that you can use the FINDPAR call
C at any point to position the file to a point after a line containing
C a parameter name.  These parameter names help to document what the
C input values are for; they also don't have to appear in the order in
C which the program uses them.  An input file might look like this:

C ! Parameters for July 16 run.

C NLAYERS
C 5

C LAYER PARAMETERS
C 5,1013,296
C .001,.001,.001,.001,.001,.001,.001,.001

C MOLECULE FLAGS
C t,t,t,t,t,t,t,t

C  By convention, ! or a space introduces a comment line.  A
C parameter name appears on a line of its own, followed by the data
C lines to be read using FORTRAN statements.  Here is a sample of
C FORTRAN code to read this parameter file:

C      OPEN (UNIT = IN_UNIT, FILE = INPUT, STATUS = 'OLD')
C      CALL FINDPAR ( IN_UNIT, 'NLAYERS' )
C      READ (IN_UNIT, *) NLAYERS
C      CALL FINDPAR ( IN_UNIT, 'LAYER PARAMETERS' )
C      READ (IN_UNIT, *) Z, PRESS, TEMP
C      READ (IN_UNIT, *) BRO, CLO, HCHO, NO2, O2, O3, OCLO, SO2
C      CALL FINDPAR ( IN_UNIT, 'MOLECULE FLAGS' )
C      READ (IN_UNIT, *) IF_BRO, IF_CLO, IF_HCHO, IF_NO2, IF_O2, IF_O3,
C     $ IF_OCLO, IF_SO2
C      CLOSE (UNIT = IN_UNIT)

C----- Details of parameter file format

C The parameter file consists of a sequence of comments and
C parameter-name/data chunks, in any order (except as noted below).

C Comments are essentially anything that doesn't happen to get matched
C as a parameter name or read as data.  However, by convention comments
C are introduced by an exclamation point, to make absolutely sure that
C they don't get mistaken for parameter names.  It is also conventional
C to ``comment out'' a parameter name that you don't want to match by
C simply indenting it one space.

C Parameter names can theoretically be almost any sequence of
C characters.  However, FINDPAR ignores trailing spaces both in the
C parameter name you give it and in lines it reads from the parameter
C file.  And, by convention, parameter names don't begin with spaces.
C Parameter names may contain embedded spaces.  Parameter name matching
C is case-insensitive but otherwise exact: if you use odd things
C like control characters in a parameter name they have to match
C exactly, too, so it's usually best to avoid control characters and
C tabs in parameter names.

C The data lines that follow a line with a parameter name are ordinary
C FORTRAN input lines.  FINDPAR has no influence on how they're handled,
C and all the ordinary rules about FORTRAN READ statements apply.

C Comments may not appear on the same lines as parameter names or data,
C and they can't appear within parameter name/data chunks.  When FINDPAR
C locates a parameter name it positions the file at the line following
C that name; your program will then use ordinary FORTRAN READ statements
C to read the data, and will probably be confused by comments that appear
C at that point.

C There's a maximum length for input lines and parameter names, set
C by the parameter MAXLINE in the code.

C By convention, a parameter name that ends in ? precedes an
C input line that contains a single logical variable (T or F).

C----- Repeated parameters

C Sometimes you want to read several parameters with the same name:
C for example, parameters for each layer of the atmosphere, where there
C might be a large number of layers.  In this case the order of
C parameters in the file does matter, since you usually just want to
C put the numbers in the file in some natural order rather than giving
C a number in each input line specifying where it appears in the input
C order.

C The normal FINDPAR approach is to organize these parameters as
C follows: begin with some parameter that appears only once in the file
C (the number of layers, for example); then follow it with the instances
C of the repeated parameter and associated data.  The only-once
C parameter can even be a dummy, with no data following; its importance
C is that reading it gets the parameter file positioned at the right
C point.

C The other problem here is knowing when a list of repeated parameters
C is over.  FINDPAR doesn't have any of the tricks for doing this that
C RDPAR has; you must either know exactly how many instances of the
C repeated parameter there are going to be (by reading some once-only
C parameter that tells you), or else you need some special numbers to
C flag the end of the list (zeros, negative numbers, etc.).  You can't
C just keep looking for the same parameter name indefinitely, because
C FINDPAR will just rewind the file for you and loop through its
C contents forever.

C----- Errors and XFINDPAR

C If FINDPAR can't find a parameter, it prints an error message
C and returns.  Your READ statements following the FINDPAR call are
C very likely to run into errors of their own at this point, since
C you'll be positioned at some random point in the file (usually at
C the beginning, in this version, but that isn't true in all cases).

C If you want to do something more intelligent about such errors,
C you can use XFINDPAR, which is a function returning a logical value:
C .true. if the parameter name was found, .false. if not; XFINDPAR
C doesn't display FINDPAR's error message when a parameter name cannot
C be found.

C Both FINDPAR and XFINDPAR can run into ordinary FORTRAN input errors
C as they read the file: no special action is taken on these---the
C system's default action, whatever that is, occurs.

C       9/10/91         John Lavagnino

C  Find a parameter name in a parameter file.

C        SUBROUTINE FINDPAR ( NUNIT, PARNAME )

C  Input arguments.

C        INTEGER         NUNIT
C        CHARACTER * (*) PARNAME

C  Local variables.

C        LOGICAL         XFINDPAR
C        EXTERNAL        XFINDPAR

C  No local variables need to be SAVEd.

C***********************************************************************

C        IF ( .NOT. XFINDPAR ( NUNIT, PARNAME ) ) THEN
C          WRITE (*, *) 'FINDPAR error'
C          WRITE (*, *) '   Unit number ', NUNIT
C          WRITE (*, *) '   Parameter name ', PARNAME
C        END IF

C        END
C
C  Find a parameter name in a paramter file: return .true. if found,
C .false. if not.

        LOGICAL FUNCTION XFINDPAR ( NUNIT, PARNAME )

C  Parameters.

        INTEGER         MAXLINE
        PARAMETER       ( MAXLINE = 132 )

C  Input arguments.

        INTEGER         NUNIT
        CHARACTER * (*) PARNAME

C  Local variables.

        INTEGER                 III
        INTEGER                 PARLEN
        INTEGER                 START
        INTEGER                 END
        LOGICAL                 ENDSEEN
        CHARACTER * (MAXLINE)   LINE
        CHARACTER * (MAXLINE)   NAME

C  No local variables need to be SAVEd.

C***********************************************************************

C  Determine the length of the parameter name.

        DO 10 III = LEN ( PARNAME ), 1, -1
          IF ( PARNAME ( III : III ) .NE. ' ' ) THEN
            PARLEN = III
            GO TO 20
          END IF
10      CONTINUE

C  If we get to here, then name contains nothing but blanks.
C  We just return, claiming success, in such a case.

        GO TO 500

C  If we get here, then there's a non-null parameter name; but it
C  might still be too long, in which case we always return
C  signaling failure, since we couldn't ever find such a name in the
C  file.

20      CONTINUE
        IF ( PARLEN  .GT.  MAXLINE ) GO TO 400

C  Convert the name to lower-case.

        NAME = PARNAME ( 1 : PARLEN )
        CALL LCSTRING ( NAME ( 1 : PARLEN ) )

C  Top of main loop.

        ENDSEEN = .FALSE.

100     CONTINUE

          LINE = ' '
          READ ( UNIT = NUNIT, FMT = '(A)', END = 200 ) LINE
          CALL LCSTRING ( LINE ( 1 : PARLEN ) )
          IF ( LINE ( 1 : PARLEN ) .NE. NAME ( 1 : PARLEN ) )
     $          GO TO 100

          START = PARLEN + 1
          END = MAXLINE
          IF ( START  .GT.  END ) GO TO 500
          IF ( LINE ( START : END ) .EQ. ' ' ) GO TO 500
          GO TO 100

C  End-of-file branch.

200       CONTINUE
          REWIND ( UNIT = NUNIT )
          IF ( ENDSEEN ) THEN
            GO TO 400
          ELSE
            ENDSEEN = .TRUE.
            GO TO 100
          END IF

C  End of loop: failure, no parameter name found.

400     CONTINUE
        XFINDPAR = .FALSE.
        RETURN

C  End of loop: successful location of parameter name.

500     CONTINUE
        XFINDPAR = .TRUE.
        RETURN

        END
C


C
C  Find a parameter name in a parameter file, with optional added
C prefix on parameter name.  Error messages are displayed but
C instead of stopping the program a logical flag is returned.
C  Returns .TRUE. if some form of the parameter name was found;
C .FALSE. if not.
C ERROR and GFINDPAR are initialised for a successful search.
C ERROR is just the negation of GFINDPAR. Error messages are removed. 

        LOGICAL FUNCTION GFINDPAR (NUNIT, PREFIX,
     $                                    ERROR, PARNAME)


C  Input arguments.  If PREFIX is something other than a bunch of
C blanks, it will be added onto PARNAME, with a blank to separate
C it, for the FINDPAR call.  If that FINDPAR call fails, then
C PARNAME without prefix is tried instead, and only if that call fails
C is there an error message.  This means that unprefixed versions of
C the parameter name act as defaults.
C  If PREFIX is a bunch of blanks, this works much as FINDPAR does.
C  (This is assuming that PREFIX has no trailing blanks.  Probably
C it should be checking for them and trimming them if necessary.)

        INTEGER         NUNIT
        CHARACTER * (*) PREFIX
        CHARACTER * (*) PARNAME

C  Modified argument.  If all goes well, this is not changed.  If
C the parameter name can't be found, it's set to .TRUE.

        LOGICAL         ERROR

C  Local variables.

        CHARACTER * 132 LINE

        LOGICAL         XFINDPAR
        EXTERNAL        XFINDPAR

C***********************************************************************

C  Initialise return values.

        GFINDPAR = .TRUE.
        ERROR    = .FALSE.

C  Compose search line for first search.

        IF ( PREFIX  .EQ.  ' ' ) THEN
          LINE = PARNAME
        ELSE
          LINE = PREFIX // ' ' // PARNAME
        END IF

C  Try the parameter name with prefix.

        IF (XFINDPAR (NUNIT, LINE)) THEN
          RETURN
        END IF

C  If that doesn't work, try it without the prefix, if the prefix
C was non-null.

        IF (PREFIX  .NE.  ' ') THEN
          IF (XFINDPAR (NUNIT, PARNAME)) THEN
            RETURN
          END IF
        END IF

C  If we got here, we just couldn't find the parameter name.

        GFINDPAR = .FALSE.
        ERROR    = .TRUE.

        END
