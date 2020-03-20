C ###########################################################
C #                                                         #
C #                    THE LIDORT FAMILY                    #
C #                                                         #
C #      (LInearized Discrete Ordinate Radiative Transfer)  #
C #       --         -        -        -         -          #
C #                                                         #
C ###########################################################

C ###########################################################
C #                                                         #
C #  Author :      Robert. J. D. Spurr                      #
C #                                                         #
C #  Address :     RT Solutions, Inc.                       #
C #                9 Channing Street                        #
C #                Cambridge, MA 02138, USA                 #
C #                                                         #
C #  Tel:          (617) 492 1183                           #
C #  Email :        rtsolutions@verizon.net                 #
C #                                                         #
C #  This Version :   3.3                                   #
C #  Release Date :   September 2007                        #
C #                                                         #
C #       NEW: THERMAL SUPPLEMENT INCLUDED    (3.2)         #
C #       NEW: OUTGOING SPHERICITY CORRECTION (3.2)         #
C #       NEW: TOTAL COLUMN JACOBIANS         (3.3)         #
C #                                                         #
C ###########################################################

C    #####################################################
C    #                                                   #
C    #   This Version of LIDORT comes with a GNU-style   #
C    #   license. Please read the license carefully.     #
C    #                                                   #
C    #####################################################

C ###############################################################
C #                                                             #
C # Subroutines in this Module                                  #
C #                                                             #
C #            LISPARSE_FUNCTION_PLUS                           #
C #            LIDENSE_FUNCTION_PLUS                            #
C #            HAPKE_FUNCTION_PLUS                              #
C #            RAHMAN_FUNCTION_PLUS                             #
C #            COXMUNK_FUNCTION_PLUS                            #
C #            COXMUNK_FUNCTION_PLUS_DB                         #
C #                                                             #
C ###############################################################

C

      SUBROUTINE HAPKE_FUNCTION_PLUS 
     I       ( MAXPARS, NPARS, PARS, DO_DERIV_PARS,
     I         XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,
     O         HAPKE_KERNEL, HAPKE_DERIVATIVES )

C  include file of constants amd dimensions

      INCLUDE '../includes/LIDORT.PARS'

C  Subroutine arguments

      INTEGER          MAXPARS, NPARS
      DOUBLE PRECISION PARS ( MAXPARS )
      LOGICAL          DO_DERIV_PARS ( MAXPARS )
      DOUBLE PRECISION XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION HAPKE_KERNEL
      DOUBLE PRECISION HAPKE_DERIVATIVES ( MAXPARS )

C  Hapke Kernel function.
C    - New version, Fresh Coding
C    - Old version uses DISORT code; for validation.

C  input variables:

c    XI, SXI  : Cosine/Sine of angle of reflection (positive)
c    XJ, SXJ  : Cosine/Sine of angle of incidence (positive)
c    XPHI     : Difference of azimuth angles of incidence and reflection
c    PARS(1)  : single scattering albedo in Hapke's BDR model
c    PARS(2)  : angular width parameter of opposition effect in Hapke's model
c    PARS(3)  : Empirical hot spot multiplier

C  local variables
c    B0_EMPIR : empirical factor to account for the finite size of
c               particles in Hapke's BDR model
c    B_HOT    : term that accounts for the opposition effect
c               (retroreflectance, hot spot) in Hapke's BDR model
c    CTHETA   : cosine of phase angle in Hapke's BDR model
c    GAMMA    : albedo factor in Hapke's BDR model
c    PHASE    : scattering phase function in Hapke's BDR model
c    THETA  : phase angle (radians); the angle between incidence and
c             reflection directions in Hapke's BDR model

C  Local variables

      INTEGER          J
      DOUBLE PRECISION CTHETA, THETA, PHASE
      DOUBLE PRECISION HOTSPOT, B0_EMPIR, HELP_HOT, B_HOT
      DOUBLE PRECISION SSALBEDO, GAMMA, REFLEC, FUNCTION
      DOUBLE PRECISION HELP_J, GHELP_J, TERM_J
      DOUBLE PRECISION HELP_I, GHELP_I, TERM_I
      DOUBLE PRECISION TI_TJ, DT1, DT2
      DOUBLE PRECISION XPHI, CKPHI

C  Initialise

      HAPKE_KERNEL       = ZERO
      DO J = 1, NPARS
        HAPKE_DERIVATIVES(J) = ZERO
      ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

C  geometrical part

C  This is the code that is in DISORT - not right, I think.
C        CTHETA = XI * XJ + DABS(SXI) *  DABS(SXJ) * CKPHI

      CTHETA = XI * XJ + SXI * SXJ * CKPHI
      THETA  = DACOS( CTHETA )
      PHASE  = ONE + HALF * CTHETA

C  hot spot parameterization

      HOTSPOT  = PARS(2)
      B0_EMPIR = PARS(3)
      HELP_HOT = HOTSPOT + DTAN ( HALF * THETA )
      B_HOT    = B0_EMPIR * HOTSPOT / HELP_HOT

C  Albedo parameterization

      SSALBEDO = PARS(1)
      GAMMA    = DSQRT ( ONE - SSALBEDO )
      HELP_J   = TWO * XJ
      GHELP_J  = ( ONE + HELP_J * GAMMA )
      TERM_J   = ( ONE + HELP_J ) / GHELP_J
      HELP_I   = TWO * XI
      GHELP_I  = ( ONE + HELP_I * GAMMA )
      TERM_I   = ( ONE + HELP_I ) / GHELP_I
      TI_TJ    = TERM_J * TERM_I

C  Function

      REFLEC       = SSALBEDO * QUARTER / ( XI + XJ )
      FUNCTION     = ( ONE + B_HOT ) * PHASE + TI_TJ - ONE
      HAPKE_KERNEL = REFLEC * FUNCTION

C  ssalbedo derivative

      IF ( DO_DERIV_PARS(1) ) THEN
        DT1 = HAPKE_KERNEL / SSALBEDO
        DT2 = ( HELP_J / GHELP_J ) + ( HELP_I / GHELP_I )
        DT2 = DT2 * TI_TJ * HALF / GAMMA
        HAPKE_DERIVATIVES(1) = DT1 + DT2 * REFLEC
      ENDIF

C  Hotspot  derivative

      IF ( DO_DERIV_PARS(2) ) THEN
        DT1 = B_HOT * ( B0_EMPIR - B_HOT ) / B0_EMPIR / HOTSPOT
        HAPKE_DERIVATIVES(2) = DT1 * REFLEC * PHASE
      ENDIF

C  empirical factor derivative

      IF ( DO_DERIV_PARS(3) ) THEN
        DT1 = B_HOT / B0_EMPIR 
        HAPKE_DERIVATIVES(3) = DT1 * REFLEC * PHASE
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE LISPARSE_FUNCTION_PLUS 
     I       ( MAXPARS, NPARS, PARS, DO_DERIV_PARS,
     I         XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,
     O         LISPARSE_KERNEL, LISPARSE_DERIVATIVES )

C  include file of constants amd dimensions

      INCLUDE '../includes/LIDORT.PARS'

C  Subroutine arguments

      INTEGER          MAXPARS, NPARS
      DOUBLE PRECISION PARS ( MAXPARS )
      LOGICAL          DO_DERIV_PARS ( MAXPARS )
      DOUBLE PRECISION XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION LISPARSE_KERNEL
      DOUBLE PRECISION LISPARSE_DERIVATIVES ( MAXPARS )

C  local variables

      INTEGER          J
      DOUBLE PRECISION X_INC, X_REF, SX_INC, SX_REF, ANG_P, TX
      DOUBLE PRECISION T_INC, T_REF, T_INC_SQ, T_REF_SQ
      DOUBLE PRECISION CKSI, DELTA, T, COST, SINT, DSQ, SINTCOST
      DOUBLE PRECISION A, B, H, R, P, Q, DT1, DT2, DT2SQ, QR
      DOUBLE PRECISION A2, R2, DX_H, DX_R, DX_Q, DX_P, DX_QR, DY_Q
      DOUBLE PRECISION XPHI, CKPHI

C  Initialise
C    -- Return for special case

      LISPARSE_KERNEL       = ZERO
      DO J = 1, NPARS
        LISPARSE_DERIVATIVES(J) = ZERO
      ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI
      IF ( ( XI .EQ. XJ ) .AND. ( CKPHI.EQ.ONE ) ) RETURN

C  Function
C  ========

C  .. incidence

      TX       = SXJ / XJ
      T_INC    = PARS(2) * TX
      T_INC_SQ = T_INC * T_INC
      ANG_P    = DATAN ( T_INC )
      X_INC    = DCOS(ANG_P)
      SX_INC   = DSIN(ANG_P)

C  .. reflection

      TX       = SXI / XI
      T_REF    = PARS(2) * TX
      T_REF_SQ = T_REF * T_REF
      ANG_P    = DATAN ( T_REF )
      X_REF    = DCOS(ANG_P)
      SX_REF   = DSIN(ANG_P)

C  ksi cosine

      CKSI = X_INC  * X_REF + SX_INC * SX_REF * CKPHI

C  contributions P and R and derivatives (if flagged)

      P = ( ONE + CKSI ) / X_REF
      A = ( ONE / X_INC )
      B = ( ONE / X_REF )
      R = A + B

      DX_R = ZERO
      DX_P = ZERO
      IF ( DO_DERIV_PARS(2) ) THEN
        DX_R = R * ( ONE - ( ONE / A / B ) ) / PARS(2)
        R2   = R * R
        A2   = A * A
        DX_P = ( P * ( ONE + A2 ) - ( R2 / B ) ) / PARS(2) / A2
      ENDIF

C  evaluate cos(t)

      DT1   = T_REF_SQ + T_INC_SQ
      DT2   = T_INC * T_REF
      DT2SQ = DT2 * DT2
      DELTA = DSQRT ( DT1 - TWO * DT2 * CKPHI )
      DSQ   = DELTA * DELTA
      H     = DSQRT ( DSQ + SKPHI * SKPHI * DT2SQ )
      COST  = PARS(1) * H / R

C  set Q function and its derivatives if flagged

      DX_Q = ZERO
      DY_Q = ZERO
      IF ( COST .GT. ONE ) THEN
        Q = ONE
      ELSE
        T        = DACOS(COST)
        SINT     = DSQRT ( ONE - COST * COST )
        SINTCOST = SINT * COST
        Q = ONE -  ( ( T - SINTCOST ) / PIE )
        IF ( DO_DERIV_PARS(2) ) THEN
          DX_H = ( TWO * H * H - DSQ ) / H / PARS(2)
          DX_Q = TWO * SINTCOST * ( (DX_H/H) - (DX_R/R) ) / PIE
        ENDIF
        IF ( DO_DERIV_PARS(1) ) THEN
          DY_Q = TWO * SINTCOST / PIE / PARS(1)
        ENDIF
      ENDIF

C  set the kernel
C  --------------

      QR = Q * R 
      LISPARSE_KERNEL = HALF * P - QR

C  Set derivatives
C  ---------------

      IF ( DO_DERIV_PARS(1) ) THEN
        LISPARSE_DERIVATIVES(1) = - R * DY_Q
      ENDIF

      IF ( DO_DERIV_PARS(2) ) THEN
        DX_QR = DX_R * Q + DX_Q * R
        LISPARSE_DERIVATIVES(2) = HALF * DX_P - DX_QR
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDENSE_FUNCTION_PLUS 
     I       ( MAXPARS, NPARS, PARS, DO_DERIV_PARS,
     I         XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,
     O         LIDENSE_KERNEL, LIDENSE_DERIVATIVES )

C  include file of constants amd dimensions

      INCLUDE '../includes/LIDORT.PARS'

C  Subroutine arguments

      INTEGER          MAXPARS, NPARS
      DOUBLE PRECISION PARS ( MAXPARS )
      LOGICAL          DO_DERIV_PARS ( MAXPARS )
      DOUBLE PRECISION XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION LIDENSE_KERNEL
      DOUBLE PRECISION LIDENSE_DERIVATIVES ( MAXPARS )

C  local variables

      INTEGER          J
      DOUBLE PRECISION X_INC, X_REF, SX_INC, SX_REF, ANG_P, TX
      DOUBLE PRECISION T_INC, T_REF, T_INC_SQ, T_REF_SQ
      DOUBLE PRECISION CKSI, DELTA, T, COST, SINT, DSQ, SINTCOST
      DOUBLE PRECISION A, B, H, R, P, Q, DT1, DT2, DT2SQ, P_QR
      DOUBLE PRECISION A2, R2, DX_H, DX_R, DX_Q, DX_P, DX_P_QR, DY_Q
      DOUBLE PRECISION XPHI, CKPHI

C  Initialise
C    -- Return for special case

      LIDENSE_KERNEL       = ZERO
      DO J = 1, NPARS
        LIDENSE_DERIVATIVES(J) = ZERO
      ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI
      IF ( ( XI .EQ. XJ ) .AND. ( CKPHI.EQ.ONE ) ) RETURN

C  Function
C  ========

C  .. incidence

      TX       = SXJ / XJ
      T_INC    = PARS(2) * TX
      T_INC_SQ = T_INC * T_INC
      ANG_P    = DATAN ( T_INC )
      X_INC    = DCOS(ANG_P)
      SX_INC   = DSIN(ANG_P)

C  .. reflection

      TX       = SXI / XI
      T_REF    = PARS(2) * TX
      T_REF_SQ = T_REF * T_REF
      ANG_P    = DATAN ( T_REF )
      X_REF    = DCOS(ANG_P)
      SX_REF   = DSIN(ANG_P)

C  ksi cosine

      CKSI = X_INC  * X_REF + SX_INC * SX_REF * CKPHI

C  contributions P and R and derivatives (if flagged)

      P = ( ONE + CKSI ) / X_REF
      A = ( ONE / X_INC )
      B = ( ONE / X_REF )
      R = A + B

      DX_R = ZERO
      DX_P = ZERO
      IF ( DO_DERIV_PARS(2) ) THEN
        DX_R = R * ( ONE - ( ONE / A / B ) ) / PARS(2)
        R2   = R * R
        A2   = A * A
        DX_P = ( P * ( ONE + A2 ) - ( R2 / B ) ) / PARS(2) / A2
      ENDIF

C  evaluate cos(t)

      DT1   = T_REF_SQ + T_INC_SQ
      DT2   = T_INC * T_REF
      DT2SQ = DT2 * DT2
      DELTA = DSQRT ( DT1 - TWO * DT2 * CKPHI )
      DSQ   = DELTA * DELTA
      H     = DSQRT ( DSQ + SKPHI * SKPHI * DT2SQ )
      COST  = PARS(1) * H / R

C  set Q function and its derivatives if flagged

      DX_Q = ZERO
      DY_Q = ZERO
      IF ( COST .GT. ONE ) THEN
        Q = ONE
      ELSE
        T        = DACOS(COST)
        SINT     = DSQRT ( ONE - COST * COST )
        SINTCOST = SINT * COST
        Q = ONE -  ( ( T - SINTCOST ) / PIE )
        IF ( DO_DERIV_PARS(2) ) THEN
          DX_H = ( TWO * H * H - DSQ ) / H / PARS(2)
          DX_Q = TWO * SINTCOST * ( (DX_H/H) - (DX_R/R) ) / PIE
        ENDIF
        IF ( DO_DERIV_PARS(1) ) THEN
          DY_Q = TWO * SINTCOST / PIE / PARS(1)
        ENDIF
      ENDIF

C  set the kernel
C  --------------

      P_QR = P / Q / R 
      LIDENSE_KERNEL = P_QR - TWO

C  Set derivatives
C  ---------------

      IF ( DO_DERIV_PARS(1) ) THEN
        LIDENSE_DERIVATIVES(1) = - P_QR * DY_Q / Q 
      ENDIF

      IF ( DO_DERIV_PARS(2) ) THEN
        DX_P_QR = ( DX_P / P ) - ( DX_R / R ) - ( DX_Q / Q )
        LIDENSE_DERIVATIVES(2) = P_QR * DX_P_QR
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE RAHMAN_FUNCTION_PLUS 
     I       ( MAXPARS, NPARS, PARS, DO_DERIV_PARS,
     I         XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,
     O         RAHMAN_KERNEL, RAHMAN_DERIVATIVES )

C  include file of constants amd dimensions

      INCLUDE '../includes/LIDORT.PARS'

C  Subroutine arguments

      INTEGER          MAXPARS, NPARS
      DOUBLE PRECISION PARS ( MAXPARS )
      LOGICAL          DO_DERIV_PARS ( MAXPARS )
      DOUBLE PRECISION XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION RAHMAN_KERNEL
      DOUBLE PRECISION RAHMAN_DERIVATIVES ( MAXPARS )

C  local variables

      INTEGER          J
      DOUBLE PRECISION T_INC, T_REF, DT1, DT2 
      DOUBLE PRECISION CXI, DELTA, K1_SQ, FACT, K0, K1, K2
      DOUBLE PRECISION HELPM, HELPR, HELPG, D_HELPM, D_FACT
      DOUBLE PRECISION GEOM, PHASE, RFAC, RFAC1, D_K0, D_K1, D_K2
      DOUBLE PRECISION XPHI, CKPHI

C  Initialise

      RAHMAN_KERNEL = ZERO
      DO J = 1, NPARS
        RAHMAN_DERIVATIVES(J) = ZERO
      ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI
      IF ( XI.EQ.ZERO .OR. XJ.EQ.ZERO ) RETURN

C  parameters

      K0 = PARS(1)
      K1 = PARS(2)
      K2 = PARS(3)

C  geometrical angle xi

        CXI = XI * XJ + SXI * SXJ * CKPHI
      IF ( CXI .GT. ONE ) CXI = ONE

C  Phase function

      K1_SQ = K1 * K1
      HELPM = ONE - K1_SQ 
      FACT  = ONE + K1_SQ + TWO * K1 * CXI
      PHASE = HELPM / ( FACT ** ONEP5 )

C  Delta and R-factor

      T_INC = SXI / XI
      T_REF = SXJ / XJ
      DT1   = T_INC*T_INC + T_REF*T_REF
      DT2   = T_INC * T_REF
      DELTA = DSQRT ( DT1 - TWO * DT2 * CKPHI )
      HELPR = ONE / ( ONE + DELTA )
      RFAC  = ( ONE - K0 ) * HELPR
      RFAC1 = ONE + RFAC

C  Geom factor and kernel

      HELPG = XI * XJ * ( XI + XJ )
      GEOM  = HELPG ** ( K2 - ONE)
      RAHMAN_KERNEL = K0 * PHASE * RFAC1 * GEOM

C  K0 derivative

      IF ( DO_DERIV_PARS(1) ) THEN
        D_K0   = ( ONE / K0 ) - ( HELPR / RFAC1 )
        RAHMAN_DERIVATIVES(1) = RAHMAN_KERNEL * D_K0
      ENDIF

C  Phase function derivative

      IF ( DO_DERIV_PARS(2) ) THEN
        D_FACT  =   TWO * K1 + TWO * CXI
        D_HELPM = - TWO * K1
        D_K1    = ( D_HELPM / HELPM ) - ONEP5 * ( D_FACT / FACT )
        RAHMAN_DERIVATIVES(2) = RAHMAN_KERNEL * D_K1
      ENDIF

C  K2 derivative

      IF ( DO_DERIV_PARS(3) ) THEN
        D_K2 = DLOG ( HELPG )
        RAHMAN_DERIVATIVES(3) = RAHMAN_KERNEL * D_K2
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE COXMUNK_FUNCTION_PLUS 
     I       ( MAXPARS, NPARS, PARS, DO_DERIV_PARS,
     I         XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,
     O         COXMUNK_KERNEL, COXMUNK_DERIVATIVES )

C  include file of constants amd dimensions

      INCLUDE '../includes/LIDORT.PARS'

C  Subroutine arguments

      INTEGER          MAXPARS, NPARS
      DOUBLE PRECISION PARS ( MAXPARS )
      LOGICAL          DO_DERIV_PARS ( MAXPARS )
      DOUBLE PRECISION XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION COXMUNK_KERNEL
      DOUBLE PRECISION COXMUNK_DERIVATIVES ( MAXPARS )

C  Critical exponent taken out

      DOUBLE PRECISION CRITEXP
      PARAMETER       ( CRITEXP = 88.0D0 )

C  Local variables

      INTEGER          J
      DOUBLE PRECISION Z, Z1, Z2, Z2_SQ_M1, H1, H2, RP, RL, XMP
      DOUBLE PRECISION A, B, TA, ARGUMENT, PROB, FAC1, FAC2
      DOUBLE PRECISION S1, S2, S3, XXI, XXJ
      DOUBLE PRECISION T1_I, T2_I, DCOT_I
      DOUBLE PRECISION T1_R, T2_R, DCOT_R
      DOUBLE PRECISION SHADOWI, SHADOWR, SHADOW
      DOUBLE PRECISION XPHI, CKPHI

      DOUBLE PRECISION H1H2, H2Z2, TA_SQ, DFAC2, DH1, DH2, DRP, DRL
      DOUBLE PRECISION D_S1, D_S2, D_T1, D_T2
      DOUBLE PRECISION D_SHADOWI, D_SHADOWR, D_SHADOW

C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C               Remark on Use of shadow effect
C               ------------------------------
C  Shadow effect is controlled by the third parameter. That is, if
C  PARS(3) not equal to then shadow effect will be included.
C    --- NPARS should always be 3 for this Kernel.
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

C  Initialise

      COXMUNK_KERNEL     = ZERO
      DO J = 1, NPARS
        COXMUNK_DERIVATIVES(J) = ZERO
      ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

C  Kernel

C  ..Scatter angles

C old   Z = - XI * XJ + SXI * SXJ * CKPHI   
C old   IF ( Z .LT. MINUS_ONE) Z = MINUS_ONE
C old   Z1 = DACOS(-Z)
C old   Z2 = DCOS(Z1*HALF)

      Z = XI * XJ + SXI * SXJ * CKPHI   
      IF ( Z .GT. ONE) Z = ONE
      Z1 = DACOS(Z)
      Z2 = DCOS(Z1*HALF)

C  .. Fresnel coefficients

      Z2_SQ_M1 = Z2 * Z2 + MINUS_ONE
      H1 = PARS(2) * Z2
      H2 = DSQRT ( PARS(2) + Z2_SQ_M1 )
      H1H2 = H1 + H2
      RP = ( H1 - H2 ) / H1H2
      H2Z2 = Z2 + H2
      RL = ( Z2 - H2 ) / H2Z2
      XMP = HALF * ( RP*RP + RL*RL )

C  Coxmunk Function

      A = TWO * Z2                 
      B = ( XI + XJ ) / A                  
      IF ( B .GT. ONE ) B = ONE           
      A  = PIO2 - DASIN(B)
      TA = DTAN(A)
      TA_SQ = TA * TA            
      ARGUMENT = TA_SQ  / PARS(1)
      IF ( ARGUMENT .LT. CRITEXP ) THEN
        PROB = DEXP ( - ARGUMENT )
        FAC1 = PROB / PARS(1)
        FAC2 = QUARTER / XI / ( B ** FOUR )
        COXMUNK_KERNEL = XMP * FAC1 * FAC2 / XJ
      ENDIF

C  inverse slope-squared derivative

      IF ( DO_DERIV_PARS(1) ) THEN
        IF ( ARGUMENT .LT. CRITEXP ) THEN
          DFAC2 = ( PARS(1) - TA_SQ ) / PARS(1) / PARS(1)
          COXMUNK_DERIVATIVES(1) = - COXMUNK_KERNEL * DFAC2
        ENDIF
      ENDIF

C  square refractive index derivative
C  --This section of code was formerly at the end of routine
C    -- Now moved here before the shadowing option
C    -- otherwise derivative will not get done
C         Bug found by V. Natraj in VLIDORT. 02 February 2007.

      IF ( DO_DERIV_PARS(2) ) THEN
        IF ( ARGUMENT .LT. CRITEXP ) THEN
          DH1 = Z2
          DH2 = HALF / H2
          DRP = ( DH1 * ( ONE - RP ) - DH2 * ( ONE + RP ) ) / H1H2
          DRL =  - DH2 * ( ONE + RL ) / H2Z2
          DFAC2 = ( RP*DRP + RL*DRL ) / XMP
          COXMUNK_DERIVATIVES(2) = COXMUNK_KERNEL * DFAC2
        ENDIF
      ENDIF

C  No Shadow code if not flagged

      IF ( PARS(3) .EQ. ZERO ) RETURN

C  Shadow code
C  -----------

      S1 = DSQRT(PARS(1)/PIE)
      S3 = ONE/(DSQRT(PARS(1)))
      S2 = S3*S3

      XXI  = XI*XI
      DCOT_I = XI/DSQRT(ONE-XXI)
      T1_I   = DEXP(-DCOT_I*DCOT_I*S2)
      T2_I   = DERFC(DCOT_I*S3)
      SHADOWI = HALF * ( S1*T1_I/DCOT_I - T2_I )

      XXJ  = XJ*XJ
      DCOT_R = XJ/DSQRT(ONE-XXJ)
      T1_R   = DEXP(-DCOT_R*DCOT_R*S2)
      T2_R   = DERFC(DCOT_R*S3)
      SHADOWR = HALF * ( S1*T1_R/DCOT_R - T2_R )

      SHADOW = ONE/(ONE+SHADOWI+SHADOWR)
      COXMUNK_KERNEL = COXMUNK_KERNEL * SHADOW

C  Update Scalar derivatives
C  -------------------------

C  add the shadow derivative to inverse slope-squared derivative

      IF ( DO_DERIV_PARS(1) ) THEN
        D_S1 = HALF / PIE / S1
        D_S2 = - S2 * S2
        D_T1 = - T1_I * DCOT_I * DCOT_I * D_S2
        D_T2 = S2 * S2 * DCOT_I * S1 * T1_I 
        D_SHADOWI = HALF * ( D_S1*T1_I/DCOT_I + S1*D_T1/DCOT_I - D_T2 )
        D_T1 = - T1_R * DCOT_R * DCOT_R * D_S2
        D_T2 = S2 * S2 * DCOT_R * S1 * T1_R
        D_SHADOWR = HALF * ( D_S1*T1_R/DCOT_R + S1*D_T1/DCOT_R - D_T2 )
        D_SHADOW = - SHADOW * SHADOW * ( D_SHADOWI + D_SHADOWR )
        COXMUNK_DERIVATIVES(1) =
     &         COXMUNK_DERIVATIVES(1) * SHADOW+
     *         COXMUNK_KERNEL * D_SHADOW / SHADOW
      ENDIF

C  Refractive index derivative, update

      IF ( DO_DERIV_PARS(2) ) THEN
        COXMUNK_DERIVATIVES(2) = COXMUNK_DERIVATIVES(2) * SHADOW 
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE COXMUNK_FUNCTION_PLUS_DB
     I       ( MAXPARS, NPARS, PARS, DO_DERIV_PARS,
     I         XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,
     O         COXMUNK_KERNEL, COXMUNK_DERIVATIVES )

C  include file of constants

      INCLUDE '../includes/LIDORT.PARS'

C  Subroutine arguments

      INTEGER          MAXPARS, NPARS
      DOUBLE PRECISION PARS ( MAXPARS )
      LOGICAL          DO_DERIV_PARS ( MAXPARS )
      DOUBLE PRECISION XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION COXMUNK_KERNEL
      DOUBLE PRECISION COXMUNK_DERIVATIVES ( MAXPARS )

C  local variables
C  ---------------

C  help variables

      integer          n, k, i, i1, q, N_phiquad_HALF
      DOUBLE PRECISION XM, SXM, sum_pr, pr, dfunc(2)
      DOUBLE PRECISION sumr, w_p, reflec_0, reflec_1
      double precision d_reflec_0(3), d_reflec_1(3)
      double precision phi_sub1, cphi_sub1, sphi_sub1
      double precision phi_sub2, cphi_sub2, sphi_sub2

C  Local quadrature stuff

      integer          max_msrs_muquad, max_msrs_phiquad
      integer          n_muquad, n_phiquad
      parameter       ( max_msrs_muquad = 40 )
      parameter       ( max_msrs_phiquad = 100 )

C  arrays

      DOUBLE PRECISION X_MUQUAD (max_msrs_muquad)
      DOUBLE PRECISION W_MUQUAD (max_msrs_muquad)
      DOUBLE PRECISION SX_MUQUAD (max_msrs_muquad)
      DOUBLE PRECISION WXX_MUQUAD(max_msrs_muquad)

      DOUBLE PRECISION X_PHIQUAD (max_msrs_phiquad)
      DOUBLE PRECISION W_PHIQUAD (max_msrs_phiquad)

      DOUBLE PRECISION R0_QUAD_IN  (max_msrs_muquad,max_msrs_phiquad)
      DOUBLE PRECISION R0_OUT_QUAD (max_msrs_muquad,max_msrs_phiquad)
      DOUBLE PRECISION
     &       D_R0_QUAD_IN  (max_msrs_muquad,max_msrs_phiquad,2),
     &       D_R0_OUT_QUAD (max_msrs_muquad,max_msrs_phiquad,2)

C  Safety first zeroing

      REFLEC_0 = ZERO
      REFLEC_1 = ZERO
      COXMUNK_KERNEL = ZERO

      DO Q = 1, NPARS
        IF ( DO_DERIV_PARS(Q) ) THEN
          COXMUNK_DERIVATIVES(Q) = ZERO
          D_REFLEC_0(Q) = ZERO
          D_REFLEC_1(Q) = ZERO
        ENDIF

      ENDDO

C  Air to water, Polar quadrature

      n_muquad = 40
      CALL GAULEG ( ZERO, ONE, X_muquad, W_muquad, n_muquad )
      DO I = 1, N_MUQUAD
        XM = X_MUQUAD(I)
        SX_MUQUAD(I) = DSQRT(ONE-XM*XM)
        WXX_MUQUAD(I) = XM * XM * W_MUQUAD(I)
      ENDDO

C  Azimuth quadrature

      n_phiquad = 100
      N_phiquad_HALF = N_PHIQUAD / 2
      CALL GAULEG ( ZERO, ONE, X_PHIQUAD, W_PHIQUAD, N_PHIQUAD_HALF )
      DO I = 1, N_PHIQUAD_HALF
        I1 = I + N_PHIQUAD_HALF
        X_PHIQUAD(I1) = - X_PHIQUAD(I)
        W_PHIQUAD(I1) =   W_PHIQUAD(I)
      ENDDO
      DO I = 1, N_PHIQUAD
        X_PHIQUAD(I)  = PIE * X_PHIQUAD(I)
      ENDDO

C  Single scattering (zero order), Phi is in degrees here!

      CALL COXMUNK_FUNCTION_PLUS
     I       ( MAXPARS, NPARS, PARS, DO_DERIV_PARS,
     I         XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,
     O         REFLEC_0, D_REFLEC_0 )

C  Quadrature output for first order R/T calculations 

      DO K = 1, n_muquad
        XM  = X_MUQUAD(K)
        SXM = SX_MUQUAD(K)
        DO N = 1, N_PHIQUAD
          PHI_SUB1  = X_PHIQUAD(N)
          CPHI_SUB1 = DCOS(PHI_SUB1)
          SPHI_SUB1 = DSIN(PHI_SUB1)
          PHI_SUB2  = PHI*DEG_TO_RAD - X_PHIQUAD(N)
          CPHI_SUB2 = DCOS(PHI_SUB2)
          SPHI_SUB2 = DSIN(PHI_SUB2)
          CALL COXMUNK_FUNCTION_PLUS
     I       ( MAXPARS, NPARS, PARS, DO_DERIV_PARS,
     I         XM, SXM, XI, SXI, PHI_SUB2, CPHI_SUB2, SPHI_SUB2,
     O         R0_OUT_QUAD(K,N), DFUNC )
          D_R0_OUT_QUAD(K,N,1) = DFUNC(1)
          D_R0_OUT_QUAD(K,N,2) = DFUNC(2)
          CALL COXMUNK_FUNCTION_PLUS
     I       ( MAXPARS, NPARS, PARS, DO_DERIV_PARS,
     I         XJ, SXJ, XM, SXM, PHI_SUB1, CPHI_SUB1, SPHI_SUB1,
     O         R0_QUAD_IN(K,N),  DFUNC )
          D_R0_QUAD_IN(K,N,1) = DFUNC(1)
          D_R0_QUAD_IN(K,N,2) = DFUNC(2)
        ENDDO
      ENDDO

C  compute the next order

      SUMR = ZERO
      DO K = 1, n_muquad
        SUM_PR = ZERO
        DO N = 1, N_PHIQUAD
          W_P  = W_PHIQUAD(N)
          SUM_PR = SUM_PR + W_P * R0_QUAD_IN(K,N) * R0_OUT_QUAD(K,N)
        ENDDO
        SUMR = SUMR + SUM_PR * WXX_MUQUAD(K)
      ENDDO
      REFLEC_1 = SUMR

C  Compute total

      COXMUNK_KERNEL = REFLEC_0 + REFLEC_1

C  Derivatives

      DO Q = 1, 2
       IF ( DO_DERIV_PARS(Q) ) THEN
        SUMR = ZERO
        DO K = 1, n_muquad
         PR = ZERO
         DO N = 1, N_PHIQUAD
          W_P  = W_PHIQUAD(N)
          PR = PR + W_P *   R0_QUAD_IN(K,N)   * D_R0_OUT_QUAD(K,N,Q)
     &            + W_P * D_R0_QUAD_IN(K,N,Q) *   R0_OUT_QUAD(K,N)
         ENDDO
         SUMR = SUMR + PR * WXX_MUQUAD(K)
        ENDDO
        D_REFLEC_1(Q) = SUMR
        COXMUNK_DERIVATIVES(Q) = D_REFLEC_0(Q) + D_REFLEC_1(Q)
       ENDIF
      ENDDO

c      write(34,'(1p6e14.5)')reflec_0, reflec_1, pars(1),
c     &   dacos(xi)/deg_to_rad, dacos(xj)/deg_to_rad, phi

C  Finish

      RETURN
      END
