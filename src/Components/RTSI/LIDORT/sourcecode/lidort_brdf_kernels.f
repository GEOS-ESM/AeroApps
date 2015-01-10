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
C #            ROSSTHIN_FUNCTION                                #
C #            ROSSTHICK_FUNCTION                               #
C #            LISPARSE_FUNCTION                                #
C #            LIDENSE_FUNCTION                                 #
C #            ROUJEAN_FUNCTION                                 #
C #            HAPKE_FUNCTION                                   #
C #            RAHMAN_FUNCTION                                  #
C #            COXMUNK_FUNCTION                                 #
C #            COXMUNK_FUNCTION_DB                              #
C #                                                             #
C ###############################################################

      SUBROUTINE ROSSTHIN_FUNCTION 
     I       ( MAXPARS, NPARS, PARS,
     I         XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,
     O         ROSSTHIN_KERNEL )

C  include file of constants

      INCLUDE '../includes/LIDORT.PARS'

C  Subroutine arguments

      INTEGER          MAXPARS, NPARS
      DOUBLE PRECISION PARS ( MAXPARS )
      DOUBLE PRECISION XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION ROSSTHIN_KERNEL

C  Local variables

      DOUBLE PRECISION DS1, DS2, CKSI, SKSI, KSI, FUNC
      DOUBLE PRECISION XPHI, CKPHI

C  Initialise

      ROSSTHIN_KERNEL = ZERO
      XPHI = PIE - PHI
      CKPHI = - CPHI

C  kernel

      DS1 = XI * XJ
      DS2 = SXI * SXJ
      CKSI = DS1 + DS2 * CKPHI
      IF ( CKSI.GT.ONE ) CKSI = ONE
      SKSI = DSQRT(ONE-CKSI*CKSI)
      KSI = DACOS(CKSI)
      FUNC = ((PIO2-KSI)*CKSI + SKSI)/DS1
      ROSSTHIN_KERNEL = FUNC - PIO2

C  Finish

      RETURN
      END

C

      SUBROUTINE ROSSTHICK_FUNCTION 
     I       ( MAXPARS, NPARS, PARS,
     I         XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,
     O         ROSSTHICK_KERNEL )

C  include file of constants

      INCLUDE '../includes/LIDORT.PARS'

C  Subroutine arguments

      INTEGER          MAXPARS, NPARS
      DOUBLE PRECISION PARS ( MAXPARS )
      DOUBLE PRECISION XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION ROSSTHICK_KERNEL

C  Local variables

      DOUBLE PRECISION DS1, DS2, DS3, CKSI, SKSI, KSI, FUNC
      DOUBLE PRECISION XPHI, CKPHI

C  Initialise

      ROSSTHICK_KERNEL = ZERO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

C  kernel

      DS1 = XI * XJ
      DS2 = SXI * SXJ
      DS3 = XI  + XJ
      CKSI = DS1 + DS2 * CKPHI
      IF ( CKSI.GT.ONE ) CKSI = ONE
      SKSI = DSQRT(ONE-CKSI*CKSI)
      KSI = DACOS(CKSI)
      FUNC = ((PIO2-KSI)*CKSI + SKSI)/DS3
      ROSSTHICK_KERNEL = FUNC - PIO4

C  Finish

      RETURN
      END

C

      SUBROUTINE ROUJEAN_FUNCTION 
     I       ( MAXPARS, NPARS, PARS,
     I         XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,
     O         ROUJEAN_KERNEL )

C  include file of constants

      INCLUDE '../includes/LIDORT.PARS'

C  Subroutine arguments

      INTEGER          MAXPARS, NPARS
      DOUBLE PRECISION PARS ( MAXPARS )
      DOUBLE PRECISION XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION ROUJEAN_KERNEL

C  Local variables

      DOUBLE PRECISION DS1, DS2, DS3, TXJ, TXI, PHIFAC, S1, S2
      DOUBLE PRECISION XPHI_R, CXPHI_R, SXPHI_R, XPHI_C
      DOUBLE PRECISION XPHI, CKPHI

C  Initialise

      ROUJEAN_KERNEL = ZERO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

C  kernel

      XPHI_C = XPHI
      IF ( XPHI .GT. PIE )  XPHI_C = TWO*PIE - XPHI
      IF ( XPHI .LT. ZERO ) XPHI_C = - XPHI

      IF ( SXI .LT. ZERO ) THEN
        XPHI_R  = ( PIE - XPHI_C )
        CXPHI_R = DCOS ( XPHI_R )
        SXPHI_R = DSIN ( XPHI_R )
        TXI =  - ( SXI / XI )
      ELSE
        TXI =   ( SXI / XI )
        XPHI_R  = XPHI_C
        CXPHI_R = DCOS ( XPHI_R )
        SXPHI_R = DSIN ( XPHI_R )
      ENDIF

      TXJ =  ( SXJ / XJ )
      DS1 = TWO * TXJ * TXI
      DS2 = TXJ + TXI
      DS3 = TXJ*TXJ  + TXI*TXI
      PHIFAC = ( ( PIE - XPHI_R ) * CXPHI_R + SXPHI_R ) / PI4
      S1 = PHIFAC * DS1
      S2 = ( DS2 + DSQRT ( DS3 - DS1 * CXPHI_R ) ) / PIE
      ROUJEAN_KERNEL = S1 - S2

C  Finish

      RETURN
      END

C

      SUBROUTINE LISPARSE_FUNCTION
     I       ( MAXPARS, NPARS, PARS,
     I         XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,
     O         LISPARSE_KERNEL )

C  include file of constants

      INCLUDE '../includes/LIDORT.PARS'

C  Subroutine arguments

      INTEGER          MAXPARS, NPARS
      DOUBLE PRECISION PARS ( MAXPARS )
      DOUBLE PRECISION XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION LISPARSE_KERNEL

C  local variables

      DOUBLE PRECISION X_INC, X_REF, SX_INC, SX_REF, ANG_P, TX
      DOUBLE PRECISION T_INC, T_REF, T_INC_SQ, T_REF_SQ
      DOUBLE PRECISION CKSI, DELTA, T, COST, SINT, DSQ, SINTCOST
      DOUBLE PRECISION A, B, H, R, P, Q, DT1, DT2, DT2SQ, QR
      DOUBLE PRECISION XPHI, CKPHI

C  Initialise
C    -- Return for special case

      LISPARSE_KERNEL = ZERO
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

C  contributions P and R

      P = ( ONE + CKSI ) / X_REF
      A = ( ONE / X_INC )
      B = ( ONE / X_REF )
      R = A + B

C  evaluate cos(t)

      DT1   = T_REF_SQ + T_INC_SQ
      DT2   = T_INC * T_REF
      DT2SQ = DT2 * DT2
      DELTA = DSQRT ( DT1 - TWO * DT2 * CKPHI )
      DSQ   = DELTA * DELTA
      H     = DSQRT ( DSQ + SKPHI * SKPHI * DT2SQ )
      COST  = PARS(1) * H / R

C  set Q function

      IF ( COST .GT. ONE ) THEN
        Q = ONE
      ELSE
        T        = DACOS(COST)
        SINT     = DSQRT ( ONE - COST * COST )
        SINTCOST = SINT * COST
        Q = ONE -  ( ( T - SINTCOST ) / PIE )
      ENDIF

C  set the kernel
C  --------------

      QR = Q * R 
      LISPARSE_KERNEL = HALF * P - QR

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDENSE_FUNCTION 
     I       ( MAXPARS, NPARS, PARS,
     I         XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,
     O         LIDENSE_KERNEL )

C  include file of constants

      INCLUDE '../includes/LIDORT.PARS'

C  Subroutine arguments

      INTEGER          MAXPARS, NPARS
      DOUBLE PRECISION PARS ( MAXPARS )
      DOUBLE PRECISION XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION LIDENSE_KERNEL

C  local variables

      DOUBLE PRECISION X_INC, X_REF, SX_INC, SX_REF, ANG_P, TX
      DOUBLE PRECISION T_INC, T_REF, T_INC_SQ, T_REF_SQ
      DOUBLE PRECISION CKSI, DELTA, T, COST, SINT, DSQ, SINTCOST
      DOUBLE PRECISION A, B, H, R, P, Q, DT1, DT2, DT2SQ, P_QR
      DOUBLE PRECISION XPHI, CKPHI

C  Initialise
C    -- Return for special case

      LIDENSE_KERNEL = ZERO
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

C  contributions P and R

      P = ( ONE + CKSI ) / X_REF
      A = ( ONE / X_INC )
      B = ( ONE / X_REF )
      R = A + B

C  evaluate cos(t)

      DT1   = T_REF_SQ + T_INC_SQ
      DT2   = T_INC * T_REF
      DT2SQ = DT2 * DT2
      DELTA = DSQRT ( DT1 - TWO * DT2 * CKPHI )
      DSQ   = DELTA * DELTA
      H     = DSQRT ( DSQ + SKPHI * SKPHI * DT2SQ )
      COST  = PARS(1) * H / R

C  set Q function

      IF ( COST .GT. ONE ) THEN
        Q = ONE
      ELSE
        T        = DACOS(COST)
        SINT     = DSQRT ( ONE - COST * COST )
        SINTCOST = SINT * COST
        Q = ONE -  ( ( T - SINTCOST ) / PIE )
      ENDIF

C  set the kernel
C  --------------

      P_QR = P / Q / R 
      LIDENSE_KERNEL = P_QR - TWO

C  Finish

      RETURN
      END

C

      SUBROUTINE HAPKE_FUNCTION 
     I       ( MAXPARS, NPARS, PARS,
     I         XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,
     O         HAPKE_KERNEL )

C  include file of constants

      INCLUDE '../includes/LIDORT.PARS'

C  Subroutine arguments

      INTEGER          MAXPARS, NPARS
      DOUBLE PRECISION PARS ( MAXPARS )
      DOUBLE PRECISION XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION HAPKE_KERNEL

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

C  local variables

      DOUBLE PRECISION CTHETA, THETA, PHASE
      DOUBLE PRECISION HOTSPOT, B0_EMPIR, HELP_HOT, B_HOT
      DOUBLE PRECISION SSALBEDO, GAMMA, REFLEC, FUNCTION
      DOUBLE PRECISION HELP_J, TERM_J, HELP_I, TERM_I
      DOUBLE PRECISION XPHI, CKPHI

C  Initialise

      HAPKE_KERNEL = ZERO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

C  kernel

C  geometrical part

C  This is the code that is in DISORT - not right, I think.
C       CTHETA = XI * XJ + DABS(SXI) *  DABS(SXJ) * CKPHI

      CTHETA = XI * XJ + SXI * SXJ * CKPHI
      IF ( CTHETA .GT. ONE ) CTHETA = ONE
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
      TERM_J   = ( ONE + HELP_J ) / ( ONE + HELP_J * GAMMA )
      HELP_I   = TWO * XI
      TERM_I   = ( ONE + HELP_I ) / ( ONE + HELP_I * GAMMA )

C  Function

      REFLEC       = SSALBEDO * QUARTER / ( XI + XJ )
      FUNCTION     = ( ONE + B_HOT ) * PHASE + TERM_J * TERM_I - ONE
      HAPKE_KERNEL = REFLEC * FUNCTION
 
C  Finish

      RETURN
      END

C

      SUBROUTINE RAHMAN_FUNCTION 
     I       ( MAXPARS,  NPARS, PARS,
     I         XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,
     O         RAHMAN_KERNEL )

C  Revision. 24 October 2007.
C  --------------------------

C    * Limiting cases and hotspot evaluation.
C    * Revision based on the DISORT_2 code
C    *  In Disort, this kernel is known as the RPV^ BRDF.

C     The RPV reference is:
C       Rahman, Pinty, Verstraete, 1993: Coupled Surface-Atmosphere 
C       Reflectance (CSAR) Model. 2. Semiempirical Surface Model Usable 
C       With NOAA Advanced Very High Resolution Radiometer Data,
C       J. Geophys. Res., 98, 20791-20801.

C  The hotspot should occur when XI = XJ and PHI = 180.

C  include file of constants

      INCLUDE '../includes/LIDORT.PARS'

C  Subroutine arguments

      INTEGER          MAXPARS, NPARS
      DOUBLE PRECISION PARS ( MAXPARS )
      DOUBLE PRECISION XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION RAHMAN_KERNEL

C  local variables

      DOUBLE PRECISION T_INC, T_REF, DT1, DT2 
      DOUBLE PRECISION CXI, DELTA, K1_SQ, FACT
      DOUBLE PRECISION GEOM, PHASE, RFAC, K0, K1, K2
      DOUBLE PRECISION XPHI, CKPHI, SMALL, HSPOT, UPPER_LIMIT
      PARAMETER        ( SMALL = 1.0d-04 )

C  Initial section
C  ---------------

C  Initialise output

      RAHMAN_KERNEL = ZERO

C  Limiting case, formerly
C      IF ( XI.EQ.ZERO .OR. XJ.EQ.ZERO ) RETURN

C  Limiting case, revised

      IF ( XJ.LT.SMALL ) RETURN

C  Azimuth convettion

      XPHI  = PIE - PHI
      CKPHI = - CPHI

C  parameters

      K0 = PARS(1)
      K1 = PARS(2)
      K2 = PARS(3)

C  Hot Spot
C  --------

C  Value of hot spot

      FACT = K0 * ( TWO - K0 )
      FACT = FACT * ( ONE - K1 ) / ( ONE + K1 ) / ( ONE + K1 )
      GEOM = ( TWO * XJ * XJ * XJ ) ** ( K2 - ONE )
      HSPOT = FACT * GEOM

C  Upper limit ( 5 times hotspot value ). Follwing comments inserted.
C     This function needs more checking; some constraints are 
C     required to avoid albedos larger than 1; in particular,
C     the BDREF is limited to 5 times the hotspot value to
C     avoid extremely large values at low polar angles

      UPPER_LIMIT = 5.0d0 * HSPOT

C  hot spot value

      IF ( DABS(PHI) .LT. SMALL .AND. XI.EQ.XJ ) THEN
        RAHMAN_KERNEL = HSPOT
        RETURN
      ENDIF

C  Use upper limit value at edges (low incidence or reflection)

      IF ( XI.LT.SMALL .OR. XJ.LT.SMALL ) THEN
        RAHMAN_KERNEL = UPPER_LIMIT
        RETURN
      ENDIF

C  Main section
C  ------------

C  geometrical angle xi

      CXI = XI * XJ + SXI * SXJ * CKPHI
      IF ( CXI .GT. ONE ) CXI = ONE

C  Phase function

      K1_SQ = K1 * K1
      FACT  = ( ONE + K1_SQ + TWO * K1 * CXI ) ** ONEP5
      PHASE = ( ONE - K1_SQ ) / FACT

C  Delta and R-factor

      T_INC = SXI / XI
      T_REF = SXJ / XJ
      DT1   = T_INC*T_INC + T_REF*T_REF
      DT2   = T_INC * T_REF
      DELTA = DSQRT ( DT1 - TWO * DT2 * CKPHI )
      RFAC = ( ONE - K0 ) / ( ONE + DELTA )

C  Geom factor and kernel

      GEOM = ( XI * XJ * ( XI + XJ ) ) ** ( K2 - ONE)
      RAHMAN_KERNEL = K0 * PHASE * ( ONE + RFAC ) * GEOM

C  Check upper limit not exceeded

      IF ( RAHMAN_KERNEL .GT. UPPER_LIMIT ) THEN
        RAHMAN_KERNEL = UPPER_LIMIT
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE HAPKE_FUNCTION_OLD 
     I       ( MAXPARS, NPARS, PARS,
     I         XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,
     O         HAPKE_KERNEL )

C  include file of constants

      INCLUDE '../includes/LIDORT.PARS'

C  Subroutine arguments

      INTEGER          MAXPARS, NPARS
      DOUBLE PRECISION PARS ( MAXPARS )
      DOUBLE PRECISION XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION HAPKE_KERNEL
      DOUBLE PRECISION XPHI, CKPHI

C  Hapke Kernel function.
C   From DISORT code; used as validation.

C  local variables

      REAL             MU, MUP, DUMMY, DPHI, HAPKEKER
      EXTERNAL         HAPKEKER

C  Initialise

      HAPKE_KERNEL = ZERO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

C  Kernel

      MUP = SNGL(XJ)
      MU = SNGL(XI)
      DPHI = SNGL(XPHI)
      DUMMY = ZERO
      HAPKE_KERNEL = HAPKEKER( DUMMY, DUMMY, MU, MUP, DPHI )

C  Finish

      RETURN
      END

C

      REAL FUNCTION  hapkeker( WVNMLO, WVNMHI, MU, MUP, DPHI )

c      Supplies surface bi-directional reflectivity.
c
c      NOTE 1: Bidirectional reflectivity in DISORT is defined
c              by Eq. 39 in STWL.
c      NOTE 2: Both MU and MU0 (cosines of reflection and incidence
c              angles) are positive.
c
c  INPUT:
c
c    WVNMLO : Lower wavenumber (inv cm) of spectral interval
c
c    WVNMHI : Upper wavenumber (inv cm) of spectral interval
c
c    MU     : Cosine of angle of reflection (positive)
c
c    MUP    : Cosine of angle of incidence (positive)
c
c    DPHI   : Difference of azimuth angles of incidence and reflection
c                (radians)
c
c  LOCAL VARIABLES:
c
c    IREF   : bidirectional reflectance options
c             1 - Hapke's BDR model
c
c    B0     : empirical factor to account for the finite size of
c             particles in Hapke's BDR model
c
c    B      : term that accounts for the opposition effect
c             (retroreflectance, hot spot) in Hapke's BDR model
c
c    CTHETA : cosine of phase angle in Hapke's BDR model
c
c    GAMMA  : albedo factor in Hapke's BDR model
c
c    H0     : H( mu0 ) in Hapke's BDR model
c
c    H      : H( mu ) in Hapke's BDR model
c
c    HH     : angular width parameter of opposition effect in Hapke's
c             BDR model
c
c    P      : scattering phase function in Hapke's BDR model
c
c    THETA  : phase angle (radians); the angle between incidence and
c             reflection directions in Hapke's BDR model
c
c    W      : single scattering albedo in Hapke's BDR model
c
c
c   Called by- DREF, SURFAC
c +-------------------------------------------------------------------+
c     .. Scalar Arguments ..

      REAL      DPHI, MU, MUP, WVNMHI, WVNMLO
c     ..
c     .. Local Scalars ..

      INTEGER   IREF
      REAL      B0, B, CTHETA, GAMMA, H0, H, HH, P, THETA, W
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC COS, SQRT
c     ..

      IREF = 1

      IF ( IREF.EQ.1 ) THEN

c                              ** Hapke's BRDF model (times Pi/Mu0)
c                              ** (Hapke, B., Theory of reflectance
c                              ** and emittance spectroscopy, Cambridge
c                              ** University Press, 1993, Eq. 8.89 on
c                              ** page 233. Parameters are from
c                              ** Fig. 8.15 on page 231, expect for w.)

         CTHETA = MU * MUP + (1.-MU**2)**.5 * (1.-MUP**2)**.5
     &            * COS( DPHI )
         THETA = ACOS( CTHETA )

         P    = 1. + 0.5 * CTHETA

         HH   = 0.06
         B0   = 1.0
         B    = B0 * HH / ( HH + TAN( THETA/2.) )

         W = 0.6
         GAMMA = SQRT( 1. - W )
         H0   = ( 1. + 2.*MUP ) / ( 1. + 2.*MUP * GAMMA )
         H    = ( 1. + 2.*MU ) / ( 1. + 2.*MU * GAMMA )

         hapkeker = W / 4. / (MU+MUP) * ( (1.+B)* P + H0 * H - 1.0 )

      END IF

      RETURN
      END

C

      SUBROUTINE COXMUNK_FUNCTION 
     I       ( MAXPARS,  NPARS, PARS,
     I         XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,
     O         COXMUNK_KERNEL )

C  include file of constants

      INCLUDE '../includes/LIDORT.PARS'

C  Subroutine arguments

      INTEGER          MAXPARS, NPARS
      DOUBLE PRECISION PARS ( MAXPARS )
      DOUBLE PRECISION XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION COXMUNK_KERNEL

C  Critical exponent taken out

      DOUBLE PRECISION CRITEXP
      PARAMETER       ( CRITEXP = 88.0D0 )

C  Local variables

      DOUBLE PRECISION Z, Z1, Z2, Z2_SQ_M1, H1, H2, RP, RL, XMP
      DOUBLE PRECISION A, B, TA, ARGUMENT, PROB, FAC1, FAC2
      DOUBLE PRECISION XPHI, CKPHI
      DOUBLE PRECISION S1, S2, S3, XXI, XXJ, T1, T2, DCOT
      DOUBLE PRECISION SHADOWI, SHADOWR, SHADOW

C  Shadow variables

C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C               Remark on Use of shadow effect
C               ------------------------------
C  Shadow effect is controlled by the third parameter. That is, if
C  PARS(3) not equal to then shadow effect will be included.
C    --- NPARS should always be 3 for this Kernel.
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

C  Initialise

      COXMUNK_KERNEL = ZERO
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
      RP = ( H1 - H2 ) / ( H1 + H2 )
      RL = ( Z2 - H2 ) / ( Z2 + H2 )
      XMP = HALF * ( RP*RP + RL*RL )

C  Coxmunk Function

      A = TWO * Z2
      B = ( XI + XJ ) / A
      IF ( B .GT. ONE ) B = ONE
      A = PIO2 - DASIN(B)
      TA = DTAN(A)
      ARGUMENT = TA * TA  / PARS(1)
      IF ( ARGUMENT .LT. CRITEXP ) THEN
        PROB = DEXP ( - ARGUMENT )
        FAC1 = PROB  / PARS(1)
        FAC2 = QUARTER / XI / ( B ** FOUR )
        COXMUNK_KERNEL = XMP * FAC1 * FAC2 / XJ
      ENDIF

C  No Shadow code if not flagged

      IF ( PARS(3) .EQ. ZERO ) RETURN

C  Shadow code

      S1 = DSQRT(PARS(1)/PIE)
      S3 = ONE/(DSQRT(PARS(1)))
      S2 = S3*S3

      XXI  = XI*XI
      DCOT = XI/DSQRT(ONE-XXI)
      T1   = DEXP(-DCOT*DCOT*S2)
      T2   = DERFC(DCOT*S3)
      SHADOWI = HALF*(S1*T1/DCOT-T2)

      XXJ  = XJ*XJ
      DCOT = XJ/DSQRT(ONE-XXJ)
      T1   = DEXP(-DCOT*DCOT*S2)
      T2   = DERFC(DCOT*S3)
      SHADOWR = HALF*(S1*T1/DCOT-T2)

      SHADOW = ONE/(ONE+SHADOWI+SHADOWR)
      COXMUNK_KERNEL = COXMUNK_KERNEL * SHADOW
 
C     Finish

      RETURN
      END

C

      SUBROUTINE COXMUNK_FUNCTION_DB
     I       ( MAXPARS, NPARS, PARS,
     I         XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,
     O         COXMUNK_KERNEL )

C  include file of constants

      INCLUDE '../includes/LIDORT.PARS'

C  Subroutine arguments

      INTEGER          MAXPARS, NPARS
      DOUBLE PRECISION PARS ( MAXPARS )
      DOUBLE PRECISION XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION COXMUNK_KERNEL

C  local variables
C  ---------------

C  help variables

      integer          n, k, i, i1, N_phiquad_HALF
      DOUBLE PRECISION XM, SXM, sum_pr
      DOUBLE PRECISION sumr, w_p, reflec_0, reflec_1
      double precision phi_sub1, cphi_sub1, sphi_sub1
      double precision phi_sub2, cphi_sub2, sphi_sub2

C  Local quadrature stuff

      integer          max_msrs_muquad, max_msrs_phiquad
      integer          n_muquad, n_phiquad
      parameter       ( max_msrs_muquad = 40 )
      parameter       ( max_msrs_phiquad = 50 )

C  arrays

      DOUBLE PRECISION X_MUQUAD (max_msrs_muquad)
      DOUBLE PRECISION W_MUQUAD (max_msrs_muquad)
      DOUBLE PRECISION SX_MUQUAD (max_msrs_muquad)
      DOUBLE PRECISION WXX_MUQUAD(max_msrs_muquad)

      DOUBLE PRECISION X_PHIQUAD (max_msrs_phiquad)
      DOUBLE PRECISION W_PHIQUAD (max_msrs_phiquad)

      DOUBLE PRECISION R0_QUAD_IN  (max_msrs_muquad,max_msrs_phiquad)
      DOUBLE PRECISION R0_OUT_QUAD (max_msrs_muquad,max_msrs_phiquad)

C  Safety first zeroing

      REFLEC_0 = ZERO
      REFLEC_1 = ZERO
      COXMUNK_KERNEL = ZERO

C  Air to water, Polar quadrature

      n_muquad = 20
      CALL GAULEG ( ZERO, ONE, X_muquad, W_muquad, n_muquad )
      DO I = 1, N_MUQUAD
        XM = X_MUQUAD(I)
        SX_MUQUAD(I) = DSQRT(ONE-XM*XM)
        WXX_MUQUAD(I) = XM * XM * W_MUQUAD(I)
      ENDDO

C  Azimuth quadrature

      n_phiquad = 40
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

      CALL COXMUNK_FUNCTION
     I       ( MAXPARS, NPARS, PARS,
     I         XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,
     O         REFLEC_0 )

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
          CALL COXMUNK_FUNCTION
     I       ( MAXPARS, NPARS, PARS,
     I         XM, SXM, XI, SXI, PHI_SUB2, CPHI_SUB2, SPHI_SUB2,
     O         R0_OUT_QUAD(K,N) )
          CALL COXMUNK_FUNCTION
     I       ( MAXPARS, NPARS, PARS,
     I         XJ, SXJ, XM, SXM, PHI_SUB1, CPHI_SUB1, SPHI_SUB1,
     O         R0_QUAD_IN(K,N) )
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

c      write(34,'(1p6e14.5)')reflec_0, reflec_1, coxmunk_kernel,
c     &   dacos(xi)/deg_to_rad, dacos(xj)/deg_to_rad, phi

C  Finish

      RETURN
      END
