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
C # Subroutines in this Module                                  #
C #                                                             #
C #            ROSSTHIN_VFUNCTION                               #
C #            ROSSTHICK_VFUNCTION                              #
C #            LISPARSE_VFUNCTION                               #
C #            LIDENSE_VFUNCTION                                #
C #            ROUJEAN_VFUNCTION                                #
C #            HAPKE_VFUNCTION                                  #
C #            RAHMAN_VFUNCTION                                 #
C #            COXMUNK_VFUNCTION                                #
C #            COXMUNK_FUNCTION_DB                              #
C #            GISSCOXMUNK_VFUNCTION                            #
C #            GISSCOXMUNK_FUNCTION_DB                          #
C #            RHERMAN_VFUNCTION                                #
C #            BREON_VFUNCTION                                  #
C #                                                             #
C #  new for Version 2.4R, introduced 30 April 2009, 6 May 2009 #
C #  2008 Veg/Soil functions based on Breon work 2008 for OCO   #
C #  2009 function is final Kernel supplied by Breon, 5/5/09    #
C #                                                             #
C #            BPDF2008VEG_VFUNCTION                            #
C #            BPDF2008SOIL_VFUNCTION                           #
C #            BPDF2009_VFUNCTION                               #
C #                                                             #
C ###############################################################

      SUBROUTINE ROSSTHIN_VFUNCTION 
     I       ( MAXPARS, N_BRDF_STOKESSQ, NPARS, PARS,
     I         XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,
     O         ROSSTHIN_VKERNEL )

C  include file of constants

      INCLUDE '../includes/VLIDORT.PARS'

C  Subroutine arguments

      INTEGER          MAXPARS, NPARS, N_BRDF_STOKESSQ
      DOUBLE PRECISION PARS ( MAXPARS )
      DOUBLE PRECISION XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION ROSSTHIN_VKERNEL(MAXSTOKES_SQ)

C  Local variables

      INTEGER          O1
      DOUBLE PRECISION DS1, DS2, CKSI, SKSI, KSI, FUNC
      DOUBLE PRECISION XPHI, CKPHI

C  Initialise

      DO O1 = 1, N_BRDF_STOKESSQ
        ROSSTHIN_VKERNEL(O1) = ZERO
      ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

C  (1,1) Function (scalar form)

      DS1 = XI * XJ
      DS2 = SXI * SXJ
      CKSI = DS1 + DS2 * CKPHI
      IF ( CKSI.GT.ONE ) CKSI = ONE
      SKSI = DSQRT(ONE-CKSI*CKSI)
      KSI = DACOS(CKSI)
      FUNC = ((PIO2-KSI)*CKSI + SKSI)/DS1
      ROSSTHIN_VKERNEL(1) = FUNC - PIO2

C  Other functions

C      Placeholder

C  Finish

      RETURN
      END

C

      SUBROUTINE ROSSTHICK_VFUNCTION 
     I       ( MAXPARS, N_BRDF_STOKESSQ, NPARS, PARS,
     I         XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,
     O         ROSSTHICK_VKERNEL )

C  include file of constants

      INCLUDE '../includes/VLIDORT.PARS'

C  Subroutine arguments

      INTEGER          MAXPARS, NPARS, N_BRDF_STOKESSQ
      DOUBLE PRECISION PARS ( MAXPARS )
      DOUBLE PRECISION XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION ROSSTHICK_VKERNEL(MAXSTOKES_SQ)

C  Local variables

      INTEGER          O1
      DOUBLE PRECISION DS1, DS2, DS3, CKSI, SKSI, KSI, FUNC
      DOUBLE PRECISION XPHI, CKPHI

C  Initialise

      DO O1 = 1, N_BRDF_STOKESSQ 
        ROSSTHICK_VKERNEL(O1) = ZERO
      ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

C  (1,1) Function (scalar form)

      DS1 = XI * XJ
      DS2 = SXI * SXJ
      DS3 = XI  + XJ
      CKSI = DS1 + DS2 * CKPHI
      IF ( CKSI.GT.ONE ) CKSI = ONE
      SKSI = DSQRT(ONE-CKSI*CKSI)
      KSI = DACOS(CKSI)
      FUNC = ((PIO2-KSI)*CKSI + SKSI)/DS3
      ROSSTHICK_VKERNEL(1) = FUNC - PIO4

C  Other functions

C      Placeholder

C  Finish

      RETURN
      END

C

      SUBROUTINE ROUJEAN_VFUNCTION 
     I       ( MAXPARS, N_BRDF_STOKESSQ, NPARS, PARS,
     I         XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,
     O         ROUJEAN_VKERNEL )

C  include file of constants

      INCLUDE '../includes/VLIDORT.PARS'

C  Subroutine arguments

      INTEGER          MAXPARS, NPARS, N_BRDF_STOKESSQ
      DOUBLE PRECISION PARS ( MAXPARS )
      DOUBLE PRECISION XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION ROUJEAN_VKERNEL(MAXSTOKES_SQ)

C  Local variables

      INTEGER          O1
      DOUBLE PRECISION DS1, DS2, DS3, TXJ, TXI, PHIFAC, S1, S2
      DOUBLE PRECISION XPHI_R, CXPHI_R, SXPHI_R, XPHI_C
      DOUBLE PRECISION XPHI, CKPHI

C  Initialise

      DO O1 = 1, N_BRDF_STOKESSQ 
        ROUJEAN_VKERNEL(O1) = ZERO
      ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

C  (1,1) Function (scalar form)

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
      ROUJEAN_VKERNEL(1) = S1 - S2

C  Other functions

C      Placeholder

C  Finish

      RETURN
      END

C

      SUBROUTINE LISPARSE_VFUNCTION
     I       ( MAXPARS, N_BRDF_STOKESSQ, NPARS, PARS,
     I         XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,
     O         LISPARSE_VKERNEL )

C  include file of constants

      INCLUDE '../includes/VLIDORT.PARS'

C  Subroutine arguments

      INTEGER          MAXPARS, NPARS, N_BRDF_STOKESSQ
      DOUBLE PRECISION PARS ( MAXPARS )
      DOUBLE PRECISION XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION LISPARSE_VKERNEL(MAXSTOKES_SQ)

C  local variables

      INTEGER          O1
      DOUBLE PRECISION X_INC, X_REF, SX_INC, SX_REF, ANG_P, TX
      DOUBLE PRECISION T_INC, T_REF, T_INC_SQ, T_REF_SQ
      DOUBLE PRECISION CKSI, DELTA, T, COST, SINT, DSQ, SINTCOST
      DOUBLE PRECISION A, B, H, R, P, Q, DT1, DT2, DT2SQ, QR
      DOUBLE PRECISION XPHI, CKPHI

C  Initialise
C    -- Return for special case

      DO O1 = 1, N_BRDF_STOKESSQ 
        LISPARSE_VKERNEL(O1) = ZERO
      ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

C  (1,1) Function (scalar form)

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
      LISPARSE_VKERNEL(1) = HALF * P - QR

C  Other functions

C      Placeholder

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDENSE_VFUNCTION 
     I       ( MAXPARS, N_BRDF_STOKESSQ, NPARS, PARS,
     I         XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,
     O         LIDENSE_VKERNEL )

C  include file of constants

      INCLUDE '../includes/VLIDORT.PARS'

C  Subroutine arguments

      INTEGER          MAXPARS, NPARS, N_BRDF_STOKESSQ
      DOUBLE PRECISION PARS ( MAXPARS )
      DOUBLE PRECISION XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION LIDENSE_VKERNEL(MAXSTOKES_SQ)

C  local variables

      INTEGER          O1
      DOUBLE PRECISION X_INC, X_REF, SX_INC, SX_REF, ANG_P, TX
      DOUBLE PRECISION T_INC, T_REF, T_INC_SQ, T_REF_SQ
      DOUBLE PRECISION CKSI, DELTA, T, COST, SINT, DSQ, SINTCOST
      DOUBLE PRECISION A, B, H, R, P, Q, DT1, DT2, DT2SQ, P_QR
      DOUBLE PRECISION XPHI, CKPHI

C  Initialise
C    -- Return for special case

      DO O1 = 1, N_BRDF_STOKESSQ 
        LIDENSE_VKERNEL(O1) = ZERO
      ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

C  (1,1) Function (scalar form)

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
      LIDENSE_VKERNEL(1) = P_QR - TWO

C  Other functions

C      Placeholder

C  Finish

      RETURN
      END

C

      SUBROUTINE HAPKE_VFUNCTION 
     I       ( MAXPARS, N_BRDF_STOKESSQ, NPARS, PARS,
     I         XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,
     O         HAPKE_VKERNEL )

C  include file of constants

      INCLUDE '../includes/VLIDORT.PARS'

C  Subroutine arguments

      INTEGER          MAXPARS, NPARS, N_BRDF_STOKESSQ
      DOUBLE PRECISION PARS ( MAXPARS )
      DOUBLE PRECISION XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION HAPKE_VKERNEL(MAXSTOKES_SQ)

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

      INTEGER          O1
      DOUBLE PRECISION CTHETA, THETA, PHASE
      DOUBLE PRECISION HOTSPOT, B0_EMPIR, HELP_HOT, B_HOT
      DOUBLE PRECISION SSALBEDO, GAMMA, REFLEC, FUNCTION
      DOUBLE PRECISION HELP_J, TERM_J, HELP_I, TERM_I
      DOUBLE PRECISION XPHI, CKPHI

C  Initialise

      DO O1 = 1, N_BRDF_STOKESSQ 
        HAPKE_VKERNEL(O1) = ZERO
      ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

C  (1,1) Function (scalar form)

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
      HAPKE_VKERNEL(1) = REFLEC * FUNCTION
 
C  Other functions

C      Placeholder

C  Finish

      RETURN
      END


C

      SUBROUTINE RAHMAN_VFUNCTION 
     I       ( MAXPARS, N_BRDF_STOKESSQ, NPARS, PARS,
     I         XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,
     O         RAHMAN_VKERNEL )

C  include file of constants

      INCLUDE '../includes/VLIDORT.PARS'

C  Subroutine arguments

      INTEGER          MAXPARS, NPARS, N_BRDF_STOKESSQ
      DOUBLE PRECISION PARS ( MAXPARS )
      DOUBLE PRECISION XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION RAHMAN_VKERNEL(MAXSTOKES_SQ)

C  local variables

      INTEGER          O1
      DOUBLE PRECISION T_INC, T_REF, DT1, DT2 
      DOUBLE PRECISION CXI, DELTA, K1_SQ, FACT
      DOUBLE PRECISION GEOM, PHASE, RFAC, K0, K1, K2
      DOUBLE PRECISION XPHI, CKPHI

C  Initialise

      DO O1 = 1, N_BRDF_STOKESSQ 
        RAHMAN_VKERNEL(O1) = ZERO
      ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

C  (1,1) Function (scalar form)

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
      RAHMAN_VKERNEL(1) = K0 * PHASE * ( ONE + RFAC ) * GEOM

C  Other functions

C      Placeholder

C  Finish

      RETURN
      END

C

      SUBROUTINE HAPKE_VFUNCTION_OLD 
     I       ( MAXPARS, N_BRDF_STOKESSQ, NPARS, PARS,
     I         XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,
     O         HAPKE_VKERNEL )

C  include file of constants

      INCLUDE '../includes/VLIDORT.PARS'

C  Subroutine arguments

      INTEGER          MAXPARS, NPARS, N_BRDF_STOKESSQ
      DOUBLE PRECISION PARS ( MAXPARS )
      DOUBLE PRECISION XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION HAPKE_VKERNEL(MAXSTOKES_SQ)
      DOUBLE PRECISION XPHI, CKPHI

C  Hapke Kernel function.
C   From DISORT code; used as validation.

C  local variables

      INTEGER          O1
      REAL             MU, MUP, DPHI, HAPKEKER
      EXTERNAL         HAPKEKER

C  Initialise

      DO O1 = 1, N_BRDF_STOKESSQ 
        HAPKE_VKERNEL(O1) = ZERO
      ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

C  (1,1) Function (scalar form)

      MUP = SNGL(XJ)
      MU = SNGL(XI)
      DPHI = SNGL(XPHI)
      HAPKE_VKERNEL(1) = HAPKEKER( MU, MUP, DPHI )

C  Other functions

C      Placeholder

C  Finish

      RETURN
      END

C

      REAL FUNCTION  hapkeker( MU, MUP, DPHI )

c      Supplies surface bi-directional reflectivity.
c
c      NOTE 1: Bidirectional reflectivity in DISORT is defined
c              by Eq. 39 in STWL.
c      NOTE 2: Both MU and MU0 (cosines of reflection and incidence
c              angles) are positive.
c
c  INPUT:
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

      REAL      DPHI, MU, MUP
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

      SUBROUTINE COXMUNK_VFUNCTION 
     I       ( MAXPARS, N_BRDF_STOKESSQ, NPARS, PARS,
     I         XJ, SXJ, XI, SXI, XPHI, CKPHI, SKPHI,
     O         COXMUNK_VKERNEL )

C  include file of constants

      INCLUDE '../includes/VLIDORT.PARS'

C  Subroutine arguments

      INTEGER          MAXPARS, NPARS, N_BRDF_STOKESSQ
      DOUBLE PRECISION PARS ( MAXPARS )
      DOUBLE PRECISION XI, SXI, XJ, SXJ, XPHI, CKPHI, SKPHI
      DOUBLE PRECISION COXMUNK_VKERNEL(MAXSTOKES_SQ)

C  Critical exponent taken out

      DOUBLE PRECISION CRITEXP
      PARAMETER       ( CRITEXP = 88.0D0 )

C  Local variables

      INTEGER          O1
      DOUBLE PRECISION Z, Z1, Z2, Z2_SQ_M1, H1, H2, RP, RL, XMP
      DOUBLE PRECISION A, B, TA, ARGUMENT, PROB, FAC1, FAC2, CKPHI_NEG
      DOUBLE PRECISION S1, S2, S3, XXI, XXJ, T1, T2, DCOT
      DOUBLE PRECISION SHADOWI, SHADOWR, SHADOW

C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C               Remark on Use of shadow effect
C               ------------------------------
C  Shadow effect is controlled by the third parameter. That is, if
C  PARS(3) not equal to then shadow effect will be included.
C    --- NPARS should always be 3 for this Kernel.
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

C  Initialise

      DO O1 = 1, N_BRDF_STOKESSQ 
        COXMUNK_VKERNEL(O1) = ZERO
      ENDDO

C  Comment. 18 January 2006.
C  We have found in comparisons with the Giss Cox-Munk code that
C  the input COSPHI (CKPHI) here is the negative of what we actually need,
C  so introduce local variable which takes care of this

C  Also removed factor of PIE in the kernel denominator
C   This makes the output exactly same as GISS model for R(1,1)

      CKPHI_NEG = - CKPHI

C  (1,1) Function (scalar form)

C  ..Scatter angles

C old   Z = - XI * XJ + SXI * SXJ * CKPHI   
C old   IF ( Z .LT. MINUS_ONE) Z = MINUS_ONE
C old   Z1 = DACOS(-Z)
C old   Z2 = DCOS(Z1*HALF)

      Z = XI * XJ + SXI * SXJ * CKPHI_NEG 
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
        FAC1 = PROB / PARS(1)
        FAC2 = QUARTER / XI / ( B ** FOUR )
        COXMUNK_VKERNEL(1) = XMP * FAC1 * FAC2 / XJ
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

      COXMUNK_VKERNEL(1) = COXMUNK_VKERNEL(1) * SHADOW

C  Finish

      RETURN
      END

C

      SUBROUTINE GISSCOXMUNK_VFUNCTION 
     I       ( MAXPARS, N_BRDF_STOKESSQ, NPARS, PARS,
     I         XI, SXI, XJ, SXJ,
     I         XPHI_REF, CKPHI_REF, SKPHI_REF,
     O         GISSCOXMUNK_VKERNEL )

C  include file of constants

      INCLUDE '../includes/VLIDORT.PARS'

C  Subroutine arguments

      INTEGER          MAXPARS, NPARS, N_BRDF_STOKESSQ
      DOUBLE PRECISION PARS ( MAXPARS )
      DOUBLE PRECISION XI, SXI, XJ, SXJ
      DOUBLE PRECISION XPHI_REF, CKPHI_REF, SKPHI_REF
      DOUBLE PRECISION GISSCOXMUNK_VKERNEL(MAXSTOKES_SQ)

C  Critical exponent taken out

      DOUBLE PRECISION CRITEXP
      PARAMETER       ( CRITEXP = 88.0D0 )

C  Local variables

      INTEGER          O1, KERNELMASK(10), I, IM
      DOUBLE PRECISION XPHI_INC, CKPHI_INC, SKPHI_INC
      DOUBLE PRECISION VI1, VI2, VI3, VR1, VR2, VR3
      DOUBLE PRECISION unit1, unit2, unit3, fact1, factor
      DOUBLE PRECISION XI1, CN1, CN2, CXI2, C2, C1, CRPER, CRPAR
      DOUBLE PRECISION TI1, TI2, TI3, TR1, TR2, TR3
      DOUBLE PRECISION PI1, PII2, PI3, PR1, PR2, PR3
      DOUBLE PRECISION PIKR, PRKI, TIKR, TRKI
      DOUBLE PRECISION E1, E2, E3, E4, SIGMA2
      DOUBLE PRECISION CF11, CF12, CF21, CF22, RDZ2, RDZ4
      DOUBLE PRECISION VP1, VP2, VP3, DMOD, DEX, DCOEFF
      DOUBLE PRECISION AF, AF11, AF12, AF21, AF22
      DOUBLE PRECISION C21, C22, CTTTP, CTTPT, CTTPP
      DOUBLE PRECISION CTPPT, CTPPP, CPTPP
      DOUBLE PRECISION S1, S2, S3, XXI, XXJ, T1, T2, DCOT
      DOUBLE PRECISION SHADOWI, SHADOWR, SHADOW, DCOEFF_0, ARGUMENT
      DATA KERNELMASK / 1,2,3,5,6,7,9,10,11,16 /

C  Initialise

      DO O1 = 1, N_BRDF_STOKESSQ 
        GISSCOXMUNK_VKERNEL(O1) = ZERO
      ENDDO

C  Transcription of the RMATR subroutine from Mishchenko/Travis code.

C   CALCULATION OF THE STOKES REFLECTION MATRIX FOR
C   ILLUMINATION FROM ABOVE FOR
C   A STATISTICALLY ROUGH SURFACE SEPARATING TWO HALF-SPACES
C   WITH REFRACTIVE INDICES OF THE UPPER AND LOWER HALF-SPACES EQUAL TO
C   CN1 AND CN2, RESPECTIVELY. THE EFFECT OF SHADOWING IS NOT
C   INCLUDED IN THIS SUBROUTINE BUT IS ADDED IN THE MAIN PROGRAM.

C   SIGMA2 = s**2 = MEAN SQUARE SURFACE SLOPE (EQ. (18) IN THE JGR PAPER) 

C   XI = ABS(COSINE OF THE INCIDENT ZENITH ANGLE)
C   XJ = ABS(COSINE OF THE REFLECTION ZENITH ANGLE).
C   SXI and SXJ are the respective SINES (input)
C   XPHI_REF = REFLECTION AZIMUTH ANGLE
C   GISSCOXMUNK_VKERNEL(16-elements) = (4X4) REFLECTION MATRIX

C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C               Remark on Use of shadow effect
C               ------------------------------
C  Shadow effect is controlled by the third parameter. That is, if
C  PARS(3) not equal to then shadow effect will be included.
C    --- NPARS should always be 3 for this Kernel.
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

C  For real case, incident azimuth taken to be zero

      XPHI_INC  = ZERO
      CKPHI_INC = ONE
      SKPHI_INC = ZERO

C  Slope square is PARS(1)

      SIGMA2 = PARS(1)

C  Check for limiting cases

      IF(DABS(XI-1D0).LT.1d-9) XI = 0.999999999999d0
      IF(DABS(XJ-1D0).LT.1d-9) XJ = 0.999999999999d0

C  help variables (coordinate transformations)

      VI1 = SXI * CKPHI_INC
      VI2 = SXI * SKPHI_INC
      VI3 = -XI
      VR1 = SXJ * CKPHI_REF
      VR2 = SXJ * SKPHI_REF
      VR3 = XJ

C    LOCAL SURFACE NORMAL FOR SPECULAR REFLECTION (normalized to 1)

      UNIT1  = VI1-VR1
      UNIT2  = VI2-VR2
      UNIT3  = VI3-VR3
      FACT1  = UNIT1*UNIT1 + UNIT2*UNIT2 + UNIT3*UNIT3
      FACTOR = DSQRT(ONE/FACT1)

C   FRESNEL REFLECTION COEFFICIENTS, assume only real for now
C   ---------------------------------------------------------

      CN1 = ONE
      CN2 = PARS(2)

C  this is the original code, but now C-variables are real

      XI1 =  FACTOR*(UNIT1*VI1+UNIT2*VI2+UNIT3*VI3)
      CXI2 = ONE - (ONE-XI1*XI1)*CN1*CN1/(CN2*CN2)
      CXI2 = DSQRT(CXI2)
      C1 = CN1*XI1
      C2 = CN2*CXI2
      CRPER = (C1-C2)/(C1+C2)
      C1 = CN2*XI1
      C2 = CN1*CXI2
      CRPAR = (C1-C2)/(C1+C2)

C  CALCULATION OF THE AMPLITUDE SCATTERING MATRIX
C  ----------------------------------------------

      TI1 = - XI * CKPHI_INC
      TI2 = - XI * SKPHI_INC
      TI3 = - SXI

      TR1 = + XJ * CKPHI_REF
      TR2 = + XJ * SKPHI_REF
      TR3 = -SXJ

      PI1  = - SKPHI_INC
      PII2 = + CKPHI_INC
      PI3  = ZERO

      PR1 = - SKPHI_REF
      PR2 = + CKPHI_REF
      PR3 = ZERO

      PIKR = PI1*VR1 + PII2*VR2 + PI3*VR3
      PRKI = PR1*VI1 + PR2*VI2  + PR3*VI3
      TIKR = TI1*VR1 + TI2*VR2  + TI3*VR3
      TRKI = TR1*VI1 + TR2*VI2  + TR3*VI3

      E1 = PIKR*PRKI
      E2 = TIKR*TRKI
      E3 = TIKR*PRKI
      E4 = PIKR*TRKI

      CF11 =  E1*CRPER+E2*CRPAR
      CF12 = -E3*CRPER+E4*CRPAR
      CF21 = -E4*CRPER+E3*CRPAR
      CF22 =  E2*CRPER+E1*CRPAR

C  CALCULATION OF THE STOKES REFLECTION MATRIX
C  -------------------------------------------

      VP1 = VI2*VR3-VI3*VR2
      VP2 = VI3*VR1-VI1*VR3
      VP3 = VI1*VR2-VI2*VR1
      DMOD = VP1*VP1+VP2*VP2+VP3*VP3
      DMOD = DMOD*DMOD

C  if DMOD = 0, that is | n x n_0 | ^ 4 = 0 in M-T formula)
C    Then we need to set the ratio CF11 / DMOD

      IF ( DMOD .EQ. ZERO ) THEN
        CF11 = CRPAR
        CF22 = CRPER
        DMOD = 1.0d0
      ENDIF

      RDZ2 = UNIT3*UNIT3
      RDZ4 = RDZ2*RDZ2

      DCOEFF_0   = ONE/(8.0D0*XI*XJ*DMOD*RDZ4*SIGMA2)
      ARGUMENT  = (UNIT1*UNIT1 + UNIT2*UNIT2) / (TWO*SIGMA2*RDZ2)
      IF ( ARGUMENT .GT. CRITEXP ) THEN
        DEX = ZERO
      ELSE
        DEX = DEXP(-ARGUMENT)
      ENDIF

      DCOEFF = DCOEFF_0 * FACT1 * FACT1 * DEX

      AF  = HALF * DCOEFF
      AF11 = DABS(CF11)
      AF12 = DABS(CF12)
      AF21 = DABS(CF21)
      AF22 = DABS(CF22)
      AF11 = AF11*AF11
      AF12 = AF12*AF12
      AF21 = AF21*AF21
      AF22 = AF22*AF22

C  original code
c      R(1,1)=(AF11+AF12+AF21+AF22)*AF
c      R(1,2)=(AF11-AF12+AF21-AF22)*AF
c      R(2,1)=(AF11-AF22+AF12-AF21)*AF
c      R(2,2)=(AF11-AF12-AF21+AF22)*AF

C  Transcribed code
      GISSCOXMUNK_VKERNEL(1) = (AF11+AF12+AF21+AF22)*AF
      GISSCOXMUNK_VKERNEL(2) = (AF11-AF12+AF21-AF22)*AF
      GISSCOXMUNK_VKERNEL(5) = (AF11-AF22+AF12-AF21)*AF
      GISSCOXMUNK_VKERNEL(6) = (AF11-AF12-AF21+AF22)*AF

C  Key Debug statement to track down DMOD
C      write(*,*)'(1,1) Fresnel = ',(AF11+AF12+AF21+AF22)*HALF/DMOD

C  Original code
C      CI=(0D0, -1D0)
C      C21=DCONJG(CF21)
C      C22=DCONJG(CF22)
C      CTTTP=CF11*DCONJG(CF12)
c      CTTPT=CF11*C21
c      CTTPP=CF11*C22
c      CTPPT=CF12*C21
c      CTPPP=CF12*C22
c      CPTPP=CF21*C22

C  replica for real variables only
      C21 = CF21
      C22 = CF22
      CTTTP=CF11*CF12
      CTTPT=CF11*C21
      CTTPP=CF11*C22
      CTPPT=CF12*C21
      CTPPP=CF12*C22
      CPTPP=CF21*C22

C  original code
c      R(1,3)=    (-CTTTP-CPTPP)*DCOEFF
c      R(1,4)=-CI*( CTTTP+CPTPP)*DCOEFF
c      R(2,3)=    (-CTTTP+CPTPP)*DCOEFF
c      R(2,4)=-CI*( CTTTP-CPTPP)*DCOEFF
c      R(3,1)=    (-CTTPT-CTPPP)*DCOEFF
c      R(3,2)=    (-CTTPT+CTPPP)*DCOEFF
c      R(3,3)=    ( CTTPP+CTPPT)*DCOEFF
c      R(3,4)= CI*( CTTPP-CTPPT)*DCOEFF
c      R(4,1)= CI*( CTTPT+CTPPP)*DCOEFF
c      R(4,2)= CI*( CTTPT-CTPPP)*DCOEFF
c      R(4,3)=-CI*( CTTPP+CTPPT)*DCOEFF
c      R(4,4)=    ( CTTPP-CTPPT)*DCOEFF

c  New code (several entries are zero)

      GISSCOXMUNK_VKERNEL(3)  =    (-CTTTP-CPTPP)*DCOEFF
      GISSCOXMUNK_VKERNEL(7)  =    (-CTTTP+CPTPP)*DCOEFF
      GISSCOXMUNK_VKERNEL(9)  =    (-CTTPT-CTPPP)*DCOEFF
      GISSCOXMUNK_VKERNEL(10) =    (-CTTPT+CTPPP)*DCOEFF
      GISSCOXMUNK_VKERNEL(11) =    ( CTTPP+CTPPT)*DCOEFF
      GISSCOXMUNK_VKERNEL(16) =    ( CTTPP-CTPPT)*DCOEFF

C  No Shadow code if not flagged

      IF ( PARS(3) .EQ. ZERO ) RETURN

C  Shadow code

      S1 = DSQRT(TWO*SIGMA2/PIE)
      S3 = ONE/(DSQRT(TWO*SIGMA2))
      S2 = S3*S3

      IF ( XI .EQ. ONE ) THEN
       SHADOWI   = ZERO
      ELSE
       XXI  = XI*XI
       DCOT = XI/DSQRT(ONE-XXI)
       T1   = DEXP(-DCOT*DCOT*S2)
       T2   = DERFC(DCOT*S3)
       SHADOWI = HALF*(S1*T1/DCOT-T2)
      ENDIF

      IF ( XJ .EQ. ONE ) THEN
       SHADOWR   = ZERO
      ELSE
       XXJ  = XJ*XJ
       DCOT = XJ/DSQRT(ONE-XXJ)
       T1   = DEXP(-DCOT*DCOT*S2)
       T2   = DERFC(DCOT*S3)
       SHADOWR = HALF*(S1*T1/DCOT-T2)
      ENDIF

      SHADOW = ONE/(ONE+SHADOWI+SHADOWR)

      DO I = 1, 10
       IM = KERNELMASK(I)
       GISSCOXMUNK_VKERNEL(IM) = GISSCOXMUNK_VKERNEL(IM) * SHADOW
      ENDDO

C  debug

c      DO M = 1, N_BRDF_STOKESSQ
c        write(33,'(I5,1p6e14.5)')m,
c     &   GISSCOXMUNK_VKERNEL(m),
c     &   dacos(xi)/deg_to_rad, dacos(xj)/deg_to_rad, phi
c      enddo

C  Finish

      RETURN
      END

C

      SUBROUTINE COXMUNK_VFUNCTION_DB
     I       ( MAXPARS, N_BRDF_STOKESSQ, NPARS, PARS,
     I         XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,
     O         COXMUNK_VKERNEL )

C  include file of constants

      INCLUDE '../includes/VLIDORT.PARS'

C  Subroutine arguments

      INTEGER          MAXPARS, NPARS, N_BRDF_STOKESSQ
      DOUBLE PRECISION PARS ( MAXPARS )
      DOUBLE PRECISION XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION COXMUNK_VKERNEL(MAXSTOKES_SQ)

C  local variables
C  ---------------

C  help variables

      integer          n, k, i, i1, N_phiquad_HALF, M, o1
      DOUBLE PRECISION XM, SXM, sum_pr, sumr, sum, w_p
      DOUBLE PRECISION reflec_0(MAXSTOKES_SQ), reflec_1(MAXSTOKES_SQ)
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

      DOUBLE PRECISION R0_QUAD_IN
     &      (MAXSTOKES_SQ,max_msrs_muquad,max_msrs_phiquad)
      DOUBLE PRECISION R0_OUT_QUAD
     &      (MAXSTOKES_SQ,max_msrs_muquad,max_msrs_phiquad)

C  Safety first zeroing

      DO O1 = 1, N_BRDF_STOKESSQ
        REFLEC_0(O1) = ZERO
        REFLEC_1(O1) = ZERO
        COXMUNK_VKERNEL(O1) = ZERO
      ENDDO

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

      CALL COXMUNK_VFUNCTION
     I       ( MAXPARS, N_BRDF_STOKESSQ, NPARS, PARS,
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
          CALL COXMUNK_VFUNCTION
     I       ( MAXPARS, N_BRDF_STOKESSQ, NPARS, PARS,
     I         XM, SXM, XI, SXI, PHI_SUB2, CPHI_SUB2, SPHI_SUB2,
     O         R0_OUT_QUAD(1,K,N) )
          CALL COXMUNK_VFUNCTION
     I       ( MAXPARS, N_BRDF_STOKESSQ, NPARS, PARS,
     I         XJ, SXJ, XM, SXM, PHI_SUB1, CPHI_SUB1, SPHI_SUB1,
     O         R0_QUAD_IN(1,K,N) )
        ENDDO
      ENDDO

C  compute the next order, (1,1) component only

      SUMR = ZERO
      DO K = 1, n_muquad
        SUM_PR = ZERO
        DO N = 1, N_PHIQUAD
          W_P  = W_PHIQUAD(N)
          SUM = R0_QUAD_IN(1,K,N) * R0_OUT_QUAD(1,K,N)
          SUM_PR = SUM_PR + W_P * SUM
        ENDDO
        SUMR =  SUMR + SUM_PR * WXX_MUQUAD(K)
      ENDDO
      REFLEC_1(1) = SUMR

C  Compute total

      DO M = 1, N_BRDF_STOKESSQ
        COXMUNK_VKERNEL(M) = REFLEC_0(M) + REFLEC_1(M)
      ENDDO
      write(*,*)'unpt',reflec_0(1),reflec_1(1)

C  debug

c      DO M = 1, N_BRDF_STOKESSQ
c        write(34,'(I5,1p6e14.5)')m,
c     &   reflec_0(m), reflec_1(m),COXMUNK_VKERNEL(m),
c     &   dacos(xi)/deg_to_rad, dacos(xj)/deg_to_rad, phi
c      enddo

C  Finish

      RETURN
      END

C

      SUBROUTINE GISSCOXMUNK_VFUNCTION_DB
     I       ( MAXPARS, N_BRDF_STOKESSQ, NPARS, PARS,
     I         XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI,
     O         GISSCOXMUNK_VKERNEL )

C  include file of constants

      INCLUDE '../includes/VLIDORT.PARS'

C  Subroutine arguments

      INTEGER          MAXPARS, NPARS, N_BRDF_STOKESSQ
      DOUBLE PRECISION PARS ( MAXPARS )
      DOUBLE PRECISION XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION GISSCOXMUNK_VKERNEL(MAXSTOKES_SQ)

C  local variables
C  ---------------

C  help variables

      integer          n, k, i, i1, N_phiquad_HALF,O1,O2,O3
      DOUBLE PRECISION XM, SXM, sum_pr(16), sumr(16), sum, w_p
      DOUBLE PRECISION reflec_0(MAXSTOKES_SQ), reflec_1(MAXSTOKES_SQ)
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

      DOUBLE PRECISION R0_QUAD_IN
     &      (MAXSTOKES_SQ,max_msrs_muquad,max_msrs_phiquad)
      DOUBLE PRECISION R0_OUT_QUAD
     &      (MAXSTOKES_SQ,max_msrs_muquad,max_msrs_phiquad)

C  Indices

      INTEGER          MASKIT(3,3),LNS,NS,NSS,M,M12,M13,M32

C  Safety first zeroing

      DO O1 = 1, N_BRDF_STOKESSQ
        REFLEC_0(O1) = ZERO
        REFLEC_1(O1) = ZERO
        GISSCOXMUNK_VKERNEL(O1) = ZERO
      ENDDO

C  Masking limits

      NSS = N_BRDF_STOKESSQ
      IF ( N_BRDF_STOKESSQ.EQ.1  ) LNS = 1
      IF ( N_BRDF_STOKESSQ.EQ.4  ) LNS = 2
      IF ( N_BRDF_STOKESSQ.EQ.9  ) LNS = 3
      NS = LNS
      IF ( N_BRDF_STOKESSQ.EQ.16 ) THEN
        LNS = 3
        NS = LNS + 1 
      ENDIF

C  masking array

      MASKIT(1,1) = 1
      MASKIT(1,2) = 2
      MASKIT(1,3) = 3
      MASKIT(2,1) = 5
      MASKIT(2,2) = 6
      MASKIT(2,3) = 7
      MASKIT(3,1) = 9
      MASKIT(3,2) = 10
      MASKIT(3,3) = 11

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

      CALL GISSCOXMUNK_VFUNCTION
     I       ( MAXPARS, N_BRDF_STOKESSQ, NPARS, PARS,
     I         XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI,
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
          CALL GISSCOXMUNK_VFUNCTION
     I       ( MAXPARS, N_BRDF_STOKESSQ, NPARS, PARS,
     I         XI, SXI, XM, SXM, PHI_SUB2, CPHI_SUB2, SPHI_SUB2,
     O         R0_OUT_QUAD(1,K,N) )
          CALL GISSCOXMUNK_VFUNCTION
     I       ( MAXPARS, N_BRDF_STOKESSQ, NPARS, PARS,
     I         XM, SXM, XJ, SXJ, PHI_SUB1, CPHI_SUB1, SPHI_SUB1,
     O         R0_QUAD_IN(1,K,N) )
c          write(67,'(2i4,1p16e12.3)')K,N,(R0_OUT_QUAD(M,K,N),M=1,16)
c          write(68,'(2i4,1p16e12.3)')K,N,(R0_QUAD_IN(M,K,N),M=1,16)
        ENDDO
      ENDDO

C  compute the next order

      DO M = 1, N_BRDF_STOKESSQ
        SUMR(M) = ZERO
      ENDDO
      DO K = 1, n_muquad
        DO M = 1, N_BRDF_STOKESSQ
          SUM_PR(M) = ZERO
        ENDDO
        DO N = 1, N_PHIQUAD
          W_P  = W_PHIQUAD(N)
          DO O1 = 1, LNS
            DO O2 = 1, LNS
              M12 = MASKIT(O1,O2)
              SUM = ZERO
              DO O3 = 1, LNS
                M13 = MASKIT(O1,O3)
                M32 = MASKIT(O3,O2)
                SUM = SUM + R0_QUAD_IN(M13,K,N) * R0_OUT_QUAD(M32,K,N)
              ENDDO  
              SUM_PR(M12) = SUM_PR(M12) + W_P * SUM
            ENDDO
          ENDDO
          IF (NS.EQ.4) THEN
            SUM = R0_QUAD_IN(NSS,K,N) * R0_OUT_QUAD(NSS,K,N)
            SUM_PR(NSS) = SUM_PR(NSS) + W_P * SUM
          ENDIF
        ENDDO
        DO M = 1, N_BRDF_STOKESSQ
          SUMR(M) =  SUMR(M) + SUM_PR(M) * WXX_MUQUAD(K)
        ENDDO
      ENDDO
      DO M = 1, N_BRDF_STOKESSQ
        REFLEC_1(M) = SUMR(M)
      ENDDO

C  Compute total

      DO M = 1, N_BRDF_STOKESSQ
        GISSCOXMUNK_VKERNEL(M) = REFLEC_0(M) + REFLEC_1(M)
      ENDDO

C  debug

c      DO M = 1, N_BRDF_STOKESSQ
c        write(34,'(I5,1p6e14.5)')m,
c     &   reflec_0(m), reflec_1(m),GISSCOXMUNK_VKERNEL(m),
c     &   dacos(xi)/deg_to_rad, dacos(xj)/deg_to_rad, phi
c      enddo

C  Finish

      RETURN
      END

C

      SUBROUTINE RHERMAN_VFUNCTION 
     I       ( MAXPARS, N_BRDF_STOKESSQ, NPARS, PARS,
     I         XI, SXI, XJ, SXJ,
     I         XPHI_REF, CKPHI_REF, SKPHI_REF,
     O         RHERMAN_VKERNEL )

C  include file of constants

      INCLUDE '../includes/VLIDORT.PARS'

C  Subroutine arguments

      INTEGER          MAXPARS, NPARS, N_BRDF_STOKESSQ
      DOUBLE PRECISION PARS ( MAXPARS )
      DOUBLE PRECISION XI, SXI, XJ, SXJ
      DOUBLE PRECISION XPHI_REF, CKPHI_REF, SKPHI_REF
      DOUBLE PRECISION RHERMAN_VKERNEL(MAXSTOKES_SQ)

C  Local variables

      INTEGER          O1

      DOUBLE PRECISION RAHMAN_VKERNEL(MAXSTOKES_SQ)
      DOUBLE PRECISION HFUNCTION

      DOUBLE PRECISION XPHI_INC, CKPHI_INC, SKPHI_INC, XPHI_RAH
      DOUBLE PRECISION VI1, VI2, VI3, VR1, VR2, VR3
      DOUBLE PRECISION VP1, VP2, VP3, DMOD
      DOUBLE PRECISION unit1, unit2, unit3, fact1, factor
      DOUBLE PRECISION XI1, CN1, CN2, CXI2, C2, C1, CRPER, CRPAR
      DOUBLE PRECISION TI1, TI2, TI3, TR1, TR2, TR3
      DOUBLE PRECISION PI1, PII2, PI3, PR1, PR2, PR3
      DOUBLE PRECISION PIKR, PRKI, TIKR, TRKI
      DOUBLE PRECISION E1, E2, E3, E4
      DOUBLE PRECISION CF11, CF12, CF21, CF22
      DOUBLE PRECISION AF11, AF12, AF21, AF22
      DOUBLE PRECISION C21, C22, CTTTP, CTTPT, CTTPP
      DOUBLE PRECISION CTPPT, CTPPP, CPTPP

C  Initialise

      DO O1 = 1, N_BRDF_STOKESSQ 
        RHERMAN_VKERNEL(O1) = ZERO
      ENDDO

C#####################################################################
C#####################################################################
C  COXMUNK scalar stuff...................
C  Also removed factor of PIE in the kernel denominator
C   This makes the output exactly same as GISS model for R(1,1)
C      CKPHI_NEG = - CKPHI
c      Z = XI * XJ + SXI * SXJ * CKPHI_NEG 
c      IF ( Z .GT. ONE) Z = ONE
c      Z1 = DACOS(Z)
c      Z2 = DCOS(Z1*HALF)
C  .. Fresnel coefficients
c      REFIDX = 1.5d0
c      REFIDX_SQ = REFIDX * REFIDX
c      Z2_SQ_M1 = Z2 * Z2 + MINUS_ONE
c      H1 = REFIDX_SQ * Z2
c      H2 = DSQRT ( REFIDX_SQ + Z2_SQ_M1 )
c      RP = ( H1 - H2 ) / ( H1 + H2 )
c      RL = ( Z2 - H2 ) / ( Z2 + H2 )
C  XMP = r(i) = (1,1) reflection coefficient for natural lights
C      XMP = HALF * ( RP*RP + RL*RL )
C  Comment. This is the Lenoble SOS code
c      xind=1.50
c      xx=cos(pi*thd/360.)
c      yy=sqrt(xind*xind+xx*xx-1.0)
c      zz=xx*xind*xind
c      rl=(zz-yy)/(zz+yy)
c      rr=(xx-yy)/(xx+yy)
C  where
C     thd  = z1 * 180 / pi  (in degrees, Z1 is in radians)
C     xx   = Z2
C     yy   = H2
C     zz   = H1
C     xind = REFIDX
C     rl   = RP
C     rr   = RL
C  Coxmunk Function
c      A = TWO * Z2
c      B = ( XI + XJ ) / A
c      IF ( B .GT. ONE ) B = ONE
c      A = PIO2 - DASIN(B)
c      TA = DTAN(A)
c      ARGUMENT = TA * TA  / PARS(1)
c      IF ( ARGUMENT .LT. CRITEXP ) THEN
c        PROB = DEXP ( - ARGUMENT )
c        FAC1 = PROB / PARS(1)
c        FAC2 = QUARTER / XI / ( B ** FOUR )
c        COXMUNK_VKERNEL(1) = XMP * FAC1 * FAC2 / XJ
c      ENDIF
C#####################################################################
C#####################################################################

C  Transcription of the RMATR subroutine from Mishchenko/Travis code.

C   CALCULATION OF THE STOKES REFLECTION MATRIX FOR
C   ILLUMINATION FROM ABOVE FOR A SURFACE SEPARATING TWO HALF-SPACES
C   WITH REFRACTIVE INDICES OF THE UPPER AND LOWER HALF-SPACES EQUAL TO
C   CN1 AND CN2, RESPECTIVELY.

C   XI = ABS(COSINE OF THE INCIDENT ZENITH ANGLE)
C   XJ = ABS(COSINE OF THE REFLECTION ZENITH ANGLE).
C   SXI and SXJ are the respective SINES (input)
C   XPHI_REF = REFLECTION AZIMUTH ANGLE
C   RHERMAN_VKERNEL(16-elements) = (4X4) REFLECTION MATRIX

C  For real case, incident azimuth taken to be zero

      XPHI_INC  = ZERO
      CKPHI_INC = ONE
      SKPHI_INC = ZERO

C  Check for limiting cases

      IF(DABS(XI-1D0).LT.1d-9) XI = 0.999999999999d0
      IF(DABS(XJ-1D0).LT.1d-9) XJ = 0.999999999999d0

C  help variables (coordinate transformations)

      VI1 = SXI * CKPHI_INC
      VI2 = SXI * SKPHI_INC
      VI3 = -XI
      VR1 = SXJ * CKPHI_REF
      VR2 = SXJ * SKPHI_REF
      VR3 = XJ

C    LOCAL SURFACE NORMAL FOR SPECULAR REFLECTION (normalized to 1)

      UNIT1  = VI1-VR1
      UNIT2  = VI2-VR2
      UNIT3  = VI3-VR3
      FACT1  = UNIT1*UNIT1 + UNIT2*UNIT2 + UNIT3*UNIT3
      FACTOR = DSQRT(ONE/FACT1)

C   FRESNEL REFLECTION COEFFICIENTS, assume only real for now
C   ---------------------------------------------------------

      CN1 = ONE
      CN2 = 1.5d0

C  this is the original code, but now C-variables are real

      XI1 =  FACTOR*(UNIT1*VI1+UNIT2*VI2+UNIT3*VI3)
      CXI2 = ONE - (ONE-XI1*XI1)*CN1*CN1/(CN2*CN2)
      CXI2 = DSQRT(CXI2)
      C1 = CN1*XI1
      C2 = CN2*CXI2
      CRPER = (C1-C2)/(C1+C2)
      C1 = CN2*XI1
      C2 = CN1*CXI2
      CRPAR = (C1-C2)/(C1+C2)

C  CALCULATION OF THE AMPLITUDE SCATTERING MATRIX
C  ----------------------------------------------

      TI1 = - XI * CKPHI_INC
      TI2 = - XI * SKPHI_INC
      TI3 = - SXI

      TR1 = + XJ * CKPHI_REF
      TR2 = + XJ * SKPHI_REF
      TR3 = -SXJ

      PI1  = - SKPHI_INC
      PII2 = + CKPHI_INC
      PI3  = ZERO

      PR1 = - SKPHI_REF
      PR2 = + CKPHI_REF
      PR3 = ZERO

      PIKR = PI1*VR1 + PII2*VR2 + PI3*VR3
      PRKI = PR1*VI1 + PR2*VI2  + PR3*VI3
      TIKR = TI1*VR1 + TI2*VR2  + TI3*VR3
      TRKI = TR1*VI1 + TR2*VI2  + TR3*VI3

      E1 = PIKR*PRKI
      E2 = TIKR*TRKI
      E3 = TIKR*PRKI
      E4 = PIKR*TRKI

C  Set the H-function

      HFUNCTION = 0.25d0 / ( XI + XJ )

C  Settting (1,1), (1,2), (2,1) and (2,2) components
C  -----------------------------------------------

      CF11 =  E1*CRPER+E2*CRPAR
      CF12 = -E3*CRPER+E4*CRPAR
      CF21 = -E4*CRPER+E3*CRPAR
      CF22 =  E2*CRPER+E1*CRPAR

C  Not to forget the normalization

      VP1 = VI2*VR3-VI3*VR2
      VP2 = VI3*VR1-VI1*VR3
      VP3 = VI1*VR2-VI2*VR1
      DMOD = VP1*VP1+VP2*VP2+VP3*VP3
      DMOD = DMOD*DMOD

C  if DMOD = 0, that is | n x n_0 | ^ 4 = 0 in M-T formula)
C    Then we need to set the ratio CF11 / DMOD

      IF ( DMOD .EQ. ZERO ) THEN
        CF11 = CRPAR
        CF22 = CRPER
        DMOD = 1.0d0
      ENDIF

!  Final H-function needs to be normalized

      HFUNCTION = HFUNCTION / DMOD

!  Compute the Electromagnetic contributions

      AF11 = DABS(CF11)
      AF12 = DABS(CF12)
      AF21 = DABS(CF21)
      AF22 = DABS(CF22)
      AF11 = AF11*AF11
      AF12 = AF12*AF12
      AF21 = AF21*AF21
      AF22 = AF22*AF22

      FACTOR = HALF * HFUNCTION
      RHERMAN_VKERNEL(1) = (AF11+AF12+AF21+AF22) * FACTOR
      RHERMAN_VKERNEL(2) = (AF11-AF12+AF21-AF22) * FACTOR
      RHERMAN_VKERNEL(5) = (AF11-AF22+AF12-AF21) * FACTOR
      RHERMAN_VKERNEL(6) = (AF11-AF12-AF21+AF22) * FACTOR

C  Setting (1,3), (2,3), (3,1), (3,2), (3,3) and (4,4) components
C  --------------------------------------------------------------

      C21 = CF21
      C22 = CF22
      CTTTP=CF11*CF12
      CTTPT=CF11*C21
      CTTPP=CF11*C22
      CTPPT=CF12*C21
      CTPPP=CF12*C22
      CPTPP=CF21*C22

      FACTOR = HFUNCTION
      RHERMAN_VKERNEL(3)  =    (-CTTTP-CPTPP) * FACTOR
      RHERMAN_VKERNEL(7)  =    (-CTTTP+CPTPP) * FACTOR
      RHERMAN_VKERNEL(9)  =    (-CTTPT-CTPPP) * FACTOR
      RHERMAN_VKERNEL(10) =    (-CTTPT+CTPPP) * FACTOR
      RHERMAN_VKERNEL(11) =    ( CTTPP+CTPPT) * FACTOR
      RHERMAN_VKERNEL(16) =    ( CTTPP-CTPPT) * FACTOR

C  Add the Diffuse term to R(1,1)
C  ------------------------------

c   IN THE PREVIOUS EQUATION, THE BRDF MODEL WE USE IS THE SINYUK-ET-AL MODEL
c   WITH FREE PARAMETERS rho,vk AND gt (see Sinyuk et al. paper).

C  This is just the Rahman Kernel.........different name !!

      XPHI_RAH = 180.0d0 - XPHI_REF
      CALL RAHMAN_VFUNCTION 
     I       ( MAXPARS, N_BRDF_STOKESSQ, NPARS, PARS,
     I         XI, SXI, XJ, SXJ,
     I         XPHI_RAH, CKPHI_REF, SKPHI_REF,
     O         RAHMAN_VKERNEL )

C  Add to the specular term

      RHERMAN_VKERNEL(1) = RHERMAN_VKERNEL(1) + RAHMAN_VKERNEL(1)

C  Here is the original code from France........
c    cs: cosinus of the solar zenith angle
c    cv: cosinus of the viewing zenith angle
c    phi: azimuth between the solar and observation vertical planes
c      xx=cs**(vk-1.)
c      yy=cs**(vk-1.)
c      zz=(cs+cv)**(1.-vk)
c      FF1=rho*xx*yy/zz
c      xx=sqrt(1-cs*cs)
c      yy=sqrt(1-cv*cv)
c      ww=cs*cv-xx*yy*cos(pi*phi/180.)
c      aa=1+gt*gt+2*gt*ww
c      FF2=(1-gt*gt)/(aa**1.5)
c      vv=xx/cs
c      ww=yy/cv
c      G=sqrt(vv*vv+ww*ww+2*vv*ww*cos(pi*phi/180))
c      FF3=1+(1-rho)/(1+G)
c      RBD=FF1*FF2*FF3
C RBD...........................IS THE REFLECTION COEFFICIENT

C  Finish

      RETURN
      END

C

      SUBROUTINE BREON_VFUNCTION 
     I       ( MAXPARS, N_BRDF_STOKESSQ, NPARS, PARS,
     I         XI, SXI, XJ, SXJ,
     I         XPHI_REF, CKPHI_REF, SKPHI_REF,
     O         BREON_VKERNEL )

C  include file of constants

      INCLUDE '../includes/VLIDORT.PARS'

C  Subroutine arguments

      INTEGER          MAXPARS, NPARS, N_BRDF_STOKESSQ
      DOUBLE PRECISION PARS ( MAXPARS )
      DOUBLE PRECISION XI, SXI, XJ, SXJ
      DOUBLE PRECISION XPHI_REF, CKPHI_REF, SKPHI_REF
      DOUBLE PRECISION BREON_VKERNEL(MAXSTOKES_SQ)

C  Local variables

      INTEGER          O1

      DOUBLE PRECISION RAHMAN_VKERNEL(MAXSTOKES_SQ)
      DOUBLE PRECISION HFUNCTION

      DOUBLE PRECISION XPHI_INC, CKPHI_INC, SKPHI_INC, XPHI_RAH
      DOUBLE PRECISION VI1, VI2, VI3, VR1, VR2, VR3
      DOUBLE PRECISION VP1, VP2, VP3, DMOD
      DOUBLE PRECISION unit1, unit2, unit3, fact1, factor
      DOUBLE PRECISION XI1, CN1, CN2, CXI2, C2, C1, CRPER, CRPAR
      DOUBLE PRECISION TI1, TI2, TI3, TR1, TR2, TR3
      DOUBLE PRECISION PI1, PII2, PI3, PR1, PR2, PR3
      DOUBLE PRECISION PIKR, PRKI, TIKR, TRKI
      DOUBLE PRECISION E1, E2, E3, E4
      DOUBLE PRECISION CF11, CF12, CF21, CF22
      DOUBLE PRECISION AF11, AF12, AF21, AF22
      DOUBLE PRECISION C21, C22, CTTTP, CTTPT, CTTPP
      DOUBLE PRECISION CTPPT, CTPPP, CPTPP

C  Initialise

      DO O1 = 1, N_BRDF_STOKESSQ 
        BREON_VKERNEL(O1) = ZERO
      ENDDO

C#####################################################################
C#####################################################################
C  COXMUNK scalar stuff...................
C  Also removed factor of PIE in the kernel denominator
C   This makes the output exactly same as GISS model for R(1,1)
C      CKPHI_NEG = - CKPHI
c      Z = XI * XJ + SXI * SXJ * CKPHI_NEG 
c      IF ( Z .GT. ONE) Z = ONE
c      Z1 = DACOS(Z)
c      Z2 = DCOS(Z1*HALF)
C  .. Fresnel coefficients
c      REFIDX = 1.5d0
c      REFIDX_SQ = REFIDX * REFIDX
c      Z2_SQ_M1 = Z2 * Z2 + MINUS_ONE
c      H1 = REFIDX_SQ * Z2
c      H2 = DSQRT ( REFIDX_SQ + Z2_SQ_M1 )
c      RP = ( H1 - H2 ) / ( H1 + H2 )
c      RL = ( Z2 - H2 ) / ( Z2 + H2 )
C  XMP = r(i) = (1,1) reflection coefficient for natural lights
C      XMP = HALF * ( RP*RP + RL*RL )
C  Comment. This is the Lenoble SOS code
c      xind=1.50
c      xx=cos(pi*thd/360.)
c      yy=sqrt(xind*xind+xx*xx-1.0)
c      zz=xx*xind*xind
c      rl=(zz-yy)/(zz+yy)
c      rr=(xx-yy)/(xx+yy)
C  where
C     thd  = z1 * 180 / pi  (in degrees, Z1 is in radians)
C     xx   = Z2
C     yy   = H2
C     zz   = H1
C     xind = REFIDX
C     rl   = RP
C     rr   = RL
C  Coxmunk Function  
c      A = TWO * Z2
c      B = ( XI + XJ ) / A
c      IF ( B .GT. ONE ) B = ONE
c      A = PIO2 - DASIN(B)
c      TA = DTAN(A)
c      ARGUMENT = TA * TA  / PARS(1)
c      IF ( ARGUMENT .LT. CRITEXP ) THEN
c        PROB = DEXP ( - ARGUMENT )
c        FAC1 = PROB / PARS(1)
c        FAC2 = QUARTER / XI / ( B ** FOUR )
c        COXMUNK_VKERNEL(1) = XMP * FAC1 * FAC2 / XJ
c      ENDIF
C#####################################################################
C#####################################################################

C  Transcription of the RMATR subroutine from Mishchenko/Travis code.

C   CALCULATION OF THE STOKES REFLECTION MATRIX FOR
C   ILLUMINATION FROM ABOVE FOR A SURFACE SEPARATING TWO HALF-SPACES
C   WITH REFRACTIVE INDICES OF THE UPPER AND LOWER HALF-SPACES EQUAL TO
C   CN1 AND CN2, RESPECTIVELY. 
C   XI = ABS(COSINE OF THE INCIDENT ZENITH ANGLE)
C   XJ = ABS(COSINE OF THE REFLECTION ZENITH ANGLE).
C   SXI and SXJ are the respective SINES (input)
C   XPHI_REF = REFLECTION AZIMUTH ANGLE
C   BREON_VKERNEL(16-elements) = (4X4) REFLECTION MATRIX

C  For real case, incident azimuth taken to be zero

      XPHI_INC  = ZERO
      CKPHI_INC = ONE
      SKPHI_INC = ZERO

C  Check for limiting cases

      IF(DABS(XI-1D0).LT.1d-9) XI = 0.999999999999d0
      IF(DABS(XJ-1D0).LT.1d-9) XJ = 0.999999999999d0

C  help variables (coordinate transformations)

      VI1 = SXI * CKPHI_INC
      VI2 = SXI * SKPHI_INC
      VI3 = -XI
      VR1 = SXJ * CKPHI_REF
      VR2 = SXJ * SKPHI_REF
      VR3 = XJ

C    LOCAL SURFACE NORMAL FOR SPECULAR REFLECTION (normalized to 1)

      UNIT1  = VI1-VR1
      UNIT2  = VI2-VR2
      UNIT3  = VI3-VR3
      FACT1  = UNIT1*UNIT1 + UNIT2*UNIT2 + UNIT3*UNIT3
      FACTOR = DSQRT(ONE/FACT1)

C   FRESNEL REFLECTION COEFFICIENTS, assume only real for now
C   ---------------------------------------------------------

      CN1 = ONE
      CN2 = 1.5d0

C  this is the original code, but now C-variables are real

      XI1 =  FACTOR*(UNIT1*VI1+UNIT2*VI2+UNIT3*VI3)
      CXI2 = ONE - (ONE-XI1*XI1)*CN1*CN1/(CN2*CN2)
      CXI2 = DSQRT(CXI2)
      C1 = CN1*XI1
      C2 = CN2*CXI2
      CRPER = (C1-C2)/(C1+C2)
      C1 = CN2*XI1
      C2 = CN1*CXI2
      CRPAR = (C1-C2)/(C1+C2)

C  CALCULATION OF THE AMPLITUDE SCATTERING MATRIX
C  ----------------------------------------------

      TI1 = - XI * CKPHI_INC
      TI2 = - XI * SKPHI_INC
      TI3 = - SXI

      TR1 = + XJ * CKPHI_REF
      TR2 = + XJ * SKPHI_REF
      TR3 = -SXJ

      PI1  = - SKPHI_INC
      PII2 = + CKPHI_INC
      PI3  = ZERO

      PR1 = - SKPHI_REF
      PR2 = + CKPHI_REF
      PR3 = ZERO

      PIKR = PI1*VR1 + PII2*VR2 + PI3*VR3
      PRKI = PR1*VI1 + PR2*VI2  + PR3*VI3
      TIKR = TI1*VR1 + TI2*VR2  + TI3*VR3
      TRKI = TR1*VI1 + TR2*VI2  + TR3*VI3

      E1 = PIKR*PRKI
      E2 = TIKR*TRKI
      E3 = TIKR*PRKI
      E4 = PIKR*TRKI

C  Set the H-function

      HFUNCTION = 0.25d0 / ( XI * XJ )

C  Settting (1,1), (1,2), (2,1) and (2,2) components
C  -----------------------------------------------

      CF11 =  E1*CRPER+E2*CRPAR
      CF12 = -E3*CRPER+E4*CRPAR
      CF21 = -E4*CRPER+E3*CRPAR
      CF22 =  E2*CRPER+E1*CRPAR

C  Not to forget the normalization

      VP1 = VI2*VR3-VI3*VR2
      VP2 = VI3*VR1-VI1*VR3
      VP3 = VI1*VR2-VI2*VR1
      DMOD = VP1*VP1+VP2*VP2+VP3*VP3
      DMOD = DMOD*DMOD

C  if DMOD = 0, that is | n x n_0 | ^ 4 = 0 in M-T formula)
C    Then we need to set the ratio CF11 / DMOD

      IF ( DMOD .EQ. ZERO ) THEN
        CF11 = CRPAR
        CF22 = CRPER
        DMOD = 1.0d0
      ENDIF

!  Final H-function needs to be normalized

      HFUNCTION = HFUNCTION / DMOD

!  Compute the Electromagnetic contributions

      AF11 = DABS(CF11)
      AF12 = DABS(CF12)
      AF21 = DABS(CF21)
      AF22 = DABS(CF22)
      AF11 = AF11*AF11
      AF12 = AF12*AF12
      AF21 = AF21*AF21
      AF22 = AF22*AF22

      FACTOR = HALF * HFUNCTION
      BREON_VKERNEL(1) = (AF11+AF12+AF21+AF22) * FACTOR
      BREON_VKERNEL(2) = (AF11-AF12+AF21-AF22) * FACTOR
      BREON_VKERNEL(5) = (AF11-AF22+AF12-AF21) * FACTOR
      BREON_VKERNEL(6) = (AF11-AF12-AF21+AF22) * FACTOR

C  Setting (1,3), (2,3), (3,1), (3,2), (3,3) and (4,4) components
C  --------------------------------------------------------------

      C21 = CF21
      C22 = CF22
      CTTTP=CF11*CF12
      CTTPT=CF11*C21
      CTTPP=CF11*C22
      CTPPT=CF12*C21
      CTPPP=CF12*C22
      CPTPP=CF21*C22

      FACTOR = HFUNCTION
      BREON_VKERNEL(3)  =    (-CTTTP-CPTPP) * FACTOR
      BREON_VKERNEL(7)  =    (-CTTTP+CPTPP) * FACTOR
      BREON_VKERNEL(9)  =    (-CTTPT-CTPPP) * FACTOR
      BREON_VKERNEL(10) =    (-CTTPT+CTPPP) * FACTOR
      BREON_VKERNEL(11) =    ( CTTPP+CTPPT) * FACTOR
      BREON_VKERNEL(16) =    ( CTTPP-CTPPT) * FACTOR

C  Add the Diffuse term to R(1,1)
C  ------------------------------

c   IN THE PREVIOUS EQUATION, THE BRDF MODEL WE USE IS THE SINYUK-ET-AL MODEL
c   WITH FREE PARAMETERS rho,vk AND gt (see Sinyuk et al. paper).

C  This is just the Rahman Kernel.........different name !!

      XPHI_RAH = 180.0d0 - XPHI_REF
      CALL RAHMAN_VFUNCTION 
     I       ( MAXPARS, N_BRDF_STOKESSQ, NPARS, PARS,
     I         XI, SXI, XJ, SXJ,
     I         XPHI_RAH, CKPHI_REF, SKPHI_REF,
     O         RAHMAN_VKERNEL )

C  Add to the specular term

      BREON_VKERNEL(1) = BREON_VKERNEL(1) + RAHMAN_VKERNEL(1)

C  Here is the original code from France........
c    cs: cosinus of the solar zenith angle
c    cv: cosinus of the viewing zenith angle
c    phi: azimuth between the solar and observation vertical planes
c      xx=cs**(vk-1.)
c      yy=cs**(vk-1.)
c      zz=(cs+cv)**(1.-vk)
c      FF1=rho*xx*yy/zz
c      xx=sqrt(1-cs*cs)
c      yy=sqrt(1-cv*cv)
c      ww=cs*cv-xx*yy*cos(pi*phi/180.)
c      aa=1+gt*gt+2*gt*ww
c      FF2=(1-gt*gt)/(aa**1.5)
c      vv=xx/cs
c      ww=yy/cv
c      G=sqrt(vv*vv+ww*ww+2*vv*ww*cos(pi*phi/180))
c      FF3=1+(1-rho)/(1+G)
c      RBD=FF1*FF2*FF3
C RBD...........................IS THE REFLECTION COEFFICIENT

C  Finish

      RETURN
      END

C

      SUBROUTINE BPDF2008VEG_VFUNCTION 
     I       ( MAXPARS, N_BRDF_STOKESSQ, NPARS, PARS,
     I         XI, SXI, XJ, SXJ,
     I         XPHI_REF, CKPHI_REF, SKPHI_REF,
     O         BPDF2008VEG_VKERNEL )

C  include file of constants

      INCLUDE '../includes/VLIDORT.PARS'

C  Subroutine arguments

      INTEGER          MAXPARS, NPARS, N_BRDF_STOKESSQ
      DOUBLE PRECISION PARS ( MAXPARS )
      DOUBLE PRECISION XI, SXI, XJ, SXJ
      DOUBLE PRECISION XPHI_REF, CKPHI_REF, SKPHI_REF
      DOUBLE PRECISION BPDF2008VEG_VKERNEL(MAXSTOKES_SQ)

C  Local variables

      INTEGER          O1

      DOUBLE PRECISION XPHI_INC, CKPHI_INC, SKPHI_INC
      DOUBLE PRECISION VI1, VI2, VI3, VR1, VR2, VR3
      DOUBLE PRECISION VP1, VP2, VP3, DMOD
      DOUBLE PRECISION unit1, unit2, unit3, fact1, factor
      DOUBLE PRECISION XI1, CN1, CN2, CXI2, C2, C1, CRPER, CRPAR
      DOUBLE PRECISION TI1, TI2, TI3, TR1, TR2, TR3
      DOUBLE PRECISION PI1, PII2, PI3, PR1, PR2, PR3
      DOUBLE PRECISION PIKR, PRKI, TIKR, TRKI
      DOUBLE PRECISION E1, E2, E3, E4
      DOUBLE PRECISION CF11, CF12, CF21, CF22
      DOUBLE PRECISION AF11, AF12, AF21, AF22
      DOUBLE PRECISION C21, C22, CTTTP, CTTPT, CTTPP
      DOUBLE PRECISION CTPPT, CTPPP, CPTPP

C  H-function variables

      DOUBLE PRECISION HFUNCTION
      DOUBLE PRECISION ATTEN, PROJECTIONS, Z, Z1, Z2
      DOUBLE PRECISION sgamma, cgamma, calpha, calpha_sq, salpha
      DOUBLE PRECISION PLEAF, GS, GV, FP0

C  Data coefficients

      DOUBLE PRECISION PLAGIOPHILE_COEFFS(4)
      DATA PLAGIOPHILE_COEFFS
     &   /0.43181098, 0.011187479, 0.043329567, 0.19262991/
   
C  F-.M. Breon vegetation model (2009)

C  Initialise

      DO O1 = 1, N_BRDF_STOKESSQ 
        BPDF2008VEG_VKERNEL(O1) = ZERO
      ENDDO

C  Transcription of the RMATR subroutine from Mishchenko/Travis code.

C   CALCULATION OF THE STOKES REFLECTION MATRIX FOR
C   ILLUMINATION FROM ABOVE FOR A SURFACE SEPARATING TWO HALF-SPACES
C   WITH REFRACTIVE INDICES OF THE UPPER AND LOWER HALF-SPACES EQUAL TO
C   CN1 AND CN2, RESPECTIVELY.

C   XI = ABS(COSINE OF THE INCIDENT ZENITH ANGLE)
C   XJ = ABS(COSINE OF THE REFLECTION ZENITH ANGLE).
C   SXI and SXJ are the respective SINES (input)
C   XPHI_REF = REFLECTION AZIMUTH ANGLE
C   BPDF2008VEG_VKERNEL(16-elements) = (4X4) REFLECTION MATRIX

C  For real case, incident azimuth taken to be zero

      XPHI_INC  = ZERO
      CKPHI_INC = ONE
      SKPHI_INC = ZERO

C  Check for limiting cases

      IF(DABS(XI-1D0).LT.1d-9) XI = 0.999999999999d0
      IF(DABS(XJ-1D0).LT.1d-9) XJ = 0.999999999999d0

C  help variables (coordinate transformations)

      VI1 = SXI * CKPHI_INC
      VI2 = SXI * SKPHI_INC
      VI3 = -XI
      VR1 = SXJ * CKPHI_REF
      VR2 = SXJ * SKPHI_REF
      VR3 = XJ

C    LOCAL SURFACE NORMAL FOR SPECULAR REFLECTION (normalized to 1)

      UNIT1  = VI1-VR1
      UNIT2  = VI2-VR2
      UNIT3  = VI3-VR3
      FACT1  = UNIT1*UNIT1 + UNIT2*UNIT2 + UNIT3*UNIT3
      FACTOR = DSQRT(ONE/FACT1)

C   FRESNEL REFLECTION COEFFICIENTS, assume only real for now
C   ---------------------------------------------------------

      CN1 = ONE
      CN2 = 1.5d0

C  this is the original code, but now C-variables are real

      XI1 =  FACTOR*(UNIT1*VI1+UNIT2*VI2+UNIT3*VI3)
      CXI2 = ONE - (ONE-XI1*XI1)*CN1*CN1/(CN2*CN2)
      CXI2 = DSQRT(CXI2)
      C1 = CN1*XI1
      C2 = CN2*CXI2
      CRPER = (C1-C2)/(C1+C2)
      C1 = CN2*XI1
      C2 = CN1*CXI2
      CRPAR = (C1-C2)/(C1+C2)

C  CALCULATION OF THE AMPLITUDE SCATTERING MATRIX
C  ----------------------------------------------

      TI1 = - XI * CKPHI_INC
      TI2 = - XI * SKPHI_INC
      TI3 = - SXI

      TR1 = + XJ * CKPHI_REF
      TR2 = + XJ * SKPHI_REF
      TR3 = -SXJ

      PI1  = - SKPHI_INC
      PII2 = + CKPHI_INC
      PI3  = ZERO

      PR1 = - SKPHI_REF
      PR2 = + CKPHI_REF
      PR3 = ZERO

      PIKR = PI1*VR1 + PII2*VR2 + PI3*VR3
      PRKI = PR1*VI1 + PR2*VI2  + PR3*VI3
      TIKR = TI1*VR1 + TI2*VR2  + TI3*VR3
      TRKI = TR1*VI1 + TR2*VI2  + TR3*VI3

      E1 = PIKR*PRKI
      E2 = TIKR*TRKI
      E3 = TIKR*PRKI
      E4 = PIKR*TRKI

C  Set the Breon Vegetation H-function
C  ===================================

C   Angle of the surface that generates specular reflection from 
C  sun to view directions (theta)
 
!      alpha = DACOS(HALF*(mus+muv)/dcos(gamma))

      Z = XI * XJ - SXI * SXJ * CKPHI_REF   
      IF ( Z .GT. ONE) Z = ONE
      Z1 = DACOS(Z)
      Z2 = DCOS(Z1*HALF)

      calpha    = HALF * (xi + xj) / Z2  
      calpha_sq = calpha*calpha
      salpha    = dsqrt(one - calpha_sq)

! Projection of leaf surface to incident direction

      gs = PLAGIOPHILE_COEFFS(1) + xi * 
     &    (PLAGIOPHILE_COEFFS(2) + xi * 
     &    (PLAGIOPHILE_COEFFS(3) + PLAGIOPHILE_COEFFS(4)*xi))

! Projection of leaf surface to reflected direction

      gv = PLAGIOPHILE_COEFFS(1) + xj * 
     &    (PLAGIOPHILE_COEFFS(2) + xj * 
     &    (PLAGIOPHILE_COEFFS(3) + PLAGIOPHILE_COEFFS(4)*xj))
      
! Probability of leaf orientation (plagiophile distr.)

      Pleaf = 16.0d0 * calpha_sq * salpha  / pie

! Polarization model for vegetation

      PROJECTIONS =  Gv/xj + Gs/xi
      Fp0 = 0.25d0 * PLEAF / xi / xj / PROJECTIONS

! attenuation factor

      cgamma = Z2
      sgamma = dsqrt ( one - cgamma * cgamma )
      atten  = one - sgamma

!  Final H-function

      HFUNCTION = Fp0 * atten

C  Setting (1,1), (1,2), (2,1) and (2,2) components
C  ------------------------------------------------

      CF11 =  E1*CRPER+E2*CRPAR
      CF12 = -E3*CRPER+E4*CRPAR
      CF21 = -E4*CRPER+E3*CRPAR
      CF22 =  E2*CRPER+E1*CRPAR

C  Not to forget the normalization

      VP1 = VI2*VR3-VI3*VR2
      VP2 = VI3*VR1-VI1*VR3
      VP3 = VI1*VR2-VI2*VR1
      DMOD = VP1*VP1+VP2*VP2+VP3*VP3
      DMOD = DMOD*DMOD

C  if DMOD = 0, that is | n x n_0 | ^ 4 = 0 in M-T formula)
C    Then we need to set the ratio CF11 / DMOD

      IF ( DMOD .EQ. ZERO ) THEN
        CF11 = CRPAR
        CF22 = CRPER
        DMOD = 1.0d0
      ENDIF

!  Final H-function needs to be normalized

      HFUNCTION = HFUNCTION / DMOD

!  Compute the Electromagnetic contributions

      AF11 = DABS(CF11)
      AF12 = DABS(CF12)
      AF21 = DABS(CF21)
      AF22 = DABS(CF22)
      AF11 = AF11*AF11
      AF12 = AF12*AF12
      AF21 = AF21*AF21
      AF22 = AF22*AF22

      FACTOR = HALF * HFUNCTION
      BPDF2008VEG_VKERNEL(1) = (AF11+AF12+AF21+AF22) * FACTOR
      BPDF2008VEG_VKERNEL(2) = (AF11-AF12+AF21-AF22) * FACTOR
      BPDF2008VEG_VKERNEL(5) = (AF11-AF22+AF12-AF21) * FACTOR
      BPDF2008VEG_VKERNEL(6) = (AF11-AF12-AF21+AF22) * FACTOR

C  Setting (1,3), (2,3), (3,1), (3,2), (3,3) and (4,4) components
C  --------------------------------------------------------------

      C21 = CF21
      C22 = CF22
      CTTTP=CF11*CF12
      CTTPT=CF11*C21
      CTTPP=CF11*C22
      CTPPT=CF12*C21
      CTPPP=CF12*C22
      CPTPP=CF21*C22

      FACTOR = HFUNCTION
      BPDF2008VEG_VKERNEL(3)  =    (-CTTTP-CPTPP) * FACTOR
      BPDF2008VEG_VKERNEL(7)  =    (-CTTTP+CPTPP) * FACTOR
      BPDF2008VEG_VKERNEL(9)  =    (-CTTPT-CTPPP) * FACTOR
      BPDF2008VEG_VKERNEL(10) =    (-CTTPT+CTPPP) * FACTOR
      BPDF2008VEG_VKERNEL(11) =    ( CTTPP+CTPPT) * FACTOR
      BPDF2008VEG_VKERNEL(16) =    ( CTTPP-CTPPT) * FACTOR

C  Finish

      RETURN
      END

C

      SUBROUTINE BPDF2008SOIL_VFUNCTION 
     I       ( MAXPARS, N_BRDF_STOKESSQ, NPARS, PARS,
     I         XI, SXI, XJ, SXJ,
     I         XPHI_REF, CKPHI_REF, SKPHI_REF,
     O         BPDF2008SOIL_VKERNEL )

C  include file of constants

      INCLUDE '../includes/VLIDORT.PARS'

C  Subroutine arguments

      INTEGER          MAXPARS, NPARS, N_BRDF_STOKESSQ
      DOUBLE PRECISION PARS ( MAXPARS )
      DOUBLE PRECISION XI, SXI, XJ, SXJ
      DOUBLE PRECISION XPHI_REF, CKPHI_REF, SKPHI_REF
      DOUBLE PRECISION BPDF2008SOIL_VKERNEL(MAXSTOKES_SQ)

C  Local variables

      INTEGER          O1

      DOUBLE PRECISION XPHI_INC, CKPHI_INC, SKPHI_INC
      DOUBLE PRECISION VI1, VI2, VI3, VR1, VR2, VR3
      DOUBLE PRECISION VP1, VP2, VP3, DMOD
      DOUBLE PRECISION unit1, unit2, unit3, fact1, factor
      DOUBLE PRECISION XI1, CN1, CN2, CXI2, C2, C1, CRPER, CRPAR
      DOUBLE PRECISION TI1, TI2, TI3, TR1, TR2, TR3
      DOUBLE PRECISION PI1, PII2, PI3, PR1, PR2, PR3
      DOUBLE PRECISION PIKR, PRKI, TIKR, TRKI
      DOUBLE PRECISION E1, E2, E3, E4
      DOUBLE PRECISION CF11, CF12, CF21, CF22
      DOUBLE PRECISION AF11, AF12, AF21, AF22
      DOUBLE PRECISION C21, C22, CTTTP, CTTPT, CTTPP
      DOUBLE PRECISION CTPPT, CTPPP, CPTPP

C  H-function variables

      DOUBLE PRECISION HFUNCTION
      DOUBLE PRECISION ATTEN, FP0, Z, Z2, Z1
      DOUBLE PRECISION sgamma, cgamma

C  Initialise

      DO O1 = 1, N_BRDF_STOKESSQ 
        BPDF2008SOIL_VKERNEL(O1) = ZERO
      ENDDO

C  Transcription of the RMATR subroutine from Mishchenko/Travis code.

C   CALCULATION OF THE STOKES REFLECTION MATRIX FOR
C   ILLUMINATION FROM ABOVE FOR A SURFACE SEPARATING TWO HALF-SPACES
C   WITH REFRACTIVE INDICES OF THE UPPER AND LOWER HALF-SPACES EQUAL TO
C   CN1 AND CN2, RESPECTIVELY. 
C   XI = ABS(COSINE OF THE INCIDENT ZENITH ANGLE)
C   XJ = ABS(COSINE OF THE REFLECTION ZENITH ANGLE).
C   SXI and SXJ are the respective SINES (input)
C   XPHI_REF = REFLECTION AZIMUTH ANGLE
C   BPDF2008_VKERNEL(16-elements) = (4X4) REFLECTION MATRIX

C  For real case, incident azimuth taken to be zero

      XPHI_INC  = ZERO
      CKPHI_INC = ONE
      SKPHI_INC = ZERO

C  Check for limiting cases

      IF(DABS(XI-1D0).LT.1d-9) XI = 0.999999999999d0
      IF(DABS(XJ-1D0).LT.1d-9) XJ = 0.999999999999d0

C  help variables (coordinate transformations)

      VI1 = SXI * CKPHI_INC
      VI2 = SXI * SKPHI_INC
      VI3 = -XI
      VR1 = SXJ * CKPHI_REF
      VR2 = SXJ * SKPHI_REF
      VR3 = XJ

C    LOCAL SURFACE NORMAL FOR SPECULAR REFLECTION (normalized to 1)

      UNIT1  = VI1-VR1
      UNIT2  = VI2-VR2
      UNIT3  = VI3-VR3
      FACT1  = UNIT1*UNIT1 + UNIT2*UNIT2 + UNIT3*UNIT3
      FACTOR = DSQRT(ONE/FACT1)

C   FRESNEL REFLECTION COEFFICIENTS, assume only real for now
C   ---------------------------------------------------------

      CN1 = ONE
      CN2 = 1.5d0

C  this is the original code, but now C-variables are real

      XI1 =  FACTOR*(UNIT1*VI1+UNIT2*VI2+UNIT3*VI3)
      CXI2 = ONE - (ONE-XI1*XI1)*CN1*CN1/(CN2*CN2)
      CXI2 = DSQRT(CXI2)
      C1 = CN1*XI1
      C2 = CN2*CXI2
      CRPER = (C1-C2)/(C1+C2)
      C1 = CN2*XI1
      C2 = CN1*CXI2
      CRPAR = (C1-C2)/(C1+C2)

C  CALCULATION OF THE AMPLITUDE SCATTERING MATRIX
C  ----------------------------------------------

      TI1 = - XI * CKPHI_INC
      TI2 = - XI * SKPHI_INC
      TI3 = - SXI

      TR1 = + XJ * CKPHI_REF
      TR2 = + XJ * SKPHI_REF
      TR3 = -SXJ

      PI1  = - SKPHI_INC
      PII2 = + CKPHI_INC
      PI3  = ZERO

      PR1 = - SKPHI_REF
      PR2 = + CKPHI_REF
      PR3 = ZERO

      PIKR = PI1*VR1 + PII2*VR2 + PI3*VR3
      PRKI = PR1*VI1 + PR2*VI2  + PR3*VI3
      TIKR = TI1*VR1 + TI2*VR2  + TI3*VR3
      TRKI = TR1*VI1 + TR2*VI2  + TR3*VI3

      E1 = PIKR*PRKI
      E2 = TIKR*TRKI
      E3 = TIKR*PRKI
      E4 = PIKR*TRKI

C  Set the Breon Soil H-function
C  =============================

C   Angle of the surface that generates specular reflection from 
C  sun to view directions (theta)
 
!      alpha = DACOS(HALF*(mus+muv)/dcos(gamma))

      Z = XI * XJ - SXI * SXJ * CKPHI_REF   
      IF ( Z .GT. ONE) Z = ONE
      Z1 = DACOS(Z)
      Z2 = DCOS(Z1*HALF)

! Polarization model for Soil

      Fp0 = 0.25d0 / xi / xj 

! attenuation factor

      cgamma = Z2
      sgamma = dsqrt ( one - cgamma * cgamma )
      atten  = one - sgamma

!  Final H-function

      HFUNCTION = Fp0 * atten

C  Settting (1,1), (1,2), (2,1) and (2,2) components
C  -----------------------------------------------

      CF11 =  E1*CRPER+E2*CRPAR
      CF12 = -E3*CRPER+E4*CRPAR
      CF21 = -E4*CRPER+E3*CRPAR
      CF22 =  E2*CRPER+E1*CRPAR

C  Not to forget the normalization

      VP1 = VI2*VR3-VI3*VR2
      VP2 = VI3*VR1-VI1*VR3
      VP3 = VI1*VR2-VI2*VR1
      DMOD = VP1*VP1+VP2*VP2+VP3*VP3
      DMOD = DMOD*DMOD

C  if DMOD = 0, that is | n x n_0 | ^ 4 = 0 in M-T formula)
C    Then we need to set the ratio CF11 / DMOD

      IF ( DMOD .EQ. ZERO ) THEN
        CF11 = CRPAR
        CF22 = CRPER
        DMOD = 1.0d0
      ENDIF

!  Final H-function needs to be normalized

      HFUNCTION = HFUNCTION / DMOD

!  Compute the Electromagnetic contributions

      AF11 = DABS(CF11)
      AF12 = DABS(CF12)
      AF21 = DABS(CF21)
      AF22 = DABS(CF22)
      AF11 = AF11*AF11
      AF12 = AF12*AF12
      AF21 = AF21*AF21
      AF22 = AF22*AF22

      FACTOR = HALF * HFUNCTION
      BPDF2008SOIL_VKERNEL(1) = (AF11+AF12+AF21+AF22) * FACTOR
      BPDF2008SOIL_VKERNEL(2) = (AF11-AF12+AF21-AF22) * FACTOR
      BPDF2008SOIL_VKERNEL(5) = (AF11-AF22+AF12-AF21) * FACTOR
      BPDF2008SOIL_VKERNEL(6) = (AF11-AF12-AF21+AF22) * FACTOR

C  Setting (1,3), (2,3), (3,1), (3,2), (3,3) and (4,4) components
C  --------------------------------------------------------------

      C21 = CF21
      C22 = CF22
      CTTTP=CF11*CF12
      CTTPT=CF11*C21
      CTTPP=CF11*C22
      CTPPT=CF12*C21
      CTPPP=CF12*C22
      CPTPP=CF21*C22

      FACTOR = HFUNCTION
      BPDF2008SOIL_VKERNEL(3)  =    (-CTTTP-CPTPP) * FACTOR
      BPDF2008SOIL_VKERNEL(7)  =    (-CTTTP+CPTPP) * FACTOR
      BPDF2008SOIL_VKERNEL(9)  =    (-CTTPT-CTPPP) * FACTOR
      BPDF2008SOIL_VKERNEL(10) =    (-CTTPT+CTPPP) * FACTOR
      BPDF2008SOIL_VKERNEL(11) =    ( CTTPP+CTPPT) * FACTOR
      BPDF2008SOIL_VKERNEL(16) =    ( CTTPP-CTPPT) * FACTOR

C  Finish

      RETURN
      END

C

      SUBROUTINE BPDF2009_VFUNCTION 
     I       ( MAXPARS, N_BRDF_STOKESSQ, NPARS, PARS,
     I         XI, SXI, XJ, SXJ,
     I         XPHI_REF, CKPHI_REF, SKPHI_REF,
     O         BPDF2009_VKERNEL )

C  include file of constants

      INCLUDE '../includes/VLIDORT.PARS'

C  Subroutine arguments

      INTEGER          MAXPARS, NPARS, N_BRDF_STOKESSQ
      DOUBLE PRECISION PARS ( MAXPARS )
      DOUBLE PRECISION XI, SXI, XJ, SXJ
      DOUBLE PRECISION XPHI_REF, CKPHI_REF, SKPHI_REF
      DOUBLE PRECISION BPDF2009_VKERNEL(MAXSTOKES_SQ)

C  Local variables

      INTEGER          O1

      DOUBLE PRECISION XPHI_INC, CKPHI_INC, SKPHI_INC
      DOUBLE PRECISION VI1, VI2, VI3, VR1, VR2, VR3
      DOUBLE PRECISION VP1, VP2, VP3, DMOD
      DOUBLE PRECISION unit1, unit2, unit3, fact1, factor
      DOUBLE PRECISION XI1, CN1, CN2, CXI2, C2, C1, CRPER, CRPAR
      DOUBLE PRECISION TI1, TI2, TI3, TR1, TR2, TR3
      DOUBLE PRECISION PI1, PII2, PI3, PR1, PR2, PR3
      DOUBLE PRECISION PIKR, PRKI, TIKR, TRKI
      DOUBLE PRECISION E1, E2, E3, E4
      DOUBLE PRECISION CF11, CF12, CF21, CF22
      DOUBLE PRECISION AF11, AF12, AF21, AF22
      DOUBLE PRECISION C21, C22, CTTTP, CTTPT, CTTPP
      DOUBLE PRECISION CTPPT, CTPPP, CPTPP

C  H-function variables

      DOUBLE PRECISION HFUNCTION, NDVI, DEXPNDVI
      DOUBLE PRECISION ATTEN, FP0, Z, Z2, Z1
      DOUBLE PRECISION sgamma, cgamma

C  Initialise

      DO O1 = 1, N_BRDF_STOKESSQ 
        BPDF2009_VKERNEL(O1) = ZERO
      ENDDO

C  Transcription of the RMATR subroutine from Mishchenko/Travis code.

C   CALCULATION OF THE STOKES REFLECTION MATRIX FOR
C   ILLUMINATION FROM ABOVE FOR A SURFACE SEPARATING TWO HALF-SPACES
C   WITH REFRACTIVE INDICES OF THE UPPER AND LOWER HALF-SPACES EQUAL TO
C   CN1 AND CN2, RESPECTIVELY. 
C   XI = ABS(COSINE OF THE INCIDENT ZENITH ANGLE)
C   XJ = ABS(COSINE OF THE REFLECTION ZENITH ANGLE).
C   SXI and SXJ are the respective SINES (input)
C   XPHI_REF = REFLECTION AZIMUTH ANGLE
C   BPDF2008_VKERNEL(16-elements) = (4X4) REFLECTION MATRIX

C  For real case, incident azimuth taken to be zero

      XPHI_INC  = ZERO
      CKPHI_INC = ONE
      SKPHI_INC = ZERO

C  Check for limiting cases

      IF(DABS(XI-1D0).LT.1d-9) XI = 0.999999999999d0
      IF(DABS(XJ-1D0).LT.1d-9) XJ = 0.999999999999d0

C  help variables (coordinate transformations)

      VI1 = SXI * CKPHI_INC
      VI2 = SXI * SKPHI_INC
      VI3 = -XI
      VR1 = SXJ * CKPHI_REF
      VR2 = SXJ * SKPHI_REF
      VR3 = XJ

C    LOCAL SURFACE NORMAL FOR SPECULAR REFLECTION (normalized to 1)

      UNIT1  = VI1-VR1
      UNIT2  = VI2-VR2
      UNIT3  = VI3-VR3
      FACT1  = UNIT1*UNIT1 + UNIT2*UNIT2 + UNIT3*UNIT3
      FACTOR = DSQRT(ONE/FACT1)

C   FRESNEL REFLECTION COEFFICIENTS, assume only real for now
C   ---------------------------------------------------------

      CN1 = ONE
      CN2 = 1.5d0

C  this is the original code, but now C-variables are real

      XI1 =  FACTOR*(UNIT1*VI1+UNIT2*VI2+UNIT3*VI3)
      CXI2 = ONE - (ONE-XI1*XI1)*CN1*CN1/(CN2*CN2)
      CXI2 = DSQRT(CXI2)
      C1 = CN1*XI1
      C2 = CN2*CXI2
      CRPER = (C1-C2)/(C1+C2)
      C1 = CN2*XI1
      C2 = CN1*CXI2
      CRPAR = (C1-C2)/(C1+C2)

C  CALCULATION OF THE AMPLITUDE SCATTERING MATRIX
C  ----------------------------------------------

      TI1 = - XI * CKPHI_INC
      TI2 = - XI * SKPHI_INC
      TI3 = - SXI

      TR1 = + XJ * CKPHI_REF
      TR2 = + XJ * SKPHI_REF
      TR3 = -SXJ

      PI1  = - SKPHI_INC
      PII2 = + CKPHI_INC
      PI3  = ZERO

      PR1 = - SKPHI_REF
      PR2 = + CKPHI_REF
      PR3 = ZERO

      PIKR = PI1*VR1 + PII2*VR2 + PI3*VR3
      PRKI = PR1*VI1 + PR2*VI2  + PR3*VI3
      TIKR = TI1*VR1 + TI2*VR2  + TI3*VR3
      TRKI = TR1*VI1 + TR2*VI2  + TR3*VI3

      E1 = PIKR*PRKI
      E2 = TIKR*TRKI
      E3 = TIKR*PRKI
      E4 = PIKR*TRKI

C  Set the Breon 2009 
C  =============================

C   Angle of the surface that generates specular reflection from 
C  sun to view directions (theta)
 
!      alpha = DACOS(HALF*(mus+muv)/dcos(gamma))

      Z = XI * XJ - SXI * SXJ * CKPHI_REF   
      IF ( Z .GT. ONE) Z = ONE
      Z1 = DACOS(Z)
      Z2 = DCOS(Z1*HALF)

!  Exponential of the NDVI
!    Out of range values default to zero

      NDVI = PARS(1)
      IF ( NDVI .GT. 1.0d0 .or. NDVI .lt. -1.0d0 ) THEN
        NDVI = 0.0d0
      ENDIF
      DEXPNDVI = DEXP ( - NDVI )

! Polarization model is isotropic

      Fp0 = 0.25d0 * DEXPNDVI / ( xi + xj )

! attenuation factor

      cgamma = Z2
      sgamma = dsqrt ( one - cgamma * cgamma )
      atten  = dexp ( - sgamma / cgamma )

!  Final H-function

      HFUNCTION = Fp0 * atten

C  Settting (1,1), (1,2), (2,1) and (2,2) components
C  -----------------------------------------------

      CF11 =  E1*CRPER+E2*CRPAR
      CF12 = -E3*CRPER+E4*CRPAR
      CF21 = -E4*CRPER+E3*CRPAR
      CF22 =  E2*CRPER+E1*CRPAR

C  Not to forget the normalization

      VP1 = VI2*VR3-VI3*VR2
      VP2 = VI3*VR1-VI1*VR3
      VP3 = VI1*VR2-VI2*VR1
      DMOD = VP1*VP1+VP2*VP2+VP3*VP3
      DMOD = DMOD*DMOD

C  if DMOD = 0, that is | n x n_0 | ^ 4 = 0 in M-T formula)
C    Then we need to set the ratio CF11 / DMOD

      IF ( DMOD .EQ. ZERO ) THEN
        CF11 = CRPAR
        CF22 = CRPER
        DMOD = 1.0d0
      ENDIF

!  Final H-function needs to be normalized

      HFUNCTION = HFUNCTION / DMOD

!  Compute the Electromagnetic contributions

      AF11 = DABS(CF11)
      AF12 = DABS(CF12)
      AF21 = DABS(CF21)
      AF22 = DABS(CF22)
      AF11 = AF11*AF11
      AF12 = AF12*AF12
      AF21 = AF21*AF21
      AF22 = AF22*AF22

      FACTOR = HALF * HFUNCTION
      BPDF2009_VKERNEL(1) = (AF11+AF12+AF21+AF22) * FACTOR
      BPDF2009_VKERNEL(2) = (AF11-AF12+AF21-AF22) * FACTOR
      BPDF2009_VKERNEL(5) = (AF11-AF22+AF12-AF21) * FACTOR
      BPDF2009_VKERNEL(6) = (AF11-AF12-AF21+AF22) * FACTOR

C  Setting (1,3), (2,3), (3,1), (3,2), (3,3) and (4,4) components
C  --------------------------------------------------------------

      C21 = CF21
      C22 = CF22
      CTTTP=CF11*CF12
      CTTPT=CF11*C21
      CTTPP=CF11*C22
      CTPPT=CF12*C21
      CTPPP=CF12*C22
      CPTPP=CF21*C22

      FACTOR = HFUNCTION
      BPDF2009_VKERNEL(3)  =    (-CTTTP-CPTPP) * FACTOR
      BPDF2009_VKERNEL(7)  =    (-CTTTP+CPTPP) * FACTOR
      BPDF2009_VKERNEL(9)  =    (-CTTPT-CTPPP) * FACTOR
      BPDF2009_VKERNEL(10) =    (-CTTPT+CTPPP) * FACTOR
      BPDF2009_VKERNEL(11) =    ( CTTPP+CTPPT) * FACTOR
      BPDF2009_VKERNEL(16) =    ( CTTPP-CTPPT) * FACTOR

C  Finish

      RETURN
      END

        double precision function derfc(x)

        double precision x

c Returns the complementary error function erfc(x) with fractional error 
c everywhere less than 1.2 * 10^7.

        double precision t,z

        z = dabs(x)
        t = 1.d0/(1.d0+0.5d0*z)
        derfc = t*dexp(-z*z-1.26551223d0+t*(1.00002368d0+t*(.37409196d0+
     * t*(.09678418d0+t*(-.18628806d0+t*(.27886807d0+t*(-1.13520398d0+
     * t*(1.48851587d0+t*(-.82215223d0+t*.17087277d0)))))))))
        if (x .lt. 0.d0) derfc = 2.d0-derfc

        return
        end
