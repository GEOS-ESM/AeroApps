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
! #                   2.5, 2.6                                  #
! #  Release Date :   December 2005  (2.0)                      #
! #  Release Date :   March 2007     (2.2)                      #
! #  Release Date :   October 2007   (2.3)                      #
! #  Release Date :   December 2008  (2.4)                      #
! #  Release Date :   April 2009     (2.4R)                     #
! #  Release Date :   July 2009      (2.4RT)                    #
! #  Release Date :   October 2010   (2.4RTC)                   #
! #  Release Date :   March 2011     (2.5)                      #
! #  Release Date :   May 2012       (2.6)                      #
! #                                                             #
! #       NEW: TOTAL COLUMN JACOBIANS         (2.4)             #
! #       NEW: BPDF Land-surface KERNELS      (2.4R)            #
! #       NEW: Thermal Emission Treatment     (2.4RT)           #
! #       Consolidated BRDF treatment         (2.4RTC)          #
! #       f77/f90 Release                     (2.5)             #
! #       External SS / New I/O Structures    (2.6)             #
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
! # Subroutines in this Module                                  #
! #                                                             #
! #            LISPARSE_VFUNCTION_PLUS                          #
! #            LIDENSE_VFUNCTION_PLUS                           #
! #            HAPKE_VFUNCTION_PLUS                             #
! #            RAHMAN_VFUNCTION_PLUS                            #
! #            COXMUNK_VFUNCTION_PLUS                           #
! #            GISSCOXMUNK_VFUNCTION_PLUS                       #
! #            COXMUNK_VFUNCTION_DB_PLUS                        #
! #            GISSCOXMUNK_VFUNCTION_DB_PLUS                    #
! #                                                             #
! ###############################################################

      MODULE vbrdf_LinSup_kernels_m

      use VLIDORT_PARS
      use vbrdf_sup_aux_m, only : derfc_e

      PRIVATE
      PUBLIC :: LISPARSE_VFUNCTION_PLUS, &
                LIDENSE_VFUNCTION_PLUS, &
                HAPKE_VFUNCTION_PLUS, &
                RAHMAN_VFUNCTION_PLUS, &
                COXMUNK_VFUNCTION_PLUS, &
                GISSCOXMUNK_VFUNCTION_PLUS, &
                COXMUNK_VFUNCTION_DB_PLUS, &
                GISSCOXMUNK_VFUNCTION_DB_PLUS

      CONTAINS


      SUBROUTINE LIDENSE_VFUNCTION_PLUS &
         ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, N_BRDF_STOKESSQ, &
           XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
           LIDENSE_VKERNEL, LIDENSE_VDERIVATIVES )

!  include file of constants

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER ::          MAXPARS, NPARS, N_BRDF_STOKESSQ
      DOUBLE PRECISION :: PARS ( MAXPARS )
      LOGICAL ::          DO_DERIV_PARS ( MAXPARS )
      DOUBLE PRECISION :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI

      DOUBLE PRECISION :: LIDENSE_VKERNEL ( MAXSTOKES_SQ )
      DOUBLE PRECISION :: LIDENSE_VDERIVATIVES ( MAXPARS, MAXSTOKES_SQ )

!  local variables

      DOUBLE PRECISION :: X_INC, X_REF, SX_INC, SX_REF, ANG_P, TX
      DOUBLE PRECISION :: T_INC, T_REF, T_INC_SQ, T_REF_SQ
      DOUBLE PRECISION :: CKSI, DELTA, T, COST, SINT, DSQ, SINTCOST
      DOUBLE PRECISION :: A, B, H, R, P, Q, DT1, DT2, DT2SQ, P_QR
      DOUBLE PRECISION :: A2, R2, DX_H, DX_R, DX_Q, DX_P, DX_P_QR, DY_Q
      DOUBLE PRECISION :: XPHI, CKPHI

!  Initialise
!    -- Return for special case

      LIDENSE_VKERNEL      = ZERO
      LIDENSE_VDERIVATIVES = ZERO
      XPHI  = PIE - PHI
      CKPHI = - CPHI
      IF ( ( XI .EQ. XJ ) .AND. ( CKPHI.EQ.ONE ) ) RETURN

!  Function
!  ========

!  .. incidence

      TX       = SXJ / XJ
      T_INC    = PARS(2) * TX
      T_INC_SQ = T_INC * T_INC
      ANG_P    = DATAN ( T_INC )
      X_INC    = DCOS(ANG_P)
      SX_INC   = DSIN(ANG_P)

!  .. reflection

      TX       = SXI / XI
      T_REF    = PARS(2) * TX
      T_REF_SQ = T_REF * T_REF
      ANG_P    = DATAN ( T_REF )
      X_REF    = DCOS(ANG_P)
      SX_REF   = DSIN(ANG_P)

!  ksi cosine

      CKSI = X_INC  * X_REF + SX_INC * SX_REF * CKPHI

!  contributions P and R and derivatives (if flagged)

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

!  evaluate cos(t)

      DT1   = T_REF_SQ + T_INC_SQ
      DT2   = T_INC * T_REF
      DT2SQ = DT2 * DT2
      DELTA = DSQRT ( DT1 - TWO * DT2 * CKPHI )
      DSQ   = DELTA * DELTA
      H     = DSQRT ( DSQ + SKPHI * SKPHI * DT2SQ )
      COST  = PARS(1) * H / R

!  set Q function and its derivatives if flagged

      Q    = ONE
      DX_Q = ZERO
      DY_Q = ZERO
      IF ( COST .LE. ONE ) THEN
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

!  set the kernel
!  --------------

      P_QR = P / Q / R
      LIDENSE_VKERNEL(1) = P_QR - TWO

!  Other vector functions

!    P L A C E H O L D E R

!  Scalar derivatives (1,1) matrix entry
!  -------------------------------------

      IF ( DO_DERIV_PARS(1) ) THEN
        LIDENSE_VDERIVATIVES(1,1) = - P_QR * DY_Q / Q
      ENDIF

      IF ( DO_DERIV_PARS(2) ) THEN
        DX_P_QR = ( DX_P / P ) - ( DX_R / R ) - ( DX_Q / Q )
        LIDENSE_VDERIVATIVES(2,1) = P_QR * DX_P_QR
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LIDENSE_VFUNCTION_PLUS

!^C

      SUBROUTINE LISPARSE_VFUNCTION_PLUS &
         ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, N_BRDF_STOKESSQ, &
           XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
           LISPARSE_VKERNEL, LISPARSE_VDERIVATIVES )

!  include file of constants

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER ::          MAXPARS, NPARS, N_BRDF_STOKESSQ
      DOUBLE PRECISION :: PARS ( MAXPARS )
      LOGICAL ::          DO_DERIV_PARS ( MAXPARS )
      DOUBLE PRECISION :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI

      DOUBLE PRECISION :: LISPARSE_VKERNEL      ( MAXSTOKES_SQ )
      DOUBLE PRECISION :: LISPARSE_VDERIVATIVES ( MAXPARS, MAXSTOKES_SQ  )

!  local variables

      DOUBLE PRECISION :: X_INC, X_REF, SX_INC, SX_REF, ANG_P, TX
      DOUBLE PRECISION :: T_INC, T_REF, T_INC_SQ, T_REF_SQ
      DOUBLE PRECISION :: CKSI, DELTA, T, COST, SINT, DSQ, SINTCOST
      DOUBLE PRECISION :: A, B, H, R, P, Q, DT1, DT2, DT2SQ, QR
      DOUBLE PRECISION :: A2, R2, DX_H, DX_R, DX_Q, DX_P, DX_QR, DY_Q
      DOUBLE PRECISION :: XPHI, CKPHI

!  Initialise
!    -- Return for special case

      LISPARSE_VKERNEL      = ZERO
      LISPARSE_VDERIVATIVES = ZERO
      XPHI  = PIE - PHI
      CKPHI = - CPHI
      IF ( ( XI .EQ. XJ ) .AND. ( CKPHI.EQ.ONE ) ) RETURN

!  Function
!  ========

!  .. incidence

      TX       = SXJ / XJ
      T_INC    = PARS(2) * TX
      T_INC_SQ = T_INC * T_INC
      ANG_P    = DATAN ( T_INC )
      X_INC    = DCOS(ANG_P)
      SX_INC   = DSIN(ANG_P)

!  .. reflection

      TX       = SXI / XI
      T_REF    = PARS(2) * TX
      T_REF_SQ = T_REF * T_REF
      ANG_P    = DATAN ( T_REF )
      X_REF    = DCOS(ANG_P)
      SX_REF   = DSIN(ANG_P)

!  ksi cosine

      CKSI = X_INC  * X_REF + SX_INC * SX_REF * CKPHI

!  contributions P and R and derivatives (if flagged)

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

!  evaluate cos(t)

      DT1   = T_REF_SQ + T_INC_SQ
      DT2   = T_INC * T_REF
      DT2SQ = DT2 * DT2
      DELTA = DSQRT ( DT1 - TWO * DT2 * CKPHI )
      DSQ   = DELTA * DELTA
      H     = DSQRT ( DSQ + SKPHI * SKPHI * DT2SQ )
      COST  = PARS(1) * H / R

!  set Q function and its derivatives if flagged

      Q = ONE
      DX_Q = ZERO
      DY_Q = ZERO
      IF ( COST .LE. ONE ) THEN
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

!  set the kernel
!  --------------

      QR = Q * R
      LISPARSE_VKERNEL(1) = HALF * P - QR

!  Other vector functions

!    P L A C E H O L D E R

!  Scalar derivatives (1,1) matrix entry
!  -------------------------------------

      IF ( DO_DERIV_PARS(1) ) THEN
        LISPARSE_VDERIVATIVES(1,1) = - R * DY_Q
      ENDIF

      IF ( DO_DERIV_PARS(2) ) THEN
        DX_QR = DX_R * Q + DX_Q * R
        LISPARSE_VDERIVATIVES(2,1) = HALF * DX_P - DX_QR
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LISPARSE_VFUNCTION_PLUS

!^C

      SUBROUTINE HAPKE_VFUNCTION_PLUS &
         ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, N_BRDF_STOKESSQ, &
           XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
           HAPKE_VKERNEL, HAPKE_VDERIVATIVES )

!  include file of constants

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Subroutine arguments

!  Subroutine arguments

      INTEGER ::          MAXPARS, NPARS, N_BRDF_STOKESSQ
      DOUBLE PRECISION :: PARS ( MAXPARS )
      LOGICAL          :: DO_DERIV_PARS ( MAXPARS )
      DOUBLE PRECISION :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI

      DOUBLE PRECISION :: HAPKE_VKERNEL      ( MAXSTOKES_SQ )
      DOUBLE PRECISION :: HAPKE_VDERIVATIVES ( MAXPARS, MAXSTOKES_SQ )

!  Hapke Kernel function.
!    - New version, Fresh Coding
!    - Old version uses DISORT code; for validation.

!  input variables:

!    XI, SXI  : Cosine/Sine of angle of reflection (positive)
!    XJ, SXJ  : Cosine/Sine of angle of incidence (positive)
!    XPHI     : Difference of azimuth angles of incidence and reflection
!    PARS(1)  : single scattering albedo in Hapke's BDR model
!    PARS(2)  : angular width parameter of opposition effect in Hapke's
!    PARS(3)  : Empirical hot spot multiplier

!  local variables
!    B0_EMPIR : empirical factor to account for the finite size of
!               particles in Hapke's BDR model
!    B_HOT    : term that accounts for the opposition effect
!               (retroreflectance, hot spot) in Hapke's BDR model
!    CTHETA   : cosine of phase angle in Hapke's BDR model
!    GAMMA    : albedo factor in Hapke's BDR model
!    PHASE    : scattering phase function in Hapke's BDR model
!    THETA  : phase angle (radians); the angle between incidence and
!             reflection directions in Hapke's BDR model

!  Local variables

      DOUBLE PRECISION :: CTHETA, THETA, PHASE
      DOUBLE PRECISION :: HOTSPOT, B0_EMPIR, HELP_HOT, B_HOT
      DOUBLE PRECISION :: SSALBEDO, GAMMA, REFLEC, FUNCT
      DOUBLE PRECISION :: HELP_J, GHELP_J, TERM_J
      DOUBLE PRECISION :: HELP_I, GHELP_I, TERM_I
      DOUBLE PRECISION :: TI_TJ, DT1, DT2
      DOUBLE PRECISION :: XPHI, CKPHI

!  Initialise

      HAPKE_VKERNEL      = ZERO
      HAPKE_VDERIVATIVES = ZERO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  (1,1) function (scalar form)

!  geometrical part

!  This is the code that is in DISORT - not right, I think.
!        CTHETA = XI * XJ + DABS(SXI) *  DABS(SXJ) * CKPHI

      CTHETA = XI * XJ + SXI * SXJ * CKPHI
      THETA  = DACOS( CTHETA )
      PHASE  = ONE + HALF * CTHETA

!  hot spot parameterization

      HOTSPOT  = PARS(2)
      B0_EMPIR = PARS(3)
      HELP_HOT = HOTSPOT + DTAN ( HALF * THETA )
      B_HOT    = B0_EMPIR * HOTSPOT / HELP_HOT

!  Albedo parameterization

      SSALBEDO = PARS(1)
      GAMMA    = DSQRT ( ONE - SSALBEDO )
      HELP_J   = TWO * XJ
      GHELP_J  = ( ONE + HELP_J * GAMMA )
      TERM_J   = ( ONE + HELP_J ) / GHELP_J
      HELP_I   = TWO * XI
      GHELP_I  = ( ONE + HELP_I * GAMMA )
      TERM_I   = ( ONE + HELP_I ) / GHELP_I
      TI_TJ    = TERM_J * TERM_I

!  Function

      REFLEC           = SSALBEDO * QUARTER / ( XI + XJ )
      FUNCT            = ( ONE + B_HOT ) * PHASE + TI_TJ - ONE
      HAPKE_VKERNEL(1) = REFLEC * FUNCT

!  Other vector functions

!    P L A C E H O L D E R

!  Scalar derivatives (1,1) matrix entry
!  -------------------------------------

!  ssalbedo derivative

      IF ( DO_DERIV_PARS(1) ) THEN
        DT1 = HAPKE_VKERNEL(1) / SSALBEDO
        DT2 = ( HELP_J / GHELP_J ) + ( HELP_I / GHELP_I )
        DT2 = DT2 * TI_TJ * HALF / GAMMA
        HAPKE_VDERIVATIVES(1,1) = DT1 + DT2 * REFLEC
      ENDIF

!  Hotspot  derivative

      IF ( DO_DERIV_PARS(2) ) THEN
        DT1 = B_HOT * ( B0_EMPIR - B_HOT ) / B0_EMPIR / HOTSPOT
        HAPKE_VDERIVATIVES(2,1) = DT1 * REFLEC * PHASE
      ENDIF

!  empirical factor derivative

      IF ( DO_DERIV_PARS(3) ) THEN
        DT1 = B_HOT / B0_EMPIR
        HAPKE_VDERIVATIVES(3,1) = DT1 * REFLEC * PHASE
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE HAPKE_VFUNCTION_PLUS


!^C

      SUBROUTINE RAHMAN_VFUNCTION_PLUS &
         ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, N_BRDF_STOKESSQ, &
           XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
           RAHMAN_VKERNEL, RAHMAN_VDERIVATIVES )

!  include file of constants

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER ::          MAXPARS, NPARS, N_BRDF_STOKESSQ
      DOUBLE PRECISION :: PARS ( MAXPARS )
      LOGICAL ::          DO_DERIV_PARS ( MAXPARS )
      DOUBLE PRECISION :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI

      DOUBLE PRECISION :: RAHMAN_VKERNEL      ( MAXSTOKES_SQ )
      DOUBLE PRECISION :: RAHMAN_VDERIVATIVES ( MAXPARS, MAXSTOKES_SQ )

!  local variables

      DOUBLE PRECISION :: T_INC, T_REF, DT1, DT2
      DOUBLE PRECISION :: CXI, DELTA, K1_SQ, FACT, K0, K1, K2
      DOUBLE PRECISION :: HELPM, HELPR, HELPG, D_HELPM, D_FACT
      DOUBLE PRECISION :: GEOM, PHASE, RFAC, RFAC1, D_K0, D_K1, D_K2
      DOUBLE PRECISION :: XPHI, CKPHI

!  Initialise

      RAHMAN_VKERNEL      = ZERO
      RAHMAN_VDERIVATIVES = ZERO
      XPHI  = PIE - PHI
      CKPHI = - CPHI
      IF ( XI.EQ.ZERO .OR. XJ.EQ.ZERO ) RETURN

!  parameters

      K0 = PARS(1)
      K1 = PARS(2)
      K2 = PARS(3)

!  geometrical angle xi

      CXI = XI * XJ + SXI * SXJ * CKPHI
      IF ( CXI .GT. ONE ) CXI = ONE

!  Phase function

      K1_SQ = K1 * K1
      HELPM = ONE - K1_SQ
      FACT  = ONE + K1_SQ + TWO * K1 * CXI
      PHASE = HELPM / ( FACT ** ONEP5 )

!  Delta and R-factor

      T_INC = SXI / XI
      T_REF = SXJ / XJ
      DT1   = T_INC*T_INC + T_REF*T_REF
      DT2   = T_INC * T_REF
      DELTA = DSQRT ( DT1 - TWO * DT2 * CKPHI )
      HELPR = ONE / ( ONE + DELTA )
      RFAC  = ( ONE - K0 ) * HELPR
      RFAC1 = ONE + RFAC

!  Geom factor and kernel

      HELPG = XI * XJ * ( XI + XJ )
      GEOM  = HELPG ** ( K2 - ONE)
      RAHMAN_VKERNEL(1) = K0 * PHASE * RFAC1 * GEOM

!  Other vector functions

!    P L A C E H O L D E R

!  Scalar derivatives (1,1) matrix entry
!  -------------------------------------

!  K0 derivative

      IF ( DO_DERIV_PARS(1) ) THEN
        D_K0   = ( ONE / K0 ) - ( HELPR / RFAC1 )
        RAHMAN_VDERIVATIVES(1,1) = RAHMAN_VKERNEL(1) * D_K0
      ENDIF

!  Phase function derivative

      IF ( DO_DERIV_PARS(2) ) THEN
        D_FACT  =   TWO * K1 + TWO * CXI
        D_HELPM = - TWO * K1
        D_K1    = ( D_HELPM / HELPM ) - ONEP5 * ( D_FACT / FACT )
        RAHMAN_VDERIVATIVES(2,1) = RAHMAN_VKERNEL(1) * D_K1
      ENDIF

!  K2 derivative

      IF ( DO_DERIV_PARS(3) ) THEN
        D_K2 = DLOG ( HELPG )
        RAHMAN_VDERIVATIVES(3,1) = RAHMAN_VKERNEL(1) * D_K2
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE RAHMAN_VFUNCTION_PLUS

!

      SUBROUTINE COXMUNK_VFUNCTION_PLUS &
         ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, N_BRDF_STOKESSQ, &
           XJ, SXJ, XI, SXI, PHI, CKPHI, SKPHI, &
           COXMUNK_VKERNEL, COXMUNK_VDERIVATIVES )

!  include file of constants

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER ::          MAXPARS, NPARS, N_BRDF_STOKESSQ
      DOUBLE PRECISION :: PARS ( MAXPARS )
      LOGICAL ::          DO_DERIV_PARS ( MAXPARS )
      DOUBLE PRECISION :: XI, SXI, XJ, SXJ, PHI, CKPHI, SKPHI

      DOUBLE PRECISION :: COXMUNK_VKERNEL      ( MAXSTOKES_SQ )
      DOUBLE PRECISION :: COXMUNK_VDERIVATIVES ( MAXPARS, MAXSTOKES_SQ )

!  Critical exponent taken out

      DOUBLE PRECISION, PARAMETER :: CRITEXP = 88.0D0

!  Local variables

      DOUBLE PRECISION :: Z, Z1, Z2, Z2_SQ_M1, H1, H2, RP, RL, XMP
      DOUBLE PRECISION :: A, B, TA, ARGUMENT, PROB, FAC1, FAC2, CKPHI_NEG
      DOUBLE PRECISION :: S1, S2, S3, XXI, XXJ
      DOUBLE PRECISION :: T1_I, T2_I, DCOT_I
      DOUBLE PRECISION :: T1_R, T2_R, DCOT_R
      DOUBLE PRECISION :: SHADOWI, SHADOWR, SHADOW

      DOUBLE PRECISION :: H1H2, H2Z2, TA_SQ, DFAC2, DH1, DH2, DRP, DRL
      DOUBLE PRECISION :: D_S1, D_S2, D_T1, D_T2
      DOUBLE PRECISION :: D_SHADOWI, D_SHADOWR, D_SHADOW

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!               Remark on Use of shadow effect
!               ------------------------------
!  Shadow effect is controlled by the third parameter. That is, if
!  PARS(3) not equal to then shadow effect will be included.
!    --- NPARS should always be 3 for this Kernel.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Initialise

      COXMUNK_VKERNEL      = ZERO
      COXMUNK_VDERIVATIVES = ZERO

!  Comment. 18 January 2006.
!  We have found in comparisons with the Giss Cox-Munk code that
!  the input COSPHI (CKPHI) here is the negative of what we actually nee
!  so introduce local variable which takes care of this

!  Also removed factor of PIE in the kernel denominator
!   This makes the output exactly same as GISS model for R(1,1)

      CKPHI_NEG = - CKPHI

!  (1,1) function (scalar form)

!  ..Scatter angles

! old   Z = - XI * XJ + SXI * SXJ * CKPHI
! old   IF ( Z .LT. MINUS_ONE) Z = MINUS_ONE
! old   Z1 = DACOS(-Z)
! old   Z2 = DCOS(Z1*HALF)

      Z = XI * XJ + SXI * SXJ * CKPHI_NEG
      IF ( Z .GT. ONE) Z = ONE
      Z1 = DACOS(Z)
      Z2 = DCOS(Z1*HALF)

!  .. Fresnel coefficients

      Z2_SQ_M1 = Z2 * Z2 + MINUS_ONE
      H1 = PARS(2) * Z2
      H2 = DSQRT ( PARS(2) + Z2_SQ_M1 )
      H1H2 = H1 + H2
      RP = ( H1 - H2 ) / H1H2
      H2Z2 = Z2 + H2
      RL = ( Z2 - H2 ) / H2Z2
      XMP = HALF * ( RP*RP + RL*RL )

!  Coxmunk Function

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
        COXMUNK_VKERNEL(1) = XMP * FAC1 * FAC2 / XJ
      ENDIF

!  inverse slope-squared derivative (regular term)

      IF ( DO_DERIV_PARS(1) ) THEN
        IF ( ARGUMENT .LT. CRITEXP ) THEN
          DFAC2 = ( PARS(1) - TA_SQ ) / PARS(1) / PARS(1)
          COXMUNK_VDERIVATIVES(1,1) = - COXMUNK_VKERNEL(1) * DFAC2
        ENDIF
      ENDIF

!  No Shadow code if not flagged

      IF ( PARS(3) .EQ. ZERO ) RETURN

!  Shadow code

      S1 = DSQRT(PARS(1)/PIE)
      S3 = ONE/(DSQRT(PARS(1)))
      S2 = S3*S3

      XXI  = XI*XI
      DCOT_I = XI/DSQRT(ONE-XXI)
      T1_I   = DEXP(-DCOT_I*DCOT_I*S2)
      T2_I   = DERFC_E(DCOT_I*S3)
      SHADOWI = HALF * ( S1*T1_I/DCOT_I - T2_I )

      XXJ  = XJ*XJ
      DCOT_R = XJ/DSQRT(ONE-XXJ)
      T1_R   = DEXP(-DCOT_R*DCOT_R*S2)
      T2_R   = DERFC_E(DCOT_R*S3)
      SHADOWR = HALF * ( S1*T1_R/DCOT_R - T2_R )

      SHADOW = ONE/(ONE+SHADOWI+SHADOWR)

      COXMUNK_VKERNEL(1) = COXMUNK_VKERNEL(1) * SHADOW

!  Scalar derivatives (1,1) matrix entry
!  -------------------------------------

!  inverse slope-squared derivative

      IF ( DO_DERIV_PARS(1) ) THEN

!  add the shadow

        D_S1 = HALF / PIE / S1
        D_S2 = - S2 * S2
        D_T1 = - T1_I * DCOT_I * DCOT_I * D_S2
        D_T2 = S2 * S2 * DCOT_I * S1 * T1_I
        D_SHADOWI = HALF * ( D_S1*T1_I/DCOT_I + S1*D_T1/DCOT_I - D_T2 )

        D_T1 = - T1_R * DCOT_R * DCOT_R * D_S2
        D_T2 = S2 * S2 * DCOT_R * S1 * T1_R
        D_SHADOWR = HALF * ( D_S1*T1_R/DCOT_R + S1*D_T1/DCOT_R - D_T2 )
        D_SHADOW = - SHADOW * SHADOW * ( D_SHADOWI + D_SHADOWR )

        COXMUNK_VDERIVATIVES(1,1) = COXMUNK_VDERIVATIVES(1,1) * SHADOW + &
                                    COXMUNK_VKERNEL(1) * D_SHADOW / &
                                    SHADOW

      ENDIF

!  square refractive index derivative

      IF ( DO_DERIV_PARS(2) ) THEN
        IF ( ARGUMENT .LT. CRITEXP ) THEN
          DH1 = Z2
          DH2 = HALF / H2
          DRP = ( DH1 * ( ONE - RP ) - DH2 * ( ONE + RP ) ) / H1H2
          DRL =  - DH2 * ( ONE + RL ) / H2Z2
          DFAC2 = ( RP*DRP + RL*DRL ) / XMP
          COXMUNK_VDERIVATIVES(2,1) = COXMUNK_VKERNEL(1) * DFAC2
        ENDIF
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE COXMUNK_VFUNCTION_PLUS

!^C

      SUBROUTINE GISSCOXMUNK_VFUNCTION_PLUS &
         ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, N_BRDF_STOKESSQ, &
           XJ, SXJ, XI, SXI, XPHI_REF, CKPHI_REF, SKPHI_REF, &
           GISSCOXMUNK_VKERNEL, GISSCOXMUNK_VDERIVATIVES )

!  include file of constants

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER ::          MAXPARS, NPARS, N_BRDF_STOKESSQ
      DOUBLE PRECISION :: PARS ( MAXPARS )
      LOGICAL ::          DO_DERIV_PARS ( MAXPARS )
      DOUBLE PRECISION :: XI, XJ
      DOUBLE PRECISION :: SXI, SXJ, XPHI_REF, CKPHI_REF, SKPHI_REF

      DOUBLE PRECISION :: GISSCOXMUNK_VKERNEL      ( MAXSTOKES_SQ )
      DOUBLE PRECISION :: GISSCOXMUNK_VDERIVATIVES ( MAXPARS, MAXSTOKES_SQ )

!  Critical exponent taken out

      DOUBLE PRECISION, PARAMETER :: CRITEXP = 88.0D0

!  Local variables

      INTEGER ::          KERNELMASK(10), I, IM
      DOUBLE PRECISION :: XPHI_INC, CKPHI_INC, SKPHI_INC
      DOUBLE PRECISION :: VI1, VI2, VI3, VR1, VR2, VR3
      DOUBLE PRECISION :: unit1, unit2, unit3, fact1, factor
      DOUBLE PRECISION :: XI1, CN1, CN2, CXI2, C2, C1, CRPER, CRPAR
      DOUBLE PRECISION :: TI1, TI2, TI3, TR1, TR2, TR3
      DOUBLE PRECISION :: PI1, PII2, PI3, PR1, PR2, PR3
      DOUBLE PRECISION :: PIKR, PRKI, TIKR, TRKI
      DOUBLE PRECISION :: E1, E2, E3, E4, SIGMA2
      DOUBLE PRECISION :: CF11, CF12, CF21, CF22, RDZ2, RDZ4
      DOUBLE PRECISION :: VP1, VP2, VP3, DMOD, DEX, DCOEFF
      DOUBLE PRECISION :: AF, AF11, AF12, AF21, AF22
      DOUBLE PRECISION :: C21, C22, CTTTP, CTTPT, CTTPP
      DOUBLE PRECISION :: CTPPT, CTPPP, CPTPP
      DOUBLE PRECISION :: S1, S2, S3, XXI, XXJ, T1, T2, DCOT
      DOUBLE PRECISION :: SHADOWI, SHADOWR, SHADOW, DCOEFF_0, ARGUMENT

      DOUBLE PRECISION :: D_S1, D_S2, D_T1, D_T2, DERFAC
      DOUBLE PRECISION :: D_SHADOWI, D_SHADOWR, D_SHADOW
      DOUBLE PRECISION :: D_DCOEFF_0, D_DCOEFF, D_DEX, D_ARGUMENT

      KERNELMASK = (/ 1,2,3,5,6,7,9,10,11,16 /)

!  Initialise

      GISSCOXMUNK_VKERNEL      = ZERO
      GISSCOXMUNK_VDERIVATIVES = ZERO

!  Transcription of the RMATR subroutine from Mishchenko/Travis code.

!   CALCULATION OF THE STOKES REFLECTION MATRIX FOR
!   ILLUMINATION FROM ABOVE FOR
!   A STATISTICALLY ROUGH SURFACE SEPARATING TWO HALF-SPACES
!   WITH REFRACTIVE INDICES OF THE UPPER AND LOWER HALF-SPACES EQUAL TO
!   CN1 AND CN2, RESPECTIVELY. THE EFFECT OF SHADOWING IS NOT
!   INCLUDED IN THIS       SUBROUTINE BUT IS ADDED IN THE MAIN PROGRAM.

!   SIGMA2 = s**2 = MEAN SQUARE SURFACE SLOPE (EQ. (18) IN THE JGR PAPER

!   XI = ABS(COSINE OF THE INCIDENT ZENITH ANGLE)
!   XJ = ABS(COSINE OF THE REFLECTION ZENITH ANGLE).
!   SXI and SXJ are the respective SINES (input)
!   XPHI_REF = REFLECTION AZIMUTH ANGLE
!   GISSCOXMUNK_VKERNEL(16-elements) = (4X4) REFLECTION MATRIX

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!               Remark on Use of shadow effect
!               ------------------------------
!  Shadow effect is controlled by the third parameter. That is, if
!  PARS(3) not equal to then shadow effect will be included.
!    --- NPARS should always be 3 for this Kernel.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  For real case, incident azimuth taken to be zero

      XPHI_INC  = ZERO
      CKPHI_INC = ONE
      SKPHI_INC = ZERO

!  Slope square is PARS(1)

      SIGMA2 = PARS(1)

!  Check for limiting cases

      IF(DABS(XI-1D0).LT.1d-9) XI = 0.999999999999d0
      IF(DABS(XJ-1D0).LT.1d-9) XJ = 0.999999999999d0

!  help variables (coordinate transformations)

      VI1 = SXI * CKPHI_INC
      VI2 = SXI * SKPHI_INC
      VI3 = -XI
      VR1 = SXJ * CKPHI_REF
      VR2 = SXJ * SKPHI_REF
      VR3 = XJ

!    LOCAL SURFACE NORMAL FOR SPECULAR REFLECTION (normalized to 1)

      UNIT1  = VI1-VR1
      UNIT2  = VI2-VR2
      UNIT3  = VI3-VR3
      FACT1  = UNIT1*UNIT1 + UNIT2*UNIT2 + UNIT3*UNIT3
      FACTOR = DSQRT(ONE/FACT1)

!   FRESNEL REFLECTION COEFFICIENTS, assume only real for now
!   ---------------------------------------------------------

      CN1 = ONE
      CN2 = PARS(2)

!  this is the original code, but now C-variables are real

      XI1 =  FACTOR*(UNIT1*VI1+UNIT2*VI2+UNIT3*VI3)
      CXI2 = ONE - (ONE-XI1*XI1)*CN1*CN1/(CN2*CN2)
      CXI2 = DSQRT(CXI2)
      C1 = CN1*XI1
      C2 = CN2*CXI2
      CRPER = (C1-C2)/(C1+C2)
      C1 = CN2*XI1
      C2 = CN1*CXI2
      CRPAR = (C1-C2)/(C1+C2)

!  CALCULATION OF THE AMPLITUDE SCATTERING MATRIX
!  ----------------------------------------------

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

!  CALCULATION OF THE STOKES REFLECTION MATRIX
!  -------------------------------------------

      VP1 = VI2*VR3-VI3*VR2
      VP2 = VI3*VR1-VI1*VR3
      VP3 = VI1*VR2-VI2*VR1
      DMOD = VP1*VP1+VP2*VP2+VP3*VP3
      DMOD = DMOD*DMOD

!  if DMOD = 0, that is | n x n_0 | ^ 4 = 0 in M-T formula)
!    Then we need to set the ratio CF11 / DMOD

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

!  Derivative of DCOEFF w.r.t. SIGMA2

      D_DCOEFF = ZERO
      IF ( DO_DERIV_PARS(1) ) THEN
       IF ( ARGUMENT .LT. CRITEXP ) THEN
        D_ARGUMENT = - ARGUMENT / SIGMA2
        D_DCOEFF_0 = - DCOEFF_0 / SIGMA2
        D_DEX      = - D_ARGUMENT * DEX
        D_DCOEFF = FACT1*FACT1 * ( DCOEFF_0 * D_DEX + D_DCOEFF_0 * DEX )
       ENDIF
      ENDIF

      AF  = HALF * DCOEFF
      AF11 = DABS(CF11)
      AF12 = DABS(CF12)
      AF21 = DABS(CF21)
      AF22 = DABS(CF22)
      AF11 = AF11*AF11
      AF12 = AF12*AF12
      AF21 = AF21*AF21
      AF22 = AF22*AF22

!  original code
!      R(1,1)=(AF11+AF12+AF21+AF22)*AF
!      R(1,2)=(AF11-AF12+AF21-AF22)*AF
!      R(2,1)=(AF11-AF22+AF12-AF21)*AF
!      R(2,2)=(AF11-AF12-AF21+AF22)*AF

!  Transcribed code
      GISSCOXMUNK_VKERNEL(1) = (AF11+AF12+AF21+AF22)*AF
      GISSCOXMUNK_VKERNEL(2) = (AF11-AF12+AF21-AF22)*AF
      GISSCOXMUNK_VKERNEL(5) = (AF11-AF22+AF12-AF21)*AF
      GISSCOXMUNK_VKERNEL(6) = (AF11-AF12-AF21+AF22)*AF

!  Original code
!      CI=(0D0, -1D0)
!      C21=DCONJG(CF21)
!      C22=DCONJG(CF22)
!      CTTTP=CF11*DCONJG(CF12)
!      CTTPT=CF11*C21
!      CTTPP=CF11*C22
!      CTPPT=CF12*C21
!      CTPPP=CF12*C22
!      CPTPP=CF21*C22

!  replica for real variables only
      C21 = CF21
      C22 = CF22
      CTTTP=CF11*CF12
      CTTPT=CF11*C21
      CTTPP=CF11*C22
      CTPPT=CF12*C21
      CTPPP=CF12*C22
      CPTPP=CF21*C22

!  original code
!      R(1,3)=    (-CTTTP-CPTPP)*DCOEFF
!      R(1,4)=-CI*( CTTTP+CPTPP)*DCOEFF
!      R(2,3)=    (-CTTTP+CPTPP)*DCOEFF
!      R(2,4)=-CI*( CTTTP-CPTPP)*DCOEFF
!      R(3,1)=    (-CTTPT-CTPPP)*DCOEFF
!      R(3,2)=    (-CTTPT+CTPPP)*DCOEFF
!      R(3,3)=    ( CTTPP+CTPPT)*DCOEFF
!      R(3,4)= CI*( CTTPP-CTPPT)*DCOEFF
!      R(4,1)= CI*( CTTPT+CTPPP)*DCOEFF
!      R(4,2)= CI*( CTTPT-CTPPP)*DCOEFF
!      R(4,3)=-CI*( CTTPP+CTPPT)*DCOEFF
!      R(4,4)=    ( CTTPP-CTPPT)*DCOEFF

!  New code (several entries are zero)

      GISSCOXMUNK_VKERNEL(3)  =    (-CTTTP-CPTPP)*DCOEFF
      GISSCOXMUNK_VKERNEL(7)  =    (-CTTTP+CPTPP)*DCOEFF
      GISSCOXMUNK_VKERNEL(9)  =    (-CTTPT-CTPPP)*DCOEFF
      GISSCOXMUNK_VKERNEL(10) =    (-CTTPT+CTPPP)*DCOEFF
      GISSCOXMUNK_VKERNEL(11) =    ( CTTPP+CTPPT)*DCOEFF
      GISSCOXMUNK_VKERNEL(16) =    ( CTTPP-CTPPT)*DCOEFF

!  Derivative before shadow effect

      DERFAC = ZERO
      IF ( DO_DERIV_PARS(1) ) THEN
       IF ( DCOEFF .NE. ZERO ) THEN
        DERFAC = D_DCOEFF / DCOEFF
       ENDIF
       DO I = 1, 10
        IM = KERNELMASK(I)
        GISSCOXMUNK_VDERIVATIVES(1,IM) = GISSCOXMUNK_VKERNEL(IM)*DERFAC
       ENDDO
      ENDIF

!  No Shadow code if not flagged

      IF ( PARS(3) .EQ. ZERO ) RETURN

!  Shadow code (includes derivative if flagged)

      S1 = DSQRT(TWO*SIGMA2/PIE)
      S3 = ONE/(DSQRT(TWO*SIGMA2))
      S2 = S3*S3
      D_S1 = ONE / PIE / S1
      D_S2 = - S2 / SIGMA2

      SHADOWI   = ZERO
      D_SHADOWI = ZERO
      IF ( XI .NE. ONE ) THEN
       XXI  = XI*XI
       DCOT = XI/DSQRT(ONE-XXI)
       T1   = DEXP(-DCOT*DCOT*S2)
       T2   = DERFC_E(DCOT*S3)
       SHADOWI = HALF*(S1*T1/DCOT-T2)
       IF ( DO_DERIV_PARS(1) ) THEN
        D_T1 = - T1 * DCOT * DCOT * D_S2
        D_T2 = TWO * S2 * DCOT * T1 / PIE / S1
        D_SHADOWI = HALF * ( D_S1*T1/DCOT + S1*D_T1/DCOT - D_T2 )
       ENDIF
      ENDIF

      SHADOWR   = ZERO
      D_SHADOWR = ZERO
      IF ( XJ .NE. ONE ) THEN
       XXJ  = XJ*XJ
       DCOT = XJ/DSQRT(ONE-XXJ)
       T1   = DEXP(-DCOT*DCOT*S2)
       T2   = DERFC_E(DCOT*S3)
       SHADOWR = HALF*(S1*T1/DCOT-T2)
       IF ( DO_DERIV_PARS(1) ) THEN
        D_T1 = - T1 * DCOT * DCOT * D_S2
        D_T2 = TWO * S2 * DCOT * T1 / PIE / S1
        D_SHADOWR = HALF * ( D_S1*T1/DCOT + S1*D_T1/DCOT - D_T2 )
       ENDIF
      ENDIF
      SHADOW = ONE/(ONE+SHADOWI+SHADOWR)

      DO I = 1, 10
       IM = KERNELMASK(I)
       GISSCOXMUNK_VKERNEL(IM) = GISSCOXMUNK_VKERNEL(IM) * SHADOW
      ENDDO

      IF ( DO_DERIV_PARS(1) ) THEN
       D_SHADOW = - SHADOW * SHADOW * ( D_SHADOWI + D_SHADOWR )
       DO I = 1, 10
        IM = KERNELMASK(I)
        GISSCOXMUNK_VDERIVATIVES(1,IM) = GISSCOXMUNK_VDERIVATIVES(1,IM) &
                                        *SHADOW +GISSCOXMUNK_VKERNEL(IM) &
                                        *D_SHADOW/SHADOW
       ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE GISSCOXMUNK_VFUNCTION_PLUS

!^C

      SUBROUTINE COXMUNK_VFUNCTION_DB_PLUS &
     ( MAXPARS, NPARS, PARS, DERIVS, ORDER, NSSQ,                        &
       n_muquad, n_phiquad, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,          & 
       X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,  &
       COXMUNK_VKERNEL, COXMUNK_VDERIVATIVES )

!  include file of constants

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER          :: MAXPARS, NPARS, NSSQ, ORDER
      DOUBLE PRECISION :: PARS   ( MAXPARS )
      LOGICAL          :: DERIVS ( MAXPARS )
      DOUBLE PRECISION :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION :: COXMUNK_VKERNEL      (MAXSTOKES_SQ)
      DOUBLE PRECISION :: COXMUNK_VDERIVATIVES ( MAXPARS, MAXSTOKES_SQ )

!  Local arrays for MSR quadrature

      INTEGER          :: n_muquad, n_phiquad
      DOUBLE PRECISION :: X_MUQUAD (max_msrs_muquad)
      DOUBLE PRECISION :: W_MUQUAD (max_msrs_muquad)
      DOUBLE PRECISION :: SX_MUQUAD (max_msrs_muquad)
      DOUBLE PRECISION :: WXX_MUQUAD (max_msrs_muquad)

      DOUBLE PRECISION :: X_PHIQUAD (max_msrs_phiquad)
      DOUBLE PRECISION :: W_PHIQUAD (max_msrs_phiquad)

!  local variables

      INTEGER         ::  s, n, k, ni, ki, nr, kr, q
      DOUBLE PRECISION :: XM, SXM, XMR, SXMR, XMI, SXMI, sum_pr, sumr, sum, w_p
      DOUBLE PRECISION :: reflec_0(16), reflec_s, R0Q, D_R0Q, reflec
      DOUBLE PRECISION :: phi_sub1, cphi_sub1, sphi_sub1
      DOUBLE PRECISION :: phi_sub2, cphi_sub2, sphi_sub2
      DOUBLE PRECISION :: d_reflec_0(3,16), d_reflec_s(3), d_reflec(3)

!  arrays

      DOUBLE PRECISION :: R0_QUAD_IN    (16,max_msrs_muquad,max_msrs_phiquad)
      DOUBLE PRECISION :: R0_OUT_QUAD   (16,max_msrs_muquad,max_msrs_phiquad)
      DOUBLE PRECISION :: D_R0_QUAD_IN  (MAXPARS,16,max_msrs_muquad,max_msrs_phiquad)
      DOUBLE PRECISION :: D_R0_OUT_QUAD (MAXPARS,16,max_msrs_muquad,max_msrs_phiquad)

      DOUBLE PRECISION :: RHOLD    (max_msrs_muquad,max_msrs_phiquad)
      DOUBLE PRECISION :: D_RHOLD  (3,max_msrs_muquad,max_msrs_phiquad)

      DOUBLE PRECISION :: R0_MSRS_QUAD   (16,max_msrs_muquad,max_msrs_phiquad)
      DOUBLE PRECISION :: D_R0_MSRS_QUAD (MAXPARS,16,max_msrs_muquad,max_msrs_phiquad)

!  Safety first zeroing

      REFLEC_0   = ZERO
      D_REFLEC_0 = ZERO
      COXMUNK_VKERNEL      = ZERO
      COXMUNK_VDERIVATIVES = ZERO

!  Only want the first element (this is a scalar routine)

!  Single scattering (zero order), Phi is in degrees here!

      CALL COXMUNK_VFUNCTION_PLUS &
         ( MAXPARS, NPARS, PARS, DERIVS, NSSQ, XJ, SXJ, XI, SXI, &
           PHI, CPHI, SKPHI, REFLEC_0, D_REFLEC_0 )

!  Higher orders scattering
!    Only want the first element (this is a scalar routine)

      REFLEC   = REFLEC_0(1)
      REFLEC_S = ZERO

      D_REFLEC(:)   = D_REFLEC_0(:,1)
      D_REFLEC_S(:) = ZERO

!  Quadrature output for first order R/T calculations
!        This will be overwritten as the orders increase

      IF ( ORDER.GE.1 ) THEN
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
               CALL COXMUNK_VFUNCTION_PLUS &
                ( MAXPARS, NPARS, PARS, DERIVS, NSSQ, XI, SXI, XM, SXM, PHI_SUB2, &
                  CPHI_SUB2, SPHI_SUB2, R0_OUT_QUAD(:,K,N), D_R0_OUT_QUAD(:,:,K,N) )
               CALL COXMUNK_VFUNCTION_PLUS &
                ( MAXPARS, NPARS, PARS, DERIVS, NSSQ, XM, SXM, XJ, SXJ, PHI_SUB1, &
                  CPHI_SUB1, SPHI_SUB1, R0_QUAD_IN(:,K,N), D_R0_QUAD_IN(:,:,K,N) )
            ENDDO
         ENDDO
      ENDIF

!  Compute the successive orders of scattering

!  compute the next order, (1,1) component only

      DO S = 1, ORDER

!  Compute result for this order

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
         REFLEC_S = SUMR

!  Derivatives

         DO Q = 1, 2
            IF ( DERIVS(Q) ) THEN
               SUMR = ZERO
               DO K = 1, n_muquad
                  SUM_PR = ZERO
                  DO N = 1, N_PHIQUAD
                     W_P  = W_PHIQUAD(N)
                     SUM_PR = SUM_PR + W_P *   R0_QUAD_IN(1,K,N)   * D_R0_OUT_QUAD(Q,1,K,N) &
                                     + W_P * D_R0_QUAD_IN(Q,1,K,N) *   R0_OUT_QUAD(1,K,N)
                  ENDDO
                  SUMR = SUMR + SUM_PR * WXX_MUQUAD(K)
               ENDDO
               D_REFLEC_S(Q) = SUMR
            ENDIF
         ENDDO

!  Finish if reached the scattering order desired

         IF ( S.EQ. ORDER ) GO TO 67

!  Compute Reflectance for next order and update
!    Quad-Quad results get computed each time: very wasteful.
!    Have to do this, as the memory is a killer.

         DO KR = 1, n_muquad
            XMR  = X_MUQUAD(KR)
            SXMR = SX_MUQUAD(KR)
            DO NR = 1, N_PHIQUAD

!  Quad-quad calculations

               DO KI = 1, n_muquad
                  XMI  = X_MUQUAD(KI)
                  SXMI = SX_MUQUAD(KI)
                  DO NI = 1, N_PHIQUAD
                     PHI_SUB1  = X_PHIQUAD(NR) - X_PHIQUAD(NI)
                     CPHI_SUB1 = DCOS(PHI_SUB1)
                     SPHI_SUB1 = DSIN(PHI_SUB1)
                     CALL COXMUNK_VFUNCTION_PLUS &
                       ( MAXPARS, NPARS, PARS, DERIVS, NSSQ, &
                         XMI, SXMI, XMR, SXMR, PHI_SUB1, CPHI_SUB1, SPHI_SUB1,&
                         R0_MSRS_QUAD(:,KI,NI), D_R0_MSRS_QUAD(:,:,KI,NI))
                  ENDDO
               ENDDO

!  Multiple reflection

               SUMR = ZERO
               DO KI = 1, N_MUQUAD
                  SUM_PR = ZERO
                  DO NI = 1, N_PHIQUAD
                     W_P  = W_PHIQUAD(NI)
                     R0Q = R0_MSRS_QUAD(1,KI,NI)
                     SUM_PR = SUM_PR + W_P * R0_QUAD_IN(1,KI,NI) * R0Q
                  ENDDO
!                  SUMR = SUMR + SUM_PR * W_MUQUAD(KI)
                  SUMR = SUMR + SUM_PR * WXX_MUQUAD(KI)
               ENDDO
               RHOLD(KR,NR) = SUMR

!  Derivatives

               DO Q = 1, 2
                  IF ( DERIVS(Q) ) THEN
                     SUMR = ZERO
                     DO KI = 1, N_MUQUAD
                        SUM_PR = ZERO
                        DO NI = 1, N_PHIQUAD
                           W_P  = W_PHIQUAD(NI)
                           R0Q   = R0_MSRS_QUAD(1,KI,NI)
                           D_R0Q = D_R0_MSRS_QUAD(Q,1,KI,NI)
                           SUM_PR = SUM_PR + W_P *   R0_QUAD_IN(1,KI,NI)   * D_R0Q &
                                           + W_P * D_R0_QUAD_IN(Q,1,KI,NI) *   R0Q
                        ENDDO
                        SUMR = SUMR +SUM_PR * WXX_MUQUAD(K)
                     ENDDO
                     D_RHOLD(Q,KR,NR) = SUMR
                  ENDIF
               ENDDO

!  End KR, NR loops

            ENDDO
         ENDDO     

!  Update

         DO KR = 1, N_MUQUAD
            DO NR = 1, N_PHIQUAD
               R0_QUAD_IN(1,KR,NR) = RHOLD(KR,NR)
            ENDDO
         ENDDO
         DO Q = 1, 2
            IF ( DERIVS(Q) ) THEN
               DO KR = 1, N_MUQUAD
                  DO NR = 1, N_PHIQUAD
                     D_R0_QUAD_IN(Q,1,KR,NR) = D_RHOLD(Q,KR,NR)
                  ENDDO
               ENDDO
            ENDIF
         ENDDO

!  Continuation point for finishing MSR
  
 67      continue

!  Add to total

          REFLEC   = REFLEC + REFLEC_S
          D_REFLEC = d_reflec + D_REFLEC_S

!  End scattering order loop

      ENDDO

!  Compute total

      COXMUNK_VKERNEL(1) = REFLEC
      DO Q = 1, 2
         IF ( DERIVS(Q) ) THEN
            COXMUNK_VDERIVATIVES(Q,1) = D_REFLEC(Q)
         ENDIF
      ENDDO

!  debug

!      DO M = 1, NSSQ
!        write(34,'(I5,1p6e14.5)')m,
!     &   reflec_0(m), reflec_1(m),COXMUNK_VKERNEL(m),
!     &   dacos(xi)/deg_to_rad, dacos(xj)/deg_to_rad, phi
!      enddo

!  Finish

      RETURN
      END SUBROUTINE COXMUNK_VFUNCTION_DB_PLUS

!

      SUBROUTINE GISSCOXMUNK_VFUNCTION_DB_PLUS &
     ( MAXPARS, NPARS, PARS, DERIVS, ORDER, NSSQ,                        &
       n_muquad, n_phiquad, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,          & 
       X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,  &
       GISSCOXMUNK_VKERNEL, GISSCOXMUNK_VDERIVATIVES )

!  include file of constants

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER ::          MAXPARS, NPARS, NSSQ, ORDER
      DOUBLE PRECISION :: PARS ( MAXPARS )
      LOGICAL          :: DERIVS ( MAXPARS )
      DOUBLE PRECISION :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION :: GISSCOXMUNK_VKERNEL(MAXSTOKES_SQ)
      DOUBLE PRECISION :: GISSCOXMUNK_VDERIVATIVES ( MAXPARS, MAXSTOKES_SQ )

!  Local arrays for MSR quadrature

      INTEGER          :: n_muquad, n_phiquad
      DOUBLE PRECISION :: X_MUQUAD (max_msrs_muquad)
      DOUBLE PRECISION :: W_MUQUAD (max_msrs_muquad)
      DOUBLE PRECISION :: SX_MUQUAD (max_msrs_muquad)
      DOUBLE PRECISION :: WXX_MUQUAD (max_msrs_muquad)

      DOUBLE PRECISION :: X_PHIQUAD (max_msrs_phiquad)
      DOUBLE PRECISION :: W_PHIQUAD (max_msrs_phiquad)

!  local variables

      INTEGER         ::  s, n, k, o1, o2, o3, ni, ki, nr, kr, q
      DOUBLE PRECISION :: XM, SXM, XMR, SXMR, XMI, SXMI, sum_pr(16), sumr(16), sum, sum1, sum2, w_p
      DOUBLE PRECISION :: reflec_0(16), reflec_s(16), R0Q, D_R0Q, reflec(16)
      DOUBLE PRECISION :: phi_sub1, cphi_sub1, sphi_sub1
      DOUBLE PRECISION :: phi_sub2, cphi_sub2, sphi_sub2
      DOUBLE PRECISION :: d_reflec_0(3,16), d_reflec_s(3,16), d_reflec(3,16)

!  arrays

      DOUBLE PRECISION :: R0_QUAD_IN    (16,max_msrs_muquad,max_msrs_phiquad)
      DOUBLE PRECISION :: R0_OUT_QUAD   (16,max_msrs_muquad,max_msrs_phiquad)
      DOUBLE PRECISION :: D_R0_QUAD_IN  (MAXPARS,16,max_msrs_muquad,max_msrs_phiquad)
      DOUBLE PRECISION :: D_R0_OUT_QUAD (MAXPARS,16,max_msrs_muquad,max_msrs_phiquad)

      DOUBLE PRECISION :: RHOLD    (16,max_msrs_muquad,max_msrs_phiquad)
      DOUBLE PRECISION :: D_RHOLD  (MAXPARS,16,max_msrs_muquad,max_msrs_phiquad)

      DOUBLE PRECISION :: R0_MSRS_QUAD   (16,max_msrs_muquad,max_msrs_phiquad)
      DOUBLE PRECISION :: D_R0_MSRS_QUAD (MAXPARS,16,max_msrs_muquad,max_msrs_phiquad)

!  Indices

      INTEGER ::          MASKIT(3,3), LNS, NS, NSS, M, M12, M13, M32

!  Safety first zeroing

      REFLEC_0   = ZERO
      D_REFLEC_0 = ZERO
      GISSCOXMUNK_VKERNEL      = ZERO
      GISSCOXMUNK_VDERIVATIVES = ZERO

!  Masking limits

      NSS = NSSQ
      IF ( NSSQ.EQ.1  ) LNS = 1
      IF ( NSSQ.EQ.4  ) LNS = 2
      IF ( NSSQ.EQ.9  ) LNS = 3
      NS = LNS
      IF ( NSSQ.EQ.16 ) THEN
        LNS = 3
        NS = LNS + 1
      ENDIF

!  masking array

      MASKIT(1,1) = 1
      MASKIT(1,2) = 2
      MASKIT(1,3) = 3
      MASKIT(2,1) = 5
      MASKIT(2,2) = 6
      MASKIT(2,3) = 7
      MASKIT(3,1) = 9
      MASKIT(3,2) = 10
      MASKIT(3,3) = 11

!  Single scattering (zero order), Phi is in degrees here!

      CALL GISSCOXMUNK_VFUNCTION_PLUS &
         ( MAXPARS, NPARS, PARS, DERIVS, NSSQ, XJ, SXJ, XI, SXI, &
           PHI, CPHI, SKPHI, REFLEC_0, D_REFLEC_0 )

!  Higher orders scattering
!    Only want the first element (this is a scalar routine)

      REFLEC   = REFLEC_0
      REFLEC_S = ZERO

      D_REFLEC   = D_REFLEC_0
      D_REFLEC_S = ZERO

!  Quadrature output for first order R/T calculations
!        This will be overwritten as the orders increase

      IF ( ORDER.GE.1 ) THEN
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
               CALL GISSCOXMUNK_VFUNCTION_PLUS &
                ( MAXPARS, NPARS, PARS, DERIVS, NSSQ, XI, SXI, XM, SXM, PHI_SUB2, &
                  CPHI_SUB2, SPHI_SUB2, R0_OUT_QUAD(:,K,N), D_R0_OUT_QUAD(:,:,K,N) )
               CALL GISSCOXMUNK_VFUNCTION_PLUS &
                ( MAXPARS, NPARS, PARS, DERIVS, NSSQ, XM, SXM, XJ, SXJ, PHI_SUB1, &
                  CPHI_SUB1, SPHI_SUB1, R0_QUAD_IN(:,K,N), D_R0_QUAD_IN(:,:,K,N) )
            ENDDO
         ENDDO
      ENDIF

!  Compute the successive orders of scattering

!  compute the next order, (1,1) component only

      DO S = 1, ORDER

!  Compute result for this order

          SUMR = ZERO
          DO K = 1, n_muquad
             SUM_PR = ZERO
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
!             SUMR(1:NSSQ) =  SUMR(1:NSSQ) + SUM_PR(1:NSSQ) * W_MUQUAD(K)  ! Warning
             SUMR(1:NSSQ) =  SUMR(1:NSSQ) + SUM_PR(1:NSSQ) * WXX_MUQUAD(K)  ! Warning
          ENDDO
          REFLEC_S(1:NSSQ) = SUMR(1:NSSQ)

!  Compute Derivative result for this order

          DO Q = 1, NPARS
            IF ( DERIVS(Q) ) THEN
              SUMR = ZERO
              DO K = 1, n_muquad
                SUM_PR = ZERO
                DO N = 1, N_PHIQUAD
                  W_P  = W_PHIQUAD(N)
                  DO O1 = 1, LNS
                    DO O2 = 1, LNS
                      M12 = MASKIT(O1,O2)
                      SUM1 = ZERO ; SUM2 = ZERO
                      DO O3 = 1, LNS
                        M13 = MASKIT(O1,O3)
                        M32 = MASKIT(O3,O2)
                        SUM1 = SUM1 +   R0_QUAD_IN(M13,K,N)   * D_R0_OUT_QUAD(Q,M32,K,N)
                        SUM2 = SUM2 + D_R0_QUAD_IN(Q,M13,K,N) *   R0_OUT_QUAD(M32,K,N)
                      ENDDO
                      SUM_PR(M12) = SUM_PR(M12) + W_P * ( SUM1 + SUM2 )
                    ENDDO
                  ENDDO
                  IF (NS.EQ.4) THEN
                    SUM =   R0_QUAD_IN(NSS,K,N)   * D_R0_OUT_QUAD(Q,NSS,K,N) + &
                          D_R0_QUAD_IN(Q,NSS,K,N) *   R0_OUT_QUAD(NSS,K,N)
                    SUM_PR(NSS) = SUM_PR(NSS) + W_P * SUM
                  ENDIF
                ENDDO
!                SUMR(1:NSSQ) =  SUMR(1:NSSQ) + SUM_PR(1:NSSQ) * W_MUQUAD(K)  ! Warning
                SUMR(1:NSSQ) =  SUMR(1:NSSQ) + SUM_PR(1:NSSQ) * WXX_MUQUAD(K)  ! Warning
              ENDDO
              D_REFLEC_S(Q,1:NSSQ) = SUMR(1:NSSQ)
            ENDIF
          ENDDO

!  Finish if reached the scattering order desired

          IF ( S.EQ. ORDER ) GO TO 67

!  Compute Reflectance for next order and update
!    Quad-Quad results get computed each time: very wasteful.
!    Have to do this, as the memory is a killer.

          DO KR = 1, n_muquad
            XMR  = X_MUQUAD(KR)
            SXMR = SX_MUQUAD(KR)
            DO NR = 1, N_PHIQUAD

!  Quad-quad calculations

               DO KI = 1, n_muquad
                  XMI  = X_MUQUAD(KI)
                  SXMI = SX_MUQUAD(KI)
                  DO NI = 1, N_PHIQUAD
                     PHI_SUB1  = X_PHIQUAD(NR) - X_PHIQUAD(NI)
                     CPHI_SUB1 = DCOS(PHI_SUB1)
                     SPHI_SUB1 = DSIN(PHI_SUB1)
                     CALL GISSCOXMUNK_VFUNCTION_PLUS &
                       ( MAXPARS, NPARS, PARS, DERIVS, NSSQ, &
                         XMI, SXMI, XMR, SXMR, PHI_SUB1, CPHI_SUB1, SPHI_SUB1,&
                         R0_MSRS_QUAD(:,KI,NI), D_R0_MSRS_QUAD(:,:,KI,NI))
                  ENDDO
               ENDDO

!  Multiple reflection

               SUMR = ZERO
               DO KI = 1, N_MUQUAD
                  SUM_PR = ZERO
                  DO NI = 1, N_PHIQUAD
                     W_P  = W_PHIQUAD(NI)
                     DO O1 = 1, LNS
                        DO O2 = 1, LNS
                           M12 = MASKIT(O1,O2)
                           SUM = ZERO
                           DO O3 = 1, LNS
                              M13 = MASKIT(O1,O3)
                              M32 = MASKIT(O3,O2)
                              R0Q = R0_MSRS_QUAD(M32,KI,NI)
                              SUM_PR(M12) = SUM_PR(M12) + W_P*R0_QUAD_IN(M13,KI,NI)*R0Q
                           ENDDO
                        ENDDO
                     ENDDO
                     IF (NS.EQ.4) THEN
                        SUM = R0_QUAD_IN(NSS,KI,NI) * R0_OUT_QUAD(NSS,KI,NI)
                        SUM_PR(NSS) = SUM_PR(NSS) + W_P * SUM
                     ENDIF
                  ENDDO
!                  SUMR(1:NSSQ) =  SUMR(1:NSSQ) + SUM_PR(1:NSSQ) * W_MUQUAD(KI) ! Warning
                  SUMR(1:NSSQ) =  SUMR(1:NSSQ) + SUM_PR(1:NSSQ) * WXX_MUQUAD(KI) ! Warning
               ENDDO
               RHOLD(1:NSSQ,KR,NR) = SUMR(1:NSSQ)

!  Derivatives Multiple reflection

               DO Q = 1, NPARS
                 IF ( DERIVS(Q) ) THEN
                   SUMR = ZERO
                   DO KI = 1, N_MUQUAD
                     SUM_PR = ZERO
                     DO NI = 1, N_PHIQUAD
                       W_P  = W_PHIQUAD(NI)
                       DO O1 = 1, LNS
                         DO O2 = 1, LNS
                           M12 = MASKIT(O1,O2)
                           SUM1 = ZERO ; SUM2 = ZERO
                           DO O3 = 1, LNS
                             M13 = MASKIT(O1,O3)
                             M32 = MASKIT(O3,O2)
                             D_R0Q = D_R0_MSRS_QUAD(Q,M32,KI,NI)
                             R0Q   = R0_MSRS_QUAD(M32,KI,NI)
                             SUM1 = D_R0_QUAD_IN(Q,M13,KI,NI) *   R0Q
                             SUM2 =   R0_QUAD_IN(M13,KI,NI)   * D_R0Q
                           ENDDO
                           SUM_PR(M12) = SUM_PR(M12) + W_P * ( SUM1 + SUM2 )
                         ENDDO
                       ENDDO
                       IF (NS.EQ.4) THEN
                         SUM =   R0_QUAD_IN(NSS,KI,NI)   * D_R0_OUT_QUAD(Q,NSS,KI,NI) + &
                               D_R0_QUAD_IN(Q,NSS,KI,NI) *   R0_OUT_QUAD(NSS,KI,NI)
                         SUM_PR(NSS) = SUM_PR(NSS) + W_P * SUM
                       ENDIF
                     ENDDO
!                     SUMR(1:NSSQ) =  SUMR(1:NSSQ) + SUM_PR(1:NSSQ) * W_MUQUAD(KI) ! Warning   
                     SUMR(1:NSSQ) =  SUMR(1:NSSQ) + SUM_PR(1:NSSQ) * WXX_MUQUAD(KI) ! Warning
                   ENDDO
                   D_RHOLD(Q,1:NSSQ,KR,NR) = SUMR(1:NSSQ)
                 ENDIF
               ENDDO

!  End loop over KR and NR

            ENDDO
         ENDDO     

!  Update

         DO KR = 1, N_MUQUAD
            DO NR = 1, N_PHIQUAD
               R0_QUAD_IN(1:NSSQ,KR,NR) = RHOLD(1:NSSQ,KR,NR)
            ENDDO
         ENDDO
         DO Q = 1, NPARS
            IF ( DERIVS(Q) ) THEN
               DO KR = 1, N_MUQUAD
                  DO NR = 1, N_PHIQUAD
                     D_R0_QUAD_IN(Q,1:NSSQ,KR,NR) = D_RHOLD(Q,1:NSSQ,KR,NR)
                  ENDDO
               ENDDO
            ENDIF
         ENDDO

!  Continuation point for finishing MSR
  
 67      continue

!  Add to total

         REFLEC(1:NSSQ) = REFLEC(1:NSSQ) + REFLEC_S(1:NSSQ)
         DO Q = 1, NPARS
           D_REFLEC(Q,1:NSSQ) = d_reflec(Q,1:NSSQ) + D_REFLEC_S(Q,1:NSSQ)
         ENDDO

!  End scattering order loop

      ENDDO

!  Compute total

      GISSCOXMUNK_VKERNEL(1:NSSQ) = REFLEC(1:NSSQ)
      DO Q = 1, NPARS
         IF ( DERIVS(Q) ) THEN
            GISSCOXMUNK_VDERIVATIVES(Q,1:NSSQ) = D_REFLEC(Q,1:NSSQ)
         ENDIF
      ENDDO

!  debug

!      DO M = 1, NSSQ
!        write(34,'(I5,1p6e14.5)')m,                      &
!        reflec_0(m), reflec_1(m),GISSCOXMUNK_VKERNEL(m), &
!        dacos(xi)/deg_to_rad, dacos(xj)/deg_to_rad, phi
!      enddo

!  Finish

      RETURN
      END SUBROUTINE GISSCOXMUNK_VFUNCTION_DB_PLUS

!  End module

      END MODULE vbrdf_LinSup_kernels_m
