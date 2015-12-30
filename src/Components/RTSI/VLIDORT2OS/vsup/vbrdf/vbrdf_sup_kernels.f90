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
! #  Release Date :   May 2014       (2.7)                      #
! #  Release Date :   August 2014    (2.7)                      #
! #                                                             #
! #       NEW: TOTAL COLUMN JACOBIANS          (2.4)            #
! #       NEW: BPDF Land-surface KERNELS       (2.4R)           #
! #       NEW: Thermal Emission Treatment      (2.4RT)          #
! #       Consolidated BRDF treatment          (2.4RTC)         #
! #       f77/f90 Release                      (2.5)            #
! #       External SS / New I/O Structures     (2.6)            #
! #                                                             #
! #       SURFACE-LEAVING / BRDF-SCALING      (2.7)             #
! #       TAYLOR Series / OMP THREADSAFE      (2.7)             #
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
! # ==========================                                  #
! #                                                             #
! # Listing for Version 2.3. (2008)                             #
! #                                                             #
! #            LAMBERTIAN_VFUNCTION                             #
! #            ROSSTHIN_VFUNCTION                               #
! #            ROSSTHICK_VFUNCTION                              #
! #            LISPARSE_VFUNCTION                               #
! #            LIDENSE_VFUNCTION                                #
! #            ROUJEAN_VFUNCTION                                #
! #            HAPKE_VFUNCTION                                  #
! #            RAHMAN_VFUNCTION                                 #
! #            COXMUNK_VFUNCTION                                #
! #            COXMUNK_VFUNCTION_MSR                            #
! #            GISSCOXMUNK_VFUNCTION                            #
! #            GISSCOXMUNK_VFUNCTION_MSR                        #
! #                                                             #
! # Additional code for complex RI Giss Cox-Munk, 15 march 2010.#
! #                                                             #
! #            GCMCRI_VFUNCTION                                 #
! #            GCMCRI_VFUNCTION_MSR                             #
! #                                                             #
! #  This kernel introduced for Version 2.4R                    #
! #  2009 function is final Kernel supplied by Breon, 5/5/09    #
! #            [ BPDF2009_VFUNCTION ] REPLACED                  #
! #                                                             #
! # New BPDF Subroutines in this Module (Version 2.7)           #
! #                                                             #
! #            BPDFVEGN_VFUNCTION                               #
! #            BPDFSOIL_VFUNCTION                               #
! #            BPDFNDVI_VFUNCTION                               #
! #            FRESNEL_VECTOR (Private, called by BPDF)         #
! #                                                             #
! # New Cox-Munk Subroutines in this Module (Version 2.7)       #
! #                                                             #
! #            VBRDF_Generalized_glint                          #
! #            VBRDF_Water_RefracIndex                          #
! #            VBRDF_WhiteCap_Reflectance                       #
! #                                                             #
! ###############################################################


      MODULE vbrdf_sup_kernels_m

!  Rob Extension 12/2/14. BPDF Kernels (replaces BPDF2009_VFUNCTION)

      use VLIDORT_PARS
      use vbrdf_sup_aux_m, only : derfc_e, VBRDF_Fresnel_Complex

      PRIVATE
      PUBLIC :: LAMBERTIAN_VFUNCTION,     &
                ROSSTHIN_VFUNCTION,       &
                ROSSTHICK_VFUNCTION,      &
                LISPARSE_VFUNCTION,       &
                LIDENSE_VFUNCTION,        &
                ROUJEAN_VFUNCTION,        &
                HAPKE_VFUNCTION,          &
                RAHMAN_VFUNCTION,         &
                COXMUNK_VFUNCTION,        &
                COXMUNK_VFUNCTION_MSR,    &
                GISSCOXMUNK_VFUNCTION,    &
                GISSCOXMUNK_VFUNCTION_MSR,&
                GCMCRI_VFUNCTION,         &
                GCMCRI_VFUNCTION_MSR,     &
                BPDFVEGN_VFUNCTION,       &
                BPDFSOIL_VFUNCTION,       &
                BPDFNDVI_VFUNCTION,       &
                VBRDF_Generalized_Glint,  &
                VBRDF_Water_RefracIndex,  &
                VBRDF_WhiteCap_Reflectance       

!  These have been replaced
!                BPDF2009_VFUNCTION,       &

      CONTAINS

      SUBROUTINE LAMBERTIAN_VFUNCTION &
         ( MAXPARS, NPARS, PARS, NSSQ, XJ, SXJ, XI, SXI, &
           PHI, CPHI, SKPHI, LAMBERTIAN_VKERNEL )

!  include file of constants

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER ::          MAXPARS, NPARS, NSSQ
      DOUBLE PRECISION :: PARS ( MAXPARS )
      DOUBLE PRECISION :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION :: LAMBERTIAN_VKERNEL(MAXSTOKES_SQ)

!  Local variables

      INTEGER ::          O1

!  Initialise

      DO O1 = 1, NSSQ
        LAMBERTIAN_VKERNEL(O1) = ZERO
      ENDDO

!  kernel: (1,1) Function (scalar form)

      LAMBERTIAN_VKERNEL(1) = ONE

!  Finish

      RETURN
      END SUBROUTINE LAMBERTIAN_VFUNCTION

!

      SUBROUTINE ROSSTHIN_VFUNCTION &
         ( MAXPARS, NPARS, PARS, NSSQ, XJ, SXJ, XI, SXI, &
           PHI, CPHI, SKPHI, ROSSTHIN_VKERNEL )

!  include file of constants

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER ::          MAXPARS, NPARS, NSSQ
      DOUBLE PRECISION :: PARS ( MAXPARS )
      DOUBLE PRECISION :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION :: ROSSTHIN_VKERNEL(MAXSTOKES_SQ)

!  Local variables

      INTEGER ::          O1
      DOUBLE PRECISION :: DS1, DS2, CKSI, SKSI, KSI, FUNC
      DOUBLE PRECISION :: XPHI, CKPHI

!  Initialise

      DO O1 = 1, NSSQ
        ROSSTHIN_VKERNEL(O1) = ZERO
      ENDDO
      XPHI = PIE - PHI
      CKPHI = - CPHI

!  kernel: (1,1) Function (scalar form)

      DS1 = XI * XJ
      DS2 = SXI * SXJ
      CKSI = DS1 + DS2 * CKPHI
      IF ( CKSI.GT.ONE ) CKSI = ONE
      SKSI = DSQRT(ONE-CKSI*CKSI)
      KSI = DACOS(CKSI)
      FUNC = ((PIO2-KSI)*CKSI + SKSI)/DS1
      ROSSTHIN_VKERNEL(1) = FUNC - PIO2

!  Finish

      RETURN
      END SUBROUTINE ROSSTHIN_VFUNCTION

!

      SUBROUTINE ROSSTHICK_VFUNCTION &
         ( MAXPARS, NPARS, PARS, NSSQ, XJ, SXJ, XI, SXI, &
           PHI, CPHI, SKPHI, ROSSTHICK_VKERNEL )

!  include file of constants

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER ::          MAXPARS, NPARS, NSSQ
      DOUBLE PRECISION :: PARS ( MAXPARS )
      DOUBLE PRECISION :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION :: ROSSTHICK_VKERNEL(MAXSTOKES_SQ)

!  Local variables

      INTEGER          :: O1
      DOUBLE PRECISION :: DS1, DS2, DS3,CKSI, SKSI, KSI, FUNC
      DOUBLE PRECISION :: XPHI, CKPHI

!  Initialise

      DO O1 = 1, NSSQ
        ROSSTHICK_VKERNEL(O1) = ZERO
      ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  (1,1) Function (scalar form)

      DS1 = XI * XJ
      DS2 = SXI * SXJ
      DS3 = XI  + XJ
      CKSI = DS1 + DS2 * CKPHI
      IF ( CKSI.GT.ONE ) CKSI = ONE
      SKSI = DSQRT(ONE-CKSI*CKSI)
      KSI = DACOS(CKSI)
      FUNC = ((PIO2-KSI)*CKSI + SKSI)/DS3
      ROSSTHICK_VKERNEL(1) = FUNC - PIO4

      RETURN
      END SUBROUTINE ROSSTHICK_VFUNCTION

!

      SUBROUTINE ROUJEAN_VFUNCTION &
         ( MAXPARS, NPARS, PARS, NSSQ, XJ, SXJ, XI, SXI, &
           PHI, CPHI, SKPHI, ROUJEAN_VKERNEL )

!  include file of constants

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER ::          MAXPARS, NPARS, NSSQ
      DOUBLE PRECISION :: PARS ( MAXPARS )
      DOUBLE PRECISION :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION :: ROUJEAN_VKERNEL(MAXSTOKES_SQ)

!  Local variables

      INTEGER          :: O1
      DOUBLE PRECISION :: DS1, DS2, DS3, TXJ, TXI, PHIFAC, S1, S2
      DOUBLE PRECISION :: XPHI_R, CXPHI_R, SXPHI_R, XPHI_C
      DOUBLE PRECISION :: XPHI, CKPHI

!  Initialise

      DO O1 = 1, NSSQ
        ROUJEAN_VKERNEL(O1) = ZERO
      ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  (1,1) Function (scalar form)

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

      RETURN
      END SUBROUTINE ROUJEAN_VFUNCTION

!

      SUBROUTINE LISPARSE_VFUNCTION &
         ( MAXPARS, NPARS, PARS, NSSQ, XJ, SXJ, XI, SXI, &
           PHI, CPHI, SKPHI, LISPARSE_VKERNEL )

!  include file of constants

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER ::          MAXPARS, NPARS, NSSQ
      DOUBLE PRECISION :: PARS ( MAXPARS )
      DOUBLE PRECISION :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION :: LISPARSE_VKERNEL(MAXSTOKES_SQ)

!  local variables

      INTEGER          :: O1
      DOUBLE PRECISION :: X_INC, X_REF, SX_INC, SX_REF, ANG_P, TX
      DOUBLE PRECISION :: T_INC, T_REF, T_INC_SQ, T_REF_SQ
      DOUBLE PRECISION :: CKSI, DELTA, T, COST, SINT, DSQ, SINTCOST
      DOUBLE PRECISION :: A, B, H, R, P, Q, DT1, DT2, DT2SQ, QR
      DOUBLE PRECISION :: XPHI, CKPHI

!  Initialise
!    -- Return for special case

      DO O1 = 1, NSSQ
        LISPARSE_VKERNEL(O1) = ZERO
      ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  (1,1) Function (scalar form)
!   9/27/14 No this is not right
!      IF ( ( XI .EQ. XJ ) .AND. ( CKPHI.EQ.ONE ) ) RETURN

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

!  contributions P and R

!    Bug found by Huan Ye (BIRA), Fixed by R. Spurr 12 October 2012
!      For the P term, Division by X_INC omitted in earlier version

! Old     P = ( ONE + CKSI ) / X_REF
      P = ( ONE + CKSI ) / X_REF / X_INC
      A = ( ONE / X_INC )
      B = ( ONE / X_REF )
      R = A + B

!  evaluate cos(t)

      DT1   = T_REF_SQ + T_INC_SQ
      DT2   = T_INC * T_REF
      DT2SQ = DT2 * DT2
      DELTA = DSQRT ( DT1 - TWO * DT2 * CKPHI )
      DSQ   = DELTA * DELTA
      H     = DSQRT ( DSQ + SKPHI * SKPHI * DT2SQ )
      COST  = PARS(1) * H / R

!  set Q function

      IF ( COST .GT. ONE ) THEN
        Q = ONE
      ELSE
        T        = DACOS(COST)
        SINT     = DSQRT ( ONE - COST * COST )
        SINTCOST = SINT * COST
        Q = ONE -  ( ( T - SINTCOST ) / PIE )
      ENDIF

!  set the kernel
!  --------------

      QR = Q * R
      LISPARSE_VKERNEL(1) = HALF * P - QR

      RETURN
      END SUBROUTINE LISPARSE_VFUNCTION

!

      SUBROUTINE LIDENSE_VFUNCTION &
         ( MAXPARS, NPARS, PARS, NSSQ, XJ, SXJ, XI, SXI, &
           PHI, CPHI, SKPHI, LIDENSE_VKERNEL )

!  include file of constants

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER ::          MAXPARS, NPARS, NSSQ
      DOUBLE PRECISION :: PARS ( MAXPARS )
      DOUBLE PRECISION :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION :: LIDENSE_VKERNEL(MAXSTOKES_SQ)

!  local variables

      INTEGER          :: O1
      DOUBLE PRECISION :: X_INC, X_REF, SX_INC, SX_REF, ANG_P, TX
      DOUBLE PRECISION :: T_INC, T_REF, T_INC_SQ, T_REF_SQ
      DOUBLE PRECISION :: CKSI, DELTA, T, COST, SINT, DSQ, SINTCOST
      DOUBLE PRECISION :: A, B, H, R, P, Q, DT1, DT2, DT2SQ, P_QR
      DOUBLE PRECISION :: XPHI, CKPHI

!  Initialise
!    -- Return for special case

      DO O1 = 1, NSSQ
        LIDENSE_VKERNEL(O1) = ZERO
      ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  (1,1) Function (scalar form)
!   No this is not right
!      IF ( ( XI .EQ. XJ ) .AND. ( CKPHI.EQ.ONE ) ) RETURN

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

!  contributions P and R

!    Bug found by Huan Ye (BIRA), Fixed by R. Spurr 12 October 2012
!      For the P term, Division by X_INC omitted in earlier version

! Old     P = ( ONE + CKSI ) / X_REF
      P = ( ONE + CKSI ) / X_REF / X_INC
      A = ( ONE / X_INC )
      B = ( ONE / X_REF )
      R = A + B

!  evaluate cos(t)

      DT1   = T_REF_SQ + T_INC_SQ
      DT2   = T_INC * T_REF
      DT2SQ = DT2 * DT2
      DELTA = DSQRT ( DT1 - TWO * DT2 * CKPHI )
      DSQ   = DELTA * DELTA
      H     = DSQRT ( DSQ + SKPHI * SKPHI * DT2SQ )
      COST  = PARS(1) * H / R

!  set Q function

      IF ( COST .GT. ONE ) THEN
        Q = ONE
      ELSE
        T        = DACOS(COST)
        SINT     = DSQRT ( ONE - COST * COST )
        SINTCOST = SINT * COST
        Q = ONE -  ( ( T - SINTCOST ) / PIE )
      ENDIF

!  set the kernel
!  --------------

      P_QR = P / Q / R
      LIDENSE_VKERNEL(1) = P_QR - TWO

      RETURN
      END SUBROUTINE LIDENSE_VFUNCTION

!

      SUBROUTINE HAPKE_VFUNCTION &
         ( MAXPARS, NPARS, PARS, NSSQ, XJ, SXJ, XI, SXI, &
           PHI, CPHI, SKPHI, HAPKE_VKERNEL )

!  include file of constants

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER ::          MAXPARS, NPARS, NSSQ
      DOUBLE PRECISION :: PARS ( MAXPARS )
      DOUBLE PRECISION :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION :: HAPKE_VKERNEL(MAXSTOKES_SQ)

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

!  local variables

      INTEGER ::          O1
      DOUBLE PRECISION :: CTHETA, THETA, PHASE
      DOUBLE PRECISION :: HOTSPOT, B0_EMPIR, HELP_HOT, B_HOT
      DOUBLE PRECISION :: SSALBEDO, GAMMA, REFLEC, FUNCT
      DOUBLE PRECISION :: HELP_J, TERM_J, HELP_I, TERM_I
      DOUBLE PRECISION :: XPHI, CKPHI

!  Initialise

      DO O1 = 1, NSSQ
        HAPKE_VKERNEL(O1) = ZERO
      ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  (1,1) Function (scalar form)

!  geometrical part

!  This is the code that is in DISORT - not right, I think.
!       CTHETA = XI * XJ + DABS(SXI) *  DABS(SXJ) * CKPHI

      CTHETA = XI * XJ + SXI * SXJ * CKPHI
      IF ( CTHETA .GT. ONE ) CTHETA = ONE
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
      TERM_J   = ( ONE + HELP_J ) / ( ONE + HELP_J * GAMMA )
      HELP_I   = TWO * XI
      TERM_I   = ( ONE + HELP_I ) / ( ONE + HELP_I * GAMMA )

!  Function

      REFLEC       = SSALBEDO * QUARTER / ( XI + XJ )
      FUNCT        = ( ONE + B_HOT ) * PHASE + TERM_J * TERM_I - ONE
      HAPKE_VKERNEL(1) = REFLEC * FUNCT

      RETURN
      END SUBROUTINE HAPKE_VFUNCTION

!

      SUBROUTINE RAHMAN_VFUNCTION &
         ( MAXPARS, NPARS, PARS, NSSQ, XJ, SXJ, XI, SXI, &
           PHI, CPHI, SKPHI, RAHMAN_VKERNEL )

!  include file of constants

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER ::          MAXPARS, NPARS, NSSQ
      DOUBLE PRECISION :: PARS ( MAXPARS )
      DOUBLE PRECISION :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION :: RAHMAN_VKERNEL(MAXSTOKES_SQ)

!  local variables

      INTEGER ::          O1
      DOUBLE PRECISION :: T_INC, T_REF, DT1, DT2
      DOUBLE PRECISION :: CXI, DELTA, K1_SQ, FACT
      DOUBLE PRECISION :: GEOM, PHASE, RFAC, K0, K1, K2
      DOUBLE PRECISION :: XPHI, CKPHI

!  Initialise

      DO O1 = 1, NSSQ
        RAHMAN_VKERNEL(O1) = ZERO
      ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  (1,1) Function (scalar form)

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
      FACT  = ( ONE + K1_SQ + TWO * K1 * CXI ) ** ONEP5
      PHASE = ( ONE - K1_SQ ) / FACT

!  Delta and R-factor

      T_INC = SXI / XI
      T_REF = SXJ / XJ
      DT1   = T_INC*T_INC + T_REF*T_REF
      DT2   = T_INC * T_REF
      DELTA = DSQRT ( DT1 - TWO * DT2 * CKPHI )
      RFAC = ( ONE - K0 ) / ( ONE + DELTA )

!  Geom factor and kernel

      GEOM = ( XI * XJ * ( XI + XJ ) ) ** ( K2 - ONE)
      RAHMAN_VKERNEL(1) = K0 * PHASE * ( ONE + RFAC ) * GEOM

!  Other functions

!      Placeholder

!  Finish

      RETURN
      END SUBROUTINE RAHMAN_VFUNCTION

!

      SUBROUTINE HAPKE_VFUNCTION_OLD &
         ( MAXPARS, NPARS, PARS, NSSQ, XJ, SXJ, XI, SXI, &
           PHI, CPHI, SKPHI, HAPKE_VKERNEL )

!  include file of constants

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER ::          MAXPARS, NPARS, NSSQ
      DOUBLE PRECISION :: PARS ( MAXPARS )
      DOUBLE PRECISION :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION :: HAPKE_VKERNEL(MAXSTOKES_SQ)

!  Hapke Kernel function.
!   From DISORT code; used as validation.

!  local variables

      INTEGER          :: O1
      DOUBLE PRECISION :: XPHI, CKPHI
      REAL             :: MU, MUP, DPHI

!  Initialise

      DO O1 = 1, NSSQ
        HAPKE_VKERNEL(O1) = ZERO
      ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  (1,1) Function (scalar form)

      MUP = SNGL(XJ)
      MU = SNGL(XI)
      DPHI = SNGL(XPHI)
      HAPKE_VKERNEL(1) = HAPKEKER( MU, MUP, DPHI )

      RETURN
      END SUBROUTINE HAPKE_VFUNCTION_OLD

!

      DOUBLE PRECISION FUNCTION hapkeker( MU, MUP, DPHI )

!      Supplies surface bi-directional reflectivity.

!      NOTE 1: Bidirectional reflectivity in DISORT is defined
!              by Eq. 39 in STWL.
!      NOTE 2: Both MU and MU0 (cosines of reflection and incidence
!              angles) are positive.

!  INPUT:

!    MU     : Cosine of angle of reflection (positive)

!    MUP    : Cosine of angle of incidence (positive)

!    DPHI   : Difference of azimuth angles of incidence and reflection
!                (radians)

!  LOCAL VARIABLES:

!    IREF   : bidirectional reflectance options
!             1 - Hapke's BDR model
!
!    B0     : empirical factor to account for the finite size of
!             particles in Hapke's BDR model

!    B      : term that accounts for the opposition effect
!             (retroreflectance, hot spot) in Hapke's BDR model

!    CTHETA : cosine of phase angle in Hapke's BDR model

!    GAMMA  : albedo factor in Hapke's BDR model

!    H0     : H( mu0 ) in Hapke's BDR model

!    H      : H( mu ) in Hapke's BDR model

!    HH     : angular width parameter of opposition effect in Hapke's
!             BDR model

!    P      : scattering phase function in Hapke's BDR model

!    THETA  : phase angle (radians); the angle between incidence and
!             reflection directions in Hapke's BDR model

!    W      : single scattering albedo in Hapke's BDR model

!   Called by- DREF, SURFAC
! +-------------------------------------------------------------------+
!     .. Scalar Arguments ..

      REAL ::    DPHI, MU, MUP
!     ..
!     .. Local Scalars ..

      INTEGER :: IREF
      REAL ::    B0, B, CTHETA, GAMMA, H0, H, HH, P, THETA, W
!     ..
!     .. Intrinsic Functions ..

      INTRINSIC COS, SQRT
!     ..

      IREF = 1

      IF ( IREF.EQ.1 ) THEN

!                              ** Hapke's BRDF model (times Pi/Mu0)
!                              ** (Hapke, B., Theory of reflectance
!                              ** and emittance spectroscopy, Cambridge
!                              ** University Press, 1993, Eq. 8.89 on
!                              ** page 233. Parameters are from
!                              ** Fig. 8.15 on page 231, expect for w.)

         CTHETA = MU * MUP + (1.-MU**2)**.5 * (1.-MUP**2)**.5 * &
                  COS( DPHI )
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
      END FUNCTION hapkeker

!

      SUBROUTINE COXMUNK_VFUNCTION &
         ( MAXPARS, NPARS, PARS, NSSQ, XJ, SXJ, XI, SXI, &
           PHI, CKPHI, SKPHI, COXMUNK_VKERNEL )

!  include file of constants

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER ::          MAXPARS, NPARS, NSSQ
      DOUBLE PRECISION :: PARS ( MAXPARS )
      DOUBLE PRECISION :: XI, SXI, XJ, SXJ, PHI, CKPHI, SKPHI
      DOUBLE PRECISION :: COXMUNK_VKERNEL(MAXSTOKES_SQ)

!  Critical exponent taken out

      DOUBLE PRECISION, PARAMETER :: CRITEXP = 88.0D0

!  local variables

      INTEGER          :: O1
      DOUBLE PRECISION :: Z, Z1, Z2, Z2_SQ_M1, H1, H2, RP, RL, XMP
      DOUBLE PRECISION :: A, B, TA, ARGUMENT, PROB, FAC1, FAC2, CKPHI_NEG
      DOUBLE PRECISION :: S1, S2, S3, XXI, XXJ, T1, T2, DCOT
      DOUBLE PRECISION :: SHADOWI, SHADOWR, SHADOW

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!               Remark on Use of shadow effect
!               ------------------------------
!  Shadow effect is controlled by the third parameter. That is, if
!  PARS(3) not equal to then shadow effect will be included.
!    --- NPARS should always be 3 for this Kernel.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Initialise

      DO O1 = 1, NSSQ
        COXMUNK_VKERNEL(O1) = ZERO
      ENDDO

!  Comment. 18 January 2006.
!  We have found in comparisons with the Giss Cox-Munk code that
!  the input COSPHI (CKPHI) here is the negative of what we actually nee
!  so introduce local variable which takes care of this

!  Also removed factor of PIE in the kernel denominator
!   This makes the output exactly same as GISS model for R(1,1)

      CKPHI_NEG = - CKPHI

!  (1,1) Function (scalar form)

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
      RP = ( H1 - H2 ) / ( H1 + H2 )
      RL = ( Z2 - H2 ) / ( Z2 + H2 )
      XMP = HALF * ( RP*RP + RL*RL )

!  Coxmunk Function

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

!  No Shadow code if not flagged

      IF ( PARS(3) .EQ. ZERO ) RETURN

!  Shadow code

      S1 = DSQRT(PARS(1)/PIE)
      S3 = ONE/(DSQRT(PARS(1)))
      S2 = S3*S3

      XXI  = XI*XI
      DCOT = XI/DSQRT(ONE-XXI)
      T1   = DEXP(-DCOT*DCOT*S2)
!      T2   = DERFC(DCOT*S3)
      T2   = DERFC_E(DCOT*S3)
      SHADOWI = HALF*(S1*T1/DCOT-T2)

      XXJ  = XJ*XJ
      DCOT = XJ/DSQRT(ONE-XXJ)
      T1   = DEXP(-DCOT*DCOT*S2)
!      T2   = DERFC(DCOT*S3)
      T2   = DERFC_E(DCOT*S3)
      SHADOWR = HALF*(S1*T1/DCOT-T2)

      SHADOW = ONE/(ONE+SHADOWI+SHADOWR)

      COXMUNK_VKERNEL(1) = COXMUNK_VKERNEL(1) * SHADOW

!  Finish

      RETURN
      END SUBROUTINE COXMUNK_VFUNCTION

!

      SUBROUTINE GISSCOXMUNK_VFUNCTION &
         ( MAXPARS, NPARS, PARS, NSSQ, XJ, SXJ, XI, SXI, &
           XPHI_REF, CKPHI_REF, SKPHI_REF, GISSCOXMUNK_VKERNEL )

!  include file of constants

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER ::          MAXPARS, NPARS, NSSQ
      DOUBLE PRECISION :: PARS ( MAXPARS )
      DOUBLE PRECISION :: XI, XJ
      DOUBLE PRECISION :: SXI, SXJ, XPHI_REF, CKPHI_REF, SKPHI_REF
      DOUBLE PRECISION :: GISSCOXMUNK_VKERNEL(MAXSTOKES_SQ)

!  Critical exponent taken out

      DOUBLE PRECISION, PARAMETER :: CRITEXP = 88.0D0

!  local variables

      INTEGER ::          O1, KERNELMASK(10), I, IM
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

      KERNELMASK = (/1,2,3,5,6,7,9,10,11,16/)

!  Initialise

      DO O1 = 1, NSSQ
        GISSCOXMUNK_VKERNEL(O1) = ZERO
      ENDDO

!  Transcription of the RMATR subroutine from Mishchenko/Travis code.

!   CALCULATION OF THE STOKES REFLECTION MATRIX FOR
!   ILLUMINATION FROM ABOVE FOR
!   A STATISTICALLY ROUGH SURFACE SEPARATING TWO HALF-SPACES
!   WITH REFRACTIVE INDICES OF THE UPPER AND LOWER HALF-SPACES EQUAL TO
!   CN1 AND CN2, RESPECTIVELY. THE EFFECT OF SHADOWING IS NOT
!   INCLUDED IN THIS SUBROUTINE BUT IS ADDED IN THE MAIN PROGRAM.

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

!  Key Debug statement to track down DMOD
!      write(*,*)'(1,1) Fresnel = ',(AF11+AF12+AF21+AF22)*HALF/DMOD

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

!  original code (Mishchenko and Travis)
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
!    We think that sine-terms need to be reversed! Entries 3, 7, 9, 10
!    Real Refractive index --> Entries 4, 8, 12, 13, 14, 15 are zero.

!  New code (Version 2.7, pending verification)

      GISSCOXMUNK_VKERNEL(3)  =    (CTTTP+CPTPP)*DCOEFF   ! Sine
      GISSCOXMUNK_VKERNEL(7)  =    (CTTTP-CPTPP)*DCOEFF   ! Sine
      GISSCOXMUNK_VKERNEL(9)  =    (CTTPT+CTPPP)*DCOEFF   ! Sine
      GISSCOXMUNK_VKERNEL(10) =    (CTTPT-CTPPP)*DCOEFF   ! Sine
      GISSCOXMUNK_VKERNEL(11) =    (CTTPP+CTPPT)*DCOEFF   ! Cos
      GISSCOXMUNK_VKERNEL(16) =    (CTTPP-CTPPT)*DCOEFF   ! Cos

!  Original VLIDORT code (Version 2.6)

!      GISSCOXMUNK_VKERNEL(3)  =    (-CTTTP-CPTPP)*DCOEFF
!      GISSCOXMUNK_VKERNEL(7)  =    (-CTTTP+CPTPP)*DCOEFF
!      GISSCOXMUNK_VKERNEL(9)  =    (-CTTPT-CTPPP)*DCOEFF
!      GISSCOXMUNK_VKERNEL(10) =    (-CTTPT+CTPPP)*DCOEFF
!      GISSCOXMUNK_VKERNEL(11) =    (-CTTPP-CTPPT)*DCOEFF ! Wrongly coded, this should be a Cos term
!      GISSCOXMUNK_VKERNEL(16) =    (-CTTPP+CTPPT)*DCOEFF ! Wrongly coded, this should be a Cos term

!  No Shadow code if not flagged

      IF ( PARS(3) .EQ. ZERO ) RETURN

!  Shadow code

      S1 = DSQRT(TWO*SIGMA2/PIE)
      S3 = ONE/(DSQRT(TWO*SIGMA2))
      S2 = S3*S3

      IF ( XI .EQ. ONE ) THEN
       SHADOWI   = ZERO
      ELSE
       XXI  = XI*XI
       DCOT = XI/DSQRT(ONE-XXI)
       T1   = DEXP(-DCOT*DCOT*S2)
!       T2   = DERFC(DCOT*S3)
       T2   = DERFC_E(DCOT*S3)
       SHADOWI = HALF*(S1*T1/DCOT-T2)
      ENDIF

      IF ( XJ .EQ. ONE ) THEN
       SHADOWR   = ZERO
      ELSE
       XXJ  = XJ*XJ
       DCOT = XJ/DSQRT(ONE-XXJ)
       T1   = DEXP(-DCOT*DCOT*S2)
!       T2   = DERFC(DCOT*S3)
       T2   = DERFC_E(DCOT*S3)
       SHADOWR = HALF*(S1*T1/DCOT-T2)
      ENDIF

      SHADOW = ONE/(ONE+SHADOWI+SHADOWR)

      DO I = 1, 10
       IM = KERNELMASK(I)
       GISSCOXMUNK_VKERNEL(IM) = GISSCOXMUNK_VKERNEL(IM) * SHADOW
      ENDDO

!  debug

!      DO M = 1, NSSQ
!        write(33,'(I5,1p6e14.5)')m,
!     &   GISSCOXMUNK_VKERNEL(m),
!     &   dacos(xi)/deg_to_rad, dacos(xj)/deg_to_rad, phi
!      enddo

!  Finish

      RETURN
      END SUBROUTINE GISSCOXMUNK_VFUNCTION

!

      SUBROUTINE GCMCRI_VFUNCTION &
         ( MAXPARS, NPARS, PARS, NSSQ, DO_SHADOW, XJ, SXJ, XI, SXI, &
           XPHI_REF, CKPHI_REF, SKPHI_REF, GISSCOXMUNK_VKERNEL )

!  include file of constants

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER ::          MAXPARS, NPARS, NSSQ
      LOGICAL ::          DO_SHADOW
      DOUBLE PRECISION :: PARS ( MAXPARS )
      DOUBLE PRECISION :: XI, XJ
      DOUBLE PRECISION :: SXI, SXJ, XPHI_REF, CKPHI_REF, SKPHI_REF
      DOUBLE PRECISION :: GISSCOXMUNK_VKERNEL(MAXSTOKES_SQ)

!  Critical exponent taken out

      DOUBLE PRECISION, PARAMETER :: CRITEXP = 88.0D0

!  local variables

      INTEGER ::          O1, KERNELMASK(16), I, IM
      DOUBLE PRECISION :: XPHI_INC, CKPHI_INC, SKPHI_INC
      DOUBLE PRECISION :: VI1, VI2, VI3, VR1, VR2, VR3
      DOUBLE PRECISION :: unit1, unit2, unit3, fact1, factor
      DOUBLE PRECISION :: XI1, TI1, TI2, TI3, TR1, TR2, TR3
      DOUBLE PRECISION :: PI1, PII2, PI3, PR1, PR2, PR3
      DOUBLE PRECISION :: PIKR, PRKI, TIKR, TRKI
      DOUBLE PRECISION :: E1, E2, E3, E4, SIGMA2
      DOUBLE PRECISION :: RDZ2, RDZ4
      DOUBLE PRECISION :: VP1, VP2, VP3, DMOD, DEX, DCOEFF
      DOUBLE PRECISION :: AF, AF11, AF12, AF21, AF22
      DOUBLE PRECISION :: S1, S2, S3, XXI, XXJ, T1, T2, DCOT
      DOUBLE PRECISION :: SHADOWI, SHADOWR, SHADOW, DCOEFF_0, ARGUMENT

      DOUBLE COMPLEX ::   CI, CF11, CF12, CF21, CF22, C21, C22
      DOUBLE COMPLEX ::   CN1, CN2, CXI2, C2, C1, CRPER, CRPAR
      DOUBLE COMPLEX ::   CTTTP, CTTPT, CTTPP, CTPPT, CTPPP, CPTPP

      KERNELMASK = (/1,2,3,5,6,7,9,10,11,4,8,12,13,14,15,16/)

!  Initialise

      DO O1 = 1, NSSQ
        GISSCOXMUNK_VKERNEL(O1) = ZERO
      ENDDO

!  Transcription of the RMATR subroutine from Mishchenko/Travis code.

!   CALCULATION OF THE STOKES REFLECTION MATRIX FOR
!   ILLUMINATION FROM ABOVE FOR
!   A STATISTICALLY ROUGH SURFACE SEPARATING TWO HALF-SPACES
!   WITH REFRACTIVE INDICES OF THE UPPER AND LOWER HALF-SPACES EQUAL TO
!   CN1 AND CN2, RESPECTIVELY. THE EFFECT OF SHADOWING IS NOT
!   INCLUDED IN THIS SUBROUTINE BUT IS ADDED IN THE MAIN PROGRAM.

!   SIGMA2 = s**2 = MEAN SQUARE SURFACE SLOPE (EQ. (18) IN THE JGR PAPER

!   XI = ABS(COSINE OF THE INCIDENT ZENITH ANGLE)
!   XJ = ABS(COSINE OF THE REFLECTION ZENITH ANGLE).
!   SXI and SXJ are the respective SINES (input)
!   XPHI_REF = REFLECTION AZIMUTH ANGLE
!   GISSCOXMUNK_VKERNEL(16-elements) = (4X4) REFLECTION MATRIX

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!               Remark on Use of shadow effect
!               ------------------------------
!  Shadow effect is controlled by the LOGICAL input flag. T
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  For real case, incident azimuth taken to be zero

      XPHI_INC  = ZERO
      CKPHI_INC = ONE
      SKPHI_INC = ZERO

!  Slope square is PARS(1)

      SIGMA2 = PARS(1)

!  refractive Index, Real + Imaginary, are PARS 2 and 3

!mick
      !CN1 = ( ONE, ZERO )
      CN1 = CMPLX ( ONE, ZERO )

      CN2 = CMPLX ( PARS(2), PARS(3) )

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

!   FRESNEL REFLECTION COEFFICIENTS
!   -------------------------------

!  this is the original code, but now C-variables are complex
!    CDSQRT is Obsolete (NAG Compiler complains)

      XI1 =  FACTOR*(UNIT1*VI1+UNIT2*VI2+UNIT3*VI3)
      CXI2 = ONE - (ONE-XI1*XI1)*CN1*CN1/(CN2*CN2)
!      CXI2 = CDSQRT(CXI2)
      CXI2 = SQRT(CXI2)
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

      AF  = HALF * DCOEFF

!  Alternative

      AF11 = real(cf11)*real(cf11) + aimag(cf11)*aimag(cf11)
      AF12 = real(cf12)*real(cf12) + aimag(cf12)*aimag(cf12)
      AF21 = real(cf21)*real(cf21) + aimag(cf21)*aimag(cf21)
      AF22 = real(cf22)*real(cf22) + aimag(cf22)*aimag(cf22)

!  original. CDABS is Obsolete  (NAG compiler complaints)

!      AF11 = CDABS(CF11)
!      AF12 = CDABS(CF12)
!      AF21 = CDABS(CF21)
!      AF22 = CDABS(CF22)
!      AF11 = AF11*AF11
!      AF12 = AF12*AF12
!      AF21 = AF21*AF21
!      AF22 = AF22*AF22

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

!  Key Debug statement to track down DMOD
!      write(*,*)'(1,1) Fresnel = ',(AF11+AF12+AF21+AF22)*HALF/DMOD

!  Original code. DCONJG is Obsolete (NAG compiler does not like)

!mick
      !CI  = ( ZERO, MINUS_ONE )
      CI  = CMPLX ( ZERO, MINUS_ONE )

!      C21 = DCONJG(CF21)
!      C22 = DCONJG(CF22)
!      CTTTP = CF11*DCONJG(CF12)
      C21 = CONJG(CF21)
      C22 = CONJG(CF22)
      CTTTP = CF11*CONJG(CF12)
      CTTPT = CF11*C21
      CTTPP = CF11*C22
      CTPPT = CF12*C21
      CTPPP = CF12*C22
      CPTPP = CF21*C22

!  replica for real variables only
!      C21 = CF21
!      C22 = CF22
!      CTTTP=CF11*CF12
!      CTTPT=CF11*C21
!      CTTPP=CF11*C22
!      CTPPT=CF12*C21
!      CTPPP=CF12*C22
!      CPTPP=CF21*C22

!  original code (Mishchenko and Travis)
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

!  Original VLIDORT code (Version 2.6)

!      GISSCOXMUNK_VKERNEL(3)  =    (-CTTTP-CPTPP)*DCOEFF
!      GISSCOXMUNK_VKERNEL(7)  =    (-CTTTP+CPTPP)*DCOEFF
!      GISSCOXMUNK_VKERNEL(9)  =    (-CTTPT-CTPPP)*DCOEFF
!      GISSCOXMUNK_VKERNEL(10) =    (-CTTPT+CTPPP)*DCOEFF
!      GISSCOXMUNK_VKERNEL(11) =    ( CTTPP+CTPPT)*DCOEFF
!        GISSCOXMUNK_VKERNEL(4)  = - CI * ( CTTTP+CPTPP)*DCOEFF
!        GISSCOXMUNK_VKERNEL(8)  = - CI * ( CTTTP-CPTPP)*DCOEFF
!        GISSCOXMUNK_VKERNEL(12) =   CI * ( CTTPP-CTPPT)*DCOEFF
!        GISSCOXMUNK_VKERNEL(13) =   CI * ( CTTPT+CTPPP)*DCOEFF
!        GISSCOXMUNK_VKERNEL(14) =   CI * ( CTTPT-CTPPP)*DCOEFF
!        GISSCOXMUNK_VKERNEL(15) = - CI * ( CTTPP+CTPPT)*DCOEFF
!        GISSCOXMUNK_VKERNEL(16) =        ( CTTPP-CTPPT)*DCOEFF

!  New code (Version 2.7, pending verification)
!    We think that sine-terms need to be reversed! Entries 3, 4, 7, 8, 9, 10, 13, 14

      GISSCOXMUNK_VKERNEL(3)  =    ( CTTTP+CPTPP)*DCOEFF   ! Sine
      GISSCOXMUNK_VKERNEL(7)  =    ( CTTTP-CPTPP)*DCOEFF   ! Sine
      GISSCOXMUNK_VKERNEL(9)  =    ( CTTPT+CTPPP)*DCOEFF   ! Sine
      GISSCOXMUNK_VKERNEL(10) =    ( CTTPT-CTPPP)*DCOEFF   ! Sine
      GISSCOXMUNK_VKERNEL(11) =    ( CTTPP+CTPPT)*DCOEFF   ! Cos
      IF ( NSSQ.EQ.16) THEN
        GISSCOXMUNK_VKERNEL(4)  = + CI * ( CTTTP+CPTPP)*DCOEFF  ! Sine
        GISSCOXMUNK_VKERNEL(8)  = + CI * ( CTTTP-CPTPP)*DCOEFF  ! Sine
        GISSCOXMUNK_VKERNEL(12) =   CI * ( CTTPP-CTPPT)*DCOEFF  ! Cos
        GISSCOXMUNK_VKERNEL(13) = - CI * ( CTTPT+CTPPP)*DCOEFF  ! Sine
        GISSCOXMUNK_VKERNEL(14) = - CI * ( CTTPT-CTPPP)*DCOEFF  ! Sine
        GISSCOXMUNK_VKERNEL(15) = - CI * ( CTTPP+CTPPT)*DCOEFF  ! Cos
        GISSCOXMUNK_VKERNEL(16) =        ( CTTPP-CTPPT)*DCOEFF  ! Cos
      ENDIF

!  No Shadow code if not flagged

      IF ( .NOT. DO_SHADOW ) RETURN

!  Shadow code

      S1 = DSQRT(TWO*SIGMA2/PIE)
      S3 = ONE/(DSQRT(TWO*SIGMA2))
      S2 = S3*S3

      IF ( XI .EQ. ONE ) THEN
       SHADOWI   = ZERO
      ELSE
       XXI  = XI*XI
       DCOT = XI/DSQRT(ONE-XXI)
       T1   = DEXP(-DCOT*DCOT*S2)
!       T2   = DERFC(DCOT*S3)
       T2   = DERFC_E(DCOT*S3)
       SHADOWI = HALF*(S1*T1/DCOT-T2)
      ENDIF

      IF ( XJ .EQ. ONE ) THEN
       SHADOWR   = ZERO
      ELSE
       XXJ  = XJ*XJ
       DCOT = XJ/DSQRT(ONE-XXJ)
       T1   = DEXP(-DCOT*DCOT*S2)
!       T2   = DERFC(DCOT*S3)
       T2   = DERFC_E(DCOT*S3)
       SHADOWR = HALF*(S1*T1/DCOT-T2)
      ENDIF

      SHADOW = ONE/(ONE+SHADOWI+SHADOWR)

      DO I = 1, NSSQ
       IM = KERNELMASK(I)
       GISSCOXMUNK_VKERNEL(IM) = GISSCOXMUNK_VKERNEL(IM) * SHADOW
      ENDDO

!  debug

!      DO M = 1, NSSQ
!        write(33,'(I5,1p6e14.5)')m, GISSCOXMUNK_VKERNEL(m), &
!         dacos(xi)/deg_to_rad, dacos(xj)/deg_to_rad, phi
!      enddo

!  Finish

      RETURN
      END SUBROUTINE GCMCRI_VFUNCTION

!

      SUBROUTINE COXMUNK_VFUNCTION_MSR &
     ( MAXPARS, NPARS, PARS, ORDER, NSSQ,                                &
       n_muquad, n_phiquad, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,          & 
       X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,  &
       COXMUNK_VKERNEL )

!  include file of constants

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER ::          MAXPARS, NPARS, NSSQ, ORDER
      DOUBLE PRECISION :: PARS ( MAXPARS )
      DOUBLE PRECISION :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION :: COXMUNK_VKERNEL(MAXSTOKES_SQ)

!  Local arrays for MSR quadrature

      INTEGER          :: n_muquad, n_phiquad
      DOUBLE PRECISION :: X_MUQUAD (max_msrs_muquad)
      DOUBLE PRECISION :: W_MUQUAD (max_msrs_muquad)
      DOUBLE PRECISION :: SX_MUQUAD (max_msrs_muquad)
      DOUBLE PRECISION :: WXX_MUQUAD (max_msrs_muquad)

      DOUBLE PRECISION :: X_PHIQUAD (max_msrs_phiquad)
      DOUBLE PRECISION :: W_PHIQUAD (max_msrs_phiquad)

!  local variables

      INTEGER         ::  s, n, k, i, i1, N_phiquad_HALF, o1, ni, ki, nr, kr
      DOUBLE PRECISION :: XM, SXM, XMR, SXMR, XMI, SXMI, sum_pr, sumr, sum, w_p
      DOUBLE PRECISION :: reflec_0(MAXSTOKES_SQ), reflec_s, R0Q, reflec
      DOUBLE PRECISION :: phi_sub1, cphi_sub1, sphi_sub1
      DOUBLE PRECISION :: phi_sub2, cphi_sub2, sphi_sub2

!  arrays

      DOUBLE PRECISION :: R0_QUAD_IN  (16,max_msrs_muquad,max_msrs_phiquad)
      DOUBLE PRECISION :: R0_OUT_QUAD (16,max_msrs_muquad,max_msrs_phiquad)
      DOUBLE PRECISION :: RHOLD       (1, max_msrs_muquad,max_msrs_phiquad)
      DOUBLE PRECISION :: R0_MSRS_QUAD(16,max_msrs_muquad,max_msrs_phiquad)

!  Safety first zeroing

      DO O1 = 1, NSSQ
        REFLEC_0(O1) = ZERO
!        REFLEC_1(O1) = ZERO
        COXMUNK_VKERNEL(O1) = ZERO
      ENDDO

!  Only want the first element (this is a scalar routine)

!  Single scattering (zero order), Phi is in degrees here!

      CALL COXMUNK_VFUNCTION &
         ( MAXPARS, NPARS, PARS, NSSQ, XJ, SXJ, XI, SXI, &
           PHI, CPHI, SKPHI, REFLEC_0 )

!  Higher orders scattering
!    Only want the first element (this is a scalar routine)

      REFLEC   = REFLEC_0(1)
      REFLEC_S = ZERO

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
               CALL COXMUNK_VFUNCTION &
                ( MAXPARS, NPARS, PARS, NSSQ, XI, SXI, XM, SXM, PHI_SUB2, &
                  CPHI_SUB2, SPHI_SUB2, R0_OUT_QUAD(1,K,N) )
               CALL COXMUNK_VFUNCTION &
                ( MAXPARS, NPARS, PARS, NSSQ, XM, SXM, XJ, SXJ, PHI_SUB1, &
                  CPHI_SUB1, SPHI_SUB1, R0_QUAD_IN(1,K,N) )
            ENDDO
         ENDDO
      ENDIF

!  Compute the successive orders of scattering.
!  Compute higher orders - (1,1) component only.

      DO S = 1, ORDER

!  Compute result for this order

         SUMR = ZERO
         DO K = 1, n_muquad
            SUM_PR = ZERO
            DO N = 1, N_PHIQUAD
               W_P = W_PHIQUAD(N)
               SUM = R0_QUAD_IN(1,K,N) * R0_OUT_QUAD(1,K,N)
               SUM_PR = SUM_PR + W_P * SUM
            ENDDO
            SUMR = SUMR + SUM_PR * WXX_MUQUAD(K)
         ENDDO
         REFLEC_S = SUMR

!  Finish if reached the scattering order desired

         IF ( S.EQ.ORDER ) GO TO 67

!  Compute reflectance for next order and update
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
                     CALL COXMUNK_VFUNCTION &
                       ( MAXPARS, NPARS, PARS, NSSQ, &
                         XMI, SXMI, XMR, SXMR, PHI_SUB1, CPHI_SUB1, SPHI_SUB1,&
                         R0_MSRS_QUAD(1,KI,NI))
                  ENDDO
               ENDDO

!  Multiple reflection

               SUMR = ZERO
               DO KI = 1, n_muquad
                  SUM_PR = ZERO
                  DO NI = 1, N_PHIQUAD
                     W_P = W_PHIQUAD(NI)
                     R0Q = R0_MSRS_QUAD(1,KI,NI)
                     SUM_PR = SUM_PR + W_P * R0_QUAD_IN(1,KI,NI) * R0Q
                  ENDDO
!                  SUMR = SUMR + SUM_PR * W_MUQUAD(KI)
                  SUMR = SUMR + SUM_PR * WXX_MUQUAD(KI)
               ENDDO
               RHOLD(1,KR,NR) = SUMR

!  End KR, NR loops

            ENDDO
         ENDDO

!  Update

         DO KR = 1, n_muquad
            DO NR = 1, N_PHIQUAD
               R0_QUAD_IN(1,KR,NR) = RHOLD(1,KR,NR)
            ENDDO
         ENDDO

!  Continuation point for finishing MSR

 67      CONTINUE

!  Add to total

          REFLEC = REFLEC + REFLEC_S

!  End scattering order loop

      ENDDO

!  Compute total

      COXMUNK_VKERNEL(1) = REFLEC

!  debug

!      DO M = 1, NSSQ
!        write(34,'(I5,1p4e14.5)')m, COXMUNK_VKERNEL(m), &
!          dacos(xi)/deg_to_rad, dacos(xj)/deg_to_rad, phi
!      ENDDO

!  Finish

      RETURN
      END SUBROUTINE COXMUNK_VFUNCTION_MSR

!

      SUBROUTINE GISSCOXMUNK_VFUNCTION_MSR &
     ( MAXPARS, NPARS, PARS, ORDER, NSSQ,                                &
       n_muquad, n_phiquad, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,          & 
       X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,  &
       GISSCOXMUNK_VKERNEL )

!  include file of constants

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER ::          MAXPARS, NPARS, NSSQ, ORDER
      DOUBLE PRECISION :: PARS ( MAXPARS )
      DOUBLE PRECISION :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION :: GISSCOXMUNK_VKERNEL(MAXSTOKES_SQ)

!  Local arrays for MSR quadrature

      INTEGER          :: n_muquad, n_phiquad
      DOUBLE PRECISION :: X_MUQUAD (max_msrs_muquad)
      DOUBLE PRECISION :: W_MUQUAD (max_msrs_muquad)
      DOUBLE PRECISION :: SX_MUQUAD (max_msrs_muquad)
      DOUBLE PRECISION :: WXX_MUQUAD (max_msrs_muquad)

      DOUBLE PRECISION :: X_PHIQUAD (max_msrs_phiquad)
      DOUBLE PRECISION :: W_PHIQUAD (max_msrs_phiquad)

!  local variables

      INTEGER          :: s, n, k, i, i1, N_phiquad_HALF, o1, o2, o3, ni, ki, nr, kr
      DOUBLE PRECISION :: XM, SXM, XMR, SXMR, XMI, SXMI, sum_pr(16), sumr(16), sum, w_p
      DOUBLE PRECISION :: reflec_0(MAXSTOKES_SQ), reflec(16), reflec_s(16),R0Q
      DOUBLE PRECISION :: phi_sub1, cphi_sub1, sphi_sub1
      DOUBLE PRECISION :: phi_sub2, cphi_sub2, sphi_sub2

!  arrays

      DOUBLE PRECISION :: R0_QUAD_IN  (16,max_msrs_muquad,max_msrs_phiquad)
      DOUBLE PRECISION :: R0_OUT_QUAD (16,max_msrs_muquad,max_msrs_phiquad)
      DOUBLE PRECISION :: RHOLD       (16,max_msrs_muquad,max_msrs_phiquad)
      DOUBLE PRECISION :: R0_MSRS_QUAD(16,max_msrs_muquad,max_msrs_phiquad)

!  Indices

      INTEGER ::          MASKIT(3,3), LNS, NS, NSS, M, M12, M13, M32

!  Safety first zeroing

      DO O1 = 1, NSSQ
        REFLEC_0(O1) = ZERO
        REFLEC(O1)   = ZERO
        GISSCOXMUNK_VKERNEL(O1) = ZERO
      ENDDO

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

      CALL GISSCOXMUNK_VFUNCTION &
         ( MAXPARS, NPARS, PARS, NSSQ, XJ, SXJ, XI, SXI, &
           PHI, CPHI, SKPHI, REFLEC_0 )

!  Higher orders scattering
!    Only want the first element (this is a scalar routine)

      REFLEC   = REFLEC_0
      REFLEC_S = ZERO

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
               CALL GISSCOXMUNK_VFUNCTION &
                ( MAXPARS, NPARS, PARS, NSSQ, XI, SXI, XM, SXM, PHI_SUB2, &
                  CPHI_SUB2, SPHI_SUB2, R0_OUT_QUAD(1,K,N) )
               CALL GISSCOXMUNK_VFUNCTION &
                ( MAXPARS, NPARS, PARS, NSSQ, XM, SXM, XJ, SXJ, PHI_SUB1, &
                  CPHI_SUB1, SPHI_SUB1, R0_QUAD_IN(1,K,N) )
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
                     CALL GISSCOXMUNK_VFUNCTION &
                       ( MAXPARS, NPARS, PARS, NSSQ, &
                         XMI, SXMI, XMR, SXMR, PHI_SUB1, CPHI_SUB1, SPHI_SUB1,&
                         R0_MSRS_QUAD(1,KI,NI))
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

!  End loop over KR and NR

            ENDDO
         ENDDO     

!  Update

         DO KR = 1, N_MUQUAD
            DO NR = 1, N_PHIQUAD
               R0_QUAD_IN(1:NSSQ,KR,NR) = RHOLD(1:NSSQ,KR,NR)
            ENDDO
         ENDDO

!  Continuation point for finishing MSR
  
 67      continue

!  Add to total

         REFLEC(1:NSSQ) = REFLEC(1:NSSQ) + REFLEC_S(1:NSSQ)

!  End scattering order loop

      ENDDO

!  Compute total

      GISSCOXMUNK_VKERNEL(1:NSSQ) = REFLEC(1:NSSQ)

!  debug

!      DO M = 1, NSSQ
!        write(34,'(I5,1p6e14.5)')m,                      &
!        reflec_0(m), reflec_1(m),GISSCOXMUNK_VKERNEL(m), &
!        dacos(xi)/deg_to_rad, dacos(xj)/deg_to_rad, phi
!      enddo

!  Finish

      RETURN
      END SUBROUTINE GISSCOXMUNK_VFUNCTION_MSR

!

      SUBROUTINE GCMCRI_VFUNCTION_MSR &
     ( MAXPARS, NPARS, PARS, ORDER, NSSQ, DO_SHADOW,                    &
       n_muquad, n_phiquad, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,          & 
       X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, X_PHIQUAD, W_PHIQUAD,  &
       GISSCOXMUNK_VKERNEL )

!  include file of constants

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER ::          MAXPARS, NPARS, NSSQ, ORDER
      LOGICAL ::          DO_SHADOW
      DOUBLE PRECISION :: PARS ( MAXPARS )
      DOUBLE PRECISION :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION :: GISSCOXMUNK_VKERNEL(MAXSTOKES_SQ)

!  Local arrays for MSR quadrature

      INTEGER          :: n_muquad, n_phiquad
      DOUBLE PRECISION :: X_MUQUAD (max_msrs_muquad)
      DOUBLE PRECISION :: W_MUQUAD (max_msrs_muquad)
      DOUBLE PRECISION :: SX_MUQUAD (max_msrs_muquad)
      DOUBLE PRECISION :: WXX_MUQUAD (max_msrs_muquad)

      DOUBLE PRECISION :: X_PHIQUAD (max_msrs_phiquad)
      DOUBLE PRECISION :: W_PHIQUAD (max_msrs_phiquad)

!  local variables

      INTEGER          :: s, n, k, i, i1, N_phiquad_HALF, o1, o2, o3, ni, ki, nr, kr
      DOUBLE PRECISION :: XM, SXM, XMR, SXMR, XMI, SXMI, sum_pr(16), sumr(16), sum, w_p
      DOUBLE PRECISION :: reflec_0(MAXSTOKES_SQ), reflec(16), reflec_S(16), R0Q
      DOUBLE PRECISION :: phi_sub1, cphi_sub1, sphi_sub1
      DOUBLE PRECISION :: phi_sub2, cphi_sub2, sphi_sub2

!  arrays

      DOUBLE PRECISION :: R0_QUAD_IN  (16,max_msrs_muquad,max_msrs_phiquad)
      DOUBLE PRECISION :: R0_OUT_QUAD (16,max_msrs_muquad,max_msrs_phiquad)
      DOUBLE PRECISION :: RHOLD       (16,max_msrs_muquad,max_msrs_phiquad)
      DOUBLE PRECISION :: R0_MSRS_QUAD(16,max_msrs_muquad,max_msrs_phiquad)

!  Indices

      INTEGER ::          MASKIT(4,4), LNS, NS, NSS, M, M12, M13, M32

!  Safety first zeroing

      DO O1 = 1, NSSQ
        REFLEC_0(O1) = ZERO
        REFLEC(O1)   = ZERO
        GISSCOXMUNK_VKERNEL(O1) = ZERO
      ENDDO

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

      MASKIT(1,4) = 4
      MASKIT(2,4) = 8
      MASKIT(3,4) = 12
      MASKIT(4,1) = 13
      MASKIT(4,2) = 14
      MASKIT(4,3) = 15
      MASKIT(4,4) = 16

!  Single scattering (zero order), Phi is in degrees here!

      CALL GCMCRI_VFUNCTION &
         ( MAXPARS, NPARS, PARS, NSSQ, DO_SHADOW, &
           XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI, &
           REFLEC_0 )

!  Higher orders scattering
!    Only want the first element (this is a scalar routine)

      REFLEC   = REFLEC_0
      REFLEC_S = ZERO

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
               CALL GCMCRI_VFUNCTION &
                ( MAXPARS, NPARS, PARS, NSSQ, DO_SHADOW, XI, SXI, XM, SXM, PHI_SUB2, &
                  CPHI_SUB2, SPHI_SUB2, R0_OUT_QUAD(1,K,N) )
               CALL GCMCRI_VFUNCTION &
                ( MAXPARS, NPARS, PARS, NSSQ, DO_SHADOW, XM, SXM, XJ, SXJ, PHI_SUB1, &
                  CPHI_SUB1, SPHI_SUB1, R0_QUAD_IN(1,K,N) )
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
                     CALL GCMCRI_VFUNCTION &
                       ( MAXPARS, NPARS, PARS, NSSQ, DO_SHADOW, &
                         XMI, SXMI, XMR, SXMR, PHI_SUB1, CPHI_SUB1, SPHI_SUB1,&
                         R0_MSRS_QUAD(1,KI,NI))
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

!  End loop over KR and NR

            ENDDO
         ENDDO     

!  Update

         DO KR = 1, N_MUQUAD
            DO NR = 1, N_PHIQUAD
               R0_QUAD_IN(1:NSSQ,KR,NR) = RHOLD(1:NSSQ,KR,NR)
            ENDDO
         ENDDO

!  Continuation point for finishing MSR
  
 67      continue

!  Add to total

          REFLEC(1:NSSQ) = REFLEC(1:NSSQ) + REFLEC_S(1:NSSQ)

!  End scattering order loop

      ENDDO

!  Compute total

      GISSCOXMUNK_VKERNEL(1:NSSQ) = REFLEC(1:NSSQ) 

!  debug

!      DO M = 1, NSSQ
!        write(34,'(I5,1p6e14.5)')m, &
!         reflec_0(m), reflec_1(m),GISSCOXMUNK_VKERNEL(m), &
!         dacos(xi)/deg_to_rad, dacos(xj)/deg_to_rad, phi
!      enddo

!  Finish

      RETURN
      END SUBROUTINE GCMCRI_VFUNCTION_MSR

!

      SUBROUTINE BPDFSOIL_VFUNCTION  &
        ( MAXPARS, NPARS, PARS, NSSQ, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, BPDFSOIL_VKERNEL )

!  module, dimensions and numbers

      USE VLIDORT_pars, only : ZERO, ONE, HALF, PIE

      implicit none

!  Subroutine arguments

      INTEGER         , intent(in)     :: MAXPARS, NPARS, NSSQ
      DOUBLE PRECISION, intent(in)     :: PARS ( MAXPARS )
      DOUBLE PRECISION, intent(inout)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION, intent(out)    :: BPDFSOIL_VKERNEL(16)

!  local variables

      DOUBLE PRECISION :: DMOD, CF11, CF12, CF21, CF22, CN1, CN2
      DOUBLE PRECISION :: AF11, AF12, AF21, AF22
      DOUBLE PRECISION :: C21, C22, CTTTP, CTTPT, CTTPP
      DOUBLE PRECISION :: CTPPT, CTPPP, CPTPP

!  H-function variables

      DOUBLE PRECISION :: HFUNCTION, FACTOR
      DOUBLE PRECISION :: ATTEN, Z, Z2, Z1
      DOUBLE PRECISION :: sgamma, cgamma, calpha, calpha_sq, salpha

!  F-.M. Breon BPDF SOIL model (2009).
!   This is the Vector model with polarization

!  Initialise

      BPDFSOIL_VKERNEL = ZERO

!  Transcription of the RMATR subroutine from Mishchenko/Travis code.
!   CALCULATION OF THE STOKES REFLECTION MATRIX FOR
!   ILLUMINATION FROM ABOVE FOR A SURFACE SEPARATING TWO HALF-SPACES
!   WITH REFRACTIVE INDICES OF THE UPPER AND LOWER HALF-SPACES EQUAL TO
!   CN1 AND CN2, RESPECTIVELY.
!   XI = ABS(COSINE OF THE INCIDENT ZENITH ANGLE)
!   XJ = ABS(COSINE OF THE REFLECTION ZENITH ANGLE).
!   SXI and SXJ are the respective SINES (input)
!   XPHI_REF = REFLECTION AZIMUTH ANGLE

!  PARS(1) = refractive index of water (real)

      CN1 = ONE
      CN2 = PARS(1)

!  Check for limiting cases

      IF(DABS(XI-1D0).LT.1d-9) XI = 0.999999999999d0
      IF(DABS(XJ-1D0).LT.1d-9) XJ = 0.999999999999d0

!   Angle of the surface that generates specular reflection from
!  sun to view directions (theta)
!      alpha = DACOS(HALF*(mus+muv)/dcos(gamma))

      Z = XI * XJ - SXI * SXJ * CPHI
      IF ( Z .GT. ONE) Z = ONE
      Z1 = DACOS(Z)
      Z2 = DCOS(Z1*HALF)

!   Angle of the surface that generates specular reflection from 
!  sun to view directions (theta)
!      alpha = DACOS(HALF*(mus+muv)/dcos(gamma))

      calpha    = HALF * (xi + xj) / Z2  
      calpha_sq = calpha*calpha
      salpha    = dsqrt(one - calpha_sq)

! attenuation factor

      cgamma = Z2
      sgamma = dsqrt ( one - cgamma * cgamma )
      atten  = one - sgamma

!  Final H-function

      HFUNCTION = 0.25d0 * atten / xi / xj

!  Call to the vector Fresnel routine

      CALL FRESNEL_VECTOR &
      ( CN1, CN2, XI, SXI, XJ, SXJ, CPHI, SKPHI, & ! Input
        CF11, CF12, CF21, CF22, DMOD )             ! Output

!  Final H-function needs to be normalized

      HFUNCTION = HFUNCTION / DMOD

!  Settting (1,1), (1,2), (2,1) and (2,2) components
!  -------------------------------------------------

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

      BPDFSOIL_VKERNEL(1) = (AF11+AF12+AF21+AF22) * FACTOR
      BPDFSOIL_VKERNEL(2) = (AF11-AF12+AF21-AF22) * FACTOR
      BPDFSOIL_VKERNEL(5) = (AF11-AF22+AF12-AF21) * FACTOR
      BPDFSOIL_VKERNEL(6) = (AF11-AF12-AF21+AF22) * FACTOR

!  Finish if NStokes = 1 or 2 (NSSQ = 1 or 4)

      if ( NSSQ .le. 4 ) return

!  Setting (1,3), (2,3), (3,1), (3,2), (3,3) and (4,4) components
!  --------------------------------------------------------------

      C21 = CF21
      C22 = CF22
      CTTTP=CF11*CF12
      CTTPT=CF11*C21
      CTTPP=CF11*C22
      CTPPT=CF12*C21
      CTPPP=CF12*C22
      CPTPP=CF21*C22

      FACTOR = HFUNCTION

!  New code (Version 2.7, pending verification)
!    We think that sine-terms need to be reversed! Entries 3, 7, 9, 10.

      BPDFSOIL_VKERNEL(3)  =    ( CTTTP+CPTPP) * FACTOR   ! Sine
      BPDFSOIL_VKERNEL(7)  =    ( CTTTP-CPTPP) * FACTOR   ! Sine
      BPDFSOIL_VKERNEL(9)  =    ( CTTPT+CTPPP) * FACTOR   ! Sine
      BPDFSOIL_VKERNEL(10) =    ( CTTPT-CTPPP) * FACTOR   ! Sine
      BPDFSOIL_VKERNEL(11) =    ( CTTPP+CTPPT) * FACTOR   ! Cos
      IF ( NSSQ.eq.16 ) then
         BPDFSOIL_VKERNEL(16) =    ( CTTPP-CTPPT) * FACTOR   ! Cos
      endif

!  Original VLIDORT code (Version 2.6)

!      BPDF2009_VKERNEL(3)  =    (-CTTTP-CPTPP) * FACTOR
!      BPDF2009_VKERNEL(7)  =    (-CTTTP+CPTPP) * FACTOR
!      BPDF2009_VKERNEL(9)  =    (-CTTPT-CTPPP) * FACTOR
!      BPDF2009_VKERNEL(10) =    (-CTTPT+CTPPP) * FACTOR
!      BPDF2009_VKERNEL(11) =    (-CTTPP-CTPPT) * FACTOR   !  Wrongly coded, these are Cosine terms
!      BPDF2009_VKERNEL(16) =    (-CTTPP+CTPPT) * FACTOR   !  Wrongly coded, these are Cosine terms

!  Finish

      RETURN
      END SUBROUTINE BPDFSOIL_VFUNCTION

!

      SUBROUTINE BPDFVEGN_VFUNCTION  &
        ( MAXPARS, NPARS, PARS, NSSQ, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, BPDFVEGN_VKERNEL )

!  module, dimensions and numbers

      USE VLIDORT_pars, only : ZERO, ONE, HALF, MINUS_ONE, PIE

      implicit none

!  Subroutine arguments

      INTEGER         , intent(in)     :: MAXPARS, NPARS, NSSQ
      DOUBLE PRECISION, intent(in)     :: PARS ( MAXPARS )
      DOUBLE PRECISION, intent(inout)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION, intent(out)    :: BPDFVEGN_VKERNEL(16)

!  local variables

      DOUBLE PRECISION :: DMOD, CF11, CF12, CF21, CF22, CN1, CN2
      DOUBLE PRECISION :: AF11, AF12, AF21, AF22
      DOUBLE PRECISION :: C21, C22, CTTTP, CTTPT, CTTPP
      DOUBLE PRECISION :: CTPPT, CTPPP, CPTPP

!  H-function variables

      DOUBLE PRECISION :: HFUNCTION, FACTOR
      DOUBLE PRECISION :: ATTEN, FP0, Z, Z2, Z1
      DOUBLE PRECISION :: sgamma, cgamma, calpha, calpha_sq, salpha
      DOUBLE PRECISION :: PLEAF, GS, GV, PROJECTIONS

!  Data coefficients

      DOUBLE PRECISION :: PLAGIOPHILE_COEFFS(4)
      DATA PLAGIOPHILE_COEFFS /0.43181098d0,  0.011187479d0, &
                               0.043329567d0, 0.19262991d0/

!  F-.M. Breon BPDF VEGN model (2009).
!   This is the Vector model with polarization

!  Initialise

      BPDFVEGN_VKERNEL = ZERO

!  Transcription of the RMATR subroutine from Mishchenko/Travis code.
!   CALCULATION OF THE STOKES REFLECTION MATRIX FOR
!   ILLUMINATION FROM ABOVE FOR A SURFACE SEPARATING TWO HALF-SPACES
!   WITH REFRACTIVE INDICES OF THE UPPER AND LOWER HALF-SPACES EQUAL TO
!   CN1 AND CN2, RESPECTIVELY.
!   XI = ABS(COSINE OF THE INCIDENT ZENITH ANGLE)
!   XJ = ABS(COSINE OF THE REFLECTION ZENITH ANGLE).
!   SXI and SXJ are the respective SINES (input)
!   XPHI_REF = REFLECTION AZIMUTH ANGLE

!  PARS(1) = refractive index of water (real)

      CN1 = ONE
      CN2 = PARS(1)

!  Check for limiting cases

      IF(DABS(XI-1D0).LT.1d-9) XI = 0.999999999999d0
      IF(DABS(XJ-1D0).LT.1d-9) XJ = 0.999999999999d0

!   Angle of the surface that generates specular reflection from
!  sun to view directions (theta)
!      alpha = DACOS(HALF*(mus+muv)/dcos(gamma))

      Z = XI * XJ - SXI * SXJ * CPHI
      IF ( Z .GT. ONE) Z = ONE
      Z1 = DACOS(Z)
      Z2 = DCOS(Z1*HALF)

!   Angle of the surface that generates specular reflection from 
!  sun to view directions (theta)
!      alpha = DACOS(HALF*(mus+muv)/dcos(gamma))

      calpha    = HALF * (xi + xj) / Z2
      calpha_sq = calpha*calpha
      salpha    = sqrt(one - calpha_sq)

! Projection of leaf surface to outgoing direction

      gv = PLAGIOPHILE_COEFFS(1) + xi * &
          (PLAGIOPHILE_COEFFS(2) + xi * &
          (PLAGIOPHILE_COEFFS(3) + PLAGIOPHILE_COEFFS(4)*xi))

! Projection of leaf surface to incident direction

      gs = PLAGIOPHILE_COEFFS(1) + xj * &
          (PLAGIOPHILE_COEFFS(2) + xj * &
          (PLAGIOPHILE_COEFFS(3) + PLAGIOPHILE_COEFFS(4)*xj))

! Probability of leaf orientation (plagiophile distr.)

      Pleaf = 16.0_fpk * calpha_sq * salpha  / pie

! Polarization model for vegetation

      PROJECTIONS =  Gv/xi + Gs/xj
      Fp0 = 0.25_fpk * PLEAF  / xi / xj / PROJECTIONS

! attenuation factor

      cgamma = Z2
      sgamma = dsqrt ( one - cgamma * cgamma )
      atten  = one - sgamma

!  Final H-function

      HFUNCTION = Fp0 * atten

!  Call to the vector Fresnel routine

      CALL FRESNEL_VECTOR &
      ( CN1, CN2, XI, SXI, XJ, SXJ, CPHI, SKPHI, & ! Input
        CF11, CF12, CF21, CF22, DMOD )             ! Output

!  Final H-function needs to be normalized

      HFUNCTION = HFUNCTION / DMOD

!  Settting (1,1), (1,2), (2,1) and (2,2) components
!  -------------------------------------------------

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

      BPDFVEGN_VKERNEL(1) = (AF11+AF12+AF21+AF22) * FACTOR
      BPDFVEGN_VKERNEL(2) = (AF11-AF12+AF21-AF22) * FACTOR
      BPDFVEGN_VKERNEL(5) = (AF11-AF22+AF12-AF21) * FACTOR
      BPDFVEGN_VKERNEL(6) = (AF11-AF12-AF21+AF22) * FACTOR

!  Finish if NStokes = 1 or 2 (NSSQ = 1 or 4)

      if ( NSSQ .le. 4 ) return

!  Setting (1,3), (2,3), (3,1), (3,2), (3,3) and (4,4) components
!  --------------------------------------------------------------

      C21 = CF21
      C22 = CF22
      CTTTP=CF11*CF12
      CTTPT=CF11*C21
      CTTPP=CF11*C22
      CTPPT=CF12*C21
      CTPPP=CF12*C22
      CPTPP=CF21*C22

      FACTOR = HFUNCTION

!  New code (Version 2.7, pending verification)
!    We think that sine-terms need to be reversed! Entries 3, 7, 9, 10.

      BPDFVEGN_VKERNEL(3)  =    ( CTTTP+CPTPP) * FACTOR   ! Sine
      BPDFVEGN_VKERNEL(7)  =    ( CTTTP-CPTPP) * FACTOR   ! Sine
      BPDFVEGN_VKERNEL(9)  =    ( CTTPT+CTPPP) * FACTOR   ! Sine
      BPDFVEGN_VKERNEL(10) =    ( CTTPT-CTPPP) * FACTOR   ! Sine
      BPDFVEGN_VKERNEL(11) =    ( CTTPP+CTPPT) * FACTOR   ! Cos
      IF ( NSSQ.eq.16 ) then
         BPDFVEGN_VKERNEL(16) =    ( CTTPP-CTPPT) * FACTOR   ! Cos
      endif

!  Original VLIDORT code (Version 2.6)

!      BPDF2009_VKERNEL(3)  =    (-CTTTP-CPTPP) * FACTOR
!      BPDF2009_VKERNEL(7)  =    (-CTTTP+CPTPP) * FACTOR
!      BPDF2009_VKERNEL(9)  =    (-CTTPT-CTPPP) * FACTOR
!      BPDF2009_VKERNEL(10) =    (-CTTPT+CTPPP) * FACTOR
!      BPDF2009_VKERNEL(11) =    (-CTTPP-CTPPT) * FACTOR   !  Wrongly coded, these are Cosine terms
!      BPDF2009_VKERNEL(16) =    (-CTTPP+CTPPT) * FACTOR   !  Wrongly coded, these are Cosine terms

!  Finish

      RETURN
      END SUBROUTINE BPDFVEGN_VFUNCTION

!

      SUBROUTINE BPDFNDVI_VFUNCTION  &
        ( MAXPARS, NPARS, PARS, NSSQ, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, BPDFNDVI_VKERNEL )

!  module, dimensions and numbers

      USE VLIDORT_pars, only : ZERO, ONE, HALF, MINUS_ONE, PIE

      implicit none

!  Subroutine arguments

      INTEGER         , intent(in)     :: MAXPARS, NPARS, NSSQ
      DOUBLE PRECISION, intent(in)     :: PARS ( MAXPARS )
      DOUBLE PRECISION, intent(inout)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION, intent(out)    :: BPDFNDVI_VKERNEL(16)

!  local variables

      DOUBLE PRECISION :: DMOD, CF11, CF12, CF21, CF22, CN1, CN2
      DOUBLE PRECISION :: AF11, AF12, AF21, AF22
      DOUBLE PRECISION :: C21, C22, CTTTP, CTTPT, CTTPP
      DOUBLE PRECISION :: CTPPT, CTPPP, CPTPP

!  H-function variables

      DOUBLE PRECISION :: HFUNCTION, NDVI, DEXPNDVI, FACTOR, C
      DOUBLE PRECISION :: ATTEN, FP0, Z, Z2, Z1
      DOUBLE PRECISION :: sgamma, cgamma

!  F-.M. Breon BPDF NDVI model (2009).
!   This is the Vector model with polarization

!  Initialise

      BPDFNDVI_VKERNEL = ZERO

!  Transcription of the RMATR subroutine from Mishchenko/Travis code.
!   CALCULATION OF THE STOKES REFLECTION MATRIX FOR
!   ILLUMINATION FROM ABOVE FOR A SURFACE SEPARATING TWO HALF-SPACES
!   WITH REFRACTIVE INDICES OF THE UPPER AND LOWER HALF-SPACES EQUAL TO
!   CN1 AND CN2, RESPECTIVELY.
!   XI = ABS(COSINE OF THE INCIDENT ZENITH ANGLE)
!   XJ = ABS(COSINE OF THE REFLECTION ZENITH ANGLE).
!   SXI and SXJ are the respective SINES (input)
!   XPHI_REF = REFLECTION AZIMUTH ANGLE

!  PARS(1) = refractive index of water (real)

      CN1 = ONE
      CN2 = PARS(1)

!  Check for limiting cases

      IF(DABS(XI-1D0).LT.1d-9) XI = 0.999999999999d0
      IF(DABS(XJ-1D0).LT.1d-9) XJ = 0.999999999999d0

!   Angle of the surface that generates specular reflection from
!  sun to view directions (theta)
!      alpha = DACOS(HALF*(mus+muv)/dcos(gamma))

      Z = XI * XJ - SXI * SXJ * CPHI
      IF ( Z .GT. ONE) Z = ONE
      Z1 = DACOS(Z)
      Z2 = DCOS(Z1*HALF)

!  PARS(2) = NDVI
!  Exponential of the NDVI.
!    Out of range values default to zero

      NDVI = PARS(2)
      IF ( NDVI .GT. ONE .or. NDVI .lt. MINUS_ONE ) THEN
        NDVI = ZERO
      ENDIF
      DEXPNDVI = DEXP ( - NDVI )
      Fp0 = 0.25d0 * DEXPNDVI / ( xi + xj )

! attenuation factor

      cgamma = Z2
      sgamma = dsqrt ( one - cgamma * cgamma )
      atten  = dexp ( - sgamma / cgamma )

!  PARS(3) = Scaling Factor

      C = PARS(3) 

!  Final H-function

      HFUNCTION = C * Fp0 * atten

!  Call to the vector Fresnel routine

      CALL FRESNEL_VECTOR &
      ( CN1, CN2, XI, SXI, XJ, SXJ, CPHI, SKPHI, & ! Input
        CF11, CF12, CF21, CF22, DMOD )             ! Output

!  Final H-function needs to be normalized

      HFUNCTION = HFUNCTION / DMOD

!  Settting (1,1), (1,2), (2,1) and (2,2) components
!  -------------------------------------------------

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

      BPDFNDVI_VKERNEL(1) = (AF11+AF12+AF21+AF22) * FACTOR
      BPDFNDVI_VKERNEL(2) = (AF11-AF12+AF21-AF22) * FACTOR
      BPDFNDVI_VKERNEL(5) = (AF11-AF22+AF12-AF21) * FACTOR
      BPDFNDVI_VKERNEL(6) = (AF11-AF12-AF21+AF22) * FACTOR

!  Finish if NStokes = 1 or 2 (NSSQ = 1 or 4)

      if ( NSSQ .le. 4 ) return

!  Setting (1,3), (2,3), (3,1), (3,2), (3,3) and (4,4) components
!  --------------------------------------------------------------

      C21 = CF21
      C22 = CF22
      CTTTP=CF11*CF12
      CTTPT=CF11*C21
      CTTPP=CF11*C22
      CTPPT=CF12*C21
      CTPPP=CF12*C22
      CPTPP=CF21*C22

      FACTOR = HFUNCTION

!  New code (Version 2.7, pending verification)
!    We think that sine-terms need to be reversed! Entries 3, 7, 9, 10.

      BPDFNDVI_VKERNEL(3)  =    ( CTTTP+CPTPP) * FACTOR   ! Sine
      BPDFNDVI_VKERNEL(7)  =    ( CTTTP-CPTPP) * FACTOR   ! Sine
      BPDFNDVI_VKERNEL(9)  =    ( CTTPT+CTPPP) * FACTOR   ! Sine
      BPDFNDVI_VKERNEL(10) =    ( CTTPT-CTPPP) * FACTOR   ! Sine
      BPDFNDVI_VKERNEL(11) =    ( CTTPP+CTPPT) * FACTOR   ! Cos
      IF ( NSSQ.eq.16 ) then
         BPDFNDVI_VKERNEL(16) =    ( CTTPP-CTPPT) * FACTOR   ! Cos
      endif

!  Original VLIDORT code (Version 2.6)

!      BPDF2009_VKERNEL(3)  =    (-CTTTP-CPTPP) * FACTOR
!      BPDF2009_VKERNEL(7)  =    (-CTTTP+CPTPP) * FACTOR
!      BPDF2009_VKERNEL(9)  =    (-CTTPT-CTPPP) * FACTOR
!      BPDF2009_VKERNEL(10) =    (-CTTPT+CTPPP) * FACTOR
!      BPDF2009_VKERNEL(11) =    (-CTTPP-CTPPT) * FACTOR   !  Wrongly coded, these are Cosine terms
!      BPDF2009_VKERNEL(16) =    (-CTTPP+CTPPT) * FACTOR   !  Wrongly coded, these are Cosine terms

!  Finish

      RETURN
      END SUBROUTINE BPDFNDVI_VFUNCTION

!

      subroutine FRESNEL_VECTOR &
      ( CN1, CN2, XI, SXI, XJ, SXJ, CKPHI_REF, SKPHI_REF, & ! Input
        CF11, CF12, CF21, CF22, DMOD )                      ! Output

!  module of constamts

      use vlidort_pars, only : zero, one

      implicit none

!  Input variables

      double precision, intent(in)    :: CN1, CN2
      double precision, intent(in)    :: XI, SXI, XJ, SXJ, CKPHI_REF, SKPHI_REF

!  Output

      double precision, intent(inout) :: CF11, CF12, CF21, CF22, DMOD

!  Local variables

      DOUBLE PRECISION :: XPHI_INC, CKPHI_INC, SKPHI_INC
      DOUBLE PRECISION :: VI1, VI2, VI3, VR1, VR2, VR3
      DOUBLE PRECISION :: VP1, VP2, VP3
      DOUBLE PRECISION :: unit1, unit2, unit3, fact1, factor
      DOUBLE PRECISION :: XI1, CXI2, C2, C1, CRPER, CRPAR
      DOUBLE PRECISION :: TI1, TI2, TI3, TR1, TR2, TR3
      DOUBLE PRECISION :: PI1, PII2, PI3, PR1, PR2, PR3
      DOUBLE PRECISION :: PIKR, PRKI, TIKR, TRKI
      DOUBLE PRECISION :: E1, E2, E3, E4

!  For real case, incident azimuth taken to be zero

      XPHI_INC  = ZERO
      CKPHI_INC = ONE
      SKPHI_INC = ZERO

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

!  Set the output

      CF11 =  E1*CRPER+E2*CRPAR
      CF12 = -E3*CRPER+E4*CRPAR
      CF21 = -E4*CRPER+E3*CRPAR
      CF22 =  E2*CRPER+E1*CRPAR

!  Not to forget the normalization

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
        DMOD = ONE
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE FRESNEL_VECTOR

!

subroutine VBRDF_WhiteCap_Reflectance &
    ( WindSpeed, Wavelength, WC_Reflectance, WC_Lambertian )

!  Stand-alone routine for computing the WhiteCap Reflectance
!   Based on 6S code, as updated by A. Sayer (2011)

!   Made compatible with VLIDORT SURFACE LEAVING code
!   renamed for VBRDF code (useful if BRDF and SLEAVE are operating together)
!   R. Spurr, 23 April 2014, 28 April 2014

   implicit none

!  Inputs
!    (Wind speed in [m/s], Wavelength in Microns)

   double precision, intent(in)  :: WindSpeed
   double precision, intent(in)  :: Wavelength

!  output

   double precision, intent(out) :: WC_Reflectance
   double precision, intent(out) :: WC_Lambertian

!  Data
!  ----

!  Single precision

   real :: Effective_WCRef(39)

! effective reflectance of the whitecaps (Koepke, 1984)
! These are the original values - superseded, A Sayer 05 Jul 2011.
!      data Effective_WCRef/ &
!     0.220,0.220,0.220,0.220,0.220,0.220,0.215,0.210,0.200,0.190,&
!     0.175,0.155,0.130,0.080,0.100,0.105,0.100,0.080,0.045,0.055,&
!     0.065,0.060,0.055,0.040,0.000,0.000,0.000,0.000,0.000,0.000,&
!     0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000/

! effective reflectance of the whitecaps (Frouin et al, 1996)
! Assume linear trends between the node points they give
! This is the spectral shape

      data Effective_WCRef/ &
     1.000,1.000,1.000,1.000,0.950,0.900,0.700,0.550,0.500,0.450,&
     0.400,0.350,0.300,0.250,0.200,0.150,0.100,0.050,0.000,0.000,&
     0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,&
     0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000/

!  Local variables
!  ---------------

!  Single precision in the original code

   integer :: iwl, iref
   real    :: Wlb, WLP, Ref(39), wspd, wl, Ref_i, Rwc

!  Initialize

   WC_Reflectance = zero
   WC_Lambertian  = zero

!  Single precision inputs in the original

   wspd = real(WindSpeed)
   wl   = real(Wavelength)

!  Scale data for value of 0.22 in the midvisible.

   DO iref = 1,39
      Ref(iref) = 0.22 * Effective_WCRef(iref)
   ENDDO

!  COMPUTE WHITECAPS REFLECTANCE (LAMBERTIAN)

   Wlb    = 0.0
   IF (wspd .le. 9.25) THEN
       Wlb = 0.01*((3.18e-03)*((wspd-3.7)**3.0))
   ELSE IF (wspd .gt. 9.25) THEN
       Wlb = 0.01*((4.82e-04)*((wspd+1.8)**3.0))
   END IF

! Original whitecap calculation - superseded, A. Sayer 05 Jul 2011.
!      W=2.95e-06*(wspd**3.52)

!  Find data point, Linearly interpolate

   iwl   = 1+int((wl-0.2)/0.1)
   wlp   = 0.5+(iwl-1)*0.1
   Ref_i = Ref(iwl+1) + ( wl-wlp)/0.1*(Ref(iwl)-Ref(iwl+1))
   Rwc   = Wlb*Ref_i

!  Final values

   WC_Lambertian  = real(Wlb,fpk)
   WC_Reflectance = real(Rwc,fpk)

!  Finish

   return
end subroutine VBRDF_WhiteCap_Reflectance

!

subroutine VBRDF_Water_RefracIndex &
     ( Wavelength, Salinity, Refrac_R, Refrac_I )

!  THIS IS FORMERLY CALLED "INDWAT" in 6S

   implicit none
   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Input/Output
!  ------------

   double precision, intent(in)   :: Wavelength
   double precision, intent(in)   :: Salinity
   double precision, intent(out)  :: Refrac_R
   double precision, intent(out)  :: Refrac_I

!  Local (Real-valued quantities)
!  -----

!subroutine indwat(wl,xsal,nr,ni)
!
! input parameters:  wl=wavelength (in micrometers)
!                    xsal=salinity (in ppt), if xsal<0 then 34.3ppt by default
! output parameters: nr=index of refraction of sea water
!                    ni=extinction coefficient of sea water

       real twl(62),tnr(62),tni(62)
       real nr,ni,wl,xwl,yr,yi,nrc,nic,xsal
       integer i

! Indices of refraction for pure water from Hale and Querry, 
! Applied Optics, March 1973, Vol. 12,  No. 3, pp. 555-563
       data twl/&
        0.250,0.275,0.300,0.325,0.345,0.375,0.400,0.425,0.445,0.475,&
        0.500,0.525,0.550,0.575,0.600,0.625,0.650,0.675,0.700,0.725,&
        0.750,0.775,0.800,0.825,0.850,0.875,0.900,0.925,0.950,0.975,&
        1.000,1.200,1.400,1.600,1.800,2.000,2.200,2.400,2.600,2.650,&
        2.700,2.750,2.800,2.850,2.900,2.950,3.000,3.050,3.100,3.150,&
        3.200,3.250,3.300,3.350,3.400,3.450,3.500,3.600,3.700,3.800,&
        3.900,4.000/
        data tnr/&
        1.362,1.354,1.349,1.346,1.343,1.341,1.339,1.338,1.337,1.336,&
        1.335,1.334,1.333,1.333,1.332,1.332,1.331,1.331,1.331,1.330,&
        1.330,1.330,1.329,1.329,1.329,1.328,1.328,1.328,1.327,1.327,&
        1.327,1.324,1.321,1.317,1.312,1.306,1.296,1.279,1.242,1.219,&
        1.188,1.157,1.142,1.149,1.201,1.292,1.371,1.426,1.467,1.483,&
        1.478,1.467,1.450,1.432,1.420,1.410,1.400,1.385,1.374,1.364,&
        1.357,1.351/
        data tni/&
        3.35E-08,2.35E-08,1.60E-08,1.08E-08,6.50E-09,&
        3.50E-09,1.86E-09,1.30E-09,1.02E-09,9.35E-10,&
        1.00E-09,1.32E-09,1.96E-09,3.60E-09,1.09E-08,&
        1.39E-08,1.64E-08,2.23E-08,3.35E-08,9.15E-08,&
        1.56E-07,1.48E-07,1.25E-07,1.82E-07,2.93E-07,&
        3.91E-07,4.86E-07,1.06E-06,2.93E-06,3.48E-06,&
        2.89E-06,9.89E-06,1.38E-04,8.55E-05,1.15E-04,&
        1.10E-03,2.89E-04,9.56E-04,3.17E-03,6.70E-03,&
        1.90E-02,5.90E-02,1.15E-01,1.85E-01,2.68E-01,&
        2.98E-01,2.72E-01,2.40E-01,1.92E-01,1.35E-01,&
        9.24E-02,6.10E-02,3.68E-02,2.61E-02,1.95E-02,&
        1.32E-02,9.40E-03,5.15E-03,3.60E-03,3.40E-03,&
        3.80E-03,4.60E-03/

!  Assign input

      wl   = real(WAVELENGTH)
      xsal = real(SALINITY)
      Refrac_R = zero
      Refrac_I = zero

!  Find wavelength point for interpolation

        i=2
 10     if (wl.lt.twl(i)) goto 20
        if (i.lt.62) then
           i=i+1
           goto 10
         endif

!  Interpolate

 20     xwl=twl(i)-twl(i-1)
        yr=tnr(i)-tnr(i-1)
        yi=tni(i)-tni(i-1)
        nr=tnr(i-1)+(wl-twl(i-1))*yr/xwl
        ni=tni(i-1)+(wl-twl(i-1))*yi/xwl
!
! Correction to be applied to the index of refraction and to the extinction 
! coefficients of the pure water to obtain the ocean water one (see for 
! example Friedman). By default, a typical sea water is assumed 
! (Salinity=34.3ppt, Chlorinity=19ppt) as reported by Sverdrup. 
! In that case there is no correction for the extinction coefficient between 
! 0.25 and 4 microns. For the index of refraction, a correction of +0.006 
! has to be applied (McLellan). For a chlorinity of 19.0ppt the correction 
! is a linear function of the salt concentration. Then, in 6S users are able 
! to enter the salt concentration (in ppt).
! REFERENCES:
! Friedman D., Applied Optics, 1969, Vol.8, No.10, pp.2073-2078.
! McLellan H.J., Elements of physical Oceanography, Pergamon Press, Inc.,
!        New-York, 1965, p 129.
! Sverdrup H.V. et al., The Oceans (Prentice-Hall, Inc., Englewood Cliffs,
!        N.J., 1942, p 173.

        nrc=0.006
        nic=0.000
        nr=nr+nrc*(xsal/34.3)
        ni=ni+nic*(xsal/34.3)

!  Assign output

    REFRAC_R = real(nr,fpk)
    REFRAC_I = real(ni,fpk)

    return
end subroutine VBRDF_Water_RefracIndex


SUBROUTINE VBRDF_Generalized_Glint &
         ( DO_ISOTROPIC, DO_SHADOW, DO_COEFFS,    &
           REFRAC_R, REFRAC_I, WINDSPEED,         &
           PHI_W, CPHI_W, SPHI_W,                 &
           XJ, SXJ, XI, SXI, PHI, CPHI, SPHI,     &
           SUNGLINT_COEFFS, SUNGLINT_REFLEC )

      implicit none

!  Subroutine Input arguments
!  --------------------------

!  Flag for using Isotropic Facet distribution

      LOGICAL  , intent(in)    :: DO_ISOTROPIC

!  Flag for including Shadow effect

      LOGICAL  , intent(in)    :: DO_SHADOW

!  Flag for Calculating Cox-Munk Coefficients
!     Only needs to be done once, so intent(inout)

      LOGICAL  , intent(inout) :: DO_COEFFS

!  Real and imaginary parts of refractive index

      double precision, intent(in)    :: REFRAC_R
      double precision, intent(in)    :: REFRAC_I

!  Windspeed m/s

      double precision, intent(in)    :: WINDSPEED

!  Azimuth between Sun and Wind directions. angle in Radians + Cosine/sine

      double precision, intent(in)    :: PHI_W, CPHI_W, SPHI_W

!  Incident and reflected ddirections: sines/cosines. Relative azimuth (angle in radians)

      double precision, intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SPHI

!  Subroutine output arguments
!  ---------------------------

!   Glitter reflectance

      double precision, intent(out)   :: SUNGLINT_REFLEC

!  Cox-Munk Coefficients. Intent(inout).

      double precision, intent(inout) :: SUNGLINT_COEFFS(7)

!  Local arguments
!  ---------------

      double precision, PARAMETER :: CRITEXP = 88.0D0
      double precision, PARAMETER :: six = two * three, twentyfour = six * four

!  Local variables

      double precision  :: B, ZX, ZY, Z, Z1, Z2, XMP
      double precision  :: TILT, TANTILT, TANTILT_SQ, COSTILT
      double precision  :: ARGUMENT, PROB, FAC2, COEFF, VAR, WSigC, WSigU
      double precision  :: XE, XN, XE_sq, XN_sq, XE_sq_1, XN_sq_1
      double precision  :: XPHI, CKPHI, SKPHI, XPHI_W, CKPHI_W, SKPHI_W
      double precision  :: S1, S2, S3, XXI, XXJ, T1, T2, DCOT
      double precision  :: SHADOWI, SHADOWR, SHADOW

!  Initialise output

      SUNGLINT_REFLEC = ZERO

!  Compute coefficients, according to 6S formulation

      IF ( DO_COEFFS ) THEN
         SUNGLINT_COEFFS = zero
         IF ( DO_ISOTROPIC ) THEN
            SUNGLINT_COEFFS(1) = 0.003d0 + 0.00512d0 * WINDSPEED
         ELSE
            SUNGLINT_COEFFS(1) = 0.003d0 + 0.00192d0 * WINDSPEED ! sigmaC
            SUNGLINT_COEFFS(2) =           0.00316d0 * WINDSPEED ! sigmaU
            SUNGLINT_COEFFS(3) = 0.010d0 - 0.00860d0 * WINDSPEED ! C21
            SUNGLINT_COEFFS(4) = 0.040d0 - 0.03300d0 * WINDSPEED ! C03
            SUNGLINT_COEFFS(5) = 0.400d0                         ! C40
            SUNGLINT_COEFFS(6) = 0.230d0                         ! C04
            SUNGLINT_COEFFS(7) = 0.120d0                         ! C22
         ENDIF
         DO_COEFFS = .false.
      ENDIF

!  Local angles

      XPHI   = PIE - PHI       ! Not used
!     CKPHI  = - CPHI          ! Original, not correct.

      CKPHI  = + CPHI
      SKPHI  = + SPHI

      XPHI_W  = PHI_W
      CKPHI_W = CPHI_W
      SKPHI_W = SPHI_W

!  Tilt angle

      B  = ONE / ( XI + XJ )
      ZX = - SXI * SKPHI * B
      ZY = ( SXJ + SXI * CKPHI ) * B
      TANTILT_SQ = ZX * ZX + ZY * ZY
      TANTILT    = SQRT ( TANTILT_SQ )
      TILT       = ATAN(TANTILT)
      COSTILT    = COS(TILT)

!  Scatter angle

      Z = XI * XJ + SXI * SXJ * CKPHI
      IF ( Z .GT. ONE) Z = ONE
      Z1 = ACOS(Z)
      Z2 = COS(Z1*HALF)

!  Fresnel
!  -------

       CALL VBRDF_Fresnel_Complex ( REFRAC_R, REFRAC_I, Z2, XMP )

!  Anisotropic
!  -----------

      IF ( .not. DO_ISOTROPIC ) THEN

!  Variance

         WSigC = ONE / Sqrt(SUNGLINT_COEFFS(1))
         WSigU = ONE / Sqrt(SUNGLINT_COEFFS(2))
         VAR   = WSigC * WSigU * HALF 

!  angles

         XE = (  CKPHI_W * ZX + SKPHI_W * ZY ) * WSigC ; XE_sq = XE * XE ; XE_sq_1 = xe_sq - one
         XN = ( -SKPHI_W * ZX + CKPHI_W * ZY ) * WSigU ; XN_sq = XN * XN ; XN_sq_1 = xn_sq - one

!  GC Coefficient

         Coeff = ONE - SUNGLINT_COEFFS(3) *      XE_sq_1      * XN * half &
                     - SUNGLINT_COEFFS(4) * ( XN_sq - three ) * XN / six  &
                     + SUNGLINT_COEFFS(5) * ( XE_sq * XE_sq - six * XE_sq + three ) / twentyfour &
                     + SUNGLINT_COEFFS(6) * ( XN_sq * XN_sq - six * XN_sq + three ) / twentyfour &
                     + SUNGLINT_COEFFS(7) * XE_sq_1 * XN_sq_1 / four

!  Probability and finish

         ARGUMENT = ( XE_sq  + XN_sq ) * HALF
         IF ( ARGUMENT .LT. CRITEXP ) THEN
            PROB = COEFF * EXP ( - ARGUMENT ) * VAR
            FAC2 = QUARTER / XI / XJ / ( COSTILT ** FOUR )
            SUNGLINT_REFLEC = XMP * PROB * FAC2
         ENDIF

      ENDIF

!  Isotropic
!  ---------

      IF ( DO_ISOTROPIC ) THEN

!  Compute Probability and finish

         COEFF = ONE
         VAR   = SUNGLINT_COEFFS(1)
         ARGUMENT = TANTILT_SQ / VAR
         IF ( ARGUMENT .LT. CRITEXP ) THEN
            PROB = COEFF * EXP ( - ARGUMENT ) / VAR
            FAC2 = QUARTER / XI / XJ / ( COSTILT ** FOUR )
            SUNGLINT_REFLEC = XMP * PROB * FAC2
         ENDIF

      ENDIF

!  No Shadow code if not flagged

      IF ( .not. DO_SHADOW  ) RETURN

!  Shadow code

      S1 = SQRT ( VAR / PIE )
      S3 = ONE / ( SQRT(VAR) )
      S2 = S3 * S3

      XXI  = XI*XI
      DCOT = XI / SQRT ( ONE - XXI )
      T1   = EXP ( - DCOT * DCOT * S2 )
      T2   = DERFC_E ( DCOT * S3 )
!      T2   = DCOT * S3 ; CALL HOMEGROWN_ERRFUNC ( T2 )     !  Error function usage in SLEAVE code
      SHADOWI = HALF * ( S1 * T1 / DCOT - T2 )

      XXJ  = XJ*XJ
      DCOT = XJ / SQRT ( ONE - XXJ )
      T1   = EXP ( - DCOT * DCOT * S2 )
      T2   = DERFC_E ( DCOT * S3 )
!      T2   = DCOT * S3 ; CALL HOMEGROWN_ERRFUNC ( T2 )     !  Error function usage in SLEAVE code
      SHADOWR = HALF * ( S1 * T1 / DCOT - T2 )

      SHADOW = ONE / ( ONE + SHADOWI + SHADOWR )
      SUNGLINT_REFLEC = SUNGLINT_REFLEC * SHADOW

!     Finish

      RETURN
      END SUBROUTINE VBRDF_Generalized_Glint

!  End module

      END MODULE vbrdf_sup_kernels_m
