MODULE RAYLEIGH
  IMPLICIT NONE

  real(kind=8), parameter ::   NA = 6.0221367D+23
  real(kind=8), parameter ::   MAIR_NOCO2=28.9595D0


interface DELTAU_RAYLEIGH
 module procedure DELTAU_RAYLEIGH_1D
 module procedure DELTAU_RAYLEIGH_2D
end interface DELTAU_RAYLEIGH

  CONTAINS

SUBROUTINE RAYLEIGH_FUNCTION &
  ( FORWARD_MAXLAMBDAS, CO2_PPMV_MIXRATIO, &
  FORWARD_NLAMBDAS,   FORWARD_LAMBDAS, &
  RAYLEIGH_XSEC, RAYLEIGH_DEPOL, NOUTM1 )
  
  !  Rayleigh cross sections and depolarization ratios
  !     Bodhaine et. al. (1999) formulae
  !     Module is stand-alone.

  !  Subroutine coded by R. Spurr, RT Solutions, Cambridge, MA.

  !  Corrected CO2 mixing ratio units for King Factor calculation, 3/2021, D. Haffner
  
  IMPLICIT NONE
  
  !  Input arguments
  !  ---------------
  !  wavelength
  
  INTEGER          FORWARD_MAXLAMBDAS, FORWARD_NLAMBDAS
  !     DOUBLE PRECISION FORWARD_LAMBDAS ( FORWARD_MAXLAMBDAS )
  ! use real to be consistent with the input from main prog
  real             FORWARD_LAMBDAS ( FORWARD_MAXLAMBDAS )
  
  !  CO2 mixing ratio
  
  DOUBLE PRECISION CO2_PPMV_MIXRATIO
  
  !  Output arguments
  !  ----------------
  
  !  cross-sections and depolarization output
  
  DOUBLE PRECISION RAYLEIGH_XSEC  ( FORWARD_MAXLAMBDAS )
  DOUBLE PRECISION RAYLEIGH_DEPOL ( FORWARD_MAXLAMBDAS )
  DOUBLE PRECISION NOUTM1 ( FORWARD_MAXLAMBDAS )
  
  !  Local variables
  !  ---------------
  
  INTEGER          W
  DOUBLE PRECISION MASS_DRYAIR
  DOUBLE PRECISION NMOL, PI, CONS
  DOUBLE PRECISION MO2,MN2,MARG,MCO2,MAIR
  DOUBLE PRECISION FO2,FN2,FARG,FCO2,FAIR
  DOUBLE PRECISION LAMBDA_C,LAMBDA_M,LPM2,LP2
  DOUBLE PRECISION N300M1,NCO2M1,NCO2
  DOUBLE PRECISION NCO2SQ, NSQM1,NSQP2,TERM
  DOUBLE PRECISION S0_A, S0_B
  DOUBLE PRECISION S1_A, S1_B, S1_C, S1_D, S1_E
  DOUBLE PRECISION S2_A
  DOUBLE PRECISION S3_A, S3_B, S3_C, S3_D, S3_E

  !  data statements and parameters
  !  ------------------------------
  
  DATA MO2  / 20.946D0 /
  DATA MN2  / 78.084D0 /
  DATA MARG / 0.934D0 /
  
  PARAMETER        ( S0_A = 15.0556D0 )
  PARAMETER        ( S0_B = 28.9595D0 )
  
  PARAMETER        ( S1_A = 8060.51D0 )
  PARAMETER        ( S1_B = 2.48099D+06 )
  PARAMETER        ( S1_C = 132.274D0 )
  PARAMETER        ( S1_D = 1.74557D+04 )
  PARAMETER        ( S1_E = 39.32957D0 )
  
  PARAMETER        ( S2_A = 0.54D0 )
  
  PARAMETER        ( S3_A = 1.034D0 )
  PARAMETER        ( S3_B = 3.17D-04 )
  PARAMETER        ( S3_C = 1.096D0 )
  PARAMETER        ( S3_D = 1.385D-03 )
  PARAMETER        ( S3_E = 1.448D-04 )
  
  !  Start of code
  !  -------------
  
  !  constants
  
  NMOL = 2.546899D19
  PI   = DATAN(1.0D0)*4.0D0
  CONS = 24.0D0 * PI * PI * PI
  
  !  convert co2
  
  MCO2 = 1.0D-06 * CO2_PPMV_MIXRATIO
  
  !  mass of dry air: Eq.(17) of BWDS
  
  MASS_DRYAIR = S0_A * MCO2 + S0_B
  
  !  start loop
  
  DO W = 1, FORWARD_NLAMBDAS
  
  !  wavelength in micrometers
  
  LAMBDA_M = 1.0D-03 * FORWARD_LAMBDAS(W)
  LAMBDA_C = 1.0D-07 * FORWARD_LAMBDAS(W)
  LPM2     = 1.0D0 / LAMBDA_M / LAMBDA_M
  
  !  step 1: Eq.(18) of BWDS

  N300M1 = S1_A + ( S1_B / ( S1_C - LPM2 ) ) + &
  ( S1_D / ( S1_E - LPM2 ) )
  N300M1 = N300M1 * 1.0D-08
  
  !  step 2: Eq.(19) of BWDS
  
  NCO2M1 = N300M1 * ( 1.0D0 + S2_A * ( MCO2  - 0.0003D0 ) )
  NCO2   = NCO2M1 + 1
  NCO2SQ = NCO2 * NCO2

  NOUTM1(W) = NCO2M1
  
  !  step 3: Eqs. (5&6) of BWDS (Bates results)
  
  FN2  = S3_A + S3_B * LPM2
  FO2  = S3_C + S3_D * LPM2 + S3_E * LPM2 * LPM2
  
  !  step 4: Eq.(23) of BWDS
  !     ---> King factor and depolarization ratio
  !
  ! convert CO2 mixing ratio to percent for King factor calc.
  
  FARG = 1.0D0
  FCO2 = 1.15D0
  MAIR = MN2 + MO2 + MARG + 100.D0*MCO2
  FAIR = MN2*FN2 + MO2*FO2 + MARG*FARG + 100.D0*MCO2*FCO2
  FAIR = FAIR / MAIR
  RAYLEIGH_DEPOL(W) = 6.0D0*(FAIR-1.0D0)/(3.0D0+7.0D0*FAIR)
  
  !  step 5: Eq.(22) of BWDS
  !     ---> Cross section
  
  LP2  = LAMBDA_C * LAMBDA_C
  NSQM1 = NCO2SQ - 1.0D0
  NSQP2 = NCO2SQ + 2.0D0
  TERM = NSQM1 / LP2 / NMOL / NSQP2
  RAYLEIGH_XSEC(W) =  CONS * TERM * TERM * FAIR
  
!  end loop
ENDDO

!  finish
!  ------

RETURN
END SUBROUTINE RAYLEIGH_FUNCTION

FUNCTION TAU_RAYLEIGH(RAYLEIGH_XSEC, PRESSURE, GEE, CO2_PPMV_MIXRATIO)
  IMPLICIT NONE
  real(kind=8), intent(in) :: RAYLEIGH_XSEC, CO2_PPMV_MIXRATIO, PRESSURE
  REAL(KIND=4), intent(in) :: GEE
  real(kind=8) :: TAU_RAYLEIGH
  real(kind=8) :: MAIR, MCO2
  MCO2 = 1.0D-6 * CO2_PPMV_MIXRATIO
  MAIR = 15.0556D0 * MCO2 + MAIR_NOCO2 
  TAU_RAYLEIGH = RAYLEIGH_XSEC * PRESSURE * NA / MAIR / GEE
END FUNCTION TAU_RAYLEIGH


FUNCTION DELTAU_RAYLEIGH_1D(RAYLEIGH_XSEC, PRESSURE, GEE, CO2_PPMV_MIXRATIO)
  IMPLICIT NONE
  real(kind=8), intent(in) :: RAYLEIGH_XSEC
  REAL(KIND=8), intent(in) :: PRESSURE(:)
  REAL (KIND=4), intent(in) :: GEE
  REAL(KIND=8), intent(in) ::  CO2_PPMV_MIXRATIO
  real(kind=8)  :: DELTAU_RAYLEIGH_1D(SIZE(PRESSURE))
  integer(kind=4) :: NLAYERS, NLAMBDAS
  real(kind=8) :: MAIR, MCO2

  integer(kind=4) :: I, J
  MCO2 = 1.0D-6 * CO2_PPMV_MIXRATIO
  MAIR = 15.0556D0 * MCO2 + MAIR_NOCO2 
  NLAYERS  = SIZE(PRESSURE)
  DO J=1, NLAYERS
     DELTAU_RAYLEIGH_1D(J) = TAU_RAYLEIGH(RAYLEIGH_XSEC, PRESSURE(J), GEE, CO2_PPMV_MIXRATIO)
  ENDDO
END FUNCTION DELTAU_RAYLEIGH_1D


FUNCTION DELTAU_RAYLEIGH_2D(RAYLEIGH_XSEC, PRESSURE, GEE, CO2_PPMV_MIXRATIO)
  IMPLICIT NONE
  real(kind=8), intent(in) :: RAYLEIGH_XSEC(:)
  REAL(KIND=8), intent(in) :: PRESSURE(:)
  REAL (KIND=4), intent(in) :: GEE
  REAL(KIND=8), intent(in) ::  CO2_PPMV_MIXRATIO
  real(kind=8)  :: DELTAU_RAYLEIGH_2D(SIZE(RAYLEIGH_XSEC), SIZE(PRESSURE))
  integer(kind=4) :: NLAYERS, NLAMBDAS
  real(kind=8) :: MAIR, MCO2
  integer(kind=4) :: I, J
  MCO2 = 1.0D-6 * CO2_PPMV_MIXRATIO
  MAIR = 15.0556D0 * MCO2 + MAIR_NOCO2 
  NLAMBDAS = SIZE(RAYLEIGH_XSEC)
  NLAYERS  = SIZE(PRESSURE)
  DO I=1, NLAMBDAS
    DO J=1, NLAYERS
     DELTAU_RAYLEIGH_2D(I,J) = TAU_RAYLEIGH(RAYLEIGH_XSEC(I), PRESSURE(J), GEE, CO2_PPMV_MIXRATIO)
    ENDDO
  ENDDO
END FUNCTION DELTAU_RAYLEIGH_2D

END MODULE RAYLEIGH
