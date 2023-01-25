MODULE GRAVITY

  USE, INTRINSIC :: ISO_FORTRAN_ENV , ONLY : REAL64, REAL32
  IMPLICIT NONE
  INTEGER, PARAMETER :: sp = REAL32
  REAL(KIND=sp), PARAMETER :: DTOR=ACOS(-1.0)/180.0

CONTAINS
  
  ! Compute gravitational acceleration at surface of the Earth as a function of latitude
  ! Eq 11 of Bodhaine (original source given in this paper)
  FUNCTION GNOT(LAT)
    IMPLICIT NONE
    REAL(KIND=sp), INTENT(IN) :: LAT
    REAL(KIND=sp) :: GNOT, COS2LAT
    COS2LAT=COS(2.0*LAT*DTOR)
    ! EQ (11)
    GNOT = 980.6160E0 * (1 - 0.0026373E0 * COS2LAT + 0.0000059E0 * COS2LAT*COS2LAT)     
  END FUNCTION GNOT


  ! Mass-weighted altitude equation from Bodhaine et al., 1999.
  ! The value for g needs to be representative of the mass weighted column of air molecules above the site, and should be calculated from Eqs. (10)–(11), modified by using a value of z_c determined from the U.S. Standard Atmosphere, as provided by List (1968). To determine z_c we used List’s (1968, p. 267) table of the density of air as a function of altitude and calculated a mass-weighted mean above each altitude value, using an average altitude and average density for each layer listed in the table. Next a least squares straight line was passed through the resulting z_c values up to 10500 m, giving the following equation:
  !   z_c = 0.73737*z + 5517.56 ! EQ (26)
  ! where z is the altitude of the observing site and zc is the effective mass-weighted altitude of the column. For example, an altitude of z 5 0 m yields an effective mass weighted column altitude of zc = 5517.56 m to use in the calculation of g. The resulting values or t R should be considered the best currently available values for the most accurate estimates of optical depth.
  FUNCTION ZC(Z)
    IMPLICIT NONE
    REAL(KIND=sp), INTENT(IN) :: Z
    REAL(KIND=sp) :: ZC
    ! EQ (26)
    ZC = 0.73737 * Z + 5517.56     
  END FUNCTION ZC

  ! Eq. (10) of Bodhaine et al., 1999 following the formula of List, 1968 for acceleration of gravity, g as a function latitude, height above sea level, z in meters and sea level acceleration of gravity, g0.
  FUNCTION GEE(GNOT, Z, LAT)
    IMPLICIT NONE
    REAL(KIND=sp), INTENT(IN) :: GNOT, Z, LAT
    REAL(KIND=sp) :: GEE, COS2LAT
    COS2LAT=COS(2.0*LAT*DTOR)
    ! EQ (10)
    GEE = GNOT - (3.085462D-4 + 2.27D-7 * COS2LAT) * Z &
         + (7.254D-11 + 1.0D-13 * COS2LAT) * Z**2 &
         - (1.517D-17 + 6.0D-20 * COS2LAT) * Z**3
  END FUNCTION GEE

END MODULE GRAVITY





