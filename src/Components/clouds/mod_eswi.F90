! file: mod_eswi.f90
! Saturation vapor pressure over water and ice
module mod_eswi

  implicit none

  private
  public esw_gcet, esi_gcet

  ! need this to back out esw, esi
  real*4, parameter :: &
    Rd = 287.04, Rv = 461.50, &
    epsilon = Rd / Rv, &
    coeff = 3.799052e2 / epsilon

contains 

  ! Saturation vapor pressure over water
  !   esw = esw_gcet(tK), tK in degK, esw in Pa
  ! GCE version based on Teten's formula
  ! (see Tao et al., MWR 117, 231--235, 1989)
  ! This one gives exactly equal esw and esi at T = 273.16

  function esw_gcet(tK) result(esw)

    implicit none
    real*4, intent(in) :: tK(:)
    real*4 :: esw(size(tK))
  
    ! GCE calculation of qsw:
    ! RP0 = 3.799052e3 / P, with P in tenths of Pa
    ! qsw = rp0 * exp(17.26939 * (tK - 273.16) / (tK - 35.86)) in [g/g]
    ! Note that P in denom implies GCE is using qsw = epsilon * esw / P,
    ! i.e. saturation vapor "content" in language of Norris et al. (2008)

    ! saturation vapor pressure in Pa over water
    esw = coeff * exp(17.26939 * (tK - 273.16) / (tK - 35.86))
    ! Note: This is consistent with the GCE qsw = epsilon * esw / P above,
    !       if P is in Pa

  end function esw_gcet

  ! Saturation vapor pressure over ice
  !   esi = esi_gcet(tK), tK in degK, esi in Pa
  ! GCE version based on Teten's formula
  ! (see Tao et al., MWR 117, 231--235, 1989)
  ! This one gives exactly equal esw and esi at T = 273.16

  function esi_gcet(tK) result(esi)

    implicit none
    real*4, intent(in) :: tK(:)
    real*4 :: esi(size(tK))
  
    ! GCE calculation of qsi:
    ! RP0 = 3.799052e3 / P, with P in tenths of Pa
    ! qsi = rp0 * exp(21.87456 * (tK - 273.16) / (tK - 7.66)) in [g/g]
    ! Note that P in denom implies GCE is using qsi = epsilon * esi / P,
    ! i.e. saturation vapor "content" in language of Norris et al. (2008)

    ! saturation vapor pressure in Pa over ice
    esi = coeff * exp(21.87456 * (tK - 273.16) / (tK - 7.66))
    ! Note: This is consistent with the GCE qsi = epsilon * esi / P above,
    !       if P is in Pa

  end function esi_gcet

end module mod_eswi
