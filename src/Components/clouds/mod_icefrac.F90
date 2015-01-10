! file: mod_icefrac.f90
! ice/liq fraction as function of temperature

module mod_icefrac

  implicit none
  private
  public icefrac

  real*4, parameter :: T0 = 273.16, T00 = T0 - 35.

contains

  subroutine icefrac (nk, T, fice, fliq)
  
    integer, intent(in) :: nk
    real*4, intent(in)  :: T(nk)
    real*4, intent(out) :: fice(nk), fliq(nk)

    fice = (T0 - T) / (T0 - T00)
    where (fice < 0.)
      fice = 0.
    elsewhere (fice > 1.)
      fice = 1.
    endwhere
    fliq = 1. - fice

  end subroutine icefrac

end module mod_icefrac
