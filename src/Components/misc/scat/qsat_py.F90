subroutine qsat ( qs, T, p, nobs )

  use MAPL_SatVaporMod, only: MAPL_EQsat

  integer, intent(in) :: nobs
  real,    intent(in) :: p(nobs) ! pressure in Pa
  real,    intent(in) :: T(nobs) ! dry-bulb temperature in K

  real,    intent(out) :: qs(nobs) ! saturation specific humidity

!               ---

  integer i

  do i = 1, nobs

     qs(i) = MAPL_EQsat(TL=t(i),PL=p(i),OverIce=.True.)

  end do

end subroutine qsat

subroutine qsatnoice ( qs, T, p, nobs )

  use MAPL_SatVaporMod, only: MAPL_EQsat

  integer, intent(in) :: nobs
  real,    intent(in) :: p(nobs) ! pressure in Pa
  real,    intent(in) :: T(nobs) ! dry-bulb temperature in K

  real,    intent(out) :: qs(nobs) ! saturation specific humidity

!               ---

  integer i

  do i = 1, nobs

     qs(i) = MAPL_EQsat(TL=t(i),PL=p(i),OverIce=.False.)

  end do

end subroutine qsatnoice




