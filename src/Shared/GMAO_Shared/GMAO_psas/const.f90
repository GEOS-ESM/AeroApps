module const

  implicit none

	! geo-physical constants

  real,parameter :: Radius_of_Earth		= 6.371*1.e6	! m
  real,parameter :: Acc_dueto_Gravity		= 9.8100	! m/s/s
  real,parameter :: Earths_angular_velocity	= 7.292e-5	! /s

  real,parameter :: R_EARTH	= Radius_of_Earth
  real,parameter :: G_EARTH	= Acc_dueto_Gravity
  real,parameter :: OMEGA	= Earths_angular_velocity

				! = 2*sin(ONEPI/4) = sqrt(2.)*OMEGA
  real,parameter :: Coriolis_45			= 1.4142136 * OMEGA
  real,parameter :: Rate_of_Coriolis_45		= Coriolis_45 / R_EARTH

  real,parameter :: BETA	= Rate_of_Coriolis_45

				! = 4*atan(1.)
  real,parameter :: ONEPI=3.1415927

  real,parameter :: T_ref	= 273.15	! Kelvin
  real,parameter :: P_ref	= 1000.		! mbar
  real,parameter :: R_dry	= 287.00    ! Gas constant of dry air

  real,parameter :: RHOBAR = 1.24 ! mean sea-level density, kg/m^3
end module const
!.
