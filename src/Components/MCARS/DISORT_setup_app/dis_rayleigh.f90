!	Created by Kerry Meyer on 2/8/08.
!
!	Parameters, variables, and equations for Rayleigh optical depth calculations.
!	Constants and formulas from:
!		Bodhaine et al., 1999: On Rayleigh Optical Depth Calculations. J. Atmos. and Ocean Tech., 16. 1854-1861.
!
!	Parameters:
!		CO2_ppm		= concentration of CO2 in ppm
!		F_Ar		= depolarization factor for argon
!		F_CO2		= depolarization factor for CO2
!		ndens_da	= number density of dry air
!		av_num		= Avogadro's number
!		g0			= gravitational constant
!
!	Variables:
!		atm_press	= atmospheric pressure at each level from atmospheric profile
!		F_O2		= depolarization factor for O2 at each channel
!		F_N2		= depolarization factor for N2 at each channel
!		F_air		= depolarization factor for dry air at each channel
!		refin_air	= refractive index of dry air at each channel
!		ray_cross	= Rayleigh scattering cross-section at each channel
!		mw_air		= molecular weight of dry air
!		g1,g2		= gravity, corrected for altitude

module dis_rayleigh

use dis_global_pars

implicit none
public

real, parameter :: CO2_ppm = 360, &
				   F_Ar = 1.0, &
				   F_CO2 = 1.15, &
				   ndens_da = 2.546899E19, &
				   av_num = 6.0221367E23, &
				   g0 = 9.832
				   
real, dimension(maxulv) :: atm_press, atm_height

!**************************************************************************************

contains

subroutine readprof (atmprof,nlyr,atm_press,atm_height)

use dis_global_pars

implicit none

!	Dummy variables
character(len=100), intent(in) :: atmprof

real, dimension(maxulv), intent(out) :: atm_press, atm_height

integer, intent(in) :: nlyr

!	Local variables
integer :: i
logical :: debug
debug = .false.
	open (50,file=atmprof)

	!print*,'nlyr: ',nlyr
	do i = 1, nlyr+1

		read (50,*) atm_height(i), atm_press(i)
        if (debug) print*,i,atm_height(i), atm_press(i)
	end do
	
	close (50)

end subroutine readprof

!**************************************************************************************

subroutine rayleigh_tau_calc (nlyr,ich,atm_press,atm_height,wvlo,wvhi,rayleigh_tau)

use dis_global_pars

implicit none

!	Dummy variables
real, dimension(:), intent(inout) :: rayleigh_tau
real, dimension(:), intent(in) :: atm_press, atm_height
real, intent(in) :: wvlo, wvhi

integer, intent(in) :: nlyr

!	Local variables
real :: wvln
real :: ray_cross_num, ray_cross_den, ray_cross
real :: F_O2, F_N2, F_air, refin_air
real :: mw_air, g1, g2

integer :: i, ich
logical :: debug
debug = .false.
!	Calculate molecular weight of air (corresponding to CO2 concentration)
	mw_air = 15.0556*(CO2_ppm/1000000.) + 28.9595

!	Calculate Rayleigh scattering cross-section
	wvln = (1. / ((wvlo + wvhi) / 2.)) * 10000.
	F_N2 = 1.034 + 3.17E-4*(1./(wvln**2))
	F_O2 = 1.096 + 1.385E-3*(1./(wvln**2)) + 1.448E-4*(1./(wvln**4))
	F_air = (78.084*F_N2 + 20.946*F_O2 + 0.934*F_Ar + CO2_ppm*F_CO2/10000.) / &
				(78.084 + 20.946 + 0.934 + CO2_ppm/10000.)
	refin_air = ((8060.51 + (2480990./(132.274-wvln**-2)) + (17455.7/(39.32957-wvln**-2)))*&
					(1. + 0.54*((CO2_ppm/1000000.)-0.0003)) / 10**8) + 1.
	ray_cross_num = 24. * (acos(-1.)**3) * ((refin_air**2 - 1.)**2)
	ray_cross_den = ((wvln*0.0001)**4) * ((refin_air**2 + 2.)**2) * ndens_da
	ray_cross_den = ray_cross_den * ndens_da
	ray_cross = (ray_cross_num / ray_cross_den) * F_air
	
	do i = 1, nlyr
		g1 = g0 * (1.-(2.*atm_height(i+1)/6.37E6))
		g2 = g0 * (1.-(2.*atm_height(i)/6.37E6))
		rayleigh_tau(i) = ((ray_cross*av_num) / mw_air) * (atm_press(i+1)/g1 - atm_press(i)/g2) * 10.
		if (debug) print*,i,rayleigh_tau(i)
	end do

end subroutine rayleigh_tau_calc

end module dis_rayleigh