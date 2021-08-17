module corrK_profile

! -=- -=- -=- -=- -=- -=- -=- -=- -=- -=- -=- -=- -=- -=- -=- -=- -=- 
! program to read climatological atmospheric profiles found in this
! directory, in preparation for building corresponding optical depth 
! profiles using Dave Kratz's correlated-k routines (e.g. modis01.f).  
! These optical depth profiles will be used as input to the discrete 
! ordinates radiative transfer code.
!
! input atmospheric profile is assumed to be from a radiosonde and
! already interpolated to the correct 35 levels
!
! NOTE: if input water vapor mixing ratio is zero at any level, a bunch
!       of NaNs will appear in the correlated k optical thicknesses
!
! Original routine by Chris Yost
! Edited by Bryan Baum September/October 2006
! adapted for HPC by Gala Wind, April 2012
! 
! Routine expect input with standard atmospheric profiles of form:
! Height     Pressure     Temp     H20 mixing ratio     O3 Density
!   km         mbar        K          g kg-1              ppmv 
!
! CALLS:  modis1.f90, modis2.f90, modis3.f90, modis4.f90, modis5a.f90, 
!         modis5b.f90, modis6.f90, modis7.f90, modis17.f90, modis18.f90,
!         modis19.f90, modis20.f90, modis22.f90, modis23.f90, modis26.f90,
!         modis27.f90, modis28.f90, modis29.f90, modis31.f90,
!         modis32.f90, modis33.f90, modis34.f90, modis35.f90, modis36.f90
!
! MODIFICATIONS: 
!    06/3/2005 - Allowed the string 'profile' to be as long as 50
!                characters instead of 30 (Yost)
!    9/7/2006 - corrected conversion of water vapor mixing ratio to g cm^-2 km^-1 (Baum)
!    9/7/2006 - worked with 24 Kratz routines, modified to f90 (baum)
!    4/10/2012 - adapted for HPC execution and converted all to proper F90
! -=- -=- -=- -=- -=- -=- -=- -=- -=- -=- -=- -=- -=- -=- -=- -=- -=- 

	use CKD_common	
	
	implicit none

	private

	type(ckd_type) :: modis_channel(25)

	real, parameter :: u_co2 = 3.80e-04


	public :: corrK_driver, init_corrK, cleanup_corrK, modis_channel


contains

	subroutine corrK_driver(x, y)

		use disort_ancillary
		
		integer, intent(in) :: x, y

		integer ::i
		
		use_levels = granule_info(x,y)%surface_level
		use_layers = use_levels-1
		

!		do i=1, use_levels
!			print*, granule_info(x,y)%heights(i), granule_info(x,y)%ozone(i)
!		
!		end do

		
	    call run_corrK(granule_info(x,y)%heights, granule_info(x,y)%pressure, granule_info(x,y)%temperature, &
						granule_info(x,y)%mix_ratio, granule_info(x,y)%ozone, u_co2)

		
	end subroutine corrK_driver



	subroutine run_corrK (zlev, p0, temp, umix_h2o, uppmv_o3, u_co2)

		use CKD_common
		use modis1
		use modis2
		use modis3
		use modis4
		use modis19
		use modis5b
		use modis5a 
		use modis26
		use modis6
		use modis7
		use modis8
		use modis9
		use modis17
		use modis18
		use modis20
		use modis22
		use modis27
		use modis28
		use modis29
		use modis31
		use modis32
		use modis33
		use modis34
		use modis35
		use modis36
	
	
		implicit none

	
		real, dimension(:), intent(in) :: zlev, p0, temp, umix_h2o, uppmv_o3
		real, intent(in) :: u_co2
		
		real, parameter :: e622 = 0.622
		real, parameter :: Rv = 461.
		real, parameter :: ATM = 1013.25
		
		integer, parameter :: nlayers = 26 !34
		integer, parameter :: nlevels = 27 !35
		
		real, dimension(nlevels) :: patm, u_h2o, u_o3
		real, dimension(nlayers) :: dzz
		real :: ppa, umix_o3

		integer :: m 

		do m=1,use_levels

         patm(m)  = p0(m)/ATM  ! convert from mb to atm
		 ppa = p0(m) * 100.0 ! convert from mb to Pa

		! 1.0e05 is number of cm in km
		! 1.0e06 is conversion from m^-3 to cm^-3
		! 28.9644 is molecular weight of dry air
		! 8.314e+07 is (R*) universal gas constant in cgs units
		! R(h2o) = gas constant = R* = 8.314e+07 bunch of units
		! mixing ratio = g(h2o) /kg(dry air) = 0.001 g(h2o)/g(dry air)
		! NOTE: need to use the molecular weight of water, not dry air
			
			
!		  u_h2o(m) = umix_h2o(m)*0.001   * 1.0e+05   *  1.0e+06      *28.9644*patm(m)/(temp(m)*8.314e+07)

		! was 1.0e-4
		  u_h2o(m) = umix_h2o(m)*0.001  *  1.0e-01   *  1000.        *(28.97)*ppa/(8.3144*temp(m))
		  umix_o3 = uppmv_o3(m)*(0.001)*(48./28.97)	! Convert ppmv to g(O3)/kg(dry air)
		  u_o3(m) = umix_o3*(0.001)*(1.0e-04)*(1000.)*(28.97)*ppa/(8.3144*temp(m))

		enddo

		do m=1, use_layers
			dzz(m) = abs(zlev(m+1) - zlev(m))
		end do

		call ch1(u_h2o, patm, temp, u_o3, dzz, modis_channel(1))
		call ch2(u_h2o, patm, temp, u_o3, dzz, modis_channel(2))
		call ch3(u_h2o, patm, temp, u_o3, dzz, modis_channel(3))
		call ch4(u_h2o, patm, temp, u_o3, dzz, modis_channel(4))
		call ch5b(u_h2o, patm, temp, u_o3, dzz, u_co2, modis_channel(5))
		call ch5a(u_h2o, patm, temp, u_o3, dzz, modis_channel(6))
		call ch06(u_h2o, patm, temp, u_o3, dzz, u_co2, modis_channel(7))
		call ch07(u_h2o, patm, temp, u_o3, dzz, u_co2, modis_channel(8))
		call ch8(u_h2o, patm, temp, u_o3, dzz, modis_channel(9))
		call ch9(u_h2o, patm, temp, u_o3, dzz, modis_channel(10))		
		call ch17(u_h2o, patm, temp, u_o3, dzz, modis_channel(11))
		call ch18(u_h2o, patm, temp, u_o3, dzz, modis_channel(12))
		call ch19(u_h2o, patm, temp, u_o3, dzz, modis_channel(13))
		call ch20(u_h2o, patm, temp, u_o3, dzz, u_co2, modis_channel(14))
		call ch22(u_h2o, patm, temp, u_o3, dzz, u_co2, modis_channel(15))
		call ch26(u_h2o, patm, temp, u_o3, dzz, modis_channel(16))
		call ch27(u_h2o, patm, temp, u_o3, dzz, modis_channel(17))
		call ch28(u_h2o, patm, temp, u_o3, dzz, modis_channel(18))
		call ch29(u_h2o, patm, temp, u_o3, dzz, modis_channel(19))
		call ch31(u_h2o, patm, temp, u_o3, dzz, u_co2, modis_channel(20))
		call ch32(u_h2o, patm, temp, u_o3, dzz, u_co2, modis_channel(21))
		call ch33(u_h2o, patm, temp, u_o3, dzz, u_co2, modis_channel(22))
		call ch34(u_h2o, patm, temp, u_o3, dzz, u_co2, modis_channel(23))  
		call ch35(u_h2o, patm, temp, u_o3, dzz, u_co2, modis_channel(24))
		call ch36(u_h2o, patm, temp, u_o3, dzz, u_co2, modis_channel(25))
	
	
	end subroutine run_corrK




	subroutine init_corrK
	
		use CKD_common
		use modis1
		use modis2
		use modis3
		use modis4
		use modis19
		use modis5b
		use modis5a 
		use modis26
		use modis6
		use modis7
		use modis8
		use modis9
		use modis17
		use modis18
		use modis20
		use modis22
		use modis27
		use modis28
		use modis29
		use modis31
		use modis32
		use modis33
		use modis34
		use modis35
		use modis36
	
	
	
		call init_ch1(modis_channel(1))
		call init_ch2(modis_channel(2))
		call init_ch3(modis_channel(3))
		call init_ch4(modis_channel(4))
		call init_ch5b(modis_channel(5))
		call init_ch5a(modis_channel(6))
		call init_ch06(modis_channel(7))
		call init_ch07(modis_channel(8))
		call init_ch8(modis_channel(9))
		call init_ch9(modis_channel(10))
		call init_ch17(modis_channel(11))
		call init_ch18(modis_channel(12))
		call init_ch19(modis_channel(13))
		call init_ch20(modis_channel(14))
		call init_ch22(modis_channel(15))
		call init_ch26(modis_channel(16))
		call init_ch27(modis_channel(17))
		call init_ch28(modis_channel(18))
		call init_ch29(modis_channel(19))
		call init_ch31(modis_channel(20))
		call init_ch32(modis_channel(21))
		call init_ch33(modis_channel(22))
		call init_ch34(modis_channel(23))
		call init_ch35(modis_channel(24))
		call init_ch36(modis_channel(25))
		
	end subroutine init_corrK
	  
	  
	subroutine cleanup_corrK
	
		integer :: i
		
		
		do i=1, 25
			deallocate(modis_channel(i)%weights, modis_channel(i)%taus)		
		end do
	

	end subroutine cleanup_corrK
	
end module corrK_profile

