! Created by Gala Wind 11/6/06
!  global variables had been moved out from the main code into this file so one doesn't have to hunt for them anymore

 module dis_global_pars

	implicit none
	public

	integer, parameter :: maxulv	= 27 ! computational levels
	integer, parameter :: maxcly	= 26  ! computational layers
!	integer, parameter :: maxcmu	= 16 ! computational polar angles
!	integer, parameter :: maxcmu	= 32 ! computational polar angles
       integer, parameter :: maxcmu    = 64 ! computational polar angles
	integer, parameter :: maxphi	= 5  ! output azimuthal angles
	integer, parameter :: maxumu	= 5 ! max number of output polar angles
	integer, parameter :: maxch		= 20  ! max number of imager channels
	integer, parameter :: maxtau	= 50  ! max number of optical thicknesses
	integer, parameter :: maxclds	= 5   ! max number of cloud layers
	integer, parameter :: maxszen	= 5  ! max number of solar zenith angles
	integer, parameter :: maxk		= 30  ! max number of k's per channel
!	integer, parameter :: maxmom	= 16  ! max number of Legendre expansion coefficients
!	integer, parameter :: maxmom	= 32  ! max number of Legendre expansion coefficients
       integer, parameter :: maxmom    = 64  ! max number of Legendre expansion coefficients

!	Out_Case:
!		3	Sixteen bands (0.65, 0.86, 0.94, 1.24b, 1.24a, 1.38, 1.64, 2.13, 3.75, 8.52, 11.0, 12.0, 13.3, 13.6, 13.9, 14.2)
!		5	Three bands (0.66, 0.86, 1.38)
!		6	Two Bands (0.66, 1.38)
!		7	Five Bands (0.66, 0.86, 1.24b, 1.24a, 1.38)
	integer, parameter :: out_case = 7 ! output case, got really sick of having to change it deep in code ...

!	real :: accur  = 1.0e-6
!	real :: c1 = 1.1910659e-8	! used to convert from radiance to brightness temperature
!	real :: c2 = 1.43878		! used to convert from radiance to brightness temperature
        real::WIND_VEC

 end module dis_global_pars
