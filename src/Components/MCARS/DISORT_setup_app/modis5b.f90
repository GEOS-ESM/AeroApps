module modis5b 

	use CKD_common

	implicit none
	private
	
	integer, parameter :: ni = 4
	integer, parameter :: nco2 = 1
	
	public :: ch5b, init_ch5b

contains

	subroutine init_ch5b(modis_channel)
	
		type(ckd_type), intent(inout) :: modis_channel
		integer :: all_weights
		
		all_weights = ni*nco2
		
		modis_channel%num_weights = all_weights
		allocate(modis_channel%weights(all_weights))
		allocate(modis_channel%taus(all_weights, max_layers))
	
	end subroutine init_ch5b

      subroutine ch5b(u0, p0, t0, ux, fac, u_co2, modis_channel)
!
!      This routine was developed for MODIS Channel 5b 
!                        8050--8210 cm^{-1}
!
!      Note that band 5 is split into two parts
!-----------------------------------------------------------------
!     Channel 5b   1.218-1.242 microns
!-----------------------------------------------------------------
!  INPUTS:
!     u0 --> water vapor amount [g/cm2/km]
!     p0 --> pressure in atmospheres
!     t0 --> temp in K
!     ux --> ozone [g/cm2/km]
!     dz --> layer thickness in km
!-----------------------------------------------------------------

		real, dimension(:), intent(in) :: u0, p0, t0, ux, fac
		real, intent(in) :: u_co2
		type(ckd_type), intent(inout) :: modis_channel
		
		integer, parameter :: nlev = max_levels
		integer, parameter :: nlay = max_layers
		integer, parameter :: mlv = nlay

		real :: p(nlay), t(nlay),u(nlay), &
				 uco2(nlay)
				
		real :: tau(ni,nlay),tauco2(nco2,nlay), &
				f(ni), fco2(nco2)
		
		integer :: i,m,ia, index


      do  m=1,use_layers
        t(m)=(t0(m)+t0(m+1))/2.0
        p(m)=(p0(m)+p0(m+1))/2.0
        u(m)=fac(m)*(u0(m)+u0(m+1))/2.0
        uco2(m)=p(m)*1.01325e+06/(1.3805e-16*t(m))*u_co2*fac(m)* &
            1.0e+05/6.023e+23*44.00995
      enddo
	  
      call ck_5b(u,f,p,t,tau)
      call ckco2_5b(uco2,fco2,p,t,tauco2)

	  
      do  m=1,use_layers
          do ia=1,nco2 ! 1 k for co2
        do i=1,ni ! 4 k's for h2o
			index = (ia-1)*ni + i
            modis_channel%taus(index,m)=tau(i,m)+tauco2(ia,m)
 			if (m == 1) modis_channel%weights(index) = f(i)*fco2(ia)
          end do
        end do
      enddo
	  
      end subroutine ch5b 


! *********************************************************************
      subroutine ck_5b(u,f,p,t,tau)
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		
		
		real :: k(ni)
		real :: coefk(ni,3,num_pressures)
		integer :: i, jp, jt

      f(1)=0.835682
      f(2)=0.067846
      f(3)=0.082925
      f(4)=0.013547
      k(1)=0.00143421052
      do i=2,ni
        k(i)=7.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00209,0.01027,0.02844,0.06026, &
     & 0.11102,0.18950,0.31197,0.50475,0.80525,1.26035, &
     &  1.89683, &
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, &
     &  0.000E+00, 0.000E+00, 4.625E-06, 4.100E-05, 1.044E-04, 2.141E-04, &
     &  3.889E-04, 6.705E-04, 1.109E-03, 1.727E-03, 2.691E-03, 4.202E-03, &
     &  6.561E-03, &
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, &
     &  0.000E+00, 0.000E+00,-8.438E-08,-1.563E-07,-3.594E-07,-5.844E-07, &
     & -9.969E-07,-1.019E-06,-1.528E-06,-2.578E-06,-4.778E-06,-9.828E-06, &
     & -1.635E-05/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00033,0.00618,0.01871,0.03807,0.06346,0.09913, &
     & 0.15404,0.23972,0.37061,0.56965,0.85077,1.17107, &
     &  1.49008, &
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, &
     & -1.175E-06, 2.912E-05, 1.066E-04, 1.872E-04, 2.869E-04, 4.284E-04, &
     &  6.340E-04, 9.260E-04, 1.355E-03, 1.958E-03, 2.987E-03, 4.232E-03, &
     &  5.651E-03, &
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, &
     & -2.687E-08,-1.719E-07,-1.719E-07,-1.750E-07,-2.187E-08, 2.344E-07, &
     &  5.813E-07, 3.312E-07,-3.281E-07,-3.606E-06,-8.787E-06,-8.500E-06, &
     & -6.112E-06/
      data ( (coefk(3,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 0.24758,0.24812,0.24897,0.25057,0.25266,0.25595, &
     & 0.26433,0.27716,0.29621,0.32389,0.36447,0.42291, &
     & 0.50184,0.60251,0.72205,0.84883,0.96324,1.04462, &
     &  1.08411, &
     &  2.888E-03, 2.888E-03, 2.888E-03, 2.888E-03, 2.888E-03, 2.877E-03, &
     &  2.852E-03, 2.849E-03, 2.841E-03, 2.827E-03, 2.828E-03, 2.850E-03, &
     &  2.947E-03, 3.192E-03, 3.605E-03, 4.108E-03, 4.537E-03, 4.762E-03, &
     &  4.703E-03, &
     &  3.947E-06, 3.944E-06, 3.941E-06, 3.775E-06, 3.787E-06, 4.062E-06, &
     &  3.863E-06, 3.884E-06, 3.788E-06, 3.681E-06, 3.562E-06, 2.781E-06, &
     &  2.516E-06, 2.138E-06, 3.719E-07,-1.250E-07, 8.125E-07, 2.937E-06, &
     &  3.369E-06/
      data ( (coefk(4,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 1.96323,1.96233,1.96089,1.96076,1.95720,1.95158, &
     & 1.94286,1.93139,1.91439,1.88623,1.84438,1.78479, &
     & 1.70211,1.59116,1.45255,1.28612,1.10238,0.91534, &
     &  0.74052, &
     &  5.494E-03, 5.494E-03, 5.495E-03, 5.497E-03, 5.499E-03, 5.503E-03, &
     &  5.483E-03, 5.490E-03, 5.527E-03, 5.511E-03, 5.478E-03, 5.418E-03, &
     &  5.274E-03, 5.006E-03, 4.545E-03, 3.905E-03, 3.253E-03, 2.637E-03, &
     &  2.155E-03, &
     & -4.625E-06,-4.628E-06,-4.625E-06,-5.963E-06,-5.981E-06,-5.978E-06, &
     & -5.400E-06,-5.369E-06,-5.894E-06,-5.112E-06,-4.353E-06,-4.431E-06, &
     & -4.994E-06,-3.103E-06,-2.500E-06,-2.169E-06,-1.728E-06,-2.262E-06, &
     & -1.850E-06/

		call apply_ckd(ni, k, u, f, p, t, tau, coefk, .false.)

	end subroutine ck_5b
!**********************************************************************
      subroutine ckco2_5b(u,f,p,t,tau)
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		
		integer :: m 

		real :: k(nco2)

      f(1)=1.000000
      do  m=1,use_layers
        k(1)=0.906589*3.91935654e-03*(1.0000+6.7076e-04*(t(m)-250.0)+ &
             2.5516e-06*(t(m)-250.0)**2)
        tau(1,m)=k(1)*u(m)
      end do

      end subroutine ckco2_5b
	  
end module modis5b
