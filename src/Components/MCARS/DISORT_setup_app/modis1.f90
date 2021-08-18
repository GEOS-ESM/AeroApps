module modis1

	use CKD_common

	implicit none
	private
	
	integer, parameter :: ni = 3
	integer, parameter :: no2 = 2
	integer, parameter :: no3 = 1
	
	public :: ch1, init_ch1

contains

	subroutine init_ch1(modis_channel)
	
		type(ckd_type), intent(inout) :: modis_channel
		integer :: all_weights
		
		all_weights = ni*no2*no3
		
		modis_channel%num_weights = all_weights
		allocate(modis_channel%weights(all_weights))
		allocate(modis_channel%taus(all_weights, max_layers))
			
	end subroutine init_ch1

	subroutine ch1(u0, p0, t0, ux, fac, modis_channel)
!
!      This routine was developed for MODIS Channel 1 
!                        14800--16200 cm^{-1}
!
!-----------------------------------------------------------------
!     Channel 1   0.62-0.67 microns
!-----------------------------------------------------------------
!  INPUTS:
!     u0 --> water vapor amount [g/cm2/km]
!     p0 --> pressure in atmospheres
!     t0 --> temp in K
!     ux --> ozone [g/cm2/km]
!     fac --> layer thickness in km
!-----------------------------------------------------------------
		real, dimension(:), intent(in) :: u0, p0, t0, ux, fac
		type(ckd_type), intent(inout) :: modis_channel
		
		integer, parameter :: nlev = max_levels
		integer, parameter :: nlay = max_layers
		integer, parameter :: mlv = nlay

		real :: p(nlay), t(nlay),u(nlay), &
				uo2(nlay),uo3(nlay)
				
		real :: tau(ni,nlay),tauo3(no3,nlay),tauo2(no2,nlay), &
				f(ni), fo2(no2), fo3(no3)
		
		integer :: i,m,ia, index
	  
		do m=1,use_layers
			t(m)=(t0(m)+t0(m+1))/2.0
			p(m)=(p0(m)+p0(m+1))/2.0
			u(m)=fac(m)*(u0(m)+u0(m+1))/2.0
			uo3(m)=fac(m)*(ux(m)+ux(m+1))/2.0
			uo2(m)=p(m)*1.01325e+06/(1.3805e-16*t(m))*0.20954*fac(m)*1.0e+05/6.023e+23*31.9988
		enddo
	  
		call ck_1(u,f,p,t,tau) 
		call cko3_1(uo3,fo3,p,t,tauo3)
		call cko2_1(uo2,fo2,p,t,tauo2)
!     **fo3=1.0000**

		do m=1,use_layers
			do ia=1,no2 ! 2 k's for o2
				do i=1,ni ! 3 k's for h2o; 1 k for o3
					index = (ia-1)*ni + i
					modis_channel%taus(index, m)=tau(i,m)+tauo3(1,m)+tauo2(ia,m)
					if (m == 1) modis_channel%weights(index) = f(i)*fo3(1)*fo2(ia)
				end do
			end do
		enddo
			  	  
				  
     end subroutine ch1

! *********************************************************************
	subroutine ck_1( u, f, p, t, tau)
	  
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		
		
		real :: k(ni)
		real :: coefk(ni,3,num_pressures)
		integer :: i, jp, jt

		f(1)=0.923520
		f(2)=0.064941
		f(3)=0.011539

		k(1)=0.000767477499
		do i=2,ni
			k(i)=14.0*k(i-1)
		end do

		data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 0.00010,0.00015,0.00021,0.00031,0.00046,0.00069, &
     & 0.00106,0.00501,0.01366,0.02855,0.05267,0.09113, &
     & 0.15111,0.24486,0.38769,0.59169,0.86087,1.17418, &
     &  1.49737, &
     &  2.063E-06, 2.663E-06, 3.525E-06, 4.813E-06, 6.737E-06, 9.588E-06, &
     &  6.375E-06,-1.487E-05,-2.263E-05,-4.325E-05,-7.950E-05,-1.385E-04, &
     & -2.294E-04,-3.700E-04,-5.548E-04,-7.496E-04,-9.725E-04,-1.093E-03, &
     & -1.084E-03, &
     &  1.406E-08, 1.469E-08, 1.562E-08, 1.594E-08, 1.469E-08, 1.219E-08, &
     &  1.906E-07, 4.063E-08, 5.938E-08, 9.375E-08, 1.438E-07, 1.437E-07, &
     &  2.219E-07, 1.313E-07,-3.500E-07,-1.072E-06,-1.994E-06,-1.744E-06, &
     & -1.875E-07/
		data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 0.42817,0.42864,0.42938,0.43055,0.43241,0.43577, &
     & 0.44078,0.45041,0.46498,0.48683,0.51860,0.56501, &
     & 0.62725,0.70676,0.79783,0.89069,0.97527,1.04583, &
     &  1.08525, &
     &  1.569E-03, 1.567E-03, 1.569E-03, 1.562E-03, 1.552E-03, 1.532E-03, &
     &  1.499E-03, 1.431E-03, 1.357E-03, 1.244E-03, 1.081E-03, 8.661E-04, &
     &  6.213E-04, 3.829E-04, 1.822E-04, 9.425E-05, 5.775E-05,-2.000E-06, &
     & -5.437E-05, &
     & -3.806E-06,-3.794E-06,-3.619E-06,-3.588E-06,-3.547E-06,-3.628E-06, &
     & -3.256E-06,-2.969E-06,-2.575E-06,-2.272E-06,-1.406E-06,-8.781E-07, &
     &  3.562E-07, 9.531E-07, 1.531E-06, 2.131E-06, 2.663E-06, 2.456E-06, &
     &  1.741E-06/
		data ( (coefk(3,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 1.66345,1.66292,1.66209,1.66078,1.65871,1.65715, &
     & 1.65372,1.64745,1.63853,1.62449,1.60054,1.56782, &
     & 1.51712,1.44758,1.35418,1.23261,1.08906,0.93147, &
     &  0.76853, &
     & -1.359E-03,-1.357E-03,-1.332E-03,-1.326E-03,-1.315E-03,-1.322E-03, &
     & -1.298E-03,-1.284E-03,-1.231E-03,-1.172E-03,-1.101E-03,-1.012E-03, &
     & -8.684E-04,-7.141E-04,-5.700E-04,-4.296E-04,-3.380E-04,-2.581E-04, &
     & -2.137E-04, &
     &  3.000E-06, 2.991E-06, 3.488E-06, 3.453E-06, 3.403E-06, 2.816E-06, &
     &  2.697E-06, 3.056E-06, 2.778E-06, 1.941E-06, 2.612E-06, 1.650E-06, &
     &  2.059E-06, 1.709E-06, 1.456E-06, 1.234E-06, 1.213E-06, 1.259E-06, &
     &  1.238E-06/

		call apply_ckd(ni, k, u, f, p, t, tau, coefk, .false.)

      end subroutine ck_1
	  
! *********************************************************************
      subroutine cko3_1(u,f,p,t,tau)
	
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		
		
		real :: k(no3)
		integer :: i, m 

		f(1)=1.000000
 
		k(1)=2.69e-21*6.023e+23/47.9982
		
		do i=1,no3
			do m=1,use_layers
				tau(i,m)=k(i)*u(m)
			end do
		end do
 
	 end subroutine cko3_1
	  
! *********************************************************************
      subroutine cko2_1( u,f,p,t,tau)

		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
				
		
		real :: k(no2)
		real :: coefk(no2,3,num_pressures)
		integer :: i, jp, jt

		f(1)=0.841000
		f(2)=0.159000
    
	    k(1)=0.0
		k(2)=0.0000347395293


      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000, &
     & 1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000, &
     & 1.0000, &
     & 0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00, &
     & 0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00, &
     & 0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00, &
     & 0.000e+00, &
     & 0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00, &
     & 0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00, &
     & 0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00, &
     & 0.000e+00/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 1.00183,1.00176,1.00164,1.00146,1.00117,1.00071, &
     & 0.99998,0.99917,0.99924,0.99929,0.99936,0.99948, &
     & 0.99965,0.99992,1.00029,1.00067,1.00018,0.99987, &
     &  0.99978, &
     & -2.022E-03,-2.022E-03,-2.021E-03,-2.020E-03,-2.018E-03,-2.014E-03, &
     & -2.009E-03,-2.016E-03,-2.017E-03,-2.018E-03,-2.018E-03,-2.019E-03, &
     & -2.020E-03,-2.022E-03,-2.024E-03,-2.020E-03,-2.015E-03,-2.018E-03, &
     & -2.018E-03, &
     &  5.991E-06, 5.984E-06, 5.984E-06, 5.978E-06, 5.966E-06, 5.953E-06, &
     &  5.928E-06, 6.028E-06, 5.978E-06, 5.981E-06, 5.987E-06, 5.988E-06, &
     &  6.000E-06, 6.009E-06, 6.013E-06, 5.906E-06, 6.006E-06, 6.022E-06, &
     &  5.994E-06/
	 

		call apply_ckd(no2, k, u, f, p, t, tau, coefk, .false.)

      end subroutine cko2_1

end module modis1
