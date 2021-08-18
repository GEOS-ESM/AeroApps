module modis6

	use CKD_common

	implicit none
	private
	
	integer, parameter :: ni = 2
	integer, parameter :: nco2 = 2
	integer, parameter :: nch4 = 2
	
	public :: ch06, init_ch06

contains

	subroutine init_ch06(modis_channel)
	
		type(ckd_type), intent(inout) :: modis_channel
		integer :: all_weights
		
		all_weights = ni*nco2*nch4
		
		modis_channel%num_weights = all_weights
		allocate(modis_channel%weights(all_weights))
		allocate(modis_channel%taus(all_weights, max_layers))
	
	end subroutine init_ch06


      subroutine ch06(u0, p0, t0, ux, fac, u_co2, modis_channel)
!
!      This routine was developed for MODIS Channel 6 
!                        6050--6150 cm^{-1}
!
!-----------------------------------------------------------------
!     Channel 6   1.628 - 1.652 microns
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
				 uco2(nlay), uch4(nlay)
				
		real :: tau(ni,nlay),tauco2(nco2,nlay), tauch4(nch4, nlay), &
				f(ni), fco2(nco2), fch4(nch4)
		
		integer :: i,m,ia, ib, index

	  
      do  m=1,use_layers
        t(m)=(t0(m)+t0(m+1))/2.0
        p(m)=(p0(m)+p0(m+1))/2.0
        u(m)=fac(m)*(u0(m)+u0(m+1))/2.0
        uco2(m)=p(m)*1.01325e+06/(1.3805e-16*t(m))*u_co2*fac(m)* &
     &       1.0e+05/6.023e+23*44.00995
        uch4(m)=p(m)*1.01325e+06/(1.3805e-16*t(m))*1.75e-06*fac(m)* &
     &       1.0e+05/6.023e+23*16.04303
	   enddo
 
      call ck_6(u,f,p,t,tau)
      call ckco2_6(uco2,fco2,p,t,tauco2)
      call ckch4_6(uch4,fch4,p,t,tauch4)

       do  m=1,use_layers
          do ib=1,nch4 ! 2 k's for ch4
         do ia=1,nco2 ! 2 k's for co2
        do i=1,ni ! 2 k's for h2o
			index = (ib-1)*nco2*ni + (ia-1)*ni + i 
            modis_channel%taus(index,m)=tau(i,m)+tauco2(ia,m)+tauch4(ib,m)
 			if (m == 1) modis_channel%weights(index) = f(i)*fco2(ia)*fch4(ib)
          end do
         end do
        end do
	  enddo
	  
      end subroutine ch06

! *********************************************************************
!     h2o 6050-6150 cm^{-1}
      subroutine ck_6(u,f,p,t,tau)
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		
		integer :: m 

		real :: k(ni)
		real :: coefk(ni,3,num_pressures)
		integer :: i, jp, jt

      f(1)=0.986574
      f(2)=0.013426
      k(1)=0.000327804908
      do i=2,ni
        k(i)=129.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/ &
      0.03474, 0.03563, 0.03707, 0.03928, 0.04282, 0.04828,&
      0.05663, 0.06935, 0.08834, 0.11607, 0.15540, 0.21171,&
      0.28868, 0.39289, 0.52722, 0.69646, 0.89905, 1.11896,&
      1.34189,&
      7.436E-04, 7.440E-04, 7.451E-04, 7.481E-04, 7.523E-04, 7.624E-04,&
      7.827E-04, 8.208E-04, 8.823E-04, 9.949E-04, 1.169E-03, 1.415E-03,&
      1.778E-03, 2.207E-03, 2.645E-03, 3.174E-03, 3.731E-03, 4.456E-03,&
      4.860E-03,&
      4.084E-06, 4.094E-06, 4.072E-06, 4.103E-06, 4.063E-06, 4.028E-06,&
      3.925E-06, 3.769E-06, 3.313E-06, 2.878E-06, 3.216E-06, 3.600E-06,&
      4.825E-06, 6.969E-06, 1.174E-05, 1.872E-05, 2.565E-05, 3.124E-05,&
      3.736E-05/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/ &
      1.56278, 1.56201, 1.56263, 1.56069, 1.55946, 1.55468,&
      1.55112, 1.54374, 1.53275, 1.51790, 1.49482, 1.46444,&
      1.41974, 1.36157, 1.28348, 1.18680, 1.07186, 0.94606,&
      0.82335,&
      3.103E-03, 3.080E-03, 3.079E-03, 3.103E-03, 3.077E-03, 3.072E-03,&
      3.058E-03, 3.058E-03, 3.016E-03, 2.932E-03, 2.827E-03, 2.694E-03,&
      2.499E-03, 2.281E-03, 2.020E-03, 1.706E-03, 1.376E-03, 9.839E-04,&
      7.570E-04,&
      2.702E-05, 2.756E-05, 2.640E-05, 2.703E-05, 2.645E-05, 2.769E-05,&
      2.664E-05, 2.751E-05, 2.802E-05, 2.738E-05, 2.795E-05, 2.687E-05,&
      2.695E-05, 2.522E-05, 2.238E-05, 1.872E-05, 1.449E-05, 1.216E-05,&
      7.525E-06/
	  
		call apply_ckd(ni, k, u, f, p, t, tau, coefk, .false.)
      end subroutine ck_6

! *********************************************************************
!     co2 6050--6150 cm^{-1}
      subroutine ckco2_6(u,f,p,t,tau)
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		
		integer :: m 

		real :: k(nco2)
		real :: coefk(nco2,3,num_pressures)
		integer :: i, jp, jt


      f(1)=0.973754
      f(2)=0.026246
      k(1)=0.00365213197
      do i=2,nco2
        k(i)=41.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 0.02755,0.02812,0.02899,0.03033,0.03354,0.03987, &
     & 0.04931,0.06348,0.08358,0.11188,0.15142,0.20406, &
     & 0.27136,0.36040,0.48384,0.65701,0.87761,1.12426, &
     &  1.36030, &
     &  3.882E-04, 3.866E-04, 3.844E-04, 3.816E-04, 3.489E-04, 3.314E-04, &
     &  3.048E-04, 2.773E-04, 2.411E-04, 1.819E-04, 9.912E-05,-5.625E-06, &
     & -1.331E-04,-2.962E-04,-5.435E-04,-8.391E-04,-1.092E-03,-1.166E-03, &
     & -1.094E-03, &
     &  1.919E-06, 1.922E-06, 1.947E-06, 1.972E-06, 2.128E-06, 2.197E-06, &
     &  2.269E-06, 2.244E-06, 2.403E-06, 2.634E-06, 2.816E-06, 3.116E-06, &
     &  3.409E-06, 4.050E-06, 4.863E-06, 5.247E-06, 5.209E-06, 3.475E-06, &
     &  3.700E-06/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 1.87906,1.87779,1.87580,1.87270,1.87004,1.86475, &
     & 1.85568,1.84316,1.82581,1.80088,1.76483,1.71762, &
     & 1.65694,1.57502,1.46330,1.30763,1.10767,0.88546, &
     &  0.67149, &
     & -7.305E-04,-7.251E-04,-7.170E-04,-7.052E-04,-7.168E-04,-6.951E-04, &
     & -6.920E-04,-6.605E-04,-6.202E-04,-5.800E-04,-5.130E-04,-3.915E-04, &
     & -2.803E-04,-1.144E-04, 8.787E-05, 3.511E-04, 5.828E-04, 6.604E-04, &
     &  5.943E-04, &
     &  7.688E-07, 7.344E-07, 6.812E-07, 6.000E-07,-1.438E-07,-2.281E-07, &
     &  3.812E-07, 4.312E-07, 2.375E-07,-7.000E-07,-2.125E-07,-5.750E-07, &
     & -1.156E-06,-9.906E-07,-1.228E-06,-2.259E-06,-2.100E-06,-1.341E-06, &
     & -8.813E-07/

		call apply_ckd(nco2, k, u, f, p, t, tau, coefk, .false.)
      end subroutine ckco2_6


! *********************************************************************
!     ch4 6050--6150 cm^{-1}
      subroutine ckch4_6(u,f,p,t,tau)
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		
		integer :: m 

		real :: k(nch4)
		real :: coefk(nch4,3,num_pressures)
		integer :: i, jp, jt
      f(1)=0.985628
      f(2)=0.014372
      k(1)=4.56504146
      do i=2,nch4
        k(i)=108.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 0.34337,0.34363,0.34404,0.34468,0.34567,0.34735, &
     & 0.35229,0.35996,0.37177,0.38976,0.41569,0.45088, &
     & 0.49796,0.56101,0.64618,0.76064,0.91324,1.09296, &
     &  1.28378, &
     &  1.099E-03, 1.097E-03, 1.095E-03, 1.092E-03, 1.088E-03, 1.064E-03, &
     &  1.030E-03, 9.907E-04, 9.285E-04, 8.651E-04, 7.820E-04, 7.284E-04, &
     &  7.160E-04, 7.654E-04, 8.562E-04, 8.891E-04, 9.259E-04, 1.099E-03, &
     &  1.221E-03, &
     & -4.241E-06,-4.228E-06,-4.203E-06,-4.175E-06,-4.122E-06,-3.697E-06, &
     & -3.794E-06,-3.537E-06,-3.106E-06,-2.966E-06,-2.700E-06,-2.297E-06, &
     & -1.756E-06,-2.303E-06,-2.737E-06,-1.384E-06,-4.747E-06,-8.116E-06, &
     & -7.209E-06/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 1.41808,1.41757,1.41677,1.41550,1.41350,1.41041, &
     & 1.40617,1.39946,1.38912,1.37525,1.35445,1.32749, &
     & 1.29584,1.25582,1.20426,1.14039,1.05227,0.94233, &
     &  0.81521, &
     &  6.508E-04, 6.535E-04, 6.579E-04, 6.645E-04, 6.746E-04, 6.855E-04, &
     &  7.056E-04, 7.397E-04, 7.717E-04, 8.509E-04, 9.021E-04, 9.507E-04, &
     &  9.903E-04, 9.651E-04, 8.978E-04, 8.664E-04, 8.205E-04, 6.340E-04, &
     &  4.632E-04, &
     & -5.619E-06,-5.637E-06,-5.672E-06,-5.713E-06,-5.778E-06,-5.788E-06, &
     & -5.991E-06,-6.206E-06,-6.025E-06,-6.716E-06,-6.391E-06,-6.762E-06, &
     & -7.856E-06,-7.603E-06,-6.525E-06,-7.497E-06,-5.900E-06,-3.794E-06, &
     & -3.687E-06/
		call apply_ckd(nco2, k, u, f, p, t, tau, coefk, .false.)
      end subroutine ckch4_6

end module modis6