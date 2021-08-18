module modis7

	use CKD_common

	implicit none
	private
	
	integer, parameter :: ni = 4
	integer, parameter :: nco2 = 1
	integer, parameter :: nch4 = 1
	integer, parameter :: nn2o = 1
	
	public :: ch07, init_ch07
	
contains

		subroutine init_ch07(modis_channel)
	
		type(ckd_type), intent(inout) :: modis_channel
		integer :: all_weights
		
		all_weights = ni*nco2*nch4*nn2o
		
		modis_channel%num_weights = all_weights
		allocate(modis_channel%weights(all_weights))
		allocate(modis_channel%taus(all_weights, max_layers))
	
	end subroutine init_ch07




      subroutine ch07(u0, p0, t0, ux, fac, u_co2, modis_channel)
!
!      This routine was developed for MODIS Channel 7 
!                        4640--4750 cm^{-1}
!
!-----------------------------------------------------------------
!     Channel 7   2.105-2.155 microns
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
				 uco2(nlay), uch4(nlay), un2o(nlay)
				
		real :: tau(ni,nlay),tauco2(nco2,nlay), tauch4(nch4, nlay), &
				taun2o(nn2o, nlay), f(ni), fco2(nco2), fch4(nch4), fn2o(nn2o)
		
		integer :: i,m,ia, index
	 
      do  m=1,use_layers
        t(m)=(t0(m)+t0(m+1))/2.0
        p(m)=(p0(m)+p0(m+1))/2.0
        u(m)=fac(m)*(u0(m)+u0(m+1))/2.0
        uco2(m)=p(m)*1.01325e+06/(1.3805e-16*t(m))*u_co2*fac(m)* &
     &       1.0e+05/6.023e+23*44.00995
        uch4(m)=p(m)*1.01325e+06/(1.3805e-16*t(m))*1.75e-06*fac(m)* &
     &       1.0e+05/6.023e+23*16.04303
        un2o(m)=p(m)*1.01325e+06/(1.3805e-16*t(m))*3.10e-07*fac(m)* &
     &       1.0e+05/6.023e+23*44.0128
      enddo
	  
      call ck_7(u,f,p,t,tau)
      call ckco2_7(uco2,fco2,p,t,tauco2)
      call ckch4_7(uch4,fch4,p,t,tauch4)
      call ckn2o_7(un2o,fn2o,p,t,taun2o)
	  
      do  m=1,use_layers
         do ia=1,nco2 ! 1 k for co2
        do i=1,ni ! 4 k's for h2o
			index = (ia-1)*ni + i 
            modis_channel%taus(index,m)=tau(i,m)+tauco2(ia,m)+tauch4(1,m)+ taun2o(1,m)
 			if (m == 1) modis_channel%weights(index) = f(i)*fco2(ia)*fch4(1)*fn2o(1)
         end do
        end do
       enddo
	   
 
	end subroutine ch07


! *********************************************************************
!     h2o 4640--4750 cm^{-1}
      subroutine ck_7(u,f,p,t,tau)
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		
		
		integer :: m 

		real :: k(ni)
		real :: coefk(ni,3,num_pressures)
		integer :: i, jp, jt

      f(1)=0.895296
      f(2)=0.031732
      f(3)=0.060290
      f(4)=0.012682
      k(1)=0.00211455525
      do i=2,ni
        k(i)=7.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     0.00069, 0.00108, 0.00171, 0.00269, 0.00425, 0.00671,&
     0.01060, 0.01675, 0.02647, 0.04179, 0.06596, 0.10392,&
     0.16230, 0.25038, 0.38134, 0.57175, 0.84015, 1.22067,&
     1.74323,&
     3.750E-07, 5.000E-07, 8.750E-07, 1.375E-06, 2.375E-06, 4.250E-06,&
     6.875E-06, 1.187E-05, 2.025E-05, 3.512E-05, 6.263E-05, 1.126E-04,&
     2.126E-04, 3.636E-04, 6.051E-04, 9.654E-04, 1.209E-03, 1.097E-03,&
     4.598E-04,&
     9.375E-09, 1.875E-08, 2.188E-08, 4.063E-08, 5.937E-08, 9.375E-08,&
     1.406E-07, 2.094E-07, 3.125E-07, 4.531E-07, 6.781E-07, 8.969E-07,&
     1.384E-06, 2.578E-06, 5.434E-06, 1.045E-05, 2.233E-05, 3.360E-05,&
     4.298E-05/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/&
     0.00123, 0.00190, 0.00295, 0.00458, 0.00712, 0.01102,&
     0.01704, 0.02636, 0.04086, 0.06325, 0.09715, 0.14691,&
     0.21613, 0.31077, 0.44054, 0.61369, 0.85680, 1.17226,&
     1.45661,&
     1.000E-06, 1.625E-06, 2.625E-06, 4.625E-06, 8.125E-06, 1.400E-05,&
     2.250E-05, 3.950E-05, 7.000E-05, 1.236E-04, 2.161E-04, 3.645E-04,&
     6.036E-04, 9.180E-04, 1.294E-03, 1.210E-03, 5.669E-04,-3.511E-04, &
     -1.525E-03, &
     6.250E-09, 1.562E-08, 1.562E-08, 9.375E-09,-9.375E-09,-1.250E-08, &
     -1.250E-08,-2.500E-08,-1.125E-07,-3.594E-07,-4.406E-07, 2.938E-07,&
     2.916E-06, 7.337E-06, 1.317E-05, 2.959E-05, 3.560E-05, 2.588E-05,&
     4.391E-05/
      data ( (coefk(3,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     0.03138, 0.03287, 0.03521, 0.03873, 0.04411, 0.05215,&
     0.06366, 0.08001, 0.10334, 0.13644, 0.18287, 0.24898,&
     0.34341, 0.46613, 0.61364, 0.78809, 0.94744, 1.06086, &
     1.12761, &
     6.776E-04, 6.746E-04, 6.709E-04, 6.678E-04, 6.681E-04, 6.753E-04, &
     6.887E-04, 7.216E-04, 7.675E-04, 8.051E-04, 8.102E-04, 7.276E-04,&
     5.464E-04, 2.079E-04,-3.499E-04,-9.739E-04,-1.483E-03,-1.799E-03,&
     -2.084E-03,&
     5.453E-06, 5.472E-06, 5.478E-06, 5.488E-06, 5.459E-06, 5.375E-06, &
     5.362E-06, 5.691E-06, 6.038E-06, 6.816E-06, 9.019E-06, 1.165E-05,&
     1.228E-05, 1.343E-05, 1.709E-05, 1.315E-05, 1.144E-05, 1.051E-05,&
     6.600E-06/
      data ( (coefk(4,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     1.93703, 1.93581, 1.93390, 1.93101, 1.92864, 1.92198,&
     1.91236, 1.90053, 1.88250, 1.85555, 1.81815, 1.76413,&
     1.68470, 1.57802, 1.44346, 1.27665, 1.09886, 0.92557,&
     0.75835,&
     -2.824E-03,-2.822E-03,-2.820E-03,-2.818E-03,-2.820E-03,-2.802E-03,&
     -2.842E-03,-2.846E-03,-2.889E-03,-2.931E-03,-2.934E-03,-2.909E-03,&
     -2.787E-03,-2.595E-03,-2.306E-03,-1.944E-03,-1.621E-03,-1.328E-03,&
     -9.416E-04,&
     1.605E-05, 1.603E-05, 1.603E-05, 1.601E-05, 1.473E-05, 1.539E-05,&
     1.605E-05, 1.510E-05, 1.479E-05, 1.415E-05, 1.301E-05, 1.082E-05,&
     1.018E-05, 8.159E-06, 5.484E-06, 5.941E-06, 4.753E-06, 2.184E-06,&
     1.416E-06/
		call apply_ckd(ni, k, u, f, p, t, tau, coefk, .true.)
      end subroutine ck_7
	  
! *********************************************************************
!     co2 4640--4750 cm^{-1}
      subroutine ckco2_7(u,f,p,t,tau)
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		
		integer :: m 

		real :: k(nco2)
		real :: coefk(nco2,3,num_pressures)
		integer :: i, jp, jt
      f(1)=1.00000
      k(1)=0.007121588379*2.0790
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     1.00131, 1.00091, 1.00028, 0.99941, 0.99941, 0.99948,&
     0.99953, 0.99963, 0.99977, 1.00000, 1.00033, 1.00071,&
     1.00061, 1.00020, 1.00005, 1.00011, 0.99999, 0.99997,&
     0.99991,&
     9.306E-03, 9.305E-03, 9.303E-03, 9.292E-03, 9.285E-03, 9.286E-03,&
     9.285E-03, 9.286E-03, 9.287E-03, 9.288E-03, 9.289E-03, 9.296E-03,&
     9.305E-03, 9.298E-03, 9.292E-03, 9.295E-03, 9.291E-03, 9.291E-03,&
     9.291E-03,&
     8.505E-05, 8.502E-05, 8.498E-05, 8.504E-05, 8.490E-05, 8.491E-05,&
     8.491E-05, 8.493E-05, 8.496E-05, 8.499E-05, 8.505E-05, 8.505E-05,&
     8.507E-05, 8.502E-05, 8.502E-05, 8.501E-05, 8.497E-05, 8.497E-05,&
     8.496E-05/

		call apply_ckd(nco2, k, u, f, p, t, tau, coefk, .true.)
      end subroutine ckco2_7

! *********************************************************************
!     ch4 4640--4750 cm^{-1}
      subroutine ckch4_7(u,f,p,t,tau)
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		
		integer :: m 

		real :: k(nch4)
		real :: coefk(nch4,3,num_pressures)
		integer :: i, jp, jt
      f(1)=1.00000
      k(1)=0.1386447515*0.80347
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/ &
      1.00224, 1.00206, 1.00179, 1.00135, 1.00066, 0.99989, &
      1.00000, 1.00001, 1.00008, 1.00020, 1.00039, 1.00067,&
      1.00106, 1.00147, 1.00102, 1.00069, 1.00024, 0.99979,&
      0.99972,&
      1.228E-02, 1.228E-02, 1.227E-02, 1.227E-02, 1.227E-02, 1.225E-02,&
      1.225E-02, 1.225E-02, 1.225E-02, 1.225E-02, 1.225E-02, 1.225E-02,&
      1.226E-02, 1.226E-02, 1.227E-02, 1.226E-02, 1.225E-02, 1.225E-02,&
      1.224E-02,&
      3.265E-05, 3.266E-05, 3.265E-05, 3.266E-05, 3.266E-05, 3.269E-05,&
      3.253E-05, 3.256E-05, 3.257E-05, 3.256E-05, 3.256E-05, 3.254E-05,&
      3.254E-05, 3.246E-05, 3.267E-05, 3.255E-05, 3.255E-05, 3.263E-05,&
      3.258E-05/
	 
		call apply_ckd(nch4, k, u, f, p, t, tau, coefk, .true.)
      end subroutine ckch4_7
! *********************************************************************
!     n2o 4640--4750 cm^{-1}
      subroutine ckn2o_7(u,f,p,t,tau)
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		
		integer :: m 

		real :: k(nn2o)
		real :: coefk(nn2o,3,num_pressures)
		integer :: i, jp, jt
      f(1)=1.00000
      k(1)=4.8762226564*0.93892
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/ &
      1.00138, 1.00096, 1.00030, 0.99950, 0.99953, 0.99959,&
      0.99965, 0.99975, 0.99990, 1.00013, 1.00047, 1.00083,&
      1.00054, 1.00036, 1.00009, 1.00012, 1.00000, 0.99997,&
      0.99993,&
     -6.287E-05,-6.050E-05,-5.688E-05,-6.338E-05,-6.575E-05,-6.663E-05,&
     -6.650E-05,-6.712E-05,-6.775E-05,-6.887E-05,-7.013E-05,-6.662E-05,&
     -5.950E-05,-6.488E-05,-6.438E-05,-6.200E-05,-6.525E-05,-6.562E-05,&
     -6.588E-05,&
     -1.359E-06,-1.369E-06,-1.384E-06,-1.259E-06,-1.325E-06,-1.334E-06,&
     -1.325E-06,-1.316E-06,-1.300E-06,-1.278E-06,-1.259E-06,-1.309E-06,&
     -1.275E-06,-1.278E-06,-1.284E-06,-1.294E-06,-1.312E-06,-1.316E-06,&
     -1.316E-06/
		call apply_ckd(nn2o, k, u, f, p, t, tau, coefk, .true.)
      end subroutine ckn2o_7

end module modis7
