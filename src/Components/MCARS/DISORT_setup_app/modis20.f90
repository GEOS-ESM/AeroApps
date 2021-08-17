module modis20

	use CKD_common

	implicit none
	private
	
	integer, parameter :: ni = 4
	integer, parameter :: nco2 = 1
	integer, parameter :: nch4 = 1
	
	public :: ch20, init_ch20
	
contains

	subroutine init_ch20(modis_channel)
	
		type(ckd_type), intent(inout) :: modis_channel
		integer :: all_weights
		
		all_weights = ni*nco2*nch4
		
		modis_channel%num_weights = all_weights
		allocate(modis_channel%weights(all_weights))
		allocate(modis_channel%taus(all_weights, max_layers))
	
	end subroutine init_ch20


      subroutine ch20(u0, p0, t0, ux, fac, u_co2, modis_channel)
!
!      This routine was developed for MODIS Channel 20 
!                        2605--2735 cm^{-1}
!
!-----------------------------------------------------------------
!     Channel 20   3.66 - 3.84 microns
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
				f(ni), fco2(nco2), fch4(nch4), xlog
		
		integer :: i,m,ia, ib, ic
	  
      do  m=1,use_layers
        t(m)=(t0(m)+t0(m+1))/2.0
        xlog=(log10(p0(m))+log10(p0(m+1)))/2.0
        p(m)=10.0**xlog
        xlog=(log10(u0(m))+log10(u0(m+1)))/2.0
        u(m)=fac(m)*10.0**xlog
        uco2(m)=p(m)*1.01325e+06/(1.3805e-16*t(m))*u_co2*fac(m)* &
      &          1.0e+05/6.023e+23*44.00995
        uch4(m)=p(m)*1.01325e+06/(1.3805e-16*t(m))*1.75e-06*fac(m)* &
      &          1.0e+05/6.023e+23*16.04303
      enddo
	  
      call ck_20(u,f,p,t,tau)
      call ckco2_20(uco2,fco2,p,t,tauco2)
      call ckch4_20(uch4,fch4,p,t,tauch4)


      do  m=1,use_layers
        do i=1,ni ! 4 k's for h2o
             modis_channel%taus(i,m)=tau(i,m)+tauco2(1,m) +tauch4(1,m)
  			if (m == 1)  modis_channel%weights(i) = f(i)*fco2(1)*fch4(1)
        end do
      enddo
	  

      end subroutine ch20

! *********************************************************************
!     h2o 2605-2735 cm^{-1}
      subroutine ck_20(u,f,p,t,tau)
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		
		

		real :: k(ni)
		real :: coefk(ni,3,num_pressures)
		integer :: i, jp, jt
      f(1)=0.728956
      f(2)=0.141171
      f(3)=0.119373
      f(4)=0.010500
      k(1)=0.0115054257
      do i=2,ni
        k(i)=6.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00012, &
     & 0.00074,0.00285,0.00892,0.02269,0.04525,0.08065, &
     & 0.13526,0.22093,0.35074,0.54211,0.81646,1.19221, &
     &  1.69219, &
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.625E-07, &
     &  1.538E-06, 4.875E-06, 9.625E-06, 1.662E-05, 3.087E-05, 5.187E-05, &
     &  7.863E-05, 1.079E-04, 1.395E-04, 1.262E-04, 6.125E-06,-2.690E-04, &
     & -6.966E-04, &
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,-5.938E-09, &
     & -3.438E-09,-4.063E-08,-8.437E-08,-2.219E-07,-3.469E-07,-5.031E-07, &
     & -7.094E-07,-8.719E-07,-1.125E-06,-1.069E-06,-1.659E-06,-2.206E-06, &
     & -5.872E-06/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 0.00000,0.00000,0.00000,0.00012,0.00094,0.00331, &
     & 0.00930,0.02075,0.03630,0.05716,0.08733,0.13215, &
     & 0.19848,0.29460,0.43142,0.61385,0.84648,1.14771, &
     &  1.47469, &
     &  0.000E+00, 0.000E+00, 0.000E+00, 1.000E-07, 3.538E-06, 1.162E-05, &
     &  2.975E-05, 4.613E-05, 5.625E-05, 7.250E-05, 9.225E-05, 1.191E-04, &
     &  1.225E-04, 1.065E-04, 6.537E-05, 1.388E-05,-1.099E-04,-1.771E-04, &
     &  1.487E-04, &
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,-1.531E-08,-5.313E-08, &
     & -8.750E-08,-2.656E-07,-2.937E-07,-4.000E-07,-2.938E-07,-1.781E-07, &
     & -1.625E-07,-3.750E-07,-1.816E-06,-2.647E-06,-3.609E-06,-6.459E-06, &
     & -6.444E-06/
      data ( (coefk(3,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 0.07742,0.07855,0.08033,0.08411,0.09114,0.10155, &
     & 0.11679,0.13925,0.17213,0.21965,0.28512,0.37317, &
     & 0.48731,0.62597,0.76722,0.88901,0.96953,0.99281, &
     &  0.95607, &
     &  6.610E-04, 6.584E-04, 6.538E-04, 6.236E-04, 6.039E-04, 5.701E-04, &
     &  5.288E-04, 4.759E-04, 4.120E-04, 3.609E-04, 3.032E-04, 2.133E-04, &
     &  1.156E-04, 2.575E-05, 1.263E-05, 3.625E-06, 5.250E-06,-3.562E-05, &
     & -4.963E-05, &
     &  9.187E-07, 9.219E-07, 9.125E-07, 9.219E-07, 8.594E-07, 8.531E-07, &
     &  9.375E-07, 8.219E-07, 5.812E-07,-2.781E-07,-5.250E-07,-6.500E-07, &
     & -1.159E-06,-2.906E-06,-4.241E-06,-4.109E-06,-3.875E-06,-2.253E-06, &
     & -2.166E-06/
      data ( (coefk(4,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 3.40794,3.40283,3.39483,3.38682,3.37222,3.35187, &
     & 3.32069,3.27570,3.20285,3.10289,2.95936,2.76430, &
     & 2.50119,2.17487,1.81644,1.45460,1.12661,0.86020, &
     &  0.64816, &
     & -1.411E-03,-1.395E-03,-1.368E-03,-1.331E-03,-1.278E-03,-1.314E-03, &
     & -1.219E-03,-1.097E-03,-1.007E-03,-9.605E-04,-8.446E-04,-6.844E-04, &
     & -5.439E-04,-3.774E-04,-3.453E-04,-3.051E-04,-2.558E-04,-5.225E-05, &
     &  3.788E-05, &
     & -1.081E-05,-1.077E-05,-1.073E-05,-1.061E-05,-1.037E-05,-1.019E-05, &
     & -1.013E-05,-1.243E-05,-9.884E-06,-8.594E-06,-7.784E-06,-7.084E-06, &
     & -4.178E-06,-8.468E-07, 1.688E-07, 1.031E-07, 1.356E-06, 3.126E-08, &
     & -2.031E-07/

		call apply_ckd(ni, k, u, f, p, t, tau, coefk, .false.)

      end subroutine ck_20
! *********************************************************************
!     co2 2605--2735 cm^{-1}
      subroutine ckco2_20(u,f,p,t,tau)
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		
		
		
		real :: k(nco2)
		integer :: i, m 
		
      f(1)=1.00000
      do  m=1,use_layers
          k(1)=0.0017949128*(1.0000-1.2487e-04*(t(m)-250.0)+ &
     &         3.4136e-06*(t(m)-250.0)**2)
          tau(1,m)=k(1)*u(m)
      end do

      end subroutine ckco2_20
! *********************************************************************
!     ch4 2605--2735 cm^{-1}
      subroutine ckch4_20(u,f,p,t,tau)
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		
		
		
		real :: k(nch4)
		integer :: i, m 
      f(1)=1.00000
      do  m=1,use_layers
          k(1)=14.524085*(1.0000+1.9377e-03*(t(m)-250.0)- &
     &         3.0119e-06*(t(m)-250.0)**2)
          tau(1,m)=k(1)*u(m)
      end do

      end subroutine ckch4_20

end module modis20
