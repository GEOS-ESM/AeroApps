module modis5a 

	use CKD_common

	implicit none
	private
	
	integer, parameter :: ni = 1
	integer, parameter :: no2 = 3
	
	public :: ch5a, init_ch5a

contains

	subroutine init_ch5a(modis_channel)
	
		type(ckd_type), intent(inout) :: modis_channel
		integer :: all_weights
		
		all_weights = ni*no2
		
		modis_channel%num_weights = all_weights
		allocate(modis_channel%weights(all_weights))
		allocate(modis_channel%taus(all_weights, max_layers))
	
	end subroutine init_ch5a

      subroutine ch5a(u0, p0, t0, ux, fac, modis_channel)
!
!      This routine was developed for MODIS Channel 5 
!                        7890--8050 cm^{-1}
!
!      Note that band 5 is split into two parts
!-----------------------------------------------------------------
!     Channel 5a   1.242-1.267 microns
!-----------------------------------------------------------------
!  INPUTS:
!     u0 --> water vapor amount [g/cm2/km]
!     p0 --> pressure in atmospheres
!     t0 --> temp in K
!     ux --> ozone [g/cm2/km]
!     dz --> layer thickness in km
!-----------------------------------------------------------------

		real, dimension(:), intent(in) :: u0, p0, t0, ux, fac
		type(ckd_type), intent(inout) :: modis_channel
		
		integer, parameter :: nlev = max_levels
		integer, parameter :: nlay = max_layers
		integer, parameter :: mlv = nlay

		real :: p(nlay), t(nlay),u(nlay), &
				 uo2(nlay)
				
		real :: tau(ni,nlay),tauo2(no2,nlay), &
				f(ni), fo2(no2)
		
		integer :: i,m,ia, index


      do  m=1,use_layers
        t(m)=(t0(m)+t0(m+1))/2.0
        p(m)=(p0(m)+p0(m+1))/2.0
        u(m)=fac(m)*(u0(m)+u0(m+1))/2.0
        uo2(m)=p(m)*1.01325e+06/(1.3805e-16*t(m))*0.20948*fac(m)* &
            1.0e+05/6.023e+23*31.9988
      enddo

      call ck_5a(u,f,p,t,tau)
      call cko2_5a(uo2,fo2,p,t,tauo2)
	  

      do  m=1,use_layers
          do ia=1,no2 ! 3 k's for o2
        do i=1,ni ! 1 k for h2o
			index = (ia-1)*ni + i
            modis_channel%taus(index,m)=tau(i,m)+tauo2(ia,m)
 			if (m == 1) modis_channel%weights(index) = f(i)*fo2(ia)
          end do
        end do
      enddo
	  

      end subroutine ch5a 


! *********************************************************************
      subroutine ck_5a(u,f,p,t,tau)
	  
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		
		
		integer :: m 

		real :: k(ni)

      f(1)=1.000000
      do  m=1,use_layers
        k(1)=0.888969*4.67179206e-04*(1.0000+1.3163e-02*(t(m)-250.0)+ &
            6.3792e-05*(t(m)-250.0)**2)
        tau(1,m)=k(1)*u(m)
      end do
      return
      end subroutine ck_5a
!**********************************************************************
      subroutine cko2_5a(u,f,p,t,tau)
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		
		
		real :: k(no2)
		real :: coefk(no2,3,num_pressures)
		integer :: i, jp, jt

      f(1)=0.889414
      f(2)=0.089339
      f(3)=0.021247
      k(1)=0.71*0.0000073895879
      do i=2,no2
        k(i)=21.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00020,0.00177,0.01554,0.04868, &
     & 0.10224,0.18791,0.31384,0.50837,0.81375,1.30088, &
     &  2.05278, &
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, &
     &  0.000E+00, 0.000E+00,-7.375E-07,-1.875E-06, 6.500E-06, 5.125E-06, &
     &  0.000E+00,-1.500E-06, 7.175E-05, 2.850E-05, 2.988E-05,-3.751E-06, &
     & -1.575E-05, &
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, &
     &  0.000E+00, 0.000E+00,-2.812E-09,-2.187E-08,-1.438E-07,-1.719E-07, &
     & -1.750E-07,-3.875E-07,-4.563E-07,-6.812E-07,-9.219E-07,-1.881E-06, &
     & -2.550E-06/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 0.00314,0.00346,0.00393,0.00457,0.00533,0.00613, &
     & 0.01133,0.02110,0.03454,0.05400,0.07995,0.11570, &
     & 0.17012,0.25380,0.38214,0.57701,0.85573,1.22220, &
     &  1.63284, &
     &  4.300E-05, 4.388E-05, 4.638E-05, 5.212E-05, 6.125E-05, 6.962E-05, &
     &  4.025E-05, 4.100E-05, 3.838E-05, 2.825E-05, 1.287E-05, 2.125E-06, &
     & -8.375E-06,-2.312E-05,-7.013E-05,-1.551E-04,-2.994E-04,-5.380E-04, &
     & -8.141E-04, &
     &  1.375E-07, 1.156E-07, 7.187E-08, 3.125E-09,-5.000E-08, 5.312E-08, &
     & -6.250E-08,-1.250E-07,-1.406E-07,-1.438E-07,-1.156E-07, 9.376E-09, &
     &  1.531E-07, 2.781E-07, 7.281E-07, 8.406E-07, 1.516E-06, 2.437E-06, &
     &  8.844E-07/
      data ( (coefk(3,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 1.33872,1.33844,1.33799,1.33595,1.33497,1.33354, &
     & 1.33155,1.32865,1.32327,1.31737,1.30873,1.29448, &
     & 1.27465,1.24794,1.20970,1.15580,1.07737,0.96672, &
     &  0.81929, &
     & -5.364E-04,-5.184E-04,-5.176E-04,-5.170E-04,-5.163E-04,-5.320E-04, &
     & -5.175E-04,-5.176E-04,-4.999E-04,-4.975E-04,-4.771E-04,-4.750E-04, &
     & -4.570E-04,-4.561E-04,-4.578E-04,-4.615E-04,-4.566E-04,-4.219E-04, &
     & -3.957E-04, &
     &  1.084E-06, 6.469E-07, 6.593E-07, 1.525E-06, 1.525E-06, 1.094E-06, &
     &  7.125E-07, 7.281E-07, 1.159E-06, 1.175E-06, 7.406E-07, 7.125E-07, &
     &  1.100E-06, 1.078E-06, 1.419E-06, 9.750E-07, 1.241E-06, 7.531E-07, &
     &  1.350E-06/

		call apply_ckd(no2, k, u, f, p, t, tau, coefk, .false.)
      end subroutine cko2_5a
	  
end module modis5a
