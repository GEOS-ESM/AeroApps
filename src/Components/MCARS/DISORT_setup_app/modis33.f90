module modis33

	use CKD_common

	implicit none
	private
	
	integer, parameter :: ni = 4
	integer, parameter :: nco2 = 4
	integer, parameter :: no3 = 1
	
	public :: ch33, init_ch33
	
contains

	subroutine init_ch33(modis_channel)
	
		type(ckd_type), intent(inout) :: modis_channel
		integer :: all_weights
		
		all_weights = ni*nco2*no3 
		
		modis_channel%num_weights = all_weights
		allocate(modis_channel%weights(all_weights))
		allocate(modis_channel%taus(all_weights, max_layers))
	
	end subroutine init_ch33

      subroutine ch33(u0, p0, t0, ux, fac, u_co2, modis_channel)
!
!      This routine was developed for MODIS Channel 33: 740--760 cm^{-1}
!
!-----------------------------------------------------------------
!     Channel 33   13.2 - 13.5 microns
!-----------------------------------------------------------------
!  INPUTS:
!     u0 --> water vapor amount [g/cm2/km]
!     p0 --> pressure in atmospheres
!     t0 --> temp in K
!	  ux --> ozone (g/cm^2/km)
!     fac --> layer thickness in km
!-----------------------------------------------------------------

		real, dimension(:), intent(in) :: u0, p0, t0, ux, fac
		real, intent(in) :: u_co2
		type(ckd_type), intent(inout) :: modis_channel
		
		integer, parameter :: nlev = max_levels
		integer, parameter :: nlay = max_layers
		integer, parameter :: mlv = nlay

		real :: p(nlay), t(nlay),u(nlay), &
				 uco2(nlay), uo3(nlay)
				
		real :: tau(ni,nlay),tauco2(nco2,nlay), tauo3(no3, nlay), &
				f(ni), fco2(nco2), fo3(no3),  xlog, twin
		
		integer :: i,m,ia,ib, index
		
		real, parameter ::  omega=750.0

	  
      do  m=1,use_layers
        t(m)=(t0(m)+t0(m+1))/2.0
        xlog=(log10(p0(m))+log10(p0(m+1)))/2.0
        p(m)=10.0**xlog
        xlog=(log10(u0(m))+log10(u0(m+1)))/2.0
        u(m)=fac(m)*10.0**xlog
        xlog=(log10(ux(m))+log10(ux(m+1)))/2.0
        uo3(m)=fac(m)*10.0**xlog
        uco2(m)=p(m)*1.01325e+06/(1.3805e-16*t(m))*u_co2*fac(m)* &
     &          1.0e+05/6.023e+23*44.00995
	  enddo
	  
      call ck_33(u,f,p,t,tau)
      call ckco2_33(uco2,fco2,p,t,tauco2)
      call cko3_33(uo3,fo3,p,t,tauo3)
	  
      do  m=1,use_layers
        call window_33(omega,u(m),twin,p(m),t(m),fac(m))
         do ia=1,nco2 ! 4 k's for CO2
        do i=1,ni ! 4 k's for H2O
			index = (ia-1)*ni + i 
           modis_channel%taus(index,m)=tau(i,m)+twin+tauco2(ia,m)+tauo3(1,m)
 			if (m == 1) modis_channel%weights(index) = f(i)*fco2(ia)*fo3(1)
         end do
	  enddo
	end do

      end subroutine ch33

! *********************************************************************
!     h2o 740--760 cm^{-1}
      subroutine ck_33(u,f,p,t,tau)
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		

		real :: k(ni)
		real :: coefk(ni,3,num_pressures)
		integer :: i, jp, jt
      f(1)=0.594100
      f(2)=0.291514
      f(3)=0.091173
      f(4)=0.023213
      k(1)=0.00789657
      do i=2,ni
        k(i)=12.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/  &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00000,0.00000,0.00030,0.00507, &
     & 0.03979,0.12412,0.25886,0.45891,0.78143,1.30000, &
     & 2.11435, &
     & 0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00, &
     & 0.000e+00,0.000e+00,0.000e+00,0.000e+00,6.625e-06,1.805e-04, &
     & 1.327e-03,3.542e-03,7.119e-03,1.190e-02,2.050e-02,3.479e-02, &
     & 5.632e-02, &
     & 0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00, &
     & 0.000e+00,0.000e+00,0.000e+00,0.000e+00,3.063e-08,2.094e-06, &
     & 1.420e-05,3.306e-05,6.478e-05,1.011e-04,1.775e-04,3.086e-04, &
     & 4.960e-04/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00018,0.00163,0.01314,0.03975,0.08007, &
     & 0.13507,0.21349,0.33541,0.51274,0.80494,1.27505, &
     & 1.99251, &
     & 0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00, &
     & 0.000e+00,4.550e-06,5.176e-05,4.037e-04,1.157e-03,2.294e-03, &
     & 3.846e-03,6.093e-03,9.576e-03,1.421e-02,2.263e-02,3.570e-02, &
     & 5.480e-02, &
     & 0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00, &
     & 0.000e+00,3.000e-08,5.309e-07,3.994e-06,1.086e-05,2.138e-05, &
     & 3.582e-05,5.713e-05,8.959e-05,1.274e-04,2.079e-04,3.297e-04, &
     & 4.884e-04/
      data ( (coefk(3,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 0.00000,0.00000,0.00000,0.00000,0.00024,0.00190, &
     & 0.00727,0.01595,0.02871,0.04618,0.07037,0.10720, &
     & 0.16346,0.25089,0.38721,0.58023,0.85328,1.17060, &
     & 1.43171, &
     & 0.000e+00,0.000e+00,0.000e+00,0.000e+00,7.288e-06,7.151e-05, &
     & 2.349e-04,4.805e-04,8.375e-04,1.323e-03,2.006e-03,3.034e-03, &
     & 4.624e-03,7.050e-03,1.075e-02,1.547e-02,2.315e-02,3.211e-02, &
     & 4.006e-02, &
     & 0.000e+00,0.000e+00,0.000e+00,0.000e+00,6.156e-08,8.809e-07, &
     & 2.516e-06,4.800e-06,8.162e-06,1.264e-05,1.895e-05,2.806e-05, &
     & 4.289e-05,6.528e-05,9.775e-05,1.318e-04,2.010e-04,2.808e-04, &
     & 3.634e-04/
      data ( (coefk(4,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 1.62416,1.62188,1.61854,1.61604,1.61215,1.60553, &
     & 1.59477,1.57903,1.55709,1.52775,1.49206,1.45074, &
     & 1.40907,1.36327,1.30740,1.19614,1.07940,0.93661, &
     & 0.77570, &
     & 4.769e-02,4.766e-02,4.763e-02,4.763e-02,4.762e-02,4.751e-02, &
     & 4.739e-02,4.715e-02,4.678e-02,4.621e-02,4.541e-02,4.432e-02, &
     & 4.307e-02,4.173e-02,4.017e-02,3.601e-02,3.337e-02,2.983e-02, &
     & 2.508e-02, &
     & 4.869e-04,4.868e-04,4.872e-04,4.872e-04,4.867e-04,4.847e-04, &
     & 4.837e-04,4.804e-04,4.762e-04,4.704e-04,4.628e-04,4.526e-04, &
     & 4.399e-04,4.269e-04,4.123e-04,3.609e-04,3.447e-04,3.213e-04, &
     & 2.728e-04/
		call apply_ckd(ni, k, u, f, p, t, tau, coefk, .false.)

      end subroutine ck_33
! *********************************************************************
!     co2 740--760 cm^{-1}
      subroutine ckco2_33(u,f,p,t,tau)
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		
		
		real :: k(nco2)
		real :: coefk(nco2,3,num_pressures)
		integer :: i, jp, jt
      f(1)=0.196902
      f(2)=0.510758
      f(3)=0.241133
      f(4)=0.051207
      k(1)=0.110955843
      do i=2,nco2
        k(i)=8.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00003, &
     & 0.00012,0.00041,0.00166,0.00512,0.01283,0.06269, &
     & 0.12158,0.20155,0.33703,0.52899,0.80912,1.24067, &
     & 1.92344, &
     & 0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,1.125E-06, &
     & 4.250E-06,1.375E-05,5.113E-05,1.606E-04,4.025E-04,1.632E-03, &
     & 3.633E-03,5.959E-03,1.016E-02,1.620E-02,2.496E-02,3.751E-02, &
     & 5.636E-02, &
     & 0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,1.562E-08, &
     & 5.000E-08,1.500E-07,4.906E-07,1.666E-06,4.362E-06,1.202E-05, &
     & 3.342E-05,5.480E-05,9.455E-05,1.542E-04,2.430E-04,3.629E-04, &
     & 5.312E-04/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 0.00000,0.00000,0.00004,0.00013,0.00036,0.00099, &
     & 0.00243,0.00581,0.01338,0.02703,0.05062,0.08767, &
     & 0.14716,0.23299,0.36725,0.57104,0.84690,1.18440, &
     & 1.64723, &
     & 0.000E+00,0.000E+00,1.500E-06,5.625E-06,1.450E-05,3.525E-05, &
     & 7.838E-05,1.774E-04,3.860E-04,7.554E-04,1.369E-03,2.308E-03, &
     & 3.765E-03,5.862E-03,8.885E-03,1.327E-02,1.946E-02,2.708E-02, &
     & 3.698E-02, &
     & 0.000E+00,0.000E+00,1.875E-08,7.188E-08,1.813E-07,4.125E-07, &
     & 8.469E-07,1.841E-06,3.869E-06,7.222E-06,1.227E-05,1.965E-05, &
     & 3.057E-05,4.842E-05,7.013E-05,9.655E-05,1.343E-04,1.894E-04, &
     & 2.462E-04/
      data ( (coefk(3,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 0.00122,0.00203,0.00337,0.00565,0.00910,0.01429, &
     & 0.02222,0.03410,0.05139,0.07954,0.12627,0.19651, &
     & 0.29222,0.40711,0.54831,0.71327,0.89947,1.10362, &
     & 1.27980, &
     & 5.925E-05,8.312E-05,1.204E-04,1.795E-04,2.605E-04,3.796E-04, &
     & 5.730E-04,8.918E-04,1.396E-03,2.138E-03,3.287E-03,4.949E-03, &
     & 7.195E-03,9.909E-03,1.323E-02,1.704E-02,2.093E-02,2.527E-02, &
     & 2.903E-02, &
     & 8.500E-07,1.078E-06,1.416E-06,1.913E-06,2.531E-06,3.409E-06, &
     & 4.906E-06,7.825E-06,1.323E-05,1.990E-05,2.849E-05,3.993E-05, &
     & 5.484E-05,7.507E-05,9.953E-05,1.273E-04,1.516E-04,1.778E-04, &
     & 2.026E-04/
      data ( (coefk(4,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 1.39285,1.29074,1.18546,1.09340,1.01971,0.97005, &
     & 0.94916,0.95801,0.99974,1.07688,1.17504,1.25839, &
     & 1.30447,1.31094,1.27437,1.19537,1.07283,0.92242, &
     & 0.76295, &
     & 2.746E-02,2.512E-02,2.301E-02,2.117E-02,1.978E-02,1.893E-02, &
     & 1.886E-02,1.968E-02,2.135E-02,2.374E-02,2.629E-02,2.801E-02, &
     & 2.839E-02,2.805E-02,2.719E-02,2.529E-02,2.244E-02,1.926E-02, &
     & 1.587E-02, &
     & 1.556E-04,1.445E-04,1.358E-04,1.249E-04,1.161E-04,1.126E-04, &
     & 1.154E-04,1.285E-04,1.491E-04,1.705E-04,1.869E-04,1.926E-04, &
     & 1.864E-04,1.798E-04,1.761E-04,1.613E-04,1.408E-04,1.228E-04, &
     & 1.002E-04/
		call apply_ckd(nco2, k, u, f, p, t, tau, coefk, .false.)
      end subroutine ckco2_33
! *********************************************************************
!     o3 740--760 cm^{-1}
      subroutine cko3_33(u,f,p,t,tau)
	
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		
		
		real :: k(no3)
		integer :: i, m 
      f(1)=1.00000
      do  m=1,use_layers
          k(1)=38.51857*(1.0000+1.0514e-03*(t(m)-250.0)- &
     &         2.5531e-06*(t(m)-250.0)**2)
          tau(1,m)=k(1)*u(m)
      end do
	  end subroutine cko3_33
	
!**********************************************************************
      subroutine window_33(omega,u,twin,pres,temp,fac)
	  
		real, intent(in) :: omega, u, pres, temp, fac
		real, intent(inout) :: twin
		
		real :: r, amnt, a1, a2, ph2o, tau 
	  
      r=8.314d+07
!     amnt=u/6.023e+23*2.687e+19*18.01534
      amnt=u
!     write(6,*)u,amnt
!     a1=4.18+5577.8*exp(-.00787*omega)
      a1=4.18+7815.6*exp(-.0083*omega)
      a1=a1*exp(1800.*(1./temp-1./296.))
      a2=0.0008*a1
      ph2o=amnt/(fac*1.0d+05)*r*temp/18.01534/1.01325d+06
      tau=amnt*(a1*ph2o+a2*(pres-ph2o))
!     tau=amnt*(a1*ph2o)
      if(tau > 88)then
		twin=88.0
      else
		twin=tau
      end if

	end subroutine window_33
	
end module modis33