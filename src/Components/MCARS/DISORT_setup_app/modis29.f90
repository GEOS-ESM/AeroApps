module modis29

	use CKD_common

	implicit none
	private
	
	integer, parameter :: ni = 5
	integer, parameter :: no3 = 1
	integer, parameter :: nch4 = 1
	integer, parameter :: nn2o = 2
	
	public :: ch29, init_ch29
	
contains

	subroutine init_ch29(modis_channel)
	
		type(ckd_type), intent(inout) :: modis_channel
		integer :: all_weights
		
		all_weights = ni*nch4*no3*nn2o 
		
		modis_channel%num_weights = all_weights
		allocate(modis_channel%weights(all_weights))
		allocate(modis_channel%taus(all_weights, max_layers))
	
	end subroutine init_ch29

      subroutine ch29(u0, p0, t0, ux, fac,modis_channel)
!
!      This routine was developed for MODIS Channel 29: 1140-1200 cm^{-1}
!
!-----------------------------------------------------------------
!     Channel 29   8.33 - 8.77 microns
!-----------------------------------------------------------------
!  INPUTS:
!     u0 --> water vapor amount [g/cm2/km]
!     p0 --> pressure in atmospheres
!     t0 --> temp in K
!	  ux --> ozone (g/cm^2/km)
!     fac --> layer thickness in km
!-----------------------------------------------------------------
		real, dimension(:), intent(in) :: u0, p0, t0, ux, fac
		type(ckd_type), intent(inout) :: modis_channel
		
		integer, parameter :: nlev = max_levels
		integer, parameter :: nlay = max_layers
		integer, parameter :: mlv = nlay

		real :: p(nlay), t(nlay),u(nlay), uo3(nlay), un2o(nlay), xlog, &
				  uch4(nlay)
				
		real :: tau(ni,nlay),tauch4(nch4, nlay), tauo3(no3, nlay), taun2o(nn2o, nlay), &
				 f(ni), fo3(no3), fn2o(nn2o), fch4(nch4), twin
		
		integer :: i,m, ic, id, index
	  
      do m=1,nlay
        t(m)=(t0(m)+t0(m+1))/2.0
        xlog=(log10(p0(m))+log10(p0(m+1)))/2.0
        p(m)=10.0**xlog
        xlog=(log10(u0(m))+log10(u0(m+1)))/2.0
        u(m)=fac(m)*10.0**xlog
        xlog=(log10(ux(m))+log10(ux(m+1)))/2.0
        uo3(m)=fac(m)*10.0**xlog
        uch4(m)=p(m)*1.01325e+06/(1.3805e-16*t(m))*1.75e-06*fac(m)* &
     &          1.0e+05/6.023e+23*16.04303
        un2o(m)=p(m)*1.01325e+06/(1.3805e-16*t(m))*3.10e-07*fac(m)* &
     &          1.0e+05/6.023e+23*44.0128
	  enddo
 
      call ck_29(u,f,p,t,tau)
      call cko3_29(uo3,fo3,p,t,tauo3)
      call ckch4_29(uch4,fch4,p,t,tauch4)
      call ckn2o_29(un2o,fn2o,p,t,taun2o)
	  
      do m=1,nlay
        call window_29(u(m),twin,p(m),t(m),fac(m))
           do id=1,nn2o   ! two k's for n2o
        do i=1,ni         ! five k's for water vapor
		   index = (id-1)*ni + i
			modis_channel%taus(index,m)=tau(i,m)+twin+tauo3(1,m) &
     &          +tauch4(1,m)+taun2o(id,m)
			if (m == 1) modis_channel%weights(index) = f(i)*fo3(1)*fch4(1)*fn2o(id)
           end do
        end do
      enddo
	  

	end subroutine ch29

! *********************************************************************
!     h2o 1140-1200 cm^{-1}
      subroutine ck_29(u,f,p,t,tau)
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		
		
		integer, parameter :: mlv = max_layers
		integer :: m 

		real :: k(ni)
		real :: coefk(ni,3,num_pressures)
		integer :: i, jp, jt
      f(1)=0.624325
      f(2)=0.253393
      f(3)=0.089834
      f(4)=0.026964
      f(5)=0.005484
      k(1)=0.0072755
      do i=2,ni
        k(i)=8.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00006,0.00065,0.00265,0.00730,0.01769, &
     & 0.04487,0.11014,0.24011,0.44539,0.77830,1.26324, &
     & 1.93516, &
     & 0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00, &
     & 0.000e+00,1.275e-06,1.470e-05,5.463e-05,1.524e-04,3.462e-04, &
     & 8.261e-04,1.945e-03,4.299e-03,7.952e-03,1.399e-02,2.282e-02, &
     & 3.441e-02, &
     & 0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00, &
     & 0.000e+00,5.625e-09,1.125e-07,3.906e-07,1.178e-06,2.606e-06, &
     & 5.284e-06,1.178e-05,2.511e-05,4.654e-05,7.976e-05,1.279e-04, &
     & 1.867e-04/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 0.00000,0.00000,0.00000,0.00000,0.00001,0.00021, &
     & 0.00100,0.00301,0.00822,0.02028,0.04442,0.09274, &
     & 0.16324,0.26366,0.39942,0.58242,0.84505,1.18790, &
     & 1.68415, &
     & 0.000e+00,0.000e+00,0.000e+00,0.000e+00,2.125e-07,4.725e-06, &
     & 2.071e-05,6.063e-05,1.574e-04,3.879e-04,8.413e-04,1.719e-03, &
     & 2.968e-03,4.783e-03,7.230e-03,1.045e-02,1.519e-02,2.140e-02, &
     & 2.962e-02, &
     & 0.000e+00,0.000e+00,0.000e+00,0.000e+00,3.125e-10,3.563e-08, &
     & 1.459e-07,4.406e-07,1.066e-06,2.516e-06,5.706e-06,1.075e-05, &
     & 1.748e-05,2.754e-05,4.141e-05,5.808e-05,8.307e-05,1.204e-04, &
     & 1.526e-04/
      data ( (coefk(3,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 0.00010,0.00012,0.00032,0.00085,0.00193,0.00458, &
     & 0.01031,0.02119,0.04015,0.06639,0.10163,0.14948, &
     & 0.21695,0.30768,0.42844,0.59514,0.84740,1.17004, &
     & 1.48726, &
     & 6.413e-06,7.337e-06,1.126e-05,2.530e-05,5.484e-05,1.101e-04, &
     & 2.316e-04,4.420e-04,7.960e-04,1.290e-03,1.956e-03,2.831e-03, &
     & 4.046e-03,5.697e-03,7.895e-03,1.074e-02,1.516e-02,2.099e-02, &
     & 2.736e-02, &
     & 1.053e-07,1.197e-07,1.322e-07,2.650e-07,5.916e-07,1.016e-06, &
     & 1.897e-06,3.206e-06,5.294e-06,8.391e-06,1.266e-05,1.748e-05, &
     & 2.388e-05,3.295e-05,4.598e-05,5.965e-05,7.798e-05,1.001e-04, &
     & 1.449e-04/
      data ( (coefk(4,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 0.02081,0.02124,0.02440,0.03128,0.03987,0.05173, &
     & 0.06697,0.08747,0.11456,0.15100,0.20218,0.27745, &
     & 0.38290,0.51788,0.65919,0.78981,0.93236,1.06208, &
     & 1.13645, &
     & 6.691e-04,6.936e-04,7.050e-04,8.089e-04,9.665e-04,1.168e-03, &
     & 1.433e-03,1.786e-03,2.247e-03,2.882e-03,3.805e-03,5.040e-03, &
     & 6.803e-03,9.243e-03,1.208e-02,1.451e-02,1.720e-02,1.948e-02, &
     & 2.062e-02, &
     & 7.641e-06,8.097e-06,7.569e-06,7.547e-06,8.419e-06,9.331e-06, &
     & 1.065e-05,1.229e-05,1.459e-05,1.789e-05,2.317e-05,2.783e-05, &
     & 3.367e-05,4.389e-05,6.509e-05,7.864e-05,9.405e-05,1.047e-04, &
     & 9.958e-05/
      data ( (coefk(5,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 1.76643,1.76732,1.77176,1.77534,1.78324,1.79220, &
     & 1.80674,1.82530,1.84127,1.84779,1.83103,1.77690, &
     & 1.69182,1.57332,1.43992,1.27381,1.10369,0.89454, &
     & 0.69900, &
     & 3.494e-02,3.492e-02,3.490e-02,3.487e-02,3.489e-02,3.491e-02, &
     & 3.488e-02,3.491e-02,3.484e-02,3.469e-02,3.423e-02,3.325e-02, &
     & 3.176e-02,2.964e-02,2.708e-02,2.392e-02,2.113e-02,1.741e-02, &
     & 1.365e-02, &
     & 1.993e-04,1.990e-04,1.970e-04,1.965e-04,1.955e-04,1.953e-04, &
     & 1.920e-04,1.894e-04,1.848e-04,1.809e-04,1.761e-04,1.714e-04, &
     & 1.643e-04,1.556e-04,1.381e-04,1.207e-04,1.094e-04,9.372e-05, &
     & 7.355e-05/
	 
		call apply_ckd(ni, k, u, f, p, t, tau, coefk, .false.)

      end subroutine ck_29
!  *********************************************************************
      subroutine window_29 (amnt,twin,patm,temp,fac)
!      Parameterization of CKD_2.4 continuum over 1140 to 1200 cm-1 band
!      Fit inaccurate for amnt < 0.01 g/cm*2;
!           however optical depths < 0.001.
!      INPUT:
!      amnt = h2O layer amount (g/cm**2)
!      patm= pressure (atm)
!      temp = temperature (k)
!      OUTPUT:
!      twin = parameterized CKD_2.4 optical depth over band
!
	  real, intent(in) :: amnt, patm, temp, fac
	  real, intent(inout) :: twin

      integer, parameter :: ncoef=7

      real ::  aa(ncoef), ph2o, tau_log
!
      data aa/7.487,1.001,-0.02105,0.02266,-0.4916,-4.798e-04,0.9770/
!
      ph2o = amnt*(8.314e+07*temp)/(fac*1.0e+05*18.01534*1.01325e+06)
!
          tau_log   = aa(1)              + &
     &                 aa(2) * log(amnt)  + &
     &                 aa(3) * temp       + &
     &                 aa(4) * patm       + &
     &                 aa(5) * ph2o       + &
     &                 aa(6) * amnt       + &
     &                 aa(7) * log(ph2o)

        twin = exp ( tau_log )
!

        end subroutine window_29
!  *********************************************************************
!      o3 1140--1200 cm^{-1}
      subroutine cko3_29(u,f,p,t,tau)
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		
		
		integer, parameter :: mlv = max_layers
		integer :: m 

		real :: k(no3)

      f(1)=1.00000
      do m=1,mlv
          k(1)=15.60821*(1.0000+4.6204e-03*(t(m)-250.0)+ &
     &          8.7844e-06*(t(m)-250.0)**2)
          tau(1,m)=k(1)*u(m)
      end do
	 end subroutine cko3_29
!  *********************************************************************
!      ch4 1140--1200 cm^{-1}
      subroutine ckch4_29(u,f,p,t,tau)
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		
		
		integer, parameter :: mlv = max_layers
		integer :: m 

		real :: k(nch4)
      f(1)=1.00000
      do m=1,mlv
          k(1)=3.44205*(1.0000+1.8357e-02*(t(m)-250.0)+ &
     &          1.4238e-04*(t(m)-250.0)**2)
          tau(1,m)=k(1)*u(m)
      end do
	end subroutine ckch4_29
!  *********************************************************************
!      n2o 1140--1200 cm^{-1}
      subroutine ckn2o_29(u,f,p,t,tau)
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		
		
		integer, parameter :: mlv = max_layers
		integer :: m 

		real :: k(nn2o)
		real :: coefk(nn2o,3,num_pressures)
		integer :: i, jp, jt
      f(1)=0.819050
      f(2)=0.180950
      k(1)=24.4321126
      do i=2,nn2o
        k(i)=11.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 0.00000,0.00000,0.00010,0.00049,0.00167,0.00505, &
     & 0.01137,0.02124,0.03650,0.05920,0.09294,0.14144, &
     & 0.20923,0.30905,0.45202,0.63756,0.86176,1.17401, &
     & 1.53221, &
     & 0.000e+00,0.000e+00,4.625e-07,2.713e-06,1.113e-05,2.500e-05, &
     & 4.000e-05,6.125e-05,9.650e-05,1.523e-04,2.298e-04,3.539e-04, &
     & 4.795e-04,6.371e-04,9.179e-04,1.314e-03,1.267e-03,1.004e-03, &
     & 1.089e-03, &
     & 0.000e+00,0.000e+00,-3.125e-10,-7.812e-09,-3.125e-09,-1.250e-08, &
     & -2.910e-13,4.375e-08,6.250e-08,1.250e-07,1.812e-07,4.344e-07, &
     & 1.725e-06,1.959e-06,9.281e-07,2.384e-06,1.004e-05,8.034e-06, &
     & 6.109e-06/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 1.52647,1.52490,1.52238,1.51850,1.51379,1.50365, &
     & 1.48864,1.46990,1.44604,1.41986,1.39347,1.36856, &
     & 1.33966,1.29786,1.23922,1.16373,1.06093,0.94282, &
     & 0.79944, &
     & 6.396e-04,6.494e-04,6.670e-04,6.926e-04,7.307e-04,7.821e-04, &
     & 8.535e-04,9.373e-04,1.033e-03,1.082e-03,1.056e-03,9.970e-04, &
     & 9.186e-04,8.666e-04,7.456e-04,5.936e-04,5.576e-04,7.041e-04, &
     & 7.050e-04, &
     & 5.103e-06,4.909e-06,4.875e-06,4.672e-06,4.537e-06,4.316e-06, &
     & 4.006e-06,3.675e-06,3.947e-06,3.956e-06,4.153e-06,3.431e-06, &
     & 3.003e-06,3.334e-06,3.466e-06,2.672e-06,-4.719e-07,-2.406e-07, &
     & 9.312e-07/
		call apply_ckd(nn2o, k, u, f, p, t, tau, coefk, .false.)

      end subroutine ckn2o_29
	  
end module modis29
