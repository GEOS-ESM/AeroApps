module modis19

	use CKD_common

	implicit none
	private
	
	integer, parameter :: ni = 4
	
	public :: ch19, init_ch19

contains

	subroutine init_ch19(modis_channel)
	
		type(ckd_type), intent(inout) :: modis_channel
		integer :: all_weights
		
		all_weights = ni
		
		modis_channel%num_weights = all_weights
		allocate(modis_channel%weights(all_weights))
		allocate(modis_channel%taus(all_weights, max_layers))
	
	end subroutine init_ch19

      subroutine ch19(u0, p0, t0, ux, fac, modis_channel)
!
!      This routine was developed for MODIS Channel 19 
!                        10430--10960 cm^{-1}
!
!-----------------------------------------------------------------
!     Channel 19   0.912-0.959 microns
!
!     revised by David Kratz on 7/26/2007
!     also accounts for the continuum (routine modis19x.f90 does not)
!
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
				uo2(nlay),uo3(nlay)
				
		real :: tau(ni,nlay), f(ni), twinx(ni), twin 
		
		integer :: m, i
		
      do  m=1,use_layers
        t(m)=(t0(m)+t0(m+1))/2.0
        p(m)=(p0(m)+p0(m+1))/2.0
        u(m)=fac(m)*(u0(m)+u0(m+1))/2.0
      enddo
	  
      call ck_19(u,f,p,t,tau)
	  
      do  m=1,use_layers
        call window_19(u(m),twin,p(m),t(m),fac(m))
        twinx(1)=0.1634974*twin*1.643519
        twinx(2)=0.3509185*twin*1.643519
        twinx(3)=0.3403536*twin*1.643519
        twinx(4)=0.1452305*twin*1.643519

         do i=1,ni ! 4 k's for h2o
          modis_channel%taus(i,m)=tau(i,m)+twinx(i)
 			if (m == 1) modis_channel%weights(i) = f(i)
         end do
      enddo
 
      end subroutine ch19

! *********************************************************************
      subroutine ck_19(u,f,p,t,tau)
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		
		
		real :: k(ni)
		real :: coefk(ni,3,num_pressures)
		integer :: i, jp, jt

      f(1)=0.537422
      f(2)=0.347788
      f(3)=0.101705
      f(4)=0.013085
      k(1)=0.040109846
      do i=2,ni
        k(i)=11.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 0.00125,0.00185,0.00273,0.00405,0.00601,0.00896, &
     & 0.01341,0.02018,0.03045,0.04619,0.06981,0.10605, &
     & 0.16153,0.24574,0.37054,0.55524,0.82447,1.20557, &
     &  1.70015, &
     &  5.125E-06, 7.000E-06, 9.250E-06, 1.262E-05, 1.737E-05, 2.362E-05, &
     &  3.250E-05, 4.562E-05, 6.550E-05, 9.500E-05, 1.390E-04, 2.065E-04, &
     &  3.081E-04, 4.695E-04, 7.158E-04, 1.017E-03, 1.432E-03, 2.051E-03, &
     &  2.800E-03, &
     &  3.125E-09,-6.250E-09,-1.250E-08,-2.813E-08,-4.063E-08,-6.563E-08, &
     & -1.063E-07,-1.844E-07,-2.750E-07,-4.500E-07,-7.125E-07,-1.019E-06, &
     & -1.491E-06,-2.012E-06,-2.431E-06,-2.937E-06,-4.234E-06,-5.688E-06, &
     & -8.038E-06/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 0.02466,0.02586,0.02771,0.03036,0.03419,0.03973, &
     & 0.04756,0.05880,0.07466,0.09742,0.12914,0.17471, &
     & 0.23979,0.33231,0.45988,0.63405,0.86191,1.14948, &
     &  1.47436, &
     &  3.060E-04, 3.040E-04, 3.026E-04, 3.030E-04, 3.066E-04, 3.120E-04, &
     &  3.224E-04, 3.389E-04, 3.619E-04, 3.895E-04, 4.204E-04, 4.584E-04, &
     &  5.066E-04, 5.650E-04, 6.431E-04, 6.987E-04, 7.815E-04, 9.155E-04, &
     &  1.079E-03, &
     &  1.156E-06, 1.131E-06, 1.066E-06, 1.019E-06, 9.781E-07, 9.000E-07, &
     &  8.656E-07, 7.469E-07, 7.156E-07, 5.938E-07, 3.656E-07, 2.344E-07, &
     &  9.688E-08,-1.875E-07,-4.219E-07,-1.025E-06,-1.862E-06,-3.950E-06, &
     & -4.584E-06/
      data ( (coefk(3,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 0.29496,0.29539,0.29641,0.29760,0.29960,0.30325, &
     & 0.30870,0.31819,0.33261,0.35556,0.38639,0.43257, &
     & 0.49728,0.58131,0.68510,0.80666,0.93435,1.06129, &
     &  1.16180, &
     &  1.732E-03, 1.730E-03, 1.721E-03, 1.712E-03, 1.702E-03, 1.675E-03, &
     &  1.641E-03, 1.588E-03, 1.513E-03, 1.405E-03, 1.256E-03, 1.092E-03, &
     &  9.172E-04, 7.450E-04, 5.989E-04, 4.624E-04, 3.570E-04, 1.972E-04, &
     &  4.012E-05, &
     & -1.287E-06,-1.269E-06,-1.356E-06,-1.316E-06,-1.147E-06,-1.206E-06, &
     & -9.031E-07,-9.906E-07,-6.344E-07,-7.406E-07,-7.344E-07,-5.875E-07, &
     & -1.125E-06,-1.019E-06,-1.459E-06,-2.522E-06,-2.181E-06,-1.275E-06, &
     & -1.416E-06/
      data ( (coefk(4,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 1.78024,1.77949,1.78023,1.77836,1.77539,1.77265, &
     & 1.76533,1.75779,1.74214,1.72111,1.67675,1.62309, &
     & 1.55539,1.47017,1.36110,1.22980,1.07615,0.91625, &
     &  0.75218, &
     & -1.576E-03,-1.573E-03,-1.594E-03,-1.586E-03,-1.551E-03,-1.556E-03, &
     & -1.528E-03,-1.487E-03,-1.429E-03,-1.370E-03,-1.275E-03,-1.171E-03, &
     & -1.064E-03,-9.365E-04,-8.485E-04,-7.939E-04,-7.222E-04,-6.522E-04, &
     & -5.941E-04, &
     &  1.056E-06, 1.047E-06, 4.562E-07, 4.312E-07, 9.875E-07, 3.687E-07, &
     &  1.506E-06, 2.250E-07, 1.319E-06, 7.688E-07, 1.501E-07, 7.750E-07, &
     &  6.250E-07, 6.250E-07, 1.694E-06, 2.278E-06, 2.963E-06, 2.788E-06, &
     &  1.797E-06/

		call apply_ckd(ni, k, u, f, p, t, tau, coefk, .false.)

      end subroutine ck_19
! *********************************************************************
      subroutine window_19 (amnt,twin,patm,temp,fac)
!     Parameterization of CKD_2.4 continuum over 10430 to 10960 cm-1 band
!     Fit inaccurate for amnt < 0.01 g/cm*2;
!          however optical depths < 0.001.
!     INPUT:
!     amnt = h2O layer amount (g/cm**2)
!     patm= pressure (atm)
!     temp = temperature (k)
!     OUTPUT:
!     twin = parameterized CKD_2.1 optical depth over band
!

	  real, intent(in) :: amnt, patm, temp, fac
	  real, intent(inout) :: twin

      integer, parameter :: ncoef=7

      real ::  aa(ncoef), ph2o, tau_log

      data aa/-2.119,1.002,-5.192e-03,0.9057,4.698,-1.894e-03,5.871e-02/

      ph2o = amnt*(8.314e+07*temp)/(fac*1.0e+05*18.01534*1.01325e+06)

          tau_log   = aa(1)              + &
                      aa(2) * log(amnt)  + &
                      aa(3) * temp       + &
                      aa(4) * patm       + &
                      aa(5) * ph2o       + &
                      aa(6) * amnt       + &
                      aa(7) * log(ph2o)

        twin = exp ( tau_log )
		

        end subroutine window_19

end module modis19
