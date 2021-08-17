module modis18

	use CKD_common

	implicit none
	private
	
	integer, parameter :: ni = 4
	
	public :: ch18, init_ch18
	
contains

	subroutine init_ch18(modis_channel)
	
		type(ckd_type), intent(inout) :: modis_channel
		integer :: all_weights
		
		all_weights = ni
		
		modis_channel%num_weights = all_weights
		allocate(modis_channel%weights(all_weights))
		allocate(modis_channel%taus(all_weights, max_layers))
	
	end subroutine init_ch18

      subroutine ch18(u0, p0, t0, ux, fac, modis_channel)
!
!      This routine was developed for MODIS Channel 18 
!                        10625--10740 cm^{-1}
!
!-----------------------------------------------------------------
!     Channel 18   0.931-0.941 microns
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

		real :: p(nlay), t(nlay),u(nlay)
				
		real :: tau(ni,nlay), f(ni), twin, twinx(ni)
		
		integer :: i,m
	  
      do  m=1,use_layers
        t(m)=(t0(m)+t0(m+1))/2.0
        p(m)=(p0(m)+p0(m+1))/2.0
        u(m)=fac(m)*(u0(m)+u0(m+1))/2.0
      enddo
	  
      call ck_18(u,f,p,t,tau)


      do  m=1,use_layers
		call window18(u(m),twin,p(m),t(m),fac(m))
        twinx(1)=0.3171613*twin*1.1551456
        twinx(2)=0.3386449*twin*1.1551456
        twinx(3)=0.2570058*twin*1.1551456
        twinx(4)=0.0871880*twin*1.1551456
        do i=1,ni ! 4 k's for h2o
          modis_channel%taus(i,m)=tau(i,m) + twinx(i)
 			if (m == 1) modis_channel%weights(i) = f(i)
        end do
      enddo
 

      end subroutine ch18

! *********************************************************************
      subroutine ck_18(u,f,p,t,tau)
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		

		real :: k(ni)
		real :: coefk(ni,3,num_pressures)
		integer :: i, jp, jt
      f(1)=0.564712
      f(2)=0.301482
      f(3)=0.114401
      f(4)=0.019405
      k(1)=0.11427365
      do i=2,ni
        k(i)=8.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     0.00178,0.00247,0.00347,0.00492,0.00705,0.01021,&
     0.01489,0.02191,0.03250,0.04846,0.07267,0.10927,&
     0.16438,0.24712,0.37051,0.55380,0.82289,1.21665,&
     1.74535,&
     9.875E-06, 1.212E-05, 1.500E-05, 1.888E-05, 2.388E-05, 3.075E-05,&
     3.963E-05, 5.213E-05, 6.950E-05, 9.300E-05, 1.253E-04, 1.663E-04,&
     2.190E-04, 2.856E-04, 3.660E-04, 4.249E-04, 4.891E-04, 5.778E-04,&
     9.047E-04,&
     9.375E-09, 3.125E-09,-6.250E-09,-1.563E-08,-3.437E-08,-6.250E-08,&
     -9.063E-08,-1.406E-07,-2.125E-07,-2.750E-07,-3.875E-07,-5.250E-07,&
     -6.187E-07,-7.469E-07,-1.150E-06,-1.772E-06,-2.659E-06,-5.100E-06,&
     -9.363E-06/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/&
     0.02158,0.02263,0.02425,0.02664,0.03007,0.03505,&
     0.04213,0.05227,0.06684,0.08765,0.11777,0.16131,&
     0.22464,0.31626,0.44584,0.62431,0.86183,1.15748,&
     1.49822,&
     2.581E-04, 2.546E-04, 2.511E-04, 2.486E-04, 2.476E-04, 2.486E-04,&
     2.507E-04, 2.555E-04, 2.618E-04, 2.651E-04, 2.680E-04, 2.665E-04,&
     2.574E-04, 2.533E-04, 2.482E-04, 1.757E-04, 1.439E-04, 6.375E-06,&
     -1.844E-04,&
     1.034E-06, 1.022E-06, 9.844E-07, 9.469E-07, 9.344E-07, 9.094E-07,&
     8.812E-07, 8.562E-07, 7.750E-07, 7.406E-07, 5.625E-07, 4.875E-07,&
     3.844E-07,-6.250E-08,-4.875E-07,-3.813E-07,-2.291E-06, 5.406E-07,&
     -7.156E-07/
      data ( (coefk(3,jt,jp), jp = 1, 19), jt = 1, 3)/&
     0.25373,0.25418,0.25492,0.25616,0.25824,0.26197,&
     0.26755,0.27693,0.29190,0.31516,0.35073,0.40191,&
     0.47139,0.56354,0.67810,0.80795,0.93901,1.05671,&
     1.11901,&
     1.636E-03, 1.633E-03, 1.627E-03, 1.618E-03, 1.603E-03, 1.581E-03,&
     1.542E-03, 1.490E-03, 1.416E-03, 1.298E-03, 1.142E-03, 9.206E-04,&
     6.606E-04, 4.164E-04, 1.916E-04,-6.375E-05,-3.354E-04,-5.551E-04,&
     -6.223E-04,&
     -8.875E-07,-8.688E-07,-8.344E-07,-7.969E-07,-7.531E-07,-6.844E-07,&
     -4.719E-07,-3.344E-07,-2.250E-07,-1.251E-08,-2.781E-07,-5.656E-07,&
     2.344E-07, 3.844E-07,-1.594E-07, 3.687E-07, 9.469E-07,-5.656E-07,&
     2.238E-06/
      data ( (coefk(4,jt,jp), jp = 1, 19), jt = 1, 3)/&
     1.85208,1.85129,1.85005,1.84809,1.84499,1.84214,&
     1.83454,1.82466,1.81037,1.78619,1.75194,1.70073,&
     1.63003,1.53482,1.41271,1.26257,1.08708,0.91236,&
     0.74049,&
     -2.769E-03,-2.765E-03,-2.759E-03,-2.751E-03,-2.737E-03,-2.720E-03,&
     -2.715E-03,-2.671E-03,-2.580E-03,-2.509E-03,-2.375E-03,-2.242E-03,&
     -2.048E-03,-1.845E-03,-1.671E-03,-1.481E-03,-1.257E-03,-1.116E-03,&
     -8.097E-04,&
     3.531E-06, 3.522E-06, 3.500E-06, 3.478E-06, 3.444E-06, 3.391E-06,&
     3.953E-06, 3.838E-06, 3.034E-06, 3.616E-06, 2.494E-06, 2.859E-06,&
     2.934E-06, 2.950E-06, 2.800E-06, 3.322E-06, 3.463E-06, 4.294E-06,&
     8.687E-07/
		call apply_ckd(ni, k, u, f, p, t, tau, coefk, .true.)
      end subroutine ck_18
	  
!c *********************************************************************
      subroutine window18 (amnt,twin,patm,temp,fac)
!c     Parameterization of CKD_2.4 continuum over 10430 to 10960 cm-1 band
!c     Fit inaccurate for amnt < 0.01 g/cm*2;
!c          however optical depths < 0.001.
!c     INPUT:
!c     amnt = h2O layer amount (g/cm**2)
!c     patm= pressure (atm)
!c     temp = temperature (k)
!c     OUTPUT:
!c     twin = parameterized CKD_2.1 optical depth over band
!c
	  real, intent(in) :: amnt, patm, temp, fac
	  real, intent(inout) :: twin

      integer, parameter :: ncoef=7

      real ::  aa(ncoef), ph2o, tau_log
!c19   data aa/-2.119,1.002,-5.192e-03,0.9057,4.698,-1.894e-03,5.871e-02/
!c18   data aa/-1.837,1.002,-4.358e-03,0.9358,4.342,-1.518e-03,3.475e-02/
      data aa/-1.837,1.002,-4.358e-03,0.9358,4.342,-1.518e-03,3.475e-02/

      ph2o = amnt*(8.314e+07*temp)/(fac*1.0e+05*18.01534*1.01325e+06)
!c
          tau_log   = aa(1)              + &
                     aa(2) * log(amnt)  +&
                     aa(3) * temp       +&
                     aa(4) * patm       +&
                     aa(5) * ph2o       +&
                     aa(6) * amnt       +&
                     aa(7) * log(ph2o)

        twin = exp ( tau_log )
!c
        end subroutine window18
	
end module modis18
