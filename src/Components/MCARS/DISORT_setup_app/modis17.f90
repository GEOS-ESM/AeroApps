module modis17

	use CKD_common

	implicit none
	private
	
	integer, parameter :: ni = 4
	
	public :: ch17, init_ch17
	
contains

	subroutine init_ch17(modis_channel)
	
		type(ckd_type), intent(inout) :: modis_channel
		integer :: all_weights
		
		all_weights = ni
		
		modis_channel%num_weights = all_weights
		allocate(modis_channel%weights(all_weights))
		allocate(modis_channel%taus(all_weights, max_layers))
	
	end subroutine init_ch17

      subroutine ch17(u0, p0, t0, ux, fac, modis_channel)
!
!      This routine was developed for MODIS Channel 17 
!                        10870--11240 cm^{-1}
!
!-----------------------------------------------------------------
!     Channel 17   0.89-0.92 microns
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

      call ck_17(u,f,p,t,tau) 

      do  m=1,use_layers
		call window17(u(m),twin,p(m),t(m),fac(m))
        twinx(1)=0.1976226*twin*1.5155475
        twinx(2)=0.2828172*twin*1.5155475
        twinx(3)=0.3781571*twin*1.5155475
        twinx(4)=0.1394031*twin*1.5155475
        do i=1,ni ! 4 k's for h2o
          modis_channel%taus(i,m)=tau(i,m) + twinx(i)
 			if (m == 1) modis_channel%weights(i) = f(i)
        end do
	  enddo
 
      end subroutine ch17

! *********************************************************************
       subroutine ck_17(u,f,p,t,tau)
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		
		
		real :: k(ni)
		real :: coefk(ni,3,num_pressures)
		integer :: i, jp, jt
      f(1)=0.599013
      f(2)=0.273002
      f(3)=0.114623
      f(4)=0.013362
      k(1)=0.00705364457
      do i=2,ni
        k(i)=10.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     0.00065,0.00102,0.00158,0.00247,0.00384,0.00598,&
     0.00931,0.01449,0.02260,0.03524,0.05513,0.08641,&
     0.13566,0.21317,0.33381,0.51919,0.80269,1.22929,&
     1.83449,&
     5.000E-07, 7.500E-07, 1.250E-06, 1.875E-06, 2.875E-06, 4.125E-06,&
     6.000E-06, 8.875E-06, 1.262E-05, 1.775E-05, 2.362E-05, 3.350E-05,&
     4.775E-05, 7.263E-05, 1.495E-04, 2.481E-04, 4.360E-04, 7.143E-04,&
     1.066E-03,&
     0.000E+00,-6.250E-09,-6.250E-09,-1.562E-08,-2.187E-08,-3.437E-08,&
     -5.625E-08,-7.812E-08,-1.406E-07,-2.125E-07,-3.406E-07,-5.938E-07,&
     -9.500E-07,-1.503E-06,-2.644E-06,-3.909E-06,-5.819E-06,-7.356E-06,&
     -8.169E-06/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/&
     0.00494,0.00590,0.00729,0.00929,0.01213,0.01624,&
     0.02214,0.03069,0.04316,0.06129,0.08804,0.12744,&
     0.18546,0.27101,0.39620,0.57700,0.83311,1.17895,&
     1.59493,&
     4.600E-05, 4.662E-05, 4.800E-05, 5.113E-05, 5.550E-05, 6.175E-05,&
     7.125E-05, 8.437E-05, 1.021E-04, 1.240E-04, 1.511E-04, 1.845E-04,&
     2.193E-04, 2.505E-04, 2.974E-04, 2.911E-04, 2.911E-04, 2.954E-04,&
     1.974E-04,&
     1.500E-07, 1.219E-07, 9.375E-08, 5.312E-08, 1.250E-08,-5.000E-08,&
     -1.125E-07,-1.844E-07,-3.094E-07,-4.125E-07,-5.156E-07,-6.125E-07,&
     -7.125E-07,-8.813E-07,-1.278E-06,-1.203E-06,-2.772E-06,-4.384E-06,&
     -5.409E-06/
      data ( (coefk(3,jt,jp), jp = 1, 19), jt = 1, 3)/&
     0.26825,0.26872,0.26949,0.27109,0.27317,0.27687,&
     0.28259,0.29150,0.30603,0.32803,0.36190,0.41119,&
     0.48045,0.57199,0.68453,0.80965,0.93566,1.04223,&
     1.10583,&
     1.832E-03, 1.829E-03, 1.821E-03, 1.817E-03, 1.804E-03, 1.780E-03,&
     1.750E-03, 1.705E-03, 1.635E-03, 1.532E-03, 1.385E-03, 1.201E-03,&
     9.851E-04, 7.405E-04, 5.110E-04, 2.854E-04, 7.988E-05,-8.175E-05,&
     -2.304E-04,&
     -9.625E-07,-9.531E-07,-8.656E-07,-9.281E-07,-8.969E-07,-9.875E-07,&
     -9.062E-07,-5.781E-07,-6.188E-07,-3.812E-07,-2.406E-07,-3.813E-07,&
     -6.281E-07,-7.625E-07,-1.488E-06,-2.034E-06,-2.597E-06,-1.694E-06,&
     -8.094E-07/
      data ( (coefk(4,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     1.85523,1.85437,1.85302,1.85304,1.84965,1.84649,&
     1.84033,1.82951,1.81582,1.79135,1.75754,1.70620,&
     1.63233,1.53246,1.40520,1.25102,1.07782,0.90958,&
     0.74666,&
     -2.242E-03,-2.238E-03,-2.261E-03,-2.226E-03,-2.212E-03,-2.218E-03,&
     -2.186E-03,-2.138E-03,-2.063E-03,-1.983E-03,-1.887E-03,-1.719E-03,&
     -1.549E-03,-1.379E-03,-1.154E-03,-9.720E-04,-7.970E-04,-6.840E-04,&
     -5.844E-04,&
     2.531E-07, 2.531E-07, 9.594E-07, 2.344E-07, 2.219E-07,-4.438E-07,&
     -5.188E-07, 7.063E-07,-7.968E-07,-2.313E-07,-3.375E-07,-9.156E-07,&
     -5.406E-07, 5.000E-07, 6.906E-07, 1.500E-06, 1.925E-06, 2.206E-06,&
     1.497E-06/
		call apply_ckd(ni, k, u, f, p, t, tau, coefk, .true.)
      end subroutine ck_17
!c *********************************************************************
      subroutine window17 (amnt,twin,patm,temp,fac)
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
      data aa/-3.000,1.002,-6.760e-03,0.7417,3.309,-1.602e-03,1.739e-01/
!c
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

        end subroutine window17
		
end module modis17
