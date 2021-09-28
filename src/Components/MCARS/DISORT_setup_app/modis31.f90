module modis31

	use CKD_common

	implicit none
	private
	
	integer, parameter :: ni = 5
	integer, parameter :: nco2 = 1
	
	public :: ch31, init_ch31
	
contains
	subroutine init_ch31(modis_channel)
	
		type(ckd_type), intent(inout) :: modis_channel
		integer :: all_weights
		
		all_weights = ni*nco2 
		
		modis_channel%num_weights = all_weights
		allocate(modis_channel%weights(all_weights))
		allocate(modis_channel%taus(all_weights, max_layers))
	
	end subroutine init_ch31

      subroutine ch31(u0, p0, t0, ux, fac, u_co2, modis_channel)
!
!      This routine was actually developed for MAS Channel 45 
!                        870--950 cm^{-1}
!
!-----------------------------------------------------------------
!     Channel 31   10.78 - 11.28 microns
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
				 uco2(nlay)
				
		real :: tau(ni,nlay),tauco2(nco2,nlay),  &
				f(ni), fco2(nco2),  xlog, twin
		
		integer :: i,m,ia

      do  m=1,use_layers
        t(m)=(t0(m)+t0(m+1))/2.0
        xlog=(log10(p0(m))+log10(p0(m+1)))/2.0
        p(m)=10.0**xlog
		xlog=(log10(u0(m))+log10(u0(m+1)))/2.0
        u(m)=fac(m)*10.0**xlog
        uco2(m)=p(m)*1.01325e+06/(1.3805e-16*t(m))*u_co2*fac(m)* &
     &          1.0e+05/6.023e+23*44.00995
      enddo
	  
      call ck_31(u,f,p,t,tau)  
      call ckco2_31(uco2,fco2,p,t,tauco2)

	  
      do  m=1,use_layers
        call window_31(u(m),twin,p(m),t(m),fac(m))
        do i=1,ni ! 5 k's for water
          modis_channel%taus(i,m)=tau(i,m)+twin+tauco2(1,m)
 			if (m == 1) modis_channel%weights(i) = f(i)*fco2(1)
         end do
      enddo

      end subroutine ch31

! *********************************************************************
      subroutine ck_31(u,f,p,t,tau)
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		
		real :: k(ni)
		real :: coefk(ni,3,num_pressures)
		integer :: i, jp, jt
      f(1)=0.822835
      f(2)=0.087911
      f(3)=0.065479
      f(4)=0.017907
      f(5)=0.005868
      k(1)=0.0011064239
      do i=2,ni
        k(i)=7.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/ & 
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, & 
     & 0.00000,0.00000,0.00029,0.00229,0.00982,0.03194, & 
     & 0.08593,0.17423,0.30684,0.50979,0.81227,1.26061, & 
     & 1.91371, & 
     & 0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00, & 
     & 0.000E+00,0.000E+00,1.060E-05,8.716E-05,4.058E-04,1.279E-03, & 
     & 2.937E-03,5.505E-03,9.504E-03,1.569E-02,2.471E-02,3.809E-02, & 
     & 5.726E-02, & 
     & 0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00, & 
     & 0.000E+00,0.000E+00,1.181E-07,1.027E-06,5.125E-06,1.608E-05, & 
     & 3.279E-05,5.805E-05,9.960E-05,1.638E-04,2.561E-04,3.921E-04, & 
     & 5.813E-04/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/ & 
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, & 
     & 0.00008,0.00107,0.00509,0.01654,0.04930,0.09949,& 
     & 0.15952,0.24791,0.38020,0.57054,0.83987,1.22580,& 
     & 1.72423,& 
     & 0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,& 
     & 3.000E-06,4.009E-05,2.075E-04,6.660E-04,1.684E-03,3.021E-03,& 
     & 4.748E-03,7.428E-03,1.137E-02,1.691E-02,2.460E-02,3.553E-02,& 
     & 5.061E-02,& 
     & 0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,& 
     & 3.187E-08,4.666E-07,2.612E-06,8.344E-06,1.861E-05,3.025E-05,& 
     & 4.803E-05,7.642E-05,1.169E-04,1.726E-04,2.472E-04,3.516E-04,& 
     & 5.155E-04/
      data ( (coefk(3,jt,jp), jp = 1, 19), jt = 1, 3)/& 
     & 0.00000,0.00000,0.00000,0.00000,0.00036,0.00167,& 
     & 0.00548,0.01590,0.03435,0.06100,0.09592,0.14384,& 
     & 0.21366,0.31323,0.45488,0.64144,0.87844,1.16403,& 
     & 1.51545,& 
     & 0.000E+00,0.000E+00,0.000E+00,0.000E+00,1.407E-05,7.285E-05,& 
     & 2.355E-04,5.983E-04,1.145E-03,1.909E-03,2.944E-03,4.403E-03,& 
     & 6.443E-03,9.303E-03,1.339E-02,1.904E-02,2.605E-02,3.426E-02,& 
     & 4.210E-02,& 
     & 0.000E+00,0.000E+00,0.000E+00,0.000E+00,1.669E-07,9.600E-07,& 
     & 3.106E-06,7.331E-06,1.270E-05,2.023E-05,3.104E-05,4.686E-05,& 
     & 6.752E-05,9.578E-05,1.363E-04,1.972E-04,2.736E-04,3.655E-04,& 
     & 4.319E-04/
      data ( (coefk(4,jt,jp), jp = 1, 19), jt = 1, 3)/ & 
     & 0.00352,0.00370,0.00570,0.01083,0.01819,0.02827,& 
     & 0.04224,0.05983,0.08285,0.11599,0.16468,0.23336,& 
     & 0.31979,0.43260,0.56679,0.73543,0.92196,1.09536,& 
     & 1.19608,& 
     & 3.024E-04,3.233E-04,3.613E-04,5.000E-04,7.402E-04,1.037E-03,& 
     & 1.424E-03,1.915E-03,2.598E-03,3.591E-03,5.009E-03,6.996E-03,& 
     & 9.664E-03,1.311E-02,1.699E-02,2.055E-02,2.367E-02,2.673E-02,& 
     & 2.822E-02,& 
     & 5.533E-06,5.955E-06,5.912E-06,6.813E-06,9.388E-06,1.242E-05,& 
     & 1.604E-05,2.085E-05,2.834E-05,3.906E-05,5.315E-05,7.239E-05,& 
     & 1.015E-04,1.405E-04,1.848E-04,2.122E-04,2.247E-04,2.290E-04,& 
     & 2.253E-04/
      data ( (coefk(5,jt,jp), jp = 1, 19), jt = 1, 3)/& 
     & 1.73163,1.73079,1.73067,1.73030,1.72960,1.72857,& 
     & 1.72888,1.72634,1.72311,1.71366,1.69383,1.65546,& 
     & 1.59832,1.51721,1.41274,1.27459,1.10668,0.91838,& 
     & 0.73202,& 
     & 4.422E-02,4.419E-02,4.415E-02,4.408E-02,4.399E-02,4.388E-02,& 
     & 4.373E-02,4.352E-02,4.313E-02,4.260E-02,4.179E-02,4.056E-02,& 
     & 3.873E-02,3.625E-02,3.319E-02,2.964E-02,2.565E-02,2.113E-02,& 
     & 1.687E-02,& 
     & 4.031E-04,4.028E-04,4.021E-04,4.013E-04,4.002E-04,3.992E-04,& 
     & 3.965E-04,3.947E-04,3.887E-04,3.814E-04,3.720E-04,3.595E-04,& 
     & 3.385E-04,3.113E-04,2.777E-04,2.438E-04,2.102E-04,1.748E-04,& 
     & 1.427E-04/
		call apply_ckd(ni, k, u, f, p, t, tau, coefk, .false.)

      end subroutine ck_31
! *********************************************************************
      subroutine ckco2_31(u,f,p,t,tau)
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		
		real :: k(nco2)
		integer :: i, m 
      f(1)=1.00000
      do  m=1,use_layers
        k(1)=0.912609*0.01047837*(1.0000+4.1069e-02*(t(m)-250.0)+ & 
     &        5.1759e-04*(t(m)-250.0)**2)
        tau(1,m)=k(1)*u(m)
      end do

      end subroutine ckco2_31
! **********************************************************************
      subroutine window_31 (amnt,twin,patm,temp,fac)
!     Parameterization of CKD_2.1 continuum over 835 to 980 cm-1 band
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
!
      data aa/9.374,0.9999,-0.02441,0.01989,-0.7838,4.191e-05,0.9918/
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

	end subroutine window_31
	
end module modis31