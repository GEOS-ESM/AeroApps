module modis32

	use CKD_common

	implicit none
	private
	
	integer, parameter :: ni = 5
	integer, parameter :: nco2 = 1
	
	public :: ch32, init_ch32
	
contains

	subroutine init_ch32(modis_channel)
	
		type(ckd_type), intent(inout) :: modis_channel
		integer :: all_weights
		
		all_weights = ni*nco2 
		
		modis_channel%num_weights = all_weights
		allocate(modis_channel%weights(all_weights))
		allocate(modis_channel%taus(all_weights, max_layers))
	
	end subroutine init_ch32

      subroutine ch32(u0, p0, t0, ux, fac, u_co2, modis_channel)
!
!      This routine was first developed for MAS Channel 46 
!                        810--870 cm^{-1}
!
!-----------------------------------------------------------------
!     Channel 32   11.49 - 12.34 microns
!-----------------------------------------------------------------
!  INPUTS:
!     u0 --> water vapor amount [g/cm2]
!     p0 --> pressure in atmospheres
!     t0 --> temp in K
!     ux --> ozone
!     dz --> layer thickness in km
!-----------------------------------------------------------------
!
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
        uco2(m)=p(m)*1.01325e+06/(1.3805e-16*t(m))*u_co2*fac(m)*  &
     &          1.0e+05/6.023e+23*44.00995
	  enddo
	  
      call ck_32(u,f,p,t,tau)  
      call ckco2_32(uco2,fco2,p,t,tauco2)


      do  m=1,use_layers
        call window_32(u(m),twin,p(m),t(m),fac(m))
        do i=1,ni ! 5 k's for water
          modis_channel%taus(i,m)=tau(i,m)+twin+tauco2(1,m)
 			if (m == 1) modis_channel%weights(i) = f(i)*fco2(1)
        end do
      enddo

      end subroutine ch32
	   
! *********************************************************************
      subroutine ck_32(u,f,p,t,tau)
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		

		real :: k(ni)
		real :: coefk(ni,3,num_pressures)
		integer :: i, jp, jt
      f(1)=0.812103
      f(2)=0.060532
      f(3)=0.082969
      f(4)=0.031873
      f(5)=0.012523
      k(1)=0.0016264818
      do i=2,ni
        k(i)=6.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00000,0.00011,0.00079,0.00349,0.01853, &
     & 0.06245,0.14125,0.27357,0.47424,0.79638,1.29226, &
     & 2.06120, &
     & 0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00, &
     & 0.000E+00,0.000E+00,5.863E-06,3.842E-05,1.654E-04,6.839E-04, &
     & 1.984E-03,4.222E-03,8.020E-03,1.386E-02,2.287E-02,3.650E-02, &
     & 5.679E-02, &
     & 0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00, &
     & 0.000E+00,0.000E+00,8.406E-08,5.481E-07,2.347E-06,7.953E-06, &
     & 2.137E-05,4.421E-05,8.254E-05,1.430E-04,2.338E-04,3.677E-04, &
     & 5.558E-04/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
     & 0.00000,0.00019,0.00120,0.00644,0.02912,0.07882, &
     & 0.13532,0.21992,0.34505,0.53277,0.82238,1.25501, &
     & 1.88556, &
     & 0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00, &
     & 0.000E+00,9.475E-06,5.435E-05,2.650E-04,1.041E-03,2.577E-03, &
     & 4.265E-03,6.990E-03,1.091E-02,1.649E-02,2.480E-02,3.706E-02, &
     & 5.445E-02, &
     & 0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00, &
     & 0.000E+00,1.356E-07,7.475E-07,3.350E-06,1.161E-05,2.808E-05, &
     & 4.588E-05,7.556E-05,1.179E-04,1.756E-04,2.554E-04,3.709E-04, &
     & 5.342E-04/
      data ( (coefk(3,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 0.00000,0.00000,0.00000,0.00000,0.00000,0.00016, &
     & 0.00091,0.00467,0.01512,0.03516,0.06444,0.10210, &
     & 0.15763,0.24078,0.36518,0.55754,0.83801,1.21862, &
     & 1.68015, &
     & 0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,6.363E-06, &
     & 3.835E-05,1.922E-04,5.534E-04,1.189E-03,2.104E-03,3.320E-03, &
     & 5.051E-03,7.595E-03,1.145E-02,1.726E-02,2.580E-02,3.791E-02, &
     & 5.114E-02, &
     & 0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,7.844E-08, &
     & 5.100E-07,2.413E-06,6.397E-06,1.338E-05,2.329E-05,3.698E-05, &
     & 5.541E-05,8.198E-05,1.231E-04,1.821E-04,2.691E-04,4.035E-04, &
     & 5.388E-04/
      data ( (coefk(4,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 0.00000,0.00000,0.00027,0.00139,0.00435,0.00983, &
     & 0.01887,0.03115,0.04795,0.07139,0.10561,0.15738, &
     & 0.23602,0.35127,0.50794,0.70236,0.91509,1.10383, &
     & 1.22070, &
     & 0.000E+00,0.000E+00,1.721E-05,7.939E-05,2.113E-04,4.257E-04, &
     & 7.283E-04,1.123E-03,1.655E-03,2.411E-03,3.494E-03,5.092E-03, &
     & 7.497E-03,1.109E-02,1.604E-02,2.186E-02,2.774E-02,3.304E-02, &
     & 3.612E-02, &
     & 0.000E+00,0.000E+00,2.834E-07,1.216E-06,2.926E-06,5.675E-06, &
     & 9.081E-06,1.347E-05,1.938E-05,2.773E-05,3.959E-05,5.655E-05, &
     & 8.115E-05,1.184E-04,1.723E-04,2.333E-04,2.869E-04,3.347E-04, &
     & 3.561E-04/
      data ( (coefk(5,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     & 1.60680,1.60598,1.60651,1.60835,1.61025,1.61230, &
     & 1.61657,1.62153,1.62441,1.62822,1.62479,1.60786, &
     & 1.56913,1.50141,1.40179,1.27122,1.10439,0.91802, &
     & 0.72779, &
     & 5.392E-02,5.385E-02,5.376E-02,5.364E-02,5.347E-02,5.319E-02, &
     & 5.284E-02,5.234E-02,5.165E-02,5.076E-02,4.976E-02,4.867E-02, &
     & 4.718E-02,4.499E-02,4.187E-02,3.793E-02,3.309E-02,2.761E-02, &
     & 2.202E-02, &
     & 5.956E-04,5.949E-04,5.935E-04,5.916E-04,5.890E-04,5.849E-04, &
     & 5.798E-04,5.721E-04,5.629E-04,5.487E-04,5.337E-04,5.195E-04, &
     & 5.018E-04,4.779E-04,4.446E-04,4.027E-04,3.557E-04,3.015E-04, &
     & 2.466E-04/
		call apply_ckd(ni, k, u, f, p, t, tau, coefk, .false.)

      end subroutine ck_32
! *********************************************************************
      subroutine ckco2_32(u,f,p,t,tau)
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		
		
		real :: k(nco2)
		integer :: i, m 
      f(1)=1.00000
      do  m=1,use_layers
        k(1)=0.942022*0.00534925618*(1.0000+4.0490e-02*(t(m)-250.0)+ &
     &       5.0470e-04*(t(m)-250.0)**2)
        tau(1,m)=k(1)*u(m)
      end do
      end subroutine ckco2_32
!**********************************************************************
      subroutine window_32(amnt,twin,patm,temp,fac)
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
      data aa/9.656,1.000,-0.02453,0.03056,-0.6620,-1.826e-04,0.9820/
!
      ph2o = amnt*(8.314e+07*temp)/(fac*1.0e+05*18.01534*1.01325e+06)
!
          tau_log   = aa(1)              + &
     &                aa(2) * log(amnt)  + &
     &                aa(3) * temp       + &
     &                aa(4) * patm       + &
     &                aa(5) * ph2o       + &
     &                aa(6) * amnt       + &
     &                aa(7) * log(ph2o)

        twin = exp ( tau_log )
!
	end subroutine window_32
	
end module modis32
