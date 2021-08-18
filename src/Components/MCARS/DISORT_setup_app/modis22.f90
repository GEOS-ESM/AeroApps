module modis22

	use CKD_common

	implicit none
	private
	
	integer, parameter :: ni = 1
	integer, parameter :: nco2 = 1
	integer, parameter :: nch4 = 1
	integer, parameter :: nn2o = 3
	
	public :: ch22, init_ch22
	
contains

	subroutine init_ch22(modis_channel)
	
		type(ckd_type), intent(inout) :: modis_channel
		integer :: all_weights
		
		all_weights = ni*nco2*nch4*nn2o
		
		modis_channel%num_weights = all_weights
		allocate(modis_channel%weights(all_weights))
		allocate(modis_channel%taus(all_weights, max_layers))
	
	end subroutine init_ch22


      subroutine ch22(u0, p0, t0, ux, fac, u_co2, modis_channel)
!
!      This routine was developed for MODIS Channel 22 
!                        2505--2545 cm^{-1}
!
!-----------------------------------------------------------------
!     Channel 22   3.929 - 3.989 microns
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

		real :: p(nlay), t(nlay),u(nlay), xlog, twin, &
				 uco2(nlay), uch4(nlay), un2o(nlay)
				
		real :: tau(ni,nlay),tauco2(nco2,nlay), tauch4(nch4, nlay), &
				taun2o(nn2o, nlay), f(ni), fco2(nco2), fch4(nch4), fn2o(nn2o)
		
		integer :: i,m,ia, ib, id
	  
	  
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
        un2o(m)=p(m)*1.01325e+06/(1.3805e-16*t(m))*3.10e-07*fac(m)* &
     &          1.0e+05/6.023e+23*44.0128
      enddo
	  
      call ck_22(u,f,p,t,tau)
      call ckco2_22(uco2,fco2,p,t,tauco2)
      call ckch4_22(uch4,fch4,p,t,tauch4)
      call ckn2o_22(un2o,fn2o,p,t,taun2o)
	  
	  
      do  m=1,use_layers
        call n2cont_22(u(m),twin,p(m),t(m),fac(m))
		   do id=1, nn2o
            modis_channel%taus(id,m)=tau(1,m)+twin+tauco2(1,m) +tauch4(1,m)+taun2o(id,m)
  			if (m == 1) modis_channel%weights(id) = f(1)*fco2(1)*fch4(1)*fn2o(id)
          end do
       enddo
 	  

      end subroutine ch22

! *********************************************************************
!     h2o 2505-2545 cm^{-1}
      subroutine ck_22(u,f,p,t,tau)
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		
		
		real :: k(ni)
		real :: coefk(ni,3,num_pressures)
		integer :: i, jp, jt
      f(1)=1.00000
      k(1)=0.001190325514*0.916330
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     1.00239,1.00239,1.00240,1.00241,1.00242,1.00245, &
     1.00250,1.00258,1.00271,1.00291,1.00319,1.00341, &
     1.00308,1.00246,1.00216,1.00159,1.00061,0.99931, &
      0.99731, &
	 2.150E-02, 2.150E-02, 2.150E-02, 2.150E-02, 2.150E-02, 2.150E-02,&
      2.150E-02, 2.150E-02, 2.150E-02, 2.151E-02, 2.151E-02, 2.152E-02,&
      2.152E-02, 2.150E-02, 2.150E-02, 2.148E-02, 2.147E-02, 2.144E-02,&
      2.141E-02, &
      1.690E-04, 1.690E-04, 1.690E-04, 1.690E-04, 1.691E-04, 1.691E-04,&
      1.691E-04, 1.691E-04, 1.691E-04, 1.691E-04, 1.691E-04, 1.692E-04,&
      1.692E-04, 1.692E-04, 1.690E-04, 1.689E-04, 1.690E-04, 1.687E-04,&
      1.686E-04/

		call apply_ckd(ni, k, u, f, p, t, tau, coefk, .true.)
      end subroutine ck_22


! *********************************************************************
      subroutine n2cont_22 (amnt,twin,patm,temp,fac)
!     under construction
!     Parameterization of CKD_2.1 continuum over 2505 to 2545 cm-1 band
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

      integer, parameter :: ncoef=8

      real ::  aa(ncoef), ph2o, tau_log
!      data aa/-13.89,0.9790,-0.004063,-2.869, &
!     &         17.16,-5.184e-03,-0.9811,9.500/
      data aa/-14.57,1.0196,-0.004022,-3.795, &
              16.21,-6.307e-03,-0.9911,11.09/
!
      ph2o = amnt*(8.314e+07*temp)/(fac*1.0e+05*18.01534*1.01325e+06)
!
          tau_log   = aa(1)              + &
     &                aa(2) * log(amnt)  + &
     &                aa(3) * temp       + &
     &                aa(4) * patm       + &
     &                aa(5) * ph2o       + &
     &                aa(6) * amnt       + &
     &                aa(7) * log(ph2o)  + &
     &                aa(8) * patm**0.5

      twin = exp ( tau_log )
!
      end subroutine n2cont_22
!*********************************************************************
!     co2 2505--2545 cm^{-1}
      subroutine ckco2_22(u,f,p,t,tau)
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		

		real :: k(nco2)
		real :: coefk(nco2,3,num_pressures)
		integer :: i, jp, jt
      f(1)=1.00000
      k(1)=0.0026355786361*1.006450
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     1.00492,1.00407,1.00339,1.00332,1.00309,1.00280,&
     1.00237,1.00184,1.00125,1.00079,1.00030,0.99958,&
     0.99910,0.99901,0.99923,0.99952,0.99977,1.00027,&
      1.00096,&
     -1.615E-03,-1.608E-03,-1.617E-03,-1.614E-03,-1.610E-03,-1.605E-03,&
     -1.597E-03,-1.587E-03,-1.576E-03,-1.567E-03,-1.552E-03,-1.547E-03,&
     -1.551E-03,-1.551E-03,-1.555E-03,-1.556E-03,-1.557E-03,-1.559E-03,&
     -1.558E-03,&
      9.406E-07, 9.187E-07, 9.781E-07, 9.063E-07, 9.031E-07, 8.438E-07,&
      7.906E-07, 7.125E-07, 6.500E-07, 5.875E-07, 4.937E-07, 6.375E-07,&
      5.969E-07, 5.937E-07, 6.125E-07, 5.656E-07, 6.094E-07, 6.063E-07,&
      5.750E-07/

		call apply_ckd(nco2, k, u, f, p, t, tau, coefk, .true.)
      end subroutine ckco2_22
! *********************************************************************
!     ch4 2505--2545 cm^{-1}
      subroutine ckch4_22(u,f,p,t,tau)
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		

		real :: k(nch4)
		real :: coefk(nch4,3,num_pressures)
		integer :: i, jp, jt
      f(1)=1.00000
      k(1)=5.24800575521*0.939533
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     1.00094,1.00056,0.99995,0.99898,0.99862,0.99864,&
     0.99870,0.99880,0.99894,0.99915,0.99946,0.99987,&
     0.99954,0.99938,0.99925,0.99947,0.99977,1.00020,&
      1.00088,&
      1.010E-03, 1.012E-03, 1.015E-03, 1.017E-03, 1.002E-03, 1.003E-03,&
      1.003E-03, 1.002E-03, 1.001E-03, 9.993E-04, 9.977E-04, 9.974E-04,&
      1.007E-03, 1.006E-03, 9.998E-04, 1.001E-03, 9.963E-04, 9.945E-04,&
      9.877E-04,&
     -6.497E-06,-6.506E-06,-6.516E-06,-6.425E-06,-6.487E-06,-6.447E-06,&
     -6.441E-06,-6.434E-06,-6.416E-06,-6.387E-06,-6.362E-06,-6.353E-06,&
     -6.225E-06,-6.306E-06,-6.256E-06,-6.259E-06,-6.250E-06,-6.250E-06,&
     -6.212E-06/

		call apply_ckd(nch4, k, u, f, p, t, tau, coefk, .true.)
      end subroutine ckch4_22
! *********************************************************************
!c     n2o 2490--2545 cm^{-1}
      subroutine ckn2o_22(u,f,p,t,tau)
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		

		real :: k(nn2o)
		real :: coefk(nn2o,3,num_pressures)
		integer :: i, jp, jt
      f(1)=0.699544
      f(2)=0.263110
      f(3)=0.037346
      k(1)=5.47256912
      do i=2,nn2o
        k(i)=12.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     0.00000,0.00000,0.00005,0.00037,0.00126,0.00301,&
     0.00642,0.01376,0.02766,0.05202,0.08901,0.14114,&
     0.21298,0.31912,0.46262,0.64385,0.86901,1.16177,&
     1.50986,&
     0.000E+00, 0.000E+00, 7.500E-07, 5.500E-06, 1.650E-05, 3.763E-05,&
     7.775E-05, 1.531E-04, 2.804E-04, 4.650E-04, 7.316E-04, 1.156E-03,&
     1.700E-03, 2.598E-03, 3.805E-03, 5.384E-03, 7.631E-03, 1.069E-02,&
     1.447E-02,&
     0.000E+00, 0.000E+00, 0.000E+00, 1.875E-08, 5.000E-08, 1.969E-07,&
     5.187E-07, 8.406E-07, 1.416E-06, 2.063E-06, 2.478E-06, 3.269E-06,&
     4.803E-06, 6.481E-06, 9.363E-06, 1.468E-05, 1.984E-05, 2.427E-05,&
     3.048E-05/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/&
     0.01935,0.02257,0.02633,0.03175,0.03841,0.04788,&
     0.06070,0.07819,0.10079,0.12961,0.16680,0.21273,&
     0.27766,0.37219,0.50555,0.67912,0.89033,1.12570,&
     1.35741,&
     2.899E-04, 4.020E-04, 4.282E-04, 4.598E-04, 5.053E-04, 5.610E-04,&
     6.435E-04, 7.509E-04, 8.804E-04, 1.037E-03, 1.218E-03, 1.472E-03,&
     1.891E-03, 2.419E-03, 3.134E-03, 4.052E-03, 4.957E-03, 5.801E-03,&
     6.521E-03,&
     5.906E-07, 2.294E-06, 2.456E-06, 2.137E-06, 2.444E-06, 2.425E-06,&
     2.506E-06, 2.353E-06, 2.622E-06, 3.500E-06, 4.162E-06, 5.762E-06,&
     6.637E-06, 6.675E-06, 5.284E-06, 3.856E-06, 1.619E-06,-1.856E-06,&
     -4.981E-06/
      data ( (coefk(3,jt,jp), jp = 1, 19), jt = 1, 3)/&
     1.76065,1.76500,1.75741,1.74600,1.73040,1.70953,&
     1.67864,1.64214,1.59801,1.55205,1.50993,1.47032,&
     1.43215,1.38185,1.31124,1.21256,1.08220,0.92054,&
     0.74651,&
     7.911E-03, 8.626E-03, 8.632E-03, 8.600E-03, 8.566E-03, 8.485E-03,&
     8.405E-03, 8.273E-03, 8.089E-03, 7.891E-03, 7.634E-03, 7.429E-03,&
     7.161E-03, 6.907E-03, 6.537E-03, 6.005E-03, 5.286E-03, 4.420E-03,&
     3.466E-03,&
     -1.910E-05,-6.116E-06,-5.409E-06,-5.422E-06,-4.409E-06,-5.228E-06,&
     -4.241E-06,-5.003E-06,-5.541E-06,-5.237E-06,-6.162E-06,-5.819E-06,&
     -7.150E-06,-7.097E-06,-7.422E-06,-6.841E-06,-6.678E-06,-6.150E-06,&
     -4.897E-06/
		call apply_ckd(nn2o, k, u, f, p, t, tau, coefk, .false.)
      end subroutine ckn2o_22

end module modis22