module modis4

	use CKD_common

	implicit none
	private
	
	integer, parameter :: ni = 1
	integer, parameter :: no3 = 1
	
	public :: ch4, init_ch4

contains

	subroutine init_ch4(modis_channel)
	
		type(ckd_type), intent(inout) :: modis_channel
		integer :: all_weights
		
		all_weights = ni*no3
		
		modis_channel%num_weights = all_weights
		allocate(modis_channel%weights(all_weights))
		allocate(modis_channel%taus(all_weights, max_layers))
	
	end subroutine init_ch4


!c     MODIS Channel 4 17700--18350 cm^{-1}


	  subroutine ch4(u0, p0, t0, ux, fac, modis_channel)
	  
		real, dimension(:), intent(in) :: u0, p0, t0, ux, fac
		type(ckd_type), intent(inout) :: modis_channel
		
		integer, parameter :: nlev = max_levels
		integer, parameter :: nlay = max_layers
		integer, parameter :: mlv = nlay

		real :: p(nlay), t(nlay),u(nlay), uo3(nlay)
		real :: zfac(nlay)		
				
		real :: tau(ni,nlay),tauo3(nlay), &
				f(ni), fo3(no3)

		integer :: m, i 
	  
      do  m=1,use_layers
        t(m)=(t0(m)+t0(m+1))/2.0
        p(m)=(p0(m)+p0(m+1))/2.0
		if (m >= 26) then 
			zfac(m) = 5.0
		else
			zfac(m) = 1.0
		endif
        u(m)=zfac(m)*(u0(m)+u0(m+1))/2.0
        uo3(m)=zfac(m)*(ux(m)+ux(m+1))/2.0
	  end do

      call ck4(u,f,p,t,tau)
      call cko34(uo3,fo3,p,t,tauo3)
!c     **fo3=1.0000**


      do  m=1,use_layers
        do i=1,ni

           modis_channel%taus(i,m)=tau(i,m)+tauo3(m)
		   if (m == 1) modis_channel%weights(i) = f(i)*fo3(1)
		
        end do
	end do
 	  

      end subroutine ch4

 
!c *********************************************************************
      subroutine ck4(u,f,p,t,tau)
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		
		real :: k(ni)
		real :: coefk(ni,3,num_pressures)
		integer :: i, jp, jt

      f(1)=1.000000
      k(1)=0.000302183032640
	  
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     1.01073,1.01073,1.01072,1.01071,1.01069,1.01066, &
     1.01062,1.01055,1.01044,1.01027,1.01000,1.00960,&
     1.00900,1.00810,1.00681,1.00473,1.00157,0.99843,&
      0.99541,&
     -6.975E-04,-6.975E-04,-6.975E-04,-6.974E-04,-6.972E-04,-6.970E-04,&
     -6.968E-04,-6.961E-04,-6.954E-04,-6.941E-04,-6.923E-04,-6.894E-04,&
     -6.849E-04,-6.784E-04,-6.689E-04,-6.486E-04,-6.356E-04,-6.261E-04,&
     -6.177E-04,&
      1.012E-06, 1.012E-06, 1.013E-06, 1.009E-06, 1.013E-06, 1.012E-06,&
      1.006E-06, 1.003E-06, 9.969E-07, 9.906E-07, 9.875E-07, 9.719E-07,&
      9.469E-07, 9.156E-07, 8.657E-07, 7.781E-07, 7.844E-07, 7.469E-07,&
      7.312E-07/

		call apply_ckd(ni, k, u, f, p, t, tau, coefk, .true.)

      end subroutine ck4
!c *********************************************************************
      subroutine cko34(u,f,p,t,tau)

		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:), intent(inout) :: tau
		
		
		
		real :: k(no3)
		integer :: m 

      f(1)=1.000000
      k(1)=3.51487e-21*6.023e+23/47.9982
      do  m=1,use_layers
          tau(m)=k(1)*u(m)
        end do



      end subroutine cko34
	  
	  
end module modis4
	  
	  
