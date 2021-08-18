module modis3

	use CKD_common

	implicit none
	private
	
	integer, parameter :: ni = 1
	integer, parameter :: no3 = 1
	
	public :: ch3, init_ch3

contains

	subroutine init_ch3(modis_channel)
	
		type(ckd_type), intent(inout) :: modis_channel
		integer :: all_weights
		
		all_weights = ni*no3
		
		modis_channel%num_weights = all_weights
		allocate(modis_channel%weights(all_weights))
		allocate(modis_channel%taus(all_weights, max_layers))
	
	end subroutine init_ch3


!c     MODIS Channel 3 20875--21785 cm^{-1}
	  
	  subroutine ch3(u0, p0, t0, ux, fac, modis_channel)

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
	  
      call ck3(u,f,p,t,tau)
      call cko3(uo3,fo3,p,t,tauo3)
!c     **fo3=1.0000**


      do  m=1,use_layers
        do i=1,ni

           modis_channel%taus(i,m)=tau(i,m)+tauo3(m)
		   if (m == 1) modis_channel%weights(i) = f(i)*fo3(1)
		
        end do
	end do


      end subroutine ch3

!c *********************************************************************
      subroutine ck3(u,f,p,t,tau)

		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		
		
		real :: k(ni)
		real :: coefk(ni,3,num_pressures)
		integer :: i, jp, jt

      f(1)=1.000000
      k(1)=0.0000835718835471
	  
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/ &
     1.00034,1.00034,1.00034,1.00034,1.00035,1.00035, &
     1.00036,1.00037,1.00038,1.00040,1.00044,1.00051, &
     1.00061,1.00076,1.00095,1.00111,1.00051,0.99955, &
      0.99852, &
     -4.221E-03,-4.221E-03,-4.221E-03,-4.221E-03,-4.221E-03,-4.221E-03, &
     -4.221E-03,-4.222E-03,-4.222E-03,-4.222E-03,-4.222E-03,-4.223E-03, &
     -4.224E-03,-4.225E-03,-4.227E-03,-4.224E-03,-4.216E-03,-4.213E-03, &
     -4.207E-03, &
      1.682E-05, 1.682E-05, 1.682E-05, 1.682E-05, 1.682E-05, 1.682E-05, &
      1.682E-05, 1.682E-05, 1.682E-05, 1.682E-05, 1.682E-05, 1.683E-05, &
      1.683E-05, 1.684E-05, 1.684E-05, 1.676E-05, 1.676E-05, 1.677E-05, &
      1.674E-05/


		call apply_ckd(ni, k, u, f, p, t, tau, coefk, .true.)

      end subroutine ck3
!c *********************************************************************
      subroutine cko3(u,f,p,t,tau)
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:), intent(inout) :: tau
		
		
		
		real :: k(no3)
		integer :: m 

      f(1)=1.000000
      k(1)=4.30347e-22*6.023e+23/47.9982
        do m=1,use_layers
          tau(m)=k(1)*u(m)
        end do

      end subroutine cko3
	  
	  
end module modis3	  
