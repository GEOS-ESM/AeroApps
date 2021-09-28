module modis2

	use CKD_common

	implicit none
	private
	
	integer, parameter :: ni = 1
	integer, parameter :: no3 = 1
	
	public :: ch2, init_ch2

contains
	subroutine init_ch2(modis_channel)
	
		type(ckd_type), intent(inout) :: modis_channel
		integer :: all_weights
		
		all_weights = ni*no3
		
		modis_channel%num_weights = all_weights
		allocate(modis_channel%weights(all_weights))
		allocate(modis_channel%taus(all_weights, max_layers))
	
	end subroutine init_ch2


      subroutine ch2(u0, p0, t0, ux, fac, modis_channel)
!
!      This routine was developed for MODIS Channel 2 
!                        11415--11890 cm^{-1}
!
!-----------------------------------------------------------------
!     Channel 2   0.841-0.876 microns
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

		real :: p(nlay), t(nlay),u(nlay), uo3(nlay)
				
		real :: tau(ni,nlay),tauo3(nlay), &
				f(ni), fo3(no3)
		real :: xlog

		integer :: m, i 

      do m=1,use_layers
        t(m)=(t0(m)+t0(m+1))/2.0
        xlog=(log10(p0(m))+log10(p0(m+1)))/2.0
        p(m)=10.0**xlog
        u(m)=fac(m)*10.0**xlog
        uo3(m)=fac(m)*(ux(m)+ux(m+1))/2.0
      enddo
	  
      call ck_2(u,f,p,t,tau)
      call cko3_2(uo3,fo3,p,t,tauo3)
!     **fo3=1.0000**


      do m=1,use_layers
        do i=1,ni ! 1 k for h2o; 1 k for o3
           modis_channel%taus(i,m)=tau(i,m)+tauo3(m)
		   if (m == 1) modis_channel%weights(i) = f(i)*fo3(1)
       end do
      enddo
	  
	  end subroutine ch2

! *********************************************************************
       subroutine ck_2(u,f,p,t,tau)
	   
		real, dimension(:), intent(in) :: p, t, u
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
	   
		real :: k(ni)
		integer :: m 

		f(1)=1.00000
		do m=1,use_layers
			k(1)=0.85608*0.001703239*(1.0000+6.5376e-04*(t(m)-250.0)-8.7657e-06*(t(m)-250.0)**2)
			tau(1,m)=k(1)*u(m)
		end do
      end subroutine ck_2
! *********************************************************************
      subroutine cko3_2(uo3,f,p,t,tau)
	  
		real, dimension(:), intent(in) :: p, t, uo3
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:), intent(inout) :: tau
		
		
		
		real :: k(no3)
		integer :: m 
	  
		f(1)=1.000000
		k(1)=0.86313*5.51e-23*6.023e+23/47.9982
		do m=1,use_layers
			tau(m)=k(1)*uo3(m)
		end do
      end subroutine cko3_2


end module modis2
