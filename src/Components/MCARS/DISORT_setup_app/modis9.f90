module modis9

	use CKD_common

	implicit none
	private
	
	integer, parameter :: no3 = 1
	
	public :: ch9, init_ch9


contains

	subroutine init_ch9(modis_channel)
	
		type(ckd_type), intent(inout) :: modis_channel
		integer :: all_weights

		all_weights = no3
		
		modis_channel%num_weights = all_weights
		allocate(modis_channel%weights(all_weights))
		allocate(modis_channel%taus(all_weights,max_layers))
	
	end subroutine init_ch9

	subroutine ch9(u0, p0, t0, ux, fac, modis_channel)
!
!      This routine was developed for MODIS Channel 1 
!                        14800--16200 cm^{-1}
!
!-----------------------------------------------------------------
!     Channel 1   0.62-0.67 microns
!-----------------------------------------------------------------
!  INPUTS:
!     u0 --> water vapor amount [g/cm2/km]
!     p0 --> pressure in atmospheres
!     t0 --> temp in K
!     ux --> ozone [g/cm2/km]
!     fac --> layer thickness in km
!-----------------------------------------------------------------
		real, dimension(:), intent(in) :: u0, p0, t0, ux, fac
		type(ckd_type), intent(inout) :: modis_channel
		
		integer, parameter :: nlev = max_levels
		integer, parameter :: nlay = max_layers
		integer, parameter :: mlv = nlay

		real :: p(nlay), t(nlay),u(nlay), &
				uo3(nlay), zfac(nlay)
				
		real :: tauo3(no3,nlay), fo3(no3)
		
		integer :: i,m


      do  m=1,use_layers
        t(m)=(t0(m)+t0(m+1))/2.0
        p(m)=(p0(m)+p0(m+1))/2.0
        if(m > 26)then
         zfac(m)=5.0
        else
         zfac(m)=1.0
        end if
        uo3(m)=zfac(m)*(ux(m)+ux(m+1))/2.0
      enddo
	  
      call cko3(uo3,fo3,p,t,tauo3)
!     **fo3=1.0000**

      do  m=1,use_layers
        do i=1,no3 ! 1 k for h2o; 1 k for o3
           modis_channel%taus(i,m)=tauo3(i,m)
		   if (m == 1) modis_channel%weights(i) = fo3(1)
       end do
      enddo
	  

	end subroutine ch9

!c *********************************************************************
      subroutine cko3(uo3,f,p,t,tau)
	  
		real, dimension(:), intent(in) :: p, t, uo3
		
		real, dimension(:), intent(inout) :: f
		real, dimension(:,:), intent(inout) :: tau
		
		
		real :: k(no3)
		integer :: m, i

      f(1)=1.000000
      k(1)=13.926e-23*6.023e+23/47.9982
      do i=1,no3
      do  m=1,use_layers
          tau(i,m)=k(i)*uo3(m)
        end do
      end do

      end subroutine cko3



end module modis9

