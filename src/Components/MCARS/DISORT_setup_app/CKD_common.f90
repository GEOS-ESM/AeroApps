module CKD_common

implicit none

integer, parameter :: num_pressures = 19
real, parameter :: stp(num_pressures) = (/  0.251, 0.398, 0.631, 1.000, 1.58, 2.51,  &
      	             3.98, 6.31, 10.0, 15.8, 25.1, 39.8, 63.1, &
      	             100.0, 158.0, 251.0, 398.0, 631.0, 1000.0 /)

integer, parameter :: max_layers = 26 ! was 34
integer, parameter :: max_levels = 27 ! was 35
integer, parameter :: max_weights = 6

integer :: use_layers
integer :: use_levels
!$omp threadprivate(use_levels, use_layers)

type ckd_type

	integer :: num_weights
	real, dimension(:), allocatable :: weights
	real, dimension(:,:), allocatable :: taus
  
end type ckd_type

contains

subroutine apply_ckd(ni, k, u, f, p, t, tau, coefk, doppler)
	integer, intent(in) :: ni
	real, dimension(:,:,:), intent(in) :: coefk
	real, dimension(:), intent(in) :: p, t, u, k

	real, dimension(:,:), intent(inout) :: tau
	real, dimension(:), intent(inout) :: f
	logical, intent(in) :: doppler
	
	real :: pmb(max_layers), x1, x2, fkg(max_weights,max_layers)
	
	integer :: i, m, ml
	integer :: mlv
	
	mlv = use_layers
	
	do i=1,ni
		do m=1,mlv
			ml=1
			pmb(m)=p(m)*1013.25
			if(pmb(m) < stp(1))then
		
				x1=coefk(i,1,1)+coefk(i,2,1)*(t(m)-250.0) + coefk(i,3,1)*(t(m)-250.0)**2
				if (doppler) then 
					fkg(i,m)=x1 ! Doppler broadening
				else 
					fkg(i,m)=x1*pmb(m)/stp(1) ! pressure broadening
				endif
			else if (pmb(m) > stp(num_pressures)) then
				x1=coefk(i,1,num_pressures-1)+coefk(i,2,num_pressures-1)*(t(m)-250.0) &
						+coefk(i,3,num_pressures-1)*(t(m)-250.0)**2
				x2=coefk(i,1,num_pressures)+coefk(i,2,num_pressures)*(t(m)-250.0) &
						+coefk(i,3,num_pressures)*(t(m)-250.0)**2
				fkg(i,m)=x1+(x2-x1)/(stp(num_pressures)-stp(num_pressures-1)) *(pmb(m)-stp(num_pressures-1))
			else
				do while(pmb(m) >= stp(ml))
					ml=ml+1
				end do
				x1=coefk(i,1,ml-1)+coefk(i,2,ml-1)*(t(m)-250.0) + coefk(i,3,ml-1)*(t(m)-250.0)**2
				x2=coefk(i,1,ml)+coefk(i,2,ml)*(t(m)-250.0)  + coefk(i,3,ml)*(t(m)-250.0)**2
				fkg(i,m)=x1+(x2-x1)/(stp(ml)-stp(ml-1)) *(pmb(m)-stp(ml-1))
			end if
		end do
	end do
	do i=1,ni
		do m=1,mlv
			tau(i,m)=k(i)*u(m)*fkg(i,m)
        end do
	end do
	

end subroutine apply_ckd
					 





end module CKD_common
