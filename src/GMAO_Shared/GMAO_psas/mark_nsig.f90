!@interface
subroutine mark_nsig(n,sigU,sigC,kl)
!@end/interface
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: mark_nsig - mark-off observations with negative sig values
!
! !INTERFACE:
!	<@interface
!
! !DESCRIPTION: (to do)
! !EXAMPLES: (to do)
! !BUGS: (to do)
! !SEE ALSO: (to do)
! !SYSTEM ROUTINES: (to do)
!
! !REVISION HISTORY:
! 	24Apr96 - J. Guo	- (to do)
!_______________________________________________________________________

implicit none
!@interface

  integer, intent(in)	:: n		! size of sigmas
  real,    intent(in)	:: sigU(n),sigC(n)	! sigmas
  logical, intent(inout):: kl(n)

!@end/interface
!-----------------------------------------------------------------------
	integer i

  do i=1,n
	! if one of sigX is less than zero, this is not a valid data
	! point

    if(kl(i)) kl(i) = .not. (sigU(i) < 0. .or. sigC(i) < 0.)
  end do

end subroutine mark_nsig
!.
