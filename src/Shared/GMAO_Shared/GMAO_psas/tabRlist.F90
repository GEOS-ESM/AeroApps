!@interface
subroutine tabRlist(res,nVal,rVal,MXrTab,nTab,rTab)
!@end/interface
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: tabRlist - merge separable real values into a table
!
! !INTERFACE:
!	<@interface
!
! !DESCRIPTION:
!	Argument res defines the resolution of separable real values.
!
! !EXAMPLES:
!
! !BUGS:
!	values in rVal and rTab are expected to be "reasonable", which
!	means when the values are scaled by the value of res, they must
!	be in the range of the integers.
!
! !SEE ALSO:
!
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
! !REVISION HISTORY:
! 	01Mar96 - J. Guo	- programmed and added the prolog
!_______________________________________________________________________
!#define	_TRACE	!

use m_die,  only : die
implicit none

!@interface
  real,    intent(in)	  :: res		! expected resolution
  integer, intent(in)	  :: nVal		! size of rVal
  real,    intent(in)	  :: rVal(nVal)		! a list of real values
  integer, intent(in)	  :: MXrTab		! possible size of rTab
  integer, intent(inout)  :: nTab		! true size of rTab
  real,    intent(inout)  :: rTab(MXrTab)	! a sorted table
!@end/interface

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Locals/workspace
integer i,ias
real ares
integer,allocatable :: iVal(:),iTab(:)

character(len=*), parameter	:: myname='tabRlist'
!-----------------------------------------------------------------------
  if(nVal.eq.0) return		! default result is null table

  allocate(iVal(nVal),iTab(MXrTab),stat=ias)
  if(ias.ne.0) call die(myname,'allocate()',ias)

#ifdef	_TRACE
	_TRACE	write(*,'(2a,1p,e12.4,0p,3i5)') myname,	&
	_TRACE	  ': entered, ',res,nVal,MXrtab,nTab
#endif

  ares=abs(res)

#ifdef	_TRACE
	_TRACE	write(*,'(2a,1p,e12.4)') myname,	&
	_TRACE	  ': ares = ',ares
#endif

  do i=1,nVal
    iVal(i)=nint(rVal(i)/ares)
  end do

#ifdef	_TRACE
	_TRACE	write(*,'(2a)') myname,': loop of iVal done'
#endif

  do i=1,nTab
    iTab(i)=nint(rTab(i)/ares)
  end do

#ifdef	_TRACE
	_TRACE	write(*,'(2a)') myname,': loop of iTab done'
	_TRACE	write(*,'(2a)') myname,': calling tabIlist()'
#endif

  call tabIlist(nVal,iVal, MXrTab,nTab,iTab)

#ifdef	_TRACE
	_TRACE	write(*,'(2a)') myname,': from tabIlist()'
#endif

  do i=1,nTab
    rTab(i)=iTab(i)*ares
  end do

#ifdef	_TRACE
	_TRACE	write(*,'(2a)') myname,': loop of rTab done'
#endif

  deallocate(iVal,iTab)
end subroutine tabRlist
!.
