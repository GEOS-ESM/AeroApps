!@interface
    subroutine tabIlist(nVal,iVal,mxiTab,nTab,iTab)
!@end/interface
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: tabIlist - merge into a table of sorted unique integers
!
! !INTERFACE:
!	<@interface
!
! !DESCRIPTION:
!
! !EXAMPLES:
!
! !BUGS:
!
! !SEE ALSO:
!
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
! !REVISION HISTORY:
! 	03Jan96 - J. Guo	- programmed and added the prolog
!_______________________________________________________________________
!#define	_TRACE	!

!@interface
	use m_die,only : die
	implicit none

	integer, intent(in)	:: nVal		! size of iVal
	integer, intent(in)	:: iVal(nVal)	! a list of integers
	integer, intent(in)	:: mxiTab	! possible size of iTab
	integer, intent(inout)	:: nTab		! true size of iTab
	integer, intent(inout)	:: iTab(mxiTab)	! a sorted table

!@end/interface
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Locals/workspace
integer it,iv,in,istat
integer, allocatable :: indx(:)		! (nVal)
integer, allocatable :: newTab(:)	! (mxiTab)
character(len=*), parameter :: myname='tabIlist'
!-----------------------------------------------------------------------
#ifdef	_TRACE
		_TRACE	write(*,'(2a,3i5)') myname,	&
		_TRACE	  ': entered, ',nVal,mxiTab,nTab
#endif

  if(nVal.eq.0) return		! default result is null table

	allocate( indx(nVal), newTab(mxiTab), stat=istat)
	if(istat.ne.0) call die(myname,'allocate()',istat)

	! sort iVal

#ifdef	_TRACE
		_TRACE	write(*,'(2a)') myname,': calling indexxi()'
#endif

  call indexxi(nVal,iVal,indx)	! indexed heap-sort Ival
  
  	! merge the sorted iVal and the sorted iTab (assumed) into a
  	! temporary array newTab

#ifdef	_TRACE
		_TRACE	write(*,'(2a)') myname,': entering merge loop'
#endif
  in=0

  it=1
  iv=1
  do while(it.le.nTab .or. iv.le.nVal)
    if(it.le.nTab .and. iv.le.nVal) then	! make a comparison
      if(iTab(it).lt.iVal(indx(iv))) then
        call getTab_	! store the current iTab entry into newTab
        call nextTab_	! go to the next different iTab entry
      elseif(iTab(it).gt.iVal(indx(iv))) then
        call getVal_	! store the current iVal entry into newTab
        call nextVal_	! go to the next different iVal entry
      elseif(iTab(it).eq.iVal(indx(iv))) then
        call getTab_	! store the current iTab entry into newTab
        call nextTab_	! go to the next different iTab entry
        call nextVal_	! go to the next different iVal entry
      endif
    elseif(it.le.nTab) then		! iv .gt. nVal
      call getTab_	! store the current iTab entry into newTab
      call nextTab_	! go to the next different iTab entry
    elseif(iv.le.nVal) then		! it .gt. nTab
      call getVal_	! store the current iVal entry into newTab
      call nextVal_	! go to the next different iVal entry
    endif
  end do
  
#ifdef	_TRACE
		_TRACE	write(*,'(2a)') myname,': merge loop done'
#endif

  nTab=min(in,mxiTab)		! a new size
  iTab(1:nTab)=newTab(1:nTab)	! a new table returned.

  deallocate(indx,newTab)
!=======================================================================
contains
  subroutine getTab_
    in=in+1
    if(in.le.mxiTab) newTab(in)=iTab(it)
  end subroutine getTab_

  subroutine getVal_
    in=in+1
    if(in.le.mxiTab) newTab(in)=iVal(indx(iv))
  end subroutine getVal_

  subroutine nextTab_
    integer curTab
    logical next
    curTab=iTab(it)
    it=it+1
    next=it.le.nTab
    if(next) next=iTab(it).eq.curTab
    do while(next)
      it=it+1
      next=it.le.nTab
      if(next) next=iTab(it).eq.curTab
    end do
  end subroutine nextTab_

  subroutine nextVal_
    integer curVal
    logical next
    curVal=iVal(indx(iv))
    iv=iv+1
    next=iv.le.nVal
    if(next) next=iVal(indx(iv)).eq.curVal
    do while(next)
      iv=iv+1
      next=iv.le.nVal
      if(next) next=iVal(indx(iv)).eq.curVal
    end do
  end subroutine nextVal_
!_______________________________________________________________________
end subroutine tabIlist
!.
