!@interface
    subroutine tabSlist(nVal,sVal,mxsTab,nTab,sTab)
!@end/interface
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: tabSlist - merge into a table of sorted unique strings
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
!@interface

      use m_die,  only : die
      implicit none

      integer, intent(in)	:: nVal		! size of sVal
      character*(*), intent(in)	:: sVal(nVal)	! a list of integers
      integer, intent(in)	:: mxsTab	! possible size of sTab
      integer, intent(inout)	:: nTab		! true size of sTab
      character*(*), intent(inout)	:: sTab(mxsTab)	! a sorted table

!@end/interface
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Locals/workspace
  integer			:: it,iv,in,istat
  integer, allocatable		:: indx(:)	! (nVal)
  character(len=len(sTab)),allocatable :: newTab(:)	! (mxsTab)

  integer, parameter	:: icap=ichar('A')-ichar('a')
  character(len=*), parameter	:: myname = 'tabSlist'
!-----------------------------------------------------------------------
  if(nVal.eq.0) return		! default result is original table

	allocate(indx(nVal), newTab(mxsTab), stat=istat)
	if(istat.ne.0) call die(myname,'allocate()',istat)

	! sort sVal array

  call indexxs(nVal,sVal,indx)	! indexed heap-sort sVal
  
  	! merge the sorted sVal and the sorted sTab (assumed) into a
  	! temporary array newTab

  in=0

  it=1
  iv=1
  do while(it.le.nTab .or. iv.le.nVal)
    if(it.le.nTab .and. iv.le.nVal) then	! make a comparison
      if(sTab(it).lt.sVal(indx(iv))) then
        call getTab_	! store the current sTab entry into newTab
        call nextTab_	! go to the next different sTab entry
      elseif(sTab(it).gt.sVal(indx(iv))) then
        call getVal_	! store the current sVal entry into newTab
        call nextVal_	! go to the next different sVal entry
      elseif(sTab(it).eq.sVal(indx(iv))) then
        call getTab_	! store the current sTab entry into newTab
        call nextTab_	! go to the next different sTab entry
        call nextVal_	! go to the next different sVal entry
      endif
    elseif(it.le.nTab) then		! iv .gt. nVal
      call getTab_	! store the current sTab entry into newTab
      call nextTab_	! go to the next different sTab entry
    elseif(iv.le.nVal) then		! it .gt. nTab
      call getVal_	! store the current sVal entry into newTab
      call nextVal_	! go to the next different sVal entry
    endif
  end do
  
  nTab=min(in,mxsTab)		! a new size
  sTab(1:nTab)=newTab(1:nTab)	! a new table returned.

  deallocate(indx,newTab)
contains
  subroutine getTab_
    in=in+1
    if(in.le.mxsTab) newTab(in)=sTab(it)
  end subroutine getTab_

  subroutine getVal_
    in=in+1
    if(in.le.mxsTab) newTab(in)=sVal(indx(iv))
  end subroutine getVal_

  subroutine nextTab_
    character(len=len(sTab))	:: curTab
    logical next
    curTab=sTab(it)		! with the current sTab(it) value
    it=it+1
    next=it.le.nTab
    if(next) next=curTab.eq.sTab(it)
    do while(next)
      it=it+1
      next=it.le.nTab
      if(next) next=curTab.eq.sTab(it)
    end do
  end subroutine nextTab_

  subroutine nextVal_
    character(len=len(sVal))	:: curVal
    logical next
    curVal=sVal(indx(iv))	! with the current sVal(indx(iv)) value
    it=it+1
    next=iv.le.nVal
    if(next) next=curVal.eq.sVal(indx(iv))
    do while(next)
      iv=iv+1
      next=iv.le.nVal
      if(next) next=curVal.eq.sVal(indx(iv))
    end do
  end subroutine nextVal_
!_______________________________________________________________________
end subroutine tabSlist
!.
