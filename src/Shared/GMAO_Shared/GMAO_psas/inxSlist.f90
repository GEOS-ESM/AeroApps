!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: inxSlist - index a list of strings respect to a table
!
! !INTERFACE:
    subroutine inxSlist(nTab,sTab,nVal,sVal,indx)
      implicit none
      integer, intent(in)	:: nTab		! size of sTab
      character*(*), intent(in)	:: Stab(nTab)	! a given string table
      integer, intent(in)	:: nVal		! size of sVal
      character*(*), intent(in)	:: sVal(nVal)	! a given string list
      integer, intent(out)	:: indx(nVal)	! pointers to the table
!   end subroutine inxSlist
!
! !DESCRIPTION:
!
! !EXAMPLES:
!
! !BUGS:
!   InxSlist() performs a simple searching with minimum consideration
!   of efficiency as a scalar algorithm.  However, it might do well on
!   a vector/multitask system, such as a Cray C90/J90.
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

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Locals
  integer iv,it
!-----------------------------------------------------------------------
  indx(1:nVal)=0
  do it=1,nTab
    do iv=1,nVal
      if(sVal(iv) == sTab(it)) indx(iv)=it
    end do
  end do
!_______________________________________________________________________
end subroutine inxSlist
!.
