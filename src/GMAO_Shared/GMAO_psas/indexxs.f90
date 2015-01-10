!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: indexxs - indexed heap-sorting of a list of strings
!
! !INTERFACE:
    subroutine indexxs(n,str,indx)
      implicit none
      integer, intent(in)	:: n		! number of strings
      character(len=*), intent(in) :: str(n)	! a list of strings
      integer, intent(out)	:: indx(n)	! a sorted index
!   end subroutine ihpsrt_str

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
! 	10Jan96 - J. Guo	- programmed and added the prolog
!_______________________________________________________________________

!   ..Locals
  integer i,j,l,ir,indxt
  character(len=len(str))	:: tStr

!-----------------------------------------------------------------------

  do i=1,n			! initialize the indices
    indx(i)=i
  end do
  if(n.le.1) return

  l = n/2 + 1
  ir = n

  do
    if( l.gt.1 ) then
      l = l - 1
      indxt = indx(l)
      tStr = str(indxt)

    else
      indxt = indx(ir)
      tStr = str(indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if( ir.eq.1 ) then
        indx(1) = indxt
        return			! <== the exit
      endif
    endif

    i = l
    j = l + l

    do while( j.le.ir )
      if( j.lt.ir ) then
	if ( str(indx(j)).lt.str(indx(j+1)) ) j = j + 1
      endif

      if( tStr.lt.str(indx(j)) ) then
        indx(i) = indx(j)
        i = j
        j = j + j

      else
        j = ir + 1

      endif

    end do

    indx(i) = indxt

  end do
!
end
!.
