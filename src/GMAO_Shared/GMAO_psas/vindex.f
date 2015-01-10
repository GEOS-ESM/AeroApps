	subroutine vindex(nt,xt,n,xn,ln)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: vindex - vector version of table lookup
!
! !SYNOPSIS:
!	subroutine vindex(nx,xt,n,xn,ln)
!	integer nt			! intent(in)
!	real	xt(nt)			! intent(in)
!	integer n			! intent(in)
!	real	xn(n)			! intent(in)
!	integer	ln(n)			! intent(out)
!
! !INPUT PARAMETERS:
!	integer nt		size of the table xt[]
!	real	xt(nt)		a table
!	integer n		size of the search list
!	real	xn(n)		a search list
!
! !OUTPUT PARAMETERS:
!	integer	ln(n)		indices to the table
!
! !DESCRIPTION: vindex search a table for a list of values and return
!		the indices of the searched values to table entries.
!		It is hoped this vectorized search process may be faster
!		than binary search
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
!	18May96 - J. Guo - Reprogrammed the searching strategy
! 	21Nov95 - J. Guo
!		- create the subroutine
!		- added the prolog
!_______________________________________________________________________

	implicit none
	integer nt	! intent(in)
	real	xt(nt)	! intent(in)
	integer n	! intent(in)
	real	xn(n)	! intent(in)
	integer	ln(n)	! intent(out)

!=======================================================================
!   ..Local variables

	integer i,it
!=======================================================================
!   ..A statement function and its arguments

	logical inorder
	real x,x1,x2
	inorder(x1,x,x2)=x.gt.min(x1,x2).and.x.lt.max(x1,x2).or.x.eq.x1
!=======================================================================


	do i=1,n
	  ln(i)=1	! you will be in trouble if nt==0
	end do
	if(nt.eq.1) return

	do i=1,n
	  do it=1,nt-1
	    if(inorder(xt(it),xn(i),xt(it+1))) then
	      ln(i)=it
	      exit
	    endif
	  end do
	  if(nt.gt.1) then
	    if(xn(i).eq.xt(nt).or.inorder(xt(1),xt(nt),xn(i)))
     &		ln(i)=nt-1
	  endif
	end do
	end
!.
