	integer function listvals(nv,vals,line)
!-----------------------------------------------------------------------
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE:  listvals
! 
! !DESCRIPTION: list data values in the most compact formats
!
! !CALLING SEQUENCE:
!	integer  listvals
!	external listvals
!	l=listvals(nv,vals,line)
!
! !INPUT PARAMETERS:
!	integer		nv		! number of vals[]
!	real		vals(nv)	! values
!
! !OUTPUT PARAMETERS:
!	character*(*)	line		! output buffer
!
! !BUGS:
!	Assumed upto 5 digits before and 2 digits after the decimal
!	point, which is reasonable for most real*4 values and some
!	physical quantities, such as pressures.
!
! !SEE ALSO:
!
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
! !REVISION HISTORY:
! 	23Aug95 - J. Guo	- added the prolog
!
!-----------------------------------------------------------------------
	implicit none
c   ..Inputs
	integer		nv		! number of vals[]
	real		vals(nv)	! values
c   ..Outputs
	character*(*)	line		! output buffer

c	..Locals
	character*8 bufr
	logical next
	integer i,j,k
	integer ln,ls,l

c   ..Checking buffer size
	ln=len(line)

c   ..For each item in vals(:), upto where the output buffer permits
	ls=0
	k=0
	do while(ls.lt.ln.and.k.lt.nv)
	  k=k+1

c	..Assumed upto [5 digits].[2 digits].

	  write(bufr,'(f8.2)') vals(k)

c	..Checking insignificant leading bytes
	  i=0
	  next=.true.
	  do while(next.and.i.lt.8)
	    i=i+1
	    next=bufr(i:i).eq.' '.or.bufr(i:i).eq.'0'
	  end do

c	..Checking insignificant trailing bytes
	  j=8+1
	  next=.true.
	  do while(next.and.j.gt.i)
	    j=j-1
	    next=bufr(j:j).eq.' ' .or. bufr(j:j).eq.'0'
	  end do

c	..Skipping the dicimal point "."?
	  if(bufr(j:j).eq.'.') then
	    if(j.gt.i) then
	      j=j-1
	    else
	      bufr(j:j)='0'
	    endif
	  endif

	  l=ls+1
	  ls=min(l + j-i+1,ln)
	  line(l:ls)=' '//bufr(i:j)

	end do

c   ..Clean up for the output
	if(ls.lt.ln) line(ls+1:ln)=' '
	listvals=ls
	end
