subroutine rdvctbl(rc_vctbl,name,desc,mxlev,nlev,plev,vctabl,istat)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: rdvctbl - read a resource of a vCor table
!
! !INTERFACE: (to do)
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
! 	11Jan96 - J. Guo	- programmed and added the prolog
!_______________________________________________________________________

use m_inpak90, only : lablin,rdnext,getwrd,getstr,str2rn
use m_stdio,   only : stderr
implicit none

!  Function of rdvctbl()
! =======================
!	Read a correlation function table in the form,

!	  > <vCorTable::>	# defined by rc_vctbl
!	  >
!	  > 	1000   1
!	  >	 850  .9  1
!	  >	 700  .6 .9  1
!	  >	 500  .5 .6 .9  1
!	  >	 400  ...
!	  >	 ...
!	  >
!	  > ::

!	Notice that following conditions are assumed,
!		1) Symmetric, C(i,j)=C(j,i).  So only the lower-half
!		   of the matrix is required and used.
!		2) C(i,i)=1, such that the diagonals are either not
!		   required, or used as a check if present.
!		3) Levels are sorted, but the order (ascending or
!		   descending) is not assumed.  It also means, if the
!		   levels are ascending (or descending) down the table,
!		   they must also be ascending (or descending) to the
!		   right (see 1).
!		4) No extrapolation is assumed.  It seems not a simple
!		   problem of positive definate is required.  Will leave
!		   it to the perent.

!  Arguments
! ===========
  character*(*),intent(in)	:: rc_vctbl	! label of a vCorTable
  character*(*),intent(out)	:: name		! name of the table
  character*(*),intent(out)	:: desc		! a short description
  integer,intent(in)		:: mxlev	! dimension of vctabl
  integer,intent(out)		:: nlev		! true size of vctabl
  real,   intent(out)		:: plev(mxlev)	! levels of vctabl
  real,   intent(out)		:: vctabl(mxlev,mxlev)	! vCorTable
  integer,intent(out)		:: istat	! input status

!  Return code of istat
! ======================
!	..For istat =	0, normal exit
!		       -1, table not found
!			1, invalid level
!			2, duplicate levels
!			3, levels not sorted
!			4, invalid correlation value
!			5, incomplete correlation table line, l<n-1
!			6, diagonal element is not 1
!			7, more in the file, but the array is full
!			8, empty table
!
!	  If istat >= 0, nlev has been assigned value to the current
!	level.
!
!  History
! =========
!	21Mar95 - Jing G. 	- Initial code
!=======================================================================

!  Local parameters and variables
! ================================
	character*7 myname
	parameter(myname='rdvctbl')

	integer lnblnk
	external lnblnk

	include "realvals.h"

	character*64 str
	integer ios	! data buffer processing status
	real val	! tmp real value
	integer n	! level counter
	integer l

	real a,b
	logical AboutEqual
	AboutEqual(a,b)=abs(a-b).le.				&
	  max(rsmlval,max(abs(a),abs(b))*rfrcval)

!  Table processing
! ==================
		! Locate and verify the existence of the label.

	call lablin(trim(rc_vctbl))
	call rdnext(ios)
	if(ios.eq.2) then
	  write(stderr,'(3a)') myname,				&
	    ': table not found, ',rc_vctbl
	  istat=-1
	  return	! bad news to the parent
	endif


		! Processe the table line by line

	istat=0		! if normal return
	name=' '
	desc=' '

	n=0
	nlev=n
	do while(ios.eq.0.and.n.lt.mxlev)
	  n=n+1
	  nlev=n	! current position if return to parent

		! Process the line token by token
	  call getwrd(ios,str)		! leading token mark the level
	  if(ios.ne.0) then
	    l=max(lnblnk(str),1)
	    write(stderr,'(4a)') myname,			&
	      ': unexpected error from getwrd(), "',str(1:l),'"'
	    istat=1
	    return	! bad news to the parent
	  endif

	  val=str2rn(str,ios)
	  if(ios.ne.0) then		! possibly a name
	    if(name.ne.' ') then	! alread defined
	      l=max(lnblnk(str),1)
	      write(stderr,'(4a)') myname,			&
		': unexpected string found, "',str(1:l),'"'
	      istat=1
	      return
	    endif

	    name=str			! it is a name for the function
	    call getstr(ios,desc)	! get the descrition string
	    if(ios.ne.0) desc='?'	! ignore the erro

	  else

	    plev(n)=val

		! Verify the levels
	    if(n.gt.1) then	! if the levels are too close
	      if( AboutEqual(plev(n),plev(n-1)) ) then
	        l=max(lnblnk(str),1)
	        write(stderr,'(4a)') myname,			&
		  ': value too close to the previous level, "',	&
		str(1:l),'"'
	        istat=2
	        return
	      endif
	    endif

		! More of verify the levels
	    if(n.gt.2) then	! if the levels are not sorted.
	      if( (plev(n)  .gt.plev(n-1)) .neqv.		&
		  (plev(n-1).gt.plev(n-2))			) then
		write(stderr,'(2a,3g10.2)') myname,		&
		  ': not sorted ,',plev(n-2),plev(n-1),plev(n)
	      
		istat=3
		return
	      endif
	    endif

		! Loop over tokens for correlation function values
	    call getwrd(ios,str)
	    l=0
	    do while(ios.eq.0.and.l.lt.n)
	      l=l+1
	      val=str2rn(str,ios)
	      if(ios.eq.0) then
	        if(val.gt.1.) then
		  l=max(lnblnk(str),1)
		  write(stderr,'(4a)') myname,			&
		    ': invalid value, "',str(1:l),'"'
		  istat=4
		  return	! bad news to the parent
	        else
	          vctabl(l,n)=val
	          vctabl(n,l)=val	! symmetry considered
	        endif
	      else
	        l=max(lnblnk(str),1)
	        write(stderr,'(4a)') myname,			&
		  ': not a number, "',str(1:l),'"'
	        istat=4
	        return	! bad news to the parent
	      endif

	      call getwrd(ios,str)
	    end do

		! If the information is complete?
	    if(l.lt.n-1) then
	      write(stderr,'(2a,i2,a,i2,a)') myname,		&
	        ': incomplete correlation function, (l=',l,	&
		') < (n=',n,')'
	      istat=5
	      return

		! If the diagonal element is missing, assign one to it
	    elseif(l.eq.n-1) then
	      vctabl(n,n)=1.

		! If the diagonal element present (l=n), verify its
		! value is expected to be "about equal" to one.
	    else
	      if(AboutEqual(vctabl(n,n),1.)) then
	        vctabl(n,n)=1.

	      else
	        write(stderr,'(2a,g10.3e1)') myname,		&
		  ': invalid diagonal correlation value, ',	&
		  vctabl(n,n)
	        istat=6
	        return	! tell the parent the bad news
	      endif
	    endif
	  endif

	  call rdnext(ios)
	end do

		! If the table is not finished, tell the parent.
	if(ios.eq.0) then
	  write(stderr,'(2a,i2,a)') myname,			&
	    ': buffer is full (mxlev=',n,'), but more in the table'
	  istat=7
	  return
	endif

		! If the table is empty, tell the parent.
	if(n.eq.0) then
	  write(stderr,'(3a)') myname,				&
	    ': empty table, ',rc_vctbl
	  istat=8
	  return
	endif

	end
