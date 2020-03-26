	subroutine rdlevels(rc_lev,mxlev,nlev,plev,levtyp,istat)

	use m_inpak90, only : lablin,getwrd,str2rn
	use m_stdio,   only : stderr
	implicit none

c  Inputs
c ========
	character*(*) rc_lev	! lable of the entry (not a table)
	integer mxlev		! size of the array

c  Outputs
c =========
	integer nlev		! returned size of plev array
	real plev(mxlev)	! returned level values
	character*(*) levtyp	! type of the level
	integer istat		! input status

c  Locals
c ========
	include "realvals.h"

	character*8 myname
	parameter  (myname='rdlevels')

	integer l
	real val
	integer ios
	character*64 str

	integer lnblnk
	external lnblnk

	real a,b
	logical AboutEqual
	AboutEqual(a,b)=abs(a-b).le.
     &	  max(rsmlval,max(abs(a),abs(b))*rfrcval)

	istat=0

c  Search for the label
c ======================
c	..Reading levels for ObsErr*<class> tables
	call lablin(trim(rc_lev))
	call getwrd(ios,str)
	if(ios.ne.0) then
	  write(stderr,'(3a)') myname,
     &	    ': label not found, ',trim(rc_lev)
	  istat=-1
	  return
	endif

c	..Save the name
	levtyp=str

c	..Loop over tokens for each level value
	nlev=0
	call getwrd(ios,str)
	do while(ios.eq.0.and.nlev.le.mxlev)

	  val=str2rn(str,ios)
	  if(ios.ne.0) then		! not a number
	    l=max(lnblnk(str),1)
	    write(stderr,'(4a)') myname,
     &		': not a number, "',str(1:l),'"'
	    istat=1
	    return
	  endif

	  nlev=nlev+1
	  plev(nlev)=val

c	  ..Verify the levels
	  if(nlev.gt.1) then	! if the levels are too close
	    if( AboutEqual(plev(nlev),plev(nlev-1)) ) then
	      l=max(lnblnk(str),1)
	      write(stderr,'(4a)') myname,
     &		': value too close to the previous level, "',
     &		str(1:l),'"'
	      istat=2
	      return
	    endif
	  endif

c	  ..Verify the levels
	  if(nlev.gt.2) then	! if the levels are not sorted.
	    if(	(plev(nlev)  .gt.plev(nlev-1)) .neqv.
     &		(plev(nlev-1).gt.plev(nlev-2)) ) then
	      write(stderr,'(2a,3g10.2)') myname,
     &		': not sorted ,',plev(nlev-2),plev(nlev-1),plev(nlev)
	      istat=3
	      return
	    endif
	  endif

	  call getwrd(ios,str)
	end do

c	..If the table is not finished, tell the parent.
	if(ios.eq.0) then
	  write(stderr,'(2a,i3)') myname,
     &	    ': more in the line, but mxlev =',mxlev
	  istat=4
	  return
	endif

c	..If the table is empty, tell the parent.
	if(nlev.eq.0) then
	  write(stderr,'(3a)') myname,
     &	    ': empty table, ',trim(rc_lev)
	  istat=5
	  return
	endif

	end
