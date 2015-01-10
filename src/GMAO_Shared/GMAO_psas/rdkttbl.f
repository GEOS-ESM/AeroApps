	subroutine rdkttbl(rc_kt,mxkt,
     &	  ktname,ktunit,ktdesc,ktmvar,istat)

	use m_inpak90, only : lablin,rdnext,getwrd,str2rn,getstr
	use m_stdio,   only : stderr

	implicit none

c  Inputs
c ========
	character*(*) rc_kt	! table label
	integer mxkt		! size of the storage arrays

c  Outputs
c =========
	character*(*) ktname(mxkt)	! variable names
	character*(*) ktunit(mxkt)	! variable units
	character*(*) ktdesc(mxkt)	! variable descriptions
	logical ktmvar(mxkt,mxkt)	! multi-variate flags
	integer istat		! input status code

c  Local Parameters and Variables
c ================================

	character*7 myname
	parameter(myname='rdkttbl')

	include "realvals.h"

c	..Locals
	character*64 str
	integer ios
	integer l,k,ls,ln
	integer ikt,jkt
	real val

	integer lnblnk
	external lnblnk

	real a
	logical NotAnInt
	NotAnInt(a)=abs(a-nint(a)).gt.rfrcval


c	..Locate the table and check the status
c	========================================
	call lablin(trim(rc_kt))
	call rdnext(ios)
	if(ios.eq.2) then
	  write(stderr,'(3a)') myname,
     &	    ': table not found, ',trim(rc_kt)
	  istat=-1
	  return	! bad news to the parent
	endif

	istat=0		! if normal return

	do ikt=1,mxkt
	  ktname(ikt)=' '
	  ktunit(ikt)=' '
	  ktdesc(ikt)=' '
	  do jkt=1,mxkt
	    ktmvar(ikt,jkt)=.false.
	  end do
	end do

c	..Processe the table line by line

	do while(ios.eq.0)

c	  ..Process the line token by token

c  #1	  ..The first token is the code number of the variable.  It is
c ====	  the number used in the input data "delta" (innovation vector)
c	  for a given variable.

	  call getwrd(ios,str)		! kt
	  if(ios.ne.0) then
	    write(stderr,'(2a,i4)') myname,
     &	      ': unexpected error from getwrd(), ios =',ios
	    istat=1
	    return	! bad news to the parent
	  endif

	  val=str2rn(str,ios)
	  if(ios.ne.0) then
	    l=max(lnblnk(str),1)
	    write(stderr,'(4a)') myname,
     &	      ': expecting a number, but "',str(1:l),'"'
	    istat=1
	    return	! bad news to the parent
	  endif

	  if(NotAnInt(val)) then
	    l=max(lnblnk(str),1)
	    write(stderr,'(4a)') myname,
     &		': expecting an integer, but "',str(1:l),'"'
	    istat=1
	    return	! bad news to the parent
	  endif

	  ikt=nint(val)
	  if(ikt.le.0 .or. ikt.gt.mxkt) then
	    l=max(lnblnk(str),1)
	    write(stderr,'(2a,i3,3a)') myname,
     &		': expecting kt value 1 to ',mxkt,', but "',
     &		str(1:l),'"'
	    istat=1
	    return	! bad news to the parent
	  endif

c	  ..Duplicate entry found
c	  =======================
	  if(ktname(ikt).ne.' ') then
	    write(stderr,'(2a,i3,2a)') myname,
     &	      ': kt =',ikt,' is already defined, ',ktname(ikt)
	    istat=1
	    return
	  endif

c  #2	  ..The second token is the name of the variable.  It should be
c ====    used in the program to relate kt number (#1) to a certain
c	  variable.

	  call getwrd(ios,str)
	  if(ios.ne.0) then
	    write(stderr,'(2a)') myname,
     &	      ': missing variable name'
	    istat=2
	    return	! bad news to the parent
	  endif
	  ktname(ikt)=str

c  #3	  ..The third token is the unit of the variable.  It may have
c ====	  not been used in the program, but for reference only.

	  call getwrd(ios,str)
	  if(ios.ne.0) then
	    write(stderr,'(2a)') myname,
     &	      ': missing variable unit'
	    istat=3
	    return	! bad news to the parent
	  endif
	  ktunit(ikt)=str

c  #4	  ..The fourth token is a description of the variable.  It is
c ====	  expected to be a "$" terminated string, or for the future,
c	  maybe a quoted string.

	  call getstr(ios,str)
	  if(ios.ne.0) str='?'
	  ktdesc(ikt)=str

c  #5	  ..Token #5 and beyond are multi-variate flags.  Each entry
c ====	  may be "0" or "1", or "F" or "T" in either lower or upper
c	  cases.  Spaces between entries are allowed.  The program will
c	  keep reading until upto nkt entries are read.

	  call getwrd(ios,str)
	  jkt=0			! for jkt=1,ikt
	  do while(ios.eq.0.and.jkt.lt.ikt)
	    ls=lnblnk(str)	! if ios=0, ls>0
	    ln=min(ls,ikt-jkt)	! only need upto the first ikt entries
	    do k=1,ln
	      jkt=jkt+1
	      if(str(k:k).eq.'0'.or.
     &		  str(k:k).eq.'F'.or.str(k:k).eq.'f') then
		ktmvar(ikt,jkt)=.false.
		ktmvar(jkt,ikt)=.false.
	      elseif(str(k:k).eq.'1'.or.
     &		  str(k:k).eq.'T'.or.str(k:k).eq.'t') then
		ktmvar(ikt,jkt)=.true.
		ktmvar(jkt,ikt)=.true.
	      else
		write(stderr,'(4a)') myname,
     &		  ': invalid multi-variate flag, "',str(1:ls),'"'
		istat=5
		return
	      endif
	    end do
	    call getwrd(ios,str)
	  end do

	  if(jkt.lt.ikt-1) then
	    write(stderr,'(2a,2(i2,a))') myname,
     &	      ': incomplete multi-variate flag, (jkt=',jkt,
     &	      ') < (ikt=',ikt,')'
	    istat=5
	    return

	  elseif(jkt.eq.ikt-1) then
	    ktmvar(ikt,ikt)=.true.

	  else if(jkt.ge.ikt) then	! use ktmvar(nkt,nkt) as a check
	    if(.not.ktmvar(ikt,ikt)) then
	      write(stderr,'(2a,l2)') myname,
     &		': invalid diagonal flag value, ',ktmvar(ikt,ikt)
	      istat=5
	      return
	    endif
	  endif

	  call rdnext(ios)
	end do

	end
