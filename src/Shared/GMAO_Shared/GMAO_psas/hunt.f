	subroutine hunt(x,n,ary,lm)
	use m_die,only : die
	use m_stdio,only : stderr

c	@(#) hunt(x,n,ary,lm): hunt the location of x in sorted ary.

c	..Hunt the location of x in ary.  The values in array "ary"
c	are assumed to be sorted in either ascending or descending
c	order.
c	  At the input, if lm is in the range from 1 to n, it is assumed
c	to be a guessed value and a hunting process will be started at
c	lm.  If lm is a value less than 1 or greater than n, a binary
c	search will be started.
c	  At the return, lm is from 1 to n-1, if x=ary(lm) or x is
c	between ary(lm) and ary(lm+1);  lm is n, if x=ary(n);  lm=0
c	if x is outside of ary(1), or lm=n+1, if x is outside of ary(n).

	implicit none
	character*4 myname
	parameter(myname='HUNT')

c	..Arguments
	real x	! a value which location in a list is to be found
	integer n	! size of the list
	real ary(n)	! a list of values
	integer lm	! Lower index where number x may locate

c	..Local variables
	integer jsm,jlg,jm
	logical guessed,ascnd,done
	integer inc

	if(n.le.0) then
	  write(stderr,'(2a,i2)') myname,
     &	    ': Zero sized array "ary", n =',n
	  call die(myname)
	endif

c	..Checking the sorted order
	ascnd=ary(n).ge.ary(1)

c	..Hunt or binary search
	guessed=lm.ge.1.and.lm.le.n

c	..Initialization
	if(ascnd) then
	  inc=1
	  jsm=1
	  jlg=n
	else
	  inc=-1
	  jsm=n
	  jlg=1
	endif

c	..Check with the boundaries first
	if(x.lt.ary(jsm)) then
	  lm=jsm-inc
	  return
	elseif(x.gt.ary(jlg)) then
	  lm=jlg+inc
	  return
	endif

c	..If a guessed value is given, start at that location.
	if(guessed) then
	  jsm=lm
	  jlg=lm
	  jm=jsm
	  done=.false.
	  do while(.not.done)
	    if(x.lt.ary(jm)) then
	      jlg=jm
	      jsm=jsm-inc
	      jm=jsm
	    else
	      jsm=jm
	      jlg=jlg+inc
	      jm=jlg
	    endif
	    done=jm.lt.1.or.jm.gt.n
	    if(.not.done) done=x.ge.ary(jsm).and.x.le.ary(jlg)
	    inc=inc+inc
	  end do
	  jsm=max(1,min(jsm,n))
	  jlg=max(1,min(jlg,n))
	endif

c	..Binary search
	do while(jlg-jsm.lt.-1.or.jlg-jsm.gt.1)
	  jm=(jsm+jlg)/2
	  if(x.lt.ary(jm)) then
	    jlg=jm
	  else
	    jsm=jm
	  endif
	end do

c	..Pick the value
	if(x.eq.ary(jlg)) then
	  lm=jlg
	elseif(x.eq.ary(jsm)) then
	  lm=jsm
	else
	  lm=min(jsm,jlg)
	endif
	end
