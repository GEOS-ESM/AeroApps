	subroutine setpindx(n,p,w,nlev,plev)
	implicit none

c	.."pindx" a given array of pressure levels against a sorted list
c	of pressure referenc values.

c	  A "pindx" is defined as a value with its integer part for the
c	index to the reference array, and its decimal part for the log-
c	linear weight off the indexed reference array element.

c  Arguments
c ===========
	integer n	! size of p[], input
	real p(n)	! pressure levels to be "pindx"ed, input
	real w(n)	! "pindx" values for the levels, output

	integer nlev	! size of plev[], input
	real plev(nlev)	! reference pressure levels, input

c  Locals
c ========
	integer l,i

	l=0
	do i=1,n
	  call hunt(p(i),nlev,plev,l)
	  if(l.lt.1) then
	    w(i)=1.
	  elseif(l.ge.nlev) then
	    w(i)=nlev
	  else
	    w(i)=l+log(p(i)/plev(l))/log(plev(l+1)/plev(l))
	  endif
	end do

	end
