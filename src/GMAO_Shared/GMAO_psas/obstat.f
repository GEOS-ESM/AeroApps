	subroutine obstat(lu,
     &	  nX,       kxX,    ktX, rlevX,
     &    Xvec,   nlevs,  plevs, header)
!-----------------------------------------------------------------------
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE:	obstat
!
! !DESCRIPTION:	print statistics of a PSAS innovation vector type data
!	array.  The statistics include mean, standard-deviation, maximum
!	and minimum of array Xvec in an orderly sequence,
!	  variable type(kt) - levels(plevs, if nlevs>0).
!
! !CALLING SEQUENCE:
!	call obstat(lu, nX,kxX,ktX,rlatX,
!    &	  Xvec,nlevs,plevs,header)
!
! !INPUT PARAMETERS:
!	integer lu		! logical output unit
!			
!	integer nX		! size of data array
!	  integer   kxX(nX)	! data sources
!	  integer   ktX(nX)	! data (variable) types
!	  real	  rlevX(nX)	! levels
!	  real	   Xvec(nX)	! data array
!
!	integer nlevs		! number of levels to check
!					=0	all levels
!					>0	values in plevs (below)
!	  real	  plevs(*)	! (pressure) level values
!
!	character*(*) header	! header message

! !OUTPUT PARAMETERS: <none>
!
! !BUGS:	Assumed the size of the data array in a 5 digit number,
!	the size of the data array of a single variable type in a single
!	region is a 4 digit number.
!
! !SEE ALSO: yet to decide
!
! !SYSTEM ROUTINES:
!
!  external procedures:
!	integer lnblnk()
!
!  external data structures:
!	include "ktmax.h"
!	include "kttabl.h"
!	include "realvals.h"
!
! !FILES USED:	logical unit lu for output
!
! !REVISION HISTORY:
! 	10Aug95 - J. Guo	- added the prolog
!				- created the procedure too
!
!-----------------------------------------------------------------------
	implicit none
	include "ktmax.h"	! ktmax :=

	integer lu		! logical output unit

c	..Data attributes
	integer nX		! size of data array
	  integer   kxX(nX)	! data sources
	  integer   ktX(nX)	! data (variable) types
	  real	  rlevX(nX)	! levels
	  real	   Xvec(nX)	! data array

c	..Table entry controls
	integer nlevs		! number of levels to check
	  real	  plevs(*)	! (pressure) level values

	character*(*) header	! header message

c	..Local variables

	integer kt,lv,ix,ln
	integer imx,imn,nav,ntt
	real xmx,xmn,xav,xdv,del,xrf,avlev
	logical fitin
	character*4 apw(ktmax)
	real amx(ktmax)
	integer nkt(ktmax)
	integer ipw
	real amag

c	..Local setting and functions

	include "realvals.h"
	include "kttabl.h"

	character*6 myname
	parameter(myname='obstat')

	integer lnblnk
	external lnblnk

	logical between
	real x1,x2,x
	between(x1,x2, x)=x.gt.min(x1,x2).and.x.lt.max(x1,x2)

c	..Write a header of the table

	ln=lnblnk(header)
	if(ln.eq.0) then
	  write(lu,'(/a)') myname,'::'
	else
	  write(lu,'(/2a)') header(1:ln),'::'
	endif
	write(lu,'(a,4x,a,5x,3(a,3x),a,4x,3(a,3x),a)')
     &	  'kt','10^','n','pres',
     &	  'mean+/-stdv','mini(imin)','maxi(imax)'

c	..Check scales of each variable types

	do kt=1,ktmax
	  amx(kt)=0.
	  nkt(kt)=0
	end do

	do ix=1,nX	! all data
	  kt=ktX(ix)
	  if(kt.ge.1.and.kt.le.ktmax) then
	    nkt(kt)=nkt(kt)+1
	    amag=abs(Xvec(ix))
	    if(amx(kt).lt.amag) amx(kt)=amag
	  endif
	end do

c	..Estimate the magnitudes for all types
	do kt=1,ktmax

c	  ..If there are some data and the value is reasonablely large
	  if(nkt(kt).gt.0.and.amx(kt).gt.rsmlval) then
	    amag=log(amx(kt))/log(10.)
	    ipw=int(amag)
	    if(ipw.gt.amag) ipw=ipw-1
	  else
	    ipw=0
	  endif
	  amx(kt)=10.**real(ipw)
	  write(apw(kt),'(a,sp,i3.2)') 'E',ipw
	end do

c	..Extimate statistics for each type-levels

	do kt=1,ktmax	! types

	  ntt=0
	  if(kt.eq.ktus.or.kt.eq.ktvs.or.kt.eq.ktslp .or.
     &	    nlevs.eq.0 ) then	! no level difference

	    nav=0
	    do ix=1,nX
	      if(kt.eq.ktX(ix)) then

	        nav=nav+1
	        if(nav.eq.1) then
	          xav=0.
	          xdv=0.
	          avlev=0.

	          xrf=Xvec(ix)

		  imn=ix
		  xmn=Xvec(ix)
		  imx=ix
		  xmx=Xvec(ix)
		endif

		if(xmn.gt.Xvec(ix)) then
		  imn=ix
		  xmn=Xvec(ix)
		endif

		if(xmx.lt.Xvec(ix)) then
		  imx=ix
		  xmx=Xvec(ix)
		endif

		del=Xvec(ix)-xrf
		xav=xav+del
		xdv=xdv+del*del
		avlev=avlev+rlevX(ix)
	      endif	! kt.eq.ktX(ix)
	    end do	! ix=1,nX

	    if(nav.gt.0) then
	      ntt=ntt+nav

	      xav=xav/nav
	      del=xdv/nav-xav*xav
	      xdv=sign(sqrt(abs(del)),del)	! permit(!) del<0
	      xav=xav+xrf

	      avlev=avlev/nav

	      write(lu,'(a8,x,a4,x,i4,x,f6.1,x,'//
     &		'f6.3,x,f5.3,2(x,f6.3,x,i5))')
     &		ktname(kt),apw(kt),nav,avlev,
     &		xav/amx(kt),xdv/amx(kt),
     &		xmn/amx(kt),imn,xmx/amx(kt),imx
	    endif

	  else		! nlevs>0

	    do lv=0,nlevs	! nlevs > 0

	      nav=0
	      do ix=1,nX

		if(kt.eq.ktX(ix)) then

c		  ..beyond level 1
		  if(lv.eq.0) then

		    fitin=plevs(1).eq.plevs(nlevs) .and.
     &		      rlevX(ix).gt.plevs(1)
		    fitin=fitin .or.
     &		      between(rlevX(ix),plevs(nlevs), plevs(1))

c		  ..beyond level nlevs, including the level
		  elseif(lv.eq.nlevs) then

		    fitin=plevs(1).eq.plevs(nlevs) .and.
     &		      rlevX(ix).le.plevs(nlevs)
		    fitin=fitin .or. rlevX(ix).eq.plevs(lv) .or.
     &		      between(rlevX(ix),plevs(1), plevs(nlevs))

c		  ..a layer between two levels, including one end.
		  else

		    fitin=rlevX(ix).eq.plevs(lv) .or.
     &		      between(plevs(lv),plevs(lv+1),rlevX(ix))

		  endif

		  if(fitin) then
		    nav=nav+1

		    if(nav.eq.1) then
		      xav=0.
		      xdv=0.
		      avlev=0.

		      xrf=Xvec(ix)

		      imn=ix
		      xmn=Xvec(ix)
		      imx=ix
		      xmx=Xvec(ix)
		    endif

		    if(xmx.lt.Xvec(ix)) then
		      imx=ix
		      xmx=Xvec(ix)
		    endif

		    if(xmn.gt.Xvec(ix)) then
		      imn=ix
		      xmn=Xvec(ix)
		    endif

		    del=Xvec(ix)-xrf	! for robust xdv
		    xav=xav+del
		    xdv=xdv+del*del

		    avlev=avlev+rlevX(ix)
		  endif	! fitin
		endif	! kt.eq.ktX(ix)
	      end do	! ix=1,nX

	      if(nav.gt.0) then
		ntt=ntt+nav

		xav=xav/nav
		del=xdv/nav-xav*xav
		xdv=sign(sqrt(abs(del)),del)
		xav=xav+xrf

		avlev=avlev/nav

		write(lu,'(a8,x,a4,x,i4,x,f6.1,x,'//
     &		  'f6.3,x,f5.3,2(x,f6.3,x,i5))')
     &		  ktname(kt),apw(kt),nav,avlev,
     &		  xav/amx(kt),xdv/amx(kt),
     &		  xmn/amx(kt),imn,xmx/amx(kt),imx

	      endif
	    end do	! lv=0,nlevs

	  endif	! kt.eq.ktsfc.or.nlevs.eq.0; else nlevs.gt.0

	  if(ntt.ne.nkt(kt)) write(lu,'(a,i6,a,i6)')
     &	    '# Total ',ntt,' out of ',nkt(kt)

	end do	! kt=1,ktmax

c	..Write a end-of-table record
	write(lu,'(a)') '::'

	end
