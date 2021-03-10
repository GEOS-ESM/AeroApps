	subroutine rdoetbl(rc_oetbl,mxtyp,ntyp,ktlist,
     &	  mxlev,loc_lev,len_lev,oetbl,ivcHHu,ivcHHc,ihcHHc,istat)

	use m_inpak90, only : lablin, rdnext, getwrd, str2rn
	use m_stdio,   only : stderr
	implicit none

c  Inputs
c ========
	character*(*) rc_oetbl
	integer mxtyp		! maximum,number of variable entries
	integer mxlev		! total number of levels for all types

c  Outputs
c =========
	integer ntyp		! number of entries
	character*(*) ktlist(mxtyp)

	integer loc_lev(mxtyp)
	integer len_lev(mxtyp)
	real oetbl(mxlev,mxtyp)	! all ObsErr (sigma-O) values

	integer ivcHHu
	integer ivcHHc,ihcHHc	! ObsErr Correlation functions of HH

	integer istat

c	..Locals
	character*7 myname
	parameter(myname='rdoetbl')

	real UNDEFoe		! Undefined value for a sigma-O
	parameter(UNDEFoe=-1.)	! Sigma-O !< 0


	integer lnblnk
	external lnblnk

	include "realvals.h"

	character*64 styp
	character*16 snum
	integer ios	! data buffer processing status
	real val	! tmp real value
	integer ival
	integer n	! type counter
	integer lv	! level counter
	integer lt,ln

	real a
	logical NotAnInt
	NotAnInt(a)=abs(a-nint(a)).gt.rfrcval

c  Table processing
c ==================
c	..Locate and verify the existence of the label.

	call lablin(trim(rc_oetbl))
	call rdnext(ios)
	if(ios.eq.2) then
	  write(stderr,'(3a)') myname,
     &	    ': table not found, ',rc_oetbl
	  istat=-1
	  return	! bad news to the parent
	endif


c	..Processe the table line by line

	ivcHHu=0
	ivcHHc=0
	ihcHHc=0

	istat=0		! if normal return

	n=0
	ntyp=0
	do while(ios.eq.0)	! read all entries

c  #1
c ====
	  call getwrd(ios,styp)
	  if(ios.ne.0) then
	    ln=max(lnblnk(styp),1)
	    write(stderr,'(4a)') myname,
     &	      ': unexpected error from getwrd(), "',styp(1:ln),'"'
	    istat=1
	    return
	  endif

c  #2:
c ====
c	  ..Correlation function class
c	  =============================
	  if(index(styp,'Cor_').gt.0) then

c	    ..Unknonw Correlation function class
	    if(	styp.ne.'vCor_HH'  .and.styp.ne.'hCor_HH'  .and.
     &		styp.ne.'vCor_HH.c'.and.styp.ne.'hCor_HH.c'.and.
     &		styp.ne.'vCor_HH.u' ) then
	      lt=max(lnblnk(styp),1)
	      write(stderr,'(4a)') myname,
     &		': unknown correlation type skipped, "',styp(1:lt),'"'

c	    ..If the name is known
	    else

	      call getwrd(ios,snum)
	      if(ios.ne.0) then
		lt=max(lnblnk(styp),1)
		write(stderr,'(4a)') myname,
     &		  ': missing an integer number, "',styp(1:lt),'"'
		istat=2
		return
	      endif

	      val=str2rn(snum,ios)
	      if(ios.ne.0) then
		lt=max(lnblnk(styp),1)
		ln=max(lnblnk(snum),1)
		write(stderr,'(6a)') myname,
     &		  ': expecting a number for "',
     &		  styp(1:lt),'", but "',snum(1:ln),'"'
		istat=2
		return
	      endif

	      if(NotAnInt(val)) then
		lt=max(lnblnk(styp),1)
		ln=max(lnblnk(snum),1)
		write(stderr,'(6a)') myname,
     &		  ': expecting an integer for "',
     &		  styp(1:lt),'", but "',snum(1:ln),'"'
		istat=2
		return
	      endif

	      ival=nint(val)
	      if(styp.eq.'vCor_HH'.or.styp.eq.'vCor_HH.c') then
		ivcHHc=ival
	      elseif(styp.eq.'hCor_HH'.or.styp.eq.'hCor_HH.c') then
		ihcHHc=ival
	      elseif(styp.eq.'vCor_HH.u') then
		ivcHHu=ival
	      endif
	    endif	! a known/unknown correlation function type

c	  ..Variable name for a sigma-O list
c	  ==================================
	  else
	    n=n+1
	    ntyp=n
	    if(n.le.mxtyp) then
	      ktlist(n)=styp
	      loc_lev(n)=0
	      len_lev(n)=0
	      lv=0

c	      ..Loop over token in an variable entry
c	      ======================================
	      call getwrd(ios,snum)
	      do while(ios.eq.0.and.lv.lt.mxlev)
		lv=lv+1

		if(snum.eq.'-') then
		  oetbl(lv,n)=UNDEFoe
		else
		  val=str2rn(snum,ios)
		  if(ios.eq.0 .and. val .lt. 0.) ios=1
		  if(ios.ne.0) then
		    lt=max(lnblnk(styp),1)
		    ln=max(lnblnk(snum),1)
		    write(stderr,'(6a)') myname,
     &			': invalid ObsErr number, "',
     &			styp(1:lt),' ',snum(1:ln),'"'
		    istat=3
		    return
		  else
		    if(loc_lev(n).eq.0) loc_lev(n)=lv
		    len_lev(n)=len_lev(n)+1
		    oetbl(lv,n)=val
		  endif
		endif

		call getwrd(ios,snum)
	      end do

c	      ..If the table is already full, return with an error code.
	      if(ios.eq.0) then
		write(stderr,'(2a,i2,a)') myname,
     &		  ': buffer is full (mxlev=',mxlev,
     &		  '), but more in the table'
		istat=4
		return
	      endif
	    endif

	  endif

	  call rdnext(ios)
	end do

	end
