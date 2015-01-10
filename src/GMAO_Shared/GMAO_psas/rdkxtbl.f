	subroutine rdkxtbl(rc_kx,mxkx,
     &	  kxclas,kxrank,kxdesc,istat)

	use m_inpak90, only : lablin,rdnext,getwrd,str2rn,getstr
	use m_stdio,   only : stderr

	implicit none

c  Inputs
c ========
	character*(*) rc_kx	! label of the table
	integer mxkx		! size of the storage arrays

c  Outputs
c =========
	character*(*) kxclas(mxkx)	! error class of the source
	integer kxrank(mxkx)		! rank of the source
	character*(*) kxdesc(mxkx)	! data source descriptions
	integer istat		! input status code

c  Local Parameters and Variables
c ================================

	character*7 myname
	parameter(myname='rdkxtbl')

	include "realvals.h"

c	..Locals
	character*64 str
	integer ios
	integer l,ikx,irank
	real val

	integer lnblnk
	external lnblnk

	real a
	logical NotAnInt
	NotAnInt(a)=abs(a-nint(a)).gt.rfrcval

c	..Locate the table and check the status
c	========================================
	call lablin(trim(rc_kx))
	call rdnext(ios)
	if(ios.eq.2) then
	  write(stderr,'(3a)') myname,
     &	    ': table not found, ',trim(rc_kx)
	  istat=-1
	  return	! bad news to the parent
	endif

	istat=0		! if normal return

	do ikx=1,mxkx
	  kxclas(ikx)='-'
	  kxrank(ikx)=-1
	  kxdesc(ikx)='<undefined>'
	end do

c	..Processe the table line by line

	do while(ios.eq.0)

c	  ..Process the line token by token

c  #1	  ..The first token is the code number of the data source
c ====
	  call getwrd(ios,str)		! kx
	  if(ios.ne.0) then
	    write(stderr,'(2a,i4)') myname,
     &		': unexpected error from getwrd() for kx, ios =',ios
	    istat=1
	    return
	  endif

	  val=str2rn(str,ios)
	  if(ios.ne.0) then
	    l=max(lnblnk(str),1)
	    write(stderr,'(4a)') myname,
     &		': expecting a number for kx, but "',str(1:l),'"'
	    istat=1
	    return	! bad news to the parent
	  endif

	  if(NotAnInt(val)) then
	    l=max(lnblnk(str),1)
	    write(stderr,'(4a)') myname,
     &		': expecting an integer for kx, but "',str(1:l),'"'
	    istat=1
	    return	! bad news to the parent
	  endif

c	  ..Check the value
c	  =================
	  ikx=nint(val)
	  if(ikx.le.0 .or. ikx.gt.mxkx ) then
	    l=max(lnblnk(str),1)
	    write(stderr,'(2a,i3,3a)') myname,
     &		': expecting kx value 1 to ',mxkx,', but "',
     &		str(1:l),'"'
	    istat=1
	    return	! bad news to the parent
	  endif

c	  ..Duplicate entry found
c	  =======================
	  if(kxclas(ikx).ne.'-') then
	    write(stderr,'(2a,i3,2a)') myname,
     &	      ': kx =',ikx,' is already defined, ',kxclas(ikx)
	    istat=1
	    return
	  endif

c  #2	  ..The second token is the ObsErr class name of the source.
c ====
	  call getwrd(ios,str)
	  if(ios.ne.0) then
	    write(stderr,'(2a)') myname,
     &	      ': missing ObsErr class name'
	    istat=2
	    return	! bad news to the parent
	  endif
	  kxclas(ikx)=str

c  #3	  ..The third token is rank of the data source.
c ====
	  call getwrd(ios,str)
	  if(ios.ne.0) then
	    write(stderr,'(2a,i4)') myname,
     &		': unexpected error from getwrd() for kxrank, ios =',ios
	    istat=3
	    return
	  endif

	  val=str2rn(str,ios)
	  if(ios.ne.0) then
	    l=max(lnblnk(str),1)
	    write(stderr,'(4a)') myname,
     &		': expecting a number for kxrank, but "',str(1:l),'"'
	    istat=3
	    return	! bad news to the parent
	  endif

	  if(NotAnInt(val)) then
	    l=max(lnblnk(str),1)
	    write(stderr,'(4a)') myname,
     &		': expecting an integer for kxrank, but "',str(1:l),'"'
	    istat=3
	    return	! bad news to the parent
	  endif

	  irank=nint(val)
	  kxrank(ikx)=irank

c  #4	  ..The fourth token is a description of the source.
c ====
	  call getstr(ios,str)
	  if(ios.ne.0) str='?'
	  kxdesc(ikx)=str

	  call rdnext(ios)
	end do

	end
