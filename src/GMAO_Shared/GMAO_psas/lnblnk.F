#ifdef sgi
	integer function dummy_lnblnk()
	  dummy_lnblnk=-1
	end
#else
	integer function lnblnk(line)

c"     ..Find the real length of the string "line".  the real length
c"     is defined as the position of the last non-blank character in
c"     the string.  a non-blank character is any character other than
c"     a space(ascii 32), a null(ascii 0), or a tab(ascii 9).

	implicit none
	character*(*) line

	integer lc,ic,ich
	integer spc,nul,tab
	parameter(spc=32,nul=0,tab=9)
	logical blank

	lc=len(line)

	ic=lc
	blank=ic.gt.0
	if(blank) then
	   ich=ichar(line(ic:ic))
	   blank=ich.eq.spc.or.ich.eq.nul.or.ich.eq.tab
	endif
	do while(blank)   
	   ic=ic-1
	   blank=ic.gt.0
	   if(blank) then
	      ich=ichar(line(ic:ic))
	      blank=ich.eq.spc.or.ich.eq.nul.or.ich.eq.tab
	   endif
	end do

	lnblnk=ic
	end
#endif
