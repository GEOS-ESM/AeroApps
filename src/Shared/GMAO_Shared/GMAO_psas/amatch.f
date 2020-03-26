	logical function amatch(str,pat,case)
	implicit none
	character*(*) str
	character*(*) pat
	logical case	! case sensitive?

c	..Match a pattern letter by letter, except "~", but abbrievia-
c	tion is allowed.

c	The form of a pattern:

c		aaa~oooo

c	..Local vars.
	integer ls,lp
	integer i,is,ip,icdif
	integer la
	logical match
	character*1 c

c	..Functions
	integer lnblnk
	external lnblnk

	ls=lnblnk(str)
	lp=lnblnk(pat)

	icdif=ichar('a')-ichar('A')

	la=index(pat,'~')-1
	if(la.lt.0) la=lp		! if no "~"
	match=(la.eq.lp.and.ls.eq.lp).or.		! w/out "~"
     &	  (la.lt.lp.and.ls.ge.la.and.ls.le.lp-1)	! with "~"

c	..Matching the first part
	i=0
	do while(i.lt.la.and.match)
	  i=i+1
	  match=i.le.ls
	  if(match) then
	    if(.not.case) then
	      if(pat(i:i).ge.'A'.and.pat(i:i).le.'Z') then
		c=char(ichar(pat(i:i))+icdif)
	      elseif(pat(i:i).ge.'a'.and.pat(i:i).le.'z') then
		c=char(ichar(pat(i:i))-icdif)
	      else
		c=pat(i:i)
	      endif
	      match=str(i:i).eq.pat(i:i).or.str(i:i).eq.c
	    else
	      match=str(i:i).eq.pat(i:i)
	    endif
	  endif
	end do

c	..Matching the second part
	is=i+1
	ip=i+2	! skip "~"
	do while(is.lt.ls.and.match)
	  is=is+1
	  ip=ip+1
	  if(.not.case) then
	    if(pat(ip:ip).ge.'A'.and.pat(ip:ip).le.'Z') then
	      c=char(ichar(pat(ip:ip))+icdif)
	    elseif(pat(ip:ip).ge.'a'.and.pat(ip:ip).le.'z') then
	      c=char(ichar(pat(ip:ip))-icdif)
	    else
	      c=pat(ip:ip)
	    endif
	    match=str(is:is).eq.pat(ip:ip).or.str(is:is).eq.c
	  else
	    match=str(is:is).eq.pat(ip:ip)
	  endif
	end do

	amatch=match
	end
