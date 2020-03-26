	integer function lstins(ne,list,item,case)

c	..Searching a item in the list which is the same as item.
c	Trailing spaces are ignored.

	implicit none

	integer ne	! the size of the list
	character*(*) list(ne)	! the list
	character*(*) item	! the item to be found
	logical case		! if case sensitive
c-----------------------------------------------------------------------
c	..Locals
	integer l
	logical found

	integer lnblnk
	external lnblnk

	logical match
	integer li,ll
	integer icshft
	character*1 ch

c-----------------------------------------------------------------------
	icshft=ichar('a')-ichar('A')
	li=lnblnk(item)

	l=0
	found=.false.
	do while(l.lt.ne.and..not.found)
	  l=l+1

	  ll=lnblnk(list(l))

	  match=li.eq.ll
	  do while(match.and.ll.gt.0)

	    ch=list(l)(ll:ll)
	    if(.not.case) then
              if(ch.ge.'A'.and.ch.le.'Z') then
                ch=char(ichar(ch)+icshft)       ! set ch to lowercase
              elseif(ch.ge.'a'.and.ch.le.'z') then
                ch=char(ichar(ch)-icshft)       ! set ch to uppercase
              endif
            endif

	    match=item(ll:ll).eq.ch .or. item(ll:ll).eq.list(l)(ll:ll)
	    ll=ll-1
	  end do

	  found=match
	end do
	if(.not.found) l=0

	lstins=l
	end
