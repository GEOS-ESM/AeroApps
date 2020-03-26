!	(ktmax.h)
!		- Total number of data type entries possible.

!	24Mar95 - Jing G.	- New version for text based data table.
!                                                     02/05/93 - ktmax.h
	integer      ktmax
	parameter  ( ktmax = 7 )

!.......................................................................

	integer      ktus,ktvs,ktslp		! surface variables
	integer	     ktHH,ktuu,ktvv,ktqq	! upper-air varibles
	parameter  ( ktus  = 1 )
	parameter  ( ktvs  = 2 )
	parameter  ( ktslp = 3 )
	parameter  ( ktuu  = 4 )
	parameter  ( ktvv  = 5 )
	parameter  ( ktHH  = 6 )
	parameter  ( ktqq  = 7 )

!.
