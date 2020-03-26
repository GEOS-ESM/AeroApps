!	(kxmax.h)
!		- Total number of data sourcse entries possible.

!                                                     02/05/93 - kxmax.h
!	24Mar95 - Jing G.	- New version for text based data table.
!       13Jul98 - Genia Brin    - changed kxmax from 113 to 145 to
!                                 accomodate Dao_tovs
!	10Feb99 - Jing Guo	- Extended to kxmax=300.  The issue is
!		a different approach must be implemented, such as use
!		an internally defined index to a given kx number.
!	18Aug00 - R. Todling	- Changed max kx to 512
 
	integer   kxmax
	parameter(kxmax= 512)
	integer   kxmod
	parameter(kxmod=1000)	! the smallest 10's power > kxmax
!.
