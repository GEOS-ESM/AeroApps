!	(ktwanted.h)
!-----------------------------------------------------------------------
!///////////////////////////////////////////////////////////////////////
!-----------------------------------------------------------------------
!	20Jan95 - Jing G.	- Created this block to store control
!				  parameters specifying analysis
!				  increment computation.

!	..Ktmax.h must be included before include this one.

	logical ktwanted
	common /ktcontrl/ ktwanted(ktmax)
