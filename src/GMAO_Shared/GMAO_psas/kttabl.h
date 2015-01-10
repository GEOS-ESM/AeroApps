!	(kttabl.h)

!	.."ktmax.h" is required before this, to define the value of
!	ktmax.

!	..Name of the table
!	===================
	character*15 RC_kt
	parameter   (RC_kt='DataTypeTable::')

	character*8  ktname	! names, [ktmax]
	character*8  ktunit	! units, [ktmax]
	character*32 ktdesc	! descriptions, [ktmax]
	logical      ktmvar	! multi-variation flags, [ktmax,ktmax]

	common/kttabli0/ktmvar(ktmax,ktmax)
	common/kttablc0/ktname(ktmax),ktunit(ktmax),ktdesc(ktmax)

!.
