!	(levtabl.h)

!	.."lvmax.h" is required before this, to define the value of
!	lvmax.
!
!	..Plev_oe array defined in here is expected to be shared with
!	both voec tables and vfec tables, as while as ObsErr tables for
!	different instrument classes.

!-	..Name of the table
!-	===================
!-	character*14 RC_plev
!-	parameter   (RC_plev='ObsErr*Levels:')

!-	integer      nlev_oe	! number of levels, <=lvmax
!-	real         plev_oe	! pressure levels (for small tables)
	integer      nveclev
	real	     pveclev	! pressure levels (for LARGE tables)

!-	common/levtabl/nlev_oe,plev_oe(lvmax),nveclev,pveclev(MXveclev)
	common/levtabl/nveclev,pveclev(MXveclev)
!.
