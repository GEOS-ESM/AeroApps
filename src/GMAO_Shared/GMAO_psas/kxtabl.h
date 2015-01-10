!	(kxtabl.h)

!	.."kxmax.h" is required before this, to define the value of
!	kxmax.

!	..Name of the table
!	===================
	character*17 RC_kx
	parameter   (RC_kx='DataSourceTable::')

	integer      kxrank	! quality ranks, [kxmax]
	integer      kxtmask	! mask data types of data sources
	integer      i_kxclas	! indices to module OEclass_tbl
	integer      i_hoecHc	! indices to module OEhcor_tbl
	integer      i_voecHc	! indices to module OEvcor_tbl
	integer      i_voecHu	! indices to module OEvcor_tbl

	character*8  kxclas	! ObsErr classes, [kxmax]
	character*32 kxdesc	! descriptions, [kxmax]

	common/kxtablc0/kxclas(kxmax),kxdesc(kxmax)

	common/kxtabli0/kxrank(kxmax),kxtmask(kxmax,ktmax)
	common/kxtabli1/i_kxclas(kxmax)
	common/kxtabli2/i_hoecHc(kxmax)
	common/kxtabli3/i_voecHc(kxmax),i_voecHu(kxmax)

! the same common block name is used in several declarations to avoid
! using continuation lines that won't be compatible with both free/F90
! and fixed/F77 formats.
!.
