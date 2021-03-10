!	(hoecHH.h)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: hoecHH.h
!
! !SYNOPSIS:	include "hoecHH.h"
!
! !DESCRIPTION: "hoecHH.h" defines the indirect factorization matrices of
!	the height-height horizontal observation error correlation
!	matrices.
!
! !EXAMPLES:
!
! !BUGS:
!
! !SEE ALSO:
!
! !REVISION HISTORY:
! 	18Sep95 - J. Guo	- Created multiple copies look-up tables
!				- Added the prolog
!	24Mar95 - Jing G.	- created for text based data table
!  5/10/94 modified to create table a fct of cos(sep angle) rather
!          than fct of distance
!_______________________________________________________________________

	integer   nHHotab	! table size (ntab_oHH)
	parameter(nHHotab = 1800 +1)
	
	integer   nOOtb1	! (ntb1_oHH)
	parameter(nOOtb1=nHHotab/2)
	integer   nOOtb2	! (ntb2_oHH)
	parameter(nOOtb2=nHHotab-1-nOOtb1)

	real hoecHH		! the packed matrix
	real Ocoslim,OObeg2	! tau section limits (tx2_oHH, tx1_oHH)
	real qxOtb1,qxOtb2	! tau resolutions (qx1_oHH, qx2_oHH)

	common /hoecHH0/hoecHH(MXveclev,nHHotab,MX_hoecH)
	common /hoecHH1/Ocoslim(MX_hoecH),OObeg2(MX_hoecH)
	common /hoecHH2/qxOtb1(MX_hoecH),qxOtb2(MX_hoecH)
!.
