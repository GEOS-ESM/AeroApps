!	(vfecQ.h)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: vfecQ.h
!
! !SYNOPSIS:	#include "vfecQ.h"
!
! !DESCRIPTION: "vfecQ.h" defines parameters and global storages for
!		vertical forecast-<Q,Q>-error-correlation function
!		tables.
!
!		"lvmax.h" is required to define the value of lvmax, and
!		"mbanks.h" is required to define the value of mbanks.
!
! !EXAMPLES:
!
! !BUGS:
!
! !SEE ALSO:
!
! !REVISION HISTORY:
!  18Sep95  - Jing G.	- Created multiple copies of statistical tables
!			- Added the prolog
!  24Mar95  - Jing G.	- created for text based data table
!_______________________________________________________________________

!-	character*17 RC_vfecQ
!-	parameter   (RC_vfecQ='FcstErr*vCor_QQ::')

!	..Vertical Q Forecast-Error correlation table
!-	real		 vfecQ,vfecQm
!-	common/vfecQccc/vfecQ(0:mbanks-1,lvmax,lvmax),
!-    &			vfecQm(0:mbanks-1,MXveclev,MXveclev)

	real		vfecQQ
	common/vfecQQ0/ vfecQQ(MXveclev,MXveclev)
!.
