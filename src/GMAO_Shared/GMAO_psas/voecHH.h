!	(voecH.h)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: voecH.h
!
! !SYNOPSIS:	#include "voecH.h"
!
! !DESCRIPTION: "voecH.h" defines parameters and global storages for
!		vertical observationa;-<H,H>-error-correlation function
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

!-	integer      MxvcHH
!-	parameter   (MxvcHH=4)

!-	character*18 RC_voecH1
!-	character*18 RC_voecH2
!-	character*18 RC_voecH3
!-	character*18 RC_voecH4
!-	parameter   (RC_voecH1='ObsErr*vCor_HH-1::')
!-	parameter   (RC_voecH2='ObsErr*vCor_HH-2::')
!-	parameter   (RC_voecH3='ObsErr*vCor_HH-3::')
!-	parameter   (RC_voecH4='ObsErr*vCor_HH-4::')

!-	character*18 rc_voecH

!	..Vertical H observation error correlation table

!-	real		 voecH,voecHm
!-	common/voecHcom/voecH(0:mbanks-1,lvmax,lvmax,MxvcHH),
!-    &			voecHm(0:mbanks-1,MXveclev,MXveclev,MxvcHH)

	real		 voecHH
	common/voecHH0/ voecHH(MXveclev,MXveclev,MX_voecH)
!.
