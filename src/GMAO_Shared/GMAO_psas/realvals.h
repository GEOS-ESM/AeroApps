!***********************************************************************
!	(realvals.h)
!***********************************************************************
!
! WRITTEN BY:	Jing Guo
!
!!DESCRIPTION:	Defining a set of practical real values for real number
!		precision checking.  The values are practical for real*4
!		representations on most operating systems.
!
! CALLED FROM:	n/a
!
! SYSTEM ROUTINES USED: n/a
!
! SUBROUTINES CALLED:	n/a
!
!!INPUT PARAMETERS:	n/a
!
!!OUTPUT PARAMETERS:	n/a
!
! FILES USED:	n/a
!
!!REVISION HISTORY:
!	30Dec94	- Jing G.	- Initial code
!
!***********************************************************************

! ..A big real*4 value for most machines and applications.

	real		 rbigval
	parameter	(rbigval=1.e+30)

! ..A small real*4 value for most machines and applications.

	real		 rsmlval
	parameter	(rsmlval=1.e-30)

! ..A siginificant real*4 value relative to 1, for most machines and
! applications.  I did not use 1.e-6, since a 16-based old real*4 did
! not have that kind of resolution.

	real		 rfrcval
	parameter	(rfrcval=1.e-7)

!	(end realvals.h)
