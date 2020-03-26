subroutine ll2qvec(	ndat,rlats,rlons,			&
			qr_x,qr_y,qr_z,qm_x,qm_y,qm_z,ql_x,ql_y	)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: ll2qvec - set unit vectors qr, qm, and ql
!
! !SYNOPSIS:
!
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!
! !EXAMPLES:
!
! !BUGS:
!
! !SEE ALSO:
!
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
! !REVISION HISTORY:
! 	08Dec95 - J. Guo	- initial code replacing ll2xyz() and
!				  qtrig0() etc.
!_______________________________________________________________________

implicit none

integer, intent(in)	:: ndat		! size of the arrays
real,    intent(in)	:: rlats(ndat)	! latitudes
real,    intent(in)	:: rlons(ndat)	! longitudes
real,    intent(out)	:: qr_x(ndat)	! x of qr
real,    intent(out)	:: qr_y(ndat)	! y of qr
real,    intent(out)	:: qr_z(ndat)	! z of qr
real,    intent(out)	:: qm_x(ndat)	! x of d(qr)/dm
real,    intent(out)	:: qm_y(ndat)	! y of d(qr)/dm
real,    intent(out)	:: qm_z(ndat)	! z of d(qr)/dm
real,    intent(out)	:: ql_x(ndat)	! x of d(qr)/dl
real,    intent(out)	:: ql_y(ndat)	! y of d(qr)/dl

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! locals
!========
  real    :: slat,clat,slon,clon
  real	  :: deg
  integer :: i

!-----------------------------------------------------------------------
	deg=4.*atan(1.)/180.
!-----------------------------------------------------------------------
do i=1,ndat
  slat=sin(rlats(i)*deg)
  clat=cos(rlats(i)*deg)
  slon=sin(rlons(i)*deg)
  clon=cos(rlons(i)*deg)

  qr_x(i)= clat*clon
  qr_y(i)= clat*slon
  qr_z(i)= slat
  qm_x(i)=-slat*clon
  qm_y(i)=-slat*slon
  qm_z(i)= clat
  ql_x(i)=-slon
  ql_y(i)= clon
end do
!_______________________________________________________________________
end
!.
