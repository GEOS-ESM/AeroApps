subroutine setpix(nobs,kx,rlat,rlon,pix)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS     !
!---------------------------------------------------------------------
!
! !ROUTINE: setpix - set profile indices of observations
!
! !INTERFACE: (to do)
!
! !DESCRIPTION:
!	setpix() sets profile indices of observations according to their
!	data sources (kx) and locations (rlat & rlon).  The defination
!	used here may be different from ks (sounding index) used in ODS.
!
! !EXAMPLES: (to do)
!
! !BUGS:
!	should be modified to be consistent to ODS when the defination
!	of ks is final in ODS.
!
! !SEE ALSO: (to do)
!
! !REVISION HISTORY:
! 	06Dec95 - J. Guo	- initial coding
!_____________________________________________________________________
implicit none
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  integer,intent(in)	:: nobs
  integer,intent(in)	:: kx(nobs)
  real,   intent(in)	:: rlat(nobs)
  real,   intent(in)	:: rlon(nobs)

  integer,intent(out)	:: pix(nobs)
!_______________________________________________________________________

  include "kxmax.h"

  integer kx_i,i
  integer knt(kxmax)
  integer pix_now
  logical begin_prof
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (nobs.le.0) return

  knt(1:kxmax)=-1		! none yet

	! this part needs to be serially processed, since it is
	! it is upto the previous entry to determine the pix
	! value of the current entry.

  kx_i=kx(1)

  knt(kx_i)=knt(kx_i)+1
  pix_now=kx_i+knt(kx_i)*kxmod

  pix(1)=pix_now

  do i=2,nobs
    kx_i=kx(i)

    begin_prof=	kx_i.ne.kx(i-1)		.or.		&
		rlat(i).ne.rlat(i-1)	.or.		&
		rlon(i).ne.rlon(i-1)

    if(begin_prof) then		! If a new profile, count it
      knt(kx_i)=knt(kx_i)+1
      pix_now=kx_i+knt(kx_i)*kxmod
    endif

    pix(i)=pix_now
  end do
!_______________________________________________________________________
end
!.
