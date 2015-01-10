!#define _TRACE	!
subroutine set_oecHH
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: set_oecHH - initialize matrix tables hoecHH and voecHH
!
! !INTERFACE: (to do)
!
! !DESCRIPTION:
!	set_oecHH() compute matrix tables hoecHH and voecHH, based on
!	the input in hoecH_tbl module and voecH_tbl module
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
! 	12Jan96 - J. Guo	- programmed and added the prolog
!_______________________________________________________________________

#ifndef _REAL_
# define	_REAL_	real*8
#endif

use config, only : MXpar_hc, lvmax_vc, MX_hoecH, MX_voecH
use const,  only : R_EARTH
use m_stdio,only : stdout
use rlev_imat,only : MXveclev,nveclev,pveclev
use hoecH_tbl
use voecH_tbl
implicit none

include "hoecHH.h"
include "voecHH.h"

character(len=*), parameter	:: myname='set_oecHH'
real,parameter :: rade = .001*R_EARTH	! convert from m to km

integer i,lv
real corlen,dlim
real HHomx1,HHomx2,OOinc1,OOinc2
_REAL_ ctaus(nHHotab)

do i=1,n_hoecH	! for each function type

	! determine region 1
  corlen=0.
  do lv=1,nlev_hoecH(i)
	! a tempory fix
    corlen=corlen+pars_hoecH(2,lv,i)
  end do
  corlen=corlen/nlev_hoecH(i)
  HHomx1=1.-cos(2.*corlen/rade)
  OOinc1=HHomx1/nOOtb1
  
  	! determine region 2
  dlim=0.
  do lv=1,nlev_hoecH(i)
    if(pars_hoecH(1,lv,i).gt.dlim) dlim=pars_hoecH(1,lv,i)
  end do
  HHomx2=1.-cos(dlim/rade)
  if(HHomx2.le.HHomx1) HHomx2=1.-cos(8.*corlen/rade)
  OOinc2=(HHomx2-HHomx1)/nOOtb2

#ifdef	_TRACE
	_TRACE write(stdout,*) corlen,HHomx1,OOinc1
	_TRACE write(stdout,*) dlim,HHomx2,OOinc2
#endif

  call intp_ctaus(nOOtb1,OOinc1,nOOtb2,OOinc2,nHHotab,ctaus)

#ifdef	_TRACE
	_TRACE write(stdout,*) ctaus(1),ctaus(nOOtb1+1),ctaus(nHHotab)
#endif

  OObeg2(i)=ctaus(nOOtb1+1)
  Ocoslim(i)=1.-ctaus(nHHotab)

#ifdef	_TRACE
	_TRACE write(stdout,*) OObeg2(i),Ocoslim(i)
#endif

  qxOtb1(i)=1./OOinc1
  qxOtb2(i)=1./OOinc2

  call intp_hCor( name_hoecH(i),				&
		  type_hoecH(i),nlev_hoecH(i),plev_hoecH(1,i),	&
		  MXpar_hc,     npar_hoecH(i),pars_hoecH(1,1,i),&
		  nveclev,pveclev, nHHotab, ctaus,		&
		  0,MXveclev,hoecHH(1,1,i),hoecHH(1,1,i),	&
		  hoecHH(1,1,i)					)

#ifdef	_TRACE
	_TRACE	write(stdout,*) OObeg2(i),Ocoslim(i)

	_TRACE	write(stdout,'(2a,x,a,4i5,1p,10e10.2)')		&
	_TRACE	  ':: ',name_hoecH(i)(1:2),type_hoecH(i),	&
	_TRACE	  nlev_hoecH(i),nHHotab,nOOtb1,nOOtb2,		&
	_TRACE	  Ocoslim(i),OObeg2(i),qxOtb1(i),qxOtb2(i),	&
	_TRACE	  OOinc1,OOinc2,				&
	_TRACE	  corlen,HHomx1,dlim,HHomx2
#endif

	! normalization for non-separable form

  hoecHH(1:nveclev,1:nHHotab,i)=.5*hoecHH(1:nveclev,1:nHHotab,i)

end do

do i=1,n_voecH
  call intp_vCor( name_voecH(i),				&
		  lvmax_vc,nlev_voecH(i),plev_voecH(1,i),	&
		  corr_voecH(1,1,i),				&
		  MXveclev,nveclev,pveclev, voecHH(1,1,i)	)
end do

end
!.
