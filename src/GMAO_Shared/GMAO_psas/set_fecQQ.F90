subroutine set_fecQQ
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: set_fecQQ - initialize matrix tables hfecQQ and vfecQQ
!
! !INTERFACE: (to do)
!
! !DESCRIPTION:
!	set_fecQQ() compute matrix tables hfecQQ and vfecQQ, based on
!	the input in hfecQ_tbl module and vfecQ_tbl module
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
!#define	_TRACE	!

#ifndef _REAL_
# define	_REAL_	real*8
#endif

use config, only : MXpar_hc, lvmax_vc
use const,  only : R_EARTH
use m_stdio,only : stdout
use rlev_imat,only : MXveclev,nveclev,pveclev
use hfecQ_tbl
use vfecQ_tbl
implicit none

include "hfecQQ.h"
include "vfecQQ.h"

character(len=*), parameter	:: myname='set_fecQQ'
real,parameter :: rade = .001*R_EARTH	! convert from m to km

real dlim
real QQmx1,QQmx2
real QQinc1,QQinc2

_REAL_ ctaus(nQQtab)
integer lv

  QQmx2=1.-cos(6030./rade)	! an earlier maximum distance in tau.
  QQmx1=.036*QQmx2		! keep the same HHmx2, so it becomes
				! function model independent
	! dlim is where the correlation vanishes

  dlim=0.
  do lv=1,nlev_hfecQ
    if(pars_hfecQ(1,lv).gt.dlim) dlim=pars_hfecQ(1,lv)
  end do
  QQmx2=1.-cos(dlim/rade)	! the maximum distance in tau
  QQmx2=max(3.*QQmx1,QQmx2)	! limit the minimum value in tau

  QQinc1=QQmx1/nQQtb1
  QQinc2=(QQmx2-QQmx1)/nQQtb2

#ifdef	_TRACE
	_TRACE	write(stdout,'(2a)') myname,': calling intp_ctaus()'
#endif

  call intp_ctaus(nQQtb1,QQinc1,nQQtb2,QQinc2,nQQtab,ctaus)

#ifdef	_TRACE
	_TRACE	write(stdout,'(2a,1p,2e10.3)') myname,		&
	_TRACE	  ':3 ',QQinc1,QQinc2
#endif

  QQbeg2=ctaus(nQQtb1+1)
  Qcoslim=1.-ctaus(nQQtab)

  qxQtb1=1./QQinc1
  qxQtb2=1./QQinc2
  
#ifdef	_TRACE
	_TRACE	write(stdout,'(2a)') myname,': calling intp_hCor()'
#endif

  call intp_hCor( name_hfecQ,				&
		  type_hfecQ,nlev_hfecQ,plev_hfecQ,	&
		  MXpar_hc,  npar_hfecQ,pars_hfecQ,	&
		  nveclev,pveclev, nQQtab, ctaus,	&
		  0,MXveclev,hfecQQ,hfecQQ,hfecQQ	)

	! normalization for non-separable form

  hfecQQ(1:nveclev,1:nQQtab)=.5*hfecQQ(1:nveclev,1:nQQtab)

#ifdef	_TRACE
	_TRACE	write(stdout,'(2a)') myname,': calling intp_vCor()'
#endif

  call intp_vCor( 'vfecQQ',					&
		  lvmax_vc,nlev_vfecQ,plev_vfecQ,corr_vfecQ,	&
		  MXveclev,nveclev,pveclev, vfecQQ		)

end
!.
