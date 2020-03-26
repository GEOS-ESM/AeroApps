!#define	_TRACE	!
!#define 	_DEBUG	!
subroutine intp_sigO(nobs,kxobs,ktobs,rlevs,sigC,sigU)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: intp_sigO - interpolate rms OE to observation locations
!
! !INTERFACE: (to do)
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
! 	21Dec95 - J. Guo	- programmed and added the prolog
!	24Sep96 - J. Guo	- add onep to avoid (real).eq.(real)
!_______________________________________________________________________
use config, only : kxmax	! Max number of kx values
use config, only : ktslp,ktus,ktvs
use m_stdio,only : stderr,stdout
use m_die,  only : die
use OEclass_tbl
implicit none
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  integer, intent(in)	:: nobs		! number of obs.
  integer, intent(in)	:: ktobs(nobs)	! data types
  integer, intent(in)	:: kxobs(nobs)	! data sources
  real,    intent(in)	:: rlevs(nobs)	! data pressure levels

  real,    intent(out)	:: sigC(nobs)	! horizontally correlated OE
  real,    intent(out)	:: sigU(nobs)	! horizontally uncorrelated OE
!_______________________________________________________________________

  include "kxtabl.h"		! contains common blocks

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Local workspace and variables
  integer,allocatable :: klev(:)	! index to plev_oe
  real   ,allocatable :: wlev(:)	! weight table of plev_oe

  integer kx,kc,kt,kl
  real    wt
  integer i,lc,ln,le,istat

#ifdef	_TRACE
	_TRACE	integer  ll
	_TRACE	integer  listvals
	_TRACE	external listvals
	_TRACE	character(len=132) line
#endif

!   Parameters

  character(len=*), parameter	:: myname='intp_sigO'
  logical,          parameter	:: nearest=.true.
  real,             parameter	:: onep = 1.+.01
!_______________________________________________________________________
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		! If no data to processing, either by the input or the
		! output.

  if(nobs.le.0) return

  if(nlev_oe.eq.0.or.nOEclas.eq.0) then
     if(nlev_oe.eq.0) write(stderr,'(2a)') myname,		&
	': yet to be defined, plev_oe'
     if(nOEclas.eq.0) write(stderr,'(2a)') myname,		&
	': yet to be defined, OEclas'
     call die(myname)
  endif

!-----------------------------------------------------------------------
	allocate(klev(nobs),wlev(nobs), stat=istat)
	if(istat.ne.0) then
	  write(stderr,'(2a,i5)') myname,		&
		': allocate() error, stat =',istat
	  call die(myname)
	endif

		! define the table indices and the weights.  Notice
		! the indices are not rounded off to the nearest tabled
		! levels (see slogtab), but those "earlier" levels.

#ifdef	_TRACE
	_TRACE	write(stdout,'(2a)') myname,': calling slogtab()'
	_TRACE	write(stdout,'(2a,i3)') myname,': ',nlev_oe
	_TRACE	if(nlev_oe.gt.0) then
	_TRACE	  ll=listvals(nlev_oe,plev_oe,line)
	_TRACE	  write(stdout,'(3a)') myname,': ',line(1:ll)
	_TRACE	endif
#endif

  call slogtab(.not.nearest,nlev_oe,plev_oe,nobs,rlevs,klev,wlev)

!		! verify if the data flow is as expected
!	_TRACE	write(stdout,'(2a)') myname,': loop verify i_kxclas'
!  do i=1,nobs
!    kx=kxobs(i)
!    kc=i_kxclas(kx)
!    if(kc.le.0) then
!      write(stdout,'(2a,i3)') myname,			&
!	': unclassified data source, kx=',kx
!      call die(myname)
!    endif
!  end do

		! loop over all observations for sigOc/sigOu
#ifdef	_TRACE
	_TRACE	write(stdout,'(2a)') myname,': loop over nobs'
#endif
  do i=1,nobs
    kx=kxobs(i)		! data source index
    kc=i_kxclas(kx)	! OE class index

    sigC(i)=-1.
    sigU(i)=-1.

	! If the instrument is supported in the user defined (psas.rc)
	! observatiton error covariance model, kc will be a positive
	! interger.

    if(kc.gt.0) then

      kt=ktobs(i)		! data type index

		! for this class (kc), the horizontally correlated OE
		! component

      lc=loc_sigOc(kt,kc)	! the level the sigO table starts
      ln=len_sigOc(kt,kc)	! the number of levels in sigO table
      le=lc+ln-1

		! This component may be not defined for this (kt,kx).
		! There is a possibility that the data type (kt) or its
		! correlated component is not valid for this instrument.

      if( lc.gt.0 .and. ln .gt.0 .and.		&
	  klev(i).ge.lc .and. klev(i).le.le	) then

	if(kt.eq.ktslp .or. kt.eq.ktus .or. kt.eq.ktvs) then

		! two-dimensional variable, only one piece of data

	  sigC(i)=sigOc(lc,kt,kc)

	elseif(rlevs(i).ge.min(plev_oe(lc),plev_oe(le))/onep .and. &
	       rlevs(i).le.max(plev_oe(lc),plev_oe(le))*onep	  ) then

		! three-dimensional variable, a log-linear interpolation
		! if the data is within the range of the table.

	  kl=klev(i)
	  sigC(i)=sigOc(kl,kt,kc)
	  if(kl.lt.le) then
	    wt=wlev(i)
	    sigC(i)=sigC(i) + wt*(sigOc(kl+1,kt,kc)-sigC(i))
	  endif
      
	endif	! kt.eq.ktslp .or. ...

      endif	! lc.gt.0 .and. ...

		! for this class (kc), the horizontally uncorrelated OE
		! component

      lc=loc_sigOu(kt,kc)	! the level the sigO table starts
      ln=len_sigOu(kt,kc)	! the number of levels in sigO table
      le=lc+ln-1


		! This component may be not defined for this (kt,kx).
		! There is a possibility that the data type (kt) or its
		! uncorrelated component is not valid for this
		! instrument.

      if( lc.gt.0 .and. ln .gt.0 .and.		&
	  klev(i).ge.lc .and. klev(i).le.le	) then

	if(kt.eq.ktslp .or. kt.eq.ktus .or. kt.eq.ktvs) then

		! two-dimensional variable, only one piece of data

	  sigU(i)=sigOu(lc,kt,kc)

	elseif(rlevs(i).ge.min(plev_oe(lc),plev_oe(le))/onep .and. &
 	       rlevs(i).le.max(plev_oe(lc),plev_oe(le))*onep	 ) then

		! three-dimensional variable, a log-linear interpolation
		! if the data is within the range of the table.

	  kl=klev(i)
	  sigU(i)=sigOu(kl,kt,kc)
	  if(kl.lt.le) then
	    wt=wlev(i)
	    sigU(i)=sigU(i) + wt*(sigOu(kl+1,kt,kc)-sigU(i))
	  endif
      
	endif	! kt.eq.ktslp .or. ...
      endif	! lc.gt.0 .and. ...

	! when both sigOu and sigOc are invalid, this data is invalid,
	! otherwise, it is valid.

      if(.not. (sigU(i) < 0. .and. sigC(i) < 0.)) then
	sigU(i)=max(0.,sigU(i))
	sigC(i)=max(0.,sigC(i))
      endif

    endif	! kc .gt. 0

  end do	! i=1,nobs; all observations

  deallocate(klev,wlev)
end
!.
