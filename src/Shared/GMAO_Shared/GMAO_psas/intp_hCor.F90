subroutine intp_hCor(	CorName,CorType, nlev_in,plev_in,	&
			mpar_in, npar_in,pars_in,		&
			nlev,   plev,   nctau,ctaus,		&
			norder, mlev, hC,hD,hE			)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: intp_hCor - compute horizontal correlations of all levels
!
! !INTERFACE: (to do)
!
! !DESCRIPTION:
!	intp_hCor() computes horizontal correlation matrix table of all
!	levels (plev) with a given parameter table (pars_in)
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
! 	18Dec95 - J. Guo	- programmed and added the prolog
!       16dec99 - da Silva     - Adde WIN_POWERLAW_3000km support
!_______________________________________________________________________
!#define _SHCOR	!
!#define _TRACE	!

#ifndef _REAL_
# define	_REAL_	real*8
#endif

use const, only : R_EARTH	! radius of Earth
use m_stdio, only : stderr	! standard error output unit
use m_die, only: die
use hcorfuns

implicit none
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  character(len=*), intent(in)	:: CorName	! Function name
  character(len=*), intent(in)	:: CorType	! Functional type
  integer, intent(in)	:: nlev_in	! size of plev_in
  real,    intent(in)	:: plev_in(nlev_in)	! levels of pars_in
  integer, intent(in)	:: mpar_in	! stride of pars_in
  integer, intent(in)	:: npar_in	! number of values in pars_in
  real,    intent(in)	:: pars_in(mpar_in,nlev_in)	! parameters

  integer, intent(in)	:: nlev		! size of plev
  real,    intent(in)	:: plev(nlev)	! values of levels
  integer, intent(in)	:: nctau	! size of ctaus
  _REAL_,  intent(in)	:: ctaus(nctau)	! values of ctaus
  integer, intent(in)	:: norder	! do derivatives?
  integer, intent(in)	:: mlev		! stride of hctbl
  real,    intent(out)	:: hC(mlev,nctau)	! function itself
  real,    intent(out)	:: hD(mlev,nctau)	! its first derivatives
  real,    intent(out)	:: hE(mlev,nctau)	! its second derivatives
!_______________________________________________________________________
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Local workspace and variables
  integer,allocatable :: klev(:)
  real,   allocatable :: wlev(:)
  _REAL_, allocatable :: htmp(:,:,:)
  integer,allocatable :: ierr_par(:)
  integer,allocatable :: klev_par(:)
  integer nlev_par

  logical valid			! if the model is valid
  logical compact		! if the model is compact-supported
  integer istat
  integer i,k,l
  integer kpar
  integer lk
  integer nk
  integer niter
  real    wk,s,r,dist,d3500,dmax
  _REAL_  dctau,dzero
  integer i_CorType

!   Parameters

  logical,	    parameter	:: roundoff=.true.
  character(len=*), parameter	:: myname="intp_hCor"
  real,parameter :: rade = .001*R_EARTH	! convert from m to km
!_______________________________________________________________________
#ifdef	_TRACE
	_TRACE	write(*,'(2a)') myname,': entered'
	_TRACE	write(*,'(5a)') myname,': ',CorName,' ',CorType
	_TRACE	write(*,*) myname,': ',	&
	_TRACE	  nlev_in,mpar_in,npar_in,nlev,nctau,mlev
#endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		! No data to processing, either by the input or the
		! output.

  if(nlev.le.0.or.nctau.le.0) return	! not an error

  if(nlev_in.le.0.or.npar_in.le.0.or.npar_in.gt.mpar_in) then
    if(nlev_in.le.0) write(stderr,'(2a)') myname,	&
	': yet to be defined, plev_in(nlev_in)'
    if(npar_in.le.0) write(stderr,'(2a)') myname,	&
	': yet to be defined, pars_in(npar_in,nlev_in)'
    if(npar_in.gt.mpar_in) write(stderr,'(a,2(a,i4))') myname,	&
	': invalid pars_in defination,  mpar_in = ',mpar_in,	&
	', npar_in = ',npar_in
    call die(myname)
  endif

  if(norder.ne.0.and.norder.ne.2) then
    write(stderr,'(2a,i3)') myname,			&
      & ': invalid parameter, norder =',norder
    call die(myname)
  endif

  if(norder.eq.2.and.nctau.lt.2) then
    write(stderr,'(2a,i3)') myname,			&
      & ': invalid parameter, nctau =',nctau
    call die(myname)
  endif
!-----------------------------------------------------------------------
		! Verify the function name and number of pars

  i_CorType=i_funcname(CorType)


  valid=.true.
  compact=.false.

		! both niter and dctau values are tuned, not by reasons
		! so far.  It appeared that niter=1 or 2 won't work.
		! They will cause a poor condition of the wind covar-
		! iance matrix.  The reason to try niter=1 or 2 is
		! because the spline function (Gaspari-Cohn) does not
		! have continues derivatives higher than second order.
  niter=3
  dzero=0.	! to ensure the arguments have proper types.
  dctau=0.
  if(norder.eq.2.and.nctau.ge.2) dctau=.1*(ctaus(2)-ctaus(1))

  select case (i_CorType)
  case ( i_DAMP_COSINE)
    kpar=6
  case (    i_GAUSSIAN)
    kpar=2
  case ( i_EXPONENTIAL)
    kpar=2
  case (   i_POWER_LAW)
    kpar=2
  case (i_GASPARI_COHN)
    kpar=2
    compact=.true.
  case (i_WIN_POWERLAW)
    kpar=2
    compact=.true.
  case (i_WIN_POWERLAW_3000km)
    kpar=2
    compact=.true.
  case default
    valid=.false.
    kpar=npar_in
  end select

  if(.not.valid) then

    l=len_trim(CorType)
    write(stderr,'(4a)') myname,': unknown function name, "',	&
	CorType(1:l),'"'
    call die(myname)

  elseif(npar_in < kpar) then

    l=len_trim(CorType)
    write(stderr,'(3a,2(a,i2))') myname,		&
	': invalid information for "',CorType(1:l),	&
	'", expecting ',kpar,					&
	' parameters, but given ',npar_in
    call die(myname)

  endif

#ifdef	_TRACE
	_TRACE	write(*,*) myname,': checked npar_in'
#endif
!-----------------------------------------------------------------------
	! index/weight plev to plev_in

	allocate(klev(nlev),wlev(nlev),			&
	  & ierr_par(1:nlev_in),klev_par(1:nlev_in),	&
	  & htmp(nctau,nlev_in,0:norder), stat=istat	)
	if(istat.ne.0) then
	  write(stderr,'(2a,i5)') myname,		&
		': allocate() error, stat =',istat
	  call die(myname)
	endif

	! Index the levels to be used for vertical interpolation
  call slogtab(.not.roundoff,nlev_in,plev_in,nlev,plev,klev,wlev)

	! Search for the levels having been indexed
  nlev_par=0
  do k=1,nlev
    lk=klev(k)
    nk=-1
    if(lk.lt.nlev_in) nk=lk+1
    
    	! Mark off existent levels
    do i=1,nlev_par
      if(lk.eq.klev_par(i)) lk=-1
      if(nk.eq.klev_par(i)) nk=-1
    end do

	! Add a level
    if(lk.ne.-1) then
      nlev_par=nlev_par+1
      klev_par(nlev_par)=lk
    endif

	! Add the level next to it
    if(nk.ne.-1) then
      nlev_par=nlev_par+1
      klev_par(nlev_par)=nk
    endif
  end do

#ifdef	_TRACE
	_TRACE	write(*,*) myname,': slogtab() done'
#endif
!_______________________________________________________________________
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  d3500=3500./rade
  dmax=acos(1.-ctaus(nctau))

!-----------------------------------------------------------------------

#ifdef _MT_
!cheung
#if(_UNICOS)
!MIC$	PARALLEL
!MIC$&		SHARED (nlev_par,klev_par,ierr_par,npar_in,pars_in,
!MIC$&			nctau,ctaus,dzero,dctau,htmp, i_CorType,
!MIC$&			norder, niter, d3500,dmax,compact,stderr)
!MIC$&		PRIVATE(k,lk,i,dist,wk,istat		)
!MIC$	DO PARALLEL
#endif
#if(sgi)
!$DOACROSS local(k,lk,i,dist,wk,istat            ),
!$&        mp_schedtype = dynamic,
!$&        share( nlev_par,klev_par,ierr_par,npar_in,pars_in,
!$&                 nctau,ctaus,dzero,dctau,htmp, i_CorType,
!$&                 norder, niter, d3500,dmax,compact,stderr)
#endif
#endif

  do k=1,nlev_par
    lk=klev_par(k)

		! Evaluate the function values

    call corfun(nctau,ctaus(1),dzero,htmp(1,lk,0),	&
      & i_CorType,npar_in,pars_in(1,lk), istat		)
    if(istat /= 0) then
      write(stderr,'(a,3(a,i3),a,i5)') myname,	&
	': corfun(k=',k, ', lk=',lk,		&
	', i_CorType=',i_CorType,		&
	') error, istat =',istat
    endif

    if(norder.eq.2.and.istat.eq.0) then

		! Evaluate the derivatives

      call richxtr(	nctau,ctaus(1),dctau,		&
	& htmp(1,lk,0),htmp(1,lk,1),htmp(1,lk,2),	&
	& i_CorType,npar_in,pars_in(1,lk),		&
	niter, istat	)
      if(istat /= 0) then
	write(stderr,'(a,3(a,i3),a,i5)') myname,	&
		': richxtr(k=',k, ', lk=',lk,		&
		', i_CorType=',i_CorType,		&
		') error, istat =',istat
      endif

    endif

    ierr_par(lk)=istat

	! apply a crude window truncation function.  A better approach
	! should be considered to replace this!

    if(istat.eq.0.and..not.compact) then
		! On the correlation functions
      do i=1,nctau
	dist=acos(1.-ctaus(i))
	if(dist.gt.d3500) then
	  wk=(dmax-dist)/(dmax-d3500)
	  htmp(i,lk,0)=htmp(i,lk,0)*wk
	endif
      end do

      if(norder.eq.2) then

		! On the derivatives of the correlation functions.  The
		! derivative of the window function is ignored.
	do i=1,nctau
	  dist=acos(1.-ctaus(i))
	  if(dist.gt.d3500) then
	    wk=(dmax-dist)/(dmax-d3500)
	    htmp(i,lk,1)=htmp(i,lk,1)*wk
	    htmp(i,lk,2)=htmp(i,lk,2)*wk
	  endif
	end do
      endif
    endif
  end do

#ifdef _MT_
#if(_UNICOS)
!MIC$	END DO
!MIC$	END PARALLEL
#endif
#endif

#ifdef	_TRACE
	_TRACE	write(*,*) myname,': richxtr() evaluation done'
	_TRACE	write(*,*) myname,': ierr =',	&
	_TRACE	  &	ierr_par(klev_par(1:nlev_par))
#endif

  istat=0
  do k=1,nlev_par
    lk=klev_par(k)
    if(ierr_par(lk).ne.0) then
      istat=lk
      exit
    endif
  end do

  if(istat.gt.0) then
    l=len_trim(CorType)
    write(stderr,'(3a,a,f8.3,a,i5)') myname,			&
      &	': error when evaluating htmp("',CorType(1:l),	&
      &	'",0:2) at level ',plev_in(istat),', status =',ierr_par(istat)
    call die(myname)
  endif
!-----------------------------------------------------------------------
	! Vertical linear interpolations

#ifdef _MT_
!cheung
#if(_UNICOS)
!MIC$	PARALLEL
!MIC$&		SHARED (nlev_in,klev,wlev,nctau,htmp,	&
!MIC$&			norder,nlev,hC,hD,hE		)
!MIC$&		PRIVATE(k,lk,wk,nk,i			)
!MIC$	DO PARALLEL
#endif
#if(sgi)
!$DOACROSS local(k,lk,wk,nk,i                    ),
!$&        mp_schedtype = dynamic,
!$&        share( nlev_in,klev,wlev,nctau,htmp,
!$&               norder,nlev,hC,hD,hE            )
#endif
#endif

  do k=1,nlev
    lk=klev(k)
    wk=max(0.,min(wlev(k),1.))	! only between [0.,1.)
    nk=lk+1

    do i=1,nctau
      hC(k,i) = htmp(i,lk,0)
    end do
    if(norder.eq.2) then
      do i=1,nctau
        hD(k,i) =-htmp(i,lk,1)
        hE(k,i) = htmp(i,lk,2)
      end do
    endif

    if(lk.lt.nlev_in) then
      do i=1,nctau
        hC(k,i) = hC(k,i) + wk*(htmp(i,nk,0)-hC(k,i))
      end do

      if(norder.eq.2) then
      	
	do i=1,nctau
          hD(k,i) = hD(k,i) + wk*(-htmp(i,nk,1)-hD(k,i))
          hE(k,i) = hE(k,i) + wk*( htmp(i,nk,2)-hE(k,i))
	end do
      endif
    endif

  end do
#ifdef _MT_
#if(_UNICOS)
!MIC$	END DO
!MIC$	END PARALLEL
#endif
#endif

#ifdef	_TRACE
	_TRACE	write(*,*) myname,': bilinear interpolation done'
#endif
!-----------------------------------------------------------------------
	!  List the interpolated table

#ifdef	_SHCOR
  _SHCOR  call slogtab(roundoff,nlev,plev,	&
  _SHCOR	& min(nlev,nlev_in),plev_in,klev,wlev)
  _SHCOR  do k=1,min(nlev,nlev_in)
  _SHCOR    lk=klev(k)
  _SHCOR    l=len_trim(CorName)
  _SHCOR    write(*,'(2a,i6.6)') '#HC{',CorName(1:l),nint(100.*plev(lk))
  _SHCOR
  _SHCOR    l=len_trim(CorType)
  _SHCOR    write(*,'(3a)')      '#HC ','# ',CorType(1:l)
  _SHCOR
  _SHCOR    if(norder.eq.2) then
  _SHCOR      write(*,'(2a,6x,a)') '#HC ','# ',	&
  _SHCOR        & '1-tau    r/a    s/a      s     hCor   dC/dT  dC2/d2T'
  _SHCOR      do i=1,nctau,1
  _SHCOR        s=acos(1.-ctaus(i))
  _SHCOR        r=2.*sin(.5*s)
  _SHCOR        write(*,'(a,i4,x,f7.6,2f7.4,f7.1,1p,3e12.4)') '#HC ', &
  _SHCOR	  & i,ctaus(i),r,s,s*rade,hC(lk,i),hD(lk,i),hE(lk,i)
  _SHCOR      end do
  _SHCOR    else
  _SHCOR      write(*,'(2a,6x,a)') '#HC ','# ',	&
  _SHCOR        & '1-tau    r/a    s/a      s     hCor'
  _SHCOR      do i=1,nctau,1
  _SHCOR        s=acos(1.-ctaus(i))
  _SHCOR        r=2.*sin(.5*s)
  _SHCOR        write(*,'(a,i4,x,f7.6,2f7.4,f7.1,1p,e12.4)') '#HC ', &
  _SHCOR	  & i,ctaus(i),r,s,s*rade,hC(lk,i)
  _SHCOR      end do
  _SHCOR    endif
  _SHCOR    write(*,'(a)') '#HC}::'
  _SHCOR  end do
  _SHCOR  write(*,'(a)') '::'
#endif
!-----------------------------------------------------------------------

	deallocate(klev,wlev,klev_par,ierr_par,htmp)
end subroutine intp_hCor
!.







