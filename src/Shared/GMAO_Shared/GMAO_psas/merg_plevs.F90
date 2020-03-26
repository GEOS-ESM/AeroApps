!@interface
subroutine merg_plevs(nolev,olevs,nalev,alevs,mxlev,nplev,plevs)
!@end/interface

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: merg_plevs - merges two given lists of pressure values
!
! !INTERFACE:
!	<@interface
!
! !INPUT PARAMETERS:
!     integer nolev	! array size of olevs
!     real    olevs(:)	! one of two input pressure level list
!     integer nalev	! array size of alevs
!     real    alevs(:)	! one of two input pressure level list
!     integer mxlev	! array size limit of plevs
!
! !OUTPUT PARAMETERS:
!     integer nplev	! actual array size of plevs
!     real    plevs(:)	! a list of all levels in olevs/alevs
!
! !DESCRIPTION:
!	Merg_plevs merges unique pressure levels given by the two input
!	lists, olevs and alevs.  The returned list in plevs is limited
!	to a given range with a given precision (see the source code for
!	information about (pmin, pmax, and pres).  For this version, the
!	values of (pmin, pmax, pres) is (.01, 1000., .01).
!
! !EXAMPLES:
!	For levels in an innovation vector (rlevs(nobs)) and an analysis
!	grid (glevs(mlev)),
!
!		call merg_plevs(nobs,rlevs,mlev,glevs,MXlev,nplev,plevs)
!
!	For more than two input arrays,
!
!		call merg_plevs(na,pa,     0,plevs, MXlev,nplev,plevs)
!		call merg_plevs(nb,pb, nplev,plevs, MXlev,nplev,plevs)
!		call merg_plevs(nc,pc, nplev,plevs, MXlev,nplev,plevs)
!		...
!
! !BUGS:
!
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
! !REVISION HISTORY:
! 	28Nov95 - J. Guo	- added the prolog
!_______________________________________________________________________
!#define	_TRACE	!
!#define	_DEBUG	!

!!$  use m_mpout, only : mpout_log
  use m_mpout, only : mpout
  use m_die,  only : die
  use m_mall, only : mall_ison,mall_mci,mall_mco

  implicit none

!-----------------------------------------------------------------------
!   Arguments:
!  ============

!@interface
  integer,intent(in)	:: nolev	! array size of olevs
  real,   intent(in)	:: olevs(nolev)	! a pressure level list
  integer,intent(in)	:: nalev	! array size of alevs
  real,   intent(in)	:: alevs(nalev)	! a pressure level list
  integer,intent(in)	:: mxlev	! size limit of plevs

  integer,intent(out)	:: nplev	! actual size of plevs
  real,   intent(out)	:: plevs(mxlev)	! a list of all levels
!@end/interface

!=======================================================================
!   Local vars.
!  =============
  integer		:: nwlev	! a size of the workspace
  integer, allocatable	:: ilevs(:)	! integer level values
  integer, allocatable	:: lndex(:)	! a sorted index
  real,    allocatable	:: all_levs(:)  ! list of all level information
  integer		:: istat

!-----------------------------------------------------------------------
!   Parameters
!  ============

  real	pmin,pmax,pres
  parameter(pmin= .01 )		! minimum pressure level range
  parameter(pmax=1200.)		! maximum pressure level range
  parameter(pres= .01 )		! a relative resolution

  character(len=*), parameter	:: myname='merg_plevs'

!-----------------------------------------------------------------------
!   Local vars.
!  ============
  integer	i,l
  real		p,rstp
  real		pflr,pfix
  real		smin,smax,sfix
  real          plast, pnext
  integer       ilev, ilev_last, ll

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		! Allocate workspace

	! set resolution for the first pass.  For pres=.01, it is a
	! resolution step (sfix) of 10 mb for p > 1000 mb.  It is
	! however limited to a value best for lower atmosphere of
	! pressure p > 100 mb.

#ifdef	_TRACE
	_TRACE	write(stdout,'(2a,$)') myname,': entered'
	_TRACE	write(stdout,*) nolev,nalev,mxlev
#endif

!-----------------------------------------------------------------------
!   Count the number of valid levels
!  ==================================

!!$	Call mpout_log(myname,'nolev =',nolev)
!!$	Call mpout_log(myname,'nalev =',nalev)
!!$        call mpout_log(myname,'mxlev =',mxlev)

	
	allocate(all_levs(nolev+nalev),ilevs(nolev+nalev), lndex(nolev+nalev), stat=istat)
	if(istat.ne.0) call die(myname,'allocate()',istat)
	if(mall_ison()) then
	  call mall_mci(all_levs,myname)
	  call mall_mci(ilevs,myname)
	  call mall_mci(lndex,myname)
	endif


	all_levs(1:nolev)             = olevs(1:nolev)
	all_levs(nolev+1:nolev+nalev) = alevs(1:nalev)

!-----------------------------------------------------------------------
!   Sort all levels
!  =========================
#ifdef	_TRACE
	_TRACE	write(stdout,'(2a)') myname,': calling index-permute'
#endif

	  call indexxr(nolev+nalev,all_levs,lndex)      	! build a heapsort index
	  call permutr(all_levs,lndex,nolev+nalev,all_levs)	! sort the array

!!$          call mpout_log(myname,'first = ',all_levs(1))
!!$          call mpout_log(myname,'last  = ',all_levs(nolev+nalev))

rstp=pres

do	! until(nwlev.le.mxlev.or.rstp.le.0.)
		! set parameters for pressure truncation/scaling

#ifdef	_TRACE
	_TRACE	write(stdout,'(2a,f8.3)') myname,	&
	_TRACE	  ': loop with rstp = ',rstp
#endif

  p=999.
  pflr=10.**floor(log10(p))		! truncate p=999 to a log10 step
  smax=rstp*pflr			! the precision for the step

  p=pmin
  pflr=10.**floor(log10(p))		! truncate pmin to a log10 step
  smin=rstp*pflr			! the precision for the step

  l=0		! set the counter
  ll = 0
  pnext = pmin
  ilev_last = -1  ! initial impossible value

  do i=1,nolev + nalev

     p = all_levs(i)
     If (p >= pnext .and. p <= pmax) Then
        ll = ll + 1


                ! truncate the pressure to 1000 100 10 1 .1 .01 mb
        pflr=10.**floor(log10(p))

		! the precision for this log10 step, but at least 1 mb
        sfix=min(smax,rstp*pflr)

		! round off the pressure to nearest sfix step
        pfix=sfix*nint(p/sfix)
        pnext=sfix*(nint(p/sfix)+0.4)

		! scale the pressure to the minimum log10 step
        ilev = nint(pfix/smin)

        If (ilev /= ilev_last) Then
           l = l + 1
           ilevs(l)= ilev
           ilev_last = ilev
           If (l > mxlev) exit
        End if

     endif

  end do
!!$  Call mpout_log(myname,' ll   = ',ll)

  nwlev=l	! set the number of levels
!!$  call mpout_log(myname,'nwlev = ',nwlev)

  if(nwlev.gt.mxlev) then
    write(mpout,'(2a,i5,$)') myname,': nwlev = ',nwlev
    write(mpout,'(a,f4.2,$)') ', readjusted rstp from ',rstp
    if (rstp .le. .01) then
      rstp=.02
    elseif (rstp .le. .02) then
      rstp=.05
    elseif (rstp .le. .05) then
      rstp=.1
    elseif (rstp .le. .1) then
      rstp=.2
    elseif (rstp .le. .2) then
      rstp=.5
    elseif (rstp .le. .5) then
      rstp=-1.
    endif
    write(mpout,'(a,f4.2)') ' to ',rstp
  endif

  if(nwlev.le.mxlev.or.rstp.le.0.) exit
enddo		! do until(nwlev.le.mxlev.or.rstp.le.0)

!-----------------------------------------------------------------------
!   Convert back to normal pressure unit
!  ======================================
  nplev=min(nwlev,mxlev)
  plevs(1:nplev)=smin*ilevs(1:nplev)	! scaled with smin

!-----------------------------------------------------------------------
!   For debugging only
!  ====================
#ifdef	_DEBUG
	_DEBUG	write(stdout,'(2a)') myname,'::'
	_DEBUG	write(stdout,'(a)') '# Sieved vertical levels'
	_DEBUG	write(stdout,'(t5,a,i4)') 'nplev ',nplev
	_DEBUG	do l=1,nplev,10
	_DEBUG	  write(stdout,'(f7.2,9f8.2)')	&
	_DEBUG	    (plevs(i),i=l,min(nplev,l+9))
	_DEBUG	end do
	_DEBUG	write(stdout,'(a)') '::'
#endif

		if(mall_ison()) then
		  call mall_mco(lndex,myname)
		  call mall_mco(all_levs,myname)
		  call mall_mco(ilevs,myname)
		endif
	deallocate(lndex,all_levs,ilevs)
!_______________________________________________________________________
end subroutine merg_plevs
!.
