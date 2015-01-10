subroutine intp_vCor( CorName,mxlev_in,nlev_in,plev_in,vCor_in,	&
		      mxlev,   nlev,   plev,   vCor		)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: intp_vCor - interpolate a vCor matrix to finer levels
!
! !INTERFACE: (to do)
!
! !DESCRIPTION:
!	intp_vCor() interpolate a vCor matrix to a matrix with finer
!	levels.  The input matrix is expected to be already modified
!	arround its diagonal line from symmetric, to assymetric about
!	1. (i.e. vCor_in is assymmetric and diag(vCor_in)=1).  The
!	returned matrix vCor, however, is symmetric with diag(vCor)=1.
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
!#define _DEBUG	!

use m_die,only : die

implicit none

	! specify the input data structure

character(len=*), intent(in)	:: CorName	! Function name
integer, intent(in)	:: mxlev_in
integer, intent(in)	:: nlev_in
real,    intent(in)	:: plev_in(nlev_in)
real,    intent(in)	:: vCor_in(mxlev_in,nlev_in)

	! specify the output data structure
integer, intent(in)	:: mxlev
integer, intent(in)	:: nlev
real,    intent(in)	:: plev(nlev)
real,    intent(out)	:: vCor(mxlev,nlev)

	! local vars.

  integer  k, l
  integer lk,ll
  real    wk,wl
  integer nk,nl
  real    vc,an,al
  integer istat

  integer, allocatable :: klev(:)
  real,    allocatable :: wlev(:)

	! parameters

  logical,	    parameter	:: nearest=.true.
  character(len=*), parameter	:: myname="intp_vCor"

!=======================================================================
	! index plev to plev_in table

	allocate(klev(nlev),wlev(nlev), stat=istat)
	if(istat.ne.0)	&
	  call die(myname,'allocate()',istat)

  call slogtab(.not.nearest,nlev_in,plev_in,nlev,plev,klev,wlev)

	! interpolation for each level

  do k=1,nlev
    lk=klev(k)
    wk=wlev(k)
    nk=lk+1

	! The diagonal elements should be 1, although it looks like a
	! linear interpolation.

    vCor(k,k)=1.

	! Off diagonal elements are bilinearly interpolated with
	! symmetric implied.

    do l=1,k-1
      ll=klev(l)
      wl=wlev(l)
      nl=ll+1

      al=vCor_in(ll,lk)+wl*(vCor_in(nl,lk)-vCor_in(ll,lk))
      an=vCor_in(ll,nk)+wl*(vCor_in(nl,nk)-vCor_in(ll,nk))

      vc=al+wk*(an-al)

      vCor(l,k)=min(vc,2.-vc)
      vCor(k,l)=vCor(l,k)
    end do
  end do

!=======================================================================
	! List the interpolated table at the input levels

#ifdef	_DEBUG
  _DEBUG call slogtab(nearest,nlev,plev,		&
  _DEBUG   & min(nlev,nlev_in),plev_in,klev,wlev	)
  _DEBUG l=max(len_trim(CorName),1)
  _DEBUG write(*,'(/4a)') myname,'*',CorName(1:l),'::'
  _DEBUG do k=1,min(nlev,nlev_in)
  _DEBUG   if(plev(klev(k)).ge.1.) then
  _DEBUG     write(*,'(i4,x,$)') nint(plev(klev(k)))
  _DEBUG   else
  _DEBUG     write(*,'(x,f3.2,x,$)') plev(klev(k))
  _DEBUG   endif
  _DEBUG   do l=1,k
  _DEBUG     write(*,'(x,f4.2,$)') vCor(klev(l),klev(k))
  _DEBUG   end do
  _DEBUG   write(*,*)
  _DEBUG end do
  _DEBUG write(*,'(a)') '::'
#endif

	deallocate(klev,wlev)
end subroutine intp_vCor
!.
