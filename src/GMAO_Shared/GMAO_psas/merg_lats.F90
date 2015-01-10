!#define	_TRACE	!
!#define	_DEBUG	!
subroutine merg_lats(ngrdlat,nobslat,obslats, MXveclat,nveclat,veclats)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: merg_lats - (to do)
!
! !INTERFACE: (to do)
! !DESCRIPTION: (to do)
! !EXAMPLES: (to do)
! !BUGS: (to do)
! !SEE ALSO: (to do)
! !SYSTEM ROUTINES: (to do)
!
! !REVISION HISTORY:
! 	04Mar96 - J. Guo	- (to do)
!	14Aug96 - J. Guo	- Fixed for cases of ngrdlat == 0.
!_______________________________________________________________________
use m_stdio,only : stderr
use m_die,  only : die
implicit none

!@interface

  integer,intent(in)	:: ngrdlat	! number of grid latitudes

  integer,intent(in)	:: nobslat	! number of obs. latitudes
     real,intent(in)	:: obslats(nobslat)	! obs. latitudes

  integer,intent(in)	:: MXveclat	! maximum size of veclats(:)
  integer,intent(inout)	:: nveclat	! actual size of veclats(:)
     real,intent(inout)	:: veclats(MXveclat)	! tabulared lat. values

!@end/interface
!=======================================================================
	! Local vars.
  integer j			! loop counter
  real    dlat			! grid interval
  real	  res			! expected resolution
  real, allocatable :: grdlats(:) ! a workspace, for no side-effect
  integer istat

  character(len=*),parameter :: myname='merg_lats'

!=======================================================================
	! consistency checking

  if(ngrdlat.gt.MXveclat) then
    write(stderr,'(2a,i3)') myname,	&
      ': not enough workspace, MXveclat, <',ngrdlat
    call die(myname)
  endif

#ifdef	_TRACE
	_TRACE	write(*,'(2a,4i6)') myname,': entered, ',	&
	_TRACE	  ngrdlat,nobslat,Mxveclat,nveclat
#endif

	! merging grided latitude values, if there is any

  if (ngrdlat > 0) then

	! grid point interval

    dlat=180./(ngrdlat-1)

	! assign grid point values to the table, grdlats(:)

	allocate(grdlats(ngrdlat), stat=istat)
	if(istat.ne.0) call die(myname,'allocate()',istat)

    do j=1,ngrdlat
      grdlats(j)=(j-1)*dlat-90.
    end do
  endif

	! an assigned resolution of the final resultant list

  res=180./(MXveclat-1)

	! a merge to make sure the table agreed with the assumed order
	! by tabRlist(), both input and the output used the same memory
	! location.

#ifdef	_TRACE
	_TRACE	write(*,'(2a,2i4)') myname,	&
	_TRACE	  ': first tabRlist() call, ',nveclat,ngrdlat
#endif

  if (ngrdlat > 0) then
    call tabRlist(res, ngrdlat,grdlats, MXveclat,nveclat,veclats)

#ifdef	_DEBUG
	_DEBUG	write(*,'(2a)') myname,':: Grid latitudes'
	_DEBUG	do j=1,nveclat,10
	_DEBUG	  write(*,'(i3,x,10f7.2)') j,veclats(j:min(j+9,nveclat))
	_DEBUG	end do
	_DEBUG	write(*,'(a)') '::'
#endif
  endif

	! merge in the observations

#ifdef	_TRACE
	_TRACE	write(*,*) 'second tabRlist() call, ',nveclat,MXveclat
#endif

  call tabRlist(res, nobslat,obslats, MXveclat,nveclat,veclats)

#ifdef	_DEBUG
	_DEBUG	write(*,'(2a)') myname,':: Sieved latitudes'
	_DEBUG	do j=1,nveclat,10
	_DEBUG	  write(*,'(i3,x,10f7.2)') j,veclats(j:min(j+9,nveclat))
	_DEBUG	end do
	_DEBUG	write(*,'(a)') '::'
#endif

  if (ngrdlat > 0) deallocate(grdlats)
end subroutine merg_lats
