!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_lookup - look up for indices
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_lookup
      implicit none
      private	! except

      public :: lookup		! The class data structure

      interface lookup; module procedure	&
	lookup_,	&
	lookupweights_,	&
        lookupgrid_,	&
        lookupgridweights_
      end interface

! !REVISION HISTORY:
! 	08Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_lookup'

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: lookupweights_ - lookup values for indices and weights
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine lookupweights_(rlev, klev,wlev, rtbl,linear)
      implicit none
      real   ,dimension(:),intent(in ) :: rlev	! values to lookup
      integer,dimension(:),intent(out) :: klev	! indices
      real   ,dimension(:),intent(out) :: wlev	! weights

      real,dimension(:),intent(in) :: rtbl	! the lookup table
      logical          ,intent(in) :: linear	! if not log-linear

! !REVISION HISTORY:
! 	08Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::lookupweights_'

  logical,parameter :: NEAREST=.true.

  select case(linear)
  case(.true.)
    call slintab(.not.NEAREST, size(rtbl),rtbl,		&
	size(rlev),rlev, klev,wlev)

  case default
    call slontab(.not.NEAREST, size(rtbl),rtbl,		&
	size(rlev),rlev, klev,wlev)
  end select

end subroutine lookupweights_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: lookup_ - lookup values for indices only
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine lookup_(rlev, klev, rtbl,linear)
      use m_die,only : die
      implicit none
      real   ,dimension(:),intent(in ) :: rlev	! values to lookup
      integer,dimension(:),intent(out) :: klev	! indices

      real,dimension(:),intent(in) :: rtbl	! the lookup table
      logical          ,intent(in) :: linear	! if not log-linear

! !REVISION HISTORY:
! 	08Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::lookup_'

  logical,parameter :: NEAREST=.true.
  real,allocatable,dimension(:) :: wlev
  integer :: nlev
  integer :: ier

  nlev=size(rlev)

	allocate(wlev(nlev),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)

  select case(linear)
  case(.true.)
    call slintab(.not.NEAREST, size(rtbl),rtbl,		&
	size(rlev),rlev, klev,wlev)

  case default
    call slontab(.not.NEAREST, size(rtbl),rtbl,		&
	size(rlev),rlev, klev,wlev)
  end select

	deallocate(wlev,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

end subroutine lookup_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: lookupgridweights_ - lookup values on a grid mesh
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine lookupgridweights_(rdat, idat,wdat,	&
	grid_min,grid_range,ngrid,periodic)

      use m_die,only : die
      implicit none

      real   ,dimension(:),intent(in ) :: rdat
      integer,dimension(:),intent(out) :: idat
      real   ,dimension(:),intent(out) :: wdat

      real   ,intent(in) :: grid_min
      real   ,intent(in) :: grid_range
      integer,intent(in) :: ngrid
      logical,optional,intent(in) :: periodic

! !REVISION HISTORY:
! 	08Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::lookupgridweights_'
  logical :: periodic_
  logical :: valid
  integer :: i
  real :: del,xdat

  periodic_=.false.
  if(present(periodic)) periodic_=periodic

  valid=ngrid>0
  if(valid .and. periodic_) valid=ngrid>1
  if(.not.valid) call die(myname_,'invalid ngrid',ngrid)

  del=grid_range/ngrid
  if(.not.periodic_) del=grid_range/(ngrid-1)

  do i=1,size(rdat)
    xdat=(rdat(i)-grid_min)/del

    idat(i)=floor(xdat)
    if(periodic_) then
		! idat isin [0,ngrid-1)
      idat(i)=modulo(idat(i),ngrid)
    else
      idat(i)=max(0,min(ngrid-2,idat(i)))
    endif

    wdat(i)=xdat-idat(i)
    idat(i)=idat(i)+1
  end do

end subroutine lookupgridweights_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: lookupgrid_ - lookup values on a grid mesh
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine lookupgrid_(rdat, idat,	&
	grid_min,grid_range,ngrid,periodic)

      use m_die,only : die
      implicit none

      real   ,dimension(:),intent(in ) :: rdat
      integer,dimension(:),intent(out) :: idat

      real   ,intent(in) :: grid_min
      real   ,intent(in) :: grid_range
      integer,intent(in) :: ngrid
      logical,optional,intent(in) :: periodic

! !REVISION HISTORY:
! 	08Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::lookupgrid_'
  logical :: periodic_
  logical :: valid
  integer :: i
  real :: del,xdat

  periodic_=.false.
  if(present(periodic)) periodic_=periodic

  valid=ngrid>0
  if(valid .and. periodic_) valid=ngrid>1
  if(.not.valid) call die(myname_,'invalid ngrid',ngrid)

  del=grid_range/ngrid
  if(.not.periodic_) del=grid_range/(ngrid-1)

  do i=1,size(rdat)
    xdat=(rdat(i)-grid_min)/del
    idat(i)=nint(xdat)

    if(periodic_) then
      idat(i)=modulo(idat(i),ngrid)
    else
      idat(i)=max(0,min(ngrid-2,idat(i)))
    endif

    idat(i)=idat(i)+1
  end do

end subroutine lookupgrid_
end module m_lookup
