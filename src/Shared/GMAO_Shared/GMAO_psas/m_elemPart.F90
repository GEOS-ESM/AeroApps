!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_elemPart - Elemntal parallel paritioner
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_elemPart
      implicit none
      private	! except

      public :: elemPart		! The class data structure

      interface elemPart; module procedure	&
	indexedPartI_,	&
	indexedPartR_,	&
	segmentPartI_,	&
	segmentPartR_
      end interface

! !REVISION HISTORY:
! 	18Mar99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_elemPart'

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexedPartR_ - an indexed parallel Cost partitioner
!
! !DESCRIPTION:
!
!   There is an assumping that all PEs have the same computation
!   efficiency for this routine.  Costs partition algorithms with
!   diffierent PE computation efficiencies may be implemented with
!   minor changes.
!
! !INTERFACE:

    subroutine indexedPartR_(Costs,iCostx,nCostx,comm,	&
	TotalCost,myCost,noSort)

      use m_SortingTools, only : indexSort
      use m_recurPart,    only : recurPart
      use m_die,      only : MP_perr_die
      use m_mpif90,only : MP_comm_rank
      use m_mpif90,only : MP_comm_size

      implicit none

      real,   dimension(:),intent(in)    :: Costs  ! a list of costs
      integer,dimension(:),intent(inout) :: iCostx ! indices to Costs
      integer,             intent(inout) :: nCostx ! <=size(iCostx)

      integer,intent(in) :: comm

      real,   optional,intent(out) :: TotalCost
      real,   optional,intent(out) :: myCost

			! The default action of indexedPart() is to
			! sort the accessing (card-dealing) order by
			! their _Costs_ values first.  The sorting
			! will be _stable_, such that the order of
			! any elements with the same _Costs_ value
			! will not be altered.  With noSort=.true.,
			! the default costs-sorting action will be
			! turned off.  The accessing order will be
			! as is.

      logical,optional,intent(in)  :: noSort

! !REVISION HISTORY:
! 	17Mar99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::indexedPartR_'

  integer :: i,l,ier
  integer :: nPE,lPE,uPE,myPE
  real    :: TotalCost_,myCost_
  logical :: noSort_

  if(nCostx <= 0) return	! No work left to _divide_

  call MP_comm_size(comm,nPE,ier)
	if(ier /= 0) call MP_perr_die(myname_,'MPI_comm_size',ier)

	! This is the initial range

  lPE=0
  uPE=nPE-1

  if(uPE <= lPE) return	! There is only one PE in this communicator

  call MP_comm_rank(comm,myPE,ier)
	if(ier /= 0) call MP_perr_die(myname_,'MPI_comm_rank',ier)

	! Set the default partitioning algorithm

  noSort_=.false.
  if(present(noSort)) noSort_=noSort

	! Sort the order according to their costs.

  if(.not.noSort_)	&
    call IndexSort(nCostx,iCostx,Costs,descend=.true.)

	! Estimate the _total_ cost to be divided

  TotalCost_=0.
  do i=1,nCostx
    l=iCostx(i)
    TotalCost_=TotalCost_+Costs(l)
  end do
  if(present(TotalCost)) TotalCost=TotalCost_

  myCost_=TotalCost_

	! Distribute the costs through an indexed recursive bi-section
	! prorcess

  call recurPart(Costs,iCostx,nCostx,myCost_, myPE,lPE,uPE)

  if(present(myCost)) myCost=myCost_

end subroutine indexedPartR_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexedPartI_ - an indexed parallel Cost partitioner
!
! !DESCRIPTION:
!
!   There is an assumping that all PEs have the same computation
!   efficiency for this routine.  Costs partition algorithms with
!   diffierent PE computation efficiencies may be implemented with
!   minor changes.
!
! !INTERFACE:

    subroutine indexedPartI_(Costs,iCostx,nCostx,comm,	&
	TotalCost,myCost,noSort)

      use m_SortingTools,only : indexSort
      use m_recurPart,   only : recurPart
      use m_die,      only : MP_perr_die
      use m_mpif90,   only : MP_comm_size
      use m_mpif90,   only : MP_comm_rank

      implicit none

      integer,dimension(:),intent(in)    :: Costs  ! a list of costs
      integer,dimension(:),intent(inout) :: iCostx ! indices to Costs
      integer,             intent(inout) :: nCostx ! <=size(iCostx)

      integer,intent(in) :: comm

      integer,optional,intent(out) :: TotalCost
      integer,optional,intent(out) :: myCost

			! The default action of indexedPart() is to
			! sort the accessing (card-dealing) order by
			! their _Costs_ values first.  The sorting
			! will be _stable_, such that the order of
			! any elements with the same _Costs_ value
			! will not be altered.  With noSort=.true.,
			! the default costs-sorting action will be
			! turned off.  The accessing order will be
			! as is.

      logical,optional,intent(in)  :: noSort

! !REVISION HISTORY:
! 	17Mar99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::indexedPartI_'

  integer :: i,l,ier
  integer :: nPE,lPE,uPE,myPE
  integer :: TotalCost_,myCost_
  logical :: noSort_

  if(nCostx <= 0) return	! No work left to _divide_

  call MP_comm_size(comm,nPE,ier)
	if(ier /= 0) call MP_perr_die(myname_,'MPI_comm_size',ier)

	! This is the initial range

  lPE=0
  uPE=nPE-1

  if(uPE <= lPE) return	! There is only one PE in this communicator

  call MP_comm_rank(comm,myPE,ier)
	if(ier /= 0) call MP_perr_die(myname_,'MPI_comm_rank',ier)

	! Set the default partitioning algorithm

  noSort_=.false.
  if(present(noSort)) noSort_=noSort

	! Sort the order according to their costs.

  if(.not.noSort_)	&
    call IndexSort(nCostx,iCostx,Costs,descend=.true.)

	! Estimate the _total_ cost to be divided

  TotalCost_=0
  do i=1,nCostx
    l=iCostx(i)
    TotalCost_=TotalCost_+Costs(l)
  end do
  if(present(TotalCost)) TotalCost=TotalCost_

  myCost_=TotalCost_

	! Distribute the costs through an indexed recursive bi-section
	! prorcess

  call recurPart(Costs,iCostx,nCostx,myCost_, myPE,lPE,uPE)

  if(present(myCost)) myCost=myCost_

end subroutine indexedPartI_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: segmentPartI_ - Segmented parallel partition
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine segmentPartI_(ncost,costs, lbound,ubound, comm,	&
	TotalCost,myCost,nozero					)

      use m_mpif90,   only : MP_comm_size,MP_comm_rank
      use m_die,      only : MP_perr_die
      use m_die,      only : perr_die
      use m_mall,     only : mall_ison,mall_mci,mall_mco
      use m_recurPart,only : recurPart

      implicit none

      integer,             intent(in)  :: ncost
      integer,dimension(:),intent(in)  :: costs ! costs
      integer,             intent(out) :: lbound,ubound
      integer,             intent(in)    :: comm

      integer,optional,intent(out) :: TotalCost
      integer,optional,intent(out) :: myCost

      logical,optional,intent(in) :: nozero

! !REVISION HISTORY:
! 	17Mar99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::segmentPartI_'

  integer :: nPE,lPE,uPE,myPE
  integer :: i,ier
  integer,allocatable,dimension(:) :: acost
  integer :: myCost_,TotalCost_

!________________________________________

		! allocate the workspace of accumulated costs array

	allocate(acost(0:ncost),stat=ier)
	if(ier /= 0) call perr_die(myname_,'allocate()',ier)
	if(mall_ison()) call mall_mci(acost,myname)

  acost(0)=0
  do i=1,ncost
    acost(i)=acost(i-1)+costs(i)
  end do
  TotalCost_=acost(ncost)

  lbound=1
  ubound=ncost

  myCost_=TotalCost_

  if(present(TotalCost)) TotalCost=TotalCost_
  if(present(   myCost))    myCost=   myCost_

!________________________________________

  call MP_comm_size(comm,nPE,ier)
	if(ier /= 0) call MP_perr_die(myname_,'MPI_comm_size',ier)

	! This is the initial range

  lPE=0
  uPE=nPE-1

  if(uPE <= lPE) then
	if(mall_ison()) call mall_mco(acost,myname)
    deallocate(acost,stat=ier)
	if(ier /=0) call perr_die(myname_,'deallocate()',ier)

    return	! There is only one PE in this communicator
  endif

  call MP_comm_rank(comm,myPE,ier)
	if(ier /= 0) call MP_perr_die(myname_,'MPI_comm_rank',ier)
	! Bisection with a workspace defined
!________________________________________

  call recurPart(acost, lbound,ubound, myCost_, myPE,lPE,uPE,	&
	nozero=nozero)
  if(present(myCost)) myCost=myCost_

		! deallocate the workspace of accumulated costs array

		if(mall_ison()) call mall_mco(acost,myname)
	deallocate(acost,stat=ier)
	if(ier /= 0) call perr_die(myname_,'deallocate()',ier)

end subroutine segmentPartI_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: segmentPartR_ - Segmented parallel partition
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine segmentPartR_(ncost,costs, lbound,ubound, comm,	&
	TotalCost,myCost,nozero					)

      use m_mpif90,   only : MP_comm_size,MP_comm_rank
      use m_die,      only : MP_perr_die
      use m_die,      only : perr_die
      use m_mall,     only : mall_ison,mall_mci,mall_mco
      use m_recurPart,only : recurPart

      implicit none

      integer,             intent(in)  :: ncost
      real   ,dimension(:),intent(in)  :: costs ! costs
      integer,             intent(out) :: lbound,ubound
      integer,             intent(in)    :: comm

      real   ,optional,intent(out) :: TotalCost
      real   ,optional,intent(out) :: myCost

      logical,optional,intent(in)  :: nozero

! !REVISION HISTORY:
! 	17Mar99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::segmentPartR_'

  integer :: nPE,lPE,uPE,myPE
  integer :: i,ier
  real,allocatable,dimension(:) :: acost
  real    :: myCost_,TotalCost_

!________________________________________

		! allocate the workspace of accumulated costs array

	allocate(acost(0:ncost),stat=ier)
	if(ier /= 0) call perr_die(myname_,'allocate()',ier)
	if(mall_ison()) call mall_mci(acost,myname)

  acost(0)=0.
  do i=1,ncost
    acost(i)=acost(i-1)+costs(i)
  end do
  TotalCost_=acost(ncost)

  lbound=1
  ubound=ncost

  myCost_=TotalCost_

  if(present(TotalCost)) TotalCost=TotalCost_
  if(present(   myCost))    myCost=   myCost_

!________________________________________

  call MP_comm_size(comm,nPE,ier)
	if(ier /= 0) call MP_perr_die(myname_,'MPI_comm_size',ier)

	! This is the initial range

  lPE=0
  uPE=nPE-1

  if(uPE <= lPE) then
	if(mall_ison()) call mall_mco(acost,myname)
    deallocate(acost,stat=ier)
	if(ier /=0) call perr_die(myname_,'deallocate()',ier)

    return	! There is only one PE in this communicator
  endif

  call MP_comm_rank(comm,myPE,ier)
	if(ier /= 0) call MP_perr_die(myname_,'MPI_comm_rank',ier)
	! Bisection with a workspace defined
!________________________________________

  call recurPart(acost, lbound,ubound, myCost_, myPE,lPE,uPE,	&
	nozero=nozero)
  if(present(myCost)) myCost=myCost_

		! deallocate the workspace of accumulated costs array

		if(mall_ison()) call mall_mco(acost,myname)
	deallocate(acost,stat=ier)
	if(ier /= 0) call perr_die(myname_,'deallocate()',ier)

end subroutine segmentPartR_

end module m_elemPart
!.
