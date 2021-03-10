!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_rectPartMatx - an order-preserving cost-weighted parallel partitioner
!        
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_rectPartMatx
      implicit none
      private	! except

      public :: PartMatx		! The class data structure

      interface PartMatx; module procedure part_; end interface

! !REVISION HISTORY:
! 	18Mar99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_rectPartMatx'

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: part_ - a cost-weighted parallel partitioner
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine part_(Matx,iMatx,nMatx,myCost,comm)

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : ptr_rAttr
      use m_AttrVect, only : AttrVect_indexRA => indexRA
      use m_SortingTools,only : IndexSort
      use m_recurPart, only : recurPart
      use m_mpif90,   only : MP_comm_size,MP_comm_rank
      use m_stdio,    only : stderr
      use m_die,      only : MP_perr_die
      use m_die,      only : die

      use m_mpout, only : mpout_log
      use m_mall,  only : mall_ison, mall_mci, mall_mco
      implicit none

      type(AttrVect),intent(in) :: Matx	! a list of matrix blocks

      integer,dimension(:),intent(inout) :: iMatx ! indices to Matx
      integer,             intent(inout) :: nMatx ! <=size(iMatx)
      real,                intent(inout) :: myCost

      integer,intent(in) :: comm

! !REVISION HISTORY:
! 	17Mar99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::part_'

  real,pointer,dimension(:,:) :: rAttr
  real,allocatable :: acost(:)

  integer :: nPE,lPE,uPE,myPE
  integer :: ier,i,l
  integer :: icost
  integer :: lb, ub
  Integer :: n_threads, i_thread

  if(nMatx <= 0) return	! No work left to _divide_

  call MP_comm_size(comm,nPE,ier)
	if(ier /= 0) call MP_perr_die(myname_,'MPI_comm_size',ier)

	! This is the initial range

  lPE=0
  uPE=nPE-1

  if(uPE <= lPE) return	! There is only one PE in this communicator

  call MP_comm_rank(comm,myPE,ier)
	if(ier /= 0) call MP_perr_die(myname_,'MPI_comm_rank',ier)

	! Check indices

  icost=AttrVect_indexRA(Matx,'cost')

	if(icost<=0) then
	  write(stderr,'(2a,i4)') myname_,': invalid, icost =',icost
	  call die(myname_)
	endif

	rAttr => ptr_rAttr(Matx)

	! Distribute the costs through an indexed recursive bi-section
	! prorcess

  lb = 1
  ub = nMatx

! Create a table of partial sums for partitioning the workload.
! -------------------------------------------------------------
  Allocate(acost(0:nmatx),STAT=ier)
      If (ier /= 0) Call die(myname_,'allocate()',ier)
      If (mall_ison()) Call mall_mci(acost,myname)

  acost(0) = 0
  Do i=1,nMatx
     acost(i) = acost(i-1) + rAttr(icost,iMatx(i))
  End Do
		nullify(rAttr) ! no longer needed

  Call recurPart(acost,lb,ub,myCost, myPE,lPE,uPE)
  nMatx = ub-lb+1
  iMatx(1:nMatx) = iMatx(lb:ub)

     If (mall_ison()) Call mall_mco(acost,myname)
     Deallocate(acost)

end subroutine part_

end module m_rectPartMatx
!.
