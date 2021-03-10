!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_scatPartMatx - a cost-weighted parallel partitioner
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_scatPartMatx
      implicit none
      private	! except

      public :: PartMatx		! The class data structure

      interface PartMatx; module procedure part_; end interface

! !REVISION HISTORY:
! 	18Mar99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_scatPartMatx'

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

  integer :: nPE,lPE,uPE,myPE
  integer :: ier,i,l
  integer :: icost

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

	! Sort the order according to their costs.

		rAttr => ptr_rAttr(Matx)

  call IndexSort(nMatx,iMatx,rAttr(icost,:),descend=.true.)

	! Distribute the costs through an indexed recursive bi-section
	! prorcess

  call recurPart(rAttr(icost,:),iMatx,nMatx,myCost, myPE,lPE,uPE)

		nullify(rAttr)

end subroutine part_

end module m_scatPartMatx
!.
