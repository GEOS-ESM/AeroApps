!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_HilbertPartMatx - an order-preserving cost-weighted parallel partitioner
!        
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_HilbertPartMatx
      implicit none
      private	! except

      public :: PartMatx		! The class data structure

      interface PartMatx; module procedure part_; end interface

! !REVISION HISTORY:
! 	18Mar99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_HilbertPartMatx'

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
      use m_AttrVect, only : ptr_rAttr, ptr_iAttr
      use m_AttrVect, only : AttrVect_indexRA => indexRA
      use m_AttrVect, only : AttrVect_indexIA => indexIA
      use m_SortingTools,only : IndexSort
      use m_recurPart, only : recurPart
      use m_mpif90,   only : MP_comm_size,MP_comm_rank
      use m_stdio,    only : stderr
      use m_die,      only : MP_perr_die
      use m_die,      only : die, assert_
#include "assert.H"
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
  Integer,pointer,dimension(:,:) :: iAttr

  integer :: nPE,lPE,uPE,myPE
  integer :: ier,i,l
  integer :: icost, lb, ub
  Integer :: ic, ir, icol, irow
  Integer, Allocatable :: idx(:)
  Real, Allocatable :: acost(:)
#ifdef _OPENMP
  Integer :: n_threads
  Integer :: chunk_size
#endif

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
		iAttr => ptr_iAttr(Matx)

  ic=AttrVect_indexIA(Matx,'icol')
  ir=AttrVect_indexIA(Matx,'irow')

  Allocate(idx(1:nMatx), acost(0:nMatx), STAT = ier)
     If (ier /= 0) Call die(myname_,'allocate()',ier)
  Do i = 1, nMatx
     icol   = iAttr(ic,i)
     irow   = iAttr(ir,i)
     idx(i) = Hilbert(ic,ir,nMatx)
  End Do
		nullify(iAttr)

  Call IndexSort(nMatx, iMatx, idx, descend=.false.)
       Deallocate(idx)

#ifdef _OPENMP
!!$	n_threads = OMP_GET_MAX_THREADS()
!!$	If (n_threads > 1) Then
!!$           chunk_size = nmatx/n_threads
!!$           call IndexSort(nMatx,iMatx,(/ (mod(i,chunk_size),i=1,nmatx) /),descend=.true.)
!!$        End If
#endif		

        
!
	! Distribute the costs through an indexed recursive bi-section
	! prorcess
  lb = 1 
  ub = nMatx

  acost(0) = 0
  Do i=1,nMatx
     acost(i) = acost(i-1) + rAttr(icost,iMatx(i))
  End Do
		nullify(rAttr)

  call recurPart(acost,lb,ub,myCost, myPE,lPE,uPE)
       Deallocate(acost)
  
  nMatx = ub - lb + 1
  iMatx(1:nMatx) = iMatx(lb:ub)

		nullify(rAttr)

end subroutine part_

Integer Function Hilbert(icol, irow, n)
  Use m_die, only : assert_

  Integer, Intent(In) :: icol, irow, n

  ! determines translation for next iteration
  Integer, Parameter :: table(3,0:1,0:1,4) = Reshape( (/ &
       & 0,0,2,     1,1,4,    0,1,1,   1,0,1,   &
       & 0,0,1,     0,1,2,    1,1,3,   1,0,2,   &
       & 1,0,3,     0,1,3,    1,1,2,   0,0,4,   &
       & 1,0,4,     1,1,1,    0,1,4,   0,0,3 /), (/ 3,2,2,4 /))
  Integer, Parameter :: ROW = 1, COL = 2, ROT = 3

  Integer :: i, qi, j, qj, nn, s
  Integer :: ir, ii, jj

  i = irow - 1 ! start from zero
  j = icol - 1 ! start from zero
  s = 0
  ALWAYS_ASSERT(n >= 2)
  nn = n/2

  ir = 1
  Do

     qi = i/nn
     qj = j/nn
     i  = mod(i,nn)
     j  = mod(j,nn)

     ii = table(ROW,qi,qj,ir)
     jj = table(COL,qi,qj,ir)
     ir = table(ROT,qi,qj,ir)

     s = 4*s + 2*ii + jj
     ASSERT(s < Huge(1)/4) ! make certain storage is available
     If (nn <= 1) Exit
     nn = (nn+1)/2

  End Do

  Hilbert = s
  
End Function Hilbert

end module m_HilbertPartMatx
!.
