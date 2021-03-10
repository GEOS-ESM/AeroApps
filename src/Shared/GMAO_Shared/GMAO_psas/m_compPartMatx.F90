!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_compPartMatx - a parallel recursive (Matx) partitioner
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_compPartMatx
      implicit none
      private	! except

      public :: PartMatx		! The class data structure

      interface PartMatx; module procedure part_; end interface

! !REVISION HISTORY:
! 	18Mar99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_compPartMatx'

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: part_ - a parallele recursive (Matx) partitioner
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine part_(Matx, iMatx,nMatx, myCost, comm)
      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_indexIA => indexIA
      use m_AttrVect, only : AttrVect_indexRA => indexRA
      use m_mpif90,   only : MP_comm_size,MP_comm_rank
      use m_stdio,    only : stderr
      use m_die,      only : MP_perr_die
      use m_die,      only : perr_die
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

  integer :: irmin,irmax
  integer :: icmin,icmax
  integer :: nPE,lPE,uPE,myPE
  integer :: ier
  integer :: irow,icol,icost
  integer :: ncost
  real,allocatable,dimension(:) :: acost

  if(nMatx <= 0) return	! No work left to _part

  call MP_comm_size(comm,nPE,ier)
	if(ier /= 0) call MP_perr_die(myname_,'MPI_comm_size',ier)

	! This is the initial range

  lPE=0
  uPE=nPE-1

  if(uPE <= lPE) return	! There is only one PE in this communicator

  call MP_comm_rank(comm,myPE,ier)
	if(ier /= 0) call MP_perr_die(myname_,'MPI_comm_rank',ier)

	! Check indices

  irow =AttrVect_indexIA(Matx,'irow')
  icol =AttrVect_indexIA(Matx,'icol')
  icost=AttrVect_indexRA(Matx,'cost')

  if(irow<=0 .or. icol<=0 .or. icost<=0) then
    if(irow <=0) write(stderr,'(2a,i4)') myname_,	&
	': invalid Matx, irow =',irow
    if(icol <=0) write(stderr,'(2a,i4)') myname_,	&
	': invalid Matx, icol =',icol
    if(icost<=0) write(stderr,'(2a,i4)') myname_,	&
	': invalid Matx, icost =',icost
    write(stderr,'(2a)') myname_,': invalid argument, Matx'
    call die(myname_)
  endif

	! Determine the row/column ranges for the cost matrix

  call MatxRanges_(irmin,irmax,icmin,icmax,	&
	iMatx,nMatx, irow,icol,Matx		)

  ncost=(irmax-irmin+1)*(icmax-icmin+1)

		! allocate the workspace of cost-matrix

	allocate(acost(-1:ncost-1),stat=ier)
	if(ier /= 0) call perr_die(myname_,'allocate()',ncost/100)

	! Bisection with a workspace defined

  call recPart_(Matx,irow,icol,icost,		&
	iMatx,nMatx, myCost, myPE,lPE,uPE,	&
	irmin,irmax,icmin,icmax, acost		)

		! deallocate the workspace of cost-matrix

	deallocate(acost,stat=ier)
	if(ier /= 0) call perr_die(myname_,'deallocate()',ier)

end subroutine part_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: MatxRanges_ - check the row/column ranges of Matx
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine MatxRanges_(irmin,irmax,icmin,icmax,	&
	iMatx,nMatx, irow,icol,Matx)

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : ptr_iAttr
      implicit none

      integer,intent(out) :: irmin,irmax
      integer,intent(out) :: icmin,icmax

      integer,dimension(:),intent(in) :: iMatx
      integer,       intent(in) :: nMatx

      integer,       intent(in) :: irow,icol
      type(AttrVect),intent(in) :: Matx

! !REVISION HISTORY:
! 	18Mar99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::MatxRanges_'
  integer :: i,l,ir,ic
  integer,pointer,dimension(:,:) :: iAttr

  l=iMatx(1)

	iAttr => ptr_iAttr(Matx)

  ir=iAttr(irow,l)
  ic=iAttr(icol,l)

  irmin=ir
  irmax=ir

  icmin=ic
  icmax=ic

  do i=2,nMatx
    l=iMatx(i)
    ir=iAttr(irow,l)
    ic=iAttr(icol,l)

    if(ir<irmin) irmin=ir
    if(ir>irmax) irmax=ir

    if(ic<icmin) icmin=ic
    if(ic>icmax) icmax=ic
  end do

	nullify(iAttr)

end subroutine MatxRanges_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: recPart_ - Parallel load balancing with a simple bisection
!
! !DESCRIPTION:
!
!	The algorithm used in this subroutine is recursive.  It 
!   starts from impirisly-serial, and ends to an empirisly-parallel
!   situation.
!
!	This algorithm also allow to be used in a communicator with
!   non-power-of-two number of PEs.
!
!	The costs of the atomic works are pre-defined and stored in
!   Matx.  The accessing order of the works is defined by iMatx.  This
!   will allow a user to experiment different partitioning approaches.
!
! !INTERFACE:

    recursive subroutine recPart_(Matx, irow,icol,icost,	&
	iMatx,nMatx,myCost, myPE,lPE,uPE,			&
	irmin,irmax,icmin,icmax, acost		)

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : ptr_iAttr
      use m_AttrVect, only : ptr_rAttr
      use m_die,      only : perr_die
      implicit none

      type(AttrVect),intent(in) :: Matx	! a list of matrix blocks
      integer,intent(in) :: irow,icol,icost

      integer,dimension(:),intent(inout) :: iMatx ! indices to Matx
      integer,             intent(inout) :: nMatx ! <=size(iMatx)
      real,                intent(inout) :: myCost

      integer,intent(in) :: myPE	! my PE ID.
      integer,intent(in) :: lPE,uPE	! The range of all PEs

      integer,intent(in) :: irmin,irmax,icmin,icmax
      real,dimension(-1:) :: acost	! workspace

! !REVISION HISTORY:
! 	17Mar99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::recPart_'
  real :: half,halfcost
  integer :: ir,ic
  integer :: irmn,irmx
  integer :: icmn,icmx
  integer :: lbPE,ubPE
  integer :: mPE
  integer :: ir_range,ic_range
  integer :: ncost
  logical :: colMajor
  integer :: i,k,l,lsm,llg,lmed,lrc
  integer :: ier

  integer,pointer,dimension(:,:) :: iAttr
  real   ,pointer,dimension(:,:) :: rAttr

	! The _end-recursive_ condition: 
	!
	!	if there is only ONE or less PE left. i.e. uPE<=lPE
	!	if there is no local work-loads left. i.e. nMatx<=0

  if(uPE <= lPE) return
  if(nMatx <= 0) return

	! Compute the middle-PE.  The range [lPE:uPE] will be
	! devided into two sub-ranges [lPE:mPE] and [mPE+1:uPE]

  mPE=lPE+(uPE-lPE)/2

	! Please verify the conditions of:
	!
	!   1) uPE==lPE+1 => mPE==lPE+0==lPE, mPE+1==uPE
	!   2) uPE==lPE   => mPE==lPE+0==lPE, mPE+1> uPE

	! Partition the work-loads to the two sub-ranges.  This
	! approach allows a communicator setting with non-power-of-
	! two number of PEs.  i.e. the meaning of _halfing_ is defined
	! by a proper definition of _half_ number of PEs.

  half=real(mPE-lPE+1)/real(uPE-lPE+1)

  ir_range=irmax-irmin+1
  ic_range=icmax-icmin+1

	! The working size of acost(:)

  ncost=ir_range * ic_range
	if(ncost-1 > ubound(acost,1))	&
	  call perr_die(myname_,'ubound(acost)',ubound(acost,1))

	! If the range of (ir) is larger than the range of (ic), 
	! the bisection wil be in the columns

  colMajor= ir_range > ic_range

	! Fill the cost-matrix

  acost(-1:ncost-1)=0.

		iAttr => ptr_iAttr(Matx)
		rAttr => ptr_rAttr(Matx)

  do i=1,nMatx
    l=iMatx(i)
    ir=iAttr(irow,l)-irmin
    ic=iAttr(icol,l)-icmin

	! The default is row-major

    lrc=ic*ir_range + ir
    if(colMajor) lrc=ir*ic_range + ic

    acost(lrc)=rAttr(icost,l)
  end do

		nullify(iAttr)
		nullify(rAttr)

	! Accumulate the cost. acost(-1:0) is already defined

  do i=1,ncost-1
    acost(i)=acost(i-1)+acost(i)
  end do

	! _halfcost_ in term of the total computation cost

  halfcost=half*acost(ncost-1)

	! Binary search for the nearest location of the breaking 
	! point, _halfcost_.

  lsm=-1
  llg=ncost-1

  do while(llg > lsm+1)		! Until the binary search leads to
    lmed=lsm+(llg-lsm)/2	! a single section.

    if(halfcost > acost(lmed)) then
      lsm=lmed
    else
      llg=lmed
    endif
  end do

	! At this point, llg==lsm+1 and acost(llg) >= acost(lsm)

  lmed=lsm
  if( lsm == -1) then
    lmed=llg

  else
	! Recenter the breaking point, if using llg may improve the
	! relative load balance between the two halves.  Note that
	! the equation used below is derived from a simple relation:
	!
	!
	!   If the relative excessive cost of the lower half with
	!   _llg_ as the tbreaking point is less than the relative
	!   excessive cost of the upper half with _lsm_ as the
	!   breaking point, use llg.
	!
	!   Relative excessive cost: (L-H)/H, where L is the actual
	!   load (total-cost), H is the expected _halfcost_ for the
	!   given _half_ (the lower-half or the upper-half).

    if(	acost(llg) * (acost(ncost-1)-halfcost)<	&
	(acost(ncost-1)-acost(lsm)) * halfcost	) lmed=llg

  endif


	! Mark all matrix blocks to two sections

		iAttr => ptr_iAttr(Matx)

  do i=1,nMatx
    l=iMatx(i)
    ir=iAttr(irow,l)-irmin
    ic=iAttr(icol,l)-icmin

	! The default is row-major

    lrc=ic*ir_range + ir
    if(colMajor) lrc=ir*ic_range + ic

    if(lrc > lmed) iMatx(i)=-l
  end do
		nullify(iAttr)

	! Bi-section the work-loads

  if(myPE <= mPE) then

	! Take the lower half matrix blocks.

    k=0
    do i=1,nMatx
      l=iMatx(i)
      if(l > 0) then
	k=k+1
	iMatx(k)=l
      endif
    end do
    nMatx=k		! The new size of iMatx

    lbPE=lPE
    ubPE=mPE

    myCost=acost(lmed)

  else
	! Take the upper half matrix blocks.

    k=0
    do i=1,nMatx
      l=-iMatx(i)
      if(l > 0) then
	k=k+1
	iMatx(k)=l
      endif
    end do
    nMatx=k		! The new size of iMatx

    lbPE=mPE+1
    ubPE=uPE

    myCost=acost(ncost-1)-acost(lmed)

  endif

	! If there is nothing to be done, don't try MatxRanges_()

  if(nMatx<=0 .or. ubPE<=lbPE) return

  call MatxRanges_(irmn,irmx,icmn,icmx,	&
	iMatx,nMatx, irow,icol,Matx	)

  call recPart_(Matx,irow,icol,icost,		&
	iMatx,nMatx, myCost, myPE,lbPE,ubPE,	&
	irmn,irmx,icmn,icmx, acost		)

end subroutine recPart_
end module m_compPartMatx
!.
