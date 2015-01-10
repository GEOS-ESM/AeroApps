!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_recurPart - Elemental parallel recursive partitioner
!
! !DESCRIPTION:
!
!	The subroutines in this module assume all preconditions have
!   been properly configured.  No checking will be done for efficiency
!   reasons.  Therefore, it may be only suitable for lower level
!   programming purposes.
!
! !INTERFACE:

    module m_recurPart
      implicit none
      private	! except

      public :: recurPart		! The class data structure

      interface recurPart; module procedure	&
	iPartR_,	&	! indexed with real arguments
	iPartI_,	&	! indexed with integer arguments
	sPartR_,	&	! segmented with real arguments
	sPartI_			! segmented with integer arguments
      end interface

! !REVISION HISTORY:
! 	18Mar99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_recurPart'

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: iPartR_ - Indexed parallel recursive partitioner
!
! !DESCRIPTION:
!
! !INTERFACE:

    recursive subroutine iPartR_(costs,iCostx,nCostx,myCost,	&
	myPE,lPE,uPE	)

      implicit none

      real,   dimension(:),intent(in)    :: costs
      integer,dimension(:),intent(inout) :: iCostx ! indices to costs
      integer,             intent(inout) :: nCostx ! <=size(iCostx)
      real,                intent(inout) :: myCost

      integer,intent(in) :: myPE	! my PE ID.
      integer,intent(in) :: lPE,uPE	! The range of all PEs

! !REVISION HISTORY:
! 	17Mar99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::iPartR_'
  real :: half,halfcost
  real :: lowerHalf,upperHalf
  real :: lowerLeft,upperLeft
  integer :: lbPE,ubPE
  integer :: mPE
  integer :: i,k,l
  integer :: ier
  logical :: lower

	! The _end-recursive_ condition: 
	!
	!	if there is only ONE or less PE left. i.e. uPE<=lPE
	!	if there is no local work-loads left. i.e. nCostx<=0

  if(uPE <= lPE) return
  if(nCostx <= 0) return

	! Compute the middle-PE.  The range [lPE:uPE] will be
	! divided into two sub-ranges [lPE:mPE] and [mPE+1:uPE]

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

	! _halfwork_ in term of the _total_ cost.  Note that myCost
	! is supposed to be the total-cost at the input

  halfcost=half*myCost

  lowerLeft=       halfcost
  upperLeft=myCost-halfcost

	! Pick the bin with a larger colume

  lower=.true.
  if(lowerLeft/= upperLeft) lower = lowerLeft > upperLeft

	! Recompute what costs are in both bins.

  LowerHalf=0.
  UpperHalf=0.

  do i=1,nCostx
    l=iCostx(i)
    if(.not.lower) iCostx(i)=-l	! Mark the elemement for the upper bin

	! Compute how much room left and the actual costs are in
	! both bins.

    if(lower) then
      lowerLeft=lowerLeft-costs(l)
      lowerHalf=lowerHalf+costs(l)
    else
      upperLeft=upperLeft-costs(l)
      upperHalf=upperHalf+costs(l)
    endif

	! Switch to the candidite bin that has a larger room.  Do not
	! switch if the room left on both bins is the same.

    if(lowerLeft /= upperLeft) lower = lowerLeft > upperLeft
  end do

	! Divide the work-loads to two parts

  if(myPE <= mPE) then

	! Take the _lower-half_ matrix blocks.

    k=0
    do i=1,nCostx
      l=iCostx(i)
      if(l > 0) then
	k=k+1
	if(k<i) iCostx(k)=l
      endif
    end do
    nCostx=k		! The new size of iCostx

    lbPE=lPE
    ubPE=mPE

    myCost=lowerHalf	! Update the total cost of the divided costs(:)

  else
	! Take the _upper-half_ matrix blocks.

    k=0
    do i=1,nCostx
      l=-iCostx(i)
      if(l > 0) then
	k=k+1
	iCostx(k)=l
      endif
    end do
    nCostx=k		! The new size of iCostx

    lbPE=mPE+1
    ubPE=uPE

    myCost=upperHalf	! Update the total cost of the divided costs(:)

  endif

  call iPartR_(costs,iCostx,nCostx,myCost,myPE,lbPE,ubPE)

end subroutine iPartR_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: iPartI_ - A wrapper of iPartR_() with integer arguments
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine iPartI_(costs,iCostx,nCostx,myCost,	&
	myPE,lPE,uPE	)

      implicit none

      integer,dimension(:),intent(in)    :: costs
      integer,dimension(:),intent(inout) :: iCostx ! indices to Costx
      integer,             intent(inout) :: nCostx ! <=size(iCostx)
      integer,             intent(inout) :: myCost

      integer,intent(in) :: myPE	! my PE ID.
      integer,intent(in) :: lPE,uPE	! The range of all PEs

! !REVISION HISTORY:
! 	09Apr99 - Jing Guo <guo@dao> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::iPartI_'

  real :: myCost_

	myCost_=myCost

  call iPartR_(real(costs),iCostx,nCostx,myCost_, myPE,lPE,uPE)

	myCost =nint(myCost_)

end subroutine iPartI_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: sPartR_ - Segmented parallel recursive partitioner
!
! !DESCRIPTION:
!
! !INTERFACE:

    recursive subroutine sPartR_(acost, lbound,ubound, myCost,	&
	myPE,lPE,uPE,nozero)
      use m_die,only : die,perr
      implicit none

      real   ,dimension(0:),intent(in)    :: acost
      integer,              intent(inout) :: lbound,ubound
      real,                 intent(inout) :: myCost

      integer,intent(in) :: myPE	! my PE ID.
      integer,intent(in) :: lPE,uPE	! The range of all PEs
      logical,optional,intent(in) :: nozero	! no empty bin

! !REVISION HISTORY:
! 	17Mar99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::sPartR_'
  real :: half,mid_cost,acost_r
  integer :: lbPE,ubPE
  integer :: mPE
  integer :: i,lsm,llg,lmd

  logical :: nozero_

  nozero_=.false.
  if(present(nozero)) nozero_=nozero

	! The _end_of_recursive_ condition: 
	!
	!	if there is only ONE or less PE left. (uPE<=lPE)
	!	if there is no local work-loads left. (lbound>ubound)

  if(uPE <= lPE) return

  if(nozero_) then
	! Check for an exception, if there are less data than PEs.

    if(ubound-lbound < uPE-lPE) then
      call perr(myname_,'lbnd',lbound,'ubnd',ubound)
      call perr(myname_,'lPE' ,lPE   ,'uPE' ,uPE)
      call die(myname_,'unexpected zero-size')
    endif
  endif

  if(lbound > ubound) return

	! Compute the middle-PE.  The range [lPE:uPE] will be divided
	! into two sub-ranges [lPE:mPE] and [mPE+1:uPE].  This approach
	! allows a communicator setting with non-power-of-two number of
	! PEs.  i.e. the meaning of _halfcost_ is proportional to the
	! definition of _half_ number of PEs.

  mPE=lPE+(uPE-lPE)/2

  half=real(mPE-lPE+1)/real(uPE-lPE+1)

	! _halfcost_ is the middle value of acost(lbound-1:ubound)

  acost_r=acost(lbound-1)
  mid_cost=acost_r + half*(acost(ubound)-acost_r)

	! Binary search for the nearest location of _mid_cost_

  lsm=lbound-1
  llg=ubound

  do while(llg > lsm+1)		! Until the binary search leads to
    lmd=lsm+(llg-lsm)/2		! a single section.

    if(mid_cost > acost(lmd)) then
      lsm=lmd
    else
      llg=lmd
    endif
  end do

	! At this point, llg==lsm+1 and acost(llg) >= acost(lsm)

  lmd=lsm
  if( lsm < lbound ) then
    lmd=llg
  else

	! Tune between the two locations: lsm or llg?  The relation
	! is derived from:
	!
	!   If the relative excessive cost of the lower half with
	!   _llg_ as the tbreaking point is less than the relative
	!   excessive cost of the upper half with _lsm_ as the
	!   breaking point, use llg.
	!
	!   Relative excessive cost: (L-H)/H, where L is the actual
	!   load (total-cost), H is the expected _halfcost_ for the
	!   given _half_ (the lower-half or the upper-half).


    if( (acost(llg)   -acost_r   )*(acost(ubound)-mid_cost) < &
	(acost(ubound)-acost(lsm))*(mid_cost     -acost_r ) ) lmd=llg

  endif

	! Final tuning

  if(nozero_) then

		! Check if any one of two sections could be zero-sized.
		! The purpose is to leave at least one entry on any
		! given PE.  It requires any one of the two sections to
		! have at least the same number of entries as the number
		! of PEs.  Note that the checking is done at the
		! beginning of the next iteration.

    if(lmd-lbound < mPE-lPE) lmd=lbound+mPE-lPE
    if(ubound-lmd < uPE-mPE) lmd=ubound-uPE+mPE

  endif

	! Bi-section the work-loads

  if(myPE <= mPE) then

	! Take the lower half section

    ubound=lmd

    lbPE=lPE
    ubPE=mPE

    myCost=acost(lmd)-acost_r

  else
	! Take the upper half section

    lbound=lmd+1

    lbPE=mPE+1
    ubPE=uPE

    myCost=acost(ubound)-acost(lmd)

  endif

  call sPartR_(acost, lbound,ubound, myCost, myPE,lbPE,ubPE,	&
	nozero=nozero)

end subroutine sPartR_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: sPartI_ - A wrapper of sPartR_ with integer arguments
!
! !DESCRIPTION:
!
! !INTERFACE:

    recursive subroutine sPartI_(acost, lbound,ubound,myCost,	&
	myPE,lPE,uPE,nozero)
      use m_die,only : die,perr
      implicit none
      integer,dimension(0:),intent(in) :: acost
      integer,intent(inout) :: lbound,ubound
      integer,intent(inout) :: myCost

      integer,intent(in) :: myPE,lPE,uPE
      logical,optional,intent(in) :: nozero

! !REVISION HISTORY:
! 	09Apr99 - Jing Guo <guo@dao> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::sPartI_'
  real    :: half
  integer :: mid_cost
  integer :: acost_0
  integer :: lbPE,ubPE
  integer :: mPE
  integer :: i,lsm,llg,lmd
  logical :: nozero_

  nozero_=.false.
  if(present(nozero)) nozero_=nozero

	! The _end_of_recursive_ condition: 
	!
	!	if there is only ONE or less PE left. (uPE<=lPE)
	!	if there is no local work-loads left. (lbound>ubound)

  if(uPE <= lPE) return

  if(nozero_) then
	! Check for an exception, if there are less data than PEs.

    if(ubound-lbound < uPE-lPE) then
      call perr(myname_,'lbnd',lbound,'ubnd',ubound)
      call perr(myname_,'lPE' ,lPE   ,'uPE' ,uPE)
      call die(myname_,'unexpected zero-size')
    endif
  endif

  if(lbound > ubound) return

	! Compute the middle-PE.  The range [lPE:uPE] will be divided
	! into two sub-ranges [lPE:mPE] and [mPE+1:uPE].  This approach
	! allows a communicator setting with non-power-of-two number of
	! PEs.  i.e. the meaning of _halfcost_ is proportional to the
	! definition of _half_ number of PEs.

  mPE=lPE+(uPE-lPE)/2

  half=real(mPE-lPE+1)/real(uPE-lPE+1)

	! _mid_cost_ is the middle value of acost(lbound-1:ubound)

  acost_0=acost(lbound-1)
  mid_cost=acost_0 + half*(acost(ubound)-acost_0)

	! Binary search for the nearest location of _mid_cost_

  lsm=lbound-1
  llg=ubound

  do while(llg > lsm+1)		! Until the binary search leads to
    lmd=lsm+(llg-lsm)/2		! a single section.

    if(mid_cost > acost(lmd)) then
      lsm=lmd
    else
      llg=lmd
    endif
  end do

	! At this point, llg==lsm+1 and acost(llg) >= acost(lsm)

  lmd=lsm
  if( lsm < lbound ) then
    lmd=llg
  else

	! Tune between the two locations: lsm or llg?  The relation
	! is derived from:
	!
	!   If the relative excessive cost of the lower half with
	!   _llg_ as the tbreaking point is less than the relative
	!   excessive cost of the upper half with _lsm_ as the
	!   breaking point, use llg.
	!
	!   Relative excessive cost: (L-H)/H, where L is the actual
	!   load (total-cost), H is the expected _halfcost_ for the
	!   given _half_ (the lower-half or the upper-half).


    if( (acost(llg)   -acost_0   )*(acost(ubound)-mid_cost) < &
	(acost(ubound)-acost(lsm))*(mid_cost     -acost_0 ) ) lmd=llg

  endif

	! Final tuning

  if(nozero_) then

		! Check if any one of two sections could be zero-sized.
		! The purpose is to leave at least one entry on any
		! given PE.  It requires any one of the two sections to
		! have at least the same number of entries as the number
		! of PEs.  Note that the checking is done at the
		! beginning of the next iteration.

    if(lmd-lbound < mPE-lPE) lmd=lbound+mPE-lPE
    if(ubound-lmd < uPE-mPE) lmd=ubound-uPE+mPE

  endif

	! Bi-section the work-loads

  if(myPE <= mPE) then

	! Take the lower half section

    ubound=lmd

    lbPE=lPE
    ubPE=mPE

    myCost=acost(lmd)-acost_0

  else
	! Take the upper half section

    lbound=lmd+1

    lbPE=mPE+1
    ubPE=uPE

    myCost=acost(ubound)-acost(lmd)

  endif

  call sPartI_(acost, lbound,ubound, myCost, myPE,lbPE,ubPE,	&
	nozero=nozero_)

end subroutine sPartI_

end module m_recurPart
!.
