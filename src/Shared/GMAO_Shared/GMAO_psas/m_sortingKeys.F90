!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_sortingKeys - Keys for a sequence of sortings
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_sortingKeys
      implicit none
      private	! except

      public :: sortingKeys		! The class data structure
      public :: sortingKeys_set
      public :: sortingKeys_get

      public :: sortingKeys_indexedMerge
      public :: sortingKeys_indexedRank
      public :: sortingKeys_indexedSort

      public :: sortingKeys_merge
      public :: sortingKeys_rank
      public :: sortingKeys_sort

      public :: ANINTSORTINGKEY
      public :: AREALSORTINGKEY

    type sortingKeys
      character(len=1)	:: typ
      integer		:: indx
      logical		:: descend
    end type sortingKeys

    interface sortingKeys_set; module procedure		&
	seti_,	&
	setr_
    end interface

    interface sortingKeys_get; module procedure get_; end interface

    interface sortingKeys_merge    ; module procedure	&
	merge_    ; end interface
    interface sortingKeys_indexedMerge; module procedure	&
	indexmerge_; end interface

    interface sortingKeys_rank     ; module procedure	&
	rank_    ; end interface
    interface sortingKeys_indexedRank ; module procedure	&
	indexRank_; end interface

    interface sortingKeys_sort     ; module procedure	&
	sort_     ; end interface
    interface sortingKeys_indexedSort ; module procedure	&
	indexSort_ ; end interface

! !REVISION HISTORY:
! 	06Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_sortingKeys'

  character(len=1),parameter :: ANINTSORTINGKEY = 'I'
  character(len=1),parameter :: AREALSORTINGKEY = 'R'

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: seti_ - set a integer sorting key
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine seti_(keys,ikey,keyType,keyIndex,descend)
      use m_die,  only : die
      implicit none
      type(SortingKeys),dimension(:),intent(out) :: keys
      integer,intent(in) :: ikey
      integer,intent(in) :: keyType
      integer,intent(in) :: keyIndex
      logical,optional,intent(in) :: descend

! !REVISION HISTORY:
! 	06Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::seti_'

  if(ikey<1 .or. ikey>size(keys)) call die(myname_,'invalid ikey',ikey)

  keys(ikey)%typ    =ANINTSORTINGKEY
  keys(ikey)%indx   =keyIndex
  keys(ikey)%descend=.false.
  if(present(descend)) keys(ikey)%descend=descend

end subroutine seti_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: setr_ - set a REAL sorting key
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine setr_(keys,ikey,keyType,keyIndex,descend)
      use m_die,  only : die
      implicit none
      type(SortingKeys),dimension(:),intent(out) :: keys
      integer,intent(in) :: ikey
      real   ,intent(in) :: keyType
      integer,intent(in) :: keyIndex
      logical,optional,intent(in) :: descend

! !REVISION HISTORY:
! 	06Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::setr_'

  if(ikey<1 .or. ikey>size(keys)) call die(myname_,'invalid ikey',ikey)

  keys(ikey)%typ    =AREALSORTINGKEY
  keys(ikey)%indx   =keyIndex
  keys(ikey)%descend=.false.
  if(present(descend)) keys(ikey)%descend=descend

end subroutine setr_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: get_ - get a sorting key set
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine get_(key,keyType,keyIndex,descend)
      use m_chars,only : uppercase
      use m_die,  only : die
      implicit none
      type(SortingKeys),intent(in) :: key
      character(len=1),optional,intent(out) :: keyType
      integer,optional,intent(out) :: keyIndex
      logical,optional,intent(out) :: descend

! !REVISION HISTORY:
! 	06Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::get_'

  if(present(keyType )) keyType =key%typ
  if(present(keyIndex)) keyIndex=key%indx
  if(present(descend )) descend =key%descend

end subroutine get_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: sort_ - sorting with keys
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine sort_(iAttr,rAttr,skeys, iOut,rOut)
      use m_SortingTools,only : indexSet
      use m_die,only : die
      implicit none
      integer,dimension(:,:),intent(in) :: iAttr
      real   ,dimension(:,:),intent(in) :: rAttr
      type(SortingKeys),dimension(:),intent(in) :: skeys
      integer,dimension(:,:),intent(out) :: iOut
      real   ,dimension(:,:),intent(out) :: rOut


! !REVISION HISTORY:
! 	13Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::sort_'
  integer :: nIA,nRA
  integer :: k,l,n
  integer :: ier
  integer,allocatable,dimension(:) :: indx

  n=min(size(iAttr,2),size(rAttr,2))
	allocate(indx(n),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)

  call indexSet(indx)

  call indexsort_(indx,iAttr,rAttr,skeys)

  nIA=min(size(iAttr,1),size(iOut,1))
  do k=1,min(size(iOut,2),n)	! k is the rank, l is the location
    l=indx(k)
    iOut(1:nIA,k)=iAttr(1:nIA,l)
  end do

  nRA=min(size(rAttr,1),size(rOut,1))
  do k=1,min(size(rOut,2),n)	! k is the rank, l is the location
    l=indx(k)
    rOut(1:nRA,k)=rAttr(1:nRA,l)
  end do

	deallocate(indx,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

end subroutine sort_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexsort_ - index sorting with keys
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine indexsort_(indx,iAttr,rAttr,skeys)
      use m_SortingTools,only : indexSort
      use m_die,only : die
      implicit none
      integer,dimension(:)  ,intent(inout) :: indx	! indices
      integer,dimension(:,:),intent(in)    :: iAttr	! int. attrib.
      real   ,dimension(:,:),intent(in)    :: rAttr	! real attrib.
      type(SortingKeys),dimension(:),intent(in) :: skeys ! sortingKeys

! !REVISION HISTORY:
! 	10Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::indexsort_'
  integer :: k
  character(len=1) :: keyType
  integer :: keyIndex
  logical :: descend

	! The order of sorting-keys is from high to low.  For example,
	! the sortingKeys for a sorting of dates (year-month-day) will
	! be specified as (/iyear,imonth,iday/).  For the actual
	! sorting process, the index-sorting will be made in the order
	! of (day, month, year), from low to high.

  do k=size(skeys),1,-1
    call get_(skeys(k),keyType,keyIndex,descend)

    select case(keyType)
    case(ANINTSORTINGKEY)
      call indexSort(indx,iAttr(keyIndex,:),descend=descend)
    case(AREALSORTINGKEY)
      call indexSort(indx,rAttr(keyIndex,:),descend=descend)
    case default
      call die(myname_,'invalid sortingKey type ("'//keyType//'")',k)
    end select
  end do

end subroutine indexsort_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: merge_ - merge two arrays with multiple keys
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine merge_(iAttr,rAttr,		&
	iAttr_i,rAttr_i,iAttr_j,rAttr_j, mkeys	)

      use m_SortingTools, only : rankSet,rankMerge
      use m_die,only : die

      implicit none

      integer,dimension(:,:),intent(out) :: iAttr
      real   ,dimension(:,:),intent(out) :: rAttr

      integer,dimension(:,:),intent(in ) :: iAttr_i
      real   ,dimension(:,:),intent(in ) :: rAttr_i
      integer,dimension(:,:),intent(in ) :: iAttr_j
      real   ,dimension(:,:),intent(in ) :: rAttr_j

      type(SortingKeys),dimension(:),intent(in) :: mkeys

! !REVISION HISTORY:
! 	13Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//':merge_'

  integer :: ni,nj
  integer,allocatable,dimension(:) :: krank_i
  integer,allocatable,dimension(:) :: krank_j
  character(len=1) :: keyType
  integer :: keyIndex
  logical :: descend
  integer :: i,j,k
  integer :: ier

  ni=size(iAttr_i)
  nj=size(iAttr_j)

	allocate(krank_i(ni),krank_j(nj),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)

  call rankSet(krank_i)
  call rankSet(krank_j)

  do k=1,size(mkeys)
    call get_(mkeys(k),keyType,keyIndex,descend)

    select case(keyType)
    case(ANINTSORTINGKEY)

      call rankMerge( iAttr_i(keyIndex,:),iAttr_j(keyIndex,:),	&
		      krank_i(:),krank_j(:), descend=descend	)

    case(AREALSORTINGKEY)

      call rankMerge( rAttr_i(keyIndex,:),rAttr_j(keyIndex,:),	&
		      krank_i(:),krank_j(:), descend=descend	)

    case default
      call die(myname_,'invalid sortingKey type ("'//keyType//'")',k)
    end select
  end do

  call rankMerge(krank_i,krank_j)

	! Physically realize the sorting, by placing the two input
	! vectors into their locations.

  do i=1,ni
    k=krank_i(i)	! k is the rank, i is the location in i-vec
    iAttr(:,k)=iAttr_i(:,i)
    rAttr(:,k)=rAttr_i(:,i)
  end do

  do j=1,nj
    k=krank_j(j)	! k is the rank, j is the location in j-vec
    iAttr(:,k)=iAttr_j(:,j)
    rAttr(:,k)=rAttr_j(:,j)
  end do

	deallocate(krank_i,krank_j,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

end subroutine merge_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexmerge_ - index-merge two segments with multiple keys
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine indexmerge_(indx,indx_i,indx_j,iAttr,rAttr, mkeys)

      use m_SortingTools, only : rankSet,rankMerge,IndexedRankMerge
      use m_die,only : die

      implicit none

      integer,dimension(:),intent(out) :: indx		! of i&j segs

      integer,dimension(:),intent(in ) :: indx_i	! of i-seg
      integer,dimension(:),intent(in ) :: indx_j	! of j-seg

      integer,dimension(:,:),intent(in ) :: iAttr	! of the full
      real   ,dimension(:,:),intent(in ) :: rAttr	! of the full

      type(SortingKeys),dimension(:),intent(in) :: mkeys

! !REVISION HISTORY:
! 	13Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//':indexmerge_'

  integer :: ni,nj
  integer,allocatable,dimension(:) :: krank_i
  integer,allocatable,dimension(:) :: krank_j
  character(len=1) :: keyType
  integer :: keyIndex
  logical :: descend
  integer :: i,j,k
  integer :: ier

  ni=size(indx_i)
  nj=size(indx_j)

	allocate(krank_i(ni),krank_j(nj),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)

  call rankSet(krank_i)
  call rankSet(krank_j)

  do k=1,size(mkeys)
    call get_(mkeys(k),keyType,keyIndex,descend)

    select case(keyType)
    case(ANINTSORTINGKEY)

      call IndexedRankMerge( indx_i(:), indx_j(:),iAttr(keyIndex,:),&
			    krank_i(:),krank_j(:),descend=descend)

    case(AREALSORTINGKEY)

      call IndexedRankMerge( indx_i(:), indx_j(:),rAttr(keyIndex,:),&
			      krank_i(:),krank_j(:),descend=descend)

    case default
      call die(myname_,'invalid sortingKey type ("'//keyType//'")',k)
    end select
  end do

  call rankMerge(krank_i,krank_j)

	! Define the index array of the merged segment of the two, by
	! placing the input indices into their new locations.  Note
	! that

  do i=1,ni
    k=krank_i(i)	! k is the rank, i is the location in i-vec
    indx(k)=indx_i(i)
  end do

  do j=1,nj
    k=krank_j(j)	! k is the rank, j is the location in j-vec
    indx(k)=indx_j(j)
  end do

	deallocate(krank_i,krank_j,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

end subroutine indexmerge_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: rank_ - rank data according the given keys
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine rank_(ranks,iAttr,rAttr,rkeys)
      use m_SortingTools,only : rankSet,rankMerge
      use m_die,only : die
      implicit none
      integer,dimension(:),  intent(out):: ranks
      integer,dimension(:,:),intent(in) :: iAttr
      real   ,dimension(:,:),intent(in) :: rAttr
      type(SortingKeys),dimension(:),intent(in) :: rkeys

! !REVISION HISTORY:
! 	23Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::rank_'

  integer,dimension(0) :: ranks_n
  integer,dimension(size(iAttr,1),0) :: iAttr_n
  real   ,dimension(size(rAttr,1),0) :: rAttr_n

  character(len=1) :: keyType
  integer :: keyIndex
  logical :: descend
  integer :: k
  integer :: ier

  call rankSet(ranks)

  do k=1,size(rkeys)
    call get_(rkeys(k),keyType,keyIndex,descend)

    select case(keyType)
    case(ANINTSORTINGKEY)

      call rankMerge( iAttr(keyIndex,:),iAttr_n(keyIndex,:),	&
		      ranks(:),ranks_n(:), descend=descend	)

    case(AREALSORTINGKEY)

      call rankMerge( rAttr(keyIndex,:),rAttr_n(keyIndex,:),	&
		      ranks(:),ranks_n(:), descend=descend	)

    case default
      call die(myname_,'invalid sortingKey type ("'//keyType//'")',k)
    end select
  end do

end subroutine rank_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexRank_ - rank data according the given indexed keys
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine indexRank_(ranks,indx,iAttr,rAttr,rkeys)
      use m_SortingTools,only : rankSet,indexedRankMerge
      use m_die,only : die
      implicit none
      integer,dimension(:),  intent(out):: ranks
      integer,dimension(:),  intent(in) :: indx
      integer,dimension(:,:),intent(in) :: iAttr
      real   ,dimension(:,:),intent(in) :: rAttr
      type(SortingKeys),dimension(:),intent(in) :: rkeys

! !REVISION HISTORY:
! 	23Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::indexRank_'

  integer,dimension(0) :: ranks_n
  integer,dimension(0) ::  indx_n

  character(len=1) :: keyType
  integer :: keyIndex
  logical :: descend
  integer :: k
  integer :: ier

  call rankSet(ranks)

  do k=1,size(rkeys)
    call get_(rkeys(k),keyType,keyIndex,descend)

    select case(keyType)
    case(ANINTSORTINGKEY)

      call IndexedRankMerge( indx(:), indx_n(:),iAttr(keyIndex,:),&
			    ranks(:),ranks_n(:),descend=descend)

    case(AREALSORTINGKEY)

      call IndexedRankMerge( indx(:), indx_n(:),rAttr(keyIndex,:),&
			    ranks(:),ranks_n(:),descend=descend)

    case default
      call die(myname_,'invalid sortingKey type ("'//keyType//'")',k)
    end select
  end do

end subroutine indexRank_

end module m_sortingKeys
