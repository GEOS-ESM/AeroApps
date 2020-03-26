!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_AttrVectSubset - subset an Attribute vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_AttrVectSubset
      implicit none
      private	! except

      public :: AttrVect_init	! The class data structure
      public :: init		! The class data structure

      interface AttrVect_init; module procedure	&
	subset_; end interface
      interface init; module procedure subset_; end interface

! !REVISION HISTORY:
! 	10Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_AttrVectSubset'

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: subset_ - create a subset AttrVect
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine subset_(subset, attr,keySet)
      use m_AttrVect,only : AttrVect
      use m_AttrVect,only : AttrVect_init
      use m_AttrVect,only : lsize
      use m_AttrVect,only : getIList
      use m_AttrVect,only : getRList
      use m_AttrVect,only : ptr_iAttr
      use m_AttrVect,only : ptr_rAttr
      use m_sortingKeys,only : sortingKeys
      use m_sortingKeys,only : sortingKeys_get
      use m_sortingKeys,only : ANINTSORTINGKEY
      use m_sortingKeys,only : AREALSORTINGKEY
      use m_String,only : String
      use m_String,only : clean
      implicit none
      type(AttrVect),intent(out) :: subset
      type(AttrVect),intent(in ) :: attr
      type(sortingKeys),dimension(:),intent(in) :: keySet

! !REVISION HISTORY:
! 	10Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::subset_'
  integer :: kIndx
  character(len=1) :: kType
  integer :: lsize_
  integer :: nIA,nRA,nkeys
  integer :: ier
  integer :: k
  type(String),allocatable,dimension(:)   :: iList
  type(String),allocatable,dimension(:)   :: rList
  integer     ,allocatable,dimension(:)   :: indxIA
  integer     ,allocatable,dimension(:)   :: indxRA
  integer     ,pointer    ,dimension(:,:) :: iAttr
  real        ,pointer    ,dimension(:,:) :: rAttr
  integer     ,pointer    ,dimension(:,:) :: isuba
  real        ,pointer    ,dimension(:,:) :: rsuba

  nkeys=size(keySet)
  allocate(iList(nkeys),rList(nkeys),	&
	indxIA(nkeys),indxRA(nkeys),stat=ier)
	if(ier/=0) call die(myname_,'allocate()',ier)

  lsize_=lsize(attr)
  nIA=0
  nRA=0
  do k=1,nkeys
    call sortingKeys_get(keySet(k),keyIndex=kIndx,keyType=kType)

    select case(kType)
    case(ANINTSORTINGKEY)

      nIA=nIA+1
      call getIList(iList(nIA),kIndx,attr)
      indxIA(nIA)=kIndx

    case(AREALSORTINGKEY)

      nRA=nRA+1
      call getRList(rList(nRA),kIndx,attr)
      indxRA(nRA)=kIndx

    case default
	call die(myname_,'unknwon keyType',kType)
    end select
  end do

  call AttrVect_init(subset,iList(1:nIA),rList(1:nRA),lsize_)

  iAttr => ptr_iAttr(attr)	! pointers as input
  rAttr => ptr_rAttr(attr)

  isuba => ptr_iAttr(subset)	! pointers as output
  rsuba => ptr_rAttr(subset)

  isuba(:,:)=iAttr(indxIA(1:nIA),:)
  rsuba(:,:)=rAttr(indxRA(1:nRA),:)

  nullify(iAttr)
  nullify(rAttr)
  nullify(iSuba)
  nullify(rSuba)

  do k=1,nIA
    call clean(iList(k))
  end do
  do k=1,nRA
    call clean(rList(k))
  end do

	deallocate(iList,rList,indxIA,indxRA,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)
end subroutine subset_
end module m_AttrVectSubset

