!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_Distribution - Distribution of a vector class
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_Distribution
      use m_Collector ,only : Collector
      use m_Transposer,only : Transposer
      implicit none
      private	! except

      public :: Distribution		! Class data structure
      public :: Distribution_init ,init		! initialize
      public :: Distribution_clean,clean	! clean

      public :: localSize		! local size
      public :: globalSize		! global size
      public :: distrSize		! distributed size
      public :: new
      public :: delete

		! Advanced interfaces for "friends"

      public :: ptr_aindx		! refering to %aindx
      public :: ptr_prank		! refering to %prank
      public :: ptr_coll		! refering to %coll
      public :: ptr_tran		! refering to %tran
      
      type Distribution
	private
				! permutation/unpermutation indices

        integer,pointer,dimension(:) :: aindx ! pre-distribution
        integer,pointer,dimension(:) :: prank ! post-distribution

				! Message-passing patterns

        type(Collector) ,pointer :: coll	! allgather
        type(Transposer),pointer :: tran	! alltoall
      end type Distribution

      interface Distribution_init; module procedure	&
	initbyAttr_,	&
	initbyDict_; end interface
      interface init; module procedure	&
	initbyAttr_,	&
	initbyDict_; end interface
        
      interface Distribution_clean; module procedure	&
	clean_; end interface
      interface clean; module procedure clean_; end interface

      interface  localSize; module procedure lsize_; end interface
      interface globalSize; module procedure gsize_; end interface
      interface  distrSize; module procedure dsize_; end interface

      interface    new; module procedure    new_; end interface
      interface delete; module procedure delete_; end interface

      interface ptr_aindx; module procedure ptr_aindx_; end interface
      interface ptr_prank; module procedure ptr_prank_; end interface

      interface ptr_coll; module procedure ptr_coll_; end interface
      interface ptr_tran; module procedure ptr_tran_; end interface
        
! !REVISION HISTORY:
! 	09Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_Distribution'

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initbyAttr_ - Initialize an object and other by-products
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine initbyAttr_(dstr,attr_recv,man,keyList,attr_send, &
	comm,wpower,nozero)

      use m_die ,only : die
      use m_mall,only : mall_ison,mall_ci,mall_mci
      use m_Attributes,only : Attributes
      use m_Attributes,only : lsize
      use m_Attributes,only : permute
      use m_Attributes,only : unpermute
      use m_Attributes,only : clean
      use m_AttributesComm,only : trans
      use m_AttributesMAN ,only : build_fromAttr

      use m_MultiAccessNavigator,only : MultiAccessNavigator

      use m_keyChain,only : keyChain
      use m_keyChain,only : keyChain_init
      use m_keyChain,only : clean

      use m_Partitioner,only : Partitioner
      use m_Partitioner,only : Partitioner_init
      use m_Partitioner,only : Partitioner_partition
      use m_Partitioner,only : Partitioner_clean

      use m_Transposer,only : recvSize
      use m_Transposer,only : new
      use m_Collector ,only : new

      implicit none

			! targetted output and by-products

      type(Distribution),intent(out) :: dstr	    ! targetted object
      type(Attributes)  ,intent(out) :: attr_recv   ! distributed data
      type(MultiAccessNavigator),intent(out) :: man ! block Navigators

      integer,dimension(:),intent(in) :: keyList   ! partition keys

			! keyList gives one or more keys to specify the
			! partition strategy of the given attributes.

      type(Attributes)    ,intent(in) :: attr_send ! being distributed
      integer             ,intent(in) :: comm	   ! communicator

			! The partitioner may be defined by balancing
			! the (wpower) power of the counts.

      real   ,optional,intent(in) :: wpower

			! The partitioner may be required to make sure
			! there is no PE left with no data.

      logical,optional,intent(in) :: nozero

! !REVISION HISTORY:
! 	09Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::initbyAttr_'
  integer :: lsend,lrecv
  integer :: ier

  type(Partitioner) :: ptnr
  type(keyChain)    :: keyCh
  type(Attributes)  :: sattr	! send buffer
  type(Attributes)  :: rattr	! recv buffer

			! Allocate new components

  dstr%tran => new(dstr%tran)
  dstr%coll => new(dstr%coll)

			! Allocate the memory for pre-distribution
			! sorting indices.

  lsend=lsize(attr_send)
  allocate(dstr%aindx(lsend), STAT = ier)
        if(ier/=0) call die(myname_,'allocate(%aindx)',ier)
        if(mall_ison()) call mall_mci(dstr%aindx,myname)

			! Build a partitioner with a given keyChain.
			! This partitioner defines how the given
			! attributes should be distributed.


	call keyChain_init(keyCh,keyList)
	call Partitioner_init(ptnr,comm,attr_send,keyCh,	&
		wpower=wpower,nozero=nozero)

			! Use the partitioner to define the distribution
			! of the given attributes.

  call Partitioner_partition(attr_send,      &
        dstr%aindx,dstr%tran,dstr%coll, ptnr,comm)

			! Clean the partitioner, and the keyChain
			! defined for one or more recursive building.

	call Partitioner_clean(ptnr)
	call clean(keyCh)

	                ! Allocate the memory for post-distribution
			! unpermutation indices.

  lrecv=recvSize(dstr%tran)
  allocate(dstr%prank(lrecv),stat=ier)
	if(ier/=0) call die(myname_,'allocate(%prank)',ier)
        if(mall_ison()) call mall_mci(dstr%prank,myname)

	! Use a m_AttributesComm interface that creates a distributed
	! version of the given Attributes.

  call permute(attr_send,sattr,dstr%aindx)

  call trans(sattr,rattr,dstr%tran,comm,stat=ier)
	if(ier/=0) call die(myname_,'trans()',ier)

	! Use a m_AttributesMAN interface to create a post-distribution
	! sorting index and a MultiAccessNavigator (MAN) of a given
	! distributed Attributes.  Note that the sorting order for the
	! sorted output Attributes is not arbitrary if one wants to
	! access them through the derived MAN.  Therefore, there is _no_
	! argument for the sorting order.
	
  call build_fromAttr(man,dstr%prank,rattr)

  call setrank_((dstr%prank(:)),dstr%prank)

	! Use a m_Attributes interface to create a sorted version of
	! the given Attributes.

  call unpermute(rattr,attr_recv,dstr%prank)

	call clean(sattr)
	call clean(rattr)

end subroutine initbyAttr_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initbyDict_ - Initialize an object and other by-products
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine initbyDict_(dstr,attr_recv,man, attr_send,	&
	dict,attr_dict, comm,wpower,nozero)

      use m_die ,only : die
      use m_mall,only : mall_ison,mall_ci,mall_mci
      use m_Attributes,only : Attributes
      use m_Attributes,only : lsize
      use m_Attributes,only : permute
      use m_Attributes,only : unpermute
      use m_Attributes,only : clean
      use m_AttributesComm,only : trans
      use m_AttributesMAN ,only : build_fromAttr

      use m_MultiAccessNavigator,only : MultiAccessNavigator

      use m_Dictionary,only : Dictionary

      use m_Partitioner,only : Partitioner
      use m_Partitioner,only : Partitioner_init
      use m_Partitioner,only : Partitioner_partition
      use m_Partitioner,only : Partitioner_clean

      use m_Transposer,only : recvSize
      use m_Transposer,only : new
      use m_Collector ,only : new

      implicit none

			! targetted output and by-products

      type(Distribution),intent(out) :: dstr	    ! targetted object
      type(Attributes)  ,intent(out) :: attr_recv   ! distributed data
      type(MultiAccessNavigator),intent(out) :: man ! block Navigators

      type(Attributes),target,intent(in) :: attr_send

			! A pre-determined Dictionary has a Dictionary
			! and an Attributes.

      type(Dictionary),target,intent(in) :: dict
      type(Attributes),target,intent(in) :: attr_dict

      integer             ,intent(in) :: comm	   ! communicator

			! The partitioner may be defined by balancing
			! the (wpower) power of the counts.

      real   ,optional,intent(in) :: wpower

			! The partitioner may be required to make sure
			! there is no PE left with no data.

      logical,optional,intent(in) :: nozero

! !REVISION HISTORY:
!	31Oct01	- Jing Guo
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::initbyDict_'
  integer :: lsend,lrecv
  integer :: ier

  type(Partitioner) :: ptnr
  type(Attributes)  :: sattr	! send buffer
  type(Attributes)  :: rattr	! recv buffer

			! Allocate new components

  dstr%tran => new(dstr%tran)
  dstr%coll => new(dstr%coll)

			! Allocate the memory for pre-distribution
			! sorting indices.

  lsend=lsize(attr_send)
  allocate(dstr%aindx(lsend), STAT = ier)
        if(ier/=0) call die(myname_,'allocate(%aindx)',ier)
        if(mall_ison()) call mall_mci(dstr%aindx,myname)

			! Build a partitioner with a pre-determined
			! Dictionary.  A partitioner defines how the
			! given attributes should be distributed.

	call Partitioner_init(ptnr,comm,dict,attr_dict,	&
		wpower=wpower,nozero=nozero)

			! Use the partitioner to define the distribution
			! of the given attributes.

  call Partitioner_partition(attr_send,      &
        dstr%aindx,dstr%tran,dstr%coll, ptnr,comm)

			! Clean the partitioner

	call Partitioner_clean(ptnr)

	                ! Allocate the memory for post-distribution
			! unpermutation indices.

  lrecv=recvSize(dstr%tran)
  allocate(dstr%prank(lrecv),stat=ier)
	if(ier/=0) call die(myname_,'allocate(%prank)',ier)
        if(mall_ison()) call mall_mci(dstr%prank,myname)

	! Use a m_AttributesComm interface that creates a distributed
	! version of the given Attributes.

  call permute(attr_send,sattr,dstr%aindx)

  call trans(sattr,rattr,dstr%tran,comm,stat=ier)
	if(ier/=0) call die(myname_,'trans()',ier)

	! Use a m_AttributesMAN interface to create a post-distribution
	! sorting index and a MultiAccessNavigator (MAN) of a given
	! distributed Attributes.  Note that the sorting order for the
	! sorted output Attributes is not arbitrary if one wants to
	! access them through the derived MAN.  Therefore, there is _no_
	! argument for the sorting order.
	
  call build_fromAttr(man,dstr%prank,rattr)

  call setrank_((dstr%prank(:)),dstr%prank)

	! Use a m_Attributes interface to create a sorted version of
	! the given Attributes.

  call unpermute(rattr,attr_recv,dstr%prank)

	call clean(sattr)
	call clean(rattr)

end subroutine initbyDict_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: setrank_ - convert indx to rank
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine setrank_(indx,rank)
      implicit none
      integer,dimension(:),intent(in) :: indx
      integer,dimension(:),intent(out) :: rank

! !REVISION HISTORY:
! 	06Nov00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::setrank_'
  integer :: i,l

  do i=1,size(indx)		! i is for the ranks of the sorted data
    l=indx(i)			! l is where the entry is
    rank(l)=i
  end do

end subroutine setrank_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - Clean an object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(dstr)
      use m_die ,only : die
      use m_mall,only : mall_ison,mall_mco
      use m_Collector ,only :  Collector_clean,delete
      use m_Transposer,only : Transposer_clean,delete

      implicit none
      type(Distribution),intent(inout) :: dstr

! !REVISION HISTORY:
! 	09Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ier
  
  call Collector_clean(dstr%coll)
  call delete(dstr%coll)

  call Transposer_clean(dstr%tran)
  call delete(dstr%tran)

        if(mall_ison()) then
          call mall_mco(dstr%aindx,myname)
          call mall_mco(dstr%prank,myname)
        endif
        
  deallocate(dstr%aindx,dstr%prank,stat=ier)
        if(ier/=0) call die(myname_,'deallocate()',ier)
end subroutine clean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: lsize_ - local size of a distribution
!
! !DESCRIPTION:
!
! !INTERFACE:

    function lsize_(dstr)
      use m_Collector,only : localSize
      implicit none
      type(Distribution),intent(in) :: dstr
      integer :: lsize_

! !REVISION HISTORY:
! 	11Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::lsize_'

  lsize_ = localSize(dstr%coll)

end function lsize_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: gsize_ - global size of a distribution
!
! !DESCRIPTION:
!
! !INTERFACE:

    function gsize_(dstr)
      use m_Collector,only : globalSize
      implicit none
      type(Distribution),intent(in) :: dstr
      integer :: gsize_

! !REVISION HISTORY:
! 	11Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::gsize_'

  gsize_=globalSize(dstr%coll)

end function gsize_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: dsize_ - distributed size of a distribution
!
! !DESCRIPTION:
!
! !INTERFACE:

    function dsize_(dstr)
      use m_Transposer,only : recvSize
      implicit none
      type(Distribution),intent(in) :: dstr
      integer :: dsize_

! !REVISION HISTORY:
! 	11Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::dsize_'

  dsize_ = recvSize(dstr%tran)

end function dsize_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_aindx_ - refering to %aindx
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_aindx_(dstr)
      implicit none
      type(Distribution),intent(in) :: dstr
      integer,pointer,dimension(:) :: ptr_aindx_

! !REVISION HISTORY:
! 	11Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_aindx_'

  ptr_aindx_ => dstr%aindx

end function ptr_aindx_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_prank_ - refering to %prank
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_prank_(dstr)
      implicit none
      type(Distribution),intent(in) :: dstr
      integer,pointer,dimension(:) :: ptr_prank_

! !REVISION HISTORY:
! 	11Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_prank_'

  ptr_prank_ => dstr%prank

end function ptr_prank_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_tran_ - refering to %tran
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_tran_(dstr)
      use m_Transposer,only : Transposer
      implicit none
      type(Distribution),intent(in) :: dstr
      type(Transposer),pointer :: ptr_tran_

! !REVISION HISTORY:
! 	11Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_tran_'

  ptr_tran_ => dstr%tran

end function ptr_tran_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_coll_ - refering to %coll
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_coll_(dstr)
      use m_Collector,only : Collector
      implicit none
      type(Distribution),intent(in) :: dstr
      type(Collector),pointer :: ptr_coll_

! !REVISION HISTORY:
! 	11Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_coll_'

  ptr_coll_ => dstr%coll

end function ptr_coll_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: new_ - create an object for a pointer
!
! !DESCRIPTION:
!
! !INTERFACE:

    function new_(mold,stat)
      use m_die, only : die
      use m_mall,only : mall_ison,mall_ci
      implicit none
      type(Distribution),pointer :: mold
      integer,optional,intent(out) :: stat
      type(Distribution),pointer :: new_

! !REVISION HISTORY:
! 	28Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::new_'
  type(Distribution),pointer :: obj
  integer :: ier

  if(present(stat)) stat=0
  allocate(obj,stat=ier)
	if(ier/=0) then
	  if(.not.present(stat)) call die(myname_,'allocate()',ier)
	  stat=ier
	  return
	endif

	if(mall_ison()) call mall_ci(1,myname)

  new_ => obj
  nullify(obj)		! to prevent the compiler touching the memory.
end function new_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: delete_ - delete an object for a pointer
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine delete_(obj,stat)
      use m_die, only : die
      use m_mall,only : mall_ison,mall_co
      implicit none
      type(Distribution),pointer :: obj
      integer,optional,intent(out) :: stat


! !REVISION HISTORY:
! 	28Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::delete_'
  integer :: ier

  if(present(stat)) stat=0

	if(mall_ison()) call mall_co(1,myname)

  deallocate(obj,stat=ier)
	if(ier/=0) then
	  if(.not.present(stat)) call die(myname_,'deallocate()',ier)
	  stat=ier
	  return
	endif

end subroutine delete_

end module m_Distribution
