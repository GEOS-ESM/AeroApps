!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_Dictionary - An attribute tree
!
! !DESCRIPTION:
!
!	A Dictionary is used in the following UseCases:
!
!   Usecase 1: build from random data entries in a container.
!
!	. Create an _insertable_ dictionary with a given key
!	definition
!	. Given a container
!	Given a container, each e
!
!   Usecase 2: build from order data entries
!
!     . If an ordered data table is available as a type(Attributes)
!	variable attr, with attr_counts(:) for all entries.
!
!	call Dictionary_init(dict,key,nsize)
!	  dict_counts => ptr_counts(dict)
!	  dict_counts(:) = attr_counts(indx(1:nsize))
!	  nullify(dict_counts)
!
!	call Dictionary_compile(dict)
!
!	do i=1,nsize
!	  call Dictionary_search(i,srched_attr,	&
!		where(i),dict,attr,verify=(i==1))
!	end do
!	[...]
!
!	call Dictionary_clean(dict)
!
! !INTERFACE:

    module m_Dictionary
      implicit none
      private	! except

      public :: Dictionary		! The class data structure
      public :: DictionaryEntry		! secondary data structure

      public :: Dictionary_init
      public :: Dictionary_clean
      public :: Dictionary_insert
      public :: Dictionary_compile

      public :: Dictionary_key
      public :: Dictionary_size
      public :: Dictionary_compiled

      public :: Dictionary_get
      public :: Dictionary_search

      public :: new
      public :: delete
      public :: ptr_counts

    type DictionaryNode
	private
	integer :: count	! number of inserts of the entry
	integer :: irank	! rank of the entry after compiled
	integer :: index	! index to the source container
	type(DictionaryNode),   pointer :: left
	type(DictionaryNode),   pointer :: right
    end type DictionaryNode

    type Dictionary
      private
      integer :: status		! status of this dictionary
      integer :: key		! key for this dictionary

		! compiled (ranked) dictionary contents

      integer,pointer,dimension(:)   :: counts	! of each entry

		! a searchable tree either as insertable or compiled

      integer :: dsize
      type(DictionaryNode),pointer :: node
    end type Dictionary

    type DictionaryEntry
      private
      type(DictionaryNode),pointer :: node
    end type DictionaryEntry

    interface Dictionary_init ; module procedure	&
	initInsertbl_,	&
	initPackable_; end interface
    interface Dictionary_clean; module procedure clean_; end interface

    interface Dictionary_compile; module procedure	&
	compile_; end interface

    interface Dictionary_insert; module procedure insert_; end interface
    interface Dictionary_search; module procedure	&
	search_,	&
	searchi_; end interface

    interface Dictionary_key   ; module procedure key_   ; end interface
    interface Dictionary_size  ; module procedure size_  ; end interface
    interface Dictionary_compiled; module procedure	&
	compiled_; end interface

    interface Dictionary_get; module procedure get_; end interface
    interface ptr_counts; module procedure ptr_counts_; end interface
    interface new   ; module procedure new_   ; end interface
    interface delete; module procedure delete_; end interface

! !REVISION HISTORY:
! 	16Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!	03Apr00	- Jing Guo
!		. Reopened the development on this module to support
!		  the parallel partitioner development
!	06Apr00	- Jing Guo
!		. After several major redesign, this module has been
!		  reimplemented and now tested.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_Dictionary'

  integer,parameter :: status_INSERTBL	= 1
  integer,parameter :: status_UNSRCHBL	= 2
  integer,parameter :: status_COMPILED	= 3

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: insert_ - insert an entry
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine insert_(dict, indx,data, count,where,verify)
      use m_die,only : die
      use m_Attributes,only : Attributes
      use m_Attributes,only : key
      use m_Attributes,only : aSuperset
      implicit none

      type(Dictionary),intent(inout) :: dict	! data dictionary

      integer,intent(in) :: indx		! where in the container
      type(Attributes),intent(in) :: data	! a container

      integer,optional,intent(in) :: count	! weighted counting
      type(DictionaryEntry),optional,intent(out) :: where ! this node

      logical,optional,intent(in) :: verify

! !REVISION HISTORY:
! 	16Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::insert_'
  integer :: kdata,kdict

  if(present(verify)) then
    if(verify) then
      if(dict%status /= status_INSERTBL)	&
	call die(myname_,'not insertable')

      kdict=key_(dict)
      kdata=key(data)

      if(.not.aSuperset(kdata,kdict)) call die(myname_,	&
	'key(dict)',kdict,'key(data)',kdata)
    endif
  endif

	! Insert recursively.  Checking on present(where) is needed
	! because not _where_, but _where%node_ is referenced.

  if(present(where)) then
    call rinsert_(dict%node,dict%key, indx,data,	&
	count=count,where=where%node)
  else
    call rinsert_(dict%node,dict%key, indx,data,	&
	count=count)
  endif

end subroutine insert_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: rinsert_ - insert an entry
!
! !DESCRIPTION:
!
! !INTERFACE:

    recursive subroutine rinsert_(node,key, indx,data, count,where)
      use m_die,only : die
      use m_mall,only : mall_ison,mall_ci
      use m_Attributes,only : Attributes
      use m_Attributes,only : order
      use m_Attributes,only : LESS,GREATER

      implicit none

      type(DictionaryNode),pointer    :: node	! (inout)
      integer,intent(in) :: key		! dictionary key set

      integer,intent(in) :: indx		! where in the container
      type(Attributes),intent(in) :: data	! a container

      integer,optional,intent(in) :: count	! weighted counting
      type(DictionaryNode),optional,pointer :: where	! (out)

! !REVISION HISTORY:
! 	16Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::rinsert_'

  integer :: ier
  integer :: count_
  integer :: iorder

	! if this is a fresh node
  if(.not.associated(node)) then

    allocate(node,stat=ier)
	if(ier/=0) call die(myname_,'allocate(node)',ier)
	if(mall_ison()) call mall_ci(1,myname)

    node%irank=indx	! %irank can not be defined until compile_().
			! However, either the dictionary is INSERTABL
			! (a source container is used) or COMPILED (a
			! compiled container is used), one can always
			! use %irank to local the actual data.

    node%index=indx	! %index stores the index to the source
			! container.
    node%count=1
    if(present(count)) node%count=count

    nullify(node%left)
    nullify(node%right)

    if(present(where)) where => node
    return
  endif

	! If there is a value at this node, compare the value and grow
	! the tree if the values are not the same.

  iorder=order(indx,data,node%irank,data,key)
  select case(iorder)
  case(LESS)	! data(indx) < data(node%irank), left-insert.

    call rinsert_(node%left ,key,indx,data,count=count,where=where)

  case(GREATER)	! data(indx) > data(node%irank), right-insert.

    call rinsert_(node%right,key,indx,data,count=count,where=where)

  case default	! data(indx) == data(node%irank), just count it in.

    if(present(where)) where => node
    count_=1
    if(present(count)) count_=count
    node%count=node%count+count_
  end select

end subroutine rinsert_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: search_ - search for a value
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine search_(indx,data, where,dict,attr,verify)
      use m_die,only : die,perr
      use m_Attributes,only : Attributes
      use m_Attributes,only : key
      use m_Attributes,only : aSuperset
      implicit none
      integer         ,intent(in) :: indx	! data index
      type(Attributes),intent(in) :: data	! data container

      type(DictionaryEntry),intent(out) :: where ! location as a node

      type(Dictionary),intent(in) :: dict	! dictionary
      type(Attributes),intent(in) :: attr	! container for dict

      logical,optional,intent(in) :: verify

! !REVISION HISTORY:
! 	16Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::search_'
  integer :: kdata,kdict,kattr
  logical :: notin_data,notin_attr

  if(present(verify)) then
    if(verify) then

      if( dict%status /= status_INSERTBL	.and.	&
	  dict%status /= status_COMPILED	)	&
		call die(myname_,'not searchable')

      kdata=key(data)
      kdict=key_(dict)
      kattr=key(attr)

      notin_data=.not.aSuperset(kdata,kdict)
      notin_attr=.not.aSuperset(kattr,kdict)

      if(notin_data .or. notin_attr) then

	if(notin_data) call perr(myname_,	&
		'key(data)',kdata,'key(dict)',kdict)
	if(notin_attr) call perr(myname_,	&
		'key(attr)',kattr,'key(dict)',kdict)
	call die(myname_)
      endif
    endif
  endif

  call rsearch_(indx,data, where%node, dict%node,attr,dict%key)

end subroutine search_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: searchi_ - search for a value
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine searchi_(indx,data, irank,dict,attr,verify)
      use m_die,only : die,perr
      use m_Attributes,only : Attributes
      use m_Attributes,only : key
      use m_Attributes,only : aSuperset
      implicit none
      integer         ,intent(in) :: indx	! data index
      type(Attributes),intent(in) :: data	! data container

      integer,intent(out) :: irank		! location as its rank

      type(Dictionary),intent(in) :: dict	! dictionary
      type(Attributes),intent(in) :: attr	! container for dict

      logical,optional,intent(in) :: verify

! !REVISION HISTORY:
! 	16Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::searchi_'
  integer :: kdata,kdict,kattr
  logical :: notin_data,notin_attr
  type(DictionaryNode),pointer :: node

  if(present(verify)) then
    if(verify) then

      if( dict%status /= status_INSERTBL	.and.	&
	  dict%status /= status_COMPILED	)	&
		call die(myname_,'not searchable')

      kdata=key(data)
      kdict=key_(dict)
      kattr=key(attr)

      notin_data=.not.aSuperset(kdata,kdict)
      notin_attr=.not.aSuperset(kattr,kdict)

      if(notin_data .or. notin_attr) then

	if(notin_data) call perr(myname_,	&
		'key(data)',kdata,'key(dict)',kdict)
	if(notin_attr) call perr(myname_,	&
		'key(attr)',kattr,'key(dict)',kdict)
	call die(myname_)
      endif
    endif
  endif

  if(dict%status == status_UNSRCHBL)	&
	call die(myname_,'dictionary not searchable')

  call rsearch_(indx,data, node, dict%node,attr,dict%key)

  irank=0
  if(associated(node)) irank=node%irank

	nullify(node)

	if(irank<=0) call die(myname_,'failed',indx)

end subroutine searchi_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: rsearch_ - locate a value on the tree
!
! !DESCRIPTION:
!
! !INTERFACE:

    recursive subroutine rsearch_(indx,data, where, node,attr,key)
      use m_die,only : die
      use m_Attributes,only : Attributes	! the container type
      use m_Attributes,only : order
      use m_Attributes,only : LESS,GREATER
      implicit none
      integer         ,intent(in) :: indx	! data index
      type(Attributes),intent(in) :: data	! data container

      type(DictionaryNode),pointer    :: where	! (out)

      type(DictionaryNode),pointer    :: node	! (in)
      type(Attributes)    ,intent(in) :: attr	! container for dict
      integer             ,intent(in) :: key		! searching key

! !REVISION HISTORY:
! 	16Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::rsearch_'
  integer :: iorder

  nullify(where)
  if(.not.associated(node)) return	! nowhere to find

  iorder=order(indx,data,node%irank,attr,key)
  select case(iorder)
  case(LESS)	! data(indx) < attr(node%irank), left-search

      call rsearch_(indx,data, where, node%left ,attr,key)

  case(GREATER)	! data(indx) > attr(node%irank), right-search

      call rsearch_(indx,data, where, node%right,attr,key)

  case default	! data(indx) == attr(node%irank), this is the one

    where => node
  end select
end subroutine rsearch_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initInsertbl_ - initialize a Dictionary as _INSERTBL
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine initInsertbl_(dict,key)
      implicit none
      type(Dictionary),intent(out) :: dict
      integer,intent(in) :: key		! a code for ordering, defined
					! by the given container

! !REVISION HISTORY:
! 	16Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::initInsertbl_'

  nullify(dict%node)
  dict%dsize=-1

  dict%status=status_INSERTBL
  dict%key=key
end subroutine initInsertbl_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initPackable_ - initialize a Dictionary with compiled data
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine initPackable_(dict,key,nsize)
      implicit none
      type(Dictionary),intent(out) :: dict
      integer,intent(in) :: key		! a code for ordering, defined
					! by the given container
      integer,intent(in) :: nsize	! thesize of the compiled data

! !REVISION HISTORY:
! 	05Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::initPackable_'
  integer :: ier

  nullify(dict%node)
  dict%dsize=-1
  dict%status=status_INSERTBL

  call alloc_(dict,nsize)

  dict%status=status_UNSRCHBL
  dict%key=key
end subroutine initPackable_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: alloc_ - allocate storate for packing
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine alloc_(dict,nsize)
      use m_mall,only : mall_ison,mall_mci
      use m_die ,only : perr,die
      implicit none
      type(Dictionary),intent(inout) :: dict
      integer,intent(in) :: nsize

! !REVISION HISTORY:
! 	05Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::alloc_'
  integer :: ier

  if(dict%status /= status_INSERTBL)	&
	call die(myname_,'unexpected object status')

  allocate(dict%counts(nsize), stat=ier)
	if(ier/=0) call die(myname_,'allocate()',ier)

	if(mall_ison()) call mall_mci(dict%counts ,myname)

end subroutine alloc_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean a Dictionary
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(dict)
      use m_die,only : die
      implicit none
      type(Dictionary),intent(inout) :: dict

! !REVISION HISTORY:
! 	03Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'

  select case(dict%status)
  case(status_INSERTBL)
        		! clean the tree
    call rcleanTree_(dict%node)

  case(status_UNSRCHBL)
        		! clean the packed data
    call cleanData_(dict)

  case(status_COMPILED)
        		! clean the packed data
    call cleanData_(dict)

        		! and clean the tree
    call rcleanTree_(dict%node)

  case default
    call die(myname_,'object undefined')
  end select

  nullify(dict%node)
  dict%dsize=-1
  dict%status=status_INSERTBL
  	! dict%key is left untouched!
end subroutine clean_
    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: cleanData_ - clean the memory for packed data
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine cleanData_(dict)
      use m_mall,only : mall_ison,mall_mco
      use m_die,only : die
      implicit none
      type(Dictionary),intent(inout) :: dict

! !REVISION HISTORY:
! 	05Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::cleanData_'
  integer :: ier

	if(mall_ison()) call mall_mco(dict%counts ,myname)
  deallocate(dict%counts, stat=ier)
	if(ier/=0) call die(myname_,'deallocate()',ier)

end subroutine cleanData_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: rcleanTree_ - clean up a whole tree
!
! !DESCRIPTION:
!
! !INTERFACE:

    recursive subroutine rcleanTree_(node)
      use m_die,only : die
      use m_mall,only : mall_ison,mall_co
      implicit none
      type(DictionaryNode),pointer :: node

! !REVISION HISTORY:
! 	16Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::rcleanTree_'

  integer :: ier

  if(.not.associated(node)) return

                                ! clean the branches
  call rcleanTree_(node%left )
  call rcleanTree_(node%right)

  nullify(node%left)
  nullify(node%right)

                                ! remove the node
	if(mall_ison()) call mall_co(1,myname)
  deallocate(node,stat=ier)
	if(ier/=0) call die(myname_,'deallocate(node)',ier)

end subroutine rcleanTree_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: size_ - get the unknown size of a Dictionary
!
! !DESCRIPTION:
!
! !INTERFACE:

    function size_(dict)
      use m_die,only : die
      implicit none
      type(Dictionary),intent(in) :: dict
      integer :: size_

! !REVISION HISTORY:
! 	16Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::size_'

  select case(dict%status)
  case(status_INSERTBL)
    size_=0
    call rsize_(dict%node,size_)
  case(status_UNSRCHBL)
    size_=size(dict%counts,1)
  case(status_COMPILED)
    size_=dict%dsize
  case default
    call die(myname_,'object undefined')
  end select

end function size_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: compiled_ - if the dictionary is status_COMPILED
!
! !DESCRIPTION:
!
! !INTERFACE:

    function compiled_(dict)
      implicit none
      type(Dictionary),intent(in) :: dict
      logical :: compiled_

! !REVISION HISTORY:
! 	24Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::compiled_'

  compiled_ = dict%status == status_COMPILED

end function compiled_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: key_ - get the ordering code of a dictionary
!
! !DESCRIPTION:
!
! !INTERFACE:

    function key_(dict)
      implicit none
      type(Dictionary),intent(in) :: dict
      integer :: key_

! !REVISION HISTORY:
! 	20Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::key_'

  key_=dict%key

end function key_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: rsize_ - count all nodes on a DictionaryNode
!
! !DESCRIPTION:
!
! !INTERFACE:

    recursive subroutine rsize_(node,n)
      implicit none
      type(DictionaryNode),pointer :: node
      integer,intent(inout) :: n

! !REVISION HISTORY:
! 	16Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::rsize_'

  if(.not.associated(node)) return

  call rsize_(node%left, n)
  n=n+1
  call rsize_(node%right,n)

end subroutine rsize_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: get_ - get ranked attributes
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine get_(dictent,irank,count)
      use m_die,only : die
      implicit none
      type(DictionaryEntry),intent(in)  :: dictent	! an entry
      integer,optional,intent(out) :: irank
      integer,optional,intent(out) :: count

! !REVISION HISTORY:
! 	03Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::get_'
  integer :: n
  type(DictionaryNode),pointer :: node

  node => dictent%node
  if(.not.associated(node)) call die(myname_,'an invalid node')

  if(present(irank)) irank=node%irank
  if(present(count)) count=node%count

  nullify(node)

end subroutine get_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_counts_ - referning to %counts (status /= INSERTBL)
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_counts_(dict)
      use m_die,only : die
      implicit none
      type(Dictionary),intent(in) :: dict
      integer,pointer,dimension(:) :: ptr_counts_

! !REVISION HISTORY:
! 	05Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_counts_'

  if( dict%status /= status_UNSRCHBL .and.	&
      dict%status /= status_COMPILED	)	&
	call die(myname_,'%status invalid')

  ptr_counts_ => dict%counts(:)

end function ptr_counts_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: compile_ - compile an uncompiled Dictionary
!
! !DESCRIPTION:
!
!	An un-compiled Dictionary is either INSERTBL or UNSRCHBL
!
! !INTERFACE:

    subroutine compile_(dict,attr,data)
      use m_die,only : die
      use m_mall,only : mall_ison,mall_mci,mall_mco
      use m_Attributes,only : Attributes
      use m_Attributes,only : subset
      implicit none
      type(Dictionary),intent(inout) :: dict
      type(Attributes),optional,intent(out) :: attr ! sample container
      type(Attributes),optional,intent(in)  :: data ! source container

! !REVISION HISTORY:
! 	05Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::compile_'
  integer,allocatable,dimension(:) :: indices
  integer :: irank
  integer :: dsize
  integer :: ier

  select case(dict%status)
  case(status_INSERTBL)

	! Expecting both attr and data

    if(.not.(present(attr).and.present(data)))	&
	call die(myname_,'missing argument(s)')


    			! Rank an INSERTBL dictionary, and get the
			! size.
    irank=0
    call rrank_(irank,dict%node)
    dsize=irank		! now the dictionary size is known

			! Get the indices as data locations to the
			! source container.  Then store the data
			! in their "Dictionary" order.  i.e.
			! attr(:) = data(indices(:)).

	allocate(indices(dsize),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)
		if(mall_ison()) call mall_mci(indices,myname)

    irank=0
    call rindx_(irank,dict%node,indices)
    call subset(attr,dict%key,indices,data)

		if(mall_ison()) call mall_mco(indices,myname)
	deallocate(indices,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

    call alloc_(dict,dsize)

			! Get %counts(:) in order, then create a tree.
    irank=0
    call rpack_(irank,dict%node,dict%counts)

			! rlink() grows a tree and leaves on existent
			! nodes.
    call rlink_(dict)
    dict%dsize=dsize

  case(status_UNSRCHBL)

    if(present(attr).or.present(data))	&
	call die(myname_,'unexpected argument(s)')

			! Data is assumed be sorted, and %counts(:)
			! has been defined.

    dsize=size(dict%counts,1)

			! rgrow_() grows a new tree from given data

    call rgrow_(dict%node,1,dsize,dict%counts)
    dict%dsize=dsize

  case(status_COMPILED)
    call die(myname_,'object already compiled')
  case default
    call die(myname_,'object undefined')
  end select

  dict%status=status_COMPILED
end subroutine compile_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: rrank_ - rank entries recursively

    recursive subroutine rrank_(irank,node)
      implicit none

      integer,intent(inout) :: irank
      type(DictionaryNode),pointer :: node

! !REVISION HISTORY:
! 	03Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::rrank_'

  if(.not.associated(node)) return

  call rrank_(irank,node%left)		! Count all on my left wing

  irank=irank+1				! Count me
  node%irank=irank

  call rrank_(irank,node%right)		! Count all on my right wing

end subroutine rrank_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: rindx_ - get %index entries, recursively

    recursive subroutine rindx_(irank,node,indices)
      implicit none

      integer,intent(inout) :: irank
      type(DictionaryNode),pointer :: node
      integer,dimension(:),intent(inout) :: indices

! !REVISION HISTORY:
! 	03Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::rindx_'

  if(.not.associated(node)) return

  call rindx_(irank,node%left ,indices)

  irank=irank+1
  indices(irank)=node%index	! 

  call rindx_(irank,node%right,indices)

end subroutine rindx_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: rpack_ - recursively pack ranked dictionary entries
!
! !DESCRIPTION:
!
! !INTERFACE:

    recursive subroutine rpack_(irank,node,counts)
      implicit none
      integer,intent(inout) :: irank
      type(DictionaryNode),pointer :: node
      integer,dimension(:),intent(inout) :: counts

! !REVISION HISTORY:
! 	05Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::rpack_'

  if(.not.associated(node)) return

  call rpack_(irank,node%left ,counts)

	! Index overflow condition should not exist, since this 
	! routine is called within pack_() with all dimension
	! values checked.

  irank=irank+1
  counts(irank)  =node%count

  call rpack_(irank,node%right,counts)

end subroutine rpack_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: rgrow_ - recursively grow a balanced tree
!
! !DESCRIPTION:
!
! !INTERFACE:

    recursive subroutine rgrow_(node,nLeft,nRight,counts)

      use m_mall,only : mall_ison,mall_ci
      use m_die, only : die
      implicit none
      type(DictionaryNode),  pointer    :: node
      integer,intent(in) :: nLeft,nRight
      integer,dimension(:)  ,intent(in) :: counts

! !REVISION HISTORY:
! 	23Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::rgrow_'
  integer :: nMiddle
  integer :: ier

  nullify(node)
  if(nLeft > nRight) return

  nMiddle=(nLeft+nRight+1)/2

	allocate(node,stat=ier)
		if(ier/=0) call die(myname_,'allocate(node)',ier)
		if(mall_ison()) call mall_ci(1,myname)

  node%irank=nMiddle
  node%count=counts(nMiddle)

  call rgrow_(node%left, nLeft,    nMiddle-1,counts)
  call rgrow_(node%right,nMiddle+1,nRight,   counts)

end subroutine rgrow_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: rlink_ - remake a compiled Dictionary
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine rlink_(dict)
      use m_die,only : die
      use m_mall,only : mall_ison,mall_ci,mall_co
      implicit none
      type(Dictionary),intent(inout) :: dict

! !REVISION HISTORY:
! 	04Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::rlink_'
  type(DictionaryEntry),allocatable,dimension(:) :: list
  type(DictionaryNode),pointer :: node
  integer :: dsize
  integer :: ier
  integer :: k,n

  dsize=size(dict%counts,1)

	allocate(list(dsize),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)
		if(mall_ison()) call mall_ci(size(list),myname)

  n=0
  call r2list_(n,list,dict%node)	! holding up the memories

  do k=1,n                              ! redefine the contents
    node => list(k)%node
    node%irank =  k
    node%count =  dict%counts(k)
  end do
  nullify(node)

  call r2tree_(dict%node,1,n,list)	! redefine the pointers

		if(mall_ison()) call mall_co(size(list),myname)
	deallocate(list,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

end subroutine rlink_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: r2list_ - recursively list sorted values
!
! !DESCRIPTION:
!
! !INTERFACE:

    recursive subroutine r2list_(n,list,node)
      implicit none
      integer,                      intent(inout) :: n
      type(DictionaryEntry),dimension(:),intent(inout) :: list
      type(DictionaryNode),         pointer       :: node

! !REVISION HISTORY:
! 	16Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::r2list_'

  if(.not.associated(node)) return

  call r2list_(n,list,node%left )

	! Index overflow condition should not exist, since this 
	! routine is called within rlink() with all dimension
	! values checked.

  n=n+1
  list(n)%node => node

  call r2list_(n,list,node%right)

end subroutine r2list_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: r2tree_ - Note only the pointers are reassgned
!
! !DESCRIPTION:
!
! !INTERFACE:

    recursive subroutine r2tree_(node,nLeft,nRight,list)
      implicit none
      type(DictionaryNode),pointer :: node
      integer,                      intent(in) :: nLeft,nRight
      type(DictionaryEntry),dimension(:),intent(in) :: list

! !REVISION HISTORY:
! 	23Aug99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::r2tree_'
  integer :: nMiddle

  nullify(node)
  if(nLeft > nRight) return

  nMiddle=(nLeft+nRight+1)/2

  node => list(nMiddle)%node
  call r2tree_(node%left, nLeft,    nMiddle-1,list)
  call r2tree_(node%right,nMiddle+1,nRight,   list)

end subroutine r2tree_
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
      type(Dictionary),pointer :: mold
      integer,optional,intent(out) :: stat
      type(Dictionary),pointer :: new_

! !REVISION HISTORY:
! 	28Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::new_'
  type(Dictionary),pointer :: obj
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
      type(Dictionary),pointer :: obj
      integer,optional,intent(out) :: stat


! !REVISION HISTORY:
! 	28Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
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
end module m_Dictionary
!.
