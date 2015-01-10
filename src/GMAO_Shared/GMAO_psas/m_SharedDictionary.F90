!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_SharedDictionary - Tools for shared Dictionary
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_SharedDictionary
      implicit none
      private	! except

      public :: Dictionary_make
      public :: Dictionary_rank
      public :: Dictionary_lookup

      interface Dictionary_make; module procedure       &
        recurDict_; end interface
      interface Dictionary_rank; module procedure       &
        rank_; end interface
      interface Dictionary_lookup; module procedure     &
        lookitup_; end interface
        
! !REVISION HISTORY:
! 	04Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_SharedDictionary'
!
!   Given a (sortingKeyList) keyList, an (AttrVect) attr, a
!   distribution-triple is returned through variables dict:attr:rank.
!
!	call Dictionary_make(data,keyList, dict,attr,comm)
!
!		allocate(rank(size(dict)))
!	call Dictionary_rank(rank,dict,comm)
!       KeySet=key(keyList)
!
!		allocate(scounts(0:nPE-1))
!		allocate(indx(lsize(data))
!	call Dictionary_lookup(data, indx,scounts,   &
!               dict,attr,rank, comm)

!
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: recurDict_ - recursively build a shared Dictionary
!
! !DESCRIPTION:
!
! !INTERFACE:

    recursive subroutine recurDict_(dict,attr, data,keyCh,	&
	comm,weights,verify)

      use m_die ,only : die
      use m_mall,only : mall_ison,mall_mci,mall_mco

      use m_keyChain,only : keyChain
      use m_keyChain,only : ptr_next
      use m_keyChain,only : key

      use m_Attributes,only : Attributes
      use m_Attributes,only : clean

      use m_Dictionary,only : Dictionary
      use m_Dictionary,only : Dictionary_size
      use m_Dictionary,only : Dictionary_clean

      use m_Transposer,only : Transposer
      use m_Transposer,only : Transposer_init
      use m_Transposer,only : clean

      implicit none

      type(Dictionary),intent(out) :: dict	! shared dictionary
      type(Attributes),intent(out) :: attr	! dictionary container

      type(Attributes),intent(in)  :: data	! distributed containers
      type(keyChain)  ,intent(in)  :: keyCh	! a list of subsets

      integer         ,intent(in)  :: comm	! communicator

      integer,dimension(:),optional,intent(in) :: weights
      logical,optional,intent(in) :: verify	! verify keys

		! If an entry given by data(:) should be weighted,
		! use optional argument weights(:)

! !REVISION HISTORY:
! 	04Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::recurDict_'

  type(keyChain),pointer :: keyChNext

  type(Dictionary) :: dict_local
  type(Attributes) :: attr_local

  type(Dictionary) :: dictNext
  type(Attributes) :: attrNext
  integer,allocatable,dimension(:) :: rankNext

  integer :: dsize_
  integer :: ier

!________________________________________
! Aquire a distribution lookup table

	! Sampling the attribute data to compile a Dictionary.  The
	! current sorting key for the dictionary is obtained from
	! key(keyCh)

  call sample_(data,key(keyCh), dict_local,attr_local,	&
	weights=weights,verify=verify)

  keyChNext => ptr_next(keyCh)  ! The keyChain without this key

  if(associated(keyChNext)) then

        !_____________________________________________________________
        ! Creating a global dictionary of attr_local in parallel
	! requires a distribution-triple (dict,rank,attr) to
	! determine how the local merging work should be distributed.

                ! .. Get the next distribution, dict:attr:rank.
		! Note since attr_local is the actual list of entries,
		! weighting is not needed.

    call recurDict_(dictNext,attrNext, attr_local,keyChNext,	&
	comm,verify=verify)
    
        !_____________________________________________________________
                ! .. Get the next distribution.rank.  This step
		! determine the "distribution".

	dsize_=Dictionary_size(dictNext)
	allocate(rankNext(dsize_),stat=ier)
		if(ier/=0) call die(myname_,'allocate(rankNext)',ier)
		if(mall_ison()) call mall_mci(rankNext,myname)

    call rank_(rankNext,dictNext,comm)

        !_____________________________________________________________
        ! Merge dict_local:attr_local to a shared Dictionary,
	! dict:attr
        
    call mergex_(dict_local,attr_local, dict,attr,	&
	comm, DictNext,attrNext,rankNext,verify=verify)
    
        !_____________________________________________________________
                ! Remove dictNext:attrNext:rankNext

    call Dictionary_clean(dictNext)
    call clean(attrNext)

		if(mall_ison()) call mall_mco(rankNext,myname)
	deallocate(rankNext,stat=ier)
		if(ier/=0) call die(myname_,'deallocate(rankNext)',ier)

  else
        !_____________________________________________________________
	! Collect all dict_local and their counts by all PEs

    call merge1_(dict_local,attr_local, dict,attr, comm,verify=verify)

  endif
  
  nullify(keyChNext)			! remove a temporary keyChain

  call Dictionary_clean(dict_local)	! remove a temporary Dictionary
  call clean(attr_local)

end subroutine recurDict_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: rank_ - partition dictionary weighted by their counts
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine rank_(rank,dict,comm,wpower,nozero)
      use m_die,only : MP_die
      use m_elemPart,only : elemPart

      use m_Dictionary,only : Dictionary
      use m_Dictionary,only : ptr_counts

      use m_Collector,only : Collector
      use m_Collector,only : Collector_init => init
      use m_Collector,only : clean

      use m_mpif90,only : MP_comm_rank
      use m_CollectorComm,only : allgatherv
      use m_mpout, only : mpout_log, mpout
      implicit none

      integer,dimension(:),intent(out) :: rank	! PE map
      type(Dictionary)    ,intent(in)  :: dict	! key list
      integer,intent(in) :: comm

      real   ,optional,intent(in) :: wpower	! specify a weight as
						! a simple power of
						! total counts.
      logical,optional,intent(in) :: nozero

! !REVISION HISTORY:
! 	06Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::rank_'
  integer,pointer,dimension(:) :: weights
  type(Collector) :: coll
  integer :: lbound,ubound
  integer :: ier
  integer :: i
  integer :: iPE
  real    :: wpower_

	! It makes a better sense to have wpower=1. as the default.
	! wpower=2. is used to make it compatable with m-r3p3.

  wpower_=1.
  if(present(wpower)) wpower_=wpower

	! Determine the local size on this PE

  weights => ptr_counts(dict)
  call elemPart(size(weights),real(weights)**wpower_,	&
	lbound,ubound,comm,nozero=nozero)
  nullify(weights)



  call MP_comm_rank(comm,iPE,ier)
	if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

		! Create a collector from all local sizes

  call Collector_init(coll,ubound-lbound+1,comm)
  call allgatherv((/(iPE,i=lbound,ubound)/),rank,coll,comm)
  call clean(coll)
end subroutine rank_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: lookitup_ - determine distribution with a distribution
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine lookitup_(data, indx,scounts, dict,attr,rank,verify)

      use m_Attributes,only : Attributes
      use m_Attributes,only : lsize
      use m_Attributes,only : key
      use m_Attributes,only : aSuperset

      use m_Dictionary,only : Dictionary
      use m_Dictionary,only : DictionaryEntry
      use m_Dictionary,only : Dictionary_search
      use m_Dictionary,only : Dictionary_key

      use m_SortingTools,only : indexSet
      use m_SortingTools,only : indexSort

      use m_die ,only : die,perr
      use m_mall,only : mall_ison,mall_mci,mall_mco
      implicit none

      type(Attributes)     ,intent(in ) :: data		! data
      integer,dimension(:) ,intent(out) :: indx		! (lsize(data))
      integer,dimension(0:),intent(out) :: scounts	! (0:nPE-1)

		! Distribution-triple

      type(Dictionary)    ,intent(in) :: dict	! dictionary
      type(Attributes)    ,intent(in) :: attr	! container
      integer,dimension(:),intent(in) :: rank	! ranking by PE

      logical,optional,intent(in) :: verify

! !REVISION HISTORY:
!	05Sep01	- Jing Guo
!		. Removed weights(:).  I believe the earlier
!		  implementition was put in a wrong place.
! 	07Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::lookitup_'
  integer,allocatable,dimension(:) :: krank
  type(DictionaryEntry) :: node
  integer :: kdata,kattr,kdict
  logical :: notin_data,notin_attr

  type(Attributes) :: subset
  integer :: lsize_
  integer :: nPE
  integer :: iPE
  integer :: i,k
  integer :: ier

	! Check for invariants

	if(present(verify)) then
	  if(verify) then
  	    kdict=Dictionary_key(dict)
	    kdata=key(data)
	    kattr=key(attr)

	    notin_data=.not.aSuperset(kdata,kdict)
	    notin_attr=.not.aSuperset(kattr,kdict)

	    if(notin_data.or.notin_attr) then
	      if(notin_data) call perr(myname_,	&
		'key(data)',kdata,'key(dict)',kdict)
	      if(notin_attr) call perr(myname_,	&
		'key(attr)',kattr,'key(dict)',kdict)
	      call die(myname_)
	    endif
	  endif
	endif

  lsize_=lsize(data)
	allocate( krank(lsize_),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)
		if(mall_ison()) call mall_mci(krank,myname)

  nPE=size(scounts)
  scounts(0:nPE-1)=0

  do i=1,lsize_
    call Dictionary_search(i,data, k ,dict,attr)

	! Counting of scounts(:) is done by bin-counting, not by
	! seg-counting.  The latter would exclude any bin with a zero
	! count.

    iPE=rank(k)
    krank(i)=iPE
    if(iPE<0.or.iPE>=nPE) call die(myname_,'unexpected rank',iPE)

    scounts(iPE)=scounts(iPE)+1
  end do

		! Index-sort the data by their PE ranks.  Note that it
		! is not required the PE-ranks to be in the same order
		! as the dictionary, although the partition pattern
		! (rank(:)) is normally determined from the dictionary.

  call indexSet(indx(1:lsize_))
  call indexSort(indx(1:lsize_),krank)

		if(mall_ison()) call mall_mco(krank,myname)
	deallocate(krank,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

end subroutine lookitup_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: mergex_ - Merge local Dictionaries to a shared version
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine mergex_(dict_local,attr_local, dict,attr, comm,	&
	dictDstr,attrDstr,rankDstr,verify)

      use m_Dictionary,only : Dictionary
      use m_Dictionary,only : Dictionary_size
      use m_Dictionary,only : Dictionary_clean

      use m_Attributes,only : Attributes
      use m_Attributes,only : clean

      use m_Transposer,only : Transposer
      use m_Transposer,only : Transposer_init
      use m_Transposer,only : Transposer_clean

      use m_die   ,only : die,MP_die
      use m_mall  ,only : mall_ison,mall_mci,mall_mco
      use m_mpif90,only : MP_comm_size

      implicit none

			! Local dictionary with its container

      type(Dictionary),intent(in) :: dict_local ! dictionary
      type(Attributes),intent(in) :: attr_local	! container

			! Merged dictionary with its container

      type(Dictionary),intent(out) :: dict	! dictionary
      type(Attributes),intent(out) :: attr	! container

      integer,intent(in) :: comm		! communicator

			! Distribution-triple

      type(Dictionary)    ,intent(in) :: dictDstr	! dictionary
      type(Attributes)    ,intent(in) :: attrDstr	! container
      integer,dimension(:),intent(in) :: rankDstr	! (PE) rank

      logical,optional,intent(in) :: verify	! verify keys

! !REVISION HISTORY:
! 	08Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::mergex_'

  integer :: lsize_
  integer :: nPE
  integer :: ier

  integer,allocatable,dimension(:) :: indx
  integer,allocatable,dimension(:) :: scounts

  type(Dictionary) :: dict_distr
  type(Attributes) :: attr_distr
  type(Transposer) :: tran

	lsize_=Dictionary_size(dict_local)
	call MP_comm_size(comm,nPE,ier)
	        if(ier/=0) call MP_die(myname_,'MP_comm_size()',ier)

	allocate(indx(lsize_),scounts(0:nPE-1),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)
		if(mall_ison()) then
		  call mall_mci(indx   ,myname)
		  call mall_mci(scounts,myname)
		endif
		
	! Look the keys of sendd up in dict

    call lookitup_(attr_local, indx,scounts,	&
	dictDstr,attrDstr,rankDstr,verify=verify)

        call Transposer_init(tran,scounts,comm)

	! Transpose the contents of local dictionary sendd to their
	! desitation PEs, and merge into a distributed dictionary.
	! Note that local dictionaries on all PEs may have multiple 
        ! entries with the same keys, while the distributed dictionary
        ! dict_distr has uniq keys globally.

    call distrmerge_(dict_local,attr_local,	&
	dict_distr,attr_distr, indx,tran,comm,verify=verify)

	! Populate the distributed dictionary to form a shared global
	! Dictionary.

    call allgather_(dict_distr,attr_distr, dict,attr, comm)

    call Dictionary_clean(dict_distr)
    call clean(attr_distr)

	call Transposer_clean(tran)

		if(mall_ison()) then
		  call mall_mco(indx   ,myname)
		  call mall_mco(scounts,myname)
		endif
	deallocate(indx,scounts,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

end subroutine mergex_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: merge1_ - Merge local Dictionaries to a shared version
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine merge1_(dict_local,attr_local, dict,attr, comm,verify)
      use m_Dictionary,only : Dictionary
      use m_Dictionary,only : Dictionary_size
      use m_Dictionary,only : Dictionary_key
      use m_Dictionary,only : ptr_counts

      use m_Collector ,only : Collector
      use m_Collector ,only : Collector_init
      use m_Collector ,only : Collector_clean
      use m_Collector ,only : globalSize
      use m_CollectorComm ,only : allgatherv

      use m_Attributes,only : Attributes
      use m_Attributes,only : clean
      use m_AttributesComm,only : allgatherv

      use m_die ,only : die
      use m_mall,only : mall_ison,mall_mci,mall_mco

      implicit none

      type(Dictionary),intent(in)  :: dict_local
      type(Attributes),intent(in)  :: attr_local

      type(Dictionary),intent(out) :: dict
      type(Attributes),intent(out) :: attr

      integer,intent(in) :: comm
      logical,optional,intent(in) :: verify

! !REVISION HISTORY:
! 	06Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//':: merge1_'

  type(Collector) :: coll

  type(Attributes) :: attr_coll
  integer,allocatable,dimension(:) :: counts_coll

  integer :: dsize_
  integer :: gsize_
  integer :: ier

		! Find the pattern of the initial distribution

    dsize_=Dictionary_size(dict_local)
    call Collector_init(coll,dsize_,comm)
    gsize_=globalSize(coll)

	allocate(counts_coll(gsize_),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)
		if(mall_ison()) call mall_mci(counts_coll,myname)

    call allgatherv(ptr_counts(dict_local),counts_coll,coll,comm)

    call allgatherv(attr_local,attr_coll,coll,comm)

		! The merged dictionary may have a size different
		! from attr_coll(:), since there may be duplicate key
		! entries in attr_coll(:).  When duplicates occur,
		! sample_() will count in their weights through
		! argument weights=counts_coll(:).

    call sample_(attr_coll,Dictionary_key(dict_local),	&
	dict,attr, weights=counts_coll,verify=verify)

	call clean(attr_coll)

		if(mall_ison()) call mall_mco(counts_coll,myname)
	deallocate(counts_coll,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

    call Collector_clean(coll)

end subroutine merge1_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: distrmerge_ - Merge local dictionaries in distribution
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine distrmerge_(dict_local,attr_local,	&
	dict_distr,attr_distr, indx,tran,comm,verify)

      use m_Transposer,only : Transposer
      use m_Transposer,only : recvSize
      use m_TransposerComm,only : indexedTrans

      use m_Dictionary,only : Dictionary
      use m_Dictionary,only : Dictionary_key
      use m_Dictionary,only : ptr_counts

      use m_Attributes,only : Attributes
      use m_Attributes,only : clean
      use m_AttributesComm,only : indexedTrans

      use m_die ,only : die
      use m_mall,only : mall_ison,mall_mci,mall_mco

      implicit none
			! local dictionary and its data container

      type(Dictionary),intent(in)  :: dict_local
      type(Attributes),intent(in)  :: attr_local

			! distributed dictionary and its data containyer

      type(Dictionary),intent(out) :: dict_distr
      type(Attributes),intent(out) :: attr_distr

      integer,dimension(:),intent(in) :: indx	! to dictionary entries
      type(Transposer)    ,intent(in) :: tran	! all-to-all pattern
      integer,intent(in) :: comm		! communicator
      logical,optional,intent(in) :: verify	! verify keys

! !REVISION HISTORY:
! 	06Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::distrmerge_'

        ! Content buffers of the distributed dictionary 

  type(Attributes) :: attr_trans
  integer,allocatable,dimension(:) :: counts_trans

        ! Contents pointers of the local dictionary

  integer,pointer,dimension(:) ::  counts_local

  integer :: psize_
  integer :: ier

	! Prepare workspace

                        ! Alias send buffers
                        
                 counts_local => ptr_counts(dict_local)

                        ! Prepare recv buffers

	psize_=recvSize(tran)

	allocate(counts_trans(psize_),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)
		if(mall_ison()) call mall_mci(counts_trans,myname)

        ! Transpose the contents of local dictionaries to their
        ! destinations.

    call indexedTrans(indx,counts_local,counts_trans,tran,comm)

    call indexedTrans(indx,  attr_local,  attr_trans,tran,comm)

	! Merge in distribution to build a distributed dictionary
	! with a container, dict_distr:attr_distr, by compiling
	! attr_trans with their weights counts_trans.  The reason of
	! using weights is that this sample_() process is a merging
	! process of locally compiled dictionaries.

    call sample_(attr_trans,Dictionary_key(dict_local),	&
	dict_distr,attr_distr, weights=counts_trans,verify=verify)

        ! Clean workspace and aliases.

	nullify(counts_local)

	call clean(attr_trans)

		if(mall_ison()) call mall_mco(counts_trans,myname)
	deallocate(counts_trans,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

end subroutine distrmerge_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: allgather_ - Gather distributed Dictionaries to all PEs
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine allgather_(dict_distr,attr_distr,dict,attr,comm)
      use m_Dictionary,only : Dictionary
      use m_Dictionary,only : Dictionary_init
      use m_Dictionary,only : Dictionary_size
      use m_Dictionary,only : Dictionary_key
      use m_Dictionary,only : Dictionary_compile
      use m_Dictionary,only : ptr_counts
      use m_Dictionary,only : Dictionary_compiled

      use m_Collector,only : Collector
      use m_Collector,only : Collector_init
      use m_Collector,only : Collector_clean
      use m_Collector,only : globalSize
      use m_CollectorComm,only : allgatherv

      use m_Attributes,only : Attributes
      use m_AttributesComm,only : allgatherv

      use m_die,only : die

      implicit none
      type(Dictionary),intent(in)  :: dict_distr
      type(Attributes),intent(in)  :: attr_distr
      type(Dictionary),intent(out) :: dict
      type(Attributes),intent(out) :: attr
      integer         ,intent(in)  :: comm

! !REVISION HISTORY:
! 	06Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::allgather_'
  integer :: lsize_
  integer :: gsize_
  type(Collector) :: coll

  integer,pointer,dimension(:) :: indx_d,cntx_d
  integer,pointer,dimension(:) :: indx_c,cntx_c
  type(Attributes) :: attr_c

! Note that local segments of a distributed dictionary are not only
! sorted and unique locally, but also sorted and unique globally.
! Therefore, one can get the shared dictionary by allgather().

  lsize_=Dictionary_size(dict_distr)

  call Collector_init(coll,lsize_,comm)
  gsize_=globalSize(coll)

  if(.not.Dictionary_compiled(dict_distr))		&
	call die(myname_,'uncompiled dict_distr')

                ! Initialize a buffered Dictionary

  call Dictionary_init(dict,Dictionary_key(dict_distr),gsize_)

		! Set aliases to the components of the output

	 cntx_c => ptr_counts(dict)

		! Set aliases to the components of the input

	 cntx_d => ptr_counts(dict_distr)

    call allgatherv(attr_distr,attr,coll,comm)

	! It is expected that attr_distr are globally sorted such that
	! collected attr is therefore properly sorted.  The dictionary
	! will be compiled based on a given %counts(:) with the same
	! sorted order.

    call allgatherv(cntx_d,cntx_c,coll,comm)

	nullify(cntx_d)
	nullify(cntx_c)

  call Collector_clean(coll)

		! Compile the output Dictionary

  call Dictionary_compile(dict)

end subroutine allgather_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: sample_ - compile a dictionary of attribute values
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine sample_(data,key, dict,attr,weights,verify)
      use m_Dictionary,only : Dictionary
      use m_Dictionary,only : Dictionary_init
      use m_Dictionary,only : Dictionary_insert
      use m_Dictionary,only : Dictionary_compile

      use m_Attributes,only : Attributes
      use m_Attributes,only : lsize
      use m_Attributes,only : Attributes_key => key
      use m_Attributes,only : aSuperset

      use m_die,only : die

      implicit none
      type(Attributes),intent(in)  :: data	! data container
      integer         ,intent(in)  :: key	! code of subsetting

      type(Dictionary),intent(out) :: dict	! as a Dictionary
      type(Attributes),intent(out) :: attr	! container for dict.

      integer,dimension(:),optional,intent(in) :: weights
      logical,optional,intent(in) :: verify	! verify keys

! !REVISION HISTORY:
! 	06Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::sample_'
  integer :: kdata
  integer :: i

		! Check invariant

	if(present(verify)) then
	  if(verify) then
	    kdata=Attributes_key(data)
	    if(.not.aSuperset(kdata,key)) call die(myname_,	&
		'key(data)',kdata,'key',key)
	  endif
	endif

		! Creating a Dictionary

  call Dictionary_init(dict,key)

  if(present(weights)) then
    do i=1,lsize(data)
      call Dictionary_insert(dict,i,data,count=weights(i))
    end do
  else
    do i=1,lsize(data)
      call Dictionary_insert(dict,i,data)
    end do
  endif

		! Compiling _dict_ results an _attr_ with
		! information from _data_.

  call Dictionary_compile(dict,attr=attr,data=data)

end subroutine sample_

end module m_SharedDictionary
