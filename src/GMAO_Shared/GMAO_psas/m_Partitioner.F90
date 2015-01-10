!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_Partitioner - Partitioner-triple (dict,attr,rank)
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_Partitioner
      use m_Dictionary,only : Dictionary
      use m_Attributes,only : Attributes
      implicit none
      private	! except

      public :: Partitioner            ! Class data structure
      public :: Partitioner_init ,init
      public :: Partitioner_clean,clean
      public :: Partitioner_partition

      public :: ptr_dict
      public :: ptr_attr
      public :: ptr_rank
      
      type Partitioner
	private
			! A partitioner-triple contains a dictionary,
			! a data container of the dictionary, and the
			! ranking list of the dictionary entries

	logical :: aliasedDict 
        type(Dictionary),pointer :: dict        ! a shared dictionary
	type(Attributes),pointer :: attr	! a container for dict.

        integer,pointer,dimension(:) :: rank ! corresponding PEs
      end type Partitioner

      interface Partitioner_init; module procedure	&
	initbyAttr_,	&
	initbyDict_; end interface
      interface init; module procedure			&
	initbyAttr_,	&
	initbyDict_; end interface

      interface Partitioner_clean; module procedure	&
	clean_; end interface
      interface clean; module procedure clean_; end interface

      interface Partitioner_partition; module procedure	&
	partition_; end interface

      interface ptr_dict; module procedure ptr_dict_; end interface
      interface ptr_attr; module procedure ptr_attr_; end interface
      interface ptr_rank; module procedure ptr_rank_; end interface
        
! !REVISION HISTORY:
! 	09Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_Partitioner'

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initbyAttr_ - Create a Partitioner-triple
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine initbyAttr_(ptnr,comm,attr,keyCh,weights,wpower,nozero)
      use m_die,only : die
      use m_mall,only : mall_ison
      use m_mall,only : mall_mci

      use m_keyChain,only : keyChain

      use m_Attributes,only : Attributes
      use m_Attributes,only : new

      use m_Dictionary,only : Dictionary
      use m_Dictionary,only : Dictionary_size
      use m_Dictionary,only : new

      use m_SharedDictionary,only : Dictionary_make
      use m_SharedDictionary,only : Dictionary_rank

      implicit none

      type(Partitioner),intent(out) :: ptnr	! partitioner-triple
      integer          ,intent(in)  :: comm	! Communicator
      type(Attributes) ,intent(in)  :: attr	! attribute pool
      type(keyChain)   ,intent(in)  :: keyCh	! key list of subsets

				! If the data in attr have weights, ...

      integer,optional,dimension(:),intent(in) :: weights

				! If one would like to weight the
				! partition schem with a simple "power-
				! of-counts" scheme, specify wpower.

      real   ,optional,intent(in) :: wpower

				! If one does not want to have any PE
				! left empty, let nozero=.true..
				! However, if there is no way to make
				! all PEs happy, an exception (die) will
				! occur.

      logical,optional,intent(in) :: nozero

! !REVISION HISTORY:
! 	09Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::initbyAttr_'
  integer :: dsize_
  integer :: ier
  
  ptnr%aliasedDict=.false.
  ptnr%dict => new(ptnr%dict)
  ptnr%attr => new(ptnr%attr)

	! Dictionary_make() constructs a "dictionary" and computes the
	! weights for all entries.  The weight computed by _make() is
	! a linear weight.  This weight is obtained either through a
	! simple counting or a weighted counting using array weights(:).

  call Dictionary_make(ptnr%dict,ptnr%attr, attr,keyCh, comm,	&
	weights=weights)
  
  dsize_=Dictionary_size(ptnr%dict)
        allocate(ptnr%rank(dsize_),stat=ier)
                if(ier/=0) call die(myname_,'allocate()',ier)
                if(mall_ison()) call mall_mci(ptnr%rank,myname)
                
	! Dictionary_rank() computes a distribution table based on a
	! "weight" table.  This "weight" table can be either the weights
	! stored in the dictionary, or derived from a simple function of
	! of that weights.  If it is useful, one can even implement a
	! new argument to specify a totally independent weight table.
	! However, I don't think it makes sense to do that.

  call Dictionary_rank(ptnr%rank,ptnr%dict,comm,	&
	wpower=wpower,nozero=nozero)

end subroutine initbyAttr_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initbyDict_ - Create a Partitioner-triple
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine initbyDict_(ptnr,comm,dict,attr,wpower,nozero)
      use m_die,only : die
      use m_mall,only : mall_ison
      use m_mall,only : mall_mci

      use m_Attributes,only : Attributes

      use m_Dictionary,only : Dictionary
      use m_Dictionary,only : Dictionary_size

      use m_SharedDictionary,only : Dictionary_rank

      implicit none

      type(Partitioner),intent(out) :: ptnr	! partitioner-triple
      integer          ,intent(in)  :: comm	! Communicator
      type(Dictionary),target,intent(in)  :: dict
      type(Attributes),target,intent(in)  :: attr	! attribute pool

				! If one would like to weight the
				! partition schem with a simple "power-
				! of-counts" scheme, specify wpower.

      real   ,optional,intent(in) :: wpower

				! If one does not want to have any PE
				! left empty, let nozero=.true..
				! However, if there is no way to make
				! all PEs happy, an exception (die) will
				! occur.

      logical,optional,intent(in) :: nozero

! !REVISION HISTORY:
! 	09Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::initbyDict_'
  integer :: dsize_
  integer :: ier
  
  ptnr%aliasedDict=.true.
  ptnr%dict => dict
  ptnr%attr => attr

  dsize_=Dictionary_size(ptnr%dict)
        allocate(ptnr%rank(dsize_),stat=ier)
                if(ier/=0) call die(myname_,'allocate()',ier)
                if(mall_ison()) call mall_mci(ptnr%rank,myname)
                
	! Dictionary_rank() computes a distribution table based on a
	! "weight" table.  This "weight" table can be either the weights
	! stored in the dictionary, or derived from a simple function of
	! of that weights.  If it is useful, one can even implement a
	! new argument to specify a totally independent weight table.
	! However, I don't think it makes sense to do that.

  call Dictionary_rank(ptnr%rank,ptnr%dict,comm,	&
	wpower=wpower,nozero=nozero)

end subroutine initbyDict_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - Clean a Partitioner-triple
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(ptnr)
      use m_die ,only : die
      use m_mall,only : mall_ison
      use m_mall,only : mall_mco

      use m_Attributes,only : clean,delete
      use m_Dictionary,only : Dictionary_clean,delete

      implicit none

      type(Partitioner) ,intent(inout) :: ptnr

! !REVISION HISTORY:
! 	09Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ier

  if(.not. ptnr%aliasedDict) then
    call Dictionary_clean(ptnr%dict)
    call clean(ptnr%attr)

    call delete(ptnr%dict)
    call delete(ptnr%attr)
  else
    nullify(ptnr%dict)
    nullify(ptnr%attr)
  endif
  ptnr%aliasedDict=.false.

        if(mall_ison()) call mall_mco(ptnr%rank,myname)
  deallocate(ptnr%rank,stat=ier)
        if(ier/=0) call die(myname_,'deallocate()',ier)


end subroutine clean_ 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: partition_ - apply partitioner to a given data
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine partition_(attr, indx,tran,coll, ptnr,comm)
      use m_die,only : die
      use m_die,only : MP_die
      use m_mall,only : mall_ison,mall_mci,mall_mco
      use m_Attributes  ,only : Attributes
      use m_Attributes  ,only : lsize
      use m_SharedDictionary,only : Dictionary_lookup
      use m_Transposer,only : Transposer
      use m_Transposer,only : Transposer_init
      use m_Collector ,only : Collector
      use m_Collector ,only : Collector_init
      use m_mpif90,only : MP_comm_size

      implicit none

      type(Attributes)    ,intent(in)  :: attr

				! partition-triple
      integer,dimension(:),intent(out) :: indx	! (lsize(attr))
      type(Transposer)    ,intent(out) :: tran
      type(Collector)     ,intent(out) :: coll

      type(Partitioner)   ,intent(in)  :: ptnr
      integer             ,intent(in)  :: comm

! !REVISION HISTORY:
! 	09Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::partition_'
  integer,allocatable,dimension(:) :: scounts
  integer :: lsize_
  integer :: ier
  integer :: nPE

  call MP_comm_size(comm,nPE,ier)
	if(ier/=0) call MP_die(myname_,'MP_comm_size()',ier)

	allocate(scounts(0:nPE-1),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)
		if(mall_ison()) call mall_mci(scounts,myname)

  call Dictionary_lookup(attr, indx, scounts,	&
	ptnr%dict,ptnr%attr,ptnr%rank)

  call Transposer_init(tran,scounts,comm)

		if(mall_ison()) call mall_mco(scounts,myname)
	deallocate(scounts,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

  lsize_=lsize(attr)

  call Collector_init(coll,lsize_,comm)
end subroutine partition_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_dict_ - referencing %dict component
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_dict_(obj)
      use m_Dictionary,only : Dictionary
      implicit none
      type(Partitioner),intent(in) :: obj
      type(Dictionary),pointer :: ptr_dict_

! !REVISION HISTORY:
! 	25Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_dict_'

  ptr_dict_ => obj%dict

end function ptr_dict_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_attr_ - referencing %attr component
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_attr_(obj)
      use m_Attributes,only : Attributes
      implicit none
      type(Partitioner),intent(in) :: obj
      type(Attributes),pointer :: ptr_attr_

! !REVISION HISTORY:
! 	25Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_attr_'

  ptr_attr_ => obj%attr

end function ptr_attr_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_rank_ - referencing %rank component
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_rank_(obj)
      implicit none
      type(Partitioner),intent(in) :: obj
      integer,pointer,dimension(:) :: ptr_rank_

! !REVISION HISTORY:
! 	25Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_rank_'

  ptr_rank_ => obj%rank

end function ptr_rank_

end module m_Partitioner
